#include <MultiGridLevel.h>

#include <Connectivity.h>
#include <Domain.h>
#include <Edge.h>

#include <set>
#include <mpi.h>

#include <GeoSource.h>

#include <DebugTools.h>

template <class Scalar>
CommPattern<Scalar>* internal_select(CommPattern<int>* a, CommPattern<double>* b) {

  std::cout << "Error: should not be called!" << std::endl;
  return NULL;
}

template <>
CommPattern<int>* internal_select<int>(CommPattern<int>* a, CommPattern<double>* b) {

  return a;
}

template <>
CommPattern<double>* internal_select<double>(CommPattern<int>* a, CommPattern<double>* b) {

  return b;
}

template <class Scalar>
const CommPattern<Scalar>* internal_select(const CommPattern<int>* a, const CommPattern<double>* b) {

  std::cout << "Error: should not be called!" << std::endl;
  return NULL;
}

template <>
const CommPattern<int>* internal_select<int>(const CommPattern<int>* a,const  CommPattern<double>* b) {

  return a;
}

template <>
const CommPattern<double>* internal_select<double>(const CommPattern<int>* a,const  CommPattern<double>* b) {

  return b;
}



template<class Scalar>
MultiGridLevel<Scalar>::MultiGridLevel(MultiGridMethod mgm,MultiGridLevel* mg, Domain& domain, DistInfo& refinedNodeDistInfo, DistInfo& refinedEdgeDistInfo)
  : nodeIdPattern(0), nodeVolPattern(0), nodePosnPattern(0),nodeVecPattern(0), domain(domain),
    nodeDistInfo(new DistInfo(refinedNodeDistInfo.numLocThreads, refinedNodeDistInfo.numLocSub, refinedEdgeDistInfo.numGlobSub,
                              refinedNodeDistInfo.locSubToGlobSub, refinedNodeDistInfo.com)),
    edgeDistInfo(new DistInfo(refinedEdgeDistInfo.numLocThreads, refinedEdgeDistInfo.numLocSub, refinedEdgeDistInfo.numGlobSub,
                              refinedEdgeDistInfo.locSubToGlobSub, refinedEdgeDistInfo.com)),
faceNormDistInfo(new DistInfo(refinedEdgeDistInfo.numLocThreads, refinedEdgeDistInfo.numLocSub, refinedEdgeDistInfo.numGlobSub,
                              refinedEdgeDistInfo.locSubToGlobSub, refinedEdgeDistInfo.com)),
    numLocSub(refinedNodeDistInfo.numLocSub), ownsData(true), connectivity(new Connectivity*[numLocSub]),
    edges(new EdgeSet*[numLocSub]),faces(new FaceSet*[numLocSub]), pNodeMapping(new DistVec<int>(refinedNodeDistInfo)), pEdgeMapping(new DistVec<int>(refinedEdgeDistInfo)),
    lineMap(new DistSVec<int,2>(refinedNodeDistInfo)),lineIDMap(new DistVec<int>(refinedNodeDistInfo)),lineLocIDMap(new DistVec<int>(refinedNodeDistInfo))
{
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    connectivity[iSub] = 0;
    edges[iSub] = 0;
    faces[iSub] = 0;
  }
  *lineIDMap = -1;
  agglomType = AgglomerationLocal;
  mgMethod = mgm;

  mesh_topology_threshold = 0.95;

  useVolumeWeightedAverage = true;

  numLines = new int[numLocSub];
  memset(numLines,0,sizeof(int)*numLocSub);
  lineids = new std::vector<int>[numLocSub];

  lineLengths = new int*[numLocSub];

  edgeNormals = NULL; 
 
  parent = mg;

  coarseNodeTopology = NULL;

  faceDistInfo = new DistInfo(refinedNodeDistInfo.numLocThreads, refinedNodeDistInfo.numLocSub, refinedEdgeDistInfo.numGlobSub,refinedNodeDistInfo.locSubToGlobSub, refinedNodeDistInfo.com);
  inletNodeDistInfo = new DistInfo(refinedNodeDistInfo.numLocThreads, refinedNodeDistInfo.numLocSub, refinedEdgeDistInfo.numGlobSub,refinedNodeDistInfo.locSubToGlobSub, refinedNodeDistInfo.com);
  agglomFaceDistInfo = new DistInfo(refinedNodeDistInfo.numLocThreads, refinedNodeDistInfo.numLocSub, refinedEdgeDistInfo.numGlobSub,refinedNodeDistInfo.locSubToGlobSub, refinedNodeDistInfo.com);

  edgeArea = NULL;

  nodeToNodeMaskILU = NULL;

  faceMapping = NULL;
  nodeVecPattern = NULL;
  nodeVecPatternEq1 = NULL;
  nodeVecPatternEq2 = NULL;

  matPattern2 = NULL;
}

struct agg_face {

  int node_id;
  int code;
  Vec3D normal;
  double area;
};

struct loc_edge {

  std::pair<int,int> ij;
  double area;
  Vec3D normal;
  bool owned;
};

// Read the multigrid level from a file
template<class Scalar>
MultiGridLevel<Scalar>::MultiGridLevel(MultiGridLevel* mg,int level_num,
                                       double ref_length,
				       Domain& domain,const char* fn_base,
				       int dim, int neq1, int neq2,
				       GeoSource& geoSource): domain(domain),
neq1(neq1), neq2(neq2) {


  connectivity = NULL;
  lineids = NULL;
  numLines = NULL;
  sharedNodes = NULL;
  agglomType = AgglomerationGlobal;
  DistInfo& dI = domain.getNodeDistInfo();

  parent = mg;
  
  faceMapping = NULL;
  pNodeMapping = new DistVec<int>(*parent->nodeDistInfo);

  int *neighCPUcount = new int[dI.numGlobSub];

  FILE* fin_base = fopen(fn_base,"rb");

  fread(neighCPUcount, sizeof(int), dI.numGlobSub, fin_base);

  subToSub = new Connectivity(dI.numGlobSub, neighCPUcount);
  
  for (int i = 0 ; i < dI.numGlobSub; ++i)
    fread((*subToSub)[i],sizeof(int), neighCPUcount[i],fin_base);

  fclose(fin_base);

  Communicator* com = domain.getCommunicator();
  levelSubDTopo = new SubDTopo(com->cpuNum(), subToSub, geoSource.getCpuToSub());

  nodeIdPattern = NULL;
  nodePosnPattern = NULL;
  offDiagMatPattern = NULL;
  edgeAreaPattern = NULL;
  edgeVecPattern = NULL;
  nodeVecPatternEq1 = nodeVecPatternEq2 = NULL;
  nodeIdPattern = new CommPattern<int>(levelSubDTopo, domain.getCommunicator(), CommPattern<int>::CopyOnSend);
  nodeVolPattern = new CommPattern<double>(levelSubDTopo, domain.getCommunicator(), CommPattern<double>::CopyOnSend);
  nodeVecPattern = new CommPattern<double>(levelSubDTopo, domain.getCommunicator(), CommPattern<double>::CopyOnSend);
  
  if (neq2 > 0) {
    nodeVecPatternEq1 = new CommPattern<double>(levelSubDTopo, domain.getCommunicator(), CommPattern<double>::CopyOnSend);
    nodeVecPatternEq2 = new CommPattern<double>(levelSubDTopo, domain.getCommunicator(), CommPattern<double>::CopyOnSend);
  }
  nodePosnPattern = new CommPattern<double>(levelSubDTopo, domain.getCommunicator(), CommPattern<double>::CopyOnSend);
  matPattern = new CommPattern<double>(levelSubDTopo, domain.getCommunicator(), CommPattern<double>::CopyOnSend);
  matPattern2 = new CommPattern<double>(levelSubDTopo, domain.getCommunicator(), CommPattern<double>::CopyOnSend);

  numLocSub = domain.getNumLocSub();
  
  nodeToNodeMaskILU = new Connectivity*[numLocSub];
  memset(nodeToNodeMaskILU, 0 , sizeof(Connectivity*)*numLocSub);
  
  nodeDistInfo = new DistInfo(dI.numLocThreads, dI.numLocSub, dI.numGlobSub,
			      dI.locSubToGlobSub, dI.com);

  DistInfo& edI = domain.getEdgeDistInfo();
  edgeDistInfo = new DistInfo(dI.numLocThreads, dI.numLocSub, dI.numGlobSub,
			      dI.locSubToGlobSub, dI.com);
  int s;

  inletNodeDistInfo = new DistInfo(dI.numLocThreads, dI.numLocSub, dI.numGlobSub,
			      dI.locSubToGlobSub, dI.com);
  faceDistInfo = new DistInfo(dI.numLocThreads, dI.numLocSub, dI.numGlobSub,
 	      dI.locSubToGlobSub, dI.com);
  agglomFaceDistInfo = new DistInfo(dI.numLocThreads, dI.numLocSub, dI.numGlobSub,
 	                           dI.locSubToGlobSub, dI.com);
  agglomeratedFaces = new AgglomeratedFaceSet*[numLocSub];
  edges = new EdgeSet*[numLocSub];

  numNodes = new int[numLocSub]; 

  int numTotalNodes = 0;

  mgMethod = MultiGridGeometric; 
 
  useVolumeWeightedAverage = true;
  mgSubdomains = new MultigridSubdomain[numLocSub];
  std::vector<loc_edge>* myEdgesArray = new std::vector<loc_edge>[numLocSub];
  nodeType = new int*[numLocSub];
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    SubDomain * S = domain.getSubDomain()[iSub];
    int glnum = S->getGlobSubNum();
    MultigridSubdomain& SD = mgSubdomains[iSub];
    
    char fn[256];
    sprintf(fn,"%s.%d",fn_base,glnum+1);
        
    FILE* fin = fopen(fn,"rb");

    int parent_nn;
    fread(&parent_nn,sizeof(int), 1, fin);
    assert(parent_nn == pNodeMapping->subSize(iSub));
    int* tmp_map = new int[parent_nn];
    fread(tmp_map, sizeof(int),parent_nn,fin);
    int* tmp_order = new int[parent_nn];
    fread(tmp_order, sizeof(int),parent_nn,fin);

    // We have to reorder the nodes based upon the actual numbering of 
    // the nodes used for this simulation
    std::map<int,int> reordering;
    if (level_num == 1) {

      for (int i = 0; i < dI.subSize(iSub); ++i) {

        reordering[S->getNodeMap()[i]] = i;
      }

    } else {

      MultigridSubdomain& SD2 = parent->mgSubdomains[iSub];
      for (int i = 0; i < SD2.ownedNodes.size(); ++i) {

        reordering[SD2.ownedNodes[i]] = i;
      }
      
      for (int i = 0; i < SD2.sharedNodes.size(); ++i) {

        reordering[SD2.sharedNodes[i]] = i+SD2.ownedNodes.size();
      }

    }

    for (int i = 0; i < parent_nn; ++i) {

      (*pNodeMapping)(iSub)[reordering[tmp_order[i]]] = tmp_map[i];

    }

    delete [] tmp_map;
    delete [] tmp_order;

    int rnk;
    MPI_Comm_rank(MPI_COMM_WORLD,&rnk);
    /*
    if (rnk == 7) {

      for (int i = 0; i < parent_nn; ++i) {

        if (pNodeMapping->v[i] == 1617679)
          std::cout << pNodeMapping->v[i] << " ";
      }
    }
    */

    fread(&s, sizeof(int), 1, fin);
    SD.ownedNodes.resize(s);
    fread(&SD.ownedNodes[0], sizeof(int), s, fin);

    numTotalNodes += s;

    fread(&s, sizeof(int), 1, fin);
    SD.sharedNodes.resize(s);
    fread(&SD.sharedNodes[0], sizeof(int), s, fin);

    SD.numNodes = SD.ownedNodes.size() + SD.sharedNodes.size();

    numNodes[iSub] = SD.numNodes;

    SD.ctrlVolCount = new int[SD.numNodes];
    SD.locToGlobMap = new int[SD.numNodes];

    SD.nodeTopology = new MultigridSubdomain::Topology[SD.numNodes];

    nodeType[iSub] = new int[SD.numNodes];

    fread(SD.ctrlVolCount , sizeof(int),SD.numNodes,fin);
    fread(SD.locToGlobMap , sizeof(int),SD.numNodes,fin);
    fread(nodeType[iSub], sizeof(int),SD.numNodes,fin);
    fread(SD.nodeTopology, sizeof(MultigridSubdomain::Topology),SD.numNodes,fin);

    for (int i = 0; i < parent_nn; ++i) {

      assert((*pNodeMapping)(iSub)[i] >= 0 && (*pNodeMapping)(iSub)[i] < SD.numNodes);
      if (SD.locToGlobMap[(*pNodeMapping)(iSub)[i] ] == 1617679)
        std::cout << "map " << level_num-1 << "->" << level_num << ": " << i  << "->" << (*pNodeMapping)(iSub)[i] << std::endl;
    }

    nodeDistInfo->setLen(iSub, SD.numNodes);
    
    double* tmp = new double[SD.numNodes*4];
    fread(tmp, sizeof(double),SD.numNodes*4,fin);
    
    SD.nodeVolumes = new double[SD.numNodes];
    SD.X = new Vec3D[SD.numNodes];
    
    for (int i = 0; i < SD.numNodes; ++i) {

      SD.nodeVolumes[i] = tmp[i*4] / (ref_length*ref_length*ref_length);
      SD.X[i] = Vec3D(tmp[i*4+1],tmp[i*4+2],tmp[i*4+3])/(ref_length);
    }

    delete [] tmp;    
    
    // Read in the faces
    std::vector<agg_face> myFacesOwned, myFacesShared;

    fread(&s, sizeof(int), 1, fin);
    myFacesOwned.resize(s);

    Vec3D tot = 0.0;
    for (int i = 0; i < s; ++i) {

      int glSub = domain.getSubDomain()[iSub]->getGlobSubNum();
      fread(&myFacesOwned[i].node_id,sizeof(int),1, fin);
      fread(&myFacesOwned[i].code,sizeof(int),1, fin);
      fread(&myFacesOwned[i].area,sizeof(double),1, fin);
      fread(&myFacesOwned[i].normal,sizeof(Vec3D),1, fin);
      myFacesOwned[i].normal /= (ref_length*ref_length);
      myFacesOwned[i].area /= (ref_length*ref_length);
      if (SD.locToGlobMap[myFacesOwned[i].node_id] == 1617679) {
        std::cout << "Owned face (" << glSub << ", " << level_num << "), code = " << myFacesOwned[i].code << ": ";
        myFacesOwned[i].normal.print();
	tot += myFacesOwned[i].normal;
      }
    }
  
    fread(&s, sizeof(int), 1, fin);
    myFacesShared.resize(s);
    for (int i = 0; i < s; ++i) {

      int glSub = domain.getSubDomain()[iSub]->getGlobSubNum();
      fread(&myFacesShared[i].node_id,sizeof(int),1, fin);
      fread(&myFacesShared[i].code,sizeof(int),1, fin);
      fread(&myFacesShared[i].area,sizeof(double),1, fin);
      fread(&myFacesShared[i].normal,sizeof(Vec3D),1, fin);
      myFacesShared[i].normal /= (ref_length*ref_length);
      myFacesShared[i].area /= (ref_length*ref_length);
      if (SD.locToGlobMap[myFacesShared[i].node_id] == 1617679) {
        std::cout << "Shared face (" << glSub << ", " << level_num << "), code = " << myFacesShared[i].code << ": ";
        myFacesShared[i].normal.print();
      }
    }

    agglomeratedFaces[iSub] = new AgglomeratedFaceSet(myFacesOwned.size() + 
						      myFacesShared.size());
    
    
    for (int i = 0; i < myFacesOwned.size(); ++i) {

      AgglomeratedFace A(myFacesOwned[i].node_id, myFacesOwned[i].code);
      A.getArea() = myFacesOwned[i].area;
      A.getNormal() = myFacesOwned[i].normal;
      (*agglomeratedFaces[iSub])[i] = A;
    }
    
    for (int i = 0; i < myFacesShared.size(); ++i) {

      AgglomeratedFace A(myFacesShared[i].node_id, myFacesShared[i].code);
      A.getArea() = myFacesShared[i].area;
      A.getNormal() = myFacesShared[i].normal;
      A.getMasterFlag() = false;
      (*agglomeratedFaces[iSub])[myFacesOwned.size()+i] = A;
    }

    edges[iSub] = new EdgeSet;
    std::vector<loc_edge>& myEdges = myEdgesArray[iSub];
    fread(&s, sizeof(int), 1, fin);
    for (int i = 0; i < s; ++i) {

      loc_edge E;
      fread(&E.ij, sizeof(E.ij),1,fin);
      fread(&E.area, sizeof(E.area),1,fin);
      fread(&E.normal, sizeof(E.normal),1,fin);
      E.normal /= (ref_length*ref_length);
      E.area /= (ref_length*ref_length);
      E.owned = true;
      myEdges.push_back(E);
      edges[iSub]->find(E.ij.first,E.ij.second);

      if (SD.locToGlobMap[E.ij.first] == 1617679 || SD.locToGlobMap[E.ij.second] == 1617679) {
        std::cout << "Owned edge (" << SD.locToGlobMap[E.ij.first] << " , "
		  <<  SD.locToGlobMap[E.ij.second] << "): ";
	if (SD.locToGlobMap[E.ij.first] == 1617679)
	  tot += E.normal;
	else
	  tot -= E.normal;
	E.normal.print();
      }
    }
      
    fread(&s, sizeof(int), 1, fin);
    for (int i = 0; i < s; ++i) {

      loc_edge E;
      fread(&E.ij, sizeof(E.ij),1,fin);
      fread(&E.area, sizeof(E.area),1,fin);
      fread(&E.normal, sizeof(E.normal),1,fin);
      E.normal /= (ref_length*ref_length);
      E.area /= (ref_length*ref_length);
      E.owned = false;
      myEdges.push_back(E);
      edges[iSub]->find(E.ij.first,E.ij.second);
      if (SD.locToGlobMap[E.ij.first] == 1617679 || SD.locToGlobMap[E.ij.second] == 1617679) {
        std::cout << "Shared edge (" << SD.locToGlobMap[E.ij.first] << " , "
		  <<  SD.locToGlobMap[E.ij.second] << "): ";
	E.normal.print();
      }
    }  

    //std::cout << "ref_length = " << ref_length << " tot = ";
    //tot.print();

    edgeDistInfo->setLen(iSub, myEdges.size());

    Vec<int> newNum(edges[iSub]->size());
    edges[iSub]->createPointers(newNum);
    edges[iSub]->setMasterFlag(new bool[edges[iSub]->size()]);

    for (int i = 0; i < myEdges.size(); ++i) {

      int l = newNum[i];
      edges[iSub]->getMasterFlag()[l] = myEdges[i].owned;
    }

    fread(&s, sizeof(int), 1, fin);

    int t,u;

    for (int i = 0; i < s; ++i) {

      fread(&t, sizeof(int), 1, fin);
      NeighborDomain* N = new NeighborDomain;
      
      fread(&N->id,sizeof(int),1,fin);
      
      fread(&u, sizeof(int), 1, fin);
      
      N->sharedNodes.resize(u);
      fread(&N->sharedNodes[0], sizeof(std::pair<int,int>),u, fin);
      mgSubdomains[iSub].neighbors[t] = N;
    }
    
    inletNodeDistInfo->setLen(iSub, 0);
    faceDistInfo->setLen(iSub,myFacesOwned.size()+myFacesShared.size());
    agglomFaceDistInfo->setLen(iSub,myFacesOwned.size()+myFacesShared.size());
  }

  com->globalSum(1,&numTotalNodes);
  com->fprintf(stdout, "Number of nodes in multigrid level %d = %d\n",level_num, numTotalNodes);
  
  nodeDistInfo->finalize(true);
  edgeDistInfo->finalize(true);
  inletNodeDistInfo->finalize(true);
  faceDistInfo->finalize(true);
  agglomFaceDistInfo->finalize(true);

  ctrlVolCount = new DistVec<double>(*nodeDistInfo);  

  fv_comp_tag = new DistSVec<int,2>(*nodeDistInfo);
  *fv_comp_tag = 0; 
  
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    memset(nodeDistInfo->getMasterFlag(iSub), 0, sizeof(bool)*nodeDistInfo->subSize(iSub));
    for (int i = 0; i < mgSubdomains[iSub].ownedNodes.size(); ++i)
      nodeDistInfo->getMasterFlag(iSub)[i] = true;

    for (int i = 0; i < nodeDistInfo->subSize(iSub); ++i)
      (*ctrlVolCount)(iSub)[i] = mgSubdomains[iSub].ctrlVolCount[i];
  }
  
  
  myGeoState = new DistGeoState(parent->myGeoState->getGeoData(),&domain,
                                *nodeDistInfo,*edgeDistInfo);

  ctrlVol = &myGeoState->getCtrlVol();//new DistVec<double>(*nodeDistInfo);
  
  Xn = &myGeoState->getXn();//new DistSVec<double,3>(*nodeDistInfo); 
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    MultigridSubdomain& SD = mgSubdomains[iSub];
    memcpy(&(*ctrlVol)(iSub)[0], SD.nodeVolumes,sizeof(double)*SD.numNodes);
    memcpy(reinterpret_cast<double*>((*Xn)(iSub).data()), SD.X,sizeof(Vec3D)*SD.numNodes);
    
  }

  edgeArea = new DistVec<double>(*edgeDistInfo);
  *edgeArea = 0.0;

  edgeNormals = &myGeoState->getEdgeNormal();//new DistVec<Vec3D>(*edgeDistInfo);  
  *edgeNormals = Vec3D(0.0);
  
  myGeoState->getEdgeNormalVel() = 0.0;
  myGeoState->getFaceNorVel() = 0.0;

  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    int glSub = domain.getSubDomain()[iSub]->getGlobSubNum();
    std::vector<loc_edge>& myEdges = myEdgesArray[iSub];
    MultigridSubdomain* D = &mgSubdomains[iSub];
    for (int i = 0; i < myEdges.size(); ++i) {

      int l = edges[iSub]->find(myEdges[i].ij.first, myEdges[i].ij.second);

      /*if (D->locToGlobMap[myEdges[i].ij.first] == 1617679 ||
          D->locToGlobMap[myEdges[i].ij.second] == 1617679) {

        std::cout << glSub << ", " << level_num << ": " << D->locToGlobMap[myEdges[i].ij.first] << " " <<
	  D->locToGlobMap[myEdges[i].ij.second] << " " << myEdges[i].owned << " "; 
        myEdges[i].normal.print();
      }
      */
      if (myEdges[i].normal.norm() < 1e-16) {

	//std::cout << "Error: normal = 0!, (" << level_num << ") " << std::endl;
        //std::cout << glSub << ", " << level_num << ": " << D->locToGlobMap[myEdges[i].ij.first] << " " <<
	//  D->locToGlobMap[myEdges[i].ij.second] << " " << myEdges[i].owned << " "; 
        //myEdges[i].normal.print();
      }
      (*edgeNormals)(iSub)[l] = myEdges[i].normal;
      (*edgeArea)(iSub)[l] = myEdges[i].area;
    }
  }
  

  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    int glSub = domain.getSubDomain()[iSub]->getGlobSubNum();
    for (std::map<int, NeighborDomain*>::iterator itr =	
	   mgSubdomains[iSub].neighbors.begin();
	 itr != mgSubdomains[iSub].neighbors.end(); ++itr) {
      int neighb = itr->first;
      nodeVecPattern->setLen(sndChannel(glSub,neighb),dim*itr->second->sharedNodes.size());
      nodeVecPattern->setLen(rcvChannel(glSub,neighb),dim*itr->second->sharedNodes.size());
      
      nodeVolPattern->setLen(sndChannel(glSub,neighb),itr->second->sharedNodes.size());
      nodeVolPattern->setLen(rcvChannel(glSub,neighb),itr->second->sharedNodes.size());

      nodeIdPattern->setLen(sndChannel(glSub,neighb),itr->second->sharedNodes.size());
      nodeIdPattern->setLen(rcvChannel(glSub,neighb),itr->second->sharedNodes.size());
      if (neq2 > 0) {
	nodeVecPatternEq1->setLen(sndChannel(glSub,neighb),neq1*itr->second->sharedNodes.size());
	nodeVecPatternEq1->setLen(rcvChannel(glSub,neighb),neq1*itr->second->sharedNodes.size());

	nodeVecPatternEq2->setLen(sndChannel(glSub,neighb),neq2*itr->second->sharedNodes.size());
	nodeVecPatternEq2->setLen(rcvChannel(glSub,neighb),neq2*itr->second->sharedNodes.size());
      }
      matPattern->setLen(sndChannel(glSub,neighb),neq1*neq1*itr->second->sharedNodes.size());
      matPattern->setLen(rcvChannel(glSub,neighb),neq1*neq1*itr->second->sharedNodes.size());

      if (neq2 > 0) {
	matPattern2->setLen(sndChannel(glSub,neighb),neq2*neq2*itr->second->sharedNodes.size());
	matPattern2->setLen(rcvChannel(glSub,neighb),neq2*neq2*itr->second->sharedNodes.size());
      }
    }
  }

  nodeIdPattern->finalize();
  nodeVolPattern->finalize();
  nodeVecPattern->finalize();
  if (neq2 > 0) {
    nodeVecPatternEq1->finalize();
    nodeVecPatternEq2->finalize();
    matPattern2->finalize();
  }
  //nodePosnPattern->finalize();
  matPattern->finalize();
}

/*
template<class Scalar>
void MultiGridLevel<Scalar>::MultiGridLevel(MultiGridLevel& level1, MultiGridLevel& level2) :
my_dim(level1.my_dim), neq1(level1.neq1), neq2(level1.neq2),
domain(level1.domain)  
  {

  nodeIdPattern = level2.nodeIdPattern;
  nodeVolPattern = level2.nodeValPattern;
  nodeVecPattern = level2.nodeVecPattern;
  nodeVecPatternEq1 = level2.nodeVecPatternEq1;
  nodeVecPatternEq2 = level2.nodeVecPatternEq2;
  nodeNormalsPattern = level2.nodeNormalsPattern;
  nodePosnPattern = level2.nodePosnPattern;
  matPattern = level2.matPattern;
  offDiagMatPattern = level2.offDiagMatPattern;
  edgeAreaPattern = level2.edgeAreaPattern;
  edgeVecPattern = level2.edgeVecPattern;

  sharedNodes = level2.sharedNodes;
  nodeToNodeMaskILU = level2.nodeToNodeMaskILU;

  myGeoState = level2.myGeoState;
  myIoData = level2.myIoData;

  useVolumeWeightedAverage = level2.useVolumeWeightedAverage;

  nodeDistInfo = level2.nodeDistInfo;
  edgeDistInfo = level2.edgeDistInfo;
  faceDistInfo = level2.faceDistInfo;
  agglomFaceDistInfo = level2.agglomFaceDistInfo;
  faceNormDistInfo = level2.faceNormDistInfo;
  inletNodeDistInfo = level2.inletNodeDistInfo;

  numLocSub = level2.numLocSub;
  connectivity = level2.connectivity;
  edges = level2.edges;
  faces = level2.faces;
  agglomeratedFaces = level2.agglomeratedFaces;
  sharedEdges = level2.sharedEdges 
}
*/
template<class Scalar>
void MultiGridLevel<Scalar>::WriteTopFile(const std::string& fileName) {

  std::ofstream outfile(fileName.c_str());
  int nodenum = 1;
  int* node_starts = new int[numLocSub];
  outfile << "Nodes FluidNodes\n";
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    node_starts[iSub] = nodenum;
    SVec<double,3>& myXn = (*Xn)(iSub);
    for (int i = 0; i < myXn.size(); ++i,++nodenum)
      outfile << nodenum << " " << myXn[i][0] << " " << myXn[i][1] << " " << myXn[i][2] << "\n";
  }

  int elem_num = 1;
/*  outfile << "Elements FluidMesh using FluidNodes\n";
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    ElemSet& elem_set = *elems[iSub];
    for (int i = 0; i < elem_set.size(); ++i,++elem_num) {

      outfile << elem_num << " 5 ";
      for (int j = 0; j < 4; ++j)
        outfile << node_starts[iSub] + elem_set[i][j] << " ";
      outfile << "\n";
    }

  }
*/
  outfile << "Elements StickMovingSurface_0 using FluidNodes\n";
  elem_num = 1;
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    FaceSet& elem_set = *faces[iSub];
    for (int i = 0; i < elem_set.size(); ++i,++elem_num) {

      outfile << elem_num << " 4 ";
      for (int j = 0; j < 3; ++j)
        outfile << node_starts[iSub] + elem_set[i][j] << " ";
      outfile << "\n";
    }
  }
  

  outfile.close();

}

template<class Scalar>
template<int dim>
void MultiGridLevel<Scalar>::writeXpostFile(const std::string& fileName,
                                            DistSVec<Scalar,dim>& val,int id) {

  std::ofstream outfile(fileName.c_str());
  outfile << "Scalar Pressure under load for FluidNodes\n";
  int nodes_inc_overlap = 0;
#pragma omp parallel for 
  for(int iSub = 0; iSub < numLocSub; ++iSub) { 
    nodes_inc_overlap += val.subSize(iSub);
  }
  outfile << nodes_inc_overlap << "\n0\n";
#pragma omp parallel for 
  for(int iSub = 0; iSub < numLocSub; ++iSub) { 
    for (int i = 0; i < val.subSize(iSub); ++i) {

      outfile << val(iSub)[i][id] << "\n";
    }
  }
}

template<class Scalar>
void MultiGridLevel<Scalar>::writeXpostFile(const std::string& fileName,
                                            DistVec<Scalar>& val) {

  std::ofstream outfile(fileName.c_str());
  outfile << "Scalar Pressure under load for FluidNodes\n";
  int nodes_inc_overlap = 0;
#pragma omp parallel for 
  for(int iSub = 0; iSub < numLocSub; ++iSub) { 
    nodes_inc_overlap += val.subSize(iSub);
  }
  outfile << nodes_inc_overlap << "\n0\n";
#pragma omp parallel for 
  for(int iSub = 0; iSub < numLocSub; ++iSub) { 
    for (int i = 0; i < val.subSize(iSub); ++i) {

      outfile << val(iSub)[i] << "\n";
    }
  }
}

//------------------------------------------------------------------------------

template<class Scalar>
MultiGridLevel<Scalar>::~MultiGridLevel()
{
  if(ownsData){
    if (nodeIdPattern)
      delete nodeIdPattern;
    delete nodeVolPattern;
    if (nodeVecPattern)
      delete nodeVecPattern;
    if (nodePosnPattern)
      delete nodePosnPattern;
    if (nodeVecPatternEq1) {
      delete nodeVecPatternEq1;
      delete nodeVecPatternEq2;
    }
    delete matPattern;
    if (matPattern2)
      delete matPattern2;
    if (offDiagMatPattern)
      delete offDiagMatPattern;
    if (edgeAreaPattern)
      delete edgeAreaPattern;
    if (edgeVecPattern)
      delete edgeVecPattern;
    
  }
  
  if (faceMapping)
    delete faceMapping;
  delete nodeDistInfo;
  delete edgeDistInfo;

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub){
    if(ownsData) {
      if (connectivity)
        delete connectivity[iSub];
      delete edges[iSub];
    }
    if (nodeToNodeMaskILU[iSub])
      delete  nodeToNodeMaskILU[iSub];
  }
  if (sharedNodes)
    delete []sharedNodes;
  if (connectivity)
    delete []connectivity;
  delete []edges;
 
  if (lineids) { 
    delete [] lineids;
    delete [] numLines;
  }

  delete [] numNodes;
  if (nodeToNodeMaskILU)
    delete [] nodeToNodeMaskILU;
  
}

//------------------------------------------------------------------------------

template<class Scalar>
void MultiGridLevel<Scalar>::mergeFinerInto(MultiGridLevel& finer) {

  DistVec<int>* pn = new DistVec<int>(*finer.pNodeMapping);
  DistVec<int>* en = new DistVec<int>(*finer.pEdgeMapping);

  Connectivity** nToN = finer.parent->getConnectivity(); 

  delete lineIDMap;
  delete lineMap;

  lineLocIDMap = new DistVec<int>(*finer.pNodeMapping);
  lineIDMap = new DistVec<int>(*finer.pNodeMapping);
 
  lineMap = new DistSVec<int,2>(*finer.parent->nodeDistInfo);

  *lineIDMap = -1;

  *lineMap = -1;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    delete [] lineLengths[iSub];
    lineLengths[iSub] = new int[(*pn)(iSub).size()];
    lineids[iSub].clear();
    numLines[iSub] = 0;

    std::vector<int>* agglom = new std::vector<int>[nodeDistInfo->subSize(iSub)];
 
    for (int i = 0; i < (*pn)(iSub).size(); ++i) {
   
      (*pn)(iSub)[i] = (*pNodeMapping)(iSub)[(*pn)(iSub)[i]];
      if (nodeDistInfo->getMasterFlag(iSub)[i])
        agglom[(*pn)(iSub)[i]].push_back(i);
    }
    for (int i = 0; i < (*en)(iSub).size(); ++i) {
  
      if ((*en)(iSub)[i] >= 0) 
        (*en)(iSub)[i] = (*pEdgeMapping)(iSub)[(*en)(iSub)[i]];
    }

    for (int i = 0; i < nodeDistInfo->subSize(iSub); ++i) {

      bool isLine = true;

      int j = 0;
      for(std::vector<int>::const_iterator iter = agglom[i].begin(); iter != agglom[i].end(); iter++,++j) {
        const int current_node = *iter;
        //nodeMapping(iSub)[current_node] = agglom_id;
        for(int n = 0; n < nToN[iSub]->num(current_node) && isLine; ++n) {
          const int neighbor = (*nToN[iSub])[current_node][n];
          if (current_node == neighbor) continue;
  
          for (int k = 0; k < agglom[i].size(); ++k) {
            if ((k > j+1 || k < j-1) && neighbor == agglom[i][k]) { 
              isLine = false;
              break;
            }
          }
        }
      }
 
      if (isLine && agglom[i].size() > 1) {
        int k;
        for (k = 0; k < agglom[i].size(); ++k) {
  
          (*lineLocIDMap)(iSub)[agglom[i][k]] = k;         
 
          (*lineIDMap)(iSub)[agglom[i][k]] = numLines[iSub];
          if (k > 0)
            (*lineMap)(iSub)[agglom[i][k]][0] = agglom[i][k-1];
          else
            (*lineMap)(iSub)[agglom[i][k]][0] = -1;
          if (k < agglom[i].size()-1)
            (*lineMap)(iSub)[agglom[i][k]][1] = agglom[i][k+1];  
          else
            (*lineMap)(iSub)[agglom[i][k]][1] = -1;
          lineids[iSub].push_back(agglom[i][k]);
        }
        lineLengths[iSub][numLines[iSub]] = k;

        for (; k < 8; ++k)
          lineids[iSub].push_back(-1);

        ++numLines[iSub];
      } else {
        for (int k = 0; k < agglom[i].size(); ++k) {
          (*lineMap)(iSub)[agglom[i][k]][0] = (*lineMap)(iSub)[agglom[i][k]][1] = -1;
        }
      }      
    }
    delete [] agglom;
  }

  delete pEdgeMapping;
  delete pNodeMapping;

  parent = finer.parent;

  pNodeMapping = pn;
  pEdgeMapping = en;

  delete nodeTopology;
  delete topologyNormal;

  computeNodeNormalClasses();
}

template<class Scalar>
void MultiGridLevel<Scalar>::copyRefinedState(const DistInfo& refinedNodeDistInfo, const DistInfo& refinedEdgeDistInfo,const DistInfo& refinedInletNodeDistInfo,const DistInfo& refinedFaceDistInfo,
                                              DistGeoState& refinedGeoState, Domain& domain,
                                              int dim,int lneq1,int lneq2)
{
  ownsData = false;
  nodeIdPattern = domain.getLevelPat();
  nodeVolPattern = domain.getVolPat();
  nodeVecPattern = domain.getVecPat();
  nodePosnPattern = domain.getVec3DPat();
  sharedNodes = new Connectivity*[numLocSub];

  sharedEdges = new EdgeDef**[numLocSub];
  numSharedEdges = new int*[numLocSub];
  numNodes = new int[numLocSub]; 
  nodeToNodeMaskILU = new Connectivity*[numLocSub];
  memset(nodeToNodeMaskILU, 0 , sizeof(Connectivity*)*numLocSub);

  my_dim = dim;
  neq1 = lneq1;
  neq2 = lneq2;

  if (neq2 > 0) {
    nodeVecPatternEq1 = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);
    nodeVecPatternEq2 = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);
  }
 
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    connectivity[iSub] = domain.getSubDomain()[iSub]->getNodeToNode();
    sharedNodes[iSub] = domain.getSubDomain()[iSub]->getSharedNodes();
    edges[iSub] = &domain.getSubDomain()[iSub]->getEdges();
    faces[iSub] = &domain.getSubDomain()[iSub]->getFaces();
    nodeDistInfo->setLen(iSub, refinedNodeDistInfo.subSize(iSub));
    edgeDistInfo->setLen(iSub, refinedEdgeDistInfo.subSize(iSub));
    faceDistInfo->setLen(iSub, refinedFaceDistInfo.subSize(iSub));
    inletNodeDistInfo->setLen(iSub, refinedInletNodeDistInfo.subSize(iSub));

    if (neq2 > 0) {

      domain.getSubDomain()[iSub]->setComLenNodes(neq1, *nodeVecPatternEq1); 
      domain.getSubDomain()[iSub]->setComLenNodes(neq2, *nodeVecPatternEq2); 
    }

    SubDomain& subD(*domain.getSubDomain()[iSub]);
    numSharedEdges[iSub] = new int[subD.getNumNeighb()];
    sharedEdges[iSub] = new EdgeDef*[subD.getNumNeighb()];
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
 
      numSharedEdges[iSub][jSub] = subD.getNumSharedEdges()[jSub];
      sharedEdges[iSub][jSub] = new EdgeDef[numSharedEdges[iSub][jSub]];
      memcpy(sharedEdges[iSub][jSub], subD.getSharedEdges()[jSub], 
             sizeof(EdgeDef)*numSharedEdges[iSub][jSub]);
    }

    numNodes[iSub] = subD.numNodes();
  }
  nodeDistInfo->finalize(true);
  edgeDistInfo->finalize(true);
  faceDistInfo->finalize(true);
  inletNodeDistInfo->finalize(true);

  if (neq2 > 0) { 

    nodeVecPatternEq1->finalize();
    nodeVecPatternEq2->finalize();
  }

  edgeNormals = &refinedGeoState.getEdgeNormal();
 
  agglomeratedFaces = NULL;

  myGeoState = &refinedGeoState;

  //ctrlVol = &refinedGeoState.getCtrlVol();
  Xn = &refinedGeoState.getXn();
  ctrlVol = new DistVec<double>(refinedNodeDistInfo);
  
  *ctrlVol = refinedGeoState.getCtrlVol();

  ctrlVolCount = new DistVec<double>(refinedNodeDistInfo);  

  *ctrlVolCount = 1.0;

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub)
    for(int i = 0; i < refinedNodeDistInfo.subSize(iSub); ++i) {
      nodeDistInfo->getMasterFlag(iSub)[i] = refinedNodeDistInfo.getMasterFlag(iSub)[i];
      nodeDistInfo->getInvWeight(iSub)[i] = refinedNodeDistInfo.getInvWeight(iSub)[i];
    }

  total_mesh_volume = myGeoState->getCtrlVol().sum();

}

//------------------------------------------------------------------------------

namespace { // MPI helper functions
template<class T,class OperType>
void assemble(Domain & domain, CommPattern<T> & commPat, Connectivity ** sharedNodes, DistVec<T> & data, const OperType & op) {
#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      SubRecInfo<T> sInfo = commPat.getSendBuffer(subD.getSndChannel()[jSub]);
      T * buffer = reinterpret_cast<T *>(sInfo.data);

      for (int iNode = 0; iNode < sharedNodes[iSub]->num(jSub); ++iNode)
        buffer[iNode] = data(iSub)[ (*sharedNodes[iSub])[jSub][iNode] ];
    }
  }

  commPat.exchange();

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      SubRecInfo<T> sInfo = commPat.recData(subD.getRcvChannel()[jSub]);
      T * buffer = reinterpret_cast<T *>(sInfo.data);

      for (int iNode = 0; iNode < sharedNodes[iSub]->num(jSub); ++iNode)
        data(iSub)[ (*sharedNodes[iSub])[jSub][iNode] ] = OperType::apply(data(iSub)[ (*sharedNodes[iSub])[jSub][iNode] ], buffer[iNode]);
    }
  }
}

template<class T,class OperType,int dim>
void assemble(Domain & domain, CommPattern<T> & commPat, Connectivity ** sharedNodes, DistSVec<T,dim> & data, const OperType & op,std::map<int,std::map<int,int> >* locToGlobMap = NULL) {
#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      SubRecInfo<T> sInfo = commPat.getSendBuffer(subD.getSndChannel()[jSub]);
      T (*buffer)[dim] = reinterpret_cast<T (*)[dim]>(sInfo.data);

      for (int iNode = 0; iNode < sharedNodes[iSub]->num(jSub); ++iNode)
        for (int j = 0; j < dim; ++j)
          buffer[iNode][j] = data(iSub)[ (*sharedNodes[iSub])[jSub][iNode] ][j];
    }
  }

  commPat.exchange();

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      SubRecInfo<T> sInfo = commPat.recData(subD.getRcvChannel()[jSub]);
      T (*buffer)[dim] = reinterpret_cast<T (*)[dim]>(sInfo.data);

      for (int iNode = 0; iNode < sharedNodes[iSub]->num(jSub); ++iNode) {
        /*if (locToGlobMap && (*locToGlobMap)[iSub][(*sharedNodes[iSub])[jSub][iNode]] == 3911)
          std::cout << iSub << " " << data(iSub)[ (*sharedNodes[iSub])[jSub][iNode] ][0] <<
                    " " << buffer[iNode][0] << std::endl;
        */
        for (int j = 0; j < dim; ++j)
          data(iSub)[ (*sharedNodes[iSub])[jSub][iNode] ][j] = OperType::apply(data(iSub)[ (*sharedNodes[iSub])[jSub][iNode]][j], buffer[iNode][j]);
      }
    }
  }
}

template<class T,int dim>
void assemble(Domain & domain, CommPattern<T> & commPatDiag, Connectivity ** sharedNodes, DistMat<T,dim> & A) {

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {

      SubRecInfo<T> sInfo = commPatDiag.getSendBuffer(subD.getSndChannel()[jSub]);
      T (*buffer)[dim*dim] = reinterpret_cast<T (*)[dim*dim]>(sInfo.data);

      for (int iNode = 0; iNode < sharedNodes[iSub]->num(jSub); ++iNode) {
        T * a = A(iSub).getElem_ii((*sharedNodes[iSub])[jSub][iNode]);
        for (int j=0; j<dim*dim; ++j) buffer[iNode][j] = a[j];
      }
    }
  }

  commPatDiag.exchange();

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      SubRecInfo<T> sInfo = commPatDiag.recData(subD.getRcvChannel()[jSub]);
      T (*buffer)[dim*dim] = reinterpret_cast<T (*)[dim*dim]>(sInfo.data);

      for (int iNode = 0; iNode < sharedNodes[iSub]->num(jSub); ++iNode) {
        T * a = A(iSub).getElem_ii((*sharedNodes[iSub])[jSub][iNode]);
        for (int j = 0; j < dim*dim; ++j) a[j] += buffer[iNode][j];
      }
    }
  }

}

template<class T,int dim>
void assemble(Domain & domain, CommPattern<T> & commPatOffDiag, int **numSharedEdges, EdgeDef ***sharedEdges, DistMat<T,dim> & A) {

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {

      SubRecInfo<T> sInfo = commPatOffDiag.getSendBuffer(subD.getSndChannel()[jSub]);
      T (*buffer)[2][dim*dim] = reinterpret_cast<T (*)[2][dim*dim]>(sInfo.data);

      for (int iEdge = 0; iEdge < numSharedEdges[iSub][jSub]/*sharedNodes[iSub]->num(jSub)*/; ++iEdge) {
       
        int edgeNum = sharedEdges[iSub][jSub][iEdge].edgeNum;

        T *aij, *aji;

        if (sharedEdges[iSub][jSub][iEdge].sign > 0) {
  	  aij = A(iSub).getElem_ij(edgeNum);
	  aji = A(iSub).getElem_ji(edgeNum);
        }
        else {
 	  aij = A(iSub).getElem_ji(edgeNum);
	  aji = A(iSub).getElem_ij(edgeNum);
        }

        if (aij && aji) {
  	  for (int k=0; k<dim*dim; ++k) {
	    buffer[iEdge][0][k] = aij[k];
	    buffer[iEdge][1][k] = aji[k];
	  }
        }
        else {
 	  for (int k=0; k<dim*dim; ++k) {
	    buffer[iEdge][0][k] = 0.0;
	    buffer[iEdge][1][k] = 0.0;
	  }
        }      
      }
    }
  }

  commPatOffDiag.exchange();

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      SubRecInfo<T> sInfo = commPatOffDiag.recData(subD.getRcvChannel()[jSub]);
      T (*buffer)[2][dim*dim] = reinterpret_cast<T (*)[2][dim*dim]>(sInfo.data);
      for (int iEdge = 0; iEdge < numSharedEdges[iSub][jSub]; ++iEdge) {

        int edgeNum = sharedEdges[iSub][jSub][iEdge].edgeNum;

        T *aij, *aji;

        if (sharedEdges[iSub][jSub][iEdge].sign > 0) {
  	  aij = A(iSub).getElem_ij(edgeNum);
	  aji = A(iSub).getElem_ji(edgeNum);
        }
        else {
	  aij = A(iSub).getElem_ji(edgeNum);
	  aji = A(iSub).getElem_ij(edgeNum);
        }

        if (aij && aji) {
	  for (int k=0; k<dim*dim; ++k) {
	    aij[k] += buffer[iEdge][0][k];
	    aji[k] += buffer[iEdge][1][k];
 	  }
        }

      }

    }
  }
}

template<class T>
void assemble(Domain & domain, CommPattern<T> & edgePat, int **numSharedEdges, EdgeDef ***sharedEdges, DistVec<T> & V) {

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {

      SubRecInfo<T> sInfo = edgePat.getSendBuffer(subD.getSndChannel()[jSub]);
      T *buffer = reinterpret_cast<T *>(sInfo.data);

      for (int iEdge = 0; iEdge < numSharedEdges[iSub][jSub]; ++iEdge) {
       
        int edgeNum = sharedEdges[iSub][jSub][iEdge].edgeNum;
        buffer[iEdge] = V(iSub)[edgeNum];
      }
    }
  }

  edgePat.exchange();

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      SubRecInfo<T> sInfo = edgePat.recData(subD.getRcvChannel()[jSub]);
      T *buffer = reinterpret_cast<T*>(sInfo.data);
      for (int iEdge = 0; iEdge < numSharedEdges[iSub][jSub]; ++iEdge) {

        int edgeNum = sharedEdges[iSub][jSub][iEdge].edgeNum;
        V(iSub)[edgeNum] += buffer[iEdge];
      }

    }
  }
}

void assemble(Domain & domain, CommPattern<double> & edgePat, int **numSharedEdges, EdgeDef ***sharedEdges, DistVec<Vec3D> & V) {

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {

      SubRecInfo<double> sInfo = edgePat.getSendBuffer(subD.getSndChannel()[jSub]);
      double (*buffer)[3] = reinterpret_cast<double (*)[3]>(sInfo.data);

      for (int iEdge = 0; iEdge < numSharedEdges[iSub][jSub]; ++iEdge) {
       
	int edgeNum = sharedEdges[iSub][jSub][iEdge].edgeNum;
        for (int k = 0; k < 3; ++k)
  	  buffer[iEdge][k] = V(iSub)[edgeNum][k] * sharedEdges[iSub][jSub][iEdge].sign;
      }
    }
  }

  edgePat.exchange();

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      SubRecInfo<double> sInfo = edgePat.recData(subD.getRcvChannel()[jSub]);
      double (*buffer)[3] = reinterpret_cast<double (*)[3]>(sInfo.data);
      for (int iEdge = 0; iEdge < numSharedEdges[iSub][jSub]; ++iEdge) {

	int edgeNum = sharedEdges[iSub][jSub][iEdge].edgeNum;
        for (int k = 0; k < 3; ++k)
  	  V(iSub)[edgeNum][k] += buffer[iEdge][k] * sharedEdges[iSub][jSub][iEdge].sign;
      }

    }
  }
}
void assemble(Domain & domain, CommPattern<double> & commPat, Connectivity ** sharedNodes, std::list<Vec3D>** data, int max_node_adj) {
#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      SubRecInfo<double> sInfo = commPat.getSendBuffer(subD.getSndChannel()[jSub]);
      double * buffer = reinterpret_cast<double *>(sInfo.data);
      memset(buffer, 0, sizeof(double)*3*max_node_adj*sharedNodes[iSub]->num(jSub));
      for (int iNode = 0; iNode < sharedNodes[iSub]->num(jSub); ++iNode) {
        
        int cnt = 0;
        std::list<Vec3D>& nl = data[iSub][ (*sharedNodes[iSub])[jSub][iNode] ];
        for (std::list<Vec3D>::iterator it = nl.begin(); it != nl.end(); ++it) {
          for (int k = 0; k < 3; ++k,++cnt)
            buffer[iNode*max_node_adj*3+cnt] = (*it)[k];
        }
      }
    }
  }

  commPat.exchange();

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      SubRecInfo<double> sInfo = commPat.recData(subD.getRcvChannel()[jSub]);
      double * buffer = reinterpret_cast<double *>(sInfo.data);

      for (int iNode = 0; iNode < sharedNodes[iSub]->num(jSub); ++iNode) {

        std::list<Vec3D>& nl = data[iSub][ (*sharedNodes[iSub])[jSub][iNode] ];
        for (int cnt = 0; cnt < max_node_adj; ++cnt) {

          double* b = buffer+iNode*max_node_adj*3+cnt*3;
          Vec3D v(b[0],b[1],b[2]);
          // A norm of 0.0 signifies the end of the list for this node
          if (v.norm() == 0.0) 
            break;
          nl.push_back(v);
        }
      }
    }
  }
}


};
   
template<class Scalar>
template<class T,int dim>
void MultiGridLevel<Scalar>::assembleInternal(DistMat<T,dim> & A) const {

  
  
#pragma omp parallel for
  for (int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {

    MultigridSubdomain* D = &mgSubdomains[iSub];
    int glSub = domain.getSubDomain()[iSub]->getGlobSubNum();

    for (std::map<int,NeighborDomain*>::iterator it = D->neighbors.begin();
	 it != D->neighbors.end(); ++it) {
      
      NeighborDomain* N = it->second;
      int snd = sndChannel(glSub, it->first);
      SubRecInfo<T> sInfo = (dim == neq1? matPattern : matPattern2)->getSendBuffer(snd);
      T (*buffer)[dim*dim] = reinterpret_cast<T (*)[dim*dim]>(sInfo.data);
    
      for (int iNode = 0; iNode < N->sharedNodes.size(); ++iNode) {
        T * a = A(iSub).getElem_ii(N->sharedNodes[iNode].second);
        for (int j=0; j<dim*dim; ++j) buffer[iNode][j] = a[j];
      }
    }
    
  }

  (dim == neq1? matPattern : matPattern2)->exchange();

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {

    MultigridSubdomain* D = &mgSubdomains[iSub];
    int glSub = domain.getSubDomain()[iSub]->getGlobSubNum();

    for (std::map<int,NeighborDomain*>::iterator it = D->neighbors.begin();
	 it != D->neighbors.end(); ++it) {
      
      NeighborDomain* N = it->second;

      int rcv = rcvChannel(glSub, it->first);
      SubRecInfo<T> sInfo = (dim == neq1? matPattern : matPattern2)->recData(rcv);
      T (*buffer)[dim*dim] = reinterpret_cast<T (*)[dim*dim]>(sInfo.data);

      for (int iNode = 0; iNode < N->sharedNodes.size(); ++iNode) {
        T * a = A(iSub).getElem_ii(N->sharedNodes[iNode].second);
        for (int j = 0; j < dim*dim; ++j) a[j] += buffer[iNode][j];
      }
    }
  }

  
}
   
template<class Scalar>
template<class T,int dim>
void MultiGridLevel<Scalar>::assembleInternal(DistSVec<T,dim> & A) const {

  
#pragma omp parallel for
  for (int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {

    MultigridSubdomain* D = &mgSubdomains[iSub];
    int glSub = domain.getSubDomain()[iSub]->getGlobSubNum();

    for (std::map<int,NeighborDomain*>::iterator it = D->neighbors.begin();
	 it != D->neighbors.end(); ++it) {
      
      NeighborDomain* N = it->second;
      int snd = sndChannel(glSub, it->first);
      SubRecInfo<T> sInfo = nodeVecPattern->getSendBuffer(snd);
      T (*buffer)[dim] = reinterpret_cast<T (*)[dim]>(sInfo.data);
    
      for (int iNode = 0; iNode < N->sharedNodes.size(); ++iNode) {

        T * a = A(iSub)[N->sharedNodes[iNode].second];
        for (int j=0; j<dim; ++j) buffer[iNode][j] = a[j];
      }
    }
    
  }

  nodeVecPattern->exchange();

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {

    MultigridSubdomain* D = &mgSubdomains[iSub];
    int glSub = domain.getSubDomain()[iSub]->getGlobSubNum();

    for (std::map<int,NeighborDomain*>::iterator it = D->neighbors.begin();
	 it != D->neighbors.end(); ++it) {
      
      NeighborDomain* N = it->second;

      int rcv = rcvChannel(glSub, it->first);
      SubRecInfo<T> sInfo = nodeVecPattern->recData(rcv);
      T (*buffer)[dim] = reinterpret_cast<T (*)[dim]>(sInfo.data);

      for (int iNode = 0; iNode < N->sharedNodes.size(); ++iNode) {
        T * a = A(iSub)[N->sharedNodes[iNode].second];
        for (int j = 0; j < dim; ++j) a[j] += buffer[iNode][j];
      }
    }
  }

  
}
  
template<class Scalar>
template<class T>
void MultiGridLevel<Scalar>::assembleInternal(DistVec<T> & A) const {

  
  CommPattern<T>* compat = internal_select<T>(nodeIdPattern, nodeVolPattern);

#pragma omp parallel for
  for (int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {

    MultigridSubdomain* D = &mgSubdomains[iSub];
    int glSub = domain.getSubDomain()[iSub]->getGlobSubNum();

    for (std::map<int,NeighborDomain*>::iterator it = D->neighbors.begin();
	 it != D->neighbors.end(); ++it) {
      
      NeighborDomain* N = it->second;
      int snd = sndChannel(glSub, it->first);
      SubRecInfo<T> sInfo = compat->getSendBuffer(snd);
      T *buffer = reinterpret_cast<T *>(sInfo.data);
    
      for (int iNode = 0; iNode < N->sharedNodes.size(); ++iNode) {
        T * a = &A(iSub)[N->sharedNodes[iNode].second];
        buffer[iNode] = *a;
      }
    }
    
  }

  compat->exchange();

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {

    MultigridSubdomain* D = &mgSubdomains[iSub];
    int glSub = domain.getSubDomain()[iSub]->getGlobSubNum();

    for (std::map<int,NeighborDomain*>::iterator it = D->neighbors.begin();
	 it != D->neighbors.end(); ++it) {
      
      NeighborDomain* N = it->second;

      int rcv = rcvChannel(glSub, it->first);
      SubRecInfo<T> sInfo = compat->recData(rcv);
      T *buffer = reinterpret_cast<T *>(sInfo.data);

      for (int iNode = 0; iNode < N->sharedNodes.size(); ++iNode) {
        T * a = &A(iSub)[N->sharedNodes[iNode].second];
        *a += buffer[iNode];
      }
    }
  }

  
}

template<class Scalar>
template<class T>
void MultiGridLevel<Scalar>::
assembleInternalMax(DistVec<T> & A) const {

  
  CommPattern<T>* compat = internal_select<T>(nodeIdPattern, nodeVolPattern);
#pragma omp parallel for
  for (int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {

    MultigridSubdomain* D = &mgSubdomains[iSub];
    int glSub = domain.getSubDomain()[iSub]->getGlobSubNum();

    for (std::map<int,NeighborDomain*>::iterator it = D->neighbors.begin();
	 it != D->neighbors.end(); ++it) {
      
      NeighborDomain* N = it->second;
      int snd = sndChannel(glSub, it->first);
      SubRecInfo<T> sInfo = compat->getSendBuffer(snd);

      T *buffer = reinterpret_cast<T *>(sInfo.data);
    
      for (int iNode = 0; iNode < N->sharedNodes.size(); ++iNode) {
        T * a = &A(iSub)[N->sharedNodes[iNode].second];
        buffer[iNode] = *a;// std::max<T>(*a,buffer[iNode]);
      }
    }
    
  }

  compat->exchange();

#pragma omp parallel for
  for(int iSub = 0; iSub < domain.getNumLocSub(); ++iSub) {

    MultigridSubdomain* D = &mgSubdomains[iSub];
    int glSub = domain.getSubDomain()[iSub]->getGlobSubNum();

    for (std::map<int,NeighborDomain*>::iterator it = D->neighbors.begin();
	 it != D->neighbors.end(); ++it) {
      
      NeighborDomain* N = it->second;

      int rcv = rcvChannel(glSub, it->first);
      SubRecInfo<T> sInfo = compat->recData(rcv);
      T *buffer = reinterpret_cast<T *>(sInfo.data);

      for (int iNode = 0; iNode < N->sharedNodes.size(); ++iNode) {
        T * a = &A(iSub)[N->sharedNodes[iNode].second];
        *a = std::max<T>(buffer[iNode],*a);
      }
    }
  }

 
}

template <class Scalar>
void MultiGridLevel<Scalar>::mapNodeList(int iSub,Aerof_unordered_set<int>::type&  l) {

  Aerof_unordered_set<int>::type tmp(l);
  l.clear();

  for (Aerof_unordered_set<int>::type::iterator it = tmp.begin();
       it != tmp.end(); ++it) {

    l.insert(parent->mapFineToCoarse(iSub,*it));
  }
}

template <class Scalar>
void MultiGridLevel<Scalar>::identifyEdges(CommPattern<int> &edgeNumPat,
                                           int mySub) {
 
  int iSub, jSub;
  SubDomain& subD(*domain.getSubDomain()[mySub]);
  
  Connectivity *nodeToNeighb = sharedNodes[mySub]->reverse();

  int numNeighb = subD.getNumNeighb();

  numSharedEdges[mySub] = new int[numNeighb];

  for (iSub = 0; iSub < numNeighb; ++iSub) numSharedEdges[mySub][iSub] = 0;

  int (*edgePtr)[2] = edges[mySub]->getPtr();

  int l;
  for (l=0; l<edges[mySub]->size(); ++l) {

    int leftNode = edgePtr[l][0];
    int rightNode = edgePtr[l][1];

    /*if (mySub == 1 && leftNode == 384 && rightNode == 429)
      std::cout << "Hello" << std::endl;
*/
    if (nodeToNeighb->num(leftNode) == 0 ||
	nodeToNeighb->num(rightNode) == 0) continue;

    for (iSub = 0; iSub < nodeToNeighb->num(leftNode); ++iSub) {

      int subI = (*nodeToNeighb)[leftNode][iSub];

      for (jSub = 0; jSub < nodeToNeighb->num(rightNode); ++jSub)
	if (subI == (*nodeToNeighb)[rightNode][jSub]) numSharedEdges[mySub][subI] += 1;

    }

  }

  sharedEdges[mySub] = new EdgeDef *[numNeighb];

  for (iSub = 0; iSub < numNeighb; ++iSub)
    sharedEdges[mySub][iSub] = new EdgeDef[ numSharedEdges[mySub][iSub] ];

  for (iSub = 0; iSub < numNeighb; ++iSub) numSharedEdges[mySub][iSub] = 0;

  for (l=0; l<edges[mySub]->size(); ++l) {

    int leftNode = edgePtr[l][0];
    int rightNode = edgePtr[l][1];

    if (nodeToNeighb->num(leftNode) == 0 ||
       nodeToNeighb->num(rightNode) == 0) continue;

    for (iSub = 0; iSub < nodeToNeighb->num(leftNode); ++iSub) {

      int subI = (*nodeToNeighb)[leftNode][iSub];

      for (jSub = 0; jSub < nodeToNeighb->num(rightNode); ++jSub) {

	if (subI == (*nodeToNeighb)[rightNode][jSub]) {

	  sharedEdges[mySub][subI][ numSharedEdges[mySub][subI] ].edgeNum = l;
	  sharedEdges[mySub][subI][ numSharedEdges[mySub][subI] ].glLeft = locToGlobMap[mySub][leftNode];
	  sharedEdges[mySub][subI][ numSharedEdges[mySub][subI] ].glRight = locToGlobMap[mySub][rightNode];
	  sharedEdges[mySub][subI][ numSharedEdges[mySub][subI] ].order();

	  numSharedEdges[mySub][subI] += 1;

	}

      }

    }
  }

  for (iSub = 0; iSub < numNeighb; ++iSub) {
#ifdef OLD_STL
    sort(sharedEdges[mySub][iSub], sharedEdges[mySub][iSub]+numSharedEdges[mySub][iSub]);
#else
    stable_sort(sharedEdges[mySub][iSub], sharedEdges[mySub][iSub]+numSharedEdges[mySub][iSub]);
#endif
    edgeNumPat.setLen(subD.getSndChannel()[iSub], 2*numSharedEdges[mySub][iSub]);
  }

  delete nodeToNeighb;

}

template <class Scalar>
void MultiGridLevel<Scalar>::sndEdgeInfo(CommPattern<int> &edgeNumPat,int mySub)
{

  int iSub, iEdge;
  SubDomain& subD(*domain.getSubDomain()[mySub]);
  
  int numNeighb = subD.getNumNeighb();

  for (iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<int> sInfo = edgeNumPat.getSendBuffer(subD.getSndChannel()[iSub]);

    for (iEdge = 0; iEdge < numSharedEdges[mySub][iSub]; ++iEdge) {
      sInfo.data[2*iEdge]   = sharedEdges[mySub][iSub][iEdge].glLeft;
      sInfo.data[2*iEdge+1] = sharedEdges[mySub][iSub][iEdge].glRight;
    }

  }
}

template <class Scalar>
void MultiGridLevel<Scalar>::rcvEdgeInfo(CommPattern<int> &edgeNumPat,int mySub,int dim)
{

  int iSub, iEdge;

  SubDomain& subD(*domain.getSubDomain()[mySub]);
  
  // We receive the neighbor's list of edges and then also build the edgeMasterFlag

  bool *edgeMasterFlag = new bool[edges[mySub]->size()];

  for (iEdge = 0; iEdge < edges[mySub]->size(); ++iEdge) edgeMasterFlag[iEdge] = true;

  for (iSub = 0; iSub < subD.getNumNeighb(); ++iSub) {

    SubRecInfo<int> sInfo = edgeNumPat.recData(subD.getRcvChannel()[iSub]);

    int nIndex = 0;
    int myIndex = 0;

    for (iEdge = 0; iEdge < numSharedEdges[mySub][iSub]; ++iEdge) {

      int glLeft  = sharedEdges[mySub][iSub][iEdge].glLeft;
      int glRight = sharedEdges[mySub][iSub][iEdge].glRight;

      while (2*nIndex < sInfo.len &&
	     (sInfo.data[2*nIndex] < glLeft ||
	      (sInfo.data[2*nIndex] == glLeft && sInfo.data[2*nIndex+1] < glRight)))
	nIndex++;

      if (2*nIndex < sInfo.len &&
	  (sInfo.data[2*nIndex] == glLeft && sInfo.data[2*nIndex+1] == glRight)) {
	sharedEdges[mySub][iSub][myIndex] = sharedEdges[mySub][iSub][iEdge];
	myIndex++;
      }

    }

    numSharedEdges[mySub][iSub] = myIndex;
    offDiagMatPattern->setLen(subD.getSndChannel()[iSub],2*dim*dim*numSharedEdges[mySub][iSub]);
    offDiagMatPattern->setLen(subD.getRcvChannel()[iSub],2*dim*dim*numSharedEdges[mySub][iSub]);
    
    edgeAreaPattern->setLen(subD.getSndChannel()[iSub],numSharedEdges[mySub][iSub]);
    edgeAreaPattern->setLen(subD.getRcvChannel()[iSub],numSharedEdges[mySub][iSub]);
    edgeVecPattern->setLen(subD.getSndChannel()[iSub],numSharedEdges[mySub][iSub]*3);
    edgeVecPattern->setLen(subD.getRcvChannel()[iSub],numSharedEdges[mySub][iSub]*3);

    if (subD.getNeighb()[iSub] < subD.getGlobSubNum()) // I cannot be the master of these edges
      for (iEdge = 0; iEdge < numSharedEdges[mySub][iSub]; ++iEdge)
	edgeMasterFlag[ sharedEdges[mySub][iSub][iEdge].edgeNum ] = false;

  }

  int num_edges = 0;
  for (iEdge = 0; iEdge < edges[mySub]->size(); ++iEdge) 
    if (edgeMasterFlag[iEdge]) num_edges++;
 // std::cout << num_edges << std::endl;
  
  edges[mySub]->setMasterFlag(edgeMasterFlag);
}

template<class Scalar>
Connectivity* MultiGridLevel<Scalar>::createEdgeBasedConnectivity(int iSub)
{

  int numNodesLocal = numNodes[iSub];

  int *numNeigh = reinterpret_cast<int *>(alloca(sizeof(int) * numNodesLocal));

  int i;
  for (i=0; i<numNodesLocal; ++i) numNeigh[i] = 1;

  int (*edgePtr)[2] = edges[iSub]->getPtr();

  int l;
  for (l=0; l<edges[iSub]->size(); ++l) {
    ++numNeigh[ edgePtr[l][0] ];
    ++numNeigh[ edgePtr[l][1] ];
  }

  int nnz = 0;
  for (i=0; i<numNodesLocal; ++i) nnz += numNeigh[i];

  if (nnz != (2*edges[iSub]->size() + numNodesLocal)) {
    fprintf(stderr,"*** Error: wrong number of nonzero blocks\n");
    exit(1);
  }

  // construction of ia

  int *ia = new int[numNodesLocal+1];

  ia[0] = 0;

  for (i=0; i<numNodesLocal; ++i) ia[i+1] = ia[i] + numNeigh[i];

  // construction of ja

  int *ja = new int[nnz];

  for (i=0; i<numNodesLocal; ++i) {
    ja[ia[i]] = i;
    numNeigh[i] = 1;
  }

  for (l=0; l<edges[iSub]->size(); ++l) {

    int i1 = edgePtr[l][0];
    int i2 = edgePtr[l][1];

    ja[ ia[i1] + numNeigh[i1] ] = i2;
    ++numNeigh[i1];

    ja[ ia[i2] + numNeigh[i2] ] = i1;
    ++numNeigh[i2];

  }

  for (i=0; i<numNodesLocal; ++i)
#ifdef OLD_STL
    sort(ja+ia[i], ja+ia[i+1]);
#else
    stable_sort(ja+ia[i], ja+ia[i+1]);
#endif

  Connectivity *nodeToNode = new Connectivity(numNodesLocal, ia, ja);

  return nodeToNode;

}

template<class Scalar>
compStruct* MultiGridLevel<Scalar>::createRenumbering(int iSub,Connectivity *nodeToNode,
					 int typeRenum, int print)
{

  int numNodesLocal = numNodes[iSub];

  compStruct *nodeRenum = nodeToNode->renumByComponent(typeRenum);

  int *order = new int[numNodesLocal];

  int i;
  for (i=0; i<numNodesLocal; ++i) order[i] = -1;

  for (i=0; i<numNodesLocal; ++i)
    if (nodeRenum->renum[i] >= 0) order[nodeRenum->renum[i]] = i;

  nodeRenum->order = order;

  int *ia = (*nodeToNode).ptr();
  int *ja = (*nodeToNode)[0];
  double aa;
  double (*a)[1] = reinterpret_cast<double (*)[1]>(&aa);

  SparseMat<double,1> A(numNodesLocal, ia[numNodesLocal], ia, ja, a, 0, 0);

  //if (print == 1)
  //  F77NAME(psplotmask)(numNodesLocal, ia, ja, 0, 10+globSubNum);

  A.permute(nodeRenum->renum);

  //if (print == 1)
  //  F77NAME(psplotmask)(numNodesLocal, ia, ja, 0, 50+globSubNum);

  return nodeRenum;

}
template<class Scalar>
template<int dim>
SparseMat<Scalar,dim>* MultiGridLevel<Scalar>::createMaskILU(int iSub,int fill, 
                                                             int renum, int *ndType)
{

  nodeToNodeMaskILU[iSub] = createEdgeBasedConnectivity(iSub);

  compStruct *nodeRenum = createRenumbering(iSub,nodeToNodeMaskILU[iSub], renum, 0);

  int *ia = (*nodeToNodeMaskILU[iSub]).ptr();
  int *ja = (*nodeToNodeMaskILU[iSub])[0];
  int n = numNodes[iSub];
  int nnz = ia[n];

  SparseMat<Scalar,dim> *A = new SparseMat<Scalar,dim>(n, nnz, ia, ja, 0, nodeRenum, ndType);

  A->symbolicILU(fill);

  A->createPointers(*edges[iSub]);

  return A;

}
    

void PrintPriorityNodes(const PriorityNodes& P) {

  for (PriorityNodes::const_iterator itr = P.begin(); itr != P.end();
       ++itr) {

    std::cout << *itr << " ";
  }
}

template<class Scalar>
int MultiGridLevel<Scalar>::
getNodeWithMostAgglomeratedNeighbors(std::vector<std::set<int> >& C,//Connectivity* C, 
                                                PriorityNodes& P,
                                                Vec<int>& nodeMapping,int iSub) {

  int max_count = -1, curr_node = -1;
 
  Topology currentTopo = TopoUnknown; 
  for (PriorityNodes::iterator it = P.begin(); 
       it != P.end();) {

    int tmp = *it;
    int tmp_count = 0;
    if (nodeMapping[tmp] >= 0) {
      PriorityNodes::iterator it2 = it;
      P.erase(it2);
      if (it == P.end()) break;
      ++it;
      continue;
    }

    for (std::set<int>::iterator it2 = C[tmp].begin();
         it2 != C[tmp].end(); ++it2) {//int n = 0; n < C->num(tmp); ++n) {
      int n = *it2;
      if (nodeMapping[n] >= 0)
        ++tmp_count;
    }
    if (tmp_count > max_count/* && currentTopo == (*nodeTopology)(iSub)[tmp]*/) {

      max_count = tmp_count;
      curr_node = tmp;
    }/* else if (currentTopo > (*nodeTopology)(iSub)[tmp]) {

      max_count = tmp_count;
      currentTopo = (*nodeTopology)(iSub)[tmp];
      curr_node = tmp;
    }*/
    ++it;
  }
  if (curr_node >= 0)
    P.erase(P.find(curr_node));  
  return curr_node;
}
                                                

template<class Scalar>
void MultiGridLevel<Scalar>::agglomerate(const DistInfo& refinedNodeDistInfo,
                                         const DistInfo& refinedEdgeDistInfo,
                                         CommPattern<int>& refinedNodeIdPattern,
                                         Connectivity** refinedSharedNodes,
                                         Connectivity ** nToN, EdgeSet ** refinedEdges,
                                         EdgeDef*** refinedSharedEdges,
                                         int** refinedNumSharedEdges,
                                         Domain& domain,int dim,int lneq1,int lneq2,
                                         DistVec<Vec3D>& refinedEdgeNormals,
                                         DistVec<double>* refinedEdgeAreas,
                                         DistVec<double>& refinedVol,
                                         DistVec<int>* finestNodeMapping_p1,
                                         double beta,int maxNumNodesPerAgglom)
{
  static int cnt = 0;
  ++cnt;
  my_dim = dim;
  neq1 = lneq1;
  neq2 = lneq2;
  numNodes = new int[numLocSub];
  nodeToNodeMaskILU = new Connectivity*[numLocSub];

  nodeIdPattern = new CommPattern<int>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<int>::CopyOnSend);
  nodeVolPattern = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);
  nodeVecPattern = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);
  
  if (neq2 > 0) {
    nodeVecPatternEq1 = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);
    nodeVecPatternEq2 = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);
  }
  nodePosnPattern = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);
  matPattern = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);
  offDiagMatPattern = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);
  edgeAreaPattern = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);
  edgeVecPattern = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);

  Connectivity* subToSub = domain.getSubToSub();

  DistVec<int>& nodeMapping = *pNodeMapping;
  DistVec<int>& edgeMapping = *pEdgeMapping;

  nodeMapping = -1;
  edgeMapping = -1;

  faceMapping = new DistVec<int>(*parent->faceDistInfo);

  *faceMapping = -1;

  *lineMap = -1;

  int myLevel = getMyLevel();  
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  
  sharedNodes = new Connectivity*[numLocSub];

  computeNodeNormalClasses();

  DistVec<int> nodeBlacklist(*parent->nodeDistInfo);
  nodeBlacklist = 0;

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    for (int i = 0; i < parent->nodeDistInfo->subSize(iSub); ++i) {

      if (!parent->nodeDistInfo->getMasterFlag(iSub)[i]) continue;
      bool ok = false;
      int myid;
      for (int j = 0; j < nToN[iSub]->num(i); ++j) {

        if ((*nToN[iSub])[i][j] == i) continue;
        myid = (*nToN[iSub])[i][j];
        if (parent->nodeDistInfo->getMasterFlag(iSub)[(*nToN[iSub])[i][j]]) {
          ok = true;
          break;
        }
      }
      if (!ok)
        nodeBlacklist(iSub)[myid] = 1;
    }
  }

  operMax<int> maxOp;
  ::assemble(domain, refinedNodeIdPattern, refinedSharedNodes, nodeBlacklist, maxOp);

  DistVec<double> topo(*parent->nodeDistInfo);
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    for (int i = 0; i < topo.subSize(iSub); ++i)
      topo(iSub)[i] = static_cast<double>((*nodeTopology)(iSub)[i]);
  }
  
//  parent->writeXpostFile("nodeClass.xpost",topo);

  std::vector<Topology>* newTopology = new std::vector<Topology>[numLocSub];
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    std::vector<std::set<int> > myNtoN;
    myNtoN.resize(parent->nodeDistInfo->subSize(iSub));

    SubDomain& subD(*domain.getSubDomain()[iSub]);
    // We have to ensure that any agglomerated nodes are shared by the
    // correct subdomains (that they can talk to each other)
    Connectivity& rSharedNodes(*refinedSharedNodes[iSub]);
    std::map<int,std::set<int> > nodeSharedData;
    for(int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      for(int i = 0; i < rSharedNodes.num(jSub); ++i) {
        nodeSharedData[rSharedNodes[jSub][i]].insert(subD.getNeighb()[jSub]);
      }
    }

    int seed_node = 0, agglom_id = 0;
    for(int i = 0; i < nodeMapping(iSub).size(); ++i) { // Disable shared slave nodes
      if(!refinedNodeDistInfo.getMasterFlag(iSub)[i]) nodeMapping(iSub)[i]=0;
      else if (nodeBlacklist(iSub)[i]) {
        nodeMapping(iSub)[i] = agglom_id++;
      }

    }

    // Perform local agglomeration
    lineLengths[iSub] = new int[nodeMapping(iSub).size()];
  
    PriorityNodes nodesOnSurfaceAndFront, nodesOnSurface,
                  nodesOnSolidAndFront, nodesOnSolid,
                  nodesOnFarFieldAndFront, nodesOnFarField,
                  nodesOnSubDomainBoundaryAndFront, nodesOnSubDomainBoundary,
                  nodesOnFront;
    subD.getSurfaceNodes(nodesOnSurface);
    subD.getSolidBoundaryNodes(nodesOnSolid);
    subD.getFarFieldBoundaryNodes(nodesOnFarField);
    subD.getSubDomainBoundaryNodes(nodesOnSubDomainBoundary);

    mapNodeList(iSub,nodesOnSolid);
    mapNodeList(iSub,nodesOnSurface);
    mapNodeList(iSub,nodesOnFarField);
    mapNodeList(iSub,nodesOnSubDomainBoundary);

    // Its possible the previous functions will give the same node twice in two
    // different lists.  But this is ok.  Before we construct a node mapping we
    // will see if the node has already been mapped

    for (int l = 0 ; l < refinedEdges[iSub]->size(); ++l) {

      myNtoN[refinedEdges[iSub]->getPtr()[l][0]].insert(refinedEdges[iSub]->getPtr()[l][1]);
      myNtoN[refinedEdges[iSub]->getPtr()[l][1]].insert(refinedEdges[iSub]->getPtr()[l][0]);
    }
 
    int tmp; 
    while(1) {
 
      seed_node = -1;
      seed_node = getNodeWithMostAgglomeratedNeighbors(myNtoN,//nToN[iSub],
                                           nodesOnSurfaceAndFront,
                                           nodeMapping(iSub),iSub);
      
      if (seed_node < 0) {
      
        seed_node = getNodeWithMostAgglomeratedNeighbors(myNtoN,//nToN[iSub],
                                           nodesOnSurface,
                                           nodeMapping(iSub),iSub);
      }
 
      
      if (seed_node < 0) {
        seed_node = getNodeWithMostAgglomeratedNeighbors(myNtoN,//nToN[iSub],
                                             nodesOnSolidAndFront,
                                             nodeMapping(iSub),iSub);
      }

      if (seed_node < 0) {
      
        seed_node = getNodeWithMostAgglomeratedNeighbors(myNtoN,//nToN[iSub],
                                           nodesOnSolid,
                                           nodeMapping(iSub),iSub);
      }
 
      if (seed_node < 0) {

        seed_node = getNodeWithMostAgglomeratedNeighbors(myNtoN,//nToN[iSub],
                                           nodesOnFarFieldAndFront,
                                           nodeMapping(iSub),iSub);
      }
      
      if (seed_node < 0) {

        seed_node = getNodeWithMostAgglomeratedNeighbors(myNtoN,//nToN[iSub],
                                           nodesOnFarField,
                                           nodeMapping(iSub),iSub);
      }
      
      if (seed_node < 0) {
        
        seed_node = getNodeWithMostAgglomeratedNeighbors(myNtoN,//nToN[iSub],
                                           nodesOnSubDomainBoundaryAndFront,
                                           nodeMapping(iSub),iSub);

      }
      
      if (seed_node < 0) {

        seed_node = getNodeWithMostAgglomeratedNeighbors(myNtoN,//nToN[iSub],
                                           nodesOnSubDomainBoundary,
                                           nodeMapping(iSub),iSub);
      }
      
      if (seed_node < 0) {

        seed_node = getNodeWithMostAgglomeratedNeighbors(myNtoN,//nToN[iSub],
                                           nodesOnFront,
                                           nodeMapping(iSub),iSub);
      }
  
      if (seed_node < 0)
        break; 

      Topology currentTopology = (*nodeTopology)(iSub)[seed_node];

      std::vector<int> agglom;
      std::map<int,double> weights;
      agglom.reserve(8);

      agglom.push_back(seed_node);

      double currentAgglomVolume = refinedVol(iSub)[seed_node];

      if (currentTopology != TopoVertex) {

        std::set<int> agglomSubDs;
        agglomSubDs.insert(subD.getGlobSubNum());
        for (std::set<int>::iterator it = nodeSharedData[seed_node].begin();
             it != nodeSharedData[seed_node].end(); ++it)
          agglomSubDs.insert(*it);

        const int current_node = seed_node;
        //while(agglom.size() < 8) {
        //  const int current_node=agglom.back();
        double max_weight = 0.0;
        std::set<int> valid_neighbors;
        Vec3D currentNormal;
        std::map<int, Vec3D> neighboring_normal_sums;
        std::map<int, Vec3D> neighboring_normal_sums_agg;
        // std::map<int,Vec3D> myAggNormal1,myAggNormal2;
        if (currentTopology == TopoFace)
          currentNormal = (*topologyNormal)(iSub)[current_node];
        for(std::set<int>::iterator s_it = myNtoN[current_node].begin(); 
            s_it != myNtoN[current_node].end(); ++s_it) {
          //int n = 0; n < nToN[iSub]->num(current_node); ++n) {
          const int neighbor = *s_it;//(*nToN[iSub])[current_node][n];
          if (neighbor == current_node)
            continue;
          if (neighbor < current_node) {
            neighboring_normal_sums[neighbor] = 
              refinedEdgeNormals(iSub)[refinedEdges[iSub]->findOnly(current_node,neighbor)];
            //myAggNormal1[neighbor] = -refinedEdgeNormals(iSub)[refinedEdges[iSub]->findOnly(current_node,neighbor)];
          } else {
            neighboring_normal_sums[neighbor] = -
              refinedEdgeNormals(iSub)[refinedEdges[iSub]->findOnly(current_node,neighbor)];
            //myAggNormal1[neighbor] = refinedEdgeNormals(iSub)[refinedEdges[iSub]->findOnly(current_node,neighbor)];
          }
          if (nodeMapping(iSub)[neighbor] >= 0) {
            if (neighbor < current_node) {
              neighboring_normal_sums_agg[nodeMapping(iSub)[neighbor]] += 
                refinedEdgeNormals(iSub)[refinedEdges[iSub]->findOnly(current_node,neighbor)];
           //   myAggNormal2[nodeMapping(iSub)[neighbor]] -= refinedEdgeNormals(iSub)[refinedEdges[iSub]->findOnly(current_node,neighbor)];
            }
            else {
              neighboring_normal_sums_agg[nodeMapping(iSub)[neighbor]] -=
                refinedEdgeNormals(iSub)[refinedEdges[iSub]->findOnly(current_node,neighbor)];
             //   myAggNormal2[nodeMapping(iSub)[neighbor]] += refinedEdgeNormals(iSub)[refinedEdges[iSub]->findOnly(current_node,neighbor)];
            }
          }

          //max_weight = std::max(max_weight, refinedEdgeNormals(iSub)[refinedEdges[iSub]->findOnly(current_node, neighbor)].norm());
          max_weight = std::max(max_weight, 1.0/refinedEdges[iSub]->viewEdgeLength()[refinedEdges[iSub]->findOnly(current_node, neighbor)]);
          const Topology& Toth = (*nodeTopology)(iSub)[neighbor];
          /*if (Toth == TopoVertex) continue;
          if (currentTopology == TopoLine) {

            if (Toth != TopoLine) continue;
            if (fabs((*topologyNormal)(iSub)[current_node]*(*topologyNormal)(iSub)[neighbor]) < 
                mesh_topology_threshold)
              continue;
            Vec3D a(parent->getXn()(iSub)[current_node][0]-parent->getXn()(iSub)[neighbor][0],
                    parent->getXn()(iSub)[current_node][1]-parent->getXn()(iSub)[neighbor][1],
                    parent->getXn()(iSub)[current_node][2]-parent->getXn()(iSub)[neighbor][2]);
            a /= a.norm();
            if (fabs(a*(*topologyNormal)(iSub)[current_node]) < mesh_topology_threshold)
              continue;
          } else if (Toth == TopoLine)
            continue; 
          else if (currentTopology == TopoFace) {

            if (Toth == TopoFace) {
 
              if (currentNormal*(*topologyNormal)(iSub)[neighbor] < mesh_topology_threshold)
                continue;
            } else
              continue;
          }*/
          if (currentTopology == TopoInterior) {

            if (Toth != TopoInterior) {
              continue;
            }
          }

          if(current_node == neighbor || nodeMapping(iSub)[neighbor] >= 0) continue;
          bool ok = true;
          bool has_new = false;
          bool is_superset = true;
          for (std::set<int>::iterator it = nodeSharedData[neighbor].begin();
               it != nodeSharedData[neighbor].end(); ++it) {
          
            if (agglomSubDs.find(*it) == agglomSubDs.end()) {

              has_new = true;
              break;
            }
          }

          for (std::set<int>::iterator it = agglomSubDs.begin();
               it != agglomSubDs.end(); ++it) {

            if (nodeSharedData[neighbor].find(*it) == nodeSharedData[neighbor].end()) {

              is_superset = false;
              break;
            }
          }
          if (has_new && !is_superset) continue;
        }
        
        for(std::set<int>::iterator s_it = myNtoN[current_node].begin(); 
            s_it != myNtoN[current_node].end(); ++s_it) {
          const int neighbor = *s_it;//(*nToN[iSub])[current_node][n];
 
        /*for(int n = 0; n < nToN[iSub]->num(current_node); ++n) {
          const int neighbor = (*nToN[iSub])[current_node][n];
        */
          if(current_node == neighbor || nodeMapping(iSub)[neighbor] >= 0) continue;

          const Topology& Toth = (*nodeTopology)(iSub)[neighbor];
          if (Toth == TopoVertex) continue;
          if (currentTopology == TopoLine) {

            if (Toth != TopoLine) continue;
            if (fabs((*topologyNormal)(iSub)[current_node]*(*topologyNormal)(iSub)[neighbor]) < 
                mesh_topology_threshold)
              continue;
            Vec3D a(parent->getXn()(iSub)[current_node][0]-parent->getXn()(iSub)[neighbor][0],
                    parent->getXn()(iSub)[current_node][1]-parent->getXn()(iSub)[neighbor][1],
                    parent->getXn()(iSub)[current_node][2]-parent->getXn()(iSub)[neighbor][2]);
            a /= a.norm();
            if (fabs(a*(*topologyNormal)(iSub)[current_node]) < mesh_topology_threshold)
              continue;
          } else if (Toth == TopoLine)
            continue; 
          else if (currentTopology == TopoFace) {

            if (Toth == TopoFace) {
 
              if (currentNormal*(*topologyNormal)(iSub)[neighbor] < mesh_topology_threshold)
                continue;
            } else
              continue;
          } else { // currentTopology == TopoInterior

            if (Toth == TopoFace) {
              continue;
              currentNormal = (*topologyNormal)(iSub)[neighbor];
              currentTopology = TopoFace;
            }
          }

          bool ok = true;
          bool has_new = false;
          bool is_superset = true;
          for (std::set<int>::iterator it = nodeSharedData[neighbor].begin();
               it != nodeSharedData[neighbor].end(); ++it) {
          
            if (agglomSubDs.find(*it) == agglomSubDs.end()) {

              has_new = true;
              break;
            }
          }

          for (std::set<int>::iterator it = agglomSubDs.begin();
               it != agglomSubDs.end(); ++it) {

            if (nodeSharedData[neighbor].find(*it) == nodeSharedData[neighbor].end()) {

              is_superset = false;
              break;
            }
          }
          if (has_new && !is_superset) continue;
          //if (beta*max_weight <  
          //    refinedEdgeNormals(iSub)[refinedEdges[iSub]->findOnly(current_node, neighbor)].norm()) {
        
          if (beta*max_weight < 1.0/refinedEdges[iSub]->viewEdgeLength()[refinedEdges[iSub]->findOnly(current_node, neighbor)]) {
            bool encloses_node = false;
            std::map<int,Vec3D> tmpSums;
            for(std::set<int>::iterator s_it2 = myNtoN[neighbor].begin(); 
               s_it2 != myNtoN[neighbor].end(); ++s_it2) {
              //for(int np = 0; np < nToN[iSub]->num(neighbor); ++np) {

              int n2 = *s_it2;//(*nToN[iSub])[neighbor][np];
              if (neighbor == n2)
                continue;
              int k;
              for (k = 0; k < agglom.size(); ++k) {
                if (agglom[k] == n2) break;
              }
              if (n2 == current_node || k < agglom.size())
                continue;
              Vec3D a = neighboring_normal_sums[n2],
                    b = refinedEdgeNormals(iSub)[refinedEdges[iSub]->findOnly(neighbor,n2)];
              if (n2 < neighbor) {
                if ((a+b).norm() < 1.0e-12)
                  encloses_node = true;
              } else  {
                if ((a-b).norm() < 1.0e-12)
                  encloses_node = true;
              } 

              if (nodeMapping(iSub)[n2] >= 0) {
                a = neighboring_normal_sums_agg[nodeMapping(iSub)[n2]],
                  b = refinedEdgeNormals(iSub)[refinedEdges[iSub]->findOnly(neighbor,n2)];
                if (tmpSums.find(nodeMapping(iSub)[n2]) == tmpSums.end())
                  tmpSums[nodeMapping(iSub)[n2]] = a;
                if (n2 < neighbor) {
                  tmpSums[nodeMapping(iSub)[n2]] += b;
                } else  {
                  tmpSums[nodeMapping(iSub)[n2]] -= b;
                }
              }
            }
            for (std::map<int,Vec3D>::iterator my_ts_it = tmpSums.begin();
                 my_ts_it != tmpSums.end(); ++my_ts_it) {
              if (my_ts_it->second.norm() < 1.0e-12)
                encloses_node = true;
            }
            if (encloses_node)
              continue;
            for(std::set<int>::iterator s_it2 = myNtoN[neighbor].begin(); 
               s_it2 != myNtoN[neighbor].end(); ++s_it2) {
 //           for(int np = 0; np < nToN[iSub]->num(neighbor); ++np) {

              int n2 = *s_it2;//(*nToN[iSub])[neighbor][np];
              if (neighbor == n2)
                continue;
              int k;
              for (k = 0; k < agglom.size(); ++k) {
                if (agglom[k] == n2) break;
              }
              if (n2 == current_node || k < agglom.size())
                continue;
              Vec3D& a = neighboring_normal_sums[n2],
                     b = refinedEdgeNormals(iSub)[refinedEdges[iSub]->findOnly(neighbor,n2)];
              if (n2 < neighbor) {
                a += b;
              } else {
                a -= b; 
              } 
 
              if (nodeMapping(iSub)[n2] >= 0) {
                Vec3D& ap = neighboring_normal_sums_agg[nodeMapping(iSub)[n2]],
                  bp = refinedEdgeNormals(iSub)[refinedEdges[iSub]->findOnly(neighbor,n2)];
                if (n2 < neighbor) {
                  ap += bp;
                } else  {
                  ap -= bp;
                }
              }
 
            }
           
            if ((currentAgglomVolume+refinedVol(iSub)[neighbor]) > 
                getTotalMeshVolume()/250.0) {

              continue;
            }
            if (agglom.size() >= maxNumNodesPerAgglom)
              continue;

            currentAgglomVolume += refinedVol(iSub)[neighbor];
            
            agglom.push_back(neighbor);
            
            for (std::set<int>::iterator it = nodeSharedData[neighbor].begin();
                 it != nodeSharedData[neighbor].end(); ++it) {
              agglomSubDs.insert(*it);
            }
           
          }
          
        }
      }
      int j = 0;
      bool isLine = true;

      for(std::vector<int>::const_iterator iter = agglom.begin(); iter != agglom.end(); iter++,++j) {
        const int current_node = *iter;
        //nodeMapping(iSub)[current_node] = agglom_id;
        for(std::set<int>::iterator s_it = myNtoN[current_node].begin(); 
            s_it != myNtoN[current_node].end() && isLine; ++s_it) {
//        for(int n = 0; n < nToN[iSub]->num(current_node) && isLine; ++n) {
          const int neighbor = *s_it;//(*nToN[iSub])[current_node][n];
          if (current_node == neighbor) continue;
  
          for (int k = 0; k < agglom.size(); ++k) {
            if ((k > j+1 || k < j-1) && neighbor == agglom[k]) { 
              isLine = false;
              break;
            }
          }
        }
   
        nodeMapping(iSub)[*iter] = agglom_id;
      }
        
      for(std::vector<int>::const_iterator iter = agglom.begin(); iter != agglom.end(); iter++) {
        const int current_node = *iter;
        for(std::set<int>::iterator s_it = myNtoN[current_node].begin(); 
            s_it != myNtoN[current_node].end(); ++s_it) {
        //for(int n = 0; n < nToN[iSub]->num(current_node); ++n) {
          const int neighbor = *s_it;//(*nToN[iSub])[current_node][n];
          if (nodeMapping(iSub)[neighbor] < 0) {

            if (nodesOnSurface.find(neighbor) != nodesOnSurface.end()) {

              nodesOnSurface.erase(nodesOnSurface.find(neighbor));
              nodesOnSurfaceAndFront.insert(neighbor);
            } else if (nodesOnSolid.find(neighbor) != nodesOnSolid.end()) {

              nodesOnSolid.erase(nodesOnSolid.find(neighbor));
              nodesOnSolidAndFront.insert(neighbor);
            } else if (nodesOnFarField.find(neighbor) != nodesOnFarField.end()) {

              nodesOnFarField.erase(nodesOnFarField.find(neighbor));
              nodesOnFarFieldAndFront.insert(neighbor);
            } else if (nodesOnSubDomainBoundary.find(neighbor) != nodesOnSubDomainBoundary.end()) {

              nodesOnSubDomainBoundary.erase(nodesOnSubDomainBoundary.find(neighbor));
              nodesOnSubDomainBoundaryAndFront.insert(neighbor);
            } else {

              nodesOnFront.insert(neighbor);
            }
          }          
        }
      }
 
      if (isLine && agglom.size() > 1) {
        int k;
        for (k = 0; k < agglom.size(); ++k) {
  
          assert(refinedNodeDistInfo.getMasterFlag(iSub)[agglom[k]]);
          (*lineLocIDMap)(iSub)[agglom[k]] = k;         
 
          (*lineIDMap)(iSub)[agglom[k]] = numLines[iSub];
          if (k > 0)
            (*lineMap)(iSub)[agglom[k]][0] = agglom[k-1];
          else
            (*lineMap)(iSub)[agglom[k]][0] = -1;
          if (k < agglom.size()-1)
            (*lineMap)(iSub)[agglom[k]][1] = agglom[k+1];  
          else
            (*lineMap)(iSub)[agglom[k]][1] = -1;
          lineids[iSub].push_back(agglom[k]);
        }
        lineLengths[iSub][numLines[iSub]] = k;

        for (; k < 8; ++k)
          lineids[iSub].push_back(-1);

        ++numLines[iSub];
      } else {
        for (int k = 0; k < agglom.size(); ++k) {
          (*lineMap)(iSub)[agglom[k]][0] = (*lineMap)(iSub)[agglom[k]][1] = -1;
        }
      }

      newTopology[iSub].push_back(currentTopology);
      ++agglom_id;
    }
 
    for (seed_node = 0; 
         seed_node < nodeMapping(iSub).size(); ++seed_node) {
      if (nodeMapping(iSub)[seed_node] < 0) {
        nodeMapping(iSub)[seed_node] = agglom_id++;
        newTopology[iSub].push_back((*nodeTopology)(iSub)[seed_node]);
      }
    }

    numNodes[iSub] = agglom_id;
  }

  // Check for nodes that are interior, but only have face neighbors.
/*
  DistVec<int> nodeStatus(refinedNodeDistInfo);
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    for (int i = 0; i < nodeStatus(iSub).size(); ++i) {
      if (newTopology[iSub][nodeMapping(iSub)[i]] != TopoInterior)
        nodeStatus(iSub)[i] = 1;
      else
        nodeStatus(iSub)[i] = 0;
    }

    int (*edgePtr)[2] = refinedEdges[iSub]->getPtr();
    for(int l = 0; l < refinedEdges[iSub]->size(); ++l) {
      int i = edgePtr[l][0], j = edgePtr[l][1];
      if (nodeMapping(iSub)[i] == nodeMapping(iSub)[j]) 
        continue;
      if (nodeMapping(iSub)[i] < 0 || nodeMapping(iSub)[j] < 0 )
        continue;
      if (newTopology[iSub][nodeMapping(iSub)[i]] == TopoInterior &&
          newTopology[iSub][nodeMapping(iSub)[j]] == TopoInterior) {
        nodeStatus(iSub)[i] = nodeStatus(iSub)[j] = 1;
      }
    }
  }

  ::assemble(domain, refinedNodeIdPattern, refinedSharedNodes, nodeStatus, maxOp);

  // If any of these nodes have nodeStatus = 0, then they are interior but have only
  // face neighbors  
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
  
    for (int i = 0; i < nodeStatus(iSub).size(); ++i) { 

      if (!nodeStatus(iSub)[i]) {
        int rnk;
        MPI_Comm_rank(MPI_COMM_WORLD,&rnk);
        
        std::cout << "Found weird node: " << rnk << " " << iSub << " " << i << " " << nodeMapping(iSub)[i] << std::endl;
      }
    }
  }
*/
// ------------------------- MPI PARALLEL CRAP -----------------------------------------------------------------
  DistVec<int> ownerSubDomain(refinedNodeDistInfo);
 
  maxNodesPerSubD = 0;
  for(int iSub = 0; iSub < numLocSub; ++iSub) maxNodesPerSubD = max(maxNodesPerSubD, nodeMapping(iSub).size());
 
  domain.getCommunicator()->globalMax(1, &maxNodesPerSubD);
 
  int tot_nodes = 0;
  for(int iSub = 0; iSub < numLocSub; ++iSub) tot_nodes += numNodes[iSub];
  domain.getCommunicator()->globalSum(1,&tot_nodes);

  domain.getCommunicator()->fprintf(stdout,"Total # of nodes on coarse grid: %d\n", tot_nodes);
 
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub)
    for(int i = 0; i < ownerSubDomain(iSub).size(); ++i)
      if(!refinedNodeDistInfo.getMasterFlag(iSub)[i])
        nodeMapping(iSub)[i] = ownerSubDomain(iSub)[i] = -1;
      else
        ownerSubDomain(iSub)[i] = domain.getSubDomain()[iSub]->getGlobSubNum();

  ::assemble(domain, refinedNodeIdPattern, refinedSharedNodes, ownerSubDomain, maxOp);
  ::assemble(domain, refinedNodeIdPattern, refinedSharedNodes, nodeMapping, maxOp);
 
  toTransfer = new std::map<int,std::set<int> >[numLocSub]; 
  std::set<int> globNodeIds;
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    Connectivity& rSharedNodes(*refinedSharedNodes[iSub]);

    std::map<int,std::map<int,int> > newlyInsertedNodes;
    std::map<int,int> newNodeIds;

    int * numNeighbors = new int[subD.getNumNeighb()];
    for(int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      for(int i = 0; i < rSharedNodes.num(jSub); ++i) {
        const int index = rSharedNodes[jSub][i];
        assert(nodeMapping(iSub)[index] >= 0);
        const int uniqueNodeID = ownerSubDomain(iSub)[index] * maxNodesPerSubD + nodeMapping(iSub)[index];
        if(refinedNodeDistInfo.getMasterFlag(iSub)[index]) {
          toTransfer[iSub][jSub].insert(nodeMapping(iSub)[index]);
          newlyInsertedNodes[jSub][uniqueNodeID] = nodeMapping(iSub)[index];
        } else {
          if(newNodeIds.find(uniqueNodeID) == newNodeIds.end())
            newNodeIds[uniqueNodeID] = numNodes[iSub]++;
          toTransfer[iSub][jSub].insert(newNodeIds[uniqueNodeID]);
          newlyInsertedNodes[jSub][uniqueNodeID] = newNodeIds[uniqueNodeID];
        }
      }
    }

    for(int i = 0; i < nodeMapping(iSub).size(); ++i) {
      const int uniqueNodeID = ownerSubDomain(iSub)[i] * maxNodesPerSubD + nodeMapping(iSub)[i];
      globNodeIds.insert(uniqueNodeID);
      if(ownerSubDomain(iSub)[i] != subD.getGlobSubNum()) {
        nodeMapping(iSub)[i] = newNodeIds[ownerSubDomain(iSub)[i] * maxNodesPerSubD + nodeMapping(iSub)[i]];
        locToGlobMap[iSub][nodeMapping(iSub)[i]] = uniqueNodeID;
      } else
        locToGlobMap[iSub][nodeMapping(iSub)[i]] = uniqueNodeID;
    }

    for(std::map<int,std::set<int> >::const_iterator iter=toTransfer[iSub].begin(); iter != toTransfer[iSub].end(); ++iter) {
      numNeighbors[iter->first] = iter->second.size();
      nodeIdPattern->setLen(subD.getSndChannel()[iter->first],iter->second.size());
      nodeIdPattern->setLen(subD.getRcvChannel()[iter->first],iter->second.size());
      nodeVolPattern->setLen(subD.getSndChannel()[iter->first],iter->second.size());
      nodeVolPattern->setLen(subD.getRcvChannel()[iter->first],iter->second.size());
      nodePosnPattern->setLen(subD.getSndChannel()[iter->first],3*iter->second.size());
      nodePosnPattern->setLen(subD.getRcvChannel()[iter->first],3*iter->second.size());
      nodeVecPattern->setLen(subD.getSndChannel()[iter->first],dim*iter->second.size());
      nodeVecPattern->setLen(subD.getRcvChannel()[iter->first],dim*iter->second.size());
      if (neq2 > 0) {
        nodeVecPatternEq1->setLen(subD.getSndChannel()[iter->first],neq1*iter->second.size());
        nodeVecPatternEq1->setLen(subD.getRcvChannel()[iter->first],neq1*iter->second.size());
        nodeVecPatternEq2->setLen(subD.getSndChannel()[iter->first],neq2*iter->second.size());
        nodeVecPatternEq2->setLen(subD.getRcvChannel()[iter->first],neq2*iter->second.size());
      }
      matPattern->setLen(subD.getSndChannel()[iter->first],neq1*neq1*iter->second.size());
      matPattern->setLen(subD.getRcvChannel()[iter->first],neq1*neq1*iter->second.size());
    }

    sharedNodes[iSub] = new Connectivity(subD.getNumNeighb(),numNeighbors);
    int* j = new int[subD.getNumNeighb()];
    memset(j,0,sizeof(int)*subD.getNumNeighb());
    for(std::map<int,std::map<int,int> >::const_iterator neighbor = newlyInsertedNodes.begin(); neighbor != newlyInsertedNodes.end(); ++neighbor) {
      for(std::map<int,int>::const_iterator iter = neighbor->second.begin(); iter != neighbor->second.end(); ++iter) {
        assert(subD.getNeighb()[neighbor->first] != subD.getGlobSubNum());
        assert(j[neighbor->first] < numNeighbors[neighbor->first]);
        (*sharedNodes[iSub])[neighbor->first][j[neighbor->first]++] = iter->second;
      }
    }
    for (int ii = 0; ii < subD.getNumNeighb(); ++ii)
      assert(j[ii] == numNeighbors[ii]);
  }
  nodeIdPattern->finalize();
  nodeVolPattern->finalize();
  nodeVecPattern->finalize();
  if (neq2 > 0) {
    nodeVecPatternEq1->finalize();
    nodeVecPatternEq2->finalize();
  }
  nodePosnPattern->finalize();
  matPattern->finalize();

  

#if 0
// #pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    fprintf(stderr,"SubDomain %d has... [%d nodes]\n", domain.getSubDomain()[iSub]->getGlobSubNum(), numNodes[iSub]);
    for(int jSub = 0; jSub < domain.getSubDomain()[iSub]->getNumNeighb(); ++jSub) {
      fprintf(stderr,"\tSubDomain %d as a neighbor with %d shared fine nodes; %d shared coarse nodes\n",
              domain.getSubDomain()[iSub]->getNeighb()[jSub],
              refinedSharedNodes[iSub]->num(jSub),
              sharedNodes[iSub]->num(jSub));
    }
  }
#endif
// ------------------------- END MPI PARALLEL CRAP -------------------------------------------------------------

  //std::cout << "gl node size " << globNodeIds.size() << std::endl;
  

  std::set< std::pair<int,int> > globalEdges;

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    // After agglomeration, allocate and build all the auxilliary data structures we'll be using
    int *numNodeNeighbors = new int[numNodes[iSub]];
    std::fill(numNodeNeighbors, numNodeNeighbors+numNodes[iSub], 1);

    int num_edges_skipped = 0;
    int (*edgePtr)[2] = refinedEdges[iSub]->getPtr();

    // Build the coarse EdgeSet, NodeToNode, and fine-to-coarse edge mapping
    edges[iSub] = new EdgeSet();
    int num_edges = 0;
    std::multimap<int,int> reverseEdgeMapping;
    for(int l = 0; l < refinedEdges[iSub]->size(); ++l) {
      int i = edgePtr[l][0], j = edgePtr[l][1];
      int new_i = nodeMapping(iSub)[i], new_j = nodeMapping(iSub)[j];
#if 1 // Debug output
      if(new_i < 0 || new_j < 0 || new_i > numNodes[iSub] || new_j > numNodes[iSub]) {
          fprintf(stderr,"ERROR; on SubD %d, i = %d, new_i = %d, j = %d, new_j = %d, but total num nodes = %d!\n",
            domain.getSubDomain()[iSub]->getGlobSubNum(), i, new_i, j, new_j, numNodes[iSub]);
          assert(false);
      }
#endif
      if(new_i != new_j) {

        int new_edge = edges[iSub]->find(new_i, new_j);
        edgeMapping(iSub)[l] = new_edge;
        reverseEdgeMapping.insert(std::pair<int,int>(new_edge,l));
        
        if(num_edges < edges[iSub]->size()) {
          numNodeNeighbors[new_i]++;
          numNodeNeighbors[new_j]++;
          num_edges = edges[iSub]->size();
        }
      } else {
        ++num_edges_skipped;
        edgeMapping(iSub)[l] = -1;
      }
    }

    Vec<int> newNum(edges[iSub]->size());
    edges[iSub]->createPointers(newNum);
    edges[iSub]->setMasterFlag(new bool[edges[iSub]->size()]);

//    Vec<int> tmp(edgeMapping(iSub));
    // Edges are renumbered by create pointers!  Update our edge mapping
    // accordingly

    for (std::multimap<int,int>::iterator itr = reverseEdgeMapping.begin();
         itr != reverseEdgeMapping.end(); ++itr) {

      edgeMapping(iSub)[(*itr).second] = newNum[(*itr).first];
    }
    
    std::vector<Face*> newFaces;
    int faceCount = 0;
    for (int i = 0; i < parent->faces[iSub]->size(); ++i) {

      bool ok;
      int newMap[4];
      Face& face((*parent->faces[iSub])[i]);
      switch (face.type()) {
        
      case Face::TRIA:
        ok = true;
        for (int k = 0; k < 3; ++k) {
          int new_node = nodeMapping(iSub)[face[k]];
          for (int j = 0; j < k; ++j) {
            if (new_node == newMap[j])
              ok = false;
          }
          newMap[k] = new_node;
        }
        if (ok) {
          Face* new_face = new(parent->faces[iSub]->getBlockAllocator()) FaceTria;
          (*faceMapping)(iSub)[i] = faceCount;
          new_face->setup(face.getCode(),newMap, faceCount++, face.getSurfaceID());
          for (int l = 0; l < 3; ++l) {
            new_face->setEdgeNum(l, edgeMapping(iSub)[face.getEdgeNum(l)]);
          }
          newFaces.push_back(new_face);
        }
      }
    }
    
    faces[iSub] = new FaceSet(newFaces.size());
    for (int i = 0; i < newFaces.size(); ++i)
      faces[iSub]->addFace(i, newFaces[i]); 

    faceNormDistInfo->setLen(iSub, faceCount);

    // Build the connectivity graph
    connectivity[iSub] = new Connectivity(numNodes[iSub], numNodeNeighbors);
    for(int i = 0; i < numNodes[iSub]; ++i) {
      (*connectivity[iSub])[i][0] = i;
      for(int n = 1; n < connectivity[iSub]->num(i);++n) (*connectivity[iSub])[i][n] = -1;
    }

    edgePtr = edges[iSub]->getPtr();
    for(int l = 0; l < edges[iSub]->size(); ++l) {
      int i=edgePtr[l][0], j=edgePtr[l][1];
      int index_i=0, index_j=0;
      while((*connectivity[iSub])[i][index_i] != -1) {
        ++index_i;
#if 1 // Debug output
        if(index_i >= connectivity[iSub]->num(i)) {
          fprintf(stderr,"ERROR: Trying to establish edge %d -> %d, but index_i = %d which is greater than %d\n",
                  i, j, index_i, connectivity[iSub]->num(i));
          assert(false);
        }
#endif
      }
      while((*connectivity[iSub])[j][index_j] != -1) {
          ++index_j;
#if 1 // Debug output
          if(index_j >= connectivity[iSub]->num(j)) {
          fprintf(stderr,"ERROR: Trying to establish edge %d -> %d, but index_j = %d which is greater than %d\n",
                  i, j, index_j, connectivity[iSub]->num(j));
          assert(false);
        }
#endif
      }
      (*connectivity[iSub])[i][index_i] = j;
      (*connectivity[iSub])[j][index_j] = i;
    }

    for (int l = 0; l < edges[iSub]->size(); ++l)
      edges[iSub]->getMasterFlag()[l] = true;


    /*fprintf(stderr, "\tAgglomeration finished with %d coarse cells (vs. %d fine cells) and %d edges vanished, %d coarse edges remaining (of %d)\n",
            numNodes[iSub], nodeMapping(iSub).size(), num_edges_skipped, edges[iSub]->size(), refinedEdges[iSub]->size());
*/
    nodeDistInfo->setLen(iSub, numNodes[iSub]);
    edgeDistInfo->setLen(iSub, edges[iSub]->size());
    faceDistInfo->setLen(iSub, faces[iSub]->size());
    inletNodeDistInfo->setLen(iSub, 0);
  }
//  fprintf(stderr, "\n");

  nodeDistInfo->finalize(true);
  faceDistInfo->finalize(true);
  faceNormDistInfo->finalize(true);
  inletNodeDistInfo->finalize(true);

  sharedEdges = new EdgeDef**[numLocSub];
  numSharedEdges = new int*[numLocSub];
  
  CommPattern<int> edgeNumPat(domain.getSubTopo(), domain.getCommunicator(), CommPattern<int>::CopyOnSend,
			      CommPattern<int>::NonSym);
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub)
    identifyEdges(edgeNumPat,iSub);

  edgeNumPat.finalize();
  
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub)
    sndEdgeInfo(edgeNumPat,iSub);

  edgeNumPat.exchange();

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    rcvEdgeInfo(edgeNumPat,iSub,dim);
  }

  edgeDistInfo->finalize(false);
  offDiagMatPattern->finalize();
  edgeAreaPattern->finalize();
  edgeVecPattern->finalize();
  
  myGeoState = new DistGeoState(parent->myGeoState->getGeoData(),&domain,
                                *nodeDistInfo,*edgeDistInfo,*faceNormDistInfo);

  edgeArea = new DistVec<double>(*edgeDistInfo);
  *edgeArea = 0.0;
  edgeNormals = &myGeoState->getEdgeNormal();//new DistVec<Vec3D>(*edgeDistInfo);  
  *edgeNormals = Vec3D(0.0);

  globalFaceNormals = new DistSVec<double,3>(*nodeDistInfo);
  *globalFaceNormals = 0.0;

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {


    (*myGeoState)(iSub).getEdgeNormalVel() = 0.0;//new DistVec<Vec3D>(*edgeDistInfo);  
    for(int l = 0; l < refinedEdges[iSub]->size(); ++l) {
      const int coarse_index = edgeMapping(iSub)[l];
      if(coarse_index < 0) continue;
      else {
        //(*edgeNormals)(iSub)[coarse_index] += refinedEdgeNormals(iSub)[l];
        if (refinedEdges[iSub]->getMasterFlag()[l]) {
          if (refinedEdgeAreas)
            (*edgeArea)(iSub)[coarse_index] += (*refinedEdgeAreas)(iSub)[l];
          else
            (*edgeArea)(iSub)[coarse_index] += refinedEdgeNormals(iSub)[l].norm();
          double sgn = 1.0;
          if (nodeMapping(iSub)[ refinedEdges[iSub]->getPtr()[l][0] ] >
              nodeMapping(iSub)[ refinedEdges[iSub]->getPtr()[l][1] ])
            sgn = -1.0;
          (*edgeNormals)(iSub)[coarse_index] += sgn*refinedEdgeNormals(iSub)[l];
        //edges[iSub]->viewEdgeLength()[coarse_index] = std::max(refinedEdges[iSub]->length(l), edges[iSub]->viewEdgeLength()[coarse_index]);
//        edges[iSub]->getMasterFlag()[coarse_index] = refinedEdges[iSub]->getMasterFlag()[l] &&
//                                                     edges[iSub]->getMasterFlag()[coarse_index];
        }
      }
    }
  } 

  ::assemble(domain, *edgeAreaPattern, numSharedEdges, sharedEdges, *edgeArea);
  ::assemble(domain, *edgeVecPattern, numSharedEdges, sharedEdges, *edgeNormals);

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {

    for (int i = 0; i < edges[iSub]->size(); ++i) {

      if (fabs((*edgeNormals)(iSub)[i].norm()) < 1.0e-15) {

        //std::cout << "Error: normal = 0!" << std::endl;
      }
    }
  }

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
   
    *nodeDistInfo->getMasterFlag(iSub) = true;
    for(int i = 0; i < nodeMapping(iSub).size(); ++i) {

      if (nodeMapping(iSub)[i] >= 0) {
        // Nodes
        nodeDistInfo->getMasterFlag(iSub)[nodeMapping(iSub)[i]] = refinedNodeDistInfo.getMasterFlag(iSub)[i];
  //      nodeDistInfo->getInvWeight(iSub)[nodeMapping(iSub)[i]] = min(nodeDistInfo->getInvWeight(iSub)[nodeMapping(iSub)[i]],
  //                                                                   refinedNodeDistInfo.getInvWeight(iSub)[i]);
      
      }
    }

  }

  finestNodeMapping = new DistVec<int>(*nodeDistInfo);
  MultiGridLevel* lvl = getFinestLevel();
  DistInfo& finestNodeInfo = lvl->getNodeDistInfo();
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    for (int i = 0; i < finestNodeInfo.subSize(iSub); ++i) {

      (*finestNodeMapping)(iSub)[mapFineToCoarse(iSub,i)] = i;
    }
  }
  
  //zero = new DistSVec<double,dim>(*nodeDistInfo);
  //*zero = 0.0;

  fv_comp_tag = new DistSVec<int,2>(*nodeDistInfo);
  *fv_comp_tag = 0; 

  ctrlVolCount = new DistVec<double>(*nodeDistInfo);
  ctrlVol = &myGeoState->getCtrlVol();//new DistVec<double>(*nodeDistInfo);
  *ctrlVol = 0.0;
  *ctrlVolCount = 0.0;
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    
    for(int i = 0; i < nodeMapping(iSub).size(); ++i) {

      if (nodeDistInfo->getMasterFlag(iSub)[nodeMapping(iSub)[i]]) {
        (*ctrlVolCount)(iSub)[nodeMapping(iSub)[i]] += 1.0;//(refinedVol)(iSub)[i];
        (*ctrlVol)(iSub)[nodeMapping(iSub)[i]] += (refinedVol)(iSub)[i];
      }
    }
  }

  operAdd<double> addOp;
  ::assemble(domain, *nodeVolPattern, sharedNodes, *ctrlVol, addOp);
  ::assemble(domain, *nodeVolPattern, sharedNodes, *ctrlVolCount, addOp);

  coarseNodeTopology = new DistVec<Topology>(*nodeDistInfo);

  coarseTopologyNormal = new DistVec<Vec3D>(*nodeDistInfo);
 
  Xn = &myGeoState->getXn();//new DistSVec<double,3>(*nodeDistInfo); 
  *Xn = 0.0;
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    SubDomain& subD(*domain.getSubDomain()[iSub]);
    PriorityNodes nodesOnSolidAndFront, nodesOnSolid,
                  nodesOnFarFieldAndFront, nodesOnFarField,
                  nodesOnSubDomainBoundaryAndFront, nodesOnSubDomainBoundary,
                  nodesOnFront;

    subD.getSolidBoundaryNodes(nodesOnSolid);
    subD.getFarFieldBoundaryNodes(nodesOnFarField);
    subD.getSubDomainBoundaryNodes(nodesOnSubDomainBoundary);

    mapNodeList(iSub,nodesOnSolid);
    mapNodeList(iSub,nodesOnFarField);
    mapNodeList(iSub,nodesOnSubDomainBoundary);
 
    std::set<int>* myBoundaryNodes = new std::set<int>[(*Xn)(iSub).size()];
    Topology* newTopo = &(*coarseNodeTopology)(iSub)[0];
    memset(newTopo,0,sizeof(Topology)*(*Xn)(iSub).size());
    for (int i = 0; i < (*Xn)(iSub).size(); ++i)
      newTopo[i] = TopoInterior;
    double* coarseVol = new double[(*Xn)(iSub).size()];
    memset(coarseVol,0,sizeof(double)*(*Xn)(iSub).size());
    for(int i = 0; i < refinedVol(iSub).size(); ++i) {
      const int coarseIndex = nodeMapping(iSub)[i];
      if (!refinedNodeDistInfo.getMasterFlag(iSub)[i])
        continue;
      if (nodesOnSolid.find(i) != nodesOnSolid.end() ||
          nodesOnFarField.find(i) != nodesOnFarField.end()) {
        myBoundaryNodes[coarseIndex].insert(i);
      }
      if ((*nodeTopology)(iSub)[i] == TopoFace) {
        newTopo[coarseIndex] = TopoFace; 
        (*coarseTopologyNormal)(iSub)[coarseIndex] = (*topologyNormal)(iSub)[i];
      }
      else if ((*nodeTopology)(iSub)[i] == TopoVertex) {
        newTopo[coarseIndex] = TopoVertex;
        (*coarseTopologyNormal)(iSub)[coarseIndex] = (*topologyNormal)(iSub)[i];
      }
      else if ((*nodeTopology)(iSub)[i] == TopoLine) {
        newTopo[coarseIndex] = TopoLine;
        (*coarseTopologyNormal)(iSub)[coarseIndex] = (*topologyNormal)(iSub)[i];
      }
    }

    for(int i = 0; i < refinedVol(iSub).size(); ++i) {
      const int coarseIndex = nodeMapping(iSub)[i];
      if (!refinedNodeDistInfo.getMasterFlag(iSub)[i])
        continue;
      if (newTopo[coarseIndex] == TopoFace &&
          (*nodeTopology)(iSub)[i] != TopoFace)
        continue;
      if (myBoundaryNodes[coarseIndex].empty() ||
          myBoundaryNodes[coarseIndex].find(i) != myBoundaryNodes[coarseIndex].end()) {
        coarseVol[coarseIndex] += (refinedVol)(iSub)[i];
        for(int j = 0; j < 3; ++j)
          (*Xn)(iSub)[coarseIndex][j] += (refinedVol)(iSub)[i] * (*parent->Xn)(iSub)[i][j];
      }
    }

    for(int i = 0; i < (*Xn)(iSub).size(); ++i) {
      if (!nodeDistInfo->getMasterFlag(iSub)[i])
        continue;
      const Scalar one_over_volume = 1.0/coarseVol[i];
      for(int j = 0; j < 3; ++j) (*Xn)(iSub)[i][j] *= one_over_volume;
    }
    delete [] myBoundaryNodes;
    delete [] coarseVol;
  }

  operAdd<double> addOp2;
  ::assemble(domain, *nodePosnPattern, sharedNodes, *Xn, addOp2);
    
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    for(int l = 0; l < edges[iSub]->size(); ++l) {

      int i = edges[iSub]->getPtr()[l][0], j = edges[iSub]->getPtr()[l][1];
      double v[3];
      for (int k = 0; k <3; ++k) 
        v[k] = (*Xn)(iSub)[i][k]-(*Xn)(iSub)[j][k];
      edges[iSub]->viewEdgeLength()[l] = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
     // assert(edge_areas(iSub)[l] > 0);
      //for (int k = 0; k <3; ++k)
      //  (*edgeNormals)(iSub)[l][k] = -edge_areas(iSub)[l]*v[k] / edges[iSub]->viewEdgeLength()[l]; 
    }
  }

   
#pragma omp parallel for 
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    for (int i = 0; i < faces[iSub]->size(); ++i) {

      (*faces[iSub])[i].computeNormal((*Xn)(iSub), myGeoState->getFaceNormal()(iSub)[i]);
    
    }
  }

  myGeoState->getFaceNorVel() = 0.0;

  agglomeratedFaces = new AgglomeratedFaceSet*[numLocSub]; 
  // Construct the agglomerated faces
  if (!parent->agglomeratedFaces) {

    typedef std::map<std::pair<int,int>, AgglomeratedFace> AggFaces;
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      AggFaces myNewFaces;
      for (int i = 0; i < parent->faces[iSub]->size(); ++i) {

        Face& F = (*parent->faces[iSub])[i];
        int code = F.getCode();
        for (int j = 0; j < F.numNodes(); ++j) {

          int node = nodeMapping(iSub)[F[j]];
          AggFaces::iterator itr = myNewFaces.find(std::pair<int,int>(node,code));
          Vec3D n = F.getNormal(parent->myGeoState->getFaceNormal()(iSub),j);
          for (int k = 0; k < 3; ++k)
            (*globalFaceNormals)(iSub)[node][k] += n[k];
          if (itr != myNewFaces.end()) {

            (itr->second).getNormal() += n;
            (itr->second).getArea() += n.norm();
          } else {

            AgglomeratedFace A(node,code);
            A.getArea() += n.norm();
            A.getNormal() += n;
            myNewFaces[std::pair<int,int>(node,code)] = A;
          }
        }
      }

      agglomeratedFaces[iSub] = new AgglomeratedFaceSet(myNewFaces.size());
      agglomFaceDistInfo->setLen(iSub, myNewFaces.size());
      int i = 0;
      for (AggFaces::iterator it = myNewFaces.begin();
           it != myNewFaces.end(); ++it) {

        agglomeratedFaces[iSub]->operator[](i) = it->second;
        ++i;
      }
    }
  } else {
    typedef std::map<std::pair<int,int>, AgglomeratedFace> AggFaces;
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      AggFaces myNewFaces;
      for (int i = 0; i < parent->agglomeratedFaces[iSub]->size(); ++i) {

        AgglomeratedFace& F = (*parent->agglomeratedFaces[iSub])[i];
        int code = F.getCode();

        int node = nodeMapping(iSub)[F.getNode()];
        AggFaces::iterator itr = myNewFaces.find(std::pair<int,int>(node,code));
        Vec3D n = F.getNormal();
        for (int k = 0; k < 3; ++k) 
          (*globalFaceNormals)(iSub)[node][k] += n[k];
        if (itr != myNewFaces.end()) {

          (itr->second).getNormal() += n;
          (itr->second).getArea() += F.getArea();
        } else {

          AgglomeratedFace A(node,code);
          A.getNormal() += n;
          A.getArea() += F.getArea();
          myNewFaces[std::pair<int,int>(node,code)] = A;
        }
      }

      agglomeratedFaces[iSub] = new AgglomeratedFaceSet(myNewFaces.size());
      agglomFaceDistInfo->setLen(iSub, myNewFaces.size());
      int i = 0;
      for (AggFaces::iterator it = myNewFaces.begin();
           it != myNewFaces.end(); ++it) {

        agglomeratedFaces[iSub]->operator[](i) = it->second;
        ++i;
      }
    } 
  }

  agglomFaceDistInfo->finalize(true);
  
  ::assemble(domain, *nodePosnPattern, sharedNodes, *globalFaceNormals, addOp2);
   
  setNodeType(); 
}

template<class Scalar>
template<int dim>
void MultiGridLevel<Scalar>::setupBcs(DistBcData<dim>& fineBcData, DistBcData<dim>& coarseBcData,
                                      DistSVec<Scalar,dim>& boundaryState) {

  Restrict(*parent, fineBcData.getUfarin(), coarseBcData.getUfarin());
  Restrict(*parent, fineBcData.getUfarout(), coarseBcData.getUfarout());
  Restrict(*parent, fineBcData.getNodeStateVector(), 
           coarseBcData.getNodeStateVector());

  DistVec<double>* d2wall_coarse = myGeoState->getd2wall(),
   *d2wall_fine = parent->myGeoState->getd2wall();

  Restrict(*parent, *d2wall_fine, *d2wall_coarse);

#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {

/*    for (int i=0; i<faces[iSub]->size(); ++i) 
      (*faces[iSub])[i].template assignFreeStreamValues2<dim>(coarseBcData.getUfarin(),
                                                     coarseBcData.getUfarout() ,
                                                     coarseBcData.getFaceStateVector()(iSub)[i]);
*/
    for (int i=0; i<agglomeratedFaces[iSub]->size(); ++i) 
      (*agglomeratedFaces[iSub])[i].template assignFreeStreamValues2<dim>(coarseBcData.getUfarin(),
                                                     coarseBcData.getUfarout() ,
                                                     boundaryState(iSub)[i]);
 
  }  
}

template<class Scalar>
void MultiGridLevel<Scalar>::setNodeType() {

  int* bcpriority = reinterpret_cast<int *>(alloca(sizeof(int) * (BC_MAX_CODE - BC_MIN_CODE + 1)));
  bcpriority -= BC_MIN_CODE;

  bcpriority[BC_ADIABATIC_WALL_MOVING     ] = 20;
  bcpriority[BC_ISOTHERMAL_WALL_MOVING    ] = 19;
  bcpriority[BC_SLIP_WALL_MOVING          ] = 18;
  bcpriority[BC_POROUS_WALL_MOVING        ] = 17;
  bcpriority[BC_MASSFLOW_INLET_MOVING     ] = 16;
  bcpriority[BC_MASSFLOW_OUTLET_MOVING    ] = 15;
  bcpriority[BC_DIRECTSTATE_INLET_MOVING  ] = 14;
  bcpriority[BC_DIRECTSTATE_OUTLET_MOVING ] = 13;
  bcpriority[BC_INLET_MOVING              ] = 12;
  bcpriority[BC_OUTLET_MOVING             ] = 11;
  bcpriority[BC_ADIABATIC_WALL_FIXED      ] = 10;
  bcpriority[BC_ISOTHERMAL_WALL_FIXED     ] =  9;
  bcpriority[BC_SLIP_WALL_FIXED           ] =  8;
  bcpriority[BC_POROUS_WALL_FIXED         ] =  7;
  bcpriority[BC_MASSFLOW_INLET_FIXED      ] =  6;
  bcpriority[BC_MASSFLOW_OUTLET_FIXED     ] =  5;
  bcpriority[BC_DIRECTSTATE_INLET_FIXED   ] =  4;
  bcpriority[BC_DIRECTSTATE_OUTLET_FIXED  ] =  3;
  bcpriority[BC_INLET_FIXED               ] =  2;
  bcpriority[BC_OUTLET_FIXED              ] =  1;
  bcpriority[BC_SYMMETRY                  ] =  0;
  bcpriority[BC_INTERNAL                  ] = -1;

  nodeType = new int*[numLocSub];

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    
    nodeType[iSub] = new int[nodeDistInfo->subSize(iSub)];
    for (int i = 0; i < nodeDistInfo->subSize(iSub); ++i) {
      nodeType[iSub][i] = BC_INTERNAL;
    }
    
    for (int i = 0; i < agglomeratedFaces[iSub]->size(); ++i)
      (*agglomeratedFaces[iSub])[i].setNodeType(bcpriority, nodeType[iSub]);

    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      SubRecInfo<int> nInfo = nodeIdPattern->getSendBuffer(subD.getSndChannel()[jSub]);
      for (int i = 0; i < sharedNodes[iSub]->num(jSub); ++i)
        nInfo.data[i] = nodeType[iSub][ (*sharedNodes[iSub])[jSub][i] ];
    }
 
  }

  nodeIdPattern->exchange();
   
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int jSub = 0; jSub < subD.getNumNeighb(); ++jSub) {
      SubRecInfo<int> nInfo = nodeIdPattern->recData(subD.getRcvChannel()[jSub]);
      for (int i = 0; i < sharedNodes[iSub]->num(jSub); ++i)
        if (bcpriority[ nInfo.data[i] ] > bcpriority[ nodeType[iSub][ (*sharedNodes[iSub])[jSub][i] ] ])
          nodeType[iSub][ (*sharedNodes[iSub])[jSub][i] ] = nInfo.data[i];
    }

  }   

}

template<class Scalar>
void MultiGridLevel<Scalar>::computeNodeNormalClasses() {

  if (parent->coarseNodeTopology) {
    nodeTopology = parent->coarseNodeTopology;
    topologyNormal = new DistVec<Vec3D>(*parent->nodeDistInfo);
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {

      for (int i = 0; i < (*parent->pNodeMapping)(iSub).size(); ++i) {

        if ((*parent->nodeTopology)(iSub)[i] != TopoInterior) {

          (*topologyNormal)(iSub)[(*parent->pNodeMapping)(iSub)[i]] = 
            (*parent->topologyNormal)(iSub)[i];
        }
      }
    }
    return;
  }

  nodeNormalCount = new DistVec<int>(parent->getNodeDistInfo());
  *nodeNormalCount = 0;

#pragma omp parallel for 
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    for (int i = 0; i < parent->faces[iSub]->size(); ++i) {

      Face& F = (*parent->faces[iSub])[i];
      for (int j = 0; j < F.numNodes(); ++j) {

        int nd = F[j];
        ++(*nodeNormalCount)(iSub)[nd];
      }
    }
  }

  operAdd<int> addOp;
  ::assemble(domain, *parent->nodeIdPattern, parent->sharedNodes, 
             *nodeNormalCount, addOp);
  
  int maxNodeAdj = nodeNormalCount->max();
  nodeNormalsPattern = new CommPattern<double>(domain.getSubTopo(), domain.getCommunicator(), CommPattern<double>::CopyOnSend);

#pragma omp parallel for 
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    SubDomain& subD(*domain.getSubDomain()[iSub]);
    for (int ch = 0; ch < parent->nodeVolPattern->numCh(); ++ch) {

      nodeNormalsPattern->setLen(ch,parent->nodeVolPattern->getLen(ch) *3*maxNodeAdj);
    }
  }

  nodeNormalsPattern->finalize();

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  nodeNormals = new std::list<Vec3D>*[numLocSub];
#pragma omp parallel for 
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    nodeNormals[iSub] = new std::list<Vec3D>[parent->getXn()(iSub).size()];
    for (int i = 0; i < parent->faces[iSub]->size(); ++i) {

      Face& F = (*parent->faces[iSub])[i];
      Vec3D N;
      F.computeNormal(parent->getXn()(iSub),N);
      for (int j = 0; j < F.numNodes(); ++j) {

        int nd = F[j];
        nodeNormals[iSub][nd].push_back(N);
      }
    }
  }
 
  ::assemble(domain, *nodeNormalsPattern,  parent->sharedNodes,
             nodeNormals,maxNodeAdj);

  double costh_cutoff = mesh_topology_threshold;

  nodeTopology = new DistVec<Topology>(*parent->nodeDistInfo);
  topologyNormal = new DistVec<Vec3D>(*parent->nodeDistInfo);
#pragma omp parallel for 
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    
    Vec<Topology>& T = (*nodeTopology)(iSub);
    for (int i = 0; i < T.size(); ++i) {

      std::list<Vec3D>& L = nodeNormals[iSub][i];
      int lsize = L.size();
      if (lsize == 0) {
        T[i] = TopoInterior;
        continue;
      }
      std::list<Vec3D> class_list;
      
      for (std::list<Vec3D>::iterator itr = L.begin(); itr != L.end(); 
           ++itr) {
        Vec3D& a = *itr;
        bool classified = false;
        for (std::list<Vec3D>::iterator itrc = class_list.begin(); 
             itrc != class_list.end(); ++itrc) {
 
          Vec3D& b = *itrc;
          if (a*b / (a.norm()*b.norm()) > costh_cutoff) {

            classified = true;
            break;
          }
        }
        if (!classified)
          class_list.push_back(a);
      }
      if (class_list.size() == 2 )
        T[i] = TopoLine;
      else if (class_list.size() > 2)
        T[i] = TopoVertex;
      else
        T[i] = TopoFace;

      if (T[i] == TopoLine) {
        Vec3D& a = *class_list.begin();
        Vec3D& b = *++class_list.begin();
        Vec3D no = a^b;
        no /= no.norm();
        (*topologyNormal)(iSub)[i] = no;
      } else if (T[i] == TopoFace) {
        Vec3D no(0,0,0);
        for (std::list<Vec3D>::iterator itr = L.begin(); itr != L.end();
             ++itr) {
          
          no += *itr;
        }
        no /= no.norm();
        (*topologyNormal)(iSub)[i] = no; 
      }
    }
  }

  delete nodeNormalsPattern;

  for (int iSub = 0; iSub < numLocSub; ++iSub)
    delete [] nodeNormals[iSub];
 
  delete [] nodeNormals;
}

//------------------------------------------------------------------------------

template<class Scalar>
void MultiGridLevel<Scalar>::computeRestrictedQuantities(const DistGeoState& refinedGeoState)
{
  /* 
  distGeoState->getXn() *= 0.0;
  distGeoState->getCtrlVol() *= 0.0;
  const DistVec<double>& refinedVolume(refinedGeoState.getCtrlVol());
  const DistSVec<double, 3>& refinedX(refinedGeoState.getXn());

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    for(int i = 0; i < refinedVolume(iSub).size(); ++i) {
      const int coarseIndex = nodeMapping(iSub)[i];
      distGeoState->getCtrlVol()(iSub)[coarseIndex] += refinedVolume(iSub)[i];
      for(int j = 0; j < 3; ++j)
        distGeoState->getXn()(iSub)[coarseIndex][j] += refinedVolume(iSub)[i] * refinedX(iSub)[i][j];
    }

    for(int i = 0; i < distGeoState->getXn()(iSub).size(); ++i) {
      const Scalar one_over_volume = 1.0/distGeoState->getCtrlVol()(iSub)[i];
      for(int j = 0; j < 3; ++j) distGeoState->getXn()(iSub)[i][j] *= one_over_volume;
    }
  }

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub)
    for(int i = 0; i < distGeoState->getXn()(iSub).size(); ++i)
      if(!nodeDistInfo->getMasterFlag(iSub)[i]) {
        for(int j = 0; j < 3; ++j) distGeoState->getXn()(iSub)[i][j] = 0.0;
        distGeoState->getCtrlVol()(iSub)[i] = 0.0;
      }

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    (*distGeoState)(iSub).getEdgeNormal() = 0.0;
    for(int l = 0; l < edgeMapping(iSub).size(); ++l) {
      const int coarse_l = edgeMapping(iSub)[l];
      if(coarse_l < 0) continue;
      else {
        (*distGeoState)(iSub).getEdgeNormal()[coarse_l] += refinedGeoState(iSub).getEdgeNormal()[l];
      }
    }
  }
  operAdd<double> addOp;
  ::assemble(domain, *nodePosnPattern, sharedNodes, distGeoState->getXn(), addOp);
  ::assemble(domain, *nodeVolPattern, sharedNodes, distGeoState->getCtrlVol(), addOp);

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    for (int i = 0; i < faces[iSub]->size(); ++i) {

      (*faces[iSub])[i].computeNormal(distGeoState->getXn()(iSub), distGeoState->getFaceNormal()(iSub));
    }
  }
*/
}

//------------------------------------------------------------------------------

template<class Scalar> template<class Scalar2, int dim>
void MultiGridLevel<Scalar>::Restrict(const MultiGridLevel<Scalar>& fineGrid, const DistSVec<Scalar2, dim>& fineData, DistSVec<Scalar2, dim>& coarseData,bool average,
                                      bool apply_relaxation) const
{
  coarseData = 0.0;

  int rnk;
  MPI_Comm_rank(MPI_COMM_WORLD,&rnk);
/*  if (rnk == 110) {

    std::cout << "Calling restrict --- " << std::endl;
  }
*/

  DistVec<int>& nodeMapping = *pNodeMapping;
  if (average) {
    if (useVolumeWeightedAverage) {
#pragma omp parallel for
      for(int iSub = 0; iSub < numLocSub; ++iSub) {
	for(int i = 0; i < fineData(iSub).size(); ++i)  {
	  /*       if (fineGrid.getMyLevel() == 0 && iSub == 0 && rnk == 110 && 
		   nodeMapping(iSub)[i] == 3136)
		   std::cout << "f = [" << fineData(iSub)[i][0] << " " <<
		   fineData(iSub)[i][1] << " " <<
		   fineData(iSub)[i][2] << " " <<
		   fineData(iSub)[i][3] << " " <<
		   fineData(iSub)[i][4] << "]" << std::endl;
	  */
	  if (!fineGrid.nodeDistInfo->getMasterFlag(iSub)[i])
	    continue;
	  /*if (rnk == 7 && nodeMapping(iSub)[i] == 2) {
	    
	    std::cout << "i = " << i << ", vol = " << fineGrid.getCtrlVol().v[i] << std::endl;
	    for(int j = 0; j < dim; ++j)
            std::cout << "V[" << j << "] = " << fineData.v[i][j] << std::endl;
	    
	    }*/
	  
	  double loc_relax_factor = 1.0;
	  if (mgSubdomains[iSub].nodeTopology[nodeMapping(iSub)[i]] != 
	      MultigridSubdomain::TopoInterior && apply_relaxation)
	    loc_relax_factor *= 0.5;

	  for(int j = 0; j < dim; ++j) {
	    //if (j < 5) {
	      coarseData(iSub)[nodeMapping(iSub)[i]][j] += loc_relax_factor*fineGrid.getCtrlVol()(iSub)[i] * fineData(iSub)[i][j];
	      //} else {
	      //  coarseData(iSub)[nodeMapping(iSub)[i]][j] += loc_relax_factor*fineGrid.getCtrlVol()(iSub)[i] * std::max(fineData(iSub)[i][j],0.0);
	      //}
	  }
	}
	
	for(int i = 0; i < coarseData(iSub).size(); ++i) {
	  const Scalar one_over_volume = 1.0 / getCtrlVol()(iSub)[i];
	  for(int j = 0; j < dim; ++j) coarseData(iSub)[i][j] *= one_over_volume;
	}
	
      }
    } else {
      
      assert(false);
#pragma omp parallel for
      for(int iSub = 0; iSub < numLocSub; ++iSub) {
	for(int i = 0; i < fineData(iSub).size(); ++i)  {
	  for(int j = 0; j < dim; ++j)
	    coarseData(iSub)[nodeMapping(iSub)[i]][j] += /*(*fineGrid.ctrlVolCount)(iSub)[i] */ fineData(iSub)[i][j];
	}
	
	for(int i = 0; i < coarseData(iSub).size(); ++i) {
	  const Scalar one_over_volume = 1.0 / (*ctrlVolCount)(iSub)[i];
	  if (nodeDistInfo->getMasterFlag(iSub)[i]) {
	    for(int j = 0; j < dim; ++j) coarseData(iSub)[i][j] *= one_over_volume;
	  } else {
	    for(int j = 0; j < dim; ++j) coarseData(iSub)[i][j] = 0.0;
	  }      
	}
	
      }
      
    }
  } else {

#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub) {
      for(int i = 0; i < fineData(iSub).size(); ++i)  {
	if (!fineGrid.nodeDistInfo->getMasterFlag(iSub)[i])
	  continue;
	
	double loc_relax_factor = 1.0;
	if (mgSubdomains[iSub].nodeTopology[nodeMapping(iSub)[i]] != 
	    MultigridSubdomain::TopoInterior && apply_relaxation)
	  loc_relax_factor *= 0.5;
	
	for(int j = 0; j < dim; ++j)
	  coarseData(iSub)[nodeMapping(iSub)[i]][j] += loc_relax_factor*fineData(iSub)[i][j];
      }
      
    }
  }
 
  operAdd<double> addOp;
  if (agglomType == AgglomerationLocal)
    ::assemble(domain, *nodeVecPattern, const_cast<Connectivity **>(sharedNodes), coarseData, addOp);
  else
    assembleInternal(coarseData);
  
  //std::cout << "Sum = " << sum << std::endl;
}

template<class Scalar> template<class Scalar2>
void MultiGridLevel<Scalar>::Restrict(const MultiGridLevel<Scalar>& fineGrid, const DistVec<Scalar2>& fineData, DistVec<Scalar2>& coarseData) const
{
  coarseData = 0.0;

  DistVec<int>& nodeMapping = *pNodeMapping;

  if (useVolumeWeightedAverage) {
#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub) {
      for(int i = 0; i < fineData(iSub).size(); ++i) {

        if (!fineGrid.nodeDistInfo->getMasterFlag(iSub)[i])
          continue;
        coarseData(iSub)[nodeMapping(iSub)[i]] += /*fineGrid.getCtrlVol()(iSub)[i] * */fineData(iSub)[i];
      }
      /*
      for(int i = 0; i < coarseData(iSub).size(); ++i) {
        const Scalar one_over_volume = 1.0 / getCtrlVol()(iSub)[i];
        coarseData(iSub)[i] *= one_over_volume;
      }
      */
    }
  } else {
    assert(false);
#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub) {
      for(int i = 0; i < fineData(iSub).size(); ++i)
        coarseData(iSub)[nodeMapping(iSub)[i]] += /*(*fineGrid.ctrlVolCount)(iSub)[i] */ fineData(iSub)[i];

      for(int i = 0; i < coarseData(iSub).size(); ++i) {
        const Scalar one_over_volume = 1.0 / (*ctrlVolCount)(iSub)[i];
        if (nodeDistInfo->getMasterFlag(iSub)[i]) {
          coarseData(iSub)[i] *= one_over_volume;
        } else {
          coarseData(iSub)[i] = 0.0;
        }      
      }

    }
 
  }
  
  operAdd<double> addOp;
  if (agglomType == AgglomerationLocal)
    ::assemble(domain, *nodeVolPattern, const_cast<Connectivity **>(sharedNodes), coarseData, addOp);
  else
    assembleInternal(coarseData);
}

template<class Scalar> template<class Scalar2, int dim>
void MultiGridLevel<Scalar>::RestrictFaceVector(const MultiGridLevel<Scalar>& fineGrid, const DistSVec<Scalar2, dim>& fineData, DistSVec<Scalar2, dim>& coarseData) const
{
  coarseData = 0.0;
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    for(int i = 0; i < fineData(iSub).size(); ++i)  {
      int fidx = (*faceMapping)(iSub)[i];
      if (fidx < 0) continue;
      for(int j = 0; j < dim; ++j)
        coarseData(iSub)[fidx][j] = fineData(iSub)[i][j];
    }
  }

}


//------------------------------------------------------------------------------

template<class Scalar> template<class Scalar2, int dim>
void MultiGridLevel<Scalar>::RestrictOperator(const MultiGridLevel<Scalar>& fineGrid,
                                              DistMat<Scalar2,dim>& fineOperator,
                                              DistMat<Scalar2,dim>& coarseOperator) {

  coarseOperator = 0.0;

  DistVec<int>& nodeMapping = *pNodeMapping;
  DistVec<int>& edgeMapping = *pEdgeMapping;

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    bool* masterFlag = fineGrid.nodeDistInfo->getMasterFlag(iSub);
    for(int i = 0; i < fineGrid.nodeDistInfo->subSize(iSub); ++i) {

      if (!masterFlag[i])
        continue;

      int m = nodeMapping(iSub)[i];
      Scalar2* Aii = fineOperator(iSub).getElem_ii(i);
      Scalar2* Aii_coarse = coarseOperator(iSub).getElem_ii(m);
      Scalar2 rat = fineGrid.getCtrlVol()(iSub)[i] / 
                   getCtrlVol()(iSub)[m];
      for (int k = 0; k < dim*dim; ++k) {
     
        Aii_coarse[k] += Aii[k] * rat /* rat*/;
      }
    }

    int (*edgePtr)[2] = fineGrid.edges[iSub]->getPtr();
    int (*edgePtrCoarse)[2] = edges[iSub]->getPtr();

    masterFlag = fineGrid.edges[iSub]->getMasterFlag();    
    for (int l = 0; l < fineGrid.edges[iSub]->size(); ++l) {

      if (!masterFlag[l])
        continue;

      int i = edgePtr[l][0],j = edgePtr[l][1];
      int ci = nodeMapping(iSub)[i], cj = nodeMapping(iSub)[j];
      Scalar2* Aij = fineOperator(iSub).getElem_ij(l),
            * Aji = fineOperator(iSub).getElem_ji(l);
      Scalar2 rati = fineGrid.getCtrlVol()(iSub)[i] /
                     getCtrlVol()(iSub)[ci],
             ratj = fineGrid.getCtrlVol()(iSub)[j] /
                    getCtrlVol()(iSub)[cj];
      if (ci == cj) {
 
        Scalar2* Aii_coarse = coarseOperator(iSub).getElem_ii(ci);
        for (int k = 0; k < dim*dim; ++k) {

          Aii_coarse[k] += Aij[k]*rati/*rati*/+Aji[k]*ratj/*ratj*/;
        }
      } else {

        int cl = edgeMapping(iSub)[l];
        assert(cl >= 0);
        Scalar2* Aij_coarse = coarseOperator(iSub).getElem_ij(cl);
        Scalar2* Aji_coarse = coarseOperator(iSub).getElem_ji(cl);
        Scalar2* tmp;
        if (edgePtrCoarse[cl][0] != ci) {

          assert(edgePtrCoarse[cl][1] == ci);
          tmp = Aij_coarse; Aij_coarse = Aji_coarse; Aji_coarse = tmp;
        }
        for (int k = 0; k < dim*dim; ++k) { 
          Aij_coarse[k] += Aij[k]*rati/*ratj*/;
          Aji_coarse[k] += Aji[k]*ratj/*rati*/;
        }
      }
    }
  }

  ::assemble(domain,*offDiagMatPattern, numSharedEdges,sharedEdges,coarseOperator);
  ::assemble(domain, *matPattern, sharedNodes, coarseOperator);
}                                              

template<class Scalar> template<class Scalar2, int dim>
void MultiGridLevel<Scalar>::Prolong(MultiGridLevel<Scalar>& fineGrid, const DistSVec<Scalar2,dim>& coarseInitialData,
                                     const DistSVec<Scalar2,dim>& coarseData, DistSVec<Scalar2,dim>& fineData,
				     DistSVec<Scalar2,dim>& fine_ref,double relax_factor,VarFcn* varFcn,
				     DistLevelSetStructure* coarseLSS, DistLevelSetStructure* fineLSS) const
{

  int rnk;
  MPI_Comm_rank(MPI_COMM_WORLD,&rnk);

  DistVec<int>& nodeMapping = *pNodeMapping;
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    
    LevelSetStructure* coarseLSS_sub = (coarseLSS ? &(*coarseLSS)(iSub) : NULL);
    LevelSetStructure* fineLSS_sub = (fineLSS ? &(*fineLSS)(iSub) : NULL);
    for(int i = 0; i < fineData(iSub).size(); ++i) {
      const int coarseIndex = nodeMapping(iSub)[i];
      //if (coarseLSS_sub && !(coarseLSS_sub->isActive(0.0, coarseIndex)))
      //	  continue;
      if (fineLSS_sub && !(fineLSS_sub->isActive(0.0, i)))
      	  continue;
      //if (coarseLSS_sub && !(coarseLSS_sub->isActive(0.0, coarseIndex)))
      //	std::cout << coarseData(iSub)[coarseIndex][0] << " " <<  coarseInitialData(iSub)[coarseIndex][0] << std::endl;

      if (mgMethod == MultiGridAlgebraic) {
        for(int j = 0; j < dim; ++j) {
          fineData(iSub)[i][j] += fineGrid.getCtrlVol()(iSub)[i]*(coarseData(iSub)[coarseIndex][j] - coarseInitialData(iSub)[coarseIndex][j]) / getCtrlVol()(iSub)[coarseIndex];
        }
      } else {
        /*if (iSub == 0 && rnk == 110 && (coarseIndex == 3136 || coarseIndex == 2426 || coarseIndex == 2428 || coarseIndex == 3502 || coarseIndex == 3612)) {
          std::cout << "Pl[" << coarseIndex << "] = ";
          for(int j = 0; j < dim; ++j)
            std::cout << coarseData(iSub)[coarseIndex][j] - coarseInitialData(iSub)[coarseIndex][j] << " ";
          std::cout << std::endl;
        }
          */
        double loc_relax_factor = relax_factor;
        //if (mgSubdomains[iSub].nodeTopology[coarseIndex] != 
        //    MultigridSubdomain::TopoInterior)
        //  loc_relax_factor *= 0.5;
	if (dim > 5 && fine_ref(iSub)[i][5] > turbRelaxCutoff) {

	  loc_relax_factor *= turbRelaxCutoff/ fine_ref(iSub)[i][5];
	  //std::cout << "loc_relax_factor  = " << loc_relax_factor << " " << fineData(iSub)[i][5] << std::endl;
	}

	/*if (mgSubdomains[iSub].nodeTopology[nodeMapping(iSub)[i]] ==
	    MultigridSubdomain::TopoLine ||
	    mgSubdomains[iSub].nodeTopology[nodeMapping(iSub)[i]] ==
	    MultigridSubdomain::TopoVertex)
	  loc_relax_factor = 0.0;
	*/
        double tmp[dim], tmpv[dim];
        for(int j = 0; j < dim; ++j) {
          tmp[j] = fineData(iSub)[i][j] + loc_relax_factor*(coarseData(iSub)[coarseIndex][j] - coarseInitialData(iSub)[coarseIndex][j]);
        }
        if (varFcn) {
          varFcn->conservativeToPrimitive(tmp, tmpv);
          if (tmpv[0] <= 0.0 || tmpv[4] <= 0.0) {
            memcpy(tmp, fineData(iSub)[i], sizeof(tmp));
          }
        }
        memcpy(fineData(iSub)[i], tmp, sizeof(tmp));

      }
    }
  }
  /*
  if (coarseLSS) {

    
  }*/
}

template<class Scalar> template<class Scalar2, int dim>
void MultiGridLevel<Scalar>::ExtrapolateProlongation(MultiGridLevel<Scalar>& fineGrid, 
						     const DistSVec<Scalar2,dim>& coarseInitialData,
						     DistSVec<Scalar2,dim>& coarseData,
						     DistLevelSetStructure* coarseLSS, 
						     DistLevelSetStructure* fineLSS) {

  DistVec<int> cnt(coarseData.info());
  cnt = 0;
  DistSVec<Scalar2,dim> extUpdate(coarseData);
  extUpdate = 0.0;

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    
    LevelSetStructure* coarseLSS_sub = (coarseLSS ? &(*coarseLSS)(iSub) : NULL);

    for (int i = 0; i < cnt.subSize(iSub); ++i) {

      bool iActive = coarseLSS_sub->isActive(0.0, i);

      if (!iActive) {
	for (int k = 0; k < dim; ++k)
	  coarseData(iSub)[i][k] = coarseInitialData(iSub)[i][k];
      }
    }
  }

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    
    LevelSetStructure* coarseLSS_sub = (coarseLSS ? &(*coarseLSS)(iSub) : NULL);
    bool* edgeFlag = edges[iSub]->getMasterFlag();
    int (*edgePtr)[2] = edges[iSub]->getPtr();
    for (int l = 0; l < edges[iSub]->size(); ++l) {

      if (!edgeFlag[l])
        continue;

      int i = edgePtr[l][0];
      int j = edgePtr[l][1];

      bool iActive = coarseLSS_sub->isActive(0.0, i),
	jActive = coarseLSS_sub->isActive(0.0, j);

      if ((iActive && jActive) ||
	  (!iActive && !jActive))
	continue;

      if (iActive) {

	++cnt(iSub)[j];
	for (int k = 0; k < dim; ++k)
	  extUpdate(iSub)[j][k] += (coarseData(iSub)[i][k]-coarseInitialData(iSub)[i][k]);
      } else {

	++cnt(iSub)[i];
	for (int k = 0; k < dim; ++k)
	  extUpdate(iSub)[i][k] += (coarseData(iSub)[j][k]-coarseInitialData(iSub)[j][k]);

      }    
      
    }
  }

  assembleInternal(extUpdate);
  assembleInternal(cnt);

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    
    for (int i = 0; i < cnt.subSize(iSub); ++i) {

      int c = cnt(iSub)[i];
      /*
      LevelSetStructure* fineLSS_sub = (fineLSS ? &(*fineLSS)(iSub) : NULL);
     
      if (!fineLSS_sub->isActive(0.0,i))
        continue;
      */
      if (c == 0)
	continue;

      //std::cout << "c = " << c << extUpdate(iSub)[i][0] << std::endl;

      for (int k = 0; k < dim; ++k)
      	coarseData(iSub)[i][k] += extUpdate(iSub)[i][k] / c;
    }
  }
  
}


template <class Scalar>
template<class Scalar2, int dim>
void MultiGridLevel<Scalar>::
ProjectResidual(DistSVec<Scalar2,dim>& r) const {

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {

    for(int i = 0; i < r(iSub).size(); ++i) {
  
      Topology myTopo = (*coarseNodeTopology)(iSub)[i];
      Vec3D topoNormal = (*coarseTopologyNormal)(iSub)[i];
      switch (nodeType[iSub][i]) {

        case BC_ADIABATIC_WALL_MOVING:
        case BC_ADIABATIC_WALL_FIXED:
        case BC_SLIP_WALL_MOVING:
        case BC_SLIP_WALL_FIXED:
        case BC_ISOTHERMAL_WALL_MOVING:
        case BC_ISOTHERMAL_WALL_FIXED:
        case BC_POROUS_WALL_MOVING:
        case BC_POROUS_WALL_FIXED:
        case BC_SYMMETRY: {

          double* pr = r(iSub)[i];
          double cmp = pr[1]*topoNormal[0]+pr[2]*topoNormal[1]+pr[3]*topoNormal[2];
          if (myTopo == TopoFace) {
 
            pr[1] -= cmp*topoNormal[0];
            pr[2] -= cmp*topoNormal[1];
            pr[3] -= cmp*topoNormal[2];
          } else if (myTopo == TopoLine) {
            pr[1] = cmp*topoNormal[0];
            pr[2] = cmp*topoNormal[1];
            pr[3] = cmp*topoNormal[2];
          }
        }
        default:
          break;
      }          
    }
  }  
}

//------------------------------------------------------------------------------

template<class Scalar>
bool MultiGridLevel<Scalar>::isLine(int iSub,int edgei,int edgej,int* lineid,int* loci,int* locj)
{

  int idm1,idp1;
  idm1 = (*lineMap)(iSub)[edgei][0];
  idp1 = (*lineMap)(iSub)[edgei][1];

  if (edgei == edgej) {

    if (idm1 >= 0 || idp1 >= 0) {
      *lineid = (*lineIDMap)(iSub)[edgei];
      *loci = (*lineLocIDMap)(iSub)[edgei];
      *locj = *loci;
      return true;
    }
  }

  if (idm1 == edgej || idp1 == edgej) {
    *lineid = (*lineIDMap)(iSub)[edgei];
    *loci = (*lineLocIDMap)(iSub)[edgei];
    *locj = (*lineLocIDMap)(iSub)[edgej];
    return true;
  } else {
    return false;
  }
  
}
//------------------------------------------------------------------------------
template<class Scalar>   
template <class Scalar2> 
void MultiGridLevel<Scalar>::assemble(DistVec<Scalar2>& V)
{

  operAdd<double> addOp;
  if (agglomType == AgglomerationLocal)
    ::assemble(domain, *(internal_select<Scalar2>(nodeIdPattern, nodeVolPattern)), sharedNodes, V, addOp);
  else
    assembleInternal(V);
}

template<class Scalar>   
template <class Scalar2,int dim> 
void MultiGridLevel<Scalar>::assemble(DistSVec<Scalar2,dim>& V)
{

  operAdd<double> addOp;
  if (agglomType == AgglomerationLocal) {
    if (dim == my_dim) {
      ::assemble(domain, *(internal_select<Scalar2>(nodeIdPattern, nodeVecPattern)) , sharedNodes, V, addOp);
    } else if (neq1 == dim) {
      ::assemble(domain, *(internal_select<Scalar2>(nodeIdPattern, nodeVecPatternEq1)), sharedNodes, V, addOp);
    } else if (neq2 == dim) {
      ::assemble(domain, *(internal_select<Scalar2>(nodeIdPattern, nodeVecPatternEq2)), sharedNodes, V, addOp);
    } else {
      std::cout << "Assemble called for dim = " << dim << std::endl;
      DebugTools::PrintBacktrace();
      exit(-1);
    }
  } else
    assembleInternal(V);    
}

template<class Scalar>   
template <class Scalar2,int dim> 
void MultiGridLevel<Scalar>::assemble(DistMat<Scalar2,dim>& A)
{
  //assert(dim == neq1);
  if (agglomType == AgglomerationLocal)
    ::assemble(domain, *matPattern, sharedNodes, A);
  else
    assembleInternal(A);
}

template<class Scalar>   
template <class Scalar2,int dim> 
void MultiGridLevel<Scalar>::assembleMax(DistSVec<Scalar2,dim>& V)
{

  operMax<Scalar2> maxOp;
//  std::cout << "Error: not used!" << std::endl;
  if (agglomType == AgglomerationLocal)
    ::assemble(domain, *nodeVecPattern, sharedNodes, V, maxOp);
}

template<class Scalar>   
template <class Scalar2> 
void MultiGridLevel<Scalar>::assembleMax(DistVec<Scalar2>& V) {

  operMax<Scalar2> maxOp;
  if (agglomType == AgglomerationLocal)
    ::assemble(domain, *(internal_select<Scalar2>(nodeIdPattern, nodeVecPattern))  ,  sharedNodes, V, maxOp);
  else
    assembleInternalMax(V);
}

template<class Scalar> 
int MultiGridLevel<Scalar>::getNumNeighborSubdomains(int mySub) {

  MultigridSubdomain* D = &mgSubdomains[mySub];

  return D->neighbors.size();
}

template<class Scalar> 
int* MultiGridLevel<Scalar>::getNeighboringSubdomains(int mySub) {

  MultigridSubdomain* D = &mgSubdomains[mySub];

  int* res = new int[D->neighbors.size()];
  
  int i = 0;
  for (std::map<int,NeighborDomain*>::iterator it = D->neighbors.begin();
       it != D->neighbors.end(); ++it, ++i) {

    res[i] = it->first;
  }
  
  return res;
}

template<class Scalar>
template <class Scalar2,int dim>
void MultiGridLevel<Scalar>::computeMatVecProd(DistMat<Scalar2,dim>& mat,
                                               DistSVec<Scalar2,dim>& p,
                                               DistSVec<Scalar2,dim>& prod) {

  prod = 0.0;

  int i,j,l;

  int num_edges = 0;

#pragma omp parallel for 
  for (int iSub = 0; iSub < p.numLocSub(); ++iSub) {
    
    bool* nodeFlag = nodeDistInfo->getMasterFlag(iSub);
    bool* edgeFlag = edges[iSub]->getMasterFlag();

    GenMat<Scalar2,dim>& A = mat(iSub);
    Scalar2 (*a)[dim*dim] = A.data(); 
    int numNodes = nodeDistInfo->subSize(iSub);
    for (i=0; i<numNodes; ++i) {

      if (!nodeFlag[i])
        continue;

      DenseMatrixOp<Scalar2,dim,dim*dim>::applyAndAddToVector(a, i, p(iSub).v, i, prod(iSub).v, i);
    }

    int (*edgePtr)[2] = edges[iSub]->getPtr();
 
    for (l = 0; l < edgeDistInfo->subSize(iSub); ++l) {

      if (!edgeFlag[l])
        continue;

      ++num_edges;

      i = edgePtr[l][0];
      j = edgePtr[l][1];

      DenseMatrixOp<Scalar2,dim,dim*dim>::applyAndAddToVector(a, numNodes + 2*l, p(iSub).v, j, prod(iSub).v, i);
      DenseMatrixOp<Scalar2,dim,dim*dim>::applyAndAddToVector(a, numNodes + 2*l + 1, p(iSub).v, i, prod(iSub).v, j);
    }

  } 

  operAdd<double> addOp;
  if (agglomType == AgglomerationLocal)
    ::assemble(domain, *nodeVecPattern, sharedNodes, prod, addOp,&locToGlobMap);
  else
    assembleInternal(prod);
}

template <class Scalar>
template <class Scalar2,int dim>
void MultiGridLevel<Scalar>::
computeGreenGaussGradient(DistSVec<Scalar2,dim>& V,
			  DistSVec<Scalar2,dim>& dX,
			  DistSVec<Scalar2,dim>& dY,
			  DistSVec<Scalar2,dim>& dZ,
                          DistLevelSetStructure* distLSS) {
  /*
  DistSVec<Scalar2,dim> Vtmp(V);
#pragma omp parallel for 
  for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {

    for (int i = 0; i < V(iSub).size(); ++i) {

      for (int k = 0; k < dim; ++k) {
	V(iSub)[i][k] = (*Xn)(iSub)[i][0] + 5.0*(*Xn)(iSub)[i][2];
      }
    }
  }
  */
  dX = dY = dZ = 0.0;
#pragma omp parallel for 
  for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {
 
    bool* edgeFlag = edges[iSub]->getMasterFlag();
    int (*edgePtr)[2] = edges[iSub]->getPtr();
    LevelSetStructure* lss = (distLSS ? &(*distLSS)(iSub) : NULL);
    for (int l = 0; l < edges[iSub]->size(); ++l) {

      if (!edgeFlag[l])
        continue;

      int i = edgePtr[l][0];
      int j = edgePtr[l][1];

      if (!lss || (lss->isActive(0.0,i) && lss->isActive(0.0,j))) {
        for (int k = 0; k < dim; ++k) {
          dX(iSub)[i][k] += 0.5*(V(iSub)[i][k]+V(iSub)[j][k])*(*edgeNormals)(iSub)[l][0];
          dY(iSub)[i][k] += 0.5*(V(iSub)[i][k]+V(iSub)[j][k])*(*edgeNormals)(iSub)[l][1];
          dZ(iSub)[i][k] += 0.5*(V(iSub)[i][k]+V(iSub)[j][k])*(*edgeNormals)(iSub)[l][2];
        
          dX(iSub)[j][k] -= 0.5*(V(iSub)[i][k]+V(iSub)[j][k])*(*edgeNormals)(iSub)[l][0];
          dY(iSub)[j][k] -= 0.5*(V(iSub)[i][k]+V(iSub)[j][k])*(*edgeNormals)(iSub)[l][1];
          dZ(iSub)[j][k] -= 0.5*(V(iSub)[i][k]+V(iSub)[j][k])*(*edgeNormals)(iSub)[l][2];
        }
      } else if (lss->isActive(0.0,i)) {
        for (int k = 0; k < dim; ++k) {
          if (k < 1 || k > 3) {
            dX(iSub)[i][k] += V(iSub)[i][k]*(*edgeNormals)(iSub)[l][0];
            dY(iSub)[i][k] += V(iSub)[i][k]*(*edgeNormals)(iSub)[l][1];
            dZ(iSub)[i][k] += V(iSub)[i][k]*(*edgeNormals)(iSub)[l][2];
	  }
        }
      } else if (lss->isActive(0.0,j)) {
        for (int k = 0; k < dim; ++k) {
          if (k < 1 || k > 3) {
            dX(iSub)[j][k] -= V(iSub)[j][k]*(*edgeNormals)(iSub)[l][0];
            dY(iSub)[j][k] -= V(iSub)[j][k]*(*edgeNormals)(iSub)[l][1];
            dZ(iSub)[j][k] -= V(iSub)[j][k]*(*edgeNormals)(iSub)[l][2];
	  }
        }
      }
    }

    for (int i = 0; i < agglomeratedFaces[iSub]->size(); ++i) {

      if (!(*agglomeratedFaces[iSub])[i].getMasterFlag())
	continue;
      int nd = (*agglomeratedFaces[iSub])[i].getNode();
      if (!lss || (lss->isActive(0.0,nd))) {
        for (int k = 0; k < dim; ++k) {
          dX(iSub)[nd][k] += V(iSub)[nd][k]*(*agglomeratedFaces[iSub])[i].getNormal()[0];
          dY(iSub)[nd][k] += V(iSub)[nd][k]*(*agglomeratedFaces[iSub])[i].getNormal()[1];
          dZ(iSub)[nd][k] += V(iSub)[nd][k]*(*agglomeratedFaces[iSub])[i].getNormal()[2];
        }
      }
    }
  }
  assemble(dX); 
  assemble(dY); 
  assemble(dZ);

#pragma omp parallel for 
  for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {

    for (int i = 0; i < V(iSub).size(); ++i) {

      double invvol = 1.0/getCtrlVol()(iSub)[i];
      for (int k = 0; k < dim; ++k) {
        dX(iSub)[i][k] *= invvol;
        dY(iSub)[i][k] *= invvol;
        dZ(iSub)[i][k] *= invvol;
      }
    }
  }   
  //dX = dY = dZ = 1.0;
  //dZ = V;

  //V = Vtmp;
}


template<class Scalar>
void MultiGridLevel<Scalar>::writePVTUFile(const char* filename) {

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if (rank == 0) {
    std::ofstream pvtufile((std::string(filename)+".pvtu").c_str());
    pvtufile << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    pvtufile << "\t<PUnstructuredGrid GhostLevel=\"0\">\n";
    pvtufile << "\t\t<PPoints>\n";
    pvtufile << "\t\t<PDataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\"/>\n";
    pvtufile << "\t\t</PPoints>\n";
    for (int i = 0; i < size; ++i)
      pvtufile << "<Piece Source=\"" << filename << i << ".vtu\"/>\n";
    pvtufile << "</PUnstructuredGrid> </VTKFile>\n";
  }

  char nfname[256];
  sprintf(nfname,"%s%d.vtu",filename,rank);
  std::ofstream vtufile(nfname);
  int numPoints = 0,numCells = 0;
  DistSVec<double,3>& X = getXn();
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    numPoints += X(iSub).size();
    numCells += faces[iSub]->size();
  }
  
  vtufile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  vtufile << "<UnstructuredGrid>";
  vtufile << "<Piece NumberOfPoints=\"" << numPoints << "\" NumberOfCells=\"" << numCells << "\">" << std::endl;
  vtufile << "<Points> <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
  int locOffsets[512];
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    for (int i = 0; i < X(iSub).size(); ++i) {

      vtufile << X(iSub)[i][0] << " " << X(iSub)[i][1] << " " << X(iSub)[i][2] << "\n";
    }
    if (iSub > 0)
      locOffsets[iSub] = locOffsets[iSub-1] + X(iSub-1).size();
    else
      locOffsets[iSub] = 0;
  }

  vtufile << "</DataArray>  </Points>\n";

  vtufile << "<Cells> <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n";
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    for (int i = 0; i < faces[iSub]->size(); ++i) {

      Face& F = (*faces[iSub])[i];
      vtufile << F[0]+locOffsets[iSub] << " " << F[1]+locOffsets[iSub]  << " " << F[2]+locOffsets[iSub] << "\n";
    }
  }
  vtufile << "</DataArray> <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">";
  int idx = 0;
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) { 
    for (int i = 0; i < faces[iSub]->size(); ++i,++idx) { 
      vtufile << idx*3+3 << " ";
    }
  }
  vtufile << "</DataArray> <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) { 
    for (int i = 0; i < faces[iSub]->size(); ++i) { 
      vtufile << 5 << " ";
    }
  }

  vtufile << "</DataArray>      </Cells>    </Piece>  </UnstructuredGrid></VTKFile>\n";

}

template<class Scalar>
void MultiGridLevel<Scalar>::writePVTUAgglomerationFile(const char* filename) {

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if (rank == 0) {
    std::ofstream pvtufile((std::string(filename)+".pvtu").c_str());
    pvtufile << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    pvtufile << "\t<PUnstructuredGrid GhostLevel=\"0\">\n";
    pvtufile << "\t\t<PPoints>\n";
    pvtufile << "\t\t<PDataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\"/>\n";
    pvtufile << "\t\t</PPoints>\n";
    for (int i = 0; i < size; ++i)
      pvtufile << "<Piece Source=\"" << filename << i << ".vtu\"/>\n";
    pvtufile << "</PUnstructuredGrid> </VTKFile>\n";
  }

  char nfname[256];
  sprintf(nfname,"%s%d.vtu",filename,rank);
  std::ofstream vtufile(nfname);
  int numPoints = 0,numCells = 0;
  DistSVec<double,3>& X = getXn();
  typedef std::map<std::pair<int,int>,int> EdgePointMap;
  EdgePointMap* edgePoints = new EdgePointMap[numLocSub];
  MultiGridLevel* finestLevel = getFinestLevel();
  std::vector<Vec3D> points;
  std::vector<std::pair<int,int> > outEdges;
  DistSVec<double,3>& Xf = finestLevel->getXn();
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    EdgeSet* E = finestLevel->edges[iSub];
    for (int i = 0; i < E->size(); ++i) {

      int I = mapFineToCoarse(iSub, E->getPtr()[i][0]);
      int J = mapFineToCoarse(iSub, E->getPtr()[i][1]);
      DistVec<Topology>* nt = getFinestTopology();
      if ((*nt)(iSub)[E->getPtr()[i][0]] == TopoInterior ||
          (*nt)(iSub)[E->getPtr()[i][1]] == TopoInterior)
        continue;
      if (I != J) {

        int tmp;
        I = E->getPtr()[i][0];
        J = E->getPtr()[i][1];
        if (I > J) { tmp = I; I = J; J = tmp; }
        edgePoints[iSub][std::pair<int,int>(I,J)] = points.size();
        Vec3D r = 0.0;
        for (int k = 0; k < 3; ++k)
          r[k] = Xf(iSub)[E->getPtr()[i][0]][k] + Xf(iSub)[E->getPtr()[i][1]][k];
        points.push_back(0.5*r);
      }
    }
  }

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    FaceSet& F = *finestLevel->faces[iSub];
    for (int i = 0; i < F.size(); ++i) {

      int I = mapFineToCoarse(iSub, F[i][0]);
      int J = mapFineToCoarse(iSub, F[i][1]);
      int K = mapFineToCoarse(iSub, F[i][2]);
      if (I == J && J == K) continue;

      Vec3D r = 0.0;
      for (int l = 0; l < 3; ++l)
        for (int k = 0; k < 3; ++k)
          r[k] += Xf(iSub)[F[i][l]][k]/3.0;
      
      int Ip = F[i][0];
      int Jp = F[i][1];
      int Kp = F[i][2];
     
      int tmp;
      if (Ip > Jp) { tmp = Ip; Ip = Jp; Jp = tmp;  tmp = I; I = J; J = tmp; }
      if (Jp > Kp) { tmp = Jp; Jp = Kp; Kp = tmp;   tmp = J; J = K; K = tmp;}
      if (Ip > Jp) { tmp = Ip; Ip = Jp; Jp = tmp;   tmp = I; I = J; J = tmp;}
      
      if (I != J) {

        outEdges.push_back(std::pair<int,int>(edgePoints[iSub][std::pair<int,int>(Ip,Jp)],
                                              points.size()));
      }
      if (J != K) {

        outEdges.push_back(std::pair<int,int>(edgePoints[iSub][std::pair<int,int>(Jp,Kp)],
                                              points.size()));
      }
      if (I != K) {

        outEdges.push_back(std::pair<int,int>(edgePoints[iSub][std::pair<int,int>(Ip,Kp)],
                                              points.size()));
      }
      points.push_back(r);
    }
  } 
/*
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    numPoints += X(iSub).size();
    numCells += faces[iSub]->size();
  }
  */
  numPoints = points.size();
  numCells = outEdges.size();
  vtufile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  vtufile << "<UnstructuredGrid>";
  vtufile << "<Piece NumberOfPoints=\"" << numPoints << "\" NumberOfCells=\"" << numCells << "\">" << std::endl;
  vtufile << "<Points> <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
  for (int i = 0; i < points.size(); ++i) {

    vtufile << points[i][0] << " " << points[i][1] << " " << points[i][2] << "\n";
  }
  
  vtufile << "</DataArray>  </Points>\n";

  vtufile << "<Cells> <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n";
  for (int i = 0; i < outEdges.size(); ++i) {

    vtufile << outEdges[i].first << " " << outEdges[i].second << "\n";
  }
  vtufile << "</DataArray> <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">";
  int idx = 0;
  for (int i = 0; i < outEdges.size(); ++i) { 
    vtufile << i*2+2 << " "; 
  }
  vtufile << "</DataArray> <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
  for (int i = 0; i < outEdges.size(); ++i) { 
    vtufile << 3 << " ";
  }

  vtufile << "</DataArray>      </Cells>    </Piece>  </UnstructuredGrid></VTKFile>\n";

  delete [] edgePoints;

}

template<class Scalar>
template<int dim>
void MultiGridLevel<Scalar>::writePVTUSolutionFile(const char* filename,
                                                   DistSVec<double,dim>& sol) {

  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if (rank == 0) {
    std::ofstream pvtufile((std::string(filename)+".pvtu").c_str());
    pvtufile << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    pvtufile << "\t<PUnstructuredGrid GhostLevel=\"0\">\n";
    pvtufile << "\t\t<PPoints>\n";
    pvtufile << "\t\t<PDataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\"/>\n";
    pvtufile << "\t\t</PPoints>\n";
    pvtufile << "<PCellData Vectors=\"Result\">\n" ;
    pvtufile << "  <PDataArray type=\"Float32\" Name=\"Result\" NumberOfComponents=\"" << 3/*dim*/ << 
            "\" format=\"appended\"/> </PCellData>\n";
    for (int i = 0; i < size; ++i)
      pvtufile << "<Piece Source=\"" << filename << i << ".vtu\"/>\n";
    pvtufile << "</PUnstructuredGrid> </VTKFile>\n";
  }

  char nfname[256];
  sprintf(nfname,"%s%d.vtu",filename,rank);
  std::ofstream vtufile(nfname);
  int numPoints = 0,numCells = 0;
  DistSVec<double,3>& X = getXn();
  typedef std::map<std::pair<int,int>,int> EdgePointMap;
  typedef std::map<int,int> VertexPointMap;
  EdgePointMap* edgePoints = new EdgePointMap[numLocSub];
  VertexPointMap* vertexPoints = new VertexPointMap[numLocSub];
  MultiGridLevel* finestLevel = getFinestLevel();
  std::vector<Vec3D> points;
  std::vector<int > outConnectivity;
  std::vector<int> outTypes;
  std::vector<int> outOffsets;
  std::vector<double> outData;
  DistSVec<double,3>& Xf = finestLevel->getXn();
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    EdgeSet* E = finestLevel->edges[iSub];
    for (int i = 0; i < E->size(); ++i) {

      int I = mapFineToCoarse(iSub, E->getPtr()[i][0]);
      int J = mapFineToCoarse(iSub, E->getPtr()[i][1]);
      //DistVec<Topology>* nt = getFinestTopology();
      //if ((*nt)(iSub)[E->getPtr()[i][0]] == TopoInterior ||
      //    (*nt)(iSub)[E->getPtr()[i][1]] == TopoInterior)
      //  continue;
      //if (I != J) {

        int tmp;
        I = E->getPtr()[i][0];
        J = E->getPtr()[i][1];
        if (I > J) { tmp = I; I = J; J = tmp; }
        edgePoints[iSub][std::pair<int,int>(I,J)] = points.size();
        Vec3D r = 0.0;
        for (int k = 0; k < 3; ++k)
          r[k] = Xf(iSub)[E->getPtr()[i][0]][k] + Xf(iSub)[E->getPtr()[i][1]][k];
        points.push_back(0.5*r);
      //}
    }

    for (int i = 0; i < finestLevel->nodeDistInfo->subSize(iSub); ++i) {

      Vec3D r(Xf(iSub)[i][0],Xf(iSub)[i][1],Xf(iSub)[i][2]);
      vertexPoints[iSub][i] = points.size();
      points.push_back(r);
    }
  }

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    FaceSet& F = *finestLevel->faces[iSub];
    for (int i = 0; i < F.size(); ++i) {

      int I = mapFineToCoarse(iSub, F[i][0]);
      int J = mapFineToCoarse(iSub, F[i][1]);
      int K = mapFineToCoarse(iSub, F[i][2]);
      int Ip = F[i][0];
      int Jp = F[i][1];
      int Kp = F[i][2];
     
      if (I == J && J == K) {
        
        outTypes.push_back(5);
        outConnectivity.push_back(vertexPoints[iSub][Ip]);
        outConnectivity.push_back(vertexPoints[iSub][Jp]);
        outConnectivity.push_back(vertexPoints[iSub][Kp]);
        outOffsets.push_back((outOffsets.size() > 0 ? 
                              outOffsets.back()+3: 3));
        for (int k = 0; k < 3/*dim*/; ++k)
          outData.push_back(sol(iSub)[I][k]); 
      }

      Vec3D r = 0.0;
      for (int l = 0; l < 3; ++l)
        for (int k = 0; k < 3; ++k)
          r[k] += Xf(iSub)[F[i][l]][k]/3.0;
      
      int tmp;

      std::pair<int,int> e1(Ip,Jp); if (Ip > Jp) { std::swap(e1.first,e1.second); }     
      std::pair<int,int> e2(Jp,Kp); if (Jp > Kp) { std::swap(e2.first,e2.second); }     
      std::pair<int,int> e3(Kp,Ip); if (Kp > Ip) { std::swap(e3.first,e3.second); }     

      outConnectivity.push_back(vertexPoints[iSub][Ip]); 
      outConnectivity.push_back(edgePoints[iSub][e1]); 
      outConnectivity.push_back(points.size()); 
      outConnectivity.push_back(edgePoints[iSub][e3]); 

      outConnectivity.push_back(vertexPoints[iSub][Jp]); 
      outConnectivity.push_back(edgePoints[iSub][e2]); 
      outConnectivity.push_back(points.size()); 
      outConnectivity.push_back(edgePoints[iSub][e1]); 

      outConnectivity.push_back(vertexPoints[iSub][Kp]); 
      outConnectivity.push_back(edgePoints[iSub][e3]); 
      outConnectivity.push_back(points.size()); 
      outConnectivity.push_back(edgePoints[iSub][e2]); 

      outOffsets.push_back((outOffsets.size() > 0 ?
                            outOffsets.back()+4: 4));      
      outOffsets.push_back((outOffsets.size() > 0 ?
                            outOffsets.back()+4: 4));      
      outOffsets.push_back((outOffsets.size() > 0 ?
                            outOffsets.back()+4: 4));      

      outTypes.push_back(9);
      outTypes.push_back(9);
      outTypes.push_back(9);
      for (int k = 0; k < 3;/*dim*/ ++k)
        outData.push_back(sol(iSub)[I][k]); 
      for (int k = 0; k < 3;/*dim*/ ++k)
        outData.push_back(sol(iSub)[J][k]); 
      for (int k = 0; k < 3/*dim*/; ++k)
        outData.push_back(sol(iSub)[K][k]); 

      points.push_back(r);
    }
  } 
/*
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    numPoints += X(iSub).size();
    numCells += faces[iSub]->size();
  }
  */
  numPoints = points.size();
  numCells = outTypes.size();
  vtufile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  vtufile << "<UnstructuredGrid>";
  vtufile << "<Piece NumberOfPoints=\"" << numPoints << "\" NumberOfCells=\"" << numCells << "\">" << std::endl;
  vtufile << "<Points> <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
  for (int i = 0; i < points.size(); ++i) {

    vtufile << points[i][0] << " " << points[i][1] << " " << points[i][2] << "\n";
  }
  
  vtufile << "</DataArray>  </Points>\n";

  vtufile << "<Cells> <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n";
  for (int i = 0; i < outConnectivity.size(); ++i) {

    vtufile << outConnectivity[i] << " " << "\n";
  }
  vtufile << "</DataArray> <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">";
  int idx = 0;
  for (int i = 0; i < outOffsets.size(); ++i) { 
    vtufile << outOffsets[i] << " "; 
  }
  vtufile << "</DataArray> <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
  for (int i = 0; i < outTypes.size(); ++i) { 
    vtufile << outTypes[i] << " ";
  }

  vtufile << "</DataArray>      </Cells>\n";

  vtufile << "<CellData Vectors=\"Result\">\n" ;
  vtufile << "  <DataArray type=\"Float32\" Name=\"Result\" NumberOfComponents=\"" << 3/*dim*/ << 
            "\" format=\"ascii\">\n";
  for (int i = 0; i < outData.size(); ++i) { 
    vtufile << (fabs(outData[i])<1.0e-10?0.0:outData[i]) << "\n";
  }
  vtufile << "</DataArray> </CellData>\n";

  vtufile <<"  </Piece>  </UnstructuredGrid></VTKFile>\n";

  delete [] edgePoints;
  delete [] vertexPoints;
}

//------------------------------------------------------------------------------

#define INSTANTIATION_HELPER(T,dim) \
  template void MultiGridLevel<T>::Restrict(const MultiGridLevel<T> &, const DistSVec<double,dim> &, DistSVec<double,dim> &,bool, bool) const; \
  template void MultiGridLevel<T>::RestrictFaceVector(const MultiGridLevel<T> &, const DistSVec<double,dim> &, DistSVec<double,dim> &) const; \
  template void MultiGridLevel<T>::Prolong(  MultiGridLevel<T> &, const DistSVec<double,dim> &, const DistSVec<double,dim> &, DistSVec<double,dim> &, \
DistSVec<double,dim>&,double, VarFcn*, DistLevelSetStructure*,DistLevelSetStructure*) const; \
template void MultiGridLevel<T>::ExtrapolateProlongation(MultiGridLevel<T>& fineGrid,  \
							 const DistSVec<double,dim>& coarseInitialData,\
							 DistSVec<double,dim>& coarseData,\
							 DistLevelSetStructure* coarseLSS, \
							 DistLevelSetStructure* fineLSS); \
  template void MultiGridLevel<T>::assemble(DistSVec<double,dim> &); \
  template void MultiGridLevel<T>::ProjectResidual(DistSVec<double,dim> &) const; \
  template void MultiGridLevel<T>::setupBcs(DistBcData<dim>&,DistBcData<dim>&,DistSVec<T,dim>&); \
  template void MultiGridLevel<T>::assemble(DistMat<double,dim> &); \
  template void MultiGridLevel<T>::assembleInternal(DistSVec<double,dim> &) const; \
  template void MultiGridLevel<T>::assembleInternal(DistMat<double,dim> &) const; \
  template void MultiGridLevel<T>::assembleMax(DistSVec<double,dim> &); \
  template void MultiGridLevel<T>::computeMatVecProd(DistMat<double,dim>& mat, \
                                                     DistSVec<double,dim>& p, \
                                                     DistSVec<double,dim>& prod); \
  template void MultiGridLevel<T>::RestrictOperator(const MultiGridLevel<T>& fineGrid, \
                                                    DistMat<double,dim>& fineOperator, \
                                                    DistMat<double,dim>& coarseOperator); \
  template void MultiGridLevel<T>::computeGreenGaussGradient(DistSVec<double,dim>& V, \
                               DistSVec<double,dim>& dX, \
                               DistSVec<double,dim>& dY, \
                               DistSVec<double,dim>& dZ, \
                               DistLevelSetStructure*); \
  template void MultiGridLevel<T>::writeXpostFile(const std::string&, \
                                                    DistSVec<T,dim>&, \
                                                    int); \
  template SparseMat<T,dim>* MultiGridLevel<T>::createMaskILU(int iSub,int fill,  \
                                           int renum, int *ndType); \
  template void MultiGridLevel<T>::writePVTUSolutionFile(const char*, \
                                                    DistSVec<double,dim>& \
                                                  ); 

#define INSTANTIATION_HELPER2(T) \
  template void MultiGridLevel<T>::assemble(DistVec<double> &); \
  template void MultiGridLevel<T>::assemble(DistVec<int> &); \
  template void MultiGridLevel<T>::assembleInternal(DistVec<double> &) const; \
  template void MultiGridLevel<T>::assembleMax(DistVec<double> &); \
  template void MultiGridLevel<T>::assembleMax(DistVec<int> &); \
  template void MultiGridLevel<T>::assembleInternalMax(DistVec<double> &) const; \
  template void MultiGridLevel<T>::Restrict(const MultiGridLevel<T> &, const DistVec<double> &, DistVec<double> &) const; 

template class MultiGridLevel<double>;
INSTANTIATION_HELPER2(double);
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
INSTANTIATION_HELPER(double,5);
INSTANTIATION_HELPER(double,6);
INSTANTIATION_HELPER(double,7);
/*template class MultiGridLevel<float>;
INSTANTIATION_HELPER2(float);
INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,5);
INSTANTIATION_HELPER(float,6);
*/

