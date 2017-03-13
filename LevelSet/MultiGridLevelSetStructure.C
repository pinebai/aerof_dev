#include "MultiGridLevelSetStructure.h"
#include "MultiGridLevel.h"

MultiGridLevelSetStructure::
MultiGridLevelSetStructure(DistMultiGridLevelSetStructure& lss,
			   SubDomain& sub,
			   Vec<int>& status,Vec<double>& distance,Vec<bool>& is_swept,Vec<bool>& is_active,Vec<bool>& is_occluded,Vec<bool>& edge_intersects,
			   Vec<Vec3D>& surfaceNormals,
			   LevelSetStructure* parent, int mySub,MultiGridLevel<double>* myLevel)
	: LevelSetStructure(status, distance, is_swept, is_active, is_occluded, edge_intersects, 
							  /*dummy from here*/  edge_intersects, distance, distance, surfaceNormals, status, 
							                                        distance, distance, surfaceNormals, status),
    distLSS(lss), parent(parent),
    subD(sub), mySub(mySub), myLevel(myLevel), edges(*myLevel->getEdges()[mySub]),
    status(status), distance(distance), is_swept(is_swept), is_active(is_active),
    surfaceNormals(surfaceNormals),
    is_occluded(is_occluded), edge_intersects(edge_intersects)
    {}

void MultiGridLevelSetStructure::
recompute() {

  status = -1;

  Vec<int>& stat = parent->getStatus();

  int N = stat.size();

  int rnk;
  MPI_Comm_rank(MPI_COMM_WORLD,&rnk);
  
  Vec<int>& nodeMapping =  myLevel->getNodeMapping()(mySub);
  for (int i = 0; i < N; ++i) {

    status[nodeMapping[i] ] = std::max<int>(status[nodeMapping[i] ] ,
					    stat[i]);
    
    //if (nodeMapping[i] == 1636 && rnk == 1)
    //  std::cout << i << " " << status[nodeMapping[i]] << " " << stat[i] << std::endl;
    if (stat[i] > 1)
      std::cout << "Error! status = " << stat[i] << std::endl;
  }

  EdgeSet& edges = *myLevel->getParent()->getEdges()[mySub];
  int (*ptr)[2] = edges.getPtr();
  for (int l = 0; l < edges.size(); ++l) {
 
    int i = ptr[l][0],j = ptr[l][1];
    if (nodeMapping[i] == nodeMapping[j] && 
        parent->edgeIntersectsStructure(0.0,l)  ) {
      status[nodeMapping[i]] = 1;
    }
  }

//  status = 0;
}

void MultiGridLevelSetStructure::
computeEdgeCrossing(SVec<double,3>& nodeNormals) {

  int N = edges.size();
  
  int (*ptr)[2] = edges.getPtr();

  for (int i = 0; i < N; ++i) {

    edge_intersects[i] = (status[ptr[i][0]] != status[ptr[i][1]]);
    
    double* v = nodeNormals[ptr[i][0]];
    if (v[0]*v[0] + v[1]*v[1]+v[2]*v[2] > 0.0) {

      surfaceNormals[i] = Vec3D(v[0],v[1],v[2]);
    } else {
      v = nodeNormals[ptr[i][1]];
      //assert(v[0]*v[0] + v[1]*v[1]+v[2]*v[2] > 0.0);
      
      surfaceNormals[i] = Vec3D(v[0],v[1],v[2]);
     
    }
    
  } 

  
}

LevelSetResult MultiGridLevelSetStructure::
getLevelSetDataAtEdgeCenter(double t, int l, bool i_less_j, double *Xr, double *Xg) {
  if (!edge_intersects[l]) {
    int (*ptr)[2] = edges.getPtr();
    int i=i_less_j ? ptr[l][0] : ptr[l][1],
        j=i_less_j ? ptr[l][1] : ptr[l][0];
    //fprintf(stderr,"%02d There is no intersection between node %d(status:%d,occluded=%d) and %d(status:%d,occluded=%d) along edge %d! Abort...\n",
    //               globIndex,locToGlobNodeMap[i]+1, status[i],is_occluded[i], locToGlobNodeMap[j]+1, status[j],is_occluded[j],l);
    fprintf(stderr,"There is no intersection between node %d(status:%d,occluded=%d) and %d(status:%d,occluded=%d) along edge %d! Abort...\n",
                   i, status[i],is_occluded[i], j, status[j],is_occluded[j],l);
    exit(-1);
  }
  
  LevelSetResult lsRes;

  lsRes.alpha = 0.5;
  lsRes.xi[0] = 1.0/3.0;
  lsRes.xi[1] = 1.0/3.0;
  lsRes.xi[2] = 1.0/3.0;
  
  lsRes.trNodes[0] = -1;
  lsRes.trNodes[1] = -1;
  lsRes.trNodes[2] = -1;

  lsRes.normVel = 0;

  lsRes.gradPhi = surfaceNormals[l];

  return lsRes;
}

int MultiGridLevelSetStructure::numOfFluids() {
  
  return distLSS.numOfFluids();
}


int
DistMultiGridLevelSetStructure::
recompute(double dtf, double dtfLeft, double dts, bool findStatus, bool retry) {

  
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    subLSS[iSub]->recompute(); 
  }
  
  //std::cout << status->size() << " " << status->v[1636] << std::endl;

  myLevel->assembleMax(*status);

  //std::cout << status->size() << " " << status->v[1636] << std::endl;

  DistSVec<double,3>* nodeNormals = 
    new DistSVec<double,3>(myLevel->getNodeDistInfo());

  *nodeNormals = 0.0;
  

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    int (*ptr)[2] = myLevel->getParent()->getEdges()[iSub]->getPtr();
    bool* masterFlag = myLevel->getParent()->getEdges()[iSub]->getMasterFlag();
    int nE = myLevel->getParent()->getEdges()[iSub]->size();

    for (int i = 0; i < nE; ++i) {

      if (!masterFlag[i])
	continue;

      if (!(*parent)(iSub).edgeIntersectsStructure(0.0, i))
	continue;
      
      LevelSetResult res = 
	(*parent)(iSub).getLevelSetDataAtEdgeCenter(0.0,i, !(*parent)(iSub).isActive(0.0, ptr[i][0]));

      int j = ( (*parent)(iSub).isActive(0.0, ptr[i][0]) ? ptr[i][1] : ptr[i][0] );

      j = myLevel->getNodeMapping()(iSub)[j];

      (*nodeNormals)(iSub)[j][0] += res.gradPhi[0];
      (*nodeNormals)(iSub)[j][1] += res.gradPhi[1];
      (*nodeNormals)(iSub)[j][2] += res.gradPhi[2];

      j = ( (*parent)(iSub).isActive(0.0, ptr[i][0]) ? ptr[i][0] : ptr[i][1] );

      j = myLevel->getNodeMapping()(iSub)[j];

      (*nodeNormals)(iSub)[j][0] += res.gradPhi[0];
      (*nodeNormals)(iSub)[j][1] += res.gradPhi[1];
      (*nodeNormals)(iSub)[j][2] += res.gradPhi[2];

    }    

  }

  myLevel->assemble(*nodeNormals);
  

  int nOfF = numOfFluids();


#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    int N = (*is_active)(iSub).size();
    for (int i = 0; i < N; ++i) {
      (*is_active)(iSub)[i] = !(*status)(iSub)[i];
      
      double* v = (*nodeNormals)(iSub)[i];
      double mag = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
      if (mag > 0.0) {
	mag = 1.0/sqrt(mag);
	v[0] *= mag;
	v[1] *= mag;
	v[2] *= mag;
      }
    }
	
	
  }

  *is_occluded = false;
/*
  DistVec<int> numActiveNeighbors(myLevel->getNodeDistInfo());
  numActiveNeighbors = 0;
  
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    int (*ptr)[2] = myLevel->getEdges()[iSub]->getPtr();
    bool* masterFlag = myLevel->getEdges()[iSub]->getMasterFlag();
    int nE = myLevel->getEdges()[iSub]->size();

    for (int i = 0; i < nE; ++i) {

      if (!masterFlag[i])
	continue;
      
      if (!(*status)(iSub)[ptr[i][0]])
	numActiveNeighbors(iSub)[ptr[i][1]]++;

      if (!(*status)(iSub)[ptr[i][1]])
	numActiveNeighbors(iSub)[ptr[i][0]]++;
    }
    
  }

  myLevel->assemble(numActiveNeighbors);

  
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    
    for (int i = 0; i < numActiveNeighbors(iSub).size(); ++i) {
      
      if (numActiveNeighbors(iSub)[i] <= 2) {
	(*status)(iSub)[i] = 1;
	(*is_active)(iSub)[i] = 0;
      }
	
    }
  }
    
  */
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    subLSS[iSub]->computeEdgeCrossing((*nodeNormals)(iSub));
  }

  delete nodeNormals;

  return 0;
}

DistMultiGridLevelSetStructure::
DistMultiGridLevelSetStructure(IoData &iod, Communicator *comm,
			       DistLevelSetStructure* parent,
			       MultiGridLevel<double>* level) : DistLevelSetStructure(),
								myLevel(level),
								parent(parent),
								com(comm) {

  

}

void DistMultiGridLevelSetStructure::
initialize(Domain * d, DistSVec<double,3> &X, DistSVec<double,3> &Xn, 
	   IoData &iod, DistVec<int> *point_based_id, 
	   DistVec<int>* oldStatus) {

  
  domain = d;
  numLocSub = d->getNumLocSub();
  subLSS = new MultiGridLevelSetStructure*[numLocSub];

  //closest = new DistVec<ClosestPoint>(myLevel->getNodeDistInfo()); //needed only for multi-phase cracking.
  status = new DistVec<int>(myLevel->getNodeDistInfo());
  distance = new DistVec<double>(myLevel->getNodeDistInfo());
  is_swept = new DistVec<bool>(myLevel->getNodeDistInfo());
  is_active = new DistVec<bool>(myLevel->getNodeDistInfo());
  is_occluded = new DistVec<bool>(myLevel->getNodeDistInfo());
  edge_intersects = new DistVec<bool>(myLevel->getEdgeDistInfo());

  surfaceNormals = new DistVec<Vec3D>(myLevel->getEdgeDistInfo());


#pragma omp parallel for
  for(int i = 0; i < numLocSub; ++i)
    subLSS[i] = new MultiGridLevelSetStructure(*this,*domain->getSubDomain()[i],
					       (*status)(i), (*distance)(i),
					       (*is_swept)(i),
					       (*is_active)(i),
					       (*is_occluded)(i),
					       (*edge_intersects)(i),
					       (*surfaceNormals)(i),
					       &(*parent)(i), i, 
					       myLevel);


  *distance=0.0;

}



LevelSetStructure &
DistMultiGridLevelSetStructure::
operator()(int subNum) const {
  return *subLSS[subNum];
}
