#include <SubDomain.h>

#include <BcDef.h>
#include <EdgeGrad.h>
#include <MacroCell.h>
#include <TimeData.h>
#include <Vector.h>
#include <Vector3D.h>
#include <SparseMatrix.h>
#include <Connectivity.h>
#include <Communicator.h>
#include <BinFileHandler.h>
#include <BCApplier.h>
#include <MatchNode.h>
#include <LinkF77.h>
#include <Extrapolation.h>

#include <cstdio>
#include <cmath>
#include <cstring>
#include <alloca.h>
#ifdef USE_EIGEN3
#include <Eigen/Sparse>
typedef Eigen::SparseMatrix<double> SpMat;
#endif


#ifdef OLD_STL
#include <defalloc.h>
#include <algo.h>
#else
#include <algorithm>
using std::stable_sort;
#endif

#include <RTree.h>

#include <TriangulatedInterface.h>

extern "C" {
  void F77NAME(psplotmask)(const int &, int *, int *, const int &, const int &);
};

//------------------------------------------------------------------------------

SubDomain::SubDomain(int locN, int clusN, int globN, int nClNd, char *clstN,
		     NodeSet *locNodes, FaceSet *locFaces, ElemSet *locElems,
		     int numLocNeighb, int *locNeighb, Connectivity *locSharedNodes,
		     int *locNodeMap, int *locFaceMap, int *locElemMap,
		     int nNdRng, int (*ndRng)[3]) :
  nodes(*locNodes), faces(*locFaces), elems(*locElems)
{

  locSubNum = locN;
  clusSubNum = clusN;
  globSubNum = globN;
  numClusNodes = nClNd;

  sampleMesh = false;

  sprintf(suffix, "%s", clstN);

  locToGlobNodeMap = locNodeMap;
  locToGlobFaceMap = locFaceMap;
  locToGlobElemMap = locElemMap;

  numNeighb = numLocNeighb;
  neighb = locNeighb;
  sharedNodes = locSharedNodes;
  numNodeRanges = nNdRng;
  nodeRanges = ndRng;

  sndChannel = 0;
  rcvChannel = 0;
  numSharedEdges = 0;
  sharedEdges = 0;
  nodeType = 0;
  nodeFaceType = 0;
  offWallNode = 0;
  mmsBCs = 0; //HB
  rotOwn = 0;
  nodesToMCNodes = 0;
  sharedInletNodes = 0;
  NodeToNode = 0;
  NodeToSubD = 0;
  NodeToElem = 0;
  ElemToElem = 0;
  nodeToNodeMaskJacobian = 0;
  nodeToNodeMaskILU      = 0;
  higherOrderMF = 0;
  higherOrderFSI = 0;

  geomTol = 1.0e-10;
  bigNum  = 1.0e16;

  int j;
  for(int i=0;i<3;i++)  {
    gradP[i] = new double[locNodes->size()];
    for (j = 0; j < locNodes->size(); j++)
      gradP[i][j] = 0.0;
  }
  dGradPSVec = new SVec<double,3>(locNodes->size()); 


// Included (MB*)
  numOffDiagEntries = 0;
  for(int i=0;i<3;i++)
   dGradP[i] = new double[locNodes->size()];

  myTree = new RTree<Elem>;  
  myTree->construct<&Elem::computeBoundingBox>(nodes, elems.getPointer(), elems.size());
}

//------------------------------------------------------------------------------

SubDomain::~SubDomain()
{

  if (locToGlobNodeMap) delete [] locToGlobNodeMap;
  if (locToGlobFaceMap) delete [] locToGlobFaceMap;
  if (locToGlobElemMap) delete [] locToGlobElemMap;
  if (neighb) delete [] neighb;
  if (sharedNodes) delete sharedNodes;
  if (sndChannel) delete [] sndChannel;
  if (rcvChannel) delete [] rcvChannel;
  if (numSharedEdges) delete [] numSharedEdges;

  if (sharedEdges) {
    for (int iSub = 0; iSub < numNeighb; ++iSub)
      if (sharedEdges[iSub])
	delete [] sharedEdges[iSub];

    delete [] sharedEdges;
  }

  if (nodeRanges) delete [] nodeRanges;
  if (nodeType) delete [] nodeType;
  if (nodeFaceType) delete [] nodeFaceType;
  if (offWallNode) delete [] offWallNode;
  if (nodesToMCNodes) delete [] nodesToMCNodes;
  if (sharedInletNodes) delete sharedInletNodes;
  if (NodeToNode) delete NodeToNode;
  if (NodeToSubD) delete NodeToSubD;
  if (NodeToElem) delete NodeToElem;
  if (ElemToElem) delete ElemToElem;
  //if (nodeToNodeMaskJacobian) delete nodeToNodeMaskJacobian;
  //if (nodeToNodeMaskILU) delete nodeToNodeMaskILU;
  delete &nodes;
  delete &faces;
  delete &elems;  
  for(int i=0;i<3;i++) 
    {
      delete[] dGradP[i];
      delete[] gradP[i];
    }
  if (rotOwn) delete[] rotOwn;
  if (mmsBCs) delete mmsBCs;
//  if (triaSurf) delete triaSurf;
//
  myTree->destruct();
  delete myTree;
  if(dGradPSVec) delete dGradPSVec;

}

//------------------------------------------------------------------------------

int SubDomain::numberEdges()
{

  int i;
  for (i=0; i<elems.size(); ++i)
    elems[i].numberEdges(edges);

  Vec<int> newNum(edges.size());

  edges.createPointers(newNum);

  for (i=0; i<elems.size(); ++i)
    elems[i].renumberEdges(newNum);

  for (i=0; i<faces.size(); ++i)
    faces[i].numberEdges(edges);

  return edges.size();
}

//------------------------------------------------------------------------------


void printMacroCells(Connectivity &cToN)

{

  // Get number of cells //

  int nCells = cToN.csize();

  // Print out the MacroCells Number and the Cells in it //

  for (int i = 0; i < nCells; ++i) {
    printf("MacroCell No %d \n", i+1);
    for (int j = 0; j < cToN.num(i); ++j) {
      printf("%d\n", cToN[i][j]);
    }
  }

}

//------------------------------------------------------------------------------

Connectivity *SubDomain::createElemBasedConnectivity()
{

  Connectivity eToN(&elems);
  Connectivity *nToE = eToN.reverse();
  Connectivity *nToN = nToE->transcon(&eToN);
  delete nToE;

  return nToN;

}

//------------------------------------------------------------------------------

Connectivity *SubDomain::createNodeToElementConnectivity()
{
  Connectivity eToN(&elems);
  Connectivity *nToE = eToN.reverse();
  return nToE;
}

//------------------------------------------------------------------------------

Connectivity *SubDomain::createElementToElementConnectivity()
{
  Connectivity eToN(&elems);
  Connectivity *nToE = eToN.reverse();
  Connectivity *eToE = eToN.transcon(nToE);
  delete nToE;

  return eToE;
}

//------------------------------------------------------------------------------

Connectivity *SubDomain::createElementToNodeConnectivity()
{
  return new Connectivity (&elems);
}
//------------------------------------------------------------------------------

Connectivity *SubDomain::createEdgeBasedConnectivity()
{

  int numNodes = nodes.size();

  int *numNeigh = reinterpret_cast<int *>(alloca(sizeof(int) * numNodes));

  int i;
  for (i=0; i<numNodes; ++i) numNeigh[i] = 1;

  int (*edgePtr)[2] = edges.getPtr();

  int l;
  for (l=0; l<edges.size(); ++l) {
    ++numNeigh[ edgePtr[l][0] ];
    ++numNeigh[ edgePtr[l][1] ];
  }

  int nnz = 0;
  for (i=0; i<numNodes; ++i) nnz += numNeigh[i];

  if (nnz != (2*edges.size() + numNodes)) {
    fprintf(stderr,"*** Error: wrong number of nonzero blocks\n");
    exit(1);
  }

  // construction of ia

  int *ia = new int[numNodes+1];

  ia[0] = 0;

  for (i=0; i<numNodes; ++i) ia[i+1] = ia[i] + numNeigh[i];

  // construction of ja

  int *ja = new int[nnz];

  for (i=0; i<numNodes; ++i) {
    ja[ia[i]] = i;
    numNeigh[i] = 1;
  }

  for (l=0; l<edges.size(); ++l) {

    int i1 = edgePtr[l][0];
    int i2 = edgePtr[l][1];

    ja[ ia[i1] + numNeigh[i1] ] = i2;
    ++numNeigh[i1];

    ja[ ia[i2] + numNeigh[i2] ] = i1;
    ++numNeigh[i2];

  }

  for (i=0; i<numNodes; ++i)
#ifdef OLD_STL
    sort(ja+ia[i], ja+ia[i+1]);
#else
    stable_sort(ja+ia[i], ja+ia[i+1]);
#endif

  Connectivity *nodeToNode = new Connectivity(numNodes, ia, ja);

  return nodeToNode;

}

//------------------------------------------------------------------------------

Connectivity *SubDomain::createNodeToConstantConnectivity()
{

  int numNodes = nodes.size();

  int *numNeigh = reinterpret_cast<int *>(alloca(sizeof(int) * numNodes));

  int i;
  for (i=0; i<numNodes; ++i) numNeigh[i] = 1;


  int l;
  int nnz = 0;
  for (i=0; i<numNodes; ++i) nnz += numNeigh[i];

  if (nnz != numNodes) {
    fprintf(stderr,"*** Error: wrong number of nonzero blocks\n");
    exit(1);
  }

  // construction of ia

  int *ia = new int[numNodes+1];

  ia[0] = 0;

  for (i=0; i<numNodes; ++i) ia[i+1] = ia[i] + numNeigh[i];

  // construction of ja

  int *ja = new int[nnz];

  for (i=0; i<numNodes; ++i) {
    ja[ia[i]] = 0;
    numNeigh[i] = 1;
  }

  for (i=0; i<numNodes; ++i)
#ifdef OLD_STL
    sort(ja+ia[i], ja+ia[i+1]);
#else
    stable_sort(ja+ia[i], ja+ia[i+1]);
#endif

  Connectivity *nodeToConstant = new Connectivity(numNodes, ia, ja);

  return nodeToConstant;

}


Connectivity *SubDomain::createConstantToNodeConnectivity()
{

  int numNodes = nodes.size();

  int *numNeigh = reinterpret_cast<int *>(alloca(sizeof(int) * 1));

  int i;
  numNeigh[0] = numNodes;


  int l;
  int nnz = numNodes;

  if (nnz != numNodes) {
    fprintf(stderr,"*** Error: wrong number of nonzero blocks\n");
    exit(1);
  }

  // construction of ia

  int *ia = new int[2];

  ia[0] = 0;
  ia[1] = numNeigh[0];

  // construction of ja

  int *ja = new int[nnz];

  for (l=0; l<numNodes; ++l) {
    ja[l] = l;
  }

  for (i=0; i<1; ++i)
#ifdef OLD_STL
    sort(ja+ia[i], ja+ia[i+1]);
#else
    stable_sort(ja+ia[i], ja+ia[i+1]);
#endif

  Connectivity *constantToNode = new Connectivity(1, ia, ja);

  return constantToNode;

}

//------------------------------------------------------------------------------

Connectivity *SubDomain::createConstantToConstantConnectivity()
{

  int *numNeigh = reinterpret_cast<int *>(alloca(sizeof(int) * 1));

  int i;
  numNeigh[0] = 1;


  int l;
  int nnz = 1;

  if (nnz != 1) {
    fprintf(stderr,"*** Error: wrong number of nonzero blocks\n");
    exit(1);
  }

  // construction of ia

  int *ia = new int[2];

  ia[0] = 0;
  ia[1] = numNeigh[0];

  // construction of ja

  int *ja = new int[nnz];

  for (l=0; l<1; ++l) {
    ja[l] = l;
  }

  for (i=0; i<1; ++i)
#ifdef OLD_STL
    sort(ja+ia[i], ja+ia[i+1]);
#else
    stable_sort(ja+ia[i], ja+ia[i+1]);
#endif

  Connectivity *constantToconstant = new Connectivity(1, ia, ja);

  return constantToconstant;

}

//------------------------------------------------------------------------------


//------------------------------------------------------------------------------

Connectivity *SubDomain::createElementBasedEdgeToNodeConnectivity()
{


  int numEdges = edges.size();
  std::set<int> *numNeigh = new std::set<int>[numEdges];

  int edgenum;
  int nodenum0, nodenum1, nodenum2, nudenum3;

  int i, l;
  for (i=0; i<elems.size(); ++i) {
    for (l=0; l<6; ++l) {
      edgenum = elems[i].edgeNum(l);
      numNeigh[edgenum].insert(elems[i].nodeNum(0));
      numNeigh[edgenum].insert(elems[i].nodeNum(1));
      numNeigh[edgenum].insert(elems[i].nodeNum(2));
      numNeigh[edgenum].insert(elems[i].nodeNum(3));
    }
  }

  // construction of ia
  int *ia = new int[numEdges+1];

  ia[0] = 0;

  for (l=0; l<numEdges; ++l) ia[l+1] = ia[l] + numNeigh[l].size();

  // construction of ja

  int nnz = ia[numEdges];
  int *ja = new int[nnz];

  int k=0;
  for (l=0; l<numEdges; ++l) {
    for (std::set<int>::iterator it=numNeigh[l].begin(); it!=numNeigh[l].end(); ++it) {
      ja[k] = *it;
      k++;
    }
  }


  for (l=0; l<numEdges; ++l)
#ifdef OLD_STL
    sort(ja+ia[l], ja+ia[l+1]);
#else
    stable_sort(ja+ia[l], ja+ia[l+1]);
#endif

  Connectivity *edgeToNode = new Connectivity(numEdges, ia, ja);
  delete [] numNeigh;

  return edgeToNode;

}

//------------------------------------------------------------------------------

Connectivity *SubDomain::createElementBasedNodeToEdgeConnectivity()
{

  int numNodes = nodes.size();

  int *numNeigh = reinterpret_cast<int *>(alloca(sizeof(int) * numNodes));

  int i;
  for (i=0; i<numNodes; ++i) numNeigh[i] = 0;

  int (*edgePtr)[2] = edges.getPtr();

  int l;
  for (l=0; l<edges.size(); ++l) {
    ++numNeigh[ edgePtr[l][0] ];
    ++numNeigh[ edgePtr[l][1] ];
  }

  int nnz = 0;
  for (i=0; i<numNodes; ++i) nnz += numNeigh[i];

  if (nnz != 2*edges.size()) {
    fprintf(stderr,"*** Error: wrong number of nonzero blocks\n");
    exit(1);
  }

  // construction of ia

  int *ia = new int[numNodes+1];

  ia[0] = 0;

  for (i=0; i<numNodes; ++i) ia[i+1] = ia[i] + numNeigh[i];

  // construction of ja

  int *ja = new int[nnz];

  for (i=0; i<numNodes; ++i) {
    numNeigh[i] = 0;
  }

  for (l=0; l<edges.size(); ++l) {

    int i1 = edgePtr[l][0];
    int i2 = edgePtr[l][1];

    ja[ ia[i1] + numNeigh[i1] ] = l;
    ++numNeigh[i1];

    ja[ ia[i2] + numNeigh[i2] ] = l;
    ++numNeigh[i2];

  }

  for (i=0; i<numNodes; ++i)
#ifdef OLD_STL
    sort(ja+ia[i], ja+ia[i+1]);
#else
    stable_sort(ja+ia[i], ja+ia[i+1]);
#endif

  Connectivity *nodeToEdge = new Connectivity(numNodes, ia, ja);

  return nodeToEdge;

}

//------------------------------------------------------------------------------

Connectivity *SubDomain::createNodeToFaceConnectivity()
{

  int numFaces = faces.size();
  int numNodes = nodes.size();

  int *numNeigh = reinterpret_cast<int *>(alloca(sizeof(int) * numNodes));

  int i;
  for (i=0; i<numNodes; ++i) numNeigh[i] = 0;

  int l;
  for (l=0; l<numFaces; ++l) {
    ++numNeigh[ faces[l][0] ];
    ++numNeigh[ faces[l][1] ];
    ++numNeigh[ faces[l][2] ];
  }

  int nnz = 0;
  for (i=0; i<numNodes; ++i) nnz += numNeigh[i];

  if (nnz != 3*faces.size()) {
    fprintf(stderr,"*** Error: wrong number of nonzero blocks\n");
    exit(1);
  }

  // construction of ia

  int *ia = new int[numNodes+1];

  ia[0] = 0;

  for (i=0; i<numNodes; ++i) ia[i+1] = ia[i] + numNeigh[i];

  // construction of ja
  int *ja = new int[nnz];

  for (i=0; i<numNodes; ++i) {
    numNeigh[i] = 0;
  }

  for (l=0; l<faces.size(); ++l) {

    int i1 = faces[l][0];
    int i2 = faces[l][1];
    int i3 = faces[l][2];

    ja[ ia[i1] + numNeigh[i1] ] = l;
    ++numNeigh[i1];

    ja[ ia[i2] + numNeigh[i2] ] = l;
    ++numNeigh[i2];

    ja[ ia[i3] + numNeigh[i3] ] = l;
    ++numNeigh[i3];
  }

  for (l=0; l<numNodes; ++l)
#ifdef OLD_STL
    sort(ja+ia[l], ja+ia[l+1]);
#else
    stable_sort(ja+ia[l], ja+ia[l+1]);
#endif

  Connectivity *nodeToFace = new Connectivity(numNodes, ia, ja);

  return nodeToFace;

}

//------------------------------------------------------------------------------

Connectivity *SubDomain::createFaceToNodeConnectivity()
{

  int numFaces = faces.size();

  int nnz = 3*numFaces;

  if (nnz != 3*faces.size()) {
    fprintf(stderr,"*** Error: wrong number of nonzero blocks\n");
    exit(1);
  }

  // construction of ia

  int *ia = new int[numFaces+1];

  ia[0] = 0;

  int l;
  for (l=0; l<numFaces; ++l) ia[l+1] = ia[l] + 3;

  // construction of ja

  int *ja = new int[nnz];

  for (l=0; l<faces.size(); ++l) {

    ja[ ia[l] ] = faces[l][0];
    ja[ ia[l] + 1 ] = faces[l][1];
    ja[ ia[l] + 2 ] = faces[l][2];

  }

  for (l=0; l<numFaces; ++l)
#ifdef OLD_STL
    sort(ja+ia[l], ja+ia[l+1]);
#else
    stable_sort(ja+ia[l], ja+ia[l+1]);
#endif

  Connectivity *faceToNode = new Connectivity(numFaces, ia, ja);

  return faceToNode;

}

//------------------------------------------------------------------------------

Connectivity *SubDomain::createNodeToSubDomainConnectivity()
{
  int numNodes = nodes.size();
  int *numOwner = reinterpret_cast<int *>(alloca(sizeof(int) * numNodes));
  for (int i=0; i<numNodes; i++)
    numOwner[i] = 1;

  for(int iSub=0; iSub<numNeighb; iSub++)
    for(int i=0; i<sharedNodes->num(iSub); i++)
      ++numOwner[(*sharedNodes)[iSub][i]];

  int nnz = 0;
  for(int i=0; i<numNodes; i++)
    nnz += numOwner[i];

  // construction of ia
  int *ia = new int[numNodes+1];
  ia[0] = 0;
  for(int i=0; i<numNodes; i++)
    ia[i+1] = ia[i] + numOwner[i];

  // construction of ja
  int *ja = new int[nnz];
  for (int i=0; i<numNodes; i++) {
    ja[ia[i]] = globSubNum;
    numOwner[i] = 1;
  }

  for(int iSub=0; iSub<numNeighb; iSub++)
    for(int i=0; i<sharedNodes->num(iSub); i++) {
      int me = (*sharedNodes)[iSub][i];
      ja[ia[me]+numOwner[me]] = neighb[iSub];
      ++numOwner[me];
    }

  // no need to sort.
  
  Connectivity *nodeToSubD = new Connectivity(numNodes, ia, ja);

  return nodeToSubD;

}

//------------------------------------------------------------------------------

void SubDomain::createSharedInletNodeConnectivity(int subd)
{        /* this function creates a connectivity for the shared inlet nodes but in the
         * indexation of the inletnodes in the subdomain we are considering/
         */

        int iSub, i, j;

//      count of the # of inlet nodes shared with each of the neighboring subdomains
        int *count = new int [numNeighb];
        for (iSub = 0; iSub<numNeighb; iSub++){
                count[iSub] = 0;
                for (i = 0; i < sharedNodes->num(iSub); i++){
                        if (nodeType[ (*sharedNodes)[iSub][i] ] == BC_INLET_MOVING ||
                            nodeType[ (*sharedNodes)[iSub][i] ] == BC_OUTLET_MOVING ||
                            nodeType[ (*sharedNodes)[iSub][i] ] == BC_INLET_FIXED ||
                            nodeType[ (*sharedNodes)[iSub][i] ] == BC_OUTLET_FIXED )
                                        count[iSub]++;
                }
        }

        /* As it is not equivalent to loop through the shared nodes and to loop
         * through the inlet nodes (there can be a redundancy of the nodes in the shared nodes
         * where as there can't be in the inlet nodes structure), we need to create two arrays:
         * each of size nodes.size(), one to store the number of subdomains a node is shared with
         * the second to store the number of the neighbouring subdomains it belongs to.
         */

        int *shareNum = new int[nodes.size()];
        int **shareSub = new int*[nodes.size()];
        for (i=0; i<nodes.size(); i++){
                shareNum[i] = 0;
                shareSub[i] = 0;
        }
//      count of the number of subDomains a node is shared with.
        for (iSub=0; iSub<numNeighb; iSub++)
                for (i = 0; i < sharedNodes->num(iSub); i++)
                        shareNum[ (*sharedNodes)[iSub][i] ]++;

//      verification and computation of the total number of nodes (with their redundancy)
        int n1=0;
        int n2=0;
        for (iSub=0; iSub<numNeighb; iSub++) n1 += count[iSub];
        for( i=0; i<nodes.size(); i++)
                if(nodeType[ i ] == BC_INLET_MOVING ||
                            nodeType[ i ] == BC_OUTLET_MOVING ||
                            nodeType[ i ] == BC_INLET_FIXED ||
                            nodeType[ i ] == BC_OUTLET_FIXED )
                        n2 += shareNum[i];


//  for each node in the subdomain, an array is created with the subdomains number
//  it is shared with.
        for (i=0; i<nodes.size(); i++)
                if(shareNum[i]>0)
                        shareSub[i] = new int[shareNum[i]];

        int *count2 = new int[nodes.size()];
        for (i=0; i<nodes.size(); i++) count2[i]=0;

        for (iSub=0; iSub<numNeighb; iSub++){
                for (i = 0; i < sharedNodes->num(iSub); i++){
                        shareSub[ (*sharedNodes)[iSub][i] ][count2[ (*sharedNodes)[iSub][i] ]] = iSub;
                        count2[ (*sharedNodes)[iSub][i] ]++;
                }
        }

//  creation of the connectivity itself with two arrays ia and ja
        int *ia = new int[numNeighb+1];
        ia[0]=0;
        for (iSub=0; iSub<numNeighb; iSub++) ia[iSub+1] = ia[iSub] + count[iSub];

        int *ja = new int[n1];

        int num;
        int *count3 = new int[numNeighb];
        for (iSub=0; iSub<numNeighb; iSub++) count3[iSub]=0;
        for (i=0; i<inletNodes.size(); i++){
                num = inletNodes[i].getNodeNum();
                if (shareNum[num]>0){
                        for (j=0; j<shareNum[num]; j++){
                                ja[ ia[ shareSub[num][j] ] + count3[ shareSub[num][j] ] ] = i;
                                count3[shareSub[num][j]]++;
                        }
                }
        }

        sharedInletNodes = new Connectivity(numNeighb, ia, ja);

        if (count) delete [] count;
        if (count2) delete [] count2;
        if (count3) delete [] count3;
        if (shareNum) delete [] shareNum;
        if (shareSub)
           for (i=0; i<nodes.size(); i++)
             if (shareSub[i]) delete [] shareSub[i];
               delete [] shareSub;

}

//------------------------------------------------------------------------------

Connectivity *SubDomain::createNodeToMacroCellNodeConnectivity(MacroCellSet *macroCells)
{

  int numNodes = nodes.size();
  int *numMacroNodes = reinterpret_cast<int *>(alloca(sizeof(int) * numNodes));
  int i, j;

  for (i=0; i<numNodes; ++i) numMacroNodes[i] = 0;

  for (i=0; i<nodes.size(); ++i)
    for (j=0; j<nodes.size(); ++j)
      if (macroCells->containing(i) == macroCells->containing(j)) ++numMacroNodes[i];

  // construction of ptr //
  int *ptr = new int[numNodes+1];
  ptr[0] = 0;

  for (i=0; i<numNodes; ++i) ptr[i+1] = ptr[i] + numMacroNodes[i];

  // construction of list //

  int nnz = 0;
  for (i=0; i<numNodes; ++i) nnz += numMacroNodes[i];

  int *list = new int[nnz];

  for (i=0; i<numNodes; ++i) numMacroNodes[i] = 0;

  for (i=0; i<nodes.size(); ++i)
    for (j=0; j<nodes.size(); ++j)
      if (macroCells->containing(i) == macroCells->containing(j)) {
        list[ptr[i]+numMacroNodes[i]] = j;
        ++numMacroNodes[i];
      }


  for (i=0; i<numNodes; ++i)
#ifdef OLD_STL
    sort(list+ptr[i], list+ptr[i+1]);
#else
    stable_sort(list+ptr[i], list+ptr[i+1]);
#endif

  Connectivity *nodeToMacroCellNode = new Connectivity(numNodes, ptr, list);

//  printMacroCells(*nodeToMacroCellNode); exit(1);

  return nodeToMacroCellNode;

}

//------------------------------------------------------------------------------

compStruct *SubDomain::createRenumbering(Connectivity *nodeToNode,
					 int typeRenum, int print)
{

  int numNodes = nodes.size();

  compStruct *nodeRenum = nodeToNode->renumByComponent(typeRenum);

  int *order = new int[numNodes];

  int i;
  for (i=0; i<numNodes; ++i) order[i] = -1;

  for (i=0; i<numNodes; ++i)
    if (nodeRenum->renum[i] >= 0) order[nodeRenum->renum[i]] = i;

  nodeRenum->order = order;

  int *ia = (*nodeToNode).ptr();
  int *ja = (*nodeToNode)[0];
  double aa;
  double (*a)[1] = reinterpret_cast<double (*)[1]>(&aa);

  SparseMat<double,1> A(numNodes, ia[numNodes], ia, ja, a, 0, 0);

  if (print == 1)
    F77NAME(psplotmask)(numNodes, ia, ja, 0, 10+globSubNum);

  A.permute(nodeRenum->renum);

  if (print == 1)
    F77NAME(psplotmask)(numNodes, ia, ja, 0, 50+globSubNum);

  return nodeRenum;

}

//------------------------------------------------------------------------------

void SubDomain::getElementStatistics(int &numNodes, int &numEdges,
				     int &numFaces, int &numElems)
{

  numNodes = nodes.size();
  numEdges = edges.size();
  numFaces = faces.size();
  numElems = elems.size();

}

// Get the local node number in the subdomain of the global node <id>
// Returns -1 if it does not exist.  Warning: This method is O(N)
int SubDomain::getLocalNodeNum(int globNodeNum) const {

  for (int i = 0; i < nodes.size(); ++i) {
    if (locToGlobNodeMap[i] == globNodeNum)
      return i;
  }
  
  return -1;
}

//------------------------------------------------------------------------------

int SubDomain::computeControlVolumes(int numInvElem, double lscale,
				     SVec<double,3> &X, Vec<double> &ctrlVol)
{

  int i, ierr = 0;

  ctrlVol = 0.0;

  for (i=0; i<elems.size(); ++i) {
    double volume = elems[i].computeControlVolumes(X, ctrlVol);

    if (volume <= 0.0) {
      fprintf(stderr,"Element %i has a negative volume of %e\n",locToGlobElemMap[i]+1,volume);
      ++ierr;
      ++numInvElem;
      if (numInvElem)
       elems[i].printInvalidElement(numInvElem, lscale, i, locToGlobNodeMap,
                                    locToGlobElemMap, nodes, X);
    }
  }
/*
#ifdef YDEBUG
  if(ierr > 0) {
    const char* output = "elementvolumecheck";
    ofstream out(output, ios::out);
    if(!out) { cerr << "Error: cannot open file" << output << endl;  exit(-1); }
    out << ierr << endl;
    out.close();
    exit(-1);
  }
#endif
*/
  return ierr;

}

//------------------------------------------------------------------------------

// Included (MB)
int SubDomain::computeDerivativeOfControlVolumes(int numInvElem, double lscale,
                                                 SVec<double,3> &X, SVec<double,3> &dX, Vec<double> &dCtrlVol)
{

  dCtrlVol = 0.0;

  for (int i=0; i<elems.size(); ++i) {
    double dVolume = elems[i].computeDerivativeOfControlVolumes(X, dX, dCtrlVol);
  }

  return 0;

}

//------------------------------------------------------------------------------

// Included (YC)
int SubDomain::computeDerivativeOfControlVolumes(RectangularSparseMat<double,3,1> &dCtrlVoldX,
                                                 SVec<double,3> &dX, Vec<double> &dCtrlVol)
{

  dCtrlVol = 0.0;
  
  dCtrlVoldX.apply(dX, dCtrlVol);

  return 0;

}

//------------------------------------------------------------------------------

// Included (YC)
int SubDomain::computeTransposeDerivativeOfControlVolumes(RectangularSparseMat<double,3,1> &dCtrlVoldX,
                                                          Vec<double> &dCtrlVol, SVec<double,3> &dX)
{
 
  SVec<double,3> dummy(dX);
  dummy = 0.0; 
  dCtrlVoldX.applyTranspose(dCtrlVol, dummy);
  dX += dummy;
  return 0;

}

//------------------------------------------------------------------------------

// Included (YC)
int SubDomain::computeDerivativeOperatorsOfControlVolumes(SVec<double,3> &X, RectangularSparseMat<double,3,1> &dCtrlVoldX)
{

  for (int i=0; i<elems.size(); ++i) {
    elems[i].computeDerivativeOperatorsOfControlVolumes(X, dCtrlVoldX);
  }

  return 0;

}

//------------------------------------------------------------------------------

void SubDomain::computeFaceNormals(SVec<double,3> &X, Vec<Vec3D> &faceNorm)
{

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeNormal(X, faceNorm);

}

//------------------------------------------------------------------------------

void SubDomain::computeFaceEdgeNormals(SVec<double,3>& X, SVec<double,6>& normals)
{

  normals = 0.0;
  for (int i=0; i<faces.size(); ++i)
    faces[i].computeEdgeNormals(X, locToGlobNodeMap, normals);

}

//------------------------------------------------------------------------------

void SubDomain::computeEdgeDihedralAngle(double threshold, SVec<double,6>& normals,
					 Vec<double>& tag)
{

  tag = 0.0;

  Vec<bool> bedges(edges.size());
  bedges = false;
  for (int i=0; i<faces.size(); ++i)
    faces[i].tagEdgesOnBoundaries(bedges);

  int (*ptr)[2] = edges.getPtr();
  for (int l=0; l<edges.size(); ++l) {
    if (bedges[l]) {
      Vec3D ni(normals[l][0], normals[l][1], normals[l][2]);
      Vec3D nj(normals[l][3], normals[l][4], normals[l][5]);
      /*
      if ( (locToGlobNodeMap[ptr[l][0]]+1 == 17 && locToGlobNodeMap[ptr[l][1]]+1 == 157) ||
	   (locToGlobNodeMap[ptr[l][1]]+1 == 17 && locToGlobNodeMap[ptr[l][0]]+1 == 157) )
	printf("%d (%d -> %d, %d) (%d -> %d, %d) %e %e %e %e %e %e %e\n",
	       l+1, ptr[l][0], locToGlobNodeMap[ptr[l][0]]+1, nodeType[ ptr[l][0] ],
	       ptr[l][1], locToGlobNodeMap[ptr[l][1]]+1, nodeType[ ptr[l][1] ],
	       ni[0], ni[1], ni[2], nj[0], nj[1], nj[2], acos(ni*nj) *180.0/acos(-1.0));
      */
      if (acos(ni*nj) >= threshold) {
	if (nodeType[ ptr[l][0] ] < BC_INTERNAL)
	  tag[ ptr[l][0] ] += 1.0;
	if (nodeType[ ptr[l][1] ] < BC_INTERNAL)
	  tag[ ptr[l][1] ] += 1.0;
      }
    }
  }

}

//------------------------------------------------------------------------------

void SubDomain::propagateInfoAlongEdges(Vec<double>& tag)
{

  int (*edgePtr)[2] = edges.getPtr();
  Vec<double> newtag(tag.size());
  newtag = 0.0;
  for (int l=0; l<edges.size(); ++l) {
    int i = edgePtr[l][0];
    int j = edgePtr[l][1];
    if (tag[i] > 0.0)
      newtag[j] += 1.0;
    if (tag[j] > 0.0)
      newtag[i] += 1.0;
  }

  tag += newtag;

}

//------------------------------------------------------------------------------

// Included (MB)
void SubDomain::computeDerivativeOfNormals(SVec<double,3> &X, SVec<double,3> &dX,
                                           Vec<Vec3D> &edgeNorm, Vec<Vec3D> &dEdgeNorm, Vec<double> &edgeNormVel, Vec<double> &dEdgeNormVel,
                                           Vec<Vec3D> &faceNorm, Vec<Vec3D> &dFaceNorm, Vec<double> &faceNormVel, Vec<double> &dFaceNormVel)
{
  dEdgeNorm = 0.0;
  dEdgeNormVel = 0.0;

  int i;
  for (i=0; i<elems.size(); ++i)
    elems[i].computeDerivativeOfEdgeNormals(X, dX, edgeNorm, dEdgeNorm, edgeNormVel, dEdgeNormVel);


  for (i=0; i<faces.size(); ++i)
    faces[i].computeDerivativeOfNormal(X, dX, faceNorm[i], dFaceNorm[i], faceNormVel[i], dFaceNormVel[i]);
}

//------------------------------------------------------------------------------

// Included (YC)
void SubDomain::computeDerivativeOfNormals(RectangularSparseMat<double,3,3> &dEdgeNormdX,
                                           RectangularSparseMat<double,3,3> &dFaceNormdX,
                                           SVec<double,3> &dX,
                                           Vec<Vec3D> &dEdgeNorm, 
                                           Vec<Vec3D> &dFaceNorm)
{

  dEdgeNorm = 0.0;

  dEdgeNormdX.apply(dX, dEdgeNorm);
  dFaceNormdX.apply(dX, dFaceNorm);

}

//------------------------------------------------------------------------------

// Included (YC)
void SubDomain::computeTransposeDerivativeOfNormals(RectangularSparseMat<double,3,3> &dEdgeNormdX,
                                                    RectangularSparseMat<double,3,3> &dFaceNormdX,
                                                    Vec<Vec3D> &dEdgeNorm, 
                                                    Vec<Vec3D> &dFaceNorm,
                                                    SVec<double,3> &dX)
{

  SVec<double,3> dummy(dX);
  dummy = 0.0;
  dEdgeNormdX.applyTranspose(dEdgeNorm, dummy);
  dX += dummy;
  dummy = 0.0;
  dFaceNormdX.applyTranspose(dFaceNorm, dummy);
  dX += dummy;
  return;

}

//------------------------------------------------------------------------------

// Included (YC)
void SubDomain::computeDerivativeOperatorsOfNormals(SVec<double,3> &X, 
                                                    RectangularSparseMat<double,3,3> &dEdgeNormdX, 
                                                    RectangularSparseMat<double,3,3> &dFaceNormdX)
{

  int i;
  for (i=0; i<elems.size(); ++i)
    elems[i].computeDerivativeOperatorsOfEdgeNormals(X, dEdgeNormdX);

  for (i=0; i<faces.size(); ++i)
    faces[i].computeDerivativeOperatorsOfNormal(i, X, dFaceNormdX);

}

//------------------------------------------------------------------------------

void SubDomain::computeNormalsConfig(SVec<double,3> &Xconfig, SVec<double,3> &Xdot,
                                     Vec<Vec3D> &edgeNorm, Vec<double> &edgeNormVel,
                                     Vec<Vec3D> &faceNorm, Vec<double> &faceNormVel)
{
  int i;
  for (i=0; i<elems.size(); ++i)
    elems[i].computeEdgeNormalsConfig(Xconfig, Xdot, edgeNorm, edgeNormVel);

  for (i=0; i<faces.size(); ++i)
    faces[i].computeNormalConfig(Xconfig, Xdot, faceNorm, faceNormVel);

}

//------------------------------------------------------------------------------

void SubDomain::computeNormalsGCL1(SVec<double,3> &Xn, SVec<double,3> &Xnp1,
				   SVec<double,3> &Xdot,
				   Vec<Vec3D> &edgeNorm, Vec<double> &edgeNormVel,
				   Vec<Vec3D> &faceNorm, Vec<double> &faceNormVel)
{

  edgeNorm = 0.0;
  edgeNormVel = 0.0;

  int i;
  for (i=0; i<elems.size(); ++i)
    elems[i].computeEdgeNormalsGCL1(Xn, Xnp1, Xdot, edgeNorm, edgeNormVel);

  for (i=0; i<faces.size(); ++i)
    faces[i].computeNormalGCL1(Xn, Xnp1, Xdot, faceNorm, faceNormVel);

}

//------------------------------------------------------------------------------

void SubDomain::computeNormalsEZGCL1(double oodt, SVec<double,3>& Xn, SVec<double,3>& Xnp1,
				     Vec<Vec3D>& edgeNorm, Vec<double>& edgeNormVel,
				     Vec<Vec3D>& faceNorm, Vec<double>& faceNormVel)
{

  edgeNorm = 0.0;
  edgeNormVel = 0.0;

  int i;
  for (i=0; i<elems.size(); ++i)
    elems[i].computeEdgeNormalsEZGCL1(oodt, Xn, Xnp1, edgeNorm, edgeNormVel);

  for (i=0; i<faces.size(); ++i)
    faces[i].computeNormalEZGCL1(oodt, Xn, Xnp1, faceNorm, faceNormVel);

}

//------------------------------------------------------------------------------

void SubDomain::computeSmoothedSensor(SVec<double,3>& X, Vec<double>& sigma,
                                      SVec<double,3>& sensor)
{

  sensor = 0.0;

  bool* edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {
    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    double dx[3];
    dx[0] = X[j][0] - X[i][0];
    dx[1] = X[j][1] - X[i][1];
    dx[2] = X[j][2] - X[i][2];
    double oolength = 1.0 / sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

    sensor[i][0] += oolength * sigma[j];
    sensor[i][1] += oolength;
    sensor[j][0] += oolength * sigma[i];
    sensor[j][1] += oolength;
  }

}

//------------------------------------------------------------------------------
/*
@ARTICLE{haselbacher-blazek-00,
  author = "Haselbacher, A. and Blazek, J.",
  title = "Accurate and efficient discretization of {N}avier-{S}tokes
           equations on mixed grids",
  journal = aiaaj,
  year = 2000,
  volume = 38,
  number = 11,
  pages = "2094--2102",
  month = nov,
}
*/
void SubDomain::computeWeightsLeastSquaresEdgePart(SVec<double,3> &X, SVec<double,6> &R)
{

  R = 0.0;

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    double dx[3];
    dx[0] = X[j][0] - X[i][0];
    dx[1] = X[j][1] - X[i][1];
    dx[2] = X[j][2] - X[i][2];

    double dxdx = dx[0] * dx[0];
    double dydy = dx[1] * dx[1];
    double dzdz = dx[2] * dx[2];
    double dxdy = dx[0] * dx[1];
    double dxdz = dx[0] * dx[2];
    double dydz = dx[1] * dx[2];

    R[i][0] += dxdx;
    R[j][0] += dxdx;

    R[i][1] += dxdy;
    R[j][1] += dxdy;

    R[i][2] += dxdz;
    R[j][2] += dxdz;

    R[i][3] += dydy;
    R[j][3] += dydy;

    R[i][4] += dydz;
    R[j][4] += dydz;

    R[i][5] += dzdz;
    R[j][5] += dzdz;

  }

}

//------------------------------------------------------------------------------

// Included (YC)
void SubDomain::computeDerivativeTransposeOfWeightsLeastSquaresEdgePart(SVec<double,3> &X, SVec<double,6> &dR, SVec<double,6> &R, SVec<double,3> &dX)
{

  R = 0.0;
  dX = 0.0;

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    double dx[3];
    dx[0] = X[j][0] - X[i][0];
    dx[1] = X[j][1] - X[i][1];
    dx[2] = X[j][2] - X[i][2];

    double dr[6];
    dr[0] = dR[j][0] + dR[i][0];
    dr[1] = dR[j][1] + dR[i][1];
    dr[2] = dR[j][2] + dR[i][2];
    dr[3] = dR[j][3] + dR[i][3];
    dr[4] = dR[j][4] + dR[i][4];
    dr[5] = dR[j][5] + dR[i][5];

    double dxdx = dx[0] * dx[0];
    double dydy = dx[1] * dx[1];
    double dzdz = dx[2] * dx[2];
    double dxdy = dx[0] * dx[1];
    double dxdz = dx[0] * dx[2];
    double dydz = dx[1] * dx[2];

    R[i][0] += dxdx;
    R[j][0] += dxdx;

    R[i][1] += dxdy;
    R[j][1] += dxdy;

    R[i][2] += dxdz;
    R[j][2] += dxdz;

    R[i][3] += dydy;
    R[j][3] += dydy;

    R[i][4] += dydz;
    R[j][4] += dydz;

    R[i][5] += dzdz;
    R[j][5] += dzdz;

    dX[j][0] += 2*dx[0]*dr[0] + dx[1]*dr[1] + dx[2]*dr[2];  
    dX[i][0] += -(2*dx[0]*dr[0] + dx[1]*dr[1] + dx[2]*dr[2]); 

    dX[j][1] += dx[0]*dr[1] + 2*dx[1]*dr[3] + dx[2]*dr[4];  
    dX[i][1] += -(dx[0]*dr[1] + 2*dx[1]*dr[3] + dx[2]*dr[4]);

    dX[j][2] += dx[0]*dr[2] + dx[1]*dr[4] + 2*dx[2]*dr[5];  
    dX[i][2] += -(dx[0]*dr[2] + dx[1]*dr[4] + 2*dx[2]*dr[5]);

  }

}

//------------------------------------------------------------------------------
// Include (YC)
void SubDomain::compute_dRdX(SVec<double,3> &X, RectangularSparseMat<double,3,6> &dRdX)
{

  // size of dRdX must be (6*numNodes, 3*numNodes)

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    double dx[3];
    dx[0] = X[j][0] - X[i][0];
    dx[1] = X[j][1] - X[i][1];
    dx[2] = X[j][2] - X[i][2];

    double C[2*6][2*3] = {0};
    C[6][0] = C[0][0] = -2*dx[0];      C[6][3] = C[0][3] = 2*dx[0];
    C[7][0] = C[1][0] = -dx[1];        C[7][3] = C[1][3] = dx[1];
    C[7][1] = C[1][1] = -dx[0];        C[7][4] = C[1][4] = dx[0];
    C[8][0] = C[2][0] = -dx[2];        C[8][3] = C[2][3] = dx[2];
    C[8][2] = C[2][2] = -dx[0];        C[8][5] = C[2][5] = dx[0];
    C[9][1] = C[3][1] = -2*dx[1];      C[9][4] = C[3][4] = 2*dx[1];
    C[10][1] = C[4][1] = -dx[2];       C[10][4] = C[4][4] = dx[2];
    C[10][2] = C[4][2] = -dx[1];       C[10][5] = C[4][5] = dx[1];
    C[11][2] = C[5][2] = -2*dx[2];     C[11][5] = C[5][5] = 2*dx[2];

    int ndList[2] = {i, j};
    dRdX.addContrib(2, ndList, C[0]);
  }

}

//------------------------------------------------------------------------------
// Included (YC)
void SubDomain::computeDerivativeOfWeightsLeastSquaresEdgePart(RectangularSparseMat<double,3,6> &dRdX, SVec<double,3> &dX, SVec<double,6> &dR)
{

  dR = 0.0;
  dRdX.apply(dX, dR);

}

//------------------------------------------------------------------------------
// Included (YC)
void SubDomain::computeTransposeDerivativeOfWeightsLeastSquaresEdgePart(RectangularSparseMat<double,3,6> &dRdX, SVec<double,6> &dR, SVec<double,3> &dX)
{

  SVec<double,3> dummy(dX);
  dummy = 0.0;
  dRdX.applyTranspose(dR, dummy);
  dX += dummy;

}

//------------------------------------------------------------------------------
// Included (MB)
// YC: if you intend to modified this routine,
// you should also modify SubDomain::computeDerivativeTransposeOfWeightsLeastSquaresEdgePart accordingly
// and                    SubDomain::compute_dRdX
void SubDomain::computeDerivativeOfWeightsLeastSquaresEdgePart(SVec<double,3> &X, SVec<double,3> &dX, SVec<double,6> &R, SVec<double,6> &dR)
{
  R = 0.0;
  dR = 0.0;

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    double dx[3];
    dx[0] = X[j][0] - X[i][0];
    dx[1] = X[j][1] - X[i][1];
    dx[2] = X[j][2] - X[i][2];

    double ddx[3];
    ddx[0] = dX[j][0] - dX[i][0];
    ddx[1] = dX[j][1] - dX[i][1];
    ddx[2] = dX[j][2] - dX[i][2];

    double dxdx = dx[0] * dx[0];
    double dydy = dx[1] * dx[1];
    double dzdz = dx[2] * dx[2];
    double dxdy = dx[0] * dx[1];
    double dxdz = dx[0] * dx[2];
    double dydz = dx[1] * dx[2];

    double ddxdx = ddx[0] * dx[0] + dx[0] * ddx[0];
    double ddydy = ddx[1] * dx[1] + dx[1] * ddx[1];
    double ddzdz = ddx[2] * dx[2] + dx[2] * ddx[2];
    double ddxdy = ddx[0] * dx[1] + dx[0] * ddx[1];
    double ddxdz = ddx[0] * dx[2] + dx[0] * ddx[2];
    double ddydz = ddx[1] * dx[2] + dx[1] * ddx[2];

    R[i][0] += dxdx;
    R[j][0] += dxdx;

    R[i][1] += dxdy;
    R[j][1] += dxdy;

    R[i][2] += dxdz;
    R[j][2] += dxdz;

    R[i][3] += dydy;
    R[j][3] += dydy;

    R[i][4] += dydz;
    R[j][4] += dydz;

    R[i][5] += dzdz;
    R[j][5] += dzdz;

    dR[i][0] += ddxdx;
    dR[j][0] += ddxdx;

    dR[i][1] += ddxdy;
    dR[j][1] += ddxdy;

    dR[i][2] += ddxdz;
    dR[j][2] += ddxdz;

    dR[i][3] += ddydy;
    dR[j][3] += ddydy;

    dR[i][4] += ddydz;
    dR[j][4] += ddydz;

    dR[i][5] += ddzdz;
    dR[j][5] += ddzdz;
  }
}

//------------------------------------------------------------------------------
// least square gradient involving only nodes of same fluid (multiphase flow)
//d2d$
void SubDomain::computeWeightsLeastSquaresEdgePart(SVec<double,3> &X, const Vec<int> &fluidId,
                                                   SVec<int,1> &count, SVec<double,6> &R, 
						   LevelSetStructure *LSS, bool includeSweptNodes)
{

  R = 0.0;
  count = 0;

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

	for(int l=0; l<edges.size(); ++l) 
	{
    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

		bool validEdge = true;

		if(fluidId[i] != fluidId[j]) validEdge = false;

		if(LSS)
		{			
			if(LSS->edgeWithSI(l) || LSS->edgeIntersectsStructure(0.0, l)) validEdge = false;
			if(!LSS->isActive(0.0, i) || !LSS->isActive(0.0, j)) validEdge = false;			
			if(!includeSweptNodes && (LSS->isSwept(0.0, i) || LSS->isSwept(0.0, j))) validEdge = false;
		}

		if(!validEdge) continue;

    count[i][0]++;
    count[j][0]++;

    double dx[3];
    dx[0] = X[j][0] - X[i][0];
    dx[1] = X[j][1] - X[i][1];
    dx[2] = X[j][2] - X[i][2];

    double dxdx = dx[0] * dx[0];
    double dydy = dx[1] * dx[1];
    double dzdz = dx[2] * dx[2];
    double dxdy = dx[0] * dx[1];
    double dxdz = dx[0] * dx[2];
    double dydz = dx[1] * dx[2];

    R[i][0] += dxdx;
		R[i][1] += dxdy;
		R[i][2] += dxdz;
		R[i][3] += dydy;
		R[i][4] += dydz;
		R[i][5] += dzdz;
		
    R[j][0] += dxdx;
    R[j][1] += dxdy;
    R[j][2] += dxdz;
    R[j][3] += dydy;
    R[j][4] += dydz;
    R[j][5] += dzdz;
  }
}

//------------------------------------------------------------------------------
// least square gradient involving only nodes of same fluid (FSI)
// with option to take into account of Riemann solution at interface
void SubDomain::computeWeightsLeastSquaresEdgePart(SVec<double,3> &X, const Vec<int> &fluidId,
                                                   SVec<int,1> &count, SVec<double,6> &R, 
												   Vec<int> &countWstarij, Vec<int> &countWstarji,
												   LevelSetStructure *LSS)
{

  R = 0.0;
  count = 0;

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {
    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

	 bool checkN = true;

    if(fluidId[i] != fluidId[j]) checkN = false;

	 if(LSS)
	 {
		 if(LSS->edgeIntersectsStructure(0.0, l)) checkN = false;

		 if(!LSS->isActive(0.0, i) && !LSS->isActive(0.0, j)) checkN = false;
	 }

	 if(checkN == false) continue; 

    double dx[3];
    dx[0] = X[j][0] - X[i][0];
    dx[1] = X[j][1] - X[i][1];
    dx[2] = X[j][2] - X[i][2];

    double dxdx = dx[0] * dx[0];
    double dydy = dx[1] * dx[1];
    double dzdz = dx[2] * dx[2];
    double dxdy = dx[0] * dx[1];
    double dxdz = dx[0] * dx[2];
    double dydz = dx[1] * dx[2];

	bool updateRi = false;
	bool updateRj = false;
	double locWeighti = 1.0;
	double locWeightj = 1.0;

	 if(checkN == true){
	  count[i][0]++;
	  count[j][0]++;
	  updateRi = true;
	  updateRj = true;
	}
	else {
	  if (LSS->isActive(0.0,i))
	    if (countWstarij[l]!=0) {
		  count[i][0]++;
		  LevelSetResult resij = LSS->getLevelSetDataAtEdgeCenter(0.0,l,true);
		  updateRi = true;
		  locWeighti = (1.0-resij.alpha)*(1.0-resij.alpha);
		}
	  if (LSS->isActive(0.0,j))
	    if (countWstarji[l]!=0) {
		  count[j][0]++;
		  LevelSetResult resji = LSS->getLevelSetDataAtEdgeCenter(0.0,l,false);
		  updateRj = true;
		  locWeightj = (1.0-resji.alpha)*(1.0-resji.alpha);
		}
	}
	if (updateRi) {
      R[i][0] += dxdx*locWeighti;
      R[i][1] += dxdy*locWeighti;
      R[i][2] += dxdz*locWeighti;
      R[i][3] += dydy*locWeighti;
      R[i][4] += dydz*locWeighti;
      R[i][5] += dzdz*locWeighti;
	}
	if (updateRj) {
      R[j][0] += dxdx*locWeightj;
      R[j][1] += dxdy*locWeightj;
      R[j][2] += dxdz*locWeightj;
      R[j][3] += dydy*locWeightj;
      R[j][4] += dydz*locWeightj;
      R[j][5] += dzdz*locWeightj;
	}
  }
}

//------------------------------------------------------------------------------

void SubDomain::computeWeightsLeastSquaresNodePart(SVec<double,6> &R)
{

  for (int i=0; i<R.size(); ++i) {

    double r11  = sqrt(R[i][0]);
    double or11 = 1.0 / r11;
    double r12  = R[i][1] * or11;
    double r13  = R[i][2] * or11;
    double r22  = sqrt(R[i][3] - r12*r12);
    double r23  = (R[i][4] - r12*r13) / r22;
    double r33  = sqrt(R[i][5] - (r13*r13 + r23*r23));

    R[i][0] = r11;
    R[i][1] = r12;
    R[i][2] = r13;
    R[i][3] = r22;
    R[i][4] = r23;
    R[i][5] = r33;

  }

}

//------------------------------------------------------------------------------
void SubDomain::computeWeightsLeastSquaresEdgePartForEmbeddedStruct(LevelSetStructure &LSS, 
								    SVec<double,3> &X, SVec<int,1> &count, 
																						  SVec<double,10> &R, Vec<int> &init, bool externalSI)
{

  R = 0.0;
  count = 0;
  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

	for(int l=0; l<edges.size(); ++l) 
	{
	if (!edgeFlag[l]) continue;
	int i = edgePtr[l][0];
	int j = edgePtr[l][1];

		bool validEdge = true;

		if(externalSI && LSS.edgeWithSI(l)) validEdge = false;

		if(LSS.edgeIntersectsStructure(0.0, l)) validEdge = false;

		if(!LSS.isActive(0.0, i) || !LSS.isActive(0.0, j)) validEdge = false;
		
		if(!validEdge) continue;
		
		double dx[3];
		dx[0] = X[j][0] - X[i][0];
		dx[1] = X[j][1] - X[i][1];
		dx[2] = X[j][2] - X[i][2];

	double dxdx = dx[0]*dx[0];
	double dydy = dx[1]*dx[1];
	double dzdz = dx[2]*dx[2];
	double dxdy = dx[0]*dx[1];
	double dxdz = dx[0]*dx[2];
	double dydz = dx[1]*dx[2];

		double dxcn = dx[0];
		double dycn = dx[1];
		double dzcn = dx[2];
		double cncn = 1.0;

		if((init[i]<1 && LSS.isActive(0.0,i)) && init[j]==1) 
		{
	  count[i][0]++;
	  R[i][0] += dxdx;  R[i][1] += dxdy;  R[i][2] += dxdz;  R[i][3] += dxcn;
	  R[i][4] += dydy;  R[i][5] += dydz;  R[i][6] += dycn;
	  R[i][7] += dzdz;  R[i][8] += dzcn;
	  R[i][9] += cncn;
	}

		if((init[j]<1 && LSS.isActive(0.0,j)) && init[i]==1) 
		{
	  count[j][0]++;
	  R[j][0] += dxdx;  R[j][1] += dxdy;  R[j][2] += dxdz;  R[j][3] -= dxcn;
	  R[j][4] += dydy;  R[j][5] += dydz;  R[j][6] -= dycn;
	  R[j][7] += dzdz;  R[j][8] -= dzcn;
	  R[j][9] += cncn;
	} 
  }
}

//------------------------------------------------------------------------------
void SubDomain::computeWeightsLeastSquaresNodePartForEmbeddedStruct(SVec<int,1> &count, SVec<double,10> &R) 
{

	  double r11, r12, r13, r14, r22, r23, r24, r33, r34, r44;
	  double or11, or22, or33;

	for(int i=0; i<R.size(); ++i) 
	{
		if (count[i][0] < 1) continue;

		if(count[i][0] > 3) 
		{
			if(R[i][0] <= 0.0) 
			{
        r11 = 0.0;  r12 = 0.0;  r13 = 0.0;  r14 = 0.0; 
		r22 = 0.0;  r23 = 0.0;  r24 = 0.0; 
		r33 = 0.0;  r34 = 0.0;
		r44 = 0.0;
      }
			else 
			{
		r11  = sqrt(R[i][0]);
		or11 = 1.0/r11;
		r12  = R[i][1]*or11;
		r13  = R[i][2]*or11;
		r14  = R[i][3]*or11;

				if((R[i][4] - r12*r12) <= 0.0)
				{
          r11 = 0.0;  r12 = 0.0;  r13 = 0.0;  r14 = 0.0; 
		  r22 = 0.0;  r23 = 0.0;  r24 = 0.0; 
		  r33 = 0.0;  r34 = 0.0;
		  r44 = 0.0;
        }
				else 
				{
		  r22  = sqrt(R[i][4]-r12*r12);
		  or22 = 1.0/r22;
		  r23  = (R[i][5]-r12*r13)*or22;
		  r24  = (R[i][6]-r12*r14)*or22;
					
					if((R[i][7]-(r13*r13 + r23*r23)) <= 0.0)
					{
            r11 = 0.0;  r12 = 0.0;  r13 = 0.0;  r14 = 0.0; 
		    r22 = 0.0;  r23 = 0.0;  r24 = 0.0; 
		    r33 = 0.0;  r34 = 0.0;
		    r44 = 0.0;
		  }
					else 
					{
			r33  = sqrt(R[i][7]-(r13*r13+r23*r23));
			or33 = 1.0/r33;
			r34  = (R[i][8]-(r13*r14+r23*r24))*or33;
			r44  = R[i][9]-(r14*r14+r24*r24+r34*r34);

						if(r44 <= 0.0)
						{
              r11 = 0.0;  r12 = 0.0;  r13 = 0.0;  r14 = 0.0; 
		      r22 = 0.0;  r23 = 0.0;  r24 = 0.0; 
		      r33 = 0.0;  r34 = 0.0;
		      r44 = 0.0;
			}
			else
			  r44 = sqrt(r44);
		  }
        }
      }
	  R[i][0] = r11;  R[i][1] = r12;  R[i][2] = r13;  R[i][3] = r14;
		  			  R[i][4] = r22;  R[i][5] = r23;  R[i][6] = r24;
	  								  R[i][7] = r33;  R[i][8] = r34;
	  												  R[i][9] = r44;
	}
		else 
		{
	  R[i][0] = 0.0;  R[i][1] = 0.0;  R[i][2] = 0.0;  R[i][3] = 0.0;
	 	 			  R[i][4] = 0.0;  R[i][5] = 0.0;  R[i][6] = 0.0;
	  								  R[i][7] = 0.0;  R[i][8] = 0.0;
	  												  R[i][9] = 0.0;
	}
  }
}

void SubDomain::computeWeightsLeastSquaresEdgePartForFF(LevelSetStructure *LSS, 
							SVec<double,3> &X, SVec<int,1> &count, SVec<double,10> &R, Vec<int> &init,
							Vec<int>& fluidId)
{
  R = 0.0;
  count = 0;
  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();
  for (int l=0; l<edges.size(); ++l) {
	if (!edgeFlag[l]) continue;
	int i = edgePtr[l][0];
	int j = edgePtr[l][1];
	if ((i==j)||(fluidId[i] != fluidId[j])) continue;
	double dx[3] = {X[j][0]-X[i][0],X[j][1]-X[i][1],X[j][2]-X[i][2]};
	double dxdx = dx[0]*dx[0];
	double dydy = dx[1]*dx[1];
	double dzdz = dx[2]*dx[2];
	double dxdy = dx[0]*dx[1];
	double dxdz = dx[0]*dx[2];
	double dydz = dx[1]*dx[2];
	double dxcn = dx[0]*1.0;
	double dycn = dx[1]*1.0;
	double dzcn = dx[2]*1.0;
	double cncn = 1.0*1.0;
	if ((init[i]<1/*&&LSS.isActive(0.0,i)*/)&&init[j]==1) {
//    bool toContinue;
//	toContinue = (((!LSS.isSwept(0.0,j))&&(LSS.isActive(0.0,j)))&&
//		((LSS.isSwept(0.0,i))&&(LSS.isActive(0.0,i))));
//	if (toContinue) {
	  count[i][0]++;
	  R[i][0] += dxdx;  R[i][1] += dxdy;  R[i][2] += dxdz;  R[i][3] += dxcn;
	  R[i][4] += dydy;  R[i][5] += dydz;  R[i][6] += dycn;
	  R[i][7] += dzdz;  R[i][8] += dzcn;
	  R[i][9] += cncn;
	}
	if ((init[j]<1/*&&LSS.isActive(0.0,j)*/)&&init[i]==1) {
//	toContinue = (((!LSS.isSwept(0.0,i))&&(LSS.isActive(0.0,i)))&&
//		((LSS.isSwept(0.0,j))&&(LSS.isActive(0.0,j))));
//	if (toContinue) {
	  count[j][0]++;
	  R[j][0] += dxdx;  R[j][1] += dxdy;  R[j][2] += dxdz;  R[j][3] -= dxcn;
	  R[j][4] += dydy;  R[j][5] += dydz;  R[j][6] -= dycn;
	  R[j][7] += dzdz;  R[j][8] -= dzcn;
	  R[j][9] += cncn;
	} 
  }
}

//------------------------------------------------------------------------------
void SubDomain::computeWeightsLeastSquaresNodePartForFF(SVec<int,1> &count, SVec<double,10> &R) {

  for (int i=0; i<R.size(); ++i) {
	if (count[i][0]<1) continue;
	if (count[i][0]>3) {
	  double r11, r12, r13, r14, r22, r23, r24, r33, r34, r44;
	  double or11, or22, or33;
      if (!(R[i][0]>0.0)) {
        r11 = 0.0;  r12 = 0.0;  r13 = 0.0;  r14 = 0.0; 
		r22 = 0.0;  r23 = 0.0;  r24 = 0.0; 
		r33 = 0.0;  r34 = 0.0;
		r44 = 0.0;
      }
	  else {
		r11  = sqrt(R[i][0]);
		or11 = 1.0/r11;
		r12  = R[i][1]*or11;
		r13  = R[i][2]*or11;
		r14  = R[i][3]*or11;
        if(!(R[i][4]-r12*r12>0.0)) {
          r11 = 0.0;  r12 = 0.0;  r13 = 0.0;  r14 = 0.0; 
		  r22 = 0.0;  r23 = 0.0;  r24 = 0.0; 
		  r33 = 0.0;  r34 = 0.0;
		  r44 = 0.0;
        }
		else {
		  r22  = sqrt(R[i][4]-r12*r12);
		  or22 = 1.0/r22;
		  r23  = (R[i][5]-r12*r13)*or22;
		  r24  = (R[i][6]-r12*r14)*or22;
		  if (!(R[i][7]-(r13*r13+r23*r23)>0.0)) {
            r11 = 0.0;  r12 = 0.0;  r13 = 0.0;  r14 = 0.0; 
		    r22 = 0.0;  r23 = 0.0;  r24 = 0.0; 
		    r33 = 0.0;  r34 = 0.0;
		    r44 = 0.0;
		  }
		  else {
			r33  = sqrt(R[i][7]-(r13*r13+r23*r23));
			or33 = 1.0/r33;
			r34  = (R[i][8]-(r13*r14+r23*r24))*or33;
			r44  = R[i][9]-(r14*r14+r24*r24+r34*r34);
			if (!(r44>0.0)) {
              r11 = 0.0;  r12 = 0.0;  r13 = 0.0;  r14 = 0.0; 
		      r22 = 0.0;  r23 = 0.0;  r24 = 0.0; 
		      r33 = 0.0;  r34 = 0.0;
		      r44 = 0.0;
			}
			else
			  r44 = sqrt(r44);
		  }
        }
      }
	  R[i][0] = r11;  R[i][1] = r12;  R[i][2] = r13;  R[i][3] = r14;
		  			  R[i][4] = r22;  R[i][5] = r23;  R[i][6] = r24;
	  								  R[i][7] = r33;  R[i][8] = r34;
	  												  R[i][9] = r44;
	}
	else {
	  R[i][0] = 0.0;  R[i][1] = 0.0;  R[i][2] = 0.0;  R[i][3] = 0.0;
	 	 			  R[i][4] = 0.0;  R[i][5] = 0.0;  R[i][6] = 0.0;
	  								  R[i][7] = 0.0;  R[i][8] = 0.0;
	  												  R[i][9] = 0.0;
	}
  }
}




//------------------------------------------------------------------------------

// Included (YC)
void SubDomain::computeDerivativeTransposeOfWeightsLeastSquaresNodePart(SVec<double,6> &R, SVec<double,6> &dR)
{
  for (int i=0; i<dR.size(); ++i) {

    double r11  = sqrt(R[i][0]);
    double or11 = 1.0 / r11;
    double r12  = R[i][1] * or11;
    double o12 = -R[i][1]*or11*or11/(2.0*r11);
    double o13 = - R[i][2]*or11*or11/(2.0*r11);
    double r13  = or11*R[i][2];
    double r22  = sqrt(R[i][3] - r12*r12);
    double o22 = -2.0*r12;
    double oo22 = 1.0/(2.0*r22);
    double r23  = (R[i][4] - r12*r13) / r22;
    double o23 = 2.0*oo22;
    double k23 = -o23*r23*oo22;
    double l23 = -o23*r12*or11;
    double q23 = -(o23*r13*or11+o23*r23*oo22*o22*or11);
    double p23 = -(o23*r13*o12+o23*r12*o13+o23*r23*oo22*o22*o12);
    double r33  = sqrt(R[i][5] - (r13*r13 + r23*r23));
    double o33 = 1.0/(2.0*r33);
    double p33 = -2.0*o33*r23;
    double k33 = p33*o23;
    double w33 = p33*k23;
    double u33 = p33*l23 - 2.0*o33*r13*or11;
    double m33 = p33*q23;
    double c33 = p33*p23 - 2.0*o33*r13*o13;

    double dr0 = dR[i][0];
    double dr1 = dR[i][1];
    double dr2 = dR[i][2];
    double dr3 = dR[i][3];
    double dr4 = dR[i][4];
    double dr5 = dR[i][5];

    dR[i][0] = 1.0/(2.0*r11)*dr0 + o12*dr1 + o13*dr2 + oo22*o22*o12*dr3 + p23*dr4 + c33*dr5;
    dR[i][1] = or11*dr1 + oo22*o22*or11*dr3 + q23*dr4 + m33*dr5; 
    dR[i][2] = or11*dr2 + l23*dr4 + u33*dr5; 
    dR[i][3] = oo22*dr3 + k23*dr4 + w33*dr5; 
    dR[i][4] = o23*dr4 + k33*dr5; 
    dR[i][5] = o33*dr5;

  }
}

//------------------------------------------------------------------------------
// Included (YC)
void SubDomain::compute_dRdR(SVec<double,6> &R, RectangularSparseMat<double,6,6> &dRdR)
{
  for (int i=0; i<R.size(); ++i) {

    double r11  = sqrt(R[i][0]);
    double or11 = 1.0 / r11;
    double r12  = R[i][1] * or11;
    double o12 = -R[i][1]*or11*or11/(2.0*r11);
    double o13 = - R[i][2]*or11*or11/(2.0*r11);
    double r13  = or11*R[i][2];
    double r22  = sqrt(R[i][3] - r12*r12);
    double o22 = -2.0*r12;
    double oo22 = 1.0/(2.0*r22);
    double r23  = (R[i][4] - r12*r13) / r22;
    double o23 = 2.0*oo22;
    double k23 = -o23*r23*oo22;
    double l23 = -o23*r12*or11;
    double q23 = -(o23*r13*or11+o23*r23*oo22*o22*or11);
    double p23 = -(o23*r13*o12+o23*r12*o13+o23*r23*oo22*o22*o12);
    double r33  = sqrt(R[i][5] - (r13*r13 + r23*r23));
    double o33 = 1.0/(2.0*r33);
    double p33 = -2.0*o33*r23;
    double k33 = p33*o23;
    double w33 = p33*k23;
    double u33 = p33*l23 - 2.0*o33*r13*or11;
    double m33 = p33*q23;
    double c33 = p33*p23 - 2.0*o33*r13*o13;
    double o00 = 1.0/(2.0*r11); 
    double o14 = oo22*o22*o12;
    double or22 = oo22*o22*or11;

    double C[6][6] = {0};
    C[0][0] = o00;
    C[1][0] = o12;  C[1][1] = or11;
    C[2][0] = o13;                   C[2][2] = or11;
    C[3][0] = o14;  C[3][1] = or22;                   C[3][3] = oo22;
    C[4][0] = p23;  C[4][1] = q23;   C[4][2] = l23;   C[4][3] = k23;  C[4][4] = o23;
    C[5][0] = c33;  C[5][1] = m33;   C[5][2] = u33;   C[5][3] = w33;  C[5][4] = k33;  C[5][5] = o33;

    int ndList[1] = {i};
    dRdR.addContrib(1, ndList, C[0]);
  }
}
//------------------------------------------------------------------------------

// Included (YC)
void SubDomain::computeDerivativeOfWeightsLeastSquaresNodePart(RectangularSparseMat<double,6,6> &dRdR, SVec<double,6> &dR)
{

  SVec<double,6> dr(dR); 
  dRdR.apply(dr, dR);

}

//------------------------------------------------------------------------------

// Included (YC)
void SubDomain::computeTransposeDerivativeOfWeightsLeastSquaresNodePart(RectangularSparseMat<double,6,6> &dRdR, SVec<double,6> &dR)
{

  SVec<double,6> dr(dR); 
  dRdR.applyTranspose(dr, dR);
}

//------------------------------------------------------------------------------

// Included (MB)
void SubDomain::computeDerivativeOfWeightsLeastSquaresNodePart(SVec<double,6> &R, SVec<double,6> &dR)
{
  for (int i=0; i<dR.size(); ++i) {
    double r11  = sqrt(R[i][0]);
    double dr11  = 1.0/(2.0*r11)*dR[i][0];
    double or11 = 1.0 / r11;
    double dor11 = -1.0 /(r11* r11)*dr11;
    double r12  = R[i][1] * or11;
    double dr12  = dR[i][1] * or11 + R[i][1] * dor11;
    double r13  = R[i][2] * or11;
    double dr13  = dR[i][2] * or11 + R[i][2] * dor11;
    double r22  = sqrt(R[i][3] - r12*r12);
    double dr22  = 1.0/(2.0*r22)*(dR[i][3] - 2.0*r12*dr12);
    double r23  = (R[i][4] - r12*r13) / r22;
    double dr23  = ( (dR[i][4] - dr12*r13 - r12*dr13)*r22 - (R[i][4] - r12*r13)*dr22 ) / (r22*r22);
    double r33  = sqrt(R[i][5] - (r13*r13 + r23*r23));
    double dr33  = 1.0/(2.0*r33)*(dR[i][5] - 2.0*(r13*dr13 + r23*dr23));

    dR[i][0] = dr11;
    dR[i][1] = dr12;
    dR[i][2] = dr13;
    dR[i][3] = dr22;
    dR[i][4] = dr23;
    dR[i][5] = dr33;
  }
}

//------------------------------------------------------------------------------
//d2d$
void SubDomain::computeWeightsLeastSquaresNodePart(SVec<int,1> &count, SVec<double,6> &R)
{

      double r11, or11, r12, r13, r22, r23, r33;

	for(int i=0; i<R.size(); ++i) 
	{		
		if(count[i][0]>3) 
		{ 

			/* -----------------------------------------------------------------*/
			if(R[i][0] <= 0.0)
			{
        r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0; r23 = 0.0; r33 = 0.0;
			}
			else
			{
        r11  = sqrt(R[i][0]);
        or11 = 1.0 / r11;
        r12  = R[i][1] * or11;
        r13  = R[i][2] * or11;

				if((R[i][3] - r12*r12) <= 0.0)
				{
          r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0; r23 = 0.0; r33 = 0.0;
				}
				else
				{
          r22  = sqrt(R[i][3] - r12*r12);
          r23  = (R[i][4] - r12*r13) / r22;

					if((R[i][5] - (r13*r13 + r23*r23)) < 0.0)
					{
            r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0; r23 = 0.0; r33 = 0.0;
					}
					else
						r33 = sqrt(R[i][5] - (r13*r13 + r23*r23));
				}
				
			}
			/* -----------------------------------------------------------------*/

      R[i][0] = r11;
      R[i][1] = r12;
      R[i][2] = r13;
      R[i][3] = r22;
      R[i][4] = r23;
      R[i][5] = r33;
		}
		else
		{
      R[i][0] = 0.0;
      R[i][1] = 0.0;
      R[i][2] = 0.0;
      R[i][3] = 0.0;
      R[i][4] = 0.0;
      R[i][5] = 0.0;
    }

  }

}

//------------------------------------------------------------------------------

void SubDomain::computeWeightsGalerkin(SVec<double,3> &X, SVec<double,3> &wii,
				       SVec<double,3> &wij, SVec<double,3> &wji)
{

  wii = 0.0;
  wij = 0.0;
  wji = 0.0;

  for (int i=0; i<elems.size(); ++i){
    elems[i].computeWeightsGalerkin(X, wii, wij, wji);
  }

}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void SubDomain::computeWeightsGalerkin(SVec<double,3> &X, const Vec<int> &fluidId,
				       SVec<double,3> &wii,
				       SVec<double,3> &wij, 
				       SVec<double,3> &wji,
				       LevelSetStructure *LSS, bool includeSweptNodes)
{

  wii = 0.0;
  wij = 0.0;
  wji = 0.0;

  //d2d
  bool e_flag = false;

  for(int i=0; i<elems.size(); ++i)
  {	  
	  for(int j=0; j<6; j++)
	  {

      int le = elems[i].edgeNum(j);
      int Nj = elems[i].nodeNum( elems[i].edgeEnd(j,0) );
      int Nk = elems[i].nodeNum( elems[i].edgeEnd(j,1) );

		if(!LSS)
		{
			if(fluidId[Nj] != fluidId[Nk])
			{
	e_flag = true;
	break;
      }
		}
		else
		{
			if( (!LSS->isActive(0.0, Nj) || !LSS->isActive(0.0, Nk)) || LSS->edgeIntersectsStructure(0.0, le) )
			{
	e_flag = true;
	break;
      }

			if(!includeSweptNodes && (LSS->isSwept(0.0,Nj) || LSS->isSwept(0.0,Nk)))
			{
				e_flag = true;
				break;
			}

		}
 
    }

    if(!e_flag) elems[i].computeWeightsGalerkin(X, wii, wij, wji);

  }

}


//------------------------------------------------------------------------------


// Included (MB)
void SubDomain::computeDerivativeOfWeightsGalerkin(SVec<double,3> &X, SVec<double,3> &dX, SVec<double,3> &dwii, SVec<double,3> &dwij, SVec<double,3> &dwji)
{

  dwii = 0.0;
  dwij = 0.0;
  dwji = 0.0;

  for (int i=0; i<elems.size(); ++i)
    elems[i].computeDerivativeOfWeightsGalerkin(X, dX, dwii, dwij, dwji);

}

//------------------------------------------------------------------------------

// Included (YC)
void SubDomain::computeDerivativeTransposeOfWeightsGalerkin(SVec<double,3> &X, SVec<double,3> &dwii, SVec<double,3> &dwij, SVec<double,3> &dwji, SVec<double,3> &dX)
{

  dX = 0.0;

  for (int i=0; i<elems.size(); ++i) {
    elems[i].computeDerivativeTransposeOfWeightsGalerkin(X, dwii, dwij, dwji, dX);
  }

}


//------------------------------------------------------------------------------

void SubDomain::computeEdgeWeightsGalerkin(SVec<double,3> &X, Vec<double> &A, SVec<double,9> &M)
{

  M = 0.0;

  for (int i=0; i<elems.size(); ++i)
    elems[i].computeEdgeWeightsGalerkin(X, M);

  M *= -1.0/9.0;

  if (globSubNum == 0) {

    Vec<bool> tmp(edges.size());
    tmp = true;

    for (int iSub = 0; iSub < numNeighb; ++iSub)
      for (int iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge)
	tmp[ sharedEdges[iSub][iEdge].edgeNum ] = false;

    Vec<double> u(nodes.size());
    SVec<double,9> duddxj(nodes.size());

    int i, j, k;
    for (i=0; i<nodes.size(); ++i)
      u[i] = 2.0*X[i][0];// *X[i][0];

    duddxj = 0.0;

    int (*edgePtr)[2] = edges.getPtr();

    for (int l=0; l<edges.size(); ++l) {

      i = edgePtr[l][0];
      j = edgePtr[l][1];

      double du = u[j] - u[i];
      for (k=0; k<9; ++k) {
	duddxj[i][k] += M[l][k] * du;
	duddxj[j][k] -= M[l][k] * du;
      }

      if (nodeType[i] == BC_INTERNAL) {
	for (k=0; k<9; ++k)
	  fprintf(stderr, "%d %e\n", k, M[l][k]);
	exit(1);
      }

      /*
      int ff = 1;

      if (!tmp[l])
	ff = 0;

      fprintf(stderr, "%d %e %e %e %d %d (%d)\n",
	      l, M[l][1], M[l][6], M[l][8], nodeType[i], nodeType[j], ff);
      */

    }

    for (i=0; i<nodes.size(); ++i) {
      double a = 1.0 / A[i];
      for (k=0; k<9; ++k)
	duddxj[i][k] *= a;

      fprintf(stderr, "%d %e %d\n", i, duddxj[i][0], nodeType[i]);

    }

    exit(1);

  }

}

//------------------------------------------------------------------------------

void SubDomain::computeFilterWidth(SVec<double,3> &X, Vec<double> &Delta)
{
  double third = 1.0/3.0;

  for (int elemNum=0; elemNum < elems.size(); ++elemNum) {
    double vol = elems[elemNum].computeVolume(X);
    double delta = pow(vol, third);
    for (int i=0; i<4; ++i)
      Delta[elems[elemNum][i]] += delta * vol;
  }
}

//------------------------------------------------------------------------------
void SubDomain::computeLocalAvg(SVec<double,3> &X, Vec<double> &Q, Vec<double> &W)
{
  double fourth = 1.0/4.0;

	for(int tetNum=0; tetNum < elems.size(); ++tetNum) 
	{
    double vol = elems[tetNum].computeVolume(X);
    int idx[4] = {elems[tetNum][0], elems[tetNum][1],
                  elems[tetNum][2], elems[tetNum][3]};
    double Qcg = fourth * (Q[idx[0]] + Q[idx[1]] + Q[idx[2]] + Q[idx[3]]);

		for(int i=0; i<4; ++i) W[idx[i]] += Qcg * vol;
  }

}

//------------------------------------------------------------------------------
void SubDomain::applySmoothing(Vec<double> &ctrlVol, Vec<double> &Q)
{

	for(int i=0; i<nodes.size(); ++i) 
	{
    double coef = 1.0 / (4.0*ctrlVol[i]);
    Q[i] = coef * Q[i];
 }
}

//------------------------------------------------------------------------------

void SubDomain::computeLocalAvg(SVec<double,3> &X, SVec<double,2> &Q, SVec<double,2> &W)
{

  double fourth = 1.0/4.0;

  for (int tetNum=0; tetNum < elems.size(); ++tetNum) {
    double vol = elems[tetNum].computeVolume(X);
    int idx[4] = {elems[tetNum][0], elems[tetNum][1],
                  elems[tetNum][2], elems[tetNum][3]};
    double Qcg[2];
    Qcg[0] = fourth * (Q[idx[0]][0] + Q[idx[1]][0] + Q[idx[2]][0] + Q[idx[3]][0]);
    Qcg[1] = fourth * (Q[idx[0]][1] + Q[idx[1]][1] + Q[idx[2]][1] + Q[idx[3]][1]);

    for (int i=0; i<4; ++i) {
//      if (nodeType[idx[i]] == BC_ADIABATIC_WALL_MOVING  || nodeType[idx[i]] == BC_ADIABATIC_WALL_FIXED  ||
//          nodeType[idx[i]] == BC_ISOTHERMAL_WALL_MOVING || nodeType[idx[i]] == BC_ISOTHERMAL_WALL_FIXED) {
//            W[idx[i]][0] += Q[idx[i]][0] * vol;
//            W[idx[i]][1] += Q[idx[i]][1] * vol;   // essentially no smoothing on the boundary
//      }
//      else  {
        W[idx[i]][0] += Qcg[0] * vol;
        W[idx[i]][1] += Qcg[1] * vol;
//      }
    }

  }

}

//------------------------------------------------------------------------------

void SubDomain::applySmoothing(Vec<double> &ctrlVol, SVec<double,2> &Q)
{


 for (int i=0; i<nodes.size(); ++i) {
    double coef = 1.0 / (4.0*ctrlVol[i]);
    Q[i][0] = coef * Q[i][0];
    Q[i][1] = coef * Q[i][1];
 }


}

//------------------------------------------------------------------------------

/*

//---------------------------------------//
// This is the old agglomeration routine //
//---------------------------------------//

MacroCellSet* SubDomain::findAgglomerateMesh(int scopeWidth,
                                             bool *masterFlag,
                                             double gam)
{

  int numNodes        = nodes.size();
  int numEdges        = edges.size();
  int (*ptrEdge)[2]   = edges.getPtr();
  int l, m;

  // numNeighbors holds number of nodes neighboring each node
  // neighbors holds node identifiers corresponding to neighbors of each node
  int *numNeighbors = new int  [numNodes];
  int **neighbors   = new int *[numNodes];

  // consideredNodes is an array of booleans with one element per node
  // True  --> node already considered for inclusion in a macro-cell
  // False --> node not yet considered for inclusion in a macro-cell
  // Note that a value of true does not imply that the node has been included in
  // a macro-cell, only that it has been considered.  Nodes for which
  // the masterFlag is false will be considered but not included.
  bool* consideredNodes = new bool [numNodes];

  // Loop over all nodes in subdomain to initialize numNeighbors
  // and specify all nodes as not yet considered for any macro-cells
  for (m=0; m<numNodes; ++m) {
    consideredNodes[m] = false;
    numNeighbors[m]  = 0;
  }

  // Count up the number of neighbors each node has
  for (l=0; l<numEdges; ++l) {
    numNeighbors[ ptrEdge[l][0] ]++;
    numNeighbors[ ptrEdge[l][1] ]++;
  }

  // Allocate the size of neighbors based upon numNeighbors
  for (m=0; m<numNodes; ++m) {
    neighbors[m]    = new int[ numNeighbors[m] ];
    numNeighbors[m] = 0;
  }
  // Fill the list of neighbors
  for (l=0; l<numEdges; ++l) {
    neighbors[ ptrEdge[l][0] ][ numNeighbors[ ptrEdge[l][0] ] ] = ptrEdge[l][1];
    neighbors[ ptrEdge[l][1] ][ numNeighbors[ ptrEdge[l][1] ] ] = ptrEdge[l][0];
    numNeighbors[ ptrEdge[l][0] ]++;
    numNeighbors[ ptrEdge[l][1] ]++;
  }

  // The j-th element of nodeToMacroCellMap contains the ID of the
  // macro-cell containing the jth node (cell).  numMacroCells will
  // hold the total number of macro-cells to be created.
  int* nodeToMacroCellMap = new int [numNodes];
  int numMacroCells       = 0;

  // Loop over all nodes to determine the location and number of the macro-cells
  for (m=0; m<numNodes; ++m) {

    // Only consider nodes that are masters of their subdomains
    if (masterFlag[m]) {

      // Continue only if node hasn't yet been considered for inclusion
      if (!consideredNodes[m]) {

        consideredNodes[m]    = true;           // This node has now been considered
        nodeToMacroCellMap[m] = numMacroCells;  // add this node to the current macro-cell

        // Fill the new macro-cell with the not-yet-included neighbors.
        // scopeWidth determines the level of recursion
        assimilateCells(m, scopeWidth, numNeighbors, neighbors, masterFlag,
                        numMacroCells, nodeToMacroCellMap, consideredNodes);

        // Go to the next macro-cell
        numMacroCells++;
      }
    }

    // If node isn't a master of its subdomain, specify as considered, but
    // assign -1 as the macro-cell to indicate node will be included via
    // a neighboring subdomain.
    else {
      nodeToMacroCellMap[m] = -1;
      consideredNodes[m]    = true;
    }
  }

  // Create set of macro-cells based upon the relevant information
  MacroCellSet* macroCells = new MacroCellSet(numMacroCells,
                                              nodeToMacroCellMap,
                                              numNodes,
                                              gam);

  // Delete temporary arrays
  for (m=0; m<numNodes; m++) {
    if (neighbors[m]) delete[] neighbors[m];
  }

  if (numNeighbors) delete [] numNeighbors;
  if (consideredNodes) delete [] consideredNodes;
  if (nodeToMacroCellMap) delete [] nodeToMacroCellMap;
  if (neighbors) delete [] neighbors;

  return macroCells;

}
*/

//------------------------------------------------------------------------------

// This function assimilates cells into a macro-cell and may be
// called recursively depending upon the algorithm

void SubDomain::assimilateCells(int baseNode,
                                int scopeWidth,
                                int *numNeighbors,
                                int **neighbors,
                                bool *masterFlag,
                                int macroCellID,
                                int *nodeToMacroCellMap,
                                bool *consideredNodes)
{

  int n;

  // Loop over nodes neighboring baseNode
  for (n=0; n < numNeighbors[baseNode]; ++n) {

    // Only consider nodes that are masters of their subdomains
    if (masterFlag[ neighbors[baseNode][n] ]) {

      // Test for inclusion of each of baseNode's neighbors in a macro-cell
      // and if not included, add the neighbor to current macro-cell and
      // specify that the node has been considered
      if (!consideredNodes[ neighbors[baseNode][n] ]) {
        consideredNodes[ neighbors[baseNode][n] ]    = true;
        nodeToMacroCellMap[ neighbors[baseNode][n] ] = macroCellID;
      }
    }

    // If node isn't a master of its subdomain, specify as considered, but
    // assign -1 as the macro-cell to indicate node will be included via
    // a neighboring subdomain.
    else {
      consideredNodes[ neighbors[baseNode][n] ]    = true;
      nodeToMacroCellMap[ neighbors[baseNode][n] ] = -1;
    }

    // If the recursion depth is greater than one, then examine the
    // neighbors of the neighbors
    if (scopeWidth>1) {
      assimilateCells(neighbors[baseNode][n], scopeWidth-1, numNeighbors,
                      neighbors, masterFlag, macroCellID, nodeToMacroCellMap,
                      consideredNodes);
    }
  }
}

//------------------------------------------------------------------------------

/*
Connectivity *SubDomain::agglomerate(Connectivity &nToN, int maxLevel, bool *masterFlag)
{
  // Get the number of "nodes"
  int nNodes = nToN.csize();
  int *list = new int[nNodes];
  // Initially we do not know how many cells we will end up with, but at most
  // we will have nNodes
  int *ptr = new int[nNodes];
  int *level = new int[nNodes];
  int i, j, k;

  for(i = 0; i < nNodes; ++i)
    level[i] = -1;

  int nCells = 0;
  int examine;
  int idx = 0;

//  printf("I am here\n");

  // Now find a seed point (random)
  for(examine = 0; examine < nNodes; ++examine) {
     if(level[examine] >= 0 || masterFlag[examine] == false)
       continue;
     // Take examine as the start point
     list[idx] = examine;
     level[examine] = 0;
     int trailer = idx;
     int header = idx+1;
     idx++;
     for(int l = 0; l < maxLevel; ++l) {
       while(trailer < header) {
          for(j = 0; j < nToN.num(list[trailer]); ++j) {
	    int nd = nToN[list[trailer]][j];
	    if(level[nd] < 0 && masterFlag[nd]) {
	      level[nd] = l+1;
	      list[idx++] = nd;
	    }
	  }
          trailer++;
       }
       header = idx;
     }
     ptr[nCells] = idx;
     nCells++;
  }

  // Now create the REAL pointers
  int *rPtr = new int[nCells+1];
  rPtr[0] = 0;
  for(i = 0; i < nCells;++i)
   rPtr[i+1] = ptr[i];
  delete [] level;
  delete [] ptr;
  return new Connectivity(nCells, rPtr, list);
}
*/

//------------------------------------------------------------------------------

// This is the driver routine for the new agglomeration routine

Connectivity *SubDomain::agglomerate(Connectivity &nToN, int maxLevel, bool *masterFlag)
{
  int nNodes = nToN.csize();     // Get the number of "nodes"
  int *list = new int[nNodes];

  // Initially we do not know how many cells we will end up with, but at most
  // we will have nNodes

  int *ptr = new int[nNodes+1];
  int *level = new int[nNodes];    // tag denoting if the node has been examined or not
  int *cellNum = new int[nNodes];  // tag that contains the cell to which a node belongs
  int i, j, k;

  for(i = 0; i < nNodes; ++i)
    level[i] = -1;  // initialize the tag that denotes if the node has been examined  to not examined

  int nCells = 0;   // counter that registers the number of macrocells
  int examine;
  int idx = 0;

  ptr[0] = 0;       // initialize the pointer to the start of the macrocell-1 to 0

  // Now find a seed point
  // (random - actually depends on the data structure and hence is not random in real sense)
  for(examine = 0; examine < nNodes; ++examine) {
     if(level[examine] >= 0 || masterFlag[examine] == false)
       continue;         // if the node has been examined (or)
                         // if the master flag of the node is false, then skip

     // Take examine as the start point for creating the next macro cell
     list[idx] = examine;
     cellNum[examine] = nCells;
     level[examine] = 0;
     int trailer = idx;
     int header = idx+1;
     idx++;

     for(int l = 0; l < maxLevel; ++l) {  // loop till the specified number of neighbours have been agglomerated
       while(trailer < header) {
          for(j = 0; j < nToN.num(list[trailer]); ++j) {  // add neighbours of the examined node
            int nd = nToN[list[trailer]][j];
            if(level[nd] < 0 && masterFlag[nd]) {
              level[nd] = l+1;
              list[idx++] = nd;
              cellNum[nd] = nCells;
            }
          }
          trailer++;
       }
       header = idx;
     }
     ptr[nCells+1] = idx;  // point to the start of the next macrocell
     nCells++;
  }

  // Examine all macrocells and deal with the ones that are too small.
  float ratio = 1/2.2;  // can be changed (just heurisic)
  int maxSize = 0;
  int iCell;

  // finding out the number of nodes in the largest macrocell
  for(iCell = 0; iCell < nCells; ++iCell)
     maxSize = std::max(maxSize, ptr[iCell+1]-ptr[iCell]);
  int acceptableNumber = int(ratio*maxSize+0.5);  // acceptable number is a percentage of the number
                                                  // of nodes in the biggest macro cell
  int nRedist = 0;
  bool *isBigEnough = new bool[nCells];

  // Build the list of nodes that need to be redistributed

  // The idea is if the number of nodes in a macrocell is smaller than a  //
  // percentage of number of nodes in the biggest macrocell then all the  //
  // nodes are taken up to be redistributed                               //

  // Phase 1: count how many there are
  for(iCell = 0; iCell < nCells; ++iCell) {
    if(ptr[iCell+1]-ptr[iCell] < acceptableNumber) {
      isBigEnough[iCell] = false;
      nRedist += ptr[iCell+1]-ptr[iCell];
    } else
      isBigEnough[iCell] = true;
  }
  int *redistList = new int[nRedist];  // create the redistribution list array
  nRedist = 0;

  // Phase 2: build the actual list
  // ie.. the resdistList array is now filled
  for(iCell = 0; iCell < nCells; ++iCell)
    if(!isBigEnough[iCell]) {
      for(i = ptr[iCell]; i < ptr[iCell+1]; ++i)
         redistList[nRedist++] = list[i];
    }
  int *cellCount = new int[nCells];
  int iter = 0;
  int threshold = 3;
  int lastNRedist = nRedist+1;
  while(nRedist > 0 && (nRedist != lastNRedist || threshold > 1)) {
    if(nRedist == lastNRedist)
      threshold--;
    lastNRedist = nRedist;
    for(i = nRedist; i-- > 0; ) {
      int nd = redistList[i];
      for(j = 0; j < nToN.num(nd); ++j)
        if (masterFlag[nToN[nd][j]])
          cellCount[cellNum[nToN[nd][j]]] = 0;
      int maxCount = 0;
      int bestCell = -1;
      for(j = 0; j < nToN.num(nd); ++j) {
        if(!masterFlag[nToN[nd][j]]) continue;
        int cell = cellNum[nToN[nd][j]];
        cellCount[cell]++;
        if(isBigEnough[cell] && cellCount[cell] > maxCount) {
           maxCount = cellCount[cell];
           bestCell = cell;
        }
      }
      if(maxCount >= threshold) {
        cellNum[nd] = bestCell;
        redistList[i] =redistList[--nRedist];
      }
    }
  }

  // Phase 3: Now create the REAL pointers
  // First remap the cell number to remove the empty cells.
  int cellMap = 0;

  // count how many nodes are in each cell.
  for(iCell = 0; iCell < nCells; ++iCell)
    ptr[iCell] = 0;
  for(i = 0; i < nNodes; ++i)
    if (masterFlag[i])
      ptr[cellNum[i]]++;
  for(iCell = 0; iCell < nCells; ++iCell)
    if(ptr[iCell] > 0)
      level[iCell] = cellMap++;
  for(i = 0; i < nNodes; ++i)
    if (masterFlag[i])
     cellNum[i] = level[cellNum[i]];

  delete [] level;
  delete [] ptr;

  nCells = cellMap;
  ptr = new int[nCells+1];
  for(i = 0; i <= nCells; ++i)
    ptr[i] = 0;
  for(i = 0; i < nNodes; ++i)
    if (masterFlag[i])
      ptr[cellNum[i]]++;
  for(i = 1; i <= nCells; ++i) {
    ptr[i] += ptr[i-1];
  }
  for(i = 0; i < nNodes; ++i)
    if (masterFlag[i])
      list[--ptr[cellNum[i]]] = i;

  return new Connectivity(nCells, ptr, list);

}

//------------------------------------------------------------------------------
// This is the new and correct agglomeration routine

MacroCellSet** SubDomain::findAgglomerateMesh(int scopeWidth, int scopeDepth,
                                              bool *masterFlag, double gam)
{
  int l, m, n;
  int numNodes        = nodes.size();

  MacroCellSet** cells = new MacroCellSet*[scopeDepth];

  // This is redoing work already done, but oh well.... //

  Connectivity *nodeToNode = createEdgeBasedConnectivity();

  // numNeighbors holds number of nodes neighboring each node //
  // neighbors holds node identifiers corresponding to neighbors of each node //

  Connectivity *cellToNode = agglomerate(*nodeToNode, scopeWidth, masterFlag);

  // Create set of macro-cells based upon the relevant information //

    cells[0] = new MacroCellSet(cellToNode->csize(),
                                              cellToNode,
                                              numNodes,
                                              gam);
    for (n = 1; n<scopeDepth; ++n) {
      Connectivity *nodeToCell = cellToNode->reverse();
      Connectivity *cellToNode2 = cellToNode->transcon(nodeToNode);
      Connectivity *cellToCell = cellToNode2->transcon(nodeToCell);
      bool *cellMastFlag = new bool[cellToCell->csize()];
      for(l = 0; l < cellToCell->csize(); ++l)
        cellMastFlag[l] = true;
      Connectivity *macroCellToCell = agglomerate(*cellToCell, scopeWidth, cellMastFlag);
      Connectivity *macroCellToNode = macroCellToCell->transcon(cellToNode);

      cells[n] = new MacroCellSet(macroCellToNode->csize(),
                                              macroCellToNode,
                                              numNodes,
                                              gam);
      cellToNode = macroCellToNode;
    }

//  printMacroCells(*cellToNode); exit(1);

  delete nodeToNode;
  delete cellToNode;  // If cell to node becomes the basis for MacroCellSet, remove this line
  return cells;

}

//------------------------------------------------------------------------------
void SubDomain::createMacroCellConnectivities(MacroCellSet **macroCells, int scopeDepth)
{
  nodesToMCNodes = new Connectivity*[scopeDepth];
  for (int i=0; i<scopeDepth; ++i)
    nodesToMCNodes[i] = createNodeToMacroCellNodeConnectivity(macroCells[i]);
}

//------------------------------------------------------------------------------
extern double orient3d(double*, double*, double*, double*);
extern double orient3dslow(double*, double*, double*, double*);

inline
int checkPiercePointInTriangle(double vol1, double vol2, double vol3)
{

  if ((vol1 > 0.0 && vol2 > 0.0 && vol3 > 0.0) ||
      (vol1 < 0.0 && vol2 < 0.0 && vol3 < 0.0))
    return 0;
  else if (vol1 == 0.0) {
    if (vol2*vol3 > 0.0)
      return 1;
    else if (vol2 == 0)
      return -2;
    else if (vol3 == 0)
      return -1;
    else
      return 14;
  }
  else if (vol2 == 0.0) {
    if (vol1*vol3 > 0.0)
      return 2;
    else if (vol1 == 0)
      return -2;
    else if (vol3 == 0)
      return -3;
    else
      return 15;
  }
  else if (vol3 == 0.0) {
    if (vol1*vol2 > 0.0)
      return 3;
    else if (vol1 == 0)
      return -1;
    else if (vol2 == 0)
      return -3;
    else
      return 16;
  }

  return 17;

}

//------------------------------------------------------------------------------

inline
int checkPiercePoint(Vec3D& a, Vec3D& b, Vec3D& c, Vec3D& d, Vec3D& e, bool adaptive)
{
        //note that d and e are most probably not on the same side of plane abc in the
        //code (even though they are physically) since we imposed a factor 1000 in the
        //findTetrahedra routine!!!!
        // the factor 1000 is an arbitrary choice, but this factor should be enough
        //if we assume the tetrahedra are "homogeneous" to a certain degree

	double tol = 1.0e-12;
  
  double (*orient)(double*, double*, double*, double*);
  if (adaptive)
    orient = &orient3d;
  else
    orient = &orient3dslow;

  double vol1 = orient(a, b, c, d);
  double vol2 = orient(a, b, c, e);
  double vol3;

	//d2d
	if(fabs(vol1) <= tol) vol1 = 0.0; 
	if(fabs(vol2) <= tol) vol2 = 0.0;
	//d2d

  if (vol1 * vol2 > 0.0) return 10;          // [e,d] doesn't cross the plane
  if (vol1 == 0.0 && vol2 == 0.0) return 11; // [e,d] is in the plane
  if (vol1 == 0.0 && vol2 != 0.0) return 12; // d touches the plane
  if (vol1 != 0.0 && vol2 == 0.0) return 13; // e touches the plane

  vol1 = orient(d, a, b, e);
  vol2 = orient(d, b, c, e);
  vol3 = orient(d, c, a, e);

   //d2d
	if(fabs(vol1) < tol) vol1 = 0.0;
	if(fabs(vol2) < tol) vol2 = 0.0;
	if(fabs(vol3) < tol) vol3 = 0.0;
   //d2d
  
  int pierce = checkPiercePointInTriangle(vol1, vol2, vol3);

  return pierce;

}

//------------------------------------------------------------------------------
// plane = c + r * (a-c) + t * (b-c)
// this routine computes the barycentric coordinates r and t.
// however, it should also be able to compute the barycentric coordinates s, r and t
//   but it seems that s is not properly implemented...?

inline
void constructPiercePoint(Vec3D& a, Vec3D& b, Vec3D& c, Vec3D& d, Vec3D& e,
                          double* s, double* r, double* t, bool test)
{

    Vec3D ac = a - c;
    Vec3D cb = c - b;
    Vec3D cd = c - d;
    Vec3D ed = e - d;

    double coef0 = ac[2]*cb[1] - ac[1]*cb[2];
    double coef1 = ac[2]*cb[0] - ac[0]*cb[2];
    double coef2 = ac[1]*cb[0] - ac[0]*cb[1];

    double den = ed[0]*coef0 - ed[1]*coef1 + ed[2]*coef2;
    if (den == 0.0) {
      if (test){
        fprintf(stderr, "*** Error: vanishing denominator\n");
        exit(1);
      }
    }else{
    double ooden = 1.0 / den;

    if (s)
      *s = ooden * (cd[0]*coef0 - cd[1]*coef1 + cd[2]*coef2);

    if (r)
      *r = ooden * (ed[0] * (cb[2]*cd[1] - cb[1]*cd[2]) -
                    ed[1] * (cb[2]*cd[0] - cb[0]*cd[2]) +
                    ed[2] * (cb[1]*cd[0] - cb[0]*cd[1]));

    if (t)
      *t = ooden * (ed[0] * (ac[2]*cd[1] - ac[1]*cd[2]) -
                    ed[1] * (ac[2]*cd[0] - ac[0]*cd[2]) +
                    ed[2] * (ac[1]*cd[0] - ac[0]*cd[1]));
    }
}

//------------------------------------------------------------------------------
#ifdef EDGE_LENGTH //HB
bool SubDomain::findTetrahedron(int i, int j, Vec<int>& count, int** list,
				SVec<double,3>& X, V6NodeData& v6data, bool adaptive, double* refLength)
#else
bool SubDomain::findTetrahedron(int i, int j, Vec<int>& count, int** list,
				SVec<double,3>& X, V6NodeData& v6data, bool adaptive)
#endif
{

  Vec3D d(X[i]);
  Vec3D e(X[j]);

/*
  double unit = sqrt ((d[0]-e[0])*(d[0]-e[0]) + (d[1]-e[1])*(d[1]-e[1]) + (d[2]-e[2])*(d[2]-e[2]));
  double length;
  if (unit == 0.0){
    fprintf(stderr, "*** Error: normal is equal to null vector\n");
    exit(1);
  }
*/

#ifdef EDGE_LENGTH //HB
  if(refLength) {
    Vec3D ed = d-e;
    d = e + (100.0*(*refLength/ed.norm()))*ed;
  }
  else { d = e + 1000.0*(d-e); }
#else
  d = e + 1000.0*(d-e);
#endif

  for (int k=0; k<count[i]; ++k) {
    int tt = list[i][k];
    Elem &elem = elems[tt];

/*
    length = elem.computeLongestEdge(X);
    d = d +5.0 * length/unit * (d-e);
*/

    if ((elem[0] != j) && (elem[1] != j) &&
	(elem[2] != j) && (elem[3] != j)) {
      // get the tetra face (should always work)
      int tf = -1;
      int kk;
      for (kk=0; kk<4; ++kk) {
	if (elem[ elem.faceDef(kk,0) ] != i &&
	    elem[ elem.faceDef(kk,1) ] != i &&
	    elem[ elem.faceDef(kk,2) ] != i) {
	  tf = kk;
	  break;
	}
      }
      if (tf == -1) {
	fprintf(stderr, "*** Error: face not found\n");
	exit(1);
      }

      // perform the intersection
      Vec3D a(X[ elem[ elem.faceDef(tf,0) ] ]);
      Vec3D b(X[ elem[ elem.faceDef(tf,1) ] ]);
      Vec3D c(X[ elem[ elem.faceDef(tf,2) ] ]);
      int pierce = checkPiercePoint(a, b, c, d, e, adaptive);
      if (pierce >= -3 && pierce <= 3) {
        double r, t;
        constructPiercePoint(a, b, c, d, e, 0, &r, &t, true);
        v6data.tet = tt;
        v6data.face = tf;
        v6data.r = r;
        v6data.t = t;
        return true;
      }
    }
  }

  return false;

}

//------------------------------------------------------------------------------

void SubDomain::findEdgeTetrahedra(SVec<double,3>& X, V6NodeData (*&v6data)[2])
{

  const int numNodes = nodes.size();
  const int numElems  = elems.size();

  // create new node type (including interfaces)
  Vec<bool> inode(numNodes);
  int i, j, l;
  for (i=0; i<numNodes; ++i) {
    inode[i] = (nodeType[i]==0) ? true : false;
  }
  for (int iSub=0; iSub<numNeighb; ++iSub)
    for (i=0; i<sharedNodes->num(iSub); ++i)
      inode[ (*sharedNodes)[iSub][i] ] = false;

  // create nodes -> tetrahedra data structure
  Vec<int> count(numNodes);
  count = 0;
  for (i=0; i<numElems; ++i)
    for (j=0; j<4; ++j)
      count[ elems[i][j] ]++;

  int** list = new int*[numNodes];
  for (i=0; i<numNodes; ++i) {
    list[i] = new int[ count[i] ];
    count[i]= 0;
  }
  for (i=0; i<numElems; ++i)
    for (j=0; j<4; ++j)
      list[  elems[i][j] ][ count[ elems[i][j] ]++ ] = i;

  // create edges -> tetrahedra data structure
  const int numEdges = edges.size();
  if (!v6data)
    v6data = new V6NodeData[numEdges][2];
  int (*ptrEdge)[2] = edges.getPtr();

#ifdef EDGE_LENGTH //HB
  // buid array maxEdgeLength: give for each node the largest length of the edges connnected to this node
  edges.updateLength(X); // compute edges'length
  Vec<double> maxEdgeLength(numNodes);
  maxEdgeLength = 0.0;
  for (l=0; l<numEdges; ++l) {
    double elength = edges.length(l);
    if(maxEdgeLength[ptrEdge[l][0]]<elength) maxEdgeLength[ptrEdge[l][0]] = elength;
    if(maxEdgeLength[ptrEdge[l][1]]<elength) maxEdgeLength[ptrEdge[l][1]] = elength;
  }
#endif
  for (l=0; l<numEdges; ++l) {
    for (int dir=0; dir<2; ++dir) {
      if (dir == 0) {
	i = ptrEdge[l][0];
	j = ptrEdge[l][1];
      }
      else {
	i = ptrEdge[l][1];
	j = ptrEdge[l][0];
      }
      double refLength = -1.0;
#ifdef EDGE_LENGTH
      refLength = maxEdgeLength[i];
      if(refLength==0.0) fprintf(stderr, " *** Warning: no refLength = 0.0 at node %6d of subd %3d\n",
                                 i,globSubNum);
      bool ans = findTetrahedron(i, j, count, list, X, v6data[l][dir], true, &refLength);
#else
      bool ans = findTetrahedron(i, j, count, list, X, v6data[l][dir]);
#endif
      if (inode[i] && !ans)
#ifdef EDGE_LENGTH
	ans = findTetrahedron(i, j, count, list, X, v6data[l][dir], false, &refLength);
#else
	ans = findTetrahedron(i, j, count, list, X, v6data[l][dir], false);
#endif
      if (inode[i] && !ans)
#ifdef EDGE_LENGTH
	fprintf(stderr, "*** Warning: no tetrahedron found for edge %d -> %d, ref length = %e\n",
		locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1,refLength);
#else
	fprintf(stderr, "*** Warning: no tetrahedron found for edge %d -> %d\n",
		locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1);
#endif
    }
  }

  /*
  for (l=0; l<numEdges; ++l) {
    i = ptrEdge[l][0];
    j = ptrEdge[l][1];
    if (l % 1000 == 0) {
      char name[500];
      sprintf(name, "edge.%d", l);
      FILE* fp = fopen(name, "w");
      if (v6data[l][0].tet != -1) {
	fprintf(fp, "Elements t0.%d using FluidNodes\n", l);
	fprintf(fp, "1 5 %d %d %d %d\n",
		locToGlobNodeMap[elems[v6data[l][0].tet][0]]+1,
		locToGlobNodeMap[elems[v6data[l][0].tet][1]]+1,
		locToGlobNodeMap[elems[v6data[l][0].tet][2]]+1,
		locToGlobNodeMap[elems[v6data[l][0].tet][3]]+1);
      }
      if (v6data[l][1].tet != -1) {
	fprintf(fp, "Elements t1.%d using FluidNodes\n", l);
	fprintf(fp, "1 5 %d %d %d %d\n",
		locToGlobNodeMap[elems[v6data[l][1].tet][0]]+1,
		locToGlobNodeMap[elems[v6data[l][1].tet][1]]+1,
		locToGlobNodeMap[elems[v6data[l][1].tet][2]]+1,
		locToGlobNodeMap[elems[v6data[l][1].tet][3]]+1);
      }
      fprintf(fp, "Elements t.%d using FluidNodes\n", l);
      int k;
      for (k=0; k<count[i]; ++k)
	fprintf(fp, "%d 5 %d %d %d %d\n", k+1,
		locToGlobNodeMap[elems[list[i][k]][0]]+1,
		locToGlobNodeMap[elems[list[i][k]][1]]+1,
		locToGlobNodeMap[elems[list[i][k]][2]]+1,
		locToGlobNodeMap[elems[list[i][k]][3]]+1);
      for (k=0; k<count[j]; ++k)
	fprintf(fp, "%d 5 %d %d %d %d\n", k+count[i]+1,
		locToGlobNodeMap[elems[list[j][k]][0]]+1,
		locToGlobNodeMap[elems[list[j][k]][1]]+1,
		locToGlobNodeMap[elems[list[j][k]][2]]+1,
		locToGlobNodeMap[elems[list[j][k]][3]]+1);
      fprintf(fp, "Elements e.%d using FluidNodes\n", l);
      fprintf(fp, "1 1 %d %d\n", locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1);
      fclose(fp);
    }
  }
  exit(1);
  */

  if (list) {
    for (i=0; i<numNodes; i++)
      if (list[i]) delete [] list[i];
    delete [] list;
  }

}

//------------------------------------------------------------------------------
bool SubDomain::findNormalTet1(Vec3D a, Vec3D b, Vec3D c, Vec3D d, Vec3D e, int tt, int tf, ExtrapolationNodeData &data, bool adaptive)
{
  int pierce = checkPiercePoint(a, b, c, d, e, adaptive);
  if (pierce >= -3 && pierce <= 3) {
    double r, t;
    constructPiercePoint(a, b, c, d, e, 0, &r, &t, true);
    data.tet = tt;
    data.face = tf;
    data.r = r;
    data.t = t;
    return true;
  }

  return false;
}
//------------------------------------------------------------------------------

inline
bool find1in3(int node, int* itet)
{

  return (node == itet[0] || node == itet[1] || node == itet[2]);

}

//------------------------------------------------------------------------------
void SubDomain::findNormalTet2( int* itet, Vec3D d, Vec3D e, int tt,
                                SVec<double,3> X,
                                ExtrapolationNodeData &data)
{

  int i, j, k, l;
  int nodesFace[3];

  data.tet = tt;

  for ( i=0; i<4; i++){
    Elem &elem = elems[tt];
    nodesFace[0] = elem[ elem.faceDef(i,0) ];
    nodesFace[1] = elem[ elem.faceDef(i,1) ];
    nodesFace[2] = elem[ elem.faceDef(i,2) ];
    if( !find1in3(itet[0], nodesFace) || !find1in3(itet[1], nodesFace) || !find1in3(itet[2], nodesFace) ){
      Vec3D a(X[nodesFace[0]]);
      Vec3D b(X[nodesFace[1]]);
      Vec3D c(X[nodesFace[2]]);

      double s,r,t;
      int pierce = checkPiercePoint(a, b, c, d, e, false);
      if (pierce >= -3 && pierce <= 3) {

        constructPiercePoint(a, b, c, d, e, 0, &r, &t, false);
	data.face = i;
	data.r = r;
	data.t = t;
	return;
      }
    }
  }
}

//------------------------------------------------------------------------------

bool SubDomain::findNormalTetrahedron(int node, Vec3D normal,
                                        int *list, int *list2, int num,
                                        SVec<double,3> X,
                                        ExtrapolationNodeData *xpoldata,
                                        bool adaptive)
{
  int nn = locToGlobNodeMap[node]+1;
  int itemp[3];
  bool ans;
  int itet[3];
  Vec3D dprime(X[node]);
  Vec3D d;
  Vec3D e(X[node][0]+normal[0], X[node][1]+normal[1], X[node][2]+normal[2]);
  double unit = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
  double length;
  if (unit == 0.0){
    fprintf(stderr, "*** Error: normal is equal to null vector\n");
    exit(1);
  }

  for (int k=0; k<num; ++k){
    int tt = list[k];
    Elem &elem = elems[tt];

    length = elem.computeLongestEdge(X);
    d = dprime + 100.0*length/unit * ( dprime-e );
    int tf = -1;
    for ( int kk=0; kk<4; kk++){
      if (elem[ elem.faceDef(kk,0) ] != node &&
	  elem[ elem.faceDef(kk,1) ] != node &&
	  elem[ elem.faceDef(kk,2) ] != node) {
	tf = kk;
	break;
      }
    }
    if (tf == -1){
      fprintf(stderr, "***Error: face not found\n");
      exit(1);
    }
    itet[0] = elem[ elem.faceDef(tf,0) ];
    itet[1] = elem[ elem.faceDef(tf,1) ];
    itet[2] = elem[ elem.faceDef(tf,2) ];

    Vec3D a(X[ itet[0]  ]);
    Vec3D b(X[ itet[1]  ]);
    Vec3D c(X[ itet[2]  ]);

    //for Constant EXTRAPOLATION
    ans = findNormalTet1(a, b, c, d, e, tt, tf, xpoldata[0], adaptive);

    //for Linear EXTRAPOLATION
    if (ans && list2){
      int tt2 = list2[k];
      if (tt2 != -1){
	findNormalTet2(itet, d, e, tt2, X, xpoldata[1]);

	Vec3D h1 = c +xpoldata[0].r*(a-c)+xpoldata[0].t*(b-c);

	elem = elems[xpoldata[1].tet];
	itemp[0] = elem[ elem.faceDef(xpoldata[1].face,0) ];
	itemp[1] = elem[ elem.faceDef(xpoldata[1].face,1) ];
	itemp[2] = elem[ elem.faceDef(xpoldata[1].face,2) ];
	Vec3D aa(X[ itemp[0] ]);
	Vec3D bb(X[ itemp[1] ]);
	Vec3D cc(X[ itemp[2] ]);
	Vec3D h2 = cc +xpoldata[1].r*(aa-cc)+xpoldata[1].t*(bb-cc);
	Vec3D d10 = h1-dprime;
	Vec3D d20 = h2-dprime;
	double tet1[3], tet2[3], tet3[3];
	//tet1 = n^X1X0, tet2 = n^X2X0, tet3 = X2x0^X1X0
	tet1[0]=normal[1]*d10[2]-normal[2]*d10[1];
	tet1[1]=normal[2]*d10[0]-normal[0]*d10[2];
	tet1[2]=normal[0]*d10[1]-normal[1]*d10[0];

	tet2[0]=normal[1]*d20[2]-normal[2]*d20[1];
	tet2[1]=normal[2]*d20[0]-normal[0]*d20[2];
	tet2[2]=normal[0]*d20[1]-normal[1]*d20[0];

	tet3[0]=d10[1]*d20[2]-d10[2]*d20[1];
	tet3[1]=d10[2]*d20[0]-d10[0]*d20[2];
	tet3[2]=d10[0]*d20[1]-d10[1]*d20[0];
	double tot1, tot2, tot3;
	tot1 = pow(tet1[0]*tet1[0]+tet1[1]*tet1[1]+tet1[2]*tet1[2], 0.5);
	tot2 = pow(tet2[0]*tet2[0]+tet2[1]*tet2[1]+tet2[2]*tet2[2], 0.5);
	tot3 = pow(tet3[0]*tet3[0]+tet3[1]*tet3[1]+tet3[2]*tet3[2], 0.5);
      }
    }
    if (ans) break;
  }
  return ans;
}

//------------------------------------------------------------------------------

void SubDomain::findNormalTetrahedra(SVec<double,3>& X, Vec<Vec3D>& normals,
				     ExtrapolationNodeData (*&xpoldata)[2])
{
  const int numBcNodes  = inletNodes.size();
  const int numNodes = nodes.size();
  int *listtets  = 0;
  int *listtets2 = 0;
  int numtets, node;

  int i;
  int counter = 0;
  int counter2 = 0;

  //inode[i] is false if node i is shared.
  Vec<bool> inode(numNodes);
  for ( i=0; i<numNodes; i++)
    inode[i] = true;
  for (int iSub=0; iSub<numNeighb; ++iSub)
    for (i=0; i<sharedNodes->num(iSub); ++i){
      inode[ (*sharedNodes)[iSub][i] ] = false;
    }

  if(!xpoldata)
    xpoldata = new ExtrapolationNodeData[numBcNodes][2];

  for ( i=0; i<numBcNodes; i++){
    node = inletNodes[i].getNodeNum();
    listtets = inletNodes[i].getTets();
    listtets2 = inletNodes[i].getTets2();
    numtets = inletNodes[i].getNumTets();
    bool ans = findNormalTetrahedron(node, normals[i], listtets, listtets2,  numtets, X, xpoldata[i]);

    if(!ans){
       ans = findNormalTetrahedron(node, normals[i], listtets, listtets2, numtets, X, xpoldata[i], false);
    }
    if(inode[node] && !ans){
      counter++;
      fprintf(stderr, "***Warning: no tetrahedron found for inletNode %d\n", locToGlobNodeMap[node]+1);
    }
  }

}

//------------------------------------------------------------------------------

void SubDomain::checkNormalTetrahedra(ExtrapolationNodeData (*&xpoldata)[2])
{

  int i, counter=0;
  int counter2=0;

  for (i=0; i<inletNodes.size(); i++){
    if ( xpoldata[i][0].tet==-1 ){
      counter++;
      //fprintf(stderr, "node %d has no corresponding tet\n", i);
    }
  }

  for (int iSub=0; iSub<numNeighb; ++iSub)
    for (i=0; i<sharedNodes->num(iSub); ++i)
      counter2++;
  fprintf(stderr, "# of nodes w/o tets is %d      and # of shared nodes for that subdomain is %d\n", counter, counter2);

}

//------------------------------------------------------------------------------

void SubDomain::setInletNodes(IoData& ioData)
{

  int i, j;

  //what we need to pass in the InletNodeSet class to store
  int numInletNodes = 0;

  //find numInletNodes and tag the nodes that are at inlet/outlet
  Vec<bool> tagNodes(nodes.size());
  tagNodes = false;
  for ( i=0; i<nodes.size(); i++){
    if ( nodeType[i] == BC_OUTLET_MOVING || nodeType[i] == BC_INLET_MOVING ||
         nodeType[i] == BC_OUTLET_FIXED  || nodeType[i] == BC_INLET_FIXED ){
      tagNodes[i] = true;
      numInletNodes++;
    }
  }

  //check that we get the same number of InletNodes counting via nodes, faces, and tets
  //via tets
  Vec<bool> tcounted(nodes.size());
  tcounted = false;
  int tnumInletNodes=0;
  int n;
  for( i=0; i<elems.size(); i++){
    for( j=0; j<4; j++){
      n=elems[i][j];
      if (!tcounted[n] && (nodeType[n] ==BC_OUTLET_MOVING || nodeType[n] ==BC_INLET_MOVING ||
			   nodeType[n] ==BC_OUTLET_FIXED ||nodeType[n] ==BC_INLET_FIXED )){
        tcounted[n] = true;
        tnumInletNodes++;
      }
    }
  }

  //via faces
  Vec<bool> fcounted(nodes.size());
  fcounted = false;
  int fnumInletNodes=0;
  for( i=0; i<faces.size(); i++){
    for( j=0; j<3; j++){
      n=faces[i][j];
      if (!fcounted[n] && (nodeType[n] ==BC_OUTLET_MOVING || nodeType[n] ==BC_INLET_MOVING ||
			   nodeType[n] ==BC_OUTLET_FIXED ||nodeType[n] ==BC_INLET_FIXED )){
	fcounted[n] = true;
	fnumInletNodes++;
      }
    }
  }

  //take care of numTets and tets for each InletNode
  //first numTets of each InletNode
  int *counttets = new int[nodes.size()];
  int **listtets = new int*[nodes.size()];

  for ( i=0; i<nodes.size(); i++)
    counttets[i]=0;

  for ( i=0; i<elems.size(); i++)
    for ( j=0; j<4; j++)
      counttets[ elems[i][j] ]++;


  //now the tets connected to each InletNode
  for ( i=0; i<nodes.size(); i++ )
    if ( tagNodes[i]==true ){
      listtets[i] = new int[counttets[i]];
      counttets[i]=0;
    }

  for( i=0; i<elems.size(); i++ )
    for( j=0; j<4; j++)
      if( tagNodes[ elems[i][j] ] == true )
        listtets[ elems[i][j] ][ counttets[elems[i][j]]++ ] = i;


  //take care of numFaces and faces for each InletNode
  //first numFaces of each InletNode
  int *countfaces = new int[nodes.size()];
  int **listfaces = new int*[nodes.size()];

  for ( i=0; i<nodes.size(); i++)
    countfaces[i]=0;

  for ( i=0; i<faces.size(); i++)
    if (faces[i].getCode() == BC_INLET_FIXED  || faces[i].getCode() == BC_INLET_MOVING ||
        faces[i].getCode() == BC_OUTLET_FIXED || faces[i].getCode() == BC_OUTLET_MOVING )
      for ( j=0; j<3; j++)
        countfaces[ faces[i][j] ]++;

  //now the faces connected to each InletNode
  for ( i=0; i<nodes.size(); i++ )
    if ( tagNodes[i] == true ){
      listfaces[i] = new int[countfaces[i]];
      countfaces[i]=0;
    }

  for( i=0; i<faces.size(); i++ )
    for( j=0; j<3; j++)
      if (tagNodes[ faces[i][j] ] == true &&
	  (faces[i].getCode() == BC_INLET_FIXED  || faces[i].getCode() == BC_INLET_MOVING ||
	   faces[i].getCode() == BC_OUTLET_FIXED || faces[i].getCode() == BC_OUTLET_MOVING ))
        listfaces[ faces[i][j] ][ countfaces[faces[i][j]]++ ] = i;

  //and create all the InletNode we need
  inletNodes.setup(locSubNum, numInletNodes, ioData);

  //finally get only the inlet nodes in the storage facility inletNodes
  int counter = 0;
  for ( i=0; i<nodes.size(); i++){
    if ( tagNodes[i] == true ){
      inletNodes[counter].setInletNode(i, countfaces[i], counttets[i], listfaces[i], listtets[i]);
      counter++;
    }
  }

  //if(countfaces) delete  countfaces;
  //if(counttets) delete  counttets;
  //if(listfaces)
  //      delete [] listfaces;
  //if(listtets)
  //      delete [] listtets;

}

//------------------------------------------------------------------------------

void SubDomain::setInletNodes2(IoData& ioData)
{

  int i, j, k, l;
  int node, count, tet;
  int inum[3];
  int *list;
  int numInletNodes = inletNodes.size();


  int** listtets = new int*[numInletNodes];

  for ( i=0; i<numInletNodes; i++){
    node  = inletNodes[i].getNodeNum();
    count = inletNodes[i].getNumTets();
    listtets[i] = new int[count];
    list  = inletNodes[i].getTets();
    for ( j=0; j<count; j++){
      tet = list[j];
      l = 0;
      for ( k=0; k<4; k++)
        if ( elems[tet][k] != node)
          inum[l++] = elems[tet][k];
      listtets[i][j] = findOppositeTet(inum, node);

    }
    inletNodes[i].addSecondTets(listtets[i]);
  }
}

//------------------------------------------------------------------------------

int SubDomain::findOppositeTet(int *itet, int node)
{

  int tet;
  int numTets = elems.size();
  bool match1 = false;
  bool match2 = false;
  bool match3 = false;
  bool match4 = false;
  bool total  = true;
  tet = 0;

  while(total){
    match1 = findNodeInTet(itet[0], tet);
    if (match1){
      match2 = findNodeInTet(itet[1], tet);
      if (match2){
        match3 = findNodeInTet(itet[2], tet);
        if (match1 && match2 && match3) match4 = !(findNodeInTet(node, tet));
      }
    }
    total = !(match1 && match2 && match3 && match4);
    if(total && (tet==numTets-1)){
      return -1;
    }
    tet++;
  }
  return tet-1;
}

//------------------------------------------------------------------------------

bool SubDomain::findNodeInTet(int node, int tet)
{
  return ( (node==elems[tet][0]) || (node==elems[tet][1]) ||
	   (node==elems[tet][2]) || (node==elems[tet][3]) );
}

//------------------------------------------------------------------------------
void SubDomain::checkInletNodes()
{
  inletNodes.checkInletNodes(locToGlobNodeMap);
}

//------------------------------------------------------------------------------

void SubDomain::sumInletNormals(Vec<Vec3D>& inletNodeNorm, Vec<Vec3D>& faceNorm, Vec<int>& numFaceNeighb)
{
  int numNodes = inletNodes.size();
  int i,j,m;
  int *listFaces;

  inletNodeNorm = 0.0;
  numFaceNeighb = 0;

  for ( i=0; i<numNodes; i++){
    m = inletNodes[i].getNodeNum();
    numFaceNeighb[i] = inletNodes[i].getNumFaces();
    listFaces = inletNodes[i].getFaces();
    for ( j=0; j<numFaceNeighb[i]; j++)
      inletNodeNorm[i] += faces[listFaces[j]].getNormal(faceNorm);
  }
}

//------------------------------------------------------------------------------

void SubDomain::numDivideNormals(Vec<Vec3D>& inletNodeNorm, Vec<int>& numFaceNeighb)
{
  for( int i=0; i<inletNodes.size(); i++)
    if (numFaceNeighb[i]==0) {
      fprintf(stderr, "division by zero!\n");
      exit(1);
    } else
      inletNodeNorm[i] /= numFaceNeighb[i];
}

//------------------------------------------------------------------------------

void SubDomain::getReferenceMeshPosition(SVec<double,3> &X)
{

  X = nodes;

}

//------------------------------------------------------------------------------

void SubDomain::computeDisplacement(SVec<double,3> &X, SVec<double,3> &dX)
{

  dX = X - nodes;

}

//------------------------------------------------------------------------------

void SubDomain::computeDisplacement(SVec<double,3> &X, double* dX,int node)
{

  for (int k = 0; k < 3; ++k)
    dX[k] = X[node][k] - nodes[node][k];

}


//------------------------------------------------------------------------------


void SubDomain::markLenKirchhoffNodes(IoData &iod, DistInfo &distInfo)
{

// This function counts the nodes on a Kirchhoff surface.
// The nodes are counted per subdomain.
// Some nodes may appear on several subdomains.
//
// The function sets the number of nodes in the object 'distInfo'.
//
// UH (07/2012)

  if (strlen(iod.output.restart.strKPtraces) == 0)
  {
    distInfo.setLen(locSubNum, 0);
    return;
  }

  for (int j = 0; j < faces.size(); ++j)
  {
    if (faces[j].getCode() == BC_KIRCHHOFF_SURFACE)
    {
      for (int k = 0; k < faces[j].numNodes(); ++k)
        kirchhoffNodesList.insert(faces[j][k]);
    }
  }

  distInfo.setLen(locSubNum, kirchhoffNodesList.size());

}

//------------------------------------------------------------------------------

void SubDomain::makeMasterFlag(DistInfo &distInfo)
{

  bool *flag    = distInfo.getMasterFlag(locSubNum);
  double *weight= distInfo.getInvWeight(locSubNum);

  if (!flag) return;

  int numNodes = nodes.size();

  int i;
  for (i = 0; i < numNodes; ++i)
    flag[i] = true;

  if (weight)
    for (i = 0; i < numNodes; ++i)
      weight[i] = 1.0;

  int iSub;
  for (iSub = 0; iSub < numNeighb; ++iSub) {

    if (weight)
      for (i = 0; i < sharedNodes->num(iSub); ++i)
        weight[ (*sharedNodes)[iSub][i] ]++;

    if (neighb[iSub] > globSubNum) continue;

    for (i = 0; i < sharedNodes->num(iSub); ++i)
      flag[ (*sharedNodes)[iSub][i] ] = false;

  }

  if (weight)
    for (i = 0; i < numNodes; ++i)
      weight[i] = sqrt(1.0 / weight[i]);

}

//------------------------------------------------------------------------------

void SubDomain::makeEdgeMasterFlag(DistInfo &distInfo, CommPattern<int> &edgeNumPat)
{

  int iSub, iEdge;

  // We receive the neighbor's list of edges and then also build the edgeMasterFlag

  bool *edgeMasterFlag = distInfo.getMasterFlag(locSubNum); 
  if(!edgeMasterFlag) return;

  for (iEdge = 0; iEdge < edges.size(); ++iEdge) edgeMasterFlag[iEdge] = true;

  for (iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<int> sInfo = edgeNumPat.recData(rcvChannel[iSub]);

    int nIndex = 0;
    int myIndex = 0;

    for (iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge) {

      int glLeft  = sharedEdges[iSub][iEdge].glLeft;
      int glRight = sharedEdges[iSub][iEdge].glRight;

      while (2*nIndex < sInfo.len && (sInfo.data[2*nIndex] < glLeft || (sInfo.data[2*nIndex] == glLeft && sInfo.data[2*nIndex+1] < glRight)))
        nIndex++;

      if (2*nIndex < sInfo.len &&
	      (sInfo.data[2*nIndex] == glLeft && sInfo.data[2*nIndex+1] == glRight)) 
      {
        sharedEdges[iSub][myIndex] = sharedEdges[iSub][iEdge];
        myIndex++;
      }

    }

    numSharedEdges[iSub] = myIndex;

    if (neighb[iSub] < globSubNum) // I cannot be the master of these edges
      for (iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge)
        edgeMasterFlag[ sharedEdges[iSub][iEdge].edgeNum ] = false;

  }

}

//------------------------------------------------------------------------------

void SubDomain::setChannelNums(SubDTopo &subTopo)
{

  rcvChannel = new int[numNeighb];
  sndChannel = new int[numNeighb];

  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    rcvChannel[iSub] = subTopo.getChannelID(neighb[iSub], globSubNum);
    sndChannel[iSub] = subTopo.getChannelID(globSubNum, neighb[iSub]);
  }

}

//------------------------------------------------------------------------------

void SubDomain::identifyEdges(CommPattern<int> &edgeNumPat)
{

  int iSub, jSub;

  Connectivity *nodeToNeighb = sharedNodes->reverse();

  numSharedEdges = new int[numNeighb];

  for (iSub = 0; iSub < numNeighb; ++iSub) numSharedEdges[iSub] = 0;

  int (*edgePtr)[2] = edges.getPtr();

  int l;
  for (l=0; l<edges.size(); ++l) {

    int leftNode = edgePtr[l][0];
    int rightNode = edgePtr[l][1];

    if (nodeToNeighb->num(leftNode) == 0 ||
	nodeToNeighb->num(rightNode) == 0) continue;

    for (iSub = 0; iSub < nodeToNeighb->num(leftNode); ++iSub) {

      int subI = (*nodeToNeighb)[leftNode][iSub];

      for (jSub = 0; jSub < nodeToNeighb->num(rightNode); ++jSub)
	if (subI == (*nodeToNeighb)[rightNode][jSub]) numSharedEdges[subI] += 1;

    }

  }

  sharedEdges = new EdgeDef *[numNeighb];

  for (iSub = 0; iSub < numNeighb; ++iSub)
    sharedEdges[iSub] = new EdgeDef[ numSharedEdges[iSub] ];

  for (iSub = 0; iSub < numNeighb; ++iSub) numSharedEdges[iSub] = 0;

  for (l=0; l<edges.size(); ++l) {

    int leftNode = edgePtr[l][0];
    int rightNode = edgePtr[l][1];

    if (nodeToNeighb->num(leftNode) == 0 ||
       nodeToNeighb->num(rightNode) == 0) continue;

    for (iSub = 0; iSub < nodeToNeighb->num(leftNode); ++iSub) {

      int subI = (*nodeToNeighb)[leftNode][iSub];

      for (jSub = 0; jSub < nodeToNeighb->num(rightNode); ++jSub) {

	if (subI == (*nodeToNeighb)[rightNode][jSub]) {

	  sharedEdges[subI][ numSharedEdges[subI] ].edgeNum = l;
	  sharedEdges[subI][ numSharedEdges[subI] ].glLeft = locToGlobNodeMap[leftNode];
	  sharedEdges[subI][ numSharedEdges[subI] ].glRight = locToGlobNodeMap[rightNode];
	  sharedEdges[subI][ numSharedEdges[subI] ].order();

	  numSharedEdges[subI] += 1;

	}

      }

    }
  }

  for (iSub = 0; iSub < numNeighb; ++iSub) {
#ifdef OLD_STL
    sort(sharedEdges[iSub], sharedEdges[iSub]+numSharedEdges[iSub]);
#else
    stable_sort(sharedEdges[iSub], sharedEdges[iSub]+numSharedEdges[iSub]);
#endif
    edgeNumPat.setLen(sndChannel[iSub], 2*numSharedEdges[iSub]);
  }

  delete nodeToNeighb;

}

//------------------------------------------------------------------------------

void SubDomain::sndEdgeInfo(CommPattern<int> &edgeNumPat)
{

  int iSub, iEdge;

  for (iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<int> sInfo = edgeNumPat.getSendBuffer(sndChannel[iSub]);

    for (iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge) {
      sInfo.data[2*iEdge]   = sharedEdges[iSub][iEdge].glLeft;
      sInfo.data[2*iEdge+1] = sharedEdges[iSub][iEdge].glRight;
    }

  }

}

//------------------------------------------------------------------------------

void SubDomain::rcvEdgeInfo(CommPattern<int> &edgeNumPat)
{

  int iSub, iEdge;

  // We receive the neighbor's list of edges and then also build the edgeMasterFlag

  bool *edgeMasterFlag = new bool[edges.size()];

  for (iEdge = 0; iEdge < edges.size(); ++iEdge) edgeMasterFlag[iEdge] = true;

  for (iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<int> sInfo = edgeNumPat.recData(rcvChannel[iSub]);

    int nIndex = 0;
    int myIndex = 0;

    for (iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge) {

      int glLeft  = sharedEdges[iSub][iEdge].glLeft;
      int glRight = sharedEdges[iSub][iEdge].glRight;

      while (2*nIndex < sInfo.len && (sInfo.data[2*nIndex] < glLeft ||
	      (sInfo.data[2*nIndex] == glLeft && sInfo.data[2*nIndex+1] < glRight)))
        nIndex++;

      if (2*nIndex < sInfo.len &&
	      (sInfo.data[2*nIndex] == glLeft && sInfo.data[2*nIndex+1] == glRight)) 
      {
        sharedEdges[iSub][myIndex] = sharedEdges[iSub][iEdge];
        myIndex++;
      }

    }

    numSharedEdges[iSub] = myIndex;

    if (neighb[iSub] < globSubNum) // I cannot be the master of these edges
      for (iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge)
        edgeMasterFlag[ sharedEdges[iSub][iEdge].edgeNum ] = false;

  }

  edges.setMasterFlag(edgeMasterFlag);
}

//------------------------------------------------------------------------------

void SubDomain::sndNormals(CommPattern<double> &edgePat, Vec3D *edgeNorm,
			   double *edgeNormVel)
{

  int iSub, iEdge;

  for (iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<double> sInfo = edgePat.getSendBuffer(sndChannel[iSub]);

    for (iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge) {
      int edgeNum = sharedEdges[iSub][iEdge].edgeNum;
      if (sharedEdges[iSub][iEdge].sign > 0) {
	sInfo.data[4*iEdge]   = edgeNorm[edgeNum][0];
	sInfo.data[4*iEdge+1] = edgeNorm[edgeNum][1];
	sInfo.data[4*iEdge+2] = edgeNorm[edgeNum][2];
	sInfo.data[4*iEdge+3] = edgeNormVel[edgeNum];
      }
      else {
	sInfo.data[4*iEdge]   = -edgeNorm[edgeNum][0];
	sInfo.data[4*iEdge+1] = -edgeNorm[edgeNum][1];
	sInfo.data[4*iEdge+2] = -edgeNorm[edgeNum][2];
	sInfo.data[4*iEdge+3] = -edgeNormVel[edgeNum];
      }
    }
  }

}

//------------------------------------------------------------------------------

void SubDomain::sndEdgeData(CommPattern<double> &edgePat, double *edgeData)
{
  int iSub, iEdge;
  for (iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<double> sInfo = edgePat.getSendBuffer(sndChannel[iSub]);
    for (iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge) {
      int edgeNum = sharedEdges[iSub][iEdge].edgeNum;
      sInfo.data[iEdge]   = edgeData[edgeNum];
    }
  }
}



//------------------------------------------------------------------------------

void SubDomain::rcvNormals(CommPattern<double> &edgePat, Vec3D *edgeNorm,
			   double *edgeNormVel)
{

  int iSub, iEdge;

  for (iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<double> sInfo = edgePat.recData(rcvChannel[iSub]);

    for (iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge) {

      int edgeNum = sharedEdges[iSub][iEdge].edgeNum;

      if (sharedEdges[iSub][iEdge].sign > 0) {
	edgeNorm[edgeNum][0] += sInfo.data[4*iEdge];
	edgeNorm[edgeNum][1] += sInfo.data[4*iEdge+1];
	edgeNorm[edgeNum][2] += sInfo.data[4*iEdge+2];
	edgeNormVel[edgeNum] += sInfo.data[4*iEdge+3];
      }
      else {
	edgeNorm[edgeNum][0] -= sInfo.data[4*iEdge];
	edgeNorm[edgeNum][1] -= sInfo.data[4*iEdge+1];
	edgeNorm[edgeNum][2] -= sInfo.data[4*iEdge+2];
	edgeNormVel[edgeNum] -= sInfo.data[4*iEdge+3];
      }

    }

  }

}

//------------------------------------------------------------------------------

void SubDomain::rcvEdgeData(CommPattern<double> &edgePat, double *edgeData)
{
  int iSub, iEdge;
  for (iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<double> sInfo = edgePat.recData(rcvChannel[iSub]);
    for (iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge) {
      int edgeNum = sharedEdges[iSub][iEdge].edgeNum;
      edgeData[edgeNum] += sInfo.data[iEdge];
    }
  }
}

//------------------------------------------------------------------------------

void SubDomain::setFaceType(int *facemap)
{

  for (int i=0; i<faces.size(); ++i)
  {
    faces[i].setType(facemap);
  }

}

//------------------------------------------------------------------------------

void SubDomain::setNodeType(int* priority, CommPattern<int> &ntP)
{

  nodeType = new int[ nodes.size() ];

  int i;
  for (i = 0; i < nodes.size(); ++i)
    nodeType[i] = BC_INTERNAL;

  for (i = 0; i < faces.size(); ++i)
    faces[i].setNodeType(priority, nodeType);

  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<int> nInfo = ntP.getSendBuffer(sndChannel[iSub]);
    for (i = 0; i < sharedNodes->num(iSub); ++i)
      nInfo.data[i] = nodeType[ (*sharedNodes)[iSub][i] ];
  }

}

//------------------------------------------------------------------------------

void SubDomain::setNodeFaceType(CommPattern<int> &ntP)
{
/* ARL: creates a new array telling if a node is
 *      connected to faces that are
 *      inlet only (1)
 *      inlet and wall (2)
 *      wall only (-1)
 *      if it belongs to no face, it is then set to 0
 */

  nodeFaceType = new int[ nodes.size() ];

  int i;
  for (i = 0; i < nodes.size(); ++i)
    nodeFaceType[i] = 0;

  for (i = 0; i < faces.size(); ++i)
    faces[i].setNodeFaceType(nodeFaceType);

  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<int> nInfo = ntP.getSendBuffer(sndChannel[iSub]);
    for (i = 0; i < sharedNodes->num(iSub); ++i)
      nInfo.data[i] = nodeFaceType[ (*sharedNodes)[iSub][i] ];
  }

}

//------------------------------------------------------------------------------

int* SubDomain::completeNodeType(int* priority, CommPattern<int> &ntP)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<int> nInfo = ntP.recData(rcvChannel[iSub]);
    for (int i = 0; i < sharedNodes->num(iSub); ++i)
      if (priority[ nInfo.data[i] ] > priority[ nodeType[ (*sharedNodes)[iSub][i] ] ])
        nodeType[ (*sharedNodes)[iSub][i] ] = nInfo.data[i];
  }

  return nodeType;

}

// -----------------------------------------------------------

int* SubDomain::completeNodeFaceType(CommPattern<int> &ntP)
{
  int type;
  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<int> nInfo = ntP.recData(rcvChannel[iSub]);
    for (int i = 0; i < sharedNodes->num(iSub); ++i){
      type = nodeFaceType[ (*sharedNodes)[iSub][i] ];
      if (type == 0)
        nodeFaceType[ (*sharedNodes)[iSub][i] ] = nInfo.data[i];
      if (type == 1)
        if (nInfo.data[i] == -1 || nInfo.data[i] == 2)
          nodeFaceType[ (*sharedNodes)[iSub][i] ] = 2;
      if (type == -1)
        if (nInfo.data[i] == 1 || nInfo.data[i] == 2)
          nodeFaceType[ (*sharedNodes)[iSub][i] ] = 2;
    }
  }

  return nodeFaceType;

}

//------------------------------------------------------------------------------
//
void SubDomain::computeOffWallNode(LevelSetStructure *LSS, CommPattern<int> &ntP)
{
  int (*edgePtr)[2] = edges.getPtr();

  if(!offWallNode) offWallNode = new int[nodes.size()];

  for (int i=0; i<nodes.size(); ++i) offWallNode[i] = 0;
  for (int l=0; l<edges.size(); ++l) {
    if(LSS->edgeIntersectsStructure(0.0,l)) { //at interface
      int i = edgePtr[l][0];
      int j = edgePtr[l][1];
      bool iIsActive = LSS->isActive(0.0,i);
      bool jIsActive = LSS->isActive(0.0,j);
      if(iIsActive) offWallNode[i] = 1;
      if(jIsActive) offWallNode[j] = 1;
    }
  }

  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<int> nInfo = ntP.getSendBuffer(sndChannel[iSub]);
    for (int i = 0; i < sharedNodes->num(iSub); ++i)
      nInfo.data[i] = offWallNode[ (*sharedNodes)[iSub][i] ];
  }
}

//------------------------------------------------------------------------------

int* SubDomain::completeOffWallNode(CommPattern<int> &ntP)
{
  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<int> nInfo = ntP.recData(rcvChannel[iSub]);
    for (int i = 0; i < sharedNodes->num(iSub); ++i)
      if ( nInfo.data[i] > offWallNode[ (*sharedNodes)[iSub][i] ] )
        offWallNode[ (*sharedNodes)[iSub][i] ] = nInfo.data[i];
  }

  return offWallNode;
}

//------------------------------------------------------------------------------
// -----------------------------------------------------------
// HB: create the dofType array using the matchNodeSet, the sliding faces & the nodeType array
// Note that the order in which the dofType is filled is crucial: it is first set to BC_FREE (i.e. all
// the dofs are assumed to be free to move), and then they are (potentially) constrained if necessary.

int*
SubDomain::getMeshMotionDofType(map<int,SurfaceData*>& surfaceMap, CommPattern<int> &ntP, MatchNodeSet* matchNodes)
{
  int* DofType = new int[3*nodes.size()];
  int (*dofType)[3] = reinterpret_cast<int (*)[3]>(DofType);

  // Step 1. Set all the dofs as BC_FREE
  int nMoving = 0;
  for (int i=0; i<nodes.size(); i++) {
    if(nodeType[i] < BC_INTERNAL)
      nMoving++; // this node is labelled as moving
    for(int l = 0; l < 3; l++)
      dofType[i][l] = BC_FREE;
  }

  // Step 2. Take into account the sliding plane constraints
  for (int i=0; i<faces.size(); i++) { // Loop over faces
    bool isSliding = false;
    map<int,SurfaceData*>::iterator it = surfaceMap.find(faces[i].getSurfaceID());
    if(it!=surfaceMap.end()) { // surface has attribute in the input file
      if(it->second->nx != 0.0 || it->second->ny != 0.0 || it->second->nz != 0.0) // it's a sliding surface
        isSliding = true;
  
    }
    if(!isSliding) { // -> constrain the nodes according to face "fluid code"
      switch(faces[i].getCode()) {
        case(BC_SYMMETRY): // by default a symmetry plane is fixed ...
        case(BC_ISOTHERMAL_WALL_FIXED):
        case(BC_ADIABATIC_WALL_FIXED):
        case(BC_POROUS_WALL_FIXED):
        case(BC_OUTLET_FIXED):
        case(BC_INLET_FIXED):
        case(BC_DIRECTSTATE_OUTLET_FIXED):
        case(BC_DIRECTSTATE_INLET_FIXED):
        case(BC_MASSFLOW_OUTLET_FIXED):
        case(BC_MASSFLOW_INLET_FIXED):
        case(BC_SLIP_WALL_FIXED):
        for(int j=0; j<faces[i].numNodes();j++)
          for(int l=0; l<3; l++) dofType[faces[i][j]][l] = BC_FIXED;
        break;
      }
    }
    else  // Sliding is there
    {
      if (it->second->nx != 0.0 ) // it's a sliding surface
        for(int j=0; j<faces[i].numNodes();j++)	 dofType[faces[i][j]][0] = BC_FIXED;

      if (it->second->ny != 0.0) // it's a sliding surface
        for(int j=0; j<faces[i].numNodes();j++)  dofType[faces[i][j]][1] = BC_FIXED;

      if (it->second->nz != 0.0) // it's a sliding surface
        for(int j=0; j<faces[i].numNodes();j++)  dofType[faces[i][j]][2] = BC_FIXED;
    }
  }


  // Step 3. Take into account the matched nodes
  Vec<bool> isMatched(nodes.size());
  isMatched = false;
  if(matchNodes) {
    for(int i=0;i<matchNodes->size();i++) {
      int inode = matchNodes->subMatchNode(i); // current matched node
      isMatched[inode] = true;
      for(int l=0; l<3; l++) dofType[inode][l] = BC_MATCHED;
    }
  } else { // Only for ForcedMeshMotion (see ForcedMeshMotionHandler) where matchNodes is not
           // explicitly created -> consider all the node labelled as moving as "matched" nodes
    for(int i=0;i<nodes.size(); i++)
      if(nodeType[i]<BC_INTERNAL) {  // this node is labelled as moving
        isMatched[i] = true;
        for(int l=0; l<3; l++) dofType[i][l] = BC_MATCHED;
      }
  }

#ifdef HB_MESHMOTION_DEBUG
  // for debugging: count number of nodes that were labeled as moving
  for(int i=0;i<nodes.size(); i++)
    if(nodeType[i]<BC_INTERNAL) { // this node is labelled as moving
      if(!isMatched[i] & (dofType[i][0]!=BC_FREE) & (dofType[i][1]!=BC_FREE) & (dofType[i][2]!=BC_FREE)) {
         fprintf(stderr," ... PROBLEM in SubDomain::getMeshMotionDofType: sub %d, non-matched moving node %d (%d) is not free ...\n",globSubNum,i,locToGlobNodeMap[i]); fflush(stderr);
      }
    }
#endif

  // Step 4. Take into account the mesh motion Dirichlet BCs from the input file (if any)
  if(mmsBCs)
    for(int i=0;i<mmsBCs->size(); i++)
      dofType[(*mmsBCs)[i].nnum][(*mmsBCs)[i].dofnum] = BC_CONSTRAINED;

  // Step 6. Fill communication buffer (to ensure same constraint on the shared nodes/dofs)
  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<int> nInfo = ntP.getSendBuffer(sndChannel[iSub]);
    int (*buffer)[3] = reinterpret_cast<int (*)[3]>(nInfo.data);
    for (int i = 0; i < sharedNodes->num(iSub); ++i)
      for(int l=0; l<3; l++)
        buffer[i][l] = dofType[(*sharedNodes)[iSub][i]][l];
  }

  return(DofType);
}

// -----------------------------------------------------------

int*
SubDomain::getEmbeddedALEMeshMotionDofType(map<int,SurfaceData*>& surfaceMap, CommPattern<int> &ntP, MatchNodeSet* matchNodes,
                                           DefoMeshMotionData::SlidingSurfaceTreatment slidingSurfaceTreatment)
{
  int* DofType = new int[3*nodes.size()];
  int (*dofType)[3] = reinterpret_cast<int (*)[3]>(DofType);

  // Step 1. Set all the dofs as BC_FREE
  int nMoving = 0;
  for (int i=0; i<nodes.size(); i++) {
    if(nodeType[i] < BC_INTERNAL)
      nMoving++; // this node is labelled as moving
    for(int l = 0; l < 3; l++)
      dofType[i][l] = BC_FREE;
  }

  // Step 2. Set appropriate fixed bc's 
  for (int i=0; i<faces.size(); i++) { // Loop over faces
    bool isSliding = false;
    map<int,SurfaceData*>::iterator it = surfaceMap.find(faces[i].getSurfaceID());
    if(it!=surfaceMap.end()) { // surface has attribute in the input file
      if(it->second->nx != 0.0 || it->second->ny != 0.0 || it->second->nz != 0.0) // it's a sliding surface
        isSliding = true;
    }

    if(slidingSurfaceTreatment == DefoMeshMotionData::PrescribedAverage) {
      switch(faces[i].getCode()) {
        case(BC_SYMMETRY):
          if(!isSliding)
            for(int j=0; j<faces[i].numNodes();j++)
              for(int l=0; l<3; l++) dofType[faces[i][j]][l] = BC_FIXED;
          else
            for(int j=0; j<faces[i].numNodes();j++)
              for(int l=0; l<3; l++) dofType[faces[i][j]][l] = (dofType[faces[i][j]][l] == BC_FIXED) ? BC_FIXED : BC_MATCHEDSLIDE;
        break;
        case(BC_ISOTHERMAL_WALL_FIXED):
        case(BC_ADIABATIC_WALL_FIXED):
        case(BC_POROUS_WALL_FIXED):
        case(BC_OUTLET_FIXED):
        case(BC_INLET_FIXED):
        case(BC_DIRECTSTATE_OUTLET_FIXED):
        case(BC_DIRECTSTATE_INLET_FIXED):
        case(BC_MASSFLOW_OUTLET_FIXED):
        case(BC_MASSFLOW_INLET_FIXED):
        case(BC_SLIP_WALL_FIXED):
          for(int j=0; j<faces[i].numNodes();j++)
            for(int l=0; l<3; l++) dofType[faces[i][j]][l] = BC_FIXED;
          break;
      }
    }
    else {
      if(!isSliding) { // -> constrain the nodes according to face "fluid code"
        switch(faces[i].getCode()) {
          case(BC_SYMMETRY): // by default a symmetry plane is fixed ...
          case(BC_ISOTHERMAL_WALL_FIXED):
          case(BC_ADIABATIC_WALL_FIXED):
          case(BC_POROUS_WALL_FIXED):
          case(BC_OUTLET_FIXED):
          case(BC_INLET_FIXED):
          case(BC_DIRECTSTATE_OUTLET_FIXED):
          case(BC_DIRECTSTATE_INLET_FIXED):
          case(BC_MASSFLOW_OUTLET_FIXED):
          case(BC_MASSFLOW_INLET_FIXED):
          case(BC_SLIP_WALL_FIXED):
          for(int j=0; j<faces[i].numNodes();j++)
            for(int l=0; l<3; l++) dofType[faces[i][j]][l] = BC_FIXED;
          break;
        }
      }
      else { // Sliding is there
        if (it->second->nx != 0.0 ) // it's a sliding surface
          for(int j=0; j<faces[i].numNodes();j++)  dofType[faces[i][j]][0] = BC_FIXED;

        if (it->second->ny != 0.0) // it's a sliding surface
          for(int j=0; j<faces[i].numNodes();j++)  dofType[faces[i][j]][1] = BC_FIXED;

        if (it->second->nz != 0.0) // it's a sliding surface
          for(int j=0; j<faces[i].numNodes();j++)  dofType[faces[i][j]][2] = BC_FIXED;
      }
    } 
  }

  // Step 3. Take into account the matched nodes
  if(matchNodes) {
    for(int i=0; i<matchNodes->size();i++) {
      int inode = matchNodes->subMatchNode(i);
      for(int l=0; l<3; l++) dofType[inode][l] = BC_MATCHED;
    }
  }

  // Step 4. Fill communication buffer (to ensure same constraint on the shared nodes/dofs)
  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<int> nInfo = ntP.getSendBuffer(sndChannel[iSub]);
    int (*buffer)[3] = reinterpret_cast<int (*)[3]>(nInfo.data);
    for (int i = 0; i < sharedNodes->num(iSub); ++i)
      for(int l=0; l<3; l++)
        buffer[i][l] = dofType[(*sharedNodes)[iSub][i]][l];
  }

  return(DofType);
}

//------------------------------------------------------------------------------
//HB: to ensure same constraint on the shared nodes/dofs (for the mesh motion)
void
SubDomain::completeMeshMotionDofType(int* DofType, CommPattern<int> &ntP)
{
  int (*dofType)[3] = reinterpret_cast<int (*)[3]>(DofType);

  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<int> nInfo = ntP.recData(rcvChannel[iSub]);
    int (*buffer)[3] = reinterpret_cast<int (*)[3]>(nInfo.data);
    for (int i = 0; i < sharedNodes->num(iSub); ++i)
      for(int l=0; l<3; l++)
        if(buffer[i][l]!=BC_FREE)
          dofType[(*sharedNodes)[iSub][i]][l] = buffer[i][l];
  }
}

//------------------------------------------------------------------------------

void
SubDomain::completeEmbeddedALEMeshMotionDofType(int* DofType, CommPattern<int> &ntP)
{
  int (*dofType)[3] = reinterpret_cast<int (*)[3]>(DofType);

  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<int> nInfo = ntP.recData(rcvChannel[iSub]);
    int (*buffer)[3] = reinterpret_cast<int (*)[3]>(nInfo.data);
    for (int i = 0; i < sharedNodes->num(iSub); ++i)
      for(int l=0; l<3; l++)
        if(buffer[i][l]!=BC_FREE) {
          if(buffer[i][l] != BC_MATCHEDSLIDE || dofType[(*sharedNodes)[iSub][i]][l] == BC_FREE)
            dofType[(*sharedNodes)[iSub][i]][l] = buffer[i][l];
        }
  }
}

//------------------------------------------------------------------------------

int SubDomain::setFaceToElementConnectivity()
{

  Vec<bool> tagNodes(nodes.size());

  tagNodes = false;

  int i;
  for (i=0; i<faces.size(); ++i) {
    faces[i].setElementNumber(-1, 0);
    faces[i].tagNodesOnBoundaries(tagNodes);
  }

  MapFaces mf;

  int (*fm)[2] = new int[faces.size()][2];
  for (i=0; i<faces.size(); ++i) {
    MaxFace f(faces[i]);
    f.reorder();
    mf[f] = i;
    fm[i][0] = -1;
  }

  int nswap = 0;
  for (i=0; i<elems.size(); ++i)
    if (elems[i].countNodesOnBoundaries(tagNodes) > 2)
      nswap += elems[i].setFaceToElementConnectivity(i, tagNodes, mf, fm);

  for (i=0; i<faces.size(); ++i) {
    faces[i].setElementNumber(fm[i][0], fm[i][1]);
    if (faces[i].getElementNumber() == -1) {
      fprintf(stderr, "*** Error: boundary face %d not connected\n",
	      locToGlobFaceMap[i]+1);
      exit(1);
    }
  }

  delete[] fm;

  return nswap;

}

//------------------------------------------------------------------------------

void SubDomain::getNdAeroLists(int &nInterfNd, int *&interfNd, int &nInfNd,
			       int *&infNd, int &nInternalNd, int *&internalNd, MatchNodeSet* matchNodes)
{
  nInterfNd  = 0;
  nInfNd     = 0;
  nInternalNd= 0;

  if(matchNodes) { //HB: new stuff ...
    nInterfNd = matchNodes->size();
    interfNd = new int[nInterfNd];
    Vec<bool> isMatched(nodes.size());
    isMatched = false;
    for(int i=0;i<matchNodes->size();i++) {
      int inode = matchNodes->subMatchNode(i); // current matched node
      isMatched[inode] = true;
      interfNd[i] = inode;
    }
    for (int i = 0; i < nodes.size(); ++i) {
      if(isMatched[i]) continue;
      if (nodeType[i] == BC_INLET_FIXED || nodeType[i] == BC_OUTLET_FIXED || 
          nodeType[i] == BC_DIRECTSTATE_INLET_FIXED || nodeType[i] == BC_DIRECTSTATE_OUTLET_FIXED ||
          nodeType[i] == BC_MASSFLOW_INLET_FIXED || nodeType[i] == BC_MASSFLOW_OUTLET_FIXED) nInfNd++;
      else nInternalNd++;
    }
    infNd      = new int[nInfNd];
    internalNd = new int[nInternalNd];
    nInfNd     = 0;
    nInternalNd= 0;
    for (int i = 0; i < nodes.size(); ++i) {
      if(isMatched[i]) continue;
      if (nodeType[i] == BC_INLET_FIXED || nodeType[i] == BC_OUTLET_FIXED || 
          nodeType[i] == BC_DIRECTSTATE_INLET_FIXED || nodeType[i] == BC_DIRECTSTATE_OUTLET_FIXED ||
          nodeType[i] == BC_MASSFLOW_INLET_FIXED || nodeType[i] == BC_MASSFLOW_OUTLET_FIXED) infNd[nInfNd++] = i;
      else internalNd[nInternalNd++] = i;
    }
  } else { // Only for ForcedMeshMotion (see ForcedMeshMotionHandler) where matchNodes is not
           // explicitly created -> consider all the node labelled as moving as "matched" nodes
    for (int i = 0; i < nodes.size(); ++i) {
      if (nodeType[i] == BC_INLET_FIXED || nodeType[i] == BC_OUTLET_FIXED || 
          nodeType[i] == BC_DIRECTSTATE_INLET_FIXED || nodeType[i] == BC_DIRECTSTATE_OUTLET_FIXED ||
          nodeType[i] == BC_MASSFLOW_INLET_FIXED || nodeType[i] == BC_MASSFLOW_OUTLET_FIXED) nInfNd++;
      else if (nodeType[i] < BC_INTERNAL) nInterfNd++;
      else nInternalNd++;
    }

    interfNd   = new int[nInterfNd];
    infNd      = new int[nInfNd];
    internalNd = new int[nInternalNd];

    nInterfNd = 0;
    nInfNd = 0;
    nInternalNd = 0;

    for (int i = 0; i < nodes.size(); ++i) {
      if (nodeType[i] == BC_INLET_FIXED || nodeType[i] == BC_OUTLET_FIXED || 
          nodeType[i] == BC_DIRECTSTATE_INLET_FIXED || nodeType[i] == BC_DIRECTSTATE_OUTLET_FIXED ||
          nodeType[i] == BC_MASSFLOW_INLET_FIXED || nodeType[i] == BC_MASSFLOW_OUTLET_FIXED) infNd[nInfNd++] = i;
      else if (nodeType[i] < BC_INTERNAL) interfNd[nInterfNd++] = i;
      else internalNd[nInternalNd++] = i;
    }
  }
}

//------------------------------------------------------------------------------

void SubDomain::testNormals(Vec<Vec3D> &edgeNorm, Vec<double> &edgeNormVel,
			    Vec<Vec3D> &faceNorm, Vec<double> &faceNormVel)
{

  /*
  double nnrm = 0.0;
  double vnrm = 0.0;

  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    int i = locToGlobNodeMap[ edgePtr[l][0] ];
    int j = locToGlobNodeMap[ edgePtr[l][1] ];

    if (locToGlobNodeMap[ edgePtr[l][0] ] > locToGlobNodeMap[ edgePtr[l][1] ]) {
      i = locToGlobNodeMap[ edgePtr[l][1] ];
      j = locToGlobNodeMap[ edgePtr[l][0] ];
    }

    nnrm += edgeNorm[l]*edgeNorm[l];
    vnrm += edgeNormVel[l]*edgeNormVel[l];

    printf("%d %d: %e %e %e %e\n", i+1, j+1,
	   edgeNorm[l][0], edgeNorm[l][1], edgeNorm[l][2], edgeNormVel[l]);

  }

  printf("%d %.10g %.10g\n", locSubNum, nnrm, vnrm);
  */

  for (int i=0; i<faces.size(); ++i) {
    int n0   = locToGlobNodeMap[ faces[i][0] ];
    int n1   = locToGlobNodeMap[ faces[i][1] ];
    int n2   = locToGlobNodeMap[ faces[i][2] ];
    Vec3D  ni    = faces[i].getNormal(faceNorm);
    double nveli = faces[i].getNormalVel(faceNormVel);

    printf("%d %d %d: %e %e %e %e\n", n0+1, n1+1, n2+1,
	   ni[0], ni[1], ni[2], nveli);
  }

}

//------------------------------------------------------------------------------
#ifdef CXFS
void SubDomain::printInfo(FILE *fp)
{

  fprintf(fp, "%d %d %s\n", globSubNum, clusSubNum, suffix);
  fprintf(fp, "%d %d %d\n", nodes.size(), numClusNodes, numNodeRanges);
  for (int i=0; i<numNodeRanges; ++i)
    fprintf(fp, "  %d %d %d\n", nodeRanges[i][0], nodeRanges[i][1], nodeRanges[i][2]);

}
#endif
//------------------------------------------------------------------------------

//HB: create & fill the sliding surface node ownership array which gives for each
//node the sliding surfaces to which this node belongs (using bit)
//For instance, surfOwn[i] = [0 0 1 0 0 1] -> node i belongs to sliding surface 2 & 5
int*
SubDomain::getSlipSurfOwnership(CommPattern<int> &cpat,
                                     map<int,SurfaceData *> &surfaceMap)
{
  int *surfOwn = new int[nodes.size()];
  for(int i = 0; i < nodes.size(); ++i) surfOwn[i] = 0;
#ifdef HB_MESHMOTION_DEBUG
  int nSlidingFaces = 0;
#endif
  for(int i = 0; i < faces.size(); ++i) {
    map<int,SurfaceData *>::iterator it = surfaceMap.find(faces[i].getSurfaceID());
    if(it != surfaceMap.end()) // surface has attribut is the input file
      if(it->second->nx != 0.0 || it->second->ny != 0.0 || it->second->nz != 0.0) { // it's a sliding surface
        int ownFlag = it->second->sBit;                                             //this test is not necessary
        surfOwn[faces[i][0]] |= ownFlag;                                            //if sBit has been properly
        surfOwn[faces[i][1]] |= ownFlag;                                            //initialized to zero for
        surfOwn[faces[i][2]] |= ownFlag;                                            //all the surface
#ifdef HB_MESHMOTION_DEBUG
        nSlidingFaces++;
#endif
      }
  }
#ifdef HB_MESHMOTION_DEBUG
  if(nSlidingFaces) fprintf(stderr," -> in SubDomain::getSlipSurfOwnership: sub %4d has %4d sliding faces.\n",globSubNum,nSlidingFaces);
#endif

  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<int> nInfo = cpat.getSendBuffer(sndChannel[iSub]);
    for (int i = 0; i < sharedNodes->num(iSub); ++i)
      nInfo.data[i] = surfOwn[ (*sharedNodes)[iSub][i] ];
  }

  return(surfOwn);
}

//------------------------------------------------------------------------------
//HB: finish the sliding surface ownership and create the projections for the sliding nodes
//For each "sliding node"i, we perform an orthogonalization (on the fly) of the normal directions
//of the sliding plane this node belongs to in order to ensure orthogonal projections (and so commutattive)
void
SubDomain::createSlipSurfProjection(int*surfOwn, CommPattern<int>&cpat,
                                    BCApplier* bcApplier, SurfaceData** surfData)
{
  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<int> nInfo = cpat.recData(rcvChannel[iSub]);
    for (int i = 0; i < sharedNodes->num(iSub); ++i)
      surfOwn[(*sharedNodes)[iSub][i]] |= nInfo.data[i];
  }
  double normals[3][3]; // no more than 3 independant projection directions at a node
  for(int i = 0; i < nodes.size(); ++i) { // loop over subdomain's nodes
    int surfNum = 0;
    int locOwn = surfOwn[i];
    int numActDir = 0; // count the number of "active" projection directions at this node
    while(locOwn != 0) { // loop over the sliding surfaces
      if((locOwn & 1) != 0) { // this node belong to sliding surface "surfNum"
        double nx = surfData[surfNum]->nx;
        double ny = surfData[surfNum]->ny;
        double nz = surfData[surfNum]->nz;
				double invlen = 1.0/sqrt(nx*nx+ny*ny+nz*nz);
				nx *= invlen; ny *= invlen; nz *= invlen;
        // orthogonalyze with respect to the previous projection directions found at this node
				for(int j = 0; j < numActDir; ++j) {
				  double dot = nx*normals[j][0]+ny*normals[j][1]+nz*normals[j][2];
				  nx -= dot*normals[j][0];
				  ny -= dot*normals[j][1];
				  nz -= dot*normals[j][2];
				}
	      double len = sqrt(nx*nx+ny*ny+nz*nz);
				if(len > 1e-4) {
				  double invlen = 1.0/len;
				  normals[numActDir][0] = invlen*nx;
				  normals[numActDir][1] = invlen*ny;
				  normals[numActDir][2] = invlen*nz;
				  bcApplier->addProj(locSubNum, i, normals[numActDir]);
				  numActDir++;
				  if(numActDir==3) // no more than 3 independant projection directions at a node
				    break;
				}
      }
      surfNum++; locOwn >>= 1; // go to next sliding surface
    }
  }
}

//------------------------------------------------------------------------------
//
int* SubDomain::getRotSurfaceOwnership(CommPattern<int> &cpat,
                                     map<int,SurfaceData *> &surfaceMap)
{
  rotOwn = new int[nodes.size()];

  for(int i = 0; i < nodes.size(); ++i) rotOwn[i] = -1;
  for(int i = 0; i < faces.size(); ++i) {
    map<int,SurfaceData *>::iterator it = surfaceMap.find(faces[i].getSurfaceID());
    if(it != surfaceMap.end()) {
       int rotID = it->second->rotationID; // = -1 if not defined in input file
       if ((rotOwn[faces[i][0]] != -1 && rotOwn[faces[i][0]] != rotID) ||
           (rotOwn[faces[i][1]] != -1 && rotOwn[faces[i][1]] != rotID) ||
           (rotOwn[faces[i][2]] != -1 && rotOwn[faces[i][2]] != rotID))  {
         //fprintf(stderr, " ... WARNING: Node associated to more than 1 Rotation ID\n");
       }
       rotOwn[faces[i][0]] = rotID;
       rotOwn[faces[i][1]] = rotID;
       rotOwn[faces[i][2]] = rotID;
    }
  }
  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<int> nInfo = cpat.getSendBuffer(sndChannel[iSub]);
    for (int i = 0; i < sharedNodes->num(iSub); ++i)
      nInfo.data[i] = rotOwn[ (*sharedNodes)[iSub][i] ];
  }

  return(rotOwn);
}

//------------------------------------------------------------------------------

void SubDomain::completeRotateSurfaceOwnership(CommPattern<int>&cpat)  {

  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<int> nInfo = cpat.recData(rcvChannel[iSub]);
    for (int i = 0; i < sharedNodes->num(iSub); ++i)  {
      if (rotOwn[(*sharedNodes)[iSub][i]] != -1 && nInfo.data[i] != -1 && rotOwn[(*sharedNodes)[iSub][i]] != nInfo.data[i])  {
        fprintf(stderr, " ... WARNING: Node %d in GlSub %d associated with more than 1 rotation\n", (*sharedNodes)[iSub][i], globSubNum);
      }
      rotOwn[(*sharedNodes)[iSub][i]] = std::max(rotOwn[(*sharedNodes)[iSub][i]], nInfo.data[i]);
    }
  }
  //int count = 0;
  //for(int i = 0; i < nodes.size(); ++i)
  //  count += (rotOwn[i]>=0) ? 1 : 0;
  //if(count) { fprintf(stderr, " -> subdomain %3d has %6d 'rotating' nodes\n",globSubNum,count); fflush(stderr); }
}

//------------------------------------------------------------------------------
void SubDomain::multiPade(bcomp *compMat, int *stepParam, double *deltaFreqCoarse, bcomp *padeMat, bcomp *padeVec, int L, int M, int nPoints, double deltaFreqFine, double freqMidFreq, bcomp *snaps, double *freqCoarse)

{
  int size = L+M+1;
  int mFreqStart;
  double mDeltaFreq[nPoints];
  int midFreq[1];
  int nSteps = stepParam[0];
  int numFreqCoarse = stepParam[1];
  int nWrittenSnaps[1];
  nWrittenSnaps[0] = 0;

  if (nPoints == numFreqCoarse) {
    buildPadeMatrix(compMat, stepParam, numFreqCoarse, deltaFreqCoarse, padeMat, L, M, -1);
    buildPadeRhs(compMat, stepParam, numFreqCoarse, padeVec, L, M, -1);
    solveLinearSystem(padeMat, padeVec, size);
    padeSolution(padeVec, stepParam, freqMidFreq, nPoints, L, M, deltaFreqFine, nSteps, snaps,nWrittenSnaps, -1, -1, -1);
  }

  else {
    int nInterval = 1;
    while (1+nInterval*(nPoints-1) <= numFreqCoarse) {
      mFreqStart = multiPointsFreq(nInterval, nPoints, freqCoarse, numFreqCoarse,1);
      multiPointsDeltaFreq(mFreqStart, freqCoarse, nPoints, mDeltaFreq,midFreq);
      buildPadeMatrix(compMat, stepParam, nPoints, mDeltaFreq, padeMat, L, M, nInterval);
      buildPadeRhs(compMat, stepParam, nPoints, padeVec, L, M, nInterval);
      solveLinearSystem(padeMat, padeVec, size);
      padeSolution(padeVec, stepParam, freqCoarse[mFreqStart+(*midFreq)], nPoints, L, M, deltaFreqFine, nSteps, snaps, nWrittenSnaps, nInterval, freqCoarse[mFreqStart], freqCoarse[mFreqStart+nPoints-1]);
      nInterval++;
    }
    if ((numFreqCoarse-1) % (nPoints-1) != 0) {
      mFreqStart = multiPointsFreq(nInterval, nPoints, freqCoarse, numFreqCoarse,0);
      multiPointsDeltaFreq(mFreqStart, freqCoarse, nPoints, mDeltaFreq,midFreq);
      buildPadeMatrix(compMat, stepParam, nPoints, mDeltaFreq, padeMat, L, M, 0);
      buildPadeRhs(compMat, stepParam, nPoints, padeVec, L, M, 0);
      solveLinearSystem(padeMat, padeVec, size);
      padeSolution(padeVec, stepParam, freqCoarse[mFreqStart+(*midFreq)], nPoints, L, M, deltaFreqFine, nSteps, snaps, nWrittenSnaps, nInterval, freqCoarse[mFreqStart], freqCoarse[mFreqStart+nPoints-1]);
    }
    if (nWrittenSnaps[0] != nSteps+1)
      fprintf(stderr," Error : %d snapshots written out of a total of %d : ****\n",nWrittenSnaps[0],nSteps+1);
  }
}
//-----------------------------------------------------------------------------
int SubDomain::multiPointsFreq(int nInterval, int nPoints, double *freqCoarse, int numFreqCoarse, int flag)
{

  if (flag == 1) {
    return (nInterval-1)*(nPoints-1);
  }
  else if (flag == 0)
    return ((numFreqCoarse -1) - (nPoints-1));
  else {
   fprintf(stderr, " Warning : flag different than 0 or 1 in multiPointsFreq\n");
   return(0);
  }
}
//------------------------------------------------------------------------------
void SubDomain::multiPointsDeltaFreq(int mFreqStart, double *freqCoarse, int nPoints, double* mDeltaFreq, int *midFreq)
{
  int i;
  if (nPoints % 2 == 1)
    midFreq[0] = int(floor(nPoints/2.0)+1)-1;
  else
    midFreq[0] = nPoints/2-1;
  for (i=0; i<nPoints; i++) {
    mDeltaFreq[i] = freqCoarse[i+mFreqStart] - freqCoarse[midFreq[0]+mFreqStart];
  }

}
//------------------------------------------------------------------------------
void SubDomain::buildPadeMatrix(bcomp *data, int *stepParam, int nPoints, double *deltaFreq, bcomp *padeMat, int L, int M, int flag)
{
  bcomp oneReal(1.0, 0.0);
  bcomp oneImag(0.0, 1.0);
  bcomp zero(0.0,0.0);

  int nFreqCoarse;
  if (flag < 0)
    nFreqCoarse = stepParam[1];
  else
    nFreqCoarse = nPoints;
  int numPadeDeriv = int(ceil((L+M+1.0)/((double)nFreqCoarse))-1.0);
  double deltaF;
  int size = L+M+1;
  int l,r,k,i;
  bcomp *ddata;

  if (flag > 0)
    ddata = data + ((flag-1)*(nPoints-1))*(numPadeDeriv+1);
  else if (flag < 0)
    ddata = data;
  else if (flag == 0)
    ddata = data + ((stepParam[1]-1) - (nPoints-1))*(numPadeDeriv+1);

  for (i=0; i<size; i++) {
    for (k=0; k<size; k++)
      *(padeMat + size*i + k) = zero;
  }

  for (i=0; i< nFreqCoarse; i++) {
    deltaF = *(deltaFreq+i);
    for (k=0; k < numPadeDeriv+1; k++) {
      if (k*(nFreqCoarse) + i < size) {
        //Construction of the coefficients for the numerator
        for (l=k; l < L+1; l++) {
          *(padeMat + i*(numPadeDeriv+1)*size+ size*k+l) = (-exp(lgamma(l+1))) / (exp(lgamma(l-k+1))) * pow(deltaF,l-k)*oneReal;
        }
        //Construction of the coefficients for the denominator
        for (l=1; l < M+1; l++) {
          *(padeMat + i*(numPadeDeriv+1)*size+ size*k+l+L) = zero;
          for (r=0; r < min(l,k)+1; r++) {
            *(padeMat + i*(numPadeDeriv+1)*size+ size*k+l+L) =  *(padeMat + i*(numPadeDeriv+1)*(L+M+1)+ size*k+l+L) +  exp(lgamma(l+1)) * exp(lgamma(k+1)) / ( exp(lgamma(k-r+1)) * exp(lgamma(l-r+1)) * exp(lgamma(r+1)) ) * pow(deltaF,l-r) * (*(ddata + i*(numPadeDeriv+1) + (k-r)));
          }
        }
      }
    }
  }
}
//--------------------------------------------------------------------------
void SubDomain::buildPadeRhs(bcomp *data, int *stepParam, int nPoints, bcomp *padeVec, int L, int M, int flag)
{
  bcomp oneReal(1.0, 0.0);
  bcomp oneImag(0.0, 1.0);
  bcomp zero(0.0,0.0);

  int nFreqCoarse;
  if (flag < 0)
    nFreqCoarse = stepParam[1];
  else
    nFreqCoarse =  nPoints;
  int numPadeDeriv = int(ceil((L+M+1.0)/((double)nFreqCoarse))-1.0);
  int i,j;
  bcomp *ddata;

  if (flag > 0)
    ddata = data + (flag-1)*(nPoints-1)*(numPadeDeriv+1);
  else if (flag < 0)
    ddata = data;
  else if (flag == 0)
    ddata = data + ((stepParam[1]-1) - (nPoints-1))*(numPadeDeriv+1);
  for (i=0; i < nFreqCoarse; i++) {

    for (j=0; j < numPadeDeriv+1; j++) {
      if (j*nFreqCoarse + i < L+M+1)
        *(padeVec + i*(numPadeDeriv+1) + j) = -(*(ddata + (numPadeDeriv+1)*i + j));
    }

  }
}
//--------------------------------------------------------------------------
void SubDomain::padeSolution(bcomp *padeVec, int *stepParam, double midFreq, int nPoints, int L, int M, double deltaFreq, int nSteps, bcomp *snaps, int* nWrittenSnaps, int flag, double freqMin, double freqMax)
{
  int i,iSteps,k,l;
  bcomp numPade(0.0, 0.0);
  bcomp denPade(1.0, 0.0);
  bcomp oneReal(1.0, 0.0);
  bcomp oneImag(0.0, 1.0);
  bcomp zero(0.0,0.0);
  double freq = 0.0;
  int nFreqCoarse = stepParam[1];
  int startSnaps = nWrittenSnaps[0];
  if (flag < 0) {

    for (iSteps=0; iSteps<nSteps+1; iSteps++) {
      freq = iSteps * deltaFreq;
      for (k=0; k<L+1; k++)
        numPade = numPade + (*(padeVec + k)) * pow(freq - midFreq,k);
      for (l=1; l<M+1; l++)
        denPade = denPade + (*(padeVec + L + l)) * pow(freq - midFreq,l);
      *(snaps + iSteps) = numPade / denPade;
      numPade = zero;
      denPade = oneReal;
    }

  }
  else {
    for (iSteps=startSnaps; iSteps<nSteps+1; iSteps++) {
      freq = iSteps * deltaFreq;
      if ((freq >= freqMin || flag==1)  && (freq < freqMax || flag==int(ceil((nFreqCoarse-1.0)/(nPoints-1.0)))) ) {
        for (k=0; k<L+1; k++)
          numPade = numPade + (*(padeVec + k)) * pow(freq - midFreq,k);
        for (l=1; l<M+1; l++)
          denPade = denPade + (*(padeVec + L + l)) * pow(freq - midFreq,l);


        *(snaps + iSteps) = numPade / denPade;
        nWrittenSnaps[0]++;
        numPade = zero;
        denPade = oneReal;
      }
    }
  }
}
//------------------------------------------------------------------------
void SubDomain::solveLinearSystem(bcomp *a, bcomp *b, int n) {


  DenseMatrixOp<bcomp, 5, 5>::lu(a,b,n);


}
//--------------------------------------------------------------------------
void SubDomain::buildDeltaFreq(double *deltaFreqCoarse,int numFreqCoarse, double *freqCoarse, int *midFreq){


  int i;
  if (numFreqCoarse % 2 == 1)
    midFreq[0] = int(floor(numFreqCoarse/2.0)+1)-1;
  else
    midFreq[0] = numFreqCoarse/2-1;
  for (i=0; i<numFreqCoarse; i++)
    deltaFreqCoarse[i] = freqCoarse[i] - freqCoarse[midFreq[0]];


}
//--------------------------------------------------------------------------
void SubDomain::extractElementsRelativeToAComponentAndAMode(double* tempMat, bcomp* compMat, int iDim, int iStrMode,int numPadeDeriv, int numFreqCoarse, int nStrMode, int nSnapsCoarse, double freq1)
{
  int i,j,start1,start2,start,startF;
  bcomp oneReal(1.0, 0.0);
  bcomp oneImag(0.0, 1.0);
  if (freq1 == 0.0) {
    start1 = iStrMode*(2*(numPadeDeriv+1)-1); //starting point
    start2 = nStrMode*(2*(numPadeDeriv+1)-1);
  }
  else {
    start1 = iStrMode*2*(numPadeDeriv+1);
    start2 = nStrMode*2*(numPadeDeriv+1);
  }
  // for each freq of the coarse grid
  for (i = 0; i<numFreqCoarse; i++) {



    if ( i == 0 ) {
      *(compMat) = oneReal*(*(tempMat +iDim*nSnapsCoarse + start1));
      if (freq1 != 0.0)
       *(compMat) = *(compMat) + oneImag*(*(tempMat + iDim*nSnapsCoarse +start1+1));



      if (numPadeDeriv > 0) {
        for (j = 1; j<numPadeDeriv+1; j++) {
          if (freq1 == 0.0)
            *(compMat+j) = oneReal*(*(tempMat +iDim*nSnapsCoarse + start1+(2*j-1))) + oneImag*(*(tempMat+iDim*nSnapsCoarse + start1+2*j));
          else
            *(compMat+j) = oneReal*(*(tempMat+iDim*nSnapsCoarse + start1+2*j)) + oneImag*(*(tempMat+ iDim*nSnapsCoarse + start1+2*j+1));



        }
      }
    }
    else if ( i > 0 ) {
      startF = nStrMode*2*(numPadeDeriv+1);
      start = iStrMode*2*(numPadeDeriv+1);
      *(compMat+ i*(numPadeDeriv+1)) = oneReal*(*(tempMat+iDim*nSnapsCoarse + start2+(i-1)*startF+start)) + oneImag*(*(tempMat+iDim*nSnapsCoarse +start2+(i-1)*startF+start+1));

     if (numPadeDeriv > 0) {
       for (j = 1; j<numPadeDeriv+1; j++)


          *(compMat+i*(numPadeDeriv+1)+j) = oneReal*(*(tempMat+iDim*nSnapsCoarse +start2+(i-1)*startF+start+2*j)) + oneImag*(*(tempMat+iDim*nSnapsCoarse +start2+(i-1)*startF+start+2*j+1));
     }
    }
  }
}
//--------------------------------------------------------------------------

void SubDomain::finalizeTags(SVec<int,2> &tag)
{

  int (*edgePtr)[2] = edges.getPtr();
  for (int l=0; l<edges.size(); ++l) {
    int i = edgePtr[l][0];
    int j = edgePtr[l][1];
    tag[i][1] = tag[i][0] +  tag[j][0];
    tag[j][1] = tag[i][1];
  }

}

//------------------------------------------------------------------------------

// Included (MB)
void SubDomain::checkVec(SVec<double,3> &V)
{

  for (int i=0; i<V.size(); ++i) 
  {

    if ((nodeType[i] != BC_ADIABATIC_WALL_MOVING)  && (nodeType[i] != BC_ISOTHERMAL_WALL_MOVING))  
    {
      if ((V[i][0] != 0.0) || ((V[i][1] != 0.0) || (V[i][2] != 0.0))) 
      {
        fprintf(stderr,"*** Error: Vector dXdsb is different from zero");
        fprintf(stderr," at a point %d in the interior of the mesh\n", i);
//        exit(1);
      }
      //----
      if (V[i][0] != 0.0)
      {
        std::cerr << "*** Overwrite the value in the x-direction to 0\n";
        V[i][0] = 0.0;
      }
      if (V[i][1] != 0.0)
      {
        std::cerr << "*** Overwrite the value in the y-direction to 0\n";
        V[i][1] = 0.0;
      }
      if (V[i][2] != 0.0)
      {
        std::cerr << "*** Overwrite the value in the z-direction to 0\n";
        V[i][2] = 0.0;
      }

    }
  }

}

//--------------------------------------------------------------------------

void SubDomain::computeTetsConnectedToNode(Vec<int> &Ni)
{

  for (int tetNum=0; tetNum < elems.size(); ++tetNum)
      for (int i=0; i<4; ++i)
        ++Ni[elems[tetNum][i]];

}

//--------------------------------------------------------------------------

void SubDomain::computeLij(double Lij[3][3], double r_u[3], double r_u_u[6], double vc[5])
{

  if (vc[0] < 1.0e-7) vc[0] = 1.0e-7;

  Lij[0][0] = r_u_u[0] - (1.0/vc[0])*(r_u[0]*r_u[0]);
  Lij[0][1] = r_u_u[1] - (1.0/vc[0])*(r_u[0]*r_u[1]);
  Lij[0][2] = r_u_u[2] - (1.0/vc[0])*(r_u[0]*r_u[2]);
  Lij[1][1] = r_u_u[3] - (1.0/vc[0])*(r_u[1]*r_u[1]);
  Lij[1][2] = r_u_u[4] - (1.0/vc[0])*(r_u[1]*r_u[2]);
  Lij[2][2] = r_u_u[5] - (1.0/vc[0])*(r_u[2]*r_u[2]);
  Lij[2][0] = Lij[0][2];
  Lij[2][1] = Lij[1][2];
  Lij[1][0] = Lij[0][1];
  Lij[0][0] = Lij[0][0] - (1.0/3.0)*(Lij[0][0]+Lij[1][1]+Lij[2][2]);
  Lij[1][1] = Lij[1][1] - (1.0/3.0)*(Lij[0][0]+Lij[1][1]+Lij[2][2]);
  Lij[2][2] = Lij[2][2] - (1.0/3.0)*(Lij[0][0]+Lij[1][1]+Lij[2][2]);

}

//-----------------------------------------------------------------------

void SubDomain::computeBij(double Bij[3][3], double r_s_p[6], double sqrt2S2,
                           double Pij[3][3], double sq_rat_delta, double vc[5])


{

   Bij[0][0] = r_s_p[0] - sq_rat_delta*vc[0]*sqrt2S2*Pij[0][0];
   Bij[0][1] = r_s_p[3] - sq_rat_delta*vc[0]*sqrt2S2*Pij[0][1];
   Bij[0][2] = r_s_p[4] - sq_rat_delta*vc[0]*sqrt2S2*Pij[0][2];
   Bij[1][1] = r_s_p[1] - sq_rat_delta*vc[0]*sqrt2S2*Pij[1][1];
   Bij[1][2] = r_s_p[5] - sq_rat_delta*vc[0]*sqrt2S2*Pij[1][2];
   Bij[2][2] = r_s_p[2] - sq_rat_delta*vc[0]*sqrt2S2*Pij[2][2];
   Bij[1][0] = Bij[0][1];
   Bij[2][0] = Bij[0][2];
   Bij[2][1] = Bij[1][2];

}

//--------------------------------------------------------------------------

void SubDomain::computeZi(double Zi[3], double sq_rat_delta, double sqrt2S2, double dtdxj[3],
			  double r_s_dtdxj[3], double vc[5], double Cp)
{

  Zi[0] = Cp * (sq_rat_delta*vc[0]*sqrt2S2*dtdxj[0] - r_s_dtdxj[0]);
  Zi[1] = Cp * (sq_rat_delta*vc[0]*sqrt2S2*dtdxj[1] - r_s_dtdxj[1]);
  Zi[2] = Cp * (sq_rat_delta*vc[0]*sqrt2S2*dtdxj[2] - r_s_dtdxj[2]);

}

//--------------------------------------------------------------------------

void SubDomain::computeLi(double Li[3], double r_e, double r_e_plus_p, double r_u[3], double vc[5])
{

  if (vc[0] < 0.0000001) vc[0] = 0.0000001;
  Li[0] = (r_e+vc[4])*(r_u[0]/vc[0]) - (r_e_plus_p)*vc[1];
  Li[1] = (r_e+vc[4])*(r_u[1]/vc[0]) - (r_e_plus_p)*vc[2];
  Li[2] = (r_e+vc[4])*(r_u[2]/vc[0]) - (r_e_plus_p)*vc[3];

}

//--------------------------------------------------------------------------

void SubDomain::findNodeBoundingBoxes(SVec<double,3>&X, SVec<double,3> &Xmin, SVec<double,3> &Xmax)
{
  for (int iNode=0; iNode<nodes.size(); iNode++) {
    Xmin[iNode][0] = Xmax[iNode][0] = X[iNode][0];
    Xmin[iNode][1] = Xmax[iNode][1] = X[iNode][1];
    Xmin[iNode][2] = Xmax[iNode][2] = X[iNode][2];
  }
  int (*ptr)[2] = edges.getPtr();
  bool* edgeMasterFlag = edges.getMasterFlag();
  int p, q;

  for (int l=0; l<edges.size(); l++) {
    if (!edgeMasterFlag[l])
      continue;
    p = ptr[l][0];
    q = ptr[l][1];
    for (int k=0; k<3; k++) {
      if (X[q][k]<Xmin[p][k])  Xmin[p][k] = X[q][k];
      if (X[p][k]<Xmin[q][k])  Xmin[q][k] = X[p][k];
    }
    for (int k=0; k<3; k++) {
      if (X[q][k]>Xmax[p][k])  Xmax[p][k] = X[q][k];
      if (X[p][k]>Xmax[q][k])  Xmax[q][k] = X[p][k];
    }
  }
}

//--------------------------------------------------------------------------

void SubDomain::findSubDomainBoundingBoxes(SVec<double,3>&X, double *Xmin, double *Xmax)
{
  SVec<double,3> BBmin(nodes.size());
  SVec<double,3> BBmax(nodes.size());
  findNodeBoundingBoxes(X, BBmin, BBmax);
  
  for (int i=0; i<3; i++) {
    Xmin[i] = BBmin[0][i];
    Xmax[i] = BBmax[0][i];
  }

  for (int i=0; i<nodes.size(); i++) 
    for (int j=0; j<3; j++) {
      if (Xmin[j]>BBmin[i][j])  Xmin[j] = BBmin[i][j];
      if (Xmax[j]<BBmax[i][j])  Xmax[j] = BBmax[i][j];
    }
}

//--------------------------------------------------------------------------

void SubDomain::outputCsDynamicLES(DynamicLESTerm *dles, SVec<double,2> &Cs,
                                   SVec<double,3> &X, Vec<double> &CsVal)
{

  for (int tetNum=0; tetNum < elems.size(); ++tetNum) {
    double dp1dxj[4][3];
    double vol = elems[tetNum].computeGradientP1Function(X, dp1dxj);
    double cs[4] = {Cs[elems[tetNum][0]][0], Cs[elems[tetNum][1]][0],
                    Cs[elems[tetNum][2]][0], Cs[elems[tetNum][3]][0]};

    double csval = dles->outputCsValues(vol, cs, X, elems[tetNum]);
    for (int i=0; i<4; ++i)
      CsVal[elems[tetNum][i]] += csval * vol;
  }

}

//--------------------------------------------------------------------------
void SubDomain::changeSurfaceType(map<int,SurfaceData*>& surfaceMap, map<int,BoundaryData*>& bcMap)  {
for (int i=0;i<faces.size(); i++) { // Loop over faces
    map<int,SurfaceData*>::iterator it = surfaceMap.find(faces[i].getSurfaceID());
    if(it!=surfaceMap.end()) { // surface has attribut in the input file
      if(it->second->type == SurfaceData::ADIABATIC) {
        if(faces[i].getCode() == BC_ISOTHERMAL_WALL_MOVING)
          faces[i].setType(BC_ADIABATIC_WALL_MOVING);
        if(faces[i].getCode() == BC_ISOTHERMAL_WALL_FIXED)
          faces[i].setType(BC_ADIABATIC_WALL_FIXED);
        }
      if(it->second->type == SurfaceData::ISOTHERMAL) {
if(faces[i].getCode()!=-1)
        if(faces[i].getCode() == BC_ADIABATIC_WALL_MOVING){
          faces[i].setType(BC_ISOTHERMAL_WALL_MOVING);
        }
        if(faces[i].getCode() == BC_ADIABATIC_WALL_FIXED){
          faces[i].setType(BC_ISOTHERMAL_WALL_FIXED);
        }
      }

      map<int,BoundaryData*>::iterator it2 = bcMap.find(it->second->bcID);
      if(it2!=bcMap.end()) {
        if(it2->second->type == BoundaryData::DIRECTSTATE) {
          if(faces[i].getCode() == BC_INLET_MOVING)
            faces[i].setType(BC_DIRECTSTATE_INLET_MOVING);
          if(faces[i].getCode() == BC_INLET_FIXED)
            faces[i].setType(BC_DIRECTSTATE_INLET_FIXED);
          if(faces[i].getCode() == BC_OUTLET_MOVING)
            faces[i].setType(BC_DIRECTSTATE_OUTLET_MOVING);
          if(faces[i].getCode() == BC_OUTLET_FIXED)
            faces[i].setType(BC_DIRECTSTATE_OUTLET_FIXED);
        }
        if(it2->second->type == BoundaryData::MASSFLOW) {
          if(faces[i].getCode() == BC_INLET_MOVING)
            faces[i].setType(BC_MASSFLOW_INLET_MOVING);
          if(faces[i].getCode() == BC_INLET_FIXED)
            faces[i].setType(BC_MASSFLOW_INLET_FIXED);
          if(faces[i].getCode() == BC_OUTLET_MOVING)
            faces[i].setType(BC_MASSFLOW_OUTLET_MOVING);
          if(faces[i].getCode() == BC_OUTLET_FIXED)
            faces[i].setType(BC_MASSFLOW_OUTLET_FIXED);
        }
        if(it2->second->type == BoundaryData::POROUSWALL) {
          if(faces[i].getCode() == BC_SLIP_WALL_MOVING)
            faces[i].setType(BC_POROUS_WALL_MOVING);
          if(faces[i].getCode() == BC_SLIP_WALL_FIXED)
            faces[i].setType(BC_POROUS_WALL_FIXED);
        }
      }
    }
  }
}

//--------------------------------------------------------------------------
void SubDomain::markFaceBelongsToSurface(Vec<int> &faceFlag, CommPattern<int> &cpat) {

  for (int i = 0; i < faces.size(); ++i) {
    for (int j=0; j < faces[i].numNodes(); ++j) {
      int surfNum = faces[i].getSurfaceID()+1;
      const Face &face = faces[i];
      faceFlag[face.nodeNum(j)] |= 1 << surfNum;
    }
  }

  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<int> nInfo = cpat.getSendBuffer(sndChannel[iSub]);
    for (int i = 0; i < sharedNodes->num(iSub); ++i)
      nInfo.data[i] = faceFlag[ (*sharedNodes)[iSub][i] ];
  }
}

//--------------------------------------------------------------------------

void SubDomain::completeFaceBelongsToSurface(Vec<int> &ndToSurfFlag, Vec<double> &nodeTemp, map<int,SurfaceData*>& surfaceMap, CommPattern<int> &cpat) {
   for (int iSub = 0; iSub < numNeighb; ++iSub) {
     SubRecInfo<int> nInfo = cpat.recData(rcvChannel[iSub]);
     for (int i = 0; i < sharedNodes->num(iSub); ++i)
       ndToSurfFlag[(*sharedNodes)[iSub][i]] |= nInfo.data[i];
   }

   for(int i = 0; i < ndToSurfFlag.size(); ++i)
     if(ndToSurfFlag[i] != 0) {
       int count = (ndToSurfFlag[i] & 1) != 0 ? 1 : 0;
       double totTemp =  (ndToSurfFlag[i] & 1) != 0 ? nodeTemp[i] : 0.0;
       for(int j = 1; j < 8*sizeof(int); j++) {
         if( ndToSurfFlag[i] & (1 << j) ) {
             map<int,SurfaceData*>::iterator it = surfaceMap.find(j-1);
             if(it == surfaceMap.end())
               continue;
             if(it->second->type == SurfaceData::ISOTHERMAL) {
               if(it->second->temp > 0){ 
               totTemp += (it->second->temp);
               count++;
             }
               }
         }
       }
       if(count != 0){
       nodeTemp[i] = totTemp/count; 
       }
    }
}

//------------------------------------------------------------------------------
//             SOLID LEVEL SET SOLUTION                                      ---
//------------------------------------------------------------------------------

//-----------------------------------------------------------------------------

void SubDomain::getTriangulatedSurfaceFromFace( SVec<double,3> &X, TriangulatedSurface *triaSurf2 )
{
   triaSurf2->addSurfaceFromFace(X, faces);
}

//----------------------------------------------------------------------------

void SubDomain::getTriangulatedSurfaceFromFace( TriangulatedSurface *triaSurf2 )
{
  triaSurf2->addSurfaceFromFace(faces);
}

//--------------------------------------------------------------------------

bool SubDomain::isINodeinITet(Vec3D nodeCoord, int iTet, SVec<double,3> &X)
{
  int* nodeNo = elems[iTet].nodeNum();
  Vec3D tetCoord[4];
  for (int i=0; i<4; i++)
    tetCoord[i] = X[nodeNo[i]];

  //a fast check.
  double min[3], max[3];
  for (int i=0; i<3; i++) {
    if (tetCoord[0][i]<tetCoord[1][i]) {min[i] = tetCoord[0][i]; max[i] = tetCoord[1][i];}
    else {min[i] = tetCoord[1][i];  max[i] = tetCoord[0][i];}
    if (min[i]>tetCoord[2][i]) min[i] = tetCoord[2][i];
    if (min[i]>tetCoord[3][i]) min[i] = tetCoord[3][i];
    if (max[i]<tetCoord[2][i]) max[i] = tetCoord[2][i];
    if (max[i]<tetCoord[3][i]) max[i] = tetCoord[3][i];
  }

  if (nodeCoord[0] < min[0] || nodeCoord[0] > max[0] || nodeCoord[1] < min[1] || nodeCoord[1] > max[1] ||
      nodeCoord[2] < min[2] || nodeCoord[2] > max[2] )  return false;

  //a careful check
  double vol1, vol2;
  vol1 = orient3dslow(tetCoord[0], tetCoord[1], tetCoord[2], tetCoord[3]);
  vol2 = orient3dslow(tetCoord[0], tetCoord[1], tetCoord[2], nodeCoord);
  if (vol1*vol2<0.0) return false;
  vol1 = orient3dslow(tetCoord[1], tetCoord[2], tetCoord[3], tetCoord[0]);
  vol2 = orient3dslow(tetCoord[1], tetCoord[2], tetCoord[3], nodeCoord);
  if (vol1*vol2<0.0) return false;
  vol1 = orient3dslow(tetCoord[2], tetCoord[3], tetCoord[0], tetCoord[1]);
  vol2 = orient3dslow(tetCoord[2], tetCoord[3], tetCoord[0], nodeCoord);
  if (vol1*vol2<0.0) return false;
  vol1 = orient3dslow(tetCoord[3], tetCoord[0], tetCoord[1], tetCoord[2]);
  vol2 = orient3dslow(tetCoord[3], tetCoord[0], tetCoord[1], nodeCoord);
  if (vol1*vol2<0.0) return false;
/*[for debug]
  fprintf(stderr, "element number = %d + 1.\n", iTet);
  fprintf(stderr, "coords of node 1: [ %.4f,  %.4f,  %.4f]\n", tetCoord[0][0], tetCoord[0][1], tetCoord[0][2]);
  fprintf(stderr, "coords of node 2: [ %.4f,  %.4f,  %.4f]\n", tetCoord[1][0], tetCoord[1][1], tetCoord[1][2]);
  fprintf(stderr, "coords of node 3: [ %.4f,  %.4f,  %.4f]\n", tetCoord[2][0], tetCoord[2][1], tetCoord[2][2]);
  fprintf(stderr, "coords of node 4: [ %.4f,  %.4f,  %.4f]\n", tetCoord[3][0], tetCoord[3][1], tetCoord[3][2]);
*/
  return true;
}

//--------------------------------------------------------------------------------

void SubDomain::localCoord(Vec3D image, int iTet, SVec<double,3> &X, Vec3D& localcoords)
{
  int* nodeNo = elems[iTet].nodeNum();
  Vec3D tetCoord[4];
  for (int i=0; i<4; i++)
    tetCoord[i] = X[nodeNo[i]];
  Vec3D coord1 = tetCoord[1] - tetCoord[0];
  Vec3D coord2 = tetCoord[2] - tetCoord[0];
  Vec3D coord3 = tetCoord[3] - tetCoord[0];

  double a11, a12, a13, a21, a22, a23, a31, a32, a33;
  double b1, b2, b3;
  double c11, c12, c13, c21, c22, c23, c31, c32, c33;
  double detA;      //info. for solving a 3 by 3 linear system.

  a11 = coord1[0];  a21 = coord1[1];  a31 = coord1[2];
  a12 = coord2[0];  a22 = coord2[1];  a32 = coord2[2];
  a13 = coord3[0];  a23 = coord3[1];  a33 = coord3[2];

  detA = a11*(a22*a33 - a32*a23) - a21*(a12*a33 - a32*a13) + a31*(a12*a23 - a22*a13);
  if (detA == 0.0) {fprintf(stderr,"error message: det(A) = 0 in SubDomain::localCoord(...) \n");  exit(1);}


  c11 = 1.0/detA*(a22*a33-a32*a23);  c21 = 1.0/detA*(a23*a31-a33*a21);  c31 = 1.0/detA*(a21*a32-a31*a22);
  c12 = 1.0/detA*(a13*a32-a33*a12);  c22 = 1.0/detA*(a11*a33-a31*a13);  c32 = 1.0/detA*(a12*a31-a32*a11);
  c13 = 1.0/detA*(a12*a23-a22*a13);  c23 = 1.0/detA*(a13*a21-a23*a11);  c33 = 1.0/detA*(a11*a22-a12*a21);

  b1 = image[0]-tetCoord[0][0];  b2 = image[1]-tetCoord[0][1];  b3 = image[2]-tetCoord[0][2];

  localcoords[0] = c11*b1 + c12*b2 + c13*b3;
  localcoords[1] = c21*b1 + c22*b2 + c23*b3;
  localcoords[2] = c31*b1 + c32*b2 + c33*b3;

//  fprintf(stderr,"Tet: %d, Coords: [%.4f, %.4f, %.4f].\n\n\n", iTet, localcoords[0], localcoords[1], localcoords[2]);
}

//-------------------------------------------------------------------------------------

int* SubDomain::getNeiElemOfNode(int iNode, int depth, int& size )
{
  if (!NodeToElem)  NodeToElem = createNodeToElementConnectivity();
  size = NodeToElem->num(iNode);

  int *list = new int[elems.size()];

  for (int i=0; i<size; i++){
    list[i] = (*NodeToElem)[iNode][i];
  }
  if (depth == 1){
    int *finalList = new int[size];
    for (int i=0; i<size; i++)  finalList[i] = list[i];

    delete[] list;
    return finalList;
  }

  int *list2 = new int[elems.size()];
  int size2 = 0;
  if (!ElemToElem) ElemToElem = createElementToElementConnectivity();

  for (int count=0; count<depth-1; count++){
    size2 = 0;
    for (int i=0; i<size; i++){
      if (size2>0){
        for (int j=0; j<ElemToElem->num(list[i]); j++){
          bool newElem = true;
          for (int listNode=0; listNode<size2; listNode++){
            if (list2[listNode] ==(* ElemToElem)[list[i]][j]) {newElem = false;  break;}
          }
          if (newElem == true)  list2[size2++] = (*ElemToElem)[list[i]][j];
        }
      }
      else {for (int l=0; l<ElemToElem->num(list[i]); l++)  list2[size2++] = (*ElemToElem)[list[i]][l]; }
    }
    int *list_temp = list;
    list = list2;
    delete[] list_temp;
    size = size2;
  }

  int *finalList = new int[size];
  for (int i=0; i<size; i++)
    finalList[i] = list[i];
  delete[] list;
  return finalList;
}

//--------------------------------------------------------------------------------
int* SubDomain::getNeiNodeOfNode(int iNode, int& size)
{

	if(!NodeToNode) NodeToNode = createEdgeBasedConnectivity();

	size = NodeToNode->num(iNode);

	int *list = new int[nodes.size()];
	
	for(int i=0; i<size; i++)	list[i] = (*NodeToNode)[iNode][i];

	int *finalList = new int[size];

   for(int i=0; i<size; i++) finalList[i] = list[i];

	delete[] list;

	return finalList;
	
}

//--------------------------------------------------------------------------------

void SubDomain::getNodeCoords(int iNode, SVec<double,3> &X, double& x, double& y, double& z)
{
  x = X[iNode][0];  y = X[iNode][1];  z = X[iNode][2];
}

//--------------------------------------------------------------------------------

double SubDomain::scalarNormalExtrap(double* value, Vec3D target, Vec3D normal, int myTet, SVec<double,3> &X, bool check = 0)
{
  if (check) { //check if target lies in iTet
    bool correct = isINodeinITet(target, myTet, X);
    if (!correct) {fprintf(stderr,"ERROR. Position lies outside the tetrahedron. Abort.\n"); exit(-1);}
  }

  double res = 0.0;
  bool found = false;

  int* myNodes = elems[myTet].nodeNum();

  for (int iFace=0; iFace<4; iFace++) {
    int nodeA = myNodes[iFace%4], nodeB = myNodes[(iFace+1)%4], nodeC = myNodes[(iFace+2)%4];
    Vec3D A(X[nodeA][0], X[nodeA][1], X[nodeA][2]);
    Vec3D B(X[nodeB][0], X[nodeB][1], X[nodeB][2]);
    Vec3D C(X[nodeC][0], X[nodeC][1], X[nodeC][2]);

    Vec3D faceNormal = (B-A)^(C-A);  faceNormal = 1.0/faceNormal.norm()*faceNormal;
    double r = ((A-target)*faceNormal)/(faceNormal*normal);
    if (r<0) continue; // opposite direction.
    Vec3D P = target + r*normal;

    //debug for P.
    if ((P-A)*faceNormal>1e-10) {fprintf(stderr,"dot = %f, error in projection calculation. Abort.\n,",(P-A)*faceNormal); exit(-1);}

    //get barycentric coords.
    double rA = 0.5*((B-P)^(C-P)).norm(),  rB = 0.5*((C-P)^(A-P)).norm(),  rC = 0.5*((A-P)^(B-P)).norm();
    double rU = 0.5*((B-A)^(C-A)).norm();
    if (rA<0.0) rA=-rA;  if (rB<0.0) rB=-rB;  if (rC<0.0) rC=-rC;  if (rU<0.0) rU=-rU;
    rA /= rU;  rB /= rU;  rC /= rU;


    if ((rA+rB+rC-1.0)>1e-10 || (rA+rB+rC-1.0)<-1e-10) continue;
    found = true;
    res = value[iFace%4]*rA + value[(iFace+1)%4]*rB + value[(iFace+2)%4]*rC;
    break;
  }

  if (!found) {fprintf(stderr,"error in normalExtrap(...). Abort.\n"); fprintf(stderr,"target: [%e, %e, %e], normal[%e,%e,%e], C1[%e,%e,%e], C2[%e,%e,%e], C3[%e,%e,%e], C4[%e,%e,%e]. \n", target[0], target[1], target[2], normal[0], normal[1], normal[2], X[myNodes[0]][0], X[myNodes[0]][1], X[myNodes[0]][2], X[myNodes[1]][0], X[myNodes[1]][1], X[myNodes[1]][2], X[myNodes[2]][0], X[myNodes[2]][1], X[myNodes[2]][2], X[myNodes[3]][0], X[myNodes[3]][1], X[myNodes[3]][2]); exit(-1);}
  return res;
}

//--------------------------------------------------------------------------------

void SubDomain::updateNodeTag(SVec<double,3> &X, LevelSetStructure &LSS, Vec<int> &nodeTag0, Vec<int> &nodeTag)
{
  // EdgeLoop(X, LSS, nodeTag0, nodeTag);
  // NodeLoop
}

//--------------------------------------------------------------------------------

void SubDomain::computeCharacteristicEdgeLength(SVec<double,3>&X, double& minLength, double& aveLength, double& maxLength, int& numInsideEdges, const double xmin, const double xmax, const double ymin, const double ymax, const double zmin, const double zmax)
#define EDGE_LENGTH
#ifdef EDGE_LENGTH
{
  edges.computeCharacteristicEdgeLength(X,minLength,aveLength,maxLength, numInsideEdges, xmin, xmax, ymin, ymax, zmin, zmax);
}
#else
{ fprintf(stderr,"EDGE_LENGTH not defined! SubDomain::computeCharacteristicEdgeLength did noting!\n");}
#endif

//---------------------------------------------------------------------------------------------

bool SubDomain::insideOutside(double *position, const double xmin, const double xmax,
                     const double ymin, const double ymax, const double zmin, const double zmax)
{
  if (position[0] < xmin || position[0] > xmax || position[1] < ymin || position[1] > ymax ||
      position[2] < zmin || position[2] > zmax )  return false;
  return true;
}

//-----------------------------------------------------------------------------------------------

void SubDomain::getMeshInBoundingBox(SVec<double,3> &X, const double xmin, const double xmax,
                                       const double ymin, const double ymax,
                                       const double zmin, const double zmax,
                                       int* nodeTag, int& numChosenNodes, int* tempNodeList,
                                       int& numChosenElems, int (*tempElemList)[4])
{
  //1. build a temporary element list that contains only elem indices.
  numChosenElems = 0;

  for (int i=0; i<numElems(); i++) {
    int* elemNodeNum = getElemNodeNum(i);
    // check if any node is (1)inside bounding box & (2) within the "band".
    if ((this->insideOutside(X[elemNodeNum[0]], xmin, xmax, ymin, ymax, zmin, zmax) ||
        this->insideOutside(X[elemNodeNum[1]], xmin, xmax, ymin, ymax, zmin, zmax)  ||
        this->insideOutside(X[elemNodeNum[2]], xmin, xmax, ymin, ymax, zmin, zmax)  ||
        this->insideOutside(X[elemNodeNum[3]], xmin, xmax, ymin, ymax, zmin, zmax) ))
    {
      tempElemList[numChosenElems][0] = i;
      numChosenElems++;
    }
  }

  //2. construct a tag on nodes. build a temporary node list.
  numChosenNodes = 0;
  for (int i=0; i<numNodes(); i++) nodeTag[i] = -1;

  for (int iElem=0; iElem<numChosenElems; iElem++) {
    int* elemNodeNum = getElemNodeNum(tempElemList[iElem][0]);
    for (int j=0; j<4; j++)
      if (nodeTag[elemNodeNum[j]] < 0) {
        nodeTag[elemNodeNum[j]] = numChosenNodes;
        tempNodeList[numChosenNodes] = elemNodeNum[j];
        numChosenNodes++;
      }
  }

  //3. fill the element list by the node indices of each element.
  for (int i=0; i<numChosenElems; i++) {
    int* elemNodeNum = getElemNodeNum(tempElemList[i][0]);
    for (int j=0; j<4; j++) tempElemList[i][j] = nodeTag[elemNodeNum[j]];
  }

  //fprintf(stderr,"numChosenElems = %d, numChosenNodes = %d.\n", numChosenElems, numChosenNodes);
}


//-----------------------------------------------------------------------------------------------
/*
void SubDomain::computeCVBasedForceLoad(int forceApp, int orderOfAccuracy, GeoState& geoState,
                                        SVec<double,3> &X, double (*Fs)[3], int sizeFs,
                                        LevelSetStructure &LSS, Vec<double> &pstarij, 
                                        Vec<double> &pstarji, double pInfty)
{
  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  bool* masterFlag = edges.getMasterFlag();
  int (*ptr)[2];
  int i,j;
  ptr = edges.getPtr();
  bool iActive, jActive, intersect;
  if (forceApp!=1) {fprintf(stderr,"ERROR: force method (%d) not recognized! Abort..\n", forceApp); exit(-1);}

  for (int l=0; l<edges.size(); l++) {
    if (!masterFlag[l]) continue;
    i = ptr[l][0];
    j = ptr[l][1];
    if(LSS.withCracking()) {
      iActive = !LSS.isOccluded(0,i);
      jActive = !LSS.isOccluded(0,j);
    } else {
      iActive = LSS.isActive(0,i);
      jActive = LSS.isActive(0,j);
    }
    intersect = LSS.edgeIntersectsStructure(0,l);

    if (!iActive && !jActive) continue; //both inside structure
    if (!intersect) continue; 
    //if (iActive && jActive && !intersect) continue;  //KW: switched to the above line on Nov.30,2010. 

    // now (i,j) must intersect the structure.
    if (iActive) {
      Vec3D flocal(0.0,0.0,0.0);
      LevelSetResult lsRes = LSS.getLevelSetDataAtEdgeCenter(0.0,l,true);
      flocal = (pstarij[l]-pInfty)*normal[l];
      sendLocalForce(flocal, lsRes, Fs);
    }
    if (jActive) {
      Vec3D flocal(0.0,0.0,0.0);
      LevelSetResult lsRes = LSS.getLevelSetDataAtEdgeCenter(0.0,l,false);
      flocal = -(pstarji[l]-pInfty)*normal[l];
      sendLocalForce(flocal, lsRes, Fs);
    }
  }
}
*/
//-----------------------------------------------------------------------------------------------
/*
void SubDomain::computeRecSurfBasedForceLoad(int forceApp, int order, SVec<double,3> &X,
                                             double (*Fs)[3], int sizeFs, LevelSetStructure &LSS,
                                             Vec<double> &pstarij, Vec<double> &pstarji, double pInfty)
{
  int T[4]; //nodes in a tet.
  double x[4][3]; //coords of nodes in a tet.
  int nPos, nNeg;
  int (*ptr)[2];
  ptr = edges.getPtr();

  int CODE;
  if (forceApp==3 && orderOfAccuracy==1)  CODE = 1;
  else if (forceApp==3 && orderOfAccuracy==2)  CODE = 3;
  else {fprintf(stderr,"ERROR: force method not recognized! Abort...\n"); exit(-1);}

  // for debugging: output the reconstructed surface.
  char nodeFileName[1024];std::sprintf(nodeFileName,"firstLayer/nodes_%d.top",globSubNum);
  char elemFileName[1024];std::sprintf(elemFileName,"firstLayer/elems_%d.top",globSubNum);
  char nodesName[1024];std::sprintf(nodesName,"InsideNodes%d",globSubNum);
  char elemsName[1024];std::sprintf(elemsName,"InsideElems%d",globSubNum);
#if 0 // Debug 
  FILE* nodeFile = fopen(nodeFileName,"w");
  FILE* elemFile = fopen(elemFileName,"w");
#else
  FILE* nodeFile = 0;
  FILE* elemFile = 0;
#endif
  int nodeCount = 1;
  if(nodeFile) fprintf(nodeFile, "Nodes %s\n",nodesName);

  int elemCount = 1;
  if(elemFile) fprintf(elemFile, "Elements %s using %s\n",elemsName,nodesName);
  // -----------------------------------------------

  for (int iElem=0; iElem<elems.size(); iElem++) {
    for (int i=0; i<4; i++) T[i] = elems[iElem][i];
    for (int i=0; i<4; i++) {
        x[i][0] = X[T[i]][0];  x[i][1] = X[T[i]][1];  x[i][2] = X[T[i]][2];}
    PolygonReconstructionData polygons[4];
    Vec3D nf;
    int numberOfPolygons = getPolygons(iElem, LSS, polygons); assert(numberOfPolygons<=4);
    int index,vertex1,vertex2;
    for(int i=0; i < numberOfPolygons; ++i){
        PolygonReconstructionData& polygon=polygons[i];
        LevelSetResult lsRes[4];
        if(!polygon.numberOfEdges) continue;
        for (int k=0; k<polygon.numberOfEdges; ++k) {
            lsRes[k] = LSS.getLevelSetDataAtEdgeCenter(0,polygon.edge[k],polygon.edgeWithVertex[k][0]<polygon.edgeWithVertex[k][1]);
            if (lsRes[k].alpha<0){fprintf(stderr,"Unable to get intersection results at edge center! Abort...\n"); exit(-1);}}

        Vec3D Xinter[4];
        Vec3D start_vertex;
        double pStar[4];
        double dist13,dist02;
        for(int k=0; k<polygon.numberOfEdges; ++k) {
            vertex1 = polygon.edgeWithVertex[k][0];
            vertex2 = polygon.edgeWithVertex[k][1];
            if(k==0){start_vertex[0]=X[vertex1][0];start_vertex[1]=X[vertex1][1];start_vertex[2]=X[vertex1][2];}
            int l = edges.find(polygon.edgeWithVertex[k][0], polygon.edgeWithVertex[k][1]);
            double alpha = lsRes[k].alpha;
            Xinter[k][0] = alpha*X[vertex1][0] + (1.0-alpha)*X[vertex2][0];
            Xinter[k][1] = alpha*X[vertex1][1] + (1.0-alpha)*X[vertex2][1];
            Xinter[k][2] = alpha*X[vertex1][2] + (1.0-alpha)*X[vertex2][2];

            if((LSS.withCracking()   && LSS.isOccluded(0,vertex1)) ||
               (!LSS.withCracking()) && LSS.isActive(0,vertex1)) {
              pStar[k] = (vertex1<vertex2) ? pstarij[l] : pstarji[l];}
            else{pStar[k] = pInfty;}}

        switch(polygon.numberOfEdges){
#if 0
        // TODO(jontg): Verify that these reconstructions are necessary
        case 1:
              index=1;
              for(int search_node=0; search_node < 4; ++search_node)
                  if(vertex1 != T[search_node] && vertex2 != T[search_node]){
                      Xinter[index][0]   = 0.333333333333333*(X[vertex1][0]+X[vertex2][0]+x[search_node][0]);
                      Xinter[index][1]   = 0.333333333333333*(X[vertex1][1]+X[vertex2][1]+x[search_node][1]);
                      Xinter[index++][2] = 0.333333333333333*(X[vertex1][2]+X[vertex2][2]+x[search_node][2]);}
              pStar[1] = pStar[2] = pStar[0];
              lsRes[1].alpha=-10.0; lsRes[2].alpha=-10.0;
              nf = (0.5*(Xinter[1]-Xinter[0])^(Xinter[2]-Xinter[0]));
              if(nf*(Xinter[1]-start_vertex) >= 0) nf *=-1;
              addLocalForce(CODE,nf,pStar[0],pStar[1],pStar[2],lsRes[0],lsRes[1],lsRes[2],Fs);
              break;
        case 2:
              Xinter[2][0] = 0.25*(x[0][0]+x[1][0]+x[2][0]+x[3][0]);
              Xinter[2][1] = 0.25*(x[0][1]+x[1][1]+x[2][1]+x[3][1]);
              Xinter[2][2] = 0.25*(x[0][2]+x[1][2]+x[2][2]+x[3][2]);
              pStar[2]=0.5*(pStar[0]+pStar[1]);
              lsRes[2].alpha=-10.0;
              nf = (0.5*(Xinter[1]-Xinter[0])^(Xinter[2]-Xinter[0]));
              if(nf*(Xinter[1]-start_vertex) >= 0) nf *=-1;
              addLocalForce(CODE,nf,pStar[0],pStar[1],pStar[2],lsRes[0],lsRes[1],lsRes[2],Fs);
              break;
#endif
        case 3:
              nf = 0.5*(Xinter[1]-Xinter[0])^(Xinter[2]-Xinter[0]);
              if(nf*(Xinter[1]-start_vertex) >= 0) nf *=-1;
              addLocalForce(CODE,nf,pStar[0],pStar[1],pStar[2],lsRes[0],lsRes[1],lsRes[2],Fs);
              break;
        case 4:
              dist02 = (Xinter[2]-Xinter[0]).norm();
              dist13 = (Xinter[3]-Xinter[1]).norm();
              if(dist02<dist13){ // connect 0,2.
                  nf = 0.5*(Xinter[1]-Xinter[0])^(Xinter[2]-Xinter[0]);
                  if(nf*(Xinter[1]-start_vertex) >= 0) nf *=-1;
                  addLocalForce(CODE,nf, pStar[0], pStar[1], pStar[2], lsRes[0], lsRes[1], lsRes[2], Fs);
                  nf = 0.5*(Xinter[2]-Xinter[0])^(Xinter[3]-Xinter[0]);
                  if(nf*(Xinter[1]-start_vertex) >= 0) nf *=-1;
                  addLocalForce(CODE,nf, pStar[0], pStar[2], pStar[3], lsRes[0], lsRes[2], lsRes[3], Fs);}
              else{
                  nf = 0.5*(Xinter[2]-Xinter[1])^(Xinter[3]-Xinter[1]);
                  if(nf*(Xinter[1]-start_vertex) >= 0) nf *=-1;
                  addLocalForce(CODE,nf, pStar[1], pStar[2], pStar[3], lsRes[1], lsRes[2], lsRes[3], Fs);
                  nf = 0.5*(Xinter[1]-Xinter[0])^(Xinter[3]-Xinter[0]);
                  if(nf*(Xinter[1]-start_vertex) >= 0) nf *=-1;
                  addLocalForce(CODE,nf, pStar[0], pStar[1], pStar[3], lsRes[0], lsRes[1], lsRes[3], Fs);}
              break;
#if 0
        default: fprintf(stderr,"ANOTHER PROBLEM!!!\n");break;
#endif
        default: break;
        }
        if(nodeFile&&elemFile) {
            fprintf(stderr,"Printing FirstLayers to Disc \"%s\" and \"%s\" %d\n",nodeFileName,elemFileName,polygon.numberOfEdges);
            fprintf(stderr,"%d %e %e %e\n", nodeCount, Xinter[0][0], Xinter[0][1], Xinter[0][2]); 
            fprintf(nodeFile,"%d %e %e %e\n", nodeCount++, Xinter[0][0], Xinter[0][1], Xinter[0][2]); 
            fprintf(stderr,"%d %e %e %e\n", nodeCount, Xinter[1][0], Xinter[1][1], Xinter[1][2]); 
            fprintf(nodeFile,"%d %e %e %e\n", nodeCount++, Xinter[1][0], Xinter[1][1], Xinter[1][2]); 
            fprintf(stderr,"%d %e %e %e\n", nodeCount, Xinter[2][0], Xinter[2][1], Xinter[2][2]); 
            fprintf(nodeFile,"%d %e %e %e\n", nodeCount++, Xinter[2][0], Xinter[2][1], Xinter[2][2]); 
            if(polygon.numberOfEdges<=3){
                fprintf(elemFile,"%d %d %d %d %d\n", elemCount++, 4, nodeCount-3, nodeCount-2, nodeCount-1);}
            else{
                fprintf(nodeFile,"%d %e %e %e\n", nodeCount++, Xinter[3][0], Xinter[3][1], Xinter[3][2]); 
                if(dist02<dist13){
                    fprintf(elemFile,"%d %d %d %d %d\n", elemCount++, 4, nodeCount-4, nodeCount-3, nodeCount-2);
                    fprintf(elemFile,"%d %d %d %d %d\n", elemCount++, 4, nodeCount-4, nodeCount-2, nodeCount-1);}
                else{
                    fprintf(elemFile,"%d %d %d %d %d\n", elemCount++, 4, nodeCount-3, nodeCount-2, nodeCount-1);
                    fprintf(elemFile,"%d %d %d %d %d\n", elemCount++, 4, nodeCount-4, nodeCount-3, nodeCount-1);}
            }
        }
    }
  }
  if(nodeFile) fclose(nodeFile);
  if(elemFile) fclose(elemFile);
}
*/

//--------------------------------------------------------------------------

void SubDomain::addLocalForce(int METHOD, Vec3D nf, double p1, double p2, double p3,
                              LevelSetResult& lsRes1, LevelSetResult& lsRes2, LevelSetResult& lsRes3,
                              double(*Fs)[3])
{
  Vec3D flocal;
  if (METHOD==1) {
    flocal = -1.0/3.0*p1*nf;
    sendLocalForce(flocal, lsRes1, Fs);
    flocal = -1.0/3.0*p2*nf;
    sendLocalForce(flocal, lsRes2, Fs);
    flocal = -1.0/3.0*p3*nf;
    sendLocalForce(flocal, lsRes3, Fs);
    return;
  }
  if (METHOD==2) {
    double c1 = 1.0/6.0, c2 = 1.0/12.0;
    flocal = -(c1*p1 + c2*(p2+p3))*nf;
    sendLocalForce(flocal, lsRes1, Fs);
    flocal = -(c1*p2 + c2*(p1+p3))*nf;
    sendLocalForce(flocal, lsRes2, Fs);
    flocal = -(c1*p3 + c2*(p1+p2))*nf;
    sendLocalForce(flocal, lsRes3, Fs);
    return;
  }
  fprintf(stderr,"METHOD=%d doesn't exist! Skipping the force calculation.\n", METHOD);
}

//-----------------------------------------------------------------------------------------------

void SubDomain::sendLocalForce(Vec3D flocal, LevelSetResult& lsRes, double(*Fs)[3])
{
  if(lsRes.alpha < 0) return;

	for(LevelSetResult::iterator it = lsRes.begin(); it != lsRes.end(); ++it) 
	{
    int n = it.nodeNum();
    double coef = it.Ni();

		for(int i=0; i<3; ++i) Fs[n][i] += coef*flocal[i];
	}

}

//-----------------------------------------------------------------------------------------------

int SubDomain::findFarfieldNode()
{
  for (int i=0; i<faces.size(); i++) {
    Face &curFace = faces[i];
    if (curFace.getCode() == BC_OUTLET_MOVING || curFace.getCode() == BC_OUTLET_FIXED || 
        curFace.getCode() == BC_INLET_MOVING || curFace.getCode() == BC_INLET_FIXED   ||
        curFace.getCode() == BC_DIRECTSTATE_OUTLET_MOVING || curFace.getCode() == BC_DIRECTSTATE_OUTLET_FIXED || 
        curFace.getCode() == BC_DIRECTSTATE_INLET_MOVING || curFace.getCode() == BC_DIRECTSTATE_INLET_FIXED ||
        curFace.getCode() == BC_MASSFLOW_OUTLET_MOVING || curFace.getCode() == BC_MASSFLOW_OUTLET_FIXED || 
        curFace.getCode() == BC_MASSFLOW_INLET_MOVING || curFace.getCode() == BC_MASSFLOW_INLET_FIXED) 
      return curFace[0];
  }
  return -1;
}

//-----------------------------------------------------------------------------------------------

void SubDomain::setupFluidIdVolumesInitialConditions(const int volid, const int myId, Vec<int> &fluidId)
{
  for (int iElem = 0; iElem < elems.size(); iElem++)  {
    if (elems[iElem].getVolumeID() == volid)  {
      int *nodeNums = elems[iElem].nodeNum();
      for (int iNode = 0; iNode < elems[iElem].numNodes(); iNode++) {
        int myNode = nodeNums[iNode];
        if(fluidId[myNode]>0 && fluidId[myNode]!=myId)
          fprintf(stderr,"WARNING: FluidId on node %d converted from %d to %d using volume Id(s)!\n", locToGlobNodeMap[myNode]+1, fluidId[myNode], myId); 
        fluidId[myNode] = myId;
      }
    }
  }
}

void SubDomain::solicitFluidIdFS(LevelSetStructure &LSS, Vec<int> &fluidId, SVec<bool,4> &poll)
{
  /* poll[0,1,2,3]|   indication
     -------------+---------------
     1  0  0  0   |     Id = 0
     0  1  0  0   |     Id = 1
     0  0  1  0   |     Id = 2
     0  0  0  1   |     Id = occluded
     0  0  0  0   |     no info available
     1  1  1  1   |     can't decide */
  
  //if(LSS.numOfFluids()!=2) {fprintf(stderr,"ERROR: #Fluid must be 2! Now it is %d\n",LSS.numOfFluids());exit(-1);}
  const Connectivity &Node2Node = *getNodeToNode();
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  for(int i=0; i<nodes.size(); i++) {
    poll[i][0] = poll[i][1] = poll[i][2] = poll[i][3] = false;
    bool swept = LSS.isSwept(0.0,i);
    bool occluded = LSS.isOccluded(0.0,i);

/*    if (i == 4498)
      fprintf(stderr,"Rank %d, i = %d, swept = %d, occluded = %d, fluidId = %d, dimLS = %d, poll = %d %d %d %d. LSS.numOfFluids = %d.\n", rank,i, swept, occluded, fluidId[i], 0xFFFF, poll[i][0], poll[i][1], poll[i][2], poll[i][3], LSS.numOfFluids());
*/

/*
    if(locToGlobNodeMap[i]+1==349859) {
      fprintf(stderr,"rank = %d, i = %d.\n", rank, i);      
    }
*/
    if(!swept){ //fluidId should not change.
      if(!occluded && fluidId[i]!=LSS.numOfFluids()) {//this "if" is false when the structural elment covering node i got deleted in Element Deletion.
        poll[i][fluidId[i]] = true;
        continue;
      }
    }

    if(occluded){ //set Id to 2
      poll[i][3] = true;
      continue;
    }

    int caught = 0;
    if(locToGlobNodeMap[i]+1==349859)
      caught = 1;

    int myId = -1;
    int iNei;
    bool consistent = false;
    for(int j=0; j<Node2Node.num(i); j++) {
      iNei = Node2Node[i][j];
      if(i==iNei)
        continue;

      if(caught)
        fprintf(stderr,"Sub %d, Nei of 349859 (myId = %d): %d, fluidId = %d, occluded = %d, swept = %d, X = %d\n",
                globSubNum, myId, locToGlobNodeMap[iNei]+1, fluidId[iNei], LSS.isOccluded(0.0,iNei), LSS.isSwept(0.0,iNei), (int)LSS.edgeIntersectsStructure(0.0,edges.findOnly(i,iNei)));

      if(LSS.isOccluded(0.0,iNei) || LSS.isSwept(0.0,iNei) || LSS.edgeIntersectsStructure(0.0,edges.findOnly(i,iNei)))
        continue;


//----------------------------------------
/*      if(myId==-1) {
        myId = fluidId[iNei];
        consistent = true;
      } else if(myId!=fluidId[iNei]) {
        consistent = false;
        break;
      }
*/
      poll[i][fluidId[iNei]] = true;

//------------------------------------------
    }
/*
    if(consistent)
      poll[i][myId] = true; //its visible neighbors have the same id
    else
      poll[i][0] = poll[i][1] = poll[i][2] = (myId!=-1); //either 'no info' or 'can't decide'.
*/


    if(caught)
      fprintf(stderr,"Sub %d, 349859, myId = %d, consistent = %d, poll = %d %d %d %d\n",
              globSubNum, myId, consistent, poll[i][0], poll[i][1], poll[i][2], poll[i][3]);

  }
}

//------------------------------------------------------------------------------

void SubDomain::computeConnectedTopology(const std::vector<int> &locSampleNodes_, const std::vector<int> &globalNeighborNodes) 
{

	sampleMesh = true;
	locSampleNodes = locSampleNodes_;
	numSampledNodes = locSampleNodes.size();
	elems.computeConnectedElems(locSampleNodes);
	edges.computeConnectedEdges(locSampleNodes);
	edges.computeGlobalConnectedEdges(globalNeighborNodes,locToGlobNodeMap);
	faces.computeConnectedFaces(locSampleNodes);

}

void SubDomain::createHigherOrderMultiFluid() {

  higherOrderMF = new HigherOrderMultiFluid();

  edges.attachHigherOrderMultiFluid(higherOrderMF);
  faces.attachHigherOrderMF(higherOrderMF);
}

void SubDomain::createHigherOrderFSI() {

  higherOrderFSI = new HigherOrderFSI();

  edges.attachHigherOrderFSI(higherOrderFSI);
  //  faces.attachHigherOrderFSI(higherOrderFSI);
}

void SubDomain::assignErrorHandler(ErrorHandler* in){
  errorHandler = in;
  edges.assignErrorHandler(errorHandler);
}

void SubDomain::getSurfaceNodes(Aerof_unordered_set<int>::type& boundaryNodes) const {

  for (int i = 0; i < faces.size(); ++i) {

    int code = faces[i].getCode();
    if (code != BC_SLIP_WALL_MOVING && code != BC_SLIP_WALL_FIXED && 
        code != BC_ADIABATIC_WALL_MOVING && code != BC_ADIABATIC_WALL_FIXED && 
        code != BC_ISOTHERMAL_WALL_MOVING && code != BC_ISOTHERMAL_WALL_FIXED &&
        code != BC_POROUS_WALL_MOVING && code != BC_POROUS_WALL_FIXED)
      continue;

    for (int k = 0; k < faces[i].numNodes(); ++k)
      boundaryNodes.insert(faces[i][k]); 
    
  }
}


void SubDomain::getSolidBoundaryNodes(Aerof_unordered_set<int>::type& boundaryNodes) const {

  for (int i = 0; i < faces.size(); ++i) {

    int code = faces[i].getCode();
    if (code != BC_SYMMETRY && 
        code != BC_SLIP_WALL_MOVING && code != BC_SLIP_WALL_FIXED && 
        code != BC_ADIABATIC_WALL_MOVING && code != BC_ADIABATIC_WALL_FIXED && 
        code != BC_ISOTHERMAL_WALL_MOVING && code != BC_ISOTHERMAL_WALL_FIXED &&
        code != BC_POROUS_WALL_MOVING && code != BC_POROUS_WALL_FIXED)
      continue;

    for (int k = 0; k < faces[i].numNodes(); ++k)
      boundaryNodes.insert(faces[i][k]); 
    
  }
}

void SubDomain::getFarFieldBoundaryNodes(Aerof_unordered_set<int>::type& boundaryNodes) const {

  for (int i = 0; i < faces.size(); ++i) {

    int code = faces[i].getCode();
    if (code != BC_OUTLET_MOVING && code != BC_OUTLET_FIXED &&
        code != BC_INLET_MOVING && code != BC_INLET_FIXED &&
        code != BC_DIRECTSTATE_OUTLET_MOVING && code != BC_DIRECTSTATE_OUTLET_FIXED &&
        code != BC_DIRECTSTATE_INLET_MOVING && code != BC_DIRECTSTATE_INLET_FIXED &&
        code != BC_MASSFLOW_OUTLET_MOVING && code != BC_MASSFLOW_OUTLET_FIXED &&
        code != BC_MASSFLOW_INLET_MOVING && code != BC_MASSFLOW_INLET_FIXED)  
      continue;

    for (int k = 0; k < faces[i].numNodes(); ++k)
      boundaryNodes.insert(faces[i][k]); 
    
  }
}

void SubDomain::getSubDomainBoundaryNodes(Aerof_unordered_set<int>::type& boundaryNodes) const {

  for(int iSub=0; iSub<numNeighb; iSub++) {
    for(int i=0; i<sharedNodes->num(iSub); i++) {
      int me = (*sharedNodes)[iSub][i];
      boundaryNodes.insert(me);
    }
  }
}

// This implementation is currently incomplete (AM)
void SubDomain::constructLines(std::vector<std::vector<int>*>& pLines, int& numLines) {

  Aerof_unordered_set<int>::type surfaceNodes;
  Aerof_unordered_set<int>::type subdBoundaryNodes;

  getSurfaceNodes(surfaceNodes);
  getSubDomainBoundaryNodes(subdBoundaryNodes);

  Vec<int> available(numNodes());
  available = 1;

  for (int i = 0; i < numNodes(); ++i) {

  }
  
}

//------------------------------------------------------------------------------

void SubDomain::setFarFieldNodes()
{

  farFieldNodes.clear();

  int i, j;

  //find numFarFieldNodes and tag the nodes that are at inlet/outlet
  Vec<bool> tagNodes(nodes.size());
  tagNodes = false;
  int numFarFieldNodes = 0;
  for ( i=0; i<nodes.size(); i++){
    if ( nodeType[i] == BC_OUTLET_MOVING || nodeType[i] == BC_INLET_MOVING ||
         nodeType[i] == BC_OUTLET_FIXED  || nodeType[i] == BC_INLET_FIXED ){
      farFieldNodes.push_back(i);
      numFarFieldNodes++;
    }
  }

  //check that we get the same number of far field nodes counting via nodes, faces, and tets
  //via tets
  /*
  Vec<bool> tcounted(nodes.size());
  tcounted = false;
  int tnumFarFieldNodes=0;
  int n;
  for( i=0; i<elems.size(); i++){
    for( j=0; j<4; j++){
      n=elems[i][j];
      if (!tcounted[n] && (nodeType[n] ==BC_OUTLET_MOVING || nodeType[n] ==BC_INLET_MOVING ||
			   nodeType[n] ==BC_OUTLET_FIXED ||nodeType[n] ==BC_INLET_FIXED )){
        tcounted[n] = true;
        tnumFarFieldNodes++;
      }
    }
  }

  //via faces
  Vec<bool> fcounted(nodes.size());
  fcounted = false;
  int fnumFarFieldNodes=0;
  for( i=0; i<faces.size(); i++){
    for( j=0; j<3; j++){
      n=faces[i][j];
      if (!fcounted[n] && (nodeType[n] ==BC_OUTLET_MOVING || nodeType[n] ==BC_INLET_MOVING ||
			   nodeType[n] ==BC_OUTLET_FIXED ||nodeType[n] ==BC_INLET_FIXED )){
	fcounted[n] = true;
	fnumFarFieldNodes++;
      }
    }
  }

  assert((numFarFieldNodes==fnumFarFieldNodes)&&(numFarFieldNodes==tnumFarFieldNodes)); 
  */
}

//------------------------------------------------------------------------------

void SubDomain::setWallNodes()
{

  wallNodes.clear();

  int i, j;

  //find numWallNodes and tag the nodes that are at inlet/outlet
  Vec<bool> tagNodes(nodes.size());
  tagNodes = false;
  int numWallNodes = 0;
  for ( i=0; i<nodes.size(); i++){
    if ( nodeType[i] == BC_ADIABATIC_WALL_MOVING || nodeType[i] == BC_ADIABATIC_WALL_FIXED ||
         nodeType[i] == BC_SLIP_WALL_MOVING || nodeType[i] == BC_SLIP_WALL_FIXED ||
         nodeType[i] == BC_ISOTHERMAL_WALL_MOVING || nodeType[i] == BC_ISOTHERMAL_WALL_FIXED ){
      wallNodes.push_back(i);
      numWallNodes++;
    }
  }

  //check that we get the same number of far field nodes counting via nodes, faces, and tets
  //via tets
  /*
  Vec<bool> tcounted(nodes.size());
  tcounted = false;
  int tnumWallNodes=0;
  int n;
  for( i=0; i<elems.size(); i++){
    for( j=0; j<4; j++){
      n=elems[i][j];
      if (!tcounted[n] && (nodeType[n] == BC_ADIABATIC_WALL_MOVING || nodeType[n] == BC_ADIABATIC_WALL_FIXED ||
         nodeType[n] == BC_SLIP_WALL_MOVING || nodeType[n] == BC_SLIP_WALL_FIXED ||
         nodeType[n] == BC_ISOTHERMAL_WALL_MOVING || nodeType[n] == BC_ISOTHERMAL_WALL_FIXED )){
        tcounted[n] = true;
        tnumWallNodes++;
      }
    }
  }

  //via faces
  Vec<bool> fcounted(nodes.size());
  fcounted = false;
  int fnumWallNodes=0;
  for( i=0; i<faces.size(); i++){
    for( j=0; j<3; j++){
      n=faces[i][j];
      if (!fcounted[n] && (nodeType[n] == BC_ADIABATIC_WALL_MOVING || nodeType[n] == BC_ADIABATIC_WALL_FIXED ||
         nodeType[n] == BC_SLIP_WALL_MOVING || nodeType[n] == BC_SLIP_WALL_FIXED ||
         nodeType[n] == BC_ISOTHERMAL_WALL_MOVING || nodeType[n] == BC_ISOTHERMAL_WALL_FIXED )){
	      fcounted[n] = true;
	      fnumWallNodes++;
      }
    }
  }

  assert((numWallNodes==fnumWallNodes)&&(numWallNodes==tnumWallNodes));  */

}

//------------------------------------------------------------------------------

void SubDomain::setFarFieldMask(Vec<double>& ffMask, Vec<double>& neighborMask){

  // Note: vectors are initialized to zero in DomainCore.C

  const Connectivity &Node2Node = *getNodeToNode();

  for (int iNode=0; iNode<farFieldNodes.size(); iNode++){
    int ffNode = farFieldNodes[iNode];
    ffMask[ffNode] = 1.0;
    for(int jNode=0; jNode<Node2Node.num(ffNode); jNode++) {
      int neighbor = Node2Node[ffNode][jNode];
      if(neighbor==ffNode) {
        continue;
      } else {
        neighborMask[neighbor] = 1.0;
      }
    }
  }

}

//------------------------------------------------------------------------------

void SubDomain::setWallMask(Vec<double>& wallMask, Vec<double>& neighborMask){

  // Note: vectors are initialized to zero in DomainCore.C

  const Connectivity &Node2Node = *getNodeToNode();

  for (int iNode=0; iNode<wallNodes.size(); iNode++){
    int wallNode = wallNodes[iNode];
    wallMask[wallNode] = 1.0;
    for(int jNode=0; jNode<Node2Node.num(wallNode); jNode++) {
      int neighbor = Node2Node[wallNode][jNode];
      if(neighbor==wallNode) {
        continue;
      } else {
        neighborMask[neighbor] = 1.0;
      }
    }
  }

}


//------------------------------------------------------------------------------


void SubDomain::maskHHVector(Vec<double>& hh) {

  for (int i = 0; i < faces.size(); ++i) {

    switch (faces[i].getCode()) {

    case BC_OUTLET_MOVING:
    case BC_OUTLET_FIXED:
    case BC_INLET_MOVING:
    case BC_INLET_FIXED:
      break;
    default:
      hh[i] = 0;
      break;
    }
  }
}

void SubDomain::attachTriangulatedInterfaceLSS(LevelSetStructure* LSS) {

  edges.attachTriangulatedInterfaceLSS(LSS);
}

class ElemSearchValid {

  public:

    bool Valid(Elem*) { return true; }
};


Elem* SubDomain::searchPoint(Vec3D Xp, SVec<double,3>& X) {

  ElemSearchValid myObj;
  Elem* E = myTree->search<&Elem::isPointInside, ElemSearchValid,
 	                   &ElemSearchValid::Valid>(&myObj, X, Xp);

  return E;
}


//-----------------------------------------------------------------------------

void SubDomain::printDistVecBool(Vec<bool> &X, int numSub, bool status)
{
  int glob;
  for (int i=0; i<nodes.size(); i++){
    glob = locToGlobNodeMap[i]+1;
    if(X[i] == status) fprintf(stdout, "X[%d] is %d. global node is %d\n", i, X[i], glob);
  }
}

bool SubDomain::getPiercedPoint(Vec3D va, Vec3D vb, Vec3D vc, Vec3D vd, Vec3D ve, double &r, double &t)
{ 

	r = 0.0; 
	t = 0.0;

	int pierce = checkPiercePoint(va, vb, vc, vd, ve, true);

	if(pierce >= -3 && pierce <= 3) 
	{
		constructPiercePoint(va, vb, vc, vd, ve, 0, &r, &t, true); 
		return true;
	}
	else 
		return false;

}

// ----------------------------------------

bool SubDomain::getSIstencil(int Ni, int Nj, SVec<double,3> &X, LevelSetStructure &LSS, Vec<int> &fluidId,
									  Vec3D &normWall, Vec3D &xWall, V6NodeData &SiStencilData)
{

	double min_dist0 = bigNum;
	double min_dist1 = bigNum;

	bool Wmode1 = false;
	bool Wmode0 = false;

	int ele1, ele0, face1, face0;
	ele1 = ele0 = face1 = face0 = -1;

	double r1, r0, t1, t0;
	r1 = r0 = t1 = t0 = 0.0;

	Vec3D X_si; 
	for(int k=0; k<3; ++k) X_si[k] = 0.5*(X[Ni][k] + X[Nj][k]);

	Vec3D dir = X_si - xWall;
	double norm = sqrt(dir * dir);
	if(norm != 0.0) dir *= 1.0 / norm;

	Vec3D ve = X_si + 1000.0*dir;

	int N_act = LSS.isActive(0.0, Ni) ? Ni : Nj;

	// Loop over the elements connected to nodes Ni and Nj
	for(int idN=0; idN<2; ++idN)
	{
		int Nnode = (idN == 0) ? Ni : Nj;

		int Esize;
		int* Elist;
		Elist = getNeiElemOfNode(Nnode, 1, Esize);

		for(int n=0; n<Esize; ++n)
		{	
			int nE = Elist[n];

			// Loop over the faces of the current element
			for(int tf=0; tf<4; ++tf)
			{				
				int n0_ = elems[nE].faceDef(tf,0); int N0 = elems[nE].nodeNum(n0_);
				int n1_ = elems[nE].faceDef(tf,1); int N1 = elems[nE].nodeNum(n1_);
				int n2_ = elems[nE].faceDef(tf,2); int N2 = elems[nE].nodeNum(n2_);
				
				// Coordinates of the 3 vertices of the current face
				Vec3D va(X[N0]); Vec3D vb(X[N1]); Vec3D vc(X[N2]);

				double r=0, t=0;
				bool newPoint = getPiercedPoint(va, vb, vc, X_si, ve, r, t);

				int mode = 0;
				if(LSS.isActive(0.0, N0) && 
					LSS.isActive(0.0, N1) && 
					LSS.isActive(0.0, N2) &&				
					fluidId[N0] == fluidId[N_act] &&	
					fluidId[N1] == fluidId[N_act] && 
					fluidId[N2] == fluidId[N_act]) mode = 1;

				if(newPoint && mode == 1)
				{
					/* 
						--- mode1 ---
						the line (Xwall - Xsi) pierces a face which
						has only active nodes; Xp is the true pierced point
					*/	

					Vec3D Xp;
					for(int k=0; k<3; ++k)
				 		Xp[k] = X[N2][k] + r*(X[N0][k] - X[N2][k])+
						                   t*(X[N1][k] - X[N2][k]);

					Vec3D vDist = Xp - X_si;
					double dist = sqrt(vDist * vDist);
					
					if(dist < min_dist1)
					{
						min_dist1 = dist;

						ele1 = nE; face1 = tf;

						r1 = r; t1 = t;

						Wmode1 = true;
					}
				}
				else if(newPoint && mode == 0)
				{
					/* 
						--- mode0 ---
						the line (Xwall - Xsi) pierces a face which
						has some inactive nodes; Xp is the active node
						closest to Xsi
					*/	
					
					Vec3D vDista = va - X_si; double dista = sqrt(vDista * vDista);
					Vec3D vDistb = vb - X_si; double distb = sqrt(vDistb * vDistb);
					Vec3D vDistc = vc - X_si; double distc = sqrt(vDistc * vDistc);

					// To exclude nodes with a different fluidId
					if(fluidId[N0] != fluidId[N_act]) dista = bigNum;
					if(fluidId[N1] != fluidId[N_act]) distb = bigNum;
					if(fluidId[N2] != fluidId[N_act]) distc = bigNum;
               //... or inactive nodes
					if(!LSS.isActive(0.0, N0)) dista = bigNum;
					if(!LSS.isActive(0.0, N1)) distb = bigNum;
					if(!LSS.isActive(0.0, N2)) distc = bigNum;

					double dist = std::min(dista, std::min(distb, distc));
					
					if(dist < min_dist0)
					{						
						min_dist0 = dist;

						ele0 = nE; face0 = tf;

						if(     fabs(dist - dista) < geomTol) {r0 = 1.0; t0 = 0.0;}
						else if(fabs(dist - distb) < geomTol) {r0 = 0.0; t0 = 1.0;}
						else                                  {r0 = 0.0; t0 = 0.0;}

						Wmode0 = true;
					}
				}
			}
		}
	}

	if(Wmode1)
	{
		SiStencilData.tet  = ele1;
		SiStencilData.face = face1;
		SiStencilData.r    = r1;
		SiStencilData.t    = t1;
	}
	else if(Wmode0 && !Wmode1) 
	{
		SiStencilData.tet  = ele0;
		SiStencilData.face = face0;
		SiStencilData.r    = r0;
		SiStencilData.t    = t0;
	}
	else
	{
		SiStencilData.tet  = -1;
		SiStencilData.face = -1;
		SiStencilData.r    = -1.0;
		SiStencilData.t    = -1.0;
	}

	return Wmode0 || Wmode1;

}

// ----------------------------------------

bool SubDomain::getFEMstencil(int Ni, SVec<double,3> &X, 
										LevelSetStructure &LSS, Vec<int> &fluidId,
										Vec3D &normWall, Vec3D &xWall, 
										V6NodeData &NodeStencilData)
{

	double min_dist0 = bigNum;
	double min_dist1 = bigNum;

	bool Wmode1 = false;
	bool Wmode0 = false;

	int ele1, ele0, face1, face0;
	ele1 = ele0 = face1 = face0 = -1;

	double r1, r0, t1, t0;
	r1 = r0 = t1 = t0 = 0.0;

	Vec3D X_i = X[Ni]; 

	Vec3D dir;
	if(LSS.isOccluded(0.0, Ni))
		dir = xWall - X_i;
	else
		dir = X_i - xWall;
		
	double norm = sqrt(dir * dir);
	if(norm != 0.0) dir *= 1.0 / norm;

	Vec3D ve = X_i + 1000.0*dir;

	int Esize;
	int* Elist;
	Elist = getNeiElemOfNode(Ni, 2, Esize);

	for(int n=0; n<Esize; ++n) 
	{	
		int nE = Elist[n];

		// Loop over the faces of the current element
		for(int tf=0; tf<4; ++tf)
		{				
			int n0_ = elems[nE].faceDef(tf,0); int N0 = elems[nE].nodeNum(n0_);
			int n1_ = elems[nE].faceDef(tf,1); int N1 = elems[nE].nodeNum(n1_);
			int n2_ = elems[nE].faceDef(tf,2); int N2 = elems[nE].nodeNum(n2_);
				
			// Coordinates of the 3 vertices of the current face
			Vec3D va(X[N0]); Vec3D vb(X[N1]); Vec3D vc(X[N2]);

			double r=0, t=0;
			bool newPoint = getPiercedPoint(va, vb, vc, X_i, ve, r, t);

			int mode = 0;
			if(LSS.isActive(0.0, N0) && 
				LSS.isActive(0.0, N1) && 
				LSS.isActive(0.0, N2) &&				
				fluidId[N0] == fluidId[Ni] &&	
				fluidId[N1] == fluidId[Ni] && 
				fluidId[N2] == fluidId[Ni]) mode = 1;

				if(newPoint && mode == 1)
				{
					/* 
						--- mode1 ---
						the line (Xwall - Xi) pierces a face which
						has only active nodes; Xp is the true pierced point
					*/	

					Vec3D Xp;
					for(int k=0; k<3; ++k)
				 		Xp[k] = X[N2][k] + r*(X[N0][k] - X[N2][k])+
						                   t*(X[N1][k] - X[N2][k]);

					Vec3D vDist = Xp - X_i;
					double dist = sqrt(vDist * vDist);
					
					if(dist < min_dist1)
					{
						min_dist1 = dist;

						ele1 = nE; face1 = tf;

						r1 = r; t1 = t;

						Wmode1 = true;
					}
				}
				else if(newPoint && mode == 0)
				{
					/* 
						--- mode0 ---
						the line (Xwall - Xi) pierces a face which
						has some inactive nodes; Xp is the active node
						closest to Xi
					*/	
					
					Vec3D vDista = va - X_i; double dista = sqrt(vDista * vDista);
					Vec3D vDistb = vb - X_i; double distb = sqrt(vDistb * vDistb);
					Vec3D vDistc = vc - X_i; double distc = sqrt(vDistc * vDistc);

					// To exclude nodes with a different fluidId...
					if(fluidId[N0] != fluidId[Ni]) dista = bigNum;
					if(fluidId[N1] != fluidId[Ni]) distb = bigNum;
					if(fluidId[N2] != fluidId[Ni]) distc = bigNum;
               //... or inactive nodes
					if(!LSS.isActive(0.0, N0)) dista = bigNum;
					if(!LSS.isActive(0.0, N1)) distb = bigNum;
					if(!LSS.isActive(0.0, N2)) distc = bigNum;

					double dist = std::min(dista, std::min(distb, distc));
					
					if(dist < min_dist0)
					{						
						min_dist0 = dist;

						ele0 = nE; face0 = tf;

						if(     fabs(dist - dista) < geomTol) {r0 = 1.0; t0 = 0.0;}
						else if(fabs(dist - distb) < geomTol) {r0 = 0.0; t0 = 1.0;}
						else                                  {r0 = 0.0; t0 = 0.0;}

						Wmode0 = true;
					}
				}
			}
		}


	if(Wmode1)
	{
		NodeStencilData.tet  = ele1;
		NodeStencilData.face = face1;
		NodeStencilData.r    = r1;
		NodeStencilData.t    = t1;
	}
	else if(Wmode0 && !Wmode1) 
	{
		NodeStencilData.tet  = ele0;
		NodeStencilData.face = face0;
		NodeStencilData.r    = r0;
		NodeStencilData.t    = t0;
	}
	else
	{
		NodeStencilData.tet  = -1;
		NodeStencilData.face = -1;
		NodeStencilData.r    = -1.0;
		NodeStencilData.t    = -1.0;
	}

	return Wmode0 || Wmode1;

}

// ------------------------------------------------------------

bool SubDomain::getFEMstencil2(int Ni, SVec<double,3> &X, 
										LevelSetStructure &LSS, Vec<int> &fluidId,
										Vec3D &normWall, Vec3D &xWall, 
										V6NodeData &NodeStencilData_p, 
										V6NodeData &NodeStencilData_m)
{

	NodeStencilData_p.tet  = -1;   	NodeStencilData_m.tet  = -1;
	NodeStencilData_p.face = -1;	   NodeStencilData_m.face = -1;
	NodeStencilData_p.r    = -1.0;   NodeStencilData_m.r    = -1.0;
	NodeStencilData_p.t    = -1.0;   NodeStencilData_m.t    = -1.0;

	double min_dist0_p = bigNum, min_dist0_m = bigNum; 
	double min_dist1_p = bigNum, min_dist1_m = bigNum;

	bool Wmode1_p = false, Wmode1_m = false;
	bool Wmode0_p = false, Wmode0_m = false;

	int ele1_p =-1,  ele1_m=-1,  ele0_p=-1,  ele0_m=-1;
	int face1_p=-1, face1_m=-1, face0_p=-1, face0_m=-1;

	double r1_p=0, r1_m=0, r0_p=0, r0_m=0;
	double t1_p=0, t1_m=0, t0_p=0, t0_m=0;

	Vec3D X_i = X[Ni];

	Vec3D dir = X_i - xWall;
		
	double norm = sqrt(dir * dir);

	if(norm > geomTol) 
		dir *= 1.0 / norm;
	else
	{
		// X_i and xWall coincide: use nWall as searching direction
	   norm = sqrt(normWall * normWall);
		dir = normWall * (1.0/norm);
	}

	Vec3D ve_p = X_i + 1000.0*dir;
	Vec3D ve_m = X_i - 1000.0*dir;

	int Esize;
	int* Elist;
	Elist = getNeiElemOfNode(Ni, 2, Esize);
	
	bool newPoint_p, newPoint_m;

	double r, t, dist;

	for(int n=0; n<Esize; ++n) 
	{	
		int nE = Elist[n];
		
		// Loop over the faces of the current element
		for(int tf=0; tf<4; ++tf)
		{
			int n0_ = elems[nE].faceDef(tf,0); int N0 = elems[nE].nodeNum(n0_);
			int n1_ = elems[nE].faceDef(tf,1); int N1 = elems[nE].nodeNum(n1_);
			int n2_ = elems[nE].faceDef(tf,2); int N2 = elems[nE].nodeNum(n2_);
			
			// Coordinates of the 3 vertices of the current face
			Vec3D va(X[N0]); Vec3D vb(X[N1]); Vec3D vc(X[N2]);
			
			double r_p = 0, t_p = 0;
			newPoint_p = getPiercedPoint(va, vb, vc, X_i, ve_p, r_p, t_p);

			double r_m = 0, t_m = 0;
			newPoint_m = getPiercedPoint(va, vb, vc, X_i, ve_m, r_m, t_m);	

			if(!newPoint_p && !newPoint_m) continue;

			if(newPoint_p && newPoint_m) 
			{
				fprintf(stderr, "   **** Error: points +/- identical\n");
				exit(-1);
			}

			int mode = 0;
			if(LSS.isActive(0.0, N0) &&
				LSS.isActive(0.0, N1) &&
				LSS.isActive(0.0, N2) &&
				fluidId[N0] == fluidId[Ni] &&
				fluidId[N1] == fluidId[Ni] &&
				fluidId[N2] == fluidId[Ni]) mode = 1;
			
			if(mode == 1)
			{
				/*
				  --- mode1 ---
				  the line (Xwall - Xi) pierces a face which
				  has only active nodes; Xp is the true pierced point
				*/

				if(newPoint_p) {r = r_p; t = t_p;}
				else           {r = r_m; t = t_m;}

				Vec3D Xp;
				for(int k=0; k<3; ++k)
					Xp[k] = X[N2][k] + r*(X[N0][k] - X[N2][k])+
						                t*(X[N1][k] - X[N2][k]);
					
				Vec3D vDist = Xp - X_i;
				dist = sqrt(vDist * vDist);
			}
			else
			{
				/* 
					--- mode0 ---
					the line (Xwall - Xi) pierces a face which
					has some inactive nodes; Xp is the active node
					closest to Xi
				*/	
					
				Vec3D vDista = va - X_i; double dista = sqrt(vDista * vDista);
				Vec3D vDistb = vb - X_i; double distb = sqrt(vDistb * vDistb);
				Vec3D vDistc = vc - X_i; double distc = sqrt(vDistc * vDistc);

				// To exclude nodes with a different fluidId...
				if(fluidId[N0] != fluidId[Ni]) dista = bigNum;
				if(fluidId[N1] != fluidId[Ni]) distb = bigNum;
				if(fluidId[N2] != fluidId[Ni]) distc = bigNum;
				//... or inactive nodes
				if(!LSS.isActive(0.0, N0)) dista = bigNum;
				if(!LSS.isActive(0.0, N1)) distb = bigNum;
				if(!LSS.isActive(0.0, N2)) distc = bigNum;
				
				dist = std::min(dista, std::min(distb, distc));

				if(     fabs(dist - dista) < geomTol) {r = 1.0; t = 0.0;}
				else if(fabs(dist - distb) < geomTol) {r = 0.0; t = 1.0;}
				else                                  {r = 0.0; t = 0.0;}
			}
			
			if(newPoint_p)
			{
				if(mode == 1 && dist < min_dist1_p)
				{
					min_dist1_p = dist;

					ele1_p = nE; face1_p = tf; 
					r1_p = r; t1_p = t;

					Wmode1_p = true;
				}
				else if(mode == 0 && dist < min_dist0_p)
				{
					min_dist0_p = dist;

					ele0_p = nE; face0_p = tf; 
					r0_p = r; t0_p = t;

					Wmode0_p = true;
				}				
			}
			else if(newPoint_m)
			{
				if(mode == 1 && dist < min_dist1_m)
				{
					min_dist1_m = dist;

					ele1_m = nE; face1_m = tf; 
					r1_m = r; t1_m = t;

					Wmode1_m = true;
				}
				else if(mode == 0 && dist < min_dist0_m)
				{
					min_dist0_m = dist;

					ele0_m = nE; face0_m = tf; 
					r0_m = r; t0_m = t;

					Wmode0_m = true;
				}
			}
		}		
	}

	// ----------------------------- //

	if(Wmode1_p)
	{
		NodeStencilData_p.tet  = ele1_p;
		NodeStencilData_p.face = face1_p;
		NodeStencilData_p.r    = r1_p;
		NodeStencilData_p.t    = t1_p;		
	}
	else if(Wmode0_p && !Wmode1_p)
	{
		NodeStencilData_p.tet  = ele0_p;
		NodeStencilData_p.face = face0_p;
		NodeStencilData_p.r    = r0_p;
		NodeStencilData_p.t    = t0_p;
	}

	if(Wmode1_m)
	{
		NodeStencilData_m.tet  = ele1_m;
		NodeStencilData_m.face = face1_m;
		NodeStencilData_m.r    = r1_m;
		NodeStencilData_m.t    = t1_m;
	}
	else if(Wmode0_m && !Wmode1_m)
	{
		NodeStencilData_m.tet  = ele0_m;
		NodeStencilData_m.face = face0_m;
		NodeStencilData_m.r    = r0_m;
		NodeStencilData_m.t    = t0_m;
	}
 
	bool Wmode = (Wmode0_m || Wmode1_m) || (Wmode0_p || Wmode1_p);

	return Wmode;

}

void SubDomain::printPoint(int Ni, SVec<double,3> &X, Vec3D &xWall, V6NodeData & NSdata, 
									int dir, int mode)
{
	int    idxTet  = NSdata.tet;
	int    idxFace = NSdata.face;
	double face_r  = NSdata.r;
	double face_t  = NSdata.t;

	int n0_loc = elems[idxTet].faceDef(idxFace, 0);
	int n1_loc = elems[idxTet].faceDef(idxFace, 1);
	int n2_loc = elems[idxTet].faceDef(idxFace, 2);

	int n0 = elems[idxTet].nodeNum(n0_loc);
	int n1 = elems[idxTet].nodeNum(n1_loc);
	int n2 = elems[idxTet].nodeNum(n2_loc);

	Vec3D Xf;
	for(int k=0; k<3; ++k)
		Xf[k] = X[n2][k] + face_r*(X[n0][k] - X[n2][k])
			              + face_t*(X[n1][k] - X[n2][k]);

	fprintf(stdout, "%d, %d, %15.6e,%15.6e,%15.6e, %15.6e,%15.6e,%15.6e, %15.6e,%15.6e,%15.6e, ", mode, dir, 
			  xWall[0], xWall[1], xWall[2],
			  X[Ni][0],X[Ni][1],X[Ni][2], 
			  Xf[0],Xf[1],Xf[2]);

}
