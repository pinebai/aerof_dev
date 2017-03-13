#include "TriangulatedInterface.h"

#include <stdlib.h>

#include <IntersectorPhysBAM/IntersectorPhysBAM.h>
 
TriangulatedInterface::TriangulatedInterface() {

  abc = NULL;
  xyz = NULL;
  xyzdot = NULL;
  constrained  = NULL;
}

TriangulatedInterface::~TriangulatedInterface() {

  if (abc)
    delete [] abc;
  
  if (xyz)
    delete [] xyz;

  if (xyzdot)
    delete [] xyzdot;

  if (constrained)
    delete [] constrained;
}

int TriangulatedInterface::NumNodes() {

  return numNodes;
}

int TriangulatedInterface::NumElems() {

  return numElems;
}
  
int (* TriangulatedInterface::getIndices())[3] {

  return abc;
}

double* TriangulatedInterface::getVertexLocations() {

  return xyz;
}

void TriangulatedInterface::initializeAsSquare(int N) {

  numNodes = (N+1)*(N+1);
  numElems = N*N*2;
  abc = new int[numElems][3];
 
  xyz = new double[numNodes*3];
  xyzdot = new double[numNodes*3];
  constrained = new int[numNodes];
  memset(constrained,0,sizeof(int)*numNodes);

  for (int i = 0; i <= N; ++i)  {

    double y = (double)i / N;
    for (int j = 0; j <= N; ++j) {

      double z = (double)j / N;
      xyz[(i*(N+1)+j)*3] = 0.499;
      xyz[(i*(N+1)+j)*3+1] = y*0.2;
      xyz[(i*(N+1)+j)*3+2] = z*0.2;  
    }
  }

  for (int i = 0; i <= N; ++i)  {

    constrained[i*(N+1)] |= 2;
    constrained[i*(N+1)+N] |= 2;
  }

  for (int i = 0; i <= N; ++i)  {

    constrained[i] |= 1;
    constrained[i+N*(N+1)] |= 1;
  }
  
  int k = 0;
  for (int i = 0; i < N; ++i)  {
    for (int j = 0; j < N; ++j, k += 2)  {
   
      abc[k][0] = i*(N+1)+j;
      abc[k][2] = (i+1)*(N+1)+j;
      abc[k][1] = i*(N+1)+j+1;

      abc[k+1][0] = i*(N+1)+j+1;
      abc[k+1][2] = (i+1)*(N+1)+j;
      abc[k+1][1] = (i+1)*(N+1)+j+1;
     
    }
  }

  memset(xyzdot, 0, sizeof(double)*3*numNodes);
}

void TriangulatedInterface::setExactSquare(double loc, double vel) {
 
  for (int i = 0; i < numNodes; ++i) {
    xyz[i*3] = loc;
    xyzdot[i*3] = vel;
  }
}
  
void TriangulatedInterface::initializeIntersector(IoData& iod,  Communicator* com, 
                                                  Domain* dom, DistSVec<double,3>& X) {

  pIntersector = new DistIntersectorPhysBAM(iod, com, numNodes,
                                            xyz, numElems, abc, NULL);

  pIntersector->initialize(dom, X, X, iod, NULL);
}

void TriangulatedInterface::update(double dt) {

  pIntersector->updateStructure(xyz, xyzdot, numNodes);

 // std::cout << "Recomputing intersections ... " << std::endl;
  pIntersector->recompute(dt,dt,dt,true); 
}

class DistLevelSetStructure* 
TriangulatedInterface::getIntersector() {

  return pIntersector;
}

class LevelSetStructure* TriangulatedInterface::
getSubLSS(int iSub) {

  return &(*pIntersector)(iSub);
}



