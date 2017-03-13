#include <PolygonReconstructionData.h>

#include <Elem.h>
#include <LevelSet/LevelSetStructure.h>
#include <Vector3D.h>

#include <cstdlib>
#include <cstdio>
//------------------------------------------------------------------------------

int getPolygons(Elem &elem, LevelSetStructure &LSS, PolygonReconstructionData* polygons)
{
    int numberOfPolygons=0,intersectedEdges=0;
    int oppositeNodes[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
    int T[4]; for (int i=0; i<4; ++i) T[i] = elem[i]; //nodes in a tet.
    int isBlocked[4][4];
    int edge[4][4];
    int parity[4] = {0,0,0,0};
    for(int i=0; i<6; ++i) {
        int l = elem.edgeNum(i);
        int ni=elem.edgeEnd(i,0), nj=elem.edgeEnd(i,1);
        isBlocked[ni][nj] = isBlocked[nj][ni] = (LSS.edgeIntersectsStructure(0, l) ? 1 : 0);
        edge[ni][nj] = edge[nj][ni] = l;
        if(isBlocked[ni][nj] == 1){
            ++parity[ni]; ++parity[nj];
            ++intersectedEdges;
        }
    }

    switch (intersectedEdges) {
    case 6:{
        numberOfPolygons = 4;
        polygons[0].AssignTriangleSingle(T[0], T[1], T[2], T[3], edge[0][1], edge[0][2], edge[0][3]);
        polygons[1].AssignTriangleSingle(T[1], T[0], T[2], T[3], edge[0][1], edge[1][2], edge[1][3]);
        polygons[2].AssignTriangleSingle(T[2], T[1], T[0], T[3], edge[2][1], edge[0][2], edge[2][3]);
        polygons[3].AssignTriangleSingle(T[3], T[1], T[2], T[0], edge[3][1], edge[3][2], edge[0][3]);
        break;}
    case 5:{
        numberOfPolygons = 3;
        int c=0,n=3,e0=0,e1=3;
        while(parity[c] != 2) ++c;
        while(parity[e0] != 3) ++e0;
        while(parity[n] != 2) --n;
        while(parity[e1] != 3) --e1;
        polygons[0].AssignTriangleSingle(T[e0], T[c], T[n], T[e1], edge[e0][c], edge[e0][n], edge[e0][e1]);
        polygons[1].AssignTriangleSingle(T[e1], T[c], T[n], T[e0], edge[e1][c], edge[e1][n], edge[e1][e0]);
        polygons[2].AssignQuadrilateral(T[c], T[n], T[e0], T[e1], edge[c][e0], edge[c][e1], edge[n][e1], edge[n][e0]);
        break;}
    case 4:{
        for(int i=0;i<4;++i) if(parity[i] == 3) {numberOfPolygons = 4;
            int e0 = oppositeNodes[i][0], e1 = oppositeNodes[i][1], e2 = oppositeNodes[i][2];
            polygons[0].AssignTriangleSingle(T[i], T[e0], T[e1], T[e2], edge[i][e0], edge[i][e1], edge[i][e2]);
            if(isBlocked[e0][e1]) {
                polygons[1].AssignQuadTriangle(T[e0], T[e2], T[i], T[e1], edge[e0][i], edge[e0][e1], edge[e2][i]);
                polygons[2].AssignQuadTriangle(T[e1], T[e2], T[i], T[e0], edge[e1][i], edge[e1][e0], edge[e2][i]);
                polygons[3].AssignConnectingQuadTriangle(T[e0], T[e1], T[e2], T[i], edge[e0][e1], edge[e2][i]);
            } else if(isBlocked[e0][e2]) {
                polygons[1].AssignQuadTriangle(T[e0], T[e1], T[i], T[e2], edge[e0][i], edge[e0][e2], edge[e1][i]);
                polygons[2].AssignQuadTriangle(T[e2], T[e1], T[i], T[e0], edge[e2][i], edge[e2][e0], edge[e1][i]);
                polygons[3].AssignConnectingQuadTriangle(T[e0], T[e2], T[e1], T[i], edge[e0][e2], edge[e1][i]);
            } else {
                polygons[1].AssignQuadTriangle(T[e1], T[e0], T[i], T[e2], edge[e1][i], edge[e1][e2], edge[e0][i]);
                polygons[2].AssignQuadTriangle(T[e2], T[e0], T[i], T[e1], edge[e2][i], edge[e2][e1], edge[e0][i]);
                polygons[3].AssignConnectingQuadTriangle(T[e1], T[e2], T[e0], T[i], edge[e1][e2], edge[e0][i]);
            }
            break;
        }
        if(numberOfPolygons == 0) {
            numberOfPolygons = 2;
            int n,e0,e1;
            for(n=1 ;  n<4; ++n) if(!isBlocked[0][n]) break;
            for(e0=1; e0<4;++e0) if(e0 != n) break;
            for(e1=3;e1>=0;--e1) if(e1 != n) break;
            polygons[0].AssignQuadrilateral(T[0],  T[n],  T[e0], T[e1], edge[0][e0], edge[0][e1], edge[n][e1], edge[n][e0]);
            polygons[1].AssignQuadrilateral(T[e0], T[e1], T[0],  T[n],  edge[e0][0], edge[e0][n], edge[e1][n], edge[e1][0]);
        }
        break;}
    case 3:{
        for(int i=0;i<4;++i) {
            if(parity[i] == 0) {numberOfPolygons = 2;
                int e0 = oppositeNodes[i][0], e1 = oppositeNodes[i][1], e2 = oppositeNodes[i][2];
                polygons[0].AssignTwoEdges(       T[e0], T[e1], T[e2], T[i], edge[e0][e1], edge[e0][e2]);
                polygons[1].AssignTwoEdgesPartTwo(T[e0], T[e1], T[e2], T[i], edge[e0][e1], edge[e0][e2], edge[e1][e2]);
                break;
            } else if(parity[i] == 3) {numberOfPolygons = 2;
                int e0 = oppositeNodes[i][0], e1 = oppositeNodes[i][1], e2 = oppositeNodes[i][2];
                polygons[0].AssignTriangleSingle(T[i], T[e0], T[e1], T[e2], edge[i][e0], edge[i][e1], edge[i][e2]);
                polygons[1].AssignTriangleMulti( T[i], T[e0], T[e1], T[e2], edge[i][e0], edge[i][e1], edge[i][e2]);
                break;
            }
        }
        if(numberOfPolygons == 0) {numberOfPolygons = 4;
            int c=0,n=3,e0=0,e1=3;
            while(parity[c] != 2) ++c;
            while(parity[e0] != 1) ++e0;
            while(parity[n] != 2) --n;
            while(parity[e1] != 1) --e1;
            if(isBlocked[c][e0]) {
                polygons[0].AssignQuadTriangle(T[c], T[e1], T[n], T[e0], edge[c][n], edge[c][e0], edge[e1][n]);
                polygons[1].AssignQuadTriangle(T[n], T[e0], T[c], T[e1], edge[n][c], edge[n][e1], edge[e0][c]);
                polygons[2].AssignConnectingQuadTriangle(T[c], T[e0] ,T[e1], T[n], edge[c][e0], edge[n][e1]);
                polygons[3].AssignConnectingQuadTriangle(T[n], T[e1] ,T[e0], T[c], edge[n][e1], edge[c][e0]);
            } else {
                polygons[0].AssignQuadTriangle(T[c], T[e0], T[n], T[e1], edge[c][n], edge[c][e1], edge[e0][n]);
                polygons[1].AssignQuadTriangle(T[n], T[e1], T[c], T[e0], edge[n][c], edge[n][e0], edge[e1][c]);
                polygons[2].AssignConnectingQuadTriangle(T[c], T[e1] ,T[e0], T[n], edge[c][e1], edge[n][e0]);
                polygons[3].AssignConnectingQuadTriangle(T[n], T[e0] ,T[e1], T[c], edge[n][e0], edge[c][e1]);
            }
        }
        break;}
    case 2:{
        for(int i=0;i<4;++i) {
            if(parity[i] == 2) {numberOfPolygons = 2;
                int e0 = oppositeNodes[i][0], e1 = oppositeNodes[i][1], e2 = oppositeNodes[i][2];
                if(!isBlocked[i][e0])      polygons[0].AssignTwoEdges(T[i], T[e1], T[e2], T[e0], edge[i][e1], edge[i][e2]);
                else if(!isBlocked[i][e1]) polygons[0].AssignTwoEdges(T[i], T[e0], T[e2], T[e1], edge[i][e0], edge[i][e2]);
                else                       polygons[0].AssignTwoEdges(T[i], T[e0], T[e1], T[e2], edge[i][e0], edge[i][e1]);
                break;
            }
        }
        break;}
    case 1:
    case 0:
    default:
        break;
    }
    return numberOfPolygons;
}


//--------------------------------------------------------------------------

void getPolygonNormal(SVec<double,3>& X, Vec3D &normal, LevelSetStructure &LSS, PolygonReconstructionData &polygon)
{
  int nEdges = polygon.numberOfEdges;
  std::vector<LevelSetResult> lsRes(nEdges);
  std::vector<Vec3D> Xinter(nEdges);
  Vec3D start_vertex;for(int m=0;m<3;++m) start_vertex[m]=X[polygon.nodeToLookFrom][m];
  for(int m=0; m<nEdges; m++) {
    lsRes[m] = LSS.getLevelSetDataAtEdgeCenter(0, polygon.edge[m], polygon.edgeWithVertex[m][0]<polygon.edgeWithVertex[m][1]);
    double alpha = lsRes[m].alpha;
    if (alpha<0) {fprintf(stderr,"Unable to get intersection results at edge center! Abort...\n"); exit(-1);}
    Xinter[m][0] = alpha*X[polygon.edgeWithVertex[m][0]][0] + (1.0-alpha)*X[polygon.edgeWithVertex[m][1]][0];
    Xinter[m][1] = alpha*X[polygon.edgeWithVertex[m][0]][1] + (1.0-alpha)*X[polygon.edgeWithVertex[m][1]][1];
    Xinter[m][2] = alpha*X[polygon.edgeWithVertex[m][0]][2] + (1.0-alpha)*X[polygon.edgeWithVertex[m][1]][2];
  }
  switch(nEdges) {
    case 3: // got a triangle.
      normal = (Xinter[1]-Xinter[0])^(Xinter[2]-Xinter[0]);
      if (normal.norm() != 0.0) {normal = 1.0/normal.norm()*normal;}
      // Then we check the orientation of the normal.
      for (int i=0; i<3; i++) {start_vertex[i] = X[polygon.edgeWithVertex[0][0]][i];}
      if (LSS.isActive(0, polygon.edgeWithVertex[0][0])) {
        if (normal*(start_vertex-Xinter[0]) < 0) {normal = -1.0*normal;}
      } else if (normal*(start_vertex-Xinter[0]) > 0) {normal = -1.0*normal;}
      break;
    case 4: // got a quadrangle... We cut it into two triangles and we average the two resulting normals. 
      Vec3D normal1, normal2;
      normal1 = (Xinter[1]-Xinter[0])^(Xinter[3]-Xinter[0]);
      if (normal1.norm() != 0.0) {normal1 = 1.0/normal1.norm()*normal1;}
      normal2 = (Xinter[2]-Xinter[1])^(Xinter[3]-Xinter[1]);
      if (normal2.norm() != 0.0) {normal2 = 1.0/normal2.norm()*normal2;}
      normal = 0.5*(normal1+normal2);
      if (normal.norm() != 0.0) {normal = 1.0/normal.norm()*normal;}
      // Then we check the orientation of the normal.
      for (int i=0; i<3; i++) {start_vertex[i] = X[polygon.edgeWithVertex[0][0]][i];}
      if (LSS.isActive(0, polygon.edgeWithVertex[0][0])) {
        if (normal*(start_vertex-Xinter[0]) < 0) {normal = -1.0*normal;}
      } else if (normal*(start_vertex-Xinter[0]) > 0) {normal = -1.0*normal;}
      break;
  }
}

