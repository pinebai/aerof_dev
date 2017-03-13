#ifndef _TRIANGULATEDSURFACE_H_
#define _TRIANGULATEDSURFACE_H_

#include <Vector.h>
#include <Face.h>

class TriangulatedSurface {
  int numTriangle;  //Node: This is used for both coords-list and node# list. 
  Vec3D (*triangleList)[3];
  int (*triangleNodeNum)[3];

public:
  TriangulatedSurface();
  ~TriangulatedSurface();

  void addTriangle(Vec3D &, Vec3D &, Vec3D &); //add to triangleLIst;
  void addTriangle(int, int, int); // add to triangleNodeNum
  void removeTriangle(int );
  void addSurfaceFromFace( SVec<double,3> &, FaceSet &); //add to triangleList
  void addSurfaceFromFace( FaceSet &);  // add to triangleNodeNum
  Vec3D (*getTriangleList() const)[3] {return triangleList; }
  int (*getTriangleNodeNum() const)[3] {return triangleNodeNum; }
  int getNumTriangle() {return numTriangle; }
};

#endif

