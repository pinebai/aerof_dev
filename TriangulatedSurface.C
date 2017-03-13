#include <Vector.h>
#include <Face.h>
#include <TriangulatedSurface.h>

//------------------------------------------------------------------------

TriangulatedSurface::TriangulatedSurface()
{
  numTriangle = 0;
  //triangleList = NULL;
  triangleList = 0;
  triangleNodeNum = 0;
}

//-------------------------------------------------------------------------

TriangulatedSurface::~TriangulatedSurface()
{
  if (triangleList) delete[] triangleList;
  if (triangleNodeNum) delete[] triangleNodeNum; //wrong! need to be modified.
}

//-------------------------------------------------------------------------

void TriangulatedSurface::addTriangle( Vec3D &Node1, Vec3D &Node2, Vec3D &Node3 ) // add a new triangle.
{
  numTriangle++;
  triangleList[ numTriangle-1 ][0] = Node1;
  triangleList[ numTriangle-1 ][1] = Node2;
  triangleList[ numTriangle-1 ][2] = Node3;
}

//--------------------------------------------------------------------------

void TriangulatedSurface::addTriangle(int Node1, int Node2, int Node3)
{
  numTriangle++;
  triangleNodeNum[ numTriangle-1][0] = Node1;
  triangleNodeNum[ numTriangle-1][1] = Node2;
  triangleNodeNum[ numTriangle-1][2] = Node3;
}


//----------------------------------------------------------------------------

void TriangulatedSurface::removeTriangle(int j)
{
  // :(
}

//-----------------------------------------------------------------------------

void TriangulatedSurface::addSurfaceFromFace( SVec<double,3> &X, FaceSet &boundaryFaces )
{
  int totTriagulatedSurfacesToCreate = 0;
  for (int i=0; i<boundaryFaces.size(); i++) {
    if (boundaryFaces[i].type() != Face::TRIA) return; //should also output an error message.
    if (boundaryFaces[i].getCode() == BC_ADIABATIC_WALL_MOVING
        || boundaryFaces[i].getCode() == BC_ADIABATIC_WALL_FIXED
        || boundaryFaces[i].getCode() == BC_SLIP_WALL_MOVING
        || boundaryFaces[i].getCode() == BC_SLIP_WALL_FIXED
        || boundaryFaces[i].getCode() == BC_ISOTHERMAL_WALL_MOVING
        || boundaryFaces[i].getCode() == BC_ISOTHERMAL_WALL_FIXED
        || boundaryFaces[i].getCode() == BC_POROUS_WALL_MOVING
        || boundaryFaces[i].getCode() == BC_POROUS_WALL_FIXED )
      totTriagulatedSurfacesToCreate++;
  }
  if(!triangleList) triangleList = new Vec3D[totTriagulatedSurfacesToCreate][3];

  for (int i=0; i<boundaryFaces.size(); i++) {
    if (boundaryFaces[i].type() != Face::TRIA) return; //should also output an error message.
    if (boundaryFaces[i].getCode() == BC_ADIABATIC_WALL_MOVING
        || boundaryFaces[i].getCode() == BC_ADIABATIC_WALL_FIXED
        || boundaryFaces[i].getCode() == BC_SLIP_WALL_MOVING
        || boundaryFaces[i].getCode() == BC_SLIP_WALL_FIXED
        || boundaryFaces[i].getCode() == BC_ISOTHERMAL_WALL_MOVING
        || boundaryFaces[i].getCode() == BC_ISOTHERMAL_WALL_FIXED
        || boundaryFaces[i].getCode() == BC_POROUS_WALL_MOVING
        || boundaryFaces[i].getCode() == BC_POROUS_WALL_FIXED ){

       Vec3D A( X[ boundaryFaces[i][0] ]);
       Vec3D B( X[ boundaryFaces[i][1] ]);
       Vec3D C( X[ boundaryFaces[i][2] ]);
       
       addTriangle ( A, B, C);
    }
  }
  
  //for debug
/*  fprintf(stderr, "# of triangles: %d\n", numTriangle);
  for (int i=0; i<numTriangle; i++) {
    double x=(triangleList[i][0])[0];
    double y=(triangleList[i][0])[1];
    double z=(triangleList[i][0])[2];
    fprintf(stderr, "node numbers: [ %.4f, %.4f, %.4f ],   ", x, y, z);

    x=(triangleList[i][1])[0];
    y=(triangleList[i][1])[1];
    z=(triangleList[i][1])[2];
    fprintf(stderr, "node numbers: [ %.4f, %.4f, %.4f ],   ", x, y, z);

    x=(triangleList[i][2])[0];
    y=(triangleList[i][2])[1];
    z=(triangleList[i][2])[2];
    fprintf(stderr, "node numbers: [ %.4f, %.4f, %.4f ]\n", x, y, z);

  }*/
}

//------------------------------------------------------------------------------

void TriangulatedSurface::addSurfaceFromFace( FaceSet &boundaryFaces )
{
  int totTriagulatedSurfacesToCreate = 0;
  for (int i=0; i<boundaryFaces.size(); i++) {
    if (boundaryFaces[i].type() != Face::TRIA) return; //should also output an error message.
    if (boundaryFaces[i].getCode() == BC_ADIABATIC_WALL_MOVING
        || boundaryFaces[i].getCode() == BC_ADIABATIC_WALL_FIXED
        || boundaryFaces[i].getCode() == BC_SLIP_WALL_MOVING
        || boundaryFaces[i].getCode() == BC_SLIP_WALL_FIXED
        || boundaryFaces[i].getCode() == BC_ISOTHERMAL_WALL_MOVING
        || boundaryFaces[i].getCode() == BC_ISOTHERMAL_WALL_FIXED
        || boundaryFaces[i].getCode() == BC_POROUS_WALL_MOVING
        || boundaryFaces[i].getCode() == BC_POROUS_WALL_FIXED )
      totTriagulatedSurfacesToCreate++;
  }
  if(!triangleNodeNum) triangleNodeNum = new int[totTriagulatedSurfacesToCreate][3];

  for (int i=0; i<boundaryFaces.size(); i++) {
    if (boundaryFaces[i].type() != Face::TRIA) return; //should also output an error message.
    if (boundaryFaces[i].getCode() == BC_ADIABATIC_WALL_MOVING
        || boundaryFaces[i].getCode() == BC_ADIABATIC_WALL_FIXED
        || boundaryFaces[i].getCode() == BC_SLIP_WALL_MOVING
        || boundaryFaces[i].getCode() == BC_SLIP_WALL_FIXED
        || boundaryFaces[i].getCode() == BC_ISOTHERMAL_WALL_MOVING
        || boundaryFaces[i].getCode() == BC_ISOTHERMAL_WALL_FIXED
        || boundaryFaces[i].getCode() == BC_POROUS_WALL_MOVING
        || boundaryFaces[i].getCode() == BC_POROUS_WALL_FIXED ){

       int A = boundaryFaces[i][0];
       int B = boundaryFaces[i][1];
       int C = boundaryFaces[i][2];

       addTriangle ( A, B, C);
    }
  }
//for debug
//  fprintf(stderr, "# of triangles: %d\n", numTriangle);
//  for (int i=0; i<numTriangle; i++)
//    fprintf(stderr, "node numbers: [ %d, %d, %d ]\n", triangleNodeNum[i][0], triangleNodeNum[i][1],
//                                                      triangleNodeNum[i][2]);
}


//-----------------------------------------------------------------------------

//Vec3D TriangulatedSurface::(*getTriangleList() const)[3] { return triangleList; }

//-------------------------------------------------------------------------------

//int TriangulatedSurface::getNumTriangle() { return numTriangle; }

//-------------------------------------------------------------------------------




















