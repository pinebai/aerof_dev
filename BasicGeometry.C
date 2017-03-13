// BasicGeometry.C

#include "BasicGeometry.h"
#include <cassert>

double BasicGeometry::ComputeTetVolume(const Vec3D x[4]) {

  Vec3D v1 = x[1] - x[0];
  Vec3D v2 = x[2] - x[0];
  Vec3D v3 = x[3] - x[0];

  return 1.0/6.0 * (v3 * (v1 ^ v2));
}

// Base of the wedge is (x[0],x[1],x[2])
// Top is (x[3],x[4],x[5])
double BasicGeometry::ComputeWedgeVolume(const Vec3D x[6]) {

  // Compute the volume by decomposing into three tets
  Vec3D T1[4] = {x[0],x[1],x[2],x[3]};
  Vec3D T2[4] = {x[1],x[2],x[3],x[4]};
  Vec3D T3[4] = {x[2],x[3],x[4],x[5]};

  double v1 = ComputeTetVolume(T1);
  double v2 = ComputeTetVolume(T2);
  double v3 = ComputeTetVolume(T3);

  //assert(v1 > 0.0 && v2 > 0.0 && v3 > 0.0);
  return v1+v2+v3;
}

