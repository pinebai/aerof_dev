// BasicGeometry.h

#include "Vector3D.h"

class BasicGeometry {

  public:

   static double ComputeTetVolume(const Vec3D x[4]);

   static double ComputeWedgeVolume(const Vec3D x[6]);

};
