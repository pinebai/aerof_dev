#ifndef _GEO_DATA_H_
#define _GEO_DATA_H_

#include <IoData.h>

//------------------------------------------------------------------------------

class GeoData {

public:

  DGCLData::Normals typeNormals;
  DGCLData::Velocities typeVelocities;

  int config;

// Included (MB)
  int configSA;

  bool use_n;
  bool use_nm1;
  bool use_nm2;
  bool use_save;

public:

  GeoData(IoData &);
  GeoData(const GeoData &);
  ~GeoData() {}

};

//------------------------------------------------------------------------------

#endif
