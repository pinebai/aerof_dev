#ifndef _DIST_GEO_STATE_H_
#define _DIST_GEO_STATE_H_

#include <GeoData.h>

class TimeData;
class Domain;
class GeoState;
class DistInfo;
class Communicator;

struct Vec3D;

template<class Scalar> class DistVec;
template<class Scalar, int dim> class DistSVec;
template<class Scalar, int dim, int dim2> class RectangularSparseMat;

//------------------------------------------------------------------------------

class DistGeoState {

  int numLocSub;
  double lscale;
  double oolscale;

  GeoData data;

  DistSVec<double,3> *Xn;	// nodal positions at time n
  DistSVec<double,3> *Xnm1;	// nodal positions at time n-1
  DistSVec<double,3> *Xnm2;	// nodal positions at time n-2
  DistSVec<double,3> *Xdot;
  DistSVec<double,3> *Xsave;

  DistVec<double> *ctrlVol_n;
  DistVec<double> *ctrlVol_nm1;
  DistVec<double> *ctrlVol_nm2;
  DistVec<double> *ctrlVol_save;

  DistVec<double> *d2wall;

  DistVec<Vec3D>  *edgeNorm;
  DistVec<double> *edgeNormVel;
  DistVec<Vec3D>  *edgeNorm_nm1;
  DistVec<double> *edgeNormVel_nm1;
  DistVec<Vec3D>  *edgeNorm_nm2;
  DistVec<double> *edgeNormVel_nm2;

  //There is one normal per face node (however for triangular nodes only one normal 
  //is stored as they are all equal). The mapping to the face normals is stored with 
  //each face (normNum in Face). The total number of face normals is stored in 
  //the faceset (numFaceNorms in FaceSet).
  DistVec<Vec3D>  *faceNorm;
  DistVec<double> *faceNormVel;
  DistVec<Vec3D>  *faceNorm_nm1;
  DistVec<double> *faceNormVel_nm1;
  DistVec<Vec3D>  *faceNorm_nm2;
  DistVec<double> *faceNormVel_nm2;

  DistVec<Vec3D> *inletNodeNorm;
  DistVec<int> *numFaceNeighb; 	   //number of faces connected to an inletnode, 
                                   //independantly of subdomains
 
  
  Domain *domain;

  GeoState **subGeoState;

  Communicator *com;

// Included (MB)
  int optFlag;
  DistSVec<double,3> *Xsa;
  DistSVec<double,3> *dXsa;
  DistVec<Vec3D> *dEdgeNorm;
  DistVec<Vec3D> *dFaceNorm;
  DistVec<double> *dEdgeNormVel;
  DistVec<double> *dFaceNormVel;

public:

  DistGeoState(IoData &, Domain *);

  // These two are used (only?) for multigrid
  DistGeoState(const GeoData &, Domain *, DistInfo& nodeDistInfo, DistInfo& edgeDistInfo,
               DistInfo& faceNormDistInfo);
  DistGeoState(const GeoData &, Domain *, DistInfo& nodeDistInfo, DistInfo& edgeDistInfo);

  ~DistGeoState();

  GeoState &operator() (int i) const { return *subGeoState[i]; }
  const GeoData& getGeoData() const { return data; }

  void setup(const char *, TimeData &, DistSVec<double,3> *, DistVec<double> *);
  void setup1(const char *, DistSVec<double,3> *, DistVec<double> *);
  void setup2(TimeData &);
  void setup3(const char *, DistSVec<double,3> *, DistVec<double> *);
  void setupInitialDisplacement(const char *, DistSVec<double,3> *, DistVec<double> *);

  void compute(TimeData &, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double> &);
  void interpolate(double, double, DistSVec<double,3> &, DistSVec<double,3> &);
  void update(DistSVec<double,3> &, DistVec<double> &);
  void writeToDisk(char *);

  int getConfig() const { return data.config; }
  DistSVec<double,3> &getXn() const { return *Xn; }
  DistVec<double> &getCtrlVol() const { return *ctrlVol_n; }
  DistVec<Vec3D> &getFaceNormal() const { return *faceNorm; }
  DistVec<Vec3D> &getEdgeNormal() const { return *edgeNorm; }
  DistVec<double> &getEdgeNormalVel() const { return *edgeNormVel; }
  DistVec<double> &getFaceNorVel() const { return *faceNormVel; }
  DistVec<double> *getd2wall() const { return d2wall; }
  DistVec<Vec3D> &getInletNodeNorm() const { return *inletNodeNorm; }

// Included (MB)
  void updateConfigSA() { data.configSA += 1; }
  void resetConfigSA() { data.configSA = 0; }
  int getConfigSA() const { return data.configSA; }
  void computeDerivatives(DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double> &);
  void computeDerivatives(RectangularSparseMat<double,3,3> **, RectangularSparseMat<double,3,3> **, RectangularSparseMat<double,3,1> **,
		                  DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double> &, DistVec<Vec3D>&, DistVec<Vec3D>&, DistVec<double>&);
  void computeTransposeDerivatives(RectangularSparseMat<double,3,3> **, RectangularSparseMat<double,3,3> **, RectangularSparseMat<double,3,1> **,
                                   DistVec<double>&, DistVec<Vec3D>&, DistVec<Vec3D>&, DistVec<double>&, DistSVec<double,3> &);
  void computeDerivativeOperators(DistSVec<double,3> &, RectangularSparseMat<double,3,3> **, RectangularSparseMat<double,3,3> **, RectangularSparseMat<double,3,1> **); 
  void reset(DistSVec<double,3> &);

  DistVec<Vec3D> &getdEdgeNormal() const { return *dEdgeNorm; }
  DistVec<double> &getdEdgeNormalVel() const { return *dEdgeNormVel; }
  

};

//------------------------------------------------------------------------------

#endif
