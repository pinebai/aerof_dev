#ifndef _WALL_FCN_H_
#define _WALL_FCN_H_

class IoData;
class VarFcn;
class ViscoFcn;

// Included (MB)
class Communicator;

struct Vec3D;

//------------------------------------------------------------------------------

class WallFcn {

  const static double third;
  const static double eleventh;

  double gam;
  double prandtl;

  VarFcn *varFcn;
  ViscoFcn *viscoFcn;

protected:

  double vkcst;

// Included (MB)
  double dRedMach;

  double reynolds;

  virtual void computeWallValues(double utau, double delta, double rho,
				 double ut, double mu, double *V) {}

// Included (MB)
  virtual void computeDerivativeOfWallValues(double utau, double dutau, double delta, double rho, double drho,
				 double ut, double dut, double mu, double dmu, double dMach, double *V, double *dV) {}

private:

  Vec3D computeTangentVector(Vec3D &, Vec3D &);
  void computeFaceValues(double [3], double *, double *[3], double &, Vec3D &, 
			 double &, Vec3D &, double &, double &);
  double computeFrictionVelocity(double, double, double, double);
  double computeFrictionTemperature(double, double, double, double, double);

// Included (MB)
  Vec3D computeDerivativeOfTangentVector(Vec3D &, Vec3D &, Vec3D &, Vec3D &);
  void computeDerivativeOfFaceValues(double [3], double *, double *, double *[3], double *[3], double &, Vec3D &,
			 double &, Vec3D &, double &, double &, double &);
  double computeDerivativeOfFrictionVelocity(Vec3D &, Vec3D &, double, double, double, Vec3D &, Vec3D &, double, double, double);
  double computeDerivativeOfFrictionTemperature(double, double, double, double, double, double, double, double, double, double);

public:

  WallFcn(IoData &, VarFcn *, ViscoFcn *);
  ~WallFcn() {}

  void computeSurfaceTerm(int, Vec3D &, double [3], double *, double *[3], double *);
  Vec3D computeForce(Vec3D &, double [3], double *, double *[3]);
  double computeInterfaceWork(Vec3D&, double [3], double*, double* [3]);
  double computeHeatPower(Vec3D &, double [3], double *, double *[3]);
  double computeDeltaPlus(Vec3D &, double [3], double *, double *[3]);

  double computedudT(double rho, double T, double ut, double dT, 
						 double d2w, double &dudn, double &dTdn);

  template<int neq>
  void computeJacobianSurfaceTerm(int, Vec3D &, double [3], double *, 
				  double *[3], double (*)[neq][neq]);


// Included (MB)
  template<int neq>
  void computeBCsJacobianWallValues(int, Vec3D &, double [3], double *, double *, double *[3]);
  void computeDerivativeOfSurfaceTerm(int, Vec3D &, Vec3D &, double [3], double *, double *, double *[3], double *[3], double, double *);
  void rstVar(IoData &, Communicator *);
  Vec3D computeDerivativeOfForce(Vec3D &, Vec3D &, double [3], double *, double *, double *[3], double *[3], double);
  double computeDerivativeOfHeatPower(Vec3D &, Vec3D &, double [3], double *, double *, double *[3], double *[3], double);

};

//------------------------------------------------------------------------------
/*
@ARTICLE{reichardt-42,
  author = "Reichardt, H.",
  title = "Gesetzm{\"a}ssigkeiten der freien {T}urbulenz",
  journal = "VDI-Forschungsheft 414, 1st edition, Berlin",
  year = 1942,
} 
*/

class WallFcnSA : public WallFcn {

  double cv1_pow3;

public:

  WallFcnSA(IoData &, VarFcn *, ViscoFcn *);
  ~WallFcnSA() {}

  void computeWallValues(double, double, double, double, double, double *);

// Included (MB)
  void computeDerivativeOfWallValues(double, double, double, double, double, double, double, double, double, double, double *, double *);

};

//------------------------------------------------------------------------------

class WallFcnKE : public WallFcn {

  double orcmu;

public:

  WallFcnKE(IoData &, VarFcn *, ViscoFcn *);
  ~WallFcnKE() {}

  void computeWallValues(double, double, double, double, double, double *);

// Included (MB)
  void computeDerivativeOfWallValues(double, double, double, double, double, double, double, double, double, double, double *, double *);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <WallFcn.C>
#endif

#endif
