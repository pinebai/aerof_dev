#ifndef _DYNAMIC_LES_TERM_H_
#define _DYNAMIC_LES_TERM_H_

#include <IoData.h>
#include <Vector.h>
#include <WallFcn.h>
#include <NavierStokesTerm.h>

#include <cmath>

//-------------------------------------------------------------------------

class DynamicLESTerm: public NavierStokesTerm {

protected:

  double gamma;
  double oogamma1;
  double R;
  double ooR;
  double coeff;
  double csprime; 
  double Prsgs;
  double csmax, ptmin, ptmax;
  double Prandtl;
  double mach;

  double onethird;
  double onefourth;
  double twothirds;
  double x0,y0,z0,x1,y1,z1;
  bool trip;
  bool delta;
  WallFcn* wallFcn;

public:

  DynamicLESTerm(IoData &, VarFcn *);
  ~DynamicLESTerm();

  void compute(double, double[4][3], double *[4], double [4], double [4], double *, SVec<double,3> &, int [4]);

  void computeDelta(double, SVec<double,3> &, int[4], double &);

  double computeEddyViscosity(double, double, double [3][3]);

  double computeFilterWidth(double, SVec<double,3> &, int [4]);

  double max(double a, double b) { return (a>b) ? a : b; }

  void computeVelocity(double *[4], double[4][3], double[3]);

  void computeVelocityGradient(double [4][3], double [4][3], double [3][3]);

  void computeStressTensor(double, double [3][3], double[3][3]);

  void computeEnergy(double *[4], double [4]);

  void computeEnergyGradient(double [4][3], double [4], double [3]);

  double computeMutOMu(double, double [4][3], double *[4], double [4], double [4], SVec<double,3> &, int [4]);

  double outputCsValues(double, double [4], SVec<double,3> &, int [4]);

};

//-------------------------------------------------------------------------

#endif
