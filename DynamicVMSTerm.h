#ifndef _DYNAMIC_VMS_TERM_H_
#define _DYNAMIC_VMS_TERM_H_

#include <IoData.h>
#include <Vector.h>
#include <WallFcn.h>
#include <NavierStokesTerm.h>

//------------------------------------------------------------------------

class DynamicVMSTerm: public NavierStokesTerm {

protected:

  double gamma;
  double oogamma1;
  double PrT;
  double coeff;
  double csprime;
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

  DynamicVMSTerm(IoData &, VarFcn *);
  ~DynamicVMSTerm();

  void compute(double [4], double [4], double, double [4][3], double *[4],
	       double *[4], double *, SVec<double,3> &, int [4], bool);

  double computeEddyViscosity(double, double[3][3]);

  double computeFilterWidth(double, SVec<double,3> &, int [4]);

  double max(double a, double b) { return (a>b) ? a : b; }

  void computeVelocityGradient(double [4][3], double [4][3], double [3][3]);

  void computeStressTensor(double, double [3][3], double[3][3]);

  void getPrimeValues(double *[4], double *[4], double [4][3], double [4]);

  void computeSmallEnergyGradient(double [4][3], double [4], double [3]);


  double computeMutOMu(double, double [4][3], double *[4], double *[4], 
                       double [4], SVec<double,3> &, int [4]);


};

//------------------------------------------------------------------------

#endif
