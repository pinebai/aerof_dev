#ifndef _WALE_LES_TERM_H_
#define _WALE_LES_TERM_H_

#include <IoData.h>
#include <Vector.h>
#include <WallFcn.h>
#include <NavierStokesTerm.h>

#include <cmath>

//-------------------------------------------------------------------------

class WaleLESTerm : public NavierStokesTerm {

protected:

  double gamma;
  double oogamma1;
  double PrT;
  double coeff;
  double Cw;

  double onethird;
  double onefourth;
  double twothirds;
  double onesixth;
  double x0,y0,z0,x1,y1,z1;
  bool trip;
  bool delta;
  WallFcn* wallFcn;

public:

  WaleLESTerm(IoData &, VarFcn *vf);
  ~WaleLESTerm();

  void compute(double, double[4][3], double *[4], double *, SVec<double,3> &, int [4]);

  double computeEddyViscosity(double, double, double [3][3]);

  double computeFilterWidth(double, SVec<double,3> &, int [4]);

  double max(double a, double b) { return (a>b) ? a : b; }

  void computeVelocity(double *[4], double[4][3], double[3]);

  void computeVelocityGradient(double [4][3], double [4][3], double [3][3]);

  void computeStressTensor(double, double [3][3], double[3][3]);

  void computeEnergy(double *[4], double [4]);

  void computeEnergyGradient(double [4][3], double [4], double [3]);

  double computeMutOMu(double , double [4][3], double *[4], SVec<double,3> &, int [4]);

};

//-------------------------------------------------------------------------

#endif
