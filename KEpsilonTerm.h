#ifndef _K_EPSILON_TERM_H_
#define _K_EPSILON_TERM_H_

#include <IoData.h>

//------------------------------------------------------------------------------
/*
@ARTICLE{jones-launder-72,
  author = "Jones, W. P. and Launder, B. E.",
  title = "The Prediction of Laminarization with a Two-Equation
           Turbulence Model",
  journal = "International Journal of Heat and Mass Transfer",
  year = 1972,
  volume = 15,
  pages = "301--314",
}
*/
class KEpsilonTerm {

private:

// Included (MB)
  double dRe_mudMach;

  double reynolds;
  double ooreynoldsKE;

protected:

  double alpha;

  double sigma_k;
  double sigma_eps;
  double sigma_eps1;
  double sigma_eps2;
  double c_mu;

public:

  KEpsilonTerm(IoData &);
  ~KEpsilonTerm() {}

  double computeTurbulentViscosity(double *[4], double &, double &);
  double computeTurbulentViscosity(double *);
  double computeSecondTurbulentViscosity(double lambdal, double mul, double mut);
  double computeDerivativeOfSecondTurbulentViscosity(double lambdal, double dlambdal, double mul, double dmul, double mut, double dmut);

  template<int neq, int shift>
  void computeJacobianVolumeTermKE(double [4][3], double, double, double, double,
				   double *[4], double (*)[3][neq][neq], double (*)[neq][neq]);

// Included (MB)
  double computeDerivativeOfTurbulentViscosity(double *[4], double *[4], double &, double &, double);
  double computeDerivativeOfTurbulentViscosity(double *, double *, double);
  void rstVarKE(IoData &);

};

//------------------------------------------------------------------------------

inline
KEpsilonTerm::KEpsilonTerm(IoData &iod)
{

// Included (MB)
  dRe_mudMach = iod.ref.dRe_mudMach;

  reynolds = iod.ref.reynolds_mu;
  ooreynoldsKE = 1.0/ reynolds;
  alpha = iod.eqs.fluidModel.gasModel.specificHeatRatio / iod.eqs.tc.prandtlTurbulent;

  sigma_k = iod.eqs.tc.tm.ke.sigma_k;
  sigma_eps = iod.eqs.tc.tm.ke.sigma_eps;
  sigma_eps1 = iod.eqs.tc.tm.ke.sigma_eps1;
  sigma_eps2 = iod.eqs.tc.tm.ke.sigma_eps2;
  c_mu = iod.eqs.tc.tm.ke.c_mu;

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void KEpsilonTerm::rstVarKE(IoData &iod)
{

  reynolds = iod.ref.reynolds_mu;
  dRe_mudMach = iod.ref.dRe_mudMach;
  ooreynoldsKE = 1.0/reynolds;

}

//------------------------------------------------------------------------------

inline
double KEpsilonTerm::computeTurbulentViscosity(double *V[4], double &rhok, double &rhoeps)
{

  static const double fourth = 1.0 / 4.0;

  rhok = fourth * (V[0][0]*V[0][5] + V[1][0]*V[1][5] +
		   V[2][0]*V[2][5] + V[3][0]*V[3][5]);
  rhoeps = fourth * (V[0][0]*V[0][6] + V[1][0]*V[1][6] +
		     V[2][0]*V[2][6] + V[3][0]*V[3][6]);

  return reynolds * c_mu * rhok*rhok / rhoeps;

}

//------------------------------------------------------------------------------

// Included (MB)
inline
double KEpsilonTerm::computeDerivativeOfTurbulentViscosity(double *V[4], double *dV[4], double &drhok, double &drhoeps, double dMach)
{

  static const double fourth = 1.0 / 4.0;

  double rhok = fourth * (V[0][0]*V[0][5] + V[1][0]*V[1][5] +
		   V[2][0]*V[2][5] + V[3][0]*V[3][5]);
  drhok = fourth * (dV[0][0]*V[0][5] + V[0][0]*dV[0][5] + dV[1][0]*V[1][5] + V[1][0]*dV[1][5] +
		   dV[2][0]*V[2][5] + V[2][0]*dV[2][5] + dV[3][0]*V[3][5] + V[3][0]*dV[3][5]);
  double rhoeps = fourth * (V[0][0]*V[0][6] + V[1][0]*V[1][6] +
		     V[2][0]*V[2][6] + V[3][0]*V[3][6]);
  drhoeps = fourth * (dV[0][0]*V[0][6] + V[0][0]*dV[0][6] + dV[1][0]*V[1][6] + V[1][0]*dV[1][6] +
		     dV[2][0]*V[2][6] + V[2][0]*dV[2][6] + dV[3][0]*V[3][6] + V[3][0]*dV[3][6]);

  return ( ( dRe_mudMach * dMach * c_mu * rhok*rhok + reynolds * c_mu *
  2.0*rhok*drhok ) * rhoeps - reynolds * c_mu * rhok*rhok * drhoeps ) / (rhoeps * rhoeps);

}

//------------------------------------------------------------------------------

inline
double KEpsilonTerm::computeTurbulentViscosity(double *V)
{

  return reynolds * c_mu * V[0] * V[5]*V[5] / V[6];

}

//------------------------------------------------------------------------------

// Included (MB)
inline
double KEpsilonTerm::computeDerivativeOfTurbulentViscosity(double *V, double *dV, double dMach)
{

  return (dRe_mudMach * dMach * c_mu * V[0] * V[5]*V[5] / V[6] +  reynolds * c_mu * dV[0] * V[5]*V[5] / V[6] + reynolds * c_mu * V[0] * 2.0*V[5]*dV[5] / V[6] - reynolds * c_mu * V[0] * V[5]*V[5] / (V[6] * V[6]) * dV[6]);

}

//------------------------------------------------------------------------------
inline
double KEpsilonTerm::computeSecondTurbulentViscosity(double lambdal, double mul, double mut)
{

  //simple model that remains true when the Stokes' hypothesis is assumed
  return -2.0*mut/3.0;

}

//------------------------------------------------------------------------------

inline
double KEpsilonTerm::computeDerivativeOfSecondTurbulentViscosity(double lambdal, double dlambdal,
    double mul, double dmul, double mut, double dmut)
{

  return -2.0*dmut/3.0;

}

//------------------------------------------------------------------------------

template<int neq, int shift>
void KEpsilonTerm::computeJacobianVolumeTermKE(double dp1dxj[4][3], double mul,
					       double mut, double rhok, double rhoeps, 
					       double *V[4], double (*dRdU)[3][neq][neq], 
					       double (*dSdU)[neq][neq])
{

  double muk = ooreynoldsKE * (mul + sigma_k * mut);
  double mueps = ooreynoldsKE * (mul + sigma_eps * mut);

  int k;
  for (k=0; k<4; ++k) {
    double oorho = 1.0 / V[k][0];
    double nuk = muk * oorho;
    double nueps = mueps * oorho;

    dRdU[k][0][0 + shift][0 + shift] = nuk * dp1dxj[k][0];
    dRdU[k][0][0 + shift][1 + shift] = 0.0;
    dRdU[k][0][1 + shift][0 + shift] = 0.0;
    dRdU[k][0][1 + shift][1 + shift] = nueps * dp1dxj[k][0];

    dRdU[k][1][0 + shift][0 + shift] = nuk * dp1dxj[k][1];
    dRdU[k][1][0 + shift][1 + shift] = 0.0;
    dRdU[k][1][1 + shift][0 + shift] = 0.0;
    dRdU[k][1][1 + shift][1 + shift] = nueps * dp1dxj[k][1];

    dRdU[k][2][0 + shift][0 + shift] = nuk * dp1dxj[k][2];
    dRdU[k][2][0 + shift][1 + shift] = 0.0;
    dRdU[k][2][1 + shift][0 + shift] = 0.0;
    dRdU[k][2][1 + shift][1 + shift] = nueps * dp1dxj[k][2];
  }

  double s00 = - rhoeps / rhok;
  double s01 = -1.0;
  double s10 = sigma_eps2 * s00*s00;
  double s11 = 2.0 * sigma_eps2 * s00;
  s00 *= 0.25;
  s01 *= 0.25;
  s10 *= 0.25;
  s11 *= 0.25;

  for (k=0; k<4; ++k) {
    //dSdU[k][0 + shift][0 + shift] = s00;
    dSdU[k][0 + shift][0 + shift] = 0.0;
    dSdU[k][0 + shift][1 + shift] = s01;
    dSdU[k][1 + shift][0 + shift] = s10;
    dSdU[k][1 + shift][1 + shift] = s11;
  }

}

//------------------------------------------------------------------------------

#endif
