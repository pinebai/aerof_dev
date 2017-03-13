#ifndef _DES_TERM_H_
#define _DES_TERM_H_

#include <IoData.h>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::max;
#endif

//------------------------------------------------------------------------------

class DESTerm {

  double oorey;
  
// Included (MB)
  double dRe_mudMach;
  
protected:

  double alpha;

  double cdes;
  double cb1;
  double cb2;
  double cw1;
  double cw2;
  double cw3_pow6;
  double opcw3_pow;
  double cv1_pow3;
  double oocv2;
  double oosigma;
  double oovkcst2;
  bool usefv3;

public:

  DESTerm(IoData &);
  ~DESTerm() {}

  double computeTurbulentViscosity(double *[4], double, double &);
  double computeTurbulentViscosity(double *, double);
  double computeSecondTurbulentViscosity(double lambdal, double mul, double mut);
  double computeDerivativeOfSecondTurbulentViscosity(double lambdal, double dlambdal, double mul, double dmul, double mut, double dmut);
  double max(double a, double b) { return (a>b) ? a : b; }
  double min(double a, double b) { return (a<b) ? a : b; }

  template<int neq, int shift>
  void computeJacobianVolumeTermDES(double [4][3], double [4], double [3][3], double, 
				   double, double *[4], double (*)[3][neq][neq],
				   double (*)[neq][neq],  SVec<double,3> &, int [4]);

// Included (MB)
  double computeDerivativeOfTurbulentViscosity(double *[4], double *[4], double, double, double &, double &);
  double computeDerivativeOfTurbulentViscosity(double *, double *, double, double);
  void rstVarDES(IoData &);

};

//------------------------------------------------------------------------------

inline
DESTerm::DESTerm(IoData &iod)
{

// Included (MB)
  dRe_mudMach = iod.ref.dRe_mudMach;

  oorey = 1.0 / iod.ref.reynolds_mu;
  alpha = iod.eqs.fluidModel.gasModel.specificHeatRatio / iod.eqs.tc.prandtlTurbulent;

  cb1 = iod.eqs.tc.tm.des.cb1;
  cb2 = iod.eqs.tc.tm.des.cb2;
  cw2 = iod.eqs.tc.tm.des.cw2;
  double cw3 = iod.eqs.tc.tm.des.cw3;
  cw3_pow6 = cw3*cw3*cw3*cw3*cw3*cw3;
  opcw3_pow = pow(1.0 + cw3_pow6, 1.0/6.0);
  double cv1 = iod.eqs.tc.tm.des.cv1;
  cv1_pow3 = cv1*cv1*cv1;
  oocv2 = 1.0 / iod.eqs.tc.tm.des.cv2;
  oosigma = 1.0 / iod.eqs.tc.tm.des.sigma;
  oovkcst2 = 1.0 / (iod.eqs.tc.tm.des.vkcst*iod.eqs.tc.tm.des.vkcst);
  cw1 = cb1*oovkcst2 + (1.0+cb2) * oosigma;
  cdes = iod.eqs.tc.tm.des.cdes;

  cw1 /= iod.ref.reynolds_mu;
  oosigma /= iod.ref.reynolds_mu;

  if (iod.eqs.tc.tm.des.form == DESModelData::FV3)
    usefv3 = true;
  else
    usefv3 = false;

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void DESTerm::rstVarDES(IoData &iod)
{

  dRe_mudMach = iod.ref.dRe_mudMach;
  oorey = 1.0 / iod.ref.reynolds_mu;
  alpha = iod.eqs.fluidModel.gasModel.specificHeatRatio / iod.eqs.tc.prandtlTurbulent;

  cb1 = iod.eqs.tc.tm.des.cb1;
  cb2 = iod.eqs.tc.tm.des.cb2;
  cw2 = iod.eqs.tc.tm.des.cw2;
  double cw3 = iod.eqs.tc.tm.des.cw3;
  cw3_pow6 = cw3*cw3*cw3*cw3*cw3*cw3;
  opcw3_pow = pow(1.0 + cw3_pow6, 1.0/6.0);
  double cv1 = iod.eqs.tc.tm.des.cv1;
  cv1_pow3 = cv1*cv1*cv1;
  oocv2 = 1.0 / iod.eqs.tc.tm.des.cv2;
  oosigma = 1.0 / iod.eqs.tc.tm.des.sigma;
  oovkcst2 = 1.0 / (iod.eqs.tc.tm.des.vkcst*iod.eqs.tc.tm.des.vkcst);
  cw1 = cb1*oovkcst2 + (1.0+cb2) * oosigma;
  cdes = iod.eqs.tc.tm.des.cdes;

  cw1 /= iod.ref.reynolds_mu;
  oosigma /= iod.ref.reynolds_mu;

}

//------------------------------------------------------------------------------

inline
double DESTerm::computeTurbulentViscosity(double *V[4], double mul, double &mutilde)
{

  mutilde = 0.25 * (V[0][0]*V[0][5] + V[1][0]*V[1][5] +
		    V[2][0]*V[2][5] + V[3][0]*V[3][5]);
  double chi = mutilde / mul;
  double chi3 = chi*chi*chi;
  double fv1 = chi3 / (chi3 + cv1_pow3);

  return mutilde*fv1;

}

//------------------------------------------------------------------------------

// Included (MB)
inline
double DESTerm::computeDerivativeOfTurbulentViscosity(double *V[4], double *dV[4], double mul, double dmul, double &mutilde, double &dmutilde)
{

  mutilde = 0.25 * (V[0][0]*V[0][5] + V[1][0]*V[1][5] +
		    V[2][0]*V[2][5] + V[3][0]*V[3][5]);

  dmutilde = 0.25 * (dV[0][0]*V[0][5] + V[0][0]*dV[0][5] + dV[1][0]*V[1][5] + V[1][0]*dV[1][5] + 
		     dV[2][0]*V[2][5] + V[2][0]*dV[2][5] + dV[3][0]*V[3][5] + V[3][0]*dV[3][5]);

  double chi = mutilde / mul;
  double dchi = dmutilde / mul - mutilde / (mul * mul) * dmul;

  double chi3 = chi*chi*chi;
  double dchi3 = 3.0*chi*chi*dchi;

  double fv1 = chi3 / (chi3 + cv1_pow3);
  double dfv1 = dchi3 / (chi3 + cv1_pow3) - chi3 / ( (chi3 + cv1_pow3) * (chi3 + cv1_pow3) ) * dchi3;

  return dmutilde*fv1 + mutilde*dfv1;

}

//------------------------------------------------------------------------------

inline
double DESTerm::computeTurbulentViscosity(double *V, double mul)
{

  double mutilde = V[0]*V[5];
  double chi = mutilde / mul;
  double chi3 = chi*chi*chi;
  double fv1 = chi3 / (chi3 + cv1_pow3);

  return mutilde*fv1;

}

//------------------------------------------------------------------------------

// Included (MB)
inline
double DESTerm::computeDerivativeOfTurbulentViscosity(double *V, double *dV, double mul, double dmul)
{

  double mutilde = V[0]*V[5];
  double dmutilde = dV[0]*V[5] + V[0]*dV[5];

  double chi = mutilde / mul;
  double dchi = dmutilde / mul - mutilde / (mul * mul) * dmul;

  double chi3 = chi*chi*chi;
  double dchi3 = 3.0*chi*chi*dchi;

  double fv1 = chi3 / (chi3 + cv1_pow3);
  double dfv1 = dchi3 / (chi3 + cv1_pow3) - chi3 / ( (chi3 + cv1_pow3) * (chi3 + cv1_pow3) ) * dchi3;

  return dmutilde*fv1 + mutilde*dfv1;

}

//------------------------------------------------------------------------------
inline
double DESTerm::computeSecondTurbulentViscosity(double lambdal, double mul, double mut)
{

  //simple model that remains true when the Stokes' hypothesis is assumed
  return -2.0*mut/3.0;

}

//------------------------------------------------------------------------------

inline
double DESTerm::computeDerivativeOfSecondTurbulentViscosity(double lambdal, double dlambdal,
    double mul, double dmul, double mut, double dmut)
{

  return -2.0*dmut/3.0;

}

//------------------------------------------------------------------------------

template<int neq, int shift>
void DESTerm::computeJacobianVolumeTermDES(double dp1dxj[4][3], double d2w[4], 
					 double dudxj[3][3], double mul, double mutilde, 
					 double *V[4], double (*dRdU)[3][neq][neq], 
					 double (*dSdU)[neq][neq], SVec<double,3> &X, int nodeNum[4])
{

  const double sixth = 1.0/6.0;

  double absmutilde = fabs(mutilde);
  double maxmutilde = max(mutilde, 0.0);
  double dabsmutilde,dmaxmutilde;
  if (mutilde != 0.0) 
    dabsmutilde = fabs(mutilde)/mutilde;
  else
    dabsmutilde = 0.0;

  if (maxmutilde == 0.0) 
    dmaxmutilde = 0.0;
  else
    dmaxmutilde = 1.0;
      
  double mu5 = oosigma * (mul + absmutilde);
  double dnutildedx = dp1dxj[0][0]*V[0][5] + dp1dxj[1][0]*V[1][5] + 
    dp1dxj[2][0]*V[2][5] + dp1dxj[3][0]*V[3][5];
  double dnutildedy = dp1dxj[0][1]*V[0][5] + dp1dxj[1][1]*V[1][5] + 
    dp1dxj[2][1]*V[2][5] + dp1dxj[3][1]*V[3][5];
  double dnutildedz = dp1dxj[0][2]*V[0][5] + dp1dxj[1][2]*V[1][5] + 
    dp1dxj[2][2]*V[2][5] + dp1dxj[3][2]*V[3][5];
  double drdx = oosigma * 0.25 * dnutildedx;
  double drdy = oosigma * 0.25 * dnutildedy;
  double drdz = oosigma * 0.25 * dnutildedz;

  int k;
  for (k=0; k<4; ++k) {
    double nu = mu5 / V[k][0];
    dRdU[k][0][shift][shift] = drdx + nu * dp1dxj[k][0];
    dRdU[k][1][shift][shift] = drdy + nu * dp1dxj[k][1];
    dRdU[k][2][shift][shift] = drdz + nu * dp1dxj[k][2];
  }

  double maxl,sidel;
  maxl=-1.0;

  for (int i=0; i<4; i++) {
    for (int j=i+1; j<4; j++){
      sidel=sqrt((X[nodeNum[i]][0]-X[nodeNum[j]][0])*(X[nodeNum[i]][0]-X[nodeNum[j]][0]) +
	       (X[nodeNum[i]][1]-X[nodeNum[j]][1])*(X[nodeNum[i]][1]-X[nodeNum[j]][1]) +
               (X[nodeNum[i]][2]-X[nodeNum[j]][2])*(X[nodeNum[i]][2]-X[nodeNum[j]][2]));
      maxl = max(maxl,sidel);
    }
  }

  double d2wall = 0.25 * (d2w[0] + d2w[1] + d2w[2] + d2w[3]);
  d2wall = min(d2wall,cdes*maxl);

  if (d2wall < 1.e-15) {
    for (k=0; k<4; ++k)
      dSdU[k][shift][shift] = 0.0;
    return;
  }
  double chi = max(mutilde/mul, 0.001);
  double chi3 = chi*chi*chi;
  double fv1 = chi3 / (chi3 + cv1_pow3);
  double fv2  = 1.-chi/(1.+chi*fv1);
  double fv3  = 1.0;
  if (usefv3) {
    fv2 = 1.0 + oocv2*chi;
    fv2 = 1.0 / (fv2*fv2*fv2);
    fv3 = (1.0 + chi*fv1) * (1.0 - fv2) / chi;
  }
  double ood2wall2 = 1.0 / (d2wall * d2wall);
  double rho = 0.25 * (V[0][0] + V[1][0] + V[2][0] + V[3][0]);
  double oorho = 1.0 / rho;
  double zz = oorey * oovkcst2 * maxmutilde * oorho * ood2wall2;
  double s12 = dudxj[0][1] - dudxj[1][0];
  double s23 = dudxj[1][2] - dudxj[2][1];
  double s31 = dudxj[2][0] - dudxj[0][2];
  double s = sqrt(s12*s12 + s23*s23 + s31*s31);
  double Stilde = max(s*fv3 + zz*fv2,1.0e-12); // To avoid possible numerical problems, the term \tilde S must never be allowed to reach zero or go negative.
  double rr = min(zz/Stilde, 2.0);
  double rr2 = rr*rr;
  double gg = rr + cw2 * (rr2*rr2*rr2 - rr);
  double gg2 = gg*gg;
  double fw = opcw3_pow * gg * pow(gg2*gg2*gg2 + cw3_pow6, -sixth);

  double chi2 = chi*chi;
  double dchi = 1.0 / mul;
  if (chi == 0.001) dchi = 0.0;
  double coef1 = 1.0 / (chi3 + cv1_pow3);
  double dfv1 = 3.0*chi2*dchi*cv1_pow3 * coef1*coef1;
  double coef2 = 1.0 / (1.0 + chi*oocv2);
  double coef3 = coef2 * coef2;
  double dfv2 = (fv2-1.)*dchi/chi+(1.-fv2)*(1-fv2)*(dfv1+fv1*dchi/chi);
  double dfv3 = 0;
  if (usefv3) {
    dfv2 = -3.0*dchi*oocv2 * coef3*coef3;
    dfv3 = ((dchi*fv1 + chi*dfv1)*(1.0 - fv2) - 
	   (1.0 + chi*fv1)*dfv2 - fv3*dchi) / chi;
  }
  double dStilde = s*dfv3 + oorey*oovkcst2*oorho*ood2wall2 * (fv2*dmaxmutilde + maxmutilde*dfv2);
  if (Stilde == 1.0e-12) dStilde = 0.0;
  double drr = oorey*oovkcst2*oorho*ood2wall2 * (Stilde*dmaxmutilde - maxmutilde*dStilde) / (Stilde*Stilde);
  if (rr == 2.0) drr = 0.0;
  double dgg = (1.0 + cw2 * (6.0*rr2*rr2*rr - 1.0)) * drr;
  double dfw = pow(gg2*gg2*gg2 + cw3_pow6, 7.0*sixth);
  dfw = cw3_pow6 * opcw3_pow * dgg / dfw;

  double P = cb1 * Stilde * dabsmutilde;
  double D = cw1 * fw * oorho * maxmutilde * dmaxmutilde * ood2wall2;
  double dP = cb1 * dStilde * absmutilde;
  double dD = cw1 * oorho * ood2wall2 * (fw * maxmutilde * dmaxmutilde + maxmutilde * maxmutilde * dfw);
  double s00 = 0.25 * (max(D - P, 0.0) + max(dD - dP, 0.0));
  //s00 = 0.25 * (D - P + (dD - dP) * maxmutilde);
  double coef4 = oosigma * cb2 * rho * 2.0;

  for (k=0; k<4; ++k)
    dSdU[k][shift][shift] = coef4 / V[k][0] * 
      (dnutildedx*dp1dxj[k][0] + dnutildedy*dp1dxj[k][1] + dnutildedz*dp1dxj[k][2]) - s00;
  
}

//------------------------------------------------------------------------------

#endif
