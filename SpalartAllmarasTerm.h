#ifndef _SA_TERM_H_
#define _SA_TERM_H_

#include <IoData.h>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::max;
using std::min;
#endif

//------------------------------------------------------------------------------
/*
@UNPUBLISHED{spalart-allmaras-92,
  author = "Spalart, P. R. and Allmaras, S. R.",
  title = "A One-Equation Turbulence Model for Aerodynamic Flows",
  month = jan,
  year = 1992,
  note = "{AIAA} paper 92-0439",
} 

Multigrid for the 2-D compressible Navier-Stokes equations
AIAA Paper 99-3336 
Steven R. Allmaras (Boeing Co., Seattle, WA)   
AIAA, Aerospace Sciences Meeting and Exhibit, 37th, Reno, NV, Jan. 11-14, 1999
*/

class SATerm {

  double oorey;
  
protected:

  double alpha;

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

  SATerm(IoData &);
  ~SATerm() {}

  double computeTurbulentViscosity(double *[4], double, double &);
  double computeTurbulentViscosity(double *, double);
  double computeSecondTurbulentViscosity(double lambdal, double mul, double mut);
  double computeDerivativeOfSecondTurbulentViscosity(double lambdal, double dlambdal, double mul, double dmul, double mut, double dmut);

  template<int neq, int shift>
  void computeJacobianVolumeTermSA(double [4][3], double [4], double [3][3], double, 
				   double, double *[4], double (*)[3][neq][neq],
				   double (*)[neq][neq]);

// Included (MB)
  double computeDerivativeOfTurbulentViscosity(double *[4], double *[4], double, double, double &);

  double computeDerivativeOfTurbulentViscosity(double *[4],double, double, double);

  double computeDerivativeOfTurbulentViscosity(double *, double *, double, double);

  template<int neq, int shift>
  void computeJacobianVolumeTermSA(double [4][3], double [4], double [3][3], double, double [4][6],
				   double, double [4][6], double *[4], double (*)[3][neq][neq],
				   double (*)[neq][neq]);

  void rstVarSA(IoData &);

};

//------------------------------------------------------------------------------

inline
SATerm::SATerm(IoData &iod)
{

  oorey = 1.0 / iod.ref.reynolds_mu;
  alpha = iod.eqs.fluidModel.gasModel.specificHeatRatio / iod.eqs.tc.prandtlTurbulent;

  cb1 = iod.eqs.tc.tm.sa.cb1;
  cb2 = iod.eqs.tc.tm.sa.cb2;
  cw2 = iod.eqs.tc.tm.sa.cw2;
  double cw3 = iod.eqs.tc.tm.sa.cw3;
  cw3_pow6 = cw3*cw3*cw3*cw3*cw3*cw3;
  opcw3_pow = pow(1.0 + cw3_pow6, 1.0/6.0);
  double cv1 = iod.eqs.tc.tm.sa.cv1;
  cv1_pow3 = cv1*cv1*cv1;
  oocv2 = 1.0 / iod.eqs.tc.tm.sa.cv2;
  oosigma = 1.0 / iod.eqs.tc.tm.sa.sigma;
  oovkcst2 = 1.0 / (iod.eqs.tc.tm.sa.vkcst*iod.eqs.tc.tm.sa.vkcst);
  cw1 = cb1*oovkcst2 + (1.0+cb2) * oosigma;

  cw1 /= iod.ref.reynolds_mu;
  oosigma /= iod.ref.reynolds_mu;

  if (iod.eqs.tc.tm.sa.form == SAModelData::FV3)
    usefv3 = true;
  else
    usefv3 = false;

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void SATerm::rstVarSA(IoData &iod)
{

  oorey = 1.0 / iod.ref.reynolds_mu;

  oosigma = 1.0 / iod.eqs.tc.tm.sa.sigma;
  cw1 = cb1*oovkcst2 + (1.0+cb2) * oosigma;

  cw1 /= iod.ref.reynolds_mu;
  oosigma /= iod.ref.reynolds_mu;

}

//------------------------------------------------------------------------------

inline
double SATerm::computeTurbulentViscosity(double *V[4], double mul, double &mutilde)
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
double SATerm::computeDerivativeOfTurbulentViscosity(double *V[4], double *dV[4], double mul, double dmul, double &dmutilde)
{

  double mutilde = 0.25 * (V[0][0]*V[0][5] + V[1][0]*V[1][5] +
		    V[2][0]*V[2][5] + V[3][0]*V[3][5]);
  dmutilde = 0.25 * (dV[0][0]*V[0][5] + V[0][0]*dV[0][5] + dV[1][0]*V[1][5] + V[1][0]*dV[1][5] +
		    dV[2][0]*V[2][5] + V[2][0]*dV[2][5] + dV[3][0]*V[3][5] + V[3][0]*dV[3][5]);
  double chi = mutilde / mul;
  double dchi = (dmutilde * mul - mutilde * dmul ) / ( mul * mul );
  double chi3 = chi*chi*chi;
  double fv1 = chi3 / (chi3 + cv1_pow3);
  double dfv1 = cv1_pow3 * 3.0 * chi * chi * dchi / ((chi3 + cv1_pow3) * (chi3 + cv1_pow3));


  return dmutilde*fv1 + mutilde*dfv1;

}

//------------------------------------------------------------------------------
inline
double SATerm::computeSecondTurbulentViscosity(double lambdal, double mul, double mut)
{

  //simple model that remains true when the Stokes' hypothesis is assumed
  return -2.0*mut/3.0;

}

//------------------------------------------------------------------------------

inline
double SATerm::computeDerivativeOfSecondTurbulentViscosity(double lambdal, double dlambdal,
    double mul, double dmul, double mut, double dmut)
{

  return -2.0*dmut/3.0;

}

//------------------------------------------------------------------------------

// Included (MB)
inline
double SATerm::computeDerivativeOfTurbulentViscosity(double *V[4], double mul, double dmul, double dmutilde)
{

  double mutilde = 0.25 * (V[0][0]*V[0][5] + V[1][0]*V[1][5] +
		    V[2][0]*V[2][5] + V[3][0]*V[3][5]);
  double chi = mutilde / mul;
  double dchi = (dmutilde * mul - mutilde * dmul ) / ( mul * mul );
  double chi3 = chi*chi*chi;
  double fv1 = chi3 / (chi3 + cv1_pow3);
  double dfv1 = cv1_pow3 * 3.0 * chi * chi * dchi / ((chi3 + cv1_pow3) * (chi3 + cv1_pow3));


  return dmutilde*fv1 + mutilde*dfv1;

}

//------------------------------------------------------------------------------

inline
double SATerm::computeTurbulentViscosity(double *V, double mul)
{

	//double mutilde = V[0]*V[5];
	double mutilde = V[0] * std::max(0.0,V[5]);
  double chi = mutilde / mul;
  double chi3 = chi*chi*chi;
  double fv1 = chi3 / (chi3 + cv1_pow3);

  return mutilde*fv1;

}

//------------------------------------------------------------------------------

// Included (MB)
inline
double SATerm::computeDerivativeOfTurbulentViscosity(double *V, double *dV, double mul, double dmul)
{

  double mutilde = V[0]*V[5];
  double dmutilde = dV[0]*V[5] + V[0]*dV[5];
  double chi = mutilde / mul;
  double dchi = (dmutilde*mul - mutilde*dmul) / (mul*mul);
  double chi3 = chi*chi*chi;
  double dchi3 = 3.0*chi*chi*dchi;
  double fv1 = chi3 / (chi3 + cv1_pow3);
  double dfv1 = (dchi3*(chi3 + cv1_pow3)-chi3*dchi3) / ((chi3 + cv1_pow3)*(chi3 + cv1_pow3));

  return dmutilde*fv1 + mutilde*dfv1;

}

//------------------------------------------------------------------------------

template<int neq, int shift>
void SATerm::computeJacobianVolumeTermSA(double dp1dxj[4][3], double d2w[4], 
					 double dudxj[3][3], double mul, double mutilde, 
					 double *V[4], double (*dRdU)[3][neq][neq], 
					 double (*dSdU)[neq][neq])
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

  double d2wall = 0.25 * (d2w[0] + d2w[1] + d2w[2] + d2w[3]);
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

// Included (MB)
template<int neq, int shift>
void SATerm::computeJacobianVolumeTermSA(double dp1dxj[4][3], double d2w[4], 
					 double dudxj[3][3], double mul, double dmul[4][6], double mutilde, 
					 double dmutilde[4][6], double *V[4], double (*dRdU)[3][neq][neq], 
					 double (*dSdU)[neq][neq])
{
  int k;

  const double sixth = 1.0/6.0;

  double absmutilde = fabs(mutilde);
  double dabsmutilde[4][6];
  if  (mutilde != 0.0)
    for (k=0; k<4; ++k) {
      dabsmutilde[k][0] = ( fabs(mutilde) / mutilde ) * dmutilde[k][0];
      dabsmutilde[k][1] = ( fabs(mutilde) / mutilde ) * dmutilde[k][1];
      dabsmutilde[k][2] = ( fabs(mutilde) / mutilde ) * dmutilde[k][2];
      dabsmutilde[k][3] = ( fabs(mutilde) / mutilde ) * dmutilde[k][3];
      dabsmutilde[k][4] = ( fabs(mutilde) / mutilde ) * dmutilde[k][4];
      dabsmutilde[k][5] = ( fabs(mutilde) / mutilde ) * dmutilde[k][5];
    }
  else {
    //  fprintf(stderr, "***** Inside the file SpalartAllmarasTerm.h the varible mutilde is zero *****\n");
    //exit(1);
    for (k=0; k<4; ++k) {
      dabsmutilde[k][0] = 0.0;
      dabsmutilde[k][1] = 0.0;
      dabsmutilde[k][2] = 0.0;
      dabsmutilde[k][3] = 0.0;
      dabsmutilde[k][4] = 0.0;
      dabsmutilde[k][5] = 0.0;
    }
  }

  double maxmutilde = max(mutilde, 0.0);
  double dmaxmutilde[4][6];
  if ( maxmutilde == 0.0 )  {
    for (k=0; k<4; ++k) {
      dmaxmutilde[k][0] = 0.0;
      dmaxmutilde[k][1] = 0.0;
      dmaxmutilde[k][2] = 0.0;
      dmaxmutilde[k][3] = 0.0;
      dmaxmutilde[k][4] = 0.0;
      dmaxmutilde[k][5] = 0.0;
    }
    //   fprintf(stderr, "***** Inside the file SpalartAllmarasTerm.h the varible mutilde is zero *****\n");
    //exit(1);
  }
  else {
    for (k=0; k<4; ++k) {
      dmaxmutilde[k][0] = dmutilde[k][0];
      dmaxmutilde[k][1] = dmutilde[k][1];
      dmaxmutilde[k][2] = dmutilde[k][2];
      dmaxmutilde[k][3] = dmutilde[k][3];
      dmaxmutilde[k][4] = dmutilde[k][4];
      dmaxmutilde[k][5] = dmutilde[k][5];
    }
  }

  double mu5 = oosigma * (mul + absmutilde);
  double dmu5[4][6];
  for (k=0; k<4; ++k) {
    dmu5[k][0] = oosigma * (dmul[k][0] + dabsmutilde[k][0]);
    dmu5[k][1] = oosigma * (dmul[k][1] + dabsmutilde[k][1]);
    dmu5[k][2] = oosigma * (dmul[k][2] + dabsmutilde[k][2]);
    dmu5[k][3] = oosigma * (dmul[k][3] + dabsmutilde[k][3]);
    dmu5[k][4] = oosigma * (dmul[k][4] + dabsmutilde[k][4]);
    dmu5[k][5] = oosigma * (dmul[k][5] + dabsmutilde[k][5]);
  }
  
  double dnutildedx = dp1dxj[0][0]*V[0][5] + dp1dxj[1][0]*V[1][5] + dp1dxj[2][0]*V[2][5] + dp1dxj[3][0]*V[3][5];
  double dnutildedy = dp1dxj[0][1]*V[0][5] + dp1dxj[1][1]*V[1][5] + dp1dxj[2][1]*V[2][5] + dp1dxj[3][1]*V[3][5];
  double dnutildedz = dp1dxj[0][2]*V[0][5] + dp1dxj[1][2]*V[1][5] + dp1dxj[2][2]*V[2][5] + dp1dxj[3][2]*V[3][5];

  for (k=0; k<4; ++k) {
    dRdU[k][0][5][0] = dmu5[k][0] * dnutildedx - mu5 * dp1dxj[k][0] * V[k][5] / V[k][0];
    dRdU[k][0][5][1] = dmu5[k][1] * dnutildedx;
    dRdU[k][0][5][2] = dmu5[k][2] * dnutildedx;
    dRdU[k][0][5][3] = dmu5[k][3] * dnutildedx;
    dRdU[k][0][5][4] = dmu5[k][4] * dnutildedx;
    dRdU[k][0][5][5] = dmu5[k][5] * dnutildedx + mu5 * dp1dxj[k][0] / V[k][0];

    dRdU[k][1][5][0] = dmu5[k][0] * dnutildedy - mu5 * dp1dxj[k][1] * V[k][5] / V[k][0];
    dRdU[k][1][5][1] = dmu5[k][1] * dnutildedy;
    dRdU[k][1][5][2] = dmu5[k][2] * dnutildedy;
    dRdU[k][1][5][3] = dmu5[k][3] * dnutildedy;
    dRdU[k][1][5][4] = dmu5[k][4] * dnutildedy;
    dRdU[k][1][5][5] = dmu5[k][5] * dnutildedy + mu5 * dp1dxj[k][1] / V[k][0];

    dRdU[k][2][5][0] = dmu5[k][0] * dnutildedz - mu5 * dp1dxj[k][2] * V[k][5] / V[k][0];
    dRdU[k][2][5][1] = dmu5[k][1] * dnutildedz;
    dRdU[k][2][5][2] = dmu5[k][2] * dnutildedz;
    dRdU[k][2][5][3] = dmu5[k][3] * dnutildedz;
    dRdU[k][2][5][4] = dmu5[k][4] * dnutildedz;
    dRdU[k][2][5][5] = dmu5[k][5] * dnutildedz + mu5 * dp1dxj[k][2] / V[k][0];
  }

  double d2wall = 0.25 * (d2w[0] + d2w[1] + d2w[2] + d2w[3]);
  if (d2wall < 1.e-15) {
    for (k=0; k<4; ++k)
      dSdU[k][5][5] = 0.0;
    return;
  }
  double chi;
  double dchi[4][6];
  double chi3;
  double fv1;
  double dfv1[4][6];
  double fv2;
  double dfv2[4][6];
  double fv3;
  double dfv3[4][6];
  double ood2wall2;
  double rho;
  double drho[4][6];
  double oorho;
  double doorho[4][6];
  double zz;
  double dzz[4][6];
  double s12;
  double ds12[4][6];
  double s23;
  double ds23[4][6];
  double s31;
  double ds31[4][6];
  double s;
  double ds[4][6];
  double Stilde;
  double dStilde[4][6];
  double rr;
  double drr[4][6];
  double rr2;
  double gg;
  double dgg[4][6];
  double gg2;
  double fw;
  double dfw[4][6];

//  double AA;
  double dAA[4][6];
//  double BB;
  double dBB[4][6];
//  double CC;
  double dCC[4][6];

  rho = 0.25 * (V[0][0] + V[1][0] + V[2][0] + V[3][0]);

  for (k=0; k<4; ++k) {

    chi = max(mutilde/mul, 0.001);
    if (chi == 0.001) {
      dchi[k][0] = 0.0;
      dchi[k][1] = 0.0;
      dchi[k][2] = 0.0;
      dchi[k][3] = 0.0;
      dchi[k][4] = 0.0;
      dchi[k][5] = 0.0;
    }
    else {
      dchi[k][0] = - mutilde / (mul * mul) * dmul[k][0] + dmutilde[k][0] / mul;
      dchi[k][1] = - mutilde / (mul * mul) * dmul[k][1] + dmutilde[k][1] / mul;
      dchi[k][2] = - mutilde / (mul * mul) * dmul[k][2] + dmutilde[k][2] / mul;
      dchi[k][3] = - mutilde / (mul * mul) * dmul[k][3] + dmutilde[k][3] / mul;
      dchi[k][4] = - mutilde / (mul * mul) * dmul[k][4] + dmutilde[k][4] / mul;
      dchi[k][5] = - mutilde / (mul * mul) * dmul[k][5] + dmutilde[k][5] / mul;
    }

    chi3 = chi*chi*chi;
    fv1 = chi3 / (chi3 + cv1_pow3);
    dfv1[k][0] = ( 3.0*chi*chi*dchi[k][0]*(chi3 + cv1_pow3) - chi3 * 3.0*chi*chi*dchi[k][0] ) / ( (chi3 + cv1_pow3) * (chi3 + cv1_pow3) );
    dfv1[k][1] = ( 3.0*chi*chi*dchi[k][1]*(chi3 + cv1_pow3) - chi3 * 3.0*chi*chi*dchi[k][1] ) / ( (chi3 + cv1_pow3) * (chi3 + cv1_pow3) );
    dfv1[k][2] = ( 3.0*chi*chi*dchi[k][2]*(chi3 + cv1_pow3) - chi3 * 3.0*chi*chi*dchi[k][2] ) / ( (chi3 + cv1_pow3) * (chi3 + cv1_pow3) );
    dfv1[k][3] = ( 3.0*chi*chi*dchi[k][3]*(chi3 + cv1_pow3) - chi3 * 3.0*chi*chi*dchi[k][3] ) / ( (chi3 + cv1_pow3) * (chi3 + cv1_pow3) );
    dfv1[k][4] = ( 3.0*chi*chi*dchi[k][4]*(chi3 + cv1_pow3) - chi3 * 3.0*chi*chi*dchi[k][4] ) / ( (chi3 + cv1_pow3) * (chi3 + cv1_pow3) );
    dfv1[k][5] = ( 3.0*chi*chi*dchi[k][5]*(chi3 + cv1_pow3) - chi3 * 3.0*chi*chi*dchi[k][5] ) / ( (chi3 + cv1_pow3) * (chi3 + cv1_pow3) );

    fv2  = 1.-chi/(1.+chi*fv1);
    dfv2[k][0] = (fv2-1.)*dchi[k][0]/chi+(1.-fv2)*(1-fv2)*(dfv1[k][0]+fv1*dchi[k][0]/chi);
    dfv2[k][1] = (fv2-1.)*dchi[k][1]/chi+(1.-fv2)*(1-fv2)*(dfv1[k][1]+fv1*dchi[k][1]/chi);
    dfv2[k][2] = (fv2-1.)*dchi[k][2]/chi+(1.-fv2)*(1-fv2)*(dfv1[k][2]+fv1*dchi[k][2]/chi);
    dfv2[k][3] = (fv2-1.)*dchi[k][3]/chi+(1.-fv2)*(1-fv2)*(dfv1[k][3]+fv1*dchi[k][3]/chi);
    dfv2[k][4] = (fv2-1.)*dchi[k][4]/chi+(1.-fv2)*(1-fv2)*(dfv1[k][4]+fv1*dchi[k][4]/chi);
    dfv2[k][5] = (fv2-1.)*dchi[k][5]/chi+(1.-fv2)*(1-fv2)*(dfv1[k][5]+fv1*dchi[k][5]/chi);
    fv3  = 1.0;
    dfv3[k][0] = 0.; 
    dfv3[k][1] = 0.; 
    dfv3[k][2] = 0.; 
    dfv3[k][3] = 0.; 
    dfv3[k][4] = 0.; 
    dfv3[k][5] = 0.; 

    if (usefv3) {
      fv2 = 1.0 + oocv2*chi;
      dfv2[k][0] = oocv2*dchi[k][0];
      dfv2[k][1] = oocv2*dchi[k][1];
      dfv2[k][2] = oocv2*dchi[k][2];
      dfv2[k][3] = oocv2*dchi[k][3];
      dfv2[k][4] = oocv2*dchi[k][4];
      dfv2[k][5] = oocv2*dchi[k][5];
      dfv2[k][0] = -3.0 / (fv2*fv2*fv2*fv2)*dfv2[k][0];
      dfv2[k][1] = -3.0 / (fv2*fv2*fv2*fv2)*dfv2[k][1];
      dfv2[k][2] = -3.0 / (fv2*fv2*fv2*fv2)*dfv2[k][2];
      dfv2[k][3] = -3.0 / (fv2*fv2*fv2*fv2)*dfv2[k][3];
      dfv2[k][4] = -3.0 / (fv2*fv2*fv2*fv2)*dfv2[k][4];
      dfv2[k][5] = -3.0 / (fv2*fv2*fv2*fv2)*dfv2[k][5];
      fv2 = 1.0 / (fv2*fv2*fv2);

      fv3 = (1.0 + chi*fv1) * (1.0 - fv2) / chi;
      dfv3[k][0] = ( ( dchi[k][0]*fv1 + chi*dfv1[k][0] ) * (1.0 - fv2) * chi + (1.0 + chi*fv1) * (- dfv2[k][0]) * chi - (1.0 + chi*fv1) * (1.0 - fv2) * dchi[k][0] ) / ( chi * chi );
      dfv3[k][1] = ( ( dchi[k][1]*fv1 + chi*dfv1[k][1] ) * (1.0 - fv2) * chi + (1.0 + chi*fv1) * (- dfv2[k][1]) * chi - (1.0 + chi*fv1) * (1.0 - fv2) * dchi[k][1] ) / ( chi * chi );
      dfv3[k][2] = ( ( dchi[k][2]*fv1 + chi*dfv1[k][2] ) * (1.0 - fv2) * chi + (1.0 + chi*fv1) * (- dfv2[k][2]) * chi - (1.0 + chi*fv1) * (1.0 - fv2) * dchi[k][2] ) / ( chi * chi );
      dfv3[k][3] = ( ( dchi[k][3]*fv1 + chi*dfv1[k][3] ) * (1.0 - fv2) * chi + (1.0 + chi*fv1) * (- dfv2[k][3]) * chi - (1.0 + chi*fv1) * (1.0 - fv2) * dchi[k][3] ) / ( chi * chi );
      dfv3[k][4] = ( ( dchi[k][4]*fv1 + chi*dfv1[k][4] ) * (1.0 - fv2) * chi + (1.0 + chi*fv1) * (- dfv2[k][4]) * chi - (1.0 + chi*fv1) * (1.0 - fv2) * dchi[k][4] ) / ( chi * chi );
      dfv3[k][5] = ( ( dchi[k][5]*fv1 + chi*dfv1[k][5] ) * (1.0 - fv2) * chi + (1.0 + chi*fv1) * (- dfv2[k][5]) * chi - (1.0 + chi*fv1) * (1.0 - fv2) * dchi[k][5] ) / ( chi * chi );
    }

    ood2wall2 = 1.0 / (d2wall * d2wall);

    drho[k][0] = 0.25;
    drho[k][1] = 0.0;
    drho[k][2] = 0.0;
    drho[k][3] = 0.0;
    drho[k][4] = 0.0;
    drho[k][5] = 0.0;
    oorho = 1.0 / rho;
    doorho[k][0] = -1.0 / ( rho * rho ) * drho[k][0];	  
    doorho[k][1] = -1.0 / ( rho * rho ) * drho[k][1];	  
    doorho[k][2] = -1.0 / ( rho * rho ) * drho[k][2];	  
    doorho[k][3] = -1.0 / ( rho * rho ) * drho[k][3];	  
    doorho[k][4] = -1.0 / ( rho * rho ) * drho[k][4];	  
    doorho[k][5] = -1.0 / ( rho * rho ) * drho[k][5];	  

    zz = oorey * oovkcst2 * maxmutilde * oorho * ood2wall2;
    dzz[k][0] = oorey*oovkcst2*dmaxmutilde[k][0]*oorho*ood2wall2 + oorey*oovkcst2*mutilde*doorho[k][0]*ood2wall2;
    dzz[k][1] = oorey*oovkcst2*dmaxmutilde[k][1]*oorho*ood2wall2 + oorey*oovkcst2*mutilde*doorho[k][1]*ood2wall2;
    dzz[k][2] = oorey*oovkcst2*dmaxmutilde[k][2]*oorho*ood2wall2 + oorey*oovkcst2*mutilde*doorho[k][2]*ood2wall2;
    dzz[k][3] = oorey*oovkcst2*dmaxmutilde[k][3]*oorho*ood2wall2 + oorey*oovkcst2*mutilde*doorho[k][3]*ood2wall2;
    dzz[k][4] = oorey*oovkcst2*dmaxmutilde[k][4]*oorho*ood2wall2 + oorey*oovkcst2*mutilde*doorho[k][4]*ood2wall2;
    dzz[k][5] = oorey*oovkcst2*dmaxmutilde[k][5]*oorho*ood2wall2 + oorey*oovkcst2*mutilde*doorho[k][5]*ood2wall2;

    s12 = dudxj[0][1] - dudxj[1][0];
    s23 = dudxj[1][2] - dudxj[2][1];
    s31 = dudxj[2][0] - dudxj[0][2];

    ds12[k][0] = - (dp1dxj[k][1]*V[k][1] - dp1dxj[k][0]*V[k][2]) / V[k][0];
    ds12[k][1] = dp1dxj[k][1] / V[k][0];
    ds12[k][2] = - dp1dxj[k][0] / V[k][0];
    ds12[k][3] = 0.0;
    ds12[k][4] = 0.0;
    ds12[k][5] = 0.0;
    
    ds23[k][0] = - (dp1dxj[k][2]*V[k][2] - dp1dxj[k][1]*V[k][3]) / V[k][0];
    ds23[k][1] = 0.0;
    ds23[k][2] = dp1dxj[k][2] / V[k][0];
    ds23[k][3] = - dp1dxj[k][1] / V[k][0];
    ds23[k][4] = 0.0;
    ds23[k][5] = 0.0;

    ds31[k][0] = - (dp1dxj[k][0]*V[k][3] - dp1dxj[k][2]*V[k][1]) / V[k][0];
    ds31[k][1] = - dp1dxj[k][2] / V[k][0];
    ds31[k][2] = 0.0;
    ds31[k][3] = dp1dxj[k][0] / V[k][0];
    ds31[k][4] = 0.0;
    ds31[k][5] = 0.0;

    s = sqrt(s12*s12 + s23*s23 + s31*s31);
    ds[k][0] = 1.0 / s * (s12*ds12[k][0] + s23*ds23[k][0] + s31*ds31[k][0]);
    ds[k][1] = 1.0 / s * (s12*ds12[k][1] + s23*ds23[k][1] + s31*ds31[k][1]);
    ds[k][2] = 1.0 / s * (s12*ds12[k][2] + s23*ds23[k][2] + s31*ds31[k][2]);
    ds[k][3] = 1.0 / s * (s12*ds12[k][3] + s23*ds23[k][3] + s31*ds31[k][3]);
    ds[k][4] = 1.0 / s * (s12*ds12[k][4] + s23*ds23[k][4] + s31*ds31[k][4]);
    ds[k][5] = 1.0 / s * (s12*ds12[k][5] + s23*ds23[k][5] + s31*ds31[k][5]);

    Stilde = max(s*fv3 + zz*fv2,1.0e-12);// To avoid possible numerical problems, the term \tilde S must never be allowed to reach zero or go negative. 
    if (Stilde == 1.0e-12) {
      dStilde[k][0] = 0.0;
      dStilde[k][1] = 0.0;
      dStilde[k][2] = 0.0;
      dStilde[k][3] = 0.0;
      dStilde[k][4] = 0.0;
      dStilde[k][5] = 0.0;
    }
    else {
      dStilde[k][0] = ds[k][0]*fv3 + s*dfv3[k][0] + dzz[k][0]*fv2 + zz*dfv2[k][0];
      dStilde[k][1] = ds[k][1]*fv3 + s*dfv3[k][1] + dzz[k][1]*fv2 + zz*dfv2[k][1];
      dStilde[k][2] = ds[k][2]*fv3 + s*dfv3[k][2] + dzz[k][2]*fv2 + zz*dfv2[k][2];
      dStilde[k][3] = ds[k][3]*fv3 + s*dfv3[k][3] + dzz[k][3]*fv2 + zz*dfv2[k][3];
      dStilde[k][4] = ds[k][4]*fv3 + s*dfv3[k][4] + dzz[k][4]*fv2 + zz*dfv2[k][4];
      dStilde[k][5] = ds[k][5]*fv3 + s*dfv3[k][5] + dzz[k][5]*fv2 + zz*dfv2[k][5];
    }

    rr = min(zz/Stilde, 2.0);
    if (rr==2.0) {
      drr[k][0] = 0.0;
      drr[k][1] = 0.0;
      drr[k][2] = 0.0;
      drr[k][3] = 0.0;
      drr[k][4] = 0.0;
      drr[k][5] = 0.0;
    }
    else {
      drr[k][0] = ( dzz[k][0] * Stilde - zz * dStilde[k][0] ) / ( Stilde * Stilde );
      drr[k][1] = ( dzz[k][1] * Stilde - zz * dStilde[k][1] ) / ( Stilde * Stilde );
      drr[k][2] = ( dzz[k][2] * Stilde - zz * dStilde[k][2] ) / ( Stilde * Stilde );
      drr[k][3] = ( dzz[k][3] * Stilde - zz * dStilde[k][3] ) / ( Stilde * Stilde );
      drr[k][4] = ( dzz[k][4] * Stilde - zz * dStilde[k][4] ) / ( Stilde * Stilde );
      drr[k][5] = ( dzz[k][5] * Stilde - zz * dStilde[k][5] ) / ( Stilde * Stilde );
    }
    rr2 = rr*rr;

    gg = rr + cw2 * (rr2*rr2*rr2 - rr);
    dgg[k][0] = drr[k][0] + cw2 * (6.0*rr*rr2*rr2*drr[k][0] - drr[k][0]);
    dgg[k][1] = drr[k][1] + cw2 * (6.0*rr*rr2*rr2*drr[k][1] - drr[k][1]);
    dgg[k][2] = drr[k][2] + cw2 * (6.0*rr*rr2*rr2*drr[k][2] - drr[k][2]);
    dgg[k][3] = drr[k][3] + cw2 * (6.0*rr*rr2*rr2*drr[k][3] - drr[k][3]);
    dgg[k][4] = drr[k][4] + cw2 * (6.0*rr*rr2*rr2*drr[k][4] - drr[k][4]);
    dgg[k][5] = drr[k][5] + cw2 * (6.0*rr*rr2*rr2*drr[k][5] - drr[k][5]);
    gg2 = gg*gg;

    fw = opcw3_pow * gg * pow(gg2*gg2*gg2 + cw3_pow6, -sixth);
    dfw[k][0] = opcw3_pow * dgg[k][0] * pow(gg2*gg2*gg2 + cw3_pow6, -sixth) + opcw3_pow * gg * (-sixth) * pow(gg2*gg2*gg2 + cw3_pow6, (-sixth - 1.0) ) * 6.0*gg*gg2*gg2*dgg[k][0];
    dfw[k][1] = opcw3_pow * dgg[k][1] * pow(gg2*gg2*gg2 + cw3_pow6, -sixth) + opcw3_pow * gg * (-sixth) * pow(gg2*gg2*gg2 + cw3_pow6, (-sixth - 1.0) ) * 6.0*gg*gg2*gg2*dgg[k][1];
    dfw[k][2] = opcw3_pow * dgg[k][2] * pow(gg2*gg2*gg2 + cw3_pow6, -sixth) + opcw3_pow * gg * (-sixth) * pow(gg2*gg2*gg2 + cw3_pow6, (-sixth - 1.0) ) * 6.0*gg*gg2*gg2*dgg[k][2];
    dfw[k][3] = opcw3_pow * dgg[k][3] * pow(gg2*gg2*gg2 + cw3_pow6, -sixth) + opcw3_pow * gg * (-sixth) * pow(gg2*gg2*gg2 + cw3_pow6, (-sixth - 1.0) ) * 6.0*gg*gg2*gg2*dgg[k][3];
    dfw[k][4] = opcw3_pow * dgg[k][4] * pow(gg2*gg2*gg2 + cw3_pow6, -sixth) + opcw3_pow * gg * (-sixth) * pow(gg2*gg2*gg2 + cw3_pow6, (-sixth - 1.0) ) * 6.0*gg*gg2*gg2*dgg[k][4];
    dfw[k][5] = opcw3_pow * dgg[k][5] * pow(gg2*gg2*gg2 + cw3_pow6, -sixth) + opcw3_pow * gg * (-sixth) * pow(gg2*gg2*gg2 + cw3_pow6, (-sixth - 1.0) ) * 6.0*gg*gg2*gg2*dgg[k][5];

//    AA = oosigma * cb2 * rho * (dnutildedx*dnutildedx + dnutildedy*dnutildedy + dnutildedz*dnutildedz);
    dAA[k][0] = oosigma * cb2 * drho[k][0] * (dnutildedx*dnutildedx + dnutildedy*dnutildedy + dnutildedz*dnutildedz) - oosigma * cb2 * rho * 2.0 * (dnutildedx*dp1dxj[k][0]*V[k][5] + dnutildedy*dp1dxj[k][1]*V[k][5] + dnutildedz*dp1dxj[k][2]*V[k][5]) / V[k][0];
    dAA[k][1] = 0.0;
    dAA[k][2] = 0.0;
    dAA[k][3] = 0.0;
    dAA[k][4] = 0.0;
    dAA[k][5] = oosigma * cb2 * rho * 2.0 * (dnutildedx*dp1dxj[k][0] + dnutildedy*dp1dxj[k][1] + dnutildedz*dp1dxj[k][2]) / V[k][0];

//    BB = cb1 * Stilde * absmutilde;
    dBB[k][0] = cb1 * dStilde[k][0] * absmutilde + cb1 * Stilde * dabsmutilde[k][0];
    dBB[k][1] = cb1 * dStilde[k][1] * absmutilde + cb1 * Stilde * dabsmutilde[k][1];
    dBB[k][2] = cb1 * dStilde[k][2] * absmutilde + cb1 * Stilde * dabsmutilde[k][2];
    dBB[k][3] = cb1 * dStilde[k][3] * absmutilde + cb1 * Stilde * dabsmutilde[k][3];
    dBB[k][4] = cb1 * dStilde[k][4] * absmutilde + cb1 * Stilde * dabsmutilde[k][4];
    dBB[k][5] = cb1 * dStilde[k][5] * absmutilde + cb1 * Stilde * dabsmutilde[k][5];

//    CC = - cw1 * fw * oorho * maxmutilde*maxmutilde * ood2wall2;
    dCC[k][0] = - cw1 * dfw[k][0] * oorho * maxmutilde*maxmutilde * ood2wall2 - cw1 * fw * doorho[k][0] * maxmutilde*maxmutilde * ood2wall2 - cw1 * fw * oorho * 2.0*maxmutilde*dmaxmutilde[k][0] * ood2wall2;
    dCC[k][1] = - cw1 * dfw[k][1] * oorho * maxmutilde*maxmutilde * ood2wall2 - cw1 * fw * doorho[k][1] * maxmutilde*maxmutilde * ood2wall2 - cw1 * fw * oorho * 2.0*maxmutilde*dmaxmutilde[k][1] * ood2wall2;
    dCC[k][2] = - cw1 * dfw[k][2] * oorho * maxmutilde*maxmutilde * ood2wall2 - cw1 * fw * doorho[k][2] * maxmutilde*maxmutilde * ood2wall2 - cw1 * fw * oorho * 2.0*maxmutilde*dmaxmutilde[k][2] * ood2wall2;
    dCC[k][3] = - cw1 * dfw[k][3] * oorho * maxmutilde*maxmutilde * ood2wall2 - cw1 * fw * doorho[k][3] * maxmutilde*maxmutilde * ood2wall2 - cw1 * fw * oorho * 2.0*maxmutilde*dmaxmutilde[k][3] * ood2wall2;
    dCC[k][4] = - cw1 * dfw[k][4] * oorho * maxmutilde*maxmutilde * ood2wall2 - cw1 * fw * doorho[k][4] * maxmutilde*maxmutilde * ood2wall2 - cw1 * fw * oorho * 2.0*maxmutilde*dmaxmutilde[k][4] * ood2wall2;
    dCC[k][5] = - cw1 * dfw[k][5] * oorho * maxmutilde*maxmutilde * ood2wall2 - cw1 * fw * doorho[k][5] * maxmutilde*maxmutilde * ood2wall2 - cw1 * fw * oorho * 2.0*maxmutilde*dmaxmutilde[k][5] * ood2wall2;

    dSdU[k][5][0] = dAA[k][0] + dBB[k][0] + dCC[k][0];
    dSdU[k][5][1] = dAA[k][1] + dBB[k][1] + dCC[k][1];
    dSdU[k][5][2] = dAA[k][2] + dBB[k][2] + dCC[k][2];
    dSdU[k][5][3] = dAA[k][3] + dBB[k][3] + dCC[k][3];
    dSdU[k][5][4] = dAA[k][4] + dBB[k][4] + dCC[k][4];
    dSdU[k][5][5] = dAA[k][5] + dBB[k][5] + dCC[k][5];

  }
  
}

//------------------------------------------------------------------------------

#endif
