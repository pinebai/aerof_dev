#include <WallFcn.h>
#include <VarFcn.h>
#include <Vector3D.h>

#include <cmath>

// Included
#include <ViscoFcn.h>
#include <BcDef.h>

//------------------------------------------------------------------------------

// This function was modified to account for the derivative of the muw with respect to the conservative variables
// Included (MB*)
template<int neq>
void WallFcn::computeJacobianSurfaceTerm(int code, Vec3D &normal, double d2wall[3],
					 double *Vwall, double *V[3], 
					 double (*dRdU)[neq][neq])
{

  double delta, dT, rhow, muw;

  double drhow[3][neq];

  double dTdu0;
  double dTdu1;
  double dTdu2;
  double dTdu3;
  double dTdu4;

  double dutau[3][neq];

  double dadu0;
  double dadu1;
  double dadu2;
  double dadu3;
  double dadu4;

  double dmuw[3][neq];

  Vec3D ddudu0;
  Vec3D ddudu1;
  Vec3D ddudu2;
  Vec3D ddudu3;
  Vec3D ddudu4;

  Vec3D dtdu0;
  Vec3D dtdu1;
  Vec3D dtdu2;
  Vec3D dtdu3;
  Vec3D dtdu4;

  Vec3D du, uw;

  computeFaceValues(d2wall, Vwall, V, delta, du, dT, uw, rhow, muw);

  for (int i=0; i<3; ++i) {
    for (int j=0; j<neq; ++j) {
      drhow[i][j] = 0.0;
      dutau[i][j] = 0.0;
      dmuw[i][j] = 0.0;
    }
  }

  drhow[0][0] = third;
  drhow[1][0] = third;
  drhow[2][0] = third;

  double norm = sqrt(normal*normal);

  Vec3D n = (1.0/norm) * normal;
  Vec3D dndu0;
  Vec3D dndu1;
  Vec3D dndu2;
  Vec3D dndu3;
  Vec3D dndu4;

  dndu0[0] = 0.0;
  dndu0[1] = 0.0;
  dndu0[2] = 0.0;
  dndu1[0] = 0.0;
  dndu1[1] = 0.0;
  dndu1[2] = 0.0;
  dndu2[0] = 0.0;
  dndu2[1] = 0.0;
  dndu2[2] = 0.0;
  dndu3[0] = 0.0;
  dndu3[1] = 0.0;
  dndu3[2] = 0.0;
  dndu4[0] = 0.0;
  dndu4[1] = 0.0;
  dndu4[2] = 0.0;

  Vec3D t = computeTangentVector(n, du);

  double dut = du * t;

  double utau = computeFrictionVelocity(dut, delta, rhow, muw);
  
  double a = - rhow * utau*utau * norm;

  double Ttau;
  double dTtaudu0;
  double dTtaudu1;
  double dTtaudu2;
  double dTtaudu3;
  double dTtaudu4;
  
  double dMach = 0.0;
  
  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED) {
    Ttau = computeFrictionTemperature(utau, delta, rhow, dT, muw);
  }

  for (int k=0; k<3; ++k) {
  
    ddudu0[0] = - third * varFcn->getVelocityX(V[k]) / varFcn->getDensity(V[k]);
    ddudu0[1] = - third * varFcn->getVelocityY(V[k]) / varFcn->getDensity(V[k]);
    ddudu0[2] = - third * varFcn->getVelocityZ(V[k]) / varFcn->getDensity(V[k]);
    ddudu1[0] = third / varFcn->getDensity(V[k]);
    ddudu1[1] = 0.0;
    ddudu1[2] = 0.0;
    ddudu2[0] = 0.0;
    ddudu2[1] = third / varFcn->getDensity(V[k]);
    ddudu2[2] = 0.0;
    ddudu3[0] = 0.0;
    ddudu3[1] = 0.0;
    ddudu3[2] = third / varFcn->getDensity(V[k]);
    ddudu4[0] = 0.0;
    ddudu4[1] = 0.0;
    ddudu4[2] = 0.0;

    dtdu0 = computeDerivativeOfTangentVector(n, dndu0, du, ddudu0);
    dtdu1 = computeDerivativeOfTangentVector(n, dndu1, du, ddudu1);
    dtdu2 = computeDerivativeOfTangentVector(n, dndu2, du, ddudu2);
    dtdu3 = computeDerivativeOfTangentVector(n, dndu3, du, ddudu3);
    dtdu4 = computeDerivativeOfTangentVector(n, dndu4, du, ddudu4);

    dTdu0 = - third / V[k][0] * (varFcn->computeTemperature(V[k]) - 0.5 * (V[k][1]*V[k][1] + V[k][2]*V[k][2] + V[k][3]*V[k][3]));
    dTdu1 = - third / V[k][0] * V[k][1];
    dTdu2 = - third / V[k][0] * V[k][2];
    dTdu3 = - third / V[k][0] * V[k][3];
    dTdu4 = third / V[k][0];

    dmuw[k][0] = viscoFcn->compute_muDerivative(varFcn->computeTemperature(V[k]), dTdu0, dMach);
    dmuw[k][1] = viscoFcn->compute_muDerivative(varFcn->computeTemperature(V[k]), dTdu1, dMach);
    dmuw[k][2] = viscoFcn->compute_muDerivative(varFcn->computeTemperature(V[k]), dTdu2, dMach);
    dmuw[k][3] = viscoFcn->compute_muDerivative(varFcn->computeTemperature(V[k]), dTdu3, dMach); 
    dmuw[k][4] = viscoFcn->compute_muDerivative(varFcn->computeTemperature(V[k]), dTdu4, dMach); 

    dutau[k][0] = computeDerivativeOfFrictionVelocity(t, dtdu0, delta, rhow, drhow[k][0], du, ddudu0, muw, dmuw[k][0], dMach);
    dutau[k][1] = computeDerivativeOfFrictionVelocity(t, dtdu1, delta, rhow, drhow[k][1], du, ddudu1, muw, dmuw[k][1], dMach);
    dutau[k][2] = computeDerivativeOfFrictionVelocity(t, dtdu2, delta, rhow, drhow[k][2], du, ddudu2, muw, dmuw[k][2], dMach);
    dutau[k][3] = computeDerivativeOfFrictionVelocity(t, dtdu3, delta, rhow, drhow[k][3], du, ddudu3, muw, dmuw[k][3], dMach);
    dutau[k][4] = computeDerivativeOfFrictionVelocity(t, dtdu4, delta, rhow, drhow[k][4], du, ddudu4, muw, dmuw[k][4], dMach);
    
    dadu0 = - drhow[k][0] * utau*utau * norm - rhow * 2.0*utau*dutau[k][0] * norm;
    dadu1 = - drhow[k][1] * utau*utau * norm - rhow * 2.0*utau*dutau[k][1] * norm;
    dadu2 = - drhow[k][2] * utau*utau * norm - rhow * 2.0*utau*dutau[k][2] * norm;
    dadu3 = - drhow[k][3] * utau*utau * norm - rhow * 2.0*utau*dutau[k][3] * norm;
    dadu4 = - drhow[k][4] * utau*utau * norm - rhow * 2.0*utau*dutau[k][4] * norm;

    dRdU[k][0][0] = 0.0;
    dRdU[k][0][1] = 0.0;
    dRdU[k][0][2] = 0.0;
    dRdU[k][0][3] = 0.0;
    dRdU[k][0][4] = 0.0;

    dRdU[k][1][0] = dadu0 * t[0] + a * dtdu0[0];
    dRdU[k][1][1] = dadu1 * t[0] + a * dtdu1[0];
    dRdU[k][1][2] = dadu2 * t[0] + a * dtdu2[0];
    dRdU[k][1][3] = dadu3 * t[0] + a * dtdu3[0];
    dRdU[k][1][4] = dadu4 * t[0] + a * dtdu4[0];

    dRdU[k][2][0] = dadu0 * t[1] + a * dtdu0[1];
    dRdU[k][2][1] = dadu1 * t[1] + a * dtdu1[1];
    dRdU[k][2][2] = dadu2 * t[1] + a * dtdu2[1];
    dRdU[k][2][3] = dadu3 * t[1] + a * dtdu3[1];
    dRdU[k][2][4] = dadu4 * t[1] + a * dtdu4[1];

    dRdU[k][3][0] = dadu0 * t[2] + a * dtdu0[2];
    dRdU[k][3][1] = dadu1 * t[2] + a * dtdu1[2];
    dRdU[k][3][2] = dadu2 * t[2] + a * dtdu2[2];
    dRdU[k][3][3] = dadu3 * t[2] + a * dtdu3[2];
    dRdU[k][3][4] = dadu4 * t[2] + a * dtdu4[2];

    dRdU[k][4][0] = dadu0 * (uw * t) + a * (uw * dtdu0);
    dRdU[k][4][1] = dadu1 * (uw * t) + a * (uw * dtdu1);
    dRdU[k][4][2] = dadu2 * (uw * t) + a * (uw * dtdu2);
    dRdU[k][4][3] = dadu3 * (uw * t) + a * (uw * dtdu3);
    dRdU[k][4][4] = dadu4 * (uw * t) + a * (uw * dtdu4);

    if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED) {
      dTtaudu0 = computeDerivativeOfFrictionTemperature(utau, dutau[k][0], delta, rhow, drhow[k][0], dT, dTdu0, muw, dmuw[k][0], dMach);
      dTtaudu1 = computeDerivativeOfFrictionTemperature(utau, dutau[k][1], delta, rhow, drhow[k][1], dT, dTdu1, muw, dmuw[k][1], dMach);
      dTtaudu2 = computeDerivativeOfFrictionTemperature(utau, dutau[k][2], delta, rhow, drhow[k][2], dT, dTdu2, muw, dmuw[k][2], dMach);
      dTtaudu3 = computeDerivativeOfFrictionTemperature(utau, dutau[k][3], delta, rhow, drhow[k][3], dT, dTdu3, muw, dmuw[k][3], dMach);
      dTtaudu4 = computeDerivativeOfFrictionTemperature(utau, dutau[k][4], delta, rhow, drhow[k][4], dT, dTdu4, muw, dmuw[k][4], dMach);
      dRdU[k][4][0] += gam * drhow[k][0] * utau * Ttau * norm + gam * rhow * dutau[k][0] * Ttau * norm + gam * rhow * utau * dTtaudu0 * norm;
      dRdU[k][4][1] += gam * drhow[k][1] * utau * Ttau * norm + gam * rhow * dutau[k][1] * Ttau * norm + gam * rhow * utau * dTtaudu1 * norm;
      dRdU[k][4][2] += gam * drhow[k][2] * utau * Ttau * norm + gam * rhow * dutau[k][2] * Ttau * norm + gam * rhow * utau * dTtaudu2 * norm;
      dRdU[k][4][3] += gam * drhow[k][3] * utau * Ttau * norm + gam * rhow * dutau[k][3] * Ttau * norm + gam * rhow * utau * dTtaudu3 * norm;
      dRdU[k][4][4] += gam * drhow[k][4] * utau * Ttau * norm + gam * rhow * dutau[k][4] * Ttau * norm + gam * rhow * utau * dTtaudu4 * norm;
    }

  }

}

//------------------------------------------------------------------------------

// Original
/*
template<int neq>
void WallFcn::computeJacobianSurfaceTerm(int code, Vec3D &normal, double d2wall[3],
					 double *Vwall, double *V[3], 
					 double (*dRdU)[neq][neq])
{

  double delta, dT, rhow, muw;

  Vec3D du, uw;

  computeFaceValues(d2wall, Vwall, V, delta, du, dT, uw, rhow, muw);

  double oomuw = 1.0 / muw;

  double norm = sqrt(normal*normal);

  Vec3D n = (1.0/norm) * normal;

  Vec3D t = computeTangentVector(n, du);

  double utau = computeFrictionVelocity(t, delta, rhow, du, muw);
  
  double dplus = reynolds * utau * delta * rhow * oomuw;

  double f = 2.5 * log(1.0 + vkcst*dplus) 
    + 7.8 * (1.0 - exp(-eleventh*dplus) - eleventh*dplus*exp(-0.33*dplus));

  double df = 2.5*vkcst/(1.0 + vkcst*dplus) + 7.8*eleventh* 
    (exp(-eleventh*dplus) - exp(-0.33*dplus) + 0.33*dplus*exp(-0.33*dplus));

  double dfdutau = df * reynolds * delta * rhow * oomuw;

  double dfdrho = df * third * reynolds * utau * delta * oomuw;

  double oodFdutau = 1.0 / (f + utau * dfdutau); 

  for (int k=0; k<3; ++k) {
    double oorho3 = third / V[k][0];
    Vec3D u = this->varFcn->getVelocity(V[k]);

    double dFdU0 = oorho3 * (u*t) + utau * dfdrho;
    double dFdU1 = - oorho3 * t[0];
    double dFdU2 = - oorho3 * t[1];
    double dFdU3 = - oorho3 * t[2];

    double dutaudU0 = - dFdU0 * oodFdutau;
    double dutaudU1 = - dFdU1 * oodFdutau;
    double dutaudU2 = - dFdU2 * oodFdutau;
    double dutaudU3 = - dFdU3 * oodFdutau;
    
    double trut = 2.0 * rhow * utau;

    double drut2du0 = third * utau*utau + trut * dutaudU0;
    double drut2du1 = trut * dutaudU1;
    double drut2du2 = trut * dutaudU2;
    double drut2du3 = trut * dutaudU3;

//    t *= norm;

    dRdU[k][0][0] = 0.0;
    dRdU[k][0][1] = 0.0;
    dRdU[k][0][2] = 0.0;
    dRdU[k][0][3] = 0.0;
    dRdU[k][0][4] = 0.0;

    dRdU[k][1][0] = - drut2du0 * t[0] * norm;
    dRdU[k][1][1] = - drut2du1 * t[0] * norm;
    dRdU[k][1][2] = - drut2du2 * t[0] * norm;
    dRdU[k][1][3] = - drut2du3 * t[0] * norm;
    dRdU[k][1][4] = 0.0;

    dRdU[k][2][0] = - drut2du0 * t[1] * norm;
    dRdU[k][2][1] = - drut2du1 * t[1] * norm;
    dRdU[k][2][2] = - drut2du2 * t[1] * norm;
    dRdU[k][2][3] = - drut2du3 * t[1] * norm;
    dRdU[k][2][4] = 0.0;

    dRdU[k][3][0] = - drut2du0 * t[2] * norm;
    dRdU[k][3][1] = - drut2du1 * t[2] * norm;
    dRdU[k][3][2] = - drut2du2 * t[2] * norm;
    dRdU[k][3][3] = - drut2du3 * t[2] * norm;
    dRdU[k][3][4] = 0.0;

    double utw = uw * t * norm;

    dRdU[k][4][0] = - drut2du0 * utw;
    dRdU[k][4][1] = - drut2du1 * utw;
    dRdU[k][4][2] = - drut2du2 * utw;
    dRdU[k][4][3] = - drut2du3 * utw;
    dRdU[k][4][4] = 0.0;
  }

}
*/

//------------------------------------------------------------------------------

// This function was included to account for the derivative of the Vwall[5] with respect to the conservative variables
// Included
template<int neq>
void WallFcn::computeBCsJacobianWallValues(int code, Vec3D &normal, double d2wall[3],
					   double *Vwall, double *dVwall, double *V[3])
{

  double delta, dT, rhow, muw;

  double dTdu0;
  double dTdu1;
  double dTdu2;
  double dTdu3;
  double dTdu4;
  
  double dMach = 0.0;

  double drhow[3][neq];

  double dutau[3][neq];

  double dmuw[3][neq];

  double ddut[3][neq];

  double dvwall[3][neq];

  Vec3D ddudu0;
  Vec3D ddudu1;
  Vec3D ddudu2;
  Vec3D ddudu3;
  Vec3D ddudu4;

  Vec3D dtdu0;
  Vec3D dtdu1;
  Vec3D dtdu2;
  Vec3D dtdu3;
  Vec3D dtdu4;

  Vec3D du, uw;

  computeFaceValues(d2wall, Vwall, V, delta, du, dT, uw, rhow, muw);

  for (int i=0; i<3; ++i) {
    for (int j=0; j<neq; ++j) {
      drhow[i][j] = 0.0;
      dutau[i][j] = 0.0;
      dmuw[i][j] = 0.0;
      ddut[i][j] = 0.0;
      dvwall[i][j] = 0.0;
    }
  }

  drhow[0][0] = third;
  drhow[1][0] = third;
  drhow[2][0] = third;

  double norm = sqrt(normal*normal);

  Vec3D n = (1.0/norm) * normal;
  Vec3D dndu0;
  Vec3D dndu1;
  Vec3D dndu2;
  Vec3D dndu3;
  Vec3D dndu4;

  dndu0[0] = 0.0;
  dndu0[1] = 0.0;
  dndu0[2] = 0.0;
  dndu1[0] = 0.0;
  dndu1[1] = 0.0;
  dndu1[2] = 0.0;
  dndu2[0] = 0.0;
  dndu2[1] = 0.0;
  dndu2[2] = 0.0;
  dndu3[0] = 0.0;
  dndu3[1] = 0.0;
  dndu3[2] = 0.0;
  dndu4[0] = 0.0;
  dndu4[1] = 0.0;
  dndu4[2] = 0.0;

  Vec3D t = computeTangentVector(n, du);

  double dut = du * t;

  double utau = computeFrictionVelocity(dut, delta, rhow, muw);
  
  for (int k=0; k<3; ++k) {
  
    ddudu0[0] = - third * varFcn->getVelocityX(V[k]) / varFcn->getDensity(V[k]);
    ddudu0[1] = - third * varFcn->getVelocityY(V[k]) / varFcn->getDensity(V[k]);
    ddudu0[2] = - third * varFcn->getVelocityZ(V[k]) / varFcn->getDensity(V[k]);
    ddudu1[0] = third / varFcn->getDensity(V[k]);
    ddudu1[1] = 0.0;
    ddudu1[2] = 0.0;
    ddudu2[0] = 0.0;
    ddudu2[1] = third / varFcn->getDensity(V[k]);
    ddudu2[2] = 0.0;
    ddudu3[0] = 0.0;
    ddudu3[1] = 0.0;
    ddudu3[2] = third / varFcn->getDensity(V[k]);
    ddudu4[0] = 0.0;
    ddudu4[1] = 0.0;
    ddudu4[2] = 0.0;

    dtdu0 = computeDerivativeOfTangentVector(n, dndu0, du, ddudu0);
    dtdu1 = computeDerivativeOfTangentVector(n, dndu1, du, ddudu1);
    dtdu2 = computeDerivativeOfTangentVector(n, dndu2, du, ddudu2);
    dtdu3 = computeDerivativeOfTangentVector(n, dndu3, du, ddudu3);
    dtdu4 = computeDerivativeOfTangentVector(n, dndu4, du, ddudu4);

    dTdu0 = - third / V[k][0] * (varFcn->computeTemperature(V[k]) - 0.5 * (V[k][1]*V[k][1] + V[k][2]*V[k][2] + V[k][3]*V[k][3]));
    dTdu1 = - third / V[k][0] * V[k][1];
    dTdu2 = - third / V[k][0] * V[k][2];
    dTdu3 = - third / V[k][0] * V[k][3];
    dTdu4 = third / V[k][0];

    dmuw[k][0] = viscoFcn->compute_muDerivative(varFcn->computeTemperature(V[k]), dTdu0, dMach);
    dmuw[k][1] = viscoFcn->compute_muDerivative(varFcn->computeTemperature(V[k]), dTdu1, dMach);
    dmuw[k][2] = viscoFcn->compute_muDerivative(varFcn->computeTemperature(V[k]), dTdu2, dMach);
    dmuw[k][3] = viscoFcn->compute_muDerivative(varFcn->computeTemperature(V[k]), dTdu3, dMach); 
    dmuw[k][4] = viscoFcn->compute_muDerivative(varFcn->computeTemperature(V[k]), dTdu4, dMach); 

    dutau[k][0] = computeDerivativeOfFrictionVelocity(t, dtdu0, delta, rhow, drhow[k][0], du, ddudu0, muw, dmuw[k][0], dMach);
    dutau[k][1] = computeDerivativeOfFrictionVelocity(t, dtdu1, delta, rhow, drhow[k][1], du, ddudu1, muw, dmuw[k][1], dMach);
    dutau[k][2] = computeDerivativeOfFrictionVelocity(t, dtdu2, delta, rhow, drhow[k][2], du, ddudu2, muw, dmuw[k][2], dMach);
    dutau[k][3] = computeDerivativeOfFrictionVelocity(t, dtdu3, delta, rhow, drhow[k][3], du, ddudu3, muw, dmuw[k][3], dMach);
    dutau[k][4] = computeDerivativeOfFrictionVelocity(t, dtdu4, delta, rhow, drhow[k][4], du, ddudu4, muw, dmuw[k][4], dMach);

    ddut[k][0] = ddudu0 * t + du * dtdu0;
    ddut[k][1] = ddudu1 * t + du * dtdu1;
    ddut[k][2] = ddudu2 * t + du * dtdu2;
    ddut[k][3] = ddudu3 * t + du * dtdu3;
    ddut[k][4] = ddudu4 * t + du * dtdu4;

    for (int l=0; l<neq; ++l) {
      computeDerivativeOfWallValues(utau, dutau[k][l], delta, rhow, drhow[k][l], dut, ddut[k][l], muw, dmuw[k][l], 0.0, Vwall, dVwall);
      dvwall[k][l] = dVwall[5]; 
    }

  }

  for (int l=0; l<neq; ++l)
    dVwall[l] = dvwall[0][l] + dvwall[1][l] + dvwall[2][l]; 

}

//------------------------------------------------------------------------------
