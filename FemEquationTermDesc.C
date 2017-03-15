#include <FemEquationTermDesc.h>

#include <cstdlib>
#include <cstdio>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::max;
using std::min;
#endif

//------------------------------------------------------------------------------
//CHANGES_FOR_WATER
// as the Stokes' law for gas does not apply anymore, we have to distinguish
// between lambda and mu (the Lame coefficients)
//------------------------------------------------------------------------------

const double NavierStokesTerm::third = 1.0/3.0;
const double NavierStokesTerm::twothird = 2.0/3.0;
const double NavierStokesTerm::fourth = 1.0/4.0;

//------------------------------------------------------------------------------

FemEquationTermNS::FemEquationTermNS(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), FemEquationTerm(iod.volumes.volumeMap.dataMap)
{

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = new WallFcn(iod, varFcn, viscoFcn);

  velocity = iod.ref.rv.velocity;
  density = iod.ref.rv.density;
  length = iod.ref.rv.length; 

// Included (MB)
   if (iod.eqs.fluidModel.fluid == FluidModelData::PERFECT_GAS)
     completeJac = true;
   else
     completeJac = false;
}

//------------------------------------------------------------------------------

void FemEquationTermNS::computeTransportCoefficients(
  const double T, double &mu, double &lambda, double &kappa)
{

  mu     = viscoFcn->compute_mu(T);
  lambda = viscoFcn->compute_lambda(T,mu);
  kappa  = thermalCondFcn->compute(T);

}

void FemEquationTermNS::computeTransportCoefficientsPublic(
  const double T, double &mu, double &lambda, double &kappa)
{

  computeTransportCoefficients(T,mu,lambda,kappa);
  
  mu     *= ooreynolds_mu;
  lambda *= ooreynolds_mu;
  kappa  *= ooreynolds_mu;

}


//------------------------------------------------------------------------------

double FemEquationTermNS::computeViscousTimeStep(double X[3], double *V)
{
  double T;
  computeTemperature(V,T);
  double mul = ooreynolds_mu * viscoFcn->compute_mu(T);
  return mul/V[0];

}

//------------------------------------------------------------------------------

// Included (MB)
double FemEquationTermNS::computeDerivativeOfViscousTimeStep(double X[3], double dX[3], double *V, double *dV, double dMach)
{

  double T, dT;
  computeTemperature(V,T);
  computeDerivativeOfTemperature(V,dV,dT);
  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;
  double mul = ooreynolds_mu * viscoFcn->compute_mu(T);
  double dmul = dooreynolds_mu * viscoFcn->compute_mu(T) + ooreynolds_mu * viscoFcn->compute_muDerivative(T,dT,dMach);
  return (dmul*V[0]-mul*dV[0])/(V[0]*V[0]);

}

//------------------------------------------------------------------------------

bool FemEquationTermNS::computeVolumeTerm(double dp1dxj[4][3], double d2w[4], 
					  double *V[4], double *r, double *S, 
                                          double *PR, double tetVol, 
                                          SVec<double,3> &X, int nodeNum[4],  
                                          int material_id)
{

  bool porousmedia = false; 

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double dTdxj[3];
  computeTemperatureGradient(dp1dxj, T, dTdxj);

  double mu, lambda, kappa;
  computeTransportCoefficients(Tcg, mu, lambda, kappa);
  mu     *= ooreynolds_mu;
  lambda *= ooreynolds_mu;
  kappa  *= ooreynolds_mu;

  double (*R)[5] = reinterpret_cast<double (*)[5]>(r);
  computeVolumeTermNS(mu, lambda, kappa, ucg, dudxj, dTdxj, R);
  
  // Initialize PR (porous media term)
  for (int j=0; j<3*4; ++j) PR[j] = 0.0; 

  if(material_id > 0)
  {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);

	  if(it != volInfo.end()) 
	  {
		  // if porous media with material_id has been defined in the input file
       porousmedia = computeVolumeTermPorousCore(tetVol, it, length, density, velocity, ucg, V, PR);
    }
  }
 
  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;

  return (porousmedia);

}

//------------------------------------------------------------------------------

// Included (MB)
bool FemEquationTermNS::computeDerivativeOfVolumeTerm(double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4],
																		double *V[4], double *dV[4], double dMach, double *dr, 
																		double *dS, double *dPR, double dtetVol, SVec<double,3> &X,
                                          int nodeNum[4], int material_id)
{

  bool porousmedia = false; 

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double du[4][3], ducg[3];
  computeDerivativeOfVelocity(dV, du, ducg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dT[4], dTcg;
  computeDerivativeOfTemperature(V, dV, dT, dTcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double ddudxj[3][3];
  computeDerivativeOfVelocityGradient(dp1dxj, ddp1dxj, u, du, ddudxj);

  double dTdxj[3];
  computeTemperatureGradient(dp1dxj, T, dTdxj);

  double ddTdxj[3];
  computeDerivativeOfTemperatureGradient(dp1dxj, ddp1dxj, T, dT, ddTdxj);

  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;

  double mu, lambda, kappa;
  computeTransportCoefficients(Tcg, mu, lambda, kappa);

  double dmu     = dooreynolds_mu * mu + ooreynolds_mu * viscoFcn->compute_muDerivative(Tcg, dTcg, dMach);
  double dlambda = dooreynolds_mu * lambda + ooreynolds_mu * viscoFcn->compute_lambdaDerivative(mu, dmu, dMach);
  double dkappa  = dooreynolds_mu * kappa + ooreynolds_mu * thermalCondFcn->computeDerivative(Tcg, dTcg, dMach);

  mu     *= ooreynolds_mu;
  lambda *= ooreynolds_mu;
  kappa  *= ooreynolds_mu;

  double (*dR)[5] = reinterpret_cast<double (*)[5]>(dr);
  computeDerivativeOfVolumeTermNS(mu, dmu, lambda, dlambda, kappa, dkappa, ucg, ducg, dudxj, ddudxj, dTdxj, ddTdxj, dR);

  // Initialize PR (porous media term)
  for (int j=0; j<3*4; ++j) dPR[j] = 0.0; 

  if(material_id > 0)
  {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
	  if(it!=  volInfo.end()) 
	  { 
		  // if porous media with material_id has been defined in the input file
       porousmedia = true;
       fprintf(stderr, "***** Inside the file FemEquationTermDesc.C the derivative related to porus media is not implemented *****\n");
       exit(1);
    } 
  }
 
  dS[0] = 0.0;
  dS[1] = 0.0;
  dS[2] = 0.0;
  dS[3] = 0.0;
  dS[4] = 0.0;

  return (porousmedia);

}


// Included (YC)
void FemEquationTermNS::computeDerivativeOperatorsOfVolumeTerm(double dp1dxj[4][3], double *V[4],
            double (*drddp1dxj)[5][4][3], double (*drdV)[5][4][5], double (*drdMach)[5])
{
  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double dudV[4][3][4][4] = {0}, ducgdV[3][4][4] = {0};
  computeDerivativeOperatorsOfVelocity(dudV, ducgdV);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dTdV[4][5] = {0}, dTcgdV[4] = {0};
  computeDerivativeOperatorsOfTemperature(V, dTdV, dTcgdV);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double ddudxj[3][3] = {0};
  double ddudxjddp1dxj[3][3][4][3] = {0}, ddudxjdu[3][3][4][3] = {0};
  computeDerivativeOperatorsOfVelocityGradient(dp1dxj, u, ddudxjddp1dxj, ddudxjdu);

  double dTdxj[3];
  computeTemperatureGradient(dp1dxj, T, dTdxj);

  double ddTdxjddp1dxj[3][4][3] = {0}, ddTdxjdT[3][4] = {0};
  computeDerivativeOperatorsOfTemperatureGradient(dp1dxj, T, ddTdxjddp1dxj, ddTdxjdT);

  double coef = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS;
  double mu, lambda, kappa;
  computeTransportCoefficients(Tcg, mu, lambda, kappa);

  double dmudTcg(0), dmudMach(0);
  viscoFcn->compute_muDerivativeOperators(Tcg, dmudTcg, dmudMach);
  double dlambdadmu(0), dlambdadMach(0);
  viscoFcn->compute_lambdaDerivativeOperators(dlambdadmu, dlambdadMach);
  double dkappadTcg(0), dkappadMach(0);
  thermalCondFcn->computeDerivativeOperators(Tcg, dkappadTcg, dkappadMach);

  mu     *= ooreynolds_mu;
  lambda *= ooreynolds_mu;
  kappa  *= ooreynolds_mu;
  double drdmu[3][5] = {0}, drdlambda[3][5] = {0}, drdkappa[3][5] = {0}, drdu[3][5][3] = {0}, drddudxj[3][5][3][3] = {0}, drddTdxj[3][5][3] = {0};
  computeDerivativeOperatorsOfVolumeTermNS(mu, lambda, kappa, ucg, dudxj, dTdxj,
                                           drdmu, drdlambda, drdkappa, drdu, drddudxj, drddTdxj);
  for(int i=0; i<3; ++i) {
    for(int j=0; j<5; ++j) {
      drdMach[i][j] += (drdkappa[i][j]*(coef*kappa+ooreynolds_mu*dkappadMach) + drdmu[i][j]*(coef*mu+ooreynolds_mu*dmudMach) + drdlambda[i][j]*(coef*lambda+ooreynolds_mu*dlambdadMach+ooreynolds_mu*dlambdadmu*ooreynolds_mu*dmudMach+ooreynolds_mu*dlambdadmu*mu*coef));
      for(int k=0; k<4; ++k) {
        for(int l=0; l<5; ++l) {
          drdV[i][j][k][l] += (drdmu[i][j]*ooreynolds_mu*dmudTcg + drdkappa[i][j]*ooreynolds_mu*dkappadTcg + drdlambda[i][j]*ooreynolds_mu*dlambdadmu*ooreynolds_mu*dmudTcg)*dTcgdV[k];
        }
      }
      for(int k=0; k<3; ++k) {
        for(int l=0; l<4; ++l) {
          for(int m=0; m<5; ++m) {
            drdV[i][j][l][m] += drddTdxj[i][j][k]*ddTdxjdT[k][l]*dTdV[l][m];//added
          }
          for(int w=0; w<3; ++w) {
            drddp1dxj[i][j][l][w] += drddTdxj[i][j][k]*ddTdxjddp1dxj[k][l][w];
          }
        }
        for(int l=0; l<4; ++l) {
          for(int m=0; m<4; ++m) {
            drdV[i][j][l][m] += drdu[i][j][k]*ducgdV[k][l][m];
          }
        }
        for(int l=0; l<3; ++l) {
          for(int m=0; m<4; ++m) {
            for(int n=0; n<3; ++n) {
              drddp1dxj[i][j][m][n] += drddudxj[i][j][k][l]*ddudxjddp1dxj[k][l][m][n];
              for(int o=0; o<4; ++o) {
                for(int p=0; p<4; ++p) {
                  drdV[i][j][o][p] += drddudxj[i][j][k][l]*ddudxjdu[k][l][m][n]*dudV[m][n][o][p];
                }
              }
            }
          }
        }
      }
    }
  }
}



//------------------------------------------------------------------------------

// This function was modified to account for the derivative of the mu with respect to the conservative variables
// Included (MB*)
bool FemEquationTermNS::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						  double *V[4], double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                  SVec<double,3> &X, int nodeNum[4], int material_id)
{

  bool porousmedia = false;
  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double mu, lambda, kappa;
  computeTransportCoefficients(Tcg, mu, lambda, kappa);
  mu     *= ooreynolds_mu;
  lambda *= ooreynolds_mu;
  kappa  *= ooreynolds_mu;

  //fprintf(stdout, "mu = %e, lambda = %e, kappa = %e\n", mu, lambda, kappa);
  //mu = 2.672612e-07;
  //lambda = -1.781742e-07;
  //kappa = 5.196746e-07;

  double (*dRdU)[3][5][5] = reinterpret_cast<double (*)[3][5][5]>(drdu);
  double (*dPdU)[4][5][5] = reinterpret_cast<double (*)[4][5][5]>(dpdu);

  double dmu[4][5];

  double dlambda[4][5];
  
  double dkappa[4][5];

  double dTcgdu0 = 0.0;
  double dTcgdu1 = 0.0;
  double dTcgdu2 = 0.0;
  double dTcgdu3 = 0.0;
  double dTcgdu4 = 0.0;
  
  double dMach = 0.0;

  for (int k=0; k<4*3*5*5; ++k)
    drdu[k] = 0.0;
  for (int k=0; k<4*5*5; ++k)
    dsdu[k] = 0.0;
  for (int k=0; k<4*4*5*5; ++k)
    dpdu[k] = 0.0;

  for (int k=0; k<4; ++k) {

    if (completeJac) { 
      dTcgdu0 = - 0.25 / V[k][0] * (T[k] - 0.5 * (V[k][1]*V[k][1] + V[k][2]*V[k][2] + V[k][3]*V[k][3]));
      dTcgdu1 = - 0.25 / V[k][0] * V[k][1];
      dTcgdu2 = - 0.25 / V[k][0] * V[k][2];
      dTcgdu3 = - 0.25 / V[k][0] * V[k][3];
      dTcgdu4 = 0.25 / V[k][0];
    }

    dkappa[k][0] = ooreynolds_mu * thermalCondFcn->computeDerivative(Tcg, dTcgdu0, dMach);
    dkappa[k][1] = ooreynolds_mu * thermalCondFcn->computeDerivative(Tcg, dTcgdu1, dMach);
    dkappa[k][2] = ooreynolds_mu * thermalCondFcn->computeDerivative(Tcg, dTcgdu2, dMach);
    dkappa[k][3] = ooreynolds_mu * thermalCondFcn->computeDerivative(Tcg, dTcgdu3, dMach); 
    dkappa[k][4] = ooreynolds_mu * thermalCondFcn->computeDerivative(Tcg, dTcgdu4, dMach); 

    dmu[k][0] = viscoFcn->compute_muDerivative(Tcg, dTcgdu0, dMach);
    dmu[k][1] = viscoFcn->compute_muDerivative(Tcg, dTcgdu1, dMach);
    dmu[k][2] = viscoFcn->compute_muDerivative(Tcg, dTcgdu2, dMach);
    dmu[k][3] = viscoFcn->compute_muDerivative(Tcg, dTcgdu3, dMach); 
    dmu[k][4] = viscoFcn->compute_muDerivative(Tcg, dTcgdu4, dMach); 

    for(int i=0; i<5; ++i){
      dlambda[k][i] = ooreynolds_mu * viscoFcn->compute_lambdaDerivative(mu/ooreynolds_mu, dmu[k][i], dMach);
      dmu[k][i] *= ooreynolds_mu;
    }
    //fprintf(stdout, "node = %d -- dmu     = (%e %e %e %e %e)\n", k, dmu[k][0], dmu[k][1], dmu[k][2], dmu[k][3], dmu[k][4]);
    //fprintf(stdout, "node = %d -- dlambda = (%e %e %e %e %e)\n", k, dlambda[k][0], dlambda[k][1], dlambda[k][2], dlambda[k][3], dlambda[k][4]);
    //fprintf(stdout, "node = %d -- dkappa  = (%e %e %e %e %e)\n", k, dkappa[k][0], dkappa[k][1], dkappa[k][2], dkappa[k][3], dkappa[k][4]);

    //for(int toto=0; toto<5; ++toto){
    //  dmu[k][toto] = 0.0;
    //  dlambda[k][toto] = 0.0;
    //  dkappa[k][toto] = 0.0;
    //}
  }

  //fprintf(stdout, "dp1dxj = %e %e %e %e %e %e %e %e %e %e %e %e\n", dp1dxj[0][0], dp1dxj[0][1], dp1dxj[0][2], dp1dxj[1][0], dp1dxj[1][1], dp1dxj[1][2], dp1dxj[2][0], dp1dxj[2][1], dp1dxj[2][2], dp1dxj[3][0], dp1dxj[3][1], dp1dxj[3][2]); 
  //T[0] *= 1.1; T[1] *= 1.2; T[2] *= 0.9; T[3] *= 0.98;
  //fprintf(stdout, " T = %e %e %e %e\n", T[0], T[1], T[2], T[3]);
  //for(int toto=0; toto<5; ++toto){
    //V[0][toto] *= 1.1;
    //V[1][toto] *= 1.2;
    //V[2][toto] *= 0.9;
    //V[3][toto] *= 0.98;
  //}
  //for(int toto=0; toto<4; ++toto)
  //fprintf(stdout, " V[%d] = %e %e %e %e %e\n", toto, V[toto][0], V[toto][1], V[toto][2], V[toto][3], V[toto][4]);
  computeJacobianVolumeTermNS(dp1dxj, mu, dmu, lambda, dlambda, kappa, dkappa, V, T, dRdU);
  //fprintf(stdout, "dRdU = \n");
  //for(int toto=0; toto<4*3*5*5; ++toto) fprintf(stdout, "   %e\n", drdu[toto]);
  //fflush(stdout);
  //exit(1);

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
      double u[4][3], ucg[3];
      computeVelocity(V, u, ucg);
      porousmedia = computeJacobianVolumeTermPorousCore(tetVol, it, length, density, velocity, ucg, V, dPdU);
    }
  }

  return (porousmedia);

}

//------------------------------------------------------------------------------

void FemEquationTermNS::computeSurfaceTerm(int code, Vec3D &n, double d2w[3], 
					   double *Vwall, double *V[3], double *R)
{
  wallFcn->computeSurfaceTerm(code, n, d2w, Vwall, V, R);
}

//------------------------------------------------------------------------------

double FemEquationTermNS::computeNormDerivWallFcn(double rho, double T, double Du1, 
																double DT1, double d2w, 
																double &dudn, double &dTdn)
{
	double ut = wallFcn->computedudT(rho, T, Du1, DT1, d2w, dudn, dTdn);

	return ut;
}

//------------------------------------------------------------------------------

// Included (MB)
void FemEquationTermNS::computeDerivativeOfSurfaceTerm(int code, Vec3D &n, Vec3D &dn, double d2w[3],
					   double *Vwall, double *dVwall, double *V[3], double *dV[3], double dMach, double *dR)
{

  wallFcn->computeDerivativeOfSurfaceTerm(code, n, dn, d2w, Vwall, dVwall, V, dV, dMach, dR);

}

//------------------------------------------------------------------------------

void FemEquationTermNS::computeJacobianSurfaceTerm(int code, Vec3D &n, 
						   double d2w[3], double *Vwall, 
						   double *V[3], double *drdu)
{

  for (int k=0; k<3*5*5; ++k)
    drdu[k] = 0.0;

  double (*dRdU)[5][5] = reinterpret_cast<double (*)[5][5]>(drdu);
  wallFcn->computeJacobianSurfaceTerm(code, n, d2w, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

void FemEquationTermNS::computeSurfaceTerm(double dp1dxj[4][3], int code,
					   Vec3D &n, double d2w[4],
					   double *Vwall, double *V[4], double *R)
{
  
  computeSurfaceTermNS(dp1dxj, n, Vwall, V, R);

}

//------------------------------------------------------------------------------

// Included (MB)
void FemEquationTermNS::computeDerivativeOfSurfaceTerm(double dp1dxj[4][3], double ddp1dxj[4][3], int code,
					   Vec3D &n, Vec3D &dn, double d2w[4],
					   double *Vwall, double *dVwall, double *V[4], double *dV[4], double dMach, double *dR)
{

  computeDerivativeOfSurfaceTermNS(dp1dxj, ddp1dxj, n, dn, Vwall, dVwall, V, dV, dMach, dR);

}

//------------------------------------------------------------------------------

void FemEquationTermNS::computeJacobianSurfaceTerm(double dp1dxj[4][3], int code,
						   Vec3D &n, double d2w[4], 
						   double *Vwall, double *V[4], 
						   double *drdu)
{

  double (*dRdU)[5][5] = reinterpret_cast<double (*)[5][5]>(drdu);
  computeJacobianSurfaceTermNS(dp1dxj, n, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

FemEquationTermSA::FemEquationTermSA(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), SATerm(iod), FemEquationTerm(iod.volumes.volumeMap.dataMap),
  turbThermalCondFcn(iod, viscoFcn, vf)
{

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = new WallFcnSA(iod, varFcn, viscoFcn);

  x0 =  iod.eqs.tc.tr.bfix.x0;
  x1 =  iod.eqs.tc.tr.bfix.x1;
  y0 =  iod.eqs.tc.tr.bfix.y0;
  y1 =  iod.eqs.tc.tr.bfix.y1;
  z0 =  iod.eqs.tc.tr.bfix.z0;
  z1 =  iod.eqs.tc.tr.bfix.z1;
  
  if (x0>x1 || y0>y1 || z0>z1) trip = 0;
  else   trip = 1;

  if (iod.ts.implicit.tmcoupling == ImplicitData::STRONG && trip==1) { 
    fprintf(stderr,"** Warning: Laminar-turbulent trip not implemented for Strongly Coupled NS-SA simulation \n");
    trip = 0;
  }

  if (iod.eqs.tc.tm.sa.form == SAModelData::FV3)
    usefv3 = true;
  else
    usefv3 = false;

  velocity = iod.ref.rv.velocity;
  density = iod.ref.rv.density;
  length = iod.ref.rv.length;

// Included (MB)
   if (iod.eqs.fluidModel.fluid == FluidModelData::PERFECT_GAS)
     completeJac = true;
   else
     completeJac = false;

}

//------------------------------------------------------------------------------

void FemEquationTermSA::computeLaminarTransportCoefficients(
  const double T, double &mul, double &lambdal, double &kappal)
{

  mul     = viscoFcn->compute_mu(T);
  lambdal = viscoFcn->compute_lambda(T,mul);
  kappal  = thermalCondFcn->compute(T);

}

//------------------------------------------------------------------------------

void FemEquationTermSA::computeTurbulentTransportCoefficients(
  double *V[], int nodeNum[], SVec<double,3> &X,
  const double mul, const double lambdal, const double kappal,
  double &mutilde, double &mut, double &lambdat, double &kappat)
{

  mut = computeTurbulentViscosity(V, mul, mutilde);

  //Applying the laminar-turbulent trip
  if(trip) {
    int in_trip = 0;
    for (int k=0; k<3; k++)
      if (X[nodeNum[k]][0]>=x0 && X[nodeNum[k]][0]<=x1 &&
          X[nodeNum[k]][1]>=y0 && X[nodeNum[k]][1]<=y1 &&
	  X[nodeNum[k]][2]>=z0 && X[nodeNum[k]][2]<=z1 )
        in_trip = 1;

    if (!in_trip) mut = 0.0;
  }
  lambdat = computeSecondTurbulentViscosity(lambdal, mul, mut);
  kappat  = turbThermalCondFcn.turbulentConductivity(mut);

}

//------------------------------------------------------------------------------

double FemEquationTermSA::computeViscousTimeStep(double X[3], double *V)
{
  double T;
  computeTemperature(V,T);
  double mul = viscoFcn->compute_mu(T);
  double mut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1)
       mut = computeTurbulentViscosity(V, mul);
    else
       mut = 0.0;
  }
  else
    mut = computeTurbulentViscosity(V, mul);

  return ooreynolds_mu * (mul+mut)/V[0];

}

//------------------------------------------------------------------------------

// Included (MB)
double FemEquationTermSA::computeDerivativeOfViscousTimeStep(double X[3], double dX[3], double *V, double *dV, double dMach)
{
  double T,dT;
  computeTemperature(V,T);
  computeDerivativeOfTemperature(V,dV,dT);
  double mul = viscoFcn->compute_mu(T);
  double dmul = viscoFcn->compute_muDerivative(T,dT,dMach);
  double mut, dmut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1) {
       mut = computeTurbulentViscosity(V, mul);
       dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
    }
    else {
       mut = 0.0;
       dmut = 0.0;
    }
  }
  else {
    mut = computeTurbulentViscosity(V, mul);
    dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
  }

  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;

  return dooreynolds_mu * (mul+mut)/V[0] + ooreynolds_mu * ((dmul+dmut)*V[0]-(mul+mut)*dV[0])/(V[0]*V[0]);

}

//------------------------------------------------------------------------------

bool FemEquationTermSA::computeVolumeTerm(double dp1dxj[4][3], double d2w[4], 
					  double *V[4], double *r, double *S, 
                                          double *PR, double tetVol,
                                          SVec<double,3> &X,
                                          int nodeNum[4], int material_id)
{

  bool porousmedia = false;

  const double sixth = 1.0/6.0;

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double dTdxj[3];
  computeTemperatureGradient(dp1dxj, T, dTdxj);

  double mul, lambdal, kappal;
  computeLaminarTransportCoefficients(Tcg, mul, lambdal, kappal);

  double mutilde;
  double mut, lambdat, kappat;
  computeTurbulentTransportCoefficients(V, nodeNum, X, mul, lambdal, kappal, mutilde, mut, lambdat, kappat);

  double (*R)[6] = reinterpret_cast<double (*)[6]>(r);

  double absmutilde = fabs(mutilde);
  double maxmutilde = max(mutilde, 0.0);
  double mu5 = oosigma * (mul + absmutilde);

  double dnutildedx = dp1dxj[0][0]*V[0][5] 
	                 + dp1dxj[1][0]*V[1][5] 
	                 + dp1dxj[2][0]*V[2][5] 
	                 + dp1dxj[3][0]*V[3][5];

  double dnutildedy = dp1dxj[0][1]*V[0][5] 
	                 + dp1dxj[1][1]*V[1][5] 
                 	  + dp1dxj[2][1]*V[2][5] 
                 	  + dp1dxj[3][1]*V[3][5];

  double dnutildedz = dp1dxj[0][2]*V[0][5] 
	                 + dp1dxj[1][2]*V[1][5] 
                    + dp1dxj[2][2]*V[2][5] 
	                 + dp1dxj[3][2]*V[3][5];

  R[0][5] = mu5 * dnutildedx;
  R[1][5] = mu5 * dnutildedy;
  R[2][5] = mu5 * dnutildedz;

  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;

  double d2wall = 0.25 * (d2w[0] + d2w[1] + d2w[2] + d2w[3]);

  if(d2wall >= 1.e-15) 
  {
    double chi = max(mutilde/mul, 0.001);
    double chi3 = chi*chi*chi;
    double fv1 = chi3 / (chi3 + cv1_pow3);
    double fv2  = 1.-chi/(1.+chi*fv1);
    double fv3  = 1.0;

    if(usefv3) 
	 {
      fv2 = 1.0 + oocv2*chi;
      fv2 = 1.0 / (fv2*fv2*fv2);
      fv3 = (1.0 + chi*fv1) * (1.0 - fv2) / chi;
    }
    double ood2wall2 = 1.0 / (d2wall * d2wall);
    double rho = 0.25 * (V[0][0] + V[1][0] + V[2][0] + V[3][0]);
    double oorho = 1.0 / rho;
    double zz = ooreynolds_mu * oovkcst2 * maxmutilde * oorho * ood2wall2;
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

    double AA = oosigma * cb2 * rho * (dnutildedx*dnutildedx + dnutildedy*dnutildedy + dnutildedz*dnutildedz);
    double BB = cb1 * Stilde * absmutilde;
    double CC = - cw1 * fw * oorho   * maxmutilde*maxmutilde * ood2wall2;
    //double CC = - cw1 * fw * oorho * oorho * maxmutilde*maxmutilde * ood2wall2;
    S[5] = AA + BB + CC;
  }
  else 
  {
    S[5] = 0.0;
  }
  
  // Initialize PR (porous media term)
  for (int j=0; j<3*4; ++j) PR[j] = 0.0; 

  if(material_id>0) 
  {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
	  if(it != volInfo.end()) 
	  {
		  // if porous media with material_id has been defined in the input file
      porousmedia = true;

      mut     = computePorousTurbulentViscosity(it, ucg, length);
      lambdat = computeSecondPorousTurbulentViscosity(lambdal, mul, mut);
      kappat  = turbThermalCondFcn.turbulentConductivity(mut);

      porousmedia = computeVolumeTermPorousCore(tetVol, it, length, density, velocity, ucg, V, PR);
    }
  }

  // In all cases, porous media or not, compute viscous term for Navier-Stokes equations
  double mu, lambda, kappa;
  mu     = ooreynolds_mu * (mul + mut);
  lambda = ooreynolds_mu * (lambdal + lambdat);
  kappa  = ooreynolds_mu * (kappal + kappat);
  computeVolumeTermNS(mu, lambda, kappa, ucg, dudxj, dTdxj, R);

  return (porousmedia);

}

void FemEquationTermSA::computeSourceTerm(double dudxj[3][3],double dnudx[3],
					  double d2wall, 
					  double *V,  double *S)
{

  bool porousmedia = false;

  const double sixth = 1.0/6.0;

  double Tcg = varFcn->computeTemperature(V);

  double mul, lambdal, kappal;
  computeLaminarTransportCoefficients(Tcg, mul, lambdal, kappal);

  double mutilde;
  double mut, lambdat, kappat;

  mutilde = V[0]*V[5];

  double chi = mutilde / mul;
  double chi3 = chi*chi*chi;
  double fv1 = chi3 / (chi3 + cv1_pow3);

  mut = mutilde*fv1;

  lambdat = computeSecondTurbulentViscosity(lambdal, mul, mut);
  kappat  = turbThermalCondFcn.turbulentConductivity(mut);

  double absmutilde = fabs(mutilde);
  double maxmutilde = max(mutilde, 0.0);
  double mu5 = oosigma * (mul + absmutilde);

  double dnutildedx = dnudx[0];
  double dnutildedy = dnudx[1];
  double dnutildedz = dnudx[2];
  
  /*double dnutildedx = dp1dxj[0][0]*V[0][5] + dp1dxj[1][0]*V[1][5] + 
    dp1dxj[2][0]*V[2][5] + dp1dxj[3][0]*V[3][5];
  double dnutildedy = dp1dxj[0][1]*V[0][5] + dp1dxj[1][1]*V[1][5] + 
    dp1dxj[2][1]*V[2][5] + dp1dxj[3][1]*V[3][5];
  double dnutildedz = dp1dxj[0][2]*V[0][5] + dp1dxj[1][2]*V[1][5] + 
    dp1dxj[2][2]*V[2][5] + dp1dxj[3][2]*V[3][5];

  R[0][5] = mu5 * dnutildedx;
  R[1][5] = mu5 * dnutildedy;
  R[2][5] = mu5 * dnutildedz;
  */
  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;


  if  (d2wall >= 1.e-15) {
    chi = max(mutilde/mul, 0.001);
    chi3 = chi*chi*chi;
    fv1 = chi3 / (chi3 + cv1_pow3);
    double fv2  = 1.-chi/(1.+chi*fv1);
    double fv3  = 1.0;
    if (usefv3) {
      fv2 = 1.0 + oocv2*chi;
      fv2 = 1.0 / (fv2*fv2*fv2);
      fv3 = (1.0 + chi*fv1) * (1.0 - fv2) / chi;
    }
    double ood2wall2 = 1.0 / (d2wall * d2wall);
    double rho = V[0];
    double oorho = 1.0 / rho;
    double zz = ooreynolds_mu * oovkcst2 * maxmutilde * oorho * ood2wall2;
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

    double AA = oosigma * cb2 * rho * 
      (dnutildedx*dnutildedx + dnutildedy*dnutildedy + dnutildedz*dnutildedz);
    double BB = cb1 * Stilde * absmutilde;
    double CC = - cw1 * fw * oorho   * maxmutilde*maxmutilde * ood2wall2;
    //double CC = - cw1 * fw * oorho * oorho * maxmutilde*maxmutilde * ood2wall2;
    S[5] = AA + BB + CC;
  }
  else {
    S[5] = 0.0;
  }
}

//------------------------------------------------------------------------------

// Included (MB)
bool FemEquationTermSA::computeDerivativeOfVolumeTerm(double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4],
					  double *V[4], double *dV[4], double dMach, double *dr, double *dS, double *dPR, double dtetVol, SVec<double,3> &X,
                                          int nodeNum[4], int material_id)
{

  bool porousmedia = false;

  const double sixth = 1.0/6.0;

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double du[4][3], ducg[3];
  computeDerivativeOfVelocity(dV, du, ducg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dT[4], dTcg;
  computeDerivativeOfTemperature(V, dV, dT, dTcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double ddudxj[3][3];
  computeDerivativeOfVelocityGradient(dp1dxj, ddp1dxj, u, du, ddudxj);

  double dTdxj[3];
  computeTemperatureGradient(dp1dxj, T, dTdxj);

  double ddTdxj[3];
  computeDerivativeOfTemperatureGradient(dp1dxj, ddp1dxj, T, dT, ddTdxj);

  double mul  = viscoFcn->compute_mu(Tcg);
  double dmul = viscoFcn->compute_muDerivative(Tcg, dTcg, dMach);
  double lambdal  = viscoFcn->compute_lambda(Tcg, mul);
  double dlambdal = viscoFcn->compute_lambdaDerivative(mul, dmul, dMach);
  double kappal  = thermalCondFcn->compute(Tcg);
  double dkappal = thermalCondFcn->computeDerivative(Tcg, dTcg, dMach);

// Test
  double mutilde = 0.0;
  double dmutilde = 0.0;

  double mut, lambdat;
  double dmut, dlambdat;
  
  // Applying the laminar-turbulent trip
  if(trip){
    if((X[nodeNum[0]][0]>=x0 && X[nodeNum[0]][0]<=x1 && X[nodeNum[0]][1]>=y0 && X[nodeNum[0]][1]<=y1 &&
    	X[nodeNum[0]][2]>=z0 && X[nodeNum[0]][2]<=z1) || (X[nodeNum[1]][0]>=x0 && X[nodeNum[1]][0]<=x1 &&
    	X[nodeNum[1]][1]>=y0 && X[nodeNum[1]][1]<=y1 && X[nodeNum[1]][2]>=z0 && X[nodeNum[1]][2]<=z1) ||
    	(X[nodeNum[2]][0]>=x0 && X[nodeNum[2]][0]<=x1 && X[nodeNum[2]][1]>=y0 && X[nodeNum[2]][1]<=y1 &&
    	X[nodeNum[2]][2]>=z0 && X[nodeNum[2]][2]<=z1) || (X[nodeNum[3]][0]>=x0 && X[nodeNum[3]][0]<=x1 &&
    	X[nodeNum[3]][1]>=y0 && X[nodeNum[3]][1]<=y1 && X[nodeNum[3]][2]>=z0 && X[nodeNum[3]][2]<=z1)) {
    	mut = computeTurbulentViscosity(V, mul, mutilde);
    	dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul, dmutilde);
    }
    else {
        computeTurbulentViscosity(V, mul, mutilde);
        computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul, dmutilde);
    	mut = 0.0; 
    	dmut = 0.0;
    }
  }
  else {
    mut = computeTurbulentViscosity(V, mul, mutilde);
    dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul, dmutilde);
  }

  lambdat  = computeSecondTurbulentViscosity(lambdal, mul, mut);
  dlambdat = computeDerivativeOfSecondTurbulentViscosity(lambdal, dlambdal, mul, dmul, mut, dmut);

  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;
  
  double mu;
  double dmu;
  double lambda;
  double dlambda;
  double kappa;
  double dkappa;

  double (*dR)[6] = reinterpret_cast<double (*)[6]>(dr);

  double absmutilde = fabs(mutilde);
  double dabsmutilde;
  if  (mutilde != 0.0)
    dabsmutilde = ( fabs(mutilde) / mutilde ) * dmutilde;
  else {
    fprintf(stderr, "***** Inside the file FemEquationTermDesc.C the varible mutilde is zero *****\n");
    //exit(1);
    dabsmutilde = 0.0;
  }

  if (mutilde == 0.0) {
    fprintf(stderr, "***** Inside the file FemEquationTermDesc.C the varibles in the function max are equal *****\n");
    //exit(1);
  }
  double maxmutilde = max(mutilde, 0.0);
  double dmaxmutilde;
  if ( maxmutilde == 0.0 )  {
    dmaxmutilde = 0.0;
  }
  else {
    dmaxmutilde = dmutilde;
  }

  double mu5 = oosigma * (mul + absmutilde);

  // These values can be non-zero.
  double d_oosigma = (SATerm::oosigma / NavierStokesTerm::ooreynolds) * dooreynolds_mu;
  double d_cw1 = 0.0;
  d_cw1 += ((1.0 + cb2) * d_oosigma) * NavierStokesTerm::ooreynolds;
  d_cw1 += (cb1*oovkcst2 + (1.0 + cb2) * SATerm::oosigma) * dooreynolds_mu;
  //----

  double dmu5 = oosigma * (dmul + dabsmutilde) + d_oosigma * (mul + absmutilde);

  double dnutildedx = dp1dxj[0][0]*V[0][5] + dp1dxj[1][0]*V[1][5] +
    dp1dxj[2][0]*V[2][5] + dp1dxj[3][0]*V[3][5];
  double ddnutildedx = ddp1dxj[0][0]*V[0][5] + dp1dxj[0][0]*dV[0][5] + ddp1dxj[1][0]*V[1][5] + dp1dxj[1][0]*dV[1][5] +
    ddp1dxj[2][0]*V[2][5] + dp1dxj[2][0]*dV[2][5] + ddp1dxj[3][0]*V[3][5] + dp1dxj[3][0]*dV[3][5];
  double dnutildedy = dp1dxj[0][1]*V[0][5] + dp1dxj[1][1]*V[1][5] +
    dp1dxj[2][1]*V[2][5] + dp1dxj[3][1]*V[3][5];
  double ddnutildedy = ddp1dxj[0][1]*V[0][5] + dp1dxj[0][1]*dV[0][5] + ddp1dxj[1][1]*V[1][5] + dp1dxj[1][1]*dV[1][5] +
    ddp1dxj[2][1]*V[2][5] + dp1dxj[2][1]*dV[2][5] + ddp1dxj[3][1]*V[3][5] + dp1dxj[3][1]*dV[3][5];
  double dnutildedz = dp1dxj[0][2]*V[0][5] + dp1dxj[1][2]*V[1][5] +
    dp1dxj[2][2]*V[2][5] + dp1dxj[3][2]*V[3][5];
  double ddnutildedz = ddp1dxj[0][2]*V[0][5] + dp1dxj[0][2]*dV[0][5] + ddp1dxj[1][2]*V[1][5] + dp1dxj[1][2]*dV[1][5] +
    ddp1dxj[2][2]*V[2][5] + dp1dxj[2][2]*dV[2][5] + ddp1dxj[3][2]*V[3][5] + dp1dxj[3][2]*dV[3][5];

  dR[0][5] = dmu5 * dnutildedx + mu5 * ddnutildedx;
  dR[1][5] = dmu5 * dnutildedy + mu5 * ddnutildedy;
  dR[2][5] = dmu5 * dnutildedz + mu5 * ddnutildedz;

  dS[0] = 0.0;
  dS[1] = 0.0;
  dS[2] = 0.0;
  dS[3] = 0.0;
  dS[4] = 0.0;
  dS[5] = 0.0;

  double d2wall = 0.25 * (d2w[0] + d2w[1] + d2w[2] + d2w[3]);
  if (d2wall >= 1.e-15) {
    if (mutilde/mul == 0.001) {
      fprintf(stderr, "***** Inside the file FemEquationTermDesc.C the varibles in the function max are equal *****\n");
      //exit(1);
    }
    double chi = max(mutilde/mul, 0.001);
    double dchi;
    if (chi == 0.001)
      dchi = 0.0;
    else
      dchi = ( dmutilde * mul - mutilde * dmul ) / ( mul * mul );
    double chi3 = chi*chi*chi;
    double fv1 = chi3 / (chi3 + cv1_pow3);
    double dfv1 = ( 3.0*chi*chi*dchi*(chi3 + cv1_pow3) - chi3 * 3.0*chi*chi*dchi ) / ( (chi3 + cv1_pow3) * (chi3 + cv1_pow3) );

    double fv2  = 1.-chi/(1.+chi*fv1);
    double dfv2 = (fv2-1.)*dchi/chi+(1.-fv2)*(1-fv2)*(dfv1+fv1*dchi/chi);
    double fv3 = 1.0;
    double dfv3 = 0.0;
    if (usefv3) {
      fv2 = 1.0 + oocv2*chi;
      dfv2 = oocv2*dchi;
      dfv2 = -3.0 / (fv2*fv2*fv2*fv2)*dfv2;
      fv2 = 1.0 / (fv2*fv2*fv2);
      fv3 = (1.0 + chi*fv1) * (1.0 - fv2) / chi;
      dfv3 = ( ( dchi*fv1 + chi*dfv1 ) * (1.0 - fv2) * chi + (1.0 + chi*fv1) * (- dfv2) * chi - (1.0 + chi*fv1) * (1.0 - fv2) * dchi ) / ( chi * chi );
    }

    double ood2wall2 = 1.0 / (d2wall * d2wall);
    double rho = 0.25 * (V[0][0] + V[1][0] + V[2][0] + V[3][0]);
    double drho = 0.25 * (dV[0][0] + dV[1][0] + dV[2][0] + dV[3][0]);
    double oorho = 1.0 / rho;
    double doorho = -1.0 / ( rho * rho ) * drho;
    double zz = ooreynolds * oovkcst2 * mutilde * oorho * ood2wall2;
    double dzz = dooreynolds_mu * oovkcst2 * mutilde * oorho * ood2wall2 + ooreynolds_mu * oovkcst2 * dmutilde * oorho * ood2wall2 + ooreynolds * oovkcst2 * mutilde * doorho * ood2wall2;
    double s12 = dudxj[0][1] - dudxj[1][0];
    double ds12 = ddudxj[0][1] - ddudxj[1][0];
    double s23 = dudxj[1][2] - dudxj[2][1];
    double ds23 = ddudxj[1][2] - ddudxj[2][1];
    double s31 = dudxj[2][0] - dudxj[0][2];
    double ds31 = ddudxj[2][0] - ddudxj[0][2];
    double s = sqrt(s12*s12 + s23*s23 + s31*s31);
    double ds = 1.0 / ( 2.0*s ) * (2.0*s12*ds12 + 2.0*s23*ds23 + 2.0*s31*ds31);
    double Stilde = s*fv3 + zz*fv2;
    double dStilde = ds*fv3 + s*dfv3 + dzz*fv2 + zz*dfv2;
    double rr = min(zz/Stilde, 2.0);
    double drr;
    if (rr==2.0)
      drr = 0.0;
    else
      drr = ( dzz * Stilde - zz*dStilde ) / ( Stilde * Stilde );
    double rr2 = rr*rr;
    double gg = rr + cw2 * (rr2*rr2*rr2 - rr);
    double dgg = drr + cw2 * (6.0*rr*rr2*rr2*drr - drr);
    double gg2 = gg*gg;
    double fw = opcw3_pow * gg * pow(gg2*gg2*gg2 + cw3_pow6, -sixth);
    double dfw = opcw3_pow * dgg * pow(gg2*gg2*gg2 + cw3_pow6, -sixth) + opcw3_pow * gg * (-sixth) * pow(gg2*gg2*gg2 + cw3_pow6, (-sixth - 1.0) ) * 6.0*gg*gg2*gg2*dgg;

//  double AA = oosigma * cb2 * rho *
//    (dnutildedx*dnutildedx + dnutildedy*dnutildedy + dnutildedz*dnutildedz);
    double dAA = 0.0;
    dAA += d_oosigma * cb2 * rho * (dnutildedx*dnutildedx + dnutildedy*dnutildedy + dnutildedz*dnutildedz);
    dAA += oosigma * cb2 * drho * (dnutildedx*dnutildedx + dnutildedy*dnutildedy + dnutildedz*dnutildedz);
    dAA += oosigma * cb2 * rho * (2.0*dnutildedx*ddnutildedx + 2.0*dnutildedy*ddnutildedy + 2.0*dnutildedz*ddnutildedz);
//  double BB = cb1 * Stilde * absmutilde;
    double dBB = cb1 * dStilde * absmutilde + cb1 * Stilde * dabsmutilde;
//  double CC = - cw1 * fw * oorho * maxmutilde*maxmutilde * ood2wall2;
    double dCC = 0.0;
    dCC -= d_cw1 * fw * oorho * maxmutilde * maxmutilde * ood2wall2;
    dCC -= cw1 * dfw * oorho * maxmutilde * maxmutilde * ood2wall2;
    dCC -= cw1 * fw * doorho * maxmutilde * maxmutilde * ood2wall2;
    dCC -= cw1 * fw * oorho * 2.0 * maxmutilde * dmaxmutilde * ood2wall2;
    //-----
    dS[5] = dAA + dBB + dCC;
  }
  else {
    dS[5] = 0.0;
  }

  // Initialize PR (porous media term)
  for (int j=0; j<3*4; ++j) dPR[j] = 0.0; 

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
      porousmedia = true;
      fprintf(stderr, "***** Inside the file FemEquationTermDesc.C the derivative related to porus media is not implemented *****\n");
      exit(1);
    }
  }

  // If element wasn't flagged as porous media, treat it as standard fluid 
  if (!porousmedia) {
    mu = ooreynolds_mu * (mul + mut);
    dmu = dooreynolds_mu * (mul + mut) + ooreynolds_mu * (dmul + dmut);

    lambda  = ooreynolds_mu * (lambdal + lambdat);
    dlambda = 
        dooreynolds_mu * (lambdal + lambdat) + ooreynolds_mu * (dlambdal + dlambdat);

    double kappat = turbThermalCondFcn.turbulentConductivity(mut);
    kappa  = ooreynolds_mu  * (kappal + kappat);
    dkappa = 
        dooreynolds_mu * (kappal + kappat) +
        ooreynolds_mu  * (dkappal + turbThermalCondFcn.turbulentConductivityDerivative(dmut));

    computeDerivativeOfVolumeTermNS(mu, dmu, lambda, dlambda, kappa, dkappa, ucg, ducg, dudxj, ddudxj, dTdxj, ddTdxj, dR);
  }

  return (porousmedia);

}

//------------------------------------------------------------------------------

// This function was modified to account for the derivative of the mu with respect to the conservative variables
// Included (MB*)
bool FemEquationTermSA::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						  double *V[4], double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                  SVec<double,3> &X, int nodeNum[4], int material_id)
{

  bool porousmedia = false;

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double mul, lambdal, kappal;
  computeLaminarTransportCoefficients(Tcg, mul, lambdal, kappal);

  double mutilde;
  double mut, lambdat, kappat;
  computeTurbulentTransportCoefficients(V, nodeNum, X, mul, lambdal, kappal, mutilde, mut, lambdat, kappat);

  double mu     = ooreynolds_mu * (mul + mut);
  double lambda = ooreynolds_mu * (lambdal + lambdat);
  double kappa  = ooreynolds_mu * (kappal + kappat);

  int k;
  for (k=0; k<4*3*6*6; ++k)
    drdu[k] = 0.0;
  for (k=0; k<4*6*6; ++k)
    dsdu[k] = 0.0;
  for (k=0; k<4*4*6*6; ++k)
    dpdu[k] = 0.0;

  double (*dRdU)[3][6][6] = reinterpret_cast<double (*)[3][6][6]>(drdu);
  double (*dSdU)[6][6] = reinterpret_cast<double (*)[6][6]>(dsdu);
  double (*dPdU)[4][6][6] = reinterpret_cast<double (*)[4][6][6]>(dpdu);

  double dmul[4][6];
  double dmutilde[4][6];
  double dmut[4][6];
  double dmu[4][6];
  double dlambda[4][6];
  double dkappa[4][6];

  double dTcgdu0 = 0.0;
  double dTcgdu1 = 0.0;
  double dTcgdu2 = 0.0;
  double dTcgdu3 = 0.0;
  double dTcgdu4 = 0.0;
  double dTcgdu5 = 0.0;
  
  double dMach = 0.0;

  for (int k=0; k<4; ++k) {

    if (completeJac) {
      dTcgdu0 = - 0.25 / V[k][0] * (T[k] - 0.5 * (V[k][1]*V[k][1] + V[k][2]*V[k][2] + V[k][3]*V[k][3]));
      dTcgdu1 = - 0.25 / V[k][0] * V[k][1];
      dTcgdu2 = - 0.25 / V[k][0] * V[k][2];
      dTcgdu3 = - 0.25 / V[k][0] * V[k][3];
      dTcgdu4 = 0.25 / V[k][0];
      dTcgdu5 = 0.0;
    }

    dmul[k][0] = viscoFcn->compute_muDerivative(Tcg, dTcgdu0, dMach);
    dmul[k][1] = viscoFcn->compute_muDerivative(Tcg, dTcgdu1, dMach);
    dmul[k][2] = viscoFcn->compute_muDerivative(Tcg, dTcgdu2, dMach);
    dmul[k][3] = viscoFcn->compute_muDerivative(Tcg, dTcgdu3, dMach); 
    dmul[k][4] = viscoFcn->compute_muDerivative(Tcg, dTcgdu4, dMach); 
    dmul[k][5] = viscoFcn->compute_muDerivative(Tcg, dTcgdu5, dMach); 

    dmutilde[k][0] = 0.0;
    dmutilde[k][1] = 0.0;
    dmutilde[k][2] = 0.0;
    dmutilde[k][3] = 0.0;
    dmutilde[k][4] = 0.0;
    if (completeJac) 
      dmutilde[k][5] = 0.25;
    else
      dmutilde[k][5] = 0.0;

    for(int i=0; i<5; ++i)
      dmut[k][i] = computeDerivativeOfTurbulentViscosity(V, mul, dmul[k][i], dmutilde[k][i]);

    dkappa[k][0] = ooreynolds_mu * (thermalCondFcn->computeDerivative(Tcg, dTcgdu0, dMach) + turbThermalCondFcn.turbulentConductivityDerivative(dmut[k][0]));
    dkappa[k][1] = ooreynolds_mu * (thermalCondFcn->computeDerivative(Tcg, dTcgdu1, dMach) + turbThermalCondFcn.turbulentConductivityDerivative(dmut[k][1]));
    dkappa[k][2] = ooreynolds_mu * (thermalCondFcn->computeDerivative(Tcg, dTcgdu2, dMach) + turbThermalCondFcn.turbulentConductivityDerivative(dmut[k][2]));
    dkappa[k][3] = ooreynolds_mu * (thermalCondFcn->computeDerivative(Tcg, dTcgdu3, dMach) + turbThermalCondFcn.turbulentConductivityDerivative(dmut[k][3]));
    dkappa[k][4] = ooreynolds_mu * (thermalCondFcn->computeDerivative(Tcg, dTcgdu4, dMach) + turbThermalCondFcn.turbulentConductivityDerivative(dmut[k][4]));
    dkappa[k][5] = ooreynolds_mu * (thermalCondFcn->computeDerivative(Tcg, dTcgdu5, dMach) + turbThermalCondFcn.turbulentConductivityDerivative(dmut[k][5]));

    for(int i=0; i<5; ++i){
      dmu[k][i] = ooreynolds_mu * (dmul[k][i] + dmut[k][i]);
      double dlambdal = viscoFcn->compute_lambdaDerivative(mu, dmu[k][i], dMach);
      double dlambdat = computeDerivativeOfSecondTurbulentViscosity(lambdal, dlambdal, mul, dmul[k][i], mut, dmut[k][i]);
      dlambda[k][i] = ooreynolds_mu * (dlambdal + dlambdat);
    }

  }
  
  double s = sqrt((dudxj[0][1] - dudxj[1][0])*(dudxj[0][1] - dudxj[1][0]) + (dudxj[1][2] - dudxj[2][1])*(dudxj[1][2] - dudxj[2][1]) + (dudxj[2][0] - dudxj[0][2])*(dudxj[2][0] - dudxj[0][2]));
  if (s != 0.0)
    computeJacobianVolumeTermSA<6,5>(dp1dxj, d2w, dudxj, mul, dmul, mutilde, dmutilde, V, dRdU, dSdU);
  else
    computeJacobianVolumeTermSA<6,5>(dp1dxj, d2w, dudxj, mul, mutilde, V, dRdU, dSdU);
     
  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {
      // if porous media with material_id has been defined in the input file
      mut     = computePorousTurbulentViscosity(it, ucg, length);
      lambdat = computeSecondPorousTurbulentViscosity(lambdal, mul, mut);
      kappat  = turbThermalCondFcn.turbulentConductivity(mut);
      mu      = ooreynolds_mu * (mul + mut);
      lambda  = ooreynolds_mu * (lambdal + lambdat);
      kappa   = ooreynolds_mu * (kappal + kappat);

      // approximate jacobian for the viscous term
      // (approximate = transport coefficients are frozen)
      computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);

      // jacobian for the porous term
      porousmedia = computeJacobianVolumeTermPorousCore(tetVol, it, length, density, velocity, ucg, V, dPdU);

    }
  }

  // If element wasn't flagged as porous media, treat it as standard fluid 
  // jacobian takes into account derivatives of mu, lambda and kappa
  if (!porousmedia)
    computeJacobianVolumeTermNS(dp1dxj, mu, dmu, 
        lambda, dlambda, kappa, dkappa, V, T, dRdU);

  return (porousmedia);


}

//------------------------------------------------------------------------------

void FemEquationTermSA::computeSurfaceTerm(int code, Vec3D &n, double d2w[3], 
					   double *Vwall, double *V[3], double *R)
{

  wallFcn->computeSurfaceTerm(code, n, d2w, Vwall, V, R);

  R[5] = 0.0;

}

//------------------------------------------------------------------------------

double FemEquationTermSA::computeNormDerivWallFcn(double rho, double T, double Du1, 
																double DT1, double d2w, 
																double &dudn, double &dTdn)
{

	double ut = wallFcn->computedudT(rho, T, Du1, DT1, d2w, dudn, dTdn);

	return ut;

}
//------------------------------------------------------------------------------


// Included (MB)
void FemEquationTermSA::computeDerivativeOfSurfaceTerm(int code, Vec3D &n, Vec3D &dn, double d2w[3],
					   double *Vwall, double *dVwall, double *V[3], double *dV[3], double dMach, double *dR)
{

  wallFcn->computeDerivativeOfSurfaceTerm(code, n, dn, d2w, Vwall, dVwall, V, dV, dMach, dR);

  dR[5] = 0.0;

}

//------------------------------------------------------------------------------

void FemEquationTermSA::computeJacobianSurfaceTerm(int code, Vec3D &n, 
						   double d2w[3], double *Vwall, 
						   double *V[3], double *drdu)
{

  for (int k=0; k<3*6*6; ++k)
    drdu[k] = 0.0;

  double (*dRdU)[6][6] = reinterpret_cast<double (*)[6][6]>(drdu);
  wallFcn->computeJacobianSurfaceTerm(code, n, d2w, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

// Included (MB)
void FemEquationTermSA::computeBCsJacobianWallValues(int code, Vec3D &n, double d2w[3], double *Vwall, double *dVwall, double *V[3])
{

  wallFcn->computeBCsJacobianWallValues<6>(code, n, d2w, Vwall, dVwall, V);

}

//------------------------------------------------------------------------------

void FemEquationTermSA::computeSurfaceTerm(double dp1dxj[4][3], int code,
					   Vec3D &n, double d2w[4],
					   double *Vwall, double *V[4], double *R)
{
  
  computeSurfaceTermNS(dp1dxj, n, Vwall, V, R);

  R[5] = 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
void FemEquationTermSA::computeDerivativeOfSurfaceTerm(double dp1dxj[4][3], double ddp1dxj[4][3], int code,
					   Vec3D &n, Vec3D &dn, double d2w[4],
					   double *Vwall, double *dVwall, double *V[4], double *dV[4], double dMach, double *dR)
{

  computeDerivativeOfSurfaceTermNS(dp1dxj, ddp1dxj, n, dn, Vwall, dVwall, V, dV, dMach, dR);

  dR[5] = 0.0;

}

//------------------------------------------------------------------------------

void FemEquationTermSA::computeJacobianSurfaceTerm(double dp1dxj[4][3], int code,
						   Vec3D &n, double d2w[4], 
						   double *Vwall, double *V[4], 
						   double *drdu)
{

  for (int k=0; k<4*6*6; ++k)
    drdu[k] = 0.0;

  double (*dRdU)[6][6] = reinterpret_cast<double (*)[6][6]>(drdu);
  computeJacobianSurfaceTermNS(dp1dxj, n, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

FemEquationTermSAmean::FemEquationTermSAmean(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), SATerm(iod), FemEquationTerm(iod.volumes.volumeMap.dataMap),
  turbThermalCondFcn(iod, viscoFcn, vf)

{

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = new WallFcn(iod, varFcn, viscoFcn);

  x0 =  iod.eqs.tc.tr.bfix.x0;
  x1 =  iod.eqs.tc.tr.bfix.x1;
  y0 =  iod.eqs.tc.tr.bfix.y0;
  y1 =  iod.eqs.tc.tr.bfix.y1;
  z0 =  iod.eqs.tc.tr.bfix.z0;
  z1 =  iod.eqs.tc.tr.bfix.z1;

  if(x0>x1 || y0>y1 || z0>z1)  trip = 0;
  else    trip = 1;

  velocity = iod.ref.rv.velocity;
  density = iod.ref.rv.density;
  length = iod.ref.rv.length;

}

//------------------------------------------------------------------------------

void FemEquationTermSAmean::computeLaminarTransportCoefficients(
  const double T, double &mul, double &lambdal, double &kappal)
{

  mul     = viscoFcn->compute_mu(T);
  lambdal = viscoFcn->compute_lambda(T,mul);
  kappal  = thermalCondFcn->compute(T);

}

//------------------------------------------------------------------------------

void FemEquationTermSAmean::computeTurbulentTransportCoefficients(
  double *V[], int nodeNum[], SVec<double,3> &X,
  const double mul, const double lambdal, const double kappal,
  double &mutilde, double &mut, double &lambdat, double &kappat)
{

  mut = computeTurbulentViscosity(V, mul, mutilde);

  //Applying the laminar-turbulent trip
  if(trip) {
    int in_trip = 0;
    for (int k=0; k<3; k++)
      if (X[nodeNum[k]][0]>=x0 && X[nodeNum[k]][0]<=x1 &&
          X[nodeNum[k]][1]>=y0 && X[nodeNum[k]][1]<=y1 &&
	  X[nodeNum[k]][2]>=z0 && X[nodeNum[k]][2]<=z1 )
        in_trip = 1;

    if (!in_trip) mut = 0.0;
  }
  lambdat = computeSecondTurbulentViscosity(lambdal, mul, mut);
  kappat  = turbThermalCondFcn.turbulentConductivity(mut);

}

//------------------------------------------------------------------------------

double FemEquationTermSAmean::computeViscousTimeStep(double X[3], double *V)
{

  double T;
  computeTemperature(V,T);
  double mul = viscoFcn->compute_mu(T);
  double mut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1)
       mut = computeTurbulentViscosity(V, mul);
    else
       mut = 0.0;
  }
  else
    mut = computeTurbulentViscosity(V, mul);

  return ooreynolds_mu * (mul+mut)/V[0];

}

//------------------------------------------------------------------------------

// Included (MB)
double FemEquationTermSAmean::computeDerivativeOfViscousTimeStep(double X[3], double dX[3], double *V, double *dV, double dMach)
{

  double T,dT;
  computeTemperature(V,T);
  computeDerivativeOfTemperature(V,dV,dT);
  double mul = viscoFcn->compute_mu(T);
  double dmul = viscoFcn->compute_muDerivative(T,dT,dMach);
  double mut, dmut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1) {
       mut = computeTurbulentViscosity(V, mul);
       dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
    }
    else {
       mut = 0.0;
       dmut = 0.0;
    }
  }
  else {
    mut = computeTurbulentViscosity(V, mul);
    dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
  }

  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;

  return dooreynolds_mu * (mul+mut)/V[0] + ooreynolds_mu * ((dmul+dmut)*V[0]-(mul+mut)*dV[0])/(V[0]*V[0]);

}

//------------------------------------------------------------------------------

bool FemEquationTermSAmean::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						      double *V[4], double *drdu, double *dSdU, double *dpdu, double tetVol,
                                                      SVec<double,3> &X, int nodeNum[4], int material_id)
{

  bool porousmedia = false;
                                                                                                                                                                   
  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);
                                                                                                                                                                   
  double T[4], Tcg;
  computeTemperature(V, T, Tcg);
                                                                                                                                                                   
  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);
                                                                                                                                                                   
  double mul, lambdal, kappal;
  computeLaminarTransportCoefficients(Tcg, mul, lambdal, kappal);

  double mutilde;
  double mut, lambdat, kappat; //may be modified later if element is flagged as a porous media
  computeTurbulentTransportCoefficients(V, nodeNum, X, mul, lambdal, kappal, mutilde, mut, lambdat, kappat);

  double (*dRdU)[3][5][5] = reinterpret_cast<double (*)[3][5][5]>(drdu);
  double (*dPdU)[4][5][5] = reinterpret_cast<double (*)[4][5][5]>(dpdu);

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {
      // if porous media with material_id has been defined in the input file
      mut     = computePorousTurbulentViscosity(it, ucg, length);
      lambdat = computeSecondPorousTurbulentViscosity(lambdal, mul, mut);
      kappat  = turbThermalCondFcn.turbulentConductivity(mut);

      // jacobian for the porous term
      porousmedia = computeJacobianVolumeTermPorousCore(tetVol, it, length, density, velocity, ucg, V, dPdU);

    }
  }

  // In all cases, porous media or not, approximate jacobian of viscous term is computed
  // (approximate = transport coefficients are frozen)
  double mu, lambda, kappa;
  mu = ooreynolds_mu * (mul + mut);
  lambda  = ooreynolds_mu * (lambdal + lambdat);
  kappa   = ooreynolds_mu * (kappal + kappat);
  computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);

  return (porousmedia);

}

//------------------------------------------------------------------------------

void FemEquationTermSAmean::computeJacobianSurfaceTerm(int code, Vec3D &n, 
						       double d2w[3], double *Vwall, 
						       double *V[3], double *drdu)
{

  double (*dRdU)[5][5] = reinterpret_cast<double (*)[5][5]>(drdu);
  wallFcn->computeJacobianSurfaceTerm(code, n, d2w, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

void FemEquationTermSAmean::computeJacobianSurfaceTerm(double dp1dxj[4][3], int code,
						       Vec3D &n, double d2w[4],
						       double *Vwall, double *V[4],
						       double *drdu)
{

  double (*dRdU)[5][5] = reinterpret_cast<double (*)[5][5]>(drdu);
  computeJacobianSurfaceTermNS(dp1dxj, n, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

FemEquationTermSAturb::FemEquationTermSAturb(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), SATerm(iod), FemEquationTerm(iod.volumes.volumeMap.dataMap)
{

}

//------------------------------------------------------------------------------

bool FemEquationTermSAturb::computeJacobianVolumeTerm(double dp1dxj[4][3], 
						      double d2w[4], double *V[4], 
						      double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                      SVec<double,3> &X, int nodeNum[4], int material_id)
{

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double mul = viscoFcn->compute_mu(Tcg);
  double mutilde;

  computeTurbulentViscosity(V, mul, mutilde);

  double (*dRdU)[3][1][1] = reinterpret_cast<double (*)[3][1][1]>(drdu);
  double (*dSdU)[1][1] = reinterpret_cast<double (*)[1][1]>(dsdu);
  computeJacobianVolumeTermSA<1,0>(dp1dxj, d2w, dudxj, mul, mutilde, V, dRdU, dSdU);

  return false;

}

//------------------------------------------------------------------------------

FemEquationTermDES::FemEquationTermDES(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), DESTerm(iod), FemEquationTerm(iod.volumes.volumeMap.dataMap),
  turbThermalCondFcn(iod, viscoFcn, vf)
{

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = new WallFcnSA(iod, varFcn, viscoFcn);

  cdes = iod.eqs.tc.tm.des.cdes;

  x0 =  iod.eqs.tc.tr.bfix.x0;
  x1 =  iod.eqs.tc.tr.bfix.x1;
  y0 =  iod.eqs.tc.tr.bfix.y0;
  y1 =  iod.eqs.tc.tr.bfix.y1;
  z0 =  iod.eqs.tc.tr.bfix.z0;
  z1 =  iod.eqs.tc.tr.bfix.z1;

  if (x0>x1 || y0>y1 || z0>z1) trip = 0;
  else   trip = 1;

  if (iod.ts.implicit.tmcoupling == ImplicitData::STRONG && trip == 1) { 
    fprintf(stderr,"** Warning: Laminar-turbulent trip not implemented for Strongly Coupled NS-DES simulation \n");
    trip = 0;
  }

  if (iod.eqs.tc.tm.des.form == DESModelData::FV3)
    usefv3 = true;
  else
    usefv3 = false;

  velocity = iod.ref.rv.velocity;
  density = iod.ref.rv.density;
  length = iod.ref.rv.length;

// Included (MB)
   if (iod.eqs.fluidModel.fluid == FluidModelData::PERFECT_GAS)
     completeJac = true;
   else
     completeJac = false;

}

//------------------------------------------------------------------------------

void FemEquationTermDES::computeLaminarTransportCoefficients(
  const double T, double &mul, double &lambdal, double &kappal)
{

  mul     = viscoFcn->compute_mu(T);
  lambdal = viscoFcn->compute_lambda(T,mul);
  kappal  = thermalCondFcn->compute(T);

}

//------------------------------------------------------------------------------

void FemEquationTermDES::computeTurbulentTransportCoefficients(
  double *V[], int nodeNum[], SVec<double,3> &X,
  const double mul, const double lambdal, const double kappal,
  double &mutilde, double &mut, double &lambdat, double &kappat)
{

  mut = computeTurbulentViscosity(V, mul, mutilde);

  //Applying the laminar-turbulent trip
  if(trip) {
    int in_trip = 0;
    for (int k=0; k<3; k++)
      if (X[nodeNum[k]][0]>=x0 && X[nodeNum[k]][0]<=x1 &&
          X[nodeNum[k]][1]>=y0 && X[nodeNum[k]][1]<=y1 &&
	  X[nodeNum[k]][2]>=z0 && X[nodeNum[k]][2]<=z1 )
        in_trip = 1;

    if (!in_trip) mut = 0.0;
  }
  lambdat = computeSecondTurbulentViscosity(lambdal, mul, mut);
  kappat  = turbThermalCondFcn.turbulentConductivity(mut);

}

//------------------------------------------------------------------------------

double FemEquationTermDES::computeViscousTimeStep(double X[3], double *V)
{

  double T;
  computeTemperature(V,T);
  double mul = viscoFcn->compute_mu(T);
  double mut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1)
       mut = computeTurbulentViscosity(V, mul);
    else
       mut = 0.0;
  }
  else
    mut = computeTurbulentViscosity(V, mul);

  return ooreynolds_mu * (mul+mut)/V[0];

}

//------------------------------------------------------------------------------

// Included (MB)
double FemEquationTermDES::computeDerivativeOfViscousTimeStep(double X[3], double dX[3], double *V, double *dV, double dMach)
{

  double T, dT;
  computeTemperature(V,T);
  computeDerivativeOfTemperature(V,dV,dT);
  double mul = viscoFcn->compute_mu(T);
  double dmul = viscoFcn->compute_muDerivative(T,dT,dMach);
  double mut, dmut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1) {
       mut = computeTurbulentViscosity(V, mul);
       dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
    }
    else {
       mut = 0.0;
       dmut = 0.0;
    }
  }
  else {
    mut = computeTurbulentViscosity(V, mul);
    dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
  }

  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;

  return dooreynolds_mu * (mul+mut)/V[0] + ooreynolds_mu * ((dmul+dmut)*V[0]-(mul+mut)*dV[0])/(V[0]*V[0]);

}

//------------------------------------------------------------------------------

bool FemEquationTermDES::computeVolumeTerm(double dp1dxj[4][3], double d2w[4], 
					  double *V[4], double *r, double *S, double *PR, double tetVol,
                                          SVec<double,3> &X,
                                          int nodeNum[4], int material_id)
{

  bool porousmedia = false;

  const double sixth = 1.0/6.0;

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double dTdxj[3];
  computeTemperatureGradient(dp1dxj, T, dTdxj);
 
  double mul, lambdal, kappal;
  computeLaminarTransportCoefficients(Tcg, mul, lambdal, kappal);

  double mutilde;
  double mut, lambdat, kappat; //may be modified later if element is flagged as a porous media
  computeTurbulentTransportCoefficients(V, nodeNum, X, mul, lambdal, kappal, mutilde, mut, lambdat, kappat);

  double (*R)[6] = reinterpret_cast<double (*)[6]>(r);

  double absmutilde = fabs(mutilde);
  double maxmutilde = max(mutilde, 0.0);
  double mu5 = oosigma * (mul + absmutilde);
  double dnutildedx = dp1dxj[0][0]*V[0][5] + dp1dxj[1][0]*V[1][5] + 
    dp1dxj[2][0]*V[2][5] + dp1dxj[3][0]*V[3][5];
  double dnutildedy = dp1dxj[0][1]*V[0][5] + dp1dxj[1][1]*V[1][5] + 
    dp1dxj[2][1]*V[2][5] + dp1dxj[3][1]*V[3][5];
  double dnutildedz = dp1dxj[0][2]*V[0][5] + dp1dxj[1][2]*V[1][5] + 
    dp1dxj[2][2]*V[2][5] + dp1dxj[3][2]*V[3][5];

  R[0][5] = mu5 * dnutildedx;
  R[1][5] = mu5 * dnutildedy;
  R[2][5] = mu5 * dnutildedz;

  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;

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

  if  (d2wall >= 1.e-15) {
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
    double zz = ooreynolds_mu * oovkcst2 * maxmutilde * oorho * ood2wall2;
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

    double AA = oosigma * cb2 * rho * 
      (dnutildedx*dnutildedx + dnutildedy*dnutildedy + dnutildedz*dnutildedz);
    double BB = cb1 * Stilde * absmutilde;
    double CC = - cw1 * fw * oorho * maxmutilde*maxmutilde * ood2wall2;
    S[5] = AA + BB + CC;
  }
  else {
    S[5] = 0.0;
  }

  // Initialize PR (porous media term)
  for (int j=0; j<3*4; ++j) PR[j] = 0.0; 

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
      porousmedia = true;

      mut     = computePorousTurbulentViscosity(it, ucg, length);
      lambdat = computeSecondPorousTurbulentViscosity(lambdal, mul, mut);
      kappat  = turbThermalCondFcn.turbulentConductivity(mut);

      porousmedia = computeVolumeTermPorousCore(tetVol, it, length, density, velocity, ucg, V, PR);
    }
  }

  // In all cases, porous media or not, compute viscous term for Navier-Stokes equations
  double mu, lambda, kappa;
  mu     = ooreynolds_mu * (mul + mut);
  lambda = ooreynolds_mu * (lambdal + lambdat);
  kappa  = ooreynolds_mu * (kappal + kappat);
  computeVolumeTermNS(mu, lambda, kappa, ucg, dudxj, dTdxj, R);

  return (porousmedia);

}

//------------------------------------------------------------------------------

// Included (MB)
bool FemEquationTermDES::computeDerivativeOfVolumeTerm(double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4],
					  double *V[4], double *dV[4], double dMach, double *dr, double *dS, double *dPR, double dtetVol, SVec<double,3> &X,
                                          int nodeNum[4], int material_id)
{

  bool porousmedia = false;
  sleep(1);//TODO delete line

  fprintf(stderr, "***** FemEquationTermDesc::computeDerivativeOfVolumeTerm");//TODO uncomment line
  fprintf(stderr, " >> Function not defined for TermDES *****\n");//TODO uncomment line

  sleep(10);//TODO delete line
  fflush(stdout);//TODO delete line
  exit(-1);

  return (porousmedia);

}

//------------------------------------------------------------------------------

bool FemEquationTermDES::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						  double *V[4], double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                  SVec<double,3> &X, int nodeNum[4], int material_id)
{

  bool porousmedia = false;

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);
                                                                                                   
  double mul, lambdal, kappal;
  computeLaminarTransportCoefficients(Tcg, mul, lambdal, kappal);

  double mutilde;
  double mut, lambdat, kappat; //may be modified later if element is flagged as a porous media
  computeTurbulentTransportCoefficients(V, nodeNum, X, mul, lambdal, kappal, mutilde, mut, lambdat, kappat);

  int k;
  for (k=0; k<4*3*6*6; ++k)
    drdu[k] = 0.0;
  for (k=0; k<4*6*6; ++k)
    dsdu[k] = 0.0;
  for (k=0; k<4*4*6*6; ++k)
    dpdu[k] = 0.0;

  double (*dRdU)[3][6][6] = reinterpret_cast<double (*)[3][6][6]>(drdu);
  double (*dSdU)[6][6] = reinterpret_cast<double (*)[6][6]>(dsdu);
  double (*dPdU)[4][6][6] = reinterpret_cast<double (*)[4][6][6]>(dpdu);

  computeJacobianVolumeTermDES<6,5>(dp1dxj, d2w, dudxj, mul, mutilde, V, dRdU, dSdU, X, nodeNum);

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {
      // if porous media with material_id has been defined in the input file

      mut     = computePorousTurbulentViscosity(it, ucg, length);
      lambdat = computeSecondPorousTurbulentViscosity(lambdal, mul, mut);
      kappat  = turbThermalCondFcn.turbulentConductivity(mut);

      // jacobian for the porous term
      porousmedia = computeJacobianVolumeTermPorousCore(tetVol, it, length, density, velocity, ucg, V, dPdU);

    }
  }

  // In all cases, porous media or not, approximate jacobian of viscous term is computed
  // (approximate = transport coefficients are frozen)
  double mu, lambda, kappa;
  mu      = ooreynolds_mu * (mul + mut);
  lambda  = ooreynolds_mu * (lambdal + lambdat);
  kappa   = ooreynolds_mu * (kappal + kappat);
  computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);

  return (porousmedia);

}

//------------------------------------------------------------------------------

void FemEquationTermDES::computeSurfaceTerm(int code, Vec3D &n, double d2w[3], 
					   double *Vwall, double *V[3], double *R)
{

  wallFcn->computeSurfaceTerm(code, n, d2w, Vwall, V, R);

  R[5] = 0.0;

}

//------------------------------------------------------------------------------

double FemEquationTermDES::computeNormDerivWallFcn(double rho, double T, double Du1, 
																	double DT1, double d2w, 
																	double &dudn, double &dTdn)
{
	//fprintf(stderr, "FemEquationTermDES::computeNormDerivWallFcn not implemented\n");
	//return 0.0;

 	double ut = wallFcn->computedudT(rho, T, Du1, DT1, d2w, dudn, dTdn);
	return ut;
}

//------------------------------------------------------------------------------

// Included (MB)
void FemEquationTermDES::computeDerivativeOfSurfaceTerm(int code, Vec3D &n, Vec3D &dn, double d2w[3],
					   double *Vwall, double *dVwall, double *V[3], double *dV[3], double dMach, double *dR)
{
  fprintf(stderr, "***** Inside the file FemEquationTermDesc.C the derivative related to FemEquationTermDES::computeDerivativeOfSurfaceTerm is not implemented *****\n");
  exit(1);
}

//------------------------------------------------------------------------------

void FemEquationTermDES::computeJacobianSurfaceTerm(int code, Vec3D &n, 
						   double d2w[3], double *Vwall, 
						   double *V[3], double *drdu)
{

  for (int k=0; k<3*6*6; ++k)
    drdu[k] = 0.0;

  double (*dRdU)[6][6] = reinterpret_cast<double (*)[6][6]>(drdu);
  wallFcn->computeJacobianSurfaceTerm(code, n, d2w, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

void FemEquationTermDES::computeSurfaceTerm(double dp1dxj[4][3], int code,
					   Vec3D &n, double d2w[4],
					   double *Vwall, double *V[4], double *R)
{
  
  computeSurfaceTermNS(dp1dxj, n, Vwall, V, R);

  R[5] = 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
void FemEquationTermDES::computeDerivativeOfSurfaceTerm(double dp1dxj[4][3], double ddp1dxj[4][3], int code, Vec3D &n, Vec3D &dn, double d2w[4], double *Vwall, double *dVwall, double *V[4], double *dV[4], double dMach, double *dR)
{
  fprintf(stderr, "***** Inside the file FemEquationTermDesc.C the derivative related to FemEquationTermDES::computeDerivativeOfSurfaceTerm is not implemented *****\n");
  exit(1);
}

//------------------------------------------------------------------------------

void FemEquationTermDES::computeJacobianSurfaceTerm(double dp1dxj[4][3], int code,
						   Vec3D &n, double d2w[4], 
						   double *Vwall, double *V[4], 
						   double *drdu)
{

  for (int k=0; k<4*6*6; ++k)
    drdu[k] = 0.0;

  double (*dRdU)[6][6] = reinterpret_cast<double (*)[6][6]>(drdu);
  computeJacobianSurfaceTermNS(dp1dxj, n, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

FemEquationTermDESmean::FemEquationTermDESmean(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), DESTerm(iod), FemEquationTerm(iod.volumes.volumeMap.dataMap),
  turbThermalCondFcn(iod, viscoFcn, vf)
{

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = new WallFcn(iod, varFcn, viscoFcn);

  x0 =  iod.eqs.tc.tr.bfix.x0;
  x1 =  iod.eqs.tc.tr.bfix.x1;
  y0 =  iod.eqs.tc.tr.bfix.y0;
  y1 =  iod.eqs.tc.tr.bfix.y1;
  z0 =  iod.eqs.tc.tr.bfix.z0;
  z1 =  iod.eqs.tc.tr.bfix.z1;

  if(x0>x1 || y0>y1 || z0>z1)  trip = 0;
  else    trip = 1;

  velocity = iod.ref.rv.velocity;
  density = iod.ref.rv.density;
  length = iod.ref.rv.length;

}

//------------------------------------------------------------------------------

void FemEquationTermDESmean::computeLaminarTransportCoefficients(
  const double T, double &mul, double &lambdal, double &kappal)
{

  mul     = viscoFcn->compute_mu(T);
  lambdal = viscoFcn->compute_lambda(T,mul);
  kappal  = thermalCondFcn->compute(T);

}

//------------------------------------------------------------------------------

void FemEquationTermDESmean::computeTurbulentTransportCoefficients(
  double *V[], int nodeNum[], SVec<double,3> &X,
  const double mul, const double lambdal, const double kappal,
  double &mutilde, double &mut, double &lambdat, double &kappat)
{

  mut = computeTurbulentViscosity(V, mul, mutilde);

  //Applying the laminar-turbulent trip
  if(trip) {
    int in_trip = 0;
    for (int k=0; k<3; k++)
      if (X[nodeNum[k]][0]>=x0 && X[nodeNum[k]][0]<=x1 &&
          X[nodeNum[k]][1]>=y0 && X[nodeNum[k]][1]<=y1 &&
	  X[nodeNum[k]][2]>=z0 && X[nodeNum[k]][2]<=z1 )
        in_trip = 1;

    if (!in_trip) mut = 0.0;
  }
  lambdat = computeSecondTurbulentViscosity(lambdal, mul, mut);
  kappat  = turbThermalCondFcn.turbulentConductivity(mut);

}

//------------------------------------------------------------------------------

double FemEquationTermDESmean::computeViscousTimeStep(double X[3], double *V)
{

  double T;
  computeTemperature(V,T);
  double mul = viscoFcn->compute_mu(T);
  double mut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1)
       mut = computeTurbulentViscosity(V, mul);
    else
       mut = 0.0;
  }
  else
    mut = computeTurbulentViscosity(V, mul);

  return ooreynolds_mu * (mul+mut)/V[0];

}

//------------------------------------------------------------------------------

// Included (MB)
double FemEquationTermDESmean::computeDerivativeOfViscousTimeStep(double X[3], double dX[3], double *V, double *dV, double dMach)
{

  double T,dT;
  computeTemperature(V,T);
  computeDerivativeOfTemperature(V,dV,dT);
  double mul = viscoFcn->compute_mu(T);
  double dmul = viscoFcn->compute_muDerivative(T,dT,dMach);
  double mut, dmut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1) {
       mut = computeTurbulentViscosity(V, mul);
       dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
    }
    else {
       mut = 0.0;
       dmut = 0.0;
    }
  }
  else {
    mut = computeTurbulentViscosity(V, mul);
    dmut = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
  }

  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;

  return dooreynolds_mu * (mul+mut)/V[0] + ooreynolds_mu * ((dmul+dmut)*V[0]-(mul+mut)*dV[0])/(V[0]*V[0]);

}

//------------------------------------------------------------------------------

bool FemEquationTermDESmean::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						      double *V[4], double *drdu, double *dSdU, double *dpdu, double tetVol,
                                                      SVec<double,3> &X, int nodeNum[4], int material_id)
{

  bool porousmedia = false;
                                                                                                                                                                   
  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);
                                                                                                                                                                   
  double T[4], Tcg;
  computeTemperature(V, T, Tcg);
                                                                                                                                                                   
  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double mul, lambdal, kappal;
  computeLaminarTransportCoefficients(Tcg, mul, lambdal, kappal);

  double mutilde;
  double mut, lambdat, kappat; //may be modified later if element is flagged as a porous media
  computeTurbulentTransportCoefficients(V, nodeNum, X, mul, lambdal, kappal, mutilde, mut, lambdat, kappat);

  double (*dRdU)[3][5][5] = reinterpret_cast<double (*)[3][5][5]>(drdu);
  double (*dPdU)[4][5][5] = reinterpret_cast<double (*)[4][5][5]>(dpdu);

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {
      // if porous media with material_id has been defined in the input file

      mut     = computePorousTurbulentViscosity(it, ucg, length);
      lambdat = computeSecondPorousTurbulentViscosity(lambdal, mul, mut);
      kappat  = turbThermalCondFcn.turbulentConductivity(mut);

      // jacobian for the porous term
      porousmedia = computeJacobianVolumeTermPorousCore(tetVol, it, length, density, velocity, ucg, V, dPdU);

    }
  }

  // In all cases, porous media or not, approximate jacobian of viscous term is computed
  // (approximate = transport coefficients are frozen)
  double mu, lambda, kappa;
  mu     = ooreynolds_mu * (mul + mut);
  lambda = ooreynolds_mu * (lambdal + lambdat);
  kappa  = ooreynolds_mu * (kappal + kappat);
  computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);

  return (porousmedia);

}

//------------------------------------------------------------------------------

void FemEquationTermDESmean::computeJacobianSurfaceTerm(int code, Vec3D &n, 
						       double d2w[3], double *Vwall, 
						       double *V[3], double *drdu)
{

  double (*dRdU)[5][5] = reinterpret_cast<double (*)[5][5]>(drdu);
  wallFcn->computeJacobianSurfaceTerm(code, n, d2w, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

void FemEquationTermDESmean::computeJacobianSurfaceTerm(double dp1dxj[4][3], int code,
						       Vec3D &n, double d2w[4],
						       double *Vwall, double *V[4],
						       double *drdu)
{

  double (*dRdU)[5][5] = reinterpret_cast<double (*)[5][5]>(drdu);
  computeJacobianSurfaceTermNS(dp1dxj, n, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

FemEquationTermDESturb::FemEquationTermDESturb(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), DESTerm(iod), FemEquationTerm(iod.volumes.volumeMap.dataMap)
{

}

//------------------------------------------------------------------------------

bool FemEquationTermDESturb::computeJacobianVolumeTerm(double dp1dxj[4][3], 
						      double d2w[4], double *V[4], 
						      double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                      SVec<double,3> &X, int nodeNum[4], int material_id)
{

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);
                                                                                                  
  double mul = viscoFcn->compute_mu(Tcg);
  double mutilde;
  computeTurbulentViscosity(V, mul, mutilde);

  double (*dRdU)[3][1][1] = reinterpret_cast<double (*)[3][1][1]>(drdu);
  double (*dSdU)[1][1] = reinterpret_cast<double (*)[1][1]>(dsdu);
  computeJacobianVolumeTermDES<1,0>(dp1dxj, d2w, dudxj, mul, mutilde, V, dRdU, dSdU, X, nodeNum);

  return false;

}

//------------------------------------------------------------------------------

FemEquationTermKE::FemEquationTermKE(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), KEpsilonTerm(iod), FemEquationTerm(iod.volumes.volumeMap.dataMap),
  turbThermalCondFcn(iod, viscoFcn, vf)
{

  wallFcn = new WallFcnKE(iod, varFcn, viscoFcn);

  x0 = iod.eqs.tc.tr.bfix.x0;
  x1 =  iod.eqs.tc.tr.bfix.x1;
  y0 =  iod.eqs.tc.tr.bfix.y0;
  y1 =  iod.eqs.tc.tr.bfix.y1;
  z0 =  iod.eqs.tc.tr.bfix.z0;
  z1 =  iod.eqs.tc.tr.bfix.z1;

  if (x0>x1 || y0>y1 || z0>z1) trip = 0;
  else   trip = 1;

  if (iod.ts.implicit.tmcoupling == ImplicitData::STRONG && trip == 1) { 
    fprintf(stderr,"** Warning: Laminar-turbulent trip not implemented for Strongly Coupled NS-KEpsilon simulation \n");
    trip = 0;
  }

  velocity = iod.ref.rv.velocity;
  density = iod.ref.rv.density;
  length = iod.ref.rv.length;

// Included (MB)
   if (iod.eqs.fluidModel.fluid == FluidModelData::PERFECT_GAS)
     completeJac = true;
   else
     completeJac = false;

}

//------------------------------------------------------------------------------

void FemEquationTermKE::computeLaminarTransportCoefficients(
  const double T, double &mul, double &lambdal, double &kappal)
{

  mul     = viscoFcn->compute_mu(T);
  lambdal = viscoFcn->compute_lambda(T,mul);
  kappal  = thermalCondFcn->compute(T);

}

//------------------------------------------------------------------------------

void FemEquationTermKE::computeTurbulentTransportCoefficients(
  double *V[], int nodeNum[], SVec<double,3> &X,
  const double mul, const double lambdal, const double kappal,
  double &rhok, double &rhoeps, double &mut, double &lambdat, double &kappat)
{

  mut = computeTurbulentViscosity(V, rhok, rhoeps);

  //Applying the laminar-turbulent trip
  if(trip) {
    int in_trip = 0;
    for (int k=0; k<3; k++)
      if (X[nodeNum[k]][0]>=x0 && X[nodeNum[k]][0]<=x1 &&
          X[nodeNum[k]][1]>=y0 && X[nodeNum[k]][1]<=y1 &&
	  X[nodeNum[k]][2]>=z0 && X[nodeNum[k]][2]<=z1 )
        in_trip = 1;

    if (!in_trip) mut = 0.0;
  }
  lambdat = computeSecondTurbulentViscosity(lambdal, mul, mut);
  kappat  = turbThermalCondFcn.turbulentConductivity(mut);

}

//------------------------------------------------------------------------------

double FemEquationTermKE::computeViscousTimeStep(double X[3], double *V)
{

  double T;
  computeTemperature(V,T);
  double mul = viscoFcn->compute_mu(T);
  double mut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1)
       mut = computeTurbulentViscosity(V);
    else
       mut = 0.0;
  }
  else
    mut = computeTurbulentViscosity(V);

  return ooreynolds_mu * (mul+mut)/V[0];

}

//------------------------------------------------------------------------------

// Included (MB)
double FemEquationTermKE::computeDerivativeOfViscousTimeStep(double X[3], double dX[3], double *V, double *dV, double dMach)
{

  double T, dT;
  computeTemperature(V,T);
  computeDerivativeOfTemperature(V,dV,dT);
  double mul = viscoFcn->compute_mu(T);
  double dmul = viscoFcn->compute_muDerivative(T,dT,dMach);
  double mut, dmut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1) {
       mut = computeTurbulentViscosity(V);
       dmut = computeDerivativeOfTurbulentViscosity(V,dV,dMach);
    }
    else {
       mut = 0.0;
       dmut = 0.0;
    }
  }
  else {
    mut = computeTurbulentViscosity(V);
    dmut = computeDerivativeOfTurbulentViscosity(V,dV,dMach);
  }

  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;

  return dooreynolds_mu * (mul+mut)/V[0] + ooreynolds_mu * ((dmul+dmut)*V[0]-(mul+mut)*dV[0])/(V[0]*V[0]);

}

//------------------------------------------------------------------------------

bool FemEquationTermKE::computeVolumeTerm(double dp1dxj[4][3], double d2w[4], 
					  double *V[4], double *r, double *S, double *PR, double tetVol,
                                          SVec<double,3> &X, int nodeNum[4], int material_id)
{
  bool porousmedia = false;

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double dTdxj[3];
  computeTemperatureGradient(dp1dxj, T, dTdxj);

  double mul, lambdal, kappal;
  computeLaminarTransportCoefficients(Tcg, mul, lambdal, kappal);

  double rhok, rhoeps;
  double mut, lambdat, kappat; //may be modified later if element is flagged as a porous media
  computeTurbulentTransportCoefficients(V, nodeNum, X, mul, lambdal, kappal, rhok, rhoeps, mut, lambdat, kappat);

  // compute viscous term for k-epsilon turbulence model equation

  double (*R)[7] = reinterpret_cast<double (*)[7]>(r);

  double muk = ooreynolds_mu * (mul + sigma_k * mut);
  double mueps = ooreynolds_mu * (mul + sigma_eps * mut);

  R[0][5] = muk * (dp1dxj[0][0]*V[0][5] + dp1dxj[1][0]*V[1][5] + 
		   dp1dxj[2][0]*V[2][5] + dp1dxj[3][0]*V[3][5]);
  R[0][6] = mueps * (dp1dxj[0][0]*V[0][6] + dp1dxj[1][0]*V[1][6] + 
		     dp1dxj[2][0]*V[2][6] + dp1dxj[3][0]*V[3][6]);

  R[1][5] = muk * (dp1dxj[0][1]*V[0][5] + dp1dxj[1][1]*V[1][5] + 
		   dp1dxj[2][1]*V[2][5] + dp1dxj[3][1]*V[3][5]);
  R[1][6] = mueps * (dp1dxj[0][1]*V[0][6] + dp1dxj[1][1]*V[1][6] + 
		     dp1dxj[2][1]*V[2][6] + dp1dxj[3][1]*V[3][6]);

  R[2][5] = muk * (dp1dxj[0][2]*V[0][5] + dp1dxj[1][2]*V[1][5] + 
		   dp1dxj[2][2]*V[2][5] + dp1dxj[3][2]*V[3][5]);
  R[2][6] = mueps * (dp1dxj[0][2]*V[0][6] + dp1dxj[1][2]*V[1][6] + 
		     dp1dxj[2][2]*V[2][6] + dp1dxj[3][2]*V[3][6]);

  double div = dudxj[0][0] + dudxj[1][1] + dudxj[2][2];
  double div2 = dudxj[0][0]*dudxj[0][0] + dudxj[1][1]*dudxj[1][1] + dudxj[2][2]*dudxj[2][2];
  double a = dudxj[0][1] + dudxj[1][0];
  double b = dudxj[0][2] + dudxj[2][0];
  double c = dudxj[1][2] + dudxj[2][1];
  double prod = ooreynolds_mu * mut * (2.0 * div2 - twothird * div*div + a*a + b*b + c*c)
    - twothird * rhok * div;

  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;
  S[5] = - rhoeps + prod;
  S[6] = (sigma_eps1 * rhoeps * prod - sigma_eps2 * rhoeps*rhoeps) / rhok;

  // initialize porous term
  for (int j=0; j<3*4; ++j) PR[j] = 0.0; 

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {
      // if porous media with material_id has been defined in the input file

      mut     = computePorousTurbulentViscosity(it, ucg, length);
      lambdat = computeSecondPorousTurbulentViscosity(lambdal, mul, mut);
      kappat  = turbThermalCondFcn.turbulentConductivity(mut);

      // porous term
      porousmedia = computeVolumeTermPorousCore(tetVol, it, length, density, velocity, ucg, V, PR);

    }
  }

  // In all cases, porous media or not, compute viscous term for Navier-Stokes equations
  double mu, lambda, kappa;
  mu      = ooreynolds_mu * (mul + mut);
  lambda  = ooreynolds_mu * (lambdal + lambdat);
  kappa   = ooreynolds_mu * (kappal + kappat);
  computeVolumeTermNS(mu, lambda, kappa, ucg, dudxj, dTdxj, R);

  return(porousmedia);

}

//------------------------------------------------------------------------------

// Included (MB)
bool FemEquationTermKE::computeDerivativeOfVolumeTerm(double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4],
					  double *V[4], double *dV[4], double dMach, double *dr, double *dS, double *dPR, double dtetVol, SVec<double,3> &X,
                                          int nodeNum[4], int material_id)
{

  bool porousmedia = false;

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double du[4][3], ducg[3];
  computeDerivativeOfVelocity(dV, du, ducg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dT[4], dTcg;
  computeDerivativeOfTemperature(V, dV, dT, dTcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double ddudxj[3][3];
  computeDerivativeOfVelocityGradient(dp1dxj, ddp1dxj, u, du, ddudxj);

  double dTdxj[3];
  computeTemperatureGradient(dp1dxj, T, dTdxj);

  double ddTdxj[3];
  computeDerivativeOfTemperatureGradient(dp1dxj, ddp1dxj, T, dT, ddTdxj);

  double mul      = viscoFcn->compute_mu(Tcg);
  double dmul     = viscoFcn->compute_muDerivative(Tcg, dTcg, dMach);
  double lambdal  = viscoFcn->compute_lambda(Tcg, mul);
  double dlambdal = viscoFcn->compute_lambdaDerivative(mul, dmul, dMach);
  double kappal   = thermalCondFcn->compute(Tcg);
  double dkappal  = thermalCondFcn->computeDerivative(Tcg, dTcg, dMach);
  double rhok, rhoeps;
  double mut, lambdat;

  double drhok, drhoeps;
  double dmut, dlambdat;

  mut  = computeTurbulentViscosity(V, rhok, rhoeps);
  dmut = computeDerivativeOfTurbulentViscosity(V, dV, drhok, drhoeps, dMach);     

  // Applying the laminar-turbulent trip
  if(trip) {
    int in_trip = 0;
    for (int k=0; k<4; k++)
      if ( X[nodeNum[k]][0]>=x0 && X[nodeNum[k]][0]<=x1 && 
          X[nodeNum[k]][1]>=y0 && X[nodeNum[k]][1]<=y1 &&
          X[nodeNum[k]][2]>=z0 && X[nodeNum[k]][2]<=z1 )
        in_trip = 1;

    if (!in_trip)
      mut = dmut = 0.0; 
  }
  lambdat  = computeSecondTurbulentViscosity(lambdal, mul, mut);
  dlambdat = computeDerivativeOfSecondTurbulentViscosity(lambdal, dlambdal, mul, dmul, mut, dmut);

  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;

  double mu, lambda;
  double dmu, dlambda;
  double kappa;
  double dkappa;

  double (*dR)[7] = reinterpret_cast<double (*)[7]>(dr);

  double muk = ooreynolds_mu * (mul + sigma_k * mut);
  double dmuk = dooreynolds_mu * (mul + sigma_k * mut) + ooreynolds_mu * (dmul + sigma_k * dmut);
  double mueps = ooreynolds_mu * (mul + sigma_eps * mut);
  double dmueps = dooreynolds_mu * (mul + sigma_eps * mut) + ooreynolds_mu * (dmul + sigma_eps * dmut);

  dR[0][5] = dmuk * (dp1dxj[0][0]*V[0][5] + dp1dxj[1][0]*V[1][5] +
		   dp1dxj[2][0]*V[2][5] + dp1dxj[3][0]*V[3][5]) +
           muk * (ddp1dxj[0][0]*V[0][5] + dp1dxj[0][0]*dV[0][5] + ddp1dxj[1][0]*V[1][5] + dp1dxj[1][0]*dV[1][5] +
		   ddp1dxj[2][0]*V[2][5] + dp1dxj[2][0]*dV[2][5] + ddp1dxj[3][0]*V[3][5] + dp1dxj[3][0]*dV[3][5]);
  dR[0][6] = dmueps * (dp1dxj[0][0]*V[0][6] + dp1dxj[1][0]*V[1][6] +
		     dp1dxj[2][0]*V[2][6] + dp1dxj[3][0]*V[3][6]) +
             mueps * (ddp1dxj[0][0]*V[0][6] + dp1dxj[0][0]*dV[0][6] + ddp1dxj[1][0]*V[1][6] + dp1dxj[1][0]*dV[1][6] +
		     ddp1dxj[2][0]*V[2][6] + dp1dxj[2][0]*dV[2][6] + ddp1dxj[3][0]*V[3][6] + dp1dxj[3][0]*dV[3][6]);

  dR[1][5] = dmuk * (dp1dxj[0][1]*V[0][5] + dp1dxj[1][1]*V[1][5] +
		   dp1dxj[2][1]*V[2][5] + dp1dxj[3][1]*V[3][5]) +
           muk * (ddp1dxj[0][1]*V[0][5] + dp1dxj[0][1]*dV[0][5] + ddp1dxj[1][1]*V[1][5] + dp1dxj[1][1]*dV[1][5] +
		   ddp1dxj[2][1]*V[2][5] + dp1dxj[2][1]*dV[2][5] + ddp1dxj[3][1]*V[3][5] + dp1dxj[3][1]*dV[3][5]);
  dR[1][6] = dmueps * (dp1dxj[0][1]*V[0][6] + dp1dxj[1][1]*V[1][6] +
		     dp1dxj[2][1]*V[2][6] + dp1dxj[3][1]*V[3][6]) +
             mueps * (ddp1dxj[0][1]*V[0][6] + dp1dxj[0][1]*dV[0][6] + ddp1dxj[1][1]*V[1][6] + dp1dxj[1][1]*dV[1][6] +
		     ddp1dxj[2][1]*V[2][6] + dp1dxj[2][1]*dV[2][6] + ddp1dxj[3][1]*V[3][6] + dp1dxj[3][1]*dV[3][6]);

  dR[2][5] = dmuk * (dp1dxj[0][2]*V[0][5] + dp1dxj[1][2]*V[1][5] +
		   dp1dxj[2][2]*V[2][5] + dp1dxj[3][2]*V[3][5]) +
           muk * (ddp1dxj[0][2]*V[0][5] + dp1dxj[0][2]*dV[0][5] + ddp1dxj[1][2]*V[1][5] + dp1dxj[1][2]*dV[1][5] +
		   ddp1dxj[2][2]*V[2][5] + dp1dxj[2][2]*dV[2][5] + ddp1dxj[3][2]*V[3][5] + dp1dxj[3][2]*dV[3][5]);
  dR[2][6] = dmueps * (dp1dxj[0][2]*V[0][6] + dp1dxj[1][2]*V[1][6] +
		     dp1dxj[2][2]*V[2][6] + dp1dxj[3][2]*V[3][6]) +
             mueps * (ddp1dxj[0][2]*V[0][6] + dp1dxj[0][2]*dV[0][6] + ddp1dxj[1][2]*V[1][6] + dp1dxj[1][2]*dV[1][6] +
		     ddp1dxj[2][2]*V[2][6] + dp1dxj[2][2]*dV[2][6] + ddp1dxj[3][2]*V[3][6] + dp1dxj[3][2]*dV[3][6]);

  double div = dudxj[0][0] + dudxj[1][1] + dudxj[2][2];
  double ddiv = ddudxj[0][0] + ddudxj[1][1] + ddudxj[2][2];
  double div2 = dudxj[0][0]*dudxj[0][0] + dudxj[1][1]*dudxj[1][1] + dudxj[2][2]*dudxj[2][2];
  double ddiv2 = 2.0*dudxj[0][0]*ddudxj[0][0] + 2.0*dudxj[1][1]*ddudxj[1][1] + 2.0*dudxj[2][2]*ddudxj[2][2];
  double a = dudxj[0][1] + dudxj[1][0];
  double da = ddudxj[0][1] + ddudxj[1][0];
  double b = dudxj[0][2] + dudxj[2][0];
  double db = ddudxj[0][2] + ddudxj[2][0];
  double c = dudxj[1][2] + dudxj[2][1];
  double dc = ddudxj[1][2] + ddudxj[2][1];
  double prod = ooreynolds * mut * (2.0 * div2 - twothird * div*div + a*a + b*b + c*c)
    - twothird * rhok * div;

  double dprod = dooreynolds_mu * mut * (2.0 * div2 - twothird * div*div + a*a + b*b + c*c)
                 - twothird * rhok * div + ooreynolds * dmut * (2.0 * div2 - twothird * div*div + a*a + b*b + c*c) +
                 ooreynolds_mu * mut * (2.0 * ddiv2 - twothird * 2.0*div*ddiv + 2.0*a*da + 2.0*b*db + 2.0*c*dc)
                 - twothird * drhok * div - twothird * rhok * ddiv;

  dS[0] = 0.0;
  dS[1] = 0.0;
  dS[2] = 0.0;
  dS[3] = 0.0;
  dS[4] = 0.0;
  dS[5] = - drhoeps + dprod;
  dS[6] = ( (sigma_eps1 * drhoeps * prod + sigma_eps1 * rhoeps * dprod  - sigma_eps2 * 2.0*rhoeps*drhoeps) * rhok -  (sigma_eps1 * rhoeps * prod - sigma_eps2 * rhoeps*rhoeps) * drhok ) / ( rhok * rhok );

  // Initialize dPR (porous media term)
  for (int j=0; j<3*4; ++j) dPR[j] = 0.0; 

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file
      porousmedia = true;
      fprintf(stderr, "***** Inside the file FemEquationTermDesc.C the derivative related to porus media is not implemented *****\n");
      exit(1);

    }

  }

  // If element wasn't flagged as porous media, treat it as standard fluid 
  if (!porousmedia) {
    mu  =  ooreynolds_mu * (mul + mut);
    dmu = dooreynolds_mu * (mul + mut) + ooreynolds_mu * (dmul + dmut);

    lambda  = ooreynolds_mu * (lambdal + lambdat);
    dlambda = dooreynolds_mu * (lambdal + lambdat) + ooreynolds_mu * (dlambdal + dlambdat);

    kappa  = ooreynolds_mu * (kappal + turbThermalCondFcn.turbulentConductivity(mut));
    dkappa =
        dooreynolds_mu * (kappal + turbThermalCondFcn.turbulentConductivity(mut)) +
        ooreynolds_mu  * (dkappal + turbThermalCondFcn.turbulentConductivityDerivative(dmut));

    computeDerivativeOfVolumeTermNS(mu, dmu, lambda, dlambda, kappa, dkappa, 
	ucg, ducg, dudxj, ddudxj, dTdxj, ddTdxj, dR);
  }

  return(porousmedia);

}

//------------------------------------------------------------------------------

bool FemEquationTermKE::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						  double *V[4], double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                  SVec<double,3> &X, int nodeNum[4], int material_id)
{

  bool porousmedia = false;

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double mul, lambdal, kappal;
  computeLaminarTransportCoefficients(Tcg, mul, lambdal, kappal);

  double rhok, rhoeps;
  double mut, lambdat, kappat; //may be modified later if element is flagged as a porous media
  computeTurbulentTransportCoefficients(V, nodeNum, X, mul, lambdal, kappal, rhok, rhoeps, mut, lambdat, kappat);

  int k;
  for (k=0; k<4*3*7*7; ++k)
    drdu[k] = 0.0;
  for (k=0; k<4*7*7; ++k)
    dsdu[k] = 0.0;
  for (k=0; k<4*4*6*6; ++k)
    dpdu[k] = 0.0;

  double (*dRdU)[3][7][7] = reinterpret_cast<double (*)[3][7][7]>(drdu);
  double (*dSdU)[7][7] = reinterpret_cast<double (*)[7][7]>(dsdu);
  double (*dPdU)[4][6][6] = reinterpret_cast<double (*)[4][6][6]>(dpdu);

  computeJacobianVolumeTermKE<7,5>(dp1dxj, mul, mut, rhok, rhoeps, V, dRdU, dSdU);

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {
      // if porous media with material_id has been defined in the input file
      //
      mut     = computePorousTurbulentViscosity(it, ucg, length);
      lambdat = computeSecondPorousTurbulentViscosity(lambdal, mul, mut);
      kappat  = turbThermalCondFcn.turbulentConductivity(mut);

      // jacobian of the porous term
      porousmedia = computeJacobianVolumeTermPorousCore(tetVol, it, length, density, velocity, ucg, V, dPdU);

    }
  }

  // In all cases, porous media or not, approximate jacobian of viscous term is computed
  // (approximate = transport coefficients are frozen)
  double mu, lambda, kappa;
  mu      = ooreynolds_mu * (mul + mut);
  lambda  = ooreynolds_mu * (lambdal + lambdat);
  kappa   = ooreynolds_mu * (kappal + kappat);
  computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);

  return (porousmedia);
}

//------------------------------------------------------------------------------

void FemEquationTermKE::computeSurfaceTerm(int code, Vec3D &n, double d2w[3], 
					   double *Vwall, double *V[3], double *R)
{
  
  wallFcn->computeSurfaceTerm(code, n, d2w, Vwall, V, R);

  R[5] = 0.0;
  R[6] = 0.0;

}
//------------------------------------------------------------------------------

double FemEquationTermKE::computeNormDerivWallFcn(double rho, double T, double Du1, 
 																double DT1, double d2w, 
																double &dudn, double &dTdn)
{
	fprintf(stderr, "FemEquationTermKE::computeNormDerivWallFcn not implemented");
	return 0.0;

 	//double ut = wallFcn->computedudT(rho, T, Du1, DT1, d2w, dudn, dTdn);
	//return ut;
}

//------------------------------------------------------------------------------

// Included (MB)
void FemEquationTermKE::computeDerivativeOfSurfaceTerm(int code, Vec3D &n, Vec3D &dn, double d2w[3],
					   double *Vwall, double *dVwall, double *V[3], double *dV[3], double dMach, double *dR)
{

  wallFcn->computeDerivativeOfSurfaceTerm(code, n, dn, d2w, Vwall, dVwall, V, dV, dMach, dR);

  dR[5] = 0.0;
  dR[6] = 0.0;

}

//------------------------------------------------------------------------------

void FemEquationTermKE::computeJacobianSurfaceTerm(int code, Vec3D &n, 
						   double d2w[3], double *Vwall, 
						   double *V[3], double *drdu)
{

  for (int k=0; k<3*7*7; ++k)
    drdu[k] = 0.0;

  double (*dRdU)[7][7] = reinterpret_cast<double (*)[7][7]>(drdu);
  wallFcn->computeJacobianSurfaceTerm(code, n, d2w, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

FemEquationTermKEmean::FemEquationTermKEmean(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), KEpsilonTerm(iod), FemEquationTerm(iod.volumes.volumeMap.dataMap),
  turbThermalCondFcn(iod, viscoFcn, vf)
{

  wallFcn = new WallFcn(iod, varFcn, viscoFcn);

  x0 =  iod.eqs.tc.tr.bfix.x0;
  x1 =  iod.eqs.tc.tr.bfix.x1;
  y0 =  iod.eqs.tc.tr.bfix.y0;
  y1 =  iod.eqs.tc.tr.bfix.y1;
  z0 =  iod.eqs.tc.tr.bfix.z0;
  z1 =  iod.eqs.tc.tr.bfix.z1;

  if(x0>x1 || y0>y1 || z0>z1)  trip = 0;
  else    trip = 1;

  velocity = iod.ref.rv.velocity;
  density = iod.ref.rv.density;
  length = iod.ref.rv.length;

}

//------------------------------------------------------------------------------

void FemEquationTermKEmean::computeLaminarTransportCoefficients(
  const double T, double &mul, double &lambdal, double &kappal)
{

  mul     = viscoFcn->compute_mu(T);
  lambdal = viscoFcn->compute_lambda(T,mul);
  kappal  = thermalCondFcn->compute(T);

}

//------------------------------------------------------------------------------

void FemEquationTermKEmean::computeTurbulentTransportCoefficients(
  double *V[], int nodeNum[], SVec<double,3> &X,
  const double mul, const double lambdal, const double kappal,
  double &rhok, double &rhoeps, double &mut, double &lambdat, double &kappat)
{

  mut = computeTurbulentViscosity(V, rhok, rhoeps);

  //Applying the laminar-turbulent trip
  if(trip) {
    int in_trip = 0;
    for (int k=0; k<3; k++)
      if (X[nodeNum[k]][0]>=x0 && X[nodeNum[k]][0]<=x1 &&
          X[nodeNum[k]][1]>=y0 && X[nodeNum[k]][1]<=y1 &&
	  X[nodeNum[k]][2]>=z0 && X[nodeNum[k]][2]<=z1 )
        in_trip = 1;

    if (!in_trip) mut = 0.0;
  }
  lambdat = computeSecondTurbulentViscosity(lambdal, mul, mut);
  kappat  = turbThermalCondFcn.turbulentConductivity(mut);

}

//------------------------------------------------------------------------------

double FemEquationTermKEmean::computeViscousTimeStep(double X[3], double *V)
{

  double T;
  computeTemperature(V,T);
  double mul = viscoFcn->compute_mu(T);
  double mut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1)
       mut = computeTurbulentViscosity(V);
    else
       mut = 0.0;
  }
  else
    mut = computeTurbulentViscosity(V);

  return ooreynolds_mu * (mul+mut)/V[0];

}

//------------------------------------------------------------------------------

// Included (MB)
double FemEquationTermKEmean::computeDerivativeOfViscousTimeStep(double X[3], double dX[3], double *V, double *dV, double dMach)
{

  double T,dT;
  computeTemperature(V,T);
  computeDerivativeOfTemperature(V,dV,dT);
  double mul = viscoFcn->compute_mu(T);
  double dmul = viscoFcn->compute_muDerivative(T,dT,dMach);
  double mut, dmut;
  if(trip){
    if(X[0]>x0 && X[0]<x1 &&
       X[1]>y0 && X[1]<y1 &&
       X[2]>z0 && X[2]<z1) {
       mut = computeTurbulentViscosity(V);
       dmut = computeDerivativeOfTurbulentViscosity(V,dV,dMach);
    }
    else {
       mut = 0.0;
       dmut = 0.0;
    }
  }
  else {
    mut = computeTurbulentViscosity(V);
    dmut = computeDerivativeOfTurbulentViscosity(V,dV,dMach);
  }

  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;

  return dooreynolds_mu * (mul+mut)/V[0] + ooreynolds_mu * ((dmul+dmut)*V[0]-(mul+mut)*dV[0])/(V[0]*V[0]);

}

//------------------------------------------------------------------------------

bool FemEquationTermKEmean::computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], 
						      double *V[4], double *drdu, double *dSdU, double *dpdu, double tetVol,
                                                      SVec<double,3> &X, int nodeNum[4], int material_id)
{

  bool porousmedia = false;
                                                                                                                                                                   
  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double mul, lambdal, kappal;
  computeLaminarTransportCoefficients(Tcg, mul, lambdal, kappal);

  double rhok, rhoeps;
  double mut, lambdat, kappat; //may be modified later if element is flagged as a porous media
  computeTurbulentTransportCoefficients(V, nodeNum, X, mul, lambdal, kappal, rhok, rhoeps, mut, lambdat, kappat);

  double (*dRdU)[3][5][5] = reinterpret_cast<double (*)[3][5][5]>(drdu);
  double (*dPdU)[4][5][5] = reinterpret_cast<double (*)[4][5][5]>(dpdu);

  if(material_id>0) {
    map<int,PorousMedia *>::iterator it = volInfo.find(material_id);
    if(it!=  volInfo.end()) {     // if porous media with material_id has been defined in the input file

      mut     = computePorousTurbulentViscosity(it, ucg, length);
      lambdat = computeSecondPorousTurbulentViscosity(lambdal, mul, mut);
      kappat  = turbThermalCondFcn.turbulentConductivity(mut);

      // jacobian for the porous term
      porousmedia = computeJacobianVolumeTermPorousCore(tetVol, it, length, density, velocity, ucg, V, dPdU);

    }
  }

  // In all cases, porous media or not, approximate jacobian of viscous term is computed
  // (approximate = transport coefficients are frozen)
  double mu, lambda, kappa;
  mu      = ooreynolds_mu * (mul + mut);
  lambda  = ooreynolds_mu * (lambdal + lambdat);
  kappa   = ooreynolds_mu * (kappal + kappat);
  computeJacobianVolumeTermNS(dp1dxj, mu, lambda, kappa, V, T, dRdU);

  return (porousmedia);

}

//------------------------------------------------------------------------------

void FemEquationTermKEmean::computeJacobianSurfaceTerm(int code, Vec3D &n, 
						       double d2w[3], double *Vwall, 
						       double *V[3], double *drdu)
{

  double (*dRdU)[5][5] = reinterpret_cast<double (*)[5][5]>(drdu);
  wallFcn->computeJacobianSurfaceTerm(code, n, d2w, Vwall, V, dRdU);

}

//------------------------------------------------------------------------------

FemEquationTermKEturb::FemEquationTermKEturb(IoData &iod, VarFcn *vf) :
  NavierStokesTerm(iod, vf), KEpsilonTerm(iod), FemEquationTerm(iod.volumes.volumeMap.dataMap)
{

}

//------------------------------------------------------------------------------

bool FemEquationTermKEturb::computeJacobianVolumeTerm(double dp1dxj[4][3], 
						      double d2w[4], double *V[4], 
						      double *drdu, double *dsdu, double *dpdu, double tetVol,
                                                      SVec<double,3> &X, int nodeNum[4], int material_id)
{

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double mul = viscoFcn->compute_mu(Tcg);

  double rhok, rhoeps;
  double mut;

  mut = computeTurbulentViscosity(V, rhok, rhoeps);

  double (*dRdU)[3][2][2] = reinterpret_cast<double (*)[3][2][2]>(drdu);
  double (*dSdU)[2][2] = reinterpret_cast<double (*)[2][2]>(dsdu);
  computeJacobianVolumeTermKE<2,0>(dp1dxj, mul, mut, rhok, rhoeps, V, dRdU, dSdU);

  return false;

}

//------------------------------------------------------------------------------
