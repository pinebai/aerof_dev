#include <WaleLESTerm.h>
#include <cmath>

//-----------------------------------------------------------------------

WaleLESTerm::WaleLESTerm(IoData& iod, VarFcn *vf): NavierStokesTerm(iod, vf)
{

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = new WallFcn(iod, varFcn, viscoFcn);

  Cw         = iod.eqs.tc.les.wale.c_w;
  gamma      = iod.eqs.fluidModel.gasModel.specificHeatRatio;
  oogamma1   = 1.0/(gamma - 1.0);
  PrT        = iod.eqs.tc.prandtlTurbulent;
  coeff      = gamma/PrT;
  onethird   = 1.0/3.0;
  onefourth  = 1.0/4.0;
  twothirds  = 2.0/3.0;
  onesixth   = 1.0/6.0;

  x0 =  iod.eqs.tc.tr.bfix.x0;
  x1 =  iod.eqs.tc.tr.bfix.x1;
  y0 =  iod.eqs.tc.tr.bfix.y0;
  y1 =  iod.eqs.tc.tr.bfix.y1;
  z0 =  iod.eqs.tc.tr.bfix.z0;
  z1 =  iod.eqs.tc.tr.bfix.z1;

  if(x0 == 0.0 && x1 == -1.0 && y0 == 0.0 && y1 == -1.0 && z0 == 0.0 && z1 == -1.0){
    trip = 0;}
  else{
    trip = 1;
  }

  if (iod.eqs.tc.les.delta == LESModelData::VOLUME)
    delta = 0; 
  else
    delta = 1;

}

//-----------------------------------------------------------------------

WaleLESTerm::~WaleLESTerm()
{

}

//-----------------------------------------------------------------------

void WaleLESTerm::compute(double tetVol, double dp1dxj[4][3],
				 double *V[4], double *r, SVec<double,3> &X, int nodeNum[4])

{

  double (*R)[5] = reinterpret_cast<double (*)[5]>(r);

  // Compute filter width
  double delta = computeFilterWidth(tetVol, X, nodeNum);

  // Compute Reynolds stress tensor and energy gradient
  double u[4][3];
  double ucg[3];
  double en[4];
  double dudxj[3][3];
  double dedxj[3];
  double tauij[3][3];
  double rhoCG = onefourth * (V[0][0] + V[1][0] +
			      V[2][0] + V[3][0]);
  computeVelocity(V, u, ucg);
  computeVelocityGradient(dp1dxj, u, dudxj);

  double muT;

  // Applying the laminar-turbulent trip
  if(trip){
    if(x0>x1 || y0>y1 || z0>z1){
       fprintf(stderr," ***Error** The values of the extremes in the box for Tripping is wrong\n");
       exit(1);
    }
    else{
       if((X[nodeNum[0]][0]>=x0 && X[nodeNum[0]][0]<=x1 && X[nodeNum[0]][1]>=y0 && X[nodeNum[0]][1]<=y1 &&
          X[nodeNum[0]][2]>=z0 && X[nodeNum[0]][2]<=z1) || (X[nodeNum[1]][0]>=x0 && X[nodeNum[1]][0]<=x1 &&
          X[nodeNum[1]][1]>=y0 && X[nodeNum[1]][1]<=y1 && X[nodeNum[1]][2]>=z0 && X[nodeNum[1]][2]<=z1) ||
          (X[nodeNum[2]][0]>=x0 && X[nodeNum[2]][0]<=x1 && X[nodeNum[2]][1]>=y0 && X[nodeNum[2]][1]<=y1 &&
          X[nodeNum[2]][2]>=z0 && X[nodeNum[2]][2]<=z1) || (X[nodeNum[3]][0]>=x0 && X[nodeNum[3]][0]<=x1 &&
          X[nodeNum[3]][1]>=y0 && X[nodeNum[3]][1]<=y1 && X[nodeNum[3]][2]>=z0 && X[nodeNum[3]][2]<=z1)){
          muT = computeEddyViscosity(rhoCG, delta, dudxj);} 
       else{
          muT = 0.0;
       }
    }
  }
  else{
    muT = computeEddyViscosity(rhoCG, delta, dudxj); 
  }

  computeStressTensor(muT, dudxj, tauij);
  computeEnergy(V, en);
  computeEnergyGradient(dp1dxj, en, dedxj);

  // Now fill R with the appropriate values
  double coeff2 = coeff * muT;
  
  R[0][0] = 0.0;
  R[0][1] = tauij[0][0];
  R[0][2] = tauij[1][0];
  R[0][3] = tauij[2][0];
  R[0][4] = coeff2*dedxj[0];

  R[1][0] = 0.0;
  R[1][1] = tauij[0][1];
  R[1][2] = tauij[1][1];
  R[1][3] = tauij[2][1];
  R[1][4] = coeff2*dedxj[1];

  R[2][0] = 0.0;
  R[2][1] = tauij[0][2];
  R[2][2] = tauij[1][2];
  R[2][3] = tauij[2][2];
  R[2][4] = coeff2*dedxj[2];

}

//-----------------------------------------------------------------------

double WaleLESTerm::computeEddyViscosity(double rho, double Delta, double duidxj[3][3])
{


  /* Reference: Subgrid Scale Stress Modelling Based on the Square of 
                the Velocity Gradient Tensor, F.Nicoud and F. Ducros
		Flow, Turbulence and Combustion 62, 183-200, 1999   */


  // strain rate tensor

  double S[3][3];

  S[0][0] = duidxj[0][0];
  S[1][1] = duidxj[1][1];
  S[2][2] = duidxj[2][2];

  S[0][1] = 0.5 * (duidxj[0][1] + duidxj[1][0]);
  S[0][2] = 0.5 * (duidxj[0][2] + duidxj[2][0]);
  S[1][2] = 0.5 * (duidxj[1][2] + duidxj[2][1]);

  S[1][0] = S[0][1];
  S[2][0] = S[0][2];
  S[2][1] = S[1][2];

  double S2 = (S[0][0]*S[0][0] + 
	       S[0][1]*S[0][1] +
	       S[0][2]*S[0][2] +

	       S[1][0]*S[1][0] +
	       S[1][1]*S[1][1] +
	       S[1][2]*S[1][2] +

	       S[2][0]*S[2][0] +
	       S[2][1]*S[2][1] +
	       S[2][2]*S[2][2]);

  double sqS2 = S2 * S2;

  // computing Aij = Sij * Sjk

  double A[3][3];

  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      A[i][j] = 0.0;

  
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      for (int k=0; k<3; ++k)
        A[i][j] += S[i][k] * S[k][j];


  // vorticity tensor

  double O[3][3];

  O[0][0] = 0.0;
  O[1][1] = 0.0;  
  O[2][2] = 0.0;  

  O[0][1] = 0.5 * (duidxj[0][1] - duidxj[1][0]);
  O[0][2] = 0.5 * (duidxj[0][2] - duidxj[2][0]);
  O[1][2] = 0.5 * (duidxj[1][2] - duidxj[2][1]);

  O[1][0] = -O[0][1];
  O[2][0] = -O[0][2];
  O[2][1] = -O[1][2];

  double O2 = (O[0][0]*O[0][0] +
               O[0][1]*O[0][1] +
               O[0][2]*O[0][2] +
                                                                                                                                               
               O[1][0]*O[1][0] +
               O[1][1]*O[1][1] +
               O[1][2]*O[1][2] +
                                                                                                                                               
               O[2][0]*O[2][0] +
               O[2][1]*O[2][1] +
               O[2][2]*O[2][2]);

  double sqO2 = O2 * O2;

  double S2O2 = S2*O2;

 // computing Bij = Oij * Ojk

  double B[3][3];

  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      B[i][j] = 0.0;

  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      for (int k=0; k<3; ++k)
        B[i][j] += O[i][k] * O[k][j];

  // computing FSO = Sij*Sjk*Okl*Oli

  double FSO = A[0][0]*B[0][0] + A[1][1]*B[1][1] + A[2][2]*B[2][2] +
               A[0][1]*B[1][0] + A[0][2]*B[2][0] +
               A[1][0]*B[0][1] + A[1][2]*B[2][1] +
               A[2][0]*B[0][2] + A[2][1]*B[1][2];

  // computing second invariant SD

  double SD = onesixth*(sqS2 + sqO2) + twothirds*S2O2 + 2*FSO;

  double num, denom;

  if (fabs(SD) < 1e-10) {
    num = 0.0;
    denom =  1.0;
  }
  else if (S2 < 1e-10) {
    num = pow(fabs(SD),1.5);
    denom = pow(fabs(SD),1.25);
  }
  else {
    num = pow(fabs(SD),1.5);
    denom = pow(S2,2.5) + pow(fabs(SD),1.25);
  }

  return rho*(Cw*Delta)*(Cw*Delta)*(num/denom);

}

//-----------------------------------------------------------------------

double WaleLESTerm::computeFilterWidth(double tetVol, SVec<double,3> &X, int nodeNum[4])

{

  if (delta){ 
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
    return maxl;
  }
  else {
    return pow(tetVol, onethird);
  }

}

//-----------------------------------------------------------------------

void WaleLESTerm::computeVelocity(double *V[4], double u[4][3],
					 double ucg[3])
{

  u[0][0] = V[0][1];
  u[0][1] = V[0][2];
  u[0][2] = V[0][3];

  u[1][0] = V[1][1];
  u[1][1] = V[1][2];
  u[1][2] = V[1][3];

  u[2][0] = V[2][1];
  u[2][1] = V[2][2];
  u[2][2] = V[2][3];

  u[3][0] = V[3][1];
  u[3][1] = V[3][2];
  u[3][2] = V[3][3];

  ucg[0] = onefourth * (u[0][0] + u[1][0] + u[2][0] + u[3][0]);
  ucg[1] = onefourth * (u[0][1] + u[1][1] + u[2][1] + u[3][1]);
  ucg[2] = onefourth * (u[0][2] + u[1][2] + u[2][2] + u[3][2]);

}

//-----------------------------------------------------------------------

void WaleLESTerm::computeVelocityGradient(double dp1dxj[4][3],
					      double u[4][3],
					      double dudxj[3][3])

{

  dudxj[0][0] = dp1dxj[0][0]*u[0][0] + dp1dxj[1][0]*u[1][0] + 
    dp1dxj[2][0]*u[2][0] + dp1dxj[3][0]*u[3][0];

  dudxj[0][1] = dp1dxj[0][1]*u[0][0] + dp1dxj[1][1]*u[1][0] + 
    dp1dxj[2][1]*u[2][0] + dp1dxj[3][1]*u[3][0];

  dudxj[0][2] = dp1dxj[0][2]*u[0][0] + dp1dxj[1][2]*u[1][0] + 
    dp1dxj[2][2]*u[2][0] + dp1dxj[3][2]*u[3][0];

  dudxj[1][0] = dp1dxj[0][0]*u[0][1] + dp1dxj[1][0]*u[1][1] + 
    dp1dxj[2][0]*u[2][1] + dp1dxj[3][0]*u[3][1];

  dudxj[1][1] = dp1dxj[0][1]*u[0][1] + dp1dxj[1][1]*u[1][1] + 
    dp1dxj[2][1]*u[2][1] + dp1dxj[3][1]*u[3][1];

  dudxj[1][2] = dp1dxj[0][2]*u[0][1] + dp1dxj[1][2]*u[1][1] + 
    dp1dxj[2][2]*u[2][1] + dp1dxj[3][2]*u[3][1];

  dudxj[2][0] = dp1dxj[0][0]*u[0][2] + dp1dxj[1][0]*u[1][2] + 
    dp1dxj[2][0]*u[2][2] + dp1dxj[3][0]*u[3][2];

  dudxj[2][1] = dp1dxj[0][1]*u[0][2] + dp1dxj[1][1]*u[1][2] + 
    dp1dxj[2][1]*u[2][2] + dp1dxj[3][1]*u[3][2];

  dudxj[2][2] = dp1dxj[0][2]*u[0][2] + dp1dxj[1][2]*u[1][2] + 
    dp1dxj[2][2]*u[2][2] + dp1dxj[3][2]*u[3][2];

}

//-----------------------------------------------------------------------

void WaleLESTerm::computeStressTensor(double mu, double dudxj[3][3],
					     double tij[3][3])

{

  tij[0][0] = mu * twothirds * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);
  tij[1][1] = mu * twothirds * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);
  tij[2][2] = mu * twothirds * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);
  tij[0][1] = mu * (dudxj[1][0] + dudxj[0][1]);
  tij[0][2] = mu * (dudxj[2][0] + dudxj[0][2]);
  tij[1][2] = mu * (dudxj[2][1] + dudxj[1][2]);
  tij[1][0] = tij[0][1];
  tij[2][0] = tij[0][2];
  tij[2][1] = tij[1][2];

}

//-----------------------------------------------------------------------

void WaleLESTerm::computeEnergy(double *V[4], double en[4])
{

  en[0] = (V[0][4] / V[0][0]) * oogamma1;
  en[1] = (V[1][4] / V[1][0]) * oogamma1;
  en[2] = (V[2][4] / V[2][0]) * oogamma1;
  en[3] = (V[3][4] / V[3][0]) * oogamma1;

}

//-----------------------------------------------------------------------

void WaleLESTerm::computeEnergyGradient(double dp1dxj[4][3],
					       double en[4],
					       double dedxj[3])

{

  dedxj[0] = (dp1dxj[0][0]*en[0] +
	      dp1dxj[1][0]*en[1] +
	      dp1dxj[2][0]*en[2] +
	      dp1dxj[3][0]*en[3]);

  dedxj[1] = (dp1dxj[0][1]*en[0]+
	      dp1dxj[1][1]*en[1]+
	      dp1dxj[2][1]*en[2]+
	      dp1dxj[3][1]*en[3]);

  dedxj[2] = (dp1dxj[0][2]*en[0]+
	      dp1dxj[1][2]*en[1]+
	      dp1dxj[2][2]*en[2]+
	      dp1dxj[3][2]*en[3]);
}

//-----------------------------------------------------------------------

double WaleLESTerm::computeMutOMu(double tetVol, double dp1dxj[4][3],
                                         double *V[4], SVec<double,3> &X, 
                                         int nodeNum[4])

{

  // Compute filter width
  double delta = computeFilterWidth(tetVol, X, nodeNum);

  // Compute Reynolds stress tensor and energy gradient
  double u[4][3];
  double ucg[3];
  double dudxj[3][3];
  double rhoCG = onefourth * (V[0][0] + V[1][0] +
                              V[2][0] + V[3][0]);
  computeVelocity(V, u, ucg);
  computeVelocityGradient(dp1dxj, u, dudxj);

  double muT;

  // Applying the laminar-turbulent trip
  if(trip){
    if(x0>x1 || y0>y1 || z0>z1){
       fprintf(stderr," ***Error** The values of the extremes in the box for Tripping is wrong\n");
       exit(1);
    }
    else{
       if((X[nodeNum[0]][0]>=x0 && X[nodeNum[0]][0]<=x1 && X[nodeNum[0]][1]>=y0 && X[nodeNum[0]][1]<=y1 &&
          X[nodeNum[0]][2]>=z0 && X[nodeNum[0]][2]<=z1) || (X[nodeNum[1]][0]>=x0 && X[nodeNum[1]][0]<=x1 &&
          X[nodeNum[1]][1]>=y0 && X[nodeNum[1]][1]<=y1 && X[nodeNum[1]][2]>=z0 && X[nodeNum[1]][2]<=z1) ||
          (X[nodeNum[2]][0]>=x0 && X[nodeNum[2]][0]<=x1 && X[nodeNum[2]][1]>=y0 && X[nodeNum[2]][1]<=y1 &&
          X[nodeNum[2]][2]>=z0 && X[nodeNum[2]][2]<=z1) || (X[nodeNum[3]][0]>=x0 && X[nodeNum[3]][0]<=x1 &&
          X[nodeNum[3]][1]>=y0 && X[nodeNum[3]][1]<=y1 && X[nodeNum[3]][2]>=z0 && X[nodeNum[3]][2]<=z1)){
          muT = computeEddyViscosity(rhoCG, delta, dudxj);}
       else{
          muT = 0.0;
       }
    }
  }
  else{
    muT = computeEddyViscosity(rhoCG, delta, dudxj);
  }

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);
  double mu = ooreynolds_mu * viscoFcn->compute_mu(Tcg);

  return (muT/mu);

}

//-----------------------------------------------------------------------
