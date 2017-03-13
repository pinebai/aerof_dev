#include <DynamicLESTerm.h>

#include <cmath>

//-----------------------------------------------------------------------

DynamicLESTerm::DynamicLESTerm(IoData& iod, VarFcn *vf): NavierStokesTerm(iod, vf)
{

                                                                                                                                                    
  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = new WallFcn(iod, varFcn, viscoFcn);

  gamma      = iod.eqs.fluidModel.gasModel.specificHeatRatio;
  R          = iod.eqs.fluidModel.gasModel.idealGasConstant; 
  oogamma1   = 1.0/(gamma - 1.0);
  ooR        = 1.0/R;
  
  onethird   = 1.0/3.0;
  onefourth  = 1.0/4.0;
  twothirds  = 2.0/3.0;

  csmax = iod.eqs.tc.les.dles.clip.cs_max;
  ptmin = iod.eqs.tc.les.dles.clip.pt_min;
  ptmax = iod.eqs.tc.les.dles.clip.pt_max;
  Prandtl = iod.eqs.tc.prandtlTurbulent;

  mach = iod.bc.inlet.mach;

  x0 =  iod.eqs.tc.tr.bfix.x0;
  x1 =  iod.eqs.tc.tr.bfix.x1;
  y0 =  iod.eqs.tc.tr.bfix.y0;
  y1 =  iod.eqs.tc.tr.bfix.y1;
  z0 =  iod.eqs.tc.tr.bfix.z0;
  z1 =  iod.eqs.tc.tr.bfix.z1;

  if(x0>x1 || y0>y1 || z0>z1)  trip = 0;
  else    trip = 1;

  if (iod.eqs.tc.les.delta == LESModelData::VOLUME){
    delta = 0; }
  else{
    delta = 1;
  }

}

//-----------------------------------------------------------------------

DynamicLESTerm::~DynamicLESTerm()
{

}

//-----------------------------------------------------------------------

void DynamicLESTerm::computeDelta(double tetVol, SVec<double,3> &X, int nodeNum[4], double &delta)

{

  delta = computeFilterWidth(tetVol, X, nodeNum);

}
                                                                                                                                                       
//-----------------------------------------------------------------------

void DynamicLESTerm::compute(double tetVol, double dp1dxj[4][3], double *V[4], 
                             double Cs[4], double Pt[4], double *r, 
									  SVec<double,3> &X, int nodeNum[4])

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

  csprime = onefourth*(Cs[0]+Cs[1]+Cs[2]+Cs[3]);
  Prsgs   = onefourth*(Pt[0]+Pt[1]+Pt[2]+Pt[3]);
    
  // Clipping for (cs*delta)^2 and pt
  if(csprime > pow(csmax*delta, 2)) 
	  csprime = csmax*delta*delta; 

  if (mach < 0.3) 
	  Prsgs = Prandtl; // for subsonic flow use constant prandtl number
  else {
    if(Prsgs < ptmin) Prsgs = ptmin;
    if(Prsgs > ptmax) Prsgs = ptmax;
  }

  computeVelocity(V, u, ucg);
  computeVelocityGradient(dp1dxj, u, dudxj);

  double muT;

// Applying the laminar-turbulent trip

  if(trip)
  {
     if((X[nodeNum[0]][0]>=x0 && X[nodeNum[0]][0]<=x1 && X[nodeNum[0]][1]>=y0 && X[nodeNum[0]][1]<=y1 &&
        X[nodeNum[0]][2]>=z0 && X[nodeNum[0]][2]<=z1) || (X[nodeNum[1]][0]>=x0 && X[nodeNum[1]][0]<=x1 &&
        X[nodeNum[1]][1]>=y0 && X[nodeNum[1]][1]<=y1 && X[nodeNum[1]][2]>=z0 && X[nodeNum[1]][2]<=z1) ||
        (X[nodeNum[2]][0]>=x0 && X[nodeNum[2]][0]<=x1 && X[nodeNum[2]][1]>=y0 && X[nodeNum[2]][1]<=y1 &&
        X[nodeNum[2]][2]>=z0 && X[nodeNum[2]][2]<=z1) || (X[nodeNum[3]][0]>=x0 && X[nodeNum[3]][0]<=x1 &&
        X[nodeNum[3]][1]>=y0 && X[nodeNum[3]][1]<=y1 && X[nodeNum[3]][2]>=z0 && X[nodeNum[3]][2]<=z1)){
        muT = computeEddyViscosity(rhoCG, delta, dudxj);} 
     else
        muT = 0.0;
     }
  else
    muT = computeEddyViscosity(rhoCG, delta, dudxj); 

// mu + muT check: if mu + muT < 0.0 then muT is set to -mu

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);
  double mu = ooreynolds_mu * viscoFcn->compute_mu(Tcg);

  if((mu + muT) < 0.0) muT = -mu;

// trubulence stress and energy gradient calculation

  computeStressTensor(muT, dudxj, tauij);
  computeEnergy(V, en);
  computeEnergyGradient(dp1dxj, en, dedxj);

// Now fill R with the appropriate values

  double coeff2 = gamma*muT/Prsgs;

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
	
double DynamicLESTerm::computeEddyViscosity(double rho, double Delta,
						double duidxj[3][3])
{

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

  double sqrt2S2 = sqrt(2.0 * S2);

  return (rho)*csprime*sqrt2S2;

}

//-----------------------------------------------------------------------

double DynamicLESTerm::computeFilterWidth(double tetVol, SVec<double,3> &X, int nodeNum[4])

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

void DynamicLESTerm::computeVelocity(double *V[4], double u[4][3],
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

void DynamicLESTerm::computeVelocityGradient(double dp1dxj[4][3],
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

void DynamicLESTerm::computeStressTensor(double mu, double dudxj[3][3],
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

void DynamicLESTerm::computeEnergy(double *V[4], double en[4])
{

  en[0] = (V[0][4] / V[0][0]) * oogamma1;
  en[1] = (V[1][4] / V[1][0]) * oogamma1;
  en[2] = (V[2][4] / V[2][0]) * oogamma1;
  en[3] = (V[3][4] / V[3][0]) * oogamma1;

}

//-----------------------------------------------------------------------

void DynamicLESTerm::computeEnergyGradient(double dp1dxj[4][3],
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

double DynamicLESTerm::computeMutOMu(double tetVol, double dp1dxj[4][3], double *V[4], 
                                   double Cs[4], double Pt[4], SVec<double,3> &X, int nodeNum[4])

{

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

  csprime = onefourth*(Cs[0]+Cs[1]+Cs[2]+Cs[3]);
  Prsgs = onefourth*(Pt[0]+Pt[1]+Pt[2]+Pt[3]);

// Clipping for (cs*delta)^2 and pt

  if (csprime > pow(csmax*delta,2)) csprime = csmax*delta*delta; 
  if (mach < 0.3) Prsgs = Prandtl; // for subsonic flow use constant prandtl number
  else {
    if(Prsgs < ptmin) Prsgs = ptmin;
    if(Prsgs > ptmax) Prsgs = ptmax;
  }

  computeVelocity(V, u, ucg);
  computeVelocityGradient(dp1dxj, u, dudxj);

  double muT;

// Applying the laminar-turbulent trip

  if(trip){
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
  else{
    muT = computeEddyViscosity(rhoCG, delta, dudxj); 
  }

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);
  double mu = ooreynolds * viscoFcn->compute_mu(Tcg);


  return (muT/mu);

}

//-----------------------------------------------------------------------

double DynamicLESTerm::outputCsValues(double tetVol, double Cs[4], SVec<double,3> &X, int nodeNum[4])
{

  double delta = computeFilterWidth(tetVol, X, nodeNum);

  csprime = onefourth*(Cs[0]+Cs[1]+Cs[2]+Cs[3]);

  return (csprime/(delta*delta));

}

//-----------------------------------------------------------------------
