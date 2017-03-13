#include <MacroCell.h>
#include <VMSLESTerm.h>

#include <cmath>
#include <cstdlib>

//------------------------------------------------------------------------

VMSLESTerm::VMSLESTerm(IoData& iod, VarFcn *vf): NavierStokesTerm(iod, vf)
{

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = new WallFcn(iod, varFcn, viscoFcn);

  csprime    = iod.eqs.tc.les.vms.c_s_prime;
  gamma      = iod.eqs.fluidModel.gasModel.specificHeatRatio;
  oogamma1   = 1.0/(gamma - 1.0);
  PrT        = iod.eqs.tc.prandtlTurbulent;
  coeff      = gamma/PrT;
  onethird   = 1.0/3.0;
  onefourth  = 1.0/4.0;
  twothirds  = 2.0/3.0;

  x0 =  iod.eqs.tc.tr.bfix.x0;
  x1 =  iod.eqs.tc.tr.bfix.x1;
  y0 =  iod.eqs.tc.tr.bfix.y0;
  y1 =  iod.eqs.tc.tr.bfix.y1;
  z0 =  iod.eqs.tc.tr.bfix.z0;
  z1 =  iod.eqs.tc.tr.bfix.z1;

  if(x0>x1 || y0>y1 || z0>z1)  trip = 0;
  else         trip = 1;

  if (iod.eqs.tc.les.delta == LESModelData::VOLUME){
    delta = 0; }
  else{
    delta = 1;
  }

}

//------------------------------------------------------------------------

VMSLESTerm::~VMSLESTerm()
{

}

//------------------------------------------------------------------------

void VMSLESTerm::compute(double tetVol, 
			 double dp1dxj[4][3], 
			 double *VBar[4],
			 double *V[4], 
			 double *r, 
			 SVec<double,3> &X,
                         int nodeNum[4])

{

  double (*R)[5] = reinterpret_cast<double (*)[5]>(r);

  // Compute filter width //
  
  double Delta = computeFilterWidth(tetVol, X, nodeNum);

  // Compute small scale Reynolds stress tensor and small scale //
  // energy gradient //
  
  double uPrime[4][3];
  double ePrime[4];
  double duPrimedxj[3][3];
  double dePrimedxj[3];
  double tauPrimeij[3][3];
  getPrimeValues(V, VBar, uPrime, ePrime);
  computeVelocityGradient(dp1dxj, uPrime, duPrimedxj);
  double rhoBarCG = onefourth * (VBar[0][0] + VBar[1][0] +
				 VBar[2][0] + VBar[3][0]);
  double muT;
  
  // Applying the laminar-turbulent trip //

  if(trip){
    if((X[nodeNum[0]][0]>=x0 && X[nodeNum[0]][0]<=x1 && X[nodeNum[0]][1]>=y0 && X[nodeNum[0]][1]<=y1 &&
        X[nodeNum[0]][2]>=z0 && X[nodeNum[0]][2]<=z1) || (X[nodeNum[1]][0]>=x0 && X[nodeNum[1]][0]<=x1 &&
        X[nodeNum[1]][1]>=y0 && X[nodeNum[1]][1]<=y1 && X[nodeNum[1]][2]>=z0 && X[nodeNum[1]][2]<=z1) ||
       (X[nodeNum[2]][0]>=x0 && X[nodeNum[2]][0]<=x1 && X[nodeNum[2]][1]>=y0 && X[nodeNum[2]][1]<=y1 &&
        X[nodeNum[2]][2]>=z0 && X[nodeNum[2]][2]<=z1) || (X[nodeNum[3]][0]>=x0 && X[nodeNum[3]][0]<=x1 &&
        X[nodeNum[3]][1]>=y0 && X[nodeNum[3]][1]<=y1 && X[nodeNum[3]][2]>=z0 && X[nodeNum[3]][2]<=z1)){
          muT = computeEddyViscosity(rhoBarCG, Delta, duPrimedxj);
    }
    else{
       muT = 0.0; 
    }
  }
  else{
    muT = computeEddyViscosity(rhoBarCG, Delta, duPrimedxj);
  }
  
  computeStressTensor(muT, duPrimedxj, tauPrimeij);
  computeSmallEnergyGradient(dp1dxj, ePrime, dePrimedxj);

  // Now fill R with the appropriate values //
  
  double coeff2 = coeff * muT;
  
  R[0][0] = 0.0;
  R[0][1] = tauPrimeij[0][0];
  R[0][2] = tauPrimeij[1][0];
  R[0][3] = tauPrimeij[2][0];
  R[0][4] = coeff2*dePrimedxj[0];

  R[1][0] = 0.0;
  R[1][1] = tauPrimeij[0][1];
  R[1][2] = tauPrimeij[1][1];
  R[1][3] = tauPrimeij[2][1];
  R[1][4] = coeff2*dePrimedxj[1];

  R[2][0] = 0.0;
  R[2][1] = tauPrimeij[0][2];
  R[2][2] = tauPrimeij[1][2];
  R[2][3] = tauPrimeij[2][2];
  R[2][4] = coeff2*dePrimedxj[2];

}

//------------------------------------------------------------------------

double VMSLESTerm::computeEddyViscosity(double rhoBar, double Delta,
					double duidxj[3][3])
{

  double Sprime[3][3];

  Sprime[0][0] = duidxj[0][0];
  Sprime[1][1] = duidxj[1][1];
  Sprime[2][2] = duidxj[2][2];

  Sprime[0][1] = 0.5 * (duidxj[0][1] + duidxj[1][0]);
  Sprime[0][2] = 0.5 * (duidxj[0][2] + duidxj[2][0]);
  Sprime[1][2] = 0.5 * (duidxj[1][2] + duidxj[2][1]);

  Sprime[1][0] = Sprime[0][1];
  Sprime[2][0] = Sprime[0][2];
  Sprime[2][1] = Sprime[1][2];

  double Sprime2 = (Sprime[0][0]*Sprime[0][0] + 
		    Sprime[0][1]*Sprime[0][1] +
		    Sprime[0][2]*Sprime[0][2] +

		    Sprime[1][0]*Sprime[1][0] +
		    Sprime[1][1]*Sprime[1][1] +
		    Sprime[1][2]*Sprime[1][2] +

		    Sprime[2][0]*Sprime[2][0] +
		    Sprime[2][1]*Sprime[2][1] +
		    Sprime[2][2]*Sprime[2][2]);

  double sqrt2Sprime2 = sqrt(2.0 * Sprime2);

  return (rhoBar)*(csprime*Delta)*(csprime*Delta)*sqrt2Sprime2;

}

//------------------------------------------------------------------------

double VMSLESTerm::computeFilterWidth(double tetVol, SVec<double,3> &X, int nodeNum[4])

{

   if(delta) {
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
   else{
     return pow(tetVol, onethird);
   }

}

//------------------------------------------------------------------------

void VMSLESTerm::computeVelocityGradient(double dp1dxj[4][3],
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

//------------------------------------------------------------------------

void VMSLESTerm::computeStressTensor(double mu, double dudxj[3][3],
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

//------------------------------------------------------------------------

void VMSLESTerm::getPrimeValues(double *V[4], double *VBar[4],
				double uPrime[4][3], double ePrime[4])

{

  double u, v, w, en;

  // Get primes for node 0 //

  u            =  V[0][1];
  v            =  V[0][2];
  w            =  V[0][3];
  en           = (V[0][4] * oogamma1) / V[0][0];

  uPrime[0][0] = u  -  VBar[0][1];
  uPrime[0][1] = v  -  VBar[0][2];
  uPrime[0][2] = w  -  VBar[0][3];
  ePrime[0]    = en - (VBar[0][4] * oogamma1) / VBar[0][0];

  // End get primes for node 0 //


  // Get primes for node 1 //

  u            =  V[1][1];
  v            =  V[1][2];
  w            =  V[1][3];
  en           = (V[1][4] * oogamma1) / V[1][0];

  uPrime[1][0] = u  -  VBar[1][1];
  uPrime[1][1] = v  -  VBar[1][2];
  uPrime[1][2] = w  -  VBar[1][3];
  ePrime[1]    = en - (VBar[1][4] * oogamma1) / VBar[1][0];


  // End get primes for node 1 //


  // Get primes for node 2 //

  u            =  V[2][1];
  v            =  V[2][2];
  w            =  V[2][3];
  en           = (V[2][4] * oogamma1) / V[2][0];

  uPrime[2][0] = u  -  VBar[2][1];
  uPrime[2][1] = v  -  VBar[2][2];
  uPrime[2][2] = w  -  VBar[2][3];
  ePrime[2]    = en - (VBar[2][4]*oogamma1) / VBar[2][0];

  // End get primes for node 2 //


  // Get primes for node 3 //
  
  u            =  V[3][1];
  v            =  V[3][2];
  w            =  V[3][3];
  en           = (V[3][4] * oogamma1) / V[3][0];

  uPrime[3][0] = u  -  VBar[3][1];
  uPrime[3][1] = v  -  VBar[3][2];
  uPrime[3][2] = w  -  VBar[3][3];
  ePrime[3]    = en - (VBar[3][4]*oogamma1) / VBar[3][0];

  // End get primes for node 3 //
  

}

//------------------------------------------------------------------------

void VMSLESTerm::computeSmallEnergyGradient(double dp1dxj[4][3],
					    double eprime[4],
					    double dePrimedxj[3])

{

  dePrimedxj[0] = (dp1dxj[0][0]*eprime[0] +
		   dp1dxj[1][0]*eprime[1] +
		   dp1dxj[2][0]*eprime[2] +
		   dp1dxj[3][0]*eprime[3]);

  dePrimedxj[1] = (dp1dxj[0][1]*eprime[0]+
		   dp1dxj[1][1]*eprime[1]+
		   dp1dxj[2][1]*eprime[2]+
		   dp1dxj[3][1]*eprime[3]);

  dePrimedxj[2] = (dp1dxj[0][2]*eprime[0]+
		   dp1dxj[1][2]*eprime[1]+
		   dp1dxj[2][2]*eprime[2]+
		   dp1dxj[3][2]*eprime[3]);
}

//------------------------------------------------------------------------

double VMSLESTerm::computeMutOMu(double tetVol, 
			         double dp1dxj[4][3], 
			         double *VBar[4],
			         double *V[4], 
			         SVec<double,3> &X,
                                 int nodeNum[4])

{

  // Compute filter width //
  
  double Delta = computeFilterWidth(tetVol, X, nodeNum);

  // Compute small scale Reynolds stress tensor and small scale //
  // energy gradient //
  
  double uPrime[4][3];
  double ePrime[4];
  double duPrimedxj[3][3];
  double dePrimedxj[3];
  double tauPrimeij[3][3];
  getPrimeValues(V, VBar, uPrime, ePrime);
  computeVelocityGradient(dp1dxj, uPrime, duPrimedxj);
  double rhoBarCG = onefourth * (VBar[0][0] + VBar[1][0] +
				 VBar[2][0] + VBar[3][0]);
  double muT;
  
  // Applying the laminar-turbulent trip //

  if(trip){
    if((X[nodeNum[0]][0]>=x0 && X[nodeNum[0]][0]<=x1 && X[nodeNum[0]][1]>=y0 && X[nodeNum[0]][1]<=y1 &&
        X[nodeNum[0]][2]>=z0 && X[nodeNum[0]][2]<=z1) || (X[nodeNum[1]][0]>=x0 && X[nodeNum[1]][0]<=x1 &&
        X[nodeNum[1]][1]>=y0 && X[nodeNum[1]][1]<=y1 && X[nodeNum[1]][2]>=z0 && X[nodeNum[1]][2]<=z1) ||
       (X[nodeNum[2]][0]>=x0 && X[nodeNum[2]][0]<=x1 && X[nodeNum[2]][1]>=y0 && X[nodeNum[2]][1]<=y1 &&
        X[nodeNum[2]][2]>=z0 && X[nodeNum[2]][2]<=z1) || (X[nodeNum[3]][0]>=x0 && X[nodeNum[3]][0]<=x1 &&
        X[nodeNum[3]][1]>=y0 && X[nodeNum[3]][1]<=y1 && X[nodeNum[3]][2]>=z0 && X[nodeNum[3]][2]<=z1)){
          muT = computeEddyViscosity(rhoBarCG, Delta, duPrimedxj);
    }
    else{
       muT = 0.0; 
    }
  }
  else{
    muT = computeEddyViscosity(rhoBarCG, Delta, duPrimedxj);
  }

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);
  double mu = ooreynolds_mu * viscoFcn->compute_mu(Tcg);

  return (muT/mu);

}

//-----------------------------------------------------------------------
