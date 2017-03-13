#include "ImplicitRiemann.h"
#include "ODEIntegrator.h"

#include <cmath>

#include <iostream>
#include <cstring>

#define max(a,b) ( (a > b) ? (a) : (b) )

// out[0] = K_rho_k
// out[1] = K_P_k
// out[2] = K_k*
// out[3] = dR/drho_k
// out[4] = dR/dPk
// out[5] = dR/dP*
void stiffenedGas_shock(double Pstar, double gamma, double Pck, double Pk,double rhok, double out[])
{

  //std::cout << "Pck = " << Pck << std::endl;
  double tmp = (Pstar+Pck) / (Pk+Pck);
  double k = (gamma-1.0)/(gamma+1.0);
  double a = 2.0/(rhok*(gamma+1.0));
  double b = k*(Pk+Pck);
  double tmp2 = (1.0-k*k) / ((tmp*k+1.0)*(tmp*k+1.0));

  double f = (Pstar-Pk)*sqrt(a/(Pstar+Pck+b));

  out[0] = -f/(2.0*rhok);
  out[1] = -(k*f*f+2.0*a*(Pstar-Pk))/(2.0*f*(Pstar+Pck+b));
  out[2] = (2.0*a*(Pstar-Pk)-f*f)/(2.0*f*(Pstar+Pck+b));

  double rhostar = rhok*(tmp+k)/(tmp*k+1.0);
  out[3] = rhostar/rhok;
  out[4] = -rhok*(Pstar+Pck)/( (Pk+Pck)*(Pk+Pck) ) * tmp2;
  out[5] = rhok/(Pk+Pck) * tmp2;

}

void stiffenedGas_rarefaction(double Pstar, double gamma, double Pck, double Pk,double rhok, double out[])
{

  //std::cout << "Pck = " << Pck << " " << "Pk = " << Pk << std::endl;
  double tmp = (Pstar+Pck) / (Pk+Pck);
  double a = sqrt(gamma*(Pk+Pck) / rhok);
  double f = 2.0*a/(gamma-1.0)*(pow( tmp, (gamma-1.0)/(2.0*gamma) ) - 1.0);

  out[0] = -f/(2.0*rhok);
  out[1] = (f-2.0*a)/(2.0*gamma*(Pk+Pck));
  out[2] = ((gamma-1.0)*f+2.0*a)/(2.0*gamma*(Pstar+Pck));

  double rhostar = rhok*pow(tmp, 1.0/gamma);
  out[3] = rhostar/rhok;
  out[4] = -rhostar/gamma/(Pk+Pck);
  out[5] = rhostar/gamma/(Pstar+Pck);

}

void tait_shock(double rhostar, double Pinf,double a, double alpha,
		double beta,double rhok, double out[])
{

  double q = pow(rhostar,beta)-pow(rhok,beta);
  double r = rhostar-rhok;
  double g = sqrt(alpha*q*r/(rhok*rhostar));
  out[0] = 0.5*g*(-1.0/rhok-beta*pow(rhok,beta-1.0)/q-1.0/r);
  out[1] = 0.5*g*(-1.0/rhostar+beta*pow(rhostar,beta-1.0)/q+1.0/r);  
}

void tait_rarefaction(double rhostar, double Pinf,double a, double alpha,
		      double beta,double rhok, double out[])
{

  double t = a*pow( rhostar/rhok, (beta-3.0)/2.0 );
  double da = a/rhok*(beta-1.0)*0.5;
  double gp = 2.0/(beta-1.0)*(pow( rhostar/rhok,(beta-1.0)*0.5) - 1.0);
  out[0] = gp*da-t*rhostar/(rhok*rhok);
  out[1] = t/rhok;
}

void computeInterfaceQuantities(const double outi[], const double outj[], 
				double jaci[], double jacj[] ) {

  const double D = outi[2]+outj[2];
  const double Dm = outj[2]-outi[2];
  const double Dinv = 1.0/D;
  const double DmDinv = Dm*Dinv;
  
  // First compute derivatives wrt i
  jaci[0] = -outi[0]*Dinv;
  jaci[1] = Dinv;
  jaci[2] = -outi[1]*Dinv;
  
  jaci[3] = -0.5*outi[0]-0.5*DmDinv*outi[0];
  jaci[4] = 0.5+0.5*DmDinv;
  jaci[5] = -0.5*outi[1]-0.5*DmDinv*outi[1];
  
  jaci[6] = outi[3]-outi[5]*outi[0]*Dinv;
  jaci[7] = Dinv*outi[5];
  jaci[8] = outi[4]-outi[5]*outi[1]*Dinv;
  
  jaci[9] = -outi[5]*outj[0]*Dinv;
  jaci[10] = -outi[5]*Dinv;	
  jaci[11] = -outi[5]*outj[1]*Dinv;
  
  // Compute derivatives wrt j
  jacj[0] = -outj[0]*Dinv;
  jacj[1] = -Dinv;
  jacj[2] = -outj[1]*Dinv;
  
  jacj[3] = 0.5*outj[0]-0.5*DmDinv*outj[0];
  jacj[4] = 0.5-0.5*DmDinv;
  jacj[5] = 0.5*outj[1]-0.5*DmDinv*outj[1];
  
  jacj[6] = outj[3]-outj[5]*outj[0]*Dinv;
  jacj[7] = -Dinv*outj[5];
  jacj[8] = outj[4]-outj[5]*outj[1]*Dinv;
  
  jacj[9] = -outj[5]*outi[0]*Dinv;
  jacj[10] = outj[5]*Dinv;	
  jacj[11] = -outj[5]*outi[1]*Dinv;

}

static void computeQuantities(double jaci[], double jacj[], double* dWidWi, double* dWidWj,
			      double* dWjdWj, double* dWjdWi) {
  
  dWidWi[0] = jaci[6]; dWidWi[1] = jaci[7]; dWidWi[2] = jaci[8];
  dWidWi[3] = jaci[3]; dWidWi[4] = jaci[4]; dWidWi[5] = jaci[5];
  dWidWi[6] = jaci[0]; dWidWi[7] = jaci[1]; dWidWi[8] = jaci[2];
  
  dWidWj[0] = jaci[9]; dWidWj[1] = jaci[10]; dWidWj[2] = jaci[11];
  dWidWj[3] = jacj[3]; dWidWj[4] = jacj[4]; dWidWj[5] = jacj[5];
  dWidWj[6] = jacj[0]; dWidWj[7] = jacj[1]; dWidWj[8] = jacj[2];
  
  
  dWjdWj[0] = jacj[6]; dWjdWj[1] = jacj[7]; dWjdWj[2] = jacj[8];
  dWjdWj[3] = jacj[3]; dWjdWj[4] = jacj[4]; dWjdWj[5] = jacj[5];
  dWjdWj[6] = jacj[0]; dWjdWj[7] = jacj[1]; dWjdWj[8] = jacj[2];
  
  dWjdWi[0] = jacj[9]; dWjdWi[1] = jacj[10]; dWjdWi[2] = jacj[11];
  dWjdWi[3] = jaci[3]; dWjdWi[4] = jaci[4]; dWjdWi[5] = jaci[5];
  dWjdWi[6] = jaci[0]; dWjdWi[7] = jaci[1]; dWjdWi[8] = jaci[2];
}

void ImplicitRiemann::computeGasGasJacobian(double Pstar, double gammai, double Pci, double Pi,double rhoi,
					    double gammaj, double Pcj, double Pj,double rhoj, double* dWidWi, 
					    double* dWidWj, double* dWjdWj, double* dWjdWi) {

  double outi[6], outj[6];
  double jaci[12], jacj[12];
  if (Pstar <= Pi) {
    
    stiffenedGas_rarefaction(Pstar, gammai, Pci, Pi, rhoi, outi);
  } else {
    
    stiffenedGas_shock(Pstar, gammai, Pci, Pi, rhoi, outi);
  }
  
  
  if (Pstar <= Pj) {
    
    stiffenedGas_rarefaction(Pstar, gammaj, Pcj, Pj, rhoj, outj);
    
  } else {
    
    stiffenedGas_shock(Pstar, gammaj, Pcj, Pj, rhoj, outj);
    
  }
  
  computeInterfaceQuantities(outi, outj, 
			     jaci, jacj );

  computeQuantities(jaci, jacj, dWidWi, dWidWj, dWjdWj, dWjdWi);
} 
								
void ImplicitRiemann::computeTaitTaitJacobian(double Pstar, double alphai, double betai,double Pinfi,double Pi,double rhoi,
					      double alphaj, double betaj, double Pinfj, double Pj,double rhoj, double* dWidWi, 
					      double* dWidWj, double* dWjdWj, double* dWjdWi) {

  double outi[6], outj[6];
  double jaci[12], jacj[12];

  double ai = sqrt(alphai*betai*pow(rhoi, betai-1.0));
  double aj = sqrt(alphaj*betaj*pow(rhoj, betaj-1.0));
  
  double rhostari = pow((Pstar-Pinfi)/alphai,1.0/betai);
  double rhostarj = pow((Pstar-Pinfj)/alphaj,1.0/betaj);

  double dpjdpi = alphai/alphaj*betai/betaj*pow(rhostari,betai-1.0)/pow(rhostarj,betaj-1.0);
   
  if (Pstar <= Pi) {
    
    tait_rarefaction(rhostari, Pinfi, ai,alphai,betai, rhoi, outi );
  } else {
    
    tait_shock(rhostari, Pinfi, ai,alphai,betai, rhoi, outi);
  }
  
  
  if (Pstar <= Pj) {
    
    tait_rarefaction(rhostarj, Pinfj, aj,alphaj,betaj, rhoj, outj);
    
  } else {
    
    tait_shock(rhostarj, Pinfj,aj,alphaj,betaj, rhoj, outj);
    
  }

  double si = outi[1]+outj[1]*dpjdpi;
  double sj = outj[1]+outi[1]/dpjdpi;
  
  memset(dWidWi,0,sizeof(double)*9);
  memset(dWidWj,0,sizeof(double)*9);
  memset(dWjdWi,0,sizeof(double)*9);
  memset(dWjdWj,0,sizeof(double)*9);

  dWidWi[0] = -outi[0]/si;
  dWidWi[1] = 1.0/si;

  dWidWj[0] = -outj[0]/si;
  dWidWj[1] = -1.0/si;

  dWjdWi[0] = -outi[0]/sj;
  dWjdWi[1] = 1.0/sj;

  dWjdWj[0] = -outj[0]/sj;
  dWjdWj[1] = -1.0/sj;

  dWidWi[3] = -0.5*(outi[0]+outi[1]*dWidWi[0])+0.5*(outj[1]*dWjdWi[0]);
  dWidWj[3] = 0.5*(outj[0]+outj[1]*dWjdWj[0])-0.5*(outi[1]*dWidWj[0]);
  dWjdWi[3] = -0.5*(outi[0]+outi[1]*dWidWi[0])+0.5*(outj[1]*dWjdWi[0]);
  dWjdWj[3] = 0.5*(outj[0]+outj[1]*dWjdWj[0])-0.5*(outi[1]*dWidWj[0]);

  dWidWi[4] = 0.5+0.5*(outj[1]*dWjdWi[1]-outi[1]*dWidWi[1]);
  dWidWj[4] = 0.5+0.5*(outj[1]*dWjdWj[1]-outi[1]*dWidWj[1]);
  dWjdWi[4] = 0.5+0.5*(outj[1]*dWjdWi[1]-outi[1]*dWidWi[1]);
  dWjdWj[4] = 0.5+0.5*(outj[1]*dWjdWj[1]-outi[1]*dWidWj[1]);
  
  dWidWi[8] = dWjdWj[8] = 1.0; 
} 

void ImplicitRiemann::computeGasTaitJacobian(double Pstar, double gammai, double Pci, double Pi,double rhoi,
					     double alphaj, double betaj, double Pinfj, double Pj,double rhoj,  
					     double* dWidWi, double* dWidWj, double* dWjdWj, double* dWjdWi,
                                             double cp, bool isBurnable) {


  double outi[6], outj[6];
  double jaci[12], jacj[12];

  //  double ai = sqrt(alphai*betai*pow(rhoi, betai-1.0));
  double aj = sqrt(alphaj*betaj*pow(rhoj, betaj-1.0));
  
  //double rhostari = pow((Pstar-Pinfi)/alphai,1.0/betai);
  double rhostarj = pow((Pstar-Pinfj)/alphaj,1.0/betaj);

  double drhodp = pow((Pstar-Pinfj)/alphaj,1.0/betaj-1.0)/(alphaj*betaj);
  if (Pstar <= Pi) {
    
    stiffenedGas_rarefaction(Pstar, gammai, Pci, Pi, rhoi, outi);
  } else {
    
    stiffenedGas_shock(Pstar, gammai, Pci, Pi, rhoi, outi);
  }
  
  if (rhostarj <= rhoj) {
    
    tait_rarefaction(rhostarj, Pinfj, aj,alphaj,betaj, rhoj, outj);
    
  } else {
    
    tait_shock(rhostarj, Pinfj,aj,alphaj,betaj, rhoj, outj);
    
  }

  for (int i = 0; i < 3; ++i) {
    outi[i] *= -1.0;
    outj[i] *= -1.0;
  }
  
  //std::cout << "Hello" << std::endl;

  memset(dWidWi,0,sizeof(double)*9);
  memset(dWidWj,0,sizeof(double)*9);
  memset(dWjdWi,0,sizeof(double)*9);
  memset(dWjdWj,0,sizeof(double)*9);

  //std::cout << "B " <<outi[2] << " " << outj[1] << " " << drhodp << std::endl;
  double s1 = (outi[2]+outj[1]*drhodp);
  dWidWi[8] = -outi[1]/s1;
  dWidWi[6] = -outi[0]/s1;
  dWidWi[7] = 1.0/s1;
  //dWjdWj[7] = 1.0/s1;
  //dWjdWi[7] = dWidWi[7];
  dWidWj[7] = -1.0/s1;
  dWidWj[6] = -outj[0]/s1;
  //std::cout << "A " << outi[2]*dWidWi[7] << " " << outj[1] << std::endl;
  dWjdWi[1] = (1.0-outi[2]*dWidWi[7])/outj[1];
  dWjdWj[1] = (-1.0-outi[2]*dWidWj[7])/outj[1];
  dWjdWi[2] = (-outi[1]-outi[2]*dWidWi[8])/outj[1];
  dWjdWi[0] = -(outi[0]+outi[2]*dWidWi[6])/outj[1];
  dWjdWj[0] = -(outi[2]*dWidWj[6]+outj[0])/outj[1];
  dWidWi[0] = outi[3]+outi[5]*dWidWi[6];
  dWidWj[0] = outi[5]*dWidWj[6];
  dWidWi[1] = outi[5]*dWidWi[7];
  dWidWj[1] = outi[5]*dWidWj[7];
  dWidWi[2] = outi[4]+outi[5]*dWidWi[8];

  dWidWi[3] = 0.5*outj[1]*dWjdWi[0]-0.5*outi[0]-0.5*outi[2]*dWidWi[6];
  dWidWi[4] = 0.5+0.5*outj[1]*dWjdWi[1]-0.5*outi[2]*dWidWi[7];
  dWidWi[5] = 0.5*outj[1]*dWjdWi[2]-0.5*outi[1]-0.5*outi[2]*dWidWi[8];
  dWjdWi[3] = dWidWi[3];
  dWjdWi[4] = dWidWi[4];
  dWjdWi[5] = dWidWi[5];

  dWidWj[3] = 0.5*outj[1]*dWjdWj[0]+0.5*outj[0]-0.5*outi[2]*dWidWj[6]; // PJSA
  dWidWj[4] = 0.5+0.5*outj[1]*dWjdWj[1]-0.5*outi[2]*dWidWj[7];
  dWidWj[5] = 0.0;
  dWjdWj[3] = dWidWj[3];
  dWjdWj[4] = dWidWj[4];
  dWjdWj[5] = dWidWj[5];

  //dWidWj[8] = dWidWj[5] = dWidWj[2] = 0.0;
  //dWjdWj[8] = 1.0;

/*double g1 = gammai-1.0;
  dWjdWj[6] = 1/(g1*rhostarj)*dWidWj[6]-(Pstar+gammai*Pci)/(g1*rhostarj*rhostarj)*dWjdWj[0];
  dWjdWj[7] = 1/(g1*rhostarj)*dWidWj[7]-(Pstar+gammai*Pci)/(g1*rhostarj*rhostarj)*dWjdWj[1];

  dWjdWi[6] = 1/(g1*rhostarj)*dWidWi[6]-(Pstar+gammai*Pci)/(g1*rhostarj*rhostarj)*dWjdWi[0];
  dWjdWi[7] = 1/(g1*rhostarj)*dWidWi[7]-(Pstar+gammai*Pci)/(g1*rhostarj*rhostarj)*dWjdWi[1];
  dWjdWi[8] = 1/(g1*rhostarj)*dWidWi[8]-(Pstar+gammai*Pci)/(g1*rhostarj*rhostarj)*dWjdWi[2];*/

  // PJSA
  if(isBurnable) {
    dWjdWj[6] = -0.5/cp*((1/rhostarj-1/rhoj)*(dWidWj[6]+alphaj*betaj*pow(rhoj,betaj-1))-(Pstar+Pj)*(1/(rhostarj*rhostarj)*dWjdWj[0]-1/(rhoj*rhoj)));
    dWjdWj[7] = -0.5/cp*((1/rhostarj-1/rhoj)*dWidWj[7]-(Pstar+Pj)/(rhostarj*rhostarj)*dWjdWj[1]);
    dWjdWj[8] = 1;

    dWjdWi[6] = -0.5/cp*((1/rhostarj-1/rhoj)*dWidWi[6]-(Pstar+Pj)/(rhostarj*rhostarj)*dWjdWi[0]);
    dWjdWi[7] = -0.5/cp*((1/rhostarj-1/rhoj)*dWidWi[7]-(Pstar+Pj)/(rhostarj*rhostarj)*dWjdWi[1]);
    dWjdWi[8] = -0.5/cp*((1/rhostarj-1/rhoj)*dWidWi[8]-(Pstar+Pj)/(rhostarj*rhostarj)*dWjdWi[2]);
  }
  else {
    dWjdWj[6] = 0.5/cp*((1/rhostarj+1/rhoj)*(dWidWj[6]-alphaj*betaj*pow(rhoj,betaj-1))-(Pstar-Pj)*(1/(rhostarj*rhostarj)*dWjdWj[0]+1/(rhoj*rhoj)));
    dWjdWj[7] = 0.5/cp*((1/rhostarj+1/rhoj)*dWidWj[7]-(Pstar-Pj)/(rhostarj*rhostarj)*dWjdWj[1]);
    dWjdWj[8] = 1;

    dWjdWi[6] = 0.5/cp*((1/rhostarj+1/rhoj)*dWidWi[6]-(Pstar-Pj)/(rhostarj*rhostarj)*dWjdWi[0]);
    dWjdWi[7] = 0.5/cp*((1/rhostarj+1/rhoj)*dWidWi[7]-(Pstar-Pj)/(rhostarj*rhostarj)*dWjdWi[1]);
    dWjdWi[8] = 0.5/cp*((1/rhostarj+1/rhoj)*dWidWi[8]-(Pstar-Pj)/(rhostarj*rhostarj)*dWjdWi[2]);
  }


  //std::cout << "Hello2" << std::endl;
}

struct JwlInfo {

  JwlInfo(double om, double ent,VarFcn* v, int f) : omega(om), vf_(v), fluidId(f),entropy(ent) {}
  double omega,entropy;
  VarFcn* vf_;
  int fluidId;
};

struct JwlKernel2 {

  double computeJwlIntegralKernel2(double rho, const JwlInfo& J) {
  
    double c = J.vf_->computeSoundSpeed(rho,J.entropy, J.fluidId);
    return 0.5*(J.omega+1.0)/(c*pow(rho,1.0-J.omega));
  }
};

double ImplicitRiemann::computeJwlIntegral2(VarFcn* vf_, int fluidId, double omega,
					    double entropy, double rhostar, double rho) {
  
  double res = 0.0;
  ODEIntegrator myIntegrator(rhostar, rho, 2000);
  
  JwlKernel2 J;
  myIntegrator.integrateMidpoint(J , res, &JwlKernel2::computeJwlIntegralKernel2,
				 JwlInfo(omega, entropy, vf_,fluidId));
  return res;
}

// out[0] = dV/dp
// out[1] = dV/drho
// out[2] = dV/drho*
// out[3] = dQ/dp;
// out[4] = dQ/drho
// out[5] = dQ/drho*
void ImplicitRiemann::computeGasDerivRarefaction2x2(double gamma,double p,double pref,
						    double rhostar, double rho, double out[6], 
						    double c, double) {

  double pp2 = pow(rhostar/rho, gamma);
  double ppref = (p+pref)*pp2;
  double cstar = sqrt(gamma*ppref/rhostar);

/*out[0] = 1.0/(gamma-1.0)*(c-cstar)/(p+pref);
  out[1] = 1.0/(gamma-1.0)*(gamma*cstar/rho-c/rho);
  out[2] = -cstar/rhostar;*/
  out[0] = -1.0/(gamma-1.0)*(c-cstar/pp2)/(p+pref); // PJSA
  out[1] = -1.0/(gamma-1.0)*(gamma*cstar/rho-c/rho); // PJSA
  out[2] = cstar/rhostar; // PJSA

  out[3] = pp2;
  out[4] = -pp2*c*c;
  out[5] = cstar*cstar;
}

void ImplicitRiemann::computeGasDerivShock2x2(double gamma,double p,
					      double pref,
					      double rhostar, double rho, 
					      double out[6]) {

  double a = (gamma)/(gamma-1.0);
  double b = 0.5*(1.0/rhostar+1.0/rho);
  double d = a/rhostar-b;
  double g = a/rho-b;
  double h = (1.0/rhostar-1.0/rho);

  double Q = g/d*(p+pref)-pref;
  double V = sqrt((p-Q)*h);

  out[3] = g/d;
  out[4] = -(p+pref)/(rho*rho)*((a-0.5)*d+0.5*g)/(d*d);
  out[5] = (p+pref)/(rhostar*rhostar)*(0.5*d+(a-0.5)*g)/(d*d);
  
  out[0] = h*(1.0-out[3])/(2.0*V);
  out[1] = 1.0/(2.0*V)*(-out[4]*h+(p-Q)/(rho*rho));
  out[2] = 1.0/(2.0*V)*(-out[5]*h-(p-Q)/(rhostar*rhostar));
}

void ImplicitRiemann::computeTaitDerivRarefaction2x2(double alpha,double beta,double pinf,double p,
						    double rhostar, double rho, double out[6], 
						    double c, double cstar) {
  
  //std::cout << "TAit rarefaction" << std::endl;
  double f = sqrt(alpha*beta);
  out[1] = -f*pow(rho,0.5*(beta-3.0) ); // PJSA
  out[2] = f*pow(rhostar,0.5*(beta-3.0) ); // PJSA
  out[0] = 0.0;

  out[3] = 0.0;
  out[4] = 0.0;
  out[5] = alpha*beta*pow(rhostar,beta-1.0);
}

void ImplicitRiemann::computeTaitDerivShock2x2(double alpha,double beta,double pinf,double p,
					      double rhostar, double rho, 
					      double out[6]) {

  //std::cout << "TAit shock" << std::endl;
  double pstar = pinf+alpha*pow(rhostar,beta);
  double f = 0.5/sqrt( max(0.0, -(pstar-p)*(1.0/rhostar-1.0/rho) ) );
  
  out[3] = 0.0;
  out[4] = 0.0;
  out[5] = alpha*beta*pow(rhostar,beta-1.0);

  out[0] = 0.0;
  out[1] = f*(-1.0/(rho*rho)*(pstar-p) + (1.0/rhostar-1.0/rho)*(alpha*beta*pow(rho,beta-1.0)) );
  out[2] = f*(-out[5]*(1.0/rhostar-1.0/rho) + (pstar-p)/(rhostar*rhostar) );
  
}

void ImplicitRiemann::computeJwlDerivRarefaction(VarFcn* vf_, int fluidId,double omega,
						 double entropy,
						 double rhostar, double rho, 
						 double out[6], double c, double cstar,double* dVdv) {

  double sign = -1.0;
  if (!dVdv) {
    double integral = computeJwlIntegral2(vf_, fluidId, omega, entropy, rhostar, rho);
  
    double pp = pow(rho, omega+1);
    out[0] = sign*1.0/pp*integral;
    out[1] = sign*(c/rho-c*c/pp*integral);
  } else {
    out[0] = sign*dVdv[0];
    out[1] = sign*dVdv[1];
  }

  out[2] = -sign*cstar/rhostar;

  double pp2 = pow(rhostar/rho, omega+1);
  out[3] = pp2;
  out[4] = -pp2*c*c;
  out[5] = cstar*cstar;
  //out[5] = entropy*(omega+1.0)*pow(rhostar,omega)+vf_->getA1(fluidId)*vf_->getR1r(fluidId)/(rhostar*rhostar)*exp(-vf_->getR1r(fluidId)/rhostar);
}

void ImplicitRiemann::computeJwlDerivShock(double omega,
					   double rhostar, double rho, 
					   double out[6], double c, double cstar,
					   double frho, double frhostar,
					   double fdrho,double fdrhostar,
					   double p) {

  double a = (omega+1.0)/omega;
  double b = 0.5*(1.0/rhostar+1.0/rho);
  double d = a/rhostar-b;
  double g = 1.0/rhostar-1.0/rho;

  double Q = ((a/rho-b)*p+1.0/omega*(frhostar/rhostar-frho/rho))/(a/rhostar-b);
  double V = sqrt((p-Q)*g);

  out[3] = (a/rho-b)/d;
//out[4] = -1.0/(rho*rho*d)*((a-0.5)*p-(fdrho*rho-frho)/omega+0.5*Q);
  out[4] = -1.0/(rho*rho*d)*((a-0.5)*p+(fdrho*rho-frho)/omega+0.5*Q); // PJSA
  out[5] = 1.0/(rhostar*rhostar*d)*(0.5*p+(fdrhostar*rhostar-frhostar)/omega+(a-0.5)*Q);
  
  out[0] = g*(1.0-out[3])/(2.0*V);
  out[1] = 1.0/(2.0*V)*(-out[4]*g+(p-Q)/(rho*rho));
  out[2] = 1.0/(2.0*V)*(-out[5]*g-(p-Q)/(rhostar*rhostar));
}

void ImplicitRiemann::computeInterfaceQuantities2x2(double outi[6], double outj[6],
						    double* jacii,double* jacij,
						    double* jacji,double* jacjj) {
  
  double delta = 1.0/(outj[5]*outi[2]+outi[5]*outj[2]);
  
  jacii[1] = delta*outj[5];
  jacij[1] = -delta*outj[5];
  jacii[0] = -delta*(outi[1]*outj[5]+outj[2]*outi[4]);
  jacij[0] = -delta*(outj[1]*outj[5]-outj[2]*outj[4]);
  jacii[2] = -delta*(outi[0]*outj[5]+outj[2]*outi[3]);
  jacij[2] = -delta*(outj[0]*outj[5]-outj[2]*outj[3]);
  
  jacji[1] = delta*outi[5];
  jacjj[1] = -delta*outi[5];
  jacji[0] = -delta*(outi[1]*outi[5]-outi[2]*outi[4]);
  jacjj[0] = -delta*(outj[1]*outi[5]+outi[2]*outj[4]);
  jacji[2] = -delta*(outi[0]*outi[5]-outi[2]*outi[3]);
  jacjj[2] = -delta*(outj[0]*outi[5]+outi[2]*outj[3]);

  jacii[4] = jacji[4] = 0.5*(1.0-outi[2]*jacii[1]+outj[2]*jacji[1]);
  jacij[4] = jacjj[4] = 0.5*(1.0-outi[2]*jacij[1]+outj[2]*jacjj[1]);
  jacii[3] = jacji[3] = 0.5*(outj[2]*jacji[0]-outi[2]*jacii[0]-outi[1]);
  jacij[3] = jacjj[3] = 0.5*(outj[2]*jacjj[0]-outi[2]*jacij[0]+outj[1]);
  jacii[5] = jacji[5] = 0.5*(outj[2]*jacji[2]-outi[2]*jacii[2]-outi[0]);
  jacij[5] = jacjj[5] = 0.5*(outj[2]*jacjj[2]+outj[0]-outi[2]*jacij[2]);

  jacii[7] = jacji[7] = 0.5*(outi[5]*jacii[1]+outj[5]*jacji[1]);
  jacij[7] = jacjj[7] = 0.5*(outi[5]*jacij[1]+outj[5]*jacjj[1]);
  jacii[6] = jacji[6] = 0.5*(outi[5]*jacii[0]+outi[4]+outj[5]*jacji[0]);
  jacij[6] = jacjj[6] = 0.5*(outi[5]*jacij[0]+outj[4]+outj[5]*jacjj[0]);
  jacii[8] = jacji[8] = 0.5*(outi[5]*jacii[2]+outi[3]+outj[5]*jacji[2]);
  jacij[8] = jacjj[8] = 0.5*(outi[5]*jacij[2]+outj[3]+outj[5]*jacjj[2]);


}

void ImplicitRiemann::computeJwlJwlJacobian(VarFcn* vf, int fluidi, int fluidj,
					    double* Vi, double* Vj,
					    double* Wi, double* Wj,
					    double* jacii,double* jacij,
					    double* jacji,double* jacjj) {

  double outi[6],outj[6];
  double ci = vf->computeSoundSpeed(Vi,fluidi);
  double cj = vf->computeSoundSpeed(Vj,fluidj);
  double cstari = vf->computeSoundSpeed(Wi,fluidi);
  double cstarj = vf->computeSoundSpeed(Wj,fluidj);  
  double omegai = vf->getOmega(fluidi);
  double omegaj = vf->getOmega(fluidj);
  double entropyi = vf->computeEntropy(Vi[0],Vi[4], fluidi);
  double entropyj = vf->computeEntropy(Vj[0],Vj[4], fluidj);
  
  if (Wi[4] <= Vi[4]) {

    computeJwlDerivRarefaction(vf,fluidi,omegai,entropyi,Wi[0],Vi[0], outi, ci,cstari);
    
  } else {

    double frhoi = vf->computeFrho(Vi,fluidi);
    double fdrhoi = vf->computeFrhop(Vi,fluidi);
    double frhostari = vf->computeFrho(Wi,fluidi);
    double fdrhostari = vf->computeFrhop(Wi,fluidi);
    
    computeJwlDerivShock(omegai,Wi[0],Vi[0], 
			 outi,ci,cstari,
			 frhoi,frhostari,
			 fdrhoi,fdrhostari,Vi[4]); 
  }
  
  if (Wi[4] <= Vj[4]) {

    computeJwlDerivRarefaction(vf,fluidj,omegaj,entropyj,Wj[0], Vj[0], outj, cj,cstarj);
    
  } else {
    
    double frhoj = vf->computeFrho(Vj,fluidj);
    double fdrhoj = vf->computeFrhop(Vj,fluidj);
    double frhostarj = vf->computeFrho(Wj,fluidj);
    double fdrhostarj = vf->computeFrhop(Wj,fluidj);
    
    computeJwlDerivShock(omegaj,Wj[0],Vj[0], 
			 outj,cj,cstarj,
			 frhoj,frhostarj,
			 fdrhoj,fdrhostarj,Vj[4]);
    
  }
  
  computeInterfaceQuantities2x2(outi, outj,
				jacii,jacij,
				jacji,jacjj);
}

void ImplicitRiemann::computeGasJwlJacobian(VarFcn* vf, int fluidi, int fluidj,
					    double* Vi, double* Vj,
					    double* Wi, double* Wj,
					    double* jacii,double* jacij,
					    double* jacjj,double* jacji,
					    double* dVdv) {

  double outi[6],outj[6];
  double ci = vf->computeSoundSpeed(Vi,fluidi);
  double cj = vf->computeSoundSpeed(Vj,fluidj);
  double cstari = vf->computeSoundSpeed(Wi,fluidi);
  double cstarj = vf->computeSoundSpeed(Wj,fluidj);  
  double gammai = vf->getGamma(fluidi);
  double omegaj = vf->getOmega(fluidj);
  double prefi = vf->getPressureConstant(fluidi);
  double entropyi = vf->computeEntropy(Vi[0],Vi[4], fluidi);
  double entropyj = vf->computeEntropy(Vj[0],Vj[4], fluidj);

  if (Wi[4] <= Vi[4]) {

    computeGasDerivRarefaction2x2(gammai, Vi[4],prefi, Wi[0], Vi[0], outi, ci, cstari);
    
  } else {

    computeGasDerivShock2x2(gammai,Vi[4],prefi,Wi[0],Vi[0],outi);
  }
  
  if (Wi[4] <= Vj[4]) {

    computeJwlDerivRarefaction(vf,fluidj,omegaj,entropyj,Wj[0], Vj[0], outj, cj,cstarj,dVdv);
    
  } else {
    
    double frhoj = vf->computeFrho(Vj,fluidj);
    double fdrhoj = vf->computeFrhop(Vj,fluidj);
    double frhostarj = vf->computeFrho(Wj,fluidj);
    double fdrhostarj = vf->computeFrhop(Wj,fluidj);
    
    computeJwlDerivShock(omegaj,Wj[0],Vj[0], 
			 outj,cj,cstarj,
			 frhoj,frhostarj,
			 fdrhoj,fdrhostarj,Vj[4]);
    
  }
  
  computeInterfaceQuantities2x2(outi, outj,
				jacii,jacij,
				jacji,jacjj);
}

void ImplicitRiemann::computeTaitJwlJacobian(VarFcn* vf, int fluidi, int fluidj,
					     double* Vi, double* Vj,
					     double* Wi, double* Wj,
					     double* jacii,double* jacij,
					     double* jacjj,double* jacji,
					     double* dVdv) {

  double outi[6],outj[6];
  double ci = vf->computeSoundSpeed(Vi,fluidi);
  double cj = vf->computeSoundSpeed(Vj,fluidj);
  double cstari = vf->computeSoundSpeed(Wi,fluidi);
  double cstarj = vf->computeSoundSpeed(Wj,fluidj);
  double alphai = vf->getAlphaWater(fluidi);
  double betai = vf->getBetaWater(fluidi);
  double pinfi = vf->getPrefWater(fluidi);
  double cp = vf->specificHeatCstPressure(fluidi);
  bool isBurnable = vf->isBurnable(fluidi);
  
  double omegaj = vf->getOmega(fluidj);
  double entropyj = vf->computeEntropy(Vj[0],Vj[4], fluidj);
  double pi = vf->getPressure(Vi,fluidi);

  if (Wj[4] <= pi) {

    computeTaitDerivRarefaction2x2(alphai,betai,pinfi,pi, Wi[0], Vi[0], outi, ci, cstari);
    
  } else {

    computeTaitDerivShock2x2(alphai,betai,pinfi,pi,Wi[0],Vi[0],outi);
  }
  
  if (Wj[4] <= Vj[4]) {

    computeJwlDerivRarefaction(vf,fluidj,omegaj,entropyj,Wj[0], Vj[0], outj, cj,cstarj,dVdv);
    
  } else {
    
    double frhoj = vf->computeFrho(Vj,fluidj);
    double fdrhoj = vf->computeFrhop(Vj,fluidj);
    double frhostarj = vf->computeFrho(Wj,fluidj);
    double fdrhostarj = vf->computeFrhop(Wj,fluidj);
    
    computeJwlDerivShock(omegaj,Wj[0],Vj[0], 
			 outj,cj,cstarj,
			 frhoj,frhostarj,
			 fdrhoj,fdrhostarj,Vj[4]);
    
  }
  
  computeInterfaceQuantities2x2(outi, outj,
				jacii,jacij,
				jacji,jacjj);

  // PJSA
  double pstari = pinfi + alphai*pow(Wi[0],betai);
  const double &rhostari = Wi[0];
  const double &rhoi = Vi[0];
  if(isBurnable) {
    jacii[6] = -0.5/cp*((1/rhostari-1/rhoi)*(jacji[6]+alphai*betai*pow(rhoi,betai-1))-(pstari+pi)*(1/(rhostari*rhostari)*jacii[0]-1/(rhoi*rhoi)));
    jacii[7] = -0.5/cp*((1/rhostari-1/rhoi)*jacji[7]-(pstari+pi)/(rhostari*rhostari)*jacii[1]);
    jacii[8] = 1;
  }
  else {
    jacii[6] = 0.5/cp*((1/rhostari+1/rhoi)*(jacji[6]-alphai*betai*pow(rhoi,betai-1))-(pstari-pi)*(1/(rhostari*rhostari)*jacii[0]+1/(rhoi*rhoi)));
    jacii[7] = 0.5/cp*((1/rhostari+1/rhoi)*jacji[7]-(pstari-pi)/(rhostari*rhostari)*jacii[1]);
    jacii[8] = 1;
  }

  if(isBurnable) {
    jacij[6] = -0.5/cp*((1/rhostari-1/rhoi)*jacjj[6]-(pstari+pi)/(rhostari*rhostari)*jacij[0]);
    jacij[7] = -0.5/cp*((1/rhostari-1/rhoi)*jacjj[7]-(pstari+pi)/(rhostari*rhostari)*jacij[1]);
    jacij[8] = -0.5/cp*((1/rhostari-1/rhoi)*jacjj[8]-(pstari+pi)/(rhostari*rhostari)*jacij[2]);
  }
  else {
    jacij[6] = 0.5/cp*((1/rhostari+1/rhoi)*jacjj[6]-(pstari-pi)/(rhostari*rhostari)*jacij[0]);
    jacij[7] = 0.5/cp*((1/rhostari+1/rhoi)*jacjj[7]-(pstari-pi)/(rhostari*rhostari)*jacij[1]);
    jacij[8] = 0.5/cp*((1/rhostari+1/rhoi)*jacjj[8]-(pstari-pi)/(rhostari*rhostari)*jacij[2]);
  }
}
