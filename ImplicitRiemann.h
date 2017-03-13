#pragma once

#include "VarFcn.h"

class ImplicitRiemann {

public:


  static void computeGasGasJacobian(double Pstar, double gammai, double Pci, double Pi,double rhoi,
				    double gammaj, double Pcj, double Pj,double rhoj, 
				    double* dWidWi, double* dWidWj, double* dWjdWj, double* dWjdWi);

  static void computeTaitTaitJacobian(double Pstar, double alphai, double betai,double Pinfi,double Pi,double rhoi,
				      double alphaj, double betaj, double Pinfj, double Pj,double rhoj, double* dWidWi, 
				      double* dWidWj, double* dWjdWj, double* dWjdWi) ;

  static void computeGasTaitJacobian(double Pstar, double gammai, double Pci, double Pi,double rhoi,
				     double alphaj, double betaj, double Pinfj, double Pj,double rhoj,  
				     double* dWidWi, double* dWidWj, double* dWjdWj, double* dWjdWi,
                                     double cp, bool isBurnable);

  static void computeJwlJwlJacobian(VarFcn* vf, int fluidi, int fluidj,
				    double* Vi, double* Vj,
				    double* Wi, double* Wj,
				    double* jacii,double* jacij,
				    double* jacji,double* jacjj);

  static void computeGasJwlJacobian(VarFcn* vf, int fluidi, int fluidj,
				    double* Vi, double* Vj,
				    double* Wi, double* Wj,
				    double* jacii,double* jacij,
				    double* jacji,double* jacjj,double* dVdv);

  static void computeTaitJwlJacobian(VarFcn* vf, int fluidi, int fluidj,
				     double* Vi, double* Vj,
				     double* Wi, double* Wj,
				     double* jacii,double* jacij,
				     double* jacjj,double* jacji,
				     double* dVdv);
 private:

  static double computeJwlIntegral2(VarFcn* vf_, int fluidId, double omega,
				    double entropy, double rhostar, double rho);

  static void computeGasDerivRarefaction2x2(double gamma,double p,double pref,
					    double rhostar, double rho, double out[6], 
					    double c, double cstar);

  static void computeGasDerivShock2x2(double gamma,double p,
				      double pref,
				      double rhostar, double rho, 
				      double out[6]) ;

  static void computeTaitDerivRarefaction2x2(double alpha,double beta,double pinf,double p,
					     double rhostar, double rho, double out[6], 
					     double c, double cstar);

  static void computeTaitDerivShock2x2(double alpha,double beta,double pinf,double p,
				       double rhostar, double rho, 
				       double out[6]);

  static void computeJwlDerivRarefaction(VarFcn* vf_, int fluidId,double omega,
					 double entropy,
					 double rhostar, double rho, 
					 double out[6], double c, double cstar,
					 double* dVdv = NULL);

  static void computeJwlDerivShock(double omega,
				   double rhostar, double rho, 
				   double out[6], double c, double cstar,
				   double frho, double frhostar,
				   double fdrho,double fdrhostar,
				   double p);

  static void computeInterfaceQuantities2x2(double outi[6], double outj[6],
					    double* jacii,double* jacij,
					    double* jacji,double* jacjj);
};
