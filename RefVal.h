#ifndef _REF_VAL_H_
#define _REF_VAL_H_

//------------------------------------------------------------------------------
/*
	for gas!
  rho_bar = rho/rho_ref, rho_ref = rho_inf
  u_bar = u/U_ref, U_ref = U_inf = Mach * sqrt(gam * p_inf / rho_inf)
  p_bar = p/p_ref, p_ref = rho_inf * U_inf * U_inf
  t_bar = t/t_ref, t_ref = L_ref/U_ref
	
	for liquid, sound speed is different!
  rho_bar = rho/rho_ref, rho_ref = rho_inf
  u_bar = u/U_ref, U_ref = U_inf = Mach * sqrt(alpha * beta * rho_inf^(beta-1)))
  p_bar = p/p_ref, p_ref = rho_inf * U_inf * U_inf
  t_bar = t/t_ref, t_ref = L_ref/U_ref

*/

class IoData;
class RefVal {

public:

  enum Mode {NON_DIMENSIONAL = 0, DIMENSIONAL = 1} mode;

  double length;
  double density;
  double velocity;
  double pressure;
  double temperature;
  double viscosity_mu;
  double nutilde;
  double kenergy;
  double epsilon;
  double time;
  double force;
  double energy;
  double power;
  double entropy;

  double tlength;
  double tvelocity;
  double tforce;
  double tpower;

  double mach;
  double moment;
  double cf2force;
  double cm2moment;
  //double surface;
  
// Included (MB)
  double dvelocitydMach;
  double dtimedMach;

public:

  RefVal();
  //RefVal(IoData &ioData);
  ~RefVal() {}

// Included (MB)
  void rstVar(IoData &);

};

//------------------------------------------------------------------------------

#endif
