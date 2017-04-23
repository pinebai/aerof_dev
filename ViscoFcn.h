#ifndef _VISCO_FCN_H_
#define _VISCO_FCN_H_

#include <IoData.h>

#include <cmath>


const double twothird = 2.0/3.0;
//------------------------------------------------------------------------------
//CHANGES_FOR_WATER
//	new behavior of Lame coefficients (lambda and mu)
//------------------------------------------------------------------------------

class ViscoFcn {

protected:
  double ooreynolds_mu;
// Included (MB)
  double dReMach;
  double dRe_muMach;

public:

  ViscoFcn(IoData &iod) { 
  ooreynolds_mu = 1.0/iod.ref.reynolds_mu; 
// Included (MB)
  dRe_muMach = iod.ref.dRe_mudMach;
  }
  virtual ~ViscoFcn() {}

  virtual double compute_mu(double) = 0;
  virtual double compute_lambda(double, double) = 0;

// Included (MB)
  virtual double compute_muDerivative(double, double, double) = 0;
  virtual void compute_muDerivativeOperators(double, double&, double &) = 0;
  virtual double compute_lambdaDerivative(double, double, double) = 0;
  virtual void compute_lambdaDerivativeOperators(double&, double&) = 0;
  virtual void rstVar(IoData &) = 0;

};

//------------------------------------------------------------------------------

class ConstantViscoFcn : public ViscoFcn {

  double bulkViscosity;
public:

  ConstantViscoFcn(IoData &iod):ViscoFcn(iod) { bulkViscosity = iod.eqs.viscosityModel.bulkViscosity;}
  ~ConstantViscoFcn() {}

  double compute_mu(double T) {
    //std::cout<<"compute_mu Constant"<<std::endl;
    return 1.0;}
  double compute_lambda(double T, double mu) { return bulkViscosity - twothird*mu; }

// Included (MB)
  double compute_muDerivative(double T, double dT, double dMach) { return 0.0; }
  void compute_muDerivativeOperators(double T, double &dmudT, double &dmudMach) { dmudT = 0.0; dmudMach = 0.0; }
  double compute_lambdaDerivative(double mu, double dmu, double dMach) { 
    return -twothird*dmu;
  }
  void compute_lambdaDerivativeOperators(double &dlambdadmu, double &dlambdadMach) { //YC
    dlambdadmu = -twothird;  dlambdadMach = 0.0;
  }

  void rstVar(IoData &iod) {  
    ooreynolds_mu = 1.0/iod.ref.reynolds_mu; 
    dRe_muMach = iod.ref.dRe_mudMach;
  }

};

//------------------------------------------------------------------------------

class SutherlandViscoFcn : public ViscoFcn {

  double alpha;
  double Ts;

// Included (MB)
  double dalpha;

public:

  SutherlandViscoFcn(IoData &iod) : ViscoFcn(iod)
  { 
    double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
    alpha = gam*(gam - 1.0) * iod.ref.mach*iod.ref.mach; 
    Ts = iod.eqs.viscosityModel.sutherlandReferenceTemperature / iod.ref.temperature;

// Included (MB)
    dalpha = 2.0*gam*(gam - 1.0) * iod.ref.mach;
  }
  ~SutherlandViscoFcn() {}

  double compute_mu(double Tadim)
  {
    //std::cout<<"compute_mu Sutherland"<<std::endl;
    double T = alpha * Tadim;
    return T * sqrt(T) * (1.0 + Ts) / (T + Ts); 
  }
  
  double compute_lambda(double Tadim, double mu) { return -twothird*mu; }

// Included (MB)
  double compute_muDerivative(double Tadim, double dTadim, double dMach)
  {
    double dAlpha = dalpha*dMach;
    double T = alpha * Tadim;
    double dT = dAlpha * Tadim + alpha * dTadim;
    return ( ( 1.5*(1.0 + Ts)*sqrt(T)*dT*(T + Ts) - T*sqrt(T)*(1.0 + Ts)*dT ) / ( (T + Ts)*(T + Ts) ) );
  }

  void compute_muDerivativeOperators(double Tadim, double &dmudTadim, double &dmudMach) //YC
  {
    double T = alpha * Tadim;
    double coef = 1.0/( (T + Ts)*(T + Ts) ), coef1 = 1.5*(1.0 + Ts)*sqrt(T)*(T+Ts), coef2 = - T*sqrt(T)*(1.0 + Ts);
    dmudTadim = coef*(coef1+coef2)*alpha;
    dmudMach = coef*(coef1+coef2)*dalpha*Tadim;
  }

  double compute_lambdaDerivative(double mu, double dmu, double dMach) { 
    return -twothird*dmu;
  }

  void compute_lambdaDerivativeOperators(double &dlambdadmu, double &dlambdadMach) {//YC
    dlambdadmu = -twothird;  dlambdadMach = 0.0;
  }

  void rstVar(IoData &iod)
  {
    ooreynolds_mu = 1.0/iod.ref.reynolds_mu; 
    dRe_muMach = iod.ref.dRe_mudMach;

    double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
    alpha = gam*(gam - 1.0) * iod.ref.mach*iod.ref.mach; 
    Ts = iod.eqs.viscosityModel.sutherlandReferenceTemperature / iod.ref.temperature;
    dalpha = 2.0*gam*(gam - 1.0) * iod.ref.mach;
  }

};

//------------------------------------------------------------------------------

class PrandtlViscoFcn : public ViscoFcn {

  double alpha;

// Included (MB)
  double dalpha;

public:

  PrandtlViscoFcn(IoData &iod) : ViscoFcn(iod)
  { 
    double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
    alpha = gam*(gam - 1.0) * iod.ref.mach*iod.ref.mach;

// Included (MB)
    dalpha = 2.0*gam*(gam - 1.0) * iod.ref.mach;
  }
  ~PrandtlViscoFcn() {}

  double compute_mu(double Tadim) {
    //std::cout<<"compute_mu Prandtl"<<std::endl
    return alpha * Tadim; }
  double compute_lambda(double Tadim, double mu) { return -twothird*mu;}

// Included (MB)
  double compute_muDerivative(double Tadim, double dTadim, double dMach)
  {
    double dAlpha = dalpha*dMach;
    return (dAlpha * Tadim + alpha * dTadim);
  }
  void compute_muDerivativeOperators(double Tadim, double &dmudTadim, double &dmudMach) //YC
  {
    dmudTadim = alpha;
    dmudMach = dalpha*Tadim;
  }
  double compute_lambdaDerivative(double mu, double dmu, double dMach) { 
    return -twothird*dmu;
  }

  void compute_lambdaDerivativeOperators(double &dlambdadmu, double &dlambdadMach) {//YC
    dlambdadmu = -twothird;  dlambdadMach = 0.0;
  }

  void rstVar(IoData &iod)
  {
    ooreynolds_mu = 1.0/iod.ref.reynolds_mu; 
    dRe_muMach = iod.ref.dRe_mudMach;

    double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
    alpha = gam*(gam - 1.0) * iod.ref.mach*iod.ref.mach;
    dalpha = 2.0*gam*(gam - 1.0) * iod.ref.mach;
  }

};

//------------------------------------------------------------------------------

#endif
