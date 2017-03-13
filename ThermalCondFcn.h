#ifndef _THERMAL_COND_FCN_H_
#define _THERMAL_COND_FCN_H_

#include "IoData.h"
#include "ViscoFcn.h"
#include "VarFcn.h"



// Pure virtual base class
class ThermalCondFcn {

public:

  ThermalCondFcn() {}
  virtual ~ThermalCondFcn() {}

  virtual double compute(double) = 0;

// Included (MB)
  virtual double computeDerivative(double, double, double) = 0;
  virtual void computeDerivativeOperators(double, double &, double &) = 0;
  virtual void rstVar(IoData&) = 0;

};
//------------------------------------------------------------------------------

class ConstantThermalCondFcn : public ThermalCondFcn {

  // non-dimensional value of thermal conductivity coefficient
  // (non-dimensionalized by reference viscosity x reference Cv, already done in IoDataCore.C)
  double thermal_conductivity_coefficient;

public:

  ConstantThermalCondFcn(IoData &iod){
    thermal_conductivity_coefficient = iod.eqs.thermalCondModel.conductivity;
  }
  ~ConstantThermalCondFcn() {}

  double compute(double Tadim) { return thermal_conductivity_coefficient; }
  double computeDerivative(double Tadim, double dTadim, double dMach) { return 0.0; }
  void computeDerivativeOperators(double Tadim, double &dkappadTaim, double &dkappadMach) { dkappadTaim=0.0; dkappadMach=0.0;}
  void rstVar(IoData &iod) { thermal_conductivity_coefficient = iod.eqs.thermalCondModel.conductivity; }

};

//------------------------------------------------------------------------------
// assumption: constant Prandtl number (as well as constant cp)
class ConstantPrandtlThermalCondFcn : public ThermalCondFcn {

  double ooPrandtl;
  double ooTurbPrandtl;

  ViscoFcn *viscoFcn;
  VarFcn *varFcn;

public:

  ConstantPrandtlThermalCondFcn(IoData &iod, ViscoFcn *visf, VarFcn *vfn) :
    varFcn(vfn)
  {
    viscoFcn = visf;

    ooPrandtl = 1.0 / iod.eqs.thermalCondModel.prandtl;
    ooTurbPrandtl = 1.0 / iod.eqs.tc.prandtlTurbulent;

  }
  ~ConstantPrandtlThermalCondFcn() { viscoFcn = 0; varFcn = 0;}

  double compute(double Tadim)
  {
    return ooPrandtl * varFcn->specificHeatCstPressure() * viscoFcn->compute_mu(Tadim);
  }

  void computeDerivativeOperators(double Tadim, double &dkappadTadim, double &dkappadMach)
  {
    double dmudTadim(0), dmudMach(0);
    viscoFcn->compute_muDerivativeOperators(Tadim, dmudTadim, dmudMach);
    dkappadTadim = ooPrandtl * varFcn->specificHeatCstPressure() * dmudTadim;
    dkappadMach = ooPrandtl * varFcn->specificHeatCstPressure() * dmudMach;
  }

  double turbulentConductivity(double mut)
  { 
    return varFcn->specificHeatCstPressure() * mut * ooTurbPrandtl;
  }

// Included (MB)
  double computeDerivative(double Tadim, double dTadim, double dMach)
  {
    return ooPrandtl * varFcn->specificHeatCstPressure() * viscoFcn->compute_muDerivative(Tadim, dTadim, dMach);
  }

  double turbulentConductivityDerivative(double dmut)
  {
    return varFcn->specificHeatCstPressure() * dmut * ooTurbPrandtl;
  }
  void rstVar(IoData &iod)
  {
    viscoFcn->rstVar(iod);
    varFcn->rstVar(iod);
  }

};

//---------------------------------------------------------------------------------

#endif
