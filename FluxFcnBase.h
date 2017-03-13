#ifndef _FLUX_FCN_BASE_H_
#define _FLUX_FCN_BASE_H_

#include <VarFcnBase.h>

//----------------------------------------------------------------------------------------

class FluxFcnBase {

public:

  enum Type {CONSERVATIVE = 0, PRIMITIVE = 1} typeJac;

protected:
  VarFcnBase *vf;

  double* hhcoeffptr;

public:
  FluxFcnBase(VarFcnBase *varFcn, Type tp);
  virtual ~FluxFcnBase() { vf = 0; }

  virtual void compute(double, double, double *, double, double *, double *, double *, bool) {}
  virtual void computeJacobian(double, double, double *, double, double *, double *, double *, bool) {}
  virtual void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool) {}
  virtual void computeJacobianRight(double length, double irey, double *normal, double normalVel,
				    double *VL, double *VR, double *jacR, bool useLimiter) {
    std::cout << "Error: computeJacobianRight() not implemented!" << std::endl;
    exit(-1);
  }

  virtual void computeJacobianFarfield(double length, double irey, double *normal, double normalVel, double * VL, double * Ub, double *jac, bool useLimiter = true) {
    
    std::cout << "Error: computeJacobianFarfield() not implemented!" << std::endl;
    exit(-1);
  }
  

  void setHHCoeffPointer(double* hh) { hhcoeffptr = hh; }
 
  VarFcnBase* getVarFcnBase() const { return vf; }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv, 
    double *vl, double *dvl, double *vr, double *dvr, 
    double dmach, double *f, double *df
    , bool useLimiter = false
  ) 
  {
    std::cout << "\n !!! FluxFcnBase::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  virtual void compute_dFluxdNormal_dFluxdNormalVel_dFluxdVL_dFluxdVR(double *n, double nv, double *vl, double *vr,
                                                                      double dmach, double *f, double dfdn[7][3],
                                                                      double *fdnv, double dfdvl[7][7], double dfdvr[7][7])
  {
    std::cout << "\n !!! FluxFcnBase::compute_dFluxdNormal_dFluxdNormalVel_dFluxdVL_dFluxdVR is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv, 
    double *v, double *ub, double *dub, double *f, double *df
  ) 
  {
    std::cout << "\n !!! FluxFcnBase::computeDerivative (11 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  virtual void computeDerivativeOperators(double *n, double nv, double *v, 
                                          double *ub, double dfdn[7][3], double dfdv[7][1], double dfdub[7][7])
  {
    std::cout << "\n !!! FluxFcnBase::computeDerivativeOperators is not implemented !!!\n\n";
    exit(1);
  }

};


//----------------------------------------------------------------------------------------

inline
FluxFcnBase::FluxFcnBase(VarFcnBase *varFcn,Type tp) : vf(varFcn) {
 
  typeJac = tp;
  
}
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

template<int dim>
class FluxFcnFD : public FluxFcnBase {

public:
  FluxFcnFD(VarFcnBase *vf,Type tp) : FluxFcnBase(vf,tp) {}
  ~FluxFcnFD() { vf = 0; }

  virtual void compute(double, double, double *, double, double *, double *, double *, bool){} 
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);
  void computeJacobianRight(double length, double irey, double *normal, double normalVel,
			    double *VL, double *VR, double *jacR, bool useLimiter);
};

//------------------------------------------------------------------------------

template<int dim>
inline
void FluxFcnFD<dim>::computeJacobian(double length, double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *jacL, bool useLimiter)
{

  const double eps0 = 1.e-6;

  double Veps[dim], flux[dim], fluxeps[dim], dfdVL[dim*dim];
 
  

  compute(length, irey, normal, normalVel, VL, VR, flux, useLimiter); 

  int k;
  for (k=0; k<dim; ++k)
    Veps[k] = VL[k];

  for (k=0; k<dim; ++k) {

    double eps;
    if (fabs(Veps[k]) < 1.e-10)
      eps = eps0;
    else
      eps = eps0 * Veps[k];

    double inveps = 1.0 / eps;

    Veps[k] += eps;

    if (k != 0)
      Veps[k-1] = VL[k-1];

    compute(length, irey, normal, normalVel, Veps, VR, fluxeps, useLimiter);

    for (int j=0; j<dim; ++j) 
      dfdVL[dim*j + k] = (fluxeps[j]-flux[j]) * inveps;

  }

  if (typeJac == CONSERVATIVE)
    vf->postMultiplyBydVdU(VL, dfdVL, jacL);
  else{
    for (k=0; k<dim*dim; ++k) 
      jacL[k] = dfdVL[k];
  }

}

template<int dim>
inline
void FluxFcnFD<dim>::computeJacobianRight(double length, double irey, double *normal, double normalVel,
					  double *VL, double *VR, double *jacR, bool useLimiter)
{

  const double eps0 = 1.e-6;

  double Veps[dim], flux[dim], fluxeps[dim], dfdVR[dim*dim];
 
  compute(length, irey, normal, normalVel, VL, VR, flux, useLimiter); 

  int k;
  for (k=0; k<dim; ++k)
    Veps[k] = VR[k];

  for (k=0; k<dim; ++k) {

    double eps;
    if (fabs(Veps[k]) < 1.e-10)
      eps = eps0;
    else
      eps = eps0 * Veps[k];

    double inveps = 1.0 / eps;

    Veps[k] += eps;

    if (k != 0)
      Veps[k-1] = VR[k-1];

    compute(length, irey, normal, normalVel, VL, Veps, fluxeps, useLimiter);

    for (int j=0; j<dim; ++j) 
      dfdVR[dim*j + k] = (fluxeps[j]-flux[j]) * inveps;

  }

  if (typeJac == CONSERVATIVE)
    vf->postMultiplyBydVdU(VR, dfdVR, jacR);
  else{
    for (k=0; k<dim*dim; ++k) 
      jacR[k] = dfdVR[k];
  }

}

//------------------------------------------------------------------------------

template<int dim>
inline
void FluxFcnFD<dim>::computeJacobians(double length, double irey, double *normal, double normalVel, double *VL,
                                      double *VR, double *jacL, double *jacR, bool useLimiter)
{

  double n[3] = {normal[0], normal[1], normal[2]};

  computeJacobian(length, irey, n, normalVel, VL, VR, jacL, useLimiter);

  n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2];
  computeJacobian(length, irey, n, -normalVel, VR, VL, jacR, useLimiter);
  for (int k=0; k<dim*dim; ++k) jacR[k] = -jacR[k];

}

//------------------------------------------------------------------------------

#endif
