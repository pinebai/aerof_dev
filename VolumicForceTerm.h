#ifndef _VOLUMIC_FORCE_TERM_H_
#define _VOLUMIC_FORCE_TERM_H_


#include<IoData.h>
#include <cmath>

//------------------------------------------------------------------------------

class VolumicForceTerm {

  double volforce[3];

public:

  VolumicForceTerm(IoData &iod){
    volforce[0] = iod.eqs.gravity_x;
    volforce[1] = iod.eqs.gravity_y;
    volforce[2] = iod.eqs.gravity_z;
  }
    
  ~VolumicForceTerm() {}

  void computeVolumeTerm(double ctrlVol, double *V, double *flux){
    double rhovol = V[0]*ctrlVol;

    flux[0] = 0.0;
    flux[1] = -rhovol*volforce[0];
    flux[2] = -rhovol*volforce[1];
    flux[3] = -rhovol*volforce[2];
    flux[4] = -rhovol*(V[1]*volforce[0]+V[2]*volforce[1]+V[3]*volforce[2]);
  }

  void computeJacobianVolumeTerm(int dim, double ctrlVol, 
                                 double *V, double *jac){
    for (int i=0; i<dim*dim; i++)
      jac[i] = 0.0;
  
    jac[dim]     = volforce[0];
    jac[2*dim]   = volforce[1];
    jac[3*dim]   = volforce[2];
    jac[4*dim+1] = jac[dim];
    jac[4*dim+2] = jac[2*dim];
    jac[4*dim+3] = jac[3*dim];
  }

// Included (MB)
  void computeDerivativeOfVolumeTerm(double ctrlVol, double dCtrlVol, double *V, double *dV, double *dFlux){
    double rhovol  = V[0]*ctrlVol;
    double drhovol = dV[0]*ctrlVol + V[0]*dCtrlVol;

    dFlux[0] = 0.0;
    dFlux[1] = -drhovol*volforce[0];
    dFlux[2] = -drhovol*volforce[1];
    dFlux[3] = -drhovol*volforce[2];
    dFlux[4] = -drhovol*(V[1]*volforce[0]+V[2]*volforce[1]+V[3]*volforce[2]) 
             - rhovol*(dV[1]*volforce[0]+dV[2]*volforce[1]+dV[3]*volforce[2]);
  }

};
  
//------------------------------------------------------------------------------

#endif
