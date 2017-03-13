#include <FluxFcnDescSG.h>
#include <AutoDiff/roeturkeljac5.h>
#include <LinkF77.h>

#include <cstdlib>
#include <cstdio>
#include <iostream>


//------------------------------------------------------------------------------
// fortran routines located in f77src folder
/*
@ARTICLE{roe-81,
  author = "Roe, P. L.",
  title = "Approximate {R}iemann Solvers, Parameter Vectors
           and Difference Schemes",
  journal = j. comp. phys.,
  year = 1981,
  volume = 43,
  pages = "357--372",
} 
*/

extern "C" {

  void F77NAME(roeflux1)(const double&, const double&, const double&, double*, 
			 const double&, double*, double*, double*,
			 const double&, const double&, const double&,
                         const double&, const int&);
  void F77NAME(hlleflux1)(const double&, const double&, const double&, double*, 
			 const double&, double*, double*, double*,
			 const double&, const double&, const double&,
                         const double&, const double&, const double&, const int&);
  void F77NAME(roeflux5)(const int&, const double&, const double&, const double&, double*,
                         const double&, double*, double*, double*, double*, double*, 
                         const double&, const double&, const double&, const double&, 
                         const double&, const double&, const int&);
  void F77NAME(roeflux6)(const int&, const double&, const double&, const double&, double*,
                         const double&, double*, double*, double*, double*, double*);
  void F77NAME(roeflux7)(const int&, const double&, const double&, double*,
                         const double&, double*, double*, double*);
  void F77NAME(roejac2)(const int&, const double&, const double&, 
			const double&, double*,
			const double&, double*, double*, double*, double*,
			const double&, const double&, const double&,
			const double&, const int&);
  void F77NAME(roejac5)(const int&, const double&, const double&, const double&, double*, 
			const double&, double*, double*, double*);
  void F77NAME(roejac6)(const int&, const double&, const double&, const double&, double*, 
			const double&, double*, double*, double*,const int&,
                        const double&, const double&, const double&, const double&, 
                        const double&, const double&, const int&);
  void F77NAME(hllejac)(const int&, const double&, const double&, const double&, double*, 
			const double&, double*, double*, double*,const int&,
                        const double&, const double&, const double&, const double&,
                        const double&, const double&, const int&);
  void F77NAME(boundflux5)(const int&, const double&, double*, const double&, 
			   double*, double*, double*);
  void F77NAME(boundjac2)(const int&, const double&, double*, const double&, 
			  double*, double*, double*);
  void F77NAME(genbcfluxgas)(const int&, const double &, const double&, double*, 
			     const double&, double*, double*, double*);
  void F77NAME(genbcfluxgas_hh)(const int&, const double &, const double&, double*, 
			     const double&, double*, double*, double*, double*, double&, double&, const double&);
  void F77NAME(hlleflux)(const int&, const double&, const double&, const double&, double*,
                         const double&, double*, double*, double*, double*, double*, 
                         const double&, const double&, const double&, const double&, 
                         const double&, const double&, const int&);
  void F77NAME(hllcflux)(const int&, const double&, const double&, const double&, double*,
                         const double&, double*, double*, double*, double*, double*,
                         const double&, const double&, const double&, const double&,
                         const double&, const double&, const int&);
  void F77NAME(hllcflux1)(const double&, const double&, const double&, double*, 
			  const double&, double*, double*, double*,
			  const double&, const double&, const double&,
			  const double&, const int&);
  void F77NAME(hllcjac)(const int&, const double&, const double&, const double&, double*,
			const double&, double*, double*, double*, double*, const double&,
			const double&, const double&, const double&, const int&);


// Included (MB)
  void F77NAME(gxroeflux5)(const int&, const double&, const double&, const double&, const double&, double*, double*, const double&, const double&, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, const double&, const double&, const double&, const double&, const double&, const double&, const int&);
  void F77NAME(gxroeflux6)(const int&, const double&, const double&, const double&, const double&, double*, double*, const double&, const double&, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
  void F77NAME(gxroeflux6temp)(const int&, const double&, const double&, const double&, const double&, double*, double*, const double&, const double&, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
  void F77NAME(computedfluxoperators)(const int&, const double&, const double&, const double&, const double&, double*, const double&, double*, double*, double*, double*, double*, double (*)[3], double*, double (*)[7], double (*)[7]);
  void F77NAME(gxboundflux5)(const int&, const double&, double*, double*, const double&, const double&, double*, double*, double*, double*, double*);
  void F77NAME(gxboundflux5temp)(const int&, const double&, double*, double*, const double&, const double&, double*, double*, double*, double*, double*);
  void F77NAME(computedboundfluxoperator)(const int&, const double&, double*, const double&, double*, double*, double (*)[3], double (*)[1], double (*)[7]);

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void FluxFcnSGFDJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                       double *VL, double *VR, double *flux, bool useLimiter)
{
  int prec = useLimiter ? sprec.getPrecTag() : 0;
  if(prec) { // roe flux with spatial low mach preconditioner and numerical viscosity
    F77NAME(roeflux5)(0, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux,
                      sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, prec);
  }
  else if(gamma != 1) { // original implementation of the roe flux with numerical viscosity
    F77NAME(roeflux6)(0, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux);
  }
  else { // faster implementation of the roe flux
    F77NAME(roeflux7)(0, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
  }
}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGFDJacRoeEuler3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, bool useLimiter)
{

  F77NAME(gxroeflux5)(0, gamma, vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VL, dVL, VR, dVR, VR, dVR, flux, dFlux, sprec.getMinMach(), 0.0*dMach, sprec.getSlope(), sprec.getCutOffMach(), irey, dIrey, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnSGApprJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                         double *VL, double *VR, double *flux, bool useLimiter)
{
  int prec = useLimiter ? sprec.getPrecTag() : 0;
  if(prec) { // roe flux with spatial low mach preconditioner, numerical viscosity and galerkin advection operator
    F77NAME(roeflux5)(0, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux,
                      sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, prec);
  }
  else if(gamma != 1 || rshift) { // original implementation of the roe flux with numerical viscosity and galerkin advection operator
    F77NAME(roeflux6)(0, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux);
  }
  else { // faster implementation of the roe flux
    F77NAME(roeflux7)(0, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
  }
}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGApprJacRoeEuler3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, bool useLimiter)
{

  F77NAME(gxroeflux5)(0, gamma, vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VL+rshift, dVL+rshift, VR, dVR, VR+rshift, dVR+rshift, flux, dFlux, sprec.getMinMach(), 0.0*dMach, sprec.getSlope(), sprec.getCutOffMach(), irey, dIrey, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

template<int dim>
inline
void roejacappr3Dgas(int type, double gamma, VarFcnBase *vf, double vfgam, double vfp, 
                  FluxFcnBase::Type localTypeJac, double* normal, double normalVel, double* VL, 
		  double* VR, double* jacL, double* jacR, double irey, double betaRef, 
		  double k1, double cmach, double shockreducer, double length,
                  int prec)
{
  const int dimm1 = dim-1;
  const int dimm2 = dim-2;
  const int dim2 = dim*dim;

  double dfdVL[dim2], dfdVR[dim2];
  double dfdUL[dim2], dfdUR[dim2];
  double n[3] = {normal[0], normal[1], normal[2]};

  F77NAME(roejac6)(type, gamma, vfgam, vfp, n, normalVel, VL, VR, dfdUL, 1, betaRef, k1, cmach, shockreducer, irey, length, prec);
  F77NAME(roejac6)(type, gamma, vfgam, vfp, n, normalVel, VR, VL, dfdUR, 2, betaRef, k1, cmach, shockreducer, irey, length, prec);

  if (type == 1 || type == 2) {
    vf->postMultiplyBydUdV(VL, dfdUL, dfdVL);
    vf->postMultiplyBydUdV(VR, dfdUR, dfdVR);

    double f1;
    F77NAME(roeflux1)(gamma, vfgam, vfp, n, normalVel, VL, VR, &f1, betaRef,k1,cmach,irey,prec);

    if (type == 1) {
      if (f1 >= 0.0) {
        dfdVL[dim2-6] = dfdVL[0]*VL[dimm1]; dfdVR[dim2-6] = dfdVR[0]*VL[dimm1];
        dfdVL[dim2-5] = dfdVL[1]*VL[dimm1]; dfdVR[dim2-5] = dfdVR[1]*VL[dimm1];
        dfdVL[dim2-4] = dfdVL[2]*VL[dimm1]; dfdVR[dim2-4] = dfdVR[2]*VL[dimm1];
        dfdVL[dim2-3] = dfdVL[3]*VL[dimm1]; dfdVR[dim2-3] = dfdVR[3]*VL[dimm1];
        dfdVL[dim2-2] = dfdVL[4]*VL[dimm1]; dfdVR[dim2-2] = dfdVR[4]*VL[dimm1];
        dfdVL[dim2-1] = f1;                 dfdVR[dim2-1] = 0.0;
      }
      else {
        dfdVL[dim2-6] = dfdVL[0]*VR[dimm1]; dfdVR[dim2-6] = dfdVR[0]*VR[dimm1];
        dfdVL[dim2-5] = dfdVL[1]*VR[dimm1]; dfdVR[dim2-5] = dfdVR[1]*VR[dimm1];
        dfdVL[dim2-4] = dfdVL[2]*VR[dimm1]; dfdVR[dim2-4] = dfdVR[2]*VR[dimm1];
        dfdVL[dim2-3] = dfdVL[3]*VR[dimm1]; dfdVR[dim2-3] = dfdVR[3]*VR[dimm1];
        dfdVL[dim2-2] = dfdVL[4]*VR[dimm1]; dfdVR[dim2-2] = dfdVR[4]*VR[dimm1];
        dfdVL[dim2-1] = 0.0;                dfdVR[dim2-1] = f1;
      }
    }
    else if (type == 2) {
      if (f1 >= 0.0) {
        dfdVL[dim2-14] = dfdVL[0]*VL[dimm2]; dfdVR[dim2-14] = dfdVR[0]*VL[dimm2];
        dfdVL[dim2-13] = dfdVL[1]*VL[dimm2]; dfdVR[dim2-13] = dfdVR[1]*VL[dimm2];
        dfdVL[dim2-12] = dfdVL[2]*VL[dimm2]; dfdVR[dim2-12] = dfdVR[2]*VL[dimm2];
        dfdVL[dim2-11] = dfdVL[3]*VL[dimm2]; dfdVR[dim2-11] = dfdVR[3]*VL[dimm2];
        dfdVL[dim2-10] = dfdVL[4]*VL[dimm2]; dfdVR[dim2-10] = dfdVR[4]*VL[dimm2];
        dfdVL[dim2-9] = f1;                  dfdVR[dim2-9] = 0.0;
        dfdVL[dim2-8] = 0.0;                 dfdVR[dim2-8] = 0.0;

        dfdVL[dim2-7] = dfdVL[0]*VL[dimm1]; dfdVR[dim2-7] = dfdVR[0]*VL[dimm1];
        dfdVL[dim2-6] = dfdVL[1]*VL[dimm1]; dfdVR[dim2-6] = dfdVR[1]*VL[dimm1];
        dfdVL[dim2-5] = dfdVL[2]*VL[dimm1]; dfdVR[dim2-5] = dfdVR[2]*VL[dimm1];
        dfdVL[dim2-4] = dfdVL[3]*VL[dimm1]; dfdVR[dim2-4] = dfdVR[3]*VL[dimm1];
        dfdVL[dim2-3] = dfdVL[4]*VL[dimm1]; dfdVR[dim2-3] = dfdVR[4]*VL[dimm1];
        dfdVL[dim2-2] = 0.0;                dfdVR[dim2-2] = 0.0;
        dfdVL[dim2-1] = f1;                 dfdVR[dim2-1] = 0.0;
      }
      else {
        dfdVL[dim2-14] = dfdVL[0]*VR[dimm2]; dfdVR[dim2-14] = dfdVR[0]*VR[dimm2];
        dfdVL[dim2-13] = dfdVL[1]*VR[dimm2]; dfdVR[dim2-13] = dfdVR[1]*VR[dimm2];
        dfdVL[dim2-12] = dfdVL[2]*VR[dimm2]; dfdVR[dim2-12] = dfdVR[2]*VR[dimm2];
        dfdVL[dim2-11] = dfdVL[3]*VR[dimm2]; dfdVR[dim2-11] = dfdVR[3]*VR[dimm2];
        dfdVL[dim2-10] = dfdVL[4]*VR[dimm2]; dfdVR[dim2-10] = dfdVR[4]*VR[dimm2];
        dfdVL[dim2-9] = 0.0;                 dfdVR[dim2-9] = f1;
        dfdVL[dim2-8] = 0.0;                 dfdVR[dim2-8] = 0.0;

        dfdVL[dim2-7] = dfdVL[0]*VR[dimm1]; dfdVR[dim2-7] = dfdVR[0]*VR[dimm1];
        dfdVL[dim2-6] = dfdVL[1]*VR[dimm1]; dfdVR[dim2-6] = dfdVR[1]*VR[dimm1];
        dfdVL[dim2-5] = dfdVL[2]*VR[dimm1]; dfdVR[dim2-5] = dfdVR[2]*VR[dimm1];
        dfdVL[dim2-4] = dfdVL[3]*VR[dimm1]; dfdVR[dim2-4] = dfdVR[3]*VR[dimm1];
        dfdVL[dim2-3] = dfdVL[4]*VR[dimm1]; dfdVR[dim2-3] = dfdVR[4]*VR[dimm1];
        dfdVL[dim2-2] = 0.0;                dfdVR[dim2-2] = 0.0;
        dfdVL[dim2-1] = 0.0;                dfdVR[dim2-1] = f1;
      }
    }
    vf->postMultiplyBydVdU(VL, dfdVL, dfdUL);
    vf->postMultiplyBydVdU(VR, dfdVR, dfdUR);
  }

  int k;

  if (localTypeJac == FluxFcnBase::CONSERVATIVE) {
    for (k=0; k<dim2; ++k) {
      jacL[k] = dfdUL[k];
      jacR[k] = dfdUR[k];
    }
  }
  else {
    vf->postMultiplyBydUdV(VL, dfdUL, jacL);
    vf->postMultiplyBydUdV(VR, dfdUR, jacR);
  }

}

//------------------------------------------------------------------------------

void FluxFcnSGApprJacRoeEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                                  double *VL, double *VR, double *jacL, double *jacR, bool useLimiter)
{
  roejacappr3Dgas<5>(0, gamma, vf, vf->getGamma(), vf->getPressureConstant(), typeJac, normal, normalVel, VL, VR, jacL, jacR,
                     irey, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), length,
                     useLimiter ? sprec.getPrecTag() : 0);
}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnSGExactJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                          double *VL, double *VR, double *flux, bool useLimiter)
{
  int prec = useLimiter ? sprec.getPrecTag() : 0;
  if(prec) { // roe flux with spatial low mach preconditioner and numerical viscosity
    F77NAME(roeflux5)(0, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux,
                      sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, prec);
  }
  else if(gamma != 1) { // original implementation of the roe flux with numerical viscosity
    F77NAME(roeflux6)(0, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux);
  }
  else { // faster implementation of the roe flux
    F77NAME(roeflux7)(0, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
  }
}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGExactJacRoeEuler3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, bool useLimiter)
{

//  fprintf(stderr," ... in FluxFcnSGExactJacRoeEuler3D::computeDerivative 001\n");
  F77NAME(gxroeflux6)(0, gamma, vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VL, dVL, VR, dVR, VR, dVR, flux, dFlux);
/*
  double dFluxdNormal[7][3] = {0}, dFluxdNormalVel[7] = {0}, dFluxdVL[7][7] = {0}, dFluxdVR[7][7] = {0};
  F77NAME(computedfluxoperators)(0, gamma, vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux, dFluxdNormal, dFluxdNormalVel, dFluxdVL, dFluxdVR);


  double dFLUX2[7] = {0};

  for(int k=0; k<5; ++k) {
    for(int p=0; p<3; ++p)
      dFLUX2[k] += dFluxdNormal[k][p]*dNormal[p];
    for(int p=0; p<7; ++p)
      dFLUX2[k] += dFluxdVL[k][p]*dVL[p] + dFluxdVR[k][p]*dVR[p];
  }

  double dFluxnorm(0), diff2(0), dFLUX2norm(0);
  for(int k=0; k<5; ++k) {
    diff2+= (dFlux[k]-dFLUX2[k])*(dFlux[k]-dFLUX2[k]);
    dFLUX2norm += dFLUX2[k]*dFLUX2[k];
    dFluxnorm += dFlux[k]*dFlux[k];
  }
  diff2 = sqrt(diff2);
  dFLUX2norm = sqrt(dFLUX2norm);
  dFluxnorm = sqrt(dFluxnorm);
  if(dFluxnorm != 0) fprintf(stderr, "... rel. diff2 = %e, dFluxnorm = %e, dFLUX2norm = %e\n", diff2/dFluxnorm, dFluxnorm, dFLUX2norm); 
  else fprintf(stderr, "... abs. diff2 = %e, dFluxnorm = %e, dFLUX2norm = %e\n", diff2, dFluxnorm, dFLUX2norm);
*/
/*
  double dFLUX[7];
  F77NAME(gxroeflux6temp)(0, gamma, vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VL, dVL, VR, dVR, VR, dVR, flux, dFLUX);
  double dpstiff = vf->getDerivativeOfPressureConstant();
 
  if(dpstiff != 0) fprintf(stderr," ... dpstiff is %e\n", vf->getDerivativeOfPressureConstant());
  double diffsq(0), dfluxsqnorm(0);
  for(int i=0; i<5; ++i) {
    diffsq += (dFlux[i] - dFLUX[i])*(dFlux[i] - dFLUX[i]);
    dfluxsqnorm += dFlux[i]*dFlux[i];
  }
  double diff = std::sqrt(diffsq);
  double dfluxnorm = std::sqrt(dfluxsqnorm);
  if(dfluxsqnorm != 0) fprintf(stderr," rel. diff is %e and dfluxnorm is %e\n", diff/dfluxnorm, dfluxnorm);
  else fprintf(stderr, " abs. diff is %e\n", diff);
*/

}

//------------------------------------------------------------------------------

// Included (YC)
void FluxFcnSGExactJacRoeEuler3D::compute_dFluxdNormal_dFluxdNormalVel_dFluxdVL_dFluxdVR(double *normal, double normalVel, double *VL, double *VR, 
                                                                                         double dMach, double *flux, double dFluxdNormal[7][3], 
                                                                                         double *dFluxdNormalVel, double dFluxdVL[7][7], double dFluxdVR[7][7])
{
  F77NAME(computedfluxoperators)(0, gamma, vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux, dFluxdNormal, dFluxdNormalVel, dFluxdVL, dFluxdVR);
}
//------------------------------------------------------------------------------

template<int dim>
inline
void roejacexact3D(int type, double gamma, VarFcnBase *vf, FluxFcnBase::Type localTypeJac, double* normal, 
		   double normalVel, double* VL, double* VR, double* jacL, double* jacR, double irey, double betaRef,
                   double k1, double cmach, double shockreducer, double length, int prec)
{

  const int dimm1 = dim-1;
  const int dimm2 = dim-2;
  const int dim2 = dim*dim;

  double dfdVL[dim2], dfdVR[dim2];
  double n[3] = {normal[0], normal[1], normal[2]};
  
  if(prec) {
    roeturkeljac5(type, gamma, vf->getGamma(), vf->getPressureConstant(), n, normalVel, VL, VR, dfdVL,
                  betaRef, k1, cmach, shockreducer, irey, length, prec);
    n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2];
    roeturkeljac5(type, gamma, vf->getGamma(), vf->getPressureConstant(), n, -normalVel, VR, VL, dfdVR,
                  betaRef, k1, cmach, shockreducer, irey, length, prec);
  }
  else {
    F77NAME(roejac5)(type, gamma, vf->getGamma(), vf->getPressureConstant(), n, normalVel, VL, VR, dfdVL);
    n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2];
    F77NAME(roejac5)(type, gamma, vf->getGamma(), vf->getPressureConstant(), n, -normalVel, VR, VL, dfdVR);
  }

  int k;
  for (k=0; k<dim2; ++k)
    dfdVR[k] = -dfdVR[k];

  if (type == 1 || type == 2) {
    double f1;
    F77NAME(roeflux1)(gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, &f1, betaRef, k1, cmach, irey, prec);

    if (type == 1) {
      if (f1 >= 0.0) {
	dfdVL[dim2-6] = dfdVL[0]*VL[dimm1]; dfdVR[dim2-6] = dfdVR[0]*VL[dimm1];
	dfdVL[dim2-5] = dfdVL[1]*VL[dimm1]; dfdVR[dim2-5] = dfdVR[1]*VL[dimm1];
	dfdVL[dim2-4] = dfdVL[2]*VL[dimm1]; dfdVR[dim2-4] = dfdVR[2]*VL[dimm1];
	dfdVL[dim2-3] = dfdVL[3]*VL[dimm1]; dfdVR[dim2-3] = dfdVR[3]*VL[dimm1];
	dfdVL[dim2-2] = dfdVL[4]*VL[dimm1]; dfdVR[dim2-2] = dfdVR[4]*VL[dimm1];
	dfdVL[dim2-1] = f1;                 dfdVR[dim2-1] = 0.0;
      }
      else {
	dfdVL[dim2-6] = dfdVL[0]*VR[dimm1]; dfdVR[dim2-6] = dfdVR[0]*VR[dimm1];
	dfdVL[dim2-5] = dfdVL[1]*VR[dimm1]; dfdVR[dim2-5] = dfdVR[1]*VR[dimm1];
	dfdVL[dim2-4] = dfdVL[2]*VR[dimm1]; dfdVR[dim2-4] = dfdVR[2]*VR[dimm1];
	dfdVL[dim2-3] = dfdVL[3]*VR[dimm1]; dfdVR[dim2-3] = dfdVR[3]*VR[dimm1];
	dfdVL[dim2-2] = dfdVL[4]*VR[dimm1]; dfdVR[dim2-2] = dfdVR[4]*VR[dimm1];
	dfdVL[dim2-1] = 0.0;                dfdVR[dim2-1] = f1;
      }
    }
    else if (type == 2) {
      if (f1 >= 0.0) {
	dfdVL[dim2-14] = dfdVL[0]*VL[dimm2]; dfdVR[dim2-14] = dfdVR[0]*VL[dimm2];
	dfdVL[dim2-13] = dfdVL[1]*VL[dimm2]; dfdVR[dim2-13] = dfdVR[1]*VL[dimm2];
	dfdVL[dim2-12] = dfdVL[2]*VL[dimm2]; dfdVR[dim2-12] = dfdVR[2]*VL[dimm2];
	dfdVL[dim2-11] = dfdVL[3]*VL[dimm2]; dfdVR[dim2-11] = dfdVR[3]*VL[dimm2];
	dfdVL[dim2-10] = dfdVL[4]*VL[dimm2]; dfdVR[dim2-10] = dfdVR[4]*VL[dimm2];
	dfdVL[dim2-9] = f1;                  dfdVR[dim2-9] = 0.0;
	dfdVL[dim2-8] = 0.0;                 dfdVR[dim2-8] = 0.0;

	dfdVL[dim2-7] = dfdVL[0]*VL[dimm1]; dfdVR[dim2-7] = dfdVR[0]*VL[dimm1];
	dfdVL[dim2-6] = dfdVL[1]*VL[dimm1]; dfdVR[dim2-6] = dfdVR[1]*VL[dimm1];
	dfdVL[dim2-5] = dfdVL[2]*VL[dimm1]; dfdVR[dim2-5] = dfdVR[2]*VL[dimm1];
	dfdVL[dim2-4] = dfdVL[3]*VL[dimm1]; dfdVR[dim2-4] = dfdVR[3]*VL[dimm1];
	dfdVL[dim2-3] = dfdVL[4]*VL[dimm1]; dfdVR[dim2-3] = dfdVR[4]*VL[dimm1];
	dfdVL[dim2-2] = 0.0;                dfdVR[dim2-2] = 0.0;
	dfdVL[dim2-1] = f1;                 dfdVR[dim2-1] = 0.0;
      }
      else {
	dfdVL[dim2-14] = dfdVL[0]*VR[dimm2]; dfdVR[dim2-14] = dfdVR[0]*VR[dimm2];
	dfdVL[dim2-13] = dfdVL[1]*VR[dimm2]; dfdVR[dim2-13] = dfdVR[1]*VR[dimm2];
	dfdVL[dim2-12] = dfdVL[2]*VR[dimm2]; dfdVR[dim2-12] = dfdVR[2]*VR[dimm2];
	dfdVL[dim2-11] = dfdVL[3]*VR[dimm2]; dfdVR[dim2-11] = dfdVR[3]*VR[dimm2];
	dfdVL[dim2-10] = dfdVL[4]*VR[dimm2]; dfdVR[dim2-10] = dfdVR[4]*VR[dimm2];
	dfdVL[dim2-9] = 0.0;                 dfdVR[dim2-9] = f1;
	dfdVL[dim2-8] = 0.0;                 dfdVR[dim2-8] = 0.0;

	dfdVL[dim2-7] = dfdVL[0]*VR[dimm1]; dfdVR[dim2-7] = dfdVR[0]*VR[dimm1];
	dfdVL[dim2-6] = dfdVL[1]*VR[dimm1]; dfdVR[dim2-6] = dfdVR[1]*VR[dimm1];
	dfdVL[dim2-5] = dfdVL[2]*VR[dimm1]; dfdVR[dim2-5] = dfdVR[2]*VR[dimm1];
	dfdVL[dim2-4] = dfdVL[3]*VR[dimm1]; dfdVR[dim2-4] = dfdVR[3]*VR[dimm1];
	dfdVL[dim2-3] = dfdVL[4]*VR[dimm1]; dfdVR[dim2-3] = dfdVR[4]*VR[dimm1];
	dfdVL[dim2-2] = 0.0;                dfdVR[dim2-2] = 0.0;
	dfdVL[dim2-1] = 0.0;                dfdVR[dim2-1] = f1;
      }
    }
  }

  if (localTypeJac == FluxFcnBase::CONSERVATIVE) {
    vf->postMultiplyBydVdU(VL, dfdVL, jacL);
    vf->postMultiplyBydVdU(VR, dfdVR, jacR);
  }
  else {
    for (k=0; k<dim2; ++k) { 
      jacL[k] = dfdVL[k]; 
      jacR[k] = dfdVR[k]; 
    }
  }

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnSGExactJacRoeEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                                   double *VL, double *VR, double *jacL, double *jacR, bool useLimiter)
{
  roejacexact3D<5>(0, gamma, vf, typeJac, normal, normalVel, VL, VR, jacL, jacR, irey, sprec.getMinMach(), sprec.getSlope(),
                   sprec.getCutOffMach(), sprec.getShockParameter(), length, useLimiter ? sprec.getPrecTag() : 0);
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

template<int dim>
inline
void hllejacappr3Dgas(int type, double gamma, VarFcnBase *vf, double vfgam, double vfp, 
                  FluxFcnBase::Type localTypeJac, double* normal, double normalVel, double* VL, 
		  double* VR, double* jacL, double* jacR, double irey, double betaRef, 
		  double k1, double cmach, int prec, double shockreducer, double length)
{
  const int dimm1 = dim-1;
  const int dimm2 = dim-2;
  const int dim2 = dim*dim;

  double dfdVL[dim2], dfdVR[dim2];
  double dfdUL[dim2], dfdUR[dim2];
  double n[3] = {normal[0], normal[1], normal[2]};

  F77NAME(hllejac)(type, gamma, vfgam, vfp, n, normalVel, VL, VR, dfdUL,1, betaRef, k1, cmach, shockreducer, irey, length, prec);
  F77NAME(hllejac)(type, gamma, vfgam, vfp, n, normalVel, VR, VL, dfdUR,2, betaRef, k1, cmach, shockreducer, irey, length, prec);

  if (type == 1 || type == 2) {
    vf->postMultiplyBydUdV(VL, dfdUL, dfdVL);
    vf->postMultiplyBydUdV(VR, dfdUR, dfdVR);

    double f1;
    F77NAME(hlleflux1)(gamma, vfgam, vfp, n, normalVel, VL, VR, &f1, betaRef, k1, cmach, shockreducer, irey, length, prec);


    if (type == 1) {
      if (f1 >= 0.0) {
        dfdVL[dim2-6] = dfdVL[0]*VL[dimm1]; dfdVR[dim2-6] = dfdVR[0]*VL[dimm1];
        dfdVL[dim2-5] = dfdVL[1]*VL[dimm1]; dfdVR[dim2-5] = dfdVR[1]*VL[dimm1];
        dfdVL[dim2-4] = dfdVL[2]*VL[dimm1]; dfdVR[dim2-4] = dfdVR[2]*VL[dimm1];
        dfdVL[dim2-3] = dfdVL[3]*VL[dimm1]; dfdVR[dim2-3] = dfdVR[3]*VL[dimm1];
        dfdVL[dim2-2] = dfdVL[4]*VL[dimm1]; dfdVR[dim2-2] = dfdVR[4]*VL[dimm1];
        dfdVL[dim2-1] = f1;                 dfdVR[dim2-1] = 0.0;
      }
      else {
        dfdVL[dim2-6] = dfdVL[0]*VR[dimm1]; dfdVR[dim2-6] = dfdVR[0]*VR[dimm1];
        dfdVL[dim2-5] = dfdVL[1]*VR[dimm1]; dfdVR[dim2-5] = dfdVR[1]*VR[dimm1];
        dfdVL[dim2-4] = dfdVL[2]*VR[dimm1]; dfdVR[dim2-4] = dfdVR[2]*VR[dimm1];
        dfdVL[dim2-3] = dfdVL[3]*VR[dimm1]; dfdVR[dim2-3] = dfdVR[3]*VR[dimm1];
        dfdVL[dim2-2] = dfdVL[4]*VR[dimm1]; dfdVR[dim2-2] = dfdVR[4]*VR[dimm1];
        dfdVL[dim2-1] = 0.0;                dfdVR[dim2-1] = f1;
      }
    }
    else if (type == 2) {
      if (f1 >= 0.0) {
        dfdVL[dim2-14] = dfdVL[0]*VL[dimm2]; dfdVR[dim2-14] = dfdVR[0]*VL[dimm2];
        dfdVL[dim2-13] = dfdVL[1]*VL[dimm2]; dfdVR[dim2-13] = dfdVR[1]*VL[dimm2];
        dfdVL[dim2-12] = dfdVL[2]*VL[dimm2]; dfdVR[dim2-12] = dfdVR[2]*VL[dimm2];
        dfdVL[dim2-11] = dfdVL[3]*VL[dimm2]; dfdVR[dim2-11] = dfdVR[3]*VL[dimm2];
        dfdVL[dim2-10] = dfdVL[4]*VL[dimm2]; dfdVR[dim2-10] = dfdVR[4]*VL[dimm2];
        dfdVL[dim2-9] = f1;                  dfdVR[dim2-9] = 0.0;
        dfdVL[dim2-8] = 0.0;                 dfdVR[dim2-8] = 0.0;

        dfdVL[dim2-7] = dfdVL[0]*VL[dimm1]; dfdVR[dim2-7] = dfdVR[0]*VL[dimm1];
        dfdVL[dim2-6] = dfdVL[1]*VL[dimm1]; dfdVR[dim2-6] = dfdVR[1]*VL[dimm1];
        dfdVL[dim2-5] = dfdVL[2]*VL[dimm1]; dfdVR[dim2-5] = dfdVR[2]*VL[dimm1];
        dfdVL[dim2-4] = dfdVL[3]*VL[dimm1]; dfdVR[dim2-4] = dfdVR[3]*VL[dimm1];
        dfdVL[dim2-3] = dfdVL[4]*VL[dimm1]; dfdVR[dim2-3] = dfdVR[4]*VL[dimm1];
        dfdVL[dim2-2] = 0.0;                dfdVR[dim2-2] = 0.0;
        dfdVL[dim2-1] = f1;                 dfdVR[dim2-1] = 0.0;
      }
      else {
        dfdVL[dim2-14] = dfdVL[0]*VR[dimm2]; dfdVR[dim2-14] = dfdVR[0]*VR[dimm2];
        dfdVL[dim2-13] = dfdVL[1]*VR[dimm2]; dfdVR[dim2-13] = dfdVR[1]*VR[dimm2];
        dfdVL[dim2-12] = dfdVL[2]*VR[dimm2]; dfdVR[dim2-12] = dfdVR[2]*VR[dimm2];
        dfdVL[dim2-11] = dfdVL[3]*VR[dimm2]; dfdVR[dim2-11] = dfdVR[3]*VR[dimm2];
        dfdVL[dim2-10] = dfdVL[4]*VR[dimm2]; dfdVR[dim2-10] = dfdVR[4]*VR[dimm2];
        dfdVL[dim2-9] = 0.0;                 dfdVR[dim2-9] = f1;
        dfdVL[dim2-8] = 0.0;                 dfdVR[dim2-8] = 0.0;

        dfdVL[dim2-7] = dfdVL[0]*VR[dimm1]; dfdVR[dim2-7] = dfdVR[0]*VR[dimm1];
        dfdVL[dim2-6] = dfdVL[1]*VR[dimm1]; dfdVR[dim2-6] = dfdVR[1]*VR[dimm1];
        dfdVL[dim2-5] = dfdVL[2]*VR[dimm1]; dfdVR[dim2-5] = dfdVR[2]*VR[dimm1];
        dfdVL[dim2-4] = dfdVL[3]*VR[dimm1]; dfdVR[dim2-4] = dfdVR[3]*VR[dimm1];
        dfdVL[dim2-3] = dfdVL[4]*VR[dimm1]; dfdVR[dim2-3] = dfdVR[4]*VR[dimm1];
        dfdVL[dim2-2] = 0.0;                dfdVR[dim2-2] = 0.0;
        dfdVL[dim2-1] = 0.0;                dfdVR[dim2-1] = f1;
      }
    }
    vf->postMultiplyBydVdU(VL, dfdVL, dfdUL);
    vf->postMultiplyBydVdU(VR, dfdVR, dfdUR);
  }

  int k;

  if (localTypeJac == FluxFcnBase::CONSERVATIVE) {
    for (k=0; k<dim2; ++k) {
      jacL[k] = dfdUL[k];
      jacR[k] = dfdUR[k];
    }
  }
  else {
    vf->postMultiplyBydUdV(VL, dfdUL, jacL);
    vf->postMultiplyBydUdV(VR, dfdUR, jacR);
  }

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void FluxFcnSGFDJacHLLEEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *flux, bool useLimiter)
{

  F77NAME(hlleflux)(0, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGApprJacHLLEEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                       double *VL, double *VR, double *flux, bool useLimiter)
{
  F77NAME(hlleflux)(0, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGApprJacHLLEEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                                double *VL, double *VR,
                                                double *jacL, double *jacR, bool useLimiter)
{
  
  hllejacappr3Dgas<5>(0, gamma, vf, vf->getGamma(), vf->getPressureConstant(), typeJac, normal, normalVel, VL, VR, jacL, jacR, irey, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), useLimiter ? sprec.getPrecTag() : 0, sprec.getShockParameter(), length);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void hllcjacappr3Dgas(int type, double gamma, VarFcnBase *vf, double vfgam, double vfp, 
                  FluxFcnBase::Type localTypeJac, double* normal, double normalVel, double* VL, 
		  double* VR, double* jacL, double* jacR, double irey, double betaRef, 
		  double k1, double cmach, int prec)
{
  const int dimm1 = dim-1;
  const int dimm2 = dim-2;
  const int dim2 = dim*dim;

  double dfdVL[dim2], dfdVR[dim2];
  double dfdUL[dim2], dfdUR[dim2];
  double n[3] = {normal[0], normal[1], normal[2]};

  F77NAME(hllcjac)(type, gamma, vfgam, vfp, n, normalVel, VL, VR, dfdUL, dfdUR, betaRef, k1, cmach, irey, prec);

  if (type == 1 || type == 2) {
    vf->postMultiplyBydUdV(VL, dfdUL, dfdVL);
    vf->postMultiplyBydUdV(VR, dfdUR, dfdVR);

    double f1;
    F77NAME(hllcflux1)(gamma, vfgam, vfp, n, normalVel, VL, VR, &f1, betaRef,k1,cmach,irey,prec);


    if (type == 1) {
      if (f1 >= 0.0) {
        dfdVL[dim2-6] = dfdVL[0]*VL[dimm1]; dfdVR[dim2-6] = dfdVR[0]*VL[dimm1];
        dfdVL[dim2-5] = dfdVL[1]*VL[dimm1]; dfdVR[dim2-5] = dfdVR[1]*VL[dimm1];
        dfdVL[dim2-4] = dfdVL[2]*VL[dimm1]; dfdVR[dim2-4] = dfdVR[2]*VL[dimm1];
        dfdVL[dim2-3] = dfdVL[3]*VL[dimm1]; dfdVR[dim2-3] = dfdVR[3]*VL[dimm1];
        dfdVL[dim2-2] = dfdVL[4]*VL[dimm1]; dfdVR[dim2-2] = dfdVR[4]*VL[dimm1];
        dfdVL[dim2-1] = f1;                 dfdVR[dim2-1] = 0.0;
      }
      else {
        dfdVL[dim2-6] = dfdVL[0]*VR[dimm1]; dfdVR[dim2-6] = dfdVR[0]*VR[dimm1];
        dfdVL[dim2-5] = dfdVL[1]*VR[dimm1]; dfdVR[dim2-5] = dfdVR[1]*VR[dimm1];
        dfdVL[dim2-4] = dfdVL[2]*VR[dimm1]; dfdVR[dim2-4] = dfdVR[2]*VR[dimm1];
        dfdVL[dim2-3] = dfdVL[3]*VR[dimm1]; dfdVR[dim2-3] = dfdVR[3]*VR[dimm1];
        dfdVL[dim2-2] = dfdVL[4]*VR[dimm1]; dfdVR[dim2-2] = dfdVR[4]*VR[dimm1];
        dfdVL[dim2-1] = 0.0;                dfdVR[dim2-1] = f1;
      }
    }
    else if (type == 2) {
      if (f1 >= 0.0) {
        dfdVL[dim2-14] = dfdVL[0]*VL[dimm2]; dfdVR[dim2-14] = dfdVR[0]*VL[dimm2];
        dfdVL[dim2-13] = dfdVL[1]*VL[dimm2]; dfdVR[dim2-13] = dfdVR[1]*VL[dimm2];
        dfdVL[dim2-12] = dfdVL[2]*VL[dimm2]; dfdVR[dim2-12] = dfdVR[2]*VL[dimm2];
        dfdVL[dim2-11] = dfdVL[3]*VL[dimm2]; dfdVR[dim2-11] = dfdVR[3]*VL[dimm2];
        dfdVL[dim2-10] = dfdVL[4]*VL[dimm2]; dfdVR[dim2-10] = dfdVR[4]*VL[dimm2];
        dfdVL[dim2-9] = f1;                  dfdVR[dim2-9] = 0.0;
        dfdVL[dim2-8] = 0.0;                 dfdVR[dim2-8] = 0.0;

        dfdVL[dim2-7] = dfdVL[0]*VL[dimm1]; dfdVR[dim2-7] = dfdVR[0]*VL[dimm1];
        dfdVL[dim2-6] = dfdVL[1]*VL[dimm1]; dfdVR[dim2-6] = dfdVR[1]*VL[dimm1];
        dfdVL[dim2-5] = dfdVL[2]*VL[dimm1]; dfdVR[dim2-5] = dfdVR[2]*VL[dimm1];
        dfdVL[dim2-4] = dfdVL[3]*VL[dimm1]; dfdVR[dim2-4] = dfdVR[3]*VL[dimm1];
        dfdVL[dim2-3] = dfdVL[4]*VL[dimm1]; dfdVR[dim2-3] = dfdVR[4]*VL[dimm1];
        dfdVL[dim2-2] = 0.0;                dfdVR[dim2-2] = 0.0;
        dfdVL[dim2-1] = f1;                 dfdVR[dim2-1] = 0.0;
      }
      else {
        dfdVL[dim2-14] = dfdVL[0]*VR[dimm2]; dfdVR[dim2-14] = dfdVR[0]*VR[dimm2];
        dfdVL[dim2-13] = dfdVL[1]*VR[dimm2]; dfdVR[dim2-13] = dfdVR[1]*VR[dimm2];
        dfdVL[dim2-12] = dfdVL[2]*VR[dimm2]; dfdVR[dim2-12] = dfdVR[2]*VR[dimm2];
        dfdVL[dim2-11] = dfdVL[3]*VR[dimm2]; dfdVR[dim2-11] = dfdVR[3]*VR[dimm2];
        dfdVL[dim2-10] = dfdVL[4]*VR[dimm2]; dfdVR[dim2-10] = dfdVR[4]*VR[dimm2];
        dfdVL[dim2-9] = 0.0;                 dfdVR[dim2-9] = f1;
        dfdVL[dim2-8] = 0.0;                 dfdVR[dim2-8] = 0.0;

        dfdVL[dim2-7] = dfdVL[0]*VR[dimm1]; dfdVR[dim2-7] = dfdVR[0]*VR[dimm1];
        dfdVL[dim2-6] = dfdVL[1]*VR[dimm1]; dfdVR[dim2-6] = dfdVR[1]*VR[dimm1];
        dfdVL[dim2-5] = dfdVL[2]*VR[dimm1]; dfdVR[dim2-5] = dfdVR[2]*VR[dimm1];
        dfdVL[dim2-4] = dfdVL[3]*VR[dimm1]; dfdVR[dim2-4] = dfdVR[3]*VR[dimm1];
        dfdVL[dim2-3] = dfdVL[4]*VR[dimm1]; dfdVR[dim2-3] = dfdVR[4]*VR[dimm1];
        dfdVL[dim2-2] = 0.0;                dfdVR[dim2-2] = 0.0;
        dfdVL[dim2-1] = 0.0;                dfdVR[dim2-1] = f1;
      }
    }
    vf->postMultiplyBydVdU(VL, dfdVL, dfdUL);
    vf->postMultiplyBydVdU(VR, dfdVR, dfdUR);
  }

  int k;

  if (localTypeJac == FluxFcnBase::CONSERVATIVE) {
    for (k=0; k<dim2; ++k) {
      jacL[k] = dfdUL[k];
      jacR[k] = dfdUR[k];
    }
  }
  else {
    vf->postMultiplyBydUdV(VL, dfdUL, jacL);
    vf->postMultiplyBydUdV(VR, dfdUR, jacR);
  }

}

//------------------------------------------------------------------------------

void FluxFcnSGFDJacHLLCEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                     double *VL, double *VR, double *flux, bool useLimiter)
{

  F77NAME(hllcflux)(0, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGApprJacHLLCEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                       double *VL, double *VR, double *flux, bool useLimiter)
{

  F77NAME(hllcflux)(0, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGApprJacHLLCEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                                double *VL, double *VR, double *jacL, double *jacR, bool useLimiter)
{
  
  hllcjacappr3Dgas<5>(0, gamma, vf, vf->getGamma(), vf->getPressureConstant(), typeJac, normal, normalVel, VL, VR, jacL, jacR, irey, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGVanLeerEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                    double *VL, double *VR, double *flux, bool useLimiter)
{

  //compute the split fluxes
  double fPlus[5], fMinus[5];
  evalFlux(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, fPlus, 1);
  evalFlux(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VR, fMinus, -1);

  for (int i = 0; i < 5; i++)
    flux[i] = fPlus[i] + fMinus[i];

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGVanLeerEuler3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, bool useLimiter)
{

  //compute the split fluxes
  double fPlus[5], fMinus[5];
  evalFlux(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, fPlus, 1);
  evalFlux(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VR, fMinus, -1);

  for (int i = 0; i < 5; i++)
    flux[i] = fPlus[i] + fMinus[i];

  //compute the split fluxes
  double dFPlus[5], dFMinus[5];
  evalDerivativeOfFlux(vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, dMach, dFPlus, 1);
  evalDerivativeOfFlux(vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VR, dVR, dMach, dFMinus, -1);

  for (int i = 0; i < 5; i++)
    dFlux[i] = dFPlus[i] + dFMinus[i];

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnSGVanLeerEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                             double *VL, double *VR,
                                             double *jacL, double *jacR, bool useLimiter)
{

  //compute the jacs for left and right states
  evalJac(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, jacL, 1);
  evalJac(vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VR, jacR, -1);

}

//------------------------------------------------------------------------------

void FluxFcnSGWallEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *V, double *Ub, double *flux, bool useLimiter)
{

  flux[0] = 0.0;
  flux[1] = V[4] * normal[0];
  flux[2] = V[4] * normal[1];
  flux[3] = V[4] * normal[2];
  flux[4] = V[4] * normalVel;

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGWallEuler3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal, 
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

//  fprintf(stderr," ... in FluxFcnSGWallEuler3D::computeDerivative 001\n");
  dFlux[0] = 0.0;
  dFlux[1] = V[4] * dNormal[0];
  dFlux[2] = V[4] * dNormal[1];
  dFlux[3] = V[4] * dNormal[2];
  dFlux[4] = V[4] * dNormalVel;

}

//------------------------------------------------------------------------------

// Included (YC)
void FluxFcnSGWallEuler3D::computeDerivativeOperators
(
  double *normal, double normalVel, double *V, double *Ub, 
  double dFluxdNormal[7][3], double dFluxdNormalVel[7][1], double dFluxdUb[7][7]
)
{

//  fprintf(stderr," ... in FluxFcnSGWallEuler3D::computeDerivative 001\n");
//  dFlux[0] = 0.0;
//  dFlux[1] = V[4] * dNormal[0];
//  dFlux[2] = V[4] * dNormal[1];
//  dFlux[3] = V[4] * dNormal[2];
//  dFlux[4] = V[4] * dNormalVel;

  dFluxdNormal[1][0] = V[4];
  dFluxdNormal[2][1] = V[4];
  dFluxdNormal[3][2] = V[4];
  dFluxdNormalVel[4][0] = V[4];

}

//------------------------------------------------------------------------------

// Included (MB*)
void FluxFcnSGWallEuler3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                   double *V, double *Ub, double *jac, bool useLimiter)
{

//  flux[0] = 0.0;
//  flux[1] = V[4] * normal[0];
//  flux[2] = V[4] * normal[1];
//  flux[3] = V[4] * normal[2];
//  flux[4] = V[4] * normalVel;

  double dfdV[25];

  dfdV[0] = 0.0;
  dfdV[1] = 0.0;
  dfdV[2] = 0.0;
  dfdV[3] = 0.0;
  dfdV[4] = 0.0;

  dfdV[5] = 0.0;
  dfdV[6] = 0.0;
  dfdV[7] = 0.0;
  dfdV[8] = 0.0;
  dfdV[9] = normal[0];

  dfdV[10] = 0.0;
  dfdV[11] = 0.0;
  dfdV[12] = 0.0;
  dfdV[13] = 0.0;
  dfdV[14] = normal[1];

  dfdV[15] = 0.0;
  dfdV[16] = 0.0;
  dfdV[17] = 0.0;
  dfdV[18] = 0.0;
  dfdV[19] = normal[2];

  dfdV[20] = 0.0;
  dfdV[21] = 0.0;
  dfdV[22] = 0.0;
  dfdV[23] = 0.0;
  dfdV[24] = normalVel;

  if (typeJac == FluxFcnBase::CONSERVATIVE)
    vf->postMultiplyBydVdU(V, dfdV, jac);
  else
    for (int k=0; k<25; ++k)
      jac[k] = dfdV[k];

}

//------------------------------------------------------------------------------

void FluxFcnSGPorousWallEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *V, double *Ub, double *flux, bool useLimiter)
{

  double gam, gam1, invgam1, pstiff;
  gam = vf->getGamma();
  pstiff = vf->getPressureConstant();
  gam1 = gam - 1.0; 
  invgam1 = 1.0 / gam1;

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;
  double porosity = Ub[1];

  if (porosity != 0 ) {
    double rho,rhoun,u,v,w,p,q;
    if ( Ub[0] <= 0.0 ) {
      rho = invgam1 * (V[4]+pstiff) / Ub[4];
      double ujet = Ub[0]/rho;
      double un = ujet/porosity;
      rhoun = rho*(ujet - nVel);
      u = un*n[0];
      v = un*n[1];
      w = un*n[2];
      q = u*u + v*v + w*w;
      p = V[4];
    }
    else if ( Ub[0] > 0.0 ) {
      rho = V[0];
      double ujet = Ub[0]/rho;
      double Un = V[1]*n[0] + V[2]*n[1] + V[3]*n[2];
      double Vt[3] = {V[1] - Un*n[0],V[2] - Un*n[1],V[3] - Un*n[2]};
      rhoun = rho*(ujet - nVel);
      u = ujet*n[0] + Vt[0];
      v = ujet*n[1] + Vt[1];
      w = ujet*n[2] + Vt[2];
      double u1 = ujet*n[0]/porosity + Vt[0];
      double v1 = ujet*n[1]/porosity + Vt[1];
      double w1 = ujet*n[2]/porosity + Vt[2];
      q = u1*u1 + v1*v1 + w1*w1;
      p = V[4];
    } 
    flux[0] = S * rhoun;
    flux[1] = S * (rhoun*u + p*n[0]);
    flux[2] = S * (rhoun*v + p*n[1]);
    flux[3] = S * (rhoun*w + p*n[2]);
    flux[4] = S * (rhoun*(gam*invgam1*(p+pstiff)/rho + 0.5*q) + p*nVel);
  }
  else {
    flux[0] = 0.0;
    flux[1] = V[4] * normal[0];
    flux[2] = V[4] * normal[1];
    flux[3] = V[4] * normal[2];
    flux[4] = V[4] * normalVel;
  }

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGPorousWallEuler3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal, 
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

// not implemented 

}

//------------------------------------------------------------------------------

void FluxFcnSGGhidagliaEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                   double *V, double *Ub, double *flux, bool useLimiter)
{

  F77NAME(genbcfluxgas)(0, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

void FluxFcnSGModifiedGhidagliaEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                   double *V, double *Ub, double *flux, bool useLimiter)
{
  const int dim0 = 5; 
  F77NAME(genbcfluxgas_hh)(0, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux, 
                          hhcoeffptr, *(hhcoeffptr+3), *(hhcoeffptr+4), *(hhcoeffptr+5));

}

void FluxFcnSGModifiedGhidagliaEuler3D::computeJacobianFarfield(double length, double irey, double *normal, double normalVel, 
	 						        double *V, double *Ub, double *jac,
								bool useLimiter)
{

  const int dim = 5; 
  
  double ff[dim],ffp[dim];
  
  const double eps0 = 1.0e-6;

  compute(length, irey, normal, normalVel, V, Ub, ff,useLimiter);

  double olds = *(hhcoeffptr+4);
  *(hhcoeffptr+4) *=(1.0+eps0);
  double news = *(hhcoeffptr+4);
  compute(length, irey, normal, normalVel, V, Ub, ffp,useLimiter);
  
  for (int k = 0; k < dim; ++k)
    jac[k] = (ffp[k] - ff[k]) / (news-olds);
  
  *(hhcoeffptr+4) = olds;
 
}


//------------------------------------------------------------------------------

void FluxFcnSGInflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                   double *V, double *Ub, double *flux, bool useLimiter)
{

  F77NAME(boundflux5)(0, vf->getGamma(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void FluxFcnSGInflowEuler3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  
  F77NAME(gxboundflux5)(0, vf->getGamma(), normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------

// Included (YC)
inline
void FluxFcnSGInflowEuler3D::computeDerivativeOperators
(
  double *normal, double normalVel, double *V, double *Ub,
  double dFluxdNormal[7][3], double dFluxdNormalVel[7][1], double dFluxdUb[7][7]
)
{
  F77NAME(computedboundfluxoperator)(0, vf->getGamma(), normal, normalVel, V, Ub, dFluxdNormal, dFluxdNormalVel, dFluxdUb);
}

//------------------------------------------------------------------------------

void FluxFcnSGOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                    double *V, double *Ub, double *flux, bool useLimiter)
{

  F77NAME(boundflux5)(0, vf->getGamma(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void FluxFcnSGOutflowEuler3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{
//  fprintf(stderr," ... in FluxFcnSGOutflowEuler3D::computeDerivative 001\n");
/*
  dNormal[0] = 1;
  dNormal[1] = 2;
  dNormal[2] = 0.5;
  dNormalVel = 3.2;
  dUb[0] = 0.3;
  dUb[1] = 1.3;
  dUb[2] = 2.3;
  dUb[3] = 3.3;
  dUb[4] = 4.3;
  dUb[5] = 0.25;
*/
  F77NAME(gxboundflux5)(0, vf->getGamma(), normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);
/*
  double dFLUX[7]={0}, dfluxdnormal[7][3]={0}, dfluxdsig[7][1]={0}, dfluxdub[7][7]={0};
  F77NAME(computedboundfluxoperator)(0, vf->getGamma(), normal, normalVel, V, Ub, dfluxdnormal, dfluxdsig, dfluxdub);
//  F77NAME(gxboundflux5temp)(0, vf->getGamma(), normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

  for(int i=0; i<5; ++i) {
    dFLUX[i] += dfluxdsig[i][0]*dNormalVel;
    for(int j=0; j<3; ++j) {
      dFLUX[i] += dfluxdnormal[i][j]*dNormal[j];
    }
    for(int j=0; j<5; ++j) {
      dFLUX[i] += dfluxdub[i][j]*dUb[j];
    }
  }


  double diff(0), dFluxnorm(0);
  for(int i=0; i<5; ++i) {
    diff += (dFlux[i] - dFLUX[i])*(dFlux[i] - dFLUX[i]);
    dFluxnorm += dFlux[i]*dFlux[i];
  }
  diff = sqrt(diff);
  dFluxnorm = sqrt(dFluxnorm);
  if(dFluxnorm != 0) fprintf(stderr, " ... rel. diff is %e\n", diff/dFluxnorm);
  else fprintf(stderr, " ... abs. diff is %e\n", diff);
*/  

}

//------------------------------------------------------------------------------

// Included (YC)
inline
void FluxFcnSGOutflowEuler3D::computeDerivativeOperators
(
  double *normal, double normalVel, double *V, double *Ub,
  double dFluxdNormal[7][3], double dFluxdNormalVel[7][1], double dFluxdUb[7][7]
)
{
  F77NAME(computedboundfluxoperator)(0, vf->getGamma(), normal, normalVel, V, Ub, dFluxdNormal, dFluxdNormalVel, dFluxdUb);
}

//------------------------------------------------------------------------------

template<int dim>
inline
void influx3D(int type, VarFcnBase *vf, double* normal, 
	      double normalVel, double* V, double* Ub, double* flux)
{
  const int dimm1 = dim-1;
  const int dimm2 = dim-2;

  double gam, gam1, invgam1, pstiff;
  gam = vf->getGamma();
  pstiff = vf->getPressureConstant();
  gam1 = gam - 1.0; 
  invgam1 = 1.0 / gam1;

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;

  double Vb[dim];
  vf->conservativeToPrimitive(Ub, Vb);

  double Vvel2 = V[1]*V[1]+V[2]*V[2]+V[3]*V[3];
 
  double rho, u, v, w, p, nut, eps, k;
  
 // if (V[1]*n[0]+V[2]*n[1]+V[3]*n[2] == 0.0){
//    rho = Vb[0];
//    u = V[1];
//    v = V[2];
//    w = V[3];
//    p = V[4];
//  }

 if (V[1]*n[0]+V[2]*n[1]+V[3]*n[2] <= 0.0){
    rho = Vb[0];
    u = Vb[1];
    v = Vb[2];
    w = Vb[3];
    // If supersonic, all characteristics come from boundary
    if(Vvel2<1) p = V[4];
    else p = Vb[4];
    if (type == 1)
      nut = Vb[dimm1];
    else if (type == 2){
      k = Vb[dimm2];
      eps = Vb[dimm2];
    }
  }else{
    rho = V[0];
    u = V[1];
    v = V[2];
    w = V[3];
    // If supersonic, all characteristics come from interior
    if (Vvel2<1) p = Vb[4];
    else p = V[4];
    if (type == 1)
      nut = V[dimm1];
    else if (type == 2){
      k = V[dimm2];
      eps = V[dimm2];
    }

  }
  double rhoun = rho * (u*n[0] + v*n[1] + w*n[2] - nVel);
  double q = u*u + v*v + w*w;

  flux[0] = S * rhoun;
  flux[1] = S * (rhoun*u + p*n[0]);
  flux[2] = S * (rhoun*v + p*n[1]);
  flux[3] = S * (rhoun*w + p*n[2]);
  flux[4] = S * (rhoun*(gam*invgam1*(p+pstiff)/rho + 0.5*q) + p*nVel);

  if (type == 1) {
    flux[5] = S * rhoun*nut;
  }
  else if (type == 2) {
    flux[5] = S * rhoun*k;
    flux[6] = S * rhoun*eps;
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
inline
void influx3DDerivative(int type, VarFcnBase *vf, double* normal, double* dNormal,
	      double normalVel, double dNormalVel, double* V, double* Ub, double* dUb, double* flux, double* dFlux)
{

  const int dimm1 = dim-1;
  const int dimm2 = dim-2;

  double gam, gam1, invgam1, pstiff;
  gam = vf->getGamma();
  pstiff = vf->getPressureConstant();
  gam1 = gam - 1.0; 
  invgam1 = 1.0 / gam1;

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double dS = 1/S*(normal[0]*dNormal[0] + normal[1]*dNormal[1] + normal[2]*dNormal[2]);

  double ooS = 1.0 / S;
  double dooS = -1.0 / (S*S) * dS;

  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double dn[3] = {dNormal[0]*ooS + normal[0]*dooS, dNormal[1]*ooS + normal[1]*dooS, dNormal[2]*ooS + normal[2]*dooS};

  double nVel = normalVel * ooS;
  double dnVel = dNormalVel * ooS + normalVel * dooS;
  double Vvel2 = V[1]*V[1]+V[2]*V[2]+V[3]*V[3];

  double Vb[dim];
  vf->conservativeToPrimitive(Ub, Vb);

  double dVb[dim];
  vf->conservativeToPrimitiveDerivative(Ub, dUb, Vb, dVb);

  double rho, u, v, w, p, nut, eps, k;
  double drho, du, dv, dw, dp, dnut, deps, dk;
  
  double dV[dim];
  for (int i=1;i<dim;i++)
    dV[i] = 0.0;
  
//  if (V[1]*n[0]+V[2]*n[1]+V[3]*n[2] == 0.0){
//    rho = Vb[0];
//    u = V[1];
//    v = V[2];
//    w = V[3];
//    p = V[4];
//    drho = dVb[0];
//    du = dV[1];
//    dv = dV[2];
//    dw = dV[3];
//    dp = dV[4];
//  }

  if (V[1]*n[0]+V[2]*n[1]+V[3]*n[2] <= 0.0){
    rho = Vb[0];
    u = Vb[1];
    v = Vb[2];
    w = Vb[3];
    if(Vvel2 < 1) p = V[4];
    else p = Vb[4];
    drho = dVb[0];
    du = dVb[1];
    dv = dVb[2];
    dw = dVb[3];
    if(Vvel2 < 1) dp = dV[4];
    else dp = dVb[4];
    if (type == 1) {
      nut = Vb[dimm1];
      dnut = dVb[dimm1];
    }
    else if (type == 2){
      k = Vb[dimm2];
      eps = Vb[dimm2];
      dk = Vb[dimm2];
      deps = Vb[dimm2];
    }
  }else{
    rho = V[0];
    u = V[1];
    v = V[2];
    w = V[3];
    if(Vvel2 < 1) p = Vb[4];
    else p=V[4];
    drho = dV[0];
    du = dV[1];
    dv = dV[2];
    dw = dV[3];
    if(Vvel2 < 1) dp = dVb[4];
    else dp = dV[4];
    if (type == 1) {
      nut = V[dimm1];
      dnut = dV[dimm1];
    }
    else if (type == 2){
      k = V[dimm2];
      eps = V[dimm2];
      dk = dV[dimm2];
      deps = dV[dimm2];
    }

  }
  double rhoun = rho * (u*n[0] + v*n[1] + w*n[2] - nVel);
  double q = u*u + v*v + w*w;
  double drhoun = drho * (u*n[0] + v*n[1] + w*n[2] - nVel) + rho * (du*n[0] + u*dn[0] + dv*n[1] + v*dn[1] + dw*n[2] + w*dn[2] - dnVel);
  double dq = 2.0*(u*du + v*dv + w*dw);

  flux[0] = S * rhoun;
  flux[1] = S * (rhoun*u + p*n[0]);
  flux[2] = S * (rhoun*v + p*n[1]);
  flux[3] = S * (rhoun*w + p*n[2]);
  flux[4] = S * rhoun*(gam*invgam1*p/rho + 0.5*q);

  dFlux[0] = dS * rhoun + S * drhoun;
  dFlux[1] = dS * (rhoun*u + p*n[0]) + S * (drhoun*u + rhoun*du + dp*n[0] + p*dn[0]);
  dFlux[2] = dS * (rhoun*v + p*n[1]) + S * (drhoun*v + rhoun*dv + dp*n[1] + p*dn[1]);
  dFlux[3] = dS * (rhoun*w + p*n[2]) + S * (drhoun*w + rhoun*dw + dp*n[2] + p*dn[2]);
  dFlux[4] = dS * rhoun*(gam*invgam1*p/rho + 0.5*q) + S * drhoun*(gam*invgam1*p/rho + 0.5*q) + S * rhoun*(gam*invgam1*(dp*rho-p*drho)/(rho*rho) + 0.5*dq);

  if (type == 1) {
    flux[5] = S * rhoun*nut;
    dFlux[5] = dS * rhoun*nut + S * drhoun*nut + S * rhoun*dnut;
  }
  else if (type == 2) {
    flux[5] = S * rhoun*k;
    flux[6] = S * rhoun*eps;
    dFlux[5] = dS * rhoun*k + S * drhoun*k + S * rhoun*dk;
    dFlux[6] = dS * rhoun*eps + S * drhoun*eps + S * rhoun*deps;
  }

}

//------------------------------------------------------------------------------

template<int dim>
inline
void jacinflux3D(int type, VarFcnBase *vf, FluxFcnBase::Type localTypeJac, 
		 double* normal, double normalVel, double* V, double* Ub, double* jac)
{

  double _dfdV[dim][dim];
  double* dfdV = reinterpret_cast<double*>(_dfdV);
  int j;
  for (j=0; j<dim*dim; ++j)
    dfdV[j] = 0.0;

  const int dimm1 = dim-1;
  const int dimm2 = dim-2;

  double gam, gam1, invgam1, pstiff;
  gam = vf->getGamma();
  pstiff = vf->getPressureConstant();
  gam1 = gam - 1.0;
  invgam1 = 1.0 / gam1;

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;
  double Vvel2 = V[1]*V[1]+V[2]*V[2]+V[3]*V[3];

  double Vb[dim];
  vf->conservativeToPrimitive(Ub, Vb);
  double rho, u, v, w, p;

  if(V[1]*n[0]+V[2]*n[1]+V[3]*n[2]<= 0.0){
    if (Vvel2 > 1){ //supersonic
      // nothing to do: jacobian is null
    } else { // subsonic
      u = Vb[1];
      v = Vb[2];
      w = Vb[3];

      _dfdV[1][4] = S * n[0];
      _dfdV[2][4] = S * n[1];
      _dfdV[3][4] = S * n[2];
      _dfdV[4][4] = S * (gam*invgam1 * (u*n[0] + v*n[1] + w*n[2] - nVel) + nVel);
    }
  }else{
    rho = V[0];
    u = V[1];
    v = V[2];
    w = V[3];
    if(Vvel2<1) p = Vb[4];
    else p = V[4];
                                                                                                                                                          
    double un = u*n[0] + v*n[1] + w*n[2] - nVel;
    double q = u*u + v*v + w*w;
    double rhoH = 0.5*rho*q + gam*invgam1*(p+pstiff);
                                                                                                                                                                                                     
    _dfdV[0][0] = S * un;
    _dfdV[0][1] = S * rho*n[0];
    _dfdV[0][2] = S * rho*n[1];
    _dfdV[0][3] = S * rho*n[2];
    _dfdV[1][0] = S * u*un;
    _dfdV[1][1] = S * (rho*un + rho*u*n[0]);
    _dfdV[1][2] = S * rho*u*n[1];
    _dfdV[1][3] = S * rho*u*n[2];
    _dfdV[2][0] = S * v*un;
    _dfdV[2][1] = S * rho*v*n[0];
    _dfdV[2][2] = S * (rho*un + rho*v*n[1]);
    _dfdV[2][3] = S * rho*v*n[2];
    _dfdV[3][0] = S * w*un;
    _dfdV[3][1] = S * rho*w*n[0];
    _dfdV[3][2] = S * rho*w*n[1];
    _dfdV[3][3] = S * (rho*un + rho*w*n[2]);
    _dfdV[4][0] = S * 0.5*q*un;
    _dfdV[4][1] = S * (rho*u*un + rhoH*n[0]);
    _dfdV[4][2] = S * (rho*v*un + rhoH*n[1]);
    _dfdV[4][3] = S * (rho*w*un + rhoH*n[2]);
                                                                                                                                                                                                     
    if (type == 1) {
      double nut = V[dimm1];
      _dfdV[dimm1][0] = S * nut*un;
      _dfdV[dimm1][1] = S * rho*nut*n[0];
      _dfdV[dimm1][2] = S * rho*nut*n[1];
      _dfdV[dimm1][3] = S * rho*nut*n[2];
      _dfdV[dimm1][dimm1] = S * rho*un;
    }
    else if (type == 2) {
      double k = V[dimm2];
      double eps = V[dimm1];
      _dfdV[dimm2][0] = S * k*un;
      _dfdV[dimm2][1] = S * rho*k*n[0];
      _dfdV[dimm2][2] = S * rho*k*n[1];
      _dfdV[dimm2][3] = S * rho*k*n[2];
      _dfdV[dimm2][dimm2] = S * rho*un;
      _dfdV[dimm1][0] = S * eps*un;
      _dfdV[dimm1][1] = S * rho*eps*n[0];
      _dfdV[dimm1][2] = S * rho*eps*n[1];
      _dfdV[dimm1][3] = S * rho*eps*n[2];
      _dfdV[dimm1][dimm1] = S * rho*un;
    }
  }


  if (localTypeJac == FluxFcnBase::CONSERVATIVE)
    vf->postMultiplyBydVdU(V, dfdV, jac);
  else
    for (j=0; j<dim*dim; ++j)
      jac[j] = dfdV[j]; 

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

template<int dim>
inline
void directstateflux3D(int type, VarFcnBase *vf, double* normal, 
	      double normalVel, double* V, double* Ub, double* flux)
{
  const int dimm1 = dim-1;
  const int dimm2 = dim-2;

  double gam, gam1, invgam1, pstiff;
  gam = vf->getGamma();
  pstiff = vf->getPressureConstant();
  gam1 = gam - 1.0; 
  invgam1 = 1.0 / gam1;

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;

  double Vb[dim];
//  vf->conservativeToPrimitive(Ub, Vb);

  double Vvel2 = V[1]*V[1]+V[2]*V[2]+V[3]*V[3];
 
  double rho, u, v, w, p, nut, eps, k;

  double un = V[1]*n[0]+V[2]*n[1]+V[3]*n[2];

  if (un <= 0.0){
    double ci = vf->computeSoundSpeed(V);
    double Rplus = un + 2.*ci*invgam1;
    double Ht = invgam1*ci*ci + 0.5*Vvel2;

    double a = 1 + 2.*invgam1;
    double b = -2*Rplus;
    double c = 0.5*gam1*(Rplus*Rplus - 2*Ht);

    double term1 = -b/(2*a);
    double term2 = sqrt(b*b - 4*a*c)/(2*a);

    double cb = std::max(term1 + term2, term1 - term2);
    double U = Rplus - 2*cb*invgam1;
    double Mb = fabs(U)/cb;
    double pb = Ub[4]/pow(1.0 + 0.5*gam1*Mb*Mb,gam*invgam1);
    double Tb = Ub[0]/(1.0 + 0.5*gam1*Mb*Mb);

  
    rho = invgam1*(pb+pstiff)/Tb;
    u = U*n[0];
    v = U*n[1];
    w = U*n[2];
    p = pb;

//  fprintf(stdout,"U,cbc,cb,Mb -> %f,%f,%f,%f\n",U,sqrt(gam*p/rho),cb,Mb);
    if (type == 1)
      nut = Ub[dimm1];
    else if (type == 2){
      k = Ub[dimm2];
      eps = Ub[dimm2];
    }
  }
  else{
    rho = V[0];
    u = V[1];
    v = V[2];
    w = V[3];
    // If supersonic, all characteristics come from interior
    if (vf->computeMachNumber(V)<1.0) p = Ub[4];
    else p = V[4];
    if (type == 1)
      nut = V[dimm1];
    else if (type == 2){
      k = V[dimm2];
      eps = V[dimm2];
    }
  }

  double rhoun = rho * (u*n[0] + v*n[1] + w*n[2] - nVel);
  double q = u*u + v*v + w*w;

  flux[0] = S * rhoun;
  flux[1] = S * (rhoun*u + p*n[0]);
  flux[2] = S * (rhoun*v + p*n[1]);
  flux[3] = S * (rhoun*w + p*n[2]);
  flux[4] = S * (rhoun*(gam*invgam1*(p+pstiff)/rho + 0.5*q) + p*nVel);

  if (type == 1) {
    flux[5] = S * rhoun*nut;
  }
  else if (type == 2) {
    flux[5] = S * rhoun*k;
    flux[6] = S * rhoun*eps;
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
inline
void directstateflux3DDerivative(int type, VarFcnBase *vf, double* normal, double* dNormal,
	      double normalVel, double dNormalVel, double* V, double* Ub, double* dUb, double* flux, double* dFlux)
{
// not implemented yet
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

template<int dim>
inline
void massflowflux3D(int type, VarFcnBase *vf, double* normal, 
	      double normalVel, double* V, double* Ub, double* flux)
{
  const int dimm1 = dim-1;
  const int dimm2 = dim-2;

  double gam, gam1, invgam1, pstiff;
  gam = vf->getGamma();
  pstiff = vf->getPressureConstant();
  gam1 = gam - 1.0; 
  invgam1 = 1.0 / gam1;

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;

  double Vb[dim];
//  vf->conservativeToPrimitive(Ub, Vb);

  double Vvel2 = V[1]*V[1]+V[2]*V[2]+V[3]*V[3];
 
  double rho, u, v, w, p, nut, eps, k;

  double un = V[1]*n[0]+V[2]*n[1]+V[3]*n[2];

  if (un <= 0.0){
    double ci = vf->computeSoundSpeed(V);
    double Rplus = un + 2.*ci*invgam1;
    double Ht = invgam1*ci*ci + 0.5*Vvel2;

    double a = 1 + 2.*invgam1;
    double b = -2*Rplus;
    double c = 0.5*gam1*(Rplus*Rplus - 2*Ht);

    double term1 = -b/(2*a);
    double term2 = sqrt(b*b - 4*a*c)/(2*a);

    double cb = std::max(term1 + term2, term1 - term2);
    double U = Rplus - 2*cb*invgam1;
    double Mb = fabs(U)/cb;
    double pb = Ub[4]/pow(1.0 + 0.5*gam1*Mb*Mb,gam*invgam1);
    double Tb = Ub[0]/(1.0 + 0.5*gam1*Mb*Mb);

  
    rho = invgam1*(pb+pstiff)/Tb;
    u = U*n[0];
    v = U*n[1];
    w = U*n[2];
    p = pb;

//  fprintf(stdout,"U,cbc,cb,Mb -> %f,%f,%f,%f\n",U,sqrt(gam*p/rho),cb,Mb);
    if (type == 1)
      nut = Ub[dimm1];
    else if (type == 2){
      k = Ub[dimm2];
      eps = Ub[dimm2];
    }
  }
  else{
    rho = V[0];
    u = V[1];
    v = V[2];
    w = V[3];
    // If supersonic, all characteristics come from interior
    if (vf->computeMachNumber(V)<1.0) p = Ub[4];
    else p = V[4];
    if (type == 1)
      nut = V[dimm1];
    else if (type == 2){
      k = V[dimm2];
      eps = V[dimm2];
    }
  }

  double rhoun = rho * (u*n[0] + v*n[1] + w*n[2] - nVel);
  double q = u*u + v*v + w*w;

  flux[0] = S * rhoun;
  flux[1] = S * (rhoun*u + p*n[0]);
  flux[2] = S * (rhoun*v + p*n[1]);
  flux[3] = S * (rhoun*w + p*n[2]);
  flux[4] = S * (rhoun*(gam*invgam1*(p+pstiff)/rho + 0.5*q) + p*nVel);

  if (type == 1) {
    flux[5] = S * rhoun*nut;
  }
  else if (type == 2) {
    flux[5] = S * rhoun*k;
    flux[6] = S * rhoun*eps;
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
inline
void massflowflux3DDerivative(int type, VarFcnBase *vf, double* normal, double* dNormal,
	      double normalVel, double dNormalVel, double* V, double* Ub, double* dUb, double* flux, double* dFlux)
{
// not implemented yet
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


void FluxFcnSGInternalInflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                           double *V, double *Ub, double *flux, bool useLimiter)
{

  influx3D<5>(0, vf, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGInternalInflowEuler3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  influx3DDerivative<5>(0, vf, normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------

void FluxFcnSGInternalInflowEuler3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                   double *V, double *Ub, double *jacL, bool useLimiter)
{

  jacinflux3D<5>(0, vf, typeJac, normal, normalVel, V, Ub, jacL);

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//
void FluxFcnSGDirectStateInflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                           double *V, double *Ub, double *flux, bool useLimiter)
{

  directstateflux3D<5>(0, vf, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGDirectStateInflowEuler3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  directstateflux3DDerivative<5>(0, vf, normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------

void FluxFcnSGMassFlowInflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                           double *V, double *Ub, double *flux, bool useLimiter)
{

  massflowflux3D<5>(0, vf, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGMassFlowInflowEuler3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  massflowflux3DDerivative<5>(0, vf, normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------
void FluxFcnSGInternalOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                            double *V, double *Ub, double *flux, bool useLimiter)
{

  influx3D<5>(0, vf, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGInternalOutflowEuler3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  influx3DDerivative<5>(0, vf, normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------

void FluxFcnSGInternalOutflowEuler3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                    double *V, double *Ub, double *jacL, bool useLimiter)
{

  jacinflux3D<5>(0, vf, typeJac, normal, normalVel, V, Ub, jacL);

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void FluxFcnSGDirectStateOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                            double *V, double *Ub, double *flux, bool useLimiter)
{

  directstateflux3D<5>(0, vf, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGDirectStateOutflowEuler3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  directstateflux3DDerivative<5>(0, vf, normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------

void FluxFcnSGMassFlowOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel,
                                            double *V, double *Ub, double *flux, bool useLimiter)
{

  massflowflux3D<5>(0, vf, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGMassFlowOutflowEuler3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  massflowflux3DDerivative<5>(0, vf, normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------
//turbulence

void FluxFcnSGFDJacRoeSA3D::compute(double length, double irey, double *normal, double normalVel,
                                    double *VL, double *VR, double *flux, bool useLimiter)
{
  int prec = useLimiter ? sprec.getPrecTag() : 0;
  if(prec) { // roe flux with spatial low mach preconditioner and numerical viscosity
    F77NAME(roeflux5)(1, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux,
                      sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, prec);
  }
  else if(gamma != 1) { // original implementation of the roe flux with numerical viscosity
    F77NAME(roeflux6)(1, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux);
  }
  else { // faster implementation of the roe flux
    F77NAME(roeflux7)(1, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
  }
}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGFDJacRoeSA3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, bool useLimiter)
{

  F77NAME(gxroeflux5)(1, gamma, vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VL, dVL, VR, dVR, VR, dVR, flux, dFlux, sprec.getMinMach(), 0.0*dMach, sprec.getSlope(), sprec.getCutOffMach(), irey, dIrey, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGApprJacRoeSA3D::compute(double length, double irey, double *normal, double normalVel,
                                      double *VL, double *VR, double *flux, bool useLimiter)
{
  int prec = useLimiter ? sprec.getPrecTag() : 0;
  if(prec) { // roe flux with spatial low mach preconditioner, numerical viscosity and galerkin advective operator
    F77NAME(roeflux5)(1, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux,
                      sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, prec);
  }
  else if(gamma != 1 || rshift) { // original implementation of the roe flux with numerical viscosity and galerkin advection operator
    F77NAME(roeflux6)(1, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux);
  }
  else { // faster implementation of the roe flux
    F77NAME(roeflux7)(1, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
  }
}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGApprJacRoeSA3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, bool useLimiter)
{

  F77NAME(gxroeflux5)(1, gamma, vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VL+rshift, dVL+rshift, VR, dVR, VR+rshift, dVR+rshift, flux, dFlux, sprec.getMinMach(), 0.0*dMach, sprec.getSlope(), sprec.getCutOffMach(), irey, dIrey, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGApprJacRoeSA3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                               double *VL, double *VR, double *jacL, double *jacR, bool useLimiter)
{
  roejacappr3Dgas<6>(1, gamma, vf, vf->getGamma(), vf->getPressureConstant(), typeJac, normal, normalVel, VL, VR, jacL, jacR,
                     irey, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), length,
                     useLimiter ? sprec.getPrecTag() : 0);
}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnSGExactJacRoeSA3D::compute(double length, double irey, double *normal, double normalVel,
                                       double *VL, double *VR, double *flux, bool useLimiter)
{
  int prec = useLimiter ? sprec.getPrecTag() : 0;
  if(prec) { // roe flux with spatial low mach preconditioner and numerical viscosity
    F77NAME(roeflux5)(1, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux,
                      sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, prec);
  }
  else if(gamma != 1) { // original implementation of the roe flux with numerical viscosity
    F77NAME(roeflux6)(1, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux);
  }
  else { // faster implementation of the roe flux
    F77NAME(roeflux7)(1, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
  }
}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGExactJacRoeSA3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, bool useLimiter)
{
  F77NAME(gxroeflux6)(1, gamma, vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VL, dVL, VR, dVR, VR, dVR, flux, dFlux);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnSGExactJacRoeSA3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                                double *VL, double *VR, double *jacL, double *jacR, bool useLimiter)
{
  roejacexact3D<6>(1, gamma, vf, typeJac, normal, normalVel, VL, VR, jacL, jacR, irey, sprec.getMinMach(), sprec.getSlope(),
                   sprec.getCutOffMach(), sprec.getShockParameter(), length, useLimiter ? sprec.getPrecTag() : 0);
}

//------------------------------------------------------------------------------

void FluxFcnSGFDJacHLLESA3D::compute(double length, double irey, double *normal, double normalVel,
                                  double *VL, double *VR, double *flux, bool useLimiter)
{

  F77NAME(hlleflux)(1, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGApprJacHLLESA3D::compute(double length, double irey, double *normal, double normalVel,
                                    double *VL, double *VR, double *flux, bool useLimiter)
{

  F77NAME(hlleflux)(1, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGApprJacHLLESA3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                             double *VL, double *VR,
                                             double *jacL, double *jacR, bool useLimiter)
{

  hllejacappr3Dgas<6>(1, gamma, vf, vf->getGamma(), vf->getPressureConstant(), typeJac, normal, normalVel, VL, VR, jacL, jacR, irey, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), useLimiter ? sprec.getPrecTag() : 0, sprec.getShockParameter(), length);

}

//------------------------------------------------------------------------------

void FluxFcnSGFDJacHLLCSA3D::compute(double length, double irey, double *normal, double normalVel,
                                  double *VL, double *VR, double *flux, bool useLimiter)
{

  F77NAME(hllcflux)(1, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGApprJacHLLCSA3D::compute(double length, double irey, double *normal, double normalVel,
                                    double *VL, double *VR, double *flux, bool useLimiter)
{

  F77NAME(hllcflux)(1, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGApprJacHLLCSA3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                             double *VL, double *VR,
                                             double *jacL, double *jacR, bool useLimiter)
{

  hllcjacappr3Dgas<6>(1, gamma, vf, vf->getGamma(), vf->getPressureConstant(), typeJac, normal, normalVel, VL, VR, jacL, jacR, irey, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGWallSA3D::compute(double length, double irey, double *normal, double normalVel,
                              double *V, double *Ub, double *flux, bool useLimiter)
{

  flux[0] = 0.0;
  flux[1] = V[4] * normal[0];
  flux[2] = V[4] * normal[1];
  flux[3] = V[4] * normal[2];
  flux[4] = V[4] * normalVel;
  flux[5] = 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGWallSA3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  dFlux[0] = 0.0;
  dFlux[1] = V[4] * dNormal[0];
  dFlux[2] = V[4] * dNormal[1];
  dFlux[3] = V[4] * dNormal[2];
  dFlux[4] = V[4] * dNormalVel;
  dFlux[5] = 0.0;

}

//------------------------------------------------------------------------------

// Included (MB*)
void FluxFcnSGWallSA3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                   double *V, double *Ub, double *jac, bool useLimiter)
{

//  flux[0] = 0.0;
//  flux[1] = V[4] * normal[0];
//  flux[2] = V[4] * normal[1];
//  flux[3] = V[4] * normal[2];
//  flux[4] = V[4] * normalVel;
//  flux[5] = 0.0;

  double dfdV[36];

  dfdV[0] = 0.0;
  dfdV[1] = 0.0;
  dfdV[2] = 0.0;
  dfdV[3] = 0.0;
  dfdV[4] = 0.0;
  dfdV[5] = 0.0;

  dfdV[6] = 0.0;
  dfdV[7] = 0.0;
  dfdV[8] = 0.0;
  dfdV[9] = 0.0;
  dfdV[10] = normal[0];
  dfdV[11] = 0.0;

  dfdV[12] = 0.0;
  dfdV[13] = 0.0;
  dfdV[14] = 0.0;
  dfdV[15] = 0.0;
  dfdV[16] = normal[1];
  dfdV[17] = 0.0;

  dfdV[18] = 0.0;
  dfdV[19] = 0.0;
  dfdV[20] = 0.0;
  dfdV[21] = 0.0;
  dfdV[22] = normal[2];
  dfdV[23] = 0.0;

  dfdV[24] = 0.0;
  dfdV[25] = 0.0;
  dfdV[26] = 0.0;
  dfdV[27] = 0.0;
  dfdV[28] = normalVel;
  dfdV[29] = 0.0;

  dfdV[30] = 0.0;
  dfdV[31] = 0.0;
  dfdV[32] = 0.0;
  dfdV[33] = 0.0;
  dfdV[34] = 0.0;
  dfdV[35] = 0.0;

  if (typeJac == FluxFcnBase::CONSERVATIVE)
    vf->postMultiplyBydVdU(V, dfdV, jac);
  else
    for (int k=0; k<36; ++k)
      jac[k] = dfdV[k];

}

//------------------------------------------------------------------------------

void FluxFcnSGPorousWallSA3D::compute(double length, double irey, double *normal, double normalVel,
                              double *V, double *Ub, double *flux, bool useLimiter)
{

  double gam, gam1, invgam1, pstiff;
  gam = vf->getGamma();
  pstiff = vf->getPressureConstant();
  gam1 = gam - 1.0; 
  invgam1 = 1.0 / gam1;

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;
  double porosity = Ub[1];

  if (porosity != 0 ) {
    double rho,rhoun,u,v,w,p,q;
    if ( Ub[0] <= 0.0 ) {
      rho = invgam1 * (V[4]+pstiff) / Ub[4];
      double ujet = Ub[0]/rho;
      double un = ujet/porosity;
      rhoun = rho*(ujet - nVel);
      u = un*n[0];
      v = un*n[1];
      w = un*n[2];
      q = u*u + v*v + w*w;
      p = V[4];
    }
    else if ( Ub[0] > 0.0 ) {
      rho = V[0];
      double ujet = Ub[0]/rho;
      double Un = V[1]*n[0] + V[2]*n[1] + V[3]*n[2];
      double Vt[3] = {V[1] - Un*n[0],V[2] - Un*n[1],V[3] - Un*n[2]};
      rhoun = rho*(ujet - nVel);
      u = ujet*n[0] + Vt[0];
      v = ujet*n[1] + Vt[1];
      w = ujet*n[2] + Vt[2];
      double u1 = ujet*n[0]/porosity + Vt[0];
      double v1 = ujet*n[1]/porosity + Vt[1];
      double w1 = ujet*n[2]/porosity + Vt[2];
      q = u1*u1 + v1*v1 + w1*w1;
      p = V[4];
    } 
    flux[0] = S * rhoun;
    flux[1] = S * (rhoun*u + p*n[0]);
    flux[2] = S * (rhoun*v + p*n[1]);
    flux[3] = S * (rhoun*w + p*n[2]);
    flux[4] = S * (rhoun*(gam*invgam1*(p+pstiff)/rho + 0.5*q) + p*nVel);
  }
  else {
    flux[0] = 0.0;
    flux[1] = V[4] * normal[0];
    flux[2] = V[4] * normal[1];
    flux[3] = V[4] * normal[2];
    flux[4] = V[4] * normalVel;
  }
  flux[5] = 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGPorousWallSA3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{
// not implemented
}

//------------------------------------------------------------------------------

void FluxFcnSGGhidagliaSA3D::compute(double length, double irey, double *normal, double normalVel,
                                   double *V, double *Ub, double *flux, bool useLimiter)
{

  F77NAME(genbcfluxgas)(1, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

void FluxFcnSGOutflowSA3D::compute(double length, double irey, double *normal, double normalVel,
                                 double *V, double *Ub, double *flux, bool useLimiter)
{

  F77NAME(boundflux5)(1, vf->getGamma(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGOutflowSA3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  F77NAME(gxboundflux5)(1, vf->getGamma(), normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------

void FluxFcnSGInternalInflowSA3D::compute(double length, double irey, double *normal, double normalVel,
                                        double *V, double *Ub, double *flux, bool useLimiter)
{

  influx3D<6>(1, vf, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGInternalInflowSA3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  influx3DDerivative<6>(1, vf, normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------

void FluxFcnSGInternalInflowSA3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                double *V, double *Ub, double *jacL, bool useLimiter)
{

  jacinflux3D<6>(1, vf, typeJac, normal, normalVel, V, Ub, jacL);

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//
void FluxFcnSGDirectStateInflowSA3D::compute(double length, double irey, double *normal, double normalVel,
                                        double *V, double *Ub, double *flux, bool useLimiter)
{

  directstateflux3D<6>(1, vf, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGDirectStateInflowSA3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  directstateflux3DDerivative<6>(1, vf, normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------

void FluxFcnSGMassFlowInflowSA3D::compute(double length, double irey, double *normal, double normalVel,
                                        double *V, double *Ub, double *flux, bool useLimiter)
{

  massflowflux3D<6>(1, vf, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGMassFlowInflowSA3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  massflowflux3DDerivative<6>(1, vf, normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void FluxFcnSGInternalOutflowSA3D::compute(double length, double irey, double *normal, double normalVel,
                                         double *V, double *Ub, double *flux, bool useLimiter)
{

  influx3D<6>(1, vf, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGInternalOutflowSA3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  influx3DDerivative<6>(1, vf, normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------

void FluxFcnSGInternalOutflowSA3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                 double *V, double *Ub, double *jacL, bool useLimiter)
{

  jacinflux3D<6>(1, vf, typeJac, normal, normalVel, V, Ub, jacL);

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void FluxFcnSGDirectStateOutflowSA3D::compute(double length, double irey, double *normal, double normalVel,
                                         double *V, double *Ub, double *flux, bool useLimiter)
{

  directstateflux3D<6>(1, vf, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGDirectStateOutflowSA3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  directstateflux3DDerivative<6>(1, vf, normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------

void FluxFcnSGMassFlowOutflowSA3D::compute(double length, double irey, double *normal, double normalVel,
                                         double *V, double *Ub, double *flux, bool useLimiter)
{

  massflowflux3D<6>(1, vf, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGMassFlowOutflowSA3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  massflowflux3DDerivative<6>(1, vf, normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------
// note: jacL = dFdUL and jacR = dFdUR
                                                                                                                  
void FluxFcnSGRoeSAturb3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                          double *VL, double *VR,
                                          double *jacL, double *jacR, bool useLimiter)
{

  F77NAME(roejac2)(1, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR, sprec.getMinMach(),sprec.getSlope(),sprec.getCutOffMach(),irey, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGWallSAturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                          double *V, double *Ub, double *jacL, bool useLimiter)
{

  jacL[0] = 0.0;

}

//------------------------------------------------------------------------------

void FluxFcnSGPorousWallSAturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                          double *V, double *Ub, double *jacL, bool useLimiter)
{

  jacL[0] = 0.0;

}

//------------------------------------------------------------------------------

void FluxFcnSGGhidagliaSAturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                             double *V, double *Ub, double *jacL, bool useLimiter)
{

  double flux[6];
  F77NAME(genbcfluxgas)(1, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);

  double updir = 1.0;
  if(flux[0]<0.0)  updir = 0.0;
  else if(flux[0]==0.0) updir = 0.5;

  jacL[0] = flux[0] * updir / V[0];

}

//------------------------------------------------------------------------------
// note: jacL = dFdUL

void FluxFcnSGOutflowSAturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                             double *V, double *Ub, double *jacL, bool useLimiter)
{

  F77NAME(boundjac2)(1, vf->getGamma(), normal, normalVel, V, Ub, jacL);

}

//------------------------------------------------------------------------------

void FluxFcnSGInternalInflowSAturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                    double *V, double *Ub, double *jacL, bool useLimiter)
{

  jacL[0] = 0.0;

}

//------------------------------------------------------------------------------

void FluxFcnSGInternalOutflowSAturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                     double *V, double *Ub, double *jacL, bool useLimiter)
{

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;

  double rho = V[0];
  double u = V[1];
  double v = V[2];
  double w = V[3];
  double un = u*n[0] + v*n[1] + w*n[2] - nVel;

  if (typeJac == FluxFcnBase::CONSERVATIVE)
    jacL[0] = S * un;
  else
    jacL[0] = S * rho*un;

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//
void FluxFcnSGDirectStateInflowSAturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                    double *V, double *Ub, double *jacL, bool useLimiter)
{

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;

  double rho = Ub[0];
  double u = Ub[1]/Ub[0];
  double v = Ub[2]/Ub[0];
  double w = Ub[3]/Ub[0];
  double un = u*n[0] + v*n[1] + w*n[2] - nVel;

  if (typeJac == FluxFcnBase::CONSERVATIVE)
    jacL[0] = S * un;
  else
    jacL[0] = S * rho*un;

}

//------------------------------------------------------------------------------

void FluxFcnSGMassFlowInflowSAturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                    double *V, double *Ub, double *jacL, bool useLimiter)
{

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;

  double rho = Ub[0];
  double u = Ub[1]/Ub[0];
  double v = Ub[2]/Ub[0];
  double w = Ub[3]/Ub[0];
  double un = u*n[0] + v*n[1] + w*n[2] - nVel;

  if (typeJac == FluxFcnBase::CONSERVATIVE)
    jacL[0] = S * un;
  else
    jacL[0] = S * rho*un;

}

//------------------------------------------------------------------------------

void FluxFcnSGDirectStateOutflowSAturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                     double *V, double *Ub, double *jacL, bool useLimiter)
{

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;

  double rho = Ub[0];
  double u = Ub[1]/Ub[0];
  double v = Ub[2]/Ub[0];
  double w = Ub[3]/Ub[0];
  double un = u*n[0] + v*n[1] + w*n[2] - nVel;

  if (typeJac == FluxFcnBase::CONSERVATIVE)
    jacL[0] = S * un;
  else
    jacL[0] = S * rho*un;

}

//------------------------------------------------------------------------------

void FluxFcnSGMassFlowOutflowSAturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                     double *V, double *Ub, double *jacL, bool useLimiter)
{

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;

  double rho = Ub[0];
  double u = Ub[1]/Ub[0];
  double v = Ub[2]/Ub[0];
  double w = Ub[3]/Ub[0];
  double un = u*n[0] + v*n[1] + w*n[2] - nVel;

  if (typeJac == FluxFcnBase::CONSERVATIVE)
    jacL[0] = S * un;
  else
    jacL[0] = S * rho*un;

}

//------------------------------------------------------------------------------

void FluxFcnSGFDJacRoeKE3D::compute(double length, double irey, double *normal, double normalVel,
                                    double *VL, double *VR, double *flux, bool useLimiter)
{
  int prec = useLimiter ? sprec.getPrecTag() : 0;
  if(prec) { // roe flux with spatial low mach preconditioner with numerical viscosity
    F77NAME(roeflux5)(2, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux,
                      sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, prec);
  }
  else if(gamma != 1) { // original implementation of the roe flux with numerical viscosity
    F77NAME(roeflux6)(2, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux);
  }
  else { // faster implementation of the roe flux
    F77NAME(roeflux7)(2, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
  }
}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGFDJacRoeKE3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, bool useLimiter)
{
  F77NAME(gxroeflux5)(2, gamma, vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VL, dVL, VR, dVR, VR, dVR, flux, dFlux, sprec.getMinMach(), 0.0*dMach, sprec.getSlope(), sprec.getCutOffMach(), irey, dIrey, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGApprJacRoeKE3D::compute(double length, double irey, double *normal, double normalVel,
                                      double *VL, double *VR, double *flux, bool useLimiter)
{
  int prec = useLimiter ? sprec.getPrecTag() : 0;
  if(prec) { // roe flux with spatial low mach preconditioner, numerical viscosity and galerkin advective operator
    F77NAME(roeflux5)(2, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux,
                      sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, prec);
  }
  else if(gamma != 1 || rshift) { // original implementation of the roe flux with numerical viscosity and galerkin advection operator
    F77NAME(roeflux6)(2, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux);
  }
  else { // faster implementation of the roe flux
    F77NAME(roeflux7)(2, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
  }
}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGApprJacRoeKE3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, bool useLimiter)
{

  F77NAME(gxroeflux5)(2, gamma, vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VL+rshift, dVL+rshift, VR, dVR, VR+rshift, dVR+rshift, flux, dFlux, sprec.getMinMach(), 0.0*dMach, sprec.getSlope(), sprec.getCutOffMach(), irey, dIrey, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGApprJacRoeKE3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                               double *VL, double *VR, double *jacL, double *jacR, bool useLimiter)
{
  roejacappr3Dgas<7>(2, gamma, vf, vf->getGamma(), vf->getPressureConstant(), typeJac, normal, normalVel, VL, VR, jacL, jacR,
                     irey, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), length,
                     useLimiter ? sprec.getPrecTag() : 0);
}

//------------------------------------------------------------------------------

void FluxFcnSGExactJacRoeKE3D::compute(double length, double irey, double *normal, double normalVel,
                                       double *VL, double *VR, double *flux, bool useLimiter)
{
  int prec = useLimiter ? sprec.getPrecTag() : 0;
  if(prec) { // roe flux with spatial low mach preconditioner and numerical viscosity
    F77NAME(roeflux5)(2, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux,
                      sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, prec);
  }
  else if(gamma != 1) { // original implementation of the roe flux with numerical viscosity
    F77NAME(roeflux6)(2, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux);
  }
  else { // faster implementation of the roe flux
    F77NAME(roeflux7)(2, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, flux);
  }
}

//------------------------------------------------------------------------------

// Included (MB)
inline
void FluxFcnSGExactJacRoeKE3D::computeDerivative(double irey, double dIrey, double *normal, double *dNormal, double normalVel, double dNormalVel,
                                          double *VL, double *dVL, double *VR, double *dVR, double dMach, double *flux, double *dFlux, bool useLimiter)
{

  F77NAME(gxroeflux6)(2, gamma, vf->getGamma(), vf->getPressureConstant(), vf->getDerivativeOfPressureConstant(), normal, dNormal, normalVel, dNormalVel, VL, dVL, VL, dVL, VR, dVR, VR, dVR, flux, dFlux);

}

//------------------------------------------------------------------------------

void FluxFcnSGExactJacRoeKE3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                                double *VL, double *VR, double *jacL, double *jacR, bool useLimiter)
{
  roejacexact3D<7>(2, gamma, vf, typeJac, normal, normalVel, VL, VR, jacL, jacR, irey, sprec.getMinMach(), sprec.getSlope(),
                   sprec.getCutOffMach(), sprec.getShockParameter(), length, useLimiter ? sprec.getPrecTag() : 0);
}

//------------------------------------------------------------------------------

void FluxFcnSGFDJacHLLEKE3D::compute(double length, double irey, double *normal, double normalVel,
                                  double *VL, double *VR, double *flux, bool useLimiter)
{

  F77NAME(hlleflux)(2, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGApprJacHLLEKE3D::compute(double length, double irey, double *normal, double normalVel,
                                    double *VL, double *VR, double *flux, bool useLimiter)
{

  F77NAME(hlleflux)(2, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGApprJacHLLEKE3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                             double *VL, double *VR,
                                             double *jacL, double *jacR, bool useLimiter)
{

  hllejacappr3Dgas<7>(2, gamma, vf, vf->getGamma(), vf->getPressureConstant(), typeJac, normal, normalVel, VL, VR, jacL, jacR, irey, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), useLimiter ? sprec.getPrecTag() : 0, sprec.getShockParameter(), length);

}

//------------------------------------------------------------------------------

void FluxFcnSGFDJacHLLCKE3D::compute(double length, double irey, double *normal, double normalVel,
                                  double *VL, double *VR, double *flux, bool useLimiter)
{

  F77NAME(hllcflux)(2, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL, VR, VR, flux, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGApprJacHLLCKE3D::compute(double length, double irey, double *normal, double normalVel,
                                    double *VL, double *VR, double *flux, bool useLimiter)
{

  F77NAME(hllcflux)(2, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), irey, length, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGApprJacHLLCKE3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                             double *VL, double *VR,
                                             double *jacL, double *jacR, bool useLimiter)
{

  hllcjacappr3Dgas<7>(2, gamma, vf, vf->getGamma(), vf->getPressureConstant(), typeJac, normal, normalVel, VL, VR, jacL, jacR, irey, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGWallKE3D::compute(double length, double irey, double *normal, double normalVel,
                              double *V, double *Ub, double *flux, bool useLimiter)
{

  flux[0] = 0.0;
  flux[1] = V[4] * normal[0];
  flux[2] = V[4] * normal[1];
  flux[3] = V[4] * normal[2];
  flux[4] = V[4] * normalVel;
  flux[5] = 0.0;
  flux[6] = 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGWallKE3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  dFlux[0] = 0.0;
  dFlux[1] = V[4] * dNormal[0];
  dFlux[2] = V[4] * dNormal[1];
  dFlux[3] = V[4] * dNormal[2];
  dFlux[4] = V[4] * dNormalVel;
  dFlux[5] = 0.0;
  dFlux[6] = 0.0;

}

//------------------------------------------------------------------------------

// Included (MB*)
void FluxFcnSGWallKE3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                   double *V, double *Ub, double *jac, bool useLimiter)
{

//  flux[0] = 0.0;
//  flux[1] = V[4] * normal[0];
//  flux[2] = V[4] * normal[1];
//  flux[3] = V[4] * normal[2];
//  flux[4] = V[4] * normalVel;
//  flux[5] = 0.0;
//  flux[6] = 0.0;

  double dfdV[49];

  dfdV[0] = 0.0;
  dfdV[1] = 0.0;
  dfdV[2] = 0.0;
  dfdV[3] = 0.0;
  dfdV[4] = 0.0;
  dfdV[5] = 0.0;
  dfdV[6] = 0.0;

  dfdV[7] = 0.0;
  dfdV[8] = 0.0;
  dfdV[9] = 0.0;
  dfdV[10] = 0.0;
  dfdV[11] = normal[0];
  dfdV[12] = 0.0;
  dfdV[13] = 0.0;

  dfdV[14] = 0.0;
  dfdV[15] = 0.0;
  dfdV[16] = 0.0;
  dfdV[17] = 0.0;
  dfdV[18] = normal[1];
  dfdV[19] = 0.0;
  dfdV[20] = 0.0;

  dfdV[21] = 0.0;
  dfdV[22] = 0.0;
  dfdV[23] = 0.0;
  dfdV[24] = 0.0;
  dfdV[25] = normal[2];
  dfdV[26] = 0.0;
  dfdV[27] = 0.0;

  dfdV[28] = 0.0;
  dfdV[29] = 0.0;
  dfdV[30] = 0.0;
  dfdV[31] = 0.0;
  dfdV[32] = normalVel;
  dfdV[33] = 0.0;
  dfdV[34] = 0.0;

  dfdV[35] = 0.0;
  dfdV[36] = 0.0;
  dfdV[37] = 0.0;
  dfdV[38] = 0.0;
  dfdV[39] = 0.0;
  dfdV[40] = 0.0;
  dfdV[41] = 0.0;

  dfdV[42] = 0.0;
  dfdV[43] = 0.0;
  dfdV[44] = 0.0;
  dfdV[45] = 0.0;
  dfdV[46] = 0.0;
  dfdV[47] = 0.0;
  dfdV[48] = 0.0;
  
  if (typeJac == FluxFcnBase::CONSERVATIVE)
    vf->postMultiplyBydVdU(V, dfdV, jac);
  else
    for (int k=0; k<49; ++k)
      jac[k] = dfdV[k];

}

//------------------------------------------------------------------------------

void FluxFcnSGPorousWallKE3D::compute(double length, double irey, double *normal, double normalVel,
                              double *V, double *Ub, double *flux, bool useLimiter)
{

  double gam, gam1, invgam1, pstiff;
  gam = vf->getGamma();
  pstiff = vf->getPressureConstant();
  gam1 = gam - 1.0; 
  invgam1 = 1.0 / gam1;

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;
  double porosity = Ub[1];

  if (porosity != 0 ) {
    double rho,rhoun,u,v,w,p,q;
    if ( Ub[0] <= 0.0 ) {
      rho = invgam1 * (V[4]+pstiff) / Ub[4];
      double ujet = Ub[0]/rho;
      double un = ujet/porosity;
      rhoun = rho*(ujet - nVel);
      u = un*n[0];
      v = un*n[1];
      w = un*n[2];
      q = u*u + v*v + w*w;
      p = V[4];
    }
    else if ( Ub[0] > 0.0 ) {
      rho = V[0];
      double ujet = Ub[0]/rho;
      double Un = V[1]*n[0] + V[2]*n[1] + V[3]*n[2];
      double Vt[3] = {V[1] - Un*n[0],V[2] - Un*n[1],V[3] - Un*n[2]};
      rhoun = rho*(ujet - nVel);
      u = ujet*n[0] + Vt[0];
      v = ujet*n[1] + Vt[1];
      w = ujet*n[2] + Vt[2];
      double u1 = ujet*n[0]/porosity + Vt[0];
      double v1 = ujet*n[1]/porosity + Vt[1];
      double w1 = ujet*n[2]/porosity + Vt[2];
      q = u1*u1 + v1*v1 + w1*w1;
      p = V[4];
    } 
    flux[0] = S * rhoun;
    flux[1] = S * (rhoun*u + p*n[0]);
    flux[2] = S * (rhoun*v + p*n[1]);
    flux[3] = S * (rhoun*w + p*n[2]);
    flux[4] = S * (rhoun*(gam*invgam1*(p+pstiff)/rho + 0.5*q) + p*nVel);
  }
  else {
    flux[0] = 0.0;
    flux[1] = V[4] * normal[0];
    flux[2] = V[4] * normal[1];
    flux[3] = V[4] * normal[2];
    flux[4] = V[4] * normalVel;
  }

  flux[5] = 0.0;
  flux[6] = 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGPorousWallKE3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{
// not implemented
}

//------------------------------------------------------------------------------

void FluxFcnSGGhidagliaKE3D::compute(double length, double irey, double *normal, double normalVel,
                                   double *V, double *Ub, double *flux, bool useLimiter)
{

  F77NAME(genbcfluxgas)(2, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------
                                                                                                                  
void FluxFcnSGOutflowKE3D::compute(double length, double irey, double *normal, double normalVel,
                                 double *V, double *Ub, double *flux, bool useLimiter)
{

  F77NAME(boundflux5)(2, vf->getGamma(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void FluxFcnSGOutflowKE3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  F77NAME(gxboundflux5)(2, vf->getGamma(), normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------
//
void FluxFcnSGDirectStateInflowKE3D::compute(double length, double irey, double *normal, double normalVel,
                                        double *V, double *Ub, double *flux, bool useLimiter)
{

  directstateflux3D<7>(2, vf, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGDirectStateInflowKE3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  directstateflux3DDerivative<7>(2, vf, normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------

void FluxFcnSGMassFlowInflowKE3D::compute(double length, double irey, double *normal, double normalVel,
                                        double *V, double *Ub, double *flux, bool useLimiter)
{

  massflowflux3D<7>(2, vf, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGMassFlowInflowKE3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  massflowflux3DDerivative<7>(2, vf, normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------

void FluxFcnSGDirectStateOutflowKE3D::compute(double length, double irey, double *normal, double normalVel,
                                         double *V, double *Ub, double *flux, bool useLimiter)
{

  directstateflux3D<7>(2, vf, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGDirectStateOutflowKE3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{

  directstateflux3DDerivative<7>(2, vf, normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------

void FluxFcnSGMassFlowOutflowKE3D::compute(double length, double irey, double *normal, double normalVel,
                                         double *V, double *Ub, double *flux, bool useLimiter)
{

  massflowflux3D<7>(2, vf, normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// Included (MB)
void FluxFcnSGMassFlowOutflowKE3D::computeDerivative
(
  double irey, double dIrey, double *normal, double *dNormal,
  double normalVel, double dNormalVel, double *V,
  double *Ub, double *dUb, double *flux, double *dFlux
)
{
 
  massflowflux3DDerivative<7>(2, vf, normal, dNormal, normalVel, dNormalVel, V, Ub, dUb, flux, dFlux);

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// note: jacL = dFdUL and jacR = dFdUR
                                                                                                                  
void FluxFcnSGRoeKEturb3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                          double *VL, double *VR,
                                          double *jacL, double *jacR, bool useLimiter)
{

  F77NAME(roejac2)(2, gamma, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, VL, VR, jacL, jacR, sprec.getMinMach(),sprec.getSlope(),sprec.getCutOffMach(),irey,useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnSGWallKEturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                          double *V, double *Ub, double *jacL, bool useLimiter)
{

  jacL[0] = 0.0;
  jacL[1] = 0.0;
  jacL[2] = 0.0;
  jacL[3] = 0.0;

}

//------------------------------------------------------------------------------

void FluxFcnSGPorousWallKEturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                          double *V, double *Ub, double *jacL, bool useLimiter)
{

  jacL[0] = 0.0;
  jacL[1] = 0.0;
  jacL[2] = 0.0;
  jacL[3] = 0.0;

}

//------------------------------------------------------------------------------

void FluxFcnSGGhidagliaKEturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                             double *V, double *Ub, double *jacL, bool useLimiter)
{

  double flux[7];
  F77NAME(genbcfluxgas)(2, vf->getGamma(), vf->getPressureConstant(), normal, normalVel, V, Ub, flux);

  double updir = 1.0;
  if(flux[0]<0.0)  updir = 0.0;
  else if(flux[0]==0.0) updir = 0.5;

  jacL[0] = flux[0] * updir / V[0];
  jacL[1] = 0.0;
  jacL[2] = 0.0;
  jacL[3] = jacL[0];

}

//------------------------------------------------------------------------------
// note: jacL = dFdUL

void FluxFcnSGOutflowKEturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                             double *V, double *Ub, double *jacL, bool useLimiter)
{

  F77NAME(boundjac2)(2, vf->getGamma(), normal, normalVel, V, Ub, jacL);

}

//------------------------------------------------------------------------------

void FluxFcnSGDirectStateInflowKEturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                    double *V, double *Ub, double *jacL, bool useLimiter)
{

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;

  double rho = Ub[0];
  double u = Ub[1]/Ub[0];
  double v = Ub[2]/Ub[0];
  double w = Ub[3]/Ub[0];
  double un = u*n[0] + v*n[1] + w*n[2] - nVel;

  if (typeJac == FluxFcnBase::CONSERVATIVE)
    jacL[0] = S * un;
  else
    jacL[0] = S * rho*un;
  jacL[1] = 0.0;
  jacL[2] = 0.0;
  jacL[3] = jacL[0];

}

//------------------------------------------------------------------------------

void FluxFcnSGMassFlowInflowKEturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                    double *V, double *Ub, double *jacL, bool useLimiter)
{

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;

  double rho = Ub[0];
  double u = Ub[1]/Ub[0];
  double v = Ub[2]/Ub[0];
  double w = Ub[3]/Ub[0];
  double un = u*n[0] + v*n[1] + w*n[2] - nVel;

  if (typeJac == FluxFcnBase::CONSERVATIVE)
    jacL[0] = S * un;
  else
    jacL[0] = S * rho*un;
  jacL[1] = 0.0;
  jacL[2] = 0.0;
  jacL[3] = jacL[0];

}

//------------------------------------------------------------------------------

void FluxFcnSGDirectStateOutflowKEturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                     double *V, double *Ub, double *jacL, bool useLimiter)
{

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;

  double rho = Ub[0];
  double u = Ub[1]/Ub[0];
  double v = Ub[2]/Ub[0];
  double w = Ub[3]/Ub[0];
  double un = u*n[0] + v*n[1] + w*n[2] - nVel;

  if (typeJac == FluxFcnBase::CONSERVATIVE)
    jacL[0] = S * un;
  else
    jacL[0] = S * rho*un;
  jacL[1] = 0.0;
  jacL[2] = 0.0;
  jacL[3] = jacL[0];

}

//------------------------------------------------------------------------------

void FluxFcnSGMassFlowOutflowKEturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                                     double *V, double *Ub, double *jacL, bool useLimiter)
{

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;

  double rho = Ub[0];
  double u = Ub[1]/Ub[0];
  double v = Ub[2]/Ub[0];
  double w = Ub[3]/Ub[0];
  double un = u*n[0] + v*n[1] + w*n[2] - nVel;

  if (typeJac == FluxFcnBase::CONSERVATIVE)
    jacL[0] = S * un;
  else
    jacL[0] = S * rho*un;
  jacL[1] = 0.0;
  jacL[2] = 0.0;
  jacL[3] = jacL[0];

}

//------------------------------------------------------------------------------
