#include <FluxFcnDescJwl.h>

#include <LinkF77.h>

#include <cstdlib>
#include <cstdio>

#ifdef USE_MPI
#include <mpi.h>
#endif

//------------------------------------------------------------------------------
// fortran routines located in f77src folder

extern "C" {

  void F77NAME(roeflux5jwl)(const int&, const double&, const double&, const double&, 
                            const double&, const double&, const double&,
                            double*, const double&,double*,double*,double*,double*,double*,
                            const double&, const double&, const double&, const double&, 
                            const double&, const double&, const int&);
  void F77NAME(roejac6jwl)(const int&, const double&, const double&, const double&, 
                           const double&, const double&, const double&, double*, 
			   const double&, double*, double*, double*,const int&,
                           const double&, const double&, const double&, const double&, 
                           const double&, const double&, const int&);
  void F77NAME(genbcfluxjwl)(const int&, const double&, const double&, 
                             const double&, const double&, const double&,
                             double*, const double&,double*,double*,double*);

};


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void FluxFcnJwlFDJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				     double *VL, double *VR, double *flux, bool useLimiter)
{

   fprintf(stderr, "*** Error: FluxFcnJwlFDJacRoeEuler3D::compute not implemented\n");
   exit(1);

}

//------------------------------------------------------------------------------

FILE* ftest = NULL;

void FluxFcnJwlApprJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				       double *VL, double *VR, double *flux, bool useLimiter)
{
  /*if (!ftest) {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    char fname[80];
    sprintf(fname,"test%d.log",myrank);
    ftest = fopen(fname,"w");
  }

  std::stringstream str;
  str << VL[0] << " " << VL[1] << " " << VL[2] << " " << VL[3] << " " << VL[4] << std::endl;
  str << VR[0] << " " << VR[1] << " " << VR[2] << " " << VR[3] << " " << VR[4] << std::endl << std::endl;
  fprintf(ftest,"%s",str.str().c_str());
  fflush(ftest);
  */
   F77NAME(roeflux5jwl)(0, gamma, vf->getOmega(), vf->getA1(), vf->getA2(), vf->getR1r(), vf->getR2r(),
                        normal, normalVel, 
                        VL, VL+rshift, VR, VR+rshift, flux, 
                        sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), 
                        sprec.getShockParameter(), irey, length, useLimiter ? sprec.getPrecTag() : 0);
 
}

//------------------------------------------------------------------------------

void FluxFcnJwlApprJacRoeEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel, 
						double *VL, double *VR, 
						double *jacL, double *jacR, bool useLimiter)
{

  const int dim = 5;
  const int dim2 = dim*dim;
  const int type = 0; //no turbulence


  double dfdUL[dim2], dfdUR[dim2];
  for (int kk = 0; kk<dim2; kk++) dfdUR[kk] = 0.0;
  double n[3] = {normal[0], normal[1], normal[2]};
  F77NAME(roejac6jwl)(type, gamma, vf->getOmega(), vf->getA1(), vf->getA2(), vf->getR1r(), vf->getR2r(), n, normalVel, 
                      VL, VR, dfdUL, 1, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), 
                      irey, length, useLimiter ? sprec.getPrecTag() : 0);
  F77NAME(roejac6jwl)(type, gamma, vf->getOmega(), vf->getA1(), vf->getA2(), vf->getR1r(), vf->getR2r(), n, normalVel,
                      VR, VL, dfdUR, 2, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), sprec.getShockParameter(), 
                      irey, length, useLimiter ? sprec.getPrecTag() : 0);
  if (typeJac == FluxFcnBase::CONSERVATIVE) {
    for (int k=0; k<dim2; ++k) {
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

void FluxFcnJwlExactJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel, 
					double *VL, double *VR, double *flux, bool useLimiter)
{

  fprintf(stderr, "*** Error: FluxFcnJwlExactJacRoeEuler3D::compute not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnJwlExactJacRoeEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel, 
						 double *VL, double *VR, 
						 double *jacL, double *jacR, bool useLimiter)
{

  fprintf(stderr, "*** Error: FluxFcnJwlExactJacRoeEuler3D::computeJacobians not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnJwlWallEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *V, double *Ub, double *flux, bool useLimiter)
{

  flux[0] = 0.0;
  flux[1] = V[4] * normal[0];
  flux[2] = V[4] * normal[1];
  flux[3] = V[4] * normal[2];
  flux[4] = V[4] * normalVel;

}

//------------------------------------------------------------------------------

void FluxFcnJwlGhidagliaEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *V, double *Ub, double *flux, bool useLimiter)
{

  F77NAME(genbcfluxjwl)(0, vf->getOmega(), vf->getA1(), vf->getA2(),
                        vf->getR1r(), vf->getR2r(),
                        normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

void FluxFcnJwlInflowEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *V, double *Ub, double *flux, bool useLimiter)
{
    fprintf(stderr, "*** Error: FluxFcnJwlInflowEuler3D::compute not implemented\n");
    exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnJwlOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				    double *V, double *Ub, double *flux, bool useLimiter)
{
    fprintf(stderr, "*** Error: FluxFcnJwlOutflowEuler3D::compute not implemented\n");
    exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnJwlInternalInflowEuler3D::compute(double length, double irey, double *normal, double normalVel, 
					   double *V, double *Ub, double *flux, bool useLimiter)
{

    fprintf(stderr, "*** Error: FluxFcnJwlInternalInflowEuler3D::compute not implemented\n");
    exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnJwlInternalInflowEuler3D::computeJacobian(double length, double irey, double *normal, double normalVel, 
						   double *V, double *Ub, double *jacL, bool useLimiter)
{

    fprintf(stderr, "*** Error: FluxFcnJwlInternalInflowEuler3D::computeJacobian not implemented\n");
    exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnJwlInternalOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel, 
					    double *V, double *Ub, double *flux, bool useLimiter)
{

    fprintf(stderr, "*** Error: FluxFcnJwlInternalOutflowEuler3D::compute not implemented\n");
    exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnJwlInternalOutflowEuler3D::computeJacobian(double length, double irey, double *normal, double normalVel, 
						    double *V, double *Ub, double *jacL, bool useLimiter)
{
 
    fprintf(stderr, "*** Error: FluxFcnJwlInternalOutflowEuler3D::computeJacobian not implemented\n");
    exit(1);

}

//------------------------------------------------------------------------------

