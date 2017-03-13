#include <FluxFcnDescTait.h>

#include <LinkF77.h>

#include <cstdlib>
#include <cstdio>


//------------------------------------------------------------------------------
// fortran routines located in f77src folder

extern "C" {
  void F77NAME(roeflux1water)(const double&, const double&, const double&,
                              const double&, const double&, double*,
                              const double&, double*, double*, double*);
  void F77NAME(roeflux5waterdissprec)(const int&, const double&, const double&, const double&,
                         const double&, const double&, double*,
                         const double&, double*, double*, double*, double*, double*,
                         const double&, const double&, const double&, const double&, const int&);
  void F77NAME(roeflux5waterburn)(const int&, const double&, const double&, const double&,
                         const double&, const double&, double*,
                         const double&, double*, double*, double*, double*, double*,
                         const double&, const double&, const double&, const double&, const int&);
  void F77NAME(roejac5waterdissprec)(const int&, const double&, const double&, const double&,
                         const double&, const double&, double*,
                         const double&, double*, double*, double*, 
			 const double&, const double&, const double&, 
			 const double&, const int&);
  void F77NAME(genbcfluxtait)(const int&, const double&, const double&, const double&, 
			 const double&,  double*, const double&, double*, double*, double*);
  void F77NAME(genbcfluxtait_hh)(const int&, const double&, const double&, const double&,
                         const double&,  double*, const double&, double*, double*, double*,double*, double&, const double&, const double&);
};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

template<int dim>
inline
void flux3Dwater(int type, VarFcnBase *vf, double *normal, double normalVel, double *V, double *Ub, double *flux){
 
  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;
  
  double Vb[dim];
  vf->conservativeToPrimitive(Ub, Vb);

  double VV[5];
  double machb = vf->computeMachNumber(Vb);
  double cb = vf->computeSoundSpeed(Vb);
  double unb = Vb[1]*n[0] + Vb[2]*n[1] + Vb[3]*n[2];//-nVel;
  double un = V[1]*n[0] + V[2]*n[1] + V[3]*n[2];
  double c = vf->computeSoundSpeed(V);
  double rhoun, p, rhoE;
  
  if(un==0.0){                   // SLIP-WALL LIKE, if boundary velocity is in the same plane as the face of the boundary
    VV[0]=vf->getDensity(Vb);
    VV[1]=0.0; VV[2]=0.0; VV[3]=0.0; VV[4]=0.0;

    rhoun=0.0;
    p = vf->getPressure(VV);
    rhoE = 0.0;
  }else{
    if(un<0.0){ 		//INLET, as the normal is going outward.
      if(-un-c>0.0){ 			//SUPERSONIC
	VV[0]=vf->getDensity(Vb);
	VV[1]=Vb[1];
	VV[2]=Vb[2];
	VV[3]=Vb[3];
	VV[4]=vf->computeTemperature(Vb);
      }else{        			//SUBSONIC
	VV[0] = vf->getDensity(V);
	VV[1] = Vb[1];
	VV[2] = Vb[2];
	VV[3] = Vb[3];
	VV[4] = vf->computeTemperature(Vb);
      }
    }
    else{                  //OUTLET
      if(un-c>0.0){                       //SUPERSONIC
	VV[0] = vf->getDensity(V);
	VV[1] = V[1];
	VV[2] = V[2];
	VV[3] = V[3];
	VV[4] = vf->computeTemperature(V);
      }else{                              //SUBSONIC
	VV[0] = vf->getDensity(Vb);
	VV[1] = V[1];
	VV[2] = V[2];
	VV[3] = V[3];
        VV[4] = vf->computeTemperature(V);
      }  
    }
    rhoun = VV[0] * ( VV[1]*n[0] + VV[2]*n[1] + VV[3]*n[2] - nVel );
    p = vf->getPressure(VV);
    rhoE = vf->computeRhoEnergy(VV);
  }
  flux[0] = S * rhoun;
  flux[1] = S * (rhoun*VV[1] + p*n[0]);
  flux[2] = S * (rhoun*VV[2] + p*n[1]);
  flux[3] = S * (rhoun*VV[3] + p*n[2]);
  flux[4] = S * ((rhoE + p) * rhoun/VV[0] + p*nVel); 

}

//------------------------------------------------------------------------------



template<int dim>
inline
void jacflux3Dwater(int type, VarFcnBase *vf, FluxFcnBase::Type localTypeJac, double *normal,
		      double normalVel, double *V, double *Ub, double *jac){

  double S = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  double ooS = 1.0 / S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};
  double nVel = normalVel * ooS;

  double VV[dim];
  double Vb[dim];
  vf->conservativeToPrimitive(Ub, Vb);

  double unb = Vb[1]*n[0] + Vb[2]*n[1] + Vb[3]*n[2]; 
  double cb = vf->computeSoundSpeed(Vb);
  double un = V[1]*n[0]+V[2]*n[1]+V[3]*n[2];
  double c  = vf->computeSoundSpeed(V);

  double _dfdV[dim][dim];
  double *dfdV = reinterpret_cast<double *>(_dfdV);
  for(int k=0; k<dim*dim; ++k)
    dfdV[k] = 0.0;


  if(un == 0.0){                 //SLIP-WALL LIKE
    
    //nothing to do, the flux has only pressure terms
    //which depend on density only, which is the one 
    //at the boundary (not inside the volume)

  }else{  
    if(un < 0.0){ 		//INLET, as the normal is going outward.
      if(-un-c > 0.0){                  //SUPERSONIC
	
	//nothing to do: jacobian is null
	
      }else{                         //SUBSONIC
	
	//derivative wrt pressure
	VV[0] = vf->getDensity(V);
	VV[1] = Vb[1];
	VV[2] = Vb[2];
	VV[3] = Vb[3];
	VV[4] = vf->computeTemperature(Vb);
	double cp = vf->computeSoundSpeed(VV);
	double cp2 = cp*cp;
	double unp = VV[1]*n[0] +VV[2]*n[1] + VV[3]*n[2] - nVel;
	_dfdV[0][0] = S * unp;
	_dfdV[1][0] = S * (VV[1]*unp + cp2*n[0]);
	_dfdV[2][0] = S * (VV[2]*unp + cp2*n[1]);
	_dfdV[3][0] = S * (VV[3]*unp + cp2*n[2]);
	_dfdV[4][0] = S * unp * (cp2 + vf->computeRhoEnergy(VV)/vf->getDensity(VV));
      }
    }else{                  //OUTLET
      if(un-c > 0.0){                       //SUPERSONIC
	VV[0] = vf->getDensity(V);
	VV[1] = V[1];
	VV[2] = V[2];
	VV[3] = V[3];
	VV[4] = vf->computeTemperature(V);
	
	double unp = VV[1]*n[0]+VV[2]*n[1]+VV[3]*n[2]-nVel;
	double cp = vf->computeSoundSpeed(VV);
	double cp2 = cp*cp;
	
	_dfdV[0][0] = S * unp;
	_dfdV[0][1] = S * VV[0]*n[0];
	_dfdV[0][2] = S * VV[0]*n[1];
	_dfdV[0][3] = S * VV[0]*n[2];
	
	_dfdV[1][0] = S * (VV[1]*unp + cp2*n[0]);
	_dfdV[1][1] = S * VV[0]*(unp + VV[1]*n[0]);
	_dfdV[1][2] = S * VV[0]*VV[1]*n[1];
	_dfdV[1][3] = S * VV[0]*VV[1]*n[2];
	
	_dfdV[2][0] = S * (VV[2]*unp +cp2*n[1]);
	_dfdV[2][1] = S * VV[0]*VV[2]*n[0];
	_dfdV[2][2] = S * VV[0]*(unp + VV[2]*n[1]);
	_dfdV[2][3] = S * VV[0]*VV[2]*n[2];
	
	_dfdV[3][0] = S * (VV[3]*unp + cp2*n[2]);
	_dfdV[3][1] = S * VV[0]*VV[3]*n[0];
	_dfdV[3][2] = S * VV[0]*VV[3]*n[1];
	_dfdV[3][3] = S * VV[0]*(unp+VV[3]*n[2]);
	
	_dfdV[4][0] = S * unp * (cp2 + vf->computeRhoEnergy(VV)/vf->getDensity(VV));
	_dfdV[4][1] = S * (VV[0]*VV[1]*unp + n[0]*(vf->computeRhoEnergy(VV) + vf->getPressure(VV)));
	_dfdV[4][2] = S * (VV[0]*VV[2]*unp + n[1]*(vf->computeRhoEnergy(VV) + vf->getPressure(VV)));
	_dfdV[4][3] = S * (VV[0]*VV[3]*unp + n[2]*(vf->computeRhoEnergy(VV) + vf->getPressure(VV)));
	_dfdV[4][4] = S * vf->computeRhoEpsilon(VV)/vf->computeTemperature(VV) * unp;
	
      }else{//subsonic outlet
	
	VV[0] = vf->getDensity(Vb);
	VV[1] = V[1];
	VV[2] = V[2];
	VV[3] = V[3];
	VV[4] = vf->computeTemperature(V);
	
	double unp = VV[1]*n[0]+VV[2]*n[1]+VV[3]*n[2]-nVel;
	
	
	_dfdV[0][1] = S * VV[0]*n[0];
	_dfdV[0][2] = S * VV[0]*n[1];
	_dfdV[0][3] = S * VV[0]*n[2];
	
	
	_dfdV[1][1] = S * VV[0]*(unp + VV[1]*n[0]);
	_dfdV[1][2] = S * VV[0]*VV[1]*n[1];
	_dfdV[1][3] = S * VV[0]*VV[1]*n[2];
	
	
	_dfdV[2][1] = S * VV[0]*VV[2]*n[0];
	_dfdV[2][2] = S * VV[0]*(unp + VV[2]*n[1]);
	_dfdV[2][3] = S * VV[0]*VV[2]*n[2];
	
	
	_dfdV[3][1] = S * VV[0]*VV[3]*n[0];
	_dfdV[3][2] = S * VV[0]*VV[3]*n[1];
	_dfdV[3][3] = S * VV[0]*(unp+VV[3]*n[2]);
	
	
	_dfdV[4][1] = S * (VV[0]*VV[1]*unp + n[0]*(vf->computeRhoEnergy(VV) + vf->getPressure(VV)));
	_dfdV[4][2] = S * (VV[0]*VV[2]*unp + n[1]*(vf->computeRhoEnergy(VV) + vf->getPressure(VV)));
	_dfdV[4][3] = S * (VV[0]*VV[3]*unp + n[2]*(vf->computeRhoEnergy(VV) + vf->getPressure(VV)));
	_dfdV[4][4] = S * vf->computeRhoEpsilon(VV)/vf->computeTemperature(VV) * unp;
      }
    }
  }
  
  
  
  if (localTypeJac == FluxFcnBase::CONSERVATIVE)
    vf->postMultiplyBydVdU(V, dfdV, jac);
  else
    for (int k=0; k<dim*dim; ++k)
      jac[k] = dfdV[k]; 
  
}

//------------------------------------------------------------------------------

template<int dim>
inline
void roejactait3D(int type, double gamma, VarFcnBase *vf, FluxFcnBase::Type localTypeJac, double *normal,
          double normalVel, double *VL, double *VR, SpatialLowMachPrec sprec, double irey,
          double *jacL, double *jacR, bool useLimiter){


  const int dimm1 = dim-1;
  const int dimm2 = dim-2;
  const int dim2 = dim*dim;

  double dfdUL[dim2], dfdUR[dim2];
  for (int kk = 0; kk<dim2; kk++) dfdUR[kk] = 0.0;
  double n[3] = {normal[0], normal[1], normal[2]};
  F77NAME(roejac5waterdissprec)(type, gamma,
                                vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(),
                                n, normalVel,
                                VL, VR, dfdUL, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), irey, useLimiter ? sprec.getPrecTag() : 0);
  n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2];
  F77NAME(roejac5waterdissprec)(type, gamma,
                                vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(),
                                n, -normalVel,
                                VR, VL, dfdUR, sprec.getMinMach(), sprec.getSlope(), sprec.getCutOffMach(), irey, useLimiter ? sprec.getPrecTag() : 0);
  for (int k=0; k<dim2; k++)
    dfdUR[k] = -dfdUR[k];

  // if type is non-zero, then a turbulence model is being considered.
  if (type == 1 || type == 2) {
    double dfdVL[dim2], dfdVR[dim2];
    vf->postMultiplyBydUdV(VL, dfdUL, dfdVL);
    vf->postMultiplyBydUdV(VR, dfdUR, dfdVR);

    double flux[dim];
    F77NAME(roeflux5waterdissprec)(type, gamma, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(),
            normal, normalVel, VL, VL, VR, VR, flux, sprec.getMinMach(),
            sprec.getSlope(), sprec.getCutOffMach(), irey, useLimiter ? sprec.getPrecTag() : 0);
    double f1 = flux[0];

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

void FluxFcnTaitFDJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				     double *VL, double *VR, double *flux, bool useLimiter)
{

  fprintf(stderr, "*** Error: FluxFcnTaitFDJacRoeEuler3D::compute not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnTaitApprJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				       double *VL, double *VR, double *flux, bool useLimiter)
{
  if (!vf->isBurnable()) {
    F77NAME(roeflux5waterdissprec)(0, gamma, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(),
            normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux, sprec.getMinMach(),
            sprec.getSlope(), sprec.getCutOffMach(), irey, useLimiter ? sprec.getPrecTag() : 0);
  } else {
    F77NAME(roeflux5waterburn)(0, gamma, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(),
            normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux, sprec.getMinMach(),
            sprec.getSlope(), sprec.getCutOffMach(), irey, useLimiter ? sprec.getPrecTag() : 0);

  } 
}

//------------------------------------------------------------------------------

void FluxFcnTaitApprJacRoeEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel, 
						double *VL, double *VR, 
						double *jacL, double *jacR, bool useLimiter)
{
 
  roejactait3D<5>(0, gamma, vf, typeJac, normal, normalVel, VL, VR, sprec, irey, jacL, jacR, useLimiter);

}

//------------------------------------------------------------------------------

void FluxFcnTaitExactJacRoeEuler3D::compute(double length, double irey, double *normal, double normalVel, 
					double *VL, double *VR, double *flux, bool useLimiter)
{

  fprintf(stderr, "*** Error: FluxFcnTaitExactJacRoeEuler3D::compute not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnTaitExactJacRoeEuler3D::computeJacobians(double length, double irey, double *normal, double normalVel, 
						 double *VL, double *VR, 
						 double *jacL, double *jacR, bool useLimiter)
{

  fprintf(stderr, "*** Error: FluxFcnTaitExactJacRoeEuler3D::computeJacobians not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void FluxFcnTaitWallEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *V, double *Ub, double *flux, bool useLimiter)
{

  double P = vf->getPrefWater() + vf->getAlphaWater()*pow(V[0], vf->getBetaWater());
  flux[0] = 0.0;
  flux[1] = P * normal[0];
  flux[2] = P * normal[1];
  flux[3] = P * normal[2];
  flux[4] = P * normalVel;

}

//------------------------------------------------------------------------------

void FluxFcnTaitGhidagliaEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *V, double *Ub, double *flux, bool useLimiter)
{

  F77NAME(genbcfluxtait)(0, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, V, Ub, flux);

}
/*
void FluxFcnTaitGhidagliaEuler3D::computeJacobian(double length, double irey, double *normal, double normalVel, 
						  double *VL, double *Ub, double *jacL, bool useLimiter) {

}
*/
void FluxFcnTaitModifiedGhidagliaEuler3D::compute(double length, double irey, double *normal, double normalVel, 
                                                  double *V, double *Ub, double *flux, bool useLimiter)
{

  int dim = 5;
 // fprintf(stderr,"faceCenter = %e %e %e;  s_ff = %e;  dt = %e.\n", flux[dim], flux[dim+1], flux[dim+2], flux[dim+3], flux[dim+4]);
  double snew;
  F77NAME(genbcfluxtait_hh)(0, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, V, Ub, flux, hhcoeffptr, snew, *(hhcoeffptr+4), *(hhcoeffptr+5));
  //hhcoeffptr[3] = snew;
 
}

void FluxFcnTaitModifiedGhidagliaEuler3D::computeJacobianFarfield(double length, double irey, double *normal, double normalVel, 
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

void FluxFcnTaitInflowEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *V, double *Ub, double *flux, bool useLimiter)
{
  fprintf(stderr, "*** Error: FluxFcnTaitInflowEuler3D::compute not implemented\n");
  exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnTaitOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel, 
				    double *V, double *Ub, double *flux, bool useLimiter)
{
  fprintf(stderr, "*** Error: FluxFcnTaitOutflowEuler3D::compute not implemented\n");
  exit(1);
}

//------------------------------------------------------------------------------

void FluxFcnTaitInternalInflowEuler3D::compute(double length, double irey, double *normal, double normalVel, 
					   double *V, double *Ub, double *flux, bool useLimiter)
{

  flux3Dwater<5>(0,vf,normal,normalVel,V,Ub,flux);

}

//------------------------------------------------------------------------------

void FluxFcnTaitInternalInflowEuler3D::computeJacobian(double length, double irey, double *normal, double normalVel, 
						   double *V, double *Ub, double *jacL, bool useLimiter)
{

  jacflux3Dwater<5>(0,vf,typeJac,normal,normalVel,V,Ub,jacL);

}

//------------------------------------------------------------------------------

void FluxFcnTaitInternalOutflowEuler3D::compute(double length, double irey, double *normal, double normalVel, 
					    double *V, double *Ub, double *flux, bool useLimiter)
{

  flux3Dwater<5>(0,vf,normal,normalVel,V,Ub,flux);

}

//------------------------------------------------------------------------------

void FluxFcnTaitInternalOutflowEuler3D::computeJacobian(double length, double irey, double *normal, double normalVel, 
						    double *V, double *Ub, double *jacL, bool useLimiter)
{
 
  jacflux3Dwater<5>(0,vf,typeJac,normal,normalVel,V,Ub,jacL);

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void FluxFcnTaitApprJacRoeSA3D::compute(double length, double irey, double *normal, double normalVel,
               double *VL, double *VR, double *flux, bool useLimiter)
{

  F77NAME(roeflux5waterdissprec)(1, gamma, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(),
          normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux, sprec.getMinMach(),
          sprec.getSlope(), sprec.getCutOffMach(), irey, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnTaitApprJacRoeSA3D::computeJacobians(double length, double irey, double *normal, double normalVel,
            double *VL, double *VR,
            double *jacL, double *jacR, bool useLimiter)
{

  roejactait3D<6>(1, gamma, vf, typeJac, normal, normalVel, VL, VR, sprec, irey, jacL, jacR, useLimiter);

}

//------------------------------------------------------------------------------

void FluxFcnTaitWallSA3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *V, double *Ub, double *flux, bool useLimiter)
{

  double P = vf->getPrefWater() + vf->getAlphaWater()*pow(V[0], vf->getBetaWater());
  flux[0] = 0.0;
  flux[1] = P * normal[0];
  flux[2] = P * normal[1];
  flux[3] = P * normal[2];
  flux[4] = P * normalVel;
  flux[5] = 0.0;

}

//------------------------------------------------------------------------------

void FluxFcnTaitGhidagliaSA3D::compute(double length, double irey, double *normal, double normalVel,
           double *V, double *Ub, double *flux, bool useLimiter)
{

  F77NAME(genbcfluxtait)(1, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------
// note: jacL = dFdUL and jacR = dFdUR

void FluxFcnTaitRoeSAturb3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                          double *VL, double *VR,
                                          double *jacL, double *jacR, bool useLimiter)
{

  double flux[6];
  F77NAME(roeflux5waterdissprec)(1, gamma, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(),
          normal, normalVel, VL, VL, VR, VR, flux, sprec.getMinMach(),
          sprec.getSlope(), sprec.getCutOffMach(), irey, useLimiter ? sprec.getPrecTag() : 0);

  double updir = 1.0;
  if(flux[0]<0.0)  updir = 0.0;
  else if(flux[0]==0.0) updir = 0.5;

  jacL[0] = flux[0] * updir / VL[0];
  jacR[0] = flux[0] * (1.0 - updir) / VR[0];

}

//------------------------------------------------------------------------------

void FluxFcnTaitWallSAturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                          double *V, double *Ub, double *jacL, bool useLimiter)
{

  jacL[0] = 0.0;

}

//------------------------------------------------------------------------------

void FluxFcnTaitGhidagliaSAturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                             double *V, double *Ub, double *jacL, bool useLimiter)
{

  double flux[6];
  F77NAME(genbcfluxtait)(1, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, V, Ub, flux);

  double updir = 1.0;
  if(flux[0]<0.0)  updir = 0.0;
  else if(flux[0]==0.0) updir = 0.5;

  jacL[0] = flux[0] * updir / V[0];

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

void FluxFcnTaitApprJacRoeKE3D::compute(double length, double irey, double *normal, double normalVel,
               double *VL, double *VR, double *flux, bool useLimiter)
{

  F77NAME(roeflux5waterdissprec)(2, gamma, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(),
          normal, normalVel, VL, VL+rshift, VR, VR+rshift, flux, sprec.getMinMach(),
          sprec.getSlope(), sprec.getCutOffMach(), irey, useLimiter ? sprec.getPrecTag() : 0);

}

//------------------------------------------------------------------------------

void FluxFcnTaitApprJacRoeKE3D::computeJacobians(double length, double irey, double *normal, double normalVel,
            double *VL, double *VR,
            double *jacL, double *jacR, bool useLimiter)
{

  roejactait3D<7>(2, gamma, vf, typeJac, normal, normalVel, VL, VR, sprec, irey, jacL, jacR, useLimiter);

}

//------------------------------------------------------------------------------

void FluxFcnTaitWallKE3D::compute(double length, double irey, double *normal, double normalVel, 
				   double *V, double *Ub, double *flux, bool useLimiter)
{

  double P = vf->getPrefWater() + vf->getAlphaWater()*pow(V[0], vf->getBetaWater());
  flux[0] = 0.0;
  flux[1] = P * normal[0];
  flux[2] = P * normal[1];
  flux[3] = P * normal[2];
  flux[4] = P * normalVel;
  flux[5] = 0.0;
  flux[6] = 0.0;

}

//------------------------------------------------------------------------------

void FluxFcnTaitGhidagliaKE3D::compute(double length, double irey, double *normal, double normalVel,
           double *V, double *Ub, double *flux, bool useLimiter)
{

  F77NAME(genbcfluxtait)(2, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, V, Ub, flux);

}

//------------------------------------------------------------------------------

// note: jacL = dFdUL and jacR = dFdUR

void FluxFcnTaitRoeKEturb3D::computeJacobians(double length, double irey, double *normal, double normalVel,
                                          double *VL, double *VR,
                                          double *jacL, double *jacR, bool useLimiter)
{

  double flux[7];
  F77NAME(roeflux5waterdissprec)(2, gamma, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(),
          normal, normalVel, VL, VL, VR, VR, flux, sprec.getMinMach(),
          sprec.getSlope(), sprec.getCutOffMach(), irey, useLimiter ? sprec.getPrecTag() : 0);

  double updir = 1.0;
  if(flux[0]<0.0)  updir = 0.0;
  else if(flux[0]==0.0) updir = 0.5;

  double jacleft = flux[0] * updir / VL[0];
  double jacright = flux[0] * (1.0 - updir) / VR[0];

  jacL[0] = jacleft;
  jacL[1] = 0.0;
  jacL[2] = 0.0;
  jacL[3] = jacL[0];

  jacR[0] = jacright;
  jacR[1] = 0.0;
  jacR[2] = 0.0;
  jacR[3] = jacR[0];

}

//------------------------------------------------------------------------------

void FluxFcnTaitWallKEturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                          double *V, double *Ub, double *jacL, bool useLimiter)
{

  jacL[0] = 0.0;
  jacL[1] = 0.0;
  jacL[2] = 0.0;
  jacL[3] = 0.0;

}

//------------------------------------------------------------------------------

void FluxFcnTaitGhidagliaKEturb3D::computeJacobian(double length, double irey, double *normal, double normalVel,
                                             double *V, double *Ub, double *jacL, bool useLimiter)
{

  double flux[7];
  F77NAME(genbcfluxtait)(2, vf->getCv(), vf->getPrefWater(), vf->getAlphaWater(), vf->getBetaWater(), normal, normalVel, V, Ub, flux);

  double updir = 1.0;
  if(flux[0]<0.0)  updir = 0.0;
  else if(flux[0]==0.0) updir = 0.5;

  jacL[0] = flux[0] * updir / V[0];
  jacL[1] = 0.0;
  jacL[2] = 0.0;
  jacL[3] = jacL[0];

}


