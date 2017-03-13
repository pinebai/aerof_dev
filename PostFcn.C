#include <PostFcn.h>

#include <WallFcn.h>
#include <Vector3D.h>

#include <cstdlib>
#include <cstdio>

#include "./Dev/devtools.h"//TODO delete line

const double PostFcnEuler::third = 1.0/3.0;
//------------------------------------------------------------------------------

PostFcn::PostFcn(VarFcn *vf)
{

  varFcn = vf;

}

//------------------------------------------------------------------------------

double PostFcn::computeNodeScalarQuantity(ScalarType type, double *V, double *X, int fluidId,double* phi)
{

  fprintf(stderr, "*** Warning: computeNodeScalarQuantity not defined\n");

  return 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcn::computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType type, double dS[3], double *X, double *dX, double *V, double *dV, double phi)
{

  fprintf(stderr, "*** Warning: computeDerivativeOfNodeScalarQuantity not defined\n");

  return 0.0;

}

//------------------------------------------------------------------------------

double PostFcn::computeFaceScalarQuantity(ScalarType type, double dp1dxj[4][3], 
					  Vec3D& n, double d2w[3], double* Vwall, 
					  double* Vface[3], double* Vtet[4])
{

  fprintf(stderr, "*** Warning: computeFaceScalarQuantity not defined\n");

  return 0.0;

}

//------------------------------------------------------------------------------

PostFcnEuler::PostFcnEuler(IoData &iod, VarFcn *vf) : PostFcn(vf)
{

  if (iod.eqs.fluidModel.fluid == FluidModelData::PERFECT_GAS ||
      iod.eqs.fluidModel.fluid == FluidModelData::STIFFENED_GAS){
    mach = iod.ref.mach;
    pinfty = iod.bc.inlet.pressure;

// Included (MB)
    dpinfty = -2.0 / (iod.eqs.fluidModel.gasModel.specificHeatRatio * mach*mach*mach);

    if (iod.problem.mode == ProblemData::DIMENSIONAL) {
      dimFlag = true;
    }
    else if (iod.problem.mode == ProblemData::NON_DIMENSIONAL) {
      dimFlag = false;
      if (iod.sa.apressFlag == false) {
        if (iod.sa.pressFlag == false) {
          dPin = dpinfty;
        }
        else {
          dPin = 0.0;
        }      
      }
      else {
        dPin = 0.0;
      }
    }

  } else if (iod.eqs.fluidModel.fluid == FluidModelData::LIQUID){
    mach = iod.ref.mach;
    double P = vf->getPrefWater();
    double a = vf->getAlphaWater();
    double b = vf->getBetaWater();
    pinfty = (P+a*pow(iod.bc.inlet.density, b));
  } else if (iod.eqs.fluidModel.fluid == FluidModelData::JWL){
    mach = iod.ref.mach;
    pinfty = iod.bc.inlet.pressure;
  }

}

//------------------------------------------------------------------------------

// Included (MB)
void PostFcnEuler::rstVar(IoData &iod, Communicator *com)
{

  mach = iod.ref.mach;
  pinfty = iod.bc.inlet.pressure;
  dpinfty = -2.0 / (iod.eqs.fluidModel.gasModel.specificHeatRatio * mach*mach*mach);

  if (iod.problem.mode == ProblemData::DIMENSIONAL) {
    dimFlag = true;
  }
  else if (iod.problem.mode == ProblemData::NON_DIMENSIONAL) {
    dimFlag = false;
    if (iod.sa.apressFlag == false) {
      if (iod.sa.pressFlag == false) {
        dPin = dpinfty;
      }
      else {
        dPin = 0.0;
      }      
    }
    else {
      dPin = 0.0;
    }
  }
  
}

//------------------------------------------------------------------------------

double PostFcnEuler::computeNodeScalarQuantity(ScalarType type, double *V, double *X, int fluidId,double* phi)
{
  double q = 0.0;
  double n[3];

  if (type == DENSITY)
    q = varFcn->getDensity(V, fluidId);
  else if (type == MACH)
    q = varFcn->computeMachNumber(V, fluidId);
  else if (type == WTMACH)
    q = varFcn->computeWtMachNumber(V, fluidId);
  else if (type == SPEED)
    q = sqrt(varFcn->computeU2(V));
  else if (type == WTSPEED)
    q = sqrt(varFcn->computeWtU2(V));
  else if (type == PRESSURE) {
    q = varFcn->getPressure(V, fluidId);
  }
  else if (type == DIFFPRESSURE)
    q = varFcn->getPressure(V, fluidId)-pinfty;
  else if (type == TEMPERATURE)
    q = varFcn->computeTemperature(V, fluidId);
  else if (type == TOTPRESSURE)
    q = varFcn->computeTotalPressure(mach, V, fluidId);
  else if (type == NUT_TURB)
    q = varFcn->getTurbulentNuTilde(V, fluidId);
  else if (type == K_TURB)
    q = varFcn->getTurbulentKineticEnergy(V, fluidId);
  else if (type == EPS_TURB)
    q = varFcn->getTurbulentDissipationRate(V,fluidId);
  else if(type == HYDROSTATICPRESSURE)
    q = varFcn->hydrostaticPressure(V[0],X);
  else if(type == HYDRODYNAMICPRESSURE)
    q = varFcn->hydrodynamicPressure(V,X,fluidId);
  else if(type == PRESSURECOEFFICIENT)
    q = varFcn->computePressureCoefficient(V, pinfty, mach, dimFlag,fluidId);
  else if(type == PHILEVEL)
    //q = static_cast<double>(fluidId);
    q = phi[0];///varFcn->getDensity(V, fluidId);
  else if(type == PHILEVEL2)
    //q = static_cast<double>(fluidId);
    q = phi[1];///varFcn->getDensity(V, fluidId);
  else if (type == FLUIDID)
    q = static_cast<double>(fluidId);
 // Included (MB)
  else if (type == VELOCITY_NORM)
    q = varFcn->getVelocityNorm(V,fluidId);

  return q;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcnEuler::computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType type,
							   double dS[3],
							   double *V,
							   double *dV,
							   double *X,
							   double *dX,
							   double phi)
{
  //Dev::Error(MPI_COMM_WORLD,"tempoarary Error to look at backtrace",true);//TODO delete line

  double q = 0.0;

  if (type == DERIVATIVE_DENSITY)
    q = varFcn->getDensity(dV);
  else if (type == DERIVATIVE_MACH)
    q = varFcn->computeDerivativeOfMachNumber(V, dV, dS[0]);
  else if (type == DERIVATIVE_PRESSURE)
    q = varFcn->getPressure(dV);
  else if (type == DERIVATIVE_TEMPERATURE)
    q = varFcn->computeDerivativeOfTemperature(V, dV);
  else if (type == DERIVATIVE_TOTPRESSURE)
    q = varFcn->computeDerivativeOfTotalPressure(mach, dS[0], V, dV, dS[0]);
  else if (type == DERIVATIVE_NUT_TURB)
    q = varFcn->getTurbulentNuTilde(dV);
  else if (type == DERIVATIVE_VELOCITY_SCALAR)
    q = varFcn->getDerivativeOfVelocityNorm(V, dV);
  else if (type == DERIVATIVE_EDDY_VISCOSITY)//not existing in Euler
    {fprintf(stderr, "*** ERROR derivative of eddy viscosity cannot be computed in an Euler simulation\n"); exit(-1);}
  else if (type == DERIVATIVE_SPATIAL_RES)//not existing in Euler
    {fprintf(stderr, "*** ERROR derivative of eddy viscosity cannot be computed in an Euler simulation\n"); exit(-1);}
  else if (type == DSSIZE)//not existing in Euler
    {fprintf(stderr, "*** ERROR derivative of eddy viscosity cannot be computed in an Euler simulation\n"); exit(-1);}
  else
  {
    // Error message
    fprintf(stderr, "*** ERROR: PostFcnEuler::computeDerivativeOfNodeScalarQuantity does not define the type %d\n", type);
    exit(-1);
  }

  return q;

}

//------------------------------------------------------------------------------

void PostFcnEuler::computeForce(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double d2w[3],
                                double *Vwall, double *Vface[3], double *Vtet[4],
				double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, double dPdx[3][3], int hydro,
				int fid)
{

  double pcg[3], p[3];
  double pcgin;
  int i;

  switch(hydro)
    {
    case 0:
      for(i=0;i<3;i++) pcg[i] = varFcn->getPressure(Vface[i],fid);
      break;
    case 1: // hydrostatic pressure
      for(i=0;i<3;i++) pcg[i] = varFcn->hydrostaticPressure(Vface[i][0],Xface[i]);
      break;
    case 2: // hydrodynamic pressure
      for(i=0;i<3;i++) pcg[i] = varFcn->hydrodynamicPressure(Vface[i],Xface[i]);
      break;
    default:
      fprintf(stderr,"hydro parameter is not correct. Pressure at the face cannot be computed. hydro = %d\n",hydro);
      exit(-1);
    }

  if (pin)
    pcgin = *pin;
  else
    if (hydro == 0)
      pcgin = pinfty;
    else
      pcgin = 0.0;

  p[0] = (pcg[0] - pcgin) ;
  p[1] = (pcg[1] - pcgin) ;
  p[2] = (pcg[2] - pcgin) ;

// ##################
 Vec3D x0, x1, x2, x;
 double temp;
 for(i = 0; i<3; i++)
 {
  x0[i] = Xface[0][i];
  x1[i] = Xface[1][i];
  x2[i] = Xface[2][i];
 }

// for node 0
 i=0; // i represents node
 x = (7.0/18.0)*(x1 + x2 - 2.0*x0);
 temp = 2.0*p[i] + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 Fi0 = (1.0/6.0*temp)*n;

// for node 1
 i=1;
 x = (7.0/18.0)*(x2 + x0 - 2.0*x1);
 temp = 2.0*p[i] + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 Fi1 = (1.0/6.0*temp)*n;

// for node 2
 i=2;
 x = (7.0/18.0)*(x0 + x1 - 2.0*x2);
 temp = 2.0*p[i] + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 Fi2 = (1.0/6.0*temp)*n;

 Fv = 0.0;

}
//------------------------------------------------------------------------------

void PostFcnEuler::computeForceEmbedded(int orderOfAccuracy, double dp1dxj[4][3], 
					double *Xface[3], Vec3D &n, double d2w[3],
					double *Vwall, double *Vface[3], double *Vtet[4], double pin, 
					Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, double dPdx[3][3], 
                                        int hydro, int* fid, bool applyRealForce)
{
  if(hydro!=0) {fprintf(stderr,"hydro parameter (%d) not supported...\n",hydro); exit(-1);}
  Vec3D p;
  if(applyRealForce)
    for(int i=0;i<3;i++) 
      p[i] = varFcn->getPressure(Vface[i],(fid?fid[i]:0)) - pin;
  else
    for(int i=0;i<3;i++)
      p[i] = Vface[i][4] - pin;

  // Computes Int_{T} (N_i * Sum_{j\in T} p_j N_j) dx
  // At first order, N_j = Chi_j (constant per control volume)
  // At second order, N_j = Phi_j (P1 lagrangian basis functions within T)
  switch (orderOfAccuracy) {
    case 1: 
      Fi0 = (p[0]/3.0)*n;
      Fi1 = (p[1]/3.0)*n;
      Fi2 = (p[2]/3.0)*n;
      break;
    case 2:
      double c1 = 1.0/6.0, c2 = 1.0/12.0;
      Fi0 = (c1*p[0]+c2*(p[1]+p[2]))*n;
      Fi1 = (c1*p[1]+c2*(p[2]+p[0]))*n;
      Fi2 = (c1*p[2]+c2*(p[0]+p[1]))*n;
      break;
  }     
  Fv = 0.0;
}

//------------------------------------------------------------------------------

// Included (MB)
void PostFcnEuler::computeDerivativeOfForce(double dp1dxj[4][3], double ddp1dxj[4][3], double *Xface[3], double *dXface[3],
                                            Vec3D &n, Vec3D &dn, double d2w[3], double *Vwall, double *dVwall,
											double *Vface[3],  double *dVface[3], double *Vtet[4], double *dVtet[4], double dS[3],
											double *pin,  Vec3D &dFi0, Vec3D &dFi1, Vec3D &dFi2, Vec3D &dFv, double dPdx[3][3], double ddPdx[3][3], int hydro)
{

  double pcg[3], p[3];
  double dPcg[3], dP[3];
  double pcgin;
  double dPcgin;
  int i;

  double dPinfty = dpinfty * dS[0];

  if (hydro == 0) {
    for(i=0;i<3;i++) {
      pcg[i] = varFcn->getPressure(Vface[i]);
      dPcg[i] = varFcn->getPressure(dVface[i]);
    }
  } 
  else if (hydro == 1){ // hydrostatic pressure
     for(i=0;i<3;i++) {
        pcg[i]  = varFcn->hydrostaticPressure(Vface[i][0],Xface[i]);
	dPcg[i] = varFcn->DerivativeHydrostaticPressure(dVface[i][0], Vface[i][0], Xface[i], dXface[i]);
     }
  }else if (hydro == 2){ // hydrodynamic pressure
    for (i=0; i<3; i++) 
    {
      pcg[i]  = varFcn->hydrodynamicPressure(Vface[i], Xface[i]);
      dPcg[i] = varFcn->getPressure(dVface[i]) 
              - varFcn->DerivativeHydrostaticPressure(dVface[i][0], Vface[i][0], Xface[i], dXface[i]);
    }

  }

  if (pin)
    pcgin = *pin;
  else
    if (hydro == 0)
      pcgin = pinfty;
    else
      pcgin = 0.0;

  if (pin) {
    if (dimFlag)
      dPcgin = (-2.0*(*pin)/mach) * dS[0];
    else
      dPcgin = dPin * dS[0];
  }
  else {
    if (hydro == 0)
      dPcgin = dPinfty;
    else
      dPcgin = 0.0;
  }

  p[0] = (pcg[0] - pcgin) ;
  p[1] = (pcg[1] - pcgin) ;
  p[2] = (pcg[2] - pcgin) ;

  dP[0] = (dPcg[0] - dPcgin) ;
  dP[1] = (dPcg[1] - dPcgin) ;
  dP[2] = (dPcg[2] - dPcgin) ;

 Vec3D x0, x1, x2, x;
 Vec3D dx0, dx1, dx2, dx;
 double temp, dTemp;
 for(int i = 0; i<3; i++)
 {
  x0[i] = Xface[0][i];
  x1[i] = Xface[1][i];
  x2[i] = Xface[2][i];
  dx0[i] = dXface[0][i];
  dx1[i] = dXface[1][i];
  dx2[i] = dXface[2][i];
 }

// for node 0
 i=0; // i represents node
 x = (7.0/18.0)*(x1 + x2 - 2.0*x0);
 dx = (7.0/18.0)*(dx1 + dx2 - 2.0*dx0);
 temp = 2.0*p[i] + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dTemp = 2.0*dP[i] + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];
 dFi0 = (1.0/6.0*dTemp)*n + (1.0/6.0*temp)*dn;

// for node 1
 i=1;
 x = (7.0/18.0)*(x2 + x0 - 2.0*x1);
 dx = (7.0/18.0)*(dx2 + dx0 - 2.0*dx1);
 temp = 2.0*p[i] + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dTemp = 2.0*dP[i] + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];
 dFi1 = (1.0/6.0*dTemp)*n + (1.0/6.0*temp)*dn;

// for node 2
 i=2;
 x = (7.0/18.0)*(x0 + x1 - 2.0*x2);
 dx = (7.0/18.0)*(dx0 + dx1 - 2.0*dx2);
 temp = 2.0*p[i] + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dTemp = 2.0*dP[i] + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];
 dFi2 = (1.0/6.0*dTemp)*n + (1.0/6.0*temp)*dn;

  dFv = 0.0;

}
//-----------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void PostFcnEuler::computeDerivativeOfForce2(double dp1dxj[4][3], double ddp1dxj[4][3], double *Xface[3], double *dXface[3],
                                             Vec3D &n, Vec3D &dn, double d2w[3], double *Vwall, double *dVwall,
                                             double *Vface[3],  double *dVface[3], double *Vtet[4], double *dVtet[4], double dS[3],
                                             double *pin,  Vec3D &dFi0, Vec3D &dFi1, Vec3D &dFi2, Vec3D &dFv, double dPdx[3][3], double ddPdx[3][3], int hydro)
{
  double pcg[3], p[3];
  double dPcg[3];
  double pcgin;
  double dPcgin;
  int i;

  double dPinfty = dpinfty * dS[0];
  double dPinftydS[3] = {0};
  dPinftydS[0] = dpinfty;
  double dPcgdVface[3][5] = {0};

  if (hydro == 0) {
    for(i=0;i<3;i++) {
      pcg[i] = varFcn->getPressure(Vface[i]);
      varFcn->computedPdV(dPcgdVface[i]);
      dPcg[i] = 0.0;
      for(int l=0; l<5; ++l) {
        dPcg[i] += dPcgdVface[i][l]*dVface[i][l];
      }
    }
  } else if (hydro == 1) { // hydrostatic pressure
    for(i=0;i<3;i++) {
      pcg[i]  = varFcn->hydrostaticPressure(Vface[i][0],Xface[i]);
      dPcg[i] = varFcn->DerivativeHydrostaticPressure(dVface[i][0], Vface[i][0], Xface[i], dXface[i]);
    }
  } else if (hydro == 2){ // hydrodynamic pressure
    for (i=0; i<3; i++) {
      pcg[i]  = varFcn->hydrodynamicPressure(Vface[i], Xface[i]);
      dPcg[i] = varFcn->getPressure(dVface[i])
              - varFcn->DerivativeHydrostaticPressure(dVface[i][0], Vface[i][0], Xface[i], dXface[i]);
    }
  }

  if (pin)
    pcgin = *pin;
  else
    if (hydro == 0)
      pcgin = pinfty;
    else
      pcgin = 0.0;

  double dPcgindS[3] = {0};
  double dPcgindPinfty(0);
  if (pin) {
    if (dimFlag) {
//      dPcgin = (-2.0*(*pin)/mach) * dS[0];
      dPcgindS[0] = -2.0*(*pin)/mach;
     } else {
//      dPcgin = dPin * dS[0];
      dPcgindS[0] = dPin;
     }
  }
  else {
    if (hydro == 0) {
      dPcgin = dPinfty;
      dPcgindPinfty = 1.0;
    } else
      dPcgin = 0.0;
  }
  for(int j=0; j<3; ++j) {
    dPcgindS[j] += dPcgindPinfty*dPinftydS[j];
    dPcgin += dPcgindS[j]*dS[j];
  }

  p[0] = (pcg[0] - pcgin) ;
  p[1] = (pcg[1] - pcgin) ;
  p[2] = (pcg[2] - pcgin) ;

 Vec3D xC, xS, xP, x[3];
 Vec3D temp = 0.0;;
 for(int i = 0; i<3; i++) {
  xC[i] = Xface[0][i];
  xS[i] = Xface[1][i];
  xP[i] = Xface[2][i];
 }

 x[0] = (7.0/18.0)*(xS + xP - 2.0*xC);
 x[1] = (7.0/18.0)*(xP + xC - 2.0*xS);
 x[2] = (7.0/18.0)*(xC + xS - 2.0*xP);

 double dx0dXface0[3] = {0};   double dx0dXface1[3] = {0};   double dx0dXface2[3] = {0};
 double dx1dXface0[3] = {0};   double dx1dXface1[3] = {0};   double dx1dXface2[3] = {0};
 double dx2dXface0[3] = {0};   double dx2dXface1[3] = {0};   double dx2dXface2[3] = {0};

 for(int i=0; i<3; ++i) {
   dx0dXface0[i] = -7.0/9.0;   dx0dXface1[i] = 7.0/18.0;   dx0dXface2[i] = 7.0/18.0;
   dx1dXface0[i] = 7.0/18.0;   dx1dXface1[i] = -7.0/9.0;   dx1dXface2[i] = 7.0/18.0;
   dx2dXface0[i] = 7.0/18.0;   dx2dXface1[i] = 7.0/18.0;   dx2dXface2[i] = -7.0/9.0;
 }

 double dTempdPcg[3] = {0}, dTempdPcgin[3] = {0};
 double dTempdx[3][3] = {0}, dTempddPdx[3][3] = {0};
 for(int i=0; i<3; ++i) {
   dTempdPcg[i] = 2.0;
   dTempdPcgin[i] = -2.0;
   temp[i] += 2.0*p[i];
   for(int j=0; j<3; ++j) {
     dTempdx[i][j] = dPdx[i][j];
     dTempddPdx[i][j] = x[i][j];
     temp[i] += dPdx[i][j]*x[i][j];
   }
 }

 double dFidTemp[3] = {0}, dFi0dn[3] = {0}, dFi1dn[3] = {0}, dFi2dn[3] = {0};

 for(int i=0; i<3; ++i) {
   dFidTemp[i] = 1.0/6.0*n[i];
   dFi0dn[i] = (1.0/6.0*temp[0]);
   dFi1dn[i] = (1.0/6.0*temp[1]);
   dFi2dn[i] = (1.0/6.0*temp[2]);
 }

 dFi0 = 0.0;  dFi1 = 0.0;  dFi2 = 0.0;

 double dFi0dPcg0[3] = {0}, dFi1dPcg1[3] = {0}, dFi2dPcg2[3] = {0};
 double dFi0dPcgin[3] = {0}, dFi1dPcgin[3] = {0}, dFi2dPcgin[3] = {0};
 double dFi0ddPdx[3][3] = {0}, dFi1ddPdx[3][3] = {0}, dFi2ddPdx[3][3] = {0};
 double dFi0dXface0[3][3] = {0}, dFi0dXface1[3][3] = {0}, dFi0dXface2[3][3] = {0};
 double dFi1dXface0[3][3] = {0}, dFi1dXface1[3][3] = {0}, dFi1dXface2[3][3] = {0};
 double dFi2dXface0[3][3] = {0}, dFi2dXface1[3][3] = {0}, dFi2dXface2[3][3] = {0};
 double dFi0dS[3][3] = {0}, dFi1dS[3][3] = {0}, dFi2dS[3][3] = {0};
 double dFi0dVface[3][5] = {0}, dFi1dVface[3][5] = {0}, dFi2dVface[3][5] = {0};
 for(int i=0; i<3; ++i) {
   dFi0dPcg0[i] = dFidTemp[i]*dTempdPcg[0];
   dFi1dPcg1[i] = dFidTemp[i]*dTempdPcg[1];
   dFi2dPcg2[i] = dFidTemp[i]*dTempdPcg[2];
   dFi0dPcgin[i] = dFidTemp[i]*dTempdPcgin[0];
   dFi1dPcgin[i] = dFidTemp[i]*dTempdPcgin[1];
   dFi2dPcgin[i] = dFidTemp[i]*dTempdPcgin[2];

   dFi0[i] += dFi0dn[i]*dn[i];
   dFi1[i] += dFi1dn[i]*dn[i];
   dFi2[i] += dFi2dn[i]*dn[i];

   for(int l=0; l<5; ++l) {
     dFi0dVface[i][l] = dFi0dPcg0[i]*dPcgdVface[0][l];
     dFi1dVface[i][l] = dFi1dPcg1[i]*dPcgdVface[1][l];
     dFi2dVface[i][l] = dFi2dPcg2[i]*dPcgdVface[2][l];

     dFi0[i] += dFi0dVface[i][l]*dVface[0][l];
     dFi1[i] += dFi1dVface[i][l]*dVface[1][l];
     dFi2[i] += dFi2dVface[i][l]*dVface[2][l];
   }
   for(int j=0; j<3; ++j) {
     dFi0dS[i][j] = dFi0dPcgin[i]*dPcgindS[j];
     dFi1dS[i][j] = dFi1dPcgin[i]*dPcgindS[j];
     dFi2dS[i][j] = dFi2dPcgin[i]*dPcgindS[j];
     dFi0ddPdx[i][j] = dFidTemp[i]*dTempddPdx[0][j];
     dFi1ddPdx[i][j] = dFidTemp[i]*dTempddPdx[1][j];
     dFi2ddPdx[i][j] = dFidTemp[i]*dTempddPdx[2][j];
     dFi0dXface0[i][j] = dFidTemp[i]*dTempdx[0][j]*dx0dXface0[j];
     dFi0dXface1[i][j] = dFidTemp[i]*dTempdx[0][j]*dx0dXface1[j];
     dFi0dXface2[i][j] = dFidTemp[i]*dTempdx[0][j]*dx0dXface2[j];
     dFi1dXface0[i][j] = dFidTemp[i]*dTempdx[1][j]*dx1dXface0[j];
     dFi1dXface1[i][j] = dFidTemp[i]*dTempdx[1][j]*dx1dXface1[j];
     dFi1dXface2[i][j] = dFidTemp[i]*dTempdx[1][j]*dx1dXface2[j];
     dFi2dXface0[i][j] = dFidTemp[i]*dTempdx[2][j]*dx2dXface0[j];
     dFi2dXface1[i][j] = dFidTemp[i]*dTempdx[2][j]*dx2dXface1[j];
     dFi2dXface2[i][j] = dFidTemp[i]*dTempdx[2][j]*dx2dXface2[j];

     dFi0[i] += dFi0ddPdx[i][j]*ddPdx[0][j] + dFi0dXface0[i][j]*dXface[0][j] + dFi0dXface1[i][j]*dXface[1][j] + dFi0dXface2[i][j]*dXface[2][j] + dFi0dS[i][j]*dS[j];
     dFi1[i] += dFi1ddPdx[i][j]*ddPdx[1][j] + dFi1dXface0[i][j]*dXface[0][j] + dFi1dXface1[i][j]*dXface[1][j] + dFi1dXface2[i][j]*dXface[2][j] + dFi1dS[i][j]*dS[j];
     dFi2[i] += dFi2ddPdx[i][j]*ddPdx[2][j] + dFi2dXface0[i][j]*dXface[0][j] + dFi2dXface1[i][j]*dXface[1][j] + dFi2dXface2[i][j]*dXface[2][j] + dFi2dS[i][j]*dS[j];
   }
 }

 dFv = 0.0;

}

//------------------------------------------------------------------------------

void PostFcnEuler::computeDerivativeOperatorsOfForce(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double *Vface[3], double *Vtet[4], double *pin,
                                                     double dPdx[3][3], int hydro, double dFi0dn[3], double dFi1dn[3], double dFi2dn[3],
                                                     double dFi0ddPdx[3][3], double dFi1ddPdx[3][3], double dFi2ddPdx[3][3],
                                                     double dFi0dXface0[3][3], double dFi0dXface1[3][3], double dFi0dXface2[3][3],
                                                     double dFi1dXface0[3][3], double dFi1dXface1[3][3], double dFi1dXface2[3][3],
                                                     double dFi2dXface0[3][3], double dFi2dXface1[3][3], double dFi2dXface2[3][3],
                                                     double dFi0dS[3][3], double dFi1dS[3][3], double dFi2dS[3][3],
                                                     double dFi0dVface[3][5], double dFi1dVface[3][5], double dFi2dVface[3][5],
                                                     double dFvddp1dxj[3][4][3], double dFvdn[3][3], double dFvdV[3][4][5])
{
  double pcg[3], p[3];
  double pcgin;
  int i;

  double dPinftydS[3] = {0};
  dPinftydS[0] = dpinfty;
  double dPcgdVface[3][5] = {0};

  if (hydro == 0) {
    for(i=0;i<3;i++) {
      pcg[i] = varFcn->getPressure(Vface[i]);
      varFcn->computedPdV(dPcgdVface[i]);
    }
  } else if (hydro == 1) { // hydrostatic pressure
  } else if (hydro == 2){ // hydrodynamic pressure
  }

  if (pin)
    pcgin = *pin;
  else
    if (hydro == 0)
      pcgin = pinfty;
    else
      pcgin = 0.0;

  double dPcgindS[3] = {0};
  double dPcgindPinfty(0);
  if (pin) {
    if (dimFlag) {
//      dPcgin = (-2.0*(*pin)/mach) * dS[0];
      dPcgindS[0] = -2.0*(*pin)/mach;
     } else {
//      dPcgin = dPin * dS[0];
      dPcgindS[0] = dPin;
     }
  }
  else {
    if (hydro == 0) {
      dPcgindPinfty = 1.0;
    }
  }
  for(int j=0; j<3; ++j) {
    dPcgindS[j] += dPcgindPinfty*dPinftydS[j];
  }

  p[0] = (pcg[0] - pcgin) ;
  p[1] = (pcg[1] - pcgin) ;
  p[2] = (pcg[2] - pcgin) ;

 Vec3D xC, xS, xP, x[3];
 Vec3D temp = 0.0;;
 for(int i = 0; i<3; i++) {
  xC[i] = Xface[0][i];
  xS[i] = Xface[1][i];
  xP[i] = Xface[2][i];
 }

 x[0] = (7.0/18.0)*(xS + xP - 2.0*xC);
 x[1] = (7.0/18.0)*(xP + xC - 2.0*xS);
 x[2] = (7.0/18.0)*(xC + xS - 2.0*xP);

 double dx0dXface0[3] = {0};   double dx0dXface1[3] = {0};   double dx0dXface2[3] = {0};
 double dx1dXface0[3] = {0};   double dx1dXface1[3] = {0};   double dx1dXface2[3] = {0};
 double dx2dXface0[3] = {0};   double dx2dXface1[3] = {0};   double dx2dXface2[3] = {0};

 for(int i=0; i<3; ++i) {
   dx0dXface0[i] = -7.0/9.0;   dx0dXface1[i] = 7.0/18.0;   dx0dXface2[i] = 7.0/18.0;
   dx1dXface0[i] = 7.0/18.0;   dx1dXface1[i] = -7.0/9.0;   dx1dXface2[i] = 7.0/18.0;
   dx2dXface0[i] = 7.0/18.0;   dx2dXface1[i] = 7.0/18.0;   dx2dXface2[i] = -7.0/9.0;
 }

 double dTempdPcg[3] = {0}, dTempdPcgin[3] = {0};
 double dTempdx[3][3] = {0}, dTempddPdx[3][3] = {0};
 for(int i=0; i<3; ++i) {
   dTempdPcg[i] = 2.0;
   dTempdPcgin[i] = -2.0;
   temp[i] += 2.0*p[i];
   for(int j=0; j<3; ++j) {
     dTempdx[i][j] = dPdx[i][j];
     dTempddPdx[i][j] = x[i][j];
     temp[i] += dPdx[i][j]*x[i][j];
   }
 }

 double dFidTemp[3] = {0};

 for(int i=0; i<3; ++i) {
   dFidTemp[i] = 1.0/6.0*n[i];
   dFi0dn[i] = (1.0/6.0*temp[0]);
   dFi1dn[i] = (1.0/6.0*temp[1]);
   dFi2dn[i] = (1.0/6.0*temp[2]);
 }

 double dFi0dPcg0[3] = {0}, dFi1dPcg1[3] = {0}, dFi2dPcg2[3] = {0};
 double dFi0dPcgin[3] = {0}, dFi1dPcgin[3] = {0}, dFi2dPcgin[3] = {0};
 for(int i=0; i<3; ++i) {
   dFi0dPcg0[i] = dFidTemp[i]*dTempdPcg[0];
   dFi1dPcg1[i] = dFidTemp[i]*dTempdPcg[1];
   dFi2dPcg2[i] = dFidTemp[i]*dTempdPcg[2];
   dFi0dPcgin[i] = dFidTemp[i]*dTempdPcgin[0];
   dFi1dPcgin[i] = dFidTemp[i]*dTempdPcgin[1];
   dFi2dPcgin[i] = dFidTemp[i]*dTempdPcgin[2];

   for(int l=0; l<5; ++l) {
     dFi0dVface[i][l] = dFi0dPcg0[i]*dPcgdVface[0][l];
     dFi1dVface[i][l] = dFi1dPcg1[i]*dPcgdVface[1][l];
     dFi2dVface[i][l] = dFi2dPcg2[i]*dPcgdVface[2][l];
   }
   for(int j=0; j<3; ++j) {
     dFi0dS[i][j] = dFi0dPcgin[i]*dPcgindS[j];
     dFi1dS[i][j] = dFi1dPcgin[i]*dPcgindS[j];
     dFi2dS[i][j] = dFi2dPcgin[i]*dPcgindS[j];
     dFi0ddPdx[i][j] = dFidTemp[i]*dTempddPdx[0][j];
     dFi1ddPdx[i][j] = dFidTemp[i]*dTempddPdx[1][j];
     dFi2ddPdx[i][j] = dFidTemp[i]*dTempddPdx[2][j];
     dFi0dXface0[i][j] = dFidTemp[i]*dTempdx[0][j]*dx0dXface0[j];
     dFi0dXface1[i][j] = dFidTemp[i]*dTempdx[0][j]*dx0dXface1[j];
     dFi0dXface2[i][j] = dFidTemp[i]*dTempdx[0][j]*dx0dXface2[j];
     dFi1dXface0[i][j] = dFidTemp[i]*dTempdx[1][j]*dx1dXface0[j];
     dFi1dXface1[i][j] = dFidTemp[i]*dTempdx[1][j]*dx1dXface1[j];
     dFi1dXface2[i][j] = dFidTemp[i]*dTempdx[1][j]*dx1dXface2[j];
     dFi2dXface0[i][j] = dFidTemp[i]*dTempdx[2][j]*dx2dXface0[j];
     dFi2dXface1[i][j] = dFidTemp[i]*dTempdx[2][j]*dx2dXface1[j];
     dFi2dXface2[i][j] = dFidTemp[i]*dTempdx[2][j]*dx2dXface2[j];
   }
 }

}


void PostFcnEuler::computeForceTransmitted(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double d2w[3],
					   double *Vwall, double *Vface[3], double *Vtet[4],
					   double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, double dPdx[3][3], int hydro,int fid)
{

  double pcg[3], p[3];
  double pcgin;
  int i;

  if (hydro == 0) {
    for(i=0;i<3;i++)
      pcg[i] = varFcn->getPressure(Vface[i],fid);
  }
  else if (hydro == 1){ // hydrostatic pressure
     for(i=0;i<3;i++)
        pcg[i] = varFcn->hydrostaticPressure(Vface[i][0],Xface[i]);
  }else if (hydro == 2){ // hydrodynamic pressure
    for (i=0; i<3; i++)
      pcg[i] = varFcn->hydrodynamicPressure(Vface[i],Xface[i]);

  }

  if (pin)
    pcgin = *pin;
  else
    if (hydro == 0)
      pcgin = pinfty;
    else
      pcgin = 0.0;

  p[0] = (pcg[0] - pcgin) ;
  p[1] = (pcg[1] - pcgin) ;
  p[2] = (pcg[2] - pcgin) ;

// ##########################

 Vec3D xC, xS, xP, x1, x2, x3, x0, x;
 double p_1C, p_1S, p_2S, p_2P, p_3P, p_3C;
 double p_C, p_S, p_P, p_0C, p_0S, p_0P;

 // C to 0, S to 1 P to 2
 for(int i = 0; i<3; i++)
 {
  xC[i] = Xface[0][i];
  xS[i] = Xface[1][i];
  xP[i] = Xface[2][i];
 }

 x1 = (xC+xS)/2.0;
 x2 = (xS+xP)/2.0;
 x3 = (xP+xC)/2.0;
 x0 = (xC+xS+xP)/3.0;

 p_C = p[0], p_S = p[1], p_P = p[2];
 // Computing p_1C, p_0C and p_3C
 i = 0;
 x = x1-xC;
 p_1C = p_C + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 x = x0-xC;
 p_0C =  p_C + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 x = x3-xC;
 p_3C =  p_C + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];

 // Computing p_1S, p_0S and p_2S
 i = 1;
 x = x1-xS;
 p_1S = p_S + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 x = x0-xS;
 p_0S =  p_S + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 x = x2-xS;
 p_2S =  p_S + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];

 // Computing p_2P, p_0P and p_3P
 i = 2;
 x = x2-xP;
 p_2P = p_P + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 x = x0-xP;
 p_0P =  p_P + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 x = x3-xP;
 p_3P =  p_P + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];

 // computation of pressure flux
 double phi_C, phi_S, phi_P, phi_1, phi_2, phi_3, phi_0;
 phi_C = (2.0*p_C + p_0C + p_3C)/12.0 + (2.0*p_C + p_0C + p_1C)/12.0 ;
 phi_S = (2.0*p_S + p_0S + p_1S)/12.0 + (2.0*p_S + p_0S + p_2S)/12.0 ;
 phi_P = (2.0*p_P + p_0P + p_2P)/12.0 + (2.0*p_P + p_0P + p_3P)/12.0 ;

 phi_1 = (2.0*p_1C + p_0C + p_C)/12.0 + (2.0*p_1S + p_0S + p_S)/12.0 ;
 phi_2 = (2.0*p_2S + p_0S + p_S)/12.0 + (2.0*p_2P + p_0P + p_P)/12.0 ;
 phi_3 = (2.0*p_3P + p_0P + p_P)/12.0 + (2.0*p_3C + p_0C + p_C)/12.0 ;

 phi_0 = (2.0*p_0C + p_1C + p_C)/12.0 + (2.0*p_0C + p_3C + p_C)/12.0 + (2.0*p_0S + p_1S + p_S)/12.0 + (2.0*p_0S + p_2S + p_S)/12.0 + (2.0*p_0P + p_2P + p_P)/12.0 + (2.0*p_0P + p_3P + p_P)/12.0 ;

 Fi0 = (phi_C + phi_1/2.0 + phi_3/2.0 + phi_0/3.0) * (1.0/6.0*n);
 Fi1 = (phi_S + phi_1/2.0 + phi_2/2.0 + phi_0/3.0) * (1.0/6.0*n);
 Fi2 = (phi_P + phi_2/2.0 + phi_3/2.0 + phi_0/3.0) * (1.0/6.0*n);

 Fv = 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
void PostFcnEuler::computeDerivativeOfForceTransmitted(double dp1dxj[4][3], double ddp1dxj[4][3], double *Xface[3], double *dXface[3],
                                            Vec3D &n, Vec3D &dn, double d2w[3], double *Vwall, double *dVwall,
					    double *Vface[3],  double *dVface[3], double *Vtet[4], double *dVtet[4], double dS[3],  
					    double *pin,  Vec3D &dFi0, Vec3D &dFi1, Vec3D &dFi2, Vec3D &dFv, double dPdx[3][3], double ddPdx[3][3], int hydro)
{

  double pcg[3], p[3];
  double dPcg[3], dP[3];
  double pcgin;
  double dPcgin;
  int i;

  double dPinfty = dpinfty * dS[0];

  if (hydro == 0) {
    for(i=0;i<3;i++) {
      pcg[i] = varFcn->getPressure(Vface[i]);
      dPcg[i] = varFcn->getPressure(dVface[i]);
    }
  } 
  else if (hydro == 1){ // hydrostatic pressure
     for(i=0;i<3;i++) {
        pcg[i]  = varFcn->hydrostaticPressure(Vface[i][0],Xface[i]);
        dPcg[i] = varFcn->DerivativeHydrostaticPressure(dVface[i][0], Vface[i][0], Xface[i], dXface[i]);
     }
  }else if (hydro == 2){ // hydrodynamic pressure
    for (i=0; i<3; i++) 
    {
      pcg[i]  = varFcn->hydrodynamicPressure(Vface[i], Xface[i]);
      dPcg[i] = varFcn->getPressure(dVface[i]) 
              - varFcn->DerivativeHydrostaticPressure(dVface[i][0], Vface[i][0], Xface[i], dXface[i]);
    }

  }

  if (pin)
    pcgin = *pin;
  else
    if (hydro == 0)
      pcgin = pinfty;
    else
      pcgin = 0.0;

  if (pin) {
    if (dimFlag)
      dPcgin = (-2.0*(*pin)/mach) * dS[0];
    else
      dPcgin = dPin * dS[0];
  }
  else {
    if (hydro == 0)
      dPcgin = dPinfty;
    else
      dPcgin = 0.0;
  }

  p[0] = (pcg[0] - pcgin) ;
  p[1] = (pcg[1] - pcgin) ;
  p[2] = (pcg[2] - pcgin) ;

  dP[0] = (dPcg[0] - dPcgin) ;
  dP[1] = (dPcg[1] - dPcgin) ;
  dP[2] = (dPcg[2] - dPcgin) ;

 Vec3D xC, xS, xP, x1, x2, x3, x0, x;
 Vec3D dxC, dxS, dxP, dx1, dx2, dx3, dx0, dx;
 double p_1C, p_1S, p_2S, p_2P, p_3P, p_3C;
 double dP_1C, dP_1S, dP_2S, dP_2P, dP_3P, dP_3C;
 double p_C, p_S, p_P, p_0C, p_0S, p_0P;
 double dP_C, dP_S, dP_P, dP_0C, dP_0S, dP_0P;

 // C to 0, S to 1 P to 2
 for(int i = 0; i<3; i++)
 {
  xC[i] = Xface[0][i];
  xS[i] = Xface[1][i];
  xP[i] = Xface[2][i];
  dxC[i] = dXface[0][i];
  dxS[i] = dXface[1][i];
  dxP[i] = dXface[2][i];
 }

 x1 = (xC+xS)/2.0;
 x2 = (xS+xP)/2.0;
 x3 = (xP+xC)/2.0;
 x0 = (xC+xS+xP)/3.0;
 dx1 = (dxC+dxS)/2.0;
 dx2 = (dxS+dxP)/2.0;
 dx3 = (dxP+dxC)/2.0;
 dx0 = (dxC+dxS+dxP)/3.0;

 p_C = p[0], p_S = p[1], p_P = p[2];
 dP_C = dP[0], dP_S = dP[1], dP_P = dP[2];
 // Computing p_1C, p_0C and p_3C, and computing dP_1C, dP_0C and dP_3C
 i = 0;
 x = x1-xC;
 dx = dx1-dxC;
 p_1C = p_C + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dP_1C = dP_C + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];
 x = x0-xC;
 dx = dx0-dxC;
 p_0C =  p_C + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dP_0C =  dP_C + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];
 x = x3-xC;
 dx = dx3-dxC;
 p_3C =  p_C + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dP_3C =  dP_C + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];

 // Computing p_1S, p_0S and p_2S, and computing dP_1S, dP_0S and dP_2S
 i = 1;
 x = x1-xS;
 dx = dx1-dxS;
 p_1S = p_S + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dP_1S = dP_S + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];
 x = x0-xS;
 dx = dx0-dxS;
 p_0S =  p_S + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dP_0S =  dP_S + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];
 x = x2-xS;
 dx = dx2-dxS;
 p_2S =  p_S + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dP_2S =  dP_S + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];

 // Computing p_2P, p_0P and p_3P, and computing dP_2P, dP_0P and dP_3P
 i = 2;
 x = x2-xP;
 dx = dx2-dxP;
 p_2P = p_P + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dP_2P = dP_P + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];
 x = x0-xP;
 dx = dx0-dxP;
 p_0P =  p_P + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dP_0P =  dP_P + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];
 x = x3-xP;
 dx = dx3-dxP;
 p_3P =  p_P + dPdx[i][0]*x[0] + dPdx[i][1]*x[1] + dPdx[i][2]*x[2];
 dP_3P =  dP_P + ddPdx[i][0]*x[0] + dPdx[i][0]*dx[0] + ddPdx[i][1]*x[1] + dPdx[i][1]*dx[1] + ddPdx[i][2]*x[2] + dPdx[i][2]*dx[2];

 // computation of pressure flux
 double phi_C, phi_S, phi_P, phi_1, phi_2, phi_3, phi_0;
 double dPhi_C, dPhi_S, dPhi_P, dPhi_1, dPhi_2, dPhi_3, dPhi_0;
 phi_C = (2.0*p_C + p_0C + p_3C)/12.0 + (2.0*p_C + p_0C + p_1C)/12.0 ;
 dPhi_C = (2.0*dP_C + dP_0C + dP_3C)/12.0 + (2.0*dP_C + dP_0C + dP_1C)/12.0 ;
 phi_S = (2.0*p_S + p_0S + p_1S)/12.0 + (2.0*p_S + p_0S + p_2S)/12.0 ;
 dPhi_S = (2.0*dP_S + dP_0S + dP_1S)/12.0 + (2.0*dP_S + dP_0S + dP_2S)/12.0 ;
 phi_P = (2.0*p_P + p_0P + p_2P)/12.0 + (2.0*p_P + p_0P + p_3P)/12.0 ;
 dPhi_P = (2.0*dP_P + dP_0P + dP_2P)/12.0 + (2.0*dP_P + dP_0P + dP_3P)/12.0 ;

 phi_1 = (2.0*p_1C + p_0C + p_C)/12.0 + (2.0*p_1S + p_0S + p_S)/12.0 ;
 dPhi_1 = (2.0*dP_1C + dP_0C + dP_C)/12.0 + (2.0*dP_1S + dP_0S + dP_S)/12.0 ;
 phi_2 = (2.0*p_2S + p_0S + p_S)/12.0 + (2.0*p_2P + p_0P + p_P)/12.0 ;
 dPhi_2 = (2.0*dP_2S + dP_0S + dP_S)/12.0 + (2.0*dP_2P + dP_0P + dP_P)/12.0 ;
 phi_3 = (2.0*p_3P + p_0P + p_P)/12.0 + (2.0*p_3C + p_0C + p_C)/12.0 ;
 dPhi_3 = (2.0*dP_3P + dP_0P + dP_P)/12.0 + (2.0*dP_3C + dP_0C + dP_C)/12.0 ;

 phi_0 = (2.0*p_0C + p_1C + p_C)/12.0 + (2.0*p_0C + p_3C + p_C)/12.0 + (2.0*p_0S + p_1S + p_S)/12.0 + (2.0*p_0S + p_2S + p_S)/12.0 + (2.0*p_0P + p_2P + p_P)/12.0 + (2.0*p_0P + p_3P + p_P)/12.0 ;
 dPhi_0 = (2.0*dP_0C + dP_1C + dP_C)/12.0 + (2.0*dP_0C + dP_3C + dP_C)/12.0 + (2.0*dP_0S + dP_1S + dP_S)/12.0 + (2.0*dP_0S + dP_2S + dP_S)/12.0 + (2.0*dP_0P + dP_2P + dP_P)/12.0 + (2.0*dP_0P + dP_3P + dP_P)/12.0 ;

 dFi0 = (dPhi_C + dPhi_1/2.0 + dPhi_3/2.0 + dPhi_0/3.0) * (1.0/6.0*n) + (phi_C + phi_1/2.0 + phi_3/2.0 + phi_0/3.0) * (1.0/6.0*dn);
 dFi1 = (dPhi_S + dPhi_1/2.0 + dPhi_2/2.0 + dPhi_0/3.0) * (1.0/6.0*n) + (phi_S + phi_1/2.0 + phi_2/2.0 + phi_0/3.0) * (1.0/6.0*dn);
 dFi2 = (dPhi_P + dPhi_2/2.0 + dPhi_3/2.0 + dPhi_0/3.0) * (1.0/6.0*n) + (phi_P + phi_2/2.0 + phi_3/2.0 + phi_0/3.0) * (1.0/6.0*dn);

 dFv = 0.0;

}

//------------------------------------------------------------------------------

// Included (YC)

// Included (YC)
void PostFcnEuler::computeDerivativeOperatorsOfForceTransmitted(double dp1dxj[4][3],
                                                                double *Xface[3], Vec3D &n,
                                                                double *Vface[3], double *Vtet[4], double *pin,
                                                                double dPdx[3][3], int hydro,
                                                                double dFi0dn[3], double dFi0dS[3][3], double dFi0dVface[3][3][5],
                                                                double dFi0ddPdx[3][3][3], double dFi0dXface[3][3][3],
                                                                double dFi1dn[3], double dFi1dS[3][3], double dFi1dVface[3][3][5],
                                                                double dFi1ddPdx[3][3][3], double dFi1dXface[3][3][3],
                                                                double dFi2dn[3], double dFi2dS[3][3], double dFi2dVface[3][3][5],
																double dFi2ddPdx[3][3][3], double dFi2dXface[3][3][3],
																double dFvddp1dxj[3][4][3], double dFvdn[3][3], double dFvdV[3][4][5]
                                                               )
{

  double pcg[3], p[3];
  double dPcg[3];
  double pcgin;
  double dPcgin;
  double cf[500];
  int i;

  double dPinftydS[3] = {0};
  dPinftydS[0] = dpinfty;
  double dPcgdVface[3][5] = {0};

  if (hydro == 0) {
    for(i=0;i<3;i++) {
      pcg[i] = varFcn->getPressure(Vface[i]);
      varFcn->computedPdV(dPcgdVface[i]);
      dPcg[i] = 0.0;
    }
  } 
  else if (hydro == 1){ // hydrostatic pressure
    fprintf(stderr, " *** Error: PostFcnEuler::computeDerivativeOperatorsOfForceTransmitted for hydro = 1 is not implemented\n");
  } else if (hydro == 2){ // hydrodynamic pressure
    fprintf(stderr, " *** Error: PostFcnEuler::computeDerivativeOperatorsOfForceTransmitted for hydro = 2 is not implemented\n");
  }

  if (pin) pcgin = *pin;
  else {
    if (hydro == 0) pcgin = pinfty;
    else pcgin = 0.0;
  }

  double dPcgindS[3] = {0};
  double dPcgindPinfty(0);
  if (pin) {
    if (dimFlag) dPcgindS[0] = -2.0*(*pin)/mach;
    else dPcgindS[0] = dPin;
  } else {
    if (hydro == 0) dPcgindPinfty = 1.0;
  }
  dPcgin = 0.0;
  for(int j=0; j<3; ++j) dPcgindS[j] += dPcgindPinfty*dPinftydS[j]; 

  p[0] = (pcg[0] - pcgin) ;
  p[1] = (pcg[1] - pcgin) ;
  p[2] = (pcg[2] - pcgin) ;

 Vec3D xC, xS, xP, x1, x2, x3, x0, x;
 double p_1C, p_1S, p_2S, p_2P, p_3P, p_3C;
 double p_C, p_S, p_P, p_0C, p_0S, p_0P;

 // C to 0, S to 1 P to 2
 for(int i = 0; i<3; i++)
 {
  xC[i] = Xface[0][i];
  xS[i] = Xface[1][i];
  xP[i] = Xface[2][i];
 }

 x1 = (xC+xS)/2.0;
 x2 = (xS+xP)/2.0;
 x3 = (xP+xC)/2.0;
 x0 = (xC+xS+xP)/3.0;

 p_C = p[0], p_S = p[1], p_P = p[2];
 // Computing p_1C, p_0C and p_3C, and computing dP_1C, dP_0C and dP_3C
 p_1C = p_C 
      + dPdx[0][0]*(x1[0]-xC[0]) 
      + dPdx[0][1]*(x1[1]-xC[1]) 
      + dPdx[0][2]*(x1[2]-xC[2]);
 p_0C =  p_C 
      + dPdx[0][0]*(x0[0]-xC[0]) 
      + dPdx[0][1]*(x0[1]-xC[1]) 
      + dPdx[0][2]*(x0[2]-xC[2]);
 p_3C =  p_C 
      + dPdx[0][0]*(x3[0]-xC[0]) 
      + dPdx[0][1]*(x3[1]-xC[1]) 
      + dPdx[0][2]*(x3[2]-xC[2]);

 // Computing p_1S, p_0S and p_2S, and computing dP_1S, dP_0S and dP_2S
 p_1S = p_S 
      + dPdx[1][0]*(x1[0]-xS[0]) 
      + dPdx[1][1]*(x1[1]-xS[1]) 
      + dPdx[1][2]*(x1[2]-xS[2]);
 p_0S =  p_S 
      + dPdx[1][0]*(x0[0]-xS[0]) 
      + dPdx[1][1]*(x0[1]-xS[1]) 
      + dPdx[1][2]*(x0[2]-xS[2]);
 p_2S =  p_S 
      + dPdx[1][0]*(x2[0]-xS[0]) 
      + dPdx[1][1]*(x2[1]-xS[1]) 
      + dPdx[1][2]*(x2[2]-xS[2]);

 // Computing p_2P, p_0P and p_3P, and computing dP_2P, dP_0P and dP_3P
 p_2P = p_P 
      + dPdx[2][0]*(x2[0]-xP[0]) 
      + dPdx[2][1]*(x2[1]-xP[1]) 
      + dPdx[2][2]*(x2[2]-xP[2]);
 p_0P =  p_P 
      + dPdx[2][0]*(x0[0]-xP[0]) 
      + dPdx[2][1]*(x0[1]-xP[1]) 
      + dPdx[2][2]*(x0[2]-xP[2]);
 p_3P =  p_P 
      + dPdx[2][0]*(x3[0]-xP[0]) 
      + dPdx[2][1]*(x3[1]-xP[1]) 
      + dPdx[2][2]*(x3[2]-xP[2]);

 // computation of pressure flux
 double phi_C, phi_S, phi_P, phi_1, phi_2, phi_3, phi_0;
 phi_C = (2.0*p_C + p_0C + p_3C)/12.0 + (2.0*p_C + p_0C + p_1C)/12.0;
 cf[0] = 2.0*(x0[0]-xC[0])+x3[0]-xC[0]+x1[0]-xC[0];
 cf[1] = 2.0*(x0[1]-xC[1])+(x3[1]-xC[1])+(x1[1]-xC[1]);
 cf[2] = 2.0*(x0[2]-xC[2])+(x3[2]-xC[2])+(x1[2]-xC[2]);
 cf[3] = 2.0/3.0*dPdx[0][0]+0.5*dPdx[0][0];
 cf[4] = 2.0/3.0*dPdx[0][1]+0.5*dPdx[0][1];
 cf[5] = 2.0/3.0*dPdx[0][2]+0.5*dPdx[0][2];
 cf[6] = 2.0/3.0*dPdx[0][0]+0.5*dPdx[0][0];
 cf[7] = 2.0/3.0*dPdx[0][1]+0.5*dPdx[0][1];
 cf[8] = 2.0/3.0*dPdx[0][2]+0.5*dPdx[0][2];
 cf[9] = -4.0/3.0*dPdx[0][0] - dPdx[0][0];
 cf[10] = -4.0/3.0*dPdx[0][1] - dPdx[0][1];
 cf[11] = - 4.0/3.0*dPdx[0][2] - dPdx[0][2];
 phi_S = (2.0*p_S + p_0S + p_1S)/12.0 + (2.0*p_S + p_0S + p_2S)/12.0;
 cf[12] = 2.0*(x0[0]-xS[0]) + (x1[0]-xS[0]) + (x2[0]-xS[0]);
 cf[13] = 2.0*(x0[1]-xS[1]) + (x1[1]-xS[1]) + (x2[1]-xS[1]);
 cf[14] = 2.0*(x0[2]-xS[2]) + (x1[2]-xS[2]) + (x2[2]-xS[2]);
 cf[15] = 2.0/3.0*dPdx[1][0] + 0.5*dPdx[1][0];
 cf[16] = 2.0/3.0*dPdx[1][1] + 0.5*dPdx[1][1];
 cf[17] = 2.0/3.0*dPdx[1][2] + 0.5*dPdx[1][2];
 cf[18] = 2.0/3.0*dPdx[1][0] + 0.5*dPdx[1][0];
 cf[19] = 2.0/3.0*dPdx[1][1] + 0.5*dPdx[1][1];
 cf[20] = 2.0/3.0*dPdx[1][2] + 0.5*dPdx[1][2];
 cf[21] = - 4.0/3.0*dPdx[1][0] - dPdx[1][0];
 cf[22] = - 4.0/3.0*dPdx[1][1] - dPdx[1][1];
 cf[23] = - 4.0/3.0*dPdx[1][2] - dPdx[1][2];
 phi_P = (2.0*p_P + p_0P + p_2P)/12.0 + (2.0*p_P + p_0P + p_3P)/12.0;
 cf[24] = 2.0*(x0[0]-xP[0]) + (x2[0]-xP[0]) + (x3[0]-xP[0]);
 cf[25] = 2.0*(x0[1]-xP[1]) + (x2[1]-xP[1]) + (x3[1]-xP[1]);
 cf[26] = (x3[2]-xP[2]) + (x2[2]-xP[2]) + 2.0*(x0[2]-xP[2]);
 cf[27] = 2.0/3.0*dPdx[2][0] + 0.5*dPdx[2][0];
 cf[28] = 2.0/3.0*dPdx[2][1] + 0.5*dPdx[2][1];
 cf[29] = 2.0/3.0*dPdx[2][2] + 0.5*dPdx[2][2];
 cf[30] = 2.0/3.0*dPdx[2][0] + 0.5*dPdx[2][0];
 cf[31] = 2.0/3.0*dPdx[2][1] + 0.5*dPdx[2][1];
 cf[32] = 2.0/3.0*dPdx[2][2] + 0.5*dPdx[2][2];
 cf[33] = - 4.0/3.0*dPdx[2][0] - dPdx[2][0];
 cf[34] = - 4.0/3.0*dPdx[2][1] - dPdx[2][1];
 cf[35] = - 4.0/3.0*dPdx[2][2] - dPdx[2][2]; 
 phi_1 = (2.0*p_1C + p_0C + p_C)/12.0 + (2.0*p_1S + p_0S + p_S)/12.0;
 cf[36] = 2.0*(x1[0]-xC[0]) + (x0[0]-xC[0]);
 cf[37] = (x0[1]-xC[1]) + 2.0*(x1[1]-xC[1]);
 cf[38] = 2.0*(x1[2]-xC[2]) + (x0[2]-xC[2]);
 cf[39] = 2.0*(x1[0]-xS[0]) + (x0[0]-xS[0]);
 cf[40] = (x0[1]-xS[1]) + 2.0*(x1[1]-xS[1]);
 cf[41] = (x0[2]-xS[2]) + 2.0*(x1[2]-xS[2]);
 cf[42] = 4.0/3.0*dPdx[0][0] - 5.0/3.0*dPdx[1][0];
 cf[43] = 4.0/3.0*dPdx[0][1] - 5.0/3.0*dPdx[1][1];
 cf[44] = 4.0/3.0*dPdx[0][2] - 5.0/3.0*dPdx[1][2];
 cf[45] = 4.0/3.0*dPdx[1][0] - 5.0/3.0*dPdx[0][0];
 cf[46] = 4.0/3.0*dPdx[1][1] - 5.0/3.0*dPdx[0][1];
 cf[47] = 4.0/3.0*dPdx[1][2] - 5.0/3.0*dPdx[0][2];
 cf[48] = 1.0/3.0*dPdx[0][0] + 1.0/3.0*dPdx[1][0];
 cf[49] = 1.0/3.0*dPdx[0][1] + 1.0/3.0*dPdx[1][1];
 cf[50] = 1.0/3.0*dPdx[0][2] + 1.0/3.0*dPdx[1][2]; 
 phi_2 = (2.0*p_2S + p_0S + p_S)/12.0 + (2.0*p_2P + p_0P + p_P)/12.0;
 cf[51] = 2.0*(x2[0]-xS[0]) + (x0[0]-xS[0]);
 cf[52] = 2.0*(x2[1]-xS[1]) + (x0[1]-xS[1]);
 cf[53] = 2.0*(x2[2]-xS[2]) + (x0[2]-xS[2]);
 cf[54] = 2.0*(x2[0]-xP[0]) + (x0[0]-xP[0]);
 cf[55] = 2.0*(x2[1]-xP[1]) + (x0[1]-xP[1]);
 cf[56] = 2.0*(x2[2]-xP[2]) + (x0[2]-xP[2]);
 cf[57] = 1.0/3.0*dPdx[1][0] + 1.0/3.0*dPdx[2][0];
 cf[58] = 1.0/3.0*dPdx[1][1] + 1.0/3.0*dPdx[2][1];
 cf[59] = 1.0/3.0*dPdx[1][2] + 1.0/3.0*dPdx[2][2];
 cf[60] = 4.0/3.0*dPdx[2][0] - 5.0/3.0*dPdx[1][0];
 cf[61] = 4.0/3.0*dPdx[2][1] - 5.0/3.0*dPdx[1][1];
 cf[62] = 4.0/3.0*dPdx[2][2] - 5.0/3.0*dPdx[1][2]; 
 cf[63] = 4.0/3.0*dPdx[1][0] - 5.0/3.0*dPdx[2][0];
 cf[64] = 4.0/3.0*dPdx[1][1] - 5.0/3.0*dPdx[2][1];
 cf[65] = 4.0/3.0*dPdx[1][2] - 5.0/3.0*dPdx[2][2];
 phi_3 = (2.0*p_3P + p_0P + p_P)/12.0 + (2.0*p_3C + p_0C + p_C)/12.0;
 cf[66] = (x0[0]-xC[0]) + 2.0*(x3[0]-xC[0]);
 cf[67] = 2.0*(x3[1]-xC[1]) + (x0[1]-xC[1]);
 cf[68] = 2.0*(x3[2]-xC[2]) + (x0[2]-xC[2]);
 cf[69] = (x0[0]-xP[0]) + 2.0*(x3[0]-xP[0]);
 cf[70] = (x0[1]-xP[1]) + 2.0*(x3[1]-xP[1]);
 cf[71] = (x0[2]-xP[2]) + 2.0*(x3[2]-xP[2]);
 cf[72] = 4.0/3.0*dPdx[2][0] - 5.0/3.0*dPdx[0][0];
 cf[73] = 4.0/3.0*dPdx[2][1] - 5.0/3.0*dPdx[0][1];
 cf[74] = 4.0/3.0*dPdx[2][2] - 5.0/3.0*dPdx[0][2];
 cf[75] = 1.0/3.0*dPdx[2][0] + 1.0/3.0*dPdx[0][0];
 cf[76] = 1.0/3.0*dPdx[0][1] + 1.0/3.0*dPdx[2][1];
 cf[77] = 1.0/3.0*dPdx[0][2] + 1.0/3.0*dPdx[2][2];
 cf[78] = 4.0/3.0*dPdx[0][0] - 5.0/3.0*dPdx[2][0];
 cf[79] = 4.0/3.0*dPdx[0][1] - 5.0/3.0*dPdx[2][1];
 cf[80] = 4.0/3.0*dPdx[0][2] - 5.0/3.0*dPdx[2][2]; 
 phi_0 = (2.0*p_0C + p_1C + p_C)/12.0 + (2.0*p_0C + p_3C + p_C)/12.0 + (2.0*p_0S + p_1S + p_S)/12.0 + (2.0*p_0S + p_2S + p_S)/12.0 + (2.0*p_0P + p_2P + p_P)/12.0 + (2.0*p_0P + p_3P + p_P)/12.0;
 cf[81] = (x3[0]-xC[0]) + (x1[0]-xC[0]) + 4.0*(x0[0]-xC[0]);
 cf[82] = 4.0*(x0[1]-xC[1]) + (x1[1]-xC[1]) + (x3[1]-xC[1]);
 cf[83] = 4.0*(x0[2]-xC[2]) + (x1[2]-xC[2]) + (x3[2]-xC[2]);
 cf[84] = (x1[0]-xS[0]) + (x2[0]-xS[0]) + 4.0*(x0[0]-xS[0]);
 cf[85] = 4.0*(x0[1]-xS[1]) + (x2[1]-xS[1]) + (x1[1]-xS[1]);
 cf[86] = 4.0*(x0[2]-xS[2]) + (x2[2]-xS[2]) + (x1[2]-xS[2]);
 cf[87] = (x2[0]-xP[0]) + (x3[0]-xP[0]) + 4.0*(x0[0]-xP[0]);
 cf[88] = 4.0*(x0[1]-xP[1]) + (x3[1]-xP[1]) + (x2[1]-xP[1]);
 cf[89] = 4.0*(x0[2]-xP[2]) + (x3[2]-xP[2]) + (x2[2]-xP[2]);
 cf[90] = - 11.0/3.0*dPdx[0][0]+ 11.0/6.0*dPdx[1][0]+ 11.0/6.0*dPdx[2][0];
 cf[91] = - 11.0/3.0*dPdx[0][1] + 11.0/6.0*dPdx[1][1] + 11.0/6.0*dPdx[2][1];
 cf[92] = 11.0/6.0*dPdx[2][2] - 11.0/3.0*dPdx[0][2] + 11.0/6.0*dPdx[1][2];
 cf[93] = - 11.0/3.0*dPdx[1][0] + 11.0/6.0*dPdx[2][0] + 11.0/6.0*dPdx[0][0];
 cf[94] = 11.0/6.0*dPdx[2][1] + 11.0/6.0*dPdx[0][1] - 11.0/3.0*dPdx[1][1];
 cf[95] = 11.0/6.0*dPdx[0][2] - 11.0/3.0*dPdx[1][2] + 11.0/6.0*dPdx[2][2];
 cf[96] = 11.0/6.0*dPdx[0][0] + 11.0/6.0*dPdx[1][0] - 11.0/3.0*dPdx[2][0];
 cf[97] = 11.0/6.0*dPdx[0][1] + 11.0/6.0*dPdx[1][1] - 11.0/3.0*dPdx[2][1];
 cf[98] = 11.0/6.0*dPdx[1][2] - 11.0/3.0*dPdx[2][2] + 11.0/6.0*dPdx[0][2]; 
 cf[99] = (phi_C + phi_1/2.0 + phi_3/2.0 + phi_0/3.0)/6.0;
 cf[100] = 1.0/12.0*cf[0] + 1.0/24.0*cf[36] + 1.0/24.0*cf[66] + 1.0/36.0*cf[81];
 cf[101] = 1.0/36.0*cf[82] + 1.0/24.0*cf[67] + 1.0/24.0*cf[37] + 1.0/12.0*cf[1];
 cf[102] = 1.0/36.0*cf[83] + 1.0/24.0*cf[68] + 1.0/24.0*cf[38] + 1.0/12.0*cf[2];
 cf[103] = 1.0/24.0*cf[39] + 1.0/36.0*cf[84];
 cf[104] = 1.0/24.0*cf[40] + 1.0/36.0*cf[85];
 cf[105] = 1.0/24.0*cf[41] + 1.0/36.0*cf[86];
 cf[106] = 1.0/36.0*cf[87] + 1.0/24.0*cf[69];
 cf[107] = 1.0/36.0*cf[88] + 1.0/24.0*cf[70];
 cf[108] = 1.0/36.0*cf[89] + 1.0/24.0*cf[71];
 cf[109] = 1.0/12.0*cf[9] + 1.0/24.0*cf[72] + 1.0/24.0*cf[45] + 1.0/36.0*cf[90];
 cf[110] = 1.0/12.0*cf[10] + 1.0/24.0*cf[46] + 1.0/24.0*cf[73] + 1.0/36.0*cf[91];
 cf[111] = 1.0/36.0*cf[92] + 1.0/24.0*cf[74] + 1.0/24.0*cf[47] + 1.0/12.0*cf[11];
 cf[112] = 1.0/12.0*cf[3] + 1.0/24.0*cf[42] + 1.0/24.0*cf[75] + 1.0/36.0*cf[93];
 cf[113] = 1.0/24.0*cf[76] + 1.0/24.0*cf[43] + 1.0/12.0*cf[4] + 1.0/36.0*cf[94];
 cf[114] = 1.0/36.0*cf[95] + 1.0/24.0*cf[77] + 1.0/24.0*cf[44] + 1.0/12.0*cf[5];
 cf[115] = 1.0/12.0*cf[6] + 1.0/24.0*cf[48] + 1.0/24.0*cf[78] + 1.0/36.0*cf[96];
 cf[116] = 1.0/24.0*cf[79] + 1.0/24.0*cf[49] + 1.0/12.0*cf[7] + 1.0/36.0*cf[97];
 cf[117] = 1.0/36.0*cf[98] + 1.0/24.0*cf[80] + 1.0/24.0*cf[50] + 1.0/12.0*cf[8];
 double dFi0dPcg[3][3] = {0};
 for(int j=0; j<3; ++j) {
   dFi0dPcg[j][0] = 11.0/54.0*n[j];
   dFi0dPcg[j][1] = 7.0/108.0*n[j];
   dFi0dPcg[j][2] = 7.0/108.0*n[j];
 }
 for(int j=0; j<3; ++j) 
   for(int k=0; k<3; ++k) 
     for(int l=0; l<5; ++l) 
       dFi0dVface[j][k][l] += dFi0dPcg[j][k]*dPcgdVface[k][l];


 double dFi0dPcgin[3] = {0};
 for(int j=0; j<3; ++j) {
   dFi0dPcgin[j] = -1.0/3.0*n[j];
 }
 for(int j=0; j<3; ++j) {
   dFi0ddPdx[j][0][0] = 1.0/6.0*n[j]*cf[100];
   dFi0ddPdx[j][0][1] = 1.0/6.0*n[j]*cf[101];
   dFi0ddPdx[j][0][2] = 1.0/6.0*n[j]*cf[102];
   dFi0ddPdx[j][1][0] = 1.0/6.0*n[j]*cf[103];
   dFi0ddPdx[j][1][1] = 1.0/6.0*n[j]*cf[104];
   dFi0ddPdx[j][1][2] = 1.0/6.0*n[j]*cf[105];
   dFi0ddPdx[j][2][0] = 1.0/6.0*n[j]*cf[106];
   dFi0ddPdx[j][2][1] = 1.0/6.0*n[j]*cf[107];
   dFi0ddPdx[j][2][2] = 1.0/6.0*n[j]*cf[108];
   dFi0dXface[j][0][0] = 1.0/6.0*n[j]*cf[109];
   dFi0dXface[j][0][1] = 1.0/6.0*n[j]*cf[110];
   dFi0dXface[j][0][2] = 1.0/6.0*n[j]*cf[111];
   dFi0dXface[j][1][0] = 1.0/6.0*n[j]*cf[112];
   dFi0dXface[j][1][1] = 1.0/6.0*n[j]*cf[113];
   dFi0dXface[j][1][2] = 1.0/6.0*n[j]*cf[114];
   dFi0dXface[j][2][0] = 1.0/6.0*n[j]*cf[115];
   dFi0dXface[j][2][1] = 1.0/6.0*n[j]*cf[116];
   dFi0dXface[j][2][2] = 1.0/6.0*n[j]*cf[117];
 }
 dFi0dn[0] = cf[99];  dFi0dn[1] = cf[99];  dFi0dn[2] = cf[99];
 for(int l=0; l<3; ++l) {
   for(int j=0; j<3; ++j) {
     dFi0dS[l][j] = dFi0dPcgin[l]*dPcgindS[j]; 
   }
 }
 cf[118] = (phi_S + phi_1/2.0 + phi_2/2.0 + phi_0/3.0)/6.0;
 cf[119] = (phi_P + phi_2/2.0 + phi_3/2.0 + phi_0/3.0)/6.0;
 cf[120] = 1.0/24.0*cf[36] + 1.0/36.0*cf[81];
 cf[121] = 1.0/24.0*cf[37] + 1.0/36.0*cf[82];
 cf[122] = 1.0/24.0*cf[38] + 1.0/36.0*cf[83];
 cf[123] = 1.0/36.0*cf[84] + 1.0/12.0*cf[12] + 1.0/24.0*cf[39] + 1.0/24.0*cf[51];
 cf[124] = 1.0/24.0*cf[52] + 1.0/36.0*cf[85] + 1.0/24.0*cf[40] + 1.0/12.0*cf[13];
 cf[125] = 1.0/12.0*cf[14] + 1.0/24.0*cf[53] + 1.0/36.0*cf[86] + 1.0/24.0*cf[41];
 cf[126] = 1.0/36.0*cf[87] + 1.0/24.0*cf[54];
 cf[127] = 1.0/36.0*cf[88] + 1.0/24.0*cf[55];
 cf[128] = 1.0/36.0*cf[89] + 1.0/24.0*cf[56];
 cf[129] = 1.0/12.0*cf[15] + 1.0/24.0*cf[45] + 1.0/24.0*cf[57] + 1.0/36.0*cf[90];
 cf[130] = 1.0/36.0*cf[91] + 1.0/24.0*cf[58] + 1.0/24.0*cf[46] + 1.0/12.0*cf[16];
 cf[131] = 1.0/36.0*cf[92] + 1.0/24.0*cf[59] + 1.0/24.0*cf[47] + 1.0/12.0*cf[17];
 cf[132] = 1.0/12.0*cf[21] + 1.0/24.0*cf[60] + 1.0/24.0*cf[42] + 1.0/36.0*cf[93];
 cf[133] = 1.0/36.0*cf[94] + 1.0/24.0*cf[61] + 1.0/24.0*cf[43] + 1.0/12.0*cf[22];
 cf[134] = 1.0/36.0*cf[95] + 1.0/24.0*cf[62] + 1.0/12.0*cf[23] + 1.0/24.0*cf[44];
 cf[135] = 1.0/24.0*cf[48] + 1.0/36.0*cf[96] + 1.0/12.0*cf[18] + 1.0/24.0*cf[63];
 cf[136] = 1.0/12.0*cf[19] + 1.0/24.0*cf[64] + 1.0/36.0*cf[97] + 1.0/24.0*cf[49];
 cf[137] = 1.0/12.0*cf[20] + 1.0/24.0*cf[50] + 1.0/24.0*cf[65] + 1.0/36.0*cf[98];

 double dFi1dPcg[3][3] = {0};
 for(int j=0; j<3; ++j) {
   dFi1dPcg[j][0] = 7.0/108.0*n[j];
   dFi1dPcg[j][1] = 11.0/54.0*n[j];
   dFi1dPcg[j][2] = 7.0/108.0*n[j];
 }
 for(int j=0; j<3; ++j) 
   for(int k=0; k<3; ++k) 
     for(int l=0; l<5; ++l) 
       dFi1dVface[j][k][l] += dFi1dPcg[j][k]*dPcgdVface[k][l];

 double dFi1dPcgin[3] = {0};
 for(int j=0; j<3; ++j) {
   dFi1dPcgin[j] = -1.0/3.0*n[j];
 }
 for(int j=0; j<3; ++j) {
   dFi1ddPdx[j][0][0] = 1.0/6.0*n[j]*cf[120];
   dFi1ddPdx[j][0][1] = 1.0/6.0*n[j]*cf[121];
   dFi1ddPdx[j][0][2] = 1.0/6.0*n[j]*cf[122];
   dFi1ddPdx[j][1][0] = 1.0/6.0*n[j]*cf[123];
   dFi1ddPdx[j][1][1] = 1.0/6.0*n[j]*cf[124];
   dFi1ddPdx[j][1][2] = 1.0/6.0*n[j]*cf[125];
   dFi1ddPdx[j][2][0] = 1.0/6.0*n[j]*cf[126];
   dFi1ddPdx[j][2][1] = 1.0/6.0*n[j]*cf[127];
   dFi1ddPdx[j][2][2] = 1.0/6.0*n[j]*cf[128];
   dFi1dXface[j][0][0] = 1.0/6.0*n[j]*cf[129];
   dFi1dXface[j][0][1] = 1.0/6.0*n[j]*cf[130];
   dFi1dXface[j][0][2] = 1.0/6.0*n[j]*cf[131];
   dFi1dXface[j][1][0] = 1.0/6.0*n[j]*cf[132];
   dFi1dXface[j][1][1] = 1.0/6.0*n[j]*cf[133];
   dFi1dXface[j][1][2] = 1.0/6.0*n[j]*cf[134];
   dFi1dXface[j][2][0] = 1.0/6.0*n[j]*cf[135];
   dFi1dXface[j][2][1] = 1.0/6.0*n[j]*cf[136];
   dFi1dXface[j][2][2] = 1.0/6.0*n[j]*cf[137];
 }
 dFi1dn[0] = cf[118];  dFi1dn[1] = cf[118];  dFi1dn[2] = cf[118];
 for(int l=0; l<3; ++l) {
   for(int j=0; j<3; ++j) {
     dFi1dS[l][j] = dFi1dPcgin[l]*dPcgindS[j];
   }
 }
 cf[138] = 1.0/24.0*cf[66] + 1.0/36.0*cf[81] ; 
 cf[139] = 1.0/36.0*cf[82] + 1.0/24.0*cf[67]; 
 cf[140] = 1.0/36.0*cf[83] + 1.0/24.0*cf[68]; 
 cf[141] = 1.0/36.0*cf[84] + 1.0/24.0*cf[51];
 cf[142] = 1.0/24.0*cf[52] + 1.0/36.0*cf[85];
 cf[143] = 1.0/24.0*cf[53] + 1.0/36.0*cf[86];
 cf[144] = 1.0/24.0*cf[54] + 1.0/36.0*cf[87]  + 1.0/12.0*cf[24] + 1.0/24.0*cf[69];
 cf[145] = 1.0/24.0*cf[55] + 1.0/36.0*cf[88] + 1.0/12.0*cf[25] + 1.0/24.0*cf[70];
 cf[146] = 1.0/24.0*cf[56] + 1.0/36.0*cf[89] + 1.0/24.0*cf[71] + 1.0/12.0*cf[26];
 cf[147] = 1.0/12.0*cf[27] + 1.0/24.0*cf[57] + 1.0/36.0*cf[90] + 1.0/24.0*cf[72];
 cf[148] = 1.0/24.0*cf[73] + 1.0/36.0*cf[91] + 1.0/24.0*cf[58] + 1.0/12.0*cf[28];
 cf[149] = 1.0/24.0*cf[74] + 1.0/36.0*cf[92] + 1.0/24.0*cf[59] + 1.0/12.0*cf[29];
 cf[150] = 1.0/12.0*cf[30] + 1.0/24.0*cf[60] + 1.0/24.0*cf[75] + 1.0/36.0*cf[93];
 cf[151] = 1.0/36.0*cf[94] + 1.0/24.0*cf[76] + 1.0/24.0*cf[61] + 1.0/12.0*cf[31];
 cf[152] = 1.0/36.0*cf[95] + 1.0/24.0*cf[77] + 1.0/24.0*cf[62] + 1.0/12.0*cf[32];
 cf[153] = 1.0/12.0*cf[33] + 1.0/24.0*cf[63] + 1.0/24.0*cf[78] + 1.0/36.0*cf[96];
 cf[154] = 1.0/12.0*cf[34] + 1.0/24.0*cf[64] + 1.0/36.0*cf[97] + 1.0/24.0*cf[79];
 cf[155] = 1.0/12.0*cf[35] + 1.0/24.0*cf[65] + 1.0/24.0*cf[80] + 1.0/36.0*cf[98]; 

 double dFi2dPcg[3][3] = {0};
 for(int j=0; j<3; ++j) {
   dFi2dPcg[j][0] = 7.0/108.0*n[j];
   dFi2dPcg[j][1] = 7.0/108.0*n[j];
   dFi2dPcg[j][2] = 11.0/54.0*n[j];
 }
 for(int j=0; j<3; ++j) 
   for(int k=0; k<3; ++k) 
     for(int l=0; l<5; ++l) 
       dFi2dVface[j][k][l] += dFi2dPcg[j][k]*dPcgdVface[k][l];

 double dFi2dPcgin[3] = {0};
 for(int j=0; j<3; ++j) {
   dFi2dPcgin[j] = -1.0/3.0*n[j];
 }
 for(int j=0; j<3; ++j) {
   dFi2ddPdx[j][0][0] = 1.0/6.0*n[j]*cf[138];
   dFi2ddPdx[j][0][1] = 1.0/6.0*n[j]*cf[139];
   dFi2ddPdx[j][0][2] = 1.0/6.0*n[j]*cf[140];
   dFi2ddPdx[j][1][0] = 1.0/6.0*n[j]*cf[141];
   dFi2ddPdx[j][1][1] = 1.0/6.0*n[j]*cf[142];
   dFi2ddPdx[j][1][2] = 1.0/6.0*n[j]*cf[143];
   dFi2ddPdx[j][2][0] = 1.0/6.0*n[j]*cf[144];
   dFi2ddPdx[j][2][1] = 1.0/6.0*n[j]*cf[145];
   dFi2ddPdx[j][2][2] = 1.0/6.0*n[j]*cf[146];
   dFi2dXface[j][0][0] = 1.0/6.0*n[j]*cf[147];
   dFi2dXface[j][0][1] = 1.0/6.0*n[j]*cf[148];
   dFi2dXface[j][0][2] = 1.0/6.0*n[j]*cf[149];
   dFi2dXface[j][1][0] = 1.0/6.0*n[j]*cf[150];
   dFi2dXface[j][1][1] = 1.0/6.0*n[j]*cf[151];
   dFi2dXface[j][1][2] = 1.0/6.0*n[j]*cf[152];
   dFi2dXface[j][2][0] = 1.0/6.0*n[j]*cf[153];
   dFi2dXface[j][2][1] = 1.0/6.0*n[j]*cf[154];
   dFi2dXface[j][2][2] = 1.0/6.0*n[j]*cf[155];
 }
 dFi2dn[0] = cf[119];  dFi2dn[1] = cf[119];  dFi2dn[2] = cf[119];
 for(int l=0; l<3; ++l) {
   for(int j=0; j<3; ++j) {
     dFi2dS[l][j] = dFi2dPcgin[l]*dPcgindS[j];
   }
 }

}

//--------------------------------------------------------------------------------------

double PostFcnEuler::computeHeatPower(double dp1dxj[4][3], Vec3D& n, double d2w[3], 
				      double* Vwall, double* Vface[3], double* Vtet[4])
{

  return 0.0;

}

//--------------------------------------------------------------------------------------

double PostFcnEuler::computeHeatFluxRelatedValues(double dp1dxj[4][3], Vec3D& n, double d2w[3],
                                               double* Vwall, double* Vface[3], double* Vtet[4], bool includeKappa)
{

  return 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcnEuler::computeDerivativeOfHeatPower(double dp1dxj[4][3], double ddp1dxj[4][3], Vec3D& n, Vec3D& dn, double d2w[3], 
				      double* Vwall, double* dVwall, double* Vface[3], double* dVface[3], double* Vtet[4], double* dVtet[4], double dS[3])
{

  return 0.0;

}

//------------------------------------------------------------------------------

double PostFcnEuler::computeInterfaceWork(double dp1dxj[4][3], Vec3D& n, double ndot, 
					  double d2w[3], double* Vwall, double* Vface[3], 
					  double* Vtet[4], double pin)
{
  double p = third * ( varFcn->getPressure(Vface[0]) + varFcn->getPressure(Vface[1]) +
		       varFcn->getPressure(Vface[2]) ) - pin;
  double W = - ndot * p;

  return W;
}

//------------------------------------------------------------------------------

PostFcnNS::PostFcnNS(IoData &iod, VarFcn *vf) 
  : PostFcnEuler(iod, vf), NavierStokesTerm(iod, vf)
{

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = new WallFcn(iod, PostFcn::varFcn, viscoFcn);
  else
    wallFcn = 0;

}

//------------------------------------------------------------------------------

PostFcnNS::~PostFcnNS()
{
  if (wallFcn) delete wallFcn;
}

//------------------------------------------------------------------------------

// Included (MB)
void PostFcnNS::rstVar(IoData &iod, Communicator *com)
{

  PostFcnEuler::rstVar(iod, com);

  NavierStokesTerm::rstVarNS(iod, com);
  
  if (wallFcn)
    wallFcn->rstVar(iod, com);

}

//------------------------------------------------------------------------------

double PostFcnNS::computeFaceScalarQuantity(ScalarType type, double dp1dxj[4][3], 
					    Vec3D& n, double d2w[3], double* Vwall, 
					    double* Vface[3], double* Vtet[4])
{

  double q = 0.0;

  if (type == DELTA_PLUS) {
#if defined(HEAT_FLUX)
    q = computeHeatPower(dp1dxj, n, d2w, Vwall, Vface, Vtet) / sqrt(n*n);
#else
    if (wallFcn)
      q = wallFcn->computeDeltaPlus(n, d2w, Vwall, Vface);
    else
      fprintf(stderr, "*** Warning: yplus computation not implemented\n");
#endif
  }
  else if (type == SKIN_FRICTION) {
    Vec3D t(1.0, 0.0, 0.0);
    Vec3D F = computeViscousForce(dp1dxj, n, d2w, Vwall, Vface, Vtet);
    q = 2.0 * t * F / sqrt(n*n);
  }

  return q;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcnNS::computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType type, double dS[3], double *V, double *dV, double *X, double *dX, double phi)
{

  double q = 0.0;

  q = PostFcnEuler::computeDerivativeOfNodeScalarQuantity(type, dS, V, dV, X, dX);

  return q;

}

//------------------------------------------------------------------------------

Vec3D PostFcnNS::computeViscousForce(double dp1dxj[4][3], Vec3D& n, double d2w[3], 
				     double* Vwall, double* Vface[3], double* Vtet[4])
{

  Vec3D Fv;
	/*
  if (wallFcn && Vwall)
	{
    Fv = wallFcn->computeForce(n, d2w, Vwall, Vface);
	}
	else 
	{
	*/
    double u[4][3], ucg[3];
    computeVelocity(Vtet, u, ucg);

    double T[4], Tcg;
    computeTemperature(Vtet, T, Tcg);

    double dudxj[3][3];
    computeVelocityGradient(dp1dxj, u, dudxj);

    double mu     = viscoFcn->compute_mu(Tcg);
    double lambda = viscoFcn->compute_lambda(Tcg, mu);
    mu     *= ooreynolds_mu;
    lambda *= ooreynolds_mu;

    double tij[3][3];
    computeStressTensor(mu, lambda, dudxj, tij);

    Fv[0] = tij[0][0] * n[0] + tij[0][1] * n[1] + tij[0][2] * n[2];
    Fv[1] = tij[1][0] * n[0] + tij[1][1] * n[1] + tij[1][2] * n[2];
    Fv[2] = tij[2][0] * n[0] + tij[2][1] * n[1] + tij[2][2] * n[2];
//	}

  return -1.0 * Fv;

}

//------------------------------------------------------------------------------

Vec3D PostFcnNS::computeViscousForceCVBoundary(Vec3D& n,  double* Vi, double dudxj[3][3])
{

  Vec3D Fv;
  // Could be useful later
  /*
  if (wallFcn)
    Fv = wallFcn->computeForce(n, d2w, Vwall, Vface);
  else {
  */
  double T;
  computeTemperature(Vi,T);
  
  double mu     = viscoFcn->compute_mu(T);
  double lambda = viscoFcn->compute_lambda(T, mu);
  mu     *= ooreynolds_mu;
  lambda *= ooreynolds_mu;
  
  double tij[3][3];
  computeStressTensor(mu, lambda, dudxj, tij);
  
  Fv[0] = tij[0][0] * n[0] + tij[0][1] * n[1] + tij[0][2] * n[2];
  Fv[1] = tij[1][0] * n[0] + tij[1][1] * n[1] + tij[1][2] * n[2];
  Fv[2] = tij[2][0] * n[0] + tij[2][1] * n[1] + tij[2][2] * n[2];
  
  return -1.0 * Fv;
}

//------------------------------------------------------------------------------

// Included (MB)
Vec3D PostFcnNS::computeDerivativeOfViscousForce(double dp1dxj[4][3], double ddp1dxj[4][3], Vec3D& n, Vec3D& dn, double d2w[3],
				     double* Vwall, double* dVwall, double* Vface[3], double* dVface[3], double* Vtet[4], double* dVtet[4], double dS[3])
{

  Vec3D dFv;

  if (wallFcn)
    dFv = wallFcn->computeDerivativeOfForce(n, dn, d2w, Vwall, dVwall, Vface, dVface, dS[0]);
  else {
    double u[4][3], ucg[3];
    computeVelocity(Vtet, u, ucg);

    double du[4][3], ducg[3];
    computeDerivativeOfVelocity(dVtet, du, ducg);

    double T[4], Tcg;
    computeTemperature(Vtet, T, Tcg);

    double dT[4], dTcg;
    computeDerivativeOfTemperature(Vtet, dVtet, dT, dTcg);

    double dudxj[3][3];
    computeVelocityGradient(dp1dxj, u, dudxj);

    double ddudxj[3][3];
    computeDerivativeOfVelocityGradient(dp1dxj, ddp1dxj, u, du, ddudxj);

    double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dS[0];

    double mu     = viscoFcn->compute_mu(Tcg);
    double lambda = viscoFcn->compute_lambda(Tcg, mu);

    double dmu = dooreynolds_mu * mu + ooreynolds_mu * viscoFcn->compute_muDerivative(Tcg, dTcg, dS[0]);
    double dlambda = dooreynolds_mu * lambda + ooreynolds_mu * viscoFcn->compute_lambdaDerivative(mu, dmu, dS[0]);

    mu     *= ooreynolds_mu;
    lambda *= ooreynolds_mu;

    double tij[3][3];
    computeStressTensor(mu, lambda, dudxj, tij);

    double dtij[3][3];
    computeDerivativeOfStressTensor(mu, dmu, lambda, dlambda, dudxj, ddudxj, dtij);

    dFv[0] = dtij[0][0] * n[0] + dtij[0][1] * n[1] + dtij[0][2] * n[2] + tij[0][0] * dn[0] + tij[0][1] * dn[1] + tij[0][2] * dn[2];
    dFv[1] = dtij[1][0] * n[0] + dtij[1][1] * n[1] + dtij[1][2] * n[2] + tij[1][0] * dn[0] + tij[1][1] * dn[1] + tij[1][2] * dn[2];
    dFv[2] = dtij[2][0] * n[0] + dtij[2][1] * n[1] + dtij[2][2] * n[2] + tij[2][0] * dn[0] + tij[2][1] * dn[1] + tij[2][2] * dn[2];

  }

  return -1.0 * dFv;

}

//------------------------------------------------------------------------------


void PostFcnNS::computeDerivativeOperatorsOfViscousForce(double dp1dxj[4][3], Vec3D& n, double* Vtet[4],
				                      double dFvddp1dxj[3][4][3], double dFvdn[3][3], double dFvdVtet[3][4][5])
{ //YC

  Vec3D dFv;
  if (wallFcn) {
    fprintf(stderr, " *** Error: PostFcnNS::computeDerivativeOperatorsOfViscousForce is not implemented\n"); exit(-1);
  } else {
    double u[4][3], ucg[3];
    computeVelocity(Vtet, u, ucg);

    double dudVtet[4][3][4][4] = {0}, ducgdVtet[3][4][4] = {0};
    computeDerivativeOperatorsOfVelocity(dudVtet, ducgdVtet);

    double T[4], Tcg;
    computeTemperature(Vtet, T, Tcg);

    double dTdVtet[4][4][5] = {0}, dTcgdVtet[4][5] = {0};
    computeDerivativeOperatorsOfTemperature2(Vtet, dTdVtet, dTcgdVtet);

    double dudxj[3][3];
    computeVelocityGradient(dp1dxj, u, dudxj);

    double ddudxjddp1dxj[3][3][4][3] = {0}, ddudxjdu[3][3][4][3] = {0};
    computeDerivativeOperatorsOfVelocityGradient(dp1dxj, u, ddudxjddp1dxj, ddudxjdu);

    double mu     = viscoFcn->compute_mu(Tcg);
    double lambda = viscoFcn->compute_lambda(Tcg, mu);

    mu     *= ooreynolds_mu;
    lambda *= ooreynolds_mu;

    double tij[3][3];
    computeStressTensor(mu, lambda, dudxj, tij);

    double dtijddudxj[3][3][3][3] = {0};
    computeDerivativeOperatorsOfStressTensor(mu, lambda, dudxj, dtijddudxj, NULL, NULL);

    double dFvdtij[3][3][3] = {0};
    dFvdtij[2][2][0] = dFvdtij[1][1][0] = dFvdtij[0][0][0] = -n[0];
    dFvdtij[2][2][1] = dFvdtij[1][1][1] = dFvdtij[0][0][1] = -n[1];
    dFvdtij[2][2][2] = dFvdtij[1][1][2] = dFvdtij[0][0][2] = -n[2];
    dFvdn[0][0] = -tij[0][0];   dFvdn[0][1] = -tij[0][1];    dFvdn[0][2] = -tij[0][2];
    dFvdn[1][0] = -tij[1][0];   dFvdn[1][1] = -tij[1][1];    dFvdn[1][2] = -tij[1][2];
    dFvdn[2][0] = -tij[2][0];   dFvdn[2][1] = -tij[2][1];    dFvdn[2][2] = -tij[2][2];

    for(int i=0; i<3; ++i)
      for(int j=0; j<3; ++j)
        for(int k=0; k<3; ++k)
          for(int l=0; l<3; ++l)
            for(int m=0; m<3; ++m)
              for(int n=0; n<4; ++n)
                for(int o=0; o<3; ++o) {
                  dFvddp1dxj[i][n][o] += dFvdtij[i][j][k]*dtijddudxj[j][k][l][m]*ddudxjddp1dxj[l][m][n][o];
                  for(int p=0; p<4; ++p)
                    for(int q=0; q<4; ++q)
                      dFvdVtet[i][p][q] += dFvdtij[i][j][k]*dtijddudxj[j][k][l][m]*ddudxjdu[l][m][n][o]*dudVtet[n][o][p][q];
                }

  }

}

//------------------------------------------------------------------------------


void PostFcnNS::computeForce(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double d2w[3], 
												     double *Vwall, double *Vface[3], double *Vtet[4], 
					                   double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, double dPdx[3][3], int hydro, int fid)
{

  PostFcnEuler::computeForce(dp1dxj, Xface, n, d2w, Vwall, Vface, Vtet, pin, Fi0, Fi1, Fi2, Fv, dPdx, hydro,fid);

  Fv = computeViscousForce(dp1dxj, n, d2w, Vwall, Vface, Vtet);

}

//------------------------------------------------------------------------------

void PostFcnNS::computeForceEmbedded(int orderOfAccuracy,double dp1dxj[4][3], double *Xface[3], Vec3D &n, double d2w[3], 
				     double *Vwall, double *Vface[3], double *Vtet[4], 
				     double pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, double dPdx[3][3], 
                                     int hydro, int* fid, bool applyRealForce)
{

  PostFcnEuler::computeForceEmbedded(orderOfAccuracy,dp1dxj, Xface, n, d2w, Vwall, Vface, Vtet, pin, Fi0, Fi1, Fi2, Fv, dPdx, 
                                     hydro, fid, applyRealForce);

  if(applyRealForce)
    Fv = computeViscousForce(dp1dxj, n, d2w, Vwall, Vface, Vtet);
  else
    Fv = 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void PostFcnNS::computeDerivativeOfForce(double dp1dxj[4][3], double ddp1dxj[4][3], double *Xface[3], double *dXface[3],
                                            Vec3D &n, Vec3D &dn, double d2w[3], double *Vwall, double *dVwall,
                                            double *Vface[3], double *dVface[3], double *Vtet[4],
                                            double *dVtet[4], double dS[3], double *pin,  Vec3D &dFi0, Vec3D &dFi1, Vec3D &dFi2,
                                            Vec3D &dFv, double dPdx[3][3], double ddPdx[3][3], int hydro)
{

  PostFcnEuler::computeDerivativeOfForce(dp1dxj, ddp1dxj, Xface, dXface, n, dn, d2w, Vwall, dVwall, Vface, dVface, Vtet, dVtet, dS, pin, dFi0, dFi1, dFi2, dFv, dPdx,  ddPdx, hydro);

  dFv = computeDerivativeOfViscousForce(dp1dxj, ddp1dxj, n, dn, d2w, Vwall, dVwall, Vface, dVface, Vtet, dVtet, dS);

}

//------------------------------------------------------------------------------

void PostFcnNS::computeForceTransmitted(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double d2w[3],
                             double *Vwall, double *Vface[3], double *Vtet[4],
                    double *pin, Vec3D &Fi0, Vec3D &Fi1, Vec3D &Fi2, Vec3D &Fv, double dPdx[3][3], int hydro, int fid)
{

  PostFcnEuler::computeForceTransmitted(dp1dxj, Xface, n, d2w, Vwall, Vface, Vtet, pin, Fi0, Fi1, Fi2, Fv, dPdx, hydro,fid);

  Fv = computeViscousForce(dp1dxj, n, d2w, Vwall, Vface, Vtet);
}

//------------------------------------------------------------------------------

// Included (MB)
inline
void PostFcnNS::computeDerivativeOfForceTransmitted(double dp1dxj[4][3], double ddp1dxj[4][3], double *Xface[3], double *dXface[3],
                                            Vec3D &n, Vec3D &dn, double d2w[3], double *Vwall, double *dVwall,
                                            double *Vface[3], double *dVface[3], double *Vtet[4],
                                            double *dVtet[4], double dS[3], double *pin,  Vec3D &dFi0, Vec3D &dFi1, Vec3D &dFi2,
                                            Vec3D &dFv, double dPdx[3][3], double ddPdx[3][3], int hydro)
{

  PostFcnEuler::computeDerivativeOfForceTransmitted(dp1dxj, ddp1dxj, Xface, dXface, n, dn, d2w, Vwall, dVwall, Vface, dVface, Vtet, dVtet, dS, pin, dFi0, dFi1, dFi2, dFv, dPdx,  ddPdx, hydro);

  dFv = computeDerivativeOfViscousForce(dp1dxj, ddp1dxj, n, dn, d2w, Vwall, dVwall, Vface, dVface, Vtet, dVtet, dS);

}

//------------------------------------------------------------------------------

void PostFcnNS::computeDerivativeOperatorsOfForce(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double *Vface[3], double *Vtet[4], double *pin,
                                                  double dPdx[3][3], int hydro, double dFi0dn[3], double dFi1dn[3], double dFi2dn[3],
                                                  double dFi0ddPdx[3][3], double dFi1ddPdx[3][3], double dFi2ddPdx[3][3],
                                                  double dFi0dXface0[3][3], double dFi0dXface1[3][3], double dFi0dXface2[3][3],
                                                  double dFi1dXface0[3][3], double dFi1dXface1[3][3], double dFi1dXface2[3][3],
                                                  double dFi2dXface0[3][3], double dFi2dXface1[3][3], double dFi2dXface2[3][3],
                                                  double dFi0dS[3][3], double dFi1dS[3][3], double dFi2dS[3][3],
                                                  double dFi0dVface[3][5], double dFi1dVface[3][5], double dFi2dVface[3][5],
                                                  double dFvddp1dxj[3][4][3], double dFvdn[3][3], double dFvdV[3][4][5])
{

  PostFcnEuler::computeDerivativeOperatorsOfForce(dp1dxj, Xface, n, Vface, Vtet, pin, dPdx, hydro,
                                                  dFi0dn, dFi1dn, dFi2dn, dFi0ddPdx, dFi1ddPdx, dFi2ddPdx,
                                                  dFi0dXface0, dFi0dXface1, dFi0dXface2, dFi1dXface0, dFi1dXface1, dFi1dXface2,
                                                  dFi2dXface0, dFi2dXface1, dFi2dXface2, dFi0dS, dFi1dS, dFi2dS,
                                                  dFi0dVface, dFi1dVface, dFi2dVface, dFvddp1dxj, dFvdn, dFvdV);

  computeDerivativeOperatorsOfViscousForce(dp1dxj, n, Vtet, dFvddp1dxj, dFvdn, dFvdV);

}

//------------------------------------------------------------------------------

void PostFcnNS::computeDerivativeOperatorsOfForceTransmitted(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double *Vface[3], double *Vtet[4], double *pin,
                                                             double dPdx[3][3], int hydro, double dFi0dn[3], double dFi0dS[3][3], double dFi0dVface[3][3][5],
                                                             double dFi0ddPdx[3][3][3], double dFi0dXface[3][3][3], double dFi1dn[3], double dFi1dS[3][3], double dFi1dVface[3][3][5],
                                                             double dFi1ddPdx[3][3][3], double dFi1dXface[3][3][3], double dFi2dn[3], double dFi2dS[3][3], double dFi2dVface[3][3][5],
                                                             double dFi2ddPdx[3][3][3], double dFi2dXface[3][3][3], double dFvddp1dxj[3][4][3], double dFvdn[3][3], double dFvdV[3][4][5])
{
  PostFcnEuler::computeDerivativeOperatorsOfForceTransmitted(dp1dxj, Xface, n, Vface, Vtet, pin, dPdx, hydro, dFi0dn, dFi0dS, dFi0dVface, dFi0ddPdx, dFi0dXface, dFi1dn, dFi1dS, dFi1dVface, dFi1ddPdx, dFi1dXface, dFi2dn, dFi2dS, dFi2dVface, dFi2ddPdx, dFi2dXface, dFvddp1dxj, dFvdn, dFvdV);

  computeDerivativeOperatorsOfViscousForce(dp1dxj, n, Vtet, dFvddp1dxj, dFvdn, dFvdV);

}


//------------------------------------------------------------------------------

double PostFcnNS::computeHeatPower(double dp1dxj[4][3], Vec3D& n, double d2w[3], 
				   double* Vwall, double* Vface[3], double* Vtet[4])
{

  double hp = 0.0;

  if (wallFcn)
    hp = wallFcn->computeHeatPower(n, d2w, Vwall, Vface);
  else {
    double T[4], Tcg;
    computeTemperature(Vtet, T, Tcg);
    double dTdxj[3];
    computeTemperatureGradient(dp1dxj, T, dTdxj);
    double kappa = ooreynolds_mu * thermalCondFcn->compute(Tcg);
    double qj[3];
    NavierStokesTerm::computeHeatFluxVector(kappa, dTdxj, qj);
    hp = qj[0]*n[0] + qj[1]*n[1] + qj[2]*n[2]; 
}

  return hp;

}
//------------------------------------------------------------------------------
double PostFcnNS::computeHeatFluxRelatedValues(double dp1dxj[4][3], Vec3D& n, double d2w[3],
                                   double* Vwall, double* Vface[3], double* Vtet[4], bool includeKappa)
{
  double hp = 0.0;

  if (wallFcn)
    hp = wallFcn->computeHeatPower(n, d2w, Vwall, Vface);
  else {
    double T[4], Tcg;
    computeTemperature(Vtet, T, Tcg);
    double dTdxj[3];
    computeTemperatureGradient(dp1dxj, T, dTdxj);  

    double kappa = -1; //The fact that it is negative balances the minus sign in NavierStokesTerm::computeHeatFluxVector
    if(includeKappa == true)
      {
        kappa = ooreynolds_mu * thermalCondFcn->compute(Tcg);
      }
     
       double qj[3];
    NavierStokesTerm::computeHeatFluxVector(kappa, dTdxj, qj);
    hp = qj[0]*n[0] + qj[1]*n[1] + qj[2]*n[2]; 
    }
  return hp;

}


//------------------------------------------------------------------------------

// Included (MB)
double PostFcnNS::computeDerivativeOfHeatPower(double dp1dxj[4][3], double ddp1dxj[4][3], Vec3D& n, Vec3D& dn, double d2w[3], 
				   double* Vwall, double* dVwall, double* Vface[3], double* dVface[3], double* Vtet[4], double* dVtet[4], double dS[3])
{

  double dhp = 0.0;

  if (wallFcn)
    dhp = wallFcn->computeDerivativeOfHeatPower(n, dn, d2w, Vwall, dVwall, Vface, dVface, dS[0]);
  else {
    double T[4], Tcg;
    computeTemperature(Vtet, T, Tcg);
    double dT[4], dTcg;
    computeDerivativeOfTemperature(Vtet, dVtet, dT, dTcg);
    double dTdxj[3];
    computeTemperatureGradient(dp1dxj, T, dTdxj);
    double ddTdxj[3];
    computeDerivativeOfTemperatureGradient(dp1dxj, ddp1dxj, T, dT, ddTdxj);
    double kappa = ooreynolds * thermalCondFcn->compute(Tcg);
    double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dS[0];
    double dkappa = dooreynolds_mu * thermalCondFcn->compute(Tcg) + ooreynolds_mu * thermalCondFcn->computeDerivative(Tcg, dTcg, dS[0]);
    double qj[3];
    NavierStokesTerm::computeHeatFluxVector(kappa, dTdxj, qj);
    double dqj[3];
    computeDerivativeOfHeatFluxVector(kappa, dkappa, dTdxj, ddTdxj, dqj);
    dhp = dqj[0]*n[0] + qj[0]*dn[0] + dqj[1]*n[1] + qj[1]*dn[1] + dqj[2]*n[2] + qj[2]*dn[2]; 
  }

  return dhp;

}

//------------------------------------------------------------------------------

double PostFcnNS::computeInterfaceWork(double dp1dxj[4][3], Vec3D& n, double ndot, 
				       double d2w[3], double* Vwall, double* Vface[3], 
				       double* Vtet[4], double pin)
{

  double W = PostFcnEuler::computeInterfaceWork(dp1dxj, n, ndot, d2w, Vwall, Vface, Vtet, pin);

  if (wallFcn)
    W += wallFcn->computeInterfaceWork(n, d2w, Vwall, Vface);
  else {
    double u[4][3], ucg[3];
    computeVelocity(Vtet, u, ucg);

    double T[4], Tcg;
    computeTemperature(Vtet, T, Tcg);

    double dudxj[3][3];
    computeVelocityGradient(dp1dxj, u, dudxj);

    double mu     = viscoFcn->compute_mu(Tcg);
    double lambda = viscoFcn->compute_lambda(Tcg, mu);
    mu     *= ooreynolds_mu;
    lambda *= ooreynolds_mu;

    double tij[3][3];
    computeStressTensor(mu, lambda, dudxj, tij);

    W += (Vwall[1] * tij[0][0] + Vwall[2] * tij[1][0] + Vwall[3] * tij[2][0]) * n[0] +
      (Vwall[1] * tij[0][1] + Vwall[2] * tij[1][1] + Vwall[3] * tij[2][1]) * n[1] +
      (Vwall[1] * tij[0][2] + Vwall[2] * tij[1][2] + Vwall[3] * tij[2][2]) * n[2];
  }

  return W;

}

//------------------------------------------------------------------------------

PostFcnSA::PostFcnSA(IoData &iod, VarFcn *vf) : PostFcnNS(iod, vf), SATerm(iod)
{
}

//------------------------------------------------------------------------------

PostFcnSA::~PostFcnSA()
{
}
//------------------------------------------------------------------------------

// Included (MB)
void PostFcnSA::rstVar(IoData &iod, Communicator *com)
{

  PostFcnNS::rstVar(iod, com);

  rstVarSA(iod);

}

//------------------------------------------------------------------------------

double PostFcnSA::computeNodeScalarQuantity(ScalarType type, double *V, double *X,  int fluidId,double* phi)
{

  double q = 0.0;

  if (type == EDDY_VISCOSITY) {
    double T = PostFcn::varFcn->computeTemperature(V);
    double mul = viscoFcn->compute_mu(T);
    q = computeTurbulentViscosity(V, mul);
  }
  else
    q = PostFcnEuler::computeNodeScalarQuantity(type, V, X, fluidId);

  return q;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcnSA::computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType type, double dS[3], double *V, double *dV, double *X, double *dX, double phi)
{

  double dq = 0.0;

  if (type == DERIVATIVE_EDDY_VISCOSITY) {
    double T = PostFcn::varFcn->computeTemperature(V);
    double dT = PostFcn::varFcn->computeDerivativeOfTemperature(V, dV);
    double mul = viscoFcn->compute_mu(T);
    double dmul = viscoFcn->compute_muDerivative(T, dT, dS[0]);
    dq = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
  }
  else
  {
    dq = PostFcnEuler::computeDerivativeOfNodeScalarQuantity
         (type, dS, V, dV, X, dX);
  }

  return dq;

}

// Included (MB)
void PostFcnDES::rstVar(IoData &iod, Communicator *com)
{

  PostFcnNS::rstVar(iod, com);

  rstVarDES(iod);

}

//------------------------------------------------------------------------------
PostFcnDES::PostFcnDES(IoData &iod, VarFcn *vf) : PostFcnNS(iod, vf), DESTerm(iod)
{

}

//------------------------------------------------------------------------------

                                                                                                
double PostFcnDES::computeNodeScalarQuantity(ScalarType type, double *V, double *X,int fluidId,double* phi)
{

  double q = 0.0;

  if (type == EDDY_VISCOSITY) {
    double T = PostFcn::varFcn->computeTemperature(V);
    double mul = viscoFcn->compute_mu(T);
    q = computeTurbulentViscosity(V, mul);
  }
  else
    q = PostFcnEuler::computeNodeScalarQuantity(type, V, X, fluidId);

  return q;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcnDES::computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType type, double dS[3], double *V, double *dV, double *X, double *dX, double phi)
{

  double dq = 0.0;

  if (type == DERIVATIVE_EDDY_VISCOSITY) {
    double T = PostFcn::varFcn->computeTemperature(V);
    double dT = PostFcn::varFcn->computeDerivativeOfTemperature(V, dV);
    double mul = viscoFcn->compute_mu(T);
    double dmul = viscoFcn->compute_muDerivative(T, dT, dS[0]);
    dq = computeDerivativeOfTurbulentViscosity(V, dV, mul, dmul);
  }
  else
    dq = PostFcnEuler::computeDerivativeOfNodeScalarQuantity(type, dS, V, dV, X, dX);

  return dq;

}

//------------------------------------------------------------------------------

PostFcnKE::PostFcnKE(IoData &iod, VarFcn *vf) : PostFcnNS(iod, vf), KEpsilonTerm(iod)
{

}

//------------------------------------------------------------------------------

// Included (MB)
void PostFcnKE::rstVar(IoData &iod, Communicator *com)
{

  PostFcnNS::rstVar(iod, com);

  rstVarKE(iod);

}

//------------------------------------------------------------------------------

double PostFcnKE::computeNodeScalarQuantity(ScalarType type, double *V, double *X, int fluidId,double* phi)
{

  double q = 0.0;

  if (type == EDDY_VISCOSITY)
    q = computeTurbulentViscosity(V);
  else
    q = PostFcnEuler::computeNodeScalarQuantity(type, V, X, fluidId);

  return q;

}

//------------------------------------------------------------------------------

// Included (MB)
double PostFcnKE::computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType type, double dS[3], double *V, double *dV, double *X, double *dX, double phi)
{

  double dq = 0.0;

  if (type == DERIVATIVE_EDDY_VISCOSITY)
    dq = computeDerivativeOfTurbulentViscosity(V, dV, dS[0]);
  else
    dq = PostFcnEuler::computeDerivativeOfNodeScalarQuantity(type, dS, V, dV, X, dX);

  return dq;

}

//------------------------------------------------------------------------------
