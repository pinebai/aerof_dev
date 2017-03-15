#include <FluidShapeOptimizationHandler.h>

#include <IoData.h>
#include <Domain.h>
#include <GeoSource.h>
#include <DistVector.h>
#include <MeshMotionSolver.h>
#include <MatVecProd.h>
#include <KspPrec.h>
#include <KspSolver.h>
#include <MemoryPool.h>

#include <cmath>

#include <iostream>
#include <string>

#include "Dev/devtools.h"
//------------------------------------------------------------------------------

template<int dim>
FluidShapeOptimizationHandler<dim>::FluidShapeOptimizationHandler
(
  IoData &ioData,
  GeoSource &geoSource,
  Domain *dom//,
) :
step(-1),
actvar(-2),
numLocSub(-1),
reynolds0(-1.0),
kenergy0(-1.0),
length(-1.0),
surface(-1.0),
xmach(-1.0),
alprad(-1.0),
teta(-1.0),
ImplicitCoupledTsDesc<dim>(ioData, geoSource, dom),
Ap(NULL),
Am(NULL),
Xp(NULL),
Xm(NULL),
Lp(NULL),
Lm(NULL),
Z(NULL),
load(NULL),
dLoad(NULL),
dLoadref(NULL),
Fp(NULL),
Fm(NULL),
Up(NULL),
Um(NULL),
outFile(NULL),
domain(dom),
dXb(dom->getNodeDistInfo()),
dXdS(dom->getNodeDistInfo()),
dXdSb(dom->getNodeDistInfo()),
Xc(dom->getNodeDistInfo()),
dAdS(dom->getNodeDistInfo()),
dFdS(dom->getNodeDistInfo()),
dFdSref(dom->getNodeDistInfo()),
dUdS(dom->getNodeDistInfo()),
dfaU(dom->getNodeDistInfo()),
dfaX(dom->getNodeDistInfo()),
lambdaU(dom->getNodeDistInfo()),
lambdaSDisp(dom->getNodeDistInfo()),
lambdaX(dom->getNodeDistInfo()),
p(dom->getNodeDistInfo()),
dPdS(dom->getNodeDistInfo()),
Flux(dom->getNodeDistInfo()),
FluxFD(dom->getNodeDistInfo()),
Pin(dom->getFaceDistInfo()),
Uc(dom->getNodeDistInfo()),
// Tests
Xplus(dom->getNodeDistInfo()),
Xminus(dom->getNodeDistInfo()),
dX(dom->getNodeDistInfo()),
dddx(dom->getNodeDistInfo()),
dddy(dom->getNodeDistInfo()),
dddz(dom->getNodeDistInfo()),
dR(dom->getNodeDistInfo()),
dEdgeNorm(dom->getEdgeDistInfoMF()),
dFaceNorm(dom->getFaceNormDistInfo()),
dFaceNormVel(dom->getFaceNormDistInfo()),
dGradP(dom->getNodeDistInfo())
{
  // Initialize
  step = 0;

  if ( ioData.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_ ) {
    load = new DistSVec<double,3>(domain->getNodeDistInfo());
    dLoad = new DistSVec<double,3>(domain->getNodeDistInfo());
    dLoadref = new DistSVec<double,3>(domain->getNodeDistInfo());
  } else {
    load = 0;
    dLoad = 0;
    dLoadref = 0;
  }

  mms = 0;

//TODO this if-else conditions should be commented back in once all verification is removed
//  if ( ioData.sa.scFlag == SensitivityAnalysis::FINITEDIFFERENCE ) {
    Xp = new DistSVec<double,3>(domain->getNodeDistInfo());
    Xm = new DistSVec<double,3>(domain->getNodeDistInfo());
    Lp = new DistSVec<double,3>(domain->getNodeDistInfo());
    Lm = new DistSVec<double,3>(domain->getNodeDistInfo());
    Ap = new DistVec<double>(domain->getNodeDistInfo());
    Am = new DistVec<double>(domain->getNodeDistInfo());
    Fp = new DistSVec<double,dim>(domain->getNodeDistInfo());
    Fm = new DistSVec<double,dim>(domain->getNodeDistInfo());
    Up = new DistSVec<double,dim>(domain->getNodeDistInfo());
    Um = new DistSVec<double,dim>(domain->getNodeDistInfo());
/*  }
  else if ( ioData.sa.scFlag == SensitivityAnalysis::SEMIANALYTICAL ) {
    Xp = new DistSVec<double,3>(domain->getNodeDistInfo());
    Xm = new DistSVec<double,3>(domain->getNodeDistInfo());
    Lp = 0;
    Lm = 0;
    Ap = new DistVec<double>(domain->getNodeDistInfo());
    Am = new DistVec<double>(domain->getNodeDistInfo());
    Fp = new DistSVec<double,dim>(domain->getNodeDistInfo());
    Fm = new DistSVec<double,dim>(domain->getNodeDistInfo());
    Up = 0;
    Um = 0;
   }
  else {
    Xp = 0;
    Xm = 0;
    Lp = 0;
    Lm = 0;
    Ap = 0;
    Am = 0;
    Fp = 0;
    Fm = 0;
  } */

  Z = new DistSVec<double,3>(domain->getNodeDistInfo());

  if (ioData.sa.homotopy == SensitivityAnalysis::ON_HOMOTOPY)
  {
    if (ioData.sa.mvp == SensitivityAnalysis::H2)
    {
      mvp = new MatVecProdH2<dim,MatScalar,dim>(ioData, this->varFcn, this->timeState, this->spaceOp, domain, this->geoState);
    }
    else
    {
      mvp = new MatVecProdFD<dim,dim>(ioData.ts.implicit, this->timeState, this->geoState, this->spaceOp, domain, ioData);
    }
  }
  else
  {
    if (ioData.sa.mvp == SensitivityAnalysis::H2)
    {
      mvp = new MatVecProdH2<dim,MatScalar,dim>(ioData, this->varFcn, 0, this->spaceOp, domain, this->geoState);
    }
    else
    {
      mvp = new MatVecProdFD<dim,dim>(ioData.ts.implicit,  0, this->geoState, this->spaceOp, domain, ioData);
    }
  }

  pc = ImplicitTsDesc<dim>::template
    createPreconditioner<PrecScalar,dim>(ioData.sa.ksp.pc, domain);

  ksp = this->createKrylovSolver(this->getVecInfo(), ioData.sa.ksp, mvp, pc, this->com);
  dRdX = new MatVecProd_dRdX<dim,double,dim>(ioData, this->varFcn, this->timeState, this->spaceOp, domain, this->geoState);

  MemoryPool mp;

  mvp->exportMemory(&mp);
  pc->exportMemory(&mp);
//  dRdX->exportMemory(&mp);

  if (ioData.sa.sensMesh == SensitivityAnalysis::ON_SENSITIVITYMESH ||
      ioData.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_) {
    mms = new TetMeshMotionSolver(ioData.dmesh, geoSource.getMatchNodes(),domain,0);
  } else mms = 0;

  length = ioData.output.transient.length;
  surface = ioData.output.transient.surface;

  numLocSub = domain->getNumLocSub();

  dXdS=0.0;
  dFdS=0.0;
  dFdSref=0.0;
  dUdS=0.0;
  dfaU=0.0;
  dfaX=0.0;
  p=0.0;
  dPdS=0.0;
  dddx=0.0;
  dddy=0.0;
  dddz=0.0;
  dR=0.0;
  dEdgeNorm=0.0;
  dFaceNorm=0.0;
  dFaceNormVel=0.0;
  dGradP=0.0;

  FluxFD = 0.0;
  Flux = 0.0;

  reynolds0 = ioData.ref.reynolds_mu;
  kenergy0 = ioData.bc.inlet.kenergy;

  steadyTol = 1.0;
}

//------------------------------------------------------------------------------

template<int dim>
FluidShapeOptimizationHandler<dim>::~FluidShapeOptimizationHandler()
{

  if (mms) delete mms;

  if (mvp) delete mvp;

  if (dRdX) delete dRdX;

  if (pc) delete pc;

  if (ksp) delete ksp;

  if (load) delete load;

  if (dLoad) delete dLoad;

  if (dLoadref) delete dLoadref;

//  if (tsSolver) delete tsSolver;

}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoRestartBcFluxs(IoData &ioData)
{

  double gamma = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
  double R = ioData.eqs.fluidModel.gasModel.idealGasConstant;
  double Pstiff = ioData.eqs.fluidModel.gasModel.pressureConstant;

// Remark: For internal flows the SA using inlet or outlet Mach, Alpha and Beta should be specified in the input file

  ioData.bc.inlet.mach = xmach;
  ioData.bc.outlet.mach = xmach;

  ioData.bc.inlet.alpha = alprad;
  ioData.bc.outlet.alpha = alprad;

  ioData.bc.inlet.beta = teta;
  ioData.bc.outlet.beta = teta;

  if (ioData.problem.mode == ProblemData::NON_DIMENSIONAL)
  {
    this->com->fprintf(stderr, "Sensitivity Analysis does not support NON-Dimensional analysis");
    exit(-1);

    ioData.ref.mach = xmach;
    ioData.ref.dRe_mudMach = 0.0;

    if (ioData.sa.densFlag == false) {
      ioData.bc.inlet.density = 1.0;
      ioData.bc.outlet.density = 1.0;
    }

    if (ioData.sa.pressFlag == false) {
      ioData.bc.inlet.pressure = 1.0/(gamma*ioData.ref.mach*ioData.ref.mach);
      ioData.bc.outlet.pressure = 1.0/(gamma*ioData.ref.mach*ioData.ref.mach);
    }

    if (ioData.sa.apressFlag == false) {
      if (ioData.sa.pressFlag == false) {
        ioData.aero.pressure = ioData.bc.inlet.pressure;
        Pin = ioData.aero.pressure;
        this->postOp->rstVarPostFcn(ioData);
        this->postOp->rstVar(ioData);
        this->spaceOp->rstFluxFcn(ioData);
        mvp->rstSpaceOp(ioData, this->varFcn, this->spaceOp, false);
        this->rstVarImplicitCoupledTsDesc(ioData);
      }
    }
  }
  else if (ioData.problem.mode == ProblemData::DIMENSIONAL) {


    // Step 1: Re-scale all the parameters

    ioData.eqs.fluidModel.pmin *= ioData.ref.rv.pressure;

    ioData.bc.inlet.density   *= ioData.ref.rv.density;
    ioData.bc.inlet.pressure  *= ioData.ref.rv.pressure;
    ioData.bc.inlet.temperature *= ioData.ref.rv.temperature;
    ioData.bc.inlet.nutilde   *= ioData.ref.rv.nutilde;
    ioData.bc.inlet.kenergy   *= ioData.ref.rv.kenergy;
    ioData.bc.inlet.eps       *= ioData.ref.rv.epsilon;
    ioData.bc.outlet.density  *= ioData.ref.rv.density;
    ioData.bc.outlet.pressure *= ioData.ref.rv.pressure;
    ioData.bc.outlet.temperature *= ioData.ref.rv.temperature;
    ioData.bc.outlet.nutilde  *= ioData.ref.rv.nutilde;
    ioData.bc.outlet.kenergy  *= ioData.ref.rv.kenergy;
    ioData.bc.outlet.eps      *= ioData.ref.rv.epsilon;

    ioData.restart.etime  *= ioData.ref.rv.time;
    ioData.restart.dt_nm1 *= ioData.ref.rv.time;
    ioData.restart.dt_nm2 *= ioData.ref.rv.time;
    ioData.restart.energy *= ioData.ref.rv.energy;
    ioData.bc.wall.temperature *= ioData.ref.rv.temperature;
    ioData.ts.timestep    *= ioData.ref.rv.time;
    ioData.ts.maxTime     *= ioData.ref.rv.time;
    ioData.rmesh.vx       *= ioData.ref.rv.velocity;
    ioData.rmesh.vy       *= ioData.ref.rv.velocity;
    ioData.rmesh.vz       *= ioData.ref.rv.velocity;
    ioData.rmesh.ax       *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.rmesh.ay       *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.rmesh.az       *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.rmesh.timestep *= ioData.ref.rv.time;

    for (int j=0; j<ioData.rmesh.num; j++){
      ioData.rmesh.vpts[j]->time      *= ioData.ref.rv.time;
      ioData.rmesh.vpts[j]->velocityX *= ioData.ref.rv.velocity;
      ioData.rmesh.vpts[j]->velocityY *= ioData.ref.rv.velocity;
      ioData.rmesh.vpts[j]->velocityZ *= ioData.ref.rv.velocity;
    }
    ioData.aero.pressure    *= ioData.ref.rv.pressure;
    ioData.forced.timestep  *= ioData.ref.rv.time;
    ioData.forced.frequency /= ioData.ref.rv.time;

    ioData.eqs.gravity_x  *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.eqs.gravity_y  *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.eqs.gravity_z  *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.bc.hydro.depth *= ioData.ref.length;

    //
    // Step 2: Restart all values based on new inlet values
    // Step 2.1: Reset values in ioData.ref
    //

    ioData.ref.mach     = ioData.bc.inlet.mach;
    ioData.ref.density  = ioData.bc.inlet.density;
    ioData.ref.pressure = ioData.bc.inlet.pressure;
    double velocity     = ioData.ref.mach * sqrt(gamma * (ioData.ref.pressure+Pstiff) / ioData.ref.density);
    ioData.ref.temperature = (ioData.ref.pressure + gamma*Pstiff)/ (ioData.ref.density * R);
    double viscosity = ioData.eqs.viscosityModel.sutherlandConstant * sqrt(ioData.ref.temperature) /
      (1.0 + ioData.eqs.viscosityModel.sutherlandReferenceTemperature/ioData.ref.temperature);
    ioData.ref.reynolds_mu = velocity * ioData.ref.length * ioData.ref.density / viscosity;

    double dvelocitydMach = sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
    ioData.ref.dRe_mudMach = dvelocitydMach * ioData.ref.length * ioData.ref.density / viscosity;


    // Step 2.2: Reset values in ioData.ref.rv
    ioData.ref.rv.mode = RefVal::DIMENSIONAL;
    ioData.ref.rv.density = ioData.ref.density;
    ioData.ref.rv.velocity = velocity;
    ioData.ref.rv.pressure = ioData.ref.density * velocity*velocity;
    ioData.ref.rv.temperature = gamma*(gamma - 1.0) * ioData.ref.mach*ioData.ref.mach * (ioData.ref.pressure+Pstiff)/(R*ioData.ref.density);
    ioData.ref.rv.viscosity_mu = viscosity;
    ioData.ref.rv.nutilde = viscosity / ioData.ref.density;
    ioData.ref.rv.kenergy = velocity*velocity;
    ioData.ref.rv.epsilon = velocity*velocity*velocity / ioData.ref.length;
    ioData.ref.rv.time = ioData.ref.length / velocity;    // Problem in RigidMeshMotionHandler
    ioData.ref.rv.force = ioData.ref.density * velocity*velocity * ioData.ref.length*ioData.ref.length;
    ioData.ref.rv.energy = ioData.ref.density * velocity*velocity * ioData.ref.length*ioData.ref.length*ioData.ref.length;
    ioData.ref.rv.power = ioData.ref.density * velocity*velocity*velocity * ioData.ref.length*ioData.ref.length;
    ioData.ref.rv.tvelocity = velocity / ioData.aero.displacementScaling;
    ioData.ref.rv.tforce = ioData.ref.rv.force / ioData.aero.forceScaling;
    ioData.ref.rv.tpower = ioData.ref.rv.power / ioData.aero.powerScaling;

    ioData.ref.rv.dvelocitydMach = dvelocitydMach;
    ioData.ref.rv.dtimedMach = - ioData.ref.length / (velocity * velocity) * dvelocitydMach;

    ioData.eqs.fluidModel.pmin /= ioData.ref.rv.pressure;

    //
    // Step 2.3: Reset values of bc.inlet
    //

    ioData.bc.inlet.density /= ioData.ref.rv.density;
    ioData.bc.inlet.pressure /= ioData.ref.rv.pressure;
    ioData.bc.inlet.temperature /= ioData.ref.rv.temperature;
    ioData.bc.inlet.nutilde /= ioData.ref.rv.nutilde;
    ioData.bc.inlet.kenergy /= ioData.ref.rv.kenergy;
    ioData.bc.inlet.eps /= ioData.ref.rv.epsilon;

    //
    // Step 2.4: Reset values of bc.outlet
    //

    ioData.bc.outlet.density /= ioData.ref.rv.density;
    ioData.bc.outlet.pressure /= ioData.ref.rv.pressure;
    ioData.bc.outlet.temperature /= ioData.ref.rv.temperature;
    ioData.bc.outlet.nutilde /= ioData.ref.rv.nutilde;
    ioData.bc.outlet.kenergy /= ioData.ref.rv.kenergy;
    ioData.bc.outlet.eps /= ioData.ref.rv.epsilon;

    ioData.restart.etime /= ioData.ref.rv.time;
    ioData.restart.dt_nm1 /= ioData.ref.rv.time;
    ioData.restart.dt_nm2 /= ioData.ref.rv.time;
    ioData.restart.energy /= ioData.ref.rv.energy;

    ioData.bc.wall.temperature /= ioData.ref.rv.temperature;
    ioData.linearizedData.stepsize = ioData.ts.timestep;
    ioData.ts.timestep /= ioData.ref.rv.time;             // Problem in RigidRollMeshMotionHandler
    ioData.ts.maxTime /= ioData.ref.rv.time;              // Problem in RigidRollMeshMotionHandler
    ioData.rmesh.vx /= ioData.ref.rv.velocity;            // Problem in RigidMeshMotionHandler
    ioData.rmesh.vy /= ioData.ref.rv.velocity;
    ioData.rmesh.vz /= ioData.ref.rv.velocity;
    ioData.rmesh.ax /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.rmesh.ay /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.rmesh.az /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.rmesh.timestep /= ioData.ref.rv.time;          // Problem in AccMeshMotionHandler

    for (int j=0; j<ioData.rmesh.num; j++){
      ioData.rmesh.vpts[j]->time     /= ioData.ref.rv.time;
      ioData.rmesh.vpts[j]->velocityX /= ioData.ref.rv.velocity;
      ioData.rmesh.vpts[j]->velocityY /= ioData.ref.rv.velocity;
      ioData.rmesh.vpts[j]->velocityZ /= ioData.ref.rv.velocity;
    }
    ioData.aero.pressure /= ioData.ref.rv.pressure;
    ioData.forced.timestep /= ioData.ref.rv.time;         // Problem in ForcedMeshMotionHandler
    ioData.forced.frequency *= ioData.ref.rv.time;        // Problem in ForcedMeshMotionHandler

    ioData.eqs.gravity_x /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.eqs.gravity_y /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.eqs.gravity_z /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.bc.hydro.depth /= ioData.ref.length;

    double theta_k = 1.0;
    double theta_w = 10.0;
    if (kenergy0 == pow(10.0, -theta_k) * theta_w / reynolds0) {
      ioData.bc.inlet.kenergy = pow(10.0, -theta_k) * theta_w / ioData.ref.reynolds_mu;
      ioData.bc.inlet.eps = ioData.eqs.tc.tm.ke.c_mu * ioData.bc.inlet.kenergy * theta_w;
    }

    Pin = ioData.aero.pressure;

    this->postOp->rstVarPostFcn(ioData);

    this->postOp->rstVar(ioData);

    if (ioData.eqs.type == EquationsData::NAVIER_STOKES)
      this->spaceOp->rstVarFet(ioData);

    this->spaceOp->rstFluxFcn(ioData);

    mvp->rstSpaceOp(ioData, this->varFcn, this->spaceOp, false);

    this->rstVarImplicitCoupledTsDesc(ioData);

    this->bcData->rstVar(ioData);

    (this->timeState->getData()).rstVar(ioData);

    this->timeState->rstVar(ioData);

    this->output->rstVar(ioData);

    this->restart->rstVar(ioData);

    this->data->rstVar(ioData);

    this->refVal->rstVar(ioData);

    this->varFcn->rstVar(ioData);

  }

// initialize boundary condition fluxes
  this->bcData->initialize(ioData, *this->X);

}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsaRestartBcFluxs(IoData &ioData)
{
  Dev::Error(this->com,"Not yet implemented",true);
}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoGetEfforts(IoData &ioData,
                                         DistSVec<double,3> &X,
                                         DistSVec<double,dim> &U,
                                         Vec3D &F,
                                         Vec3D &M,
                                         Vec3D &L)
{

  int nSurfs = this->postOp->getNumSurf();

  Vec3D x0, force, moment;

  Vec3D *Fi = new Vec3D[nSurfs];//internal forces
  Vec3D *Mi = new Vec3D[nSurfs];//internal moments
  Vec3D *Fv = new Vec3D[nSurfs];//transmitted forces
  Vec3D *Mv = new Vec3D[nSurfs];//transmitted moments

  x0[0] = ioData.output.transient.x0;
  x0[1] = ioData.output.transient.y0;
  x0[2] = ioData.output.transient.z0;

  this->spaceOp->computeGradP(X, *this->A, U);

  this->postOp->computeForceAndMoment(x0, X, U, 0, Fi, Mi, Fv, Mv);

  F = 0.0;
  M = 0.0;

  F = Fi[0] + Fv[0];
  M = Mi[0] + Mv[0];

  if (this->refVal->mode == RefVal::NON_DIMENSIONAL) {
    F *= 2.0 * this->refVal->length*this->refVal->length / surface;
    M *= 2.0 * this->refVal->length*this->refVal->length*this->refVal->length / (surface * length);
  }
  else {
    F *= this->refVal->force;
    M *= this->refVal->energy;
  }

  L[0] = F[0]*cos(ioData.bc.inlet.alpha)*cos(ioData.bc.inlet.beta) +
         F[1]*cos(ioData.bc.inlet.alpha)*sin(ioData.bc.inlet.beta) +
         F[2]*sin(ioData.bc.inlet.alpha);

  L[1] = -F[0]*sin(ioData.bc.inlet.beta) + F[1]*cos(ioData.bc.inlet.beta);

  L[2] = -F[0]*sin(ioData.bc.inlet.alpha)*cos(ioData.bc.inlet.beta) -
          F[1]*sin(ioData.bc.inlet.alpha)*sin(ioData.bc.inlet.beta) +
          F[2]*cos(ioData.bc.inlet.alpha);

}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoGetDerivativeOfEffortsFiniteDifference(
                                         IoData &ioData,
                                         DistSVec<double,3> &X,
                                         DistSVec<double,3> &dX,
                                         DistVec<double> &A,
                                         DistSVec<double,dim> &U,
                                         DistSVec<double,dim> &dU,
                                         Vec3D &dForces,
                                         Vec3D &dMoments,
                                         Vec3D &dL)
{

  //
  //Remark: Error mesage for pointers
  //

  if (Up == 0) {
    fprintf(stderr, "*** Error: Variable Up does not exist!\n");
    exit(1);
  }
  if (Um == 0) {
    fprintf(stderr, "*** Error: Variable Um does not exist!\n");
    exit(1);
  }
  if (Xp == 0) {
    fprintf(stderr, "*** Error: Variable Xp does not exist!\n");
    exit(1);
  }
  if (Xm == 0) {
    fprintf(stderr, "*** Error: Variable Xm does not exist!\n");
    exit(1);
  }
  if (Ap == 0) {
    fprintf(stderr, "*** Error: Variable Ap does not exist!\n");
    exit(1);
  }
  if (Am == 0) {
    fprintf(stderr, "*** Error: Variable Am does not exist!\n");
    exit(1);
  }

  double eps=ioData.sa.eps;
  double xmachc;
  double alpradc;
  double tetac;
  double gamma    = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
  double velocity = ioData.ref.mach * sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
  double dVelocity= sqrt(gamma * ioData.ref.pressure / ioData.ref.density)*DFSPAR[0];
  double Force    = ioData.ref.density*ioData.ref.length*ioData.ref.length*ioData.ref.length*velocity*ioData.ref.length*velocity;
  double dForce   = 2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*velocity*dVelocity;
  double dEnergy  = 2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*ioData.ref.length*velocity*dVelocity;

  int nSurfs = this->postOp->getNumSurf();

  Vec3D x0, F, Fplus, Fminus, dF, M, Mp, Mm, dM;

  Vec3D *Fi = new Vec3D[nSurfs];
  Vec3D *Mi = new Vec3D[nSurfs];
  Vec3D *Fv = new Vec3D[nSurfs];
  Vec3D *Mv = new Vec3D[nSurfs];

  Vec3D *Fip = new Vec3D[nSurfs];
  Vec3D *Mip = new Vec3D[nSurfs];
  Vec3D *Fvp = new Vec3D[nSurfs];
  Vec3D *Mvp = new Vec3D[nSurfs];

  Vec3D *Fim = new Vec3D[nSurfs];
  Vec3D *Mim = new Vec3D[nSurfs];
  Vec3D *Fvm = new Vec3D[nSurfs];
  Vec3D *Mvm = new Vec3D[nSurfs];

  Vec3D *dFi2 = new Vec3D[nSurfs];
  Vec3D *dMi2 = new Vec3D[nSurfs];
  Vec3D *dFv2 = new Vec3D[nSurfs];
  Vec3D *dMv2 = new Vec3D[nSurfs];

  x0[0] = ioData.output.transient.x0;
  x0[1] = ioData.output.transient.y0;
  x0[2] = ioData.output.transient.z0;

  *Up = 0.0;
  *Um = 0.0;
  *Xp = 0.0;
  *Xm = 0.0;
  *Ap = 0.0;
  *Am = 0.0;

  this->spaceOp->computeGradP(X, A, U);

  this->postOp->computeForceAndMoment(x0, X, U, 0, Fi, Mi, Fv, Mv);

  F = 0.0;
  M = 0.0;

  F = Fi[0] + Fv[0];
  M = Mi[0] + Mv[0];

  xmachc=xmach;
  alpradc=alprad;
  tetac=teta;

  //
  // Evaluate efforts at U + eps * dU
  //

  Xc=X;

  *Xp=X+eps*dX;

  X=*Xp;

  *Up=U+eps*dU;

  xmach=xmachc+eps*DFSPAR[0];
  alprad=alpradc+eps*DFSPAR[1];
  teta=tetac+eps*DFSPAR[2];

  fsoRestartBcFluxs(ioData);

  this->geoState->reset(*Xp);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), *Xp, *Ap);
  this->bcData->update(*Xp);

  A=*Ap;

  this->spaceOp->computeGradP(*Xp, *Ap, *Up);

  this->postOp->computeForceAndMoment(x0, *Xp, *Up, 0, Fip, Mip, Fvp, Mvp);

  Fplus = 0.0;
  Mp = 0.0;

  Fplus = Fip[0] + Fvp[0];
  Mp = Mip[0] + Mvp[0];

  //
  // Evaluate efforts at U - eps * dU
  //

  X=Xc;

  *Xm=X-eps*dX;

  X=*Xm;

  *Um=U-eps*dU;

  xmach=xmachc-eps*DFSPAR[0];
  alprad=alpradc-eps*DFSPAR[1];
  teta=tetac-eps*DFSPAR[2];

  fsoRestartBcFluxs(ioData);

  this->geoState->reset(*Xm);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), *Xm, *Am);
  this->bcData->update(*Xm);

  A=*Am;

  this->spaceOp->computeGradP(*Xm, *Am, *Um);

  this->postOp->computeForceAndMoment(x0, *Xm, *Um, 0, Fim, Mim, Fvm, Mvm);

  Fminus = 0.0;
  Mm = 0.0;

  Fminus = Fim[0] + Fvm[0];
  Mm = Mim[0] + Mvm[0];

  //
  // Compute the derivatives
  //

  dF=1.0/(2.0*eps)*(Fplus-Fminus);
  dM=1.0/(2.0*eps)*(Mp-Mm);

  //
  // Reset values to the steady state
  //

  X=Xc;

  xmach=xmachc;
  alprad=alpradc;
  teta=tetac;

  fsoRestartBcFluxs(ioData);

  this->geoState->reset(X);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), X, A);
  this->bcData->update(X);

  this->spaceOp->computeGradP(X, A, U);

  if (this->refVal->mode == RefVal::NON_DIMENSIONAL) {
    dF *= 2.0 * this->refVal->length*this->refVal->length / surface;
    dM *= 2.0 * this->refVal->length*this->refVal->length*this->refVal->length / (surface * length);
    dForces=dF;
    dMoments=dM;
  }
  else {
    dF *= this->refVal->force;
    dM *= this->refVal->energy;
    std::cout<<"===force reference value is: "<<this->refVal->force<<std::endl;//TODO delete line
    std::cout<<"===Force value is: "<<Force<<std::endl;//TODO delete line
    std::cout<<"===Force deriv part 1: "<<dF[0]<<" "<<dF[1]<<" "<<dF[2]<<std::endl;//TODO delete line
    Vec3D t(F*dForce);
    std::cout<<"===Force deriv part 2: "<<t[0]<<" "<<t[1]<<" "<<t[2]<<std::endl;//TODO delete line
    F*=dForce;
    M*=dEnergy;
    dForces=dF+F;
    dMoments=dM+M;
  }

  double sin_a = sin(ioData.bc.inlet.alpha);
  double cos_a = cos(ioData.bc.inlet.alpha);
  double sin_b = sin(ioData.bc.inlet.beta);
  double cos_b = cos(ioData.bc.inlet.beta);

  dL[0] =  dF[0]*cos_a*cos_b + dF[1]*cos_a*sin_b + dF[2]*sin_a;
  dL[1] = -dF[0]*sin_b       + dF[1]*cos_b;
  dL[2] = -dF[0]*sin_a*cos_b - dF[1]*sin_a*sin_b + dF[2]*cos_a;

  double dsin_a = cos_a*DFSPAR[1], dcos_a = -sin_a*DFSPAR[1];
  double dsin_b = cos_b*DFSPAR[2], dcos_b = -sin_b*DFSPAR[2];

  dL[0] += F[0]*(dcos_a*cos_b + cos_a*dcos_b) +
           F[1]*(dcos_a*sin_b + cos_a*dsin_b) +
           F[2]*dsin_a;

  dL[1] += -F[0]*dsin_b + F[1]*dcos_b;

  dL[2] += -F[0]*(dsin_a*cos_b + sin_a*dcos_b) -
            F[1]*(dsin_a*sin_b + sin_a*dsin_b) +
            F[2]*dcos_a;

}


/*****************************************************************************************
 * Takes the derivative of the state vector as input and calculated the                  *
 * derivative of the desired integral quantity out of it. E.g. derivative of             *
 * Force or derivative of Lift.                                                          *
 *****************************************************************************************/
template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoGetDerivativeOfEffortsAnalytical(
                                         bool isSparse,
                                         IoData &ioData,
                                         DistSVec<double,3> &X,   //mesh motion
                                         DistSVec<double,3> &dX,  //derivative of mesh motion
                                         DistSVec<double,dim> &U, //state vector
                                         DistSVec<double,dim> &dU,//derivative of state vector
                                         Vec3D &dForces,          //derivative of forces
                                         Vec3D &dMoments,         //derivative of moments
                                         Vec3D &dL)               //derivative of lift/drag
{


  double gamma = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
  double velocity = ioData.ref.mach * sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
  double dVelocity= sqrt(gamma * ioData.ref.pressure / ioData.ref.density)*DFSPAR[0];
  double dForce =2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*velocity*dVelocity;
  double dEnergy=2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*ioData.ref.length*velocity*dVelocity;

  int nSurfs = this->postOp->getNumSurf();

  Vec3D x0, F, dF, M, dM;

  Vec3D *Fi = new Vec3D[nSurfs];
  Vec3D *Mi = new Vec3D[nSurfs];
  Vec3D *Fv = new Vec3D[nSurfs];
  Vec3D *Mv = new Vec3D[nSurfs];

  Vec3D *dFi = new Vec3D[nSurfs];
  Vec3D *dMi = new Vec3D[nSurfs];
  Vec3D *dFv = new Vec3D[nSurfs];
  Vec3D *dMv = new Vec3D[nSurfs];

  x0[0] = ioData.output.transient.x0; //x-coordinate of the point around which the moments are computed
  x0[1] = ioData.output.transient.y0; //y-coordinate of the point around which the moments are computed
  x0[2] = ioData.output.transient.z0; //z-coordinate of the point around which the moments are computed


  //it seems, that the computed grad P is only stored at the sub-domain itself and not returned
  this->spaceOp->computeGradP(X, *this->A, U);

  //computes Fi, Mi, Fv, Mv
  this->postOp->computeForceAndMoment(x0,  //reference point for moment calculation
                                  X,   //mesh motion
                    U,   //fluid state vector
                    0,   //fluid ID
                    Fi,  //internal forces
                    Mi,  //internal moments
                    Fv,  //transmitted forces
                    Mv); //transmitted moments

  F = 0.0;
  M = 0.0;

  F = Fi[0] + Fv[0];
  M = Mi[0] + Mv[0];

  dRdXoperators<dim> *dRdXop = dRdX->getdRdXop();
//  this->spaceOp->computeDerivativeOfGradP(dRdXop, dX, dAdS, dU, dddx, dddy, dddz, dR, dGradP);
//  this->postOp->computeDerivativeOfForceAndMoment(x0, X, dX, U, dU, DFSPAR, dFi, dMi, dFv, dMv);


  if(isSparse) this->spaceOp->computeDerivativeOfGradP(dRdXop, dX, dAdS, dU, dddx, dddy, dddz, dR, dGradP);
  else this->spaceOp->computeDerivativeOfGradP(X, dX, *this->A, dAdS, U, dU);

  if(isSparse) {
    this->postOp->computeDerivativeOfForceAndMoment(dRdXop, dX, dU, DFSPAR, dGradP, dFi, dMi, dFv, dMv);
  }
  else this->postOp->computeDerivativeOfForceAndMoment(x0, X, dX, U, dU, DFSPAR, dFi, dMi, dFv, dMv);


  dF = 0.0;
  dM = 0.0;

  dF = dFi[0] + dFv[0];
  dM = dMi[0] + dMv[0];

  if (this->refVal->mode == RefVal::NON_DIMENSIONAL) {
    this->com->fprintf(stderr, "Sensitivity Analysis does not support NON-Dimensional analysis");
    exit(-1);
    dF *= 2.0 * this->refVal->length*this->refVal->length / surface;
    dM *= 2.0 * this->refVal->length*this->refVal->length*this->refVal->length / (surface * length);
    dForces=dF;
    dMoments=dM;
  }
  else {
    dF *= this->refVal->force;
    dM *= this->refVal->energy;
    std::cout<<"===Force deriv part 1: "<<dF[0]<<" "<<dF[1]<<" "<<dF[2]<<std::endl;//TODO delete line
    Vec3D t(F*dForce);
    std::cout<<"===Force deriv part 2: "<<t[0]<<" "<<t[1]<<" "<<t[2]<<std::endl;//TODO delete line
    dForces = dF+F*dForce;//product rule
    dMoments = dM+M*dEnergy;//product rule
  }

//  dForces2dLifts(ioData,dForces,dL);//new version that is not finished yet
  dL = 0;
  double sin_a = sin(ioData.bc.inlet.alpha);
  double cos_a = cos(ioData.bc.inlet.alpha);
  double sin_b = sin(ioData.bc.inlet.beta);
  double cos_b = cos(ioData.bc.inlet.beta);

  dL[0] =  dF[0]*cos_a*cos_b + dF[1]*cos_a*sin_b + dF[2]*sin_a;
  dL[1] = -dF[0]*sin_b       + dF[1]*cos_b;
  dL[2] = -dF[0]*sin_a*cos_b - dF[1]*sin_a*sin_b + dF[2]*cos_a;

  double dsin_a = cos_a*DFSPAR[1], dcos_a = -sin_a*DFSPAR[1];
  double dsin_b = cos_b*DFSPAR[2], dcos_b = -sin_b*DFSPAR[2];

  double convfac = ((ioData.sa.angleRad == ioData.sa.OFF_ANGLERAD) && (DFSPAR[1] || DFSPAR[2])) ? perRad2perDeg : 1.0;

  dL[0] += (F[0]*(dcos_a*cos_b + cos_a*dcos_b) +
            F[1]*(dcos_a*sin_b + cos_a*dsin_b) +
            F[2]*dsin_a                          )*convfac;

  dL[1] += (-F[0]*dsin_b + F[1]*dcos_b           )*convfac;

  dL[2] += (-F[0]*(dsin_a*cos_b + sin_a*dcos_b) -
             F[1]*(dsin_a*sin_b + sin_a*dsin_b) +
             F[2]*dcos_a                         )*convfac;

}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoGetDerivativeOfEffortsWRTStateAndMeshPositionAnalytical(
                                         IoData &ioData,
                                         Vec3D &dForces,
                                         Vec3D &dMoments,
                                         Vec3D &dL,
                                         DistSVec<double,3> &X,
                                         DistSVec<double,dim> &U,
                                         DistSVec<double,3> &dQdX,
                                         DistSVec<double,dim> &dQdU)
{

// Q is a quantity of your interest.
// If you want Q to be first component of lift, that is L[0],
// then set dL = [1 0 0]
// If you want Q to be second component of forces,
// then set dForces = [0 1 0]

  int nSurfs = this->postOp->getNumSurf();
  if(nSurfs != 1) { this->com->fprintf(stderr, " *** Error : Sparse format supports only nSurfs = 1\n");  exit(-1); }

  Vec3D x0, F, dF, M, dM;
  Vec3D *Fi = new Vec3D[nSurfs];
  Vec3D *Mi = new Vec3D[nSurfs];
  Vec3D *Fv = new Vec3D[nSurfs];
  Vec3D *Mv = new Vec3D[nSurfs];

  x0[0] = ioData.output.transient.x0;
  x0[1] = ioData.output.transient.y0;
  x0[2] = ioData.output.transient.z0;

  this->spaceOp->computeGradP(X, *this->A, U);
  this->postOp->computeForceAndMoment(x0, X, U, 0, Fi, Mi, Fv, Mv);

  F = 0.0;
  M = 0.0;
  F = Fi[0] + Fv[0];
  M = Mi[0] + Mv[0];

  double sin_a = sin(ioData.bc.inlet.alpha);
  double cos_a = cos(ioData.bc.inlet.alpha);
  double sin_b = sin(ioData.bc.inlet.beta);
  double cos_b = cos(ioData.bc.inlet.beta);

  SVec<double,3> dFiSVec(1), dMiSVec(1), dFvSVec(1), dMvSVec(1), dSSVec(1);
  dddx = 0.0;  dddy = 0.0;  dddz = 0.0;  dR = 0.0;  dAdS = 0.0;
  // transpose computation.
  dQdX = 0.0;  dGradP = 0.0;  dQdU = 0.0;
  dF = 0.0;  dSSVec = 0.0;  dM = 0.0;

  double dLdF[3][3] = {0}, dLdS[3][3] = {0};
  dLdF[0][0] = cos_a*cos_b;  dLdF[0][1] = cos_a*sin_b;  dLdF[0][2] = sin_a;
  dLdF[1][0] = -sin_b;       dLdF[1][1] = cos_b;
  dLdF[2][0] = -sin_a*cos_b; dLdF[2][1] =-sin_a*sin_b;  dLdF[2][2] = cos_a;

  dLdS[0][1] = (F[2]*cos_a - F[0]*sin_a*cos_b - F[1]*sin_a*sin_b);
  dLdS[0][2] = (F[0]*cos_a*sin_b + F[1]*cos_a*cos_b);
  dLdS[1][2] = -(F[0]*cos_b + F[1]*sin_b);
  dLdS[2][1] = -(F[0]*cos_a*cos_b + F[1]*cos_a*sin_b + F[2]*sin_a);
  dLdS[2][2] = (F[0]*sin_a*sin_b - F[1]*sin_a*cos_b);

  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j) {
      dF[i] += dLdF[j][i]*dL[j];
      dSSVec[0][i] += dLdS[j][i]*dL[j];
    }


  dF += dForces;
  dM += dMoments;
  if(this->refVal->mode == RefVal::NON_DIMENSIONAL) {
    dF *= 2.0 * this->refVal->length*this->refVal->length / surface;
    dM *= 2.0 * this->refVal->length*this->refVal->length*this->refVal->length / (surface * length);
  }
  else {
    dF *= this->refVal->force;
    dM *= this->refVal->energy;
  }

  dFiSVec[0][0] = dF[0];    dFiSVec[0][1] = dF[1];   dFiSVec[0][2] = dF[2];
  dFvSVec[0][0] = dF[0];    dFvSVec[0][1] = dF[1];   dFvSVec[0][2] = dF[2];
  dMiSVec[0][0] = dM[0];    dMiSVec[0][1] = dM[1];   dMiSVec[0][2] = dM[2];
  dMvSVec[0][0] = dM[0];    dMvSVec[0][1] = dM[1];   dMvSVec[0][2] = dM[2];
  dRdXoperators<dim> *dRdXop = dRdX->getdRdXop();
  this->postOp->computeTransposeDerivativeOfForceAndMoment(dRdXop, dFiSVec, dMiSVec, dFvSVec, dMvSVec, dQdX, dQdU, dSSVec, dGradP);
  this->spaceOp->computeTransposeDerivativeOfGradP(dRdXop, dGradP, dddx, dddy, dddz, dR, dAdS, dQdX, dQdU);

}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoGetDerivativeOfLoadFiniteDifference(
                                         IoData &ioData,
                                         DistSVec<double,3> &X,   DistSVec<double,3> &dX,
                                         DistVec<double> &A,
                                         DistSVec<double,dim> &U,  DistSVec<double,dim> &dU,
                                         DistSVec<double,3> &load, DistSVec<double,3>   &dLoad)
{

//Remark: Error mesage for pointers
  if (Up == 0) {
    fprintf(stderr, "*** Error: Variable Up does not exist!\n");
    exit(1);
  }
  if (Um == 0) {
    fprintf(stderr, "*** Error: Variable Um does not exist!\n");
    exit(1);
  }
  if (Xp == 0) {
    fprintf(stderr, "*** Error: Variable Xp does not exist!\n");
    exit(1);
  }
  if (Xm == 0) {
    fprintf(stderr, "*** Error: Variable Xm does not exist!\n");
    exit(1);
  }
  if (Lp == 0) {
    fprintf(stderr, "*** Error: Variable Lp does not exist!\n");
    exit(1);
  }
  if (Lm == 0) {
    fprintf(stderr, "*** Error: Variable Lm does not exist!\n");
    exit(1);
  }
  if (Ap == 0) {
    fprintf(stderr, "*** Error: Variable Ap does not exist!\n");
    exit(1);
  }
  if (Am == 0) {

    fprintf(stderr, "*** Error: Variable Am does not exist!\n");
    exit(1);
  }

  double eps=ioData.sa.eps;

  double xmachc;
  double alpradc;
  double tetac;
  double gamma = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
  double velocity = ioData.ref.mach * sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
  double dVelocity= sqrt(gamma * ioData.ref.pressure / ioData.ref.density)*DFSPAR[0];
  double dForce=2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*velocity*dVelocity;

  *Up = 0.0;
  *Um = 0.0;
  *Xp = 0.0;
  *Xm = 0.0;
  *Lp = 0.0;
  *Lm = 0.0;
  *Ap = 0.0;
  *Am = 0.0;

  load=0.0;
  dLoad=0.0;

  this->spaceOp->computeGradP(X, A, U);
  this->postOp->computeNodalForce(X, U, Pin, load);

  xmachc=xmach;
  alpradc=alprad;
  tetac=teta;

  Xc=X;

  *Xp=X+eps*dX;

  X=*Xp;

  *Up=U+eps*dU;

  xmach=xmachc+eps*DFSPAR[0];
  alprad=alpradc+eps*DFSPAR[1];
  teta=tetac+eps*DFSPAR[2];

  fsoRestartBcFluxs(ioData);

  this->geoState->reset(*Xp);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), *Xp, *Ap);
  this->bcData->update(*Xp);

  A=*Ap;

  this->spaceOp->computeGradP(*Xp, *Ap, *Up);

  this->postOp->computeNodalForce(*Xp, *Up, Pin, *Lp);

  X=Xc;

  *Xm=X-eps*dX;

  X=*Xm;

  *Um=U-eps*dU;

  xmach=xmachc-eps*DFSPAR[0];
  alprad=alpradc-eps*DFSPAR[1];
  teta=tetac-eps*DFSPAR[2];

  fsoRestartBcFluxs(ioData);

  this->geoState->reset(*Xm);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), *Xm, *Am);
  this->bcData->update(*Xm);

  A=*Am;

  this->spaceOp->computeGradP(*Xm, *Am, *Um);

  this->postOp->computeNodalForce(*Xm, *Um, Pin, *Lm);

  dLoad=1.0/(2.0*eps)*((*Lp)-(*Lm));

  X=Xc;

  xmach=xmachc;
  alprad=alpradc;
  teta=tetac;

  fsoRestartBcFluxs(ioData);

  this->geoState->reset(X);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), X, A);
  this->bcData->update(X);

  this->spaceOp->computeGradP(X, A, U);

  if (this->refVal->mode == RefVal::NON_DIMENSIONAL) {
    dLoad *= 2.0 * this->refVal->length*this->refVal->length / surface;
  }
  else {
    dLoad += (dForce / this->refVal->force) * load;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoGetDerivativeOfLoadAnalytical(
                                         bool isSparse,
                                         IoData &ioData,
                                         DistSVec<double,3> &X,    DistSVec<double,3> &dX,
                                         DistSVec<double,dim> &U,  DistSVec<double,dim> &dU,
                                         DistSVec<double,3> &load, DistSVec<double,3> &dLoad)
{

  double gamma = ioData.eqs.fluidModel.gasModel.specificHeatRatio;

  double velocity = ioData.ref.mach * sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
  double dVelocity= sqrt(gamma * ioData.ref.pressure / ioData.ref.density)*DFSPAR[0];
  double dForce=2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*velocity*dVelocity;

  load=0.0;
  dLoad=0.0;

  this->spaceOp->computeGradP(X, *this->A, U);
  this->postOp->computeNodalForce(X, U, Pin, load);

  dRdXoperators<dim> *dRdXop = dRdX->getdRdXop();
  if(isSparse) {
    this->spaceOp->computeDerivativeOfGradP(dRdXop, dX, dAdS, dU, dddx, dddy, dddz, dR, dGradP);
  } else this->spaceOp->computeDerivativeOfGradP(X, dX, *this->A, dAdS, U, dU);

  //TODO: must treat dS2 better in case that dS is not zero.
  double dS2[3] = {0};

  if(isSparse) {
//    DistSVec<double,3> dLoad2(dLoad), diff(dLoad);
    this->postOp->computeDerivativeOfNodalForce(dRdXop->dForcedX, dRdXop->dForcedGradP, dRdXop->dForcedV, dRdXop->dForcedS,
                                                dRdXop->dVdU, dX, dGradP, dU, DFSPAR, dLoad);
//    this->postOp->computeDerivativeOfNodalForce(X, dX, U, dU, Pin, DFSPAR, dLoad2);
//    diff = dLoad2 - dLoad;
//    this->com->fprintf(stderr, " difference between dLoad and dLoad2 is %e and dLoad = %e, dLoad2 = %e\n", diff.norm()/dLoad.norm(), dLoad.norm(), dLoad2.norm());
  }
  else this->postOp->computeDerivativeOfNodalForce(X, dX, U, dU, Pin, DFSPAR, dLoad);

  if (this->refVal->mode == RefVal::NON_DIMENSIONAL) {
    dLoad *= 2.0 * this->refVal->length*this->refVal->length / surface;
  }
  else {
    dLoad += (dForce / this->refVal->force) * load;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoGetTransposeDerivativeOfLoadAnalytical(IoData &ioData,
                                     DistSVec<double,3> &dLoad, DistSVec<double,3> &dX, DistSVec<double,dim> &dU)
{

  double gamma = ioData.eqs.fluidModel.gasModel.specificHeatRatio;

  double velocity = ioData.ref.mach * sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
  double dVelocity= sqrt(gamma * ioData.ref.pressure / ioData.ref.density)*DFSPAR[0];
  double dForce=2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*velocity*dVelocity;


  dRdXoperators<dim> *dRdXop = dRdX->getdRdXop();

  //TODO: must treat DFSPAR2 better in case that it is not zero.
  double DFSPAR2[3] = {0};

  DistSVec<double,3> dGradP2(dGradP);
  DistVec<double> dAdS2(dAdS);

  dddx = 0.0; dddy = 0.0;  dddz = 0.0;  dR = 0.0;   dGradP2 = 0.0;    dU = 0.0;   dX = 0.0;   dAdS2 = 0.0;
  this->postOp->computeTransposeDerivativeOfNodalForce(dRdXop->dForcedX,
                                                       dRdXop->dForcedGradP,
                                                       dRdXop->dForcedV,
                                                       dRdXop->dForcedS,
                                                       dRdXop->dVdU,
                                                       dLoad, dX, dGradP2,
                                                       dU, DFSPAR2);

  this->spaceOp->computeTransposeDerivativeOfGradP(dRdXop, dGradP2, dddx, dddy, dddz, dR, dAdS2, dX, dU, false);
  dEdgeNorm = 0;  dFaceNorm = 0;  dFaceNormVel = 0;
  this->geoState->computeTransposeDerivatives(dRdXop->dEdgeNormdX, dRdXop->dFaceNormdX, dRdXop->dCtrlVoldX, dAdS2, dEdgeNorm, dFaceNorm, dFaceNormVel, dX);

  if (this->refVal->mode == RefVal::NON_DIMENSIONAL) {
    dX *= 2.0 * this->refVal->length*this->refVal->length / surface;
    dU *= 2.0 * this->refVal->length*this->refVal->length / surface;
  }
  else {
    this->com->fprintf(stderr, " Sensitivity Analysis doed not support non-dimensional calculation");  exit(-1);
    //TODO: needs to add the term below if Mach number is used as a sensitivity variable.
//  dLoad += (dForce / this->refVal->force) * load;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoSemiAnalytical
(
  IoData &ioData,
  DistSVec<double,3> &X,
  DistVec<double> &A,
  DistSVec<double,dim> &U,
  DistSVec<double,dim> &dF
)
{

  //
  // Error mesage for pointers
  //

  if (Fp == 0) {
    fprintf(stderr, "*** Error: Variable Fp does not exist!\n");
    exit(1);
  }
  if (Fm == 0) {
    fprintf(stderr, "*** Error: Variable Fm does not exist!\n");
    exit(1);
  }
  if (Xp == 0) {
    fprintf(stderr, "*** Error: Variable Xp does not exist!\n");

    exit(1);
  }
  if (Xm == 0) {
    fprintf(stderr, "*** Error: Variable Xm does not exist!\n");
    exit(1);
  }
  if (Ap == 0) {
    fprintf(stderr, "*** Error: Variable Ap does not exist!\n");
    exit(1);
  }
  if (Am == 0) {
    fprintf(stderr, "*** Error: Variable Am does not exist!\n");
    exit(1);
  }

  double dtLeft;
  double eps=ioData.sa.eps;
  double xmachc;
  double alpradc;
  double tetac;

  xmachc=xmach;
  alpradc=alprad;
  tetac=teta;

  *Fp = 0.0;
  *Fm = 0.0;
   dF = 0.0;
  *Xp = 0.0;
  *Xm = 0.0;
  *Ap = 0.0;
  *Am = 0.0;

  //
  // Compute the first flux for the FD approach
  //

  Xc=X;

  *Xp=X+eps*dXdS;

  X=*Xp;

  xmach=xmachc+eps*DFSPAR[0];
  alprad=alpradc+eps*DFSPAR[1];
  teta=tetac+eps*DFSPAR[2];

  fsoRestartBcFluxs(ioData);

  this->geoState->reset(*Xp);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), *Xp, *Ap);
  this->bcData->update(*Xp);

  A=*Ap;

  dtLeft = 0.0;
  this->computeTimeStep(1, &dtLeft, U);

  this->spaceOp->computeResidual(*Xp, *Ap, U, *Fp, this->timeState);

  this->spaceOp->applyBCsToResidual(U, *Fp);

  //
  // Compute the second flux for the FD approach
  //

  X=Xc;

  *Xm=X-eps*dXdS;

  X=*Xm;

  xmach=xmachc-eps*DFSPAR[0];
  alprad=alpradc-eps*DFSPAR[1];
  teta=tetac-eps*DFSPAR[2];

  fsoRestartBcFluxs(ioData);

  this->geoState->reset(*Xm);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), *Xm, *Am);
  this->bcData->update(*Xm);

  A=*Am;

  dtLeft = 0.0;
  this->computeTimeStep(1, &dtLeft, U);

  this->spaceOp->computeResidual(*Xm, *Am, U, *Fm, this->timeState);

  this->spaceOp->applyBCsToResidual(U, *Fm);

  dAdS=1.0/(2.0*eps)*((*Ap)-(*Am));

  dF=1.0/(2.0*eps)*((*Fp)-(*Fm));
  if((ioData.sa.angleRad == ioData.sa.OFF_ANGLERAD) && (DFSPAR[1] || DFSPAR[2]))
     dF *= perRad2perDeg;
  //
  // Reset the steady state
  //

  X=Xc;

  xmach=xmachc;
  alprad=alpradc;
  teta=tetac;

  fsoRestartBcFluxs(ioData);

  this->geoState->reset(X);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), X, A);
  this->bcData->update(X);

  dtLeft = 0.0;
  this->computeTimeStep(1, &dtLeft, U);

  this->spaceOp->computeGradP(X, A, U);

  if ( ioData.sa.scFlag == SensitivityAnalysis::SEMIANALYTICAL ) {
    this->geoState->updateConfigSA();
    this->bcData->initializeSA(ioData, X, dXdS, DFSPAR[0], DFSPAR[1], DFSPAR[2]);
  }

}


//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoAnalytical
(bool isSparse, IoData &ioData, DistSVec<double,3> &X, DistSVec<double,3> &dXdS, DistVec<double> &A, DistSVec<double,dim> &U, DistSVec<double,dim> &dFdS)
{
  //
  // Computing the normal, derivative of the normal and of the control volume
  //
  DistSVec<double,dim> dddx2(dddx), dddy2(dddy), dddz2(dddz);
  DistSVec<double,6> dR2(dR);
  DistVec<double> dFaceNormVel(domain->getFaceNormDistInfo());
  dEdgeNorm = 0.0;  dFaceNorm = 0.0;  dddx = 0.0;  dddy = 0.0;  dddz = 0.0;  dFaceNormVel = 0.0;  dR = 0.0;


  dRdXoperators<dim> *dRdXop = dRdX->getdRdXop();
  isSparse=false;//TODO HACK
  if(isSparse) {

    this->geoState->computeDerivatives(dRdXop->dEdgeNormdX,
                                       dRdXop->dFaceNormdX,
                                       dRdXop->dCtrlVoldX,
                                       X,
                                       dXdS,
                                       dAdS,
                                       dEdgeNorm,
                                       dFaceNorm,
                                       dFaceNormVel);
  }
  else{//non-sparse routine
    std::cout<<"geoState->computeDerivatives X: "<<X.norm()<<"  dXdS: "<<dXdS.norm()<<"  dAdS: "<<dAdS.norm()<<std::endl;//TODO delete line
    std::cout<<"\033[94mFULL this->bcData->getVelocityVector() \033[00m"<<this->bcData->getVelocityVector().norm()<<std::endl;//TODO delete line
    std::cout<<"\033[94mFULL this->bcData->getDerivativeOfVelocityVector() \033[00m"<<this->bcData->getDerivativeOfVelocityVector().norm()<<std::endl;//TODO delete line
    this->geoState->computeDerivatives(X,
                                       dXdS,
                                       this->bcData->getVelocityVector(),            //velocity of the boundary nodes
                                       this->bcData->getDerivativeOfVelocityVector(),//derivative of the boundary nodes velocity
                                       dAdS);
  }
  //
  // Computing the derivatives of the boundary fluxes
  //
  this->bcData->initializeSA(ioData, X, dXdS, DFSPAR[0], DFSPAR[1], DFSPAR[2]);

  //
  // Computing the partial derivative of the flux with respect to the variables
  //
  if(isSparse) {
  //      DistSVec<double,dim> dFdS2(dFdS), diff(dFdS);
  //      this->spaceOp->computeDerivativeOfResidual(X, dXdS, A, dAdS, U, DFSPAR[0], Flux, dFdS2, this->timeState);
	  this->spaceOp->computeDerivativeOfResidual(dRdXop, X, dXdS, A, dAdS, dEdgeNorm, dFaceNorm, dFaceNormVel, Flux, dFdS, dR, dddx, dddy, dddz, U, DFSPAR[0], this->timeState);
  //      diff = dFdS2 - dFdS;
  //      this->com->fprintf(stderr, " diff between sparse_dFdS and sparse_dFdS is %e, dFdS is %e, nonsparse_dFdS2 is %e\n", diff.norm()/dFdS.norm(), dFdS.norm(), dFdS2.norm());
  }
  else//non-sparse version
  {
      this->spaceOp->computeDerivativeOfResidual(X, dXdS, A, dAdS, U, DFSPAR[0], Flux, dFdS, this->timeState);
  }

  this->spaceOp->applyBCsToDerivativeOfResidual(U, dFdS);
  //if(dFdS.norm()<1e-10){std::cout<<"ERRROR"<<std::endl; exit(-1);}//TODO delete line
  if((ioData.sa.angleRad == ioData.sa.OFF_ANGLERAD)&& (DFSPAR[1] || DFSPAR[2]))
    dFdS *= Deg2Rad;

}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoApply_dFdXtranspose
(DistVec<double> &A, DistSVec<double,dim> &lambdaU, DistSVec<double,3> &rhs)
{

  DistVec<double> dAdS2(dAdS);
  DistSVec<double,dim> lambdaU1(lambdaU);

  dEdgeNorm = 0.0;  dFaceNorm = 0.0;  dddx = 0.0;  dddy = 0.0;  dddz = 0.0;  dFaceNormVel = 0.0;  dR = 0.0;   rhs = 0.0;  dAdS2 = 0.0;  dFaceNormVel = 0.0;
  dRdXoperators<dim> *dRdXop = dRdX->getdRdXop();

  this->spaceOp->computeTransposeDerivativeOfResidual(dRdXop, Flux, lambdaU1, A, dAdS2, rhs, dddx, dddy, dddz, dEdgeNorm, dFaceNorm, dFaceNormVel, dR);
  this->geoState->computeTransposeDerivatives(dRdXop->dEdgeNormdX, dRdXop->dFaceNormdX, dRdXop->dCtrlVoldX, dAdS2, dEdgeNorm, dFaceNorm, dFaceNormVel, rhs);

}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoSetUpLinearSolver(IoData &ioData, DistSVec<double,3> &X, DistVec<double> &A,
                                                              DistSVec<double,dim> &U, DistSVec<double,dim> &dFdS)
{

// Preparing the linear solver
  fsoRestartBcFluxs(ioData);

  this->geoState->reset(X);

  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), X, A);

  this->bcData->update(X);

  this->spaceOp->computeResidual(X, A, U, FluxFD, this->timeState);

  if (ioData.sa.homotopy == SensitivityAnalysis::ON_HOMOTOPY)
    this->timeState->add_dAW_dt(1, *this->geoState, A, U, FluxFD);

  this->spaceOp->applyBCsToResidual(U, FluxFD);

  mvp->evaluate(0, X, A, U, FluxFD);

  DistMat<PrecScalar,dim> *_pc = dynamic_cast<DistMat<PrecScalar,dim> *>(pc);

  if (_pc) {
    MatVecProdFD<dim,dim> *mvpfd = dynamic_cast<MatVecProdFD<dim,dim> *>(mvp);
    MatVecProdH2<dim,MatScalar,dim> *mvph2 = dynamic_cast<MatVecProdH2<dim,MatScalar,dim> *>(mvp);

    if (mvpfd || mvph2)
    {
      this->spaceOp->computeJacobian(X, A, U, *_pc, this->timeState);
      if (ioData.sa.homotopy == SensitivityAnalysis::ON_HOMOTOPY)
        this->timeState->addToJacobian(A, *_pc, U);
      this->spaceOp->applyBCsToJacobian(U, *_pc);
    }

  } // END if (_pc)

  pc->setup();

  // Computing flux for compatibility correction of the derivative of the flux
  this->spaceOp->computeResidual(X, A, U, Flux, this->timeState, false);

}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoSetUpAdjointLinearSolver(
                                         IoData &ioData, DistSVec<double,3> &X,
                                         DistVec<double> &A,
                                         DistSVec<double,dim> &U,
                                         DistSVec<double,dim> &dFdS)
{

// Preparing the linear solver

  fsoRestartBcFluxs(ioData);

  this->geoState->reset(X);

  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), X, A);

  this->bcData->update(X);

  this->spaceOp->computeResidual(X, A, U, FluxFD, this->timeState);

  if (ioData.sa.homotopy == SensitivityAnalysis::ON_HOMOTOPY)
    this->timeState->add_dAW_dt(1, *this->geoState, A, U, FluxFD);

  this->spaceOp->applyBCsToResidual(U, FluxFD);

  mvp->evaluate(0, X, A, U, FluxFD);

  DistMat<PrecScalar,dim> *_pc = dynamic_cast<DistMat<PrecScalar,dim> *>(pc);

  if (_pc) {
    MatVecProdFD<dim,dim> *mvpfd = dynamic_cast<MatVecProdFD<dim,dim> *>(mvp);
    MatVecProdH2<dim,MatScalar,dim> *mvph2 = dynamic_cast<MatVecProdH2<dim,MatScalar,dim> *>(mvp);

    if (mvpfd || mvph2){
      this->spaceOp->computeJacobian(X, A, U, *_pc, this->timeState);
      if (ioData.sa.homotopy == SensitivityAnalysis::ON_HOMOTOPY)
        this->timeState->addToJacobian(A, *_pc, U);
      this->spaceOp->applyBCsToJacobian(U, *_pc);
    }

  } // END if (_pc)

  pc->setup();
  pc->setupTR();

  // Computing flux for compatibility correction of the derivative of the flux
  this->spaceOp->computeResidual(X, A, U, Flux, this->timeState, false);
}

//------------------------------------------------------------------------------


/******************************************************************************
 * Solves the equation [dFdU]*[dUdS]=[dFdS] for dUdS by an iterative solver   *
 * The matrix dFdU, that is the Jacobian is impicitly stored within the       *
 * KspSolver object                                                           *
 ******************************************************************************/
template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoLinearSolver(
                                         IoData &ioData,
                                         DistSVec<double,dim> &dFdS,
                                         DistSVec<double,dim> &dUdS,
                                         bool isFSI
)
{
  fsoPrintTextOnScreen("Starting LinearSolver");

//  dUdS = 0.0;

  dFdS *= (-1.0);
  if(!isFSI) ksp->setup(0, 1, dFdS);

  int numberIteration;
  bool istop = false;
  int iter = 0;

  while ((istop == false) && (iter < 100))
  {
    numberIteration = ksp->solve(dFdS, dUdS);
    if ((!ioData.sa.excsol) || (numberIteration < ioData.sa.ksp.maxIts))
      istop = true;
    iter += 1;

  }

  dFdS *= (-1.0);

  fsoPrintTextOnScreen("Finished LinearSolver");

}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoAdjointLinearSolver
(
  IoData &ioData,
  DistSVec<double,dim> &dQdU, DistSVec<double,dim> &lambdaU,
  bool isFSI
)
{

  DistSVec<double,dim> rhs(dQdU);
  if(isFSI) {
    dfaX = 0;  dfaU = 0;
    fsoGetTransposeDerivativeOfLoadAnalytical(ioData, lambdaSDisp, dfaX, dfaU);
    rhs += dfaU;
  }

  if(!isFSI) ksp->setup(0, 1, rhs);

  int numberIteration;
  bool istop = false;
  int iter = 0;

  while ((istop == false) && (iter < 100))
  {
    numberIteration = ksp->solveT(rhs, lambdaU);
    //ioData.sa.excsol is an important parameter. If it is set to 0, only one iteration will be carried out
    if ((!ioData.sa.excsol) || (numberIteration < ioData.sa.ksp.maxIts))
      istop = true;
    iter += 1;
  }

}

//------------------------------------------------------------------------------


template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoPrintTextOnScreen(const char *Text)
{
   this->com->fprintf(stderr, Text);
}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsaPrintTextOnScreen(const char *Text)
{
   this->com->fprintf(stderr, Text);
}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoOutput1D(const char *fileName, DistVec<double> &V)
{

outFile = fopen(fileName,"w");


#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    for (int j=0; j<V(iSub).size(); ++j)
        fprintf(outFile," %20.17e \n",V(iSub)[j]);


fclose(outFile);

}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoOutput3D(const char *fileName, DistSVec<double,3> &V)
{

outFile = fopen(fileName,"w");

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    for (int j=0; j<V(iSub).size(); ++j) {
      for (int k=0; k<3; ++k)
        fprintf(outFile," %20.17e ",V(iSub)[j][k]);
      fprintf(outFile,"\n");

    }
  }

fclose(outFile);

}

//------------------------------------------------------------------------------


template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoOutputDimD(const char *fileName, DistSVec<double,dim> &V)
{


outFile = fopen(fileName,"w");

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub ; ++iSub) {
    for (int j=0; j<V(iSub).size(); ++j) {
      for (int k=0; k<dim; ++k)
        fprintf(outFile," %20.17e ",V(iSub)[j][k]);
      fprintf(outFile,"\n");
    }

  }

fclose(outFile);

}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoInitialize(IoData &ioData, DistSVec<double,dim> &U)
{
  this->output->openAsciiFiles();

  // Reseting the configuration control of the geometry datas
  this->geoState->resetConfigSA();

  if (this->com->cpuNum() == 0) {
    outFile = fopen(ioData.sa.sensoutput,"w");
    if (outFile)
      fclose(outFile);
  }

  xmach = ioData.sa.machref;
  alprad = ioData.sa.alpharef;
  teta = ioData.sa.betaref;

  double dtLeft = 0.0;
  this->computeTimeStep(1, &dtLeft, U);

  this->computeMeshMetrics();
  this->updateStateVectors(U);

  // Setting up the linear solver.
  fsoSetUpLinearSolver(ioData, *this->X, *this->A, U, dFdS);
}

//------------------------------------------------------------------------------

template<int dim>
int FluidShapeOptimizationHandler<dim>::fsoHandler(IoData &ioData, DistSVec<double,dim> &U)
{
  // xmach      -  Mach number
  // alpha      -  pitch angle
  // teta       -  yaw angle
  // DFSPAR[1]  -  This is 1.0 if the sensitivity with respect ot the mach number should be caluclated
  // DFSPAR[2]  -  This is 1.0 if the sensitivity with respect ot the AoA (alpha) should be caluclated
  // DFSPAR[3]  -  This is 1.0 if the sensitivity with respect ot the yq-angle (beta) should be caluclated

  // Start basic timer
  double MyLocalTimer = -this->timer->getTime();

  double dtLeft = 0.0;
  this->computeTimeStep(1, &dtLeft, U);
  this->computeMeshMetrics();
  this->updateStateVectors(U);


  bool isSparse       = bool(ioData.sa.sparseFlag);
  if(ioData.sa.method == SensitivityAnalysis::ADJOINT) isSparse = true;

  if(isSparse) {
    Vec3D x0;
    x0[0] = ioData.output.transient.x0;
    x0[1] = ioData.output.transient.y0;
    x0[2] = ioData.output.transient.z0;
    dRdX->constructOperators(x0, *this->X, *this->A, U, DFSPAR[0], Flux, Pin, this->timeState, this->postOp);
  }

  //Adjoint routines
  if (ioData.sa.method    == SensitivityAnalysis::ADJOINT &&
       (ioData.sa.sensMach  == SensitivityAnalysis::ON_SENSITIVITYMACH  ||
        ioData.sa.sensAlpha == SensitivityAnalysis::ON_SENSITIVITYALPHA ||
        ioData.sa.sensBeta  == SensitivityAnalysis::ON_SENSITIVITYBETA)    ){

      this->com->fprintf(stderr, "\033[91m\nERROR: Adjoint SA currently only supports mesh-sensitivity\n\n\033[00m");
      exit(-1);
  }
  //Direct routines
  else if (ioData.sa.method == SensitivityAnalysis::DIRECT) {
    fsoSetUpLinearSolver(ioData, *this->X, *this->A, U, dFdS);
    if (ioData.sa.sensMesh  == SensitivityAnalysis::ON_SENSITIVITYMESH)  fso_on_sensitivityMesh(isSparse, ioData, U);
    if (ioData.sa.sensMach  == SensitivityAnalysis::ON_SENSITIVITYMACH)  fso_on_sensitivityMach(isSparse, ioData, U);
    if (ioData.sa.sensAlpha == SensitivityAnalysis::ON_SENSITIVITYALPHA) fso_on_sensitivityAlpha(isSparse, ioData, U);
    if (ioData.sa.sensBeta  == SensitivityAnalysis::ON_SENSITIVITYBETA)  fso_on_sensitivityBeta(isSparse, ioData, U);
  }
  else if(ioData.sa.method  == SensitivityAnalysis::ADJOINT) {
    fsoSetUpAdjointLinearSolver(ioData, *this->X, *this->A, U, dFdS);
    if (ioData.sa.sensMesh  == SensitivityAnalysis::ON_SENSITIVITYMESH)
      fso_on_AdjointSensitivityMesh(ioData, U);
  }

  bool lastIt = true;

  this->output->closeAsciiFiles();


  MyLocalTimer += this->timer->getTime();
  if (this->com->cpuNum() == 0)
  {
    std::cout << "\n *** FluidShapeOptimizationHandler::fsoHandler >> Exit";
    std::cout << " (" << MyLocalTimer << " s)";
    std::cout << "\n\n";
  }

  return -1;

}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::setDFSPAR(IoData &ioData)
{
  switch (actvar) {
    case 1:
      if(ioData.sa.sensMesh == SensitivityAnalysis::OFF_SENSITIVITYMESH) {
        fsoPrintTextOnScreen(" ***** Error: SensitivityMesh in fluid input must be on\n"); exit(-1); }
      DFSPAR[0] = 0;  DFSPAR[1] = 0;  DFSPAR[2] = 0;
      break;
    case 2:
      if(ioData.sa.sensMach == SensitivityAnalysis::OFF_SENSITIVITYMACH) {
        fsoPrintTextOnScreen(" ***** Error: SensitivityMach in fluid input must be on\n"); exit(-1); }
      DFSPAR[0] = 1;  DFSPAR[1] = 0;  DFSPAR[2] = 0;
      break;
    case 3:
      if(ioData.sa.sensAlpha == SensitivityAnalysis::OFF_SENSITIVITYALPHA) {
        fsoPrintTextOnScreen(" ***** Error: SensitivityAlpha in fluid input must be on\n"); exit(-1); }
      DFSPAR[0] = 0;  DFSPAR[1] = 1;  DFSPAR[2] = 0;
      break;
    case 4:
      if(ioData.sa.sensBeta == SensitivityAnalysis::OFF_SENSITIVITYBETA) {
        fsoPrintTextOnScreen(" ***** Error: SensitivityBeta in fluid input must be on\n"); exit(-1); }
      DFSPAR[0] = 0;  DFSPAR[1] = 0;  DFSPAR[2] = 1;
      break;
    case 5:
      if(ioData.sa.sensFSI == SensitivityAnalysis::OFF_SENSITIVITYFSI) {
        fsoPrintTextOnScreen(" ***** Error: SensitivityFSI in fluid input must be on\n"); exit(-1); }
      DFSPAR[0] = 0;  DFSPAR[1] = 0;  DFSPAR[2] = 0;
      break;
    default:
      fsoPrintTextOnScreen(" ***** Error: invalid value for active sensitivity variable!\n"); exit(-1);
      break;
  }

  //New implementations that will replace DFSPAR
  switch (actvar) {
    case 1:
      if(ioData.sa.sensMesh == SensitivityAnalysis::OFF_SENSITIVITYMESH) { fsoPrintTextOnScreen(" ***** Error: SensitivityMesh in fluid input must be on\n"); exit(-1); }
      senstype.mach = false;  senstype.alpha = false;  senstype.beta = false; senstype.mesh=true;
      break;
    case 2:
      if(ioData.sa.sensMach == SensitivityAnalysis::OFF_SENSITIVITYMACH) { fsoPrintTextOnScreen(" ***** Error: SensitivityMach in fluid input must be on\n"); exit(-1); }
      senstype.mach = true;  senstype.alpha = false;  senstype.beta = false; senstype.mesh=false;
      break;
    case 3:
      if(ioData.sa.sensAlpha == SensitivityAnalysis::OFF_SENSITIVITYALPHA) { fsoPrintTextOnScreen(" ***** Error: SensitivityAlpha in fluid input must be on\n"); exit(-1); }
      senstype.mach = false;  senstype.alpha = true;  senstype.beta = false; senstype.mesh=false;
      break;
    case 4:
      if(ioData.sa.sensBeta == SensitivityAnalysis::OFF_SENSITIVITYBETA) { fsoPrintTextOnScreen(" ***** Error: SensitivityBeta in fluid input must be on\n"); exit(-1); }
      senstype.mach = false;  senstype.alpha = false;  senstype.beta = true; senstype.mesh=false;
      break;
    case 5:
      if(ioData.sa.sensFSI == SensitivityAnalysis::OFF_SENSITIVITYFSI) { fsoPrintTextOnScreen(" ***** Error: SensitivityFSI in fluid input must be on\n"); exit(-1); }
      senstype.mach = false;  senstype.alpha = false;  senstype.beta = false; senstype.mesh=false;
      break;
    default:
      fsoPrintTextOnScreen(" ***** Error: invalid value for active sensitivity variable!\n"); exit(-1);
      break;
  }
}

//------------------------------------------------------------------------------

template<int dim>
int FluidShapeOptimizationHandler<dim>::fsoAeroelasticHandler(IoData &ioData, DistSVec<double,dim> &U)
{

  // xmach      -  Mach number
  // alpha      -  pitch angle
  // teta       -  yaw angle
  // DFSPAR[1]  -  Mach number differential
  // DFSPAR[2]  -  angle of attack differential
  // DFSPAR[3]  -  yaw angle differential

  // Start basic timer
  double MyLocalTimer = -this->timer->getTime();

  bool isSparse = bool(ioData.sa.sparseFlag);
  if(ioData.sa.method == SensitivityAnalysis::ADJOINT) isSparse = true;
  double dtLeft = 0.0;
  this->computeTimeStep(1, &dtLeft, U);
  this->computeMeshMetrics();
  this->updateStateVectors(U);


  if(ioData.sa.method == SensitivityAnalysis::DIRECT) {
    fsoSetUpLinearSolver(ioData, *this->X, *this->A, U, dFdS);
    int totalNumParamTypes;
    this->getNumParam(totalNumParamTypes,actvar,steadyTol);
    if(ioData.sa.sensMach  == SensitivityAnalysis::ON_SENSITIVITYMACH) { totalNumParamTypes++; }
    if(ioData.sa.sensAlpha == SensitivityAnalysis::ON_SENSITIVITYALPHA){ totalNumParamTypes++; }
    if(ioData.sa.sensBeta  == SensitivityAnalysis::ON_SENSITIVITYBETA) { totalNumParamTypes++; }

    if(isSparse) {
      Vec3D x0;
      x0[0] = ioData.output.transient.x0;
      x0[1] = ioData.output.transient.y0;
      x0[2] = ioData.output.transient.z0;
      dRdX->constructOperators(x0, *this->X, *this->A, U, DFSPAR[0], Flux, Pin, this->timeState, this->postOp);
    }

    for(int iparam=0; iparam<totalNumParamTypes; ++iparam) {
      int numParam;
      this->getNumParam(numParam,actvar,steadyTol);
      setDFSPAR(ioData);
      for(int i=0; i<numParam; ++i) fso_on_aeroelasticSensitivityFSI(isSparse, ioData, U);
    }
    bool lastIt = true;

    this->output->closeAsciiFiles();

  } else if(ioData.sa.method == SensitivityAnalysis::ADJOINT) {
    fsoSetUpAdjointLinearSolver(ioData, *this->X, *this->A, U, dFdS);
    int totalNumStructureQuantities(0), totalNumFluidQuantities(0);
    this->getNumParam(totalNumStructureQuantities,actvar,steadyTol);
    if(ioData.output.transient.dLiftx[0] != 0) { totalNumFluidQuantities++; }
    if(ioData.output.transient.dLifty[0] != 0) { totalNumFluidQuantities++; }
    if(ioData.output.transient.dLiftz[0] != 0) { totalNumFluidQuantities++; }
    this->sendNumParam(totalNumFluidQuantities);
    this->com->fprintf(stderr, " ... In fluid, totalNumStructureQuantities = %d, totalNumFluidQuantities = %d\n",
                                totalNumStructureQuantities, totalNumFluidQuantities);
    if(isSparse) {

      Vec3D x0;
      x0[0] = ioData.output.transient.x0;
      x0[1] = ioData.output.transient.y0;
      x0[2] = ioData.output.transient.z0;

      dRdX->constructOperators(x0, *this->X, *this->A, U, DFSPAR[0], Flux, Pin, this->timeState, this->postOp);
    }

    // structure sensitivity quantities (e.g., aggregated von Mises stress, tip displacement of structure, etc) are computed first
    for(int iquan=0; iquan<totalNumStructureQuantities; ++iquan) {
      int numQuans;
      this->getNumParam(numQuans,actvar,steadyTol);  // numQuans is number of quantities related to a sensitivity quantity type
                                                     // actvar = 1: shape variable
                                                     // actvar = 2: Mach number variable
                                                     // actvar = 3: angle of attack variable
                                                     // actvar = 4: yaw angle variable
                                                     // actvar = 5: FSI variable such as structure thickness
                                                     // actvar = 6: aggregated von Mises stress quantity
                                                     // actvar = 7: von Mises stress quantity
                                                     // actvar = 8: structure displacement quantity
                                                     // actvar = 9: liftx quantity
                                                     // actvar = 10: lifty quantity
                                                     // actvar = 11: liftz quantity
      this->com->fprintf(stderr, " ... In fluid, numQuans = %d, actvar = %d, steadyTol = %e\n", numQuans, actvar, steadyTol);
      for(int i=0; i<numQuans; ++i) fso_on_aeroelasticAdjointSensitivityFSI(ioData, U);
    }
    int numStructParamTypes(0);
    this->getNumParam(numStructParamTypes, actvar, steadyTol);
    this->com->fprintf(stderr, " ... numStructParamTypes = %d\n", numStructParamTypes);
    if(ioData.output.transient.dLiftx[0] != 0) {
      this->com->fprintf(stderr, "x-direction lift sensitivity will be computed\n"); actvar = 9;
      fso_on_aeroelasticAdjointSensitivityFSI(ioData, U);
      for(int iStParam=0; iStParam<numStructParamTypes; ++iStParam) {
        double dlift(0);
        this->getRelResidual(dlift);
        this->output->writeDerivativeOfLiftxToDisk(dlift);
      }
    }
    if(ioData.output.transient.dLifty[0] != 0) {
      this->com->fprintf(stderr, "y-direction lift sensitivity will be computed\n"); actvar = 10;
      fso_on_aeroelasticAdjointSensitivityFSI(ioData, U);
      for(int iStParam=0; iStParam<numStructParamTypes; ++iStParam) {
        double dlift(0);
        this->getRelResidual(dlift);
        this->output->writeDerivativeOfLiftyToDisk(dlift);
      }
    }
    if(ioData.output.transient.dLiftz[0] != 0) {
      this->com->fprintf(stderr, "z-direction lift sensitivity will be computed\n"); actvar = 11;
      fso_on_aeroelasticAdjointSensitivityFSI(ioData, U);
      for(int iStParam=0; iStParam<numStructParamTypes; ++iStParam) {
        double dlift(0);
        this->getRelResidual(dlift);
        this->output->writeDerivativeOfLiftzToDisk(dlift);
      }
    }



    bool lastIt = true;

    this->output->closeAsciiFiles();

  }


  MyLocalTimer += this->timer->getTime();
  if (this->com->cpuNum() == 0)
  {
    std::cout << "\n *** FluidShapeOptimizationHandler::fsoAeroelasticHandler >> Exit";
    std::cout << " (" << MyLocalTimer << " s)";
    std::cout << "\n\n";
  }

  return -1;

}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fso_on_sensitivityBeta(bool isSparse, IoData &ioData, DistSVec<double,dim> &U)
{
    dXdS = 0.0;
    dAdS = 0.0;
    DFSPAR[0] = 0.0;
    DFSPAR[1] = 0.0;
    DFSPAR[2] = 1.0;
    actvar = 4;

    if (!ioData.sa.angleRad)
      ioData.sa.eps *= perRad2perDeg;

    fsoComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U, false, isSparse);

    fsoComputeSensitivities(isSparse, ioData, "Derivatives with respect to the yaw angle:", ioData.sa.sensoutput, *this->X, U);

    fsoPrintTextOnScreen("\n ***** Derivatives with respect to the yaw angle were computed! \n");

    step = step + 1;

    if (!ioData.sa.angleRad)
      ioData.sa.eps *= perDeg2perRad;
}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fso_on_sensitivityAlpha(bool isSparse, IoData &ioData, DistSVec<double,dim> &U)
{
    dXdS = 0.0;
    dAdS = 0.0;
    DFSPAR[0] = 0.0;
    DFSPAR[1] = 1.0;
    DFSPAR[2] = 0.0;
    actvar = 3;

    if (!ioData.sa.angleRad)
      ioData.sa.eps *= perRad2perDeg;

    fsoComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U, false, isSparse);

    fsoComputeSensitivities(isSparse, ioData, "Derivatives with respect to the angle of attack:", ioData.sa.sensoutput, *this->X, U);

    fsoPrintTextOnScreen("\n ***** Derivatives with respect to the angle of attack were computed! \n");

    step = step + 1;

    if (!ioData.sa.angleRad)
      ioData.sa.eps *= perDeg2perRad;
}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fso_on_sensitivityMach(bool isSparse, IoData &ioData, DistSVec<double,dim> &U)
{
    dXdS = 0.0;
    dAdS = 0.0;
    DFSPAR[0] = 1.0;
    DFSPAR[1] = 0.0;
    DFSPAR[2] = 0.0;
    actvar = 2;

    fsoComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U, false, isSparse);

    fsoComputeSensitivities(isSparse, ioData, "Derivatives with respect to the Mach number:", ioData.sa.sensoutput, *this->X, U);

    fsoPrintTextOnScreen("\n ***** Derivatives with respect to the Mach number were computed! \n");

    step = step + 1;
}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fso_on_aeroelasticSensitivityFSI(bool isSparse, IoData &ioData, DistSVec<double,dim> &U)
{

    double tag = 0.0;
    bool lastIt = false;
    int iter = 0;
    dXdS = 0.0;
    dXdSb = 0.0;
    dAdS = 0.0;

    while (!lastIt) {

      this->cmdCom(&lastIt);
      if(lastIt) { dXdSb = 0.0; break; }
      this->com->fprintf(stderr, "fso_aeroelatic_sensitivity Iteration\t");
      if(ioData.sa.adaptiveEpsFSI) {
        double relres;
        this->getRelResidual(relres);
        ksp->setEps(relres);
      } else ksp->setEps(steadyTol);
      // Reading derivative of the overall deformation
      this->receiveBoundaryPositionSensitivityVector(dXdSb); // [F] receive boundary displacement sensitivity from structure ...

      // Checking if dXdSb has entries different from zero at the interior of the mesh
      this->postOp->checkVec(dXdSb);

      if (dXdSb.norm() == 0.0)
      {
        this->com->fprintf(stderr, "\n *** WARNING *** No Surface Mesh Sensitivity Perturbation \n\n");
        if(!ioData.sa.fsiFlag) exit(1);
      }

      // Updating the mesh
      dXdS = *this->X;
      mms->solve(dXdSb, dXdS);
      dXdS -= *this->X;

      // Check that the mesh perturbation is propagated
      if (dXdS.norm() == 0.0) this->com->fprintf(stderr, "\n !!! WARNING !!! No Mesh Sensitivity Perturbation !!!\n\n");

      fsoComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U, true, isSparse);

      fsoComputeAndSendForceSensitivities(isSparse, ioData, ioData.sa.sensoutput, *this->X, U);

      dXdSb = 0.0;
      iter++;
    }

    fsoComputeSensitivities(isSparse, ioData, "Derivatives with respect to the FSI parameter:", ioData.sa.sensoutput, *this->X, U);
    fsoPrintTextOnScreen("\n ***** Derivatives of mesh position and state were computed! \n");
    step++;

}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fso_on_aeroelasticAdjointSensitivityFSI(IoData &ioData, DistSVec<double,dim> &U)
{

  double tag = 0.0;
  bool lastIt = false;
  int iter = 0;
  dXdS = 0.0;
  dAdS = 0.0;
  lambdaSDisp = 0.0;
  DistSVec<double,3> dQdX(*this->X);
  DistSVec<double,dim> dQdU(U);
  Vec3D dForces(0.0), dMoments(0.0), dL(0.0);
  double relres, relres_p(1.0), relres_pp(1.0);

  while (!lastIt) {

    this->cmdCom(&lastIt);
    if(lastIt) { dXdSb = 0.0; break; }
    this->com->fprintf(stderr, "fso_aeroelatic_adjoint_sensitivity Iteration\t");
    this->getRelResidual(relres);
    if(ioData.sa.adaptiveEpsFSI && relres < relres_p && relres < relres_pp && relres < 1.0) {
      ksp->setEps(relres);
      relres_pp = relres_p;
      relres_p = relres;
    } else ksp->setEps(steadyTol);
    // Reading derivative of the overall deformation
    this->receiveBoundaryPositionSensitivityVector(lambdaSDisp, true); // [F] receive dual variable for structure displacement TODO: need to specify second argument (true)

    // Checking if dXdSb has entries different from zero at the interior of the mesh
    this->postOp->checkVec(lambdaSDisp);

    if (lambdaSDisp.norm() == 0.0)
    {
      this->com->fprintf(stderr, "\n *** WARNING *** zero structure displacement dual variable \n\n");
      if(!ioData.sa.fsiFlag) exit(1);
    }
    if(actvar < 9) { dQdX = 0; dQdU =0; }
    else if(actvar == 9) { // Liftx
      dL = 0;  dForces = 0;  dMoments = 0;  dL[0] = 1;  dQdX = 0;  dQdU = 0;
      fsoGetDerivativeOfEffortsWRTStateAndMeshPositionAnalytical(ioData, dForces, dMoments, dL, *this->X, U, dQdX, dQdU);
    } else if(actvar == 10) { // Lifty
      dL = 0;  dForces = 0;  dMoments = 0;  dL[1] = 1;  dQdX = 0;  dQdU = 0;
      fsoGetDerivativeOfEffortsWRTStateAndMeshPositionAnalytical(ioData, dForces, dMoments, dL, *this->X, U, dQdX, dQdU);
    } else if(actvar == 11) { // Liftz
      dL = 0;  dForces = 0;  dMoments = 0;  dL[2] = 1;  dQdX = 0;  dQdU = 0;
      fsoGetDerivativeOfEffortsWRTStateAndMeshPositionAnalytical(ioData, dForces, dMoments, dL, *this->X, U, dQdX, dQdU);
    }
    fsoComputeAdjoint(ioData, *this->A, dQdX, dQdU, true);
    this->sendForceSensitivity(&lambdaX, false);

    lambdaSDisp = 0;

    iter++;
  }
  step++;
}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fso_on_AdjointSensitivityMesh(IoData &ioData, DistSVec<double,dim> &U)
{
  // Computing efforts (F: force, M: moment, L:LiftAndDrag)
  //    Vec3D F, M, L;
  //    fsoGetEfforts(ioData, X, U, F, M, L);

  double tag = 0.0;

  step = 0;
  dXdS = 0.0;
  dXdSb = 0.0;
  dAdS = 0.0;
  DFSPAR[0] = 0.0;
  DFSPAR[1] = 0.0;
  DFSPAR[2] = 0.0;
  actvar = 1;

  int numShapeVars = 100; // maximum shape vairbales is set 100.
  DistSVec<double,3> **dDdS = new DistSVec<double,3>*[numShapeVars];
  while(true) {

    // Reading derivative of the overall deformation
    bool readOK = domain->readVectorFromFile(this->input->shapederivatives, step, &tag, dXdSb);
    if(!readOK) break;
    this->com->fprintf(stderr, "\n ***** Surface derivatives of shape variable %d were read\n", step);

    // Checking if dXdSb has entries different from zero at the interior of the mesh
    this->postOp->checkVec(dXdSb);

    if (dXdSb.norm() == 0.0) {
      this->com->fprintf(stderr, "\n *** WARNING *** No Mesh Perturbation \n\n");
      if(!ioData.sa.fsiFlag) exit(1);
    }
    dDdS[step] = new DistSVec<double,3>(dXdSb);
    dXdSb = 0.0;

    step = step + 1;
  }
  if(step < numShapeVars) numShapeVars = step;

  Vec3D dForces(0.0), dMoments(0.0), dL(0.0);
  DistSVec<double,3> dQdX(*this->X);
  DistSVec<double,dim> dQdU(U);

  ////////////////////////////////////////////////////////////////////
  if(ioData.output.transient.dLiftDrag != NULL) {


    dQdX = 0;    dQdU = 0;     dL = 0;  dMoments = 0;  dForces = 0;

    double liftresults[3][numShapeVars];
    double dforceresults[3][numShapeVars];



    dL[0] = 1.0; dL[1] = 0.0;  dL[2] = 0.0;
    //probably computes dQdx and dQdU
    fsoGetDerivativeOfEffortsWRTStateAndMeshPositionAnalytical(
        ioData, dForces, dMoments, dL, *this->X, U, dQdX, dQdU);
    // computes lambdaX
    fsoComputeAdjoint(ioData, *this->A, dQdX, dQdU, false);
    for(step = 0; step<numShapeVars; ++step) {
      liftresults[0][step] = -1.0*(*dDdS[step]*lambdaX);
    }



    dL[0] = 0.0; dL[1] = 1.0;  dL[2] = 0.0;
    //probably computes dQdy and dQdU
    fsoGetDerivativeOfEffortsWRTStateAndMeshPositionAnalytical(
        ioData, dForces, dMoments, dL, *this->X, U, dQdX, dQdU);
    // computes lambdaX
    fsoComputeAdjoint(ioData, *this->A, dQdX, dQdU, false);
    for(step = 0; step<numShapeVars; ++step) {
      liftresults[1][step] = -1.0*(*dDdS[step]*lambdaX);
    }


    dL[0] = 0.0; dL[1] = 0.0;  dL[2] = 1.0;
    //probably computes dQdz and dQdU
    fsoGetDerivativeOfEffortsWRTStateAndMeshPositionAnalytical(
        ioData, dForces, dMoments, dL, *this->X, U, dQdX, dQdU);
    // computes lambdaX
    fsoComputeAdjoint(ioData, *this->A, dQdX, dQdU, false);
    for(step = 0; step<numShapeVars; ++step) {
      liftresults[2][step] = -1.0*(*dDdS[step]*lambdaX);
    }



    for(step = 0; step<numShapeVars; ++step) {
      Vec3D dL, L;
      L[0]=0;  L[1]=0;  L[2]=0;
      dL[0]=liftresults[0][step];
      dL[1]=liftresults[1][step];
      dL[2]=liftresults[2][step];
      this->output->writeDerivativeOfLiftDragToDisk(step,1,L,dL);


      //TODO HACK, fix this so that it actually outputs adjoint forces
      double sboom =0;
      Vec3D F, M;
      this->output->writeDerivativeOfForcesToDisk(step,1,F,dForces,M,dMoments,sboom,sboom);

    }






  }
  ////////////////////////////////////////////////////////////////////


  if(ioData.output.transient.dLiftx[0] != 0) {
    dQdX = 0;  dQdU = 0;  dL = 0;  dMoments = 0;  dForces = 0;  dL[0] = 1.0;
    fsoGetDerivativeOfEffortsWRTStateAndMeshPositionAnalytical(
        ioData, dForces, dMoments, dL, *this->X, U, dQdX, dQdU);
    fsoComputeAdjoint(ioData, *this->A, dQdX, dQdU, false);
    for(step = 0; step<numShapeVars; ++step) {
      double dlift = -1.0*(*dDdS[step]*lambdaX);
      this->output->writeDerivativeOfLiftxToDisk(dlift);
    }
  }
  if(ioData.output.transient.dLifty[0] != 0) {
    dQdX = 0;  dQdU = 0;  dL = 0;  dMoments = 0;  dForces = 0;  dL[1] = 1.0;
    fsoGetDerivativeOfEffortsWRTStateAndMeshPositionAnalytical(
        ioData, dForces, dMoments, dL, *this->X, U, dQdX, dQdU);
    fsoComputeAdjoint(ioData, *this->A, dQdX, dQdU, false);
    for(step = 0; step<numShapeVars; ++step) {
      double dlift = -1.0*(*dDdS[step]*lambdaX);
      this->output->writeDerivativeOfLiftyToDisk(dlift);
    }
  }
  if(ioData.output.transient.dLiftz[0] != 0) {
    dQdX = 0;  dQdU = 0;  dL = 0;  dMoments = 0;  dForces = 0;  dL[2] = 1.0;
    fsoGetDerivativeOfEffortsWRTStateAndMeshPositionAnalytical(
        ioData, dForces, dMoments, dL, *this->X, U, dQdX, dQdU);
    fsoComputeAdjoint(ioData, *this->A, dQdX, dQdU, false);
    for(step = 0; step<numShapeVars; ++step) {
      double dlift = -1.0*(*dDdS[step]*lambdaX);
      this->output->writeDerivativeOfLiftzToDisk(dlift);
    }
  }

  fsoPrintTextOnScreen("\n ***** Derivatives of mesh position were computed! \n");

}


//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fso_on_sensitivityMesh(bool isSparse, IoData &ioData, DistSVec<double,dim> &U)
{
    this->com->fprintf(stderr, "\n\033[93m fso_on_sensitivityMesh entered\033[00m \n\n");

    double tag = 0.0;

    step = 0;
    dXdS = 0.0;
    dXdSb = 0.0;
    dAdS = 0.0;
    DFSPAR[0] = 0.0;
    DFSPAR[1] = 0.0;
    DFSPAR[2] = 0.0;
    actvar = 1;

    while (true) {

      // Reading derivative of the overall deformation
      bool readOK = domain->readVectorFromFile(this->input->shapederivatives, step, &tag, dXdSb);
      if(!readOK) break;

      // Checking if dXdSb has entries different from zero at the interior of the mesh
      this->postOp->checkVec(dXdSb);

      if (dXdSb.norm() == 0.0) {
        this->com->fprintf(stderr, "\n *** WARNING *** No Mesh Perturbation \n\n");
        if(!ioData.sa.fsiFlag) exit(1);
      }

      this->com->fprintf(stderr, "\n ***** Shape variable %d\n", step);

      // Updating the mesh
      dXdS = *this->X;
      mms->solve(dXdSb, dXdS);
      dXdS -= *this->X;

      // Check that the mesh perturbation is propagated
      if (dXdS.norm() == 0.0) this->com->fprintf(stderr, "\n !!! WARNING !!! No Mesh Sensitivity Perturbation !!!\n\n");

      fsoComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U, false, isSparse);
      fsoComputeSensitivities(isSparse, ioData, "Derivatives with respect to the mesh position:", ioData.sa.sensoutput, *this->X, U);

      dXdSb = 0.0;
      step = step + 1;
    }

    fsoPrintTextOnScreen("\n ***** Derivatives of mesh position were computed! \n");
    this->com->fprintf(stderr, "\n\033[93m fso_on_sensitivityMeh exited\033[00m \n\n");

}


//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoComputeDerivativesOfFluxAndSolution(IoData &ioData,
                                         DistSVec<double,3> &X,
                                         DistVec<double> &A,
                                         DistSVec<double,dim> &U,
                                         bool isFSI,
                                         bool isSparse)
{

  dFdS = 0.0;

  // Derivative of the Flux, either analytical or semi-analytical
  if ( ioData.sa.scFlag == SensitivityAnalysis::ANALYTICAL ) {

    DistSVec<double,dim> dFdS2(dFdS), diff(dFdS);

    fsoAnalytical(isSparse, ioData, X, dXdS, A, U, dFdS);

  } else {
    fsoSemiAnalytical(ioData, X, A, U, dFdS);
  }

  //TODO BUGHUNT writing the linear solver right hand side to disk,
  //this is can than be postprocessed with sower and xp2exo
  if (ioData.sa.linsolverhs != NULL)
	  this->output->writeAnyVectorToDisk(ioData.sa.linsolverhs,step,step,dFdS);

  // Computing the derivative of the fluid variables
  // with respect to the optimization variables


  //get dUdS by Solving the following equation:  [dFdU]*[dUdS]=[dFdS]
  fsoLinearSolver(ioData, dFdS, dUdS,isFSI);

}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoComputeAdjoint(IoData &ioData, DistVec<double> &A, DistSVec<double,3> &dQdX, DistSVec<double,dim> &dQdU, bool isFSI)
{

//  lambdaU = 0.0;
  lambdaX = 0.0;
  DistSVec<double,3> rhs(dXdS);
  rhs = 0.0;

  // Derivative of the Flux, either analytical or semi-analytical
  if ( ioData.sa.scFlag != SensitivityAnalysis::ANALYTICAL ) {
    this->com->fprintf(stderr, " --- WARNING : only analytical adjoint sensitivities are available\n");
  }

  fsoAdjointLinearSolver(ioData, dQdU, lambdaU, isFSI);
  fsoApply_dFdXtranspose(A, lambdaU, rhs);
  rhs -= dQdX;
  if(isFSI) rhs -= dfaX;

  mms->setAdjointFlagOn();
  // solve for lambdaX
  lambdaX = *this->X;
  mms->solveAdjoint(rhs, lambdaX);
  lambdaX -= *this->X;

  mms->applyProjectorTranspose(lambdaX);

}

//------------------------------------------------------------------------------


template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoComputeAndSendForceSensitivities(
                                         bool isSparse,
                                         IoData &ioData, const char *fileName,
                                         DistSVec<double,3> &X, DistSVec<double,dim> &U)
{

    if (ioData.sa.scFlag == SensitivityAnalysis::FINITEDIFFERENCE ) {
      fsoGetDerivativeOfLoadFiniteDifference(ioData, X, dXdS, *this->A, U, dUdS, *load, *dLoad);
    } else {
      fsoGetDerivativeOfLoadAnalytical(isSparse, ioData, X, dXdS, U, dUdS, *load, *dLoad);
    }

  this->sendForceSensitivity(dLoad);

}

//------------------------------------------------------------------------------

template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoComputeSensitivities(
                                         bool isSparse,
                                         IoData &ioData,
                                         const char *mesage,
                                         const char *fileName,
                                         DistSVec<double,3> &X,
                                         DistSVec<double,dim> &U)
{

// Computing efforts (F: force, M: moment, L:LiftAndDrag)
  Vec3D F, M, L;
  fsoGetEfforts(ioData, X, U, F, M, L);

// Computing derivative of the efforts
  Vec3D dFds, dMds, dLdS;

  if ( ioData.sa.scFlag == SensitivityAnalysis::FINITEDIFFERENCE ){
    fsoGetDerivativeOfEffortsFiniteDifference(ioData, X, dXdS, *this->A, U, dUdS, dFds, dMds,dLdS);
  }
  else {
    fsoGetDerivativeOfEffortsAnalytical(isSparse, ioData, X, dXdS, U, dUdS, dFds, dMds, dLdS);//TODO uncomments
  }


  if (this->com->cpuNum() == 0) {
    outFile = fopen(fileName,"a+");
    if (outFile) {
      this->com->fprintf(outFile,mesage);
      this->com->fprintf(outFile,"\n");
      this->com->fprintf(outFile,"Fx= %16.13e \n",F[0]);
      this->com->fprintf(outFile,"Fy= %16.13e \n",F[1]);
      this->com->fprintf(outFile,"Fz= %16.13e \n",F[2]);
      this->com->fprintf(outFile,"dFx/ds= %16.13e \n",dFds[0]);
      this->com->fprintf(outFile,"dFy/ds= %16.13e \n",dFds[1]);
      this->com->fprintf(outFile,"dFz/ds= %16.13e \n",dFds[2]);
      this->com->fprintf(outFile,"Mx= %16.13e \n",M[0]);
      this->com->fprintf(outFile,"My= %16.13e \n",M[1]);
      this->com->fprintf(outFile,"Mz= %16.13e \n",M[2]);
      this->com->fprintf(outFile,"dMx/ds= %16.13e \n",dMds[0]);
      this->com->fprintf(outFile,"dMy/ds= %16.13e \n",dMds[1]);
      this->com->fprintf(outFile,"dMz/ds= %16.13e \n",dMds[2]);
      this->com->fprintf(outFile,"\n");
      fclose(outFile);
    }
  }

  double sboom = 0.0;
  double dSboom = 0.0;

/*  // Compute Flux norm and derivative
  if (ioData.sa.homotopy == SensitivityAnalysis::ON_HOMOTOPY)
     this->computeFunction(0,U,Flux,true);
  else
     this->computeFunction(0,U,Flux,false);

  double normF = Flux.norm();
  double normF2 = 0.5*normF*normF;
  double dnormF2;

  DFluxDs = 0;
  for (int i = 0; i < this->nPod; ++i)
    DFluxDs += (this->AJ[i])*dYdS[i];
  dnormF2 = (this->F)*DFluxDs;*/


  //
  // This function is simply writing to the disk.
  //
  //this->output->writeDerivativeOfFluxNormToDisk(step, actvar, normF2, dnormF2);
  this->output->writeDerivativeOfForcesToDisk(step, actvar, F, dFds, M, dMds, sboom, dSboom);
  this->output->writeDerivativeOfLiftDragToDisk(step, actvar, L, dLdS);

  //
  // This function is writing to the disk quantities of interest in binary files.
  // The possible quantities of interest include
  // - dUdS
  // - Derivative of Scalar Quantities: Density, Mach, Pressure, Temperature, TotPressure,
  //   NutTurb, EddyViscosity, VelocityScalar
  // - Derivative of Vector Quantities: VelocityVector, Displacement
  //
  //


  if (ioData.sa.dFdS_final != NULL)
	  this->output->writeAnyVectorToDisk(ioData.sa.dFdS_final,1,1,dFdS);

  this->output->writeBinaryDerivativeOfVectorsToDisk(
                step+1,         //iteration index
                actvar,         //variable type
                DFSPAR,
                *this->X,
                dXdS,
                U,               //state vector
                dUdS,            //derivative of state vector w.r.t shape variable
                this->timeState,
                this->A);

}

//------------------------------------------------------------------------------
