#include <FluidRomShapeOptimizationHandler.h>

#include <IoData.h>
#include <Domain.h>
#include <GeoSource.h>
#include <DistVector.h>
#include <MeshMotionSolver.h>
//#include <StructExc.h>
#include <MatVecProd.h>
//#include <KspPrec.h>
//#include <KspSolver.h>
#include <MemoryPool.h>

#include <cmath>

#include <iostream>
#include <string>

//------------------------------------------------------------------------------

template<int dim>
FluidRomShapeOptimizationHandler<dim>::FluidRomShapeOptimizationHandler
(
  IoData &ioData,
  GeoSource &geoSource,
  Domain *dom//,
) :
ImplicitPGTsDesc<dim>(ioData, geoSource, dom),
domain(dom),
dXb(dom->getNodeDistInfo()),
dXdS(dom->getNodeDistInfo()),
dXdSb(dom->getNodeDistInfo()),
Xc(dom->getNodeDistInfo()),
dAdS(dom->getNodeDistInfo()),
dFdS(dom->getNodeDistInfo()),
dFrdS(this->nPod),
//Frp(this->nPod),
//Frm(this->nPod),
dYdS(this->nPod),
dUdS(dom->getNodeDistInfo()),
DFluxDs(dom->getNodeDistInfo()),
p(dom->getNodeDistInfo()),
dPdS(dom->getNodeDistInfo()),
Flux(dom->getNodeDistInfo()),
FluxFD(dom->getNodeDistInfo()),
Pin(dom->getFaceDistInfo()),
Uc(dom->getNodeDistInfo()),
// Tests
Xplus(dom->getNodeDistInfo()),
Xminus(dom->getNodeDistInfo()),
dX(dom->getNodeDistInfo())
{

  // Initialize
  step = 0;

  if ( ioData.sa.scFlag == SensitivityAnalysis::FINITEDIFFERENCE ) {
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
  }  
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
  }

  Z = new DistSVec<double,3>(domain->getNodeDistInfo());

//  this->com->fprintf(stderr,"ioData.ts.implicit.mvp = %i\n",ioData.ts.implicit.mvp);
//  this->com->fprintf(stderr,"ioData.sa.mvp = %i\n",ioData.sa.mvp);
//  this->com->fprintf(stderr,"ioData.sa.lsSolver = %i\n",ioData.sa.lsSolver);

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

  MemoryPool mp;

  mvp->exportMemory(&mp);

  if (ioData.sa.sensMesh == SensitivityAnalysis::ON_SENSITIVITYMESH) {
    mms = new TetMeshMotionSolver(ioData.dmesh, geoSource.getMatchNodes(),domain,0);
  }

  length = ioData.output.transient.length;
  surface = ioData.output.transient.surface;

  numLocSub = domain->getNumLocSub();

  dXdS=0.0;
  dFdS=0.0;
  dFrdS=0.0;
  dUdS=0.0;
  dYdS=0.0;
  p=0.0;
  dPdS=0.0;

  FluxFD = 0.0;
  Flux = 0.0;

  reynolds0 = ioData.ref.reynolds_mu;
  kenergy0 = ioData.bc.inlet.kenergy;

}

//------------------------------------------------------------------------------

template<int dim>
FluidRomShapeOptimizationHandler<dim>::~FluidRomShapeOptimizationHandler()
{

  if (mms) delete mms;

  if (mvp) delete mvp;

}

//------------------------------------------------------------------------------

template<int dim>
void FluidRomShapeOptimizationHandler<dim>::fsoRestartBcFluxs(IoData &ioData)
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

//  this->com->fprintf(stderr, "\n\n FluidSensitivityAnalysis values: \n\n");
//  this->com->fprintf(stderr, "\n\n Inlet Alpha = %20.17e \n\n",ioData.bc.inlet.alpha);
//  this->com->fprintf(stderr, "\n\n Inlet Beta = %20.17e \n\n",ioData.bc.inlet.beta);
//  this->com->fprintf(stderr, "\n\n Outlet Alpha = %20.17e \n\n",ioData.bc.outlet.alpha);
//  this->com->fprintf(stderr, "\n\n Outlet Beta = %20.17e \n\n",ioData.bc.outlet.beta);
      
  if (ioData.problem.mode == ProblemData::NON_DIMENSIONAL) 
  {

    //
    // UH (08/10)
    // From IoDataCore::resetInputValues, the code does not allow NON_DIMENSIONAL.
    // The following lines will not be executed.
    //

    ioData.ref.mach = xmach;
    ioData.ref.dRe_mudMach = 0.0;

    if (ioData.sa.densFlag == false) {
      ioData.bc.inlet.density = 1.0;
      ioData.bc.outlet.density = 1.0;
      this->com->fprintf(stderr, "\n\n NonDim Density = %e \n\n",ioData.bc.inlet.density);
    }

    if (ioData.sa.pressFlag == false) {
      ioData.bc.inlet.pressure = 1.0/(gamma*ioData.ref.mach*ioData.ref.mach);
      ioData.bc.outlet.pressure = 1.0/(gamma*ioData.ref.mach*ioData.ref.mach);

      this->com->fprintf(stderr, "\n\n NonDim Pressure = %e \n\n",ioData.bc.inlet.pressure);
    }

    if (ioData.sa.apressFlag == false) {
      if (ioData.sa.pressFlag == false) {
        ioData.aero.pressure = ioData.bc.inlet.pressure;
        Pin = ioData.aero.pressure;
//        this->com->fprintf(stderr, "\n\n NonDim Internal Pressure = %e \n\n",ioData.aero.pressure);

        this->postOp->rstVarPostFcn(ioData);
        this->postOp->rstVar(ioData);

        this->spaceOp->rstFluxFcn(ioData);
 
        mvp->rstSpaceOp(ioData, this->varFcn, this->spaceOp, false);

        this->rstVarImplicitRomTsDesc(ioData);
      }
    }
  }
  else if (ioData.problem.mode == ProblemData::DIMENSIONAL) {
    
    //
    // Step 1: Re-scale all the parameters
    //

    ioData.eqs.fluidModel.pmin *= ioData.ref.rv.pressure;

    ioData.bc.inlet.density *= ioData.ref.rv.density;
    ioData.bc.inlet.pressure *= ioData.ref.rv.pressure;
    ioData.bc.inlet.temperature *= ioData.ref.rv.temperature;
    ioData.bc.inlet.nutilde *= ioData.ref.rv.nutilde;
    ioData.bc.inlet.kenergy *= ioData.ref.rv.kenergy;
    ioData.bc.inlet.eps *= ioData.ref.rv.epsilon;
    ioData.bc.outlet.density *= ioData.ref.rv.density;
    ioData.bc.outlet.pressure *= ioData.ref.rv.pressure;
    ioData.bc.outlet.temperature *= ioData.ref.rv.temperature;
    ioData.bc.outlet.nutilde *= ioData.ref.rv.nutilde;
    ioData.bc.outlet.kenergy *= ioData.ref.rv.kenergy;
    ioData.bc.outlet.eps *= ioData.ref.rv.epsilon;

    ioData.restart.etime *= ioData.ref.rv.time;
    ioData.restart.dt_nm1 *= ioData.ref.rv.time;
    ioData.restart.dt_nm2 *= ioData.ref.rv.time;
    ioData.restart.energy *= ioData.ref.rv.energy;
    ioData.bc.wall.temperature *= ioData.ref.rv.temperature;
    ioData.ts.timestep *= ioData.ref.rv.time;
    ioData.ts.maxTime *= ioData.ref.rv.time;
    ioData.rmesh.vx *= ioData.ref.rv.velocity;
    ioData.rmesh.vy *= ioData.ref.rv.velocity;
    ioData.rmesh.vz *= ioData.ref.rv.velocity;
    ioData.rmesh.ax *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.rmesh.ay *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.rmesh.az *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.rmesh.timestep *= ioData.ref.rv.time;

    for (int j=0; j<ioData.rmesh.num; j++){
      ioData.rmesh.vpts[j]->time     *= ioData.ref.rv.time;
      ioData.rmesh.vpts[j]->velocityX *= ioData.ref.rv.velocity;
      ioData.rmesh.vpts[j]->velocityY *= ioData.ref.rv.velocity;
      ioData.rmesh.vpts[j]->velocityZ *= ioData.ref.rv.velocity;
    }
    ioData.aero.pressure *= ioData.ref.rv.pressure;
    ioData.forced.timestep *= ioData.ref.rv.time;
    ioData.forced.frequency /= ioData.ref.rv.time;

    ioData.eqs.gravity_x *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.eqs.gravity_y *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.eqs.gravity_z *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.bc.hydro.depth *= ioData.ref.length;

    //
    // Step 2: Restart all values based on new inlet values
    // Step 2.1: Reset values in ioData.ref
    //

    ioData.ref.mach = ioData.bc.inlet.mach;
    ioData.ref.density = ioData.bc.inlet.density;
    ioData.ref.pressure = ioData.bc.inlet.pressure;
    double velocity = ioData.ref.mach * sqrt(gamma * (ioData.ref.pressure+Pstiff) / ioData.ref.density);
    ioData.ref.temperature = (ioData.ref.pressure + gamma*Pstiff)/ (ioData.ref.density * R);
    double viscosity = ioData.eqs.viscosityModel.sutherlandConstant * sqrt(ioData.ref.temperature) /
      (1.0 + ioData.eqs.viscosityModel.sutherlandReferenceTemperature/ioData.ref.temperature);
    ioData.ref.reynolds_mu = velocity * ioData.ref.length * ioData.ref.density / viscosity;

    if (ioData.eqs.type == EquationsData::NAVIER_STOKES)
      this->com->fprintf(stderr, "\n\n Reynolds = %e \n\n",ioData.ref.reynolds_mu);

    double dvelocitydMach = sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
    ioData.ref.dRe_mudMach = dvelocitydMach * ioData.ref.length * ioData.ref.density / viscosity;

    //
    // Step 2.2: Reset values in ioData.ref.rv
    //

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

    this->rstVarImplicitRomTsDesc(ioData);
    
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
void FluidRomShapeOptimizationHandler<dim>::fsoGetEfforts(IoData &ioData, 
                                  DistSVec<double,3> &X, DistSVec<double,dim> &U, Vec3D &F, Vec3D &M, Vec3D &L)
{

  int nSurfs = this->postOp->getNumSurf();

  Vec3D x0, force, moment;

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
void FluidRomShapeOptimizationHandler<dim>::fsoGetDerivativeOfEffortsFiniteDifference(IoData &ioData,
                                                          DistSVec<double,3> &X, DistSVec<double,3> &dX, DistVec<double> &A,
                                                          DistSVec<double,dim> &U, DistSVec<double,dim> &dU,
                                                          Vec3D &dForces, Vec3D &dMoments)
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
  double gamma = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
  double velocity = ioData.ref.mach * sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
  double dVelocity= sqrt(gamma * ioData.ref.pressure / ioData.ref.density)*DFSPAR[0];
  double dForce=2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*velocity*dVelocity;
  double dEnergy=2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*ioData.ref.length*velocity*dVelocity;

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
    F*=dForce;
    M*=dEnergy;
    dForces=dF+F;
    dMoments=dM+M;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void FluidRomShapeOptimizationHandler<dim>::fsoGetDerivativeOfEffortsAnalytical(IoData &ioData,
                                                          DistSVec<double,3> &X, DistSVec<double,3> &dX,
                                                          DistSVec<double,dim> &U, DistSVec<double,dim> &dU,
                                                          Vec3D &dForces, Vec3D &dMoments, Vec3D &dL)
{


  double gamma = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
  double velocity = ioData.ref.mach * sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
  double dVelocity= sqrt(gamma * ioData.ref.pressure / ioData.ref.density)*DFSPAR[0];
  double dForce=2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*velocity*dVelocity;
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

  x0[0] = ioData.output.transient.x0;
  x0[1] = ioData.output.transient.y0;
  x0[2] = ioData.output.transient.z0;

  this->spaceOp->computeGradP(X, *this->A, U);

  this->postOp->computeForceAndMoment(x0, X, U, 0, Fi, Mi, Fv, Mv);

  F = 0.0;
  M = 0.0;

  F = Fi[0] + Fv[0];
  M = Mi[0] + Mv[0];

  this->spaceOp->computeDerivativeOfGradP(X, dX, *this->A, dAdS, U, dU);

  this->postOp->computeDerivativeOfForceAndMoment(x0, X, dX, U, dU, DFSPAR, dFi, dMi, dFv, dMv);

  dF = 0.0;
  dM = 0.0;

  dF = dFi[0] + dFv[0];
  dM = dMi[0] + dMv[0];

  if (this->refVal->mode == RefVal::NON_DIMENSIONAL) {
    dF *= 2.0 * this->refVal->length*this->refVal->length / surface;
    dM *= 2.0 * this->refVal->length*this->refVal->length*this->refVal->length / (surface * length);
    dForces=dF;
    dMoments=dM;
  }
  else {
    dF *= this->refVal->force;  
    dM *= this->refVal->energy;
    F *= dForce;
    M *= dEnergy;
    dForces = dF+F;
    dMoments = dM+M;
  }

  dL[0] = dF[0]*cos(ioData.bc.inlet.alpha)*cos(ioData.bc.inlet.beta) +
          dF[1]*cos(ioData.bc.inlet.alpha)*sin(ioData.bc.inlet.beta) +
          dF[2]*sin(ioData.bc.inlet.alpha);

  dL[1] = -dF[0]*sin(ioData.bc.inlet.beta) + dF[1]*cos(ioData.bc.inlet.beta);

  dL[2] = -dF[0]*sin(ioData.bc.inlet.alpha)*cos(ioData.bc.inlet.beta) -
           dF[1]*sin(ioData.bc.inlet.alpha)*sin(ioData.bc.inlet.beta) +
           dF[2]*cos(ioData.bc.inlet.alpha);

}

//------------------------------------------------------------------------------

template<int dim>
void FluidRomShapeOptimizationHandler<dim>::fsoGetDerivativeOfLoadFiniteDifference(IoData &ioData, DistSVec<double,3> &X, DistSVec<double,3> &dX, DistVec<double> &A,
                                                                                                                            DistSVec<double,dim> &U, DistSVec<double,dim> &dU, DistSVec<double,3> &load, DistSVec<double,3> &dLoad)
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
void FluidRomShapeOptimizationHandler<dim>::fsoGetDerivativeOfLoadAnalytical(IoData &ioData, DistSVec<double,3> &X, DistSVec<double,3> &dX, 
                                             DistSVec<double,dim> &U, DistSVec<double,dim> &dU, DistSVec<double,3> &load, DistSVec<double,3> &dLoad)
{

  double gamma = ioData.eqs.fluidModel.gasModel.specificHeatRatio;

  double velocity = ioData.ref.mach * sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
  double dVelocity= sqrt(gamma * ioData.ref.pressure / ioData.ref.density)*DFSPAR[0];
  double dForce=2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*velocity*dVelocity;

  load=0.0;
  dLoad=0.0;

  this->spaceOp->computeGradP(X, *this->A, U);

  this->postOp->computeNodalForce(X, U, Pin, load);

  this->spaceOp->computeDerivativeOfGradP(X, dX, *this->A, dAdS, U, dU);

  this->postOp->computeDerivativeOfNodalForce(X, dX, U, dU, Pin, DFSPAR, dLoad);

  if (this->refVal->mode == RefVal::NON_DIMENSIONAL) {
    dLoad *= 2.0 * this->refVal->length*this->refVal->length / surface;
  }
  else {
    dLoad += (dForce / this->refVal->force) * load;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void FluidRomShapeOptimizationHandler<dim>::fsoSemiAnalytical
(
  IoData &ioData, 
  DistSVec<double,3> &X,
  DistVec<double> &A,
  DistSVec<double,dim> &U,
  Vec<double> &dF
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

  dFdS=1.0/(2.0*eps)*((*Fp)-(*Fm));

  this->projectVector(this->AJ, dFdS, dF);

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
void FluidRomShapeOptimizationHandler<dim>::fsoAnalytical
(IoData &ioData, DistSVec<double,3> &X, DistVec<double> &A, DistSVec<double,dim> &U, Vec<double> &dFrds)
{
 
  //
  // Computing the normal, derivative of the normal and of the control volume
  //
  this->geoState->computeDerivatives(X, dXdS, this->bcData->getVelocityVector(), this->bcData->getDerivativeOfVelocityVector(), dAdS);

  //
  // Computing the derivatives of the boundary fluxes
  //
  this->bcData->initializeSA(ioData, X, dXdS, DFSPAR[0], DFSPAR[1], DFSPAR[2]);

  //
  // Computing the partial derivative of the flux with respect to the variables
  //
  this->spaceOp->computeDerivativeOfResidual(X, dXdS, A, dAdS, U, DFSPAR[0], Flux, dFdS, this->timeState);

  this->spaceOp->applyBCsToDerivativeOfResidual(U, dFdS);

  // Compute sensitivity of ROM residual (from HDM residual)
  this->projectVector(this->AJ, dFdS, dFrds);

  //
  // ROM contribution to dFdS
  //
  // TODO: dFrds_ij = d2Rdwds_kpj*Phi_pi*R_k + (dRdw)_kp*Phi_pi*(dRdS)_kj
  // Start with second term ONLY!
  
}

//------------------------------------------------------------------------------

template<int dim>
void FluidRomShapeOptimizationHandler<dim>::computeAJ(int it, DistSVec<double, dim> &Q)  {

  mvp->evaluate(it, *this->X, *this->A, Q, this->F);

  for (int iPod = 0; iPod < this->nPod; iPod++)
    mvp->apply(this->pod[iPod], this->AJ[iPod]);

}

//------------------------------------------------------------------------------

template<int dim>
void FluidRomShapeOptimizationHandler<dim>::fsoSetUpLinearSolver(IoData &ioData, DistSVec<double,3> &X, DistVec<double> &A, DistSVec<double,dim> &U)
{

  fsoRestartBcFluxs(ioData);
  this->geoState->reset(X);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), X, A);
  this->bcData->update(X);

  if (ioData.sa.homotopy == SensitivityAnalysis::ON_HOMOTOPY)
    this->computeFullResidual(0,U,false,NULL,true);
  else
    this->computeFullResidual(0,U,false,NULL,false);
  //this->computeAJ(0,U);
  computeAJ(0,U);

  // Set up reduced linear system (for PG: Psi = dFdw*Phi)
  // d(Psi^T*F)_i/dy_j = (Psi^T*(dFdw)*Phi)_ij + F_k*d(Psi_ki)/dy_j
  // TODO: Add second term via finite differences
  double* jactmp = new double [this->nPod * this->nPod];
  transMatMatProd(this->AJ,this->AJ,jactmp);
  this->jac.zero();
  for (int iRow = 0; iRow < this->nPod; ++iRow)
    for (int iCol = 0; iCol < this->nPod; ++iCol)
      this->jac[iRow][iCol] = jactmp[iRow + iCol * this->nPod];
  delete[] jactmp;

}

//------------------------------------------------------------------------------

template<int dim>
void FluidRomShapeOptimizationHandler<dim>::fsoPrintTextOnScreen(const char *Text)
{
   this->com->fprintf(stderr, Text);
}

//------------------------------------------------------------------------------

template<int dim>
void FluidRomShapeOptimizationHandler<dim>::fsoOutput1D(const char *fileName, DistVec<double> &V)
{

outFile = fopen(fileName,"w");


#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    for (int j=0; j<V(iSub).size(); ++j)
        fprintf(outFile," %20.17e \n",V(iSub)[j]);
//        fprintf(outFile," %9.6e \n",V(iSub)[j]);


fclose(outFile);

}

//------------------------------------------------------------------------------

template<int dim>
void FluidRomShapeOptimizationHandler<dim>::fsoOutput3D(const char *fileName, DistSVec<double,3> &V)
{

outFile = fopen(fileName,"w");

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    for (int j=0; j<V(iSub).size(); ++j) {
      for (int k=0; k<3; ++k)
        fprintf(outFile," %20.17e ",V(iSub)[j][k]);
//        fprintf(outFile," %9.6e ",V(iSub)[j][k]);
      fprintf(outFile,"\n");

    }
  }

fclose(outFile);

}

//------------------------------------------------------------------------------

template<int dim>
void FluidRomShapeOptimizationHandler<dim>::fsoOutputDimD(const char *fileName, DistSVec<double,dim> &V)
{


outFile = fopen(fileName,"w");

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub ; ++iSub) {
    for (int j=0; j<V(iSub).size(); ++j) {
      for (int k=0; k<dim; ++k)
        fprintf(outFile," %20.17e ",V(iSub)[j][k]);
//        fprintf(outFile," %9.6e ",V(iSub)[j][k]);
      fprintf(outFile,"\n");
    }

  }

fclose(outFile);

}

//------------------------------------------------------------------------------

template<int dim>
void FluidRomShapeOptimizationHandler<dim>::fsoInitialize(IoData &ioData, DistSVec<double,dim> &U)
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

  this->checkLocalRomStatus(U,0);
  this->jac.setNewSize(this->nPod,this->nPod);
  dFrdS.resize(this->nPod);
  dYdS.resize(this->nPod);
  fsoSetUpLinearSolver(ioData,*this->X, *this->A, U);
}

//------------------------------------------------------------------------------

template<int dim>
int FluidRomShapeOptimizationHandler<dim>::fsoHandler(IoData &ioData, DistSVec<double,dim> &U)
{

  // xmach      -  Mach number
  // alpha      -  pitch angle
  // teta       -  yaw angle
  // DFSPAR(1)  -  Mach number differential
  // DFSPAR(2)  -  angle of attack differential
  // DFSPAR(3)  -  yaw angle differential

  // Start basic timer
  double MyLocalTimer = -this->timer->getTime();

  double dtLeft = 0.0;
  this->computeTimeStep(1, &dtLeft, U);
  this->computeMeshMetrics();
  this->updateStateVectors(U);

  fsoSetUpLinearSolver(ioData,*this->X, *this->A, U);

  if (ioData.sa.sensMesh == SensitivityAnalysis::ON_SENSITIVITYMESH) {

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

      if (dXdSb.norm() == 0.0)
      {
        this->com->fprintf(stderr, "\n *** ERROR *** No Mesh Perturbation \n\n");
        exit(1);
      }

      this->com->fprintf(stderr, "\n ***** Shape variable %d\n", step);

      // Updating the mesh
      dXdS = *this->X;
      mms->solve(dXdSb, dXdS);
      dXdS -= *this->X;

      // Check that the mesh perturbation is propagated
      if (dXdS.norm() == 0.0)
      {
        this->com->fprintf(stderr, "\n !!! WARNING !!! No Mesh Perturbation !!!\n\n");
      }

      fsoComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U);
  
      fsoComputeSensitivities(ioData, "Derivatives with respect to the mesh position:", ioData.sa.sensoutput, *this->X, U);

      dXdSb = 0.0;

      step = step + 1;
    }

    fsoPrintTextOnScreen("\n ***** Derivatives with respect to the mesh position were computed! \n");

  }

  if (ioData.sa.sensMach == SensitivityAnalysis::ON_SENSITIVITYMACH) {

    dXdS = 0.0;
    dAdS = 0.0;
    DFSPAR[0] = 1.0;
    DFSPAR[1] = 0.0;
    DFSPAR[2] = 0.0;
    actvar = 2;

    fsoComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U);

    fsoComputeSensitivities(ioData, "Derivatives with respect to the Mach number:", ioData.sa.sensoutput, *this->X, U);

    fsoPrintTextOnScreen("\n ***** Derivatives with respect to the Mach number were computed! \n");

    step = step + 1;
  }

  if (ioData.sa.sensAlpha == SensitivityAnalysis::ON_SENSITIVITYALPHA) {

    dXdS = 0.0;
    dAdS = 0.0;
    DFSPAR[0] = 0.0;
    DFSPAR[1] = 1.0;
    DFSPAR[2] = 0.0;
    actvar = 3;

    if (!ioData.sa.angleRad) 
      ioData.sa.eps *= acos(-1.0) / 180.0;

    fsoComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U);

    fsoComputeSensitivities(ioData, "Derivatives with respect to the angle of attack:", ioData.sa.sensoutput, *this->X, U);

    fsoPrintTextOnScreen("\n ***** Derivatives with respect to the angle of attack were computed! \n");

    step = step + 1;

    if (!ioData.sa.angleRad)
      ioData.sa.eps /= acos(-1.0) / 180.0;
  }

  if (ioData.sa.sensBeta == SensitivityAnalysis::ON_SENSITIVITYBETA) {

    dXdS = 0.0;
    dAdS = 0.0;
    DFSPAR[0] = 0.0;
    DFSPAR[1] = 0.0;
    DFSPAR[2] = 1.0;
    actvar = 4;

    if (!ioData.sa.angleRad)
      ioData.sa.eps *= acos(-1.0) / 180.0;

    fsoComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U);

    fsoComputeSensitivities(ioData, "Derivatives with respect to the yaw angle:", ioData.sa.sensoutput, *this->X, U);

    fsoPrintTextOnScreen("\n ***** Derivatives with respect to the yaw angle were computed! \n");

    step = step + 1;

    if (!ioData.sa.angleRad)
      ioData.sa.eps /= acos(-1.0) / 180.0;
  }

  bool lastIt = true;
//  this->outputToDisk(ioData, &lastIt, 0, 0, 0, 0, dtLeft, U); 
  this->outputPositionVectorToDisk(U);

  this->output->closeAsciiFiles();

  this->com->barrier();
  MyLocalTimer += this->timer->getTime();
  if (this->com->cpuNum() == 0)
  {
    std::cout << "\n *** FluidSensityAnalysisHandler::fsoHandler >> Exit";
    std::cout << " (" << MyLocalTimer << " s)";
    std::cout << "\n\n";
  }

  return -1;

}

//------------------------------------------------------------------------------

template<int dim>
void FluidRomShapeOptimizationHandler<dim>::fsoComputeDerivativesOfFluxAndSolution(IoData &ioData, DistSVec<double,3> &X, DistVec<double> &A, DistSVec<double,dim> &U)
{

  dFdS = 0.0;
  dFrdS = 0.0;

  // Derivative of the Flux, either analytical or semi-analytical
  if ( ioData.sa.scFlag == SensitivityAnalysis::ANALYTICAL )
    fsoAnalytical(ioData, X, A, U, dFrdS);
  else
    fsoSemiAnalytical(ioData, X, A, U, dFrdS);

  // Computing the derivative of the fluid variables 
  // with respect to the fsoimization variables
  if (ioData.sa.lsSolver == SensitivityAnalysis::QR) {
     double** lsCoefficients = new double*[1];
     lsCoefficients[0] = new double[this->nPod];
     RefVec<DistSVec<double,dim> > residualRef2(dFdS);
     this->parallelRom->parallelLSMultiRHS(this->AJ,residualRef2,this->nPod,1,lsCoefficients);
     for (int i=0; i< this->nPod; ++i)
        dYdS[i] = -lsCoefficients[0][i];
     delete[] lsCoefficients[0];
     delete[] lsCoefficients;
  } else if (ioData.sa.lsSolver == SensitivityAnalysis::NORMAL_EQUATIONS) {
     this->solveLinearSystem(0, dFrdS, dYdS);
     dYdS*=(-1.0);
  }

//  dYdS*=0.0;
//  dYdS[1]=1.0;

  this->com->fprintf(stderr,"dYdS (before) = \n");
  for (int i=0; i<dYdS.len; ++i)
    this->com->fprintf(stderr,"%20.16f  \n",dYdS[i]);

  DistSVec<double,dim> tmp(domain->getNodeDistInfo());
  tmp = 0;
  for (int i = 0; i < this->nPod; ++i)
    tmp += (this->AJ[i])*dYdS[i];
  tmp += dFdS;
  this->com->fprintf(stderr,"||J*Phi*dYdS - dRdS|| = %20.16f\n",tmp.norm());
}

//------------------------------------------------------------------------------

template<int dim>
void FluidRomShapeOptimizationHandler<dim>::fsoComputeSensitivities(IoData &ioData, const char *mesage, const char *fileName, DistSVec<double,3> &X, DistSVec<double,dim> &U)
{

// Computing efforts (F: force, M: moment, L:LiftAndDrag)
  Vec3D F, M, L;
  fsoGetEfforts(ioData, X, U, F, M, L);

// Computing derivative of the efforts
  Vec3D dFds, dMds, dLds;

  this->expandVector(dYdS,dUdS);
  if ( ioData.sa.scFlag == SensitivityAnalysis::FINITEDIFFERENCE ){
    fsoGetDerivativeOfEffortsFiniteDifference(ioData, X, dXdS, *this->A, U, dUdS, dFds, dMds);
    dLds[0] = dFds[0]*cos(ioData.bc.inlet.alpha)*cos(ioData.bc.inlet.beta) +
              dFds[1]*cos(ioData.bc.inlet.alpha)*sin(ioData.bc.inlet.beta) +
              dFds[2]*sin(ioData.bc.inlet.alpha);

    dLds[1] = -dFds[0]*sin(ioData.bc.inlet.beta) + dFds[1]*cos(ioData.bc.inlet.beta);

    dLds[2] = -dFds[0]*sin(ioData.bc.inlet.alpha)*cos(ioData.bc.inlet.beta) -
               dFds[1]*sin(ioData.bc.inlet.alpha)*sin(ioData.bc.inlet.beta) +
               dFds[2]*cos(ioData.bc.inlet.alpha);
    this->com->fprintf(stderr,"dLds = %20.16f, %20.16f, %20.16f\n",dLds[0],dLds[1],dLds[2]);

  } else
    fsoGetDerivativeOfEffortsAnalytical(ioData, X, dXdS, U, dUdS, dFds, dMds, dLds);

  if ((!ioData.sa.angleRad) && (DFSPAR[1] || DFSPAR[2])) {
    dFds *= acos(-1.0) / 180.0;
    dMds *= acos(-1.0) / 180.0;
    dLds *= acos(-1.0) / 180.0;
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

  // Compute Flux norm and derivative
  double normF = (this->F).norm();
  double normF2 = 0.5*normF*normF;
  double dnormF2;

  this->com->fprintf(stderr,"normF = %20.16f\n",normF);
  this->com->fprintf(stderr,"normF2 = %20.16f\n",normF2);

  DFluxDs = 0;
  for (int i = 0; i < this->nPod; ++i)
    DFluxDs += (this->AJ[i])*dYdS[i];
  dnormF2 = (this->F)*DFluxDs; 

  //
  // This function is simply writing to the disk.
  //
  this->output->writeDerivativeOfFluxNormToDisk(step, actvar, normF2, dnormF2);
  this->output->writeDerivativeOfForcesToDisk(step, actvar, F, dFds, M, dMds, sboom, dSboom);
  this->output->writeDerivativeOfLiftDragToDisk(step, actvar, L, dLds); 
 
  //
  // This function is writing to the disk quantities of interest in binary files.
  // The possible quantities of interest include
  // - dUdS
  // - Derivative of Scalar Quantities: Density, Mach, Pressure, Temperature, TotPressure,
  //   NutTurb, EddyViscosity, VelocityScalar
  // - Derivative of Vector Quantities: VelocityVector, Displacement
  //
  //
  this->output->writeBinaryDerivativeOfVectorsToDisk(step+1, actvar, DFSPAR, *this->X, dXdS, U, dUdS, this->timeState, this->A);

}

//------------------------------------------------------------------------------


