#include <EmbeddedFluidShapeOptimizationHandler.h>

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

//------------------------------------------------------------------------------

template<int dim>
EmbeddedFluidShapeOptimizationHandler<dim>::EmbeddedFluidShapeOptimizationHandler
( IoData &ioData,
  GeoSource &geoSource,
  Domain *dom 
) :
  ImplicitEmbeddedCoupledTsDesc<dim>(ioData, geoSource, dom),
  dXdS(dom->getNodeDistInfo()),
  Flux(dom->getNodeDistInfo()),
FluxFD(dom->getNodeDistInfo()),
  dFdS(dom->getNodeDistInfo()), 
  dUdS(dom->getNodeDistInfo()),
   Pin(dom->getFaceDistInfo()),
//debug
dFdS_debug(dom->getNodeDistInfo()),
difference(dom->getNodeDistInfo())
{

  step = 0;

  if ( ioData.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_ ) {
        load = new DistSVec<double,3>(dom->getNodeDistInfo());
       dLoad = new DistSVec<double,3>(dom->getNodeDistInfo());
    dLoadref = new DistSVec<double,3>(dom->getNodeDistInfo());
  } else {
        load = 0;
       dLoad = 0;
    dLoadref = 0;
  }

  X_ = new DistSVec<double,3>(dom->getNodeDistInfo());
  A_ = new DistVec<double>(dom->getNodeDistInfo());   

  Lp = new DistSVec<double,3>(dom->getNodeDistInfo());
  Lm = new DistSVec<double,3>(dom->getNodeDistInfo());

  Fp = new DistSVec<double,dim>(dom->getNodeDistInfo());
  Fm = new DistSVec<double,dim>(dom->getNodeDistInfo());

  Up = new DistSVec<double,dim>(dom->getNodeDistInfo());
  Um = new DistSVec<double,dim>(dom->getNodeDistInfo());   

  if (ioData.sa.homotopy == SensitivityAnalysis::ON_HOMOTOPY) {

    if (ioData.ts.implicit.mvp == ImplicitData::H2){
      mvp = new MatVecProdH2<dim,double,dim>(ioData, this->varFcn, this->timeState, this->spaceOp, dom, this->geoState);
    } else {
      mvp = new MatVecProdFD<dim,dim>(ioData.ts.implicit, this->timeState, this->geoState, this->spaceOp, dom, ioData);
    }

  } else {

    if (ioData.ts.implicit.mvp == ImplicitData::H2) {
      mvp = new MatVecProdH2<dim,double,dim>(ioData, this->varFcn, 0, this->spaceOp, dom, this->geoState);
    } else {
      mvp = new MatVecProdFD<dim,dim>(ioData.ts.implicit,  0, this->geoState, this->spaceOp, dom, ioData);
    }

  } 

  pc = ImplicitEmbeddedTsDesc<dim>::template 
    createPreconditioner<dim>(ioData.sa.ksp.pc, this->domain);

  ksp = this->createKrylovSolver(this->getVecInfo(), ioData.sa.ksp, mvp, pc, this->com);

  typename MatVecProd<dim,dim>::_fsi fsi = {
    this->distLSS,
   &this->nodeTag,
    this->riemann,
    this->linRecAtInterface,
    this->viscSecOrder,
   &this->Wtemp,
    this->riemannNormal,
    this->ghostPoints,
  };
  mvp->AttachStructure(fsi);

   length = ioData.output.transient.length;
  surface = ioData.output.transient.surface;

  numLocSub = dom->getNumLocSub();

  dFdS = 0.0;

    Flux = 0.0;
  FluxFD = 0.0;

  reynolds0 = ioData.ref.reynolds_mu;
   kenergy0 = ioData.bc.inlet.kenergy;

  steadyTol = 1.0;
  
  int sp = strlen(ioData.input.prefix) + 1;
  dXdSb_file = new char[sp + strlen(ioData.input.shapederivatives)];
  sprintf(dXdSb_file, "%s%s", ioData.input.prefix, ioData.input.shapederivatives);

}

//------------------------------------------------------------------------------

template<int dim>
EmbeddedFluidShapeOptimizationHandler<dim>::~EmbeddedFluidShapeOptimizationHandler()
{

  if (mvp)      delete mvp;
  if (pc)       delete pc;
  if (ksp)      delete ksp;
  if (load)     delete load;
  if (dLoad)    delete dLoad;
  if (dLoadref) delete dLoadref;

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidShapeOptimizationHandler<dim>::fsoInitialize(IoData &ioData, DistSVec<double,dim> &U)
{

  this->output->openAsciiFiles();

  // Reseting the configuration control of the geometry datas
  this->geoState->resetConfigSA();

  // ?????
  if (this->com->cpuNum() == 0) {
    outFile = fopen(ioData.sa.sensoutput, "w");
    if (outFile) 
      fclose(outFile);
  }
  // ?????

  xmach  = ioData.sa.machref;
  alprad = ioData.sa.alpharef;
  teta   = ioData.sa.betaref;

  double dtLeft = 0.0;
  this->computeTimeStep(1, &dtLeft, U); //?

  int it = this->getInitialIteration();
  double dt = this->computeTimeStep(it, &dtLeft, U, -2.0);

  this->updateFarfieldCoeffs(dt);
  this->interpolatePositionVector(dt, dtLeft);
  this->computeMeshMetrics();
  this->updateStateVectors(U);  

  // Setting up the linear solver
  fsoSetUpLinearSolver(ioData, *this->X, *this->A, U, dFdS);

}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidShapeOptimizationHandler<dim>::fsoSetUpLinearSolver
(
 IoData &ioData, 
 DistSVec<double,3> &X, 
 DistVec<double> &A, 
 DistSVec<double,dim> &U,
 DistSVec<double,dim> &dFdS
){

  fsoRestartBcFluxs(ioData);

  this->geoState->reset(X);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), X, A);
  this->bcData->update(X);

  this->spaceOp->computeResidual(X, A, U, 
											*this->Wstarij, *this->Wstarji, *this->Wextij, this->distLSS,
                                 this->linRecAtInterface, this->viscSecOrder, 
				 this->nodeTag, 
				 FluxFD, 
				 this->riemann, this->riemannNormal, 1, this->ghostPoints);

  if (ioData.sa.homotopy == SensitivityAnalysis::ON_HOMOTOPY)
    this->timeState->add_dAW_dt(1, *this->geoState, A, U, FluxFD, this->distLSS); 

  this->spaceOp->applyBCsToResidual(U, FluxFD, this->distLSS);

  mvp->evaluate(0, X, A, U, FluxFD);  

  DistMat<double,dim> *_pc = dynamic_cast<DistMat<double,dim> *>(pc);

  if (_pc) {

    MatVecProdFD<dim,dim>           *mvpfd = dynamic_cast<MatVecProdFD<dim,dim> *>(mvp);
    MatVecProdH2<dim,MatScalar,dim> *mvph2 = dynamic_cast<MatVecProdH2<dim,double,dim> *>(mvp);

    if (mvpfd || mvph2) {

      this->spaceOp->computeJacobian(X, A, U, 
				     this->distLSS, this->nodeTag, 
				     this->riemann, this->riemannNormal,
				     this->ghostPoints, *_pc, this->timeState);

      if (ioData.sa.homotopy == SensitivityAnalysis::ON_HOMOTOPY){
	this->timeState->addToJacobian(A, *_pc, U);
      }

      this->spaceOp->applyBCsToJacobian(U, *_pc, this->distLSS);      
    }

  }

  pc->setup();

  this->spaceOp->computeResidual(X, A, U, 
											*this->Wstarij, *this->Wstarji, *this->Wextij, this->distLSS,
                                 this->linRecAtInterface, this->viscSecOrder, 
				 this->nodeTag, 
				 Flux, 
				 this->riemann, this->riemannNormal, 1, this->ghostPoints, false);
}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidShapeOptimizationHandler<dim>::fsoRestartBcFluxs(IoData &ioData)
{

  double gamma  = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
  double R      = ioData.eqs.fluidModel.gasModel.idealGasConstant;
  double Pstiff = ioData.eqs.fluidModel.gasModel.pressureConstant;

  ioData.bc.inlet.mach  = xmach;
  ioData.bc.inlet.alpha = alprad;
  ioData.bc.inlet.beta  = teta;
 
  ioData.bc.outlet.mach  = xmach;
  ioData.bc.outlet.alpha = alprad;
  ioData.bc.outlet.beta  = teta;
  
  if (ioData.problem.mode == ProblemData::NON_DIMENSIONAL) {

    this->com->fprintf(stderr, "\n NON-DIMENSIONAL simulation not allowed \n");

  } else {

    // Step 1: From NonDimensial to Dimensional values
    ioData.eqs.fluidModel.pmin *= ioData.ref.rv.pressure;

    ioData.bc.inlet.density     *= ioData.ref.rv.density;
    ioData.bc.inlet.pressure    *= ioData.ref.rv.pressure;
    ioData.bc.inlet.temperature *= ioData.ref.rv.temperature;
    ioData.bc.inlet.nutilde     *= ioData.ref.rv.nutilde;
    ioData.bc.inlet.kenergy     *= ioData.ref.rv.kenergy;
    ioData.bc.inlet.eps         *= ioData.ref.rv.epsilon;

    ioData.bc.outlet.density     *= ioData.ref.rv.density;
    ioData.bc.outlet.pressure    *= ioData.ref.rv.pressure;
    ioData.bc.outlet.temperature *= ioData.ref.rv.temperature;
    ioData.bc.outlet.nutilde     *= ioData.ref.rv.nutilde;
    ioData.bc.outlet.kenergy     *= ioData.ref.rv.kenergy;
    ioData.bc.outlet.eps         *= ioData.ref.rv.epsilon;

    ioData.restart.etime  *= ioData.ref.rv.time;
    ioData.restart.dt_nm1 *= ioData.ref.rv.time;
    ioData.restart.dt_nm2 *= ioData.ref.rv.time;
    ioData.restart.energy *= ioData.ref.rv.energy;

    ioData.bc.wall.temperature *= ioData.ref.rv.temperature;

    ioData.ts.timestep *= ioData.ref.rv.time;
    ioData.ts.maxTime  *= ioData.ref.rv.time;

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
    ioData.aero.pressure *= ioData.ref.rv.pressure;

    ioData.forced.timestep  *= ioData.ref.rv.time;
    ioData.forced.frequency /= ioData.ref.rv.time;

    ioData.eqs.gravity_x *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.eqs.gravity_y *= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.eqs.gravity_z *= ioData.ref.rv.velocity / ioData.ref.rv.time;

    ioData.bc.hydro.depth *= ioData.ref.length;

    // Step 2: Restart all values based on new inlet values
    //
    // Step 2.1: Reset values in ioData.ref
    ioData.ref.mach     = ioData.bc.inlet.mach;
    ioData.ref.density  = ioData.bc.inlet.density;
    ioData.ref.pressure = ioData.bc.inlet.pressure;

    double velocity = ioData.ref.mach * sqrt(gamma * (ioData.ref.pressure+Pstiff) / ioData.ref.density);

    ioData.ref.temperature = (ioData.ref.pressure + gamma*Pstiff)/ (ioData.ref.density * R);

    double viscosity = ioData.eqs.viscosityModel.sutherlandConstant * sqrt(ioData.ref.temperature) /
                       (1.0 + ioData.eqs.viscosityModel.sutherlandReferenceTemperature/ioData.ref.temperature);

    ioData.ref.reynolds_mu = velocity * ioData.ref.length * ioData.ref.density / viscosity;

    double dvelocitydMach = sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
    ioData.ref.dRe_mudMach = dvelocitydMach * ioData.ref.length * ioData.ref.density / viscosity;

    // Step 2.2: Reset values in ioData.ref.rv
    ioData.ref.rv.mode           = RefVal::DIMENSIONAL;
    ioData.ref.rv.density        = ioData.ref.density;
    ioData.ref.rv.velocity       = velocity;
    ioData.ref.rv.pressure       = ioData.ref.density * velocity*velocity;
    ioData.ref.rv.temperature    = gamma*(gamma - 1.0) * ioData.ref.mach*ioData.ref.mach * (ioData.ref.pressure+Pstiff)/(R*ioData.ref.density); //?
    ioData.ref.rv.viscosity_mu   = viscosity;
    ioData.ref.rv.nutilde        = viscosity / ioData.ref.density;
    ioData.ref.rv.kenergy        = velocity*velocity;
    ioData.ref.rv.epsilon        = velocity*velocity*velocity / ioData.ref.length;
    ioData.ref.rv.time           = ioData.ref.length / velocity;
    ioData.ref.rv.force          = ioData.ref.density * velocity*velocity * ioData.ref.length*ioData.ref.length;
    ioData.ref.rv.energy         = ioData.ref.density * velocity*velocity * ioData.ref.length*ioData.ref.length*ioData.ref.length;
    ioData.ref.rv.power          = ioData.ref.density * velocity*velocity*velocity * ioData.ref.length*ioData.ref.length;
    ioData.ref.rv.tvelocity      = velocity / ioData.aero.displacementScaling;
    ioData.ref.rv.tforce         = ioData.ref.rv.force / ioData.aero.forceScaling;
    ioData.ref.rv.tpower         = ioData.ref.rv.power / ioData.aero.powerScaling;
    ioData.ref.rv.dvelocitydMach = dvelocitydMach;
    ioData.ref.rv.dtimedMach     = -ioData.ref.length / (velocity * velocity) * dvelocitydMach;

    ioData.eqs.fluidModel.pmin /= ioData.ref.rv.pressure;

    // Step 2.3: Back to dimensionless bc.inlet
    ioData.bc.inlet.density     /= ioData.ref.rv.density;
    ioData.bc.inlet.pressure    /= ioData.ref.rv.pressure;
    ioData.bc.inlet.temperature /= ioData.ref.rv.temperature;
    ioData.bc.inlet.nutilde     /= ioData.ref.rv.nutilde;
    ioData.bc.inlet.kenergy     /= ioData.ref.rv.kenergy;
    ioData.bc.inlet.eps         /= ioData.ref.rv.epsilon;

    // Step 2.4: Back to dimensionless bc.outlet
    ioData.bc.outlet.density     /= ioData.ref.rv.density;
    ioData.bc.outlet.pressure    /= ioData.ref.rv.pressure;
    ioData.bc.outlet.temperature /= ioData.ref.rv.temperature;
    ioData.bc.outlet.nutilde     /= ioData.ref.rv.nutilde;
    ioData.bc.outlet.kenergy     /= ioData.ref.rv.kenergy;
    ioData.bc.outlet.eps         /= ioData.ref.rv.epsilon;

    ioData.restart.etime  /= ioData.ref.rv.time;
    ioData.restart.dt_nm1 /= ioData.ref.rv.time;
    ioData.restart.dt_nm2 /= ioData.ref.rv.time;
    ioData.restart.energy /= ioData.ref.rv.energy;

    ioData.bc.wall.temperature /= ioData.ref.rv.temperature;

    ioData.linearizedData.stepsize = ioData.ts.timestep;

    ioData.ts.timestep /= ioData.ref.rv.time;             // Problem in RigidRollMeshMotionHandler
    ioData.ts.maxTime  /= ioData.ref.rv.time;             // Problem in RigidRollMeshMotionHandler

    ioData.rmesh.vx       /= ioData.ref.rv.velocity;      // Problem in RigidMeshMotionHandler
    ioData.rmesh.vy       /= ioData.ref.rv.velocity;
    ioData.rmesh.vz       /= ioData.ref.rv.velocity;
    ioData.rmesh.ax       /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.rmesh.ay       /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.rmesh.az       /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.rmesh.timestep /= ioData.ref.rv.time;          // Problem in AccMeshMotionHandler
    
    for (int j=0; j<ioData.rmesh.num; j++){
      ioData.rmesh.vpts[j]->time      /= ioData.ref.rv.time;
      ioData.rmesh.vpts[j]->velocityX /= ioData.ref.rv.velocity;
      ioData.rmesh.vpts[j]->velocityY /= ioData.ref.rv.velocity;
      ioData.rmesh.vpts[j]->velocityZ /= ioData.ref.rv.velocity;
    }

    ioData.aero.pressure    /= ioData.ref.rv.pressure;
    ioData.forced.timestep  /= ioData.ref.rv.time;        // Problem in ForcedMeshMotionHandler
    ioData.forced.frequency *= ioData.ref.rv.time;        // Problem in ForcedMeshMotionHandler

    ioData.eqs.gravity_x /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.eqs.gravity_y /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.eqs.gravity_z /= ioData.ref.rv.velocity / ioData.ref.rv.time;

    ioData.bc.hydro.depth /= ioData.ref.length;

    double theta_k = 1.0;
    double theta_w = 10.0;
    if (kenergy0 == pow(10.0, -theta_k) * theta_w / reynolds0) {
      ioData.bc.inlet.kenergy = pow(10.0, -theta_k) * theta_w / ioData.ref.reynolds_mu;
      ioData.bc.inlet.eps     = ioData.eqs.tc.tm.ke.c_mu * ioData.bc.inlet.kenergy * theta_w;
    }

    Pin = ioData.aero.pressure;

    //////////////////////////

    this->postOp->rstVarPostFcn(ioData);
    this->postOp->rstVar(ioData);

    if (ioData.eqs.type == EquationsData::NAVIER_STOKES)
      this->spaceOp->rstVarFet(ioData);

    this->spaceOp->rstFluxFcn(ioData);

    mvp->rstSpaceOp(ioData, this->varFcn, this->spaceOp, false);

    this->rstVarImplicitEmbeddedCoupledTsDesc(ioData);

    this->bcData->rstVar(ioData);

    (this->timeState->getData()).rstVar(ioData);
 
    this->timeState->rstVar(ioData);

    this->output->rstVar(ioData);

    this->restart->rstVar(ioData);
 
    this->data->rstVar(ioData);

    this->refVal->rstVar(ioData);
       
    //OFF this->varFcn->rstVar(ioData);
    
  }

  // initialize boundary condition fluxes
  this->bcData->initialize(ioData, *this->X);

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

template<int dim>
int EmbeddedFluidShapeOptimizationHandler<dim>::fsoHandler(IoData &ioData, DistSVec<double,dim> &U)
{

  // xmach      -  Mach number
  // alpha      -  pitch angle
  // teta       -  yaw angle
  // DFSPAR(1)  -  Mach number differential
  // DFSPAR(2)  -  angle of attack differential
  // DFSPAR(3)  -  yaw angle differential

  double MyLocalTimer = -this->timer->getTime();

  
  double dtLeft = 0.0;
  this->computeTimeStep(1, &dtLeft, U);
  this->computeMeshMetrics();
  this->updateStateVectors(U);
  fsoSetUpLinearSolver(ioData, *this->X, *this->A, U, dFdS);
 
  if(ioData.sa.sensMesh  == SensitivityAnalysis::ON_SENSITIVITYMESH)  fso_on_sensitivityMesh(ioData,  U);
  if(ioData.sa.sensMach  == SensitivityAnalysis::ON_SENSITIVITYMACH)  fso_on_sensitivityMach(ioData,  U);
  if(ioData.sa.sensAlpha == SensitivityAnalysis::ON_SENSITIVITYALPHA) fso_on_sensitivityAlpha(ioData, U); 
  if(ioData.sa.sensBeta  == SensitivityAnalysis::ON_SENSITIVITYBETA)  fso_on_sensitivityBeta(ioData,  U); 

  bool lastIt = true;

  //this->output->closeAsciiFiles(); ??
  
  MyLocalTimer += this->timer->getTime();

  if (this->com->cpuNum() == 0)  {
    std::cout << "\n *** FluidShapeOptimizationHandler::fsoHandler >> Exit";
    std::cout << " (" << MyLocalTimer << " s)";
    std::cout << "\n\n";
  }
  
  return -1;

}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidShapeOptimizationHandler<dim>::fso_on_sensitivityMesh
(
IoData &ioData, 
DistSVec<double,dim> &U
){

  double tag = 0.0;

  step = 0;

  DFSPAR[0] = 0.0;
  DFSPAR[1] = 0.0;
  DFSPAR[2] = 0.0;

  actvar = 1;
  
  while (true) {

    // Reading derivative of the overall deformation
    bool readOK = getdXdSb(step);
    if(!readOK) break;

    this->com->fprintf(stderr, "\n ***** Shape variable #%d\n", step);

    fsoComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U);

    fsoComputeSensitivities(ioData, "Derivatives with respect to the mesh position:", ioData.sa.sensoutput, *this->X, U);

    step = step + 1;
  }

  fsoPrintTextOnScreen("\n ***** Derivatives of mesh position were computed! \n");
  
}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidShapeOptimizationHandler<dim>::fso_on_sensitivityMach
(
 IoData &ioData, 
 DistSVec<double,dim> &U
){

  DFSPAR[0] = 1.0;
  DFSPAR[1] = 0.0;
  DFSPAR[2] = 0.0;
  actvar = 2;

  fsoComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U);
  
  fsoComputeSensitivities(ioData, "Derivatives with respect to the Mach number:", ioData.sa.sensoutput, *this->X, U);

  fsoPrintTextOnScreen("\n ***** Derivatives with respect to the Mach number were computed! \n");

  step = step + 1;

}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidShapeOptimizationHandler<dim>::fso_on_sensitivityBeta
(
 IoData &ioData, 
 DistSVec<double,dim> &U
){

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

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidShapeOptimizationHandler<dim>::fso_on_sensitivityAlpha
(
 IoData &ioData, 
 DistSVec<double,dim> &U
)
{

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

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidShapeOptimizationHandler<dim>::fsoComputeDerivativesOfFluxAndSolution
(
 IoData &ioData, 
 DistSVec<double,3> &X, 
 DistVec<double> &A, 
 DistSVec<double,dim> &U, 
 bool isFSI
){

  if ( ioData.sa.scFlag == SensitivityAnalysis::ANALYTICAL ) {

    fsoAnalytical(ioData, X, A, U, dFdS);
       
    /*
    //debug----------------------------------------
    dFdS_debug = 0.0;
    fsoSemiAnalytical(ioData, X, A, U, dFdS_debug);
    difference = dFdS - dFdS_debug;
    this->com->fprintf(stderr, "!!! dFdS and dFdSref do not match. The relative difference norm is %e\n", difference.norm()/dFdS.norm());
    this->com->fprintf(stderr, "!!! dFdSref norm is %e\n", dFdS.norm()); 
    //exit(-1);
    //debug----------------------------------------
    */

  } else {
    fsoSemiAnalytical(ioData, X, A, U, dFdS);
  }

  fsoLinearSolver(ioData, dFdS, dUdS, isFSI);

}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidShapeOptimizationHandler<dim>::fsoAnalytical
(
  IoData &ioData, 
  DistSVec<double,3> &X,
  DistVec<double> &A,
  DistSVec<double,dim> &U,
  DistSVec<double,dim> &dFdS
 ){

  // Computing the derivatives of the boundary fluxes
  dXdS = 0.0;
  this->bcData->initializeSA(ioData, X, dXdS, DFSPAR[0], DFSPAR[1], DFSPAR[2]);
 
  // Computing the partial derivative of the flux with respect to the variables
  this->spaceOp->computeDerivativeOfResidual(X, A, U, this->distLSS, 
					     this->linRecAtInterface, this->viscSecOrder, 
					     this->nodeTag, this->riemann, this->riemannNormal,
					     this->ghostPoints, DFSPAR[0], 
					     Flux, dFdS, this->timeState);

  this->spaceOp->applyBCsToDerivativeOfResidual(U, dFdS);

  // dF/dS_rad ---> dF/dS_deg 
  if(DFSPAR[1] || DFSPAR[2]) dFdS *= (acos(-1.0)/180.0);
  
}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidShapeOptimizationHandler<dim>::fsoSemiAnalytical
(
  IoData &ioData, 
  DistSVec<double,3> &X,
  DistVec<double> &A,
  DistSVec<double,dim> &U,
  DistSVec<double,dim> &dF
)
{

  if (Fp == 0) {
    fprintf(stderr, "*** Error: Variable Fp does not exist!\n");
    exit(1);
  }
  if (Fm == 0) {
    fprintf(stderr, "*** Error: Variable Fm does not exist!\n");
    exit(1);
  }
  if (X_ == 0) {
    fprintf(stderr, "*** Error: Variable X_ does not exist!\n");
    exit(1);
  }
  if (A_ == 0) {
    fprintf(stderr, "*** Error: Variable A_ does not exist!\n");
    exit(1);
  }

  double dtLeft;
  double xmachc;
  double alpradc;
  double tetac;

  double eps = ioData.sa.eps;

  xmachc  = xmach;
  alpradc = alprad;
  tetac   = teta;

  *Fp = 0.0;
  *Fm = 0.0;
   dF = 0.0;
  *X_ = 0.0;
  *A_ = 0.0;

  DistVec<int> pb_dummy(this->domain->getNodeDistInfo());
  *X_ = X;

  //------
  // Plus
  //---------------------------------------------------------------------
  if(ioData.sa.sensMesh  == SensitivityAnalysis::ON_SENSITIVITYMESH){
    this->distLSS->updateXb(eps);
    this->distLSS->initialize(this->domain, X, this->geoState->getXn(), ioData, &pb_dummy, &this->nodeTag);
  }

  xmach  = xmachc  + eps*DFSPAR[0];
  alprad = alpradc + eps*DFSPAR[1];
  teta   = tetac   + eps*DFSPAR[2];

  fsoRestartBcFluxs(ioData);

  this->geoState->reset(*X_);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), *X_, *A_);
  this->bcData->update(*X_);
  
  dtLeft = 0.0;
  this->computeTimeStep(1, &dtLeft, U);
  
  this->spaceOp->computeResidual(*X_, *A_, U, 
											*this->Wstarij, *this->Wstarji, *this->Wextij, this->distLSS,
                                 this->linRecAtInterface, this->viscSecOrder, 
				 this->nodeTag, 
				 *Fp, 
				 this->riemann, this->riemannNormal, 1, this->ghostPoints);
  
  this->spaceOp->applyBCsToResidual(U, *Fp, this->distLSS);
  //---------------------------------------------------------------------
  
  //-------
  // Minus
  //-------------------------------
  if(ioData.sa.sensMesh  == SensitivityAnalysis::ON_SENSITIVITYMESH){
    this->distLSS->updateXb(-1.0*eps);
    this->distLSS->initialize(this->domain, X, this->geoState->getXn(), ioData, &pb_dummy, &this->nodeTag);
  }

  xmach  = xmachc  - eps*DFSPAR[0];
  alprad = alpradc - eps*DFSPAR[1];
  teta   = tetac   - eps*DFSPAR[2];

  fsoRestartBcFluxs(ioData);

  this->geoState->reset(*X_);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), *X_, *A_);
  this->bcData->update(*X_);

  dtLeft = 0.0;
  this->computeTimeStep(1, &dtLeft, U);
  
  this->spaceOp->computeResidual(*X_, *A_, U, 
											*this->Wstarij, *this->Wstarji, *this->Wextij, this->distLSS,
                                 this->linRecAtInterface, this->viscSecOrder, 
				 this->nodeTag, 
				 *Fm, 
				 this->riemann, this->riemannNormal, 1, this->ghostPoints);
  
  this->spaceOp->applyBCsToResidual(U, *Fm, this->distLSS);
  //---------------------------------------------------------------------

  dF = (1.0/(2.0*eps))*((*Fp) - (*Fm));

  // dF/dS_rad ---> dF/dS_deg 
  if(DFSPAR[1] || DFSPAR[2]) dF *= (acos(-1.0)/180.0);
  
  //-----------------------
  // Reset the steady state
  //-------------------------------
  xmach  = xmachc;
  alprad = alpradc;
  teta   = tetac;

  this->distLSS->updateXb(0.0);
  this->distLSS->initialize(this->domain, X, this->geoState->getXn(), ioData, &pb_dummy, &this->nodeTag);

  fsoRestartBcFluxs(ioData);
  
  this->geoState->reset(X);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), X, A);
  this->bcData->update(X);

  dtLeft = 0.0;
  this->computeTimeStep(1, &dtLeft, U);

  //// ?????
  if ( ioData.sa.scFlag == SensitivityAnalysis::SEMIANALYTICAL ) {
    this->geoState->updateConfigSA();
    this->bcData->initializeSA(ioData, X, dXdS, DFSPAR[0], DFSPAR[1], DFSPAR[2]);
  }

}
//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidShapeOptimizationHandler<dim>::fsoLinearSolver
(
  IoData &ioData, 
  DistSVec<double,dim> &dFdS, 
  DistSVec<double,dim> &dUdS,
  bool isFSI
)
{

  int numberIteration;
  bool istop = false;
  int iter = 0;

  dFdS *= (-1.0);
  
  this->embeddeddQ = 0.0;

  this->embeddedB.ghost() =  0.0;
  this->embeddedB.real()  = dFdS;
  
  ksp->setup(0, 1, this->embeddedB);

  while ((istop == false) && (iter < 100)) {

    numberIteration = ksp->solve(this->embeddedB, this->embeddeddQ);

    if ((!ioData.sa.excsol) || (numberIteration < ioData.sa.ksp.maxIts))
      istop = true; 

    iter += 1;

  }

  dUdS = this->embeddeddQ.real();
  //this->embeddedU.ghost() += this->embeddeddQ.ghost();

  dFdS *= (-1.0);

}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidShapeOptimizationHandler<dim>::fsoComputeSensitivities
(
 IoData &ioData, 
 const char *mesage, 
 const char *fileName, 
 DistSVec<double,3> &X, 
 DistSVec<double,dim> &U
){

  // Computing efforts (F: force, M: moment, L: LiftAndDrag)
  Vec3D F = 0.0, M = 0.0, L = 0.0;
  fsoGetEfforts(ioData, X, U, F, M, L);

  // Computing derivative of the efforts
  Vec3D dFds = 0.0, dMds = 0.0, dLds = 0.0;
  
  if ( ioData.sa.scFlag == SensitivityAnalysis::FINITEDIFFERENCE ){
    
    fsoGetDerivativeOfEffortsFiniteDifference(ioData, X, U, dUdS, dFds, dMds, dLds);

  } else {

    fsoGetDerivativeOfEffortsAnalytical(ioData, X, U, dUdS, dFds, dMds, dLds);

    /*
    //debug-----------------------------------------------------------------------------------
    Vec3D dFds_fd = 0.0, dMds_fd = 0.0, dLs_fd = 0.0;
    fsoGetDerivativeOfEffortsFiniteDifference(ioData, X, U, dUdS, dFds_fd, dMds_fd, dLs_fd);
    //
    this->com->fprintf(stdout, "F     : %+15.5e %+15.5e %+15.5e\n", F[0], F[1], F[2]);
    this->com->fprintf(stdout, "df    : %+15.5e %+15.5e %+15.5e\n", dFds[0], dFds[1], dFds[2]);
    this->com->fprintf(stdout, "df_fd : %+15.5e %+15.5e %+15.5e\n", dFds_fd[0], dFds_fd[1], dFds_fd[2]);
    this->com->fprintf(stdout, "err_F : %+15.5e %+15.5e %+15.5e\n\n", abs(dFds[0]-dFds_fd[0]), 
                                                	              abs(dFds[1]-dFds_fd[1]), 
                                                                      abs(dFds[2]-dFds_fd[2]));
    //
    this->com->fprintf(stdout, "dm    : %+15.5e %+15.5e %+15.5e\n", dMds[0], dMds[1], dMds[2]);
    this->com->fprintf(stdout, "dm_fd : %+15.5e %+15.5e %+15.5e\n", dMds_fd[0], dMds_fd[1], dMds_fd[2]);
    this->com->fprintf(stdout, "err_M : %+15.5e %+15.5e %+15.5e\n\n", abs(dMds[0]-dMds_fd[0]), 
                                 	                              abs(dMds[1]-dMds_fd[1]), 
                                                                      abs(dMds[2]-dMds_fd[2]));

    // 
    this->com->fprintf(stdout, "L     : %+15.5e %+15.5e %+15.5e\n",      L[0],      L[1],      L[2]);
    this->com->fprintf(stdout, "dL    : %+15.5e %+15.5e %+15.5e\n",   dLds[0],   dLds[1],   dLds[2]);
    this->com->fprintf(stdout, "dL_fd : %+15.5e %+15.5e %+15.5e\n", dLs_fd[0], dLs_fd[1], dLs_fd[2]);
    this->com->fprintf(stdout, "err_L : %+15.5e %+15.5e %+15.5e\n\n", abs(dLds[0]-dLs_fd[0]), 
                                                	              abs(dLds[1]-dLs_fd[1]), 
                                                                      abs(dLds[2]-dLs_fd[2]));
    //debug-----------------------------------------------------------------------------------
    */  
  }

  if ((!ioData.sa.angleRad) && (DFSPAR[1] || DFSPAR[2])) {
    dFds *= acos(-1.0) / 180.0;
    dMds *= acos(-1.0) / 180.0;
  }  
  
  if (this->com->cpuNum() == 0) {
    outFile = fopen(fileName,"a+");
    if (outFile) {
      this->com->fprintf(outFile,mesage);
      this->com->fprintf(outFile,"\n");
      this->com->fprintf(outFile,"Fx= %24.13e \n",F[0]);
      this->com->fprintf(outFile,"Fy= %24.13e \n",F[1]);
      this->com->fprintf(outFile,"Fz= %24.13e \n",F[2]);
      this->com->fprintf(outFile,"dFx/ds= %24.13e \n",dFds[0]);
      this->com->fprintf(outFile,"dFy/ds= %24.13e \n",dFds[1]);
      this->com->fprintf(outFile,"dFz/ds= %24.13e \n",dFds[2]);
      this->com->fprintf(outFile,"Mx= %24.13e \n",M[0]);
      this->com->fprintf(outFile,"My= %24.13e \n",M[1]);
      this->com->fprintf(outFile,"Mz= %24.13e \n",M[2]);
      this->com->fprintf(outFile,"dMx/ds= %24.13e \n",dMds[0]);
      this->com->fprintf(outFile,"dMy/ds= %24.13e \n",dMds[1]);
      this->com->fprintf(outFile,"dMz/ds= %24.13e \n",dMds[2]);
      this->com->fprintf(outFile,"\n");
      fclose(outFile);
    }
  }

  double  sboom = 0.0;
  double dSboom = 0.0;


  // This function is simply writing to the disk.
  this->output->writeDerivativeOfForcesToDisk(step, actvar, F, dFds, M, dMds, sboom, dSboom);
  this->output->writeDerivativeOfLiftDragToDisk(step, actvar, L, dLds); 
 
  // This function is writing to the disk quantities of interest in binary files.
  // The possible quantities of interest include
  // - dUdS
  // - Derivative of Scalar Quantities: Density, Mach, Pressure, Temperature, TotPressure,
  //   NutTurb, EddyViscosity, VelocityScalar
  // - Derivative of Vector Quantities: VelocityVector, Displacement
  this->output->writeBinaryDerivativeOfVectorsToDisk(step+1, actvar, DFSPAR, *this->X, dXdS, U, dUdS, this->timeState);


}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidShapeOptimizationHandler<dim>::fsoGetEfforts
(IoData &ioData, 
 DistSVec<double,3> &X, 
 DistSVec<double,dim> &U, 
 Vec3D &F, Vec3D &M, Vec3D &L
){

  int nSurfs = this->postOp->getNumSurf();

  Vec3D x0, force, moment;

  Vec3D *Fi = new Vec3D[nSurfs];
  Vec3D *Mi = new Vec3D[nSurfs];
  Vec3D *Fv = new Vec3D[nSurfs];
  Vec3D *Mv = new Vec3D[nSurfs];

  x0[0] = ioData.output.transient.x0;
  x0[1] = ioData.output.transient.y0;
  x0[2] = ioData.output.transient.z0;

  this->postOp->computeForceAndMoment(x0, X, U, &this->nodeTag, Fi, Mi, Fv, Mv);

  F = 0.0; F = Fi[0] + Fv[0];
  M = 0.0; M = Mi[0] + Mv[0];

  // DIMENSIONAL
  F *= this->refVal->force;
  M *= this->refVal->energy;
  //

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
void EmbeddedFluidShapeOptimizationHandler<dim>::fsoGetDerivativeOfEffortsAnalytical
(IoData &ioData,  DistSVec<double,3> &X,
 DistSVec<double,dim> &U, DistSVec<double,dim> &dU,
 Vec3D &dForces, Vec3D &dMoments, Vec3D &dL
){

  double gamma     = ioData.eqs.fluidModel.gasModel.specificHeatRatio;

  double velocity  = ioData.ref.mach * sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
  double dVelocity = sqrt(gamma * ioData.ref.pressure / ioData.ref.density)*DFSPAR[0];

  double dForce    = 2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*velocity*dVelocity;
  double dEnergy   = 2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*ioData.ref.length*velocity*dVelocity;

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

  this->postOp->computeForceAndMoment(x0, X, U, &this->nodeTag, Fi, Mi, Fv, Mv);

  F = 0.0;  M = 0.0;
  F = Fi[0] + Fv[0];
  M = Mi[0] + Mv[0];
 
  this->postOp->computeDerivativeOfForceAndMoment(x0, X, U, dU, &this->nodeTag, DFSPAR, dFi, dMi, dFv, dMv);

  dF = 0.0; dM = 0.0;
  dF = dFi[0] + dFv[0];
  dM = dMi[0] + dMv[0];

  // DIMENSIONAL
  dF *= this->refVal->force;
  dM *= this->refVal->energy;

  F *= dForce;
  M *= dEnergy;

  dForces  = dF + F;
  dMoments = dM + M;
  //

  dL = 0.0;
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

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidShapeOptimizationHandler<dim>::fsoGetDerivativeOfEffortsFiniteDifference
(IoData &ioData,  DistSVec<double,3> &X,
 DistSVec<double,dim> &U, DistSVec<double,dim> &dU,
 Vec3D &dForces, Vec3D &dMoments,  Vec3D &dL){

  //Remark: Error mesage for pointers
  if (Up == 0) {
    fprintf(stderr, "*** Error: Variable Up does not exist!\n");
    exit(1);
  }
  if (Um == 0) {
    fprintf(stderr, "*** Error: Variable Um does not exist!\n");
    exit(1);
  }
  if (X_ == 0) {
    fprintf(stderr, "*** Error: Variable Xp does not exist!\n");
    exit(1);
  }
  if (A_ == 0) {
    fprintf(stderr, "*** Error: Variable Ap does not exist!\n");
    exit(1);
  }

  double eps = ioData.sa.eps;

  double gamma     = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
  double velocity  = ioData.ref.mach * sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
  double dVelocity = sqrt(gamma * ioData.ref.pressure / ioData.ref.density)*DFSPAR[0];

  double dForce    = 2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*velocity*dVelocity;
  double dEnergy   = 2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*ioData.ref.length*velocity*dVelocity;

  Vec3D x0, F, Fplus, Fminus, dF, M, Mplus, Mminus, dM;

  double xmachc;
  double alpradc;
  double tetac;

  int nSurfs = this->postOp->getNumSurf();

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
  *X_ = 0.0;
  *A_ = 0.0;

  DistVec<int> pb_dummy(this->domain->getNodeDistInfo());

  this->postOp->computeForceAndMoment(x0, X, U, &this->nodeTag, Fi, Mi, Fv, Mv); 

  F = 0.0;  M = 0.0;
  F = Fi[0] + Fv[0];
  M = Mi[0] + Mv[0];

  xmachc  = xmach; 
  alpradc = alprad;
  tetac   = teta;

  *X_ = X;

  //------
  // Plus
  //---------------------------------------------------------------------
  *Up = U + eps*dU;
 
  if(ioData.sa.sensMesh == SensitivityAnalysis::ON_SENSITIVITYMESH) {
    this->distLSS->updateXb(eps);
    this->distLSS->initialize(this->domain, X, this->geoState->getXn(), ioData, &pb_dummy, &this->nodeTag);
  }

  xmach  = xmachc  + eps*DFSPAR[0];
  alprad = alpradc + eps*DFSPAR[1];
  teta   = tetac   + eps*DFSPAR[2];

  fsoRestartBcFluxs(ioData);

  this->geoState->reset(*X_);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), *X_, *A_);
  this->bcData->update(*X_);

  this->spaceOp->computeResidual(*X_, *A_, *Up, 
											*this->Wstarij, *this->Wstarji, *this->Wextij, this->distLSS,
                                 this->linRecAtInterface, this->viscSecOrder, 
				 this->nodeTag, 
				 *Fp, //dummy 
				 this->riemann, this->riemannNormal, 1, this->ghostPoints);

  this->postOp->computeForceAndMoment(x0, *X_, *Up, &this->nodeTag, Fip, Mip, Fvp, Mvp);

  Fplus = 0.0; Mplus = 0.0;
  Fplus = Fip[0] + Fvp[0];
  Mplus = Mip[0] + Mvp[0];

  this->com->fprintf(stdout, "F_plus: %+15.5e %+15.5e %+15.5e\n", Fplus[0]*this->refVal->force, 
		                                                  Fplus[1]*this->refVal->force, 
		                                                  Fplus[2]*this->refVal->force);

  //-------
  // Minus
  //---------------------------------------------------------------------
  *Um = U - eps*dU;

  if(ioData.sa.sensMesh  == SensitivityAnalysis::ON_SENSITIVITYMESH) {
    this->distLSS->updateXb(-1.0*eps);
    this->distLSS->initialize(this->domain, X, this->geoState->getXn(), ioData, &pb_dummy, &this->nodeTag);
  }

  xmach  = xmachc  - eps*DFSPAR[0];
  alprad = alpradc - eps*DFSPAR[1];
  teta   = tetac   - eps*DFSPAR[2];

  fsoRestartBcFluxs(ioData);

  this->geoState->reset(*X_);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), *X_, *A_);
  this->bcData->update(*X_);

  this->spaceOp->computeResidual(*X_, *A_, *Um, 
											*this->Wstarij, *this->Wstarji, *this->Wextij, this->distLSS,
                                 this->linRecAtInterface, this->viscSecOrder, 
				 this->nodeTag, 
				 *Fm, //dummy 
				 this->riemann, this->riemannNormal, 1, this->ghostPoints);

  this->postOp->computeForceAndMoment(x0, *X_, *Um, &this->nodeTag, Fim, Mim, Fvm, Mvm);

  Fminus = 0.0;  Mminus = 0.0;
  Fminus = Fim[0] + Fvm[0];
  Mminus = Mim[0] + Mvm[0];

  this->com->fprintf(stdout, "F_minu: %+15.5e %+15.5e %+15.5e\n", Fminus[0]*this->refVal->force, 
		                                                  Fminus[1]*this->refVal->force, 
		                                                  Fminus[2]*this->refVal->force);

  // Compute the derivatives
  dF = 1.0/(2.0*eps)*(Fplus - Fminus);
  dM = 1.0/(2.0*eps)*(Mplus - Mminus);
  //---------------------------------------------------------------------

  // Reset values to the steady state
  xmach  = xmachc;
  alprad = alpradc;
  teta   = tetac;

  if(ioData.sa.sensMesh  == SensitivityAnalysis::ON_SENSITIVITYMESH){
    this->distLSS->updateXb(0.0);
    this->distLSS->initialize(this->domain, X, this->geoState->getXn(), ioData, &pb_dummy, &this->nodeTag);
  }

  fsoRestartBcFluxs(ioData);

  this->geoState->reset(X);
  this->geoState->compute(this->timeState->getData(), this->bcData->getVelocityVector(), X, *this->A);
  this->bcData->update(X);

  //Dimensional
  dF *= this->refVal->force;
  dM *= this->refVal->energy;

   F *= dForce;
   M *= dEnergy;

  dForces  = dF + F;
  dMoments = dM + M;

  ///
 
  dL = 0.0;
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

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidShapeOptimizationHandler<dim>::fsoPrintTextOnScreen(const char *Text){
   this->com->fprintf(stderr, Text);
}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedFluidShapeOptimizationHandler<dim>::fsaPrintTextOnScreen(const char *Text){
   this->com->fprintf(stderr, Text);
}

//------------------------------------------------------------------------------

template<int dim>
bool EmbeddedFluidShapeOptimizationHandler<dim>::getdXdSb(int istep){

  bool found = false;

  FILE *dFile;
  dFile = fopen(dXdSb_file, "r");

  if(dFile == NULL) {
    this->com->fprintf(stderr, "Embedded surface sensitivity file  (%s) doesn't exist.\n", dXdSb_file);
    exit(-1);
  }
  
  const int MAXrLINE = 500;
  char line[MAXrLINE];

  fgets(line, MAXrLINE, dFile);

  int NumNodes;
  fgets(line, MAXrLINE, dFile);
  sscanf(line, "%d", &NumNodes);

  if(NumNodes < 1) {
    this->com->fprintf(stderr, "Error: total num nodes in Embedded surface sensitivity file is %i \n", NumNodes);
    exit(-1);
  }

  double* dxdSb = new double[NumNodes];
  double* dydSb = new double[NumNodes];
  double* dzdSb = new double[NumNodes];

  int key, id;

  while (fgets(line, MAXrLINE, dFile) != 0) {

    sscanf(line, "%d", &key);
    if(key == istep){

      found = true;
      for(int i=0; i<NumNodes; ++i){
        fgets(line, MAXrLINE, dFile);
        //sscanf(line, "%d %lf %lf %lf", &id, &dxdSb[i], &dydSb[i], &dzdSb[i]);
	sscanf(line, "%lf %lf %lf", &dxdSb[i], &dydSb[i], &dzdSb[i]);
      }
      break;

    } else {

      found = false;
      for(int i=0; i<NumNodes; ++i)
	fgets(line, MAXrLINE, dFile);

    }

  }
  
  fclose(dFile);

  if(found)
    this->distLSS->setdXdSb(NumNodes, dxdSb, dydSb, dzdSb);
 
  delete [] dxdSb;
  delete [] dydSb;
  delete [] dzdSb;

  return found;

}

//------------------------------------------------------------------------------
