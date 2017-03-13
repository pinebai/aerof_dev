/*
 * FluidSegShapeOptimizationHandler.C
 *
 *  Created on: Dec 7, 2016
 *      Author: lscheuch
 */

#include <FluidSegShapeOptimizationHandler.h>

#include <IoData.h>
#include <Domain.h>
#include <GeoSource.h>
#include <DistVector.h>
#include <MeshMotionSolver.h>
//#include <StructExc.h>
#include <MatVecProd.h>
#include <KspPrec.h>
#include <KspSolver.h>
#include <MemoryPool.h>

#include <cmath>

#include <iostream>
#include <string>

#include "Dev/devtools.h"

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
FluidSegShapeOptimizationHandler<dim,neq1,neq2>::FluidSegShapeOptimizationHandler
(
  IoData &ioData,
  GeoSource &geoSource,
  Domain *dom//,
) :
ImplicitSegTsDesc<dim,neq1,neq2>(ioData, geoSource, dom),
domain(dom),
dXb(dom->getNodeDistInfo()),
dXdS(dom->getNodeDistInfo()),
dXdSb(dom->getNodeDistInfo()),
Xc(dom->getNodeDistInfo()),
dAdS(dom->getNodeDistInfo()),
dFdS(dom->getNodeDistInfo()),
dFdSref(dom->getNodeDistInfo()),
dUdS(dom->getNodeDistInfo()),
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
      mvp1 = new MatVecProdH2<dim,MatScalar,neq1>(ioData, this->varFcn, this->timeState, this->spaceOp1, domain, this->geoState);
      mvp2 = new MatVecProdH2<dim,MatScalar,neq2>(ioData, this->varFcn, this->timeState, this->spaceOp2, domain, this->geoState);
    }
    else
    {
      mvp1 = new MatVecProdFD<dim,neq1>(ioData.ts.implicit, this->timeState, this->geoState, this->spaceOp1, domain, ioData);
      mvp2 = new MatVecProdFD<dim,neq2>(ioData.ts.implicit, this->timeState, this->geoState, this->spaceOp2, domain, ioData);
    }
  }
  else
  {
    if (ioData.sa.mvp == SensitivityAnalysis::H2)
    {
      mvp1 = new MatVecProdH2<dim,MatScalar,neq1>(ioData, this->varFcn, 0, this->spaceOp1, domain, this->geoState);
      mvp2 = new MatVecProdH2<dim,MatScalar,neq2>(ioData, this->varFcn, 0, this->spaceOp2, domain, this->geoState);
    }
    else
    {
      mvp1 = new MatVecProdFD<dim,neq1>(ioData.ts.implicit,  0, this->geoState, this->spaceOp1, domain, ioData);
      mvp2 = new MatVecProdFD<dim,neq2>(ioData.ts.implicit,  0, this->geoState, this->spaceOp2, domain, ioData);
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

  if (ioData.sa.sensMesh == SensitivityAnalysis::ON_SENSITIVITYMESH || ioData.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_) {
    mms = new TetMeshMotionSolver(ioData.dmesh, geoSource.getMatchNodes(),domain,0);
  } else mms = 0;

  length = ioData.output.transient.length;
  surface = ioData.output.transient.surface;

  numLocSub = domain->getNumLocSub();

  dXdS=0.0;
  dFdS=0.0;
  dFdSref=0.0;
  dUdS=0.0;
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

template<int dim, int neq1, int neq2>
FluidSegShapeOptimizationHandler<dim,neq1,neq2>::~FluidSegShapeOptimizationHandler()
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

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoRestartBcFluxs(IoData &ioData)
{
  //TODO HACK
  //ioData.ref.mach=ioData.sa.machref;

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

        std::cout<<"\033[96m"<<__FILE__<<__LINE__<<"\033[00m"<<std::endl;//TODO delete line
        this->rstVarImplicitSegTsDesc(ioData);//TODO HACK
        std::cout<<"\033[96m"<<__FILE__<<__LINE__<<"\033[00m"<<std::endl;//TODO delete line
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

    //ioData.ref.mach = ioData.bc.inlet.mach;//TODO HACK
    ioData.ref.density = ioData.bc.inlet.density;
    ioData.ref.pressure = ioData.bc.inlet.pressure;
    std::cout<<"        ioData.sa"<<ioData.sa.machref<<std::endl;//TODO delelte line
    std::cout<<"        ioData.bc.inlet.mach"<<ioData.bc.inlet.mach<<std::endl;//TODO delelte line
    std::cout<<"        ioData.ref.mach"<<ioData.ref.mach<<std::endl;//TODO delelte line
    std::cout<<"        gamma"<<gamma<<std::endl;//TODO delelte line
    std::cout<<"        ioData.ref.pressure"<<ioData.ref.pressure<<std::endl;//TODO delelte line
    std::cout<<"        Pstiff"<<Pstiff<<std::endl;//TODO delelte line
    std::cout<<"        ioData.ref.density"<<ioData.ref.density<<std::endl;//TODO delelte line
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

    //std::cout<<"XXXXXXX check 7"<<std::endl;
    ioData.ref.rv.mode = RefVal::DIMENSIONAL;
    //std::cout<<"XXXXXXX check 7.1"<<std::endl;
    ioData.ref.rv.density = ioData.ref.density;
    //std::cout<<"XXXXXXX check 7.2"<<std::endl;
    ioData.ref.rv.velocity = velocity;
    //std::cout<<"XXXXXXX Density: "<<ioData.ref.density<<std::endl;
    //std::cout<<"XXXXXXX Velocity: "<<velocity<<std::endl;
    ioData.ref.rv.pressure = ioData.ref.density * velocity*velocity;
    //std::cout<<"XXXXXXX check 7.3"<<std::endl;
    ioData.ref.rv.temperature = gamma*(gamma - 1.0) * ioData.ref.mach*ioData.ref.mach * (ioData.ref.pressure+Pstiff)/(R*ioData.ref.density);
    ioData.ref.rv.viscosity_mu = viscosity;
    //std::cout<<"XXXXXXX check 7.4"<<std::endl;
    ioData.ref.rv.nutilde = viscosity / ioData.ref.density;
    //std::cout<<"XXXXXXX check 7.5"<<std::endl;
    ioData.ref.rv.kenergy = velocity*velocity;
    //std::cout<<"XXXXXXX check 7.6"<<std::endl;
    ioData.ref.rv.epsilon = velocity*velocity*velocity / ioData.ref.length;
    //std::cout<<"XXXXXXX check 7.7"<<std::endl;
    ioData.ref.rv.time = ioData.ref.length / velocity;    // Problem in RigidMeshMotionHandler
    //std::cout<<"XXXXXXX check 7.8"<<std::endl;
    ioData.ref.rv.force = ioData.ref.density * velocity*velocity * ioData.ref.length*ioData.ref.length;
    //std::cout<<"XXXXXXX check 7.9"<<std::endl;
    ioData.ref.rv.energy = ioData.ref.density * velocity*velocity * ioData.ref.length*ioData.ref.length*ioData.ref.length;
    //std::cout<<"XXXXXXX check 7.10"<<std::endl;
    ioData.ref.rv.power = ioData.ref.density * velocity*velocity*velocity * ioData.ref.length*ioData.ref.length;
    ioData.ref.rv.tvelocity = velocity / ioData.aero.displacementScaling;
    //std::cout<<"XXXXXXX check 7.11"<<std::endl;
    ioData.ref.rv.tforce = ioData.ref.rv.force / ioData.aero.forceScaling;
    //std::cout<<"XXXXXXX check 7.12"<<std::endl;
    ioData.ref.rv.tpower = ioData.ref.rv.power / ioData.aero.powerScaling;

    ioData.ref.rv.dvelocitydMach = dvelocitydMach;
    //std::cout<<"XXXXXXX check 7.13"<<std::endl;
    ioData.ref.rv.dtimedMach = - ioData.ref.length / (velocity * velocity) * dvelocitydMach;

    ioData.eqs.fluidModel.pmin /= ioData.ref.rv.pressure;

    //
    // Step 2.3: Reset values of bc.inlet
    //

    //std::cout<<"XXXXXXX check 8"<<std::endl;
    ioData.bc.inlet.density /= ioData.ref.rv.density;
    ioData.bc.inlet.pressure /= ioData.ref.rv.pressure;
    ioData.bc.inlet.temperature /= ioData.ref.rv.temperature;
    ioData.bc.inlet.nutilde /= ioData.ref.rv.nutilde;
    ioData.bc.inlet.kenergy /= ioData.ref.rv.kenergy;
    ioData.bc.inlet.eps /= ioData.ref.rv.epsilon;

    //
    // Step 2.4: Reset values of bc.outlet
    //

    //std::cout<<"XXXXXXX check 9"<<std::endl;
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

    //std::cout<<"XXXXXXX check 10"<<std::endl;
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
    //std::cout<<"XXXXXXX check 11"<<std::endl;
    ioData.aero.pressure /= ioData.ref.rv.pressure;
    ioData.forced.timestep /= ioData.ref.rv.time;         // Problem in ForcedMeshMotionHandler
    ioData.forced.frequency *= ioData.ref.rv.time;        // Problem in ForcedMeshMotionHandler

    ioData.eqs.gravity_x /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.eqs.gravity_y /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.eqs.gravity_z /= ioData.ref.rv.velocity / ioData.ref.rv.time;
    ioData.bc.hydro.depth /= ioData.ref.length;

    //std::cout<<"XXXXXXX check 12"<<std::endl;
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

    std::cout<<"\033[96m"<<__FILE__<<__LINE__<<"\033[00m"<<std::endl;//TODO delete line
    this->rstVarImplicitSegTsDesc(ioData);//TODO HACK
    std::cout<<"\033[96m"<<__FILE__<<__LINE__<<"\033[00m"<<std::endl;//TODO delete line

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






template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsaRestartBcFluxs(IoData &ioData)
{
  Dev::Error(this->com,"Not correctly implemented yet",true);
//  //TODO HACK
//  ioData.ref.mach=ioData.sa.machref;
//
//  double gamma = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
//  double R = ioData.eqs.fluidModel.gasModel.idealGasConstant;
//  double Pstiff = ioData.eqs.fluidModel.gasModel.pressureConstant;
//
//// Remark: For internal flows the SA using inlet or outlet Mach, Alpha and Beta should be specified in the input file
//
//  ioData.bc.inlet.mach = xmach;//TODO this is faulty
//  ioData.bc.outlet.mach = xmach;
//
//  ioData.bc.inlet.alpha = alprad;
//  ioData.bc.outlet.alpha = alprad;
//
//  ioData.bc.inlet.beta = teta;
//  ioData.bc.outlet.beta = teta;
//
////  this->com->fprintf(stderr, "\n\n FluidSensitivityAnalysis values: \n\n");
////  this->com->fprintf(stderr, "\n\n Inlet Alpha = %20.17e \n\n",ioData.bc.inlet.alpha);
////  this->com->fprintf(stderr, "\n\n Inlet Beta = %20.17e \n\n",ioData.bc.inlet.beta);
////  this->com->fprintf(stderr, "\n\n Outlet Alpha = %20.17e \n\n",ioData.bc.outlet.alpha);
////  this->com->fprintf(stderr, "\n\n Outlet Beta = %20.17e \n\n",ioData.bc.outlet.beta);
//
//  if (ioData.problem.mode == ProblemData::NON_DIMENSIONAL)
//  {
//    Dev::Error(this->com,"NON-Dimensional currently not supported for SA");
//  }
//  else if (ioData.problem.mode == ProblemData::DIMENSIONAL) {
//
//    //
//    // Step 1: Re-scale all the parameters
//    //
//
//    ioData.eqs.fluidModel.pmin *= ioData.ref.rv.pressure;
//
//    ioData.bc.inlet.density *= ioData.ref.rv.density;
//    ioData.bc.inlet.pressure *= ioData.ref.rv.pressure;
//    ioData.bc.inlet.temperature *= ioData.ref.rv.temperature;
//    ioData.bc.inlet.nutilde *= ioData.ref.rv.nutilde;
//    ioData.bc.inlet.kenergy *= ioData.ref.rv.kenergy;
//    ioData.bc.inlet.eps *= ioData.ref.rv.epsilon;
//    ioData.bc.outlet.density *= ioData.ref.rv.density;
//    ioData.bc.outlet.pressure *= ioData.ref.rv.pressure;
//    ioData.bc.outlet.temperature *= ioData.ref.rv.temperature;
//    ioData.bc.outlet.nutilde *= ioData.ref.rv.nutilde;
//    ioData.bc.outlet.kenergy *= ioData.ref.rv.kenergy;
//    ioData.bc.outlet.eps *= ioData.ref.rv.epsilon;
//
//    ioData.restart.etime *= ioData.ref.rv.time;
//    ioData.restart.dt_nm1 *= ioData.ref.rv.time;
//    ioData.restart.dt_nm2 *= ioData.ref.rv.time;
//    ioData.restart.energy *= ioData.ref.rv.energy;
//    ioData.bc.wall.temperature *= ioData.ref.rv.temperature;
//    ioData.ts.timestep *= ioData.ref.rv.time;
//    ioData.ts.maxTime *= ioData.ref.rv.time;
//    ioData.rmesh.vx *= ioData.ref.rv.velocity;
//    ioData.rmesh.vy *= ioData.ref.rv.velocity;
//    ioData.rmesh.vz *= ioData.ref.rv.velocity;
//    std::cout<<"XXXXXXX check 1"<<std::endl;
//    ioData.rmesh.ax *= ioData.ref.rv.velocity / ioData.ref.rv.time;
//    ioData.rmesh.ay *= ioData.ref.rv.velocity / ioData.ref.rv.time;
//    ioData.rmesh.az *= ioData.ref.rv.velocity / ioData.ref.rv.time;
//    std::cout<<"XXXXXXX check 2"<<std::endl;
//    ioData.rmesh.timestep *= ioData.ref.rv.time;
//
//    for (int j=0; j<ioData.rmesh.num; j++){
//      ioData.rmesh.vpts[j]->time     *= ioData.ref.rv.time;
//      ioData.rmesh.vpts[j]->velocityX *= ioData.ref.rv.velocity;
//      ioData.rmesh.vpts[j]->velocityY *= ioData.ref.rv.velocity;
//      ioData.rmesh.vpts[j]->velocityZ *= ioData.ref.rv.velocity;
//    }
//    std::cout<<"XXXXXXX check 3"<<std::endl;
//    ioData.aero.pressure *= ioData.ref.rv.pressure;
//    ioData.forced.timestep *= ioData.ref.rv.time;
//    ioData.forced.frequency /= ioData.ref.rv.time;
//
//    std::cout<<"XXXXXXX check 4"<<std::endl;
//    ioData.eqs.gravity_x *= ioData.ref.rv.velocity / ioData.ref.rv.time;
//    ioData.eqs.gravity_y *= ioData.ref.rv.velocity / ioData.ref.rv.time;
//    ioData.eqs.gravity_z *= ioData.ref.rv.velocity / ioData.ref.rv.time;
//    ioData.bc.hydro.depth *= ioData.ref.length;
//
//    //
//    // Step 2: Restart all values based on new inlet values
//    // Step 2.1: Reset values in ioData.ref
//    //
//
//    std::cout<<"XXXXXXX check 5"<<std::endl;
//    //ioData.ref.mach = ioData.bc.inlet.mach;//TODO HACK
//    ioData.ref.density = ioData.bc.inlet.density;
//    ioData.ref.pressure = ioData.bc.inlet.pressure;
//    std::cout<<"        ioData.sa"<<ioData.sa.machref<<std::endl;//TODO delelte line
//    std::cout<<"        ioData.bc.inlet.mach"<<ioData.bc.inlet.mach<<std::endl;//TODO delelte line
//    std::cout<<"        ioData.ref.mach"<<ioData.ref.mach<<std::endl;//TODO delelte line
//    std::cout<<"        gamma"<<gamma<<std::endl;//TODO delelte line
//    std::cout<<"        ioData.ref.pressure"<<ioData.ref.pressure<<std::endl;//TODO delelte line
//    std::cout<<"        Pstiff"<<Pstiff<<std::endl;//TODO delelte line
//    std::cout<<"        ioData.ref.density"<<ioData.ref.density<<std::endl;//TODO delelte line
//    double velocity = ioData.sa.machref * sqrt(gamma * (ioData.ref.pressure+Pstiff) / ioData.ref.density);
//
//    ioData.ref.temperature = (ioData.ref.pressure + gamma*Pstiff)/ (ioData.ref.density * R);
//    double viscosity = ioData.eqs.viscosityModel.sutherlandConstant * sqrt(ioData.ref.temperature) /
//      (1.0 + ioData.eqs.viscosityModel.sutherlandReferenceTemperature/ioData.ref.temperature);
//    ioData.ref.reynolds_mu = velocity * ioData.ref.length * ioData.ref.density / viscosity;
//
//    std::cout<<"XXXXXXX check 6"<<std::endl;
//    if (ioData.eqs.type == EquationsData::NAVIER_STOKES)
//      this->com->fprintf(stderr, "\n\n Reynolds = %e \n\n",ioData.ref.reynolds_mu);
//
//    double dvelocitydMach = sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
//    ioData.ref.dRe_mudMach = dvelocitydMach * ioData.ref.length * ioData.ref.density / viscosity;
//
//    //
//    // Step 2.2: Reset values in ioData.ref.rv
//    //
//
//    std::cout<<"XXXXXXX check 7"<<std::endl;
//    ioData.ref.rv.mode = RefVal::DIMENSIONAL;
//    std::cout<<"XXXXXXX check 7.1"<<std::endl;
//    ioData.ref.rv.density = ioData.ref.density;
//    std::cout<<"XXXXXXX check 7.2"<<std::endl;
//    ioData.ref.rv.velocity = velocity;
//    std::cout<<"XXXXXXX Density: "<<ioData.ref.density<<std::endl;
//    std::cout<<"XXXXXXX Velocity: "<<velocity<<std::endl;
//    ioData.ref.rv.pressure = ioData.ref.density * velocity*velocity;
//    std::cout<<"XXXXXXX check 7.3"<<std::endl;
//    ioData.ref.rv.temperature = gamma*(gamma - 1.0) * ioData.sa.machref*ioData.sa.machref * (ioData.ref.pressure+Pstiff)/(R*ioData.ref.density);
//    ioData.ref.rv.viscosity_mu = viscosity;
//    std::cout<<"XXXXXXX check 7.4"<<std::endl;
//    ioData.ref.rv.nutilde = viscosity / ioData.ref.density;
//    std::cout<<"XXXXXXX check 7.5"<<std::endl;
//    ioData.ref.rv.kenergy = velocity*velocity;
//    std::cout<<"XXXXXXX check 7.6"<<std::endl;
//    ioData.ref.rv.epsilon = velocity*velocity*velocity / ioData.ref.length;
//    std::cout<<"XXXXXXX check 7.7"<<std::endl;
//    ioData.ref.rv.time = ioData.ref.length / velocity;    // Problem in RigidMeshMotionHandler
//    std::cout<<"XXXXXXX check 7.8"<<std::endl;
//    ioData.ref.rv.force = ioData.ref.density * velocity*velocity * ioData.ref.length*ioData.ref.length;
//    std::cout<<"XXXXXXX check 7.9"<<std::endl;
//    ioData.ref.rv.energy = ioData.ref.density * velocity*velocity * ioData.ref.length*ioData.ref.length*ioData.ref.length;
//    std::cout<<"XXXXXXX check 7.10"<<std::endl;
//    ioData.ref.rv.power = ioData.ref.density * velocity*velocity*velocity * ioData.ref.length*ioData.ref.length;
//    ioData.ref.rv.tvelocity = velocity / ioData.aero.displacementScaling;
//    std::cout<<"XXXXXXX check 7.11"<<std::endl;
//    ioData.ref.rv.tforce = ioData.ref.rv.force / ioData.aero.forceScaling;
//    std::cout<<"XXXXXXX check 7.12"<<std::endl;
//    ioData.ref.rv.tpower = ioData.ref.rv.power / ioData.aero.powerScaling;
//
//    ioData.ref.rv.dvelocitydMach = dvelocitydMach;
//    std::cout<<"XXXXXXX check 7.13"<<std::endl;
//    ioData.ref.rv.dtimedMach = - ioData.ref.length / (velocity * velocity) * dvelocitydMach;
//
//    ioData.eqs.fluidModel.pmin /= ioData.ref.rv.pressure;
//
//    //
//    // Step 2.3: Reset values of bc.inlet
//    //
//
//    std::cout<<"XXXXXXX check 8"<<std::endl;
//    ioData.bc.inlet.density /= ioData.ref.rv.density;
//    ioData.bc.inlet.pressure /= ioData.ref.rv.pressure;
//    ioData.bc.inlet.temperature /= ioData.ref.rv.temperature;
//    ioData.bc.inlet.nutilde /= ioData.ref.rv.nutilde;
//    ioData.bc.inlet.kenergy /= ioData.ref.rv.kenergy;
//    ioData.bc.inlet.eps /= ioData.ref.rv.epsilon;
//
//    //
//    // Step 2.4: Reset values of bc.outlet
//    //
//
//    std::cout<<"XXXXXXX check 9"<<std::endl;
//    ioData.bc.outlet.density /= ioData.ref.rv.density;
//    ioData.bc.outlet.pressure /= ioData.ref.rv.pressure;
//    ioData.bc.outlet.temperature /= ioData.ref.rv.temperature;
//    ioData.bc.outlet.nutilde /= ioData.ref.rv.nutilde;
//    ioData.bc.outlet.kenergy /= ioData.ref.rv.kenergy;
//    ioData.bc.outlet.eps /= ioData.ref.rv.epsilon;
//
//    ioData.restart.etime /= ioData.ref.rv.time;
//    ioData.restart.dt_nm1 /= ioData.ref.rv.time;
//    ioData.restart.dt_nm2 /= ioData.ref.rv.time;
//    ioData.restart.energy /= ioData.ref.rv.energy;
//
//    std::cout<<"XXXXXXX check 10"<<std::endl;
//    ioData.bc.wall.temperature /= ioData.ref.rv.temperature;
//    ioData.linearizedData.stepsize = ioData.ts.timestep;
//    ioData.ts.timestep /= ioData.ref.rv.time;             // Problem in RigidRollMeshMotionHandler
//    ioData.ts.maxTime /= ioData.ref.rv.time;              // Problem in RigidRollMeshMotionHandler
//    ioData.rmesh.vx /= ioData.ref.rv.velocity;            // Problem in RigidMeshMotionHandler
//    ioData.rmesh.vy /= ioData.ref.rv.velocity;
//    ioData.rmesh.vz /= ioData.ref.rv.velocity;
//    ioData.rmesh.ax /= ioData.ref.rv.velocity / ioData.ref.rv.time;
//    ioData.rmesh.ay /= ioData.ref.rv.velocity / ioData.ref.rv.time;
//    ioData.rmesh.az /= ioData.ref.rv.velocity / ioData.ref.rv.time;
//    ioData.rmesh.timestep /= ioData.ref.rv.time;          // Problem in AccMeshMotionHandler
//
//    for (int j=0; j<ioData.rmesh.num; j++){
//      ioData.rmesh.vpts[j]->time     /= ioData.ref.rv.time;
//      ioData.rmesh.vpts[j]->velocityX /= ioData.ref.rv.velocity;
//      ioData.rmesh.vpts[j]->velocityY /= ioData.ref.rv.velocity;
//      ioData.rmesh.vpts[j]->velocityZ /= ioData.ref.rv.velocity;
//    }
//    std::cout<<"XXXXXXX check 11"<<std::endl;
//    ioData.aero.pressure /= ioData.ref.rv.pressure;
//    ioData.forced.timestep /= ioData.ref.rv.time;         // Problem in ForcedMeshMotionHandler
//    ioData.forced.frequency *= ioData.ref.rv.time;        // Problem in ForcedMeshMotionHandler
//
//    ioData.eqs.gravity_x /= ioData.ref.rv.velocity / ioData.ref.rv.time;
//    ioData.eqs.gravity_y /= ioData.ref.rv.velocity / ioData.ref.rv.time;
//    ioData.eqs.gravity_z /= ioData.ref.rv.velocity / ioData.ref.rv.time;
//    ioData.bc.hydro.depth /= ioData.ref.length;
//
//    std::cout<<"XXXXXXX check 12"<<std::endl;
//    double theta_k = 1.0;
//    double theta_w = 10.0;
//    if (kenergy0 == pow(10.0, -theta_k) * theta_w / reynolds0) {
//      ioData.bc.inlet.kenergy = pow(10.0, -theta_k) * theta_w / ioData.ref.reynolds_mu;
//      ioData.bc.inlet.eps = ioData.eqs.tc.tm.ke.c_mu * ioData.bc.inlet.kenergy * theta_w;
//    }
//
//    Pin = ioData.aero.pressure;
//
//    this->postOp->rstVarPostFcn(ioData);
//
//    this->postOp->rstVar(ioData);
//
//    if (ioData.eqs.type == EquationsData::NAVIER_STOKES)
//      this->spaceOp->rstVarFet(ioData);
//
//    this->spaceOp->rstFluxFcn(ioData);
//
//    mvp->rstSpaceOp(ioData, this->varFcn, this->spaceOp, false);
//
//    std::cout<<"\033[96m"<<__FILE__<<__LINE__<<"\033[00m"<<std::endl;//TODO delete line
//    this->rstVarImplicitSegTsDesc(ioData);//TODO HACK
//    std::cout<<"\033[96m"<<__FILE__<<__LINE__<<"\033[00m"<<std::endl;//TODO delete line
//
//    this->bcData->rstVar(ioData);
//
//    (this->timeState->getData()).rstVar(ioData);
//
//    this->timeState->rstVar(ioData);
//
//    this->output->rstVar(ioData);
//
//    this->restart->rstVar(ioData);
//
//    this->data->rstVar(ioData);
//
//    this->refVal->rstVar(ioData);
//
//    this->varFcn->rstVar(ioData);
//
//  }
//
//// initialize boundary condition fluxes
//  this->bcData->initialize(ioData, *this->X);

}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoGetEfforts(IoData &ioData,
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

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoGetDerivativeOfEffortsFiniteDifference(IoData &ioData,
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

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoGetDerivativeOfEffortsAnalytical(bool isSparse, IoData &ioData,
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

  dRdXoperators<dim> *dRdXop = dRdX->getdRdXop();
  this->spaceOp->computeDerivativeOfGradP(dRdXop, dX, dAdS, dU, dddx, dddy, dddz, dR, dGradP);

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

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoGetDerivativeOfLoadFiniteDifference(IoData &ioData, DistSVec<double,3> &X, DistSVec<double,3> &dX, DistVec<double> &A,
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

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoGetDerivativeOfLoadAnalytical(bool isSparse, IoData &ioData, DistSVec<double,3> &X, DistSVec<double,3> &dX,
                                             DistSVec<double,dim> &U, DistSVec<double,dim> &dU, DistSVec<double,3> &load, DistSVec<double,3> &dLoad)
{

  double gamma = ioData.eqs.fluidModel.gasModel.specificHeatRatio;

  double velocity = ioData.ref.mach * sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
  double dVelocity= sqrt(gamma * ioData.ref.pressure / ioData.ref.density)*DFSPAR[0];
  double dForce=2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*velocity*dVelocity;

  load=0.0;
  dLoad=0.0;
//  DistSVec<double,3> dLoad2(dLoad), diff(dLoad);
  this->spaceOp->computeGradP(X, *this->A, U);
  this->postOp->computeNodalForce(X, U, Pin, load);

  dRdXoperators<dim> *dRdXop = dRdX->getdRdXop();
  if(isSparse)
    this->spaceOp->computeDerivativeOfGradP(dRdXop, dX, dAdS, dU, dddx, dddy, dddz, dR, dGradP);
  else this->spaceOp->computeDerivativeOfGradP(X, dX, *this->A, dAdS, U, dU);
/*/////////////// checking spaceOp->computeDerivativeOfGradP & spaceOp->computeTransposeDerivativeOfGradP
 //
 //
  DistSVec<double,3> dGradP2(dGradP), dX2(dX);
  DistSVec<double,dim> dddx2(dddx), dddy2(dddy), dddz2(dddz), dU2(dU);
  DistSVec<double,6> dR2(dR);
  DistVec<double> dAdS2(dAdS);
  dGradP2 = 0;  dX2 = 0;  dddx2 = 0;  dddy2 = 0;  dddz2 = 0;  dU2 = 0;  dR2 = 0;  dAdS2 = 0;
  this->spaceOp->computeDerivativeOfGradP(dRdXop, dX, dAdS, dU, dddx2, dddy2, dddz2, dR2, dGradP2);
  double aa = dGradP2*dGradP; // + dR2*dR;

  dddx2 = 0;  dddy2 = 0;   dddz2 = 0;
  this->spaceOp->computeTransposeDerivativeOfGradP(dRdXop, dGradP, dddx2, dddy2, dddz2, dR, dAdS2, dX2, dU2);
  double bb = dAdS2*dAdS + dX2*dX + dU2*dU;

  double diffnorm = sqrt((aa-bb)*(aa-bb));
  if(aa != 0) this->com->fprintf(stderr, " ... rel. diff is %e\n", diffnorm/std::abs(aa));
  else this->com->fprintf(stderr, " ... abs. diff is %e\n", diffnorm);
*/

  //TODO: must treat dS2 better in case that dS is not zero.
  double dS2[3] = {0};

  if(isSparse) this->postOp->computeDerivativeOfNodalForce(dRdXop->dForcedX, dRdXop->dForcedGradP, dRdXop->dForcedV, dRdXop->dForcedS,
                                                           dRdXop->dVdU, dX, dGradP, dU, DFSPAR, dLoad);
  else this->postOp->computeDerivativeOfNodalForce(X, dX, U, dU, Pin, DFSPAR, dLoad);

/* ////////////////// checking computeDerivativeOfNodalForce && computeTransposeDerivativeOfNodalForce
 //
 //
  DistSVec<double,dim> dU2(dU);
  DistSVec<double,3> dGradP2(dGradP), dX2(dX), dLoad2(dLoad);
  DistVec<double> dAdS2(dAdS);
  dGradP2 = 0.0;  dX2 = 0.0;  dU2 = 0.0;  dAdS2 = 0.0;  dLoad2 = 0.0;

  this->postOp->computeDerivativeOfNodalForce(dRdXop->dForcedX, dRdXop->dForcedGradP, dRdXop->dForcedV, dRdXop->dForcedS,
                                              dRdXop->dVdU, dX, dGradP, dU, DFSPAR, dLoad2);
  double aa = dLoad2*dLoad;

  this->postOp->computeTransposeDerivativeOfNodalForce(dRdXop->dForcedX,
                                                       dRdXop->dForcedGradP,
                                                       dRdXop->dForcedV,
                                                       dRdXop->dForcedS,
                                                       dRdXop->dVdU,
                                                       dLoad, dX2, dGradP2,
                                                       dU2, dS2);

  double bb = dU2*dU + dX2*dX + dAdS2*dAdS + dGradP2*dGradP;
  double diffnorm = sqrt((aa-bb)*(aa-bb));
  if( aa != 0.0 ) this->com->fprintf(stderr, " ... final rel. diff = %e\n", diffnorm/sqrt(aa*aa));
  else this->com->fprintf(stderr, " ... final abs. diff = %e\n", diffnorm);

  this->com->fprintf(stderr, " ... dS2[0] = %e, dS2[1] = %e, dS2[2] = %e\n", dS2[0], dS2[1], dS2[2]);
*/



/* //////////////////// checking spaceOp->computeDerivativeOfGradP & spaceOp->computeTransposeDerivativeOfGradP
 ////////////////////// checking computeDerivativeOfNodalForce && computeTransposeDerivativeOfNodalForce
 //
 //
  DistSVec<double,dim> dU2(dU);
  DistSVec<double,3> dGradP2(dGradP), dX2(dX), dLoad2(dLoad);
  DistVec<double> dAdS2(dAdS);

  dGradP2 = 0.0;  dLoad2 = 0.0;  dddx = 0.0; dddy = 0.0;  dddz = 0.0;  dR = 0.0;
  this->spaceOp->computeDerivativeOfGradP(dRdXop, dX, dAdS, dU, dddx, dddy, dddz, dR, dGradP2);
  this->postOp->computeDerivativeOfNodalForce(dRdXop->dForcedX, dRdXop->dForcedGradP, dRdXop->dForcedV, dRdXop->dForcedS,
                                              dRdXop->dVdU, dX, dGradP2, dU, DFSPAR, dLoad2);
  double aa = dLoad2*dLoad;

  dddx = 0.0; dddy = 0.0;  dddz = 0.0;  dR = 0.0;   dGradP2 = 0.0;    dU2 = 0.0;   dX2 = 0.0;   dAdS2 = 0.0;
  this->postOp->computeTransposeDerivativeOfNodalForce(dRdXop->dForcedX,
                                                       dRdXop->dForcedGradP,
                                                       dRdXop->dForcedV,
                                                       dRdXop->dForcedS,
                                                       dRdXop->dVdU,
                                                       dLoad, dX2, dGradP2,
                                                       dU2, dS2);

  this->spaceOp->computeTransposeDerivativeOfGradP(dRdXop, dGradP2, dddx, dddy, dddz, dR, dAdS2, dX2, dU2);
  double bb = dU2*dU + dX2*dX + dAdS2*dAdS;
  double diffnorm = sqrt((aa-bb)*(aa-bb));
  if( aa != 0.0 ) this->com->fprintf(stderr, " ... final rel. diff = %e\n", diffnorm/sqrt(aa*aa));
  else this->com->fprintf(stderr, " ... final abs. diff = %e\n", diffnorm);

  this->com->fprintf(stderr, " ... dS2[0] = %e, dS2[1] = %e, dS2[2] = %e\n", dS2[0], dS2[1], dS2[2]);
*/


/*
  diff = dLoad2 - dLoad;
  double diffnorm = diff.norm();
  double dLoadnorm = dLoad.norm();
  double dLoad2norm = dLoad2.norm();
  if(dLoadnorm != 0) this->com->fprintf(stderr, " ... rel. diff is %e, dLoadnorm = %e, dLoad2norm = %e\n", diffnorm/dLoadnorm, dLoadnorm, dLoad2norm);
  else this->com->fprintf(stderr, " ... abs. diff is %e\n", diffnorm);
*/

  if (this->refVal->mode == RefVal::NON_DIMENSIONAL) {
    dLoad *= 2.0 * this->refVal->length*this->refVal->length / surface;
  }
  else {
    dLoad += (dForce / this->refVal->force) * load;
  }

}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoGetTransposeDerivativeOfLoadAnalytical(IoData &ioData,
                                     DistSVec<double,3> &dLoad, DistSVec<double,3> &dX, DistSVec<double,dim> &dU)
{

  double gamma = ioData.eqs.fluidModel.gasModel.specificHeatRatio;

  double velocity = ioData.ref.mach * sqrt(gamma * ioData.ref.pressure / ioData.ref.density);
  double dVelocity= sqrt(gamma * ioData.ref.pressure / ioData.ref.density)*DFSPAR[0];
  double dForce=2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*velocity*dVelocity;

  dX = 0.0;
  dU = 0.0;

  dRdXoperators<dim> *dRdXop = dRdX->getdRdXop();

  //TODO: must treat DFSPAR2 better in case that it is not zero.
  double DFSPAR2[3] = {0};

  DistSVec<double,dim> dU2(dU);
  DistSVec<double,3> dGradP2(dGradP), dX2(dX), dLoad2(dLoad);
  DistVec<double> dAdS2(dAdS);

  dddx = 0.0; dddy = 0.0;  dddz = 0.0;  dR = 0.0;   dGradP = 0.0;    dAdS2 = 0.0;
  this->postOp->computeTransposeDerivativeOfNodalForce(dRdXop->dForcedX,
                                                       dRdXop->dForcedGradP,
                                                       dRdXop->dForcedV,
                                                       dRdXop->dForcedS,
                                                       dRdXop->dVdU,
                                                       dLoad, dX, dGradP2,
                                                       dU, DFSPAR2);

  this->spaceOp->computeTransposeDerivativeOfGradP(dRdXop, dGradP2, dddx, dddy, dddz, dR, dAdS2, dX, dU);

  if (this->refVal->mode == RefVal::NON_DIMENSIONAL) {
    dX *= 2.0 * this->refVal->length*this->refVal->length / surface;
    dU *= 2.0 * this->refVal->length*this->refVal->length / surface;
  }
  else {
      //TODO: needs to add the term below if Mach number is used as a sensitivity variable.
//    dLoad += (dForce / this->refVal->force) * load;
  }

}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoSemiAnalytical
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
  if(DFSPAR[1] || DFSPAR[2]) dF *= 0.017453292519943295769236907684886127134428718885417254560;  // convert radian to degree
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

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoAnalytical
(bool isSparse, IoData &ioData, DistSVec<double,3> &X, DistVec<double> &A, DistSVec<double,dim> &U, DistSVec<double,dim> &dFdS)
{

  //
  // Computing the normal, derivative of the normal and of the control volume
  //


  DistSVec<double,dim> dddx2(dddx), dddy2(dddy), dddz2(dddz);
  DistSVec<double,6> dR2(dR);
  DistVec<double> dFaceNormVel(domain->getFaceNormDistInfo());
  dEdgeNorm = 0.0;  dFaceNorm = 0.0;  dddx = 0.0;  dddy = 0.0;  dddz = 0.0;  dFaceNormVel = 0.0;  dR = 0.0;


  dRdXoperators<dim> *dRdXop = dRdX->getdRdXop();
  if(isSparse) {
    this->geoState->computeDerivatives(dRdXop->dEdgeNormdX, dRdXop->dFaceNormdX, dRdXop->dCtrlVoldX, X, dXdS, dAdS, dEdgeNorm, dFaceNorm, dFaceNormVel);
  } else
    this->geoState->computeDerivatives(X, dXdS, this->bcData->getVelocityVector(), this->bcData->getDerivativeOfVelocityVector(), dAdS);


  //
  // Computing the derivatives of the boundary fluxes
  //
// TODO:: uncomment this!
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
  } else
       this->spaceOp->computeDerivativeOfResidual(X, dXdS, A, dAdS, U, DFSPAR[0], Flux, dFdS, this->timeState);

/*  checking   BOTH   geoState->computeDerivatives   &&   spaceOp->computeDerivativeOfResidual
//
    DistSVec<double,3> dXdS2(dXdS);
    DistSVec<double,dim> dFdS2(dFdS);
    DistVec<double> dFaceNormVel2(dFaceNormVel), dAdS2(dAdS);
    DistVec<Vec3D> dEdgeNorm2(dEdgeNorm), dFaceNorm2(dFaceNorm);

    dFdS2 = 0.0;   dR2 = 0.0;   dddx2 = 0.0;  dddy2 = 0.0;   dddz2 = 0.0;  dAdS2 = 0.0;    dEdgeNorm2 = 0.0;   dFaceNorm2 = 0.0;   dFaceNormVel2 = 0.0;
    this->geoState->computeDerivatives(dRdXop->dEdgeNormdX, dRdXop->dFaceNormdX, dRdXop->dCtrlVoldX, dXdS, dAdS2, dEdgeNorm2, dFaceNorm2, dFaceNormVel2);
    this->spaceOp->computeDerivativeOfResidual(dRdXop, dXdS, dAdS2, dEdgeNorm2, dFaceNorm2, dFaceNormVel2, dFdS2, dR2, dddx2, dddy2, dddz2);
    double aa = dFdS2*dFdS;

    dXdS2 = 0.0;   dR2 = 0.0;   dddx2 = 0.0;  dddy2 = 0.0;   dddz2 = 0.0;  dAdS2 = 0.0;    dEdgeNorm2 = 0.0;   dFaceNorm2 = 0.0;   dFaceNormVel2 = 0.0;
    this->spaceOp->computeTransposeDerivativeOfResidual(dRdXop, dFdS, dAdS2, dXdS2, dddx2, dddy2, dddz2, dEdgeNorm2, dFaceNorm2, dFaceNormVel2, dR2);
    this->geoState->computeTransposeDerivatives(dRdXop->dEdgeNormdX, dRdXop->dFaceNormdX, dRdXop->dCtrlVoldX, dAdS2, dEdgeNorm2, dFaceNorm2, dFaceNormVel2, dXdS2);

    double bb = dXdS2*dXdS;
    double diff = sqrt((aa-bb)*(aa-bb));
    if(aa != 0) this->com->fprintf(stderr, " ... dFlux/dX relative error = %e, aa = %e, bb = %e\n", diff/abs(aa), aa, bb);
    else this->com->fprintf(stderr, " ... dFlux/dX absolute error = %e, aa = %e, bb = %e\n", diff, aa, bb);
*/

/*  checking   spaceOp->computeDerivativeOfResidual
    DistSVec<double,3> dXdS2(dXdS);
    DistSVec<double,dim> dFdS2(dFdS);
    DistVec<double> dFaceNormVel2(dFaceNormVel), dAdS2(dAdS);
    DistVec<Vec3D> dEdgeNorm2(dEdgeNorm), dFaceNorm2(dFaceNorm);

    dFdS2 = 0.0;   dR2 = 0.0;   dddx2 = 0.0;  dddy2 = 0.0;   dddz2 = 0.0;
    this->spaceOp->computeDerivativeOfResidual(dRdXop, dXdS, dAdS, dEdgeNorm, dFaceNorm, dFaceNormVel, dFdS2, dR2, dddx2, dddy2, dddz2);
    double aa = dFdS2*dFdS;

    dXdS2 = 0.0;   dR2 = 0.0;   dddx2 = 0.0;  dddy2 = 0.0;   dddz2 = 0.0;  dAdS2 = 0.0;    dEdgeNorm2 = 0.0;   dFaceNorm2 = 0.0;   dFaceNormVel2 = 0.0;
    this->spaceOp->computeTransposeDerivativeOfResidual(dRdXop, dFdS, dAdS2, dXdS2, dddx2, dddy2, dddz2, dEdgeNorm2, dFaceNorm2, dFaceNormVel2, dR2);


    DistSVec<double,3> dEdgeNormSVec(dEdgeNorm.info()), dEdgeNorm2SVec(dEdgeNorm2.info());
    DistSVec<double,3> dFaceNormSVec(dFaceNorm.info()), dFaceNorm2SVec(dFaceNorm2.info());
    for(int iSub=0; iSub< dEdgeNorm.info().numLocThreads; iSub++)
      for(int i=0; i<dEdgeNorm[iSub]->size(); ++i)
        for(int j=0; j<3; ++j) {
          dEdgeNormSVec(iSub)[i][j] = dEdgeNorm(iSub)[i][j];
          dEdgeNorm2SVec(iSub)[i][j] = dEdgeNorm2(iSub)[i][j];
        }
    for(int iSub=0; iSub< dFaceNorm.info().numLocThreads; iSub++)
      for(int i=0; i<dFaceNorm[iSub]->size(); ++i)
        for(int j=0; j<3; ++j) {
          dFaceNormSVec(iSub)[i][j] = dFaceNorm(iSub)[i][j];
          dFaceNorm2SVec(iSub)[i][j] = dFaceNorm2(iSub)[i][j];
        }

    double bb = dXdS2*dXdS + dAdS2*dAdS + dEdgeNormSVec*dEdgeNorm2SVec + dFaceNormSVec*dFaceNorm2SVec + dFaceNormVel2*dFaceNormVel;
    double diff = sqrt((aa-bb)*(aa-bb));
    if(aa != 0) this->com->fprintf(stderr, " ... computeTransposeDerivativeOfResidual ... relative error = %e, aa = %e, bb = %e\n", diff/abs(aa), aa, bb);
    else this->com->fprintf(stderr, " ... computeTransposeDerivativeOfResidual ... absolute error = %e, aa = %e, bb = %e\n", diff, aa, bb);
*/

/*  checking   geoState->computeDerivatives
    DistSVec<double,3> dXdS2(dXdS);
    DistSVec<double,dim> dFdS2(dFdS);
    DistVec<double> dFaceNormVel2(dFaceNormVel), dAdS2(dAdS);
    DistVec<Vec3D> dEdgeNorm2(dEdgeNorm), dFaceNorm2(dFaceNorm);

    dAdS2 = 0.0;    dEdgeNorm2 = 0.0;   dFaceNorm2 = 0.0;   dFaceNormVel2 = 0.0;
    this->geoState->computeDerivatives(dRdXop->dEdgeNormdX, dRdXop->dFaceNormdX, dRdXop->dCtrlVoldX, dXdS, dAdS2, dEdgeNorm2, dFaceNorm2, dFaceNormVel2);
    DistSVec<double,3> dEdgeNormSVec(dEdgeNorm.info()), dEdgeNorm2SVec(dEdgeNorm2.info());
    DistSVec<double,3> dFaceNormSVec(dFaceNorm.info()), dFaceNorm2SVec(dFaceNorm2.info());
    for(int iSub=0; iSub< dEdgeNorm.info().numLocThreads; iSub++)
      for(int i=0; i<dEdgeNorm[iSub]->size(); ++i)
        for(int j=0; j<3; ++j) {
          dEdgeNormSVec(iSub)[i][j] = dEdgeNorm(iSub)[i][j];
          dEdgeNorm2SVec(iSub)[i][j] = dEdgeNorm2(iSub)[i][j];
        }
    for(int iSub=0; iSub< dFaceNorm.info().numLocThreads; iSub++)
      for(int i=0; i<dFaceNorm[iSub]->size(); ++i)
        for(int j=0; j<3; ++j) {
          dFaceNormSVec(iSub)[i][j] = dFaceNorm(iSub)[i][j];
          dFaceNorm2SVec(iSub)[i][j] = dFaceNorm2(iSub)[i][j];
        }
    double aa = dAdS2*dAdS + dEdgeNorm2SVec*dEdgeNormSVec + dFaceNorm2SVec*dFaceNormSVec + dFaceNormVel2*dFaceNormVel;

    dXdS2 = 0.0;
    this->geoState->computeTransposeDerivatives(dRdXop->dEdgeNormdX, dRdXop->dFaceNormdX, dRdXop->dCtrlVoldX, dAdS, dEdgeNorm, dFaceNorm, dFaceNormVel, dXdS2);
    double bb = dXdS2*dXdS;
    double diff = sqrt((aa-bb)*(aa-bb));
    if(aa != 0) this->com->fprintf(stderr, " ... relative error = %e, aa = %e, bb = %e\n", diff/abs(aa), aa, bb);
    else this->com->fprintf(stderr, " ... absolute error = %e, aa = %e, bb = %e\n", diff, aa, bb);
*/


/*
  DistSVec<double,dim> dFdS2(dFdS), diff(dFdS);
  this->spaceOp->computeDerivativeOfResidual(X, dXdS, A, dAdS, U, DFSPAR[0], Flux, dFdS2, this->timeState, false, dRdXop);

  diff = dFdS-dFdS2;
  double dFdSnorm(0), dFdS2norm(0), diffnorm(0);
  dFdSnorm = dFdS.norm();
  dFdS2norm = dFdS2.norm();
  diffnorm = diff.norm();
  if(dFdSnorm != 0) this->com->fprintf(stderr, "... rel. error = %e\n", diffnorm/dFdSnorm);
  else this->com->fprintf(stderr, "... abs. error = %e\n", diffnorm);
*/

  this->spaceOp->applyBCsToDerivativeOfResidual(U, dFdS);
  if(DFSPAR[1] || DFSPAR[2]) dFdS *= 0.017453292519943295769236907684886127134428718885417254560;  // convert radian to degree
}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoAnalyticalTranspose
(DistSVec<double,dim> &dFdS, DistSVec<double,3> &dXdS)
{

  dEdgeNorm = 0.0;  dFaceNorm = 0.0;  dddx = 0.0;  dddy = 0.0;  dddz = 0.0;  dFaceNormVel = 0.0;  dR = 0.0;
  dRdXoperators<dim> *dRdXop = dRdX->getdRdXop();
  this->spaceOp->computeTransposeDerivativeOfResidual(dRdXop, dFdS, dAdS, dXdS, dddx, dddy, dddz, dEdgeNorm, dFaceNorm, dFaceNormVel, dR);

  //
  // Computing the derivatives of the boundary fluxes
  //
// TODO:: uncomment this!
//  this->bcData->initializeSA(ioData, X, dXdS, DFSPAR[0], DFSPAR[1], DFSPAR[2]);

  //
  // Computing the partial derivative of the flux with respect to the variables
  //

  this->geoState->computeTransposeDerivatives(dRdXop->dEdgeNormdX, dRdXop->dFaceNormdX, dRdXop->dCtrlVoldX, dAdS, dEdgeNorm, dFaceNorm, dFaceNormVel, dXdS);

//  this->spaceOp->applyBCsToDerivativeOfResidual(U, dFdS);
//  if(DFSPAR[1] || DFSPAR[2]) dFdS *= 0.017453292519943295769236907684886127134428718885417254560;  // convert radian to degree
}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoSetUpLinearSolver(IoData &ioData, DistSVec<double,3> &X, DistVec<double> &A,
                                                              DistSVec<double,dim> &U, DistSVec<double,dim> &dFdS)
{

  //int* debugtest = NULL;
  //if(debugtest==NULL)//TODO delete lines
  //  Dev::Error(this->com,"Variable 'debugtest' not initialized");

// Preparing the linear solver
  fsaRestartBcFluxs(ioData);//TODO HACK

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


template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoLinearSolver
(
  IoData &ioData,
  DistSVec<double,dim> &dFdS, DistSVec<double,dim> &dUdS,
  bool isFSI
)
{

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

}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoPrintTextOnScreen(const char *Text)
{
   this->com->fprintf(stderr, Text);
}


//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsaPrintTextOnScreen(const char *Text)
{
   this->com->fprintf(stderr, Text);
}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoOutput1D(const char *fileName, DistVec<double> &V)
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

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoOutput3D(const char *fileName, DistSVec<double,3> &V)
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


template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoOutputDimD(const char *fileName, DistSVec<double,dim> &V)
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

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoInitialize(IoData &ioData, DistSVec<double,dim> &U)
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

  // Setting up the linear solver
  fsoSetUpLinearSolver(ioData, *this->X, *this->A, U, dFdS);
}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
int FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoHandler(IoData &ioData, DistSVec<double,dim> &U)
{
  std::cout<<"========"<<__FILE__<<":"<<__LINE__<<std::endl;//TODO delete line

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
  bool isSparse = false;

  fsoSetUpLinearSolver(ioData, *this->X, *this->A, U, dFdS);

  if (ioData.sa.sensMesh == SensitivityAnalysis::ON_SENSITIVITYMESH) fso_on_sensitivityMesh(isSparse, ioData, U);
  if (ioData.sa.sensMach == SensitivityAnalysis::ON_SENSITIVITYMACH) fso_on_sensitivityMach(isSparse, ioData, U);
  if (ioData.sa.sensAlpha == SensitivityAnalysis::ON_SENSITIVITYALPHA) fso_on_sensitivityAlpha(isSparse, ioData, U);
  if (ioData.sa.sensBeta == SensitivityAnalysis::ON_SENSITIVITYBETA) fso_on_sensitivityBeta(isSparse, ioData, U);

  bool lastIt = true;
//  this->outputToDisk(ioData, &lastIt, 0, 0, 0, 0, dtLeft, U);
//  this->outputPositionVectorToDisk(U);

  this->output->closeAsciiFiles();


//  this->com->barrier();
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

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::setDFSPAR(IoData &ioData)
{
  switch (actvar) {
    case 1:
      if(ioData.sa.sensMesh == SensitivityAnalysis::OFF_SENSITIVITYMESH) { fsoPrintTextOnScreen(" ***** Error: SensitivityMesh in fluid input must be on\n"); exit(-1); }
      DFSPAR[0] = 0;  DFSPAR[1] = 0;  DFSPAR[2] = 0;
      break;
    case 2:
      if(ioData.sa.sensMach == SensitivityAnalysis::OFF_SENSITIVITYMACH) { fsoPrintTextOnScreen(" ***** Error: SensitivityMach in fluid input must be on\n"); exit(-1); }
      DFSPAR[0] = 1;  DFSPAR[1] = 0;  DFSPAR[2] = 0;
      break;
    case 3:
      if(ioData.sa.sensAlpha == SensitivityAnalysis::OFF_SENSITIVITYALPHA) { fsoPrintTextOnScreen(" ***** Error: SensitivityAlpha in fluid input must be on\n"); exit(-1); }
      DFSPAR[0] = 0;  DFSPAR[1] = 1;  DFSPAR[2] = 0;
      break;
    case 4:
      if(ioData.sa.sensBeta == SensitivityAnalysis::OFF_SENSITIVITYBETA) { fsoPrintTextOnScreen(" ***** Error: SensitivityBeta in fluid input must be on\n"); exit(-1); }
      DFSPAR[0] = 0;  DFSPAR[1] = 0;  DFSPAR[2] = 1;
      break;
    case 5:
      if(ioData.sa.sensFSI == SensitivityAnalysis::OFF_SENSITIVITYFSI) { fsoPrintTextOnScreen(" ***** Error: SensitivityFSI in fluid input must be on\n"); exit(-1); }
      DFSPAR[0] = 0;  DFSPAR[1] = 0;  DFSPAR[2] = 0;
      break;
    default:
      fsoPrintTextOnScreen(" ***** Error: invalid value for active sensitivity variable!\n"); exit(-1);
      break;
  }
}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
int FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoAeroelasticHandler(IoData &ioData, DistSVec<double,dim> &U)
{

  // xmach      -  Mach number
  // alpha      -  pitch angle
  // teta       -  yaw angle
  // DFSPAR(1)  -  Mach number differential
  // DFSPAR(2)  -  angle of attack differential
  // DFSPAR(3)  -  yaw angle differential

  // Start basic timer
  double MyLocalTimer = -this->timer->getTime();

  bool isSparse = false;
  double dtLeft = 0.0;
  this->computeTimeStep(1, &dtLeft, U);
  this->computeMeshMetrics();
  this->updateStateVectors(U);

  fsoSetUpLinearSolver(ioData, *this->X, *this->A, U, dFdS);
  int totalNumParamTypes;
  this->getNumParam(totalNumParamTypes,actvar,steadyTol);
  if(ioData.sa.sensMach == SensitivityAnalysis::ON_SENSITIVITYMACH) {  totalNumParamTypes++; }
  if(ioData.sa.sensAlpha == SensitivityAnalysis::ON_SENSITIVITYALPHA) { totalNumParamTypes++; }
  if(ioData.sa.sensBeta == SensitivityAnalysis::ON_SENSITIVITYBETA) { totalNumParamTypes++; }

  if(isSparse)  dRdX->constructOperators(*this->X, *this->A, U, DFSPAR[0], Flux, Pin, this->timeState, this->postOp);

  for(int iparam=0; iparam<totalNumParamTypes; ++iparam) {
    int numParam;
    this->getNumParam(numParam,actvar,steadyTol);
    setDFSPAR(ioData);
    for(int i=0; i<numParam; ++i) fso_on_aeroelasticSensitivityFSI(isSparse, ioData, U);
  }
  bool lastIt = true;
//  this->outputToDisk(ioData, &lastIt, 0, 0, 0, 0, dtLeft, U);
//  this->outputPositionVectorToDisk(U);

  this->output->closeAsciiFiles();


//  this->com->barrier();
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

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fso_on_sensitivityBeta(bool isSparse, IoData &ioData, DistSVec<double,dim> &U)
{
    dXdS = 0.0;
    dAdS = 0.0;
    DFSPAR[0] = 0.0;
    DFSPAR[1] = 0.0;
    DFSPAR[2] = 1.0;
    actvar = 4;

    if (!ioData.sa.angleRad)
      ioData.sa.eps *= acos(-1.0) / 180.0;

    fsoComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U, false, isSparse);

    fsoComputeSensitivities(isSparse, ioData, "Derivatives with respect to the yaw angle:", ioData.sa.sensoutput, *this->X, U);

    fsoPrintTextOnScreen("\n ***** Derivatives with respect to the yaw angle were computed! \n");

    step = step + 1;

    if (!ioData.sa.angleRad)
      ioData.sa.eps /= acos(-1.0) / 180.0;
}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fso_on_sensitivityAlpha(bool isSparse, IoData &ioData, DistSVec<double,dim> &U)
{
    dXdS = 0.0;
    dAdS = 0.0;
    DFSPAR[0] = 0.0;
    DFSPAR[1] = 1.0;
    DFSPAR[2] = 0.0;
    actvar = 3;

    if (!ioData.sa.angleRad)
      ioData.sa.eps *= acos(-1.0) / 180.0;

    fsoComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U, false, isSparse);

    fsoComputeSensitivities(isSparse, ioData, "Derivatives with respect to the angle of attack:", ioData.sa.sensoutput, *this->X, U);

    fsoPrintTextOnScreen("\n ***** Derivatives with respect to the angle of attack were computed! \n");

    step = step + 1;

    if (!ioData.sa.angleRad)
      ioData.sa.eps /= acos(-1.0) / 180.0;
}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fso_on_sensitivityMach(bool isSparse, IoData &ioData, DistSVec<double,dim> &U)
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

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fso_on_aeroelasticSensitivityFSI(bool isSparse, IoData &ioData, DistSVec<double,dim> &U)
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
//      this->com->fprintf(stderr, "norm of dXdSb is %e\n",dXdSb.norm());

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
//      else this->com->fprintf(stderr, "\n norm of dXdS is %e\n", dXdS.norm());

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

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fso_on_sensitivityMesh(bool isSparse, IoData &ioData, DistSVec<double,dim> &U)
{

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

      if ( ioData.input.shapederivativesType == InputData::WALL) {
        // Reading derivative of the overall deformation
        bool readOK = domain->readVectorFromFile(this->input->shapederivatives, step, &tag, dXdSb);
        if(!readOK) break;

        // Checking if dXdSb has entries different from zero at the interior of the mesh
        this->postOp->checkVec(dXdSb);

      if (dXdSb.norm() == 0.0) {
          this->com->fprintf(stderr, "\n *** WARNING *** No Mesh Perturbation \n\n");
          if(!ioData.sa.fsiFlag) exit(1);
        }


        // Updating the mesh
        dXdS = *this->X;
        mms->solve(dXdSb, dXdS);
        dXdS -= *this->X;
      } else if ( ioData.input.shapederivativesType == InputData::VOLUME) {
        // Reading derivative of the overall deformation
        bool readOK = domain->readVectorFromFile(this->input->shapederivatives, step, &tag, dXdS);
        if(!readOK) break;
      }

      // Check that the mesh perturbation is propagated
      if (dXdS.norm() == 0.0) this->com->fprintf(stderr, "\n !!! WARNING !!! No Mesh Sensitivity Perturbation !!!\n\n");

      this->com->fprintf(stderr, "\n ***** Shape variable %d\n", step);

      fsoComputeDerivativesOfFluxAndSolution(ioData, *this->X, *this->A, U, false, isSparse);
      fsoComputeSensitivities(isSparse, ioData, "Derivatives with respect to the mesh position:", ioData.sa.sensoutput, *this->X, U);

      dXdSb = 0.0;
      step = step + 1;
    }

    fsoPrintTextOnScreen("\n ***** Derivatives of mesh position were computed! \n");

}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoComputeDerivativesOfFluxAndSolution(IoData &ioData, DistSVec<double,3> &X, DistVec<double> &A, DistSVec<double,dim> &U, bool isFSI, bool isSparse)
{

  dFdS = 0.0;

  // Derivative of the Flux, either analytical or semi-analytical
  if ( ioData.sa.scFlag == SensitivityAnalysis::ANALYTICAL ) {
    fsoAnalytical(isSparse, ioData, X, A, U, dFdS);
/*    dFdSref = 0.0;
    fsoSemiAnalytical(ioData, X, A, U, dFdSref);
    DistSVec<double,dim> difference(domain->getNodeDistInfo());
    difference = dFdS - dFdSref;
    this->com->fprintf(stderr, "!!! dFdS and dFdSref do not match. The relative difference norm is %e\n", difference.norm()/dFdS.norm());
    this->com->fprintf(stderr, "!!! dFdSref norm is %e\n", dFdSref.norm());  */
  } else {
    fsoSemiAnalytical(ioData, X, A, U, dFdS);
  }

  // Computing the derivative of the fluid variables
  // with respect to the fsoimization variables
  fsoLinearSolver(ioData, dFdS, dUdS,isFSI);

}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoComputeAndSendForceSensitivities(bool isSparse, IoData &ioData, const char *fileName,
                                                                             DistSVec<double,3> &X, DistSVec<double,dim> &U)
{

//  if ( ioData.sa.sensFSI == SensitivityAnalysis::ON_SENSITIVITYFSI ) {
    if (ioData.sa.scFlag == SensitivityAnalysis::FINITEDIFFERENCE ) {
      fsoGetDerivativeOfLoadFiniteDifference(ioData, X, dXdS, *this->A, U, dUdS, *load, *dLoad);
    } else {
      fsoGetDerivativeOfLoadAnalytical(isSparse, ioData, X, dXdS, U, dUdS, *load, *dLoad);
/*      *dLoadref = 0.0;
      fsoGetDerivativeOfLoadFiniteDifference(ioData, X, dXdS, *this->A, U, dUdS, *load, *dLoadref);
      DistSVec<double,3> difference(domain->getNodeDistInfo());
      difference = *dLoad - *dLoadref;
      this->com->fprintf(stderr, "!!! dLoad and dLoadref do not match. The relative difference norm is %e\n", difference.norm()/dLoad->norm());
      this->com->fprintf(stderr, "!!! dLoadref norm is %e\n", dLoadref->norm()); */
    }
//  }

  this->sendForceSensitivity(dLoad);

}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
void FluidSegShapeOptimizationHandler<dim,neq1,neq2>::fsoComputeSensitivities(bool isSparse,
                                                                 IoData &ioData, const char *mesage, const char *fileName,
                                                                 DistSVec<double,3> &X, DistSVec<double,dim> &U)
{

// Computing efforts (F: force, M: moment, L:LiftAndDrag)
  Vec3D F, M, L;
  fsoGetEfforts(ioData, X, U, F, M, L);

// Computing derivative of the efforts
  Vec3D dFds, dMds, dLds;

  if ( ioData.sa.scFlag == SensitivityAnalysis::FINITEDIFFERENCE )
    fsoGetDerivativeOfEffortsFiniteDifference(ioData, X, dXdS, *this->A, U, dUdS, dFds, dMds);
  else
    fsoGetDerivativeOfEffortsAnalytical(isSparse, ioData, X, dXdS, U, dUdS, dFds, dMds, dLds);

  if ((!ioData.sa.angleRad) && (DFSPAR[1] || DFSPAR[2])) {
    dFds *= acos(-1.0) / 180.0;
    dMds *= acos(-1.0) / 180.0;
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
