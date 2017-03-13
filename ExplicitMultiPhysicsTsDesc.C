#include <IoData.h>
#include <Domain.h>
#include <GeoSource.h>
#include <DistTimeState.h>
#include <SpaceOperator.h>
#include <LevelSet.h>
#include <FluidSelector.h>

#ifdef TYPE_MAT
#define MatScalar TYPE_MAT
#else
#define MatScalar double
#endif

#ifdef TYPE_PREC
#define PrecScalar TYPE_PREC
#else
#define PrecScalar double
#endif

//------------------------------------------------------------------------------

template<int dim, int dimLS>
ExplicitMultiPhysicsTsDesc<dim,dimLS>::
ExplicitMultiPhysicsTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  MultiPhysicsTsDesc<dim,dimLS>(ioData, geoSource, dom),
  k1(this->getVecInfo()), k2(this->getVecInfo()),
  k3(this->getVecInfo()), k4(this->getVecInfo()), Ubc(this->getVecInfo()), 
  p1(this->getVecInfo()), p2(this->getVecInfo()),
  p3(this->getVecInfo()), p4(this->getVecInfo()), 
  U0(this->getVecInfo()), Phi0(this->getVecInfo()),
  fluidId0(this->getVecInfo())
{
  timeType = ioData.ts.expl.type;
  
  //initialize mmh (EmbeddedMeshMotionHandler).
  if(this->dynNodalTransfer) {
    MeshMotionHandler *_mmh = 0;
    _mmh = new EmbeddedMeshMotionHandler(ioData, dom, this->dynNodalTransfer, this->distLSS);
    this->mmh = _mmh;
  } else this->mmh = 0;


  if (this->modifiedGhidaglia) {

    hh1 = new DistVec<double>(*this->bcData->getBoundaryStateHH());
    hh2 = new DistVec<double>(*this->bcData->getBoundaryStateHH());
    hh3 = new DistVec<double>(*this->bcData->getBoundaryStateHH());
    hh4 = new DistVec<double>(*this->bcData->getBoundaryStateHH());
    hhorig = new DistVec<double>(*this->bcData->getBoundaryStateHH());
  }

}

//------------------------------------------------------------------------------
  
template<int dim, int dimLS>
ExplicitMultiPhysicsTsDesc<dim,dimLS>::~ExplicitMultiPhysicsTsDesc()
{ /* nothing to destroy */ } 
  
//------------------------------------------------------------------------------

template<int dim, int dimLS>
int ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNonLinearSystem(DistSVec<double,dim> &U, int)
{
  solveNLSystemTwoBlocks(U);
  return 1;
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNLSystemTwoBlocks(DistSVec<double,dim> &U)
{
  //Vec3D p(0,0,0);
  //DebugTools::PrintFluidId("Printing Id 1: ", *(this->fluidSelector.fluidId), *(this->X), p, this->com->cpuNum());
  DistVec<int> fluidId_copy(*this->fluidSelector.fluidId);

  if(this->mmh && (!this->inSubCycling || this->inRedoTimestep)) {
    // get structural time-step and recompute FS intersections.
    if(!this->inSubCycling) recomputeIntersections();

    //DebugTools::PrintFluidId("Printing Id 2: ", *(this->fluidSelector.fluidId), *(this->X), p, this->com->cpuNum());

    // update fluidId.
    updateFluidIdFS(U);

    //DebugTools::PrintFluidId("Printing Id 3: ", *(this->fluidSelector.fluidId), *(this->X), p, this->com->cpuNum());

    // update the phase-change (U & Phi) caused by the motion of FS interface
    updatePhaseChangeFS(U);

    //DebugTools::PrintFluidId("Printing Id 4: ", *(this->fluidSelector.fluidId), *(this->X), p, this->com->cpuNum());
  }
//  // populate ghost nodes (only for Navier-Stokes.)
//  populateGhostPointsForNavierStokes(U);
  // evolve the fluid equation using FE, RK2, or RK4
  solveNLNavierStokes(U);

  this->errorHandler->reduceError();
  this->data->resolveErrors();
  if(this->errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP]) {
    *this->fluidSelector.fluidId = fluidId_copy;
    this->inRedoTimestep = true;
    return;
  }
  else {
    this->inRedoTimestep = false;
  }

  //DebugTools::PrintFluidId("Printing Id 5: ", *(this->fluidSelector.fluidId), *(this->X), p, this->com->cpuNum());

  // evolve the level-set equation using FE, RK2, or RK4.
  solveNLLevelSet(U);

  //DebugTools::PrintFluidId("Printing Id 6: ", *(this->fluidSelector.fluidId), *(this->X), p, this->com->cpuNum());

  // update fluidId (fluidId0 = fluidId, fluidId = new).
  fluidId0 = *(this->fluidSelector.fluidId); // used in updatePhaseChangeFF
//  DebugTools::PrintElement("Phi_post_update",this->Phi,63,0,503);

  //DebugTools::PrintFluidId("Printing Id 7: ", *(this->fluidSelector.fluidId), *(this->X), p, this->com->cpuNum());

  if(this->withCracking && this->withMixedLS) {
//    this->com->fprintf(stderr,"calling updateFluidIdFF2!\n");
    this->fluidSelector.updateFluidIdFF2(this->distLSS, this->Phi);
  } else {
//    this->com->fprintf(stderr,"calling updateFluidIdFF!\n");
    this->fluidSelector.updateFluidIdFF(this->distLSS, this->Phi);
  }
  // update the phase-change (only U) caused by the movement of FF interface

  //DebugTools::PrintFluidId("Printing Id 8: ", *(this->fluidSelector.fluidId), *(this->X), p, this->com->cpuNum());

  updatePhaseChangeFF(U);
  //this->com->fprintf(stderr,"DONE with updatePhaseChangeFF!\n");
  // check the consistency of Phi and FluidId. Can be removed for better efficiency! 
  if (this->lsMethod == 0)
    this->LS->conservativeToPrimitive(this->Phi, this->PhiV, U);
//  if(this->withCracking && this->withMixedLS)
//    this->domain->debugMultiPhysics(*(this->distLSS), this->PhiV, *(this->fluidSelector.fluidId), U);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::recomputeIntersections()
{
  // just to get structural time-step (dts).
  this->dts = this->mmh->update(0, 0, 0, this->bcData->getVelocityVector(), *this->Xs);

  this->com->barrier();
  double tw = this->timer->getTime();
  if(this->withCracking && this->withMixedLS) // no need for the intersector to determine fluidId.
    this->distLSS->recompute(this->dtf, this->dtfLeft, this->dts, false, TsDesc<dim>::failSafeFlag); 
  else
    this->distLSS->recompute(this->dtf, this->dtfLeft, this->dts, true, TsDesc<dim>::failSafeFlag); 

  this->timer->addIntersectionTime(tw);
  this->com->barrier();
  this->timer->removeIntersAndPhaseChange(tw);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::updatePhaseChangeFS(DistSVec<double,dim> &U)
{
  double tw = this->timer->getTime();
  this->multiPhaseSpaceOp->updateSweptNodes(*this->X, this->phaseChangeChoice, U, this->Vtemp,
                                            *this->Weights, *this->VWeights, this->Phi, this->PhiWeights,
                                            *this->Wstarij, *this->Wstarji, this->distLSS, this->vfar,
                                            (this->withCracking && this->withMixedLS),
                                            this->fluidSelector.fluidIdn, this->fluidSelector.fluidId);
  this->timer->addEmbedPhaseChangeTime(tw);
  this->timer->removeIntersAndPhaseChange(tw);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::updateFluidIdFS(DistSVec<double,dim> &U)
{
  if (this->lsMethod == 0)
    this->LS->conservativeToPrimitive(this->Phi, this->PhiV, U);
  else
    this->PhiV = this->Phi;

  if(this->withCracking && this->withMixedLS) {
    //this->multiPhaseSpaceOp->extrapolatePhiV2(this->distLSS, this->PhiV);
    //this->fluidSelector.updateFluidIdFS2(this->distLSS, this->PhiV);
//    this->com->fprintf(stderr,"calling updateFluidIdFS2!\n");
    DistSVec<bool,4> poll(this->domain->getNodeDistInfo());
    this->domain->updateFluidIdFS2Prep(*(this->distLSS), this->PhiV, *(this->fluidSelector.fluidId), poll);
    this->fluidSelector.updateFluidIdFS2(this->distLSS, this->PhiV, poll);
  } else {
    //this->multiPhaseSpaceOp->extrapolatePhiV(this->distLSS, this->PhiV);
//    this->com->fprintf(stderr,"calling fluidSelector.updateFluidIdFS!\n");
    this->fluidSelector.updateFluidIdFS(this->distLSS, this->PhiV);
  }
  this->PhiV = 0.0; //PhiV is no longer a distance function now. Only its sign (+/-)
                    //  is meaningful. We destroy it so people wouldn't use it
                    //  by mistake later on.
}

//------------------------------------------------------------------------------
//TODO: check with Adam to make sure this function works for Multi-Phase Fluid-Structure
template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::populateGhostPointsForNavierStokes(DistSVec<double,dim> &U)
{
  if(this->eqsType == MultiPhysicsTsDesc<dim,dimLS>::NAVIER_STOKES) {
    this->ghostPoints->deletePointers();
    this->multiPhaseSpaceOp->populateGhostPoints(this->ghostPoints,*this->X,U,this->varFcn,this->distLSS,this->viscSecOrder,*(this->fluidSelector.fluidId));
  }
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNLNavierStokes(DistSVec<double,dim> &U)
{
  double t0 = this->timer->getTime();
  switch (timeType) {
    case ExplicitData::FORWARD_EULER :
      solveNLNavierStokesFE(U); break;
    case ExplicitData::RUNGE_KUTTA_2 :
      solveNLNavierStokesRK2(U); break;
    case ExplicitData::RUNGE_KUTTA_4 :
      solveNLNavierStokesRK4(U); break;
    default:
      this->com->fprintf(stderr,"ERROR: Choose time-integrator from ForwardEuler, RungeKutta2, and RungeKutta4!\n");
  }

  this->updateBoundaryExternalState();

  // for FF phase-change update using extrapolation
  this->varFcn->conservativeToPrimitive(U,this->V0,this->fluidSelector.fluidId);
  this->riemann->storePreviousPrimitive(this->V0, *this->fluidSelector.fluidId, *this->X);

  this->timer->addFluidSolutionTime(t0);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNLNavierStokesFE(DistSVec<double,dim> &U)
{
  if (this->lsMethod == 0)
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);

  computeRKUpdate(U, k1, 1);
  this->multiPhaseSpaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  this->multiPhaseSpaceOp->applyBCsToSolutionVector(U0); //(?)for Navier-Stokes only

  if (this->modifiedGhidaglia) {

    computeRKUpdateHH(U, *hh1);
    *(this->bcData->getBoundaryStateHH()) -= *hh1;
    //std::cout << "Updating mod. ghidaglia" << std::endl;
  }

  U = U0;
  this->checkSolution(U);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNLNavierStokesRK2(DistSVec<double,dim> &U)
{
  if (this->lsMethod == 0)
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);

  computeRKUpdate(U, k1, 1);
  this->multiPhaseSpaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;
//  this->com->barrier();
//  this->com->fprintf(stderr,"Half done.\n");
//  this->com->barrier();
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  this->checkSolution(U0);

  if (this->modifiedGhidaglia) {
    *hhorig = *(this->bcData->getBoundaryStateHH());
    computeRKUpdateHH(U, *hh1);
    *(this->bcData->getBoundaryStateHH()) -= *hh1;
  }

//  // Ghost-Points Population
//  if(this->eqsType == MultiPhysicsTsDesc<dim,dimLS>::NAVIER_STOKES)
//    {
//      this->ghostPoints->deletePointers();
//      this->multiPhaseSpaceOp->populateGhostPoints(this->ghostPoints,*this->X,U,this->varFcn,this->distLSS,this->viscSecOrder, *this->fluidSelector.fluidId);
//    }

  computeRKUpdate(U0, k2, 2);
  this->multiPhaseSpaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  
  if (this->modifiedGhidaglia) {
    computeRKUpdateHH(U0, *hh2);
    *(this->bcData->getBoundaryStateHH()) = (*hhorig)- 0.5*((*hh1)+(*hh2));
  }

  U = U - 1.0/2.0 * (k1 + k2);
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U, Ubc);
  this->multiPhaseSpaceOp->applyBCsToSolutionVector(U);
  this->checkSolution(U);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNLNavierStokesRK4(DistSVec<double,dim> &U)
{ //TODO: no Ghost-Points Population ???

  if (this->lsMethod == 0)
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);

  computeRKUpdate(U, k1, 1);
  this->multiPhaseSpaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;

  if (this->modifiedGhidaglia) {
    *hhorig = *(this->bcData->getBoundaryStateHH());
    computeRKUpdateHH(U, *hh1);
    *(this->bcData->getBoundaryStateHH()) -= 0.5*(*hh1);
  }

  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  this->checkSolution(U0);

  computeRKUpdate(U0, k2, 2);

  if (this->modifiedGhidaglia) {
    computeRKUpdateHH(U0, *hh2);
    *(this->bcData->getBoundaryStateHH()) = *hhorig - 0.5*(*hh2);
  }

  this->multiPhaseSpaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U0 = U - 0.5 * k2;
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  this->checkSolution(U0);

  computeRKUpdate(U0, k3, 3);
  if (this->modifiedGhidaglia) {
    computeRKUpdateHH(U0, *hh3);
    *(this->bcData->getBoundaryStateHH()) = *hhorig - (*hh3);
  }
  this->multiPhaseSpaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U0 = U - k3;
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  this->checkSolution(U0);

  computeRKUpdate(U0, k4, 4);

  if (this->modifiedGhidaglia) {
    computeRKUpdateHH(U0, *hh4);
    *(this->bcData->getBoundaryStateHH()) = *hhorig - 
       1.0/6.0*(*hh1+2.0*(*hh2+*hh3)+*hh4);
  }
  this->multiPhaseSpaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U = U - 1.0/6.0 * (k1 + 2.0 * (k2 + k3) + k4);
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U, Ubc);
  this->multiPhaseSpaceOp->applyBCsToSolutionVector(U);
  this->checkSolution(U);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::computeRKUpdate(DistSVec<double,dim>& Ulocal,
                                                            DistSVec<double,dim>& dU, int it)
{
  DistSVec<double,dimLS>* locphi = &this->Phi;
  if (this->lsMethod == 0)
    locphi = &this->PhiV;
  
  this->multiPhaseSpaceOp->applyBCsToSolutionVector(Ulocal);
  this->multiPhaseSpaceOp->computeResidual(*this->X, *this->A, Ulocal, *this->Wstarij, *this->Wstarji,
                                           this->distLSS, this->linRecAtInterface, this->viscSecOrder, this->riemann, 
                                           this->riemannNormal, *locphi, this->fluidSelector,
                                           dU, it, this->ghostPoints);
                                           //Q: why send PhiV?
                                           //A: Riemann solver needs gradPhi.
                                           //Note: PhiV should be pre-computed. 
  this->timeState->multiplyByTimeStep(dU);

  if(this->numFluid==1&&!this->mmh)
    this->timeState->multiplyByPreconditioner(Ulocal,dU); //low-mach preconditioner
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNLLevelSet(DistSVec<double,dim>& U)
{
  double t0 = this->timer->getTime();
  switch (timeType) {
    case ExplicitData::FORWARD_EULER :
      solveNLLevelSetFE(U); break;
    case ExplicitData::RUNGE_KUTTA_2 :
      solveNLLevelSetRK2(U); break;
    case ExplicitData::RUNGE_KUTTA_4 :
      solveNLLevelSetRK4(U); break;
    default:
      this->com->fprintf(stderr,"ERROR: Choose time-integrator from ForwardEuler, RungeKutta2, and RungeKutta4!\n");
  }
  this->timer->addLevelSetSolutionTime(t0);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNLLevelSetFE(DistSVec<double,dim> &U)
{
  computeRKUpdateLS(this->Phi, *this->fluidSelector.fluidId, p1, U);
  Phi0 = this->Phi;
  this->Phi = this->Phi - p1;
  if(this->withCracking && this->withMixedLS)
    this->riemann->avoidNewPhaseCreation(this->Phi, Phi0);
  else
    this->riemann->avoidNewPhaseCreation(this->Phi, Phi0, this->distLSS);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNLLevelSetRK2(DistSVec<double,dim> &U)
{
  computeRKUpdateLS(this->Phi, *this->fluidSelector.fluidId, p1, U);
  Phi0 = this->Phi - p1;
/* Kevin thinks that it's not necessary to get an updated fluidId0.
  OLD
  this->fluidSelector.getFluidId(fluidId0,Phi0,&(this->distLSS->getStatus()));
  computeRKUpdateLS(Phi0, fluidId0, p2, U);
*/
  computeRKUpdateLS(Phi0, *this->fluidSelector.fluidId, p2, U);
  Phi0 = this->Phi;
  this->Phi = this->Phi - 1.0/2.0 * (p1+p2);
  if(this->withCracking && this->withMixedLS)
    this->riemann->avoidNewPhaseCreation(this->Phi, Phi0);
  else
    this->riemann->avoidNewPhaseCreation(this->Phi, Phi0, this->distLSS);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::solveNLLevelSetRK4(DistSVec<double,dim> &U)
{
  computeRKUpdateLS(this->Phi, *this->fluidSelector.fluidId, p1, U);
  Phi0 = this->Phi - 0.5 * p1;
//  this->fluidSelector.getFluidId(fluidId0,Phi0,&(this->distLSS->getStatus()));

  computeRKUpdateLS(Phi0, *this->fluidSelector.fluidId, p2, U);
  Phi0 = this->Phi - 0.5 * p2;
//  this->fluidSelector.getFluidId(fluidId0,Phi0,&(this->distLSS->getStatus()));

  computeRKUpdateLS(Phi0, *this->fluidSelector.fluidId, p3, U);
  Phi0 = this->Phi - p3;
//  this->fluidSelector.getFluidId(fluidId0,Phi0,&(this->distLSS->getStatus()));

  computeRKUpdateLS(Phi0, *this->fluidSelector.fluidId, p4, U);
  Phi0 = this->Phi;
  this->Phi -= 1.0/6.0 * (p1 + 2.0 * (p2 + p3) + p4);
  if(this->withCracking && this->withMixedLS)
    this->riemann->avoidNewPhaseCreation(this->Phi, Phi0);
  else
    this->riemann->avoidNewPhaseCreation(this->Phi, Phi0, this->distLSS);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::computeRKUpdateLS(DistSVec<double,dimLS> &Philocal,
                                            DistVec<int> &localFluidId,
                                            DistSVec<double,dimLS> &dPhi, DistSVec<double,dim> &U)
{
  this->multiPhaseSpaceOp->computeResidualLS(*this->X, *this->A, Philocal, localFluidId, U, dPhi, this->distLSS, this->linRecAtInterface, this->lsMethod);
  this->timeState->multiplyByTimeStep(dPhi);
  this->LS->checkTrueLevelSetUpdate(dPhi);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::updatePhaseChangeFF(DistSVec<double,dim> &U)
{
  if (this->phaseChangeType == 0) {
    this->riemann->updatePhaseChange(this->V0, *this->fluidSelector.fluidId, fluidId0);
  }
  else
    this->multiPhaseSpaceOp->extrapolatePhaseChange(*this->X,*this->A, this->multiFluidInterfaceOrder-1,
						    U, this->V0,
						    *this->Weights,*this->VWeights,
						    NULL, *this->fluidSelector.fluidId,fluidId0,
						    this->limitHigherOrderExtrapolation);
  
  this->varFcn->primitiveToConservative(this->V0,U,this->fluidSelector.fluidId);
  this->checkSolution(U);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitMultiPhysicsTsDesc<dim,dimLS>::
computeRKUpdateHH(DistSVec<double,dim>& Ulocal,
		  DistVec<double>& dHH) {

  //this->varFcn->conservativeToPrimitive(Ulocal,*(this->V),&this->nodeTag);
  //*(this->bcData->getBoundaryStateHH()) = HHlocal;
  this->domain->computeHHBoundaryTermResidual(*this->bcData,Ulocal,dHH, this->varFcn);
  this->timeState->multiplyByTimeStep(dHH);

}
