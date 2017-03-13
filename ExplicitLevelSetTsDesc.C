#include <ExplicitLevelSetTsDesc.h>

#include "IoData.h"
#include "Domain.h"
#include <GeoSource.h>
#include <DistTimeState.h>
#include <SpaceOperator.h>
#include "LevelSet.h"

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
ExplicitLevelSetTsDesc<dim,dimLS>::
ExplicitLevelSetTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  LevelSetTsDesc<dim,dimLS>(ioData, geoSource, dom),
  k1(this->getVecInfo()), k2(this->getVecInfo()), 
  k3(this->getVecInfo()), k4(this->getVecInfo()), 
  p1(this->getVecInfo()), p2(this->getVecInfo()), 
  p3(this->getVecInfo()), p4(this->getVecInfo()), 
  U0(this->getVecInfo()), Phi0(this->getVecInfo()),
  ratioTimesU(this->getVecInfo()), 
  ratioTimesPhi(this->getVecInfo()),
  fluidId0(this->getVecInfo())
{
  this->mmh = this->createMeshMotionHandler(ioData, geoSource, 0);

  timeType = ioData.ts.expl.type;

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
ExplicitLevelSetTsDesc<dim,dimLS>::~ExplicitLevelSetTsDesc()
{
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
int ExplicitLevelSetTsDesc<dim,dimLS>::solveNonLinearSystem(DistSVec<double,dim> &U, int)
{

  if(timeType == ExplicitData::FORWARD_EULER || 
     timeType == ExplicitData::ONE_BLOCK_RK2 ||
     timeType == ExplicitData::ONE_BLOCK_RK2bis )
    solveNLSystemOneBlock(U);
  else if(timeType == ExplicitData::RUNGE_KUTTA_2 ||
          timeType == ExplicitData::RUNGE_KUTTA_4 )
    solveNLSystemTwoBlocks(U);

  return 1;

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLSystemOneBlock(DistSVec<double,dim> &U)
{

  if(timeType == ExplicitData::FORWARD_EULER)
    solveNLAllFE(U);
  else if(timeType == ExplicitData::ONE_BLOCK_RK2)
    solveNLAllRK2(U);
  else if(timeType == ExplicitData::ONE_BLOCK_RK2bis)
    solveNLAllRK2bis(U);

}
//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLAllFE(DistSVec<double,dim> &U)
{

  double t0 = this->timer->getTime();

  DistSVec<double,dim> Ubc(this->getVecInfo());
  if (this->lsMethod == 0) {
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);
  }

  computeRKUpdate(U, k1,1);

  if (this->modifiedGhidaglia) {

    computeRKUpdateHH(U, *hh1);
    *(this->bcData->getBoundaryStateHH()) -= *hh1;
  }

  this->multiPhaseSpaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);

  this->checkSolution(U);
  this->errorHandler->reduceError();
  this->data->resolveErrors();
  if(this->errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP]) return;

  this->timer->addFluidSolutionTime(t0);

  if(!(this->interfaceType==MultiFluidData::FSF)){
    this->varFcn->conservativeToPrimitive(U0,this->V0,this->fluidSelector.fluidId);
    this->riemann->storePreviousPrimitive(this->V0, *this->fluidSelector.fluidId, *this->X);
    t0 = this->timer->getTime();

    computeRKUpdateLS(this->Phi, *this->fluidSelector.fluidId, p1, U);

    this->Phi = this->Phi - p1;
    this->riemann->avoidNewPhaseCreation(this->Phi, this->LS->Phin);
    
    // Riemann overwrite using the value of Phi_{n+1} 
    //
    this->setPhiExact();
     
    if (this->myTriangulatedInterface)
      this->myTriangulatedInterface->update(this->currentTimeStep);

    if (this->lsMethod == 0 ||
        this->lsMethod == 1)
      (this->fluidSelector).getFluidId(this->Phi); //update fluidId accordingly
    else 
      (this->fluidSelector).getFluidId(this->myTriangulatedInterface);

    this->timer->addLevelSetSolutionTime(t0);

    if (this->phaseChangeType == 0)
      this->riemann->updatePhaseChange(this->V0, *this->fluidSelector.fluidId, *this->fluidSelector.fluidIdn);
    else {
      this->multiPhaseSpaceOp->setLastPhaseChangeValues(this->riemann);
      this->multiPhaseSpaceOp->extrapolatePhaseChange(*this->X,*this->A, this->interfaceOrder-1,
						      U0, this->V0,
						      this->Weights,this->VWeights,
						      NULL, *this->fluidSelector.fluidId, *this->fluidSelector.fluidIdn,
                                                      this->limitHigherOrderExtrapolation);

    }
    this->varFcn->primitiveToConservative(this->V0,U0,this->fluidSelector.fluidId);
   
  }

  U = U0;

  this->checkSolution(U);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLAllRK2(DistSVec<double,dim> &U)
{

  double t0 = this->timer->getTime();
  this->domain->computePrdtWCtrlVolRatio(ratioTimesU, U, *this->A, *this->geoState);
  this->domain->computePrdtPhiCtrlVolRatio(ratioTimesPhi, this->Phi, *this->A, *this->geoState);

  DistSVec<double,dim> Ubc(this->getVecInfo());
  if (this->lsMethod == 0) {
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);
  }
  // *** prediction step ***
  computeRKUpdate(U, k1,1);

  if (this->modifiedGhidaglia) {
    *hhorig = *(this->bcData->getBoundaryStateHH());
    computeRKUpdateHH(U, *hh1);
    *(this->bcData->getBoundaryStateHH()) -= *hh1;
  }

  this->multiPhaseSpaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = ratioTimesU - k1;
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);

  this->checkSolution(U);
  this->errorHandler->reduceError();
  this->data->resolveErrors();
  if(this->errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP]) return;

  this->timer->addFluidSolutionTime(t0);

  if(!(this->interfaceType==MultiFluidData::FSF)){
    this->varFcn->conservativeToPrimitive(U0,this->V0,this->fluidSelector.fluidId);
    this->riemann->storePreviousPrimitive(this->V0, *this->fluidSelector.fluidId, *this->X);

    t0 = this->timer->getTime();

    /*computeRKUpdateLS(this->Phi, *this->fluidSelector.fluidId, p1, U);
    Phi0 = ratioTimesPhi - p1;
    this->fluidSelector.getFluidId(fluidId0,Phi0);
*/

    this->timer->addLevelSetSolutionTime(t0);

    // Riemann overwrite on U0 using the value of Phi0
  }

  t0 = this->timer->getTime();

  // *** corrector step ***
  computeRKUpdate(U0, k2,2);
  
  if (this->modifiedGhidaglia) {
    computeRKUpdateHH(U0, *hh2);
    *(this->bcData->getBoundaryStateHH()) = (*hhorig)- 0.5*((*hh1)+(*hh2));
  }

  this->multiPhaseSpaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U = ratioTimesU - 0.5 * (k1 + k2);
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U, Ubc);
  this->multiPhaseSpaceOp->applyBCsToSolutionVector(U);

  this->checkSolution(U);
  this->errorHandler->reduceError();
  this->data->resolveErrors();
  if(this->errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP]) return;

  this->timer->addFluidSolutionTime(t0);

  if(!(this->interfaceType==MultiFluidData::FSF)){

    t0 = this->timer->getTime();

    computeRKUpdateLS(Phi0, fluidId0, p2, U0);
    this->Phi = ratioTimesPhi - 0.5 * (p1 + p2);
    this->riemann->avoidNewPhaseCreation(this->Phi, this->LS->Phin);
    // Riemann overwrite on U_{n+1} using the value of Phi_{n+1}
    this->setPhiExact();
    if (this->myTriangulatedInterface)
      this->myTriangulatedInterface->update(this->currentTimeStep);

    if (this->lsMethod == 0 ||
        this->lsMethod == 1)
      (this->fluidSelector).getFluidId(this->Phi); //update fluidId accordingly
    else 
      (this->fluidSelector).getFluidId(this->myTriangulatedInterface);


    this->timer->addLevelSetSolutionTime(t0);

    if (this->phaseChangeType == 0)
      this->riemann->updatePhaseChange(this->V0, *this->fluidSelector.fluidId, *this->fluidSelector.fluidIdn);
    else
      this->multiPhaseSpaceOp->extrapolatePhaseChange(*this->X,*this->A, this->interfaceOrder-1,
						      U, this->V0,
						      this->Weights,this->VWeights,
						      NULL, *this->fluidSelector.fluidId, *this->fluidSelector.fluidIdn,
                                                      this->limitHigherOrderExtrapolation);

 
    this->varFcn->primitiveToConservative(this->V0,U,this->fluidSelector.fluidId);
  }
  this->checkSolution(U);

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLAllRK2bis(DistSVec<double,dim> &U)
{

  double t0 = this->timer->getTime();

  DistSVec<double,dim> Ubc(this->getVecInfo());
  if (this->lsMethod == 0) {
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);
  }
  // *** prediction step ***
  computeRKUpdate(U, k1,1);
  this->multiPhaseSpaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - k1;
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);

  this->checkSolution(U);
  this->errorHandler->reduceError();
  this->data->resolveErrors();
  if(this->errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP]) return;

  this->timer->addFluidSolutionTime(t0);

  if(!(this->interfaceType==MultiFluidData::FSF)){
    this->varFcn->conservativeToPrimitive(U0,this->V0,this->fluidSelector.fluidId);
    this->riemann->storePreviousPrimitive(this->V0, *this->fluidSelector.fluidId, *this->X);

    t0 = this->timer->getTime();

    computeRKUpdateLS(this->Phi, *this->fluidSelector.fluidId, p1, U0);
    Phi0 = this->Phi - p1;
    
    this->fluidSelector.getFluidId(fluidId0,Phi0);

    this->timer->addLevelSetSolutionTime(t0);

    // Riemann overwrite on U0 using the value of Phi0
    //this->multiPhaseSpaceOp->updatePhaseChange(this->Vg, U0, fluidId0,
    //                                 this->fluidSelector.fluidIdn, this->Vgf,
    //                                 this->Vgfweight, this->riemann);
    //this->riemann->updatePhaseChange();
  }

  t0 = this->timer->getTime();

  // *** corrector step ***
  computeRKUpdate(U0, k2,2);
  this->multiPhaseSpaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U -= 0.5 * (k1 + k2);
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U, Ubc);
  this->multiPhaseSpaceOp->applyBCsToSolutionVector(U);

  this->checkSolution(U);
  this->errorHandler->reduceError();
  this->data->resolveErrors();
  if(this->errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP]) return;

  this->timer->addFluidSolutionTime(t0);

  if(!(this->interfaceType==MultiFluidData::FSF)){
    //this->multiPhaseSpaceOp->storePreviousPrimitive(U, this->Vg, this->Phi,
    //                                      this->Vgf, this->Vgfweight);
    //this->riemann->storePreviousPrimitive();

    t0 = this->timer->getTime();

    computeRKUpdateLS(Phi0, fluidId0, p2, U);
    this->Phi = this->Phi - 0.5 * (p1 + p2);
    this->riemann->avoidNewPhaseCreation(this->Phi, this->LS->Phin);
    (this->fluidSelector).getFluidId(this->Phi);

    this->timer->addLevelSetSolutionTime(t0);

    // Riemann overwrite on U_{n+1} using the value of Phi_{n+1}
    this->setPhiExact();
    if (this->phaseChangeType == 0)
      this->riemann->updatePhaseChange(this->V0, *this->fluidSelector.fluidId, *this->fluidSelector.fluidIdn);
    else
      this->multiPhaseSpaceOp->extrapolatePhaseChange(*this->X,*this->A, this->interfaceOrder-1,
						      U, this->V0,
						      this->Weights,this->VWeights,
						      NULL, *this->fluidSelector.fluidId, *this->fluidSelector.fluidIdn,
                                                      this->limitHigherOrderExtrapolation);

    this->varFcn->primitiveToConservative(this->V0,U,this->fluidSelector.fluidId);
  }
  this->checkSolution(U);

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
int ExplicitLevelSetTsDesc<dim,dimLS>::solveNLSystemTwoBlocks(DistSVec<double,dim> &U)
{

  solveNLEuler(U);

  if(!(this->interfaceType==MultiFluidData::FSF)){

    this->varFcn->conservativeToPrimitive(U,this->V0,this->fluidSelector.fluidId);
    this->riemann->storePreviousPrimitive(this->V0, *this->fluidSelector.fluidId, *this->X);

    if (this->myExactInterface.size() == 0)
      solveNLLevelSet(U);
    else {
      this->setPhiExact();
    }
   
    if (this->myTriangulatedInterface) {
      this->varFcn->conservativeToPrimitive(U,this->V0,this->fluidSelector.fluidId);
      if (this->myExactInterface.size() == 0) {
        DistNodalGrad<dim, double> * G = this->multiPhaseSpaceOp->getDistNodalGrad(U);
        G->compute(this->geoState->getConfig(), *(this->X), 
                   *this->A, *this->fluidSelector.fluidId, (this->V0));
        this->myTriangulatedInterface->integraterk2(this->domain, *this->X,
                                                    this->V0,*G,
                                                    *this->fluidSelector.fluidId,
                                                    this->currentTimeStep);
      }
      this->myTriangulatedInterface->update(this->currentTimeStep);
    }

    if (this->lsMethod == 0 ||
      this->lsMethod == 1)
      (this->fluidSelector).getFluidId(this->Phi); //update fluidId accordingly
    else 
      (this->fluidSelector).getFluidId(this->myTriangulatedInterface);
   
    if (this->phaseChangeType == 0)
      this->riemann->updatePhaseChange(this->V0, *this->fluidSelector.fluidId, *this->fluidSelector.fluidIdn);
    else {
      this->multiPhaseSpaceOp->setLastPhaseChangeValues(this->riemann);
      this->multiPhaseSpaceOp->extrapolatePhaseChange(*this->X, *this->A,this->interfaceOrder-1,
						      U, this->V0,
						      this->Weights,this->VWeights,
						      NULL, *this->fluidSelector.fluidId, *this->fluidSelector.fluidIdn,
                                                      this->limitHigherOrderExtrapolation);
    }

    this->varFcn->primitiveToConservative(this->V0,U,this->fluidSelector.fluidId);
  }

  return 0;
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLEuler(DistSVec<double,dim> &U)
{
  double t0 = this->timer->getTime();

  if (timeType == ExplicitData::RUNGE_KUTTA_4) solveNLEulerRK4(U);
  else                                         solveNLEulerRK2(U);

  this->checkSolution(U);
  this->errorHandler->reduceError();
  this->data->resolveErrors();
  if(this->errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP]) return;
 
  this->timer->addFluidSolutionTime(t0);

}
//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLEulerRK2(DistSVec<double,dim> &U)
{
  this->domain->computePrdtWCtrlVolRatio(ratioTimesU, U, *this->A, *this->geoState);
  DistSVec<double,dim> Ubc(this->getVecInfo());
  if (this->lsMethod == 0) {
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);
  }

  computeRKUpdate(U, k1,1);

  if (this->modifiedGhidaglia) {
    *hhorig = *(this->bcData->getBoundaryStateHH());
    computeRKUpdateHH(U, *hh1);
    *(this->bcData->getBoundaryStateHH()) -= *hh1;
  }

  this->multiPhaseSpaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = ratioTimesU - k1;

  this->domain->setExactBoundaryValues(U0, *this->X, this->ioData, 
				       this->currentTime + this->currentTimeStep,
				       this->spaceOp->getVarFcn());

  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  this->checkSolution(U0);
    
  /*if (this->lsMethod == 1) {
    computeRKUpdateLS(this->Phi, *this->fluidSelector.fluidId, p1, U);
    Phi0 = this->Phi - 0.5*p1;
  }*/

  computeRKUpdate(U0, k2,2);
  
  if (this->modifiedGhidaglia) {
    computeRKUpdateHH(U0, *hh2);
    *(this->bcData->getBoundaryStateHH()) = (*hhorig)- 0.5*((*hh1)+(*hh2));
  }


  this->multiPhaseSpaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U = ratioTimesU - 1.0/2.0 * (k1 + k2);

  this->domain->setExactBoundaryValues(U, *this->X, this->ioData, 
				       this->currentTime + this->currentTimeStep,
				       this->spaceOp->getVarFcn());

  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U, Ubc);
  this->multiPhaseSpaceOp->applyBCsToSolutionVector(U);
  this->checkSolution(U);

  this->updateBoundaryExternalState();
}
//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLEulerRK4(DistSVec<double,dim> &U)
{

  DistSVec<double,dim> Ubc(this->getVecInfo());
  if (this->lsMethod == 0) {
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);
  }

  computeRKUpdate(U, k1, 1);
 
  if (this->modifiedGhidaglia) {
    *hhorig = *(this->bcData->getBoundaryStateHH());
    computeRKUpdateHH(U, *hh1);
    *(this->bcData->getBoundaryStateHH()) -= 0.5*(*hh1);
  }


  this->multiPhaseSpaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U0 = U - 0.5 * k1;
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  this->checkSolution(this->U0);

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

  computeRKUpdate(U0, k4, 1);

  if (this->modifiedGhidaglia) {
    computeRKUpdateHH(U0, *hh4);
    *(this->bcData->getBoundaryStateHH()) = *hhorig - 
       1.0/6.0*(*hh1+2.0*(*hh2+*hh3)+*hh4);
  }


  this->multiPhaseSpaceOp->getExtrapolationValue(U0, Ubc, *this->X);
  U -= 1.0/6.0 * (k1 + 2.0 * (k2 + k3) + k4);
  this->multiPhaseSpaceOp->applyExtrapolationToSolutionVector(U, Ubc);
  this->multiPhaseSpaceOp->applyBCsToSolutionVector(U);
  this->checkSolution(U);

  this->updateBoundaryExternalState();
}

//------------------------------------------------------------------------------
template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLLevelSet(DistSVec<double,dim> &U)
{

  double t0 = this->timer->getTime();

  if (timeType == ExplicitData::RUNGE_KUTTA_4) solveNLLevelSetRK4(U);
  else                                         solveNLLevelSetRK2(U);

  this->timer->addLevelSetSolutionTime(t0);

}
//------------------------------------------------------------------------------
template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLLevelSetRK2(DistSVec<double,dim> &U)
{

  this->domain->computePrdtPhiCtrlVolRatio(ratioTimesPhi, this->Phi, *this->A, *this->geoState);
  computeRKUpdateLS(this->Phi, *this->fluidSelector.fluidId, p1, U);
  Phi0 = ratioTimesPhi - p1;
  this->riemann->avoidNewPhaseCreation(this->Phi0, this->LS->Phin);
  this->fluidSelector.getFluidId(fluidId0,Phi0);

  computeRKUpdateLS(Phi0, fluidId0, p2, U);
  this->Phi = ratioTimesPhi - 1.0/2.0 * (p1+p2);
  this->riemann->avoidNewPhaseCreation(this->Phi, this->LS->Phin);
  // Riemann overwrite on U_{n+1} using the value of Phi_{n+1}
  this->setPhiExact();
}
//------------------------------------------------------------------------------
template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::solveNLLevelSetRK4(DistSVec<double,dim> &U)
{
  computeRKUpdateLS(this->Phi, *this->fluidSelector.fluidId, p1, U);
  Phi0 = this->Phi - 0.5 * p1;
  this->riemann->avoidNewPhaseCreation(this->Phi0, this->LS->Phin);
  this->fluidSelector.getFluidId(fluidId0,Phi0);

  computeRKUpdateLS(Phi0, fluidId0, p2, U);
  Phi0 = this->Phi - 0.5 * p2;
  this->riemann->avoidNewPhaseCreation(this->Phi0, this->LS->Phin);
  this->fluidSelector.getFluidId(fluidId0,Phi0);

  computeRKUpdateLS(Phi0, fluidId0, p3, U);
  Phi0 = this->Phi - p3;
  this->riemann->avoidNewPhaseCreation(this->Phi0, this->LS->Phin);
  this->fluidSelector.getFluidId(fluidId0,Phi0);

  computeRKUpdateLS(Phi0, fluidId0, p4, U);
  this->Phi -= 1.0/6.0 * (p1 + 2.0 * (p2 + p3) + p4);
  this->riemann->avoidNewPhaseCreation(this->Phi, this->LS->Phin);
  // Riemann overwrite on U_{n+1} using the value of Phi_{n+1}
  this->setPhiExact();
}
//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::computeRKUpdate(DistSVec<double,dim>& Ulocal,
							DistSVec<double,dim>& dU, int it)
{
  this->multiPhaseSpaceOp->applyBCsToSolutionVector(Ulocal);

  DistSVec<double,dimLS>* locphi = &this->Phi;
  if (this->lsMethod == 0)
    locphi = &this->PhiV;
  /*else if (it > 1 &&  timeType == ExplicitData::RUNGE_KUTTA_2)
    locphi = this->Phi0;
*/

  if (this->interfaceOrder == 2) {
    this->varFcn->conservativeToPrimitive(Ulocal,this->V0,this->fluidSelector.fluidId);
    this->varFcn->primitiveToConservative(this->V0,Ulocal,this->fluidSelector.fluidId);
  }

  this->multiPhaseSpaceOp->computeResidual(*this->X, *this->A, Ulocal, *locphi, this->fluidSelector, 
					   dU, this->riemann,this->timeState,it);
                                 //Q: why send PhiV?
                                 //A: Riemann solver needs gradPhi.
  // for RK2 on moving grids
  this->timeState->multiplyByTimeStep(dU);
}

//------------------------------------------------------------------------------
template<int dim, int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::computeRKUpdateLS(DistSVec<double,dimLS> &Philocal,
                                  DistVec<int> &localFluidId,
                                  DistSVec<double,dimLS> &dPhi, DistSVec<double,dim> &U)
{

  this->multiPhaseSpaceOp->computeResidualLS(*this->X, *this->A, Philocal, localFluidId, U, dPhi,0,0,this->lsMethod, this->interfaceOrder);
  // for RK2 on moving grids
  this->timeState->multiplyByTimeStep(dPhi);
  this->LS->checkTrueLevelSetUpdate(dPhi);

}

template<int dim,int dimLS>
void ExplicitLevelSetTsDesc<dim,dimLS>::computeRKUpdateHH(DistSVec<double,dim>& Ulocal,
				 		          DistVec<double>& dHH) {

  //this->varFcn->conservativeToPrimitive(Ulocal,*(this->V),&this->nodeTag);
  //*(this->bcData->getBoundaryStateHH()) = HHlocal;
  this->domain->computeHHBoundaryTermResidual(*this->bcData,Ulocal,dHH, this->varFcn);
  this->timeState->multiplyByTimeStep(dHH);

}

