#include <ExplicitTsDesc.h>

#include <DistTimeState.h>
#include <GeoSource.h>
#include <SpaceOperator.h>
#include <Domain.h>

//------------------------------------------------------------------------------

template<int dim>
ExplicitTsDesc<dim>::ExplicitTsDesc(IoData& ioData, GeoSource& geoSource, Domain* dom) 
  : TsDesc<dim>(ioData, geoSource, dom), k1(this->getVecInfo()), k2(this->getVecInfo()), 
  k3(this->getVecInfo()), k4(this->getVecInfo()), U0(this->getVecInfo()), Ubc(this->getInletVecInfo()),ratioTimesU(this->getVecInfo())
{
  this->timeState = new DistTimeState<dim>(ioData, this->spaceOp, this->varFcn, this->domain, this->V);
  this->mmh = this->createMeshMotionHandler(ioData, geoSource, 0);
  
  if(ioData.ts.expl.type == ExplicitData::RUNGE_KUTTA_4) RK4 = true;
  else {
    RK4 = false;
    if(ioData.ts.expl.type == ExplicitData::FORWARD_EULER) FE = true;
    else FE = false;
  }
}

//------------------------------------------------------------------------------

template<int dim>
ExplicitTsDesc<dim>::~ExplicitTsDesc()
{

}

//------------------------------------------------------------------------------

template<int dim>
int ExplicitTsDesc<dim>::solveNonLinearSystem(DistSVec<double,dim>& U, int)
{

  double t0 = this->timer->getTime();

  if (RK4) computeRKFourthOrder(U);   // Runge-Kutta 4
  else if(FE) computeForwardEuler(U); // Forward Euler
  else computeRKSecondOrder(U);       // Runge-Kutta 2 

  this->updateBoundaryExternalState();

  int ierr = this->checkSolution(U);  // checking solution for negative pressure or density
  if (ierr > 0) exit(1);

  this->timer->addFluidSolutionTime(t0);

  return 1;

}

//------------------------------------------------------------------------------

template<int dim>
void ExplicitTsDesc<dim>::computeRKFourthOrder(DistSVec<double,dim>& U)
{

// 1st step of Runge Kutta
  computeRKUpdate(U, k1);
  this->spaceOp->getExtrapolationValue(U, Ubc, *this->X);

  U0 = U - 0.5 * k1;
  this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);

//2nd step
  computeRKUpdate(U0, k2);
  this->spaceOp->getExtrapolationValue(U0, Ubc, *this->X);

  U0 = U - 0.5 * k2;
  this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  
//3rd step
  computeRKUpdate(U0, k3);
  this->spaceOp->getExtrapolationValue(U0, Ubc, *this->X);

  U0 = U - k3;
  this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);

//4th step
  computeRKUpdate(U0, k4);
  this->spaceOp->getExtrapolationValue(U0, Ubc, *this->X);

  U -= 1.0/6.0 * (k1 + 2.0 * (k2 + k3) + k4);
  this->spaceOp->applyExtrapolationToSolutionVector(U, Ubc);


  this->spaceOp->applyBCsToSolutionVector(U);

}

//------------------------------------------------------------------------------

template<int dim>
void ExplicitTsDesc<dim>::computeRKSecondOrder(DistSVec<double,dim>& U)
{
  this->domain->computePrdtWCtrlVolRatio(ratioTimesU, U, *this->A, *this->geoState);

  computeRKUpdate(U, k1);
 
  this->spaceOp->getExtrapolationValue(U, Ubc, *this->X);

  U0 = ratioTimesU - k1;
  this->spaceOp->applyExtrapolationToSolutionVector(U0, Ubc);
  computeRKUpdate(U0, k2);
  this->spaceOp->getExtrapolationValue(U0, Ubc, *this->X);

  U = ratioTimesU - 1.0/2.0 * (k1 + k2);

  this->spaceOp->applyExtrapolationToSolutionVector(U, Ubc);
  this->spaceOp->applyBCsToSolutionVector(U);

}

//------------------------------------------------------------------------------

template<int dim>
void ExplicitTsDesc<dim>::computeForwardEuler(DistSVec<double,dim>& U)
{
  this->domain->computePrdtWCtrlVolRatio(ratioTimesU, U, *this->A, *this->geoState);

  computeRKUpdate(U, k1);

  this->spaceOp->getExtrapolationValue(U, Ubc, *this->X);
  U = ratioTimesU - k1;
  this->spaceOp->applyExtrapolationToSolutionVector(U, Ubc);
  this->spaceOp->applyBCsToSolutionVector(U);

}

//------------------------------------------------------------------------------

template<int dim>
void ExplicitTsDesc<dim>::computeRKUpdate(DistSVec<double,dim>& U,
				DistSVec<double,dim>& dU)
{

  this->spaceOp->applyBCsToSolutionVector(U);

  if(this->wallRecType==BcsWallData::CONSTANT) 
    this->spaceOp->computeResidual(*this->X, *this->A, U, dU, this->timeState);
  else
    this->spaceOp->computeResidual(this->riemann1, *this->X, *this->A, U, dU, this->timeState);

  this->timeState->multiplyByTimeStep(dU);
  this->timeState->multiplyByPreconditioner(U,dU);
}

//------------------------------------------------------------------------------
