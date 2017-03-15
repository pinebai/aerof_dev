#include <ImplicitTsDesc.h>

#include <DistTimeState.h>
#include <GeoSource.h>
#include <SpaceOperator.h>
#include <Domain.h>
#include <MatVecProd.h>
#include <KspSolver.h>
#include <NewtonSolver.h>
#include <MultiGridPrec.h>
#include <cstring>

//------------------------------------------------------------------------------

template<int dim>
ImplicitTsDesc<dim>::ImplicitTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  TsDesc<dim>(ioData, geoSource, dom)
{

  tag = 0;

  failSafeNewton = ioData.ts.implicit.newton.failsafe;
  maxItsNewton = ioData.ts.implicit.newton.maxIts;
  epsNewton = ioData.ts.implicit.newton.eps;  
  epsAbsResNewton = ioData.ts.implicit.newton.epsAbsRes;
  epsAbsIncNewton = ioData.ts.implicit.newton.epsAbsInc;
  maxItsLS = ioData.ts.implicit.newton.lineSearch.maxIts;
  contractionLS = ioData.ts.implicit.newton.lineSearch.rho;
  sufficDecreaseLS = ioData.ts.implicit.newton.lineSearch.c1;
  if (strcmp(ioData.ts.implicit.newton.output, "") == 0)
    outputNewton = 0;
  else if (strcmp(ioData.ts.implicit.newton.output, "stdout") == 0)
    outputNewton = stdout;
  else if (strcmp(ioData.ts.implicit.newton.output, "stderr") == 0)
    outputNewton = stderr;
  else {
    outputNewton = fopen(ioData.ts.implicit.newton.output, "w");
    if (!outputNewton) {
      this->com->fprintf(stderr, "*** Error: could not open \'%s\'\n", ioData.ts.implicit.newton.output);
      exit(1);
    }
  }

  this->timeState = new DistTimeState<dim>(ioData, this->spaceOp, this->varFcn, this->domain, this->V);
  ns = new NewtonSolver<ImplicitTsDesc<dim> >(this);

  myIoDataPtr = &ioData;
  
  if (strcmp(ioData.output.rom.krylovVector,"")==0) {
    kspBinaryOutput = NULL;
  } else {
    kspBinaryOutput = new KspBinaryOutput<DistSVec<double,dim> >(this->domain->getCommunicator(), &ioData, this->domain);
  }

  prevOutputState = NULL; 
  stateIncrement = NULL;

}

//------------------------------------------------------------------------------

template<int dim>
ImplicitTsDesc<dim>::~ImplicitTsDesc()
{
  if (tag) { delete tag; tag = 0; }
  if (ns) { delete ns; ns = 0; }
  if (kspBinaryOutput) delete kspBinaryOutput;
  if (prevOutputState) delete prevOutputState;
  if (stateIncrement) delete stateIncrement;
}

//------------------------------------------------------------------------------

template<int dim>
int ImplicitTsDesc<dim>::solveNonLinearSystem(DistSVec<double,dim> &U, const int timeStep)
{

  int its = 0;
  double t0 = this->timer->getTime();

  TsDesc<dim>::setFailSafe(false);

  its = this->ns->solve(U, timeStep);

  this->errorHandler->reduceError();
  this->data->resolveErrors();
  if(this->errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP]) return its;

  if(its<0) { // failSafe
    U = this->timeState->getUn();
    return its;
  }

  if(TsDesc<dim>::timeState->getData().typeIntegrator == ImplicitData::THREE_POINT_BDF &&
               TsDesc<dim>::timeStepCalculation == TsData::ERRORESTIMATION )
    doErrorEstimation(U);

  this->updateBoundaryExternalState();

  this->timer->addFluidSolutionTime(t0);
  return its;

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitTsDesc<dim>::doErrorEstimation(DistSVec<double,dim> &U)
{
  DistSVec<double,dim> *flux = new DistSVec<double,dim>(TsDesc<dim>::domain->getNodeDistInfo());

  if(this->wallRecType==BcsWallData::CONSTANT)
    this->spaceOp->computeResidual(*this->X, *this->A, this->timeState->getUn(), *flux, this->timeState);
  else
    this->spaceOp->computeResidual(this->riemann1, *this->X, *this->A, this->timeState->getUn(), *flux, this->timeState);

  this->timeState->calculateErrorEstiNorm(U, *flux);

  delete flux;

}

//------------------------------------------------------------------------------
// this function evaluates (Aw),t + F(w,x,v)
template<int dim>
void ImplicitTsDesc<dim>::computeFunction(int it, DistSVec<double,dim> &Q, 
					  DistSVec<double,dim> &F)
{
  // XML
  //spaceOp->applyBCsToSolutionVector(Q);

  if(this->wallRecType==BcsWallData::CONSTANT) {
    this->spaceOp->computeResidual(*this->X, *this->A, Q, F, this->timeState);
  } else {
    this->spaceOp->computeResidual(this->riemann1, *this->X, *this->A, Q, F, this->timeState);
  }

  this->timeState->add_dAW_dt(it, *this->geoState, *this->A, Q, F);

  this->spaceOp->applyBCsToResidual(Q, F);

}

//------------------------------------------------------------------------------
template<int dim>
void ImplicitTsDesc<dim>::recomputeFunction(DistSVec<double,dim> &Q,
                                            DistSVec<double,dim> &rhs)
{
  this->spaceOp->recomputeRHS(*this->X, Q, rhs);
}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
KspPrec<neq> *ImplicitTsDesc<dim>::createPreconditioner(PcData &pcdata, Domain *dom)
{

  KspPrec<neq> *_pc = 0;

  if (pcdata.type == PcData::IDENTITY)
    _pc = new IdentityPrec<neq>();
  else if (pcdata.type == PcData::JACOBI)
    _pc = new JacobiPrec<Scalar,neq>(DiagMat<Scalar,neq>::DENSE, dom);
  else if (pcdata.type == PcData::AS || 
	   pcdata.type == PcData::RAS || 
	   pcdata.type == PcData::ASH || 
	   pcdata.type == PcData::AAS ||
           (pcdata.type == PcData::MG && neq < 5))
    _pc = new IluPrec<Scalar,neq>(pcdata, dom);
  else if (pcdata.type == PcData::MG)  {
    // I just need something to compile
    _pc = new MultiGridPrec<Scalar,neq>(dom,*this->geoState, 
                                        *myIoDataPtr);
  }

  return _pc;

}

//------------------------------------------------------------------------------

template<int dim>
template<int neq>
KspSolver<DistSVec<double,neq>, MatVecProd<dim,neq>, KspPrec<neq>, Communicator> *
ImplicitTsDesc<dim>::createKrylovSolver(const DistInfo &info, KspData &kspdata, 
					MatVecProd<dim,neq> *_mvp, KspPrec<neq> *_pc, 
					Communicator *_com)
{

  KspSolver<DistSVec<double,neq>, MatVecProd<dim,neq>, 
    KspPrec<neq>, Communicator> *_ksp = 0;

  if (kspdata.type == KspData::RICHARDSON)
    _ksp = new RichardsonSolver<DistSVec<double,neq>, MatVecProd<dim,neq>, 
      KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  else if (kspdata.type == KspData::CG)
    _ksp = new CgSolver<DistSVec<double,neq>, MatVecProd<dim,neq>, 
      KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  else if (kspdata.type == KspData::GMRES)
    _ksp = new GmresSolver<DistSVec<double,neq>, MatVecProd<dim,neq>, 
      KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  else if (kspdata.type == KspData::GCR)
     _ksp = new GcrSolver<DistSVec<double,neq>, MatVecProd<dim,neq>,
       KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  return _ksp;

}

//------------------------------------------------------------------------------
template<int dim>
int ImplicitTsDesc<dim>::checkFailSafe(DistSVec<double,dim>& U)
{

  if (!failSafeNewton) return 0;

  if (!tag)
    tag = new DistSVec<bool,2>(this->getVecInfo());

  this->domain->checkFailSafe(this->varFcn, U, *tag);
  this->spaceOp->fix(*tag);

  return failSafeNewton;

}

//------------------------------------------------------------------------------
template<int dim>
void ImplicitTsDesc<dim>::resetFixesTag()
{

  this->spaceOp->resetTag();

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitTsDesc<dim>::writeBinaryVectorsToDiskRom(bool lastNewtonIt, int timeStep, int newtonIt, 
                                                        DistSVec<double,dim> *state, DistSVec<double,dim> *residual)
{
  if (myIoDataPtr->output.rom.avgStateIncrements) {
    if (prevOutputState==NULL) {
      prevOutputState = new DistSVec<double,dim>(this->domain->getNodeDistInfo());
      stateIncrement = new DistSVec<double,dim>(this->domain->getNodeDistInfo());
      int status = this->output->writeBinaryVectorsToDiskRom(lastNewtonIt, timeStep, newtonIt, NULL, residual);
      *prevOutputState = *state;
    } else {
      *stateIncrement = *state - *prevOutputState;
      *stateIncrement *= (1.0 / ((double) myIoDataPtr->output.rom.stateOutputFreqTime));
      int status = this->output->writeBinaryVectorsToDiskRom(lastNewtonIt, timeStep, newtonIt, stateIncrement, residual);
      if (status==1 || status==3) { // state vector was written to disk
        *prevOutputState = *state;
      }
    }
  } else {
    int status = this->output->writeBinaryVectorsToDiskRom(lastNewtonIt, timeStep, newtonIt, state, residual); 
  }
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitTsDesc<dim>::calculateSpatialResidual(DistSVec<double,dim> &Q, DistSVec<double,dim> &spatialRes) {

  if(this->wallRecType==BcsWallData::CONSTANT) {
    this->spaceOp->computeResidual(*this->X, *this->A, Q, spatialRes, this->timeState);
  } else {
    this->spaceOp->computeResidual(this->riemann1, *this->X, *this->A, Q, spatialRes, this->timeState);
  }

  this->spaceOp->applyBCsToResidual(Q, spatialRes);

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitTsDesc<dim>::incrementNewtonOutputTag()
{
    ++(*(this->domain->getNewtonTag()));

}


