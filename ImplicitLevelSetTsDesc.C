#include <ImplicitLevelSetTsDesc.h>

#include "IoData.h"
#include "GeoSource.h"
#include "Domain.h"
#include "LevelSet.h"
#include "DistTimeState.h"

#include <MatVecProd.h>
#include <KspSolver.h>
#include <SpaceOperator.h>
#include <NewtonSolver.h>


#ifdef TYPE_PREC
#define PrecScalar TYPE_PREC
#else
#define PrecScalar double
#endif

//------------------------------------------------------------------------------

template<int dim, int dimLS>
ImplicitLevelSetTsDesc<dim,dimLS>::
ImplicitLevelSetTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  LevelSetTsDesc<dim,dimLS>(ioData, geoSource, dom),U0(this->getVecInfo()),Fold(this->getVecInfo()), embeddedU(dom->getNodeDistInfo()), embeddedB(dom->getNodeDistInfo()),embeddeddQ(dom->getNodeDistInfo())

{
  tag = 0;
  ImplicitData &implicitData = ioData.ts.implicit;
  
  // NewtonSolver
  ns = new NewtonSolver<ImplicitLevelSetTsDesc<dim,dimLS> >(this);
  failSafeNewton = implicitData.newton.failsafe;
  maxItsNewton = implicitData.newton.maxIts;
  epsNewton = implicitData.newton.eps;
  epsAbsResNewton = implicitData.newton.epsAbsRes;
  epsAbsIncNewton = implicitData.newton.epsAbsInc;
  maxItsLS = implicitData.newton.lineSearch.maxIts;
  contractionLS = implicitData.newton.lineSearch.rho;
  sufficDecreaseLS = implicitData.newton.lineSearch.c1;
  if (strcmp(implicitData.newton.output, "") == 0)
    outputNewton = 0;
  else if (strcmp(implicitData.newton.output, "stdout") == 0)
    outputNewton = stdout;
  else if (strcmp(implicitData.newton.output, "stderr") == 0)
    outputNewton = stderr;
  else {
    outputNewton = fopen(implicitData.newton.output, "w");
    if (!outputNewton) {
      this->com->fprintf(stderr, "*** Error: could not open \'%s\'\n", implicitData.newton.output);
      exit(1);
    }
  }

  // MatVecProd, Prec and Krylov solver for Euler equations
  if (implicitData.mvp == ImplicitData::FD)
    mvp = new MatVecProdFDMultiPhase<dim,dimLS>(this->timeState, this->geoState,
                                                this->multiPhaseSpaceOp, this->riemann,
                                                &this->fluidSelector, this->domain,
						ioData);
  else if (implicitData.mvp == ImplicitData::H1)
    mvp = new MatVecProdH1MultiPhase<dim,dimLS>(this->timeState, this->multiPhaseSpaceOp, this->riemann, &this->fluidSelector, this->domain);
  else{
    this->com->fprintf(stdout, "*** Error: MatVecProd H2 for MultiPhase is not available\n");
    exit(1);
  }
  
  pc = createPreconditioner<PrecScalar,dim>(implicitData.newton.ksp.ns.pc, this->domain);
  
  ksp = createKrylovSolver(this->getVecInfo(), implicitData.newton.ksp.ns, mvp, pc, this->com);
  

  // MatVecProd, Prec and Krylov solver for LevelSet equation
  mvpLS  = new MatVecProdLS<dim,dimLS>(this->timeState, this->geoState, this->multiPhaseSpaceOp, this->domain, this->LS);

  pcLS = createPreconditioner<PrecScalar,dimLS>(implicitData.newton.ksp.lsi.pc, this->domain);
  kspLS = createKrylovSolverLS(this->getVecInfo(), implicitData.newton.ksp.lsi, mvpLS, pcLS, this->com);


  // meshmotion
  MemoryPool mp;
  
  mvp->exportMemory(&mp);
  pc->exportMemory(&mp);
  
  this->mmh = this->createMeshMotionHandler(ioData, geoSource, &mp);

  if (this->modifiedGhidaglia) {
    embeddedU.addHHBoundaryTerm(dom->getFaceDistInfo());
    embeddeddQ.addHHBoundaryTerm(dom->getFaceDistInfo());
    embeddedB.addHHBoundaryTerm(dom->getFaceDistInfo());
    
    hhResidual = new DistVec<double>(dom->getFaceDistInfo());
  }
 
  if (this->modifiedGhidaglia)
    mvp->attachHH(this->embeddedU);

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
ImplicitLevelSetTsDesc<dim,dimLS>::~ImplicitLevelSetTsDesc()
{
  if (tag)   delete tag;
  if (mvp)   delete mvp;
  if (pc)    delete pc;
  if (ksp)   delete ksp;
  if (mvpLS) delete mvpLS;
  if (kspLS) delete kspLS;
  if (pcLS)  delete pcLS;
  if (ns)    delete ns;

}

//------------------------------------------------------------------------------
//  Internal routines to setup the class (called in constructor)
//------------------------------------------------------------------------------

template<int dim, int dimLS>
template<class Scalar, int neq>
KspPrec<neq> *ImplicitLevelSetTsDesc<dim,dimLS>::createPreconditioner(PcData &pcdata, Domain *dom)
{
  
  KspPrec<neq> *_pc = 0;
  
  if (pcdata.type == PcData::IDENTITY)
    _pc = new IdentityPrec<neq>();
  else if (pcdata.type == PcData::JACOBI)
    _pc = new JacobiPrec<Scalar,neq>(DiagMat<Scalar,neq>::DENSE, dom);
  else if (pcdata.type == PcData::AS ||
	   pcdata.type == PcData::RAS ||
	   pcdata.type == PcData::ASH ||
	   pcdata.type == PcData::AAS)
    _pc = new IluPrec<Scalar,neq>(pcdata, dom);
  
  return _pc;
  
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
template<int neq, class MatVecProdOp>
KspSolver<DistEmbeddedVec<double,neq>, MatVecProdOp, KspPrec<neq>, Communicator> *
ImplicitLevelSetTsDesc<dim,dimLS>::createKrylovSolver(
                               const DistInfo &info, KspData &kspdata,
                               MatVecProdOp *_mvp, KspPrec<neq> *_pc,
                               Communicator *_com)
{
  
  KspSolver<DistEmbeddedVec<double,neq>, MatVecProdOp, KspPrec<neq>, Communicator> *_ksp = 0;
  
  if (kspdata.type == KspData::RICHARDSON)
    _ksp = new RichardsonSolver<DistEmbeddedVec<double,neq>, MatVecProdOp,
                 KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  else if (kspdata.type == KspData::CG)
    _ksp = new CgSolver<DistEmbeddedVec<double,neq>, MatVecProdOp,
                 KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  else if (kspdata.type == KspData::GMRES)
    _ksp = new GmresSolver<DistEmbeddedVec<double,neq>, MatVecProdOp,
                 KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  
  return _ksp;
  
}

template<int dim, int dimLS>
template<int neq, class MatVecProdOp>
KspSolver<DistSVec<double,neq>, MatVecProdOp, KspPrec<neq>, Communicator> *
ImplicitLevelSetTsDesc<dim,dimLS>::createKrylovSolverLS(
                               const DistInfo &info, KspData &kspdata,
                               MatVecProdOp *_mvp, KspPrec<neq> *_pc,
                               Communicator *_com)
{
  
  KspSolver<DistSVec<double,neq>, MatVecProdOp, KspPrec<neq>, Communicator> *_ksp = 0;
  
  if (kspdata.type == KspData::RICHARDSON)
    _ksp = new RichardsonSolver<DistSVec<double,neq>, MatVecProdOp,
                 KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  else if (kspdata.type == KspData::CG)
    _ksp = new CgSolver<DistSVec<double,neq>, MatVecProdOp,
                 KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  else if (kspdata.type == KspData::GMRES)
    _ksp = new GmresSolver<DistSVec<double,neq>, MatVecProdOp,
                 KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  
  return _ksp;
  
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// External routine to solve problem (called by TsSolver)
// It calls for the NewtonSolver ns, which in turn will
// call routines below from this same file or from LevelSetTsDesc
//------------------------------------------------------------------------------
template<int dim, int dimLS>
int ImplicitLevelSetTsDesc<dim,dimLS>::solveNonLinearSystem(DistSVec<double,dim> &U, int)
{
  
  int its;

  double t0 = this->timer->getTime();

  if (this->timeState->useNm1() && this->timeState->existsNm1() &&
      this->phaseChangeType == 1) {
    DistSVec<double,dim>& Unm1 = this->timeState->getUnm1();
  
    this->varFcn->conservativeToPrimitive(U,this->V0,this->fluidSelector.fluidId);
    this->multiPhaseSpaceOp->extrapolatePhaseChange(*this->X, *this->A, this->interfaceOrder-1,
 					            Unm1, this->V0,
						    this->Weights,this->VWeights,
						    NULL, *this->fluidSelector.fluidId, *this->fluidSelector.fluidIdnm1,
						    false);
    this->varFcn->primitiveToConservative(this->V0,Unm1,this->fluidSelector.fluidId);
  }
    
 
/*  if (this->interfaceOrder == 2)
    int itsLS = this->ns->solveLS(this->Phi, U);
*/

  if (this->modifiedGhidaglia)
    embeddedU.hh() = *this->bcData->getBoundaryStateHH();
  
  this->domain->setExactBoundaryValues(U, *this->X, this->ioData, 
				       this->currentTime + this->currentTimeStep,
				       this->spaceOp->getVarFcn());
  

  its = this->ns->solve(U);

  this->errorHandler->reduceError();
  this->data->resolveErrors();
  if(this->errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP]) return its;

  this->timer->addFluidSolutionTime(t0);
  this->Utilde = U;

  if(TsDesc<dim>::timeState->getData().typeIntegrator == ImplicitData::THREE_POINT_BDF &&
     TsDesc<dim>::timeStepCalculation == TsData::ERRORESTIMATION )
    doErrorEstimation(U);
  
  if(!(this->interfaceType==MultiFluidData::FSF)){
    this->varFcn->conservativeToPrimitive(U,this->V0,this->fluidSelector.fluidId);
    this->riemann->storePreviousPrimitive(this->V0, *this->fluidSelector.fluidId, *this->X);
    
    double t1 = this->timer->getTime();
    int itsLS = this->ns->solveLS(this->Phi, U);
    this->riemann->storeOldV(U);
    this->riemann->avoidNewPhaseCreation(this->Phi, this->LS->Phin);
    (this->fluidSelector).getFluidId(this->Phi);
    this->timer->addLevelSetSolutionTime(t1);
    
    // Riemann overwrite using the value of Phi_{n+1}
    this->setPhiExact();

    if (this->phaseChangeType == 0)
      this->riemann->updatePhaseChange(this->V0, *this->fluidSelector.fluidId, *this->fluidSelector.fluidIdn);
    else
      this->multiPhaseSpaceOp->extrapolatePhaseChange(*this->X, *this->A,this->interfaceOrder-1,
						      U, this->V0,
						      this->Weights,this->VWeights,
						      NULL, *this->fluidSelector.fluidId, *this->fluidSelector.fluidIdn,
                                                      this->limitHigherOrderExtrapolation);

    this->varFcn->primitiveToConservative(this->V0,U,this->fluidSelector.fluidId);
  }

  this->checkSolution(U);

 if (this->modifiedGhidaglia) {
    *this->bcData->getBoundaryStateHH() = embeddedU.hh();
    this->timeState->updateHH(embeddedU.hh());
  }

  this->domain->setExactBoundaryValues(U, *this->X, this->ioData, 
				       this->currentTime + this->currentTimeStep,
				       this->spaceOp->getVarFcn());

  this->updateBoundaryExternalState();
  return its;
  
}

template<int dim,int dimLS>
void ImplicitLevelSetTsDesc<dim,dimLS>::doErrorEstimation(DistSVec<double,dim> &U) 
{
  DistSVec<double,dim> *F = ns->GetResidual();

  if (this->lsMethod == 0) {
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,U);
    this->multiPhaseSpaceOp->computeResidual(*this->X, *this->A, this->timeState->getUn(), this->PhiV, 
					     this->fluidSelector, *F, this->riemann, 
                                             this->timeState,0);
  } else {

    this->multiPhaseSpaceOp->computeResidual(*this->X, *this->A, this->timeState->getUn(), this->Phi, 
					     this->fluidSelector, *F, this->riemann, 
                                             this->timeState,0);
  }

  this->timeState->calculateErrorEstiNorm(U, *F); 
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// External routines to solve Euler equations implicitly (called by NewtonSolver)
//------------------------------------------------------------------------------

// this function evaluates (Aw),t + F(w,x,v)
template<int dim, int dimLS>
void ImplicitLevelSetTsDesc<dim,dimLS>::computeFunction(int it, DistSVec<double,dim> &Q,
                                                  DistSVec<double,dim> &F)
{
  // phi is obtained once and for all for this iteration
  // no need to recompute it before computation of jacobian.

  DistSVec<double,dimLS>* locphi = &this->Phi;
  if (this->lsMethod == 0)
    locphi = &this->PhiV;
  //else if (it > 1 &&  timeType == ExplicitData::RUNGE_KUTTA_2)
  //  locphi = this->Phi0;

  this->domain->setExactBoundaryValues(Q, *this->X, this->ioData, 
				       this->currentTime + this->currentTimeStep,
				       this->spaceOp->getVarFcn());
  
  if (this->lsMethod == 0) {
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,Q);
    this->multiPhaseSpaceOp->computeResidual(*this->X, *this->A, Q, this->PhiV, 
					     this->fluidSelector, F, this->riemann,
                                             this->timeState, 1);
  } else {

    this->multiPhaseSpaceOp->computeResidual(*this->X, *this->A, Q, this->Phi, 
					     this->fluidSelector, F, this->riemann, 
                                             this->timeState,1);
  }
  
  this->timeState->add_dAW_dt(it, *this->geoState, *this->A, Q, F);
  this->multiPhaseSpaceOp->applyBCsToResidual(Q, F);
  this->domain->setExactBoundaryResidual(F, *this->X, this->ioData, 
					 this->currentTime + this->currentTimeStep,
					 this->spaceOp->getVarFcn());

 
  if (this->modifiedGhidaglia) {

    *hhResidual = 0.0;
    this->domain->
      computeHHBoundaryTermResidual(*this->bcData,Q,*hhResidual, this->varFcn);
       
    this->timeState->add_dAW_dt_HH(-1, *this->geoState, *this->A,*this->bcData->getBoundaryStateHH()
    			     , *hhResidual);
  }

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void ImplicitLevelSetTsDesc<dim,dimLS>::recomputeFunction(DistSVec<double,dim> &Q,
                                            DistSVec<double,dim> &rhs)
{
  this->multiPhaseSpaceOp->recomputeRHS(*this->X, Q, *this->fluidSelector.fluidId, rhs);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
int ImplicitLevelSetTsDesc<dim,dimLS>::checkFailSafe(DistSVec<double,dim>& U)
{
  //fprintf(stdout, "CheckFailSafe for ImplicitLevelSetTsDesc to be rewritten\n");

  if (!this->failSafeNewton) return 0;

  if (!this->tag)
    this->tag = new DistSVec<bool,2>(this->getVecInfo());

  this->domain->checkFailSafe(this->varFcn, U, *this->tag, this->fluidSelector.fluidId);
  this->multiPhaseSpaceOp->fix(*this->tag);

  return this->failSafeNewton;

}

//------------------------------------------------------------------------------
template<int dim, int dimLS> 
void ImplicitLevelSetTsDesc<dim,dimLS>::resetFixesTag()
{

  this->multiPhaseSpaceOp->resetTag();

}

//------------------------------------------------------------------------------
template<int dim, int dimLS>
void ImplicitLevelSetTsDesc<dim,dimLS>::computeJacobian(int it, DistSVec<double,dim> &Q,
							DistSVec<double,dim> &F)
{

  MatVecProdH1<dim,double,dim> *mvph1 = dynamic_cast<MatVecProdH1<dim,double,dim> *>(mvp);
  if (this->modifiedGhidaglia)
    mvp->evaluateHH(*this->hhResidual, *this->bcData->getBoundaryStateHH());
  
  if (this->lsMethod == 0) {
    mvp->evaluate(it, *this->X, *this->A, Q, this->PhiV, F);
  } else
    mvp->evaluate(it, *this->X, *this->A, Q, this->Phi, F);
  
  if (mvph1) 
    this->domain->setExactBoundaryJacobian(Q, *this->X, this->ioData, 
					   this->currentTime + this->currentTimeStep,
					   this->spaceOp->getVarFcn(), *mvph1);
}
//------------------------------------------------------------------------------
template<int dim, int dimLS>
void ImplicitLevelSetTsDesc<dim,dimLS>::setOperators(DistSVec<double,dim> &Q)
{
  
  DistMat<PrecScalar,dim> *_pc = dynamic_cast<DistMat<PrecScalar,dim> *>(pc);
  
  if (_pc) {
    
    MatVecProdFDMultiPhase<dim,dimLS> *mvpfd = dynamic_cast<MatVecProdFDMultiPhase<dim,dimLS> *>(mvp);
    MatVecProdH1MultiPhase<dim,dimLS> *mvph1 = dynamic_cast<MatVecProdH1MultiPhase<dim,dimLS> *>(mvp);
    
    if (mvpfd) {

      this->multiPhaseSpaceOp->computeJacobian(*this->X, *this->A, Q, *_pc, this->fluidSelector, this->riemann,this->timeState);
      this->timeState->addToJacobian(*this->A, *_pc, Q);
      this->multiPhaseSpaceOp->applyBCsToJacobian(Q, *_pc);
    }
    else if (mvph1) {
      JacobiPrec<PrecScalar,dim> *jac = dynamic_cast<JacobiPrec<PrecScalar,dim> *>(pc);
      IluPrec<PrecScalar,dim> *ilu = dynamic_cast<IluPrec<PrecScalar,dim> *>(pc);
      
      if (jac)
	jac->getData(*mvph1);
      else if (ilu)
	ilu->getData(*mvph1);
    }
    
  }
  
  double t0 = this->timer->getTime();
  
  pc->setup();
  
  double t = this->timer->addPrecSetupTime(t0);
  
  this->com->printf(6, "Fluid preconditioner computation: %f s\n", t);
  
}

//------------------------------------------------------------------------------
template<int dim, int dimLS>
int ImplicitLevelSetTsDesc<dim,dimLS>::solveLinearSystem(int it, DistSVec<double,dim> &b,
						  DistSVec<double,dim> &dQ)
{
  
  double t0 = this->timer->getTime();

  dQ = 0.0;

  this->embeddeddQ = 0.0;
  this->embeddedB.ghost() = 0.0;
  this->embeddedB.real() = b;
   
  if (this->modifiedGhidaglia)
    this->embeddedB.hh() = -1.0*(*this->hhResidual); 

  ksp->setup(it, this->maxItsNewton, this->embeddedB);

  int lits = ksp->solve(this->embeddedB, this->embeddeddQ);

  if(lits==kspLS->maxits) this->errorHandler->localErrors[ErrorHandler::SATURATED_LS] += 1;
 
  dQ = this->embeddeddQ.real();
//  this->embeddedU.ghost() += this->embeddeddQ.ghost();
  if (this->modifiedGhidaglia) {
    this->embeddedU.hh() += this->embeddeddQ.hh();

    *this->bcData->getBoundaryStateHH() = this->embeddedU.hh();
  }
 
  //mvpLS->apply(dQ, fnew);
 
  this->timer->addLSKspTime(t0);

  return lits;

}

//------------------------------------------------------------------------------
// External routines to solve LevelSet equation implicitly (called by NewtonSolver)
//------------------------------------------------------------------------------
// this function evaluates (Aw),t + F(w,x,v)
template<int dim, int dimLS>
void ImplicitLevelSetTsDesc<dim,dimLS>::computeFunctionLS(int it,
                                                    DistSVec<double,dim> &U,
                                                    DistSVec<double,dimLS> &PhiF)
{
  this->multiPhaseSpaceOp->computeResidualLS(*this->X, *this->A, this->Phi, *this->fluidSelector.fluidId, U, PhiF,0,0,this->lsMethod);

  this->timeState->add_dAW_dtLS(it, *this->geoState, *this->A, this->Phi,
			        this->LS->Phin, this->LS->Phinm1, 
      				this->LS->Phinm2, PhiF,this->requireSpecialBDF);

}

//------------------------------------------------------------------------------
template<int dim, int dimLS>
void ImplicitLevelSetTsDesc<dim,dimLS>::computeJacobianLS(int it,
                                                    DistSVec<double,dim> &U,
                                                    DistSVec<double,dimLS> &PhiF)
{
  mvpLS->evaluate(it, *this->X, *this->A, this->Phi,
		  U,this->V0, PhiF, *this->fluidSelector.fluidId,this->requireSpecialBDF, 0,this->lsMethod);
}

//------------------------------------------------------------------------------
template<int dim, int dimLS>
void ImplicitLevelSetTsDesc<dim,dimLS>::setOperatorsLS(DistSVec<double,dimLS> &Q)
{
  
  DistMat<PrecScalar,dimLS> *_pc = dynamic_cast<DistMat<PrecScalar,dimLS> *>(pcLS);
  
  if (_pc) {
    
      JacobiPrec<PrecScalar,dimLS> *jac = dynamic_cast<JacobiPrec<PrecScalar,dimLS> *>(pcLS);
      IluPrec<PrecScalar,dimLS> *ilu = dynamic_cast<IluPrec<PrecScalar,dimLS> *>(pcLS);
      
      if (jac)
	jac->getData(*mvpLS);
      else if (ilu)
	ilu->getData(*mvpLS);
    
  }
  
  double t0 = this->timer->getTime();
  
  pcLS->setup();
  
  double t = this->timer->addLSPrecSetupTime(t0);
  
  this->com->printf(6, "Fluid preconditioner computation: %f s\n", t);
  
}
//------------------------------------------------------------------------------

template<int dim, int dimLS>
int ImplicitLevelSetTsDesc<dim,dimLS>::solveLinearSystemLS(int it, DistSVec<double,dimLS> &b,
                                                  DistSVec<double,dimLS> &dQ)
{
  double t0 = this->timer->getTime();

  dQ = 0.0;

  kspLS->setup(it, this->maxItsNewton, b);

  int lits = kspLS->solve(b, dQ);

  //mvpLS->apply(dQ, fnew);
  
  //  this->domain->getCommunicator()->fprintf(stdout,"%e\n",dQ.norm());
  
  this->timer->addLSKspTime(t0);
  
  return lits;
}

//------------------------------------------------------------------------------

