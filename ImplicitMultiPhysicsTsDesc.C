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

template<int dim,int dimLS>
ImplicitMultiPhysicsTsDesc<dim,dimLS>::
ImplicitMultiPhysicsTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  MultiPhysicsTsDesc<dim,dimLS>(ioData, geoSource, dom) , embeddedU(dom->getNodeDistInfo()), embeddedB(dom->getNodeDistInfo()),embeddeddQ(dom->getNodeDistInfo())
{
  tag = 0;
  ImplicitData &implicitData = ioData.ts.implicit;
  
  // NewtonSolver
  ns = new NewtonSolver<ImplicitMultiPhysicsTsDesc<dim,dimLS> >(this);
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
                                      this->multiPhaseSpaceOp,this->riemann,&(this->fluidSelector),
                                      this->domain,ioData);
  else if (implicitData.mvp == ImplicitData::H1)
    mvp = new MatVecProdH1MultiPhase<dim,dimLS>(this->timeState, this->multiPhaseSpaceOp, this->riemann,&(this->fluidSelector),this->domain);
  else{
    this->com->fprintf(stdout, "*** Error: MatVecProdH2 is not available\n");
    exit(1);
  }
  
  pc = createPreconditioner<dim>(implicitData.newton.ksp.ns.pc, this->domain);
  
  ksp = createKrylovSolver(this->getVecInfo(), implicitData.newton.ksp.ns, mvp, pc, this->com);

  // MatVecProd, Prec and Krylov solver for LevelSet equation
  mvpLS  = new MatVecProdLS<dim,dimLS>(this->timeState, this->geoState, this->multiPhaseSpaceOp, this->domain, this->LS);

  pcLS = createPreconditioner<dimLS>(implicitData.newton.ksp.lsi.pc, this->domain);
  kspLS = createKrylovSolverLS(this->getVecInfo(), implicitData.newton.ksp.lsi, mvpLS, pcLS, this->com);
 
 
  //initialize mmh (EmbeddedMeshMotionHandler).
  if(this->dynNodalTransfer) 
    {
      /*
      MeshMotionHandler *_mmh = 0;
      _mmh = new EmbeddedMeshMotionHandler(ioData, dom, this->dynNodalTransfer, this->distLSS);
      this->mmh = _mmh;
      */
      this->mmh = new EmbeddedMeshMotionHandler(ioData, dom, this->dynNodalTransfer, this->distLSS);
    } 
  else
    { 
      this->mmh = 0;
    }

  typename MatVecProdMultiPhase<dim,dimLS>::_fsi fsi = {

    this->distLSS,
    this->fluidSelector.fluidId,
    this->riemann,
    this->linRecAtInterface,
    this->viscSecOrder,
    &this->Wtemp,
    this->riemannNormal,
    this->ghostPoints,
  };

  mvp->AttachStructure(fsi);
  
  this->existsWstarnm1 = false;

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

template<int dim,int dimLS>
ImplicitMultiPhysicsTsDesc<dim,dimLS>::~ImplicitMultiPhysicsTsDesc()
{
  if (tag)   delete tag;
  if (mvp)   delete mvp;
  if (pc)    delete pc;
  if (ksp)   delete ksp;
  if (ns)    delete ns;

}

//------------------------------------------------------------------------------
//  Internal routines to setup the class (called in constructor)
//------------------------------------------------------------------------------

template<int dim,int dimLS>
template <int neq>
KspPrec<neq> *ImplicitMultiPhysicsTsDesc<dim,dimLS>::createPreconditioner(PcData &pcdata, Domain *dom)
{
  
  KspPrec<neq> *_pc = 0;
  
  if (pcdata.type == PcData::IDENTITY)
    _pc = new IdentityPrec<neq>();
  else if (pcdata.type == PcData::JACOBI)
    _pc = new JacobiPrec<double,neq>(DiagMat<double,neq>::DENSE, dom);
  else if (pcdata.type == PcData::AS ||
	   pcdata.type == PcData::RAS ||
	   pcdata.type == PcData::ASH ||
	   pcdata.type == PcData::AAS)
    _pc = new IluPrec<double,neq>(pcdata, dom);
  
  return _pc;
  
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
template<int neq, class MatVecProdOp>
KspSolver<DistEmbeddedVec<double,neq>, MatVecProdOp, KspPrec<neq>, Communicator> *
ImplicitMultiPhysicsTsDesc<dim,dimLS>::createKrylovSolver(
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
ImplicitMultiPhysicsTsDesc<dim,dimLS>::createKrylovSolverLS(							    const DistInfo &info, KspData &kspdata,
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
template<int dim,int dimLS>
void ImplicitMultiPhysicsTsDesc<dim,dimLS>::commonPart(DistSVec<double,dim> &U)
{
  // Adam 04/06/10: Took everything in common in solveNLAllFE and solveNLAllRK2 and put it here. Added Ghost-Points treatment for viscous flows.

  if(this->mmh && !this->inSubCycling) {
    //get structure timestep dts
    this->dts = this->mmh->update(0, 0, 0, this->bcData->getVelocityVector(), *this->Xs);

    //recompute intersections
    double tw = this->timer->getTime();
    if(this->withCracking && this->withMixedLS) // no need for the intersector to determine fluidId.
      this->distLSS->recompute(this->dtf, this->dtfLeft, this->dts, false, TsDesc<dim>::failSafeFlag); 
    else
      this->distLSS->recompute(this->dtf, this->dtfLeft, this->dts, true, TsDesc<dim>::failSafeFlag); 

    this->timer->addIntersectionTime(tw);
    this->com->barrier();
    this->timer->removeIntersAndPhaseChange(tw);
   
    if (this->lsMethod == 0)
      this->LS->conservativeToPrimitive(this->Phi, this->PhiV, U);
    else
      this->PhiV = this->Phi;

    //this->multiPhaseSpaceOp->extrapolatePhiV(this->distLSS, this->PhiV);
    if(this->withCracking && this->withMixedLS) {
      //this->multiPhaseSpaceOp->extrapolatePhiV2(this->distLSS, this->PhiV);
      //this->fluidSelector.updateFluidIdFS2(this->distLSS, this->PhiV);
      DistSVec<bool,4> poll(this->domain->getNodeDistInfo());
      this->domain->updateFluidIdFS2Prep(*(this->distLSS), this->PhiV, *(this->fluidSelector.fluidId), poll);
      this->fluidSelector.updateFluidIdFS2(this->distLSS, this->PhiV, poll);
    } else {
      this->fluidSelector.updateFluidIdFS(this->distLSS, this->PhiV);
    }
    this->PhiV = 0.0; //PhiV is no longer a distance function now. Only its sign (+/-)
                       //  is meaningful. We destroy it so people wouldn't use it
                       //  by mistake later on.

    //update phase-change
    tw = this->timer->getTime();
    this->multiPhaseSpaceOp->updateSweptNodes(*this->X, this->phaseChangeChoice, U, this->Vtemp, *this->Weights, *this->VWeights,
                                              this->LS->Phin, this->PhiWeights, *this->Wstarij, *this->Wstarji,
                                              this->distLSS, this->vfar, (this->withCracking && this->withMixedLS), this->fluidSelector.fluidIdn, this->fluidSelector.fluidId);

    this->timeState->getUn() = U;

    // BDF update (Unm1)
    if (this->timeState->useNm1() && this->timeState->existsNm1()) {
      tw = this->timer->getTime();
      DistSVec<double,dim>& Unm1 = this->timeState->getUnm1();
  
      if (!this->existsWstarnm1) {
        fprintf(stderr,"*** Error: I ignored this case!\n");
        exit(1);  
      }

      this->multiPhaseSpaceOp->updateSweptNodes(*this->X, this->phaseChangeChoice, Unm1, this->Vtemp, *this->Weights, *this->VWeights,
                                                this->LS->Phinm1, this->PhiWeights, *this->Wstarij_nm1, *this->Wstarji_nm1,
                                                this->distLSS, this->vfar, (this->withCracking && this->withMixedLS), this->fluidSelector.fluidIdn, this->fluidSelector.fluidId);
      this->timer->addEmbedPhaseChangeTime(tw);
      this->timer->removeIntersAndPhaseChange(tw);
    }

    if (this->timeState->useNm1()) {
      *this->Wstarij_nm1 = *this->Wstarij;
      *this->Wstarji_nm1 = *this->Wstarji;
      this->existsWstarnm1 = true;
    }
 
  }
  else if(this->mmh && this->inSubCycling) {
     // PJSA: reset after ErrorHandler::REDO_TIMESTEP is raised
     U = this->timeState->getUn();
     *this->fluidSelector.fluidId = *this->fluidSelector.fluidIdn;
  }

  *this->fluidSelector.fluidIdn = *this->fluidSelector.fluidId;

  // Ghost-Points Population
/*  if(this->eqsType == TsDesc<dim>::NAVIER_STOKES)
    {
      this->ghostPoints->deletePointers();
      this->mulitPhaseSpaceOp->populateGhostPoints(this->ghostPoints,*this->X,U,this->varFcn,this->distLSS,this->viscSecOrder,this->nodeTag);
    }
*/
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// External routine to solve problem (called by TsSolver)
// It calls for the NewtonSolver ns, which in turn will
// call routines below from this same file or from LevelSetTsDesc
//------------------------------------------------------------------------------
template<int dim,int dimLS>
int ImplicitMultiPhysicsTsDesc<dim,dimLS>::solveNonLinearSystem(DistSVec<double,dim> &U, int)
{ 
  double t0 = this->timer->getTime();
  DistVec<int> fluidId_copy(*this->fluidSelector.fluidId);
  commonPart(U);

  if (this->modifiedGhidaglia)
    embeddedU.hh() = *this->bcData->getBoundaryStateHH();

 
  if (this->timeState->useNm1() && this->timeState->existsNm1() &&
      this->phaseChangeType == 1) {
    DistSVec<double,dim>& Unm1 = this->timeState->getUnm1();
  
    this->varFcn->conservativeToPrimitive(U,this->V0,this->fluidSelector.fluidId);
    this->multiPhaseSpaceOp->extrapolatePhaseChange(*this->X, *this->A, this->multiFluidInterfaceOrder-1,
 					            Unm1, this->V0,
						    *this->Weights,*this->VWeights,
						    NULL, *this->fluidSelector.fluidId, *this->fluidSelector.fluidIdnm1,
						    false);
    this->varFcn->primitiveToConservative(this->V0,Unm1,this->fluidSelector.fluidId);
  }
    

 
  int its = this->ns->solve(U);
 
  this->errorHandler->reduceError();
  this->data->resolveErrors();

  this->timer->addFluidSolutionTime(t0);

  if(this->errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP]) {
    *this->fluidSelector.fluidId = fluidId_copy;
    return its;
  }
 
  this->varFcn->conservativeToPrimitive(U,this->V0,this->fluidSelector.fluidId);
  this->riemann->storePreviousPrimitive(this->V0, *this->fluidSelector.fluidId, *this->X);
    
  double t1 = this->timer->getTime();
  int itsLS = this->ns->solveLS(this->Phi, U);
  this->riemann->storeOldV(U);
  if(this->withCracking && this->withMixedLS)
    this->riemann->avoidNewPhaseCreation(this->Phi, this->LS->Phin);
  else
    this->riemann->avoidNewPhaseCreation(this->Phi, this->LS->Phin,this->distLSS);
  
  
//  this->fluidSelector.getFluidId(this->Phi,&(this->distLSS->getStatus()));
  DistVec<int> fluidId0(*this->fluidSelector.fluidId);
  if(this->withCracking && this->withMixedLS) {
    this->fluidSelector.updateFluidIdFF2(this->distLSS, this->Phi);
  } else {
    this->fluidSelector.updateFluidIdFF(this->distLSS, this->Phi);
  }

  this->timer->addLevelSetSolutionTime(t1);

  if (this->phaseChangeType == 0)
    this->riemann->updatePhaseChange(this->V0, *this->fluidSelector.fluidId, fluidId0);
  else
    this->multiPhaseSpaceOp->extrapolatePhaseChange(*this->X,*this->A, this->multiFluidInterfaceOrder-1,
						    U, this->V0,
						    *this->Weights,*this->VWeights,
						    NULL, *this->fluidSelector.fluidId,fluidId0,
						    this->limitHigherOrderExtrapolation);

  this->varFcn->primitiveToConservative(this->V0,U,this->fluidSelector.fluidId);

  this->checkSolution(U);

  this->updateBoundaryExternalState();
  
  return its;
}
//------------------------------------------------------------------------------
// External routines to solve Euler equations implicitly (called by NewtonSolver)
//------------------------------------------------------------------------------

// this function evaluates (Aw),t + F(w,x,v)
template<int dim,int dimLS>
void ImplicitMultiPhysicsTsDesc<dim,dimLS>::computeFunction(int it, DistSVec<double,dim> &Q,
                                                  DistSVec<double,dim> &F)
{
  // phi is obtained once and for all for this iteration
  // no need to recompute it before computation of jacobian.
  if (this->lsMethod == 0) {  
    this->LS->conservativeToPrimitive(this->Phi,this->PhiV,Q);
    this->multiPhaseSpaceOp->computeResidual(*this->X, *this->A, Q, *this->Wstarij, *this->Wstarji, this->distLSS,
                                   this->linRecAtInterface, this->viscSecOrder, this->riemann,  
                                   this->riemannNormal, this->PhiV, this->fluidSelector,F, 1, this->ghostPoints);
  } else {

    this->multiPhaseSpaceOp->computeResidual(*this->X, *this->A, Q, *this->Wstarij, *this->Wstarji, this->distLSS,
                                   this->linRecAtInterface, this->viscSecOrder, this->riemann,  
                                   this->riemannNormal, this->Phi, this->fluidSelector,F, 1, this->ghostPoints);
  }
  this->timeState->add_dAW_dt(it, *this->geoState, *this->A, Q, F);
  this->multiPhaseSpaceOp->applyBCsToResidual(Q, F);

  if (this->modifiedGhidaglia) {

    *hhResidual = 0.0;
    this->domain->
      computeHHBoundaryTermResidual(*this->bcData,Q,*hhResidual, this->varFcn);
       
    this->timeState->add_dAW_dt_HH(-1, *this->geoState, *this->A,*this->bcData->getBoundaryStateHH()
    			     , *hhResidual);
  }

}

//------------------------------------------------------------------------------

template<int dim,int dimLS>
void ImplicitMultiPhysicsTsDesc<dim,dimLS>::recomputeFunction(DistSVec<double,dim> &Q,
                                            DistSVec<double,dim> &rhs)
{
  this->multiPhaseSpaceOp->recomputeRHS(*this->X, Q, rhs);
}

//------------------------------------------------------------------------------

template<int dim,int dimLS>
int ImplicitMultiPhysicsTsDesc<dim,dimLS>::checkFailSafe(DistSVec<double,dim>& U)
{
//  this->com->fprintf(stdout, "WARNING: At the moment CheckFailSafe is not supported by the embedded framework with an implicit time-integrator!\n");

  if (!this->failSafeNewton) return 0;

  if (!this->tag)
    this->tag = new DistSVec<bool,2>(this->getVecInfo());

  this->domain->checkFailSafe(this->varFcn, U, *this->tag, this->fluidSelector.fluidId);
  this->multiPhaseSpaceOp->fix(*this->tag);

  return 1;

}

//------------------------------------------------------------------------------
template<int dim,int dimLS> 
void ImplicitMultiPhysicsTsDesc<dim,dimLS>::resetFixesTag()
{

  this->multiPhaseSpaceOp->resetTag();

}

//------------------------------------------------------------------------------
template<int dim,int dimLS>
void ImplicitMultiPhysicsTsDesc<dim,dimLS>::computeJacobian(int it, DistSVec<double,dim> &Q,
							DistSVec<double,dim> &F)
{

  if (this->modifiedGhidaglia)
    mvp->evaluateHH(*this->hhResidual, *this->bcData->getBoundaryStateHH());
  
  if (this->lsMethod == 0)
    mvp->evaluate(it,*(this->X) ,*(this->A), Q,this->PhiV, F);
  else
    mvp->evaluate(it,*(this->X) ,*(this->A), Q,this->Phi, F);

}
//------------------------------------------------------------------------------
template<int dim,int dimLS>
void ImplicitMultiPhysicsTsDesc<dim,dimLS>::setOperators(DistSVec<double,dim> &Q)
{
  
  DistMat<double,dim> *_pc = dynamic_cast<DistMat<double,dim> *>(pc);
  
  if (_pc) {
    
    MatVecProdFDMultiPhase<dim,dimLS> *mvpfd = dynamic_cast<MatVecProdFDMultiPhase<dim,dimLS> *>(mvp);
    MatVecProdH1MultiPhase<dim,dimLS> *mvph1 = dynamic_cast<MatVecProdH1MultiPhase<dim,dimLS> *>(mvp);
    
    if (mvpfd) {

      this->multiPhaseSpaceOp->computeJacobian(this->riemann, *this->X, Q,*this->A,this->distLSS,
                                   this->riemannNormal, (this->fluidSelector),*_pc,this->timeState);
      this->timeState->addToJacobian(*this->A, *_pc, Q);
      this->multiPhaseSpaceOp->applyBCsToJacobian(Q, *_pc);
    }
    else if (mvph1) {
      JacobiPrec<double,dim> *jac = dynamic_cast<JacobiPrec<double,dim> *>(pc);
      IluPrec<double,dim> *ilu = dynamic_cast<IluPrec<double,dim> *>(pc);
      
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
template<int dim,int dimLS>
int ImplicitMultiPhysicsTsDesc<dim,dimLS>::solveLinearSystem(int it, DistSVec<double,dim> &b,
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
  
  if(lits==ksp->maxits) this->errorHandler->localErrors[ErrorHandler::SATURATED_LS] += 1;
  
  dQ = this->embeddeddQ.real();
//  this->embeddedU.ghost() += this->embeddeddQ.ghost();
  if (this->modifiedGhidaglia) {
    this->embeddedU.hh() += this->embeddeddQ.hh();

    *this->bcData->getBoundaryStateHH() = this->embeddedU.hh();
  }
 
  this->timer->addKspTime(t0);
  
  return lits;
  
}

//------------------------------------------------------------------------------
// External routines to solve LevelSet equation implicitly (called by NewtonSolver)
//------------------------------------------------------------------------------
// this function evaluates (Aw),t + F(w,x,v)
template<int dim, int dimLS>
void ImplicitMultiPhysicsTsDesc<dim,dimLS>::computeFunctionLS(int it,
                                                    DistSVec<double,dim> &U,
                                                    DistSVec<double,dimLS> &PhiF)
{
  this->multiPhaseSpaceOp->computeResidualLS(*this->X, *this->A, this->Phi, *this->fluidSelector.fluidId, U, PhiF,this->distLSS,this->linRecAtInterface,this->lsMethod);

  this->timeState->add_dAW_dtLS(it, *this->geoState, *this->A, this->Phi,
                                this->LS->Phin, this->LS->Phinm1,
                                this->LS->Phinm2, PhiF,this->requireSpecialBDF);

}

//------------------------------------------------------------------------------
template<int dim, int dimLS>
void ImplicitMultiPhysicsTsDesc<dim,dimLS>::computeJacobianLS(int it,
                                                    DistSVec<double,dim> &U,
                                                    DistSVec<double,dimLS> &PhiF)
{
  mvpLS->evaluate(it, *this->X, *this->A, this->Phi,
                  U,this->V0, PhiF, *this->fluidSelector.fluidId,this->requireSpecialBDF,this->distLSS,this->lsMethod);
}

//------------------------------------------------------------------------------
template<int dim, int dimLS>
void ImplicitMultiPhysicsTsDesc<dim,dimLS>::setOperatorsLS(DistSVec<double,dimLS> &Q)
{

  DistMat<double,dimLS> *_pc = dynamic_cast<DistMat<double,dimLS> *>(pcLS);

  if (_pc) {

      JacobiPrec<double,dimLS> *jac = dynamic_cast<JacobiPrec<double,dimLS> *>(pcLS);
      IluPrec<double,dimLS> *ilu = dynamic_cast<IluPrec<double,dimLS> *>(pcLS);

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
int ImplicitMultiPhysicsTsDesc<dim,dimLS>::solveLinearSystemLS(int it, DistSVec<double,dimLS> &b,
                                                  DistSVec<double,dimLS> &dQ)
{

  double t0 = this->timer->getTime();

  dQ = 0.0;
  
  kspLS->setup(it, this->maxItsNewton, b);
  
  int lits = kspLS->solve(b, dQ);

  if(lits==kspLS->maxits) this->errorHandler->localErrors[ErrorHandler::SATURATED_LS] += 1;

  //mvpLS->apply(dQ, fnew);

  this->timer->addLSKspTime(t0);

  return lits;

}
  
//------------------------------------------------------------------------------

