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

template<int dim>
ImplicitEmbeddedTsDesc<dim>::
ImplicitEmbeddedTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  EmbeddedTsDesc<dim>(ioData, geoSource, dom), embeddedU(dom->getNodeDistInfo()), embeddedB(dom->getNodeDistInfo()),embeddeddQ(dom->getNodeDistInfo())
{
  tag = 0;
  ImplicitData &implicitData = ioData.ts.implicit;
  
  // NewtonSolver
  ns = new NewtonSolver<ImplicitEmbeddedTsDesc<dim> >(this);
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

  //initialize emmh (EmbeddedMeshMotionHandler).
  if(this->dynNodalTransfer) 
    {
      /*
      MeshMotionHandler *_mmh = 0;
      _mmh = new EmbeddedMeshMotionHandler(ioData, dom, this->dynNodalTransfer, this->distLSS);
      this->emmh = _mmh;
      */
      this->emmh = new EmbeddedMeshMotionHandler(ioData, dom, this->dynNodalTransfer, this->distLSS);
    } 
  else
    { 
      this->emmh = 0;
    }

  this->existsWstarnm1 = false;

  if (ioData.problem.framework==ProblemData::EMBEDDEDALE && this->emmh) 
    this->mmh = this->createEmbeddedALEMeshMotionHandler(ioData, geoSource, this->distLSS);
  else
    this->mmh = 0;

  if (this->modifiedGhidaglia) {
    embeddedU.addHHBoundaryTerm(dom->getFaceDistInfo());
    embeddeddQ.addHHBoundaryTerm(dom->getFaceDistInfo());
    embeddedB.addHHBoundaryTerm(dom->getFaceDistInfo());
    
    hhResidual = new DistVec<double>(dom->getFaceDistInfo());
  }

    
}

//------------------------------------------------------------------------------

template<int dim>
ImplicitEmbeddedTsDesc<dim>::~ImplicitEmbeddedTsDesc()
{
  if (tag)   delete tag;
  if (ns)    delete ns;

}

//------------------------------------------------------------------------------
//  Internal routines to setup the class (called in constructor)
//------------------------------------------------------------------------------

template<int dim>
template <int neq>
KspPrec<neq> *ImplicitEmbeddedTsDesc<dim>::createPreconditioner(PcData &pcdata, Domain *dom)
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

template<int dim>
template<int neq, class MatVecProdOp>
KspSolver<DistEmbeddedVec<double,neq>, MatVecProdOp, KspPrec<neq>, Communicator> *
ImplicitEmbeddedTsDesc<dim>::createKrylovSolver(
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

template<int dim>
int ImplicitEmbeddedTsDesc<dim>::commonPart(DistSVec<double,dim> &U)
{
  // Adam 04/06/10: Took everything in common in solveNLAllFE and solveNLAllRK2 and put it here. Added Ghost-Points treatment for viscous flows.

  //if(this->emmh && !this->inSubCycling) {  //subcycling is allowed from now on
	if(this->emmh) 
	{
    int failSafe=0;

		if(TsDesc<dim>::failSafeFlag == false)
		{
      *this->WstarijCopy = *this->Wstarij;
      *this->WstarjiCopy = *this->Wstarji;

			if(this->timeState->useNm1()) 
			{
        *this->Wstarij_nm1Copy = *this->Wstarij_nm1;
        *this->Wstarji_nm1Copy = *this->Wstarji_nm1;
      }
      *this->nodeTagCopy = this->distLSS->getStatus();
      *EmbeddedTsDesc<dim>::UCopy = U;
    }
		else
		{
      *this->Wstarij = *this->WstarijCopy;
      *this->Wstarji = *this->WstarjiCopy;

			if(this->timeState->useNm1()) 
			{
        *this->Wstarij_nm1 = *this->Wstarij_nm1Copy;
        *this->Wstarji_nm1 = *this->Wstarji_nm1Copy;
      }
      this->distLSS->setStatus(*this->nodeTagCopy);
      this->nodeTag = *this->nodeTagCopy;
    }

    //get structure timestep dts
    this->dts = this->emmh->update(0, 0, 0, this->bcData->getVelocityVector(), *this->Xs);


    //recompute intersections
    double tw = this->timer->getTime();

		failSafe = this->distLSS->recompute(this->dtf, this->dtfLeft, 
														this->dts, true, TsDesc<dim>::failSafeFlag);
		
    this->com->globalMin(1, &failSafe);
    if(failSafe<0) //in case of intersection failure -1 is returned by recompute 
       return failSafe;

    this->timer->addIntersectionTime(tw);
    this->com->barrier();
    this->timer->removeIntersAndPhaseChange(tw);

    //update nodeTags (only for numFluid>1)
 //   if(this->numFluid>1) {
      this->nodeTag0 = this->nodeTag;
      this->nodeTag = this->distLSS->getStatus();
 //   }

    //store previous states for phase-change update
    tw = this->timer->getTime();

		this->spaceOp->updateSweptNodes(*this->X, *this->A, 
												  this->phaseChangeChoice, this->phaseChangeAlg, 
												  U, this->Vtemp,  *this->Weights, *this->VWeights, 
												  *this->Wstarij, *this->Wstarji, this->distLSS, 
												  (double*)this->vfar, 
												  this->ioData.embed.interfaceLimiter == EmbeddedFramework::LIMITERALEX1, 
												  &this->nodeTag);

    this->timer->addEmbedPhaseChangeTime(tw);
    this->timer->removeIntersAndPhaseChange(tw);

    //this->timeState->update(U); 
    this->timeState->getUn() = U;

    // BDF update (Unm1)
		if (this->timeState->useNm1() && this->timeState->existsNm1()) 
		{
      tw = this->timer->getTime();
      DistSVec<double,dim>& Unm1 = this->timeState->getUnm1();
  
			if (!this->existsWstarnm1) 
			{        
				this->spaceOp->computeResidual(*this->X, *this->A, Unm1, 
														 *this->Wstarij, *this->Wstarji, *this->Wextij, this->distLSS,
														 this->linRecAtInterface, this->viscSecOrder, 
														 this->nodeTag, this->Vtemp, this->riemann, 
                                 this->riemannNormal, 1, this->ghostPoints);
      }

      tw = this->timer->getTime();

			this->spaceOp->updateSweptNodes(*this->X,*this->A, 
													  this->phaseChangeChoice, this->phaseChangeAlg, 
													  Unm1, this->Vtemp, *this->Weights, *this->VWeights, 
													  *this->Wstarij_nm1, *this->Wstarji_nm1,
				      this->distLSS, (double*)this->vfar,
													  this->ioData.embed.interfaceLimiter == EmbeddedFramework::LIMITERALEX1,
													  &this->nodeTag);

      this->timer->addEmbedPhaseChangeTime(tw);
      this->timer->removeIntersAndPhaseChange(tw);

      this->timer->addEmbedPhaseChangeTime(tw);
      this->timer->removeIntersAndPhaseChange(tw);
    }

		if (this->timeState->useNm1()) 
		{
      *this->Wstarij_nm1 = *this->Wstarij;
      *this->Wstarji_nm1 = *this->Wstarji;
      this->existsWstarnm1 = true;
    }
  
  }

  if (this->modifiedGhidaglia)
    embeddedU.hh() = *this->bcData->getBoundaryStateHH();

  return 0;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// External routine to solve problem (called by TsSolver)
// It calls for the NewtonSolver ns, which in turn will
// call routines below from this same file or from LevelSetTsDesc
//------------------------------------------------------------------------------
template<int dim>
int ImplicitEmbeddedTsDesc<dim>::solveNonLinearSystem(DistSVec<double,dim> &U, int)
{ 
  double t0 = this->timer->getTime();
  DistSVec<double,dim> Ubc(this->getVecInfo());
  Ubc = U;

  int its = 0;
  its = commonPart(U); //failure gives negative values
  projectStateOntoROB(U); // Lei Lei, 11 Oct 2016, for Rom projection
	if(its<0) return its; //failSafe

  TsDesc<dim>::setFailSafe(false);

  its = this->ns->solve(U);

  this->errorHandler->reduceError();
  this->data->resolveErrors();
  if(this->errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP]) return its;

	if(its < 0)
	{  
      //failSafe
    U = *EmbeddedTsDesc<dim>::UCopy;
    return its;
  }
  
  this->timer->addFluidSolutionTime(t0);
   
  int ierr = this->checkSolution(U);

	if(ierr > 0)
	{  
      //failSafe
    U = *EmbeddedTsDesc<dim>::UCopy;
    return (-ierr);
  }

  if(TsDesc<dim>::timeState->getData().typeIntegrator == ImplicitData::THREE_POINT_BDF &&
                 TsDesc<dim>::timeStepCalculation == TsData::ERRORESTIMATION )
    doErrorEstimation(U);

	if (this->modifiedGhidaglia) 
	{
    *this->bcData->getBoundaryStateHH() = embeddedU.hh();
    this->timeState->updateHH(embeddedU.hh());
  }
  
  return its;
}
//------------------------------------------------------------------------------
// External routines to solve Euler equations implicitly (called by NewtonSolver)
//------------------------------------------------------------------------------

// this function evaluates (Aw),t + F(w,x,v)
template<int dim>
void ImplicitEmbeddedTsDesc<dim>::computeFunction(int it, DistSVec<double,dim> &Q,
                                                  DistSVec<double,dim> &F)
{

  // Included for test with twilight zone problems (AM)
  // (Usually does nothing)
  this->domain->setExactBoundaryValues(Q, *this->X, this->ioData, 
				       this->currentTime + this->currentTimeStep,
				       this->spaceOp->getVarFcn());

	this->spaceOp->computeResidual(*this->X, *this->A, Q, 
											 *this->Wstarij, *this->Wstarji, *this->Wextij, this->distLSS,
											 this->linRecAtInterface, this->viscSecOrder, 
											 this->nodeTag, F, this->riemann, 
                                 this->riemannNormal, 1, this->ghostPoints);

	this->timeState->add_dAW_dt(it, *this->geoState, *this->A, Q, F,this->distLSS);

  this->spaceOp->applyBCsToResidual(Q, F,this->distLSS);

  this->domain->setExactBoundaryResidual(F, *this->X, this->ioData, 
					 this->currentTime + this->currentTimeStep,
					 this->spaceOp->getVarFcn());

 
	if (this->modifiedGhidaglia) 
	{
		*hhResidual = 0.0;
		this->domain->computeHHBoundaryTermResidual(*this->bcData,Q,*hhResidual, this->varFcn);

		this->timeState->add_dAW_dt_HH(-1, *this->geoState, *this->A,
												 *this->bcData->getBoundaryStateHH(),
												 *hhResidual);
  }

	
/*	int isubd;
#pragma omp parallel for
	for (isubd=0; isubd < this->domain->getNumLocSub(); ++isubd) 
	{
		int lnu = (*this->Wextij)(isubd).size();
		for(int k_=0; k_<lnu; ++k_)
		{
			if( !(*this->distLSS)(isubd).isActive(0.0,k_) ){				
				Vec3D XX = (*this->X)(isubd)[k_];
				if(XX[1] > 0.0 && XX[2] > 0.02) 
					fprintf(stdout, "%d, %f,%f,%f, %f\n", isubd, XX[0],XX[1],XX[2], (*this->Wextij)(isubd)[k_][4]);
			}
		} 
	}
*/
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitEmbeddedTsDesc<dim>::doErrorEstimation(DistSVec<double,dim> &U) 
{
  DistSVec<double,dim> *flux = new DistSVec<double,dim>(TsDesc<dim>::domain->getNodeDistInfo());

  this->spaceOp->computeResidual(*this->X, *this->A, this->timeState->getUn(), *this->Wstarij, *this->Wstarji, *this->Wextij, 
											this->distLSS, this->linRecAtInterface, this->viscSecOrder, this->nodeTag, *flux, this->riemann, 
                                 this->riemannNormal, 1, this->ghostPoints);

  this->timeState->calculateErrorEstiNorm(U, *flux); 

  delete flux;
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitEmbeddedTsDesc<dim>::recomputeFunction(DistSVec<double,dim> &Q,
                                            DistSVec<double,dim> &rhs)
{
  this->spaceOp->recomputeRHS(*this->X, Q, rhs);
}

//------------------------------------------------------------------------------

template<int dim>
int ImplicitEmbeddedTsDesc<dim>::checkFailSafe(DistSVec<double,dim>& U)
{
//  this->com->fprintf(stdout, "WARNING: At the moment CheckFailSafe is not supported by the embedded framework with an implicit time-integrator!\n");

  if (!this->failSafeNewton) return 0;

  if (!this->tag)
    this->tag = new DistSVec<bool,2>(this->getVecInfo());

  this->domain->checkFailSafe(this->varFcn, U, *this->tag);
  this->spaceOp->fix(*this->tag);

  return 1;

}

//------------------------------------------------------------------------------
template<int dim> 
void ImplicitEmbeddedTsDesc<dim>::resetFixesTag()
{

  this->spaceOp->resetTag();

}

//------------------------------------------------------------------------------
