#include <cmath>

#include <MatVecProd.h>
#include <IoData.h>

#include <BcDef.h>
#include <RecFcnDesc.h>
#include <FluxFcn.h>
#include <DistTimeState.h>
#include <DistGeoState.h>
#include <SpaceOperator.h>
#include <Domain.h>
#include <MemoryPool.h>

// Included (MB)
#include <FemEquationTermDesc.h>


// included for rom (lei lei, Sep 27 2016)
#include <VectorSet.h>
//------------------------------------------------------------------------------

template<int dim, int neq>
MatVecProdFD<dim, neq>::MatVecProdFD
(
  ImplicitData &data, 
  DistTimeState<dim> *ts,
  DistGeoState *gs, 
  SpaceOperator<dim> *spo, 
  Domain *domain, 
  IoData &ioData
)
  : geoState(gs), Qeps(domain->getNodeDistInfo()), Feps(domain->getNodeDistInfo())
                , Qepstmp(domain->getNodeDistInfo()), Fepstmp(domain->getNodeDistInfo())
                , Q(domain->getNodeDistInfo()), F(domain->getNodeDistInfo())
                , Ftmp(domain->getNodeDistInfo()), iod(&ioData)
{

  com = domain->getCommunicator();

  if (ts)
    timeState = new DistTimeState<dim>(*ts, false, ioData);
  else
    timeState = 0;

  spaceOp = new SpaceOperator<dim>(*spo, false);

  recFcnCon = 0;
  Rn = 0;

  if (data.mvp == ImplicitData::H1FD) {
    recFcnCon = new RecFcnConstant<dim>;
    spaceOp->setRecFcn(recFcnCon);
    if (data.type == ImplicitData::CRANK_NICOLSON) {
      Rn = new DistSVec<double,dim>(domain->getNodeDistInfo());
      timeState->setResidual(Rn);
    }
  }

  fdOrder = ioData.ts.implicit.fdOrder; 

  output = NULL;

}

//------------------------------------------------------------------------------

template<int dim, int neq>
MatVecProdFD<dim, neq>::~MatVecProdFD()
{ 

  if (spaceOp) delete spaceOp;
  if (timeState) delete timeState;
  recFcnCon = 0; // deleted by spaceOp
  Rn = 0; // deleted by timeState

}

template<int dim, int neq>
void MatVecProdFD<dim, neq>::attachHH(DistEmbeddedVec<double,dim>& v) {

  hhRes = new DistVec<double>(v.hhInfo());
  hhEps = new DistVec<double>(v.hhInfo());
  hhVal = new DistVec<double>(v.hhInfo());
}

//------------------------------------------------------------------------------

template<int dim, int neq>
void MatVecProdFD<dim, neq>::evaluate(int it, DistSVec<double,3> &x, DistVec<double> &cv,
				 DistSVec<double,dim> &q, DistSVec<double,dim> &f)
{
	// ****
  X = &x;
  ctrlVol = &cv;
  Qeps = q;

	if(recFcnCon && !this->isFSI) 
	{
    spaceOp->computeResidual(*X, *ctrlVol, Qeps, Feps, timeState);

		if (timeState) timeState->add_dAW_dt(it, *geoState, *ctrlVol, Qeps, Feps);

    spaceOp->applyBCsToResidual(Qeps, Feps);
  }
	else Feps = f;

	if(timeState && iod->ts.dualtimestepping == TsData::ON) 
	{
		timeState->add_dAW_dtau(it, *geoState, *ctrlVol, Qeps, Feps);

		spaceOp->applyBCsToResidual(Qeps, Feps);
	}

  Qeps.strip(Q);
  Feps.strip(F);
  
}

//------------------------------------------------------------------------------

template<int dim, int neq>
void MatVecProdFD<dim, neq>::evaluateHH(DistVec<double> &hhterm,
					DistVec<double> &bcVal ) {

  *hhRes = hhterm;
  *hhVal = bcVal;
}

//------------------------------------------------------------------------------

template<int dim, int neq>
void MatVecProdFD<dim, neq>::evaluate(DistExactRiemannSolver<dim> &riemann, int it, DistSVec<double,3> &x, DistVec<double> &cv,
                                 DistSVec<double,dim> &q, DistSVec<double,dim> &f)
{

  X = &x;
  ctrlVol = &cv;
  Qeps = q;

  if (recFcnCon && !this->isFSI) {
    spaceOp->computeResidual(&riemann, *X, *ctrlVol, Qeps, Feps, timeState);

    if (timeState)
      timeState->add_dAW_dt(it, *geoState, *ctrlVol, Qeps, Feps);

    spaceOp->applyBCsToResidual(Qeps, Feps);

  }
  else  {
    Feps = f;
  }

  if (timeState && iod->ts.dualtimestepping == TsData::ON) {
    timeState->add_dAW_dtau(it, *geoState, *ctrlVol, Qeps, Feps);
    spaceOp->applyBCsToResidual(Qeps, Feps);
  }

  Qeps.strip(Q);
  Feps.strip(F);

}

//------------------------------------------------------------------------------

template<int dim, int neq>
void MatVecProdFD<dim, neq>::evaluateRestrict(int it, DistSVec<double,3> &x,
		DistVec<double> &cv, DistSVec<double,dim> &q, DistSVec<double,dim> &f,
		RestrictionMapping<dim> & restrictionMapping,
                TsDesc<dim>* probDesc,
                int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &))
{

	std::vector<std::vector<int> > sampledLocNodes =
		restrictionMapping.getRestrictedToOriginLocNode() ;
  X = &x;
  ctrlVol = &cv;
  Qeps = q;
  Qeps.strip(Q);

  // clipping
  if (dim>5 && probDesc!=NULL && checkSolution!=NULL) {
    (probDesc->*checkSolution)(Qeps);
  }

  if (recFcnCon) {
    spaceOp->computeResidualRestrict(*X, *ctrlVol, Qeps, Feps, timeState, sampledLocNodes);

    if (timeState)
      timeState->add_dAW_dtRestrict(it, *geoState, *ctrlVol, Qeps, Feps, sampledLocNodes);

    spaceOp->applyBCsToResidual(Qeps, Feps);

  }
  else  {
    Feps = f;
  }
  Feps.strip(F);
  
}
//------------------------------------------------------------------------------

// Included (MB)
template<int dim, int neq>
void MatVecProdFD<dim, neq>::evaluateInviscid(int it, DistSVec<double,3> &x, DistVec<double> &cv,
				 DistSVec<double,dim> &q, DistSVec<double,dim> &f)
{

  X = &x;
  ctrlVol = &cv;
  Qeps = q;

  if (recFcnCon) {
    spaceOp->computeInviscidResidual(*X, *ctrlVol, Qeps, Feps, timeState);

    if (timeState)
      timeState->add_dAW_dt(it, *geoState, *ctrlVol, Qeps, Feps);

    spaceOp->applyBCsToResidual(Qeps, Feps);
  }
  else  {
    Feps = f;
  }

  if (timeState && iod->ts.dualtimestepping == TsData::ON) {
    timeState->add_dAW_dtau(it, *geoState, *ctrlVol, Qeps, Feps);
    spaceOp->applyBCsToResidual(Qeps, Feps);
  }
  
  Qeps.strip(Q);
  Feps.strip(F);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, int neq>
void MatVecProdFD<dim, neq>::evaluateViscous(int it, DistSVec<double,3> &x, DistVec<double> &cv,
				 DistSVec<double,dim> &q, DistSVec<double,dim> &f)
{

  X = &x;
  ctrlVol = &cv;
  Qeps = q;

  if (recFcnCon) {
    spaceOp->computeViscousResidual(*X, *ctrlVol, Qeps, Feps, timeState);

    spaceOp->applyBCsToResidual(Qeps, Feps);
  }
  else  {
    Feps = f;
  }

  Qeps.strip(Q);
  Feps.strip(F);

}
//------------------------------------------------------------------------------
template<int dim, int neq>
void MatVecProdFD<dim, neq>::apply(DistSVec<double,neq> &p, DistSVec<double,neq> &prod)
{


  double eps = computeEpsilon(Q, p);

// Included (MB)
  Qepstmp = Q + eps * p;

  Qepstmp.pad(Qeps);
  
  if (!this->isFSI)
    spaceOp->computeResidual(*X, *ctrlVol, Qeps, Feps, timeState);
  else
    spaceOp->computeResidual(*X,*ctrlVol, Qeps, *(this->fsi.Wtemp), *(this->fsi.Wtemp), *(this->fsi.Wtemp),
                             this->fsi.LSS, this->fsi.linRecAtInterface, this->fsi.viscSecOrder, *(this->fsi.fluidId),
                             Feps, this->fsi.riemann, this->fsi.Nriemann, 0, this->fsi.ghostPoints);

  if (timeState) {
    timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, Feps);
    if (iod->ts.dualtimestepping == TsData::ON)
      timeState->add_dAW_dtau(-1, *geoState, *ctrlVol, Qeps, Feps);
  }

  if (this->isFSI)
    spaceOp->applyBCsToResidual(Qeps, Feps, this->fsi.LSS);
  else
    spaceOp->applyBCsToResidual(Qeps, Feps);

  Feps.strip(Fepstmp);
  if (output) int status = output->writeBinaryVectorsToDiskRom(false, 0, 0, NULL, &Feps);

  if (fdOrder == 1) {

    prod = (1.0/eps) * (Fepstmp - F);
 
  }
  else if (fdOrder == 2) {

    Qepstmp = Q - eps * p;
    
    Qepstmp.pad(Qeps);
  
    if (!this->isFSI)
      spaceOp->computeResidual(*X, *ctrlVol, Qeps, Feps, timeState);
    else
      spaceOp->computeResidual(*X,*ctrlVol, Qeps, *(this->fsi.Wtemp),*(this->fsi.Wtemp),*(this->fsi.Wtemp),
                               this->fsi.LSS, this->fsi.linRecAtInterface, this->fsi.viscSecOrder, *(this->fsi.fluidId),
                               Feps, this->fsi.riemann, this->fsi.Nriemann, 0, this->fsi.ghostPoints);
 
    if (timeState) {
      timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, Feps);
      if (iod->ts.dualtimestepping == TsData::ON)
        timeState->add_dAW_dtau(-1, *geoState, *ctrlVol, Qeps, Feps);
    }

    if (this->isFSI)
      spaceOp->applyBCsToResidual(Qeps, Feps, this->fsi.LSS);
    else
      spaceOp->applyBCsToResidual(Qeps, Feps);

    Feps.strip(Ftmp);

    prod = (0.5/eps) * (Fepstmp - Ftmp);

    if (output) int status = output->writeBinaryVectorsToDiskRom(false, 0, 0, NULL, &Feps);

  }

// Original
/*

  //#define MVP_CHECK_ONE_EQ 5
#if MVP_CHECK_ONE_EQ
  int numLocSub = p.numLocSub();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    double (*qeps)[dim] = Qeps.subData(iSub);
    double (*q)[dim] = Q->subData(iSub);
    double (*lp)[dim] = p.subData(iSub);
    for (int i=0; i<p.subSize(iSub); ++i) {
      for (int j=0; j<dim; ++j)
	qeps[i][j] = q[i][j];
      qeps[i][MVP_CHECK_ONE_EQ] = q[i][MVP_CHECK_ONE_EQ] + eps * lp[i][MVP_CHECK_ONE_EQ];
    }
  }
#else
  Qepstmp = Q + eps * p;
  Qepstmp.pad(Qeps);
#endif

  if(Phi)
    spaceOp->computeResidual(*X, *ctrlVol, Qeps, *Phi, Feps, riemann, 0);
  else
    spaceOp->computeResidual(*X, *ctrlVol, Qeps, Feps, timeState);


  if (timeState) {
    timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, Feps);
    if (iod->ts.dualtimestepping == TsData::ON)
        timeState->add_dAW_dtau(-1, *geoState, *ctrlVol, Qeps, Feps);
  }

  spaceOp->applyBCsToResidual(Qeps, Feps);

  Feps.strip(Fepstmp);
  
  prod = (1.0/eps) * (Fepstmp - F);
*/

}
//------------------------------------------------------------------------------
template<int dim, int neq>
void MatVecProdFD<dim, neq>::applyRestrict(DistSVec<double,neq> &p,
		DistSVec<double,neq> &prod, RestrictionMapping<neq> & restrictionMapping,
                TsDesc<dim>* probDesc,
                int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &))
{
  std::vector<std::vector<int> > sampledLocNodes = restrictionMapping.getRestrictedToOriginLocNode() ;
  DistSVec<double, neq> pRestricted(restrictionMapping.restrictedDistInfo());
  restrictionMapping.restriction(p,pRestricted);

  double eps = computeEpsilon(Q, pRestricted);

// Included (MB)
  Qepstmp = Q + eps * p;

  Qepstmp.pad(Qeps);

  // clipping
  if (dim>5 && probDesc!=NULL && checkSolution!=NULL) {
    (probDesc->*checkSolution)(Qeps);
  }

  spaceOp->computeResidualRestrict(*X, *ctrlVol, Qeps, Feps, timeState, sampledLocNodes);

  if (timeState)
    timeState->add_dAW_dtRestrict(-1, *geoState, *ctrlVol, Qeps, Feps, sampledLocNodes);

  spaceOp->applyBCsToResidual(Qeps, Feps);

  Feps.strip(Fepstmp);

  if (fdOrder == 1) {

    prod = (1.0/eps) * (Fepstmp - F);
 
  }
  else if (fdOrder == 2) {

    Qepstmp = Q - eps * p;
    
    Qepstmp.pad(Qeps);

    spaceOp->computeResidualRestrict(*X, *ctrlVol, Qeps, Feps, timeState, sampledLocNodes);

    if (timeState)
      timeState->add_dAW_dtRestrict(-1, *geoState, *ctrlVol, Qeps, Feps, sampledLocNodes);

    spaceOp->applyBCsToResidual(Qeps, Feps);

    Feps.strip(Ftmp);

    prod = (0.5/eps) * (Fepstmp - Ftmp);

  }

}

//-----------------------------------------------------

  // ROMs minimize the residual, so the weighting of the residual becomes very important.
  // These functions allow for Jacobians of residuals with non-constant weights.

template<int dim, int neq>
void MatVecProdFD<dim, neq>::evaluateWeighted(int it, DistSVec<double,3> &x, DistVec<double> &cv,
                            DistSVec<double,dim> &q, DistSVec<double,dim> &f, VarFcn *varFcn)
{
  X = &x;
  ctrlVol = &cv;
  Qeps = q;

  if (recFcnCon && !this->isFSI) {
    spaceOp->computeResidual(*X, *ctrlVol, Qeps, Feps, timeState);
    if (timeState) {
      timeState->add_dAW_dt(it, *geoState, *ctrlVol, Qeps, Feps);
    }
    spaceOp->applyBCsToResidual(Qeps, Feps);
  } else {
    Feps = f;
  }

  if (timeState && iod->ts.dualtimestepping == TsData::ON) {
    timeState->add_dAW_dtau(it, *geoState, *ctrlVol, Qeps, Feps);
    spaceOp->applyBCsToResidual(Qeps, Feps);
  }

  DistSVec<double, dim>* componentwiseScalingVec = varFcn->computeScalingVec(*iod,Qeps);
  Feps *= *componentwiseScalingVec;
  delete componentwiseScalingVec;

  Qeps.strip(Q);
  Feps.strip(F);

}

template<int dim, int neq>
void MatVecProdFD<dim, neq>::applyWeighted(DistSVec<double,neq> &p, DistSVec<double,neq> &prod, VarFcn *varFcn)
{
  double eps = computeEpsilon(Q, p);

  Qepstmp = Q + eps * p;

  Qepstmp.pad(Qeps);

  if (!this->isFSI)
    spaceOp->computeResidual(*X, *ctrlVol, Qeps, Feps, timeState);
  else
    spaceOp->computeResidual(*X,*ctrlVol, Qeps, *(this->fsi.Wtemp),*(this->fsi.Wtemp),*(this->fsi.Wtemp),
                             this->fsi.LSS, this->fsi.linRecAtInterface, this->fsi.viscSecOrder, *(this->fsi.fluidId),
                             Feps, this->fsi.riemann, this->fsi.Nriemann, 0, this->fsi.ghostPoints);

  if (timeState) {
    timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, Feps);
    if (iod->ts.dualtimestepping == TsData::ON)
      timeState->add_dAW_dtau(-1, *geoState, *ctrlVol, Qeps, Feps);
  }

  if (this->isFSI)
    spaceOp->applyBCsToResidual(Qeps, Feps, this->fsi.LSS);
  else
    spaceOp->applyBCsToResidual(Qeps, Feps);

  DistSVec<double, dim>* componentwiseScalingVec = varFcn->computeScalingVec(*iod,Qeps);
  Feps *= *componentwiseScalingVec;
  delete componentwiseScalingVec;

  Feps.strip(Fepstmp);
  if (output) int status = output->writeBinaryVectorsToDiskRom(false, 0, 0, NULL, &Feps);

  if (fdOrder == 1) {
    prod = (1.0/eps) * (Fepstmp - F);
  } else if (fdOrder == 2) {

    Qepstmp = Q - eps * p;
    
    Qepstmp.pad(Qeps);
  
    if (!this->isFSI)
      spaceOp->computeResidual(*X, *ctrlVol, Qeps, Feps, timeState);
    else
      spaceOp->computeResidual(*X,*ctrlVol, Qeps, *(this->fsi.Wtemp),*(this->fsi.Wtemp),*(this->fsi.Wtemp),
                               this->fsi.LSS, this->fsi.linRecAtInterface, this->fsi.viscSecOrder, *(this->fsi.fluidId),
                               Feps, this->fsi.riemann, this->fsi.Nriemann, 0, this->fsi.ghostPoints);

 
    if (timeState) {
      timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, Feps);
      if (iod->ts.dualtimestepping == TsData::ON)
        timeState->add_dAW_dtau(-1, *geoState, *ctrlVol, Qeps, Feps);
    }

    if (this->isFSI)
      spaceOp->applyBCsToResidual(Qeps, Feps, this->fsi.LSS);
    else
      spaceOp->applyBCsToResidual(Qeps, Feps);

    DistSVec<double, dim>* componentwiseScalingVec = varFcn->computeScalingVec(*iod,Qeps);
    Feps *= *componentwiseScalingVec;
    delete componentwiseScalingVec;

    Feps.strip(Ftmp);

    prod = (0.5/eps) * (Fepstmp - Ftmp);

    if (output) int status = output->writeBinaryVectorsToDiskRom(false, 0, 0, NULL, &Feps);

  }
}

//--------------------------------------------------

template<int dim, int neq>
void MatVecProdFD<dim, neq>::evaluateWeightedRestrict(int it, DistSVec<double,3> &x,
                                  DistVec<double> &cv,
                                  DistSVec<double,dim> &q, DistSVec<double,dim> &f,
                                  RestrictionMapping<dim> & restrictionMapping,
                                  VarFcn *varFcn, TsDesc<dim>* probDesc,
                                  int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &)){

  std::vector<std::vector<int> > sampledLocNodes = restrictionMapping.getRestrictedToOriginLocNode() ;
  X = &x;
  ctrlVol = &cv;
  Qeps = q;

  // clipping
  if (dim>5 && probDesc!=NULL && checkSolution!=NULL) {
    (probDesc->*checkSolution)(Qeps);
  }

  if (recFcnCon) {
    spaceOp->computeResidualRestrict(*X, *ctrlVol, Qeps, Feps, timeState, sampledLocNodes);

    if (timeState)
      timeState->add_dAW_dtRestrict(it, *geoState, *ctrlVol, Qeps, Feps, sampledLocNodes);

    spaceOp->applyBCsToResidual(Qeps, Feps);

  }
  else  {
    Feps = f;
  }

  DistSVec<double, dim>* componentwiseScalingVec = varFcn->computeScalingVec(*iod,Qeps);
  Feps *= *componentwiseScalingVec;
  delete componentwiseScalingVec;

  Qeps.strip(Q);
  Feps.strip(F);

}

//--------------------------------------------------

template<int dim, int neq>
void MatVecProdFD<dim, neq>::applyWeightedRestrict(DistSVec<double,neq> &p,
                                         DistSVec<double,neq> &prod,
                                         RestrictionMapping<neq> & restrictionMapping,
                                         VarFcn *varFcn,
                                         TsDesc<dim>* probDesc,
                                         int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &))
{

  std::vector<std::vector<int> > sampledLocNodes = restrictionMapping.getRestrictedToOriginLocNode() ;
  DistSVec<double, neq> pRestricted(restrictionMapping.restrictedDistInfo());
  restrictionMapping.restriction(p,pRestricted);

  double eps = computeEpsilon(Q, pRestricted);

  Qepstmp = Q + eps * p;

  Qepstmp.pad(Qeps);

  // clipping
  if (dim>5 && probDesc!=NULL && checkSolution!=NULL) {
    (probDesc->*checkSolution)(Qeps);
  }

  spaceOp->computeResidualRestrict(*X, *ctrlVol, Qeps, Feps, timeState, sampledLocNodes);

  if (timeState)
    timeState->add_dAW_dtRestrict(-1, *geoState, *ctrlVol, Qeps, Feps, sampledLocNodes);

  spaceOp->applyBCsToResidual(Qeps, Feps);

  DistSVec<double, dim>* componentwiseScalingVec = varFcn->computeScalingVec(*iod,Qeps);
  Feps *= *componentwiseScalingVec;
  delete componentwiseScalingVec;

  Feps.strip(Fepstmp);

  if (fdOrder == 1) {
    prod = (1.0/eps) * (Fepstmp - F);
  }
  else if (fdOrder == 2) {
    Qepstmp = Q - eps * p;
    Qepstmp.pad(Qeps);
    spaceOp->computeResidualRestrict(*X, *ctrlVol, Qeps, Feps, timeState, sampledLocNodes);

    if (timeState)
      timeState->add_dAW_dtRestrict(-1, *geoState, *ctrlVol, Qeps, Feps, sampledLocNodes);

    spaceOp->applyBCsToResidual(Qeps, Feps);

    DistSVec<double, dim>* componentwiseScalingVec = varFcn->computeScalingVec(*iod,Qeps);
    Feps *= *componentwiseScalingVec;
    delete componentwiseScalingVec;

    Feps.strip(Ftmp);
    prod = (0.5/eps) * (Fepstmp - Ftmp);
  }

}

//--------------------------------------------------

template<int dim, int neq>
void MatVecProdFD<dim,neq>::apply(DistEmbeddedVec<double,neq> & p, DistEmbeddedVec<double,neq> & prod) 
{
	// ***
  double eps = computeEpsilon(Q, p.real());

// Included (MB)
  Qepstmp = Q + eps * p.real();

  Qepstmp.pad(Qeps);
  
	if(p.hasHHBoundaryTerm()) 
	{
    *hhEps = *hhVal + eps*p.hh();
    *spaceOp->getDistBcData()->getBoundaryStateHH() = *hhEps;
  }

  if (!this->isFSI)
    spaceOp->computeResidual(*X, *ctrlVol, Qeps, Feps, timeState);
  else
	  spaceOp->computeResidual(*X, *ctrlVol, Qeps, 
										*(this->fsi.Wtemp), *(this->fsi.Wtemp), *(this->fsi.Wtemp), this->fsi.LSS, 
										this->fsi.linRecAtInterface, this->fsi.viscSecOrder, 
										*(this->fsi.fluidId), Feps, this->fsi.riemann, 
										this->fsi.Nriemann, 0, this->fsi.ghostPoints);

	if(p.hasHHBoundaryTerm()) 
	{
    *hhEps = 0.0;

		spaceOp->getDomain()->computeHHBoundaryTermResidual(*spaceOp->getDistBcData(),Qeps,*hhEps, spaceOp->getVarFcn());
       
		timeState->add_dAW_dt_HH(-1, *geoState, *ctrlVol,*spaceOp->getDistBcData()->getBoundaryStateHH(), *hhEps);
  }

	if(timeState) 
	{
		timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, Feps, this->fsi.LSS);
		if(iod->ts.dualtimestepping == TsData::ON) timeState->add_dAW_dtau(-1, *geoState, *ctrlVol, Qeps, Feps, this->fsi.LSS);
  }

	if(this->isFSI) spaceOp->applyBCsToResidual(Qeps, Feps, this->fsi.LSS);
	else            spaceOp->applyBCsToResidual(Qeps, Feps);

  Feps.strip(Fepstmp);

	if(fdOrder == 1) 
	{
    prod.real() = (1.0/eps) * (Fepstmp - F);

		if(p.hasHHBoundaryTerm())
		{
      *hhEps = (1.0/eps) * (*hhEps - *hhRes);
      prod.setHH(*hhEps);
    }
 
  }
	else if(fdOrder == 2) 
	{
    Qepstmp = Q - eps * p.real();
    
    Qepstmp.pad(Qeps);
  
		if(!this->isFSI)	spaceOp->computeResidual(*X, *ctrlVol, Qeps, Feps, timeState);
		else              spaceOp->computeResidual(*X,*ctrlVol, Qeps, *(this->fsi.Wtemp), *(this->fsi.Wtemp),*(this->fsi.Wtemp),
																 this->fsi.LSS, this->fsi.linRecAtInterface, 
																 this->fsi.viscSecOrder, *(this->fsi.fluidId),
                               Feps, this->fsi.riemann, this->fsi.Nriemann, 0, this->fsi.ghostPoints);
 
		if(timeState) 
		{
			timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, Feps, this->fsi.LSS);

			if (iod->ts.dualtimestepping == TsData::ON) timeState->add_dAW_dtau(-1, *geoState, *ctrlVol, Qeps, Feps, this->fsi.LSS);
    }

		if(p.hasHHBoundaryTerm()) 
		{
      prod.setHH(*hhEps);
      *hhEps = *hhVal - eps*p.hh();
      *spaceOp->getDistBcData()->getBoundaryStateHH() = *hhEps;
      *hhEps = 0.0;
			
			spaceOp->getDomain()->computeHHBoundaryTermResidual(*spaceOp->getDistBcData(),Qeps,*hhEps, spaceOp->getVarFcn());
       
			timeState->add_dAW_dt_HH(-1, *geoState, *ctrlVol,*spaceOp->getDistBcData()->getBoundaryStateHH(), *hhEps);

      prod.hh() = (0.5/eps)*(prod.hh()-*hhEps);
    }

		if(this->isFSI) spaceOp->applyBCsToResidual(Qeps, Feps, this->fsi.LSS);
		else            spaceOp->applyBCsToResidual(Qeps, Feps);

    Feps.strip(Ftmp);

    prod.real() = (0.5/eps) * (Fepstmp - Ftmp);
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, int neq>
void MatVecProdFD<dim, neq>::applyInviscid(DistSVec<double,neq> &p, DistSVec<double,neq> &prod)
{

  double eps = computeEpsilon(Q, p);

  Qepstmp = Q + eps * p;

  Qepstmp.pad(Qeps);

  spaceOp->computeInviscidResidual(*X, *ctrlVol, Qeps, Feps, timeState);

  if (timeState) {
    timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, Feps);
    if (iod->ts.dualtimestepping == TsData::ON)
      timeState->add_dAW_dtau(-1, *geoState, *ctrlVol, Qeps, Feps);
  }

  Q.pad(Qeps);

  spaceOp->applyBCsToResidual(Qeps, Feps);

  Feps.strip(Fepstmp);

  if (output) int status = output->writeBinaryVectorsToDiskRom(false, 0, 0, NULL, &Feps);
 
  if (fdOrder == 1) {

    spaceOp->computeInviscidResidual(*X, *ctrlVol, Qeps, Feps, timeState);

    if (timeState) {
      timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, Feps);
      if (iod->ts.dualtimestepping == TsData::ON)
        timeState->add_dAW_dtau(-1, *geoState, *ctrlVol, Qeps, Feps);
    }

    spaceOp->applyBCsToResidual(Qeps, Feps);

    Feps.strip(Ftmp);
    
    prod += (1.0/eps) * (Fepstmp - Ftmp);
 
  }
  else if (fdOrder == 2) {

    Qepstmp = Q - eps * p;

    Qepstmp.pad(Qeps);

    spaceOp->computeInviscidResidual(*X, *ctrlVol, Qeps, Feps, timeState);

    if (timeState) {
      timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, Feps);
      if (iod->ts.dualtimestepping == TsData::ON)
        timeState->add_dAW_dtau(-1, *geoState, *ctrlVol, Qeps, Feps);
    }

    Q.pad(Qeps);

    spaceOp->applyBCsToResidual(Qeps, Feps);

    Feps.strip(Ftmp);

    prod += (0.5/eps) * (Fepstmp - Ftmp);

    if (output) int status = output->writeBinaryVectorsToDiskRom(false, 0, 0, NULL, &Feps);

  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, int neq>
void MatVecProdFD<dim, neq>::applyViscous(DistSVec<double,neq> &p, DistSVec<double,neq> &prod)
{

  double eps = computeEpsilon(Q, p);

  Qepstmp = Q + eps * p;

  Qepstmp.pad(Qeps);

  spaceOp->computeViscousResidual(*X, *ctrlVol, Qeps, Feps, timeState);

  spaceOp->applyBCsToResidual(Qeps, Feps);

  Feps.strip(Fepstmp);
  
  if (output) int status = output->writeBinaryVectorsToDiskRom(false, 0, 0, NULL, &Feps);

  if (fdOrder == 1) {

    Q.pad(Qeps);

    spaceOp->computeViscousResidual(*X, *ctrlVol, Qeps, Feps, timeState);

    spaceOp->applyBCsToResidual(Qeps, Feps);

    Feps.strip(Ftmp);

    prod += (1.0/eps) * (Fepstmp - Ftmp);
 
  }
  else if (fdOrder == 2) {

    Qepstmp = Q - eps * p;

    Qepstmp.pad(Qeps);

    spaceOp->computeViscousResidual(*X, *ctrlVol, Qeps, Feps, timeState);

    spaceOp->applyBCsToResidual(Qeps, Feps);

    Feps.strip(Ftmp);

    prod += (0.5/eps) * (Fepstmp - Ftmp);

    if (output) int status = output->writeBinaryVectorsToDiskRom(false, 0, 0, NULL, &Feps);

  }

}

//------------------------------------------------------------------------------

template<int dim, int neq>
double MatVecProdFD<dim, neq>::computeEpsilon(DistSVec<double,neq> &U, DistSVec<double,neq> &p)
{

  int iSub, size = 0;
  double eps0 = 1.e-6;

  double norm = sqrt(p*p);
  if (norm < 1.0e-14)
    return eps0;

  const DistInfo &distInfo = U.info();

  double *alleps = reinterpret_cast<double *>(alloca(sizeof(double) * distInfo.numGlobSub));

  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) alleps[iSub] = 0.0;

#pragma omp parallel for reduction(+: size)
  for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

    int locOffset = distInfo.subOffset[iSub];

    int locsize = 0;
    double loceps = 0.0;

    for (int i=0; i<distInfo.subLen[iSub]; ++i) {

      if (distInfo.masterFlag[locOffset+i]) {
	for (int j=0; j<neq; ++j) {
	  ++locsize;
	  loceps += eps0*fabs(U[locOffset+i][j]) + eps0;
	}
      }

    }

    size += locsize;
    alleps[distInfo.locSubToGlobSub[iSub]] = loceps;

  }

  this->com->globalSum(1, &size);
  this->com->globalSum(distInfo.numGlobSub, alleps);

  double eps = 0.0;
  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) eps += alleps[iSub];

  eps /= double(size) * norm;
 
  return eps;

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, int neq>
void MatVecProdFD<dim, neq>::rstSpaceOp(IoData & ioData, VarFcn *varFcn, SpaceOperator<dim> *spo, bool typeAlloc, SpaceOperator<dim> *spofd)
{

  // UH (09/10) -> Check for memory leak?
  // No FluxFcn pointer is deleted.

  spaceOp->rstFluxFcn(ioData);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
MatVecProdH1<dim,Scalar,neq>::MatVecProdH1(DistTimeState<dim> *ts, SpaceOperator<dim> *spo,
					   Domain *domain) : 
  DistMat<Scalar,neq>(domain), timeState(ts), spaceOp(spo)
{

#ifdef _OPENMP 
  this->numLocSub = DistMat<Scalar,neq>::numLocSub; //BUG omp
#endif

  A = new MvpMat<Scalar,neq>*[this->numLocSub];

  double size = 0.0;

#pragma omp parallel for reduction (+: size)
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {

    A[iSub] = this->subDomain[iSub]->template createMaskMatVecProd<Scalar,neq>();

    size += double(A[iSub]->numNonZeroBlocks()*neq*neq*sizeof(Scalar)) / (1024.*1024.);

  }

  areHHTermsActive = false;

  this->com->globalSum(1, &size);
  
  this->com->printf(2, "Memory required for matvec with H1 (dim=%d): %3.2f MB\n", neq, size);

}
//------------------------------------------------------------------------------
                                                                                                  
template<int dim, class Scalar, int neq>
MatVecProdH1<dim,Scalar,neq>::MatVecProdH1(DistTimeState<dim> *ts, SpaceOperator<dim> *spo,
                                           Domain *domain, IoData &ioData) :
  DistMat<Scalar,neq>(domain), timeState(ts), spaceOp(spo)
{

#ifdef _OPENMP
  this->numLocSub = DistMat<Scalar,neq>::numLocSub; //BUG omp
#endif

  A = new MvpMat<Scalar,neq>*[this->numLocSub];

  double size = 0.0;

#pragma omp parallel for reduction (+: size)
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {
    A[iSub] = this->subDomain[iSub]->template createMaskMatVecProd<Scalar,neq>();
    size += double(A[iSub]->numNonZeroBlocks()*neq*neq*sizeof(Scalar)) / (1024.*1024.);
  }

  areHHTermsActive = false;

  this->com->globalSum(1, &size);

  this->com->printf(2, "Memory required for matvec with H1 (dim=%d): %3.2f MB\n", neq, size);

}

template<int dim, class Scalar,  int neq>
void MatVecProdH1<dim,Scalar, neq>::attachHH(DistEmbeddedVec<double,dim>& v) {

  areHHTermsActive = true;

#pragma omp parallel for reduction (+: size)
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {

    A[iSub]->enableHHTerms(v.hh()(iSub).size());
  }

  hhVal = new DistVec<double>(v.hhInfo());
}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
MatVecProdH1<dim,Scalar,neq>::~MatVecProdH1()
{

  if (A) {
#pragma omp parallel for
    for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
      if (A[iSub]) delete A[iSub];

    delete [] A;
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
DistMat<Scalar,neq> &MatVecProdH1<dim,Scalar,neq>::operator= (const Scalar x)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub)
    *A[iSub] = x;

  return *this;

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH1<dim,Scalar,neq>::exportMemory(MemoryPool *mp)
{

  if (!mp) return;

  for (int iSub = 0; iSub < this->numLocSub; ++iSub)
    mp->set(A[iSub]->numNonZeroBlocks() * neq*neq * sizeof(Scalar), 
	    reinterpret_cast<void *>(A[iSub]->data()));

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH1<dim,Scalar,neq>::evaluateHH(DistVec<double> &hhterm,
					      DistVec<double> &bcVal ) {

  *hhVal = bcVal;
}

template<int dim, class Scalar, int neq>
void MatVecProdH1<dim,Scalar,neq>::evaluate(int it, DistSVec<double,3> &X, DistVec<double> &ctrlVol, 
					    DistSVec<double,dim> &Q, DistSVec<double,dim> &F)
{

  if (!this->isFSI)
    spaceOp->computeJacobian(X, ctrlVol, Q, *this, timeState);
  else
    spaceOp->computeJacobian(X,ctrlVol, Q, this->fsi.LSS, *(this->fsi.fluidId), 
                             this->fsi.riemann, this->fsi.Nriemann,
                             this->fsi.ghostPoints, *this,timeState);

  if (timeState)
    timeState->addToJacobian(ctrlVol, *this, Q);

  if (this->isFSI)
    spaceOp->applyBCsToJacobian(Q, *this, this->fsi.LSS);
  else
    spaceOp->applyBCsToJacobian(Q, *this);

  if (areHHTermsActive) {

#pragma omp parallel for
    for (int iSub = 0; iSub < this->numLocSub; ++iSub) {

      this->subDomain[iSub]->computeJacobianFiniteVolumeTermHH(spaceOp->getFluxFcn(),
							       (*spaceOp->getDistBcData())(iSub) ,
							       (*spaceOp->getGeoState())(iSub),
							       ctrlVol(iSub), Q(iSub), *A[iSub], spaceOp->getVarFcn());
    }
    
    this->timeState->addToHHJacobian(ctrlVol, *this, *hhVal);
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH1<dim,Scalar,neq>::evaluateRestrict(int it, DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                            DistSVec<double,dim> &Q, DistSVec<double,dim> &F, 
                                            RestrictionMapping<dim> & restrictionMapping,
                                            TsDesc<dim>* probDesc,
                                            int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &)) {
  //TODO
  MatVecProdH1<dim,Scalar,neq>::evaluate(it, X, ctrlVol, Q, F);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH1<dim,Scalar,neq>::evaluate(DistExactRiemannSolver<dim> &riemann,
                                            int it, DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                            DistSVec<double,dim> &Q, DistSVec<double,dim> &F)
{

  spaceOp->computeJacobian(&riemann, X, ctrlVol, Q, *this, timeState);

  if (timeState)
    timeState->addToJacobian(ctrlVol, *this, Q);

  spaceOp->applyBCsToJacobian(Q, *this);

  if (areHHTermsActive) {

#pragma omp parallel for
    for (int iSub = 0; iSub < this->numLocSub; ++iSub) {

      this->subDomain[iSub]->computeJacobianFiniteVolumeTermHH(spaceOp->getFluxFcn(),
							       (*spaceOp->getDistBcData())(iSub),
							       (*spaceOp->getGeoState())(iSub),
							       ctrlVol(iSub), Q(iSub), *A[iSub], spaceOp->getVarFcn());
    }
    this->timeState->addToHHJacobian(ctrlVol, *this, *hhVal);
    
  }
}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH1<dim,Scalar,neq>::evaluateViscous(int it, DistSVec<double,3> &X,
                                   DistVec<double> &cv)  {

  // Compute the Jacobian of viscous terms
  spaceOp->computeViscousJacobian(X, cv, *this);

}

//------------------------------------------------------------------------------
// note: this can be done in another way (but less efficient) !!
// (1) compute off-diag products (2) assemble (3) compute diag products (->redundancy)

template<int dim, class Scalar, int neq>
void MatVecProdH1<dim,Scalar,neq>::apply(DistSVec<double,neq> &p, DistSVec<double,neq> &prod)
{

  int iSub;
  
#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub) {
    this->subDomain[iSub]->computeMatVecProdH1(p.getMasterFlag(iSub), *A[iSub],
					 p(iSub), prod(iSub));
    this->subDomain[iSub]->sndData(*this->vecPat, prod.subData(iSub));
  }

  this->vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    this->subDomain[iSub]->addRcvData(*this->vecPat, prod.subData(iSub));

}

template<int dim, class Scalar, int neq>
void MatVecProdH1<dim,Scalar,neq>::applyTranspose(DistSVec<double,neq> &p, DistSVec<double,neq> &prod)
{

  int iSub;
  
#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub) {
    this->subDomain[iSub]->computeMatVecProdH1transpose(p.getMasterFlag(iSub), *A[iSub],
                                                        p(iSub), prod(iSub));
    this->subDomain[iSub]->sndData(*this->vecPat, prod.subData(iSub));
  }

  this->vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    this->subDomain[iSub]->addRcvData(*this->vecPat, prod.subData(iSub));

}

template<int dim, class Scalar, int neq>
void MatVecProdH1<dim,Scalar,neq>::applyRestrict(DistSVec<double,neq> &p, DistSVec<double,neq> &prod,
                                         RestrictionMapping<neq> & restrictionMapping,
                                         TsDesc<dim>* probDesc,
                                         int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &))
{ //TODO
  MatVecProdH1<dim,Scalar,neq>::apply(p, prod);
}

//------------------------------------------------------------------------------


template<int dim, class Scalar, int neq>
void MatVecProdH1<dim,Scalar,neq>::apply(DistEmbeddedVec<double,neq> &p, DistEmbeddedVec<double,neq> &prod)
{

  int iSub;
  
#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub) {
    this->subDomain[iSub]->computeMatVecProdH1(p.real().getMasterFlag(iSub), *A[iSub],
					       p.real()(iSub), prod.real()(iSub), 
                                               p.ghost()(iSub), prod.ghost()(iSub) );

    if (p.hasHHBoundaryTerm()) {
      if (!prod.hasHHBoundaryTerm()) {
    
        prod.setHH(p.hh());
      }
      prod.hh() = 0.0;
      this->subDomain[iSub]->
	computeMatVecProdH1FarFieldHH(p.real().getMasterFlag(iSub),
				      *A[iSub],p.real()(iSub), prod.real()(iSub), 
				      p.hh()(iSub), prod.hh()(iSub));
    }

    this->subDomain[iSub]->sndData(*this->vecPat, prod.real().subData(iSub));
  }


  this->vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    this->subDomain[iSub]->addRcvData(*this->vecPat, prod.real().subData(iSub));

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub) {
    this->subDomain[iSub]->sndData(*this->vecPat, prod.ghost().subData(iSub));
  }
  this->vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    this->subDomain[iSub]->addRcvData(*this->vecPat, prod.ghost().subData(iSub));
}

template<int dim, class Scalar, int neq>
void MatVecProdH1<dim,Scalar,neq>::clearGhost()
{
#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {
    A[iSub]->clearGhost();
  }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar, int neq>
void MatVecProdH1<dim,Scalar,neq>::rstSpaceOp(IoData & ioData, VarFcn *varFcn, SpaceOperator<dim> *spo, bool typeAlloc, SpaceOperator<dim> *spofd)
{

  // UH (09/10) -> Check for memory leak
  // No FluxFcn pointer is deleted.

  spaceOp = spo;

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
MatVecProdH2<dim,Scalar,neq>::MatVecProdH2
(
  IoData &ioData, VarFcn *varFcn, DistTimeState<dim> *ts,
  SpaceOperator<dim> *spo, Domain *domain, DistGeoState *gs
) :
  MatVecProd<dim,neq>(),
  DistMat<Scalar,dim>(domain),
  A(0),
  aij(domain->getEdgeDistInfo()), aji(domain->getEdgeDistInfo()),
  bij(domain->getEdgeDistInfo()), bji(domain->getEdgeDistInfo()),
  betaij(domain->getEdgeDistInfo()), betaji(domain->getEdgeDistInfo()),
  timeState(ts),
  spaceOp(0),
  fluxFcn(0),
  X(0),
  ctrlVol(0),
  Q(0),
  F(0)
  , R(0)
  , RFD(0)
  , vProd(0)
{

#ifdef _OPENMP 
  this->numLocSub = DistMat<Scalar,dim>::numLocSub; //BUG omp
#endif

  A = new MvpMat<Scalar,dim>*[this->numLocSub];

  double size = 0.0;
  double coefsize = double(4*aij.size()*dim*sizeof(double)) / (1024.*1024.);

  // allocate for viscous flux jacobian term
  bool nsFlag = false;
  spaceOp = new SpaceOperator<dim>(*spo, false);

// Included (MB*)
  if ((ioData.eqs.type == EquationsData::NAVIER_STOKES) && (ioData.bc.wall.integration != BcsWallData::WALL_FUNCTION))
    nsFlag = true;

/*
//
// Original
//
//  if (ioData.eqs.type == EquationsData::NAVIER_STOKES)  {
//    R = new MatVecProdH1<dim, Scalar ,dim>(ts, spo, domain);
//    nsFlag = true;
//  }
//  else
//    R = 0;
*/

#pragma omp parallel for reduction (+: size)
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {

    A[iSub] = this->subDomain[iSub]->template createMaskMatVecProd<Scalar,dim>(nsFlag);

    size += double(A[iSub]->numNonZeroBlocks()*dim*dim*sizeof(Scalar)) / (1024.*1024.);

  }

  this->com->globalSum(1, &size);
  this->com->globalSum(1, &coefsize);
  
  this->com->printf(2, "Memory required for matvec with H2: ");
  this->com->printf(2, "%3.2f+%3.2f=%3.2f MB\n", size, coefsize, size+coefsize);
  

// Included (MB)
  fluxFcn = new FluxFcn*[BC_MAX_CODE - BC_MIN_CODE + 1]; 
  fluxFcn -= BC_MIN_CODE;
  if(BC_MAX_CODE-BC_MIN_CODE+1 < 22)
    fprintf(stderr,"Be prepared to see a segmentation fault shortly...\n");
  fluxFcn[BC_SYMMETRY] = new FluxFcn(0,BC_SYMMETRY,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_MASSFLOW_OUTLET_MOVING] = new FluxFcn(0,BC_MASSFLOW_OUTLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_MASSFLOW_OUTLET_FIXED] = new FluxFcn(0,BC_MASSFLOW_OUTLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_MASSFLOW_INLET_MOVING] = new FluxFcn(0,BC_MASSFLOW_INLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_MASSFLOW_INLET_FIXED] = new FluxFcn(0,BC_MASSFLOW_INLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_DIRECTSTATE_OUTLET_MOVING] = new FluxFcn(0,BC_DIRECTSTATE_OUTLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_DIRECTSTATE_OUTLET_FIXED] = new FluxFcn(0,BC_DIRECTSTATE_OUTLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_DIRECTSTATE_INLET_MOVING] = new FluxFcn(0,BC_DIRECTSTATE_INLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_DIRECTSTATE_INLET_FIXED] = new FluxFcn(0,BC_DIRECTSTATE_INLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_OUTLET_MOVING] = new FluxFcn(0,BC_OUTLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_OUTLET_FIXED] = new FluxFcn(0,BC_OUTLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_INLET_MOVING] = new FluxFcn(0,BC_INLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_INLET_FIXED] = new FluxFcn(0,BC_INLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_POROUS_WALL_MOVING] = new FluxFcn(0,BC_POROUS_WALL_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_POROUS_WALL_FIXED] = new FluxFcn(0,BC_POROUS_WALL_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ADIABATIC_WALL_MOVING] = new FluxFcn(0,BC_ADIABATIC_WALL_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ADIABATIC_WALL_FIXED] = new FluxFcn(0,BC_ADIABATIC_WALL_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_SLIP_WALL_MOVING] = new FluxFcn(0,BC_SLIP_WALL_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_SLIP_WALL_FIXED] = new FluxFcn(0,BC_SLIP_WALL_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ISOTHERMAL_WALL_MOVING] = new FluxFcn(0,BC_ISOTHERMAL_WALL_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ISOTHERMAL_WALL_FIXED] = new FluxFcn(0,BC_ISOTHERMAL_WALL_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_INTERNAL] = new FluxFcn(0,BC_INTERNAL,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  spaceOp->setFluxFcn(fluxFcn); //TODO: should avoid doing this!

  // Included (MB)
  // (09/10) UH >> Modification
  // The flag viscJacContrib was previously controlled by the sensitivity module.
  // It governs the computation of the Jacobian for the viscous term
  // (laminar flux and, possibly, turbulent flux).
  // The case viscJacContrib == 1 gives an "exact" computation
  // (exact for the laminar part and, possibly, approximate for the turbulent flux).
  // The case viscJacContrib == 2 gives a finite difference approximation.
  // It is not activated but kept for reference.
  int viscJacContrib = 1;
  if (ioData.eqs.type == EquationsData::NAVIER_STOKES)
  {
    vProd = new DistSVec<double,neq>(domain->getNodeDistInfo());
    if (viscJacContrib == 1)
    {
      R = new MatVecProdH1<dim, Scalar, neq>(ts, spo, domain);
    }
    else if ((viscJacContrib == 2) && (gs))
    {
      RFD = new MatVecProdFD<dim,neq>(ioData.ts.implicit, ts, gs, spo, domain, ioData);
    }
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
MatVecProdH2<dim,Scalar,neq>::~MatVecProdH2()
{ 

  if (spaceOp) delete spaceOp;
  fluxFcn = 0; // deleted by spaceOperator.

  if (A) {
#pragma omp parallel for
    for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
      if (A[iSub]) delete A[iSub];

    delete [] A;
  }

  if (vProd)
    delete vProd;
  vProd = 0;

  if (R)
    delete R;
  R = 0;

  if (RFD)
    delete RFD;
  RFD = 0;

}

//------------------------------------------------------------------------------
template<int dim, class Scalar, int neq>
DistMat<Scalar,dim> &MatVecProdH2<dim,Scalar,neq>::operator= (const Scalar x)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub)
    *A[iSub] = x;

  return *this;

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH2<dim,Scalar,neq>::evaluate(int it, DistSVec<double,3> &x, DistVec<double> &cv, 
					DistSVec<double,dim> &q, DistSVec<double,dim> &f)
{

// Included (MB)
  evaluateInviscid(it, x, cv, q, f);
  evaluateViscous( it, x, cv, q, f);

// Original
/*
  X = &x;
  ctrlVol = &cv;
  Q = &q;

  spaceOp->computeH2(*X, *ctrlVol, *Q, *this, aij, aji, bij, bji);

  if (timeState)
    timeState->addToH2(*ctrlVol, *Q, *this);
  


  // compute viscous flux jacobian
  if (R)  {
    spaceOp->applyBCsToH2Jacobian(*Q, *this);
    R->evaluateViscous(it, *X, *ctrlVol);
    spaceOp->applyBCsToJacobian(*Q, *R);
  }
*/

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH2<dim,Scalar,neq>::evaluateRestrict(int it, DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                            DistSVec<double,dim> &Q, DistSVec<double,dim> &F, 
                                            RestrictionMapping<dim> & restrictionMapping,
                                            TsDesc<dim>* probDesc,
                                            int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &)) {

  // TODO
  MatVecProdH2<dim,Scalar,neq>::evaluate(it, X, ctrlVol, Q, F);

}


//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar, int neq>
void MatVecProdH2<dim,Scalar,neq>::evaluateInviscid(int it, DistSVec<double,3> &x, DistVec<double> &cv, 
					DistSVec<double,dim> &q, DistSVec<double,dim> &f)
{

  X = &x;
  ctrlVol = &cv;
  Q = &q;

  if (!this->isFSI){
    spaceOp->computeH2(*X, *ctrlVol, *Q, *this, aij, aji, bij, bji);
  }else{
    spaceOp->computeH2(*X, *ctrlVol, *Q, this->fsi.LSS, *(this->fsi.fluidId), 
		       this->fsi.riemann, this->fsi.Nriemann,
		       this->fsi.ghostPoints, *this, aij, aji, bij, bji, betaij, betaji);
  }

  if (timeState)
    timeState->addToH2(*ctrlVol, *Q, *this);

  spaceOp->applyBCsToH2Jacobian(*Q, *this);


}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar, int neq>
void MatVecProdH2<dim,Scalar,neq>::evaluateViscous(int it, DistSVec<double,3> &x, DistVec<double> &cv, 
					DistSVec<double,dim> &q, DistSVec<double,dim> &f)
{

  // compute viscous flux jacobian
  if (R)  {
    R->evaluateViscous(it, *X, *ctrlVol);
    spaceOp->applyBCsToJacobian(*Q, *R);
  }

  if (RFD) {
    F = &f;
    RFD->evaluateViscous(it, *X, *ctrlVol, *Q, *F);
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH2<dim,Scalar,neq>::evaluate(int it, DistSVec<double,3> &x, 
                               DistVec<double> &cv, DistSVec<double,dim> &q, 
                               DistSVec<double,dim> &F, Scalar shift)
{

  X = &x;
  ctrlVol = &cv;
  Q = &q;

  spaceOp->computeH2(*X, *ctrlVol, *Q, *this, aij, aji, bij, bji);

  if (timeState) {
     switch (it)  {
       //case for the construction of the POD
       case 0:
         timeState->addToH2(*ctrlVol, *Q, *this, shift, 1.0);
         break;   
       case 1:
         timeState->addToH2(*ctrlVol, *Q, *this, Scalar(2.0), 1.0);
         break;
       case 2:
         timeState->addToH2(*ctrlVol, *Q, *this, Scalar(2.0), -1.0);
         break;
       case 3:
         timeState->addToH2(*ctrlVol, *Q, *this, Scalar(3.0), 2.0);
         break;
       case 4:
         timeState->addToH2(*ctrlVol, *Q, *this, Scalar(4.0), -2.0);
         break;
       case 5:
         timeState->addToH2(*ctrlVol, *Q, *this, Scalar(2.0), -2.0);
         break;
       case 6:
         timeState->addToH2(*ctrlVol, *Q, *this, (Scalar(4.0)*shift+Scalar(2.0))/(shift*(shift+Scalar(1.0))), 2.0);
         break;
       case 7:
         timeState->addToH2(*ctrlVol, *Q, *this, Scalar(8.0/3.0), 2.0);
         break;
       case 8:
         timeState->addToH2(*ctrlVol, *Q, *this, Scalar(16.0/3.0), 2.0);
         break;
       case 9:
         timeState->addToH2(*ctrlVol, *Q, *this, Scalar(10.0/3.0), 2.0);
         break;

    }
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH2<dim,Scalar,neq>::evaluate2(int it, DistSVec<double,3> &x, DistVec<double> &cv, 
                                         DistSVec<double,dim> &q, DistSVec<double,dim> &F)
{

  X = &x;
  ctrlVol = &cv;
  Q = &q;

  spaceOp->computeH2(*X, *ctrlVol, *Q, *this, aij, aji, bij, bji);

  if (timeState) {
    timeState->addToH2Minus(*ctrlVol, *Q, *this);
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH2<dim,Scalar,neq>::evalH(int it, DistSVec<double,3> &x,
                               DistVec<double> &cv, DistSVec<double,dim> &q)  
{

  X = &x;
  ctrlVol = &cv;
  Q = &q;

  spaceOp->computeH2(*X, *ctrlVol, *Q, *this, aij, aji, bij, bji);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH2<dim,Scalar,neq>::apply(DistSVec<double,neq> &p, DistSVec<double,neq> &prod)
{

// Original
/*
  spaceOp->applyH2(*X, *ctrlVol, *Q, *this, aij, aji, bij, bji, p, prod);

  if (R)  {
    DistSVec<double, dim> vProd(p);
    vProd = 0.0;
    R->apply(p, vProd);
    prod += vProd;
  }
*/

  //std::cout << "$$$$$ IN MatVecProc H2 applyXXXX\n";

  Multiplier<dim,neq,Scalar,double> Operator;
  Operator.Apply(spaceOp, *X, *ctrlVol, *Q, *this, aij, aji, bij, bji, p, prod,
                 R, RFD, vProd);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH2<dim,Scalar,neq>::applyRestrict(DistSVec<double,neq> &p, DistSVec<double,neq> &prod,
                                         RestrictionMapping<neq> & restrictionMapping,
                                         TsDesc<dim>* probDesc,
                                         int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &))
{
  //TODO
  MatVecProdH2<dim,Scalar,neq>::apply(p, prod);
}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH2<dim,Scalar,neq>::applyTranspose(DistSVec<double,neq> &p, DistSVec<double,neq> &prod)
{

// Original
/*
  spaceOp->applyH2(*X, *ctrlVol, *Q, *this, aij, aji, bij, bji, p, prod);

  if (R)  {
    DistSVec<double, dim> vProd(p);
    vProd = 0.0;
    R->apply(p, vProd);
    prod += vProd;
  }
*/

  Multiplier<dim,neq,Scalar,double> Operator;
  Operator.ApplyTranspose(spaceOp, *X, *ctrlVol, *Q, *this, aij, aji, bij, bji, p, prod,
                          R, RFD, vProd);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH2<dim,Scalar,neq>::apply(DistEmbeddedVec<double,dim> &p, DistEmbeddedVec<double,dim> &prod)
{

  //std::cout << "$$$$$ IN MatVecProc EMB H2 applyXXXX \n";

  spaceOp->applyH2(*X, *ctrlVol, *Q, 
		   this->fsi.LSS, *(this->fsi.fluidId), 
		   this->fsi.linRecAtInterface, this->fsi.viscSecOrder,
		   this->fsi.riemann, this->fsi.Nriemann,
		   this->fsi.ghostPoints, 
		   *this, aij, aji, bij, bji, betaij, betaji,
		   p.real(), prod.real());

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
template<int dd, int nn, class Scalar1, class Scalar2>
void MatVecProdH2<dim,Scalar,neq>::Multiplier<dd,nn,Scalar1,Scalar2>::Apply
(
  SpaceOperator<dd> *spaceOp
  , DistSVec<double,3> &X
  , DistVec<double> &ctrlVol
  , DistSVec<double,dd> &U
  , DistMat<Scalar1,dd> &H2
  , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
  , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
  , DistSVec<Scalar2,nn> &p, DistSVec<Scalar2,nn> &prod
  , MatVecProdH1<dd, Scalar1, nn> *R
  , MatVecProdFD<dd, nn> *RFD
  , DistSVec<Scalar2, nn> *vProd
)
{

  if (nn > dd)
  {
    std::cout << "\n !!! Apply Not Implemented for dd = " << dd;
    std::cout << " nn = " << nn << std::endl;
    exit(-1);
  }

  DistSVec<Scalar2,dd> pExt(p.info());
  pExt = (Scalar2) 0;

  int numLocSub = p.numLocSub();
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
  {
    Scalar2 (*locp)[nn] = p.subData(iSub);
    Scalar2 (*locExt)[dd] = pExt.subData(iSub);
    for (int i = 0; i < p.subSize(iSub); ++i)
    {
      for (int jj = 0; jj < nn; ++jj)
        locExt[i][jj] = locp[i][jj];
    }
  } // for (int iSub = 0; iSub < numLocSub; ++iSub)

  spaceOp->applyH2(X, ctrlVol, U, H2, aij, aji, bij, bji, pExt, pExt);

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
  {
    Scalar2 (*locp)[nn] = prod.subData(iSub);
    Scalar2 (*locExt)[dd] = pExt.subData(iSub);
    for (int i = 0; i < p.subSize(iSub); ++i)
    {
      for (int jj = 0; jj < nn; ++jj)
        locp[i][jj] = locExt[i][jj];
    }
  } // for (int iSub = 0; iSub < numLocSub; ++iSub)

  if (R)
  {
    *vProd = (Scalar2) 0;
    R->apply(p, *vProd);
    prod += *vProd;
  }
  else if (RFD)
  {
    *vProd = (Scalar2) 0;
    RFD->applyViscous(p, *vProd);
    prod += *vProd;
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
template<int dd, int nn, class Scalar1, class Scalar2>
void MatVecProdH2<dim,Scalar,neq>::Multiplier<dd,nn,Scalar1,Scalar2>::ApplyTranspose
(
  SpaceOperator<dd> *spaceOp
  , DistSVec<double,3> &X
  , DistVec<double> &ctrlVol
  , DistSVec<double,dd> &U
  , DistMat<Scalar1,dd> &H2
  , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
  , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
  , DistSVec<Scalar2,nn> &p, DistSVec<Scalar2,nn> &prod
  , MatVecProdH1<dd, Scalar1, nn> *R
  , MatVecProdFD<dd, nn> *RFD
  , DistSVec<Scalar2, nn> *vProd
)
{

  if (nn > dd)
  {
    std::cout << "\n !!! Apply Not Implemented for dd = " << dd;
    std::cout << " nn = " << nn << std::endl;
    exit(-1);
  }

  DistSVec<Scalar2,dd> pExt(p.info());
  pExt = (Scalar2) 0;

  int numLocSub = p.numLocSub();
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
  {
    Scalar2 (*locp)[nn] = p.subData(iSub);
    Scalar2 (*locExt)[dd] = pExt.subData(iSub);
    for (int i = 0; i < p.subSize(iSub); ++i)
    {
      for (int jj = 0; jj < nn; ++jj)
        locExt[i][jj] = locp[i][jj];
    }
  } // for (int iSub = 0; iSub < numLocSub; ++iSub)

  spaceOp->applyH2transposeNew(X, ctrlVol, U, H2, aij, aji, bij, bji, pExt, pExt);

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
  {
    Scalar2 (*locp)[nn] = prod.subData(iSub);
    Scalar2 (*locExt)[dd] = pExt.subData(iSub);
    for (int i = 0; i < p.subSize(iSub); ++i)
    {
      for (int jj = 0; jj < nn; ++jj)
        locp[i][jj] = locExt[i][jj];
    }
  } // for (int iSub = 0; iSub < numLocSub; ++iSub)

  if (R)
  {
    *vProd = (Scalar2) 0;
    R->apply(p, *vProd);
    prod += *vProd;
  }
  else if (RFD)
  {
    *vProd = (Scalar2) 0;
    RFD->applyViscous(p, *vProd);
    prod += *vProd;
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
template<int dd, class Scalar1, class Scalar2>
void MatVecProdH2<dim,Scalar,neq>::Multiplier<dd,dd,Scalar1,Scalar2>::Apply
(
  SpaceOperator<dd> *spaceOp
  , DistSVec<double,3> &X
  , DistVec<double> &ctrlVol
  , DistSVec<double,dd> &U
  , DistMat<Scalar1,dd> &H2
  , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
  , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
  , DistSVec<Scalar2,dd> &p, DistSVec<Scalar2,dd> &prod
  , MatVecProdH1<dd, Scalar1, dd> *R
  , MatVecProdFD<dd, dd> *RFD
  , DistSVec<Scalar2, dd> *vProd
)
{

  spaceOp->applyH2(X, ctrlVol, U, H2, aij, aji, bij, bji, p, prod);

  if (R)
  {
    *vProd = (Scalar2) 0;
    R->apply(p, *vProd);
    prod += *vProd;
  }
  else if (RFD)
  {
    *vProd = (Scalar2) 0;
    RFD->applyViscous(p, *vProd);
    prod += *vProd;
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
template<int dd, class Scalar1, class Scalar2>
void MatVecProdH2<dim,Scalar,neq>::Multiplier<dd,dd,Scalar1,Scalar2>::ApplyTranspose
(
  SpaceOperator<dd> *spaceOp
  , DistSVec<double,3> &X
  , DistVec<double> &ctrlVol
  , DistSVec<double,dd> &U
  , DistMat<Scalar1,dd> &H2
  , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
  , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
  , DistSVec<Scalar2,dd> &p, DistSVec<Scalar2,dd> &prod
  , MatVecProdH1<dd, Scalar1, dd> *R
  , MatVecProdFD<dd, dd> *RFD
  , DistSVec<Scalar2, dd> *vProd
)
{

  spaceOp->applyH2transposeNew(X, ctrlVol, U, H2, aij, aji, bij, bji, p, prod);

  if (R)
  {
//    std::cout << "\n !!! R is being added !!\n";
    *vProd = (Scalar2) 0;
    R->applyTranspose(p, *vProd);
    prod += *vProd;
  }
  else if (RFD) //TODO: YC: have not checked this works
  {
    *vProd = (Scalar2) 0;
    RFD->applyViscous(p, *vProd);
    prod += *vProd;
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH2<dim,Scalar,neq>::apply(DistSVec<bcomp,neq> &p,
                DistSVec<bcomp,neq> &prod)
{

  //std::cout << "$$$$$ IN MatVecProc H2 C applyapply\n";

  Multiplier<dim,neq,Scalar,bcomp> Operator;
  Operator.Apply(spaceOp, *X, *ctrlVol, *Q, *this, aij, aji, bij, bji, p, prod,
                 0, 0, 0);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH2<dim,Scalar,neq>::applyT(DistSVec<double,neq> &p,
        DistSVec<double,neq> &prod)
{

  Multiplier<dim,neq,Scalar,double> Operator;
  Operator.ApplyT(spaceOp, *X, *ctrlVol, *Q, *this, aij, aji, bij, bji, p, prod);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void MatVecProdH2<dim,Scalar,neq>::applyT(DistSVec<bcomp,neq> &p,
        DistSVec<bcomp,neq> &prod)
{

  Multiplier<dim,neq,Scalar,bcomp> Operator;
  Operator.ApplyT(spaceOp, *X, *ctrlVol, *Q, *this, aij, aji, bij, bji, p, prod);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
template<int dd, int nn, class Scalar1, class Scalar2>
void MatVecProdH2<dim,Scalar,neq>::Multiplier<dd,nn,Scalar1,Scalar2>::ApplyT
(
  SpaceOperator<dd> *spaceOp
  , DistSVec<double,3> &X
  , DistVec<double> &ctrlVol
  , DistSVec<double,dd> &U
  , DistMat<Scalar1,dd> &H2
  , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
  , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
  , DistSVec<Scalar2,nn> &p, DistSVec<Scalar2,nn> &prod
)
{

  std::cout << "\n !!! ApplyT Not Implemented !!\n";
  exit(-1);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
template<int dd, class Scalar1, class Scalar2>
void MatVecProdH2<dim,Scalar,neq>::Multiplier<dd,dd,Scalar1,Scalar2>::ApplyT
(
  SpaceOperator<dd> *spaceOp
  , DistSVec<double,3> &X
  , DistVec<double> &ctrlVol
  , DistSVec<double,dd> &U
  , DistMat<Scalar1,dd> &H2
  , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
  , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
  , DistSVec<Scalar2,dd> &p, DistSVec<Scalar2,dd> &prod
)
{

  spaceOp->applyH2T(X, ctrlVol, U, H2, aij, aji, bij, bji, p, prod);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar, int neq>
void MatVecProdH2<dim,Scalar,neq>::rstSpaceOp
(
  IoData & ioData, VarFcn *varFcn, SpaceOperator<dim> *spo, 
  bool typeAlloc, SpaceOperator<dim> *spofd
)
{

  if (dim != neq)
  {
    // UH (09/10) This function is only called from the sensitivity module.
    // The sensitivity module assumes a strong turbulence model coupling
    // (i.e. dim == neq).
    this->com->fprintf(stderr, "\n *** MatVecProdH2<dim,Scalar,neq>::rstSpaceOp");
    this->com->fprintf(stderr, " is not verified for weakly coupled systems.\n\n");
    exit(1);
  }

  // UH (09/10) -> Check for memory leak
  // No FluxFcn pointer is deleted.

  fluxFcn = new FluxFcn*[BC_MAX_CODE - BC_MIN_CODE + 1]; 
  fluxFcn -= BC_MIN_CODE;
  if(BC_MAX_CODE-BC_MIN_CODE+1 < 22)
    fprintf(stderr,"Be prepared to see a segmentation fault shortly...\n");

  fluxFcn[BC_SYMMETRY] = new FluxFcn(0,BC_SYMMETRY,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_MASSFLOW_OUTLET_MOVING] = new FluxFcn(0,BC_MASSFLOW_OUTLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_MASSFLOW_OUTLET_FIXED] = new FluxFcn(0,BC_MASSFLOW_OUTLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_MASSFLOW_INLET_MOVING] = new FluxFcn(0,BC_MASSFLOW_INLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_MASSFLOW_INLET_FIXED] = new FluxFcn(0,BC_MASSFLOW_INLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_DIRECTSTATE_OUTLET_MOVING] = new FluxFcn(0,BC_DIRECTSTATE_OUTLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_DIRECTSTATE_OUTLET_FIXED] = new FluxFcn(0,BC_DIRECTSTATE_OUTLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_DIRECTSTATE_INLET_MOVING] = new FluxFcn(0,BC_DIRECTSTATE_INLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_DIRECTSTATE_INLET_FIXED] = new FluxFcn(0,BC_DIRECTSTATE_INLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_OUTLET_MOVING] = new FluxFcn(0,BC_OUTLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_OUTLET_FIXED] = new FluxFcn(0,BC_OUTLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_INLET_MOVING] = new FluxFcn(0,BC_INLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_INLET_FIXED] = new FluxFcn(0,BC_INLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_POROUS_WALL_MOVING] = new FluxFcn(0,BC_POROUS_WALL_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_POROUS_WALL_FIXED] = new FluxFcn(0,BC_POROUS_WALL_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ADIABATIC_WALL_MOVING] = new FluxFcn(0,BC_ADIABATIC_WALL_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ADIABATIC_WALL_FIXED] = new FluxFcn(0,BC_ADIABATIC_WALL_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_SLIP_WALL_MOVING] = new FluxFcn(0,BC_SLIP_WALL_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_SLIP_WALL_FIXED] = new FluxFcn(0,BC_SLIP_WALL_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ISOTHERMAL_WALL_MOVING] = new FluxFcn(0,BC_ISOTHERMAL_WALL_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ISOTHERMAL_WALL_FIXED] = new FluxFcn(0,BC_ISOTHERMAL_WALL_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_INTERNAL] = new FluxFcn(0,BC_INTERNAL,ioData,varFcn,FluxFcnBase::PRIMITIVE);

  spaceOp->setFluxFcn(fluxFcn);

}

//----------------------------------------------------------------------------//
//                MatVecProd for Multiphase Euler equations                   //
//----------------------------------------------------------------------------//
// Finite Difference method for Multiphase Euler equations
template<int dim, int dimLS>
MatVecProdFDMultiPhase<dim,dimLS>::MatVecProdFDMultiPhase(
                                   DistTimeState<dim> *ts, DistGeoState *gs,
                                   MultiPhaseSpaceOperator<dim,dimLS> * spo,
                                   DistExactRiemannSolver<dim> * rsolver,
                                   FluidSelector *fs, Domain *dom,
				   IoData& ioData) :
  MatVecProdMultiPhase<dim,dimLS>(ts,spo,rsolver,fs), geoState(gs),
  Qeps(dom->getNodeDistInfo()), Feps(dom->getNodeDistInfo()),
  Q(dom->getNodeDistInfo()), F(dom->getNodeDistInfo()), iod(&ioData)
{

  X = 0;
  ctrlVol = 0;
  Phi = 0;
  com = dom->getCommunicator();

  fdOrder = ioData.ts.implicit.fdOrder;
}

//----------------------------------------------------------------------------//

template<int dim, int dimLS>
MatVecProdFDMultiPhase<dim,dimLS>::~MatVecProdFDMultiPhase()
{

  this->spaceOp = 0;
  this->timeState = 0;
  this->riemann = 0;
  this->fluidSelector = 0;

  geoState = 0;
  X = 0;
  ctrlVol = 0;
  Phi = 0;
  com = 0;

}

template<int dim, int dimLS>
void MatVecProdFDMultiPhase<dim, dimLS>::attachHH(DistEmbeddedVec<double,dim>& v) {

  hhRes = new DistVec<double>(v.hhInfo());
  hhEps = new DistVec<double>(v.hhInfo());
  hhVal = new DistVec<double>(v.hhInfo());
}


template<int dim, int dimLS>
void MatVecProdFDMultiPhase<dim, dimLS>::evaluateHH(DistVec<double> &hhterm,
			   		DistVec<double> &bcVal ) {

  *hhRes = hhterm;
  *hhVal = bcVal;
}


//----------------------------------------------------------------------------//

template<int dim, int dimLS>
void MatVecProdFDMultiPhase<dim, dimLS>::evaluate(int it,
                                 DistSVec<double,3> &x, DistVec<double> &cv,
                                 DistSVec<double,dim> &q, DistSVec<double,dimLS> &phi,
                                 DistSVec<double,dim> &f)
{
  
  X = &x;
  ctrlVol = &cv;
  Qeps = q;
  Phi = &phi;

/*  this->spaceOp->computeResidual(*X, *ctrlVol, Qeps, *Phi, *this->fluidSelector, Feps, this->riemann, 1);

  if (this->timeState)
    this->timeState->add_dAW_dt(it, *geoState, *ctrlVol, Qeps, Feps);

  this->spaceOp->applyBCsToResidual(Qeps, Feps);
*/
  Q = Qeps;
  F = f;//Feps;
  
  if (this->timeState && iod->ts.dualtimestepping == TsData::ON) {
    this->timeState->add_dAW_dtau(it, *geoState, *ctrlVol, Q, F);
    this->spaceOp->applyBCsToResidual(Q, F);
  }
  
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MatVecProdFDMultiPhase<dim, dimLS>::apply(DistSVec<double,dim> &p,
                                               DistSVec<double,dim> &prod)
{

  double eps = computeEpsilon(Q, p);

// Included (MB)
  Qeps = Q + eps * p;

  //com->fprintf(stderr,"Computed eps = %e; p.p = %e; Q.Q = %e\n",eps,p*p, Q*Q);  
  if (!this->isFSI)
    this->spaceOp->computeResidual(*X, *ctrlVol, Qeps, *Phi, *this->fluidSelector, Feps, this->riemann,this->timeState, -1);
  else
    this->spaceOp->computeResidual(*X, *ctrlVol, Qeps, *this->fsi.Wtemp, *this->fsi.Wtemp,
                                   this->fsi.LSS, this->fsi.linRecAtInterface, this->fsi.viscSecOrder,
                                   this->riemann, this->fsi.Nriemann, *Phi, 
                                   *this->fluidSelector, Feps, 0, this->fsi.ghostPoints);    

  if (this->timeState)
    this->timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, Feps);
    if (iod->ts.dualtimestepping == TsData::ON)
      this->timeState->add_dAW_dtau(-1, *geoState, *ctrlVol, Qeps, Feps);

  this->spaceOp->applyBCsToResidual(Qeps, Feps);

  if (fdOrder == 1) {
    
    prod = (1.0/eps) * (Feps - F);
 
  }
  else if (fdOrder == 2) {

    Qeps = Q - eps * p;
    
    if (!this->isFSI)
      this->spaceOp->computeResidual(*X, *ctrlVol, Qeps, *Phi, *this->fluidSelector, F, this->riemann,this->timeState, -1);
    else
      this->spaceOp->computeResidual(*X, *ctrlVol, Qeps, *this->fsi.Wtemp, *this->fsi.Wtemp,
                                   this->fsi.LSS, this->fsi.linRecAtInterface, this->fsi.viscSecOrder,
                                   this->riemann,this->fsi.Nriemann, *Phi, 
                                   *this->fluidSelector, F, 0, this->fsi.ghostPoints);    

    if (this->timeState)
      this->timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, F);
      if (iod->ts.dualtimestepping == TsData::ON)
        this->timeState->add_dAW_dtau(-1, *geoState, *ctrlVol, Qeps, F);

    this->spaceOp->applyBCsToResidual(Qeps, F);

    prod = (0.5/eps) * (Feps - F);

  }

}

template<int dim, int dimLS>
void MatVecProdFDMultiPhase<dim, dimLS>::apply(DistEmbeddedVec<double,dim> &p,
                                               DistEmbeddedVec<double,dim> &prod)
{

  double eps = computeEpsilon(Q, p.real());

// Included (MB)
  Qeps = Q + eps * p.real();

  //com->fprintf(stderr,"Computed eps = %e; p.p = %e; Q.Q = %e\n",eps,p*p, Q*Q);  
  if (p.hasHHBoundaryTerm()) {
    *hhEps = *hhVal + eps*p.hh();
    *this->spaceOp->getDistBcData()->getBoundaryStateHH() = *hhEps;
  }

  if (!this->isFSI)
    this->spaceOp->computeResidual(*X, *ctrlVol, Qeps, *Phi, *this->fluidSelector, Feps, this->riemann,this->timeState, -1);
  else
    this->spaceOp->computeResidual(*X, *ctrlVol, Qeps, *this->fsi.Wtemp, *this->fsi.Wtemp,
                                   this->fsi.LSS, this->fsi.linRecAtInterface, this->fsi.viscSecOrder,
                                   this->riemann, this->fsi.Nriemann, *Phi, 
                                   *this->fluidSelector, Feps, 0, this->fsi.ghostPoints);    

  //std::cout << iod->ts.dualtimestepping << std::endl;
  if (this->timeState)
    this->timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, Feps);
    if (iod->ts.dualtimestepping == TsData::ON)
      this->timeState->add_dAW_dtau(-1, *geoState, *ctrlVol, Qeps, Feps);

  this->spaceOp->applyBCsToResidual(Qeps, Feps);

  if (p.hasHHBoundaryTerm()) {
    *hhEps = 0.0;
    this->spaceOp->getDomain()->
      computeHHBoundaryTermResidual(*this->spaceOp->getDistBcData(),Qeps,*hhEps, this->spaceOp->getVarFcn());
       
    this->timeState->add_dAW_dt_HH(-1, *geoState, *ctrlVol,*this->spaceOp->getDistBcData()->getBoundaryStateHH()
			     , *hhEps);
  }



  if (fdOrder == 1) {
    
    prod.real() = (1.0/eps) * (Feps - F);
 
    if (p.hasHHBoundaryTerm()) {
      *hhEps = (1.0/eps) * (*hhEps - *hhRes);
      prod.setHH(*hhEps);
    }

  }
  else if (fdOrder == 2) {

    Qeps = Q - eps * p.real();
    
    if (!this->isFSI)
      this->spaceOp->computeResidual(*X, *ctrlVol, Qeps, *Phi, *this->fluidSelector, F, this->riemann,this->timeState, -1);
    else
      this->spaceOp->computeResidual(*X, *ctrlVol, Qeps, *this->fsi.Wtemp, *this->fsi.Wtemp,
                                   this->fsi.LSS, this->fsi.linRecAtInterface, this->fsi.viscSecOrder,
                                   this->riemann,this->fsi.Nriemann, *Phi, 
                                   *this->fluidSelector, F, 0, this->fsi.ghostPoints);    

    if (this->timeState)
      this->timeState->add_dAW_dt(-1, *geoState, *ctrlVol, Qeps, F);
      if (iod->ts.dualtimestepping == TsData::ON)
        this->timeState->add_dAW_dtau(-1, *geoState, *ctrlVol, Qeps, F);

    if (p.hasHHBoundaryTerm()) {

      prod.setHH(*hhEps);
      *hhEps = *hhVal - eps*p.hh();
      *this->spaceOp->getDistBcData()->getBoundaryStateHH() = *hhEps;
      *hhEps = 0.0;
      this->spaceOp->getDomain()->
        computeHHBoundaryTermResidual(*this->spaceOp->getDistBcData(),Qeps,*hhEps, this->spaceOp->getVarFcn());
       
      this->timeState->add_dAW_dt_HH(-1, *geoState, *ctrlVol,*this->spaceOp->getDistBcData()->getBoundaryStateHH()
  			       , *hhEps);
      prod.hh() = (0.5/eps)*(prod.hh()-*hhEps);
    }

    this->spaceOp->applyBCsToResidual(Qeps, F);

    prod.real() = (0.5/eps) * (Feps - F);

  }

}


//------------------------------------------------------------------------------

template<int dim, int dimLS>
double MatVecProdFDMultiPhase<dim, dimLS>::computeEpsilon(DistSVec<double,dim> &U, DistSVec<double,dim> &p)
{

  int iSub, size = 0;
  double eps0 = 1.e-8;

  const DistInfo &distInfo = U.info();

  //double *alleps = reinterpret_cast<double *>(alloca(sizeof(double) * distInfo.numGlobSub));
 
  double* alleps = new double[distInfo.numGlobSub];

  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) alleps[iSub] = 0.0;

#pragma omp parallel for reduction(+: size)
  for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

    int locOffset = distInfo.subOffset[iSub];

    int locsize = 0;
    double loceps = 0.0;

    for (int i=0; i<distInfo.subLen[iSub]; ++i) {

      if (distInfo.masterFlag[locOffset+i]) {
	for (int j=0; j<dim; ++j) {
	  ++locsize;
	  loceps += eps0*fabs(U[locOffset+i][j]) + eps0;
	}
      }

    }

    size += locsize;
    alleps[distInfo.locSubToGlobSub[iSub]] = loceps;

  }

  this->com->globalSum(1, &size);
  this->com->globalSum(distInfo.numGlobSub, alleps);

  double eps = 0.0;
  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) eps += alleps[iSub];

  double norm = sqrt(p*p);

  if (norm > 1.e-14) eps /= double(size) * norm;
  else eps = eps0;

  delete [] alleps;

  return eps;
 
/*
  DistSVec<double,dim> Up(U);
  Up *= eps0;

  double totsum = Up.norm(); 
  
  double norm = sqrt(p*p);

  if (norm > 1.0e-14)
    totsum /= (double)(Up.size()*dim)*norm;
  else
    totsum = eps0;

  return totsum;
*/
}

//------------------------------------------------------------------------------
// H1 MatVecProd for Multiphase Euler equations

template<int dim, int dimLS>
MatVecProdH1MultiPhase<dim,dimLS>::MatVecProdH1MultiPhase(DistTimeState<dim> *ts,
                                   MultiPhaseSpaceOperator<dim,dimLS> *spo,
                                   DistExactRiemannSolver<dim> *rsolver,
                                   FluidSelector *fs, Domain *dom) :
  MatVecProdMultiPhase<dim,dimLS>(ts,spo,rsolver,fs), DistMat<double,dim>(dom)
{

  A = new MvpMat<double,dim>*[this->numLocSub];

  double size = 0.0;

#pragma omp parallel for reduction (+: size)
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {
    A[iSub] = this->subDomain[iSub]->template createMaskMatVecProd<double,dim>();
    size += double(A[iSub]->numNonZeroBlocks()*dim*dim*sizeof(double)) / (1024.*1024.);
  }

  this->com->globalSum(1, &size);

  areHHTermsActive = false;


  this->com->printf(2, "Memory required for MultiPhaseMatVec with H1 (dim=%d): %3.2f MB\n", dim, size);


}

template<int dim, int dimLS>
void MatVecProdH1MultiPhase<dim,dimLS>::attachHH(DistEmbeddedVec<double,dim>& v) {

  areHHTermsActive = true;

#pragma omp parallel for reduction (+: size)
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {

    A[iSub]->enableHHTerms(v.hh()(iSub).size());
  }

  hhVal = new DistVec<double>(v.hhInfo());
}


//------------------------------------------------------------------------------

template<int dim, int dimLS>
MatVecProdH1MultiPhase<dim,dimLS>::~MatVecProdH1MultiPhase()
{

  if (A) {
#pragma omp parallel for
    for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
      if (A[iSub]) delete A[iSub];

    delete [] A;
  }

  this->timeState = 0;
  this->spaceOp   = 0;
  this->riemann   = 0;
  this->fluidSelector = 0;
  
}
//------------------------------------------------------------------------------

template<int dim, int dimLS>
DistMat<double,dim> &MatVecProdH1MultiPhase<dim,dimLS>::operator= (const double x)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub)
    *A[iSub] = x;

  return *this;

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MatVecProdH1MultiPhase<dim,dimLS>::exportMemory(MemoryPool *mp)
{

  if (!mp) return;

  for (int iSub = 0; iSub < this->numLocSub; ++iSub)
    mp->set(A[iSub]->numNonZeroBlocks() * dim*dim * sizeof(double), 
	    reinterpret_cast<void *>(A[iSub]->data()));

}

template<int dim, int dimLS>
void MatVecProdH1MultiPhase<dim,dimLS>::evaluateHH(DistVec<double> &hhterm,
					      DistVec<double> &bcVal ) {

  *hhVal = bcVal;
}


//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MatVecProdH1MultiPhase<dim,dimLS>::evaluate(int it, DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                                            DistSVec<double,dim> &Q, DistSVec<double,dimLS> &Phi,
                                            DistSVec<double,dim> &F)
{
  if (!this->isFSI)
    this->spaceOp->computeJacobian(X, ctrlVol, Q, *this, *this->fluidSelector, this->riemann,this->timeState);
  else
    this->spaceOp->computeJacobian(this->riemann, X, Q,ctrlVol,this->fsi.LSS,
                                   this->fsi.Nriemann, *this->fluidSelector,*this,this->timeState);

  if (this->timeState)
    this->timeState->addToJacobian(ctrlVol, *this, Q);

  this->spaceOp->applyBCsToJacobian(Q, *this);

  if (areHHTermsActive) {

#pragma omp parallel for
    for (int iSub = 0; iSub < this->numLocSub; ++iSub) {

      this->subDomain[iSub]->computeJacobianFiniteVolumeTermHH(this->spaceOp->getFluxFcn(),
							       (*this->spaceOp->getDistBcData())(iSub) ,
							       (*this->spaceOp->getGeoState())(iSub),
							       ctrlVol, Q, *A[iSub],this->spaceOp->getVarFcn());
							      
    }
    
    this->timeState->addToHHJacobian(ctrlVol, *this, *hhVal);
  }


}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MatVecProdH1MultiPhase<dim,dimLS>::apply(DistSVec<double,dim> &p, DistSVec<double,dim> &prod)
{

  int iSub;
  
#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub) {
    this->subDomain[iSub]->computeMatVecProdH1(p.getMasterFlag(iSub), *A[iSub],
					 p(iSub), prod(iSub));
    this->subDomain[iSub]->sndData(*this->vecPat, prod.subData(iSub));
  }

  this->vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    this->subDomain[iSub]->addRcvData(*this->vecPat, prod.subData(iSub));

}

template<int dim, int dimLS>
void MatVecProdH1MultiPhase<dim,dimLS>::apply(DistEmbeddedVec<double,dim> &p, DistEmbeddedVec<double,dim> &prod)
{

  int iSub;
  
#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub) {
    this->subDomain[iSub]->computeMatVecProdH1(p.real().getMasterFlag(iSub), *A[iSub],
					       p.real()(iSub), prod.real()(iSub)/*, 
                                               p.ghost()(iSub), prod.ghost()(iSub)*/ );

    if (p.hasHHBoundaryTerm()) {
      if (!prod.hasHHBoundaryTerm()) {
    
        prod.setHH(p.hh());
      }
      prod.hh() = 0.0;
      this->subDomain[iSub]->
	computeMatVecProdH1FarFieldHH(p.real().getMasterFlag(iSub),
				      *A[iSub],p.real()(iSub), prod.real()(iSub), 
				      p.hh()(iSub), prod.hh()(iSub));
    }

    this->subDomain[iSub]->sndData(*this->vecPat, prod.real().subData(iSub));
  }


  this->vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    this->subDomain[iSub]->addRcvData(*this->vecPat, prod.real().subData(iSub));
/*
#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub) {
    this->subDomain[iSub]->sndData(*this->vecPat, prod.ghost().subData(iSub));
  }
  this->vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    this->subDomain[iSub]->addRcvData(*this->vecPat, prod.ghost().subData(iSub));
*/
}

//------------------------------------------------------------------------------

//----------------------------------------------------------------------------//
//                MatVecProd for Level Set equations                          //
//----------------------------------------------------------------------------//

template<int dim, int dimLS>
MatVecProdLS<dim,dimLS>::MatVecProdLS(DistTimeState<dim> *ts, DistGeoState *gs,
                                      MultiPhaseSpaceOperator<dim,dimLS> *spo, 
                                      Domain *dom, LevelSet<dimLS> *ls) :
  timeState(ts), spaceOp(spo), levelSet(ls), geoState(gs), 
  Qeps(dom->getNodeDistInfo()), Feps(dom->getNodeDistInfo()), DistMat<double,dimLS>(dom)
{
  X = 0; ctrlVol = 0; Q = 0; U = 0; F = 0; FluidId = 0;
  com = dom->getCommunicator();


  // H1 stuff
  A = new MvpMat<double,dimLS>*[this->numLocSub];

  double size = 0.0;

#pragma omp parallel for reduction (+: size)
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {
    A[iSub] = this->subDomain[iSub]->template createMaskMatVecProd<double,dimLS>();
    size += double(A[iSub]->numNonZeroBlocks()*dimLS*dimLS*sizeof(double)) / (1024.*1024.);
  }

  this->com->globalSum(1, &size);

  this->com->printf(2, "Memory required for MultiPhaseMatVec with H1 (dim=%d): %3.2f MB\n", dim, size);

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
MatVecProdLS<dim,dimLS>::~MatVecProdLS()
{
  timeState = 0;
  spaceOp   = 0;
  X = 0; ctrlVol = 0; Q = 0; U = 0; F = 0; FluidId = 0;
  com = 0;
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MatVecProdLS<dim,dimLS>::apply(DistSVec<double,dimLS> &p, DistSVec<double,dimLS> &prod)
{  

  int iSub;
  
#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub) {
    this->subDomain[iSub]->computeMatVecProdH1(p.getMasterFlag(iSub), *A[iSub],
					 p(iSub), prod(iSub));
    this->subDomain[iSub]->sndData(*this->vecPat, prod.subData(iSub));
  }

  this->vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    this->subDomain[iSub]->addRcvData(*this->vecPat, prod.subData(iSub));
  
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MatVecProdLS<dim,dimLS>::evaluate(int it, DistSVec<double,3> &x, DistVec<double> &cv,
				       DistSVec<double,dimLS> &q, DistSVec<double,dim> &u,
				       DistSVec<double,dim> &v, DistSVec<double,dimLS> &f, 
				       DistVec<int> &fluidId,bool requireSpecialBDF,DistLevelSetStructure* distLSS,
				       int lsMethod)
{

  X       = &x;
  ctrlVol = &cv;
  Q       = &q;
  U       = &u;
  F       = &f;
  V = &v;
  FluidId = &fluidId;

  spaceOp->computeJacobianLS(x,v,*ctrlVol, *Q,*this, fluidId,distLSS,lsMethod);

  if (timeState)
    timeState->addToJacobianLS(cv, *this, u,requireSpecialBDF);

}

//------------------------------------------------------------------------------
// WARNING: this routine is very different from MatVecProdFD<dim,neq>::computeEpsilon(...)
template<int dim, int dimLS>
double MatVecProdLS<dim, dimLS>::computeEpsilon(DistSVec<double,dimLS> &U, DistSVec<double,dimLS> &p)
{

  int iSub, size = 0;
  double eps0 = 1.e-6;

  const DistInfo &distInfo = U.info();

  double *alleps = reinterpret_cast<double *>(alloca(sizeof(double) * distInfo.numGlobSub));

  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) alleps[iSub] = 0.0;

#pragma omp parallel for reduction(+: size)
  for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

    int locOffset = distInfo.subOffset[iSub];

    int locsize = 0;
    double loceps = 0.0;

    for (int i=0; i<distInfo.subLen[iSub]; ++i) {

      if (distInfo.masterFlag[locOffset+i]) {
	for (int j=0; j<dimLS; ++j) {
	  ++locsize;
	  loceps += eps0*fabs(U[locOffset+i][j]) + eps0;
	}
      }

    }

    size += locsize;
    alleps[distInfo.locSubToGlobSub[iSub]] = loceps;

  }

  this->com->globalSum(1, &size);
  this->com->globalSum(distInfo.numGlobSub, alleps);

  double eps = 0.0;
  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) eps += alleps[iSub];

  double norm = sqrt(p*p);

  if (norm > 1.e-14) eps /= double(size) * norm;
  else eps = eps0;
 
  return eps;

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
DistMat<double,dimLS> &MatVecProdLS<dim,dimLS>::operator= (const double x)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub)
    *A[iSub] = x;

  return *this;

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
MatVecProd_dRdX<dim,Scalar,neq>::MatVecProd_dRdX
(
  IoData &ioData, VarFcn *varFcn, DistTimeState<dim> *ts,
  SpaceOperator<dim> *spo, Domain *domain, DistGeoState *gs
) :
  MatVecProd<dim,neq>(),
  timeState(ts),
  iod(&ioData),
  spaceOp(0),
  fluxFcn(0),
  X(0),
  ctrlVol(0),
  Q(0),
  F(0)
{

  numLocSub = domain->getNumLocSub();
  subDomain = domain->getSubDomain();
  com = domain->getCommunicator();
  dRdXop = new dRdXoperators<dim>();
  dRdXop->setNumLocSub(numLocSub);
  dRdXop->initialize();

  double size = 0.0;

  // allocate for viscous flux jacobian term
  bool nsFlag = false;
  spaceOp = new SpaceOperator<dim>(*spo, false);

// Included (MB*)
  if ((ioData.eqs.type == EquationsData::NAVIER_STOKES) && (ioData.bc.wall.integration != BcsWallData::WALL_FUNCTION))
    nsFlag = true;

#pragma omp parallel for reduction (+: size)
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    dRdXop->dEdgeNormdX[iSub] = subDomain[iSub]->template create_EdgeBaseddRdXoperators<3,3>();
    dRdXop->dFluxdEdgeNorm[iSub] = subDomain[iSub]->template create_NodeToEdgeBaseddRdXoperators<3,dim>();
    dRdXop->dFluxdFaceNormal[iSub] = subDomain[iSub]->template create_NodeToFaceBaseddRdXoperators<3,dim>();
    dRdXop->dFluxdFaceNormalVel[iSub] = subDomain[iSub]->template create_NodeToFaceBaseddRdXoperators<1,dim>();
    dRdXop->dFluxdUb[iSub] = subDomain[iSub]->template create_NodeToFaceBaseddRdXoperators<dim,dim>();
    dRdXop->dFaceNormdX[iSub] = subDomain[iSub]->template create_FaceBaseddRdXoperators<3,3>();
    dRdXop->dFidGradP[iSub] = subDomain[iSub]->template create_ConstantToNodeBaseddRdXoperators<3,3>();
    dRdXop->dFidX[iSub] = subDomain[iSub]->template create_ConstantToNodeBaseddRdXoperators<3,3>();
    dRdXop->dFvdX[iSub] = subDomain[iSub]->template create_ConstantToNodeBaseddRdXoperators<3,3>();
    dRdXop->dFidV[iSub] = subDomain[iSub]->template create_ConstantToNodeBaseddRdXoperators<dim,3>();
    dRdXop->dFvdV[iSub] = subDomain[iSub]->template create_ConstantToNodeBaseddRdXoperators<dim,3>();
    dRdXop->dFidS[iSub] = subDomain[iSub]->template create_ConstantToConstantBaseddRdXoperators<3,3>();
    dRdXop->dMidGradP[iSub] = subDomain[iSub]->template create_ConstantToNodeBaseddRdXoperators<3,3>();
    dRdXop->dMidX[iSub] = subDomain[iSub]->template create_ConstantToNodeBaseddRdXoperators<3,3>();
    dRdXop->dMvdX[iSub] = subDomain[iSub]->template create_ConstantToNodeBaseddRdXoperators<3,3>();
    dRdXop->dMvdV[iSub] = subDomain[iSub]->template create_ConstantToNodeBaseddRdXoperators<dim,3>();
    dRdXop->dMidV[iSub] = subDomain[iSub]->template create_ConstantToNodeBaseddRdXoperators<dim,3>();
    dRdXop->dMidS[iSub] = subDomain[iSub]->template create_ConstantToConstantBaseddRdXoperators<3,3>();
    dRdXop->dCtrlVoldX[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<3,1>();
    dRdXop->dRdX[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<3,6>();
    dRdXop->dRdR[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<6,6>();
    dRdXop->dddxdX[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<3,dim>();
    dRdXop->dddydX[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<3,dim>();
    dRdXop->dddzdX[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<3,dim>();
    dRdXop->dddxdV[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<dim,dim>();
    dRdXop->dddydV[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<dim,dim>();
    dRdXop->dddzdV[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<dim,dim>();
    dRdXop->dddxdR[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<6,dim>();
    dRdXop->dddydR[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<6,dim>();
    dRdXop->dddzdR[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<6,dim>();
    dRdXop->dFluxdddx[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<dim,dim>();
    dRdXop->dFluxdddy[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<dim,dim>();
    dRdXop->dFluxdddz[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<dim,dim>();
    dRdXop->dFluxdX[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<3,dim>();
    dRdXop->dViscousFluxdX[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<3,dim>();
    dRdXop->dGradPdddx[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<dim,3>();
    dRdXop->dGradPdddy[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<dim,3>();
    dRdXop->dGradPdddz[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<dim,3>();
    dRdXop->dForcedGradP[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<3,3>();
    dRdXop->dForcedX[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<3,3>();
    dRdXop->dForcedV[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<dim,3>();
    dRdXop->dForcedS[iSub] = subDomain[iSub]->template create_NodeToConstantBaseddRdXoperators<3,3>();
    dRdXop->dVdU[iSub] = subDomain[iSub]->template create_NodeBaseddRdXoperators<dim,dim>();
    dRdXop->dVdPstiff[iSub] = subDomain[iSub]->template create_NodeToConstantBaseddRdXoperators<1,dim>();
    size += double(dRdXop->dMidGradP[iSub]->numNonZeroBlocks()*3*3*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dMidX[iSub]->numNonZeroBlocks()*3*3*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dMvdX[iSub]->numNonZeroBlocks()*3*3*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dMvdV[iSub]->numNonZeroBlocks()*dim*3*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dMidV[iSub]->numNonZeroBlocks()*dim*3*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dMidS[iSub]->numNonZeroBlocks()*3*3*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dFidGradP[iSub]->numNonZeroBlocks()*3*3*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dFidX[iSub]->numNonZeroBlocks()*3*3*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dFvdX[iSub]->numNonZeroBlocks()*3*3*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dFidV[iSub]->numNonZeroBlocks()*dim*3*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dFvdV[iSub]->numNonZeroBlocks()*dim*3*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dFidS[iSub]->numNonZeroBlocks()*3*3*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dCtrlVoldX[iSub]->numNonZeroBlocks()*3*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dEdgeNormdX[iSub]->numNonZeroBlocks()*3*3*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dFaceNormdX[iSub]->numNonZeroBlocks()*3*3*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dFluxdFaceNormal[iSub]->numNonZeroBlocks()*3*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dFluxdFaceNormalVel[iSub]->numNonZeroBlocks()*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dFluxdUb[iSub]->numNonZeroBlocks()*dim*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dRdX[iSub]->numNonZeroBlocks()*3*6*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dRdR[iSub]->numNonZeroBlocks()*6*6*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dddxdX[iSub]->numNonZeroBlocks()*3*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dddydX[iSub]->numNonZeroBlocks()*3*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dddzdX[iSub]->numNonZeroBlocks()*3*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dddxdV[iSub]->numNonZeroBlocks()*dim*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dddydV[iSub]->numNonZeroBlocks()*dim*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dddzdV[iSub]->numNonZeroBlocks()*dim*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dddxdR[iSub]->numNonZeroBlocks()*6*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dddydR[iSub]->numNonZeroBlocks()*6*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dddzdR[iSub]->numNonZeroBlocks()*6*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dFluxdddx[iSub]->numNonZeroBlocks()*dim*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dFluxdddy[iSub]->numNonZeroBlocks()*dim*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dFluxdddz[iSub]->numNonZeroBlocks()*dim*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dFluxdX[iSub]->numNonZeroBlocks()*3*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dViscousFluxdX[iSub]->numNonZeroBlocks()*3*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dGradPdddx[iSub]->numNonZeroBlocks()*3*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dGradPdddy[iSub]->numNonZeroBlocks()*3*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dGradPdddz[iSub]->numNonZeroBlocks()*3*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dFluxdEdgeNorm[iSub]->numNonZeroBlocks()*3*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dForcedGradP[iSub]->numNonZeroBlocks()*3*3*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dForcedX[iSub]->numNonZeroBlocks()*3*3*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dForcedV[iSub]->numNonZeroBlocks()*dim*3*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dForcedS[iSub]->numNonZeroBlocks()*3*3*sizeof(double)) / (1024.*1024.);//TODO BUGHUNY bug here
    size += double(dRdXop->dVdU[iSub]->numNonZeroBlocks()*dim*dim*sizeof(double)) / (1024.*1024.);
    size += double(dRdXop->dVdPstiff[iSub]->numNonZeroBlocks()*dim*sizeof(double)) / (1024.*1024.);
  }

  com->globalSum(1, &size);
  
  com->printf(2, "Memory required for matvec with dRdX: ");
  com->printf(2, "%3.2f MB\n", size);

  fluxFcn = new FluxFcn*[BC_MAX_CODE - BC_MIN_CODE + 1]; 
  fluxFcn -= BC_MIN_CODE;
  if(BC_MAX_CODE-BC_MIN_CODE+1 < 22)
    fprintf(stderr,"Be prepared to see a segmentation fault shortly...\n");
  fluxFcn[BC_SYMMETRY] = new FluxFcn(0,BC_SYMMETRY,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_MASSFLOW_OUTLET_MOVING] = new FluxFcn(0,BC_MASSFLOW_OUTLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_MASSFLOW_OUTLET_FIXED] = new FluxFcn(0,BC_MASSFLOW_OUTLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_MASSFLOW_INLET_MOVING] = new FluxFcn(0,BC_MASSFLOW_INLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_MASSFLOW_INLET_FIXED] = new FluxFcn(0,BC_MASSFLOW_INLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_DIRECTSTATE_OUTLET_MOVING] = new FluxFcn(0,BC_DIRECTSTATE_OUTLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_DIRECTSTATE_OUTLET_FIXED] = new FluxFcn(0,BC_DIRECTSTATE_OUTLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_DIRECTSTATE_INLET_MOVING] = new FluxFcn(0,BC_DIRECTSTATE_INLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_DIRECTSTATE_INLET_FIXED] = new FluxFcn(0,BC_DIRECTSTATE_INLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_OUTLET_MOVING] = new FluxFcn(0,BC_OUTLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_OUTLET_FIXED] = new FluxFcn(0,BC_OUTLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_INLET_MOVING] = new FluxFcn(0,BC_INLET_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_INLET_FIXED] = new FluxFcn(0,BC_INLET_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_POROUS_WALL_MOVING] = new FluxFcn(0,BC_POROUS_WALL_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_POROUS_WALL_FIXED] = new FluxFcn(0,BC_POROUS_WALL_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ADIABATIC_WALL_MOVING] = new FluxFcn(0,BC_ADIABATIC_WALL_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ADIABATIC_WALL_FIXED] = new FluxFcn(0,BC_ADIABATIC_WALL_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_SLIP_WALL_MOVING] = new FluxFcn(0,BC_SLIP_WALL_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_SLIP_WALL_FIXED] = new FluxFcn(0,BC_SLIP_WALL_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ISOTHERMAL_WALL_MOVING] = new FluxFcn(0,BC_ISOTHERMAL_WALL_MOVING,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_ISOTHERMAL_WALL_FIXED] = new FluxFcn(0,BC_ISOTHERMAL_WALL_FIXED,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  fluxFcn[BC_INTERNAL] = new FluxFcn(0,BC_INTERNAL,ioData,varFcn,FluxFcnBase::PRIMITIVE);
  spaceOp->setFluxFcn(fluxFcn); //TODO: should avoid doing this!

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
MatVecProd_dRdX<dim,Scalar,neq>::~MatVecProd_dRdX()
{ 

  if (spaceOp) delete spaceOp;
  fluxFcn = 0; // deleted by spaceOperator.
  if (dRdXop) delete dRdXop;

}

//------------------------------------------------------------------------------
template<int dim, class Scalar, int neq>
void MatVecProd_dRdX<dim,Scalar,neq>::initializeOperators(double x)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
	*dRdXop->dMidX[iSub] = x;
	*dRdXop->dMvdX[iSub] = x;
	*dRdXop->dMvdV[iSub] = x;
	*dRdXop->dMidV[iSub] = x;
	*dRdXop->dMidS[iSub] = x;
	*dRdXop->dMidGradP[iSub] = x;
	*dRdXop->dFidX[iSub] = x;
	*dRdXop->dFvdX[iSub] = x;
	*dRdXop->dFidV[iSub] = x;
	*dRdXop->dFvdV[iSub] = x;
	*dRdXop->dFidS[iSub] = x;
	*dRdXop->dFidGradP[iSub] = x;
    *dRdXop->dCtrlVoldX[iSub] = x;
    *dRdXop->dEdgeNormdX[iSub] = x;
    *dRdXop->dFaceNormdX[iSub] = x;
    *dRdXop->dRdX[iSub] = x;
    *dRdXop->dRdR[iSub] = x;
    *dRdXop->dddxdX[iSub] = x;
    *dRdXop->dddydX[iSub] = x;
    *dRdXop->dddzdX[iSub] = x;
    *dRdXop->dddxdV[iSub] = x;
    *dRdXop->dddydV[iSub] = x;
    *dRdXop->dddzdV[iSub] = x;
    *dRdXop->dddxdR[iSub] = x;
    *dRdXop->dddydR[iSub] = x;
    *dRdXop->dddzdR[iSub] = x;
    *dRdXop->dFluxdddx[iSub] = x;
    *dRdXop->dFluxdddy[iSub] = x;
    *dRdXop->dFluxdddz[iSub] = x;
    *dRdXop->dFluxdX[iSub] = x;
    *dRdXop->dViscousFluxdX[iSub] = x;
    *dRdXop->dGradPdddx[iSub] = x;
    *dRdXop->dGradPdddy[iSub] = x;
    *dRdXop->dGradPdddz[iSub] = x;
    *dRdXop->dForcedX[iSub] = x;
    *dRdXop->dForcedS[iSub] = x;
    *dRdXop->dForcedV[iSub] = x;
    *dRdXop->dForcedGradP[iSub] = x;
    *dRdXop->dFluxdEdgeNorm[iSub] = x;
    *dRdXop->dFluxdFaceNormal[iSub] = x;
    *dRdXop->dFluxdFaceNormalVel[iSub] = x;
    *dRdXop->dFluxdUb[iSub] = x;
    *dRdXop->dVdU[iSub] = x;
    *dRdXop->dVdPstiff[iSub] = x;
  }

}

//------------------------------------------------------------------------------
template<int dim, class Scalar, int neq>
void MatVecProd_dRdX<dim,Scalar,neq>::constructOperators(Vec3D &x0,
                                                         DistSVec<double,3> &X,
                                                         DistVec<double> &ctrlVol,
                                                         DistSVec<double,dim> &U,
                                                         double dMach,
                                                         DistSVec<double,dim> &R,
                                                         DistVec<double> &Pin,
                                                         DistTimeState<dim> *timeState,
                                                         PostOperator<dim> *postOp)
{
  com->printf(5," ... in MatVecProd_dRdX<dim,Scalar,neq>::constructOperators, norm of X is %e\n", X.norm());
  initializeOperators(0.0);
  spaceOp->computeDerivativeOperators(x0, X, ctrlVol, U, dMach, R, Pin, timeState, postOp, dRdXop);
}


