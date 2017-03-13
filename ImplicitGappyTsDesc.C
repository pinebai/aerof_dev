#include <Communicator.h>

#include <cmath>
#include <VecSetOp.h>

//------------------------------------------------------------------------------

template<int dim>
ImplicitGappyTsDesc<dim>::ImplicitGappyTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom)
{

  jactmp = NULL;
  column = NULL;
  
}

//------------------------------------------------------------------------------

template<int dim>
ImplicitGappyTsDesc<dim>::~ImplicitGappyTsDesc() 
{

  if (jactmp) delete [] jactmp;
  if (column) delete [] column;
  
  ResRestrict.reset();
  AJRestrict.reset();

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGappyTsDesc<dim>::deleteRestrictedQuantities() {

  ResRestrict.reset();
  AJRestrict.reset();
  
  this->rom->deleteRestrictedQuantities();

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGappyTsDesc<dim>::computeFullResidual(int it, DistSVec<double, dim> &Q, bool applyWeighting,  DistSVec<double, dim> *R, bool includeHomotopy)
{
  // Evaluate residual on full mesh

  if (R==NULL) R=&(this->F);

  //clipping
  DistSVec<double, dim> Qeval(Q); 
  if (dim>5) {
    this->checkSolution(Qeval);
  } 

  double tRes = this->timer->getTime(); 

  this->spaceOp->computeResidualRestrict(*this->X, *this->A, Qeval, *R, this->timeState, (this->rom->restrictMapping())->getRestrictedToOriginLocNode());

  if (includeHomotopy) {
    this->timeState->add_dAW_dtRestrict(it, *this->geoState, *this->A, Qeval, *R, (this->rom->restrictMapping())->getRestrictedToOriginLocNode());
  }

  this->spaceOp->applyBCsToResidual(Qeval, *R);

  if (applyWeighting && (this->ioData->romOnline.residualScaling != NonlinearRomOnlineData::SCALING_OFF)) {
    if (this->componentwiseScalingVec) delete this->componentwiseScalingVec;
    this->componentwiseScalingVec = this->varFcn->computeScalingVec(*(this->ioData),Q,R);
    *R *= *this->componentwiseScalingVec;
  }

  this->timer->addResidualTime(tRes); 

  double t0 = this->timer->getTime();
  (this->rom->restrictMapping())->restriction(*R, *ResRestrict);
  this->timer->addRestrictionTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGappyTsDesc<dim>::computeAJ(int it, DistSVec<double, dim> &Q, bool applyWeighting, DistSVec<double, dim> *R)  {

  // Evaluate action of Jacobian on full mesh
  if (R==NULL) R = &this->F;

  bool componentScaling = (applyWeighting && (this->ioData->romOnline.residualScaling != NonlinearRomOnlineData::SCALING_OFF));

  double t0 = this->timer->getTime();
  if (componentScaling) {
    this->mvp->evaluateWeightedRestrict(it, *this->X, *this->A, Q, *R, 
                                        *(this->rom->restrictMapping()), this->varFcn,this, &TsDesc<dim>::checkSolution);
  } else { 
    this->mvp->evaluateRestrict(it, *this->X, *this->A, Q, *R, *(this->rom->restrictMapping()), this, &TsDesc<dim>::checkSolution);
  }
  this->timer->addJacEvaluateTime(t0);

  if (this->AJ.numVectors()!=this->nPod) this->AJ.resize(this->nPod);

  t0 = this->timer->getTime();
  for (int iPod = 0; iPod < this->nPod; iPod++) {
    if (componentScaling) {
      this->mvp->applyWeightedRestrict(this->pod[iPod], this->AJ[iPod], *(this->rom->restrictMapping()),this->varFcn, 
                             this, &TsDesc<dim>::checkSolution);
    } else {
      this->mvp->applyRestrict(this->pod[iPod], this->AJ[iPod], *(this->rom->restrictMapping()),
                             this, &TsDesc<dim>::checkSolution);
    }
  }
  this->timer->addJacApplyTime(t0);

  t0 = this->timer->getTime();
  for (int iPod = 0; iPod < this->nPod; iPod++) { // TODO only on local pod
    (this->rom->restrictMapping())->restriction(this->AJ[iPod], (*AJRestrict)[iPod]);
  }
  this->timer->addRestrictionTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
bool ImplicitGappyTsDesc<dim>::breakloop1(const bool breakloop) {

	return false;
	
}

//------------------------------------------------------------------------------

template<int dim>
bool ImplicitGappyTsDesc<dim>::breakloop2(const bool breakloop) {

	return breakloop;

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGappyTsDesc<dim>::monitorInitialState(int it, DistSVec<double,dim> &U)
{

  this->com->printf(2, "State vector norm = %.12e\n", sqrt(U*U));  

  if (!this->unsteady) {
    this->com->printf(2, "\nNOTE: For steady gappy simulations the reported residual is calculated using only the sample mesh,\n");
    this->com->printf(2, "      and is relative to the residual of the initial condition calculated on the same sample mesh.\n");
    this->com->printf(2, "      (This reference residual is re-restricted after every cluster switch for consistency).\n");
 
    this->Uinit = new DistSVec<double, dim>(this->domain->getNodeDistInfo());
    *(this->Uinit) = U;  // needed for computing the restricted residual after each cluster switch
  }

  this->com->printf(2, "\n");

}

//------------------------------------------------------------------------------

template<int dim>
bool ImplicitGappyTsDesc<dim>::checkForLastIteration(IoData &ioData, int it, double t, double dt, DistSVec<double,dim> &U)
{

  if (!this->unsteady && monitorConvergence(it, U))
    return true;

  if (!this->problemType[ProblemData::AERO] && !this->problemType[ProblemData::THERMO] && it >= this->data->maxIts) return true;

  if (this->unsteady)
    if(t >= this->data->maxTime - 0.01 * dt)
      return true;

  return false;

}

//------------------------------------------------------------------------------

template<int dim>
bool ImplicitGappyTsDesc<dim>::monitorConvergence(int it, DistSVec<double,dim> &U)
{// only called for steady simulations

  this->data->residual = computeGappyResidualNorm(U);

  if (this->data->residual == 0.0 || this->data->residual < this->data->eps * this->restart->residual || this->data->residual < this->data->epsabs)
    return true;
  else
    return false;

}

//------------------------------------------------------------------------------

template<int dim>
double ImplicitGappyTsDesc<dim>::computeGappyResidualNorm(DistSVec<double,dim>& Q)
{ // spatial only

  this->spaceOp->computeResidualRestrict(*this->X, *this->A, Q, this->F, this->timeState, (this->rom->restrictMapping())->getRestrictedToOriginLocNode());

  this->spaceOp->applyBCsToResidual(Q, this->F);

  double t0 = this->timer->getTime();

  (this->rom->restrictMapping())->restriction(this->F, *ResRestrict); 

  this->timer->addRestrictionTime(t0);

  double res = 0.0;
  res = (*ResRestrict) * (*ResRestrict);

  return sqrt(res);

}


//------------------------------------------------------------------------------

template<int dim>
void ImplicitGappyTsDesc<dim>::setReferenceResidual()
{
  if (this->Uinit) this->restart->residual = computeGappyResidualNorm(*(this->Uinit));

  this->com->printf(2, "Norm of restricted reference residual = %.12e\n", this->restart->residual);

}




