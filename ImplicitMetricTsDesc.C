#include <Communicator.h>

#include <cmath>
#include <VecSetOp.h>

//------------------------------------------------------------------------------

template<int dim>
ImplicitMetricTsDesc<dim>::ImplicitMetricTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitGappyTsDesc<dim>(ioData, geoSource, dom)
	//, leastSquaresSolver(this->com, this->com->size(), 1)// all cpus along rows
{

  this->projVectorTmp = NULL;
  
}

//------------------------------------------------------------------------------

template<int dim>
ImplicitMetricTsDesc<dim>::~ImplicitMetricTsDesc() {}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitMetricTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop, DistSVec<double, dim> &U, const int& totalTimeSteps)  {

  // Form the normal equations
  double t0 = this->timer->getTime();

  DistSVec<double, dim> weightedResRestrict(*this->ResRestrict);
  if (this->rom->metric->numVectors() == 1) {
    weightedResRestrict *= (*this->rom->metric)[0];
  } else {
    this->com->fprintf(stderr, "... the metric has %d vectors...\n", this->rom->metric->numVectors());
    this->com->fprintf(stderr, "... this hasn't been coded yet ...\n");
    sleep(1);
    exit(-1);
  }
  this->projectVector(*this->AJRestrict, weightedResRestrict, From);     // different from PG
  rhs = -1.0 * From;
  
  res = rhs*rhs;
  
  if (res < 0.0){
    fprintf(stderr, "*** negative residual: %e\n", res);
    exit(1);
  }
  res = sqrt(res);

  // Check onvergence criterion
  if (it == 0) {
    this->res0 = res;
    this->target = this->epsNewton * this->res0;
    this->com->fprintf(stderr, " ... ROM residual = %e (target = %e)\n", res, this->target);
  } else {
    this->com->fprintf(stderr, " ... ROM residual = %e\n", res);
  }
  if (breakloop = (res == 0.0) || (res <= this->target)) return;

 
  VecSet< DistSVec<double, dim> > weightedAJRestrict(*this->AJRestrict);
  if (this->rom->metric->numVectors() == 1) {
    for (int iVec=0; iVec<weightedAJRestrict.numVectors(); ++iVec)
      weightedAJRestrict[iVec] *= (*this->rom->metric)[0];
  } else {
    this->com->fprintf(stderr, "... the metric has %d vectors...\n", this->rom->metric->numVectors());
    this->com->fprintf(stderr, "... this hasn't been coded yet ...\n");
    sleep(1);
    exit(-1);
  }

  transMatMatProd(*this->AJRestrict,weightedAJRestrict,this->jactmp);
  for (int iRow = 0; iRow < this->nPod; ++iRow) {
    for (int iCol = 0; iCol < this->nPod; ++iCol) { // different from PG
      this->jac[iRow][iCol] = this->jactmp[iRow + iCol * this->pod.numVectors()];
    }
  }

  // homotopy on reduced-coordinates for spatial-only problems
  if (this->spatialOnlyWithHomotopy) {
    double homotopyStep = min(this->homotopyStepInitial*pow(this->homotopyStepGrowthRate,totalTimeSteps), this->homotopyStepMax);
    double invHomotopyStep = 1/homotopyStep;
    Vec<double> dUrom(this->dUromTimeIt);
    dUrom *= invHomotopyStep;
    rhs -= dUrom;
    for (int iDiag = 0; iDiag < this->nPod; ++iDiag) {
      this->jac[iDiag][iDiag] += invHomotopyStep;
    }
  }

  this->timer->addLinearSystemFormTime(t0);

  // Solve the normal equations
  t0 = this->timer->getTime();
  this->solveLinearSystem(it, rhs, this->dUromNewtonIt);
  this->timer->addLinearSystemSolveTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitMetricTsDesc<dim>::setProblemSize(DistSVec<double, dim> &U) {

  this->jac.setNewSize(this->nPod,this->nPod);
  From.resize(this->nPod);
  rhs.resize(this->nPod);

  if (this->projVectorTmp) delete [] (this->projVectorTmp);
  this->projVectorTmp = new double [this->nPod];

  //nSampleNodes = this->rom->sampleNodes.size();
  //leastSquaresSolver.problemSizeIs(nSampleNodes*dim, this->nPod);

  if (this->jactmp) delete [] this->jactmp;
  if (this->column) delete [] this->column;
  this->jactmp = new double [this->nPod * this->nPod];
  this->column = new double [this->nPod];

  this->AJRestrict.reset(new VecSet<DistSVec<double, dim> >(this->nPod, this->rom->getRestrictedDistInfo()));
  this->ResRestrict.reset(new DistSVec<double, dim> (this->rom->getRestrictedDistInfo())); 

}

//------------------------------------------------------------------------------

template<int dim>
double ImplicitMetricTsDesc<dim>::meritFunction(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> &dQ, DistSVec<double, dim> &F, double stepLength)  {
  // merit function: norm of the residual (want to minimize residual)

  DistSVec<double, dim> newQ(this->domain->getNodeDistInfo());
  newQ = Q + stepLength*dQ;
  this->checkSolution(newQ);
  this->computeFullResidual(it,newQ,true,&F);

  DistSVec<double, dim> weightedResRestrict(*this->ResRestrict);
  if (this->rom->metric->numVectors() == 1) {
    weightedResRestrict *= (*this->rom->metric)[0];
  } else {
    this->com->fprintf(stderr, "... the metric has %d vectors...\n", this->rom->metric->numVectors());
    this->com->fprintf(stderr, "... this hasn't been coded yet ...\n");
    sleep(1);
    exit(-1);
  }

  double merit = 0.0;

  if (this->ioData->romOnline.meritFunction == NonlinearRomOnlineData::HDM_RESIDUAL) {
    // merit function = 1/2 * (norm of full-order residual)^2
    merit = *this->ResRestrict * weightedResRestrict;	
    merit *= 0.5;
  } else if (this->ioData->romOnline.meritFunction == NonlinearRomOnlineData::ROM_RESIDUAL) {
    // Form the ROM residual
    Vec<double> romResidual;
    romResidual.resize(this->nPod);
    this->computeAJ(it, newQ, true);
    this->projectVector(*this->AJRestrict, weightedResRestrict, romResidual);
    merit = romResidual*romResidual;
    merit *= 0.5;
  } else {
    fprintf(stderr,"*** Error: unrecognized choice of merit function!\n");
    exit(-1);
  }

  return merit;

}

//------------------------------------------------------------------------------


