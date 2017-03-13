#include <Communicator.h>

#include <cmath>
#include <iostream>

//------------------------------------------------------------------------------

template<int dim>
ImplicitRomPostproTsDesc<dim>::ImplicitRomPostproTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom), Uinitial(dom->getNodeDistInfo())
{
	this->maxItsNewton = 0;	// never do iterations

  reducedCoordsFile = NULL;

  if (!(strcmp(this->ioData->input.reducedCoords,"")==0)) {
    char *fullReducedCoordsName = new char[strlen(this->ioData->input.prefix) + 1 + strlen(this->ioData->input.reducedCoords) + 1];
    sprintf(fullReducedCoordsName, "%s%s", this->ioData->input.prefix, this->ioData->input.reducedCoords);
    reducedCoordsFile = fopen(fullReducedCoordsName, "r");
    if (!reducedCoordsFile)  {
      this->com->fprintf(stderr, "*** Error opening reduced coordinates file: %s\n", fullReducedCoordsName);
      fflush(stderr);
      exit (-1);
    }
    delete [] fullReducedCoordsName;
  }

  // avoid calling computeTimeStep for steady
  dt = 0.0;

  if (this->ioData->output.rom.resjacfrequency >= 1)
    this->rom->initializeClusteredOutputs();

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomPostproTsDesc<dim>::computeFullResidual(int it, DistSVec<double, dim> &Q, bool applyWeighting,  DistSVec<double, dim> *R, bool includeHomotopy) {


}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomPostproTsDesc<dim>::computeAJ(int it, DistSVec<double, dim> &Q, bool applyWeighting,  DistSVec<double, dim> *R) {

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomPostproTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop, DistSVec<double, dim> &U, const int& totalTimeSteps)  {

  breakloop = true;	// after loop, exit will occur because maxItsNewton = 1
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomPostproTsDesc<dim>::checkLocalRomStatus(DistSVec<double, dim> &U, const int totalTimeSteps)  {

  // checks whether the local ROM needs to be modified

  int tmp, _n, closestCluster, nCoords, prevCluster; 
  char switchStr[50], updateStr[50];

  if (this->currentCluster == -1) { // first iteration
    if (this->ioData->romOnline.distanceComparisons)
      this->rom->initializeDistanceComparisons(U);
    if (this->ioData->romOnline.basisUpdates==NonlinearRomOnlineData::UPDATES_FAST_EXACT) {
      this->rom->initializeFastExactUpdatesQuantities(U);
    } else if (this->ioData->romOnline.projectSwitchStateOntoAffineSubspace!=NonlinearRomOnlineData::PROJECT_OFF) {
      this->rom->initializeProjectionQuantities(U);
    }
  }

  _n = fscanf(reducedCoordsFile, "%d %s %s %d %d %d %d %d", &tmp, switchStr, updateStr, &closestCluster, &nCoords, &tmp, &tmp, &tmp);
  if (_n == EOF) {
    this->com->fprintf(stdout, "*** Error: reached end of reduced coordinates file before final time step.  Exiting...\n\n\n");
    exit(-1);  
  }

  this->com->fprintf(stdout, "%s %s %d %d\n", switchStr, updateStr, closestCluster, nCoords);

  this->clusterSwitch = (strcmp(switchStr,"switch")==0) ? true : false;
  this->updatePerformed = (strcmp(updateStr,"update")==0) ? true : false;

  if (this->updatePerformed || this->clusterSwitch) this->loadCluster(closestCluster, this->clusterSwitch, U);

  if (this->nPod != nCoords) {
    this->com->fprintf(stderr, "*** Error: dimension of reduced coordinates (%d) does not match dimension of state ROB (%d)\n"
                       ,nCoords, this->nPod);
    exit(-1);
  }

  double tmp2;
  for (int iPod = 0; iPod < this->nPod; ++iPod) {
    _n = fscanf(reducedCoordsFile, "%le", &tmp2);
    this->dUromTimeIt[iPod] = tmp2;
    this->dUromNewtonIt[iPod] = tmp2;
  }
  this->dUromCurrentROB += this->dUromTimeIt;
  if (this->UromCurrentROB.size() == this->dUromTimeIt.size()) 
    this->UromCurrentROB += this->dUromTimeIt;
  if (this->ioData->romOnline.distanceComparisons)
    this->rom->advanceDistanceComparisons(this->currentCluster, this->dUromTimeIt, this->UromCurrentROB);

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomPostproTsDesc<dim>::postProStep(DistSVec<double, dim> &U, int totalTimeSteps)  {


  // option to output residual and action-of-jacobian during post processing (before incementing)
  int freq = this->ioData->output.rom.resjacfrequency;
  if ((freq >= 1) && (totalTimeSteps%freq == 0) ) {
    ImplicitRomTsDesc<dim>::computeFullResidual(0, U, false);

    if (this->rom->jacActionSnapsFileNameSpecified) { 
      double tAJ = this->timer->getTime();
      ImplicitRomTsDesc<dim>::computeAJ(0, U, true);     // skipped some times for Broyden
      this->timer->addAJTime(tAJ);
    }

    this->saveNewtonSystemVectorsAction(totalTimeSteps);
  }

  // increment solution
  DistSVec<double, dim> dU(this->domain->getNodeDistInfo());
	this->expandVector(this->dUromTimeIt, dU); // solution increment in full coordinates
	U += dU;

  // option to output residual and action-of-jacobian during post processing (after incrementing)
  if ( (this->ioData->ts.implicit.type != ImplicitData::SPATIAL_ONLY)
       && !(this->rom->jacActionSnapsFileNameSpecified) // can't output jacAction without the appropriate dUromNewtonIt 
       && ((freq >= 1) && (totalTimeSteps%freq == 0))) {

    ImplicitRomTsDesc<dim>::computeFullResidual(1, U, false);

    this->saveNewtonSystemVectorsAction(totalTimeSteps);
  }

}

//------------------------------------------------------------------------------

template<int dim>
bool ImplicitRomPostproTsDesc<dim>::monitorConvergence(int it, DistSVec<double,dim> &U)
{// avoid residual calculation

  int freq = this->ioData->output.rom.resjacfrequency;
  if ((freq >= 1) && (it%freq == 0) ) {
    this->data->residual = this->F.norm();
  }

  return false;

}

//------------------------------------------------------------------------------

template<int dim>
double ImplicitRomPostproTsDesc<dim>::computeTimeStep(int it, double *dtLeft, DistSVec<double,dim> &U, double angle)
{
  double t0 = this->timer->getTime();
  this->data->computeCflNumber(it - 1, this->data->residual / this->restart->residual, angle);
  int numSubCycles = 1;

  if(it==1) {
    if(this->failSafeFlag == false){
      dt = this->timeState->computeTimeStep(this->data->cfl, this->data->dualtimecfl, dtLeft, &numSubCycles, *this->geoState, *this->X, *this->A, U);
    } else {
      dt = this->timeState->computeTimeStepFailSafe(dtLeft, &numSubCycles);
    }
  }

  if (this->problemType[ProblemData::UNSTEADY])
    this->com->printf(5, "Global dt: %g (remaining subcycles = %d)\n", dt*this->refVal->time, numSubCycles);
  this->timer->addFluidSolutionTime(t0);
  this->timer->addTimeStepTime(t0);

  return dt;
}


