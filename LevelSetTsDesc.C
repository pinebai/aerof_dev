#include <LevelSetTsDesc.h>

#include "IoData.h"
#include "GeoSource.h"
#include "Domain.h"
#include "LevelSet.h"
#include <DistExactRiemannSolver.h>
#include <FluidSelector.h>

#include <cmath>
#include <limits>
                                                                                                        
#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
#endif


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
LevelSetTsDesc<dim,dimLS>::
LevelSetTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  TsDesc<dim>(ioData, geoSource, dom), Phi(this->getVecInfo()), V0(this->getVecInfo()),
  PhiV(this->getVecInfo()),
  fluidSelector(ioData.eqs.numPhase, ioData, dom),umax(this->getVecInfo()), programmedBurn(NULL),Utilde(this->getVecInfo()),
  VWeights(this->getVecInfo()), Weights(this->getVecInfo()),ioData(ioData)

{
  multiPhaseSpaceOp = new MultiPhaseSpaceOperator<dim,dimLS>(ioData, this->varFcn, this->bcData, this->geoState, 
                                                             this->domain, this->V);
  this->timeState = new DistTimeState<dim>(ioData, multiPhaseSpaceOp, this->varFcn, this->domain, this->V);

  if (ioData.mf.levelSetMethod == MultiFluidData::TRIANGULATED) 
    lsMethod = 2;
  else if (ioData.mf.levelSetMethod == MultiFluidData::PRIMITIVE)
    lsMethod = 1;
  else
    lsMethod = 0;

  LS = new LevelSet<dimLS>(ioData, this->domain);
  riemann = new DistExactRiemannSolver<dim>(ioData,this->domain,this->varFcn);

  frequencyLS = ioData.mf.frequency;
  interfaceType = ioData.mf.interfaceType;

  Prate = ioData.implosion.Prate;
  Pinit = ioData.implosion.Pinit;
  tmax = (Prate == 0) ? std::numeric_limits<double>::max() : (ioData.bc.inlet.pressure - Pinit)/Prate;

  requireSpecialBDF = false;

  int numBurnableFluids = ProgrammedBurn::countBurnableFluids(ioData);
  //std::cout << "Num burnable fluids = " << numBurnableFluids << std::endl;
  if (numBurnableFluids > 0) {
    programmedBurn = new ProgrammedBurn(ioData,this->X);
    this->fluidSelector.attachProgrammedBurn(programmedBurn);
  }

  phaseChangeType = 0;
  
  if (ioData.mf.typePhaseChange == MultiFluidData::EXTRAPOLATION) {

    phaseChangeType = 1;
  }

  interfaceOrder = 1;

  if (ioData.mf.interfaceTreatment == MultiFluidData::SECONDORDER) {

    dom->createHigherOrderMultiFluid();

    interfaceOrder = 2;

#pragma omp parallel for
    for (int iSub = 0; iSub < dom->getNumLocSub(); ++iSub) {
      V6NodeData (*v6data)[2];
      v6data = 0;
      dom->getSubDomain()[iSub]->findEdgeTetrahedra((*this->X)(iSub), v6data);
      dom->getSubDomain()[iSub]->getHigherOrderMF()->
	initialize<dim>(dom->getNodeDistInfo().subSize(iSub),
			dom->getSubDomain()[iSub]->getElems(),
			v6data); 

      if (ioData.mf.interfaceLimiter == MultiFluidData::LIMITERALEX1)
        dom->getSubDomain()[iSub]->getHigherOrderMF()->setLimitedExtrapolation();
    }
  }

  if (strcmp(ioData.input.exactInterfaceLocation,"") != 0) {
    
    loadExactInterfaceFile(ioData, ioData.input.exactInterfaceLocation);
  }

  myTriangulatedInterface = NULL;
  if (ioData.mf.testCase == 1 && lsMethod == 2) {

    myTriangulatedInterface = new TriangulatedInterface();
    myTriangulatedInterface->initializeAsSquare(11);

  }

  limitHigherOrderExtrapolation = false;
  if (ioData.mf.interfaceLimiter == MultiFluidData::LIMITERALEX1)
    limitHigherOrderExtrapolation = true;
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
LevelSetTsDesc<dim,dimLS>::~LevelSetTsDesc()
{

  if (LS) delete LS;
  if (riemann) delete riemann;
}

//------------------------------------------------------------------------------
template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::setupTimeStepping(DistSVec<double,dim> *U, IoData &ioData)
{

  this->geoState->setup2(this->timeState->getData());

  // initalize solution
  this->timeState->setup(this->input->solutions, *this->X, this->bcData->getInletBoundaryVector(), *U, ioData);
  LS->setup(this->input->levelsets, *this->X, *U, Phi, ioData,&fluidSelector,this->varFcn,0,0,lsMethod);
  fluidSelector.initializeFluidIds(Phi, LS->Phinm1, LS->Phinm2); //setup fluidId in fluidSelector

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(this->mmh);
  if (_mmh) 
    _mmh->setup(&this->restart->frequency, &this->data->maxTime, this->postOp, *this->X, *U, fluidSelector.fluidId);
 
  *this->Xs = *this->X;

  this->initializeFarfieldCoeffs();

  this->timer->setSetupTime();
  
  if (this->modifiedGhidaglia)
    this->timeState->attachHH(*this->bcData->getBoundaryStateHH());

  if (lsMethod == 2) {

    myTriangulatedInterface->initializeIntersector(ioData,this->domain->getCommunicator(),
                                                   this->domain, *this->X);
    this->multiPhaseSpaceOp->attachTriangulatedInterface(myTriangulatedInterface);
  }
  
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
double LevelSetTsDesc<dim,dimLS>::computeTimeStep(int it, double *dtLeft,
                                            DistSVec<double,dim> &U, double angle)
{
  double t0 = this->timer->getTime();
  this->data->allowstop = this->timeState->allowcflstop;

  this->data->computeCflNumber(it - 1, this->data->residual / this->restart->residual, angle);

  umax = 0.0;
  int numSubCycles = 1;
  double dt,dtmax;
  if(TsDesc<dim>::timeStepCalculation == TsData::CFL || it==1){
    dt = this->timeState->computeTimeStep(this->data->cfl, this->data->dualtimecfl, dtLeft,
                              &numSubCycles, *this->geoState, *this->A, U, *(fluidSelector.fluidId),&umax);
  }
  else { //time step size with error estimation
    // First compute the maximum possible time step (based on the speed of the interface)
    // To do this, we simply pass in a cfl=10000 to the previous routine
    dtmax = this->timeState->computeTimeStep(10000, this->data->dualtimecfl, dtLeft,
                                             &numSubCycles, *this->geoState, *this->A,
                                             U, *(fluidSelector.fluidId),&umax);
    dt = this->timeState->computeTimeStep(it, dtLeft, &numSubCycles);
    dt = min(dt,dtmax);
    this->timeState->computeCoefficients(dt);
  }

  if (this->problemType[ProblemData::UNSTEADY])
    this->com->printf(5, "Global dt: %g (remaining subcycles = %d)\n",
                      dt*this->refVal->time, numSubCycles);
  
  if(TsDesc<dim>::timeStepCalculation == TsData::ERRORESTIMATION && it == 1)
    this->timeState->setDtMin(dt * TsDesc<dim>::data->getCflMinOverCfl0());

  this->timer->addFluidSolutionTime(t0);
  this->timer->addTimeStepTime(t0);

  // for algNum 22 the maxTime can be set to less than the current time, to trigger termination
  if (dt + this->currentTime > this->data->maxTime && (this->mmh && this->mmh->getAlgNum() != 22))
    dt = this->data->maxTime - this->currentTime; 

  currentTimeStep = dt;

  return dt;

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::updateStateVectors(DistSVec<double,dim> &U, int it)
{

  this->geoState->update(*this->X, *this->A);

  if(frequencyLS > 0 && it%frequencyLS == 0){
//    this->com->printf(5, "LevelSet norm before reinitialization = %e\n", Phi.norm());
    if (this->lsMethod == 0)
      LS->conservativeToPrimitive(Phi,PhiV,U);
    else
      PhiV = Phi;
    LS->reinitializeLevelSet(*this->X, PhiV);
    if (this->lsMethod == 0)
      LS->primitiveToConservative(PhiV,Phi,U);
    else
      Phi = PhiV;
//    this->com->printf(5, "LevelSet norm after reinitialization = %e\n", Phi.norm());
    LS->update(Phi);

    // If we are doing 3BDF, reinitialization destroys unm1
    // Create a new version using F(U) = dU/dt
    if (this->timeState->useNm1()) {
      DistSVec<double,dimLS>& Phinm1 = LS->getPhinm1();
      this->multiPhaseSpaceOp->computeResidualLS(*this->X, *this->A, Phi, *fluidSelector.fluidId, U, Phinm1);
      Phinm1 = -1.0*Phinm1;
      requireSpecialBDF = true;
    }      
  } else {
    requireSpecialBDF = false;
    LS->update(Phi);
  }

  if (programmedBurn) {

    programmedBurn->setFluidIds(currentTime, *fluidSelector.fluidId,U);
  }

  fluidSelector.update();



  this->timeState->update(U,Utilde,  *(fluidSelector.fluidIdn), fluidSelector.fluidIdnm1, riemann);

  this->spaceOp->updateFixes();
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
int LevelSetTsDesc<dim,dimLS>::checkSolution(DistSVec<double,dim> &U)
{

  int ierr = this->domain->checkSolution(this->varFcn, *this->A, U, *fluidSelector.fluidId, *fluidSelector.fluidIdn);

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::setupOutputToDisk(IoData &ioData, bool *lastIt,
                          int it, double t, DistSVec<double,dim> &U)
{
  if (it == this->data->maxIts)
    *lastIt = true;
  else
    this->monitorInitialState(it, U); // Phi?

  this->output->setMeshMotionHandler(ioData, this->mmh);
  this->output->openAsciiFiles();
  this->output->cleanProbesFile();

  this->timer->setSetupTime();
  double fluxNorm = 0.5*(this->data->residual)*(this->data->residual);

  if (it == 0) {
    this->output->writeForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, fluidSelector.fluidId);
    this->output->writeLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, fluidSelector.fluidId);
    this->output->writeMatchPressureToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, *this->A, U, this->timeState, fluidSelector.fluidId);
    this->output->writeFluxNormToDisk(it, 0, 0, t, fluxNorm);
    this->output->writeHydroForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, fluidSelector.fluidId);
    this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, fluidSelector.fluidId);
    this->output->writeResidualsToDisk(it, 0.0, 1.0, this->data->cfl);
    this->output->writeMaterialVolumesToDisk(it, 0.0, *this->A, fluidSelector.fluidId);
    this->output->writeMaterialMassEnergyToDisk(it, 0.0, U, *this->A, fluidSelector.fluidId);
    this->output->writeCPUTimingToDisk(*lastIt, it, t, this->timer);
    this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, *fluidSelector.fluidId,&Phi);
    this->output->writeHeatFluxesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, fluidSelector.fluidId);
  }

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::outputToDisk(IoData &ioData, bool* lastIt, int it,
                                               int itSc, int itNl,
                                               double t, double dt, DistSVec<double,dim> &U)
{
  this->com->globalSum(1, &interruptCode);
  if (interruptCode)
    *lastIt = true;

  double cpu = this->timer->getRunTime();
  double res = this->data->residual / this->restart->residual;
  double fluxNorm = 0.5*(this->data->residual)*(this->data->residual);
                                                                                                      
  this->output->writeLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, fluidSelector.fluidId);
  this->output->writeMatchPressureToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, *this->A, U, this->timeState, fluidSelector.fluidId);
  this->output->writeFluxNormToDisk(it, itSc, itNl, t, fluxNorm);
  this->output->writeHydroForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, fluidSelector.fluidId);
  this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, fluidSelector.fluidId);
  this->output->writeResidualsToDisk(it, cpu, res, this->data->cfl);
  this->output->writeMaterialVolumesToDisk(it, t, *this->A, fluidSelector.fluidId);
  this->output->writeMaterialMassEnergyToDisk(it, t, U, *this->A, fluidSelector.fluidId);
  this->output->writeCPUTimingToDisk(*lastIt, it, t, this->timer);
  this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState,*fluidSelector.fluidId,&Phi);
  this->output->writeProbesToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState,*fluidSelector.fluidId,&Phi);
  this->restart->writeToDisk(this->com->cpuNum(), *lastIt, it, t, dt, *this->timeState, *this->geoState, LS,NULL, &fluidSelector);

  this->output->updatePrtout(t);
  this->restart->updatePrtout(t);
  if (*lastIt) {
    this->timer->setRunTime();
    if (this->com->getMaxVerbose() >= 2)
      this->timer->print(this->domain->getStrTimer());
    this->output->closeAsciiFiles();

    this->output->template writeLinePlotsToDisk<1>(true, it, t, *this->X,
						   *this->A, U, 
						   this->timeState, *fluidSelector.fluidId);

    if (strcmp(ioData.input.convergence_file,"") != 0) {

      this->varFcn->conservativeToPrimitive(U, *this->V, fluidSelector.fluidId);
      computeConvergenceInformation(ioData,ioData.input.convergence_file,*this->V);
    }
    else if (ioData.mf.testCase == 2 ||
	     ioData.mf.testCase == 3) {
      
      //this->varFcn->conservativeToPrimitive(U, *this->V, fluidSelector.fluidId);
      computeConvergenceInformation(ioData,ioData.input.convergence_file,U);//*this->V);
    }
  }

}
//------------------------------------------------------------------------------

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::outputForces(IoData &ioData, bool* lastIt, int it, int itSc, int itNl,
                               double t, double dt, DistSVec<double,dim> &U)  {

  double cpu = this->timer->getRunTime();
  this->output->writeForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, fluidSelector.fluidId);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::resetOutputToStructure(DistSVec<double,dim> &U)
{
  this->com->printf(5,"LevelSetTsDesc<dim,dimLS>::resetOutputToStructure\n");

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(this->mmh);
  if (_mmh) 
    _mmh->resetOutputToStructure(this->postOp, *this->X, U, fluidSelector.fluidId);

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::updateOutputToStructure(double dt, double dtLeft,
					  DistSVec<double,dim> &U)
{

  this->com->printf(5,"LevelSetTsDesc<dim,dimLS>::resetOutputToStructure\n");
  if (this->mmh) {
    double work[2];
    this->mmh->computeInterfaceWork(dt, this->postOp, this->geoState->getXn(), 
                                    this->timeState->getUn(), *this->X, U, 
                                    work, fluidSelector.fluidIdn, fluidSelector.fluidId);
    this->restart->energy[0] += work[0];
    this->restart->energy[1] += work[1];
  }

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(this->mmh);
  if (_mmh)
    _mmh->updateOutputToStructure(dt, dtLeft, this->postOp, *this->X, U, fluidSelector.fluidId);

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
bool LevelSetTsDesc<dim,dimLS>::IncreasePressure(int it, double dt, double t, DistSVec<double,dim> &U)
{
  if(Pinit<0.0 || Prate<0.0) return true; // no setup for increasing pressure

  if(t>tmax && t-dt>tmax) return true; // max pressure was reached, so now we solve
  else{ // max pressure not reached, so we do not solve and we increase pressure and let structure react
    this->com->fprintf(stdout, "about to increase pressure to value of %e\n", Pinit+t*Prate);
    this->domain->IncreasePressure(Pinit+t*Prate, this->varFcn, U);
    return false;
  }


}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::avoidNewPhaseCreation(DistSVec<double,dimLS> &localPhi)
{

  this->domain->avoidNewPhaseCreation(localPhi, LS->Phin);

}

//------------------------------------------------------------------------------

template<int dim,int dimLS>
void LevelSetTsDesc<dim,dimLS>::fixSolution(DistSVec<double,dim>& U,DistSVec<double,dim>& dU) {

  //this->domain->fixSolution2(this->varFcn,U,dU,fluidSelector.fluidId);
  if (this->fixSol == 1)
    this->domain->fixSolution(this->varFcn,U,dU,fluidSelector.fluidId);
}


template<int dim,int dimLS>
void LevelSetTsDesc<dim,dimLS>::setCurrentTime(double t,DistSVec<double,dim>& U) { 

  currentTime = t;

  if (programmedBurn)
    programmedBurn->setCurrentTime(t,multiPhaseSpaceOp->getVarFcn(), U,*(fluidSelector.fluidId),*(fluidSelector.fluidIdn));
}

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::computeConvergenceInformation(IoData &ioData, const char* file, DistSVec<double,dim>& U) {

  DistSVec<double,dim> Uexact(U.info());

  std::cout << "current time = " << currentTime*ioData.ref.rv.time << std::endl;

  if (ioData.mf.testCase == 2 || 
      ioData.mf.testCase == 3 ) {

    
    if (ioData.mf.testCase == 2) {
      
      DistSVec<double,1> dummy1(U.info());
      DistVec<int> dummy2(U.info());
      
      ExactSolution::Fill<&ExactSolution::CylindricalBubble,
	dim, 1>(Uexact,dummy2,
		dummy1, *this->X, ioData,currentTime+currentTimeStep,
		this->spaceOp->getVarFcn());
      
    } else if (ioData.mf.testCase == 3) {
      
      DistSVec<double,1> dummy1(U.info());
      DistVec<int> dummy2(U.info());
      
      ExactSolution::Fill<&ExactSolution::AcousticTwoFluid,
	dim, 1>(Uexact,dummy2,
		dummy1, *this->X, ioData,currentTime+currentTimeStep,
		this->spaceOp->getVarFcn());
      
    }
  
  } else {
    
    OneDimensional::read1DSolution(ioData,file, Uexact, 
				   (DistSVec<double,dimLS>*)0,
				   &fluidSelector,
				   multiPhaseSpaceOp->getVarFcn(),
				   *this->X,
				   *this->domain,
				   OneDimensional::ModeU,
				   false) ;
    
  }
  
  double error[dim];
  double refs[dim] = {ioData.ref.rv.density, ioData.ref.rv.velocity,
		      ioData.ref.rv.velocity, ioData.ref.rv.velocity,
		      ioData.ref.rv.pressure};

  //std::cout << "ref pressure = " <<  ioData.ref.rv.pressure << " ref length =  " <<  ioData.ref.rv.length << std::endl;
  double tot_error = 0.0;
  DistVec<int> nnv(U.info());
  nnv = 1;
  int nNodes = nnv.sum();

  this->domain->computeL1Error(U,Uexact,*this->A,error);
  for (int k = 0; k < dim; ++k) {
    tot_error += error[k];
    this->domain->getCommunicator()->fprintf(stdout,"L1 error [%d]: %e\n", k, error[k]*refs[k]);
  }
  this->domain->getCommunicator()->fprintf(stdout,"L1 error (total): %e\n", tot_error);

  tot_error = 0.0;
  this->domain->computeL2Error(U,Uexact,*this->A,error);
  for (int k = 0; k < dim; ++k) {
    tot_error += error[k];
    this->domain->getCommunicator()->fprintf(stdout,"L2 error [%d]: %e\n", k, error[k]*refs[k]);
  }
  this->domain->getCommunicator()->fprintf(stdout,"L2 error (total): %e\n", tot_error);


  tot_error = 0.0;
  this->domain->computeLInfError(U,Uexact,error);
  for (int k = 0; k < dim; ++k) {
    tot_error = max(error[k],tot_error);
    this->domain->getCommunicator()->fprintf(stdout,"Linf error [%d]: %e\n", k, error[k]*refs[k]);
  }
  this->domain->getCommunicator()->fprintf(stdout,"Linf error (total): %e\n", tot_error);

  
}

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::loadExactInterfaceFile(IoData& ioData, const char* file) {

  FILE* f = fopen(file,"r");
  if (!f) {

    this->domain->getCommunicator()->fprintf(stderr,"Cannot open exact interface location file: %s\n",file);
    exit(-1);
  }
  
  while (!feof(f)) {

    exactInterfacePoint l;
    int nm = fscanf(f,"%lf %lf", &l.time,&l.loc);
    l.time /= ioData.ref.rv.time;
    l.loc /= ioData.ref.rv.length;
    myExactInterface.push_back(l);
  }

  fclose(f);
}

template<int dim, int dimLS>
void LevelSetTsDesc<dim,dimLS>::setPhiExact() {

  if (this->myExactInterface.size() > 0) {
    
    double il = 0.0, vl;
    for (int i = 0; i < this->myExactInterface.size(); ++i) {
      
      if (this->myExactInterface[i].time > this->currentTime) {
	
	il = -(this->myExactInterface[i-1].time-this->currentTime) * this->myExactInterface[i].loc / (this->myExactInterface[i].time-this->myExactInterface[i-1].time) + 
	  (this->myExactInterface[i].time-this->currentTime) * this->myExactInterface[i-1].loc / (this->myExactInterface[i].time-this->myExactInterface[i-1].time);
        vl = (this->myExactInterface[i].loc-this->myExactInterface[i-1].loc)/(this->myExactInterface[i].time-this->myExactInterface[i-1].time) ;
	break;
      }
    }

    if (myTriangulatedInterface)
      myTriangulatedInterface->setExactSquare(il,  vl);   
 
#pragma omp parallel for
    for (int iSub = 0; iSub < this->Phi.numLocSub(); ++iSub) {
      
      SVec<double,dimLS>& phi(this->Phi(iSub));
      for (int i = 0 ; i < phi.size(); ++i) {
	
	phi[i][0] = (*this->X)[i][0]-il;
      }
    }
  }

}

