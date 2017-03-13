#include <DistExactRiemannSolver.h>
#include <Domain.h>
#include <FluidSelector.h>
#include <FSI/DynamicNodalTransfer.h>
#include <GeoSource.h>
#include <LevelSet.h>
#include <MultiPhysicsTsDesc.h>

#ifdef DO_EMBEDDED
#include <IntersectorFRG/IntersectorFRG.h>
#include <IntersectorPhysBAM/IntersectorPhysBAM.h>
#endif

#include <cmath>
#include <limits>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
#endif

//------------------------------------------------------------------------------

template<int dim, int dimLS>
MultiPhysicsTsDesc<dim,dimLS>::
MultiPhysicsTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  TsDesc<dim>(ioData, geoSource, dom), Phi(this->getVecInfo()), V0(this->getVecInfo()),
  PhiV(this->getVecInfo()), PhiWeights(this->getVecInfo()), InterfaceTag(this->getVecInfo()), 
  fluidSelector(ioData.eqs.numPhase, ioData, dom), //memory allocated for fluidIds
  Vtemp(this->getVecInfo()), numFluid(ioData.eqs.numPhase),Wtemp(this->getVecInfo()),umax(this->getVecInfo()), programmedBurn(NULL)
{
  simType = (ioData.problem.type[ProblemData::UNSTEADY]) ? 1 : 0;
  orderOfAccuracy = (ioData.schemes.ns.reconstruction == SchemeData::CONSTANT) ? 1 : 2;

  if (ioData.mf.levelSetMethod == MultiFluidData::PRIMITIVE)
    lsMethod = 1;
  else
    lsMethod = 0;

  phaseChangeChoice  = (ioData.embed.eosChange==EmbeddedFramework::RIEMANN_SOLUTION) ? 1 : 0;
  phaseChangeAlg	 = (ioData.embed.phaseChangeAlg==EmbeddedFramework::LEAST_SQUARES) ? 1 : 0;
  interfaceAlg		 = (ioData.embed.interfaceAlg==EmbeddedFramework::INTERSECTION) ? 1 : 0;
  intersectAlpha	 = ioData.embed.alpha;

  multiPhaseSpaceOp = new MultiPhaseSpaceOperator<dim,dimLS>(ioData, this->varFcn, this->bcData, this->geoState, 
                                                             this->domain, this->V);
  this->timeState = new DistTimeState<dim>(ioData, multiPhaseSpaceOp, this->varFcn, this->domain, this->V);
  riemann = new DistExactRiemannSolver<dim>(ioData,this->domain,this->varFcn);

  int numBurnableFluids = ProgrammedBurn::countBurnableFluids(ioData);
  this->com->fprintf(stderr,"numBurnableFluids = %d.\n", numBurnableFluids);
  if (numBurnableFluids > 0) {
    programmedBurn = new ProgrammedBurn(ioData,this->X);
    this->fluidSelector.attachProgrammedBurn(programmedBurn);
  }

  setupEmbeddedFSISolver(ioData);
  setupMultiPhaseFlowSolver(ioData);

  // only for IncreasePressure
  if(ioData.implosion.type==ImplosionSetup::LINEAR)
    implosionSetupType = MultiPhysicsTsDesc<dim,dimLS>::LINEAR;
  else if(ioData.implosion.type==ImplosionSetup::SMOOTHSTEP)
    implosionSetupType = MultiPhysicsTsDesc<dim,dimLS>::SMOOTHSTEP;
  else {
    this->com->fprintf(stderr,"ERROR: ImplosionSetup::Type = %d. Code not supported!\n", ioData.implosion.type);
    implosionSetupType = MultiPhysicsTsDesc<dim,dimLS>::NONE;
    exit(-1);
  }

  Prate = ioData.implosion.Prate;
  Pinit = ioData.implosion.Pinit;
  Pscale = ioData.ref.rv.pressure;
  intersector_freq = ioData.implosion.intersector_freq;
  if(ioData.implosion.type==ImplosionSetup::LINEAR)
    tmax = (Prate == 0) ? std::numeric_limits<double>::max() : (ioData.bc.inlet.pressure - Pinit)/Prate;
  else if(ioData.implosion.type==ImplosionSetup::SMOOTHSTEP) {
    tmax = ioData.implosion.tmax;
    if(tmax<=0.0) {
      this->com->fprintf(stderr,"ERROR: ImplosionSetup::Tmax = %e!\n", tmax);
      exit(-1);
    }
  }

  if(intersector_freq<1) {
    this->com->fprintf(stderr,"ERROR: InterfaceTrackingFrequency must be larger than 0. Currently it is %d.\n", intersector_freq);
    exit(-1);
  }

  // miscellaneous
  globIt = -1;
  inSubCycling = false;
  inRedoTimestep = false;
  double *Vin = this->bcData->getInletPrimitiveState();
  for(int i=0; i<dim; i++)
    vfar[i] =Vin[i]; //for phase-change update only
  
  requireSpecialBDF = false;
  increasingPressure = false;

//  ioData.printDebug();


  multiFluidInterfaceOrder = 1;

  lastLSUpdateIteration = -1;

  if (ioData.mf.interfaceTreatment == MultiFluidData::SECONDORDER) {

    dom->createHigherOrderMultiFluid();

    multiFluidInterfaceOrder = 2;

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


  if (interfaceAlg == 1) {

    dom->createHigherOrderFSI();

#pragma omp parallel for
    for (int iSub = 0; iSub < dom->getNumLocSub(); ++iSub) {
      V6NodeData (*v6data)[2];
      v6data = 0;
      dom->getSubDomain()[iSub]->findEdgeTetrahedra((*this->X)(iSub), v6data);
      dom->getSubDomain()[iSub]->getHigherOrderFSI()->
	initialize<dim>(dom->getNodeDistInfo().subSize(iSub),
			dom->getSubDomain()[iSub]->getElems(),
			v6data); 

      //if (ioData.mf.interfaceLimiter == MultiFluidData::LIMITERALEX1)
      dom->getSubDomain()[iSub]->getHigherOrderFSI()->setLimitedExtrapolation();
    }
    
  }


  limitHigherOrderExtrapolation = false;
  if (ioData.mf.interfaceLimiter == MultiFluidData::LIMITERALEX1)
    limitHigherOrderExtrapolation = true;

  phaseChangeType = 0;
  
  if (ioData.mf.typePhaseChange == MultiFluidData::EXTRAPOLATION) {

    phaseChangeType = 1;
  }

}


//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::setupEmbeddedFSISolver(IoData &ioData)
{

  this->postOp->setForceGenerator(this);

  phaseChangeChoice  = (ioData.embed.eosChange==EmbeddedFramework::RIEMANN_SOLUTION) ? 1 : 0;
//  forceApp           = (ioData.embed.forceAlg==EmbeddedFramework::RECONSTRUCTED_SURFACE) ? 3 : 1;
  switch (ioData.embed.forceAlg) {
    case EmbeddedFramework::CONTROL_VOLUME_BOUNDARY :
      forceApp = 1;
      break;
    case EmbeddedFramework::EMBEDDED_SURFACE :
      forceApp = 2;
      break;
    case EmbeddedFramework::RECONSTRUCTED_SURFACE :
      forceApp = 3;
      break;
    default:
      this->com->fprintf(stderr,"ERROR: force approach not specified correctly! Abort...\n"); 
      exit(-1);
  }
  linRecAtInterface  = (ioData.embed.reconstruct==EmbeddedFramework::LINEAR) ? true : false;
  viscSecOrder  = (ioData.embed.viscousinterfaceorder==EmbeddedFramework::SECOND) ? true : false;
  riemannNormal = (int)ioData.embed.riemannNormal;

  if(orderOfAccuracy==1) { //first-order everywhere...
    linRecAtInterface = false;
    viscSecOrder = false;
  }

  //for phase-change update
  Weights = 0;
  VWeights = 0;

//------------- For Fluid-Structure Interaction ---------------
  withCracking = false;
  if(ioData.problem.type[ProblemData::FORCED] || ioData.problem.type[ProblemData::AERO]) {
    dynNodalTransfer = new DynamicNodalTransfer(ioData, *this->domain->getCommunicator(), *this->domain->getStrCommunicator(),
                                                this->domain->getTimer());
    withCracking = dynNodalTransfer->cracking();

    //for updating phase change
    Weights  = new DistVec<double>(this->getVecInfo());
    VWeights = new DistSVec<double,dim>(this->getVecInfo());
  } else
    dynNodalTransfer = 0;
//-------------------------------------------------------------

//-------------------- Setup Intersector ----------------------
#ifdef DO_EMBEDDED
  switch (ioData.embed.intersectorName) {
    case EmbeddedFramework::FRG : // use IntersectorFRG
      if(dynNodalTransfer && dynNodalTransfer->embeddedMeshByFEM()) {
//        this->com->fprintf(stderr,"Using dynamic nodal transfer to get embedded surface data.\n");
        int nNodes = dynNodalTransfer->numStNodes();
        int nElems = dynNodalTransfer->numStElems();

        if(withCracking) {
          this->com->fprintf(stderr,"ERROR: IntersectorFRG is not capable of handling cracking structures!\n");
          this->com->fprintf(stderr,"       Try IntersectorPhysBAM instead!\n");
          exit(-1);
        }

        double *xyz = dynNodalTransfer->getStNodes();
        int (*abc)[3] = dynNodalTransfer->getStElems();
        distLSS = new DistIntersectorFRG(ioData, this->com, nNodes, xyz, nElems, abc);
      } else
        distLSS = new DistIntersectorFRG(ioData, this->com);
      break;
    case EmbeddedFramework::PHYSBAM : // use IntersectorPhysBAM
      if(dynNodalTransfer && dynNodalTransfer->embeddedMeshByFEM()) {
//        this->com->fprintf(stderr,"Using dynamic nodal transfer to get embedded surface data.\n");
        int nNodes = dynNodalTransfer->numStNodes();
        int nElems = dynNodalTransfer->numStElems();
        double *xyz = dynNodalTransfer->getStNodes();
        int (*abc)[3] = dynNodalTransfer->getStElems();

        if(withCracking) {
          this->com->fprintf(stderr,"Note: Topology change of the embedded surface will be considered...\n");
          distLSS = new DistIntersectorPhysBAM(ioData, this->com, nNodes, xyz, nElems, abc, dynNodalTransfer->getCrackingSurface());
        } else
          distLSS = new DistIntersectorPhysBAM(ioData, this->com, nNodes, xyz, nElems, abc);

      } else
        distLSS = new DistIntersectorPhysBAM(ioData, this->com);
      break;
    default:
      this->com->fprintf(stderr,"ERROR: No valid intersector specified! Check input file\n");
      exit(-1);
  }
#else
  this->com->fprintf(stderr,"ERROR: Embedded framework is NOT compiled! Check your makefile.\n");
  exit(-1);
#endif
//-------------------------------------------------------------

  //Riemann solution stored on edges
  Wstarij = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  Wstarji = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  *Wstarij = 0.0;
  *Wstarji = 0.0;

  if (this->timeState->useNm1()) {
    Wstarij_nm1 = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
    Wstarji_nm1 = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
    *Wstarij_nm1 = 0.0;
    *Wstarji_nm1 = 0.0;
  } else {
    Wstarij_nm1 = 0;
    Wstarji_nm1 = 0;
  }

//------ load structure mesh information ----------------------
  Fs = 0;
  numStructNodes = distLSS->getNumStructNodes();
  totStructNodes = dynNodalTransfer ? dynNodalTransfer->totStNodes() : numStructNodes;

  if (numStructNodes>0) {
    //this->com->fprintf(stderr,"- Embedded Structure Surface: %d (%d) nodes.\n", numStructNodes, totStructNodes);
    // We allocate Fs from memory that allows fast one-sided MPI communication
    Fs = new (*this->com) double[totStructNodes][3];
  } else
    this->com->fprintf(stderr,"Warning: failed loading structure mesh information!\n");

  FsComputed = false;
//-------------------------------------------------------------

//------ Setup ghostPoints for Viscous Flows ------------------
  ghostPoints = 0;
  switch (ioData.eqs.type)
    {
    case EquationsData::EULER:
      eqsType = MultiPhysicsTsDesc<dim,dimLS>::EULER;
      break;
    case EquationsData::NAVIER_STOKES:
      eqsType = MultiPhysicsTsDesc<dim,dimLS>::NAVIER_STOKES;
      ghostPoints  = new DistVec<GhostPoint<dim> *>(this->getVecInfo());
      ghostPoints->nullifyPointers();
      break;
    }
//-------------------------------------------------------------

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::setupMultiPhaseFlowSolver(IoData &ioData)
{
  LS = new LevelSet<dimLS>(ioData, this->domain);
  frequencyLS = ioData.mf.frequency;
  withMixedLS = (ioData.embed.crackingWithLevelset==EmbeddedFramework::ON);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
MultiPhysicsTsDesc<dim,dimLS>::~MultiPhysicsTsDesc()
{
  if (riemann) delete riemann;

  // Embedded FSI
  if (distLSS) delete distLSS;
  if (Wstarij) delete Wstarij;
  if (Wstarji) delete Wstarji;
  if (Wstarij_nm1) delete Wstarij_nm1;
  if (Wstarji_nm1) delete Wstarji_nm1;
  if (Weights) delete Weights;
  if (VWeights) delete VWeights;
  if (dynNodalTransfer) delete dynNodalTransfer;
  if (Fs) operator delete[] (Fs, *this->com);
  if(ghostPoints) {
    ghostPoints->deletePointers();
    delete ghostPoints;
  }

  // Level-set
  if (LS) delete LS;
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::setupTimeStepping(DistSVec<double,dim> *U, IoData &ioData)
{
  // Setup fluid mesh geometry
  this->geoState->setup2(this->timeState->getData());
  // Initialize intersector and compute intersections
  DistVec<int> point_based_id(this->domain->getNodeDistInfo());

  // If we have the fluid ID from a restart file, use that.  Otherwise,
  // it will be set in initialize()
  // Added by Alex Main (September 2013)
  //
  if (ioData.input.fluidId[0] != 0 || ioData.input.restart_file_package[0] != 0) {
    distLSS->initialize(this->domain,*this->X, this->geoState->getXn(), ioData, &point_based_id,
                        fluidSelector.fluidId);
  } else
    distLSS->initialize(this->domain,*this->X, this->geoState->getXn(), ioData, &point_based_id);

  // Initialize fluid state vector
  this->timeState->setup(this->input->solutions, *this->X, this->bcData->getInletBoundaryVector(),
                         *U, ioData, &point_based_id); //populate U by i.c. or restart data.

  *(this->Xs) = *(this->X);

  this->initializeFarfieldCoeffs();

  // Initialize level-sets 
  if(withCracking && withMixedLS) 
    LS->setup(this->input->levelsets, *this->X, *U, Phi, ioData, &fluidSelector, this->varFcn, 
              &(distLSS->getClosestPoints()), &(distLSS->getStatus()),lsMethod);
  else
    LS->setup(this->input->levelsets, *this->X, *U, Phi, ioData, &fluidSelector, this->varFcn, 0, 0, lsMethod);

  // Initialize or reinitialize (i.e. at the beginning of a restart) fluid Ids
  //if(this->input->levelsets[0] == 0) // init
  //  fluidSelector.initializeFluidIds(distLSS->getStatus(), *this->X, ioData);
  //else //restart
  /*if(withCracking && withMixedLS) 
    fluidSelector.reinitializeFluidIdsWithCracking(distLSS->getStatus(), Phi); 
  else 
    fluidSelector.reinitializeFluidIds(distLSS->getStatus(), Phi);
  */

  // If we have the fluid ID from a restart file, use that.
  // Added by Alex Main (September 2013)
  //
  if(ioData.input.fluidId[0] == 0 && ioData.input.restart_file_package[0] == 0) // init
    fluidSelector.reinitializeFluidIds(distLSS->getStatus(), Phi);
 
  if (programmedBurn)
    programmedBurn->setFluidIds(this->getInitialTime(), *(fluidSelector.fluidId), *U);

  // Initialize the embedded FSI handler
  EmbeddedMeshMotionHandler* _mmh = dynamic_cast<EmbeddedMeshMotionHandler*>(this->mmh);
  if(_mmh) {
    double *tMax = &(this->data)->maxTime;
    _mmh->setup(tMax); //obtain maxTime from structure
  }

  // If 'IncreasePressure' is activated, re-initialize the fluid state
  if(Pinit>=0.0 && Prate>=0.0 && this->getInitialTime()<tmax) {
    increasingPressure = true;
    this->domain->IncreasePressure(currentPressure(this->getInitialTime()), this->varFcn, *U, *(this->fluidSelector.fluidId));
  }

  //compute fluid indueced force on the structural surface
  // construct Wij, Wji from U. Then use them for force calculation.
  DistSVec<double,dim> *Wij = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  DistSVec<double,dim> *Wji = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  DistSVec<double,dim> VV(this->getVecInfo());
  *Wij = 0.0;
  *Wji = 0.0;
  this->varFcn->conservativeToPrimitive(*U,VV,fluidSelector.fluidId);
  SubDomain **subD = this->domain->getSubDomain();

  int iSub;
#pragma omp parallel for
  for (iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
    SVec<double,dim> &subWij = (*Wij)(iSub);
    SVec<double,dim> &subWji = (*Wji)(iSub);
    SVec<double,dim> &subVV = VV(iSub);
    int (*ptr)[2] =  (subD[iSub]->getEdges()).getPtr();
    if (subWij.size()!=(subD[iSub]->getEdges()).size()) {fprintf(stderr,"WRONG!!!\n"); exit(-1);}

    for (int l=0; l<subWij.size(); l++) {
      int i = ptr[l][0];
      int j = ptr[l][1];
      for (int k=0; k<dim; k++) {
        subWij[l][k] = subVV[i][k];
        subWji[l][k] = subVV[j][k];
      }
    }
  }

  this->spaceOp->conservativeToPrimitive(*U,fluidSelector.fluidId);
  computeForceLoad(Wij, Wji);
  delete Wij;
  delete Wji;

  FsComputed = true;
  // Now "accumulate" the force for the embedded structure
  if(dynNodalTransfer){
    numStructNodes = dynNodalTransfer->numStNodes();
    SVec<double,3> v(numStructNodes, Fs);
    dynNodalTransfer->updateOutputToStructure(0.0, 0.0, v); //dt=dtLeft=0.0-->They are not used!
  }

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
double MultiPhysicsTsDesc<dim,dimLS>::computeTimeStep(int it, double *dtLeft,
                                                  DistSVec<double,dim> &U, double angle)
{
  if(!FsComputed&&dynNodalTransfer) this->com->fprintf(stderr,"WARNING: FSI force not computed!\n");
  FsComputed = false; //reset FsComputed at the beginning of a fluid iteration

  //check if it's in subcycling with iCycle>1.
  if(globIt==it)
    inSubCycling = true;
  else {
    globIt = it;
    inSubCycling = false;
  }

  double t0 = this->timer->getTime();
  this->data->computeCflNumber(it - 1, this->data->residual / this->restart->residual, angle);
  int numSubCycles = 1;

  double dt;
  umax = 0.0; 
  
  dt = this->timeState->computeTimeStep(this->data->cfl, this->data->dualtimecfl, dtLeft,
                          &numSubCycles, *this->geoState, *this->A, U, *(fluidSelector.fluidId),&umax);

  if (this->problemType[ProblemData::UNSTEADY])
    this->com->printf(5, "Global dt: %g (remaining subcycles = %d)\n",
                      dt*this->refVal->time, numSubCycles);

  dtf = dt;
  dtfLeft = *dtLeft + dt;

  this->timer->addFluidSolutionTime(t0);
  this->timer->addTimeStepTime(t0);

  return dt;
}

//---------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::updateStateVectors(DistSVec<double,dim> &U, int it)
{
  this->geoState->update(*this->X, *this->A);
  
//  if(withCracking && withMixedLS)
//    fluidSelector.checkLSConsistency(Phi);

	  
  // Speed improvement.  Reinitialization is relatively expensive.
  // so only do it if we haven't already done it this iteration.
  // This is effective when we are doing fluid subcycling.
  // Added by Alex Main (October 2013)
  //
  if(frequencyLS > 0 && it%frequencyLS == 0  && lastLSUpdateIteration != it){
    lastLSUpdateIteration = it;
    if (this->lsMethod == 0)
      LS->conservativeToPrimitive(Phi,PhiV,U);
    else
      PhiV = Phi;
    if(withCracking && withMixedLS) {
      this->multiPhaseSpaceOp->resetFirstLayerLevelSetFS(PhiV, this->distLSS, *fluidSelector.fluidId, InterfaceTag);
      for (int k = 0; k < dimLS; ++k)
      	LS->reinitializeLevelSet(*this->X, PhiV, false,k);
    } else {
      for (int k = 0; k < dimLS; ++k)
	LS->reinitializeLevelSet(*this->X, PhiV,true, k);
    }
    if (this->lsMethod == 0)
      LS->primitiveToConservative(PhiV,Phi,U);
    else
      Phi = PhiV;
    LS->update(Phi);
    if (this->timeState->useNm1()) {
      DistSVec<double,dimLS>& Phinm1 = LS->getPhinm1();
      this->multiPhaseSpaceOp->computeResidualLS(*this->X, *this->A, Phi, *fluidSelector.fluidId, U, Phinm1, distLSS, linRecAtInterface);
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

  this->timeState->update(U, U, *(fluidSelector.fluidIdn), fluidSelector.fluidIdnm1, riemann,distLSS,increasingPressure);
                            //fluidIdn, fluidIdnm1 and riemann are used only for implicit time-integrators

  this->spaceOp->updateFixes();
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
int MultiPhysicsTsDesc<dim,dimLS>::checkSolution(DistSVec<double,dim> &U)
{
  int ierr = this->domain->checkSolution(this->varFcn, *this->A, U, *fluidSelector.fluidId, *fluidSelector.fluidIdn);
                             // fluidIdn is only used for screen output when an error is found
  return ierr;
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::setupOutputToDisk(IoData &ioData, bool *lastIt,
                          int it, double t, DistSVec<double,dim> &U)
{
//  conservationErrors(U,it); /* not checking conservation errors for the moment...*/

  if (it == this->data->maxIts)
    *lastIt = true;
  else
    monitorInitialState(it, U); // Phi?

//  this->output->setMeshMotionHandler(ioData, this->mmh); /*currently not supported by embedded framework.*/

  this->output->openAsciiFiles();
  this->timer->setSetupTime();
  this->output->cleanProbesFile();

  double fluxNorm = 0.5*(this->data->residual)*(this->data->residual);

  if (it == 0) {
    this->multiPhaseSpaceOp->computeGradP(*this->X, *this->A, U); /*really used???*/
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
    this->output->writeEmbeddedSurfaceToDisk(*lastIt, it, t, distLSS->getStructPosition_n(), distLSS->getStructPosition_0());
    this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, *fluidSelector.fluidId, &Phi);
//    this->output->writeHeatFluxesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, fluidSelector.fluidId);
//    this->output->writeAvgVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);
  }
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::outputToDisk(IoData &ioData, bool* lastIt, int it,
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
  this->output->writeEmbeddedSurfaceToDisk(*lastIt, it, t, distLSS->getStructPosition_n(), distLSS->getStructPosition_0());
  this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, *fluidSelector.fluidId, &Phi);
//  this->output->writeAvgVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);

  this->output->writeProbesToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState,*fluidSelector.fluidId,&Phi);
  this->restart->writeToDisk(this->com->cpuNum(), *lastIt, it, t, dt, *this->timeState, *this->geoState, LS, dynNodalTransfer,&fluidSelector);
  this->restart->writeStructPosToDisk(this->com->cpuNum(), *lastIt, distLSS->getStructPosition_n()); //KW: must be after writeToDisk

  this->output->updatePrtout(t);
  this->restart->updatePrtout(t);
  if (*lastIt) {
    this->timer->setRunTime();
    if (this->com->getMaxVerbose() >= 2)
      this->timer->print(this->domain->getStrTimer());
    this->output->closeAsciiFiles();
  }
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::outputForces(IoData &ioData, bool* lastIt, int it, int itSc, int itNl,
                               double t, double dt, DistSVec<double,dim> &U)
{
  double cpu = this->timer->getRunTime();
  this->output->writeForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, fluidSelector.fluidId);

  // PJSA this should written before getting the new cracking data
  this->restart->writeCrackingDataToDisk(this->com->cpuNum(), *lastIt, it, t, dynNodalTransfer);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::outputPositionVectorToDisk(DistSVec<double,dim> &U)
{/* currently only consider the embedded framework. No mesh motion. */}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::resetOutputToStructure(DistSVec<double,dim> &U)
{
  /* currently only consider the embedded framework. No mesh motion. */
/*  this->com->printf(5,"LevelSetTsDesc<dim,dimLS>::resetOutputToStructure\n");

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(this->mmh);
  if (_mmh)
    _mmh->resetOutputToStructure(this->postOp, *this->X, U, fluidSelector.fluidId);
*/
}

//------------------------------------------------------------------------------

template <int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::updateOutputToStructure(double dt, double dtLeft, DistSVec<double,dim> &U)
{
  if (this->mmh) {
    double work[2];
    this->mmh->computeInterfaceWork(dt, this->postOp, this->geoState->getXn(),
                                    this->timeState->getUn(), *this->X, U,
                                    work, fluidSelector.fluidIdn, fluidSelector.fluidId); //FF interface
    this->restart->energy[0] += work[0];
    this->restart->energy[1] += work[1];
  }

  if(dynNodalTransfer) {
    computeForceLoad(this->Wstarij, this->Wstarji);
    FsComputed = true; //to avoid redundant computation of Fs.

    // Now "accumulate" the force for the embedded structure
    numStructNodes = dynNodalTransfer->numStNodes();
    SVec<double,3> v(numStructNodes, Fs);
    dynNodalTransfer->updateOutputToStructure(dt, dtLeft, v);
  }
}

//-------------------------------------------------------------------------------

template<int dim, int dimLS>
double MultiPhysicsTsDesc<dim,dimLS>::computeResidualNorm(DistSVec<double,dim>& U)
{
//  // Ghost-Points Population (for Navier-Stokes only)
//  if(this->eqsType == MultiPhysicsTsDesc<dim,dimLS>::NAVIER_STOKES) {
//    this->ghostPoints->deletePointers();
//    this->multiPhaseSpaceOp->populateGhostPoints(this->ghostPoints,*this->X,U,this->varFcn,this->distLSS,this->viscSecOrder,*(this->fluidSelector.fluidId));
//  }
  if (this->lsMethod == 0)
    LS->conservativeToPrimitive(Phi,PhiV,U);
  else
    PhiV = Phi;
  this->multiPhaseSpaceOp->computeResidual(*this->X, *this->A, U, *Wstarij, *Wstarji, distLSS, linRecAtInterface, viscSecOrder, this->riemann, riemannNormal, PhiV, fluidSelector, *this->R, 0, 0);

  this->multiPhaseSpaceOp->applyBCsToResidual(U, *this->R);

  double res = 0.0;
  if(this->numFluid==1)
    res = this->multiPhaseSpaceOp->computeRealFluidResidual(*this->R, *this->Rreal, *distLSS);
  else {
    this->com->fprintf(stderr,"WARNING: MultiPhysicsTsDesc::computeResidualNorm is not implemented for numFluid > 1\n");
    res = 1.0;
  }

  return sqrt(res);
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::monitorInitialState(int it, DistSVec<double,dim> &U)
{
  /* only used for steady simulations. Is it meaningful when a ff interface is present */
  //this->com->printf(2, "State vector norm = %.12e\n", sqrt(U*U));
  if (!this->problemType[ProblemData::UNSTEADY]) {
    double trhs = this->timer->getTimeSyncro();
    this->data->residual = computeResidualNorm(U);
    trhs = this->timer->getTimeSyncro() - trhs;
    if (it == 0)
      this->restart->residual = this->data->residual;
    if (this->data->resType == -1)
      this->com->printf(2, "Spatial residual norm = %.12e\n", this->data->residual);
    else
      this->com->printf(2, "Spatial residual norm[%d] = %.12e\n", this->data->resType, this->data->residual);
    this->com->printf(2, "Time for one residual evaluation: %f s\n", trhs);
  }

  this->com->printf(2, "\n");

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::computeForceLoad(DistSVec<double,dim> *Wij, DistSVec<double,dim> *Wji)
{
  if (!Fs) {fprintf(stderr,"computeForceLoad: Fs not initialized! Cannot compute the load!\n"); return;}
  double t0 = this->timer->getTime();
  if(dynNodalTransfer)
    numStructNodes = dynNodalTransfer->numStNodes();
  for (int i=0; i<numStructNodes; i++)
    Fs[i][0] = Fs[i][1] = Fs[i][2] = 0.0;
  this->multiPhaseSpaceOp->computeForceLoad(forceApp, orderOfAccuracy, *this->X, *this->A, Fs, numStructNodes, distLSS, *Wij, *Wji, 0, 0, 
														  this->postOp->getPostFcn(),fluidSelector.fluidId); 
// 0 is the pointer to the GhostPoint structure. May be needed in the future. 
  this->timer->addEmbeddedForceTime(t0);
  //at this stage Fs is NOT globally assembled!
}

//-------------------------------------------------------------------------------

template <int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::getForcesAndMoments(map<int,int> & surfOutMap, DistSVec<double,dim> &U, DistSVec<double,3> &X,
                                           Vec3D *Fi, Vec3D *Mi)
{
  int idx;
  if (!FsComputed) 
    computeForceLoad(this->Wstarij, this->Wstarji);

  if(dynNodalTransfer)
    numStructNodes = dynNodalTransfer->numStNodes();

  Vec<Vec3D>& Xstruc = distLSS->getStructPosition();

  for (int i=0; i<numStructNodes; i++) {
    map<int,int>::iterator it = surfOutMap.find(distLSS->getSurfaceID(i));
    if(it != surfOutMap.end() && it->second != -2)
      idx = it->second;
    else {
      idx = 0;
    }

    Fi[idx][0] += Fs[i][0]; 
    Fi[idx][1] += Fs[i][1]; 
    Fi[idx][2] += Fs[i][2];

    Mi[idx][0] += Xstruc[i][1]*Fs[i][2]-Xstruc[i][2]*Fs[i][1];
    Mi[idx][1] += Xstruc[i][2]*Fs[i][0]-Xstruc[i][0]*Fs[i][2];
    Mi[idx][2] += Xstruc[i][0]*Fs[i][1]-Xstruc[i][1]*Fs[i][0];
  }
}

//-------------------------------------------------------------------------------
template <int dim, int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::getderivativeOfForcesAndMoments(map<int,int> & surfOutMap, 
								    DistSVec<double,dim> &U, DistSVec<double,dim> &dU, 
								    DistSVec<double,3> &X, double dS[3],
								    Vec3D *dFi, Vec3D *dMi) 
{
  
  fprintf(stderr, "ERROR: MultiPhys getderivativeOfForcesAndMoments not implemented \n");

}
//-------------------------------------------------------------------------------
template<int dim, int dimLS>
bool MultiPhysicsTsDesc<dim,dimLS>::IncreasePressure(int it, double dt, double t, DistSVec<double,dim> &U)
{

  increasingPressure = false;
  if(Pinit<0.0 || Prate<0.0) return true; // no setup for increasing pressure

  if(t>tmax && t-dt>tmax) {// max pressure was reached, so now we solve
    return true;
  }

  increasingPressure = true;

  // max pressure not reached, so we do not solve and we increase pressure and let structure react

  if(this->mmh && !inSubCycling) {
    //get structure timestep dts
    this->dts = this->mmh->update(0, 0, 0, this->bcData->getVelocityVector(), *this->Xs);
    //recompute intersections
    double tw = this->timer->getTime();
    if(intersector_freq==1||((it-1)%intersector_freq==0&&(it>1))) {
      this->com->fprintf(stderr,"recomputing fluid-structure intersections.\n");
      this->distLSS->recompute(this->dtf, this->dtfLeft, this->dts, true, TsDesc<dim>::failSafeFlag); 
    }

    this->timer->addIntersectionTime(tw);
    this->timer->removeIntersAndPhaseChange(tw);
    //updateFluidIdFS
    if (this->lsMethod == 0)
      LS->conservativeToPrimitive(Phi,PhiV,U);
    else
      PhiV = Phi;
    //this->multiPhaseSpaceOp->extrapolatePhiV(this->distLSS, this->PhiV);
    this->fluidSelector.updateFluidIdFS(this->distLSS, this->PhiV);
    this->PhiV = 0.0; //PhiV is no longer a distance function now. Only its sign (+/-)
                      //  is meaningful. We destroy it so people wouldn't use it
                      //  by mistake later on.
    //phase-change update
    tw = this->timer->getTime();
    this->multiPhaseSpaceOp->updateSweptNodes(*this->X, this->phaseChangeChoice, U, this->Vtemp, *this->Weights, *this->VWeights,
                                              this->Phi, this->PhiWeights, *this->Wstarij, *this->Wstarji,
                                              this->distLSS, this->vfar, false, this->fluidSelector.fluidIdn, this->fluidSelector.fluidId);
 
    this->timer->addEmbedPhaseChangeTime(tw);
    this->timer->removeIntersAndPhaseChange(tw);
  }

  // construct Wij, Wji from U. 
  DistSVec<double,dim> VV(this->getVecInfo());
  this->varFcn->conservativeToPrimitive(U,VV,this->fluidSelector.fluidId);
  SubDomain **subD = this->domain->getSubDomain();

  int iSub;
#pragma omp parallel for
  for (iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {
    SVec<double,dim> &subWstarij = (*Wstarij)(iSub);
    SVec<double,dim> &subWstarji = (*Wstarji)(iSub);
    SVec<double,dim> &subVV = VV(iSub);
    int (*ptr)[2] =  (subD[iSub]->getEdges()).getPtr();

    for (int l=0; l<subWstarij.size(); l++) {
      int i = ptr[l][0];
      int j = ptr[l][1];
      for (int k=0; k<dim; k++) {
        subWstarij[l][k] = subVV[i][k];
        subWstarji[l][k] = subVV[j][k];
      }
    }
  }

  // Population of spaceOp->V for the force computation
  this->spaceOp->conservativeToPrimitive(U, this->fluidSelector.fluidId); // PJSA

  double pnow = currentPressure(t);
  this->com->fprintf(stdout, "about to increase pressure to %e\n", pnow*Pscale);
  this->domain->IncreasePressure(pnow, this->varFcn, U, *(this->fluidSelector.fluidId));

  return false;

}
//-------------------------------------------------------------------------------

template<int dim,int dimLS>
void MultiPhysicsTsDesc<dim,dimLS>::setCurrentTime(double t,DistSVec<double,dim>& U) { 

  currentTime = t;

  if (programmedBurn)
    programmedBurn->setCurrentTime(t,multiPhaseSpaceOp->getVarFcn(), U,*(fluidSelector.fluidId),*(fluidSelector.fluidIdn));
}

//-------------------------------------------------------------------------------

template<int dim,int dimLS>
double MultiPhysicsTsDesc<dim,dimLS>::currentPressure(double t)
{
  double p;
  if(implosionSetupType==MultiPhysicsTsDesc<dim,dimLS>::LINEAR)
    p = Pinit + t*Prate;
  else if(implosionSetupType==MultiPhysicsTsDesc<dim,dimLS>::SMOOTHSTEP) {
    double tbar = t/tmax;
    p = Pinit + (Pfinal - Pinit)*tbar*tbar*tbar*(10.0-15.0*tbar+6.0*tbar*tbar);
  } else {
    this->com->fprintf(stderr,"ERROR! ImplosionSetup::Type = %d NOT recognized!\n", implosionSetupType);
    exit(-1);
  }
  return p;
}
