#include <EmbeddedTsDesc.h>
#include <DistExactRiemannSolver.h>
#include <FSI/DynamicNodalTransfer.h>
#include <FSI/CrackingSurface.h>
#include <Domain.h>

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

template<int dim>
EmbeddedTsDesc<dim>::
EmbeddedTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  TsDesc<dim>(ioData, geoSource, dom), nodeTag(this->getVecInfo()), nodeTag0(this->getVecInfo()),
  Vtemp(this->getVecInfo()), numFluid(ioData.eqs.numPhase), Wtemp(this->getVecInfo()),
  ioData(ioData)
{

  currentTime = 0.0;

  simType         = (ioData.problem.type[ProblemData::UNSTEADY]) ? 1 : 0;
  orderOfAccuracy = (ioData.schemes.ns.reconstruction == SchemeData::CONSTANT) ? 1 : 2;

  this->postOp->setForceGenerator(this);

  phaseChangeChoice  = (ioData.embed.eosChange==EmbeddedFramework::RIEMANN_SOLUTION) ? 1 : 0;
  phaseChangeAlg	 = (ioData.embed.phaseChangeAlg==EmbeddedFramework::LEAST_SQUARES) ? 1 : 0;
  interfaceAlg		 = (ioData.embed.interfaceAlg==EmbeddedFramework::INTERSECTION) ? 1 : 0;
  intersectAlpha	 = ioData.embed.alpha;
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
      
	//first-order everywhere... //d2d
	if(orderOfAccuracy == 1) 
	{
    linRecAtInterface = false; 
    viscSecOrder = false; 
  }

	if(interfaceAlg == 1 || ioData.embed.surrogateinterface == EmbeddedFramework::EXTERNAL) 
	{
    dom->createHigherOrderFSI();

		if(ioData.embed.surrogateinterface == EmbeddedFramework::HYBRID) 
		{
#pragma omp parallel for
			for(int iSub = 0; iSub < dom->getNumLocSub(); ++iSub) 
			{				
				V6NodeData (*v6data)[2]; v6data = 0;

      dom->getSubDomain()[iSub]->findEdgeTetrahedra((*this->X)(iSub), v6data);
			 
				dom->getSubDomain()[iSub]->getHigherOrderFSI()->initialize<dim>(dom->getNodeDistInfo().subSize(iSub),
			dom->getSubDomain()[iSub]->getElems(),
			v6data); 

      if (ioData.embed.interfaceLimiter == EmbeddedFramework::LIMITERALEX1)
        dom->getSubDomain()[iSub]->getHigherOrderFSI()->setLimitedExtrapolation();
    }
    
  }
		else if(ioData.embed.surrogateinterface == EmbeddedFramework::EXTERNAL) 
		{		 
			bool hmode = (interfaceAlg == 1) ? true : false;

#pragma omp parallel for
			for (int iSub = 0; iSub < dom->getNumLocSub(); ++iSub)
				dom->getSubDomain()[iSub]->getHigherOrderFSI()->initialize<dim>(ioData, this->com, dom->getNodeDistInfo().subSize(iSub),
																									 dom->getSubDomain()[iSub]->getElems(),
																									 hmode, viscSecOrder);
		}
		
	}

  this->timeState = new DistTimeState<dim>(ioData, this->spaceOp, this->varFcn, this->domain, this->V);

  riemann = new DistExactRiemannSolver<dim>(ioData,this->domain,this->varFcn);

  //for phase-change update
  Weights = 0;
  VWeights = 0;

//------------- For Fluid-Structure Interaction ---------------
  withCracking = false;
	if(ioData.problem.type[ProblemData::FORCED] || ioData.problem.type[ProblemData::AERO]) 
	{
    dynNodalTransfer = new DynamicNodalTransfer(ioData, *this->domain->getCommunicator(), *this->domain->getStrCommunicator(),
                                                this->domain->getTimer());
    withCracking = dynNodalTransfer->cracking(); 

    //for updating phase change
    Weights  = new DistVec<double>(this->getVecInfo());
    VWeights = new DistSVec<double,dim>(this->getVecInfo());
	} 
	else
    dynNodalTransfer = 0;

  emmh = 0;
//-------------------------------------------------------------

#ifdef DO_EMBEDDED
	switch (ioData.embed.intersectorName) 
	{
    case EmbeddedFramework::FRG :
      if(dynNodalTransfer && dynNodalTransfer->embeddedMeshByFEM()) {
      
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

    case EmbeddedFramework::PHYSBAM : 
      if(dynNodalTransfer && dynNodalTransfer->embeddedMeshByFEM()) {
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
  wall_computer=new ReinitializeDistanceToWall<1>(ioData, *this->domain);

#else
  this->com->fprintf(stderr,"ERROR: Embedded framework is NOT compiled! Check your makefile.\n");
  exit(-1);
#endif

  //Riemann solution stored on edges
  Wstarij = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  Wstarji = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  *Wstarij = 0.0;
  *Wstarji = 0.0;

  //linRecAtInterface  = (ioData.embed.reconstruct==EmbeddedFramework::LINEAR) ? true : false;
  if (interfaceAlg) {
    countWstarij = new DistVec<int>(this->domain->getEdgeDistInfo());
    countWstarji = new DistVec<int>(this->domain->getEdgeDistInfo());
    *countWstarij = 0;
    *countWstarji = 0;
  }
  else {
    countWstarij = NULL;
    countWstarji = NULL;
  }

  //copies for fail safe
  WstarijCopy = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  WstarjiCopy = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  *WstarijCopy = 0.0;
  *WstarjiCopy = 0.0;



  if (this->timeState->useNm1()) {
    Wstarij_nm1 = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
    Wstarji_nm1 = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
    *Wstarij_nm1 = 0.0;
    *Wstarji_nm1 = 0.0;

    Wstarij_nm1Copy = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
    Wstarji_nm1Copy = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
    *Wstarij_nm1Copy = 0.0;
    *Wstarji_nm1Copy = 0.0;
  } else {
    Wstarij_nm1 = 0;
    Wstarji_nm1 = 0;

    Wstarij_nm1Copy = 0;
    Wstarji_nm1Copy = 0;
  }
  UCopy = new DistSVec<double,dim>(this->domain->getNodeDistInfo());

	Wextij = new DistSVec<double,dim>(this->domain->getNodeDistInfo());
	*Wextij = 0.0;

  //TODO: should be merged with fluidId in TsDesc
  nodeTag0 = 0;
  nodeTag = 0;

 nodeTagCopy = new DistVec<int>(this->getVecInfo());
 *nodeTagCopy = 0;
 nodeTag0Copy = new DistVec<int>(this->getVecInfo());
 *nodeTag0Copy = 0;


//---------------------------------------------------------------------
  // for IncreasePressure
  if(ioData.implosion.type==ImplosionSetup::LINEAR)
    implosionSetupType = EmbeddedTsDesc<dim>::LINEAR;
  else if(ioData.implosion.type==ImplosionSetup::SMOOTHSTEP)
    implosionSetupType = EmbeddedTsDesc<dim>::SMOOTHSTEP;
  else {
    this->com->fprintf(stderr,"ERROR: ImplosionSetup::Type = %d. Code not supported!\n", ioData.implosion.type);
    implosionSetupType = EmbeddedTsDesc<dim>::NONE;
    exit(-1);
  }

  Prate = ioData.implosion.Prate;
  Pinit = ioData.implosion.Pinit;
  Pfinal = ioData.bc.inlet.pressure;

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

  recomputeIntersections = true;
  unifPressure[0] = unifPressure[1] = Pinit;
//---------------------------------------------------------------------

  globIt = -1;
  inSubCycling = false;

  //store farfield state for phase-change update for fluid-fullbody
  double *Vin = this->bcData->getInletPrimitiveState();
  for(int i=0; i<dim; i++)
    vfar[i] =Vin[i];

//------ load structure mesh information ----------------------
  Fs = 0;
  numStructNodes = distLSS->getNumStructNodes();
  totStructNodes = dynNodalTransfer ? dynNodalTransfer->totStNodes() : numStructNodes;

  if (numStructNodes>0) {
    //this->com->fprintf(stderr,"- Embedded Structure Surface: %d (%d) nodes\n", numStructNodes, totStructNodes);
    // We allocate Fs from memory that allows fast one-sided MPI communication
    Fs = new (*this->com) double[totStructNodes][3];
  } else 
    this->com->fprintf(stderr,"Warning: failed loading structure mesh information!\n");

  FsComputed = false;

  dFs = 0;
  if(ioData.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_ ||
     ioData.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_ ||
	 ioData.problem.alltype == ProblemData::_SENSITIVITY_ANALYSIS_)
    dFs = new (*this->com) double[totStructNodes][3];

//-------------------------------------------------------------

  // Adam 04/06/2010
  ghostPoints = 0;
  switch (ioData.eqs.type)
    {
    case EquationsData::EULER:
      eqsType = EmbeddedTsDesc<dim>::EULER;
      break;
    case EquationsData::NAVIER_STOKES:
      eqsType = EmbeddedTsDesc<dim>::NAVIER_STOKES;
      ghostPoints  = new DistVec<GhostPoint<dim> *>(this->getVecInfo());
      ghostPoints->nullifyPointers();
      break;
    }

  increasingPressure = false;

}


//------------------------------------------------------------------------------

template<int dim>
EmbeddedTsDesc<dim>::~EmbeddedTsDesc()
{
  if (distLSS) delete distLSS;
  if (riemann) delete riemann;
  if (Wstarij) delete Wstarij;
  if (Wstarji) delete Wstarji;
  if (countWstarij) delete countWstarij;
  if (countWstarji) delete countWstarji;
  if (Wstarij_nm1) delete Wstarij_nm1;
  if (Wstarji_nm1) delete Wstarji_nm1;
  if (WstarijCopy) delete WstarijCopy;
  if (WstarjiCopy) delete WstarjiCopy;
  if (Wstarij_nm1Copy) delete Wstarij_nm1Copy;
  if (Wstarji_nm1Copy) delete Wstarji_nm1Copy;
  if (UCopy) delete UCopy;
  if (Weights) delete Weights;
  if (VWeights) delete VWeights;
  delete wall_computer;

  if (dynNodalTransfer) delete dynNodalTransfer;
  if (emmh) delete emmh;
  if (Fs) operator delete[] (Fs, *this->com);
  if (dFs) delete [] dFs;
  
  if(ghostPoints) 
    {
      ghostPoints->deletePointers();
      delete ghostPoints;
    }

  if (nodeTagCopy) delete nodeTagCopy;
  if (nodeTag0Copy) delete nodeTag0Copy;

  if(Wextij) delete Wextij;
}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::setupTimeStepping(DistSVec<double,dim> *U, IoData &ioData)
{
  // Setup fluid mesh geometry
  this->geoState->setup2(this->timeState->getData());
  // Initialize intersector and compute intersections
  DistVec<int> point_based_id(this->domain->getNodeDistInfo());

	if(ioData.input.fluidId[0] != 0 || ioData.input.restart_file_package[0] != 0) 
	{
    FluidSelector f(2, ioData,this->domain);
    nodeTag0 = nodeTag = *f.fluidId;
    distLSS->initialize(this->domain,*this->X, this->geoState->getXn(), ioData, &point_based_id, &nodeTag);
	} 
	else 
    distLSS->initialize(this->domain,*this->X, this->geoState->getXn(), ioData, &point_based_id);

   //d2d
	if(ioData.embed.surrogateinterface == EmbeddedFramework::EXTERNAL) 
	{
		this->spaceOp->setSIstencil(*this->X, distLSS, this->nodeTag, *U); 

		if(eqsType == EmbeddedTsDesc<dim>::NAVIER_STOKES)
			this->spaceOp->setFEMstencil(*this->X, distLSS, this->nodeTag, *U);
	}

  // Initialize fluid state vector
  this->timeState->setup(this->input->solutions, *this->X, this->bcData->getInletBoundaryVector(),
                         *U, ioData, &point_based_id); //populate U by i.c. or restart data.

/*////////////////////////// debug
  int isubd;
#pragma omp parallel for
  for (isubd=0; isubd<this->domain->getNumLocSub(); ++isubd) 
  {
	  int lnu = (*U)(isubd).size();
	  for(int k_=0; k_<lnu; ++k_)
	  {
		  if( !(*distLSS)(isubd).isActive(0.0, k_) ) (*U)(isubd)[k_][1] = 0.0;
		  // if( fabs( (*this->X)(isubd)[k_][1] ) <= 0.0105 && 
		  // 		(*this->X)(isubd)[k_][0] >= 0.1 && 
		  // 		(*this->X)(isubd)[k_][0] <= 0.97  )
		  // {
		  // 	  (*U)(isubd)[k_][1] = 0.0;
		  // }
	  } 
  }
/////////////////////////// test */

  this->spaceOp->applyBCsToTurbSolutionVector(*U,distLSS);
  // Initialize fluid Ids (not on restart)
  if (ioData.input.fluidId[0] == 0 && ioData.input.restart_file_package[0] == 0)
    nodeTag0 = nodeTag = distLSS->getStatus();
  // Initialize the embedded FSI handler
  EmbeddedMeshMotionHandler* _emmh = dynamic_cast<EmbeddedMeshMotionHandler*>(this->emmh);
	if(_emmh) 
	{
    double *tMax = &(this->data)->maxTime;
    _emmh->setup(tMax); //obtain maxTime from structure
  }

  EmbeddedALEMeshMotionHandler* _mmh = dynamic_cast<EmbeddedALEMeshMotionHandler*>(this->mmh);

	if(_mmh) _mmh->setup(*(this->X),this->bcData->getVelocityVector());

	if(this->hth) this->hth->setup(&(this->restart)->frequency, &(this->data)->maxTime);
  
  *(this->Xs) = *(this->X);

  this->initializeFarfieldCoeffs();

  // If 'IncreasePressure' is activated, re-initialize the fluid state
	if(Pinit>=0.0 && Prate>=0.0 && this->getInitialTime()<tmax) 
	{
    increasingPressure = true;
    this->domain->IncreasePressure(currentPressure(this->getInitialTime()), this->varFcn, *U, nodeTag);
  }
 
  //compute force
  DistSVec<double,dim> *Wij = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  DistSVec<double,dim> *Wji = new DistSVec<double,dim>(this->domain->getEdgeDistInfo());
  DistSVec<double,dim> VV(this->getVecInfo());
  *Wij = 0.0;
  *Wji = 0.0;

   this->varFcn->conservativeToPrimitive(*U, VV, distLSS, &nodeTag);

  SubDomain **subD = this->domain->getSubDomain();

  int iSub;
#pragma omp parallel for
	for(iSub=0; iSub<this->domain->getNumLocSub(); iSub++) 
	{
    SVec<double,dim> &subWij = (*Wij)(iSub);
    SVec<double,dim> &subWji = (*Wji)(iSub);
    SVec<double,dim> &subVV = VV(iSub);
    int (*ptr)[2] =  (subD[iSub]->getEdges()).getPtr();

		if(subWij.size()!=(subD[iSub]->getEdges()).size()) 
		{
			fprintf(stderr,"WRONG!!!\n"); 
			exit(-1);
		}

		for(int l=0; l<subWij.size(); l++) 
		{
      int i = ptr[l][0];
      int j = ptr[l][1];

			for(int k=0; k<dim; k++) 
			{
        subWij[l][k] = subVV[i][k];
        subWji[l][k] = subVV[j][k];
      }
    }
  }

	this->spaceOp->computeNodalGrad(*this->X, *this->A, *U, &nodeTag, this->distLSS);

  // Ghost-Points Population
  if(this->eqsType == EmbeddedTsDesc<dim>::NAVIER_STOKES)
      this->spaceOp->populateGhostPoints(this->ghostPoints,*this->X,*U,this->varFcn,this->distLSS,this->viscSecOrder,this->nodeTag);

  // Population of spaceOp->V for the force computation
	this->spaceOp->conservativeToPrimitive(*U, this->distLSS, &this->nodeTag);

  computeForceLoad(Wij, Wji);
  delete Wij;
  delete Wji;

	FsComputed = true;

  // Now "accumulate" the force for the embedded structure
	if(dynNodalTransfer)
	{
    numStructNodes = dynNodalTransfer->numStNodes();
    SVec<double,3> v(numStructNodes, Fs);
    dynNodalTransfer->updateOutputToStructure(0.0, 0.0, v); //dt=dtLeft=0.0-->They are not used!
  }

	if(this->modifiedGhidaglia) this->timeState->attachHH(*this->bcData->getBoundaryStateHH());

}

//------------------------------------------------------------------------------

template<int dim>
double EmbeddedTsDesc<dim>::computeTimeStep(int it, double *dtLeft,
                                             DistSVec<double,dim> &U, double angle)
{
  if(!FsComputed&&dynNodalTransfer) this->com->fprintf(stderr,"WARNING: FSI force not computed!\n");

   //reset FsComputed at the beginning of a fluid iteration
	FsComputed = false; 

  //check if it's in subcycling with iCycle>1.
  if(globIt==it)
    inSubCycling = true;
	else 
	{
    globIt = it;
    inSubCycling = false;
  }

	this->com->barrier();

  double t0 = this->timer->getTime();
  this->data->allowstop = this->timeState->allowcflstop; 
  int numSubCycles = 1;
  double dt=0.0;

	if(TsDesc<dim>::timeStepCalculation == TsData::CFL || it==1)
	{
      this->data->computeCflNumber(it - 1, this->data->residual / this->restart->residual, angle);
		
      if(numFluid==1) 
		{
        dt = this->timeState->computeTimeStep(this->data->cfl, this->data->dualtimecfl, dtLeft,
															  &numSubCycles, *this->geoState, *this->X, *this->A, U, this->distLSS);
      } 
		else 
		{
         //numFLuid>1
        dt = this->timeState->computeTimeStep(this->data->cfl, this->data->dualtimecfl, dtLeft,
                                &numSubCycles, *this->geoState, *this->A, U, nodeTag);
      }
    }
	else 
	{ 
      //time step size with error estimation
      dt = this->timeState->computeTimeStep(it, dtLeft, &numSubCycles);
   }
 
  if(TsDesc<dim>::timeStepCalculation == TsData::ERRORESTIMATION && it == 1)
    this->timeState->setDtMin(dt * TsDesc<dim>::data->getCflMinOverCfl0());

  if (this->problemType[ProblemData::UNSTEADY])
    this->com->printf(5, "Global dt: %g (remaining subcycles = %d)\n",
                      dt*this->refVal->time, numSubCycles);
 
  this->timer->addFluidSolutionTime(t0);
  this->timer->addTimeStepTime(t0);

  dtf = dt;
  dtfLeft = *dtLeft + dt;

  return dt;
}
//------------------------------------------------------------------------------

template<int dim>
double EmbeddedTsDesc<dim>::computePositionVector(bool *lastIt, int it, double t, DistSVec<double,dim> &U)
{
  double dt = 0.0;

  if (this->emmh && this->emmh->structureSubcycling()) {
    double dtleft = 0.0;
    double dtf = this->computeTimeStep(it+1, &dtleft, U); 
    this->emmh->storeFluidSuggestedTimestep(dtf);
  }

  if (this->emmh) {
    double t0 = this->timer->getTime();
    dt = this->emmh->updateStep1(lastIt, it, t, this->bcData->getVelocityVector(), *(this->Xs), &(this->data)->maxTime);
    this->timer->addMeshSolutionTime(t0);
  }
    
  if (this->emmh) {
    double t0 = this->timer->getTime();
    this->emmh->updateStep2(lastIt, it, t, this->bcData->getVelocityVector(), *(this->Xs));
    this->timer->addMeshSolutionTime(t0);
  }

  if (this->mmh) {
    double t0 = this->timer->getTime();
    this->mmh->updateStep1(lastIt, it, t, this->bcData->getVelocityVector(), *(this->Xs), &(this->data)->maxTime);
    this->timer->addMeshSolutionTime(t0);
  }

  if (this->hth) {
    double dth = this->hth->updateStep1(lastIt, it, this->bcData->getTemperatureVector());
    if (!this->mmh && !this->emmh)
      dt = dth;
  }

  if (this->mmh) {
    double t0 = this->timer->getTime();
    this->mmh->updateStep2(lastIt, it, t, this->bcData->getVelocityVector(), *(this->Xs));
    this->timer->addMeshSolutionTime(t0);
  }

  if (this->hth) {
    this->hth->updateStep2(lastIt, it, this->bcData->getTemperatureVector());
  }


  // Once we know the time step and the time, if we are doing an exact solution problem,
  // set up Unm1.  This is not the best place to do this, but oh well.
  
  // In the case of an exact solution
  if (ioData.embed.testCase == 1 && it == 0) {

#pragma omp parallel for
    for (int iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {

      int lsize = U(iSub).size();
      for (int i = 0; i < lsize; ++i) {

	double* x = (*this->X)(iSub)[i];
	double V[5];
	ExactSolution::AcousticBeam(ioData,x[0],x[1],x[2],-dt, V);

	this->varFcn->primitiveToConservative(V, this->timeState->getUnm1()(iSub)[i], 0); //*************************
	
      }
    }
    this->timeState->setExistsNm1();
    this->timeState->setDtNm1(dt);
  } else if (ioData.embed.testCase == 2 && it == 0) {

#pragma omp parallel for
    for (int iSub=0; iSub<this->domain->getNumLocSub(); iSub++) {

      int lsize = U(iSub).size();
      for (int i = 0; i < lsize; ++i) {

	double* x = (*this->X)(iSub)[i];
	double V[5];
	ExactSolution::AcousticViscousBeam(ioData,x[0],x[1],x[2],-dt, V);

			this->varFcn->primitiveToConservative(V, this->timeState->getUnm1()(iSub)[i], 0); //*************************
	
      }
    }
    this->timeState->setExistsNm1();
    this->timeState->setDtNm1(dt);
  }



  return dt;

}

//---------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::updateStateVectors(DistSVec<double,dim> &U, int it)
{
  this->geoState->update(*this->X, *this->A);
  this->timeState->update(U,increasingPressure); 
  this->spaceOp->updateFixes();

}

//-----------------------------------------------------------------------------

template<int dim>
int EmbeddedTsDesc<dim>::checkSolution(DistSVec<double,dim> &U)
{
  int ierr = 0;
  if(numFluid==1) {
    if (dim == 6)
      ierr = this->domain->template 
        clipSolution<dim,1>(this->clippingType, this->wallType, this->varFcn, this->bcData->getInletConservativeState(), U);
    else if (dim == 7)
      ierr = this->domain->template 
        clipSolution<dim,2>(this->clippingType, this->wallType, this->varFcn, this->bcData->getInletConservativeState(), U);
    else
		 ierr = this->domain->checkSolution(this->varFcn, U, this->distLSS); //also check ghost nodes.
  }
  else {
	  ierr = this->domain->checkSolution(this->varFcn, U, nodeTag, this->distLSS);
  }

  return ierr;
}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::setupOutputToDisk(IoData &ioData, bool *lastIt, int it, double t,
                                                  DistSVec<double,dim> &U)
{

  if (it == this->data->maxIts)
    *lastIt = true;
  else 
    monitorInitialState(it, U);

  this->output->setMeshMotionHandler(ioData, this->mmh);
  this->output->openAsciiFiles();
  this->timer->setSetupTime();
  this->output->cleanProbesFile();
  double fluxNorm = 0.5*(this->data->residual)*(this->data->residual);

  if (it == 0) {
    // First time step: compute GradP before computing forces
	 this->spaceOp->computeGradP(*this->X, *this->A, U, &nodeTag, this->distLSS);  
    this->output->writeForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &nodeTag);
    this->output->writeLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &nodeTag);
    this->output->writeMatchPressureToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, *this->A, U, this->timeState, &nodeTag);
    this->output->writeFluxNormToDisk(it, 0, 0, t, fluxNorm);
    this->output->writeHydroForcesToDisk(*lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &nodeTag);
    this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, this->restart->energy, *this->X, U, &nodeTag);
    this->output->writeResidualsToDisk(it, 0.0, 1.0, this->data->cfl);
    // Lei Lei, 02/01/2015: added for embedded ROM training purpose
    DistSVec<char, dim> temp(this->domain->getNodeDistInfo());
    Bool2Char(distLSS->getIsActive(), temp);
    this->output->writeStateMaskVectorsToDiskRom(it, U, temp);
    //
    this->output->writeMaterialVolumesToDisk(it, 0.0, *this->A, &nodeTag);
    this->output->writeMaterialMassEnergyToDisk(it, 0.0, U, *this->A, &nodeTag);
    this->output->writeCPUTimingToDisk(*lastIt, it, t, this->timer);
    this->output->writeEmbeddedSurfaceToDisk(*lastIt, it, t, distLSS->getStructPosition_n(), distLSS->getStructPosition_0());
    this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, nodeTag, this->Wextij, this->distLSS, ghostPoints);
    this->output->writeAvgVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::outputToDisk(IoData &ioData, bool* lastIt, int it, int itSc, int itNl,
                                             double t, double dt, DistSVec<double,dim> &U)
{

  this->com->globalSum(1, &interruptCode);
  if (interruptCode)
    *lastIt = true;

  double cpu = this->timer->getRunTime();
  double res = this->data->residual / this->restart->residual;
  double fluxNorm = 0.5*(this->data->residual)*(this->data->residual);

  this->output->writeLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &nodeTag);
  this->output->writeMatchPressureToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, *this->A, U, this->timeState, &nodeTag);
  this->output->writeFluxNormToDisk(it, itSc, itNl, t, fluxNorm);
  this->output->writeHydroForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &nodeTag);
  this->output->writeHydroLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &nodeTag);
  this->output->writeResidualsToDisk(it, cpu, res, this->data->cfl);
  // Lei Lei, 02/01/2015: added for embedded ROM training purpose
  /*
  DistSVec<char, dim> temp(this->domain->getNodeDistInfo());
  Bool2Char(distLSS->getIsActive(), temp);
  this->output->writeStateMaskVectorsToDiskRom(it, U, temp);
   */
  //
  this->output->writeMaterialVolumesToDisk(it, t, *this->A, &nodeTag);
  this->output->writeMaterialMassEnergyToDisk(it, t, U, *this->A, &nodeTag);
  this->output->writeCPUTimingToDisk(*lastIt, it, t, this->timer);
  if (*lastIt)
    this->output->writeEmbeddedSurfaceToDisk(*lastIt, it, t, distLSS->getStructPosition(), distLSS->getStructPosition_0());
  else
    this->output->writeEmbeddedSurfaceToDisk(*lastIt, it, t, distLSS->getStructPosition_n(), distLSS->getStructPosition_0());
  this->output->writeBinaryVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState, nodeTag, this->Wextij, this->distLSS, ghostPoints); //d2d
  this->output->writeProbesToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState,nodeTag, this->distLSS, ghostPoints);
  this->output->writeAvgVectorsToDisk(*lastIt, it, t, *this->X, *this->A, U, this->timeState);

  TsRestart *restart2 = this->restart; // Bug: compiler does not accept this->restart->writeToDisk<dim,1>(...)
                                       //      it does not seem to understand the template

  // temporary fluid selector
  FluidSelector fluidSelector(nodeTag,this->domain);
  restart2->template writeToDisk<dim,1>(this->com->cpuNum(), *lastIt, it, t, dt, *this->timeState, *this->geoState, NULL, NULL, &fluidSelector);
  if (*lastIt)
    this->restart->writeStructPosToDisk(this->com->cpuNum(), *lastIt, this->distLSS->getStructPosition()); //KW: must be after writeToDisk
  else
    this->restart->writeStructPosToDisk(this->com->cpuNum(), *lastIt, this->distLSS->getStructPosition_n()); //KW: must be after writeToDisk

  this->output->updatePrtout(t);
  this->restart->updatePrtout(t);
  if (*lastIt) {

    this->output->template writeLinePlotsToDisk<1>(true, it, t, *this->X,
						   *this->A, U, 
						   this->timeState, nodeTag);
    this->timer->setRunTime();
    if (this->com->getMaxVerbose() >= 2)
      this->timer->print(this->domain->getStrTimer());
    this->output->closeAsciiFiles();

    if (strcmp(ioData.input.convergence_file,"") != 0) {
      
      this->varFcn->conservativeToPrimitive(U, *this->V, &nodeTag); //*************************
      computeConvergenceInformation(ioData,ioData.input.convergence_file,*this->V);
    }

    if (ioData.embed.testCase == 1 || ioData.embed.testCase == 2) {

      DistSVec<double,dim> Uexact(U);
      
      if (ioData.embed.testCase == 1) {

	ExactSolution::Fill<&ExactSolution::AcousticBeam, dim>(Uexact, *this->X,
							     ioData, t, 
							     this->spaceOp->getVarFcn());

      } else if (ioData.embed.testCase == 2) {

	ExactSolution::Fill<&ExactSolution::AcousticViscousBeam, dim>(Uexact, *this->X,
							     ioData, t, 
							     this->spaceOp->getVarFcn());

      }


      double error[dim];
      double refs[dim] = {ioData.ref.rv.density, ioData.ref.rv.velocity,
			  ioData.ref.rv.velocity, ioData.ref.rv.velocity,
			  ioData.ref.rv.pressure};
      double tot_error = 0.0;
      this->domain->computeL1Error(U,Uexact,*this->A,error, this->distLSS);
      for (int k = 0; k < dim; ++k) {
	tot_error += error[k];
	this->domain->getCommunicator()->fprintf(stdout,"L1 error [%d]: %e\n", k, error[k]*refs[k]);
      }
      this->domain->getCommunicator()->fprintf(stdout,"L1 error (total): %e\n", tot_error);
 
      tot_error = 0.0;
      this->domain->computeL2Error(U,Uexact,*this->A,error, this->distLSS);
      for (int k = 0; k < dim; ++k) {
	tot_error += error[k];
	this->domain->getCommunicator()->fprintf(stdout,"L2 error [%d]: %e\n", k, error[k]*refs[k]);
      }
      this->domain->getCommunicator()->fprintf(stdout,"L2 error (total): %e\n", tot_error);
      
      tot_error = 0.0;
      this->domain->computeLInfError(U,Uexact,error, this->distLSS);
      for (int k = 0; k < dim; ++k) {
	tot_error = max(error[k],tot_error);
	this->domain->getCommunicator()->fprintf(stdout,"Linf error [%d]: %e\n", k, error[k]*refs[k]);
      }
      this->domain->getCommunicator()->fprintf(stdout,"Linf error (total): %e\n", tot_error);

  
    }

  }

}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::outputForces(IoData &ioData, bool* lastIt, int it, int itSc, int itNl,
                               double t, double dt, DistSVec<double,dim> &U)
{ 

  double cpu = this->timer->getRunTime();
  if(this->numFluid==1)
    this->output->writeForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U);
  else
    this->output->writeForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, this->restart->energy, *this->X, U, &this->nodeTag);
}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::outputPositionVectorToDisk(DistSVec<double,dim> &U) 
{
  TsDesc<dim>::outputPositionVectorToDisk(U);
  
  if(emmh && emmh->getAlgNum() == 1)
    this->restart->writeStructPosToDisk(this->com->cpuNum(), true, this->distLSS->getStructPosition());
}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::resetOutputToStructure(DistSVec<double,dim> &U)
{}

//------------------------------------------------------------------------------

template<int dim>
double EmbeddedTsDesc<dim>::computeResidualNorm(DistSVec<double,dim>& U)
{  

  this->spaceOp->computeResidual(*this->X, *this->A, U, *Wstarij, *Wstarji, *Wextij, distLSS, 
											linRecAtInterface,  viscSecOrder, nodeTag, *this->R, 
											this->riemann, riemannNormal, 0, ghostPoints);

  this->spaceOp->applyBCsToResidual(U, *this->R, distLSS);

  double res = 0.0;
  if(this->numFluid==1)
    res = this->spaceOp->computeRealFluidResidual(*this->R, *this->Rreal, *distLSS);
  else {
    this->com->fprintf(stderr,"WARNING: EmbeddedTsDesc::computeResidualNorm is not implemented for numFluid > 1\n");
    res = 1.0;
  }

  return sqrt(res);
}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::monitorInitialState(int it, DistSVec<double,dim> &U)
{

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

template<int dim>
void EmbeddedTsDesc<dim>::computeForceLoad(DistSVec<double,dim> *Wij, DistSVec<double,dim> *Wji)
{
	
	if(!Fs)
	{
		fprintf(stderr,"computeForceLoad: Fs not initialized! Cannot compute the load!\n"); 
		return;
	}

  double t0 = this->timer->getTime();

	if(dynNodalTransfer) numStructNodes = dynNodalTransfer->numStNodes();
  
	if(!increasingPressure || recomputeIntersections) 
	{
		for(int i=0; i<numStructNodes; i++)  Fs[i][0] = Fs[i][1] = Fs[i][2] = 0.0;

		this->spaceOp->computeForceLoad(forceApp, orderOfAccuracy, *this->X,*this->A, Fs, 
												  numStructNodes, distLSS, *Wij, *Wji,  this->Wextij,
                                    ghostPoints, this->postOp->getPostFcn(), &nodeTag);
	} 
	else 
	{
		if(unifPressure[0]==0) 
		{
      this->com->fprintf(stderr,"ERROR: Detected pressure p = %e in Implosion Setup.\n", unifPressure[0]);
      exit(-1);
    }
		for(int i=0; i<numStructNodes; i++) 
		{
      Fs[i][0] *= unifPressure[1]/unifPressure[0]; 
      Fs[i][1] *= unifPressure[1]/unifPressure[0]; 
      Fs[i][2] *= unifPressure[1]/unifPressure[0]; 
    }
  }

  this->timer->addEmbeddedForceTime(t0);
  //at this stage Fs is NOT globally assembled!
}

//-------------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::computederivativeOfForceLoad(DistSVec<double,dim> *Wij, 
						       DistSVec<double,dim> *Wji,
						       double dS[3],
						       DistSVec<double,dim> &dV) {

  if (!dFs) {fprintf(stderr,"dFs not initialized!"); exit(-1);}

  double t0 = this->timer->getTime();
 
  for (int i=0; i<numStructNodes; i++) 
    dFs[i][0] = dFs[i][1] = dFs[i][2] = 0.0;

  this->spaceOp->computederivativeOfForceLoad(forceApp, orderOfAccuracy, *this->X, *this->A,
					      dFs, numStructNodes, distLSS, *Wij, *Wji, dV, dS,
					      ghostPoints, this->postOp->getPostFcn(), &nodeTag);

  this->timer->addEmbeddedForceTime(t0);

}

//-------------------------------------------------------------------------------

template <int dim>
void EmbeddedTsDesc<dim>::getForcesAndMoments(map<int,int> & surfOutMap, DistSVec<double,dim> &U, DistSVec<double,3> &X,
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

template <int dim>
void EmbeddedTsDesc<dim>::getderivativeOfForcesAndMoments(map<int,int> & surfOutMap, 
							  DistSVec<double,dim> &V, DistSVec<double,dim> &dV, 
							  DistSVec<double,3> &X, double dS[3],
							  Vec3D *dFi, Vec3D *dMi) 
{

  int idx;
  computederivativeOfForceLoad(this->Wstarij, this->Wstarji, dS, dV);

  Vec<Vec3D>& Xstruc = distLSS->getStructPosition();

  for (int i=0; i<numStructNodes; i++) {

     map<int,int>::iterator it = surfOutMap.find(distLSS->getSurfaceID(i));

     if(it != surfOutMap.end() && it->second != -2)
       idx = it->second;
     else {
       idx = 0;
     }

     dFi[idx][0] += dFs[i][0]; 
     dFi[idx][1] += dFs[i][1]; 
     dFi[idx][2] += dFs[i][2];

     Vec<Vec3D>& dXstruc = distLSS->getStructDerivative();
     
     dMi[idx][0] += (dXstruc[i][1]* Fs[i][2] -  dXstruc[i][2]* Fs[i][1] +
                      Xstruc[i][2]*dFs[i][2] -   Xstruc[i][2]*dFs[i][1]);

     dMi[idx][1] += (dXstruc[i][2]* Fs[i][0] - dXstruc[i][0]* Fs[i][2] +
       	              Xstruc[i][2]*dFs[i][0] -  Xstruc[i][0]*dFs[i][2]);

     dMi[idx][2] += (dXstruc[i][0]* Fs[i][1] - dXstruc[i][1]* Fs[i][0] +
		      Xstruc[i][0]*dFs[i][1] -  Xstruc[i][1]*dFs[i][0]);
  
  }

}

//-------------------------------------------------------------------------------

template <int dim>
void EmbeddedTsDesc<dim>::updateOutputToStructure(double dt, double dtLeft, DistSVec<double,dim> &U)
{
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

template<int dim>
bool EmbeddedTsDesc<dim>::IncreasePressure(int it, double dt, double t, DistSVec<double,dim> &U)
{

  increasingPressure = false;
  if(Pinit<0.0 || Prate<0.0) return true; // no setup for increasing pressure

  if(t>tmax && t-dt>tmax) {// max pressure was reached, so now we solve
    return true;
  } 

  increasingPressure = true;

  // max pressure not reached, so we do not solve and we increase pressure and let structure react
  
  if(this->emmh && !inSubCycling) {
    //get structure timestep dts
    this->dts = this->emmh->update(0, 0, 0, this->bcData->getVelocityVector(), *this->Xs);
    //recompute intersections
    double tw = this->timer->getTime();
    if(intersector_freq==1||((it-1)%intersector_freq==0&&(it>1))) {
      this->com->fprintf(stderr,"recomputing fluid-structure intersections.\n");
      recomputeIntersections = true;
      this->distLSS->recompute(this->dtf, this->dtfLeft, this->dts, true, TsDesc<dim>::failSafeFlag); 
    } else
      recomputeIntersections = false;

    this->timer->addIntersectionTime(tw);
    this->timer->removeIntersAndPhaseChange(tw);

    nodeTag0 = this->nodeTag;
    nodeTag = this->distLSS->getStatus();

    //store previous states for phase-change update
    tw = this->timer->getTime();
    if(recomputeIntersections)
      this->spaceOp->updateSweptNodes(*this->X,*this->A, this->phaseChangeChoice, this->phaseChangeAlg, U, this->Vtemp,
            *this->Weights, *this->VWeights, *this->Wstarij, *this->Wstarji,
				      this->distLSS, (double*)this->vfar, this->ioData.embed.interfaceLimiter == EmbeddedFramework::LIMITERALEX1,(this->numFluid == 1 ? (DistVec<int>*)0 : &this->nodeTag));
    this->timer->addEmbedPhaseChangeTime(tw);
    this->timer->removeIntersAndPhaseChange(tw);
  } 
  
  // construct Wij, Wji from U. 
  DistSVec<double,dim> VV(this->getVecInfo());
  this->varFcn->conservativeToPrimitive(U,VV,&nodeTag); //*************************
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
  this->spaceOp->conservativeToPrimitive(U, &this->nodeTag); // PJSA  //*************************

  double pnow = currentPressure(t);
  this->com->fprintf(stdout, "about to increase pressure to %e\n", pnow*Pscale);
  this->domain->IncreasePressure(pnow, this->varFcn, U, nodeTag);
  unifPressure[0] = unifPressure[1];
  unifPressure[1] = pnow;

  return false;

}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::fixSolution(DistSVec<double,dim>& U,DistSVec<double,dim>& dU) {

  if (this->fixSol == 1)
    this->domain->fixSolution(this->varFcn,U,dU,&this->nodeTag);
}

//------------------------------------------------------------------------------

template<int dim>
double EmbeddedTsDesc<dim>::currentPressure(double t)
{
  double p;
  if(implosionSetupType==EmbeddedTsDesc<dim>::LINEAR)
    p = Pinit + t*Prate;
  else if(implosionSetupType==EmbeddedTsDesc<dim>::SMOOTHSTEP) {
    double tbar = t/tmax;
    //fprintf(stderr,"t = %e, tmax = %e, tbar = %e.\n", t, tmax, tbar);
    p = Pinit + (Pfinal - Pinit)*tbar*tbar*tbar*(10.0-15.0*tbar+6.0*tbar*tbar);
  } else {
    this->com->fprintf(stderr,"ERROR! ImplosionSetup::Type = %d NOT recognized!\n", implosionSetupType);
    exit(-1);
  }
  return p;
}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::computeDistanceToWall(IoData &ioData)
{ 
	if(ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY)
	{
    if (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
         ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES) 
		{
      double t0 = this->timer->getTime();

      wall_computer->ComputeWallFunction(*this->distLSS,*this->X,*this->geoState);

      this->timer->addWallDistanceTime(t0);

      this->domain->computeOffWallNode(this->distLSS);
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
MeshMotionHandler *EmbeddedTsDesc<dim>::
createEmbeddedALEMeshMotionHandler(IoData &ioData, GeoSource &geoSource, DistLevelSetStructure *distLSS)
{

  MeshMotionHandler *_mmh = 0;

  if (ioData.problem.type[ProblemData::AERO]) {
    _mmh = new EmbeddedALEMeshMotionHandler(ioData, this->domain, geoSource.getMatchNodes(), distLSS);
    //check that algorithm number is consistent with simulation in special case RK2-CD
    // if C0 and RK2 then RK2DGCL is needed!
    if(_mmh->getAlgNum() == 20 || _mmh->getAlgNum() == 21 || _mmh->getAlgNum() == 22){
      if(ioData.ts.type == TsData::EXPLICIT &&
         (ioData.ts.expl.type == ExplicitData::RUNGE_KUTTA_2 ||
          ioData.ts.expl.type == ExplicitData::ONE_BLOCK_RK2 ||
          ioData.ts.expl.type == ExplicitData::ONE_BLOCK_RK2bis )){
        if(!((ioData.dgcl.normals    == DGCLData::EXPLICIT_RK2     || ioData.dgcl.normals == DGCLData::AUTO) &&
             (ioData.dgcl.velocities == DGCLData::EXPLICIT_RK2_VEL || ioData.dgcl.velocities == DGCLData::AUTO_VEL))){
          this->com->fprintf(stderr, "***Error: Computation of the normals or velocities (%d,%d)\n", ioData.dgcl.normals, ioData.dgcl.velocities);
          this->com->fprintf(stderr, "***       is not consistent with Aeroelastic algorithm\n");
          exit(1);
        }
      }
    }
  }
  else if (ioData.problem.type[ProblemData::FORCED]) {
    _mmh = new EmbeddedALEMeshMotionHandler(ioData, this->domain, geoSource.getMatchNodes(), distLSS);
  }
  else if (ioData.problem.type[ProblemData::ACCELERATED])
    _mmh = new AccMeshMotionHandler(ioData, this->varFcn, this->bcData->getInletPrimitiveState(), this->domain);
  else if (ioData.problem.type[ProblemData::ROLL])
    _mmh = new RigidRollMeshMotionHandler(ioData, this->bcData->getInletAngles(), this->domain);
  else if (ioData.problem.type[ProblemData::RBM])
    _mmh = new RbmExtractor(ioData, this->domain);

  return _mmh;

}


template<int dim>
void EmbeddedTsDesc<dim>::computeConvergenceInformation(IoData &ioData, const char* file, DistSVec<double,dim>& U) {

  DistSVec<double,dim> Uexact(U);
  FluidSelector F2(2, ioData, this->domain);
  OneDimensional::read1DSolution(ioData,file, Uexact, 
				 (DistSVec<double,1>*)0,
				 &F2,
				 this->spaceOp->getVarFcn(),
				 *this->X,
				 *this->domain,
				 OneDimensional::ModeU,
				 true) ;

  double error[dim];
  double refs[dim] = {ioData.ref.rv.density, ioData.ref.rv.velocity,
		       ioData.ref.rv.velocity, ioData.ref.rv.velocity,
		      ioData.ref.rv.pressure};
  double tot_error = 0.0;
  this->domain->computeL1Error(U,Uexact,*this->A,error, this->distLSS);
  for (int k = 0; k < dim; ++k) {
    tot_error += error[k];
    this->domain->getCommunicator()->fprintf(stdout,"L1 error [%d]: %lf\n", k, error[k]*refs[k]);
  }
  this->domain->getCommunicator()->fprintf(stdout,"L1 error (total): %lf\n", tot_error);

  tot_error = 0.0;
  this->domain->computeLInfError(U,Uexact,error, this->distLSS);
  for (int k = 0; k < dim; ++k) {
    tot_error = max(error[k],tot_error);
    this->domain->getCommunicator()->fprintf(stdout,"Linf error [%d]: %lf\n", k, error[k]*refs[k]);
  }
  this->domain->getCommunicator()->fprintf(stdout,"Linf error (total): %lf\n", tot_error);

  
}

//------------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::setCurrentTime(double t,DistSVec<double,dim>& U) { 

  currentTime = t;
}

template<int dim>
void EmbeddedTsDesc<dim>::setCurrentTimeStep(double dt) { 

  currentTimeStep = dt;
}

//-------------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::Bool2Char(DistVec<bool> &X, DistSVec<char, dim> &Y) {
  for(int iSub = 0; iSub < X.numLocSub(); iSub++){
    bool *buffer = X.subData(iSub);
    char (*res)[dim] = Y.subData(iSub);
    for(int iNode = 0; iNode < X.subSize(iSub); iNode++){
      for(int k = 0; k < dim; k++){
        res[iNode][k] = buffer[iNode] ? 1 : 0;
      }
    }
  }
}

//-------------------------------------------------------------------------------

template<int dim>
void EmbeddedTsDesc<dim>::writeBinaryVectorsToDiskRom(bool lastNewtonIt, int timeStep, int newtonIter,
DistSVec<double, dim> *state, DistSVec<double, dim> *residual)
{
  if (newtonIter > 0) return;
  //Lei Lei: 04 July 2016, does not capture residual during Newton Iterations for now.
  DistSVec<char, dim> temp(this->domain->getNodeDistInfo());
  Bool2Char(distLSS->getIsActive(), temp);
  this->output->writeStateMaskVectorsToDiskRom(timeStep, *state, temp);
}
