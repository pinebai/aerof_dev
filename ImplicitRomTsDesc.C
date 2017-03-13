#include <DistTimeState.h>
#include <GeoSource.h>
#include <SpaceOperator.h>
#include <Domain.h>
#include <MatVecProd.h>
#include <NewtonSolver.h>

#include <TsOutput.h>

#ifdef TYPE_MAT
#define MatScalar TYPE_MAT
#else
#define MatScalar double
#endif

//------------------------------------------------------------------------------

template<int dim>
ImplicitRomTsDesc<dim>::ImplicitRomTsDesc(IoData &_ioData, GeoSource &geoSource, Domain *dom) :
  TsDesc<dim>(_ioData, geoSource, dom), pod(0, dom->getNodeDistInfo()), F(dom->getNodeDistInfo()), AJ(0, dom->getNodeDistInfo()),
  rom(0) {

	ioData = &_ioData;
  tag = 0;

  maxItsNewton = ioData->ts.implicit.newton.maxIts;
  epsNewton = ioData->ts.implicit.newton.eps;  
  epsAbsResNewton = ioData->ts.implicit.newton.epsAbsRes;  
  epsAbsIncNewton = ioData->ts.implicit.newton.epsAbsInc;  
  maxItsLS = ioData->ts.implicit.newton.lineSearch.maxIts;
  rho = ioData->ts.implicit.newton.lineSearch.rho;
  c1 = ioData->ts.implicit.newton.lineSearch.c1;

  this->timeState = new DistTimeState<dim>(*ioData, this->spaceOp, this->varFcn, this->domain, this->V);

  ImplicitData &implicitData = ioData->ts.implicit;
  if (implicitData.mvp == ImplicitData::FD || implicitData.mvp == ImplicitData::H1FD) {
    mvp = new MatVecProdFD<dim,dim>(implicitData, this->timeState, this->geoState, this->spaceOp, this->domain, *ioData);
  } else if (implicitData.mvp == ImplicitData::H1) {
    //mvp = new MatVecProdH1<dim,MatScalar,dim>(this->timeState, this->spaceOp, this->domain, *ioData);
    this->com->fprintf(stderr, "*** Warning: Switching to finite difference matrix vector products (approximate is not currently supported for ROMs)");
    mvp = new MatVecProdFD<dim,dim>(implicitData, this->timeState, this->geoState, this->spaceOp, this->domain, *ioData);
  } else if (implicitData.mvp == ImplicitData::H2) {
    //mvp = new MatVecProdH2<dim,MatScalar,dim>(*ioData, this->varFcn, this->timeState, this->spaceOp, this->domain, this->geoState);
    this->com->fprintf(stderr, "*** Warning: Switching to finite difference matrix vector products (exact is not currently supported for ROMs)");
    mvp = new MatVecProdFD<dim,dim>(implicitData, this->timeState, this->geoState, this->spaceOp, this->domain, *ioData);
  }

  std::vector<double>* weights = NULL;
  if (strcmp(ioData->input.multiSolutionsParams,"")!=0) {
    this->formInterpolationWeights(*ioData);
    weights = &this->interpolatedICWeights; 
  }

  if ((ioData->problem.alltype == ProblemData::_STEADY_NONLINEAR_ROM_POST_
           || ioData->problem.alltype == ProblemData::_UNSTEADY_NONLINEAR_ROM_POST_)
           && (strcmp(ioData->romDatabase.files.surfacePrefix,"")!=0 
               || strcmp(ioData->romDatabase.files.surfaceStateBasisName,"")!=0)) {
    rom = new NonlinearRomOnlineIII<dim>(dom->getCommunicator(), _ioData, *dom, weights);
    systemApprox=true;
  } else { 
    switch (ioData->romOnline.systemApproximation) {
      case (NonlinearRomOnlineData::SYSTEM_APPROXIMATION_NONE):
        rom = new NonlinearRomOnlineII<dim>(dom->getCommunicator(), _ioData, *dom, weights);
        systemApprox=false;
        break;
      case (NonlinearRomOnlineData::GNAT):
        rom = new NonlinearRomOnlineIII<dim>(dom->getCommunicator(), _ioData, *dom, weights);
        systemApprox=true;
        break;
      case (NonlinearRomOnlineData::COLLOCATION):
        rom = new NonlinearRomOnlineIII<dim>(dom->getCommunicator(), _ioData, *dom, weights);
        systemApprox=true;
        break;
      case (NonlinearRomOnlineData::APPROX_METRIC_NL):
        rom = new NonlinearRomOnlineIII<dim>(dom->getCommunicator(), _ioData, *dom, weights);
        systemApprox=true;
        break;
      default:
        this->com->fprintf(stderr, "*** Error:  Unexpected system approximation type\n");
        exit (-1);
    }
  }

  currentCluster = -1;

  // nPod = 0 ?
  nPod = (ioData->romOnline.maxDimension > 0) ? ioData->romOnline.maxDimension : 0;

  dUromNewtonIt = 0.0;
  dUromTimeIt = 0.0;
  dUromCurrentROB = 0.0;
 
  MemoryPool mp;
  this->mmh = this->createMeshMotionHandler(*ioData, geoSource, &mp);
  mvp->exportMemory(&mp);
 
  basisUpdateFreq = ioData->romOnline.basisUpdateFreq;
  tryAllFreq = ioData->romOnline.tryAllFreq;

  updateFreq = false;
  clusterSwitch = false;
  updatePerformed = false;

  componentwiseScalingVec = NULL; 

  levenbergMarquardtWeight = this->ioData->romOnline.levenbergMarquardtWeight;

  Uinit = NULL;

  useIncrements = (ioData->romDatabase.avgIncrementalStates==NonlinearRomFileSystemData::AVG_INCREMENTAL_STATES_TRUE) ? true : false;

  if ((!systemApprox) && useIncrements) {
    Uprev = new DistSVec<double, dim>(this->domain->getNodeDistInfo());
  } else {
    Uprev = NULL;
  }

  unsteady = this->problemType[ProblemData::UNSTEADY];

  tryingAllClusters = false;

  checkSolutionInNewton = false;

  // reduced coordinate homotopy for Spatial-Only problems
  homotopyStepInitial = ioData->romOnline.romSpatialOnlyInitialHomotomyStep;
  homotopyStepMax = ioData->romOnline.romSpatialOnlyMaxHomotomyStep;
  homotopyStepGrowthRate = max(ioData->romOnline.romSpatialOnlyHomotomyStepExpGrowthRate,1.0);

  if ((ioData->ts.implicit.type == ImplicitData::SPATIAL_ONLY) 
      && (homotopyStepInitial > 0) && (homotopyStepInitial <= homotopyStepMax)) {
    if ((ioData->romOnline.systemApproximation == NonlinearRomOnlineData::SYSTEM_APPROXIMATION_NONE) 
        && (ioData->romOnline.lsSolver!=NonlinearRomOnlineData::NORMAL_EQUATIONS)) {
      this->com->fprintf(stderr, "*** Error:  Spatial-Only ROM homotopy is currently only implemented for the normal equations\n");
      exit (-1); 
    }  
    spatialOnlyWithHomotopy = true;
  } else {
    spatialOnlyWithHomotopy = false;
  }

}

//------------------------------------------------------------------------------

template<int dim>
ImplicitRomTsDesc<dim>::~ImplicitRomTsDesc()
{

  if (tag) delete tag;
  if (componentwiseScalingVec) delete componentwiseScalingVec;
  if (Uinit) delete Uinit;
  if (Uprev) delete Uprev;

}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
KspPrec<neq> *ImplicitRomTsDesc<dim>::createPreconditioner(PcData &pcdata, Domain *dom)
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
//  else if (pcdata.type == PcData::MG)
//    _pc = new MultiGridPrec<Scalar,neq>(dom, *this->geoState);

  return _pc;

}



//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomTsDesc<dim>::tryAllClusters(DistSVec<double, dim> &U, const int totalTimeSteps, int* bestCluster) {


  if (rom->nClusters == 1) {
    *bestCluster=0;
    return;
  }

  this->com->fprintf(stdout, " ... trying to reduce residual with all clusters ...\n");

  double bestResReduction = -1.0;
  
  DistSVec<double, dim> Ubackup(this->domain->getNodeDistInfo());
  Ubackup = U;

  int dummyInt = 0;
  bool dummyTrue = true;
  bool dummyFalse = false;

  tryingAllClusters=true;

  for (int iCluster=0; iCluster<(rom->nClusters); ++iCluster) {

    loadCluster(iCluster, dummyTrue, U);

    computeFullResidual(dummyInt, U, true);
    double resInit = (this->F)*(this->F);

    solveNonLinearSystem(U, totalTimeSteps);

    computeFullResidual(dummyInt, U, true);
    double resFinal = (this->F)*(this->F);

    if (resFinal==0.0) {
      *bestCluster = iCluster;
      return;
    }

    double resReduction = resFinal / resInit;

    this->com->fprintf(stdout, " ... reduced residual from %e to %e (ratio of %e) using cluster %d\n\n",
                              resInit, resFinal, resReduction, iCluster);
  
    if ((resReduction<bestResReduction) || (bestResReduction<0.0)) {
      *bestCluster = iCluster;
      bestResReduction = resReduction;
    }

    U = Ubackup;

  }

  tryingAllClusters=false;
  currentCluster = -1;

}

//------------------------------------------------------------------------------
template<int dim>
void ImplicitRomTsDesc<dim>::checkLocalRomStatus(DistSVec<double, dim> &U, const int totalTimeSteps)  {

  // checks whether the local ROM needs to be modified

  if ((rom->nClusters > 1) || (basisUpdateFreq>0) || (currentCluster == -1)) {

    if (currentCluster == -1) { // first iteration
      if (ioData->romOnline.distanceComparisons)
        rom->initializeDistanceComparisons(U);
      if (ioData->romOnline.basisUpdates==NonlinearRomOnlineData::UPDATES_FAST_EXACT) {
        rom->initializeFastExactUpdatesQuantities(U);
      } else if (ioData->romOnline.projectSwitchStateOntoAffineSubspace!=NonlinearRomOnlineData::PROJECT_OFF) {
        rom->initializeProjectionQuantities(U);
      }
    }

    int closestCluster;

    if (rom->nClusters > 1) {
      if (useIncrements) {
        bool tryAllNow = ((tryAllFreq > 0) && (totalTimeSteps%tryAllFreq == 0)) ? true : false;
        if ((currentCluster==-1) || (tryAllNow) || (dUromTimeIt.norm()<ioData->romOnline.incrementCoordsTol)) { 
          tryAllClusters(U, totalTimeSteps, &closestCluster);
        } else {
          if (systemApprox) {
            rom->closestCenter(U, &closestCluster); // U isn't actually used (and Uprev is NULL)
          } else {
            DistSVec<double, dim> Uincrement(this->domain->getNodeDistInfo());
            Uincrement = U - *Uprev; 
            rom->closestCenter(Uincrement, &closestCluster);
          }
        }
      } else {
        rom->closestCenter(U, &closestCluster);
      }
      this->com->fprintf(stdout, " ... using cluster %d\n", closestCluster);
    } else {
      closestCluster = 0;
    }

    updateFreq = ((basisUpdateFreq > 0) && (totalTimeSteps%basisUpdateFreq == 0)) ? true : false;
    clusterSwitch = (currentCluster != closestCluster) ? true : false;
    updatePerformed = false;

    if (updateFreq || clusterSwitch) loadCluster(closestCluster, clusterSwitch, U);

  } else { // single basis, no updateFreq
    updateFreq = false;
    clusterSwitch = false;
    updatePerformed = false;
  }

  dUromTimeIt = 0.0;

}

//-----------------------------------------------------------------------

template<int dim>
void ImplicitRomTsDesc<dim>::loadCluster(int closestCluster, bool clusterSwitch, DistSVec<double, dim> &U) {

  if (clusterSwitch) {
    deleteRestrictedQuantities(); 
    int prevCluster = currentCluster;
    currentCluster = closestCluster;
    rom->readClusteredOnlineQuantities(currentCluster);  // read state basis, update info, and (if applicable) gappy matrices
    if (this->ioData->romOnline.projectSwitchStateOntoAffineSubspace!=NonlinearRomOnlineData::PROJECT_OFF) {
      if (this->ioData->romOnline.basisUpdates==NonlinearRomOnlineData::UPDATES_OFF) {
        rom->projectSwitchStateOntoAffineSubspace(currentCluster, prevCluster, U, UromCurrentROB);
      } else {
        this->com->fprintf(stderr, "*** Warning: Updates and projection were both specified; not performing projection\n");
      }
    }
  }

  if (this->ioData->romOnline.basisUpdates!=NonlinearRomOnlineData::UPDATES_OFF) {
    updatePerformed = rom->updateBasis(currentCluster, U, &dUromCurrentROB);
    if (this->ioData->romOnline.bufferEnergy!=0.0) rom->truncateBufferedBasis();
    if (this->ioData->romOnline.krylov.include) rom->appendNonStateDataToBasis(currentCluster,"krylov");
    if (this->ioData->romOnline.sensitivity.include) rom->appendNonStateDataToBasis(currentCluster,"sensitivity");
  }
  
  nPod = rom->basis->numVectors();
  pod.resize(nPod);
  for (int iVec=0; iVec<nPod; ++iVec) {
    pod[iVec] = (*(rom->basis))[iVec];
  }
  AJ.resize(nPod);
  dUromNewtonIt.resize(nPod);
  dUromTimeIt.resize(nPod);
  dUromCurrentROB.resize(nPod);
  dUromCurrentROB = 0.0;
  setProblemSize(U);  // defined in derived classes
  // TODO also set new reference residual if the weighting changes
  if (clusterSwitch && !unsteady) setReferenceResidual(); // for steady gappy simulations (reference residual is restricted to currently active nodes)

}


//------------------------------------------------------------------------------
template<int dim>
int ImplicitRomTsDesc<dim>::solveNonLinearSystem(DistSVec<double, dim> &U, const int totalTimeSteps)  {

  if (Uprev) *Uprev = U;

  // initializations 

  double t0 = this->timer->getTime();

  int it = 0;
  int fsIt = 0;

  DistSVec<double, dim> dUfull(this->domain->getNodeDistInfo());	// solution increment at EACH NEWTON ITERATION in full coordinates
  dUfull = 0.0;	// initial zero increment

  double res;
  bool breakloop = false;
  bool breakloopNow = false;

  // line search variables
  double alpha;
  bool convergeFlag=0;

  postProStep(U,totalTimeSteps);

  if (totalTimeSteps == 1 && this->ioData->romOnline.residualScaling != NonlinearRomOnlineData::SCALING_OFF) {
    computeFullResidual(it, U, true); 
  }

  if (unsteady && (this->ioData->romOnline.residualScaling == NonlinearRomOnlineData::SCALING_BALANCED)) {  
    this->com->fprintf(stdout, " ... using the initial weighting vector for all timesteps\n");
    // reset the balancing scaling at the start of every timestep for unsteady?
     //this->com->fprintf(stdout, "resetting the weighting vector\n");
     //(this->varFcn->weights)[0]=-1;
  }

  double maxStep = -1;
  double newtonStepThreshold = this->ioData->romOnline.newtonStepThreshold;

  for (it = 0; it < maxItsNewton; it++)  {

    computeFullResidual(it, U, false);

    double tAJ = this->timer->getTime();
    computeAJ(it, U, true);	// skipped some times for Broyden
    this->timer->addAJTime(tAJ);

    if (this->ioData->romOnline.residualScaling != NonlinearRomOnlineData::SCALING_OFF) {
      computeFullResidual(it, U, true);
    }

    solveNewtonSystem(it, res, breakloop, U, totalTimeSteps);	// 1) check if residual small enough, 2) solve 
      // INPUTS: AJ, F
      // OUTPUTS: dUromNewtonIt, res, breakloop

    breakloopNow = breakloop1(breakloop);
    if (breakloopNow) break;

    if (this->ioData->romOnline.lineSearch==NonlinearRomOnlineData::LINE_SEARCH_WOLF) { 
      // do line search (linesearch exits with alpha=0 and convergenceFlag if convergence criteria is satisfied)
      alpha = lineSearch(U,dUromNewtonIt,it,AJ,epsNewton, convergeFlag);
      if (it > 0 && convergeFlag == 1) break;
      dUromNewtonIt *= alpha;
    } else if (this->ioData->romOnline.lineSearch==NonlinearRomOnlineData::LINE_SEARCH_BACKTRACKING) {
      double baselineMerit = this->meritFunction(it, U, dUfull, this->F, 0);
      alpha = 1.0;
      expandVector(dUromNewtonIt, dUfull);
      for (int itLS=0; itLS<maxItsLS; ++itLS) {
        if (itLS>0) alpha *= rho;
        double testMerit = this->meritFunction(it, U, dUfull, this->F, alpha);
        //this->com->fprintf(stderr, "*** baselineMerit=%e, baslineMerit*coeff=%e, testMerit=%e\n", baselineMerit, sqrt(1-2.0*alpha*c1)*baselineMerit, testMerit);
        if (testMerit < sqrt(1-2.0*alpha*c1)*baselineMerit) {// || dQ.norm() <= epsAbsInc)
          dUromNewtonIt *= alpha;
          break;
        }
        if (itLS == maxItsLS-1 && maxItsLS != 1) {
          this->com->printf(1, "*** Warning: Line Search reached %d its ***\n", maxItsLS);
          dUromNewtonIt *= alpha;
        }
      }
    }

    double tSol = this->timer->getTime();    
    dUromTimeIt += dUromNewtonIt; // solution increment in reduced coordinates (initialized to zero in checkLocalRomStatus)
    expandVector(dUromNewtonIt, dUfull); // solution increment in full coordinates
    U += dUfull;
    this->timer->addSolutionIncrementTime(tSol);

    saveNewtonSystemVectors(totalTimeSteps); // only implemeted for PG rom

    // verify that the solution is physical (also calls clipping)
    if (checkSolutionInNewton && this->checkSolution(U)) {
      if (checkFailSafe(U) && fsIt < 5) {
        this->com->fprintf(stderr, "*** Warning: Not yet implemented\n");
        //fprintf(stderr,"*** Warning: Newton solver redoing iteration %d\n", it+1);
        //Q = rhs;
        //--it;
        //++fsIt;
      }
      else{
        this->com->fprintf(stderr, "*** Exiting\n");
        exit(-1);
      }
    }

    double step = dUromNewtonIt.norm();
    this->com->fprintf(stdout, " ... step length = %e\n", step);

    if (step > maxStep) {
      maxStep = step;
    } else if (step < newtonStepThreshold) {
      break;
    }
 
    breakloopNow = breakloop2(breakloop);
    if (breakloopNow) break;

  } // end Newton loop
	
  //savedUnormAccum();
  if (fsIt > 0 && checkFailSafe(U) == 1)
  resetFixesTag();

  if (it == maxItsNewton && maxItsNewton!=1 && maxItsNewton!=0) {
    this->com->fprintf(stderr, "*** Warning: ROM Newton solver reached %d its", maxItsNewton);
    this->com->fprintf(stderr, " (Residual: initial=%.2e, reached=%.2e, target=%.2e)\n", res0, res, target);
  }

  if (it==0 && (ioData->problem.alltype != ProblemData::_STEADY_NONLINEAR_ROM_POST_) 
     && (ioData->problem.alltype != ProblemData::_UNSTEADY_NONLINEAR_ROM_POST_)) {
    this->com->fprintf(stderr, "*** Warning: ROM converged on first iteration");
    this->com->fprintf(stderr, " (Residual: initial=%.2e, reached=%.2e, target=%.2e)\n", res0, res, target);
  }

  this->timer->addFluidSolutionTime(t0);

  this->com->fprintf(stdout, "\n");

  // If trying all clusters, exit before writing reduced coordinates
  if (tryingAllClusters) return (maxItsNewton == 0 || it==0) ? 1 : it;

  // output POD coordinates
  dUromCurrentROB += dUromTimeIt;
  if (UromCurrentROB.size() == dUromTimeIt.size()) UromCurrentROB += dUromTimeIt;

  rom->writeReducedCoords(totalTimeSteps, clusterSwitch, updatePerformed, currentCluster, dUromTimeIt); 

  if (ioData->romOnline.distanceComparisons)
    rom->advanceDistanceComparisons(currentCluster, dUromTimeIt, UromCurrentROB);

  return (maxItsNewton == 0 || it==0) ? 1 : it;

}
                                                                                                           
//------------------------------------------------------------------------------

template<int dim>
int ImplicitRomTsDesc<dim>::solveLinearSystem(int it , Vec<double> &rhs, Vec<double> &sol)
{

  double *x = rhs.data();
  FullM myjac(jac);
  myjac.factor();
  myjac.reSolve(x);
  sol = rhs;  

  return 0;

}

//------------------------------------------------------------------------------

// this function evaluates (Aw),t + F(w,x,v)
template<int dim>
void ImplicitRomTsDesc<dim>::computeFullResidual(int it, DistSVec<double, dim> &Q, bool applyWeighting,  DistSVec<double, dim> *R, bool includeHomotopy)
{
  if (R==NULL) R = &F;  // make R an alias for F

  double tRes = this->timer->getTime();

  this->spaceOp->computeResidual(*this->X, *this->A, Q, *R, this->timeState);

  if (includeHomotopy)
      this->timeState->add_dAW_dt(it, *this->geoState, *this->A, Q, *R);

  this->spaceOp->applyBCsToResidual(Q, *R);  // wall BCs only

  if (applyWeighting && (this->ioData->romOnline.residualScaling != NonlinearRomOnlineData::SCALING_OFF)) {
    if (componentwiseScalingVec) delete componentwiseScalingVec;
    componentwiseScalingVec = this->varFcn->computeScalingVec(*(this->ioData),Q,R);
    *R *= *componentwiseScalingVec;
  }

  this->timer->addResidualTime(tRes);
 
}



//------------------------------------------------------------------------------
template<int dim>
double ImplicitRomTsDesc<dim>::meritFunctionDeriv(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> &p, DistSVec<double, dim> &R, double currentMerit)  { 

  MatVecProdFD<dim,dim>* mvpfd = dynamic_cast<MatVecProdFD<dim,dim>*>(mvp);
  if (mvpfd==NULL) {
    this->com->fprintf(stderr,"*** ERROR: Wolfe line search only implemented for MatVecProd=FD! \n"); 
    exit(-1);
  }

  double eps = mvpfd->computeEpsilon(Q,p);
  DistSVec<double, dim> newR(this->domain->getNodeDistInfo());

  double newMerit = this->meritFunction(it, Q, p, newR, eps);

  double meritDeriv = (newMerit - currentMerit) / eps;

  return meritDeriv;

/*
  double eps = mvpfd->computeEpsilon(Q,p);
  DistSVec<double, dim> newF(this->domain->getNodeDistInfo());

  DistSVec<double, dim> newQ(this->domain->getNodeDistInfo());
  newQ = Q + eps*p;
  computeFullResidual(it,newQ,&newF);

  newF -= F;  // overwrite new Flux with finite difference
  newF *= (1.0/eps);

  double meritDeriv = 0.0;
  meritDeriv = newF*F; // Take inner product
    // merit function = norm of full-order residual
  return meritDeriv;
*/
}

//------------------------------------------------------------------------------
template<int dim>
double ImplicitRomTsDesc<dim>::lineSearch(DistSVec<double, dim> &Q, Vec<double> &prom, int it, VecSet<DistSVec<double, dim> > &leftProj, double eps, bool &convergeFlag){
  // This function (along with "zoom"), finds a steplength satisfying the Strong Wolfe conditions in two steps:
  // 1) Find a valid bracket containing a point that satisfies these conditions
  // 2) Perform interpolation/bisection within the bracket to locate this valid point
  // See p. 61 of Numerical Optimization, 2nd ed by Nocedal and Wright for details

  // Parameters for Wolfe conditions

  double c1Wolfe = 0.00001;  // for Armijo (sufficient decrease) condition
  double c2 = 0.0001; // relative gradient condition.
  //NOTE: We require 0<c1<c2<1 
  //NOTE: the smaller the value of c2, the closer the steplength will be to a local minimizer

  // Parameters for bracketing step (step 1)

  int maxIter = 100;  // max iterations for finding a bracket
  DistSVec<double, dim> Qnew(this->domain->getNodeDistInfo()); // new state vector initialization
  double beta = 1.2;  // amount to increment alpha each time
  double alpha_max = 100000000.0;    // bound on step size alpha  (this was about 50)
  int maxedFlag = 0;  // flag to determine whether or not alpha_max has been reached
  int count = 0;  // counter for number of iterations
  double alpha = 1.0; // initial alpha is just the Newton step
  double merit, meritZero, meritOld;
  double meritDeriv, meritDerivZero,targetMeritDerivZero, meritDerivOld;

  DistSVec<double,dim> dQ(this->domain->getNodeDistInfo());
  DistSVec<double,dim> R(this->domain->getNodeDistInfo());
  Vec<double> From(nPod);
  expandVector(prom, dQ);

  //----- Obtain an initial bracket -----//

  meritZero = this->meritFunction(it, Q, dQ, R, 0.0);
  meritOld = meritZero;
  meritDerivZero = this->meritFunctionDeriv(it, Q, dQ, R, meritZero);

  // Reverse the search direction if it is not a descent direction

  if (meritDerivZero > 0) {
   this->com->fprintf(stderr,"Reversing the search direction: original direction NOT a descent direction! \n");
    dQ = -1.0*dQ;
    prom = -1.0*prom;
    meritDerivZero = -meritDerivZero;
  }
  this->com->fprintf(stderr,"meritDerivZero = %e, meritZero = %e \n",meritDerivZero, meritZero);
 
  targetMeritDerivZero = eps*fabs(meritZero);  //Compare targetMeritDeriv to meritZero based on Taylor series explanation
 
  if (fabs(meritDerivZero) <= targetMeritDerivZero){ // convergence criterion based on stationarity of full-order residual norm
     convergeFlag = 1;
     return 0;
  }

  meritDerivOld = meritDerivZero;
  double alphaOld = 0.0;  //Left endpoint is zero
  
  while (count<maxIter){
    merit = this->meritFunction(it, Q, dQ, R, alpha); // evaluate merit function at current alpha
    if (merit>meritZero+c1Wolfe*alpha*meritDerivZero || (merit>=meritOld && count>0)){
        // alphaLo=alphaOld, alphaHi=alpha
        this->com->fprintf(stderr,"Entering zoom from location 1, alpha = %e \n",alpha);
        alpha = zoom(alphaOld, alpha, meritOld, merit, meritDerivOld, meritDerivZero, meritZero, c1Wolfe, c2, Q, dQ, R, it);
        return alpha;
    }
    Qnew = Q+alpha*dQ;  // Q at current alpha
    meritDeriv = this->meritFunctionDeriv(it, Qnew, dQ, R, merit);

    if (fabs(meritDeriv)<=-c2*meritDerivZero){
       this->com->fprintf(stderr,"Condition one is satisfied: %d. Condition two is satisfied: %d.\n",fabs(meritDeriv)<=-c2*meritDerivZero,merit<meritZero+c1Wolfe*alpha*meritDerivZero);
       this->com->fprintf(stderr,"Returning alpha without zoom, alpha = %e \n",alpha);
       return alpha;  // a valid solution has been found
    }
    if (meritDeriv>=0){
       // alphaLo=alpha, alphaHi=alphaOld
       this->com->fprintf(stderr,"Entering zoom from location 2.  The alphas are: alpha = %e, alphaOld = %e \n",alpha,alphaOld);
       alpha = zoom(alpha, alphaOld, merit, meritOld, meritDeriv, meritDerivZero, meritZero, c1Wolfe, c2, Q, dQ, R, it);
       return alpha;
    }

    //----- update values of alpha, merit function, and derivative -----//

    alphaOld = alpha;
    meritOld = merit;
    meritDerivOld = meritDeriv;

    //----- choose alpha between current alpha and alpha_max -----//

    alpha = beta*alphaOld;   //increment alpha by a constant factor

    //----- correct if alpha exceeds alpha_max -----

    while (alpha>alpha_max){
      beta=(1.0+beta)/2.0;    //decrease beta
      alpha=beta*alphaOld; //try a new value of alpha
      this->com->fprintf(stderr,"Decreasing beta because alpha_max violated \n");
    }
    ++count;
    
  }

  this->com->fprintf(stderr,"Leaving linesearch because max iterations was exceeded. Conditions NOT satisfied! alpha = %e\n",alpha);  
  return alpha;

}
//-----------------------------------------------------------------------------
template<int dim>
double ImplicitRomTsDesc<dim>::zoom(double alphaLo, double alphaHi, double meritLo, double meritHi, double meritDerivLo, double meritDerivZero, double meritZero, double c1Wolfe, double c2, DistSVec<double, dim> Q, DistSVec<double, dim> dQ, DistSVec<double, dim> R, int it){

    // Given a valid bracket, this function finds an alpha within the bracket satisfying the Strong Wolfe conditions
    // See p. 61 of Numerical Optimization, 2nd ed by Nocedal and Wright for details

    // initialize parameters
    int count = 0;
    int maxIter = 40; // max number of zooming steps
    double merit, meritDeriv, alpha;
    DistSVec<double, dim> Qnew(this->domain->getNodeDistInfo());
    double epsilon = 0.1; // Percent reduction needed for quadratic approx before bisection is done
    double tolInterp = 0.05; // Determines if interpolation will be executed-we have had problems with the matrix being ill-conditioned
 
    while (count<maxIter){

      // Perform a quadratic interpolation using slope of the low point IF matrix will be well-conditioned using a heuristic

      if (false) { //(fabs(alphaHi-alphaLo)>tolInterp){
      
            FullM Vdm(3); // Vandermonde matrix
   
            Vdm[0][0] = alphaLo*alphaLo; Vdm[0][1] = alphaLo; Vdm[0][2] = 1.0;
            Vdm[1][0] = 2.0*alphaLo; Vdm[1][1] = 1.0; Vdm[1][2] = 0.0;
            Vdm[2][0] = alphaHi*alphaHi; Vdm[2][1] = alphaHi; Vdm[2][2] = 1.0;
            double *VdmRHS = new double[3];
            VdmRHS[2] = meritHi; VdmRHS[1] = meritDerivLo; VdmRHS[0] = meritLo;
            Vdm.Factor();
            Vdm.ReSolve(VdmRHS); //obtain interpolant coefficients
            alpha = -VdmRHS[1]/(2.0*VdmRHS[0]);//new value from minimizing quadratic model
   
            delete[] VdmRHS;
      }

      // ----- Use bisection if interpolation was bypassed or if there was an insufficient reduction-----//

      if (fabs(alphaHi-alphaLo)<=tolInterp || (alpha-alphaLo)/(alphaHi-alphaLo)<epsilon || (alphaHi-alpha)/(alphaHi-alphaLo)<epsilon)
        alpha=(alphaLo+alphaHi)/2.0; // bisection

      // Evaluate merit function

      merit = this->meritFunction(it, Q, dQ, R, alpha);

      if (merit>meritZero+c1Wolfe*alpha*meritDerivZero || merit>=meritLo){
           alphaHi = alpha;
           meritHi = merit;
      }
      else{
           Qnew = Q+alpha*dQ;
           meritDeriv = this->meritFunctionDeriv(it, Qnew, dQ, R, merit);
           if (fabs(meritDeriv)<=-c2*meritDerivZero){
                this->com->fprintf(stderr,"Condition one is satisfied: %d. Condition two is satisfied: %d.\n",fabs(meritDeriv)<=-c2*meritDerivZero,merit<meritZero+c1Wolfe*alpha*meritDerivZero);
                this->com->fprintf(stderr,"alpha = %e, # zoom iter = %d\n",alpha, count);

                return alpha;  // alpha satisfies Strong Wolfe conditions: exit loop
           }
           if (meritDeriv*(alphaHi-alphaLo)>=0){
                alphaHi = alphaLo;
                meritHi = meritLo;
           }
           alphaLo = alpha;
           meritLo = merit;
           meritDerivLo = meritDeriv; 
      }
      ++count;
    }
    this->com->fprintf(stderr,"Leaving zoom because max iterations was exceeded. Conditions NOT satisfied! Count = %d, alpha = %e\n",count, alpha);
    return alpha;

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomTsDesc<dim>::expandVector(Vec<double> &romV, DistSVec<double, dim> &fullV)  {

  fullV = 0.0;
  for (int iVec = 0; iVec < nPod; iVec++)
    fullV += romV[iVec]*pod[iVec];

}


//------------------------------------------------------------------------------


template<int dim>
void ImplicitRomTsDesc<dim>::projectVector(VecSet<DistSVec<double, dim> > &leftProj, DistSVec<double,dim> &fullV, Vec<double> &romV)  {

	transMatVecProd(leftProj,fullV,projVectorTmp);

  for (int iVec = 0; iVec < pod.numVectors(); iVec++)
    romV[iVec] = projVectorTmp[iVec];

}

//------------------------------------------------------------------------------

template<int dim>
int ImplicitRomTsDesc<dim>::checkFailSafe(DistSVec<double,dim>& U)
{

  this->com->fprintf(stderr, "*** Warning: Checkfailsafe Not yet implemented\n");
  int failSafeNewton = 0;
  return failSafeNewton;

}

//------------------------------------------------------------------------------
template<int dim>
void ImplicitRomTsDesc<dim>::resetFixesTag()
{
  this->spaceOp->resetTag();
}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomTsDesc<dim>::computeRedHessianSums(int it, DistSVec<double, dim> &Q)  {

  //mvpfd->evaluate(it, *this->X, *this->A, Q, F);
  
  //for (int iPod = 0; iPod < nPod; iPod++)
  //  mvpfd->apply(pod[iPod], AJ[iPod]);

}

//------------------------------------------------------------------------------


template<int dim>
void ImplicitRomTsDesc<dim>::computeAJ(int it, DistSVec<double, dim> &Q, bool applyWeighting, DistSVec<double, dim> *R)  {

  if (R==NULL) R = &F;

  bool componentScaling = (applyWeighting && (this->ioData->romOnline.residualScaling != NonlinearRomOnlineData::SCALING_OFF));

  double t0 = this->timer->getTime();
  if (componentScaling) {
    mvp->evaluateWeighted(it, *this->X, *this->A, Q, *R, this->varFcn);
  } else {
    mvp->evaluate(it, *this->X, *this->A, Q, *R);
  }

  this->timer->addJacEvaluateTime(t0);

  if (AJ.numVectors()!=nPod) AJ.resize(nPod);  
  t0 = this->timer->getTime();
  for (int iPod = 0; iPod < nPod; iPod++) {
    if (componentScaling) {
      mvp->applyWeighted(pod[iPod], AJ[iPod], this->varFcn);
    } else {
      mvp->apply(pod[iPod], AJ[iPod]);
    }
  }
  this->timer->addJacApplyTime(t0);


}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomTsDesc<dim>::saveNewtonSystemVectorsAction(const int totalTimeSteps) {
	// only do for modelII online and modelIII postpro

  int freq = ioData->output.rom.resjacfrequency;

  if ((freq >= 1) && (totalTimeSteps%freq == 0)) {
    if (rom->jacActionSnapsFileNameSpecified) {
      DistSVec<double, dim> AJsol(this->domain->getNodeDistInfo()); 
      AJsol = 0.0;
      for (int i=0; i<this->nPod; ++i)
         AJsol += this->AJ[i] * this->dUromNewtonIt[i]; 

      // saving 1) residual and 2) this->AJ * this->dUromNewtonIt (for GappyPOD)
      rom->writeClusteredBinaryVectors(currentCluster, &(this->F), &AJsol);
    } else {
      // just residual
      rom->writeClusteredBinaryVectors(currentCluster, &(this->F)); 
    }
  }
}


//------------------------------------------------------------------------------

//template<int dim>
//void ImplicitRomTsDesc<dim>::savedUnormAccum() {

//	for (int iPod = 0; iPod < nPod; ++iPod) {
//		dUnormAccum[0][iPod] += fabs(dUromNewtonIt[iPod]);	// 1 norm
//		dUnormAccum[1][iPod] += dUromNewtonIt[iPod] * dUromNewtonIt[iPod];	// 2 norm
//	}

//}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitRomTsDesc<dim>::rstVarImplicitRomTsDesc(IoData &ioData)
{

  mvp->rstSpaceOp(ioData, this->varFcn, this->spaceOp, false);

}

//------------------------------------------------------------------------------

template<int dim>
bool ImplicitRomTsDesc<dim>::breakloop1(const bool breakloop) {

	return breakloop;
	
}

//------------------------------------------------------------------------------

template<int dim>
bool ImplicitRomTsDesc<dim>::breakloop2(const bool breakloop) {

	return false;

}

//------------------------------------------------------------------------------

