#include <TsDesc.h>
#include <Domain.h>
#include <cmath>
#include <RefVal.h>
#include <GeoSource.h>
#include <DistBcData.h>
#include <DistTimeState.h>
#include <DistGeoState.h>
#include <SpaceOperator.h>
#include <PostOperator.h>
#include <MeshMotionHandler.h>
#include <HeatTransferHandler.h>
#include <DistVector.h>
#include <MemoryPool.h>
#include <Timer.h>
#include <alloca.h>
#include <DistExactRiemannSolver.h>
//#include <RBFInterpND.h>
//#include <r8lib.h>

extern int interruptCode;

//------------------------------------------------------------------------------

template<int dim>
TsDesc<dim>::TsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) : domain(dom),fluidIdDummy(dom->getNodeDistInfo())
{

  X = new DistSVec<double,3>(getVecInfo());
  A = new DistVec<double>(getVecInfo());
  Xs = new DistSVec<double,3>(getVecInfo());

  // Initialize the values
  *X = 0.0;
  *A = 0.0;
  *Xs = 0.0;

  fluidIdDummy = 0;

  Uic = new DistSVec<double,dim>(getVecInfo());
  V = new DistSVec<double,dim>(getVecInfo());
  F = new DistSVec<double,dim>(getVecInfo());
  R = new DistSVec<double,dim>(getVecInfo());
  Rinlet = new DistSVec<double,dim>(getVecInfo());
  Rreal = new DistSVec<double,dim>(getVecInfo());
  timer = domain->getTimer();
  com = domain->getCommunicator();
  errorHandler = domain->getErrorHandler();

  problemType = ioData.problem.type;
  clippingType = ioData.ts.typeClipping;
  wallType = ioData.bc.wall.integration;
  wallRecType = ioData.bc.wall.reconstruction;
  timeStepCalculation = ioData.ts.timeStepCalculation;

  refVal = new RefVal(ioData.ref.rv);

  varFcn = new VarFcn(ioData);

  input = new TsInput(ioData);
  geoState = new DistGeoState(ioData, domain);

  // restart the geoState (positions of the mesh) At return X contains the last
  // position of the mesh
  //
  // .
  if ((ioData.problem.framework==ProblemData::BODYFITTED || ioData.problem.framework==ProblemData::EMBEDDEDALE) &&
      (ioData.problem.type[ProblemData::AERO] ||
       ioData.problem.type[ProblemData::ACCELERATED] ||
       ioData.problem.type[ProblemData::FORCED] ||
       ioData.problem.type[ProblemData::ROLL] ||
       ioData.problem.type[ProblemData::RBM] ||
       ioData.problem.alltype==ProblemData::_STEADY_)) {
    if (strcmp(input->displacements,"")!=0 && strcmp(input->positions,"")==0) {
      // read initial displacement of full mesh
      geoState->setupInitialDisplacement(input->displacements, X, A);
      mems = 0;
    } else {
      geoState->setup1(input->positions, X, A);
      moveMesh(ioData, geoSource);
    }
  } else {
    if ((input->displacements && strcmp(input->displacements,"")!=0) && (!input->positions || strcmp(input->positions,"")==0)) {
      // read initial displacement of full mesh
      geoState->setupInitialDisplacement(input->displacements, X, A);
      mems = 0;
    } else {
      char temp[1]; temp[0] = '\0';
      geoState->setup1(temp, X, A);
      moveMesh(ioData, geoSource);
    }
  }
  bcData = createBcData(ioData);

  spaceOp = new SpaceOperator<dim>(ioData, varFcn, bcData, geoState, domain, V);

  postOp = new PostOperator<dim>(ioData, varFcn, bcData, geoState, domain, V);

  data = new TsParameters(ioData);
  data->assignErrorHandler(dom->getErrorHandler());
  output = new TsOutput<dim>(ioData, refVal, domain, postOp);
  restart = new TsRestart(ioData, refVal);

  hth = createHeatTransferHandler(ioData, geoSource);

  riemann1 = new DistExactRiemannSolver<dim>(ioData, domain, varFcn);
// Included (MB)
  forceNorm = 0.0;
  failSafeFlag = false;
  if (ioData.sa.avgsIt) {
    forceNorms = new double[ioData.sa.avgsIt];
  }
  else {
    forceNorms = 0;
  }

  iForce = 0;
  iTotal = 0;

  modifiedGhidaglia = (ioData.schemes.bc.type==BoundarySchemeData::MODIFIED_GHIDAGLIA);

  if (ioData.sa.fixsol == 0)
    fixSol = 0;
  else if (ioData.sa.fixsol == 1)
    fixSol = 1;

  timeState = 0;
  mmh = 0;

  isMultigridTsDesc = false;
  outputOnlySpatialResidualBool = ioData.output.rom.outputOnlySpatialResidual==ROMOutputData::OUTPUT_ONLY_SPATIAL_RES_ON;

  interpolatedICWeights.clear();
}

//------------------------------------------------------------------------------

template<int dim>
TsDesc<dim>::~TsDesc()
{
  if (X) delete X;
  if (Xs) delete Xs;
  if (A) delete A;
  if (V) delete V;
  if (R) delete R;
  if (Rinlet) delete Rinlet;
  if (Rreal) delete Rreal;
  if (data) delete data;
  if (input) delete input;
  if (output) delete output;
  if (restart) delete restart;
  if (refVal) delete refVal;
  if (varFcn) delete varFcn;
  if (timeState) delete timeState;
  if (bcData) delete bcData;
  if (geoState) delete geoState;
  if (spaceOp) delete spaceOp;
  if (postOp) delete postOp;
  if (mmh) delete mmh;
  if (mems) delete mems;
  if (hth) delete hth;
  if (forceNorms) delete forceNorms;
  if (riemann1) delete riemann1;
}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::moveMesh(IoData &ioData, GeoSource &geoSource)
{
    if (strcmp(input->wallsurfacedisplac,"") != 0 && strcmp(input->positions,"") == 0) {
      PosVecType dXb(getVecInfo());
      mems = new TetMeshMotionSolver(ioData.dmesh, geoSource.getMatchNodes(), domain, 0);
//      mems = new TetMeshMotionSolver(ioData.dmesh, 0, domain, 0);
      domain->readVectorFromFile(input->wallsurfacedisplac, 0, 0, dXb);
      mems->solve(dXb, *X);
      *Xs = *X;
      if(X->norm() == 0.0)
      {
        this->com->fprintf(stderr, "\n *** ERROR *** No Mesh Perturbation \n\n");
        exit(1);
      }
      com->fprintf(stderr," ... mesh has been moved.\n");

/*      double tag = 0.0;
      for(int i=0; i<1; ++i) {
        com->fprintf(stderr," *** mesh sensitivity has been computed 0.\n");
        bool readOK = domain->readVectorFromFile(this->input->shapederivatives, 0, &tag, *dXdSb0);
        com->fprintf(stderr," *** mesh sensitivity has been computed 1.\n");
        if(readOK) {
          // Checking if dXdSb0 has entries different from zero at the interior of the mesh
//          this->postOp->checkVec(*dXdSb0);
          com->fprintf(stderr," *** mesh sensitivity has been computed 2.\n");

          if (dXdSb0->norm() == 0.0)
          {
            this->com->fprintf(stderr, "\n *** WARNING *** No Mesh Perturbation \n\n");
            if(!ioData.sa.fsiFlag) exit(1);
          }
        } else *dXdSb0 = 0.0;
        *dXdS0 = *this->X;
        com->fprintf(stderr," *** mesh sensitivity has been computed 3.\n");
        mems->solve(*dXdSb0, *dXdS0);
        com->fprintf(stderr," *** mesh sensitivity has been computed 4.\n");
        *dXdS0 -= *this->X;
      }
*/
      char temp[1]; temp[0] = '\0';
      geoState->setup3(temp, X, A);
      if(ioData.output.restart.positions[0] != 0) {
        int sp = strlen(ioData.output.restart.prefix) + 1;
        char *posFile = new char[sp + strlen(ioData.output.restart.positions)];
        sprintf(posFile, "%s%s", ioData.output.restart.prefix, ioData.output.restart.positions);
        com->fprintf(stderr, " ... Writing Fluid mesh positions to %s\n", posFile);
        domain->writeVectorToFile(posFile, 0, 0.0, *X);
        delete [] posFile;
      }
    } else {
      mems = 0;
    }
}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::printf(int verbose, const char *format, ...)
{

  if (com->cpuNum() == 0 && verbose <= com->getMaxVerbose()) {
    va_list args;
    va_start(args, format);
    vfprintf(stdout, format, args);
    ::fflush(stdout);
    va_end(args);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::fprintf(FILE *fp, const char *format, ...)
{

  if (com->cpuNum() == 0) {
    va_list args;
    va_start(args, format);
    vfprintf(fp, format, args);
    ::fflush(fp);
    va_end(args);
  }

}

//------------------------------------------------------------------------------

template<int dim>
VarFcn *TsDesc<dim>::createVarFcn(IoData &ioData)
{

  VarFcn *vf = 0;
  fprintf(stderr,"ERROR: obsolete function createVarFcn(...) called!\n");
  exit(-1);
  return vf;

}

//------------------------------------------------------------------------------

template<int dim>
DistBcData<dim> *TsDesc<dim>::createBcData(IoData &ioData)
{

  DistBcData<dim> *bc = 0;

  if (ioData.eqs.type == EquationsData::NAVIER_STOKES &&
      ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
    if (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
        ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES)
      bc = new DistBcDataSA<dim>(ioData, varFcn, domain, *X);
    else if (ioData.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE)
      bc = new DistBcDataKE<dim>(ioData, varFcn, domain, *X);
  }
  else
    bc = new DistBcDataEuler<dim>(ioData, varFcn, domain, *X);

  if (!bc) {
    com->fprintf(stderr, "*** Error: no valid choice for BCs\n");
    exit(1);
  }

  return bc;

}

//------------------------------------------------------------------------------

template<int dim>
MeshMotionHandler *TsDesc<dim>::
createMeshMotionHandler(IoData &ioData, GeoSource &geoSource, MemoryPool *mp)
{

  MeshMotionHandler *_mmh = 0;

  if (ioData.problem.type[ProblemData::AERO]) {
    if (ioData.problem.type[ProblemData::ACCELERATED])
      _mmh = new AccAeroMeshMotionHandler(ioData, varFcn, bcData->getInletPrimitiveState(),
                                          geoSource.getMatchNodes(), domain, mp);
    else
      _mmh = new AeroMeshMotionHandler(ioData, varFcn, bcData->getInletPrimitiveState(),
                                       geoSource.getMatchNodes(), domain, mp);
    //check that algorithm number is consistent with simulation in special case RK2-CD
    // if C0 and RK2 then RK2DGCL is needed!
    if(_mmh->getAlgNum() == 20 || _mmh->getAlgNum() == 21 || _mmh->getAlgNum() == 22){
      if(ioData.ts.type == TsData::EXPLICIT &&
         (ioData.ts.expl.type == ExplicitData::RUNGE_KUTTA_2 ||
          ioData.ts.expl.type == ExplicitData::ONE_BLOCK_RK2 ||
          ioData.ts.expl.type == ExplicitData::ONE_BLOCK_RK2bis )){
        if(!((ioData.dgcl.normals    == DGCLData::EXPLICIT_RK2     || ioData.dgcl.normals == DGCLData::AUTO) &&
             (ioData.dgcl.velocities == DGCLData::EXPLICIT_RK2_VEL || ioData.dgcl.velocities == DGCLData::AUTO_VEL))){
          com->fprintf(stderr, "***Error: Computation of the normals or velocities (%d,%d)\n", ioData.dgcl.normals, ioData.dgcl.velocities);
          com->fprintf(stderr, "***       is not consistent with Aeroelastic algorithm\n");
          exit(1);
        }
      }
    }
  }
  else if (ioData.problem.type[ProblemData::FORCED]) {
    if (ioData.forced.type == ForcedData::HEAVING) {
      if (ioData.problem.type[ProblemData::ACCELERATED]){
        _mmh = new AccHeavingMeshMotionHandler(ioData, varFcn, bcData->getInletPrimitiveState(), domain);
      } else {
        _mmh = new HeavingMeshMotionHandler(ioData, domain);
      }
    } else if (ioData.forced.type  == ForcedData::PITCHING){
      if (ioData.problem.type[ProblemData::ACCELERATED]){
        _mmh = new AccPitchingMeshMotionHandler(ioData, varFcn, bcData->getInletPrimitiveState(), domain);
      } else {
        _mmh = new PitchingMeshMotionHandler(ioData, domain);
      }
    } else if (ioData.forced.type  == ForcedData::DEFORMING){
      if (ioData.problem.type[ProblemData::ACCELERATED]){
        _mmh = new AccDeformingMeshMotionHandler(ioData, varFcn, bcData->getInletPrimitiveState(), domain);
      } else {
        _mmh = new DeformingMeshMotionHandler(ioData, domain);
      }
    } else if (ioData.forced.type == ForcedData::SPIRALING){
      if (ioData.problem.type[ProblemData::ACCELERATED]){
        com->fprintf(stderr,"***Error: Accelerated Spiraling is not currently a supported problem type\n");
        exit(-1);
      } else {
        _mmh = new SpiralingMeshMotionHandler(ioData, domain);
      }
    }
  }
  else if (ioData.problem.type[ProblemData::ACCELERATED])
    _mmh = new AccMeshMotionHandler(ioData, varFcn, bcData->getInletPrimitiveState(), domain);
  else if (ioData.problem.type[ProblemData::ROLL])
    _mmh = new RigidRollMeshMotionHandler(ioData, bcData->getInletAngles(), domain);
  else if (ioData.problem.type[ProblemData::RBM])
    _mmh = new RbmExtractor(ioData, domain);


  return _mmh;

}

//------------------------------------------------------------------------------

template<int dim>
HeatTransferHandler* TsDesc<dim>::createHeatTransferHandler(IoData& iod, GeoSource& gs)

{

  HeatTransferHandler* _hth = 0;

  if (iod.problem.type[ProblemData::THERMO])
    _hth = new HeatTransferHandler(iod, gs.getMatchNodes(), domain);

  return _hth;

}

//------------------------------------------------------------------------------
template<int dim>
double TsDesc<dim>::recomputeResidual(DistSVec<double,dim> &F, DistSVec<double,dim> &Finlet)
{

  return spaceOp->recomputeResidual(F,Finlet);

}

//------------------------------------------------------------------------------
template<int dim>
void TsDesc<dim>::evaluateFluxAtMultipleSolutions(IoData &iod, char* best_soln)
{
  com->fprintf(stderr," ... In TsDesc<dim>::evaluateFluxAtMultipleSolutions ...\n");

  FILE *inFP = fopen(input->multiSolutions,"r");
  if (!inFP)  {
    com->fprintf(stderr, "*** Error: No solution data FILES in %s\n", input->multiSolutions);
    exit (-1);
  }
  int nData, _n;
  _n = fscanf(inFP, "%d",&nData);
  com->fprintf(stdout, "Reading %d Solutions for Flux Evaluation\n",nData);

  if (nData == 0) {
     strcpy(best_soln,"");
     return;
  }

  char solnFile1[500];
  char** solnFile = new char*[nData];
  for (int iData=0; iData < nData; ++iData)
    solnFile[iData] = new char[500];
  double* normF = new double[nData];
  double tmp, bestSoFar;

  for (int i=0; i < nData; ++i) {
    _n = fscanf(inFP, "%s", solnFile1);
    com->fprintf(stderr,"     solnFile = %s\n",solnFile1);
    strcpy(solnFile[i],solnFile1);
    domain->readVectorFromFile(solnFile[i], 0, 0, *Uic);

    spaceOp->computeResidual(*X, *A, *Uic, *F, timeState);
    spaceOp->applyBCsToResidual(*Uic, *F);
    tmp = (*F).norm();
    normF[i] = 0.5*tmp*tmp;
    if (i == 0 || tmp < bestSoFar){
       bestSoFar = tmp;
       strcpy(best_soln,solnFile[i]);
    }
  }
  fclose(inFP);

  if (com->cpuNum() == 0) {
    if (iod.output.transient.multiSolnFluxNorm[0] != 0){
      int dsp = strlen(iod.output.transient.prefix)+1;
      char* MultiSolnFluxNorm = new char[dsp + strlen(iod.output.transient.multiSolnFluxNorm)];
      sprintf(MultiSolnFluxNorm,"%s%s",iod.output.transient.prefix,iod.output.transient.multiSolnFluxNorm);

      FILE *fpMultiSolnFluxNorm = fopen(MultiSolnFluxNorm,"w");
      if (!fpMultiSolnFluxNorm) {
        fprintf(stderr,"*** Error: could not open \'%s'\n",MultiSolnFluxNorm);
      }
      fprintf(fpMultiSolnFluxNorm,"Solution FluxNorm Name\n");
      for (int iData=0; iData < nData; ++iData)
         fprintf(fpMultiSolnFluxNorm,"%d %e %s\n",iData,normF[iData],solnFile[iData]);
      delete[] MultiSolnFluxNorm;

       if (fpMultiSolnFluxNorm) fclose(fpMultiSolnFluxNorm);
    }
  }
}
//------------------------------------------------------------------------------
template<int dim>
void TsDesc<dim>::formInterpolationWeights(IoData &iod) {

  com->fprintf(stdout, " ... determining the initial condition via interpolation of the solutions from %s\n", iod.input.multiSolutionsParams);

  // read file, store solutions and parameters
  FILE *paramsFile = fopen(iod.input.multiSolutionsParams,"r");
  if (!paramsFile)  {
    com->fprintf(stderr, "*** Error: No solution data FILES in %s\n", iod.input.multiSolutionsParams);
    exit (-1);
  }
  int nData, _n, nParams;
  _n = fscanf(paramsFile, "%d",&nData);
  com->fprintf(stdout, " ... forming interpolated initial condition from %d precomputed states\n",nData);
  _n = fscanf(paramsFile, "%d",&nParams);
  com->fprintf(stdout, " ... %d-dimensional parameter space\n",nParams);
  std::vector<std::vector<double> > solutionsParameters;
  solutionsParameters.resize(nData);

  std::vector<double> minParamValues(nParams,1e16);
  std::vector<double> maxParamValues(nParams,-1e16);

  char solnFile[500];
  for (int iData=0; iData < nData; ++iData) {
    // skip solution name
    _n = fscanf(paramsFile, "%s", solnFile);
    // read parameters
    solutionsParameters[iData].resize(nParams);
    for (int iParam=0; iParam < nParams; ++iParam) {
      _n = fscanf(paramsFile, "%lf", &(solutionsParameters[iData][iParam]));
      minParamValues[iParam] = (solutionsParameters[iData][iParam] < minParamValues[iParam]) ? solutionsParameters[iData][iParam]: minParamValues[iParam];
      maxParamValues[iParam] = (solutionsParameters[iData][iParam] > maxParamValues[iParam]) ? solutionsParameters[iData][iParam]: minParamValues[iParam];
    }
  }
  fclose(paramsFile);

  for (int iData=0; iData < nData; ++iData) {
    for (int iParam=0; iParam < nParams; ++iParam) {
      double paramRange = maxParamValues[iParam] - minParamValues[iParam];
      solutionsParameters[iData][iParam] = (paramRange>1e-15) ? (solutionsParameters[iData][iParam] - minParamValues[iParam]) / paramRange : 0;
    }
  }

  // read parameters for current simulation
  std::vector<double> parameters;
  parameters.resize(nParams);

  FILE *inParams = fopen(iod.input.parameters,"r");
  if (!inParams)  {
    com->fprintf(stderr, "*** Error: could not open parameter file (%s)\n", iod.input.parameters);
    exit (-1);
  }
  int tmpNParams;
  _n = fscanf(inParams, "%d",&tmpNParams);
  if (tmpNParams!=nParams) {
    com->fprintf(stderr, "*** Error: mismatch in number of parameters (%d vs %d)\n",nParams, tmpNParams);
    exit(-1);
  }

  com->fprintf(stdout,"\n");
  for (int iParam=0; iParam < nParams; ++iParam) {
    _n = fscanf(inParams, "%lf", &(parameters[iParam]));
    com->fprintf(stdout," ... current operating point: param #%d = %e\n",iParam,parameters[iParam]);
    double paramRange = maxParamValues[iParam] - minParamValues[iParam];
    parameters[iParam] = (paramRange>1e-15) ? (parameters[iParam] - minParamValues[iParam]) / paramRange : 0;
  }

  // determine weighting
  // strategy: use a convex combination of the (primitive) states, with interpolatedICWeights proportional to
  //           the inverse of the distance btw. training parameters and current operating point
  com->fprintf(stdout, " ... calculating interpolation weights for initial condition\n");
  std::vector<double> distances;
  distances.resize(nData);
  double distanceExponent = iod.input.parametricDistanceExponent;
  if (distanceExponent!=1.0) {
    com->fprintf(stdout, " ... using exponent of %e for all distances\n", distanceExponent);
  }
  int maxInterpolatedSolutions = iod.input.maxInterpolatedSolutions;
  if (maxInterpolatedSolutions>=0) {
    com->fprintf(stdout, " ... using only the %d nearest solutions for interpolation\n", maxInterpolatedSolutions);
  }

  int reproductiveOperatingPoint = -1;
  for (int iData=0; iData<nData; ++iData) {
    distances[iData] = 0.0;
    for (int iParam=0; iParam < nParams; ++iParam) {
      distances[iData] += pow(parameters[iParam] - solutionsParameters[iData][iParam],2.0);
    }
    distances[iData] = pow(distances[iData],distanceExponent/2.0);
    //distances[iData] = pow(distances[iData],0.5); //use squared distances
    if (distances[iData] <= 1e-10) {
      if (reproductiveOperatingPoint>=0) {
        com->fprintf(stderr, "*** ERROR: More than one training solution matches the current operating point... exiting.\n");
        exit(-1);
      }
      reproductiveOperatingPoint = iData;
    }
  }


  if (maxInterpolatedSolutions>0 && maxInterpolatedSolutions<nData) {
    std::vector<double> sortedDistances(distances);
    std::stable_sort(sortedDistances.begin(), sortedDistances.end());
    int count = nData;
    for (int iData=0; iData<nData; ++iData) {
      if (distances[iData] > sortedDistances[maxInterpolatedSolutions-1]) {
        --count;
        distances[iData] = 1e16;
      }
    }
    if (count != maxInterpolatedSolutions)
      com->fprintf(stderr, "*** Warning: using %d interpolated solutions instead of %d\n", count, maxInterpolatedSolutions);
  }

  interpolatedICWeights.resize(nData, 0.0);

  if (reproductiveOperatingPoint<0) {
    //rbf_interp_nd_test04( ); // can use radial basis functions if necessary
    double weightSum = 0.0;
    for (int iData=0; iData<nData; ++iData) {
      interpolatedICWeights[iData] = 1.0/distances[iData];
      weightSum += interpolatedICWeights[iData];
    }
    for (int iData=0; iData<nData; ++iData) {
      interpolatedICWeights[iData] = interpolatedICWeights[iData]/weightSum;
    }
  } else {
    com->fprintf(stdout, " ... using the training solution corresponding to the current operating point...\n");
    interpolatedICWeights[reproductiveOperatingPoint]=1.0;
  }
  for (int iData=0; iData<nData; ++iData) {
    com->fprintf(stdout, " ... weight[%d]=%e...\n", iData, interpolatedICWeights[iData]);
  }

  //this->setInterpWeightsForMultiIC(interpolatedICWeights);

 }

//------------------------------------------------------------------------------
template<int dim>
void TsDesc<dim>::formInterpolatedInitialCondition(DistSVec<double,dim> *U, IoData &iod)  {
  // overloaded for ImplicitGnatTsDesc, which needs to handle this a bit differently

  com->fprintf(stdout, " ... entering formInterpolatedIC, reading %s\n", iod.input.multiSolutionsParams);
  FILE *paramsFile = fopen(iod.input.multiSolutionsParams,"r");
  if (!paramsFile)  {
    com->fprintf(stderr, "*** Error: No solution data FILES in %s\n", iod.input.multiSolutionsParams);
    exit (-1);
  }
  int nData, _n, nParams;
  double tmp;
  _n = fscanf(paramsFile, "%d", &nData);
  _n = fscanf(paramsFile, "%d", &nParams);
  com->fprintf(stdout, " ... reading %d solutions for interpolation (%d params)\n", nData, nParams);

  char solnFile[500];
  DistSVec<double,dim> Utmp(domain->getNodeDistInfo());
  *U = 0.0;
  for (int iData=0; iData < nData; ++iData) {
    // read solution
    _n = fscanf(paramsFile, "%s", solnFile);
    domain->readVectorFromFile(solnFile, 0, 0, Utmp);
    // add this vector's contribution to U
    *U += interpolatedICWeights[iData]*Utmp;
    // skip through parameters
    for (int iParam=0; iParam < nParams; ++iParam) {
      _n = fscanf(paramsFile, "%lf", &tmp);
    }
  }
  fclose(paramsFile);

 // DistSVec<double,dim> Vtmp(domain->getNodeDistInfo());
 // DistSVec<double,dim> V(domain->getNodeDistInfo());
 // U = weighted sum of solutions in primitive variables (to avoid negative pressure issues)
 // V = 0.0;
 // for (int iData=0; iData<nData; ++iData) {
 //   // add this vector's contribution to U
 //   _n = fscanf(inFP, "%s", solnFile);
 //   domain->readVectorFromFile(solnFile, 0, 0, Utmp);
 //   varFcn->conservativeToPrimitive(Utmp, Vtmp);
 //   V += interpolatedICWeights[iData]*Vtmp;
 // }
 // varFcn->primitiveToConservative(V, *U);

}

//------------------------------------------------------------------------------
template<int dim>
void TsDesc<dim>::setupTimeStepping(DistSVec<double,dim> *U, IoData &iod)
{
  char * name = new char[500];

  geoState->setup2(timeState->getData());
  if (iod.input.solutions[0] == 0 && iod.input.multiSolutionsParams[0] != 0) {
    // interpolate stored solutions to find an appropriate initial condition
    double icInterpTime = timer->getTime();
    if (interpolatedICWeights.size()==0) formInterpolationWeights(iod);
    this->formInterpolatedInitialCondition(U, iod);// modifies U
    timer->addICInterpTime(icInterpTime);
    timeState->setup("", *X, *U, *U, iod); // pass modified U instead of Ufarfield, name=NULL ensures this state is used
  } else if ( iod.input.solutions[0] == 0 && iod.input.multiSolutions[0] != 0 && iod.problem.solveWithMultipleICs) {
    // U has already been set to the appropriate initial condition (in TsSolver)
    timeState->setup("", *X, *U, *U, iod);
  } else { // initial condition will be read from a file
     if ( iod.input.solutions[0] == 0 && iod.input.multiSolutions[0] != 0){
       evaluateFluxAtMultipleSolutions(iod,name);
     } else {
       strcpy(name,input->solutions);
     }
    timeState->setup(name, *X, bcData->getInletBoundaryVector(), *U, iod);
  }

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(mmh);
  DeformingMeshMotionHandler* _dmmh = dynamic_cast<DeformingMeshMotionHandler*>(mmh);
  HeavingMeshMotionHandler* _hmmh = dynamic_cast<HeavingMeshMotionHandler*>(mmh);
  SpiralingMeshMotionHandler* _smmh = dynamic_cast<SpiralingMeshMotionHandler*>(mmh);
  PitchingMeshMotionHandler* _pmmh = dynamic_cast<PitchingMeshMotionHandler*>(mmh);

  if (_mmh)
    _mmh->setup(&restart->frequency, &data->maxTime, postOp, *X, *U);
  else if (_dmmh)
    _dmmh->setup(*X);
  else if (_hmmh)
    _hmmh->setup(*X);
  else if (_smmh)
    _smmh->setup(*X);
  else if (_pmmh)
    _pmmh->setup(*X);

  if (hth)
    hth->setup(&restart->frequency, &data->maxTime);

  *Xs = *X;

  initializeFarfieldCoeffs();
}

//------------------------------------------------------------------------------

template<int dim>
double TsDesc<dim>::computeTimeStep(int it, double *dtLeft, DistSVec<double,dim> &U, double angle)
{
  double t0 = timer->getTime();

  com->barrier();
  this->data->allowstop = this->timeState->allowcflstop;
  data->computeCflNumber(it - 1, data->residual / restart->residual, angle);
  int numSubCycles = 1;

  double dt = 0.0;
  if(failSafeFlag == false){
    if(timeStepCalculation == TsData::CFL || it==1) {
      dt = timeState->computeTimeStep(data->cfl, data->dualtimecfl, dtLeft, &numSubCycles, *geoState, *X, *A, U);
    } else { //time step size with error estimation
      dt = timeState->computeTimeStep(it, dtLeft, &numSubCycles);
    }
  }
  else { //if time step is repeated

    dt = this->timeState->computeTimeStepFailSafe(dtLeft, &numSubCycles);
  }

  if(timeStepCalculation == TsData::ERRORESTIMATION && it == 1) {
    this->timeState->setDtMin(dt * data->getCflMinOverCfl0());
  }

  if (problemType[ProblemData::UNSTEADY])
    com->printf(5, "Global dt: %g (remaining subcycles = %d)\n", dt*refVal->time, numSubCycles);

  timer->addFluidSolutionTime(t0);
  timer->addTimeStepTime(t0);

  return dt;
}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::getNumParam(int &numParam, int &actvar, double &steadyTol)
{
  if (mmh) mmh->getNumParam(numParam, actvar, steadyTol);
}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::sendNumParam(int numParam)
{
  if (mmh) mmh->sendNumParam(numParam);
}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::getRelResidual(double &relres)
{
  if (mmh) mmh->getRelResidual(relres);
}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::cmdCom(bool *lastIt)
{
  if (mmh) mmh->cmdCom(lastIt);
}

//------------------------------------------------------------------------------

template<int dim>
double TsDesc<dim>::computePositionVector(bool *lastIt, int it, double t, DistSVec<double,dim> &U)
{
  double dt = 0.0;
  //this->com->printf(9, "deubgging: entering computePositionVector()\n");
  if (mmh && mmh->structureSubcycling()) {
    double dtleft = 0.0;
    double dtf = this->computeTimeStep(it+1, &dtleft, U);
    mmh->storeFluidSuggestedTimestep(dtf);
  }

  if (mmh) {
    double t0 = timer->getTime();
    dt = mmh->updateStep1(lastIt, it, t, bcData->getVelocityVector(), *Xs, &data->maxTime);
    timer->addMeshSolutionTime(t0);
  }

  if (hth) {
    double dth = hth->updateStep1(lastIt, it, bcData->getTemperatureVector());
    if (!mmh)
      dt = dth;
  }

  if (mmh) {
    double t0 = timer->getTime();
    //this->com->printf(9, "debuggin: entering mmh->updateStep2()\n");
    mmh->updateStep2(lastIt, it, t, bcData->getVelocityVector(), *Xs);
    timer->addMeshSolutionTime(t0);
  }

  if (hth) {
    hth->updateStep2(lastIt, it, bcData->getTemperatureVector());
  }

  return dt;

}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::setMeshSensitivitySolverPositionVector()
{
  if(mmh) {
    mmh->setPositionVector(*Xs);
  }
}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::receiveBoundaryPositionSensitivityVector(DistSVec<double,3> &dXdSb, bool applyScale)
{
  if (mmh) {
	mmh->updateDStep2(*Xs,dXdSb, applyScale);
  }
}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::negotiate()
{
  if (mmh)  mmh->negotiate();
}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::sendForceSensitivity(DistSVec<double,3> *dFdS, bool applyScale)
{
  if (mmh)  mmh->sendForceSensitivity(dFdS, applyScale);
}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::interpolatePositionVector(double dt, double dtLeft)
{
  if (!mmh) return;

//  EmbeddedMeshMotionHandler* _mmh = dynamic_cast<EmbeddedMeshMotionHandler*>(mmh);
//  if (_mmh) return;

  geoState->interpolate(dt, dtLeft, *Xs, *X);

}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::computeMeshMetrics(int it)
{
//  EmbeddedMeshMotionHandler* _mmh = dynamic_cast<EmbeddedMeshMotionHandler*>(mmh);

  if (mmh) {
    if (it >= 0) com->fprintf(stderr, "GeoState Computing for it %d\n", it);
    double t0 = timer->getTime();
    geoState->compute(timeState->getData(), bcData->getVelocityVector(), *X, *A);
    timer->addMeshMetricsTime(t0);
    timer->addFluidSolutionTime(t0);
  }

  if (mmh || hth) {
    bcData->update(*X);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::updateStateVectors(DistSVec<double,dim> &U, int it)
{

  geoState->update(*X, *A);
  timeState->update(U);

  spaceOp->updateFixes();
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
bool TsDesc<dim>::checkForLastIteration(IoData &ioData, int it, double t, double dt, DistSVec<double,dim> &U)
{

  if ((ioData.eqs.type == EquationsData::NAVIER_STOKES) && (ioData.sa.fres)) {
    if (!problemType[ProblemData::UNSTEADY]) {
      bool forceconv = false;
      if (ioData.sa.avgsIt)
        forceconv = monitorAvgForceConvergence(ioData, it, U);
      else
        forceconv = monitorForceConvergence(ioData, it, U);
      bool fluidconv = monitorConvergence(it, U);
      if ((forceconv || fluidconv) || it >= data->maxIts) {
        if (forceconv)
          com->fprintf(stderr,"\n***** Residual of the aerodynamic force is satisfied \n\n");
        else if (fluidconv)
          com->fprintf(stderr,"\n***** Residual of the fluid solution is satisfied \n\n");
        else if (it >= data->maxIts)
          com->fprintf(stderr,"\n***** Maximum number of iteration is reached \n\n");
        return true;
      }
    }
  }
  else {
    if (!problemType[ProblemData::UNSTEADY] && monitorConvergence(it, U))
      return true;
  }

  if (!problemType[ProblemData::AERO] && !problemType[ProblemData::THERMO] && it >= data->maxIts) return true;

  if (problemType[ProblemData::UNSTEADY] )
    if(t >= data->maxTime - 0.01 * dt)
      return true;

  return false;

}

//------------------------------------------------------------------------------

template<int dim>
int TsDesc<dim>::checkSolution(DistSVec<double,dim> &U)
{

  int ierr = 0;

  if (dim == 6)
    ierr = domain->template
      clipSolution<dim,1>(clippingType, wallType, varFcn, bcData->getInletConservativeState(), U);
  else if (dim == 7)
    ierr = domain->template
      clipSolution<dim,2>(clippingType, wallType, varFcn, bcData->getInletConservativeState(), U);
  else
    ierr = domain->checkSolution(varFcn, U);

  return ierr;

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void TsDesc<dim>::fixSolution(DistSVec<double,dim> &U, DistSVec<double,dim> &dU)
{

  if (fixSol == 1)
    domain->fixSolution(varFcn, U, dU);

}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::setupOutputToDisk(IoData &ioData, bool *lastIt, int it, double t,
                                    DistSVec<double,dim> &U)
{
  if (it == data->maxIts) {
    *lastIt = true;
  } else if (!ioData.sa.fsiFlag) {
    monitorInitialState(it, U);
  }

  output->setMeshMotionHandler(ioData, mmh);
  output->openAsciiFiles();
  timer->setSetupTime();
  output->cleanProbesFile();

  if (it == 0) {
    // First time step: compute GradP before computing forces
    spaceOp->computeGradP(*X, *A, U);

    if (wallRecType==BcsWallData::CONSTANT) {
      output->writeForcesToDisk(*lastIt, it, 0, 0, t, 0.0, restart->energy, *X, U);
    } else { //wallRecType == EXACT_RIEMANN
      output->writeForcesToDisk(*riemann1, *lastIt, it, 0, 0, t, 0.0, restart->energy, *X, U);
    }

    double fluxNorm = 0.5*(data->residual)*(data->residual);

    output->writeLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, restart->energy, *X, U);
    output->writeMatchPressureToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, restart->energy, *X, *A, U, timeState);
    output->writeMatchStateToDisk(ioData, it, 0.0, 0.0, U, *A);
    output->writeFluxNormToDisk(it, 0, 0, t, fluxNorm);
    output->writeHydroForcesToDisk(*lastIt, it, 0, 0, t, 0.0, restart->energy, *X, U);
    output->writeHydroLiftsToDisk(ioData, *lastIt, it, 0, 0, t, 0.0, restart->energy, *X, U);
    output->writeResidualsToDisk(it, 0.0, 1.0, data->cfl);
    output->writeMaterialVolumesToDisk(it, 0.0, *A);
    output->writeMaterialConservationScalarsToDisk(it, 0.0, U,*A);
    output->writeCPUTimingToDisk(*lastIt, it, t, timer);
    output->writeBinaryVectorsToDisk(*lastIt, it, t, *X, *A, U, timeState);
    output->writeAvgVectorsToDisk(*lastIt, it, t, *X, *A, U, timeState);
    output->writeHeatFluxesToDisk(*lastIt, it, 0, 0, t, 0.0, restart->energy, *X, U);
    restart->writeKPtracesToDisk(ioData, *lastIt, it, t, *X, *A, U, timeState, domain, postOp);
    writeStateRomToDisk(it, 0.0);
    writeErrorToDisk(it, 0.0);
  }
}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::outputToDisk(IoData &ioData, bool* lastIt, int it, int itSc, int itNl,
                               double t, double dt, DistSVec<double,dim> &U)
{

  com->globalSum(1, &interruptCode);
  if (interruptCode) {
    *lastIt = true;
  }

  double cpu = timer->getRunTime();
  double res = data->residual / restart->residual;
  double fluxNorm = 0.5*(data->residual)*(data->residual);

  output->writeLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, restart->energy, *X, U);
  output->writeMatchPressureToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, restart->energy, *X, *A, U, timeState);
  output->writeMatchStateToDisk(ioData, it, t, cpu, U, *A);
  output->writeFluxNormToDisk(it, itSc, itNl, t, fluxNorm);
  output->writeHydroForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, restart->energy, *X, U);
  output->writeHydroLiftsToDisk(ioData, *lastIt, it, itSc, itNl, t, cpu, restart->energy, *X, U);
  output->writeResidualsToDisk(it, cpu, res, data->cfl);
  writeStateRomToDisk(it, cpu);
  output->writeMaterialVolumesToDisk(it, t, *A);
  output->writeMaterialConservationScalarsToDisk(it, 0.0, U,*A);
  output->writeCPUTimingToDisk(*lastIt, it, t, timer);
  writeErrorToDisk(it, cpu);
  output->writeBinaryVectorsToDisk(*lastIt, it, t, *X, *A, U, timeState);
  output->writeAvgVectorsToDisk(*lastIt, it, t, *X, *A, U, timeState);
  output->writeProbesToDisk(*lastIt, it, t, *X, *A, U, timeState,fluidIdDummy);
  restart->writeToDisk<dim,1>(com->cpuNum(), *lastIt, it, t, dt, *timeState, *geoState);
  output->writeHeatFluxesToDisk(*lastIt, it, itSc, itNl, t, cpu, restart->energy, *X, U);

  restart->writeKPtracesToDisk(ioData, *lastIt, it, t, *X, *A, U, timeState, domain, postOp);

  this->output->updatePrtout(t);
  if (*lastIt) {

    if (strcmp(ioData.input.convergence_file,"") != 0) {

      this->varFcn->conservativeToPrimitive(U, *this->V);
      computeConvergenceInformation(ioData,ioData.input.convergence_file,*this->V);
    }


    timer->setRunTime();
    if (com->getMaxVerbose() >= 2) {
      timer->print(domain->getStrTimer());
    }

    if(ioData.problem.alltype != ProblemData::_SHAPE_OPTIMIZATION_ &&
       ioData.problem.alltype != ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_ &&
       ioData.problem.alltype != ProblemData::_ROM_SHAPE_OPTIMIZATION_ &&
       ioData.problem.alltype != ProblemData::_SENSITIVITY_ANALYSIS_ ) { //TODO CHECK if really needed
      output->closeAsciiFiles();
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::outputForces(IoData &ioData, bool* lastIt, int it, int itSc, int itNl,
                               double t, double dt, DistSVec<double,dim> &U)  {

  double cpu = timer->getRunTime();
  if (wallRecType==BcsWallData::CONSTANT)
    output->writeForcesToDisk(*lastIt, it, itSc, itNl, t, cpu, restart->energy, *X, U);
  else //wallRecType==EXACT_RIEMANN
    output->writeForcesToDisk(*riemann1, *lastIt, it, itSc, itNl, t, cpu, restart->energy, *X, U);
}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::outputPositionVectorToDisk(DistSVec<double,dim> &U)
{

  *X = *Xs;

  if (mmh)  {
    int algNum = mmh->getAlgNum();
    if (algNum == 8)  return;
  }

  domain->writeVectorToFile(restart->positions[0], 0, 0.0, *Xs, &(refVal->tlength));

  if(mmh && mmh->getAlgNum() == 1) {
    output->writeDisplacementVectorToDisk(1, 1.0, *X, U);
  }

  timer->setRunTime();
  if (com->getMaxVerbose() >= 2)
    timer->print(domain->getStrTimer());

  DistVec<double> As(getVecInfo());
  int ierr = domain->computeControlVolumes(refVal->tlength, *Xs, As);
#ifdef YDEBUG
  if(ierr) {
    const char* output = "elementvolumecheck";
    ofstream out(output, ios::out);
    if(!out) { cerr << "Error: cannot open file" << output << endl;  exit(-1); }
    out << ierr << endl;
    out.close();
    exit(-1);
  }
#endif

}

//------------------------------------------------------------------------------
/*
template<int dim>
void TsDesc<dim>::outputPositionSensitivityVectorToDisk(DistSVec<double,dim> &dUds)
{

  if (mmh)  {
    int algNum = mmh->getAlgNum();
    if (algNum == 8)  return;
  }

  output->writePositionSensitivityVectorToDisk(0, 0.0, *dXds);

  timer->setRunTime();
  if (com->getMaxVerbose() >= 2)
    timer->print(domain->getStrTimer());

}
*/
//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::resetOutputToStructure(DistSVec<double,dim> &U)
{

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(mmh);
  if (_mmh)
    _mmh->resetOutputToStructure(postOp, *X, U);

}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::updateOutputToStructure(double dt, double dtLeft,
                                          DistSVec<double,dim> &U)
{

  if (mmh) {
    double work[2];
    mmh->computeInterfaceWork(dt, postOp, geoState->getXn(), timeState->getUn(), *X, U, work);
    restart->energy[0] += work[0];
    restart->energy[1] += work[1];
  }

  AeroMeshMotionHandler* _mmh = dynamic_cast<AeroMeshMotionHandler*>(mmh);
  if (_mmh)
    _mmh->updateOutputToStructure(dt, dtLeft, postOp, *X, U);

  if (hth)
    hth->updateOutputToStructure(dt, dtLeft, postOp, *X, U);

}

//------------------------------------------------------------------------------

template<int dim>
double TsDesc<dim>::computeResidualNorm(DistSVec<double,dim>& U)
{
  if (wallRecType==BcsWallData::CONSTANT)
    spaceOp->computeResidual(*X, *A, U, *R, timeState);
  else //wallRecTyp == ExactRiemann
    spaceOp->computeResidual(riemann1, *X, *A, U, *R, timeState);

  spaceOp->applyBCsToResidual(U, *R);

  double res = 0.0;
  if (data->resType == -1){
    res = (*R)*(*R);


  }else{
    int iSub;
    const DistInfo& distInfo = R->info();
#ifndef MPI_OMP_REDUCTION
    double* allres = reinterpret_cast<double*>(alloca(sizeof(double) * distInfo.numGlobSub));
    for (iSub=0; iSub<distInfo.numGlobSub; ++iSub)
      allres[iSub] = 0.0;
#endif

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
    for (iSub=0; iSub<distInfo.numLocSub; ++iSub) {
      double (*r)[dim] = R->subData(iSub);
      bool* flag = R->getMasterFlag(iSub);
      double locres = 0.0;
      for (int i=0; i<R->subSize(iSub); ++i) {
        if (flag[i])
          locres += r[i][data->resType]*r[i][data->resType];
      }
#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif
    }

#ifdef MPI_OMP_REDUCTION
    distInfo.com->globalSum(1, &res);
#else
    distInfo.com->globalSum(distInfo.numGlobSub, allres);
    res = 0.0;
    for (iSub=0; iSub<distInfo.numGlobSub; ++iSub)
      res += allres[iSub];
#endif
  }

  return sqrt(res);

}

//------------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::monitorInitialState(int it, DistSVec<double,dim> &U)
{

  if (!problemType[ProblemData::UNSTEADY]) {
    double trhs = timer->getTimeSyncro();
    data->residual = computeResidualNorm(U);
    trhs = timer->getTimeSyncro() - trhs;
    if (it == 0)
      restart->residual = data->residual;
    if (data->resType == -1)
      com->printf(2, "Spatial residual norm = %.12e\n", data->residual);
    else
      com->printf(2, "Spatial residual norm[%d] = %.12e\n", data->resType, data->residual);
    com->printf(2, "Time for one residual evaluation: %f s\n", trhs);
  }

  com->printf(2, "\n");

}

//------------------------------------------------------------------------------

template<int dim>
bool TsDesc<dim>::monitorConvergence(int it, DistSVec<double,dim> &U)
{

  // For multigrid, monitorConvergence() was already called in the
  // smoothing process.  It was called with it == 0
  // in this case we do not need to recompute it
  //
  if (!isMultigridTsDesc || it == 0)
    data->residual = computeResidualNorm(U);

  if ((problemType[ProblemData::AERO] || problemType[ProblemData::THERMO]) && (it == 1 || it == 2))
    restart->residual = data->residual;

  if (data->residual == 0.0 || data->residual < data->eps * restart->residual || data->residual < data->epsabs)
    return true;
  else
    return false;

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
bool TsDesc<dim>::monitorForceConvergence(IoData &ioData, int it, DistSVec<double,dim> &U)
{

  double forceNorm_n;
  double resForce;

  int nSurfs = postOp->getNumSurf();

  Vec3D x0, F, M;

  Vec3D *Fi = new Vec3D[nSurfs];
  Vec3D *Mi = new Vec3D[nSurfs];
  Vec3D *Fv = new Vec3D[nSurfs];
  Vec3D *Mv = new Vec3D[nSurfs];

  x0[0] = ioData.output.transient.x0;
  x0[1] = ioData.output.transient.y0;
  x0[2] = ioData.output.transient.z0;

  postOp->computeForceAndMoment(x0, *this->X, U, 0, Fi, Mi, Fv, Mv);

  F = 0.0;
  M = 0.0;

  F = Fi[0] + Fv[0];
  M = Mi[0] + Mv[0];

  forceNorm_n = sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2]);

  resForce = fabs(forceNorm_n - forceNorm) / forceNorm_n;

  forceNorm = forceNorm_n;

  com->fprintf(stderr,"\n***** It = %d, Force residual = %e, Target = %e\n\n", it, resForce, ioData.sa.fres);

  if ((resForce == 0.0) || (resForce < ioData.sa.fres))
    return true;
  else
    return false;

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
bool TsDesc<dim>::monitorAvgForceConvergence(IoData &ioData, int it, DistSVec<double,dim> &U)
{

  double avgForceNorm;
  double forceNorm_n;
  double resForce;

  if (forceNorms == 0)
    forceNorms = new double[ioData.sa.avgsIt];

  int nSurfs = postOp->getNumSurf();

  Vec3D x0, F, M;

  Vec3D *Fi = new Vec3D[nSurfs];
  Vec3D *Mi = new Vec3D[nSurfs];
  Vec3D *Fv = new Vec3D[nSurfs];
  Vec3D *Mv = new Vec3D[nSurfs];

  x0[0] = ioData.output.transient.x0;
  x0[1] = ioData.output.transient.y0;
  x0[2] = ioData.output.transient.z0;

  postOp->computeForceAndMoment(x0, *this->X, U, 0, Fi, Mi, Fv, Mv);

  F = 0.0;
  M = 0.0;

  F = Fi[0] + Fv[0];
  M = Mi[0] + Mv[0];

  forceNorm_n = sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2]);

  if (iTotal < ioData.sa.avgsIt) {
    forceNorms[iForce] = forceNorm_n;
    ++iForce;
    avgForceNorm = 0.0;
    for (int i = 0; i < iForce; ++i)
      avgForceNorm += forceNorms[i];
    avgForceNorm *= 1.0/double(iForce);
  }
  else {
    if (iForce == ioData.sa.avgsIt)
      iForce = 0;
    forceNorms[iForce] = forceNorm_n;
    ++iForce;
    avgForceNorm = 0.0;
    for (int i = 0; i < ioData.sa.avgsIt; ++i)
      avgForceNorm += forceNorms[i];
    avgForceNorm *= 1.0/double(ioData.sa.avgsIt);
  }

  if (iTotal)
    resForce = fabs(forceNorm_n - avgForceNorm) / avgForceNorm;
  else
    resForce = 1.0;

  ++iTotal;

  com->fprintf(stderr,"\n***** It = %d, Force residual = %e, Target = %e\n\n", it, resForce, ioData.sa.fres);

  if ((resForce == 0.0) || (resForce < ioData.sa.fres))
    return true;
  else
    return false;

}

//-----------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::updateGhostFluid(DistSVec<double,dim> &U, Vec3D& totalForce, double dt)
{
/*
  if (eulerFSI)
    eulerFSI->updateGhostFluid(X, U, totalForce, dt);
*/
}

//-----------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::printNodalDebug(int globNodeId, int identifier, DistSVec<double,dim> *U, DistVec<int> *Id, DistVec<int> *Id0)
{ //Kevin:For debug only!
  int nSub = domain->getNumLocSub();
  SubDomain **sub = domain->getSubDomain();
  for(int iSub=0; iSub<nSub; iSub++) {
    int* locToGlob = sub[iSub]->getNodeMap();
    for(int i=0; i<(*U)(iSub).size(); i++)
      if(locToGlob[i]+1==globNodeId) {
        fprintf(stderr,"*** %d Node %d: U = %e %e %e %e %e ", identifier, globNodeId,
                (*U)(iSub)[i][0], (*U)(iSub)[i][1], (*U)(iSub)[i][2], (*U)(iSub)[i][3], (*U)(iSub)[i][4]);
        if(Id)
          fprintf(stderr,", Id = %d ", (*Id)(iSub)[i]);
        if(Id0)
          fprintf(stderr,", Id0 = %d ", (*Id0)(iSub)[i]);
        fprintf(stderr,"\n");
      }
  }
}

//----------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::computeDistanceToWall(IoData &ioData)
{
  // Nothing to do here by default.
}

//----------------------------------------------------------------------------


template<int dim>
void TsDesc<dim>::updateFarfieldCoeffs(double dt)
{
  if(!modifiedGhidaglia) return;
  int nSub = domain->getNumLocSub();
  SubDomain **sub = domain->getSubDomain();
  int iSub;
}

//----------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::updateBoundaryExternalState()
{
  if(!modifiedGhidaglia) return;
  int nSub = domain->getNumLocSub();
  SubDomain **sub = domain->getSubDomain();
  int iSub;
}

//----------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::initializeFarfieldCoeffs()
{
  if(!modifiedGhidaglia) return;
  double *Vin = bcData->getInletPrimitiveState();
  double soundspeed = varFcn->computeSoundSpeed(Vin);
  double gamma;
  if (varFcn->getType() == VarFcnBase::STIFFENEDGAS ||
      varFcn->getType() == VarFcnBase::PERFECTGAS)
    gamma = varFcn->getGamma();
  else
    gamma = varFcn->getBetaWater();

  double HH_init = -2.0*soundspeed/(gamma - 1.0);

  *(bcData->getBoundaryStateHH()) = HH_init;

  this->domain->maskHHVector(*(bcData->getBoundaryStateHH()));
}

template<int dim>
void TsDesc<dim>::computeConvergenceInformation(IoData &ioData, const char* file, DistSVec<double,dim>& U) {

  DistSVec<double,dim> Uexact(U);
  OneDimensional::read1DSolution(ioData,file, Uexact,
				 (DistSVec<double,1>*)0,
				 NULL,//&fluidSelector,
				 spaceOp->getVarFcn(),
				 *this->X,
				 *this->domain,
				 OneDimensional::ModeU,
				 false) ;

  double error[dim];
  double refs[dim] = {ioData.ref.rv.density, ioData.ref.rv.velocity,
		       ioData.ref.rv.velocity, ioData.ref.rv.velocity,
		      ioData.ref.rv.pressure};

  double tot_error = 0.0;
  DistVec<int> nnv(U.info());
  nnv = 1;
  int nNodes = nnv.sum();

  this->domain->computeL1Error(U,Uexact,*this->A,error);
  for (int k = 0; k < dim; ++k) {
    tot_error += error[k] / nNodes;
    this->domain->getCommunicator()->fprintf(stdout,"L1 error [%d]: %lf\n", k, error[k]*refs[k] / nNodes);
  }
  this->domain->getCommunicator()->fprintf(stdout,"L1 error (total): %lf\n", tot_error);

  tot_error = 0.0;
  this->domain->computeLInfError(U,Uexact,error);
  for (int k = 0; k < dim; ++k) {
    tot_error = max(error[k],tot_error);
    this->domain->getCommunicator()->fprintf(stdout,"Linf error [%d]: %lf\n", k, error[k]*refs[k]);
  }
  this->domain->getCommunicator()->fprintf(stdout,"Linf error (total): %lf\n", tot_error);


}

//----------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::performPostProForState(DistSVec<double,dim> &outVec, int tmpIt)
{ // public function that performs post processing on a state vector. Used during Nonlinear ROM preprocessing
  bool tmpLastIt = false;
  output->writeBinaryVectorsToDisk(tmpLastIt, tmpIt, double(tmpIt), *X, *A, outVec, timeState);
}

//----------------------------------------------------------------------------

template<int dim>
void TsDesc<dim>::readICFromDisk(char * solnFile, int iData, int nData, DistSVec<double,dim> &U) {
    com->fprintf(stdout,"\n\n\nInitial condition %d of %d (%s)\n", iData+1, nData, solnFile);
    domain->readVectorFromFile(solnFile, 0, 0, U);
}
