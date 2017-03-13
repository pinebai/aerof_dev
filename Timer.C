#include <cstdio>
#include <cstdlib>
#include <alloca.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <Timer.h>
#include <Communicator.h>

#include <DebugTools.h>

//------------------------------------------------------------------------------

Timer::Timer(Communicator *communicator) : com(communicator)
{

  ioData = 0;
  initialTime = getTime();

  numTimings = (int)NUMTIMINGS;
  
  counter = new int[numTimings];
  data = new double[numTimings];

  for (int i=0; i<numTimings; ++i) {
    counter[i] = 0;
    data[i] = 0.0;
  }

}

//------------------------------------------------------------------------------

Timer::~Timer()
{

  if (counter) delete [] counter;
  if (data) delete [] data;

}

//------------------------------------------------------------------------------

double Timer::getTime()
{

  static double micro = 1.e-6;

  timeval tp;
  struct timezone tz;
  
  gettimeofday(&tp, &tz);

  // return 1000.0*tp.tv_sec + tp.tv_usec/1000.0;
  return double(tp.tv_sec) + double(tp.tv_usec) * micro;

}

//------------------------------------------------------------------------------

double Timer::getTimeSyncro()
{

  if (com) com->barrier();

  return getTime();

}

//------------------------------------------------------------------------------

double Timer::getRunTime()
{

  return getTime() - data[setup] - initialTime;

}

//------------------------------------------------------------------------------
/*
double Timer::getCpuTime()
{

  static struct rusage r;
  static double micro = (1./1.e6);

  double t;

  getrusage(RUSAGE_SELF,&r);

  // maximum resident set size utilized (in kilobytes)
  // size = r.ru_maxrss;

  // amount of time spent executing in user mode (in seconds)
  t = double(r.ru_utime.tv_sec) + double(r.ru_utime.tv_usec) * micro;

  // amount of time spent in the system executing on behalf of the process(es)
  t += double(r.ru_stime.tv_sec) + double(r.ru_stime.tv_usec) * micro;

  return t;

}
*/

//------------------------------------------------------------------------------
void Timer::setIoData(IoData &_ioData)
{

ioData = &_ioData;

/*if (ioData->problem.alltype == ProblemData::_ROB_CONSTRUCTION_)
  numTimings += 4;
else if (ioData->problem.alltype == ProblemData::_ROM_AEROELASTIC_)
  numTimings += 3;  
*/
}

//------------------------------------------------------------------------------

void Timer::setSetupTime() 
{ 

  counter[setup]++;
  data[setup] = getTime() - initialTime; 

}

//------------------------------------------------------------------------------

void Timer::setRunTime() 
{ 

  counter[run]++;
  data[run] = getRunTime(); 

}

//------------------------------------------------------------------------------

double Timer::addTimeStepTime(double t0) 
{ 

  double t = getTime() - t0;
  
  counter[timeStep]++;
  data[timeStep] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addNodalWeightsTime(double t0) 
{ 

  double t = getTime() - t0;
  
  counter[nodalWeights]++;
  data[nodalWeights] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addNodalGradTime(double t0) 
{ 
  
  double t = getTime() - t0;

  counter[nodalGrad]++;
  data[nodalGrad] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addFiniteVolumeTermTime(double t0) 
{ 

  double t = getTime() - t0;

//  if (com->cpuNum() == 0)
 //   DebugTools::PrintBacktrace();

  counter[fvTerm]++;
  data[fvTerm] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addFiniteElementTermTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[feTerm]++;
  data[feTerm] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addFiniteVolumeJacTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[fvJac]++;
  data[fvJac] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addFiniteElementJacTime(double t0) 
{ 
  
  double t = getTime() - t0;

  counter[feJac]++;
  data[feJac] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addVMSLESTime(double t0)
{

  double t = getTime() - t0;

  counter[vms]++;
  data[vms] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addDynamicVMSLESTime(double t0)
{
  double t = getTime() - t0;
  counter[dvms]++;
  data[dvms] += t;
  return t;
}

//------------------------------------------------------------------------------
	   
double Timer::addH2SetupTime(double t0) 
{ 
  
  double t = getTime() - t0;

  counter[h2Assembly]++;
  data[h2Assembly] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addPrecSetupTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[fluidPrecSetup]++;
  data[fluidPrecSetup] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addKspTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[fluidKsp]++;
  data[fluidKsp] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addMeshMetricsTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[meshMetrics]++;
  data[meshMetrics] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addEmbeddedForceTime(double t0)
{

  double t = getTime() - t0;

  counter[embedforce]++;
  data[embedforce] += t;
  data[eulerFSI] += t;

  return t;

}

//------------------------------------------------------------------------------
                                                                                                                             
double Timer::addStructUpdTime(double t0)
{

  double t = getTime() - t0;

  counter[structUpd]++;
  data[structUpd] += t;

  return t;

}

//------------------------------------------------------------------------------
                                                                                                                             
double Timer::addICInterpTime(double t0)
{

  double t = getTime() - t0;

  counter[icInterp]++;
  data[icInterp] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addFluidSolutionTime(double t0)
{

  double t = getTime() - t0;

  counter[fluid]++;
  data[fluid] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addMeshSolutionTime(double t0)
{

  double t = getTime() - t0;

  counter[mesh]++;
  data[mesh] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addMeshAssemblyTime(double t0)
{

  double t = getTime() - t0;

  counter[meshAssembly]++;
  data[meshAssembly] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addMeshPrecSetupTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[meshPrecSetup]++;
  data[meshPrecSetup] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addMeshKspTime(double t0)
{

  double t = getTime() - t0;

  counter[meshKsp]++;
  data[meshKsp] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::removeForceAndDispComm(double t0)
{

  double t = getTime() - t0;

  data[mesh] -= t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addPodConstrTime(double t0)
{

  double t = getTime() - t0;

  counter[podConstr]++;
  data[podConstr] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addSnapsLinSolvTime(double t0)
{

  double t = getTime() - t0;

  counter[snapsLinSolv]++;
  data[snapsLinSolv] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addPadeReconstrTime(double t0)
{

  double t = getTime() - t0;

  counter[padeReconstr]++;
  data[padeReconstr] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addCorrelMatrixTime(double t0)
{

  double t = getTime() - t0;

  counter[correlMatrix]++;
  data[correlMatrix] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addEigSolvTime(double t0)
{

  double t = getTime() - t0;

  counter[eigSolv]++;
  data[eigSolv] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addAJTime(double t0)
{

  double t = getTime() - t0;

  counter[aj]++;
  data[aj] += t;

  return t;

}


//------------------------------------------------------------------------------

double Timer::addJacEvaluateTime(double t0)
{

  double t = getTime() - t0;

  counter[jacEvaluate]++;
  data[jacEvaluate] += t;

  return t;

}


//------------------------------------------------------------------------------

double Timer::addJacApplyTime(double t0)
{

  double t = getTime() - t0;

  counter[jacApply]++;
  data[jacApply] += t;

  return t;

}


//------------------------------------------------------------------------------

double Timer::addResidualTime(double t0)
{

  double t = getTime() - t0;

  counter[residual]++;
  data[residual] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addRestrictionTime(double t0)
{

  double t = getTime() - t0;

  counter[restriction]++;
  data[restriction] += t;

  return t;

}
//------------------------------------------------------------------------------

double Timer::addSolutionIncrementTime(double t0)
{

  double t = getTime() - t0;

  counter[solutionIncrement]++;
  data[solutionIncrement] += t;

  return t;

}
//------------------------------------------------------------------------------

double Timer::addLinearSystemFormTime(double t0)
{

  double t = getTime() - t0;

  counter[linearSystemForm]++;
  data[linearSystemForm] += t;

  return t;

}
//------------------------------------------------------------------------------

double Timer::addLinearSystemSolveTime(double t0)
{

  double t = getTime() - t0;

  counter[linearSystemSolve]++;
  data[linearSystemSolve] += t;

  return t;

}
//------------------------------------------------------------------------------

double Timer::addCheckConvergenceTime(double t0)
{

  double t = getTime() - t0;

  counter[checkConvergence]++;
  data[checkConvergence] += t;

  return t;

}
//------------------------------------------------------------------------------

double Timer::addGramSchmidtTime(double t0)
{

  double t = getTime() - t0;

  counter[gramSchmidt]++;
  data[gramSchmidt] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addRomSolTime(double t0)
{

  double t = getTime() - t0;

  counter[romSol]++;
  data[romSol] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addRomConstrTime(double t0)
{

  double t = getTime() - t0;

  counter[romConstr]++;
  data[romConstr] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addRomTimeIntegTime(double t0)
{

  double t = getTime() - t0;

  counter[romTimeInteg]++;
  data[romTimeInteg] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addLocalComTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[localCom]++;
  data[localCom] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addGlobalComTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[globalCom]++;
  data[globalCom] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addInterComTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[interCom]++;
  data[interCom] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addBinaryReadTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[binread]++;
  data[binread] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addBinaryWriteTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[binwrite]++;
  data[binwrite] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addLevelSetSolutionTime(double t0)
{

  double t = getTime() - t0;

  counter[levelSet]++;
  data[levelSet] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addLSNodalWeightsAndGradTime(double t0)  {

  double t = getTime() - t0;

  counter[lsNodalWeightsAndGrad]++;
  data[lsNodalWeightsAndGrad] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addLSFiniteVolumeTermTime(double t0)
{

  double t = getTime() - t0;

  counter[lsFvTerm]++;
  data[lsFvTerm] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addLSKspTime(double t0)
{

  double t = getTime() - t0;

  counter[lsKsp]++;
  data[lsKsp] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addLSPrecSetupTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[lsPrecSetup]++;
  data[lsPrecSetup] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addLSFiniteVolumeJacTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[lsJac]++;
  data[lsJac] += t; 

  return t;

}

double Timer::addLSReinitializationTime(double t0) 
{ 

  double t = getTime() - t0;

  counter[lsreinitialization]++;
  data[lsreinitialization] += t; 

  return t;

}

//------------------------------------------------------------------------------

double Timer::addWaitAndReceiveDisp(double t0)
{

  double t = getTime() - t0;

  counter[waitrec]++;
  data[waitrec] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addIntersectionTime(double t0)
{
  double t = getTime() - t0;
  
  counter[intersect]++;
  data[intersect] += t;
  data[eulerFSI] += t;

  return t;
}

//------------------------------------------------------------------------------

double Timer::addEmbedPhaseChangeTime(double t0)
{
  double t = getTime() - t0;

  counter[embedPhaseChange]++;
  data[embedPhaseChange] += t;
  data[eulerFSI] += t;

  return t;
}

//------------------------------------------------------------------------------

double Timer::addRMAComTime(double t0)
{
  double t = getTime() - t0;

//  counter[rmaCom]++;
  data[rmaCom] += t;

  return t;
}

//------------------------------------------------------------------------------

double Timer::removeIntersAndPhaseChange(double t0)  //removed from "Fluid Solution"
{

  double t = getTime() - t0;

  data[fluid] -= t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addWallDistanceTime(double t0)
{

  double t = getTime() - t0;

  counter[walldistance]++;
  data[walldistance] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addReadSnapshotFileTime(double t0) {

  double t = getTime() - t0;
  data[readSnapshotFile] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addClusteringTime(double t0) {

  double t = getTime() - t0;
  data[clustering] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addPODTime(double t0) {

  double t = getTime() - t0;
  data[pod] += t;

  return t;

}

//------------------------------------------------------------------------------

double Timer::addDistCalcsPreproTime(double t0) {

  double t = getTime() - t0;
  data[distCalcsPrepro] += t;

  return t;

}

//------------------------------------------------------------------------------
                
double Timer::addExactUpdatesPreproTime(double t0) {
                
  double t = getTime() - t0;
  data[exactUpdatesPrepro] += t;

  return t;

} 

//------------------------------------------------------------------------------
                
double Timer::addProjErrorTime(double t0) {
                
  double t = getTime() - t0;
  data[projError] += t;

  return t;

} 

//------------------------------------------------------------------------------

double Timer::addMDSTime(double t0) {
  
  double t = getTime() - t0;
  data[mds] += t;

  return t;

}

//------------------------------------------------------------------------------
                
double Timer::addApproxMetricPreproTime(double t0) {
                
  double t = getTime() - t0;
  data[approxMetricPrepro] += t;

  return t;

} 

//------------------------------------------------------------------------------
                
double Timer::addSurfaceMeshConstructionTime(double t0) {
                
  double t = getTime() - t0;
  data[surfaceMeshConstruction] += t;

  return t;

} 

//------------------------------------------------------------------------------
                
double Timer::addSurfaceOutputTime(double t0) {
                
  double t = getTime() - t0;
  data[surfaceOutput] += t;

  return t;

} 


//------------------------------------------------------------------------------
                
double Timer::addSampledMeshConstructionTime(double t0) {
                
  double t = getTime() - t0;
  data[sampledMeshConstruction] += t;

  return t;

} 

//------------------------------------------------------------------------------
                
double Timer::addSampledOutputTime(double t0) {
                
  double t = getTime() - t0;
  data[sampledOutput] += t;

  return t;

} 

//------------------------------------------------------------------------------
                
double Timer::addPseudoInvTime(double t0) {
                
  double t = getTime() - t0;
  data[pseudoInv] += t;

  return t;

} 

//------------------------------------------------------------------------------
                
double Timer::addTotalGappyOfflineTime(double t0) {
                
  double t = getTime() - t0;
  data[gappyOffline] += t;

  return t;

} 

//------------------------------------------------------------------------------

double Timer::addTotalOfflineTime(double t0) {

  double t = getTime() - t0;
  data[romOffline] += t;

  return t;

}

//------------------------------------------------------------------------------
// note: the timings of both fluid and mesh parts contain their communication
void Timer::print(Timer *str, FILE *fp)
{

  if (!com) return;

  double *tmin = reinterpret_cast<double *>(alloca(numTimings * sizeof(double)));
  double *tmax = reinterpret_cast<double *>(alloca(numTimings * sizeof(double)));
  double *tavg = reinterpret_cast<double *>(alloca(numTimings * sizeof(double)));

  if (str) {
    counter[interCom] += str->counter[interCom];
    data[interCom] += str->data[interCom];
  }

  data[total] = data[setup] + data[run];


  data[comm] = data[localCom] + data[globalCom] + data[rmaCom] + data[interCom];
  data[io] = data[binread] + data[binwrite];
  
  if (ioData->problem.alltype == ProblemData::_POD_CONSTRUCTION_)
    data[podConstr] -= data[io];

  int i;

  for (i=0; i<numTimings ; ++i) {
    tmin[i] = data[i];
    tmax[i] = data[i];
    tavg[i] = data[i];
  }

  com->globalMin(numTimings, tmin);
  com->globalMax(numTimings, tmax);
  com->globalSum(numTimings, tavg);

  int numCPU = com->size();

  for (i=0; i<numTimings ; ++i)
    tavg[i] /= numCPU;

  com->fprintf(fp, "\n");
  com->fprintf(fp, "----------------------------------------------------------------------\n");
  com->fprintf(fp, "Elapsed Time Report (s)       :        Min        Max        Avg   # Calls\n");
  com->fprintf(fp, "\n");
  com->fprintf(fp, "Problem Setup                 : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[setup], tmax[setup], tavg[setup], 
	       counter[setup]);
  com->fprintf(fp, "\n");
  com->fprintf(fp, "Fluid Solution                : %10.2f %10.2f %10.2f         -\n", 
	       tmin[fluid], tmax[fluid], tavg[fluid]);
  com->fprintf(fp, "  Time Steps                  : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[timeStep], tmax[timeStep], tavg[timeStep], 
	       counter[timeStep]);
  com->fprintf(fp, "  Nodal Weights and Gradients : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[nodalGrad], tmax[nodalGrad], tavg[nodalGrad], 
	       counter[nodalGrad]);
  com->fprintf(fp, "  FV Fluxes                   : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[fvTerm], tmax[fvTerm], tavg[fvTerm], 
	       counter[fvTerm]);
  com->fprintf(fp, "  FE Fluxes                   : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[feTerm], tmax[feTerm], tavg[feTerm], 
	       counter[feTerm]);
  com->fprintf(fp, "  FV Jacobian                 : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[fvJac], tmax[fvJac], tavg[fvJac], 
	       counter[fvJac]);
  com->fprintf(fp, "  FE Jacobian                 : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[feJac], tmax[feJac], tavg[feJac], 
	       counter[feJac]);
  com->fprintf(fp, "  VMS-LES Modeling            : %10.2f %10.2f %10.2f %9d\n",
               tmin[vms], tmax[vms], tavg[vms],
               counter[vms]);
  com->fprintf(fp, "  Dynamic VMS-LES Modeling    : %10.2f %10.2f %10.2f %9d\n",
               tmin[dvms], tmax[dvms], tavg[dvms],
               counter[dvms]);
  com->fprintf(fp, "  H2 Matrix Assembly          : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[h2Assembly], tmax[h2Assembly], tavg[h2Assembly], 
	       counter[h2Assembly]);
  com->fprintf(fp, "  Preconditioner Setup        : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[fluidPrecSetup], tmax[fluidPrecSetup], tavg[fluidPrecSetup], 
	       counter[fluidPrecSetup]);
  com->fprintf(fp, "  Linear Solver               : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[fluidKsp], tmax[fluidKsp], tavg[fluidKsp], 
	       counter[fluidKsp]);
  com->fprintf(fp, "  Mesh Metrics Update         : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[meshMetrics], tmax[meshMetrics], tavg[meshMetrics], 
	       counter[meshMetrics]);
  if (ioData->problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_AEROELASTIC_)  {
    com->fprintf(fp, "  Structural Update           : %10.2f %10.2f %10.2f %9d\n",
               tmin[structUpd], tmax[structUpd], tavg[structUpd], counter[structUpd]);
  }
  if (ioData->input.multiSolutionsParams[0] != 0)  {
    com->fprintf(fp, "  Interpolated IC             : %10.2f %10.2f %10.2f %9d\n",
               tmin[icInterp], tmax[icInterp], tavg[icInterp], counter[icInterp]);
  }

  if(ioData->problem.framework == ProblemData::EMBEDDED || ioData->problem.framework == ProblemData::EMBEDDEDALE) {
    if (ioData->eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
      if (ioData->eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
          ioData->eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES) {
        com->fprintf(fp, "  Wall Distance Computation   : %10.2f %10.2f %10.2f %9d\n",
                   tmin[walldistance], tmax[walldistance], tavg[walldistance],
                   counter[walldistance]);
      }
    }
  }
  com->fprintf(fp, "\n");

  if ((ioData->problem.alltype == ProblemData::_STEADY_NONLINEAR_ROM_ ) || 
      (ioData->problem.alltype == ProblemData::_UNSTEADY_NONLINEAR_ROM_ ) ||
      (ioData->problem.alltype == ProblemData::_ACC_UNSTEADY_NONLINEAR_ROM_ ) ||
      (ioData->problem.alltype == ProblemData::_FORCED_NONLINEAR_ROM_ ))  {
    com->fprintf(fp, "  Residual evaluation         : %10.2f %10.2f %10.2f %9d\n",
              tmin[residual], tmax[residual], tavg[residual], counter[residual]);
    com->fprintf(fp, "  Action of Jacobian on ROB   : %10.2f %10.2f %10.2f %9d\n",
              tmin[aj], tmax[aj], tavg[aj], counter[aj]);
    com->fprintf(fp, "    Jacobian evaluation       : %10.2f %10.2f %10.2f %9d\n",
              tmin[jacEvaluate], tmax[jacEvaluate], tavg[jacEvaluate], counter[jacEvaluate]);
    com->fprintf(fp, "    Jacobian vector product   : %10.2f %10.2f %10.2f %9d\n",
              tmin[jacApply], tmax[jacApply], tavg[jacApply], counter[jacApply]);
    com->fprintf(fp, "  Solution increment          : %10.2f %10.2f %10.2f %9d\n",
              tmin[solutionIncrement], tmax[solutionIncrement], tavg[solutionIncrement], counter[solutionIncrement]);
    com->fprintf(fp, "  Linear system form          : %10.2f %10.2f %10.2f %9d\n",
              tmin[linearSystemForm], tmax[linearSystemForm], tavg[linearSystemForm], counter[linearSystemForm]);
    com->fprintf(fp, "  Linear system solve         : %10.2f %10.2f %10.2f %9d\n",
              tmin[linearSystemSolve], tmax[linearSystemSolve], tavg[linearSystemSolve], counter[linearSystemSolve]);
    com->fprintf(fp, "  Restriction                 : %10.2f %10.2f %10.2f %9d\n",
              tmin[restriction], tmax[restriction], tavg[restriction], counter[restriction]);
    com->fprintf(fp, "  Check convergence           : %10.2f %10.2f %10.2f %9d\n",
              tmin[checkConvergence], tmax[checkConvergence], tavg[checkConvergence], counter[checkConvergence]);
	}

  // Output Mesh solution time (except for Euler FSI)
  if(ioData->problem.framework != ProblemData::EMBEDDED) {
    com->fprintf(fp, "Mesh Solution                 : %10.2f %10.2f %10.2f         -\n", 
                 tmin[mesh], tmax[mesh], tavg[mesh]);
    com->fprintf(fp, "  K Matrix Assembly           : %10.2f %10.2f %10.2f %9d\n", 
                 tmin[meshAssembly], tmax[meshAssembly], tavg[meshAssembly], 
                 counter[meshAssembly]);
    com->fprintf(fp, "  Preconditioner Setup        : %10.2f %10.2f %10.2f %9d\n", 
                 tmin[meshPrecSetup], tmax[meshPrecSetup], tavg[meshPrecSetup], 
	         counter[meshPrecSetup]);
    com->fprintf(fp, "  Linear Solver               : %10.2f %10.2f %10.2f %9d\n", 
                 tmin[meshKsp], tmax[meshKsp], tavg[meshKsp], 
	         counter[meshKsp]);
    com->fprintf(fp, "\n");
  }

  // Output Level-Set Timers
  if (ioData->eqs.numPhase > 1) {
    com->fprintf(fp, "LevelSet Solution             : %10.2f %10.2f %10.2f         -\n", tmin[levelSet], tmax[levelSet], tavg[levelSet]);

    com->fprintf(fp, "  Nodal Weights and Grad      : %10.2f %10.2f %10.2f %9d\n",
               tmin[lsNodalWeightsAndGrad], tmax[lsNodalWeightsAndGrad], tavg[lsNodalWeightsAndGrad],
               counter[lsNodalWeightsAndGrad]);
    com->fprintf(fp, "  FV Fluxes                   : %10.2f %10.2f %10.2f %9d\n",
               tmin[lsFvTerm], tmax[lsFvTerm], tavg[lsFvTerm], counter[lsFvTerm]);
    com->fprintf(fp, "  FV Jacobian                 : %10.2f %10.2f %10.2f %9d\n", 
		 tmin[lsJac], tmax[lsJac], tavg[lsJac], 
		 counter[lsJac]);
    com->fprintf(fp, "  Reinitialization            : %10.2f %10.2f %10.2f %9d\n", 
		 tmin[lsreinitialization], tmax[lsreinitialization],
                 tavg[lsreinitialization], 
		 counter[lsreinitialization]);
    com->fprintf(fp, "  Preconditioner Setup        : %10.2f %10.2f %10.2f %9d\n", 
		 tmin[lsPrecSetup], tmax[lsPrecSetup], tavg[lsPrecSetup], 
		 counter[lsPrecSetup]);
    com->fprintf(fp, "  Linear Solver               : %10.2f %10.2f %10.2f %9d\n", tmin[lsKsp], tmax[lsKsp], tavg[lsKsp], counter[lsKsp]);
    com->fprintf(fp, "\n");
  }


  // Output POD Timers
  if (ioData->problem.alltype == ProblemData::_NONLINEAR_ROM_OFFLINE_) {
    com->fprintf(fp, "Offline ROM Precomputations   : %10.2f %10.2f %10.2f         -\n",
               tmin[romOffline], tmax[romOffline], tavg[romOffline]);
    com->fprintf(fp, "  Read State Snapshot Files   : %10.2f %10.2f %10.2f         -\n",
               tmin[readSnapshotFile], tmax[readSnapshotFile], tavg[readSnapshotFile]);
    com->fprintf(fp, "  K-Means Clustering          : %10.2f %10.2f %10.2f         -\n",
               tmin[clustering], tmax[clustering], tavg[clustering]);
    com->fprintf(fp, "  SVDs                        : %10.2f %10.2f %10.2f         -\n",
               tmin[pod], tmax[pod], tavg[pod]);
    com->fprintf(fp, "  Fast Distance Calcs Prepro  : %10.2f %10.2f %10.2f         -\n",
               tmin[distCalcsPrepro], tmax[distCalcsPrepro], tavg[distCalcsPrepro]);
    com->fprintf(fp, "  Fast Exact Updates Prepro   : %10.2f %10.2f %10.2f         -\n",
               tmin[exactUpdatesPrepro], tmax[exactUpdatesPrepro], tavg[exactUpdatesPrepro]);
    com->fprintf(fp, "  Relative Projection Error   : %10.2f %10.2f %10.2f         -\n",
               tmin[projError], tmax[projError], tavg[projError]);
    com->fprintf(fp, "  Multi-Dimensional Scaling   : %10.2f %10.2f %10.2f         -\n",
               tmin[mds], tmax[mds], tavg[mds]);
    com->fprintf(fp, "  Offline Gappy Prepro        : %10.2f %10.2f %10.2f         -\n",
               tmin[gappyOffline], tmax[gappyOffline], tavg[gappyOffline]);
    com->fprintf(fp, "    Approx Metric Prepro      : %10.2f %10.2f %10.2f         -\n",
               tmin[approxMetricPrepro], tmax[approxMetricPrepro], tavg[approxMetricPrepro]);
    com->fprintf(fp, "    Sampled Mesh Construction : %10.2f %10.2f %10.2f         -\n",
               tmin[sampledMeshConstruction], tmax[sampledMeshConstruction], tavg[sampledMeshConstruction]);
    com->fprintf(fp, "    Sampled Quantity Output   : %10.2f %10.2f %10.2f         -\n",
               tmin[sampledOutput], tmax[sampledOutput], tavg[sampledOutput]);
    com->fprintf(fp, "    Pseudo-Inverse            : %10.2f %10.2f %10.2f         -\n",
               tmin[pseudoInv], tmax[pseudoInv], tavg[pseudoInv]);
    com->fprintf(fp, "    Surface Mesh Construction : %10.2f %10.2f %10.2f         -\n",
               tmin[surfaceMeshConstruction], tmax[surfaceMeshConstruction], tavg[surfaceMeshConstruction]);
    com->fprintf(fp, "    Surface Quantity Output   : %10.2f %10.2f %10.2f         -\n",
               tmin[surfaceOutput], tmax[surfaceOutput], tavg[surfaceOutput]);
    com->fprintf(fp, "\n");
  }
  else if (ioData->problem.alltype == ProblemData::_POD_CONSTRUCTION_) {
    com->fprintf(fp, "POD Basis Construction        : %10.2f %10.2f %10.2f         -\n",
               tmin[podConstr], tmax[podConstr], tavg[podConstr]);
    com->fprintf(fp, "  Snapshot Linear Solver      : %10.2f %10.2f %10.2f %9d\n",
               tmin[snapsLinSolv], tmax[snapsLinSolv], tavg[snapsLinSolv], counter[snapsLinSolv]);
    if (ioData->linearizedData.padeReconst == LinearizedData::TRUE) {
    com->fprintf(fp, "  Pade Reconstruction         : %10.2f %10.2f %10.2f %9d\n",
               tmin[padeReconstr], tmax[padeReconstr], tavg[padeReconstr], counter[padeReconstr]);
    }
    com->fprintf(fp, "  Correlation Matrix          : %10.2f %10.2f %10.2f %9d\n",
               tmin[correlMatrix], tmax[correlMatrix], tavg[correlMatrix], counter[correlMatrix]);
    com->fprintf(fp, "  SVD Solver                  : %10.2f %10.2f %10.2f %9d\n",
               tmin[eigSolv], tmax[eigSolv], tavg[eigSolv], counter[eigSolv]);
    com->fprintf(fp, "  Gram-Schmidt                : %10.2f %10.2f %10.2f %9d\n",
               tmin[gramSchmidt], tmax[gramSchmidt], tavg[gramSchmidt], counter[gramSchmidt]);
    com->fprintf(fp, "\n");
  }
  else if (ioData->problem.alltype == ProblemData::_ROM_AEROELASTIC_)  {
    com->fprintf(fp, "ROM Solution                  : %10.2f %10.2f %10.2f         -\n",
              tmin[romSol], tmax[romSol], tavg[romSol]);
    com->fprintf(fp, "  ROM Construction            : %10.2f %10.2f %10.2f %9d\n",
               tmin[romConstr], tmax[romConstr], tavg[romConstr], counter[romConstr]);
    com->fprintf(fp, "  ROM Time Integration        : %10.2f %10.2f %10.2f %9d\n",
               tmin[romTimeInteg], tmax[romTimeInteg], tavg[romTimeInteg], counter[romTimeInteg]);
    com->fprintf(fp, "\n");
  }

  if(ioData->problem.framework == ProblemData::EMBEDDED) {
    com->fprintf(fp, "Eulerian FSI                  : %10.2f %10.2f %10.2f         -\n",
                 tmin[eulerFSI], tmax[eulerFSI], tavg[eulerFSI]);
    com->fprintf(fp, "  F-S Intersections           : %10.2f %10.2f %10.2f %9d\n", 
  	         tmin[intersect], tmax[intersect], tavg[intersect], 
  	         counter[intersect]);
    com->fprintf(fp, "  Force calculation           : %10.2f %10.2f %10.2f %9d\n",
                 tmin[embedforce], tmax[embedforce], tavg[embedforce],
                 counter[embedforce]);
    com->fprintf(fp,"\n");
  }

  com->fprintf(fp, "Communication/Synchronization : %10.2f %10.2f %10.2f         -\n", 
	       tmin[comm], tmax[comm], tavg[comm]);
  com->fprintf(fp, "  Local                       : %10.2f %10.2f %10.2f         -\n",
	       tmin[localCom], tmax[localCom], tavg[localCom]);
  com->fprintf(fp, "  Global                      : %10.2f %10.2f %10.2f         -\n", 
	       tmin[globalCom], tmax[globalCom], tavg[globalCom]);
  com->fprintf(fp, "  RMA                         : %10.2f %10.2f %10.2f         -\n", 
	       tmin[rmaCom], tmax[rmaCom], tavg[rmaCom]);
  com->fprintf(fp, "  Inter                       : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[interCom], tmax[interCom], tavg[interCom], 
	       counter[interCom]);
  com->fprintf(fp, "\n");
  com->fprintf(fp, "I/O                           : %10.2f %10.2f %10.2f         -\n", 
	       tmin[io], tmax[io], tavg[io]);
  com->fprintf(fp, "  Binary Read                 : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[binread], tmax[binread], tavg[binread], 
	       counter[binread]);
  com->fprintf(fp, "  Binary Write                : %10.2f %10.2f %10.2f %9d\n", 
	       tmin[binwrite], tmax[binwrite], tavg[binwrite], 
	       counter[binwrite]);

  com->fprintf(fp, "\n");
  com->fprintf(fp, "Total Simulation              : %10.2f %10.2f %10.2f         -\n", 
	       tmin[total], tmax[total], tavg[total]);
  com->fprintf(fp, "----------------------------------------------------------------------\n");

}

//------------------------------------------------------------------------------
