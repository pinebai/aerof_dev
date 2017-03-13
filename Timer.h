#ifndef _TIMER_H_
#define _TIMER_H_

#include <IoData.h>

class Communicator;
class IoData;
//------------------------------------------------------------------------------

class Timer {

  enum TimerIndex {
		setup, run, total, fluid, nodalWeights, nodalGrad, fvTerm, feTerm, fvJac,
		feJac, vms, dvms, h2Assembly, fluidPrecSetup, fluidKsp, meshMetrics,
		structUpd, mesh, meshAssembly, meshPrecSetup, meshKsp, podConstr,
		snapsLinSolv, padeReconstr, correlMatrix, eigSolv, aj, jacEvaluate, jacApply, residual, restriction,
		solutionIncrement, linearSystemForm, linearSystemSolve, checkConvergence, gramSchmidt, romSol,
		romConstr, romTimeInteg, comm, localCom, globalCom, interCom, rmaCom, io,
		binread, binwrite, levelSet, lsNodalWeightsAndGrad, lsFvTerm,
		lsKsp,lsPrecSetup,lsJac, waitrec, timeStep, intersect, embedPhaseChange,
		eulerFSI, embedforce, walldistance, lsreinitialization, readSnapshotFile,
		clustering, pod, distCalcsPrepro, exactUpdatesPrepro, projError, mds, 	
		approxMetricPrepro, surfaceMeshConstruction, surfaceOutput, sampledMeshConstruction, sampledOutput, pseudoInv,
                gappyOffline, romOffline, icInterp, NUMTIMINGS
  };

  int numTimings;

  double initialTime;

  int *counter;
  double *data;
  
  IoData *ioData;
  
  Communicator *com;

public:

  Timer(Communicator *);
  ~Timer();

  double getTime();
  double getTimeSyncro();
  double getRunTime();

  void setIoData(IoData &_ioData);
  void setSetupTime();
  void setRunTime();

  double addTimeStepTime(double);
  double addNodalWeightsTime(double);
  double addNodalGradTime(double);
  double addFiniteVolumeTermTime(double);
  double addFiniteElementTermTime(double);
  double addFiniteVolumeJacTime(double);
  double addFiniteElementJacTime(double);
  double addVMSLESTime(double);
  double addDynamicVMSLESTime(double);
  double addH2SetupTime(double);
  double addPrecSetupTime(double);
  double addKspTime(double);
  double addMeshMetricsTime(double);
  double addEmbeddedForceTime(double);
  double addStructUpdTime(double);
  double addICInterpTime(double);
  double addMeshSolutionTime(double);
  double addMeshAssemblyTime(double);
  double addMeshPrecSetupTime(double);
  double addMeshKspTime(double);
  double removeForceAndDispComm(double);
  double addPodConstrTime(double);
  double addProjectTime(double); //CBM--check
  double addSnapsLinSolvTime(double);
  double addPadeReconstrTime(double);
  double addCorrelMatrixTime(double);
  double addEigSolvTime(double);
  double addAJTime(double);
  double addJacEvaluateTime(double);
  double addJacApplyTime(double);
  double addResidualTime(double);
  double addRestrictionTime(double);
  double addSolutionIncrementTime(double);
  double addLinearSystemFormTime(double);
  double addLinearSystemSolveTime(double);
  double addCheckConvergenceTime(double);
  double addGramSchmidtTime(double);
  double addRomSolTime(double);
  double addRomConstrTime(double);
  double addRomTimeIntegTime(double);
  double addLocalComTime(double);
  double addGlobalComTime(double);
  double addInterComTime(double);
  double addRMAComTime(double);
  double addBinaryReadTime(double);
  double addBinaryWriteTime(double);
  double addFluidSolutionTime(double);
  double addWaitAndReceiveDisp(double);

  // Level-Set Timer Functions
  double addLevelSetSolutionTime(double);
  double addLSNodalWeightsAndGradTime(double);
  double addLSFiniteVolumeTermTime(double);
  double addLSKspTime(double);
  double addLSPrecSetupTime(double);
  double addLSFiniteVolumeJacTime(double);
  double addLSReinitializationTime(double);

  // Embedded FSI Timer Functions
  double addIntersectionTime(double);
  double addEmbedPhaseChangeTime(double);
  double removeIntersAndPhaseChange(double);
  double addWallDistanceTime(double);

  // Nonlinear ROM Offline Timer Functions
  double addReadSnapshotFileTime(double);
  double addClusteringTime(double);
  double addPODTime(double);
  double addDistCalcsPreproTime(double);
  double addExactUpdatesPreproTime(double);
  double addProjErrorTime(double);
  double addMDSTime(double);
  double addApproxMetricPreproTime(double);
  double addSampledMeshConstructionTime(double);
  double addSampledOutputTime(double);
  double addSurfaceMeshConstructionTime(double);
  double addSurfaceOutputTime(double);
  double addPseudoInvTime(double);
  double addTotalGappyOfflineTime(double);
  double addTotalOfflineTime(double);

  void print(Timer *, FILE * = stdout);

};

//------------------------------------------------------------------------------

#endif
