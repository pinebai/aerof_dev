#ifndef _TS_DESC_H_
#define _TS_DESC_H_

#include <IoData.h>
#include <TsInput.h>
#include <TsOutput.h>
#include <TsRestart.h>
#include <TsParameters.h>
//#include <Domain.h>
#include <DistVector.h>
#include <MultiGridKernel.h>
#include <ErrorHandler.h>

class RefVal;
class VarFcn;
class GeoSource;
class DistGeoState;
class MeshMotionHandler;
class MeshMotionSolver;
class HeatTransferHandler;
class MemoryPool;
class Timer;

template<int dim> class DistBcData;
template<int dim> class DistTimeState;
template<int dim> class SpaceOperator;
template<int dim> class PostOperator;
template<int dim> class DistExactRiemannSolver;

//------------------------------------------------------------------------------
/** Base class for time integration descriptors */
template<int dim>
class TsDesc {

public:

  typedef DistSVec<double,dim> SolVecType;
  typedef DistSVec<double,3> PosVecType;
  typedef DistVec<double> VolVecType;
  typedef DistSVec<double,1> PhiVecType;

  TsParameters *data;
  virtual double computeResidualNorm(DistSVec<double,dim>&);
  virtual void monitorInitialState(int, DistSVec<double,dim> &);

  Timer *timer;

  void performPostProForState(DistSVec<double,dim> &, int);

protected:
  DistExactRiemannSolver<dim> *riemann1;

  PosVecType *X; //vector of positions
  VolVecType *A; //vector of cell volumes
  PosVecType *Xs;

  bool* problemType;
  int numPhase;
  TsData::Clipping clippingType;
  BcsWallData::Integration wallType;
  BcsWallData::Reconstruction wallRecType;
  TsData::TimeStepCalculation timeStepCalculation;

//  TsParameters *data;
  TsInput *input;
  TsOutput<dim> *output;
  TsRestart *restart;

  DistSVec<double,dim> *Uic;
  DistSVec<double,dim> *V;
  DistSVec<double,dim> *F;
  DistSVec<double,dim> *R;
  DistSVec<double,dim> *Rinlet;
  DistSVec<double,dim> *Rreal;
  RefVal *refVal;
  VarFcn *varFcn;

  DistVec<int> fluidIdDummy;
  DistTimeState<dim> *timeState;
  DistBcData<dim> *bcData;
  DistGeoState *geoState;
  SpaceOperator<dim> *spaceOp;
  PostOperator<dim> *postOp;

  MeshMotionHandler* mmh;
  MeshMotionSolver* mems;
  HeatTransferHandler* hth;

  Domain *domain;

//  Timer *timer;
  Communicator *com;

  ErrorHandler *errorHandler;

// Included (MB)
  int fixSol;
  int iForce;
  int iTotal;

  double forceNorm;
  double *forceNorms;

  bool failSafeFlag;

  //Modified Ghidaglia scheme with 'external state' estimated using the Hagstrom b.c.
  bool modifiedGhidaglia;

  bool isMultigridTsDesc;

protected:

//  void monitorInitialState(int, DistSVec<double,dim> &);

  void moveMesh(IoData &ioData, GeoSource &geoSource);  // YC

// Included (MB)
  bool monitorForceConvergence(IoData &, int, DistSVec<double,dim> &);
  bool monitorAvgForceConvergence(IoData &, int, DistSVec<double,dim> &);

public:

  Communicator* gComm(){return this->com;};//TODO delete line HACK

  TsDesc(IoData &, GeoSource &, Domain *);
  virtual ~TsDesc();

  void printf(int, const char *, ...);
  void fprintf(FILE *, const char *, ...);
  VarFcn *createVarFcn(IoData &);
  DistBcData<dim> *createBcData(IoData &);
  MeshMotionHandler *createMeshMotionHandler(IoData &, GeoSource &, MemoryPool *);

  HeatTransferHandler* createHeatTransferHandler(IoData&, GeoSource&);

  SpaceOperator<dim>* getSpaceOperator() { return spaceOp; }
  Communicator* getCommunicator() { return com; }

  virtual bool monitorConvergence(int, DistSVec<double,dim> &);

  double recomputeResidual(DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluateFluxAtMultipleSolutions(IoData &iod, char* best_soln);
  void setRestartIterationAndTime(int it, double t) { restart->etime = t; restart->iteration = it;}
  void formInterpolationWeights(IoData&);
  virtual void formInterpolatedInitialCondition(DistSVec<double,dim> *U, IoData &iod);
  virtual void setInterpWeightsForMultiIC(std::vector<double> ) {}
  virtual void setupTimeStepping(DistSVec<double,dim> *, IoData &);
  virtual double computeTimeStep(int, double *, DistSVec<double,dim> &, double);
  virtual double computeTimeStep(int a, double *b, DistSVec<double,dim> &c){ return computeTimeStep(a,b,c,-2); }
  virtual void cmdCom(bool *);
  virtual void getNumParam(int &, int &, double &);
  virtual void sendNumParam(int);
  virtual void getRelResidual(double &);
  virtual double computePositionVector(bool *, int, double, DistSVec<double,dim> &);
//  virtual double computePositionSensitivityVector(bool *, int, double);
  virtual void setMeshSensitivitySolverPositionVector();
  void negotiate();
  void sendForceSensitivity(DistSVec<double,3> *, bool applyScale = true);
  void interpolatePositionVector(double, double);
  void computeMeshMetrics(int it = -1);
  virtual void updateStateVectors(DistSVec<double,dim> &, int = 0);
  bool checkForLastIteration(int, double, double, DistSVec<double,dim> &); //KW: not used?

  void setFailSafe(bool flag){ failSafeFlag = flag; }

  DistSVec<double,dim>& getCurrentResidual() { return *R; }

// Modified (MB)
  virtual bool checkForLastIteration(IoData &, int, double, double, DistSVec<double,dim> &);

  virtual void setupOutputToDisk(IoData &, bool *, int, double,
                                 DistSVec<double,dim> &);
  virtual void outputToDisk(IoData &, bool*, int, int, int, double, double,
                            DistSVec<double,dim> &);
  virtual void writeStateRomToDisk(int it, double cpu) {};
  virtual void writeErrorToDisk(int it, double cpu) {};

  virtual void outputForces(IoData &, bool*, int, int, int, double, double,
                            DistSVec<double,dim> &);

  virtual void outputPositionVectorToDisk(DistSVec<double,dim> &U);
//  void outputPositionSensitivityVectorToDisk(DistSVec<double,dim> &dUds);
  virtual void resetOutputToStructure(DistSVec<double,dim> &);
  virtual void updateOutputToStructure(double, double, DistSVec<double,dim> &);

  virtual bool IncreasePressure(int it, double dt, double t, DistSVec<double,dim> &U){return true;}
  virtual int solveNonLinearSystem(DistSVec<double,dim> &U, int) { return 0; }
  virtual int checkSolution(DistSVec<double,dim> &);

  int getInitialIteration() const { return restart->iteration; }
  double getInitialTime() const { return restart->etime; }
  DistInfo &getVecInfo() const { return domain->getNodeDistInfo(); }
  DistInfo &getInletVecInfo() const {return domain->getInletNodeDistInfo(); }
  int getNumPhase() { return numPhase; }
  int structureSubcycling() {return (mmh ? mmh->structureSubcycling() : 0);}
  virtual bool willNotSolve(double dts, double t) {return false;}

// Included (MB)
  virtual void fixSolution(DistSVec<double,dim> &, DistSVec<double,dim> &);

  virtual void setCurrentTime(double t,DistSVec<double,dim>& U) { }
  virtual void setCurrentTimeStep(double dt) { }

  virtual void setFluidSubcycling(bool inSub) { }

  void updateGhostFluid(DistSVec<double,dim> &, Vec3D&, double);

  void updateFarfieldCoeffs(double dt);
  void updateBoundaryExternalState();
  void initializeFarfieldCoeffs();

  void printNodalDebug(int globNodeId, int identifier, DistSVec<double,dim> *U, DistVec<int> *Id=0, DistVec<int> *Id0=0);

  void computeDistanceToWall(IoData &ioData);

  TsParameters* getTsParams() {return data;}
  ErrorHandler* getErrorHandler() {return errorHandler;}
  void computeConvergenceInformation(IoData &ioData, const char* file, DistSVec<double,dim>&);
  void receiveBoundaryPositionSensitivityVector(DistSVec<double,3> &, bool applyScale = false);

  virtual void checkLocalRomStatus(DistSVec<double, dim> &, const int) {}
  virtual void writeBinaryVectorsToDiskRom(bool, int, int, DistSVec<double,dim> *, DistSVec<double,dim> *) {}  // state, residual
  virtual void incrementNewtonOutputTag() {}
  int *getTimeIt() { return domain->getTimeIt(); }
  int *getNewtonIt() { return domain->getNewtonIt(); }
  int *getNumResidualsOutputCurrentNewtonIt() { return domain->getNumResidualsOutputCurrentNewtonIt(); }
    // Lei Lei, 26 Sep 2016
    virtual void test() {} // a placeholder for derived class to run tests/diagnostics
  void readICFromDisk(char *, int, int, DistSVec<double,dim> &);
  bool outputOnlySpatialResidualBool;
  bool outputOnlySpatialResidual() {return outputOnlySpatialResidualBool;}
  virtual void calculateSpatialResidual(DistSVec<double,dim> &, DistSVec<double,dim> &) {}

  std::vector<double> interpolatedICWeights;
};

//------------------------------------------------------------------------------


#ifdef TEMPLATE_FIX
#include <TsDesc.C>
#endif

#endif
