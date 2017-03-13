#ifndef _MESH_MOTION_SOLVER_H_
#define _MESH_MOTION_SOLVER_H_

#include <IoData.h>
#include <Domain.h>
#include <DistVector.h>
#include <KspPrec.h>
#include <BCApplier.h>
#include <TsParameters.h>
#include <ErrorHandler.h>

class MatchNodeSet;
class CorotSolver;
class Communicator;
class MemoryPool;
class Timer;

template<class Scalar, int dim> class StiffMat;
template<class ProbDesc> class NewtonSolver;
#ifndef _KSPSLVR_TMPL_
#define _KSPSLVR_TMPL_
template<class VecType, class MvpOp, class PrecOp, class IoOp, class ScalarT = double> class KspSolver;
#endif


//------------------------------------------------------------------------------

class MeshMotionSolver {

public:

  MeshMotionSolver() {};
  virtual ~MeshMotionSolver() {}

  virtual int solve(DistSVec<double,3> &, DistSVec<double,3> &) = 0;
  virtual int solveAdjoint(DistSVec<double,3> &, DistSVec<double,3> &) = 0;
  virtual int solveSensitivity(DistSVec<double,3> &, DistSVec<double,3> &) = 0;
  virtual void setup(DistSVec<double,3> &) = 0;
  virtual void applyProjector(DistSVec<double,3> &dX) {}
  virtual void applyProjectorTranspose(DistSVec<double,3> &dX) {}
  virtual void setAdjointFlagOn() {}
  virtual void set_dX0(DistSVec<double,3> &dX) {}
  virtual void computeStiffnessMatrix(DistSVec<double,3> &) {}
  virtual int solveLinearSystem(int, DistSVec<double,3> &, DistSVec<double,3> &) { return 0; }
  virtual void setOperators(DistSVec<double,3> &) {}
  virtual void apply(DistSVec<double,3> &, DistSVec<double,3> &) {}

};

//------------------------------------------------------------------------------

class TetMeshMotionSolver : public MeshMotionSolver {

public:

  typedef DistSVec<double,3> SolVecType;
  typedef DistSVec<double,1> PhiVecType;
  typedef DistVec<double> VolVecType;

protected:

  int maxItsNewton;
  double epsNewton;
  double epsAbsResNewton, epsAbsIncNewton;
  FILE *outputNewton;
  int maxItsLS;
  double contractionLS, sufficDecreaseLS;

  DefoMeshMotionData::Element typeElement;

  DistSVec<double,3> *F0;
  DistSVec<double,3> *dX0;

  CorotSolver *cs;
  StiffMat<double,3> *mvp;
  KspPrec<3> *pc;
  KspSolver<DistSVec<double,3>, StiffMat<double,3>, KspPrec<3>, Communicator> *ksp;
  NewtonSolver<TetMeshMotionSolver> *ns;

  Domain *domain;

  Communicator *com;
  Timer *timer;

  double volStiff;
  bool adjointFlag;
  bool stiffFlag;
  bool sensitivityFlag;

  BCApplier* meshMotionBCs; //HB

public:

  TetMeshMotionSolver(DefoMeshMotionData &, MatchNodeSet **, Domain *, MemoryPool *);

  //Constructor for Embedded ALE case
  TetMeshMotionSolver(Domain *dom) :
	domain(dom),
	adjointFlag(false),
	stiffFlag(false),
	sensitivityFlag(false),
	maxItsNewton(0),
	epsNewton(0.0),
	epsAbsResNewton(0.0),
	epsAbsIncNewton(0.0),
	outputNewton(NULL),
	maxItsLS(0),
	contractionLS(0.0),
	sufficDecreaseLS(0.0),
	F0(NULL),
	dX0(NULL),
	cs(NULL),
	mvp(NULL),
	pc(NULL),
	ksp(NULL),
	ns(NULL),
	com(NULL),
	timer(NULL),
	volStiff(0.0),
	meshMotionBCs(NULL),
	typeElement(DefoMeshMotionData::LINEAR_FE)
	{};
  virtual ~TetMeshMotionSolver();

  virtual int solve(DistSVec<double,3> &, DistSVec<double,3> &);
  virtual int solveAdjoint(DistSVec<double,3> &, DistSVec<double,3> &);
  int solveSensitivity(DistSVec<double,3> &, DistSVec<double,3> &) { return 0; };

  void applyProjector(DistSVec<double,3> &X);
  void applyProjectorTranspose(DistSVec<double,3> &X);
  void apply(DistSVec<double,3> &, DistSVec<double,3> &);
 
  void printf(int, const char *, ...);
  void fprintf(FILE *, const char *, ...);
  virtual void computeFunction(int, DistSVec<double,3> &, DistSVec<double,3> &);
  void recomputeFunction(DistSVec<double,3> &, DistSVec<double,3> &) {}
  virtual void computeStiffnessMatrix(DistSVec<double,3> &);
  double recomputeResidual(DistSVec<double,3> &, DistSVec<double,3> &) { return 0.0; }
  void computeJacobian(int, DistSVec<double,3> &, DistSVec<double,3> &);
  void setOperators(DistSVec<double,3> &);
  void setAdjointFlagOn() { adjointFlag = true; }
  int solveLinearSystem(int, DistSVec<double,3> &, DistSVec<double,3> &);

  int checkSolution(DistSVec<double,3> &X) { return 0; }
  int checkFailSafe(DistSVec<double,3>& X) { return 0; }
  void resetFixesTag() { return;}
  int getMaxItsNewton() const { return maxItsNewton; }
  double getEpsNewton() const { return epsNewton; }
  double getEpsAbsResNewton() const { return epsAbsResNewton; }
  double getEpsAbsIncNewton() const { return epsAbsIncNewton; }
  FILE* getOutputNewton() const { return outputNewton; }
  int getLineSearch() const { return (maxItsLS>0); }
  int getMaxItsLineSearch() const { return maxItsLS; }
  double getContractionLineSearch() const { return contractionLS; }
  double getSufficientDecreaseLineSearch() const { return sufficDecreaseLS; }
  DistInfo &getVecInfo() const { return domain->getNodeDistInfo(); }
  
  virtual void setup(DistSVec<double,3> &X);
  virtual void set_dX0(DistSVec<double,3> &dX);
  TsParameters* getTsParams() { return NULL; }
  ErrorHandler* getErrorHandler() {return NULL; }

  // Included (MB)
  int fixSolution(DistSVec<double,3> &X, DistSVec<double,3> &dX) { return 0; }

  void printNodalDebug(int globNodeId, int identifier, DistSVec<double,3> *U, DistVec<int> *Id=0, DistVec<int> *Id0=0) {}

  void checkLocalRomStatus(DistSVec<double, 3> &, const int) {}
  void calculateSpatialResidual(DistSVec<double, 3> &, DistSVec<double,3> &) {}
  bool outputOnlySpatialResidual() {return false;}

  void writeBinaryVectorsToDiskRom(bool, int, int, DistSVec<double,3> *, DistSVec<double,3> *) {}
  void incrementNewtonOutputTag() {}
  int *getTimeIt() { return domain->getTimeIt(); }
  int *getNewtonIt() { return domain->getNewtonIt(); }
  int *getNumResidualsOutputCurrentNewtonIt() { return domain->getNumResidualsOutputCurrentNewtonIt(); }
  void setCurrentStateForKspBinaryOutput(DistSVec<double,3>&) {}
};

//------------------------------------------------------------------------------
//
class EmbeddedALETetMeshMotionSolver : public TetMeshMotionSolver {

public:

  EmbeddedALETetMeshMotionSolver(DefoMeshMotionData &, MatchNodeSet **, Domain *, MemoryPool *);
  ~EmbeddedALETetMeshMotionSolver(){};
  int solve(DistSVec<double,3> &, DistSVec<double,3> &);
};

//------------------------------------------------------------------------------

#endif
