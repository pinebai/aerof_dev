#ifndef _MULTIPHYSICS_TS_DESC_H_
#define _MULTIPHYSICS_TS_DESC_H_

#include <GhostPoint.h>
#include <TsDesc.h>

struct DistInfo;
class Domain;
class DynamicNodalTransfer;
class GeoSource;
class IoData;

template<int dim> class DistExactRiemannSolver;
template<class Scalar, int dim> class DistSVec;
template<int dimLS> class LevelSet;

//------------------------------------------------------------------------

template<int dim, int dimLS>
class MultiPhysicsTsDesc : public TsDesc<dim> , ForceGenerator<dim> {

 protected:

  double vfar[dim]; //farfield state
  DistExactRiemannSolver<dim> *riemann; //Riemann solver -- used at both FF and FS interfaces
  enum Type {EULER = 0, NAVIER_STOKES = 1} eqsType;
  const int numFluid;

  // time-step and related.
  double dtf;     //<! fluid time-step
  double dtfLeft; //<! time until next structure time-step is reached.
  double dts;     //<! structure time-step
  int globIt;         //<! current global(i.e. structure) iteration
  bool inSubCycling;  //<! is it in subcyling (i.e. itSc>1)
  bool inRedoTimestep;

  bool requireSpecialBDF;
  bool increasingPressure;

  int multiFluidInterfaceOrder;

  int limitHigherOrderExtrapolation;

  int phaseChangeType;

  //------------------------------------------------------------------------
  // EulerFSI: basic parameters
  bool withCracking;

  int orderOfAccuracy; // consistent with the reconstruction type for space
  bool linRecAtInterface, viscSecOrder;
  int simType;        // 0: steady-state    1: unsteady
  int riemannNormal;  // 0: struct normal;  1: fluid normal (w.r.t. control volume face)
  int numStructNodes;
  int totStructNodes;
  int forceApp; // now have four options.
                // = 1 : on GammaF, formula 1;
                // = 2 : on GammaF, formula 2; (not used)
                // = 3 : on Gamma*, formula 1;
                // = 4 : on Gamma*, formula 2; (not used)
  int phaseChangeChoice; // = 0. use nodal values.
                         // = 1. use solutions of Riemann problems.

  // EulerFSI: force calculation
  double (*Fs)[3]; //force distribution on the structure surfac3
  bool FsComputed; //whether Fs has been computed for this (fluid-)time step.

  bool existsWstarnm1;

  // EulerFSI: interface tracking
  DistLevelSetStructure *distLSS; //<! tool for FS tracking (not necessarily using level-sets) 

  DistSVec<double,dim> *Wstarij,*Wstarij_nm1;  //<! stores the FS Riemann solution (i->j) along edges
  DistSVec<double,dim> *Wstarji,*Wstarji_nm1;  //<! stores the FS Riemann solution (j->i) along edges
  DistSVec<double,dim> Vtemp;     //<! the primitive variables.
  DistSVec<double,dim> *VWeights; //<! stores U*Weights for each node. Used in updating phase change.
  DistVec<double> *Weights;       //<! weights for each node. Used in updating phase change.
  DistVec<GhostPoint<dim>*> *ghostPoints;

  DistVec<double> umax;

  int phaseChangeAlg;	 // = 0. use averaged value, given phaseChangeChocie==0
  						 // = 1. use least-squares, given phaseChangeChoice==0
  int interfaceAlg;		 // = 0. do not use information of intersection at surrogate interface
  						 // = 1. use information of intersection at surrogate interface
  double intersectAlpha; //	relevant only if interfaceAlg==1

  // EulerFSI: FS communication
  DistSVec<double,dim> Wtemp;
  DynamicNodalTransfer *dynNodalTransfer;

  //------------------------------------------------------------------------
  // MultiPhaseFlow: basics
  bool withMixedLS;
  MultiPhaseSpaceOperator<dim,dimLS> *multiPhaseSpaceOp;
  int frequencyLS; // frequency for reinitialization of level set

  // MultiPhaseFlow: level-sets and fluidIds
  LevelSet<dimLS> *LS;
  FluidSelector fluidSelector;
  DistSVec<double,dimLS> Phi;           //conservative variables
  DistSVec<double,dimLS> PhiWeights;    //<! stores Phi*Weights for each ndoe. Used in updating phase change
  DistSVec<double,dimLS> PhiV;          //primitive variables
  DistSVec<double,dim> V0;
  DistSVec<bool,2> InterfaceTag;

  //------------------------------------------------------------------------
  // buckling cylinder parameters
  // pressure is increased in the fluid at rate Prate from
  // initial pressure Pinit until it reaches the pressure
  // given by boundary conditions which happens at tmax.
  enum ImplosionSetupType {LINEAR = 0, SMOOTHSTEP = 1, NONE = 2} implosionSetupType;
  double tmax;
  double Prate;
  double Pinit;
  double Pfinal;
  double Pscale;
  int intersector_freq;
 
  int lsMethod;

  ProgrammedBurn* programmedBurn;

  double currentTime;

  int lastLSUpdateIteration;

 protected:
  void setupEmbeddedFSISolver(IoData &ioData);
  void setupMultiPhaseFlowSolver(IoData &ioData);
  /** computes the force load. Wij and Wji must be edge-based primitive state vectors. */ 
  void computeForceLoad(DistSVec<double,dim> *Wij, DistSVec<double,dim> *Wji);

 public:
  MultiPhysicsTsDesc(IoData &, GeoSource &, Domain *);
  ~MultiPhysicsTsDesc();

  //-- overrides the functions implemented in TsDesc.
  void setupTimeStepping(DistSVec<double,dim> *, IoData &);
  double computeTimeStep(int, double *, DistSVec<double,dim> &, double);
  double computeTimeStep(int a, double *b, DistSVec<double,dim> &c){ return computeTimeStep(a,b,c,-2.0);}
  void updateStateVectors(DistSVec<double,dim> &, int = 0);
  int checkSolution(DistSVec<double,dim> &);
  void setupOutputToDisk(IoData &, bool *, int, double,
                        DistSVec<double,dim> &);
  void outputToDisk(IoData &, bool*, int, int, int, double, double,
                        DistSVec<double,dim> &);
  void outputForces(IoData &, bool*, int, int, int, double, double,
                    DistSVec<double,dim> &);
  void outputPositionVectorToDisk(DistSVec<double,dim>&);
  void resetOutputToStructure(DistSVec<double,dim> &);
  /** Override the TsDesc routine because forces are sent to the structure
   * in a different way than the general case */
  void updateOutputToStructure(double, double, DistSVec<double,dim> &);

  double computeResidualNorm(DistSVec<double,dim>& );
  void monitorInitialState(int, DistSVec<double,dim>& );

  void getForcesAndMoments(map<int,int> & surfOutMap, DistSVec<double,dim> &U, DistSVec<double,3> &X,
                                           Vec3D* Fi, Vec3D* Mi);

  void getderivativeOfForcesAndMoments(map<int,int> & surfOutMap, 
				       DistSVec<double,dim> &U, DistSVec<double,dim> &dU, 
				       DistSVec<double,3> &X, double dS[3],
				       Vec3D *dFi, Vec3D *dMi);

  bool IncreasePressure(int it, double dt, double t, DistSVec<double,dim> &U);

  virtual int solveNonLinearSystem(DistSVec<double,dim> &, int)=0;
  virtual bool willNotSolve(double dts, double t) {return (t+dts*2)<tmax;}
  virtual void setFluidSubcycling(bool inSub) {inSubCycling = inSub;}


  void setCurrentTime(double t,DistSVec<double,dim>& U);
  double currentPressure(double t);

};

#ifdef TEMPLATE_FIX
#include <MultiPhysicsTsDesc.C>
#endif

#endif
