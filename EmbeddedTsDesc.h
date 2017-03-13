#ifndef _EMBEDDED_TS_DESC_H_
#define _EMBEDDED_TS_DESC_H_

#include <TsDesc.h>

#include <IoData.h>
#include <PostOperator.h>
#include <ReinitializeDistanceToWall.h>
#include <GhostPoint.h>

struct DistInfo;
class DynamicNodalTransfer;
class GeoSource;
class Domain;
template<class Scalar, int dim> class DistSVec;
template<int dim> class DistExactRiemannSolver;

//------------------------------------------------------------------------

template<int dim>
class EmbeddedTsDesc : public TsDesc<dim> , ForceGenerator<dim> {

 protected:

  IoData& ioData;

  DistExactRiemannSolver<dim> *riemann; //Riemann solver -- used at both FF and FS interfaces
  double vfar[dim]; //farfield state

  DistVec<int> nodeTag; // = 1 for fluid #1; = -1 for fluid #2.
  DistVec<int> nodeTag0; // node tag for the previous time-step.

  double timeStep;

  bool withCracking;

  double (*Fs)[3]; //force distribution on the structure surfac3
  bool FsComputed; //whether Fs has been computed for this (fluid-)time step.
  int numStructNodes;
  int totStructNodes;
  bool linRecAtInterface, viscSecOrder;
  int simType;        // 0: steady-state    1: unsteady
  int riemannNormal;  // 0: struct normal;  1: fluid normal (w.r.t. control volume face)

  double (*dFs)[3]; // derivative of force distribution on the structure surface

  bool increasingPressure;
  bool recomputeIntersections;
  double unifPressure[2];
 
  // ----------- time steps -----------------------------------------------------------
  double dtf;     //<! fluid time-step
  double dtfLeft; //<! time until next structure time-step is reached.
  double dts;     //<! structure time-step
  int globIt;         //<! current global(i.e. structure) iteration
  bool inSubCycling;  //<! is it in subcyling (i.e. itSc>1)
  // ----------------------------------------------------------------------------------

  bool existsWstarnm1;

  // ----------- components for Fluid-Structure interface. -----------------------------
  DistLevelSetStructure *distLSS; //<! tool for FS tracking (not necessarily a  "levelset solver".)
  DistVec<int> *countWstarij, *countWstarji;   //<! only used if ioData.embed.interfaceAlg==INTERSECTION
  DistSVec<double,dim> *Wstarij,*Wstarij_nm1;  //<! stores the FS Riemann solution (i->j) along edges
  DistSVec<double,dim> *Wstarji,*Wstarji_nm1;  //<! stores the FS Riemann solution (j->i) along edges
  DistSVec<double,dim> Vtemp;     //<! the primitive variables.
  DistSVec<double,dim> *VWeights; //<! stores U*Weights for each node. Used in updating phase change.
  DistVec<double> *Weights;       //<! weights for each node. Used in updating phase change.

  DistSVec<double,dim> *Wextij;

  ReinitializeDistanceToWall<1> *wall_computer;
  // ------------------------------------------------------------------------------------

  // Copies for fail safe ----- -----------------------------
  DistSVec<double,dim> *WstarijCopy,*Wstarij_nm1Copy;  //<! stores the FS Riemann solution (i->j) along edges
  DistSVec<double,dim> *WstarjiCopy,*Wstarji_nm1Copy;  //<! stores the FS Riemann solution (j->i) along edges
  DistVec<int> *nodeTagCopy; // = 1 for fluid #1; = -1 for fluid #2.
  DistVec<int> *nodeTag0Copy; // node tag for the previous time-step.
  DistSVec<double,dim> *UCopy;     //<! the primitive variables.
  // ------------------------------------------------------------------------------------

  DistSVec<double,dim> Wtemp;
  DynamicNodalTransfer *dynNodalTransfer;
  MeshMotionHandler* emmh;

  //buckling cylinder parameters
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

  double currentTime;
  double currentTimeStep;

  // Adam 04/06/2010
  enum Type {EULER = 0, NAVIER_STOKES = 1} eqsType;

  DistVec<GhostPoint<dim>*> *ghostPoints;
  // ghostPoints is a pointer on a vector of pointer on GhostPoint<dim> objects. Each subdomain can be accessed separatly.

 public:
  int orderOfAccuracy; // consistent with the reconstruction type for space
  int forceApp; // now have four options.
                // = 1 : on GammaF, formula 1;
                // = 2 : on GammaF, formula 2;
                // = 3 : on Gamma*, formula 1;
                // = 4 : on Gamma*, formula 2;
  int phaseChangeChoice; // = 0. use nodal values.
                         // = 1. use solutions of Riemann problems.
  int phaseChangeAlg;	 // = 0. use averaged value, given phaseChangeChocie==0
  						 // = 1. use least-squares, given phaseChangeChoice==0
  int interfaceAlg;		 // = 0. do not use information of intersection at surrogate interface
  						 // = 1. use information of intersection at surrogate interface
  double intersectAlpha; //	relevant only if interfaceAlg==1
  const int numFluid;   //numFluid = 1 (for fluid-fullbody)
                            //     = 2 (for fluid-shell-fluid)

  EmbeddedTsDesc(IoData &, GeoSource &, Domain *);
  ~EmbeddedTsDesc();


  //-- overrides the functions implemented in TsDesc.
  void setupTimeStepping(DistSVec<double,dim> *, IoData &);
  double computeTimeStep(int, double *, DistSVec<double,dim> &, double);
  double computeTimeStep(int a, double * b, DistSVec<double,dim> & c){ return computeTimeStep(a,b,c,-2.0);}
  double computePositionVector(bool *, int, double, DistSVec<double,dim> &);
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

  void computeForceLoad(DistSVec<double,dim> *Wij, DistSVec<double,dim> *Wji);
  /** computes the force load. Wij and Wji must be edge-based primitive state vectors. */ 

  void computederivativeOfForceLoad(DistSVec<double,dim> *Wij, 
				    DistSVec<double,dim> *Wji, 
				    double dS[3], 
				    DistSVec<double,dim> &dV);

  virtual int solveNonLinearSystem(DistSVec<double,dim> &, int)=0;

  void getForcesAndMoments(map<int,int> & surfOutMap, DistSVec<double,dim> &U, DistSVec<double,3> &X,
                                           Vec3D* Fi, Vec3D* Mi);

  void getderivativeOfForcesAndMoments(map<int,int> & surfOutMap, 
				       DistSVec<double,dim> &V, DistSVec<double,dim> &dV, 
				       DistSVec<double,3> &X, double dS[3],
				       Vec3D *dFi, Vec3D *dMi);

  bool IncreasePressure(int it, double dt, double t, DistSVec<double,dim> &U);
  virtual bool willNotSolve(double dts, double t) {return (t+dts*2)<tmax;}
  virtual void setFluidSubcycling(bool inSub) {inSubCycling = inSub;}

  void fixSolution(DistSVec<double,dim>& U,DistSVec<double,dim>& dU);
  double currentPressure(double t);

  void computeDistanceToWall(IoData &ioData);

  MeshMotionHandler *createEmbeddedALEMeshMotionHandler(IoData &, GeoSource &, DistLevelSetStructure *);

  void computeConvergenceInformation(IoData &ioData, const char* file, DistSVec<double,dim>& U);

  void setCurrentTime(double t,DistSVec<double,dim>& U);
  void setCurrentTimeStep(double dt);
    void Bool2Char(DistVec<bool> &X, DistSVec<char, dim> &Y); //<! Lei Lei, 24 March 2016
    void writeBinaryVectorsToDiskRom(bool lastNewtonIt, int timeStep, int newtonIter,
                                     DistSVec<double, dim> *state, DistSVec<double, dim> *residual); //<! Lei Lei, 04 July 2016
};


#ifdef TEMPLATE_FIX
#include <EmbeddedTsDesc.C>
#endif

#endif

