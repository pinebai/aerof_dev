#ifndef _LEVELSET_TS_DESC_H_
#define _LEVELSET_TS_DESC_H_

#include <TsDesc.h>

#include <ProgrammedBurn.h>

#include <HigherOrderMultiFluid.h>

#include <OneDimensionalSolver.h>

#include <TriangulatedInterface.h>

class IoData;
class GeoSource;
class Domain;
class FluidSelector;
template<int dimLS> class LevelSet;
template<int dim> class DistExactRiemannSolver;

struct DistInfo;
template<class Scalar, int dim> class DistSVec;

//------------------------------------------------------------------------

template<int dim, int dimLS>
class LevelSetTsDesc : public TsDesc<dim> {

 protected:

  IoData& ioData;

  MultiPhaseSpaceOperator<dim,dimLS> *multiPhaseSpaceOp;
  FluidSelector fluidSelector;
  LevelSet<dimLS> *LS;

  TriangulatedInterface* myTriangulatedInterface;

  DistExactRiemannSolver<dim> *riemann;
  DistSVec<double,dimLS> Phi;           //conservative variables
  DistSVec<double,dimLS> PhiV;          //primitive variables
  DistSVec<double,dim> V0;
  DistSVec<double,dim> Utilde;

  DistVec<double> umax;

  // frequency for reinitialization of level set
  int frequencyLS;

  MultiFluidData::InterfaceType interfaceType; //to advance levelset or not
  
  //buckling cylinder parameters
  // pressure is increased in the fluid at rate Prate from
  // initial pressure Pinit until it reaches the pressure
  // given by boundary conditions which happens at tmax.
  double tmax;
  double Prate;
  double Pinit;

  bool requireSpecialBDF;

  int lsMethod;

  int interfaceOrder;

  int phaseChangeType;
  
  DistVec<double> Weights;
  DistSVec<double,dim> VWeights;

  struct exactInterfacePoint {

    double time,loc;
  };

  std::vector< exactInterfacePoint > myExactInterface;

  bool limitHigherOrderExtrapolation;

 public:
  LevelSetTsDesc(IoData &, GeoSource &, Domain *);
  ~LevelSetTsDesc();

  //-- overrides the functions implemented in TsDesc.
  void setupTimeStepping(DistSVec<double,dim> *, IoData &);
  double computeTimeStep(int, double *, DistSVec<double,dim> &, double);
  double computeTimeStep(int a, double * b, DistSVec<double,dim> & c){ return computeTimeStep(a,b,c,-2.0);}
  void updateStateVectors(DistSVec<double,dim> &, int = 0);
  int checkSolution(DistSVec<double,dim> &);
  void setupOutputToDisk(IoData &, bool *, int, double, 
                        DistSVec<double,dim> &);
  void outputToDisk(IoData &, bool*, int, int, int, double, double, 
		    	DistSVec<double,dim> &);
  void outputForces(IoData &, bool*, int, int, int, double, double,
                    DistSVec<double,dim> &);
  void resetOutputToStructure(DistSVec<double,dim> &);
  void updateOutputToStructure(double, double, DistSVec<double,dim> &);

  bool IncreasePressure(int it, double dt, double t, DistSVec<double,dim> &U);
  
  void fixSolution(DistSVec<double,dim>& U, DistSVec<double,dim>& dU);

  virtual int solveNonLinearSystem(DistSVec<double,dim> &, int)=0;
  
  void setCurrentTime(double t,DistSVec<double,dim>& U);

  void computeConvergenceInformation(IoData &ioData, const char* file, DistSVec<double,dim>& U);

  void loadExactInterfaceFile(IoData& ioData, const char* file);

  void setPhiExact();
 protected:
  void avoidNewPhaseCreation(DistSVec<double,dimLS> &localPhi);

  double currentTime,progBurnIgnitionTime;
  bool progBurnIgnited;

  ProgrammedBurn* programmedBurn;

  double currentTimeStep;
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <LevelSetTsDesc.C>
#endif

#endif
