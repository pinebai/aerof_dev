#ifndef _DIST_TIME_STATE_H_
#define _DIST_TIME_STATE_H_

#include <TimeData.h>
#include <DistVector.h>
#include <DistMacroCell.h>
#include <LowMachPrec.h>
#include <LevelSet/LevelSetStructure.h>
#include <ErrorHandler.h>
#include <RefVal.h>

struct FluidModelData;
struct InitialConditions;

class VarFcn;
class Domain;
class DistGeoState;
class FemEquationTerm;
class DistLevelSetStructure;

template<int dim> class TimeState;
template<class Scalar, int dim> class DistMat;
template<int dim> class SpaceOperator;
template<int dim> class DistExactRiemannSolver;

//------------------------------------------------------------------------------

template<int dim>
class DistTimeState {

private:

  VarFcn *varFcn;
  RefVal *refVal;
  FemEquationTerm *fet;
  DistSVec<double,dim> *V;
  DistSVec<double,dim> *Vn;
  DistSVec<double,dim> *VnBar;
  DistSVec<double,dim> *QBar;

  DistVec<int> *firstOrderNodes;
  
  bool locAlloc;
  int numLocSub;

  double gam;
  double pstiff;

  double refTime;

  TimeLowMachPrec    tprec;
  SpatialLowMachPrec sprec; //only for computation of irey
  /*bool prec;
  double beta;
  double mach;
  double cmach;
  double k1;
  double betav;
*/
  double viscousCst;

  TimeData *data;

  DistVec<double> *dt;			//actual   time stepping
  DistVec<double> *idti;		//inverse inviscid time stepping
  DistVec<double> *idtv;		//inverse viscous  time stepping
  DistVec<double> *dtau;		//dual time stepping
  DistVec<double> *irey;
  DistSVec<double,dim> *Un;
  DistSVec<double,dim> *Unm1;
  DistSVec<double,dim> *Unm2;
  DistSVec<double,dim> *UnBar;
  DistSVec<double,dim> *Unm1Bar;
  DistSVec<double,dim> *Unm2Bar;
  DistSVec<double,dim> *Rn;
  double errorEstiNorm;                 //norm of estimated error
  double dtMin;

  Domain *domain;

  TimeState<dim> **subTimeState;

// Included (MB)
  DistVec<double> *dIrey;
  DistVec<double> *dIdti;
  DistVec<double> *dIdtv;

  bool isGFMPAR;
  int fvmers_3pbdf ;

  int mf_phase_change_type;

  double dt_coeff;
  int dt_coeff_count;

  double checkForRapidlyChangingPressure;
  double checkForRapidlyChangingDensity;

  DistVec<double>* hhn,*hhnm1;
  ErrorHandler* errorHandler;

  bool checkForRapidlyChangingValues;

public:
  bool allowcflstop;
  bool allowdtstop;

private:
  void computeInitialState(InitialConditions &ic, FluidModelData &fm, double UU[dim]);
public:

  TimeLowMachPrec& getTimeLowMachPrec() { return tprec; }
  SpatialLowMachPrec& getSpatialLowMachPrec() { return sprec; }

  DistTimeState(IoData &, SpaceOperator<dim> *, VarFcn *, Domain *, DistSVec<double,dim> * = 0);
  DistTimeState(IoData &, SpaceOperator<dim> *, VarFcn *, Domain *,
                DistInfo& dI, DistSVec<double,dim> * = 0);
  DistTimeState(const DistTimeState<dim> &, bool, IoData &);
  ~DistTimeState();

  void attachHH(DistVec<double>&);

  void copyTimeData(DistTimeState<dim>* oth);

  void disableRapidlyChangingValueCheck() { checkForRapidlyChangingValues = false; }

  DistVec<double>& getDt() const { return *dt; }

  void createSubStates();

  void initialize(IoData &ioData, SpaceOperator<dim> *spo, VarFcn *vf,
		  Domain *dom, DistSVec<double,dim> *v, DistInfo& dI) ;

  TimeState<dim> &operator() (int i) const { return *subTimeState[i]; }

  void setResidual(DistSVec<double,dim> *rn) { if (Rn != 0) delete Rn; Rn = rn; }

  void setGlobalTimeStep (double t) { *dt = t; }

  void setup(const char *name, DistSVec<double,3> &X, DistSVec<double,dim> &Ufar,
             DistSVec<double,dim> &U, IoData &iod, DistVec<int> *fluidId = 0); 
  void setupUVolumesInitialConditions(IoData &iod);
  void setupUOneDimensionalSolution(IoData &iod, DistSVec<double,3> &X);
  void setupUMultiFluidInitialConditions(IoData &iod, DistSVec<double,3> &X);
  void setupUFluidIdInitialConditions(IoData &iod, DistVec<int> &fluidId);
  void setupUExactSolutionInitialConditions(IoData &iod, DistSVec<double,3> &X);

  void update(DistSVec<double,dim> &,bool increasingPressure = false);
  void update(DistSVec<double,dim> &Q,  DistSVec<double,dim> &Qtilde,DistVec<int> &fluidId, DistVec<int> *fluidIdnm1, 
              DistExactRiemannSolver<dim> *riemann,class DistLevelSetStructure* = 0, bool increasingPressure = false);

  void updateHH(DistVec<double> & hh);

  void writeToDisk(char *);

  double computeTimeStep(double, double, double*, int*, DistGeoState &, 
								 DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistLevelSetStructure *distLSS=0);

  double computeTimeStep(double, double, double*, int*, DistGeoState &,
                         DistVec<double> &, DistSVec<double,dim> &, DistVec<int> &, 
								 DistVec<double>* = NULL);

  //computes time step size when time step failed
  double computeTimeStepFailSafe(double* , int*);

  //computes time step with error estimation
  //see: "Adaptive Time Stepping for Non-Linear Hyperbolic Problems"
  double computeTimeStep(int, double* , int*);

  //calculates the norm of the error estimator
  void calculateErrorEstiNorm(DistSVec<double,dim> &, DistSVec<double,dim> &);

  //minimum time step size is set for error estimation
  void setDtMin(double dt){dtMin = dt;}

  void computeCoefficients(double);

  void add_dAW_dt(int, DistGeoState &, DistVec<double> &, 
		  DistSVec<double,dim> &, DistSVec<double,dim> &, DistLevelSetStructure *distLSS=0);

  void add_dAW_dt_HH(int, DistGeoState &, DistVec<double> &, 
		     DistVec<double> &, DistVec<double> &);

  void add_dAW_dtRestrict(int, DistGeoState &, DistVec<double> &, 
			  DistSVec<double,dim> &, DistSVec<double,dim> &, const std::vector<std::vector<int> > &);
  template<int dimLS>
  void add_dAW_dtLS(int, DistGeoState &, DistVec<double> &, 
			 DistSVec<double,dimLS> &, DistSVec<double,dimLS> &, DistSVec<double,dimLS> &, 
			 DistSVec<double,dimLS> &, DistSVec<double,dimLS> &,bool requireSpecialBDF = false);
  void add_dAW_dtau(int, DistGeoState &, DistVec<double> &, 
		  DistSVec<double,dim> &, DistSVec<double,dim> &, DistLevelSetStructure *distLSS=0);

  template<class Scalar, int neq>
  void addToJacobian(DistVec<double> &, DistMat<Scalar,neq> &, DistSVec<double,dim> &);

  template<class Scalar, int neq>
  void addToHHJacobian(DistVec<double> &, DistMat<Scalar,neq> &, DistVec<double> &);

  template<class Scalar, int neq>
  void addToJacobianNoPrec(DistVec<double> &, DistMat<Scalar,neq> &, DistSVec<double,dim> &);

  template<class Scalar, int neq>
  void addToJacobianLS(DistVec<double> &, DistMat<Scalar,neq> &, DistSVec<double,dim> &,bool);

  template<class Scalar, int neq>
  void addToJacobianGasPrec(DistVec<double> &, DistMat<Scalar,neq> &, DistSVec<double,dim> &);

  template<class Scalar, int neq>
  void addToJacobianLiquidPrec(DistVec<double> &, DistMat<Scalar,neq> &, DistSVec<double,dim> &);

  template<class Scalar, int neq>
  void addToH1(DistVec<double> &, DistMat<Scalar,neq> &);

  template<class Scalar, int neq>
  void addToH1(DistVec<double> &, DistMat<Scalar,neq> &, Scalar);

  template<class Scalar, int neq>
  void addToH2(DistVec<double> &, DistSVec<double,dim> &, DistMat<Scalar,neq> &);

  template<class Scalar, int neq>
  void addToH2(DistVec<double> &, DistSVec<double,dim> &, DistMat<Scalar,neq> &, Scalar);

  template<class Scalar, int neq>
  void addToH2(DistVec<double> &, DistSVec<double,dim> &, DistMat<Scalar,neq> &, Scalar, double);
  
  template<class Scalar, int neq>
  void addToH2Minus(DistVec<double> &, DistSVec<double,dim> &, DistMat<Scalar,neq> &);

  void computeBar(bool, DistMacroCellSet *, DistGeoState &, int);
                                                                                                                          
  void get_dW_dt(bool, DistGeoState &, DistVec<double> &, DistSVec<double,dim> &,
                 DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
                 DistMacroCellSet *, DistSVec<double,1> **, int);

  DistVec<double>* getInvReynolds(){ return irey; }
  
  void multiplyByTimeStep(DistVec<double>&);                                                                        
  void multiplyByTimeStep(DistSVec<double,dim>&);
  template<int dimLS>
  void multiplyByTimeStep(DistSVec<double,dimLS>&);
  void multiplyByPreconditioner(DistSVec<double,dim> &, DistSVec<double,dim>&);
  void multiplyByPreconditionerPerfectGas(DistSVec<double,dim> &, DistSVec<double,dim>&);
  void multiplyByPreconditionerLiquid(DistSVec<double,dim> &, DistSVec<double,dim>&);

  TimeData &getData() { return *data; }
  DistSVec<double,dim> &getUn() const { return *Un; }
  DistSVec<double,dim> &getUnm1() const { return *Unm1; }

  inline bool existsNm1() const { return data->exist_nm1; }
  inline bool useNm1() const { return data->use_nm1; }

  void setExistsNm1();
  void setDtNm1(double dt);
    void setDt(double dt){ *(this->dt) = dt; }

  double getTime()  { return data->dt_n; }
  
  void updateDtCoeff();

  void rstVar(IoData &);

  DistVec<double>* getDerivativeOfInvReynolds(DistGeoState &, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &, double);

	double getNewtonTag() const;
  int getNewtonStateStep() const;
  int getNewtonResidualStep() const;
  int getKrylovStep() const;

  DistVec<int> * getFirstOrderNodeSet() const { return firstOrderNodes; }
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <DistTimeState.C>
#endif

#endif
