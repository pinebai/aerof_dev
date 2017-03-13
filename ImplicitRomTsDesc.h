#ifndef _IMPLICIT_ROM_TS_DESC_H_
#define _IMPLICIT_ROM_TS_DESC_H_

#include <IoData.h>
#include <TsDesc.h>
#include <KspPrec.h>
#include <NonlinearRom.h>
#include <NonlinearRomOnlineII.h>
#include <NonlinearRomOnlineIII.h>

struct DistInfo;

class GeoSource;
class Domain;
class Communicator;

template<class Scalar, int dim> class DistSVec;
template<int dim, int neq> class MatVecProd;


//------------------------------------------------------------------------------

template<int dim>
class ImplicitRomTsDesc : public TsDesc<dim> {

protected:

  IoData *ioData;

  int maxItsNewton;
  double epsNewton;
  double epsAbsResNewton;
  double epsAbsIncNewton;

  MatVecProd<dim,dim> *mvp;

  DistSVec<bool,2> *tag;

  VecSet<DistSVec<double, dim> > pod;

  NonlinearRom<dim>* rom;

  int currentCluster;

  int basisUpdateFreq;  
  int tryAllFreq;

  FullM jac;
  
  int nPod;
  
  bool unsteady;
  bool systemApprox;
  bool useIncrements;
  bool tryingAllClusters;

  DistSVec<double, dim> F;	// residual
  VecSet<DistSVec<double, dim> > AJ; // Action of Jacobian (AJ) on reduced-order basis

  DistSVec<double, dim>* componentwiseScalingVec; // Intended to fix the scaling of the residual.  Gives energy interpretation to all equations.

  double levenbergMarquardtWeight;

  // backtracking line search
  double rho;
  double c1;
  int maxItsLS; 

  Vec<double> dUromNewtonIt;    // set to zero before each newton iteration
  Vec<double> dUromTimeIt;      // set to zero before each time iteration
  Vec<double> dUromCurrentROB;  // set to zero after each cluster switch
  Vec<double> UromCurrentROB;   // for projection only: initialized at each cluster switch

  double target, res0;	// for Newton convergence

  virtual void computeAJ(int, DistSVec<double, dim> &, bool applyWeighting = false, DistSVec<double, dim> *R = NULL);
  virtual void computeRedHessianSums(int, DistSVec<double, dim> &);	// Broyden doesn't do this every time 
  virtual void computeFullResidual(int, DistSVec<double, dim> &, bool applyWeighting = false, DistSVec<double, dim> *R = NULL, bool includeHomotopy = true); 
  virtual void saveNewtonSystemVectors(const int _it) {};	// only implemented for PG/Galerkin
  void saveNewtonSystemVectorsAction(const int);	// implementation for PG/Galerkin
  virtual void solveNewtonSystem(const int &it, double &res, bool &breakloop, DistSVec<double, dim> &, const int &totalTimeSteps = 0) = 0;
	// each ROM has a different way of solving the Newton system
  int solveLinearSystem(int, Vec<double> &, Vec<double> &);
  virtual double meritFunction(int, DistSVec<double, dim> &, DistSVec<double, dim> &, DistSVec<double, dim> &, double) = 0; // depends on MOR approach
  double meritFunctionDeriv(int, DistSVec<double, dim> &, DistSVec<double, dim> &, DistSVec<double, dim> &, double);
  double lineSearch(DistSVec<double, dim> &, Vec<double> &, int, VecSet<DistSVec<double, dim> > &,double, bool &);
  double lineSearchBacktrack(DistSVec<double, dim> &, Vec<double> &, int, VecSet<DistSVec<double, dim> > &,double, bool &);
  double zoom(double, double, double, double, double, double, double, double, double, DistSVec<double,dim>,DistSVec<double,dim>, DistSVec<double,dim>, int);
  int checkFailSafe(DistSVec<double,dim>&);
  void resetFixesTag();
  void projectVector(VecSet<DistSVec<double, dim> >&, DistSVec<double, dim> &, Vec<double> &);
  void expandVector(Vec<double> &, DistSVec<double, dim> &);

  void loadCluster(int, bool, DistSVec<double, dim> &);

  virtual void postProStep(DistSVec<double,dim> &, int) {};	// by default, do not do post processing
  virtual bool breakloop1(const bool);
  virtual bool breakloop2(const bool);

  virtual void setReferenceResidual() {};
  virtual void setProblemSize(DistSVec<double, dim> &) {};
  virtual void deleteRestrictedQuantities() {};

  double *projVectorTmp; // temporary vector for projectVector

  bool updateFreq;
  bool clusterSwitch;
  bool updatePerformed;

  bool checkSolutionInNewton;

  DistSVec<double, dim>* Uinit;  // initial condition of the steady state simulation, 
                                 // stored to recalculate reference residual after 
                                 // cluster switch (new sampled mesh for GNAT) or after
                                 // changing the residual weighting 

  DistSVec<double, dim>* Uprev;  // solution at the beginning of the previous time step (needed for model II incremental bases) 

  double homotopyStepInitial;  // Reduced coordinate pseudo-time-stepping for spatial-only simulations. Only implemented for normal equations.
  double homotopyStepMax;  
  double homotopyStepGrowthRate;
  bool spatialOnlyWithHomotopy;

  void tryAllClusters(DistSVec<double, dim>&, const int totalTimeSteps, int*);

protected:
  template<class Scalar, int neq>
  KspPrec<neq> *createPreconditioner(PcData &, Domain *);

public:
  
  ImplicitRomTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitRomTsDesc();

  int solveNonLinearSystem(DistSVec<double, dim> &, const int _it);
  void rstVarImplicitRomTsDesc(IoData &);
  void checkLocalRomStatus(DistSVec<double, dim> &, const int);
  //void setInterpWeightsForMultiIC(std::vector<double> vec) {rom->interpWeightsForMultiIC = vec;}
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitRomTsDesc.C>
#endif

#endif
