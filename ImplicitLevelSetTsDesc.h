#ifndef _IMPLICIT_LEVELSET_TS_DESC_H_
#define _IMPLICIT_LEVELSET_TS_DESC_H_

#include <LevelSetTsDesc.h>

class IoData;
class GeoSource;
class Domain;

struct DistInfo;
template<class Scalar, int dim> class DistSVec;

template<int dim, int dimLS> class MatVecProdMultiPhase;
template<int dim, int dimLS> class MatVecProdLS;
template<int dim, class Scalar2> class KspPrec;
#ifndef _KSPSLVR_TMPL_
#define _KSPSLVR_TMPL_
template<class VecType, class MvpOp, class PrecOp, class IoOp, class ScalarT = double> class KspSolver;
#endif


//------------------------------------------------------------------------

template<int dim, int dimLS>
class ImplicitLevelSetTsDesc : public LevelSetTsDesc<dim, dimLS> {


public:
  typedef DistSVec<double,dimLS> PhiVecType;
 protected: 

  DistSVec<double,dim> U0,Fold;
  
  DistEmbeddedVec<double,dim> embeddedU,embeddedB,embeddeddQ;
  int numBlur;
  
  DistSVec<bool,2> *tag;

  MatVecProdMultiPhase<dim,dimLS> *mvp;
  KspPrec<dim> *pc;
  KspSolver<DistEmbeddedVec<double,dim>, MatVecProdMultiPhase<dim,dimLS>, 
            KspPrec<dim>, Communicator> *ksp;
	
  MatVecProdLS<dim,dimLS> *mvpLS;
  KspPrec<dimLS> *pcLS;
  KspSolver<DistSVec<double,dimLS>, MatVecProdLS<dim,dimLS>, KspPrec<dimLS>, Communicator> *kspLS;

  NewtonSolver<ImplicitLevelSetTsDesc<dim,dimLS> > *ns;

  int failSafeNewton;
  int maxItsNewton;
  double epsNewton;
  double epsAbsResNewton, epsAbsIncNewton;
  FILE *outputNewton;
  int maxItsLS;
  double contractionLS, sufficDecreaseLS;

  DistVec<double>* hhResidual;
 public:
  ImplicitLevelSetTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitLevelSetTsDesc();
  
  int getMaxItsNewton() const { return maxItsNewton; }
  double getEpsNewton() const { return epsNewton; }
  double getEpsAbsResNewton() const { return epsAbsResNewton; }
  double getEpsAbsIncNewton() const { return epsAbsIncNewton; }
  FILE* getOutputNewton() const { return outputNewton; }
  int getLineSearch() const { return (maxItsLS>0); }
  int getMaxItsLineSearch() const { return maxItsLS; }
  double getContractionLineSearch() const { return contractionLS; }
  double getSufficientDecreaseLineSearch() const { return sufficDecreaseLS; }


  //-- functions for solving Euler equations
  int solveNonLinearSystem(DistSVec<double,dim> &, int);
  void computeFunction(int, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void recomputeFunction(DistSVec<double,dim> &, DistSVec<double,dim> &);
  
  int checkFailSafe(DistSVec<double,dim>&);
  void resetFixesTag();

  void computeJacobian(int, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void setOperators(DistSVec<double,dim> &);
  int solveLinearSystem(int, DistSVec<double,dim> &, DistSVec<double,dim> &);

  //-- new functions for solving LevelSet equation
  void computeFunctionLS(int, DistSVec<double,dim> &,DistSVec<double,dimLS> &);
  void computeJacobianLS(int, DistSVec<double,dim> &,DistSVec<double,dimLS> &);
  void setOperatorsLS(DistSVec<double,dimLS> &Q);

  int solveLinearSystemLS(int, DistSVec<double,dimLS> &, DistSVec<double,dimLS> &);

  void doErrorEstimation(DistSVec<double,dim> &U);

  void setCurrentStateForKspBinaryOutput(DistSVec<double,dim> &Q) {}

 protected:
  template<class Scalar, int neq>
  KspPrec<neq> *createPreconditioner(PcData &, Domain *);

  template<int neq, class MatVecProdOp>
  KspSolver<DistEmbeddedVec<double,neq>, MatVecProdOp, KspPrec<neq>,
  Communicator> *createKrylovSolver(const DistInfo &, KspData &, 
                                    MatVecProdOp *, KspPrec<neq> *, Communicator *);

  template<int neq, class MatVecProdOp>
  KspSolver<DistSVec<double,neq>, MatVecProdOp, KspPrec<neq>,
  Communicator> *createKrylovSolverLS(const DistInfo &, KspData &, 
                                      MatVecProdOp *, KspPrec<neq> *, Communicator *);
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitLevelSetTsDesc.C>
#endif

#endif
