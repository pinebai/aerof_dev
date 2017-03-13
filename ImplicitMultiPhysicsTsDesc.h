#ifndef _IMPLICIT_MULTIPHYSICS_TS_DESC_H_
#define _IMPLICIT_MULTIPHYSICS_TS_DESC_H_

#include <MultiPhysicsTsDesc.h>
#include <DistVector.h>
#include <MatVecProd.h>

class IoData;
class Domain;
class GeoSource;
//template<class Scalar, int dim> class DistVec;

//------------------------------------------------------------------------

template<int dim, int dimLS>
class ImplicitMultiPhysicsTsDesc : public MultiPhysicsTsDesc<dim,dimLS> {

 public:
  typedef DistSVec<double,dimLS> PhiVecType;
 protected:
  DistSVec<bool,2> *tag;

  MatVecProdMultiPhase<dim,dimLS> *mvp;
  KspPrec<dim> *pc;
  KspSolver<DistEmbeddedVec<double,dim>, MatVecProdMultiPhase<dim,dimLS>, 
            KspPrec<dim>, Communicator> *ksp;

  MatVecProdLS<dim,dimLS> *mvpLS;
  KspPrec<dimLS> *pcLS;
  KspSolver<DistSVec<double,dimLS>, MatVecProdLS<dim,dimLS>, KspPrec<dimLS>, Communicator> *kspLS;

  NewtonSolver<ImplicitMultiPhysicsTsDesc<dim,dimLS> > *ns;

  DistEmbeddedVec<double,dim> embeddedU,embeddedB,embeddeddQ;

  DistVec<double>* hhResidual;

  int failSafeNewton;
  int maxItsNewton;
  double epsNewton;
  double epsAbsResNewton, epsAbsIncNewton;
  FILE *outputNewton;
  int maxItsLS;
  double contractionLS, sufficDecreaseLS;

 public:
  ImplicitMultiPhysicsTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitMultiPhysicsTsDesc();

  int solveNonLinearSystem(DistSVec<double,dim> &U, int);

  int getMaxItsNewton() const { return maxItsNewton; }
  double getEpsNewton() const { return epsNewton; }
  double getEpsAbsResNewton() const { return epsAbsResNewton; }
  double getEpsAbsIncNewton() const { return epsAbsIncNewton; }
  FILE* getOutputNewton() const { return outputNewton; }
  int getLineSearch() const { return (maxItsLS>0); }
  int getMaxItsLineSearch() const { return maxItsLS; }
  double getContractionLineSearch() const { return contractionLS; }
  double getSufficientDecreaseLineSearch() const { return sufficDecreaseLS; }

  void computeFunction(int, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void recomputeFunction(DistSVec<double,dim> &, DistSVec<double,dim> &);

  int checkFailSafe(DistSVec<double,dim>&);
  void resetFixesTag();

  void computeJacobian(int, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void setOperators(DistSVec<double,dim> &);
  int solveLinearSystem(int, DistSVec<double,dim> &, DistSVec<double,dim> &);

  void commonPart(DistSVec<double,dim> &U);

  //-- new functions for solving LevelSet equation
  void computeFunctionLS(int, DistSVec<double,dim> &,DistSVec<double,dimLS> &);
  void computeJacobianLS(int, DistSVec<double,dim> &,DistSVec<double,dimLS> &);
  void setOperatorsLS(DistSVec<double,dimLS> &Q);

  int solveLinearSystemLS(int, DistSVec<double,dimLS> &, DistSVec<double,dimLS> &);

  template<int neq>
  KspPrec<neq> *createPreconditioner(PcData &, Domain *);

  template<int neq, class MatVecProdOp>
  KspSolver<DistEmbeddedVec<double,neq>, MatVecProdOp, KspPrec<neq>, Communicator> *
  createKrylovSolver(const DistInfo &info, KspData &kspdata,
                     MatVecProdOp *_mvp, KspPrec<neq> *_pc,
                     Communicator *_com);

  template<int neq, class MatVecProdOp>
  KspSolver<DistSVec<double,neq>, MatVecProdOp, KspPrec<neq>,
  Communicator> *createKrylovSolverLS(const DistInfo &, KspData &, 
                                      MatVecProdOp *, KspPrec<neq> *, Communicator *);

  void setCurrentStateForKspBinaryOutput(DistSVec<double,dim> &Q) {}
 private:

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitMultiPhysicsTsDesc.C>
#endif

#endif

