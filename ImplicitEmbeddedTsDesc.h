#pragma once

#include <EmbeddedTsDesc.h>

#include "IoData.h"
#include "GeoSource.h"
#include "Domain.h"
#include "LevelSet.h"
#include "DistTimeState.h"

#include <MatVecProd.h>
#include <KspSolver.h>
#include <SpaceOperator.h>
#include <NewtonSolver.h>
#include <DistEmbeddedVector.h>

#ifdef TYPE_PREC
#define PrecScalar TYPE_PREC
#else
#define PrecScalar double
#endif

//------------------------------------------------------------------------------

template<int dim>
class ImplicitEmbeddedTsDesc : public EmbeddedTsDesc<dim> {

protected:

  DistSVec<bool,2> *tag;

  NewtonSolver<ImplicitEmbeddedTsDesc<dim> > *ns;

  int failSafeNewton;
  int maxItsNewton;
  double epsNewton;
  double epsAbsResNewton, epsAbsIncNewton;
  FILE *outputNewton;
  
  int maxItsLS;
  double contractionLS, sufficDecreaseLS;

  template <int neq>
  KspPrec<neq> *createPreconditioner(PcData &pcdata, Domain *dom);
 
  DistEmbeddedVec<double,dim> embeddedU,embeddedB,embeddeddQ;

  DistVec<double>* hhResidual;

public:
  
  virtual DistMat<double,dim>* GetJacobian() { 
    return NULL;
  }

  ImplicitEmbeddedTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom);

  ~ImplicitEmbeddedTsDesc();

  template<int neq, class MatVecProdOp>
  KspSolver<DistEmbeddedVec<double,neq>, MatVecProdOp, KspPrec<neq>, Communicator> *
  createKrylovSolver(const DistInfo &info, KspData &kspdata,
                     MatVecProdOp *_mvp, KspPrec<neq> *_pc,
                     Communicator *_com);

  int commonPart(DistSVec<double,dim> &U);

  int solveNonLinearSystem(DistSVec<double,dim> &U, int);
  
  void computeFunction(int it, DistSVec<double,dim> &Q,
                       DistSVec<double,dim> &F);

  void recomputeFunction(DistSVec<double,dim> &Q,
                         DistSVec<double,dim> &rhs);

  void doErrorEstimation(DistSVec<double,dim> &U);

  int checkFailSafe(DistSVec<double,dim>& U);

  void resetFixesTag();

  virtual void computeJacobian(int it, DistSVec<double,dim> &Q,
		       DistSVec<double,dim> &F) = 0;

  virtual void setOperators(DistSVec<double,dim> &Q) = 0;

  virtual int solveLinearSystem(int it, DistSVec<double,dim> &b,
			DistSVec<double,dim> &dQ) = 0;
    /**
     * (Lei Lei) only overridden in subclass ImplicitEmbeddedRomTsDesc
     * similiar in function to projectSwitchStateOntoAffineSubspace
     */
    virtual void projectStateOntoROB(DistSVec<double, dim> &U) {}
    /**
     * (Lei Lei) only overridden in subclass ImplicitEmbeddedRomTsDesc
     * testing if ghost point phase change is done
     */
    virtual void printGhostPoint() {}

  int getMaxItsNewton() const { return maxItsNewton; }
  double getEpsNewton() const { return epsNewton; }
  double getEpsAbsResNewton() const { return epsAbsResNewton; }
  double getEpsAbsIncNewton() const { return epsAbsIncNewton; }
  FILE* getOutputNewton() const { return outputNewton; }
  int getLineSearch() const { return (maxItsLS>0); }
  int getMaxItsLineSearch() const { return maxItsLS; }
  double getContractionLineSearch() const { return contractionLS; }
  double getSufficientDecreaseLineSearch() const { return sufficDecreaseLS; }
  
  void setCurrentStateForKspBinaryOutput(DistSVec<double,dim> &Q) {}
};


#ifdef TEMPLATE_FIX
#include <ImplicitEmbeddedTsDesc.C>
#endif
