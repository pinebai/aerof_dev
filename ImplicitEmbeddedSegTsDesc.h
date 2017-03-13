#pragma once

#include <ImplicitEmbeddedTsDesc.h>

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

template<int dim, int neq1, int neq2>
class ImplicitEmbeddedSegTsDesc : public ImplicitEmbeddedTsDesc<dim> {
   
  DistEmbeddedVec<double,neq1> embeddedB1,embeddeddQ1;
  DistEmbeddedVec<double,neq2> embeddedB2,embeddeddQ2;

  MatVecProd<dim,neq1> *mvp1;
  MatVecProd<dim,neq2> *mvp2;

  KspPrec<neq1> *pc1;
  KspPrec<neq2> *pc2;

  KspSolver<DistEmbeddedVec<double,neq1>, MatVecProd<dim,neq1>, KspPrec<neq1>, Communicator> *ksp1;
  KspSolver<DistEmbeddedVec<double,neq2>, MatVecProd<dim,neq2>, KspPrec<neq2>, Communicator> *ksp2;
	
private:

  SpaceOperator<dim> *createSpaceOperator1(IoData &, SpaceOperator<dim> *);
  SpaceOperator<dim> *createSpaceOperator2(IoData &, SpaceOperator<dim> *);

  template<int neq>
  void setOperator(MatVecProd<dim,neq> *, KspPrec<neq> *, DistSVec<double,dim> &);

protected:

  SpaceOperator<dim> *spaceOp1;
  SpaceOperator<dim> *spaceOp2;

public:
  
  ImplicitEmbeddedSegTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom);

  ~ImplicitEmbeddedSegTsDesc();

  void computeJacobian(int it, DistSVec<double,dim> &Q,
		       DistSVec<double,dim> &F);

  void setOperators(DistSVec<double,dim> &Q);

  int solveLinearSystem(int it, DistSVec<double,dim> &b,
			DistSVec<double,dim> &dQ);
};


#ifdef TEMPLATE_FIX
#include <ImplicitEmbeddedSegTsDesc.C>
#endif
