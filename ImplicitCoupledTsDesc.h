#ifndef _IMPLICIT_COUPLED_TS_DESC_H_
#define _IMPLICIT_COUPLED_TS_DESC_H_

#include <ImplicitTsDesc.h>
#include <KspPrec.h>
#include <IoData.h>

class Domain;
template<int dim, int neq> class MatVecProd;

#ifndef _KSPSLVR_TMPL_
#define _KSPSLVR_TMPL_
template<class VecType, class MvpOp, class PrecOp, class IoOp, class ScalarT = double> class KspSolver;
#endif

//------------------------------------------------------------------------------

template<int dim>
class ImplicitCoupledTsDesc : public ImplicitTsDesc<dim> {

protected:
  MatVecProd<dim,dim> *mvp;
  KspPrec<dim> *pc;
  KspSolver<DistSVec<double,dim>, MatVecProd<dim,dim>, KspPrec<dim>, Communicator> *ksp;
  
#ifdef MVP_CHECK
  MatVecProd<dim,dim> *mvpfd1;
#endif

  template<int neq>
  KspSolver<DistSVec<double,neq>, MatVecProd<dim,neq>, KspPrec<neq>,
    Communicator> *createKrylovSolver(const DistInfo &, KspData &, MatVecProd<dim,neq> *,
              KspPrec<neq> *, Communicator *);



public:

  ImplicitCoupledTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitCoupledTsDesc();

  void computeJacobian(int, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void setOperators(DistSVec<double,dim> &);
  int solveLinearSystem(int, DistSVec<double,dim> &, DistSVec<double,dim> &);
  
// Included (MB)
  void rstVarImplicitCoupledTsDesc(IoData &);

  DistMat<double,dim>* GetJacobian() { 
    MatVecProdH1<dim,double,dim>* mvph1 = 
      dynamic_cast<MatVecProdH1<dim,double,dim>*>(this->mvp);
    if (mvph1) 
      return dynamic_cast<DistMat<double,dim>*>(mvph1);
    MatVecProdH2<dim,double,dim>* mvph2 = 
      dynamic_cast<MatVecProdH2<dim,double,dim>*>(this->mvp);
    if (mvph2) 
      return dynamic_cast<DistMat<double,dim>*>(mvph2);
    return NULL;
  }
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitCoupledTsDesc.C>
#endif

#endif
