#ifndef _IMPLICIT_SEG_TS_DESC_H_
#define _IMPLICIT_SEG_TS_DESC_H_

#include <ImplicitTsDesc.h>
#include <KspPrec.h>

template<int dim, int neq> class MatVecProd;

#ifndef _KSPSLVR_TMPL_
#define _KSPSLVR_TMPL_
template<class VecType, class MvpOp, class PrecOp, class IoOp, class ScalarT = double> class KspSolver;
#endif

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
class ImplicitSegTsDesc : public ImplicitTsDesc<dim> {

  DistSVec<double,neq1> b1, dQ1;
  DistSVec<double,neq2> b2, dQ2;


  MatVecProd<dim,neq1> *mvp1;
  MatVecProd<dim,neq2> *mvp2;

#ifdef MVP_CHECK
  MatVecProd<dim,dim> *mvpfd;
#endif

  KspPrec<neq1> *pc1;
  KspPrec<neq2> *pc2;

  KspSolver<DistSVec<double,neq1>, MatVecProd<dim,neq1>, KspPrec<neq1>, Communicator> *ksp1;
  KspSolver<DistSVec<double,neq2>, MatVecProd<dim,neq2>, KspPrec<neq2>, Communicator> *ksp2;

private:

  SpaceOperator<dim> *createSpaceOperator1(IoData &, SpaceOperator<dim> *);
  SpaceOperator<dim> *createSpaceOperator2(IoData &, SpaceOperator<dim> *);

  template<int neq>
  void setOperator(MatVecProd<dim,neq> *, KspPrec<neq> *, DistSVec<double,dim> &, SpaceOperator<dim> *);

protected:

  SpaceOperator<dim> *spaceOp1;
  SpaceOperator<dim> *spaceOp2;

public:

  ImplicitSegTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitSegTsDesc();

  void computeJacobian(int, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void setOperators(DistSVec<double,dim> &);

  void rstVarImplicitSegTsDesc(IoData &);

  int solveLinearSystem(int, DistSVec<double,dim> &, DistSVec<double,dim> &);

  DistMat<double,neq1>* GetJacobian1() {
    MatVecProdH1<dim,double,neq1>* mvph1 =
      dynamic_cast<MatVecProdH1<dim,double,neq1>*>(this->mvp1);
    if (mvph1)
      return dynamic_cast<DistMat<double,neq1>*>(mvph1);
    MatVecProdH2<dim,double,neq1>* mvph2 =
      dynamic_cast<MatVecProdH2<dim,double,neq1>*>(this->mvp1);
    if (mvph2)
      return dynamic_cast<DistMat<double,neq1>*>(mvph2);
    return NULL;
  }

  DistMat<double,neq2>* GetJacobian2() {
    MatVecProdH1<dim,double,neq2>* mvph1 =
      dynamic_cast<MatVecProdH1<dim,double,neq2>*>(this->mvp2);
    if (mvph1)
      return dynamic_cast<DistMat<double,neq2>*>(mvph1);
    MatVecProdH2<dim,double,neq2>* mvph2 =
      dynamic_cast<MatVecProdH2<dim,double,neq2>*>(this->mvp2);
    if (mvph2)
      return dynamic_cast<DistMat<double,neq2>*>(mvph2);
    return NULL;
  }

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitSegTsDesc.C>
#endif

#endif
