#include <MultiGridPrec.h>

#include <Domain.h>
#include <DistGeoState.h>
#include <MultiGridLevel.h>
#include <MatVecProd.h>
#include <DistMvpMatrix.h>

#include <IoData.h>

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
MultiGridPrec<Scalar,dim,Scalar2>::MultiGridPrec(Domain *dom, DistGeoState& distGeoState, IoData& ioData) : domain(dom), DistMat<Scalar2,dim>(dom)
{

  mgKernel = new MultiGridKernel<Scalar2>(dom, distGeoState, ioData,
                                         ioData.mg.num_multigrid_levels);

}

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::initialize() {

  mgKernel->initialize(dim,dim,0);

  smoothingMatrices = 
    new MultiGridSmoothingMatrices<Scalar2, dim>(mgKernel, 
              MultiGridSmoothingMatrix<Scalar2,dim>::RAS);

  myX.init(mgKernel);
  myR.init(mgKernel);
  myPx.init(mgKernel);
  xold.init(mgKernel);

  mvpMatrices = new MultiGridMvpMatrix<Scalar2,dim>(domain, mgKernel);

  mvpMat = new DistMvpMatrix<Scalar2,dim>(domain, mgKernel->getLevel(0)->getNodeDistInfo().subLen,
                               mgKernel->getLevel(0)->getEdges(),0);
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
MultiGridPrec<Scalar,dim,Scalar2>::~MultiGridPrec()
{
  delete mgKernel;
}

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::setup(DistSVec<Scalar2,dim>& U)
{
}

template<class Scalar, int dim, class Scalar2>
bool MultiGridPrec<Scalar,dim,Scalar2>::isInitialized() {

  return mgKernel->isInitialized();
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::
computeResidual(int lvl, MultiGridDistSVec<double,dim>& b,
                MultiGridDistSVec<double,dim>& x,
                MultiGridDistSVec<double,dim>& R) {

  if (lvl == 0)
    mgKernel->getLevel(lvl)->computeMatVecProd(*myFineMat, x(lvl), R(lvl));
  else
    mgKernel->getLevel(lvl)->computeMatVecProd((*mvpMatrices)(lvl), x(lvl), R(lvl));
  R(lvl) = b(lvl)-R(lvl);  
}

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::cycle(int lvl, MultiGridDistSVec<double,dim>& b,
                                              MultiGridDistSVec<double,dim>& x) {

  static double lvl_times[] = {0.0,0.0,0.0,0.0,0.0};
  const static int num_presmooth[] = {0,0,0,0,0,0,0};
  const static int num_postsmooth[] = {1,2,3,3,3,3,3};//3,3,3,3,3,3,3};
  static int zero_cnt = 0;
  if (lvl == 0)
    zero_cnt++;

  for (int i = 0; i < num_presmooth[lvl]; ++i) {
    computeResidual(lvl,b,x,myR);
    smoothingMatrices->apply(lvl, x, myR);
    x(lvl) += xold(lvl); 
  }
 
  if (lvl < mgKernel->numLevels()-1) {

    if ( num_presmooth[lvl] > 0)
      computeResidual(lvl,b,x,myR);
    else
      myR(lvl) = b(lvl);
    //mgKernel->Restrict(lvl+1, x(lvl), x(lvl+1)); 
    mgKernel->Restrict(lvl+1, myR(lvl), b(lvl+1));
    xold(lvl+1) = 0.0;//x(lvl+1); 
    x(lvl+1) = 0.0;
    cycle(lvl+1, b,x);

    mgKernel->Prolong(lvl+1, xold(lvl+1), x(lvl+1), x(lvl),x(lvl), 1.0, NULL);
  }

  double t = domain->getTimer()->getTime();
  for (int i = 0; i < num_postsmooth[lvl]; ++i) {
    computeResidual(lvl,b,x,myR);
    smoothingMatrices->apply(lvl, xold, myR);
    x(lvl) += xold(lvl);
  }
  lvl_times[lvl] += domain->getTimer()->getTime() - t;

  if (zero_cnt % 50 == 0)
    domain->getCommunicator()->fprintf(stdout, "Time[%d] = %e\n",lvl, lvl_times[lvl]);
  
  //computeResidual(lvl,b,x,myR);
}

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::apply(DistSVec<Scalar2,dim> & x, DistSVec<Scalar2,dim> & Px)
{
  static int cnt = 0;
  static double prec_time = 0.0;
  ++cnt; 
  double t = domain->getTimer()->getTime();
  myX(0) = x;
  myPx(0) = 0.0;//Px;
  //for (int i = 0; i < 20; ++i) {
    cycle(0, myX, myPx);
    
  //  domain->getCommunicator()->fprintf(stdout,"%e\n",myR(0).norm());;
  //}
  Px = myPx(0);
  prec_time += domain->getTimer()->getTime() - t;
  if (cnt % 50 == 0)
    domain->getCommunicator()->fprintf(stdout, "Prec time = %e\n",prec_time);
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void MultiGridPrec<Scalar,dim,Scalar2>::getData(DistMat<Scalar2,dim>& mat)/* ,
                                                DistSVec<Scalar2,dim>& V,
                                                SpaceOperator<dim>& spo,
                                                DistTimeState<dim>* timeState)*/
{

  static int cnt = 0;
  static double prec_time = 0.0;
  ++cnt; 
  double t = domain->getTimer()->getTime();
  myFineMat = &mat;
//  mgKernel->getData(mat);
  smoothingMatrices->acquire(mat);
  for (int i = 1; i < mgKernel->numLevels(); ++i) {

    if (i == 1)
      mgKernel->getLevel(i)->RestrictOperator(*mgKernel->getLevel(i-1),
                                              mat,
                                              (*mvpMatrices)(i));
   else 
      mgKernel->getLevel(i)->RestrictOperator(*mgKernel->getLevel(i-1),
                                              (*mvpMatrices)(i-1),
                                              (*mvpMatrices)(i));
    smoothingMatrices->acquire(i, (*mvpMatrices));
  }  
  prec_time += domain->getTimer()->getTime() - t;
  if (cnt % 10 == 0)
    domain->getCommunicator()->fprintf(stdout, "Prec setup time = %e\n",prec_time);
}

template<class Scalar, int dim, class Scalar2>
DistMat<Scalar2,dim>& MultiGridPrec<Scalar,dim,Scalar2>::
operator= (const Scalar2 s) {

  return (*mvpMat = s);
}

template<class Scalar, int dim, class Scalar2>
GenMat<Scalar2,dim>& MultiGridPrec<Scalar,dim,Scalar2>::
operator() (int i) {

  return (*mvpMat)(i);//mgKernel->getFineMatrix()(i);
}

//------------------------------------------------------------------------------
#define INSTANTIATION_HELPER(T,dim) \
    template class MultiGridPrec<T, dim, double>;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,5);
INSTANTIATION_HELPER(float,6);
INSTANTIATION_HELPER(float,7);

INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,5);
INSTANTIATION_HELPER(double,6);
INSTANTIATION_HELPER(double,7);
