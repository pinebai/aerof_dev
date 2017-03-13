#ifndef _KSP_MG_PREC_H_
#define _KSP_MG_PREC_H_

#include <DistVector.h>
#include <KspPrec.h>
#include <KspSolver.h>
#include <MultiGridLevel.h>
#include <MvpMatrix.h>
#include <MultiGridSmoothingMatrix.h>
#include <SpaceOperator.h>
#include <MultiGridKernel.h>
#include <MultiGridMvpMatrix.h>
#include <MultiGridSmoothingMatrices.h>
#include <MultiGridDistSVec.h>

class Domain;
class DistGeoState;

template<class Scalar, int dim,class Scalar2 = double>
class MultiGridPrec : public KspPrec<dim, Scalar2>, public DistMat<Scalar2,dim> {

public:

  MultiGridPrec(Domain *, DistGeoState &, IoData&);
  ~MultiGridPrec();

  void initialize();

  bool isInitialized(); 

  void setup() { } 

  void setup(DistSVec<Scalar2,dim>&);

  void apply(DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);

  void getData(DistMat<Scalar2,dim>& mat);/*,
               DistSVec<Scalar2,dim>& V,
               SpaceOperator<dim>& spo,
               DistTimeState<dim>* timeState);*/

  DistMat<Scalar2,dim> &operator= (const Scalar2 s);// { return macroA[0]->operator = (s); }

  GenMat<Scalar2,dim> &operator() (int i);// { return macroA[0]->operator()(i); }

  void computeResidual(int lvl, MultiGridDistSVec<double,dim>& b,
                       MultiGridDistSVec<double,dim>& x,
                       MultiGridDistSVec<double,dim>& R);

 private:

  void cycle(int lvl, MultiGridDistSVec<double,dim>& b,
             MultiGridDistSVec<double,dim>& x);

  DistMat<Scalar2,dim>* myFineMat;

  MultiGridKernel<Scalar2>* mgKernel;

  MultiGridMvpMatrix<Scalar2,dim>* mvpMatrices;

  MultiGridSmoothingMatrices<Scalar2, dim>* smoothingMatrices;

  MultiGridDistSVec<Scalar2,dim> xold, myX, myPx, myR;  

  DistMvpMatrix<Scalar2,dim>* mvpMat;

  Domain* domain;
};

#endif
