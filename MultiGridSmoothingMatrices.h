#pragma once

#include <MultiGridKernel.h>

#include <MultiGridMvpMatrix.h>

#include <MultiGridDistSVec.h>

template <class Scalar, int dim>
class MultiGridSmoothingMatrices {

 public:

  MultiGridSmoothingMatrices(MultiGridKernel<Scalar>* pKernel,
                             typename MultiGridSmoothingMatrix<Scalar,dim>::SmoothingMode);

  ~MultiGridSmoothingMatrices();

  void acquire(DistMat<Scalar,dim>&);

  void acquire(int lvl, MultiGridMvpMatrix<Scalar,dim>& mvp);
  
  void apply(int lvl, MultiGridDistSVec<Scalar,dim>& dx, 
             MultiGridDistSVec<Scalar,dim>& R);

 private:
  
  MultiGridKernel<Scalar>* pKernel;

  int nLevels;

  int numLocSub;

  MultiGridSmoothingMatrix<Scalar,dim>*** pMatrices;
};
