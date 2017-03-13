#pragma once

#include <MultiGridKernel.h>

template <class Scalar,int neq>
class MultiGridMvpMatrix {

 public:

  MultiGridMvpMatrix(Domain*,MultiGridKernel<Scalar>* pKernel);

  ~MultiGridMvpMatrix();

  DistMvpMatrix<Scalar,neq>& operator()(int i) { return *myA[i]; }

 private:

  DistMvpMatrix<Scalar,neq>** myA;

  int nLevels;

  MultiGridKernel<Scalar>* pKernel;
};
