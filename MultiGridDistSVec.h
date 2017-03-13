#pragma once

#include <MultiGridKernel.h>

template <class Scalar,int dim>
class MultiGridDistSVec {

 public:

  MultiGridDistSVec();

  ~MultiGridDistSVec();

  void init(MultiGridKernel<Scalar>* pKernel);

  void free();

  DistSVec<Scalar,dim>& operator()(int i) { return *myVecs[i]; }

 private:

  DistSVec<Scalar,dim>** myVecs;

  int nLevels;

  MultiGridKernel<Scalar>* pKernel;
};
