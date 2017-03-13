#include "MultiGridDistSVec.h"

template <class Scalar,int dim>
MultiGridDistSVec<Scalar,dim>::
MultiGridDistSVec() {

  myVecs = NULL;
  nLevels = 0;
}

template <class Scalar,int dim>
MultiGridDistSVec<Scalar,dim>::
~MultiGridDistSVec() {

  
}

template <class Scalar,int dim>
void MultiGridDistSVec<Scalar,dim>::
free() {

  for (int i = 0; i < nLevels; ++i) {
  
    delete myVecs[i];
  }

  if (myVecs)
    delete [] myVecs;
}

template <class Scalar,int dim>
void MultiGridDistSVec<Scalar,dim>::init(MultiGridKernel<Scalar>* K) {

  pKernel = K;
  nLevels = pKernel->numLevels();

  myVecs = new DistSVec<Scalar,dim>*[nLevels];

  for (int i = 0; i < nLevels; ++i) {

    myVecs[i] = new DistSVec<Scalar,dim>(pKernel->getLevel(i)->getNodeDistInfo());
  }
}

template class MultiGridDistSVec<double, 1>;
template class MultiGridDistSVec<double, 2>;
template class MultiGridDistSVec<double, 5>;
template class MultiGridDistSVec<double, 6>;
template class MultiGridDistSVec<double, 7>;
