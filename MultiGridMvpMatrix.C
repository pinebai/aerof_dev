#include "MultiGridMvpMatrix.h"

template <class Scalar,int neq>
MultiGridMvpMatrix<Scalar,neq>::
MultiGridMvpMatrix(Domain* domain,
                   MultiGridKernel<Scalar>* pKernel) : pKernel(pKernel) {

  nLevels = pKernel->numLevels();

  myA = new DistMvpMatrix<Scalar,neq>*[nLevels];

  for (int i = 1; i < nLevels; ++i) {

    myA[i] = new DistMvpMatrix<Scalar,neq>(domain, 
        pKernel->getLevel(i)->getNodeDistInfo().subLen,
        pKernel->getLevel(i)->getEdges(),(FaceSet**)NULL);
  }
}

template <class Scalar,int dim>
MultiGridMvpMatrix<Scalar,dim>::
~MultiGridMvpMatrix() {

  for (int i = 1; i < nLevels; ++i) {
  
    delete myA[i];
  }

  delete [] myA;
}

template class MultiGridMvpMatrix<double,1>;
template class MultiGridMvpMatrix<double,2>;
template class MultiGridMvpMatrix<double,5>;
template class MultiGridMvpMatrix<double,6>;
template class MultiGridMvpMatrix<double,7>;
