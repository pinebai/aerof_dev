#include <MultiGridSmoothingMatrices.h>

template <class Scalar,int dim>
MultiGridSmoothingMatrices<Scalar,dim>::
MultiGridSmoothingMatrices(MultiGridKernel<Scalar>* pKernel,
                           typename
MultiGridSmoothingMatrix<Scalar,dim>::SmoothingMode smoothingMode)  :
pKernel(pKernel) {

  nLevels = pKernel->numLevels();

  DistInfo& finestInfo = pKernel->getLevel(0)->getNodeDistInfo();

  numLocSub = finestInfo.numLocSub;

  pMatrices = new MultiGridSmoothingMatrix<Scalar,dim>**[nLevels];

  for (int i = 0; i < nLevels; ++i) {

    pMatrices[i] = new MultiGridSmoothingMatrix<Scalar,dim>*[finestInfo.numLocSub];
    
    MultiGridLevel<Scalar>* pL = NULL;
    MultiGridLevel<Scalar>* pR = pKernel->getLevel(i);
    if (i < nLevels-1)
      pL = pKernel->getLevel(i+1);    

#pragma omp parallel for
    for (int iSub = 0; iSub < finestInfo.numLocSub; ++iSub) {

      int nn = pKernel->getLevel(i)->getNodeDistInfo().subSize(iSub);
      int ne = pKernel->getLevel(i)->getEdges()[iSub]->size();

      typename MultiGridSmoothingMatrix<Scalar,dim>::SmoothingMode s = 
        smoothingMode;
  
      if (i == nLevels-1 && s == MultiGridSmoothingMatrix<Scalar,dim>::LineJacobi)
        s = MultiGridSmoothingMatrix<Scalar,dim>::RAS;    
  
      pMatrices[i][iSub] = 
        new MultiGridSmoothingMatrix<Scalar,dim>(s,iSub, nn,ne,0,
                                                 pL, pR);
                                                 
    }
  } 
}

template <class Scalar, int dim>
MultiGridSmoothingMatrices<Scalar,dim>::~MultiGridSmoothingMatrices() {

  DistInfo& finestInfo = pKernel->getLevel(0)->getNodeDistInfo();
  for (int i = 0; i < nLevels; ++i) {

#pragma omp parallel for
    for (int iSub = 0; iSub < finestInfo.numLocSub; ++iSub) {
      delete pMatrices[i][iSub];
    }
    delete [] pMatrices[i];
  }
  
  delete [] pMatrices;
}

template <class Scalar, int dim>
void MultiGridSmoothingMatrices<Scalar,dim>::
acquire(DistMat<Scalar,dim>& mvp) {

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    pMatrices[0][iSub]->getData(mvp(iSub));
  }

}


template <class Scalar, int dim>
void MultiGridSmoothingMatrices<Scalar,dim>::
acquire(int lvl, MultiGridMvpMatrix<Scalar,dim>& mvp) {

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    pMatrices[lvl][iSub]->getData(mvp(lvl)(iSub));
  }

}

template <class Scalar, int dim>
void MultiGridSmoothingMatrices<Scalar,dim>::
apply(int lvl, MultiGridDistSVec<Scalar,dim>& dx, 
      MultiGridDistSVec<Scalar,dim>& R) {

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    pMatrices[lvl][iSub]->smooth(R(lvl)(iSub), dx(lvl)(iSub));
  }

  pKernel->getLevel(lvl)->assemble(dx(lvl));
}


template class MultiGridSmoothingMatrices<double,1>;
template class MultiGridSmoothingMatrices<double,2>;
template class MultiGridSmoothingMatrices<double,5>;
template class MultiGridSmoothingMatrices<double,6>;
template class MultiGridSmoothingMatrices<double,7>;

