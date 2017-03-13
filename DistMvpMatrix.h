#ifndef _DIST_MVP_MATRIX_H_
#define _DIST_MVP_MATRIX_H_

#include <DistMatrix.h>
#include <GenMatrix.h>
#include <Communicator.h>
#include <MvpMatrix.h>

class Domain;

template <class Scalar,int dim>
class DistMvpMatrix : public DistMat<Scalar,dim> {

 public:

  DistMvpMatrix(Domain*, const int* nn,
                EdgeSet** edges, FaceSet** faces);

  ~DistMvpMatrix();

  virtual DistMat<Scalar,dim> &operator= (const Scalar);

  virtual GenMat<Scalar,dim> &operator() (int);

//  void apply(DistSVec<Scalar,dim>& p, DistSVec<Scalar,dim>& prod);

 private:

  Domain* domain;
 
  MvpMat<Scalar,dim> **A;

  EdgeSet** edges;

  FaceSet** faces;

  int* numNodes;

  int numLocSub;
};

template <class Scalar,int dim>
DistMvpMatrix<Scalar,dim>::DistMvpMatrix(Domain* d, const int* nn,
                                         EdgeSet** edges,
                                         FaceSet** faces) : DistMat<Scalar,dim>(d) , domain(d),
                                                            edges(edges), faces(faces) {
 
  numLocSub = d->getNumLocSub();
 
  A = new MvpMat<Scalar,dim> *[d->getNumLocSub()];

  numNodes = new int[d->getNumLocSub()];

  memcpy(numNodes,nn,sizeof(int)*d->getNumLocSub());

#pragma omp parallel for
  for (int iSub = 0; iSub < d->getNumLocSub(); ++iSub) {
  
//    std::cout << "Num nodes = " << numNodes[iSub] << " num edges = " << edges[iSub]->size() << std::endl;

    A[iSub] = new MvpMat<Scalar,dim>(numNodes[iSub],edges[iSub]->size());
  }
}

template <class Scalar,int dim>
DistMvpMatrix<Scalar,dim>::~DistMvpMatrix() {

  delete [] numNodes;

#pragma omp parallel for
  for (int iSub = 0; iSub < domain->getNumLocSub(); ++iSub) { 
    
    delete A[iSub];
  }

  delete [] A;
}

template <class Scalar,int dim>
DistMat<Scalar,dim> & DistMvpMatrix<Scalar,dim>::operator= (const Scalar s) {

#pragma omp parallel for 
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
 
    (*A[iSub]) = s;
  }

  return *this;
}

template <class Scalar,int dim>
GenMat<Scalar,dim>& DistMvpMatrix<Scalar,dim>::operator() (int iSub) {

  return *A[iSub];
}
/*
template <class Scalar,int dim> 
void DistMvpMatrix<Scalar,dim>::apply(DistSVec<Scalar,dim>& p, DistSVec<Scalar,dim>& prod,
                                      CommPattern* vecPat) {


}
*/
#endif // _DIST_MVP_MATRIX_H_
