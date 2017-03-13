#ifndef _MULTIGRID_SMOOTHING_MATRIX_H_
#define _MULTIGRID_SMOOTHING_MATRIX_H_

#include <Vector.h>
#include <GenMatrix.h>
#include <BlockTridiagonalMatrix.h>
#include <SparseMatrix.h>

template <class Scalar>
class MultiGridLevel;

//------------------------------------------------------------------------------

template <class Scalar, int dim>
class MultiGridSmoothingMatrix : public GenMat<Scalar,dim> {

 public:

  enum SmoothingMode { BlockJacobi, LineJacobi, RAS };

  MultiGridSmoothingMatrix(SmoothingMode m,int iSub,
                           int nn, int ne, int nBC,
                           MultiGridLevel<Scalar>*,
                           MultiGridLevel<Scalar>*);

  ~MultiGridSmoothingMatrix();
 
  MultiGridSmoothingMatrix<Scalar,dim> &operator= (const Scalar x) { a = x; return *this; }
  MultiGridSmoothingMatrix<Scalar,dim> &operator*= (const Scalar x) { a *= x; return *this; }

  double norm() {return a.norm();}

  int numNonZeroBlocks() const { return a.size(); }
  Scalar (*data())[dim*dim] { return a.data(); }

  Scalar *getElem_ii(int i) { return *(a.v + i); }
  Scalar *getElem_ij(int l) { return *(a.v + n + 2 * l); }
  Scalar *getElem_ji(int l) { return *(a.v + n + 2 * l + 1); }

  Scalar *getBcElem_ij(int l) { return *(a.v + n + 2 * numEdges + 2 * l); }
  Scalar *getBcElem_ji(int l) { return *(a.v + n + 2 * numEdges + 2 * l + 1); }
  void addContrib(int nnd, int *nd, double *K) {}

  void getData(GenMat<Scalar,dim>& mat);

  void smooth(SVec<Scalar,dim>& r, SVec<Scalar,dim>& du);

 private:

  int n;
 
  int iSub;

  int numEdges;

  SVec<Scalar,dim*dim> a;

  // Permutation indices for an LU decomposition.
  SVec<int,dim>* indices;

  BlockTridiagonalMatrix<Scalar,dim>** lineJacobiMatrices;

  SmoothingMode mySmoothingMode;

  void getDataBlockJacobi(GenMat<Scalar,dim>& mat);
  
  void smoothBlockJacobi(SVec<Scalar,dim>& r, SVec<Scalar,dim>& du);
  
  void getDataLineJacobi(GenMat<Scalar,dim>& mat);
  
  void smoothLineJacobi(SVec<Scalar,dim>& r, SVec<Scalar,dim>& du);
  
  void smoothRAS(SVec<Scalar,dim>& r, SVec<Scalar,dim>& du);

  MultiGridLevel<Scalar>* mgLevel,*mgLevelRefined;

  int numLines;

  SparseMat<Scalar,dim> *iluA;
  Vec<int>* iluJW;
};

#endif
