#include <BlockTridiagonalMatrix.h>

#include <memory.h>

template <class Scalar,int dim>
BlockTridiagonalMatrix<Scalar,dim>::BlockTridiagonalMatrix(int n) : N(n) {

  diag = new Scalar[N][dim*dim];
  sub = new Scalar[N-1][dim*dim];
  super = new Scalar[N-1][dim*dim];

}

template <class Scalar,int dim>
BlockTridiagonalMatrix<Scalar,dim>::~BlockTridiagonalMatrix() {
 
  delete [] diag;
  delete [] sub;
  delete [] super;
}

template <class Scalar,int dim>
void BlockTridiagonalMatrix<Scalar,dim>::
computeLU(const BlockTridiagonalMatrix& B) {

  indices = new int[N][dim];

  Scalar tmp1[dim],tmp2[dim*dim];
  for (int i = 0; i < N; ++i) {

    if (this != &B)
      memcpy(diag[i], B.diag[i], sizeof(Scalar)*dim*dim);
    if (i > 0) {
     
      DenseMatrixOp<Scalar,dim,dim*dim>::applyToDenseMatrix(sub,i-1, super, i-1, &tmp2,0);
      for (int j = 0; j < dim*dim; ++j)
        diag[i][j] -= tmp2[j];
    }

    DenseMatrixOp<Scalar,0,0>::ludec(diag[i], indices[i], 1.0, dim);
    if (i < N-1) { 
    
      for (int j = 0; j < dim; ++j) {
        for (int k = 0; k < dim; ++k)
          tmp1[k] = B.super[i][k*dim+j];
        DenseMatrixOp<Scalar,0,0>::ludfdbksb(diag[i], indices[i], tmp1, dim); 
        for (int k = 0; k < dim; ++k)
          super[i][k*dim+j] = tmp1[k];
      }
    }
   
    if (i > 0) {

      if (this != &B)
        memcpy(sub[i-1], B.sub[i-1], sizeof(Scalar)*dim*dim);
    }
  }
}

template <class Scalar,int dim>
void BlockTridiagonalMatrix<Scalar,dim>::computeLU() {

  computeLU(*this);
}

template <class Scalar,int dim>
Scalar (*BlockTridiagonalMatrix<Scalar,dim>::get(int i,int j))[dim*dim] {

  if (i == j) 
    return &diag[i];
  else if (i == j+1)
    return &sub[j];
  else if (i+1 == j)
    return &super[i];
  else
    return NULL;
}

template <class Scalar,int dim>
void BlockTridiagonalMatrix<Scalar,dim>::solveLU(Scalar (*b)[dim],Scalar (*x)[dim]) {

  Scalar tmp[dim];

  // forward substitution
  for (int i = 0; i < N; ++i) {

    if (i == 0)
      memcpy(x[i],b[i],sizeof(Scalar)*dim);
    else {
      DenseMatrixOp<Scalar,dim,dim*dim>::applyToVector(sub,i-1,x,i-1,x,i);
      for (int j = 0; j < dim; ++j)
        x[i][j] = b[i][j] - x[i][j];
    }
        
    DenseMatrixOp<Scalar,0,0>::ludfdbksb(diag[i], indices[i], x[i], dim);
  }

  // back substitute
  for (int i = N-2; i >= 0; --i) {
 
    DenseMatrixOp<Scalar,dim,dim*dim>::applyToVector(super,i, x, i+1, &tmp,0);
    for (int j = 0; j < dim; ++j)
      x[i][j] -= tmp[j];
  }
}

#define INST_HELPER(Sc,d) template class BlockTridiagonalMatrix<Sc,d>;

INST_HELPER(double,1)
INST_HELPER(double,2)
INST_HELPER(double,3)
INST_HELPER(double,5)
INST_HELPER(double,6)
INST_HELPER(double,7)

