#pragma once

#include <DenseMatrixOps.h>

template <class Scalar,int dim>
class BlockTridiagonalMatrix {

 public: 

  BlockTridiagonalMatrix(int n);

  ~BlockTridiagonalMatrix();

  void computeLU(const BlockTridiagonalMatrix&);

  // Compute the LU decomposition in place
  void computeLU();

  Scalar (*get(int i,int j))[dim*dim];

  void solveLU(Scalar (*b)[dim],Scalar (*x)[dim]);

 private:

  Scalar (*diag)[dim*dim];
  Scalar (*super)[dim*dim];
  Scalar (*sub)[dim*dim];

  int (*indices)[dim];

  int N;
};
