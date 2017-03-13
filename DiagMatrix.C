#include <DiagMatrix.h>

#include <BcDef.h>
#include <SparseMatrix.h>
#include <DenseMatrixOps.h>

#include <cstdlib>

//------------------------------------------------------------------------------

template<class Scalar, int dim>
DiagMat<Scalar,dim>::DiagMat(Type t, int nn, int *ndType) : a(nn) 
{ 

  type = t; 
  n = nn; 
  nodeType = ndType; 

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void DiagMat<Scalar,dim>::addContrib(int nnd, int *nd, double *K)
{

  int colSize = dim*nnd;

  if (type == DENSE) {
    for (int i=0; i<nnd; ++i) //HB: need testing ...
      for (int k=0; k<dim; ++k)
        for (int l=0; l<dim; ++l)
          a[nd[i]][k*dim+l] += K[(dim*i+k)*colSize + dim*i + l];
  }
  else if (type == DIAGONAL) {
    for (int i = 0; i < nnd; ++i)
      for (int k=0; k<dim; ++k) 
	a[nd[i]][k*dim + k] += K[(dim*i+k)*colSize + dim*i + k];
  }

}
//------------------------------------------------------------------------------

template<class Scalar, int dim>
void DiagMat<Scalar,dim>::invert()
{

  bool invert = true;
  int ierr;
  if (type == DENSE) {
    for (int i=0; i<n; ++i) {

      invert = true;
      for(int j = 0; j < dim; j++)  {
        int k = dim*j + j;
        if (a[i][k] == 0)
          invert = false;
      }
      
      if (invert)
        ierr = invertDenseMatrix<Scalar,dim>(a[i]);
      else
        ierr = 0;

      if (ierr) {
	fprintf(stderr, "*** Error: zero pivot at row %d during diagonal inversion\n", i);
	exit(1);
      }
    }
  }
  else if (type == DIAGONAL) {
    for (int i=0; i<n; ++i) {
      if (nodeType && nodeType[i] != BC_INTERNAL)
        for (int k=0; k<dim; ++k)
          a[i][k*dim + k] = 0.0;
      else
        for (int k=0; k<dim; ++k)
          a[i][k*dim + k] = 1.0 / a[i][k*dim + k];
    }
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void DiagMat<Scalar,dim>::apply(SVec<double,dim> &y, SVec<double,dim> &x)
{

  if (type == DENSE) {
    for (int i=0; i<n; ++i)
      DenseMatrixOp<Scalar,dim,dim*dim>::applyToVector(a.v, i, y.v, i, x.v, i);
  }
  else if (type == DIAGONAL) {
    for (int i=0; i<n; ++i)
      for (int k=0; k<dim; ++k)
	x[i][k] = a[i][k*dim + k] * y[i][k];
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<class MatScal>
void DiagMat<Scalar,dim>::getData(GenMat<MatScal,dim> &B)
{

  for (int i=0; i<n; ++i) {

    MatScal *b = B.getElem_ii(i);

    for (int k=0; k<dim*dim; ++k) 
      a[i][k] = b[k];

  }

}

//------------------------------------------------------------------------------
