#include <cstdlib>

#ifdef OLD_STL
#include <defalloc.h>
#include <algo.h>
#else
#include <algorithm>
using std::stable_sort;
#endif

#include <RectangularSparseMatrix.h>

#include <BcDef.h>
#include <Edge.h>
#include <Vector.h>
#include <ResizeArray.h>
#include <DenseMatrixOps.h>

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
RectangularSparseMat<Scalar,dim,dim2>::~RectangularSparseMat()
{ 

  if (nodeRenum) delete nodeRenum;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
int RectangularSparseMat<Scalar,dim,dim2>::find(int i, int j) 
{

  if (i >= n || i < 0) {
    fprintf(stderr, "*** Error: out of bound tet: %d vs %d\n", i, n);
    exit(-1);
  }

  for (int ii = ia[i]; ii < ia[i+1]; ++ii)
    if (ja[ii] == j) return ii;

  fprintf(stderr, "*** Error: connection not found in RectangularSparseMat: %d %d %d %d\n", i, j, ia[i], ja[ia[i]]);

  return -1;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::createPointers(EdgeSet &edges)
{

  ptr_ii.resize(n);
  ptr_ij.resize(edges.size());
  ptr_ji.resize(edges.size());

  // construction of ptr_ii

  int i, j, k;
  for (i=0; i<n; ++i) {

    int inew = (nodeRenum) ? nodeRenum->renum[i] : i;

    for (k=ia[inew]; k<ia[inew+1]; ++k) {
      if (ja[k] == inew) { ptr_ii[i] = k; break; }
    }

  }

  // construction of ptr_ij and ptr_ji

  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    i = (nodeRenum) ? nodeRenum->renum[ edgePtr[l][0] ] : edgePtr[l][0];
    j = (nodeRenum) ? nodeRenum->renum[ edgePtr[l][1] ] : edgePtr[l][1];

    // ptr_ji[l] points to a_ji

    for (k=ia[j]; k<ia[j+1]; ++k) {
      if (ja[k] == i) { ptr_ji[l] = k; break; }
    }

    // ptr_ij[l] points to a_ij

    for (k=ia[i]; k<ia[i+1]; ++k) {
      if (ja[k] == j) { ptr_ij[l] = k; break; }
    }
 
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::addContrib(int nNd, int *ndList, double *C)
{

  int i, j, k, l;
  int cSize = nNd*dim;

  for (i = 0; i < nNd; ++i) {
    int iNode = (nodeRenum) ? nodeRenum->renum[ ndList[i] ] : ndList[i];
    if (nodeType && nodeType[ ndList[i] ] != BC_INTERNAL) {
      int indx = find(iNode, iNode);
      for (k = 0; k < dim2; ++k)
        for (l = 0; l < dim; ++l)
          a[indx][dim*k+l] = (k == l) ? 1.0 : 0.0;
      continue;
    }

    for (j = 0; j < nNd; ++j) {
      if (nodeType && nodeType[ ndList[j] ] != BC_INTERNAL) continue;
      int jNode = (nodeRenum) ? nodeRenum->renum[ ndList[j] ] : ndList[j];
      int indx = find(iNode, jNode);
      for (k = 0; k < dim2; ++k)
        for (l = 0; l < dim; ++l) {
          a[indx][dim*k+l] += C[(i*dim2+k)*cSize + j*dim + l];
        }
    }
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::addContrib(int row, int col, double *C)
{

  int j, k, l;

  int iNode = (nodeRenum) ? nodeRenum->renum[ row ] : row;
//  fprintf(stderr, "iNode = %d, row = %d\n", iNode, row);
  if (nodeType && nodeType[ row ] != BC_INTERNAL) {
    int indx = find(iNode, iNode);
    for (k = 0; k < dim2; ++k)
      for (l = 0; l < dim; ++l)
        a[indx][dim*k+l] = (k == l) ? 1.0 : 0.0;
    return;
  }

  if (nodeType && nodeType[ col ] != BC_INTERNAL) return;
  int jNode = (nodeRenum) ? nodeRenum->renum[ col ] : col;
//  fprintf(stderr, "jNode = %d, col = %d\n", jNode, col);
  int indx = find(iNode, jNode);
  for (k = 0; k < dim2; ++k)
    for (l = 0; l < dim; ++l) {
      a[indx][dim*k+l] += C[k*dim + l];
    }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::addContrib(int row, int col, double C)
{

  int j, k, l;

  int iNode = (nodeRenum) ? nodeRenum->renum[ row ] : row;
//  fprintf(stderr, "iNode = %d, row = %d\n", iNode, row);
  if (nodeType && nodeType[ row ] != BC_INTERNAL) {
    int indx = find(iNode, iNode);
    for (k = 0; k < dim2; ++k)
      for (l = 0; l < dim; ++l)
        a[indx][dim*k+l] = (k == l) ? 1.0 : 0.0;
    return;
  }

  if (nodeType && nodeType[ col ] != BC_INTERNAL) return;
  int jNode = (nodeRenum) ? nodeRenum->renum[ col ] : col;
//  fprintf(stderr, "jNode = %d, col = %d\n", jNode, col);
  int indx = find(iNode, jNode);
  for (k = 0; k < dim2; ++k)
    for (l = 0; l < dim; ++l) {
      a[indx][dim*k+l] += (k == l) ? C : 0.0;
    }

}

//------------------------------------------------------------------------------


template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::printFirstElementIn_a()
{
  fprintf(stderr, "a[0] = %e\n", a.v[0][0]);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::apply(SVec<double,dim> &x, SVec<double,dim2> &Ax, int *ndType)
{
  Ax = 0.0;
  applyAndAdd(x,Ax,ndType);
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::applyAndAdd(SVec<double,dim> &x, SVec<double,dim2> &Ax, int *ndType)
{
  int i, l;
  for (i=0; i<n; ++i) {
    for (int k=ia[i]; k<ia[i+1]; ++k)
      RectangularDenseMatrixOp<Scalar,dim,dim2,dim*dim2>::applyAndAddToVector(a.v, k, x.v, ja[k], Ax.v, i);  // Gives Ax.v = a.v*x.v   
  }

  if (ndType) {
    for (i=0; i<n; ++i)
      if (ndType[i] != BC_INTERNAL)
        for (l=0; l<dim2; ++l) Ax[i][l] = 0.0;
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::applyTranspose(SVec<double,dim2> &x, SVec<double,dim> &ATx, int *ndType)
{
  ATx = 0.0;
  applyTransposeAndAdd(x,ATx,ndType);
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::applyTransposeAndAdd(SVec<double,dim2> &x, SVec<double,dim> &ATx, int *ndType)
{
  int i, l;
  for (i=0; i<n; ++i) {
    for (int k=ia[i]; k<ia[i+1]; ++k)
      RectangularDenseMatrixOp<Scalar,dim,dim2,dim*dim2>::applyTransposeAndAddToVector(a.v, k, x.v, i, ATx.v, ja[k]);  // Gives ATx.v = transpose(a).v*x.v   
  }

  if (ndType) {
    for (i=0; i<n; ++i)
      if (ndType[i] != BC_INTERNAL)
        for (l=0; l<dim; ++l) ATx[i][l] = 0.0;
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::apply(SVec<double,dim> &x, Vec<Vec3D> &Ax, int *ndType)
{
  Ax = 0.0;
  applyAndAdd(x,Ax,ndType);
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::applyAndAdd(SVec<double,dim> &x, Vec<Vec3D> &Ax, int *ndType)
{
  if(dim2 != 3) { fprintf(stderr, " *** Error: dim2 must be 3\n");  exit(-1); }
  int i, l;
  for (i=0; i<n; ++i) {
    for (int k=ia[i]; k<ia[i+1]; ++k) 
      RectangularDenseMatrixOp<Scalar,dim,dim2,dim*dim2>::applyAndAddToVector(a.v, k, x.v, ja[k], Ax.v, i);  // Gives Ax.v = a.v*x.v   
  }

  if (ndType) {
    for (i=0; i<n; ++i)
      if (ndType[i] != BC_INTERNAL)
        for (l=0; l<dim2; ++l) Ax[i][l] = 0.0;
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::apply(Vec<Vec3D> &x, SVec<double,dim2> &Ax, int *ndType)
{
  Ax = 0.0;
  applyAndAdd(x,Ax,ndType);
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::apply(Vec<double> &x, SVec<double,dim2> &Ax, int *ndType)
{
  Ax = 0.0;
  applyAndAdd(x,Ax,ndType);
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::applyAndAdd(Vec<double> &x, SVec<double,dim2> &Ax, int *ndType)
{
  if(dim != 1) { fprintf(stderr, " *** Error: dim must be 1\n");  exit(-1); }
  int i, l;
  for (i=0; i<n; ++i) {
    for (int k=ia[i]; k<ia[i+1]; ++k)
      RectangularDenseMatrixOp<Scalar,dim,dim2,dim*dim2>::applyAndAddToVector(a.v, k, x.v, ja[k], Ax.v, i);  // Gives Ax.v = a.v*x.v   
  }
  if (ndType) {
    for (i=0; i<n; ++i)
      if (ndType[i] != BC_INTERNAL)
        for (l=0; l<dim2; ++l) Ax[i][l] = 0.0;
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::applyAndAdd(Vec<Vec3D> &x, SVec<double,dim2> &Ax, int *ndType)
{
  if(dim != 3) { fprintf(stderr, " *** Error: dim must be 3\n");  exit(-1); }
  int i, l;
  for (i=0; i<n; ++i) {
    for (int k=ia[i]; k<ia[i+1]; ++k) 
      RectangularDenseMatrixOp<Scalar,dim,dim2,dim*dim2>::applyAndAddToVector(a.v, k, x.v, ja[k], Ax.v, i);  // Gives Ax.v = a.v*x.v   
  }
  if (ndType) {
    for (i=0; i<n; ++i)
      if (ndType[i] != BC_INTERNAL)
        for (l=0; l<dim2; ++l) Ax[i][l] = 0.0;
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::applyTranspose(SVec<double,dim2> &x, Vec<double> &ATx, int *ndType)
{
  ATx = 0.0;
  applyTransposeAndAdd(x,ATx,ndType);
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::applyTransposeAndAdd(SVec<double,dim2> &x, Vec<double> &ATx, int *ndType)
{
  if(dim != 1) { fprintf(stderr, " *** Error: dim must be 1\n");  exit(-1); }
  int i, l;
  for (i=0; i<n; ++i) {
    for (int k=ia[i]; k<ia[i+1]; ++k)
      RectangularDenseMatrixOp<Scalar,dim,dim2,dim*dim2>::applyTransposeAndAddToVector(a.v, k, x.v, i, ATx.v, ja[k]);  // Gives Ax.v = a.v*x.v   
  }
  if (ndType) {
    for (i=0; i<n; ++i)
      if (ndType[i] != BC_INTERNAL)
        ATx[i] = 0.0;
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::applyTranspose(SVec<double,dim2> &x, Vec<Vec3D> &ATx, int *ndType)
{
  ATx = 0.0;
  applyTransposeAndAdd(x,ATx,ndType);
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::applyTransposeAndAdd(SVec<double,dim2> &x, Vec<Vec3D> &ATx, int *ndType)
{

  if(dim != 3) { fprintf(stderr, " *** Error: dim must be 3\n");  exit(-1); }
  int i, l;
  for (i=0; i<n; ++i) {
    for (int k=ia[i]; k<ia[i+1]; ++k)
      RectangularDenseMatrixOp<Scalar,dim,dim2,dim*dim2>::applyTransposeAndAddToVector(a.v, k, x.v, i, ATx.v, ja[k]);  // Gives ATx.v = transpose(a).v*x.v   
  }
  if (ndType) {
    for (i=0; i<n; ++i)
      if (ndType[i] != BC_INTERNAL)
        for (l=0; l<dim; ++l) ATx[i][l] = 0.0;
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::applyTranspose(Vec<Vec3D> &x, SVec<double,dim> &ATx, int *ndType)
{
  ATx = 0.0;
  applyTransposeAndAdd(x,ATx,ndType);
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::applyTransposeAndAdd(Vec<Vec3D> &x, SVec<double,dim> &ATx, int *ndType)
{

  if(dim2 != 3) { fprintf(stderr, " *** Error: dim2 must be 3\n");  exit(-1); }
  int i, l;
  for (i=0; i<n; ++i) {
    for (int k=ia[i]; k<ia[i+1]; ++k)
      RectangularDenseMatrixOp<Scalar,dim,dim2,dim*dim2>::applyTransposeAndAddToVector(a.v, k, x.v, i, ATx.v, ja[k]);  // Gives ATx.v = transpose(a).v*x.v   
  }

  if (ndType) {
    for (i=0; i<n; ++i)
      if (ndType[i] != BC_INTERNAL)
        for (l=0; l<dim; ++l) ATx[i][l] = 0.0;
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::apply(SVec<double,dim> &x, Vec<double> &Ax, int *ndType)
{
  Ax = 0.0;
  applyAndAdd(x,Ax,ndType);
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::applyAndAdd(SVec<double,dim> &x, Vec<double> &Ax, int *ndType)
{
  if(dim2 != 1) { fprintf(stderr, " *** Error: dim2 must be 1\n");  exit(-1); }
  int i, l;
  for (i=0; i<n; ++i) {
    for (int k=ia[i]; k<ia[i+1]; ++k) 
      RectangularDenseMatrixOp<Scalar,dim,dim2,dim*dim2>::applyAndAddToVector(a.v, k, x.v, ja[k], Ax.v, i);  // Gives Ax.v = a.v*x.v   
  }

  if (ndType) {
    for (i=0; i<n; ++i)
      if (ndType[i] != BC_INTERNAL)
        Ax[i] = 0.0;
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::applyTranspose(Vec<double> &x, SVec<double,dim> &ATx, int *ndType)
{
  ATx = 0.0;
  applyTransposeAndAdd(x,ATx,ndType);
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::applyTransposeAndAdd(Vec<double> &x, SVec<double,dim> &ATx, int *ndType)
{

  if(dim2 != 1) { fprintf(stderr, " *** Error: dim2 must be 1\n");  exit(-1); }
  int i, l;
  for (i=0; i<n; ++i) {
    for (int k=ia[i]; k<ia[i+1]; ++k)
      RectangularDenseMatrixOp<Scalar,dim,dim2,dim*dim2>::applyTransposeAndAddToVector(a.v, k, x.v, i, ATx.v, ja[k]);  // Gives ATx.v = transpose(a).v*x.v   
  }

  if (ndType) {
    for (i=0; i<n; ++i)
      if (ndType[i] != BC_INTERNAL)
        for (l=0; l<dim; ++l) ATx[i][l] = 0.0;
  }
}

//------------------------------------------------------------------------------
// a(i,j) in the original matrix becomes a(perm(i),perm(j))

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::permute(int *perm)
{

  rperm(perm);

  cperm(perm);

  for (int i=0; i<n; ++i) 
#ifdef OLD_STL
    sort(ja.v+ia[i], ja.v+ia[i+1]);
#else
    stable_sort(ja.v+ia[i], ja.v+ia[i+1]);
#endif
  
}

//------------------------------------------------------------------------------
// a(i,j) in the original matrix becomes a(perm(i),j)

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::rperm(int *perm)
{

  int *iao = new int[n+1];
  int *jao = new int[nnz];

  // determine pointers for output matrix. 

  int j;
  for (j=0; j<n; ++j) {
    int i = perm[j];
    iao[i+1] = ia[j+1] - ia[j];
  }

  // get pointers from lengths

  iao[0] = 0;
  for (j=0; j<n; ++j) iao[j+1] += iao[j];

  // copying 

  for (int ii=0; ii<n; ++ii) {

    // old row = ii  -- new row = perm(ii) -- ko = new pointer

    int ko = iao[ perm[ii] ];

    for (int k=ia[ii]; k<ia[ii+1]; ++k) { 
      jao[ko] = ja[k];
      ++ko;
    }

  }

  int i;
  for (i=0; i<nnz; ++i) ja[i] = jao[i];
  for (i=0; i<n+1; ++i) ia[i] = iao[i];
  
  if (iao) delete [] iao;
  if (jao) delete [] jao;

}

//------------------------------------------------------------------------------
// a(i,j) in the original matrix becomes a(i,perm(j))

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::cperm(int *perm) 
{

  for (int k=0; k<nnz; ++k) ja[k] = perm[ ja[k] ];

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
template<class Scalar2, int dim3>
inline
void RectangularSparseMat<Scalar,dim,dim2>::renumVector(SVec<Scalar2, dim3> &x)  {

  if (!nodeRenum) return;

  SVec<Scalar2, dim3> tmpVec(x);

  for (int i=0; i<n; ++i) {
    int inew = nodeRenum->renum[i];
    for (int j=0; j<dim3; ++j)
      x[inew][j] = tmpVec[i][j];
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
template<class Scalar2, int dim3>
inline
void RectangularSparseMat<Scalar,dim,dim2>::orderVector(SVec<Scalar2,dim3> &x)
{

  if (!nodeRenum) return;

  //*tmp = x;
  SVec<Scalar2,dim3> tmpVec(x);

  for (int i=0; i<n; ++i) {
    int iold = nodeRenum->order[i];
    for (int j=0; j<dim3; ++j) 
      x[iold][j] = tmpVec[i][j];
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
inline
void RectangularSparseMat<Scalar,dim,dim2>::convert2fortran()
{

  if (fortran == 1) return;

  fortran = 1;

  int i;
  for (i=0; i<nnz; ++i) ja[i] += fortran;
  for (i=0; i<n+1; ++i) ia[i] += fortran;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
inline
void RectangularSparseMat<Scalar,dim,dim2>::convert2cplusplus()
{

  if (fortran == 0) return;

  int i;
  for (i=0; i<nnz; ++i) ja[i] -= fortran;
  for (i=0; i<n+1; ++i) ia[i] -= fortran;

  fortran = 0;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::print(FILE *fp)
{

  if (fortran == 0) fprintf(fp, "*** Warning: C++ numbering\n");
  else fprintf(fp, "*** Warning: f77 numbering\n");

  int iold, icol;
  for (int i=0; i<n; ++i) {

    if (nodeRenum) iold = nodeRenum->order[i];
    else iold = i;

    fprintf(fp, "%d(%d):", i+fortran, iold+fortran);

    for (int j=ia[i]-fortran; j<ia[i+1]-fortran; ++j) {

      if (nodeRenum) icol = nodeRenum->order[ja[j]-fortran];
      else icol = ja[j]-fortran;

      fprintf(fp, " %d(%d)", ja[j], icol+fortran);

    }

    fprintf(fp, "\n");

  }
  
  exit(1);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
void RectangularSparseMat<Scalar,dim,dim2>::printRow(int i, int *locToGlobNodeMap, FILE *fp)
{

  int iglob, inew, icol;

  if (locToGlobNodeMap) iglob = locToGlobNodeMap[i];
  else iglob = -1;

  if (nodeRenum) inew = nodeRenum->renum[i];
  else inew = i;

  fprintf(fp, "line %d[%d](%d):\n\n", i, iglob, inew);

  for (int j=ia[inew]-fortran; j<ia[inew+1]-fortran; ++j) {

    if (nodeRenum) icol = nodeRenum->order[ja[j]-fortran];
    else icol = ja[j]-fortran;

    if (locToGlobNodeMap) iglob = locToGlobNodeMap[icol];

    fprintf(fp, "%d[%d](%d):\n", icol, iglob, ja[j]-fortran);

    for (int k=0; k<dim2; ++k) {

      for (int l=0; l<dim; ++l) {
	
        fprintf(fp, " %e", a[j][dim*k + l]);
	
      }

      fprintf(fp, "\n");

    }

  }

  fprintf(fp, "\n");

}

//------------------------------------------------------------------------------
