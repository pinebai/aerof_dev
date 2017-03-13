#include <cstdlib>

#ifdef OLD_STL
#include <defalloc.h>
#include <algo.h>
#else
#include <algorithm>
using std::stable_sort;
#endif

#include <SparseMatrix.h>

#include <BcDef.h>
#include <Edge.h>
#include <Vector.h>
#include <ResizeArray.h>
#include <DenseMatrixOps.h>

//------------------------------------------------------------------------------

template<class Scalar, int dim>
SparseMat<Scalar,dim>::~SparseMat()
{ 

  if (nodeRenum) delete nodeRenum;
  //if (tmp) delete tmp;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
int SparseMat<Scalar,dim>::find(int i, int j) 
{

  if (i >= n || i < 0) {
    fprintf(stderr, "*** Error: out of bound tet: %d vs %d\n", i, n);
    exit(-1);
  }

  for (int ii = ia[i]; ii < ia[i+1]; ++ii)
    if (ja[ii] == j) return ii;

  fprintf(stderr, "*** Error: connection not found in sparseMat: %d %d\n", i, j);

  return -1;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SparseMat<Scalar,dim>::createPointers(EdgeSet &edges)
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

template<class Scalar, int dim>
void SparseMat<Scalar,dim>::addContrib(int nNd, int *ndList, double *C)
{

  int i, j, k, l;
  int cSize = nNd*dim;

  for (i = 0; i < nNd; ++i) {
    int iNode = (nodeRenum) ? nodeRenum->renum[ ndList[i] ] : ndList[i];
    if (nodeType && nodeType[ ndList[i] ] != BC_INTERNAL) {
      int indx = find(iNode, iNode);
      for (k = 0; k < dim; ++k)
	for (l = 0; l < dim; ++l)
	  a[indx][dim*k+l] = (k == l) ? 1.0 : 0.0;
      continue;
    }

    for (j = 0; j < nNd; ++j) {
      if (nodeType && nodeType[ ndList[j] ] != BC_INTERNAL) continue;
      int jNode = (nodeRenum) ? nodeRenum->renum[ ndList[j] ] : ndList[j];
      int indx = find(iNode, jNode);
      for (k = 0; k < dim; ++k)
	for (l = 0; l < dim; ++l) {
	  a[indx][dim*k+l] += C[(i*dim+k)*cSize + j*dim + l];
	}
    }
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<class MatScal>
void SparseMat<Scalar,dim>::getData(GenMat<MatScal,dim> &B)
{

  int k;

  Scalar *a1;
  MatScal *a2;

  for (int l=0; l<ptr_ij.size(); ++l) {
    a1 = getElem_ij(l);
    a2 = B.getElem_ij(l);
    for (k=0; k<dim*dim; ++k) 
      a1[k] = a2[k];

    a1 = getElem_ji(l);
    a2 = B.getElem_ji(l);
    for (k=0; k<dim*dim; ++k) 
      a1[k] = a2[k];
  }

  for (int i=0; i<ptr_ii.size(); ++i) {
    a1 = getElem_ii(i);
    a2 = B.getElem_ii(i);
    for (k=0; k<dim*dim; ++k) 
      a1[k] = a2[k];
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SparseMat<Scalar,dim>::symbolicILU(const int levfill)
{

  int i, j, k;

  ju.resize(n);

  if (levfill == 0) {

    int ndiag = 0;

    for (i=0; i<n; ++i) {

#ifdef OLD_STL
      sort(ja.v+ia[i], ja.v+ia[i+1]);
#else
      stable_sort(ja.v+ia[i], ja.v+ia[i+1]);
#endif

      for (k=ia[i]; k<ia[i+1]; ++k)
	if (ja[k] == i) { ju[i] = k; ++ndiag; break; }

    }

    if (ndiag == n) return;
    else fprintf(stderr, "*** Warning: adding %d zero diagonal terms\n", n-ndiag);

  }

  Vec<int> ial(n+1);
  Vec<int> iau(n+1);
  Vec<int> lnklst(n);
  Vec<int> curlev(n);
  Vec<int> currow(n);

  ResizeArray<int> jal(0, nnz);
  ResizeArray<int> jau(0, nnz);
  ResizeArray<int> levels(0, nnz);

  int nnzl = 0;
  int nnzu = 0;

  ial[0] = 0;
  iau[0] = 0;

  for (i=0; i<n; ++i) {

    // copy column indices of row into workspace and sort them

    int len = ia[i+1] - ia[i];

    k = 0;
    for (j=ia[i]; j<ia[i+1]; ++j) currow[k++] = ja[j];

#ifdef OLD_STL
    sort(currow.v, currow.v+len);
#else
    stable_sort(currow.v, currow.v+len);
#endif

    // construct implied linked list for row

    for (j=0; j<len-1; ++j) {
      lnklst[ currow[j] ] = currow[j+1];
      curlev[ currow[j] ] = 0;
    }

    lnklst[ currow[len-1] ] = n;
    curlev[ currow[len-1] ] = 0;

    // merge with rows in U

    int first = currow[0];
    int next = first;
    while (next < i) {

      int oldlst = next;
      int nxtlst = lnklst[next];
      int row = next;

      // scan row

      for (int ii=iau[row]+1; ii<iau[row+1]; /*nop*/) {

	int newlev;
	if (jau[ii] < nxtlst) {

	  // new fill-in

	  newlev = curlev[row] + levels[ii] + 1;
	  if (newlev <= levfill) {
	    lnklst[oldlst] = jau[ii];
	    lnklst[ jau[ii] ] = nxtlst;
	    oldlst = jau[ii];
	    curlev[ jau[ii] ] = newlev;
	  }
	  ii++;

	} 
	else if (jau[ii] == nxtlst) {

	  // update fill-in

	  oldlst = nxtlst;
	  nxtlst = lnklst[oldlst];
	  newlev = curlev[row] + levels[ii] + 1;
	  curlev[ jau[ii] ] = curlev[jau[ii]] < newlev ? curlev[jau[ii]] : newlev;
	  ii++;

	} 
	else {
	  
	  oldlst = nxtlst;
	  nxtlst = lnklst[oldlst];

	}

      }

      next = lnklst[next];

    }
        
    // gather the pattern into L and U

    next = first;
    while (next < i) {

      jal[nnzl++] = next;
      next = lnklst[next];

    }

    ial[i+1] = nnzl;

    if (next != i) {
      fprintf(stderr, "*** Warning: adding diagonal term on line %d\n", i);
      levels[nnzu] = 2*n; // infinity
      jau[nnzu++] = i;
    }

    while (next < n) {

      levels[nnzu] = curlev[next];
      jau[nnzu++] = next;
      next = lnklst[next];

    }

    iau[i+1] = nnzu;

  }

  // create new ia and ja vectors from ial,jal and iau,jau

  nnz = nnzl + nnzu;
  ja.resize(nnz);
  a.resize(nnz);

  int nnzlu = 0;
  ia[0] = 0;
  
  for (i=0; i<n; ++i) {

    for (k=ial[i]; k<ial[i+1]; ++k) ja[nnzlu++] = jal[k];

    ju[i] = nnzlu;

    for (k=iau[i]; k<iau[i+1]; ++k) ja[nnzlu++] = jau[k];
    
    ia[i+1] = nnzlu;

  }

  if (nnzlu != nnz) {
    fprintf(stderr, "*** Error: nnzlu incorrect\n");
    exit(1);
  }
  
}

//------------------------------------------------------------------------------
// a stores both L and U since L is a unit lower triangular matrix
// ju points to the diagonal elements
// the diagonal elements are already inverted

template<class Scalar, int dim>
void SparseMat<Scalar,dim>::numericILU(int *marker)
{
  int i, j, k, kk, l, m, ierr;

  Scalar mult[1][dim*dim], res[1][dim*dim];

  for (i=0; i<n; ++i) 
    marker[i] = -1;

  for (i=0; i<n; ++i) {

    for (k=ia[i]; k<ia[i+1]; ++k) 
      marker[ ja[k] ] = k;

    for (k=ia[i]; k<ju[i]; ++k) {
      
      l = ja[k];

      //Scalar mult = a[k] * a[ ju[l] ];
      //a[k] = mult;
      DenseMatrixOp<Scalar,dim,dim*dim>::applyToDenseMatrix(a.v, k, a.v, ju[l], mult, 0);
      for (kk=0; kk<dim*dim; ++kk) 
	a[k][kk] = mult[0][kk];

      for (j=ju[l]+1; j<ia[l+1]; ++j) {
	m = marker[ ja[j] ];
	if (m != -1) {
	  //a[m] -= mult * a[j];
	  DenseMatrixOp<Scalar,dim,dim*dim>::applyToDenseMatrix(mult, 0, a.v, j, res, 0);
	  for (kk=0; kk<dim*dim; ++kk) a[m][kk] -= res[0][kk];
	}
      }

    }

    for (k=ia[i]; k<ia[i+1]; ++k) 
      marker[ ja[k] ] = -1;

    ierr = invertDenseMatrix<Scalar,dim>(a[ ju[i] ]);

    if (ierr) {
      fprintf(stderr, "*** Error: zero pivot at row %d in numeric ILU\n", i);
      exit(1);
    }

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<class Scalar2>
void SparseMat<Scalar,dim>::lusol(SVec<Scalar2,dim> &y, SVec<Scalar2,dim> &x)
{

  int i, k, k1, k2, l;

  Scalar2 res[1][dim];

  renumVector(y);

  // forward solve

  for (i=0; i<n; ++i) {

    for (l=0; l<dim; ++l) x[i][l] = y[i][l];

    k1 = ia[i];
    k2 = ju[i];

#pragma ivdep
    for (k=k1; k<k2; ++k)
      DenseMatrixOp<Scalar,dim,dim*dim>::applyAndSubToVector(a.v, k, x.v, ja[k], x.v, i);

  }

  // backward solve

  for (i=n-1; i>=0; --i) {

    k1 = ju[i]+1;
    k2 = ia[i+1];

#pragma ivdep
    for (k=k1; k<k2; ++k)
      DenseMatrixOp<Scalar,dim,dim*dim>::applyAndSubToVector(a.v, k, x.v, ja[k], x.v, i);

    DenseMatrixOp<Scalar,dim,dim*dim>::applyToVector(a.v, ju[i], x.v, i, res, 0);

    for (l=0; l<dim; ++l) x[i][l] = res[0][l];

  }

  orderVector(y);
  orderVector(x);

  if (nodeType) {
    for (i=0; i<n; ++i)
      if (nodeType[i] != BC_INTERNAL) 
	for (l=0; l<dim; ++l) x[i][l] = 0.0;
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
void SparseMat<Scalar,dim>::ILUTR()
{

  int i, k, k1, k2, l;


  int *tc;
 
  kc = new int[n+1];
  ku = new int[n];
  tc = new int[n];
  for (i=0; i<n; ++i) tc[i] = 0;
  kk = new int[nnz];
  kr = new int[nnz];
  for (i=0; i<nnz; ++i) kr[i] = 0;
 
  //Develop a transformation from CSR format to CSC format
  //first build column pointers
  for (i=0; i<n; ++i) {
    for (k=ia[i]; k<ia[i+1]; ++k)
      tc[ ja[k] ] += 1;
  }
  kc[0] = 0;
 
  for (i=1; i<=n; ++i)
      kc[ i ] = kc[i-1] + tc[i-1];

  //next build row pointers
  k = 0;
  for (i=0; i<n; ++i) {
      for (l=0; l<n; ++l) {
        for (k1=ia[l]; k1<ia[l+1]; ++k1) {
          if (ja[k1] == i) {
            kr[k] = l;
            kk[k] = k1;
            k +=1;
          }
        }
      }
  }

  //Finally build L pointers
  ku[0] = 0;
  for (i=1; i<n; ++i) {
    k = 0;
    for (k1=kc[i]; k1<kc[i+1]; ++k1) {
      if (kr[k1] < i) k += 1;
      if(kr[k1] == i) ku[i] = k + kc[i];
    }
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<class Scalar2>
inline
void SparseMat<Scalar,dim>::lusolTR(SVec<Scalar2,dim> &y, SVec<Scalar2,dim> &x)
{
  int i, k, k1, k2, l;

  Scalar2 tmp[1][dim];

  renumVector(y);
 
  // forward solve

  for (i=0; i<n; ++i) {

    for (l=0; l<dim; ++l) x[i][l] = y[i][l];

    k1 = kc[i];
    k2 = ku[i];

#pragma ivdep
    for (k=k1; k<k2; ++k){
      DenseMatrixOp<Scalar,dim,dim*dim>::applyTransAndSubToVector(a.v, kk[k], x.v, kr[k], x.v, i);
      //subDenseMatrixTransTimesVector(a.v, kk[k], x.v, ja[k], x.v, i);
    }

    //denseMatrixTransTimesVector(a.v, ju[i], x.v, i, tmp, 0);
    DenseMatrixOp<Scalar,dim,dim*dim>::applyTransToVector(a.v, ju[i], x.v, i, tmp, 0);

    for (l=0; l<dim; ++l) x[i][l] = tmp[0][l];

  }

  // backward solve

  for (i=n-1; i>=0; --i) {

    k1 = ku[i]+1;
    k2 = kc[i+1];
#pragma ivdep
    for (k=k1; k<k2; ++k)
    	for (k=k1; k<k2; ++k) {
    	  DenseMatrixOp<Scalar,dim,dim*dim>::applyTransAndSubToVector(a.v, kk[k], x.v, kr[k], x.v, i);
    	  //subDenseMatrixTransTimesVector(a.v, kk[k], x.v, kr[k], x.v, i);
    	}
      //subDenseMatrixTransTimesVector(a.v, kk[k], x.v, kr[k], x.v, i);

  }

  orderVector(y);
  orderVector(x);

  if (nodeType) {
    for (i=0; i<n; ++i)
      if (nodeType[i] != 0)
        for (l=0; l<dim; ++l) x[i][l] = 0.0;
  }


}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SparseMat<Scalar,dim>::apply(SVec<double,dim> &x, SVec<double,dim> &Ax, int *ndType)
{

  int i, l;
  for (i=0; i<n; ++i) {

    for (l=0; l<dim; ++l) Ax[i][l] = 0.0;

    for (int k=ia[i]; k<ia[i+1]; ++k)
      DenseMatrixOp<Scalar,dim,dim*dim>::applyAndAddToVector(a.v, k, x.v, ja[k], Ax.v, i);  // Gives Ax.v = a.v*x.v   

  }

  if (ndType) {
    for (i=0; i<n; ++i)
      if (ndType[i] != BC_INTERNAL)
	for (l=0; l<dim; ++l) Ax[i][l] = 0.0;
  }

}

//------------------------------------------------------------------------------
// a(i,j) in the original matrix becomes a(perm(i),perm(j))

template<class Scalar, int dim>
void SparseMat<Scalar,dim>::permute(int *perm)
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

template<class Scalar, int dim>
void SparseMat<Scalar,dim>::rperm(int *perm)
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

template<class Scalar, int dim>
void SparseMat<Scalar,dim>::cperm(int *perm) 
{

  for (int k=0; k<nnz; ++k) ja[k] = perm[ ja[k] ];

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<class Scalar2>
inline
void SparseMat<Scalar,dim>::renumVector(SVec<Scalar2, dim> &x)  {

  if (!nodeRenum) return;

  SVec<Scalar2, dim> tmpVec(x);

  for (int i=0; i<n; ++i) {
    int inew = nodeRenum->renum[i];
    for (int j=0; j<dim; ++j)
      x[inew][j] = tmpVec[i][j];
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<class Scalar2>
inline
void SparseMat<Scalar,dim>::orderVector(SVec<Scalar2,dim> &x)
{

  if (!nodeRenum) return;

  //*tmp = x;
  SVec<Scalar2,dim> tmpVec(x);

  for (int i=0; i<n; ++i) {
    int iold = nodeRenum->order[i];
    for (int j=0; j<dim; ++j) 
      x[iold][j] = tmpVec[i][j];
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
void SparseMat<Scalar,dim>::convert2fortran()
{

  if (fortran == 1) return;

  fortran = 1;

  int i;
  for (i=0; i<nnz; ++i) ja[i] += fortran;
  for (i=0; i<n+1; ++i) ia[i] += fortran;
  for (i=0; i<n; ++i) ju[i] += fortran;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
void SparseMat<Scalar,dim>::convert2cplusplus()
{

  if (fortran == 0) return;

  int i;
  for (i=0; i<nnz; ++i) ja[i] -= fortran;
  for (i=0; i<n+1; ++i) ia[i] -= fortran;
  for (i=0; i<n; ++i) ju[i] -= fortran;

  fortran = 0;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SparseMat<Scalar,dim>::print(FILE *fp)
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

template<class Scalar, int dim>
void SparseMat<Scalar,dim>::printRow(int i, int *locToGlobNodeMap, FILE *fp)
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

    for (int k=0; k<dim; ++k) {

      for (int l=0; l<dim; ++l) {
	
	fprintf(fp, " %e", a[j][dim*k + l]);
	
      }

      fprintf(fp, "\n");

    }

  }

  fprintf(fp, "\n");

}

//------------------------------------------------------------------------------
