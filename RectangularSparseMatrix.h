#ifndef _RECTANGULAR_SPARSE_MATRIX_H_
#define _RECTANGULAR_SPARSE_MATRIX_H_

#include <cstdio>
#include <Vector.h>
#include <Vector3D.h>
#include <Connectivity.h>

#include <complex>
typedef std::complex<double> bcomp;

class EdgeSet;

//------------------------------------------------------------------------------

template<class Scalar, int dim, int dim2>
class RectangularSparseMat {

  int n;                    // number of unknown blocks
  int nnz;                  // number of non zero blocks
  int fortran;              // fortran=1 if fortran style numbering

  Vec<int> ia;              // pointer to lines
  Vec<int> ja;              // pointer to columns

  SVec<Scalar,dim*dim2> a;   // data stored in CSR format

  compStruct *nodeRenum;    // unknowns renumbering structure
  //SVec<double,dim> *tmp;    // tmp vector used for renumbering
  //SVec<bcomp,dim> *tmpC;    // tmp vector used for renumbering

  Vec<int> ptr_ii;          // pointer to diagonal elmenents
  Vec<int> ptr_ij;          // pointer to off diagonal elements
  Vec<int> ptr_ji;          // pointer to off diagonal elements

  int *nodeType;

private:

  int find(int, int);

  void rperm(int *);
  void cperm(int *);
  template<class Scalar2, int dim3>
  void renumVector(SVec<Scalar2,dim3> &);
  template<class Scalar2, int dim3>
  void orderVector(SVec<Scalar2,dim3> &);

  void convert2fortran();
  void convert2cplusplus();

public:

  RectangularSparseMat(int _n, int _nnz, int *_ia, int *_ja, Scalar (*_a)[dim*dim2],
	    compStruct *ndRenum, int *ndType) : 
    ia(_n+1, _ia), ja(_nnz, _ja), a(_nnz, _a), ptr_ii(0), ptr_ij(0), ptr_ji(0) {
    n = _n;
    nnz = _nnz;
    fortran = 0;
    nodeRenum = ndRenum;
    nodeType = ndType;
  }
  ~RectangularSparseMat();

  RectangularSparseMat<Scalar,dim,dim2> &operator= (const Scalar x) { a = x; return *this; }
  RectangularSparseMat<Scalar,dim,dim2> &operator*= (const Scalar x) { a *= x; return *this; }

  int numBlocks() const { return n; }
  int numNonZeroBlocks() const { return nnz; }
  int *colind() { return ja.v; }
  int *rowptr() { return ia.v; }
  Scalar (*data())[dim*dim2] { return a.v; }

  double norm() {return a.norm();}

  Scalar *getElem_ii(int i) { /*fprintf(stdout, "getElem_ii(%d) = %d\n", i, ptr_ii[i]);*/ return *(a.v + ptr_ii[i]); }
  Scalar *getElem_ij(int edgeNumber) { /*fprintf(stdout, "getElem_ij(%d) = %d\n", edgeNumber, ptr_ij[edgeNumber]);*/ return *(a.v + ptr_ij[edgeNumber]); }
  Scalar *getElem_ji(int edgeNumber) { /*fprintf(stdout, "getElem_ji(%d) = %d\n", edgeNumber, ptr_ji[edgeNumber]);*/ return *(a.v + ptr_ji[edgeNumber]); }

  void addContrib(int, int *, double *);
  void addContrib(int, int, double *);
  void addContrib(int, int, double);

  void createPointers(EdgeSet &);

  void apply(SVec<double,dim> &, SVec<double,dim2> &, int * = 0);
  void apply(SVec<double,dim> &, Vec<Vec3D> &, int * = 0);
  void apply(Vec<Vec3D> &, SVec<double,dim2> &, int * = 0);
  void apply(Vec<double> &, SVec<double,dim2> &, int * = 0);
  void apply(SVec<double,dim> &, Vec<double> &, int * = 0);

  void applyAndAdd(SVec<double,dim> &, SVec<double,dim2> &, int * = 0);
  void applyAndAdd(SVec<double,dim> &, Vec<Vec3D> &, int * = 0);
  void applyAndAdd(Vec<Vec3D> &, SVec<double,dim2> &, int * = 0);
  void applyAndAdd(Vec<double> &, SVec<double,dim2> &, int * = 0);
  void applyAndAdd(SVec<double,dim> &, Vec<double> &, int * = 0);

  void applyTranspose(SVec<double,dim2> &, SVec<double,dim> &, int * = 0);
  void applyTranspose(Vec<Vec3D> &, SVec<double,dim> &, int * = 0);
  void applyTranspose(SVec<double,dim2> &, Vec<Vec3D> &, int * = 0);
  void applyTranspose(Vec<double> &, SVec<double,dim> &, int * = 0);
  void applyTranspose(SVec<double,dim2> &, Vec<double> &, int * = 0);

  void applyTransposeAndAdd(SVec<double,dim2> &, SVec<double,dim> &, int * = 0);
  void applyTransposeAndAdd(Vec<Vec3D> &, SVec<double,dim> &, int * = 0);
  void applyTransposeAndAdd(SVec<double,dim2> &, Vec<Vec3D> &, int * = 0);
  void applyTransposeAndAdd(Vec<double> &, SVec<double,dim> &, int * = 0);
  void applyTransposeAndAdd(SVec<double,dim2> &, Vec<double> &, int * = 0);




  void permute(int *);

  void zeroData() { a = 0.0; }

  void print(FILE * = stderr);
  void printRow(int, int * =0, FILE * = stderr);
  void printFirstElementIn_a();
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <RectangularSparseMatrix.C>
#endif

#endif
