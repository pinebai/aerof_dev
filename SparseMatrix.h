#ifndef _SPARSE_MATRIX_H_
#define _SPARSE_MATRIX_H_

#include <cstdio>
#include <Vector.h>
#include <GenMatrix.h>
#include <Connectivity.h>

#include <complex>
typedef std::complex<double> bcomp;

class EdgeSet;

//------------------------------------------------------------------------------

template<class Scalar, int dim>
class SparseMat : public GenMat<Scalar,dim> {

  int n;                    // number of unknown blocks
  int nnz;                  // number of non zero blocks
  int fortran;              // fortran=1 if fortran style numbering

  Vec<int> ia;              // pointer to lines
  Vec<int> ja;              // pointer to columns
  Vec<int> ju;              // pointer to U factors in LU decomposition

  SVec<Scalar,dim*dim> a;   // data stored in CSR format

  compStruct *nodeRenum;    // unknowns renumbering structure
  //SVec<double,dim> *tmp;    // tmp vector used for renumbering
  //SVec<bcomp,dim> *tmpC;    // tmp vector used for renumbering

  Vec<int> ptr_ii;          // pointer to diagonal elmenents
  Vec<int> ptr_ij;          // pointer to off diagonal elements
  Vec<int> ptr_ji;          // pointer to off diagonal elements

  int *kc;                  // pointer to cols, transposed
  int *kk;                  // pointer to data, transposed
  int *kr;                  // pointer to rows, transposed
  int *ku;                  // pointer to L factors, transposed

  int *nodeType;

private:

  int find(int, int);

  void rperm(int *);
  void cperm(int *);
  template<class Scalar2>
  void renumVector(SVec<Scalar2,dim> &);
  template<class Scalar2>
  void orderVector(SVec<Scalar2,dim> &);

  void convert2fortran();
  void convert2cplusplus();

public:

  SparseMat(int _n, int _nnz, int *_ia, int *_ja, Scalar (*_a)[dim*dim],
	    compStruct *ndRenum, int *ndType) : 
    ia(_n+1, _ia), ja(_nnz, _ja), ju(0), a(_nnz, _a), ptr_ii(0), ptr_ij(0), ptr_ji(0) {
    n = _n;
    nnz = _nnz;
    fortran = 0;
    nodeRenum = ndRenum;
    nodeType = ndType;
    //if (nodeRenum) tmp = new SVec<double,dim>(n);
    //else tmp = 0;
  }
  ~SparseMat();

  SparseMat<Scalar,dim> &operator= (const Scalar x) { a = x; return *this; }
  SparseMat<Scalar,dim> &operator*= (const Scalar x) { a *= x; return *this; }

  int numBlocks() const { return n; }
  int numNonZeroBlocks() const { return nnz; }
  int *colind() { return ja.v; }
  int *rowptr() { return ia.v; }
  Scalar (*data())[dim*dim] { return a.v; }

  double norm() {return a.norm();}

  Scalar *getElem_ii(int i) { /*fprintf(stdout, "getElem_ii(%d) = %d\n", i, ptr_ii[i]);*/ return *(a.v + ptr_ii[i]); }
  Scalar *getElem_ij(int edgeNumber) { /*fprintf(stdout, "getElem_ij(%d) = %d\n", edgeNumber, ptr_ij[edgeNumber]);*/ return *(a.v + ptr_ij[edgeNumber]); }
  Scalar *getElem_ji(int edgeNumber) { /*fprintf(stdout, "getElem_ji(%d) = %d\n", edgeNumber, ptr_ji[edgeNumber]);*/ return *(a.v + ptr_ji[edgeNumber]); }

  void addContrib(int, int *, double *);

  void createPointers(EdgeSet &);

  void symbolicILU(const int);

  void numericILU(int *);

  void ILUTR();

  template<class Scalar2>
  void lusol(SVec<Scalar2,dim> &, SVec<Scalar2,dim> &);

  template<class Scalar2>
  void lusolTR(SVec<Scalar2,dim> &, SVec<Scalar2,dim> &);

  void apply(SVec<double,dim> &, SVec<double,dim> &, int * = 0);

  void permute(int *);

  template<class MatScal>
  void getData(GenMat<MatScal,dim> &);

  void zeroData() { a = 0.0; }

  void print(FILE * = stderr);
  void printRow(int, int * =0, FILE * = stderr);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <SparseMatrix.C>
#endif

#endif
