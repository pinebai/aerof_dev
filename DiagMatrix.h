#ifndef _DIAG_MATRIX_H_
#define _DIAG_MATRIX_H_

#include <Vector.h>
#include <GenMatrix.h>

template<class Scalar, int dim> class SparseMat;

//------------------------------------------------------------------------------

template<class Scalar, int dim>
class DiagMat : public GenMat<Scalar,dim> {

public:

  enum Type {DENSE = 0, DIAGONAL = 1} type;

private:

  int n;

  int *nodeType;

  SVec<Scalar,dim*dim> a;

public:

  DiagMat(Type, int, int *);
  ~DiagMat() {}

  DiagMat<Scalar,dim> &operator= (const Scalar x) { a = x; return *this; }
  DiagMat<Scalar,dim> &operator*= (const Scalar x) { a *= x; return *this; }

  int numNonZeroBlocks() const { return n; }
  Scalar (*data())[dim*dim] { return a.v; }

  Scalar *getElem_ii(int i) { return *(a.v + i); }
  Scalar *getElem_ij(int l) { return 0; }
  Scalar *getElem_ji(int l) { return 0; }

  void addContrib(int, int *, double *);

  void invert();

  void apply(SVec<double,dim> &, SVec<double,dim> &);

  double norm() {return a.norm();}

  template<class MatScal>
  void getData(GenMat<MatScal,dim> &);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <DiagMatrix.C>
#endif

#endif
