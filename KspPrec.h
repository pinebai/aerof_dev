#ifndef _KSP_PREC_H_
#define _KSP_PREC_H_

#include <IoData.h>
#include <DistVector.h>
#include <DistMatrix.h>
#include <DiagMatrix.h>
#include <SparseMatrix.h>
#include <DistEmbeddedVector.h>

#include <complex>
typedef std::complex<double> bcomp;
class MemoryPool;
class BCApplier;

//------------------------------------------------------------------------------

template<int dim, class Scalar2 = double>
class KspPrec {

public:

  KspPrec() {}
  virtual ~KspPrec() {}
  
  virtual void exportMemory(MemoryPool *mp) {}

  virtual void setup() = 0;
  virtual void setupTR() { std::cout << "setupTR not implemented for this Preconditioner" << endl; exit(-1);}

  // all KspPrec derived class must have the following functions
  virtual void apply(DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &)  { std::cout << "  ERROR: Using default apply function in KspPrec" << endl;  exit(-1);}
  virtual void apply(DistVec<Scalar2> &, DistVec<Scalar2> &)  { std::cout << "  ERROR: Using default apply function in KspPrec" << endl;  exit(-1);}
  virtual void applyT(DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &)  { std::cout << "  ERROR: Using default applyT function in KspPrec" << endl;  exit(-1);}
  virtual void applyT(DistEmbeddedVec<Scalar2,dim>& x, DistEmbeddedVec<Scalar2,dim>& Px) { std::cout<< " ERROR: Using default applyT function in KspPrec" << endl;  exit(-1);}
  //virtual void applyTranspose(DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &)  { std::cout << "  ERROR: Using default applyTranspose function in KspPrec" << endl;  exit(-1);}
  //virtual void applyTranspose(DistEmbeddedVec<Scalar2,dim>& x, DistEmbeddedVec<Scalar2,dim>& Px) { std::cout<< " ERROR: Using default applyTranspose function in KspPrec" << endl;  exit(-1);}

  void apply(DistEmbeddedVec<Scalar2,dim>& x, DistEmbeddedVec<Scalar2,dim>& Px) { 
    apply(x.real(), Px.real());
    Px.ghost() = x.ghost();
    if (x.hasHHBoundaryTerm())
      Px.hh() = x.hh();
  }
};

//------------------------------------------------------------------------------

template<int dim, class Scalar2 = double>
class IdentityPrec : public KspPrec<dim, Scalar2> {

  BCApplier* BCs;
public:

  IdentityPrec(BCApplier *bcs = 0) : BCs(bcs) {}
  ~IdentityPrec() {}
  
  void setup() {}
  void setupTR() {}

  void apply(DistSVec<Scalar2,dim> &x, DistSVec<Scalar2,dim> &Ix); 
  void applyT(DistSVec<Scalar2,dim> &x, DistSVec<Scalar2,dim> &Ix) { Ix = x; }
  void applyTranspose(DistSVec<Scalar2,dim> &x, DistSVec<Scalar2,dim> &Ix) { Ix = x; }

};

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2 = double>
class JacobiPrec : public KspPrec<dim, Scalar2>, public DistMat<Scalar,dim> {

#ifdef _OPENMP
  int numLocSub; //BUG omp
#endif

  DiagMat<Scalar,dim> **A;
  BCApplier* BCs; //HB
public:

  JacobiPrec(typename DiagMat<Scalar,dim>::Type, Domain *, int ** = 0, BCApplier* =0);
  ~JacobiPrec();

  DistMat<Scalar,dim> &operator= (const Scalar);

  GenMat<Scalar,dim> &operator() (int i) { return *A[i]; }
  
  void setup();

  void apply(DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);
  template<class MatScal>
  void getData(DistMat<MatScal,dim> &);

};

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2 = double>
class IluPrec : public KspPrec<dim, Scalar2>, public DistMat<Scalar,dim> {

#ifdef _OPENMP
  int numLocSub; //BUG omp
#endif

  PcData::Type type;

  SparseMat<Scalar,dim> **A;

  DistVec<int> jw;

  //DistSVec<double,dim> tmp;
  int fill;

  Domain* domain;

public:

  IluPrec(PcData &, Domain *, int ** = 0);
  ~IluPrec();

  DistMat<Scalar,dim> &operator= (const Scalar);
  
  GenMat<Scalar,dim> &operator() (int i) { return *A[i]; }

  void exportMemory(MemoryPool *);

  void setup();
  void setupTR();

  void apply(DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);
  void applyT(DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);
  void applyTranspose(DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);

  template<class MatScal>
  void getData(DistMat<MatScal,dim> &);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <KspPrec.C>
#endif

#endif
