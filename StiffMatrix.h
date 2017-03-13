#ifndef _STIFF_MAT_H_
#define _STIFF_MAT_H_

#include <SparseMatrix.h>
#include <DistVector.h>
#include <DistMatrix.h>

class MemoryPool;
class BCApplier;

//------------------------------------------------------------------------------

template<class Scalar, int dim>
class StiffMat : public DistMat<Scalar,dim> {

#ifdef _OPENMP
  int numLocSub; //BUG omp
#endif

  int **ndType;
//  BCApplier* BCs; //HB

  SparseMat<Scalar,dim> **A;

public:
  BCApplier* BCs; //PJSA
  
  StiffMat(Domain *, int **, MemoryPool *, BCApplier* bcs = NULL);
  ~StiffMat();

  DistMat<Scalar,dim> &operator= (const Scalar);

  GenMat<Scalar,dim> &operator() (int i) { return *A[i]; }

  void apply(DistSVec<double,dim> &, DistSVec<double,dim> &);

  void applyTranspose(DistSVec<double,dim> &, DistSVec<double,dim> &)
  {std::cout<<__FILE__<<":"<<__LINE__<<" empty header declaration called"; exit(-1);}

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <StiffMatrix.C>
#endif

#endif
