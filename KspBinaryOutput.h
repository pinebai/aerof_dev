#ifndef _KSP_BINARY_OUTPUT_H_
#define _KSP_BINARY_OUTPUT_H_

#include <cstdio>
#include <Vector.h>
#include <DenseMatrix.h>
#include <VectorSet.h>
#include <DistVector.h>
#include <DistEmbeddedVector.h>
#include <cstring>
template <class Scalar, int dim> class DistSVec;

template <class VecType>
class KspBinaryOutput {

//This class is used to output the krylov vectors from the linear solver (GMRES only) 

  protected:
  Domain* domain;
  Communicator* com; 
  IoData* ioData;

  int krylovFreqTime;
  int krylovFreqNewton;
  double krylovEnergy;
  int* timeIt;
  int* newtonIt;
  char* fileName;

  VecType* U;

  public:

  KspBinaryOutput(Communicator *, IoData *, Domain *);
  ~KspBinaryOutput();

  void setCurrentState(VecType&);

  template <int dim>
  void writeKrylovVectors(VecSet<DistSVec<double, dim> >&, Vec<double>, int);

  template <int dim>
  void writeKrylovVectors(VecSet<DistSVec<bcomp, dim> >&, Vec<bcomp>, int);

  template <class Scalar, int dim>
  void writeKrylovVectors(VecSet<DistEmbeddedVec<Scalar, dim> >&, Vec<Scalar>, int);
};

#include "KspBinaryOutput.C"
#endif
