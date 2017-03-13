#include "DistLeastSquareSolver.h"

#include <Communicator.h>

#include "LinkF77.h"

#include <stdexcept>
#include <algorithm>
#include <cstdlib>
#include <cassert>

#ifdef DO_SCALAPACK

#include <mpi.h>

extern "C" {
  // Context & cpu topology management
  int Csys2blacs_handle(MPI_Comm comm);
  void Cfree_blacs_system_handle(int handle);
  void Cblacs_gridinit(int * ictxt, char * order, int nprow, int npcol);
  void Cblacs_gridinfo(int ictxt, int * nprow, int * npcol, int * myrow, int * mycol);
  void Cblacs_gridexit(int ictxt);

  // Block-cycling
  int F77NAME(numroc)(const int * n, const int * nb, const int * iproc, const int * isrcproc, const int * nprocs);
  void F77NAME(descinit)(int * desc, const int * m, const int * n, const int * mb, const int * nb,
                         const int * irsrc, const int * icsrc, const int * ictxt, const int * lld, int * info);

  // Index mapping 
  int F77NAME(indxl2g)(const int * indxloc, const int * nb, const int * iproc, const int * isrcproc, const int * nprocs);
  int F77NAME(indxg2l)(const int * indxglob, const int * nb, const int * iproc_dummy, const int * isrcproc_dummy, const int * nprocs);
  int F77NAME(indxg2p)(const int * indxglob, const int * nb, const int * iproc_dummy, const int * isrcproc, const int * nprocs);

  // Computations
  void F77NAME(pdgels)(const char * trans, const int * m, const int * n, const int * nrhs,
                       double * a, const int * ia, const int * ja, const int desca[9],
                       double * b, const int * ib, const int * jb, const int descb[9],
                       double * work, const int * lwork, int * info);
}

#endif /* DO_SCALAPACK */

const int DistLeastSquareSolver::DEFAULT_BLOCK_SIZE = 32;

const int DistLeastSquareSolver::INT_ZERO = 0;
const int DistLeastSquareSolver::INT_ONE = 1;
const int DistLeastSquareSolver::INT_MINUS_ONE = -1;

DistLeastSquareSolver::DistLeastSquareSolver(Communicator * comm, int rowCpus, int colCpus) : 
  communicator_(comm),
  bufferResizePolicy_(TIGHT)
{
#ifdef DO_SCALAPACK
  if (rowCpus * colCpus > communicator_->size()) {
    throw std::invalid_argument("Not enough cpus"); 
  }
  
  if (rowCpus * colCpus < communicator_->size()) {
    throw std::invalid_argument("Too many cpus"); 
  }

  blacsHandle_ = Csys2blacs_handle(communicator_->comm);
  context_ = blacsHandle_;
  char order[] = "R";
  Cblacs_gridinit(&context_, order, rowCpus, colCpus);
  Cblacs_gridinfo(context_, &rowCpus_, &colCpus_, &localCpuRow_, &localCpuCol_);
 
  assert(rowCpus_ == rowCpus);
  assert(colCpus_ == colCpus);

  // Empty problem
  equationCount_ = 0;
  unknownCount_ = 0;
  largestDimension_ = 0;
  rhsCount_ = 0;

  blockSizeIs(DEFAULT_BLOCK_SIZE);

#else /* DO_SCALAPACK */
  throw std::runtime_error("ScaLAPACK not available");
#endif /* DO_SCALAPACK */
}

DistLeastSquareSolver::~DistLeastSquareSolver() {
#ifdef DO_SCALAPACK
  Cblacs_gridexit(context_);
  Cfree_blacs_system_handle(blacsHandle_);
#endif /* DO_SCALAPACK */
}

void
DistLeastSquareSolver::blockSizeIs(int size) {
  rowBlockSize_ = size;
  colBlockSize_ = size;
  rhsRowBlockSize_ = size;
  rhsColBlockSize_ = size;

  reset(); 
}

void
DistLeastSquareSolver::blockSizeIs(int row, int col, int rhsRow, int rhsRank) {
  rowBlockSize_ = row;
  colBlockSize_ = col;
  rhsRowBlockSize_ = rhsRow;
  rhsColBlockSize_ = rhsRank;

  reset(); 
}

void
DistLeastSquareSolver::bufferResizePolicyIs(BufferResizePolicy newPolicy) {
  if (newPolicy == bufferResizePolicy()) {
    return;
  }

  bufferResizePolicy_ = newPolicy;

  if (bufferResizePolicy() == TIGHT) {
    reset();
  }
}

void
DistLeastSquareSolver::problemSizeIs(int eqnCount, int unknownCount, int rhsCount) {
  equationCount_    = eqnCount;
  unknownCount_     = unknownCount;
  largestDimension_ = std::max(eqnCount, unknownCount);
  rhsCount_         = rhsCount;

  reset();
}

void
DistLeastSquareSolver::reset() {
#ifdef DO_SCALAPACK
  localRows_ = F77NAME(numroc)(&equationCount_, &rowBlockSize_, &localCpuRow_, &INT_ZERO, &rowCpus_);
  localCols_ = F77NAME(numroc)(&unknownCount_,  &colBlockSize_, &localCpuCol_, &INT_ZERO, &colCpus_);

  localRhsRows_      = F77NAME(numroc)(&equationCount_, &rhsRowBlockSize_, &localCpuRow_, &INT_ZERO, &rowCpus_);
  localSolutionRows_ = F77NAME(numroc)(&unknownCount_,  &rhsRowBlockSize_, &localCpuRow_, &INT_ZERO, &rowCpus_);
  localRhsCount_     = F77NAME(numroc)(&rhsCount_,      &rhsColBlockSize_, &localCpuCol_, &INT_ZERO, &colCpus_);

  const int reqMatrixBufferSize = localCols_ * localRows_;
  const int reqRhsColBufferSize = std::max(localRhsRows_, localSolutionRows_);
  const int reqRhsBufferSize    = localRhsCount_ * reqRhsColBufferSize;

  if (bufferResizePolicy() == TIGHT || reqMatrixBufferSize > matrixBuffer_.size()) {
    matrixBuffer_.sizeIs(reqMatrixBufferSize);
  }

  if (bufferResizePolicy() == TIGHT || reqRhsBufferSize > rhsBuffer_.size()) {
    rhsBuffer_.sizeIs(reqRhsBufferSize);
  }
    
  localMatrixLeadDim_ = std::max(localRows_, 1);
  localRhsLeadDim_    = std::max(reqRhsColBufferSize, 1);

  int info;
  F77NAME(descinit)(matrixDesc_, &equationCount_, &unknownCount_, &rowBlockSize_, &colBlockSize_,
                    &context_, &INT_ZERO, &INT_ZERO, &localMatrixLeadDim_, &info);
  assert(info == 0);
  
  F77NAME(descinit)(rhsDesc_, &largestDimension_, &rhsCount_, &rhsRowBlockSize_, &rhsColBlockSize_,
                    &context_, &INT_ZERO, &INT_ZERO, &localRhsLeadDim_, &info);
  assert(info == 0);
#endif /* DO_SCALAPACK */
}

int
DistLeastSquareSolver::globalRowIdx(int localRowIdx) const {
#ifdef DO_SCALAPACK
  const int localRowIdx_f = localRowIdx + 1;
  return F77NAME(indxl2g)(&localRowIdx_f, &rowBlockSize_, &localCpuRow_, &INT_ZERO, &rowCpus_) - 1;
#else /* DO_SCALAPACK */
  return 0;
#endif /* DO_SCALAPACK */
}

int
DistLeastSquareSolver::globalColIdx(int localColIdx) const {
#ifdef DO_SCALAPACK
  const int localColIdx_f = localColIdx + 1;
  return F77NAME(indxl2g)(&localColIdx_f, &colBlockSize_, &localCpuCol_, &INT_ZERO, &colCpus_) - 1;
#else /* DO_SCALAPACK */
  return 0;
#endif /* DO_SCALAPACK */
}

int
DistLeastSquareSolver::globalRhsRowIdx(int localRowIdx) const {
#ifdef DO_SCALAPACK
  const int localRowIdx_f = localRowIdx + 1;
  return F77NAME(indxl2g)(&localRowIdx_f, &rhsRowBlockSize_, &localCpuRow_, &INT_ZERO, &rowCpus_) - 1;
#else /* DO_SCALAPACK */
  return 0;
#endif /* DO_SCALAPACK */
}

int
DistLeastSquareSolver::globalRhsRankIdx(int localRankIdx) const {
#ifdef DO_SCALAPACK
  const int localRankIdx_f = localRankIdx + 1;
  return F77NAME(indxl2g)(&localRankIdx_f, &rhsColBlockSize_, &localCpuCol_, &INT_ZERO, &colCpus_) - 1;
#else /* DO_SCALAPACK */
  return 0;
#endif /* DO_SCALAPACK */
}

int
DistLeastSquareSolver::rowHostCpu(int globalRowIdx) const {
#ifdef DO_SCALAPACK
  const int globalRowIdx_f = globalRowIdx + 1;
  return F77NAME(indxg2p)(&globalRowIdx_f, &rowBlockSize_, NULL, &INT_ZERO, &rowCpus_);
#else /* DO_SCALAPACK */
  return 0;
#endif /* DO_SCALAPACK */
}

int
DistLeastSquareSolver::colHostCpu(int globalColIdx) const {
#ifdef DO_SCALAPACK
  const int globalColIdx_f = globalColIdx + 1;
  return F77NAME(indxg2p)(&globalColIdx_f, &colBlockSize_, NULL, &INT_ZERO, &colCpus_);
#else /* DO_SCALAPACK */
  return 0;
#endif /* DO_SCALAPACK */
}

int
DistLeastSquareSolver::rhsRowHostCpu(int globalRowIdx) const {
#ifdef DO_SCALAPACK
  const int globalRowIdx_f = globalRowIdx + 1;
  return F77NAME(indxg2p)(&globalRowIdx_f, &rhsRowBlockSize_, NULL, &INT_ZERO, &rowCpus_);
#else /* DO_SCALAPACK */
  return 0;
#endif /* DO_SCALAPACK */
}

int
DistLeastSquareSolver::rhsRankHostCpu(int globalRankIdx) const {
#ifdef DO_SCALAPACK
  const int globalRankIdx_f = globalRankIdx + 1;
  return F77NAME(indxg2p)(&globalRankIdx_f, &rhsColBlockSize_, NULL, &INT_ZERO, &colCpus_);
#else /* DO_SCALAPACK */
  return 0;
#endif /* DO_SCALAPACK */
}

int
DistLeastSquareSolver::localRowIdx(int globalRowIdx) const {
#ifdef DO_SCALAPACK
  const int globalRowIdx_f = globalRowIdx + 1;
  return F77NAME(indxg2l)(&globalRowIdx_f, &rowBlockSize_, NULL, NULL, &rowCpus_) - 1;
#else /* DO_SCALAPACK */
  return 0;
#endif /* DO_SCALAPACK */
}

int
DistLeastSquareSolver::localColIdx(int globalColIdx) const {
#ifdef DO_SCALAPACK
  const int globalColIdx_f = globalColIdx + 1;
  return F77NAME(indxg2l)(&globalColIdx_f, &colBlockSize_, NULL, NULL, &colCpus_) - 1;
#else /* DO_SCALAPACK */
  return 0;
#endif /* DO_SCALAPACK */
}

int
DistLeastSquareSolver::localRhsRowIdx(int globalRowIdx) const {
#ifdef DO_SCALAPACK
  const int globalRowIdx_f = globalRowIdx + 1;
  return F77NAME(indxg2l)(&globalRowIdx_f, &rhsRowBlockSize_, NULL, NULL, &rowCpus_) - 1;
#else /* DO_SCALAPACK */
  return 0;
#endif /* DO_SCALAPACK */
}

int
DistLeastSquareSolver::localRhsRankIdx(int globalRankIdx) const {
#ifdef DO_SCALAPACK
  const int globalRankIdx_f = globalRankIdx + 1;
  return F77NAME(indxg2l)(&globalRankIdx_f, &rhsColBlockSize_, NULL, NULL, &colCpus_) - 1;
#else /* DO_SCALAPACK */
  return 0;
#endif /* DO_SCALAPACK */
}

void
DistLeastSquareSolver::solve() {
#ifdef DO_SCALAPACK
  const char * trans = "N";
  int info;

  Scalar workspaceQuery;
  F77NAME(pdgels)(trans, &equationCount_, &unknownCount_, &rhsCount_,
                  matrixBuffer_.array(), &INT_ONE, &INT_ONE, matrixDesc_,
                  rhsBuffer_.array(),    &INT_ONE, &INT_ONE, rhsDesc_,
                  &workspaceQuery, &INT_MINUS_ONE, &info);
  assert(info == 0);

  const int lwork = static_cast<int>(workspaceQuery);
  SimpleBuffer<Scalar> workspace(lwork);

  F77NAME(pdgels)(trans, &equationCount_, &unknownCount_, &rhsCount_,
                  matrixBuffer_.array(), &INT_ONE, &INT_ONE, matrixDesc_,
                  rhsBuffer_.array(),    &INT_ONE, &INT_ONE, rhsDesc_,
                  workspace.array(), &lwork, &info);
  assert(info == 0);
#endif /* DO_SCALAPACK */
}
