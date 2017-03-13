#ifndef DIST_LEAST_SQUARE_SOLVER_H
#define DIST_LEAST_SQUARE_SOLVER_H

class Communicator;

#include "SimpleBuffer.h"

// ScaLAPACK-based distributed least-square solver
class DistLeastSquareSolver {
public:
  typedef double Scalar;
 
  // Underlying communicator
  const Communicator * communicator() const { return communicator_; }

  // Cpu topology
  int rowCpus() const { return rowCpus_; }
  int colCpus() const { return colCpus_; }
  int localCpuRow() const { return localCpuRow_; }
  int localCpuCol() const { return localCpuCol_; }

  // Blocking parameters (Solution has same parameters as rhs)
  int matrixRowBlockSize() const { return rowBlockSize_;    }
  int matrixColBlockSize() const { return colBlockSize_;    }
  int rhsRowBlockSize()    const { return rhsRowBlockSize_; }
  int rhsColBlockSize()    const { return rhsColBlockSize_; }

  void blockSizeIs(int size);
  void blockSizeIs(int row, int col, int rhsRow, int rhsCol);

  // Buffer resizing policy
  enum BufferResizePolicy { TIGHT, LOOSE };
  BufferResizePolicy bufferResizePolicy() const { return bufferResizePolicy_; }
  void bufferResizePolicyIs(BufferResizePolicy policy);

  // Problem size
  int equationCount()    const { return equationCount_;    }
  int unknownCount()     const { return unknownCount_;     }
  int largestDimension() const { return largestDimension_; } // = max(equationCount, unknownCount)
  int rhsCount()         const { return rhsCount_;         }

  void problemSizeIs(int eqnCount, int unknownCount, int rhsCount = 1);

  // Local data distribution
  int localRows()         const { return localRows_;         }
  int localCols()         const { return localCols_;         }
  int localRhsRows()      const { return localRhsRows_;      }
  int localSolutionRows() const { return localSolutionRows_; }
  int localRhsCount()     const { return localRhsCount_;     }

  // Global/Local (zero-based) index mapping
  // Local to global
  int globalRowIdx(int localIdx) const;
  int globalColIdx(int localIdx) const;

  int globalRhsRowIdx(int localIdx) const; // Also valid for solution
  int globalRhsRankIdx(int localIdx) const;

  // Global to local: host cpu
  int rowHostCpu(int globalIdx) const;
  int colHostCpu(int globalIdx) const;

  int rhsRowHostCpu(int globalIdx) const; // Also valid for solution
  int rhsRankHostCpu(int globalIdx) const;

  // Global to local: local index on host cpu
  int localRowIdx(int globalIdx) const;
  int localColIdx(int localIdx) const;

  int localRhsRowIdx(int globalIdx) const; // Also valid for solution
  int localRhsRankIdx(int globalIdx) const;

  // Local buffers: Internal column-major ordering, zero-based indexing
  // Local matrix buffer: [localRows by localCols]
  Scalar matrixEntry(int row, int col) const;
  const Scalar * matrixColBuffer(int col) const;
  const Scalar * matrixBuffer() const;
  
  Scalar & matrixEntry(int row, int col);
  Scalar * matrixColBuffer(int col);
  Scalar * matrixBuffer();
  
  int localMatrixLeadDim() const { return localMatrixLeadDim_; }
  
  // Local Rhs/Solution buffer: [max(localRhsRows, localSolutionRows) by localRhsCount]
  Scalar rhsEntry(int row) const; // First right-hand side
  Scalar rhsEntry(int rank, int row) const;
  const Scalar * rhsBuffer(int rank) const;
  const Scalar * rhsBuffer() const;
  
  Scalar & rhsEntry(int row); // First right-hand side
  Scalar & rhsEntry(int rank, int row);
  Scalar * rhsBuffer(int rank);
  Scalar * rhsBuffer();

  int localRhsLeadDim() const { return localRhsLeadDim_; }
  
  // Ctor and dtor
  DistLeastSquareSolver(Communicator * comm, int rowCpus, int colCpus);
  ~DistLeastSquareSolver();

  // Solve
  void solve();

private:
  Communicator * communicator_;
  int blacsHandle_;

  typedef int Context;
  
  int rowCpus_, colCpus_;
  int localCpuRow_, localCpuCol_;
  Context context_;

  int rowBlockSize_, colBlockSize_;
  int rhsRowBlockSize_, rhsColBlockSize_;
  
  typedef int ArrayDesc[9];
 
  BufferResizePolicy bufferResizePolicy_;

  int equationCount_, unknownCount_;
  int largestDimension_, rhsCount_;
  int localRows_, localCols_;
  int localRhsRows_, localSolutionRows_, localRhsCount_;
  int localMatrixLeadDim_, localRhsLeadDim_;
  ArrayDesc matrixDesc_;
  ArrayDesc rhsDesc_;

  SimpleBuffer<Scalar> matrixBuffer_;
  SimpleBuffer<Scalar> rhsBuffer_;

  // Private functions
  void reset();

  static const int DEFAULT_BLOCK_SIZE;
  
  // Adressable constants
  static const int INT_ZERO;
  static const int INT_ONE;
  static const int INT_MINUS_ONE;

  // Disallow copy & assignment
  DistLeastSquareSolver(const DistLeastSquareSolver &);
  DistLeastSquareSolver & operator=(const DistLeastSquareSolver &);
};


/* Helper functions for buffer access */

inline
const DistLeastSquareSolver::Scalar *
DistLeastSquareSolver::matrixBuffer() const {
  return matrixBuffer_.array();
} 

inline
const DistLeastSquareSolver::Scalar *
DistLeastSquareSolver::matrixColBuffer(int col) const {
  return matrixBuffer() + (col * localMatrixLeadDim());
}

inline
DistLeastSquareSolver::Scalar
DistLeastSquareSolver::matrixEntry(int row, int col) const {
  return matrixColBuffer(col)[row];
}

inline
DistLeastSquareSolver::Scalar *
DistLeastSquareSolver::matrixBuffer() {
  return matrixBuffer_.array();
} 

inline
DistLeastSquareSolver::Scalar *
DistLeastSquareSolver::matrixColBuffer(int col) {
  return matrixBuffer() + (col * localMatrixLeadDim());
}

inline
DistLeastSquareSolver::Scalar &
DistLeastSquareSolver::matrixEntry(int row, int col) {
  return matrixColBuffer(col)[row];
}

inline
const DistLeastSquareSolver::Scalar *
DistLeastSquareSolver::rhsBuffer() const {
  return rhsBuffer_.array();
} 

inline
const DistLeastSquareSolver::Scalar *
DistLeastSquareSolver::rhsBuffer(int rank) const {
  return rhsBuffer() + (rank * localRhsLeadDim());
}

inline
DistLeastSquareSolver::Scalar
DistLeastSquareSolver::rhsEntry(int rank, int row) const {
  return rhsBuffer(rank)[row];
}

inline
DistLeastSquareSolver::Scalar
DistLeastSquareSolver::rhsEntry(int row) const {
  return rhsEntry(0, row);
}

inline
DistLeastSquareSolver::Scalar *
DistLeastSquareSolver::rhsBuffer() {
  return rhsBuffer_.array();
} 

inline
DistLeastSquareSolver::Scalar *
DistLeastSquareSolver::rhsBuffer(int rank) {
  return rhsBuffer() + (rank * localRhsLeadDim());
}

inline
DistLeastSquareSolver::Scalar &
DistLeastSquareSolver::rhsEntry(int rank, int row) {
  return rhsBuffer(rank)[row];
}

inline
DistLeastSquareSolver::Scalar &
DistLeastSquareSolver::rhsEntry(int row) {
  return rhsEntry(0, row);
}

#endif /* DIST_LEAST_SQUARE_SOLVER_H */
