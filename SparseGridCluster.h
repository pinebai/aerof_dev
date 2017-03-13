#ifndef SPARSEGRIDCLUSTER_H_
#define SPARSEGRIDCLUSTER_H_

class SparseGrid;
class Communicator;
struct SparseGridData;

class SparseGridCluster
{
private:
  int dim_;
  int numSparseGrids_;
  SparseGrid *sparseGrids_;

public:
  SparseGridCluster();
  virtual ~SparseGridCluster();
  
  template <typename T>
  void generate(SparseGridData &data, double *param,
                void (T::*fn)(double *, double *, double *), T &object,
                const char *filename,
                const double *refIn, const double *refOut,
                const Communicator *com);

  void readFromFile(const int numFiles, const double *refIn, const double *refOut, 
                    const char *filename, int outputRangeFlag=0);
  
  int interpolate(const int numRes, double **coord, double **res);
  int interpolateGradient(const int numRes, double **coord, double **res);
};
#ifdef TEMPLATE_FIX
#include "SparseGridCluster.C"
#endif

#endif /*SPARSEGRIDCLUSTER_H_*/
