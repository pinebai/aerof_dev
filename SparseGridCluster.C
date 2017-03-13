#include "SparseGridCluster.h"

#include "SparseGrid.h"
#include "IoData.h"
#include "Communicator.h"
#include <cstring>
#include <cstdio>
#include <cmath>

template<typename T>
void SparseGridCluster::generate(SparseGridData &data, double *param,
                                 void (T::*fn)(double *, double *, double *), T &object,
                                 const char *filename,
                                 const double *refIn, const double *refOut,
                                 const Communicator *com)
{
  bool debugInfo = false;
  dim_ = data.numInputs;

  int *nGridsDim    = new int[dim_];
  int *nCumGridsDim = new int[dim_];
  int *nPtsDim      = new int[dim_];
  for(int idim=0; idim<dim_; idim++){
    nGridsDim[idim]    = data.numDomainDim[idim];
    nCumGridsDim[idim] = (idim==0)? nGridsDim[idim] : nCumGridsDim[idim-1]*nGridsDim[idim];
    nPtsDim[idim]      = data.numDomainDim[idim]+1;
  }
  numSparseGrids_ = nCumGridsDim[dim_-1];

  double **coordDim = new double  *[dim_];
  for(int idim=0; idim<dim_; idim++){
    coordDim[idim] = new double[nPtsDim[idim]];
    double spacing = (data.range[idim][1]-data.range[idim][0])/data.numDomainDim[idim];
    for(int i=0; i<nPtsDim[idim]; i++){
      if(i==0) coordDim[idim][0] = data.range[idim][0];
      else if(i==nPtsDim[idim]) coordDim[idim][i] = data.range[idim][1];
      else coordDim[idim][i] = coordDim[idim][i-1] + spacing;
    }
  }
  
  SparseGridData subData(data); // check what it copies since default copy constructor
  int filenameLength = strlen(filename);
  char ifilename[filenameLength+10];
  
  int *index = new int[dim_];
  for(int idim=0; idim<dim_; idim++)
    index[idim] = 0;

  int sgPerCpu = numSparseGrids_/com->size();
  if(debugInfo) 
    fprintf(stdout, "%d cpus for %d sparse grids => %d or %d sparse grids per cpu\n", com->size(), numSparseGrids_, sgPerCpu, sgPerCpu+1);
  for(int k=0; k<sgPerCpu+1; k++){
    int igrid = com->cpuNum() + k*com->size();
    if(igrid>=numSparseGrids_) break;
    //fprintf(stdout, "taking care of sparse grid # %d\n", igrid);
    int moveForward = 0;
    if(k==0) moveForward = com->cpuNum();
    else     moveForward = com->size();

    for(int move=0; move<moveForward; move++){
      index[0]++;
      for(int idim=0; idim<dim_; idim++)
        if(index[idim]==nGridsDim[idim]){
          index[idim+1]++;
          index[idim] = 0;
        }
    }

    if(debugInfo) {
      fprintf(stdout, "index = [ ");
      for(int idim=0; idim<dim_; idim++)
        fprintf(stdout, "%d ", index[idim]);
      fprintf(stdout, "]\n");
    }

    //fprintf(stdout, "subData.range = ");
    for(int idim=0; idim<dim_; idim++){
      subData.range[idim][0] = coordDim[idim][index[idim]];
      subData.range[idim][1] = coordDim[idim][index[idim]+1];
      //fprintf(stdout, "[%e %e] ", subData.range[idim][0],subData.range[idim][1]);
    }
    //fprintf(stdout, "\n");

    SparseGrid sparseGrid(subData,param,refIn,refOut);
    sparseGrid.tabulate(fn, object);
    sprintf(ifilename, "%s%d", filename, igrid+1);
    sparseGrid.printToFile(refIn, refOut, ifilename);
  }


  /*for(int igrid=0; igrid<numSparseGrids_; igrid++){
    index[0]++;
    for(int idim=0; idim<dim_; idim++)
      if(index[idim]==nGridsDim[idim]){
        index[idim+1]++;
        index[idim] = 0;
      }

    for(int idim=0; idim<dim_; idim++){
      subData.range[idim][0] = coordDim[idim][index[idim]];
      subData.range[idim][1] = coordDim[idim][index[idim]+1];
    }

    sparseGrids_[igrid].initialize(subData, param, refIn, refOut);
    sparseGrids_[igrid].tabulate(fn, object);
    sprintf(ifilename, "%s%d", filename, igrid+1);
    sparseGrids_[igrid].printToFile(refIn, refOut, ifilename);
  }*/
  
}
