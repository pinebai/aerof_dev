#include "SparseGridCluster.h"

#include "SparseGrid.h"
#include <cstring>
#include <cstdio>

//------------------------------------------------------------------------------
SparseGridCluster::SparseGridCluster()
{
  dim_            = 0;
  numSparseGrids_ = 0;
  sparseGrids_    = 0;
}
//------------------------------------------------------------------------------
SparseGridCluster::~SparseGridCluster()
{
  delete [] sparseGrids_;
}
//------------------------------------------------------------------------------
void SparseGridCluster::readFromFile(const int numFiles, 
                                     const double *refIn, const double *refOut,
                                     const char *filename,
                                     int outputRangeFlag)
{
  numSparseGrids_ = numFiles;
  sparseGrids_ = new SparseGrid[numSparseGrids_];
  
  int filenameLength = strlen(filename);
  char ifilename[filenameLength+10];
  for(int i=0; i<numSparseGrids_; i++){
    sprintf(ifilename, "%s%d", filename, i+1);
    sparseGrids_[i].readFromFile(refIn, refOut, ifilename, outputRangeFlag);
  }	
  dim_ = sparseGrids_[0].getDim();
}
//------------------------------------------------------------------------------
int SparseGridCluster::interpolate(const int numRes, double **coord, double **res)
{
  bool inSparseGrid = false;
  for(int iRes=0; iRes<numRes; iRes++){
    for(int iGrid=0; iGrid<numSparseGrids_; iGrid++){
      inSparseGrid = sparseGrids_[iGrid].contains(coord[iRes]);
      if(inSparseGrid){
        sparseGrids_[iGrid].interpolate(1, &(coord[iRes]), &(res[iRes]));
        break;
      }
    }
    if(!inSparseGrid){
      fprintf(stdout, "*** Warning: coord[%d] = ( ", iRes);
      for(int idim=0; idim<dim_; idim++)
        fprintf(stdout, "%e ", coord[iRes][idim]);
      fprintf(stdout, ") is out of range of all SparseGrids.\n");
      //exit(1);
      return 0;
    }
  }
  return 1;
	
}
//------------------------------------------------------------------------------
int SparseGridCluster::interpolateGradient(const int numRes, double **coord, double **res)
{
  bool inSparseGrid = false;
  for(int iRes=0; iRes<numRes; iRes++){
    for(int iGrid=0; iGrid<numSparseGrids_; iGrid++){
      inSparseGrid = sparseGrids_[iGrid].contains(coord[iRes]);
      if(inSparseGrid){
        sparseGrids_[iGrid].interpolateGradient(1, &(coord[iRes]), &(res[iRes]));
        break;
      }
    }
    if(!inSparseGrid){
      fprintf(stdout, "coord[%d] = ( ", iRes);
      for(int idim=0; idim<dim_; idim++)
        fprintf(stdout, "%e ", coord[iRes][idim]);
      fprintf(stdout, ") is out of range of all SparseGrids.\n");
      return 1;
    }
  }
 
  return 0;
	
}
//------------------------------------------------------------------------------
