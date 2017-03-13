#include "SparseGrid.h"
#include <cmath>
#include <cstdlib>
using namespace std;

//------------------------------------------------------------------------------

template<typename T>
void SparseGrid::tabulate(void (T::*fn)(double *, double *, double *), 
                          T &object, bool restart){

    tabulate(Functor<T>(fn,object), restart);

}

//------------------------------------------------------------------------------

template<typename FnType>
void SparseGrid::tabulate(FnType fn, bool restart){

  messages(1);

  if(!restart)
    initialize(fn);

  bool success = false;
  bool adaptivity = false;
  int currentMultiIndex, addedSubGrids;

  while(nPoints<maxPoints && !success){

    messages(2);
    
    // get next active subgrid
    currentMultiIndex = currentSubGrid(adaptivity);
    if(currentMultiIndex<0){
      fprintf(stdout, "SparseGrid: no more subgrids available\n");
      exit(1);
    }
    active[currentMultiIndex] = false;
    messages(3,currentMultiIndex);

    // resize arrays if necessary
    while(nSubGrids+dim>sizeMultiIndex) resizeMultiIndex();

    addedSubGrids = 0;
    // for each forward neighbour(one in each direction/dimension),
    //   check if it is admissible
    for(int idim=0; idim<dim; idim++){
      if(admissible(currentMultiIndex, idim)){
        messages(4,currentMultiIndex);
        // find its neighbours and index it as a neighbour of its own neighbours
        findNeighbours(currentMultiIndex, idim, addedSubGrids);
        // the forward admissible neighbour is now active
        active[neighbour[currentMultiIndex][idim]] = true;
        addedSubGrids++;
      }
    }
    messages(5,currentMultiIndex);

    // get the hierarchical surplus of the points on the admissible 
    // forward neighbouring subgrids (which become active subgrids)
    for(int forwardAdmissible=0; forwardAdmissible<addedSubGrids; forwardAdmissible++){
      int numNewPoints = integrateForwardAdmissible(nSubGrids, fn);
      nPoints += numNewPoints;
      nSubGrids++;

      if(adaptivity) nAdaptivePoints += numNewPoints;
      messages(6,forwardAdmissible);
    }

    if(nPoints>minPoints)
      success = checkAccuracy();

  }
  messages(10);

}

//------------------------------------------------------------------------------

template<typename FnType>
void SparseGrid::initialize(FnType fn){
// initialization of the data structures for the first subgrid
// in the Clenshaw-Curtis grid, which contains only one point located
// at the center of the hypercube (this is level of refinement 0 in all directions).

  nPoints   = 1;
  nSubGrids = 1;
  for(int idim=0; idim<dim; idim++) multiIndex[0][idim] = 0; 

  double firstPoint[dim]; double scaledCoord[dim]; double res[out];
  for(int idim=0; idim<dim; idim++) 
    firstPoint[idim] = 0.5;
  scale(firstPoint,scaledCoord, 0);
  fn(scaledCoord,res,parameters);
  for(int iout=0; iout<out; iout++){
    surplus[0][iout] = res[iout];
    fnmin[iout]      = res[iout];
    fnmax[iout]      = res[iout];
  }
  
  activeError[0] = surplus[0][0]; //max(surplus[0]);
  for(int iout=1; iout<out; iout++)
    if(activeError[0]<surplus[0][iout]) activeError[0] = surplus[0][iout];
  activeCost[0]  = 1.0;
  activeHeapError.insert(0,activeError);
  activeHeapCost.insert(0,activeCost);
  active[0] = true;

  for(int iout=0; iout<out; iout++)
    error[0][iout] = surplus[0][iout];

  // no backward neighbour for subGrid multiIndex[0]
  for(int isize=0; isize<sizeMultiIndex; isize++)
    for(int neigh=0; neigh<2*dim; neigh++)
      neighbour[isize][neigh] = -1;
}

//------------------------------------------------------------------------------

template<typename FnType>
int SparseGrid::integrateForwardAdmissible(const int newSubGridIndex,
                            FnType fn){

  int nPointsSubGrid;
  double ** subGrid = generateSubGrid(newSubGridIndex,nPointsSubGrid);
  evaluateFunctionOnGrid(subGrid, nPointsSubGrid, fn);
  evaluatePreviousInterpolation(subGrid, nPointsSubGrid);

  updateError(nPointsSubGrid);
  updateCost();
  activeHeapError.insert(newSubGridIndex,activeError);
  activeHeapCost.insert(newSubGridIndex,activeCost);

  for (int i=0; i<nPointsSubGrid; i++) delete [] subGrid[i];
  delete [] subGrid;
  return nPointsSubGrid;

}

//------------------------------------------------------------------------------

template<typename FnType>
void SparseGrid::evaluateFunctionOnGrid(double **subGrid,
                                        const int nPointsSubGrid, FnType fn){

  while(nPoints+nPointsSubGrid > sizeSurplus) resizeSurplus();

  double scaledCoord[dim];
  double res[out];
  for(int iPts=0; iPts<nPointsSubGrid; iPts++){
    scale(subGrid[iPts],scaledCoord, 0);
    fn(scaledCoord,res,parameters);
    for(int iout=0; iout<out; iout++)
      surplus[nPoints+iPts][iout] = res[iout];
  }
  
  bounds(nPointsSubGrid); //find the min and max values of the target function
                          //on the grid points

}

//------------------------------------------------------------------------------

template<typename T>
void SparseGrid::test(void (T::*fn)(double *, double *, double *), T &object,
                      int type, int *number, double *param){

    test(Functor<T>(fn,object), type, number, param);

}

//------------------------------------------------------------------------------

template<typename FnType>
void SparseGrid::test(FnType fn, int type, int *number, double *param){
	
  if(nPoints == 0) return;
  int numTestPoints = 0;
  int *numPointsDim = new int[dim];
  int *numCumPointsDim = new int[dim];
  double **coordDim = new double *[dim];
  
  switch(type){
//0 single point testing is done by the user from outside this class

//1 random testing in the domain ('number' random values in each dimension):
    case 1:{
      fprintf(stdout, "\n# The sparse grid is now being tested at random points...\n");
      numTestPoints = static_cast<int>(pow(static_cast<double>(number[0]),dim)+0.1);
      for(int idim=0; idim<dim; idim++){
      	numPointsDim[idim] = number[0];
        coordDim[idim] = new double[number[0]];
        numCumPointsDim[idim] = static_cast<int>(pow(static_cast<double>(number[0]),idim+1)+0.1);
      }
      for(int i=0; i<number[0]; i++)
      	for(int idim=0; idim<dim; idim++)
          coordDim[idim][i] = range[idim][0] + static_cast<double>(rand())/RAND_MAX * (range[idim][1]-range[idim][0]);
       break;
    }

//2 linear testing in the domain ('number' evenly spaces between values in each dimension):
    case 2:{
      fprintf(stdout, "\n# The sparse grid is now being tested at uniformly spaced points...\n");
      numTestPoints = static_cast<int>(pow(static_cast<double>(number[0]),dim)+0.1);
      double spacing[dim];
      for(int idim=0; idim<dim; idim++){
      	numPointsDim[idim] = number[0];
        coordDim[idim] = new double[number[0]];
        spacing[idim] = (range[idim][1]-range[idim][0])/number[0];
        numCumPointsDim[idim] = static_cast<int>(pow(static_cast<double>(number[0]),idim+1)+0.1);
        for(int i=0; i<number[0]; i++)
          coordDim[idim][i] = range[idim][0] + 0.5*spacing[idim] + i*spacing[idim];
      }
      break;
    }
    
//3 grid testing:
    case 3:{
      numTestPoints = 1;
      for(int idim=0; idim<dim; idim++){
      	double spacing = (range[idim][1]-range[idim][0])*pow(2.0, -number[idim]);
      	
      	if(number[idim] == 0) numPointsDim[idim] = 1;
      	else                  numPointsDim[idim] = static_cast<int>(pow(2.0,number[idim])+1.1);
      	numTestPoints *= numPointsDim[idim];
      	numCumPointsDim[idim] = numTestPoints;
      	
      	coordDim[idim] = new double[numPointsDim[idim]];
      	
      	if(numPointsDim[idim] == 1) coordDim[idim][0] = 0.5*(range[idim][0]+range[idim][1]);
      	else{
          coordDim[idim][0] = range[idim][0];
          for(int i=1; i<numPointsDim[idim]-1; i++)
            coordDim[idim][i] = coordDim[idim][i-1] + spacing;
          coordDim[idim][numPointsDim[idim]-1] = range[idim][1];
      	}
      }
      break;
    }
    // coordinates in each dimension have been set up.
  }

  double **coord = new double *[numTestPoints];
  for(int i=0; i<numTestPoints; i++)
    coord[i] = new double[dim];
  tensorize(coord,coordDim,numPointsDim,numCumPointsDim);
  
  // computing results with present tabulation and comparing with exact results.
  double **exact  = new double *[numTestPoints];
  double **tabul  = new double *[numTestPoints];
  for(int i=0; i<numTestPoints; i++){
    exact[i]  = new double[out];
    tabul[i]  = new double[out];
  }

  for(int i=0; i<numTestPoints; i++){
  	fn(coord[i],exact[i],param);
  }

  interpolate(numTestPoints,coord,tabul);
  
  // print to screen results
  fprintf(stdout, "# For each test point, the approximated value obtained by the sparse grid is compared with the actual computed value.\n");
  fprintf(stdout, "# The relative tolerance has been hard-coded to 1e-6.\n# Below are the pairs of test points and values that do not satisfy this criterion.\n");
  int numErrors = 0;
  for(int i=0; i<numTestPoints; i++){
    double error = (exact[i][0] - tabul[i][0])/exact[i][0];
    if(fabs(error) > 1e-6){
      numErrors++;
      fprintf(stdout, "test point %d ( ", i);
      for(int idim=0; idim<dim; idim++)
        fprintf(stdout, "%e ", coord[i][idim]);
      fflush(stdout);
      fprintf(stdout, ") \n    -- exact = %e -- approx = %e ", exact[i][0], tabul[i][0]);
      if(exact[i][0] != 0.0)
        fprintf(stdout, "(err = %e)\n", (exact[i][0] - tabul[i][0])/exact[i][0]);
      else fprintf(stdout, "\n");
        fprintf(stdout, "\n");
    }
  }
  if(numErrors == 0)
    fprintf(stdout, "# All test points satisfied the tolerance\n");

}
