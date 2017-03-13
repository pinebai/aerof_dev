#include "SparseGrid.h"
#include <cmath>
#include <iostream>


//------------------------------------------------------------------------------

SparseGrid::SparseGrid(){

  parameters = 0;

  verbose = 0;

  dim = 0;
  out = 0;

  nPoints     = 0;
  sizeSurplus = 0;
  surplus     = 0;

  nSubGrids      = 0;
  sizeMultiIndex = 0;
  multiIndex     = 0;

  maxPoints = 100;
  minPoints = 100;
  absAccuracy = -1.0;
  relAccuracy = -1.0;
  range = 0;
  dimAdaptDegree = 0.0;

  nAdaptivePoints = 0;
  active = 0;

  activeError = 0;
  activeCost  = 0;

  neighbour = 0;
  error = 0;
  fnmin = 0;
  fnmax = 0;

  logMap = 0;

}

//------------------------------------------------------------------------------

SparseGrid::~SparseGrid(){

  parameters = 0;

  for(int i=0; i<sizeSurplus; i++) delete [] surplus[i];
  delete [] surplus;
  for(int i=0; i<sizeMultiIndex; i++) delete [] multiIndex[i];
  delete [] multiIndex;
  delete [] range;
  delete [] active;
  delete [] activeError;
  delete [] activeCost;
  if(neighbour){ // not initialized if sparse grid read from a file
    for(int i=0; i<sizeMultiIndex; i++) delete [] neighbour[i];
    delete [] neighbour;
  }
  if(error){ // not initialized if sparse grid read from a file
    for(int i=0; i<sizeMultiIndex; i++) delete [] error[i];
    delete [] error;
  }
  delete [] fnmin;
  delete [] fnmax;
  
  for(int i=0; i<dim; i++) delete logMap[i];
  delete [] logMap;

}

//------------------------------------------------------------------------------

SparseGrid::SparseGrid(SparseGridData &data, double *param, 
                       const double *refIn, const double *refOut){
  initialize(data,param,refIn,refOut);
}

//------------------------------------------------------------------------------

void SparseGrid::initialize(SparseGridData &data, double *param, 
                       const double *refIn, const double *refOut){

  parameters = param;

  dim = data.numInputs;
  out = data.numOutputs;
  verbose = data.verbose;

  maxPoints = data.maxPoints;
  minPoints = data.minPoints;
  absAccuracy = data.absAccuracy;
  relAccuracy = data.relAccuracy;
  
  logMap = new LogarithmicMapping *[dim];
  range = new Range[dim];
  for(int idim=0; idim<dim; idim++){
    if(data.mapBaseValue[idim] > 1.0){
      logMap[idim] = new LogarithmicMapping(data.mapBaseValue[idim]);
    }
    else
      logMap[idim] = 0;
  	
    range[idim][0] = data.range[idim][0];
    range[idim][1] = data.range[idim][1];
    // WARNING: logMap is applied later, in scaleGrid(...)
  }
  dimAdaptDegree = data.dimAdaptDegree;

  nPoints         = 0;
  sizeSurplus     = maxPoints;
  surplus         = new double *[sizeSurplus];
  for(int i=0; i<sizeSurplus; i++) surplus[i] = new double[out];

  nSubGrids       = 0;
  sizeMultiIndex  = 10;
  multiIndex      = new int *[sizeMultiIndex];
  for(int i=0; i<sizeMultiIndex; i++) multiIndex[i] = new int[dim];

  nAdaptivePoints = 0;
  active          = new bool[sizeMultiIndex];
  activeError  = new double[sizeMultiIndex];
  activeCost   = new double[sizeMultiIndex];

  neighbour       = new int *[sizeMultiIndex];
  for(int i=0; i<sizeMultiIndex; i++) neighbour[i] = new int[2*dim];
  
  error           = new double *[sizeMultiIndex];
  for(int i=0; i<sizeMultiIndex; i++) error[i] = new double[out];
  
  fnmin = new double[out];
  fnmax = new double[out];
  
  scaleGrid(refIn, refOut, 1);

}

//------------------------------------------------------------------------------

void SparseGrid::newParameters(SparseGridData &data, double *param){

  parameters = param;

  dim = data.numInputs;
  out = data.numOutputs;
  verbose = data.verbose;

  maxPoints = data.maxPoints;
  minPoints = data.minPoints;
  absAccuracy = data.absAccuracy;
  relAccuracy = data.relAccuracy;
  
}
  
//------------------------------------------------------------------------------

int SparseGrid::currentSubGrid(bool &adaptivity){
// pick next subgrid: choose from the two available heaps

  bool done = false;
  int current  = -1;

  adaptivity = false;

  while(!done){
    current = -1; adaptivity = false;
    if(nAdaptivePoints>=dimAdaptDegree*nPoints) // non-adaptive
      current = activeHeapCost.pop(activeCost);

    if(current<0){ // adaptive
      current = activeHeapError.pop(activeError);
      if(current>=0) adaptivity = true;
    }

    if(current<0)
      current = activeHeapCost.pop(activeCost);

    if(current<0) return current;

    // make sure it is really active
    if(active[current]) done = true;

  }

  return current;

}

//------------------------------------------------------------------------------

bool SparseGrid::admissible(const int currentMultiIndex,
                            const int forwardDir){

  // currentMultiIndex has multi-index multiIndex[currentMultiIndex]
  // To check the admissibility of the newMultiIndex 
  // (defined by the forward neighbour of currentMultiIndex in the 
  //  direction forwardDir, see sketch SparseGrid::findNeighbours),
  // all its backward neighbours must be inactive (integrated) subgrids.
  // In particular, its backward neighbour currentMultiIndex
  //   refers to an inactive subgrid. 

  int backwardOfCurrent, backwardOfNew;

  for(int idim=0; idim<dim; idim++){
    if(forwardDir == idim) continue;
    backwardOfCurrent = neighbour[currentMultiIndex][dim+idim];
    if(backwardOfCurrent<0){ // already at lowest level in direction idim
      if(multiIndex[currentMultiIndex][idim]!=0){
        fprintf(stdout, "*** Error: contradiction detected (no backward neighbours but not at lowest level) ...\n");
        exit(1);
      }
      continue;
    }
    backwardOfNew = neighbour[backwardOfCurrent][forwardDir];
    if(backwardOfNew<0)       return false; // not part of sparse grid 
    if(active[backwardOfNew]) return false; // part of sparse grid, but active
    
  }

  return true;

}

//------------------------------------------------------------------------------

void SparseGrid::findNeighbours(const int currentMultiIndex,
                                const int forwardDir, const int addedSubGrids){
  int backwardOfCurrent, backwardOfNew;
  const int newMultiIndex = nSubGrids+addedSubGrids;

//    ____________________
//   |         |          |
//   | Current |   New    |  idim
//   |_________|__________|   ^
//   |         |          |   |
//   |  bOfC   |  bOfN    |   |
//   |_________|__________|   ----> forwardDir
//

  // find the neighbours
  for(int idim=0; idim<dim; idim++){
    if(forwardDir == idim){
      neighbour[currentMultiIndex][forwardDir    ] = newMultiIndex;
      neighbour[newMultiIndex    ][forwardDir+dim] = currentMultiIndex;
    }else{
      backwardOfCurrent = neighbour[currentMultiIndex][dim+idim];
      if(!(backwardOfCurrent<0)){ 
        backwardOfNew = neighbour[backwardOfCurrent][forwardDir];
        neighbour[backwardOfNew][idim    ] = newMultiIndex;
        neighbour[newMultiIndex][idim+dim] = backwardOfNew;
      }
    }
  }

  // get the levels of refinement of the new active subgrids
  for(int idim=0; idim<dim; idim++){
    backwardOfNew = neighbour[newMultiIndex][idim+dim];
    if(backwardOfNew>=0)
      multiIndex[newMultiIndex][idim] = multiIndex[backwardOfNew][idim]+1;
    else multiIndex[newMultiIndex][idim] = 0;
  }

}

//------------------------------------------------------------------------------

double **SparseGrid::generateSubGrid(const int newSubGrid,
                                     int &nPointsSubGrid) const{

  // find number of points in subgrid (Clenshaw-Curtis type)
  // and coordinates of the points in each direction.
  nPointsSubGrid = 1;
  int nPointsDim[dim]; // number of points in each dimension
  int nCumulatedPointsDim[dim]; // number of cumulated points in each dimension
  double **coordDim = new double*[dim];
  for(int i=0; i<dim; i++){
    if(multiIndex[newSubGrid][i] == 0){
      nPointsDim[i] = 1;
      coordDim[i] = new double[1];
      coordDim[i][0] = 0.5;
    }
    else if(multiIndex[newSubGrid][i] < 3){
      nPointsDim[i] = 2;
      coordDim[i] = new double[2];
      if(multiIndex[newSubGrid][i] == 1){
        coordDim[i][0] = 0.0;
        coordDim[i][1] = 1.0;
      }else{
        coordDim[i][0] = 0.25;
        coordDim[i][1] = 0.75;
      }
    }
    else{
      nPointsDim[i] = 1<<(multiIndex[newSubGrid][i]-1); //static_cast<int>(pow(2.0,multiIndex[newSubGrid][i]-1)+0.1);
      coordDim[i] = new double[nPointsDim[i]];
      for(int j=0; j<nPointsDim[i]; j++)
        coordDim[i][j] = (2.0*j+1.0)/(2.0*nPointsDim[i]);
    }

    nCumulatedPointsDim[i] = nPointsDim[i];
    nPointsSubGrid *= nPointsDim[i];

    if(i>0) nCumulatedPointsDim[i] *= nCumulatedPointsDim[i-1];
    if(i==dim-1 && nCumulatedPointsDim[i]!=nPointsSubGrid)
      fprintf(stdout, "*** Error: the number of points in this subgrid is not correct\n");
  }

  // tensorization of the coordinates is returned in res
  double **res = new double *[nPointsSubGrid];
  for(int i=0; i<nPointsSubGrid; i++)
    res[i] = new double[dim];
  tensorize(res,coordDim,nPointsDim,nCumulatedPointsDim);

  for(int i=0; i<dim; i++) delete [] coordDim[i];
  delete [] coordDim;

  return res;

}

//------------------------------------------------------------------------------

void SparseGrid::tensorize(double **res, double ** coordDim,
                           const int *nPtsDim, const int *nCumulatedPtsDim) const{

  for(int newdim=0; newdim<dim; newdim++){

    if(newdim==0){
      for(int iPts=0; iPts<nPtsDim[newdim]; iPts++)
        res[iPts][newdim] = coordDim[newdim][iPts];

    }else{
      for(int iPts=0; iPts<nPtsDim[newdim]; iPts++)
        for(int iPtsRepeat=0; iPtsRepeat<nCumulatedPtsDim[newdim-1]; iPtsRepeat++)
          res[nCumulatedPtsDim[newdim-1]*iPts+iPtsRepeat][newdim] = coordDim[newdim][iPts];
  
      for(int iPts=0; iPts<nCumulatedPtsDim[newdim-1]; iPts++)
        for(int iPtsRepeat=0; iPtsRepeat<nPtsDim[newdim]; iPtsRepeat++)
          for(int otherdim=newdim-1; otherdim>=0; otherdim--)
            res[iPts+iPtsRepeat*nCumulatedPtsDim[newdim-1]][otherdim] = res[iPts][otherdim];
    }
  }

}

//------------------------------------------------------------------------------

void SparseGrid::scale(const double *subGrid, double *scaledCoord, const int op) const{

  if(op==0) // from SparseGrid scale (0,1) to physical non-dimensional scale
    for(int i=0; i<dim; i++){
      scaledCoord[i] = subGrid[i]*(range[i][1]-range[i][0])+range[i][0];
      if(logMap[i]) scaledCoord[i] = logMap[i]->invMap(scaledCoord[i]);
    }
  else if(op==1) // from log scale (rangemin,rangemax) to SparseGrid scale (0,1)
    for(int i=0; i<dim; i++)
      scaledCoord[i] = (subGrid[i]-range[i][0])/(range[i][1]-range[i][0]);
  else if(op==2) // from non-dimensional physical scale to SparseGrid scale (0,1)
    for(int i=0; i<dim; i++){
      double temp = subGrid[i];
      if(logMap[i]) temp = logMap[i]->map(subGrid[i]);
      scaledCoord[i] = (temp-range[i][0])/(range[i][1]-range[i][0]);
    }
    
  /*if(op==0)
    for(int i=0; i<dim; i++)
      if(range[i][1]>range[i][0])
        scaledCoord[i] = subGrid[i]*(range[i][1]-range[i][0])+range[i][0];
      else scaledCoord[i] = range[i][0];
  else if(op==1)
    for(int i=0; i<dim; i++)
      if(range[i][1]>range[i][0])
        scaledCoord[i] = (subGrid[i]-range[i][0])/(range[i][1]-range[i][0]);
      else scaledCoord[i] = 0.0;
  */
}

//------------------------------------------------------------------------------

void SparseGrid::bounds(const int nPointsSubGrid){

  double temp;

  for(int iout=0; iout<out; iout++){
    for(int iPts=0; iPts<nPointsSubGrid; iPts++){
      temp = surplus[nPoints+iPts][iout]; //surplus data struct does not contain the
                                          //surplus of the function yet, it only contains
                                          //the function values on the new grid.
      if(temp<fnmin[iout]) fnmin[iout] = temp;
      if(temp>fnmax[iout]) fnmax[iout] = temp;
    }
  }

}

//------------------------------------------------------------------------------

void SparseGrid::evaluatePreviousInterpolation(double **subGrid,
                                               const int nPointsSubGrid){

  double temp[out];
  for(int iPts=0; iPts<nPointsSubGrid; iPts++){
    singleInterpolation(subGrid[iPts], temp);
    for(int iout=0; iout<out; iout++) {
      surplus[nPoints+iPts][iout] -= temp[iout];
    }
  }

}

//------------------------------------------------------------------------------

void SparseGrid::updateError(const int nPointsSubGrid){

  double indicator[out]; double temp;

  // computes the "errors" given by the surpluses for a subgrid in each output
  for(int iout=0; iout<out; iout++){
    indicator[iout] = 0.0;
    for(int iPts=0; iPts<nPointsSubGrid; iPts++){
      temp = surplus[nPoints+iPts][iout];
      indicator[iout] += ((temp>0) ? temp : -temp);
    }
    indicator[iout] /= nPointsSubGrid;
  }

  activeError[nSubGrids] = 0.0; // = max(indicator)
  for(int iout=0; iout<out; iout++)
    if(indicator[iout] > activeError[nSubGrids])
      activeError[nSubGrids] = indicator[iout];

  // computes the "errors" for the accuracy check
  for(int iout=0; iout<out; iout++){
    error[nSubGrids][iout] = 0.0;
    for(int iPts=0; iPts<nPointsSubGrid; iPts++){
      temp = surplus[nPoints+iPts][iout];
      temp = (temp>0) ? temp : -temp;
      if(temp>error[nSubGrids][iout]) 
        error[nSubGrids][iout] = temp;
    }
  }
      
}

//------------------------------------------------------------------------------

void SparseGrid::updateCost(){

  activeCost[nSubGrids] = 0.0;
  for(int idim=0; idim<dim; idim++)
    activeCost[nSubGrids] += static_cast<double>(multiIndex[nSubGrids][idim]);

  activeCost[nSubGrids] = 1.0/activeCost[nSubGrids];

}

//------------------------------------------------------------------------------

bool SparseGrid::checkAccuracy(){

  double temp;
  double maxsurplus[out];
  for(int iout=0; iout<out; iout++) maxsurplus[iout] = 0.0;
  for(int iout=0; iout<out; iout++){
    for(int iactiveErr=0; iactiveErr<activeHeapError.size(); iactiveErr++){
      if(active[activeHeapError[iactiveErr]]){
        temp = error[activeHeapError[iactiveErr]][iout];
        if(temp>maxsurplus[iout]) maxsurplus[iout] = temp;
      }
    }

    for(int iactiveCost=0; iactiveCost<activeHeapCost.size(); iactiveCost++){
      if(active[activeHeapCost[iactiveCost]]){
        temp = error[activeHeapCost[iactiveCost]][iout];
        if(temp>maxsurplus[iout]) maxsurplus[iout] = temp;
      }
    }
  }

  actualAbsAcc = 0.0;
  actualRelAcc = 0.0;
  for(int iout=0; iout<out; iout++){
    if(actualAbsAcc<maxsurplus[iout]) actualAbsAcc = maxsurplus[iout];
    if(fnmax[iout]>fnmin[iout])
      if(actualRelAcc<maxsurplus[iout]/(fnmax[iout]-fnmin[iout])) 
        actualRelAcc = maxsurplus[iout]/(fnmax[iout]-fnmin[iout]);
  }

  if(actualAbsAcc<=absAccuracy || actualRelAcc<=relAccuracy){
    messages(7);
    return true;
  }
  else return false;

}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//assumes that each coord is between 0 and 1.
void SparseGrid::singleInterpolation(const double *coord, double *output) const{

  int firstSurplus = 0; // points to first surplus of the considered subgrid
  int nPointsSubGrid;
  int surplusLocalCoord[dim], contributingSurplus;
  for(int iout=0; iout<out; iout++) output[iout] = 0.0;

  messages(8);

  for(int subGrid=0; subGrid<nSubGrids; subGrid++){

    // compute subGrid data structures
    nPointsSubGrid = 1;
    int nPointsDim[dim];
    int nCumulatedPointsDim[dim];
    for(int idim=0; idim<dim; idim++){
      if(multiIndex[subGrid][idim] == 0)
        nPointsDim[idim] = 1;
      else if(multiIndex[subGrid][idim] < 3)
        nPointsDim[idim] = 2;
      else
        nPointsDim[idim] = 1<<(multiIndex[subGrid][idim]-1); // static_cast<int>(pow(2.0,multiIndex[subGrid][idim]-1)+0.1);
  
      nCumulatedPointsDim[idim] = nPointsDim[idim];
      nPointsSubGrid *= nPointsDim[idim];
      if(idim>0) nCumulatedPointsDim[idim] *= nCumulatedPointsDim[idim-1];
    }

    // interpolation: value of the basis function at the considered subgrid (basisFnVal)
    //              : detection of the corresponding surplus in numbering of the
    //                  considered subgrid (surplusLocalCoord).
    double basisFnVal = 1.0;
    for(int idim=0; idim<dim; idim++){
      if(multiIndex[subGrid][idim] == 1){
        if(coord[idim] == 1.0) surplusLocalCoord[idim] = 1;
        else{
          int xp = static_cast<int>(floor(coord[idim]*2));
          if(xp == 0) basisFnVal *= 2.0*(0.5 - coord[idim]);
          else        basisFnVal *= 2.0*(coord[idim] - 0.5);
          surplusLocalCoord[idim] = xp;
        }
      }else if(multiIndex[subGrid][idim] == 0) surplusLocalCoord[idim] = 0;
      else if(coord[idim] == 0.0){
        basisFnVal = 0.0;
        break;
      }else{
        double scale = 1<<multiIndex[subGrid][idim]; // pow(2.0,multiIndex[subGrid][idim]);
        int xp = static_cast<int>(floor(coord[idim] * scale / 2.0));
        basisFnVal *= (1.0 - scale*fabs(coord[idim]-(2.0*static_cast<double>(xp)+1.0)/scale));
        surplusLocalCoord[idim] = xp;
      }

      if(basisFnVal == 0.0) break;

    }

    // contributing surplus in the considered subgrid is added to result
    if(basisFnVal > 0.0){
      // position of the contributing surplus in the array of surpluses
      contributingSurplus = firstSurplus + surplusLocalCoord[0];
      for(int idim=1; idim<dim; idim++)
        contributingSurplus += nCumulatedPointsDim[idim-1]*surplusLocalCoord[idim];
      // contribution added
      for(int iout=0; iout<out; iout++)
        output[iout] += basisFnVal*surplus[contributingSurplus][iout];
    }
    firstSurplus += nPointsSubGrid;

  }

}

//assumes that each coord is between 0 and 1.
void SparseGrid::singleInterpolationGradient(const double *coord, double *output) const{

  int firstSurplus = 0; // points to first surplus of the considered subgrid
  int nPointsSubGrid;
  int surplusLocalCoord[dim], contributingSurplus;
  for(int iout=0; iout<out; iout++) output[iout] = 0.0;

  messages(8);

  for(int subGrid=0; subGrid<nSubGrids; subGrid++){

    // compute subGrid data structures
    nPointsSubGrid = 1;
    int nPointsDim[dim];
    int map[dim];
    int nCumulatedPointsDim[dim];
    for(int idim=0; idim<dim; idim++){
      if(multiIndex[subGrid][idim] == 0)
        nPointsDim[idim] = 1;
      else if(multiIndex[subGrid][idim] < 3)
        nPointsDim[idim] = 2;
      else
        nPointsDim[idim] = 1<<(multiIndex[subGrid][idim]-1); //static_cast<int>(pow(2.0,multiIndex[subGrid][idim]-1)+0.1);

      nCumulatedPointsDim[idim] = nPointsDim[idim];
      nPointsSubGrid *= nPointsDim[idim];
      if(idim>0) nCumulatedPointsDim[idim] *= nCumulatedPointsDim[idim-1];
      map[idim] = 0;
    }

    // interpolation: value of the basis function at the considered subgrid (basisFnVal)
    //              : detection of the corresponding surplus in numbering of the
    //                  considered subgrid (surplusLocalCoord).
    double basisFnVal = 1.0,basisFnGrad[dim], basisStored[dim];
    for(int idim=0; idim<dim; idim++){
      if(multiIndex[subGrid][idim] == 1){
        if(coord[idim] == 1.0) surplusLocalCoord[idim] = 1;
        else{
          int xp = static_cast<int>(floor(coord[idim]*2));
          if(xp == 0) { 
            basisStored[idim] = 2.0*(0.5 - coord[idim]);//basisFnVal *= 2.0*(0.5 - coord[idim]);
            basisFnGrad[idim] = -2.0;
          }
          else {
            basisStored[idim] = 2.0*(coord[idim] - 0.5);
            basisFnGrad[idim] = 2.0;
          }
          surplusLocalCoord[idim] = xp;
	  map[idim] = 1;
	  if (basisStored[idim] == 0.0)
	    break;
        }
      }else if(multiIndex[subGrid][idim] == 0) surplusLocalCoord[idim] = 0;
      else if(coord[idim] == 0.0){
        basisFnVal = 0.0;
        break;
      }else{
        double scale = 1<<multiIndex[subGrid][idim]; //pow(2.0,multiIndex[subGrid][idim]);
        int xp = static_cast<int>(floor(coord[idim] * scale / 2.0));
        basisStored[idim] = (1.0 - scale*fabs(coord[idim]-(2.0*static_cast<double>(xp)+1.0)/scale));
	if (basisStored[idim] == 0.0)
	  break;

	double f = coord[idim]-(2.0*static_cast<double>(xp)+1.0)/scale;
	if (fabs(f) > 1.0e-8)
	  basisFnGrad[idim] = -scale*fabs(f)/(f); 
	else
	  basisFnGrad[idim] = -scale;
	  
        surplusLocalCoord[idim] = xp;
	map[idim] = 1;
      }

      if(basisFnVal == 0.0) break;

    }

    // contributing surplus in the considered subgrid is added to result
    if(basisFnVal > 0.0){
      // position of the contributing surplus in the array of surpluses
      contributingSurplus = firstSurplus + surplusLocalCoord[0];
      for(int idim=1; idim<dim; idim++)
        contributingSurplus += nCumulatedPointsDim[idim-1]*surplusLocalCoord[idim];
      // contribution added
      basisFnVal = 1.0;
      for (int idim=0; idim<dim; ++idim)
	basisFnVal *= basisStored[idim];
      
      for(int iout=0; iout<out; iout++) {
        //output[iout] += basisFnVal*surplus[contributingSurplus][iout];
        for (int idim=0; idim < dim; ++idim) {
          output[iout*dim + idim] = 0.0;
          for (int idim2 = 0; idim2 < dim; ++idim2) {
	    if (map[idim2])
	      output[iout*dim + idim] += basisFnVal/basisStored[idim2]*basisFnGrad[idim2];
	  }
          output[iout*dim+idim] *= surplus[contributingSurplus][iout];
        }
      }
    }
    firstSurplus += nPointsSubGrid;

  }

}

//------------------------------------------------------------------------------

bool SparseGrid::contains(double *coord){

  for(int idim=0; idim<dim; idim++)
    if(coord[idim] < range[idim][0] || coord[idim] > range[idim][1])
      return false;
  return true;
}

//------------------------------------------------------------------------------

void SparseGrid::interpolate(const int numRes, double **coord, double **res){

  double scaledCoord[dim];
  for(int iPts=0; iPts<numRes; iPts++){
    scale(coord[iPts],scaledCoord, 2);
    if(outOfRange(scaledCoord)) closestPointInRange(scaledCoord);
    singleInterpolation(scaledCoord, res[iPts]);
  }

}

void SparseGrid::interpolateGradient(const int numRes, double **coord, double **res){

  double scaledCoord[dim];
  for(int iPts=0; iPts<numRes; iPts++){
    scale(coord[iPts],scaledCoord, 2);
    if(outOfRange(scaledCoord)) closestPointInRange(scaledCoord);
    singleInterpolationGradient(scaledCoord, res[iPts]);
  }

}

//------------------------------------------------------------------------------

bool SparseGrid::outOfRange(double *coord) const{

  for(int idim=0; idim<dim; idim++)
    if(coord[idim]<0 || coord[idim]>1){
      fprintf(stdout, "*** Warning: coordinates(%e) in dimension %d is out of bounds [0,1]\n", coord[idim], idim);
      return true;
    }

  return false;

}

//------------------------------------------------------------------------------

void SparseGrid::closestPointInRange(double *coord){

  for(int idim=0; idim<dim; idim++){
    if(coord[idim]<0) coord[idim] = 0.0;
    if(coord[idim]>1) coord[idim] = 1.0;
  }

}

//------------------------------------------------------------------------------
// SparseGrid::print and read sparse grid in a file
//------------------------------------------------------------------------------

void SparseGrid::scaleGrid(const double *refIn, const double *refOut, int outputRangeFlag){

  if(false)
  if(refOut && out>1){
  	fprintf(stdout, "*** Warning: it is not recommended to create SparseGrids with more than 1 output and with output scaling factors\n");
  	fprintf(stdout, "***        : if each dimension represents a different physical quantity,\n");
  	fprintf(stdout, "***        : it is not obvious how to treat each one w.r.t. the others\n");
  	fprintf(stdout, "***        : in particular, what does the scalar absolute accuracy mean?\n");
  }
  
  if(refIn && range)
    for(int idim=0; idim<dim; idim++){
      if(refIn[idim] == 0.0){
        fprintf(stdout, "*** Error: SparseGrid is being rescaled with factor 0!\n");
        exit(1);
      }
      range[idim][0] /= refIn[idim]; // min of range
      range[idim][1] /= refIn[idim]; // max of range
    }
    
  for(int idim=0; idim<dim; idim++){
  	if(logMap[idim]){
  	  fprintf(stdout, "logMap for dim %d\n", idim);
      if(range[idim][0] <= 0.0){
        fprintf(stdout, "*** Error: SparseGrid is being remapped (log)\n");
        exit(1);
      }
      range[idim][0] = logMap[idim]->map(range[idim][0]);
      range[idim][1] = logMap[idim]->map(range[idim][1]);
  	}
  }

  if(refOut){
    if(refOut[0]>0) absAccuracy /= refOut[0];

    if(surplus){
      for(int iPts=0; iPts<nPoints; iPts++){
        for(int iout=0; iout<out; iout++){
          if(refOut[iout] == 0.0){
            fprintf(stdout, "*** Error: SparseGrid is being rescaled with factor 0!\n");
            exit(1);
          }
          surplus[iPts][iout] /= refOut[iout];
        }
      }
    }
    
    if(error){
      for(int subGrid=0; subGrid<nSubGrids; subGrid++)
      	for(int iout=0; iout<out; iout++)
      	  error[subGrid][iout]/= refOut[iout];
    }
    
    if(fnmax)
      for(int iout=0; iout<out; iout++)
        fnmax[iout] /= refOut[iout];
        
    if(fnmin)
      for(int iout=0; iout<out; iout++)
        fnmin[iout] /= refOut[iout];
        
  }

  if(outputRangeFlag==1){
    fprintf(stdout, "# The ranges of the sparse grid in each dimension are\n");
    for(int idim=0; idim<dim; idim++)
      fprintf(stdout, "#      [% e % e]\n", range[idim][0],range[idim][1]);
    fprintf(stdout, "\n");
  }

}

//------------------------------------------------------------------------------

void SparseGrid::printToFile(const double *refIn, const double *refOut, const char *filename) const{

  // first check scaling values are strictly positive!
  double *refIn_ = new double[dim];
  for(int idim=0; idim<dim; idim++)
    refIn_[idim] = 1.0;
    
  if(refIn){
    for(int idim=0; idim<dim; idim++){
      if(refIn[idim]<=0.0){
        fprintf(stderr, "*** Error: Sparse Grid file was not written!\n");
        return;
      }
      refIn_[idim] = refIn[idim];
    }
  }
  
  double *refOut_ = new double[out];
  for(int iout=0; iout<out; iout++)
    refOut_[iout] = 1.0;
  
  if(refOut){
    for(int iout=0; iout<out; iout++){
      if(refOut[iout]<=0.0){
        fprintf(stderr, "*** Error: Sparse Grid file was not written!\n");
        return;
      }
      refOut_[iout] = refOut[iout];
    }
  }

  // print to file the number of subgrids and the corresponding multiIndex array
  //           and the number of points and the surpluses.
  FILE *fpSPARSEGRID = fopen(filename, "w");
  if (!fpSPARSEGRID){
    fprintf(stderr, "*** Error: could not open SparseGrid file\n");
    exit(1);
  }

  // header
  fprintf(fpSPARSEGRID, "# sparse grid: hierarchical subgrids and the surpluses\n");

  // characteristics used to create the sparse grid
  fprintf(fpSPARSEGRID, "%d %d\n", dim, out);
  
  for(int idim=0; idim<dim; idim++){
  	if(logMap[idim])
      fprintf(fpSPARSEGRID, "%e ", logMap[idim]->invMap(1.0));
    else
      fprintf(fpSPARSEGRID, "%e ", 0.0);
  }
  fprintf(fpSPARSEGRID, "\n");
  
  fprintf(fpSPARSEGRID, "%d %d\n", minPoints, maxPoints);
  fprintf(fpSPARSEGRID, "%.15e %.15e %.15e %d\n", refOut_[0]*absAccuracy, relAccuracy, dimAdaptDegree, nAdaptivePoints);
  fprintf(fpSPARSEGRID, "%.15e %.15e\n", actualAbsAcc, actualRelAcc);
  for(int idim=0; idim<dim; idim++){
  	double temp[2];
  	if(logMap[idim]){
      temp[0] = logMap[idim]->invMap(range[idim][0]);
      temp[1] = logMap[idim]->invMap(range[idim][1]);
  	}else{
      temp[0] = range[idim][0];
      temp[1] = range[idim][1];
  	}

    fprintf(fpSPARSEGRID, "%.15e %.15e ", temp[0]*refIn_[idim], temp[1]*refIn_[idim]);
  }
  fprintf(fpSPARSEGRID, "\n");

  // the sparse grid data
  fprintf(fpSPARSEGRID, "%d\n", nSubGrids);
  for(int isubgrid=0; isubgrid<nSubGrids; isubgrid++){
    for(int idim=0; idim<dim; idim++)
      fprintf(fpSPARSEGRID, "%d ", multiIndex[isubgrid][idim]);
    fprintf(fpSPARSEGRID, "\n"); 
  }
  
  fprintf(fpSPARSEGRID, "%d\n", nPoints);
  for(int ipts=0; ipts<nPoints; ipts++){
    for(int iout=0; iout<out; iout++){
      fprintf(fpSPARSEGRID, "%.15e ", surplus[ipts][iout]*refOut_[iout]);
    }
    fprintf(fpSPARSEGRID, "\n"); 
  }
  
  // the data to later increase the number of points
  for(int iout=0; iout<out; iout++)
    fprintf(fpSPARSEGRID, "%.15e ", fnmin[iout]*refOut_[iout]);
  fprintf(fpSPARSEGRID, "\n");
  for(int iout=0; iout<out; iout++)
    fprintf(fpSPARSEGRID, "%.15e ", fnmax[iout]*refOut_[iout]);
  fprintf(fpSPARSEGRID, "\n");
  
  // error
  for(int isubgrid=0; isubgrid<nSubGrids; isubgrid++){
  	for(int iout=0; iout<out; iout++){
      fprintf(fpSPARSEGRID, "%.15e ", error[isubgrid][iout]*refOut_[iout]);
  	}
  	fprintf(fpSPARSEGRID, "\n");
  }
  
  // active - activeError - activeCost
  for(int isubgrid=0; isubgrid<nSubGrids; isubgrid++){
  	fprintf(fpSPARSEGRID, "%d %.15e %.15e\n", active[isubgrid], activeError[isubgrid], activeCost[isubgrid]);
  }
  
  // activeHeapError
  activeHeapError.printToFile(fpSPARSEGRID);
  // activeHeapCost
  activeHeapCost.printToFile(fpSPARSEGRID);
  
  // neighbour
  for(int isubgrid=0; isubgrid<nSubGrids; isubgrid++){
  	for(int idim=0; idim<dim; idim++) // forward neighbour
      fprintf(fpSPARSEGRID, "%d ", neighbour[isubgrid][idim]);
    for(int idim=0; idim<dim; idim++) // backward neighbour
      fprintf(fpSPARSEGRID, "%d ", neighbour[isubgrid][idim+dim]);
  	fprintf(fpSPARSEGRID, "\n");
  }
  
  fflush(fpSPARSEGRID);
  fclose(fpSPARSEGRID);

}

//------------------------------------------------------------------------------

void SparseGrid::readFromFile(const double *refIn, const double *refOut, 
                              const char* filename, int outputRangeFlag){

  char mystring [100];
  FILE *fpSPARSEGRID = fopen(filename, "r");
  if (!fpSPARSEGRID){
    fprintf(stderr, "*** Error: could not open SparseGrid file '%s' to read\n", filename);
    exit(1);
  }

  // header
  char* toto = fgets(mystring, 100, fpSPARSEGRID);

  // characteristics used to create the sparse grid
  int toto2 = fscanf(fpSPARSEGRID, "%d %d\n", &dim, &out);
  
  if(logMap)
    for(int idim=0; idim<dim; idim++)
      delete logMap[idim];
  delete [] logMap;
  	
  logMap = new LogarithmicMapping *[dim];
  for(int idim=0; idim<dim; idim++){
  	double temp;
  	toto2 = fscanf(fpSPARSEGRID, "%lf ", &temp); 
  	if(temp > 0.0)
      logMap[idim] = new LogarithmicMapping(temp);
    else if(temp == 0.0)
      logMap[idim] = 0;
    else{
      fprintf(stdout, "*** Error: negative base for logarithmic mapping\n");
      exit(1);
    }
  }

  toto2 = fscanf(fpSPARSEGRID, "%d %d", &minPoints, &maxPoints);
  toto2 = fscanf(fpSPARSEGRID, "%lf %lf %lf %d", &absAccuracy, &relAccuracy, &dimAdaptDegree, &nAdaptivePoints);
  toto2 = fscanf(fpSPARSEGRID, "%lf %lf", &actualAbsAcc, &actualRelAcc);
  delete [] range;
  range = new Range[dim];
  for(int idim=0; idim<dim; idim++)
    toto2 = fscanf(fpSPARSEGRID, "%lf %lf ", &(range[idim][0]), &(range[idim][1]));
  // WARNING: just like for constructor, logMap is applied later, in scaleGrid(...)

  // the sparse grid data
  toto2 = fscanf(fpSPARSEGRID, "%d", &nSubGrids);
  multiIndex = new int *[nSubGrids];
  for(int isubgrid=0; isubgrid<nSubGrids; isubgrid++){
    multiIndex[isubgrid] = new int[dim];
    for(int idim=0; idim<dim; idim++)
      toto2 = fscanf(fpSPARSEGRID, "%d ", &(multiIndex[isubgrid][idim]));
  }

  toto2 = fscanf(fpSPARSEGRID, "%d", &nPoints);
  surplus = new double *[nPoints];
  for(int ipts=0; ipts<nPoints; ipts++){
    surplus[ipts] = new double[out];
    for(int iout=0; iout<out; iout++)
      toto2 = fscanf(fpSPARSEGRID, "%lf ", &(surplus[ipts][iout]));
  }

  // the data to later increase the number of points
  delete [] fnmin; delete [] fnmax;
  fnmin = new double[out];
  fnmax = new double[out];
  for(int iout=0; iout<out; iout++)
    toto2 = fscanf(fpSPARSEGRID, "%lf ", &(fnmin[iout]));
  for(int iout=0; iout<out; iout++)
    toto2 = fscanf(fpSPARSEGRID, "%lf ", &(fnmax[iout]));
    
  // error
  error = new double *[nSubGrids];
  for(int isubgrid=0; isubgrid<nSubGrids; isubgrid++){
  	error[isubgrid] = new double[out];
  	for(int iout=0; iout<out; iout++){
      toto2 = fscanf(fpSPARSEGRID, "%lf ", &(error[isubgrid][iout]));
  	}
  }
  
  // active - activeError - activeCost
  active      = new bool[nSubGrids];
  activeError = new double[nSubGrids];
  activeCost  = new double[nSubGrids];
  int intToBool;
  for(int isubgrid=0; isubgrid<nSubGrids; isubgrid++){
  	toto2 = fscanf(fpSPARSEGRID, "%d %lf %lf", &intToBool, &(activeError[isubgrid]), &(activeCost[isubgrid]));
  	active[isubgrid] = intToBool;
  }
  
  // activeHeapError
  activeHeapError.readFromFile(fpSPARSEGRID);
  // activeHeapCost
  activeHeapCost.readFromFile(fpSPARSEGRID);
  
  // neighbour
  neighbour = new int *[nSubGrids];
  for(int isubgrid=0; isubgrid<nSubGrids; isubgrid++){
  	neighbour[isubgrid] = new int[2*dim];
  	for(int idim=0; idim<dim; idim++) // forward neighbour
      toto2 = fscanf(fpSPARSEGRID, "%d ", &(neighbour[isubgrid][idim]));
    for(int idim=0; idim<dim; idim++) // backward neighbour
      toto2 = fscanf(fpSPARSEGRID, "%d ", &(neighbour[isubgrid][idim+dim]));
  }
  
  fclose(fpSPARSEGRID);

  //fill in the rest
  sizeSurplus = nPoints;
  sizeMultiIndex = nSubGrids;
  
  scaleGrid(refIn, refOut, outputRangeFlag);

}

//------------------------------------------------------------------------------
// Auxiliary functions
//------------------------------------------------------------------------------

void SparseGrid::resizeSurplus(){

  double **larger = new double *[2*sizeSurplus];
  for(int i=0; i<2*sizeSurplus; i++)
  	larger[i] = new double[out];
  for(int i=0; i<sizeSurplus; i++)
    for(int iout=0; iout<out; iout++)
      larger[i][iout] = surplus[i][iout];

  for(int i=0; i<sizeSurplus; i++) delete [] surplus[i];
  delete [] surplus;

  sizeSurplus *= 2;
  surplus = larger;

}

//------------------------------------------------------------------------------

void SparseGrid::resizeMultiIndex(){

  int **largerMultiIndex = new int*[2*sizeMultiIndex];
  for(int i=0; i<2*sizeMultiIndex; i++)
    largerMultiIndex[i] = new int[dim];
  for(int i=0; i<sizeMultiIndex; i++)
    for(int idim=0; idim<dim; idim++)
      largerMultiIndex[i][idim] = multiIndex[i][idim];

  for(int i=0; i<sizeMultiIndex; i++) delete [] multiIndex[i];
  delete [] multiIndex;
  multiIndex = largerMultiIndex;

  double *largerErrorActiveSet = new double[2*sizeMultiIndex];
  double *largerCostActiveSet = new double[2*sizeMultiIndex];
  bool *largerActive= new bool[2*sizeMultiIndex];
  for(int i=0; i<sizeMultiIndex; i++){
    largerErrorActiveSet[i]= activeError[i];
    largerCostActiveSet[i]= activeCost[i];
    largerActive[i]= active[i];
  }

  delete [] activeError; delete [] activeCost; delete [] active;
  activeError = largerErrorActiveSet;
  activeCost = largerCostActiveSet;
  active= largerActive;

  double **largerError = new double *[2*sizeMultiIndex];
  for(int i=0; i<2*sizeMultiIndex; i++)
    largerError[i] = new double[out];
  for(int i=0; i<sizeMultiIndex; i++)
    for(int iout=0; iout<out; iout++)
      largerError[i][iout] = error[i][iout];

  for(int i=0; i<sizeMultiIndex; i++) delete [] error[i];
  delete [] error;
  error = largerError;

  int **largerNeighbour = new int *[2*sizeMultiIndex];
  for(int i=0; i<2*sizeMultiIndex; i++)
    largerNeighbour[i] = new int[2*dim];
  for(int i=0; i<2*sizeMultiIndex; i++)
    for(int ineigh=0; ineigh<2*dim; ineigh++)
      largerNeighbour[i][ineigh] = -1;
  for(int i=0; i<nSubGrids; i++)
    for(int ineigh=0; ineigh<2*dim; ineigh++)
      largerNeighbour[i][ineigh] = neighbour[i][ineigh];

  for(int i=0; i<sizeMultiIndex; i++) delete [] neighbour[i];
  delete [] neighbour;
  neighbour = largerNeighbour;

  sizeMultiIndex *= 2;

}

//------------------------------------------------------------------------------
// SparseGrid::Heap implementation
//------------------------------------------------------------------------------

SparseGrid::Heap::Heap(){

  size_    = 0;
  numElem_ = 0;
  elem_    = 0;

}

//------------------------------------------------------------------------------

SparseGrid::Heap::~Heap(){

  delete [] elem_;

}

//------------------------------------------------------------------------------

SparseGrid::Heap::Heap(const Heap &heap){

  size_    = heap.size_;
  numElem_ = heap.numElem_;
  elem_    = new int[size_];
  for(int i=0; i<size_; i++)
    elem_[i] = heap.elem_[i];

}

//------------------------------------------------------------------------------

void SparseGrid::Heap::sort(int index, const double *value){

  if(index>numElem_-1 || index<0) return;

  int parent, temp;
  while(index>0){
    parent = static_cast<int>(floor((index-1)/2.0));
    if(value[elem_[index]]>value[elem_[parent]]){
      temp = elem_[parent];
      elem_[parent] = elem_[index];
      elem_[index] = temp;
    }else 
      break;
    index = parent;
  }

}

//------------------------------------------------------------------------------

void SparseGrid::Heap::insert(const int newElem, const double *value){

  if(numElem_+1>size_){//resize
    if(size_>0) size_ *= 2;
    else size_++;
    int *temp = new int[size_];
    for(int i=0; i<numElem_; i++) temp[i] = elem_[i];
    delete [] elem_;
    elem_ = temp;
  }

  elem_[numElem_] = newElem;
  numElem_++;
  sort(numElem_-1,value);

}

//------------------------------------------------------------------------------

int SparseGrid::Heap::pop(const double *value){

  if(numElem_==0) return -1;

  int res = elem_[0];

  // reorder the heap 
  // first, fillup empty spot by children
  int i = 0;
  int firstChild, secondChild;
  while(2*i+1 < numElem_-1){ // condition for i to have to 2 children
    firstChild = 2*i+1;
    secondChild = firstChild + 1;
    if(value[elem_[firstChild]]>value[elem_[secondChild]]){
      elem_[i] = elem_[firstChild];
      i = firstChild;
    }else{
      elem_[i] = elem_[secondChild];
      i = secondChild;
    }
  }
  // second, put last element in last empty spot
  elem_[i] = elem_[numElem_-1];
  // third, make sure elem at position i is well positioned
  sort(i,value);

  numElem_--;

  return res;


}

//------------------------------------------------------------------------------

void SparseGrid::Heap::printToFile(FILE *file) const{

  fprintf(file, "%d\n", numElem_);
  for(int ielem=0; ielem<numElem_; ielem++){
  	fprintf(file, "%d\n", elem_[ielem]);
  }
}

//------------------------------------------------------------------------------

void SparseGrid::Heap::readFromFile(FILE *file){
	
  int toto = fscanf(file, "%d", &numElem_);
  size_ = numElem_;
  elem_ = new int[numElem_];
  for(int ielem=0; ielem<numElem_; ielem++){
    toto = fscanf(file, "%d", &(elem_[ielem]));
  }
}

//------------------------------------------------------------------------------
// functions for DEBUG PURPOSES
//------------------------------------------------------------------------------

void SparseGrid::messages(const int flag, const int arg) const{


  if(flag==10){
  	fprintf(stdout, "##############################################################\n");
  	fprintf(stdout, "# tabulation done\n");
  	fprintf(stdout, "# number of subgrids = %d\n", nSubGrids);
  	fprintf(stdout, "# number of points   = %d\n", nPoints);
  	fprintf(stdout, "# max number of points   = %d\n", maxPoints);
  	fprintf(stdout, "# relative and absolute accuracy = %e %e\n", actualRelAcc, actualAbsAcc);
  	fprintf(stdout, "##############################################################\n");
  }
  
  if(flag==1 || verbose>0){

  if(flag==1){
    fprintf(stdout, "###########################################\n");
    fprintf(stdout, "#     INITIALIZATION\n");
    fprintf(stdout, "###########################################\n");
  }

  if(flag==2 || flag==10){
    fprintf(stdout, "#     BEGIN ITERATION\n");
    fprintf(stdout, "###########################################\n");
    if(verbose>1){
      fprintf(stdout, "#     STATE AT BEGINNING OF ITERATION\n");
      fprintf(stdout, "# activity of grids is:");
      for(int kk = 0; kk<nSubGrids; kk++)
        fprintf(stdout, "  %d", active[kk]);
      fprintf(stdout, "\n");
      fprintf(stdout, "# multiIndex is:");
      for(int idim = 0; idim<dim; idim++){
        for(int kk = 0; kk<nSubGrids; kk++)
          fprintf(stdout, "  %d", multiIndex[kk][idim]);
        fprintf(stdout, "\n#               ");
      }
      fprintf(stdout, "\n# activeError is:");
      for(int kk = 0; kk<nSubGrids; kk++)
        fprintf(stdout, "  %e", activeError[kk]);
      fprintf(stdout, "\n");
      fprintf(stdout, "# activeCost is:");
      for(int kk = 0; kk<nSubGrids; kk++)
        fprintf(stdout, "  %e", activeCost[kk]);
      fprintf(stdout, "\n");
      fprintf(stdout, "# errors are:");
      for(int iout = 0; iout<out; iout++){
        for(int kk = 0; kk<nSubGrids; kk++)
          fprintf(stdout, "  %e", error[kk][iout]);
        fprintf(stdout, "\n#               ");
      }
      fprintf(stdout, "\n");
      fprintf(stdout, "# surpluses are:");
      for(int iout = 0; iout<out; iout++){
        for(int kk = 0; kk<nPoints; kk++)
          fprintf(stdout, "  %e", surplus[kk][iout]);
        fprintf(stdout, "\n#               ");
      }
      fprintf(stdout, "\n");
      fprintf(stdout, "# active grids (error)");
      for(int i=0; i<activeHeapError.size(); i++)
        fprintf(stdout, "  %d", activeHeapError[i]);
      fprintf(stdout, "\n");
      fprintf(stdout, "# active grids (cost)");
      for(int i=0; i<activeHeapCost.size(); i++)
        fprintf(stdout, "  %d", activeHeapCost[i]);
      fprintf(stdout, "\n###########################################\n");
    }
  }

  if(flag==3 && verbose>3){
    fprintf(stdout, "# adding grid %d to inactive grid set\n", arg);
    fprintf(stdout, "###########################################\n");
  }
  
  if(flag==4 && verbose>3){
    fprintf(stdout, "# neighbour of %d is admissible\n", arg);
    fprintf(stdout, "###########################################\n");
  }

  if(flag==5 && verbose>3){
    fprintf(stdout, "# neighbours of %d are:", arg);
    for(int ineigh=0; ineigh<2*dim; ineigh++)
      fprintf(stdout, "  %d", neighbour[arg][ineigh]);
    fprintf(stdout, "\n");
    fprintf(stdout, "###########################################\n");
  }

  if(flag==6 && verbose>3){
    fprintf(stdout, "# forward admissible grid %d leads to %d points and %d subgrids\n", arg, nPoints, nSubGrids);
    fprintf(stdout, "###########################################\n");
  }

  if(flag==7){
    fprintf(stdout, "#    desired absolute accuracy (%e < %e) has been reached\n", actualAbsAcc, absAccuracy);
    fprintf(stdout, "# or desired relative accuracy (%e < %e) has been reached\n", actualRelAcc, relAccuracy);
    fprintf(stdout, "###########################################\n");
  }

  if(flag==8 && verbose>2){
    fprintf(stdout, "### singleInterpolation -- nSubGrids = %d\n", nSubGrids);
    for(int subGrid=0; subGrid<nSubGrids; subGrid++){
      for(int debugDim=0; debugDim<dim; debugDim++)
        if(verbose>4) fprintf(stdout, "###      for multiIndex[%d][%d] = %d\n", subGrid, debugDim, multiIndex[subGrid][debugDim]);
    }
  }

  if(flag==9){
    fprintf(stdout, "###########################################\n");
    fprintf(stdout, "# After reading sparse grid file,");
    fprintf(stdout, "# min/maxPoints = %d / %d\n", minPoints, maxPoints);
    fprintf(stdout, "# abs/relAccuracy = %e / %e\n", absAccuracy, relAccuracy);
    fprintf(stdout, "# degree of dimensional adaptivity = %e\n", dimAdaptDegree);
    for(int idim=0; idim<dim; idim++)
      fprintf(stdout, "# range[%d] is [%e %e]\n", idim, range[idim][0], range[idim][1]);
    fprintf(stdout, "###########################################\n");
  }
  }
  
}

//------------------------------------------------------------------------------
