#include <Communicator.h>
#include <cmath>
#include <VecSetOp.h>

//------------------------------------------------------------------------------

template<int dim>
ImplicitGnatTsDesc<dim>::ImplicitGnatTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitGappyTsDesc<dim>(ioData, geoSource, dom),
	leastSquaresSolver(this->com, this->com->size(), 1)// all cpus along rows
{

  nPodJac = 0;  // set when reading the online matrices
  numResJacMat = this->rom->getNumResJacMat();

  leastSquaresSolver.blockSizeIs(32);
  
}

//------------------------------------------------------------------------------

template<int dim>
ImplicitGnatTsDesc<dim>::~ImplicitGnatTsDesc() {}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGnatTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop, DistSVec<double, dim> &U, const int& totalTimeSteps)  {

  double t0 = this->timer->getTime();
  // homotopy on reduced-coordinates for spatial-only problems
  double invSqrtHomotopyStep = 0.0;
  Vec<double> dUrom(this->dUromTimeIt);
  if (this->spatialOnlyWithHomotopy) {
    double homotopyStep = min(this->homotopyStepInitial*pow(this->homotopyStepGrowthRate,totalTimeSteps),this->homotopyStepMax);
    if (it==0) this->com->fprintf(stdout, " ... homotopy step %e\n", homotopyStep);
    invSqrtHomotopyStep = pow(homotopyStep,-0.5);
    dUrom *= invSqrtHomotopyStep;
  }

  // Form A * of and distribute
  transMatMatProd(*(this->rom->getJacMat()), *this->AJRestrict,this->jactmp);

  for (int iCol = 0; iCol < leastSquaresSolver.localCols(); ++iCol) {
    const int globalColIdx = leastSquaresSolver.globalColIdx(iCol);
    const int colOffset = globalColIdx * leastSquaresSolver.equationCount();
    for (int iRow = 0; iRow < leastSquaresSolver.localRows(); ++iRow) {
      const int globalRowIdx = leastSquaresSolver.globalRowIdx(iRow);
      int homotopyRowIdx = globalRowIdx - nPodJac;
      if (this->spatialOnlyWithHomotopy && (homotopyRowIdx >= 0)) {
        leastSquaresSolver.matrixEntry(iRow, iCol) = (globalColIdx==homotopyRowIdx) ? invSqrtHomotopyStep : 0.0;
      } else {
        leastSquaresSolver.matrixEntry(iRow, iCol) = this->jactmp[globalRowIdx + colOffset];
      }
    }
  }

  // Form B * ResRestrict and distribute
  // NOTE: do not check if current CPU has rhs

  transMatVecProd(*(this->rom->getResMat()), *this->ResRestrict, this->column);

  for (int iRow = 0; iRow < leastSquaresSolver.localRows(); ++iRow) {
    const int globalRowIdx = leastSquaresSolver.globalRowIdx(iRow);
    if (this->spatialOnlyWithHomotopy && (globalRowIdx >= nPodJac)) {
      leastSquaresSolver.rhsEntry(iRow) = -dUrom[globalRowIdx-nPodJac];
    } else {
      leastSquaresSolver.rhsEntry(iRow) = -this->column[globalRowIdx];
    }
  }
  this->timer->addLinearSystemFormTime(t0);

  // Check convergence 
  t0 = this->timer->getTime();
  double romResidual[this->nPod];

  for (int iCol = 0; iCol < this->nPod; ++iCol) {
    romResidual[iCol] = 0.0;
  }

  for (int iCol = 0; iCol < leastSquaresSolver.localCols(); ++iCol) {
    const int globalColIdx = leastSquaresSolver.globalColIdx(iCol);
    for (int iRow = 0; iRow < leastSquaresSolver.localRows(); ++iRow) {
      romResidual[globalColIdx] += leastSquaresSolver.matrixEntry(iRow, iCol) * leastSquaresSolver.rhsEntry(iRow);   
    }
  }

  this->com->globalSum(this->nPod,romResidual);
  res = 0.0;
  for (int iCol = 0; iCol < this->nPod; ++iCol) {
    res += romResidual[iCol] * romResidual[iCol];
  }

  if (it == 0) {
    this->res0 = res;
    this->target = this->epsNewton * this->res0;
    this->com->fprintf(stderr, " ... ROM residual = %e (target = %e)\n", res, this->target);
  } else {
    this->com->fprintf(stderr, " ... ROM residual = %e\n", res);
  }
  if (breakloop = (res == 0.0) || (res <= this->target)) return;
  this->timer->addCheckConvergenceTime(t0);


  // Solve least squares problem
  t0 = this->timer->getTime();
  leastSquaresSolver.solve();

  // Update vector: The first nPod rows give the components in the pod basis
  t0 = this->timer->getTime();
  this->dUromNewtonIt = 0.0;
  for (int localIRow = 0; localIRow < leastSquaresSolver.localSolutionRows(); ++localIRow) {
    const int iRow = leastSquaresSolver.globalRhsRowIdx(localIRow);
    this->dUromNewtonIt[iRow] = leastSquaresSolver.rhsEntry(localIRow);
  }
 
  // Consolidate across the cpus
  this->com->globalSum(this->nPod, this->dUromNewtonIt.data());
  this->timer->addLinearSystemSolveTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitGnatTsDesc<dim>::setProblemSize(DistSVec<double, dim> &U) {
 
  nPodJac = this->rom->getJacMat()->numVectors();

  if (this->spatialOnlyWithHomotopy) {
    leastSquaresSolver.problemSizeIs(nPodJac + this->nPod, this->nPod);
  } else {
    leastSquaresSolver.problemSizeIs(nPodJac, this->nPod);
  }

  if (this->jactmp) delete [] this->jactmp;
  if (this->column) delete [] this->column;
  this->jactmp = new double [nPodJac * this->nPod];
  this->column = new double [nPodJac];

  this->AJRestrict.reset(new VecSet<DistSVec<double, dim> >(this->nPod, this->rom->getRestrictedDistInfo()));
  this->ResRestrict.reset(new DistSVec<double, dim> (this->rom->getRestrictedDistInfo())); 
}

//------------------------------------------------------------------------------

template<int dim>
double ImplicitGnatTsDesc<dim>::meritFunction(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> &dQ, DistSVec<double, dim> &F, double stepLength)  {
  // merit function: norm of the residual (want to minimize residual)

  DistSVec<double, dim> newQ(this->domain->getNodeDistInfo());
  newQ = Q + stepLength*dQ;
  this->checkSolution(newQ);
  this->computeFullResidual(it,newQ,true,&F);

  double merit = 0.0;

  if (this->ioData->romOnline.meritFunction == NonlinearRomOnlineData::HDM_RESIDUAL) {
    // merit function = 1/2 * (norm of full-order residual)^2
    transMatVecProd(*(this->rom->getResMat()), *this->ResRestrict, this->column);
    merit = 0.0;
    for (int iRow=0; iRow<nPodJac; ++iRow) {
      merit += this->column[iRow] * this->column[iRow];
    }
    merit *= 0.5;
  } else if (this->ioData->romOnline.meritFunction == NonlinearRomOnlineData::ROM_RESIDUAL) {
    // Form the ROM residual
    Vec<double> romResidual;
    romResidual.resize(this->nPod);
    this->computeAJ(it, newQ, true);
    transMatMatProd(*(this->rom->getJacMat()), *this->AJRestrict,this->jactmp);
    transMatVecProd(*(this->rom->getResMat()), *this->ResRestrict, this->column);
  
    // Check convergence 
    romResidual = 0.0;
    for (int iCol = 0; iCol < this->nPod; ++iCol) {
      for (int iRow = 0; iRow < nPodJac; ++iRow) {
        romResidual[iCol] += this->jactmp[iRow + nPodJac*iCol] * this->column[iRow];   
      }
    }
    merit = romResidual*romResidual;
    merit *= 0.5;
  } else {
    fprintf(stderr,"*** Error: unrecognized choice of merit function!\n");
    exit(-1);
  }


  return merit;

}


