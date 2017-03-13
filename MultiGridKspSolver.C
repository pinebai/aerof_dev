#include <MultiGridKspSolver.h>

template<class Scalar,int neq,class Scalar2>
MultiGridKspSolver<Scalar,neq,Scalar2>::
MultiGridKspSolver(Domain* domain,KspData& coarseSolverData,
                   MultiGridKernel<Scalar>* pKernel)  : pKernel(pKernel) {

  nLevels = pKernel->numLevels();
  
  numLocSub = domain->getNodeDistInfo().numLocSub;

  coarseMvps = new MultiGridMatVecProd<Scalar2,neq>*[nLevels];
  coarsePrecs = new MultiGridRASPrec<Scalar2,neq>*[nLevels];
  coarseSolvers = new  KspSolver<DistSVec<Scalar2,neq>, MultiGridMatVecProd<Scalar2,neq>,
                                   MultiGridRASPrec<Scalar2,neq>, Communicator>*[nLevels];

  typename MultiGridSmoothingMatrix<Scalar2,neq>::SmoothingMode smoothing_mode;
/*  if (ioData.mg.mg_smoother == PcData::MGJACOBI)
    smoothing_mode = MultiGridSmoothingMatrix<Scalar2,neq>::BlockJacobi;
  else if (ioData.mg.mg_smoother == PcData::MGLINEJACOBI)
    smoothing_mode = MultiGridSmoothingMatrix<Scalar2,neq>::LineJacobi;
  else //if (pcData.mg_smoother == PcData::MGRAS) 
*/
  smoothing_mode = MultiGridSmoothingMatrix<Scalar2,neq>::RAS;
    
  smoothingMatrices = new MultiGridSmoothingMatrix<Scalar2,neq>**[nLevels];
  for(int level = 1; level < nLevels; ++level) {
    smoothingMatrices[level] = new MultiGridSmoothingMatrix<Scalar2,neq>*[numLocSub];
#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub)
      smoothingMatrices[level][iSub] = new
        MultiGridSmoothingMatrix<Scalar2,neq>(MultiGridSmoothingMatrix<Scalar2,neq>::RAS,
                                              iSub,
                                              pKernel->getLevel(level)->getNodeDistInfo().subSize(iSub),
                                              pKernel->getLevel(level)->getEdgeDistInfo().subSize(iSub),
                                              0,
                                              ((level < nLevels-1) ? pKernel->getLevel(level+1): NULL),
                                              pKernel->getLevel(level));
  } 


  for (int lvl = 1; lvl < nLevels; ++lvl) {
    coarseMvps[lvl] = new MultiGridMatVecProd<Scalar2,neq>(NULL,
                                                           pKernel->getLevel(lvl));

    coarsePrecs[lvl] = new MultiGridRASPrec<Scalar2,neq>(smoothingMatrices[lvl],
                                                         pKernel->getLevel(lvl));

    coarseSolvers[lvl] = new GmresSolver<DistSVec<Scalar2,neq>, MultiGridMatVecProd<Scalar2,neq>,
                                   MultiGridRASPrec<Scalar2,neq>, Communicator>(
                       pKernel->getLevel(lvl)->getNodeDistInfo(), coarseSolverData, coarseMvps[lvl],
                       coarsePrecs[lvl], domain->getCommunicator());

  }

}

template<class Scalar,int neq,class Scalar2>
MultiGridKspSolver<Scalar,neq,Scalar2>::
~MultiGridKspSolver() {
  
  for (int i = 1; i < nLevels; ++i) {

    delete coarseMvps[i];
    delete coarsePrecs[i];
    delete coarseSolvers[i];

#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub) {

      delete smoothingMatrices[i][iSub];
    }
  
    delete [] smoothingMatrices[i];

  }

  delete [] smoothingMatrices;

  delete [] coarseMvps;
  delete [] coarsePrecs;
  delete [] coarseSolvers;

}

template<class Scalar,int neq,class Scalar2>
void MultiGridKspSolver<Scalar,neq,Scalar2>::
solve(int lvl, MultiGridMvpMatrix<Scalar,neq>& M,
      MultiGridDistSVec<Scalar,neq>& f,
      MultiGridDistSVec<Scalar,neq>& x) {

  coarseMvps[lvl]->setA(&M(lvl));
  
#pragma omp parallel for1
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    smoothingMatrices[lvl][iSub]->getData(M(lvl)(iSub));
  }

  coarseSolvers[lvl]->setup(0,1,f(lvl));
  
  coarseSolvers[lvl]->solve(f(lvl),x(lvl));
}

template class MultiGridKspSolver<double,1,double>;
template class MultiGridKspSolver<double,2,double>;
template class MultiGridKspSolver<double,3,double>;

template class MultiGridKspSolver<double,5,double>;
template class MultiGridKspSolver<double,6,double>;
template class MultiGridKspSolver<double,7,double>;
