#pragma once

#include <DistMatrix.h>

#include <MultiGridSmoothingMatrix.h>
#include <MultiGridKernel.h>
#include <MultiGridMvpMatrix.h>
#include <MultiGridDistSVec.h>

template <class Scalar,int neq>
class MultiGridMatVecProd {

  DistMat<Scalar,neq>* macroA;
                                                                                                      
  MultiGridLevel<Scalar>* multiGridLevel;
  
 public:
  
  MultiGridMatVecProd(DistMat<Scalar,neq>* macroA, MultiGridLevel<Scalar>* multiGridLevel) :      
                      macroA(macroA), multiGridLevel(multiGridLevel) { }                          
  
  ~MultiGridMatVecProd() { }                                                                      

  void setA(DistMat<Scalar,neq>* A) { macroA = A; } 
 
  void apply(DistSVec<Scalar,neq>& x, DistSVec<Scalar,neq>& b) {                                      
  
    b = 0.0;                                                                                          
  
    multiGridLevel->computeMatVecProd(*dynamic_cast< DistMvpMatrix<Scalar,neq>*>(macroA),             
                                      x, b);                                                          
                                                                                                      
  }

  void applyTranspose(DistSVec<Scalar,neq>& x, DistSVec<Scalar,neq>& b) {}
  void applyT(DistSVec<Scalar,neq>& x, DistSVec<Scalar,neq>& b) {}

};

template<class Scalar,int neq,class Scalar2>                                                          
class MultiGridJacobiPrec {
  
  MultiGridSmoothingMatrix<Scalar,neq>** smoothingMatrices;                                           
  MultiGridLevel<Scalar>* level;                                                                      
                                                                                                      
 public:
                                                                                                      
  MultiGridJacobiPrec(MultiGridSmoothingMatrix<Scalar,neq>** smoothingMatrices,                   
                      MultiGridLevel<Scalar>* level) :                                            
                      smoothingMatrices(smoothingMatrices),                                        
                      level(level) { }                                                             
  
  ~MultiGridJacobiPrec() { }                                                                      
  
  void apply(DistSVec<Scalar2,neq>& x, DistSVec<Scalar2,neq>& b) {

    b = 0.0;
    int numLocSub = x.info().numLocSub;
#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub) {
      smoothingMatrices[iSub]->smooth(x(iSub),b(iSub));
    }
    level->assemble(b);
  }
};

template<class Scalar,int neq,class Scalar2 = double>
class MultiGridRASPrec {

  MultiGridSmoothingMatrix<Scalar,neq>** smoothingMatrices;
  MultiGridLevel<Scalar>* level;

 public:

  MultiGridRASPrec(MultiGridSmoothingMatrix<Scalar,neq>** smoothingMatrices,
                   MultiGridLevel<Scalar>* level) :
                   smoothingMatrices(smoothingMatrices),
                   level(level) { }

  ~MultiGridRASPrec() { }

  void apply(DistSVec<Scalar2,neq>& x, DistSVec<Scalar2,neq>& b) {

    b = 0.0;
    int numLocSub = x.info().numLocSub;
#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub) {
      smoothingMatrices[iSub]->smooth(x(iSub),b(iSub));
    }
    level->assemble(b);
  }

  void applyTranspose(DistSVec<Scalar2,neq>& x, DistSVec<Scalar2,neq>& b) {}
  void applyT(DistSVec<Scalar2,neq>& x, DistSVec<Scalar2,neq>& b) {}
};

template<class Scalar,int neq,class Scalar2>
class MultiGridKspSolver {

 public:

  MultiGridKspSolver(Domain*,KspData&,MultiGridKernel<Scalar>* pKernel);

  ~MultiGridKspSolver();

  void solve(int lvl, MultiGridMvpMatrix<Scalar,neq>&,
             MultiGridDistSVec<Scalar,neq>& f,
             MultiGridDistSVec<Scalar,neq>& x);

 private:

  KspSolver<DistSVec<Scalar2,neq>, 
            MultiGridMatVecProd<Scalar2,neq> ,
            MultiGridRASPrec<Scalar2,neq>, Communicator>** coarseSolvers;                         
    
  MultiGridMatVecProd<Scalar2,neq>** coarseMvps;
  MultiGridRASPrec<Scalar2,neq>** coarsePrecs;

  int nLevels;

  int numLocSub;

  MultiGridKernel<Scalar>* pKernel;

  MultiGridSmoothingMatrix<Scalar,neq>*** smoothingMatrices;
};

