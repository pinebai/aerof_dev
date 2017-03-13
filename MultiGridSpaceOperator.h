#pragma once

#include <MultiGridKernel.h>

#include <MultiGridMvpMatrix.h>

#include <MultiGridDistSVec.h>

#include <LevelSet/MultiGridLevelSetStructure.h>

template <class Scalar,int dim>
class MultiGridSpaceOperator {

 public:

  MultiGridSpaceOperator(IoData&,Domain*, SpaceOperator<dim>*,
                         MultiGridKernel<Scalar>* pKernel,
                         SpaceOperator<dim>* = NULL,
                         SpaceOperator<dim>* = NULL);

  ~MultiGridSpaceOperator();

  void setupBcs(DistBcData<dim>* bcd);

  void computeTimeStep(int level, double cfl, 
                       MultiGridDistSVec<Scalar,dim>& V);

  void computeResidual(int level, MultiGridDistSVec<Scalar,dim>& U,
                       MultiGridDistSVec<Scalar,dim>& V,
                       MultiGridDistSVec<Scalar,dim>& res,
                       bool addDWdt = true);

  void computeResidualEmbedded(DistExactRiemannSolver<dim>&,
			       int level, MultiGridDistSVec<Scalar,dim>& U,
			       MultiGridDistSVec<Scalar,dim>& V,
			       MultiGridDistSVec<Scalar,dim>& res,
			       DistMultiGridLevelSetStructure*,
			       bool addDWdt = true);

  void add_dAW_dtEmbedded(int level,MultiGridDistSVec<Scalar,dim>& U,
			  MultiGridDistSVec<Scalar,dim>& res,
			  DistMultiGridLevelSetStructure*);

  void updateStateVectors(int lvl, MultiGridDistSVec<Scalar,dim>& U) ;

  template <int neq>
  void computeJacobian(int level, MultiGridDistSVec<Scalar,dim>& U,
                       MultiGridDistSVec<Scalar,dim>& V,
                       MultiGridMvpMatrix<Scalar,neq>& mvp);

  template <int neq>
    void computeTurbulentJacobian(int level, MultiGridDistSVec<Scalar,dim>& U,
				  MultiGridDistSVec<Scalar,dim>& V,
				  MultiGridMvpMatrix<Scalar,neq>& mvp);

  template <int neq>
    void computeJacobianEmbedded(DistExactRiemannSolver<dim>&,
				 int level, MultiGridDistSVec<Scalar,dim>& U,
				 MultiGridDistSVec<Scalar,dim>& V,
				 MultiGridMvpMatrix<Scalar,neq>& mvp,
				 DistMultiGridLevelSetStructure*);


  MultiGridOperator<Scalar,dim>* getOperator(int i) {

    return myOperators[i];
  }

 private:

  int nLevels;

  MultiGridKernel<Scalar>* pKernel;

  MultiGridOperator<Scalar,dim>** myOperators;

  RecFcnConstant<dim> recConstant;

  FemEquationTerm* fet,*fet1, *fet2; 

  FluxFcn** fluxFcn, **fluxFcn1, **fluxFcn2;

  VarFcn* varFcn;
};
