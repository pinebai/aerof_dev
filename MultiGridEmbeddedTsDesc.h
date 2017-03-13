#pragma once

#include <ImplicitEmbeddedCoupledTsDesc.h>

#include <MultiGridKernel.h>

#include <MultiGridDistSVec.h>

#include <MultiGridSpaceOperator.h>

#include <MultiGridSmoothingMatrices.h>

#include <MultiGridKspSolver.h>

#include <LevelSet/MultiGridLevelSetStructure.h>

template <int dim>
class MultiGridEmbeddedTsDesc : public ImplicitEmbeddedCoupledTsDesc<dim> {

 public:

  MultiGridEmbeddedTsDesc(IoData &, GeoSource &, Domain *); 

  ~MultiGridEmbeddedTsDesc();

  void smooth0(DistSVec<double,dim>& x,int steps);

  void smooth(int lvl, MultiGridDistSVec<double,dim>& x,
	      DistSVec<double,dim>& f,int steps,
	      bool postsmooth);

  void cycle(int lvl, DistSVec<double,dim>& f,
             MultiGridDistSVec<double,dim>& x);

  void cycle(DistSVec<double,dim>& x);
  
  void setupTimeStepping(DistSVec<double,dim> *U, IoData &iod);

  DistSVec<double,dim>* getSmoothedVec() { return U_smoothed; }

 private:

  int numSmooths_pre[10],numSmooths_post[10];

  bool smoothWithGMRES;

  int mc;

  double prolong_relax_factor,restrict_relax_factor;

  MultiGridKernel<double>* pKernel;

  MultiGridDistSVec<double,dim> V, res, R, F, U,Uold,dx,Forig;

  MultiGridSpaceOperator<double,dim>* mgSpaceOp;

  MultiGridSmoothingMatrices<double,dim>* smoothingMatrices;

  MultiGridMvpMatrix<double,dim>* mgMvp;

  MultiGridKspSolver<double,dim,double>* mgKspSolver;

  DistMultiGridLevelSetStructure** mgLSS;

  int globalIt;

  DistSVec<double,dim>* U_smoothed;

  double densityMin, densityMax;
};
