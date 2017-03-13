#pragma once

#include <ImplicitCoupledTsDesc.h>

#include <MultiGridKernel.h>

#include <MultiGridDistSVec.h>

#include <MultiGridSpaceOperator.h>

#include <MultiGridSmoothingMatrices.h>

#include <MultiGridKspSolver.h>

template <int dim>
class MultiGridCoupledTsDesc : public ImplicitCoupledTsDesc<dim> {

 public:

  MultiGridCoupledTsDesc(IoData &, GeoSource &, Domain *); 

  ~MultiGridCoupledTsDesc();

  void smooth0(DistSVec<double,dim>& x,int steps);

  void smooth(int lvl, MultiGridDistSVec<double,dim>& x,
       DistSVec<double,dim>& f,int steps);

  void cycle(int lvl, DistSVec<double,dim>& f,
             MultiGridDistSVec<double,dim>& x);

  void cycle(DistSVec<double,dim>& x);
  
  void setupTimeStepping(DistSVec<double,dim> *U, IoData &iod);

  DistSVec<double,dim>* getSmoothedVec() { return NULL; }

 private:

  int numSmooths_pre[6],numSmooths_post[6];

  bool smoothWithGMRES;

  int mc;

  double prolong_relax_factor,restrict_relax_factor;

  MultiGridKernel<double>* pKernel;

  MultiGridDistSVec<double,dim> V, res, R, F, U,Uold,dx,Forig;

  MultiGridSpaceOperator<double,dim>* mgSpaceOp;

  MultiGridSmoothingMatrices<double,dim>* smoothingMatrices;

  MultiGridMvpMatrix<double,dim>* mgMvp;

  MultiGridKspSolver<double,dim,double>* mgKspSolver;

  int globalIt;
};
