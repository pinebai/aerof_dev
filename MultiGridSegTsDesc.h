#pragma once

#include <ImplicitSegTsDesc.h>

#include <MultiGridKernel.h>

#include <MultiGridDistSVec.h>

#include <MultiGridSpaceOperator.h>

#include <MultiGridSmoothingMatrices.h>

#include <MultiGridKspSolver.h>


template <int dim,int neq1,int neq2>
class MultiGridSegTsDesc : public ImplicitSegTsDesc<dim,neq1,neq2> {

 public:

  MultiGridSegTsDesc(IoData &, GeoSource &, Domain *); 

  ~MultiGridSegTsDesc();

  void smooth0(DistSVec<double,dim>& x,int steps);

  void smooth(int lvl, MultiGridDistSVec<double,dim>& x,
	      DistSVec<double,dim>& f,int steps,
	      bool postsmooth);

  void cycle(int lvl, DistSVec<double,dim>& f,
             MultiGridDistSVec<double,dim>& x);

  void cycle(DistSVec<double,dim>& x);
  
  void setupTimeStepping(DistSVec<double,dim> *U, IoData &iod);

  DistSVec<double,dim>* getSmoothedVec() { return NULL; }

 private:

  int numSmooths_pre[6],numSmooths_post[6];

  bool smoothWithGMRES;

  int mc;

  int globalIt;

  int addTurbulenceTerms;

  double prolong_relax_factor,restrict_relax_factor;
  
  MultiGridKernel<double>* pKernel;

  MultiGridDistSVec<double,dim> V, res, R, F, U,Uold,dx,Forig,update_tmp;
  MultiGridDistSVec<double,neq1> R1,dx1;
  MultiGridDistSVec<double,neq2> R2,dx2;

  MultiGridSpaceOperator<double,dim>* mgSpaceOp;

  MultiGridSmoothingMatrices<double,neq1>* smoothingMatrices1;
  MultiGridSmoothingMatrices<double,neq2>* smoothingMatrices2;

  MultiGridMvpMatrix<double,neq1>* mgMvp1;
  MultiGridMvpMatrix<double,neq2>* mgMvp2;

  MultiGridKspSolver<double,neq1,double>* mgKspSolver1;
  MultiGridKspSolver<double,neq2,double>* mgKspSolver2;
};
