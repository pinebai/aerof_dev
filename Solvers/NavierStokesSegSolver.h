#ifndef _NAVIERSTOKESSEGSOLVER_H_
#define _NAVIERSTOKESSEGSOLVER_H_
#include <IoData.h>
#include <GeoSource.h>
#include <Domain.h>
#include <TsSolver.h>
#include <ImplicitSegTsDesc.h>
#include <MultiGridSegTsDesc.h>
#include <MultiGridSolver.h>
#include <FluidSegShapeOptimizationHandler.h> //lscheuch //TODO delete line if unused

template<int dim, int neq1, int neq2>
void startNavierStokesSegSolver(IoData &ioData, GeoSource &geoSource, Domain &domain)
{

  domain.createVecPat(dim);
  domain.createRhsPat(dim, ioData);

  if (ioData.problem.solutionMethod == ProblemData::TIMESTEPPING) {
    ImplicitSegTsDesc<dim,neq1,neq2> tsDesc(ioData, geoSource, &domain);

    //TODO HACK
    if (ioData.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_){
      FluidSegShapeOptimizationHandler<dim,neq1,neq2> fsoh(ioData, geoSource, &domain);
      TsSolver<FluidSegShapeOptimizationHandler<dim,neq1,neq2> > tsSolver(&fsoh);
      tsSolver.fsoSolve(ioData);
    }
    //TODO HACK
    else if (ioData.problem.alltype == ProblemData::_SENSITIVITY_ANALYSIS_){
      FluidSegShapeOptimizationHandler<dim,neq1,neq2> fsoh(ioData, geoSource, &domain);
      TsSolver<FluidSegShapeOptimizationHandler<dim,neq1,neq2> > tsSolver(&fsoh);
      tsSolver.fsaSolve(ioData);
    }
    else{
      TsSolver<ImplicitSegTsDesc<dim,neq1,neq2> > tsSolver(&tsDesc);
      tsSolver.solve(ioData);
    }
  } else {

    MultiGridSegTsDesc<dim,neq1,neq2> tsDesc(ioData, geoSource, &domain);
    MultiGridSolver<MultiGridSegTsDesc<dim,neq1,neq2> > mgSolver(&tsDesc);
    mgSolver.solve(ioData);

  }
}

#endif
