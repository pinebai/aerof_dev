#ifndef _NAVIER_STOKES_EMBEDDED_COUPLEDSOLVER_H_
#define _NAVIER_STOKES_EMBEDDED_COUPLEDSOLVER_H_
#include <IoData.h>
#include <GeoSource.h>
#include <Domain.h>
#include <TsSolver.h>
#include <ExplicitEmbeddedTsDesc.h>
#include <MultiGridEmbeddedTsDesc.h>
#include <ImplicitEmbeddedCoupledTsDesc.h>
#include <MultiGridSolver.h>
#include <EmbeddedFluidShapeOptimizationHandler.h>

template<int dim>
void startNavierStokesEmbeddedCoupledSolver(IoData &ioData, GeoSource &geoSource, Domain &domain)
{

  Communicator *com = domain.getCommunicator();

  domain.createVecPat(dim, &ioData);
  domain.createRhsPat(dim, ioData);

  //Combined calculations of Steady state and sensitivities
  if (ioData.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_) {
    EmbeddedFluidShapeOptimizationHandler<dim> fsoh(ioData, geoSource, &domain);
    TsSolver<EmbeddedFluidShapeOptimizationHandler<dim> > tsSolver(&fsoh);
    tsSolver.fsoSolve(ioData);
  }
  //Sensitivity analysis on a provided steady state solution
  else if (ioData.problem.alltype == ProblemData::_SENSITIVITY_ANALYSIS_) {
    EmbeddedFluidShapeOptimizationHandler<dim> fsoh(ioData, geoSource, &domain);
    TsSolver<EmbeddedFluidShapeOptimizationHandler<dim> > tsSolver(&fsoh);
    tsSolver.fsaSolve(ioData);
  }

  else if (ioData.problem.solutionMethod == ProblemData::TIMESTEPPING) {

    if (ioData.ts.type == TsData::IMPLICIT) {

      ImplicitEmbeddedCoupledTsDesc<dim> tsDesc(ioData, geoSource, &domain);
      TsSolver<ImplicitEmbeddedCoupledTsDesc<dim> > tsSolver(&tsDesc);
      tsSolver.solve(ioData);
    }
    else{
      ExplicitEmbeddedTsDesc<dim> tsDesc(ioData, geoSource, &domain);
      TsSolver<ExplicitEmbeddedTsDesc<dim> > tsSolver(&tsDesc);
      tsSolver.solve(ioData);
    }
  } else {

    MultiGridEmbeddedTsDesc<dim> tsDesc(ioData, geoSource, &domain);
    MultiGridSolver<MultiGridEmbeddedTsDesc<dim> > mgSolver(&tsDesc);
    mgSolver.solve(ioData);
  }

}

#endif
