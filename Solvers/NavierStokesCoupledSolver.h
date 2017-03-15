#ifndef _NAVIERSTOKESCOUPLEDSOLVER_H_
#define _NAVIERSTOKESCOUPLEDSOLVER_H_

#include <IoData.h>
#include <GeoSource.h>
#include <Domain.h>
#include <TsSolver.h>
#include <ExplicitTsDesc.h>
#include <ImplicitCoupledTsDesc.h>
#include <ImplicitPGTsDesc.h>
//#include <ImplicitGalerkinTsDesc.h>
#include <ImplicitGnatTsDesc.h>
#include <ImplicitCollocationTsDesc.h>
#include <ImplicitMetricTsDesc.h>
#include <ImplicitRomPostproTsDesc.h>
#include <MultiGridSolver.h>
#include <MultiGridCoupledTsDesc.h>
#include <FluidShapeOptimizationHandler.h>
#include <FluidRomShapeOptimizationHandler.h>  // MZ
#include <ImplicitEmbeddedRomTsDesc.h>          // Lei Lei
#include <FluidGnatShapeOptimizationHandler.h>  // MZ
#include <FluidMetricShapeOptimizationHandler.h>  // MZ
#include <FluidCollocationShapeOptimizationHandler.h>  // MZ

template<int dim>
void startNavierStokesCoupledSolver(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  Communicator* com = domain.getCommunicator();

  domain.createVecPat(dim, &ioData);
  domain.createRhsPat(dim, ioData);

  //Combined calculations of Steady state and Sensitivities
  if (ioData.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_) { // YC
      FluidShapeOptimizationHandler<dim> fsoh(ioData, geoSource, &domain);
      TsSolver<FluidShapeOptimizationHandler<dim> > tsSolver(&fsoh);
      tsSolver.fsoSolve(ioData);
  }
  //Sensitivity analysis on a provided steady state solution
  else if (ioData.problem.alltype == ProblemData::_SENSITIVITY_ANALYSIS_) {
      FluidShapeOptimizationHandler<dim> fsoh(ioData, geoSource, &domain);

//Ori
      TsSolver<FluidShapeOptimizationHandler<dim> > tsSolver(&fsoh);
      tsSolver.fsaSolve(ioData);

      //TODO HACK
      //fsoh.fsaOnlySolve(ioData);
  }
  else if (ioData.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_) { // YC
      FluidShapeOptimizationHandler<dim> fsisoh(ioData, geoSource, &domain);
      TsSolver<FluidShapeOptimizationHandler<dim> > tsSolver(&fsisoh);
      tsSolver.fsisoSolve(ioData);
  }
  else if (ioData.problem.alltype == ProblemData::_ROM_SHAPE_OPTIMIZATION_) { // MZ
    if (ioData.romOnline.projection == NonlinearRomOnlineData::PETROV_GALERKIN && ioData.romOnline.systemApproximation == NonlinearRomOnlineData::SYSTEM_APPROXIMATION_NONE) {
      FluidRomShapeOptimizationHandler<dim> fsoh(ioData, geoSource, &domain);
      TsSolver<FluidRomShapeOptimizationHandler<dim> > tsSolver(&fsoh);
      tsSolver.fsoSolve(ioData);
    } else if (ioData.romOnline.projection == NonlinearRomOnlineData::PETROV_GALERKIN && ioData.romOnline.systemApproximation == NonlinearRomOnlineData::GNAT) {
      FluidGnatShapeOptimizationHandler<dim> fsoh(ioData, geoSource, &domain);
      TsSolver<FluidGnatShapeOptimizationHandler<dim> > tsSolver(&fsoh);
      tsSolver.fsoSolve(ioData);
    } else if (ioData.romOnline.projection == NonlinearRomOnlineData::PETROV_GALERKIN && ioData.romOnline.systemApproximation == NonlinearRomOnlineData::COLLOCATION) {
      FluidCollocationShapeOptimizationHandler<dim> fsoh(ioData, geoSource, &domain);
      TsSolver<FluidCollocationShapeOptimizationHandler<dim> > tsSolver(&fsoh);
      tsSolver.fsoSolve(ioData);
    } else if (ioData.romOnline.projection == NonlinearRomOnlineData::PETROV_GALERKIN && ioData.romOnline.systemApproximation == NonlinearRomOnlineData::APPROX_METRIC_NL) {
      FluidMetricShapeOptimizationHandler<dim> fsoh(ioData, geoSource, &domain);
      TsSolver<FluidMetricShapeOptimizationHandler<dim> > tsSolver(&fsoh);
      tsSolver.fsoSolve(ioData);
    } else {
      com->fprintf(stderr, "*** Error: this system approximation method is not currently supported forROM Shape optimization\n");
      exit(-1);
    }
  }
  else if ((ioData.problem.alltype == ProblemData::_STEADY_NONLINEAR_ROM_) ||
           (ioData.problem.alltype == ProblemData::_UNSTEADY_NONLINEAR_ROM_) ||
           (ioData.problem.alltype == ProblemData::_ACC_UNSTEADY_NONLINEAR_ROM_) ||
           (ioData.problem.alltype == ProblemData::_FORCED_NONLINEAR_ROM_)) {
    if (ioData.romOnline.projection == 0 && ioData.romOnline.systemApproximation == NonlinearRomOnlineData::SYSTEM_APPROXIMATION_NONE) {
        ImplicitPGTsDesc<dim> tsDesc(ioData, geoSource, &domain);
        TsSolver<ImplicitPGTsDesc<dim> > tsSolver(&tsDesc);
        tsSolver.solve(ioData);
    }
    else if (ioData.romOnline.projection == 0 && ioData.romOnline.systemApproximation == NonlinearRomOnlineData::GNAT) {
        ImplicitGnatTsDesc<dim> tsDesc(ioData, geoSource, &domain);
        TsSolver<ImplicitGnatTsDesc<dim> > tsSolver(&tsDesc);
        tsSolver.solve(ioData);
    }
    else if (ioData.romOnline.projection == 0 && ioData.romOnline.systemApproximation == NonlinearRomOnlineData::COLLOCATION) {
        ImplicitCollocationTsDesc<dim> tsDesc(ioData, geoSource, &domain);
        TsSolver<ImplicitCollocationTsDesc<dim> > tsSolver(&tsDesc);
        tsSolver.solve(ioData);
    }
    else if (ioData.romOnline.projection == 0 && ioData.romOnline.systemApproximation == NonlinearRomOnlineData::APPROX_METRIC_NL) {
        ImplicitMetricTsDesc<dim> tsDesc(ioData, geoSource, &domain);
        TsSolver<ImplicitMetricTsDesc<dim> > tsSolver(&tsDesc);
        tsSolver.solve(ioData);
    }
			/*else if (ioData.rom.projection == 1 && ioData.rom.systemApproximation == 0) {
      ImplicitGalerkinTsDesc<dim> tsDesc(ioData, geoSource, &domain);
      TsSolver<ImplicitGalerkinTsDesc<dim> > tsSolver(&tsDesc);
      tsSolver.solve(ioData);
			}*/
    else
        com->fprintf(stderr, "*** Error: this type of nonlinear ROM simulation is not currently supported\n");
  }
  else if(ioData.problem.alltype == ProblemData::_EMBEDDED_ALS_ROM_ONLINE_){
      ImplicitEmbeddedRomTsDesc<dim> tsDesc(ioData, geoSource, &domain);
      TsSolver<ImplicitEmbeddedRomTsDesc<dim> > tsSolver(&tsDesc);
      // todo: testing functionality of tsDesc
      //tsDesc.test();
      tsSolver.solve(ioData);
  }
  else if (ioData.problem.alltype == ProblemData::_UNSTEADY_NONLINEAR_ROM_POST_ ||
           ioData.problem.alltype == ProblemData::_STEADY_NONLINEAR_ROM_POST_) {
      ImplicitRomPostproTsDesc <dim> tsDesc(ioData, geoSource, &domain);
      TsSolver<ImplicitRomPostproTsDesc<dim> > tsSolver(&tsDesc);
      tsSolver.solve(ioData);
  }
  else if (ioData.ts.type == TsData::IMPLICIT) {
    if (ioData.problem.solutionMethod == ProblemData::TIMESTEPPING) {
      //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;//TODO delete line
      ImplicitCoupledTsDesc<dim> tsDesc(ioData, geoSource, &domain);
      //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;//TODO delete line
      TsSolver<ImplicitCoupledTsDesc<dim> > tsSolver(&tsDesc);
      //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;//TODO delete line
      tsSolver.solve(ioData);
      //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;//TODO delete line
    } else {
      //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;//TODO delete line
      MultiGridCoupledTsDesc<dim> tsDesc(ioData, geoSource, &domain);
      MultiGridSolver<MultiGridCoupledTsDesc<dim> > mgSolver(&tsDesc);
      mgSolver.solve(ioData);
    }
  }
  else if (ioData.ts.type == TsData::EXPLICIT) {
      std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
    ExplicitTsDesc<dim> tsDesc(ioData, geoSource, &domain);
    TsSolver<ExplicitTsDesc<dim> > tsSolver(&tsDesc);
    tsSolver.solve(ioData);
  }
  else
    com->fprintf(stderr, "*** Error: wrong time-integrator\n");

}

#endif
