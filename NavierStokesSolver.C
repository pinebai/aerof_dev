#include <IoData.h>
#include <Domain.h>
#include "Solvers/Solvers.h"
#include <SparseGridGeneratorDesc.h>
#include <ImplicitRomTsDesc.h>

//-----------------------------------------------------------------------------

void startNavierStokesSolver(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  Communicator* com = domain.getCommunicator();
  int numBurnableFluids = ProgrammedBurn::countBurnableFluids(ioData);
  if (ioData.problem.framework==ProblemData::EMBEDDED || ioData.problem.framework==ProblemData::EMBEDDEDALE) { //Trigger the embedded framework
    if (ioData.eqs.type == EquationsData::EULER) {
		//temporary testing for embeddedALSonline
		if(ioData.problem.alltype == ProblemData::_EMBEDDED_ALS_ROM_ONLINE_){
			com->fprintf(stderr, "*** Testing Embedded ALS Online ***\n");
			NavierStokesCoupledSolver<5>::solve(ioData, geoSource, domain);
		}
      com->fprintf(stderr, "*** Running an Embedded Inviscid %d Phase Fluid-Structure Simulation with %d Level-Set(s) ***\n", ioData.eqs.numPhase, ioData.embed.nLevelset);
      switch(ioData.embed.nLevelset) {
        case 0 : NavierStokesEmbeddedCoupledSolver<5>::solve(ioData, geoSource, domain); break;
        case 1 : NavierStokesMultiPhysicsEmbedded<5,1>::solve(ioData,geoSource,domain); break;
        case 2 : NavierStokesMultiPhysicsEmbedded<5,2>::solve(ioData,geoSource,domain); break;
        case 3 : NavierStokesMultiPhysicsEmbedded<5,3>::solve(ioData,geoSource,domain); break;
        // Feel free to add more here. e.g. case 4 : NavierStokesMultiPhysicsEmbedded<5,4>::solve(ioData,geoSource,domain); break;
        default: 
          com->fprintf(stderr,"*** Error: %d level-sets detected. Only support 0 ~ 3 for now although it can be extended quickly.\n", 
                       ioData.embed.nLevelset); 
          exit(-1);
      }
    }
    else if (ioData.eqs.type == EquationsData::NAVIER_STOKES)
      {
        com->fprintf(stderr, "*** Running an Embedded Viscous %d Phase Fluid-Structure simulation ***\n", ioData.eqs.numPhase);
	if(ioData.eqs.tc.type == TurbulenceClosureData::NONE)
	  {
	    com->fprintf(stderr,"--- No Turbulent Model Used ***\n");
	    NavierStokesEmbeddedCoupledSolver<5>::solve(ioData, geoSource, domain);
	  }
	else if(ioData.eqs.tc.type == TurbulenceClosureData::LES)
	  {
	    if(ioData.eqs.tc.les.type == LESModelData::SMAGORINSKY)
	    {
	      com->fprintf(stderr,"--- Smagorinsky LES Model Used ***\n");
	      NavierStokesEmbeddedCoupledSolver<5>::solve(ioData, geoSource, domain);
	    }
	    if(ioData.eqs.tc.les.type == LESModelData::DYNAMIC)
	    {
	      com->fprintf(stderr,"--- Dynamic LES Model Used ***\n");
	      NavierStokesEmbeddedCoupledSolver<5>::solve(ioData, geoSource, domain);
	    }
	    if(ioData.eqs.tc.les.type == LESModelData::WALE)
	    {
	      com->fprintf(stderr,"--- Wale LES Model Used ***\n");
	      NavierStokesEmbeddedCoupledSolver<5>::solve(ioData, geoSource, domain);
	    }
	  }
	else if(ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY)
	  {
	    if(ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS)
	      {
		com->fprintf(stderr,"--- Spalart-Allmaras Turbulent Model Used ***\n");
		if (ioData.ts.type == TsData::IMPLICIT &&
				ioData.ts.implicit.tmcoupling == ImplicitData::WEAK)
		  NavierStokesEmbeddedSegSolver<6,5,1>::solve(ioData, geoSource, domain);
		else
		  NavierStokesEmbeddedCoupledSolver<6>::solve(ioData, geoSource, domain);
	      }
	    if(ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES)
	      {
		com->fprintf(stderr,"--- DES Turbulent Model Used ***\n");
					//com->fprintf(stderr, "*** Error: This option should work but has never been tested. Please use carefully. ***\n");
					//exit(1);
		if (ioData.ts.type == TsData::IMPLICIT &&
				ioData.ts.implicit.tmcoupling == ImplicitData::WEAK)
		  NavierStokesEmbeddedSegSolver<6,5,1>::solve(ioData, geoSource, domain);
		else
		  NavierStokesEmbeddedCoupledSolver<6>::solve(ioData, geoSource, domain);
	      }
	    if(ioData.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE)
	      {
		com->fprintf(stderr,"--- K-Epsilon Turbulent Model Used ***\n");
		if (ioData.ts.type == TsData::IMPLICIT &&
				ioData.ts.implicit.tmcoupling == ImplicitData::WEAK)
		  NavierStokesEmbeddedSegSolver<7,5,2>::solve(ioData, geoSource, domain);
		else
		  NavierStokesEmbeddedCoupledSolver<7>::solve(ioData, geoSource, domain);
	      }
	  }
	else
	  {
	    com->fprintf(stderr, "*** Error: wrong turbulence closure type\n");
	    exit(1);
	  }
      }
    else
      {
	com->fprintf(stderr, "*** Error: wrong equation type\n");
	exit(1);
      }
  } 
  else if (ioData.eqs.numPhase == 1){
    if (ioData.eqs.type == EquationsData::EULER) 
      NavierStokesCoupledSolver<5>::solve(ioData, geoSource, domain);
    else if (ioData.eqs.type == EquationsData::NAVIER_STOKES) {
      if (ioData.eqs.tc.type == TurbulenceClosureData::NONE ||
	  ioData.eqs.tc.type == TurbulenceClosureData::LES)
	//startNavierStokesCoupledSolver<5>(ioData, geoSource, domain);
	NavierStokesCoupledSolver<5>::solve(ioData, geoSource, domain);
      else if (ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
	if (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
	    ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES) {
		if (ioData.ts.type == TsData::IMPLICIT &&
				ioData.ts.implicit.tmcoupling == ImplicitData::WEAK)
			if (ioData.problem.alltype == ProblemData::_STEADY_NONLINEAR_ROM_ ||
                            ioData.problem.alltype == ProblemData::_UNSTEADY_NONLINEAR_ROM_ || 
                            ioData.problem.alltype == ProblemData::_STEADY_NONLINEAR_ROM_POST_ ||
                            ioData.problem.alltype == ProblemData::_UNSTEADY_NONLINEAR_ROM_POST_ ||
                            ioData.problem.alltype == ProblemData::_ACC_UNSTEADY_NONLINEAR_ROM_ ||
                            ioData.problem.alltype == ProblemData::_FORCED_NONLINEAR_ROM_) {
				com->fprintf(stderr,"*** WARNING: Seg solver not implemented for UnsteadyROM, starting the coupled solver\n"); //CBM
				NavierStokesCoupledSolver<6>::solve(ioData, geoSource, domain);
			}
			else
	   // startNavierStokesSegSolver<6,5,1>(ioData, geoSource, domain);
	    NavierStokesSegSolver<6,5,1>::solve(ioData, geoSource, domain);
	  else
	    NavierStokesCoupledSolver<6>::solve(ioData, geoSource, domain);
	}
	else if (ioData.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE) {
	  if (ioData.ts.type == TsData::IMPLICIT &&
	      ioData.ts.implicit.tmcoupling == ImplicitData::WEAK)
	    NavierStokesSegSolver<7,5,2>::solve(ioData, geoSource, domain);
	  else
	    NavierStokesCoupledSolver<7>::solve(ioData, geoSource, domain);
	}
	else {
	  com->fprintf(stderr, "*** Error: wrong turbulence model type\n");
	  exit(1);
	}
      }
      else {
	com->fprintf(stderr, "*** Error: wrong turbulence closure type\n");
	exit(1);
      }
    }
    else {
      com->fprintf(stderr, "*** Error: wrong equation type\n");
      exit(1);
    }
  }
  else if (ioData.eqs.numPhase == 2) {
    com->fprintf(stdout, "*** Running a Multi(%d)-Phase Flow simulation ***\n", ioData.eqs.numPhase);
    LevelSetSolver<5,1>::solve(ioData, geoSource, domain);
  }
  else if (ioData.eqs.numPhase == 3) {
    com->fprintf(stdout, "*** Running a Multi(%d)-Phase Flow simulation ***\n", ioData.eqs.numPhase);
    if (numBurnableFluids == 1)
      LevelSetSolver<5,1>::solve(ioData, geoSource, domain);
    else if (numBurnableFluids == 0)
      LevelSetSolver<5,2>::solve(ioData, geoSource, domain);
    else {
      fprintf(stderr,"Num burnable fluids = %d???\n", numBurnableFluids);
    }
  }
  else if (ioData.eqs.numPhase == 4) {
    com->fprintf(stdout, "*** Running a Multi(%d)-Phase Flow simulation ***\n", ioData.eqs.numPhase);
    if (numBurnableFluids == 1)
      LevelSetSolver<5,2>::solve(ioData, geoSource, domain);
    else if (numBurnableFluids == 2)
      LevelSetSolver<5,1>::solve(ioData, geoSource, domain);
    else if (numBurnableFluids == 0)
      LevelSetSolver<5,3>::solve(ioData, geoSource, domain);
    else {
      fprintf(stderr,"Num burnable fluids = %d???\n", numBurnableFluids);
    }
  }
  else{
    com->fprintf(stdout, "*** Error: a Multi(%d)-Phase Flow simulation cannot be done yet\n", ioData.eqs.numPhase);
  }
}

//------------------------------------------------------------------------------

void startSparseGridGeneration(IoData &ioData, Domain &domain)
{

  Communicator* com = domain.getCommunicator();

  fprintf(stdout, "*** Generating a sparse grid\n");
  SparseGridGeneratorDesc sgDesc(ioData, com);
  sgDesc.tabulate(ioData);
  fprintf(stdout, "*** The sparse grid was generated and the simulation is finished\n");

}

//------------------------------------------------------------------------------
