#include <IoData.h>
#include <Communicator.h>
#include <Domain.h>
#include <Modal.h>

//------------------------------------------------------------------------------
void startModalSolver(Communicator *com, IoData &ioData, Domain &domain)
{

   if (ioData.eqs.type == EquationsData::EULER ) {
    ModalSolver<5> mSolver(com, ioData, domain);
    domain.createVecPat(5, &ioData);
    mSolver.solve();
  }
  else if (ioData.eqs.type == EquationsData::NAVIER_STOKES) {
    if (ioData.eqs.tc.type == TurbulenceClosureData::NONE ||
        ioData.eqs.tc.type == TurbulenceClosureData::LES) {
      ModalSolver<5> mSolver(com, ioData, domain);
      domain.createVecPat(5, &ioData);
      mSolver.solve();
    }
    else if (ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
      if (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
            ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES) {
        if (false/*ioData.ts.type == TsData::IMPLICIT &&
            ioData.ts.implicit.coupling == ImplicitData::WEAK*/) {
          com->fprintf(stderr, "*** Error: weak coupling not supported in model reduction\n");
          exit(1);
        }
        else {
          ModalSolver<6> mSolver(com, ioData, domain);
          domain.createVecPat(6, &ioData);
          mSolver.solve();
        }
      }
    }
  }
/*
      else if (ioData.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE) {
        if (ioData.ts.implicit.coupling == ImplicitData::WEAK) {
          com->fprintf(stderr, "*** Error: weak coupling not supported in model reduction\n");
          exit(1);
        }
          else {
            ModalSolver<7> mSolver(com, ioData, domain);
            domain.createVecPat(7, &ioData);
            mSolver.solve();
          }
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
*/
 
}

//------------------------------------------------------------------------------
