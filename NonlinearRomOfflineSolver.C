#include <IoData.h>
#include <Communicator.h>
#include <Domain.h>
#include <NonlinearRomOffline.h>

//------------------------------------------------------------------------------
void startNonlinearRomOfflineSolver(Communicator *com, IoData &ioData, Domain &domain, GeoSource &geoSource)
{
   if (ioData.eqs.type == EquationsData::EULER ) {
    NonlinearRomOfflineSolver<5> offlineSolver(com, ioData, domain, geoSource);
    domain.createVecPat(5, &ioData);
    offlineSolver.solve();
  }
  else if (ioData.eqs.type == EquationsData::NAVIER_STOKES) {
    if (ioData.eqs.tc.type == TurbulenceClosureData::NONE ||
        ioData.eqs.tc.type == TurbulenceClosureData::LES) {
      NonlinearRomOfflineSolver<5> offlineSolver(com, ioData, domain, geoSource);
      domain.createVecPat(5, &ioData);
      offlineSolver.solve();
    }
    else if (ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
      if (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
            ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES) {
        NonlinearRomOfflineSolver<6> offlineSolver(com, ioData, domain, geoSource);
        domain.createVecPat(6, &ioData);
        offlineSolver.solve();
      }
    }
  }
}

//------------------------------------------------------------------------------

