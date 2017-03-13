#ifndef _NAVIER_STOKES_EMBEDDED_H_
#define _NAVIER_STOKES_EMBEDDED_H_
#include <IoData.h>
#include <GeoSource.h>
#include <Domain.h>
#include <TsSolver.h>
#include <ExplicitMultiPhysicsTsDesc.h>
#include <ImplicitMultiPhysicsTsDesc.h>
template<int dim, int dimLS>
void startNavierStokesMultiPhysicsEmbedded(IoData &ioData, GeoSource &geoSource, Domain &domain)
{

  Communicator *com = domain.getCommunicator();

  domain.createVecPat(dim, &ioData);
  domain.createPhiVecPat(dimLS, &ioData);
  domain.createRhsPat(dim, ioData);

  if (ioData.ts.type == TsData::IMPLICIT) {
    ImplicitMultiPhysicsTsDesc<dim,dimLS> tsDesc(ioData, geoSource, &domain);
    TsSolver<ImplicitMultiPhysicsTsDesc<dim,dimLS> > tsSolver(&tsDesc);
    tsSolver.solve(ioData);

  }
  else{
    ExplicitMultiPhysicsTsDesc<dim,dimLS> tsDesc(ioData, geoSource, &domain);
    TsSolver<ExplicitMultiPhysicsTsDesc<dim,dimLS> > tsSolver(&tsDesc);
    tsSolver.solve(ioData);
  }

}

#endif
