#ifndef _NAVIER_STOKES_EMBEDDED_SEGSOLVER_H_
#define _NAVIER_STOKES_EMBEDDED_SEGSOLVER_H_
#include <IoData.h>
#include <GeoSource.h>
#include <Domain.h>
#include <TsSolver.h>
#include <ImplicitEmbeddedSegTsDesc.h>
template<int dim, int neq1, int neq2>
void startNavierStokesEmbeddedSegSolver(IoData &ioData, GeoSource &geoSource, Domain &domain)
{

  Communicator *com = domain.getCommunicator();

  domain.createVecPat(dim, &ioData);
  domain.createRhsPat(dim, ioData);

  ImplicitEmbeddedSegTsDesc<dim,neq1,neq2> tsDesc(ioData, geoSource, &domain);
  TsSolver<ImplicitEmbeddedSegTsDesc<dim,neq1,neq2> > tsSolver(&tsDesc);
  tsSolver.solve(ioData);

}

#endif
