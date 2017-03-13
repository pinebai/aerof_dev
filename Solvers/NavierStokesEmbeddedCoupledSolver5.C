#include "Solvers.h"
#include "NavierStokesEmbeddedCoupledSolver.h"

template <>
void
NavierStokesEmbeddedCoupledSolver<5>::
  solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  startNavierStokesEmbeddedCoupledSolver<5>(ioData, geoSource, domain);
}
