#include "Solvers.h"
#include "NavierStokesEmbeddedCoupledSolver.h"

template <>
void
NavierStokesEmbeddedCoupledSolver<7>::
  solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  startNavierStokesEmbeddedCoupledSolver<7>(ioData, geoSource, domain);
}
