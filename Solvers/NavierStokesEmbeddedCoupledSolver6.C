#include "Solvers.h"
#include "NavierStokesEmbeddedCoupledSolver.h"

template <>
void
NavierStokesEmbeddedCoupledSolver<6>::
  solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  startNavierStokesEmbeddedCoupledSolver<6>(ioData, geoSource, domain);
}
