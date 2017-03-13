#include "Solvers.h"
#include "NavierStokesEmbeddedSegSolver.h"

template<>
void
NavierStokesEmbeddedSegSolver<6,5,1>::
  solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  startNavierStokesEmbeddedSegSolver<6,5,1>(ioData, geoSource, domain);
}
