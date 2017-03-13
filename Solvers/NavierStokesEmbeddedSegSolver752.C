#include "Solvers.h"
#include "NavierStokesEmbeddedSegSolver.h"

template<>
void
NavierStokesEmbeddedSegSolver<7,5,2>::
  solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  startNavierStokesEmbeddedSegSolver<7,5,2>(ioData, geoSource, domain);
}
