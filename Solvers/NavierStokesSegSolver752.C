#include "Solvers.h"
#include "NavierStokesSegSolver.h"

template<>
void
NavierStokesSegSolver<7,5,2>::
  solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  startNavierStokesSegSolver<7,5,2>(ioData, geoSource, domain);
}
