#include "Solvers.h"
#include "NavierStokesSegSolver.h"

template<>
void
NavierStokesSegSolver<6,5,1>::
  solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  startNavierStokesSegSolver<6,5,1>(ioData, geoSource, domain);
}
