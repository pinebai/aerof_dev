#include "Solvers/Solvers.h"
#include "Solvers/NavierStokesCoupledSolver.h"

void 
NavierStokesCoupledSolver<7>
  ::solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
      startNavierStokesCoupledSolver<7>(ioData, geoSource, domain);
}
