#include "Solvers/Solvers.h"
#include "Solvers/NavierStokesCoupledSolver.h"

void 
NavierStokesCoupledSolver<5>
  ::solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
      startNavierStokesCoupledSolver<5>(ioData, geoSource, domain);
}
