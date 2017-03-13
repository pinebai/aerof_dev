#include "Solvers/Solvers.h"
#include "Solvers/NavierStokesCoupledSolver.h"

void 
NavierStokesCoupledSolver<6>
  ::solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
      startNavierStokesCoupledSolver<6>(ioData, geoSource, domain);
}
