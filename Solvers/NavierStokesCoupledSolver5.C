#include "Solvers/Solvers.h"
#include "Solvers/NavierStokesCoupledSolver.h"

void 
NavierStokesCoupledSolver<5>
  ::solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
       std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
      startNavierStokesCoupledSolver<5>(ioData, geoSource, domain);
}
