#include "Solvers.h"
#include "LevelSetSolver.h"

template <>
void
LevelSetSolver<5,1>::
  solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  startLevelSetSolver<5,1>(ioData, geoSource, domain);
}

template <>
void
LevelSetSolver<5,2>::
  solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  startLevelSetSolver<5,2>(ioData, geoSource, domain);
}

template <>
void
LevelSetSolver<5,3>::
  solve(IoData &ioData, GeoSource &geoSource, Domain &domain)
{
  startLevelSetSolver<5,3>(ioData, geoSource, domain);
}
