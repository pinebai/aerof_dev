#ifndef _LEVELSETSOLVER_H_
#define _LEVELSETSOLVER_H_
#include <IoData.h>
#include <GeoSource.h>
#include <Domain.h>
#include <TsSolver.h>
#include <ImplicitLevelSetTsDesc.h>
#include <ExplicitLevelSetTsDesc.h>


template<int dim, int dimLS>
void startLevelSetSolver(IoData &ioData, GeoSource &geoSource, Domain &domain)
{

  Communicator *com = domain.getCommunicator();

  domain.createVecPat(dim, &ioData);
  domain.createPhiVecPat(dimLS, &ioData);
  domain.createRhsPat(dim, ioData);

  if (ioData.ts.type == TsData::IMPLICIT) {
    ImplicitLevelSetTsDesc<dim,dimLS> tsDesc(ioData, geoSource, &domain);
    TsSolver<ImplicitLevelSetTsDesc<dim,dimLS> > tsSolver(&tsDesc);
    tsSolver.solve(ioData);
  }
  else{
    ExplicitLevelSetTsDesc<dim,dimLS> tsDesc(ioData, geoSource, &domain);
    TsSolver<ExplicitLevelSetTsDesc<dim,dimLS> > tsSolver(&tsDesc);
    tsSolver.solve(ioData);
  }

}

#endif
