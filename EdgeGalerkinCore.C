#include <EdgeGalerkin.h>

#include <SubDomain.h>
#include <Domain.h>
#include <Timer.h>

//------------------------------------------------------------------------------

EdgeGalerkin::EdgeGalerkin(IoData &ioData, Domain *domain) : 
  M(domain->getEdgeDistInfo())
{

  numLocSub = domain->getNumLocSub();
  subDomain = domain->getSubDomain();
  timer = domain->getTimer();

  cp = new CommPattern<double>(domain->getSubTopo(), domain->getCommunicator(), 
			       CommPattern<double>::CopyOnSend);

#pragma omp parallel for
  for (int iSub = 0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->setComLenEdges(9, *cp);

  cp->finalize();

  lastConfig = -1;

}

//------------------------------------------------------------------------------

EdgeGalerkin::~EdgeGalerkin()
{

  if (cp) delete cp;
  
}

//------------------------------------------------------------------------------

