#include <DistExtrapolation.h>

#include <IoData.h>
#include <Extrapolation.h>
#include <Domain.h>
#include <VarFcn.h>
#include <DistVector.h>
#include <MacroCell.h>

#include <Communicator.h>

#include <cstdio>

//------------------------------------------------------------------------------

template<int dim>
DistExtrapolation<dim>::DistExtrapolation(IoData& iod, Domain* domain, VarFcn* vf)
{

  numLocSub  = domain->getNumLocSub();
  subDomain  = domain->getSubDomain();
  it0        = iod.restart.iteration;
  lastIt     = it0;
  lastConfig = -1;

  subExtrapolation = new Extrapolation<dim>*[numLocSub];
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subExtrapolation[iSub] = new Extrapolation<dim>(iod, vf);


  myCom = domain->getCommunicator();
}

//------------------------------------------------------------------------------

template<int dim>
DistExtrapolation<dim>::~DistExtrapolation()
{

  if (subExtrapolation) { 
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) 
      if (subExtrapolation[iSub]) 
	delete subExtrapolation[iSub];
    delete [] subExtrapolation;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void DistExtrapolation<dim>::compute(int config, DistVec<Vec3D> &norm, DistSVec<double,3> &X)
{
  if ((config != lastConfig) || (lastIt == it0)) {
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->findNormalTetrahedra( X(iSub), norm(iSub), subExtrapolation[iSub]->getExtrapolationData());
    if (lastConfig != config)
      lastConfig = config;
    if (lastIt == it0)
      lastIt = -1;
  }
  
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void DistExtrapolation<dim>::computeDerivative(int config, DistVec<Vec3D> &norm, DistSVec<double,3> &X)
{

  fprintf(stderr, "***** DistExtrapolation<dim>::computeDerivative is not implemented!\n");
  exit(1);
  
}

//------------------------------------------------------------------------------
