#include <EdgeGalerkin.h>

#include <RefVal.h>
#include <ViscoFcn.h>
#include <ThermalCondFcn.h>
#include <DistBcData.h>
#include <SubDomain.h>
#include <DistVector.h>

//------------------------------------------------------------------------------

template<int dim>
void EdgeGalerkin::compute(int config, RefVal *refVal, VarFcn *varFcn, ViscoFcn *viscoFcn, 
			   ThermalCondFcn *thermalCondFcn, DistBcData<dim> &bCond, 
			   DistSVec<double,3> &X, DistVec<double> &ctrlVol, 
			   DistSVec<double,dim> &V, DistSVec<double,dim> &F)
{

  int iSub;

  if (config != lastConfig) {
    double t0 = timer->getTime();

#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub) {
      subDomain[iSub]->computeEdgeWeightsGalerkin(X(iSub), ctrlVol(iSub), M(iSub));
      subDomain[iSub]->sndEdgeData(*cp, M.subData(iSub));
    }

    cp->exchange();

#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub)
      subDomain[iSub]->addRcvEdgeData(*cp, M.subData(iSub));

    //timer->addEdgeWeightsTime(t0);

    lastConfig = config;
  }

  double t0 = timer->getTime();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeViscousFluxes2(refVal, varFcn, viscoFcn, thermalCondFcn,
					   bCond(iSub), M(iSub), X(iSub), V(iSub), F(iSub));

  timer->addVisFluxesTime(t0);

}

//------------------------------------------------------------------------------

