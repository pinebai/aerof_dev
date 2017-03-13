#include <DistExactRiemannSolver.h>


#include <ExactRiemannSolver.h>
#include <IoData.h>
#include <Domain.h>
#include <Communicator.h>
#include "SparseGridCluster.h"

#include <cstdio>

//------------------------------------------------------------------------------

template<int dim>
DistExactRiemannSolver<dim>::DistExactRiemannSolver(IoData &ioData, Domain *dom,
                                                    VarFcn *vf) :
domain(dom)
{

  phaseChangeType_ = ioData.mf.typePhaseChange;
  algorithmType_   = ioData.mf.method;

  numLocSub = dom->getNumLocSub();

  oldV = new DistSVec<double,dim>(dom->getNodeDistInfo());
  riemannupdate = new DistSVec<double,dim>(dom->getNodeDistInfo());
  weight        = new DistVec<double>(dom->getNodeDistInfo());

  interfacialWi = new DistSVec<double,dim-2>(dom->getEdgeDistInfo());
  interfacialWj = new DistSVec<double,dim-2>(dom->getEdgeDistInfo());

  fluidIdToSet = new DistVec<int>(dom->getNodeDistInfo());

  if(ioData.mf.riemannComputation == MultiFluidData::TABULATION2){
    // only the ioData.eqs.fluidModel is considered since only the Riemann invariant of one EOS is tabulated!
    // (no need to specify two different EOS)

    double *refIn  = new double[2];
    double *refOut = new double[1];
    refIn[0] = ioData.ref.rv.density;
    refIn[1] = pow(ioData.ref.rv.density,-ioData.eqs.fluidModel.jwlModel.omega)*ioData.ref.rv.velocity*ioData.ref.rv.velocity;
    refOut[0] = ioData.ref.rv.velocity;

    tabulationC = new SparseGridCluster;
    tabulationC->readFromFile(ioData.mf.sparseGrid.numberOfTabulations, refIn, refOut, ioData.mf.sparseGrid.tabulationFileName, dom->getCommunicator()->cpuNum());

  }else if(ioData.mf.riemannComputation == MultiFluidData::TABULATION5){

    double *refIn = new double[5]; double *refOut = new double[2];
    refIn[0] = ioData.ref.rv.density;
    refIn[1] = ioData.ref.rv.pressure;
    refIn[2] = ioData.ref.rv.density;
    refIn[3] = ioData.ref.rv.pressure;
    refIn[4] = ioData.ref.rv.velocity;
    refOut[0] = ioData.ref.rv.density;
    refOut[1] = ioData.ref.rv.density;

    tabulationC = new SparseGridCluster;
    tabulationC->readFromFile(ioData.mf.sparseGrid.numberOfTabulations, refIn, refOut, ioData.mf.sparseGrid.tabulationFileName, dom->getCommunicator()->cpuNum());
  }else{
    tabulationC = 0;
  }

  subExactRiemannSolver = new ExactRiemannSolver<dim>*[numLocSub];
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subExactRiemannSolver[iSub] =
      new ExactRiemannSolver<dim>(ioData, (*riemannupdate)(iSub),
                         	(*weight)(iSub), (*interfacialWi)(iSub), (*interfacialWj)(iSub), vf, tabulationC, (*fluidIdToSet)(iSub));

}
//------------------------------------------------------------------------------
template<int dim>
DistExactRiemannSolver<dim>::~DistExactRiemannSolver()
{

  delete riemannupdate;
  delete weight;
  delete interfacialWi;
  delete interfacialWj;
  delete oldV;

  delete fluidIdToSet;

  if (subExactRiemannSolver) {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      if (subExactRiemannSolver[iSub]) 
        delete subExactRiemannSolver[iSub];

    delete [] subExactRiemannSolver;
  }

  delete tabulationC;

}
//------------------------------------------------------------------------------
template<int dim>
void DistExactRiemannSolver<dim>::updatePhaseChange(DistSVec<double,dim> &V,
                                               DistVec<int> &fluidId, DistVec<int> &fluidIdn)
{

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
   
    int res = subExactRiemannSolver[iSub]->updatePhaseChange(V(iSub), fluidId(iSub), fluidIdn(iSub),
						   domain->getSubDomain()[iSub]->getHigherOrderMF());
    if (res >= 0) {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
      std::cout << "Phase change update failed at node " << 
        domain->getSubDomain()[iSub]->getNodeMap()[res]+1 << "; local id = " << res << "; rank = " << rank << std::endl;
      exit(-1);
    }
  }

}
//------------------------------------------------------------------------------
template<int dim>
void DistExactRiemannSolver<dim>::storeOldV(DistSVec<double,dim> &V) {

  *oldV = V;
}

//------------------------------------------------------------------------------
template<int dim>
void DistExactRiemannSolver<dim>::storePreviousPrimitive(DistSVec<double,dim> &V,
                                                    DistVec<int> &fluidId,
                                                    DistSVec<double,3> &X)
{
  /*
  if(phaseChangeType_ == MultiFluidData::EXTRAPOLATION){
    // *riemannupdate = 0.0;
    // *weight = 0.0;
    if(algorithmType_ == MultiFluidData::GHOSTFLUID_FOR_POOR)
      fprintf(stdout, "*** Error: not supposed to be here for GFMP\n");
    domain->storePreviousPrimitive(V,fluidId,X,*riemannupdate, *weight);
  }
  */
}
//------------------------------------------------------------------------------
template<int dim>
template<int dimLS>
void DistExactRiemannSolver<dim>::avoidNewPhaseCreation(DistSVec<double,dimLS> &Phi, DistSVec<double,dimLS> &Phin, DistLevelSetStructure *distLSS)
{
  if(algorithmType_ != MultiFluidData::GHOSTFLUID_FOR_POOR)
    domain->avoidNewPhaseCreation(Phi,Phin,*weight, distLSS, fluidIdToSet);
}
//------------------------------------------------------------------------------
