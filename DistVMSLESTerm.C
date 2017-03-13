#include <DistVMSLESTerm.h>

#include <VMSLESTerm.h>
#include <DistVector.h>
#include <Domain.h>

#include <cstdio>

//------------------------------------------------------------------------

template<int dim>
DistVMSLESTerm<dim>::DistVMSLESTerm(VarFcn *vf, IoData &iod, Domain *dom) : domain(dom), varFcn(vf)

{

  // initialization of variables that matter in computing the Bar terms only once //

  it0        = iod.restart.iteration;
  lastIt     = it0;
  lastConfig = -1;
  
  numLocSub  = domain->getNumLocSub();
  scopeWidth = iod.eqs.tc.les.vms.agglomeration_width;
  scopeDepth = iod.eqs.tc.les.vms.agglomeration_depth;

  VtBar      = new DistSVec<double,dim>(domain->getNodeDistInfo());
  volRatio   = new DistSVec<double,1>(domain->getNodeDistInfo());
  masterFlag = new bool *[numLocSub];
 
  // initialization of VBar and volRatio //

  *VtBar = 0.0;
  *volRatio = 0.0;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; iSub++) {
    masterFlag[iSub] = VtBar->getMasterFlag(iSub);
  }

  // creation of macrocells //

  macroCells = new DistMacroCellSet(domain, iod.eqs.fluidModel.gasModel.specificHeatRatio, masterFlag, scopeWidth, scopeDepth);
  vmst       = new VMSLESTerm(iod, varFcn);

}

//------------------------------------------------------------------------

template<int dim>
DistVMSLESTerm<dim>::~DistVMSLESTerm()

{

  if (macroCells) delete macroCells;
  if (vmst) delete vmst;

  if (masterFlag) { delete masterFlag; }

  if (VtBar) delete VtBar;
  if (volRatio) delete volRatio;

}

//------------------------------------------------------------------------

template<int dim>
void DistVMSLESTerm<dim>::compute(int config, DistVec<double> &ctrlVol,
				  DistSVec<double,3> &X,
				  DistSVec<double,dim> &V,
				  DistSVec<double,dim> &R)

{

  bool doInitialTasks = false;

  if ((lastConfig != config) || (lastIt == it0))
    doInitialTasks = true;

  domain->computeVMSLESTerm(vmst, macroCells, doInitialTasks,
			    ctrlVol, *VtBar, *volRatio, X, V, R, scopeDepth);

  if (lastConfig != config)
    lastConfig = config;

  if (lastIt == it0)
    lastIt = -1;

}

//------------------------------------------------------------------------

template<int dim>
void DistVMSLESTerm<dim>::computeMutOMu(DistVec<double> &ctrlVol,
                                        DistSVec<double,3> &X,
                                        DistSVec<double,dim> &V,
                                        DistVec<double> &mutOmu)
{

  bool doInitialTasks = true;

  domain->computeMutOMuVMS(vmst, macroCells, ctrlVol, 
                           doInitialTasks, *VtBar, *volRatio,
                           X, V, scopeDepth, mutOmu);

}

//------------------------------------------------------------------------
// Included (MB)
template<int dim>
void DistVMSLESTerm<dim>::computeDerivative(int config, DistVec<double> &ctrlVol,
				  DistSVec<double,3> &X,
				  DistSVec<double,dim> &V,
				  DistSVec<double,dim> &R)

{

  fprintf(stderr, "***** DistVMSLESTerm<dim>::computeDerivative is not implemented!\n");
  exit(1);

}

//------------------------------------------------------------------------
