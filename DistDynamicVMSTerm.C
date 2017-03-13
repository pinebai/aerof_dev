#include <DistDynamicVMSTerm.h>

#include <FemEquationTerm.h>
#include <DistVector.h>
#include <DistGeoState.h>
#include <DistTimeState.h>
#include <DistBcData.h>
#include <Domain.h>
#include <DistEdgeGrad.h>
#include <DistNodalGrad.h>
#include <FluxFcn.h>
#include <RecFcnDesc.h>
#include <DynamicVMSTerm.h>

#include <Communicator.h>

#include <cstdio>

//------------------------------------------------------------------------

template<int dim>
DistDynamicVMSTerm<dim>::DistDynamicVMSTerm(VarFcn *vf, IoData &iod, Domain *dom) : domain(dom),varFcn(vf)

{

  int i;
  it0        = iod.restart.iteration;
  lastIt     = it0;
  lastConfig = -1;
  
  numLocSub  = domain->getNumLocSub();

  dvmst       = new DynamicVMSTerm(iod, varFcn);
  com = domain->getCommunicator();

  scopeWidth = iod.eqs.tc.les.dvms.agglomeration_width;
  scopeDepth1 = iod.eqs.tc.les.dvms.agglomeration_depth1;
  scopeDepth2 = iod.eqs.tc.les.dvms.agglomeration_depth2;

  // following takes care of the requirement that scopeWidth2 should be > scopeWidth1 //
  
  if(scopeDepth2 <= scopeDepth1){
    fprintf(stderr, "*** Warning: AgglomerationDepth2 = %d is invalid, Using AgglomerationDepth2 =  %d\n" ,
                         scopeDepth2, scopeDepth1+1);
    scopeDepth2 = scopeDepth1 + 1;
  }

   // to switch between D1-VMS-LES or D2-VMS-LES or D3-VMS-LES //
   
   // method = 0 is the solution by germano's identity 
   // method = 1 is direct least squares solution
   // method = 2 uses a special clipping procedure to avoid -ve Cs
  
  if(iod.eqs.tc.les.dvms.type == DynamicVMSData::D1VMSLES) method = 0;
  else if(iod.eqs.tc.les.dvms.type == DynamicVMSData::D2VMSLES) method = 1;
  else if(iod.eqs.tc.les.dvms.type == DynamicVMSData::D3VMSLES) method = 2;


  // allocation of all vectors involved //

  VBar = new DistSVec<double,dim>*[2];     // contains the primitive bar terms
  volRatio   = new DistSVec<double,1>*[2]; // contains the volume ratio ie.. vol(cell)/vol(macrocell)

  for (i = 0; i < 2; ++i) {
   VBar[i] = new DistSVec<double,dim>(domain->getNodeDistInfo());
   volRatio[i] = new  DistSVec<double,1>(domain->getNodeDistInfo());
  }

  RBar = new DistSVec<double,dim>(domain->getNodeDistInfo());    // coarse mesh residual
  dWBardt = new DistSVec<double,dim>(domain->getNodeDistInfo()); // coarse mesh temporal term
  dWdt = new DistSVec<double,dim>(domain->getNodeDistInfo());    // fine mesh temporal term
  MBar = new DistSVec<double,dim>(domain->getNodeDistInfo());    // coarse mesh model term
  M = new DistSVec<double,dim>(domain->getNodeDistInfo());       // fine mesh model term

  S = new DistSVec<double,dim>(domain->getNodeDistInfo());       // flux due to reynolds stress term

  // allocation of gradients //

  ngrad = new DistNodalGrad<dim, double>(iod, domain); // nodal gradients
  egrad = 0;                                           // edge based gradients 

  if (iod.schemes.ns.dissipation == SchemeData::SIXTH_ORDER)
    egrad = new DistEdgeGrad<dim>(iod, domain);

  masterFlag = new bool *[numLocSub];  // flag that controls inclusion of cells 
                                       // on the boundary of two subdomains to exactly one of them 

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; iSub++) {
    masterFlag[iSub] = VBar[0]->getMasterFlag(iSub);
  }

  macroCells = new DistMacroCellSet(domain, iod.eqs.fluidModel.gasModel.specificHeatRatio, masterFlag, scopeWidth, scopeDepth2); // creating macrocells

}

//------------------------------------------------------------------------

template<int dim>
DistDynamicVMSTerm<dim>::~DistDynamicVMSTerm()

{

  if (dvmst) delete dvmst;
  if (masterFlag) { delete masterFlag;}
  if (macroCells) delete macroCells;
  if (VBar) delete [] VBar;
  if (volRatio) delete [] volRatio;
  if (RBar) delete RBar;;
  if (S) delete S;
  if (dWdt) delete dWdt;
  if (dWBardt) delete dWBardt;
  if (M) delete M;
  if (MBar) delete MBar;
  if (ngrad) delete ngrad;
  if (egrad) delete egrad;

}

//------------------------------------------------------------------------

template<int dim>
void DistDynamicVMSTerm<dim>::compute(FluxFcn** fluxFcn, 
                                      RecFcn* recFcn, 
				      FemEquationTerm *fet, 
				      int config, 
				      DistVec<double> &ctrlVol,
				      DistBcData<dim> &bcData,
                                      DistGeoState &geoState,
                                      DistTimeState<dim> *timeState,
				      DistSVec<double,3> &X,
				      DistSVec<double,dim> &U,
                                      DistSVec<double,dim> &V, 
				      DistSVec<double,dim> &R, 
                                      int failsafe, int rshift) 
				  
{

  bool doInitialTasks = false;  // boolean that controls one-time operations 

  if ((lastConfig != config) || (lastIt == it0))
    doInitialTasks = true;

  // print volume ratio statistics //

  if(doInitialTasks) printDelRatios(ctrlVol);
 
  // Initialization of RBar, dWdt. dWBardt,M, MBar & S //

  *RBar = 0.0; 
  *dWdt = 0.0;
  *dWBardt = 0.0;
  *M = 0.0;
  *MBar = 0.0;
  *S = 0.0;    

  // compute Bar-terms required for computing dWBardt //
  
  timeState->computeBar(doInitialTasks, macroCells, geoState, scopeDepth1);

  // computing VBar, volRatio //

  domain->computeBarTerm(macroCells, doInitialTasks, ctrlVol, VBar, volRatio, X, V, scopeDepth1, scopeDepth2);

  // computing the time derivatives of conservative variables (i.e. dWdt & dWBardt) //
  
  timeState->get_dW_dt(doInitialTasks, geoState, ctrlVol, U, *(VBar[0]), 
                       *dWdt, *dWBardt, macroCells, volRatio, scopeDepth1);
  

  // computing coarse mesh residual (RBar) //

  if (egrad)
    egrad->compute(geoState.getConfig(), X);  // computing edge gradient (if specified)

  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)  
    ngrad->compute(geoState.getConfig(), X, ctrlVol, *(VBar[0])); // nodal gradient computation
  
  
  if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)
    ngrad->limit(recFcn, X, ctrlVol, *(VBar[0]));    // applying limiter and spheres/boxes

  domain->computeGalerkinBarTerm(doInitialTasks, fet, bcData, geoState, X, macroCells, *(VBar[0]), 
                                 *(volRatio[0]), *RBar, scopeDepth1, scopeDepth2); // residual for diffusion part 
                                                                                   // with coarse mesh shape function

  
  DistVec<double> *irey = timeState->getInvReynolds(); // for local low-mach preconditioning in viscous cases (ARL)
  domain->computeFiniteVolumeBarTerm(ctrlVol, *irey, fluxFcn, recFcn, bcData, geoState, X, macroCells, *(VBar[0]), 
                                     *(volRatio[0]), *ngrad, egrad, *RBar, scopeDepth1, scopeDepth2, failsafe, rshift);
                                                                                   // residual for convective part
                                                                                   // with coarse mesh characteristic function


  domain->computeDynamicVMSTerm(dvmst, macroCells, doInitialTasks, ctrlVol, VBar, volRatio, X, V, *S,
                                R, *RBar, *dWdt, *dWBardt, *M, *MBar, scopeDepth1, scopeDepth2, method);
	                                                                         // computing (Cs*Delta)^2 and PrT dynamically 							
  R += *S;  // adding the reynolds stress flux term to residual R

  if (lastConfig != config)
    lastConfig = config;

  if (lastIt == it0)
    lastIt = -1;

  irey = 0;

}

//------------------------------------------------------------------------

// Included (MB)
template<int dim>
void DistDynamicVMSTerm<dim>::computeDerivative(FluxFcn** fluxFcn, 
                                      RecFcn* recFcn, 
				      FemEquationTerm *fet, 
				      int config, 
				      DistVec<double> &ctrlVol,
				      DistBcData<dim> &bcData,
                                      DistGeoState &geoState,
                                      DistTimeState<dim> *timeState,
				      DistSVec<double,3> &X,
				      DistSVec<double,dim> &U,
                                      DistSVec<double,dim> &V, 
				      DistSVec<double,dim> &R, 
                                      int failsafe, int rshift) 
				  
{

  fprintf(stderr, "***** DistDynamicVMSTerm<dim>::computeDerivative is not implemented!\n");
  exit(1);

}

//------------------------------------------------------------------------------------------------------------

// This routine is used to retrive the value of Smag Coeff (Cs) for post processing //

template<int dim>
void DistDynamicVMSTerm<dim>::computeCs(FluxFcn** fluxFcn,
                                        RecFcn* recFcn,
                                        FemEquationTerm *fet,
                                        int config,
                                        DistVec<double> &ctrlVol,
                                        DistBcData<dim> &bcData,
                                        DistGeoState &geoState,
                                        DistTimeState<dim> *timeState,
                                        DistSVec<double,3> &X,
                                        DistSVec<double,dim> &U,
                                        DistSVec<double,dim> &V,
                                        DistSVec<double,dim> &R,
                                        DistVec<double> *Cs, 
                                        int failsafe, int rshift)

{

    *RBar = 0.0;
    *dWdt = 0.0;
    *dWBardt = 0.0;
    *M = 0.0;
    *MBar = 0.0;
    *S = 0.0;
    
    bool doInitialTasks = true;

    if((lastConfig == config) || (lastIt != it0)) { 
        timeState->computeBar(doInitialTasks, macroCells, geoState, scopeDepth1);
    }

    domain->computeBarTerm(macroCells, doInitialTasks, ctrlVol, VBar, volRatio, X, V, scopeDepth1, scopeDepth2);

    if((lastConfig == config) || (lastIt != it0)) {
      timeState->get_dW_dt(doInitialTasks, geoState, ctrlVol, U, *(VBar[0]), 
                           *dWdt, *dWBardt, macroCells, volRatio, scopeDepth1);
    }

    if (egrad)
      egrad->compute(geoState.getConfig(), X);

    if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)
      ngrad->compute(geoState.getConfig(), X, ctrlVol, *(VBar[0]));

    if (dynamic_cast<RecFcnConstant<dim> *>(recFcn) == 0)
      ngrad->limit(recFcn, X, ctrlVol, *(VBar[0]));

    domain->computeGalerkinBarTerm(doInitialTasks, fet, bcData, geoState, X, macroCells, *(VBar[0]),
                                   *(volRatio[0]), *RBar, scopeDepth1, scopeDepth2);

    DistVec<double> *irey = timeState->getInvReynolds(); // for local low-mach preconditioning in viscous cases (ARL)
    domain->computeFiniteVolumeBarTerm(ctrlVol, *irey, fluxFcn, recFcn, bcData, geoState, X, macroCells, *(VBar[0]), 
                                       *(volRatio[0]), *ngrad, egrad, *RBar, scopeDepth1, scopeDepth2, failsafe, rshift);

    domain->computeDynamicVMSTerm(dvmst, macroCells, doInitialTasks, ctrlVol, VBar, volRatio, X, V, *S,
                                  R, *RBar, *dWdt, *dWBardt, *M, *MBar, scopeDepth1, scopeDepth2, method, Cs);

    if (lastConfig != config)
      lastConfig = config;

    if (lastIt == it0)
      lastIt = -1;

    irey = 0;

}

//------------------------------------------------------------------------------------------------------------

template<int dim>
void DistDynamicVMSTerm<dim>::computeMutOMu(DistVec<double> &ctrlVol,
                                            DistSVec<double,3> &X,
                                            DistSVec<double,dim> &V,
                                            DistVec<double> &Cs,
                                            DistVec<double> &mutOmu)
{

  bool doInitialTasks = true;

  domain->computeBarTerm(macroCells, doInitialTasks, ctrlVol, VBar, volRatio, X, V, scopeDepth1, scopeDepth2);

  domain->computeMutOMuDynamicVMS(dvmst, ctrlVol, *(VBar[0]), X, V,  Cs, mutOmu);

}

//------------------------------------------------------------------------

// This routine prints out the ratios of volumes between different levels of agglomeration //

template<int dim>
void DistDynamicVMSTerm<dim>:: printDelRatios(DistVec<double> &ctrlVol)
{

  int numavgs = scopeDepth2*(scopeDepth2+1)/2;
  double *rmax = new double[scopeDepth2];
  double *rmin = new double[scopeDepth2];
  double *rsum = new double[scopeDepth2];
  int *numNodes = new int[scopeDepth2];
  double *ravg = new double[scopeDepth2];
  double *rravg = new double[numavgs];
  double *rrsum = new double[numavgs];

  domain->computeDelRatios(macroCells, ctrlVol, scopeDepth2, rmax, rmin, rsum, rrsum, numNodes);

  com->globalMin(1, rmax);
  com->globalMax(1, rmin);
  com->globalSum(1, rsum);
  com->globalSum(1, rrsum);
  com->globalSum(1, numNodes);

  int numCPU = com->size();

  for (int i=0; i<scopeDepth2; ++i) {
    ravg[i] = rsum[i]/numNodes[i];
  }

  int k = 0;
  for (int i=0; i<scopeDepth2-1; ++i) {
    for (int j=i+1; j<scopeDepth2; ++j) {
      rravg[k] = rrsum[k]/numNodes[j-1];
      k +=1;
    }
  }

  com->printf(2, "              Volume Ratio Statistics            \n");     
  com->printf(2, "\n");     
  com->printf(2, " Agg. Depth     :   MaxRatio     AvgRatio    MinRatio\n");     
  com->printf(2, "\n");     
               
  for (int i = 0; i<scopeDepth2; ++i) {
    com->printf(2, "      %d         :   %4.3e    %4.3e   %4.3e\n", i+1, rmax[i]/100, ravg[i]/100, rmin[i]/100);
  }

  com->printf(2, "\n");     
  com->printf(2, " Level1  Level2   :   MaxRatio     AvgRatio    MinRatio\n");
  com->printf(2, " Depth   Depth \n");
  com->printf(2, "\n");

  
  k = 0;
  for (int i = 0; i<scopeDepth2-1; ++i) {
    for (int j = i+1; j<scopeDepth2; ++j) {
      com->printf(2, "   %d        %d     :   %4.3e    %4.3e   %4.3e\n", i+1, j+1, rmax[j]/rmax[i], 
                                                                         rravg[k]/100, rmin[j]/rmin[i]);
      k = k+1;
    }
  }
  com->printf(2, "\n");

  // Delete all the allocated values //

  delete [] rmax; delete [] rmin; delete [] rsum; 
  delete [] numNodes; delete [] ravg; delete [] rravg; delete [] rrsum;

}

//--------------------------------------------------------------------------------------------------------
