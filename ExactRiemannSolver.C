#include <ExactRiemannSolver.h>

#include <IoData.h>
#include <Vector.h>
#include <Vector3D.h>
#include <LocalRiemannDesc.h>

#include <cmath>


//------------------------------------------------------------------------------
template<int dim>
ExactRiemannSolver<dim>::ExactRiemannSolver(IoData &iod, SVec<double,dim> &_rupdate,
                                            Vec<double> &_weight, SVec<double,dim-2> &_interfacialWi,
                                            SVec<double,dim-2> &_interfacialWj, VarFcn *vf,
                                            SparseGridCluster *sgCluster,
                                            Vec<int>& fluidIdToSet):
                                            rupdate(_rupdate), weight(_weight),
                                            interfacialWi(_interfacialWi),
                                            interfacialWj(_interfacialWj),
                                            fluidIdToSet(fluidIdToSet)
{

  iteration = -1;
  lriemann = 0;
  numLriemann = 0;
  fsiRiemann = 0;

// FSI Riemann problem
  if(iod.problem.framework==ProblemData::EMBEDDED ||
     iod.problem.framework==ProblemData::EMBEDDEDALE ||
     iod.bc.wall.reconstruction==BcsWallData::EXACT_RIEMANN ||
     iod.oneDimensionalInfo.problemMode == OneDimensionalInfo::FSI) {
    fsiRiemann = new LocalRiemannFluidStructure<dim>();
    dynamic_cast<LocalRiemannFluidStructure<dim> *>(fsiRiemann)->setStabilAlpha(iod.embed.stabil_alpha);

    if (iod.embed.prec == EmbeddedFramework::PRECONDITIONED)
      dynamic_cast<LocalRiemannFluidStructure<dim> *>(fsiRiemann)->setPreconditioner(iod.prec.mach);
	 
	 if (iod.eqs.type == EquationsData::NAVIER_STOKES) {
		 if(iod.embed.viscousboundarycondition == EmbeddedFramework::STRONG){
			 dynamic_cast<LocalRiemannFluidStructure<dim> *>(fsiRiemann)->setViscousSwitch(1.0);
		 } else {
			 dynamic_cast<LocalRiemannFluidStructure<dim> *>(fsiRiemann)->setViscousSwitch(0.0);
		 }  
	 }
  }

    // Riemann solver for embedded constrains todo: Daniel, I allocate these riemann solvers in all cases, but they are only needed when special constrains are there
    actuatorDiskRiemann = new LocalRiemannActuatorDisk<dim>();
    //dynamic_cast<LocalRiemannActuatorDisk<dim> *>(actuatorDiskRiemann);
    symmetryplaneRiemann = new LocalRiemannFluidStructure<dim>();
    //dynamic_cast<LocalRiemannFluidStructure<dim> *>(symmetryplaneRiemann)->setStabilAlpha(0.0);


    for (int i = 0; i < 10; ++i) {
    for (int j = 0; j < 10; ++j) {
      levelSetMap[i][j] = -1; //KW: for debugging
      levelSetSign[i][j]=1.0;
    }
  }

// Multiphase Riemann problem
  if(iod.eqs.numPhase > 1){

    numLriemann = (iod.eqs.numPhase-1)*(iod.eqs.numPhase)/2;
    lriemann = new LocalRiemann*[numLriemann];
    memset(lriemann, 0, sizeof(LocalRiemann*)*numLriemann) ;
    int iRiemann = 0;

    for (int iPhase = 0; iPhase < iod.eqs.numPhase; ++iPhase) {
      for (int jPhase = iPhase + 1; jPhase < iod.eqs.numPhase; ++jPhase) {
    
	int fluidI = iPhase;
	int fluidJ = jPhase;
	int fluid1,fluid2;
	
	for (int k = 0; k < 1; ++k,++iRiemann) {
	  if (k == 0) {
	    fluid1 = fluidI;
	    fluid2 = fluidJ;
	  } else {
	    fluid1 = fluidJ;
	    fluid2 = fluidI;
	  }
	  levelSetMap[fluid1][fluid2] = iRiemann;
	  levelSetMap[fluid2][fluid1] = iRiemann;
	  
	  map<int, FluidModelData *>::iterator it1 = iod.eqs.fluidModelMap.dataMap.find(fluid1);
	  if(it1 == iod.eqs.fluidModelMap.dataMap.end()){
	    fprintf(stderr, "*** Error: no FluidModel[%d] was specified\n", fluid1);
	    exit(1);
	  }
	  map<int, FluidModelData *>::iterator it2 = iod.eqs.fluidModelMap.dataMap.find(fluid2);
	  if(it2 == iod.eqs.fluidModelMap.dataMap.end()){
	    fprintf(stderr, "*** Error: no FluidModel[%d] was specified\n", fluid2);
	    exit(1);
	  }
	  bool fluid1IsAGas = (it1->second->fluid  == FluidModelData::PERFECT_GAS || it1->second->fluid  == FluidModelData::STIFFENED_GAS);
	  bool fluid2IsAGas = (it2->second->fluid  == FluidModelData::PERFECT_GAS || it2->second->fluid  == FluidModelData::STIFFENED_GAS);
	  if(iod.mf.method == MultiFluidData::GHOSTFLUID_FOR_POOR){
	    if(fluid1IsAGas && fluid2IsAGas)
	      lriemann[iRiemann] = new LocalRiemannGfmpGasGas(vf,fluid1,fluid2);
	    else if(it1->second->fluid  == FluidModelData::LIQUID &&
		    it2->second->fluid == FluidModelData::LIQUID)
	      lriemann[iRiemann] = new LocalRiemannGfmpTaitTait(vf,fluid1,fluid2);
	    else if(it1->second->fluid  == FluidModelData::JWL &&
		    it2->second->fluid == FluidModelData::JWL)
	      lriemann[iRiemann] = new LocalRiemannGfmpJWLJWL(vf,fluid1,fluid2);
	    else if(fluid1IsAGas &&
		    it2->second->fluid == FluidModelData::JWL)
	      lriemann[iRiemann] = new LocalRiemannGfmpGasJWL(vf,fluid1,fluid2);
	/*    else{
	      // fprintf(stdout, "*** Warning: No GFMP possible between fluid models %i and %i\n",fluid1, fluid2);
	    }*/
	  }else if(iod.mf.method == MultiFluidData::GHOSTFLUID_WITH_RIEMANN){

	    
	    if (iod.mf.prec == MultiFluidData::PRECONDITIONED) {

	      lriemann[iRiemann] = new LocalRiemannLowMach(vf,fluid1,fluid2,
							   iod.prec.mach, dim);
	    } else {

	      if(fluid1IsAGas && fluid2IsAGas){
		lriemann[iRiemann] = new LocalRiemannGfmparGasGas(vf,fluid1,fluid2, iod.mf.typePhaseChange,
                                                                  iod.mf.riemannEps, iod.mf.riemannMaxIts);
	      }
	      else if(it1->second->fluid  == FluidModelData::LIQUID &&
		      fluid2IsAGas){
		levelSetSign[fluid1][fluid2] = levelSetSign[fluid2][fluid1] = -1.0;
		lriemann[iRiemann] = new LocalRiemannGfmparGasTait(vf,fluid2,fluid1, iod.mf.typePhaseChange);
	      }
	      else if(it1->second->fluid  == FluidModelData::LIQUID &&
		      it2->second->fluid == FluidModelData::LIQUID){
		lriemann[iRiemann] = new LocalRiemannGfmparTaitTait(vf,fluid1,fluid2, iod.mf.typePhaseChange,
                                                                    iod.mf.riemannEps, iod.mf.riemannMaxIts);
	      }
	      else if(fluid1IsAGas && 
		      it2->second->fluid == FluidModelData::LIQUID){
		lriemann[iRiemann] = new LocalRiemannGfmparGasTait(vf,fluid1,fluid2, iod.mf.typePhaseChange);
	      }
	      else if(fluid1IsAGas && 
		      it2->second->fluid == FluidModelData::JWL){
		lriemann[iRiemann] = new LocalRiemannGfmparGasJWL(vf,fluid1,fluid2,sgCluster,iod.mf.riemannComputation,
								  iod.mf.jwlRelaxationFactor,
								  iod.ref.rv.density,iod.ref.rv.entropy,iod.mf.typePhaseChange);
	      }
	      else if(it1->second->fluid  == FluidModelData::JWL &&
		      it2->second->fluid == FluidModelData::JWL){
		lriemann[iRiemann] = new LocalRiemannGfmparJWLJWL(vf,fluid1,fluid2, iod.mf.typePhaseChange);
	      }
	      else if(it1->second->fluid  == FluidModelData::LIQUID &&
		      it2->second->fluid == FluidModelData::JWL){
		lriemann[iRiemann] = new LocalRiemannGfmparTaitJWL(vf,fluid1,fluid2,sgCluster,iod.mf.riemannComputation,
								   iod.mf.jwlRelaxationFactor, iod.ref.rv.density,iod.ref.rv.entropy,iod.mf.typePhaseChange);
	      }/* else if(it1->second->fluid  == FluidModelData::JWL &&
		  it2->second->fluid == FluidModelData::LIQUID){
		  lriemann[iRiemann] = new LocalRiemannGfmparTaitJWL(vf,fluid2,fluid1,sgCluster,iod.mf.riemannComputation, iod.mf.typePhaseChange,-1.0);
		  }*//* else{
	      
		     // fprintf(stdout, "*** Warning: No GFMP possible between fluid models %i and %i\n",fluid1, fluid2);
		     }*/
	    }
	  } 
	}
      }
    }
  }

}

//------------------------------------------------------------------------------
template<int dim>
ExactRiemannSolver<dim>::~ExactRiemannSolver() 
{

  for(int iLriemann=0; iLriemann<numLriemann; iLriemann++) {
    if (lriemann[iLriemann])
      delete lriemann[iLriemann];
  }
  delete [] lriemann;
  delete fsiRiemann;

}

//------------------------------------------------------------------------------
template<int dim>
int ExactRiemannSolver<dim>::updatePhaseChange(SVec<double,dim> &V, Vec<int> &fluidId,
                                                Vec<int> &fluidIdn,HigherOrderMultiFluid* higherOrderMF)
{

  for(int i=0; i<V.size(); i++){ 
//    if(fluidId[i]==2) {fprintf(stderr,"OK. you are 2, you were %d.\n", fluidIdn[i]);}
//    if(fluidId[i]==2 && fluidIdn[i]==1) {fprintf(stderr,"I caught you! weight = %e.\n", weight[i]);}
    if (lriemann[0]->updatePhaseChange(V[i],fluidId[i],fluidIdn[i],rupdate[i],weight[i],false) != 0) {
      return i;
    }
    if (higherOrderMF && fluidId[i] != fluidIdn[i]) {
      higherOrderMF->setLastPhaseChangeValue<dim>(i, V[i]);
    }
    // lriemann[0] can be used for all interfaces, because  this routine does
    // not need to consider which interface has traversed that node (this
    // was done previously when computing the riemann problem)
  }

  return -1;
}
//------------------------------------------------------------------------------
template<int dim>
int ExactRiemannSolver<dim>::fluid1(int IDi, int IDj)
{
  int riemannId = levelSetMap[IDi][IDj];
  if(riemannId<0) {
    fprintf(stderr,"ERROR: There is no Riemann solver for IDi = %d and IDj = %d!\n", IDi, IDj);
    exit(-1);
  }
  return (levelSetSign[IDi][IDj] > 0) ? lriemann[riemannId]->fluid1 : lriemann[riemannId]->fluid2;
}
//------------------------------------------------------------------------------
template<int dim>
int ExactRiemannSolver<dim>::fluid2(int IDi, int IDj)
{
  int riemannId = levelSetMap[IDi][IDj];
  if(riemannId<0) {
    fprintf(stderr,"ERROR: There is no Riemann solver for IDi = %d and IDj = %d!\n", IDi, IDj);
    exit(-1);
  }
  return (levelSetSign[IDi][IDj] > 0) ? lriemann[riemannId]->fluid2 : lriemann[riemannId]->fluid1;
}
//------------------------------------------------------------------------------
template<int dim>
int ExactRiemannSolver<dim>::computeRiemannSolution(double *Vi, double *Vj,
						    int IDi, int IDj, double *nphi, VarFcn *vf,
						    double *Wi, double *Wj, int i, int j, int edgeNum,
						    double dx[3], int lsdim,
						    bool isHigherOrder)
{

  //fprintf(stdout, "Debug: calling computeRiemannSolution with IDi = %d - IDj = %d for LocalRiemann[%d]\n", IDi, IDj, IDi+IDj-1);
  int riemannId = levelSetMap[IDi][IDj];
  if(riemannId<0) {
    fprintf(stderr,"ERROR: There is no Riemann solver for IDi = %d and IDj = %d!\n", IDi, IDj);
    exit(-1);
  }
  double lssign = levelSetSign[IDi][IDj];
  double nphilss[3];
  for (int k=0; k < 3; ++k) {
    nphilss[k] = nphi[k]*lssign;
  }
  fluidIdToSet[i] = fluidIdToSet[j] = lsdim;
  return lriemann[riemannId]->computeRiemannSolution(Vi,Vj,IDi,IDj,nphilss,interfacialWi[edgeNum],interfacialWj[edgeNum],
					      Wi,Wj,rupdate[i],rupdate[j],weight[i],weight[j],
					      dx,iteration,isHigherOrder);


}
//------------------------------------------------------------------------------
template<int dim>
void ExactRiemannSolver<dim>::computeRiemannJacobian(double *Vi, double *Vj,
						     int IDi, int IDj, double *nphi, VarFcn *vf,
						     double *Wi, double *Wj,
						     int i, int j, int edgeNum, double dx[3],
						     double* dWidUi,double*  dWidUj,double* dWjdUi,double*  dWjdUj) {

  int riemannId = levelSetMap[IDi][IDj];
  if(riemannId<0) {
    fprintf(stderr,"ERROR: There is no Riemann solver for IDi = %d and IDj = %d!\n", IDi, IDj);
    exit(-1);
  }
  double lssign = levelSetSign[IDi][IDj];
  double nphilss[3];
  for (int k=0; k < 3; ++k) {
    nphilss[k] = nphi[k]*lssign;
  }
  lriemann[riemannId]->computeRiemannJacobian(Vi,Vj,IDi,IDj,nphilss,
          Wi,Wj,
          dx,iteration, dWidUi, dWidUj,dWjdUi, dWjdUj);
}

//------------------------------------------------------------------------------
template<int dim>
int ExactRiemannSolver<dim>::computeFSIRiemannSolution(double *Vi, double *Vstar,
      double *nphi, VarFcn *vf, double *Wstar, int nodej, int Id)

{
  return fsiRiemann->computeRiemannSolution(Vi,Vstar,nphi,vf,
         Wstar,rupdate[nodej],weight[nodej],iteration, Id);
}
//------------------------------------------------------------------------------
template<int dim>
void ExactRiemannSolver<dim>::computeFSIRiemannJacobian(double *Vi, double *Vstar,
      double *nphi, VarFcn *vf, double *Wstar, int nodej, double* dWdW,int Id)

{
  fsiRiemann->computeRiemannJacobian(Vi,Vstar,nphi,vf,
         Wstar,rupdate[nodej],weight[nodej],iteration, dWdW,Id);
}
//------------------------------------------------------------------------------
template<int dim>
void ExactRiemannSolver<dim>::computeFSIRiemannderivative(double *Vi, double *Vstar,
      double *nphi, VarFcn *vf, double *Wstar, int nodej, double* dWstardn, int Id)

{
  fsiRiemann->computeRiemannderivative(Vi, Vstar, nphi, vf, Wstar, dWstardn, Id);
}
//------------------------------------------------------------------------------
template<int dim>
void ExactRiemannSolver<dim>::computeFSIRiemannSolution(int tag, double *Vi, double *Vstar,
      double *nphi, VarFcn *vf, double *Wstar, int nodej)

{
  // Adam 2010.08.18
  // This function doesn't seem to be used anymore.
  // To be removed in a couple of months
  fprintf(stderr,"Oh Sorry ! Please uncomment the function (ExactRiemannSolver.C:159). I thought it wasn't needed anymore\n");
  exit(-1);
  /*
  fsiRiemann->computeRiemannSolution(tag, Vi,Vstar,nphi,vf,
         Wstar,rupdate[nodej],weight[nodej],iteration);
  */
}
//------------------------------------------------------------------------------
template<int dim>
void ExactRiemannSolver<dim>::reset(int it)
{

  iteration = it;
  if(iteration==1){
    rupdate = 0.0;
    weight  = 0.0;
    fluidIdToSet = -1;
  }

}
//------------------------------------------------------------------------------
template<int dim>
void ExactRiemannSolver<dim>::resetInterfacialW(int edgeNum)
{

  //dim-2 since interfacialW store only rho, u, p
  for(int idim=0; idim<dim-2; idim++){
    interfacialWi[edgeNum][idim] = 0.0;
    interfacialWj[edgeNum][idim] = 0.0;
  }

}
//------------------------------------------------------------------------------
template <int dim>
int ExactRiemannSolver<dim>::getRiemannSolverId(int i, int j) const {
  
  if (i == j) {
    fprintf(stderr, "Error: getInterfaceId called with fluididi = fluididj!\n");
    return 0;
  }

  return levelSetMap[i][j];
}
//------------------------------------------------------------------------------
template<int dim>
int ExactRiemannSolver<dim>::computeActuatorDiskRiemannSolution(double *Vi, double *Vj, double *Vstar, double dp,double *n_s,
                                                                double *n_f, VarFcn *vf,
                                                                double *Wi, double *Wj,int Id)

{   return actuatorDiskRiemann->computeRiemannSolution(Vi, Vj, Vstar,dp, n_s, n_f, vf, iteration, Wi, Wj,Id);

}

template<int dim>
void ExactRiemannSolver<dim>::computeActuatorDiskSourceTerm(double *Vi, double *Vj,double dp,double *n_s,
                                                                double *n_f, VarFcn *vf,
                                                                double *flux,bool method, int Id)

{    actuatorDiskRiemann->computeSourceTerm(Vi, Vj,dp, n_s, n_f, vf, flux,method, Id);

}

//------------------------------------------------------------------------------
template<int dim>
int ExactRiemannSolver<dim>::computeSymmetryPlaneRiemannSolution(double *Vi, double *Vstar,
                                                                 double *nphi, VarFcn *vf, double *Wstar, int nodej, int Id)

{
    return symmetryplaneRiemann->computeRiemannSolution(Vi,Vstar,nphi,vf,
                                                        Wstar,rupdate[nodej],weight[nodej],iteration, Id);
}