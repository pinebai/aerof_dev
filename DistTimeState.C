#include <DistTimeState.h>

#include <IoData.h>
#include <VarFcn.h>
#include <RefVal.h>
#include <TimeState.h>
#include <Domain.h>
#include <DistGeoState.h>
#include <DistMatrix.h>
#include <DistVectorOp.h>
#include <cmath>
#include <fstream>
#include <OneDimensionalSolver.h>
#include <ErrorHandler.h>
#include <ExactSolution.h>

using namespace std;
//------------------------------------------------------------------------------
template<int dim>
DistTimeState<dim>::DistTimeState(IoData &ioData, SpaceOperator<dim> *spo, VarFcn *vf,
				  Domain *dom, DistSVec<double,dim> *v) 
  : varFcn(vf), domain(dom) {

  initialize(ioData,spo,vf,dom,v,dom->getNodeDistInfo());
  
  refTime = ioData.ref.rv.time;

}

template<int dim>
DistTimeState<dim>::DistTimeState(IoData &ioData, SpaceOperator<dim> *spo, VarFcn *vf,
				  Domain *dom, DistInfo& dI,DistSVec<double,dim> *v) 
  : varFcn(vf), domain(dom) {

  initialize(ioData,spo,vf,dom,v,dI);
  
  refTime = ioData.ref.rv.time;

}

template<int dim>
void DistTimeState<dim>::initialize(IoData &ioData, SpaceOperator<dim> *spo, VarFcn *vf,
				  Domain *dom, DistSVec<double,dim> *v, DistInfo& dI) 
{
  locAlloc = true;

  checkForRapidlyChangingValues = true;

  if (v) V = v->alias();
  else V = new DistSVec<double,dim>(dI);

  numLocSub = domain->getNumLocSub();

  data = new TimeData(ioData);  

  dt  = new DistVec<double>(dI);
  idti = new DistVec<double>(dI);
  idtv = new DistVec<double>(dI);
  dtau  = new DistVec<double>(dI);
  irey  = new DistVec<double>(dI);
  viscousCst = ioData.ts.viscousCst;
  Un  = new DistSVec<double,dim>(dI);
  Vn = new DistSVec<double,dim>(dI);
  firstOrderNodes = new DistVec<int>(dI);

  *firstOrderNodes = 0;

  if (data->use_nm1)
    Unm1 = new DistSVec<double,dim>(dI);
  else
    Unm1 = Un->alias();

  if (data->use_nm2)
    Unm2 = new DistSVec<double,dim>(dI);
  else
    Unm2 = Unm1->alias();

  if (ioData.eqs.tc.les.type == LESModelData::DYNAMICVMS) {
    QBar = new DistSVec<double,dim>(dI);
    VnBar = new DistSVec<double,dim>(dI);
    UnBar = new DistSVec<double,dim>(dI);
    if (data->use_nm1)
      Unm1Bar = new DistSVec<double,dim>(dI);
    else
      Unm1Bar = UnBar->alias();
    if (data->use_nm2)
      Unm2Bar = new DistSVec<double,dim>(dI);
    else
      Unm2Bar = Unm1Bar->alias();
  }
  else {
    QBar = 0;
    VnBar = 0;
    //Vn = 0;
    UnBar = 0;
    Unm1Bar = 0;
    Unm2Bar = 0;
  }
                                                                                                                          
  if (data->typeIntegrator == ImplicitData::CRANK_NICOLSON)
    Rn = new DistSVec<double,dim>(dI);
  else
    Rn = Un->alias();

  errorEstiNorm = 0.0;
  dtMin = 1.e-10;

  gam = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
  pstiff = ioData.eqs.fluidModel.gasModel.pressureConstant;

  if (spo)
    fet = spo->getFemEquationTerm();
  else
    fet = NULL;

  //preconditioner setup
  tprec.setup(ioData);
  sprec.setup(ioData);

  subTimeState = new TimeState<dim>*[numLocSub];
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) 
    subTimeState[iSub] = 0;

// Included (MB)
  if (ioData.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_ ||
      ioData.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_ ||
	  ioData.problem.alltype == ProblemData::_ROM_SHAPE_OPTIMIZATION_ ||
	  ioData.problem.alltype == ProblemData::_SENSITIVITY_ANALYSIS_ ) {
    dIdti = new DistVec<double>(dI);
    dIdtv = new DistVec<double>(dI);
    dIrey = new DistVec<double>(dI);
    *dIdti = 0.0; 
    *dIdtv = 0.0; 
    *dIrey = 0.0;
  }
  else {
    dIdti = 0; 
    dIdtv = 0; 
    dIrey = 0;
        //each volume (volIt->first) is setup using
  }

  double* outputNewtonTag = domain->getNewtonTag();
  *outputNewtonTag = data->getNewtonTag();

  int* outputNewtonStateStep = domain->getNewtonStateStep();
  *outputNewtonStateStep = data->getNewtonStateStep();

  int* outputNewtonResidualStep = domain->getNewtonResidualStep(); 
  *outputNewtonResidualStep = data->getNewtonResidualStep();

  int* outputKrylovStep = domain->getKrylovStep();
  *outputKrylovStep = data->getKrylovStep();

  isGFMPAR = (ioData.eqs.numPhase > 1 &&
              ioData.mf.method == MultiFluidData::GHOSTFLUID_WITH_RIEMANN);
  
  fvmers_3pbdf = ioData.ts.implicit.fvmers_3pbdf;

  mf_phase_change_type = (ioData.mf.typePhaseChange != MultiFluidData::EXTRAPOLATION);

  *dtau = 1.0;
  dt_coeff = 1.0;
  dt_coeff_count = 0;
  allowcflstop = true;
  allowdtstop = true;

  *irey = 0.0;

  checkForRapidlyChangingPressure = ioData.ts.rapidPressureThreshold;
  checkForRapidlyChangingDensity = ioData.ts.rapidDensityThreshold;

  errorHandler = dom->getErrorHandler();

  refTime = ioData.ref.rv.time;
}

//------------------------------------------------------------------------------

template<int dim>
DistTimeState<dim>::DistTimeState(const DistTimeState<dim> &ts, bool typeAlloc, IoData &ioData) 
{

  locAlloc = typeAlloc;

  gam = ioData.eqs.fluidModel.gasModel.specificHeatRatio;
  pstiff = ioData.eqs.fluidModel.gasModel.pressureConstant;

  //preconditioner setup
  tprec.setup(ioData);
  sprec.setup(ioData);

  varFcn = ts.varFcn;
  fet = ts.fet;
  V = ts.V;

  numLocSub = ts.numLocSub;

  data = ts.data;

  dt  = ts.dt;
  idti = ts.idti;
  idtv = ts.idtv;
  dtau = ts.dtau;
  irey = ts.irey;
  viscousCst = ts.viscousCst;
  Un  = ts.Un;
  Unm1 = ts.Unm1;
  Unm2 = ts.Unm2;
  Rn = ts.Rn;

  QBar = ts.QBar;
  VnBar = ts.VnBar;
  Vn = ts.Vn;
  UnBar = ts.UnBar;
  Unm1Bar = ts.Unm1Bar;
  Unm2Bar = ts.Unm2Bar;

  domain = ts.domain;

  subTimeState = ts.subTimeState;

  refTime = ioData.ref.rv.time;

}

//------------------------------------------------------------------------------

template<int dim>
DistTimeState<dim>::~DistTimeState()
{

  if (!locAlloc) return;

  if (data) delete data;
  if (V) delete V;
  if (Un) delete Un;
  if (Unm1) delete Unm1;
  if (Unm2) delete Unm2;
  if (Rn) delete Rn;
  if (QBar) delete QBar;
  if (Vn) delete Vn;
  if (VnBar) delete VnBar;
  if (UnBar) delete UnBar;
  if (Unm1Bar) delete Unm1Bar;
  if (Unm2Bar) delete Unm2Bar;
                                                                                                                          
  delete dt;
  delete idti;
  delete idtv; 
  delete dtau;
  delete irey;
                                                                                                                          
  if (subTimeState) {
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub)
      if (subTimeState[iSub]) 
	delete subTimeState[iSub]; 
    
    delete [] subTimeState;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::attachHH(DistVec<double>& hh) {

  hhn = new DistVec<double>(hh);
  hhnm1 = new DistVec<double>(hh);
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    if (subTimeState[iSub]) 
      subTimeState[iSub]->attachHH(&(*hhn)(iSub),&(*hhnm1)(iSub));
  }
}

template<int dim>
void DistTimeState<dim>::setup(const char *name, DistSVec<double,3> &X,
                               DistSVec<double,dim> &Ufar,
                               DistSVec<double,dim> &U, IoData &iod,
                               DistVec<int> *point_based_id)
{
  *Un = Ufar;
 
  // convention:
  // initialization : U = Ufar
  // first,  setup U for volume with specific volume Ids
  // optional, setup U with one dimensional solution
  // second, setup U for multiphase geometric conditions (planes, then spheres)
  // third,  setup U for embedded structures (points)
  // NOTE: each new setup overwrites the previous ones.
  // CHANGED by Alexander the Great.  Now the one dimensional solution is setup first, followed by the point ic's

  setupUOneDimensionalSolution(iod,X);
  if(point_based_id)
    setupUFluidIdInitialConditions(iod, *point_based_id);

  setupUVolumesInitialConditions(iod);
  setupUMultiFluidInitialConditions(iod,X);
  setupUExactSolutionInitialConditions(iod,X);

  if (name[0] != 0) {
    domain->readVectorFromFile(name, 0, 0, *Un);
    if (data->use_nm1) {
      data->exist_nm1 = domain->readVectorFromFile(name, 1, 0, *Unm1);
      if (isGFMPAR || ( iod.problem.framework == ProblemData::EMBEDDED ||  iod.problem.framework == ProblemData::EMBEDDEDALE )) {
        //domain->getCommunicator()->fprintf(stderr,"*** Warning: Backward Euler being used instead of 3BDF for multiphase flows on first step after restart\n");
        data->exist_nm1 = false;
      }
    }
    if (data->use_nm2)
      data->exist_nm2 = domain->readVectorFromFile(name, 2, 0, *Unm2);
  }

  U = *Un;
  *Vn = *Un;
  if (data->use_nm1 && !data->exist_nm1)
    *Unm1 = *Un;
  if (data->use_nm2 && !data->exist_nm2)
    *Unm2 = *Unm1;
  
  createSubStates();
}

template<int dim>
void DistTimeState<dim>::createSubStates() {

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub)
    if (!subTimeState[iSub])
      subTimeState[iSub] = new TimeState<dim>(*data, (*dt)(iSub), (*idti)(iSub), (*idtv)(iSub), (*dtau)(iSub),
                                              (*Un)(iSub), (*Unm1)(iSub), (*Unm2)(iSub), (*Rn)(iSub));

}

//------------------------------------------------------------------------------
template<int dim>
void DistTimeState<dim>::copyTimeData(DistTimeState<dim>* oth) {

  data->copy(*oth->data);
}

template<int dim>
void DistTimeState<dim>::computeInitialState(InitialConditions &ic,
                                             FluidModelData &fm, double UU[dim])
{

  if(fm.fluid == FluidModelData::PERFECT_GAS || fm.fluid == FluidModelData::STIFFENED_GAS){
    double gam = fm.gasModel.specificHeatRatio;
    double ps = fm.gasModel.pressureConstant;

    double rho = ic.density;
    double p   = ic.pressure;
    double vel = 0.0;
    if(ic.mach>=0.0) vel = ic.mach*sqrt(gam*(p+ps)/rho);
    else vel = ic.velocity;
    double u   = vel*cos(ic.alpha)*cos(ic.beta);
    double v   = vel*cos(ic.alpha)*sin(ic.beta);
    double w   = vel*sin(ic.alpha);

    UU[0] = rho;
    UU[1] = rho*u;
    UU[2] = rho*v;
    UU[3] = rho*w;
    UU[4] = (p+gam*ps)/(gam-1.0) + 0.5 *rho*vel*vel;

  }else if(fm.fluid == FluidModelData::JWL){
    double omega = fm.jwlModel.omega;
    double A1    = fm.jwlModel.A1;
    double A2    = fm.jwlModel.A2;
    double R1    = fm.jwlModel.R1;
    double R2    = fm.jwlModel.R2;
    double rhor  = fm.jwlModel.rhoref;
    double R1r = R1*rhor; double R2r = R2*rhor;

    double rho = ic.density;
    double p   = ic.pressure;

    double frho  = A1*(1-omega*rho/R1r)*exp(-R1r/rho) + A2*(1-omega*rho/R2r)*exp(-R2r/rho);
    double frhop = A1*(-omega/R1r + (1-omega*rho/R1r)*R1r/(rho*rho)) *exp(-R1r/rho)
                 + A2*(-omega/R2r + (1-omega*rho/R2r)*R2r/(rho*rho)) *exp(-R2r/rho);

    double vel = 0.0;
    if(ic.mach>=0.0) vel = ic.mach*sqrt(((omega+1.0)*p-frho)/rho + frhop);
    else vel = ic.velocity;
    double u   = vel*cos(ic.alpha)*cos(ic.beta);
    double v   = vel*cos(ic.alpha)*sin(ic.beta);
    double w   = vel*sin(ic.alpha);

    UU[0] = rho;
    UU[1] = rho*u;
    UU[2] = rho*v;
    UU[3] = rho*w;
    UU[4] = (p-frho)/omega + 0.5 *rho*vel*vel;

  }else if(fm.fluid == FluidModelData::LIQUID){
    double pref  = fm.liquidModel.Pref;
    double alpha = fm.liquidModel.alpha;
    double beta  = fm.liquidModel.beta;
    double cv    = fm.liquidModel.specificHeat;

    double rho = ic.density;
    double temperature = ic.temperature;
    double vel = 0.0;
    if(ic.mach>=0.0) vel = ic.mach*sqrt(alpha*beta*pow(rho,beta-1.0));
    else vel = ic.velocity;
    double u   = vel*cos(ic.alpha)*cos(ic.beta);
    double v   = vel*cos(ic.alpha)*sin(ic.beta);
    double w   = vel*sin(ic.alpha);

    UU[0] = rho;
    UU[1] = rho*u;
    UU[2] = rho*v;
    UU[3] = rho*w;
    UU[4] = rho*(cv*temperature + 0.5*vel*vel) - (pref + alpha*pow(rho, beta));

  }else{
    fprintf(stderr, "*** Error: no initial state could be computed\n");
    exit(1);
  }

}

//------------------------------------------------------------------------------
 
template<int dim>
void DistTimeState<dim>::setupUVolumesInitialConditions(IoData &iod)
{

  // loop on all Volumes to setup U0
  if(!iod.volumes.volumeMap.dataMap.empty()){
    map<int, VolumeData *>::iterator volIt;
    for (volIt=iod.volumes.volumeMap.dataMap.begin(); volIt!=iod.volumes.volumeMap.dataMap.end();volIt++)
      if(volIt->second->type==VolumeData::FLUID || volIt->second->type==VolumeData::POROUS){
        //each volume (volIt->first) is setup using Input variables 'volumeInitialConditions'
        //                                 and equation of state 'fluidModel'
        map<int, FluidModelData *>::iterator fluidIt = iod.eqs.fluidModelMap.dataMap.find(volIt->second->fluidModelID);
        if(fluidIt == iod.eqs.fluidModelMap.dataMap.end()){
          fprintf(stderr, "*** Error: fluidModelData[%d] could not be found\n", volIt->second->fluidModelID);
          exit(1);
        }
        double UU[dim];
        computeInitialState(volIt->second->initialConditions, *fluidIt->second, UU);
/*        domain->getCommunicator()->fprintf(stdout, "- Initializing volume %d(EOS=%d) with \n", volIt->first, volIt->second->fluidModelID);
			 domain->getCommunicator()->fprintf(stdout, "    non-dimensionalized conservative state vector: (%g %g %g %g %g).\n", UU[0],UU[1],UU[2],UU[3],UU[4]);*/

        domain->setupUVolumesInitialConditions(volIt->first, UU, *Un);
      }
  }

}
//------------------------------------------------------------------------------
template<int dim>
void DistTimeState<dim>::setupUOneDimensionalSolution(IoData &iod, DistSVec<double,3> &X)
{
  OneDimensional::read1DSolution(iod, *Un, 
				 (DistSVec<double,1>*)NULL,NULL,
				 varFcn, X,*domain,
				 OneDimensional::ModeU);
}
//------------------------------------------------------------------------------
template<int dim>
void DistTimeState<dim>::setupUExactSolutionInitialConditions(IoData &iod, DistSVec<double,3> &X)
{
  // Test case one:
  if (iod.embed.testCase == 1) {

#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      SVec<double,dim> &u((*Un)(iSub));
      SVec<double, 3> &x(X(iSub));
      for(int i=0; i<X.subSize(iSub); i++) {

	double v[dim];
	
	ExactSolution::AcousticBeam(iod,x[i][0],x[i][1],x[i][2], 0.0, v);

	varFcn->primitiveToConservative(v, u[i], 0);
      }
    }
  } else if (iod.embed.testCase == 2) {

#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      SVec<double,dim> &u((*Un)(iSub));
      SVec<double, 3> &x(X(iSub));
      for(int i=0; i<X.subSize(iSub); i++) {

	double v[dim];
	
	ExactSolution::AcousticViscousBeam(iod,x[i][0],x[i][1],x[i][2], 0.0, v);

	varFcn->primitiveToConservative(v, u[i], 0);
      }
    }
  }

  // Test case two:
  if (iod.mf.testCase == 2) {

    DistSVec<double,1> dummy1(X.info());
    DistVec<int> dummy2(X.info());
    
    ExactSolution::Fill<&ExactSolution::CylindricalBubble,
      dim, 1>(*Un,dummy2,
	      dummy1, X, iod,0.0,
	      varFcn);
  } else if (iod.mf.testCase == 3) {

    DistSVec<double,1> dummy1(X.info());
    DistVec<int> dummy2(X.info());
    
    ExactSolution::Fill<&ExactSolution::AcousticTwoFluid,
      dim, 1>(*Un,dummy2,
	      dummy1, X, iod,0.0,
	      varFcn);
  }
}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::setupUMultiFluidInitialConditions(IoData &iod, DistSVec<double,3> &X)
{
  // Note that Un was already initialized either using the far-field
  // or using the definition of volumes. Therefore the only initialization
  // left to do is for the planes and spheres (in that order).

  map<int, FluidModelData *>::iterator fluidIt;

  if(!iod.mf.multiInitialConditions.planeMap.dataMap.empty()){
    map<int, PlaneData *>::iterator planeIt;
    for(planeIt  = iod.mf.multiInitialConditions.planeMap.dataMap.begin();
        planeIt != iod.mf.multiInitialConditions.planeMap.dataMap.end();
        planeIt++){
      fluidIt = iod.eqs.fluidModelMap.dataMap.find(planeIt->second->fluidModelID);
      if(fluidIt == iod.eqs.fluidModelMap.dataMap.end()){
        fprintf(stderr, "*** Error: fluidModelData[%d] could not be found\n", planeIt->second->fluidModelID);
        exit(1);
      }
      double UU[dim];
      computeInitialState(planeIt->second->initialConditions, *fluidIt->second, UU);
/*      domain->getCommunicator()->fprintf(stdout, "- Initializing PlaneData[%d] = (%g %g %g), (%g %g %g) with \n", planeIt->first, planeIt->second->cen_x, planeIt->second->cen_y,planeIt->second->cen_z,planeIt->second->nx,planeIt->second->ny,planeIt->second->nz);
		  domain->getCommunicator()->fprintf(stdout, "    EOS %d and non-dimensionalized conservative state vector: (%g %g %g %g %g).\n",  planeIt->second->fluidModelID, UU[0],UU[1],UU[2],UU[3],UU[4]);*/

#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        SVec<double,dim> &u((*Un)(iSub));
        SVec<double, 3> &x(X(iSub));

        double scalar = 0.0;
        for(int i=0; i<u.size(); i++) {
          scalar = planeIt->second->nx*(x[i][0] - planeIt->second->cen_x)+planeIt->second->ny*(x[i][1] - planeIt->second->cen_y)+planeIt->second->nz*(x[i][2] - planeIt->second->cen_z);
          if(scalar > 0.0) //node is on the same side indicated by vector
            for (int idim=0; idim<dim; idim++)
              u[i][idim] = UU[idim];
        }
      }
    }
  }

  if(!iod.mf.multiInitialConditions.cylinderMap.dataMap.empty()){
    map<int, CylinderData *>::iterator cylinderIt;
    for(cylinderIt  = iod.mf.multiInitialConditions.cylinderMap.dataMap.begin();
        cylinderIt != iod.mf.multiInitialConditions.cylinderMap.dataMap.end();
        cylinderIt++){
      fluidIt = iod.eqs.fluidModelMap.dataMap.find(cylinderIt->second->fluidModelID);
      if(fluidIt == iod.eqs.fluidModelMap.dataMap.end()){
        fprintf(stderr, "*** Error: fluidModelData[%d] could not be found\n", cylinderIt->second->fluidModelID);
        exit(1);
      }
      double UU[dim];
      computeInitialState(cylinderIt->second->initialConditions, *fluidIt->second, UU);
/*      domain->getCommunicator()->fprintf(stdout, "- Initializing PlaneData[%d] = (%g %g %g), (%g %g %g) with \n", planeIt->first, planeIt->second->cen_x, planeIt->second->cen_y,planeIt->second->cen_z,planeIt->second->nx,planeIt->second->ny,planeIt->second->nz);
		  domain->getCommunicator()->fprintf(stdout, "    EOS %d and non-dimensionalized conservative state vector: (%g %g %g %g %g).\n",  planeIt->second->fluidModelID, UU[0],UU[1],UU[2],UU[3],UU[4]);*/

#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        SVec<double,dim> &u((*Un)(iSub));
        SVec<double, 3> &x(X(iSub));

        double scalar = 0.0;
        for(int i=0; i<u.size(); i++) {
          scalar = cylinderIt->second->nx*(x[i][0] - cylinderIt->second->cen_x)+cylinderIt->second->ny*(x[i][1] - cylinderIt->second->cen_y)+cylinderIt->second->nz*(x[i][2] - cylinderIt->second->cen_z);
	  if (scalar < 0.0 || scalar > cylinderIt->second->L)
	    continue;
	  
	  double q[3] = {(x[i][0] - cylinderIt->second->cen_x) - scalar*cylinderIt->second->nx,
			 (x[i][1] - cylinderIt->second->cen_y) - scalar*cylinderIt->second->ny,
			 (x[i][2] - cylinderIt->second->cen_z) - scalar*cylinderIt->second->nz};
	  
          if(sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]) < 
	     cylinderIt->second->r) //node is on the same side indicated by vector
            for (int idim=0; idim<dim; idim++)
              u[i][idim] = UU[idim];
        }
      }
    }
  }

  if(!iod.mf.multiInitialConditions.sphereMap.dataMap.empty()){
    map<int, SphereData *>::iterator sphereIt;
    for(sphereIt  = iod.mf.multiInitialConditions.sphereMap.dataMap.begin();
        sphereIt != iod.mf.multiInitialConditions.sphereMap.dataMap.end();
        sphereIt++){
      fluidIt = iod.eqs.fluidModelMap.dataMap.find(sphereIt->second->fluidModelID);
      if(fluidIt == iod.eqs.fluidModelMap.dataMap.end()){
        fprintf(stderr, "*** Error: fluidModelData[%d] could not be found\n", sphereIt->second->fluidModelID);
        exit(1);
      }
      double UU[dim];
      computeInitialState(sphereIt->second->initialConditions, *fluidIt->second, UU);
/*      domain->getCommunicator()->fprintf(stdout, "- Initializing SphereData[%d] = (%g %g %g), %g with EOS %d\n", sphereIt->first, sphereIt->second->cen_x, sphereIt->second->cen_y,sphereIt->second->cen_z,sphereIt->second->radius, sphereIt->second->fluidModelID);      domain->getCommunicator()->fprintf(stdout, "    and non-dimensionalized conservative state vector: (%g %g %g %g %g).\n", UU[0],UU[1],UU[2],UU[3],UU[4]);*/


#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        SVec<double,dim> &u((*Un)(iSub));
        SVec<double, 3> &x(X(iSub));

        double dist = 0.0;
        for(int i=0; i<X.subSize(iSub); i++) {
          dist = (x[i][0] - sphereIt->second->cen_x)*(x[i][0] - sphereIt->second->cen_x) +
                 (x[i][1] - sphereIt->second->cen_y)*(x[i][1] - sphereIt->second->cen_y) +
                 (x[i][2] - sphereIt->second->cen_z)*(x[i][2] - sphereIt->second->cen_z);
          if(dist < sphereIt->second->radius*sphereIt->second->radius){ //node is inside the sphere
            for (int idim=0; idim<dim; idim++)
              u[i][idim] = UU[idim];
          }
        }
      }
    }
  }

  if(!iod.mf.multiInitialConditions.prismMap.dataMap.empty()){
    map<int, PrismData *>::iterator prismIt;
    for(prismIt  = iod.mf.multiInitialConditions.prismMap.dataMap.begin();
        prismIt != iod.mf.multiInitialConditions.prismMap.dataMap.end();
        prismIt++){
      fluidIt = iod.eqs.fluidModelMap.dataMap.find(prismIt->second->fluidModelID);
      if(fluidIt == iod.eqs.fluidModelMap.dataMap.end()){
        fprintf(stderr, "*** Error: fluidModelData[%d] could not be found\n", prismIt->second->fluidModelID);
        exit(1);
      }
      double UU[dim];
      computeInitialState(prismIt->second->initialConditions, *fluidIt->second, UU);
/*      domain->getCommunicator()->fprintf(stdout, "- Initializing PrismData[%d] = (%g %g %g with EOS %d\n", prismIt->first, prismIt->second->cen_x, prismIt->second->cen_y,prismIt->second->cen_z, prismIt->second->fluidModelID);
		  domain->getCommunicator()->fprintf(stderr, "    and non-dimensionalized primitive state vector: (%g %g %g %g %g).\n", UU[0],UU[1],UU[2],UU[3],UU[4]);*/
#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        SVec<double,dim> &u((*Un)(iSub));
        SVec<double, 3> &x(X(iSub));

        double dist = 0.0;
        for(int i=0; i<X.subSize(iSub); i++) {
          if(prismIt->second->inside(x[i][0],x[i][1],x[i][2])) { //node is inside the sphere
	    for (int idim=0; idim<dim; idim++)
	      u[i][idim] = UU[idim];
	  }
        }
      }
    }
  }

  // Test case one:
  // p = 1+10^8(x-0.2)^4*(x-0.6)^4, [x = (0.2,0.6)]
  if (iod.mf.testCase == 1) {

    int oth = (iod.eqs.fluidModelMap.dataMap.size() > 1 ? 1 : 0);   
 
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {
      SVec<double,dim> &u((*Un)(iSub));
      SVec<double, 3> &x(X(iSub));
      for(int i=0; i<X.subSize(iSub); i++) {
	double press = 1.0;
	if (x[i][0] > 0.499-0.2 && x[i][0] < 0.499+0.2)
	  press += 1.0e6*pow(x[i][0]-(0.499-0.2),4.0)*pow(x[i][0]-(0.499+0.2),4.0);
	double v[dim] = {0,0,0,0,0};
	v[0] = u[i][0];
	v[4] = press / iod.ref.rv.pressure;
	varFcn->primitiveToConservative(v, u[i], (x[i][0] < 0.499 ? 0 : oth) );
      }
    }
  }
}
//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::setupUFluidIdInitialConditions(IoData &iod, DistVec<int> &pointId)
{
  map<int, FluidModelData *>::iterator fluidIt;

  if(!iod.embed.embedIC.pointMap.dataMap.empty()){
    map<int, PointData *>::iterator pointIt;
    int count = 0;
    for(pointIt  = iod.embed.embedIC.pointMap.dataMap.begin();
        pointIt != iod.embed.embedIC.pointMap.dataMap.end();
        pointIt ++){

      count++;
      fluidIt = iod.eqs.fluidModelMap.dataMap.find(pointIt->second->fluidModelID);
      if(fluidIt == iod.eqs.fluidModelMap.dataMap.end()){
        fprintf(stderr, "*** Error: fluidModelData[%d] could not be found\n", pointIt->second->fluidModelID);
        exit(-1);
      }

      double UU[dim];
      computeInitialState(pointIt->second->initialConditions, *fluidIt->second, UU);
/*      domain->getCommunicator()->fprintf(stdout, "- Initializing PointData[%d] = (%g %g %g), with EOS %d\n", pointIt->first, pointIt->second->x, pointIt->second->y,pointIt->second->z, pointIt->second->fluidModelID);
		  domain->getCommunicator()->fprintf(stdout, "    and non-dimensionalized conservative state vector: (%g %g %g %g %g).\n", UU[0],UU[1],UU[2],UU[3],UU[4]);*/

#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        SVec<double,dim> &subUn((*Un)(iSub));
        Vec<int> &subId(pointId(iSub));

        for(int i=0; i<subUn.size(); i++)
          if(subId[i]==count)
            for(int iDim=0; iDim<dim; iDim++)
              subUn[i][iDim] = UU[iDim];
      }
    }
  } else {
    domain->getCommunicator()->fprintf(stderr, "NOTE: FluidId-based initial conditions not specified.\n");
  }
} 

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void DistTimeState<dim>::rstVar(IoData & ioData) 
{

  //mach = ioData.bc.inlet.mach;

}

//------------------------------------------------------------------------------

template<int dim>
double DistTimeState<dim>::computeTimeStep(double cfl, double dualtimecfl, double* dtLeft, int* numSubCycles,
					   DistGeoState &geoState, DistSVec<double,3> &X,
														 DistVec<double> &ctrlVol, DistSVec<double,dim> &U, 
														 DistLevelSetStructure *distLSS)
{

	if(distLSS)
	{
		varFcn->conservativeToPrimitive(U, *V, distLSS);

		domain->computeTimeStep(cfl, dualtimecfl, viscousCst, fet, varFcn, geoState, X, ctrlVol, 
										*V, *dt, *idti, *idtv, *dtau, *irey, tprec, sprec, distLSS);
	}
	else
{
  varFcn->conservativeToPrimitive(U, *V);

		domain->computeTimeStep(cfl, dualtimecfl, viscousCst, fet, varFcn, geoState, X, ctrlVol, 
										*V, *dt, *idti, *idtv, *dtau, *irey, tprec, sprec);
	}

  double dt_glob;
  updateDtCoeff();

	if(data->dt_imposed > 0.0) 
	{
    dt_glob = data->dt_imposed;
    allowcflstop = false; 
    dt_glob *= dt_coeff;
  }
	else 
	{
    dt_glob = dt->min();
    allowdtstop = false;
  }

  if (data->typeStartup == ImplicitData::MODIFIED && 
      ((data->typeIntegrator == ImplicitData::THREE_POINT_BDF && !data->exist_nm1) ||
	    (data->typeIntegrator == ImplicitData::FOUR_POINT_BDF  && (!data->exist_nm2 || !data->exist_nm1)))) 
	{
		if(*dtLeft != 0.0 && dt_glob > *dtLeft)  dt_glob = *dtLeft / 1000.0;
		else dt_glob = min(dt_glob, dt_glob*dt_glob); //dt_glob /= 1000.0;
  }

	if(*dtLeft != 0.0) 
	{
    *numSubCycles = int(*dtLeft / dt_glob);
    if (*numSubCycles == 0 || (*numSubCycles)*dt_glob != *dtLeft) ++(*numSubCycles);
    dt_glob = *dtLeft / double(*numSubCycles);
    *dtLeft -= dt_glob;
  }
  data->computeCoefficients(*dt, dt_glob);

  return dt_glob;

}

//------------------------------------------------------------------------------

template<int dim>
double DistTimeState<dim>::computeTimeStepFailSafe(double* dtLeft, int* numSubCycles)
{
  // time step is repeated with half the time step size
  double dt_glob;
  *dtLeft += dt->min();
  dt_glob = dt->min() / 2.0;

  *dtLeft -= dt_glob;

  data->computeCoefficients(*dt, dt_glob);

  return dt_glob;

}

//------------------------------------------------------------------------------

template<int dim> 
double DistTimeState<dim>::computeTimeStep(int it, double* dtLeft, int* numSubCycles)
{
  if (data->dt_imposed > 0.0 && *dtLeft == 0.0) 
  //Allows for use of subcycling in fluid only simulations with imposed time step
    *dtLeft = data->dt_imposed;

  double incfac = 1.25 + (1.15 * pow((2.71828),(- double(it-2) / 3.0)));
  double decfac = max(0.2 , (0.75 + (-1.25 * pow((2.71828),(- double(it-2) / 3.0)))));

  double factor;
  if (errorEstiNorm == 0.0)
    printf("Error: errorEstiNorm equals zero! Division by zero! \n");
  else   
    factor = min( max( pow((data->errorTol / errorEstiNorm),(1.0 / 2.0)) , decfac ), incfac );

  if(factor==decfac && it==2)
    printf("WARNING: Cfl0 is chosen too big!! \n");

  double dt_glob;
  updateDtCoeff();
  if (data->dt_imposed > 0.0) {
    dt_glob = data->dt_imposed;
    allowcflstop = false; 
    dt_glob *= dt_coeff;
  }
  else{ 
    allowdtstop = false;
    dt_glob = max(dtMin, (factor * data->dt_nm1));
  }

  if (data->typeStartup == ImplicitData::MODIFIED && 
      ((data->typeIntegrator == ImplicitData::THREE_POINT_BDF && !data->exist_nm1) ||
       (data->typeIntegrator == ImplicitData::FOUR_POINT_BDF && (!data->exist_nm2 || !data->exist_nm1)))) {
    if (*dtLeft != 0.0 && dt_glob > *dtLeft)
      dt_glob = *dtLeft / 1000.0;
    else 
      dt_glob /= 1000.0;
  }

  if (*dtLeft != 0.0) {
    *numSubCycles = int(*dtLeft / dt_glob);
    if (*numSubCycles == 0 || (*numSubCycles)*dt_glob != *dtLeft) ++(*numSubCycles);
    dt_glob = *dtLeft / double(*numSubCycles);
    *dtLeft -= dt_glob;
  }

  data->computeCoefficients(*dt, dt_glob);

  return dt_glob;

}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::updateDtCoeff()
{

  if(errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP_TIME]) {
    errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP_TIME] = 0;
    dt_coeff_count=0;
    dt_coeff /= 2.0;
    if(dt_coeff<0.001 && allowdtstop) {
      domain->getCommunicator()->fprintf(stdout, "Could not resolve unphysicality by reducing timestep. Aborting.");
      exit(-1);
    }
  }
  dt_coeff_count++;
  
  if(dt_coeff_count>4) {
    dt_coeff *= 2.0;
    dt_coeff=min(dt_coeff,1.0);
    dt_coeff_count = 2;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::calculateErrorEstiNorm(DistSVec<double,dim> &U, DistSVec<double,dim> &F)
{
  //linear extrapolation for error estimation
  //F = *Un + ( data->dt_n / data->dt_nm1 )*( *Un-*Unm1 );

  //Forward Euler step for error estimation
  F = *Un - ( data->dt_n * F);

  //common for both approaches
  F -= U;
  errorEstiNorm = F.norm() / U.norm();
}

//------------------------------------------------------------------------------
                 
//dd                                               
template<int dim>
double DistTimeState<dim>::computeTimeStep(double cfl, double dualtimecfl, double* dtLeft, int* numSubCycles,
                                           DistGeoState &geoState, DistVec<double> &ctrlVol,
                                           DistSVec<double,dim> &U, DistVec<int> &fluidId,
					   DistVec<double>* umax)
{
  varFcn->conservativeToPrimitive(U, *V, &fluidId);

	domain->computeTimeStep(cfl, dualtimecfl, viscousCst, fet, varFcn, geoState, ctrlVol, 
									*V, *dt, *idti, *idtv, *dtau, tprec, fluidId, umax);
                                                                                                         
  double dt_glob;
  updateDtCoeff();
  
  if(data->dt_imposed > 0.0)
  {
    dt_glob = data->dt_imposed;
    allowcflstop = false; 
    dt_glob *= dt_coeff;
  }
  else
  {
    dt_glob = dt->min();
    allowdtstop = false;
  }
                                                                                                         
  if(umax && isGFMPAR) 
  {
    double udt = umax->min();

	  if(udt < dt_glob) 
	  {
      domain->getCommunicator()->fprintf(stdout, "Clamped new dt %lf (old = %lf)\n", udt*refTime, dt_glob*refTime);
      domain->getCommunicator()->fprintf(stdout, "*** Warning: Cfl for this multi-phase algorithm has been clamped to \n");
      domain->getCommunicator()->fprintf(stdout, "             %lf (user specified %lf)\n", udt/dt_glob*cfl, cfl);
      dt_glob = udt;
    }
  }

  if (data->typeStartup == ImplicitData::MODIFIED &&
      ((data->typeIntegrator == ImplicitData::THREE_POINT_BDF && !data->exist_nm1) ||
      (data->typeIntegrator == ImplicitData::FOUR_POINT_BDF && (!data->exist_nm2 || !data->exist_nm1)))) 
  {
	  if (*dtLeft != 0.0 && dt_glob > *dtLeft) dt_glob = *dtLeft / 1000.0;
	  else dt_glob = min(dt_glob,dt_glob*dt_glob);//dt_glob /= 1000.0;
  }
  
  if(*dtLeft != 0.0) 
  {
    *numSubCycles = int(*dtLeft / dt_glob);
    if (*numSubCycles == 0 || (*numSubCycles)*dt_glob != *dtLeft) ++(*numSubCycles);
    dt_glob = *dtLeft / double(*numSubCycles);
    *dtLeft -= dt_glob;
  }
                                                                                            
  data->computeCoefficients(*dt, dt_glob);
  
  return dt_glob;
                                                                                                         
}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::computeCoefficients(double dt_glob)
{

  data->computeCoefficients(*dt, dt_glob);

}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::add_dAW_dt(int it, DistGeoState &geoState, 
				    DistVec<double> &ctrlVol,
				    DistSVec<double,dim> &Q, 
				    DistSVec<double,dim> &R, DistLevelSetStructure *distLSS)
{
  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  if (data->typeIntegrator == ImplicitData::CRANK_NICOLSON && it == 0) *Rn = R;

#pragma omp parallel for
	for(int iSub = 0; iSub < numLocSub; ++iSub) 
	{
    LevelSetStructure *LSS = distLSS ? &((*distLSS)(iSub)) : 0;

		if(tprec.timePreconditioner() && data->typeTimeStep == TsData::LOCAL) 
		{
      if(varFcn->getType() == VarFcnBase::PERFECTGAS || varFcn->getType() == VarFcnBase::STIFFENEDGAS)
        subTimeState[iSub]->add_GASPrec_dAW_dt(Q.getMasterFlag(iSub), geoState(iSub), ctrlVol(iSub), Q(iSub), R(iSub), gam, pstiff, (*irey)(iSub), tprec, LSS);

      else if(varFcn->getType() == VarFcnBase::TAIT)
        subTimeState[iSub]->add_LiquidPrec_dAW_dt(Q.getMasterFlag(iSub), geoState(iSub), ctrlVol(iSub), varFcn, Q(iSub), R(iSub), (*irey)(iSub), tprec, LSS);

			else 
			{
        fprintf(stderr, "*** Error: no time preconditioner for this EOS  *** EXITING\n");
        exit(1);
      }
    }
		else 
		{
      subTimeState[iSub]->add_dAW_dt(Q.getMasterFlag(iSub), geoState(iSub), ctrlVol(iSub), Q(iSub), R(iSub), LSS);
    }
  }
}   

template<int dim>
void DistTimeState<dim>::add_dAW_dt_HH(int it, DistGeoState &geoState, 
		   DistVec<double> &ctrlVol,
		   DistVec<double> &Q, 
		   DistVec<double> &R) {

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    subTimeState[iSub]->add_dAW_dt_HH(Q.getMasterFlag(iSub), geoState(iSub), 
				      ctrlVol(iSub), Q(iSub), R(iSub));
  }
  
}                                                                                                                 
//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::add_dAW_dtRestrict(int it, DistGeoState &geoState, 
					    DistVec<double> &ctrlVol,
					    DistSVec<double,dim> &Q, 
					    DistSVec<double,dim> &R, const std::vector<std::vector<int> > &sampledLocNodes)
{

  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  if (data->typeIntegrator == ImplicitData::CRANK_NICOLSON && it == 0) *Rn = R;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->add_dAW_dtRestrict(Q.getMasterFlag(iSub), geoState(iSub), 
				   ctrlVol(iSub), Q(iSub), R(iSub), sampledLocNodes[iSub]);

}
//------------------------------------------------------------------------------
template<int dim>
template<int dimLS>
void DistTimeState<dim>::add_dAW_dtLS(int it, DistGeoState &geoState,
                                      DistVec<double> &ctrlVol,
                                      DistSVec<double,dimLS> &Q,
                                      DistSVec<double,dimLS> &Qn,
                                      DistSVec<double,dimLS> &Qnm1,
                                      DistSVec<double,dimLS> &Qnm2,
                                      DistSVec<double,dimLS> &R,bool requireSpecialBDF)
{

  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  //if (data->typeIntegrator == ImplicitData::CRANK_NICOLSON && it == 0) *Rn = R;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->add_dAW_dtLS(Q.getMasterFlag(iSub), geoState(iSub),
                                     ctrlVol(iSub), Q(iSub), Qn(iSub), Qnm1(iSub),
				                             Qnm2(iSub), R(iSub),requireSpecialBDF);
}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::add_dAW_dtau(int it, DistGeoState &geoState, 
				    DistVec<double> &ctrlVol,
				    DistSVec<double,dim> &Q, 
				    DistSVec<double,dim> &R, DistLevelSetStructure *distLSS)
{

  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

//  if (data->typeIntegrator == ImplicitData::CRANK_NICOLSON && it == 0) *Rn = R;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    LevelSetStructure *LSS = distLSS ? &((*distLSS)(iSub)) : 0;

    if(tprec.timePreconditioner()) {
      if(varFcn->getType() == VarFcnBase::PERFECTGAS || varFcn->getType() == VarFcnBase::STIFFENEDGAS)
        subTimeState[iSub]->add_GASPrec_dAW_dtau(Q.getMasterFlag(iSub), geoState(iSub), ctrlVol(iSub), Q(iSub), R(iSub), gam, pstiff, (*irey)(iSub), tprec, LSS);

      else if(varFcn->getType() == VarFcnBase::TAIT)
        subTimeState[iSub]->add_LiquidPrec_dAW_dtau(Q.getMasterFlag(iSub), geoState(iSub), ctrlVol(iSub), varFcn, Q(iSub), R(iSub), (*irey)(iSub), tprec, LSS);

      else {
        fprintf(stderr, "*** Error: no time preconditioner for this EOS  *** EXITING\n");
        exit(1);
      }
    }
    else {
      subTimeState[iSub]->add_dAW_dtau(Q.getMasterFlag(iSub), geoState(iSub), ctrlVol(iSub), Q(iSub), R(iSub), LSS);
    }
  }
}                                                                                                                      
//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void DistTimeState<dim>::addToJacobian(DistVec<double> &ctrlVol, DistMat<Scalar,neq> &A,
                                       DistSVec<double,dim> &U)
{
  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  if(tprec.timePreconditioner()){
    if(varFcn->getType() == VarFcnBase::PERFECTGAS || varFcn->getType() == VarFcnBase::STIFFENEDGAS)
      addToJacobianGasPrec(ctrlVol, A, U);
    else if(varFcn->getType() == VarFcnBase::TAIT)
      addToJacobianLiquidPrec(ctrlVol, A, U);
    else{
      fprintf(stderr, "*** Error: no time preconditioner for this EOS  *** EXITING\n");
      exit(1);
    }
  }else{
     addToJacobianNoPrec(ctrlVol, A, U);
  }

}

template<int dim>
template<class Scalar, int neq>
void DistTimeState<dim>::addToHHJacobian(DistVec<double> &ctrlVol, DistMat<Scalar,neq> &A,
                                       DistVec<double> &hh)
{
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->addToJacobianHH(ctrlVol(iSub), A(iSub), hh(iSub));
}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void DistTimeState<dim>::addToJacobianLS(DistVec<double> &ctrlVol, DistMat<Scalar,neq> &A,
                                       DistSVec<double,dim> &U,bool requireSpecialBDF)
{
  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->addToJacobianLS(V->getMasterFlag(iSub), ctrlVol(iSub), A(iSub), U(iSub),
                              requireSpecialBDF);
}

//------------------------------------------------------------------------------
template<int dim>
template<class Scalar, int neq>
void DistTimeState<dim>::addToJacobianNoPrec(DistVec<double> &ctrlVol, DistMat<Scalar,neq> &A,
                                       DistSVec<double,dim> &U)
{
  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  int** nodeType = domain->getNodeTypeExtrapolation();
  if (nodeType){
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subTimeState[iSub]->addToJacobianNoPrec(V->getMasterFlag(iSub), ctrlVol(iSub), A(iSub), U(iSub),
                                varFcn, nodeType[iSub]);
  }else{
    int* empty = 0;
    #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subTimeState[iSub]->addToJacobianNoPrec(V->getMasterFlag(iSub), ctrlVol(iSub), A(iSub), U(iSub),
                                varFcn,empty);
  }
}

//------------------------------------------------------------------------------
template<int dim>
template<class Scalar, int neq>
void DistTimeState<dim>::addToJacobianGasPrec(DistVec<double> &ctrlVol, DistMat<Scalar,neq> &A,
                                       DistSVec<double,dim> &U)
{
  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  int** nodeType = domain->getNodeTypeExtrapolation();
  if (nodeType){
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subTimeState[iSub]->addToJacobianGasPrec(V->getMasterFlag(iSub), ctrlVol(iSub), A(iSub), U(iSub),
                                varFcn, gam, pstiff, tprec, (*irey)(iSub), nodeType[iSub]);
  }else{
    int* empty = 0;
    #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subTimeState[iSub]->addToJacobianGasPrec(V->getMasterFlag(iSub), ctrlVol(iSub), A(iSub), U(iSub),
                                varFcn, gam, pstiff, tprec, (*irey)(iSub), empty);
  }
}

//------------------------------------------------------------------------------
template<int dim>
template<class Scalar, int neq>
void DistTimeState<dim>::addToJacobianLiquidPrec(DistVec<double> &ctrlVol, DistMat<Scalar,neq> &A,
                                       DistSVec<double,dim> &U)
{
  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  int** nodeType = domain->getNodeTypeExtrapolation();
  if (nodeType){
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subTimeState[iSub]->addToJacobianLiquidPrec(V->getMasterFlag(iSub), ctrlVol(iSub), A(iSub), U(iSub),
                                varFcn, tprec, (*irey)(iSub), nodeType[iSub]);
  }else{
    int* empty = 0;
    #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subTimeState[iSub]->addToJacobianLiquidPrec(V->getMasterFlag(iSub), ctrlVol(iSub), A(iSub), U(iSub),
                                varFcn, tprec, (*irey)(iSub), empty);
  }
}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void DistTimeState<dim>::addToH1(DistVec<double> &ctrlVol, DistMat<Scalar,neq> &A)
{
  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->addToH1(V->getMasterFlag(iSub), ctrlVol(iSub), A(iSub));

}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void DistTimeState<dim>::addToH1(DistVec<double> &ctrlVol,
                DistMat<Scalar,neq> &A, Scalar shift)
{
  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->addToH1(V->getMasterFlag(iSub), ctrlVol(iSub),
                A(iSub), shift);

}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void DistTimeState<dim>::addToH2(DistVec<double> &ctrlVol, DistSVec<double,dim> &U,
				 DistMat<Scalar,neq> &A)
{
  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

#ifdef DOUBLE_CHECK
  varFcn->conservativeToPrimitive(U, *V);
#endif

  if(tprec.timePreconditioner()) {

    if(varFcn->getType() == VarFcnBase::PERFECTGAS || 
       varFcn->getType() == VarFcnBase::STIFFENEDGAS)

#pragma omp parallel for
      for (int iSub = 0; iSub < numLocSub; ++iSub)
        subTimeState[iSub]->addToH2GasPrec(V->getMasterFlag(iSub), varFcn, 
                    ctrlVol(iSub), (*V)(iSub), A(iSub), gam, pstiff, 
                    (*irey)(iSub), tprec);

    else if(varFcn->getType() == VarFcnBase::TAIT)

#pragma omp parallel for
      for (int iSub = 0; iSub < numLocSub; ++iSub)
        subTimeState[iSub]->addToH2LiquidPrec(V->getMasterFlag(iSub), varFcn, 
                    ctrlVol(iSub), (*V)(iSub), A(iSub), (*irey)(iSub), tprec);

    else {
      fprintf(stderr, "*** Error: no time preconditioner for this EOS  *** EXITING\n");
      exit(1);
    }
  }
  else {

#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subTimeState[iSub]->addToH2NoPrec(V->getMasterFlag(iSub), varFcn, 
                    ctrlVol(iSub), (*V)(iSub), A(iSub));

  }

//#pragma omp parallel for
//  for (int iSub = 0; iSub < numLocSub; ++iSub)
//    subTimeState[iSub]->addToH2(V->getMasterFlag(iSub), varFcn, ctrlVol(iSub), 
//				(*V)(iSub), A(iSub)); 

}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar, int neq>
void DistTimeState<dim>::addToH2(DistVec<double> &ctrlVol,
                DistSVec<double,dim> &U, DistMat<Scalar,neq> &A, Scalar shift)
{
  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

#ifdef DOUBLE_CHECK
  varFcn->conservativeToPrimitive(U, *V);
#endif

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->addToH2(V->getMasterFlag(iSub), varFcn, ctrlVol(iSub),
                                (*V)(iSub), A(iSub), shift);

}
//------------------------------------------------------------------------------


template<int dim>
template<class Scalar, int neq>
void DistTimeState<dim>::addToH2(DistVec<double> &ctrlVol,
                DistSVec<double,dim> &U, DistMat<Scalar,neq> &A, Scalar coefVol, double coefA)
{
  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

#ifdef DOUBLE_CHECK
  varFcn->conservativeToPrimitive(U, *V);
#endif

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->addToH2(V->getMasterFlag(iSub), varFcn, ctrlVol(iSub),
                                (*V)(iSub), A(iSub), coefVol, coefA);

}
//------------------------------------------------------------------------------




template<int dim>
template<class Scalar,int neq>
void DistTimeState<dim>::addToH2Minus(DistVec<double> &ctrlVol, DistSVec<double,dim> &U,
                                      DistMat<Scalar,neq> &A)
{
  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

#ifdef DOUBLE_CHECK
  varFcn->conservativeToPrimitive(U, *V);
#endif

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->addToH2Minus(V->getMasterFlag(iSub), varFcn, ctrlVol(iSub),
                                     (*V)(iSub), A(iSub));

}

//------------------------------------------------------------------------------
template<int dim>
void DistTimeState<dim>::multiplyByTimeStep(DistVec<double>& dU)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    double *du = dU.subData(iSub);
    double* _dt = dt->subData(iSub);
    for (int i=0; i<dU.subSize(iSub); ++i)
      du[i] *= _dt[i];
  }

}

template<int dim>
void DistTimeState<dim>::multiplyByTimeStep(DistSVec<double,dim>& dU)
{
  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    double (*du)[dim] = dU.subData(iSub);
    double* _dt = dt->subData(iSub);
    for (int i=0; i<dU.subSize(iSub); ++i)
      for (int j=0; j<dim; ++j)
	du[i][j] *= _dt[i];
  }

}

//------------------------------------------------------------------------------

template<int dim>
template<int dimLS>
void DistTimeState<dim>::multiplyByTimeStep(DistSVec<double,dimLS>& dPhi)
{
  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    double (*dphi)[dimLS] = dPhi.subData(iSub);
    double* _dt = dt->subData(iSub);
    for (int i=0; i<dPhi.subSize(iSub); ++i)
      for (int idim=0; idim<dimLS; idim++)
	      dphi[i][idim] *= _dt[i];
  }

}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::multiplyByPreconditioner(DistSVec<double,dim>& U0, DistSVec<double,dim>& dU)
{
  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

  if (tprec.timePreconditioner()){
    if (varFcn->getType() == VarFcnBase::PERFECTGAS || varFcn->getType() == VarFcnBase::STIFFENEDGAS)
      multiplyByPreconditionerPerfectGas(U0,dU);
    else if (varFcn->getType() == VarFcnBase::TAIT)
      multiplyByPreconditionerLiquid(U0,dU);
    else{
      fprintf(stderr, "*** Error: no time preconditioner for this EOS  *** EXITING\n");
      exit(1);
    }
  }
}

//------------------------------------------------------------------------------
template<int dim>
void DistTimeState<dim>::multiplyByPreconditionerPerfectGas(DistSVec<double,dim>& U0, DistSVec<double,dim>& dU)
{
  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    double* _irey = irey->subData(iSub);
    double (*u0)[dim] = U0.subData(iSub);
    double (*du)[dim] = dU.subData(iSub);
    for (int i=0; i<dU.subSize(iSub); ++i) {
        double ro = u0[i][0];
        double invRho = 1.0/ro;
        double u  = u0[i][1] * invRho;
        double v  = u0[i][2] * invRho;
        double w  = u0[i][3] * invRho;
        double u2 = u*u;
        double v2 = v*v;
        double w2 = w*w;
        double q2 = u2 + v2 + w2;
        double gam1 = gam - 1.0;
        double p  = gam1 * (u0[i][4] - 0.5 * ro * q2)-gam*pstiff;
        double c2 = gam*(p+pstiff)/ro;

        double locMach = sqrt(q2/c2); //local Preconditioning (ARL)
        double locbeta = tprec.getBeta(locMach,_irey[i]);
        //double locbeta = fmax(k1*locMach, beta);
        //locbeta = fmin((1.0+sqrt(_irey[i]))*locbeta,cmach);

        double beta2 =  locbeta * locbeta;
        double qhat2 = (q2 * gam1)/2.0;
                                                                                
        double nu = qhat2/c2;
        double mu = beta2 - 1.0;
        double H = (c2/gam1) + 0.5*q2;
        
 	// initialize output
	double temp[dim];
        for (int j=0; j<dim; j++)
          temp[j] = 0.0;

        // Euler or Navier-Stokes preconditioning
        double P[5][5] = {  {1.0 + nu*mu,  -u*mu*gam1/c2,      -v*mu*gam1/c2,      -w*mu*gam1/c2,       mu*gam1/c2   },
                            {u*nu*mu,     1.0 - u2*mu*gam1/c2, -u*v*mu*gam1/c2,      -u*w*mu*gam1/c2,     u*mu*gam1/c2 },
                            {v*nu*mu,     -u*v*mu*gam1/c2 ,    1.0 - v2*mu*gam1/c2,  -v*w*mu*gam1/c2,     v*mu*gam1/c2 },
                            {w*nu*mu,     -u*w*mu*gam1/c2 ,    -v*w*mu*gam1/c2,      1.0 - w2*mu*gam1/c2, w*mu*gam1/c2 },
                            {qhat2*H*mu/c2,  -u*mu*(1+nu),      -v*mu*(1+nu),       -w*mu*(1+nu), 1.0 + mu*gam1*H/c2 } };

       for (int j=0; j<5; ++j)
          for (int k=0; k<5; ++k)
             temp[j] += P[j][k]*du[i][k];

	//turbulence preconditioning
        if(dim==6){
          double t1 = u0[i][5] * invRho;
          double mup = mu*t1*gam1/c2;
          double Pt[6] = {mu*nu*t1, -mup*u, -mup*v, -mup*w, mup, 1.0};
          for (int k=0; k<6; k++)
            temp[5] += Pt[k]*du[i][k];

        }else if(dim==7){
          double t1 = u0[i][5] * invRho;
	  double t2 = u0[i][6] * invRho;
          double mup1 = mu*t1*gam1/c2;
          double mup2 = mu*t2*gam1/c2;
          double Pt[2][7] = { {mu*nu*t1, -mup1*u, -mup1*v, -mup1*w, mup1, 1.0, 0.0},
			      {mu*nu*t2, -mup2*u, -mup2*v, -mup2*w, mup2, 0.0, 1.0} };
          for (int k=0; k<7; k++){
            temp[5] += Pt[0][k]*du[i][k];
            temp[6] += Pt[1][k]*du[i][k];
          }
	}

	//get output
	for(int j=0; j<dim; ++j)
	   du[i][j] = temp[j];
    
    }
    
  }
}
                                                                                
//------------------------------------------------------------------------------
template<int dim>
void DistTimeState<dim>::multiplyByPreconditionerLiquid(DistSVec<double,dim> &U, DistSVec<double,dim> &dU)
{
  if (data->typeIntegrator == ImplicitData::SPATIAL_ONLY) return;

//ARL : turbulence preconditioner never tested...
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      double* _irey = irey->subData(iSub);
      double (*u0)[dim] = U.subData(iSub);
      double (*du)[dim] = dU.subData(iSub);
      for (int i=0; i<dU.subSize(iSub); ++i) {
        double V[dim];
        varFcn->conservativeToPrimitive(u0[i],V);
        double pressure = varFcn->getPressure(V);
        double e = varFcn->computeRhoEnergy(V)/V[0];
        double locMach = varFcn->computeMachNumber(V); //local Preconditioning (ARL)
        double locbeta = tprec.getBeta(locMach,_irey[i]);
        double beta2 = (locbeta*locbeta);
        double beta2m1 = beta2 - 1.0;
        
	// initialize output        
        double temp[dim];
        for (int j=0; j<dim; j++)
          temp[j] = 0.0;

	// preconditioning matrix
        double P[dim][dim];
 	// first initialize to identity
        for (int j=0; j<dim; j++)
          for (int k=0; k<dim; k++)
	     P[j][k] = 0.0;
        for (int j=0; j<dim; j++)
           P[j][j] = 1.0;
	// second get first column terms
        for (int j=0; j<dim; j++)
	  P[j][0] = beta2m1*V[j];
        P[0][0] = beta2;
        double c = varFcn->computeSoundSpeed(V);
        P[4][0] = beta2m1*(e+pressure/V[0]-c*c);
        
	/* The preconditioning matrix is:
         * P[dim][dim] = { { beta2,           0.0, 0.0, 0.0, 0.0 , 0.0, 0.0 },
         *                 {(beta2-1.0)*V[1], 1.0, 0.0, 0.0, 0.0 , 0.0, 0.0 },
         *                 {(beta2-1.0)*V[2], 0.0, 1.0, 0.0, 0.0 , 0.0, 0.0 },
         *                 {(beta2-1.0)*V[3], 0.0, 0.0, 1.0, 0.0 , 0.0, 0.0 },
         *                 {(beta2-1.0)*(h-c2),0.0, 0.0, 0.0, 1.0 , 0.0, 0.0 },
	 *		   {(beta2-1.0)*V[5], 0.0, 0.0, 0.0, 0.0 , 1.0, 0.0 },
	 *		   {(beta2-1.0)*V[6], 0.0, 0.0, 0.0, 0.0 , 0.0, 1.0 } };
	 * Take the first 5-by-5 matrix to get the Euler preconditioner
	 *      the first 6-by-6 matrix to get the "SA"  preconditioner
	 *      the whole 7-by-7 matrix to get the "k-e" preconditioner
         */
        for (int j=0; j<dim; ++j)
          for (int k=0; k<j+1; ++k)
            temp[j] += P[j][k]*du[i][k];

        for(int j=0; j<dim; ++j)
          du[i][j] = temp[j];
      }                                                                                                                                        
    }
}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::updateHH(DistVec<double> & hh) {

  *hhnm1 = *hhn;

  *hhn = hh;
}

struct SetFirstOrderNodes {

  VarFcn* varFcn;

  double threshold;
  ErrorHandler* errorHandler;

  int isDensity;

  SetFirstOrderNodes(VarFcn* varFcn,double t, ErrorHandler* errorHandlerIn,
		     int isDensity=0) : varFcn(varFcn), isDensity(isDensity),
      threshold(t) { errorHandler = errorHandlerIn; }

  void Perform(double* uold, double* unew, int& status,int id = 0) const {

    double vold[7], vnew[7];
    
    if (!isDensity) {
      varFcn->conservativeToPrimitive(uold,vold, id);
      varFcn->conservativeToPrimitive(unew,vnew, id);
      
      double pold = varFcn->getPressure(vold,id);
      double pnew = varFcn->getPressure(vnew,id);
      if (fabs(pnew-pold)/pold > threshold) {
	status = 1;
	errorHandler->localErrors[ErrorHandler::RAPIDLY_CHANGING_PRESSURE]++;
      }
    } else {
      varFcn->conservativeToPrimitive(uold,vold, id);
      varFcn->conservativeToPrimitive(unew,vnew, id);
      
      if (fabs(vnew[0]-vold[0])/vold[0] > threshold) {
	status = 1;
	errorHandler->localErrors[ErrorHandler::RAPIDLY_CHANGING_DENSITY]++;
      }
    }
  }
};

template<int dim>
void DistTimeState<dim>::update(DistSVec<double,dim> &Q,bool increasingPressure)
{

  data->update();

  *firstOrderNodes = 0;

  if (data->use_nm2 && data->exist_nm1) {
    *Unm2 = *Unm1;
    data->exist_nm2 = true;
  }
  // If we are increasing the fluid pressure there is
  // no relation between Un and Unm1; thus data->exist_nm1
  // should still be false.  This way when we start the fluid
  // solution startup will be done as necessary.
  if (data->use_nm1 && !increasingPressure) {
    *Unm1 = *Un;
    data->exist_nm1 = true;
  }

  if (checkForRapidlyChangingValues) {
    if (checkForRapidlyChangingPressure > 0.0)
      DistVectorOp::Op(*Un, Q,*firstOrderNodes,  
  		     SetFirstOrderNodes(varFcn,checkForRapidlyChangingPressure,errorHandler)); 
  
    if (checkForRapidlyChangingDensity > 0.0)
      DistVectorOp::Op(*Un, Q,*firstOrderNodes,
  		       SetFirstOrderNodes(varFcn,checkForRapidlyChangingDensity,errorHandler,1)); 
  }
  *Un = Q;
}

//------------------------------------------------------------------------------

struct MultiphaseRiemannCopy {
  int dim;
  MultiphaseRiemannCopy(int _dim) : dim(_dim) { }
  void Perform(double* vin, double* vout,int id1,int id2,int swept) const  {

    if (id1 != id2 && !swept) {
      for (int i = 0; i < dim; ++i) vout[i] = vin[i];
    }
  }
};


//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::update(DistSVec<double,dim> &Q, DistSVec<double,dim> &Qtilde,
				DistVec<int> &fluidId,
                                DistVec<int> *fluidIdnm1,
                                DistExactRiemannSolver<dim> *riemann,
                                DistLevelSetStructure* distLSS,bool increasingPressure)
{
  data->update();

  if (data->use_nm2 && data->exist_nm1) {
    fprintf(stdout, "4pt-BDF has not been studied for 2-phase flow\n");
    exit(1);
  }
  *firstOrderNodes = 0;
  if (data->use_nm1 && !increasingPressure) {
    if (fvmers_3pbdf == ImplicitData::BDF_SCHEME2) {
      DistVec<int> tempInt(fluidId);
      tempInt = 0;
      if (distLSS) {
	distLSS->getSwept(tempInt);
      }    
      DistVec<int> minus1(fluidId);
      minus1 = -1;
      if(!data->exist_nm1) riemann->updatePhaseChange(*Vn,fluidId,minus1);
      
      varFcn->conservativeToPrimitive(*Un, *Unm1, fluidIdnm1);

      if (checkForRapidlyChangingValues) {
        if (checkForRapidlyChangingPressure > 0.0)
          DistVectorOp::Op(*Un, Qtilde,*firstOrderNodes, *fluidIdnm1, 
                           SetFirstOrderNodes(varFcn,checkForRapidlyChangingPressure,errorHandler)); 

        if (checkForRapidlyChangingDensity > 0.0)
          DistVectorOp::Op(*Un, Qtilde,*firstOrderNodes, *fluidIdnm1, 
                           SetFirstOrderNodes(varFcn,checkForRapidlyChangingDensity,errorHandler,1)); 
      }

      int numFirstOrderNodes = firstOrderNodes->sum(); 
      if (numFirstOrderNodes > 0)
        this->domain->getCommunicator()->fprintf(stdout,"%d nodes set to first order accuracy\n",numFirstOrderNodes);

      DistVectorOp::Op(*Vn,*Unm1, fluidId, *fluidIdnm1, tempInt,MultiphaseRiemannCopy(dim) );
      varFcn->primitiveToConservative(*Unm1,*Vn,&fluidId);
      *Unm1 = *Vn;
      
      riemann->updatePhaseChange(*Vn,fluidId,minus1);
      *Un = Q;
    } else {
      *Unm1 = *Vn;
      *Vn = Q;
      if (checkForRapidlyChangingValues) {
        if (checkForRapidlyChangingPressure > 0.0)
          DistVectorOp::Op(*Un, Qtilde,*firstOrderNodes, *fluidIdnm1, 
                           SetFirstOrderNodes(varFcn,checkForRapidlyChangingPressure,errorHandler)); 

        if (checkForRapidlyChangingDensity > 0.0)
          DistVectorOp::Op(*Un, Qtilde,*firstOrderNodes, *fluidIdnm1, 
                           SetFirstOrderNodes(varFcn,checkForRapidlyChangingDensity,errorHandler,1)); 
      }
      double tau = data->getTauN();
      double beta = (1.0+2.0*tau)/((1.0+tau)*(1.0+tau));
      *Un = beta*Q+(1.0-beta)*Qtilde;
    }

    data->exist_nm1 = true;
      
  } else {
    *Vn = *Un = Q;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::computeBar(bool doInitialTasks, DistMacroCellSet *macroCells,
                                    DistGeoState &geoState, int scopeDepth)
{
                                                                                                                          
  if (doInitialTasks) {
    if (data->use_nm1){
      varFcn->conservativeToPrimitive(*Unm1, *Vn);
      domain->computeVBar(macroCells, doInitialTasks, geoState, *VnBar, *Vn, scopeDepth, 2);
      varFcn->primitiveToConservative(*VnBar, *Unm1Bar);
    }
    if (data->use_nm2){
      varFcn->conservativeToPrimitive(*Unm2, *Vn);
      domain->computeVBar(macroCells, doInitialTasks, geoState, *VnBar, *Vn, scopeDepth, 3);
      varFcn->primitiveToConservative(*VnBar, *Unm2Bar);
    }
  }
  else {
    if (data->use_nm2) *Unm2Bar = *Unm1Bar;
    if (data->use_nm1) *Unm1Bar = *UnBar;
  }
                                                                                                                          
  varFcn->conservativeToPrimitive(*Un, *Vn);
  domain->computeVBar(macroCells, doInitialTasks, geoState, *VnBar, *Vn, scopeDepth, 1);
  varFcn->primitiveToConservative(*VnBar, *UnBar);
                                                                                                                          
}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::writeToDisk(char *name)
{

  domain->writeVectorToFile(name, 0, 0.0, *Un);
  if (data->use_nm1) 
    domain->writeVectorToFile(name, 1, 0.0, *Unm1);
  if (data->use_nm2) 
    domain->writeVectorToFile(name, 2, 0.0, *Unm2);

}

//------------------------------------------------------------------------------

template<int dim>
void DistTimeState<dim>::get_dW_dt(bool doInitialTasks,
                                   DistGeoState &geoState,
                                   DistVec<double> &ctrlVol,
                                   DistSVec<double,dim> &Q,
                                   DistSVec<double,dim> &VBar,
                                   DistSVec<double,dim> &dWdt,
                                   DistSVec<double,dim> &dWBardt,
                                   DistMacroCellSet *macroCells,
                                   DistSVec<double,1> **volRatio,
                                   int scopeDepth)
{
                                                                                                                          
  DistSVec<double,dim>* Sigma = new DistSVec<double,dim>(domain->getNodeDistInfo());
                                                                                                                          
  *Sigma = 0.0;
                                                                                                                          
  varFcn->primitiveToConservative(VBar, *QBar);
                                                                                                                          
  // compute dWdt //
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->get_dW_dt(Q.getMasterFlag(iSub), geoState(iSub),
                                   ctrlVol(iSub), Q(iSub), dWdt(iSub));
                                                                                                                          
  // compute dWBardt //
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subTimeState[iSub]->get_dWBar_dt((*QBar).getMasterFlag(iSub), geoState(iSub),
                                      ctrlVol(iSub), (*QBar)(iSub), (*UnBar)(iSub),
                                      (*Unm1Bar)(iSub), (*Unm2Bar)(iSub), (*Sigma)(iSub));
                                                                                                                          
  domain->assemble_dWdt(dWdt, *Sigma);
                                                                                                                          
  domain->computedWBar_dt(dWBardt, *Sigma, macroCells, volRatio, scopeDepth);
                                                                                                                          
                                                                                                                          
  delete (Sigma);
                                                                                                                          
}
                                                                                                                          
//------------------------------------------------------------------------------

// Included (MB)
template<int dim> DistVec<double>*
DistTimeState<dim>::getDerivativeOfInvReynolds(DistGeoState &geoState,
		DistSVec<double,3> &X, DistSVec<double,3> &dX, DistVec<double> &ctrlVol,
		DistVec<double> &dCtrlVol, DistSVec<double,dim> &V, DistSVec<double,dim>
		&dV, double dMach)
{

//Remark: Error mesage for pointers
  if (dIdti == 0) {
    fprintf(stderr, "*** Error: Variable dIdti does not exist!\n");
    exit(1);
  }
  if (dIdtv == 0) {
    fprintf(stderr, "*** Error: Variable dIdtv does not exist!\n");
    exit(1);
  }
  if (dIrey == 0) {
    fprintf(stderr, "*** Error: Variable dIrey does not exist!\n");
    exit(1);
  }

  //domain->computeDerivativeOfInvReynolds(fet, varFcn, geoState, X, dX, ctrlVol, dCtrlVol, V, dV, *idti, *dIdti, *idtv, *dIdtv, *dIrey, dMach, betav, beta, dbeta, k1, cmach);
  domain->computeDerivativeOfInvReynolds(fet, varFcn, geoState, X, dX, ctrlVol, dCtrlVol, V, dV, *idti, *dIdti, *idtv, *dIdtv, *dIrey, dMach, tprec, sprec);

  return dIrey;

}

//------------------------------------------------------------------------------

template<int dim> 
double DistTimeState<dim>::getNewtonTag() const {
  return *(domain->getNewtonTag()); 
}

//------------------------------------------------------------------------------

template<int dim>
int DistTimeState<dim>::getNewtonStateStep() const {
  return *(domain->getNewtonStateStep());
}

//------------------------------------------------------------------------------

template<int dim>
int DistTimeState<dim>::getNewtonResidualStep() const {
  return *(domain->getNewtonResidualStep());
}

//------------------------------------------------------------------------------

template<int dim>
int DistTimeState<dim>::getKrylovStep() const {
  return *(domain->getKrylovStep());
}

//------------------------------------------------------------------------------

template<int dim> 
void DistTimeState<dim>::setExistsNm1() {

  data->exist_nm1 = true;
}

//------------------------------------------------------------------------------

template<int dim> 
void DistTimeState<dim>::setDtNm1(double dt) {

  data->dt_nm1 = dt;
}

