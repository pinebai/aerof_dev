#include "OneDimensionalSolver.h"

#include <fstream>
#include <iostream>
using namespace std;
#include <cmath>

#include "FluxFcn.h"
#include "VarFcn.h"
#include "LocalRiemannDesc.h"
#include "FluidSelector.h"
#include "OneDimensionalSourceTerm.h"
#include "OneDimensionalInterpolator.h"
#include "RefVal.h"

#include "IoData.h"
#include "Domain.h"

//------------------------------------------------------------------------------
OneDimensional::OneDimensional(int np,double* mesh,IoData &ioData, Domain *domain) : 
  numPoints(np),
  ctrlVol(numPoints), ctrlSurf(numPoints+1), X(numPoints), Y(numPoints+1), cutCellStatus(numPoints+1),
  U(numPoints), V(numPoints), R(numPoints),
  gradV(numPoints), Phi(numPoints), Rphi(numPoints), gradPhi(numPoints),
  fluidId(numPoints),fluidIdn(numPoints),
  fluidSelector(ioData.eqs.numPhase,ioData), refVal(ioData.ref.rv),
  Wr(numPoints),Vslope(numPoints),Phislope(numPoints),
  rupdate(numPoints), weight(numPoints), interfacialWi(numPoints),
  interfacialWj(numPoints), riemannStatus(numPoints), Phin(numPoints),
  programmedBurn(NULL), fidToSet(numPoints),lastPhaseChangeValue(numPoints)
{
  // equation modelling
  coordType  = ioData.oneDimensionalInfo.coordType;
  volumeType = ioData.oneDimensionalInfo.volumeType;
  interfaceTreatment = 0;

  // time and space domain definition
  maxDistance = mesh[np-1];
  finalTime = ioData.ts.maxTime;
  cfl = ioData.ts.cfl.cfl0;

  // Copy 1D mesh to X
  X = 0.0;
  for (int i = 0; i < np; ++i)
    X[i][0] = mesh[i];
  delete [] mesh;

  bubbleRadiusFile = new char[256];
 
  if (ioData.oneDimensionalInfo.problemMode == OneDimensionalInfo::MULTIFLUID)
    problemMode = MultiFluid;
  else
    problemMode = FSI; 


  // output
  frequency = ioData.output.transient.frequency;

  // necessary for computation: varFcn to compute different state quantities
  //                            fluxFcn to compute fluxes at interfaces
  //                            riemann to solve riemann problem between two fluids
  //                            source to compute source term(s) for spherical and cylindrical problems
  varFcn = new VarFcn(ioData);

  fluxFcn = new FluxFcn *[3];
  fluxFcn[0] = new FluxFcn(0,BC_INTERNAL,ioData,varFcn);
  fluxFcn[1] = new FluxFcn(0,BC_SYMMETRY,ioData,varFcn);
  fluxFcn[2] = new FluxFcn(0,BC_OUTLET_FIXED,ioData,varFcn);

  //riemann = new LocalRiemannGfmparGasJWL(varFcn,0,1,0,MultiFluidData::RK2);
  //riemann = new LocalRiemannGfmpGasJWL(varFcn,0,1);

  source = 0;
  if(volumeType == OneDimensionalInfo::CONSTANT_VOLUME){
    if(coordType == OneDimensionalInfo::CYLINDRICAL ||
       coordType == OneDimensionalInfo::SPHERICAL)
      source = new OneDimensionalSourceTerm();
  }else if(volumeType == OneDimensionalInfo::REAL_VOLUME){
  }

  strcpy(bubbleRadiusFile,"");

  recFcn = createRecFcn(ioData);
  recFcnLS = createRecFcnLS(ioData);

  if (ioData.ts.type != TsData::EXPLICIT) {
    fprintf(stderr,"Only explicit integration available for the 1D solver!\n");
    exit(1);
  }

  switch (ioData.ts.expl.type) {

  case ExplicitData::RUNGE_KUTTA_4:
    Vintegrator = new RKIntegrator< SVec<double,5> >(RKIntegrator< SVec<double,5> >::RK4, numPoints);
    Phiintegrator = new RKIntegrator< SVec<double,1> >(RKIntegrator< SVec<double,1> >::RK4, numPoints);
    break;
  case ExplicitData::RUNGE_KUTTA_2:
    Vintegrator = new RKIntegrator< SVec<double,5> >(RKIntegrator< SVec<double,5> >::RK2, numPoints);
    Phiintegrator = new RKIntegrator< SVec<double,1> >(RKIntegrator< SVec<double,1> >::RK2, numPoints);
    break;
  case ExplicitData::FORWARD_EULER:
    Vintegrator = new RKIntegrator< SVec<double,5> >(RKIntegrator< SVec<double,5> >::FE, numPoints);
    Phiintegrator = new RKIntegrator< SVec<double,1> >(RKIntegrator< SVec<double,1> >::FE, numPoints);
    break;

  default:
    fprintf(stderr,"Unavailable explicit integration for the 1D solver!\n");
    exit(1);
    break;
  }

  loadSparseGrid(ioData);

  fidToSet = 0;

  riemann = new ExactRiemannSolver<5>(ioData,rupdate,weight, interfacialWi,
				      interfacialWj, varFcn,
				      tabulationC, fidToSet);

  if (ioData.oneDimensionalInfo.programmedBurn.unburnedEOS >= 0) {
    programmedBurn = new ProgrammedBurn(ioData,&this->X);
    this->fluidSelector.attachProgrammedBurn(programmedBurn);
  }
  programmedBurnStopPercentDistance = ioData.ts.programmedBurnShockSensor ;
  if ( ioData.oneDimensionalInfo.programmedBurn.unburnedEOS >= 0 ){
    programmedBurnIsUsed = true ;
  }
  else{
    programmedBurnIsUsed = false ;
  }

  setupOutputFiles(ioData);
  setupProbes(ioData);

  setupFixes(ioData);

  cutCellStatus = 0;

  if (ioData.schemes.ns.dissipation == SchemeData::SIXTH_ORDER) {
    isSixthOrder = true;
  } else
    isSixthOrder = false;

  beta = ioData.schemes.ns.beta;

  if (ioData.mf.interfaceTreatment == MultiFluidData::SECONDORDER)
    interfaceTreatment = 1;

  levelSetMethod = 0;
 
  interfaceExtrapolation = 0;
  if (ioData.mf.interfaceExtrapolation == MultiFluidData::EXTRAPOLATIONSECONDORDER)
    interfaceExtrapolation = 1;


  if (ioData.mf.interfaceLimiter == MultiFluidData::LIMITERALEX1) {

    limiterLeft = limiterRight = 2; 
  } else {

    limiterLeft = limiterRight = 1; 
  }

  //limiterLeft = 0;
  
  
  if (ioData.mf.levelSetMethod == MultiFluidData::HJWENO)
    levelSetMethod = 1;
  else if (ioData.mf.levelSetMethod == MultiFluidData::SCALAR/* || ioData.mf.levelSetMethod == MultiFluidData::CONSERVATIVE*/)
    levelSetMethod = 2;
  else if (ioData.mf.levelSetMethod == MultiFluidData::PRIMITIVE)
    levelSetMethod = 3;
 
  int sto = ioData.oneDimensionalInfo.sourceTermOrder; 
  if(volumeType == OneDimensionalInfo::CONSTANT_VOLUME){
    if(coordType == OneDimensionalInfo::CYLINDRICAL)
      source->initialize(1.0,sto,X);
    else if (coordType == OneDimensionalInfo::SPHERICAL)
      source->initialize(2.0,sto,X);
  }

  typePhaseChange = ioData.mf.typePhaseChange;

  myTimer = domain->getTimer();

}
//------------------------------------------------------------------------------
OneDimensional::~OneDimensional(){

  delete varFcn;
  delete riemann;
  for(int i=0; i<3; i++) delete fluxFcn[i];
  delete [] fluxFcn;
  
  if (source)
    delete source;

  if (tabulationC)
    delete tabulationC;

}

void OneDimensional::setupFixes(IoData& ioData) {


  double spheres[20][4];
  double boxes[20][2][3];
  double cones[20][2][4];
  int j, nspheres = 0, nboxes = 0, ncones = 0;
  for (j=0; j<ioData.schemes.fixes.num; ++j) {
    if (ioData.schemes.fixes.spheres[j]->r > 0.0) {
      spheres[nspheres][0] = ioData.schemes.fixes.spheres[j]->x0;
      spheres[nspheres][1] = ioData.schemes.fixes.spheres[j]->y0;
      spheres[nspheres][2] = ioData.schemes.fixes.spheres[j]->z0;
      spheres[nspheres][3] = ioData.schemes.fixes.spheres[j]->r;
      ++nspheres;
      if (ioData.schemes.fixes.symmetry == SchemeFixData::X) {
	spheres[nspheres][0] = - ioData.schemes.fixes.spheres[j]->x0;
	spheres[nspheres][1] = ioData.schemes.fixes.spheres[j]->y0;
	spheres[nspheres][2] = ioData.schemes.fixes.spheres[j]->z0;
	spheres[nspheres][3] = ioData.schemes.fixes.spheres[j]->r;
	++nspheres;
      }
      else if (ioData.schemes.fixes.symmetry == SchemeFixData::Y) {
	spheres[nspheres][0] = ioData.schemes.fixes.spheres[j]->x0;
	spheres[nspheres][1] = - ioData.schemes.fixes.spheres[j]->y0;
	spheres[nspheres][2] = ioData.schemes.fixes.spheres[j]->z0;
	spheres[nspheres][3] = ioData.schemes.fixes.spheres[j]->r;
	++nspheres;
      }
      else if (ioData.schemes.fixes.symmetry == SchemeFixData::Z) {
	spheres[nspheres][0] = ioData.schemes.fixes.spheres[j]->x0;
	spheres[nspheres][1] = ioData.schemes.fixes.spheres[j]->y0;
	spheres[nspheres][2] = - ioData.schemes.fixes.spheres[j]->z0;
	spheres[nspheres][3] = ioData.schemes.fixes.spheres[j]->r;
	++nspheres;
      }
    }
    if (ioData.schemes.fixes.boxes[j]->x0 < ioData.schemes.fixes.boxes[j]->x1) {
      boxes[nboxes][0][0] = ioData.schemes.fixes.boxes[j]->x0;
      boxes[nboxes][0][1] = ioData.schemes.fixes.boxes[j]->y0;
      boxes[nboxes][0][2] = ioData.schemes.fixes.boxes[j]->z0;
      boxes[nboxes][1][0] = ioData.schemes.fixes.boxes[j]->x1;
      boxes[nboxes][1][1] = ioData.schemes.fixes.boxes[j]->y1;
      boxes[nboxes][1][2] = ioData.schemes.fixes.boxes[j]->z1;
      ++nboxes;
      if (ioData.schemes.fixes.symmetry == SchemeFixData::X) {
	boxes[nboxes][0][0] = -ioData.schemes.fixes.boxes[j]->x1;
	boxes[nboxes][0][1] = ioData.schemes.fixes.boxes[j]->y0;
	boxes[nboxes][0][2] = ioData.schemes.fixes.boxes[j]->z0;
	boxes[nboxes][1][0] = -ioData.schemes.fixes.boxes[j]->x0;
	boxes[nboxes][1][1] = ioData.schemes.fixes.boxes[j]->y1;
	boxes[nboxes][1][2] = ioData.schemes.fixes.boxes[j]->z1;
	++nboxes;
      }
      if (ioData.schemes.fixes.symmetry == SchemeFixData::Y) {
	boxes[nboxes][0][0] = ioData.schemes.fixes.boxes[j]->x0;
	boxes[nboxes][0][1] = -ioData.schemes.fixes.boxes[j]->y1;
	boxes[nboxes][0][2] = ioData.schemes.fixes.boxes[j]->z0;
	boxes[nboxes][1][0] = ioData.schemes.fixes.boxes[j]->x1;
	boxes[nboxes][1][1] = -ioData.schemes.fixes.boxes[j]->y0;
	boxes[nboxes][1][2] = ioData.schemes.fixes.boxes[j]->z1;
	++nboxes;
      }
     if (ioData.schemes.fixes.symmetry == SchemeFixData::Z) {
	boxes[nboxes][0][0] = ioData.schemes.fixes.boxes[j]->x0;
	boxes[nboxes][0][1] = ioData.schemes.fixes.boxes[j]->y0;
	boxes[nboxes][0][2] = -ioData.schemes.fixes.boxes[j]->z1;
	boxes[nboxes][1][0] = ioData.schemes.fixes.boxes[j]->x1;
	boxes[nboxes][1][1] = ioData.schemes.fixes.boxes[j]->y1;
	boxes[nboxes][1][2] = -ioData.schemes.fixes.boxes[j]->z0;
	++nboxes;
      }
    }
    if (ioData.schemes.fixes.cones[j]->r0 >= 0.0 && ioData.schemes.fixes.cones[j]->r1 >= 0.0) {
      cones[ncones][0][0] = ioData.schemes.fixes.cones[j]->x0;
      cones[ncones][0][1] = ioData.schemes.fixes.cones[j]->y0;
      cones[ncones][0][2] = ioData.schemes.fixes.cones[j]->z0;
      cones[ncones][0][3] = ioData.schemes.fixes.cones[j]->r0;
      cones[ncones][1][0] = ioData.schemes.fixes.cones[j]->x1;
      cones[ncones][1][1] = ioData.schemes.fixes.cones[j]->y1;
      cones[ncones][1][2] = ioData.schemes.fixes.cones[j]->z1;
      cones[ncones][1][3] = ioData.schemes.fixes.cones[j]->r1;
      ++ncones;
      if (ioData.schemes.fixes.symmetry == SchemeFixData::X) {
        cones[ncones][0][0] = -ioData.schemes.fixes.cones[j]->x0;
        cones[ncones][0][1] = ioData.schemes.fixes.cones[j]->y0;
        cones[ncones][0][2] = ioData.schemes.fixes.cones[j]->z0;
        cones[ncones][0][3] = ioData.schemes.fixes.cones[j]->r0;
        cones[ncones][1][0] = -ioData.schemes.fixes.cones[j]->x1;
        cones[ncones][1][1] = ioData.schemes.fixes.cones[j]->y1;
        cones[ncones][1][2] = ioData.schemes.fixes.cones[j]->z1;
        cones[ncones][1][3] = ioData.schemes.fixes.cones[j]->r1;
        ++ncones;
      }
      if (ioData.schemes.fixes.symmetry == SchemeFixData::Y) {
        cones[ncones][0][0] = ioData.schemes.fixes.cones[j]->x0;
        cones[ncones][0][1] = -ioData.schemes.fixes.cones[j]->y0;
        cones[ncones][0][2] = ioData.schemes.fixes.cones[j]->z0;
        cones[ncones][0][3] = ioData.schemes.fixes.cones[j]->r0;
        cones[ncones][1][0] = ioData.schemes.fixes.cones[j]->x1;
        cones[ncones][1][1] = -ioData.schemes.fixes.cones[j]->y1;
        cones[ncones][1][2] = ioData.schemes.fixes.cones[j]->z1;
        cones[ncones][1][3] = ioData.schemes.fixes.cones[j]->r1;
        ++ncones;
      }
      if (ioData.schemes.fixes.symmetry == SchemeFixData::Z) {
        cones[ncones][0][0] = ioData.schemes.fixes.cones[j]->x0;
        cones[ncones][0][1] = ioData.schemes.fixes.cones[j]->y0;
        cones[ncones][0][2] = -ioData.schemes.fixes.cones[j]->z0;
        cones[ncones][0][3] = ioData.schemes.fixes.cones[j]->r0;
        cones[ncones][1][0] = ioData.schemes.fixes.cones[j]->x1;
        cones[ncones][1][1] = ioData.schemes.fixes.cones[j]->y1;
        cones[ncones][1][2] = -ioData.schemes.fixes.cones[j]->z1;
        cones[ncones][1][3] = ioData.schemes.fixes.cones[j]->r1;
        ++ncones;
      }
    }
  }

  loctag = new int[numPoints];
  memset(loctag,0,sizeof(int)*numPoints);

  if (nspheres > 0 || nboxes > 0 || ncones > 0) {
    for (j=0; j<nspheres; ++j)
      printf( "*** Warning: set the gradients to zero in [(%g, %g, %g), %g]\n",
		  spheres[j][0], spheres[j][1], spheres[j][2], spheres[j][3]);
    for (j=0; j<nboxes; ++j)
      printf( "*** Warning: set the gradients to zero in [(%g, %g, %g), (%g, %g, %g)]\n",
		  boxes[j][0][0], boxes[j][0][1], boxes[j][0][2],
		  boxes[j][1][0], boxes[j][1][1], boxes[j][1][2]);

    for (j=0; j<ncones; ++j)
      printf( "*** Warning: set the gradients to zero in cone [(%g, %g, %g), %g; (%g, %g, %g), %g]\n",
                  cones[j][0][0], cones[j][0][1], cones[j][0][2], cones[j][0][3],
                  cones[j][1][0], cones[j][1][1], cones[j][1][2], cones[j][1][3]);

    for (int i=0; i<numPoints; ++i) {
      double x0[3] = {X[i][0],0.0,0.0};
      loctag[i] = 0;
      for (j=0; j<nspheres; ++j) {
	double r = sqrt( (x0[0] - spheres[j][0])*(x0[0] - spheres[j][0]) +
			 (x0[1] - spheres[j][1])*(x0[1] - spheres[j][1]) +
			 (x0[2] - spheres[j][2])*(x0[2] - spheres[j][2]) );
	if (r <= spheres[j][3])
	  loctag[i] = 1;
      }
      for (j=0; j<nboxes; ++j) {
	if ((x0[0] >= boxes[j][0][0]) && (x0[0] <= boxes[j][1][0]) &&
	    (x0[1] >= boxes[j][0][1]) && (x0[1] <= boxes[j][1][1]) &&
	    (x0[2] >= boxes[j][0][2]) && (x0[2] <= boxes[j][1][2]))
	    loctag[i] = 1;
      }
      for (j=0; j<ncones; ++j)  {
	Vec3D dr(cones[j][1][0]-cones[j][0][0], cones[j][1][1]-cones[j][0][1], cones[j][1][2]-cones[j][0][2]);
	double height = dr.norm();
	dr /= height;
	Vec3D xp;
	Vec3D pr0(x0[0]-cones[j][0][0], x0[1]-cones[j][0][1], x0[2]-cones[j][0][2]);
	double h = pr0*dr;
	if (h >= 0.0 && h <= height)  {
	  xp = pr0 - (h*dr);
	  double r = cones[j][0][3] + (cones[j][1][3]-cones[j][0][3]) * h / height;
	  if (xp.norm() < r)
	    loctag[i] = 1;
          }
      }
    }
  }
}

void OneDimensional::setupOutputFiles(IoData& iod) {

  int i;


  int sp = strlen(iod.output.transient.prefix) + 1;
  int spr = strlen(iod.output.restart.prefix) + 1;


  if (iod.output.restart.solutions[0] != 0) {
    outfile = new char[spr + strlen(iod.output.restart.solutions)];
    sprintf(outfile, "%s%s", 
	    iod.output.restart.prefix, iod.output.restart.solutions);
  }

  for (i=0; i<PostFcn::SSIZE; ++i) {
    sscale[i] = 1.0;
    scalars[i] = 0;
  }

  for (i=0; i<PostFcn::VSIZE; ++i) {
    vscale[i] = 1.0;
    vectors[i] = 0;
  }

  if (iod.output.transient.density[0] != 0) {
    sscale[PostFcn::DENSITY] = iod.ref.rv.density;
    scalars[PostFcn::DENSITY] = new char[sp + strlen(iod.output.transient.density)];
    sprintf(scalars[PostFcn::DENSITY], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.density);
  }
  if (iod.output.transient.mach[0] != 0) {
    scalars[PostFcn::MACH] = new char[sp + strlen(iod.output.transient.mach)];
    sprintf(scalars[PostFcn::MACH], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.mach);
  }
  if (iod.output.transient.pressure[0] != 0) {
    sscale[PostFcn::PRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::PRESSURE] = new char[sp + strlen(iod.output.transient.pressure)];
    sprintf(scalars[PostFcn::PRESSURE], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.pressure);
  }
  if (iod.output.transient.temperature[0] != 0) {
    sscale[PostFcn::TEMPERATURE] = iod.ref.rv.temperature;
    scalars[PostFcn::TEMPERATURE] = new char[sp + strlen(iod.output.transient.temperature)];
    sprintf(scalars[PostFcn::TEMPERATURE], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.temperature);
  }
  if (iod.output.transient.philevel[0] != 0) {
    sscale[PostFcn::PHILEVEL] = 1.0;
    scalars[PostFcn::PHILEVEL] = new char[sp + strlen(iod.output.transient.philevel)];
    sprintf(scalars[PostFcn::PHILEVEL], "%s%s",
            iod.output.transient.prefix, iod.output.transient.philevel);
  }
  if (iod.output.transient.fluidid[0] != 0) {
    sscale[PostFcn::FLUIDID] = 1.0;
    scalars[PostFcn::FLUIDID] = new char[sp + strlen(iod.output.transient.fluidid)];
    sprintf(scalars[PostFcn::FLUIDID], "%s%s",
            iod.output.transient.prefix, iod.output.transient.fluidid);
  }
  if (iod.output.transient.velocity[0] != 0) {
    vscale[PostFcn::VELOCITY] = iod.ref.rv.velocity;
    vectors[PostFcn::VELOCITY] = new char[sp + strlen(iod.output.transient.velocity)];
    sprintf(vectors[PostFcn::VELOCITY], "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.velocity);
  }
  if (iod.output.transient.bubbleRadius[0] != 0) {
    bubbleRadiusFile = new char[sp + strlen(iod.output.transient.bubbleRadius)];
    sprintf(bubbleRadiusFile, "%s%s", 
	    iod.output.transient.prefix, iod.output.transient.bubbleRadius);
  }
}

//------------------------------------------------------------------------------
void OneDimensional::load1DMesh(IoData& ioData,int& numPts,double* &meshPoints) {

  if (ioData.input.geometry[0] != 0) {
    char mesh1d[256];
    
    sprintf(mesh1d,"%s%s",ioData.input.prefix,ioData.input.geometry);
    FILE* fin = fopen(mesh1d,"r");
    
    if (!fin) {

      std::cout << "*** Error: cannot open mesh file " << mesh1d << std::endl;
      exit(-1);
    }
    
    int n = fscanf(fin, "%i",&numPts);
    meshPoints = new double[numPts];
    for (int i = 0; i < numPts; ++i) {

      int m = fscanf(fin,"%lf",&meshPoints[i]);
    }
  } else {

    numPts = ioData.oneDimensionalInfo.numPoints;
    meshPoints = new double[numPts];
    for (int i = 0; i < numPts; ++i)
      meshPoints[i] = (double)i / (numPts-1)*ioData.oneDimensionalInfo.maxDistance;
  }
}
//------------------------------------------------------------------------------
void OneDimensional::spatialSetup(){

  Y = 0.0;
  Y[0][0] = X[0][0];
  for(int i=1; i<numPoints; i++) Y[i][0] = 0.5*(X[i-1][0]+X[i][0]);
  Y[numPoints][0] = X[numPoints-1][0];

  ctrlVol[0][0] = 0.5*(X[1][0]-X[0][0]);
  for(int i=1; i<numPoints-1; i++)
    ctrlVol[i][0] = 0.5*(X[i+1][0]-X[i-1][0]);
  ctrlVol[numPoints-1][0] = 0.5*(X[numPoints-1][0]-X[numPoints-2][0]);

  for(int i=0; i<numPoints+1; i++) ctrlSurf[i][0] = 1.0;

  // in case a real Finite Volume approach is considered
  // with real cylindrical/spherical volumes, the control volumes and
  // surfaces must be computed differently as follows:
  if(volumeType == OneDimensionalInfo::REAL_VOLUME){

    if(coordType == OneDimensionalInfo::SPHERICAL){

      for(int i=0; i<numPoints; i++)
        ctrlVol[i][0] = (Y[i+1][0]-Y[i][0])*
                        (Y[i+1][0]*Y[i+1][0]+Y[i+1][0]*Y[i][0]+Y[i][0]*Y[i][0])/3.0;
      for(int i=0; i<numPoints+1; i++) ctrlSurf[i][0] = Y[i][0]*Y[i][0];

    }else if(coordType == OneDimensionalInfo::CYLINDRICAL){

      for(int i=0; i<numPoints; i++)
        ctrlVol[i][0] = 0.5*(Y[i+1][0]-Y[i][0])*(Y[i+1][0]+Y[i][0]);
      for(int i=0; i<numPoints+1; i++) ctrlSurf[i][0] = Y[i][0];
    }
  }

   
}
//------------------------------------------------------------------------------
void OneDimensional::temporalSetup(){

}
//------------------------------------------------------------------------------
void OneDimensional::stateInitialization(OneDimensionalInfo &data){
  
  V = 0.0;
  // initialize V
  if (data.mode == OneDimensionalInfo::CONVTEST1) {
    data.pressure1 = data.pressure2;
  }
  interfaceLocation = data.interfacePosition;
  for(int i=0; i<numPoints; i++){
    if(X[i][0]<data.interfacePosition){
      if (varFcn->getType(1) != VarFcnBase::TAIT) {
	V[i][0] = data.density1;
	V[i][1] = data.velocity1;
	V[i][4] = data.pressure1;
      } else {
	V[i][0] = data.density1;
	V[i][1] = data.velocity1;
	V[i][4] = data.temperature1;
      }
    }else{
      if (varFcn->getType(0) != VarFcnBase::TAIT) {
	V[i][0] = data.density2;
	V[i][1] = data.velocity2;
	V[i][4] = data.pressure2;
      } else {
	V[i][0] = data.density2;
	V[i][1] = data.velocity2;
	V[i][4] = data.temperature2;
      }
    }

    if (data.mode == OneDimensionalInfo::CONVTEST1) {
      V[i][0] = data.density2;
      if (fabs(X[i][0]-data.interfacePosition) < 0.2)
	V[i][4] = 
	  (1000000.0*pow((X[i][0]-data.interfacePosition)+0.2,4.0)*pow((X[i][0]-data.interfacePosition)-0.2,4.0)+1.0)*data.pressure2;
      else
	V[i][4] = data.pressure2;
      V[i][1] = 0.0;
    }

    if (data.mode == OneDimensionalInfo::CONVTEST2) {
      V[i][0] = data.density2;
      V[i][4] = (atan(-(X[i][0]-data.interfacePosition)*20.0)+3.14159265358979323846*0.5+1.0)*data.pressure2;
      V[i][1] = 0.0;
    }

  }

  isSinglePhase = varFcn->getVarFcnBase(0)->equal(varFcn->getVarFcnBase(1));

  // initialize Phi
  if (levelSetMethod == 0) {
    for(int i=0; i<numPoints; i++)
      Phi[i][0]  = -V[i][0]*(X[i][0]-data.interfacePosition);
  } else {
    for(int i=0; i<numPoints; i++)
      Phi[i][0]  = -(X[i][0]-data.interfacePosition);
  }

  // initialize fluidId and U
  fluidSelector.getFluidId(fluidId,Phi);
  varFcn->primitiveToConservative(V,U,&fluidId);

  fluidIdn = fluidId;

  // compute boundary states
  double temp0[5] = {data.density1, data.velocity1, 0.0, 0.0, data.pressure1};
  double temp1[5] = {data.density2, data.velocity2, 0.0, 0.0, data.pressure2};
  varFcn->primitiveToConservative(temp0,BC[0],fluidId[0]);
  varFcn->primitiveToConservative(temp1,BC[1],fluidId[numPoints-1]);
  BCphi[0] = Phi[0][0]/V[0][0];
  BCphi[1] = Phi[numPoints-1][0]/V[numPoints-1][0];
  
}
//------------------------------------------------------------------------------
void OneDimensional::totalTimeIntegration(){

  // messages
  SVec<double,1> timeSteps(numPoints);

  time = 0.0;
  double dt   = 0.0;
  int iteration = 0;

  lastPhaseChangeValue = -1;//V;

  resultsOutput(time,iteration);
  double cpu_time = myTimer->getTime();
  while(time<finalTime){
    // determine how much to advance in time
    dt = cfl*computeMaxTimeStep();

    if (programmedBurn){
      programmedBurn->setCurrentTime(time,varFcn, U,fluidId,fluidIdn);
    }
    if(time+dt>finalTime) dt = finalTime-time;
    if(frequency > 0 && iteration % frequency == 0)
      cout <<"*** Iteration " << iteration <<": Time = "<<time*refVal.time<<", and dt = "<<dt*refVal.time<<endl;
    iteration++;

    // advance one iteration
    singleTimeIntegration(dt);
    time += dt;

    // Programmed burn shock sensor:
    // If programmed burn is used, if the shock reaches a specified maximum location, stop the simulation:
    if (programmedBurnIsUsed == true && iteration > 10 ){
      // Obtain the shock sensor node numbers:
      int shockSensorNode3 = (int)( numPoints * programmedBurnStopPercentDistance ) ;
      int shockSensorNode2 = shockSensorNode3 - 1 ;
      int shockSensorNode1 = shockSensorNode2 - 1 ;
      // Obtain the x-velocity at the shock sensor nodes:
      double velocityAtShockSensorNode1 = std::abs(V[shockSensorNode1][1]) ;
      double velocityAtShockSensorNode2 = std::abs(V[shockSensorNode2][1]) ;
      double velocityAtShockSensorNode3 = std::abs(V[shockSensorNode3][1]) ;
      // Calculate the shock sensor value (this is similar to how a flux limiter works):
      double shockSensorValue = 1.0 ;
      if (velocityAtShockSensorNode1 > 0.0 && velocityAtShockSensorNode2 > 0.0 && velocityAtShockSensorNode3 > 0.0 ) {
        shockSensorValue = ( velocityAtShockSensorNode2 - velocityAtShockSensorNode1 )/( velocityAtShockSensorNode3 - velocityAtShockSensorNode2 ) ;
      }
      // The shock sensor value tends to 0 on a strong shock, and tends to 1 otherwise.
      // If the shock sensor is less than a threshold value (0.5), the detonation shock has reached the sensor, and the simulation is stopped:
      if (shockSensorValue < 0.5 ){
        double currentTime = time ;
        finalTime = time ;
        std::cout << "*** The shockwave has reached " << programmedBurnStopPercentDistance*100.0 << " percent of the mesh distance at the time of " << currentTime*refVal.time << " seconds." << std::endl ;
        std::cout << "*** The simulation has stopped automatically." << std::endl ;
      }
    }

    outputProbes(time,iteration-1);

    if(frequency > 0 && iteration % frequency == 0)
      resultsOutput(time,iteration);
  }
  resultsOutput(time,iteration);
  restartOutput(time,iteration);
  std::cout << "Total time: " << myTimer->getTime()-cpu_time << " s" << std::endl;
}
//------------------------------------------------------------------------------
double OneDimensional::computeMaxTimeStep(){

  // very crude CFL law
  double c = 0.0;
  double maxTimeStep = 1e20;
  varFcn->conservativeToPrimitive(U,V,&fluidId);

  for(int i=0; i<numPoints; i++){
    c = varFcn->computeSoundSpeed(V[i],fluidId[i]);
    if (c == 0.0) {
      std::cout << "*** Error: zero sound speed detected!" << std::endl;
      std::cout << "c = " << c <<  " " << fluidId[i] <<  std::endl;
    }
    maxTimeStep = min(maxTimeStep, 0.5*(Y[i+1][0]-Y[i][0])/c);
  }

  return maxTimeStep;
}

void OneDimensional::levelSetDerivative(double t0, Vec<double>& phi, Vec<double>& k) {
  
  double u;
  if (interfaceTreatment == 1) {
   
    /*for(int i=0; i<numPoints-1; i++){
      if (cutCellStatus[i] == 1)
	V[i][1] = 0.5*(V[i-1][1]+V[i+1][1]);
	}*/
    for(int i=0; i<numPoints-1; i++){
      if (cutCellStatus[i] == 1) {
	double xi = phi[0];
	double ul = V[i-1][1]*(xi-X[i-2][0])/(X[i-1][0]-X[i-2][0])+V[i-2][1]*(xi-X[i-1][0])/(X[i-2][0]-X[i-1][0]);
	double ur = V[i+1][1]*(xi-X[i+2][0])/(X[i+1][0]-X[i+2][0])+V[i+2][1]*(xi-X[i+1][0])/(X[i+2][0]-X[i+1][0]);
	k[0] = (ul+ur)*0.5;
      }
      
    }
  } else {
    OneDimensionalInterpolator::Interpolate<5>(V,u , 1,
					       X, phi[0] ,2);
    
    k[0] = u;
    //std::cout << u << std::endl;
  }
}

//------------------------------------------------------------------------------
void OneDimensional::singleTimeIntegration(double dt){
// for now, assume forward Euler

  double Vtemp[5];
  

  fluidSelector.getFluidId(fluidId,Phi);

  if (interfaceTreatment == 1) {
    //varFcn->conservativeToPrimitive(U,V,&fluidId);
    for(int i=0; i<numPoints-1; i++){
      
      if (Y[i][0] < interfaceLocation && Y[i+1][0] >= interfaceLocation) {
	for (int k = 0; k < 5; ++k) {
	  if (interfaceLocation > X[i][0])
	    V[i][k] = (-(X[i][0]-X[i-1][0])*V[i-2][k]+(X[i][0]-X[i-2][0])*V[i-1][k])/(X[i-1][0]-X[i-2][0]);//(2.0*V[i-2][k]-V[i-1][k]);
	  else
	    V[i][k] = (-(X[i][0]-X[i+1][0])*V[i+2][k]+(X[i][0]-X[i+2][0])*V[i+1][k])/(X[i+1][0]-X[i+2][0]);//(2.0*V[i-2][k]-V[i-1][k]);
	}
	cutCellStatus[i] = 1;
      }
      else {
	if (cutCellStatus[i] == 1) {
	  if (interfaceExtrapolation == 0) {
	    int j = i+1;
	    if (fluidId[i-1] == fluidId[i])
	      j = i-1;
	    memcpy(V[i],Wr[j],sizeof(double)*5);
	    //memcpy(V[i],V[j],sizeof(double)*5);
            //std::cout << "\nConstant extrapolated value: " << std::endl;
	    for (int k = 0; k < 5; ++k) {
	      //std::cout << V[i][k] << " ";
	    }
	  }
	  else if (interfaceExtrapolation == 1) {

	    // Linear extrapolation to populate the ghost state
	    int j = i+1;
	    if (fluidId[i-1] == fluidId[i])
	      j = i-1;
            int l = (j-i)*2+i;
	    //double xi = (X[j][0]+X[i][0])*0.5;
	    //std::cout << std::endl << std::endl;	
            double alpha = 1.0;
            if (limiterLeft == 0)
              alpha = 0.0;
            else if (limiterLeft == 2) {
              for (int k = 0; k < 5; ++k) {
                if (k != 2 && k != 3) {
                  alpha = std::min<double>(alpha, fabs(V[j][k]-lastPhaseChangeValue[j][k])/
                                           std::max<double>(1e-8,fabs(V[j][k]-V[l][k])));
                 // alpha = std::min<double>(alpha,fabs(2.0*V[j][1]-V[l][1])/fabs(lastPhaseChangeValue[j][1]));
                }

              }
            }

	    for (int k = 0; k < 5; ++k) {
	      //std::cout << V[j][k] << " ";
	      //std::cout << fabs(V[j][k]-lastPhaseChangeValue[j][k]) << " " <<
              //  std::max<double>(1e-8,fabs(V[j][k]-V[l][k])) << std::endl;
              //alpha = std::min<double>(alpha,fabs(2.0*V[j][1]-V[l][1])/lastPhaseChangeValue[j][1]);
              //std::cout << "alpha = " << alpha << std::endl;
              if (typePhaseChange == MultiFluidData::EXTRAPOLATION) {
  	        V[i][k] = alpha*((X[i][0]-X[j][0])/(X[l][0]-X[j][0])*V[l][k]+
		  (X[i][0]-X[l][0])/(X[j][0]-X[l][0])*V[j][k]) + 
                          (1.0-alpha)*V[j][k];
              } else {
	        V[i][k] = alpha*((X[i][0]-X[j][0])/(interfaceLocation-X[j][0])*Wr[j][k]-
		  (X[i][0]-interfaceLocation)/(interfaceLocation-X[j][0])*V[j][k]) + (1.0-alpha)*Wr[j][k];
              }
	      //std::cout << V[i][k] << " ";
              //std::cout << std::endl;
              lastPhaseChangeValue[i][k] = -1.0;
	    }
	    //memcpy(V[i],V[j],sizeof(double)*5);
	  }
	  //std::cout << std::endl << std::endl;	
	}
	cutCellStatus[i] = 0;
      }
    }
    varFcn->primitiveToConservative(V,U,&fluidId);
  }
  
  riemannStatus = 0;
  //std::cout << "U0 = " << U*U << std::endl;

  Vintegrator->integrate(this,&OneDimensional::EulerF,
			 U,time,dt);
 
  if (levelSetMethod == 1 || levelSetMethod == 0 || levelSetMethod == 3) {
    Phin = Phi;
    Phiintegrator->integrate(this,&OneDimensional::PhiF,
			     Phi,time,dt);
    //Rphi = 0.0;
    //computeLevelSetFluxes();
    
    //preliminary update of U and Phi
    for(int i=0; i<numPoints; i++){
      if (Phi[i][0]*Phin[i][0] < 0.0 &&
	  !riemannStatus[i])
	Phi[i][0] = Phin[i][0];
    }

    // store previous primitive with old fluidId
    varFcn->conservativeToPrimitive(U,V,&fluidId);
    
    // update fluidId
    fluidIdn = fluidId;
    fluidSelector.getFluidId(fluidId,Phi);
  } else {

    fluidIdn = fluidId;

    if (problemMode == MultiFluid) {
      Vec<double> phi(1);
      phi[0] = interfaceLocation;
      varFcn->conservativeToPrimitive(U,V,&fluidId);
      RKIntegrator<Vec<double> > phiI( RKIntegrator< Vec<double> >::RK4, 1);
      phiI.integrate(this, &OneDimensional::levelSetDerivative,
  		     phi,time,dt);
      interfaceLocation = phi[0];
    } else {

      varFcn->conservativeToPrimitive(U,V,&fluidId);
      interfaceLocation = 1.0-pow(time*refVal.time,4.0)/8.0;
    }
    
    for(int i=0; i<numPoints; i++){
      Phi[i][0] = -X[i][0] + interfaceLocation;
    }
    fluidIdn = fluidId;
    fluidSelector.getFluidId(fluidId,Phi);
    
  }

  for(int i=0; i<numPoints; i++){
  
    if (problemMode == FSI) {
      if (fluidId[i] != fluidIdn[i]) {

	memcpy(V[i],V[i+1],sizeof(double)*5);
      }
      continue;
    }
    if (fluidId[i] != fluidIdn[i]&& !varFcn->getVarFcnBase(0)->equal(varFcn->getVarFcnBase(1))) { // Phase change
      
      if (!riemannStatus[i])
	std::cout << "Have a problem!" <<  " " << fluidId[i] << " " << fluidIdn[i] << std::endl;
      if (interfaceTreatment == 0) {
	if (interfaceExtrapolation == 0)
	  memcpy(V[i],Wr[i],sizeof(double)*5);
	else {

	  // Linear extrapolation to populate the ghost state
	  int j = i+1;
	  if (fluidId[i-1] == fluidId[i])
	    j = i-1;
	  double xi = (X[j][0]+X[i][0])*0.5;
	  for (int k = 0; k < 5; ++k) {
	    V[i][k] = (X[i][0]-X[j][0])/(xi-X[j][0])*Wr[i][k]-
	      (X[i][0]-xi)/(xi-X[j][0])*V[j][k];
	  }
	}
      }	
    }
  }

  // initialize Phi
  for(int i=0; i<numPoints-1; i++) {

    if (Phi[i+1][0] < 0.0) {
      if (levelSetMethod == 0)
	interfaceLocation = (X[i+1][0]*Phi[i][0]/V[i][0]-X[i][0]*Phi[i+1][0]/V[i+1][0])/(Phi[i][0]/V[i][0]-Phi[i+1][0]/V[i+1][0]);
      else if (levelSetMethod == 1 || levelSetMethod == 3)
	interfaceLocation = (X[i+1][0]*Phi[i][0]-X[i][0]*Phi[i+1][0])/(Phi[i][0]-Phi[i+1][0]);
      break;
    }
  }

    
  if (levelSetMethod == 0) {
  } else if (levelSetMethod == 1 || levelSetMethod == 3) {
    for(int i=0; i<numPoints-1; i++) {
      Phi[i][0] = (interfaceLocation-X[i][0]);
    }

  }
  //update PhaseChange with new fluidId
  varFcn->primitiveToConservative(V,U,&fluidId);

  // Check solution (clip pressure, that is), if necessary
  if (varFcn->doVerification()) {
    for(int i=0; i<numPoints; i++){
      varFcn->conservativeToPrimitiveVerification(i+1, U[i], Vtemp, fluidId[i]);
    }
  }
   
  if (programmedBurn) {

    programmedBurn->setFluidIds(time, fluidId,U);
  }
}   

void OneDimensional::EulerF(double t, SVec<double,5>& y,SVec<double,5>& k) {

  R = 0.0;
  computeEulerFluxes(y);
  for(int i=0; i<numPoints; i++){
    double* kp = k[i],*Rp = R[i],icv = 1.0/ctrlVol[i][0];
    for(int idim=0; idim<dim; idim++)
      kp[idim] = -Rp[idim] * icv ;
  }

}

class EulerSource {

public:

  EulerSource(VarFcn* _vf,Vec<int>& _fid) :  vf(_vf), fid(_fid) { }

  ~EulerSource() { }

  void compute(double* U, double* f,int i) {

    double V[5];
    vf->conservativeToPrimitive(U,V,fid[i]);
    f[0] = U[1];
    f[1] = U[1]*U[1]/U[0];
    f[4] = (U[4]+vf->getPressure(V, fid[i]))*(U[1]/U[0]);
    f[2] = f[3] = 0.0;
  }

  VarFcn* vf;
  Vec<int>& fid;
};

//------------------------------------------------------------------------------
void OneDimensional::computeEulerFluxes(SVec<double,5>& y){

  double normal[3] = {1.0, 0.0, 0.0};
  double length = 1.0;
  double normalVel = 0.0; // no ALE
  double flux[dim];
  double Udummy[dim],Vtemp[dim];
  int i,j,k;
  
  for(int i=0; i<numPoints; i++){
    if (y[i][0] < 0.0)
      std::cout << "*** Error: node " << i << " has negative density " <<
            y[i][0] << "; fid = " << fluidId[i] << std::endl;
  }

  // Check solution (clip pressure, that is), if necessary
  if (varFcn->doVerification()) {
    for(int i=0; i<numPoints; i++){
      varFcn->conservativeToPrimitiveVerification(i+1, y[i], Vtemp, fluidId[i]);
    }
  } 

  varFcn->conservativeToPrimitive(y,V,&fluidId);
  if (problemMode == MultiFluid)
    computeSlopes(V,Vslope,fluidId,!varFcn->getVarFcnBase(0)->equal(varFcn->getVarFcnBase(1)));
  else
    computeSlopes(V,Vslope,fluidId,false);
  
  for(int iEdge=0; iEdge<numPoints-1; iEdge++){
    i = iEdge;
    j = iEdge+1;
    int fidi = fluidId[i],fidj = fluidId[j];
  
    double Xi = X[i][0],Xj = X[j][0]; 
    double* Vi_ptr = V[i], *Vj_ptr = V[j];

    if (problemMode == FSI && Xi < interfaceLocation && Xj < interfaceLocation)
      continue; 
 
    double Vi[dim*2],Vj[dim*2],Vsi[dim],Vsj[dim],VslopeI[dim],VslopeJ[dim];
    if (!isSixthOrder || !(i > 0 && i < numPoints-3 && ((fidi == fluidId[i+1] && 
				     fidi == fluidId[i+2] && fidi == fluidId[i-1]) || isSinglePhase))) {

      double *Vsi = Vslope[i];      
      double *Vsj = Vslope[j];      
      for (k = 0; k < dim; ++k) {
	VslopeI[k] = Vsi[k];
	VslopeJ[k] = Vsj[k];
      }
    } else {

      double *Vsi = Vslope[i];      
      double *Vsim1 = Vslope[i-1];      
      double *Vsj = Vslope[j];      
      double *Vsjp1 = Vslope[j+1];      
      for (k = 0; k < dim; ++k) {
	VslopeI[k] = (V[j][k]-V[i][k])/(X[j][0]-X[i][0])+(V[i][k]-V[i-1][k])/(X[i][0]-X[i-1][0])+
	  (-1.0/30.0)*((V[j+1][k]-V[j][k])/(X[j+1][0]-X[j][0])-2.0*(V[j][k]-V[i][k])/(X[j][0]-X[i][0])+(V[i][k]-V[i-1][k])/(X[i][0]-X[i-1][0]))/beta+
	  (-2.0/15.0)*(Vsim1[k]-2.0*Vsi[k]+Vsj[k])/beta;
	
	  VslopeJ[k] = (V[j][k]-V[i][k])/(X[j][0]-X[i][0])+(V[j+1][k]-V[j][k])/(X[j+1][0]-X[j][0])+
	  (-1.0/30.0)*((V[j+1][k]-V[j][k])/(X[j+1][0]-X[j][0])-2.0*(V[j][k]-V[i][k])/(X[j][0]-X[i][0])+(V[i][k]-V[i-1][k])/(X[i][0]-X[i-1][0]))/beta+
	  (-2.0/15.0)*(Vsjp1[k]-2.0*Vsj[k]+Vsi[k])/beta;
	
	VslopeI[k] *= 0.5;
	VslopeJ[k] *= 0.5;
	
      }
      
    }

    for (k = 0; k < dim; ++k) {
      Vsi[k] = VslopeI[k]*(Xj-Xi);
      Vsj[k] = VslopeJ[k]*(Xj-Xi);
			  
    }
    
    // One dimensional (source) term
    /*if(volumeType == OneDimensionalInfo::REAL_VOLUME){
      
      if(coordType == OneDimensionalInfo::SPHERICAL) {

	double pi = varFcn->getPressure(Vi,fluidId[i]),pj = varFcn->getPressure(Vj,fluidId[j]);
	R[i][1] -= pi*ctrlSurf[iEdge+1][0];
	R[j][1] += pj*ctrlSurf[iEdge+1][0];
      }
      }*/

    length = (Xj-Xi);

    //if (fidi == fidj || isSinglePhase)
      recFcn->compute(Vi_ptr, Vsi,Vj_ptr, Vsj, Vi, Vj);
    /*else {
      for (k = 0; k < dim; ++k) {
        Vi[k] = Vi_ptr[k];
        Vj[k] = Vj_ptr[k];
      }
    }*/
    
    //std::cout << "Hello" << std::endl;
    varFcn->getVarFcnBase(fidi)->verification(0,Udummy,Vi);
    varFcn->getVarFcnBase(fidj)->verification(0,Udummy,Vj);

    if (interfaceTreatment == 0) {
      
      double* Ri_ptr = R[i],*Rj_ptr = R[j];
      double Cs = ctrlSurf[iEdge+1][0];
       
      if (problemMode == FSI &&
          Xi < interfaceLocation && 
          Xj > interfaceLocation) {
        
        double gradphi[3] = {1.0, 0.0, 0.0};
        double Wi[2*dim], Wj[2*dim];
        Vec3D normalDir(1.0,0.0,0.0);
        double iv[3] = {-0.5*pow(time*refVal.time,3.0)*refVal.time,0.0,0.0};
        for (k = 0; k < dim; ++k) {
          Vi[k+dim] = Vi_ptr[k];
          Vj[k+dim] = Vj_ptr[k];
        }
        riemann->computeFSIRiemannSolution(Vj,iv,normalDir,varFcn,Wi,j,fidi);
  
        fluxFcn[0]->compute(length, 0.0, normal, normalVel, Wi, Vj, flux, fidi);
        for(int k=0; k<dim; ++k) {
          Rj_ptr[k] -= Cs*flux[k];
        }
      }
      else if(fidi == fidj|| isSinglePhase){
  
        fluxFcn[0]->compute(length, 0.0, normal, normalVel, Vi, Vj, flux, fidi);
        for(int k=0; k<dim; ++k) {
          Ri_ptr[k] += Cs*flux[k];
          Rj_ptr[k] -= Cs*flux[k];
        }
      } else{
        double gradphi[3] = {1.0, 0.0, 0.0};
        double Wi[2*dim], Wj[2*dim];
        double Wir[2*dim], Wjr[2*dim];
        double Vir[dim],Vjr[dim];
        double wi, wj;
        double dx[3] = {Xj-Xi, 0.0, 0.0};
        int iteration = 0;
        double fluxi[dim], fluxj[dim];
        
	memcpy(Vir, Vi, sizeof(double)*dim);
	memcpy(Vjr, Vj, sizeof(double)*dim);
        	
/*        if (programmedBurn && fidj == 1) {

          for (int k = 0; k < dim; ++k) {
            Vir[k] = Vi_ptr[k] + Vsi[k]*0.5;
            Vjr[k] = Vj_ptr[k] - Vsj[k]*0.5;
          }
        }
*/
        riemann->computeRiemannSolution(Vir,Vjr,fidi,fidj,gradphi,varFcn,
					Wir,Wjr,i,j,i,dx,0,false);
/*
        if (fidi == 2) {

          std::cout << Vir[0]*refVal.density << " " << Vir[1]*refVal.velocity << " " << Vir[4]*refVal.pressure << std::endl;
          std::cout << Wir[0]*refVal.density << " " << Wir[1]*refVal.velocity << " " << Wir[4]*refVal.pressure << std::endl;
          std::cout << Vjr[0]*refVal.density << " " << Vjr[1]*refVal.velocity << " " << Vjr[4]*refVal.temperature << std::endl;
          std::cout << Wjr[0]*refVal.density << " " << Wjr[1]*refVal.velocity << " " << Wjr[4]*refVal.temperature << std::endl;
        }
*/
	memcpy(Wi, Wir, sizeof(double)*dim);
	memcpy(Wj, Wjr, sizeof(double)*dim);
      
        memcpy(Wr[i], Wj, sizeof(double)*5);
        memcpy(Wr[j], Wi, sizeof(double)*5);
        riemannStatus[i] = riemannStatus[j] = 1;

        //fprintf(stderr,"Edge %d-%d crosses interface\n",i,j);

        fluxFcn[0]->compute(length, 0.0, normal, normalVel, Vi, Wi, fluxi, fidi);
        fluxFcn[0]->compute(length, 0.0, normal, normalVel, Wj, Vj, fluxj, fidj);
        for (int k=0; k<dim; k++){
          Ri_ptr[k] += Cs*fluxi[k];
          Rj_ptr[k] -= Cs*fluxj[k];
        }

      }
    } else if (interfaceTreatment == 1) {

      if(cutCellStatus[i] == 0 && cutCellStatus[j] == 0){

        fluxFcn[0]->compute(length, 0.0, normal, normalVel, Vi, Vj, flux, fidi);
        for(int k=0; k<dim; ++k) {
          R[i][k] += ctrlSurf[iEdge+1][0]*flux[k];
          R[j][k] -= ctrlSurf[iEdge+1][0]*flux[k];
        }
      }else{
        double gradphi[3] = {1.0, 0.0, 0.0};
        double Wi[2*dim], Wj[2*dim];
        double Wir[2*dim], Wjr[2*dim];
        double Vir[dim],Vjr[dim];
        double wi, wj;
        double dx[3] = {X[j][0]-X[i][0], 0.0, 0.0};
        int iteration = 0;
        double fluxi[dim], fluxj[dim];
        int I,J;
        double betapr = 1.0;
        double betapl = 1.0;

        double betam[5] = {0,0,0,0,0};
        if (limiterRight == 0) {
          betapr = 0.0;
          betapl = 0.0;
        }
        if (cutCellStatus[j] == 1) {

          if (interfaceExtrapolation == 1) {
	    //std::cout << "\t" << X[i-1][0] << " " << X[i][0] << " " << X[j][0] << " " << X[j+1][0] << " " << X[j+2][0] << " " << interfaceLocation << std::endl;
	    for (int k = 0; k < dim; ++k) {

	      Vir[k] = (interfaceLocation-X[i-1][0])/(X[i][0]-X[i-1][0])*V[i][k] - (interfaceLocation-X[i][0])/(X[i][0]-X[i-1][0])*V[i-1][k];
	      Vjr[k] = (interfaceLocation-X[j+2][0])/(X[j+1][0]-X[j+2][0])*V[j+1][k] - (interfaceLocation-X[j+1][0])/(X[j+1][0]-X[j+2][0])*V[j+2][k];
/*              betapr = 1.0;
              betapl = 1.0;
              if (limiterRight == 0) {
                betapr = 0.0;
                betapl = 0.0;
              }
 */
              if (k != 2 && k != 3 && limiterRight == 2) {
                if (Vjr[1] > 0.0) {
		  double rr = (V[j+3][k]-V[j+2][k]);
                  if (fabs(V[j+2][k]-V[j+1][k]) > 1.0e-8)
                    rr /= (V[j+2][k]-V[j+1][k]);
		  rr = std::max<double>(rr,0.0);
                  betapr = std::min<double>(betapr, rr);
		}
                if (Vir[1] < 0.0) {
		  double rr = (V[i-1][k]-V[i-2][k]);
                  if (fabs(V[i][k]-V[i-1][k]) > 1.0e-8)
                    rr /= (V[i][k]-V[i-1][k]);
                  else
                    rr = 1.0;
		  rr = std::max<double>(rr,0.0);
                  betapl = std::min<double>(betapl, rr);
		}
		betam[k] = std::min<double>(betapl,betapr);
              } else
                betam[k] = 1.0;
	      //std::cout << Vir[k] << " " << Vjr[k] << " " << V[i][k] << " " << V[j+1][k] << " " << V[i-1][k] << " " << V[j+2][k] << std::endl;
	    }
            //double beta = std::min<double>(betapl,betapr);
            //beta = sqrt(beta);
            /*std::cout << "beta = " << beta << " " << betapl << " " << betapr << std::endl;
            std::cout << "V[i-2]  = [" << V[i-2][0] << " " << V[i-2][1] << " " << V[i-2][4] << std::endl;
            std::cout << "V[i-1]  = [" << V[i-1][0] << " " << V[i-1][1] << " " << V[i-1][4] << std::endl;
            std::cout << "V[i]  = [" << V[i][0] << " " << V[i][1] << " " << V[i][4] << std::endl;
            std::cout << "V[j+1]  = [" << V[j+1][0] << " " << V[j+1][1] << " " << V[j+1][4] << std::endl;
            std::cout << "V[j+2]  = [" << V[j+2][0] << " " << V[j+2][1] << " " << V[j+2][4] << std::endl;
            std::cout << "V[j+3]  = [" << V[j+3][0] << " " << V[j+3][1] << " " << V[j+3][4] << std::endl;
*/
            double beta = betam[0];
/*            for (int k = 0; k < dim; ++k)
              beta = std::min<double>(beta,betam[k]);
            for (int k = 0; k < dim; ++k)
              betam[k] = beta;
	    for (int k = 0; k < dim; ++k) {

              std::cout << "beta[" << k << "] = " << betam[k] << std::endl;
            }
*/
            //betapl = sqrt(betapl);
            //betapr = sqrt(betapr);
            //betapl = betapr = beta;
	    for (int k = 0; k < dim; ++k) {
              Vir[k] = betapl*Vir[k]+(1.0-betapl)*V[i][k];
              Vjr[k] = betapr*Vjr[k]+(1.0-betapr)*V[j+1][k];
              //Vir[k] = betam[k]*Vir[k]+(1.0-betam[k])*V[i][k];
              //Vjr[k] = betam[k]*Vjr[k]+(1.0-betam[k])*V[j+1][k];
            }
	    //memcpy(Vir, V[i], sizeof(double)*dim);
	    //memcpy(Vjr, V[j+1], sizeof(double)*dim);
          } else {
	    memcpy(Vir, V[i], sizeof(double)*dim);
	    memcpy(Vjr, V[j+1], sizeof(double)*dim);
	    for (int k = 0; k < dim; ++k) {
	      //std::cout << i << " " << Vir[k] << " " << Vjr[k] << " " << V[i][k] << " " << V[j][k] << " " << V[i-1][k] << " " << V[j+1][k] << std::endl;
	    }
	    //std::cout << "p = " << varFcn->getPressure(Vjr,fluidId[j+1]) << std::endl;
	  }
        } else { // cutCellStatus[i] == 1
          
          if (interfaceExtrapolation == 1) {
	    for (int k = 0; k < dim; ++k) {
	      
	      Vir[k] = (interfaceLocation-X[i-2][0])/(X[i-1][0]-X[i-2][0])*V[i-1][k] - (interfaceLocation-X[i-1][0])/(X[i-1][0]-X[i-2][0])*V[i-2][k];
	      Vjr[k] = (interfaceLocation-X[j+1][0])/(X[j][0]-X[j+1][0])*V[j][k] - (interfaceLocation-X[j][0])/(X[j][0]-X[j+1][0])*V[j+1][k];
           /* 
              betapr = 1.0;
              betapl = 1.0;
              if (limiterRight == 0) {
                betapr = 0.0;
                betapl = 0.0;
              }
*/
              if (k != 2 && k != 3 && limiterRight == 2) {
                if (Vjr[1] > 0.0) {
		  double rr = (V[j+2][k]-V[j+1][k]);
                  if (fabs(V[j+1][k]-V[j][k]) > 1.0e-8)
                    rr /= (V[j+1][k]-V[j][k]);
                  else
                    rr = 1.0;
		  rr = std::max<double>(rr,0.0);
                  betapr = std::min<double>(betapr, rr);
		}
                if (Vir[1] < 0.0) {
		  double rr = (V[i-2][k]-V[i-3][k]);
                  if (fabs(V[i-1][k]-V[i-2][k]) > 1.0e-8)
                    rr /= (V[i-1][k]-V[i-2][k]);
                  else
                    rr = 1.0;
		  rr = std::max<double>(rr,0.0);
                  betapl = std::min<double>(betapl, rr);
		}
              }
              betam[k] = std::min<double>(betapl,betapr);
	      //std::cout << " " << Vir[k] << " " << Vjr[k] << " " << V[i][k] << " " << V[j][k] << " " << V[i-1][k] << " " << V[j+1][k] << std::endl;
	    }
            double beta = betam[0];
            for (int k = 0; k < dim; ++k)
              beta = std::min<double>(beta,betam[k]);
            for (int k = 0; k < dim; ++k)
              betam[k] = beta;
            //betapl = sqrt(betapl);
            //betapr = sqrt(betapr);
	    beta = sqrt(beta);
	    std::cout << betapl << " " << betapr << std::endl;
            //betapl = betapr = beta;
	    for (int k = 0; k < dim; ++k) {
              Vir[k] = betapl*Vir[k]+(1.0-betapl)*V[i-1][k];
              Vjr[k] = betapr*Vjr[k]+(1.0-betapr)*V[j][k];
//              Vir[k] = betam[k]*Vir[k]+(1.0-betam[k])*V[i-1][k];
//              Vjr[k] = betam[k]*Vjr[k]+(1.0-betam[k])*V[j][k];
            }
	    //memcpy(Vir, V[i-1], sizeof(double)*dim);
	    //memcpy(Vjr, V[j], sizeof(double)*dim);
          } else {
	    memcpy(Vir, V[i-1], sizeof(double)*dim);
	    memcpy(Vjr, V[j], sizeof(double)*dim);
	    for (int k = 0; k < dim; ++k) {
	      //std::cout << i << " " << Vir[k] << " " << Vjr[k] << " " << V[i][k] << " " << V[j][k] << " " << V[i-1][k] << " " << V[j+1][k] << std::endl;
	    }
          }
        }

        memset(fluxi,0,sizeof(double)*dim);
        memset(fluxj,0,sizeof(double)*dim);
	varFcn->getVarFcnBase(fluidId[i-1])->verification(0,Udummy,Vir);
	varFcn->getVarFcnBase(fluidId[j+1])->verification(0,Udummy,Vjr);

        riemann->computeRiemannSolution(Vir,Vjr,fluidId[i-1],fluidId[j+1],gradphi,varFcn,
	  			        Wir,Wjr,i,j,i,dx,0,false);

        if (lastPhaseChangeValue[i][0] < 0.0) {

          memcpy(lastPhaseChangeValue[i], Wir, sizeof(double)*5);

          //for (int mm = 0; mm < 5; ++mm)
  	  //  std::cout << lastPhaseChangeValue[i][mm] << " ";
        }
        if (lastPhaseChangeValue[j][0] < 0.0) {

          memcpy(lastPhaseChangeValue[j], Wjr, sizeof(double)*5);
          //for (int mm = 0; mm < 5; ++mm)
  	  //  std::cout << lastPhaseChangeValue[j][mm] << " ";
        }

        if (interfaceExtrapolation == 1) {
          if (cutCellStatus[j] == 1) {
	    for (int k = 0; k < dim; ++k) {
	      //Wi[k] = (Y[i+1][0]-X[i][0])/(interfaceLocation-X[i][0])*Wir[k] + (interfaceLocation-Y[i+1][0])/(interfaceLocation-X[i][0])*V[i][k];
	      Wi[k] = (Y[i+1][0]-X[i-1][0])/(interfaceLocation-X[i-1][0])*Wir[k] + (interfaceLocation-Y[i+1][0])/(interfaceLocation-X[i-1][0])*V[i-1][k];
              Wi[k] = betapl*Wi[k]+(1.0-betapl)*Wir[k];
              //Wi[k] = betam[k]*Wi[k]+(1.0-betam[k])*Wir[k];
	    }
          } else {
	    for (int k = 0; k < dim; ++k) {
	      //Wj[k] = (Y[i+1][0]-X[j][0])/(interfaceLocation-X[j][0])*Wjr[k] + (interfaceLocation-Y[i+1][0])/(interfaceLocation-X[j][0])*V[j][k];
	      Wj[k] = (Y[i+1][0]-X[j+1][0])/(interfaceLocation-X[j+1][0])*Wjr[k] + (interfaceLocation-Y[i+1][0])/(interfaceLocation-X[j+1][0])*V[j+1][k];
              Wj[k] = betapr*Wj[k]+(1.0-betapr)*Wjr[k];
              //Wj[k] = betam[k]*Wj[k]+(1.0-betam[k])*Wjr[k];
            }
          }
          //memcpy(Wi, Wir, sizeof(double)*dim);
	  //memcpy(Wj, Wjr, sizeof(double)*dim);
	} else {
          memcpy(Wi, Wir, sizeof(double)*dim);
	  memcpy(Wj, Wjr, sizeof(double)*dim);
        }

	if (interfaceTreatment == 0) {
	  memcpy(Wr[i], Wj, sizeof(double)*5);
	  memcpy(Wr[j], Wi, sizeof(double)*5);
	} else {
	  if (cutCellStatus[j] == 1)
	    memcpy(Wr[i], Wi, sizeof(double)*5);
	  else
	    memcpy(Wr[j], Wj, sizeof(double)*5);
	  
	}

	  
	riemannStatus[i] = riemannStatus[j] = 1;

        //fprintf(stderr,"Edge %d-%d crosses interface\n",i,j);

        if (cutCellStatus[j] == 1) {
          fluxFcn[0]->compute(length, 0.0, normal, normalVel, Wi, Wi, fluxi, fidi);
        } else {

          fluxFcn[0]->compute(length, 0.0, normal, normalVel, Wj, Wj, fluxj, fidj);
        }

        for (int k=0; k<dim; k++){
          R[i][k] += ctrlSurf[iEdge+1][0]*fluxi[k];
          R[j][k] -= ctrlSurf[iEdge+1][0]*fluxj[k];
        }
      }
    }
  }
	
  double dummy[dim];
  // flux at maxDistance
  fluxFcn[2]->compute(0.0, 0.0, normal, normalVel, V[numPoints-1], BC[1], flux, fluidId[numPoints-1]);
  for (int k=0; k<dim; ++k)
    R[numPoints-1][k] += ctrlSurf[numPoints][0]*flux[k];

  //R[numPoints-1][1] -= varFcn->getPressure(V[numPoints-1],0)*ctrlSurf[numPoints][0];

  // flux at left (for cartesian) - use of non-reflecting BC
  // flux at center (for cylindrical and spherical) - use of wall=symmetry
  normal[0] = ctrlSurf[0][0];
  //if(coordType == OneDimensionalInfo::CARTESIAN)
  //   fluxFcn[2]->compute(0.0, 0.0, normal, normalVel, V[0], BC[0], flux, fluidId[0]);
  /*else*/ if(coordType == OneDimensionalInfo::CYLINDRICAL)
    fluxFcn[1]->compute(0.0, 0.0, normal, normalVel, V[0], dummy, flux, fluidId[0]);
  else if(coordType == OneDimensionalInfo::SPHERICAL || coordType == OneDimensionalInfo::CARTESIAN)
    fluxFcn[1]->compute(0.0, 0.0, normal, normalVel, V[0], dummy, flux, fluidId[0]);

  //std::cout << "fnew = " << flux[0] << " " << flux[1] << " " << flux[2] << " " << flux[3] << " " << flux[4] << std::endl;

  for (int k=0; k<dim; ++k)
    R[0][k] -= ctrlSurf[0][0]*flux[k];
  
  
  if (source) {
    EulerSource E(varFcn,fluidId);
    source->compute(E,y, R, X, Y,fluidId);
  }

  if (interfaceTreatment == 1) {
    for (int i = 0; i < numPoints; ++i) {

      if (cutCellStatus[i])
        for (int k=0; k<dim; ++k)
          R[i][k] = 0.0;
    }
  }

  //std::cout << R*R << std::endl;

  // for debug
  //cout<<"flux[0] = "<<flux[0]<<" "<<flux[1]<<" "<<flux[2]<<" "<<flux[3]<<" "<<flux[4]<<endl;
  //cout<<"ctrlSurf[0] = "<<ctrlSurf[0][0]<<endl;
}
//------------------------------------------------------------------------------
void OneDimensional::PhiF(double t, SVec<double,1>& y,SVec<double,1>& k) {

  Rphi = 0.0;
  computeLevelSetFluxes(y);
  for(int i=0; i<numPoints; i++){
    k[i][0] = -Rphi[i][0] / ctrlVol[i][0];
  }
}

class LevelSetSource {

public:

  LevelSetSource(VarFcn* _vf,SVec<double,5>& _lu) :  vf(_vf), locU(_lu) { }

  ~LevelSetSource() { }

  void compute(double* U, double* f,int i) {

    f[0] = locU[i][1]*U[0]/locU[i][0];
  }

  VarFcn* vf;
  SVec<double,5>& locU;
};

double max(double e, double f) {
  return (e > f ? e : f);
}

double max(double d, double e, double f) {
  double q = max(e,f);
  return (d > q ? d : q);
}

double max(double c,double d, double e, double f) {
  double q = max(d,e,f);
  return (c > q ? c : q);
}

double max(double b,double c,double d, double e, double f) {
  double q = max(c,d,e,f);
  return (b > q ? b : q);
}

double max(double a,double b,double c,double d, double e, double f) {
  double q = max(b,c,d,e,f);
  return (a > q ? a : q);
}

void OneDimensional::computeLevelSetFluxes(SVec<double,1>& y){

  double uroe,flux;

  varFcn->conservativeToPrimitive(U,V,&fluidId);
  
  computeSlopes(y,Phislope,fluidId,false);

  if (levelSetMethod == 0) {
    double Phii[1],Phij[1],Phisi[1],Phisj[1];
    for(int iEdge=0; iEdge<numPoints-1; iEdge++){
      int i = iEdge;
      int j = iEdge+1;
      
      if (!isSixthOrder || !(i > 0)) {
	
	Phisi[0] = Phislope[i][0]*(X[j][0]-X[i][0]);
	Phisj[0] = Phislope[j][0]*(X[j][0]-X[i][0]);
      } else {
	
	Phisi[0] = (y[j][0]-y[i][0])/(X[j][0]-X[i][0])+(y[i][0]-y[i-1][0])/(X[i][0]-X[i-1][0])+
	  (-1.0/30.0)*((y[j+1][0]-y[j][0])/(X[j+1][0]-X[j][0])-2.0*(y[j][0]-y[i][0])/(X[j][0]-X[i][0])+(y[i][0]-y[i-1][0])/(X[i][0]-X[i-1][0]))/beta+
	  (-2.0/15.0)*(Phislope[i-1][0]-2.0*Phislope[i][0]+Phislope[j][0])/beta;
	
	Phisj[0] = (y[j][0]-y[i][0])/(X[j][0]-X[i][0])+(y[j+1][0]-y[j][0])/(X[j+1][0]-X[j][0])+
	  (-1.0/30.0)*((y[j+1][0]-y[j][0])/(X[j+1][0]-X[j][0])-2.0*(y[j][0]-y[i][0])/(X[j][0]-X[i][0])+(y[i][0]-y[i-1][0])/(X[i][0]-X[i-1][0]))/beta+
	  (-2.0/15.0)*(Phislope[j+1][0]-2.0*Phislope[j][0]+Phislope[i][0])/beta;
	
	Phisi[0] *= 0.5*(X[j][0]-X[i][0]);
	Phisj[0] *= 0.5*(X[j][0]-X[i][0]);
	
      }
      
      recFcnLS->compute(y[i], Phisi,y[j], Phisj, Phii, Phij);
      
      double uav = 0.5*(V[i][1]+V[j][1]);
      
      if(uav > 0.0){
	flux = Phii[0]*uav;
      }
      else{
	flux = Phij[0]*uav;
      }
      Rphi[i][0] += flux*ctrlSurf[iEdge+1][0];
      Rphi[j][0] -= flux*ctrlSurf[iEdge+1][0];
    }
    
    // flux at maxDistance
    flux = (V[numPoints-1][1] > 0) ? Phi[numPoints-1][0]*V[numPoints-1][1] : V[numPoints-1][0]*BCphi[1]*V[numPoints-1][1];
    Rphi[numPoints-1][0] += flux*ctrlSurf[numPoints][0];
    
    // flux at center is zero since u=0 for radial and spherical
    // but for cartesian, flux needs to be computed
    if(coordType == OneDimensionalInfo::CARTESIAN){
      flux = (V[0][1] < 0) ? Phi[0][0]*V[0][1] : V[0][0]*BCphi[0]*V[0][1];
      Rphi[0][0] -= flux*ctrlSurf[0][0];
    }
    
    // source term
    
    if (source) {
      for (int i = 0; i < numPoints; ++i)
        Rphi[i][0] += 2.0*U[i][1]/U[i][0]*Phi[i][0]*ctrlVol[i][0]/(X[i][0] > 0 ? X[i][0] : 1e8);
        //LevelSetSource E(varFcn,U);
      
        //source->compute(E,y, Rphi, X, Y,fluidId);
    }
  } else if (levelSetMethod == 1) // H-J WENO 
    {
      for(int i=0; i<numPoints; i++){
	double u = V[i][1];
	
	double v1 = 0.0,v2=0.0,v3=0.0,v4=0.0,v5=0.0;
	if (u > 0) {
	  if (i > 2)
	    v1 = (y[i-2][0]-y[i-3][0])/(X[i-2][0]-X[i-3][0]);
	  if (i > 1)
	    v2 = (y[i-1][0]-y[i-2][0])/(X[i-1][0]-X[i-2][0]);
	  if (i > 0)
	    v3 = (y[i][0]-y[i-1][0])/(X[i][0]-X[i-1][0]);
	  if (i < numPoints-1)
	    v4 = (y[i+1][0]-y[i][0])/(X[i+1][0]-X[i][0]);
	  if (i < numPoints-2)
	    v5 = (y[i+2][0]-y[i+1][0])/(X[i+2][0]-X[i+1][0]);
	} else {
	  if (i > 1)
	    v5 = (y[i-1][0]-y[i-2][0])/(X[i-1][0]-X[i-2][0]);
	  if (i > 0)
	    v4 = (y[i][0]-y[i-1][0])/(X[i][0]-X[i-1][0]);
	  if (i < numPoints-1)
	    v3 = (y[i+1][0]-y[i][0])/(X[i+1][0]-X[i][0]);
	  if (i < numPoints-2)
	    v2 = (y[i+2][0]-y[i+1][0])/(X[i+2][0]-X[i+1][0]);
	  if (i < numPoints-3)
	    v1 = (y[i+3][0]-y[i+2][0])/(X[i+3][0]-X[i+2][0]);
	}

	double S1 = 13.0/12.0*(v1-2.0*v2+v3)*(v1-2.0*v2+v3)+0.25*(v1-4.0*v2+3.0*v3)*(v1-4.0*v2+3.0*v3);
	double S2 = 13.0/12.0*(v2-2.0*v3+v4)*(v2-2.0*v3+v4)+0.25*(v2-v4)*(v2-v4);
	double S3 = 13.0/12.0*(v3-2.0*v4+v5)*(v3-2.0*v4+v5)+0.25*(3.0*v3-4.0*v4+v5)*(3.0*v3-4.0*v4+v5);

	double eps = 1e-6*max(v1*v1,v2*v2,v3*v3,v4*v4,v5*v5)+1.0e-30;
	double alpha1 = 0.1/((S1+eps)*(S1+eps));
	double alpha2 = 0.6/((S2+eps)*(S2+eps));
	double alpha3 = 0.3/((S3+eps)*(S3+eps));
	double asum = (alpha1+alpha2+alpha3);
	double om1 = alpha1/asum, om2 = alpha2/asum, om3 = alpha3/asum;
	
	double phi1 = v1/3.0-7.0/6.0*v2+11.0/6.0*v3;
	double phi2 = -v2/6.0+5.0/6.0*v3+1.0/3.0*v4;
	double phi3 = v3/3.0+5.0/6.0*v4-v5/6.0;
	double phix = om1*phi1+om2*phi2+om3*phi3;
	/*double phix = 0.0;
	if (u > 0) {
	  if (i > 1) {
	    double chi = (X[i][0]-X[i-1][0])/(X[i-1][0]-X[i-2][0]);
	    phix = (1.0+2.0*chi)/(1.0+chi)*y[i][0] + (-1.0-chi)*y[i-1][0]+(-(1.0+2.0*chi)/(1.0+chi)+1.0+chi)*y[i-2][0];
	    phix /= (X[i][0]-X[i-1][0]);
	  } else if (i > 0)
	    phix = (y[i][0]-y[i-1][0])/(X[i][0]-X[i-1][0]);
	} else if (u < 0) {
	  if (i < numPoints-2) {
	    double chi = (X[i][0]-X[i+1][0])/(X[i+1][0]-X[i+2][0]);
	    phix = (1.0+2.0*chi)/(1.0+chi)*y[i][0] + (-1.0-chi)*y[i+1][0]+(-(1.0+2.0*chi)/(1.0+chi)+1.0+chi)*y[i+2][0];
	    phix /= (X[i][0]-X[i+1][0]);
	  } else if (i < numPoints-1)
	    phix = (y[i+1][0]-y[i][0])/(X[i+1][0]-X[i][0]);
	    }*/
	Rphi[i][0] = u*phix * ctrlVol[i][0];
	 
      }   
      
      /*if (source) {
	LevelSetSource E(varFcn,U);
	
	source->compute(E,y, Rphi, X, Y,fluidId);
	}*/
    } else if (levelSetMethod == 3) {

      double Phii[1],Phij[1],Phisi[1],Phisj[1];
      for(int iEdge=0; iEdge<numPoints-1; iEdge++){
        int i = iEdge;
        int j = iEdge+1;
      
        if (!isSixthOrder || !(i > 0)) {
	
  	  Phisi[0] = Phislope[i][0]*(X[j][0]-X[i][0]);
	  Phisj[0] = Phislope[j][0]*(X[j][0]-X[i][0]);
        } else {
	
	  Phisi[0] = (y[j][0]-y[i][0])/(X[j][0]-X[i][0])+(y[i][0]-y[i-1][0])/(X[i][0]-X[i-1][0])+
	    (-1.0/30.0)*((y[j+1][0]-y[j][0])/(X[j+1][0]-X[j][0])-2.0*(y[j][0]-y[i][0])/(X[j][0]-X[i][0])+(y[i][0]-y[i-1][0])/(X[i][0]-X[i-1][0]))/beta+
 	    (-2.0/15.0)*(Phislope[i-1][0]-2.0*Phislope[i][0]+Phislope[j][0])/beta;
	
	  Phisj[0] = (y[j][0]-y[i][0])/(X[j][0]-X[i][0])+(y[j+1][0]-y[j][0])/(X[j+1][0]-X[j][0])+
	    (-1.0/30.0)*((y[j+1][0]-y[j][0])/(X[j+1][0]-X[j][0])-2.0*(y[j][0]-y[i][0])/(X[j][0]-X[i][0])+(y[i][0]-y[i-1][0])/(X[i][0]-X[i-1][0]))/beta+
	    (-2.0/15.0)*(Phislope[j+1][0]-2.0*Phislope[j][0]+Phislope[i][0])/beta;
	
	  Phisi[0] *= 0.5*(X[j][0]-X[i][0]);
	  Phisj[0] *= 0.5*(X[j][0]-X[i][0]);
	
        }  
      
        recFcnLS->compute(y[i], Phisi,y[j], Phisj, Phii, Phij);
      
        double uav = 0.5*(V[i][1]+V[j][1]);
      
        if(uav > 0.0){
	  flux = Phii[0]*uav;
        }
        else{
	  flux = Phij[0]*uav;
        }
        Rphi[i][0] += flux*ctrlSurf[iEdge+1][0];
        Rphi[j][0] -= flux*ctrlSurf[iEdge+1][0];
      }
    
      // flux at maxDistance
      flux = (V[numPoints-1][1] > 0) ? Phi[numPoints-1][0]*V[numPoints-1][1] : V[numPoints-1][0]*BCphi[1]*V[numPoints-1][1];
      Rphi[numPoints-1][0] += flux*ctrlSurf[numPoints][0];
    
      // flux at center is zero since u=0 for radial and spherical
      // but for cartesian, flux needs to be computed
      if(coordType == OneDimensionalInfo::CARTESIAN){
        flux = (V[0][1] < 0) ? Phi[0][0]*V[0][1] : V[0][0]*BCphi[0]*V[0][1];
        Rphi[0][0] -= flux*ctrlSurf[0][0];
      }
    
      // source term
    
      if (source) {
        for (int i = 0; i < numPoints; ++i)
          Rphi[i][0] += 2.0*U[i][1]/U[i][0]*Phi[i][0]*ctrlVol[i][0]/(X[i][0] > 0 ? X[i][0] : 1e8);
          //LevelSetSource E(varFcn,U);
      
          //source->compute(E,y, Rphi, X, Y,fluidId);
      }


      double integrand;      
      for (int i = 1; i < numPoints-1; ++i) {

        for (int q = 0; q < 10; ++q) {

          double s = ((double)q + 0.5) / 10.0;
          double phi;
          if (s < 0.5) {
            phi = Phi[i-1][0] * (0.5-s) + Phi[i][0]*(s+0.5);
            if (Phi[i][0]*Phi[i-1][0] > 0.0) {
              integrand = phi*(Vslope[i][1]);
            }  else {
              double isect = -Phi[i-1][0] /(Phi[i][0]-Phi[i-1][0]);
              if (isect - 0.5 > s)
                integrand = phi*Vslope[i-1][1];
              else
                integrand = phi*(Vslope[i][1]);
            }
          } else  {
            phi = Phi[i+1][0] * (s-0.5) + Phi[i][0]*(1.5-s);
            if (Phi[i][0]*Phi[i+1][0] > 0.0) {
              integrand = phi*(Vslope[i][1]);
            }  else {
              double isect = -Phi[i+1][0] /(Phi[i][0]-Phi[i+1][0]);
              if (1.5-isect < s)
                integrand = phi*Vslope[i+1][1];
              else
                integrand = phi*(Vslope[i][1]);
            }
          }
          Rphi[i][0] -= integrand*1.0/10.0*(Y[i+1][0]-Y[i][0]);
        }
      }
    }
   
}
//------------------------------------------------------------------------------
void OneDimensional::resultsOutput(double time, int iteration){

  fstream output;
  int i;

  varFcn->conservativeToPrimitive(U,V,&fluidId);

  //output2DVTK();

  for (i=0; i<PostFcn::SSIZE; ++i) {
    if (scalars[i]) {

      if (iteration == 0)
	output.open(scalars[i], fstream::out);
      else
	output.open(scalars[i], fstream::out | fstream::app);

      if (!output.good()) {
	std::cout << "*** Error: cannot open output file " << scalars[i] << std::endl;
	exit(-1);
      }

      output << time*refVal.time  << endl;
      for (int j = 0; j < numPoints; ++j) {
	switch ( (PostFcn::ScalarType)i ) {
	case PostFcn::DENSITY:
	  output << V[j][0]*sscale[i]; break;
	case PostFcn::MACH:
	  output << fabs(V[j][1]) / varFcn->computeSoundSpeed(V[j], fluidId[j])*sscale[i]; break;
	case PostFcn::PRESSURE:
	  output << varFcn->getPressure(V[j], fluidId[j])*sscale[i]; break;
	case PostFcn::TEMPERATURE:
	  output << varFcn->computeTemperature(V[j], fluidId[j])*sscale[i]; break;
	case PostFcn::PHILEVEL:
	  output << Phi[j][0] / V[j][0]*sscale[i]; break;
	case PostFcn::VELOCITY_NORM:
	  output << fabs(V[j][1])*sscale[i]; break;
	case PostFcn::FLUIDID:
	  output << fluidId[j]*sscale[i]; break;
	default:
	  break;	
	}  
	output << endl;
      }
      output << endl;
      output.close();
    }
  }
  for (i=0; i<PostFcn::VSIZE; ++i) {
    if (vectors[i]) {
      if (iteration == 0)
	output.open(vectors[i], fstream::out);	
      else
	output.open(vectors[i], fstream::out | fstream::app);

      if (!output.good()) {
	std::cout << "*** Error: cannot open output file " << vectors[i] << std::endl;
	exit(-1);
      }

      output << time*refVal.time << endl;
      for (int j = 0; j < numPoints; ++j) {
	switch ( (PostFcn::VectorType)i ) {
	case PostFcn::VELOCITY:
	  output << V[j][1]*vscale[i] << " " << 0.0 << " " << 0.0; break;
	default:
	  break;
	} 
        output << endl; 
      }
      output << endl;
      output.close();
    }
  }

  if (bubbleRadiusFile[0] != 0) {
    
    if (iteration == 0)
      output.open(bubbleRadiusFile, fstream::out);
    else
      output.open(bubbleRadiusFile, fstream::out | fstream::app);

    if (!output.good()) {
      std::cout << "*** Error: cannot open output file " << bubbleRadiusFile[0] << std::endl;
      exit(-1);
    }
    double rad = 0.0;
    for (i=0; i<numPoints; ++i) {
      if (fluidId[i+1] == 0) {
	if (levelSetMethod == 0)
	  rad = (X[i+1][0]*Phi[i][0]/V[i][0]-X[i][0]*Phi[i+1][0]/V[i+1][0])/(Phi[i][0]/V[i][0]-Phi[i+1][0]/V[i+1][0]);
	else if (levelSetMethod == 1 || levelSetMethod == 3)
	  rad = (X[i+1][0]*Phi[i][0]-X[i][0]*Phi[i+1][0])/(Phi[i][0]-Phi[i+1][0]);
	else
	  rad = interfaceLocation;
	break;
      }
    }
    output << time*refVal.time << " " <<  rad*refVal.length << endl;
    output.close();
  }
}

void OneDimensional::restartOutput(double time, int iteration){

  //if(iteration==0) cout << "outputting results in file "<<outfile<<endl;
  fstream output;
  int sp = strlen(outfile)+1;
  //char str[10];
  //sprintf(str,"%d",iteration);
  //char *currentfile = new char[sp + strlen(str) ];
  output.open(outfile, fstream::out | fstream::trunc);
  output.setf(ios::scientific);
  output.precision(20);

  if (!output.good()) {

    std::cout << "*** Error: Cannot open restart file " << outfile << std::endl;
    exit(-1);
  }

  cout << "outputting last results in file "<<outfile<<endl;

  varFcn->conservativeToPrimitive(U,V,&fluidId);

  double rad = 0.0;
  for (int i=0; i<numPoints; ++i) {
    if (fluidId[i+1] == 0) {
      if (levelSetMethod == 0)
	rad = (X[i+1][0]*Phi[i][0]/V[i][0]-X[i][0]*Phi[i+1][0]/V[i+1][0])/(Phi[i][0]/V[i][0]-Phi[i+1][0]/V[i+1][0]);
      else if (levelSetMethod == 1 || levelSetMethod == 3)
	rad = (X[i+1][0]*Phi[i][0]-X[i][0]*Phi[i+1][0])/(Phi[i][0]-Phi[i+1][0]);
      else
	rad = interfaceLocation;
      break;
    }
  }
  
  output << "# time = " << time*refVal.time << endl;
  output << "# " << numPoints << endl;
  for(int i=0; i<numPoints; i++) {
    output << X[i][0]*refVal.length <<" "<< V[i][0]*refVal.density <<" "<<
      V[i][1]*refVal.velocity <<" "<<
      varFcn->getPressure(V[i],fluidId[i])*refVal.pressure <<" ";
    if (levelSetMethod == 0)
      output << Phi[i][0]/V[i][0];
    else if (levelSetMethod == 1 || levelSetMethod == 3)
      output << Phi[i][0];
    else if (levelSetMethod == 2)
      output << -(X[i][0]-rad)*refVal.length;
    output << " " << fluidId[i] << " " << 
      varFcn->computeTemperature(V[i],fluidId[i])*refVal.temperature << endl;
  }
  output.close();

}

//------------------------------------------------------------------------------

template <int neq>
void OneDimensional::computeSlopes(SVec<double,neq>& VV, SVec<double,neq>& slopes,
				   Vec<int>& fid,bool crossInterface) {
  
  int i,j;
  double r,sig;
  int stat = 0;
  slopes = 0.0;
  for(i=0; i<numPoints; ++i) {

    //double A[4] = { X[i+1][0]*X[i+1][0]+X[i-1][0]*X[i-1][0]+X[i][0]*X[i][0], X[i+1][0]+X[i-1][0]+X[i][0],
    //		   X[i+1][0]+X[i-1][0]+X[i][0], 3};
  //double det = A[0]*A[3]-A[1]*A[2];
    double Xi = X[i][0], Xim1 = ((i > 0) ? X[i-1][0] : 0), Xip1 = ((i < numPoints-1) ? X[i+1][0] : 0);
    double* Vi = VV[i], *Vip1 = ((i < numPoints-1) ? VV[i+1] : 0), *Vim1 = ((i > 0) ? VV[i-1] : 0);
    double dx1 = Xip1-Xi, dx2 = Xim1-Xi;
    double dxsq = dx1*dx1+dx2*dx2;
    double* pSlopes = slopes[i];
    if (loctag[i])
      continue;
    
    //if (crossInterface) {
      if (i > 0 && i < numPoints-1 &&
	  ((interfaceTreatment == 0 && fid[i] == fid[i+1] && fid[i] == fid[i-1]) ||
	   (interfaceTreatment == 1 && cutCellStatus[i] == cutCellStatus[i-1] && cutCellStatus[i] == cutCellStatus[i+1])
	   || !crossInterface))
	stat = 0;
      else if (i < numPoints-1 && 
	       ((interfaceTreatment == 0 && fid[i] == fid[i+1]) ||
		(interfaceTreatment == 1 && cutCellStatus[i] == cutCellStatus[i+1]) ||
		!crossInterface))
	stat = 1;
      else if (i > 0 && 
	       ((interfaceTreatment == 0 && fid[i] == fid[i-1]) ||
		(interfaceTreatment == 1 && cutCellStatus[i] == cutCellStatus[i-1]) ||
		!crossInterface) )
	stat = 2;
      else
	stat = 3;
      //}
    //std::cout << stat << std::endl;
    for (j = 0; j < neq; ++j) {

      if (stat == 0) {   
	pSlopes[j] = ((Vip1[j]-Vi[j])*dx1+(Vim1[j]-Vi[j])*dx2)/(dxsq); 
      }
      else if (stat == 1) {
        pSlopes[j] = (Vip1[j]-Vi[j])/(Xip1-Xi);
      }
      else if (stat == 2)
	pSlopes[j] = (Vi[j]-Vim1[j])/(Xi-Xim1);
      else
	pSlopes[j] = 0.0;
    }
  }

}

RecFcn* OneDimensional::createRecFcn(IoData &ioData)
{

  RecFcn *rf = 0;

  double beta = ioData.schemes.ns.beta;
  double eps = ioData.schemes.ns.eps;

  if (ioData.eqs.type == EquationsData::NAVIER_STOKES &&
      ioData.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
  } else {
    if (ioData.schemes.ns.reconstruction == SchemeData::CONSTANT)
      rf = new RecFcnConstant<dim>;
    else if (ioData.schemes.ns.reconstruction == SchemeData::LINEAR) {
      if (ioData.schemes.ns.limiter == SchemeData::NONE)
	rf = new RecFcnLinear<dim>(beta, eps);
      else if (ioData.schemes.ns.limiter == SchemeData::VANALBADA)
	rf = new RecFcnVanAlbada<dim>(beta, eps);
      else if (ioData.schemes.ns.limiter == SchemeData::BARTH)
	rf = new RecFcnBarth<dim>(beta, eps);
      else if (ioData.schemes.ns.limiter == SchemeData::VENKAT)
	rf = new RecFcnVenkat<dim>(beta, eps);
      else if (ioData.schemes.ns.limiter == SchemeData::P_SENSOR)
	rf = new RecFcnLtdLinear<dim>(beta, eps);
    }
  }

  if (!rf) {
    fprintf(stderr, "*** Error: no valid choice for the reconstruction\n");
    exit(1);
  }

  return rf;

}

RecFcn* OneDimensional::createRecFcnLS(IoData &ioData)
{
  RecFcn *rf = 0;

  double beta = ioData.schemes.ls.beta;
  double eps = ioData.schemes.ls.eps;

  if (ioData.schemes.ls.reconstruction == SchemeData::CONSTANT)
    rf = new RecFcnConstant<1>;
  else if (ioData.schemes.ls.reconstruction == SchemeData::LINEAR) {
    if (ioData.schemes.ls.limiter == SchemeData::NONE)
      rf = new RecFcnLinear<1>(beta, eps);
    else if (ioData.schemes.ls.limiter == SchemeData::VANALBADA)
      rf = new RecFcnVanAlbada<1>(beta, eps);
    else if (ioData.schemes.ls.limiter == SchemeData::BARTH)
      rf = new RecFcnBarth<1>(beta, eps);
    else if (ioData.schemes.ls.limiter == SchemeData::VENKAT)
      rf = new RecFcnVenkat<1>(beta, eps);
    else if (ioData.schemes.ls.limiter == SchemeData::P_SENSOR)
      rf = new RecFcnLtdLinear<1>(beta, eps);
  }

  return rf;
}

void OneDimensional::loadSparseGrid(IoData& ioData) {

  if(ioData.mf.riemannComputation == MultiFluidData::TABULATION2){
    // only the ioData.eqs.fluidModel is considered since only the Riemann invariant of one EOS is tabulated!
    // (no need to specify two different EOS)

    double *refIn  = new double[2];
    double *refOut = new double[1];
    refIn[0] = ioData.ref.rv.density;
    refIn[1] = pow(ioData.ref.rv.density,-ioData.eqs.fluidModel.jwlModel.omega)*ioData.ref.rv.velocity*ioData.ref.rv.velocity;
    refOut[0] = ioData.ref.rv.velocity;

    tabulationC = new SparseGridCluster;
    tabulationC->readFromFile(ioData.mf.sparseGrid.numberOfTabulations, refIn, refOut, ioData.mf.sparseGrid.tabulationFileName, 0);

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
    tabulationC->readFromFile(ioData.mf.sparseGrid.numberOfTabulations, refIn, refOut, ioData.mf.sparseGrid.tabulationFileName, 0);
  }else{
    tabulationC = 0;
  }
}

void OneDimensional::setupProbes(IoData& iod) {

  int i;
  for (i=0; i<PostFcn::SSIZE; ++i) {
    nodal_scalars[i] = 0; 
  }
  for (i=0; i<PostFcn::VSIZE; ++i) {
    nodal_vectors[i] = 0; 
  }
 
  int spn = strlen(iod.output.transient.probes.prefix) + 1;
  Probes& myProbes = iod.output.transient.probes;
  nodal_output.ids.resize(Probes::MAXNODES);
  nodal_output.alpha.resize(Probes::MAXNODES);
  nodal_output.locations.resize(Probes::MAXNODES);
  for (i = 0; i < Probes::MAXNODES; ++i) {
    nodal_output.locations[i] = Vec3D(myProbes.myNodes[i].locationX,
                                      myProbes.myNodes[i].locationY,
                                      myProbes.myNodes[i].locationZ);
    nodal_output.ids[i] = myProbes.myNodes[i].id;
    if (myProbes.myNodes[i].id >= 0) {
  
    } else if (myProbes.myNodes[i].locationX < -1.0e10)
      break;
    else {

      for (int j = 0; j < numPoints; ++j) {
        if (X[j][0] > myProbes.myNodes[i].locationX) {
          nodal_output.ids[i] = j-1;
          nodal_output.alpha[i] = -(myProbes.myNodes[i].locationX-X[j][0])/(X[j][0]-X[j-1][0]);
          break;
        }
      }
    }
  }

  nodal_output.numNodes = i;

  if (iod.output.transient.probes.density[0] != 0) {
    nodal_scalars[PostFcn::DENSITY] = new char[spn + strlen(iod.output.transient.probes.density)];
    sprintf(nodal_scalars[PostFcn::DENSITY], "%s%s",
            iod.output.transient.probes.prefix, iod.output.transient.probes.density);
  }
  if (iod.output.transient.probes.pressure[0] != 0) {
    nodal_scalars[PostFcn::PRESSURE] = new char[spn + strlen(iod.output.transient.probes.pressure)];
    sprintf(nodal_scalars[PostFcn::PRESSURE], "%s%s",
            iod.output.transient.probes.prefix, iod.output.transient.probes.pressure);
  }       
  if (iod.output.transient.probes.temperature[0] != 0) {
    nodal_scalars[PostFcn::TEMPERATURE] = new char[spn + strlen(iod.output.transient.probes.temperature)];
    sprintf(nodal_scalars[PostFcn::TEMPERATURE], "%s%s",
            iod.output.transient.probes.prefix, iod.output.transient.probes.temperature);
  } 
  if (iod.output.transient.probes.velocity[0] != 0) {
    nodal_vectors[PostFcn::VELOCITY] = new char[spn + strlen(iod.output.transient.probes.velocity)];
    sprintf(nodal_vectors[PostFcn::VELOCITY], "%s%s",
            iod.output.transient.probes.prefix, iod.output.transient.probes.velocity);
  }

}

void OneDimensional::outputProbes(double time,int iteration) {

  if (nodal_output.numNodes == 0)
    return;
  
  int i,j;
  const char* mode = iteration ? "a" : "w";

  for (i=0; i<PostFcn::SSIZE; ++i) {
    if (nodal_scalars[i]) {
      FILE* out = fopen(nodal_scalars[i],mode);

      if (!out) {

	std::cout << "*** Error: cannot open probes file " << nodal_scalars[i] << std::endl;
	exit(-1);
      }

      fprintf(out,"%d %e ",iteration,time*refVal.time);
      for (j = 0; j < nodal_output.numNodes; ++j) {
        if (nodal_output.locations[i].v[0] >= 0.0) {
          switch ((PostFcn::ScalarType)i) {

            case PostFcn::DENSITY:
              fprintf(out,"%e ", (V[nodal_output.ids[j]][0]*nodal_output.alpha[j] + V[nodal_output.ids[j]+1][0]*(1.0-nodal_output.alpha[j]))*refVal.density);
              break;
            // case PostFcn::VELOCITY:
            //  fprintf(out,"%lf ", V[nodal_output.ids[j]][1]*nodal_output.alpha[j] + V[nodal_output.ids[j]+1][1]*(1.0-nodal_output.alpha[j]));
            //  break;
            case PostFcn::TEMPERATURE:
              fprintf(out,"%e ", (varFcn->computeTemperature(V[nodal_output.ids[j]],fluidId[nodal_output.ids[j]])*nodal_output.alpha[j] + 
                                  varFcn->computeTemperature(V[nodal_output.ids[j]+1],fluidId[nodal_output.ids[j]+1])*(1.0-nodal_output.alpha[j]))*refVal.temperature);
              break;
            case PostFcn::PRESSURE:
              fprintf(out,"%e ", (varFcn->getPressure(V[nodal_output.ids[j]],fluidId[nodal_output.ids[j]])*nodal_output.alpha[j] +
                                  varFcn->getPressure(V[nodal_output.ids[j]+1],fluidId[nodal_output.ids[j]+1])*(1.0-nodal_output.alpha[j]))*refVal.pressure);
              break;
          }
        } else {
          switch ((PostFcn::ScalarType)i) {
  
            case PostFcn::DENSITY:
              fprintf(out,"%e ", V[nodal_output.ids[j]][0]*refVal.density);
              break;
            // case PostFcn::VELOCITY:
            //   fprintf(out,"%lf ", V[nodal_output.ids[j]][1]*nodal_output.alpha[j] + V[nodal_output.ids[j]+1][1]*(1.0-nodal_output.alpha[j]));
            //   break;
            case PostFcn::TEMPERATURE:
              fprintf(out,"%e ", varFcn->computeTemperature(V[nodal_output.ids[j]],fluidId[nodal_output.ids[j]])*refVal.temperature);
              break;
            case PostFcn::PRESSURE:
              fprintf(out,"%e ", varFcn->getPressure(V[nodal_output.ids[j]],fluidId[nodal_output.ids[j]])*refVal.pressure);
              break;
          }
        }
      }
      fprintf(out,"\n");
      fclose(out);
    }
  }

  for (i=0; i<PostFcn::VSIZE; ++i) {
    if (nodal_vectors[i]) {
      FILE* out = fopen(nodal_vectors[i],mode); 
      if (!out) {

	std::cout << "*** Error: cannot open probes file " << nodal_vectors[i] << std::endl;
	exit(-1);
      }
      fprintf(out,"%e ",time*refVal.time);
      for (j = 0; j < nodal_output.numNodes; ++j) {
        if (nodal_output.locations[i].v[0] >= 0.0) {
          switch ((PostFcn::VectorType)i) {
  
            case PostFcn::VELOCITY:
              fprintf(out,"%e ", (V[nodal_output.ids[j]][1]*nodal_output.alpha[j] + V[nodal_output.ids[j]+1][1]*(1.0-nodal_output.alpha[j]))*refVal.velocity);
              break;
          }
        } else {
          switch ((PostFcn::VectorType)i) {

            case PostFcn::VELOCITY:
              fprintf(out,"%e ", V[nodal_output.ids[j]][1]*refVal.velocity);
              break;
          }
        }
      }
      fprintf(out,"\n");
      fclose(out);
    }
  }

}

void OneDimensional::output2DVTK() {

  static int cnt = 0;

  int N = numPoints*2-1;

  N/=2;
 
  double x1 = -X[numPoints-2][0];
  double x2 = X[numPoints-2][0];
  
  double y1 = -X[numPoints-2][0];
  double y2 = X[numPoints-2][0];

  char fname[256];
  sprintf(fname,"pressure2d%d.vti",cnt++);
  std::ofstream out(fname);

  double x0 = 0.5*(x1+x2);
  double y0 = 0.5*(y1+y2);

  double dx = (x2-x1) / (N-1);
  double dy = (y2-y1) / (N-1);

  // Write out in VTK rectilinear grid format 
  out << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;

  out << "<ImageData WholeExtent=\"" << 0 << " " << N-1 << " " << 0 << " " << N-1 << " 0 0\"" << std::endl;

  out << "\tOrigin=\"" << x0 << " " << y0 << " 0.0 \" Spacing = \"" <<
	dx << " " << dy << " 0.0\">" << std::endl;

  out << "<Piece Extent=\"" << 0 << " " << N-1 << " " << 0 << " " << N-1 << " 0 0\">" << std::endl;

  out << "<PointData Scalars = \"" << "pressure" << "\">" << std::endl;
  
  out << "<DataArray type = \"Float64\" Name = \"" << "pressure" << "\" NumberOfComponents=\"1\"" << std::endl;
  out << "format=\"ascii\">" << std::endl;

  for (int j = 0; j < N; ++j) {
    for (int i = 0; i < N; ++i) {
      double r = sqrt(pow(x1+(double)i*dx,2.0)+pow(y1+(double)j*dx,2.0)); 
      int id = numPoints-2;
      double alpha = 0.0;
      for (int k = 0; k < numPoints-2; ++k) {
        if (X[k][0] > r) {
          id = k-1;
          alpha = -(r-X[k][0])/(X[k][0]-X[k-1][0]);
          break;
        }
      }
      double p = (varFcn->getPressure(V[id],fluidId[id])*alpha +
                  varFcn->getPressure(V[id+1],fluidId[id+1])*(1.0-alpha))*refVal.pressure;
      out << p << " ";
    }
  }
  
  out << "\n</DataArray>\n</PointData>\n</Piece>\n";

  out << "</ImageData>" << std::endl;

  out << "</VTKFile>" << std::endl;

  out.close();
  
}



template <int dimp,int dimLS2>
void OneDimensional::read1DSolution(IoData& iod, DistSVec<double,dimp>& Up, 
			     DistSVec<double,dimLS2>* Phi,
			     FluidSelector* fluidSelector,
			     VarFcn* varFcn,
			     DistSVec<double,3>& X,
			     Domain& dom,
			     ReadMode mode,
			     bool spherical) 
{
    // read 1D solution
    DistSVec<double,dimp> ut(Up);
    for (map<int, OneDimensionalInputData *>::iterator itr = iod.input.oneDimensionalInput.dataMap.begin();
            itr != iod.input.oneDimensionalInput.dataMap.end(); ++itr) {
      
      std::fstream input;
      int sp = strlen(iod.input.prefix) + 1;
      char *filename = new char[sp + strlen(itr->second->file)+5];
      sprintf(filename, "%s%s", iod.input.prefix, itr->second->file);
      input.open(filename, fstream::in);
      //cout << filename << endl;
      if (!input.is_open()) {
        cout<<"*** Error: could not open 1D solution file "<<filename<<endl;
        exit(1);
      }
      
      input.ignore(256,'\n');
      input.ignore(2,' ');
      int numPoints = 0;
      input >> numPoints;
      //cout << " Num 1d points = " << numPoints << endl;
      double* x_1D = new double[numPoints];
      double* v_1D = new double[numPoints*5];// rho, u, p, phi, T
      int* fids = new int[numPoints];
      
      double rad = 0;
      for(int i=0; i<numPoints; i++) {
        input >> x_1D[i] >> v_1D[i*5] >> v_1D[i*5+1] >> v_1D[i*5+2] >> v_1D[i*5+3]>> fids[i] >> v_1D[i*5+4];

        fids[i] = std::min<int>(fids[i], varFcn->size()-1);
        x_1D[i]    /= iod.ref.rv.length;
        v_1D[i*5] /= iod.ref.rv.density;
        v_1D[i*5+1] /= iod.ref.rv.velocity;
        v_1D[i*5+2] /= iod.ref.rv.pressure;
        //v_1D[i][3] /= iod.ref.rv.length;
        v_1D[i*5+4] /= iod.ref.rv.temperature;
        //std::cout << v_1D[i*5+2] << std::endl;
        if (rad == 0 && fids[i] == 0) {
          rad = (x_1D[i-1]*v_1D[i*5+3]-x_1D[i]*v_1D[(i-1)*5+3])/(v_1D[i*5+3]-v_1D[(i-1)*5+3]);
        }
      }
      
      input.close();

      int lsdim = 0;
      
      // interpolation assuming 1D solution is centered on bubble_coord0
      double bubble_x0 = itr->second->x0;
      double bubble_y0 = itr->second->y0;
      double bubble_z0 = itr->second->z0;
      double max_distance = x_1D[numPoints-1]; 
      double localRadius;
      ut = 1.0;
#pragma omp parallel for
      for (int iSub=0; iSub<Up.numLocSub(); ++iSub) {
        SVec<double,dimp> &u(Up(iSub));
        SVec<double, 3> &x(X(iSub));
        for(int i=0; i<u.size(); i++) {
          localRadius = sqrt((x[i][0]-bubble_x0)*(x[i][0]-bubble_x0)+(x[i][1]-bubble_y0)*(x[i][1]-bubble_y0)+(x[i][2]-bubble_z0)*(x[i][2]-bubble_z0));
          for (int k = 0; k < dim; ++k)
            ut.subData(iSub)[i][k] = 0.0;
	      }

        Veval veval(varFcn,itr->second,fluidSelector,x_1D,v_1D,fids,&x,numPoints,bubble_x0,bubble_y0,bubble_z0,spherical);

        double localRadius; int np;
        double localAlpha, velocity_r;
        double localV[5];
        if (mode == ModeU) {
          for(int i=0; i<u.size(); i++) {
            double v[5];
            localRadius = sqrt((x[i][0]-bubble_x0)*(x[i][0]-bubble_x0)+(x[i][1]-bubble_y0)*(x[i][1]-bubble_y0)+(x[i][2]-bubble_z0)*(x[i][2]-bubble_z0));
	    //if (localRadius < 1.0) {
	      //varFcn->conservativeToPrimitive(ut(iSub)[i],v,1);
	    //  std::cout << "2: " << localRadius << " " << ut(iSub)[i][4] << std::endl;
	    //} 
          }
          dom.getSubDomain()[iSub]->integrateFunction(&veval, x, ut(iSub), &Veval::Eval, 4);
          for(int i=0; i<u.size(); i++) {
            double v[5];
            localRadius = sqrt((x[i][0]-bubble_x0)*(x[i][0]-bubble_x0)+(x[i][1]-bubble_y0)*(x[i][1]-bubble_y0)+(x[i][2]-bubble_z0)*(x[i][2]-bubble_z0));
	    /*if (localRadius < 1.0) {
	      varFcn->conservativeToPrimitive(ut(iSub)[i],v,1);
	      //std::cout << "2: " << localRadius << " " << ut(iSub)[i][4] << std::endl;
	      }*/
          }
          dom.getSubDomain()[iSub]->sndData(*dom.getVecPat(), ut.subData(iSub) );
        } else {
          int lsdim=0;
          for(int i=0; i<u.size(); i++) {

            localRadius = sqrt((x[i][0]-bubble_x0)*(x[i][0]-bubble_x0)+(x[i][1]-bubble_y0)*(x[i][1]-bubble_y0)+(x[i][2]-bubble_z0)*(x[i][2]-bubble_z0));
	 
            int fid_new = fids[0];
            if (itr->second->fluidRemap.dataMap.find(fids[0]) != itr->second->fluidRemap.dataMap.end())
              fid_new = itr->second->fluidRemap.dataMap.find(fids[0])->second->newID;
              lsdim = fluidSelector->getLevelSetDim(0,fid_new);
              if (fluidSelector) {
                (*Phi)(iSub)[i][lsdim] = rad-localRadius;
              }
          }
        }
      }
      
      if (mode == ModeU) {
        dom.getVecPat()->exchange();
#pragma omp parallel for
        for (int iSub=0; iSub<Up.numLocSub(); ++iSub) {
          dom.getSubDomain()[iSub]->addRcvData(*dom.getVecPat(), ut.subData(iSub));
        }
    
        DistVec<double> A(Up.info());
        dom.getCommunicator()->barrier();
        dom.computeControlVolumes(1.0, X, A);
	
#pragma omp parallel for
        for (int iSub=0; iSub<Up.numLocSub(); ++iSub) {
          SVec<double,dimp> &u(Up(iSub));
          SVec<double,dimp> &utl(ut(iSub));
          SVec<double, 3> &x(X(iSub));
          for(int i=0; i<u.size(); i++) {
            assert(A(iSub)[i] > 0 && utl[i][0] > 0);
	    
            localRadius = sqrt((x[i][0]-bubble_x0)*(x[i][0]-bubble_x0)+(x[i][1]-bubble_y0)*(x[i][1]-bubble_y0)+(x[i][2]-bubble_z0)*(x[i][2]-bubble_z0));
	    
            if (localRadius < max_distance || !spherical) {
              double v[5];
              for (int k = 0; k <dimp; ++k)
                v[k] = utl[i][k] / A(iSub)[i];

              if (spherical)
                varFcn->primitiveToConservative(v,u[i],localRadius < rad ? 1 : 0);
		
	      
	      /*if (localRadius*iod.ref.rv.length < 0.5) {
		double v[5];
		varFcn->conservativeToPrimitive(u[i],v,1);
		std::cout << localRadius << " " <<v[4]*iod.ref.rv.pressure << std::endl;
		}*/ 
            }
          }
        }
      }
    }
}



template 
void OneDimensional::read1DSolution<5,1>(IoData& , DistSVec<double,5>& , DistSVec<double,1>* , FluidSelector* , VarFcn* , DistSVec<double,3>& , Domain& , ReadMode , bool );
template 
void OneDimensional::read1DSolution<5,2>(IoData& , DistSVec<double,5>& , DistSVec<double,2>* , FluidSelector* , VarFcn* , DistSVec<double,3>& , Domain& , ReadMode , bool );
template 
void OneDimensional::read1DSolution<5,3>(IoData& , DistSVec<double,5>& , DistSVec<double,3>* , FluidSelector* , VarFcn* , DistSVec<double,3>& , Domain& , ReadMode , bool );
template 
void OneDimensional::read1DSolution<6,1>(IoData& , DistSVec<double,6>& , DistSVec<double,1>* , FluidSelector* , VarFcn* , DistSVec<double,3>& , Domain& , ReadMode , bool );
template 
void OneDimensional::read1DSolution<6,2>(IoData& , DistSVec<double,6>& , DistSVec<double,2>* , FluidSelector* , VarFcn* , DistSVec<double,3>& , Domain& , ReadMode , bool );
template 
void OneDimensional::read1DSolution<6,3>(IoData& , DistSVec<double,6>& , DistSVec<double,3>* , FluidSelector* , VarFcn* , DistSVec<double,3>& , Domain& , ReadMode , bool );
template 
void OneDimensional::read1DSolution<7,1>(IoData& , DistSVec<double,7>& , DistSVec<double,1>* , FluidSelector* , VarFcn* , DistSVec<double,3>& , Domain& , ReadMode , bool );
template 
void OneDimensional::read1DSolution<7,2>(IoData& , DistSVec<double,7>& , DistSVec<double,2>* , FluidSelector* , VarFcn* , DistSVec<double,3>& , Domain& , ReadMode , bool );
template 
void OneDimensional::read1DSolution<7,3>(IoData& , DistSVec<double,7>& , DistSVec<double,3>* , FluidSelector* , VarFcn* , DistSVec<double,3>& , Domain& , ReadMode , bool );
