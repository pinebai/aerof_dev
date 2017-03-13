#include "FluidSelector.h"
#include "Domain.h"
#include "TriangulatedInterface.h"
#include <iostream>
#include <fstream>
using namespace std;

#include <TsInput.h>

#include "TsRestart.h"

//------------------------------------------------------------------------------

FluidSelector::FluidSelector(const int nPhases, IoData &ioData, Domain *dom) : iodp(&ioData)
{ 
  numPhases = nPhases;
  domain = dom ? dom : 0;

  char* fidpath;
  if (ioData.input.restart_file_package[0] == 0) {
    const std::string prefix(ioData.input.prefix); 
    fidpath = TsInput::absolutePath(ioData.input.fluidId, prefix);
  } else {

    char dummy[6][256];
    fidpath = new char[256];
    const std::string prefix(ioData.input.prefix); 
    char* tmp = TsInput::absolutePath(ioData.input.restart_file_package, prefix);

    TsRestart::readRestartFileNames(tmp,
				    dummy[0],
				    dummy[1],
				    dummy[2],
				    dummy[3],
				    fidpath,
				    dummy[4],
				    dummy[5], dom->getCommunicator());

//   std::cout << std::string(tmp) << std::endl;
    //std::cout << std::string(fidpath) << std::endl;
  }
  
  fluidId  = 0;
  fluidIdn = 0;
  if(domain){
    fluidId  = new DistVec<int>(domain->getNodeDistInfo());
    fluidIdn = new DistVec<int>(domain->getNodeDistInfo());
    *fluidId  = 0;
    *fluidIdn = 0;
  }

  fluidIdnm1 = 0;
  fluidIdnm2 = 0;
  if(ioData.ts.implicit.type == ImplicitData::THREE_POINT_BDF){
    fluidIdnm1 = new DistVec<int>(domain->getNodeDistInfo());
    *fluidIdnm1 = 0;
  }
  else if(ioData.ts.implicit.type == ImplicitData::FOUR_POINT_BDF){
    fluidIdnm1 = new DistVec<int>(domain->getNodeDistInfo());
    *fluidIdnm1 = 0;
    fluidIdnm2 = new DistVec<int>(domain->getNodeDistInfo());
    *fluidIdnm2 = 0;
  }

  if (ioData.input.fluidId[0] != 0 ||
      ioData.input.restart_file_package[0] != 0) {
    
    //std::cout << "Reading solution " << std::string(fidpath) << std::endl;
    double scl;
    DistSVec<int,1> i1(fluidId->info(), reinterpret_cast<int(*)[1]>(fluidId->data()));
    domain->readVectorFromFile(fidpath, 0, &scl, i1);

    if (fluidIdnm1){
      DistSVec<int,1> i2(fluidIdnm1->info(), reinterpret_cast<int(*)[1]>(fluidIdnm1->data()));
      domain->readVectorFromFile(fidpath, 1,  &scl, i2);
    }

    if (fluidIdnm2){
      DistSVec<int,1> i3(fluidIdnm2->info(), reinterpret_cast<int(*)[1]>(fluidIdnm2->data()));
      domain->readVectorFromFile(fidpath, 2,  &scl,i3);
    }

  }


  programmedBurn = 0;

  ownsData = true;
}

FluidSelector::FluidSelector(DistVec<int>& nodeTag, Domain* dom) : domain(dom) {

  fluidId  = 0;
  fluidIdn = 0;

  fluidIdnm1 = 0;
  fluidIdnm2 = 0;

  ownsData = false;

  fluidId = &nodeTag;

}

//------------------------------------------------------------------------------
#define SAFE_DELETE(f) if (f) delete f;
FluidSelector::~FluidSelector()
{
  if (ownsData) {
    SAFE_DELETE(fluidId);
    SAFE_DELETE(fluidIdn);
    SAFE_DELETE(fluidIdnm1);
    SAFE_DELETE(fluidIdnm2);
  }
}

//------------------------------------------------------------------------------
// takes in the fluidId computed by intersector. Complete it with user-specified FF interfaces.
/*void FluidSelector::initializeFluidIds(DistVec<int> &fsId, DistSVec<double,3> &X, IoData &ioData)
{
  *fluidId = fsId;
  setupFluidIdVolumesInitialConditions(ioData);
  if(ioData.input.oneDimensionalSolution[0] != 0)
    setupFluidIdOneDimensionalSolution(ioData,X);
  setupFluidIdMultiFluidInitialConditions(ioData,X);

  *fluidIdn = *fluidId;
  if(fluidIdnm1) *fluidIdnm1 = *fluidId;
  if(fluidIdnm2) *fluidIdnm2 = *fluidId;
}
*/
//------------------------------------------------------------------------------

void FluidSelector::getFluidId(int &tag, double *phi){
  tag = 0;
  for(int i=0; i<numPhases-1; i++)
    if(phi[i]>0.0) {tag = i+1; return; }
}

//------------------------------------------------------------------------------

void FluidSelector::getFluidId(Vec<int> &tag, Vec<double> &phi){
  for(int iNode=0; iNode<phi.size(); iNode++)
    tag[iNode] = (phi[iNode]<0.0) ? 0 : 1;
}

//------------------------------------------------------------------------------

void FluidSelector::getFluidId(DistVec<double> &Phi){
  int numLocSub = Phi.numLocSub();
  int iSub;
#pragma omp parallel for
  for(iSub=0; iSub<numLocSub; ++iSub) {
    double *phi = Phi.subData(iSub);
    int    *tag = fluidId->subData(iSub);
    for(int iNode=0; iNode<Phi.subSize(iSub); iNode++)
      tag[iNode] = (phi[iNode]<0.0) ? 0 : 1; 
  }
}

//------------------------------------------------------------------------------

void FluidSelector::getFluidId(DistVec<int> &Tag, DistVec<double> &Phi){
  int numLocSub = Phi.numLocSub();
  int iSub;
#pragma omp parallel for
  for(iSub=0; iSub<numLocSub; ++iSub) {
    double *phi = Phi.subData(iSub);
    int    *tag = Tag.subData(iSub);
    for(int iNode=0; iNode<Phi.subSize(iSub); iNode++)
      tag[iNode] = (phi[iNode]<0.0) ? 0 : 1;
  }
}

//------------------------------------------------------------------------------

int FluidSelector::getLevelSetDim(int fluidId1, int fluidId2, int node1, int node2){
  if(fluidId1 == fluidId2){
    fprintf(stdout, "*** Error: getLevelSetDim should not be called when there is no interface.\n");
    exit(1);
  }
  if(fluidId1 < 0 || fluidId2 < 0){
    fprintf(stdout, "*** Error: arguments of getLevelSetDim (%d %d)should be positive\n", fluidId1, fluidId2);
    exit(1);
  }
  int burnTag;
  if (programmedBurn && programmedBurn->isBurnedEOS(fluidId1,burnTag)) {
    return programmedBurn->getUnburnedEOS(burnTag)-1;
  }
  if (programmedBurn && programmedBurn->isBurnedEOS(fluidId2,burnTag)) {
    return programmedBurn->getUnburnedEOS(burnTag)-1;
  }
    
  if(fluidId1 * fluidId2 != 0){
    fprintf(stdout, "*** Error: it is assumed that all interfaces are between any fluid and fluid 0\n");
    fprintf(stdout, "***      : here fluidIds are %d and %d for nodes %d and %d\n", fluidId1, fluidId2, node1, node2);
#ifdef AEROF_MPI_DEBUG
    DebugTools::SpitRank();
#endif 
    exit(1);
  }
  int fId = fluidId1 + fluidId2;
  return fId-1;
}

//------------------------------------------------------------------------------
/*
void FluidSelector::setupFluidIdVolumesInitialConditions(IoData &iod)
{
  // loop on all Volumes to setup U0
  if(!iod.volumes.volumeMap.dataMap.empty()){
    if(!domain)
      fprintf(stderr,"ERROR: Unable to setup fluidId using user-specified volume ID(s).\n");
    map<int, VolumeData *>::iterator volIt;
    for (volIt=iod.volumes.volumeMap.dataMap.begin(); volIt!=iod.volumes.volumeMap.dataMap.end();volIt++)
      if(volIt->second->type==VolumeData::FLUID){
        //each volume (volIt->first) is setup using Input variables 'volumeInitialConditions'
        //                                 and equation of state 'fluidModel'
        map<int, FluidModelData *>::iterator fluidIt = iod.eqs.fluidModelMap.dataMap.find(volIt->second->fluidModelID);
        if(fluidIt == iod.eqs.fluidModelMap.dataMap.end()){
          fprintf(stderr, "*** Error: fluidModelData[%d] could not be found\n", volIt->second->fluidModelID);
          exit(1);
        }
        // initialize fluidId.
        domain->setupFluidIdVolumesInitialConditions(volIt->first, volIt->second->fluidModelID, *fluidId);
      }
  }
}

//------------------------------------------------------------------------------

void FluidSelector::setupFluidIdOneDimensionalSolution(IoData &iod, DistSVec<double,3> &X)
{
  // read 1D solution
  fstream input;
  int sp = strlen(iod.input.prefix) + 1;
  char *filename = new char[sp + strlen(iod.input.oneDimensionalSolution)];
  sprintf(filename, "%s%s", iod.input.prefix, iod.input.oneDimensionalSolution);
  input.open(filename, fstream::in);
  if (!input.is_open()) {
    cout<<"*** Error: could not open 1D solution file "<<filename<<endl;
    exit(1);
  }

  input.ignore(256,'\n');
  input.ignore(2,' ');
  int numPoints = 0;
  input >> numPoints;
  cout <<"number of points in 1D solution is " << numPoints <<endl;
  double x_1D[numPoints];
  double v_1D[numPoints][4];// rho, u, p, phi

  for(int i=0; i<numPoints; i++) {
    input >> x_1D[i] >> v_1D[i][0] >> v_1D[i][1] >> v_1D[i][2] >> v_1D[i][3];
    x_1D[i]    /= iod.ref.rv.length;
    v_1D[i][0] /= iod.ref.rv.density;
    v_1D[i][1] /= iod.ref.rv.velocity;
    v_1D[i][2] /= iod.ref.rv.pressure;
    v_1D[i][3] /= iod.ref.rv.length;
  }

  input.close();

  // interpolation assuming 1D solution is centered on bubble_coord0
  double bubble_x0 = 0.0;
  double bubble_y0 = 0.0;
  double bubble_z0 = 0.0;
  double phi       = 0.0;
#pragma omp parallel for
  for (int iSub=0; iSub<domain->getNumLocSub(); ++iSub) {
    Vec<int> &fId((*fluidId)(iSub));
    SVec<double, 3> &x(X(iSub));

    double localRadius; int np;
    double localAlpha;
    for(int i=0; i<fluidId->size(); i++) {
      localRadius = sqrt((x[i][0]-bubble_x0)*(x[i][0]-bubble_x0)+(x[i][1]-bubble_y0)*(x[i][1]-bubble_y0)+(x[i][2]-bubble_z0)*(x[i][2]-bubble_z0));
      np = static_cast<int>(floor(localRadius/x_1D[numPoints-1]*(numPoints-1)));

      if(np>numPoints-1){
      //further away from the max radius of the 1D simulation, take last value
      //<==> constant extrapolation
        np = numPoints-1;
        phi = v_1D[np][3];
      }else{
      //linear interpolation
        localAlpha = (localRadius-x_1D[np])/(x_1D[np+1]-x_1D[np]);
        phi = localAlpha*v_1D[np][3] + (1.0-localAlpha)*v_1D[np+1][3];
      }
      // assumes that the fluidId==0 is outside flow
      // and that the fluidId==1 is the inside flow!!!
      int myId = phi>0.0 ? 1 : 0;
      if(fId[i]>0 && fId[i]!=myId)
        fprintf(stderr,"WARNING: Overwriting determined fluidId (%d) by one-dimensional result (%d).\n", fId[i], myId);
      fId[i] = myId;
    }
  }
}

//------------------------------------------------------------------------------

void FluidSelector::setupFluidIdMultiFluidInitialConditions(IoData &iod, DistSVec<double,3> &X)
{
  // Assumption: one geometric object surface does not cross another geometric
  //             object surface. In other words, there cannot exist a point
  //             in space where three or more fluids are in contact.
  // Assumption: a geometric surface always separates the main fluid
  //             given by varFcn[0] = iod.eqs.fluidModel =
  //             iod.eqs.fluidModelMap.dataMap.begin() from
  //             another fluid, at least one of the fluid involved is the main fluid

  if(!iod.mf.multiInitialConditions.planeMap.dataMap.empty()){
    map<int, PlaneData *>::iterator planeIt;
    for(planeIt  = iod.mf.multiInitialConditions.planeMap.dataMap.begin();
        planeIt != iod.mf.multiInitialConditions.planeMap.dataMap.end();
        planeIt++){
//      com->fprintf(stdout, "Processing initialization of LevelSet for plane %d paired with EOS %d\n", planeIt->first, planeIt->second->fluidModelID);
      double nx = planeIt->second->nx;
      double ny = planeIt->second->ny;
      double nz = planeIt->second->nz;
      double norm = sqrt(nx*nx+ny*ny+nz*nz);
      nx /= norm; ny /= norm; nz /= norm;

#pragma omp parallel for
      for (int iSub=0; iSub<domain->getNumLocSub(); ++iSub) {
        SVec<double,3> &x(X(iSub));
        Vec<int> &fId((*fluidId)(iSub));
        for (int i=0; i<x.size(); i++){
          double scalar = nx*(x[i][0]-planeIt->second->cen_x)+
                          ny*(x[i][1]-planeIt->second->cen_y)+
                          nz*(x[i][2]-planeIt->second->cen_z); // positive in the direction of n
          if(planeIt->second->fluidModelID > 0)
            if(scalar>0.0) {
              int myId = planeIt->second->fluidModelID;
              if(fId[i]>0 && fId[i]!=myId)
                continue; //avoid overwriting.
                //fprintf(stderr,"WARNING: Overwriting determined fluid Id (%d) by plane data (%d).\n", fId[i], myId);
              fId[i] = myId;
            }
        }
      }
    }
  }

  if(!iod.mf.multiInitialConditions.sphereMap.dataMap.empty()){
    map<int, SphereData *>::iterator sphereIt;
    for(sphereIt  = iod.mf.multiInitialConditions.sphereMap.dataMap.begin();
        sphereIt != iod.mf.multiInitialConditions.sphereMap.dataMap.end();
        sphereIt++){
//      com->fprintf(stdout, "Processing initialization of LevelSet for sphere %d paired with EOS %d\n", sphereIt->first, sphereIt->second->fluidModelID);
      double x0 = sphereIt->second->cen_x;
      double y0 = sphereIt->second->cen_y;
      double z0 = sphereIt->second->cen_z;
      double r0 = sphereIt->second->radius;
#pragma omp parallel for
      for (int iSub=0; iSub<domain->getNumLocSub(); ++iSub) {
        SVec<double,3>     &x  (X(iSub));
        Vec<int> &fId((*fluidId)(iSub));
        for (int i=0; i<x.size(); i++){
          double scalar = sqrt((x[i][0] - x0)*(x[i][0] - x0) +
                          (x[i][1] - y0)*(x[i][1] - y0) +
                          (x[i][2] - z0)*(x[i][2] - z0)) - r0;//positive outside the sphere
          if(sphereIt->second->fluidModelID > 0)
            if(scalar<0.0) {//inside the sphere
              int myId = sphereIt->second->fluidModelID;
              if(fId[i]>0 && fId[i]!=myId)
                continue; //avoid overwriting.
                //fprintf(stderr,"WARNING: Overwriting determined fluid Id (%d) by sphere data (%d).\n", fId[i], myId);
              fId[i] = myId;
            }
        }
      }
    }
  }

  if(!iod.mf.multiInitialConditions.cylinderMap.dataMap.empty()){
    map<int, CylinderData *>::iterator cylinderIt;
    for(cylinderIt  = iod.mf.multiInitialConditions.cylinderMap.dataMap.begin();
        cylinderIt != iod.mf.multiInitialConditions.cylinderMap.dataMap.end();
        cylinderIt++){
     if(cylinderIt->second->idModelID <= 0)
	continue;

#pragma omp parallel for
      for (int iSub=0; iSub<numLocSub; ++iSub) {
        SVec<double, 3> &x(X(iSub));

        Vec<int> &fId((*fluidId)(iSub));

        double scalar = 0.0;
        for(int i=0; i<fId.size(); i++) {
          scalar = cylinderIt->second->nx*(x[i][0] - cylinderIt->second->cen_x)+cylinderIt->second->ny*(x[i][1] - cylinderIt->second->cen_y)+cylinderIt->second->nz*(x[i][2] - cylinderIt->second->cen_z);
	  if (scalar < 0.0 || scalar > cylinderIt->second->L)
	    continue;
	  
	  double q[3] = {(x[i][0] - planeIt->second->cen_x) - scalar*cylinderIt->second->nx,
			 (x[i][1] - planeIt->second->cen_y) - scalar*cylinderIt->second->ny
			 (x[i][2] - planeIt->second->cen_z) - scalar*cylinderIt->second->nz};
	  
          if(sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]) < 
	     cylinderIt->second->r) { //node is on the same side indicated by vector

	    int myId = cylinderIt->second->fluidModelID;
	    if(fId[i]>0 && fId[i]!=myId)
	      continue; //avoid overwriting.

	    fId[i] = myId;
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

      double x0 = prismIt->second->cen_x;
      double y0 = prismIt->second->cen_y;
      double z0 = prismIt->second->cen_z;
      double wx = prismIt->second->w_x;
      double wy = prismIt->second->w_y;
      double wz = prismIt->second->w_z;
      
#pragma omp parallel for
      for (int iSub=0; iSub<domain->getNumLocSub(); ++iSub) {
        SVec<double,3>     &x  (X(iSub));
        Vec<int> &fId((*fluidId)(iSub));
        for (int i=0; i<x.size(); i++){

          if(prismIt->second->fluidModelID > 0)
            if(prismIt->second->inside(x[i][0],x[i][1],x[i][2])) {
              int myId = prismIt->second->fluidModelID;
              if(fId[i]>0 && fId[i]!=myId)
                continue; //avoid overwriting.
                //fprintf(stderr,"WARNING: Overwriting determined fluid Id (%d) by sphere data (%d).\n", fId[i], myId);
              fId[i] = myId;
            }
        }
      }
    }
  }
}
*/
//------------------------------------------------------------------------------

void FluidSelector::printFluidId(){
  int numLocSub = fluidId->numLocSub();
  int iSub;
#pragma omp parallel for
  for(iSub=0; iSub<numLocSub; ++iSub) {
    int    *tag = fluidId->subData(iSub);
    for (int i=0; i<fluidId->subSize(iSub); i++)
      fprintf(stdout, "fluidId[%d] = %d\n", i, tag[i]);
  }
}

//------------------------------------------------------------------------------

//template<int dim>
void FluidSelector::getFluidId(TriangulatedInterface* T){
//  assert(dim<=numPhases-1);
  int numLocSub = fluidId->numLocSub();
  int iSub;
#pragma omp parallel for
  for(iSub=0; iSub<numLocSub; ++iSub) {
    LevelSetStructure* LSS = T->getSubLSS(iSub);
    int     *tag       = fluidId->subData(iSub);
    int burnTag;
    for(int iNode=0; iNode<fluidId->subSize(iSub); iNode++){
      //if (programmedBurn && programmedBurn->isBurnedEOS(tag[iNode],burnTag))
      //	continue;
      tag[iNode] = 0;
/*      for(int i=0; i<dim; i++) {
        if(phi[iNode][i]>0.0) {
	  if (programmedBurn && (programmedBurn->isUnburnedEOS(i+1,burnTag) ||
				 programmedBurn->isBurnedEOS(i+1,burnTag)) ) {
	    if (programmedBurn->nodeInside(burnTag,iSub,iNode) || 
                programmedBurn->isFinished(burnTag))
	      tag[iNode] = programmedBurn->getBurnedEOS(burnTag);
	    else
	      tag[iNode] = programmedBurn->getUnburnedEOS(burnTag);
	    break;
	  }
	  else 
	    {  tag[iNode] = i+1; break; }

          
	}
      }
*/
      if (!LSS->isOccluded(0.0,iNode))
        tag[iNode] = LSS->fluidModel(0.0,iNode);
    }
  }
}

void FluidSelector::writeToDisk(const char* name) {

  DistSVec<int,1> i1(fluidId->info(), reinterpret_cast<int(*)[1]>(fluidId->data()));
  domain->writeVectorToFile(name, 0, 0.0, i1);

  if (fluidIdnm1){
    DistSVec<int,1> i2(fluidIdnm1->info(), reinterpret_cast<int(*)[1]>(fluidIdnm1->data()));
    domain->writeVectorToFile(name, 1, 0.0, i2);
  }

  if (fluidIdnm2){
    DistSVec<int,1> i3(fluidIdnm2->info(), reinterpret_cast<int(*)[1]>(fluidIdnm2->data()));
    domain->writeVectorToFile(name, 2, 0.0,i3);
  }

}
