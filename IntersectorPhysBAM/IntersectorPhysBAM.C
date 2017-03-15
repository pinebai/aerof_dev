#include <iostream>
#include <cstdio>
#include "IntersectorPhysBAM.h"
#include "FloodFill.h"
#include "Mpi_Utilities.h"
#include "Vector3D.h"
#include "Communicator.h"
#include <Domain.h>
#include <IoData.h>
#include <Vector.h>
#include <DistVector.h>
#include <Timer.h>
#include "FSI/CrackingSurface.h"
#include "LevelSet/LevelSetStructure.h"
#include "parser/Assigner.h"
#include <Connectivity.h>
#include <queue>

#include "TsRestart.h"

#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>

using std::pair;
using std::map;
using std::list;

using PhysBAM::TRIPLE;
using PhysBAM::GLOBAL_SUBD_ID;

typedef pair<int, int> iipair;
typedef pair<int, bool> ibpair;
typedef pair<iipair, ibpair> EdgePair;

const int IntersectorPhysBAM::UNDECIDED, IntersectorPhysBAM::INSIDE, IntersectorPhysBAM::OUTSIDE;
int IntersectorPhysBAM::OUTSIDECOLOR;

//int debug_PhysBAM_count = 0;
//----------------------------------------------------------------------------

DistIntersectorPhysBAM::DistIntersectorPhysBAM(IoData &iodata, 
															  Communicator *comm, 
															  int nNodes, double *xyz, 
															  int nElems, int (*abc)[3], 
															  CrackingSurface *cs) : DistLevelSetStructure(), iod(iodata)
{

  interface_thickness = iod.embed.interfaceThickness;

  externalSI = (iod.embed.surrogateinterface == EmbeddedFramework::EXTERNAL) ? true : false;

  if(externalSI) comm->fprintf(stdout, " +++ Using external-based surrogate interface +++ \n");

  withViscousTerms = false;
  if(iod.eqs.type == EquationsData::NAVIER_STOKES)
  {
	  withViscousTerms = true;

	  if(iod.bc.wall.type == BcsWallData::ISOTHERMAL) 
	  {
		  isoThermalwall = true;
		  Twall = iod.bc.wall.temperature;
	  }
	  else
	  {
		  isoThermalwall = false;
		  Twall = -1.0; // dummy
	  }
  }

  this->numFluid = iod.eqs.numPhase;
  floodFill=new FloodFill();
  com = comm;

  //get embedded structure surface mesh and restart pos
  char *struct_mesh, *struct_restart_pos;
  int sp = strlen(iod.input.prefix) + 1;

  struct_mesh        = new char[sp + strlen(iod.input.embeddedSurface)];
  sprintf(struct_mesh,"%s%s", iod.input.prefix, iod.input.embeddedSurface);

  if (iod.input.restart_file_package[0] == 0) {
    struct_restart_pos = new char[sp + strlen(iod.input.embeddedpositions)];
    if(iod.input.embeddedpositions[0] != 0)
      sprintf(struct_restart_pos,"%s%s", iod.input.prefix, iod.input.embeddedpositions);
    else //no restart position file provided
      strcpy(struct_restart_pos,""); 
  } else {
    struct_restart_pos = new char[256];
    char dummy[256], tmp[256];
    sprintf(tmp,"%s%s", iod.input.prefix, iod.input.restart_file_package);
    TsRestart::readRestartFileNames(tmp,
				    dummy,
				    dummy,
				    dummy,
				    dummy,
				    dummy,
				    dummy,
				    struct_restart_pos, comm);
				   
  }

  interpolatedNormal = (iod.embed.structNormal==EmbeddedFramework::NODE_BASED) ? true : false;
  
  //initialize the following to 0(NULL)
  physInterface = 0;
  triNorms = 0;
  triSize = 0;
  nodalNormal = 0;
  status0 = 0;
  closest = 0;
  occluded_node0 = 0;
  boxMin = 0;
  boxMax = 0;
  is_swept_helper = 0;

  rotOwn = 0;

  faceID = NULL;
  surfaceID = NULL;

  cracking = cs;
  gotNewCracking = false;

  with_sensitivity = false;

  //Load files. Compute structure normals. Initialize PhysBAM Interface
  if(nNodes && xyz && nElems && abc)
    init(nNodes, xyz, nElems, abc, struct_restart_pos);
  else {
    double XScale = (iod.problem.mode==ProblemData::NON_DIMENSIONAL) ? 1.0 : iod.ref.rv.length;
    init(struct_mesh, struct_restart_pos, XScale);
  }
  setStructureType();
  setPorosity();
  setActuatorDisk();
  makerotationownership();
  updatebc();

  comm->barrier();

  delete[] struct_mesh;
  if(iod.input.embeddedpositions[0] != 0) delete[] struct_restart_pos;
}

//----------------------------------------------------------------------------

DistIntersectorPhysBAM::~DistIntersectorPhysBAM() 
{
  if(Xs)          delete[] Xs;
  if(Xs0)         delete[] Xs0;
  if(Xs_n)        delete[] Xs_n;
  if(Xs_np1)      delete[] Xs_np1;
  if(Xsdot)       delete[] Xsdot;
  if(status0)     delete   status0;
  if(closest)     delete   closest;
  delete is_swept_helper;
  if(triSize)     delete[] triSize;
  if(triNorms)    delete[] triNorms;
  if(nodalNormal) delete[] nodalNormal;
  if(boxMax)      delete   boxMax;
  if(boxMin)      delete   boxMin;
  delete occluded_node0;

  for(int i = 0; i < numLocSub; ++i) {
      delete intersector[i];}
  delete []intersector;
  delete physInterface;
  delete floodFill;

  if(stElem) delete[] stElem;
  if(solidX) delete solidX;
  if(solidX0) delete solidX0;
  if(solidXn) delete solidXn;
  if(solidXnp1) delete solidXnp1;
  if(solidXdS) delete solidXdS;
  if(faceID) delete[] faceID;
  if(surfaceID) delete[] surfaceID;
  if(porosity) delete[] porosity;
  if(rotOwn) delete[] rotOwn;

  if(dXdSb) delete[] dXdSb;

}

//----------------------------------------------------------------------------

LevelSetStructure &
DistIntersectorPhysBAM::operator()(int subNum) const {
  return *intersector[subNum];
}

//----------------------------------------------------------------------------

/** Intersector initialization method
*
* \param dataTree the data read from the input file for this intersector.
*/
void DistIntersectorPhysBAM::init(char *solidSurface, char *restartSolidSurface, double XScale) {

  // Read data from the solid surface input file.
  FILE *topFile;
  topFile = fopen(solidSurface, "r");
  if (topFile == NULL) {
	  com->fprintf(stderr, "Embedded structure surface mesh doesn't exist: %s, %s\n",solidSurface,restartSolidSurface); 
	  exit(1); 
  }

  char line[MAXLINE],key1[MAXLINE],key2[MAXLINE];

  // load the nodes and initialize all node-based variables.

  // load solid nodes at t=0
  int num0 = 0; 
  int num1 = 0;

  double x1,x2,x3;
  int node1, node2, node3;

  int type_read = 0;
  int surfaceid = 0;

  std::list<int> indexList;
  std::list<Vec3D> nodeList;
  std::list<int> elemIdList;    
  std::list<int> elemList1;
  std::list<int> elemList2;
  std::list<int> elemList3;
  std::list<int> surfaceIDList;

  int maxIndex = 0, maxElem = -1;

  while (fgets(line,MAXLINE,topFile) != 0) {
    sscanf(line, "%s", key1);
    bool skip = false;
    if (strcmp(key1, "Nodes") == 0) {
      sscanf(line, "%*s %s", key2);
      skip = true;
      type_read = 1;
    }
    if (strcmp(key1, "Elements") == 0) {
      sscanf(line, "%*s %s", key2);
      skip = true;
      type_read = 2;
      int underscore_pos = -1;
      int k = 0;
      while((key2[k] != '\0') && (k<MAXLINE)) {
        if(key2[k] == '_') underscore_pos = k;
        k++;
      }
      if(underscore_pos > -1)
        sscanf(key2+(underscore_pos+1),"%d",&surfaceid); 
    }
    if (!skip) {
      if (type_read == 1) {
	sscanf(line, "%d %lf %lf %lf", &num1, &x1, &x2, &x3);
			 if(num1<1) {
				 com->fprintf(stderr,"ERROR: detected a node with index %d in the embedded surface file!\n",num1); 
				 exit(-1);
			 }

        x1 /= XScale;
        x2 /= XScale;
        x3 /= XScale;
        indexList.push_back(num1);
        if(num1>maxIndex) maxIndex = num1;
        nodeList.push_back(Vec3D(x1,x2,x3));
      }
      if (type_read == 2) {
	sscanf(line,"%d %d %d %d %d", &num0, &num1, &node1, &node2, &node3);
        elemList1.push_back(node1-1);
        elemList2.push_back(node2-1);
        elemList3.push_back(node3-1);
	surfaceIDList.push_back(surfaceid);
      }
    }
  }

  numStNodes = nodeList.size();
  if(numStNodes != maxIndex) {
    com->fprintf(stderr,"ERROR: The node set of the embedded surface have gap(s). \n");
    com->fprintf(stderr,"       Detected max index = %d, number of nodes = %d\n", maxIndex, numStNodes);
    com->fprintf(stderr,"NOTE: Currently the node set of the embedded surface cannot have gaps. Moreover, the index must start from 1.\n");
    exit(-1);
  }

  numStElems = elemList1.size();

  // feed data to Xss. 
  Xs      = new Vec3D[numStNodes];
  Xs0     = new Vec3D[numStNodes];
  Xs_n    = new Vec3D[numStNodes];
  Xs_np1  = new Vec3D[numStNodes];
  Xsdot   = new Vec3D[numStNodes];
  solidX  = new Vec<Vec3D>(numStNodes, Xs);
  solidX0 = new Vec<Vec3D>(numStNodes, Xs0);
  solidXn = new Vec<Vec3D>(numStNodes, Xs_n);
  solidXnp1 = new Vec<Vec3D>(numStNodes, Xs_np1);

     dXdSb = new Vec3D[numStNodes];
  solidXdS = new Vec<Vec3D>(numStNodes, dXdSb);

  faceID = new int[numStElems];

  surfaceID = new int[numStNodes];

  std::list<Vec3D>::iterator it1;
  std::list<int>::iterator it2;
  std::list<int>::iterator it_1;
  std::list<int>::iterator it_2;
  std::list<int>::iterator it_3;
  std::list<int>::iterator it_4;

  it2 = indexList.begin();
  for (it1=nodeList.begin(); it1!=nodeList.end(); it1++) {
    Xs[(*it2)-1] = *it1;
    it2++;
  }

  for (int k=0; k<numStNodes; k++) {
    Xs0[k]    = Xs[k];
    Xs_n[k]   = Xs[k];
    Xs_np1[k] = Xs[k];
    Xsdot[k]  = Vec3D(0.0, 0.0, 0.0);
    surfaceID[k] = 0;
    dXdSb[k]     = 0; //init
  }

  stElem = new int[numStElems][3];

  it_1 = elemList1.begin();
  it_2 = elemList2.begin();
  it_3 = elemList3.begin();
  it_4 = surfaceIDList.begin();
  for (int i=0; i<numStElems; i++) {
    stElem[i][0] = *it_1;
    stElem[i][1] = *it_2;
    stElem[i][2] = *it_3;
    faceID[i] = *it_4;
    surfaceID[*it_1] = *it_4;
    surfaceID[*it_2] = *it_4;
    surfaceID[*it_3] = *it_4;
    it_1++;
    it_2++;
    it_3++;
    it_4++;
  }

  fclose(topFile);

  int nInputs;
  char c1[200], c2[200], c3[200];
  // load solid nodes at restart time.
  if (restartSolidSurface[0] != 0) {
    FILE* resTopFile = fopen(restartSolidSurface, "r");

    if(resTopFile==NULL) {
		 com->fprintf(stderr, "restart topFile doesn't exist: \"%s\".\n",restartSolidSurface); 
		 exit(1);
	 }

    int ndMax = 0, ndMax2 = 0;
    std::list<std::pair<int,Vec3D> > nodeList2;
    std::list<std::pair<int,Vec3D> >::iterator it2;

    while(1) {
      nInputs = fscanf(resTopFile,"%s", c1);
      if(nInputs!=1) break;    
      char *endptr;
      num1 = strtol(c1, &endptr, 10);
      if(endptr == c1) break;

      int toto = fscanf(resTopFile,"%lf %lf %lf\n", &x1, &x2, &x3);
      nodeList2.push_back(std::pair<int,Vec3D>(num1,Vec3D(x1,x2,x3)));
      ndMax = std::max(num1, ndMax);
    }
    if ((!cracking && ndMax!=numStNodes) || (cracking && ndMax != cracking->totNodes())) {
      com->fprintf(stderr,"ERROR: number of nodes in restart top-file is wrong; ndMax = %d; numStNodes = %d\n", ndMax , numStNodes);
      exit(1);
    }

    for (int k=0; k<numStNodes; k++)
      Xs[k] = Vec3D(0,0,0);
    
    for (it2=nodeList2.begin(); it2!=nodeList2.end(); it2++)
      Xs[it2->first-1] = it2->second;

    for (int k=0; k<numStNodes; k++) {
      Xs_n[k]         = Xs[k];
      Xs_np1[k]    = Xs[k];
    }
    fclose(resTopFile);
  }

  totStNodes = numStNodes;
  totStElems = numStElems;

  // Verify (1)triangulated surface is closed (2) normal's of all triangles point outward.
//  com->fprintf(stderr,"- IntersectorPhysBAM: Checking the embedded structure surface...   ");
//  if (checkTriangulatedSurface()) 
//    com->fprintf(stderr,"Ok.\n");
//  else {
//    com->fprintf(stderr,"\n"); 
//    exit(-1); 
//  }

  initializePhysBAM();
}

//----------------------------------------------------------------------------
/** Intersector initialization method
*
* \param dataTree the data read from the input file for this intersector.
*/
void DistIntersectorPhysBAM::init(int nNodes, double *xyz, int nElems, int (*abc)[3], char *restartSolidSurface) {

  // node set
  if(cracking && nNodes!=cracking->usedNodes()) { 
	  com->fprintf(stderr,"SOFTWARE BUG (nNodes = %d; usedNodes = %d)!\n", nNodes, cracking->usedNodes());
	  exit(-1);
  }

  numStNodes = nNodes;
  totStNodes = cracking ? cracking->totNodes() : numStNodes;

  // feed data to Xss. 
  Xs      = new Vec3D[totStNodes];
  Xs0     = new Vec3D[totStNodes];
  Xs_n    = new Vec3D[totStNodes];
  Xs_np1  = new Vec3D[totStNodes];
  Xsdot   = new Vec3D[totStNodes];
  solidX  = new Vec<Vec3D>(totStNodes, Xs);
  solidX0 = new Vec<Vec3D>(totStNodes, Xs0);
  solidXn = new Vec<Vec3D>(totStNodes, Xs_n);
  solidXnp1 = new Vec<Vec3D>(totStNodes, Xs_np1);

     dXdSb = new Vec3D[numStNodes];
  solidXdS = new Vec<Vec3D>(numStNodes, dXdSb);

  for (int k=0; k<numStNodes; k++) {
    Xs[k]     = Vec3D(xyz[3*k], xyz[3*k+1], xyz[3*k+2]);
    Xs0[k]    = Xs[k];
    Xs_n[k]   = Xs[k];
    Xs_np1[k] = Xs[k];
    Xsdot[k]  = Vec3D(0.0, 0.0, 0.0);    
    dXdSb[k]  = 0; //init
  }

  // elem set
  if(cracking && nElems!=cracking->usedTrias()) {com->fprintf(stderr,"SOFTWARE BUG 2!\n");exit(-1);}

  numStElems = nElems;
  totStElems = cracking ? cracking->totTrias() : numStElems;

  stElem = new int[totStElems][3];
  for (int i=0; i<numStElems; i++)
    for (int k=0; k<3; k++)
      stElem[i][k] = abc[i][k];

  // load solid nodes at restart time.
  if (restartSolidSurface[0] != 0) {
    FILE* resTopFile = fopen(restartSolidSurface, "r");
    if(resTopFile==NULL) {com->fprintf(stderr, "restart topFile doesn't exist.\n"); exit(1);}
    int ndMax2 = 0;
    std::list<std::pair<int,Vec3D> > nodeList2;
    std::list<std::pair<int,Vec3D> >::iterator it2;

    int ndMax = 0;
    while(1) {
      int nInputs, num1;
      double x1, x2, x3;
      char *endptr, c1[200];

      nInputs = fscanf(resTopFile,"%s", c1);
      if(nInputs!=1) break;
      num1 = strtol(c1, &endptr, 10);
      if(endptr == c1) break;

      int toto = fscanf(resTopFile,"%lf %lf %lf\n", &x1, &x2, &x3);
      nodeList2.push_back(std::pair<int,Vec3D>(num1,Vec3D(x1,x2,x3)));
      ndMax = std::max(num1, ndMax);
    }
    if ((!cracking && ndMax!=numStNodes) || (cracking && ndMax != cracking->totNodes())) {
      com->fprintf(stderr,"ERROR: number of nodes in restart top-file is wrong; ndMax = %d; numStNodes = %d\n", ndMax , numStNodes);
      exit(1);
    }

    for (int k=0; k<numStNodes; k++)
      Xs[k] = Vec3D(0,0,0);

    for (it2=nodeList2.begin(); it2!=nodeList2.end(); it2++)
      Xs[it2->first-1] = it2->second;

    for (int k=0; k<numStNodes; k++) {
      Xs_n[k]   = Xs[k];
      Xs_np1[k] = Xs[k];
    }
    fclose(resTopFile);
  }

  // Verify (1)triangulated surface is closed (2) normal's of all triangles point outward.
//  com->fprintf(stderr,"- IntersectorPhysBAM: Checking the embedded structure surface...   ");
//  if (checkTriangulatedSurface())
//    com->fprintf(stderr,"Ok.\n");
//  else {
//    com->fprintf(stderr,"\n"); 
//    exit(-1); 
//  }

  initializePhysBAM();
}

//----------------------------------------------------------------------------
void DistIntersectorPhysBAM::setPorosity() {
  map<int,SurfaceData *> &surfaceMap = iod.surfaces.surfaceMap.dataMap;
  map<int,BoundaryData *> &bcMap = iod.bc.bcMap.dataMap;

  porosity = new double[numStElems];
  for(int i=0; i<numStElems; i++) {
    porosity[i] = 0.;
  }

  if(faceID) {
    for(int i=0; i<numStElems; i++) {
      map<int,SurfaceData*>::iterator it = surfaceMap.find(faceID[i]);
      if (it != surfaceMap.end()) {
        map<int,BoundaryData *>::iterator it2 = bcMap.find(it->second->bcID);
        if(it2 != bcMap.end()) { // the bc data have been defined
          if(it2->second->type == BoundaryData::POROUSWALL ) {
            porosity[i] = it2->second->porosity;
          }
        }
      }
    }
  }
}
//----------------------------------------------------------------------------
void DistIntersectorPhysBAM::setStructureType() {
  map<int,SurfaceData *> &surfaceMap = iod.surfaces.surfaceMap.dataMap;
  map<int,BoundaryData *> &bcMap = iod.bc.bcMap.dataMap;

  structureType = new int[numStElems];
  for(int i=0; i<numStElems; i++) {
    structureType[i] = 1;
  }

  if(faceID) {
    for(int i=0; i<numStElems; i++) {
      map<int,SurfaceData*>::iterator it = surfaceMap.find(faceID[i]);
      if (it != surfaceMap.end()) {
        map<int,BoundaryData *>::iterator it2 = bcMap.find(it->second->bcID);
        if(it2 != bcMap.end()) { // the bc data have been defined
          structureType[i] = it2->second->type;
        }
      }
    }
  }
}

//----------------------------------------------------------------------------
void DistIntersectorPhysBAM::setActuatorDisk() {
    map<int,SurfaceData *> &surfaceMap = iod.surfaces.surfaceMap.dataMap;
    map<int,BoundaryData *> &bcMap = iod.bc.bcMap.dataMap;
    actuatorDiskMethod = new int[numStElems];
    actuatorDiskPressureJump = new double[numStElems];
    actuatorDiskReconstructionMethod = new int[numStElems];
    isCorrectedMethod = new bool[numStElems];
    for(int i=0; i<numStElems; i++) {
        actuatorDiskMethod[i] = 1;
        actuatorDiskPressureJump[i] = 0.0;
        actuatorDiskReconstructionMethod[i] = -1;
        isCorrectedMethod[i] = false;
    }
    gamma = iod.eqs.fluidModel.gasModel.specificHeatRatio;
    if(faceID) {
        for(int i=0; i<numStElems; i++) {
            map<int,SurfaceData*>::iterator it = surfaceMap.find(faceID[i]);
            if (it != surfaceMap.end()) {
                map<int,BoundaryData *>::iterator it2 = bcMap.find(it->second->bcID);
                if(it2 != bcMap.end()) { // the bc data have been defined
                    if(it2->second->type == BoundaryData::ACTUATORDISK ) {
                        actuatorDiskPressureJump[i] = it2->second->pressureJump;
                        if (iod.problem.mode == ProblemData::DIMENSIONAL)
                            actuatorDiskPressureJump[i]/= iod.ref.rv.pressure;

                        //--------------------------------
                        if(it2->second->actuatorDiskMethod == BoundaryData::SOURCETERM){
                            actuatorDiskMethod[i] = 1;
                        }
                        else if(it2->second->actuatorDiskMethod == BoundaryData::RIEMANNSOLVER){
                            actuatorDiskMethod[i] = 2;
                        }else{
                            com->fprintf(stderr, "!!! WARNING: no actuator disk method specified, defaulting to SOURCETERM\n\n");
                            actuatorDiskMethod[i] = 1;
                        }

                        //--------------------------------
                        if(it2->second->velocityReconstructionMethod == BoundaryData::AVERAGE){
                            actuatorDiskReconstructionMethod[i] = 1;
                        }
                        else if(it2->second->velocityReconstructionMethod == BoundaryData::FIRSTORDER){
                            actuatorDiskReconstructionMethod[i] = 2;
                        }
                        else if(it2->second->velocityReconstructionMethod == BoundaryData::SECONDORDER){
                            actuatorDiskReconstructionMethod[i] = 3;
                        }
                            //--------------------------------
                        else{
                            com->fprintf(stderr, "!!! WARNING: no actuator disk method specified, defaulting to Average\n\n");
                            actuatorDiskReconstructionMethod[i] = 1;
                        }
                        if(it2->second->sourceTermExpression==BoundaryData::OLD){
                            isCorrectedMethod[i] = false;
                        }
                        else if(it2->second->sourceTermExpression==BoundaryData::CORRECTED){
                            isCorrectedMethod[i] = true;
                        }
                        else{
                            com->fprintf(stderr, "!!! WARNING: no Source term method specified for the actuator Disk, defaulting to Corrected\n\n");
                            isCorrectedMethod[i] = true;
                        }
                    }
                }
            }
        }
    }
}
//----------------------------------------------------------------------------



void DistIntersectorPhysBAM::setdXdSb(int N, double* dxdS, double* dydS, double* dzdS){

  if(N != numStNodes) {
    fprintf(stderr, "Error num nodes: %d and %d\n", N,  numStNodes);
    exit(-1);
  }

  with_sensitivity = true;
  //dXdSb = new Vec3D[numStNodes];

  double Normx=0.0, Normy=0.0, Normz=0.0;

  for(int i=0; i<N; ++i){
    dXdSb[i][0] = dxdS[i];
    dXdSb[i][1] = dydS[i];
    dXdSb[i][2] = dzdS[i];

    Normx += fabs(dxdS[i]);
    Normy += fabs(dydS[i]);
    Normz += fabs(dzdS[i]);
  }

  if(Normx==0 && Normy==0 && Normz==0)
    com->fprintf(stderr, "!!! WARNING: No surface Sensitivity Perturbation\n\n");

}

//----------------------------------------------------------------------------
void DistIntersectorPhysBAM::makerotationownership() {
  map<int,SurfaceData *> &surfaceMap = iod.surfaces.surfaceMap.dataMap;
  map<int,RotationData*> &rotationMap= iod.rotations.rotationMap.dataMap;
  map<int,SurfaceData *>::iterator it = surfaceMap.begin();
  rotOwn = 0;
  
  int numRotSurfs = 0;
  int numTransWalls = 0;
  while(it != surfaceMap.end()) {
    map<int,RotationData*>::iterator it1 = rotationMap.find(it->second->rotationID);
    if(it1!=rotationMap.end()) {
      if(it1->second->infRadius) {
        numTransWalls++;
        com->fprintf(stderr," ... surface %2d is ``translating''\n",it->first, it1->first);
        com->fprintf(stderr,"     -> uniform velocity V = %3.2e in direction %3.2e %3.2e %3.2e\n",
                     it1->second->omega, it1->second->nx,it1->second->ny,it1->second->nz);
      } else {
        numRotSurfs++;
        com->fprintf(stderr," ... surface %2d is ``rotating'' using rotation data %2d\n",it->first, it1->first);
        com->fprintf(stderr,"     -> omega = %3.2e, rotation axis = %3.2e %3.2e %3.2e\n",
                     it1->second->omega, it1->second->nx,it1->second->ny,it1->second->nz);
      }
    }
    it++;
  }

  if(numRotSurfs || numTransWalls) {
    rotOwn = new int[numStNodes];

    for (int k=0; k<numStNodes; k++) {
      rotOwn[k] = -1;

      map<int,SurfaceData *>::iterator it = surfaceMap.find(surfaceID[k]);
      if(it != surfaceMap.end()) {
         int rotID = it->second->rotationID; // = -1 if not defined in input file
         if ((rotOwn[k] != -1 && rotOwn[k] != rotID) ||
             (rotOwn[k] != -1 && rotOwn[k] != rotID) ||
             (rotOwn[k] != -1 && rotOwn[k] != rotID))  {
           com->fprintf(stderr, " ... WARNING: Embedded Node %d associated to more than 1 Rotation ID\n",k);
         }
         rotOwn[k] = rotID;
         rotOwn[k] = rotID;
         rotOwn[k] = rotID;
      }
    }
  }
}
//----------------------------------------------------------------------------
void DistIntersectorPhysBAM::updatebc() {
  map<int,RotationData*> &rotationMap= iod.rotations.rotationMap.dataMap;
  if (rotOwn)  { 
    for (int k=0; k<numStNodes; k++) {
      if (rotOwn[k]>=0) {    // node belongs to a (potential) "rotating" surface
        map<int,RotationData *>::iterator it =  rotationMap.find(rotOwn[k]);
        if(it != rotationMap.end()) { // the rotation data have been defined
	  if(it->second->infRadius == RotationData::TRUE) {
            double vel = it->second->omega / iod.ref.rv.velocity;
	    Xsdot[k][0] = vel*it->second->nx + Xsdot[k][0];
	    Xsdot[k][1] = vel*it->second->ny + Xsdot[k][1];
	    Xsdot[k][2] = vel*it->second->nz + Xsdot[k][2];
	  } 
          else {
            double XScale = (iod.problem.mode==ProblemData::NON_DIMENSIONAL) ? 1.0 : iod.ref.rv.length;
	    double tref = iod.ref.rv.time;
	    double ox = tref*it->second->omega*it->second->nx;
	    double oy = tref*it->second->omega*it->second->ny;
	    double oz = tref*it->second->omega*it->second->nz;
	    double xd = Xs[k][0] - it->second->x0 / XScale;
	    double yd = Xs[k][1] - it->second->y0 / XScale;
	    double zd = Xs[k][2] - it->second->z0 / XScale;
	    Xsdot[k][0] = oy*zd-oz*yd + Xsdot[k][0];
	    Xsdot[k][1] = oz*xd-ox*zd + Xsdot[k][1];
	    Xsdot[k][2] = ox*yd-oy*xd + Xsdot[k][2];
	  }
	} 
        else  { // no rotation data -> use velocity from mesh motion if any
          Xsdot[k][0] = Xsdot[k][0];
          Xsdot[k][1] = Xsdot[k][1];
          Xsdot[k][2] = Xsdot[k][2];
	}
      }
      else  { // no rotation data -> use velocity from mesh motion if any
        Xsdot[k][0] = Xsdot[k][0];
        Xsdot[k][1] = Xsdot[k][1];
        Xsdot[k][2] = Xsdot[k][2];
      }
    }
  }
}
//----------------------------------------------------------------------------
//NOTE: In PhysBAM array index starts from 1 instead of 0
void DistIntersectorPhysBAM::initializePhysBAM() { 

  // Initialize the Particles list
  PhysBAM::GEOMETRY_PARTICLES<PhysBAM::VECTOR<double,3> > *physbam_solids_particle=new PhysBAM::GEOMETRY_PARTICLES<PhysBAM::VECTOR<double,3> >();
  physbam_solids_particle->array_collection.Resize(numStNodes);

	for (int i=0; i<numStNodes; i++) {
		physbam_solids_particle->X(i+1) = PhysBAM::VECTOR<double,3>(Xs[i][0],Xs[i][1], Xs[i][2]);
	}
  
  // Initialize the Triangle list.

  // Speed improvement.  When we are doing cracking,
  // we do not add phantom triangles to the triangle hierarchy.
  // Added by Alex Main (October 2013)
  //
  int np = (cracking?cracking->numberRealTriangles():numStElems);
  PhysBAM::ARRAY<PhysBAM::VECTOR<int,3> > physbam_stElem(np);
  if (cracking) {
    for (int i=0; i<np; i++)  {
      int id = cracking->mapTriangleID(i);
      physbam_stElem(i+1) = PhysBAM::VECTOR<int,3>(stElem[id][0]+1,stElem[id][1]+1,stElem[id][2]+1);
    }
  } else {
    for (int i=0; i<np; i++)  {
      physbam_stElem(i+1) = PhysBAM::VECTOR<int,3>(stElem[i][0]+1,stElem[i][1]+1,stElem[i][2]+1);
    }
  }

  // Initialize the mesh.
  PhysBAM::TRIANGLE_MESH *mesh = new PhysBAM::TRIANGLE_MESH(numStNodes,physbam_stElem);
	mesh->Initialize_Adjacent_Elements();
	mesh->Set_Number_Nodes(numStNodes);

  // Construct TRIANGULATED_SURFACE.
  if(physInterface) delete physInterface;
  physInterface = new PhysBAMInterface<double>(*mesh,*physbam_solids_particle,cracking);
  com->fprintf(stderr,"Setting interface thickness to %e.\n", interface_thickness);
  physInterface->SetThickness(interface_thickness);

  with_sensitivity = false;

}

//----------------------------------------------------------------------------

EdgePair DistIntersectorPhysBAM::makeEdgePair(int node1, int node2, int triangleNumber) {
if(node1 < node2)
 return EdgePair(iipair(node1, node2), ibpair(triangleNumber, true));
else
 return EdgePair(iipair(node2, node1), ibpair(triangleNumber, false));
}

//----------------------------------------------------------------------------

bool DistIntersectorPhysBAM::checkTriangulatedSurface()
{
  map<iipair, ibpair> edgeMap;

  for (int iTriangle=0; iTriangle<numStElems; iTriangle++) {
    int from1, to1, from2, to2, from3, to3;
    bool found1, found2, found3;
    from1 = stElem[iTriangle][0];  to1 = stElem[iTriangle][1];  found1 = false;
    from2 = stElem[iTriangle][1];  to2 = stElem[iTriangle][2];  found2 = false;
    from3 = stElem[iTriangle][2];  to3 = stElem[iTriangle][0];  found3 = false;

    EdgePair ep[3];
    ep[0] = makeEdgePair(from1, to1, iTriangle);
    ep[1] = makeEdgePair(from2, to2, iTriangle);
    ep[2] = makeEdgePair(from3, to3, iTriangle);

    for(int i=0; i < 3; ++i) {
      map<iipair, ibpair>::iterator it = edgeMap.find(ep[i].first);
      if(it != edgeMap.end()) { // we found this edge
         if(it->second.second == ep[i].second.second)
           {com->fprintf(stderr,"ERROR: surface is not closed or a triangle orientation problem. exit.\n"); return false;}
         else {
             int oTriangle = it->second.first;
             int n1 = it->second.second ? ep[i].first.first : ep[i].first.second;
             int edgeIndex;
             if(stElem[oTriangle][0] == n1)
               edgeIndex = 0;
             else if(stElem[oTriangle][1] == n1)
               edgeIndex = 1;
             else
               edgeIndex = 2;
         }
      } else
        edgeMap[ep[i].first] = ep[i].second;
    }
  }
  return true;
}

//----------------------------------------------------------------------------

void DistIntersectorPhysBAM::buildSolidNormals() {

  if(!triNorms) triNorms = new Vec3D[totStElems];
  if(!triSize)  triSize = new double[totStElems];
  if(interpolatedNormal) {
    if(!nodalNormal) nodalNormal = new Vec3D[totStNodes];
    for(int i=0; i<numStNodes; i++) nodalNormal[i] = 0.0;
  }

  // Also look to determine a point inside the solid but away from the structure.
  double nrmMax = 0;
  int trMaxNorm = -1;
  for (int iTriangle=0; iTriangle<numStElems; iTriangle++) {
    int n1 = stElem[iTriangle][0];
    int n2 = stElem[iTriangle][1];
    int n3 = stElem[iTriangle][2];
    double x1 = Xs[n1][0];
    double y1 = Xs[n1][1];
    double z1 = Xs[n1][2];
    double dx2 = Xs[n2][0]-x1;
    double dy2 = Xs[n2][1]-y1;
    double dz2 = Xs[n2][2]-z1;
    double dx3 = Xs[n3][0]-x1;
    double dy3 = Xs[n3][1]-y1;
    double dz3 = Xs[n3][2]-z1;

   // first calculate the "size" of the triangle.
    double dx23 = Xs[n3][0] - Xs[n2][0];
    double dy23 = Xs[n3][1] - Xs[n2][1];
    double dz23 = Xs[n3][2] - Xs[n2][2];
    double size12 = Vec3D(dx2,dy2,dz2).norm();
    double size13 = Vec3D(dx3,dy3,dz3).norm();
    double size23 = Vec3D(dx23,dy23,dz23).norm();
    triSize[iTriangle] = min(size12,min(size13,size23));

    // now calculate the normal.
    triNorms[iTriangle] = Vec3D(dx2, dy2, dz2)^Vec3D(dx3,dy3,dz3);
    
    if(interpolatedNormal){ // compute nodal normal (weighted by area)
      nodalNormal[n1] += triNorms[iTriangle];
      nodalNormal[n2] += triNorms[iTriangle];
      nodalNormal[n3] += triNorms[iTriangle];
    }
    
    double nrm = triNorms[iTriangle].norm();
    if(nrm > nrmMax) {
      nrmMax = nrm;
      trMaxNorm = iTriangle;
    }
    // normalize the normal.
    if(nrm > 0.0) triNorms[iTriangle] /= nrm;
    else fprintf(stderr,"ERROR: Tri %d: (%d(%e,%e,%e), %d(%e,%e,%e), %d(%e,%e,%e)\n", iTriangle+1, n1, x1,y1,z1, n2, dx2,dy2,dz2, n3, dx3,dy3,dz3);

  }

  if(interpolatedNormal) //normalize nodal normals.
    for(int i=0; i<numStNodes; i++) {
      nodalNormal[i] /= nodalNormal[i].norm();
    }
}

//----------------------------------------------------------------------------


/** compute the intersections, node statuses and normals for the initial geometry */
void DistIntersectorPhysBAM::initialize(Domain *d, DistSVec<double,3> &X, 
													 DistSVec<double,3> &Xn, IoData &iod, 
													 DistVec<int> *point_based_id, 
													 DistVec<int>* oldStatus /*for restart*/) {

  if(this->numFluid<1) {
    fprintf(stderr,"ERROR: numFluid = %d!\n", this->numFluid);
    exit(-1);
  }
  this->X = &X;
  this->Xn = &Xn;
  domain = d;
  numLocSub = d->getNumLocSub();
  intersector = new IntersectorPhysBAM*[numLocSub];

  status0 = new DistVec<int>(domain->getNodeDistInfo());
  occluded_node0 = new DistVec<bool>(domain->getNodeDistInfo());
  is_swept_helper = new DistVec<int>(domain->getNodeDistInfo());
  boxMin = new DistSVec<double,3>(domain->getNodeDistInfo());
  boxMax = new DistSVec<double,3>(domain->getNodeDistInfo());
  closest = new DistVec<ClosestPoint>(domain->getNodeDistInfo()); //needed only for multi-phase cracking.
  status = new DistVec<int>(domain->getNodeDistInfo());
  distance = new DistVec<double>(domain->getNodeDistInfo());
  is_swept = new DistVec<bool>(domain->getNodeDistInfo());
  is_active = new DistVec<bool>(domain->getNodeDistInfo());
  is_occluded = new DistVec<bool>(domain->getNodeDistInfo());
  edge_intersects = new DistVec<bool>(domain->getEdgeDistInfo());

	edge_SI  = new DistVec<bool>(domain->getEdgeDistInfo());

	   xi_SI = new DistVec<double>(domain->getEdgeDistInfo());
	  eta_SI = new DistVec<double>(domain->getEdgeDistInfo());
	TriID_SI = new DistVec<int>(domain->getEdgeDistInfo());
	nWall_SI = new DistVec<Vec3D>(domain->getEdgeDistInfo());

	   xi_node = new DistVec<double>(domain->getNodeDistInfo());
     eta_node = new DistVec<double>(domain->getNodeDistInfo());
	TriID_node = new DistVec<int>(domain->getNodeDistInfo());
	nWall_node = new DistVec<Vec3D>(domain->getNodeDistInfo());

  // for hasCloseTriangle
  DistVec<bool> tId(domain->getNodeDistInfo());

#pragma omp parallel for
  for(int i = 0; i < numLocSub; ++i)
		intersector[i] = new IntersectorPhysBAM( *(d->getSubDomain()[i]), X(i), 
									  					    (*status0)(i), (*closest)(i),
                                            (*occluded_node0)(i), *this);

  updatePhysBAMInterface(Xs, numStNodes,X,true,false);

  buildSolidNormals();
  d->findNodeBoundingBoxes(X,*boxMin,*boxMax);

  int numIntersectedEdges=0;
#pragma omp parallel for
  for(int i = 0; i < numLocSub; ++i) {
    intersector[i]->hasCloseTriangle(X(i), Xn(i), (*boxMin)(i), (*boxMax)(i), tId(i));

    numIntersectedEdges += intersector[i]->findIntersections(X(i), tId(i), *com);
  }

  list< pair<Vec3D,int> > points; //pair points with fluid model ID (or point id if needed).
  map<int,int> pid2id; //maps point id to id

  if(!iod.embed.embedIC.pointMap.dataMap.empty()){
    map<int, PointData *>::iterator pointIt;
    int count = 0;
    for(pointIt  = iod.embed.embedIC.pointMap.dataMap.begin();
      pointIt != iod.embed.embedIC.pointMap.dataMap.end();
      pointIt ++){
      int myID = pointIt->second->fluidModelID;
      int myPID = ++count;
      Vec3D xyz(pointIt->second->x, pointIt->second->y,pointIt->second->z);

      if(point_based_id) {
        points.push_back(pair<Vec3D,int>(xyz, myPID));
        pid2id[myPID] = myID;
		  } else {
        points.push_back(pair<Vec3D,int>(xyz, myID));
		  }

      if(myID>=numFluid) { //myID should start from 0
        com->fprintf(stderr,"ERROR:FluidModel %d doesn't exist! NumPhase = %d\n", myID, numFluid);
			  exit(-1);
		  }
	  }

  } else if (iod.embed.embedIC.dummyPointMap.dataMap.empty()) {
    com->fprintf(stderr, "Point-based initial conditions could not be found.  Assuming single-phase flow\n");
  }

  if (iod.mf.levelSetMethod == MultiFluidData::TRIANGULATED) {

    points.push_back(pair<Vec3D,int>(Vec3D(0.0,0.0,0.0),0));
    points.push_back(pair<Vec3D,int>(Vec3D(1.0,0.0,0.0),1));
  }

  // If we have the fluid ID from a restart file, use that.  Otherwise,
  // it will be set by a flood fill
  // Added by Alex Main (September 2013)
  //
  if (!oldStatus)
    findActiveNodesUsingFloodFill(tId,points);
  else
    *status = *oldStatus;

  if(point_based_id) {
    *point_based_id = -999;
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++)
      for(int i=0; i<(*status)(iSub).size(); i++) {
        int myPID = (*status)(iSub)[i];
        map<int,int>::iterator pit = pid2id.find(myPID);
        if(pit!=pid2id.end()) {
          (*point_based_id)(iSub)[i] = myPID;
          (*status)(iSub)[i] = pit->second;
        }
      }
  }

  *distance=0.0;
  *status0=*status;
  *occluded_node0=*is_occluded;
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub)
	  for(int i = 0; i < (*is_active)(iSub).size(); ++i){
      (*is_active)(iSub)[i] = intersector[iSub]->testIsActive(0.0, i);
}


  ///////////////////////////////////////////////////////
  if(externalSI)
  {
	  DistVec<int> intISactive(domain->getNodeDistInfo());

#pragma omp parallel for
	  for(int iSub = 0; iSub < numLocSub; ++iSub)
	  {
		  intersector[iSub]->reFlagRealNodes(X(iSub));

		  for(int i=0; i<X(iSub).size(); ++i)
		  {
			  intISactive(iSub)[i] = (*is_active)(iSub)[i] ? 1 : 0;
			  (*TriID_node)(iSub)[i] = -1;
		  }
	  }

	  operMin<int> minOp;
	  domain->assemble(domain->getLevelPat(), intISactive, minOp);

#pragma omp parallel for
	  for(int iSub = 0; iSub < numLocSub; ++iSub)
		   for(int i=0; i<X(iSub).size(); ++i)
				(*is_active)(iSub)[i] = (intISactive(iSub)[i] == 1) ? true : false;

#pragma omp parallel for
	  for(int iSub = 0; iSub < numLocSub; ++iSub)
		  intersector[iSub]->ComputeSIbasedIntersections(iSub, X(iSub), (*boxMin)(iSub), (*boxMax)(iSub), withViscousTerms);
	  
  }
  ///////////////////////////////////////////////////////

}

//----------------------------------------------------------------------------

void DistIntersectorPhysBAM::findActiveNodesUsingFloodFill(const DistVec<bool>& tId, 
																			  const list<pair<Vec3D,int> >& points) {

  DistVec<int> nodeColors(X->info());

// Perform local floodFill
  int localColorCount[numLocSub];
#pragma omp parallel for
  for(int iSub=0;iSub<numLocSub;++iSub) {
		localColorCount[iSub] = FloodFill::floodFillSubDomain(*domain->getSubDomain()[iSub],
										 									   intersector[iSub]->edges,
                                                        intersector[iSub]->edge_intersects,
																				intersector[iSub]->is_occluded, 
																				nodeColors(iSub));
	}

  floodFill->generateConnectionsSet(*domain,*com,nodeColors);

   // only contains valid data for local SubDomains.
	map<pair<GLOBAL_SUBD_ID,int>,int> localToGlobalColorMap; 
  floodFill->unionColors(*domain,*com,numLocSub,localColorCount,localToGlobalColorMap);

// Determine the status of local colors
  map<int,int> globalColorToGlobalStatus;

#if 0 // Debug output
  for(list<pair<Vec3D,int> >::const_iterator iter = points.begin(); iter!=points.end(); iter++)
    com->fprintf(stderr,"found Point %d (%e %e %e) with specified fluid model and initial conditions.\n", 
                 iter->second, (iter->first)[0], (iter->first)[1], (iter->first)[2]);
#endif

  for(int iSub=0;iSub<numLocSub;++iSub){
    int ffNode=domain->getSubDomain()[iSub]->findFarfieldNode();
    if(ffNode >= 0){
#if 0 // Debug output
      fprintf(stderr,"%02d: Setting global color %d [from %d, %d] to Fluid Model %d BASED ON FF NODE\n",com->cpuNum(),
					 localToGlobalColorMap[pair<GLOBAL_SUBD_ID,int>(PhysBAM::getGlobSubNum(*domain->getSubDomain()[iSub]),
																					nodeColors(iSub)[ffNode])],
					 PhysBAM::getGlobSubNum(*domain->getSubDomain()[iSub]).Value(),
					 nodeColors(iSub)[ffNode],IntersectorPhysBAM::INSIDE);
#endif
      globalColorToGlobalStatus[localToGlobalColorMap[pair<GLOBAL_SUBD_ID,int>(PhysBAM::getGlobSubNum(*domain->getSubDomain()[iSub]),nodeColors(iSub)[ffNode])]]=IntersectorPhysBAM::INSIDE;
    }
  }

  for(int iSub=0;iSub<numLocSub;++iSub){
    SubDomain& sub=intersector[iSub]->subD;
    for(int iElem=0; iElem<sub.numElems(); iElem++)
      for(list<pair<Vec3D,int> >::const_iterator iP=points.begin(); iP!=points.end(); iP++){

			  if(sub.isINodeinITet(iP->first, iElem, (*X)(iSub))){ 
			  // TODO(jontg): Use a robust implementation of this routine

#if 0 // Debug output
          fprintf(stderr,"%02d: Setting global color %d [from %d, %d] to Fluid Model %d BASED ON POINT DATA\n",com->cpuNum(),
			  localToGlobalColorMap[pair<GLOBAL_SUBD_ID,int>(PhysBAM::getGlobSubNum(sub),nodeColors(iSub)[sub.getElemNodeNum(iElem)[0]])],
			  PhysBAM::getGlobSubNum(sub).Value(),nodeColors(iSub)[sub.getElemNodeNum(iElem)[0]], iP->second);
#endif
          globalColorToGlobalStatus[localToGlobalColorMap[pair<GLOBAL_SUBD_ID,int>(PhysBAM::getGlobSubNum(sub),nodeColors(iSub)[sub.getElemNodeNum(iElem)[0]])]]=iP->second;
        }
      }
  }

  PHYSBAM_MPI_UTILITIES::syncMap(*domain,*com,globalColorToGlobalStatus);

// Compute node status (occluded nodes are OUTSIDE the fluid regime)
#pragma omp parallel for
  for(int iSub=0;iSub<numLocSub;++iSub){
    SubDomain& sub=intersector[iSub]->subD;
    for(int i=0;i<(*status)(iSub).size();++i){
      int color=localToGlobalColorMap[pair<GLOBAL_SUBD_ID,int>(PhysBAM::getGlobSubNum(sub),nodeColors(iSub)[i])];
      if((*is_occluded)(iSub)[i] || globalColorToGlobalStatus.find(color)==globalColorToGlobalStatus.end()){
#if 0 // Debug output
        fprintf(stderr,"Flagging node %d as OUTSIDE COLOR %d, based on occluded = %d, global_color found = %d\n",
			intersector[iSub]->locToGlobNodeMap[i]+1, IntersectorPhysBAM::OUTSIDECOLOR, (*is_occluded)(iSub)[i],
			globalColorToGlobalStatus.find(color)!=globalColorToGlobalStatus.end());
        fprintf(stderr,"\t\tGLOBAL COLOR = %d, local color given as (%d, %d)\n", color, PhysBAM::getGlobSubNum(sub).Value(), nodeColors(iSub)[i]);
#endif
        (*is_occluded)(iSub)[i]=true;
        (*status)(iSub)[i]=IntersectorPhysBAM::OUTSIDECOLOR;

		} else { 
			(*status)(iSub)[i]=globalColorToGlobalStatus[color]; 
    }
  }
}
}

//----------------------------------------------------------------------------

void
DistIntersectorPhysBAM::findActiveNodes(const DistVec<bool>& tId) 
{

	int OUT       = IntersectorPhysBAM::OUTSIDECOLOR;
	int UNDECIDED = IntersectorPhysBAM::UNDECIDED;

    // Easy stuff first
#pragma omp parallel for
	for(int iSub=0; iSub<numLocSub; ++iSub)
	{
		for(int i=0; i<(*status)(iSub).size(); ++i)
		{
			if((*is_occluded)(iSub)[i]) 
				(*status)(iSub)[i] = OUT;
			else if(!(*is_swept)(iSub)[i] && !(*occluded_node0)(iSub)[i]) 
				(*status)(iSub)[i] = (*status0)(iSub)[i];
		}
	}

    operMax<int> maxOp;
    domain->assemble(domain->getLevelPat(),*status,maxOp);

    // Next handle swept nodes
    int iteration_count=0;
    int flags[2] = {1,1}; // flags = {needs_iteration, detected_change}

	while(flags[0]>0 && flags[1]>0)
	{
		flags[0] = 0; 
		flags[1] = 0; 

		++iteration_count;

#pragma omp parallel for
		for(int iSub=0; iSub<numLocSub; ++iSub)
		{
            SubDomain& sub=intersector[iSub]->subD;
            Connectivity &nToN = *(sub.getNodeToNode());

			for(int i=0; i<(*status)(iSub).size(); ++i) 
			{
				if((*status)(iSub)[i] == UNDECIDED)
				{
                int stat=UNDECIDED;

					for(int n=0; n<nToN.num(i); ++n)
					{
						int j = nToN[i][n];

                    if(i==j) continue;
                    if(!(*edge_intersects)(iSub)[intersector[iSub]->edges.findOnly(i,j)]
							&& (*status)(iSub)[j] != UNDECIDED && (*status)(iSub)[j] != OUT)
						{
							stat = (*status)(iSub)[j];
							flags[1] = 1;
							break;
            }
        }

					if(stat == UNDECIDED) 
						flags[0] = 1;
					else
					{
						(*status)(iSub)[i] = stat;
						(*is_swept)(iSub)[i] = true;
					}
            }
			}
		}

        com->globalOp(2,(int*)flags,MPI_SUM);

		if(flags[1]>0)
		{
            domain->assemble(domain->getLevelPat(),*status,maxOp);

#pragma omp parallel for
			for(int iSub=0; iSub<numLocSub; ++iSub) 
				for(int i=0; i<(*is_swept)(iSub).size(); ++i) 
					(*is_swept_helper)(iSub)[i] = (*is_swept)(iSub)[i] ? 1 : 0;

            domain->assemble(domain->getLevelPat(),*is_swept_helper,maxOp);
#pragma omp parallel for
			for(int iSub=0; iSub<numLocSub; ++iSub) 
				for(int i=0; i<(*is_swept)(iSub).size(); ++i) 
					(*is_swept)(iSub)[i] = (*is_swept_helper)(iSub)[i] >= 1 ? true : false;
        }
    }

    // Finish any remaining untouched nodes
#pragma omp parallel for
	for(int iSub=0; iSub<numLocSub; ++iSub)
	{
		SubDomain& sub = intersector[iSub]->subD;

        for(int i=0;i<(*status)(iSub).size();++i)
		{
			if((*status)(iSub)[i] == UNDECIDED)
			{
                (*status)(iSub)[i]=OUT;
				(*is_occluded)(iSub)[i] = true;
    }
}
	}

}

//----------------------------------------------------------------------------

void
DistIntersectorPhysBAM::updateStructure(double *xs, double *Vs, int nNodes, int (*abc)[3])
{
  int previous = numStNodes;
  if(cracking)
    gotNewCracking = cracking->getNewCrackingFlag();
  if(gotNewCracking) {
    cracking->setNewCrackingFlag(false);
  }
  if((nNodes!=numStNodes) && !gotNewCracking){
    fprintf(stderr,"ERROR: AERO-F is not sure if there is a new cracking... (Could be a software bug.) nNodes = %d, numStNodes = %d, gotNewCracking = %d\n", nNodes, numStNodes, gotNewCracking);
    exit(-1);
  }

  if(gotNewCracking) {
    if(!cracking || nNodes<numStNodes) {
      com->fprintf(stderr,"ERROR: Number of structure nodes has changed! Cracking is not considered in this simulation.\n");
      exit(-1);}
    // need to get the topology change caused by cracking.
    updateCracking(abc);
    if(nNodes!=numStNodes) {com->fprintf(stderr,"ERROR: horrible! %d v.s. %d\n", nNodes, numStNodes);exit(-1);}
  }

  for (int i=0; i<nNodes; i++)
    for(int j=0; j<3; j++) {
      Xs_n[i][j] = (i<previous) ? Xs_np1[i][j] : xs[3*i+j];
      Xs_np1[i][j] = xs[3*i+j];
      Xs[i][j] = Xs_np1[i][j];
      Xsdot[i][j] = Vs[3*i+j];}

  if(gotNewCracking) {
    std::map<int,int> newNodes = cracking->getLatestPhantomNodes();
    if(numStNodes-previous!=newNodes.size()) {
      com->fprintf(stderr,"SOFTWARE BUG: How many new phantom nodes, %d or %d ?!\n", numStNodes-previous, newNodes.size());
      exit(-1);
    }
    for(std::map<int,int>::iterator it=newNodes.begin(); it!=newNodes.end(); it++)
      Xs[it->first] = Xs[it->second]; //add new phantom nodes but keep the old node corrdinates.
    initializePhysBAM(); //delete the old interface and create a new one. Use the modified Xs.
  }

// add surface velocity contribution
  updatebc();
}

//----------------------------------------------------------------------------

void DistIntersectorPhysBAM::updateCracking(int (*abc)[3])
{ //update the topology, but not the node coordinates
  numStNodes = cracking->usedNodes();
  numStElems = cracking->usedTrias(); 
  std::set<int> newQuads = cracking->getLatestPhantomQuads();
  for(std::set<int>::iterator it=newQuads.begin(); it!=newQuads.end(); it++) {
    int trId1,trId2;
    cracking->getQuad2Tria(*it,trId1,trId2); //obtain trId1 and trId2;
    for(int j=0; j<3; j++) {
      stElem[trId1][j] = abc[trId1][j];
      if (trId2 >= 0)
	stElem[trId2][j] = abc[trId2][j];
    } 
  }
}

//----------------------------------------------------------------------------

void DistIntersectorPhysBAM::expandScope()
{
  SubDomain **subs = domain->getSubDomain();

  // 1. setup communication pattern
  int numLocSub = domain->getNumLocSub();
  SubDTopo *subTopo = domain->getSubTopo();
  CommPattern<int> trader(subTopo, com, CommPattern<int>::CopyOnSend, CommPattern<int>::NonSym);
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++) {
    int *sndChannel = subs[iSub]->getSndChannel();
    for(int iNei=0; iNei<subs[iSub]->getNumNeighb(); iNei++)
      trader.setLen(sndChannel[iNei], 1+intersector[iSub]->package[iNei].size());     
  }    
  trader.finalize();

  // 2. send packages to neighbour subdomains.
  set<int>::iterator it;
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++) {
    int *sndChannel = subs[iSub]->getSndChannel();
    for(int iNei=0; iNei<subs[iSub]->getNumNeighb(); iNei++) {
      SubRecInfo<int> sInfo = trader.getSendBuffer(sndChannel[iNei]);
      int *buffer = reinterpret_cast<int*>(sInfo.data);
      buffer[0] = intersector[iSub]->package[iNei].size(); 
      int count = 0;
      for(it=intersector[iSub]->package[iNei].begin(); it!= intersector[iSub]->package[iNei].end(); it++)
        buffer[++count] = *it;
    }
  }

  // 3. exchange
  trader.exchange();

  // 4. expand scope
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++) {

    int *rcvChannel = subs[iSub]->getRcvChannel();
    std::set<int>& scope=physInterface->getScope(iSub+1);
    for(int iNei=0; iNei<subs[iSub]->getNumNeighb(); iNei++) {
      SubRecInfo<int> sInfo = trader.recData(rcvChannel[iNei]);
      int *buffer = reinterpret_cast<int*>(sInfo.data);
      for(int j=0; j<buffer[0]; j++) scope.insert(buffer[j+1]);
    }
  }
}

//----------------------------------------------------------------------------

void
DistIntersectorPhysBAM::updatePhysBAMInterface(Vec3D *particles, int size, 
															  const DistSVec<double,3>& fluid_nodes,
															  const bool fill_scope, const bool retry) {

  if(!retry) physInterface->SaveOldState();
  for (int i=0; i<size; i++)
    physInterface->particles.X(i+1) = PhysBAM::VECTOR<double,3>(particles[i][0], particles[i][1], particles[i][2]);

  if(fill_scope){
    physInterface->Update(numLocSub,true);
    for(int iSub=0;iSub<numLocSub;++iSub){
      std::set<int>& scope=physInterface->getScope(iSub+1);
      if (!cracking)
	for(int i=1;i<=numStElems;++i) scope.insert(i);
      else
	for(int i=1;i<=cracking->numberRealTriangles();++i) scope.insert(i);	

    }

	} else {
		expandScope();
  }

	for(int iSub=0;iSub<numLocSub;++iSub) physInterface->UpdateScope(iSub+1);

}

//----------------------------------------------------------------------------

/** compute the intersections, node statuses and normals for the initial geometry */
int DistIntersectorPhysBAM::recompute(double dtf, double dtfLeft, double dts, bool findStatus, bool retry) 
{
	
	if(dtfLeft < -1.0e-6) 
	{
    fprintf(stderr,"There is a bug in time-step!\n");
    exit(-1);
  }
  //get current struct coordinates.
  double alpha = (dts - dtfLeft + dtf)/dts;
  for (int i=0; i<numStNodes; i++) 
    Xs[i] = (1.0-alpha)*Xs_n[i] + alpha*Xs_np1[i];

  // for hasCloseTriangle
  DistVec<bool> tId(X->info());
  
  updatePhysBAMInterface(Xs, numStNodes,*X, gotNewCracking,retry);

  buildSolidNormals();
  domain->findNodeBoundingBoxes(*X,*boxMin,*boxMax);

	DistVec<bool> bk_ISactive(domain->getNodeDistInfo());

	if(externalSI)
	{
#pragma omp parallel for
		for(int iSub = 0; iSub < numLocSub; ++iSub)
		{
			for(int i = 0; i < (*is_active)(iSub).size(); ++i)
				bk_ISactive(iSub)[i] = (*is_active)(iSub)[i];
		}
	}

#pragma omp parallel for
	for(int iSub = 0; iSub < numLocSub; ++iSub)
	{
      intersector[iSub]->reset(findStatus,retry);
      intersector[iSub]->hasCloseTriangle((*X)(iSub), (*Xn)(iSub), (*boxMin)(iSub), (*boxMax)(iSub), tId(iSub));
      intersector[iSub]->findIntersections((*X)(iSub),tId(iSub),*com);
      intersector[iSub]->computeSweptNodes((*X)(iSub), tId(iSub), *com, dtf);
	}

	if(findStatus) 
	{
    findActiveNodes(tId);
#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub)
		{
      for(int i = 0; i < (*is_active)(iSub).size(); ++i)
        (*is_active)(iSub)[i] = intersector[iSub]->testIsActive(0.0, i);
  }

		///////////////////////////////////////////////////////
		if(externalSI)
		{
			DistVec<int> intISactive(domain->getNodeDistInfo());
			DistVec<int>    intSwept(domain->getNodeDistInfo());

#pragma omp parallel for
			for(int iSub = 0; iSub < numLocSub; ++iSub)
			{
				intersector[iSub]->reFlagRealNodes((*X)(iSub), &bk_ISactive(iSub));

				for(int i=0; i<(*X)(iSub).size(); ++i)
				{
					intISactive(iSub)[i] = (*is_active)(iSub)[i] ? 1 : 0;
					   intSwept(iSub)[i] = (*is_swept)(iSub)[i]  ? 1 : 0;
					(*TriID_node)(iSub)[i] = -1;
				}
			}

			operMin<int> minOp;
			domain->assemble(domain->getLevelPat(), intISactive, minOp);

			operMax<int> maxOp;
			domain->assemble(domain->getLevelPat(), intSwept, maxOp);

#pragma omp parallel for
			for(int iSub = 0; iSub < numLocSub; ++iSub)
				for(int i=0; i<(*X)(iSub).size(); ++i) 
				{
					(*is_active)(iSub)[i] = (intISactive(iSub)[i] == 1) ? true : false;
					(*is_swept)(iSub)[i]  = (intSwept(iSub)[i]    == 1) ? true : false;
				}

#pragma omp parallel for
			for(int iSub = 0; iSub < numLocSub; ++iSub)
				intersector[iSub]->ComputeSIbasedIntersections(iSub, (*X)(iSub), (*boxMin)(iSub), (*boxMax)(iSub), withViscousTerms);	  
		}
		///////////////////////////////////////////////////////
	}

  return 0;
}

//----------------------------------------------------------------------------
void IntersectorPhysBAM::reFlagRealNodes(SVec<double,3>& X, Vec<bool> *bk_isActive)
{

	int (*ptr)[2] = edges.getPtr();

	bool recompute = false;
	if(bk_isActive) recompute = true;

	for(int l=0; l<edges.size(); l++) 
	{
		if(edge_intersects[l])
		{
			int i = ptr[l][0];
			int j = ptr[l][1];
		
			/*
			bool isOK = false;
			if( 
				((fabs(X[i][0]-0.165724)<1.0e-4 && fabs(X[i][1]-0.074293)<1.0e-4 && fabs(X[i][2]-0.03)<1.0e-4) &&
				 (fabs(X[j][0]-0.167208)<1.0e-4 && fabs(X[j][1]-0.034521)<1.0e-4 && fabs(X[j][2]-0.03)<1.0e-4)) ||
				((fabs(X[j][0]-0.165724)<1.0e-4 && fabs(X[j][1]-0.074293)<1.0e-4 && fabs(X[j][2]-0.03)<1.0e-4) &&
				 (fabs(X[i][0]-0.167208)<1.0e-4 && fabs(X[i][1]-0.034521)<1.0e-4 && fabs(X[i][2]-0.03)<1.0e-4)) 
				) isOK = true;
				
			if(isOK && recompute)
			{
				std::cout << X[i][0] << " " << X[i][1] << " " << X[i][2] << " | "  << X[j][0] << " " << X[j][1] << " " << X[j][2] << " || "
					       << std::boolalpha << is_active[i] << " " << is_active[j] << " ~~~ " << (*bk_isActive)[i] << " " << (*bk_isActive)[j] << " ++++ ";
			}
			else if(isOK && !recompute)
			{
				std::cout << X[i][0] << " " << X[i][1] << " " << X[i][2] << " | "  << X[j][0] << " " << X[j][1] << " " << X[j][2] << " || "
					       << std::boolalpha << is_active[i] << " " << is_active[j] << " ~~~ ";
			}
			*/

			// Intersection I ---> J
			const IntersectionResult<double>& Res_ij = CrossingEdgeRes[l];

			double alpha_1 = Res_ij.alpha;
			if(alpha_1 > 0.5 && is_active[i]) is_active[i] = false; 

			if(recompute)
			{
				if(is_active[i] != (*bk_isActive)[i]) is_swept[i] = true;
			}

			// Intersection J ---> I
			const IntersectionResult<double>& Res_ji = ReverseCrossingEdgeRes[l];
			
			double alpha_2 = Res_ji.alpha;
			if(alpha_2 > 0.5 && is_active[j]) is_active[j] = false; 

			if(recompute)
			{
				if( is_active[j] != (*bk_isActive)[j] ) is_swept[j] = true;
			}

			//if(isOK) std::cout << " " << std::boolalpha << is_active[i] << " " << is_active[j] << " " << is_swept[i] << " " << is_swept[j] << " \n";

		}

	}

}

//----------------------------------------------------------------------------
void IntersectorPhysBAM::ComputeSIbasedIntersections(int iSub, SVec<double,3>& X,  
																	  SVec<double,3> &boxMin, SVec<double,3> &boxMax, 
																	  bool withViscousTerms)
{
	int shrunk_index;

	bool occluded;

	double l1, l2, l1_i, l1_j, l2_i, l2_j, lambda[3];

	int Tri_, Tri_i, Tri_j;

	Vec3D *structX = (distIntersector.getStructPosition()).data();
	int (*triNodes)[3] = distIntersector.stElem;

	int (*ptr)[2] = edges.getPtr();

	for(int l=0; l<edges.size(); l++)
	{
		int i = ptr[l][0];
		int j = ptr[l][1];

		edge_SI[l] = false;
	  
		if(!is_active[i] || !is_active[j])
		{			
			if(is_active[i] == is_active[j]) continue;

			Vec3D X_ij;
			for(int k=0; k<3; ++k) X_ij[k] = 0.5*(X[j][k] + X[i][k]);

			edge_SI[l] = true;

			PhysBAMInterface<double>& physbam_interface=*distIntersector.physInterface;

			ARRAY<int> candidates_i;
			VECTOR<double,3> min_corner_i(boxMin[i][0],boxMin[i][1],boxMin[i][2]);
			VECTOR<double,3> max_corner_i(boxMax[i][0],boxMax[i][1],boxMax[i][2]);		
			bool t_i = physbam_interface.HasCloseTriangle(locIndex+1, VECTOR<double,3>(X[i][0],X[i][1],X[i][2]), 
																		 min_corner_i, max_corner_i, &shrunk_index, &occluded, &candidates_i);
		
			ARRAY<int> candidates_j;
			VECTOR<double,3>min_corner_j(boxMin[j][0],boxMin[j][1],boxMin[j][2]);
			VECTOR<double,3>max_corner_j(boxMax[j][0],boxMax[j][1],boxMax[j][2]);
			bool t_j= physbam_interface.HasCloseTriangle(locIndex+1, VECTOR<double,3>(X[i][0],X[i][1],X[i][2]), 
																		min_corner_j, max_corner_j, &shrunk_index, &occluded, &candidates_j);

			int trId;		
			double dist_si, dist_i, dist_j, xi[3];
			
			double min_dist_si = 1.0e16;
			double min_dist_i  = 1.0e16;
			double min_dist_j  = 1.0e16;

			Tri_ =-1;
			l1 = l2 = -1.0;

			Tri_i = Tri_j =-1;
			l1_i = l1_j = l2_i = l2_j = -1.0;
			
			if(t_i) 
			{
				for(int iArray=1; iArray <= candidates_i.Size(); iArray++) 
				{					
					trId = candidates_i(iArray) - 1;
			
					dist_si = piercing(X_ij, trId, lambda);

					if(dist_si < min_dist_si) 
					{
						min_dist_si = dist_si;

						Tri_ = trId;
						l1   = lambda[0];
						l2   = lambda[1];
					}

					if(withViscousTerms && !is_active[i])
					{
						dist_i = piercing(X[i], trId, lambda);

						if(dist_i < min_dist_i)
						{
							min_dist_i = dist_i;

							Tri_i = trId;
							l1_i  = lambda[0];
							l2_i  = lambda[1];
						}
					}
					
				}
			}

			if(t_j) 
			{		  
				for(int iArray=1; iArray <= candidates_j.Size(); iArray++) 
				{
					trId = candidates_j(iArray) - 1;

					dist_si = piercing(X_ij, trId, lambda);

					if(dist_si < min_dist_si) 
					{
						min_dist_si = dist_si;

						Tri_ = trId;
						l1   = lambda[0];
						l2   = lambda[1];
					}

					if(withViscousTerms && !is_active[j])
					{
						dist_j = piercing(X[j], trId, lambda);

						if(dist_j < min_dist_j)
						{
							min_dist_j = dist_j;

							Tri_j = trId;
							l1_j  = lambda[0];
							l2_j  = lambda[1];
						}
					}

				}
			}
			
			Vec3D d2, d3, nn;

			for(int id=0; id<3; ++id)
			{
				d2[id] = structX[triNodes[Tri_][1]][id] - structX[triNodes[Tri_][0]][id];
				d3[id] = structX[triNodes[Tri_][2]][id] - structX[triNodes[Tri_][0]][id];
			}

			nn = d2^d3;

			double norm = sqrt(nn * nn);
			if(norm != 0.0) nn *= (1.0 / norm);
			
			   xi_SI[l] = l1;
		     eta_SI[l] = l2;
			nWall_SI[l] = nn;
			TriID_SI[l] = Tri_;

			if(withViscousTerms && !is_active[i])
			{
				for(int id=0; id<3; ++id)
				{
					d2[id] = structX[triNodes[Tri_i][1]][id] - structX[triNodes[Tri_i][0]][id];
					d3[id] = structX[triNodes[Tri_i][2]][id] - structX[triNodes[Tri_i][0]][id];
				}

				nn = d2^d3;
				
				double norm = sqrt(nn * nn);
				if(norm != 0.0) nn *= 1.0 / norm;

				   xi_node[i] = l1_i;
 				  eta_node[i] = l2_i;
 				nWall_node[i] = nn;
				TriID_node[i] = Tri_i;
			}

			if(withViscousTerms && !is_active[j])
			{
				for(int id=0; id<3; ++id)
				{
					d2[id] = structX[triNodes[Tri_j][1]][id] - structX[triNodes[Tri_j][0]][id];
					d3[id] = structX[triNodes[Tri_j][2]][id] - structX[triNodes[Tri_j][0]][id];
				}

				nn = d2^d3;
				
				double norm = sqrt(nn * nn);
				if(norm != 0.0) nn *= 1.0 / norm;

				   xi_node[j] = l1_j;
 				  eta_node[j] = l2_j;
 				nWall_node[j] = nn;
				TriID_node[j] = Tri_j;
			}
		
		}

	}

/* --------- debug --------- */
#if 0
	for (int l=0; l<edges.size(); l++)
	{
		int i = ptr[l][0]; int j = ptr[l][1];

		if(edge_SI[l])
		{			
			if( (!is_active[i] && !is_active[j]) ||  
				 ( is_active[i] &&  is_active[j]) )
			{		
				std::cout << "error SI edge and " << std::boolalpha << is_active[i] << " " << is_active[j] << " " << edge_intersects[l] << endl;
				std::cout << i << " Xi = " << X[i][0]<<" "<<X[i][1]<<" "<<X[i][2] 
							 << " " << j << " Xj = " << X[j][0]<<" "<<X[j][1]<<" "<<X[j][2] <<endl;
				exit(-1);
			}
		}
		else
		{
			if( is_active[i] != is_active[j] )
			{	
				std::cout << "error  edge without SI and " << std::boolalpha << is_active[i] << " " << is_active[j] 
							 << " " << edge_intersects[l] << endl;
				std::cout << i << " Xi = " << X[i][0]<<" "<<X[i][1]<<" "<<X[i][2] << " " 
							 << j << " Xj = " << X[j][0]<<" "<<X[j][1]<<" "<<X[j][2] <<endl;
				exit(-1);
			}
		}

	}
#endif
}

//----------------------------------------------------------------------------

void IntersectorPhysBAM::printFirstLayer(SubDomain& sub, SVec<double,3>&X, int TYPE) {

  int mySub = sub.getGlobSubNum();
  int myLocSub = sub.getLocSubNum();
  int (*ptr)[2] = edges.getPtr();
  Connectivity &nToN = *(sub.getNodeToNode());
  std::string fileName = PhysBAM::STRING_UTILITIES::string_sprintf("firstLayer/%d/firstLayer%d.top",mySub);
  std::string nodesName = PhysBAM::STRING_UTILITIES::string_sprintf("InsideNodes%d\n",mySub);

  FILE* firstLayer = fopen(fileName.c_str(),"w");
  fprintf(firstLayer, "Nodes %s\n", nodesName.c_str());
  for (int i=0; i<sub.numNodes(); i++) fprintf(firstLayer,"%d %e %e %e\n", i+1, X[i][0], X[i][1], X[i][2]);
  fprintf(firstLayer, "Elements FirstLayer%d using %s\n", mySub, nodesName.c_str());
  for (int l=0; l<edges.size(); l++){
    int x1 = ptr[l][0], x2 = ptr[l][1];
    if(edge_intersects[l]) fprintf(firstLayer,"%d %d %d %d\n", l+1, (int)1, x1+1, x2+1);}
  fclose(firstLayer);
}

//----------------------------------------------------------------------------

IntersectorPhysBAM::IntersectorPhysBAM(SubDomain &sub,SVec<double,3> &X, 
													Vec<int> &stat0, Vec<ClosestPoint> &clo,
				       Vec<bool>& occ_node0, DistIntersectorPhysBAM &distInt) :
  LevelSetStructure((*distInt.status)(sub.getLocSubNum()),(*distInt.distance)(sub.getLocSubNum()),
		    (*distInt.is_swept)(sub.getLocSubNum()), (*distInt.is_active)(sub.getLocSubNum()),
							(*distInt.is_occluded)(sub.getLocSubNum()), (*distInt.edge_intersects)(sub.getLocSubNum()),
							(*distInt.edge_SI)(sub.getLocSubNum()),
							(*distInt.xi_SI)(sub.getLocSubNum()),(*distInt.eta_SI)(sub.getLocSubNum()),
							(*distInt.nWall_SI)(sub.getLocSubNum()),
							(*distInt.TriID_SI)(sub.getLocSubNum()),
							(*distInt.xi_node)(sub.getLocSubNum()),(*distInt.eta_node)(sub.getLocSubNum()),
		               (*distInt.nWall_node)(sub.getLocSubNum()),
							(*distInt.TriID_node)(sub.getLocSubNum()) ),
  subD(sub),edges(sub.getEdges()),distIntersector(distInt),package(0),
  nodeToSubD(*sub.getNodeToSubD()), status0(stat0), closest(clo),
  occluded_node0(occ_node0),locIndex(sub.getLocSubNum()),globIndex(sub.getGlobSubNum())
{
  int numEdges = edges.size();

  status0 = UNDECIDED;
  occluded_node0 = false;
  OUTSIDECOLOR = distInt.numOfFluids();

  locToGlobNodeMap = sub.getNodeMap();

  package = new set<int>[sub.getNumNeighb()];
  int *neighb = sub.getNeighb();

  for(int i=0; i<sub.getNumNeighb(); i++)  sub2pack[neighb[i]] = i;

  nFirstLayer = 0;
  status = UNDECIDED;
  distance = -1.0;
  is_swept = false;
  is_active = false;
  is_occluded = false;
  edge_intersects = false;
  edge_SI = false;
}


//----------------------------------------------------------------------------

IntersectorPhysBAM::~IntersectorPhysBAM()
{
    delete []package;
}

//----------------------------------------------------------------------------

void IntersectorPhysBAM::reset(const bool findStatus,const bool retry)
{
  for(int i=0; i<subD.getNumNeighb(); i++) package[i].clear();

  if(!retry){
    status0 = status;
    occluded_node0 = is_occluded;
  }
  status = UNDECIDED;
  distance = -1.0;
  is_swept = false;
  is_active = false;
  is_occluded = false;
  nFirstLayer = 0;
  edge_intersects = false;
  edge_SI = false;

  TriID_SI = -1;
  TriID_node = -1;

  reverse_mapping.Remove_All();
  forward_mapping.Remove_All();
  xyz.Remove_All();
  xyz_n.Remove_All();
}

//----------------------------------------------------------------------------

/*
  Find the closest structural triangle for each node. 
  If no triangle intersect the bounding box of the node,
  no closest triangle exists
*/
int IntersectorPhysBAM::hasCloseTriangle(SVec<double,3> &X, SVec<double,3> &Xn, 
													  SVec<double,3> &boxMin, SVec<double,3> &boxMax, 
													  Vec<bool> &tId) 
{
  PhysBAMInterface<double>& physbam_interface=*distIntersector.physInterface;

  int numCloseNodes=0;
  forward_mapping.Resize(X.size()+1);  PhysBAM::ARRAYS_COMPUTATIONS::Fill(forward_mapping,-1);
  reverse_mapping.Resize(X.size()+1);  PhysBAM::ARRAYS_COMPUTATIONS::Fill(reverse_mapping,-1);

  xyz.Resize(X.size()+1);  xyz_n.Resize(Xn.size()+1);

  for(int i=0;i<boxMin.size();++i){
    ARRAY<int> candidates;
    VECTOR<double,3> min_corner(boxMin[i][0],boxMin[i][1],boxMin[i][2]), max_corner(boxMax[i][0],boxMax[i][1],boxMax[i][2]);

    int shrunk_index; 
	 bool occluded;

    tId[i] = physbam_interface.HasCloseTriangle(locIndex+1, VECTOR<double,3>(X[i][0],X[i][1],X[i][2]), 
																min_corner, max_corner, &shrunk_index, &occluded, &candidates);
	 
    closest[i].mode = -2; //set to "unknown"
    if(tId[i]){
      forward_mapping(i+1)=shrunk_index;
      is_occluded[i]=occluded;
      xyz(forward_mapping(i+1))=VECTOR<double,3>(X[i][0],X[i][1],X[i][2]);
      xyz_n(forward_mapping(i+1))=VECTOR<double,3>(Xn[i][0],Xn[i][1],Xn[i][2]);
      reverse_mapping(forward_mapping(i+1))=i;
      ++numCloseNodes;
      for(int j=1;j<=candidates.Size();++j) addToPackage(i,candidates(j));
      Vec3D x0(X[i][0], X[i][1], X[i][2]);

      if(distIntersector.cracking || 1 /*need a flag for 'multi-phase'*/) {
        findNodeClosestPoint(i,x0,candidates); //fill closest[i]
    } 

    } else {

      is_occluded[i]=false;
      closest[i].mode = -1; // set to "far"
    }
  }

  nFirstLayer += numCloseNodes;

  // Trim the fat
  reverse_mapping.Resize(nFirstLayer);
  xyz.Resize(nFirstLayer); 
  xyz_n.Resize(nFirstLayer);

  return numCloseNodes;
}

//----------------------------------------------------------------------------

int IntersectorPhysBAM::findIntersections(SVec<double,3>&X,Vec<bool>& tId,Communicator& com)
{
  int (*ptr)[2] = edges.getPtr();
  ARRAY<bool> occludedNode(nFirstLayer); PhysBAM::ARRAYS_COMPUTATIONS::Fill(occludedNode,false);

	for(int i=1;i<=nFirstLayer;++i) {
		occludedNode(i) = is_occluded[reverse_mapping(i)];
	}

#if 0 // Debug output
  for(int i=0;i<tId.size();++i){
      if(tId[i] && (i != reverse_mapping(forward_mapping(i+1))))
			fprintf(stderr,"BAH %d, %d, %d\n",i+1,forward_mapping(i+1),reverse_mapping(forward_mapping(i+1)));
		else if(!tId[i] && forward_mapping(i+1) != -1)
			fprintf(stderr,"forward_mapping %d should be -1, is %d (tId = %d)\n",i,forward_mapping(i+1),tId[i]);
	}

	for(int i=1;i<=nFirstLayer;++i){
		if(i != forward_mapping(reverse_mapping(i)+1)) 
			fprintf(stderr,"BLAH %d, %d, %d\n",i,reverse_mapping(i)+1,forward_mapping(reverse_mapping(i)+1));
	}

	for(int i=0;i<tId.size();++i) {
		if(tId[i] && occludedNode(forward_mapping(i+1)) != is_occluded[i]) 
			fprintf(stderr,"OCC ERR, %d -> %d, %d -> %d\n",i,forward_mapping(i+1),is_occluded[i],occludedNode(forward_mapping(i+1)));
		else if(!tId[i] && is_occluded[i]) 
			fprintf(stderr,"OCC ERR!!!\n");
	}
#endif

  ARRAY<TRIPLE<VECTOR<int,3>,IntersectionResult<double>,IntersectionResult<double> > > edgeRes;
  for (int l=0; l<edges.size(); l++) {
    int p = ptr[l][0], q = ptr[l][1];
    if(tId[p] && tId[q])
			edgeRes.Append(TRIPLE<VECTOR<int,3>, 
								IntersectionResult<double>, 
								IntersectionResult<double> >(VECTOR<int,3>(forward_mapping(p+1),forward_mapping(q+1),l),
																	  IntersectionResult<double>(),IntersectionResult<double>()));
	}

  distIntersector.getInterface().Intersect(locIndex+1,xyz,occludedNode,edgeRes);

  int intersectedEdgeCount=0;
  for(int i=1; i<=edgeRes.Size(); ++i){
      if(edgeRes(i).y.triangleID > 0 && edgeRes(i).z.triangleID > 0) {
#if 0 // Debug output
	  int* tr0=distIntersector.stElem[edgeRes(i).y.triangleID-1];
	  int* tr1=distIntersector.stElem[edgeRes(i).z.triangleID-1];
          if(edgeRes(i).y.triangleID > distIntersector.getNumStructElems() || tr0[0] > distIntersector.getNumStructNodes() || 
	     tr0[1] > distIntersector.getNumStructNodes() || tr0[2] > distIntersector.getNumStructNodes() || 
             edgeRes(i).z.triangleID > distIntersector.getNumStructElems() || tr1[0] > distIntersector.getNumStructNodes() || 
	     tr1[1] > distIntersector.getNumStructNodes() || tr1[2] > distIntersector.getNumStructNodes()){
            fprintf(stderr,"Detected a WEIRD intersection case: [%d, %e] (%e,%e,%e) (%d,%d,%d) AND [%d, %e] (%e,%e,%e) (%d,%d,%d)\n",
			    edgeRes(i).y.triangleID,edgeRes(i).y.alpha,edgeRes(i).y.zeta[0],edgeRes(i).y.zeta[1],edgeRes(i).y.zeta[2],tr0[0],tr0[1],tr0[2],
			    edgeRes(i).z.triangleID,edgeRes(i).z.alpha,edgeRes(i).z.zeta[0],edgeRes(i).z.zeta[1],edgeRes(i).z.zeta[2],tr1[0],tr1[1],tr1[2]);
	  }
#endif
          int l=edgeRes(i).x.z;
          edge_intersects[l] = true;
          CrossingEdgeRes[l] = edgeRes(i).y;
          ReverseCrossingEdgeRes[l] = edgeRes(i).z;

	  if (distIntersector.cracking) {

	    CrossingEdgeRes[l].triangleID  = 
	      distIntersector.cracking->mapTriangleID(CrossingEdgeRes[l].triangleID-1) + 1;
	    
	    ReverseCrossingEdgeRes[l].triangleID = 
	      distIntersector.cracking->mapTriangleID(ReverseCrossingEdgeRes[l].triangleID-1) + 1;

	  }

          ++intersectedEdgeCount;

		}
	}

#if 0 // Debug output
  bool should_quit=false;
  for(int l=0;l<edges.size();++l){
      int p = ptr[l][0], q = ptr[l][1];
      if((!tId[p] || !tId[q]) && (edge_intersects[l])){
			fprintf(stderr, "%03d bad edge, case 1, on edge number %d (%d) between nodes %d (%d,%d,%d) and %d (%d,%d,%d).\n", 
					  globIndex, l, edge_intersects[l], p, locToGlobNodeMap[p]+1, tId[p], is_occluded[p], q, locToGlobNodeMap[q]+1, 
					  tId[q],is_occluded[q]);
			should_quit=true;
		}
      if((is_occluded[p] || is_occluded[q]) && (!edge_intersects[l])){
			fprintf(stderr, "%03d encountered a bad edge, case 2, on edge number %d (%d) between nodes %d (%d,%d,%d) and %d (%d,%d,%d).\n", 
					  globIndex,l,edge_intersects[l], p,locToGlobNodeMap[p]+1,tId[p],is_occluded[p],q,locToGlobNodeMap[q]+1, 
					  tId[q],is_occluded[q]);
			should_quit=true;
		}
      if(edge_intersects[l] && (CrossingEdgeRes[l].triangleID == -1 || ReverseCrossingEdgeRes[l].triangleID == -1)){
			fprintf(stderr, "%03d encountered a bad edge, case 3, on edge number %d (%d) between nodes %d (%d,%d,%d) and %d (%d,%d,%d).\n",
					  globIndex,l,edge_intersects[l], p,locToGlobNodeMap[p]+1,tId[p],is_occluded[p],q,locToGlobNodeMap[q]+1,
					  tId[q],is_occluded[q]);
			should_quit=true;
		}
  }

  for(int i=1; i<=edgeRes.Size(); ++i){
	  if(!edge_intersects[edgeRes(i).x.z] && (is_occluded[reverse_mapping(edgeRes(i).x.x)] || is_occluded[reverse_mapping(edgeRes(i).x.y)])) {
		  fprintf(stderr, 
					 "Detected a bad edge, case 2, with the following results: %d (%d,%d, [%d])\n\t%d -> %d (%d), %e, (%e, %e, %e)\n\t%d -> %d (%d), %e, (%e, %e, %e)\n", edgeRes(i).x.z, reverse_mapping(edgeRes(i).x.x),reverse_mapping(edgeRes(i).x.y),
					 edges.find(reverse_mapping(edgeRes(i).x.x),reverse_mapping(edgeRes(i).x.y)),
					 ptr[edgeRes(i).x.z][0],ptr[edgeRes(i).x.z][1],edgeRes(i).y.triangleID, 
					 edgeRes(i).y.alpha, edgeRes(i).y.zeta[0], edgeRes(i).y.zeta[1], edgeRes(i).y.zeta[2],
					 ptr[edgeRes(i).x.z][1],ptr[edgeRes(i).x.z][0],edgeRes(i).z.triangleID, 
					 edgeRes(i).z.alpha, edgeRes(i).z.zeta[0], edgeRes(i).z.zeta[1], edgeRes(i).z.zeta[2]);
  }

  }

  if(should_quit) exit(-1);
#endif
  return intersectedEdgeCount;
}

//----------------------------------------------------------------------------

int IntersectorPhysBAM::computeSweptNodes(SVec<double,3>& X, Vec<bool>& tId,Communicator& com,const double dt)
{
  ARRAY<bool> swept;

  swept.Resize(nFirstLayer); PhysBAM::ARRAYS_COMPUTATIONS::Fill(swept,false);
  distIntersector.getInterface().computeSweptNodes(locIndex+1,xyz,xyz_n,swept,(double)1.0);

  int numSweptNodes=0;

  for(int i=1;i<=nFirstLayer;++i)
  {
	  is_swept[reverse_mapping(i)] = swept(i);
	  if(swept(i)) ++numSweptNodes;
  }

  return numSweptNodes;
}

//----------------------------------------------------------------------------

LevelSetResult
IntersectorPhysBAM::getLevelSetDataAtEdgeCenter(double t, int l, bool i_less_j, double *Xr, double *Xg) {

  if (!edge_intersects[l]) {

    int (*ptr)[2] = edges.getPtr();
    int i=i_less_j ? ptr[l][0] : ptr[l][1],
        j=i_less_j ? ptr[l][1] : ptr[l][0];

    fprintf(stderr,"%02d There is no intersection between node %d(status:%d,occluded=%d) and %d(status:%d,occluded=%d) along edge %d! Abort...\n",
                   globIndex,locToGlobNodeMap[i]+1, status[i],is_occluded[i], locToGlobNodeMap[j]+1, status[j],is_occluded[j],l);

    PHYSBAM_MPI_UTILITIES::dump_stack_trace();

    exit(-1);
  }

  const IntersectionResult<double>& result = i_less_j ? CrossingEdgeRes[l] : ReverseCrossingEdgeRes[l];
  
  double alpha0      = result.alpha;
  int trueTriangleID = result.triangleID-1;

  LevelSetResult lsRes;
  lsRes.alpha      = alpha0;

  lsRes.xi[0]      = result.zeta[0];
  lsRes.xi[1]      = result.zeta[1];
  lsRes.xi[2]      = 1-result.zeta[0]-result.zeta[1];

  lsRes.trNodes[0] = distIntersector.stElem[trueTriangleID][0];
  lsRes.trNodes[1] = distIntersector.stElem[trueTriangleID][1];
  lsRes.trNodes[2] = distIntersector.stElem[trueTriangleID][2];

  lsRes.normVel    = lsRes.xi[0]*distIntersector.Xsdot[lsRes.trNodes[0]]
                   + lsRes.xi[1]*distIntersector.Xsdot[lsRes.trNodes[1]]
                   + lsRes.xi[2]*distIntersector.Xsdot[lsRes.trNodes[2]]; 

  lsRes.porosity   = distIntersector.porosity[trueTriangleID];
  lsRes.structureType   = distIntersector.structureType[trueTriangleID];
  lsRes.actuatorDiskMethod = distIntersector.actuatorDiskMethod[trueTriangleID];

  lsRes.actuatorDiskPressureJump = distIntersector.actuatorDiskPressureJump[trueTriangleID];
  lsRes.isCorrectedMethod = distIntersector.isCorrectedMethod[trueTriangleID];
  lsRes.gamma = distIntersector.gamma;
  lsRes.actuatorDiskReconstructionMethod = distIntersector.actuatorDiskReconstructionMethod[trueTriangleID];



  if(!distIntersector.interpolatedNormal){

    lsRes.gradPhi = distIntersector.getSurfaceNorm(trueTriangleID);

    if(distIntersector.with_sensitivity && Xr != 0 && Xg != 0){
      Vec3D *XX  = (distIntersector.getStructPosition()).data();
      
      Vec3D dNdS;  
      derivativeOFnormal(XX[lsRes.trNodes[0]], 
			 XX[lsRes.trNodes[1]], 
			 XX[lsRes.trNodes[2]],
			 distIntersector.dXdSb[lsRes.trNodes[0]],
			 distIntersector.dXdSb[lsRes.trNodes[1]],
			 distIntersector.dXdSb[lsRes.trNodes[2]], dNdS);
      lsRes.dnds = dNdS;

      double dads;
      dads = derivativeOFalpha(XX[lsRes.trNodes[0]], 
                               XX[lsRes.trNodes[1]], 
       			       XX[lsRes.trNodes[2]],
       			       distIntersector.dXdSb[lsRes.trNodes[0]],
			       distIntersector.dXdSb[lsRes.trNodes[1]],
			       distIntersector.dXdSb[lsRes.trNodes[2]],
			       Xr, Xg);
      lsRes.dads = dads;
    }

  } else { //use nodal normals.    
    Vec3D ns0 = distIntersector.getNodalNorm(lsRes.trNodes[0]);
    Vec3D ns1 = distIntersector.getNodalNorm(lsRes.trNodes[1]);
    Vec3D ns2 = distIntersector.getNodalNorm(lsRes.trNodes[2]);
    lsRes.gradPhi = lsRes.xi[0]*ns0 + lsRes.xi[1]*ns1 + lsRes.xi[2]*ns2;
    lsRes.gradPhi /= lsRes.gradPhi.norm();
  }


  return lsRes;
}

//----------------------------------------------------------------------------

double IntersectorPhysBAM::isPointOnSurface(Vec3D pt, int N1, int N2, int N3) 
{
  Vec<Vec3D> &solidX = distIntersector.getStructPosition();
  Vec3D X1 = solidX[N1];
  Vec3D X2 = solidX[N2];
  Vec3D X3 = solidX[N3];

  Vec3D normal = (X2-X1)^(X3-X1);
  normal /=  normal.norm();

  return fabs((pt-X1)*normal);
}

//----------------------------------------------------------------------------

void IntersectorPhysBAM::addToPackage(const int i, const int candidate)
{
  int nSub=nodeToSubD.num(i);
  for(int n=0;n<nSub;n++){
    if(nodeToSubD[i][n]==globIndex) continue;
    package[sub2pack[nodeToSubD[i][n]]].insert(candidate);}
}

//----------------------------------------------------------------------------

void IntersectorPhysBAM::findNodeClosestPoint(const int nodeId, Vec3D& x0, ARRAY<int> &cand)
{
  int (*triNodes)[3] = distIntersector.stElem;
  Vec3D *structX = (distIntersector.getStructPosition()).data();
  double dist, xi[3];
  int trId, n1, n2; 
  int mod = -2;
  const double eps = 0;

  for(int iArray=1; iArray<=cand.Size(); iArray++) {
	  
    // Speed improvement.  When we are doing cracking,
    // we do not add phantom triangles to the triangle hierarchy.
    // Added by Alex Main (October 2013)
    //
    if (!distIntersector.cracking)
      trId = cand(iArray)-1;
    else
      trId = distIntersector.cracking->mapTriangleID(cand(iArray)-1);
      

    dist = std::abs(project(x0, trId, xi[0], xi[1])); // project on plane
    xi[2] = 1.0-xi[0]-xi[1];
    if(xi[0] >= -eps && xi[1] >= -eps && xi[2] >= -eps) 
      mod = 0;
    else { 
      dist = 1.0e16;
      for (int i=0; i<3; i++) {
        if(xi[i]<-eps) {
          int p1 = triNodes[trId][(i+1)%3], p2 = triNodes[trId][(i+2)%3];
          double alpha;
          double d2l = std::abs(edgeProject(x0, p1, p2, alpha)); // project on line
          if(alpha>=-eps && alpha<=1.0+eps) {
            if(dist>d2l) {dist = d2l; mod = 1; n1 = p1; n2 = p2; xi[i] = 0.0; xi[(i+1)%3] = 1.0-alpha; xi[(i+2)%3] = alpha;}
          } else {
            if(alpha<-eps) {
              double d2p = (x0-structX[p1]).norm();
              if(dist>d2p) {dist = d2p; mod = 2; n1 = n2 = p1; xi[i] = xi[(i+2)%3] = 0.0; xi[(i+1)%3] = 1.0;}
            }
            if(alpha>1.0+eps) {
              double d2p = (x0-structX[p2]).norm();
              if(dist>d2p) {dist = d2p; mod = 2; n1 = n2 = p2; xi[i] = xi[(i+1)%3] = 0.0; xi[(i+2)%3] = 1.0;}
            }
          }
        }
      }
    }

    // if the point is in the "phantom" region...
    if(distIntersector.cracking && distIntersector.cracking->getPhi(trId, xi[0], xi[1])<0.0) {
      double phi[3];
      phi[0] = distIntersector.cracking->getPhi(trId, 1.0, 0.0);
      phi[1] = distIntersector.cracking->getPhi(trId, 0.0, 1.0);
      phi[2] = distIntersector.cracking->getPhi(trId, 0.0, 0.0);
      if(phi[0]<0.0 && phi[1]<0.0 && phi[2]<0.0)
        continue; //this is purely a phantom...

      Vec3D p[2];
      int edgeNo[2];
      int count = 0;
      for(int i=0; i<3; i++) 
        if(phi[i]*phi[(i+1)%3]<0.0) {
          p[count] = std::abs(phi[(i+1)%3])/(std::abs(phi[i])+std::abs(phi[(i+1)%3]))*structX[triNodes[trId][i]] +
                       std::abs(phi[i])/(std::abs(phi[i])+std::abs(phi[(i+1)%3]))*structX[triNodes[trId][(i+1)%3]];
          edgeNo[count++] = i;
        }

      if(count!=2) {fprintf(stderr,"Debug: this is not right! count = %d.\n", count); exit(-1);}
      double a0;
      double d2l = edgeProject(x0, p[0], p[1], a0); 
      if(a0>=-eps && a0<=1.0+eps) {dist = std::abs(d2l); mod = 0;}
      else if(a0<-eps) {a0 = 0.0; dist = (x0-p[0]).norm(); mod = 1; n1 = triNodes[trId][edgeNo[0]]; n2 = triNodes[trId][(edgeNo[0]+1)%3];}
      else {a0 = 1.0; dist = (x0-p[1]).norm(); mod = 1; n1 = triNodes[trId][edgeNo[1]]; n2 = triNodes[trId][(edgeNo[1]+1)%3];} 
      double dd = project((1.0-a0)*p[0]+a0*p[1], trId, xi[0], xi[1]);
      if(std::abs(dd)>1.0e-6) {fprintf(stderr,"Debug: interpolation error in IntersectorPhysBAM... dd = %e\n", dd);exit(-1);}
      xi[2] = 1.0-xi[0]-xi[1];
    }

    // Now we got for this triangle: trId, dist, mod, xi[0,1,2], n1, n2.
    if(closest[nodeId].mode==-2 || closest[nodeId].dist>dist) {
      closest[nodeId].mode = mod;
      closest[nodeId].dist = dist;
      distance[nodeId] = dist;
      closest[nodeId].xi1 = xi[0]; 
      closest[nodeId].xi2 = xi[1]; 
      if(mod==0) {closest[nodeId].tracker[0] = trId;}
      else if(mod==1) {closest[nodeId].tracker[0] = n1; closest[nodeId].tracker[1] = n2;}
      else if(mod==2) {closest[nodeId].tracker[0] = n1;}
      else {fprintf(stderr,"Debug: error in mode! mod = %d. Should be 0, 1, or 2.\n", mod); exit(-1);}
    } 
  }
}

//----------------------------------------------------------------------------

double IntersectorPhysBAM::edgeProject(Vec3D x0, Vec3D& xA, Vec3D& xB, double &alpha) const
{
  Vec3D AB= xB-xA;
  Vec3D AX = x0-xA;
  alpha = AB*AX/(AB*AB);
  Vec3D P = xA + alpha*AB;
  return (P-x0).norm();
}

//----------------------------------------------------------------------------

double IntersectorPhysBAM::edgeProject(Vec3D x0, int n1, int n2, double &alpha) const
{
  Vec3D *structX = (distIntersector.getStructPosition()).data();

  Vec3D xA = structX[n1];
  Vec3D xB = structX[n2];
  return edgeProject(x0, xA, xB, alpha);
}
//----------------------------------------------------------------------------

double IntersectorPhysBAM::project(Vec3D x0, int tria, double& xi1, double& xi2) const
{
  int (*triNodes)[3] = distIntersector.stElem;
  Vec3D *structX = (distIntersector.getStructPosition()).data();

  int iA = triNodes[tria][0];
  int iB = triNodes[tria][1];
  int iC = triNodes[tria][2];
   Vec3D xA = structX[iA];
  Vec3D xB = structX[iB];
  Vec3D xC = structX[iC];

  Vec3D ABC = (xB-xA)^(xC-xA);
  double areaABC = ABC.norm();
  Vec3D dir = 1.0/areaABC*ABC;

  //calculate the projection.
  double dist = (x0-xA)*dir;
  Vec3D xp = x0 - dist*dir;

  //calculate barycentric coords.
  double areaPBC = (((xB-xp)^(xC-xp))*dir);
  double areaPCA = (((xC-xp)^(xA-xp))*dir);
  xi1 = areaPBC/areaABC;
  xi2 = areaPCA/areaABC;

  return dist;
}

//----------------------------------------------------------------------------

double IntersectorPhysBAM::piercing(Vec3D x0, int tria, double xi[3])
{

   const double tol = 1.0e-12;

	int (*triNodes)[3] = distIntersector.stElem;
	Vec3D *structX = (distIntersector.getStructPosition()).data();

	double dist;

	dist = std::abs( project(x0, tria, xi[0], xi[1]) );

	xi[2] = 1.0 - xi[0] - xi[1];

	if(xi[0] >= tol && xi[1] >= tol && xi[2] >= tol)
	{		
		return dist;
	} 
	else
	{
		dist = 1.0e16;

		for(int i=0; i<3; ++i)
		{
			if(xi[i] < tol)
			{
				int p1 = triNodes[tria][(i+1)%3];
				int p2 = triNodes[tria][(i+2)%3];
				
				double alpha;
				double d2l = std::abs(edgeProject(x0, p1, p2, alpha));

				if(alpha >= tol && alpha <= 1.0 + tol) 
				{
					if(dist > d2l) 
					{
						dist        = d2l; 
						xi[i]       = 0.0; 
						xi[(i+1)%3] = 1.0-alpha; 
						xi[(i+2)%3] = alpha;
					}
				} 
				else 
				{
					if(alpha < tol) 
					{
						double d2p = (x0 - structX[p1]).norm();

						if(dist > d2p) 
						{
							dist        = d2p; 
							xi[i]       = 0.0;
							xi[(i+2)%3] = 0.0; 
							xi[(i+1)%3] = 1.0;
						}
					}

					if(alpha > 1.0+tol) 
					{
						double d2p = (x0 - structX[p2]).norm();

						if(dist > d2p) 
						{
							dist        = d2p; 
							xi[i]       = 0.0;
							xi[(i+1)%3] = 0.0; 
							xi[(i+2)%3] = 1.0;
						}
					}
				}
			}
		}

		return dist;
	}

}

//----------------------------------------------------------------------------

void IntersectorPhysBAM::xWallWithSI(int n, Vec3D &xWall)
{

	double l1 =  xi_SI[n];
	double l2 = eta_SI[n];

	int tria  = TriID_SI[n];

	Vec3D X0 = distIntersector.Xs[distIntersector.stElem[tria][0]];
	Vec3D X1 = distIntersector.Xs[distIntersector.stElem[tria][1]];
	Vec3D X2 = distIntersector.Xs[distIntersector.stElem[tria][2]];

	for(int k=0; k<3; ++k) 
		xWall[k] = X2[k] + l1*(X0[k] - X2[k]) + l2*(X1[k] - X2[k]);

}

//----------------------------------------------------------------------------

void IntersectorPhysBAM::vWallWithSI(int n, Vec3D &vWall)
{
  
	double l1 =  xi_SI[n];
	double l2 = eta_SI[n];

	int tria  = TriID_SI[n];

	Vec3D v0 = distIntersector.Xsdot[distIntersector.stElem[tria][0]];
	Vec3D v1 = distIntersector.Xsdot[distIntersector.stElem[tria][1]];
	Vec3D v2 = distIntersector.Xsdot[distIntersector.stElem[tria][2]];

	for(int k=0; k<3; ++k) 
		vWall[k] = v2[k] + l1*(v0[k] - v2[k]) + l2*(v1[k] - v2[k]);

}

//----------------------------------------------------------------------------

bool IntersectorPhysBAM::xWallNode(int i, Vec3D &xWall)
{

	double l1 =  xi_node[i];
	double l2 = eta_node[i];

	int tria  = TriID_node[i];

	if(tria < 0) return false;
	
	Vec3D X0 = distIntersector.Xs[distIntersector.stElem[tria][0]];
	Vec3D X1 = distIntersector.Xs[distIntersector.stElem[tria][1]];
	Vec3D X2 = distIntersector.Xs[distIntersector.stElem[tria][2]];

	for(int k=0; k<3; ++k) 
		xWall[k] = X2[k] + l1*(X0[k] - X2[k]) + l2*(X1[k] - X2[k]);

	return true;

}

//----------------------------------------------------------------------------

bool IntersectorPhysBAM::vWallNode(int i, Vec3D &vWall)
{
  
	double l1 =  xi_node[i];
	double l2 = eta_node[i];

	int tria  = TriID_node[i];

	if(tria < 0) return false;

	Vec3D v0 = distIntersector.Xsdot[distIntersector.stElem[tria][0]];
	Vec3D v1 = distIntersector.Xsdot[distIntersector.stElem[tria][1]];
	Vec3D v2 = distIntersector.Xsdot[distIntersector.stElem[tria][2]];

	for(int k=0; k<3; ++k) 
		vWall[k] = v2[k] + l1*(v0[k] - v2[k]) + l2*(v1[k] - v2[k]);

	return true;
}

//----------------------------------------------------------------------------

void IntersectorPhysBAM::derivativeOFnormal(Vec3D  xA, Vec3D  xB, Vec3D  xC, 
					    Vec3D dxA, Vec3D dxB, Vec3D dxC, 
					    Vec3D &dnds){

  Vec3D V1, V2, Vp, dV1, dV2, dVp; 
  double sp, dsp;

  V1 = (xC - xA);
  V2 = (xB - xA);
    
  dV1 = dxC - dxA;
  dV2 = dxB - dxA;

  Vp = -(V1^V2);

  dVp = -(dV1^V2) - (V1^dV2);

  double d = Vp*Vp;

  double tmp = dVp*Vp;

  dnds = (1.0/sqrt(d)) * (dVp) - (Vp*tmp)/pow(d, 3.0/2.0);

}

//----------------------------------------------------------------------------

void DistIntersectorPhysBAM::updateXb(double epsilon){

  for(int i=0; i<numStNodes; ++i){

    for(int j=0; j<3; ++j){
      //Xs0[i][j]
      Xs[i][j]     = Xs0[i][j] + dXdSb[i][j]*epsilon;
      Xs_n[i][j]   = Xs[i][j];
      Xs_np1[i][j] = Xs_n[i][j];
      Xsdot[i][j]  = 0.0;
    }
  }

  initializePhysBAM();
  updatebc();

}

//----------------------------------------------------------------------------

double IntersectorPhysBAM::derivativeOFalpha(Vec3D  xA, Vec3D  xB, Vec3D  xC, 
					    Vec3D dxA, Vec3D dxB, Vec3D dxC, 
					    Vec3D  X1, Vec3D  X2){

  Vec3D  u =  xB -  xA;
  Vec3D du = dxB - dxA;

  Vec3D  v =  xC -  xA;
  Vec3D dv = dxC - dxA;

  Vec3D  n = u ^ v;
  Vec3D dn = (du ^ v) + (u ^ dv);
 
  Vec3D d = X1 - X2;

  Vec3D  w = X2 - xA;
  Vec3D dw =    -dxA;

  double tmp1 = (dn*w) + (n*dw);
  double tmp2 = n*d;
  double tmp3 = n*w;
  double tmp4 = dn*d;
  
  double dalpha_ds = -(tmp1/tmp2) + tmp3*tmp4/(tmp2*tmp2);

  return dalpha_ds;

}

//----------------------------------------------------------------------------

// double
// IntersectorPhysBAM::testAlpha(Vec3D Xi, Vec3D Xj, int n0, int n1, int n2){
//   Vec3D *XX  = (distIntersector.getStructPosition()).data(); 
//   Vec3D u = XX[n1] - XX[n0];
//   Vec3D v = XX[n2] - XX[n0];
//   Vec3D n = u ^ v;
//   Vec3D d = Xi - Xj;
//   Vec3D w = Xj - XX[n0];
//   double a = -n*w;
//   double b =  n*d;
//   double g = a/b;
//   return g; 
//}
//----------------------------------------------------------------------------

double 
IntersectorPhysBAM::testAlpha(Vec3D Xi, Vec3D Xj, Vec3D X0, Vec3D X1, Vec3D X2){

   Vec3D u = X1 - X0;
   Vec3D v = X2 - X0;
   Vec3D n = u ^ v;

   Vec3D d = Xi - Xj;
   Vec3D w = Xj - X0;

   double a = -n*w;
   double b =  n*d;

   double g = a/b;

   return g;

 }


//----------------------------------------------------------------------------
