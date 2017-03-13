/*
 * DynamicNodalTransfer.cpp
 *
 *  Created on: May 12, 2009
 *      Author: michel, kevin
 */
#include <iostream>
#include <FSI/DynamicNodalTransfer.h>
#include <FSI/CrackingSurface.h>
#include <IoData.h>
#include <Vector3D.h>
#include <MatchNode.h>
#include <StructExc.h>
#include <DistVector.h>
#include <assert.h>
#include <list>
#include <map>
#include <fstream>

#include "TsRestart.h"

//------------------------------------------------------------------------------

DynamicNodalTransfer::DynamicNodalTransfer(IoData& iod, Communicator &c, Communicator &sc, Timer *tim): com(c) , F(1),
                           fScale(iod.ref.rv.tforce), XScale(iod.ref.rv.tlength), UScale(iod.ref.rv.tvelocity),
                           tScale(iod.ref.rv.time), structure(iod,c,sc,tim), iod(iod)

{
  timer = tim;

  if (cracking()) {

    if (strcmp(iod.input.cracking,"") != 0 || 
	iod.input.restart_file_package[0] != 0) {

      char* fn;
      if (iod.input.restart_file_package[0] == 0) {
	int lns = strlen(iod.input.prefix)+strlen(iod.input.cracking)+1;
	fn = new char[lns];
	sprintf(fn,"%s%s",iod.input.prefix, iod.input.cracking);
      } else {

	fn = new char[256];
	char dummy[256];
	int lns = strlen(iod.input.prefix)+strlen(iod.input.restart_file_package)+1;
	char* fn2 = new char[lns];
	sprintf(fn2,"%s%s",iod.input.prefix, iod.input.restart_file_package);
	
	TsRestart::readRestartFileNames(fn2, dummy, dummy, dummy,
					fn, dummy, dummy, dummy,&c);
	delete [] fn2;
      }
      std::ifstream infile(fn, std::ios::binary);
      readCrackingData(infile);
      delete [] fn;
    }
  }

//  com.fprintf(stderr,"fscale = %e, XScale = %e, tScale = %e.\n", fScale, XScale, tScale);

{  Communication::Window<double> window(com, 1*sizeof(double), &dts);
  window.fence(true);
  structure.sendTimeStep(&window);
  window.fence(false);
}

{  Communication::Window<double> window2(com, 1*sizeof(double), &tMax);
  window2.fence(true);
  structure.sendMaxTime(&window2);
  window2.fence(false);
}

  algNum = structure.getAlgorithmNumber();
//  com.fprintf(stderr,"--- Received initial structure Time-step: %e, Final Time: %e\n", dts, tMax);
  dts /= tScale;
  tMax /= tScale;
  com.barrier();

  //initialize windows
  dt_tmax = new double[2];
  wintime = new Communication::Window<double> (com, 2*sizeof(double), (double*)dt_tmax);

  XandUdot = new double[2*3*structure.totalNodes];
  winDisp = new  Communication::Window<double> (com, 2*3*structure.nNodes*sizeof(double), (double *)XandUdot);

  std::pair<double *, int> embedded = structure.getTargetData();
  double *embeddedData = embedded.first;
  int length = embedded.second;
  winForce = new Communication::Window<double> (com, 3*length*sizeof(double), embeddedData);

  //get structure position
  com.barrier(); //for timing purpose
  winDisp->fence(true);
  structure.sendInitialPosition(winDisp);
  winDisp->fence(false);

  int N = structure.nNodes;
  for(int i=0; i<N; i++)
    for(int j=0; j<3; j++) {
      XandUdot[i*3+j] *= 1.0/XScale;
      XandUdot[(N+i)*3+j] *= 1.0/UScale;
    }

  structureSubcycling = (algNum==22) ? getStructSubcyclingInfo() : 0;

}

//------------------------------------------------------------------------------

DynamicNodalTransfer::~DynamicNodalTransfer() {
  if(wintime)  delete   wintime;
  if(winForce) delete   winForce;
  if(winDisp)  delete   winDisp;
  if(XandUdot) delete[] XandUdot;
  if(dt_tmax)  delete[] dt_tmax;
}

//------------------------------------------------------------------------------

void
DynamicNodalTransfer::sendForce()
{
  structure.clearForceVector();

  com.barrier(); //for timing purpose
  winForce->fence(true);
  winForce->accumulate((double *)F.data(), 0, 3*F.size(), 0, 0, Communication::Window<double>::Add);
  winForce->fence(false);

  structure.processReceivedForce();
}

//------------------------------------------------------------------------------

void
DynamicNodalTransfer::sendFluidSuggestedTimestep(double dtf0)
{
  dtf0 *= tScale;
  structure.sendFluidSuggestedTimestep(dtf0);
}

//------------------------------------------------------------------------------

void
DynamicNodalTransfer::updateInfo() {
  com.barrier(); //for timing purpose
  wintime->fence(true);
  structure.sendInfo(wintime);
  wintime->fence(false);
  dts   = dt_tmax[0];
  tMax  = dt_tmax[1];
//  fprintf(stderr,"*** CPU %d received dts = %e, tMax = %e (dimensional)\n", com.cpuNum(), dts, tMax);
  dts  /= tScale;
  tMax /= tScale;
}

//------------------------------------------------------------------------------

int
DynamicNodalTransfer::getNewCracking()
{
  if(!cracking()) {
    com.fprintf(stderr,"WARNING: Cracking is not considered in the structure simulation!\n");
    return 0;
  }
  int numConnUpdate = structure.getNewCracking();

  if(numConnUpdate) { //update "windows".
    int N = structure.nNodes;
//    if((N*2*3*sizeof(double) == winDisp->size())) {com.fprintf(stderr,"WEIRD!\n");exit(-1);}

    delete winDisp;
    winDisp = new Communication::Window<double> (com, 2*3*N*sizeof(double), (double *)XandUdot);
   
    delete winForce;
    std::pair<double *, int> embedded = structure.getTargetData();
    double *embeddedData = embedded.first;
    int length = embedded.second;
    if(length!=N) {com.fprintf(stderr,"WEIRD TOO!\n");exit(-1);}
    winForce = new Communication::Window<double> (com, 3*length*sizeof(double), embeddedData);
  }

  return numConnUpdate;
}

//------------------------------------------------------------------------------

void
DynamicNodalTransfer::getDisplacement()
{
  int N = numStNodes();
  if(N*2*3*sizeof(double) != winDisp->size()) {
    com.fprintf(stderr,"SOFTWARE BUG: length of winDisp (i.e. # nodes) is %d (should be %d)!\n", 
                        winDisp->size()/(2*3*sizeof(double)), N);
    exit(-1);
  }

  com.barrier(); //for timing purpose
  winDisp->fence(true);
  structure.sendDisplacement(winDisp);
  winDisp->fence(false);

  for(int i=0; i<N; i++) { 
    for(int j=0; j<3; j++) {
      XandUdot[i*3+j] *= 1.0/XScale;
      XandUdot[(N+i)*3+j] *= 1.0/UScale;
    }
  }
}

//------------------------------------------------------------------------------

int DynamicNodalTransfer::getStructSubcyclingInfo()
{
  int subcyc = structure.sendSubcyclingInfo();
  return subcyc;
}

//------------------------------------------------------------------------------

void
DynamicNodalTransfer::updateOutputToStructure(double dt, double dtLeft, SVec<double,3> &fs)
{
  if(F.size() != fs.size()) {
   // com.fprintf(stderr,"force vector in DynamicNodalTransfer resized (from %d to %d)!\n", F.size(), fs.size());
    F.resize(fs.size());
  }
  F = fs;
}

//------------------------------------------------------------------------------

EmbeddedStructure::EmbeddedStructure(IoData& iod, Communicator &comm, Communicator &strCom, Timer *tim) : com(comm), 
                                             tScale(iod.ref.rv.time), XScale(iod.ref.rv.tlength),
                                             UScale(iod.ref.rv.tvelocity), nNodes(0), nElems(0), elemType(3),
                                             cracking(0), X(0), X0(0), Tria(0), U(0), Udot(0), XandUdot(0), F(0), it(0), 
                                             structExc(0), mns(0), algNum(6), iod(iod) //A6
{
  timer = tim;

  // ----------------------------------
  //       User-Provided Info
  // ----------------------------------
  if(iod.problem.type[ProblemData::AERO])
    coupled = true;
  else if(iod.problem.type[ProblemData::FORCED])
    coupled = false;
  else {
    com.fprintf(stderr,"ERROR: Simulation type is not supported by the embedded framework.\n");
    exit(-1);
  }

  // ---- input files ----
  int sp = strlen(iod.input.prefix) + 1;
  meshFile = new char[sp + strlen(iod.input.embeddedSurface)];
  sprintf(meshFile,"%s%s", iod.input.prefix, iod.input.embeddedSurface);

  restartmeshFile = new char[sp + strlen(iod.input.embeddedpositions)];
  if(iod.input.embeddedpositions[0] != 0)
    sprintf(restartmeshFile,"%s%s", iod.input.prefix, iod.input.embeddedpositions);
  else //no restart position file provided
    restartmeshFile[0] = '\0'; 

  matcherFile = new char[sp + strlen(iod.input.match)];
  sprintf(matcherFile,"%s%s", iod.input.prefix, iod.input.match); 

  getSurfFromFEM = coupled && (!strlen(iod.input.embeddedSurface));
  if(getSurfFromFEM)
    com.fprintf(stderr,"- Using the embedded surface provided by structure code.\n");

  surfaceID = NULL;

  // ---- for forced-motion only ------
  if(!coupled) {
    if(iod.forced.type==ForcedData::HEAVING)
      mode = 1;
    else if(iod.forced.type==ForcedData::PITCHING)
      mode = 2;
    else if(iod.forced.type==ForcedData::VELOCITY)
      mode = 3;
    else if (iod.forced.type==ForcedData::DEFORMING)
      mode = 4;
    else if (iod.forced.type==ForcedData::SPIRALING)
      mode = 5;
    else if (iod.forced.type==ForcedData::ACOUSTICVISCOUSBEAM)
      mode = 97;
    else if (iod.forced.type==ForcedData::ACOUSTICBEAM)
      mode = 98;
    else if (iod.forced.type==ForcedData::DEBUGDEFORMING)
      mode = 99;
    else {
      com.fprintf(stderr,"ERROR: Forced motion type is not supported by the embedded framework.\n");
      exit(-1);
    }
  } else mode = -1;

  // NOTE: All variables stored in EmbeddedStructure must be dimensional! (unless the simulation itself is non-dim)
  tMax  = tScale*iod.ts.maxTime; //iod.ts.maxTime is already non-dimensionalized in IoData
  dt    = tScale*iod.forced.timestep;
  t0    = tScale*iod.restart.etime;
  omega = 2.0*acos(-1.0)*(1.0/tScale)*iod.forced.frequency;

  // for heaving
  if (mode!=99) {				
    dx    = iod.ref.rv.length*iod.forced.hv.ax;
    dy    = iod.ref.rv.length*iod.forced.hv.ay;
    dz    = iod.ref.rv.length*iod.forced.hv.az;  
  }																	
  else {															
	dx = dy = dz = iod.ref.rv.length*iod.forced.df.amplification;	
  }																	

  // for pitching
  alpha_in  = (acos(-1.0)*iod.forced.pt.alpha_in) / 180.0;  // initial angle of rotation
  alpha_max = (acos(-1.0)*iod.forced.pt.alpha_max) / 180.0;  // maximum angle of rotation
  x1[0] = iod.forced.pt.x11;  x1[1] = iod.forced.pt.y11;  x1[2] = iod.forced.pt.z11;
  x2[0] = iod.forced.pt.x21;  x2[1] = iod.forced.pt.y21;  x2[2] = iod.forced.pt.z21;

  beta_in  = (acos(-1.0)*iod.forced.pt.beta_in) / 180.0;  // initial angle of rotation
  beta_max = (acos(-1.0)*iod.forced.pt.beta_max) / 180.0;  // maximum angle of rotation
  y1[0] = iod.forced.pt.x12;  y1[1] = iod.forced.pt.y12;  y1[2] = iod.forced.pt.z12;
  y2[0] = iod.forced.pt.x22;  y2[1] = iod.forced.pt.y22;  y2[2] = iod.forced.pt.z22;

  for(int i=0; i<3; i++) { //get back to user-specified coordinates.
    x1[i] *= iod.ref.rv.length;
    x2[i] *= iod.ref.rv.length;
    y1[i] *= iod.ref.rv.length;
    y2[i] *= iod.ref.rv.length;
  }
  u = x2[0]-x1[0];  v = x2[1]-x1[1];  w = x2[2]-x1[2];
  // unit normals of axis of rotation //
  ix = u/sqrt(u*u+v*v+w*w);  iy = v/sqrt(u*u+v*v+w*w);  iz = w/sqrt(u*u+v*v+w*w);

  // for deforming data
  deformMeshFile = new char[strlen(iod.forced.df.positions) + 1];
  if(iod.forced.df.positions[0] != 0)
    sprintf(deformMeshFile,"%s", iod.forced.df.positions);
  else //no deforming data position file provided
    deformMeshFile[0] = '\0'; 
  Xd = 0;
  dXmax = 0;

  // for spiraling
  cableLen = iod.ref.rv.length*iod.forced.sp.xL;
  xbeg = iod.ref.rv.length*iod.forced.sp.x0;

  // ----------------------------------
  //               End
  // ----------------------------------

  // ---------------------------------
  //       F-S Communication
  // ---------------------------------
  if(coupled) {
    mns = new MatchNodeSet *[1];
    if(com.cpuNum() == 0 && !getSurfFromFEM)
      mns[0] = new MatchNodeSet(matcherFile);
    else 
      mns[0] = new MatchNodeSet();
    structExc = new StructExc(iod, mns, 6, &strCom, &com, 1);

    if(getSurfFromFEM) {//receive embedded surface from the structure code.
      bool crack;
      int  nStNodes, nStElems, totalStNodes, totalStElems;
      structExc->getEmbeddedWetSurfaceInfo(elemType, crack, nStNodes, nStElems); 

      // initialize cracking information
      if(crack) {
        structExc->getInitialCrackingSetup(totalStNodes, totalStElems);
        cracking = new CrackingSurface(elemType, nStElems, totalStElems, nStNodes, totalStNodes);
      } else {
        totalStNodes = nStNodes;
        totalStElems = nStElems;
      }

      // allocate memory for the node list
      totalNodes = totalStNodes;
      nNodes     = nStNodes;
      X = new (com) double[totalNodes][3];
      X0 = new (com) double[totalNodes][3];
      for(int i=0; i<totalNodes; i++)
        X[i][0] = X[i][1] = X[i][2] = X0[i][0] = X0[i][1] = X0[i][2] = 0.0;

      // allocate memory for the element topology list
      int tmpTopo[nStElems][4];
      switch (elemType) {
        case 3: //all triangles
          totalElems = totalStElems;
          nElems     = nStElems; 
          Tria = new (com) int[totalElems][3];
          structExc->getEmbeddedWetSurface(nNodes, (double*)X0, nElems, (int*)Tria, elemType);
          break;
        case 4: //quadrangles include triangles represented as degenerated quadrangles.
          structExc->getEmbeddedWetSurface(nNodes, (double*)X0, nStElems, (int*)tmpTopo, elemType);
          if(cracking) {
            totalElems = totalStElems*2;
            Tria = new (com) int[totalElems][3];
            nElems = cracking->splitQuads((int*)tmpTopo, nStElems, Tria);
          } else
            splitQuads((int*)tmpTopo, nStElems); //memory for Tria will be allocated
          break;
        default:
          com.fprintf(stderr,"ERROR: Element type (%d) of the wet surface not recognized! Must be 3 or 4.\n", elemType);
          exit(-1);
      }

      for(int i=0; i<nNodes; i++)
        for(int j=0; j<3; j++)
          X[i][j] = X0[i][j];

      if(com.cpuNum()==0) {
        mns[0]->autoInit(totalNodes); //in case of cracking, match all the nodes including inactive ones.
        structExc->updateMNS(mns);
      }
    }

    structExc->negotiate();
    structExc->getInfo();
    dt = tScale*structExc->getTimeStep();
    tMax = tScale*structExc->getMaxTime();
    algNum = structExc->getAlgorithmNumber();


    if(cracking) {
      getInitialCrack();
      cracking->setNewCrackingFlag(false);
    }
  }

  // ----------------------------------
  //               End
  // ----------------------------------


  // ----------------------------------
  //    Load Structure Mesh (at t=0)
  // ----------------------------------
  if(!getSurfFromFEM) {//otherwise mesh is already loaded.
    // load structure nodes and elements (from file).
    FILE *topFile = 0;
    topFile = fopen(meshFile,"r");
    if(!topFile) com.fprintf(stderr,"ERROR: Embedded surface mesh file could not be found!\n");

    char line[MAXLINE],key1[MAXLINE],key2[MAXLINE];

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
  	  if(num1<1) {com.fprintf(stderr,"ERROR: detected a node with index %d in the embedded surface file!\n",num1); exit(-1);}
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

    nNodes = totalNodes = nodeList.size();
    nElems = totalElems = elemList1.size();

    if(nNodes != maxIndex) {
      com.fprintf(stderr,"ERROR: The node set of the embedded surface have gap(s). \n");
      com.fprintf(stderr,"       Detected max index = %d, number of nodes = %d\n", maxIndex, nNodes);
      com.fprintf(stderr,"NOTE: Currently the node set of the embedded surface cannot have gaps. Moreover, the index must start from 1.\n");
      exit(-1);
    }

    X0 = new (com) double[totalNodes][3];
    X  = new (com) double[totalNodes][3];
    surfaceID = new int[totalNodes];

    std::list<Vec3D>::iterator it1;
    std::list<int>::iterator it2;
    std::list<int>::iterator it_1;
    std::list<int>::iterator it_2;
    std::list<int>::iterator it_3;
    std::list<int>::iterator it_4;

    it2=indexList.begin();
    for (it1=nodeList.begin(); it1!=nodeList.end(); it1++) {
      X[(*it2)-1][0] = X0[(*it2)-1][0] = (*it1)[0]; 
      X[(*it2)-1][1] = X0[(*it2)-1][1] = (*it1)[1]; 
      X[(*it2)-1][2] = X0[(*it2)-1][2] = (*it1)[2];
      it2++; 
    }

    for (int k=0; k<totalNodes; k++) {
      surfaceID[k] = 0;
    }

    Tria = new (com) int[totalElems][3];

    it_1 = elemList1.begin();
    it_2 = elemList2.begin();
    it_3 = elemList3.begin();
    it_4 = surfaceIDList.begin();
    for (int i=0; i<nElems; i++) {
      Tria[i][0] = *it_1;
      Tria[i][1] = *it_2;
      Tria[i][2] = *it_3;
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
    if (restartmeshFile[0] != 0) {
      FILE* resTopFile = fopen(restartmeshFile, "r");
      if(resTopFile==NULL) {com.fprintf(stderr, "restart topFile doesn't exist.\n"); exit(1);}
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
      if (ndMax!=totalNodes) {
        com.fprintf(stderr,"ERROR: number of nodes in restart top-file is wrong.\n");
        exit(1);
      }

      for(int i=0; i<totalNodes; i++)
        X[i][0] = X[i][1] = X[i][2] = 0.0;
      
      for (it2=nodeList2.begin(); it2!=nodeList2.end(); it2++) {
        X[it2->first-1][0] = (it2->second)[0];
        X[it2->first-1][1] = (it2->second)[1];
        X[it2->first-1][2] = (it2->second)[2];
      }

      fclose(resTopFile);
    }

    if (mode==4) { //deforming data
      int nInputs;
      char c1[200], c2[200], c3[200];
      // load deforming data
      if (deformMeshFile[0] != 0) {
        FILE* defoTopFile = fopen(deformMeshFile, "r");
        if(defoTopFile==NULL) {com.fprintf(stderr, "Deforming data topFile doesn't exist.\n"); exit(1);}
        int ndMax = 0, ndMax2 = 0;
        std::list<std::pair<int,Vec3D> > nodeList2;
        std::list<std::pair<int,Vec3D> >::iterator it2;

        while(1) {
          nInputs = fscanf(defoTopFile,"%s", c1);
          if(nInputs!=1) break;    
          char *endptr;
          num1 = strtol(c1, &endptr, 10);
          if(endptr == c1) break;

          int toto = fscanf(defoTopFile,"%lf %lf %lf\n", &x1, &x2, &x3);
          nodeList2.push_back(std::pair<int,Vec3D>(num1,Vec3D(x1,x2,x3)));
          ndMax = std::max(num1, ndMax);
        }
        if (ndMax!=totalNodes) {
          com.fprintf(stderr,"ERROR: number of nodes in Deforming data top-file is wrong.\n");
          exit(1);
        }

        Xd  = new (com) double[totalNodes][3];
        dXmax  = new (com) double[totalNodes][3];

        for(int i=0; i<totalNodes; i++)
          Xd[i][0] = Xd[i][1] = Xd[i][2] = 0.0;
        
        for (it2=nodeList2.begin(); it2!=nodeList2.end(); it2++) {
          Xd[it2->first-1][0] = (it2->second)[0];
          Xd[it2->first-1][1] = (it2->second)[1];
          Xd[it2->first-1][2] = (it2->second)[2];
          dXmax[it2->first-1][0] = iod.forced.df.amplification*(Xd[it2->first-1][0] - X0[it2->first-1][0]);
          dXmax[it2->first-1][1] = iod.forced.df.amplification*(Xd[it2->first-1][1] - X0[it2->first-1][1]);
          dXmax[it2->first-1][2] = iod.forced.df.amplification*(Xd[it2->first-1][2] - X0[it2->first-1][2]);
        }

        fclose(defoTopFile);
      }
    }

  }

  makerotationownership(iod);

  // ----------------------------------
  //               End
  // ----------------------------------

  // allocate memory for other stuff...
  U = new (com) double[totalNodes][3];
  Udot = new (com) double[totalNodes][3];
  XandUdot = new (com) double[totalNodes*2][3];
  F = new (com) double[totalNodes][3];

  // prepare distinfo for struct nodes
  if(coupled) {
    int *locToGlob = new int[1];
    locToGlob[0] = 0;
    di = new DistInfo(1, 1, 1, locToGlob, &com);
    di->setLen(0,totalNodes);
    di->finalize(false);
  }

  timeStepOffset = iod.forced.tsoffset;
}

//------------------------------------------------------------------------------

EmbeddedStructure::~EmbeddedStructure()
{
  if(X0)       operator delete[] (X0, com);
  if(X)        operator delete[] (X, com);
  if(Xd)       operator delete[] (Xd, com);
  if(dXmax)    operator delete[] (dXmax, com);
  if(Tria)     operator delete[] (Tria, com);
  if(U)        operator delete[] (U, com);
  if(Udot)     operator delete[] (Udot, com);
  if(XandUdot) operator delete[] (XandUdot, com);
  if(F)        operator delete[] (F, com);
  if(surfaceID) delete[] surfaceID;
  if(rotOwn) delete[] rotOwn;

  if(structExc) delete structExc;  
  if(mns)      {delete mns[0]; delete [] mns;}
  delete[] meshFile;
  if(restartmeshFile) delete[] restartmeshFile;
  if(deformMeshFile) delete[] deformMeshFile;
  delete[] matcherFile;

  if(coupled) delete di;
}

//------------------------------------------------------------------------------

pair<double*, int>
EmbeddedStructure::getTargetData() 
{
  //clear the force.
  for (int i=0; i<nNodes; i++) {
    F[i][0] = 0.0;
    F[i][1] = 0.0;
    F[i][2] = 0.0;
  }

  return pair<double*,int>((double*)F,nNodes);
}

//------------------------------------------------------------------------------
void EmbeddedStructure::makerotationownership(IoData &iod) {
  map<int,SurfaceData *> &surfaceMap = iod.surfaces.surfaceMap.dataMap;
  map<int,SurfaceData *>::iterator it = surfaceMap.begin();
  rotationMap = &(iod.forced.vel.rotationMap.dataMap);
  rotOwn = 0;
  
  int numRotSurfs = 0;
  int numTransWalls = 0;
  while(it != surfaceMap.end()) {
    map<int,RotationData*>::iterator it1 = rotationMap->find(it->second->forceID);
    if(it1!=rotationMap->end()) {
      if(it1->second->infRadius) {
        numTransWalls++;
        com.fprintf(stderr," ... surface %2d is ``translating''\n",it->first, it1->first);
        com.fprintf(stderr,"     -> uniform velocity V = %3.2e in direction %3.2e %3.2e %3.2e\n",
                     it1->second->omega, it1->second->nx,it1->second->ny,it1->second->nz);
      } else {
        numRotSurfs++;
        com.fprintf(stderr," ... surface %2d is ``rotating'' in DynamicNodalTranser using rotation data %2d\n",it->first, it1->first);
        com.fprintf(stderr,"     -> omega = %3.2e, rotation axis = %3.2e %3.2e %3.2e\n",
                     it1->second->omega, it1->second->nx,it1->second->ny,it1->second->nz);
      }
    }
    it++;
  }

  if(numRotSurfs || numTransWalls) {
    rotOwn = new int[nNodes];

    for (int k=0; k<nNodes; k++) {
      rotOwn[k] = -1;

      map<int,SurfaceData *>::iterator it = surfaceMap.find(surfaceID[k]);
      if(it != surfaceMap.end()) {
         int rotID = it->second->forceID; // = -1 if not defined in input file
         if ((rotOwn[k] != -1 && rotOwn[k] != rotID) ||
             (rotOwn[k] != -1 && rotOwn[k] != rotID) ||
             (rotOwn[k] != -1 && rotOwn[k] != rotID))  {
           com.fprintf(stderr, " ... WARNING: Embedded Node %d associated to more than 1 Rotation ID\n",k);
         }
         rotOwn[k] = rotID;
         rotOwn[k] = rotID;
         rotOwn[k] = rotID;
      }
    }
  }
}
//------------------------------------------------------------------------------
void EmbeddedStructure::updaterotation(double time) {
  if (rotOwn)  { 
    for (int k=0; k<nNodes; k++) {
      if (rotOwn[k]>=0) {    // node belongs to a (potential) "rotating" surface
        map<int,RotationData *>::iterator it =  rotationMap->find(rotOwn[k]);
        if(it != rotationMap->end()) { // the rotation data have been defined
	  if(it->second->infRadius == RotationData::TRUE) {
            double vel = it->second->omega;
	    U[k][0] = vel*it->second->nx*time;
	    U[k][1] = vel*it->second->ny*time;
	    U[k][2] = vel*it->second->nz*time;
	    Udot[k][0] = vel*it->second->nx;
	    Udot[k][1] = vel*it->second->ny;
	    Udot[k][2] = vel*it->second->nz;
	  } 
          else {
	    double xd = X0[k][0] - it->second->x0;
	    double yd = X0[k][1] - it->second->y0;
	    double zd = X0[k][2] - it->second->z0;
	    double theta = it->second->omega*time;
	    double costheta = cos(theta);
	    double sintheta = sin(theta);
	    double ix = it->second->nx;
	    double iy = it->second->ny;
	    double iz = it->second->nz;

	    double dx[3] = {0.,0.,0.};
            dx[0] += (costheta + (1 - costheta) * ix * ix) * xd;
            dx[0] += ((1 - costheta) * ix * iy - iz * sintheta) * yd;
            dx[0] += ((1 - costheta) * ix * iz + iy * sintheta) * zd;

            dx[1] += ((1 - costheta) * ix * iy + iz * sintheta) * xd;
            dx[1] += (costheta + (1 - costheta) * iy * iy) * yd;
            dx[1] += ((1 - costheta) * iy * iz - ix * sintheta) * zd;

            dx[2] += ((1 - costheta) * ix * iz - iy * sintheta) * xd;
            dx[2] += ((1 - costheta) * iy * iz + ix * sintheta) * yd;
            dx[2] += (costheta + (1 - costheta) * iz * iz) * zd;

            dx[0] += it->second->x0;
            dx[1] += it->second->y0;
            dx[2] += it->second->z0;

	    U[k][0] = dx[0] - X0[k][0]; 
	    U[k][1] = dx[1] - X0[k][1];
	    U[k][2] = dx[2] - X0[k][2];

	    double ox = it->second->omega*it->second->nx;
	    double oy = it->second->omega*it->second->ny;
	    double oz = it->second->omega*it->second->nz;
	    xd = dx[0] - it->second->x0;
	    yd = dx[1] - it->second->y0;
	    zd = dx[2] - it->second->z0;
	    Udot[k][0] = oy*zd-oz*yd;
	    Udot[k][1] = oz*xd-ox*zd;
	    Udot[k][2] = ox*yd-oy*xd;
	  }
	} 
        else  { // no rotation data
          U[k][0] = 0.;
          U[k][1] = 0.;
          U[k][2] = 0.;
          Udot[k][0] = 0.;
          Udot[k][1] = 0.;
          Udot[k][2] = 0.;
	}
      }
      else  { // no rotation data
        U[k][0] = 0.;
        U[k][1] = 0.;
        U[k][2] = 0.;
        Udot[k][0] = 0.;
        Udot[k][1] = 0.;
        Udot[k][2] = 0.;
      }
    }
  }
}
//------------------------------------------------------------------------------

void
EmbeddedStructure::clearForceVector()
{
  //clear the force.
  for (int i=0; i<nNodes; i++) {
    F[i][0] = 0.0;
    F[i][1] = 0.0;
    F[i][2] = 0.0;
  }
}

//------------------------------------------------------------------------------

void
EmbeddedStructure::sendTimeStep(Communication::Window<double> *window)
{
  if(com.cpuNum()>0) return; // only proc #1 sends the time.
  //std::cout << "Sending the timestep (" << dt << ") to fluid " << std::endl;
{
  for(int i = 0; i < com.size(); ++i) {
    window->put(&dt, 0, 1, i, 0);
  }
}
}

//------------------------------------------------------------------------------

void
EmbeddedStructure::sendMaxTime(Communication::Window<double> *window)
{
  if(com.cpuNum()>0) return; // only proc #1 sends the time.
  //std::cout << "Sending the max time (" << tMax << ") to fluid " << std::endl;
{
  for(int i = 0; i < com.size(); ++i) {
    window->put(&tMax, 0, 1, i, 0);
  }
}
}

//------------------------------------------------------------------------------

void
EmbeddedStructure::sendInfo(Communication::Window<double> *window)
{
  if(coupled)
    structExc->getInfo();
  if(com.cpuNum()>0) return; // only proc. #1 will send.

  if(coupled) {
    dt = tScale*structExc->getTimeStep();
    tMax = tScale*structExc->getMaxTime();
  } /* else: nothing to be done */

  dt_tmax[0] = dt; dt_tmax[1] = tMax;
  for(int i = 0; i < com.size(); ++i)
    window->put((double*)dt_tmax, 0, 2, i, 0);
}

//------------------------------------------------------------------------------

void
EmbeddedStructure::sendInitialPosition(Communication::Window<double> *window)
{
  if(com.cpuNum()>0) return; // only proc. #1 will send.

  for(int i=0; i<nNodes; i++)
    for(int j=0; j<3; j++) {
      XandUdot[i][j] = X[i][j];
      XandUdot[(i+nNodes)][j] = Udot[i][j];
    }

  for(int i = 0; i < com.size(); ++i)
    window->put((double*)XandUdot, 0, 2*3*nNodes, i, 0);
}

//------------------------------------------------------------------------------

void
EmbeddedStructure::sendDisplacement(Communication::Window<double> *window)
{
  if(coupled) {
     DistSVec<double,3> Y0(*di, X);
     DistSVec<double,3> V(*di, U);
     DistSVec<double,3> Ydot(*di, Udot);
     DistSVec<double,3> Y(*di);
     Y = Y0; //KW: as long as Y = Y0, it doesn't matter if Y0 is X or X0, or anything else...
     structExc->getDisplacement(Y0, Y, Ydot, V, true);
     V = XScale*V;
     Ydot = UScale*Ydot;
  }
  if(com.cpuNum()>0) return; // only proc. #1 will send.

  it++;
  if(!coupled) {
    double time;
    time = t0 + dt*((double)it + timeStepOffset);
     
    if (mode==1) //heaving
      for(int i=0; i < nNodes; ++i) {
        U[i][0] = (1-cos(omega*time))*dx;
        U[i][1] = (1-cos(omega*time))*dy;
        U[i][2] = (1-cos(omega*time))*dz;
        Udot[i][0] = dx*omega*sin(omega*time);
        Udot[i][1] = dy*omega*sin(omega*time); 
        Udot[i][2] = dz*omega*sin(omega*time); 
      }
    else if (mode==2) { //pitching
      double theta = alpha_in + alpha_max*sin(omega*time);
      double costheta = cos(theta);
      double sintheta = sin(theta);

      double phi = beta_in + beta_max*sin(omega*time); 
      double cosphi = cos(phi);
      double sinphi = sin(phi);

// Rotate the axis of 2nd rotation about the first rotation axis
      double p[3], yy1[3], yy2[3];

      p[0] = y1[0] - x1[0];                                         
      p[1] = y1[1] - x1[1];                                         
      p[2] = y1[2] - x1[2];                                         

      yy1[0] = 0.0; yy1[1] = 0.0; yy1[2] = 0.0;           
      yy1[0] += (costheta + (1 - costheta) * ix * ix) * p[0];         
      yy1[0] += ((1 - costheta) * ix * iy - iz * sintheta) * p[1];    
      yy1[0] += ((1 - costheta) * ix * iz + iy * sintheta) * p[2];    

      yy1[1] += ((1 - costheta) * ix * iy + iz * sintheta) * p[0];   
      yy1[1] += (costheta + (1 - costheta) * iy * iy) * p[1];         
      yy1[1] += ((1 - costheta) * iy * iz - ix * sintheta) * p[2];    

      yy1[2] += ((1 - costheta) * ix * iz - iy * sintheta) * p[0];    
      yy1[2] += ((1 - costheta) * iy * iz + ix * sintheta) * p[1];    
      yy1[2] += (costheta + (1 - costheta) * iz * iz) * p[2];         

      yy1[0] += x1[0];                                                  
      yy1[1] += x1[1];                                                   
      yy1[2] += x1[2];                                                 

      p[0] = y2[0] - x1[0];                                    
      p[1] = y2[1] - x1[1];                                    
      p[2] = y2[2] - x1[2];                                    

      yy2[0] = 0.0; yy2[1] = 0.0; yy2[2] = 0.0;                         
      yy2[0] += (costheta + (1 - costheta) * ix * ix) * p[0];         
      yy2[0] += ((1 - costheta) * ix * iy - iz * sintheta) * p[1];   
      yy2[0] += ((1 - costheta) * ix * iz + iy * sintheta) * p[2];   

      yy2[1] += ((1 - costheta) * ix * iy + iz * sintheta) * p[0];    
      yy2[1] += (costheta + (1 - costheta) * iy * iy) * p[1];        
      yy2[1] += ((1 - costheta) * iy * iz - ix * sintheta) * p[2];  

      yy2[2] += ((1 - costheta) * ix * iz - iy * sintheta) * p[0];    
      yy2[2] += ((1 - costheta) * iy * iz + ix * sintheta) * p[1];   
      yy2[2] += (costheta + (1 - costheta) * iz * iz) * p[2];      

      yy2[0] += x1[0];                                                
      yy2[1] += x1[1];                                             
      yy2[2] += x1[2];                                       

// unit normals of axis of 2nd rotation //

      u = yy2[0]-yy1[0];                                      
      v = yy2[1]-yy1[1];                                  
      w = yy2[2]-yy1[2];                            

      double jx = u/sqrt(u*u+v*v+w*w);
      double jy = v/sqrt(u*u+v*v+w*w);
      double jz = w/sqrt(u*u+v*v+w*w);

      for(int i=0; i<nNodes; ++i) {
        p[0] = X0[i][0] - x1[0];
        p[1] = X0[i][1] - x1[1];
        p[2] = X0[i][2] - x1[2];

        U[i][0] = 0.0; U[i][1] = 0.0; U[i][2] = 0.0;              
        U[i][0] += (costheta + (1 - costheta) * ix * ix) * p[0];
        U[i][0] += ((1 - costheta) * ix * iy - iz * sintheta) * p[1];
        U[i][0] += ((1 - costheta) * ix * iz + iy * sintheta) * p[2];

        U[i][1] += ((1 - costheta) * ix * iy + iz * sintheta) * p[0];
        U[i][1] += (costheta + (1 - costheta) * iy * iy) * p[1];
        U[i][1] += ((1 - costheta) * iy * iz - ix * sintheta) * p[2];

        U[i][2] += ((1 - costheta) * ix * iz - iy * sintheta) * p[0];
        U[i][2] += ((1 - costheta) * iy * iz + ix * sintheta) * p[1];
        U[i][2] += (costheta + (1 - costheta) * iz * iz) * p[2];

        U[i][0] += x1[0];
        U[i][1] += x1[1];
        U[i][2] += x1[2];

        p[0] = U[i][0] -  yy1[0];                               
        p[1] = U[i][1] -  yy1[1];                               
        p[2] = U[i][2] -  yy1[2];                               

        U[i][0] = 0.0; U[i][1] = 0.0; U[i][2] = 0.0;              
        U[i][0] += (cosphi + (1 - cosphi) * jx * jx) * p[0];            
        U[i][0] += ((1 - cosphi) * jx * jy - jz * sinphi) * p[1];      
        U[i][0] += ((1 - cosphi) * jx * jz + jy * sinphi) * p[2];    

        U[i][1] += ((1 - cosphi) * jx * jy + jz * sinphi) * p[0];       
        U[i][1] += (cosphi + (1 - cosphi) * jy * jy) * p[1];           
        U[i][1] += ((1 - cosphi) * jy * jz - jx * sinphi) * p[2];     

        U[i][2] += ((1 - cosphi) * jx * jz - jy * sinphi) * p[0];   
        U[i][2] += ((1 - cosphi) * jy * jz + jx * sinphi) * p[1];  
        U[i][2] += (cosphi + (1 - cosphi) * jz * jz) * p[2];           

        U[i][0] += yy1[0];                                            
        U[i][1] += yy1[1];                                        
        U[i][2] += yy1[2];                                     

        for(int j=0; j<3; j++)
          Udot[i][j] = (U[i][j]-X[i][j])/dt;

        for(int j=0; j<3; j++)
          U[i][j] -= X0[i][j];
      }
    }
    else if (mode==3) //heaving with a constant velocity (in this case dx dy dz are velocity
//      for(int i=0; i < nNodes; ++i) {
//        U[i][0] = time*dx;
//        U[i][1] = time*dy;
//        U[i][2] = time*dz;
//        Udot[i][0] = dx;
//        Udot[i][1] = dy;
//        Udot[i][2] = dz;
//      }
      updaterotation(time);
    else if (mode==4) //deforming data
    {
      for(int i=0; i<nNodes; ++i) {
        for(int j=0; j<3; j++) {
          U[i][j] = sin(omega*time)*dXmax[i][j];
          Udot[i][j] = omega*cos(omega*time)*dXmax[i][j];
	}

      }
    }
    else if (mode==5) //spiraling data
    {
      double Rcurv = cableLen/(omega*time);
      for(int i=0; i<nNodes; ++i) {
        if ( X0[i][0] >= xbeg ) {
          double xloc = X0[i][0] - xbeg;
          double theta = xloc/Rcurv;
          U[i][0] = (Rcurv + X0[i][1])*sin(theta) + xbeg;
          U[i][1] = (Rcurv + X0[i][1])*cos(theta) - Rcurv;
          U[i][2] = X0[i][2];
        }
        else {
          U[i][0] = X0[i][0];
          U[i][1] = X0[i][1];
          U[i][2] = X0[i][2];
        }

        for(int j=0; j<3; j++)
          Udot[i][j] = (U[i][j]-X[i][j])/dt;

        for(int j=0; j<3; j++)
          U[i][j] -= X0[i][j];
      }
    }
    else if (mode==97) //deforming data
    {
      for(int i=0; i<nNodes; ++i) {
        U[i][0] = U[i][2] = 0.0;
        Udot[i][0] = Udot[i][2] = 0.0;

        ExactSolution::AcousticViscousBeamStructure(iod,X0[i][0],X0[i][1],
                                                    X0[i][2], time/tScale,
                                                    U[i][1], Udot[i][1]);
        Udot[i][1] /= tScale;
      }
    }


    else if (mode==98) //deforming data
    {
      for(int i=0; i<nNodes; ++i) {
        U[i][0] = U[i][2] = 0.0;
        Udot[i][0] = Udot[i][2] = 0.0;

        ExactSolution::AcousticBeamStructure(iod,X0[i][0],X0[i][1],
                                             X0[i][2], time/tScale,
                                             U[i][1], Udot[i][1]);
        Udot[i][1] /= tScale;
      }
    }

    else if (mode==99) { // for debugging use.
	  bool shrinking_sphere = true;	
	  if (shrinking_sphere) {
//		for (int i=0; i<nNodes; ++i) {
//		  for (int j=0; j<3; j++) {
//			U[i][j] = 0.0;
//			Udot[i][j] = 0.0;
//		  }
//		}
		double tt = time;//t0+dt*(double)(it-1);
		double dist, dir[3];
        for(int i=0; i < nNodes; ++i) {
          dist = sqrt(X0[i][0]*X0[i][0]+X0[i][1]*X0[i][1]+X0[i][2]*X0[i][2]);
          for(int j=0; j<3; j++) {
            dir[j] = X0[i][j]/dist; 
            U[i][j] = -tt*tt*tt*tt/8.0*dir[j]*dx; 
            Udot[i][j] = -tt*tt*tt/2.0*dir[j]*dx;
          }
        }
	  }	
	  else { 
        for(int i=0; i < nNodes; ++i) { // expand / shrink the structure in y-z plane w.r.t. the origin
          double cosTheta = X[i][1]/sqrt(X[i][1]*X[i][1]+X[i][2]*X[i][2]);
          double sinTheta = X[i][2]/sqrt(X[i][1]*X[i][1]+X[i][2]*X[i][2]);
          U[i][0] = 0.0;
          U[i][1] = dy*time*cosTheta;
          U[i][2] = dz*time*sinTheta;
          Udot[i][0] = 0.0;
          Udot[i][1] = dy*cosTheta;
          Udot[i][2] = dz*sinTheta;
        }
	  }
	}
  }

  for(int i=0; i<nNodes; i++) {
//    fprintf(stderr,"U %d %e %e %e\n", i+1, U[i][0], U[i][1], U[i][2]);
//    fprintf(stderr,"Udot %d %e %e %e\n", i+1, Udot[i][0], Udot[i][1], Udot[i][2]);
    for(int j=0; j<3; j++) {
      X[i][j] = X0[i][j] + U[i][j];
      XandUdot[i][j] = X[i][j];
      XandUdot[(i+nNodes)][j] = Udot[i][j];
    }
  }

  for(int i = 0; i < com.size(); ++i) 
    window->put((double*)XandUdot, 0, 2*3*nNodes, i, 0);
}

//------------------------------------------------------------------------------

int
EmbeddedStructure::sendSubcyclingInfo(/*Communication::Window<int> *window*/)
{
/*  int subcyc = structExc->getSubcyclingInfo();
  if(com.cpuNum()>0) return; // only proc. #1 will send.

  for(int i = 0; i < com.size(); ++i)
    window->put((int*)&subcyc, 0, 1, i, 0);
*/
  return structExc->getSubcyclingInfo();
}

//------------------------------------------------------------------------------

void
EmbeddedStructure::processReceivedForce()
{ 
  if(coupled) {

    DistSVec<double,3> f(*di, F);
  
    structExc->sendForce(f);
  }
}

//------------------------------------------------------------------------------

void
EmbeddedStructure::sendFluidSuggestedTimestep(double dtf0)
{
  structExc->sendFluidSuggestedTimestep(dtf0);
}

//------------------------------------------------------------------------------

void
EmbeddedStructure::splitQuads(int* quads, int nStElems)
{ 
  int nTrias = 0;
  for(int i=0; i<nStElems; i++)
    if(quads[i*4+2]==quads[i*4+3])
      nTrias += 1;
    else 
      nTrias += 2;

  Tria = new (com) int[nTrias][3];

  int count = 0;
  for(int i=0; i<nStElems; i++) { 
    Tria[count][0] = quads[i*4];
    Tria[count][1] = quads[i*4+1];
    Tria[count][2] = quads[i*4+2];
    count++;

    if(quads[i*4+2]==quads[i*4+3])
      continue;

    Tria[count][0] = quads[i*4];
    Tria[count][1] = quads[i*4+2];
    Tria[count][2] = quads[i*4+3];
    count++;
  }

  if(count!=nTrias) {com.fprintf(stderr,"Software bug in FSI/EmbeddedStructure/splitQuad!\n");exit(-1);}
  nElems = totalElems = nTrias;
}

//------------------------------------------------------------------------------

void
EmbeddedStructure::getInitialCrack()
{
  int newNodes, numConnUpdate, numLSUpdate;
  bool need2update = structExc->getNewCrackingStats(numConnUpdate, numLSUpdate, newNodes); //inputs will be modified.
  if(!need2update) return; //Nothing new :)

  // get initial phantom nodes.
  structExc->getInitialPhantomNodes(newNodes,X,nNodes);
    //NOTE: nNodes will be updated in "getNewCracking"

  // get initial phantom elements (topo change).
  getNewCracking(numConnUpdate, numLSUpdate, newNodes);
}

//------------------------------------------------------------------------------

void
EmbeddedStructure::getNewCracking(int numConnUpdate, int numLSUpdate, int newNodes)
{
  if(numConnUpdate<1)  return;

  int phantElems[5*numConnUpdate]; // elem.id and node id.
  double phi[4*numLSUpdate];
  int phiIndex[numLSUpdate];
  int new2old[std::max(1,newNodes*2)]; //in the case of Element Deletion, newNodes might be 0
 
  structExc->getNewCracking(numConnUpdate, numLSUpdate, phantElems, phi, phiIndex, new2old, newNodes); 

  if(elemType!=4) {com.fprintf(stderr,"ERROR: only support quadrangles for cracking!\n");exit(1);} 
  nNodes += newNodes;
  nElems += cracking->updateCracking(numConnUpdate, numLSUpdate, phantElems, phi, phiIndex, Tria, nNodes, new2old, newNodes);
  if(nElems!=cracking->usedTrias()) {
    com.fprintf(stderr,"ERROR: inconsistency in the number of used triangles. (Software bug.)\n");exit(-1);}

  if(com.cpuNum()==0 && newNodes)
    structExc->updateNumStrNodes(nNodes); //set numStrNodes to nNodes in structExc and mns
}

//------------------------------------------------------------------------------

int
EmbeddedStructure::getNewCracking()
{
  int newNodes, numConnUpdate, numLSUpdate;
  bool need2update = structExc->getNewCrackingStats(numConnUpdate, numLSUpdate, newNodes); //inputs will be modified.
  if(!need2update) {assert(numConnUpdate==0); return 0;} //Nothing new :)
  getNewCracking(numConnUpdate, numLSUpdate, newNodes);
  return numConnUpdate;
}

//------------------------------------------------------------------------------

void DynamicNodalTransfer::writeCrackingData(std::ofstream& restart_file) const {

  structure.writeCrackingData(restart_file);
}

void DynamicNodalTransfer::readCrackingData(std::ifstream& restart_file) {

  structure.readCrackingData(restart_file);
}

void EmbeddedStructure::writeCrackingData(std::ofstream& restart_file) const {

  if (!cracking)
    return;

  restart_file.write(reinterpret_cast<const char*>(Tria), sizeof(int)*totalElems*3);
  cracking->writeCrackingData(restart_file);
}

void EmbeddedStructure::readCrackingData(std::ifstream& restart_file){

  if (!cracking)
    return;

  restart_file.read(reinterpret_cast<char*>(Tria), sizeof(int)*totalElems*3);
  cracking->readCrackingData(restart_file);
}
