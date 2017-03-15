#include <DistGeoState.h>

#include <DistInfo.h>
#include <IoData.h>
#include <TimeData.h>
#include <Domain.h>
#include <GeoData.h>
#include <GeoState.h>
#include <Vector3D.h>
#include <DistVector.h>
#include <Communicator.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>

//------------------------------------------------------------------------------

DistGeoState::DistGeoState(IoData &ioData, Domain *dom) : data(ioData), domain(dom)
{

  numLocSub = domain->getNumLocSub();
  com = domain->getCommunicator();
  lscale = ioData.ref.rv.tlength;
  oolscale = 1.0 / lscale;

  if (data.use_n) {
    Xn = new DistSVec<double,3>(domain->getNodeDistInfo());
    ctrlVol_n = new DistVec<double>(domain->getNodeDistInfo());
    // Initialize the values
    *Xn = 0.0;
    *ctrlVol_n = 0.0;
  }
  else {
    Xn = 0;
    ctrlVol_n = 0;
  }

  if (data.use_nm1) {
    Xnm1 = new DistSVec<double,3>(domain->getNodeDistInfo());
    ctrlVol_nm1 = new DistVec<double>(domain->getNodeDistInfo());
    // Initialize the values
    *Xnm1 = 0.0;
    *ctrlVol_nm1 = 0.0;
  }
  else {
    Xnm1 = 0;
    ctrlVol_nm1 = 0;
  }

  if (data.use_nm2) {
    Xnm2 = new DistSVec<double,3>(domain->getNodeDistInfo());
    ctrlVol_nm2 = new DistVec<double>(domain->getNodeDistInfo());
    // Initialize the values
    *Xnm2 = 0.0;
    *ctrlVol_nm2 = 0.0;
  }
  else {
    Xnm2 = 0;
    ctrlVol_nm2 = 0;
  }

  Xdot = new DistSVec<double,3>(domain->getNodeDistInfo());

// for mesh motion (ALE) with RK2 only
  if (data.use_save) { //true only if explicit with moving mesh
    ctrlVol_save = new DistVec<double>(domain->getNodeDistInfo());
    Xsave = new DistSVec<double,3>(domain->getNodeDistInfo());
  }else{
    ctrlVol_save = 0;
    Xsave = 0;
  }
// end ALE-RK2

  d2wall = new DistVec<double>(domain->getNodeDistInfo());
  if (ioData.input.d2wall[0] != 0) {
    char* name = new char[strlen(ioData.input.prefix) + strlen(ioData.input.d2wall) + 1];
    sprintf(name, "%s%s", ioData.input.prefix, ioData.input.d2wall);
    DistSVec<double,1> d2w(d2wall->info(), reinterpret_cast<double (*)[1]>(d2wall->data()));
    domain->readVectorFromFile(name, 0, 0, d2w, &oolscale);
    delete [] name;
  }
  else
    *d2wall = 0.0;

  edgeNorm    = new DistVec<Vec3D>(domain->getEdgeDistInfo());
  edgeNormVel = new DistVec<double>(domain->getEdgeDistInfo());
  edgeNorm_nm1    = 0;
  edgeNormVel_nm1 = 0;
  edgeNorm_nm2    = 0;
  edgeNormVel_nm2 = 0;

  faceNorm    = new DistVec<Vec3D>(domain->getFaceNormDistInfo());
  faceNormVel = new DistVec<double>(domain->getFaceNormDistInfo());
  faceNorm_nm1    = 0;
  faceNormVel_nm1 = 0;
  faceNorm_nm2    = 0;
  faceNormVel_nm2 = 0;
  
  inletNodeNorm = new DistVec<Vec3D>(domain->getInletNodeDistInfo());
  numFaceNeighb = new DistVec<int>(domain->getInletNodeDistInfo());


// Included (MB)
  if (ioData.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_ ||
      ioData.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_ ||
	  ioData.problem.alltype == ProblemData::_ROM_SHAPE_OPTIMIZATION_ ||
	  ioData.problem.alltype == ProblemData::_SENSITIVITY_ANALYSIS_){
    optFlag = 1;
    Xsa = new DistSVec<double,3>(domain->getNodeDistInfo());
    dXsa = new DistSVec<double,3>(domain->getNodeDistInfo());
    dEdgeNorm = new DistVec<Vec3D>(domain->getEdgeDistInfo());
    dFaceNorm = new DistVec<Vec3D>(domain->getFaceNormDistInfo());
    dEdgeNormVel = new DistVec<double>(domain->getEdgeDistInfo());
    dFaceNormVel = new DistVec<double>(domain->getFaceNormDistInfo());
    *dFaceNormVel = 0.0;
  }
  else {
    optFlag = 0;
    Xsa = 0;
    dXsa = 0;
    dEdgeNorm = 0;
    dFaceNorm = 0;
    dEdgeNormVel = 0;
    dFaceNormVel = 0;
  }


  if (data.typeNormals == DGCLData::IMPLICIT_SECOND_ORDER_GCL) {
    edgeNorm_nm1    = new DistVec<Vec3D>(domain->getEdgeDistInfo());
    edgeNormVel_nm1 = new DistVec<double>(domain->getEdgeDistInfo());
    faceNorm_nm1    = new DistVec<Vec3D>(domain->getFaceNormDistInfo());
    faceNormVel_nm1 = new DistVec<double>(domain->getFaceNormDistInfo());
  } 
  else if (data.typeNormals == DGCLData::IMPLICIT_SECOND_ORDER_EZGCL) {
    edgeNormVel_nm1 = new DistVec<double>(domain->getEdgeDistInfo());
    faceNormVel_nm1 = new DistVec<double>(domain->getFaceNormDistInfo());
  } 
  else if (data.typeNormals == DGCLData::IMPLICIT_THIRD_ORDER_EZGCL) {
    edgeNormVel_nm1 = new DistVec<double>(domain->getEdgeDistInfo());
    edgeNormVel_nm2 = new DistVec<double>(domain->getEdgeDistInfo());
    faceNormVel_nm1 = new DistVec<double>(domain->getFaceNormDistInfo());
    faceNormVel_nm2 = new DistVec<double>(domain->getFaceNormDistInfo());
  } 

  subGeoState = new GeoState*[numLocSub];

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    subGeoState[iSub] = 0;
    double* d2w = d2wall->subData(iSub);
    for (int i=0; i<d2wall->subSize(iSub); ++i){
      d2w[i] += ioData.bc.wall.delta;
  }
  }

}

//------------------------------------------------------------------------------

DistGeoState::DistGeoState(const GeoData& data, Domain *dom, DistInfo& nodeDistInfo, DistInfo& edgeDistInfo,DistInfo& faceNormDistInfo)
: data(data), domain(dom), Xnm1(0), Xnm2(0), Xdot(0), Xsave(0), edgeNorm_nm1(0), edgeNormVel_nm1(0), edgeNorm_nm2(0), edgeNormVel_nm2(0),
    faceNorm_nm1(0), faceNormVel_nm1(0), faceNorm_nm2(0), faceNormVel_nm2(0), ctrlVol_save(0), Xsa(0), dXsa(0), dEdgeNorm(0), dFaceNorm(0),
    dEdgeNormVel(0), dFaceNormVel(0)
{
  numLocSub = domain->getNumLocSub();
  com = domain->getCommunicator();
  oolscale = 1.0;

  Xn = new DistSVec<double,3>(nodeDistInfo);
  ctrlVol_n = new DistVec<double>(nodeDistInfo);
  ctrlVol_nm1 = new DistVec<double>(nodeDistInfo);
  ctrlVol_nm2 = new DistVec<double>(nodeDistInfo);
  // Initialize the values
  *Xn = 0.0;
  *ctrlVol_n = 0.0;
  *ctrlVol_nm1 = 0.0;
  *ctrlVol_nm2 = 0.0;

  d2wall = new DistVec<double>(nodeDistInfo);

  edgeNorm    = new DistVec<Vec3D>(edgeDistInfo);
  edgeNormVel = new DistVec<double>(edgeDistInfo);

  faceNorm    = new DistVec<Vec3D>(faceNormDistInfo);
  faceNormVel = new DistVec<double>(faceNormDistInfo);
  
  inletNodeNorm = new DistVec<Vec3D>(domain->getInletNodeDistInfo());
  numFaceNeighb = new DistVec<int>(domain->getInletNodeDistInfo());

  subGeoState = new GeoState*[numLocSub];

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub)
    subGeoState[iSub] = new GeoState(data, (*ctrlVol_n)(iSub), (*ctrlVol_nm1)(iSub),
                                     (*ctrlVol_nm2)(iSub),
                                     (*d2wall)(iSub),
                                     (*edgeNorm)(iSub), (*faceNorm)(iSub),
                                     (*edgeNormVel)(iSub), (*faceNormVel)(iSub),
                                     (*inletNodeNorm)(iSub), (*numFaceNeighb)(iSub));
}

DistGeoState::DistGeoState(const GeoData& data, Domain *dom, DistInfo& nodeDistInfo, DistInfo& edgeDistInfo)
: data(data), domain(dom), Xnm1(0), Xnm2(0), Xdot(0), Xsave(0), edgeNorm_nm1(0), edgeNormVel_nm1(0), edgeNorm_nm2(0), edgeNormVel_nm2(0),
    faceNorm_nm1(0), faceNormVel_nm1(0), faceNorm_nm2(0), faceNormVel_nm2(0), ctrlVol_save(0), Xsa(0), dXsa(0), dEdgeNorm(0), dFaceNorm(0),
    dEdgeNormVel(0), dFaceNormVel(0)
{
  numLocSub = domain->getNumLocSub();
  com = domain->getCommunicator();
  oolscale = 1.0;

  Xn = new DistSVec<double,3>(nodeDistInfo);
  ctrlVol_n = new DistVec<double>(nodeDistInfo);
  ctrlVol_nm1 = new DistVec<double>(nodeDistInfo);
  ctrlVol_nm2 = new DistVec<double>(nodeDistInfo);
  // Initialize the values
  *Xn = 0.0;
  *ctrlVol_n = 0.0;
  *ctrlVol_nm1 = 0.0;
  *ctrlVol_nm2 = 0.0;

  d2wall = new DistVec<double>(nodeDistInfo);

  edgeNorm    = new DistVec<Vec3D>(edgeDistInfo);
  edgeNormVel = new DistVec<double>(edgeDistInfo);

  faceNorm    = new DistVec<Vec3D>(domain->getFaceDistInfo());
  faceNormVel = new DistVec<double>(domain->getFaceDistInfo());
  
  inletNodeNorm = new DistVec<Vec3D>(domain->getInletNodeDistInfo());
  numFaceNeighb = new DistVec<int>(domain->getInletNodeDistInfo());

  subGeoState = new GeoState*[numLocSub];

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub)
    subGeoState[iSub] = new GeoState(data, (*ctrlVol_n)(iSub), (*ctrlVol_nm1)(iSub),
                                     (*ctrlVol_nm2)(iSub),
                                     (*d2wall)(iSub),
                                     (*edgeNorm)(iSub), (*faceNorm)(iSub),
                                     (*edgeNormVel)(iSub), (*faceNormVel)(iSub),
                                     (*inletNodeNorm)(iSub), (*numFaceNeighb)(iSub));
}

//------------------------------------------------------------------------------

DistGeoState::~DistGeoState()
{

  if (Xn) delete Xn;
  if (Xnm1) delete Xnm1;
  if (Xnm2) delete Xnm2;
  if (Xdot) delete Xdot;
  if (Xsave) delete Xsave;

  if (ctrlVol_n) delete ctrlVol_n;
  if (ctrlVol_nm1) delete ctrlVol_nm1;
  if (ctrlVol_nm2) delete ctrlVol_nm2;
  if (ctrlVol_save) delete ctrlVol_save;

  if (d2wall) delete d2wall;

  if (edgeNorm) delete edgeNorm;
  if (edgeNormVel) delete edgeNormVel;
  if (edgeNorm_nm1) delete edgeNorm_nm1;
  if (edgeNormVel_nm1) delete edgeNormVel_nm1;
  if (edgeNorm_nm2) delete edgeNorm_nm2;
  if (edgeNormVel_nm2) delete edgeNormVel_nm2;

  if (faceNorm) delete faceNorm;
  if (faceNormVel) delete faceNormVel;
  if (faceNorm_nm1) delete faceNorm_nm1;
  if (faceNormVel_nm1) delete faceNormVel_nm1;
  if (faceNorm_nm2) delete faceNorm_nm2;
  if (faceNormVel_nm2) delete faceNormVel_nm2;
  
  if (inletNodeNorm) delete inletNodeNorm;
  if (numFaceNeighb) delete numFaceNeighb;
 

  if (subGeoState) {
#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub)
      if (subGeoState[iSub]) delete subGeoState[iSub];
    
    delete [] subGeoState;
  }

//Included
  if (Xsa) delete Xsa;
  if (dXsa) delete dXsa;
  if (dEdgeNorm) delete dEdgeNorm;
  if (dFaceNorm) delete dFaceNorm;
  if (dEdgeNormVel) delete dEdgeNormVel;

}

//------------------------------------------------------------------------------

void DistGeoState::setup(const char *name, TimeData &timeData,
			 DistSVec<double,3> *X, DistVec<double> *ctrlVol)
{
  setup1(name, X, ctrlVol);
  setup2(timeData);
}

//------------------------------------------------------------------------------

void DistGeoState::setup3(const char *name, DistSVec<double,3> *X, DistVec<double> *ctrlVol)
{


  if (!data.use_n) {
    Xn = X->alias();
    ctrlVol_n = ctrlVol->alias();
  }
  if (!data.use_nm1) {
    Xnm1 = Xn->alias();
    ctrlVol_nm1 = ctrlVol_n->alias();
  }
  if (!data.use_nm2) {
    Xnm2 = Xnm1->alias();
    ctrlVol_nm2 = ctrlVol_nm1->alias();
  }

  bool read_n = false;
  bool read_nm1 = false;
  bool read_nm2 = false;

  // This effectively `restarts' the geometry of the mesh, recovering the
  // last 3 positions.
  if (name[0] != 0) {
    read_n = domain->readVectorFromFile(name, 0, 0, *Xn, &oolscale);
    if (data.use_nm1)
      read_nm1 = domain->readVectorFromFile(name, 1, 0, *Xnm1, &oolscale);
    if (data.use_nm2) 
      read_nm2 = domain->readVectorFromFile(name, 2, 0, *Xnm2, &oolscale);
  }

  if (data.use_nm1)
    *Xnm1 = *Xn;
  if (data.use_nm2)
    *Xnm2 = *Xnm1;

  data.config = 0;

  int ierr = domain->computeControlVolumes(lscale, *Xn, *ctrlVol_n);
#ifdef YDEBUG
  if(ierr) {
    const char* output = "elementvolumecheck";
    ofstream out(output, ios::out);
    if(!out) { cerr << "Error: cannot open file" << output << endl;  exit(-1); }
    out << ierr << endl;
    out.close();
    exit(-1);
  }
#endif
  if (data.use_nm1)
    ierr = domain->computeControlVolumes(lscale, *Xnm1, *ctrlVol_nm1);
  if (data.use_nm2)
    ierr = domain->computeControlVolumes(lscale, *Xnm2, *ctrlVol_nm2);
#ifdef YDEBUG
  if(ierr) {
    const char* output = "elementvolumecheck";
    ofstream out(output, ios::out);
    if(!out) { cerr << "Error: cannot open file" << output << endl;  exit(-1); }
    out << ierr << endl;
    out.close();
    exit(-1);
  }
#endif

	*Xn = *X;
  *ctrlVol = *ctrlVol_n;

  double bbMin[3], bbMax[3];
  X->min(bbMin); X->max(bbMax);

  com->printf(2,
        "Control volume statistics: min=%.3e, max=%.3e, total=%.3e\n"
        "Mesh bounding box: (Xmin,Ymin,Zmin) = (%.3e %.3e %.3e)\n"
        "                   (Xmax,Ymax,Zmax) = (%.3e %.3e %.3e)\n",
        ctrlVol_n->min(), ctrlVol_n->max(), ctrlVol_n->sum(),
              bbMin[0], bbMin[1], bbMin[2], bbMax[0], bbMax[1], bbMax[2]);

}

//------------------------------------------------------------------------------

void DistGeoState::setup1(const char *name, DistSVec<double,3> *X, DistVec<double> *ctrlVol)
{

  if (!data.use_n) {
    Xn = X->alias();
    ctrlVol_n = ctrlVol->alias();
  } 
  if (!data.use_nm1) {
    Xnm1 = Xn->alias();
    ctrlVol_nm1 = ctrlVol_n->alias();
  } 
  if (!data.use_nm2) {
    Xnm2 = Xnm1->alias();
    ctrlVol_nm2 = ctrlVol_nm1->alias();
  } 

  bool read_n = false;
  bool read_nm1 = false;
  bool read_nm2 = false;

  // This effectively `restarts' the geometry of the mesh, recovering the
  // last 3 positions.
  if (name[0] != 0) {
    read_n = domain->readVectorFromFile(name, 0, 0, *Xn, &oolscale);
    if (data.use_nm1)
      read_nm1 = domain->readVectorFromFile(name, 1, 0, *Xnm1, &oolscale);
    if (data.use_nm2) 
      read_nm2 = domain->readVectorFromFile(name, 2, 0, *Xnm2, &oolscale);
  }

  if (!read_n)
    domain->getReferenceMeshPosition(*Xn);
  if (data.use_nm1 && !read_nm1)
    *Xnm1 = *Xn;
  if (data.use_nm2 && !read_nm2)
    *Xnm2 = *Xnm1;

  data.config = 0;
    
  int ierr = domain->computeControlVolumes(lscale, *Xn, *ctrlVol_n);
#ifdef YDEBUG
  if(ierr) {
    const char* output = "elementvolumecheck";
    ofstream out(output, ios::out);
    if(!out) { cerr << "Error: cannot open file" << output << endl;  exit(-1); }
    out << ierr << endl;
    out.close();
    exit(-1);
  }
#endif
  if (data.use_nm1)
    ierr = domain->computeControlVolumes(lscale, *Xnm1, *ctrlVol_nm1);
  if (data.use_nm2)
    ierr = domain->computeControlVolumes(lscale, *Xnm2, *ctrlVol_nm2);
#ifdef YDEBUG
  if(ierr) {
    const char* output = "elementvolumecheck";
    ofstream out(output, ios::out);
    if(!out) { cerr << "Error: cannot open file" << output << endl;  exit(-1); }
    out << ierr << endl;
    out.close();
    exit(-1);
  }
#endif

  *X = *Xn;
  *ctrlVol = *ctrlVol_n;

  double bbMin[3], bbMax[3];
  X->min(bbMin); X->max(bbMax);

  com->printf(2, 
	      "Control volume statistics: min=%.3e, max=%.3e, total=%.3e\n"
	      "Mesh bounding box: (Xmin,Ymin,Zmin) = (%.3e %.3e %.3e)\n"
	      "                   (Xmax,Ymax,Zmax) = (%.3e %.3e %.3e)\n",
	      ctrlVol_n->min(), ctrlVol_n->max(), ctrlVol_n->sum(),
              bbMin[0], bbMin[1], bbMin[2], bbMax[0], bbMax[1], bbMax[2]);

}

//-----------------------------------------------------------------------------

void DistGeoState::setup2(TimeData &timeData)
{

  double oodtnm1 = 1.0 / timeData.dt_nm1;
  double oodtnm2 = 1.0 / timeData.dt_nm2;
  if (data.typeVelocities == DGCLData::IMPLICIT_ZERO ||
      data.typeVelocities == DGCLData::EXPLICIT_RK2_VEL)
    *Xdot = 0.0;
  else
    *Xdot = oodtnm1 * (*Xn - *Xnm1);

// for mesh motion (ALE) with RK2 only
  if (ctrlVol_save)  *ctrlVol_save = *ctrlVol_n;
  if (Xsave)         *Xsave = *Xn;
// end ALE-RK2
  
  /*if(data.typeNormals    == DGCLData::EXPLICIT_RK2 &&
     data.typeVelocities == DGCLData::EXPLICIT_RK2_VEL &&
     data.use_save){
    *edgeNorm = 0.0;
    *edgeNormVel = 0.0;
    *faceNorm = 0.0;
    *faceNormVel = 0.0;

    domain->computeNormalsConfig(*Xn, *Xdot, *edgeNorm, *edgeNormVel, *faceNorm, *faceNormVel);

  }else*/
    domain->computeNormalsGCL1(*Xnm1, *Xn, *Xdot, *edgeNorm, *edgeNormVel, *faceNorm, *faceNormVel);

  domain->computeInletNormals(*inletNodeNorm, *faceNorm, *numFaceNeighb);

  if (data.typeNormals == DGCLData::IMPLICIT_SECOND_ORDER_GCL)
    domain->computeNormalsGCL1(*Xnm1, *Xn, *Xdot, *edgeNorm_nm1, *edgeNormVel_nm1, 
			       *faceNorm_nm1, *faceNormVel_nm1);
  else if (data.typeNormals == DGCLData::IMPLICIT_SECOND_ORDER_EZGCL)
    domain->computeNormalsEZGCL1(oodtnm1, *Xnm1, *Xn, *edgeNorm, *edgeNormVel_nm1, 
				 *faceNorm, *faceNormVel_nm1);
  else if (data.typeNormals == DGCLData::IMPLICIT_THIRD_ORDER_EZGCL) {
    domain->computeNormalsEZGCL1(oodtnm2, *Xnm2, *Xnm1, *edgeNorm, *edgeNormVel_nm2, 
				 *faceNorm, *faceNormVel_nm2);
    domain->computeNormalsEZGCL1(oodtnm1, *Xnm1, *Xn, *edgeNorm, *edgeNormVel_nm1, 
				 *faceNorm, *faceNormVel_nm1);
  }

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub)
    if (!subGeoState[iSub]) {
// Included (MB)
      if (optFlag)
      subGeoState[iSub] = new GeoState(data, (*ctrlVol_n)(iSub), (*ctrlVol_nm1)(iSub),
                                       (*ctrlVol_nm2)(iSub),
                                       (*d2wall)(iSub),
                                       (*edgeNorm)(iSub), (*faceNorm)(iSub),
                                       (*edgeNormVel)(iSub), (*faceNormVel)(iSub),
                                       (*inletNodeNorm)(iSub), (*numFaceNeighb)(iSub),
                                       (*dEdgeNorm)(iSub), (*dFaceNorm)(iSub),
                                       (*dEdgeNormVel)(iSub), (*dFaceNormVel)(iSub),
                                       (*Xsa)(iSub), (*dXsa)(iSub));
      else
      subGeoState[iSub] = new GeoState(data, (*ctrlVol_n)(iSub), (*ctrlVol_nm1)(iSub),
				       (*ctrlVol_nm2)(iSub),
                                       (*d2wall)(iSub), 
				       (*edgeNorm)(iSub), (*faceNorm)(iSub),
				       (*edgeNormVel)(iSub), (*faceNormVel)(iSub),
				       (*inletNodeNorm)(iSub), (*numFaceNeighb)(iSub));
    }
}

//------------------------------------------------------------------------------

void DistGeoState::compute(TimeData &timeData, DistSVec<double,3> &Xsdot, 
			   DistSVec<double,3> &X, DistVec<double> &ctrlVol)
{

  data.config += 1;
    
  //ctrlVol has the control volumes of Xnp1
  int ierr = domain->computeControlVolumes(lscale, X, ctrlVol);
#ifdef YDEBUG
  if(ierr) {
    const char* output = "elementvolumecheck";
    ofstream out(output, ios::out);
    if(!out) { cerr << "Error: cannot open file" << output << endl;  exit(-1); }
    out << ierr << endl;
    out.close();
    exit(-1);
  }
#endif

  //Xdot
  domain->computeVelocities(data.typeVelocities, timeData, Xsdot, *Xnm1, *Xn, X, *Xdot);

  //normals -> edgeNorm, edgeNormVel, faceNorm, faceNormVel
  if (data.typeNormals == DGCLData::IMPLICIT_FIRST_ORDER_GCL || 
      data.typeNormals == DGCLData::IMPLICIT_SECOND_ORDER_GCL) {
    domain->computeNormalsGCL1(*Xn, X, *Xdot, *edgeNorm, *edgeNormVel, 
			       *faceNorm, *faceNormVel);
    if (data.typeNormals == DGCLData::IMPLICIT_SECOND_ORDER_GCL)
      domain->computeNormalsGCL2(timeData, *edgeNorm, *edgeNorm_nm1, *edgeNormVel, 
				 *edgeNormVel_nm1, *faceNorm, *faceNorm_nm1,
				 *faceNormVel, *faceNormVel_nm1);
  }
  else if (data.typeNormals == DGCLData::IMPLICIT_FIRST_ORDER_EZGCL || 
	   data.typeNormals == DGCLData::IMPLICIT_SECOND_ORDER_EZGCL ||
	   data.typeNormals == DGCLData::IMPLICIT_THIRD_ORDER_EZGCL) {
    domain->computeNormalsEZGCL1(1.0/timeData.dt_n, *Xn, X, *edgeNorm, *edgeNormVel, 
				 *faceNorm, *faceNormVel);
    if (data.typeNormals == DGCLData::IMPLICIT_SECOND_ORDER_EZGCL)
      domain->computeNormalsEZGCL2(timeData, *edgeNormVel, *edgeNormVel_nm1, 
				   *faceNormVel, *faceNormVel_nm1);
    else if (data.typeNormals == DGCLData::IMPLICIT_THIRD_ORDER_EZGCL)
      domain->computeNormalsEZGCL3(timeData, *edgeNormVel, *edgeNormVel_nm1, *edgeNormVel_nm2,
				   *faceNormVel, *faceNormVel_nm1, *faceNormVel_nm2);
  }
  else if (data.typeNormals == DGCLData::IMPLICIT_CURRENT_CFG) {
    domain->computeNormalsGCL1(*Xn, *Xn, *Xdot, *edgeNorm, *edgeNormVel, 
			       *faceNorm, *faceNormVel);
  }
  else if (data.typeNormals == DGCLData::IMPLICIT_LATEST_CFG) {
    domain->computeNormalsGCL1(X, X, *Xdot, *edgeNorm, *edgeNormVel, 
			       *faceNorm, *faceNormVel);
  }
  else if (data.typeNormals == DGCLData::EXPLICIT_RK2 && data.use_save){


    // Xsave is temporarily used to compute the different configurations on which to compute the normals
    //       and the normal velocities
    *edgeNorm = 0.0;
    *edgeNormVel = 0.0;
    *faceNorm = 0.0;
    *faceNormVel = 0.0;

    *Xsave = 0.5*(1.0-1.0/sqrt(3.0))*X + 0.5*(1.0+1.0/sqrt(3.0))*(*Xn);
    domain->computeNormalsConfig(*Xsave, *Xdot, *edgeNorm, *edgeNormVel, *faceNorm, *faceNormVel, false);

    *Xsave = 0.5*(1.0+1.0/sqrt(3.0))*X + 0.5*(1.0-1.0/sqrt(3.0))*(*Xn);
    domain->computeNormalsConfig(*Xsave, *Xdot, *edgeNorm, *edgeNormVel, *faceNorm, *faceNormVel, true);

    *edgeNorm *= 0.5;
    *edgeNormVel *= 0.5;
    *faceNorm *= 0.5;
    *faceNormVel *= 0.5;

    //save values of actual configuration Xnp1

    *Xsave = X;
    *ctrlVol_save = ctrlVol;

    // return values X and ctrlVol on which to compute the fluxes!
    // for Eulerian fluxes, implicit/explicit formulation/coding requires normals and normal velocities
    // for Viscous terms, it is always required to have Xconfig which is not necessarily Xnp1!
    
    //X = 0.5*(*Xn) + 0.5*X;
    X += *Xn;
    X *= 0.5;

    // now X has the configuration on which to compute the fluxes 
    // but does not correspond to Xnp1 which will need to be restored
    // at end of iteration!

  }
  else{
    com->fprintf(stderr, "*** Error: Running a mesh moving simulation\n");
    com->fprintf(stderr, "***        but no update of the normals!!!!\n");
    com->fprintf(stderr, "*** Exiting\n");
    exit(1);
  }

  domain->computeInletNormals(*inletNodeNorm, *faceNorm, *numFaceNeighb);

}

//------------------------------------------------------------------------------

// Included (MB)
void DistGeoState::computeDerivatives(DistSVec<double,3> &X,       //mesh position
                                      DistSVec<double,3> &dX,      //derivative of mesh positions
                                      DistSVec<double,3> &Xsdot,   //boundary velocity
                                      DistSVec<double,3> &dXsdot,  //derivative of boundary velocity
                                      DistVec<double> &dCtrlVol)   //derivative of control volumes
{

//Remark: Error mesage for pointers
  if (Xsa == 0) {
    fprintf(stderr, "*** Error: Variable Xsa does not exist!\n");
    exit(1);
  }
  if (dXsa == 0) {
    fprintf(stderr, "*** Error: Variable dXsa does not exist!\n");
    exit(1);
  }
  if (dEdgeNorm == 0) {
    fprintf(stderr, "*** Error: Variable dEdgeNorm does not exist!\n");
    exit(1);
  }
  if (dFaceNorm == 0) {
    fprintf(stderr, "*** Error: Variable dFaceNorm does not exist!\n");
    exit(1);
  }
  if (dEdgeNormVel == 0) {
    fprintf(stderr, "*** Error: Variable dEdgeNormVel does not exist!\n");
    exit(1);
  }
  if (dFaceNormVel == 0) {
    fprintf(stderr, "*** Error: Variable dFaceNormVel does not exist!\n");
    exit(1);
  }

  data.configSA += 1;

  *Xsa=X;
  *dXsa=dX;

  if (data.typeNormals == DGCLData::IMPLICIT_FIRST_ORDER_GCL) {
    domain->computeDerivativeOfNormals(*Xsa, *dXsa, *edgeNorm, *dEdgeNorm, *edgeNormVel, 
                                       *dEdgeNormVel, *faceNorm, *dFaceNorm, *faceNormVel, *dFaceNormVel);
  }
  else {
    fprintf(stderr, "*******************************************************************\n");
    fprintf(stderr, "*** Warning: The normal and the derivative can not be computed, ***\n");
    fprintf(stderr, "*** please check the function type chosen in the class domain!  ***\n");
    fprintf(stderr, "********************************************************************\n");
    fprintf(stderr, "%d %d\n", data.typeNormals, DGCLData::IMPLICIT_FIRST_ORDER_GCL);
    exit(1);
  }

  dCtrlVol = 0.0;
  domain->computeDerivativeOfControlVolumes(lscale, *Xsa, *dXsa, dCtrlVol);

}

//------------------------------------------------------------------------------

// Included (YC)
void DistGeoState::computeDerivatives(RectangularSparseMat<double,3,3> **dEdgeNormdX,
                                      RectangularSparseMat<double,3,3> **dFaceNormdX,
                                      RectangularSparseMat<double,3,1> **dCtrlVoldX,
									  DistSVec<double,3> &X,
									  DistSVec<double,3> &dX,
                                      DistVec<double> &dCtrlVol,
                                      DistVec<Vec3D>& dEdgeNormal,
                                      DistVec<Vec3D>& dFaceNormal, 
                                      DistVec<double>& dFaceNormalVel)
{

//Remark: Error mesage for pointers
  if (Xsa == 0) {
	fprintf(stderr, "*** Error: Variable Xsa does not exist!\n");
	exit(1);
  }
  if (dXsa == 0) {
    fprintf(stderr, "*** Error: Variable dXsa does not exist!\n");
    exit(1);
  }
  if (dEdgeNorm == 0) {
    fprintf(stderr, "*** Error: Variable dEdgeNorm does not exist!\n");
    exit(1);
  }
  if (dFaceNorm == 0) {
    fprintf(stderr, "*** Error: Variable dFaceNorm does not exist!\n");
    exit(1);
  }
  if (dEdgeNormVel == 0) {
    fprintf(stderr, "*** Error: Variable dEdgeNormVel does not exist!\n");
    exit(1);
  }
  if (dFaceNormVel == 0) {
    fprintf(stderr, "*** Error: Variable dFaceNormVel does not exist!\n");
    exit(1);
  }

  data.configSA += 1;

  *Xsa=X;
  *dXsa=dX;

  if (data.typeNormals == DGCLData::IMPLICIT_FIRST_ORDER_GCL) {
    domain->computeDerivativeOfNormals(dEdgeNormdX, dFaceNormdX, *dXsa, dEdgeNormal, *dEdgeNormVel, dFaceNormal, dFaceNormalVel);
    *dEdgeNorm = dEdgeNormal;   *dFaceNorm = dFaceNormal;    *dFaceNormVel = dFaceNormalVel;
    //    dEdgeNormal = *dEdgeNorm;   dFaceNormal = *dFaceNorm;   dFaceNormalVel = *dFaceNormVel;
/*    
    DistSVec<double,3> dX2(dX);
    DistVec<Vec3D> dEdgeNorm2(*dEdgeNorm), dFaceNorm2(*dFaceNorm);
    DistVec<double> dEdgeNormVel2(*dEdgeNormVel), dFaceNormVel2(*dFaceNormVel);
    dX2 = 0.0;
    dEdgeNorm2 = 0.0;     dFaceNorm2 = 0.0;
    dEdgeNormVel2 = 0.0;  dFaceNormVel2 = 0.0;

    domain->computeDerivativeOfNormals(dEdgeNormdX, dFaceNormdX, *dXsa, dEdgeNorm2, dEdgeNormVel2, dFaceNorm2, dFaceNormVel2);
    DistSVec<double,3> dEdgeNorm3(dEdgeNorm2.info()), dFaceNorm3(dFaceNorm2.info());
    DistSVec<double,3> dEdgeNorm4(dEdgeNorm2.info()), dFaceNorm4(dFaceNorm2.info());
//    fprintf(stderr, " .... %d      %d     \n", dEdgeNorm3.info().getMasterFlag(0)[0], dFaceNorm3.info().getMasterFlag(0)[0]);
    for(int iSub=0; iSub < dEdgeNorm2.info().numLocThreads; iSub++) {
      for(int i=0; i<dEdgeNorm2[iSub]->size(); ++i)
        for(int j=0; j<3; ++j) {
          dEdgeNorm3(iSub)[i][j] = dEdgeNorm2(iSub)[i][j]; 
          dEdgeNorm4(iSub)[i][j] = (*dEdgeNorm)(iSub)[i][j]; 
        }
      for(int i=0; i<dFaceNorm2[iSub]->size(); ++i)
        for(int j=0; j<3; ++j) {
          dFaceNorm3(iSub)[i][j] = dFaceNorm2(iSub)[i][j]; 
          dFaceNorm4(iSub)[i][j] = (*dFaceNorm)(iSub)[i][j]; 
        }
//      fprintf(stderr, " norm of dEdgeNorm3 = %e, norm of dEdgeNorm2 = %e\n", dEdgeNorm3(iSub).norm(), dEdgeNorm2(iSub).norm());
      fprintf(stderr, " norm of dFaceNorm3 = %e, norm of dFaceNorm2 = %e\n", dFaceNorm3(iSub).norm(), dFaceNorm2(iSub).norm());
    }

//    fprintf(stderr, " norm of dEdgeNormVel2 is %e\n", dEdgeNormVel2.norm());
//    fprintf(stderr, " norm of dFaceNormVel2 is %e\n", dFaceNormVel2.norm());
    double aa = dEdgeNorm3*dEdgeNorm4 + dFaceNorm3*dFaceNorm4; // + dEdgeNormVel2*(*dEdgeNormVel) + dFaceNormVel2*(*dFaceNormVel);


    domain->computeTransposeDerivativeOfNormals(dEdgeNormdX, dFaceNormdX, *dEdgeNorm, *dFaceNorm, dX2);
    double bb = dX2*dX;
    double diff = sqrt((aa-bb)*(aa-bb));
    if(aa != 0.0) fprintf(stderr, " ... rel. diff = %e\n", diff/abs(aa));
    else fprintf(stderr, " ... abs. diff = %e\n", diff);
*/
  }
  else {
    fprintf(stderr, "*******************************************************************\n");
    fprintf(stderr, "*** Warning: The normal and the derivative can not be computed, ***\n");
    fprintf(stderr, "*** please check the function type chosen in the class domain!  ***\n");
    fprintf(stderr, "********************************************************************\n");
    fprintf(stderr, "%d %d\n", data.typeNormals, DGCLData::IMPLICIT_FIRST_ORDER_GCL);
    exit(1);
  }

  domain->computeDerivativeOfControlVolumes(dCtrlVoldX, *dXsa, dCtrlVol);
/*
  DistSVec<double,3> dX2(*dXsa);
  DistVec<double> dCtrlVol2(dCtrlVol);
  dX2 = 0.0;   dCtrlVol2 = 0.0;

  domain->computeDerivativeOfControlVolumes(dCtrlVoldX, *dXsa, dCtrlVol2);
  double aa = dCtrlVol2*dCtrlVol;

  domain->computeTransposeDerivativeOfControlVolumes(dCtrlVoldX, dCtrlVol, dX2);
  double bb = dX2*(*dXsa);
  double diff = sqrt((aa-bb)*(aa-bb));
  if(aa != 0.0) fprintf(stderr, " ... rel. diff = %e\n", diff/abs(aa));
  else fprintf(stderr, " ... abs. diff = %e\n", diff);
*/
}

//------------------------------------------------------------------------------

// Included (YC)
void DistGeoState::computeTransposeDerivatives(RectangularSparseMat<double,3,3> **dEdgeNormdX,
                                               RectangularSparseMat<double,3,3> **dFaceNormdX,
                                               RectangularSparseMat<double,3,1> **dCtrlVoldX,
                                               DistVec<double>& dCtrlVol,
                                               DistVec<Vec3D>& dNormal,
                                               DistVec<Vec3D>& dn,
                                               DistVec<double>& dndot,
                                               DistSVec<double,3>& dX) 
{

//Remark: Error mesage for pointers
  if (dXsa == 0) {
    fprintf(stderr, "*** Error: Variable dXsa does not exist!\n");
    exit(1);
  }
/*  if (dEdgeNorm == 0) {
    fprintf(stderr, "*** Error: Variable dEdgeNorm does not exist!\n");
    exit(1);
  }
  if (dFaceNorm == 0) {
    fprintf(stderr, "*** Error: Variable dFaceNorm does not exist!\n");
    exit(1);
  }
*/
  data.configSA += 1;

  domain->computeTransposeDerivativeOfControlVolumes(dCtrlVoldX, dCtrlVol, dX);

  if (data.typeNormals == DGCLData::IMPLICIT_FIRST_ORDER_GCL) {
    domain->computeTransposeDerivativeOfNormals(dEdgeNormdX, dFaceNormdX, dNormal, dn, dX);
  }
  else {
    fprintf(stderr, "*******************************************************************\n");
    fprintf(stderr, "*** Warning: The normal and the derivative can not be computed, ***\n");
    fprintf(stderr, "*** please check the function type chosen in the class domain!  ***\n");
    fprintf(stderr, "********************************************************************\n");
    fprintf(stderr, "%d %d\n", data.typeNormals, DGCLData::IMPLICIT_FIRST_ORDER_GCL);
    exit(1);
  }

}

//------------------------------------------------------------------------------

// Included (YC)
void DistGeoState::computeDerivativeOperators(DistSVec<double,3> &X, 
                                              RectangularSparseMat<double,3,3> **dEdgeNormdX,
                                              RectangularSparseMat<double,3,3> **dFaceNormdX,
                                              RectangularSparseMat<double,3,1> **dCtrlVoldX)
{

//Remark: Error mesage for pointers
  if (Xsa == 0) {
    fprintf(stderr, "*** Error: Variable Xsa does not exist!\n");
    exit(1);
  }

  *Xsa=X;

  if (data.typeNormals == DGCLData::IMPLICIT_FIRST_ORDER_GCL) {
    domain->computeDerivativeOperatorsOfNormals(*Xsa, dEdgeNormdX, dFaceNormdX); 
  }
  else {
    fprintf(stderr, "*******************************************************************\n");
    fprintf(stderr, "*** Warning: The normal and the derivative can not be computed, ***\n");
    fprintf(stderr, "*** please check the function type chosen in the class domain!  ***\n");
    fprintf(stderr, "********************************************************************\n");
    fprintf(stderr, "%d %d\n", data.typeNormals, DGCLData::IMPLICIT_FIRST_ORDER_GCL);
    exit(1);
  }

  domain->computeDerivativeOperatorsOfControlVolumes(lscale, *Xsa, dCtrlVoldX);

}

//------------------------------------------------------------------------------

// Included (MB)
void DistGeoState::reset(DistSVec<double,3> &Xmod)
{
  *Xn=Xmod;
  *Xnm1=Xmod;
}

//------------------------------------------------------------------------------

void DistGeoState::interpolate(double dt, double dtLeft, 
			       DistSVec<double,3> &Xs, DistSVec<double,3> &X)
{

  double alpha = dt / (dt + dtLeft) - 1.0;

  X = Xs + alpha * (Xs - *Xn);

}

//------------------------------------------------------------------------------

void DistGeoState::update(DistSVec<double,3> &X, DistVec<double> &ctrlVol)
{
  if (data.use_save){
    // X and ctrlVol had configuration on which fluxes were computed
    // Xnp1 and ctrlVol_np1 are restored to X and ctrlVol
    X = *Xsave;
    ctrlVol = *ctrlVol_save;
  }

  if (data.use_nm2) {
    *Xnm2 = *Xnm1;
    *ctrlVol_nm2 = *ctrlVol_nm1;
  }
  if (data.use_nm1) {
    *Xnm1 = *Xn;
    *ctrlVol_nm1 = *ctrlVol_n;
  }
  if (data.use_n) {
    *Xn = X;
    *ctrlVol_n = ctrlVol;
  }

}

//------------------------------------------------------------------------------

void DistGeoState::writeToDisk(char *name)  {

  if (data.use_n)
    domain->writeVectorToFile(name, 0, 0.0, *Xn, &lscale);
  if (data.use_nm1)
    domain->writeVectorToFile(name, 1, 0.0, *Xnm1, &lscale);
  if (data.use_nm2)  
    domain->writeVectorToFile(name, 2, 0.0, *Xnm2, &lscale);

}

//------------------------------------------------------------------------------

void DistGeoState::setupInitialDisplacement(const char *name, DistSVec<double,3> *X, DistVec<double> *ctrlVol)
{

  if (!data.use_n) {
    Xn = X->alias();
    ctrlVol_n = ctrlVol->alias();
  } 
  if (!data.use_nm1) {
    Xnm1 = Xn->alias();
    ctrlVol_nm1 = ctrlVol_n->alias();
  } 
  if (!data.use_nm2) {
    Xnm2 = Xnm1->alias();
    ctrlVol_nm2 = ctrlVol_nm1->alias();
  } 

  bool read_n = false;
  bool read_nm1 = false;
  bool read_nm2 = false;

  DistSVec<double,3> refPosition(domain->getNodeDistInfo());

  // use relative displacements rather than absolute positions.
  com->printf(2,"Note: Assuming that InitialDisplacement vector uses the same dimensions as the mesh.\n"
                " ... Non-Dimensionalizing InitialDisplacement vector by %e\n", lscale);


  oolscale = (1.0/lscale); // one over lscale
  domain->getReferenceMeshPosition(refPosition);
  read_n = domain->readVectorFromFile(name, 0, 0, *Xn, &oolscale);
  double tmp = Xn->norm();
  double tmp2 = refPosition.norm();
  *Xn = *Xn + refPosition;
  com->printf(2,"... disp = %e, orig = %e, new = %e, diff = %e\n", tmp, tmp2, Xn->norm(), Xn->norm()-tmp2);

  if (data.use_nm1) {
    read_nm1 = domain->readVectorFromFile(name, 1, 0, *Xnm1, &oolscale);
    if (read_nm1) {
      *Xnm1 = *Xnm1+refPosition;
    } else {
      *Xnm1 = *Xn ;
    }
  }
  if (data.use_nm2) {
    read_nm2 = domain->readVectorFromFile(name, 2, 0, *Xnm2, &oolscale);
    if (read_nm2) {
      *Xnm2 = *Xnm2+refPosition;
    } else {
      *Xnm2 = *Xnm1 ;
    }
  } 

  data.config = 0;
    
  int ierr = domain->computeControlVolumes(lscale, *Xn, *ctrlVol_n);
  if (data.use_nm1)
    ierr = domain->computeControlVolumes(lscale, *Xnm1, *ctrlVol_nm1);
  if (data.use_nm2)
    ierr = domain->computeControlVolumes(lscale, *Xnm2, *ctrlVol_nm2);

  *X = *Xn;
  *ctrlVol = *ctrlVol_n;

  double bbMin[3], bbMax[3];
  X->min(bbMin); X->max(bbMax);

  com->printf(2, 
	      "Control volume statistics: min=%.12e, max=%.12e, total=%.12e\n"
	      "Mesh bounding box: (Xmin,Ymin,Zmin) = (%.3e %.3e %.3e)\n"
	      "                   (Xmax,Ymax,Zmax) = (%.3e %.3e %.3e)\n",
	      ctrlVol_n->min(), ctrlVol_n->max(), ctrlVol_n->sum(),
              bbMin[0], bbMin[1], bbMin[2], bbMax[0], bbMax[1], bbMax[2]);

}

//------------------------------------------------------------------------------
