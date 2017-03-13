#include <cstring>
#include <cmath>

#include <MeshMotionHandler.h>

#include <BcDef.h>
#include <VarFcn.h>
#include <MatchNode.h>
#include <StructExc.h>
#include <CorotSolver.h>
#include <MeshMotionSolver.h>
#include <MemoryPool.h>
#include <Timer.h>
#include <Domain.h>

#include <DistVector.h> 
//------------------------------------------------------------------------------

MeshMotionHandler::MeshMotionHandler(IoData &iod, Domain *dom) : 
  domain(dom), F(dom->getNodeDistInfo()), Pin(dom->getFaceDistInfo()), 
  X0(dom->getNodeDistInfo()), dX(dom->getNodeDistInfo()), dXds(dom->getNodeDistInfo()), F0(dom->getNodeDistInfo())
{

  tscale = iod.ref.rv.time;
  oolscale = 1.0 / iod.ref.rv.tlength;
  Wn = iod.restart.energy;
  dtf0 = 0.0;

  Pin = iod.aero.pressure;

  com = dom->getCommunicator();
  domain->getReferenceMeshPosition(X0);

}

//------------------------------------------------------------------------------

RigidMeshMotionHandler::RigidMeshMotionHandler(IoData &ioData, VarFcn *vf, 
					       double *Vin, Domain *dom) :
  Xr(dom->getNodeDistInfo())
{

  typeTag = ioData.rmesh.tag;
  typeLaw = ioData.rmesh.lawtype;
  tscale = ioData.ref.rv.time;

  v[0] = ioData.rmesh.vx;
  v[1] = ioData.rmesh.vy;
  v[2] = ioData.rmesh.vz;
 
  a[0] = ioData.rmesh.ax;
  a[1] = ioData.rmesh.ay;
  a[2] = ioData.rmesh.az;

  vin = vf->getVelocity(Vin);
  cin = vf->computeSoundSpeed(Vin);

  if (typeLaw == RigidMeshMotionData::CONSTANTACCELERATION){
    nvpts = 0;
    timeVpts = 0;
    velVpts  = 0;
  }else
    setupVelocityPoints(ioData);

}

//------------------------------------------------------------------------------
                                                                                                                   
RigidMeshMotionHandler::~RigidMeshMotionHandler()
{
  if(timeVpts) delete timeVpts;
  if(velVpts) delete velVpts;
}

//------------------------------------------------------------------------------

void RigidMeshMotionHandler::setupVelocityPoints(IoData &ioData)
{

  // compute body velocity using Mach number and AoA
  double Pstiff = ioData.eqs.fluidModel.gasModel.pressureConstant;
  double velBody2 =  ioData.eqs.fluidModel.gasModel.specificHeatRatio * (ioData.bc.inlet.pressure + Pstiff)
                  * ioData.bc.inlet.mach*ioData.bc.inlet.mach / ioData.bc.inlet.density;
  double velBody = sqrt(velBody2);

  double velBodyX = velBody * cos(ioData.bc.inlet.alpha) * cos(ioData.bc.inlet.beta);
  double velBodyY = velBody * cos(ioData.bc.inlet.alpha) * sin(ioData.bc.inlet.beta);
  double velBodyZ = velBody * sin(ioData.bc.inlet.alpha);

  // find how many (time-velocity) points have been specified.
  if (ioData.rmesh.vpts[0]->time == 0.0)  {
    fprintf(stderr, " ***ERROR: Velocity Point specified for t=0.  Exiting...\n");
    exit(-1);
  }

  nvpts = 1;
  for (int i = 0; i < RigidMeshMotionData::num; i++)
    if (ioData.rmesh.vpts[i]->time > 0.0)
      nvpts++;
  

  //create the pointers
  timeVpts = new double [nvpts];
  velVpts  = new Vec3D [nvpts];

  // initial time-velocity point is at the freestream
  timeVpts[0] = 0.0;
  velVpts[0][0] = 0.0;
  velVpts[0][1] = 0.0;
  velVpts[0][2] = 0.0;

  // transform the user-specified ground-relative velocities into mesh velocities
  for (int i=1; i < nvpts; i++){
    timeVpts[i] = ioData.rmesh.vpts[i-1]->time;
    velVpts[i][0] = velBodyX - ioData.rmesh.vpts[i-1]->velocityX;
    velVpts[i][1] = velBodyY - ioData.rmesh.vpts[i-1]->velocityY;
    velVpts[i][2] = velBodyZ - ioData.rmesh.vpts[i-1]->velocityZ;
  }

  //we want the points to be ordered with increasing time
  int ierrlocal = 0;
  for (int i=0; i<nvpts-1; i++) {
    if(timeVpts[i] >= timeVpts[i+1])
      ierrlocal++;
  }
  if(ierrlocal){
    fprintf(stderr, "***Error: in under Accelerated, (time-velocity) points are not ordered... Aborting\n");
    exit(1);
  }

}

//------------------------------------------------------------------------------

DistSVec<double,3> &RigidMeshMotionHandler::getFlightPositionVector(double t, DistSVec<double,3> &X)  {

  double dx[3];

  for (int k=0; k<3; ++k)
    dx[k] = vin[k] * t;

  int numLocSub = X.numLocSub();

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {

    double (*x)[3] = X.subData(iSub);
    double (*xr)[3] = Xr.subData(iSub);

    for (int i=0; i<X.subSize(iSub); ++i)
      for (int k=0; k<3; ++k)
        xr[i][k] = x[i][k] - dx[k];

  }

  return Xr;
}

//------------------------------------------------------------------------------

void RigidMeshMotionHandler::addRigidMotion(double t, DistSVec<double,3> &Xrel,
					    DistSVec<double,3> &Xdot, DistSVec<double,3> &X, Communicator *com)
{

  double dx[3], vel[3], dtime;
  int k;
  if(typeLaw == RigidMeshMotionData::CONSTANTACCELERATION)  {
    for (k=0; k<3; ++k) {
      dx[k] = v[k] * t + 0.5 * a[k] * t*t;
      vel[k] = v[k] + a[k] * t;
    }
  }
  else if (typeLaw == RigidMeshMotionData::VELOCITYPOINTS)  {
    // find the current timeInterval: intervals are labeled from 0 to nvpts-1

    int timeInt = nvpts-1;
    for (int i=0; i<nvpts-1; i++)
      if (t >= timeVpts[i] && t < timeVpts[i+1])  {
        timeInt = i;
        break;
      }

    //velocity = interpolation within interval and last prescribed value if out of bound
    //position = integration of velocity on [0,t] passing through all intervals
    for (k=0; k<3; ++k) {
      if (timeInt < nvpts-1)  {
        vel[k] = (timeVpts[timeInt+1]-t)*velVpts[timeInt][k]+(t-timeVpts[timeInt])*velVpts[timeInt+1][k];
        vel[k] /= (timeVpts[timeInt+1]-timeVpts[timeInt]);
      }
      else
        vel[k] = velVpts[timeInt][k];

      dx[k] = 0.0;
      for (int nInt = 0; nInt < timeInt; nInt++)
        dx[k] += (velVpts[nInt][k] + velVpts[nInt+1][k]) * (timeVpts[nInt+1] - timeVpts[nInt]);
       
      dx[k] += (vel[k] + velVpts[timeInt][k]) * (t - timeVpts[timeInt]);
      dx[k] *= 0.5;

    } 
  }
  if (com)
    com->fprintf(stderr, "Body pos/vel at time %f: %f / %f with acc = %f\n", t, dx[0], vel[0], a[0]);

  int numLocSub = X.numLocSub();

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {

    double (*xrel)[3] = Xrel.subData(iSub);
    double (*xdot)[3] = Xdot.subData(iSub);
    double (*x)[3] = X.subData(iSub);

    for (int i=0; i<X.subSize(iSub); ++i) {
      for (k=0; k<3; ++k) {
	x[i][k] = xrel[i][k] + dx[k];
	xdot[i][k] += vel[k];
      }
    }

  }

  deltaX[0] = dx[0];
  deltaX[1] = dx[1];
  deltaX[2] = dx[2];

}

//------------------------------------------------------------------------------

void RigidMeshMotionHandler::updateMomentArm(Vec3D &x0)  {

  x0 += deltaX;

}

//------------------------------------------------------------------------------

DistSVec<double,3> &RigidMeshMotionHandler::getRelativePositionVector(double t, DistSVec<double,3> &X)
{

/*
  double dx[3];

  for (int k=0; k<3; ++k)
    dx[k] = v[k] * t + 0.5 * a[k] * t*t;
*/
  int numLocSub = X.numLocSub();

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {

    double (*x)[3] = X.subData(iSub);
    double (*xr)[3] = Xr.subData(iSub);

    for (int i=0; i<X.subSize(iSub); ++i)
      for (int k=0; k<3; ++k)
	xr[i][k] = x[i][k] - deltaX[k];

  }

  return Xr;

}

//------------------------------------------------------------------------------

const char* RigidMeshMotionHandler::getTagName()
{

  const char* tag = "None";

  if (typeTag == RigidMeshMotionData::MACH)
    tag = "Mach";
  else if(typeTag == RigidMeshMotionData::TIME)
    tag = "Time";
  else
    tag = "Velocity";

  return tag;

}

//------------------------------------------------------------------------------

double RigidMeshMotionHandler::getTagValue(double t)
{

  double tag;

  if (typeTag == RigidMeshMotionData::MACH) {
    Vec3D vnew = vin - (v + t * a);
    tag = sqrt(vnew*vnew) / cin;
  } else if(typeTag == RigidMeshMotionData::TIME)
    tag = t * tscale;
  else {
    Vec3D vel = vin - (v + t * a);
    tag = sqrt(vel*vel);
  }
    

  return tag;

}

//------------------------------------------------------------------------------

AccMeshMotionHandler::AccMeshMotionHandler(IoData &iod, VarFcn *vf, 
					   double *Vin, Domain *dom) :
  MeshMotionHandler(iod, dom), RigidMeshMotionHandler(iod, vf, Vin, dom)
{

  dt = iod.rmesh.timestep;
}

//------------------------------------------------------------------------------

double AccMeshMotionHandler::update(bool *lastIt, int it, double t,
                                    DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{

  if (*lastIt) return dt;

  Xdot = 0.0;

  addRigidMotion(t + dt, X0, Xdot, X, com);

  return dt;

}
//------------------------------------------------------------------------------

double AccMeshMotionHandler::updateStep2(bool *lastIt, int it, double t,
                                         DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{

  if (*lastIt) return dt;

  Xdot = 0.0;

  addRigidMotion(t + dt, X0, Xdot, X, com);

  return dt;

}

//------------------------------------------------------------------------------

double AccMeshMotionHandler::updateStep1(bool *lastIt, int it, double t,
                                         DistSVec<double,3> &Xdot, DistSVec<double,3> &X, double *tmax)
{

 return dt;
}

//------------------------------------------------------------------------------

AeroMeshMotionHandler::AeroMeshMotionHandler(IoData &ioData, VarFcn *varFcn, 
                       double *Vin, MatchNodeSet **matchNodes, 
		       Domain *dom, MemoryPool *mp) :
                       MeshMotionHandler(ioData, dom) {

  steady = !ioData.problem.type[ProblemData::UNSTEADY];
  forceComputation = ioData.aero.force;
  if (ioData.ts.type == TsData::IMPLICIT) {
    if (ioData.ts.implicit.type == ImplicitData::BACKWARD_EULER)
      timeIntegrator = IMPLICIT_FIRST_ORDER;
    else if (ioData.ts.implicit.type == ImplicitData::CRANK_NICOLSON ||
	     ioData.ts.implicit.type == ImplicitData::THREE_POINT_BDF)
      timeIntegrator = IMPLICIT_SECOND_ORDER;
  }
  else if (ioData.ts.type == TsData::EXPLICIT)
    timeIntegrator = IMPLICIT_SECOND_ORDER;

  Fn = 0;
  Fnp1 = 0;
  Favg = 0;
  if (forceComputation == AeroelasticData::AVERAGED) {
    Favg = new DistSVec<double,3>(dom->getNodeDistInfo());
    if (timeIntegrator == IMPLICIT_SECOND_ORDER) {
      Fn = new DistSVec<double,3>(dom->getNodeDistInfo());
      Fnp1 = new DistSVec<double,3>(dom->getNodeDistInfo());
    }
  }
  else if (forceComputation == AeroelasticData::LAST_KRIS) {
    Fn = new DistSVec<double,3>(dom->getNodeDistInfo());
    *Fn = 0.0;
  }

  strExc = new StructExc(ioData, matchNodes, 6, domain->getStrCommunicator(), domain->getCommunicator(),
           domain->getNumLocSub());
  strExc->negotiate();
  mppFactor = strExc->getInfo();

  //com->fprintf(stderr, " ... Starting Struct Solver\n");

  mms = new TetMeshMotionSolver(ioData.dmesh, matchNodes, domain, mp);

  it0 = ioData.restart.iteration;

  int sp = strlen(ioData.output.restart.prefix) + 1;
  posFile = new char[sp + strlen(ioData.output.restart.positions)];
  if (ioData.output.restart.positions[0] != 0)
    sprintf(posFile, "%s%s", ioData.output.restart.prefix, ioData.output.restart.positions);
  else
    sprintf(posFile, "");

}

//------------------------------------------------------------------------------

void AeroMeshMotionHandler::negotiate()
{
  strExc->negotiate();
}

//------------------------------------------------------------------------------

AeroMeshMotionHandler::~AeroMeshMotionHandler()
{

  if (Fn) delete Fn;
  if (Fnp1) delete Fnp1;
  if (Favg) delete Favg;
  if (strExc) delete strExc;
  if (mms) delete mms;
}

//------------------------------------------------------------------------------

void AeroMeshMotionHandler::sendForceSensitivity(DistSVec<double,3> *dFdS, bool applyScale)
{
	strExc->sendForce(*dFdS, applyScale);
}

//------------------------------------------------------------------------------

void AeroMeshMotionHandler::getNumParam(int &numParam, int &actvar, double &steadyTol)
{
  strExc->getNumParam(numParam, actvar, steadyTol);
}

//------------------------------------------------------------------------------

void AeroMeshMotionHandler::sendNumParam(int numParam)
{
  strExc->sendNumParam(numParam);
}

//------------------------------------------------------------------------------

void AeroMeshMotionHandler::getRelResidual(double &relres)
{
  strExc->getRelResidual(relres);
}

//------------------------------------------------------------------------------

void AeroMeshMotionHandler::cmdCom(bool *lastIt)
{
  strExc->negotiateStopping(lastIt); 
}

//------------------------------------------------------------------------------
/*
  X0 = reference configuration
  X = current configuration
  dX = relative displacement (i.e. with respect to X) of the boundaries
*/

double AeroMeshMotionHandler::update(bool *lastIt, int it, double t,
                                     DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{
  int algNum = strExc->getAlgorithmNumber();
  double dt = strExc->getTimeStep();
  if (steady)
    dt = 0.0;

  if (algNum == 6 && it == 0) 
    dt *= 0.5;
  if ((algNum == 20 || algNum == 21 || algNum == 22) && it == 0) 
    dt *= 0.5;

  if (algNum == 8) {
    getModalMotion(X);
    *lastIt = true;
  }
  else  {
    if (algNum == 4) {
      strExc->getDisplacement(X0, X, Xdot, dX, false);
      if (*lastIt) 
        return 0.0;
      strExc->sendForce(F);
    }
    else {
      if (it > it0 && algNum != 10) {
        strExc->sendForce(F);
        if (steady) {
    	  strExc->negotiateStopping(lastIt);
        if (*lastIt) 
          return 0.0;
        }
      }
      strExc->getDisplacement(X0, X, Xdot, dX, false);
      if (*lastIt) 
        return 0.0;
    }

   mms->applyProjector(Xdot); //HB: make sure Xdot satisfies the sliding conditions
   mms->solve(dX, X); //HB: the sliding conditions are also applied to dX inside the solve method

   if (algNum == 1) //Ping-Pong 
     *lastIt = true;
  }

  return dt;

}

//------------------------------------------------------------------------------
/*
  X0 = reference configuration
  X = current configuration
  dX = relative displacement (i.e. with respect to X) of the boundaries
*/

double AeroMeshMotionHandler::updateStep1(bool *lastIt, int it, double t,
                                          DistSVec<double,3> &Xdot, DistSVec<double,3> &X, double *tmax)
{
  Timer *timer;
  timer = domain->getTimer();
  double t0 = timer->getTime();

  int algNum = strExc->getAlgorithmNumber();
  double dt = strExc->getTimeStep();
  if (steady)
    dt = 0.0;

  if (algNum == 6 && it == 0)
    dt *= 0.5;
  if ((algNum == 20 || algNum == 21 || algNum == 22) && it == 0)
    dt *= 0.5;

  if (algNum == 20) { // RK2-CD algorithm with Aero-S C0_legacy
    if(it==0) {strExc->getDisplacement(X0,X,Xdot,dX, false);} //for proper restart
    else if(it==it0) {/*nothing to do*/}
    else if(it!=1) {strExc->sendForce(F);}
  }
  else if (algNum == 21) { // RK2-CD algorithm with XED3D
    if(it==0) {strExc->getDisplacement(X0,X,Xdot,dX, false);} //for proper restart
    else if(it==it0) {strExc->sendForce(F);}
    else if(it!=1) {strExc->sendForce(F);}
  }
  else if (algNum == 22) { // RK2-CD algorithm with Dyna3D-XFEM or Aero-S C0
    if(it==0) strExc->getDisplacement(X0,X,Xdot,dX, false); //for proper restart
    if(it!=1) strExc->sendForce(F);
    if(it==0 || it!=(*lastIt)) {strExc->getInfo(); dt = strExc->getTimeStep(); if(it==0) dt *= 0.5;}
    *tmax = strExc->getMaxTime();
  }
  else if (algNum == 8) {
    getModalMotion(X);
    *lastIt = true;
  }
  else if (algNum == 4 || algNum == 5) {
    if (*lastIt)
      return 0.0;
    strExc->getDisplacement(X0, X, Xdot, dX, false);
  }
  else if (algNum == 10 && (/*it == 0 ||*/ *lastIt) ) // change to align B0 with manual description
    strExc->sendForce(F);
  else if (it > it0 && algNum != 10) {
    strExc->sendForce(F);
  }
  timer->removeForceAndDispComm(t0); // do not count the communication time with the
                                     // structure in the mesh solution

  //com->fprintf(stderr, "Aero F sent Force norm = %e\n", F.norm());
  return dt;
}

//------------------------------------------------------------------------------
/*
  X0 = reference configuration
  X = current configuration
  dX = relative displacement (i.e. with respect to X) of the boundaries
*/

double AeroMeshMotionHandler::updateStep2(bool *lastIt, int it, double t,
                                          DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{
  int verbose = com->getMaxVerbose();
  Timer *timer;
  timer = domain->getTimer();
  double t0 = timer->getTime();

  int algNum = strExc->getAlgorithmNumber();
  double dt = strExc->getTimeStep();


  if (algNum == 6 && it == 0) {
    dt *= 0.5;
  }
  if ((algNum == 20 || algNum == 21 || algNum == 22) && it == 0)
    dt *= 0.5;

  if (algNum == 20 || algNum == 21) {
    if(it==0) {strExc->sendForce(F);}
    else if(!*lastIt) {strExc->getDisplacement(X0, X, Xdot, dX, false);}
    else return 0.0; // last iteration!
  }
  else if (algNum == 22) {
    if(it==0) {/*nothing to do*/}
    else if(!*lastIt) {strExc->getDisplacement(X0, X, Xdot, dX, false);}
    else return 0.0; // last iteration!
  }
  else if (algNum == 8) {
    return dt;
  }


  if (algNum == 4 || algNum == 5) {
    if (*lastIt)
      return 0.0;
    strExc->sendForce(F);
  }
  else if(!(algNum == 20 || algNum == 21 || algNum == 22)) {
    if (it > it0 && algNum != 10) {
      if (steady) {
        strExc->negotiateStopping(lastIt);  // [F] receive iSteady from structure
        if(verbose>4) com->fprintf(stderr," ... [F] received iSteady from structure\n");
        if (*lastIt) {
          return 0.0;
        }
      }
    }
    if (algNum != 10 && !*lastIt)
      strExc->getDisplacement(X0, X, Xdot, dX, false);  // [F] receive displacement from structure
    else
      if (it == it0)
        strExc->getDisplacement(X0, X, Xdot, dX, false);
      
  }
  timer->removeForceAndDispComm(t0); // do not count the communication time with the
                                     // structure in the mesh solution

  if ((algNum != 10 && !*lastIt) || it == it0)  {

    com->fprintf(stdout, "... It %5d: Received Incr. Disp. and Vel. ==> %e and %e (%e) \n", it, dX.norm(), Xdot.norm(), X0.norm());

    mms->applyProjector(Xdot); //HB: make sure Xdot satisfies the sliding conditions
    com->sync();
    mms->solve(dX, X); //HB: the sliding conditions are also applied to dX inside the solve method
    com->sync();
  }

  if (algNum == 1)
    *lastIt = true;

  return dt;
}

//------------------------------------------------------------------------------

void AeroMeshMotionHandler::setPositionVector(DistSVec<double,3> &X)
{
  if(mms1)
    mms1->setup(X);

}

//------------------------------------------------------------------------------
/*
  dXdSb = boundary displacement sensitivity
*/

void AeroMeshMotionHandler::updateDStep2(DistSVec<double,3> &X, DistSVec<double,3> &dXdSb, bool applyScale)
{
  int verbose = com->getMaxVerbose();
  Timer *timer;
  timer = domain->getTimer();
  double t0 = timer->getTime();

  int algNum = strExc->getAlgorithmNumber();

  strExc->getDisplacementSensitivity(X,dXdSb);
}

//------------------------------------------------------------------------------
int AeroMeshMotionHandler::getAlgNum()  {

  return strExc->getAlgorithmNumber();

}

//------------------------------------------------------------------------------

int AeroMeshMotionHandler::getModalMotion(DistSVec<double,3> &X)
{

  int nf, ierr(0);
  double *f;
  strExc->getMdFreq(nf, f);
  com->fprintf(stderr, " ... Reading %d modal frequencies\n", nf);
  com->fprintf(stderr, " ... Writing Fluid mesh modes to %s\n", posFile);
  com->fprintf(stderr, " ... Using MppFactor: %f\n", mppFactor);
  domain->writeVectorToFile(posFile, 0, nf, X);
  DistVec<double> cv(domain->getNodeDistInfo());
  double lscale = 1.0/oolscale;

  for (int im = 0; im < nf; ++im) {

    double rscale=1;
    X = X0;
    strExc->getMdStrDisp(im, X0, X, dX);
    for(int i=0; i<4; ++i) {
      if(i > 0) rscale /= 2;
      X = X0;
      dX *= rscale;    
      mms->solve(dX, X);

      // verify mesh integrity
      ierr = domain->computeControlVolumes(lscale, X, cv);
      if(ierr > 0) {
        com->fprintf(stderr, " ... WARNING: negative volume element is detected. turn modal displacement upside down\n");
        dX *= -1;
        X = X0;
        mms->solve(dX, X);
        ierr = domain->computeControlVolumes(lscale, X, cv);
      }

      if(ierr == 0) {
        dX = X - X0;
        dX *= (mppFactor/rscale);
        domain->writeVectorToFile(posFile, im+1, f[im], dX);
        break;
      } else com->fprintf(stderr, " ... WARNING: negative volume element is detected. attempt to increase mppFactor to %e\n", mppFactor/(rscale/2));
    }

#ifdef YDEBUG
    if(ierr > 0) {
      const char* output = "elementvolumecheck";
      ofstream out(output, ios::out);
      if(!out) { cerr << "Error: cannot open file" << output << endl;  exit(-1); }
      out << ierr << endl;
      out.close();
      exit(-1);
    }
#endif

  }
  return 1;

}

//------------------------------------------------------------------------------

AccAeroMeshMotionHandler::AccAeroMeshMotionHandler(IoData &iod, VarFcn *vf, double *Vin, 
						   MatchNodeSet **mns, Domain *dom, MemoryPool *mp) :
  AeroMeshMotionHandler(iod, vf, Vin, mns, dom, mp), 
  RigidMeshMotionHandler(iod, vf, Vin, dom)
{

}

//------------------------------------------------------------------------------

double AccAeroMeshMotionHandler::update(bool *lastIt, int it, double t,
                                        DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{

  DistSVec<double,3> &Xrel = getRelativePositionVector(t, X);

  double dt = AeroMeshMotionHandler::update(lastIt, it, t, Xdot, Xrel);

  if (*lastIt) return dt;

  addRigidMotion(t + dt, Xrel, Xdot, X, com);

  return dt;

}

//------------------------------------------------------------------------------

double AccAeroMeshMotionHandler::updateStep1(bool *lastIt, int it, double t,
                                             DistSVec<double,3> &Xdot, DistSVec<double,3> &X, double *tmax)
{

  DistSVec<double,3> &Xrel = getRelativePositionVector(t, X);

  double dt = AeroMeshMotionHandler::updateStep1(lastIt, it, t, Xdot, Xrel);

  return dt;

}

//------------------------------------------------------------------------------

double AccAeroMeshMotionHandler::updateStep2(bool *lastIt, int it, double t,
                                             DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{

  DistSVec<double,3> &Xrel = getRelativePositionVector(t, X);

  double dt = AeroMeshMotionHandler::updateStep2(lastIt, it, t, Xdot, Xrel);

  if (*lastIt) return dt;

  addRigidMotion(t + dt, Xrel, Xdot, X, com);

  return dt;

}

//------------------------------------------------------------------------------

DeformingMeshMotionHandler::DeformingMeshMotionHandler(IoData &iod, Domain *dom) :
                            MeshMotionHandler(iod, dom), dXmax(0)
{

  dt = iod.forced.timestep;
  omega = 2.0 * acos(-1.0) * iod.forced.frequency;

  if (iod.forced.df.positions[0] != 0) {
    dXmax = new DistSVec<double,3>(dom->getNodeDistInfo());
    domain->readVectorFromFile(iod.forced.df.positions, 0, 0, *dXmax, &oolscale);
    *dXmax = iod.forced.df.amplification * (*dXmax - X0);
  }

  if (iod.forced.df.domain == DeformingData::VOLUME)
    mms = 0;
  else if (iod.forced.df.domain == DeformingData::SURFACE)  {
    mms = new TetMeshMotionSolver(iod.dmesh, 0, domain, 0); //HB: as we passed not MatchNodeSet, it will create the dofType array
                                                            //(see SubDomain::getMeshMotionDofType) considering the nodes labelled
                                                            //as moving as "matched" nodes (i.e. nodes with given (non-zero) displacement)
                                                            //Another option, would have been to create an MatchNodeSet using the "moving" nodes
  }


}

//------------------------------------------------------------------------------

DeformingMeshMotionHandler::~DeformingMeshMotionHandler()
{

  if (dXmax) delete dXmax;
  if (mms) delete mms;

}

//------------------------------------------------------------------------------

void DeformingMeshMotionHandler::setup(DistSVec<double,3> &X)
{

  if(mms) mms->setup(X);

}

//------------------------------------------------------------------------------

DistSVec<double,3> DeformingMeshMotionHandler::getModes()
{

  return(*dXmax);

}

//------------------------------------------------------------------------------

double DeformingMeshMotionHandler::update(bool *lastIt, int it, double t,
                                          DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{

  Timer *timer;
  timer = domain->getTimer();


  if (*lastIt) return dt;

  if (dXmax) {
    dX = sin(omega * (t + dt)) * (*dXmax) - (X - X0);
    Xdot = omega * cos(omega * (t + dt)) * (*dXmax);
  }
  else
  {
    static const double pi = acos(-1.0);
    double dxmax = 0.05;
    double dx = -cos(2.0*pi*(t + dt)) * (t + dt)*(t + dt)*(t + dt) * 0.125*dxmax;
    double dy = -sin(2.0*pi*(t + dt)) * (t + dt)*(t + dt)*(t + dt) * 0.125*dxmax;

    com->printf(5, "dx=%e and dy=% at time = %e\n", dx, dy, (t + dt) * tscale);

    int numLocSub = domain->getNumLocSub();
    #pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub)
    {
      double (*x)[3] = X.subData(iSub);
      for (int i=0; i<X.subSize(iSub); ++i)
      {
        x[i][0] += dx * cos(pi*x[i][0]);
        x[i][1] += dy * cos(pi*x[i][1]);
      }
    }
  }

  if (mms)
  { //HB: changed to use dofType array instead of nodeType array
    int numLocSub = domain->getNumLocSub();
    BCApplier* meshMotionBCs = domain->getMeshMotionBCs();
    int** DofType = (meshMotionBCs) ? meshMotionBCs->getDofType(): 0;

    if(DofType)
    {
      #pragma omp parallel for
      for(int iSub=0; iSub<numLocSub; ++iSub)
      {
        double (*dx)[3] = dX.subData(iSub);
        int (*dofType)[3] = reinterpret_cast<int (*)[3]>(DofType[iSub]);
        for(int i=0;i<dX.subSize(iSub); i++)
          for(int l=0; l<3; l++)
            if(dofType[i][l]!=BC_MATCHED) dx[i][l] = 0.0;
      }
    }

    mms->applyProjector(Xdot); //HB: make sure Xdot satisfies the sliding conditions
    mms->solve(dX, X); //HB: the sliding conditions are also applied to dX inside the solve method
		com->fprintf(stderr,"at iteration %d, mesh is updated. ", it);

  }
  else {
    X += dX;
  }	

  return dt;

}

//------------------------------------------------------------------------------

double DeformingMeshMotionHandler::updateStep1(bool *lastIt, int it, double t,
                                               DistSVec<double,3> &Xdot, DistSVec<double,3> &X, double *tmax)
{

  return dt;

}

//------------------------------------------------------------------------------

double DeformingMeshMotionHandler::updateStep2(bool *lastIt, int it, double t,
                                               DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{

 return update(lastIt, it, t, Xdot, X);

}

//------------------------------------------------------------------------------

PitchingMeshMotionHandler::PitchingMeshMotionHandler(IoData &iod, Domain *dom) :
                            MeshMotionHandler(iod, dom)
{

  dt = iod.forced.timestep;
  omega = 2.0 * acos(-1.0) * iod.forced.frequency;

  alpha_in    = (acos(-1.0)*iod.forced.pt.alpha_in) / 180.0;    // initial angle of rotation
  alpha_max   = (acos(-1.0)*iod.forced.pt.alpha_max) / 180.0;   // maximum angle of rotation
  alpha_slope = (acos(-1.0)*iod.forced.pt.alpha_slope) / 180.0; // (for a linear pitching law)

  beta_in    = (acos(-1.0)*iod.forced.pt.beta_in) / 180.0;    // initial angle of second rotation
  beta_max   = (acos(-1.0)*iod.forced.pt.beta_max) / 180.0;   // maximum angle of second rotation
  beta_slope = (acos(-1.0)*iod.forced.pt.beta_slope) / 180.0; // (for a linear pitching law)

  x1[0] = iod.forced.pt.x11;
  x1[1] = iod.forced.pt.y11;
  x1[2] = iod.forced.pt.z11;

  x2[0] = iod.forced.pt.x21;
  x2[1] = iod.forced.pt.y21;
  x2[2] = iod.forced.pt.z21;

  y1[0] = iod.forced.pt.x12;
  y1[1] = iod.forced.pt.y12;
  y1[2] = iod.forced.pt.z12;

  y2[0] = iod.forced.pt.x22;
  y2[1] = iod.forced.pt.y22;
  y2[2] = iod.forced.pt.z22;

  u = x2[0]-x1[0];
  v = x2[1]-x1[1];
  w = x2[2]-x1[2];

  // unit normals of axis of rotation //

  ix = u/sqrt(u*u+v*v+w*w);
  iy = v/sqrt(u*u+v*v+w*w);
  iz = w/sqrt(u*u+v*v+w*w);

  if (iod.forced.pt.domain == PitchingData::VOLUME)
    mms = 0;
  else if (iod.forced.pt.domain == PitchingData::SURFACE)
    mms = new TetMeshMotionSolver(iod.dmesh, 0, domain, 0);

}

//------------------------------------------------------------------------------

PitchingMeshMotionHandler::~PitchingMeshMotionHandler()
{

  if (mms) delete mms;

}

//------------------------------------------------------------------------------

void PitchingMeshMotionHandler::setup(DistSVec<double,3> &X)
{

  if(mms) mms->setup(X);

}

//------------------------------------------------------------------------------

DistSVec<double,3> PitchingMeshMotionHandler::getModes()
{

  int numLocSub = domain->getNumLocSub();

  double theta = alpha_in + alpha_max;
  double costheta = cos(theta);
  double sintheta = sin(theta);

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {

    double (*dx)[3] = dX.subData(iSub);
    double (*x0)[3] = X0.subData(iSub);

    for (int i=0; i<dX.subSize(iSub); ++i) {

      dx[i][0] = 0.0; dx[i][1] = 0.0; dx[i][2] = 0.0;

      double p[3];
      p[0] = x0[i][0] -  x1[0];
      p[1] = x0[i][1] -  x1[1];
      p[2] = x0[i][2] -  x1[2];

      dx[i][0] += (costheta + (1 - costheta) * ix * ix) * p[0];
      dx[i][0] += ((1 - costheta) * ix * iy - iz * sintheta) * p[1];
      dx[i][0] += ((1 - costheta) * ix * iz + iy * sintheta) * p[2];

      dx[i][1] += ((1 - costheta) * ix * iy + iz * sintheta) * p[0];
      dx[i][1] += (costheta + (1 - costheta) * iy * iy) * p[1];
      dx[i][1] += ((1 - costheta) * iy * iz - ix * sintheta) * p[2];

      dx[i][2] += ((1 - costheta) * ix * iz - iy * sintheta) * p[0];
      dx[i][2] += ((1 - costheta) * iy * iz + ix * sintheta) * p[1];
      dx[i][2] += (costheta + (1 - costheta) * iz * iz) * p[2];

      dx[i][0] += x1[0];
      dx[i][1] += x1[1];
      dx[i][2] += x1[2];

    }
  }

  dX -= X0;

  return(dX);

}

//------------------------------------------------------------------------------

double PitchingMeshMotionHandler::update(bool *lastIt, int it, double t,
                                         DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{

  if (*lastIt) return dt;

  double theta = (omega==0.0 && alpha_slope!=0.0) ? alpha_in + alpha_slope * (t + dt) : alpha_in + alpha_max * sin(omega * (t + dt));
  double costheta = cos(theta);
  double sintheta = sin(theta);

  double phi = (omega==0.0 && beta_slope!=0.0) ? beta_in + beta_slope * (t + dt) : beta_in + beta_max * sin(omega * (t + dt)); 
  double cosphi = cos(phi);
  double sinphi = sin(phi);

  int numLocSub = domain->getNumLocSub();

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

  jx = u/sqrt(u*u+v*v+w*w);
  jy = v/sqrt(u*u+v*v+w*w);
  jz = w/sqrt(u*u+v*v+w*w);

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {

    double (*dx)[3] = dX.subData(iSub);
    double (*x0)[3] = X0.subData(iSub);

    for (int i=0; i<dX.subSize(iSub); ++i) {

      p[0] = x0[i][0] - x1[0];
      p[1] = x0[i][1] - x1[1];
      p[2] = x0[i][2] - x1[2];

      dx[i][0] = 0.0; dx[i][1] = 0.0; dx[i][2] = 0.0;
      dx[i][0] += (costheta + (1 - costheta) * ix * ix) * p[0];
      dx[i][0] += ((1 - costheta) * ix * iy - iz * sintheta) * p[1];
      dx[i][0] += ((1 - costheta) * ix * iz + iy * sintheta) * p[2];

      dx[i][1] += ((1 - costheta) * ix * iy + iz * sintheta) * p[0];
      dx[i][1] += (costheta + (1 - costheta) * iy * iy) * p[1];
      dx[i][1] += ((1 - costheta) * iy * iz - ix * sintheta) * p[2];

      dx[i][2] += ((1 - costheta) * ix * iz - iy * sintheta) * p[0];
      dx[i][2] += ((1 - costheta) * iy * iz + ix * sintheta) * p[1];
      dx[i][2] += (costheta + (1 - costheta) * iz * iz) * p[2];

      dx[i][0] += x1[0];
      dx[i][1] += x1[1];
      dx[i][2] += x1[2];

      p[0] = dx[i][0] - yy1[0];                               
      p[1] = dx[i][1] - yy1[1];                             
      p[2] = dx[i][2] - yy1[2];                         

      dx[i][0] = 0.0; dx[i][1] = 0.0; dx[i][2] = 0.0;              
      dx[i][0] += (cosphi + (1 - cosphi) * jx * jx) * p[0];            
      dx[i][0] += ((1 - cosphi) * jx * jy - jz * sinphi) * p[1];      
      dx[i][0] += ((1 - cosphi) * jx * jz + jy * sinphi) * p[2];    

      dx[i][1] += ((1 - cosphi) * jx * jy + jz * sinphi) * p[0];       
      dx[i][1] += (cosphi + (1 - cosphi) * jy * jy) * p[1];           
      dx[i][1] += ((1 - cosphi) * jy * jz - jx * sinphi) * p[2];     

      dx[i][2] += ((1 - cosphi) * jx * jz - jy * sinphi) * p[0];   
      dx[i][2] += ((1 - cosphi) * jy * jz + jx * sinphi) * p[1];  
      dx[i][2] += (cosphi + (1 - cosphi) * jz * jz) * p[2];           

      dx[i][0] += yy1[0];                                            
      dx[i][1] += yy1[1];                                        
      dx[i][2] += yy1[2];                                     
    }
  }

  dX -= X;

  if (mms)
  { //HB: changed to use dofType array instead of nodeType array
    int numLocSub = domain->getNumLocSub();
    BCApplier* meshMotionBCs = domain->getMeshMotionBCs();
    int** DofType = (meshMotionBCs) ? meshMotionBCs->getDofType(): 0;

    if(DofType)
    {
      #pragma omp parallel for
      for(int iSub=0; iSub<numLocSub; ++iSub)
      {
        double (*dx)[3] = dX.subData(iSub);
        int (*dofType)[3] = reinterpret_cast<int (*)[3]>(DofType[iSub]);
        for(int i=0;i<dX.subSize(iSub); i++)
          for(int l=0; l<3; l++)
            if(dofType[i][l]!=BC_MATCHED) dx[i][l] = 0.0;
      }
    }


    mms->applyProjector(Xdot); //HB: make sure Xdot satisfies the sliding conditions
    mms->solve(dX, X); //HB: the sliding conditions are also applied to dX inside the solve method

  }
  else
    X += dX;

  return dt;

}

//------------------------------------------------------------------------------

double PitchingMeshMotionHandler::updateStep1(bool *lastIt, int it, double t,
                                              DistSVec<double,3> &Xdot, DistSVec<double,3> &X, double *tmax)
{

  return dt;

}

//------------------------------------------------------------------------------

double PitchingMeshMotionHandler::updateStep2(bool *lastIt, int it, double t,
                                              DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{

  return update(lastIt, it, t, Xdot, X);

}

//------------------------------------------------------------------------------

HeavingMeshMotionHandler::HeavingMeshMotionHandler(IoData &iod, Domain *dom) :
                            MeshMotionHandler(iod, dom)
{

  dt = iod.forced.timestep;
  omega = 2.0 * acos(-1.0) * iod.forced.frequency;
  delta[0] = iod.forced.hv.ax;
  delta[1] = iod.forced.hv.ay;
  delta[2] = iod.forced.hv.az;

  if (iod.forced.hv.domain == HeavingData::VOLUME)
    mms = 0;
  else if (iod.forced.hv.domain == HeavingData::SURFACE)
    mms = new TetMeshMotionSolver(iod.dmesh, 0, domain, 0);

}

//------------------------------------------------------------------------------

HeavingMeshMotionHandler::~HeavingMeshMotionHandler()
{

  if (mms) delete mms;

}

//------------------------------------------------------------------------------

void HeavingMeshMotionHandler::setup(DistSVec<double,3> &X)
{

  if(mms) mms->setup(X);

}


//------------------------------------------------------------------------------

DistSVec<double,3> HeavingMeshMotionHandler::getModes()
{

  int numLocSub = domain->getNumLocSub();

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    double (*dx)[3] = dX.subData(iSub);
    for (int i=0; i<dX.subSize(iSub); ++i) {

      dx[i][0] = delta[0];
      dx[i][1] = delta[1];
      dx[i][2] = delta[2];

    }
  }

 return(dX);

}

//------------------------------------------------------------------------------

double HeavingMeshMotionHandler::update(bool *lastIt, int it, double t,
                                        DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{
  
 if (*lastIt) return dt;

 double hsintheta = 1.0 - cos(omega * (t + dt));
   // double hsintheta = t + dt;

 int numLocSub = domain->getNumLocSub();
 this->com->printf(9, "debugging: inside Heaving Mesh Motion Handler\n");
 #pragma omp parallel for
   for (int iSub=0; iSub<numLocSub; ++iSub) {

      double (*dx)[3] = dX.subData(iSub);
      double (*x0)[3] = X0.subData(iSub);

      for (int i=0; i<dX.subSize(iSub); ++i) {

           dx[i][0] = x0[i][0] + delta[0]*hsintheta;
           dx[i][1] = x0[i][1] + delta[1]*hsintheta;
           dx[i][2] = x0[i][2] + delta[2]*hsintheta;

        }

    }

  dX -= X;

  if (mms)
  { //HB: changed to use dofType array instead of nodeType array
    int numLocSub = domain->getNumLocSub();
    BCApplier* meshMotionBCs = domain->getMeshMotionBCs();
    int** DofType = (meshMotionBCs) ? meshMotionBCs->getDofType(): 0;

    if(DofType)
    {
      #pragma omp parallel for
      for(int iSub=0; iSub<numLocSub; ++iSub)
      {
        double (*dx)[3] = dX.subData(iSub);
        int (*dofType)[3] = reinterpret_cast<int (*)[3]>(DofType[iSub]);
        for(int i=0;i<dX.subSize(iSub); i++)
          for(int l=0; l<3; l++)
            if(dofType[i][l]!=BC_MATCHED) dx[i][l] = 0.0;
      }
    }


    mms->applyProjector(Xdot); //HB: make sure Xdot satisfies the sliding conditions
    mms->solve(dX, X); //HB: the sliding conditions are also applied to dX inside the solve method

  }
  else
    X += dX;

  return dt;

}

//------------------------------------------------------------------------------

double HeavingMeshMotionHandler::updateStep1(bool *lastIt, int it, double t,
                                             DistSVec<double,3> &Xdot, DistSVec<double,3> &X, double *tmax)
{

  return dt;

}

//------------------------------------------------------------------------------

double HeavingMeshMotionHandler::updateStep2(bool *lastIt, int it, double t,
                                             DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{

  return update(lastIt, it, t, Xdot, X);

}
//------------------------------------------------------------------------------
//
SpiralingMeshMotionHandler::SpiralingMeshMotionHandler(IoData &iod, Domain *dom) :
                            MeshMotionHandler(iod, dom)
{

  dt = iod.forced.timestep;
  omega = 2.0 * acos(-1.0) * iod.forced.frequency;
  delta[0] = iod.forced.sp.xL;
  delta[1] = iod.forced.sp.x0;

  if (iod.forced.sp.domain == SpiralingData::VOLUME)
    mms = 0;
  else if (iod.forced.sp.domain == SpiralingData::SURFACE)
    mms = new TetMeshMotionSolver(iod.dmesh, 0, domain, 0);

}

//------------------------------------------------------------------------------

SpiralingMeshMotionHandler::~SpiralingMeshMotionHandler()
{

  if (mms) delete mms;

}

//------------------------------------------------------------------------------

void SpiralingMeshMotionHandler::setup(DistSVec<double,3> &X)
{

  if(mms) mms->setup(X);

}

//------------------------------------------------------------------------------

DistSVec<double,3> SpiralingMeshMotionHandler::getModes()
{
// This is not doing the correct thing...Vinod

  int numLocSub = domain->getNumLocSub();

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    double (*dx)[3] = dX.subData(iSub);
    for (int i=0; i<dX.subSize(iSub); ++i) {

      dx[i][0] = 0.0;
      dx[i][1] = 0.0;
      dx[i][2] = 0.0;

    }
  }

 return(dX);

}

//------------------------------------------------------------------------------

double SpiralingMeshMotionHandler::update(bool *lastIt, int it, double t,
                                          DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{
  
 if (*lastIt) return dt;

 int numLocSub = domain->getNumLocSub();

 double Rcurv = delta[0]/( omega*(t+dt) );

 #pragma omp parallel for
   for (int iSub=0; iSub<numLocSub; ++iSub) {

      double (*dx)[3] = dX.subData(iSub);
      double (*x0)[3] = X0.subData(iSub);

      for (int i=0; i<dX.subSize(iSub); ++i) {
        if ( x0[i][0] >= delta[1] ) {
          double xloc = x0[i][0] - delta[1];
          double theta = xloc/Rcurv;
          dx[i][0] = (Rcurv + x0[i][1])*sin(theta) + delta[1];
          dx[i][1] = (Rcurv + x0[i][1])*cos(theta) - Rcurv;
          dx[i][2] = x0[i][2];
        }
        else {
          dx[i][0] = x0[i][0];
          dx[i][1] = x0[i][1];
          dx[i][2] = x0[i][2];
        }
      }
  }

  dX -= X;

  if (mms)
  { //HB: changed to use dofType array instead of nodeType array
    int numLocSub = domain->getNumLocSub();
    BCApplier* meshMotionBCs = domain->getMeshMotionBCs();
    int** DofType = (meshMotionBCs) ? meshMotionBCs->getDofType(): 0;

    if(DofType)
    {
      #pragma omp parallel for
      for(int iSub=0; iSub<numLocSub; ++iSub)
      {
        double (*dx)[3] = dX.subData(iSub);
        int (*dofType)[3] = reinterpret_cast<int (*)[3]>(DofType[iSub]);
        for(int i=0;i<dX.subSize(iSub); i++)
          for(int l=0; l<3; l++)
            if(dofType[i][l]!=BC_MATCHED) dx[i][l] = 0.0;
      }
    }


    mms->applyProjector(Xdot); //HB: make sure Xdot satisfies the sliding conditions
    mms->solve(dX, X); //HB: the sliding conditions are also applied to dX inside the solve method

  }
  else
    X += dX;

  return dt;

}

//------------------------------------------------------------------------------

double SpiralingMeshMotionHandler::updateStep1(bool *lastIt, int it, double t,
                                               DistSVec<double,3> &Xdot, DistSVec<double,3> &X, double *tmax)
{

  return dt;

}

//------------------------------------------------------------------------------

double SpiralingMeshMotionHandler::updateStep2(bool *lastIt, int it, double t,
                                               DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{

  return update(lastIt, it, t, Xdot, X);

}

//------------------------------------------------------------------------------

AccPitchingMeshMotionHandler::
AccPitchingMeshMotionHandler(IoData &iod, VarFcn *vf, double *Vin, Domain *dom) :
  PitchingMeshMotionHandler(iod, dom), RigidMeshMotionHandler(iod, vf, Vin, dom)
{

}

//------------------------------------------------------------------------------

double AccPitchingMeshMotionHandler::update(bool *lastIt, int it, double t,
                                          DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{

  DistSVec<double,3> &Xrel = getRelativePositionVector(t, X);

  double _dt = PitchingMeshMotionHandler::update(lastIt, it, t, Xdot, Xrel);

  if (*lastIt) return _dt;

  addRigidMotion(t + _dt, Xrel, Xdot, X, com);

  return _dt;

}

//------------------------------------------------------------------------------

double AccPitchingMeshMotionHandler::updateStep1(bool *lastIt, int it, double t,
                                               DistSVec<double,3> &Xdot, DistSVec<double,3> &X, double *tmax)
{

 DistSVec<double,3> &Xrel = getRelativePositionVector(t, X);

 double _dt = PitchingMeshMotionHandler::updateStep1(lastIt, it, t, Xdot, Xrel);
 
 return _dt;
 
}

//------------------------------------------------------------------------------

double AccPitchingMeshMotionHandler::updateStep2(bool *lastIt, int it, double t,
                                          DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{

  //DistSVec<double,3> &Xrel = getRelativePositionVector(t, X);

  //double _dt = PitchingMeshMotionHandler::updateStep2(lastIt, it, t, Xdot, Xrel);

  //if (*lastIt) return _dt;

  //addRigidMotion(t + _dt, Xrel, Xdot, X, com);

  //return _dt;

  double _dt = update(lastIt, it, t, Xdot, X);

  return _dt;
}
//------------------------------------------------------------------------------

AccHeavingMeshMotionHandler::
AccHeavingMeshMotionHandler(IoData &iod, VarFcn *vf, double *Vin, Domain *dom) :
  HeavingMeshMotionHandler(iod, dom), RigidMeshMotionHandler(iod, vf, Vin, dom)
{

}

//------------------------------------------------------------------------------

double AccHeavingMeshMotionHandler::update(bool *lastIt, int it, double t,
                                             DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{

  DistSVec<double,3> &Xrel = getRelativePositionVector(t, X);

  double _dt = HeavingMeshMotionHandler::update(lastIt, it, t, Xdot, Xrel);

  if (*lastIt) return _dt;

  addRigidMotion(t + _dt, Xrel, Xdot, X, com);

  return _dt;

}

//------------------------------------------------------------------------------

double AccHeavingMeshMotionHandler::updateStep1(bool *lastIt, int it, double t,
                                          DistSVec<double,3> &Xdot, DistSVec<double,3> &X, double *tmax)
{

 DistSVec<double,3> &Xrel = getRelativePositionVector(t, X);
 
 double _dt = HeavingMeshMotionHandler::updateStep1(lastIt, it, t, Xdot, Xrel);

 return _dt;

}

//------------------------------------------------------------------------------

double AccHeavingMeshMotionHandler::updateStep2(bool *lastIt, int it, double t,
                                          DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{

 // DistSVec<double,3> &Xrel = getRelativePositionVector(t, X);

 // double _dt = HeavingMeshMotionHandler::updateStep2(lastIt, it, t, Xdot, Xrel);

 // if (*lastIt) return _dt;

 // addRigidMotion(t + _dt, Xrel, Xdot, X, com);

 // return _dt;

  double  _dt = update(lastIt, it, t, Xdot, X);

  return _dt;

}
//------------------------------------------------------------------------------

AccDeformingMeshMotionHandler::
AccDeformingMeshMotionHandler(IoData &iod, VarFcn *vf, double *Vin, Domain *dom) :
  DeformingMeshMotionHandler(iod, dom), RigidMeshMotionHandler(iod, vf, Vin, dom)
{

}

//------------------------------------------------------------------------------

double AccDeformingMeshMotionHandler::update(bool *lastIt, int it, double t,
					  DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{

  DistSVec<double,3> &Xrel = getRelativePositionVector(t, X);

  double _dt = DeformingMeshMotionHandler::update(lastIt, it, t, Xdot, Xrel);

  if (*lastIt) return _dt;

  addRigidMotion(t + _dt, Xrel, Xdot, X, com);

  return _dt;

}

//------------------------------------------------------------------------------

double AccDeformingMeshMotionHandler::updateStep1(bool *lastIt, int it, double t,
                                          DistSVec<double,3> &Xdot, DistSVec<double,3> &X, double *tmax)
{

  DistSVec<double,3> &Xrel = getRelativePositionVector(t, X);

  double dt = DeformingMeshMotionHandler::updateStep1(lastIt, it, t, Xdot, Xrel);

  return dt;
 
}

//------------------------------------------------------------------------------

double AccDeformingMeshMotionHandler::updateStep2(bool *lastIt, int it, double t,
                                          DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{

  //DistSVec<double,3> &Xrel = getRelativePositionVector(t, X);

 // double _dt = DeformingMeshMotionHandler::updateStep2(lastIt, it, t, Xdot, Xrel);

  //if (*lastIt) return _dt;

  //addRigidMotion(t + _dt, Xrel, Xdot, X, com);

  //return _dt;

  return -1e-15;

}

//------------------------------------------------------------------------------

RigidRollMeshMotionHandler::
RigidRollMeshMotionHandler(IoData &ioData, double *angles, Domain *dom) : 
  MeshMotionHandler(ioData, dom)
{

  alpha_in = angles[0];
  beta_in = angles[1];

  xo[0] = 19.764697;
  xo[1] = 0.000000;
  xo[2] = 0.041863;

  /*
  xo[0] = 0.0;
  xo[1] = 0.0;
  xo[2] = 0.0;
  */

  u[0] = cos(alpha_in) * cos(beta_in);
  u[1] = sin(beta_in);
  u[2] = sin(alpha_in) * cos(beta_in);

  dt = ioData.ts.timestep;
  maxtime = ioData.ts.maxTime;

  mms = 0;//new TetMeshMotionSolver(ioData.dmesh, 0, domain, 0);

}

//------------------------------------------------------------------------------

double RigidRollMeshMotionHandler::computeRotationAngle(double time)
{

  double pi                = acos(-1.0);

  double t_switch1         = 0.0;
  double t_switch2         = maxtime;
  double alphaMaxInDegrees = 360.0;
  double degreesToRadians  = pi/180.0;
  double alphaMaxInRadians = alphaMaxInDegrees*degreesToRadians;
  double alpha1inDegrees   = 0.0;
  double alpha1            = alpha1inDegrees*degreesToRadians;

#ifdef CONSTANT_ROLL_RATE
  double a = (alphaMaxInRadians - alpha1) / (t_switch2 - t_switch1);
  double b = alpha1 - a * t_switch1;
#endif

  double alpha;

  if (time < t_switch1)
    alpha = alpha1;
  else if ( (time >= t_switch1) && (time <= t_switch2) ) {
#ifdef CONSTANT_ROLL_RATE
    alpha = a * time + b;
#else
    alpha = 0.5*alphaMaxInRadians * ( 1.0 - cos(time*pi/t_switch2) );
#endif
  }
  else
    alpha = alphaMaxInRadians;

  return alpha;

}

//------------------------------------------------------------------------------

double RigidRollMeshMotionHandler::update(bool *lastIt, int it, double t, 
                                          DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{

  if (*lastIt) return dt;

  double theta = computeRotationAngle(t + dt);

  com->printf(5, "Rotation angle = %e degrees at time = %e\n", 
	      theta * 180.0/acos(-1.0), (t + dt) * tscale);

  double R[3][3];

  R[0][0] = cos(alpha_in)*cos(alpha_in) + sin(alpha_in)*sin(alpha_in)*cos(theta);
  R[0][1] = -sin(alpha_in)*sin(theta);
  R[0][2] = (1.0 - cos(theta)) * cos(alpha_in)*sin(alpha_in);

  R[1][0] = sin(theta)*sin(alpha_in);
  R[1][1] = cos(theta);
  R[1][2] = -sin(theta)*cos(alpha_in);

  R[2][0] = (1.0 - cos(theta)) * cos(alpha_in)*sin(alpha_in);
  R[2][1] = cos(alpha_in)*sin(theta);
  R[2][2] = sin(alpha_in)*sin(alpha_in) + cos(alpha_in)*cos(alpha_in)*cos(theta);

  R[0][0] -= 1.0;
  R[1][1] -= 1.0;
  R[2][2] -= 1.0;

  int numLocSub = domain->getNumLocSub();

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {

    double (*dx)[3] = dX.subData(iSub);
    double (*x0)[3] = X0.subData(iSub);

    for (int i=0; i<dX.subSize(iSub); ++i) {

      Vec3D p0(x0[i]);
      double lambda = (p0 - xo) * u / (u * u);
      Vec3D xr = xo + lambda * u;
      Vec3D dp = p0 - xr;

      dx[i][0] = R[0][0]*dp[0] + R[0][1]*dp[1] + R[0][2]*dp[2];
      dx[i][1] = R[1][0]*dp[0] + R[1][1]*dp[1] + R[1][2]*dp[2];
      dx[i][2] = R[2][0]*dp[0] + R[2][1]*dp[1] + R[2][2]*dp[2];

    }

  }

  dX -= X - X0;

  if (mms) { //HB: need to be changed -> used dofType array instead of nodeType 
             // (similar to ForcedMeshMotionHandler::update)
    int **nodeType = domain->getNodeType();

#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; ++iSub) {

      double (*dx)[3] = dX.subData(iSub);

      for (int i=0; i<dX.subSize(iSub); ++i) {
	if (nodeType[iSub][i] >= BC_INTERNAL) {
	  dx[i][0] = 0.0;
	  dx[i][1] = 0.0;
	  dx[i][2] = 0.0;
	}
      }

    }

    mms->solve(dX, X);
  }
  else
    X += dX;

  return dt;

}

//------------------------------------------------------------------------------

double RigidRollMeshMotionHandler::updateStep2(bool *lastIt, int it, double t,
                                               DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{
 return update(lastIt, it, t, Xdot, X);
}

//------------------------------------------------------------------------------

EmbeddedMeshMotionHandler::EmbeddedMeshMotionHandler(IoData &iod, Domain *dom, DynamicNodalTransfer *dnTran,
  DistLevelSetStructure *distlss) : MeshMotionHandler(iod, dom)
{
  steady = !iod.problem.type[ProblemData::UNSTEADY];
  dynNodalTransfer = dnTran;
  distLSS = distlss;
  dts = 0.0;
  it0 = iod.restart.iteration; //restart time-step
}
 
//------------------------------------------------------------------------------

EmbeddedMeshMotionHandler::~EmbeddedMeshMotionHandler()
{ /*nothing*/ }

//------------------------------------------------------------------------------

void EmbeddedMeshMotionHandler::setup(double *maxTime)
{
  if(!dynNodalTransfer) {
    fprintf(stderr,"EmbeddedMeshMotionHandler is not initialized correctly!\n");
    exit(-1);
  }

  *maxTime = dynNodalTransfer->getStructureMaxTime();

}
//------------------------------------------------------------------------------

double EmbeddedMeshMotionHandler::updateStep1(bool *lastIt, int it, double t,
                                              DistSVec<double,3> &Xdot, DistSVec<double,3> &X, double *tmax)
{
  Timer *timer;
  timer = domain->getTimer();
  double ttt = timer->getTime();
  if(!dynNodalTransfer) {
    fprintf(stderr,"EmbeddedMeshMotionHandler is not initialized correctly!\n");
    exit(-1);
  }
  int algNum = dynNodalTransfer->getAlgorithmNumber();
  switch (algNum) {
    case 6: //A6 with Aero-S
      step1ForA6(lastIt,it,t,Xdot,X);
      break;
    case 20: //C0_legacy with Aero-S
      step1ForC0FEM(lastIt,it,t,Xdot,X);
      break;
    case 21: //C0 with XED3D
      if(dynNodalTransfer->cracking()) {
        fprintf(stderr,"XED3D is not supported for FSI w/ cracking!\n"); exit(-1);}
      step1ForC0XFEM(lastIt,it,t,Xdot,X);
      break;
    case 22: //C0 with Dyna3D-XFEM or Aero-S
      step1ForC0XFEM3D(lastIt,it,t,Xdot,X);
      *tmax = dynNodalTransfer->getStructureMaxTime();
      break;
    case 0: case 4: case 5: case 7: // A0/A4/A5/A7 with Aero-S
      fprintf(stderr,"***Error: Aeroelastic algorithm A%1d is not supported for embedded framework.\n", algNum);
      exit(-1);
      break;
  }
  timer->removeForceAndDispComm(ttt); // do not count the communication time with the
                                      // structure in the mesh solution

  return (steady) ? 0 : dts;
}

//------------------------------------------------------------------------------

void EmbeddedMeshMotionHandler::step1ForA6(bool *lastIt, int it, double t,
                                           DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{
  dts = dynNodalTransfer->getStructureTimeStep();

  if (it==0&& dynNodalTransfer->isCoupled()) {
    dts *= 0.5;	

    int numStructNodes = dynNodalTransfer->numStNodes();
    if(numStructNodes != distLSS->getNumStructNodes()) {
      fprintf(stderr,"SOFTWARE BUG: numStructNodes = %d (in dynNodalTransfer) and %d (in intersector)!\n",
              numStructNodes, distLSS->getNumStructNodes());
      exit(-1);
    }
  }


  // send force
  if (it>0 && it>it0)
    dynNodalTransfer->sendForce(); //send force to structure
}

//------------------------------------------------------------------------------

void EmbeddedMeshMotionHandler::step1ForC0FEM(bool *lastIt, int it, double t,
                                              DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{
  dts = dynNodalTransfer->getStructureTimeStep();

  if (it==0) {
    dts *= 0.5;

    int numStructNodes = dynNodalTransfer->numStNodes();
    if(numStructNodes != distLSS->getNumStructNodes()) {
      fprintf(stderr,"SOFTWARE BUG: numStructNodes = %d (in dynNodalTransfer) and %d (in intersector)!\n",
              numStructNodes, distLSS->getNumStructNodes());
      exit(-1);
    }

    //get displacement.
    dynNodalTransfer->getDisplacement(); //receive displacement and velocity from structure.
    distLSS->updateStructure(dynNodalTransfer->getStNodes(), dynNodalTransfer->getStVelocity(), numStructNodes, dynNodalTransfer->getStElems());

//    X = X0;
  } 

  else if (it>1 && it>it0)
    dynNodalTransfer->sendForce(); //send force to structure

}

//------------------------------------------------------------------------------

void EmbeddedMeshMotionHandler::step1ForC0XFEM(bool *lastIt, int it, double t,
                                               DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{
  dts = dynNodalTransfer->getStructureTimeStep();

  if (it==0) {
    dts *= 0.5;
    dynNodalTransfer->getDisplacement(); //receive displacement and velocity from structure.
    distLSS->updateStructure(dynNodalTransfer->getStNodes(), dynNodalTransfer->getStVelocity(), 
                             dynNodalTransfer->numStNodes());

//    X = X0;
  }

  else if (it>1 || it==it0)
    dynNodalTransfer->sendForce(); //send force to structure

}

//------------------------------------------------------------------------------

void EmbeddedMeshMotionHandler::step1ForC0XFEM3D(bool *lastIt, int it, double t,
                                                 DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{
  if (it==0) {
    //NOTE: this is the first iteration. No need to update cracking -- it has been done in the constructor
    //      of EmbeddedStructure inside the constructor of dynNodalTransfer.
    int numStructNodes = dynNodalTransfer->numStNodes();
    if(numStructNodes != distLSS->getNumStructNodes()) {
      fprintf(stderr,"SOFTWARE BUG: numStructNodes = %d (in dynNodalTransfer) and %d (in intersector)!\n",
              numStructNodes, distLSS->getNumStructNodes());
      exit(-1);
    }

    dynNodalTransfer->getDisplacement(); //receive displacement and velocity from structure.
    distLSS->updateStructure(dynNodalTransfer->getStNodes(), dynNodalTransfer->getStVelocity(), numStructNodes, dynNodalTransfer->getStElems());

    dynNodalTransfer->sendForce(); //send force to structure
    if(dynNodalTransfer->structSubcycling()) {
      dynNodalTransfer->sendFluidSuggestedTimestep(this->dtf0);
    }
    dynNodalTransfer->updateInfo();
    dts = dynNodalTransfer->getStructureTimeStep(); //dts obtained at initialization (getInfo)

    if(dts<=0.0) fprintf(stderr,"WARNING: Obtained a non-positive structural timestep (%e)!\n", dts);
    dts *= 0.5;

//    X = X0;
  }

  else if(it==it0) {
    dynNodalTransfer->sendForce(); //send force to structure
    if(dynNodalTransfer->structSubcycling())
      dynNodalTransfer->sendFluidSuggestedTimestep(this->dtf0);
    dynNodalTransfer->updateInfo(); // PJSA 11/24/2010
    dts = dynNodalTransfer->getStructureTimeStep();
  }

  else if(it==1) {
    dynNodalTransfer->updateInfo();
    dts = dynNodalTransfer->getStructureTimeStep(); 
  }

  else if (it>1 && !*lastIt) {
    dynNodalTransfer->sendForce(); //send force to structure
    if(dynNodalTransfer->structSubcycling())
      dynNodalTransfer->sendFluidSuggestedTimestep(this->dtf0);
    dynNodalTransfer->updateInfo();
    dts = dynNodalTransfer->getStructureTimeStep(); //dts obtained at initialization (getInfo)
  }

  else { // lastIt
    dynNodalTransfer->sendForce(); //send force to structure
    if(dynNodalTransfer->structSubcycling())
      dynNodalTransfer->sendFluidSuggestedTimestep(this->dtf0);
  }
}

//------------------------------------------------------------------------------

double EmbeddedMeshMotionHandler::updateStep2(bool *lastIt, int it, double t,
                                              DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{
  Timer *timer;
  timer = domain->getTimer();
  double ttt = timer->getTime();
  if(!dynNodalTransfer || !distLSS) {
    fprintf(stderr,"EmbeddedMeshMotionHandler is not initialized correctly!\n");
    exit(-1);
  }

  int algNum = dynNodalTransfer->getAlgorithmNumber();

  switch (algNum) {
    case 1: //PP with Aero-S
      step2ForPP(lastIt,it,t,Xdot,X);
      break;
    case 6: //A6 with Aero-S
      step2ForA6(lastIt,it,t,Xdot,X);
      break;
    case 20: //C0_legacy with Aero-S
    case 21: //C0 with XED3D
      step2ForC0(lastIt,it,t,Xdot,X);
      break;
    case 22: //C0 with Dyna3D-XFEM or Aero-S
      step2ForC0XFEM3D(lastIt,it,t,Xdot,X);
      break;
    case 0: case 4: case 5: case 7: // A0/A4/A5/A7 with Aero-S
      fprintf(stderr,"***Error: Aeroelastic algorithm A%1d is not supported for embedded framework.\n", algNum);
      exit(-1);
      break;
  }
  timer->removeForceAndDispComm(ttt); // do not count the communication time with the
                                      // structure in the mesh solution

  return dts;
}

//------------------------------------------------------------------------------

void EmbeddedMeshMotionHandler::step2ForPP(bool *lastIt, int it, double t,
                                           DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{
  // get displacement
  if(it==it0) {
    dynNodalTransfer->getDisplacement(); //receive displacement and velocity from structure.
    distLSS->updateStructure(dynNodalTransfer->getStNodes(), dynNodalTransfer->getStVelocity(), 
                             dynNodalTransfer->numStNodes(), dynNodalTransfer->getStElems());
  }
  *lastIt = true;
}

//------------------------------------------------------------------------------
void EmbeddedMeshMotionHandler::step2ForA6(bool *lastIt, int it, double t,
                                           DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{
  dts = dynNodalTransfer->getStructureTimeStep();

  if(it==0 && dynNodalTransfer->isCoupled() )
    dts *= 0.5;

  if(it > it0 && steady) {
    dynNodalTransfer->getStructExc()->negotiateStopping(lastIt);  // [F] receive iSteady from structure
    com->printf(4," ... [F] received iSteady from structure\n");
    if (*lastIt) {
      return;
    }
  }

  // get displacement
  if(it==it0 || !*lastIt) {
    if(it>0 && dynNodalTransfer->cracking())
      dynNodalTransfer->getNewCracking();

    dynNodalTransfer->getDisplacement(); //receive displacement and velocity from structure.
    distLSS->updateStructure(dynNodalTransfer->getStNodes(), dynNodalTransfer->getStVelocity(), 
                             dynNodalTransfer->numStNodes(), dynNodalTransfer->getStElems());

//    X = X0;
  }
}

//------------------------------------------------------------------------------

void EmbeddedMeshMotionHandler::step2ForC0(bool *lastIt, int it, double t,
                                           DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{
  dts = dynNodalTransfer->getStructureTimeStep();

  if(it==0) {
    dts *= 0.5;
    dynNodalTransfer->sendForce(); //send force to structure
  } 
  else if (!*lastIt) { //get displacement
    if(dynNodalTransfer->cracking())
      dynNodalTransfer->getNewCracking();
    dynNodalTransfer->getDisplacement(); //receive displacement and velocity from structure.
    distLSS->updateStructure(dynNodalTransfer->getStNodes(), dynNodalTransfer->getStVelocity(), 
                             dynNodalTransfer->numStNodes(), dynNodalTransfer->getStElems());

//    X = X0;
  } 
  else // last iteration 
    dts = 0.0; 
}

//------------------------------------------------------------------------------

void EmbeddedMeshMotionHandler::step2ForC0XFEM3D(bool *lastIt, int it, double t,
                                                 DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{
  if(it==0) {
    dts = dynNodalTransfer->getStructureTimeStep();
    dts *= 0.5;
  } 
  else if (!*lastIt) { //get displacement
    if(dynNodalTransfer->cracking()) 
      dynNodalTransfer->getNewCracking();
    dynNodalTransfer->getDisplacement(); //receive displacement and velocity from structure.
    distLSS->updateStructure(dynNodalTransfer->getStNodes(), dynNodalTransfer->getStVelocity(), 
                             dynNodalTransfer->numStNodes(), dynNodalTransfer->getStElems());

//    X = X0;
  } 
  else // last iteration 
    dts = 0.0; 
}
//------------------------------------------------------------------------------

EmbeddedALEMeshMotionHandler::EmbeddedALEMeshMotionHandler(IoData &iod, Domain *dom,
                              MatchNodeSet **matchNodes,
                              DistLevelSetStructure *distlss) : MeshMotionHandler(iod, dom)
{
  dt = iod.ts.timestep;
  distLSS = distlss;

  Vec<Vec3D>& Xstruct = distLSS->getStructPosition_0();
  Xs0 = new double [3*distLSS->getNumStructNodes()];

  for (int i=0; i<Xstruct.size(); i++)
    for (int j=0; j<3; j++)
      Xs0[3*i + j] = Xstruct[i][j];

  cs = new EmbeddedCorotSolver(iod, matchNodes, domain, Xs0, distLSS->getNumStructNodes(),distLSS->getStructElems());

  if (iod.dmesh.type == DefoMeshMotionData::COROTATIONAL)
    mms = new EmbeddedALETetMeshMotionSolver(iod.dmesh, matchNodes, domain, 0);
  else
    mms = 0;

}

//------------------------------------------------------------------------------

EmbeddedALEMeshMotionHandler::~EmbeddedALEMeshMotionHandler()
{

  if (cs) delete cs;
  if (mms) delete mms;
  if (Xs0) delete Xs0;

}

//------------------------------------------------------------------------------

void EmbeddedALEMeshMotionHandler::setup(DistSVec<double,3> &X, DistSVec<double,3> &Xdot)
{
  Vec<Vec3D>& Xstruct = distLSS->getStructPosition_n();

  double *Xs = new double [3*distLSS->getNumStructNodes()];

  for (int i=0; i<Xstruct.size(); i++)
    for (int j=0; j<3; j++)
      Xs[3*i + j] = Xstruct[i][j];

  cs->setup(Xs, distLSS->getNumStructNodes());
  cs->applyProjector(Xdot);

  delete [] Xs;
}

//------------------------------------------------------------------------------

double EmbeddedALEMeshMotionHandler::update(bool *lastIt, int it, double t,
                                            DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{
  
  if (*lastIt) return dt;

  Vec<Vec3D>& Xstruct = distLSS->getStructPosition_np1();

  double *Xs = new double [3*distLSS->getNumStructNodes()];

  for (int i=0; i<Xstruct.size(); i++)
    for (int j=0; j<3; j++)
      Xs[3*i + j] = Xstruct[i][j];

  cs->solve(Xs, distLSS->getNumStructNodes(), X, dX);
  cs->applyProjector(Xdot);

  if (mms) {

//--------------------------------------------------------------
// The following is only done for safety. 
// The dX coming out of cs solve should already satisfy this condition 
    BCApplier* meshMotionBCs = domain->getMeshMotionBCs();
    int** DofType = (meshMotionBCs) ? meshMotionBCs->getDofType(): 0;

    if(DofType)
    {
      int numLocSub = domain->getNumLocSub();
      #pragma omp parallel for
      for(int iSub=0; iSub<numLocSub; ++iSub)
      {
        double (*dx)[3] = dX.subData(iSub);
        int (*dofType)[3] = reinterpret_cast<int (*)[3]>(DofType[iSub]);
        for(int i=0;i<dX.subSize(iSub); i++)
          for(int l=0; l<3; l++)
            if(dofType[i][l]!=BC_FIXED && dofType[i][l]!=BC_MATCHED && dofType[i][l]!=BC_MATCHEDSLIDE) dx[i][l] = 0.0;
      }
    }
//--------------------------------------------------------------

    mms->applyProjector(Xdot); 
    mms->solve(dX, X);

  }

  delete [] Xs;

  return dt;

}

//------------------------------------------------------------------------------

double EmbeddedALEMeshMotionHandler::updateStep1(bool *lastIt, int it, double t,
                                                 DistSVec<double,3> &Xdot, DistSVec<double,3> &X, double *tmax)
{

  return dt;

}

//------------------------------------------------------------------------------

double EmbeddedALEMeshMotionHandler::updateStep2(bool *lastIt, int it, double t,
                                                 DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{

  return update(lastIt, it, t, Xdot, X);

}

//------------------------------------------------------------------------------
RbmExtractor::RbmExtractor(IoData &iod, Domain *dom) 
  : MeshMotionHandler(iod, dom)
{

  name_in = new char[strlen(iod.input.prefix) + strlen(iod.input.positions) + 1];
  sprintf(name_in, "%s%s", iod.input.prefix, iod.input.positions);

  name_out1 = new char[strlen(iod.output.transient.prefix) + 
		       strlen(iod.output.transient.displacement) + 1];
  sprintf(name_out1, "%s%s", iod.output.transient.prefix, iod.output.transient.displacement);

  const char *suffix = ".fixed";
  name_out2 = new char[strlen(name_out1) + strlen(suffix) + 1];
  sprintf(name_out2, "%s%s", name_out1, suffix);

  cs = new CorotSolver(iod.dmesh, 0, domain);

}

//------------------------------------------------------------------------------

double RbmExtractor::update(bool *lastIt, int it,double t,
                            DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{

  double amp = 10.0;

  int i = 0;

  // Initialize the tag value
  double tag = 0.0;

  DistSVec<double,3> Xrigid(domain->getNodeDistInfo());

  Xrigid = X0;

  while (domain->readVectorFromFile(name_in, i, &tag, X)) {
    dX = X - Xrigid;
    cs->solve(dX, Xrigid);
    dX = amp * (X - Xrigid);
    X = Xrigid + dX;
    domain->writeVectorToFile(name_out1, i, tag, X);
    X = X0 + dX;
    domain->writeVectorToFile(name_out2, i, tag, X);
    ++i;
  }

  *lastIt = true;

  return 0.0;

}

//------------------------------------------------------------------------------

double RbmExtractor::updateStep2(bool *lastIt, int it,double t,
                                 DistSVec<double,3> &Xdot, DistSVec<double,3> &X)
{
  return update(lastIt, it, t, Xdot, X);
}

//------------------------------------------------------------------------------
