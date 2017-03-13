#include <MeshMotionSolver.h>

#include <MatchNode.h>
#include <StiffMatrix.h>
#include <CorotSolver.h>
#include <KspSolver.h>
#include <NewtonSolver.h>
#include <MeshMotionHandler.h>
#include <Communicator.h>
#include <MemoryPool.h>
#include <Timer.h>
#include <BCApplier.h> 

#include <cstdio>

#ifdef TYPE_PREC_MESH
#define PrecScalar TYPE_PREC_MESH
#else
#define PrecScalar double
#endif

//#define HB_MESHMOTION_DEBUG

//------------------------------------------------------------------------------

TetMeshMotionSolver::TetMeshMotionSolver
(
  DefoMeshMotionData &data, MatchNodeSet **matchNodes, 
  Domain *dom, MemoryPool *mp
) 
: domain(dom), adjointFlag(false), stiffFlag(false), sensitivityFlag(false)
{

  com = domain->getCommunicator();

  KspData &kspData = data.newton.ksp;
  PcData &pcData = data.newton.ksp.pc;

  typeElement = data.element;
  maxItsNewton = data.newton.maxIts;
  if (data.element == DefoMeshMotionData::TORSIONAL_SPRINGS || data.element == DefoMeshMotionData::BALL_VERTEX)
    maxItsNewton = 1;
  epsNewton = data.newton.eps;
  epsAbsResNewton = data.newton.epsAbsRes;
  epsAbsIncNewton = data.newton.epsAbsInc;
  maxItsLS = data.newton.lineSearch.maxIts;
  contractionLS = data.newton.lineSearch.rho; 
  sufficDecreaseLS = data.newton.lineSearch.c1;
  if (strcmp(data.newton.output, "") == 0)
    outputNewton = 0;
  else if (strcmp(data.newton.output, "stdout") == 0)
    outputNewton = stdout;
  else if (strcmp(data.newton.output, "stderr") == 0)
    outputNewton = stderr;
  else {
    outputNewton = fopen(data.newton.output, "w");
    if (!outputNewton) {
      this->com->fprintf(stderr, "*** Error: could not open \'%s\'\n", data.newton.output);
      exit(1);
    }
  }

  timer = domain->getTimer();

  F0 = new DistSVec<double,3>(domain->getNodeDistInfo());

  if (data.type == DefoMeshMotionData::COROTATIONAL)
    cs = new CorotSolver(data, matchNodes, domain);
  else
    cs = 0;

  //int **ndType = domain->getNodeType();
  int **ndType = 0;

  meshMotionBCs = domain->getMeshMotionBCs(); //HB

  if (meshMotionBCs)   {
    meshMotionBCs->setDofType(matchNodes);

  }

  mvp = new StiffMat<double,3>(domain, ndType, mp, meshMotionBCs);

  if (pcData.type == PcData::IDENTITY)
    pc = new IdentityPrec<3>(meshMotionBCs);
  else if (pcData.type == PcData::JACOBI)
    //pc = new JacobiPrec<PrecScalar,3>(DiagMat<PrecScalar,3>::DIAGONAL, domain, ndType, meshMotionBCs);
    pc = new JacobiPrec<PrecScalar,3>(DiagMat<PrecScalar,3>::DENSE, domain, ndType, meshMotionBCs);

  else if (pcData.type == PcData::AS || pcData.type == PcData::RAS || pcData.type == PcData::ASH || pcData.type == PcData::AAS)
    pc = new IluPrec<PrecScalar,3>(pcData, domain, ndType);

  if (kspData.type == KspData::RICHARDSON)
    ksp = new RichardsonSolver<DistSVec<double,3>, StiffMat<double,3>, KspPrec<3>, Communicator>
      (domain->getNodeDistInfo(), kspData, mvp, pc, com);
  else if (kspData.type == KspData::CG)
    ksp = new CgSolver<DistSVec<double,3>, StiffMat<double,3>, KspPrec<3>, Communicator>
      (domain->getNodeDistInfo(), kspData, mvp, pc, com);
  else if (kspData.type == KspData::GMRES)
    ksp = new GmresSolver<DistSVec<double,3>, StiffMat<double,3>, KspPrec<3>, Communicator>
      (domain->getNodeDistInfo(), kspData, mvp, pc, com);

  ns = new NewtonSolver<TetMeshMotionSolver>(this);

  volStiff = data.volStiff;

}  


//------------------------------------------------------------------------------

TetMeshMotionSolver::~TetMeshMotionSolver()
{
  if (F0) delete F0;

  if (cs) delete cs;
  if (ns) delete ns;
  if (mvp) delete mvp;
  if (pc) delete pc;
}  

void TetMeshMotionSolver::applyProjectorTranspose(DistSVec<double,3> &X)
{
   if(meshMotionBCs) meshMotionBCs->applyPt(X);
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//HB: X <- P.X where P is the projector onto the sliding type of constraints
void
TetMeshMotionSolver::applyProjector(DistSVec<double,3> &X)
{
   if(meshMotionBCs) meshMotionBCs->applyP(X);
}

//------------------------------------------------------------------------------

/*
  X = current configuration
  dX = relative displacement (i.e. with respect to X) of the boundaries
*/
int TetMeshMotionSolver::solve(DistSVec<double,3> &dX, DistSVec<double,3> &X)  {

  // HB: dX <- P.dX where P is the projector onto  the sliding type of constraints
  // WARNING: assume the homogeneous Dirichlet BCs have already been applied to dX
  //com->fprintf(stdout, "Received 'unprojected' incr disp = %e\n",dX.norm());

  applyProjector(dX); 

  //com->fprintf(stdout, "Received 'projected' incr disp = %e\n",dX.norm());

  dX0 = &dX;

  if (cs) cs->solve(dX, X, meshMotionBCs); //HB: we may also want to apply the rotation to 
  //the normal directions used in the projections in BCApplier ...Not currently done, because 
  //it is assumed that in most practical cases (to date) there would be only one sliding plane,
  //and this sliding plane would also be the symmetry plane so that the plane normal wouldn't be
  // affected by the rotation of the fluid mesh around the symmetry plane normal ...

  ns->solve(X);

  return 0;
}

//------------------------------------------------------------------------------


//------------------------------------------------------------------------------

int TetMeshMotionSolver::solveAdjoint(DistSVec<double,3> &rhs, DistSVec<double,3> &lambdaX)  {

  dX0 = &rhs;
  ns->solve(lambdaX);

  return 0;
}

//------------------------------------------------------------------------------

void TetMeshMotionSolver::set_dX0(DistSVec<double,3> &dX)
{
  dX0 = &dX;
}

//------------------------------------------------------------------------------


void TetMeshMotionSolver::setup(DistSVec<double,3> &X)
{
  if(cs) cs->setup(X);

}

//------------------------------------------------------------------------------

void TetMeshMotionSolver::printf(int verbose, const char *format, ...)
{
  if (com->cpuNum() == 0 && verbose <= com->getMaxVerbose()) {
    va_list args;
    va_start(args, format);
    vfprintf(stdout, format, args);
    ::fflush(stdout);
    va_end(args);
  }
}

//------------------------------------------------------------------------------

void TetMeshMotionSolver::fprintf(FILE *fp, const char *format, ...)
{
  if (com->cpuNum() == 0) {
    va_list args;
    va_start(args, format);
    vfprintf(fp, format, args);
    ::fflush(fp);
    va_end(args);
  }
}

//------------------------------------------------------------------------------

void TetMeshMotionSolver::computeFunction(int it, DistSVec<double,3> &X, 
					  DistSVec<double,3> &F) 
{

  DistMat<PrecScalar,3> *_pc = dynamic_cast<DistMat<PrecScalar,3> *>(pc);

  if(it == 0 && (typeElement == DefoMeshMotionData::NON_LINEAR_FE
     || typeElement == DefoMeshMotionData::NL_BALL_VERTEX)
     && adjointFlag) {
    com->fprintf(stderr, " *** WARNING: Currently, only linear mesh motion solver is allowed for adjoint sensitivity method\n");
    exit(-1);
  }

  // PJSA FIX
  if(it == 0 && (typeElement == DefoMeshMotionData::NON_LINEAR_FE 
     || typeElement == DefoMeshMotionData::NL_BALL_VERTEX)) {
    X += *dX0; 
  }

  if((sensitivityFlag && !stiffFlag) || !sensitivityFlag) {
    domain->computeStiffAndForce(typeElement, X, F, *mvp, _pc, volStiff);
    stiffFlag = true;
  }
  if(adjointFlag) F = *dX0;  // overwrite F for adjoint method

  // PJSA FIX 
    if (it == 0) {
      if(!(typeElement == DefoMeshMotionData::NON_LINEAR_FE
        || typeElement == DefoMeshMotionData::NL_BALL_VERTEX
        || adjointFlag) ) { // compute F0 <- F0 + [Kib*dXb,0] & X <- X + [0,dXb]
      mvp->BCs = 0;
      mvp->apply(*dX0, *F0);
      mvp->BCs = meshMotionBCs;
      F += *F0;
      X += *dX0;
      } else if(adjointFlag) {
         X += *dX0;
    }
  }

  // PJSA FIX
  if(meshMotionBCs) {
//		meshMotionBCs->applyD(F);
//		com->fprintf(stderr,"F*F in TetMeshMotionSolver::computeFunction after applyD in final PJSA FIX is %e.\n", F*F);
//		meshMotionBCs->applyP(F);
	    if(adjointFlag) meshMotionBCs->applyPDt(F);
	    else meshMotionBCs->applyPD(F);
	}

}


void TetMeshMotionSolver::computeStiffnessMatrix(DistSVec<double,3> &X)
{

  DistMat<PrecScalar,3> *_pc = dynamic_cast<DistMat<PrecScalar,3> *>(pc);
  DistSVec<double,3> F(X);

  if(!stiffFlag) {
    domain->computeStiffAndForce(typeElement, X, F, *mvp, _pc, volStiff);
    stiffFlag = true;
  }
}

//------------------------------------------------------------------------------




//------------------------------------------------------------------------------

void TetMeshMotionSolver::computeJacobian(int it, DistSVec<double,3> &X, 
					  DistSVec<double,3> &F) 
{

}

//------------------------------------------------------------------------------

void TetMeshMotionSolver::setOperators(DistSVec<double,3> &X)
{

  double t0 = timer->getTime();

  pc->setup();
  
  double t = timer->addMeshPrecSetupTime(t0);

  com->printf(6, "Mesh preconditioner computation: %f s\n", t);

}

//------------------------------------------------------------------------------

int TetMeshMotionSolver::solveLinearSystem(int it, DistSVec<double,3> &rhs, 
					   DistSVec<double,3> &dX) 
{

  double t0 = timer->getTime();

  dX = 0.0;

  ksp->setup(it, maxItsNewton, rhs);

  int lits = ksp->solve(rhs, dX);

  // PJSA FIX (note rhs has already been projected in computeFunction)
  if(meshMotionBCs) {
    if(adjointFlag) meshMotionBCs->applyPDt(dX);
    else meshMotionBCs->applyPD(dX);
  }

  if(adjointFlag) {
    mvp->BCs = 0;
    DistSVec<double,3> dummy(dX);
    dummy = 0.0;
    mvp->apply(dX, dummy);
    dX = dummy;
    mvp->BCs = meshMotionBCs;
  }

  timer->addMeshKspTime(t0);
  
  return lits;

}


//------------------------------------------------------------------------------

void TetMeshMotionSolver::apply(DistSVec<double,3> &b,
                                DistSVec<double,3> &Ab)
{
  mvp->apply(b, Ab);
}

//--------------------------------------------------------------------------------------------------

EmbeddedALETetMeshMotionSolver::EmbeddedALETetMeshMotionSolver
(
  DefoMeshMotionData &data, MatchNodeSet **matchNodes, 
  Domain *dom, MemoryPool *mp
) 
: TetMeshMotionSolver(dom)
{

  com = domain->getCommunicator();

  KspData &kspData = data.newton.ksp;
  PcData &pcData = data.newton.ksp.pc;

  typeElement = data.element;
  maxItsNewton = data.newton.maxIts;
  if (data.element == DefoMeshMotionData::TORSIONAL_SPRINGS || data.element == DefoMeshMotionData::BALL_VERTEX)
    maxItsNewton = 1;
  epsNewton = data.newton.eps;
  epsAbsResNewton = data.newton.epsAbsRes;
  epsAbsIncNewton = data.newton.epsAbsInc;
  maxItsLS = data.newton.lineSearch.maxIts;
  contractionLS = data.newton.lineSearch.rho;
  sufficDecreaseLS = data.newton.lineSearch.c1;
  if (strcmp(data.newton.output, "") == 0)
    outputNewton = 0;
  else if (strcmp(data.newton.output, "stdout") == 0)
    outputNewton = stdout;
  else if (strcmp(data.newton.output, "stderr") == 0)
    outputNewton = stderr;
  else {
    outputNewton = fopen(data.newton.output, "w");
    if (!outputNewton) {
      this->com->fprintf(stderr, "*** Error: could not open \'%s\'\n", data.newton.output);
      exit(1);
    }
  }


  timer = domain->getTimer();

  F0 = new DistSVec<double,3>(domain->getNodeDistInfo());

  cs = 0;

  //int **ndType = domain->getNodeType();
  int **ndType = 0;

  meshMotionBCs = domain->getMeshMotionBCs(); //HB

  if (meshMotionBCs)   {
    meshMotionBCs->setEmbeddedALEDofType(matchNodes);
  }

  mvp = new StiffMat<double,3>(domain, ndType, mp, meshMotionBCs);

  if (pcData.type == PcData::IDENTITY)
    pc = new IdentityPrec<3>(meshMotionBCs);
  else if (pcData.type == PcData::JACOBI)
    //pc = new JacobiPrec<PrecScalar,3>(DiagMat<PrecScalar,3>::DIAGONAL, domain, ndType, meshMotionBCs);
    pc = new JacobiPrec<PrecScalar,3>(DiagMat<PrecScalar,3>::DENSE, domain, ndType, meshMotionBCs);

  else if (pcData.type == PcData::AS || pcData.type == PcData::RAS || pcData.type == PcData::ASH || pcData.type == PcData::AAS)
    pc = new IluPrec<PrecScalar,3>(pcData, domain, ndType);

  if (kspData.type == KspData::RICHARDSON)
    ksp = new RichardsonSolver<DistSVec<double,3>, StiffMat<double,3>, KspPrec<3>, Communicator>
      (domain->getNodeDistInfo(), kspData, mvp, pc, com);
  else if (kspData.type == KspData::CG)
    ksp = new CgSolver<DistSVec<double,3>, StiffMat<double,3>, KspPrec<3>, Communicator>
      (domain->getNodeDistInfo(), kspData, mvp, pc, com);
  else if (kspData.type == KspData::GMRES)
    ksp = new GmresSolver<DistSVec<double,3>, StiffMat<double,3>, KspPrec<3>, Communicator>
      (domain->getNodeDistInfo(), kspData, mvp, pc, com);

  ns = new NewtonSolver<TetMeshMotionSolver>(this);

  volStiff = data.volStiff;

}  

//------------------------------------------------------------------------------
/*
  X = current configuration
  dX = relative displacement (i.e. with respect to X) of the boundaries
*/
int EmbeddedALETetMeshMotionSolver::solve(DistSVec<double,3> &dX, DistSVec<double,3> &X)  {

  dX0 = &dX;

  ns->solve(X);

  return 0;
}

//------------------------------------------------------------------------------
