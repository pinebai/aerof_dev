#include <ImplicitCoupledTsDesc.h>
#include <Domain.h>
#include <GeoSource.h>
#include <DistTimeState.h>
#include <MatVecProd.h>
#include <KspSolver.h>
#include <MemoryPool.h>
//#include <MultiGridPrec.h>
#include <cstring>

#ifdef TYPE_MAT
#define MatScalar TYPE_MAT
#else
#define MatScalar double
#endif

#ifdef TYPE_PREC
#define PrecScalar TYPE_PREC
#else
#define PrecScalar double
#endif

//------------------------------------------------------------------------------

template<int dim>
ImplicitCoupledTsDesc<dim>::
ImplicitCoupledTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitTsDesc<dim>(ioData, geoSource, dom)
{

  ImplicitData &implicitData = ioData.ts.implicit;

  if (implicitData.mvp == ImplicitData::FD || implicitData.mvp == ImplicitData::H1FD)
  {
    mvp = new MatVecProdFD<dim,dim>(implicitData, this->timeState, this->geoState, this->spaceOp, this->domain, ioData);
    if (ioData.output.rom.fdResiduals==ROMOutputData::FD_RESIDUALS_ON) this->mvp->setTsOutput(this->output);
  }
  else if (implicitData.mvp == ImplicitData::H1)
  {
    mvp = new MatVecProdH1<dim,MatScalar,dim>(this->timeState, this->spaceOp, this->domain, ioData);
  }
  else if (implicitData.mvp == ImplicitData::H2)
  {
    mvp = new MatVecProdH2<dim,MatScalar,dim>(ioData, this->varFcn, this->timeState, this->spaceOp, this->domain, this->geoState);
  }

#ifdef MVP_CHECK
  ImplicitData fddata;
  fddata.mvp = ImplicitData::FD;
  mvpfd1 = new MatVecProdFD<dim,dim>(fddata, this->timeState, this->geoState, this->spaceOp, this->domain, ioData);
#endif

  pc = ImplicitTsDesc<dim>::template 
    createPreconditioner<PrecScalar,dim>(implicitData.newton.ksp.ns.pc, this->domain);

  ksp = this->createKrylovSolver(this->getVecInfo(), implicitData.newton.ksp.ns, mvp, pc, this->com);

  MemoryPool mp;

  mvp->exportMemory(&mp);
  pc->exportMemory(&mp);

  this->mmh = this->createMeshMotionHandler(ioData, geoSource, &mp);

}

//------------------------------------------------------------------------------

template<int dim>
ImplicitCoupledTsDesc<dim>::~ImplicitCoupledTsDesc()
{

  if (mvp) delete mvp;
  if (pc) delete pc;
  if (ksp) delete ksp;

#ifdef MVP_CHECK
  if (mvpfd1) delete mvpfd1;
#endif

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitCoupledTsDesc<dim>::computeJacobian(int it, DistSVec<double,dim> &Q,
						 DistSVec<double,dim> &F)
{

  if(this->wallRecType==BcsWallData::CONSTANT) {
    mvp->evaluate(it, *this->X, *this->A, Q, F);
  } else {
    mvp->evaluate(*this->riemann1, it, *this->X, *this->A, Q, F);
  }

#ifdef MVP_CHECK
  DistSVec<double,dim> p(this->getVecInfo());
  DistSVec<double,dim> prod(this->getVecInfo());

  p = 1.e-2;
  mvp->apply(p, prod);
  this->domain->checkMatVecProd(prod, "mvp");

  computeFunction(it, Q, F);
  mvpfd1->evaluate(it, *this->X, *this->A, Q, F);
  mvpfd1->apply(p, prod);
  this->domain->checkMatVecProd(prod, "mvpfd1");

  this->com->barrier();
  exit(1);

#endif

#ifdef DD_check
//------------------------------------------------
  std::cout << "Testing MVP \n";

  DistSVec<double,dim>       p_x(this->getVecInfo());
  DistSVec<double,dim>    prod_x(this->getVecInfo());

  DistSVec<double,dim> prod_fd(this->getVecInfo());
  DistSVec<double,dim>err_prod(this->getVecInfo());
  DistSVec<double,dim>     Q_p(this->getVecInfo());
  DistSVec<double,dim>     Q_m(this->getVecInfo());
  DistSVec<double,dim>     F_p(this->getVecInfo());
  DistSVec<double,dim>     F_m(this->getVecInfo());

for(int dd=1; dd<=10; ++dd){

  p_x = pow(10.0, -dd);

  mvp->apply(p_x, prod_x);

  Q_p = Q + p_x;
  computeFunction(it, Q_p, F_p);

  Q_m = Q - p_x;
  computeFunction(it, Q_m, F_m);

  prod_fd = 0.5*(F_p - F_m);
 
  err_prod = prod_x - prod_fd;

  double err_mvp = err_prod.norm();

  fprintf(stderr, "********* %10.5e %24.16e\n ", pow(10.0, -dd), err_mvp);
  
  
 }
 exit(-1);
//------------------------------------------------
#endif


}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitCoupledTsDesc<dim>::setOperators(DistSVec<double,dim> &Q)
{

  double t0 = this->timer->getTime();
  
  DistMat<PrecScalar,dim>  *_pc = dynamic_cast<DistMat<PrecScalar,dim> *>(pc);
  DistMat<double,    dim> *_pc2 = dynamic_cast<DistMat<double,    dim> *>(pc);

  MultiGridPrec<PrecScalar,dim> *pmg = dynamic_cast<MultiGridPrec<PrecScalar,dim> *>(pc);

  if (_pc || _pc2) {

    MatVecProdFD<dim, dim> *mvpfd = dynamic_cast<MatVecProdFD<dim, dim> *>(mvp);
    MatVecProdH1<dim,MatScalar,dim> *mvph1 = dynamic_cast<MatVecProdH1<dim,MatScalar,dim> *>(mvp);
    MatVecProdH2<dim,MatScalar,dim> *mvph2 = dynamic_cast<MatVecProdH2<dim,MatScalar,dim> *>(mvp);

    if (mvpfd || mvph2) {
      if (_pc) {
        this->spaceOp->computeJacobian(*this->X, *this->A, Q, *_pc, this->timeState);
        this->timeState->addToJacobian(*this->A, *_pc, Q);
        this->spaceOp->applyBCsToJacobian(Q, *_pc);
      } else {
        this->spaceOp->computeJacobian(*this->X, *this->A, Q, *_pc2, this->timeState);
        this->timeState->addToJacobian(*this->A, *_pc2, Q);
        this->spaceOp->applyBCsToJacobian(Q, *_pc2);
        if (pmg) {
          if (!pmg->isInitialized())
            pmg->initialize();
          pmg->getData(*_pc2);
        }
      }
    }
    else if (mvph1) {
      JacobiPrec<PrecScalar,dim> *jac = dynamic_cast<JacobiPrec<PrecScalar,dim> *>(pc);
      IluPrec<PrecScalar,dim> *ilu = dynamic_cast<IluPrec<PrecScalar,dim> *>(pc);
      MultiGridPrec<PrecScalar,dim> *pmg = dynamic_cast<MultiGridPrec<PrecScalar,dim> *>(pc);
      
      if (jac) {
        jac->getData(*mvph1);
      } else if (ilu) {
        ilu->getData(*mvph1);
      } else if (pmg) {
        if (!pmg->isInitialized())
          pmg->initialize();
        pmg->getData(*mvph1);
      }

    }

  }
    
  pc->setup();
  
  double t = this->timer->addPrecSetupTime(t0);

  this->com->printf(6, "Fluid preconditioner computation: %f s\n", t);

}

//------------------------------------------------------------------------------

template<int dim>
int ImplicitCoupledTsDesc<dim>::solveLinearSystem(int it, DistSVec<double,dim> &b, 
						  DistSVec<double,dim> &dQ)
{

  double t0 = this->timer->getTime();

//*****************    Verification of Adjointness of Jacobian dF/dW     ******************************
/*
  DistSVec<double,dim> w(this->domain->getNodeDistInfo()),p(this->domain->getNodeDistInfo()),q(this->domain->getNodeDistInfo()),z(this->domain->getNodeDistInfo());
  DistSVec<double,dim> one(this->domain->getNodeDistInfo());
  DistSVec<double,dim> bp(this->domain->getNodeDistInfo());
  one = b;
  bp = b;

  double direcTime(0), adjoinTime(0);
  for(int i=0; i<100; ++i) {
    double s_seconds(0), e_seconds(0);
    s_seconds = this->timer->getTime();
    mvp->apply(one,q);
    double r = q*q;
    e_seconds = this->timer->getTime();
    direcTime += (e_seconds-s_seconds);
    double qq = q*q;

    s_seconds = this->timer->getTime();
    mvp->applyTranspose(q,z);
    e_seconds = this->timer->getTime();
    adjoinTime += (e_seconds-s_seconds);
    double m = z*one;
    double zz = z*z;
    if (this->com->cpuNum() == 0) {
      fprintf(stderr,"q norm is %20.15e and z norm is %20.15e\n", qq, zz);
      fprintf(stderr,"r = %20.15e,\tm= %20.15e, and rel.error is %20.15e", r, m, (r-m)/r);
      if(abs((r-m)/r)>1.0e-6) fprintf(stderr,"WARNING!\n");
      else fprintf(stderr,"\n");
    }
    one = (1.0/q.norm())*q;
  }
  if(this->com->cpuNum() == 0) {
    fprintf(stderr,"direct Time is %9.6f and adjoint Time is %9.6f\n", direcTime, adjoinTime);
    exit(-1);
  }

  if(this->com->cpuNum() == 0) { fprintf(stderr,"\n\ncheckint xAx and xA'x\n"); }
  for(int i=0; i<10; ++i) {
    mvp->apply(bp,w);
    double a = bp*w;

    mvp->applyTranspose(bp,p);
    double v = bp*p;

    double ww = w*w;
    double pp = p*p;
    if (this->com->cpuNum() == 0) {
      fprintf(stderr,"w norm is %20.15e and p norm is %20.15e\n", ww, pp);
      fprintf(stderr,"a = %20.15e\tv= %20.15e, and rel.error is %20.15e", a, v, (a-v)/a);
      if(abs((a-v)/a)>1.0e-6) fprintf(stderr,"WARNING!\n");
      else fprintf(stderr,"\n");
    }
    bp = (1.0/p.norm())*p;

  }
*/
/////////////////////////////////////////

  dQ = 0.0;

  ksp->setup(it, this->maxItsNewton, b);



  int lits = ksp->solve(b, dQ);

  if (lits == ksp->maxits) this->errorHandler->localErrors[ErrorHandler::SATURATED_LS] += 1;

  this->timer->addKspTime(t0);

  return lits;

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void ImplicitCoupledTsDesc<dim>::rstVarImplicitCoupledTsDesc(IoData &ioData)
{

#ifdef MVP_CHECK
    mvpfd1->rstSpaceOp(ioData, this->varFcn, this->spaceOp, false);
#endif

    mvp->rstSpaceOp(ioData, this->varFcn, this->spaceOp, false);

}

//------------------------------------------------------------------------------


template<int dim>
template<int neq>
KspSolver<DistSVec<double,neq>, MatVecProd<dim,neq>, KspPrec<neq>, Communicator> *
ImplicitCoupledTsDesc<dim>::createKrylovSolver(const DistInfo &info, KspData &kspdata,
          MatVecProd<dim,neq> *_mvp, KspPrec<neq> *_pc,
          Communicator *_com)
{

  KspSolver<DistSVec<double,neq>, MatVecProd<dim,neq>,
    KspPrec<neq>, Communicator> *_ksp = 0;

  if (kspdata.type == KspData::RICHARDSON) {
    _ksp = new RichardsonSolver<DistSVec<double,neq>, MatVecProd<dim,neq>,
      KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  } else if (kspdata.type == KspData::CG) {
    _ksp = new CgSolver<DistSVec<double,neq>, MatVecProd<dim,neq>,
      KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  } else if (kspdata.type == KspData::GMRES) {
    _ksp = new GmresSolver<DistSVec<double,neq>, MatVecProd<dim,neq>,
      KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
    if (this->kspBinaryOutput)  _ksp->setKspBinaryOutput(this->kspBinaryOutput);
  } else if (kspdata.type == KspData::GCR) {
     _ksp = new GcrSolver<DistSVec<double,neq>, MatVecProd<dim,neq>,
       KspPrec<neq>, Communicator>(info, kspdata, _mvp, _pc, _com);
  }
  return _ksp;

}


//------------------------------------------------------------------------------
