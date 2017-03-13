#include <ImplicitSegTsDesc.h>

#include <BcFcn.h>
#include <BcDef.h>
#include <FluxFcn.h>
#include <FemEquationTermDesc.h>
#include <DistTimeState.h>
#include <GeoSource.h>
#include <DistGeoState.h>
#include <SpaceOperator.h>
#include <MatVecProd.h>
#include <KspSolver.h>
#include <MemoryPool.h>

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

template<int dim, int neq1, int neq2>
ImplicitSegTsDesc<dim,neq1,neq2>::
ImplicitSegTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  ImplicitTsDesc<dim>(ioData, geoSource, dom), b1(this->getVecInfo()), dQ1(this->getVecInfo()),
  b2(this->getVecInfo()), dQ2(this->getVecInfo())
{

  ImplicitData &implicitData = ioData.ts.implicit;

  spaceOp1 = createSpaceOperator1(ioData, this->spaceOp);
  spaceOp2 = createSpaceOperator2(ioData, this->spaceOp);

#ifdef MVP_CHECK
  ImplicitData fddata;
  fddata.mvp = ImplicitData::FD;
  mvpfd = new MatVecProdFD<dim, dim>(fddata, this->timeState, this->geoState, this->spaceOp, this->domain, ioData);
#endif

  switch (implicitData.mvp)
  {
    case ImplicitData::H2 :
      mvp1 = new MatVecProdH2<dim,MatScalar,neq1>(ioData, this->varFcn, this->timeState, spaceOp1, this->domain, this->geoState);
      mvp2 = new MatVecProdH1<dim,MatScalar,neq2>(this->timeState, spaceOp2, this->domain, ioData);
      break;
    //---------------------
    case ImplicitData::H1 :
      mvp1 = new MatVecProdH1<dim,MatScalar,neq1>(this->timeState, spaceOp1, this->domain, ioData);
      mvp2 = new MatVecProdH1<dim,MatScalar,neq2>(this->timeState, spaceOp2, this->domain, ioData);
      break;
    //---------------------
    case ImplicitData::FD :
    case ImplicitData::H1FD :
    default :
      mvp1 = new MatVecProdFD<dim,neq1>(implicitData, this->timeState, this->geoState, this->spaceOp, this->domain, ioData);
      mvp2 = new MatVecProdFD<dim,neq2>(implicitData, this->timeState, this->geoState, this->spaceOp, this->domain, ioData);
      break;
  }

  pc1 = ImplicitTsDesc<dim>::template
    createPreconditioner<PrecScalar,neq1>(implicitData.newton.ksp.ns.pc, this->domain);
  ksp1 = this->createKrylovSolver(this->getVecInfo(), implicitData.newton.ksp.ns, mvp1, pc1, this->com);
  pc2 = ImplicitTsDesc<dim>::template
    createPreconditioner<PrecScalar,neq2>(implicitData.newton.ksp.tm.pc, this->domain);
  ksp2 = this->createKrylovSolver(this->getVecInfo(), implicitData.newton.ksp.tm, mvp2, pc2, this->com);

  MemoryPool mp;

  mvp1->exportMemory(&mp);
  pc1->exportMemory(&mp);

  this->mmh = this->createMeshMotionHandler(ioData, geoSource, &mp);

}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
ImplicitSegTsDesc<dim,neq1,neq2>::~ImplicitSegTsDesc()
{

  if (spaceOp1) delete spaceOp1;
  if (spaceOp2) delete spaceOp2;
  if (mvp1) delete mvp1;
  if (mvp2) delete mvp2;
  if (pc1) delete pc1;
  if (pc2) delete pc2;
  if (ksp1) delete ksp1;
  if (ksp2) delete ksp2;

#ifdef MVP_CHECK
  if (mvpfd) delete mvpfd;
#endif

}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
SpaceOperator<dim> *ImplicitSegTsDesc<dim,neq1,neq2>::
createSpaceOperator1(IoData &ioData, SpaceOperator<dim> *spo)
{

  SpaceOperator<dim> *spo1 = new SpaceOperator<dim>(*spo, false);

  double gamma = ioData.schemes.ns.gamma;

  FluxFcn **ff1 = new FluxFcn*[BC_MAX_CODE - BC_MIN_CODE + 1];
  ff1 -= BC_MIN_CODE;
  if(BC_MAX_CODE-BC_MIN_CODE+1 < 22)
    fprintf(stderr,"Be prepared to see a segmentation fault shortly...\n");

  ff1[BC_SYMMETRY] = new FluxFcn(0,BC_SYMMETRY,ioData,this->varFcn,1);
  ff1[BC_MASSFLOW_OUTLET_MOVING] = new FluxFcn(0,BC_MASSFLOW_OUTLET_MOVING,ioData,this->varFcn,1);
  ff1[BC_MASSFLOW_OUTLET_FIXED] = new FluxFcn(0,BC_MASSFLOW_OUTLET_FIXED,ioData,this->varFcn,1);
  ff1[BC_MASSFLOW_INLET_MOVING] = new FluxFcn(0,BC_MASSFLOW_INLET_MOVING,ioData,this->varFcn,1);
  ff1[BC_MASSFLOW_INLET_FIXED] = new FluxFcn(0,BC_MASSFLOW_INLET_FIXED,ioData,this->varFcn,1);
  ff1[BC_DIRECTSTATE_OUTLET_MOVING] = new FluxFcn(0,BC_DIRECTSTATE_OUTLET_MOVING,ioData,this->varFcn,1);
  ff1[BC_DIRECTSTATE_OUTLET_FIXED] = new FluxFcn(0,BC_DIRECTSTATE_OUTLET_FIXED,ioData,this->varFcn,1);
  ff1[BC_DIRECTSTATE_INLET_MOVING] = new FluxFcn(0,BC_DIRECTSTATE_INLET_MOVING,ioData,this->varFcn,1);
  ff1[BC_DIRECTSTATE_INLET_FIXED] = new FluxFcn(0,BC_DIRECTSTATE_INLET_FIXED,ioData,this->varFcn,1);
  ff1[BC_OUTLET_MOVING] = new FluxFcn(0,BC_OUTLET_MOVING,ioData,this->varFcn,1);
  ff1[BC_OUTLET_FIXED] = new FluxFcn(0,BC_OUTLET_FIXED,ioData,this->varFcn,1);
  ff1[BC_INLET_MOVING] = new FluxFcn(0,BC_INLET_MOVING,ioData,this->varFcn,1);
  ff1[BC_INLET_FIXED] = new FluxFcn(0,BC_INLET_FIXED,ioData,this->varFcn,1);
  ff1[BC_POROUS_WALL_MOVING] = new FluxFcn(0,BC_POROUS_WALL_MOVING,ioData,this->varFcn,1);
  ff1[BC_POROUS_WALL_FIXED] = new FluxFcn(0,BC_POROUS_WALL_FIXED,ioData,this->varFcn,1);
  ff1[BC_ADIABATIC_WALL_MOVING] = new FluxFcn(0,BC_ADIABATIC_WALL_MOVING,ioData,this->varFcn,1);
  ff1[BC_ADIABATIC_WALL_FIXED] = new FluxFcn(0,BC_ADIABATIC_WALL_FIXED,ioData,this->varFcn,1);
  ff1[BC_SLIP_WALL_MOVING] = new FluxFcn(0,BC_SLIP_WALL_MOVING,ioData,this->varFcn,1);
  ff1[BC_SLIP_WALL_FIXED] = new FluxFcn(0,BC_SLIP_WALL_FIXED,ioData,this->varFcn,1);
  ff1[BC_ISOTHERMAL_WALL_MOVING] = new FluxFcn(0,BC_ISOTHERMAL_WALL_MOVING,ioData,this->varFcn,1);
  ff1[BC_ISOTHERMAL_WALL_FIXED] = new FluxFcn(0,BC_ISOTHERMAL_WALL_FIXED,ioData,this->varFcn,1);
  ff1[BC_INTERNAL] = new FluxFcn(0,BC_INTERNAL,ioData,this->varFcn,1);

  BcFcn *bf1 = 0;
  if (ioData.bc.wall.integration != BcsWallData::WALL_FUNCTION)
    bf1 = new BcFcnNS;

  FemEquationTerm *fet1 = 0;
  if (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS)
    fet1 = new FemEquationTermSAmean(ioData, this->varFcn);
  else if (ioData.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES)
    fet1 = new FemEquationTermDESmean(ioData, this->varFcn);
  else if (ioData.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE)
    fet1 = new FemEquationTermKEmean(ioData, this->varFcn);

  spo1->setFluxFcn(ff1);
  spo1->setBcFcn(bf1);
  spo1->setFemEquationTerm(fet1);

  return spo1;

}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
SpaceOperator<dim> *ImplicitSegTsDesc<dim,neq1,neq2>::
createSpaceOperator2(IoData &ioData, SpaceOperator<dim> *spo)
{

  SpaceOperator<dim> *spo2 = new SpaceOperator<dim>(*spo, false);

  double gamma = ioData.schemes.ns.gamma;

  FluxFcn** ff2 = new FluxFcn*[BC_MAX_CODE - BC_MIN_CODE + 1];
  ff2 -= BC_MIN_CODE;
  if(BC_MAX_CODE-BC_MIN_CODE+1 < 22)
    fprintf(stderr,"Be prepared to see a segmentation fault shortly...\n");

  ff2[BC_SYMMETRY] = new FluxFcn(0,BC_SYMMETRY,ioData,this->varFcn,2);
  ff2[BC_MASSFLOW_OUTLET_MOVING] = new FluxFcn(0,BC_MASSFLOW_OUTLET_MOVING,ioData,this->varFcn,2);
  ff2[BC_MASSFLOW_OUTLET_FIXED] = new FluxFcn(0,BC_MASSFLOW_OUTLET_FIXED,ioData,this->varFcn,2);
  ff2[BC_MASSFLOW_INLET_MOVING] = new FluxFcn(0,BC_MASSFLOW_INLET_MOVING,ioData,this->varFcn,2);
  ff2[BC_MASSFLOW_INLET_FIXED] = new FluxFcn(0,BC_MASSFLOW_INLET_FIXED,ioData,this->varFcn,2);
  ff2[BC_DIRECTSTATE_OUTLET_MOVING] = new FluxFcn(0,BC_DIRECTSTATE_OUTLET_MOVING,ioData,this->varFcn,2);
  ff2[BC_DIRECTSTATE_OUTLET_FIXED] = new FluxFcn(0,BC_DIRECTSTATE_OUTLET_FIXED,ioData,this->varFcn,2);
  ff2[BC_DIRECTSTATE_INLET_MOVING] = new FluxFcn(0,BC_DIRECTSTATE_INLET_MOVING,ioData,this->varFcn,2);
  ff2[BC_DIRECTSTATE_INLET_FIXED] = new FluxFcn(0,BC_DIRECTSTATE_INLET_FIXED,ioData,this->varFcn,2);
  ff2[BC_OUTLET_MOVING] = new FluxFcn(0,BC_OUTLET_MOVING,ioData,this->varFcn,2);
  ff2[BC_OUTLET_FIXED] = new FluxFcn(0,BC_OUTLET_FIXED,ioData,this->varFcn,2);
  ff2[BC_INLET_MOVING] = new FluxFcn(0,BC_INLET_MOVING,ioData,this->varFcn,2);
  ff2[BC_INLET_FIXED] = new FluxFcn(0,BC_INLET_FIXED,ioData,this->varFcn,2);
  ff2[BC_POROUS_WALL_MOVING] = new FluxFcn(0,BC_POROUS_WALL_MOVING,ioData,this->varFcn,2);
  ff2[BC_POROUS_WALL_FIXED] = new FluxFcn(0,BC_POROUS_WALL_FIXED,ioData,this->varFcn,2);
  ff2[BC_ADIABATIC_WALL_MOVING] = new FluxFcn(0,BC_ADIABATIC_WALL_MOVING,ioData,this->varFcn,2);
  ff2[BC_ADIABATIC_WALL_FIXED] = new FluxFcn(0,BC_ADIABATIC_WALL_FIXED,ioData,this->varFcn,2);
  ff2[BC_SLIP_WALL_MOVING] = new FluxFcn(0,BC_SLIP_WALL_MOVING,ioData,this->varFcn,2);
  ff2[BC_SLIP_WALL_FIXED] = new FluxFcn(0,BC_SLIP_WALL_FIXED,ioData,this->varFcn,2);
  ff2[BC_ISOTHERMAL_WALL_MOVING] = new FluxFcn(0,BC_ISOTHERMAL_WALL_MOVING,ioData,this->varFcn,2);
  ff2[BC_ISOTHERMAL_WALL_FIXED] = new FluxFcn(0,BC_ISOTHERMAL_WALL_FIXED,ioData,this->varFcn,2);
  ff2[BC_INTERNAL] = new FluxFcn(0,BC_INTERNAL,ioData,this->varFcn,2);

  BcFcn* bf2 = 0;
  FemEquationTerm* fet2 = 0;

  switch(ioData.eqs.tc.tm.type){
    case TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS:
      bf2 = new BcFcnSAturb;
      fet2 = new FemEquationTermSAturb(ioData, this->varFcn);
      break;
    case TurbulenceModelData::ONE_EQUATION_DES:
      bf2 = new BcFcnSAturb;
      fet2 = new FemEquationTermDESturb(ioData, this->varFcn);
      break;
    case TurbulenceModelData::TWO_EQUATION_KE:
      bf2 = new BcFcnKEturb;
      fet2 = new FemEquationTermKEturb(ioData, this->varFcn);
      break;
    default:
      fprintf(stderr,"Error: Seg. solver not implemented for this type of simulation!\n");
      exit(1);
      break;
  }

  spo2->setFluxFcn(ff2);
  spo2->setBcFcn(bf2);
  spo2->setFemEquationTerm(fet2);

  return spo2;

}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
void ImplicitSegTsDesc<dim,neq1,neq2>::computeJacobian(int it, DistSVec<double,dim> &Q,
						       DistSVec<double,dim> &F)
{
  if(this->wallRecType==BcsWallData::CONSTANT)
    mvp1->evaluate(it, *this->X, *this->A, Q, F);
  else
    mvp1->evaluate(*this->riemann1, it, *this->X, *this->A, Q, F);

  mvp2->evaluate(it, *this->X, *this->A, Q, F);

#ifdef MVP_CHECK
  DistSVec<double,dim> p(this->getVecInfo());
  DistSVec<double,dim> prod(this->getVecInfo());

  DistSVec<double,neq1> p1(this->getVecInfo());
  DistSVec<double,neq1> prod1(this->getVecInfo());

  DistSVec<double,neq2> p2(this->getVecInfo());
  DistSVec<double,neq2> prod2(this->getVecInfo());

  p1 = 0.0;
  p2 = 1.e-2;
  p.merge(p1, p2);

  mvp2->apply(p2, prod2);

  this->domain->checkMatVecProd(prod2, "mvp");

  computeFunction(it, Q, F);

  mvpfd->evaluate(it, *this->X, *this->A, Q, F);

  mvpfd->apply(p, prod);

  prod.split(prod1, prod2);

  this->domain->checkMatVecProd(prod2, "mvpfd");

  this->com->barrier();
  exit(1);
#endif

}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
template<int neq>
void ImplicitSegTsDesc<dim,neq1,neq2>::setOperator(MatVecProd<dim,neq> *mvp, KspPrec<neq> *pc,
						   DistSVec<double,dim> &Q, SpaceOperator<dim> *spo)
{

  DistMat<PrecScalar,neq> *_pc = dynamic_cast<DistMat<PrecScalar,neq> *>(pc);
  DistMat<double,neq> *_pc2 = dynamic_cast<DistMat<double,neq> *>(pc);

  MultiGridPrec<PrecScalar,neq> *pmg = dynamic_cast<MultiGridPrec<PrecScalar,neq> *>(pc);

  if (pmg) {
    if (!pmg->isInitialized())
      pmg->initialize();
  }

  if (_pc || _pc2) {

    MatVecProdFD<dim, neq> *mvpfdtmp = dynamic_cast<MatVecProdFD<dim, neq> *>(mvp);
    MatVecProdH1<dim,MatScalar,neq> *mvph1 = dynamic_cast<MatVecProdH1<dim,MatScalar,neq> *>(mvp);
    MatVecProdH2<dim,MatScalar,neq> *mvph2 = dynamic_cast<MatVecProdH2<dim,MatScalar,neq> *>(mvp);

    if ((mvpfdtmp) || (mvph2))  {
      if (_pc) {
        if (neq > 2)  {
          spaceOp1->computeJacobian(*this->X, *this->A, Q, *_pc, this->timeState);
          this->timeState->addToJacobian(*this->A, *_pc, Q);
          spaceOp1->applyBCsToJacobian(Q, *_pc);
        }
        else  {
          spaceOp2->computeJacobian(*this->X, *this->A, Q, *_pc, this->timeState);
          this->timeState->addToJacobian(*this->A, *_pc, Q);
          spaceOp2->applyBCsToJacobian(Q, *_pc);
        }
      } else {
        if (neq > 2)  {
          spaceOp1->computeJacobian(*this->X, *this->A, Q, *_pc2, this->timeState);
          this->timeState->addToJacobian(*this->A, *_pc2, Q);
          spaceOp1->applyBCsToJacobian(Q, *_pc2);
          if (pmg) {
            if (!pmg->isInitialized())
              pmg->initialize();
            pmg->getData(*_pc2);
          }

        }
        else  {
          spaceOp2->computeJacobian(*this->X, *this->A, Q, *_pc2, this->timeState);
          this->timeState->addToJacobian(*this->A, *_pc2, Q);
          spaceOp2->applyBCsToJacobian(Q, *_pc2);
          if (pmg) {
            if (!pmg->isInitialized())
              pmg->initialize();
            pmg->getData(*_pc2);
          }
        }
      }
    }
    else if (mvph1)
    {
      JacobiPrec<PrecScalar,neq> *jac = dynamic_cast<JacobiPrec<PrecScalar,neq> *>(pc);
      IluPrec<PrecScalar,neq> *ilu = dynamic_cast<IluPrec<PrecScalar,neq> *>(pc);
      //MultiGridPrec<PrecScalar,neq> *pmg = dynamic_cast<MultiGridPrec<PrecScalar,neq> *>(pc);
      if (jac)
        jac->getData(*mvph1);
      else if (ilu)
        ilu->getData(*mvph1);
      else if (pmg) {
        if (!pmg->isInitialized())
          pmg->initialize();
        pmg->getData(*mvph1);
      }
    }

  }

}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
void ImplicitSegTsDesc<dim,neq1,neq2>::setOperators(DistSVec<double,dim> &Q)
{

  setOperator(mvp1, pc1, Q, this->spaceOp);
  setOperator(mvp2, pc2, Q, this->spaceOp);

  double t0 = this->timer->getTime();

  pc1->setup();
  pc2->setup();

  double t = this->timer->addPrecSetupTime(t0);

  this->com->printf(6, "Fluid preconditioner computation: %f s\n", t);

}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
int ImplicitSegTsDesc<dim,neq1,neq2>::solveLinearSystem(int it, DistSVec<double,dim> &b,
							DistSVec<double,dim> &dQ)
{

  b.split(b1, b2);

  double t0 = this->timer->getTime();

  dQ1 = 0.0;

  ksp1->setup(it, this->maxItsNewton, b1);

  int lits1 = ksp1->solve(b1, dQ1);

  this->timer->addKspTime(t0);

  t0 = this->timer->getTime();

  dQ2 = 0.0;

  ksp2->setup(it,this->maxItsNewton, b2);

  int lits2 = ksp2->solve(b2, dQ2);

  if(lits1 == ksp1->maxits || lits2 == ksp2->maxits) this->errorHandler->localErrors[ErrorHandler::SATURATED_LS] += 1;

  this->timer->addKspTime(t0);

  dQ.merge(dQ1, dQ2);

  return 0;

}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
void ImplicitSegTsDesc<dim,neq1,neq2>::rstVarImplicitSegTsDesc(IoData &ioData)
{
  std::cout<<"Not implemented yet"<<std::endl; sleep(1); exit(-1);
}

//------------------------------------------------------------------------------
