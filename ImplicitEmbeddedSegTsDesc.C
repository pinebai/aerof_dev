#include "IoData.h"
#include "GeoSource.h"
#include "Domain.h"
#include "LevelSet.h"
#include "DistTimeState.h"

#include <MatVecProd.h>
#include <KspSolver.h>
#include <SpaceOperator.h>
#include <NewtonSolver.h>


#ifdef TYPE_PREC
#define PrecScalar TYPE_PREC
#else
#define PrecScalar double
#endif

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
ImplicitEmbeddedSegTsDesc<dim,neq1,neq2>::
ImplicitEmbeddedSegTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom):
  ImplicitEmbeddedTsDesc<dim>(ioData, geoSource, dom), embeddedB1(dom->getNodeDistInfo()), embeddeddQ1(dom->getNodeDistInfo()), embeddedB2(dom->getNodeDistInfo()), embeddeddQ2(dom->getNodeDistInfo())
{

  ImplicitData &implicitData = ioData.ts.implicit;

  spaceOp1 = createSpaceOperator1(ioData, this->spaceOp);
  spaceOp2 = createSpaceOperator2(ioData, this->spaceOp);
  
  // MatVecProd, Prec and Krylov solver for Euler equations
  if (implicitData.mvp == ImplicitData::FD) 
  {
    mvp1 = new MatVecProdFD<dim,neq1>(implicitData,this->timeState, this->geoState,
                                      this->spaceOp,this->domain,ioData);
    mvp2 = new MatVecProdFD<dim,neq2>(implicitData,this->timeState, this->geoState,
                                      this->spaceOp,this->domain,ioData);
  }
  else if(implicitData.mvp == ImplicitData::H1) 
  {
    mvp1 = new MatVecProdH1<dim,double,neq1>(this->timeState, this->spaceOp1, this->domain,ioData);
    mvp2 = new MatVecProdH1<dim,double,neq2>(this->timeState, this->spaceOp2, this->domain,ioData);
  }
  else
  {
    this->com->fprintf(stdout, "*** Error: MatVecProdH2 is not available\n");
    exit(1);
  }
  
  pc1 = ImplicitEmbeddedTsDesc<dim>::template 
    createPreconditioner<neq1>(implicitData.newton.ksp.ns.pc, this->domain);

  ksp1 = this->createKrylovSolver(this->getVecInfo(), implicitData.newton.ksp.ns, mvp1, pc1, this->com);
  
  pc2 = ImplicitEmbeddedTsDesc<dim>::template 
    createPreconditioner<neq2>(implicitData.newton.ksp.tm.pc, this->domain);

  ksp2 = this->createKrylovSolver(this->getVecInfo(), implicitData.newton.ksp.tm, mvp2, pc2, this->com);
  
  typename MatVecProd<dim,neq1>::_fsi fsi1 = {

    this->distLSS,
    &this->nodeTag,
    this->riemann,
    this->linRecAtInterface,
    this->viscSecOrder,
    &this->Wtemp,
    this->riemannNormal,
    this->ghostPoints,
  };

  typename MatVecProd<dim,neq2>::_fsi fsi2 = {

    this->distLSS,
    &this->nodeTag,
    this->riemann,
    this->linRecAtInterface,
    this->viscSecOrder,
    &this->Wtemp,
    this->riemannNormal,
    this->ghostPoints,
  };

  mvp1->AttachStructure(fsi1); 
  mvp2->AttachStructure(fsi2); 
  
}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
ImplicitEmbeddedSegTsDesc<dim,neq1,neq2>::~ImplicitEmbeddedSegTsDesc()
{
  if (spaceOp1) delete spaceOp1;
  if (spaceOp2) delete spaceOp2;
  if (mvp1)   delete mvp1;
  if (mvp2)   delete mvp2;
  if (pc1)    delete pc1;
  if (pc2)    delete pc2;
  if (ksp1)   delete ksp1;
  if (ksp2)   delete ksp2;

}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
SpaceOperator<dim> *ImplicitEmbeddedSegTsDesc<dim,neq1,neq2>::
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
SpaceOperator<dim> *ImplicitEmbeddedSegTsDesc<dim,neq1,neq2>::
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
void ImplicitEmbeddedSegTsDesc<dim,neq1,neq2>::computeJacobian(int it, DistSVec<double,dim> &Q,
							DistSVec<double,dim> &F)
{

  MatVecProdH1<dim,double,neq1> *mvph11 = dynamic_cast<MatVecProdH1<dim,double,neq1> *>(mvp1);
  if (mvph11)  {
    mvph11->clearGhost(); 
  }

  MatVecProdH1<dim,double,neq2> *mvph12 = dynamic_cast<MatVecProdH1<dim,double,neq2> *>(mvp2);
  if (mvph12)  {
    mvph12->clearGhost(); 
  }

  mvp1->evaluate(it,*(this->X) ,*(this->A), Q, F);
  mvp2->evaluate(it,*(this->X) ,*(this->A), Q, F);

//  mvph1 = dynamic_cast<MatVecProdH1<dim,double,dim> *>(mvp);
//  if (mvph1 && this->ghostPoints) 
//    this->domain->populateGhostJacobian(*this->ghostPoints,Q, this->varFcn, *this->distLSS, this->nodeTag,*mvph1);

}

//------------------------------------------------------------------------------

template<int dim, int neq1, int neq2>
template<int neq>
void ImplicitEmbeddedSegTsDesc<dim,neq1,neq2>::setOperator(MatVecProd<dim,neq> *mvp, KspPrec<neq> *pc, 
						   DistSVec<double,dim> &Q)
{

  DistMat<double,neq> *_pc = dynamic_cast<DistMat<double,neq> *>(pc);

  if (_pc) {

    MatVecProdFD<dim, neq> *mvpfd = dynamic_cast<MatVecProdFD<dim, neq> *>(mvp);
    MatVecProdH1<dim,double,neq> *mvph1 = dynamic_cast<MatVecProdH1<dim,double,neq> *>(mvp);

    if (mvpfd)  {
      if (neq > 2)  {
        spaceOp1->computeJacobian(*this->X, *this->A, Q,this->distLSS, this->nodeTag, this->riemann,
                                   this->riemannNormal, this->ghostPoints, *_pc,this->timeState);
        this->timeState->addToJacobian(*this->A, *_pc, Q);
        spaceOp1->applyBCsToJacobian(Q, *_pc, this->distLSS);
      }
      else  {
        spaceOp2->computeJacobian(*this->X, *this->A, Q,this->distLSS, this->nodeTag, this->riemann,
                                   this->riemannNormal, this->ghostPoints, *_pc,this->timeState);
        this->timeState->addToJacobian(*this->A, *_pc, Q);
        spaceOp2->applyBCsToJacobian(Q, *_pc, this->distLSS);
      }
    }
    else if (mvph1) {
      JacobiPrec<double,neq> *jac = dynamic_cast<JacobiPrec<double,neq> *>(pc);
      IluPrec<double,neq> *ilu = dynamic_cast<IluPrec<double,neq> *>(pc);

      if (jac) 
        jac->getData(*mvph1);
      else if (ilu) 
        ilu->getData(*mvph1);
    }

  }

}

//------------------------------------------------------------------------------
template<int dim, int neq1, int neq2>
void ImplicitEmbeddedSegTsDesc<dim,neq1,neq2>::setOperators(DistSVec<double,dim> &Q)
{
  
  setOperator<neq1>(mvp1, pc1, Q);
  setOperator<neq2>(mvp2, pc2, Q);
  
  double t0 = this->timer->getTime();

  pc1->setup();
  pc2->setup();
  
  double t = this->timer->addPrecSetupTime(t0);
  
  this->com->printf(6, "Fluid preconditioner computation: %f s\n", t);
  
}

//------------------------------------------------------------------------------
template<int dim, int neq1, int neq2>
int ImplicitEmbeddedSegTsDesc<dim,neq1,neq2>::solveLinearSystem(int it, DistSVec<double,dim> &b,
				                   DistSVec<double,dim> &dQ)
{
 
  double t0 = this->timer->getTime();  

  this->embeddeddQ = 0.0;
  this->embeddedB.ghost() = 0.0;
  this->embeddedB.real() = b;

  this->embeddedB.split(embeddedB1,embeddedB2);

  embeddeddQ1 = 0.0;

  ksp1->setup(it, this->maxItsNewton, embeddedB1);
  
  int lits1 = ksp1->solve(embeddedB1, embeddeddQ1);

  this->timer->addKspTime(t0);

  t0 = this->timer->getTime();  

  embeddeddQ2 = 0.0;

  ksp2->setup(it, this->maxItsNewton, embeddedB2);
  
  int lits2 = ksp2->solve(embeddedB2, embeddeddQ2);

  if(lits1 == ksp1->maxits || lits2 == ksp2->maxits) this->errorHandler->localErrors[ErrorHandler::SATURATED_LS] += 1;

  this->embeddeddQ.merge(embeddeddQ1,embeddeddQ2);
 
  dQ = this->embeddeddQ.real();
  this->embeddedU.ghost() += this->embeddeddQ.ghost();
  
  this->timer->addKspTime(t0);
  
  return 0;
  
}

