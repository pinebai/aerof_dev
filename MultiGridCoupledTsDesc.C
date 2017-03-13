#include <VectorSet.h>
#include <LevelSet/LevelSetStructure.h>
#include <MultiGridCoupledTsDesc.h>

template <int dim>
MultiGridCoupledTsDesc<dim>::
MultiGridCoupledTsDesc(IoData & iod, GeoSource & gs,  Domain * dom) :
  ImplicitCoupledTsDesc<dim>(iod,gs,dom) {

  memset(numSmooths_post,0,sizeof(numSmooths_post));
  numSmooths_pre[0] = 1;
  numSmooths_pre[1] = 2;
  numSmooths_pre[2] = 3;
  numSmooths_pre[3] = 3;
  numSmooths_pre[4] = 3;
  numSmooths_pre[5] = 3;

  smoothWithGMRES = (iod.mg.mg_smoother == MultiGridData::MGGMRES);

  prolong_relax_factor = iod.mg.prolong_relax_factor;
  restrict_relax_factor = iod.mg.restrict_relax_factor;

  if (iod.mg.cycle_scheme == MultiGridData::VCYCLE)
    mc = 1;
  else if (iod.mg.cycle_scheme == MultiGridData::WCYCLE)
    mc = 2;     

  globalIt = 0;

  mgMvp = NULL;
  pKernel = NULL;
  mgSpaceOp = NULL;
  mgKspSolver = NULL;
  smoothingMatrices = NULL;

  TsDesc<dim>::isMultigridTsDesc = true;
}

template <int dim>
MultiGridCoupledTsDesc<dim>::
~MultiGridCoupledTsDesc() {

  if (mgMvp)
    delete mgMvp;
  if (mgSpaceOp)
    delete mgSpaceOp;
  if (mgKspSolver)
    delete mgKspSolver;
  if (smoothingMatrices)
    delete smoothingMatrices;

  V.free();
  res.free();
  R.free();
  F.free();
  Forig.free();
  U.free();
  Uold.free(); 
  dx.free();


 if (pKernel)
    delete pKernel;

}


template <int dim>
void MultiGridCoupledTsDesc<dim>::
setupTimeStepping(DistSVec<double,dim> *U0, IoData &iod) {

  ImplicitCoupledTsDesc<dim>::setupTimeStepping(U0,iod);

  pKernel = new MultiGridKernel<double>(this->domain, *this->geoState,
                                        iod,iod.mg.num_multigrid_levels);

  pKernel->initialize(dim,dim,0);

  mgSpaceOp = 
    new MultiGridSpaceOperator<double,dim>(iod, this->domain, this->spaceOp, pKernel);

  mgMvp = new MultiGridMvpMatrix<double,dim>(this->domain,pKernel);

  mgKspSolver = new MultiGridKspSolver<double,dim,double>(this->domain, iod.ts.implicit.newton.ksp.ns,
                                                          pKernel);

  if (!smoothWithGMRES) {
    typename MultiGridSmoothingMatrix<double,dim>::SmoothingMode s;
    if (iod.mg.mg_smoother == MultiGridData::MGRAS)
      s = MultiGridSmoothingMatrix<double,dim>::RAS; 
    if (iod.mg.mg_smoother == MultiGridData::MGJACOBI)
      s = MultiGridSmoothingMatrix<double,dim>::BlockJacobi;
    if (iod.mg.mg_smoother == MultiGridData::MGLINEJACOBI)
      s = MultiGridSmoothingMatrix<double,dim>::LineJacobi;
    smoothingMatrices = new MultiGridSmoothingMatrices<double,dim>(pKernel, s);

  }

  V.init(pKernel);
  res.init(pKernel);
  R.init(pKernel);
  F.init(pKernel);
  Forig.init(pKernel);
  U.init(pKernel);
  Uold.init(pKernel); 
  dx.init(pKernel);
}

template <int dim>
void MultiGridCoupledTsDesc<dim>::
smooth0(DistSVec<double,dim>& x,int steps) {

  int i;
  double dummy = 0.0;
  this->updateStateVectors(x, 0);
  for (i = 0; i < steps; ++i) {

    this->computeTimeStep(globalIt,&dummy, x);
    ++globalIt;
    this->computeFunction(0, x, R(0));
    this->computeJacobian(0, x, R(0));
    if (smoothWithGMRES)
      this->setOperators(x);
    else
      smoothingMatrices->acquire( *this->GetJacobian());
  
    R(0) *= -1.0;
    if (smoothWithGMRES)
      this->solveLinearSystem(0, R(0),dx(0));
    else
      smoothingMatrices->apply(0, dx, R);
    
    x += dx(0);
    
    this->updateStateVectors(x, 0);
    this->monitorConvergence(0, x);
    R(0) = -1.0*this->getCurrentResidual();
    
  }
  //if (i == 0) {
  //  this->monitorConvergence(0, x);
  //  R(0) = -1.0*this->getCurrentResidual();
 // }
  double one = 1.0;
  if (globalIt%10 == 1)  {
    //this->domain->writeVectorToFile("myResidual", globalIt/10, globalIt, R(0), &one);

    DistVec<double> rmag(R(0).info());
    for (int iSub = 0; iSub < rmag.numLocSub(); ++iSub) {

      for (int i = 0; i < rmag.subSize(iSub); ++i) {

	rmag(iSub)[i] = 0.0;
	for (int k = 0; k < dim; ++k)
	  rmag(iSub)[i] += R(0)(iSub)[i][k]*R(0)(iSub)[i][k];
      }
    }
    
    double max_res = rmag.max();
    
    for (int iSub = 0; iSub < rmag.numLocSub(); ++iSub) {
      
      for (int i = 0; i < rmag.subSize(iSub); ++i) {

	if (rmag(iSub)[i] == max_res) {

	  std::cout << "Maximum residual = " << rmag(iSub)[i] << " at node " <<
	    this->domain->getSubDomain()[iSub]->getNodeMap()[i]+1 << "; " << std::endl;
	}
      }
    }
  }
}

template <int dim>
void MultiGridCoupledTsDesc<dim>::
smooth(int lvl, MultiGridDistSVec<double,dim>& x,
       DistSVec<double,dim>& f,int steps) {

  int i;
  for (i = 0; i < steps; ++i) {

    this->varFcn->conservativeToPrimitive(x(lvl), V(lvl));

    mgSpaceOp->updateStateVectors(lvl,x);

    mgSpaceOp->computeTimeStep(lvl,this->data->cfl*pow(0.75,lvl),
                               V);
 
    //mgSpaceOp->computeResidual(lvl, x, V, res);
    if (i == 0)
      mgSpaceOp->computeJacobian(lvl, x, V, *mgMvp);
    R(lvl) = f-res(lvl);
    //R(lvl) = -1.0*res(lvl);
    if (smoothWithGMRES)
      mgKspSolver->solve(lvl, *mgMvp, R, dx);  
    else {
      if (i == 0)
        smoothingMatrices->acquire(lvl, *mgMvp);
      smoothingMatrices->apply(lvl, dx, R); 
    }
    
    x(lvl) += dx(lvl);

    //x(lvl) += 1e-4*R(lvl);
   
    // char fn[64];
    // if ( i % 100 == 0) {
    //   double rnorm = R(lvl).norm();
      
    //   int rnk;
    //   MPI_Comm_rank(MPI_COMM_WORLD, &rnk);
    
    //   if (rnk == 0)
    // 	std::cout << "i = " << i << " r = " << rnorm  << std::endl;
    //   sprintf(fn, "solutionresf%i_%i", globalIt,i);
    //   pKernel->getLevel(lvl)->writePVTUSolutionFile(fn,
    // 						  f); 
      
    //   sprintf(fn, "solutionres%i_%i", globalIt,i);
    //   pKernel->getLevel(lvl)->writePVTUSolutionFile(fn,
    // 						  res(lvl)); 
    
    //   sprintf(fn, "solutionV_%i",i);
    //          pKernel->getLevel(lvl)->writePVTUSolutionFile(fn,
    //     						  V(lvl));

    // }

    this->varFcn->conservativeToPrimitive(x(lvl), V(lvl));
    pKernel->fixNegativeValues(lvl,V(lvl), x(lvl), dx(lvl), f,Forig(lvl), this->varFcn,
                               mgSpaceOp->getOperator(lvl));
    mgSpaceOp->computeResidual(lvl, x, V, res, false);
    R(lvl) = f-res(lvl);
  }
}

template <int dim>
void MultiGridCoupledTsDesc<dim>::cycle(int lvl, DistSVec<double,dim>& f,
                                        MultiGridDistSVec<double,dim>& x) {

  if (lvl == 0) { 
    smooth0(x(lvl), numSmooths_pre[0]);
    mgSpaceOp->setupBcs(this->getSpaceOperator()->getDistBcData());
  }
  else
    smooth(lvl,x, f,  numSmooths_pre[lvl]);

  if (lvl < pKernel->numLevels()-1) {

    pKernel->Restrict(lvl+1, x(lvl), U(lvl+1));
    pKernel->Restrict(lvl+1, R(lvl), R(lvl+1));
    Uold(lvl+1) = U(lvl+1);
    
    this->varFcn->conservativeToPrimitive(U(lvl+1), V(lvl+1));
    //pKernel->fixNegativeValues(lvl+1,V(lvl+1), x(lvl+1), dx(lvl+1), R(lvl+1),Forig(lvl+1), this->varFcn);
    mgSpaceOp->computeResidual(lvl+1, U, V, res, false);
    //pKernel->fixNegativeValues(lvl+1,V(lvl+1), x(lvl+1), dx(lvl+1), F(lvl+1),Forig(lvl+1), this->varFcn);
    if (lvl == 0) {

    //  pKernel->getLevel(lvl+1)->writePVTUSolutionFile("myR",R(lvl+1));
    }
    pKernel->applyFixes(lvl+1, R(lvl+1));
    F(lvl+1) = res(lvl+1) + R(lvl+1)*restrict_relax_factor;
    for (int i = 0; i < mc; ++i)
      cycle(lvl+1, F(lvl+1), U);
    
    pKernel->Prolong(lvl+1, Uold(lvl+1), U(lvl+1), x(lvl), x(lvl), prolong_relax_factor, this->varFcn);
  }
  
  if (lvl == 0) 
    smooth0(x(lvl), numSmooths_post[0]);
  else
    smooth(lvl,x, f,  numSmooths_post[lvl]);
}

template <int dim>
void MultiGridCoupledTsDesc<dim>::cycle(DistSVec<double,dim>& x) {

  F(0) = 0.0;

  U(0) = x;

  cycle(0, F(0), U);  

  x = U(0);
}

template class MultiGridCoupledTsDesc<5>;
template class MultiGridCoupledTsDesc<6>;
template class MultiGridCoupledTsDesc<7>;

