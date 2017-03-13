#include <VectorSet.h>
#include <LevelSet/LevelSetStructure.h>
#include <MultiGridSegTsDesc.h>

template <int dim,int neq1,int neq2>
MultiGridSegTsDesc<dim,neq1,neq2>::
MultiGridSegTsDesc(IoData & iod, GeoSource & gs,  Domain * dom) :
  ImplicitSegTsDesc<dim,neq1,neq2>(iod,gs,dom) {

  memset(numSmooths_post,0,sizeof(numSmooths_post));
  memset(numSmooths_pre,0,sizeof(numSmooths_pre));
  numSmooths_pre[0] = 1;
  numSmooths_pre[1] = 2;
  numSmooths_pre[2] = 3;
  numSmooths_pre[3] = 4;
  numSmooths_pre[4] = 5;
  numSmooths_pre[5] = 6;
  
  numSmooths_post[0] = 0;
  numSmooths_post[1] = 0;
  numSmooths_post[2] = 0;
  numSmooths_post[3] = 0;
  numSmooths_post[4] = 0;
  numSmooths_post[5] = 0;
  
  prolong_relax_factor = iod.mg.prolong_relax_factor;
  restrict_relax_factor = iod.mg.restrict_relax_factor;

  addTurbulenceTerms = iod.mg.addTurbulenceTerms;

  if (iod.mg.cycle_scheme == MultiGridData::VCYCLE)
    mc = 1;
  else if (iod.mg.cycle_scheme == MultiGridData::WCYCLE)
    mc = 2;     

  smoothWithGMRES = (iod.mg.mg_smoother == MultiGridData::MGGMRES);

  globalIt = 0;

  mgMvp1 = NULL;
  pKernel = NULL;
  mgSpaceOp = NULL;
  mgKspSolver1 = NULL;
  mgKspSolver2 = NULL;
  smoothingMatrices1 = NULL;
  smoothingMatrices2 = NULL;

  TsDesc<dim>::isMultigridTsDesc = true;
}

template <int dim,int neq1,int neq2>
MultiGridSegTsDesc<dim,neq1,neq2>::
~MultiGridSegTsDesc() {

  if (mgMvp1)
    delete mgMvp1;
  if (pKernel)
    delete pKernel;
  if (mgSpaceOp)
    delete mgSpaceOp;
  if (mgKspSolver1)
    delete mgKspSolver1;
  if (mgKspSolver2)
    delete mgKspSolver2;
  if (smoothingMatrices1)
    delete smoothingMatrices1;
  if (smoothingMatrices2)
    delete smoothingMatrices2;
}


template <int dim,int neq1,int neq2>
void MultiGridSegTsDesc<dim,neq1,neq2>::
setupTimeStepping(DistSVec<double,dim> *U0, IoData &iod) {

  ImplicitSegTsDesc<dim,neq1,neq2>::setupTimeStepping(U0,iod);

  pKernel = new MultiGridKernel<double>(this->domain, *this->geoState,
                                        iod,iod.mg.num_multigrid_levels);

  pKernel->initialize(dim,neq1,neq2);

  mgSpaceOp = 
    new MultiGridSpaceOperator<double,dim>(iod, this->domain, this->spaceOp, pKernel,
                                           this->spaceOp1,this->spaceOp2);

  mgMvp1 = new MultiGridMvpMatrix<double,neq1>(this->domain,pKernel);
  mgMvp2 = new MultiGridMvpMatrix<double,neq2>(this->domain,pKernel);

  mgKspSolver1 = new MultiGridKspSolver<double,neq1,double>(this->domain, iod.ts.implicit.newton.ksp.ns,
                                                          pKernel);

  mgKspSolver2 = new MultiGridKspSolver<double,neq2,double>(this->domain, iod.ts.implicit.newton.ksp.ns,
                                                          pKernel);

  pKernel->setUseVolumeWeightedAverage(iod.mg.restrictMethod == MultiGridData::VOLUME_WEIGHTED);

  if (!smoothWithGMRES) {
    typename MultiGridSmoothingMatrix<double,neq1>::SmoothingMode s;
    if (iod.mg.mg_smoother == MultiGridData::MGRAS)
      s = MultiGridSmoothingMatrix<double,neq1>::RAS; 
    if (iod.mg.mg_smoother == MultiGridData::MGJACOBI)
      s = MultiGridSmoothingMatrix<double,neq1>::BlockJacobi;
    if (iod.mg.mg_smoother == MultiGridData::MGLINEJACOBI)
      s = MultiGridSmoothingMatrix<double,neq1>::LineJacobi;
    smoothingMatrices1 = new MultiGridSmoothingMatrices<double,neq1>(pKernel, s);
    
    typename MultiGridSmoothingMatrix<double,neq2>::SmoothingMode s2;
    if (iod.mg.mg_smoother == MultiGridData::MGRAS)
      s2 = MultiGridSmoothingMatrix<double,neq2>::RAS; 
    if (iod.mg.mg_smoother == MultiGridData::MGJACOBI)
      s2 = MultiGridSmoothingMatrix<double,neq2>::BlockJacobi;
    if (iod.mg.mg_smoother == MultiGridData::MGLINEJACOBI)
      s2 = MultiGridSmoothingMatrix<double,neq2>::LineJacobi;
    smoothingMatrices2 = new MultiGridSmoothingMatrices<double,neq2>(pKernel,s2);

  }

  V.init(pKernel);
  res.init(pKernel);
  R.init(pKernel);
  R1.init(pKernel);
  R2.init(pKernel);
  F.init(pKernel);
  Forig.init(pKernel);
  U.init(pKernel);
  Uold.init(pKernel); 
  dx.init(pKernel);
  dx1.init(pKernel);
  dx2.init(pKernel);

  update_tmp.init(pKernel);  
}

template <int dim,int neq1,int neq2>
void MultiGridSegTsDesc<dim,neq1,neq2>::
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
    else {
      smoothingMatrices1->acquire(*this->GetJacobian1());
      smoothingMatrices2->acquire(*this->GetJacobian2());
    }
  
    R(0) *= -1.0;
    if (smoothWithGMRES)
      this->solveLinearSystem(0, R(0),dx(0));
    else {
      R(0).split(R1(0), R2(0));
      smoothingMatrices1->apply(0, dx1, R1);
      smoothingMatrices2->apply(0, dx2, R2);
      dx(0).merge(dx1(0),dx2(0));
    }
    
    x += dx(0);
    
    this->updateStateVectors(x, 0);
    this->monitorConvergence(0, x);
    R(0) = -1.0*this->getCurrentResidual();
    
  }
  //if (i == 0) {
  //  this->monitorConvergence(0, x);
  //  R(0) = -1.0*this->getCurrentResidual();
  //}

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

template <int dim,int neq1,int neq2>
void MultiGridSegTsDesc<dim,neq1,neq2>::
smooth(int lvl, MultiGridDistSVec<double,dim>& x,
       DistSVec<double,dim>& f,int steps, bool postsmooth) {

  int i;
  double norm;
    int rnk;
    MPI_Comm_rank(MPI_COMM_WORLD,&rnk);
    
  for (i = 0; i < /*100000000*/steps; ++i) {

    this->varFcn->conservativeToPrimitive(x(lvl), V(lvl));

    mgSpaceOp->updateStateVectors(lvl,x);

    mgSpaceOp->computeTimeStep(lvl,this->data->cfl*pow(0.75,lvl),
                               V);
 
    if (i == 0 && postsmooth)
      mgSpaceOp->computeResidual(lvl, x, V, res);

    if (i == 0) {
      mgSpaceOp->computeJacobian(lvl, x, V, *mgMvp1);
      if (addTurbulenceTerms)
	mgSpaceOp->computeTurbulentJacobian(lvl, x, V, *mgMvp2);
    }
    R(lvl) = f-1.0*res(lvl);
    //R(lvl) = -1.0*res(lvl);
    R(lvl).split(R1(lvl), R2(lvl));

    norm = R1(lvl).norm();
    if (rnk == 0)
      std::cout << "i = " << i << " Rnorm = " << norm << std::endl;

    if (smoothWithGMRES) {
      mgKspSolver1->solve(lvl, *mgMvp1, R1, dx1);
      if (addTurbulenceTerms)
	mgKspSolver2->solve(lvl, *mgMvp2, R2, dx2);
    }
    else {
      if (i == 0) {
	smoothingMatrices1->acquire(lvl, *mgMvp1);
	if (addTurbulenceTerms)
	  smoothingMatrices2->acquire(lvl, *mgMvp2);
      }
      smoothingMatrices1->apply(lvl, dx1, R1);
      if (addTurbulenceTerms)
	smoothingMatrices2->apply(lvl, dx2, R2); 
    }
    if (!addTurbulenceTerms)
      dx2(lvl) = 0.0;
    dx(lvl).merge(dx1(lvl), dx2(lvl));
    
    x(lvl) += dx(lvl);

    
    norm = dx(lvl).norm();
    if (rnk == 0)
      std::cout << "i = " << i << " dx norm = " << norm << std::endl;
   
    this->varFcn->conservativeToPrimitive(x(lvl), V(lvl));

    pKernel->fixNegativeValues(lvl,V(lvl), x(lvl), dx(lvl), f,Forig(lvl), this->varFcn,
                               mgSpaceOp->getOperator(lvl));

    mgSpaceOp->computeResidual(lvl, x, V, res, false);
    R(lvl) = f-res(lvl);

    /*
    if (lvl == 1 && i%1000==0) {
      pKernel->getLevel(lvl)->writePVTUSolutionFile("r.sol",x(lvl));
      MPI_Barrier(MPI_COMM_WORLD);
    }
    */
  }
 
/*
  if (lvl == 1) {
    pKernel->getLevel(lvl)->writePVTUSolutionFile("r.sol",x(lvl));
    MPI_Barrier(MPI_COMM_WORLD);
    exit(0);
  }
*/
}

template <int dim,int neq1,int neq2>
void MultiGridSegTsDesc<dim,neq1,neq2>::cycle(int lvl, DistSVec<double,dim>& f,
                                    MultiGridDistSVec<double,dim>& x) {

  int rnk;
    MPI_Comm_rank(MPI_COMM_WORLD,&rnk);
    double norm;
  if (lvl == 0) 
    for (int iSub = 0; iSub < this->domain->getNumLocSub(); ++iSub) {

      double Vloc[dim];
      for (int l = 0; l < x(lvl)(iSub).size(); ++l) {
        if (this->domain->getSubDomain()[iSub]->getNodeMap()[l]+1 == 9 ||
	    this->domain->getSubDomain()[iSub]->getNodeMap()[l]+1 == 127182 ||
            this->domain->getSubDomain()[iSub]->getNodeMap()[l]+1 == 1892729 || 
            this->domain->getSubDomain()[iSub]->getNodeMap()[l]+1 == 11148134) {   

          this->varFcn->conservativeToPrimitive(x(lvl)(iSub)[l],Vloc);
          std::cout << "V[" << this->domain->getSubDomain()[iSub]->getNodeMap()[l]+1  << "] = ";
          for (int k = 0; k < dim; ++k)
            std::cout << Vloc[k] << " ";
          std::cout << std::endl;
        }
     }
    }
  if (lvl == 0) { 
    smooth0(x(lvl), numSmooths_pre[0]);
    mgSpaceOp->setupBcs(this->getSpaceOperator()->getDistBcData());
  }
  else
    smooth(lvl,x, f,  numSmooths_pre[lvl],false);
  
  // Check for negative turbulence values
  
  if (lvl == 0)  {
    for (int iSub = 0; iSub < this->domain->getNumLocSub(); ++iSub) {

      double Vloc[dim];
      for (int l = 0; l < x(lvl)(iSub).size(); ++l) {

	if (x(lvl)(iSub)[l][5] < 1.0e-10)
	  x(lvl)(iSub)[l][5] = 1.0e-10;
      }
    }
  }

  if (lvl == 0) 
    for (int iSub = 0; iSub < this->domain->getNumLocSub(); ++iSub) {

      double Vloc[dim];
      for (int l = 0; l < x(lvl)(iSub).size(); ++l) {
        if (this->domain->getSubDomain()[iSub]->getNodeMap()[l]+1 == 9 ||
	    this->domain->getSubDomain()[iSub]->getNodeMap()[l]+1 == 127182 ||
            this->domain->getSubDomain()[iSub]->getNodeMap()[l]+1 == 1892729 || 
            this->domain->getSubDomain()[iSub]->getNodeMap()[l]+1 == 11148134) {   

          this->varFcn->conservativeToPrimitive(x(lvl)(iSub)[l],Vloc);
          std::cout << "V[" << this->domain->getSubDomain()[iSub]->getNodeMap()[l]+1  << "] = ";
          for (int k = 0; k < dim; ++k)
            std::cout << Vloc[k] << " ";
          std::cout << std::endl;
        }
     }
    }
  if (lvl < pKernel->numLevels()-1/* && globalIt > 50*/) {
    /*
    for (int iSub = 0; iSub < this->domain->getNumLocSub(); ++iSub) {

      double Vloc[dim];
      for (int l = 0; l < x(lvl)(iSub).size(); ++l) {
        if (this->domain->getSubDomain()[iSub]->getNodeMap()[l] == 220559) {   

          this->varFcn->conservativeToPrimitive(x(lvl)(iSub)[l],Vloc);
          std::cout << "V[220559] = ";
          for (int k = 0; k < dim; ++k)
            std::cout << Vloc[k] << " ";
          std::cout << std::endl;
        }
     }
   }
    */
    pKernel->Restrict(lvl+1, x(lvl), U(lvl+1));
    //this->domain->getCommunicator()->fprintf(stderr,"Restricting residual...\n");
    //fflush(stderr);
    //MPI_Barrier(MPI_COMM_WORLD);
    pKernel->Restrict(lvl+1, R(lvl), R(lvl+1));
    Uold(lvl+1) = U(lvl+1);
    
    this->varFcn->conservativeToPrimitive(U(lvl+1), V(lvl+1));
    /*
    int nodes[] = {7713 ,7718,7924, 7925 ,7927 , 8029,8032, 8034, 310130 };
    for (int iSub = 0; iSub < this->domain->getNumLocSub(); ++iSub) {

      for (int i = 0; i < pKernel->getLevel(lvl+1)->getNodeDistInfo().subSize(iSub); ++i) {

        for (int k = 0; k < 9; ++k) {
          if (pKernel->getLevel(lvl+1)->getMgSubDomains()[iSub].locToGlobMap[i] == 
              nodes[k]) {
            
            std::cout << "V[" <<nodes[k] <<  "] = ";
            for (int l = 0; l < dim; ++l)
              std::cout << V(lvl+1)(iSub)[i][l] << " ";
            std::cout << std::endl;
          }
        }
      }
    }
    */
    mgSpaceOp->computeResidual(lvl+1, U, V, res, false);

    pKernel->applyFixes(lvl+1, R(lvl+1));

    //pKernel->fixNegativeValues(lvl+1,V(lvl+1), U(lvl+1), dx(lvl+1), F(lvl+1), this->varFcn);
    /*if (lvl == 0 && globalIt % 25 == 0) {

      pKernel->getLevel(lvl+1)->writePVTUSolutionFile("myR",R(lvl+1));
    }*/
    //Forig(lvl+1) = F(lvl+1);
    if (rnk == 0)
      std::cout << "Current Relaxation Factor = " << restrict_relax_factor*std::min((double)(globalIt-50)/200,1.0) << std::endl;

    double loc_fac = restrict_relax_factor*std::min((double)(globalIt-50)/200,1.0);
    F(lvl+1) = res(lvl+1) + R(lvl+1)*loc_fac;
    Forig(lvl+1) = res(lvl+1);

    for (int i = 0; i < mc; ++i)
      cycle(lvl+1, F(lvl+1), U);
    
    update_tmp(lvl) = 0.0;
    loc_fac = prolong_relax_factor*std::min((double)(globalIt-50)/200,1.0);
    pKernel->Prolong(lvl+1, Uold(lvl+1), U(lvl+1), update_tmp(lvl),x(lvl), loc_fac, this->varFcn);

    pKernel->applyFixes(lvl,update_tmp(lvl));
    x(lvl) += update_tmp(lvl);
    this->varFcn->conservativeToPrimitive(x(lvl), V(lvl));
    pKernel->fixNegativeValues(lvl,V(lvl), x(lvl), update_tmp(lvl), F(lvl),F(lvl),
			       this->varFcn,
			       mgSpaceOp->getOperator(lvl));
  }
 
  if (lvl == 0) 
    for (int iSub = 0; iSub < this->domain->getNumLocSub(); ++iSub) {

      double Vloc[dim];
      for (int l = 0; l < x(lvl)(iSub).size(); ++l) {
        if (this->domain->getSubDomain()[iSub]->getNodeMap()[l]+1 == 9 ||
	    this->domain->getSubDomain()[iSub]->getNodeMap()[l]+1 == 127182 ||
            this->domain->getSubDomain()[iSub]->getNodeMap()[l]+1 == 1892729 || 
            this->domain->getSubDomain()[iSub]->getNodeMap()[l]+1 == 11148134) {   

          this->varFcn->conservativeToPrimitive(x(lvl)(iSub)[l],Vloc);
          std::cout << "V[" << this->domain->getSubDomain()[iSub]->getNodeMap()[l]+1  << "] = ";
          for (int k = 0; k < dim; ++k)
            std::cout << Vloc[k] << " ";
          std::cout << std::endl;
        }
     }
    }
  if (lvl == 0) 
    smooth0(x(lvl), numSmooths_post[0]);
  else
    smooth(lvl,x, f,  numSmooths_post[lvl], true);
}

template <int dim,int neq1,int neq2>
void MultiGridSegTsDesc<dim,neq1,neq2>::cycle(DistSVec<double,dim>& x) {

  F(0) = 0.0;

  U(0) = x;
  
  cycle(0, F(0), U);

  x = U(0); 
}

template class MultiGridSegTsDesc<6,5,1>;
template class MultiGridSegTsDesc<7,5,2>;

