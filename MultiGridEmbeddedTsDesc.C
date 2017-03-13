#include <VectorSet.h>
#include <LevelSet/LevelSetStructure.h>
#include <MultiGridEmbeddedTsDesc.h>

#include <mpi.h>

template <int dim>
MultiGridEmbeddedTsDesc<dim>::
MultiGridEmbeddedTsDesc(IoData & iod, GeoSource & gs,  Domain * dom) :
  ImplicitEmbeddedCoupledTsDesc<dim>(iod,gs,dom) {

  memset(numSmooths_pre,0,sizeof(numSmooths_pre));
  memset(numSmooths_post,0,sizeof(numSmooths_post));
  numSmooths_pre[0] = 1;
  numSmooths_pre[1] = 2;
  numSmooths_pre[2] = 3;
  numSmooths_pre[3] = 4;
  numSmooths_pre[4] = 5;
  numSmooths_pre[5] = 6;
  numSmooths_pre[6] = 7;
  numSmooths_pre[7] = 8;
  numSmooths_pre[8] = 9;
  
  numSmooths_post[0] = 0;
  numSmooths_post[1] = 0;
  numSmooths_post[2] = 0;
  numSmooths_post[3] = 0;
  numSmooths_post[4] = 0;
  numSmooths_post[5] = 0;
  numSmooths_post[6] = 0;
  numSmooths_post[7] = 0;
  numSmooths_post[8] = 0;
  
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

  mgLSS  = NULL;

  TsDesc<dim>::isMultigridTsDesc = true;

  densityMin = iod.mg.densityMin/iod.ref.rv.density;
  densityMax = iod.mg.densityMax/iod.ref.rv.density;

}

template <int dim>
MultiGridEmbeddedTsDesc<dim>::
~MultiGridEmbeddedTsDesc() {

  if (mgLSS) {

    for (int i = 1; i < pKernel->numLevels(); ++i)
      delete mgLSS[i];
    
    delete [] mgLSS;
  }

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
void MultiGridEmbeddedTsDesc<dim>::
setupTimeStepping(DistSVec<double,dim> *U0, IoData &iod) {

  EmbeddedTsDesc<dim>::setupTimeStepping(U0, iod);

  U_smoothed = new DistSVec<double,dim>(*U0);

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

  mgLSS = new DistMultiGridLevelSetStructure*[pKernel->numLevels()];
  mgLSS[0] = NULL;

  for (int i = 1; i < pKernel->numLevels(); ++i) {
    DistLevelSetStructure* parent = (i == 1 ? this->distLSS : mgLSS[i-1]);
    mgLSS[i] = new DistMultiGridLevelSetStructure(iod, this->domain->getCommunicator(),
						  parent, pKernel->getLevel(i));

    mgLSS[i]->initialize(this->domain,pKernel->getLevel(i)->getXn(),
			 pKernel->getLevel(i)->getXn(), iod);
  }

}

template <int dim>
void MultiGridEmbeddedTsDesc<dim>::
smooth0(DistSVec<double,dim>& x,int steps) {

  int i;
  double dummy = 0.0;
  this->updateStateVectors(x, 0);

  int rnk;
  MPI_Comm_rank(MPI_COMM_WORLD, &rnk);


  for (i = 0; i < steps; ++i) {

    this->domain->checkSolution(this->varFcn, x);

    this->computeTimeStep(globalIt,&dummy, x);
    ++globalIt;
    this->computeFunction(0, x, R(0));
    this->computeJacobian(0, x, R(0));
    if (smoothWithGMRES)
      this->setOperators(x);
    else
      smoothingMatrices->acquire( *this->GetJacobian());
  
    R(0) *= -1.0;

    double mag = R(0).norm();
    if (rnk == 0) {
      std::cout << "norm = " << mag << std::endl;
    }

    if (smoothWithGMRES)
      this->solveLinearSystem(0, R(0),dx(0));
    else
      smoothingMatrices->apply(0, dx, R);
    
    x += dx(0);
    
    this->updateStateVectors(x, 0);
    this->monitorConvergence(0, x);
    R(0) = -1.0*this->getCurrentResidual();

    *U_smoothed = x;

  }

  double mag = R(0).norm();
  if (rnk == 0) {
    std::cout << "norm = " << mag << std::endl;
  }

//  if (i == 0) {
//    this->monitorConvergence(0, x);
//    R(0) = -1.0*this->getCurrentResidual();
//  }
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
void MultiGridEmbeddedTsDesc<dim>::
smooth(int lvl, MultiGridDistSVec<double,dim>& x,
       DistSVec<double,dim>& f,int steps, 
       bool postsmooth) {

  int i;
  //if (globalIt > 200 && lvl == 3 && !postsmooth)
  //  steps += std::min((globalIt-200)/10,30);
      

  int rnk;
  MPI_Comm_rank(MPI_COMM_WORLD, &rnk);

//  for (i = 0; i < 100000/*steps*/; ++i) {
  for (i = 0; i < steps; ++i) {

    
    //if (rnk == 0)
    //  std::cout << "i = " << i << std::endl;
    
    dx(lvl) = 0.0;
    pKernel->fixNegativeValues(lvl,V(lvl), x(lvl), dx(lvl), f,Forig(lvl), this->varFcn,
                               mgSpaceOp->getOperator(lvl));

    this->varFcn->conservativeToPrimitive(x(lvl), V(lvl));

    mgSpaceOp->updateStateVectors(lvl,x);

    mgSpaceOp->computeTimeStep(lvl,this->data->cfl*pow(0.75,lvl),
                               V);
//    mgSpaceOp->computeTimeStep(lvl,std::min<double>(5.0 + 0.5*(i+1),1000.0),
//                               V);
 
    if (i == 0 && postsmooth)
      mgSpaceOp->computeResidualEmbedded(*this->riemann,lvl, x, V, res , mgLSS[lvl], false );

//    mgSpaceOp->add_dAW_dtEmbedded(res, mgLSS[lvl]);

    if (i %25 == 0)
      mgSpaceOp->computeJacobianEmbedded(*this->riemann,lvl, x, V, *mgMvp,  mgLSS[lvl] );
    
  //  double r = res(lvl).norm();
  //  if (rnk == 0)
  //    std::cout << r << std::endl;
    

    char fn[64];

    R(lvl) = f-1.0*res(lvl);

    for (int iSub = 0; iSub < R(lvl).numLocSub(); ++iSub) {

      for (int j = 0; j < R(lvl).subSize(iSub); ++j) {

	if (!(*mgLSS[lvl])(iSub).isActive(0.0,j)) {

	  memset(R(lvl)(iSub)[j], 0, sizeof(double)*dim);
	}
      }
    }


    double mag = R(lvl).norm();
    if (rnk == 0)
      std::cout << "i = " <<  i << " res = " << mag << std::endl;


    
    if ( /*globalIt % 25 == 0 &&*/ i == 0) {
      //sprintf(fn, "solution_dz_%i_%i", globalIt,lvl);
      //  pKernel->getLevel(lvl)->writePVTUSolutionFile(fn, 
      //					    *mgSpaceOp->getOperator(lvl)->getDerivative(2)); 
      /*
      sprintf(fn, "solutionres_%i_%i", lvl, globalIt);
      pKernel->getLevel(lvl)->writePVTUSolutionFile(fn,
      						  R(lvl)); 
   
      sprintf(fn, "solutionV_%i_%i",lvl,globalIt);
              pKernel->getLevel(lvl)->writePVTUSolutionFile(fn,
       						  V(lvl));
      */
    }

    //if (lvl == 3) {

      // double mag = R(lvl).norm();
      // if (rnk == 0)
      // 	std::cout << "i = " << i << " res = " << mag << std::endl;
      //}
    // R(lvl) = -1.0*res(lvl);
      /*
    V(lvl) = 1;
    for (int iSub = 0; iSub < R(lvl).numLocSub(); ++iSub) {

      for (int i = 0; i < R(lvl).subSize(iSub); ++i) {

	if (!(*mgLSS[lvl])(iSub).isActive(0.0,i)) {
	  memset(R(lvl)[i],0, sizeof(double)*dim);
	  memset(V(lvl)(iSub)[i],0, sizeof(double)*dim);
	}
      }
    }
      */
    
//    sprintf(fn, "isActive%i", i);
//    pKernel->getLevel(lvl)->writePVTUSolutionFile(fn,
//						  V(lvl));


    dx(lvl) = 0.0;
    if (smoothWithGMRES)
      mgKspSolver->solve(lvl, *mgMvp, R, dx);  
    else {
      if (i % 25 == 0)
        smoothingMatrices->acquire(lvl, *mgMvp);
      smoothingMatrices->apply(lvl, dx, R); 
    }

    //dx(lvl) *= 0.5;

    //dx(lvl) = 0.0;
    /*
    for (int iSub = 0; iSub < R(lvl).numLocSub(); ++iSub) {

      for (int i = 0; i < R(lvl).subSize(iSub); ++i) {

	if (!(*mgLSS[lvl])(iSub).isActive(0.0,i)) {

	  memset(dx(lvl)(iSub)[i], 0, sizeof(double)*dim);
	}
      }
    }

    std::cout << dx(lvl).norm() << std::endl;
    */

    mag = dx(lvl).norm();
    if (rnk == 0)
      std::cout << "res = " << mag << std::endl;

    x(lvl) += dx(lvl); 

    this->varFcn->conservativeToPrimitive(x(lvl), V(lvl));

    /*
    pKernel->getLevel(lvl)->writePVTUSolutionFile("solutionV2",
    V(lvl));*/
    /*
    MPI_Barrier(MPI_COMM_WORLD);

    if (globalIt == 4)
      exit(-1);
    */
    pKernel->fixNegativeValues(lvl,V(lvl), x(lvl), dx(lvl), f,Forig(lvl), this->varFcn,
                               mgSpaceOp->getOperator(lvl));
    mgSpaceOp->computeResidualEmbedded(*this->riemann,lvl, x, V, res,mgLSS[lvl], false);
    R(lvl) = f-res(lvl);
  }

  // double mag = R(lvl).norm();
  // if (rnk == 0)
  //   std::cout << " res = " << mag << std::endl;
}

template <int dim>
void MultiGridEmbeddedTsDesc<dim>::cycle(int lvl, DistSVec<double,dim>& f,
                                        MultiGridDistSVec<double,dim>& x) {

  int rnk;
  MPI_Comm_rank(MPI_COMM_WORLD, &rnk);

  // double fnorm = f.norm();
  // if (rnk == 0)
  //   std::cout << fnorm << std::endl;

  if (lvl == 0) { 
    smooth0(x(lvl), numSmooths_pre[0]);
    mgSpaceOp->setupBcs(this->getSpaceOperator()->getDistBcData());
  }
  else
    smooth(lvl,x, f,  numSmooths_pre[lvl],false);

  // fnorm = f.norm();
  // if (rnk == 0)
  //   std::cout << fnorm << std::endl;
  
  if ( lvl < pKernel->numLevels()-1) {

    pKernel->Restrict(lvl+1, x(lvl), U(lvl+1));
    pKernel->Restrict(lvl+1, R(lvl), R(lvl+1));
    Uold(lvl+1) = U(lvl+1);
    
    this->varFcn->conservativeToPrimitive(U(lvl+1), V(lvl+1));
    //pKernel->fixNegativeValues(lvl+1,V(lvl+1), x(lvl+1), dx(lvl+1), R(lvl+1),Forig(lvl+1), this->varFcn);
    mgSpaceOp->computeResidualEmbedded(*this->riemann,lvl+1, U, V, res, mgLSS[lvl+1], false);
    //pKernel->fixNegativeValues(lvl+1,V(lvl+1), x(lvl+1), dx(lvl+1), F(lvl+1),Forig(lvl+1), this->varFcn);
    if (lvl == 0 && globalIt % 100 == 0) {

    //  pKernel->getLevel(lvl)->writePVTUSolutionFile("myR",R(lvl));
    }
    pKernel->applyFixes(lvl+1, R(lvl+1));

    char fn[64];
    
    for (int iSub = 0; iSub < R(lvl+1).numLocSub(); ++iSub) {
      
      for (int i = 0; i < R(lvl+1).subSize(iSub); ++i) {

	if (V(lvl+1)(iSub)[i][0] <= densityMin ||
	    V(lvl+1)(iSub)[i][0] >= densityMax) {
	  for (int k = 0; k < dim; ++k) {
	    R(lvl+1)(iSub)[i][k] = 0.0;
	  }
	}
      }
    }
    
    if (globalIt % 100 == 0) {
      sprintf(fn, "solutionresF%i", globalIt);
     // pKernel->getLevel(lvl+1)->writePVTUSolutionFile(fn,
//						      F(lvl+1)); 

    }
    

   
    /*
    for (int iSub = 0; iSub < R(lvl+1).numLocSub(); ++iSub) {
      
      for (int i = 0; i < R(lvl+1).subSize(iSub); ++i) {

	for (int k = 0; k < dim; ++k) {
	  if (fabs(F(lvl+1)(iSub)[i][k]) > 1.0e-8 && 
	      R(lvl+1)(iSub)[i][k] / F(lvl+1)(iSub)[i][k] > 0.0)
	    F(lvl+1)(iSub)[i][k] += R(lvl+1)(iSub)[i][k]*restrict_relax_factor;
	}
      }
    }
    */
    F(lvl+1) = res(lvl+1) + R(lvl+1)*restrict_relax_factor;
    /*
    for (int iSub = 0; iSub < R(lvl+1).numLocSub(); ++iSub) {
      
      for (int i = 0; i < R(lvl+1).subSize(iSub); ++i) {

	if (!(*mgLSS[lvl+1])(iSub).isActive(0.0,i)) {
	  for (int k = 0; k < dim; ++k) {
	    F(lvl+1)(iSub)[i][k] = 0.0;
	  }
	}
      }
    }
    */
    for (int i = 0; i < mc; ++i)
      cycle(lvl+1, F(lvl+1), U);

    DistLevelSetStructure* parent = (lvl == 0 ? this->distLSS : mgLSS[lvl]);
    pKernel->ExtrapolateProlongation(lvl+1,Uold(lvl+1), U(lvl+1), 
    				     mgLSS[lvl+1], parent);

    pKernel->Prolong(lvl+1, Uold(lvl+1), U(lvl+1), x(lvl), x(lvl), prolong_relax_factor,
                     this->varFcn,
		     mgLSS[lvl+1], parent);
  }

  // fnorm = f.norm();
  // if (rnk == 0)
  //   std::cout << fnorm << std::endl;

  if (lvl == 0) 
    smooth0(x(lvl), numSmooths_post[0]);
  else
    smooth(lvl,x, f,  numSmooths_post[lvl],true);

  // fnorm = f.norm();
  // if (rnk == 0)
  //   std::cout << fnorm << std::endl;
}

template <int dim>
void MultiGridEmbeddedTsDesc<dim>::cycle(DistSVec<double,dim>& x) {

  if (globalIt % 5 == 0) {
    ImplicitEmbeddedTsDesc<dim>::commonPart(x);

    for (int i = 1; i < pKernel->numLevels(); ++i) {

      mgLSS[i]->recompute(0.0,0.0,0.0, true);
    }
  }

  F(0) = 0.0;

  U(0) = x;

  cycle(0, F(0), U);  

  x = U(0);
}

template class MultiGridEmbeddedTsDesc<5>;
template class MultiGridEmbeddedTsDesc<6>;
template class MultiGridEmbeddedTsDesc<7>;

