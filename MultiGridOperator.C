/* MultiGridOperator.C

 */

#include <MultiGridOperator.h>
#include <MultiGridLevel.h>

template <class Scalar,int dim>
MultiGridOperator<Scalar,dim>::MultiGridOperator(MultiGridLevel<Scalar>* mg_lvl,
                                                 IoData& ioData, VarFcn* varFcn, 
                                                 Domain* domain, BcFcn* bcFcn,
                                                 BcFcn* bcFcn1,
                                                 BcFcn* bcFcn2) : bcFcn(bcFcn),
                                                 bcFcn1(bcFcn1),
                                                 bcFcn2(bcFcn2) {

  mgLevel = mg_lvl;
  myBcData = new DistBcData<dim>(ioData, varFcn, domain, &mg_lvl->getNodeDistInfo(),
                                 &mg_lvl->getInletNodeDistInfo(), &mg_lvl->getFaceDistInfo());
  timeState = new DistTimeState<dim>(ioData, NULL, varFcn, domain, mg_lvl->getNodeDistInfo(),
                                     NULL);

  boundaryState = new DistSVec<Scalar,dim>(mg_lvl->getAgglomFaceDistInfo());

  timeState->createSubStates();

  zero = new DistSVec<Scalar,dim>(mg_lvl->getNodeDistInfo());
  scalar_zero = new DistVec<Scalar>(mg_lvl->getNodeDistInfo());

  idti = new DistVec<Scalar>(mg_lvl->getNodeDistInfo()); 
  idtv = new DistVec<Scalar>(mg_lvl->getNodeDistInfo()); 

  DX[0] = new DistSVec<Scalar,dim>(mg_lvl->getNodeDistInfo());
  DX[1] = new DistSVec<Scalar,dim>(mg_lvl->getNodeDistInfo());
  DX[2] = new DistSVec<Scalar,dim>(mg_lvl->getNodeDistInfo());

  Wstarij = new DistSVec<Scalar,dim>(mg_lvl->getNodeDistInfo());
  Wstarji = new DistSVec<Scalar,dim>(mg_lvl->getNodeDistInfo());
  
  myVarFcn = varFcn;

  addViscousTerms = ioData.mg.addViscousTerms;

  timeState->disableRapidlyChangingValueCheck();
}

template <class Scalar,int dim>
MultiGridOperator<Scalar,dim>::~MultiGridOperator() {

  delete myBcData;
  delete timeState;
  delete boundaryState;
  delete zero;
  delete scalar_zero;
  delete idti;
  delete idtv;
  delete DX[0];
  delete DX[1];
  delete DX[2];

  delete Wstarij;
  delete Wstarji;
}

template<class Scalar,int dim>
template <class Scalar2,int neq>
void MultiGridOperator<Scalar,dim>::
computeJacobian(DistSVec<Scalar2,dim>& U, DistSVec<Scalar2,dim>& V,
		//                                         DistVec<Scalar2>& irey,
		FluxFcn **fluxFcn,
		FemEquationTerm* fet,
		DistMvpMatrix<Scalar2,neq>& matrices) {
  
  DistVec<Scalar2>& irey = *scalar_zero; 
#pragma omp parallel for
  for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {

    matrices(iSub) = 0.0;
  }

#pragma omp parallel for
  for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {

    mgLevel->getEdges()[iSub]->computeJacobianFiniteVolumeTerm(fluxFcn, (mgLevel->getGeoState())(iSub),irey(iSub),
                                                 mgLevel->getXn()(iSub),
                                                 mgLevel->getCtrlVol()(iSub), 
                                                 V(iSub), matrices(iSub));
    //mgLevel->getFaces()[iSub]->computeJacobianFiniteVolumeTerm(fluxFcn,(*myBcData)(iSub), (mgLevel->getGeoState())(iSub),
    //                                             V(iSub), matrices(iSub));
    mgLevel->getAgglomeratedFaces()[iSub]->computeJacobianFiniteVolumeTerm(fluxFcn,V(iSub), 
                                                 getBoundaryState()(iSub),
                                                 matrices(iSub));

    if (fet && addViscousTerms) {
      /*
      mgLevel->getEdges()[iSub]->computeJacobianThinLayerViscousFiniteVolumeTerm(
        NULL, myVarFcn, fet, mgLevel->getGeoState()(iSub), 
        mgLevel->getXn()(iSub),  V(iSub), mgLevel->getCtrlVol()(iSub), matrices(iSub));

      mgLevel->getAgglomeratedFaces()[iSub]->computeJacobianThinLayerViscousFiniteVolumeTerm(
         fet,myVarFcn, V(iSub), (*DX[0])(iSub), 
         (*DX[1])(iSub), (*DX[2])(iSub), mgLevel->getCtrlVol()(iSub),
         mgLevel->getGeoState()(iSub).getDistanceToWall(),
         (*myBcData)(iSub).getFaceStateVector(),matrices(iSub));
      */
    }
 
    Vec<double>& ctrlVol = mgLevel->getCtrlVol()(iSub); 
    for (int i=0; i<ctrlVol.size(); ++i) {
      Scalar voli = 1.0 / ctrlVol[i];
      Scalar2 *Aii = matrices(iSub).getElem_ii(i);
      for (int k=0; k<neq*neq; ++k)
        Aii[k] *= voli;
    }

  }

  mgLevel->assemble(matrices);

  if (timeState) {

    timeState->addToJacobian(mgLevel->getCtrlVol(), matrices, U);
  }

  applyBCsToJacobian(U, matrices); 
}

template<class Scalar,int dim>
template <class Scalar2,int neq>
void MultiGridOperator<Scalar,dim>::
computeJacobianEmbedded(DistExactRiemannSolver<dim>& riemann,
			DistSVec<Scalar2,dim>& U, 
			DistSVec<Scalar2,dim>& V,
			//                                         DistVec<Scalar2>& irey,
			FluxFcn **fluxFcn,
			FemEquationTerm* fet,
			DistMvpMatrix<Scalar2,neq>& matrices,	
			DistMultiGridLevelSetStructure* mgLSS) {
  
  DistVec<Scalar2>& irey = *scalar_zero; 
#pragma omp parallel for
  for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {

    matrices(iSub) = 0.0;
  }

#pragma omp parallel for
  for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {

    mgLevel->getEdges()[iSub]->computeJacobianFiniteVolumeTerm(riemann(iSub),
							       fluxFcn,
							       (mgLevel->getGeoState())(iSub),
							       mgLevel->getXn()(iSub),
							       V(iSub),
							       mgLevel->getCtrlVol()(iSub),
							       (*mgLSS)(iSub),
							       (*mgLSS).getStatus()(iSub) ,
							       1,
							       matrices(iSub),
							       irey(iSub));
    //mgLevel->getFaces()[iSub]->computeJacobianFiniteVolumeTerm(fluxFcn,(*myBcData)(iSub), (mgLevel->getGeoState())(iSub),
    //                                             V(iSub), matrices(iSub));
    mgLevel->getAgglomeratedFaces()[iSub]->
      computeJacobianFiniteVolumeTerm(fluxFcn,V(iSub), 
				      getBoundaryState()(iSub),
				      matrices(iSub), (*mgLSS)(iSub));
    /*
    if (fet && addViscousTerms) {
     
      mgLevel->getEdges()[iSub]->computeJacobianThinLayerViscousFiniteVolumeTerm(
        NULL, myVarFcn, fet, mgLevel->getGeoState()(iSub), 
        mgLevel->getXn()(iSub),  V(iSub), mgLevel->getCtrlVol()(iSub), matrices(iSub));

      mgLevel->getAgglomeratedFaces()[iSub]->computeJacobianThinLayerViscousFiniteVolumeTerm(
         fet,myVarFcn, V(iSub), (*DX[0])(iSub), 
         (*DX[1])(iSub), (*DX[2])(iSub), mgLevel->getCtrlVol()(iSub),
         mgLevel->getGeoState()(iSub).getDistanceToWall(),
         (*myBcData)(iSub).getFaceStateVector(),matrices(iSub));
          
	 }*/
 
    Vec<double>& ctrlVol = mgLevel->getCtrlVol()(iSub); 
    for (int i=0; i<ctrlVol.size(); ++i) {
      Scalar voli = 1.0 / ctrlVol[i];
      Scalar2 *Aii = matrices(iSub).getElem_ii(i);
      for (int k=0; k<neq*neq; ++k)
        Aii[k] *= voli;
      
      
    }

  }

  

  mgLevel->assemble(matrices);
  
  if (timeState) {

    timeState->addToJacobian(mgLevel->getCtrlVol(), matrices, U);
  }
  
  //applyBCsToJacobian(U, matrices); 

}

template<class Scalar,int dim>
template <class Scalar2>
void MultiGridOperator<Scalar,dim>::computeResidual(DistSVec<Scalar2,dim>& V,
                                                DistSVec<Scalar2,dim>& U,
  //                                           DistVec<Scalar2>& irey,
                                             FluxFcn** fluxFcn,
                                             RecFcn* recFcn,
                                             FemEquationTerm* fet,
                                             DistSVec<Scalar2,dim>& res,
                                             bool addDWdt) {


  ElemSet dummy;
  res = 0.0;
  DistVec<Scalar2>& irey = *scalar_zero; 
  
  NavierStokesTerm* nsterm = NULL;
  FemEquationTermSA* saterm = dynamic_cast<FemEquationTermSA*>(fet);
  if (fet) {

    FemEquationTermNS* nst = dynamic_cast<FemEquationTermNS*>(fet);
    if (nst)
      nsterm = static_cast<NavierStokesTerm*>(nst);
    else if (saterm) {
      nsterm = static_cast<NavierStokesTerm*>(saterm);

    }
 /*   else {

      fprintf(stderr, "Cannot create a NavierStokesTerm from FemEquationTerm!");
      exit(-1);
    }
*/
  }

  if (nsterm) {
    //std::cout << "Computing GG gradient" << std::endl;
    mgLevel->computeGreenGaussGradient(V, *DX[0],*DX[1],*DX[2]);
  }

//float dxn[3] = {(*DX[0])*(*DX[0]), (*DX[1])*(*DX[1]), (*DX[2])*(*DX[2])};
//std::cout << dxn[0] << " " << dxn[1] << " " << dxn[2] << std::endl;

#pragma omp parallel for
  for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {

    NodalGrad<dim,Scalar2> ngrad((*zero)(iSub), (*zero)(iSub), (*zero)(iSub),
	                         (*scalar_zero)(iSub), (*scalar_zero)(iSub), (*scalar_zero)(iSub));
    mgLevel->getEdges()[iSub]->template computeFiniteVolumeTerm<dim>(NULL, irey(iSub), fluxFcn,
                                         recFcn, dummy, (mgLevel->getGeoState())(iSub),
                                         (mgLevel->getGeoState().getXn())(iSub), V(iSub), ngrad, NULL, res(iSub),
                                         (mgLevel->getFVCompTag())(iSub), 0,0); 


    //mgLevel->getFaces()[iSub]->computeFiniteVolumeTerm(fluxFcn, (*myBcData)(iSub), 
    //                                     (mgLevel->getGeoState())(iSub), V(iSub),
    //                                     res(iSub));
    mgLevel->getAgglomeratedFaces()[iSub]->computeFiniteVolumeTerm(fluxFcn, 
                                         V(iSub),getBoundaryState()(iSub),
                                         res(iSub));

    if (fet && addViscousTerms) {
 
      /* mgLevel->getEdges()[iSub]->template 
        computeThinLayerViscousFiniteVolumeTerm<dim>(NULL, myVarFcn,
                                                     fet, (mgLevel->getGeoState())(iSub),
                                                     (mgLevel->getGeoState().getXn())(iSub),
                                                     V(iSub),res(iSub));
      
      */
      mgLevel->getEdges()[iSub]->template 
        computeViscousFiniteVolumeTerm<dim>(NULL, myVarFcn,
					    fet, (mgLevel->getGeoState())(iSub),
					    (mgLevel->getGeoState().getXn())(iSub),
					    V(iSub),
					    (*DX[0])(iSub), 
					    (*DX[1])(iSub), (*DX[2])(iSub) ,
					    res(iSub));
/*
      mgLevel->getAgglomeratedFaces()[iSub]->computeThinLayerViscousFiniteVolumeTerm(
                                      fet, myVarFcn,V(iSub), (*DX[0])(iSub), 
                                      (*DX[1])(iSub), (*DX[2])(iSub) ,
                                      mgLevel->getGeoState()(iSub).getDistanceToWall(),
                                      (*myBcData)(iSub).getFaceStateVector(),
                                      res(iSub)); 
  */    
    }    

  }
  
  mgLevel->assemble(res);

#pragma omp parallel for
  for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {

    for (int i = 0; i < V.subSize(iSub); ++i) {

      for (int j = 0; j < dim; ++j)
        res(iSub)[i][j] /= (mgLevel->getCtrlVol())(iSub)[i];
    }
  }

  // Add the S-A source term.
  /*  if (saterm) {

    double S[6];
#pragma omp parallel for
    for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {
      
      for (int i = 0; i < V.subSize(iSub); ++i) {

	double d2wall = mgLevel->getGeoState()(iSub).getDistanceToWall()[i];
	
	double dudxj[3][3] = {(*DX[0])(iSub)[i][1], (*DX[1])(iSub)[i][1], (*DX[2])(iSub)[i][1],
				      (*DX[0])(iSub)[i][2], (*DX[1])(iSub)[i][2], (*DX[2])(iSub)[i][2],
			      (*DX[0])(iSub)[i][3], (*DX[1])(iSub)[i][3], (*DX[2])(iSub)[i][3]};
	
	double dnudx[3] = {(*DX[0])(iSub)[i][5], (*DX[1])(iSub)[i][5], (*DX[2])(iSub)[i][5]};

	saterm->computeSourceTerm(dudxj, dnudx, d2wall, V(iSub)[i], S);
	//if (i == 0)
	//  std::cout << S[5] << std::endl;
	for (int k = 0; k < dim; ++k) {

	  res(iSub)[i][k] -= S[k];
	}

      }
    }
  }
  */

  /*
  if (timeState && addDWdt) {

    timeState->add_dAW_dt(0, mgLevel->getGeoState(),
                          mgLevel->getCtrlVol(), U, res);
  } 
  */

  applyBCsToResidual(U, res);
  
} 

template<class Scalar,int dim>
template <class Scalar2>
void MultiGridOperator<Scalar,dim>::
computeResidualEmbedded(DistExactRiemannSolver<dim>& riemann,
			DistSVec<Scalar2,dim>& V,
			DistSVec<Scalar2,dim>& U,
			//                                           DistVec<Scalar2>& irey,
			FluxFcn** fluxFcn,
			RecFcn* recFcn,
			FemEquationTerm* fet,
			DistSVec<Scalar2,dim>& res,	
			DistMultiGridLevelSetStructure* mgLSS,
			bool addDWdt) {
  

  ElemSet dummy;
  res = 0.0;
  DistVec<Scalar2>& irey = *scalar_zero; 
  
  NavierStokesTerm* nsterm = NULL;
  FemEquationTermSA* saterm = dynamic_cast<FemEquationTermSA*>(fet);
  if (fet) {

    FemEquationTermNS* nst = dynamic_cast<FemEquationTermNS*>(fet);
    if (nst)
      nsterm = static_cast<NavierStokesTerm*>(nst);
    else if (saterm) {
      nsterm = static_cast<NavierStokesTerm*>(saterm);

    } 
  }

  if (nsterm)
    mgLevel->computeGreenGaussGradient(V, *DX[0],*DX[1],*DX[2],mgLSS );

#pragma omp parallel for
  for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {

    int* locToGlobNodeMap = mgLevel->getMgSubDomains()[iSub].locToGlobMap;

    NodalGrad<dim,Scalar2> ngrad((*zero)(iSub), (*zero)(iSub), (*zero)(iSub),
	                         (*scalar_zero)(iSub), (*scalar_zero)(iSub), (*scalar_zero)(iSub));

    mgLevel->getEdges()[iSub]->
      template computeFiniteVolumeTerm<dim>(riemann(iSub),
 					    locToGlobNodeMap , fluxFcn,
					    recFcn, dummy, (mgLevel->getGeoState())(iSub),
					    (mgLevel->getGeoState().getXn())(iSub), V(iSub), 
					    (*Wstarij)(iSub), (*Wstarji)(iSub),
					    (*mgLSS)(iSub), false,
					    (*mgLSS).getStatus()(iSub), 1,
					    ngrad, NULL, res(iSub),0,
					    (mgLevel->getFVCompTag())(iSub), 0,0); 
    

    //mgLevel->getFaces()[iSub]->computeFiniteVolumeTerm(fluxFcn, (*myBcData)(iSub), 
    //                                     (mgLevel->getGeoState())(iSub), V(iSub),
    //                                     res(iSub));
        mgLevel->getAgglomeratedFaces()[iSub]->computeFiniteVolumeTerm(fluxFcn, 
								   V(iSub),getBoundaryState()(iSub),
								   res(iSub), (*mgLSS)(iSub));
    


    if (fet && addViscousTerms) {
 
//      mgLevel->getEdges()[iSub]->template 
//        computeThinLayerViscousFiniteVolumeTerm<dim>(NULL, myVarFcn,
//                                                     fet, (mgLevel->getGeoState())(iSub),
//                                                     (mgLevel->getGeoState().getXn())(iSub),
//                                                     V(iSub),res(iSub));
//
      mgLevel->getEdges()[iSub]->template 
        computeViscousFiniteVolumeTerm<dim>(NULL, myVarFcn,
					    fet, (mgLevel->getGeoState())(iSub),
					    (mgLevel->getGeoState().getXn())(iSub),
					    V(iSub),
					    (*DX[0])(iSub), 
					    (*DX[1])(iSub), (*DX[2])(iSub) ,
					    res(iSub), &(*mgLSS)(iSub));

//      mgLevel->getAgglomeratedFaces()[iSub]->computeThinLayerViscousFiniteVolumeTerm(
//                                      fet, myVarFcn,V(iSub), (*DX[0])(iSub), 
//                                      (*DX[1])(iSub), (*DX[2])(iSub) ,
//                                      mgLevel->getGeoState()(iSub).getDistanceToWall(),
//                                      (*myBcData)(iSub).getFaceStateVector(),
//                                      res(iSub)); 
                                                    
    }    
  
  }
  
  mgLevel->assemble(res);

#pragma omp parallel for
  for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {

    for (int i = 0; i < V.subSize(iSub); ++i) {

      for (int j = 0; j < dim; ++j)
        res(iSub)[i][j] /= (mgLevel->getCtrlVol())(iSub)[i];
    }
  }
  
  if (timeState && addDWdt) {

    timeState->add_dAW_dt(0, mgLevel->getGeoState(),
                          mgLevel->getCtrlVol(), U, res);
  } 
  
  //applyBCsToResidual(U, res);
  
} 

template<class Scalar,int dim>
template <class Scalar2>
void MultiGridOperator<Scalar,dim>::
add_dAW_dtEmbedded(DistSVec<Scalar2,dim>& U,
                   DistSVec<Scalar2,dim>& res,	
                   DistMultiGridLevelSetStructure* mgLSS) {

  timeState->add_dAW_dt(0, mgLevel->getGeoState(),
			mgLevel->getCtrlVol(), U, res);  
}

template<class Scalar,int dim>
void MultiGridOperator<Scalar,dim>::computeGradientsLeastSquares(DistSVec<Scalar,dim>& V) {

  EdgeSet** edges = mgLevel->getEdges();  

#pragma omp parallel for
  for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {
    
    bool *edgeFlag = edges[iSub]->getMasterFlag();
    int (*edgePtr)[2] = edges[iSub]->getPtr();

    SVec<Scalar,dim>& ddx = myNodalGrad->getX()(iSub);
    SVec<Scalar,dim>& ddy = myNodalGrad->getY()(iSub);
    SVec<Scalar,dim>& ddz = myNodalGrad->getZ()(iSub);
    SVec<double,6>& R = myNodalGrad->getR()(iSub); 
    SVec<double,3>& X = (mgLevel->getGeoState().getXn())(iSub);
    SVec<Scalar,dim>& var = V(iSub);
    for (int l=0; l<edges[iSub]->size(); ++l) {

      if (!edgeFlag[l])
        continue;

      int i = edgePtr[l][0];
      int j = edgePtr[l][1];

      double Wi[3], Wj[3];
      Scalar deltaVar;

      double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
      computeLocalWeightsLeastSquares(dx, R[i], Wi);

      dx[0] = -dx[0]; dx[1] = -dx[1]; dx[2] = -dx[2];
      computeLocalWeightsLeastSquares(dx, R[j], Wj);

      for (int k=0; k<dim; ++k) {
        deltaVar = var[j][k] - var[i][k];

        ddx[i][k] += Wi[0] * deltaVar;
	ddy[i][k] += Wi[1] * deltaVar;
	ddz[i][k] += Wi[2] * deltaVar;
	ddx[j][k] -= Wj[0] * deltaVar;
	ddy[j][k] -= Wj[1] * deltaVar;
	ddz[j][k] -= Wj[2] * deltaVar;
      }
    }
  }

  mgLevel->assemble(myNodalGrad->getX());
  mgLevel->assemble(myNodalGrad->getY());
  mgLevel->assemble(myNodalGrad->getZ());

}

template<class Scalar,int dim>
void MultiGridOperator<Scalar,dim>::updateStateVectors(DistSVec<Scalar,dim>& U) {

  timeState->update(U);
}

template<class Scalar,int dim>
void MultiGridOperator<Scalar,dim>::computeTimeStep(double cfl, VarFcn *varFcn,
                                                    DistSVec<double,dim> &V) {

  int iSub;
  DistVec<double>& dt = timeState->getDt();
  dt = 0.0;
  *idti = 0.0;
  *idtv = 0.0;

#pragma omp parallel for
  for (iSub = 0; iSub < V.numLocSub(); ++iSub) {
    mgLevel->getEdges()[iSub]->computeTimeStep2(NULL, varFcn,  (mgLevel->getGeoState())(iSub),
                                 (mgLevel->getGeoState().getXn())(iSub),
                                 V(iSub), (*idti)(iSub), (*idtv)(iSub), 
                                 timeState->getTimeLowMachPrec(),
                                 mgLevel->getEdgeArea()(iSub));
    mgLevel->getAgglomeratedFaces()[iSub]->computeTimeStep(NULL, varFcn,  
                                 (mgLevel->getGeoState().getXn())(iSub),
                                 V(iSub), (*idti)(iSub), (*idtv)(iSub), 
                                 timeState->getTimeLowMachPrec());
  }

  mgLevel->assemble(*idti);
  *idtv = 0.0;

#pragma omp parallel for
  for (iSub = 0; iSub < V.numLocSub(); ++iSub) {
    double (*idtimev) = idtv->subData(iSub);
    double (*idtimei) = idti->subData(iSub);
    double (*dtime) = dt.subData(iSub);
    //double (*ireynolds) = irey.subData(iSub);
    double (*volume) = mgLevel->getCtrlVol().subData(iSub);
    for (int i = 0; i < mgLevel->getCtrlVol().subSize(iSub); ++i) {
      //   idtimev[i] = idtimev[i] / volume[i];
      dtime[i] = cfl *volume[i]/(1.0*idtimei[i]/* + viscous*idtimev[i]*/);
      //ireynolds[i] = -sprec.getViscousRatio()*idtimev[i] / idtimei[i];
    }
  }

  double dt_glob = dt.min();
  timeState->computeCoefficients(dt_glob);
}

template<class Scalar,int dim>
template <class Scalar2>
void MultiGridOperator<Scalar,dim>::applyBCsToResidual(DistSVec<Scalar2,dim>& U,
                                                       DistSVec<Scalar2,dim>& R) {

  if (!bcFcn)
    return;
 
#pragma omp parallel for
  for (int iSub = 0; iSub < U.numLocSub(); ++iSub) {
    SVec<double,dim> &Vwall = (*myBcData)(iSub).getNodeStateVector();

    int* nodeType = mgLevel->getNodeType(iSub);
    for (int i = 0; i < U(iSub).size(); ++i) {
 
      if (nodeType[i] != BC_INTERNAL)
        bcFcn->applyToResidualTerm(nodeType[i], Vwall[i], U(iSub)[i], R(iSub)[i]);
    }
  }
}

template<class Scalar,int dim>
template <class Scalar2,int neq>
void MultiGridOperator<Scalar,dim>::applyBCsToJacobian(DistSVec<Scalar2,dim>& U,
                                                       DistMvpMatrix<Scalar2,neq>& A) {
 
  BcFcn* B = NULL;
  if (neq > 2) {
    B = bcFcn1;
  } else {
    
    B = bcFcn2;
  }
  
  if (!B)
    return;

  EdgeSet** edges = mgLevel->getEdges();

#pragma omp parallel for
  for (int iSub = 0; iSub < U.numLocSub(); ++iSub) {

    SVec<double,dim> &Vwall = (*myBcData)(iSub).getNodeStateVector();

    int (*edgePtr)[2] = edges[iSub]->getPtr();
    int* nodeType = mgLevel->getNodeType(iSub);

    for (int l=0; l<edges[iSub]->size(); ++l) {
      int i = edgePtr[l][0];
      int j = edgePtr[l][1];

      if (nodeType[i] != BC_INTERNAL)  {
        Scalar *Aij = A(iSub).getElem_ij(l);
        if (Aij)
          B->applyToOffDiagonalTerm(nodeType[i], Aij);
      }

      if (nodeType[j] != BC_INTERNAL) {
        Scalar *Aji = A(iSub).getElem_ji(l);
        if (Aji)
          B->applyToOffDiagonalTerm(nodeType[j], Aji);
      }
    }

    for (int i=0; i<U(iSub).size(); ++i) {
      if (nodeType[i] != BC_INTERNAL) {
        Scalar *Aii = A(iSub).getElem_ii(i);
        if (Aii)
          B->applyToDiagonalTerm(nodeType[i], Vwall[i], U(iSub)[i], Aii);
      }
    }

  }

}

template<class Scalar,int dim>
double MultiGridOperator<Scalar,dim>::queryTimeStep(int iSub, int i) {

  return timeState->getDt()(iSub)[i];
}

#define INST_HELPER(S,D) \
template \
void MultiGridOperator<S,D>::computeJacobian(DistSVec<double,D>& U, DistSVec<double,D>& V,\
                                             FluxFcn **fluxFcn, \
                                             FemEquationTerm*,DistMvpMatrix<double,D>& matrices); \
 template void MultiGridOperator<S,D>::computeJacobianEmbedded(DistExactRiemannSolver<D>&,DistSVec<double,D>& U, DistSVec<double,D>& V, \
FluxFcn **fluxFcn, \
							       FemEquationTerm*,DistMvpMatrix<double,D>& matrices, DistMultiGridLevelSetStructure*); \
template void MultiGridOperator<S,D>::applyBCsToJacobian(DistSVec<double,D>& U, \
                                                       DistMvpMatrix<double,D>& A); \
template void MultiGridOperator<S,D>::computeResidual(DistSVec<double,D>& V, \
                                                DistSVec<double,D>& U, \
                                             FluxFcn** fluxFcn,\
                                             RecFcn* recFcn,\
                                             FemEquationTerm*, \
                                             DistSVec<double,D>& res,bool); \
template void MultiGridOperator<S,D>::computeResidualEmbedded(DistExactRiemannSolver<D>&,DistSVec<double,D>& V, \
                                                DistSVec<double,D>& U, \
                                             FluxFcn** fluxFcn,\
                                             RecFcn* recFcn,\
                                             FemEquationTerm*, \
                                             DistSVec<double,D>& res, DistMultiGridLevelSetStructure*,bool); \
template void MultiGridOperator<S,D>::applyBCsToResidual(DistSVec<double,D>& U, \
                                                       DistSVec<double,D>& R); \
template \
void MultiGridOperator<S,D>:: \
add_dAW_dtEmbedded(DistSVec<S,D>& U, \
                   DistSVec<S,D>& res, \
                   DistMultiGridLevelSetStructure* mgLSS);

#define INST_HELPER2(S,D,E) \
template void MultiGridOperator<S,D>::computeJacobian(DistSVec<double,D>& U, DistSVec<double,D>& V,\
FluxFcn **fluxFcn, \
FemEquationTerm*,DistMvpMatrix<double,E>& matrices); \
 template void MultiGridOperator<S,D>::computeJacobianEmbedded(DistExactRiemannSolver<D>&,DistSVec<double,D>& U, DistSVec<double,D>& V, \
FluxFcn **fluxFcn, \
							       FemEquationTerm*,DistMvpMatrix<double,E>& matrices, DistMultiGridLevelSetStructure*); \
template void MultiGridOperator<S,D>::applyBCsToJacobian(DistSVec<double,D>& U, \
							 DistMvpMatrix<double,E>& A); 

template class MultiGridOperator<double,1>;
template class MultiGridOperator<double,2>;
template class MultiGridOperator<double,5>;
template class MultiGridOperator<double,6>;
template class MultiGridOperator<double,7>;

INST_HELPER(double,1);
INST_HELPER(double,2);
INST_HELPER(double,5);
INST_HELPER(double,6);
INST_HELPER(double,7);

INST_HELPER2(double,6,1);
INST_HELPER2(double,7,2);

INST_HELPER2(double,6,5);
INST_HELPER2(double,7,5);
