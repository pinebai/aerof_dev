//------------------------------------------------------------------------------

template<int dim>
ImplicitPGTsDesc<dim>::ImplicitPGTsDesc(IoData &ioData, GeoSource &geoSource, Domain *dom) :
  residualRef(this->F),
  ImplicitRomTsDesc<dim>(ioData, geoSource, dom), From(this->nPod), rhs(this->nPod) {

  pc = ImplicitRomTsDesc<dim>::template 
  createPreconditioner<PrecScalar,dim>(this->ioData->ts.implicit.newton.ksp.ns.pc, this->domain);

  // TODO necessary?
  currentProblemSize = this->nPod;
  lsCoeff = new double*[1];
  lsCoeff[0] = new double[this->nPod];
  this->projVectorTmp = new double [this->nPod];

  parallelRom = NULL;
  jactmp = NULL;

  lsSolver = this->ioData->romOnline.lsSolver;
  if ((lsSolver==NonlinearRomOnlineData::QR) || (lsSolver==NonlinearRomOnlineData::LEVENBERG_MARQUARDT_SVD)) {
    parallelRom = new ParallelRom<dim>(*dom,this->com,dom->getNodeDistInfo());
  } else if (lsSolver==NonlinearRomOnlineData::NORMAL_EQUATIONS) { 
    jactmp = new double [this->nPod * this->nPod];
    this->jac.setNewSize(this->nPod,this->nPod);
  }

  this->rom->initializeClusteredOutputs();

}

//------------------------------------------------------------------------------
template<int dim>
ImplicitPGTsDesc<dim>::~ImplicitPGTsDesc(){

    delete [] lsCoeff[0];
    delete [] lsCoeff;  
    if (this->projVectorTmp) delete [] this->projVectorTmp;
    if (jactmp) delete [] jactmp;
    if (pc) delete pc;
    if (parallelRom) delete parallelRom;
}


//-----------------------------------------------------------------------------
template<int dim>
void ImplicitPGTsDesc<dim>::solveNewtonSystem(const int &it, double &res, bool &breakloop, DistSVec<double, dim> &U, const int& totalTimeSteps)  {

  this->projectVector(this->AJ, this->F, From);
  Vec<double> rhs(this->nPod);
  rhs = -1.0 * From;

  res = rhs*rhs;
  double resReg = 0;

  // print some residual information
  this->com->fprintf(stdout, " ... component-wise residual norms are");
  double (*iDimMask)[dim] = new double[this->domain->getNodeDistInfo().totLen][dim];
  for (int iDim=0; iDim<dim; ++iDim) {
    for (int iNode=0; iNode<this->domain->getNodeDistInfo().totLen; ++iNode) {
      for (int jDim=0; jDim<dim; ++jDim) {
        iDimMask[iNode][jDim] = 0.0;
      }
      iDimMask[iNode][iDim] = 1.0;
    }
    DistSVec<double, dim> iDimMaskedRes(this->domain->getNodeDistInfo(), iDimMask);
    iDimMaskedRes *= this->F;
    this->com->fprintf(stdout, " %e", iDimMaskedRes.norm());
  }
  this->com->fprintf(stdout, "\n");
  delete [] iDimMask;

  if (res < 0.0){
    this->com->fprintf(stderr, "*** negative residual: %e\n", res);
    exit(1);
  }
  res = sqrt(res);

  if (it == 0) {
    this->target = this->epsNewton*res;
    this->res0 = res;
  }

  if (res == 0.0 || res <= this->target || res <= this->epsAbsResNewton) {
    breakloop = true;
    return;  // do not solve the system
  }

  double t0 = this->timer->getTime();
 
  if (lsSolver==NonlinearRomOnlineData::QR) {  // ScaLAPACK least-squares

    RefVec<DistSVec<double, dim> > residualRef2(this->F);
    parallelRom->parallelLSMultiRHS(this->AJ,residualRef2,this->nPod,1,lsCoeff);
    double dUromNewtonItNormSquared = 0;
    for (int iPod=0; iPod<this->nPod; ++iPod) {
      this->dUromNewtonIt[iPod] = -lsCoeff[0][iPod];
      dUromNewtonItNormSquared += pow(this->dUromNewtonIt[iPod],2);
    }

  } else if (lsSolver==NonlinearRomOnlineData::PROBABILISTIC_SVD) {

    VecSet<DistSVec<double, dim> > resTmp(1, this->domain->getNodeDistInfo());
    resTmp[0] = this->F;

    VecSet<DistSVec<double, dim> > ajTmp(this->AJ.numVectors(), this->domain->getNodeDistInfo());
    for (int iVec=0; iVec<this->AJ.numVectors(); ++iVec) ajTmp[iVec] = this->AJ[iVec];

    std::vector<std::vector<double> > lsCoeffVec;
    lsCoeffVec.resize(1);
    lsCoeffVec[0].resize(this->nPod,0.0);

    this->rom->probabilisticLSMultiRHS(ajTmp, resTmp, lsCoeffVec, this->ioData->romOnline.randMatDimension * (it+1), 0, false);

    double dUromNewtonItNormSquared = 0;
    for (int iPod=0; iPod<this->nPod; ++iPod) {
      lsCoeff[0][iPod] = lsCoeffVec[0][iPod];
      this->dUromNewtonIt[iPod] = -lsCoeff[0][iPod];
      dUromNewtonItNormSquared += pow(this->dUromNewtonIt[iPod],2);
    }
    this->com->fprintf(stdout, " ... || dUromNewtonIt ||^2 = %1.12e \n", dUromNewtonItNormSquared);

  } else if (lsSolver==NonlinearRomOnlineData::LEVENBERG_MARQUARDT_SVD){  // Solve Levenberg-Marquardt regularized LS via ScaLAPACK SVD

    // SVD quantities
    VecSet< DistSVec<double, dim> >* U_AJ = new VecSet< DistSVec<double, dim> >(this->nPod, this->domain->getNodeDistInfo());
    Vec<double>* sVals_AJ = new Vec<double>(this->nPod);
    FullM* V_AJ = new FullM(this->nPod);

    parallelRom->parallelSVD(this->AJ, *U_AJ, sVals_AJ->data(), *V_AJ, this->nPod, true);
    // Note: V, not V_transpose

   /* double maxErr = 0.0;
    double avgErr = 0.0;
    for (int iPod = 0; iPod < this->nPod; ++iPod) {
      DistSVec<double, dim> error(this->domain->getNodeDistInfo());
      error  = this->AJ[iPod];
      for (int jPod = 0; jPod < this->nPod; ++jPod)
        error = error - (((*singVals)[jPod]*(*Vtrue)[iPod][jPod])*(*Utrue)[jPod]);
      double errorNorm = error.norm()/((this->AJ[iPod]).norm());
      avgErr += errorNorm;
      if (errorNorm > maxErr)
        maxErr = errorNorm;   
    }
    avgErr /= this->nPod;
  
    this->com->fprintf(stderr, " ... Average error on AJ after SVD = %e\n", avgErr);  
    this->com->fprintf(stderr, " ... Maximum error on AJ after SVD = %e\n", maxErr); */
    
    Vec<double> tmpVec(this->nPod);
    for (int iVec=0; iVec<U_AJ->numVectors(); ++iVec)
      tmpVec[iVec] = (*U_AJ)[iVec] * (-1.0*this->F);

    delete U_AJ;

    if ((*sVals_AJ)[this->nPod-1]>0) {
      this->com->fprintf(stdout, " ... Singular value ratio for AJ = %e \n", (*sVals_AJ)[0]/(*sVals_AJ)[this->nPod-1]);
    } else {
      this->com->fprintf(stdout, " ... AJ is rank deficient! \n");
    }

    double lambdaSquared = pow(this->levenbergMarquardtWeight,2);

    double firstTermNormSquared = 0.0;
    for (int iVec=0; iVec<this->nPod; ++iVec) {
      if ((*sVals_AJ)[iVec]>0) {
        firstTermNormSquared += pow(tmpVec[iVec],2)*pow(lambdaSquared/(pow((*sVals_AJ)[iVec],2) + lambdaSquared),2);
      } else {
        firstTermNormSquared += pow(tmpVec[iVec],2);
      }
    }

    double secondTermNormSquared = 0.0;
    for (int iVec=0; iVec<this->nPod; ++iVec) {
      tmpVec[iVec] = ((*sVals_AJ)[iVec]>0) ? tmpVec[iVec]*(*sVals_AJ)[iVec]/(pow((*sVals_AJ)[iVec],2) + lambdaSquared) : 0.0;
      secondTermNormSquared += pow(tmpVec[iVec],2);
    }

    delete sVals_AJ;

    double dUromNewtonItNormSquared = 0;
    for (int iPod=0; iPod<this->nPod; ++iPod) {
      this->dUromNewtonIt[iPod] = 0.0;
      for (int jPod=0; jPod<this->nPod; ++jPod) {
        this->dUromNewtonIt[iPod] += (*V_AJ)[iPod][jPod] * tmpVec[jPod];
      }
      dUromNewtonItNormSquared += pow(this->dUromNewtonIt[iPod],2);
      //this->com->fprintf(stdout, " ... dUromNewtonIt[%d] = %e \n", iPod, this->dUromNewtonIt[iPod]);
    }
    this->com->fprintf(stdout, " ... || dUromNewtonIt ||^2 = %e \n", dUromNewtonItNormSquared);
    this->com->fprintf(stdout, " ... || Ax - b ||^2 = %1.12e \n", firstTermNormSquared);
    this->com->fprintf(stdout, " ... || x ||^2 = %1.12e \n", secondTermNormSquared);
 
    delete V_AJ;

  } else if (lsSolver==NonlinearRomOnlineData::NORMAL_EQUATIONS)  {    // normal equations

    transMatMatSymProd(this->AJ,jactmp);  
    for (int iRow = 0; iRow < this->nPod; ++iRow) {
      for (int iCol = 0; iCol < this->nPod; ++iCol) {
        this->jac[iRow][iCol] = jactmp[iRow + iCol * this->nPod];
      }
    } 

    // homotopy on reduced-coordinates for spatial-only problems
    if (this->spatialOnlyWithHomotopy) {
      double homotopyStep = min(this->homotopyStepInitial*pow(this->homotopyStepGrowthRate,totalTimeSteps),this->homotopyStepMax);
      this->com->fprintf(stdout, " ... homotopy step = %1.12e \n", homotopyStep);
      double invHomotopyStep = 1/homotopyStep;
      Vec<double> dUrom(this->dUromTimeIt);
      dUrom *= invHomotopyStep;
      rhs -= dUrom;
      for (int iDiag = 0; iDiag < this->nPod; ++iDiag) {
        this->jac[iDiag][iDiag] += invHomotopyStep;
      }
    }

    this->solveLinearSystem(it, rhs, this->dUromNewtonIt);

    double dUromNewtonItNormSquared = 0;
    for (int iPod=0; iPod<this->nPod; ++iPod) {
      //this->com->fprintf(stdout, " ... dUromNewtonIt[%d] = %e \n", iPod, this->dUromNewtonIt[iPod]);
      dUromNewtonItNormSquared += pow(this->dUromNewtonIt[iPod],2);
    }
    //this->com->fprintf(stdout, " ... || dUromNewtonIt ||^2 = %e \n", dUromNewtonItNormSquared);

  } 
  
  this->timer->addLinearSystemSolveTime(t0); 

}

//-----------------------------------------------------------------------------

template<int dim>
void ImplicitPGTsDesc<dim>::setProblemSize(DistSVec<double, dim> &U) {

  if (currentProblemSize != this->nPod){
    currentProblemSize = this->nPod;

    if (lsSolver==NonlinearRomOnlineData::QR) parallelRom->parallelLSMultiRHSInit(this->AJ,residualRef);

    if (lsCoeff) delete [] lsCoeff[0];
    lsCoeff[0] = new double[this->nPod];

    if (this->projVectorTmp) delete [] (this->projVectorTmp);
    this->projVectorTmp = new double [this->nPod];

    if (lsSolver==NonlinearRomOnlineData::NORMAL_EQUATIONS) {
      if (jactmp) delete [] jactmp;
      jactmp = new double [this->nPod * this->nPod];
      this->jac.setNewSize(this->nPod,this->nPod);
    }

    From.resize(this->nPod);
    rhs.resize(this->nPod);
  }


}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitPGTsDesc<dim>::monitorInitialState(int it, DistSVec<double,dim> &U)
{

  this->com->printf(2, "State vector norm = %.12e\n", sqrt(U*U));

  if (!this->problemType[ProblemData::UNSTEADY]) {
    this->com->printf(2, "\nNOTE: For weighted ROM simulations the reported residual is calculated using a weighted norm,\n");
    this->com->printf(2, "      and is relative to the residual of the initial condition calculated using the same norm.\n");

    this->Uinit = new DistSVec<double, dim>(this->domain->getNodeDistInfo());
    *(this->Uinit) = U;  // needed for computing the restricted residual after each cluster switch
  }

  this->com->printf(2, "\n");
  
}   

//------------------------------------------------------------------------------

template<int dim>
bool ImplicitPGTsDesc<dim>::checkForLastIteration(IoData &ioData, int it, double t, double dt, DistSVec<double,dim> &U)
{

  if (!this->problemType[ProblemData::UNSTEADY] && monitorConvergence(it, U))
    return true;

  if (!this->problemType[ProblemData::AERO] && !this->problemType[ProblemData::THERMO] && it >= this->data->maxIts) return true;

  if (this->problemType[ProblemData::UNSTEADY] )
    if(t >= this->data->maxTime - 0.01 * dt)
      return true;

  return false;

}


//------------------------------------------------------------------------------
 
template<int dim>
bool ImplicitPGTsDesc<dim>::monitorConvergence(int it, DistSVec<double,dim> &U)
{// only called for steady simulations

  this->data->residual = computePGResidualNorm(U);

  if (this->data->residual == 0.0 || this->data->residual < this->data->eps * this->restart->residual || this->data->residual < this->data->epsabs)
    return true;
  else
    return false;

}

//------------------------------------------------------------------------------

template<int dim>
double ImplicitPGTsDesc<dim>::computePGResidualNorm(DistSVec<double,dim>& Q)
{ // spatial only

  this->spaceOp->computeResidual(*this->X, *this->A, Q, this->F, this->timeState);

  this->spaceOp->applyBCsToResidual(Q, this->F);

  double res = 0.0;
  res = (this->F) * (this->F);
  return sqrt(res);

}

//------------------------------------------------------------------------------

template<int dim>
void ImplicitPGTsDesc<dim>::setReferenceResidual()
{
  if (this->Uinit) this->restart->residual = computePGResidualNorm(*(this->Uinit));

  this->com->printf(2, "Norm of reference residual = %.12e\n", this->restart->residual);

}

//------------------------------------------------------------------------------
template<int dim>
double ImplicitPGTsDesc<dim>::meritFunction(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> &dQ, DistSVec<double, dim> &R, double stepLength)  {
	// merit function: norm of the residual (want to minimize residual)

  DistSVec<double, dim> newQ(this->domain->getNodeDistInfo());
  newQ = Q + stepLength*dQ;
  this->computeFullResidual(it,newQ,true,&R);

  double merit = 0.0;

  if (this->ioData->romOnline.meritFunction == NonlinearRomOnlineData::HDM_RESIDUAL) {
    // merit function = 1/2 * (norm of full-order residual)^2
    merit = R.norm();	
    merit *= merit;
    merit *= 0.5;
  } else if (this->ioData->romOnline.meritFunction == NonlinearRomOnlineData::ROM_RESIDUAL) {
    // Form the ROM residual
    Vec<double> romResidual;
    romResidual.resize(this->nPod);
    this->computeAJ(it, newQ, true);
    this->projectVector(this->AJ, R, romResidual);
    merit = romResidual*romResidual;
    merit *= 0.5;
  } else {
    fprintf(stderr,"*** Error: unrecognized choice of merit function!\n");
    exit(-1);
  }

  return merit;

}
