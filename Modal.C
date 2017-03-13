#include <cstdio>
#include <cmath>
#include <sys/time.h>
#include <algorithm>
#include <cstdlib>
using std::sort;

#include <Domain.h>
#include <Modal.h>
#include <IoData.h>
#include <DistVector.h>
#include <VectorSet.h>
#include <MatVecProd.h>
#include <Timer.h>
//#include <VarFcnDesc.h>
#include <DistBcData.h>
#include <DistGeoState.h>
#include <DistTimeState.h>
#include <PostOperator.h>
#include <ParallelRom.h>
#ifdef USE_EIGEN3
#include <Eigen/Eigenvalues>
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
#endif
#include <complex>

#ifdef DO_MODAL
 #include <arpack++/include/ardsmat.h>
 #include <arpack++/include/ardnsmat.h>
 #include <arpack++/include/ardssym.h>
 #include <arpack++/include/ardsnsym.h>
#endif

#include <sstream>
#include <string>
#include <cstring>
extern "C"      {

	void F77NAME(dsvdc)(double *, int &, int &, int&, double *,
				               double *, double *, int &, double *, int &,
											 double *, const int &, int &);

  void F77NAME(thinsvd)(int &, int &, int &, int &, int&, int&, int&, int&, int&, double *, int &,
						           int &, int &, int &, int &, int &, int &, double *U, double *S, double *V,
											int &, double *, int &, int &);
	void F77NAME(lworksize)(int &, int &, int &, int &, int&, int&, int&, int&, int&, int &);

}

//----------------------------------------------------------------------------------

template <int dim>
ModalSolver<dim>::ModalSolver(Communicator *_com, IoData &_ioData, Domain &dom) : 
          domain(dom), mX(0, dom.getNodeDistInfo() ), Xref(dom.getNodeDistInfo()), 
          Uref(dom.getNodeDistInfo()), DX(0, dom.getNodeDistInfo()), 
          DE(0, dom.getNodeDistInfo()), controlVol(dom.getNodeDistInfo()), 
					controlVolComp(dom.getNodeDistInfo())  {

 com = _com;
 double f = 0;
 pi = 3.14159265358979;
 ioData = &_ioData; 
 const char *modeFile = ioData->linearizedData.strModesFile;

 DistSVec<double, dim> tmpVec(domain.getNodeDistInfo());
 DistSVec<double, 3> Xtmp(domain.getNodeDistInfo());

 if (strcmp(modeFile, "") != 0)  {
   com->fprintf(stderr, " ... Reading Modefile %s\n", modeFile);
   modeFile = ioData->linearizedData.strModesFile;
   domain.readVectorFromFile(modeFile, 0, &f, Xtmp);
   nStrMode = int(f);
   if (ioData->linearizedData.numStrModes > nStrMode)  {
     com->fprintf(stderr, " *** WARNING: Setting number of structural modes to number in file: %d\n",
                nStrMode);
   }
   else
     nStrMode = ioData->linearizedData.numStrModes;

 }
 else  {
   com->fprintf(stderr, " ... Running Fluid Alone without structural mode file\n");
   nStrMode = 0;
 }
 K = new double[nStrMode];

 // We read the modal deformations
 mX.resize(nStrMode);
 com->fprintf(stderr, " ... Number of modes in file: %d\n", nStrMode);

 for(int iMode = 0; iMode < nStrMode; ++iMode) {
   domain.readVectorFromFile(modeFile, iMode+1, &f, mX[iMode]);
   com->fprintf(stderr, "Mode %d f = %e disp^2 = %e\n", iMode+1, f, mX[iMode]*mX[iMode]);
   K[iMode] = f*f*pi*pi*4.0;
 }

 tInput = new TsInput(*ioData);

 tOutput = 0;
 tRestart = 0;
 bcData = 0;
 geoState = 0;
 tState = 0;
 spaceOp = 0;
 postOp = 0;
 HOp = 0;
 HOp2 = 0;
 HOp2step1 = 0;
 HOpstep2 = 0;
 HOpstep3 = 0;
 pc = 0;
 ksp = 0;
 ksp2 = 0;
 ksp3 = 0;
 kspComp = 0;
 kspCompGcr = 0;
 totalEnergy = 0.0;	// for POD basis construction

}

//-------------------------------------------------------------------------------
template <int dim> 
void ModalSolver<dim>::solve()  {

 // set up Timers
 Timer *modalTimer = domain.getTimer();
 double t0;

 modalTimer->setSetupTime();//CBM

 if (ioData->problem.alltype == ProblemData::_INTERPOLATION_)
   interpolatePOD();
 else if (ioData->problem.alltype == ProblemData::_ROB_INNER_PRODUCT_)
   ROBInnerProducts();
 else  {
   preProcess();
   if (ioData->problem.alltype == ProblemData::_POD_CONSTRUCTION_) {
     t0 = modalTimer->getTime();
     podMethod = 0;
     constructPOD();
     modalTimer->addPodConstrTime(t0);
   }
   else if (ioData->problem.alltype == ProblemData::_AEROELASTIC_ANALYSIS_) {
     // t0 = modalTimer->getTime();
     computeDampingRatios();
   }
#ifdef USE_EIGEN3
   else if (ioData->problem.alltype == ProblemData::_NONLINEAR_EIGEN_ERROR_INDICATOR_) {
     double sReal, sImag, normalizationTerm;
     Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic> rightEigenVector(2*nStrMode,1), leftEigenVector(2*nStrMode,1), residual(nStrMode,1);
     rightEigenVector.setZero();    leftEigenVector.setZero();    residual.setZero();
     char eigFile [100];  int iEV;
     sprintf(eigFile, "%s%s", ioData->input.prefix, ioData->input.reducedEigState);
     complex<double>* rEigenVector = rightEigenVector.data();
     complex<double>* lEigenVector = leftEigenVector.data();
     domain.readEigenValuesAndVectors(eigFile, sReal, sImag, normalizationTerm, iEV, rEigenVector, lEigenVector);
//     cerr << rightEigenVector;
//     cerr << leftEigenVector;
//     computeEigenvectorsAndResidual(sReal,sImag,iEV,rightEigenVector,leftEigenVector,residual);

     com->fprintf(stderr,"ioData->linearizedData.errorIndicator is %d\n", ioData->linearizedData.errorIndicator);
     switch (ioData->linearizedData.errorIndicator) {
       case LinearizedData::OIBEI :
         com->fprintf(stderr," ... computing OIBEI\n");
         computeOIBEI(sReal,sImag,iEV);
         break;
       case LinearizedData::RBEI1 :
         com->fprintf(stderr," ... computing RBEI1\n");
         computeREigenvector(sReal, sImag, iEV, rightEigenVector);
         computeNonlinearEigenResidual(sReal, sImag, rightEigenVector, residual, "RBEI1");
         break;
       case LinearizedData::RBEI2 :
         com->fprintf(stderr," ... computing RBEI2\n");
         computeNonlinearEigenResidual(sReal, sImag, rightEigenVector, residual, "RBEI2");
         break;
       case LinearizedData::RBEI3 :
         com->fprintf(stderr," ... computing RBEI3\n");
         computeREigenvector(sReal, sImag, iEV, rightEigenVector);
         computeNonlinearEigenResidualNormalized(sReal, sImag, rightEigenVector, normalizationTerm, residual, "RBEI3");
         break;
       case LinearizedData::RBEI4 :
         com->fprintf(stderr," ... computing RBEI4\n");
         computeNonlinearEigenResidualNormalized(sReal, sImag, rightEigenVector, normalizationTerm, residual, "RBEI4");
         break;
       default:
         cout << " ... Error! invalid error indicator type\n";
     }

/*
     if (isFirstErrorIndicator) {
       double residualnorm = residual.norm(); 
       double resDenominator = computeResidualDenominator(sReal, sImag, rightEigenVector, leftEigenVector);
       com->fprintf(stderr,"residualnorm is %e\n", residualnorm); 
       com->fprintf(stderr,"resDenominator is %e\n", resDenominator); 
       com->fprintf(stderr,"errorEstimate is %e\n", residualnorm/resDenominator); 
       const char* output = "errorIndicator1";
       ofstream out(output, ios::out);
       if(!out) {
         cerr << "Error: cannot open file" << output << endl;
         exit(-1);
       } 
       out << residualnorm/resDenominator << endl;
       out.close();
//       delete [] rEigenVector;
     }
*/
   }
#endif
   else if (ioData->problem.alltype == ProblemData::_GAM_CONSTRUCTION_) {
     // t0 = modalTimer->getTime();
     computeGenAeroForceMat();
   }
   else {
     t0 = modalTimer->getTime();
     solveInTimeDomain();
     //modalTimer->addRomSolTime(t0);
     modalTimer->addFluidSolutionTime(t0);
   }

 }

 modalTimer->setRunTime();
 modalTimer->print(modalTimer);

/*
 // set up Timers
 Timer *modalTimer = domain.getTimer();
 double t0;

 if (ioData->problem.alltype == ProblemData::_INTERPOLATION_)
   interpolatePOD();
 else  {
   preProcess();
   modalTimer->setSetupTime();
   if (ioData->problem.alltype == ProblemData::_ROB_CONSTRUCTION_) {
     t0 = modalTimer->getTime();
     constructPOD();
     modalTimer->addPodConstrTime(t0);
   }
   else {
     t0 = modalTimer->getTime();
     solveInTimeDomain();
     //modalTimer->addRomSolTime(t0);
     modalTimer->addFluidSolutionTime(t0);
   }
 }

 modalTimer->setRunTime();
 modalTimer->print(modalTimer);
*/
}

//-------------------------------------------------------------------------------
template <int dim> 
void ModalSolver<dim>::solveInTimeDomain()  {

 Timer *modalTimer = domain.getTimer();

 double sdt = ioData->linearizedData.stepsize;

 // Initial Conditions
 double *delY = new double [nStrMode];;
 double *delU = new double [nStrMode];;

 //Initial Conditions
 int i;
 for (i=0; i<nStrMode; ++i) {
   delY[i] = 0.0;
   delU[i] = 0.0;
 }

 int nSteps = ioData->ts.maxIts;

 if (ioData->problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_AEROELASTIC_ || 
     ioData->problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_ ||
     ioData->linearizedData.type == LinearizedData::FORCED)  {
   VecSet<DistSVec<double, dim> > snaps(0, domain.getNodeDistInfo());

   int modeNum = ioData->linearizedData.modeNumber - 1;

   if (ioData->linearizedData.type == LinearizedData::FORCED)
     com->fprintf(stderr, " ... Running Forced Oscillations w/no initializations\n");
   else if (ioData->problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_)
     com->fprintf(stderr, " ... Running Unsteady Linearized w/no structural initializations\n");
   else if (ioData->linearizedData.initCond == LinearizedData::DISPLACEMENT)  {
     com->fprintf(stderr, " ... Initializing with Displacement Mode: %d\n", ioData->linearizedData.modeNumber);
     delU[modeNum] = ioData->linearizedData.amplification;
   }
   else if (ioData->linearizedData.initCond == LinearizedData::VELOCITY)  {
     com->fprintf(stderr, " ... Initializing with Velocity Mode: %d\n", ioData->linearizedData.modeNumber);
     delY[modeNum] = ioData->linearizedData.amplification;
   }
   else  {
     com->fprintf(stderr, "*** Error: Bad Initial Condition: %d\n", ioData->linearizedData.initCond);
     exit (-1);
   }

   int dummySnap = 0;
   timeIntegrate(snaps, nSteps, 0, delU, delY, dummySnap, sdt);
 }

 if (ioData->problem.alltype == ProblemData::_ROM_ || 
     ioData->problem.alltype == ProblemData::_ROM_AEROELASTIC_)  {

   int nPodVecs = ioData->linearizedData.numPOD;
   VecSet<DistSVec<double, dim> > podVecs(nPodVecs, domain.getNodeDistInfo());
   readPodVecs(podVecs, nPodVecs);

   // print out ROM if requested
   if (ioData->output.transient.romFile[0] != 0)  {

     if (ioData->problem.alltype == ProblemData::_ROM_AEROELASTIC_)  {
       VecSet<Vec<double> > outputRom(nPodVecs, nStrMode);
       formOutputRom(outputRom, podVecs, nPodVecs);

       evalAeroSys(outputRom, podVecs, nPodVecs);
     }
     else
       evalFluidSys(podVecs, nPodVecs);
   }

   // timeIntegrate ROM
   int modeNum = ioData->linearizedData.modeNumber - 1;
   if (ioData->problem.alltype == ProblemData::_ROM_)
     com->fprintf(stderr, " ... Running Unsteady Linearized ROM w/no structural initializations\n");
   else if (ioData->linearizedData.initCond == LinearizedData::DISPLACEMENT)  {
     com->fprintf(stderr, " ... Initializing with Displacement Mode: %d\n", ioData->linearizedData.modeNumber);
     delU[modeNum] = ioData->linearizedData.amplification;
   }
   else if (ioData->linearizedData.initCond == LinearizedData::VELOCITY)  {
     com->fprintf(stderr, " ... Initializing with Velocity Mode: %d\n", ioData->linearizedData.modeNumber);
     delY[modeNum] = ioData->linearizedData.amplification;
   }
   else  {
     com->fprintf(stderr, "*** Error: Bad Initial Condition: %d\n", ioData->linearizedData.initCond);
     exit (-1);
   }

   // Form ROM operators
   VecSet<Vec<double> > ecVecs(nStrMode, nPodVecs);
   VecSet<Vec<double> > gVecs(nStrMode, nPodVecs);

   VecSet<Vec<double> > romOperator0(nPodVecs, nPodVecs);
   double *romOperator = new double[nPodVecs*nPodVecs];
   double *romOperator1 = new double[nPodVecs*nPodVecs];
   double *romOperator2 = new double[nPodVecs*nPodVecs];

   double t0 = modalTimer->getTime();
   constructROM2(romOperator, romOperator0, romOperator1, romOperator2, ecVecs, gVecs, podVecs, nPodVecs);
   modalTimer->addRomConstrTime(t0);

   modalTimer->getTime(); 
   timeIntegrateROM(romOperator, romOperator0, romOperator1, romOperator2, ecVecs, gVecs, podVecs, nSteps, nPodVecs, delU, delY, sdt);
   modalTimer->addRomTimeIntegTime(t0);

    delete [] romOperator;
    delete [] romOperator1;
    delete [] romOperator2;

 }
  delete [] delY;
  delete [] delU;
}

//----------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::createPODInTime()  {

 if (nStrMode == 0) return;

 Timer *modalTimer = domain.getTimer();
 double t0 = modalTimer->getTime();
#ifdef DO_MODAL
 double sdt = ioData->linearizedData.stepsize;
 int nSteps = ioData->ts.maxIts;
 int numSnapsTot = 0;

 int iSnap = 0;

 VecSet<DistSVec<double, dim> > snaps(numSnapsTot, domain.getNodeDistInfo());

 //  set up Snapshot Matrix
 int nModes = nStrMode;
 int numSnapsPerImpulse = nSteps;
 int snapFreq = nSteps / numSnapsPerImpulse;
 //char *suffix = "mode";
 numSnapsTot = 2*nModes*numSnapsPerImpulse;
 snaps.resize(numSnapsTot);
 com->fprintf(stderr, " ... Allocating for %d total snapshots\n", numSnapsTot);

 // Time Integration Loop for EC impulses
 int i;

 double *delU = new double[nModes];
 double *delY = new double[nModes];

 for (i = 0; i < nModes; i++)  {
   com->fprintf(stderr, " ... Impulsing mode %d\n", i);
   //delY[i] = sqrt(K[i]) * ioData->linearizedData.amplification;
   delY[i] = ioData->linearizedData.amplification;
   timeIntegrate(snaps, nSteps, snapFreq, delU, delY, iSnap, sdt);
   com->fprintf(stderr, " ... Completed Snapshot %d\n", iSnap);
   delY[i] = 0.0;

   delU[i] = ioData->linearizedData.amplification;
   timeIntegrate(snaps, nSteps, snapFreq, delU, delY, iSnap, sdt);
   com->fprintf(stderr, " ... Completed Snapshot %d\n", iSnap);
   delU[i] = 0.0;
 }
 t0 = modalTimer->getTime();

 int nPOD = ioData->linearizedData.numPOD;
 int nconv = 0;

 com->fprintf(stderr, " ... Forming Correlation Matrix\n");
 // allocate for upper half of sym. eigprob
 double *rVals = new double[numSnapsTot*(numSnapsTot+1)/2];
 for (i = 0; i < numSnapsTot; i++)
   for (int j = 0; j <= i; j++)
     rVals[(i+1)*i/2 + j] = snaps[j] * snaps[i];

 double tolerance = ioData->linearizedData.tolerance;

 ARdsSymMatrix<double> pod(numSnapsTot, rVals, 'U');

 com->fprintf(stderr, " ... Factoring Correlation Matrix\n");

 pod.FactorA();
 int ncv = 2*nPOD;
 if (ncv > numSnapsTot)  ncv = numSnapsTot-1;
 ARluSymStdEig<double> podEigProb(nPOD, pod, "LM", ncv, tolerance, 300*nPOD);

 com->fprintf(stderr, " ... Solving EigenProblem\n");
 nconv = podEigProb.FindEigenvectors();

 com->fprintf(stderr, " ... Got %d converged eigenvectors out of %d snaps\n", nconv, numSnapsTot);
 outputPODVectors(podEigProb, snaps, nPOD, numSnapsTot);

#else
 com->fprintf(stderr, "  ... ERROR: REQUIRES COMPILATION WITH ARPACK and DO_MODAL Flag\n");
 exit(-1);
#endif

 modalTimer->addMeshSolutionTime(t0);
}

//------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::timeIntegrate(VecSet<DistSVec<double, dim> > &snaps, 
                      int nSteps, int snapFreq, double *delU, double *delY, 
                      int &iSnap, double sdt, char *snapFile)  {

 Timer *modalTimer = domain.getTimer();
 double t0;

 DistSVec<double,dim> FF(domain.getNodeDistInfo());
 VarFcn *varFcn = new VarFcn(*ioData);  
 spaceOp->computeResidual(Xref, controlVol, Uref, FF, tState);

 // basic initializations
 DistSVec<double,3> deltmp(domain.getNodeDistInfo());
 DistSVec<double,dim> delW(domain.getNodeDistInfo());
 DistSVec<double,dim> delWtmp(domain.getNodeDistInfo());
 DistSVec<double,dim> rhs(domain.getNodeDistInfo());
 DistSVec<double,dim> rhsA(domain.getNodeDistInfo());
 DistSVec<double,dim> rhsB(domain.getNodeDistInfo());
 DistSVec<double,dim> rhsC(domain.getNodeDistInfo());
 DistSVec<double,dim> rhsD(domain.getNodeDistInfo());
 DistSVec<double,dim> delWnm1(domain.getNodeDistInfo());
 DistSVec<double,dim> delWint(domain.getNodeDistInfo());

 // Uncomment for forced oscillations
 double freq = 2.0*pi*ioData->linearizedData.frequency;
 double dMax = ioData->linearizedData.amplification;

 int i;
 deltmp = 0.0;

 double sdt0 = sdt*dt0/dt;
 double *prevU = new double[nStrMode];
 double *prevY = new double[nStrMode];
 double *prevA = new double[nStrMode];

 for (i = 0; i < nStrMode; ++i)  {
   deltmp += (delU[i]+sdt*delY[i])*mX[i];
   prevU[i] = delU[i];
   prevY[i] = delY[i];
 }

 deltmp += Xref;

 // Time Loop
 ksp->printParam();
 ksp2->printParam();
 ksp3->printParam(); 

 int printFreq = nSteps / 20;
 if (nSteps < 20)
   printFreq = 1;

 com->fprintf(stderr, " ... Doing Time Integration with refVals: F= %e, rho = %e, P=%e, M=%e\n", ioData->ref.rv.force, ioData->ref.density, ioData->ref.pressure, ioData->ref.mach); 

 // Init delW
 char *nlSolFile = 0;
 if (ioData->input.perturbed[0] == 0)
   nlSolFile = tInput->solutions;
 else  {
   int sp = strlen(ioData->input.prefix) + 1;
   nlSolFile = new char[sp + strlen(ioData->input.perturbed)];
   sprintf(nlSolFile, "%s%s", ioData->input.prefix, ioData->input.perturbed);
 }
 domain.readVectorFromFile(nlSolFile, 0, 0, delW);
 com->fprintf(stderr, " ... Read Perturbed solution: W = %e\n", delW.norm());

 tOutput->writeForcesToDisk(0, 0, 0, 0, 0.0, com->cpuNum(), tRestart->energy, deltmp, delW);
 tOutput->writeBinaryVectorsToDisk(false, 0, 0, deltmp, controlVol, delW, tState);

 delW -= Uref;
 delWnm1 = delW;
 delWint = delW;

 // Compute Reference Modal Force
 DistSVec<double,3> refNodalForce(domain.getNodeDistInfo());
 DistVec<double> Pin(domain.getFaceDistInfo());
 Vec<double> refModalF(nStrMode);
 com->fprintf(stderr, " ... Computing Ref Nodal Force w/pressure: %e\n", ioData->aero.pressure);
 Pin = ioData->aero.pressure;
 postOp->computeNodalForce(Xref, Uref, Pin, refNodalForce);

 CommPattern<double> *vpat = domain.getVecPat();
 domain.assemble(vpat, refNodalForce);

 for (i= 0; i < nStrMode; i++)
   refModalF[i] = mX[i]*refNodalForce;

  //Compute initial acceleration 
  Vec<double> modalF(nStrMode);
  postOp->computeForceDerivs(Xref, Uref, delW, modalF, mX);
  modalF += refModalF;
  modalF *= ioData->ref.rv.force;

  for (i = 0; i < nStrMode; i++)
    prevA[i] = modalF[i] - K[i]*prevU[i];

  FILE *dispFP;
  if (ioData->output.transient.gendispFile[0] != 0)  {
    // generalized displacements output
    int sp = strlen(ioData->output.transient.prefix);
    char *dispFile = new char[sp + strlen(ioData->output.transient.gendispFile)+1];
    sprintf(dispFile, "%s%s", ioData->output.transient.prefix, ioData->output.transient.gendispFile);
    dispFP = fopen(dispFile, "w");
    com->barrier();
    if (ioData->problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_AEROELASTIC_)  {
   /* if (!dispFP)  {
      com->fprintf(stderr, "*** Warning: Cannot create generalized displacement FILE in %s\n", dispFile);
      //exit (-1);
    }*/
      com->fprintf(dispFP,"%d %f ",0,0.0);
      for (i=0; i < nStrMode; ++i) {
        com->fprintf(dispFP,"%.16e ",delU[i]);
      }
      com->fprintf(dispFP,"\n");
    }
  }

 int cntr = 0;
 int cntp1;

 // every step corresponds to solving for 2x each member of the BDF scheme //
 
 FF *= (-2.0); 

 for (int cnt = 0; cnt < nSteps; ++cnt) {

   cntp1 = cnt+1;

   t0 = modalTimer->getTime();  
   if (cnt == 0) {

     rhsA = 0.0;
     rhsB = 0.0;
     rhsC = 0.0;
     rhsD = 0.0;
      rhs  = 0.0;
     for (i = 0; i < nStrMode; ++i)  {
       rhsA += DX[i]*delU[i];
       rhsB += DX[i]*delY[i];
       rhsC += DE[i]*delY[i];
       rhsD += DE[i]*prevA[i];
     }

     //Predictor step
      HOp2step1->apply(delWint,rhs);
      rhs+= rhsA + rhsC - 0.5*FF;
      rhs*= 0.5*dt;
      rhs /= controlVol;
      delW = delWint - rhs;
     
      rhs = 0.0;
      //Corrector step 
      HOp2step1->apply(delW,rhs);
      rhs+= rhsA + 0.5*sdt*rhsB + rhsC + 0.5*sdt*rhsD - 0.5*FF;
      rhs*= dt;
      rhs /= controlVol;
      delW = delWint - rhs;
   }
   else if (cnt == 1) {
      
     rhs = delW*(6.0/dt) - delWnm1*(2.0/(3.0*dt));
     rhs *= controlVol;

     rhsA = 0.0;
     rhsB = 0.0;
     rhsC = 0.0;
     rhsD = 0.0;

     for (i = 0; i < nStrMode; ++i)  {
        rhsA += DX[i]*delU[i];
        rhsB += DX[i]*delY[i];
        rhsC += DE[i]*delY[i];
        rhsD += DE[i]*prevY[i];
     }
 
     rhsA *= 2.0;
     rhsB *= sdt;
     rhsC *= 3.0;

     rhs -= rhsA;
     rhs -= rhsB;
     rhs -= rhsC;
     rhs += rhsD;
      
     delWnm1 = delW;
   
     rhs += FF;

     ksp2->solve(rhs, delW);
     modalTimer->addKspTime(t0);
   }
   else if (cnt == 2) {
  
     rhs = delW*(6.0/dt) - delWnm1*(8.0/(3.0*dt));
     rhs *= controlVol;

     rhsA = 0.0;
     rhsB = 0.0;
     rhsC = 0.0;
     rhsD = 0.0;

     for (i = 0; i < nStrMode; ++i)  {
        rhsA += DX[i]*delU[i];
        rhsB += DX[i]*delY[i];
        rhsC += DE[i]*delY[i];
        rhsD += DE[i]*prevY[i];
     }
 
     rhsA *= 2.0;
     rhsB *= sdt;
     rhsC *= 3.0;

     rhs -= rhsA;
     rhs -= rhsB;
     rhs -= rhsC;
     rhs += rhsD;
      
     delWnm1 = delW;
   
     rhs += FF;

     ksp3->solve(rhs, delW);
     modalTimer->addKspTime(t0);
   }

   else {
     
      rhs = delW*(4.0/dt) - delWnm1*(1.0/dt);
      rhs *= controlVol;

      rhsA = 0.0;     
      rhsB = 0.0;     
      rhsC = 0.0;     
      rhsD = 0.0;     

     for (i = 0; i < nStrMode; ++i)  {
         rhsA += DX[i]*delU[i];
         rhsB += DX[i]*delY[i];
         rhsC += DE[i]*delY[i];
         rhsD += DE[i]*prevY[i];
     }
 
     rhsA *= 2.0;
     rhsB *= sdt;
     rhsC *= 3.0;

     rhs -= rhsA;
     rhs -= rhsB;
     rhs -= rhsC;
     rhs += rhsD;
      
     delWnm1 = delW;
   
     rhs += FF;

     ksp->solve(rhs, delW);
     modalTimer->addKspTime(t0);
   }

   // for forced oscillations          
   if (ioData->linearizedData.type == LinearizedData::FORCED)  {
     if (cnt <= 1) {
        delU[ioData->linearizedData.modeNumber-1] = sin(freq*(cnt)*sdt)*dMax;   
        delY[ioData->linearizedData.modeNumber-1] = freq*cos(freq*(cnt)*sdt)*dMax;
     } 
     else {
        delU[ioData->linearizedData.modeNumber-1] = sin(freq*((cnt)*sdt))*dMax;
        delY[ioData->linearizedData.modeNumber-1] = freq*cos(freq*((cnt)*sdt))*dMax;  
     }
   }

   if (ioData->problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_AEROELASTIC_)  {

     for (int kk = 0; kk < nStrMode; kk++)  {
       prevU[kk] = delU[kk];
       prevY[kk] = delY[kk];
     }
     t0 = modalTimer->getTime();

     computeModalDisp(sdt, deltmp, delW, delU, delY, refModalF, cnt);
     modalTimer->addStructUpdTime(t0);

     if (ioData->output.transient.gendispFile[0] != 0)  {
       // output generalized displacements
       com->fprintf(dispFP, "%d %f ",cntp1, (cnt+1)*sdt);
       for (i = 0; i < nStrMode; ++i) {
         com->fprintf(dispFP, "%.16e ", delU[i]);
       }
       com->fprintf(dispFP, "\n");
     }
   }
   // compute updated position
   deltmp = 0.0;
   for (i = 0; i < nStrMode; ++i) {
       deltmp += (delU[i]+sdt/2*delY[i])*mX[i];
   }
   deltmp += Xref;

   if (ioData->problem.alltype == ProblemData::_POD_CONSTRUCTION_)  {
     if (cnt % snapFreq == 0 || cnt == 0) {
       snaps[iSnap] = delW;
       if (snapFile) domain.writeVectorToFile(snapFile, iSnap, 0, snaps[iSnap]);
       iSnap++;
     }
   }

   delWtmp = Uref + delW;
   int ierr = domain.checkSolution(varFcn, delWtmp);

    tOutput->writeForcesToDisk(cntp1, cntp1, cntp1, cntp1, (cnt+1)*dt, com->cpuNum(), tRestart->energy, deltmp, delWtmp);
    tOutput->writeBinaryVectorsToDisk(false, cntp1, (cnt+1)*dt, deltmp, controlVol, delWtmp, tState);

   if (ierr)  {
     com->fprintf(stderr, " ... WARNING: %d nodes have neg. rho/P \n", ierr);
     cntr++;
     /*if (cntr > 10)  {
       com->barrier();
       exit(-1);
     }*/
   }

    if (cnt % 20 == 0)
      com->fprintf(stderr, " ... Iteration: %d   Time: %f\n", cnt, cnt*sdt);   

   if (ioData->problem.alltype == ProblemData::_UNSTEADY_LINEARIZED_ || 
       ioData->problem.alltype == ProblemData::_POD_CONSTRUCTION_)  {
     for (i = 0; i < nStrMode; ++i) {
       delU[i] = 0.0;
       delY[i] = 0.0;
     }
   }
 }

  for (i = 0; i < nStrMode; ++i) {
    delU[i] = 0.5*(delU[i] + prevU[i]);
    delY[i] = 0.5*(delY[i] + prevY[i]);
  }
  if (!(ioData->input.perturbed[0] == 0))
    delete [] nlSolFile;
  delete [] prevU;
  delete [] prevY;
  delete [] prevA;
}

//----------------------------------------------------------------------------------

template<int dim>
void
ModalSolver<dim>::timeIntegrateROM(double *romOp, VecSet<Vec<double> > &romOp0, double *romOp1, double *romOp2, VecSet<Vec<double> > &ecMat, VecSet<Vec<double> > &gMat, VecSet<DistSVec<double, dim> > &podVecs, int nSteps, int nPodVecs, double *delU, double *delY, double sdt)  {

#ifdef DO_MODAL
  VarFcn *varFcn = new VarFcn(*ioData);  

  // basic initializations
  DistSVec<double,3> deltmp(domain.getNodeDistInfo());
  DistSVec<double, dim> delWFull(domain.getNodeDistInfo());
  Vec<double> delWRom(nPodVecs);
  Vec<double> delWRomTemp(nPodVecs);
  Vec<double> prevWRom(nPodVecs);
  Vec<double> pprevWRom(nPodVecs);
  Vec<double> modalF(nStrMode);
  modalF = 0.0;

  delWRom = 0.0;
  deltmp = 0.0;

  double *prevU = new double[nStrMode];
  double *prevY = new double[nStrMode];

  int i;
  for (i = 0; i < nStrMode; ++i)
    deltmp += (delU[i]+0.5*sdt*delY[i])*mX[i];

  deltmp += Xref;

  // Init delW
  char *nlSolFile = 0;
  if (ioData->input.perturbed[0] == 0)
    nlSolFile = tInput->solutions;
  else  {
    int sp = strlen(ioData->input.prefix) + 1;
    nlSolFile = new char[sp + strlen(ioData->input.perturbed)+1];
    sprintf(nlSolFile, "%s%s", ioData->input.prefix, ioData->input.perturbed);
  }
  domain.readVectorFromFile(nlSolFile, 0, 0, delWFull);
  com->fprintf(stderr, " ... Read Perturbed solution: W = %e, Uref = %e\n", delWFull.norm(), Uref.norm());
  delWFull -= Uref;

  // construct initial delWRom
  if (ioData->ts.form == TsData::DESCRIPTOR) {
    DistSVec<double, dim> temp(domain.getNodeDistInfo());
    temp = delWFull;
    temp *= controlVol;
    for (i = 0; i < nPodVecs; i++)
      delWRom[i] = podVecs[i] * temp;
  }
  else {
    for (i = 0; i < nPodVecs; i++)
      delWRom[i] = podVecs[i] * delWFull;
  }

  //ROM Initial Condition Output
  FILE *romICFP;
  if (ioData->output.transient.romInitialConditionFile[0] != 0)  {
    int sp = strlen(ioData->output.transient.prefix);
    char *romICFile = new char[sp + strlen(ioData->output.transient.romInitialConditionFile)+1];
    sprintf(romICFile, "%s%s", ioData->output.transient.prefix, ioData->output.transient.romInitialConditionFile);
    romICFP = fopen(romICFile, "w");
    com->barrier();
    com->fprintf(romICFP,"%d\n",nPodVecs+2*nStrMode);
    for (i=0; i < nPodVecs; ++i)
      com->fprintf(romICFP,"%.16e\n", delWRom[i]);
    for (i=0; i < nStrMode; ++i)
      com->fprintf(romICFP,"%.16e\n", delY[i]);
    for (i=0; i < nStrMode; ++i)
      com->fprintf(romICFP,"%.16e\n", delU[i]);
    fclose(romICFP);
  }

#ifndef NDEBUG
  com->fprintf(stderr, " ... Initial Condition for Fluid ROM:\n");
  com->fprintf(stderr," q0 = [");
  for (i=0; i < nPodVecs; i++)
    com->fprintf(stderr,"%e ",delWRom[i]); 
  com->fprintf(stderr,"];\n");
#endif

 // Compute Reference Modal Force
 DistSVec<double,3> refNodalForce(domain.getNodeDistInfo());
 DistVec<double> Pin(domain.getFaceDistInfo());
 Vec<double> refModalF(nStrMode);
 com->fprintf(stderr, " ... Computing Ref Nodal Force w/pressure: %e\n", ioData->aero.pressure);
 Pin = ioData->aero.pressure;
 postOp->computeNodalForce(Xref, Uref, Pin, refNodalForce);

 CommPattern<double> *vpat = domain.getVecPat();
 domain.assemble(vpat,refNodalForce);

 // Form ROM operators
 VecSet<Vec<double> > ecOpMat(nStrMode, nPodVecs);
 VecSet<Vec<double> > gOpMat(nStrMode, nPodVecs);
 VecSet<Vec<double> > ecOpMat1(nStrMode, nPodVecs);
 VecSet<Vec<double> > gOpMat1(nStrMode, nPodVecs);
 VecSet<Vec<double> > ecOpMat2(nStrMode, nPodVecs);
 VecSet<Vec<double> > gOpMat2(nStrMode, nPodVecs);

 VecSet<Vec<double> > romOpMinus(nPodVecs, nPodVecs);
 VecSet<Vec<double> > romOperator(nPodVecs, nPodVecs);
 VecSet<Vec<double> > romOperator1(nPodVecs, nPodVecs);
 VecSet<Vec<double> > romOperator2(nPodVecs, nPodVecs);

 int iVec;
 for (iVec = 0; iVec < nPodVecs; iVec++)  {
   romOpMinus[iVec] = 0.0;
   romOpMinus[iVec][iVec] = 1.0;
 } 

 ARdsNonSymMatrix<double, double> romOpPlus(nPodVecs, romOp);
 romOpPlus.FactorA();
 ARdsNonSymMatrix<double, double> romOpPlus1(nPodVecs, romOp1);
 romOpPlus1.FactorA();
 ARdsNonSymMatrix<double, double> romOpPlus2(nPodVecs, romOp2);
 romOpPlus2.FactorA();
 com->fprintf(stderr, " ... Factored Matrix\n");

 com->fprintf(stderr, " ... Forming Reduced ROM Operator\n");

 for (iVec = 0; iVec < nPodVecs; iVec++) {
   romOpPlus.MultInvv(romOpMinus[iVec].data(), romOperator[iVec].data());
   romOpPlus1.MultInvv(romOpMinus[iVec].data(), romOperator1[iVec].data());
   romOpPlus2.MultInvv(romOpMinus[iVec].data(), romOperator2[iVec].data());
 }

 for (iVec = 0; iVec < nStrMode; iVec++)  {
   romOpPlus.MultInvv(ecMat[iVec].data(), ecOpMat[iVec].data());
   romOpPlus.MultInvv(gMat[iVec].data(), gOpMat[iVec].data());
   romOpPlus1.MultInvv(ecMat[iVec].data(), ecOpMat1[iVec].data());
   romOpPlus1.MultInvv(gMat[iVec].data(), gOpMat1[iVec].data());
   romOpPlus2.MultInvv(ecMat[iVec].data(), ecOpMat2[iVec].data());
   romOpPlus2.MultInvv(gMat[iVec].data(), gOpMat2[iVec].data());
 }

 for (i = 0; i < nStrMode; i++)
   refModalF[i] = mX[i]*refNodalForce;
 com->fprintf(stderr, " ... Norm RefModalF = %e\n", refModalF.norm());
 com->fprintf(stderr, " refModalF = [");
 for (i = 0; i < nStrMode; i++)
   com->fprintf(stderr, "%e ",ioData->ref.rv.force*refModalF[i]);
 com->fprintf(stderr, "];\n");


 com->fprintf(stderr, " ... Initial delWRom: %e\n", delWRom.norm());

 int cntp1;
 delWFull = 0.0;
 for (i = 0; i < nPodVecs; ++i)
   delWFull += podVecs[i] * delWRom[i];
 delWFull += Uref;
 tOutput->writeForcesToDisk(0, 0, 0, 0, 0.0, com->cpuNum(), tRestart->energy, deltmp, delWFull);
 tOutput->writeBinaryVectorsToDisk(false, 0, 0, deltmp, controlVol, delWFull, tState);

  pprevWRom = delWRom; 
  prevWRom = delWRom;

 
  //Precompute P*PHI
  VecSet<Vec<double> > PtimesPhi(nPodVecs, nStrMode);

  for (int iVec = 0; iVec < nPodVecs; iVec++)  {
    postOp->computeForceDerivs(Xref, Uref, podVecs[iVec], modalF, mX);
    PtimesPhi[iVec] = ioData->ref.rv.force*modalF;
  }

  // generalized displacements output
  FILE *dispFP;
  if (ioData->output.transient.gendispFile[0] != 0)  {

    int sp = strlen(ioData->output.transient.prefix);
    char *dispFile = new char[sp + strlen(ioData->output.transient.gendispFile)+1];
    sprintf(dispFile, "%s%s", ioData->output.transient.prefix, ioData->output.transient.gendispFile);
    dispFP = fopen(dispFile, "w");
    com->barrier();
    com->fprintf(dispFP,"%d %f ",0,0.0);
    for (i=0; i < nStrMode; ++i) {
      com->fprintf(dispFP,"%.16e ", delU[i]);
    }
    com->fprintf(dispFP,"\n");
  }

  //Time integration loop
  for (int cnt = 0; cnt < nSteps+1; ++cnt) {

    cntp1 = cnt+1;

    if (cnt == 0){
	    delWRomTemp = prevWRom;
	    for (i = 0; i < nPodVecs; ++i)
        delWRomTemp += 0.5*dt*delWRom[i]*romOp0[i];

      for (i = 0; i < nStrMode; ++i){
        delWRomTemp -= 0.5*dt*(delU[i]*gMat[i] + delY[i]*ecMat[i]);
        delWRom -= dt*((delU[i] + 0.5*sdt*delY[i])*gMat[i] + delY[i]*ecMat[i]);
      }

     for (i = 0; i < nPodVecs; ++i)
        delWRom += dt*delWRomTemp[i]*romOp0[i];
   }
  
    else if (cnt == 1) {
     delWRom = 0.0;

     for (i = 0; i < nPodVecs; ++i)
      delWRom += (1.0/dt)*(3.0*prevWRom[i] - (1.0/3.0)*pprevWRom[i])*romOperator[i];

      for (i = 0; i < nStrMode; ++i)
        delWRom -= ( (delU[i] + 0.5*sdt*delY[i])*gOpMat[i] + (1.5*delY[i] - 0.5*prevY[i])*ecOpMat[i] );
    } 
 
    else if (cnt == 2) {
      delWRom = 0.0;
     for (i = 0; i < nPodVecs; ++i)
        delWRom += (1.0/dt)*(3.0*prevWRom[i] - (4.0/3.0)*pprevWRom[i])*romOperator1[i];
 
     for (i = 0; i < nStrMode; ++i)
        delWRom -= ( (delU[i] + 0.5*sdt*delY[i])*gOpMat1[i] + (1.5*delY[i] - 0.5*prevY[i])*ecOpMat1[i] );
   }

   else {
      delWRom = 0.0;
      for (i = 0; i < nPodVecs; ++i)
        delWRom += (1.0/dt)*(2.0*prevWRom[i] - 0.5*pprevWRom[i])*romOperator2[i];
 
      for (i = 0; i < nStrMode; ++i)
        delWRom -= ( (delU[i] + 0.5*sdt*delY[i])*gOpMat2[i] + (1.5*delY[i] - 0.5*prevY[i])*ecOpMat2[i] );
    }

    pprevWRom = prevWRom;
    prevWRom = delWRom;

    // project soltn into full space
    delWFull = 0.0;
    for (i = 0; i < nPodVecs; ++i)
      delWFull += podVecs[i] * delWRom[i];


   if (ioData->problem.alltype == ProblemData::_ROM_AEROELASTIC_){
     for (int kk=0; kk < nStrMode; kk++){
       prevU[kk] = delU[kk];
       prevY[kk] = delY[kk];
     }
   }

    // Output the reduced order vector
    computeModalDisp(sdt, delWRom, delU, delY, refModalF, PtimesPhi, nPodVecs, cnt);

   if (ioData->output.transient.gendispFile[0] != 0)  {
     // output generalized displacements
     com->fprintf(dispFP, "%d %f ",cntp1, (cnt+1)*sdt);
     for (i = 0; i < nStrMode; ++i) {
       com->fprintf(dispFP, "%.16e ", delU[i]);
     }
     com->fprintf(dispFP, "\n");
   }
   
   // compute Cl, Cm
   deltmp = 0.0;
   for (i = 0; i < nStrMode; ++i)
     deltmp += (delU[i]+sdt/2.0*delY[i])*mX[i];

   deltmp += Xref;
   delWFull += Uref;

   int ierr = domain.checkSolution(varFcn, delWFull);

    tOutput->writeForcesToDisk(cntp1, cntp1, cntp1, cntp1, (cnt+1)*dt, com->cpuNum(), tRestart->energy, deltmp, delWFull);
    tOutput->writeBinaryVectorsToDisk(false, cntp1, (cnt+1)*dt, deltmp, controlVol, delWFull, tState);

   /*if (ierr)  {
     com->fprintf(stderr, " ... WARNING: %d nodes have neg. rho/P \n", ierr);
//     exit(-1);
   }*/

   if (cnt % 20 == 0)
     com->fprintf(stderr, " ... Iteration: %d   Time: %f\n", cnt, cnt*sdt);
  }

  for (i = 0; i < nStrMode; ++i) {
    delU[i] = 0.5*(prevU[i] + delU[i]);
    delY[i] = 0.5*(prevY[i] + delY[i]);
  } 
  delete [] prevU;
  delete [] prevY;
  if (!(ioData->input.perturbed[0] == 0))
    delete [] nlSolFile;

#else
 com->fprintf(stderr, "*** Error: REQUIRES COMPILATION WITH ARPACK and DO_MODAL Flag\n");
 exit(-1);

#endif
}

//--------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::preProcess()  {

 // setup solvers
 VarFcn *varFcn = new VarFcn(*ioData);  
 geoState = new DistGeoState(*ioData, &domain);
 geoState->setup1(tInput->positions, &Xref, &controlVol);
 bcData = new DistBcDataEuler<dim>(*ioData, varFcn, &domain, Xref);

 spaceOp = new SpaceOperator<dim>(*ioData, varFcn, bcData, geoState, &domain);
 postOp = new PostOperator<dim> (*ioData, varFcn, bcData, geoState, &domain);

 // Temporal operator contains Uref for us
 tState = new DistTimeState<dim>(*ioData, spaceOp, varFcn, &domain);
 tState->setup(tInput->solutions,  Xref, bcData->getInletBoundaryVector(), Uref, *ioData); 
 RefVal *refVal = new RefVal(ioData->ref.rv);
 tOutput = new TsOutput<dim>(*ioData, refVal, &domain, postOp);
 tRestart = new TsRestart(*ioData, refVal);

 HOp = new MatVecProdH2<dim,double, dim>(*ioData,  varFcn, tState, spaceOp, &domain); 
 HOp2 = new MatVecProdH2<dim,double,dim>(*ioData,  varFcn, tState, spaceOp, &domain); 
 HOp2step1 = new MatVecProdH2<dim,double, dim>(*ioData,  varFcn, tState, spaceOp, &domain); 
 HOpstep2 = new MatVecProdH2<dim,double,dim>(*ioData,  varFcn, tState, spaceOp, &domain); 
 HOpstep3 = new MatVecProdH2<dim,double,dim>(*ioData,  varFcn, tState, spaceOp, &domain); 

 // Stuff for solver
 PcData pcData = ioData->ts.implicit.newton.ksp.ns.pc;
 KspData &kspData = ioData->ts.implicit.newton.ksp.ns;

 if (pcData.type == PcData::IDENTITY)
   pc = new IdentityPrec<dim, double>();
 else if (pcData.type == PcData::JACOBI)
   pc = new JacobiPrec<PreScalar,dim>(DiagMat<PreScalar,dim>::DENSE, &domain);
 else if (pcData.type == PcData::AS || pcData.type == PcData::RAS ||
          pcData.type == PcData::ASH || pcData.type == PcData::AAS)
   pc = new IluPrec<PreScalar,dim>(pcData, &domain);


 if (kspData.type == KspData::GMRES) {
   ksp = new GmresSolver<DistSVec<double,dim>, MatVecProd<dim, dim>,
       KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOp, pc, com); 
   ksp2 = new GmresSolver<DistSVec<double,dim>, MatVecProd<dim, dim>,
       KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOpstep2, pc, com); 
   ksp3 = new GmresSolver<DistSVec<double,dim>, MatVecProd<dim, dim>, 
       KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOpstep3, pc, com); 
   }
 else if (kspData.type == KspData::GCR) {
   ksp = new GcrSolver<DistSVec<double,dim>, MatVecProd<dim, dim>,
       KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOp, pc, com); 
   ksp2 = new GcrSolver<DistSVec<double,dim>, MatVecProd<dim, dim>,
       KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOpstep2, pc, com); 
   ksp3 = new GcrSolver<DistSVec<double,dim>, MatVecProd<dim, dim>,
       KspPrec<dim, double>, Communicator>((&domain)->getNodeDistInfo(), kspData, HOpstep3, pc, com); 
 }

 if (ioData->linearizedData.domain == LinearizedData::FREQUENCY)  {
   HOpC = new MatVecProdH2<dim,bcomp,dim>(*ioData,  varFcn, tState, spaceOp, &domain); 

   pcComplex = new IluPrec<bcomp ,dim, bcomp>(pcData, &domain);
  if (ioData->linearizedData.padeReconst == LinearizedData::TRUE) {


    kspCompGcr = new GcrSolver<DistSVec<bcomp,dim>, MatVecProd<dim,dim>,
                   KspPrec<dim, bcomp>, Communicator, bcomp>
                   ((&domain)->getNodeDistInfo(), kspData, HOpC, pcComplex, com); 
   }
   else {

     if (kspData.type == KspData::GMRES)

       kspComp = new GmresSolver<DistSVec<bcomp,dim>, MatVecProd<dim, dim>,
                     KspPrec<dim, bcomp>, Communicator, bcomp>
                    ((&domain)->getNodeDistInfo(), kspData, HOpC, pcComplex, com); 
     else if (kspData.type == KspData::GCR)

       kspComp = new GcrSolver<DistSVec<bcomp,dim>, MatVecProd<dim, dim>,
                     KspPrec<dim, bcomp>, Communicator, bcomp>
                    ((&domain)->getNodeDistInfo(), kspData, HOpC, pcComplex, com); 


   }
 }

 geoState->setup2(tState->getData());

 //Setup Time state
 double dummyTime = 0.0;
 int dummySubCycles = 1;
 tState->computeTimeStep(1.0, 1.0, &dummyTime, &dummySubCycles, *geoState, Xref, controlVol, Uref);
 dt = ioData->linearizedData.stepsize/ioData->ref.rv.time;    //ts.timestep/  (ref.length/velocity)

 //time step for the first iteration (Forward Euler) has to be much smaller : here sdt0 = sdt^2

 dt0 = (ioData->linearizedData.stepsize)*(ioData->linearizedData.stepsize)/ioData->ref.rv.time; 
 com->fprintf(stderr, " ... Running Modal Solver w/%d struct modes, fluid step: %f, RefTime = %f (%f/%f)\n", nStrMode, dt, ioData->ref.rv.time, ioData->ref.rv.length, ioData->ref.rv.velocity);

 geoState->compute(tState->getData(), bcData->getVelocityVector(), Xref, controlVol);
 geoState->update(Xref, controlVol);
 geoState->compute(tState->getData(), bcData->getVelocityVector(), Xref, controlVol);

 // F is not used
 DistSVec<double,dim> FF(domain.getNodeDistInfo());
 spaceOp->computeResidual(Xref, controlVol, Uref, FF, tState);
 FF /= controlVol;

 com->fprintf(stderr, " ... Norm FF: %e\n", FF.norm());
 if (FF.norm() > 1e-4) {
   com->fprintf(stderr, " *** WARNING: Check Steady State Convergence for this Mach Number ***\n");
   if (FF.norm() > 1e-2) {
     com->fprintf(stderr, " *** ERROR: Residual is too high. Exiting ***\n");
     exit(-1);
   }
 }
 // Compute Derivatives***********************************
 //***Variable Initializations***
 DistSVec<double,3> Xnp1(domain.getNodeDistInfo());
 DistVec<double> CV1(domain.getNodeDistInfo());
 DistSVec<double,dim> F1(domain.getNodeDistInfo());
 DistVec<double> CV2(domain.getNodeDistInfo());
 DistSVec<double,dim> F2(domain.getNodeDistInfo());
 VecSet< DistVec<double> > DA(nStrMode, domain.getNodeDistInfo());
 DistSVec<double,dim> E(domain.getNodeDistInfo());
 VecSet< DistSVec<double,dim> > DV(nStrMode, domain.getNodeDistInfo());

 DX.resize(nStrMode);
 DE.resize(nStrMode);

 double eps = ioData->linearizedData.eps;
// double alpha = dt/ioData->linearizedData.eps2;
// double eps2 = eps*alpha;
// com->fprintf(stderr, "Alpha = %f\n", alpha);

 // Loop over the modes

 for (int ic=0;ic<nStrMode;++ic) {

   //First DFDX & DADX***********************************
   //***Computing F(X1,V=0)
   Xnp1 = Xref -  eps*mX[ic]; //1.0*Xref;

   geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV1);
   geoState->update(Xnp1, CV1);
   geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV1);

   spaceOp->computeResidual(Xnp1, CV1, Uref, F1, tState);
   spaceOp->applyBCsToResidual(Uref, F1);

   //***Computing F(X2,V=0)
   Xnp1 = Xref + eps*(mX[ic]);
   geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV2);
   geoState->update(Xnp1, CV2);
   geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV2);

   spaceOp->computeResidual(Xnp1, CV2, Uref, F2, tState);
   spaceOp->applyBCsToResidual(Uref, F2);

   DX[ic] = 1.0/(2.0*eps)*(F2 - F1);
   DA[ic] = 1.0/(2.0*eps)*(CV2 - CV1);
   com->fprintf(stderr, " ... Norm DX, Mode %i: %e\n",ic, DX[ic].norm());
   com->fprintf(stderr, " ... Norm DA, Mode %i: %e\n",ic, DA[ic].norm());

   // Now Compute w*dA/dX***********************************
   int numLocSub = Uref.numLocSub();

   #pragma omp parallel for
   for (int iSub=0; iSub<numLocSub; ++iSub) {
     double *cv = DA[ic].subData(iSub);
     double (*r)[dim] = Uref.subData(iSub);
     double (*ans)[dim] = E.subData(iSub);

     for (int i=0; i<DA[ic].subSize(iSub); ++i) {
       for (int j=0; j<dim; ++j)
         ans[i][j] = cv[i]*r[i][j];
     }
   }

   // need to mult. by time to maintain a dimensional vel. (conversion from adim time to dim time)
   DE[ic] = ioData->ref.rv.time*E;

   com->fprintf(stderr, " ... Norm DE, Mode %i: %e\n",ic, DE[ic].norm());

   // Then DFDXdot*****************************************
   // Computing F(X0,V1) 
   double Ti;
   if (K[ic]>0.0)
     Ti = (pi*2.0/sqrt(K[ic]) ) / ioData->ref.rv.time; 
   else
     Ti = 1.0;
   
   double eps2 = 0.5*dt/Ti*eps;

   Xnp1 = Xref - eps2*(mX[ic]);

   geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV1);
   geoState->update(Xnp1, CV1);

   Xnp1 = Xref + eps2*(mX[ic]);
   geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV1);

   //spaceOp->computeResidual(Xnp1, CV1, Uref, F1, tState);
   spaceOp->computeResidual(Xref, controlVol, Uref, F1, tState);
   spaceOp->applyBCsToResidual(Uref, F1);

   // Computing F(X0, -V1)
   Xnp1 = Xref + eps2*(mX[ic]);

   geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV1);
   geoState->update(Xnp1, CV1);

   Xnp1 = Xref - eps2*(mX[ic]);
   geoState->compute(tState->getData(), bcData->getVelocityVector(), Xnp1, CV1);

   //spaceOp->computeResidual(Xnp1, CV1, Uref, F2, tState);
   spaceOp->computeResidual(Xref, controlVol, Uref, F2, tState);
   spaceOp->applyBCsToResidual(Uref, F2);


   // no dt here because the spatial computation will divide the
   // given displacement change by dt to define its velocity.
   // Thus with Xnp1 = X0+dt*eps*Xm and Xn = X0-dt*eps*Xm: V_n+1/2 = 2*eps*Xm

   DV[ic] = (Ti/(2.0*eps))*(F1-F2);
   DV[ic] *= ioData->ref.rv.time;

   com->fprintf(stderr, " ... Norm DV, Mode %i: %e\n",ic, DV[ic].norm());

   // create E+C
   DE[ic] += DV[ic];
   com->fprintf(stderr, " ... Norm DE, Mode %i: %e\n\n",ic, DE[ic].norm());
 }

 // need to reset geoState after finite differences
 geoState->compute(tState->getData(), bcData->getVelocityVector(), Xref, controlVol);
 geoState->update(Xref, controlVol);
 geoState->compute(tState->getData(), bcData->getVelocityVector(), Xref, controlVol);



 // Now setup H matrices

 double delt0;
 delt0 = dt0;

 double r = dt/(2.0*delt0);
 // ***This is (c1*A+c2*dt*H)
 HOp->evaluate(3, Xref, controlVol, Uref, FF,0.0);

 //***This is (c1*A-c2*dt*H)
 HOp2->evaluate2(0, Xref, controlVol, Uref, FF);

  tState->setGlobalTimeStep(dt);
  HOp2step1->evalH(0, Xref, controlVol, Uref); 

  tState->setGlobalTimeStep(dt);
  HOpstep2->evaluate(8, Xref, controlVol, Uref, FF,0.0);

 tState->setGlobalTimeStep(dt);
  HOpstep3->evaluate(9, Xref, controlVol, Uref, FF,0.0);

 //set up preconditioner, GMRES Solver
 DistMat<PreScalar,dim> *_pc = dynamic_cast<DistMat<PreScalar,dim> *>(pc);

 if (_pc) {
   spaceOp->computeH1(Xref, controlVol, Uref, *_pc);
   tState->addToH1(controlVol, *_pc);
   spaceOp->applyBCsToJacobian(Uref, *_pc);
 }

 pc->setup();
 ksp->setup(1, ioData->ts.implicit.newton.ksp.ns.maxIts, Uref);
 ksp2->setup(1, ioData->ts.implicit.newton.ksp.ns.maxIts, Uref);
 ksp3->setup(1, ioData->ts.implicit.newton.ksp.ns.maxIts, Uref);

 tOutput->openAsciiFiles();
}

//-------------------------------------------------------------------------------
/*
template<int dim>
void ModalSolver<dim>::constructROM(VecSet<Vec<double> > &romOperator,
                      VecSet<Vec<double> > &ecVecs, VecSet<Vec<double> > &gVecs,
                      VecSet<DistSVec<double, dim> > &podVecs, int nPodVecs)  {

 com->fprintf(stderr, " ... Forming Reduced Integrator \n");
 int iVec;

 // form H ROM op
 DistSVec<double,dim> tmpVec(domain.getNodeDistInfo());
 DistSVec<double, dim> tmpRomVec(domain.getNodeDistInfo());
 ksp->printParam();
 for (iVec = 0; iVec < nPodVecs; iVec++)  {
   HOp2->apply(podVecs[iVec], tmpVec);
   ksp->solve(tmpVec, tmpRomVec);
   for (int jVec = 0; jVec < nPodVecs; jVec++)
     romOperator[iVec][jVec] = podVecs[jVec] * tmpRomVec;
 }
 com->fprintf(stderr, " ... Computed Rom Operator\n");

 // form coupling matrix ROMs
 DistSVec<double,dim> tmpECvec(domain.getNodeDistInfo());
 DistSVec<double,dim> tmpGvec(domain.getNodeDistInfo());
 for (iVec = 0; iVec < nStrMode; iVec++)  {

   tmpECvec = 0.0; 
   tmpGvec = 0.0;
   // form (A+.5*dt*H)^-1 B
   ksp->solve(DE[iVec], tmpECvec);
   ksp->solve(DX[iVec], tmpGvec);

   for (int jVec = 0; jVec < nPodVecs; jVec++)  {
     ecVecs[iVec][jVec] = podVecs[jVec] * tmpECvec;
     gVecs[iVec][jVec] = podVecs[jVec] * tmpGvec;
   }
 }
 com->fprintf(stderr, " ... Computed Coupling Rom Vectors\n");

}
*/
//-------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::constructROM2(double *romOpPlusVals, VecSet<Vec<double> > &romOperator0, double *romOpPlusVals1, double *romOpPlusVals2,
                      VecSet<Vec<double> > &ecVecs, VecSet<Vec<double> > &gVecs,
                      VecSet<DistSVec<double, dim> > &podVecs, int nPodVecs)  {
#ifdef DO_MODAL
 com->fprintf(stderr, " ... Forming Reduced System Directly \n");

 // form fluid ROM
 DistSVec<double,dim> FF(domain.getNodeDistInfo());

 // Allocate ROM operators
 VarFcn *varFcn = new VarFcn(*ioData); 
  //MatVecProdH2<dim, double, dim> *onlyHOp = new MatVecProdH2<5,double,5>(*ioData,  varFcn, tState, spaceOp, &domain);
  MatVecProdH2<dim, double, dim> *onlyHOp = new MatVecProdH2<dim,double,dim>(*ioData,  varFcn, tState, spaceOp, &domain); // PJSA 11/04/2010
 onlyHOp->evalH(0, Xref, controlVol, Uref);

 int iVec, jVec;

 // form H ROM op
 DistSVec<double,dim> tmpVec(domain.getNodeDistInfo());
 
 double romVal;
 com->barrier();
 for (iVec = 0; iVec < nPodVecs; iVec++)  {
   onlyHOp->apply(podVecs[iVec], tmpVec);
   if (ioData->ts.form == TsData::NONDESCRIPTOR)
     tmpVec /= controlVol;
   for (jVec = 0; jVec < nPodVecs; jVec++)  {
     romVal = podVecs[jVec] * tmpVec;
      romOpPlusVals[iVec*nPodVecs+jVec] = romVal;
      romOpPlusVals1[iVec*nPodVecs+jVec] = romVal;
      romOpPlusVals2[iVec*nPodVecs+jVec] = romVal;
     romOperator0[iVec][jVec] = -romVal;
   }

    romOpPlusVals[iVec*nPodVecs+iVec] += 8.0/(3.0*dt);
    romOpPlusVals1[iVec*nPodVecs+iVec] += 5.0/(3.0*dt);
    romOpPlusVals2[iVec*nPodVecs+iVec] += 3.0/(2.0*dt);

   //spaceOp has to be reset because it has been modified by the apply function
   DistSVec<double,dim> FF(domain.getNodeDistInfo());
   spaceOp->computeResidual(Xref, controlVol, Uref, FF, tState);  
 }
#ifdef DO_MODAL
// checkFluidRomStability(romOperator0, nPodVecs);
#else
 com->fprintf(stderr, "  ... ERROR: REQUIRES COMPILATION WITH ARPACK and DO_MODAL Flag\n");
 exit(-1);

#endif
 // form coupling matrix ROMs
 DistSVec<double,dim> tmpECvec(domain.getNodeDistInfo());
 DistSVec<double,dim> tmpGvec(domain.getNodeDistInfo());

 Vec<double> tmpECrom(nPodVecs);
 Vec<double> tmpGrom(nPodVecs);

 for (iVec = 0; iVec < nStrMode; iVec++)  {

   tmpECvec = DE[iVec]; 
   tmpGvec = DX[iVec];
   
   if (ioData->ts.form == TsData::NONDESCRIPTOR) {
     tmpECvec /= controlVol;
     tmpGvec /= controlVol;
   }
   for (jVec = 0; jVec < nPodVecs; jVec++)  {
     tmpECrom[jVec] = podVecs[jVec] * tmpECvec;
     tmpGrom[jVec] = podVecs[jVec] * tmpGvec;
   }

   ecVecs[iVec] =  tmpECrom;
   gVecs[iVec] = tmpGrom;

 }
 com->fprintf(stderr, " ... Computed Coupling Rom Vectors\n");

#else
 com->fprintf(stderr, "  ... ERROR: REQUIRES COMPILATION WITH ARPACK and DO_MODAL Flag\n");
 exit(-1);

#endif

}

//---------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::computeModalDisp(double sdt, Vec<double> &delWRom, double *delU, double *delY, Vec<double> &refModalF, VecSet<Vec<double> > &PtimesPhi, int nPodVecs, int timeIt) {
//For ROM
Vec<double> modalF(nStrMode);
//modalF = 0.0*ioData->ref.rv.force*refModalF;
modalF = ioData->ref.rv.force*refModalF;
for (int i = 0; i < nPodVecs; i++)
  modalF += delWRom[i]*PtimesPhi[i];


  updateModalValues(sdt, delU, delY, modalF, timeIt);
}

//---------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::computeModalDisp(double sdt, DistSVec<double, 3> &xPos, DistSVec<double, dim> &delW, double *delU, double *delY, Vec<double> &refModalF, int timeIt) {
//For FOM

 Vec<double> modalF(nStrMode);
 postOp->computeForceDerivs(xPos, Uref, delW, modalF, mX);
 //modalF += 0.0*refModalF; // DJA: DEBUG only
 modalF += refModalF;
 modalF *= ioData->ref.rv.force;

  updateModalValues(sdt, delU, delY, modalF, timeIt);
}

//---------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::updateModalValues(double sdt, double *delU, double *delY, Vec<double> &modalF, int timeIt){
//For both ROM and FOM
  double prevU, prevY;
 double *srhs = new double[nStrMode];
double sdt2 = sdt*sdt;

 if (timeIt == 0){
    for (int i = 0; i < nStrMode; i++) {
       prevU = delU[i];
       prevY = delY[i];
       srhs[i] =  modalF[i] - K[i]*prevU + (1.0/sdt - 0.5*sdt*K[i])*prevY;
       delY[i] = srhs[i]/(1.0/sdt + 0.5*sdt*K[i]);
       delU[i] = 0.5*sdt*(prevY + delY[i]) + prevU;
    }
  }else {
    for (int i = 0; i < nStrMode; i++)  {
       prevU = delU[i];
       prevY = delY[i];
       srhs[i] =  modalF[i] - K[i]*prevU + (1.0/sdt - 0.25*sdt*K[i])*prevY;
       delY[i] = srhs[i]/(1.0/sdt + 0.25*sdt*K[i]);
       delU[i] = 0.5*sdt*(prevY + delY[i]) + prevU;
     }
  }
} 
//---------------------------------------------------------------------------------

template <int dim>
void ModalSolver<dim>::constructPOD()  {

 Timer *modalTimer = domain.getTimer();
 double t0;
 
 // Initial Conditions
 VecSet<DistSVec<bcomp,dim> > prevW(nStrMode, domain.getNodeDistInfo());

 int i;
 for (i = 0; i < nStrMode; ++i)
   prevW[i] = 0.0;

 int nSteps = ioData->ts.maxIts;

 double refLength = ioData->linearizedData.refLength;

 if (ioData->linearizedData.domain == LinearizedData::FREQUENCY)  {

   // adjust (E+C) matrices for freq. domain
   double invT = 1.0/ioData->ref.rv.time;
   double refConst = ioData->ref.rv.velocity/refLength;
   for (int iMode = 0; iMode < nStrMode; iMode++)
     DE[iMode] *= invT;

   com->fprintf(stderr, " ... Adjusting E+C with vel/length: %e (%e)\n", invT, refConst);

   DistSVec<double,dim> FF(domain.getNodeDistInfo());
   DistMat<bcomp,dim> *_pcC = dynamic_cast<DistMat<bcomp,dim> *>(pcComplex);
   int iSnap = 0;

   int nSteps = ioData->ts.maxIts;
   int nSnaps = (2*nSteps+1)*nStrMode;
   double deltaFreq = ioData->linearizedData.freqStep;
   int nSnapsCoarse = 0;

   int *stepParam = new int[4];                                                               
   stepParam[0] = nSteps;
   stepParam[2] = nSnaps;
   nPadeDeriv = -1;

   com->fprintf(stderr, " ... Running POD from reduced freq. 0.0 using %d steps with refLength %f\n", nSteps, refLength);

   VecSet<DistSVec<double, dim> > snaps(nSnaps, domain.getNodeDistInfo());

   int nStepsCoarse = 0;
   double *coarseFreq = new double[11];
   int L, M, nPoints;

   if (ioData->linearizedData.padeReconst == LinearizedData::TRUE)  {
     L = ioData->linearizedData.pade.degNum;
     M = ioData->linearizedData.pade.degDen;
     nPoints = ioData->linearizedData.pade.nPoints;
     if (ioData->linearizedData.pade.freq[0] >= 0.0) coarseFreq[nStepsCoarse++] = ioData->linearizedData.pade.freq[0];
     for (i=1; i< ioData->linearizedData.pade.num; i++) 
       if (ioData->linearizedData.pade.freq[i] > 0.0) coarseFreq[nStepsCoarse++] = ioData->linearizedData.pade.freq[i];
     nPadeDeriv = int(ceil((L+M+1.0)/((double)nPoints))-1.0);
     nSnapsCoarse = nStepsCoarse*2*(nPadeDeriv+1)*nStrMode;
     if (coarseFreq[0] == 0.0)
       nSnapsCoarse -= nStrMode; // No imaginary part for the snapshot at k=0.0
     nSteps = nStepsCoarse-1;
     com->fprintf(stderr, " ... Computing the POD and their derivatives at %d coarse freq.\n", nStepsCoarse);                                 
   }
   VecSet<DistSVec<bcomp,dim> > prevDerivW(nStrMode*(nPadeDeriv+1), domain.getNodeDistInfo());
   if (nPadeDeriv+1 > 0) {
     for (i = 0; i < nStrMode*(nPadeDeriv+1); ++i)
       prevDerivW[i] = 0.0;
   }
   stepParam[1] = nStepsCoarse;
   stepParam[3] = nSnapsCoarse;
   VecSet<DistSVec<double, dim> > snapsCoarse(nSnapsCoarse, domain.getNodeDistInfo());
   double kFreq;
   for (i = 0; i <= nSteps; i++)  {

     if (ioData->linearizedData.padeReconst == LinearizedData::TRUE)
       kFreq = coarseFreq[i];
     else
       kFreq = i*deltaFreq;
     double hzFreq = kFreq*invT/6.283185307; 
     bcomp shift(0.0, kFreq);
     HOpC->evaluate(0, Xref, controlVol, Uref, FF, shift);
     com->fprintf(stderr, " ... Shifting by  k = %f(%f Hz) i\n", kFreq, hzFreq);
     if (_pcC) {
       spaceOp->computeH1(Xref, controlVol, Uref, *_pcC);
       tState->addToH1(controlVol, *_pcC, shift);
       spaceOp->applyBCsToJacobian(Uref, *_pcC);
     }
     pcComplex->setup();

     //t0 = modalTimer->getTime();
     if (ioData->linearizedData.padeReconst == LinearizedData::TRUE) {
       freqIntegrateMultipleRhs(snapsCoarse, kFreq, iSnap, prevDerivW);
     }

     else
       freqIntegrate(snaps, kFreq, iSnap, prevW);
     //modalTimer->addSnapsLinSolvTime(t0); 

   }

   if (ioData->linearizedData.padeReconst == LinearizedData::TRUE)  {
     com->fprintf(stderr, " ... Computing the Pade Reconstruction\n");
     nSteps = ioData->ts.maxIts;
     t0 = modalTimer->getTime();
     domain.padeReconstruction(snapsCoarse, snaps, stepParam, coarseFreq, deltaFreq, nStrMode, L, M, nPoints);
     modalTimer->addPadeReconstrTime(t0);
   }                                                          
   delete[] coarseFreq;
   delete[] stepParam;
   makeFreqPOD(snaps, nSnaps);
 }
 else
   createPODInTime();

}

//------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::freqIntegrate(VecSet<DistSVec<double, dim> >&snaps,
                      double kFreq, int &iSnap,
                      VecSet<DistSVec<bcomp,dim> > &prevW)  {

 Timer *modalTimer = domain.getTimer();
 double t0;

 // loop over modes
 bcomp oneReal(1.0, 0.0);
 bcomp oneImag(0.0, 1.0);
 bcomp kImag(0.0, kFreq);

 DistSVec<bcomp, dim> rhs(domain.getNodeDistInfo());
 DistSVec<bcomp, dim> delW(domain.getNodeDistInfo());
 Vec3D x0(0.0, 0.0, 0.0);

 rhs = oneReal*DX[0] + kImag*DE[0];
 t0 = modalTimer->getTime();
 if (ioData->linearizedData.padeReconst == LinearizedData::TRUE) {
   kspCompGcr->setup(1, 40, rhs);
   kspCompGcr->printParam();
   kspCompGcr->numCalcVec = 0;
 }
 else {
   kspComp->setup(1, 40, rhs);
   kspComp->printParam();
 }

 int iMode;
 for (iMode = 0; iMode < nStrMode; iMode++)  {

   // form rhs
   rhs = -oneReal*DX[iMode] - kImag*DE[iMode];

   // solve
   com->fprintf(stderr, " ... Solving for mode %d, w/rhs norm %e\n", iMode+1, rhs.norm());

   delW = prevW[iMode];
   if (ioData->linearizedData.padeReconst == LinearizedData::TRUE) {
     kspCompGcr->solve(rhs, delW);
   }
   else {
     kspComp->solve(rhs, delW);
   }
   if (ioData->problem.alltype == ProblemData::_POD_CONSTRUCTION_)  {
     snaps[iSnap++].getReal(delW);

     if (kFreq != 0.0)
       snaps[iSnap++].getImag(delW);
   }

   prevW[iMode] = delW;

 }
 modalTimer->addSnapsLinSolvTime(t0);
}

//-----------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::freqIntegrateMultipleRhs(VecSet<DistSVec<double, dim> >&snapsCoarse,
                      double kFreq, int &iSnap,
                      VecSet<DistSVec<bcomp,dim> > &prevDerivW)  {


 Timer *modalTimer = domain.getTimer();
 double t0;

 // loop over modes
 bcomp oneReal(1.0, 0.0);
 bcomp oneImag(0.0, 1.0);
 bcomp kImag(0.0, kFreq);

 int nRhsCalc = 0;
 int nRecycleVect = 10;

 DistSVec<bcomp, dim> rhs(domain.getNodeDistInfo());
 DistSVec<bcomp, dim> delW(domain.getNodeDistInfo());
 DistSVec<double, dim> delWReal(domain.getNodeDistInfo());
 DistSVec<bcomp, dim> prevWtemp(domain.getNodeDistInfo());
 Vec3D x0(0.0, 0.0, 0.0);


 rhs = oneReal*DX[0] + kImag*DE[0];

 t0 = modalTimer->getTime();
 kspCompGcr->setup(1, 40, rhs);
 kspCompGcr->printParam();


 kspCompGcr->numCalcVec = 0;
 int iMode;
 int iPadeDeriv;
 for (iMode = 0; iMode < nStrMode; iMode++)  {
   // form rhs
   kspCompGcr->numCalcVec = 0;
   rhs = -oneReal*DX[iMode] - kImag*DE[iMode];


   // solve
   com->fprintf(stderr, " ... Solving for mode %d, w/rhs norm %e\n", iMode+1, rhs.norm());

   delW = prevDerivW[iMode*(nPadeDeriv +1)];
   kspCompGcr->solveMRhs(rhs, delW);
   nRhsCalc += 1;
   if (ioData->problem.alltype == ProblemData::_POD_CONSTRUCTION_)  {
     snapsCoarse[iSnap++].getReal(delW);


     if (kFreq != 0.0)
       snapsCoarse[iSnap++].getImag(delW);
   }


   prevDerivW[iMode*(nPadeDeriv +1)] = delW;

   if (nPadeDeriv > 0 ) {

     controlVolComp = oneImag*controlVol;
     prevWtemp = prevDerivW[iMode*(nPadeDeriv +1)];
     prevWtemp *= controlVolComp;

     rhs = -oneImag*DE[iMode] - prevWtemp;

     // solve for the derivative
     com->fprintf(stderr, " ... Solving for the first derivative for mode %d\n", iMode+1);
     delW = prevDerivW[iMode*(nPadeDeriv +1)+1];
     kspCompGcr->solveMRhs(rhs, delW);
     nRhsCalc += 1;
     if (ioData->problem.alltype == ProblemData::_POD_CONSTRUCTION_)  {
       snapsCoarse[iSnap++].getReal(delW);


       snapsCoarse[iSnap++].getImag(delW);
     }

     prevDerivW[iMode*(nPadeDeriv +1)+1] = delW;

     if (nPadeDeriv > 1) {
       for (iPadeDeriv = 2; iPadeDeriv < nPadeDeriv+1; iPadeDeriv++) {

         bcomp iPadeDerivComp(0.0,iPadeDeriv);

         prevWtemp = prevDerivW[iMode*(nPadeDeriv +1)+(iPadeDeriv-1)];
         prevWtemp *= controlVolComp;


         rhs = -iPadeDerivComp*prevWtemp;

         // solve for the derivative
         if (iPadeDeriv == 2)
           com->fprintf(stderr, " ... Solving for the second derivative for mode %d\n", iMode+1);
         else if (iPadeDeriv == 3)
           com->fprintf(stderr, " ... Solving for the third derivative for mode %d\n", iMode+1);
         else
           com->fprintf(stderr, " ... Solving for the %dth derivative for mode %d\n", iPadeDeriv, iMode+1);
         delW = prevDerivW[iMode*(nPadeDeriv +1)+iPadeDeriv];
         kspCompGcr->solveMRhs(rhs, delW);
         nRhsCalc += 1;

         if (ioData->problem.alltype == ProblemData::_POD_CONSTRUCTION_)  {
           snapsCoarse[iSnap++].getReal(delW);

           snapsCoarse[iSnap++].getImag(delW);
         }

         prevDerivW[iMode*(nPadeDeriv +1)+iPadeDeriv] = delW;
       }
     }
   }

 }
 modalTimer->addSnapsLinSolvTime(t0);     

}


//-----------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::makeFreqPOD(VecSet<DistSVec<double, dim> > &snaps, int nSnaps, int nPOD, bool outputToDisk){

 Timer *modalTimer = domain.getTimer();
 DistVec<double> controlVolSqrt(domain.getNodeDistInfo());
 geoState = new DistGeoState(*ioData, &domain);
 geoState->setup1(tInput->positions, &Xref, &controlVol);


 if (nPOD == 0)
	 nPOD = (ioData->linearizedData.numPOD) ? ioData->linearizedData.numPOD : nSnaps;  //KMW

 int nconv = 0;
 if (nPOD > nSnaps) {
	 com->fprintf(stderr, " ... WARNING: nPOD > nSnaps. Changing nPOD from %d to %d\n", nPOD, nSnaps);
	 nPOD = nSnaps;
 }

 com->fprintf(stderr, " ... Computing %d POD vectors from %d snapshots\n", nPOD, nSnaps);

 if (podMethod == 0) {	// svd
#ifdef DO_SCALAPACK
  com->fprintf(stderr, "Inside DO_SCALAPACK \n");
  if (ioData->ts.form == TsData::DESCRIPTOR) {
    controlVolSqrt.pow(controlVol,0.5);
    for (int iSnap = 0; iSnap < nSnaps; ++iSnap)
      snaps[iSnap] *= controlVolSqrt;
  }

  VecSet<DistSVec<double, dim> > Utrue(nSnaps, domain.getNodeDistInfo());
  Vec<double> singVals(nSnaps);
  FullM VtrueDummy(1);	// do not need VtrueDummy

  double t0 = modalTimer->getTime();
  ParallelRom<dim> parallelRom(domain,com,domain.getNodeDistInfo());
  parallelRom.parallelSVD(snaps, Utrue, singVals.data(), VtrueDummy, nSnaps, false);

  if (ioData->ts.form == TsData::DESCRIPTOR) {
    for (int iSnap = 0; iSnap < nSnaps; iSnap++)
      Utrue[iSnap] /= controlVolSqrt;
  }
  modalTimer->addEigSolvTime(t0);  
  if (outputToDisk) {
    outputPODVectors(Utrue, singVals, nPOD);
  }
  else {	// overwrite snaps with Utrue and return
    for (int i = 0; i < nPOD; ++i) {
      snaps[i] = Utrue[i]*singVals[i];
    }
  }
#else
 com->fprintf(stderr, "*** Error: REQUIRES COMPILATION WITH SCALAPACK \n");
 exit(-1);

#endif
	}
 else if (podMethod == 1) {	// eig
#ifdef DO_MODAL

  // allocate for upper half of sym. eigprob
  double *rVals = new double[nSnaps*(nSnaps+1)/2];
  DistSVec<double, dim> CVsnap(domain.getNodeDistInfo());

  for (int i = 0; i < nSnaps; i++){
    com->fprintf(stderr," ... processing snap %d\n",i);//CBM
    if (ioData->ts.form == TsData::DESCRIPTOR) {
      CVsnap = snaps[i];
      CVsnap *= controlVol;
      for (int j = 0; j <= i; ++j) 
        rVals[(i+1)*i/2 + j] = snaps[j] * CVsnap;
     }
     else {
      for (int j = 0; j <= i; j++)
        rVals[(i+1)*i/2 + j] = snaps[j] * snaps[i];
    }
  }

  double tolerance = ioData->linearizedData.tolerance;

  com->barrier();
  ARdsSymMatrix<double> pod(nSnaps, rVals, 'U');
  com->fprintf(stderr, " ... Factoring Correlation Matrix\n");
  
  double t0 = modalTimer->getTime();
  int iSnap;
  pod.FactorA();
  ARluSymStdEig<double> podEigProb(nPOD, pod, "LM", nSnaps-1, tolerance, 300*nPOD);
  modalTimer->addCorrelMatrixTime(t0);

  t0 = modalTimer->getTime();
  nconv = podEigProb.FindEigenvectors();
  modalTimer->addEigSolvTime(t0);

	com->fprintf(stderr, " ... Got %d converged eigenvectors out of %d snaps\n", nconv, nSnaps);
	outputPODVectors(podEigProb, snaps, nPOD, nSnaps);
	delete [] rVals;

#else
 com->fprintf(stderr, "  ... ERROR: REQUIRES COMPILATION WITH ARPACK and DO_MODAL Flag\n");
 exit(-1);

#endif
 }
}

//---------------------------------------------------------------------------------------

// template<int dim>
// void ModalSolver<dim>::projectFullSoltn() {		// this function is no longer supported
// 
//  // read in POD Basis
//  int nPodVecs = ioData->rom.dimension;
//  VecSet<DistSVec<double, dim> > podVecs(nPodVecs, domain.getNodeDistInfo());
//  readPodVecs(podVecs, nPodVecs);
// 
//  // setup solvers
//  VarFcn *varFcn = new VarFcn(*ioData); 
//  //VarFcn *varFcn = new VarFcnPerfectGasEuler3D(*ioData);
//  geoState = new DistGeoState(*ioData, &domain);
//  geoState->setup1(tInput->positions, &Xref, &controlVol);
//  bcData = new DistBcDataSA<dim>(*ioData, varFcn, &domain, Xref);
// 
//  spaceOp = new SpaceOperator<dim>(*ioData, varFcn, bcData, geoState, &domain);
//  postOp = new PostOperator<dim> (*ioData, varFcn, bcData, geoState, &domain);
// 
//  // Temporal operator contains Uref for us
//  tState = new DistTimeState<dim>(*ioData, spaceOp, varFcn, &domain);
//  tState->setup(tInput->solutions,  bcData->getInletBoundaryVector(), Xref, Uref);
//  geoState->setup2(tState->getData());
// 
//  RefVal *refVal = new RefVal(ioData->ref.rv);
//  tOutput = new TsOutput<dim>(*ioData, refVal, &domain, postOp);
//  tRestart = new TsRestart(*ioData, refVal);
// 
//  tOutput->openAsciiFiles();
// 
//  // initialization
//  double* Wr = new double[nPodVecs];
//  DistSVec<double, dim> WFull(domain.getNodeDistInfo());
// 
// 
//  // read in State Vectors
//  char *snapFile = 0;
//  int sp = strlen(ioData->input.prefix) + 1;
//  snapFile = new char[sp + strlen(ioData->input.stateVecFile)];
//  sprintf(snapFile, "%s%s", ioData->input.prefix, ioData->input.stateVecFile);
// 
//  VecSet< DistSVec<double, dim> > snap(1, domain.getNodeDistInfo());
//  //double *tag = new double[1];
//  double tag = 0.0;
// 
//  // project and output forces
//  int iSnap = 0;
//  int j,ierr;
//  bool endOfFile = true;
//  while(endOfFile) {
//    endOfFile = domain.readVectorFromFile(snapFile, iSnap, &tag, snap[0]);
//    if(endOfFile) { //CBM -- BAD fix
//      for (j = 0; j < nPodVecs; ++j) {
//        Wr[j] = podVecs[j] * snap[0];
//      }
//      WFull = 0.0;
//      for (j = 0; j < nPodVecs; ++j) {
//        WFull += Wr[j] * podVecs[j];
//      }
// 
//      ierr = domain.checkSolution(varFcn, WFull);
//      if (ierr)
//        com->fprintf(stderr, " ... WARNING: %d nodes have neg. rho/P \n", ierr);
// 
//      tOutput->writeForcesToDisk(0, iSnap, 0, 0, tag, com->cpuNum(), tRestart->energy, Xref, WFull);
//      tOutput->writeLiftsToDisk((*ioData), 0, iSnap, 0, 0, tag, com->cpuNum(), tRestart->energy, Xref, WFull);
//      iSnap++;
//    }
//  }
// }

//---------------------------------------------------------------------------------------
template<int dim>
void
ModalSolver<dim>::outputModalDisp(double *delU, double *delY, double sdt, int cnt,
                                 int nStrMode, FILE *fp) {

 for (int j = 0; j < nStrMode; j++)
   com->fprintf(fp, "%e  ", delU[j]);
 com->fprintf(fp, "\n ");

}

//----------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::interpolatePOD()  {


	Timer *modalTimer = domain.getTimer();
	com->fprintf(stderr, " ... Interpolating POD on a tangent space to the Grassmann manifold \n");

	char *vecFile = tInput->podFile;
	if (!vecFile)
	{
		string str = "podFiles.in";
		vecFile    =  new char [str.size()+1];
		strcpy (vecFile, str.c_str());
		// vecFile now contains a c-string copy of str. Adam 2010.08.23 : g++4.3 was complaining
	}

	FILE *inFP = fopen(vecFile, "r");
	if (!inFP)  {
		com->fprintf(stderr, "*** Warning: No POD FILES in %s\n", vecFile);
		exit (-1);
	}
	int nData, _n;
	_n = fscanf(inFP, "%d",&nData);

	char **podFile = new char *[nData];

	for (int iData = 0; iData < nData; ++iData){
		podFile[iData] = new char[500];
		//char *podFile1 = new char[500];
		_n = fscanf(inFP, "%s", podFile[iData]);
		//podFile[iData] = podFile1;
		com->fprintf(stderr, " ... Reading POD from %s \n", podFile[iData]);
	}


	if (ioData->output.transient.podFile[0] == 0)  {
		com->fprintf(stderr, "*** ERROR: POD Basis File not specified\n");
		exit (-1);
	}
	int sp = strlen(ioData->output.transient.prefix);
	char *outFile = new char[sp + strlen(ioData->output.transient.podFile)+1];
	sprintf(outFile, "%s%s", ioData->output.transient.prefix, ioData->output.transient.podFile);

	int numPod = ioData->linearizedData.numPOD;

	double *mach = new double[nData];
	double *angle =new double[nData];
	double newMach;
	double newAngle;
	for (int iData = 0; iData < nData; ++iData)
		_n = fscanf(inFP, "%lf", mach+iData);
	_n = fscanf(inFP, "%lf", &newMach);
	for (int iData = 0; iData < nData; ++iData)
		_n = fscanf(inFP, "%lf", angle+iData);
	_n = fscanf(inFP, "%lf", &newAngle);

	com->fprintf(stderr, " ... Interpolating new POD basis at Mach = %f and angle of attack = %f from :\n",newMach, newAngle);
	for (int iData = 0; iData < nData; ++iData)
		com->fprintf(stderr,"      -> Mach = %f and angle of attack = %f\n", mach[iData], angle[iData]);


	VecSet< DistSVec<double, dim> > **pod = new VecSet< DistSVec<double, dim> >*[nData];
	for (int iData=0; iData < nData; ++iData){
		pod[iData]= new VecSet< DistSVec<double, dim> >(numPod, domain.getNodeDistInfo());
		// read in Pod Vectors
		double *eig = new double[numPod];
		domain.readVectorFromFile(podFile[iData], 0, &eig[0], (*pod[iData])[0] );
		if (numPod > eig[0])  {
			com->fprintf(stderr, "*** Warning: Resetting number of interpolated POD vectors from %d to %d\n", numPod, (int) eig[0]);
			numPod = (int) eig[0];
		}

		for (int iPod = 0; iPod < numPod; ++iPod)
			domain.readVectorFromFile(podFile[iData], iPod+1, &eig[iPod], (*pod[iData])[iPod]);
		delete [] eig;
	}

	double t0 = modalTimer->getTime();

	//build the reduced coordinates
	double *reducedMach = new double[nData];
	double *reducedAngle = new double[nData];
	if (newAngle!=0.0) {
		for (int iData=0; iData < nData; ++iData) {
			reducedMach[iData] = mach[iData]/newMach-1;
			reducedAngle[iData] = angle[iData]/newAngle-1;
		}
	}
	else {
		for (int iData=0; iData < nData; ++iData) {
			reducedMach[iData] = mach[iData]-newMach;
			reducedAngle[iData] = angle[iData]-newAngle;
		}
	}
	//find the central point
	int iDataMin = 0;
	double radMinSq = reducedMach[0]*reducedMach[0] + reducedAngle[0]*reducedAngle[0];
	double radSq = 0.0;
	for (int iData=1; iData < nData; ++iData) {
		radSq = reducedMach[iData]*reducedMach[iData] + reducedAngle[iData]*reducedAngle[iData];
		if (radSq < radMinSq){
			iDataMin = iData;
			radMinSq = radSq;
		}
	}

	//store the reference pod basis
	VecSet< DistSVec<double, dim> > podRef(*(pod[iDataMin]));

	//Logarithmic mappings
	VecSet< DistSVec<double, dim> > **projMap = new VecSet< DistSVec<double, dim> >*[nData];
	FullM matVals(numPod);

	//compute the matrices (Phi0'*Phi)^(-1)
	for (int iData=0; iData < nData; ++iData){
		if (iData!=iDataMin) {
			for (int j = 0; j < numPod; ++j) {
				for (int k = 0; k < numPod; ++k) {
					matVals[j][ k] = podRef[j] * ((*pod[iData])[k]);
				}
			}	 
			com->barrier();

			matVals.invert();

			//compute the projection mapping = Phi*(Phi0'*Phi)^(-1)-Phi0
			projMap[iData]= new VecSet< DistSVec<double, dim> >(numPod, domain.getNodeDistInfo());
			com->barrier();
			for (int iPod = 0; iPod < numPod; ++iPod) {
				(*projMap[iData])[iPod] = 0.0;
				for (int jPod = 0; jPod < numPod; ++jPod) {
					(*projMap[iData])[iPod] += (*pod[iData])[jPod] * matVals[jPod][iPod];
				}
			}
			com->barrier();
			for (int iPod = 0; iPod < numPod; ++iPod)
				(*projMap[iData])[iPod] = (*projMap[iData])[iPod] - podRef[iPod];
			com->barrier();
		}
	}

	com->barrier();
	for (int iData=0; iData<nData;++iData){
		delete pod[iData];
		com->barrier();
	}

	com->barrier();
	delete [] pod;
	//SVD decomposition of projMap



	VecSet<DistSVec<double, dim> > U(numPod, domain.getNodeDistInfo());
	double *Sigma = new double[numPod];
	FullM V(numPod);
	double *Theta = new double[numPod];
	VecSet<DistSVec<double, dim> > U2(numPod, domain.getNodeDistInfo());
	double *Sigma2 = new double[numPod];
	FullM V2(numPod);
	//Logarithmic mapping
	VecSet< DistSVec<double, dim> > **logMap = new VecSet< DistSVec<double, dim> >*[nData];
	for (int iData=0; iData < nData; ++iData) {
		if (iData!=iDataMin) {
			//compute SVD
			ParallelRom<dim> parallelRom(domain,com,domain.getNodeDistInfo()); 
			parallelRom.parallelSVD(*projMap[iData],U,Sigma,V,numPod);//call SVD here 

			com->barrier();
			delete projMap[iData];

			//compute tan^(-1) Sigma
			for (int iPod = 0; iPod < numPod; ++iPod)
				Theta[iPod] = atan(Sigma[iPod]);
			//compute the logarithmic mapping
			//build logMap
			logMap[iData]= new VecSet< DistSVec<double, dim> >(numPod, domain.getNodeDistInfo());
			for (int iPod = 0; iPod < numPod; ++iPod) {
				(*logMap[iData])[iPod] = 0.0;
				for (int jPod = 0; jPod < numPod; ++jPod){
					(*logMap[iData])[iPod] += ((Theta[jPod]*V[iPod][jPod])*U[jPod]);  // logMap = U tan^(-1) Sig V^T
				}
			}
		}
	}
	delete [] Theta;
	delete [] Sigma;

	VecSet< DistSVec<double, dim> > logMapInterp(numPod, domain.getNodeDistInfo());
	// Interpolation (Hardy's multiquadrics method is here used)
	// see "Two Dimensional Spline Interpolation Algorithms", H. Spath, A.K. Peters Ltd  
	double q = 0.25;
	double Rsq;

	//find Rsq
	double minReducedM = reducedMach[0];
	double maxReducedM = reducedMach[0];
	double minReducedA = reducedAngle[0];
	double maxReducedA = reducedAngle[0];
	for (int iData=1; iData < nData; ++iData) {
		if (reducedMach[iData] < minReducedM)
			minReducedM = reducedMach[iData];
		if (reducedMach[iData] > maxReducedM)
			maxReducedM = reducedMach[iData];
		if (reducedAngle[iData] < minReducedA)
			minReducedA = reducedAngle[iData];
		if (reducedAngle[iData] > maxReducedA)
			maxReducedA = reducedAngle[iData];
	}
	double maxDelM = maxReducedM - minReducedM;
	double maxDelA = maxReducedA - minReducedA;
	double maxConst = maxDelM > maxDelA ? maxDelM : maxDelA;
	Rsq = 0.1*maxConst;

	FullM B(nData), b(nData,1);


	double rsq;
	for (int j=0; j < nData; ++j) {
		for (int k=0; k < nData; ++k) {
			rsq = pow(reducedMach[j]-reducedMach[k],2)+pow(reducedAngle[j]-reducedAngle[k],2);
			B[j][k] = pow(rsq+Rsq,q);
		}
		rsq = pow(reducedMach[j],2)+pow(reducedAngle[j],2);
		b[j][0] = pow(rsq+Rsq,q);
	}
	B.invert();

	//compute the interpolated POD basis
	domain.hardyInterpolationLogMap(logMap,logMapInterp,nData,numPod,iDataMin,B,b);
	//Exponential mapping
	VecSet<DistSVec<double, dim> > UInterp(numPod, domain.getNodeDistInfo());
	double *SigInterp = new double[numPod];
	FullM VInterp(numPod);
	ParallelRom<dim> parallelRom(domain,com,domain.getNodeDistInfo());
	parallelRom.parallelSVD(logMapInterp,UInterp,SigInterp,VInterp,numPod);//call SVD here

	//compute podRefVInterp
	VecSet<DistSVec<double, dim> > podRefVInterp(numPod, domain.getNodeDistInfo());
	for (int iPod = 0; iPod < numPod; ++iPod) {
		podRefVInterp[iPod] = 0.0;
		for (int jPod = 0; jPod < numPod; ++jPod)
			podRefVInterp[iPod] += VInterp[jPod][iPod]*podRef[jPod];
	}

	VecSet<DistSVec<double, dim> >newPOD(numPod,domain.getNodeDistInfo());

	for (int iPod = 0; iPod < numPod; ++iPod) {
		newPOD[iPod] = cos(SigInterp[iPod])*podRefVInterp[iPod]+ sin(SigInterp[iPod])*UInterp[iPod];
        }
        double *Rmatrix =  new double[numPod*numPod];
        modifiedGramSchmidt(newPOD, Rmatrix, numPod);
		//Do a Gramm-Schmidt to finalize the POD reconstruction
		//for (int j = 0; j < iPod; ++j)
		//	newPOD[iPod] -= (newPOD[iPod] * newPOD[j])*newPOD[j];
		//newPOD[iPod] *= 1.0/newPOD[iPod].norm();

	com->fprintf(stderr, " ... Writing Interpolated POD to %s\n", outFile);
	// output interpolated pod vectors
	domain.writeVectorToFile(outFile, 0, numPod, newPOD[0]);
	double newEig=-1.0;

	for (int jPod = 0; jPod < numPod; ++jPod)  {
		com->fprintf(stderr, " ... Writing pod vec %d\n", jPod+1);
		domain.writeVectorToFile(outFile, jPod+1, newEig, newPOD[jPod]);
	}
	modalTimer->addMeshSolutionTime(t0);

	delete[] SigInterp;
	for (int iData=0; iData< nData; ++iData)
		delete [] podFile[iData];
	delete[] podFile;
	delete[] outFile;
        delete[] Rmatrix;

	com->fprintf(stderr,"End of InterpolatePOD\n'");
	modalTimer->setRunTime();


}
//----------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::evalFluidSys(VecSet<DistSVec<double, dim> > &podVecs, int nPodVecs)  {

 // Allocate ROM operators
 VarFcn *varFcn = new VarFcn(*ioData); 
 MatVecProdH2<dim, double, dim> *onlyHOp = new MatVecProdH2<dim,double,dim>(*ioData,  varFcn, tState, spaceOp, &domain);
 onlyHOp->evalH(0, Xref, controlVol, Uref);

 VecSet<Vec<double> > romOperator(nPodVecs, nPodVecs);

 // form H ROM op
 DistSVec<double,dim> tmpVec(domain.getNodeDistInfo());

 double gamma = ioData->eqs.fluidModel.gasModel.specificHeatRatio;

 double timeDimConst =  ioData->ref.mach * sqrt(gamma) / ioData->ref.length;
 com->fprintf(stderr, " ... Using M = %f, Gamma = %f, rho = %e, L = %f for dim constant: %16.10e\n",  
      ioData->ref.mach, gamma, ioData->ref.density, ioData->ref.length, timeDimConst);

 int iVec, jVec;

 for (iVec = 0; iVec < nPodVecs; iVec++)  {
   onlyHOp->apply(podVecs[iVec], tmpVec);
   if (ioData->ts.form == TsData::NONDESCRIPTOR) 
     tmpVec /= controlVol;
   for (jVec = 0; jVec < nPodVecs; jVec++)
     romOperator[iVec][jVec] = (podVecs[jVec] * tmpVec) * timeDimConst;
 }

 com->fprintf(stderr, " ... Computed Rom Operator\n");

 com->fprintf(stderr, " ... Created Fluid System\n");
 int sp = strlen(ioData->output.transient.prefix);
 char *romFile = new char[sp + strlen(ioData->output.transient.romFile)+1];
 sprintf(romFile, "%s%s", ioData->output.transient.prefix, ioData->output.transient.romFile);
 FILE *romFP = fopen(romFile, "w");
 com->barrier();

 com->fprintf(romFP, "%d 0\n", nPodVecs);
 for (iVec = 0; iVec < nPodVecs; iVec++)  {
   for (jVec = 0; jVec < nPodVecs; jVec++)
     com->fprintf(romFP, "%.16e ", romOperator[jVec][iVec]);
   com->fprintf(romFP, "\n");
 }

 delete[] romFile;

 //spaceOp has to be reset because it has been modified by the apply function
 DistSVec<double,dim> FF(domain.getNodeDistInfo());
 spaceOp->computeResidual(Xref, controlVol, Uref, FF, tState);

}

//----------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::evalAeroSys(VecSet<Vec<double> > &outRom, 
                      VecSet<DistSVec<double, dim> > &podVecs, int nPodVecs)  {

 // Allocate ROM operators
 VarFcn *varFcn = new VarFcn(*ioData);  
  MatVecProdH2<dim, double, dim> *onlyHOp = new MatVecProdH2<dim,double,dim>(*ioData,  varFcn, tState, spaceOp, &domain);
 onlyHOp->evalH(0, Xref, controlVol, Uref);

 VecSet<Vec<double> > ecVecs(nStrMode, nPodVecs);
 VecSet<Vec<double> > gVecs(nStrMode, nPodVecs);
 VecSet<Vec<double> > romOperator(nPodVecs, nPodVecs);

 // form H ROM op
 DistSVec<double,dim> tmpVec(domain.getNodeDistInfo());
 DistSVec<double,dim> tmpVec2(domain.getNodeDistInfo());

 double gamma = ioData->eqs.fluidModel.gasModel.specificHeatRatio;

 double timeDimConst = -ioData->ref.mach * sqrt(gamma) / ioData->ref.length;
 com->fprintf(stderr, " ... Using M = %f, Gamma = %f, rho = %e, L = %f for dim constant: %16.10e\n",  ioData->ref.mach, gamma, ioData->ref.density, ioData->ref.length, timeDimConst);

 int iVec, jVec;

 for (iVec = 0; iVec < nPodVecs; iVec++)  {
   onlyHOp->apply(podVecs[iVec], tmpVec);
   if (ioData->ts.form == TsData::NONDESCRIPTOR) 
     tmpVec /= controlVol;
   
   for (jVec = 0; jVec < nPodVecs; jVec++)
     romOperator[iVec][jVec] = (podVecs[jVec] * tmpVec) * timeDimConst;
 }

 com->fprintf(stderr, " ... Computed Rom Operator\n");

 // form coupling matrix ROMs
 // form Structural Component: K
 double *structSys = new double[nStrMode];

 double invDt = -1.0/ioData->ref.rv.time;

 for (iVec = 0; iVec < nStrMode; iVec++)  {
   tmpVec = DE[iVec] * invDt;
   tmpVec2 = DX[iVec];
   if (ioData->ts.form == TsData::NONDESCRIPTOR) {
     tmpVec /= controlVol;
     tmpVec2 /= controlVol;
   }

   for (jVec = 0; jVec < nPodVecs; jVec++)  {
     ecVecs[iVec][jVec] = podVecs[jVec] * tmpVec;
     gVecs[iVec][jVec] = podVecs[jVec] * tmpVec2 * timeDimConst;
   }

   structSys[iVec] = -K[iVec];
 }

 com->fprintf(stderr, " ... Computed Coupling Input Rom Vectors\n");

 int sysSize = 2*nStrMode + nPodVecs;
 double *sysVals = new double[sysSize*sysSize];
 for (int j = 0; j < sysSize*sysSize; j++)
   sysVals[j] = 0.0;
 com->fprintf(stderr, " ... Forming Aeroelastic System of size: %d\n", sysSize);


 // populate 1,1 & 2,1 block
 for (iVec = 0; iVec < nPodVecs; iVec++)  {
   for (jVec = 0; jVec < nPodVecs; jVec++)
     sysVals[iVec*sysSize+jVec] = romOperator[iVec][jVec];
   for (jVec = 0; jVec < nStrMode; jVec++)
     sysVals[iVec*sysSize + nPodVecs + jVec] = outRom[iVec][jVec];
 }

 // populate 1,2 and 3,2 block

 for (iVec = 0; iVec < nStrMode; iVec++)  {
   for (jVec = 0; jVec < nPodVecs; jVec++)
     sysVals[sysSize*nPodVecs + iVec*sysSize + jVec] = ecVecs[iVec][jVec]; 
   for (jVec = 0; jVec < nStrMode; jVec++) {
     if (jVec == iVec)
       sysVals[sysSize*nPodVecs + iVec*sysSize + nPodVecs+nStrMode + jVec] = 1.0;
     else
       sysVals[sysSize*nPodVecs + iVec*sysSize + nPodVecs+nStrMode + jVec] = 0.0;
   }
 }

 // populate 1,3 and 2,3 block
 for (iVec = 0; iVec < nStrMode; iVec++)  {
   for (jVec = 0; jVec < nPodVecs; jVec++)
     sysVals[sysSize*(nPodVecs+nStrMode) + iVec*sysSize + jVec] = gVecs[iVec][jVec];
   for (jVec = 0; jVec < nStrMode; jVec++) {
     if (jVec == iVec)
       sysVals[sysSize*(nPodVecs+nStrMode) + iVec*sysSize + nPodVecs + jVec] = structSys[iVec];
     else
       sysVals[sysSize*(nPodVecs+nStrMode) + iVec*sysSize + nPodVecs + jVec] = 0.0;
   }
 }

 com->fprintf(stderr, " ... Created Aeroelastic System\n");
 int sp = strlen(ioData->output.transient.prefix);
 char *romFile = new char[sp + strlen(ioData->output.transient.romFile)+1];
 sprintf(romFile, "%s%s", ioData->output.transient.prefix, ioData->output.transient.romFile);
 FILE *romFP = fopen(romFile, "w");
 com->barrier();

 com->fprintf(romFP, "%d %d\n", nPodVecs, nStrMode);
 for (iVec = 0; iVec < sysSize; iVec++)  {
   for (jVec = 0; jVec < sysSize; jVec++) {
     com->fprintf(romFP, "%.16e ", sysVals[jVec*sysSize+iVec]);
	 }	 
   com->fprintf(romFP, "\n");
 }

 delete[] romFile;

 //spaceOp has to be reset because it has been modified by the apply function
 DistSVec<double,dim> FF(domain.getNodeDistInfo());
 spaceOp->computeResidual(Xref, controlVol, Uref, FF, tState);

}

//------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::formOutputRom(VecSet<Vec<double> > &outputRom,
                      VecSet<DistSVec<double, dim> > &podVecs, int nPodVecs)  {

 Vec<double> modalF(nStrMode);
 modalF = 0.0;

 double gamma = ioData->eqs.fluidModel.gasModel.specificHeatRatio;
 double machSquare = ioData->ref.mach * ioData->ref.mach;
 double lengthSquare = ioData->ref.length * ioData->ref.length;

 for (int iVec = 0; iVec < nPodVecs; iVec++)  {
   postOp->computeForceDerivs(Xref, Uref, podVecs[iVec], modalF, mX);
   outputRom[iVec] = machSquare * gamma * lengthSquare * modalF;
 }
}

//------------------------------------------------------------------------------

template<int dim>
template<class Scalar>
void ModalSolver<dim>::readPodVecs(VecSet<DistSVec<Scalar, dim> > &podVecs,
                      int &nPod)  {

 // read in POD Vectors
 char *vecFile = tInput->podFile;

 double eigValue; 
 int nPodVecs;

 // read number of vecs
 DistSVec<Scalar,dim> tmpVec(domain.getNodeDistInfo());
 domain.readVectorFromFile(vecFile, 0, &eigValue, tmpVec);

 nPodVecs = (int) eigValue;
 com->fprintf(stderr, " ... There are %d total podVecs \n", nPodVecs);

 if (nPod > nPodVecs)  {
   com->fprintf(stderr, " ... There are only %d POD Vectors \n", nPodVecs);
   nPod = nPodVecs;
 }
 else
   nPodVecs = nPod;

 com->fprintf(stderr, " ... Reading %d POD Vectors from file %s\n", nPodVecs, vecFile);

 double firstEig;
 domain.readVectorFromFile(vecFile, 1, &firstEig, podVecs[0]);
 for (int iVec = 1; iVec < nPodVecs; iVec++)
   domain.readVectorFromFile(vecFile, iVec+1, &eigValue, podVecs[iVec]);

 com->fprintf(stderr, " ... Eigenvalue Ratio: (%e/%e) = %e\n", eigValue, firstEig, eigValue/firstEig);

 checkROBType(podVecs, nPod);

}

//------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::checkROBType(VecSet<DistSVec<double, dim> > &podVecs, int nPodVecs) {

  DistSVec<double, dim> temp(domain.getNodeDistInfo());  
  FullM B(nPodVecs);
  if (ioData->ts.form == TsData::DESCRIPTOR) {
    for (int i = 0; i < nPodVecs; ++i){
      temp = podVecs[i];
      temp *= controlVol;
      for (int j = 0; j < nPodVecs; ++j)
        B[j][i] = podVecs[j]*temp;
      B[i][i] -= 1.0;
    } 
  }
  else {
    for (int i = 0; i < nPodVecs; ++i){
      for (int j = 0; j < nPodVecs; ++j)
        B[j][i] = podVecs[j]*podVecs[i];
      B[i][i] -= 1.0; 
    }
  }
  double Bnorm = B.norm();

  if (ioData->ts.form == TsData::DESCRIPTOR)
    com->fprintf(stderr, " ... Norm (Phi'*A*Phi - I) = %e\n",Bnorm);
  else 
    com->fprintf(stderr, " ... Norm (Phi'*Phi - I) = %e\n",Bnorm);

  if (Bnorm > 1e-4) {
    com->fprintf(stderr, " *** ERROR: Wrong Reduced Basis Type (Descriptor vs Non-descriptor) ***\n");
    exit(-1);
  }

}
//------------------------------------------------------------------------------

template<int dim>
void ModalSolver<dim>::outputPODVectors(VecSet<DistSVec<double, dim> > &podVecs, 
                      Vec<double> &sVals, int nPOD)  {

  Timer *modalTimer = domain.getTimer();

 // write header
 int sp = strlen(ioData->output.transient.prefix) + 1;
 const char *sValExtension = ".singularVals";

 char *podFileName = new char[sp + strlen(ioData->output.transient.podFile)];
 char *sValsFileName = new char[sp + strlen(ioData->output.transient.podFile)+strlen(sValExtension)];
 sprintf(podFileName, "%s%s", ioData->output.transient.prefix, ioData->output.transient.podFile);
 if (com->cpuNum() == 0)
	 sprintf(sValsFileName, "%s%s%s", ioData->output.transient.prefix, ioData->output.transient.podFile, sValExtension);
 FILE *sValsFile;
 if (com->cpuNum() == 0)	// only open with cpu 0
   sValsFile = fopen(sValsFileName, "wt");

 com->fprintf(sValsFile,"%d\n", nPOD);
 com->fprintf(stderr, " ... Writing %d (%f)POD vectors to File\n", nPOD, (double) nPOD);

 com->fprintf(sValsFile,"Singular values\n");
 domain.writeVectorToFile(podFileName, 0, (double) nPOD, podVecs[0]);

 for (int jj = 0; jj < nPOD; ++jj)
	 com->fprintf(sValsFile,"%e ", sVals[jj]);
 com->fprintf(sValsFile,"\n");
 //computeRelativeEnergy(sValsFile, sVals, nPOD);

 //const int waitTime = 60;

 for (int jj = 0; jj < nPOD; ++jj) {
   com->fprintf(stderr, "%d %e\n", jj, sVals[jj]);
	 //com->fprintf(stderr, " ... waiting %d seconds ...\n", waitTime);
	 com->barrier();
	 //wait(waitTime);	// avoid file system crash due to excessive I/O
   domain.writeVectorToFile(podFileName, jj+1, sVals[jj]*sVals[jj], podVecs[jj]);
 }

}

#ifdef DO_MODAL
//------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::outputPODVectors(ARluSymStdEig<double> &podEigProb,
                      VecSet<DistSVec<double, dim> > &snaps, int nPOD, int numSnapsTot)  {

  Timer *modalTimer = domain.getTimer();

 DistSVec<double, dim> podVec(domain.getNodeDistInfo());
 podVec = 0.0;

 // write header
 int sp = strlen(ioData->output.transient.prefix) + 1;
 const char *sValExtension = ".singularVals";

 char *podFileName = new char[sp + strlen(ioData->output.transient.podFile)];
 char *sValsFileName = new char[sp + strlen(ioData->output.transient.podFile)+strlen(sValExtension)];
 sprintf(podFileName, "%s%s", ioData->output.transient.prefix, ioData->output.transient.podFile);
 if (com->cpuNum() == 0)
	 sprintf(sValsFileName, "%s%s%s", ioData->output.transient.prefix, ioData->output.transient.podFile, sValExtension);
 FILE *sValsFile;
 if (com->cpuNum() == 0)	// only open with cpu 0
	 sValsFile = fopen(sValsFileName, "wt");
 delete [] sValsFileName;

 com->fprintf(sValsFile,"%d\n", nPOD);
 com->fprintf(stderr, " ... Writing %d (%f) POD vectors to File\n", nPOD, (double) nPOD);
 //const char *podFileName = ioData->output.transient.podFile;

 com->fprintf(sValsFile,"Singular values\n");
 domain.writeVectorToFile(podFileName, 0, (double) nPOD, podVec);	// dummy vector

 double t0;
 double *rawEigVec = new double[numSnapsTot];
 int jj, kk;
 Vec<double> sVals(nPOD);
 //const int waitTime = 60;
 for (jj = 0; jj < nPOD; jj++)  {
   t0 = modalTimer->getTime();
   for (kk = 0; kk < numSnapsTot; kk++)
     rawEigVec[kk] = podEigProb.Eigenvector(nPOD-jj-1, kk);

   Vec<double> pVec(numSnapsTot, rawEigVec);
   podVec = 0.0;
   for (kk = 0; kk < numSnapsTot; kk++)
     podVec += snaps[kk] * pVec[kk];

   double eig = podEigProb.Eigenvalue(nPOD-jj-1);
   if (eig < 0.0) {
     com->fprintf(stderr,"Negative eigenvalue: %e\n",eig);
     break;
   }
   podVec *= 1.0/sqrt(eig);

   // do gram-schmidt on podVec
   //for (int jVec = 0; jVec < jj; jVec++)
     //podVec[jj] -= podVec[jVec] * (podVec[jj]*podVec[jVec]);

   //norm = podVec.norm();
   //podVec *= (1.0 / norm);
   modalTimer->addGramSchmidtTime(t0);
	 //com->barrier();
	 //com->fprintf(stderr, " ... waiting %d seconds ...\n", waitTime);
	 //wait(waitTime);	// avoid file system crash
	 domain.writeVectorToFile(podFileName, jj+1, eig, podVec);
	 com->fprintf(stderr, "%d %e\n", jj, sqrt(eig));
	 sVals[jj] = sqrt(eig);
 }
 delete [] rawEigVec;
 delete [] podFileName;

 for (int jj = 0; jj < nPOD; ++jj) {
	 com->fprintf(sValsFile,"%e ", sVals[jj]);
 }
 com->fprintf(sValsFile,"\n");
 computeRelativeEnergy(sValsFile, sVals, nPOD);

}
#endif

//------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::computeRelativeEnergy(FILE *sValsFile, const Vec<double> &sVals, const int nPod){

  // TODO: optimize!

  com->fprintf(sValsFile,"Relative energy: s(i)^2/sum(s(1:end).^2)\n");
  std::vector<double> relEnergy;
  int nSnap = sVals.size();

  if (totalEnergy == 0.0) {
    for (int i = 0; i < nSnap; ++i)
      totalEnergy += pow(sVals[i],2);
  }

  for (int i = 0; i < nSnap; ++i) {
    double currentRelEnergy = pow(sVals[i],2)/totalEnergy;
    com->fprintf(sValsFile,"%d %e\n", i+1, currentRelEnergy);
    relEnergy.push_back(currentRelEnergy);
  }
  com->fprintf(sValsFile,"Cumulative energy: sum(s(1:k).^2)/sum(s(1:end).^2)\n");
  double cumulativeEnergy = 0.0;
  double criteria [10] = {0.9, 0.95, 0.975, 0.99, 0.995, 0.999, 0.9995, 0.9999, 0.99995, 0.99999};
  int energyIndex [10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int handledCriteria  [10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  int critCounter = 0;
  for (int i = 0; i < nSnap; ++i) {
    cumulativeEnergy+=relEnergy[i];
    com->fprintf(sValsFile,"%d %e\n", i+1, cumulativeEnergy);

    int critCounterTmp = 0;
    for (int j = critCounter; j < 10;++j) {
      if (cumulativeEnergy >= criteria[j] && handledCriteria[j] == 0) {
        energyIndex[j] = i;
        handledCriteria[j] = 1;
        ++critCounterTmp;
      }
    }
    critCounter +=critCounterTmp;
  }
  com->fprintf(sValsFile,"Cumulative energy indices\n");

  for (int i = 0; i < 10; ++i) {
    com->fprintf(sValsFile,"%e: %d\n", criteria[i], energyIndex[i]+1);
  }
}

//------------------------------------------------------------------------------
/*  // an interesting idea, but not currently implemented (KMW)
template<int dim>
void ModalSolver<dim>::normalizeSnap(DistSVec<double, dim> &snap, const int iSnap, const int nSnaps){

	// PURPOSE: normalize snaphots

	double scalingFactor, magnitude;

	const double smallestScaling = 0.1;	// sets smallest scaling of any snaphot
	const double tau = pow((double)nSnaps-1,2)/log(smallestScaling);	

	//com->fprintf(stderr, " ... Distance weights:\n");
	magnitude = snap.norm();
	if (magnitude == 0.0)
		magnitude = 1.0;
	scalingFactor = 1.0/magnitude;
	if (ioData->snapshots.snapshotWeights == SnapshotsData::RBF) {
		double distanceWeight = exp(pow((double)iSnap,2)/tau);
		scalingFactor *= distanceWeight;
		com->fprintf(stderr, "%d: %e\n", iSnap, distanceWeight);
	}
	snap *= scalingFactor;
}
*/
//------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::wait(const int seconds )
{
	clock_t endwait;
	endwait = clock () + seconds * CLOCKS_PER_SEC ;
	while (clock() < endwait) {}
}
//------------------------------------------------------------------------------
template<int dim>
int ModalSolver<dim>::ROBInnerProductSteps(int n, int Nmax)
{

  int nPass;
  int nSteps;
  int i, j;

  if (Nmax >= n){
    nSteps = 1;
  } else{
    nPass = 1 + (int) (ceil(double(n-Nmax)/double(Nmax-1)));
    nSteps = nPass + n*(nPass-1) - Nmax*(nPass-1)*nPass/2 + (nPass-2)*(nPass-1)/2;
    if ((n - nPass*(Nmax+1) - 1) > 0)
      nSteps += n - nPass*(Nmax+1) - 1;
  }

  return nSteps;
}
//------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::ROBInnerProductSchedule(int** cache, int n, int Nmax, int nSteps )
{
  int i, j;
  int cnt = 1;
  //Initialize every column of the Cache array
  for (i = 0; i < nSteps+1; ++i)
    cache[i] = new (nothrow) int[Nmax];

  //Fill all columns with  zeros to begin
  for (i = 0; i < nSteps+1; ++i){
     for (j = 0; j < Nmax; ++j){
        cache[i][j] = 0;
     }
  }

  if (Nmax >= n){
   //Fill the second column with 1:n (only need 1 cache since it can fit everything)
   for (int j = 0; j < n; ++j)
     cache[1][j] = j+1;
   return;
  }
   //Fill the second column with 1:Nmax
   for (int j = 0; j < Nmax; ++j)
     cache[1][j] = j+1;

  //Initialize variable to be used for determine the status of elements.
  //0 -> element was already in the cache on this pass, but it is not completely exhausted
  //1 -> element currently in cache
  //2 -> element in exhausted (will never be in another cache)
  //3 -> element available
  int *SU = new int[n];

  //Initially, the cache is the first Nmax elements.  Everything else is available.
  for (i = 0; i < Nmax; ++i)
    SU[i] = 1;
  for (i = Nmax; i < n; ++i)
    SU[i] = 3;

  //Loop over each available element-> first pass
  int nAv = (int) count (SU, SU+n, 3);
  int M_Cache, M_Avail;

  for (i = 0; i < nAv; ++i){

    M_Cache = cache[cnt][Nmax-1];
    for (j = n-1; j > -1; --j){
      if (SU[j] == 3)
        break;
    }
    M_Avail = j+1;

    for (j = 0; j < Nmax-1; ++j)
      cache[cnt+1][j] = cache[cnt][j];
    cache[cnt+1][Nmax-1] = M_Avail;
    SU[M_Cache-1] = 0;
    SU[M_Avail-1] = 1;
    ++cnt;
  }

  //Update the exhausted elements and set the zero elements to 3
  for (i = 0; i < n; ++i){
     if (i < Nmax-1){
        //The first Nmax-1 elements are exhausted after 1st pass
        SU[i] = 2;
     }else if (SU[i] == 0){
        //Reset unused but available elements to available for the next pass
        SU[i] = 3;
     }
  }

/////////////////////////////////////////////////////////////////////////
  //Update the cache for the next pass.  Take the first Nmax-1 available
  //there is already 1 element in the cache we will need again.
  int inc = 0;
  for (i = 0; i < n; ++i){
     if (SU[i] == 3 && inc < Nmax-1){
        //Put the Nmax-1 smallest available elements in cache
        SU[i] = 1;
        ++inc;
     }
     //If i is in the previous cache, put it in the same position
     if (SU[i] == 1){
        for (j = 0; j < Nmax; ++j){
           if (cache[cnt][j] == i+1)
              cache[cnt+1][j] = i+1;
        }
     }
  }
  //Otherwise, put it in the first available position.
  for (i = 0; i < n; ++i){
     if (SU[i] == 1 && count(cache[cnt+1], cache[cnt+1]+Nmax, i+1) == 0){
        for (j = 0; j < Nmax; ++j){
           if (cache[cnt+1][j] == 0){
             cache[cnt+1][j] = i+1;
             break;
           }
        }
     }
  }
///////////////////////////////////////////////////////////////////////////
/*  //Update the cache for the next pass.  Take the first Nmax-1 available
  //there is already 1 element in the cache we will need again.
  int inc = 0;
  int inc2 = 0;
  for (i = 0; i < n; ++i){
     if (SU[i] == 3 && inc < Nmax-1){
        //Put the Nmax-1 smallest available elements in cache
        SU[i] = 1;
        ++inc;
     }

     if (SU[i] == 1){
        cache[cnt+1][inc2] = i+1;
        ++inc2;
     }
   }*/

   ++cnt;

  //Loop over each available element-> all other passes
  int maxLoc = 0;
  inc = 0;
  while (1){
     ++inc;
     nAv = (int) count (SU, SU+n, 3);
     for (i = 0; i < nAv; ++i){

///////////////////////////////////////////////////////////////////////////
      M_Cache = 0;
       //Determine max element in cache and its location
       for (j = 0; j < Nmax; ++j){
          if (cache[cnt][j] > M_Cache){
             M_Cache = cache[cnt][j];
             maxLoc  = j;
          }
       }
///////////////////////////////////////////////////////////////////////////
     // M_Cache = cache[cnt][Nmax-1];

      //Depending on which pass we are on, take either the maximum
      //or minimum element from the available pile
      if (inc % 2 == 0){
          //Take the maximum available element
          for (j = n-1; j > -1; --j){
            if (SU[j] == 3)
              break;
          }
       } else {
           //Take the minimum available element
           for (j = 0; j < n; ++j){
             if (SU[j] == 3)
               break;
             }
       }

       M_Avail = j+1;

///////////////////////////////////////////////////////////////////////////
       //Update cache
       for (int k = 0; k < Nmax; ++k)
         cache[cnt+1][k] = cache[cnt][k];
       cache[cnt+1][maxLoc] = M_Avail;
///////////////////////////////////////////////////////////////////////////
/*       //Update cache
       for (int k = 0; k < Nmax-1; ++k)
         cache[cnt+1][k] = cache[cnt][k];
       cache[cnt+1][Nmax-1] = M_Avail; */

       SU[M_Cache-1] = 0;
       SU[M_Avail-1] = 1;
       ++cnt;
     }
     //Update the exhausted elements and set the zero elements to 3
     for (i = 0; i < n; ++i){
        if (i < (inc+1)*Nmax - inc - 1){
           //The first Nmax-1 elements are exhausted after 1st pass
           SU[i] = 2;
        }else if (SU[i] == 0){
           //Reset unused but available elements to available for the next pass
           SU[i] = 3;
        }
     }
///////////////////////////////////////////////////////////////////////////
     //Update the cache for the next pass.  Take the first Nmax-1 available
     //there is already 1 element in the cache we will need again.
     int inc1 = 0;
     for (i = 0; i < n; ++i){
        if (SU[i] == 3 && inc1 < Nmax-1){
           //Put the Nmax-1 smallest available elements in cache
           SU[i] = 1;
           ++inc1;
        }
        //If i is in the previous cache, put it in the same position
        if (SU[i] == 1){
           for (j = 0; j < Nmax; ++j){
              if (cache[cnt][j] == i+1)
                 cache[cnt+1][j] = i+1;
           }
        }
     }
     //Otherwise, put it in the first available position.
     for (i = 0; i < n; ++i){
        if (SU[i] == 1 && count(cache[cnt+1], cache[cnt+1]+Nmax, i+1) == 0){
           for (j = 0; j < Nmax; ++j){
              if (cache[cnt+1][j] == 0){
                cache[cnt+1][j] = i+1;
                break;
              }
           }
        }
     }
///////////////////////////////////////////////////////////////////////////
/*
     //Update the cache for the next pass.  Take the first Nmax-1 available
     //there is already 1 element in the cache we will need again.
     int inc1 = 0;
     inc2 = 0;
     for (i = 0; i < n; ++i){
        if (SU[i] == 3 && inc1 < Nmax-1){
           //Put the Nmax-1 smallest available elements in cache
           SU[i] = 1;
           ++inc1;
        }

        if (SU[i] == 1){
           cache[cnt+1][inc2] = i+1;
           ++inc2;
        }
     }
*/
     ++cnt;

     if (count(SU,SU+n,2) == n || count(SU,SU+n,3) == 0)
       break;
  }
  delete [] SU;
}
//------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::ROBInnerProducts()
{

  com->fprintf(stderr, " ... Computing inner products \n");
  int numPod = ioData->linearizedData.numPOD;
  double *matVals = new double[numPod*numPod]; //will contain inner products
  double *eig = new double[numPod];

  //open POD file
  const char *vecFile = (tInput->podFile) ? tInput->podFile : "podFiles.in";
  FILE *inFP = fopen(vecFile, "r");
  if (!inFP)  {
    com->fprintf(stderr, "*** Warning: No POD FILES in %s\n", vecFile);
    exit (-1);
  }

  int nROB, nLoadMax;
  fscanf(inFP, "%d",&nROB);
  fscanf(inFP, "%d",&nLoadMax);

  if (nLoadMax>=nROB)
    nLoadMax = nROB;
 
  char **ROBFile = new char *[nROB];
  for (int iROB = 0; iROB < nROB; ++iROB) {
    ROBFile[iROB] = new char[500];
    fscanf(inFP, "%s", ROBFile[iROB]);
  }

  //inner products output file
  if (ioData->output.transient.robProductFile[0] == 0)  {
    com->fprintf(stderr, "*** ERROR: ROB Inner Products Output File not specified\n");
    exit (-1);
  }  

  int sp = strlen(ioData->output.transient.prefix);
  char *outputFile = new char[sp + strlen(ioData->output.transient.robProductFile)+1];
  sprintf(outputFile, "%s%s", ioData->output.transient.prefix, ioData->output.transient.robProductFile);
  FILE *outFP = fopen(outputFile, "w");
  if (!outFP)  {     com->fprintf(stderr, "*** Warning: No output file: %s\n", outputFile);
    exit (-1);
  }

  //R matrix from Modified Gram Schmidt output file
  FILE *outFPR;
  if (ioData->output.transient.rMatrixFile[0] != 0)  {
    int spR = strlen(ioData->output.transient.prefix);
    char *outputFileR = new char[spR + strlen(ioData->output.transient.rMatrixFile)+1];
    sprintf(outputFileR, "%s%s", ioData->output.transient.prefix, ioData->output.transient.rMatrixFile);
    outFPR = fopen(outputFileR, "w");
    if (!outFPR)  {     com->fprintf(stderr, "*** Warning: No output file: %s\n", outputFileR);
      exit (-1);
    }
  }



  //allocate memory for ROBs
  VecSet< DistSVec<double, dim> > **rob = new VecSet< DistSVec<double, dim> >*[nLoadMax];
  for (int iROB = 0; iROB < nLoadMax; ++iROB)
    rob[iROB]= new VecSet< DistSVec<double, dim> >(numPod, domain.getNodeDistInfo());

 // array to keep track of computed products
 int **computedProds = new int*[nROB];
 for (int iROB = 0; iROB < nROB; ++iROB) 
    computedProds[iROB] = new int[nROB];

 for (int iROB = 0; iROB < nROB; ++iROB) {
   for (int jROB = 0; jROB < iROB; ++jROB) {
     computedProds[iROB][jROB] = 0;
     computedProds[jROB][iROB] = 0;
   }
   computedProds[iROB][iROB] = 1;
 }
 // array to keep track of outputed R matrices 
 int *outputedRmatrix = new int[nROB];
 for (int iROB=0; iROB < nROB; ++iROB)
   outputedRmatrix[iROB] = 0;

 
 int nSteps = ROBInnerProductSteps(nROB, nLoadMax); //number of steps
 com->fprintf(stderr,"Inner Products Computation in %d Steps\n",nSteps+1);
 int **cache = new int *[nSteps+1]; 
 
 ROBInnerProductSchedule(cache, nROB, nLoadMax, nSteps);

  for (int i = 0; i < nLoadMax; ++i){
     for (int j = 0; j < nSteps+1; ++j){
        com->fprintf(stderr,"%d \t",cache[j][i]);
     }
     com->fprintf(stderr,"\n");
  }

 int iROB1, iROB2; 
  // setup solvers
 VarFcn *varFcn = new VarFcn(*ioData);
 geoState = new DistGeoState(*ioData, &domain);
 geoState->setup1(tInput->positions, &Xref, &controlVol);
 DistSVec<double, dim> temp(domain.getNodeDistInfo());
 double *Rmatrix = new double[numPod*numPod];

 for (int iStep = 0; iStep < nSteps; ++iStep) {

   // read ROBs
   for (int iData = 0; iData < nLoadMax; ++iData) {
     
     iROB1 = cache[iStep+1][iData];
     iROB2 = cache[iStep][iData];
     if (iROB1 > 0 && iROB1 != iROB2) { // need to load ROB
       domain.readVectorFromFile(ROBFile[iROB1-1], 0, &eig[0], (*rob[iData])[0] );
       if (numPod > eig[0])  {
         com->fprintf(stderr, "*** Warning: Resetting number of loaded POD vectors from %d to %d\n", numPod, (int) eig[0]);        
         numPod = (int) eig[0];
       }

       for (int iPod = 0; iPod < numPod; ++iPod)
         domain.readVectorFromFile(ROBFile[iROB1-1], iPod+1, &eig[iPod], (*rob[iData])[iPod]);
       if (ioData->linearizedData.doGramSchmidt == LinearizedData::TRUE_GS) {
         modifiedGramSchmidt(*rob[iData],Rmatrix,numPod);
         // write R matrix in output file
         if (ioData->output.transient.rMatrixFile[0] != 0)  {
           if (!outputedRmatrix[iROB1]) {
             outputedRmatrix[iROB1] = 1;
             com->fprintf(outFPR, "%d\n", iROB1);
             for (int iPod=0; iPod <numPod; ++iPod){
               for (int jPod=0; jPod <numPod; ++jPod)
                 com->fprintf(outFPR, "%.16e ", Rmatrix[iPod*numPod+jPod]);
               com->fprintf(outFPR, "\n");
             }
           }  
         } 
       } 
     }
   }
 
   for (int iData1 = 0; iData1 < nLoadMax; ++iData1) {
     for (int iData2 = 0; iData2 < iData1; ++iData2) {
       iROB1 = cache[iStep+1][iData1];
       iROB2 = cache[iStep+1][iData2];
       if (iROB1 > 0 && iROB2 > 0 && !computedProds[iROB1-1][iROB2-1]) {
        // compute inner product
         com->fprintf(stderr,"computing inner product between ROBs #%d and #%d\n",iROB1,iROB2);

         switch (ioData->ts.form) {
           case TsData::DESCRIPTOR: {
             for (int j = 0; j < numPod; j++) {
               temp = (*rob[iData2])[j];
               if (ioData->linearizedData.doGramSchmidt == LinearizedData::FALSE_GS)
                 temp *= controlVol;
               for (int k = 0; k < numPod; k++) {
                 matVals[j*numPod + k] = ((*rob[iData1])[k]) * temp;
               }
             }
             break; }
           case TsData::NONDESCRIPTOR: {
             for (int j = 0; j < numPod; j++) {
               for (int k = 0; k < numPod; k++) { 
                 matVals[j*numPod + k] = ((*rob[iData1])[k]) * ((*rob[iData2])[j]);
               }  
             }
             break; }
         }                   
         computedProds[iROB1-1][iROB2-1] = 1;
         computedProds[iROB2-1][iROB1-1] = 1;

         // write ROB inner product in output file
         com->fprintf(outFP, "%d %d\n", iROB1, iROB2);
         for (int iPod=0; iPod <numPod; ++iPod){
           for (int jPod=0; jPod <numPod; ++jPod)
             com->fprintf(outFP, "%.16e ", matVals[jPod*numPod+iPod]);
           com->fprintf(outFP, "\n");
         }
       }
     }
   }
 }


 delete [] matVals;
 delete [] eig;
 delete [] Rmatrix;
 for (int iROB = 0; iROB < nROB; ++iROB) {
   delete [] computedProds[iROB];
   delete [] ROBFile[iROB];
 }
 delete [] computedProds;
 delete [] outputedRmatrix;
 delete [] ROBFile;
 delete [] outputFile;
 for (int iStep=0; iStep < nSteps; ++iStep)
   delete [] cache[iStep];
 delete [] cache;
  
}
//------------------------------------------------------------------------------
#ifdef DO_MODAL
template<int dim>
void ModalSolver<dim>::checkFluidRomStability(VecSet<Vec<double> > &romOperator, int nPodVecs)
{

  double *rom = new double[nPodVecs*nPodVecs];
  for (int iVec = 0; iVec < nPodVecs; ++iVec) {
    for (int jVec = 0; jVec < nPodVecs; ++jVec) {
      rom[iVec*(nPodVecs)+jVec] = romOperator[iVec][jVec];
    }
  }

  // ARPACK cannot compute whole spectrum as once. hence compute first the larger magnitude eigenvalues and then the smaller magnitude ones.

  ARdsNonSymMatrix<double, double> romMat(nPodVecs, rom);
  romMat.FactorA();
  int nEv = int(floor(4*nPodVecs/5));
  // larger magnitude eigenvalues
  ARluNonSymStdEig<double> romEigProb(nEv, romMat, "LM", nPodVecs-1, 1e-8, 300*nEv);
  romEigProb.FindEigenvalues();
  // smaller magnitude eigenvalues
  ARluNonSymStdEig<double> romEigProb2(nEv, romMat, "SM", nPodVecs-1, 1e-8, 300*nEv);
  romEigProb2.FindEigenvalues();

  int stability = 1;
  for (int iPod = 0; iPod < nEv; ++iPod) {
    if (romEigProb.EigenvalueReal(iPod) > 0.0 || romEigProb2.EigenvalueReal(iPod) > 0.0){
      com->fprintf(stderr, "*** Warning: the Fluid Rom of Dimension %d has at Least One Unstable Eigenvalue\n",nPodVecs);
      stability = 0;
      break;
    }
  }
  if (stability)
    com->fprintf(stderr,"... The Fluid Rom of Dimension %d is stable\n",nPodVecs);
  delete [] rom;
}
#endif
//------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::modifiedGramSchmidt(VecSet<DistSVec<double,dim> > &vectors, double *Rmatrix, int numVecs)
{
  for (int iVec = 0; iVec < numVecs; ++iVec) {
    for (int jVec = 0; jVec < numVecs; ++jVec) 
      Rmatrix[iVec*numVecs+jVec] = 0.0;
  }
  Rmatrix[0] = vectors[0].norm();
  if (Rmatrix[0] == 0.0) {
    com->fprintf(stderr, "*** Error: Break down in Modified Gram-Schmidt: Vector #1 is zero\n");
    exit(-1);
  }
  else
    vectors[0] *= 1.0/Rmatrix[0];

  for (int jVec = 1; jVec < numVecs; ++jVec) {
    for (int iVec = 0; iVec < jVec; ++iVec) {
      Rmatrix[iVec*numVecs+jVec] = vectors[iVec] * vectors[jVec];
      vectors[iVec] -= Rmatrix[iVec*numVecs+jVec]*vectors[jVec];
    }
    Rmatrix[jVec*numVecs+jVec] = vectors[jVec].norm();
    if (Rmatrix[jVec*numVecs+jVec] == 0.0) {
      com->fprintf(stderr, "*** Error: Break down in Modified Gram-Schmidt: Vector at Step %d is zero\n",jVec);
      exit(-1);
    }
    else
      vectors[jVec] *= 1.0/Rmatrix[jVec*numVecs+jVec];
  }
}

//-------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::computeDampingRatios()
{
  int nIt = 0;
  int nMaxIt = ioData->linearizedData.maxItEV;
  int flagMaxIt[nStrMode]; 
  double sReal = 0.0;
  double sImag = 0.0;
  double absLambda = 0.0;
  double epsEV = ioData->linearizedData.epsEV;
  double dampRatio = 0.0;
  double dryModes[nStrMode];
  
  VecSet<Vec<bcomp> > compQ(nStrMode, nStrMode);
  VecSet<Vec<double> > Omega(nStrMode, nStrMode);
  VecSet<Vec<double> > eigMatA(2*nStrMode, 2*nStrMode);
  Vec<bcomp> evList(nStrMode);
  Vec<bcomp> sortEV(nStrMode);
  Vec<bcomp> sortEV_c(nStrMode);
  Vec<bcomp> sortEV_Eigen(2*nStrMode);
 
  bcomp sEVnew(0.0,0.0);
  bcomp sEVold(0.0,0.0);

  int sp = strlen(ioData->output.transient.prefix) + 1;
  char *outFile = new char[sp + strlen(ioData->output.transient.aeroelasticEigenvalues)];
  sprintf(outFile, "%s%s", ioData->output.transient.prefix, ioData->output.transient.aeroelasticEigenvalues);
  FILE *outEV = fopen(outFile, "w");
  com->barrier();
  //delete [] outFile;

  if (!outEV)  {
    com->fprintf(stderr, " *** Error: DampRatio output file not specified\n");
    exit(-1);
  }

  for (int i = 0; i < nStrMode; ++i) {
    dryModes[i] = 0.0;
    evList[i] = sEVnew;
    sortEV[i] = sEVnew;
    flagMaxIt[i] = 0; 
  }

//dry modes
  for (int i = 0; i < nStrMode; ++i) {
    for (int j = 0; j < nStrMode; ++j) {   
      if (i == j) Omega[i][j] = K[i];
      else  Omega[i][j] = 0.0;
    }
    dryModes[i] = sqrt(K[i]);
  }

  sort(dryModes, dryModes + nStrMode);
 
  for (int i = 0; i < nStrMode; ++i) {
    evList[i].imag(dryModes[i]); 
  }
 
//loop over eigenvalues/modes
  com->fprintf(stderr, "Compute eigenvalues and damping ratios ...\n");
  for (int iEV = 0; iEV < nStrMode; ++iEV) {
    sEVold = evList[iEV];
    sEVnew = evList[iEV];
    nIt = 0;
    com->fprintf(stderr, "\n ... solving for eigenmode #%i \n", iEV+1);
    int Index;

    while ((nIt == 0) || (sqrt (pow ((sEVnew.real()- sEVold.real()), 2.0)) + pow ((sEVnew.imag() - sEVold.imag()), 2.0)) / (sqrt (pow (sEVold.real(), 2.0) + pow (sEVold.imag(), 2.0))) > epsEV)  {
      
      com->fprintf(stderr, " ........ iteration #%i \n", nIt+1);
      if (nIt >= nMaxIt) {
        com->fprintf(stderr, "\n***WARNING: max. number of eigenvalue iterations (%i) is reached for mode = %i\n\n", nMaxIt, iEV+1); 
        flagMaxIt[iEV] = 1;
        break;
      }
      nIt++;
      sEVold = sEVnew;
      sReal = sEVnew.real();
      sImag = sEVnew.imag();
//call for assembly of A  
      evalMatForEvProblem(sReal, sImag, compQ, eigMatA, Omega);
//Eigen
      for (int i = 0; i < 2* nStrMode; ++i) {
          sortEV_Eigen[i] = 0.0;
      }

#ifdef USE_EIGEN3
      MatrixXd evProb(2*nStrMode,2*nStrMode);
      for (int i = 0; i < 2*nStrMode; i++) {
        for (int j = 0; j < 2*nStrMode; j++) {
          evProb(i,j) = eigMatA[j][i];
        }
      }
      Eigen::EigenSolver<MatrixXd> eigSolv(evProb);
      for (int i = 0; i < 2*nStrMode; ++i) {
        sortEV_Eigen[i] = eigSolv.eigenvalues()[i];
      }
#else
      com->fprintf(stderr, " ***  ERROR: ModalSolver<dim>::computeDampingRatios() needs EIGEN library\n");
      exit(-1); 
#endif

//sort
      for (int i = 0; i < nStrMode; ++i) {
          sortEV[i] = 0.0;
      }
      int i_pos2 = 0;
      int ipos3 [nStrMode];
      for (int i = 0; i < 2*nStrMode; ++i) {
        if ((sortEV_Eigen[i].imag() > 0.0) && (i_pos2 < nStrMode)) {
          sortEV[i_pos2] = sortEV_Eigen[i];
          ipos3[i_pos2] = i;
          i_pos2++;
        }
      }
      bcomp tmpSort(0.0,-1.0); 
      int i_addr = 0;
      for (int i = 0; i < nStrMode; i++) {
        sortEV_c[i] = sortEV[i];
      }
      for (int i = nStrMode-1; i >= 0; i--) {
        tmpSort.imag(-1.0);
        for (int j = 0; j < nStrMode; j++) {
          if (sortEV_c[j].imag() > tmpSort.imag()) {
            tmpSort = sortEV_c[j];
            i_addr = j;
          }
        }
        sortEV_c[i_addr].imag(-1.0);
        if(iEV == i)
          Index = i_addr;
        sortEV[i] = tmpSort;
      }
      sEVnew = sortEV[iEV];
      evList[iEV] = sEVnew;
      if (std::abs(sEVnew.imag()) <= 1.0e-16) {
        com->fprintf(stderr, "***WARNING: sImag is zero\n");
        break;
      }

   }
  }
  
//output damp ratio and eigenvalues
  com->fprintf(stderr, "\nWrite solution to '%s'\n\n", outFile);
 
  com->fprintf(outEV, "ModeIDnumber RealPartEigenvalue ImaginaryPartEigenvalue DampingRatio ConvergenceToSpecifiedPrecision (1 = yes, 0 = no)\n");
  
  for (int iEV = 0; iEV < nStrMode; ++iEV) {
    absLambda = sqrt (pow (evList[iEV].real(), 2.0) + pow (evList[iEV].imag(), 2.0));
    //com->fprintf(stderr, "Eigenvalue # %i ...  %e + i%e", iEV+1, evList[iEV].real(), evList[iEV].imag());

//    if (absLambda >=  1.0e-16) { 
      dampRatio = -evList[iEV].real() / absLambda;
      //com->fprintf(stderr, "\tDamping Ratio ... %e\n", dampRatio);
      //Output in File
      com->fprintf(outEV,"%i %e %e %e",iEV+1, evList[iEV].real(), evList[iEV].imag(), dampRatio);
//    } else {
//      com->fprintf(stderr, "\t***WARNING: Eigenvalue is 0!!\n");
//      com->fprintf(outEV, "%i 0.0 0.0 0.0 \tWARNING: Eigenvalue is 0", iEV+1);
//    }
    if (flagMaxIt[iEV]==1) {
//      com->fprintf(outEV, " \tWARNING: max. number of eigenvalue iterations (%i) is reached for mode = %i\n", nMaxIt, iEV+1);
      com->fprintf(outEV, " 0\n");
    } else {
      com->fprintf(outEV, " 1\n");
    }
  } 

}
//-------------------------------------------------------------------------------
#ifdef USE_EIGEN3
template<int dim>
void ModalSolver<dim>::computeREigenvector(double sReal, double sImag, int iEV,
                                           Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic> &revector)
{
  iEV--;  // subtract 1 because iEV starts with zero
  Vec<bcomp> sortEV_Eigen(2*nStrMode);
  Vec<bcomp> sortEV(nStrMode);
  Vec<bcomp> sortEV_c(nStrMode);
  
  VecSet<Vec<bcomp> > compQ(nStrMode, nStrMode);
  VecSet<Vec<double> > Omega(nStrMode, nStrMode);
  VecSet<Vec<double> > eigMatA(2*nStrMode, 2*nStrMode);

//dry modes
  for (int i = 0; i < nStrMode; ++i) {
    for (int j = 0; j < nStrMode; ++j) {   
      if (i == j) Omega[i][j] = K[i];
      else  Omega[i][j] = 0.0;
    }
  }

//loop over eigenvalues/modes
  Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> eigenvector(2*nStrMode,1);
      
//call for assembly of A  
  evalMatForEvProblem(sReal, sImag, compQ, eigMatA, Omega);

  MatrixXd evProb(2*nStrMode,2*nStrMode);
  for (int i = 0; i < 2*nStrMode; i++) {
    for (int j = 0; j < 2*nStrMode; j++) {
      evProb(i,j) = eigMatA[j][i];
    }
  }
  Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> Qmatrix(nStrMode,nStrMode);
  for (int i = 0; i < nStrMode; i++) {
    for (int j = 0; j < nStrMode; j++) {
      Qmatrix(i,j) = compQ[j][i];
    }
  }

  Eigen::EigenSolver<MatrixXd> eigSolv(evProb);

  for (int i = 0; i < 2*nStrMode; ++i) {
    sortEV_Eigen[i] = eigSolv.eigenvalues()[i];
  }


// sort
  for (int i = 0; i < nStrMode; ++i) {
    sortEV[i] = 0.0;
  }
  int i_pos2 = 0;
  int ipos3 [nStrMode]; int Index;
  for (int i = 0; i < 2*nStrMode; ++i) {
    if ((sortEV_Eigen[i].imag() > 0.0) && (i_pos2 < nStrMode)) {
      sortEV[i_pos2] = sortEV_Eigen[i];
      ipos3[i_pos2] = i;
      i_pos2++;
     }
   }
   bcomp tmpSort(0.0,-1.0); 
   int i_addr = 0;
   for (int i = 0; i < nStrMode; i++) {
     sortEV_c[i] = sortEV[i];
   }
   for (int i = nStrMode-1; i >= 0; i--) {
     tmpSort.imag(-1.0);
     for (int j = 0; j < nStrMode; j++) {
       if (sortEV_c[j].imag() > tmpSort.imag()) {
         tmpSort = sortEV_c[j];
         i_addr = j;
       }
     }
     sortEV_c[i_addr].imag(-1.0);
     if(iEV == i)
       Index = i_addr;
     sortEV[i] = tmpSort;
   }
// end sort

   eigenvector = eigSolv.eigenvectors().col(ipos3[Index]);
   for(int i=0; i<2*nStrMode; ++i) {
     revector(i,0) = eigenvector(i,0);
   }
}
#endif
//-------------------------------------------------------------------------------
#ifdef USE_EIGEN3
template<int dim>
void ModalSolver<dim>::computeOIBEI(double sReal, double sImag, int iEV)
{
  iEV--;  // subtract 1 because iEV starts with zero
  Vec<bcomp> sortEV_Eigen(2*nStrMode);
  Vec<bcomp> sortEV(nStrMode);
  Vec<bcomp> sortEV_c(nStrMode);
  
  VecSet<Vec<bcomp> > compQ(nStrMode, nStrMode);
  VecSet<Vec<double> > Omega(nStrMode, nStrMode);
  VecSet<Vec<double> > eigMatA(2*nStrMode, 2*nStrMode);
 
  bcomp sEVnew(0.0,0.0);
  bcomp sEVold(0.0,0.0);

//dry modes
  for (int i = 0; i < nStrMode; ++i) {
    for (int j = 0; j < nStrMode; ++j) {   
      if (i == j) Omega[i][j] = K[i];
      else  Omega[i][j] = 0.0;
    }
  }

//loop over eigenvalues/modes
  Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> eigenvector(2*nStrMode,1);
      
//call for assembly of A  
  evalMatForEvProblem(sReal, sImag, compQ, eigMatA, Omega);

  MatrixXd evProb(2*nStrMode,2*nStrMode);
  for (int i = 0; i < 2*nStrMode; i++) {
    for (int j = 0; j < 2*nStrMode; j++) {
      evProb(i,j) = eigMatA[j][i];
    }
  }
  Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> Qmatrix(nStrMode,nStrMode);
  for (int i = 0; i < nStrMode; i++) {
    for (int j = 0; j < nStrMode; j++) {
      Qmatrix(i,j) = compQ[j][i];
    }
  }

  Eigen::EigenSolver<MatrixXd> eigSolv(evProb);

  for (int i = 0; i < 2*nStrMode; ++i) {
    sortEV_Eigen[i] = eigSolv.eigenvalues()[i];
  }


// sort
  for (int i = 0; i < nStrMode; ++i) {
    sortEV[i] = 0.0;
  }
  int i_pos2 = 0;
  int ipos3 [nStrMode]; int Index;
  for (int i = 0; i < 2*nStrMode; ++i) {
    if ((sortEV_Eigen[i].imag() > 0.0) && (i_pos2 < nStrMode)) {
      sortEV[i_pos2] = sortEV_Eigen[i];
      ipos3[i_pos2] = i;
      i_pos2++;
     }
   }
   bcomp tmpSort(0.0,-1.0); 
   int i_addr = 0;
   for (int i = 0; i < nStrMode; i++) {
     sortEV_c[i] = sortEV[i];
   }
   for (int i = nStrMode-1; i >= 0; i--) {
     tmpSort.imag(-1.0);
     for (int j = 0; j < nStrMode; j++) {
       if (sortEV_c[j].imag() > tmpSort.imag()) {
         tmpSort = sortEV_c[j];
         i_addr = j;
       }
     }
     sortEV_c[i_addr].imag(-1.0);
     if(iEV == i)
       Index = i_addr;
     sortEV[i] = tmpSort;
   }
   sEVnew = sortEV[iEV];
// end sort
   com->fprintf(stderr, "iEV = %d\n", iEV);
   com->fprintf(stderr, "sEVnew = %e + i %e\n", sEVnew.real(), sEVnew.imag());
   com->fprintf(stderr, "sEVnew.imag() - sImag = %e\n", sEVnew.imag() - sImag);
   com->fprintf(stderr, "print sortEVs\n");
   for(int i=0; i<nStrMode; ++i) {
     com->fprintf(stderr, "%e + i %e\n", sortEV[i].real(), sortEV[i].imag());
   }
   
   Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> dRatios(nStrMode,1);
   for(int i=0; i<nStrMode; ++i) {
     dRatios(i,0) = -sortEV[i].real()/std::norm(sortEV[i]);
   }
   double minDRatio = -sEVnew.real()/sqrt(sEVnew.real()*sEVnew.real() + sEVnew.imag()*sEVnew.imag());
   double minDRatio_candidate = -sReal/sqrt(sReal*sReal+sImag*sImag);

   double OIBEI = std::fabs(minDRatio - minDRatio_candidate); 
//   double errorIndicator2 = (sqrt (pow ((sEVnew.real()- sReal), 2.0)) + pow ((sEVnew.imag() - sImag), 2.0)) / (sqrt (pow (sEVnew.real(), 2.0) + pow (sEVnew.imag(), 2.0)));
   com->fprintf(stderr, "OIBEI = %e\n", OIBEI);   
   printErrorIndicatorOutput(OIBEI,"OIBEI");
 
// find imaginary part of an eigenvalue that is closest to sImag
/*  double cur_dist;  double dist=10000000;  
  for (int i = 0; i < 2*nStrMode; ++i) {
    cur_dist = abs(sortEV_Eigen[i].imag() - sImag);
    if (cur_dist < dist) { 
      dist = cur_dist;
      Index = i;
    }
  } 
  com->fprintf(stderr, "closest distance = %e\n", dist);
  com->fprintf(stderr, "sortEV_Eigen[%d] = %e + i %e\n", Index, sortEV_Eigen[Index].real(), sortEV_Eigen[Index].imag());
*/

}
#endif
//-------------------------------------------------------------------------------
#ifdef USE_EIGEN3
template<int dim>
void ModalSolver<dim>::computeEigenvectorsAndResidual(double sReal, double sImag, int iEV, 
                                                      Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic> &revector,
                                                      Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic> &levector,
                                                      Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic> &residual)
{
  iEV--;  // subtract 1 because iEV starts with zero
  Vec<bcomp> sortEV_Eigen(2*nStrMode);
  Vec<bcomp> sortEV(nStrMode);
  Vec<bcomp> sortEV_c(nStrMode);
  
  VecSet<Vec<bcomp> > compQ(nStrMode, nStrMode);
  VecSet<Vec<double> > Omega(nStrMode, nStrMode);
  VecSet<Vec<double> > eigMatA(2*nStrMode, 2*nStrMode);
 
  bcomp sEVnew(0.0,0.0);
  bcomp sEVold(0.0,0.0);

//dry modes
  for (int i = 0; i < nStrMode; ++i) {
    for (int j = 0; j < nStrMode; ++j) {   
      if (i == j) Omega[i][j] = K[i];
      else  Omega[i][j] = 0.0;
    }
  }

//loop over eigenvalues/modes
  Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> eigenvector(2*nStrMode,1);
      
//call for assembly of A  
  evalMatForEvProblem(sReal, sImag, compQ, eigMatA, Omega);

  MatrixXd evProb(2*nStrMode,2*nStrMode);
  for (int i = 0; i < 2*nStrMode; i++) {
    for (int j = 0; j < 2*nStrMode; j++) {
      evProb(i,j) = eigMatA[j][i];
    }
  }
  Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> Qmatrix(nStrMode,nStrMode);
  for (int i = 0; i < nStrMode; i++) {
    for (int j = 0; j < nStrMode; j++) {
      Qmatrix(i,j) = compQ[j][i];
    }
  }

  Eigen::EigenSolver<MatrixXd> eigSolv(evProb);

  for (int i = 0; i < 2*nStrMode; ++i) {
    sortEV_Eigen[i] = eigSolv.eigenvalues()[i];
  }


// sort
  for (int i = 0; i < nStrMode; ++i) {
    sortEV[i] = 0.0;
  }
  int i_pos2 = 0;
  int ipos3 [nStrMode]; int Index;
  for (int i = 0; i < 2*nStrMode; ++i) {
    if ((sortEV_Eigen[i].imag() > 0.0) && (i_pos2 < nStrMode)) {
      sortEV[i_pos2] = sortEV_Eigen[i];
      ipos3[i_pos2] = i;
      i_pos2++;
     }
   }
   bcomp tmpSort(0.0,-1.0); 
   int i_addr = 0;
   for (int i = 0; i < nStrMode; i++) {
     sortEV_c[i] = sortEV[i];
   }
   for (int i = nStrMode-1; i >= 0; i--) {
     tmpSort.imag(-1.0);
     for (int j = 0; j < nStrMode; j++) {
       if (sortEV_c[j].imag() > tmpSort.imag()) {
         tmpSort = sortEV_c[j];
         i_addr = j;
       }
     }
     sortEV_c[i_addr].imag(-1.0);
     if(iEV == i)
       Index = i_addr;
     sortEV[i] = tmpSort;
   }
   sEVnew = sortEV[iEV];
// end sort
   com->fprintf(stderr, "iEV = %d\n", iEV);
   com->fprintf(stderr, "sEVnew = %e + i %e\n", sEVnew.real(), sEVnew.imag());
   com->fprintf(stderr, "sEVnew.imag() - sImag = %e\n", sEVnew.imag() - sImag);
   com->fprintf(stderr, "print sortEVs\n");
   for(int i=0; i<nStrMode; ++i) {
     com->fprintf(stderr, "%e + i %e\n", sortEV[i].real(), sortEV[i].imag());
   }
   
   Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> dRatios(nStrMode,1);
   for(int i=0; i<nStrMode; ++i) {
     dRatios(i,0) = -sortEV[i].real()/std::norm(sortEV[i]);
   }
   double minDRatio = -sEVnew.real()/sqrt(sEVnew.real()*sEVnew.real() + sEVnew.imag()*sEVnew.imag());
   double minDRatio_candidate = -sReal/sqrt(sReal*sReal+sImag*sImag);

   double errorIndicator2 = std::fabs(minDRatio - minDRatio_candidate); 
//   double errorIndicator2 = (sqrt (pow ((sEVnew.real()- sReal), 2.0)) + pow ((sEVnew.imag() - sImag), 2.0)) / (sqrt (pow (sEVnew.real(), 2.0) + pow (sEVnew.imag(), 2.0)));
   com->fprintf(stderr, "error indicator 2 = %e\n", errorIndicator2);   

   const char* output = "errorIndicator2";
   ofstream out(output, ios::out);
   if(!out) {
     cerr << "Error: cannot open file" << output << endl;
     exit(-1);
   } 
   out << errorIndicator2 << endl;
   out.close();
 
// find imaginary part of an eigenvalue that is closest to sImag
/*  double cur_dist;  double dist=10000000;  
  for (int i = 0; i < 2*nStrMode; ++i) {
    cur_dist = std::abs(sortEV_Eigen[i].imag() - sImag);
    if (cur_dist < dist) { 
      dist = cur_dist;
      Index = i;
    }
  } 
  com->fprintf(stderr, "closest distance = %e\n", dist);
  com->fprintf(stderr, "sortEV_Eigen[%d] = %e + i %e\n", Index, sortEV_Eigen[Index].real(), sortEV_Eigen[Index].imag());
*/

  com->fprintf(stderr,"nStrMode is %d\n",nStrMode);
  for(int i=0; i<2*nStrMode; ++i) {
    com->fprintf(stderr, " ... revector(%d,0) = (%e,%e)\n",i,revector(i,0).real(),revector(i,0).imag());
  }
/*
  if(isFirstErrorIndicator) {
    eigenvector = eigSolv.eigenvectors().col(ipos3[Index]);
    com->fprintf(stderr, "Eigenvalue ... %e + i %e\n", sReal, sImag);
    double normeigenvector = eigenvector.norm();
    com->fprintf(stderr, "size of eigenvector is %d\n", eigenvector.size());
    com->fprintf(stderr, "size of revector is %d\n", revector.size());
//    revector -= eigenvector;
//    com->fprintf(stderr, "relative difference: %e\n", revector.norm()/normeigenvector);
    for(int i=0; i<2*nStrMode; ++i) {
      revector(i,0) = eigenvector(i,0);
    }

    for (int i = 0; i < 2*nStrMode; i++) {
      for (int j = 0; j < 2*nStrMode; j++) {
        evProb(j,i) = eigMatA[j][i];
      }
    }
    Eigen::EigenSolver<MatrixXd> leigSolv(evProb);

    for (int i = 0; i < 2*nStrMode; ++i) {
      sortEV_Eigen[i] = leigSolv.eigenvalues()[i];
    }
 
// find imaginary part of an eigenvalue that is closest to sImag
//  dist=10000000;  
    double cur_dist;  double dist=10000000;  
    for (int i = 0; i < 2*nStrMode; ++i) {
      cur_dist = std::abs(sortEV_Eigen[i].imag() - sEVnew.imag());
      if (cur_dist < dist) { 
        dist = cur_dist;
        Index = i;
      }
    }
    eigenvector = leigSolv.eigenvectors().col(Index);
    for(int i=0; i<2*nStrMode; ++i) levector(i,0) = eigenvector(i,0);

// compute residual
    complex<double> lambda(sReal,sImag);
    Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic> reduced_revector(nStrMode,1);

    for (int jj = 0; jj < nStrMode; jj++) {
      residual(jj,0) = (K[jj] + lambda*lambda)*revector(jj,0);
      reduced_revector(jj,0) = revector(jj,0);
    }

    residual += Qmatrix*reduced_revector;

  }
*/
}
//------------------------------------------------------------------------------
#ifdef USE_EIGEN3
template<int dim>
void ModalSolver<dim>::computeNonlinearEigenResidualNormalized(double sReal, double sImag, 
                                                     Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic> rEigenVector,
                                                     double normalizationTerm, 
                                                     Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic> &residual,
                                                     const char* output) 
{

  complex<double> lambda(sReal,sImag);
  Vec<bcomp> GAMrEigenVector(nStrMode);

  for (int jj = 0; jj < nStrMode; jj++) {
    residual(jj,0) = (K[jj] + lambda*lambda)*rEigenVector(jj,0);
  }

  multiplyGAM(sReal, sImag, rEigenVector.data(), GAMrEigenVector);
  for(int jj=0; jj<nStrMode; ++jj) {
    residual(jj,0) += GAMrEigenVector[jj];
  }
  double RBEI = residual.norm()/normalizationTerm;
  printErrorIndicatorOutput(RBEI, output);
  com->fprintf(stderr, "RBEI = %e\n", RBEI);   

}
#endif
//------------------------------------------------------------------------------
#ifdef USE_EIGEN3
template<int dim>
void ModalSolver<dim>::computeNonlinearEigenResidual(double sReal, double sImag, 
                                                     Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic> rEigenVector, 
                                                     Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic> &residual,
                                                     const char* output) 
{

  complex<double> lambda(sReal,sImag);
  Vec<bcomp> GAMrEigenVector(nStrMode);

  for (int jj = 0; jj < nStrMode; jj++) {
    residual(jj,0) = (K[jj] + lambda*lambda)*rEigenVector(jj,0);
  }

  multiplyGAM(sReal, sImag, rEigenVector.data(), GAMrEigenVector);
  for(int jj=0; jj<nStrMode; ++jj) {
    residual(jj,0) += GAMrEigenVector[jj];
  }
  double RBEI = residual.norm();
  printErrorIndicatorOutput(RBEI, output);
  com->fprintf(stderr, "RBEI = %e\n", RBEI);   

}
#endif
//------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::printErrorIndicatorOutput(double errorIndicator, const char* output)
{
  ofstream out(output, ios::out);
  if(!out) {
    cerr << "Error: cannot open file" << output << endl;
    exit(-1);
  } 
  out << errorIndicator << endl;
  out.close();
} 
//------------------------------------------------------------------------------
template<int dim>
double ModalSolver<dim>::computeResidualDenominator(double sReal, double sImag, 
                                                   Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic> rEigenVector, 
                                                   Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic> lEigenVector) 
{
  complex<double> lambda(sReal,sImag);
  Vec<bcomp> GAMrEigenVector(nStrMode);
  Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic> residual(nStrMode,1);

  for (int jj = 0; jj < nStrMode; jj++) {
    residual(jj,0) = 2.0*lambda*rEigenVector(jj,0);
  }

  multiply_dQdLambda(sReal, sImag, rEigenVector.data(), GAMrEigenVector);
  Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic> reduced_lEigenVector(nStrMode,1);
  for(int jj=0; jj<nStrMode; ++jj) {
    reduced_lEigenVector(jj,0) = lEigenVector(jj,0);
    residual(jj,0) += GAMrEigenVector[jj];
  }
  return (reduced_lEigenVector.adjoint()*residual).norm();

}
#endif 
//------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::evalMatForEvProblem(double sReal, double sImag, VecSet<Vec<bcomp> > &compGAM,
                      VecSet<Vec<double> > &eigMatA, VecSet<Vec<double> > &Omega)
{
  VecSet<Vec<double> > realQ(nStrMode, nStrMode);
  VecSet<Vec<double> > imagQ(nStrMode, nStrMode);
  VecSet<Vec<double> > Aaa(nStrMode, nStrMode);
  VecSet<Vec<double> > Aab(nStrMode, nStrMode);
  VecSet<Vec<double> > Aba(nStrMode, nStrMode);
  VecSet<Vec<double> > Abb(nStrMode, nStrMode);
  computeGAM(sReal, sImag, compGAM);


// assemble Matrix for EV problem
  for (int iMode = 0; iMode < nStrMode; iMode++)  {
    for (int jj = 0; jj < nStrMode; jj++)  {

      realQ[iMode][jj] = compGAM[iMode][jj].real();
      imagQ[iMode][jj] = compGAM[iMode][jj].imag();

      Aaa[iMode][jj] = 0.0;
     
      if (iMode == jj)  {
        Aab[iMode][jj] = 1.0;
      }
      else {
        Aab[iMode][jj] = 0.0;
      }
    }
    Aba[iMode] = (-1.0)*Omega[iMode] - realQ[iMode] + (sReal/sImag)*imagQ[iMode];

    Abb[iMode] = (-1.0/sImag)*imagQ[iMode];
  }

  for (int iMode = 0; iMode < 2*nStrMode; iMode++)  {
    for (int jj = 0; jj < 2*nStrMode; jj++)  {
      if (iMode < nStrMode && jj < nStrMode) {
        eigMatA[iMode][jj] = Aaa[iMode][jj];
      }
      else if (iMode >= nStrMode && jj < nStrMode){
        eigMatA[iMode][jj] = Aab[iMode-nStrMode][jj];
      }
      else if (iMode < nStrMode && jj >= nStrMode){
        eigMatA[iMode][jj] = Aba[iMode][jj-nStrMode];
      }
      else {
        eigMatA[iMode][jj] = Abb[iMode-nStrMode][jj-nStrMode];
      }
    }
  }

}
//------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::computeGAM(double sReal, double sImag, VecSet<Vec<bcomp> > &compGAM)
{
// get aero operators and form GAM
  bcomp oneReal(1.0, 0.0);
  bcomp oneImag(0.0, 1.0);
  bcomp sCompl(sReal, sImag);
  Vec<double> modalFr(nStrMode);
  Vec<double> modalFi(nStrMode);
  modalFr = 0.0;
  modalFi = 0.0;
 
  DistSVec<bcomp, dim> rhs(domain.getNodeDistInfo());
  DistSVec<bcomp, dim> delW(domain.getNodeDistInfo());
  DistSVec<double, dim> FF(domain.getNodeDistInfo());
  DistSVec<double, dim> delWreal(domain.getNodeDistInfo());
  DistSVec<double, dim> delWimag(domain.getNodeDistInfo());
  DistMat<bcomp,dim> *_pcC = dynamic_cast<DistMat<bcomp,dim> *>(pcComplex);
  
  double invT = 1.0/ioData->ref.rv.time;
  
  Timer *modalTimer = domain.getTimer();
  double t0;

  rhs = sCompl*DE[0] + oneReal*DX[0];
  if (ioData->linearizedData.padeReconst == LinearizedData::TRUE) {
    kspCompGcr->setup(1, 40, rhs);
    kspCompGcr->numCalcVec = 0;
  }
  else {
    kspComp->setup(1, 40, rhs);
  }
//form [s/invT*A+H]
  sCompl /= invT; 
  HOpC->evaluate(0, Xref, controlVol, Uref, FF, sCompl);
  if (_pcC) {
    spaceOp->computeH1(Xref, controlVol, Uref, *_pcC);
    tState->addToH1(controlVol, *_pcC, sCompl);
    spaceOp->applyBCsToJacobian(Uref, *_pcC);
  }
  pcComplex->setup();
//loop over modes
  sCompl *= invT;
  for (int iMode = 0; iMode < nStrMode; iMode++)  {
 
    // form [s(E+C)+G]
    rhs = 0.0;
    rhs = sCompl*DE[iMode] + oneReal*DX[iMode];
   
    // solve for [sA+H]^(-1) * [s(E+C)+G]
    delW = 0.0;
    t0 = modalTimer->getTime();

    if (ioData->linearizedData.padeReconst == LinearizedData::TRUE) {
      kspCompGcr->solve(rhs, delW);
    }
    else {
      kspComp->solve(rhs, delW);
    }
    modalTimer->addKspTime(t0);
 
    delWreal.getReal(delW);
    delWimag.getImag(delW);
    
  // multiply by P
    postOp->computeForceDerivs(Xref, Uref, delWreal, modalFr, mX);
    postOp->computeForceDerivs(Xref, Uref, delWimag, modalFi, mX);
    
  // build GAM Q
    compGAM[iMode] = (ioData->ref.rv.force/ioData->ref.length)*oneReal*modalFr + (ioData->ref.rv.force/ioData->ref.length)*oneImag*modalFi;
  }
 
}

//------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::multiplyGAM(double sReal, double sImag, complex<double> *reigvector, Vec<bcomp> &GAMreigenvector)
{
// get aero operators and form GAM
  bcomp oneReal(1.0, 0.0);
  bcomp oneImag(0.0, 1.0);
  bcomp sCompl(sReal, sImag);
  Vec<double> modalFr(nStrMode);
  Vec<double> modalFi(nStrMode);
  modalFr = 0.0;
  modalFi = 0.0;
 
  DistSVec<bcomp, dim> rhs(domain.getNodeDistInfo());
  DistSVec<bcomp, dim> delW(domain.getNodeDistInfo());
  DistSVec<double, dim> FF(domain.getNodeDistInfo());
  DistSVec<double, dim> delWreal(domain.getNodeDistInfo());
  DistSVec<double, dim> delWimag(domain.getNodeDistInfo());
  DistMat<bcomp,dim> *_pcC = dynamic_cast<DistMat<bcomp,dim> *>(pcComplex);
  
  double invT = 1.0/ioData->ref.rv.time;
  
  Timer *modalTimer = domain.getTimer();
  double t0;

  rhs = sCompl*DE[0] + oneReal*DX[0];
  if (ioData->linearizedData.padeReconst == LinearizedData::TRUE) {
    kspCompGcr->setup(1, 40, rhs);
    kspCompGcr->numCalcVec = 0;
  }
  else {
    kspComp->setup(1, 40, rhs);
  }
//form [s/invT*A+H]
  sCompl /= invT; 
  HOpC->evaluate(0, Xref, controlVol, Uref, FF, sCompl);
  if (_pcC) {
    spaceOp->computeH1(Xref, controlVol, Uref, *_pcC);
    tState->addToH1(controlVol, *_pcC, sCompl);
    spaceOp->applyBCsToJacobian(Uref, *_pcC);
  }
  pcComplex->setup();
//loop over modes
  sCompl *= invT;
   
  rhs = 0.0;
  for (int iMode = 0; iMode < nStrMode; iMode++)  {
    rhs += reigvector[iMode]*(sCompl*DE[iMode] + oneReal*DX[iMode]);
  }
   
  // solve for [sA+H]^(-1) * [s(E+C)+G]
  delW = 0.0;
  t0 = modalTimer->getTime();

  if (ioData->linearizedData.padeReconst == LinearizedData::TRUE) {
    kspCompGcr->solve(rhs, delW);
  }
  else {
    kspComp->solve(rhs, delW);
  }
  modalTimer->addKspTime(t0);
 
  delWreal.getReal(delW);
  delWimag.getImag(delW);
    
  // multiply by P
  postOp->computeForceDerivs(Xref, Uref, delWreal, modalFr, mX);
  postOp->computeForceDerivs(Xref, Uref, delWimag, modalFi, mX);
    
  // build GAM Q
  GAMreigenvector = (ioData->ref.rv.force/ioData->ref.length)*oneReal*modalFr + (ioData->ref.rv.force/ioData->ref.length)*oneImag*modalFi;
}

//------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::multiply_dQdLambda(double sReal, double sImag, complex<double> *reigvector, Vec<bcomp> &GAMreigenvector)
{
// get aero operators and form GAM
  bcomp oneReal(1.0, 0.0);
  bcomp oneImag(0.0, 1.0);
  bcomp sCompl(sReal, sImag);
  Vec<double> modalFr(nStrMode);
  Vec<double> modalFi(nStrMode);
  modalFr = 0.0;
  modalFi = 0.0;
 
  DistSVec<bcomp, dim> rhs(domain.getNodeDistInfo());
  DistSVec<bcomp, dim> delW(domain.getNodeDistInfo());
  DistSVec<bcomp, dim> delW2(domain.getNodeDistInfo());
  DistSVec<double, dim> FF(domain.getNodeDistInfo());
  DistSVec<double, dim> delWreal(domain.getNodeDistInfo());
  DistSVec<double, dim> delWimag(domain.getNodeDistInfo());
  DistMat<bcomp,dim> *_pcC = dynamic_cast<DistMat<bcomp,dim> *>(pcComplex);
  
  double invT = 1.0/ioData->ref.rv.time;
  
  Timer *modalTimer = domain.getTimer();
  double t0;

  //com->fprintf(stderr, "Output dry modes********************************\n");
  //for (int i = 0; i < nStrMode; ++i) {
  //  com->fprintf(stderr, "%e\t", K[i]);
  //}
  //com->fprintf(stderr, "\nEnd Output dry modes****************************\n");

  rhs = 0.0;
  for (int iMode = 0; iMode < nStrMode; iMode++)  {
    rhs += -reigvector[iMode]*(sCompl*DE[iMode] + oneReal*DX[iMode]);
  }
   
  if (ioData->linearizedData.padeReconst == LinearizedData::TRUE) {
    kspCompGcr->setup(1, 40, rhs);
//    kspCompGcr->printParam();
    kspCompGcr->numCalcVec = 0;
  }
  else {
    kspComp->setup(1, 40, rhs);
//    kspComp->printParam();
  }
//form [s/invT*A+H]
  sCompl /= invT; 
  HOpC->evaluate(0, Xref, controlVol, Uref, FF, sCompl);
  if (_pcC) {
    spaceOp->computeH1(Xref, controlVol, Uref, *_pcC);
    tState->addToH1(controlVol, *_pcC, sCompl);
    spaceOp->applyBCsToJacobian(Uref, *_pcC);
  }
  pcComplex->setup();
//loop over modes
  sCompl *= invT;
   
  // solve for [sA+H]^(-1) * [s(E+C)+G] * reigvector
  delW = 0.0;
  t0 = modalTimer->getTime();

  if (ioData->linearizedData.padeReconst == LinearizedData::TRUE) {
    kspCompGcr->solve(rhs, delW);
  }
  else {
    kspComp->solve(rhs, delW);
  }
  modalTimer->addKspTime(t0);

  for (int iMode = 0; iMode < nStrMode; iMode++)  {
    delW += reigvector[iMode]*DE[iMode];
  }

  if (ioData->linearizedData.padeReconst == LinearizedData::TRUE) {
    kspCompGcr->setup(1, 40, delW);
//    kspCompGcr->printParam();
    kspCompGcr->numCalcVec = 0;
  }
  else {
    kspComp->setup(1, 40, delW);
//    kspComp->printParam();
  }
  delW2 = 0.0;
  t0 = modalTimer->getTime();

  if (ioData->linearizedData.padeReconst == LinearizedData::TRUE) {
    kspCompGcr->solve(delW, delW2);
  }
  else {
    kspComp->solve(delW, delW2);
  }
  modalTimer->addKspTime(t0);
 
  delWreal.getReal(delW2);
  delWimag.getImag(delW2);
    
  // multiply by P
  postOp->computeForceDerivs(Xref, Uref, delWreal, modalFr, mX);
  postOp->computeForceDerivs(Xref, Uref, delWimag, modalFi, mX);
    
  // build GAM Q
  GAMreigenvector = (ioData->ref.rv.force/ioData->ref.length)*oneReal*modalFr + (ioData->ref.rv.force/ioData->ref.length)*oneImag*modalFi;
}

//------------------------------------------------------------------------------
template<int dim>
void ModalSolver<dim>::computeGenAeroForceMat()
{
  int sp = strlen(ioData->output.transient.prefix) + 1;
   
  char *outFileGAM = new char[sp + strlen(ioData->output.transient.gamData)];
  sprintf(outFileGAM, "%s%s", ioData->output.transient.prefix, ioData->output.transient.gamData);
  FILE *outGAM = fopen(outFileGAM, "w");
  com->barrier();
  delete [] outFileGAM;

  char *outFileGAMF = new char[sp + strlen(ioData->output.transient.gamFData)];
  sprintf(outFileGAMF, "%s%s", ioData->output.transient.prefix, ioData->output.transient.gamFData);
  FILE *outGAMF = fopen(outFileGAMF, "w");
  com->barrier();
  delete [] outFileGAMF;
  

  for (int numF = 0; numF < ioData->linearizedData.numFreq; numF++) {
    if (ioData->linearizedData.gamFreq[numF] >= 0.0) {
      double GAMfreq = ioData->linearizedData.gamFreq[numF];
      com->fprintf(stderr, " ... GAMfreq%i = %e\n", numF, GAMfreq);

      VecSet<Vec<bcomp> > compQ(nStrMode, nStrMode);
      VecSet<Vec<bcomp> > compQstar(nStrMode, nStrMode);

      double invT = 1.0/ioData->ref.rv.time;
      GAMfreq *= invT;

      computeGAM(0.0, GAMfreq, compQ);

      for (int i = 0; i < nStrMode; ++i) {
        for (int j = 0; j < nStrMode; ++j) {
          compQ[i][j] *= -1.0;
          compQstar[i][j] = compQ[i][j]*(2.0/((pow(ioData->ref.rv.velocity,2))*ioData->ref.rv.density));
        }
      }

      if (strlen(ioData->output.transient.gamData) > 0) {
        if (!outGAM)  {
          com->fprintf(stderr, " *** Error: GAMData output file not specified\n");
          exit(-1);
        }

        com->fprintf(stderr, " ... write  Qstar (generalized aerodynamic matrix) for GAMFrequency%i = %e to file\n", numF+1, ioData->linearizedData.gamFreq[numF]);//gen aero matrix
        com->fprintf(outGAM, "\n%i %i %i 4QHH 1P,5E16.9\n", nStrMode, nStrMode, 1);//header

        for (int i = 0; i < nStrMode; ++i) {
          com->fprintf(outGAM, "%i %i\n", i+1, 2*nStrMode);
          for (int j = 0; j < nStrMode; ++j) {
            com->fprintf(outGAM, "%e %e ", compQstar[i][j].real(), compQstar[i][j].imag());
          }
          com->fprintf(outGAM, "\n");
        }
        com->fprintf(outGAM, "%i %i %i\n%e\n", nStrMode+1, 1, 1, ioData->linearizedData.gamFreq[numF]);
      }

      if (strlen(ioData->output.transient.gamFData) > 0) {
        if (!outGAMF)  {
          com->fprintf(stderr, " *** Error: GAMFData output file not specified\n");
          exit(-1);
        }

       com->fprintf(stderr, " ... write Q (generalized aerodynamic force matrix) for GAMFrequency%i = %e to file\n", numF+1, ioData->linearizedData.gamFreq[numF]);//gen aero FORCE matrix
       com->fprintf(outGAMF, "\n%i %i %i 4QHH 1P,5E16.9\n", nStrMode, nStrMode, 1);//header

       for (int i = 0; i < nStrMode; ++i) {
         com->fprintf(outGAMF, "%i %i\n", i+1, 2*nStrMode);
         for (int j = 0; j < nStrMode; ++j) {
           com->fprintf(outGAMF, "%e %e ", compQ[i][j].real(), compQ[i][j].imag());
         }
         com->fprintf(outGAMF, "\n");
       }
       com->fprintf(outGAMF, "%i %i %i\n%e\n", nStrMode+1, 1, 1, ioData->linearizedData.gamFreq[numF]);
       
   
      }
    }
  }

}
//------------------------------------------------------------------------------
