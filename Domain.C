#include <BcFcn.h>
#include <FluxFcn.h>
#include <RecFcnDesc.h>
#include <DistNodalGrad.h>
#include <DistEdgeGrad.h>
#include <DistExtrapolation.h>
#include <DistMacroCell.h>
#include <DistBcData.h>
#include <DistExactRiemannSolver.h>
#include <SubDomain.h>
#include <DistGeoState.h>
#include <DistVector.h>
#include <DistMatrix.h>
#include <Communicator.h>
#include <ErrorHandler.h>
#include <PostFcn.h>
#include <LowMachPrec.h>
#include <GeoState.h>
#include <NodalGrad.h>
#include <FluidSelector.h>
#include <MatVecProd.h>

#include <cstdio>
#include <cmath>
#include <unistd.h>

#include <ExactSolution.h>

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeTimeStep(double cfl, double dualtimecfl, double viscous, FemEquationTerm *fet, VarFcn *varFcn, 
			     DistGeoState &geoState,
			     DistSVec<double,3> &X, DistVec<double> &ctrlVol, DistSVec<double,dim> &V,
			     DistVec<double> &dt, DistVec<double> &idti, DistVec<double> &idtv, DistVec<double> &dtau,
									  DistVec<double> &irey, TimeLowMachPrec &tprec, SpatialLowMachPrec &sprec, 
									  DistLevelSetStructure *distLSS)
{

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
	  if(distLSS)
		  subDomain[iSub]->computeTimeStep(fet, varFcn, geoState(iSub), X(iSub), V(iSub), 
													  dt(iSub), idti(iSub), idtv(iSub), dtau(iSub), tprec, &((*distLSS)(iSub)) );
	  else
		  subDomain[iSub]->computeTimeStep(fet, varFcn, geoState(iSub), X(iSub), V(iSub), 
													  dt(iSub), idti(iSub), idtv(iSub), dtau(iSub), tprec );

    subDomain[iSub]->sndData(*volPat, reinterpret_cast<double (*)[1]>(idti.subData(iSub)));
  }

  volPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*volPat, reinterpret_cast<double (*)[1]>(idti.subData(iSub)));

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->sndData(*volPat, reinterpret_cast<double (*)[1]>(idtv.subData(iSub)));

  volPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*volPat, reinterpret_cast<double (*)[1]>(idtv.subData(iSub)));

// Included (MB)
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    double (*idtimev) = idtv.subData(iSub);
    double (*idtimei) = idti.subData(iSub);
    double (*dtime) = dt.subData(iSub);
    double (*dualtime) = dtau.subData(iSub);
    double (*ireynolds) = irey.subData(iSub);
    double (*volume) = ctrlVol.subData(iSub);

    for (int i = 0; i < ctrlVol.subSize(iSub); ++i) 
	 {
/*
		 if( distLSS && !((*distLSS)(iSub).isActive(0.0,i)))
		 {
			 dtime[i] = cfl*volume[i];
			 dualtime[i] = dualtimecfl *volume[i];
			 ireynolds[i] = -sprec.getViscousRatio();
		 }
		 else
		 {
*/
      //   idtimev[i] = idtimev[i] / volume[i];
      dtime[i] = cfl *volume[i]/(-1.0*idtimei[i] + viscous*idtimev[i]);
      dualtime[i] = dualtimecfl *volume[i]/(-1.0*idtimei[i] + viscous*idtimev[i]);
      ireynolds[i] = abs(-sprec.getViscousRatio()*idtimev[i] / idtimei[i]);
/*
    }
*/
  }
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Domain::computeDerivativeOfInvReynolds(FemEquationTerm *fet, VarFcn *varFcn, DistGeoState &geoState,
			     DistSVec<double,3> &X, DistSVec<double,3> &dX, DistVec<double> &ctrlVol,
			     DistVec<double> &dCtrlVol, DistSVec<double,dim> &V, DistSVec<double,dim> &dV,
			     DistVec<double> &idti, DistVec<double> &dIdti, DistVec<double> &idtv, DistVec<double> &dIdtv,
			     DistVec<double> &dIrey, double dMach, TimeLowMachPrec&tprec, SpatialLowMachPrec &sprec)
{

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeDerivativeOfTimeStep(fet, varFcn, geoState(iSub), X(iSub), dX(iSub), V(iSub), dV(iSub), dIdti(iSub), dIdtv(iSub), dMach, tprec);
    subDomain[iSub]->sndData(*volPat, reinterpret_cast<double (*)[1]>(dIdti.subData(iSub)));
  }

  volPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*volPat, reinterpret_cast<double (*)[1]>(dIdti.subData(iSub)));

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->sndData(*volPat, reinterpret_cast<double (*)[1]>(dIdtv.subData(iSub)));

  volPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*volPat, reinterpret_cast<double (*)[1]>(dIdtv.subData(iSub)));

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    double (*idtimev) = idtv.subData(iSub);
    double (*dIdtimev) = dIdtv.subData(iSub);
    double (*idtimei) = idti.subData(iSub);
    double (*dIdtimei) = dIdti.subData(iSub);
    double (*dIreynolds) = dIrey.subData(iSub);
    double (*volume) = ctrlVol.subData(iSub);
    double (*dVolume) = dCtrlVol.subData(iSub);
    for (int i = 0; i < ctrlVol.subSize(iSub); ++i) {
      dIdtimev[i] = (dIdtimev[i]*volume[i] - (idtimev[i]*volume[i])*dVolume[i]) / (volume[i]*volume[i]);
      dIreynolds[i] = -sprec.getViscousRatio()*(dIdtimev[i]*idtimei[i] - idtimev[i]*dIdtimei[i]) / (idtimei[i]*idtimei[i]);
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeTimeStep(double cfl, double dualtimecfl, double viscous, FemEquationTerm *fet, VarFcn *varFcn, DistGeoState &geoState,
                             DistVec<double> &ctrlVol, DistSVec<double,dim> &V,
                             DistVec<double> &dt, DistVec<double> &idti, DistVec<double> &idtv, DistVec<double> &dtau,
			     TimeLowMachPrec &tprec, DistVec<int> &fluidId, DistVec<double>* umax)
{

  int iSub;

  if (umax)
    *umax = 0.0;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    if (umax)
      subDomain[iSub]->computeTimeStep(fet, varFcn, geoState(iSub), V(iSub), dt(iSub), idti(iSub), idtv(iSub), dtau(iSub), tprec, fluidId(iSub),&(*umax)(iSub));
    else
      subDomain[iSub]->computeTimeStep(fet, varFcn, geoState(iSub), V(iSub), dt(iSub), idti(iSub), idtv(iSub), dtau(iSub), tprec, fluidId(iSub));
    subDomain[iSub]->sndData(*volPat, reinterpret_cast<double (*)[1]>(dt.subData(iSub)));
  }

  volPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*volPat, reinterpret_cast<double (*)[1]>(dt.subData(iSub)));

  if (umax) {
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {  
      subDomain[iSub]->sndData(*volPat, reinterpret_cast<double (*)[1]>(umax->subData(iSub)));
    }

    volPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*volPat, reinterpret_cast<double (*)[1]>(umax->subData(iSub)));

  }

  dtau = -dualtimecfl * ctrlVol / dt;
  dt = -cfl * ctrlVol / dt;

  if (umax)
    *umax = -0.333f * ctrlVol / (*umax);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
void Domain::computeGradientsLeastSquares(DistSVec<double,3> &X,
					  DistSVec<double,6> &R,
					  DistSVec<Scalar,dim> &var,
					  DistSVec<Scalar,dim> &ddx,
					  DistSVec<Scalar,dim> &ddy,
					  DistSVec<Scalar,dim> &ddz)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeGradientsLeastSquares(X(iSub), R(iSub), var(iSub),
						  ddx(iSub), ddy(iSub), ddz(iSub));

  CommPattern<Scalar> *vPat = getCommPat(var);
  assemble(vPat, ddx);
  assemble(vPat, ddy);
  assemble(vPat, ddz);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar>
void Domain::computeDerivativeOfGradientsLeastSquares(DistSVec<double,3> &X, DistSVec<double,3> &dX,
					  DistSVec<double,6> &R, DistSVec<double,6> &dR,
					  DistSVec<Scalar,dim> &var, DistSVec<Scalar,dim> &dvar, DistSVec<Scalar,dim> &dddx,
					  DistSVec<Scalar,dim> &dddy, DistSVec<Scalar,dim> &dddz)
{

  double t0 = timer->getTime();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    SVec<Scalar,dim> dummy(dddx(iSub));
    subDomain[iSub]->computeDerivativeOfGradientsLeastSquares(X(iSub), dX(iSub), R(iSub), dR(iSub), var(iSub), dvar(iSub),
                                                              dddx(iSub), dddy(iSub), dddz(iSub));
  }


  CommPattern<Scalar> *vPat = getCommPat(var);
  assemble(vPat, dddx);
  assemble(vPat, dddy);
  assemble(vPat, dddz);

  timer->addNodalGradTime(t0);
}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim, class Scalar>
void Domain::computeDerivativeOfGradientsLeastSquares(dRdXoperators<dim> &dRdXop, DistSVec<double,3> &dX,
                                                      DistSVec<double,6> &dR, DistSVec<double,dim> &dV,
                                                      DistSVec<Scalar,dim> &dddx,
                                                      DistSVec<Scalar,dim> &dddy, DistSVec<Scalar,dim> &dddz)
{

  double t0 = timer->getTime();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    SVec<Scalar,dim> dummy(dddx(iSub));
    subDomain[iSub]->computeDerivativeOfGradientsLeastSquares(dRdXop.dddxdX[iSub], dRdXop.dddydX[iSub], dRdXop.dddzdX[iSub],
                                                              dRdXop.dddxdR[iSub], dRdXop.dddydR[iSub], dRdXop.dddzdR[iSub],
                                                              dRdXop.dddxdV[iSub], dRdXop.dddydV[iSub], dRdXop.dddzdV[iSub],
                                                              dX(iSub), dR(iSub), dV(iSub), dddx(iSub), dddy(iSub), dddz(iSub));
  }

  timer->addNodalGradTime(t0);

  CommPattern<Scalar> *vPat = getCommPat(dV);
  assemble(vPat, dddx);
  assemble(vPat, dddy);
  assemble(vPat, dddz);

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim, class Scalar>
void Domain::computeTransposeDerivativeOfGradientsLeastSquares(dRdXoperators<dim> &dRdXop, 
                                                               DistSVec<Scalar,dim> &dddx,
                                                               DistSVec<Scalar,dim> &dddy, 
                                                               DistSVec<Scalar,dim> &dddz,
                                                               DistSVec<double,3> &dX2,
                                                               DistSVec<double,6> &dR2,
                                                               DistSVec<double,dim> &dV2)
{

  double t0 = timer->getTime();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeTransposeDerivativeOfGradientsLeastSquares(dRdXop.dddxdX[iSub], dRdXop.dddydX[iSub], dRdXop.dddzdX[iSub],
                                                                       dRdXop.dddxdR[iSub], dRdXop.dddydR[iSub], dRdXop.dddzdR[iSub],
                                                                       dRdXop.dddxdV[iSub], dRdXop.dddydV[iSub], dRdXop.dddzdV[iSub],
                                                                       dddx(iSub), dddy(iSub), dddz(iSub), dX2(iSub), dR2(iSub), dV2(iSub));
  }

  timer->addNodalGradTime(t0);

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim, class Scalar>
void Domain::computeDerivativeOperatorsOfGradientsLeastSquares(DistSVec<double,3> &X, DistSVec<double,6> &R, DistSVec<Scalar,dim> &var, 
                                                               dRdXoperators<dim> &dRdXop)
{

  double t0 = timer->getTime();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeDerivativeOperatorsOfGradientsLeastSquares(X(iSub), R(iSub), var(iSub), 
                                  *dRdXop.dddxdX[iSub], *dRdXop.dddydX[iSub], *dRdXop.dddzdX[iSub], 
                                  *dRdXop.dddxdR[iSub], *dRdXop.dddydR[iSub], *dRdXop.dddzdR[iSub],
                                  *dRdXop.dddxdV[iSub], *dRdXop.dddydV[iSub], *dRdXop.dddzdV[iSub]);

  timer->addNodalGradTime(t0);

}

//------------------------------------------------------------------------------
// least square gradient involving only nodes of same fluid (multiphase flow and FSI)
// d2d$ 
template<int dim, class Scalar>
void Domain::computeGradientsLeastSquares(DistSVec<double,3> &X,
                                          DistVec<int> &fluidId,
                                          DistSVec<double,6> &R,
                                          DistSVec<Scalar,dim> &var,
                                          DistSVec<Scalar,dim> &ddx,
                                          DistSVec<Scalar,dim> &ddy,
                                          DistSVec<Scalar,dim> &ddz,
                                          bool linFSI, DistLevelSetStructure *distLSS,
                                          bool includeSweptNodes)
{

  if(distLSS) {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->computeGradientsLeastSquares(X(iSub), fluidId(iSub), R(iSub), var(iSub),
                                                    ddx(iSub), ddy(iSub), ddz(iSub), linFSI, &((*distLSS)(iSub)),includeSweptNodes);
  } else {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->computeGradientsLeastSquares(X(iSub), fluidId(iSub), R(iSub), var(iSub),
                                                    ddx(iSub), ddy(iSub), ddz(iSub), linFSI, 0,
                                                    includeSweptNodes);
  }

  CommPattern<Scalar> *vPat = getCommPat(var);
  assemble(vPat, ddx);
  assemble(vPat, ddy);
  assemble(vPat, ddz);

}

//------------------------------------------------------------------------------
// least square gradient of single variable involving only nodes of same fluid (multiphase flow and FSI)
template<class Scalar>
void Domain::computeGradientLeastSquares(DistSVec<double,3> &X,
                                         DistVec<int> &fluidId,
                                         DistSVec<double,6> &R,
                                         DistVec<Scalar> &var,
                                         DistVec<Scalar> &ddx,
                                         DistVec<Scalar> &ddy,
                                         DistVec<Scalar> &ddz,
                                         DistLevelSetStructure *distLSS)
{

  if(distLSS) {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->computeGradientLeastSquares(X(iSub), fluidId(iSub), R(iSub), var(iSub),
                                                    ddx(iSub), ddy(iSub), ddz(iSub), &((*distLSS)(iSub)));
  } else {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->computeGradientLeastSquares(X(iSub), fluidId(iSub), R(iSub), var(iSub),
                                                    ddx(iSub), ddy(iSub), ddz(iSub), 0);
  }

  CommPattern<Scalar> *vPat = getCommPat(var);
  assemble(vPat, ddx);
  assemble(vPat, ddy);
  assemble(vPat, ddz);

}

//------------------------------------------------------------------------------
// least square gradient involving only nodes of fluid (FSI)
// Wstar is involved in gradient computation
template<int dim, class Scalar>
void Domain::computeGradientsLeastSquares(DistSVec<double,3> &X,
                                          DistVec<int> &fluidId,
                                          DistSVec<double,6> &R,
                                          DistSVec<Scalar,dim> &var,
					  DistSVec<Scalar,dim> &Wstarij,
					  DistSVec<Scalar,dim> &Wstarji,
					  DistVec<int> &countWstarij, DistVec<int> &countWstarji,
                                          DistSVec<Scalar,dim> &ddx, DistSVec<Scalar,dim> &ddy,
                                          DistSVec<Scalar,dim> &ddz,
					  bool linFSI, DistLevelSetStructure *distLSS)
{

  if(distLSS) {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->computeGradientsLeastSquares(X(iSub), fluidId(iSub), R(iSub), var(iSub),
						    Wstarij(iSub), Wstarji(iSub), 
						    countWstarij(iSub), countWstarji(iSub),
                                                    ddx(iSub), ddy(iSub), ddz(iSub), linFSI, &((*distLSS)(iSub)));
  } else {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->computeGradientsLeastSquares(X(iSub), fluidId(iSub), R(iSub), var(iSub),
                                                    ddx(iSub), ddy(iSub), ddz(iSub), linFSI, 0);
  }

  CommPattern<Scalar> *vPat = getCommPat(var);
  assemble(vPat, ddx);
  assemble(vPat, ddy);
  assemble(vPat, ddz);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
void Domain::computeGradientsGalerkin(DistVec<double> &ctrlVol, DistSVec<double,3> &wii,
				      DistSVec<double,3> &wij, DistSVec<double,3> &wji,
				      DistSVec<Scalar,dim> &var, DistSVec<Scalar,dim> &ddx,
				      DistSVec<Scalar,dim> &ddy, DistSVec<Scalar,dim> &ddz)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeGradientsGalerkin(ctrlVol(iSub), wii(iSub), wij(iSub), wji(iSub),
					      var(iSub), ddx(iSub), ddy(iSub), ddz(iSub));

  CommPattern<Scalar> *vPat = getCommPat(var);
  assemble(vPat, ddx);
  assemble(vPat, ddy);
  assemble(vPat, ddz);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar>
void Domain::computeDerivativeOfGradientsGalerkin(DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol,
                      DistSVec<double,3> &wii, DistSVec<double,3> &wij,
                      DistSVec<double,3> &wji, DistSVec<double,3> &dwii,
                      DistSVec<double,3> &dwij, DistSVec<double,3> &dwji,
                      DistSVec<Scalar,dim> &var, DistSVec<Scalar,dim> &dvar, DistSVec<Scalar,dim> &dddx,
                      DistSVec<Scalar,dim> &dddy, DistSVec<Scalar,dim> &dddz)
{

  DistSVec<Scalar,dim> ddx(getNodeDistInfo());
  DistSVec<Scalar,dim> ddy(getNodeDistInfo());
  DistSVec<Scalar,dim> ddz(getNodeDistInfo());

  double t0 = timer->getTime();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeDerivativeOfGradientsGalerkin(ctrlVol(iSub), dCtrlVol(iSub), wii(iSub), wij(iSub), wji(iSub),
                                       dwii(iSub), dwij(iSub), dwji(iSub), var(iSub), dvar(iSub), ddx(iSub), ddy(iSub), ddz(iSub), dddx(iSub), dddy(iSub), dddz(iSub));

  timer->addNodalGradTime(t0);

  CommPattern<Scalar> *vPat = getCommPat(var);
  assemble(vPat, dddx);
  assemble(vPat, dddy);
  assemble(vPat, dddz);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
void Domain::computeGradientsGalerkinT(DistVec<double> &ctrlVol,
                DistSVec<double,3> &wii, DistSVec<double,3> &wij,
                DistSVec<double,3> &wji, DistSVec<Scalar,dim> &var,
                DistSVec<Scalar,dim> &var1, DistSVec<Scalar,dim> &var2,
                DistSVec<Scalar,dim> &ddx, DistSVec<Scalar,dim> &ddy,
                DistSVec<Scalar,dim> &ddz)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeGradientsGalerkinT(ctrlVol(iSub), wii(iSub),
                wij(iSub), wji(iSub), var(iSub), var1(iSub), var2(iSub),
                ddx(iSub), ddy(iSub), ddz(iSub));


  CommPattern<Scalar> *vPat = getCommPat(var);
  assemble(vPat, ddx);
  assemble(vPat, ddy);
  assemble(vPat, ddz);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeMultiDimLimiter(RecFcnLtdMultiDim<dim> *recFcn, DistSVec<double,3> &X,
				    DistVec<double> &ctrlVol, DistSVec<double,dim> &V,
												DistSVec<double,dim> &dVdx, 
												DistSVec<double,dim> &dVdy,
												DistSVec<double,dim> &dVdz, 
												DistSVec<double,dim> &Vmin, DistSVec<double,dim> &Vmax, 
												DistSVec<double,dim> &phi)
{

  int iSub;

#pragma omp parallel for
	for(iSub = 0; iSub < numLocSub; ++iSub) 
	{
    subDomain[iSub]->computeMinMaxStencilValues(V(iSub), Vmin(iSub), Vmax(iSub));
    subDomain[iSub]->sndData(*vecPat, Vmin.subData(iSub));
  }

  vecPat->exchange();

#pragma omp parallel for
  for(iSub = 0; iSub < numLocSub; ++iSub) 
  {
    subDomain[iSub]->minRcvData(*vecPat, Vmin.subData(iSub));
    subDomain[iSub]->sndData(*vecPat, Vmax.subData(iSub));
  }

  vecPat->exchange();

#pragma omp parallel for
  for(iSub = 0; iSub < numLocSub; ++iSub) 
  {	  
    subDomain[iSub]->maxRcvData(*vecPat, Vmax.subData(iSub));
    subDomain[iSub]->computeMultiDimLimiter(recFcn, X(iSub), ctrlVol(iSub),
					    V(iSub), dVdx(iSub), dVdy(iSub), dVdz(iSub),
					    Vmin(iSub), Vmax(iSub), phi(iSub));
    subDomain[iSub]->sndData(*vecPat, phi.subData(iSub));
  }

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->minRcvData(*vecPat, phi.subData(iSub));

    double (*locphi)[dim] = phi.subData(iSub);
    double (*locdVdx)[dim] = dVdx.subData(iSub);
    double (*locdVdy)[dim] = dVdy.subData(iSub);
    double (*locdVdz)[dim] = dVdz.subData(iSub);

    for (int i=0; i<phi.subSize(iSub); ++i) {
      for (int k=0; k<dim; ++k) {
	locdVdx[i][k] *= locphi[i][k];
	locdVdy[i][k] *= locphi[i][k];
	locdVdz[i][k] *= locphi[i][k];
      }
    }
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Domain::computeDerivativeOfMultiDimLimiter(RecFcnLtdMultiDim<dim> *recFcn, DistSVec<double,3> &X, DistSVec<double,3> &dX,
				    DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol, DistSVec<double,dim> &V, DistSVec<double,dim> &dV,
				    DistSVec<double,dim> &dVdx, DistSVec<double,dim> &dVdy, DistSVec<double,dim> &dVdz,
				    DistSVec<double,dim> &ddVdx, DistSVec<double,dim> &ddVdy, DistSVec<double,dim> &ddVdz,
                    DistSVec<double,dim> &Vmin, DistSVec<double,dim> &dVmin, DistSVec<double,dim> &Vmax,
                    DistSVec<double,dim> &dVmax, DistSVec<double,dim> &phi, DistSVec<double,dim> &dphi)
{

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeDerivativeOfMinMaxStencilValues(V(iSub), dV(iSub), Vmin(iSub), dVmin(iSub), Vmax(iSub), dVmax(iSub));
    subDomain[iSub]->sndData(*vecPat, Vmin.subData(iSub));
  }

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->minRcvData(*vecPat, Vmin.subData(iSub));
    subDomain[iSub]->sndData(*vecPat, Vmax.subData(iSub));
  }

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->maxRcvData(*vecPat, Vmax.subData(iSub));
    subDomain[iSub]->sndData(*vecPat, dVmin.subData(iSub));
  }

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->minRcvData(*vecPat, dVmin.subData(iSub));
    subDomain[iSub]->sndData(*vecPat, dVmax.subData(iSub));
  }

  vecPat->exchange();


#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->maxRcvData(*vecPat, dVmax.subData(iSub));
    subDomain[iSub]->computeDerivativeOfMultiDimLimiter(recFcn, X(iSub), dX(iSub), ctrlVol(iSub), dCtrlVol(iSub),
					    V(iSub), dV(iSub), dVdx(iSub), dVdy(iSub), dVdz(iSub), ddVdx(iSub), ddVdy(iSub), ddVdz(iSub),
					    Vmin(iSub), dVmin(iSub), Vmax(iSub), dVmax(iSub), phi(iSub), dphi(iSub));
    subDomain[iSub]->sndData(*vecPat, phi.subData(iSub));
  }

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->minRcvData(*vecPat, phi.subData(iSub));
    subDomain[iSub]->sndData(*vecPat, dphi.subData(iSub));
  }

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->minRcvData(*vecPat, dphi.subData(iSub));

    double (*locphi)[dim] = phi.subData(iSub);
    double (*locdVdx)[dim] = dVdx.subData(iSub);
    double (*locdVdy)[dim] = dVdy.subData(iSub);
    double (*locdVdz)[dim] = dVdz.subData(iSub);

    double (*dlocphi)[dim] = dphi.subData(iSub);
    double (*dlocdVdx)[dim] = ddVdx.subData(iSub);
    double (*dlocdVdy)[dim] = ddVdy.subData(iSub);
    double (*dlocdVdz)[dim] = ddVdz.subData(iSub);

    for (int i=0; i<dphi.subSize(iSub); ++i) {
      for (int k=0; k<dim; ++k) {
	dlocdVdx[i][k] *= locphi[i][k];
	dlocdVdy[i][k] *= locphi[i][k];
	dlocdVdz[i][k] *= locphi[i][k];
	dlocdVdx[i][k] += locdVdx[i][k] * dlocphi[i][k];
	dlocdVdy[i][k] += locdVdy[i][k] * dlocphi[i][k];
	dlocdVdz[i][k] += locdVdz[i][k] * dlocphi[i][k];
      }
    }
  }

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    double (*locphi)[dim] = phi.subData(iSub);
    double (*locdVdx)[dim] = dVdx.subData(iSub);
    double (*locdVdy)[dim] = dVdy.subData(iSub);
    double (*locdVdz)[dim] = dVdz.subData(iSub);

    for (int i=0; i<phi.subSize(iSub); ++i) {
      for (int k=0; k<dim; ++k) {
	locdVdx[i][k] *= locphi[i][k];
	locdVdy[i][k] *= locphi[i][k];
	locdVdz[i][k] *= locphi[i][k];
      }
    }
  }

}

//------------------------------------------------------------------------------
template<int dim>
void Domain::computeMultiDimLimiter(DistSVec<double,3> &X,
												DistVec<double> &A, DistSVec<double,dim> &V,
												DistSVec<double,dim> &dVdx, 
												DistSVec<double,dim> &dVdy,
												DistSVec<double,dim> &dVdz,
												DistLevelSetStructure *distLSS)
{

	DistSVec<double,dim> *Vmin = new DistSVec<double,dim>(getNodeDistInfo());
	DistSVec<double,dim> *Vmax = new DistSVec<double,dim>(getNodeDistInfo());
	DistSVec<double,dim> *phi  = new DistSVec<double,dim>(getNodeDistInfo());

	int iSub;

#pragma omp parallel for
	for(iSub = 0; iSub < numLocSub; ++iSub) 
	{
		subDomain[iSub]->computeMinMaxStencilValues(V(iSub), (*Vmin)(iSub), (*Vmax)(iSub), &((*distLSS)(iSub)));

		subDomain[iSub]->sndData(*vecPat, (*Vmin).subData(iSub));
	}

 	vecPat->exchange();

#pragma omp parallel for
 	for(iSub = 0; iSub < numLocSub; ++iSub) 
	{
 		subDomain[iSub]->minRcvData(*vecPat, (*Vmin).subData(iSub));
		
 		subDomain[iSub]->sndData(*vecPat, (*Vmax).subData(iSub));
 	}
	
 	vecPat->exchange();

	RecFcnLtdMultiDim<dim> *rf = new RecFcnVenkat<dim>(1.0, 0.1);

#pragma omp parallel for
 	for(iSub = 0; iSub < numLocSub; ++iSub) 
	{	  
 		subDomain[iSub]->maxRcvData(*vecPat, (*Vmax).subData(iSub));

 		subDomain[iSub]->computeMultiDimLimiter(rf, X(iSub), A(iSub), V(iSub), 
 															 dVdx(iSub), dVdy(iSub), dVdz(iSub),
 															 (*Vmin)(iSub), (*Vmax)(iSub), (*phi)(iSub), &((*distLSS)(iSub)) );

 		subDomain[iSub]->sndData(*vecPat, (*phi).subData(iSub));
 	}

 	vecPat->exchange();

#pragma omp parallel for
  	for(iSub = 0; iSub < numLocSub; ++iSub) 
  	{
 		subDomain[iSub]->minRcvData(*vecPat, (*phi).subData(iSub));

 		double (*locphi)[dim]  = (*phi).subData(iSub);
 		double (*locdVdx)[dim] = dVdx.subData(iSub);
 		double (*locdVdy)[dim] = dVdy.subData(iSub);
 		double (*locdVdz)[dim] = dVdz.subData(iSub);

 		for(int i=0; i<(*phi).subSize(iSub); ++i) 
 		{			
 			for(int k=0; k<dim; ++k) 
 			{
 				locdVdx[i][k] *= locphi[i][k];
 				locdVdy[i][k] *= locphi[i][k];
 				locdVdz[i][k] *= locphi[i][k];
 			}
 		}
 	}

	delete Vmin;
	delete Vmax;
	delete phi;
	delete rf;

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computePressureSensor(double threshold, DistSVec<double,3>& X,
				   DistSVec<double,dim>& V, DistSVec<double,dim>& dVdx,
				   DistSVec<double,dim>& dVdy, DistSVec<double,dim>& dVdz,
				   DistSVec<double,3>& sensor, DistVec<double>& sigma)
{

  const int nsmooth = 2;
  const double gamma = 0.01;
  const double eps = 0.5;
  const double omeps = 1.0 - eps;

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computePressureSensor(X(iSub), V(iSub), dVdx(iSub), dVdy(iSub),
					   dVdz(iSub), sensor(iSub));
    subDomain[iSub]->sndData(*vec3DPat, sensor.subData(iSub));
  }
  vec3DPat->exchange();
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->addRcvData(*vec3DPat, sensor.subData(iSub));
    double (*s)[3] = sensor.subData(iSub);
    double* sig = sigma.subData(iSub);
    for (int i=0; i<sigma.subSize(iSub); ++i) {
      double alpha = gamma * s[i][2] / (s[i][0] + s[i][1] + s[i][2]);
      sig[i] = s[i][0] / (s[i][1] + alpha * s[i][2]);
    }
  }

  for (int k=0; k<nsmooth; ++k) {
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->computeSmoothedSensor(X(iSub), sigma(iSub), sensor(iSub));
      subDomain[iSub]->sndData(*vec3DPat, sensor.subData(iSub));
    }
    vec3DPat->exchange();
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->addRcvData(*vec3DPat, sensor.subData(iSub));
      double (*s)[3] = sensor.subData(iSub);
      double* sig = sigma.subData(iSub);
      for (int i=0; i<sigma.subSize(iSub); ++i)
	sig[i] = eps * sig[i] + omeps * s[i][0] / s[i][1];
    }
  }

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    double (*locdVdx)[dim] = dVdx.subData(iSub);
    double (*locdVdy)[dim] = dVdy.subData(iSub);
    double (*locdVdz)[dim] = dVdz.subData(iSub);
    double* sig = sigma.subData(iSub);
    for (int i=0; i<sigma.subSize(iSub); ++i) {
      if (sig[i] >= threshold) {
	sig[i] = 1.0;
	for (int k=0; k<dim; ++k) {
	  locdVdx[i][k] = 0.0;
	  locdVdy[i][k] = 0.0;
	  locdVdz[i][k] = 0.0;
	}
      }
      else
	sig[i] = 0.0;
    }
  }

}
//------------------------------------------------------------------------------

//d2d BF
template<int dim>
void Domain::computeFiniteVolumeTerm(DistVec<double> &ctrlVol, DistVec<double>& irey,
                                     FluxFcn** fluxFcn, RecFcn* recFcn,
                                     DistBcData<dim>& bcData, DistGeoState& geoState,
                                     DistSVec<double,3>& X, DistSVec<double,dim>& V,
                                     DistNodalGrad<dim>& ngrad, DistEdgeGrad<dim>* egrad,
                                     DistSVec<double,dim>& R, int failsafe, int rshift)
{

  double t0 = timer->getTime();
  int ierr = 0;

  if (!tag) {
     tag = new DistSVec<int,2>(getNodeDistInfo());
     *tag = 0;
  }

  //KW&AM: TODO: should add RR as a member.
  DistSVec<double,dim>* RR = new DistSVec<double,dim>(getNodeDistInfo());
  *RR = R; // initialize temp residual

  int iSub;
#pragma omp parallel for reduction(+: ierr)
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
    ierr += subDomain[iSub]->computeFiniteVolumeTerm(irey(iSub), fluxFcn, recFcn, bcData(iSub),
                                                     geoState(iSub), X(iSub), V(iSub), ngrad(iSub),
                                                     legrad, (*RR)(iSub), (*tag)(iSub), failsafe, rshift);
  }

  com->globalSum(1, &ierr);

  if (ierr) {
    if (!failsafe) {
      com->fprintf(stderr," ... Error: some reconstructed pressure & density are negative. Aborting....\n");
      MPI_Abort(com->comm,-1);
      exit(1);
    }
    else {   // If failsafe option is Yes or Always

      *RR = R; // reinitialize temp residual

#pragma omp parallel for
      for(iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->sndData(*fsPat,  (*tag).subData(iSub));

      fsPat->exchange();

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->addRcvData(*fsPat, (*tag).subData(iSub));

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->finalizeTags((*tag)(iSub));

      ngrad.fix(*tag);
      ngrad.compute(geoState.getConfig(), X, ctrlVol, V);
      ngrad.limit(recFcn, X, ctrlVol, V);

      if (egrad) egrad->fix(*tag);

#pragma omp parallel for reduction(+: ierr)
      for (iSub = 0; iSub < numLocSub; ++iSub) {
        EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
        subDomain[iSub]->computeFiniteVolumeTerm(irey(iSub), fluxFcn, recFcn, bcData(iSub),
                                                 geoState(iSub), X(iSub), V(iSub), ngrad(iSub),
                                                 legrad, (*RR)(iSub), (*tag)(iSub), 0, rshift);
      }

      if (failsafe == 1) *tag = 0;
    }
  }

#pragma omp parallel for reduction(+: ierr)
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->sndData(*vecPat, (*RR).subData(iSub));

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, (*RR).subData(iSub));

  R = *RR;

  timer->addFiniteVolumeTermTime(t0);

  if (RR) delete(RR); // delete temp residual

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeFiniteVolumeTerm(DistExactRiemannSolver<dim>& riemann,
                                     DistVec<double> &ctrlVol, DistVec<double>& irey,
                                     FluxFcn** fluxFcn, RecFcn* recFcn,
                                     DistBcData<dim>& bcData, DistGeoState& geoState,
                                     DistSVec<double,3>& X, DistSVec<double,dim>& V,
                                     DistNodalGrad<dim>& ngrad, DistEdgeGrad<dim>* egrad,
                                     DistSVec<double,dim>& R, int failsafe, int rshift)
{
  double t0 = timer->getTime();
  int ierr = 0;

  if (!tag) {
     tag = new DistSVec<int,2>(getNodeDistInfo());
     *tag = 0;
  }

  DistSVec<double,dim>* RR = new DistSVec<double,dim>(getNodeDistInfo());
  *RR = R; // initialize temp residual

  int iSub;
#pragma omp parallel for reduction(+: ierr)
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
    ierr += subDomain[iSub]->computeFiniteVolumeTerm(riemann(iSub), irey(iSub), fluxFcn, recFcn, bcData(iSub),
                                                     geoState(iSub), X(iSub), V(iSub), ngrad(iSub),
                                                     legrad, (*RR)(iSub), (*tag)(iSub), failsafe, rshift);
  }

  com->globalSum(1, &ierr);

  if (ierr) {
    if (!failsafe) {
      com->fprintf(stderr," ... Error: some reconstructed pressure & density are negative. Aborting....\n");
      MPI_Abort(com->comm,-1);
      exit(1);
    }
    else {   // If failsafe option is Yes or Always

      *RR = R; // reinitialize temp residual

#pragma omp parallel for
      for(iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->sndData(*fsPat,  (*tag).subData(iSub));

      fsPat->exchange();

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->addRcvData(*fsPat, (*tag).subData(iSub));

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->finalizeTags((*tag)(iSub));

      ngrad.fix(*tag);
      ngrad.compute(geoState.getConfig(), X, ctrlVol, V); // bug: this is for one-phase flow!
      //ngrad.compute(geoState.getConfig(), X, ctrlVol, fluidId, V); //where is fluidId?
      ngrad.limit(recFcn, X, ctrlVol, V);

      if (egrad) egrad->fix(*tag);

#pragma omp parallel for reduction(+: ierr)
      for (iSub = 0; iSub < numLocSub; ++iSub) {
        EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
        subDomain[iSub]->computeFiniteVolumeTerm(riemann(iSub), irey(iSub), fluxFcn, recFcn, bcData(iSub),
                                                 geoState(iSub), X(iSub), V(iSub), ngrad(iSub),
                                                 legrad, (*RR)(iSub), (*tag)(iSub), 0, rshift);
      }

      if (failsafe == 1) *tag = 0;
    }
  }

#pragma omp parallel for reduction(+: ierr)
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->sndData(*vecPat, (*RR).subData(iSub));

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, (*RR).subData(iSub));

  R = *RR;

  timer->addFiniteVolumeTermTime(t0);

  if (RR) delete(RR); // delete temp residual

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Domain::computeDerivativeOfFiniteVolumeTerm(DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol,
						 DistVec<double>& irey, DistVec<double>& dIrey,
						 FluxFcn** fluxFcn, RecFcn* recFcn,
						 DistBcData<dim>& bcData, DistGeoState& geoState,
						 DistSVec<double,3>& X, DistSVec<double,3>& dX, DistSVec<double,dim>& V, DistSVec<double,dim>& dV,
						 DistNodalGrad<dim>& ngrad, DistEdgeGrad<dim>* egrad, double dMach,
						 DistSVec<double,dim>& dF)
{

  double t0 = timer->getTime();

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
    subDomain[iSub]->computeDerivativeOfFiniteVolumeTerm(irey(iSub), dIrey(iSub), fluxFcn, recFcn, bcData(iSub), geoState(iSub),
                                                         X(iSub), dX(iSub), V(iSub), dV(iSub), ngrad(iSub), legrad, dMach, dF(iSub));
    subDomain[iSub]->sndData(*vecPat, dF.subData(iSub));

  }

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, dF.subData(iSub));

  timer->addFiniteVolumeTermTime(t0);

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void Domain::computeDerivativeOfFiniteVolumeTerm(dRdXoperators<dim> &dRdXop, 
						 DistBcData<dim>& bcData, DistGeoState& geoState, DistSVec<double,3>& dX, 
						 DistNodalGrad<dim>& ngrad, DistEdgeGrad<dim>* egrad, 
             DistSVec<double,dim>& dddx,
             DistSVec<double,dim>& dddy, 
             DistSVec<double,dim>& dddz,
             DistVec<Vec3D>& dEdgeNormal,
             DistVec<Vec3D>& dFaceNormal,
             DistVec<double>& dFaceNormalVel,
             DistSVec<double,dim>& dF)
{

  double t0 = timer->getTime();

  int iSub;
//  com->fprintf(stderr, " computeDerivativeOfFiniteVolumeTerm received dddx, dddy, dddz norms are %e %e %e\n", dddx.norm(), dddy.norm(), dddz.norm());
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
    subDomain[iSub]->computeDerivativeOfFiniteVolumeTerm(dRdXop.dFluxdddx[iSub], dRdXop.dFluxdddy[iSub], dRdXop.dFluxdddz[iSub],
                                                         dRdXop.dFluxdEdgeNorm[iSub], dRdXop.dFluxdX[iSub],
                                                         dRdXop.dFluxdFaceNormal[iSub], dRdXop.dFluxdFaceNormalVel[iSub], dRdXop.dFluxdUb[iSub],
                                                         bcData(iSub), geoState(iSub),
		                                                     dX(iSub), ngrad(iSub), legrad, 
                                                         dddx(iSub), dddy(iSub), dddz(iSub), dEdgeNormal(iSub), dFaceNormal(iSub), dFaceNormalVel(iSub), dF(iSub));
  }
  assemble(vecPat, dF);

  timer->addFiniteVolumeTermTime(t0);

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void Domain::computeTransposeDerivativeOfFiniteVolumeTerm(dRdXoperators<dim> &dRdXop, 
						 DistBcData<dim>& bcData, DistGeoState& geoState,
						 DistSVec<double,dim>& dF, 
						 DistNodalGrad<dim>& ngrad, DistEdgeGrad<dim>* egrad, 
						 DistSVec<double,3>& dX2,
             DistSVec<double,dim>& dddx2,
             DistSVec<double,dim>& dddy2,
             DistSVec<double,dim>& dddz2,
             DistVec<Vec3D>& dEdgeNormal2,
             DistVec<Vec3D>& dFaceNormal2,
             DistVec<double>& dFaceNormalVel2)
{

  double t0 = timer->getTime();
  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
    subDomain[iSub]->computeTransposeDerivativeOfFiniteVolumeTerm(dRdXop.dFluxdddx[iSub], dRdXop.dFluxdddy[iSub], dRdXop.dFluxdddz[iSub],
                                                                  dRdXop.dFluxdEdgeNorm[iSub], dRdXop.dFluxdX[iSub],
                                                                  dRdXop.dFluxdFaceNormal[iSub], dRdXop.dFluxdFaceNormalVel[iSub], dRdXop.dFluxdUb[iSub],
                                                                  bcData(iSub), geoState(iSub), dF(iSub),
                                                                  ngrad(iSub), legrad, dX2(iSub), dddx2(iSub), dddy2(iSub), dddz2(iSub), dEdgeNormal2(iSub), dFaceNormal2(iSub), dFaceNormalVel2(iSub));


  }

  DistVec<double> &dEdgeNormVel = geoState.getdEdgeNormalVel();
  DistVec<double> dEdgeNormVel2(dEdgeNormVel);
  dEdgeNormVel2 = 0.0;
  assemble(vecPat, dddx2);
  assemble(vecPat, dddy2);
  assemble(vecPat, dddz2);

#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) 
    subDomain[iSub]->sndNormals(*edgePat, dEdgeNormal2.subData(iSub), dEdgeNormVel2.subData(iSub));

  edgePat->exchange();

#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->rcvNormals(*edgePat, dEdgeNormal2.subData(iSub), dEdgeNormVel2.subData(iSub));

  timer->addFiniteVolumeTermTime(t0);

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void Domain::computeDerivativeOperatorsOfFiniteVolumeTerm(DistVec<double>& irey, DistVec<double>& dIrey,
                                                          FluxFcn** fluxFcn, RecFcn* recFcn,
                                                          DistBcData<dim>& bcData, DistGeoState& geoState,
                                                          DistSVec<double,3>& X, DistSVec<double,dim>& V,
                                                          DistNodalGrad<dim>& ngrad, DistEdgeGrad<dim>* egrad, double dMach,
                                                          dRdXoperators<dim> &dRdXop)
{

  double t0 = timer->getTime();

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
    subDomain[iSub]->computeDerivativeOperatorsOfFiniteVolumeTerm(irey(iSub), dIrey(iSub), fluxFcn, recFcn, bcData(iSub), geoState(iSub),
                                                                  X(iSub), V(iSub), ngrad(iSub), legrad, dMach, 
                                                                  *dRdXop.dFluxdEdgeNorm[iSub], 
                                                                  *dRdXop.dFluxdX[iSub],
                                                                  *dRdXop.dFluxdddx[iSub],
                                                                  *dRdXop.dFluxdddy[iSub],
                                                                  *dRdXop.dFluxdddz[iSub],
                                                                  *dRdXop.dFluxdFaceNormal[iSub],
                                                                  *dRdXop.dFluxdFaceNormalVel[iSub],
                                                                  *dRdXop.dFluxdUb[iSub]);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeDerivativeOfFiniteVolumeTerm(FluxFcn** fluxFcn, RecFcn* recFcn,
						 DistBcData<dim>& bcData, DistGeoState& geoState,
						 DistSVec<double,3> &X,
						 DistLevelSetStructure *distLSS,
						 bool linRecAtInterface, bool viscSecOrder, 
						 DistVec<int> &fluidId, 
						 DistExactRiemannSolver<dim> &riemann, 
						 int Nriemann,
						 DistNodalGrad<dim>& ngrad, 
						 DistEdgeGrad<dim>* egrad,
						 double dMach,
						 DistSVec<double,dim>& V,
						 DistSVec<double,dim>& dF)
{

  double t0 = timer->getTime();

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {

    EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;

    Vec<int> &FluidId = fluidId(iSub);

    subDomain[iSub]->computeDerivativeOfFiniteVolumeTerm(fluxFcn, recFcn, bcData(iSub), geoState(iSub),
							 X(iSub), (*distLSS)(iSub), 
							 linRecAtInterface,  viscSecOrder, 
							 FluidId, riemann(iSub), Nriemann,
							 ngrad(iSub), legrad, dMach,
							 V(iSub), dF(iSub));

    subDomain[iSub]->sndData(*vecPat, dF.subData(iSub));

  }

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, dF.subData(iSub));

  timer->addFiniteVolumeTermTime(t0);

}
//------------------------------------------------------------------------------


template<int dim, int dimLS>
void Domain::computeFiniteVolumeTerm(DistVec<double> &ctrlVol,
                                     DistExactRiemannSolver<dim> &riemann,
                                     FluxFcn** fluxFcn, RecFcn* recFcn,
                                     DistBcData<dim>& bcData, DistGeoState& geoState,
                                     DistSVec<double,3>& X, DistSVec<double,dim>& V,
                                     FluidSelector &fluidSelector,
                                     DistNodalGrad<dim>& ngrad, DistEdgeGrad<dim>* egrad,
				     DistSVec<double,dimLS>& phi,
                                     DistNodalGrad<dimLS>& ngradLS,
				     DistEdgeGrad<dimLS>* egradLS,
                                     DistSVec<double,dim>& R, int it,
                                     int failsafe, int rshift)
{

  double t0 = timer->getTime();
  int ierr = 0;

  if (!tag) {
     tag = new DistSVec<int,2>(getNodeDistInfo());
     *tag = 0;
  }

  DistSVec<double,dim>* RR = new DistSVec<double,dim>(getNodeDistInfo());
  *RR = R; // initialize temp residual
  DistVec<int> &FluidId(*(fluidSelector.fluidId));

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {

    EdgeGrad<dim>*   legrad   = (egrad)   ? &((*egrad)(iSub))   : 0;
    EdgeGrad<dimLS>* legradLS = (egradLS) ? &((*egradLS)(iSub)) : 0;

    Vec<int> &fluidId = FluidId(iSub);
    ierr = subDomain[iSub]->computeFiniteVolumeTerm(riemann(iSub),
						    fluxFcn, recFcn, bcData(iSub), geoState(iSub),
						    X(iSub), V(iSub), fluidId,
						    fluidSelector, 
						    ngrad(iSub),   legrad, phi(iSub),  
						    ngradLS(iSub), legradLS,
						    (*RR)(iSub), it,
						    (*tag)(iSub), failsafe, rshift);
  }
  com->globalSum(1, &ierr);

  if (ierr) {
    if (!failsafe) {
      com->fprintf(stderr," ... Error: some reconstructed pressure & density are negative. Aborting....\n");
      MPI_Abort(com->comm,-1);
      exit(1);
    }
    else {   // If failsafe option is Yes or Always

      *RR = R; // reinitialize temp residual

#pragma omp parallel for
      for(iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->sndData(*fsPat,  (*tag).subData(iSub));

      fsPat->exchange();

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->addRcvData(*fsPat, (*tag).subData(iSub));

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->finalizeTags((*tag)(iSub));

      ngrad.fix(*tag);
      ngrad.compute(geoState.getConfig(), X, ctrlVol, FluidId, V);
      ngrad.limit(recFcn, X, ctrlVol, V);

      if (egrad) egrad->fix(*tag);
//    if (egradLS) egradLS->fix(*tag);

#pragma omp parallel for reduction(+: ierr)
      for (iSub = 0; iSub < numLocSub; ++iSub) {

        EdgeGrad<dim>*   legrad   = (egrad)   ? &((*egrad)(iSub))   : 0;
        EdgeGrad<dimLS>* legradLS = (egradLS) ? &((*egradLS)(iSub)) : 0;

        Vec<int> &fluidId = FluidId(iSub);
        ierr = subDomain[iSub]->computeFiniteVolumeTerm(riemann(iSub),
							fluxFcn, recFcn, bcData(iSub), geoState(iSub),
							X(iSub), V(iSub), fluidId,
							fluidSelector, 
							ngrad(iSub),   legrad, phi(iSub), 
							ngradLS(iSub), legradLS, 
							(*RR)(iSub), it,
							(*tag)(iSub), 0, rshift);
      }

      if (failsafe == 1) *tag = 0;
    }
  }

#pragma omp parallel for reduction(+: ierr)
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->sndData(*vecPat, (*RR).subData(iSub));

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, (*RR).subData(iSub));

  R = *RR;

  timer->addFiniteVolumeTermTime(t0);

  if (RR) delete(RR); // delete temp residual

// subdomain communication for riemann update values (cf ExactRiemannSolver.h)
  if(it == 1){
    DistSVec<double,dim> *rupdate = riemann.getRiemannUpdate();
    DistVec<double> *weight= riemann.getRiemannWeight();
    DistVec<int>* fid = riemann.getFluidIdToSet();
    assemble(vecPat,*rupdate);
    assemble(volPat,*weight);
    operMax<int> opMax;
    assemble(levelPat, *fid, opMax);
  }

}

//------------------------------------------------------------------------------
//d2d embedded
template<int dim, int dimLS>
void Domain::computeFiniteVolumeTerm(DistVec<double> &ctrlVol, 
				     DistExactRiemannSolver<dim> &riemann,
                                     FluxFcn** fluxFcn, RecFcn* recFcn, 
				     DistBcData<dim>& bcData, 
				     DistGeoState& geoState,
                                     DistSVec<double,3>& X, 
				     DistSVec<double,dim>& V, 
				     DistSVec<double,dim>& Wstarij, DistSVec<double,dim>& Wstarji,
                                     DistLevelSetStructure *distLSS, bool linRecAtInterface, 
				     FluidSelector &fluidSelector, int Nriemann,
				     DistNodalGrad<dim>& ngrad, DistEdgeGrad<dim>* egrad,
				     DistSVec<double,dimLS>& phi,
                                     DistNodalGrad<dimLS>& ngradLS, DistEdgeGrad<dimLS>* egradLS, 
				     DistSVec<double,dim>& R, int it, int failsafe, int rshift)
{

 double t0 = timer->getTime();
  int ierr = 0;

  if (!tag) {
     tag = new DistSVec<int,2>(getNodeDistInfo());
     *tag = 0;
  }

  DistSVec<double,dim>* RR = new DistSVec<double,dim>(getNodeDistInfo());
  *RR = R; // initialize temp residual
  DistVec<int> &FluidId(*(fluidSelector.fluidId));

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {

    EdgeGrad<dim>*   legrad   = (egrad)   ? &((*egrad)(iSub))   : 0;
    EdgeGrad<dimLS>* legradLS = (egradLS) ? &((*egradLS)(iSub)) : 0;

    Vec<int> &fluidId = FluidId(iSub);
    ierr = subDomain[iSub]->computeFiniteVolumeTerm(riemann(iSub),
						    fluxFcn, recFcn, bcData(iSub), geoState(iSub),
						    X(iSub), V(iSub), Wstarij(iSub), Wstarji(iSub), (*distLSS)(iSub), 
						    linRecAtInterface, fluidId, Nriemann,
						    fluidSelector, 
						    ngrad(iSub),   legrad, phi(iSub), 
						    ngradLS(iSub), legradLS, 
						    (*RR)(iSub), it,
						    (*tag)(iSub), failsafe, rshift);
  }
  com->globalSum(1, &ierr);

  // failsafe
  if (ierr) {
    if (!failsafe) {
      com->fprintf(stderr," ... Error: some reconstructed pressure & density are negative. Aborting....\n");
      MPI_Abort(com->comm,-1);
      exit(1);
    }
    else {   // If failsafe option is Yes or Always

      *RR = R; // reinitialize temp residual

#pragma omp parallel for
      for(iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->sndData(*fsPat,  (*tag).subData(iSub));

      fsPat->exchange();

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->addRcvData(*fsPat, (*tag).subData(iSub));

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->finalizeTags((*tag)(iSub));

      ngrad.fix(*tag);
      ngrad.compute(geoState.getConfig(), X, ctrlVol, FluidId, V);
      ngrad.limit(recFcn, X, ctrlVol, V);

      if (egrad) egrad->fix(*tag);
    //if (egradLS) egradLS->fix(*tag);

#pragma omp parallel for reduction(+: ierr)
      for (iSub = 0; iSub < numLocSub; ++iSub) {

        EdgeGrad<dim>*   legrad   = (egrad)   ? &((*egrad)(iSub))   : 0;
        EdgeGrad<dimLS>* legradLS = (egradLS) ? &((*egradLS)(iSub)) : 0;

        Vec<int> &fluidId = FluidId(iSub);

        ierr = subDomain[iSub]->computeFiniteVolumeTerm(riemann(iSub),
							fluxFcn, recFcn, bcData(iSub), geoState(iSub), 
							X(iSub), V(iSub), Wstarij(iSub), Wstarji(iSub),(*distLSS)(iSub), 
							linRecAtInterface, fluidId, Nriemann,
							fluidSelector, 
							ngrad(iSub),   legrad,  phi(iSub), 
							ngradLS(iSub), legradLS,
							(*RR)(iSub), it,
							(*tag)(iSub), 0, rshift);
      }

      if (failsafe == 1) *tag = 0;
    }
  }

#pragma omp parallel for reduction(+: ierr)
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->sndData(*vecPat, (*RR).subData(iSub));

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, (*RR).subData(iSub));

  R = *RR;
  if (RR) delete(RR); // delete temp residual

// subdomain communication for riemann update values (cf ExactRiemannSolver.h)
  if(it == 1){
    DistSVec<double,dim> *rupdate = riemann.getRiemannUpdate();
    DistVec<double> *weight= riemann.getRiemannWeight();
    DistVec<int>* fid = riemann.getFluidIdToSet();
    assemble(vecPat,*rupdate);
    assemble(volPat,*weight);
    operMax<int> opMax;
    assemble(levelPat, *fid, opMax);
  }

  timer->addFiniteVolumeTermTime(t0);
}

//------------------------------------------------------------------------------

//d2d embedded structure
template<int dim>
void Domain::computeFiniteVolumeTerm(DistVec<double> &ctrlVol,
                                     DistExactRiemannSolver<dim> &riemann,
                                     FluxFcn** fluxFcn, RecFcn* recFcn,
                                     DistBcData<dim>& bcData, 
				     DistGeoState& geoState,
                                     DistSVec<double,3>& X, 
				     DistSVec<double,dim>& V,
                                     DistSVec<double,dim>& Wstarij, DistSVec<double,dim>& Wstarji,
												 DistSVec<double,dim>& Wext,
                                     DistLevelSetStructure *LSS, bool linRecAtInterface, DistVec<int> &fluidId, 
                                     int Nriemann, 
				     DistNodalGrad<dim>& ngrad, DistEdgeGrad<dim>* egrad,
                                     DistSVec<double,dim>& R, 
												 int it, int failsafe, int rshift, bool externalSI)
{
 
  double t0 = timer->getTime();
  int ierr = 0;

  if (!tag) {
     tag = new DistSVec<int,2>(getNodeDistInfo());
     *tag = 0;
  }

  DistSVec<double,dim>* RR = new DistSVec<double,dim>(getNodeDistInfo());
  *RR = R; // initialize temp residual

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) 
  {
    EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
    ierr = subDomain[iSub]->computeFiniteVolumeTerm(riemann(iSub),
						    fluxFcn, recFcn, 
						    bcData(iSub), geoState(iSub),
						    X(iSub), V(iSub), 
																	  Wstarij(iSub), Wstarji(iSub), Wext(iSub),
						    (*LSS)(iSub),
						    linRecAtInterface, fluidId(iSub), Nriemann, 
						    ngrad(iSub), legrad, 
						    (*RR)(iSub), 
																	  it, (*tag)(iSub), failsafe, rshift, externalSI);
  }

  com->globalSum(1, &ierr);

  if(ierr) 
  {	  
	  if (!failsafe) 
	  {
      com->fprintf(stderr," ... Error: some reconstructed pressure & density are negative. Aborting....\n");
      MPI_Abort(com->comm,-1);
      exit(1);
	  } 
	  else 
	  {
		  
		  // If failsafe option is Yes or Always

      *RR = R; // reinitialize temp residual

#pragma omp parallel for
      for(iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->sndData(*fsPat,  (*tag).subData(iSub));

      fsPat->exchange();

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->addRcvData(*fsPat, (*tag).subData(iSub));

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->finalizeTags((*tag)(iSub));

      ngrad.fix(*tag);
      ngrad.compute(geoState.getConfig(), X, ctrlVol, fluidId, V);
      ngrad.limit(recFcn, X, ctrlVol, V);

      if (egrad) egrad->fix(*tag);

#pragma omp parallel for reduction(+: ierr)
		  for (iSub = 0; iSub < numLocSub; ++iSub) 
		  {

        EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
        ierr = subDomain[iSub]->computeFiniteVolumeTerm(riemann(iSub),
							fluxFcn, recFcn, 
							bcData(iSub), geoState(iSub),
							X(iSub), V(iSub), 
																			  Wstarij(iSub), Wstarji(iSub), Wext(iSub),
							(*LSS)(iSub), 
							linRecAtInterface, fluidId(iSub), Nriemann, 
							ngrad(iSub), legrad, (*RR)(iSub), 
																			  it, (*tag)(iSub), 0, rshift, externalSI);
		  }
      if (failsafe == 1) *tag = 0;
    }
  }

#pragma omp parallel for reduction(+: ierr)
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->sndData(*vecPat, (*RR).subData(iSub));

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, (*RR).subData(iSub));

  R = *RR;

#pragma omp parallel for reduction(+: ierr)
  for (iSub = 0; iSub < numLocSub; ++iSub)
	  subDomain[iSub]->sndData(*vecPat, Wext.subData(iSub));

  vecPat->exchange();

#pragma omp parallel for
   for (iSub = 0; iSub < numLocSub; ++iSub)
		subDomain[iSub]->RcvData(*vecPat, Wext.subData(iSub));

  timer->addFiniteVolumeTermTime(t0);

  if (RR) delete(RR); // delete temp residual
}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeFiniteVolumeTerm(DistVec<double> &ctrlVol,
                                     DistExactRiemannSolver<dim> &riemann,
                                     FluxFcn** fluxFcn, RecFcn* recFcn,
                                     DistBcData<dim>& bcData, DistGeoState& geoState,
                                     DistSVec<double,3>& X, DistSVec<double,dim>& V,
                                     DistSVec<double,dim>& Wstarij, 
				     DistSVec<double,dim>& Wstarji,
				     DistVec<int>& countWstarij, 
				     DistVec<int>& countWstarji,
                                     DistLevelSetStructure *LSS, bool linRecAtInterface, 
				     DistVec<int> &fluidId, int Nriemann, 
				     double dt, double alpha, 
				     DistNodalGrad<dim>& ngrad, DistEdgeGrad<dim>* egrad, 
				     DistSVec<double,dim>& R, int it, int failsafe, int rshift)
{

  double t0 = timer->getTime();
  int ierr = 0;

  if (!tag) {
     tag = new DistSVec<int,2>(getNodeDistInfo());
     *tag = 0;
  }

  DistSVec<double,dim>* RR = new DistSVec<double,dim>(getNodeDistInfo());
  *RR = R; // initialize temp residual

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
    ierr = subDomain[iSub]->computeFiniteVolumeTerm(riemann(iSub),
						    fluxFcn, recFcn, bcData(iSub), geoState(iSub),
						    X(iSub), V(iSub), Wstarij(iSub), Wstarji(iSub), 
						    countWstarij(iSub), countWstarji(iSub), 
						    (*LSS)(iSub), linRecAtInterface, fluidId(iSub), 
						    Nriemann,
						    dt, alpha, ngrad(iSub), legrad, (*RR)(iSub), it,
						    (*tag)(iSub), failsafe, rshift);
  }
  com->globalSum(1, &ierr);

  if (ierr) {
    if (!failsafe) {
      com->fprintf(stderr," ... Error: some reconstructed pressure & density are negative. Aborting....\n");
      MPI_Abort(com->comm,-1);
      exit(1);
    }
    else {   // If failsafe option is Yes or Always

      *RR = R; // reinitialize temp residual

#pragma omp parallel for
      for(iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->sndData(*fsPat,  (*tag).subData(iSub));

      fsPat->exchange();

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->addRcvData(*fsPat, (*tag).subData(iSub));

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->finalizeTags((*tag)(iSub));

      ngrad.fix(*tag);
      ngrad.compute(geoState.getConfig(), X, ctrlVol, fluidId, V);
      ngrad.limit(recFcn, X, ctrlVol, V);

      if (egrad) egrad->fix(*tag);

#pragma omp parallel for reduction(+: ierr)
      for (iSub = 0; iSub < numLocSub; ++iSub) {
        EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
        ierr = subDomain[iSub]->computeFiniteVolumeTerm(riemann(iSub),fluxFcn, recFcn, 
							bcData(iSub), geoState(iSub), X(iSub), V(iSub), 
							Wstarij(iSub), Wstarji(iSub), 
							countWstarij(iSub), countWstarji(iSub), 
							(*LSS)(iSub), linRecAtInterface, fluidId(iSub), Nriemann, 
							dt, alpha, 
							ngrad(iSub), legrad, (*RR)(iSub), it, 
							(*tag)(iSub), 0, rshift); 
      }

      if (failsafe == 1) *tag = 0;
    }
  }

#pragma omp parallel for reduction(+: ierr)
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->sndData(*vecPat, (*RR).subData(iSub));

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, (*RR).subData(iSub));

  R = *RR;

  timer->addFiniteVolumeTermTime(t0);

  if (RR) delete(RR); // delete temp residual
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void Domain::computeFiniteVolumeTermLS(FluxFcn** fluxFcn, RecFcn* recFcn, RecFcn* recFcnLS,
				       DistBcData<dim>& bcData, DistGeoState& geoState,
				       DistSVec<double,3>& X, DistSVec<double,dim>& V,
                                       DistVec<int>& fluidId,
				       DistNodalGrad<dim>& ngrad,     DistEdgeGrad<dim>* egrad,
				       DistNodalGrad<dimLS>& ngradLS, DistEdgeGrad<dimLS>* egradLS,
				       DistSVec<double,dimLS>& Phi, DistSVec<double,dimLS> &PhiF,
				       DistLevelSetStructure *distLSS, int ls_order)
{

  double t0 = timer->getTime();
  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {

    EdgeGrad<dim>*   legrad   = (egrad)   ? &((*egrad)(iSub))   : 0;
    EdgeGrad<dimLS>* legradLS = (egradLS) ? &((*egradLS)(iSub)) : 0;

    LevelSetStructure* LSS = (distLSS) ? &((*distLSS)(iSub)) : 0;

    subDomain[iSub]->computeFiniteVolumeTermLS(fluxFcn, recFcn, recFcnLS, bcData(iSub),
                                               geoState(iSub),
                                               X(iSub), V(iSub),fluidId(iSub), 
					       ngrad(iSub),   legrad,
					       ngradLS(iSub), legradLS,
                                               Phi(iSub),PhiF(iSub), LSS, ls_order);
    subDomain[iSub]->sndData(*phiVecPat, PhiF.subData(iSub));
  }

  phiVecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*phiVecPat, PhiF.subData(iSub));

  timer->addLSFiniteVolumeTermTime(t0);
}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeFiniteVolumeBarTerm(DistVec<double> &ctrlVol,
                                        DistVec<double> &irey, FluxFcn** fluxFcn, RecFcn* recFcn,
                                        DistBcData<dim>& bcData, DistGeoState& geoState,
                                        DistSVec<double,3>& X, DistMacroCellSet *macroCells,
                                        DistSVec<double,dim> &VBar, DistSVec<double,1> &volRatio,
                                        DistNodalGrad<dim>& ngrad, DistEdgeGrad<dim>* egrad,
                                        DistSVec<double,dim>& RBar, int scopeDepth1, int scopeDepth2,
                                        int failsafe, int rshift)
{

  double t0 = timer->getTime();
  int ierr = 0;

  if (!tagBar) {
     tagBar = new DistSVec<int,2>(getNodeDistInfo());
     *tagBar = 0;
  }

  DistSVec<double,dim>* Sigma = new DistSVec<double,dim>(getNodeDistInfo());
  *Sigma = 0.0;
  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
    ierr += subDomain[iSub]->computeFiniteVolumeBar_Step1(irey(iSub), fluxFcn, recFcn, bcData(iSub), geoState(iSub),
                                                          X(iSub), VBar(iSub), ngrad(iSub), legrad, (*Sigma)(iSub),
                                                          (*tagBar)(iSub), failsafe, rshift);
  }

  com->globalSum(1, &ierr);

  if (ierr) {
    if (!failsafe) {
      com->fprintf(stderr," ... Error: some reconstructed pressure & density are negative. Aborting....\n");
      MPI_Abort(com->comm,-1);
      exit(1);
    }
    else {   // If failsafe option is Yes or Always

      *Sigma = 0.0; // reinitialize temp residual

#pragma omp parallel for
      for(iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->sndData(*fsPat,  (*tagBar).subData(iSub));

      fsPat->exchange();

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->addRcvData(*fsPat, (*tagBar).subData(iSub));

#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->finalizeTags((*tagBar)(iSub));

      ngrad.fix(*tagBar);
      ngrad.compute(geoState.getConfig(), X, ctrlVol, VBar);
      ngrad.limit(recFcn, X, ctrlVol, VBar);

      if (egrad) egrad->fix(*tagBar);

#pragma omp parallel for reduction(+: ierr)
      for (iSub = 0; iSub < numLocSub; ++iSub) {
        EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
        subDomain[iSub]->computeFiniteVolumeBar_Step1(irey(iSub), fluxFcn, recFcn, bcData(iSub), geoState(iSub),
                                                      X(iSub), VBar(iSub), ngrad(iSub), legrad, (*Sigma)(iSub),
                                                      (*tagBar)(iSub), 0, rshift);
      }

      if (failsafe == 1) *tagBar = 0;
    }
  }

#pragma omp parallel for reduction(+: ierr)
  for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->sndData(*vecPat, (*Sigma).subData(iSub));

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, (*Sigma).subData(iSub));


#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    MacroCellSet** macCells = new MacroCellSet*[scopeDepth2];

    for (int i = 0; i<scopeDepth2; ++i)
        macCells[i] =   macroCells->obtainMacroCell(iSub, i);

     subDomain[iSub]->computeFiniteVolumeBar_Step2(macCells, volRatio(iSub),
                                   (*Sigma)(iSub), RBar(iSub), scopeDepth1);
     subDomain[iSub]->sndData(*vecPat, RBar.subData(iSub));
     delete [] macCells;
  }

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, RBar.subData(iSub));

  timer->addFiniteVolumeTermTime(t0);

  delete (Sigma);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void Domain::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, DistBcData<dim> &bcData,
                                             DistGeoState &geoState, DistVec<double> &irey,
                                             DistSVec<double,3> &X,
                                             DistVec<double> &ctrlVol,
                                             DistSVec<double,dim> &V, DistMat<Scalar,neq> &A)
{
  int iSub;
  double t0 = timer->getTime();
  CommPattern<Scalar> *matPat = A.getDiagMatPat();

  if(inletRhsPat){
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->computeJacobianFiniteVolumeTerm(fluxFcn, bcData(iSub), geoState(iSub), irey(iSub),
                                                     X(iSub), ctrlVol(iSub), V(iSub), A(iSub), inletRhsPat);
      subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
    }

    double t = timer->addFiniteVolumeJacTime(t0);
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvDiagInletBlocks(*matPat, A(iSub));
    com->printf(6, "FV Jacobian matrix computation: %f s\n", t);


  }else{
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->computeJacobianFiniteVolumeTerm(fluxFcn, bcData(iSub), geoState(iSub), irey(iSub),
                                                     X(iSub), ctrlVol(iSub), V(iSub), A(iSub), inletRhsPat);
      subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
    }

    double t = timer->addFiniteVolumeJacTime(t0);
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvDiagBlocks(*matPat, A(iSub));
    com->printf(6, "FV Jacobian matrix computation: %f s\n", t);

  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void Domain::computeJacobianFiniteVolumeTerm(DistExactRiemannSolver<dim> &riemann,
                                             FluxFcn **fluxFcn, DistBcData<dim> &bcData,
                                             DistGeoState &geoState, DistVec<double> &irey,
                                             DistSVec<double,3> &X,
                                             DistVec<double> &ctrlVol,
                                             DistSVec<double,dim> &V, DistMat<Scalar,neq> &A)
{
  int iSub;
  double t0 = timer->getTime();
  CommPattern<Scalar> *matPat = A.getDiagMatPat();

  if(inletRhsPat){
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->computeJacobianFiniteVolumeTerm(riemann(iSub), fluxFcn, bcData(iSub), geoState(iSub), irey(iSub),
                                                     X(iSub), ctrlVol(iSub), V(iSub), A(iSub), inletRhsPat);
      subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
    }

    double t = timer->addFiniteVolumeJacTime(t0);
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvDiagInletBlocks(*matPat, A(iSub));
    com->printf(6, "FV Jacobian matrix computation: %f s\n", t);


  }else{
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->computeJacobianFiniteVolumeTerm(riemann(iSub), fluxFcn, bcData(iSub), geoState(iSub), irey(iSub),
                                                     X(iSub), ctrlVol(iSub), V(iSub), A(iSub), inletRhsPat);
      subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
    }

    double t = timer->addFiniteVolumeJacTime(t0);
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvDiagBlocks(*matPat, A(iSub));
    com->printf(6, "FV Jacobian matrix computation: %f s\n", t);

  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq, int dimLS>
void Domain::computeJacobianFiniteVolumeTerm(DistExactRiemannSolver<dim> &riemann,
                                             FluxFcn **fluxFcn, DistBcData<dim> &bcData,
                                             DistGeoState &geoState,
                                             DistNodalGrad<dim> &ngrad, DistNodalGrad<dimLS> &ngradLS,
                                             DistSVec<double,3> &X,
                                             DistVec<double> &ctrlVol,
                                             DistSVec<double,dim> &V, DistMat<Scalar,neq> &A,
                                             FluidSelector &fluidSelector)
{

  int iSub;
  double t0 = timer->getTime();
  CommPattern<Scalar> *matPat = A.getDiagMatPat();

  if(inletRhsPat){
    fprintf(stdout, "with inletRhsPat\n");
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->computeJacobianFiniteVolumeTerm(riemann(iSub), fluxFcn,
                                                     bcData(iSub), geoState(iSub),
                                                     ngrad(iSub), ngradLS(iSub),
                                                     X(iSub),
                                                     ctrlVol(iSub), V(iSub), A(iSub),
                                                     fluidSelector, 
                                                     (*(fluidSelector.fluidId))(iSub), inletRhsPat);
      subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
    }
    double t = timer->addFiniteVolumeJacTime(t0);
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvDiagInletBlocks(*matPat, A(iSub));
    com->printf(6, "FV Jacobian matrix computation: %f s\n", t);

  }else{
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->computeJacobianFiniteVolumeTerm(riemann(iSub), fluxFcn,
                                                     bcData(iSub), geoState(iSub),
                                                     ngrad(iSub), ngradLS(iSub),
                                                     X(iSub),
                                                     ctrlVol(iSub), V(iSub), A(iSub),
                                                     fluidSelector, 
                                                     (*(fluidSelector.fluidId))(iSub), inletRhsPat);
      subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
    }
    double t = timer->addFiniteVolumeJacTime(t0);
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvDiagBlocks(*matPat, A(iSub));
    com->printf(6, "FV Jacobian matrix computation: %f s\n", t);
  }
}

//d2d embedded $$
template<class Scalar,int dim,int neq>
void Domain::computeJacobianFiniteVolumeTerm(DistVec<double> &ctrlVol,
                                             DistExactRiemannSolver<dim> &riemann,
                                             FluxFcn** fluxFcn,
                                             DistBcData<dim>& bcData, DistGeoState& geoState,
                                             DistSVec<double,3>& X, DistSVec<double,dim>& V,
                                             DistLevelSetStructure *LSS, DistVec<int> &fluidId, 
                                             int Nriemann,
                                             DistMat<Scalar,neq>& A, DistVec<double>& irey, bool externalSI) 
{

  int iSub;
  double t0 = timer->getTime();
  CommPattern<Scalar> *matPat = A.getDiagMatPat();
 
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->computeJacobianFiniteVolumeTerm(riemann(iSub), fluxFcn,
                                                     bcData(iSub), geoState(iSub),
                                                     X(iSub), V(iSub),ctrlVol(iSub),
                                                     (*LSS)(iSub),
                                                     fluidId(iSub),Nriemann,
																		  A(iSub),irey(iSub), externalSI);
      subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
    }
    double t = timer->addFiniteVolumeJacTime(t0);
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvDiagBlocks(*matPat, A(iSub));
}
  

template<int dim, class Scalar, int neq, int dimLS>
void Domain::computeJacobianFiniteVolumeTerm(DistExactRiemannSolver<dim>& riemann,
                                             FluxFcn** fluxFcn, 
                                             DistBcData<dim>& bcData, DistGeoState& geoState,
                                             DistSVec<double,3>& X, DistSVec<double,dim>& V,DistVec<double>& ctrlVol,
                                             DistNodalGrad<dimLS> &ngradLS,
                                             DistLevelSetStructure *LSS,
                                             int Nriemann,
                                             FluidSelector &fluidSelector,
                                             DistMat<Scalar,neq>& A) {

  int iSub;
  double t0 = timer->getTime();
  CommPattern<Scalar> *matPat = A.getDiagMatPat();
 
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->computeJacobianFiniteVolumeTerm(riemann(iSub), fluxFcn,
                                                     bcData(iSub), geoState(iSub),
                                                     X(iSub), V(iSub),ctrlVol(iSub),
                                                     ngradLS(iSub),(*LSS)(iSub),(*(fluidSelector.fluidId))(iSub),
                                                     Nriemann, fluidSelector,
                                                     A(iSub));
      subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
    }
    double t = timer->addFiniteVolumeJacTime(t0);
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvDiagBlocks(*matPat, A(iSub));

}

//------------------------------------------------------------------------------
template<int dim, class Scalar, int dimLS>
void Domain::computeJacobianFiniteVolumeTermLS(RecFcn* recFcn, RecFcn* recFcnLS,
					   DistGeoState &geoState,DistSVec<double,3>& X,DistSVec<double,dim> &V,
					   DistNodalGrad<dim>& ngrad,DistNodalGrad<dimLS> &ngradLS,
					   DistEdgeGrad<dim>* egrad,
					   DistVec<double> &ctrlVol,DistSVec<double,dimLS>& Phi,
					   DistMat<Scalar,dimLS> &A,DistLevelSetStructure* distLSS)
{

  int iSub;
  double t0 = timer->getTime();
  CommPattern<Scalar> *matPat = A.getDiagMatPat();

  if(inletRhsPat){
    fprintf(stdout, "with inletRhsPat\n");
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
      LevelSetStructure* lss = (distLSS) ? &((*distLSS)(iSub)) : 0;
      subDomain[iSub]->computeJacobianFiniteVolumeTermLS(recFcn,recFcnLS, 
							 geoState(iSub),
							 X(iSub),V(iSub),ngrad(iSub),
							 ngradLS(iSub),
							 legrad,
							 ctrlVol(iSub),
							 Phi(iSub), 
							 A(iSub),lss,
							 inletRhsPat);
      subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
    }
    double t = timer->addLSFiniteVolumeJacTime(t0);
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvDiagInletBlocks(*matPat, A(iSub));
    com->printf(6, "FV Jacobian matrix computation: %f s\n", t);

  }else{
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub) {
      EdgeGrad<dim>* legrad = (egrad) ? &((*egrad)(iSub)) : 0;
      LevelSetStructure* lss = (distLSS) ? &((*distLSS)(iSub)) : 0;
      subDomain[iSub]->computeJacobianFiniteVolumeTermLS(recFcn,recFcnLS, 
							 geoState(iSub),
							 X(iSub),V(iSub),ngrad(iSub),
							 ngradLS(iSub),
							 legrad,
							 ctrlVol(iSub),
							 Phi(iSub), 
							 A(iSub),lss,
							 inletRhsPat);
      subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
    }
    double t = timer->addLSFiniteVolumeJacTime(t0);
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvDiagBlocks(*matPat, A(iSub));
    com->printf(6, "FV Jacobian matrix computation: %f s\n", t);
  }
}

template<int dim>
void Domain::recomputeRHS(VarFcn* vf, DistSVec<double,dim> &V, DistSVec<double,dim> &rhs,
                          DistExtrapolation<dim>* xpol, DistBcData<dim>& bcData,
                          DistGeoState& geoState, DistSVec<double,3> &X)
{
  int iSub;

#pragma omp parallel for
  for ( iSub = 0; iSub < numLocSub; iSub++){
    Extrapolation<dim>* lxpol = (xpol) ? &((*xpol)(iSub)) : 0;
    subDomain[iSub]->recomputeRHS(vf, V(iSub), rhs(iSub), lxpol,
   				 bcData(iSub), geoState(iSub), X(iSub));
    subDomain[iSub]->sndInletRhsData(*inletRhsPat, rhs.subData(iSub));
  }

  inletRhsPat->exchange();

#pragma omp parallel for
  for ( iSub = 0; iSub < numLocSub; iSub++)
    subDomain[iSub]->addRcvInletRhsData(*inletRhsPat, rhs.subData(iSub));


}
//------------------------------------------------------------------------------

template<int dim>
void Domain::recomputeRHS(VarFcn* vf, DistSVec<double,dim> &V, DistVec<int> &fluidId,
                          DistSVec<double,dim> &rhs, DistExtrapolation<dim>* xpol,
                          DistBcData<dim>& bcData, DistGeoState& geoState, DistSVec<double,3> &X)
{
  int iSub;
#pragma omp parallel for
  for ( iSub = 0; iSub < numLocSub; iSub++){
    Extrapolation<dim>* lxpol = (xpol) ? &((*xpol)(iSub)) : 0;
    subDomain[iSub]->recomputeRHS(vf, V(iSub), fluidId(iSub), rhs(iSub), lxpol,
                                  bcData(iSub), geoState(iSub), X(iSub));
    subDomain[iSub]->sndInletRhsData(*inletRhsPat, rhs.subData(iSub));
  }

  inletRhsPat->exchange();

#pragma omp parallel for
  for ( iSub = 0; iSub < numLocSub; iSub++)
    subDomain[iSub]->addRcvInletRhsData(*inletRhsPat, rhs.subData(iSub));

}

//------------------------------------------------------------------------------

template<int dim>
double Domain::recomputeResidual(DistSVec<double,dim> &F, DistSVec<double,dim> &Finlet)
{

#pragma omp parallel for
  for (int iSub=0; iSub < numLocSub; iSub++)
    subDomain[iSub]->recomputeResidual(F(iSub), Finlet(iSub));

  return Finlet*Finlet;

}

//------------------------------------------------------------------------------

template<int dim>
double Domain::computeRealFluidResidual(DistSVec<double, dim> &F, DistSVec<double,dim> &Freal,
                                        DistLevelSetStructure &dlss)
{
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->computeRealFluidResidual(F(iSub), Freal(iSub), dlss(iSub));

  return Freal*Freal;
}



//------------------------------------------------------------------------------

template<class Scalar, int neq>
void Domain::finishJacobianGalerkinTerm(DistVec<double> &ctrlVol, DistMat<Scalar,neq> &A)  {

  int iSub;

  double t0 = timer->getTime();

  CommPattern<Scalar> *matPat = A.getDiagMatPat();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->finishJacobianGalerkinTerm(ctrlVol(iSub), A(iSub));
    subDomain[iSub]->sndDiagBlocks(*matPat, A(iSub));
  }

  double t = timer->addFiniteVolumeJacTime(t0);

  matPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvDiagBlocks(*matPat, A(iSub));

  com->printf(6, "FV Jacobian matrix computation: %f s\n", t);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeGalerkinTerm(FemEquationTerm *fet, DistBcData<dim> &bcData,
				 DistGeoState &geoState, DistSVec<double,3> &X,
				 DistSVec<double,dim> &V, DistSVec<double,dim> &R,
											DistVec<GhostPoint<dim>*> *ghostPoints,
											DistLevelSetStructure *LSS, 
											bool externalSI)
{

  double t0 = timer->getTime();

  if(ghostPoints)
  {
    if(!LSS) 
    {
			std::cout<<"LSS has to be provided in the case of a viscous simulation\n";
      exit(1);
    }

#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
    {
      subDomain[iSub]->computeGalerkinTerm(fet, bcData(iSub), geoState(iSub),
															 X(iSub), V(iSub), R(iSub),
															 ghostPoints->operator[](iSub),
															 &(LSS->operator()(iSub)), 
															 externalSI);
    }
  }
  else
  {
#pragma omp parallel for
		for (int iSub = 0; iSub < numLocSub; ++iSub) 
		{
      subDomain[iSub]->computeGalerkinTerm(fet, bcData(iSub), geoState(iSub),
                                           X(iSub), V(iSub), R(iSub));
    }
  }
  timer->addFiniteElementTermTime(t0);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Domain::computeDerivativeOfGalerkinTerm(FemEquationTerm *fet, DistBcData<dim> &bcData,
				 DistGeoState &geoState, DistSVec<double,3> &X, DistSVec<double,3> &dX,
				 DistSVec<double,dim> &V, DistSVec<double,dim> &dV, double dMach, DistSVec<double,dim> &dR)
{

  double t0 = timer->getTime();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeDerivativeOfGalerkinTerm(fet, bcData(iSub), geoState(iSub),
					 X(iSub), dX(iSub), V(iSub), dV(iSub), dMach, dR(iSub));

  timer->addFiniteElementTermTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeDerivativeOfGalerkinTerm(dRdXoperators<dim> &dRdXop, FemEquationTerm *fet, DistBcData<dim> &bcData,
				 DistGeoState &geoState, DistSVec<double,3> &X, DistSVec<double,3> &dX,
				 DistSVec<double,dim> &V, DistSVec<double,dim> &dV, double dMach, DistSVec<double,dim> &dR)
{ // YC

  double t0 = timer->getTime();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeDerivativeOfGalerkinTerm(dRdXop.dViscousFluxdX[iSub], fet, bcData(iSub), geoState(iSub),
					 X(iSub), dX(iSub), V(iSub), dV(iSub), dMach, dR(iSub));

  timer->addFiniteElementTermTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeTransposeDerivativeOfGalerkinTerm(dRdXoperators<dim> &dRdXop, DistSVec<double,dim> &dR, DistSVec<double,3> &dX)
{ // YC

  double t0 = timer->getTime();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeTransposeDerivativeOfGalerkinTerm(dRdXop.dViscousFluxdX[iSub], dR(iSub), dX(iSub));

  timer->addFiniteElementTermTime(t0);

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void Domain::computeDerivativeOperatorsOfGalerkinTerm(FemEquationTerm *fet, DistBcData<dim> &bcData,
				 DistGeoState &geoState, DistSVec<double,3> &X,
				 DistSVec<double,dim> &V, RectangularSparseMat<double,3,dim> **dViscousFluxdX)
{

  double t0 = timer->getTime();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeDerivativeOperatorsOfGalerkinTerm(fet, bcData(iSub), geoState(iSub),
					 X(iSub), V(iSub), *dViscousFluxdX[iSub]);

  timer->addFiniteElementTermTime(t0);

}

//------------------------------------------------------------------------------


// Included (MB)
template<int dim>
void Domain::computeOnlyGalerkinTerm(FemEquationTerm *fet, DistBcData<dim> &bcData,
				 DistGeoState &geoState, DistSVec<double,3> &X,
				 DistSVec<double,dim> &V, DistSVec<double,dim> &R)
{

  double t0 = timer->getTime();

  int iSub;
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeGalerkinTerm(fet, bcData(iSub), geoState(iSub),
					 X(iSub), V(iSub), R(iSub));

    subDomain[iSub]->sndData(*vecPat, R.subData(iSub));
  }

  timer->addFiniteElementTermTime(t0);

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, R.subData(iSub));

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeVolumicForceTerm(VolumicForceTerm *volForce, DistVec<double> &ctrlVol,
                               DistSVec<double,dim> &V, DistSVec<double,dim> &R)
{

#pragma omp parallel for
  for (int iSub = 0; iSub <numLocSub; iSub++)
    subDomain[iSub]->computeVolumicForceTerm(volForce, ctrlVol(iSub), V(iSub), R(iSub));

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Domain::computeDerivativeOfVolumicForceTerm(VolumicForceTerm *volForce, DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol,
                               DistSVec<double,dim> &V, DistSVec<double,dim> &dV, DistSVec<double,dim> &dR)
{

#pragma omp parallel for
  for (int iSub = 0; iSub <numLocSub; iSub++)
    subDomain[iSub]->computeDerivativeOfVolumicForceTerm(volForce, ctrlVol(iSub), dCtrlVol(iSub), V(iSub), dV(iSub), dR(iSub));

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeGalerkinBarTerm(bool doInitialTasks,
                                    FemEquationTerm *fet, DistBcData<dim> &bcData,
                                    DistGeoState &geoState, DistSVec<double,3> &X,
                                    DistMacroCellSet *macroCells,
                                    DistSVec<double,dim> &VBar,
                                    DistSVec<double,1> &volRatio,
                                    DistSVec<double,dim> &RBar,
                                    int scopeDepth1, int scopeDepth2)
{

  double t0 = timer->getTime();

  DistSVec<double,dim>* Sigma = new DistSVec<double,dim>(getNodeDistInfo());
  *Sigma = 0.0;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeGalerkinBar_Step1(fet, bcData(iSub), geoState(iSub), X(iSub), VBar(iSub), (*Sigma)(iSub));
    subDomain[iSub]->sndData(*vecPat, (*Sigma).subData(iSub));
  }

  vecPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vecPat, (*Sigma).subData(iSub));


#pragma omp parallel for
  for (int iSub=0; iSub < numLocSub; ++iSub) {
    MacroCellSet** macCells = new MacroCellSet*[scopeDepth2];
    for (int i = 0; i<scopeDepth2; ++i)
        macCells[i] =   macroCells->obtainMacroCell(iSub, i);

    subDomain[iSub]->computeGalerkinBar_Step2(macCells, volRatio(iSub), (*Sigma)(iSub), RBar(iSub), scopeDepth1);
    delete [] macCells;
  }

  delete (Sigma);

  timer->addFiniteElementTermTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computePointWiseSourceTerm(DistGeoState &geoState, DistVec<double> &ctrlVol,
					DistNodalGrad<dim> &ngrad, DistSVec<double,dim> &VV,
					DistSVec<double,dim> &RR)
{

  const double sixth = 1.0/6.;
  const double cb1 = 0.1355;
  const double cb2 = 0.622;
  const double cw2 = 0.3;
  const double cw3 = 2.0;
  const double cv1 = 7.1;
  const double cv2 = 5.0;
  const double sigma = 2.0/3.0;
  const double vkcst = 0.41;
  const double reynolds = 2.91e6;

  const double cw3_pow6 = cw3*cw3*cw3*cw3*cw3*cw3;
  const double opcw3_pow = pow(1.0 + cw3_pow6, 1.0/6.0);
  const double cv1_pow3 = cv1*cv1*cv1;
  const double oocv2 = 1.0 / cv2;
  double oosigma = 1.0 / sigma;
  const double oovkcst2 = 1.0 / (vkcst*vkcst);
  double cw1 = cb1*oovkcst2 + (1.0+cb2) * oosigma;
  const double ooreynolds = 1.0/ reynolds;

  cw1 /= reynolds;
  oosigma /= reynolds;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    Vec<double>& d2w = geoState(iSub).getDistanceToWall();
    double* cvol = ctrlVol.subData(iSub);
    SVec<double,dim>& dVdx = ngrad(iSub).getX();
    SVec<double,dim>& dVdy = ngrad(iSub).getY();
    SVec<double,dim>& dVdz = ngrad(iSub).getZ();
    double (*V)[dim] = VV.subData(iSub);
    double (*R)[dim] = RR.subData(iSub);
    for (int i=0; i<RR.subSize(iSub); ++i) {
      if (d2w[i] < 1.e-10) continue;
      double mutilde = V[i][0]*V[i][5];
      double mul = 1.0;
      double chi = max(mutilde/mul, 0.001);
      double chi3 = chi*chi*chi;
      double fv1 = chi3 / (chi3 + cv1_pow3);
      double fv2 = 1.0 + oocv2*chi;
      fv2 = 1.0 / (fv2*fv2*fv2);
      double fv3 = (1.0 + chi*fv1) * (1.0 - fv2) / chi;
      double d2wall = d2w[i];
      double ood2wall2 = 1.0 / (d2wall * d2wall);
      double rho = V[i][0];
      double oorho = 1.0 / rho;
      double zz = ooreynolds * oovkcst2 * mutilde * oorho * ood2wall2;
      double s12 = dVdy[i][1] - dVdx[i][2];
      double s23 = dVdz[i][2] - dVdy[i][3];
      double s31 = dVdx[i][3] - dVdz[i][1];
      double s = sqrt(s12*s12 + s23*s23 + s31*s31);
      double Stilde = s*fv3 + zz*fv2;
      double rr = min(zz/Stilde, 2.0);
      double rr2 = rr*rr;
      double gg = rr + cw2 * (rr2*rr2*rr2 - rr);
      double gg2 = gg*gg;
      double fw = opcw3_pow * gg * pow(1.0/(gg2*gg2*gg2 + cw3_pow6), sixth);

      double S = cb1 * Stilde * mutilde;
      //S += oosigma * cb2 * rho * (dVdx[i][5]*dVdx[i][5] + dVdy[i][5]*dVdy[i][5] + dVdz[i][5]*dVdz[i][5]);
      S -= cw1 * fw * oorho * mutilde*mutilde * ood2wall2;
      R[i][5] -= cvol[i] * S;
    }
  }

}

//------------------------------------------------------------------------------
//----------------- All LES Models start  here

template<int dim>
void Domain::computeSmagorinskyLESTerm(SmagorinskyLESTerm *smag, DistSVec<double,3> &X,
				       DistSVec<double,dim> &V, DistSVec<double,dim> &R,
													DistVec<GhostPoint<dim>*> *ghostPoints, 
													DistLevelSetStructure *LSS, bool externalSI)

{

  if (ghostPoints)
  {
    if (!LSS) 
    {
      com->fprintf(stderr, "***** Domain::LSS has to be provided in the case of an LES simulation\n");
      exit(1);
    }

#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->computeSmagorinskyLESTerm(smag, X(iSub), V(iSub), R(iSub), 
																	 ghostPoints->operator[](iSub), 
																	 &(LSS->operator()(iSub)), externalSI);
	}
	else 
	{
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->computeSmagorinskyLESTerm(smag, X(iSub), V(iSub), R(iSub));
  }

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeDynamicLESTerm(DynamicLESTerm *dles, DistSVec<double,2> &Cs,
                                   DistSVec<double,3> &X, DistSVec<double,dim> &V,
				   DistSVec<double,dim> &R,
											  DistVec<GhostPoint<dim>*> *ghostPoints, 
											  DistLevelSetStructure *LSS,
											  bool externalSI)
{

  if (ghostPoints)
  {
    if (!LSS) 
    {
      com->fprintf(stderr, "***** Domain::LSS has to be provided in the case of an LES simulation\n");
      exit(1);
    }

#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->computeDynamicLESTerm(dles, Cs(iSub), X(iSub), V(iSub), R(iSub),
															ghostPoints->operator[](iSub), &(LSS->operator()(iSub)), 
															externalSI);
  }
  else 
  {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->computeDynamicLESTerm(dles, Cs(iSub), X(iSub), V(iSub), R(iSub));
  }

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeVMSLESTerm(VMSLESTerm *vmst, DistMacroCellSet *macroCells,
                               bool doInitialTasks, DistVec<double> &ctrlVol,
                               DistSVec<double,dim> &VBar, DistSVec<double,1> &volRatio,
                               DistSVec<double,3> &X, DistSVec<double,dim> &V,
                               DistSVec<double,dim> &R, int scopeDepth)
{

  double t0 = timer->getTime();

  DistSVec<double,dim>* Sigma = new DistSVec<double,dim>(getNodeDistInfo());
  *Sigma = 0.0;

  // Compute the large-scale component of V (VBar) and //
  // the volume ratios (Vi/VIi) with the macro-cells   //

  macroCells->computeVMS(doInitialTasks, ctrlVol, X, V, VBar, volRatio, scopeDepth);

  // Exchange VBar and volRatio values //

  assemble(vecPat, VBar);
  assemble(vecPat, volRatio);

  // Compute the contribution to R from the tetrahedra contained within each subdomain //
#pragma omp parallel for
  for (int iSub=0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeVMSLES_Step1(vmst, VBar(iSub),X(iSub), V(iSub), (*Sigma)(iSub));

  assemble(vecPat,*Sigma);

  // Add the subgrid scale viscosity to the residual vector R //

#pragma omp parallel for
  for (int iSub=0; iSub < numLocSub; ++iSub) {
    MacroCellSet* macCells;
    macCells = macroCells->obtainMacroCell(iSub, scopeDepth-1);
    subDomain[iSub]->computeVMSLES_Step2(volRatio(iSub),macCells,(*Sigma)(iSub), R(iSub), scopeDepth);
  }

  delete (Sigma);

  double t = timer->addVMSLESTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeDynamicVMSTerm(DynamicVMSTerm *dvmst, DistMacroCellSet *macroCells,
                                   bool doInitialTasks, DistVec<double> &ctrlVol, DistSVec<double,dim> **VBar,
                                   DistSVec<double,1> **volRatio, DistSVec<double,3> &X, DistSVec<double,dim> &V,
                                   DistSVec<double,dim> &S, DistSVec<double,dim> &R, DistSVec<double,dim> &RBar,
                                   DistSVec<double,dim> &dWdt, DistSVec<double,dim> &dWBardt,
                                   DistSVec<double,dim> &M, DistSVec<double,dim> &MBar,
                                   int scopeDepth1, int scopeDepth2, int method, DistVec<double> *Cs)
{

  double t0 = timer->getTime();

  DistSVec<double,dim>* Sigma = new DistSVec<double,dim>(getNodeDistInfo());
  DistSVec<double,dim>* SigmaBar = new DistSVec<double,dim>(getNodeDistInfo());

  *Sigma = 0.0;
  *SigmaBar = 0.0;

  // compute volume averaged filter width for each node (used in post processing) //

  if (doInitialTasks) {
    Delta = new DistVec<double>(getNodeDistInfo());
    CsDelSq = new DistVec<double>(getNodeDistInfo());
    PrT = new DistVec<double>(getNodeDistInfo());
    WCsDelSq = new DistVec<double>(getNodeDistInfo());
    WPrT = new DistVec<double>(getNodeDistInfo());

   *Delta = 0.0;

    #pragma omp parallel for
      for (int iSub=0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->computeFilterWidth(X(iSub),(*Delta)(iSub));

    applySmoothing(ctrlVol, *Delta);
  }


  *WCsDelSq = 0.0; *WPrT = 0.0;  // variables that store the smoothed CsDelSq and PrT values


  // compute MBar and M //

#pragma omp parallel for
    for (int iSub=0; iSub < numLocSub; ++iSub) {
      SVec<double,dim> **vBar = new SVec<double,dim>*[2];
      SVec<double,1> **volR = new SVec<double,1>*[2];

      for (int i = 0; i<2; ++i) {
        vBar[i]  =  &(*VBar[i])(iSub);
        volR[i]  =  &(*volRatio[i])(iSub);
      }

      subDomain[iSub]->computeMBarAndM_Step1(dvmst, vBar, volR, X(iSub), V(iSub), (*SigmaBar)(iSub), (*Sigma)(iSub));

      delete [] vBar;
      delete [] volR;
    }


    assemble(vecPat, (*SigmaBar));
    assemble(vecPat, (*Sigma));

#pragma omp parallel for
    for (int iSub=0; iSub < numLocSub; ++iSub) {
      SVec<double,1> **volR = new SVec<double,1>*[2];
      MacroCellSet** macCells = new MacroCellSet*[scopeDepth2];

      for (int i = 0; i<2; ++i) {
        volR[i]  =  &(*volRatio[i])(iSub);
      }

      for (int i = 0; i<scopeDepth2; ++i)
        macCells[i] =   macroCells->obtainMacroCell(iSub, i);

      subDomain[iSub]->computeMBarAndM_Step2(macCells, volR, MBar(iSub), M(iSub), (*SigmaBar)(iSub),
                                             (*Sigma)(iSub), scopeDepth1, scopeDepth2);

      delete [] volR;
      delete [] macCells;
    }


    assemble(vecPat, MBar);
    assemble(vecPat, M);

  // compute (Cs*Delta)^2 and PrT //

#pragma omp parallel for
   for (int iSub=0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->computeCsDeltaSq(R(iSub), RBar(iSub), M(iSub), MBar(iSub),
                                        dWdt(iSub), dWBardt(iSub), (*CsDelSq)(iSub),
                                        (*PrT)(iSub), method);
   }

  // local smoothing of (Cs*Delta)^2 and PrT //

#pragma omp parallel for
   for (int iSub=0; iSub < numLocSub; ++iSub) {
       subDomain[iSub]->computeLocalAvg(X(iSub), (*CsDelSq)(iSub), (*WCsDelSq)(iSub));
       subDomain[iSub]->computeLocalAvg(X(iSub), (*PrT)(iSub), (*WPrT)(iSub));
   }

  applySmoothing(ctrlVol, *WCsDelSq);
  applySmoothing(ctrlVol, *WPrT);


  *Sigma = 0.0;

  // computing the reynold stress flux "S" //

#pragma omp parallel for
    for (int iSub=0; iSub < numLocSub; ++iSub) {
      SVec<double,dim> **vBar = new SVec<double,dim>*[2];
      Vec<double>* cs;

      for (int i = 0; i<2; ++i) {
        vBar[i]  =  &(*VBar[i])(iSub);
      }

      if (Cs) cs = &(*Cs)(iSub);
      else cs = 0;

      subDomain[iSub]->computeDynamicVMSTerm_Step1(dvmst, vBar, X(iSub), V(iSub), (*Sigma)(iSub),
                                                   (*WCsDelSq)(iSub), (*WPrT)(iSub), cs, (*Delta)(iSub));

      delete [] vBar;
    }
    assemble(vecPat, *Sigma);

#pragma omp parallel for
    for (int iSub=0; iSub < numLocSub; ++iSub) {
      SVec<double,1> **volR = new SVec<double,1>*[2];
      MacroCellSet** macCells = new MacroCellSet*[scopeDepth2];
      Vec<double>* cs;

      for (int i = 0; i<2; ++i) {
        volR[i]  =  &(*volRatio[i])(iSub);
      }

      for (int i = 0; i<scopeDepth2; ++i)
        macCells[i] =   macroCells->obtainMacroCell(iSub, i);

      if (Cs) cs = &(*Cs)(iSub);
      else cs = 0;

      subDomain[iSub]->computeDynamicVMSTerm_Step2(macCells, volR, (*Sigma)(iSub), S(iSub), scopeDepth1);

      delete [] volR;
      delete [] macCells;
    }

    assemble(vecPat, S);

  delete (Sigma);
  delete (SigmaBar);

  double t = timer->addDynamicVMSLESTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeWaleLESTerm(WaleLESTerm *wale, DistSVec<double,3> &X,
				DistSVec<double,dim> &V, DistSVec<double,dim> &R,
										  DistVec<GhostPoint<dim>*> *ghostPoints, 
										  DistLevelSetStructure *LSS, bool externalSI)

{

  if (ghostPoints)
  {
    if (!LSS) 
    {
      com->fprintf(stderr, "***** Domain::LSS has to be provided in the case of an LES simulation\n");
      exit(1);
    }

#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->computeWaleLESTerm(wale, X(iSub), V(iSub), R(iSub),
															ghostPoints->operator[](iSub), &(LSS->operator()(iSub)), externalSI);
  }
  else {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->computeWaleLESTerm(wale, X(iSub), V(iSub), R(iSub));
  }

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeVBar(DistMacroCellSet *macroCells, bool doInitialTasks, DistGeoState &geoState,
                         DistSVec<double,dim> &VBar, DistSVec<double,dim> &V, int scopeDepth, int n)

{

  double t0 = timer->getTime();

  macroCells->computeVBar(doInitialTasks, geoState, V, VBar, scopeDepth, n);

  #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->sndData(*vecPat, VBar.subData(iSub));

    vecPat->exchange();

  #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvData(*vecPat, VBar.subData(iSub));

  double t = timer->addDynamicVMSLESTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeBarTerm(DistMacroCellSet *macroCells, bool doInitialTasks, DistVec<double> &ctrlVol,
                            DistSVec<double,dim> **VBar, DistSVec<double,1> **volRatio,
                            DistSVec<double,3> &X, DistSVec<double,dim> &V, int scopeDepth1, int scopeDepth2)

{

  double t0 = timer->getTime();

  // Computing the large-scale components VBar and volRatio //

  macroCells->computeDVMS(doInitialTasks, ctrlVol, X, V, VBar, volRatio, scopeDepth1, scopeDepth2);

   // Exchange of information between subdomains for VBar and volRatio //

   for (int i =0; i < 2; ++i) {

    #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->sndData(*vecPat, (*VBar[i]).subData(iSub));

    vecPat->exchange();

    #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvData(*vecPat, (*VBar[i]).subData(iSub));

    #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
     subDomain[iSub]->sndData(*vecPat, (*volRatio[i]).subData(iSub));

    vecPat->exchange();

    #pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
     subDomain[iSub]->addRcvData(*vecPat, (*volRatio[i]).subData(iSub));

   }

  double t = timer->addDynamicVMSLESTime(t0);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeTestFilterValues(DistVec<double> &ctrlVol,
                                     DistSVec<double,dim> &VCap,
                                     DistSVec<double,16> &Mom_Test,
                                     DistSVec<double,6> &Sij_Test,
                                     DistVec<double> &modS_Test,
                                     DistSVec<double,8> &Eng_Test,
                                     DistSVec<double,2> &Cs,
				     DistVec<int> &Ni,
                                     DistBcData<dim> &bcData,
												 DistSVec<double,3> &X, 
												 DistSVec<double,dim> &V,
				     double gam, double R,
				     DistVec<GhostPoint<dim>*> *ghostPoints, 
                                     DistLevelSetStructure *LSS, bool externalSI)
{

// computing test filtered values for all the required flow variables //

  if (ghostPoints)
  {
    if (!LSS) 
    {
      com->fprintf(stderr, "***** Domain::LSS has to be provided in the case of an LES simulation\n");
      exit(1);
    }

#pragma omp parallel for
		for(int iSub = 0; iSub < numLocSub; ++iSub)
		{
      subDomain[iSub]->computeTestFilterAvgs(VCap(iSub), Mom_Test(iSub), Sij_Test(iSub), modS_Test(iSub),
                                             Eng_Test(iSub), X(iSub), V(iSub), gam, R,
																ghostPoints->operator[](iSub), &(LSS->operator()(iSub)), externalSI);
      subDomain[iSub]->sndData(*vecPat, VCap.subData(iSub));
    }
  }
  else
  {
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub){
      subDomain[iSub]->computeTestFilterAvgs(VCap(iSub), Mom_Test(iSub), Sij_Test(iSub), modS_Test(iSub),
                                             Eng_Test(iSub), X(iSub), V(iSub), gam, R);
      subDomain[iSub]->sndData(*vecPat, VCap.subData(iSub));
    }
  }

   // START OF ALL EXCHANGES

  vecPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->addRcvData(*vecPat, VCap.subData(iSub));
    subDomain[iSub]->sndData(*momPat, Mom_Test.subData(iSub));
  }

  momPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->addRcvData(*momPat, Mom_Test.subData(iSub));
    subDomain[iSub]->sndData(*engPat, Eng_Test.subData(iSub));
  }

  engPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->addRcvData(*engPat, Eng_Test.subData(iSub));
    subDomain[iSub]->sndData(*weightPat, Sij_Test.subData(iSub));
  }

  weightPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->addRcvData(*weightPat, Sij_Test.subData(iSub));
    subDomain[iSub]->sndData(*volPat, reinterpret_cast<double (*)[1]>(modS_Test.subData(iSub)));
  }

  volPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*volPat, reinterpret_cast<double (*)[1]>(modS_Test.subData(iSub)));

   // END OF ALL EXCHANGES


// computing unsmoothed Cs and Pt values //

  if (ghostPoints) {
#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->computeCsValues(VCap(iSub), Mom_Test(iSub), Sij_Test(iSub),
                                       modS_Test(iSub), Eng_Test(iSub), Cs(iSub),
         			       Ni(iSub), X(iSub), gam, R,
                                       &(LSS->operator()(iSub)));
    }
  }
  else {
#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub) {
      subDomain[iSub]->computeCsValues(VCap(iSub), Mom_Test(iSub), Sij_Test(iSub),
                                       modS_Test(iSub), Eng_Test(iSub), Cs(iSub),
         			       Ni(iSub), X(iSub), gam, R);
    }
  }


// smoothing Cs and Pt values //

   DistSVec<double,2>* W = new DistSVec<double,2>(getNodeDistInfo());
   *W = 0.0;

#pragma omp parallel for
   for (int iSub=0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->computeLocalAvg(X(iSub), Cs(iSub), (*W)(iSub));

   applySmoothing(ctrlVol, *W);
   Cs = *W;

   delete (W);

}

//------------------------------------------------------------------------------
// Included (MB)
template<int dim>
void Domain::computeDerivativeOfSmagorinskyLESTerm(SmagorinskyLESTerm *smag, DistSVec<double,3> &X,
				       DistSVec<double,dim> &V, DistSVec<double,dim> &R)

{

  com->fprintf(stderr, "***** Domain::computeDerivativeOfSmagorinskyLESTerm is not implemented!\n");
  exit(1);

}

//----------------------- All LES Models End Here
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//--------Start of routines that compute MutOMu values

template<int dim>
void Domain::computeMutOMuSmag(SmagorinskyLESTerm *smag, DistVec<double> &ctrlVol,
                               DistSVec<double,3> &X, DistSVec<double,dim> &V,
                               DistVec<double> &mutOmu)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeMutOMuSmag(smag, X(iSub), V(iSub), mutOmu(iSub));

  applySmoothing(ctrlVol, mutOmu);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeMutOMuVMS(VMSLESTerm *vms, DistMacroCellSet *macroCells, DistVec<double> &ctrlVol,
                              bool doInitialTasks, DistSVec<double,dim> &VBar, DistSVec<double,1> &volRatio,
                              DistSVec<double,3> &X, DistSVec<double,dim> &V, int scopeDepth,
                              DistVec<double> &mutOmu)
{

  macroCells->computeVMS(doInitialTasks, ctrlVol, X, V, VBar, volRatio, scopeDepth);
  assemble(vecPat, VBar);

#pragma omp parallel for
  for (int iSub=0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeMutOMuVMS(vms, VBar(iSub),X(iSub), V(iSub), mutOmu(iSub));

  applySmoothing(ctrlVol, mutOmu);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeMutOMuDynamicVMS(DynamicVMSTerm *dvms, DistVec<double> &ctrlVol,
                                     DistSVec<double,dim> &VBar, DistSVec<double,3> &X,
                                     DistSVec<double,dim> &V, DistVec<double> &Cs, DistVec<double> &mutOmu)
{

#pragma omp parallel for
  for (int iSub=0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeMutOMuDynamicVMS(dvms, VBar(iSub),X(iSub), V(iSub), Cs(iSub), mutOmu(iSub));

  applySmoothing(ctrlVol, mutOmu);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeMutOMuWale(WaleLESTerm *wale, DistVec<double> &ctrlVol,
                               DistSVec<double,3> &X, DistSVec<double,dim> &V,
                               DistVec<double> &mutOmu)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeMutOMuWale(wale, X(iSub), V(iSub), mutOmu(iSub));

  applySmoothing(ctrlVol, mutOmu);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeMutOMuDynamicLES(DynamicLESTerm *dles, DistVec<double> &ctrlVol,
                                     DistSVec<double,2> &Cs, DistSVec<double,3> &X,
                                     DistSVec<double,dim> &V, DistVec<double> &mutOmu)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeMutOMuDynamicLES(dles, Cs(iSub), X(iSub), V(iSub), mutOmu(iSub));

  applySmoothing(ctrlVol, mutOmu);

}


//--------End of routines that compute MutOMu values
//------------------------------------------------------------------------------
template<int dim, class Scalar, int neq>
void Domain::computeJacobianGalerkinTerm(FemEquationTerm *fet, DistBcData<dim> &bcData,
					 DistGeoState &geoState, DistSVec<double,3> &X,
					 DistVec<double> &ctrlVol, DistSVec<double,dim> &V,
					 DistMat<Scalar,neq> &A,
                                         DistVec<GhostPoint<dim>*> *ghostPoints,
													  DistLevelSetStructure *distLSS, bool externalSI)
{

  int iSub;

  double t0 = timer->getTime();

  CommPattern<Scalar> *matPat = A.getOffDiagMatPat();

#pragma omp parallel for
	for(iSub = 0; iSub < numLocSub; ++iSub) 
	{
    Vec<GhostPoint<dim>*>* gp = (ghostPoints? &(*ghostPoints)(iSub) :0);
    LevelSetStructure *LSS = distLSS ? &(distLSS->operator()(iSub)) : 0;
    subDomain[iSub]->computeJacobianGalerkinTerm(fet, bcData(iSub), geoState(iSub), X(iSub),
																	ctrlVol(iSub), V(iSub), A(iSub), gp, LSS, externalSI);
		
    subDomain[iSub]->sndOffDiagBlocks(*matPat, A(iSub));
  }

  double t = timer->addFiniteElementJacTime(t0);

  matPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvOffDiagBlocks(*matPat, A(iSub));
  
	if(ghostPoints && !externalSI) 
	{
#pragma omp parallel for
		for(iSub = 0; iSub < numLocSub; ++iSub) 
      subDomain[iSub]->sndGhostOffDiagBlocks(*matPat, A(iSub));
    
    matPat->exchange();

#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->addRcvGhostOffDiagBlocks(*matPat, A(iSub));
  }

  com->printf(6, "FE Jacobian matrix computation: %f s\n", t);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Domain::computeBCsJacobianWallValues(FemEquationTerm *fet, DistBcData<dim> &bcData,
					  DistGeoState &geoState, DistSVec<double,3> &X,
					  DistSVec<double,dim> &V)
{

  int iSub;

  double t0 = timer->getTime();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeBCsJacobianWallValues(fet, bcData(iSub), geoState(iSub), X(iSub), V(iSub));
  }

  double t = timer->addFiniteElementJacTime(t0);

  com->printf(6, "FE wall BC computation: %f s\n", t);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void Domain::computeJacobianVolumicForceTerm(VolumicForceTerm *volForce,
                               DistVec<double> &ctrlVol,
                               DistSVec<double,dim> &V, DistMat<Scalar,neq> &A)
{

#pragma omp parallel for
  for (int iSub = 0; iSub <numLocSub; iSub++)
    subDomain[iSub]->computeJacobianVolumicForceTerm(volForce,
                                           ctrlVol(iSub), V(iSub), A(iSub));

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::getExtrapolationValue(DistExtrapolation<dim> *xpol, DistSVec<double,dim> &V,
                                      DistSVec<double,dim> &Ubc, VarFcn* vf,
                                 DistBcData<dim>& bcData, DistGeoState& geoState, DistSVec<double,3>& X)
{
  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub){
    Extrapolation<dim>* lxpol = (xpol) ? &((*xpol)(iSub)) : 0;
    subDomain[iSub]->getExtrapolationValue(lxpol, V(iSub), Ubc(iSub), vf, bcData(iSub), geoState(iSub), X(iSub));
    subDomain[iSub]->sndInletData(*inletRhsPat, Ubc.subData(iSub));
  }
 inletRhsPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub){
   subDomain[iSub]->addRcvInletData(*inletRhsPat, Ubc.subData(iSub), true);
  }


}
//------------------------------------------------------------------------------

template<int dim>
void Domain::applyExtrapolationToSolutionVector(DistExtrapolation<dim> *xpol, DistSVec<double,dim> &U,
                                      DistSVec<double,dim> &Ubc)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub){
    Extrapolation<dim>* lxpol = (xpol) ? &((*xpol)(iSub)) : 0;
    subDomain[iSub]->applyExtrapolationToSolutionVector(lxpol, U(iSub), Ubc(iSub));
  }

}
//------------------------------------------------------------------------------

template<int dim>
void Domain::applyBCsToSolutionVector(BcFcn *bcFcn, DistBcData<dim> &bcData,
                                      DistSVec<double,dim> &U, DistLevelSetStructure *distLSS)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    { 
      LevelSetStructure *LSS = distLSS ? &((*distLSS)(iSub)) : 0;
      subDomain[iSub]->applyBCsToSolutionVector(bcFcn, bcData(iSub), U(iSub), LSS);
    }
}

//------------------------------------------------------------------------------

template<int dim>
void Domain::applyBCsToTurbSolutionVector(BcFcn *bcFcn, DistBcData<dim> &bcData,
                                      DistSVec<double,dim> &U, DistLevelSetStructure *distLSS)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    { 
      LevelSetStructure *LSS = distLSS ? &((*distLSS)(iSub)) : 0;
      subDomain[iSub]->applyBCsToTurbSolutionVector(bcFcn, bcData(iSub), U(iSub), LSS);
    }
}

//------------------------------------------------------------------------------

template<int dim>
void Domain::applyBCsToResidual(BcFcn *bcFcn, DistBcData<dim> &bcData,
										  DistSVec<double,dim> &U, DistSVec<double,dim> &F, 
										  DistLevelSetStructure *distLSS)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    {
      LevelSetStructure *LSS = distLSS ? &((*distLSS)(iSub)) : 0;
      subDomain[iSub]->applyBCsToResidual(bcFcn, bcData(iSub), U(iSub), F(iSub), LSS);
    }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Domain::applyBCsToDerivativeOfResidual(BcFcn *bcFcn, DistBcData<dim> &bcData,
				DistSVec<double,dim> &U, DistSVec<double,dim> &dU, DistSVec<double,dim> &dF)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->applyBCsToDerivativeOfResidual(bcFcn, bcData(iSub), U(iSub), dU(iSub), dF(iSub));

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void Domain::applyBCsToJacobian(BcFcn *bcFcn, DistBcData<dim> &bcData,
				DistSVec<double,dim> &U, DistMat<Scalar,neq> &A, DistLevelSetStructure *distLSS)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
  {
    LevelSetStructure *LSS = distLSS ? &((*distLSS)(iSub)) : 0;
    subDomain[iSub]->applyBCsToJacobian(bcFcn, bcData(iSub), U(iSub), A(iSub), LSS);
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void Domain::applyBCsToH2Jacobian(BcFcn *bcFcn, DistBcData<dim> &bcData,
	                          DistSVec<double,dim> &U, DistMat<Scalar,neq> &A)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->applyBCsToH2Jacobian(bcFcn, bcData(iSub), U(iSub), A(iSub));

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar, int neq>
void Domain::applyBCsToJacobianWallValues(BcFcn *bcFcn, DistBcData<dim> &bcData,
				DistSVec<double,dim> &U, DistMat<Scalar,neq> &A)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->applyBCsToJacobianWallValues(bcFcn, bcData(iSub), U(iSub), A(iSub));

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar2>
void Domain::applyBCsToProduct(BcFcn *bcFcn, DistBcData<dim> &bcData, DistSVec<double,dim> &U, DistSVec<Scalar2,dim> &Prod)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->applyBCsToProduct(bcFcn, bcData(iSub), U(iSub), Prod(iSub));

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void Domain::computeH1(FluxFcn **fluxFcn, DistBcData<dim> &bcData,
                       DistGeoState &geoState, DistVec<double> &ctrlVol,
                       DistSVec<double,dim> &V, DistMat<Scalar,dim> &H1)
{

  int iSub;

  double t0 = timer->getTime();

  CommPattern<Scalar> *matPat = H1.getDiagMatPat();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeH1(fluxFcn, bcData(iSub), geoState(iSub),
                               ctrlVol(iSub), V(iSub), H1(iSub));
    subDomain[iSub]->sndDiagBlocks(*matPat, H1(iSub));
  }

  //double t = timer->addH1SetupTime(t0);

  matPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvDiagBlocks(*matPat, H1(iSub));

  //com->printf("Time for computing H1 matrix: %f s\n", t);

}

//------------------------------------------------------------------------------


template<int dim, class Scalar, int neq>
void Domain::computeH2transpose(FluxFcn **fluxFcn, RecFcn *recFcn,
           DistBcData<dim> &bcData, DistGeoState &geoState,
           DistSVec<double,3> &X, DistSVec<double,dim> &V,
           DistNodalGrad<dim, double> &ngrad, DistMat<Scalar,neq> &H2transpose,
           DistSVec<double,dim> &aij, DistSVec<double,dim> &aji,
           DistSVec<double,dim> &bij, DistSVec<double,dim> &bji)
{

  double t0 = timer->getTime();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeH2transpose(fluxFcn, recFcn, bcData(iSub), geoState(iSub),
             X(iSub), V(iSub), ngrad(iSub), H2transpose(iSub));
    subDomain[iSub]->precomputeRec(recFcn, X(iSub), V(iSub), ngrad(iSub),
           aij(iSub), aji(iSub), bij(iSub), bji(iSub));
  }

  double t = timer->addH2SetupTime(t0);

  com->printf(6, "H2 matrix transpose computation: %f s\n", t);

}

//------------------------------------------------------------------------------

template<class Scalar1, class Scalar2, int dim>
void Domain::computeMatVecProdH2transposeNew(IoData& iod, DistSVec<double,3> &X,
                                             DistVec<double> &ctrlVol, DistMat<Scalar1,dim> &H2,
                                             DistSVec<double,dim> &aij, DistSVec<double,dim> &aji,
                                             DistSVec<double,dim> &bij, DistSVec<double,dim> &bji,
                                             DistNodalGrad<dim, Scalar2> &dpdxj,
                                             DistSVec<Scalar2,dim> &p, DistSVec<Scalar2,dim> &prod)  {

  int iSub;

  CommPattern<Scalar2> *vPat = getCommPat(p);
  DistSVec<Scalar2,dim> ddxt(p.info());
  DistSVec<Scalar2,dim> ddyt(p.info());
  DistSVec<Scalar2,dim> ddzt(p.info());
  DistSVec<Scalar2,dim> cij(aij.info()), cji(aij.info()), dij(aij.info()), dji(aij.info());

  if(iod.schemes.ns.reconstruction == SchemeData::LINEAR) {
#pragma omp parallel for
    for (iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->computeGradientsTransposeNew(X(iSub), ctrlVol(iSub), H2(iSub),
                                                    bij(iSub), bji(iSub), cij(iSub), cji(iSub), dij(iSub), dji(iSub),
                                                    p(iSub), ddxt(iSub), ddyt(iSub), ddzt(iSub));

    assemble(vPat, ddxt);
    assemble(vPat, ddyt);
    assemble(vPat, ddzt);
  }

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeMatVecProdH2transposeNew(iod, X(iSub), ctrlVol(iSub), H2(iSub),
                                                     aij(iSub), aji(iSub), bij(iSub), bji(iSub),
                                                     cij(iSub), cji(iSub), dij(iSub), dji(iSub),
                                                     dpdxj(iSub), p(iSub), prod(iSub));
  if(iod.schemes.ns.reconstruction == SchemeData::LINEAR) {
    if(iod.schemes.ns.gradient == SchemeData::LEAST_SQUARES) {
      DistSVec<double,6> R = dpdxj.getR();
#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->addToMatVecProdH2transposeLeastSquareNew(X(iSub), R(iSub), ddxt(iSub), ddyt(iSub), ddzt(iSub),
                                                                  dpdxj(iSub), prod(iSub));
    }
    else if(iod.schemes.ns.gradient == SchemeData::GALERKIN || iod.schemes.ns.gradient == SchemeData::NON_NODAL) {
#pragma omp parallel for
      for (iSub = 0; iSub < numLocSub; ++iSub)
        subDomain[iSub]->addToMatVecProdH2transposeGalerkinNew(ctrlVol(iSub), ddxt(iSub), ddyt(iSub), ddzt(iSub),
                                                               dpdxj(iSub), prod(iSub));
    }
  }

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addDiagonalInMatVecProdH2transpose(ctrlVol(iSub), H2(iSub), p(iSub), prod(iSub));

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->sndData(*vPat, prod.subData(iSub));

  vPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vPat, prod.subData(iSub));

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void Domain::computeH2(FluxFcn **fluxFcn, RecFcn *recFcn,
		       DistBcData<dim> &bcData, DistGeoState &geoState,
		       DistSVec<double,3> &X, DistSVec<double,dim> &V,
		       DistNodalGrad<dim, double> &ngrad, DistMat<Scalar,neq> &H2,
		       DistSVec<double,dim> &aij, DistSVec<double,dim> &aji,
		       DistSVec<double,dim> &bij, DistSVec<double,dim> &bji)
{

  double t0 = timer->getTime();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeH2(fluxFcn, recFcn, bcData(iSub), geoState(iSub),
			       X(iSub), V(iSub), ngrad(iSub), H2(iSub));

    subDomain[iSub]->precomputeRec(recFcn, X(iSub), V(iSub), ngrad(iSub),
				   aij(iSub), aji(iSub), bij(iSub), bji(iSub));
  }

  double t = timer->addH2SetupTime(t0);

  com->printf(6, "H2 matrix computation: %f s\n", t);

}
//------------------------------------------------------------------------------


template<int dim, class Scalar, int neq>
void Domain::computeH2(FluxFcn **fluxFcn, RecFcn *recFcn,
		       DistBcData<dim> &bcData, DistGeoState &geoState,
		       DistSVec<double,3> &X, DistSVec<double,dim> &V,
		       DistNodalGrad<dim, double> &ngrad,
		       DistExactRiemannSolver<dim> &riemann,
		       DistLevelSetStructure *distLSS,
		       DistVec<int> &fluidId, 		      
		       int Nriemann,
		       DistMat<Scalar,neq> &H2,
		       DistSVec<double,dim> &aij, DistSVec<double,dim> &aji,
		       DistSVec<double,dim> &bij, DistSVec<double,dim> &bji,
		       DistSVec<double,dim> &betaij, DistSVec<double,dim> &betaji)
{

  double t0 = timer->getTime();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {

    Vec<int> &FluidId = fluidId(iSub);

    subDomain[iSub]->computeH2(fluxFcn, recFcn, bcData(iSub), geoState(iSub),
			       X(iSub), V(iSub), ngrad(iSub), riemann(iSub),
			       (*distLSS)(iSub), FluidId, Nriemann,
			       H2(iSub), aij(iSub), aji(iSub), bij(iSub), bji(iSub), 
			       betaij(iSub), betaji(iSub));

    subDomain[iSub]->precomputeRec(recFcn, X(iSub), V(iSub), ngrad(iSub), (*distLSS)(iSub), FluidId,
    				   aij(iSub), aji(iSub), bij(iSub), bji(iSub));

  }

  double t = timer->addH2SetupTime(t0);

  com->printf(6, "H2 matrix computation: %f s\n", t);

}


//------------------------------------------------------------------------------

template<class Scalar1, class Scalar2, int dim>
void Domain::computeMatVecProdH2(RecFcn *recFcn, DistSVec<double,3> &X,
				 DistVec<double> &ctrlVol, DistMat<Scalar1,dim> &H2,
				 DistSVec<double,dim> &aij, DistSVec<double,dim> &aji,
				 DistSVec<double,dim> &bij, DistSVec<double,dim> &bji,
				 DistSVec<Scalar2,dim> &p, DistNodalGrad<dim, Scalar2> &dpdxj,
				 DistSVec<Scalar2,dim> &prod)  {

  int iSub;

  CommPattern<Scalar2> *vPat = getCommPat(p);

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeMatVecProdH2(recFcn, X(iSub), ctrlVol(iSub), H2(iSub),
					 aij(iSub), aji(iSub), bij(iSub), bji(iSub),
					 p(iSub), dpdxj(iSub), prod(iSub));
    subDomain[iSub]->sndData(*vPat, prod.subData(iSub));
  }

  vPat->exchange();
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vPat, prod.subData(iSub));

}

//------------------------------------------------------------------------------

template<class Scalar1, class Scalar2, int dim>
void Domain::computeMatVecProdH2(FluxFcn **fluxFcn, RecFcn *recFcn, DistGeoState &geoState,
				 DistSVec<double,3> &X, DistVec<double> &ctrlVol, 
				 DistExactRiemannSolver<dim> &riemann,
				 DistLevelSetStructure *distLSS,
				 DistVec<int> &fluidId,
				 int Nriemann,
				 DistMat<Scalar1,dim> &H2,
				 DistSVec<double,dim> &aij, DistSVec<double,dim> &aji,
				 DistSVec<double,dim> &bij, DistSVec<double,dim> &bji,
				 DistSVec<double,dim> &betaij, DistSVec<double,dim> &betaji,
				 DistSVec<Scalar2,dim> &p, DistNodalGrad<dim, Scalar2> &dpdxj,
				 DistSVec<Scalar2,dim> &prod)  {

  int iSub;

  CommPattern<Scalar2> *vPat = getCommPat(p);

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {

        Vec<int> &FluidId = fluidId(iSub);

        subDomain[iSub]->computeMatVecProdH2(fluxFcn, recFcn, geoState(iSub), 
					     X(iSub), ctrlVol(iSub),
					     riemann(iSub), (*distLSS)(iSub), FluidId,
					     Nriemann, H2(iSub),
					     aij(iSub), aji(iSub), bij(iSub), bji(iSub),
					     betaij(iSub), betaji(iSub),
					     p(iSub), dpdxj(iSub), prod(iSub));
	subDomain[iSub]->sndData(*vPat, prod.subData(iSub));

  }

  
    vPat->exchange();
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vPat, prod.subData(iSub));
  
}

//------------------------------------------------------------------------------

template<class Scalar1, class Scalar2, int dim>
void Domain::computeMatVecProdH2T(RecFcn *recFcn, DistSVec<double,3> &X,
                DistVec<double> &ctrlVol, DistMat<Scalar1,dim> &H2,
                DistSVec<double,dim> &aij, DistSVec<double,dim> &aji,
                DistSVec<double,dim> &bij, DistSVec<double,dim> &bji,
                DistSVec<Scalar2,dim> &p, DistSVec<Scalar2,dim> &prod,
                DistSVec<Scalar2,dim> &prod2, DistSVec<Scalar2,dim> &prod3,
                DistSVec<Scalar2,dim> &prod4)  {

  int iSub;

  CommPattern<Scalar2> *vPat = getCommPat(p);
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeMatVecProdH2T(recFcn, X(iSub), ctrlVol(iSub),
        H2(iSub), aij(iSub), aji(iSub), bij(iSub), bji(iSub), p(iSub),
        prod(iSub), prod2(iSub), prod3(iSub), prod4(iSub));

    subDomain[iSub]->sndData(*vPat, prod.subData(iSub));
  }

  vPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*vPat, prod.subData(iSub));

 assemble(vPat, prod2);
 assemble(vPat, prod3);
 assemble(vPat, prod4);

}

//------------------------------------------------------------------------------

template<class Scalar1, class Scalar2, int dim>
void Domain::computeMatVecProdH2Tb(RecFcn *recFcn, DistSVec<double,3> &X,
                DistVec<double> &ctrlVol, DistMat<Scalar1,dim> &H2,
                DistNodalGrad<dim, Scalar2> &dpdxj, DistSVec<Scalar2,dim> &p,
                DistSVec<Scalar2,dim> &prod, DistSVec<Scalar2,dim> &prod2)
{

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->computeMatVecProdH2Tb(recFcn, X(iSub), ctrlVol(iSub),
        H2(iSub), dpdxj(iSub), p(iSub), prod(iSub), prod2(iSub) );
  }
// No assemble????

}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
void Domain::assemble(CommPattern<Scalar> *commPat, DistSVec<Scalar,dim> &W)
{

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->sndData(*commPat, W.subData(iSub));

  commPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->addRcvData(*commPat, W.subData(iSub));

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, class OpType >
void Domain::assemble(CommPattern<Scalar> *commPat, DistSVec<Scalar,dim> &W, const OpType &oper)
{

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->sndData(*commPat, W.subData(iSub));

  commPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->operateRcvData(*commPat, W.subData(iSub), oper);

}

//------------------------------------------------------------------------------

template<class Scalar>
void Domain::assemble(CommPattern<Scalar> *commPat, DistVec<Scalar> &W)
{

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) 
    {
      subDomain[iSub]->sndData(*commPat, reinterpret_cast<Scalar (*)[1]>(W.subData(iSub)));
    }

  commPat->exchange();
  
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    {
      subDomain[iSub]->addRcvData(*commPat, reinterpret_cast<Scalar (*)[1]>(W.subData(iSub)));
    }

}
//------------------------------------------------------------------------------

template<class Scalar, class OpType >
void Domain::assemble(CommPattern<Scalar> *commPat, DistVec<Scalar> &W, const OpType& oper)
{

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) 
    {
      subDomain[iSub]->sndData(*commPat, reinterpret_cast<Scalar (*)[1]>(W.subData(iSub)));
    }

  commPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    {
      subDomain[iSub]->operateRcvData(*commPat, reinterpret_cast<Scalar (*)[1]>(W.subData(iSub)),oper);
    }

}

//------------------------------------------------------------------------------
template<int dim>
void Domain::assembleGhostPoints(DistVec<GhostPoint<dim>*> &ghostPoints, VarFcn *varFcn)
{
  int iSub;
  // Adam 2010.10.27
  // Caution, the order of the calls matters, because a ghost point can lie on a domain boundary, 
  // in which case we may want to create its state after during the exchange. The num ghost 
  // states is going to be used as a parameter.
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) 
    {
      subDomain[iSub]->sndNumGhostStates(*levelPat, ghostPoints(iSub));
	  subDomain[iSub]->sndGhostStates(*vecPat, ghostPoints(iSub), 0);
    }

  levelPat->exchange();
  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    {
      subDomain[iSub]->rcvNumGhostStates(*levelPat, ghostPoints(iSub), varFcn);
	  subDomain[iSub]->rcvGhostStates(*vecPat, ghostPoints(iSub), 0);
    }

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) 
    {
	  subDomain[iSub]->sndGhostWeights(*vecPat, ghostPoints(iSub), 0);
      subDomain[iSub]->sndGhostTags(*levelPat, ghostPoints(iSub));
    }

  vecPat->exchange();
  levelPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    {
	  subDomain[iSub]->rcvGhostWeights(*vecPat, ghostPoints(iSub), 0);
      subDomain[iSub]->rcvGhostTags(*levelPat, ghostPoints(iSub));
    }

  ////

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) 
	  subDomain[iSub]->sndGhostStates(*vecPat, ghostPoints(iSub), dim);

  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
	  subDomain[iSub]->rcvGhostStates(*vecPat, ghostPoints(iSub), dim);
  
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) 
	  subDomain[iSub]->sndGhostWeights(*vecPat, ghostPoints(iSub), dim);
  
  vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
	  subDomain[iSub]->rcvGhostWeights(*vecPat, ghostPoints(iSub), dim);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
bool Domain::readTagFromFile(const char *prefix, int step, double *tag, int *numSteps)
{// unlike readVectorFromFile and writeVectorToFile, this function will not call exit() if the file does not exist -- it simply returns false

  bool fileExists = subDomain[0]->template checkIfFileExists<Scalar,dim>(prefix);

  if (fileExists) {
    int neq;
    double t = subDomain[0]->template readTagFromFile<Scalar,dim>(prefix, step, &neq, numSteps);
    if (tag) *tag = t;
    if (neq != dim) com->printf(1, "*** Warning: mismatch in dim for \'%s\' (%d vs %d)\n", prefix, neq, dim);
    if (step >= *numSteps)  return false;
    return true;
  } else {
    com->fprintf(stdout, "File [%s] does not exist\n", prefix);
    *numSteps = 0;
    *tag = 0;
    return false;
  }

}


//------------------------------------------------------------------------------

template<class Scalar, int dim>
bool Domain::readVectorFromFile(const char *prefix, int step, double *tag,
				DistSVec<Scalar,dim> &U, Scalar* scale)
{
  int neq, numSteps;
  double t = subDomain[0]->template readTagFromFile<Scalar,dim>(prefix, step, &neq, &numSteps);
  if (tag) *tag = t;

  if (neq != dim)
    com->printf(1, "*** Warning: mismatch in dim for \'%s\' (%d vs %d)\n", prefix, neq, dim);

  if (step >= numSteps)
    return false;

//  com->barrier(); //For timing (of i/o) purpose.
  double t0 = timer->getTime();


#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->readVectorFromFile(prefix, step, neq, U(iSub), scale);

  timer->addBinaryReadTime(t0);

  com->fprintf(stdout, "Read solution %d from \'%s\'\n", step, prefix);

  return true;

}


//------------------------------------------------------------------------------

template<class Scalar>
bool Domain::readVectorFromFile(const char *prefix, int step, double *tag,
                                DistVec<Scalar> &U, Scalar* scale)
{
  int neq, numSteps;
  double t = subDomain[0]->template readTagFromFile<Scalar, 1>(prefix, step, &neq, &numSteps);
  if (tag) *tag = t;
  com->fprintf(stdout, "from [%s], neq = %d, numSteps = %d, step = %d\n", prefix, neq, numSteps, step);
  if (neq != 1)
    com->fprintf(stdout, "*** Warning: mismatch in dim for \'%s\' (%d vs 1)\n", prefix, neq);

  if (step >= numSteps)
    return false;

//  com->barrier(); //For timing (of i/o) purpose.
  double t0 = timer->getTime();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->readVectorFromFile(prefix, step, U(iSub));

  timer->addBinaryReadTime(t0);

  com->fprintf(stdout, "Read solution %d from \'%s\'\n", step, prefix);

  return true;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void Domain::writeVectorToFile(const char *prefix, int step, double tag,
			       DistSVec<Scalar,dim> &U, Scalar* scale)
{

  int iSub;

  com->barrier(); //For timing (of i/o) purpose.
  double t0 = timer->getTime();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->template openFileForWriting<Scalar,dim>(prefix, step);

  if (step == 0)
    com->barrier();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->writeVectorToFile(prefix, step, U(iSub), scale);

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->template writeTagToFile<Scalar,dim>(prefix, step, tag);

#ifndef SYNCHRO_WRITE
  sync();
#endif

  timer->addBinaryWriteTime(t0);

  com->printf(1, "Wrote solution %d to \'%s\'\n", step, prefix);

}


//------------------------------------------------------------------------------

template<class Scalar>
void Domain::writeVectorToFile(const char *prefix, int step, double tag,
                               DistVec<Scalar> &U, Scalar* scale)
{

  int iSub;

  com->barrier(); //For timing (of i/o) purpose.
  double t0 = timer->getTime();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->template openFileForWriting<Scalar, 1>(prefix, step);

  if (step == 0)
    com->barrier();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->writeVectorToFile(prefix, step, U(iSub), scale);

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->template writeTagToFile<Scalar, 1>(prefix, step, tag);

#ifndef SYNCHRO_WRITE
  sync();
#endif

  timer->addBinaryWriteTime(t0);

  com->printf(1, "Wrote solution %d to \'%s\'\n", step, prefix);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void Domain::scaleSolution(DistSVec<Scalar,dim> &data, RefVal* refVal)  {

  int iSub;

  double scale[dim];

  scale[0] = refVal->density;
  scale[1] = refVal->density*refVal->velocity;
  scale[2] = refVal->density*refVal->velocity;
  scale[3] = refVal->density*refVal->velocity;
  scale[4] = refVal->energy;
  
  if (dim == 6) {
    scale[5] = refVal->density*refVal->nutilde;
  }

  if (dim == 7) {
    scale[5] = refVal->density*refVal->kenergy;
    scale[6] = refVal->density*refVal->epsilon;
  }

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)  {

    SVec<Scalar, dim> &subData = data(iSub);

    for (int i = 0; i < subData.size(); ++i)
      for (int j = 0; j < dim; ++j)
        subData[i][j] *= scale[j];
  }

}

//------------------------------------------------------------------------------

template<class S1, class S2>
void Domain::computeStiffAndForce(DefoMeshMotionData::Element type, DistSVec<double,3>& X,
                                  DistSVec<double,3>& F, DistMat<S1,3>& K,
                                  DistMat<S2,3>* P, double volStiff, int** ndType)  {

  double t0 = timer->getTime();

  CommPattern<S2>* matPat = 0;
  if (P) matPat = P->getDiagMatPat();

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    GenMat<S2,3>* p = (P) ? &((*P)(iSub)) : 0;
    int* subNdType = (ndType) ? ndType[iSub] : 0;
    subDomain[iSub]->computeStiffAndForce(type, X(iSub), F(iSub), K(iSub), p, volStiff, subNdType);
    subDomain[iSub]->sndData(*vec3DPat, F.subData(iSub));
    if (P) subDomain[iSub]->sndDiagBlocks(*matPat, (*P)(iSub));
  }

  double t = timer->addMeshAssemblyTime(t0);

  vec3DPat->exchange();
  if (P) matPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->addRcvData(*vec3DPat, F.subData(iSub));
    if (P) subDomain[iSub]->addRcvDiagBlocks(*matPat, (*P)(iSub));
  }

  F *= -1.0;

  com->printf(6, "K matrix computation: %f s\n", t);

}

//------------------------------------------------------------------------------

template<int dim>
int Domain::checkSolution(VarFcn *varFcn, DistSVec<double,dim> &U, DistLevelSetStructure *distLSS)
{

  int ierr = 0;

  if(distLSS)
  {
#pragma omp parallel for reduction(+: ierr)
	  for (int iSub = 0; iSub < numLocSub; ++iSub)
		  ierr += subDomain[iSub]->checkSolution(varFcn, U(iSub), &((*distLSS)(iSub)));
  }
  else
  {
#pragma omp parallel for reduction(+: ierr)
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    ierr += subDomain[iSub]->checkSolution(varFcn, U(iSub));
  }

  com->globalSum(1, &ierr);

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim>
int Domain::checkSolution(VarFcn *varFcn, DistSVec<double,dim> &U, DistVec<int> &fluidId, DistLevelSetStructure *distLSS)
{

  int ierr = 0;

  if(distLSS)
  {
#pragma omp parallel for reduction(+: ierr)
	  for (int iSub = 0; iSub < numLocSub; ++iSub)
		  ierr += subDomain[iSub]->checkSolution(varFcn, U(iSub), fluidId(iSub), &((*distLSS)(iSub)));
  }
  else
  {
#pragma omp parallel for reduction(+: ierr)
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    ierr += subDomain[iSub]->checkSolution(varFcn, U(iSub), fluidId(iSub));

  }

  com->globalSum(1, &ierr);

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::fixSolution(VarFcn *varFcn, DistSVec<double,dim> &U, DistSVec<double,dim> &dU,DistVec<int>* fluidId)
{

  int verboseFlag = com->getMaxVerbose();

//#pragma omp parallel for reduction(+: ierr)
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    if (fluidId) {
      subDomain[iSub]->fixSolution(varFcn, U(iSub), dU(iSub), &((*fluidId)(iSub)), verboseFlag);
    } else {
  
      subDomain[iSub]->fixSolution(varFcn,U(iSub),dU(iSub),NULL,verboseFlag);
    }
  }
}

template<int dim>
void Domain::fixSolution2(VarFcn *varFcn, DistSVec<double,dim> &U, DistSVec<double,dim> &dU,DistVec<int>* fluidId)
{

  int verboseFlag = com->getMaxVerbose();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    if (fluidId) {
      subDomain[iSub]->fixSolution2(varFcn, U(iSub), dU(iSub), &((*fluidId)(iSub)), verboseFlag);
    } else {
      subDomain[iSub]->fixSolution2(varFcn,U(iSub),dU(iSub),NULL,verboseFlag);
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
int Domain::checkSolution(VarFcn *varFcn, DistVec<double> &ctrlVol,
                          DistSVec<double,dim> &U, DistVec<int> &fluidId,
                          DistVec<int> &fluidIdn)
{

  int ierr = 0;

#pragma omp parallel for reduction(+: ierr)
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    ierr += subDomain[iSub]->checkSolution(varFcn, ctrlVol(iSub), U(iSub), fluidId(iSub), fluidIdn(iSub));

  com->globalSum(1, &ierr);

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim, int neq>
int Domain::clipSolution(TsData::Clipping ctype, BcsWallData::Integration wtype,
			 VarFcn* varFcn, double* Uin, DistSVec<double,dim>& U)
{

  const DistInfo& distInfo = U.info();

  int size = neq * distInfo.numGlobSub;
  int sizeint = 2 * size;
  int sizedouble = size;

  int* allint = reinterpret_cast<int *>(alloca(sizeof(int) * (sizeint + 1)));
  double* alldouble = reinterpret_cast<double *>(alloca(sizeof(double) * sizedouble));

  int k;
  for (k=0; k<sizeint; ++k)
    allint[k] = 0;
  for (k=0; k<sizedouble; ++k)
    alldouble[k] = 0.0;

  int (*allcmin)[neq] = reinterpret_cast<int (*)[neq]>(allint);
  int (*allpmin)[neq] = reinterpret_cast<int (*)[neq]>(allint + 1 * size);
  double (*allvmin)[neq] = reinterpret_cast<double (*)[neq]>(alldouble);

  int iSub;
  int ierr = 0;
#pragma omp parallel for reduction(+: ierr)
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    int gSub = distInfo.locSubToGlobSub[iSub];
    ierr += subDomain[iSub]->template
      clipSolution<dim,neq>(ctype, wtype, varFcn, Uin, U.getMasterFlag(iSub), U(iSub),
			    allcmin[gSub], allpmin[gSub], allvmin[gSub]);
  }

  allint[sizeint] = ierr;
  com->globalSum(sizeint + 1, allint);
  ierr = allint[sizeint];
  com->globalSum(sizedouble, alldouble);

  int cmin[neq];
  int pmin[neq];
  double vmin[neq];

  for (k=0; k<neq; ++k) {
    cmin[k] = 0;
    pmin[k] = allpmin[0][k];
    vmin[k] = allvmin[0][k];
  }

  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) {
    for (k=0; k<neq; ++k) {
      cmin[k] += allcmin[iSub][k];
      if (allvmin[iSub][k] < vmin[k]) {
	pmin[k] = allpmin[iSub][k];
	vmin[k] = allvmin[iSub][k];
      }
    }
  }

  for (k=0; k<neq; ++k) {
    if (cmin[k] > 0) {
      if (ctype == TsData::NONE)
	com->printf(1, "*** Warning: %d negative %s value%s (min=%e at %d)\n",
		    cmin[k], varFcn->pname(dim-neq+k), cmin[k]>1? "s":"", vmin[k], pmin[k]);
      else if (ctype == TsData::ABS_VALUE)
	com->printf(1, "*** Warning: %d %s value%s clipped with abs (min=%e at %d)\n",
		    cmin[k], varFcn->pname(dim-neq+k), cmin[k]>1? "s":"", vmin[k], pmin[k]);
      else if (ctype == TsData::FREESTREAM)
	com->printf(1, "*** Warning: %d %s value%s clipped at freestream (min=%e at %d)\n",
		    cmin[k], varFcn->pname(dim-neq+k), cmin[k]>1? "s":"", vmin[k], pmin[k]);
      else if (ctype == TsData::CUTOFF)
	com->printf(1, "*** Warning: %d %s value%s clipped at cutoff (min=%e at %d)\n",
		    cmin[k], varFcn->pname(dim-neq+k), cmin[k]>1? "s":"", vmin[k], pmin[k]);
    }
  }

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::checkFailSafe(VarFcn* varFcn, DistSVec<double,dim>& U,
               DistSVec<bool,2>& tag, DistVec<int> *fluidId)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub){
    if(!fluidId)
      subDomain[iSub]->checkFailSafe(varFcn, U(iSub), tag(iSub));
    else
      subDomain[iSub]->checkFailSafe(varFcn, U(iSub), tag(iSub), &(*fluidId)(iSub));
  }

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::checkGradients(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
			    DistSVec<double,dim> &V, DistNodalGrad<dim> &ngrad)
{

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->checkGradientsSetUp(X(iSub), V(iSub));

  ngrad.compute(0, ctrlVol, V);

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->checkGradientsWrite(X(iSub), ngrad(iSub));

  com->barrier();

  exit(1);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::checkMatVecProd(DistSVec<double,dim> &prod,
			     const char *msg)
{


#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->checkMatVecProd(prod(iSub), msg);
}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeForceDerivs(VarFcn *varFcn, DistSVec<double,3> &X, DistSVec<double,dim> &V,
                                DistSVec<double,dim> &deltaU, Vec<double> &modalF,
                                VecSet< DistSVec<double,3> > &mX)  {

  modalF = 0.0;

  int nStrModes = modalF.len;

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    Vec<double> locModalF(nStrModes);

    SVec<double, 3> **locModes = new SVec<double, 3> *[nStrModes];
    for (int iMode = 0; iMode < nStrModes; iMode++)
      locModes[iMode] = &mX[iMode](iSub);

    subDomain[iSub]->computeForceDerivs(varFcn, X(iSub), V(iSub), deltaU(iSub), locModalF, locModes);

    int iMode;
    for(iMode = 0; iMode < nStrModes; ++iMode) {
      #pragma omp critical
      modalF[iMode] += locModalF[iMode];
    }
  }

  com->globalSum(nStrModes, modalF.v);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::zeroInternalVals(DistSVec<double, dim> &v)  {

//#pragma omp parallel for reduction(+: ierr)
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->zeroInternalVals(v(iSub));
}

//------------------------------------------------------------------------------
template<int dim>
void Domain::printVariable(DistSVec<double,dim>&V, VarFcn *vf)
{
  com->barrier();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->printVariable(V(iSub), vf);
  com->barrier();
}
//------------------------------------------------------------------------------
template<int dim>
void Domain::printInletVariable(DistSVec<double,dim>&V)
{
  com->barrier();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->printInletVariable(V(iSub));
  com->barrier();
}
//------------------------------------------------------------------------------
template<int dim>
void Domain::printAllVariable(DistVec<int> &X, DistSVec<double,dim>&V, int it)
{
  com->barrier();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->printAllVariable(X(iSub), V(iSub), numLocSub, it);
  com->barrier();
}
//------------------------------------------------------------------------------
template<int dim>
void Domain::checkExtrapolationValue(DistSVec<double,dim>&U, VarFcn* vf,
                                     DistBcData<dim>& bcData, DistGeoState& geoState)
{
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->checkExtrapolationValue(U(iSub),  vf, bcData(iSub), geoState(iSub));
  com->barrier();
}

//------------------------------------------------------------------------------
template<class Scalar, int neq>
void Domain::printAllMatrix(DistMat<Scalar, neq> &A, int it)
{

  com->barrier();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->printAllMatrix(A(iSub), it);
  com->barrier();
  exit(1);
}
//------------------------------------------------------------------------------
template<int dim>
void Domain::assemble_dWdt(DistSVec<double, dim> &dWdt, DistSVec<double, dim> &Sigma)

{

  assemble(vecPat, dWdt);
  assemble(vecPat, Sigma);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computedWBar_dt(DistSVec<double, dim> &dWBardt, DistSVec<double, dim> &Sigma,
                             DistMacroCellSet *macroCells, DistSVec<double,1> **volRatio,
                             int scopeDepth)
{

#pragma omp parallel for
    for (int iSub=0; iSub < numLocSub; ++iSub) {
      SVec<double,1> **volR = new SVec<double,1>*[2];
      MacroCellSet** macCells = new MacroCellSet*[scopeDepth];

      for (int i = 0; i<2; ++i)
        volR[i]  =  &(*volRatio[i])(iSub);

      for (int i = 0; i<scopeDepth; ++i)
        macCells[i] =   macroCells->obtainMacroCell(iSub, i);

      subDomain[iSub]->computedWBar_dt(macCells, volR, Sigma(iSub), dWBardt(iSub), scopeDepth);

      delete [] volR;
      delete [] macCells;
    }

   assemble(vecPat, dWBardt);

}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeWeightsForEmbeddedStruct(DistSVec<double,3> &X, DistSVec<double,dim> &V, 
               DistVec<double> &Weights, DistSVec<double,dim> &VWeights, DistVec<int> &init,
															DistVec<int> &next_init, DistLevelSetStructure *distLSS, bool externalSI)
{

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) 
    subDomain[iSub]->computeWeightsForEmbeddedStruct(V(iSub),VWeights(iSub),Weights(iSub),
                                                     (*distLSS)(iSub),X(iSub),
																		 init(iSub), next_init(iSub), externalSI);

  assemble(vecPat, VWeights);
  assemble(volPat, Weights);
  operMax<int> opMax;
  assemble(levelPat, next_init, opMax); // PJSA: added opMax here to fix integer overflow
}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeWeightsLeastSquaresForEmbeddedStruct(DistSVec<double,3> &X, DistSVec<double,dim> &V, 
																			DistVec<double> &Weights, DistSVec<double,dim> &VWeights, 
																			DistVec<int> &init, DistVec<int> &next_init, DistLevelSetStructure *distLSS,
																			DistNodalGrad<dim>& DX, bool limit,	DistVec<int>* fid, bool externalSI)
{

  int iSub;
  DistSVec<double,10> *R = new DistSVec<double,10>(getNodeDistInfo());
  DistSVec<int,1> *count = new DistSVec<int,1>(getNodeDistInfo());

#pragma omp parallel for
	for(iSub = 0; iSub < numLocSub; ++iSub) 
	{ 
    subDomain[iSub]->computeWeightsLeastSquaresEdgePartForEmbeddedStruct((*distLSS)(iSub),X(iSub),
																									(*count)(iSub), (*R)(iSub), init(iSub), externalSI);

	subDomain[iSub]->sndData(*weightPhaseChangePat,(*R).subData(iSub));
	subDomain[iSub]->sndData(*levelPat,(*count).subData(iSub));
  }
  weightPhaseChangePat->exchange();
  levelPat->exchange();
 
#pragma omp parallel for
	for(iSub = 0; iSub < numLocSub; ++iSub) 
	{
	  subDomain[iSub]->addRcvData(*(weightPhaseChangePat),(*R).subData(iSub));
	  subDomain[iSub]->addRcvData(*levelPat,(*count).subData(iSub));
		subDomain[iSub]->computeWeightsLeastSquaresNodePartForEmbeddedStruct((*count)(iSub), (*R)(iSub));
	}

#pragma omp parallel for
	for(int iSub = 0; iSub < numLocSub; ++iSub) 
	{
    Vec<int>* subFid = (fid ? (&(*fid)(iSub)) : 0);
	subDomain[iSub]->computeWeightsLeastSquaresForEmbeddedStruct(X(iSub),(*R)(iSub),V(iSub),
																						 Weights(iSub), VWeights(iSub), (*distLSS)(iSub), 
																						 init(iSub), next_init(iSub), DX(iSub), limit, subFid, externalSI);
  }

  assemble(vecPat, VWeights);
  assemble(volPat, Weights);
  operMax<int> opMax;
  assemble(levelPat, next_init, opMax); // PJSA: added opMax here to fix integer overflow
  if (R) delete R;
  if (count) delete count;
}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeWeightsForFluidFluid(DistSVec<double,3> &X, DistSVec<double,dim> &V, 
					 DistVec<double> &Weights, DistSVec<double,dim> &VWeights, DistVec<int> &init,
					 DistVec<int> &next_init, DistLevelSetStructure *distLSS,
					 DistVec<int> &fluidId)
{
  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) 
    subDomain[iSub]->computeWeightsForFluidFluid(V(iSub),VWeights(iSub),Weights(iSub),
						 (distLSS ? &(*distLSS)(iSub) : NULL),
						 X(iSub),
						 init(iSub),next_init(iSub),
						 fluidId(iSub));

  assemble(vecPat, VWeights);
  assemble(volPat, Weights);
  operMax<int> opMax;
  assemble(levelPat, next_init, opMax); // PJSA: added opMax here to fix integer overflow
}

//------------------------------------------------------------------------------

template<int dim>
void Domain::computeWeightsLeastSquaresForFluidFluid(
			   DistSVec<double,3> &X, DistSVec<double,dim> &V, 
			   DistVec<double> &Weights, DistSVec<double,dim> &VWeights, DistVec<int> &init,
			   DistVec<int> &next_init, DistLevelSetStructure *distLSS,
			   DistVec<int> &fluidId,DistNodalGrad<dim>& nodalgrad,
			   bool limit)
{
  int iSub;
  DistSVec<double,10> *R = new DistSVec<double,10>(getNodeDistInfo());
  DistSVec<int,1> *count = new DistSVec<int,1>(getNodeDistInfo());
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) { 
    subDomain[iSub]->computeWeightsLeastSquaresEdgePartForFF((distLSS ? &(*distLSS)(iSub): NULL),X(iSub),
							     (*count)(iSub),(*R)(iSub),init(iSub),fluidId(iSub));
	subDomain[iSub]->sndData(*weightPhaseChangePat,(*R).subData(iSub));
	subDomain[iSub]->sndData(*levelPat,(*count).subData(iSub));
  }
  weightPhaseChangePat->exchange();
  levelPat->exchange();
 
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
	  subDomain[iSub]->addRcvData(*(weightPhaseChangePat),(*R).subData(iSub));
	  subDomain[iSub]->addRcvData(*levelPat,(*count).subData(iSub));
	  subDomain[iSub]->computeWeightsLeastSquaresNodePartForFF((*count)(iSub),(*R)(iSub));
  }

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
	subDomain[iSub]->computeWeightsLeastSquaresForFluidFluid(X(iSub),(*R)(iSub),V(iSub),
			Weights(iSub),VWeights(iSub),(distLSS ? &(*distLSS)(iSub): NULL),
								 init(iSub),next_init(iSub),
								 fluidId(iSub),nodalgrad(iSub),
								 limit);
  }

  assemble(vecPat, VWeights);
  assemble(volPat, Weights);
  operMax<int> opMax;
  assemble(levelPat, next_init, opMax); // PJSA: added opMax here to fix integer overflow
  if (R) delete R;
  if (count) delete count;
}
//------------------------------------------------------------------------------

template<int dim, int dimLS>
void Domain::computeWeightsForEmbeddedStruct(DistSVec<double,3> &X, DistSVec<double,dim> &V,
                                             DistVec<double> &Weights, DistSVec<double,dim> &VWeights, 
                                             DistSVec<double,dimLS> &Phi, DistSVec<double,dimLS> &PhiWeights,
                                             DistVec<int>& init, DistVec<int> &next_init,
                                             DistLevelSetStructure *distLSS, DistVec<int> *fluidId)
{
  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeWeightsForEmbeddedStruct(V(iSub),VWeights(iSub),Phi(iSub), PhiWeights(iSub),
                                                     Weights(iSub), (*distLSS)(iSub), X(iSub),
                                                     init(iSub), next_init(iSub), (*fluidId)(iSub));

  assemble(vecPat, VWeights);
  assemble(phiVecPat, PhiWeights);
  assemble(volPat, Weights);
  operMax<int> opMax;
  assemble(levelPat, next_init, opMax); // PJSA: added opMax here to fix integer overflow
}

//------------------------------------------------------------------------------

template<int dimLS>
void Domain::extrapolatePhiV(DistLevelSetStructure *distLSS, DistSVec<double,dimLS> &PhiV)
{
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->extrapolatePhiV((*distLSS)(iSub), PhiV(iSub));

  assemble(phiVecPat, PhiV);
  // Note: here PhiV is not a distance function. It will be used only as an indicator of the
  //       fluidId on newly uncovered nodes.
}

//------------------------------------------------------------------------------

template<int dim>
void Domain::populateGhostPoints(DistVec<GhostPoint<dim>*> *ghostPoints, DistSVec<double,3> &X, 
											DistSVec<double,dim> &U, DistNodalGrad<dim, double> *ngrad, 
											VarFcn *varFcn,DistLevelSetStructure *distLSS,bool linRecAtInterface, 
											DistVec<int> &tag, bool externalSI, FemEquationTerm *fet)
{

  int iSub;

  if(!externalSI)
  {
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) 
		  subDomain[iSub]->populateGhostPoints((*ghostPoints)(iSub), X(iSub), U(iSub), (*ngrad)(iSub),
															varFcn, (*distLSS)(iSub), linRecAtInterface, tag(iSub));
  
  assembleGhostPoints(*ghostPoints,varFcn);

  for (iSub = 0; iSub < numLocSub; ++iSub) 
		  subDomain[iSub]->reduceGhostPoints((*ghostPoints)(iSub), X(iSub));
  }
  else
  {
#pragma omp parallel for
	  for (iSub = 0; iSub < numLocSub; ++iSub)
		  subDomain[iSub]->populateGhostPoints((*ghostPoints)(iSub), X(iSub), U(iSub), (*ngrad)(iSub),
															varFcn, (*distLSS)(iSub), tag(iSub), fet);
	  
	  assembleGhostPoints(*ghostPoints, varFcn);

	  for (iSub = 0; iSub < numLocSub; ++iSub) 
		  subDomain[iSub]->reduceGhostPoints((*ghostPoints)(iSub), X(iSub));

	  for (iSub = 0; iSub < numLocSub; ++iSub) 
		  subDomain[iSub]->checkGhostPoints((*ghostPoints)(iSub), X(iSub), U(iSub), (*ngrad)(iSub),
														varFcn, (*distLSS)(iSub), tag(iSub));
  }

}

//------------------------------------------------------------------------------
//
template<int dim, class Scalar, int neq>
void Domain::populateGhostJacobian(DistVec<GhostPoint<dim>*> *ghostPoints, 
											  DistSVec<double,dim> &U, FluxFcn** fluxFcn, VarFcn *varFcn, 
											  DistLevelSetStructure *distLSS, DistVec<int> &tag, DistMat<Scalar,neq>& A) 
{

  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) 
    subDomain[iSub]->populateGhostJacobian((*ghostPoints)(iSub), U(iSub), fluxFcn, varFcn, (*distLSS)(iSub), tag(iSub), A(iSub));
 
}

//------------------------------------------------------------------------------

template<int dim>
void Domain::setSIstencil(DistSVec<double,3> &X, DistLevelSetStructure *distLSS, DistVec<int> &fluidId, DistSVec<double,dim> &U)
{

	int iSub;

#pragma omp parallel for
	for (iSub = 0; iSub < numLocSub; ++iSub) 
		subDomain[iSub]->setSIstencil(X(iSub), (*distLSS)(iSub), fluidId(iSub), U(iSub));
			
}

//------------------------------------------------------------------------------

template<int dim>
void Domain::setFEMstencil(DistSVec<double,3> &X, DistLevelSetStructure *distLSS, DistVec<int> &fluidId, DistSVec<double,dim> &U)
{

	int iSub;

#pragma omp parallel for
	for (iSub = 0; iSub < numLocSub; ++iSub)
		subDomain[iSub]->setFEMstencil(X(iSub), (*distLSS)(iSub), fluidId(iSub), U(iSub));			
}


//------------------------------------------------------------------------------

template<int dim>
void Domain::computeRiemannWeightsForEmbeddedStruct(DistSVec<double,3> &X, DistSVec<double,dim> &V,
               DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji, 
               DistVec<double> &Weights, DistSVec<double,dim> &VWeights, DistLevelSetStructure *distLSS)
{
  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeRiemannWeightsForEmbeddedStruct(V(iSub), Wstarij(iSub), Wstarji(iSub), 
                                                     VWeights(iSub),Weights(iSub), (*distLSS)(iSub),X(iSub));
  assemble(vecPat, VWeights);
  assemble(volPat, Weights);
}

//-------------------------------------------------------------------------------

template<int dim, int dimLS>
void Domain::computeRiemannWeightsForEmbeddedStruct(DistSVec<double,3> &X, DistSVec<double,dim> &V,
               DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji, 
               DistVec<double> &Weights, DistSVec<double,dim> &VWeights,
               DistSVec<double,dimLS> &Phi, DistSVec<double,dimLS> &PhiWeights,
               DistLevelSetStructure *distLSS, DistVec<int> *fluidId0, DistVec<int> *fluidId)
{
  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->computeRiemannWeightsForEmbeddedStruct(V(iSub), Wstarij(iSub), Wstarji(iSub),
                                                     VWeights(iSub),Weights(iSub), Phi(iSub), PhiWeights(iSub),
                                                     (*distLSS)(iSub),X(iSub), (*fluidId0)(iSub), (*fluidId)(iSub));
  assemble(vecPat, VWeights);
  assemble(phiVecPat, PhiWeights);
  assemble(volPat, Weights);
}

//-------------------------------------------------------------------------------
template<int dimLS>
void Domain::checkNodePhaseChange(DistSVec<double,dimLS> &PhiProduct)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->checkNodePhaseChange(PhiProduct(iSub));

}
//-------------------------------------------------------------------------------
template<int dim>
void Domain::storePreviousPrimitive(DistSVec<double,dim> &V, DistVec<int> &fluidId, 
                                    DistSVec<double,3> &X, DistSVec<double,dim> &Vupdate, 
                                    DistVec<double> &weight){
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->storePreviousPrimitive(V(iSub), fluidId(iSub), X(iSub), Vupdate(iSub), weight(iSub));

  assemble(vecPat, Vupdate);
  assemble(volPat, weight);

}
//-------------------------------------------------------------------------------
template<int dim>
void Domain::IncreasePressure(double p, VarFcn *vf,  DistSVec<double,dim> &U){

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->IncreasePressure(p,vf,U(iSub));

}
//-------------------------------------------------------------------------------
template<int dim>
void Domain::IncreasePressure(double p, VarFcn *vf,  DistSVec<double,dim> &U, DistVec<int> &fluidId){

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->IncreasePressure(p,vf,U(iSub),fluidId(iSub));

}
//-------------------------------------------------------------------------------
template<int dim>
void Domain::padeReconstruction(VecSet<DistSVec<double, dim> >&snapsCoarse, VecSet<DistSVec<double, dim> >&snaps, int *stepParam, double *freqCoarse, double deltaFreq, int nStrMode, int L, int M, int nPoints)
{
  int nSteps = stepParam[0];
  int nSnaps = stepParam[2];
  int nStepsCoarse = stepParam[1];
  int nSnapsCoarse = stepParam[3];


  int iSub;

  for (iSub = 0; iSub < numLocSub; ++iSub)  {

     SVec<double, dim> **locVecSetCoarse = new SVec<double, dim> *[nSnapsCoarse];
     SVec<double, dim> **locVecSet = new SVec<double, dim> *[nSnaps];
      for (int iSnapsCoarse = 0; iSnapsCoarse < nSnapsCoarse; iSnapsCoarse++){
        locVecSetCoarse[iSnapsCoarse] = &snapsCoarse[iSnapsCoarse](iSub);

      }
      for (int iSnaps = 0; iSnaps < nSnaps; iSnaps++) {

        locVecSet[iSnaps] = &snaps[iSnaps](iSub);


      }


      subDomain[iSub]->padeReconstruction(locVecSetCoarse, locVecSet, stepParam, freqCoarse, deltaFreq, nStrMode, L, M, nPoints);

  }

}

//------------------------------------------------------------------------------
template<int dim>
void Domain::hardyInterpolationLogMap(VecSet<DistSVec<double, dim> >**logMap, VecSet<DistSVec<double, dim> >&logMapInterp, int nData, int numPod, int iDataMin, FullM &B, FullM &b)
{

  SVec<double, dim> ***locVecSet = new SVec<double, dim> **[nData];
  for (int iSub = 0; iSub < numLocSub; ++iSub)  {
    for (int iData=0; iData < nData; ++iData) {
      if (iData != iDataMin) {
        locVecSet[iData] = new SVec<double, dim> *[numPod];
        for (int iPod = 0; iPod < numPod; ++iPod)
          (locVecSet[iData])[iPod] = &((*(logMap[iData]))[iPod])(iSub);
      }
    }
    SVec<double, dim> **locLogMapInterp = new SVec<double, dim> *[numPod];
    for (int iPod = 0; iPod < numPod; ++iPod)
      locLogMapInterp[iPod] = &logMapInterp[iPod](iSub);
    subDomain[iSub]->hardyInterpolationLogMap(locVecSet,locLogMapInterp,nData,numPod,iDataMin,B,b);
    for (int iPod = 0; iPod < numPod; ++iPod)
      locLogMapInterp[iPod] = 0;
    delete[] locLogMapInterp;
  }
  for (int iData=0; iData < nData; ++iData) {
    if (iData != iDataMin){
      for (int iPod=0; iPod<numPod;++iPod)
        (locVecSet[iData])[iPod] = 0;
      delete [] locVecSet[iData];
    }
  }
  delete [] locVecSet;
}
//------------------------------------------------------------------------------
// Included (MB)
template<int dim>
void Domain::getGradP(DistNodalGrad<dim>& ngrad)
{

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->getGradP(ngrad(iSub));

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Domain::getDerivativeOfGradP(DistNodalGrad<dim>& ngrad)
{

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->getDerivativeOfGradP(ngrad(iSub));

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void Domain::getDerivativeOfGradP(RectangularSparseMat<double,dim,3> **dGradPdddx,
                                  RectangularSparseMat<double,dim,3> **dGradPdddy,
                                  RectangularSparseMat<double,dim,3> **dGradPdddz,
                                  DistSVec<double,dim> &dddx,
                                  DistSVec<double,dim> &dddy,
                                  DistSVec<double,dim> &dddz,
                                  DistSVec<double,3> &dGradP)
{

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->getDerivativeOfGradP(*dGradPdddx[iSub], *dGradPdddy[iSub], *dGradPdddz[iSub], dddx(iSub), dddy(iSub), dddz(iSub), dGradP(iSub));

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void Domain::getTransposeDerivativeOfGradP(RectangularSparseMat<double,dim,3> **dGradPdddx,
                                           RectangularSparseMat<double,dim,3> **dGradPdddy,
                                           RectangularSparseMat<double,dim,3> **dGradPdddz,
                                           DistSVec<double,3> &dGradP,
                                           DistSVec<double,dim> &dddx,
                                           DistSVec<double,dim> &dddy,
                                           DistSVec<double,dim> &dddz)
{

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->getTransposeDerivativeOfGradP(*dGradPdddx[iSub], *dGradPdddy[iSub], *dGradPdddz[iSub], dGradP(iSub), dddx(iSub), dddy(iSub), dddz(iSub));

  CommPattern<double> *vPat = getCommPat(dddx);
  assemble(vPat, dddx);
  assemble(vPat, dddy);
  assemble(vPat, dddz);

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void Domain::computeDerivativeOperatorOfGradP(RectangularSparseMat<double,dim,3> **dGradPdddx,
                                              RectangularSparseMat<double,dim,3> **dGradPdddy,
                                              RectangularSparseMat<double,dim,3> **dGradPdddz)
{

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->computeDerivativeOperatorOfGradP(*dGradPdddx[iSub], *dGradPdddy[iSub], *dGradPdddz[iSub]);

}

//-----------------------------------------------------------------------------

template<int dim>
void Domain::computeCVBasedForceLoad(int forceApp, int orderOfAccuracy, DistGeoState& geoState,
                                     DistSVec<double,3> &X, double (*Fs)[3], int sizeFs,
                                     DistLevelSetStructure *distLSS, double pInfty,
                                     DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji, 
                                     DistSVec<double,dim> &V, DistVec<GhostPoint<dim>*> *ghostPoints, 
                                     PostFcn *postFcn, DistNodalGrad<dim, double> *ngrad, VarFcn* vf, DistVec<int> *fid)
{
  typedef double array3d[3];
  array3d **subFs = new array3d * [numLocSub];
  for(int i=0; i<numLocSub; ++i) subFs[i] = new array3d[sizeFs];

  Vec<GhostPoint<dim>*> *gp=0;

#pragma omp parallel for
	for(int iSub=0; iSub<numLocSub; iSub++) 
	{
    for (int is=0; is<sizeFs; is++) subFs[iSub][is][0] = subFs[iSub][is][1] = subFs[iSub][is][2] = 0.0;
    if(ghostPoints) gp = ghostPoints->operator[](iSub);
    subDomain[iSub]->computeCVBasedForceLoad(forceApp, orderOfAccuracy, geoState(iSub), X(iSub), subFs[iSub],
                                             sizeFs, (*distLSS)(iSub), pInfty, Wstarij(iSub), Wstarji(iSub),
                                             V(iSub),gp,postFcn,(*ngrad)(iSub),vf,fid ? (&(*fid)(iSub)) : 0);
  }
	
	for(int is=0; is<sizeFs; is++) 
	{
    Fs[is][0] = subFs[0][is][0];
    Fs[is][1] = subFs[0][is][1];
    Fs[is][2] = subFs[0][is][2];
  }
  for (int iSub=1; iSub<numLocSub; iSub++)
	{
		for (int is=0; is<sizeFs; is++) 
		{
      Fs[is][0] += subFs[iSub][is][0];
      Fs[is][1] += subFs[iSub][is][1];
      Fs[is][2] += subFs[iSub][is][2];
    }
	}

  for(int i=0; i<numLocSub; ++i) delete [] subFs[i];
  delete [] subFs;
}

//-------------------------------------------------------------------------------
template<int dim>
void Domain::computeEmbSurfBasedForceLoad(IoData &iod, int forceApp, int orderOfAccuracy, 
														DistSVec<double,3> &X, double (*Fs)[3], int sizeFs, 
														DistLevelSetStructure *distLSS, double pInfty, 
                                          DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji, 
                                          DistSVec<double,dim> *Wextij,
                                          DistSVec<double,dim> &V, DistVec<GhostPoint<dim>*> *ghostPoints, 
														PostFcn *postFcn, DistNodalGrad<dim, double> *ngrad, 
														VarFcn* vf, DistVec<int> *fid, bool externalSI)
{

  typedef double array3d[3];
  array3d **subFs = new array3d * [numLocSub];
  for(int i=0; i<numLocSub; ++i) subFs[i] = new array3d[sizeFs];

  int numStructElems = distLSS->getNumStructElems();
  int (*stElem)[3] = distLSS->getStructElems();
  Vec<Vec3D>& Xstruct = distLSS->getStructPosition();

  Vec<GhostPoint<dim>*> *gp=0;

	int**  StNodeDir;

	double** StX1;
	double** StX2;

	double* tmp1;
	double* tmp2;

	if(externalSI)
	{
		StNodeDir = new int*    [numStructElems];
		StX1      = new double* [numStructElems];
		StX2      = new double* [numStructElems];

		for(int i=0; i<numStructElems; ++i)
		{
			StNodeDir[i] = new int    [2];
			StX1[i]      = new double [3];
			StX2[i]      = new double [3];

			StNodeDir[i][0] = StNodeDir[i][1] = 0;
			StX1[i][0] = StX1[i][1] = StX1[i][2] = 0.0;
			StX2[i][0] = StX2[i][1] = StX2[i][2] = 0.0;
		}
	}

	bool rebuildTree = false;
	if(iod.problem.framework == ProblemData::EMBEDDEDALE) rebuildTree = true;

#pragma omp parallel for
	for(int iSub=0; iSub<numLocSub; iSub++) 
	{
		for(int is=0; is<sizeFs; is++) 
		{
			subFs[iSub][is][0] = 0.0;
			subFs[iSub][is][1] = 0.0;
			subFs[iSub][is][2] = 0.0;
		}
		
		if(!externalSI)
		{
    if(ghostPoints) gp = ghostPoints->operator[](iSub);
    subDomain[iSub]->computeEmbSurfBasedForceLoad(iod, forceApp, orderOfAccuracy, 
																		 X(iSub), subFs[iSub], sizeFs, 
																		 numStructElems, stElem, Xstruct,
						  (*distLSS)(iSub), pInfty, 
																		 Wstarij(iSub), Wstarji(iSub), V(iSub), 
																		 gp, postFcn, 
						  (*ngrad)(iSub), vf, fid?&((*fid)(iSub)):0);
  }
		else
			 subDomain[iSub]->computeEMBNodeScalarQuantity_step1(X(iSub), V(iSub),
																				  numStructElems, stElem, Xstruct, 
																				  (*distLSS)(iSub), 
																				  StNodeDir, StX1, StX2, rebuildTree);
	}

	if(externalSI)
	{
		tmp1 = new double[numStructElems];
		tmp2 = new double[numStructElems];

		for(int i=0; i<numStructElems; ++i) 
		{ 
			tmp1[i] = (double) StNodeDir[i][0];
			tmp2[i] = (double) StNodeDir[i][1];
		}

		com->globalSum(numStructElems, tmp1);
		com->globalSum(numStructElems, tmp2);

		for(int i=0; i<numStructElems; ++i) 
		{
			StNodeDir[i][0] = (int) tmp1[i];
			StNodeDir[i][1] = (int) tmp2[i];
		}
	  
		for(int k=0; k<3; ++k)
		{
			for(int i=0; i<numStructElems; ++i) 
			{
				tmp1[i] = StX1[i][k];
				tmp2[i] = StX2[i][k];
			}
		  
			com->globalSum(numStructElems, tmp1);
			com->globalSum(numStructElems, tmp2);
		  
			for(int i=0; i<numStructElems; ++i) 
			{
				StX1[i][k] = tmp1[i];
				StX2[i][k] = tmp2[i];
			}
		}
	}

#pragma omp parallel for
	for (int iSub=0; iSub<numLocSub; iSub++) 
	{
		if(externalSI)
		{
			if(ghostPoints) gp = ghostPoints->operator[](iSub);
			
			 subDomain[iSub]->computeEmbSurfBasedForceLoad_e(iod, forceApp, orderOfAccuracy, 
			 																X(iSub), subFs[iSub], sizeFs, 
			 																numStructElems, stElem, Xstruct,
			 																(*distLSS)(iSub), pInfty, V(iSub), 
			 																gp, postFcn, 
			 																(*ngrad)(iSub), vf, fid?&((*fid)(iSub)):0,
			 																StNodeDir, StX1, StX2);
		}
	}
	
  double res = 0.0;

	for(int is=0; is<sizeFs; is++) 
	{
    Fs[is][0] = subFs[0][is][0];
    Fs[is][1] = subFs[0][is][1];
    Fs[is][2] = subFs[0][is][2];
    res += Fs[is][0]*Fs[is][0] +  Fs[is][1]*Fs[is][1] +  Fs[is][2]*Fs[is][2];
  }
  for (int iSub=1; iSub<numLocSub; iSub++)
	{
		for (int is=0; is<sizeFs; is++) 
		{
      Fs[is][0] += subFs[iSub][is][0];
      Fs[is][1] += subFs[iSub][is][1];
      Fs[is][2] += subFs[iSub][is][2];
    }
	}

  for(int i=0; i<numLocSub; ++i) delete [] subFs[i];
  delete [] subFs;

	if(externalSI)
	{
		for(int i=0; i<numStructElems; ++i)
		{
			delete [] StNodeDir[i];
			delete [] StX1[i];
			delete [] StX2[i];
		}

		delete [] StNodeDir;
		delete [] StX1;
		delete [] StX2;
		delete [] tmp1;
		delete [] tmp2;
	}

}

//-------------------------------------------------------------------------------
template<int dim>
void Domain::computederivativeEmbSurfBasedForceLoad(IoData &iod, int forceApp, int orderOfAccuracy, DistSVec<double,3> &X, 
						    double (*dFs)[3], int sizedFs, DistLevelSetStructure *distLSS, 
						    double pInfty, double dpInfty,
						    DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji, 
						    DistSVec<double,dim> &V, DistSVec<double,dim> &dV_,
						    DistVec<GhostPoint<dim>*> *ghostPoints, PostFcn *postFcn, 
						    DistNodalGrad<dim, double> *gradV, DistNodalGrad<dim, double> *graddV,
						    VarFcn* vf, DistVec<int> *fid){

  typedef double array3d[3];
  array3d **subdFs = new array3d * [numLocSub];
  for(int i=0; i<numLocSub; ++i) subdFs[i] = new array3d[sizedFs];

  int numStructElems   = distLSS->getNumStructElems();
  int (*stElem)[3]     = distLSS->getStructElems();
  Vec<Vec3D>&  Xstruct = distLSS->getStructPosition();
  Vec<Vec3D>& dXstruct = distLSS->getStructDerivative();

  Vec<GhostPoint<dim>*> *gp=0;

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++) {

    for (int is=0; is<sizedFs; is++) subdFs[iSub][is][0] = subdFs[iSub][is][1] = subdFs[iSub][is][2] = 0.0;

    subDomain[iSub]->computederivativeEmbSurfBasedForceLoad(iod, forceApp, orderOfAccuracy, X(iSub), 
							    subdFs[iSub], sizedFs, 
							    numStructElems, stElem, Xstruct, dXstruct, (*distLSS)(iSub), 
							    pInfty, dpInfty,
							    Wstarij(iSub), Wstarji(iSub), V(iSub), dV_(iSub), gp, postFcn, 
							    (*gradV)(iSub), (*graddV)(iSub), vf, fid?&((*fid)(iSub)):0);
  }

  for (int is=0; is<sizedFs; is++) {
    dFs[is][0] = subdFs[0][is][0];
    dFs[is][1] = subdFs[0][is][1];
    dFs[is][2] = subdFs[0][is][2];
  }

#pragma omp parallel for
  for (int iSub=1; iSub<numLocSub; iSub++)
    for (int is=0; is<sizedFs; is++) {
      dFs[is][0] += subdFs[iSub][is][0];
      dFs[is][1] += subdFs[iSub][is][1];
      dFs[is][2] += subdFs[iSub][is][2];
    }

  for(int i=0; i<numLocSub; ++i) delete [] subdFs[i];
  delete [] subdFs;

}
//-------------------------------------------------------------------------------

template<int dim>
void Domain::computeRecSurfBasedForceLoad(int forceApp, int orderOfAccuracy, DistSVec<double,3> &X, 
                                          double (*Fs)[3], int sizeFs, DistLevelSetStructure *distLSS, double pInfty, 
                                          DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji, 
                                          DistSVec<double,dim> &V, 
                                          DistVec<GhostPoint<dim>*> *ghostPoints, PostFcn *postFcn, VarFcn* vf, DistVec<int> *fid)
{
  typedef double array3d[3];
  array3d **subFs = new array3d * [numLocSub];
  for(int i=0; i<numLocSub; ++i) subFs[i] = new array3d[sizeFs];

  Vec<GhostPoint<dim>*> *gp=0;
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++) {
    for (int is=0; is<sizeFs; is++) subFs[iSub][is][0] = subFs[iSub][is][1] = subFs[iSub][is][2] = 0.0;
    if(ghostPoints) gp = ghostPoints->operator[](iSub);
    subDomain[iSub]->computeRecSurfBasedForceLoad(forceApp, orderOfAccuracy, X(iSub), subFs[iSub], sizeFs,
						  (*distLSS)(iSub), pInfty, 
						  Wstarij(iSub), Wstarji(iSub), V(iSub), 
						  gp, postFcn, vf, fid?&((*fid)(iSub)):0);
  }
  for (int is=0; is<sizeFs; is++) {
    Fs[is][0] = subFs[0][is][0];
    Fs[is][1] = subFs[0][is][1];
    Fs[is][2] = subFs[0][is][2];
  }
  for (int iSub=1; iSub<numLocSub; iSub++)
    for (int is=0; is<sizeFs; is++) {
      Fs[is][0] += subFs[iSub][is][0];
      Fs[is][1] += subFs[iSub][is][1];
      Fs[is][2] += subFs[iSub][is][2];
    }

  for(int i=0; i<numLocSub; ++i) delete [] subFs[i];
  delete [] subFs;
}

//-------------------------------------------------------------------------------

template<int dim>
void Domain::computePrdtWCtrlVolRatio(DistSVec<double,dim> &ratioTimesU, DistSVec<double,dim> &U, DistVec<double> &ctrlVol, DistGeoState &geoState) {
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->computePrdtWCtrlVolRatio(ratioTimesU(iSub), U(iSub), ctrlVol(iSub), geoState(iSub));

}

//-------------------------------------------------------------------------------
//---------------------- LEVEL SET (PHI) ----------------------------------------
//-------------------------------------------------------------------------------


template<int dimLS>
void Domain::avoidNewPhaseCreation(DistSVec<double,dimLS> &Phi, DistSVec<double,dimLS> &Phin){

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->avoidNewPhaseCreation(Phi(iSub), Phin(iSub));

}

//------------------------------------------------------------------------------
template<int dimLS>
void Domain::avoidNewPhaseCreation(DistSVec<double,dimLS> &Phi, DistSVec<double,dimLS> &Phin, DistVec<double> &weight, DistLevelSetStructure *distLSS, 
                                   DistVec<int>* fluidIdToSet){

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->avoidNewPhaseCreation(Phi(iSub), Phin(iSub),weight(iSub), distLSS ? &((*distLSS)(iSub)) : 0,
                fluidIdToSet? &((*fluidIdToSet)(iSub)):0);
}
//------------------------------------------------------------------------------

template<int dimLS>
void Domain::setupPhiVolumesInitialConditions(const int volid, const int fluidId, DistSVec<double,dimLS> &Phi){

  // It is assumed that the initialization using volumes is only
  // called to distinguish nodes that are separated by a material
  // interface (structure). Thus one node cannot be at
  // the boundary of two fluids. A fluid node then gets its
  // id from the element id and there cannot be any problem
  // for parallelization.
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->setupPhiVolumesInitialConditions(volid, fluidId, Phi(iSub));

}

//------------------------------------------------------------------------------

template<int dimLS>
void Domain::TagInterfaceNodes(int lsdim, DistVec<int> &Tag, DistSVec<double,dimLS> &Phi,
                               int level, DistLevelSetStructure *distLSS)
{
  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub){
    if (distLSS)
      subDomain[iSub]->TagInterfaceNodes(lsdim, Tag(iSub),Phi(iSub),level,&((*distLSS)(iSub)));
    else
      subDomain[iSub]->TagInterfaceNodes(lsdim, Tag(iSub),Phi(iSub),level,NULL);

    subDomain[iSub]->sndData(*levelPat, reinterpret_cast<int (*)[1]>(Tag.subData(iSub)));
  }

  levelPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->maxRcvData(*levelPat, reinterpret_cast<int (*)[1]>(Tag.subData(iSub)));
}

//------------------------------------------------------------------------------

template<int dimLS>
void Domain::pseudoFastMarchingMethod(DistVec<int> &Tag, DistSVec<double,3> &X, 
				DistSVec<double,dimLS> &d2wall, int level,  int iterativeLevel,
				DistVec<int> &sortedNodes, int *nSortedNodes,
				int *firstCheckedNode, DistLevelSetStructure *distLSS)
{
  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub){
    subDomain[iSub]->pseudoFastMarchingMethod<dimLS>(Tag(iSub),X(iSub),d2wall(iSub),level,iterativeLevel,sortedNodes(iSub),*(nSortedNodes+iSub),*(firstCheckedNode+iSub),distLSS?&((*distLSS)(iSub)):NULL);
    subDomain[iSub]->sndData(*levelPat, reinterpret_cast<int (*)[1]>(Tag.subData(iSub)));
    subDomain[iSub]->sndData(*volPat, reinterpret_cast<double (*)[dimLS]>(d2wall.subData(iSub)));
  }

  levelPat->exchange();
  volPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->maxRcvDataAndCountUpdates(*levelPat, reinterpret_cast<int (*)[1]>(Tag.subData(iSub)),nSortedNodes[iSub],sortedNodes(iSub));
    subDomain[iSub]->minRcvData(*volPat, reinterpret_cast<double (*)[dimLS]>(d2wall.subData(iSub)));
  }
}

//------------------------------------------------------------------------------

template<int dimLS>
void Domain::TagInterfaceNodes(int lsdim, DistSVec<bool,2> &Tag, DistSVec<double,dimLS> &Phi,
                               DistLevelSetStructure *distLSS)
{
  int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->TagInterfaceNodes(lsdim, Tag(iSub), Phi(iSub), &((*distLSS)(iSub)));
     
    subDomain[iSub]->sndData(*bool2Pat, reinterpret_cast<bool (*)[2]>(Tag.subData(iSub)));
  }
  bool2Pat->exchange();
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->maxRcvData(*bool2Pat, reinterpret_cast<bool (*)[2]>(Tag.subData(iSub)));
}

//------------------------------------------------------------------------------
/*template<int dimLS>
void Domain::FinishReinitialization(DistVec<int> &Tag, DistSVec<double,dimLS> &Psi,
                                    int level)
{

	int iSub;
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub){
    subDomain[iSub]->FinishReinitialization(Tag(iSub),Psi(iSub),level);
    subDomain[iSub]->sndData(*levelPat, reinterpret_cast<int (*)[1]>(Tag.subData(iSub)));
    subDomain[iSub]->sndData(*volPat, Psi.subData(iSub));
  }

  levelPat->exchange();
  volPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->TagPsiExchangeData(*levelPat, reinterpret_cast<int (*)[1]>(Tag.subData(iSub)),
                                        *volPat, Psi.subData(iSub));

}*/
//------------------------------------------------------------------------------

template<int dimLS>
void Domain::printPhi(DistSVec<double,3> &X, DistSVec<double,dimLS> &Phi, int it)
{
  com->barrier();
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->printPhi(X(iSub), Phi(iSub), numLocSub);
  com->barrier();
}

//-------------------------------------------------------------------------------

template<int dimLS>
void Domain::computePrdtPhiCtrlVolRatio(DistSVec<double,dimLS> &ratioTimesPhi, DistSVec<double,dimLS> &Phi, DistVec<double> &ctrlVol, DistGeoState &geoState) {

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->computePrdtPhiCtrlVolRatio(ratioTimesPhi(iSub), Phi(iSub), ctrlVol(iSub), geoState(iSub));

}

//-------------------------------------------------------------------------------

template<int dim>
void Domain::restrictionOnPhi(DistSVec<double,dim> &initial, DistVec<int> &fluidId,
                              DistSVec<double,dim> &restriction, int fluidIdTarget){

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->restrictionOnPhi(initial(iSub),fluidId(iSub),
                                      restriction(iSub), fluidIdTarget);

}

//-------------------------------------------------------------------------------
template<int dimLS>
void Domain::getSignedDistance(int lsdim, DistSVec<double,1> &Psi, DistSVec<double,dimLS> &Phi)
{
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->getSignedDistance(lsdim, Psi(iSub),Phi(iSub));
}
//------------------------------------------------------------------------------
template<int dimLS>
void Domain::computeDistanceCloseNodes(int lsdim, DistVec<int> &Tag, DistSVec<double,3> &X,
                                       DistNodalGrad<dimLS> &lsgrad,
                                       DistSVec<double,dimLS> &Phi,DistSVec<double,1> &Psi,
                                       MultiFluidData::CopyCloseNodes copy)
{
  if(copy==MultiFluidData::FALSE){
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub){
    subDomain[iSub]->computeDistanceCloseNodes(lsdim, Tag(iSub), X(iSub), lsgrad(iSub), Phi(iSub), Psi(iSub));
    //subDomain[iSub]->sndData(*phiVecPat, Psi.subData(iSub));
    subDomain[iSub]->sndData(*volPat, Psi.subData(iSub));
  }

  volPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub){
    //subDomain[iSub]->minRcvData(*phiVecPat, Psi.subData(iSub));
    subDomain[iSub]->minRcvData(*volPat, Psi.subData(iSub));
    subDomain[iSub]->recomputeDistanceCloseNodes(lsdim, Tag(iSub), X(iSub), lsgrad(iSub), Phi(iSub), Psi(iSub));
    //subDomain[iSub]->sndData(*phiVecPat, Psi.subData(iSub));
    subDomain[iSub]->sndData(*volPat, Psi.subData(iSub));
  }

  volPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    //subDomain[iSub]->minRcvData(*phiVecPat, Psi.subData(iSub));
    subDomain[iSub]->minRcvData(*volPat, Psi.subData(iSub));
  }else{
#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    subDomain[iSub]->copyCloseNodes(lsdim, 1,Tag(iSub),Phi(iSub),Psi(iSub));
  }

}
//-------------------------------------------------------------------------------
template<int dimLS>
void Domain::computeDistanceLevelNodes(int lsdim, DistVec<int> &Tag, int level,
                                       DistSVec<double,3> &X,DistSVec<double,1> &Psi,
                                       double &_res, DistSVec<double,dimLS> &Phi,
                                       MultiFluidData::CopyCloseNodes copy)
{
  double res(_res);

	if(copy == MultiFluidData::TRUE && level==2)
	{ //KW: why?
#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subDomain[iSub]->copyCloseNodes(lsdim, 2,Tag(iSub),Phi(iSub),Psi(iSub));
    return;
  }

double locRes = 0.0;
#pragma omp parallel for reduction(+: locRes)
	for(int iSub = 0; iSub < numLocSub; ++iSub)
	{
    locRes += subDomain[iSub]->computeDistanceLevelNodes(lsdim, Tag(iSub), level, X(iSub), Psi(iSub),Phi(iSub));
    //subDomain[iSub]->sndData(*phiVecPat, Psi.subData(iSub));
    subDomain[iSub]->sndData(*volPat, Psi.subData(iSub));
  }
  res += locRes;

  volPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub)
    //subDomain[iSub]->minRcvData(*phiVecPat, Psi.subData(iSub));
    subDomain[iSub]->minRcvData(*volPat, Psi.subData(iSub));

  com->globalSum(1, &res);
  _res = sqrt(res);
}
//------------------------------------------------------------------------------

template<int dim>
void Domain::setupUVolumesInitialConditions(const int volid, double UU[dim],
                                            DistSVec<double,dim> &U)
{
  #pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->setupUVolumesInitialConditions_Step1(volid, UU, U(iSub), *vecPat);
  }
  vecPat->exchange();
  #pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->setupUVolumesInitialConditions_Step2(*vecPat, U(iSub));
  }
}

template<int dim>
void Domain::blur(DistSVec<double,dim> &U, DistSVec<double,dim> &U0)
{
  int iSub;

  DistVec<double> loc_weight(getNodeDistInfo());

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) 
    subDomain[iSub]->blur(U(iSub),U0(iSub),loc_weight(iSub));
  
  assemble(vecPat, U0);
  assemble(volPat,loc_weight);
  
#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {

    for (int i = 0; i < U.subSize(iSub); ++i) {
      for (int k = 0; k < dim; ++k)
	U0(iSub)[i][k] = 0.5*U0(iSub)[i][k]/loc_weight(iSub)[i] + 0.5*U(iSub)[i][k];
    }

  }

}

//------------------------------------------------------------------------------

template<int dimLS>
void Domain::updateFluidIdFS2Prep(DistLevelSetStructure &distLSS, DistSVec<double,dimLS> &PhiV, DistVec<int> &fluidId, DistSVec<bool,4> &poll)
{
  int iSub;
#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub) {
    subDomain[iSub]->solicitFluidIdFS(distLSS(iSub), fluidId(iSub), poll(iSub));
    subDomain[iSub]->sndData(*bool4Pat, reinterpret_cast<bool (*)[4]>(poll.subData(iSub)));
    //subDomain[iSub]->sndData(*bool3Pat, reinterpret_cast<bool (*)[3]>(poll.subData(iSub)));
  }

  bool4Pat->exchange();
  //bool3Pat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < numLocSub; ++iSub) {
    subDomain[iSub]->maxRcvData(*bool4Pat, reinterpret_cast<bool (*)[4]>(poll.subData(iSub)));
    //subDomain[iSub]->maxRcvData(*bool3Pat, reinterpret_cast<bool (*)[3]>(poll.subData(iSub)));
//    subDomain[iSub]->updateFluidIdFS2(distLSS(iSub), PhiV(iSub), poll(iSub), fluidId(iSub), (PhiV.info()).getMasterFlag(iSub));
  }
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void Domain::debugMultiPhysics(DistLevelSetStructure &distLSS, DistSVec<double,dimLS> &PhiV, 
                               DistVec<int> &fluidId, DistSVec<double,dim> &U)
{
  int iSub;
#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub)
    subDomain[iSub]->debugMultiPhysics(distLSS(iSub), PhiV(iSub), fluidId(iSub), U(iSub)); 
}

//------------------------------------------------------------------------------

template<int dim, class Obj>
void Domain::integrateFunction(Obj* obj,DistSVec<double,3> &X,DistSVec<double,dim>& V, void (Obj::*F)(int node, const double* loc,double* f),
			       int npt) {

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub){
    subDomain[iSub]->integrateFunction(obj, X(iSub), V(iSub), F, npt);
    subDomain[iSub]->sndData(*volPat, V(iSub) );
  }

  volPat->exchange();

#pragma omp parallel for
  for (int iSub = 0; iSub < numLocSub; ++iSub) {
    //subDomain[iSub]->minRcvData(*phiVecPat, Psi.subData(iSub));
    subDomain[iSub]->addRcvData(*volPat, V(iSub));
  }
}

//------------------------------------------------------------------------------
/*
template<int dim>
void Domain::readMultiPodBasis(const char *multiPodFile,VecSet< DistSVec<double,dim> > *(pod[2]), int nPod [2], int nBasesNeeded, int *whichFiles) {	

	//	multiPodFile: file containing names of bases
	//	pod: array of pointers to POD bases. Each one is individually uninitialized
	//	nPod: vector of number of required POD basis vectors
	//	nBasesNeeded: number of bases needed for the problem (default is zero, meaning ignore)
	//	whichFiles: indicates which files should be read. if [1 -1], it means
		//	you should read only second file in the list of files. if [0 1], it
		//	means you should read both


	if (whichFiles == NULL) 	// by default, just take the files in the order prescribed
		for (int i = 0; i < nBasesNeeded; ++i)
			whichFiles[i] = i;

	const char *vecFile = multiPodFile;	// already read into the function
	if (!vecFile)
		vecFile = "multiPodFile.in";	// default filename
	FILE *inFP = fopen(vecFile, "r");
	if (!inFP)  {
		com->fprintf(stderr, "*** Warning: No POD FILES in %s\n", vecFile);
		exit (-1);
	}

	// =====================================================
	// read in file containing filenames for multiple bases
	// =====================================================

	int nData; // number of available bases in the file
	fscanf(inFP, "%d",&nData);	// first entry is the number of bases in the file

	if (nBasesNeeded > 0 && nData < nBasesNeeded) {
		com->fprintf(stderr," ... ERROR: expecting %d POD files, found %d\n",nBasesNeeded,nData);
		com->fprintf(stderr," ... Exiting\n");
		exit(-1); 
	}
	com->fprintf(stderr," ... reading in %d POD files\n",nBasesNeeded);

	char **podFile = new char *[nData];	// files as they appear in the list input file

	for (int iData = 0; iData < nData; ++iData){
		podFile[iData] = new char[1000];
		fscanf(inFP, "%s", podFile[iData]);
	}

	// ================================
	// read individual POD bases
	// ================================

	for (int iData=0; iData < nBasesNeeded; ++iData){	// loop over bases
		com->fprintf(stderr, " ... Reading POD from %s \n", podFile[whichFiles[iData]]);
		readPodBasis(podFile[whichFiles[iData]], nPod[iData],*(pod[iData]));
	}
}

//------------------------------------------------------------------------------

template<int dim>
void Domain::readPodBasis(const char *podFile, int &nPod,
		VecSet<DistSVec<double, dim> > &podVecs, bool useSnaps) {

  // read in POD Vectors
  const char *vecFile = podFile;

  double eigValue;
  int nPodVecs;

  // read number of vecs
  DistSVec<double,dim> tmpVec(getNodeDistInfo());

	if (useSnaps)	// reading snapshots, not a pod basis
		nPodVecs = nPod;
	else {
		readVectorFromFile(vecFile, 0, &eigValue, tmpVec);
		nPodVecs = (int) eigValue;
		com->fprintf(stderr, " ... There are %d total podVecs \n", nPodVecs);	// unique POD vectors (first one repeated)
	}

  if (nPod > nPodVecs)  {
    com->fprintf(stderr, " ... WARNING: there are only %d POD Vectors \n", nPodVecs);
    nPod = nPodVecs;
  }
	else if (nPod < 0)  {	// if negative value specified, read in all vectors
		nPod = nPodVecs;
	}
  else
    nPodVecs = nPod;

  com->fprintf(stderr, " ... Reading %d POD Vectors from file %s\n", nPodVecs, vecFile);

  podVecs.resize(nPodVecs);

  int iVec;
  double firstEig;
  for (iVec = 0; iVec < nPodVecs; iVec++) { 
    readVectorFromFile(vecFile, iVec+1, &eigValue, podVecs[iVec]);	// read in one more
		if (iVec == 0)
			firstEig = eigValue;
	}

	double firstEigDisplayed;
	if (firstEig == 0){
		firstEigDisplayed = 1.0;
	}
	else{
		firstEigDisplayed = firstEig;
	}
  com->fprintf(stderr, " ... Eigenvalue Ratio: (%e/%e) = %e\n", eigValue, firstEigDisplayed, eigValue/firstEigDisplayed);
}

*/

template<typename Scalar>
void Domain::communicateMesh(std::vector <Scalar> * nodeOrEle, int arraySize,
		int *alreadyCommunicatedArray){	
	
	// loop over iIslands
		// figure out how many total entries each cpu has for the iIsland
		// 	numNeigh = {0 9 15 58} means that the first cpu has 9, second has 6, etc.
		// initiate memory for the total number of globalNodes (using last entry in above vector)
		// fill out entries [iCpu] to [iCpu+1] in above array using vector
		// do a global sum
		// overwrite node and element vectors with this global data (for each cpu)
	
	int* numNeigh = new int [com->size() + 1];
	int offset;	// initial values to not communicate (if already communicated some)
	for (int iArraySize = 0; iArraySize < arraySize; ++iArraySize) {
		if (com->cpuNum() > 0 && alreadyCommunicatedArray != NULL)
			offset = alreadyCommunicatedArray[iArraySize];
		else
			offset = 0;
		for (int i = 0; i <=com->size(); ++i)	// initialize
			numNeigh[i] = 0;
		numNeigh[com->cpuNum()+1] = nodeOrEle[iArraySize].size() - offset;	// number of entries on this cpu
		com->globalSum(com->size()+1,numNeigh);
		for (int i = 1; i <=com->size(); ++i)	// accumulate
			numNeigh[i] += numNeigh[i-1];
		int totalNodeOrEle = numNeigh[com->size()];	// total across all processors

		Scalar *nodeOrEleArray = new Scalar [totalNodeOrEle];
		for (int iNeighbor = 0; iNeighbor < totalNodeOrEle; ++iNeighbor) {
			if (iNeighbor >= numNeigh[com->cpuNum()] && iNeighbor < numNeigh[com->cpuNum()+1]) 
				nodeOrEleArray[iNeighbor] = nodeOrEle[iArraySize][iNeighbor - numNeigh[com->cpuNum()] + offset];	// fill in this cpu's contribution
			else
				nodeOrEleArray[iNeighbor] = 0;
		}
		com->globalSum(totalNodeOrEle,nodeOrEleArray);

		// fill in the array with all global entries
		nodeOrEle[iArraySize].clear();
                nodeOrEle[iArraySize].resize(totalNodeOrEle,0.0);
		for (int iNeighbor = 0; iNeighbor < totalNodeOrEle; ++iNeighbor) 
			nodeOrEle[iArraySize][iNeighbor]=nodeOrEleArray[iNeighbor];

		delete [] nodeOrEleArray;
		
	}
	delete [] numNeigh;
}

template<typename Scalar>
void Domain::makeUnique( std::vector <Scalar> * nodeOrEle, int length) {

	// remove redundant entries from a vector <int> nodeOrEle *
	// apply to nodeOrEle and elements

	vector<int>::iterator it;

	// put all on one vector
	for (int iIsland = 1; iIsland < length; ++iIsland) {
		for (int iEntry = 0; iEntry < nodeOrEle[iIsland].size(); ++iEntry) {
			nodeOrEle[0].push_back(nodeOrEle[iIsland][iEntry]);
		}
		nodeOrEle[iIsland].erase(nodeOrEle[iIsland].begin(),nodeOrEle[iIsland].end());	// no longer need that vector
	}

	sort(nodeOrEle[0].begin(), nodeOrEle[0].end());	// sort: puts in order
	it = unique(nodeOrEle[0].begin(), nodeOrEle[0].end()); // remove duplicate consecutive elements (reason for sort)
	nodeOrEle[0].resize(it - nodeOrEle[0].begin());	// remove extra entries
}
// Functions to compute the error (that is, the difference between two state vectors)
template <int dim>
void Domain::computeL1Error(DistSVec<double,dim>& U, DistSVec<double,dim>& Uexact, 
			    DistVec<double>& vol, double error[dim],
                            DistLevelSetStructure* distLSS) {

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++) {
    LevelSetStructure*  LSS = (distLSS ? &(*distLSS)(iSub) : NULL);
    subDomain[iSub]->computeL1Error(U.getMasterFlag(iSub),U(iSub), Uexact(iSub), 
				    vol(iSub),error, LSS);
  }

  com->globalSum(dim, error);
}

// Functions to compute the error (that is, the difference between two state vectors)
template <int dim>
void Domain::computeL2Error(DistSVec<double,dim>& U, DistSVec<double,dim>& Uexact, 
			    DistVec<double>& vol, double error[dim],
                            DistLevelSetStructure* distLSS) {

#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++) {
    LevelSetStructure*  LSS = (distLSS ? &(*distLSS)(iSub) : NULL);
    subDomain[iSub]->computeL2Error(U.getMasterFlag(iSub),U(iSub), Uexact(iSub), 
				    vol(iSub),error, LSS);
  }

  com->globalSum(dim, error);

  for (int k = 0; k < dim; ++k) {

    error[k] = sqrt(error[k]);
  }
}



template <int dim>
void Domain::computeLInfError(DistSVec<double,dim>& U, DistSVec<double,dim>& Uexact, double error[dim],DistLevelSetStructure* distLSS) {
  
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++) {
    LevelSetStructure*  LSS = (distLSS ? &(*distLSS)(iSub) : NULL);
    subDomain[iSub]->computeLInfError(U.getMasterFlag(iSub),U(iSub), Uexact(iSub), error,
                                      LSS);
  }

  com->globalMax(dim, error);
}

template <int dim>
void Domain::computeHHBoundaryTermResidual(DistBcData<dim> &bcData,DistSVec<double,dim> &U,DistVec<double>& res,
					   VarFcn* vf) {
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; iSub++)
    subDomain[iSub]->computeHHBoundaryTermResidual(bcData(iSub),U(iSub), res(iSub), vf);
}

template <int dim>
void Domain::setExactBoundaryValues(DistSVec<double,dim>& U, DistSVec<double,3>& X,
				    IoData& iod,double t, VarFcn* varFcn) {

  if (iod.embed.testCase == 1) {

#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; iSub++) {

      int lsize = U(iSub).size();
      for (int i = 0; i < lsize; ++i) {

	double* x = X(iSub)[i];
	if (x[0] == 0.0 || fabs(x[0]-1.0) < 1.0e-12/* || x[1] > 0.98*/) {
	  
	  double V[5];
	  ExactSolution::AcousticBeam(iod,x[0],x[1],x[2],t, V);

	  varFcn->primitiveToConservative(V, U(iSub)[i], 0);
	  
	}
      }
    }
  } else if (iod.embed.testCase == 2) {

#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; iSub++) {

      int lsize = U(iSub).size();
      for (int i = 0; i < lsize; ++i) {

	double* x = X(iSub)[i];
	if (x[0] == 0.0 || fabs(x[0]-1.0) < 1.0e-12/* || x[1] > 0.98*/) {
	  
	  double V[5];
	  ExactSolution::AcousticViscousBeam(iod,x[0],x[1],x[2],t, V);

	  varFcn->primitiveToConservative(V, U(iSub)[i], 0);
	  
	}
      }
    }
  } else if (iod.mf.testCase == 3) {

#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; iSub++) {

      int lsize = U(iSub).size();
      for (int i = 0; i < lsize; ++i) {

	double* x = X(iSub)[i];
	if (x[0] == 0.0 || fabs(x[0]-1.0) < 1.0e-12 ||
	    x[1] == 0.0 || fabs(x[1]-1.0) < 1.0e-12/* || x[1] > 0.98*/) {
	  
	  double V[5];
	  double dummy;
	  int fid;
	  ExactSolution::AcousticTwoFluid(iod,x[0],x[1],x[2],t, V,&dummy, fid);

	  varFcn->primitiveToConservative(V, U(iSub)[i], fid);
	  
	}
      }
    }
  }
}

template <int dim>
void Domain::setExactBoundaryResidual(DistSVec<double,dim>& U, DistSVec<double,3>& X,
				      IoData& iod,double t, VarFcn* varFcn) {

  if (iod.embed.testCase == 1 || iod.embed.testCase == 2 || iod.mf.testCase == 3) {

#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; iSub++) {

      int lsize = U(iSub).size();
      for (int i = 0; i < lsize; ++i) {

	double* x = X(iSub)[i];
	if (x[0] == 0.0 || fabs(x[0]-1.0) < 1.0e-12/* || x[1] > 0.98*/) {
	  
	  memset(U(iSub)[i], 0, sizeof(double)*dim);
	  
	}
      }
    }
  }
}

template <int dim,int neq,class Scalar>
void Domain::setExactBoundaryJacobian(DistSVec<double,dim>& U, DistSVec<double,3>& X,
				      IoData& iod,double t, VarFcn* varFcn,
				      DistMat<Scalar,neq>& A) {

  if (iod.embed.testCase == 1 || iod.embed.testCase == 2 || iod.mf.testCase == 3) {

#pragma omp parallel for
    for (int iSub=0; iSub<numLocSub; iSub++) {

      int (*edgePtr)[2] = subDomain[iSub]->getEdges().getPtr();
      
      GenMat<Scalar,neq> &Asub = A(iSub);
      for (int l=0; l<subDomain[iSub]->getEdges().size(); ++l) {
	int i = edgePtr[l][0];
	int j = edgePtr[l][1];
	
	double* xi = X(iSub)[i];
	double* xj = X(iSub)[j];
	bool isi = (xi[0] == 0.0 || fabs(xi[0]-1.0) < 1.0e-12);
	bool isj = (xj[0] == 0.0 || fabs(xj[0]-1.0) < 1.0e-12);
	
	if (isi)  {
	  Scalar *Aij = Asub.getElem_ij(l);
	  memset(Aij, 0, sizeof(double)*neq*neq);
	}
	
	if (isj)  {
	  Scalar *Aji = Asub.getElem_ji(l);
	  memset(Aji, 0, sizeof(double)*neq*neq);
	}
      }
    }
  }
}

//-------------------------------------------------------------------------------
/**   Function Domain::computeMaterailConservationScalars
   *  Mass is the mass vector to save mass of different fluid materials
   *  size is the numFluidPhases + 1(ghost node)
   *  U is the conservative state variables
   *  A is the volume of fluid control volume
   *  fluidId is the fluid Id vector for mutiphase problem, and NULL for single phase problem
   */
template<int dim>
void Domain::computeMaterialConservationScalars(double *Mass,  double *MomentumX, double *MomentumY, double *MomentumZ, double *Energy,
                                       int size, DistSVec<double,dim> &U, DistVec<double> &A, DistVec<int> *fluidId)
{
  double subMass[numLocSub][size];
  double subMomentumX[numLocSub][size];
    double subMomentumY[numLocSub][size];
    double subMomentumZ[numLocSub][size];
    double subEnergy[numLocSub][size];
#pragma omp parallel for
  for (int iSub=0; iSub<numLocSub; ++iSub) {
    for(int i=0; i<size; i++) {
      subMass[iSub][i] = 0.0;
      subEnergy[iSub][i] = 0.0;
        subMomentumX[iSub][i] = 0.0;
        subMomentumY[iSub][i] = 0.0;
        subMomentumZ[iSub][i] = 0.0;
    }
    double *subA = A.subData(iSub);
    double (*subU)[dim] = U.subData(iSub);
    int *subId = fluidId ? fluidId->subData(iSub) : 0;
    bool *subMasterFlag = A.getMasterFlag(iSub);

    for(int i=0; i<A.subSize(iSub); i++) {
      if(!subMasterFlag[i])
        continue;
      int myId = subId ? subId[i] : 0;
      if(myId>=size) {
        fprintf(stderr,"ERROR: Detected FluidId = %d. Maximum should be %d (or less). \n", myId, size-1);
        exit(-1);
      }
      //subU[i][0] is the density
        subMass[iSub][myId] += subA[i]*subU[i][0];
        subMomentumX[iSub][myId] += subA[i]*subU[i][1];
        subMomentumY[iSub][myId] += subA[i]*subU[i][2];
        subMomentumZ[iSub][myId] += subA[i]*subU[i][3];
        subEnergy[iSub][myId] += subA[i]*subU[i][4];
    }
  }

  for(int iSub=0; iSub<numLocSub; ++iSub)
    for(int i=0; i<size; i++) {
        Mass[i] += subMass[iSub][i];
        MomentumX[i] += subMomentumX[iSub][i];
        MomentumY[i] += subMomentumY[iSub][i];

        MomentumZ[i] += subMomentumZ[iSub][i];
        Energy[i] += subEnergy[iSub][i];
    }

    com->globalSum(size, Mass);
    com->globalSum(size, Energy);
    com->globalSum(size, MomentumX);
    com->globalSum(size, MomentumY);
    com->globalSum(size, MomentumZ);
}
