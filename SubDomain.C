#include <cstdio>
#include <cmath>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
using std::max;
#endif
#ifdef USE_EIGEN3
#include <Eigen/Core>
#include <Eigen/Sparse>
typedef Eigen::SparseMatrix<double> SpMat;
#endif

#include <FluxFcn.h>
#include <RecFcn.h>
#include <MacroCell.h>
#include <VMSLESTerm.h>
#include <DynamicVMSTerm.h>
#include <SmagorinskyLESTerm.h>
#include <WaleLESTerm.h>
#include <DynamicLESTerm.h>
#include <FemEquationTerm.h>
#include <VolumicForceTerm.h>
#include <PostFcn.h>
#include <BcFcn.h>
#include <BcDef.h>
#include <NodalGrad.h>
#include <EdgeGrad.h>
#include <Extrapolation.h>
#include <ExactRiemannSolver.h>
#include <BcData.h>
#include <GeoState.h>
#include <Vector.h>
#include <MvpMatrix.h>
#include <SparseMatrix.h>
#include <DenseMatrixOps.h>
#include <Connectivity.h>
#include <MemoryPool.h>
#include <Communicator.h>
#include <BinFileHandler.h>
#include <VectorSet.h>
#include <LinkF77.h>
#include <LowMachPrec.h>
#include "FluidSelector.h"
#include <GhostPoint.h>
#include <DenseMatrixOps.h>
#include <limits>
#include <PolygonReconstructionData.h> 
#include <Quadrature.h>
#include <sys/stat.h>
#include <RTree.h>
#include <ProgrammedBurn.h>
#include "LevelSet/LevelSetStructure.h"

#include "HigherOrderFSI.h"

#include "FSI/CrackingSurface.h"

extern "C" {
  void F77NAME(mvp5d)(const int &, const int &, int *, int *, int (*)[2],
		      double (*)[25], double (*)[5], double (*)[5]);
  void F77NAME(torsionspring)(double (*)[3], int [4], double (*)[12], double &invCoef);
  void F77NAME(ballvertex)(double (*)[3], double(*)[3], int [4], double (*)[12], double &invCoef);
};

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn, GeoState &geoState,
                                SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &dt,
                                Vec<double> &idti, Vec<double> &idtv, Vec<double> &dtau,
                                TimeLowMachPrec &tprec, LevelSetStructure *LSS)
{

  dt = 0.0;
  idti = 0.0;
  idtv = 0.0;
  dtau = 0.0;

	if(LSS)
	{
		edges.computeTimeStep(fet, varFcn, geoState, X, V, idti, idtv, tprec, LSS);
		if(fet) elems.computeTimeStep(fet, X, V, idtv, LSS);
		faces.computeTimeStep(fet, varFcn, geoState, X, V, idti, idtv, tprec, LSS);
	}
	else
	{
  edges.computeTimeStep(fet, varFcn, geoState, X, V, idti, idtv, tprec);
		if(fet) elems.computeTimeStep(fet, X, V, idtv, LSS);
  faces.computeTimeStep(fet, varFcn, geoState, X, V, idti, idtv, tprec);
	}

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeDerivativeOfTimeStep(FemEquationTerm *fet, VarFcn *varFcn, GeoState &geoState,
                                SVec<double,3> &X, SVec<double,3> &dX, SVec<double,dim> &V, SVec<double,dim> &dV,
                                Vec<double> &dIdti, Vec<double> &dIdtv, double dMach,
                                TimeLowMachPrec &tprec)
{

  dIdti = 0.0;
  dIdtv = 0.0;

  edges.computeDerivativeOfTimeStep(fet, varFcn, geoState, X,  dX, V, dV, dIdti, dIdtv, dMach, tprec);
  faces.computeDerivativeOfTimeStep(fet, varFcn, geoState, X,  dX, V, dV, dIdti, dIdtv, dMach, tprec);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn, GeoState &geoState,
                                SVec<double,dim> &V, Vec<double> &dt,
				Vec<double> &idti, Vec<double> &idtv, Vec<double> &dtau,
                                TimeLowMachPrec &tprec,
				Vec<int> &fluidId, Vec<double>* umax)
{

  dt = 0.0;
  dtau = 0.0;

  edges.computeTimeStep(varFcn, geoState, V, dt, tprec, fluidId, globSubNum,umax);
  faces.computeTimeStep(varFcn, geoState, V, dt, tprec, fluidId);

}

//------------------------------------------------------------------------------

inline
void computeLocalWeightsLeastSquares(double dx[3], double *R, double *W)
{

  if(R[0]*R[3]*R[5] == 0.0) fprintf(stderr, "Going to divide by 0 %e %e %e\n", 
         R[0], R[3], R[5]);
  double or11 = 1.0 / R[0];
  double or22 = 1.0 / R[3];
  double or33 = 1.0 / R[5];

  double r12or11 = R[1] * or11;
  double r23or22 = R[4] * or22;

  double psi = (R[1]*R[4] - R[2]*R[3]) * or11* or22;

  double alpha1 = dx[0] * or11 * or11;
  double alpha2 = (dx[1] - r12or11*dx[0]) * or22 * or22;
  double alpha3 = (dx[2] - r23or22*dx[1] + psi*dx[0]) * or33 * or33;

  W[0] = alpha1 - r12or11*alpha2 + psi*alpha3;
  W[1] = alpha2 - r23or22*alpha3;
  W[2] = alpha3;

}

//------------------------------------------------------------------------------
inline
void computeLocalWeightsLeastSquaresForEmbeddedStruct(double dx[3], double *R, double *W)
{
  if (R[0]*R[4]*R[7]*R[9]==0.0) 
	fprintf(stderr, "Going to be divided by 0 %e %e %e %e\n",R[0],R[4],R[7],R[9]);
  double or11 = 1.0/R[0];
  double or22 = 1.0/R[4];
  double or33 = 1.0/R[7];
  double or44 = 1.0/R[9];

  double r12or11 = R[1]*or11;
  double r23or22 = R[5]*or22;
  double r34or33 = R[8]*or33;

  double psi13 = (R[1]*R[5]-R[2]*R[4])*or11*or22;
  double psi24 = (R[5]*R[8]-R[6]*R[7])*or22*or33;
  double psi14 = (-R[2]*R[8]+R[3]*R[7]+R[1]*R[7]*psi24)*or11*or33;

  double alpha1 = dx[0]*or11*or11;
  double alpha2 = (dx[1]-r12or11*dx[0])*or22*or22;
  double alpha3 = (dx[2]-r23or22*dx[1]+psi13*dx[0])*or33*or33;
  double alpha4 = (1.0-r34or33*dx[2]+psi24*dx[1]-psi14*dx[0])*or44*or44;

  W[0] = alpha1 - r12or11*alpha2 + psi13*alpha3 - psi14*alpha4;
  W[1] = alpha2 - r23or22*alpha3 + psi24*alpha4;
  W[2] = alpha3 - r34or33*alpha4;
  W[3] = alpha4;
}

//------------------------------------------------------------------------------
// Included (YC)
inline 
void compute_dWdXAnddWdR(int pm, double dx[3], double *R, double *W, double dWdX[][6], double dWdR[][6])
{
  if(R[0]*R[3]*R[5] == 0.0) fprintf(stderr, "Going to divide by 0 %e %e %e\n", R[0], R[3], R[5]);

  double or11 = 1.0 / R[0];
  double or22 = 1.0 / R[3];
  double or33 = 1.0 / R[5];

  double r12or11 = R[1] * or11;
  double r23or22 = R[4] * or22;

  double psi = (R[1]*R[4] - R[2]*R[3]) * or11 * or22;
  double alpha1 = dx[0] * or11 * or11;
  double alpha2 = (dx[1] - r12or11*dx[0]) * or22 * or22;
  double alpha3 = (dx[2] - r23or22*dx[1] + psi*dx[0]) * or33 * or33;

  W[0] = alpha1 - r12or11*alpha2 + psi*alpha3;
  W[1] = alpha2 - r23or22*alpha3;
  W[2] = alpha3;

  double co11 = (or11*or11 + psi*or33*or33*psi + r12or11*or22*or22*r12or11);
  double co12 = (-psi*or33*or33*r23or22 - r12or11*or22*or22);
  double co13 = psi*or33*or33;
  double ro11 = (alpha2*R[1] / ( R[0]*R[0] ) - 2.0*dx[0] * or11 / ( R[0]*R[0] ) - psi*or33*or33*dx[0]*or22*(R[1]*R[4] - R[2]*R[3]) / ( R[0]*R[0] ) - r12or11*or22*or22*dx[0]*R[1] / ( R[0]*R[0] ) - alpha3*or22*(R[1]*R[4] - R[2]*R[3]) / ( R[0]*R[0] ));
  double ro12 = (r12or11*or22*or22*dx[0]*or11 + psi*or33*or33*dx[0]*or11*or22*R[4] + alpha3*or11*or22*R[4] - alpha2*or11);
  double ro13 = (-alpha3*or11*or22*R[3] - psi*or33*or33*dx[0]*or11*or22*R[3]);
  double ro14 = (r12or11*2.0*(dx[1] - r12or11*dx[0]) * or22 / ( R[3]*R[3] ) - psi*or33*or33*dx[0]*(or11*or22*R[2] + (R[1]*R[4] - R[2]*R[3]) * or11 / ( R[3]*R[3] )) + psi*or33*or33*dx[1]*R[4] /( R[3]*R[3] ) - alpha3*(or11*or22*R[2] + (R[1]*R[4] - R[2]*R[3]) * or11 / ( R[3]*R[3] )));
  double ro15 = (alpha3*or11*or22*R[1] + psi*or33*or33*dx[0]*or11*or22*R[1] - psi*or33*or33*dx[1]*or22);
  double ro16 = -psi*2.0*(dx[2] - r23or22*dx[1] + psi*dx[0]) * or33 / ( R[5]*R[5] );
  double co21 =  (-or22*or22*r12or11 - r23or22*or33*or33*psi);
  double co22 = (or22*or22 + r23or22*or33*or33*r23or22);
  double co23 = - r23or22*or33*or33;
  double ro21 = (or22*or22*dx[0]*R[1] / ( R[0]*R[0] ) + r23or22*or33*or33*dx[0]*or22*(R[1]*R[4] - R[2]*R[3]) / ( R[0]*R[0] ));
  double ro22 = -(r23or22*or33*or33*dx[0]*or11*or22*R[4] + or22*or22*dx[0]*or11);
  double ro23 = r23or22*or33*or33*dx[0]*or11*or22*R[3];
  double ro24 = (alpha3*R[4] /( R[3]*R[3] ) - 2.0*(dx[1] - r12or11*dx[0]) * or22 / ( R[3]*R[3] ) + r23or22*or33*or33*dx[0]*(or11*or22*R[2] + (R[1]*R[4] - R[2]*R[3]) * or11 / ( R[3]*R[3] )) - r23or22*or33*or33*dx[1]*R[4] /( R[3]*R[3] ));
  double ro25 = (-alpha3*or22 + r23or22*or33*or33*dx[1]*or22 - r23or22*or33*or33*dx[0]*or11*or22*R[1]);
  double ro26 = r23or22*2.0*(dx[2] - r23or22*dx[1] + psi*dx[0]) * or33 / ( R[5]*R[5] );
  double co31 = or33*or33*psi;
  double co32 = -or33*or33*r23or22;
  double co33 = or33*or33;
  double ro31 = -or33*or33*dx[0]*or22*(R[1]*R[4] - R[2]*R[3]) / ( R[0]*R[0] );
  double ro32 = or33*or33*dx[0]*or11*or22*R[4];
  double ro33 = - or33*or33*dx[0]*or11*or22*R[3];
  double ro34 = (or33*or33*dx[1]*R[4] /( R[3]*R[3] ) - or33*or33*dx[0]*(or11*or22*R[2] + (R[1]*R[4] - R[2]*R[3]) * or11 / ( R[3]*R[3] )));
  double ro35 = (or33*or33*dx[0]*or11*or22*R[1] - or33*or33*dx[1]*or22);
  double ro36 = - 2.0*(dx[2] - r23or22*dx[1] + psi*dx[0]) * or33 / ( R[5]*R[5] );

// size of dWdddx is 3x6
  dWdX[0][0] =-pm*co11;  dWdX[0][1] =-pm*co12;  dWdX[0][2] =-pm*co13;  dWdX[0][3] = pm*co11;  dWdX[0][4] = pm*co12;  dWdX[0][5] = pm*co13;
  dWdX[1][0] =-pm*co21;  dWdX[1][1] =-pm*co22;  dWdX[1][2] =-pm*co23;  dWdX[1][3] = pm*co21;  dWdX[1][4] = pm*co22;  dWdX[1][5] = pm*co23;
  dWdX[2][0] =-pm*co31;  dWdX[2][1] =-pm*co32;  dWdX[2][2] =-pm*co33;  dWdX[2][3] = pm*co31;  dWdX[2][4] = pm*co32;  dWdX[2][5] = pm*co33;

// size of dWdR is 3x6
  dWdR[0][0] = ro11;  dWdR[0][1] = ro12;  dWdR[0][2] = ro13;  dWdR[0][3] = ro14;  dWdR[0][4] = ro15;  dWdR[0][5] = ro16; 
  dWdR[1][0] = ro21;  dWdR[1][1] = ro22;  dWdR[1][2] = ro23;  dWdR[1][3] = ro24;  dWdR[1][4] = ro25;  dWdR[1][5] = ro26; 
  dWdR[2][0] = ro31;  dWdR[2][1] = ro32;  dWdR[2][2] = ro33;  dWdR[2][3] = ro34;  dWdR[2][4] = ro35;  dWdR[2][5] = ro36; 

}

//------------------------------------------------------------------------------
// Included (MB)
inline
void computeDerivativeOfLocalWeightsLeastSquares(double dx[3], double ddx[3], double *R, double *dR, double *W, double *dW)
{

  if(R[0]*R[3]*R[5] == 0.0) fprintf(stderr, "Going to divide by 0 %e %e %e\n",
         R[0], R[3], R[5]);
  double or11 = 1.0 / R[0];
  double dor11 = -1.0 / ( R[0]*R[0] )*dR[0];
  double or22 = 1.0 / R[3];
  double dor22 = -1.0 / ( R[3]*R[3] )*dR[3];
  double or33 = 1.0 / R[5];
  double dor33 = -1.0 / ( R[5]*R[5] )*dR[5];

  double r12or11 = R[1] * or11;
  double dr12or11 = dR[1] * or11 + R[1] * dor11;
  double r23or22 = R[4] * or22;
  double dr23or22 = dR[4] * or22 + R[4] * dor22;

  double psi = (R[1]*R[4] - R[2]*R[3]) * or11 * or22;
  double dpsi = (dR[1]*R[4] + R[1]*dR[4] - dR[2]*R[3] - R[2]*dR[3]) * or11 * or22 + (R[1]*R[4] - R[2]*R[3]) * dor11 * or22 + (R[1]*R[4] - R[2]*R[3]) * or11 * dor22;

  double alpha1 = dx[0] * or11 * or11;
  double dalpha1 = ddx[0] * or11 * or11 + dx[0] * dor11 * or11 + dx[0] * or11 * dor11;
  double alpha2 = (dx[1] - r12or11*dx[0]) * or22 * or22;
  double dalpha2 = (ddx[1] - dr12or11*dx[0] - r12or11*ddx[0]) * or22 * or22 + (dx[1] - r12or11*dx[0]) * dor22 * or22 + (dx[1] - r12or11*dx[0]) * or22 * dor22;
  double alpha3 = (dx[2] - r23or22*dx[1] + psi*dx[0]) * or33 * or33;
  double dalpha3 = (ddx[2] - dr23or22*dx[1] - r23or22*ddx[1] + dpsi*dx[0] + psi*ddx[0]) * or33 * or33 + (dx[2] - r23or22*dx[1] + psi*dx[0]) * dor33 * or33 + (dx[2] - r23or22*dx[1] + psi*dx[0]) * or33 * dor33;

  W[0] = alpha1 - r12or11*alpha2 + psi*alpha3;
  W[1] = alpha2 - r23or22*alpha3;
  W[2] = alpha3;

  dW[0] = dalpha1 - dr12or11*alpha2 - r12or11*dalpha2 + dpsi*alpha3  + psi*dalpha3;
  dW[1] = dalpha2 - dr23or22*alpha3 - r23or22*dalpha3;
  dW[2] = dalpha3;

}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
void SubDomain::computeGradientsLeastSquares(SVec<double,3> &X, SVec<double,6> &R,
					     SVec<Scalar,dim> &var, SVec<Scalar,dim> &ddx,
					     SVec<Scalar,dim> &ddy, SVec<Scalar,dim> &ddz)  {	//KTC: CANNOT MODIFY

  ddx = (Scalar) 0.0;
  ddy = (Scalar) 0.0;
  ddz = (Scalar) 0.0;

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

	if (sampleMesh) {

		for (int iEdge=0; iEdge<edges.getNumTwoLayersEdges(); ++iEdge) {

			int l = edges.edgesTwoLayersSampleNode[iEdge];

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

	else {
		for (int l=0; l<edges.size(); ++l) {

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

}

//------------------------------------------------------------------------------
// least square gradient involving only nodes of same fluid (multiphase flow and FSI)
// d2d
template<int dim, class Scalar>
void SubDomain::computeGradientsLeastSquares(SVec<double,3> &X,
					     const Vec<int> &fluidId, SVec<double,6> &R,
					     SVec<Scalar,dim> &var, SVec<Scalar,dim> &ddx,
					     SVec<Scalar,dim> &ddy, SVec<Scalar,dim> &ddz,
					     bool linRecFSI, LevelSetStructure *LSS,
					     bool includeSweptNodes)  {

  ddx = (Scalar) 0.0;
  ddy = (Scalar) 0.0;
  ddz = (Scalar) 0.0;

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

	for(int l=0; l<edges.size(); ++l) 
	{
    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

		bool validEdge = true;

		if(fluidId[i] != fluidId[j]) validEdge = false;

		if(LSS)
		{			
			if(LSS->edgeWithSI(l) || LSS->edgeIntersectsStructure(0.0, l)) validEdge = false;
			if(!LSS->isActive(0.0, i) || !LSS->isActive(0.0, j)) validEdge = false;
			if(!includeSweptNodes && (LSS->isSwept(0.0, i) || LSS->isSwept(0.0, j))) validEdge = false;
		}

		if(!validEdge) continue;

    double Wi[3], Wj[3];
    Scalar deltaVar;

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};

		// should be positive for a well posed least square problem
		if(R[i][0] > 0.0 && fabs(R[i][0]*R[i][3]*R[i][5]) > 1.0e-10) 
      computeLocalWeightsLeastSquares(dx, R[i], Wi);
		else
		{ 
      Wi[0] = 0.0;
      Wi[1] = 0.0;
      Wi[2] = 0.0;
    }

    dx[0] = -dx[0]; dx[1] = -dx[1]; dx[2] = -dx[2];

		// should be positive for a well posed least square problem
		if(R[j][0]>0.0 && fabs(R[j][0]*R[j][3]*R[j][5]) > 1.0e-10)
      computeLocalWeightsLeastSquares(dx, R[j], Wj);
		else
		{
      Wj[0] = 0.0;
      Wj[1] = 0.0;
      Wj[2] = 0.0;
    }

		for(int k=0; k<dim; ++k) 
		{
      deltaVar = var[j][k] - var[i][k];

      ddx[i][k] += Wi[0] * deltaVar;
      ddy[i][k] += Wi[1] * deltaVar;
      ddz[i][k] += Wi[2] * deltaVar;
      ddx[j][k] -= Wj[0] * deltaVar;
      ddy[j][k] -= Wj[1] * deltaVar;
      ddz[j][k] -= Wj[2] * deltaVar;
    }
  }

  if(!linRecFSI)
	{
	 	for(int l=0; l<edges.size(); ++l) 
		{
      int i = edgePtr[l][0];
      int j = edgePtr[l][1];

			if( fluidId[i] != fluidId[j] || (LSS && (LSS->edgeIntersectsStructure(0.0,l) || !LSS->isActive(0.0,i) || !LSS->isActive(0.0,j))) ) {
        for (int k=0; k<dim; ++k)
				{
          ddx[i][k] = ddy[i][k] = ddz[i][k] = ddx[j][k] = ddy[j][k] = ddz[j][k] = 0.0;
      }
    }
		}
	}

}

//------------------------------------------------------------------------------
// least square gradient of single variable involving only nodes of same fluid (multiphase flow and FSI)
template<class Scalar>
void SubDomain::computeGradientLeastSquares(SVec<double,3> &X,
                const Vec<int> &fluidId, SVec<double,6> &R,
                Vec<Scalar> &var, Vec<Scalar> &ddx,
                Vec<Scalar> &ddy, Vec<Scalar> &ddz,
                LevelSetStructure *LSS)  {

  ddx = (Scalar) 0.0;
  ddy = (Scalar) 0.0;
  ddz = (Scalar) 0.0;

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

	 /*	 
	 bool isValid = true;

	 if(fluidId[i] != fluidId[j]) isValid = false;
		
	 if(LSS)
	 {
		 if(LSS->edgeIntersectsStructure(0.0, l)) isValid = false;
		 
		 if(!LSS->isActive(0.0, i) || !LSS->isActive(0.0, j)) isValid = false;
	 }
	 
	 if(!isValid) continue;
	 */
    if(higherOrderMF)  continue;

    double Wi[3], Wj[3];
    Scalar deltaVar;

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    if(R[i][0]>0.0 && fabs(R[i][0]*R[i][3]*R[i][5]) > 1.0e-10) // should be positive for a well posed least square problem
      computeLocalWeightsLeastSquares(dx, R[i], Wi);
    else{ // gradient is set to 0.0
      Wi[0] = 0.0;
      Wi[1] = 0.0;
      Wi[2] = 0.0;
    }

    dx[0] = -dx[0]; dx[1] = -dx[1]; dx[2] = -dx[2];
    if(R[j][0]>0.0 && fabs(R[j][0]*R[j][3]*R[j][5]) > 1.0e-10) // should be positive for a well posed least square problem
      computeLocalWeightsLeastSquares(dx, R[j], Wj);
    else{ // gradient is set to 0.0
      Wj[0] = 0.0;
      Wj[1] = 0.0;
      Wj[2] = 0.0;
    }
    deltaVar = var[j] - var[i];

    ddx[i] += Wi[0] * deltaVar;
    ddy[i] += Wi[1] * deltaVar;
    ddz[i] += Wi[2] * deltaVar;
    ddx[j] -= Wj[0] * deltaVar;
    ddy[j] -= Wj[1] * deltaVar;
    ddz[j] -= Wj[2] * deltaVar;
  }
}

//------------------------------------------------------------------------------
// least square gradient involving only nodes of fluid (FSI)
// if Wstar is available, they are also involved in gradient computation
template<int dim, class Scalar>
void SubDomain::computeGradientsLeastSquares(SVec<double,3> &X,
                const Vec<int> &fluidId, SVec<double,6> &R,
                SVec<Scalar,dim> &var, SVec<Scalar,dim> &Wstarij,
				SVec<Scalar,dim> &Wstarji, Vec<int> &countWstarij,
				Vec<int> &countWstarji, SVec<Scalar,dim> &ddx,
                SVec<Scalar,dim> &ddy, SVec<Scalar,dim> &ddz,
                bool linRecFSI, LevelSetStructure *LSS)  {

  ddx = (Scalar) 0.0;
  ddy = (Scalar) 0.0;
  ddz = (Scalar) 0.0;

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];
	 /*
	 bool isValid = true;

	 if(fluidId[i] != fluidId[j]) isValid = false;
		
	 if(LSS)
	 {
		 if(LSS->edgeIntersectsStructure(0.0, l)) isValid = false;
		 
		 if(!LSS->isActive(0.0, i) || !LSS->isActive(0.0, j)) isValid = false;		 
	 }
	 
	 if(!isValid) continue;
	 */
    double Wi[3], Wj[3];
    Scalar deltaVari[dim];
	Scalar deltaVarj[dim];
	bool updatei = false;
	bool updatej = false;

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};

	if (((!LSS)&&(fluidId[i]==fluidId[j]))||(LSS && LSS->isActive(0.0,i) && LSS->isActive(0.0,j) && !LSS->edgeIntersectsStructure(0.0,l))) {
      if(R[i][0]>0.0 && fabs(R[i][0]*R[i][3]*R[i][5]) > 1.0e-10) // should be positive for a well posed least square problem
        computeLocalWeightsLeastSquares(dx, R[i], Wi);
      else{ // gradient is set to 0.0
        Wi[0] = 0.0;
        Wi[1] = 0.0;
        Wi[2] = 0.0;
      }
      dx[0] = -dx[0]; dx[1] = -dx[1]; dx[2] = -dx[2];
      if(R[j][0]>0.0 && fabs(R[j][0]*R[j][3]*R[j][5]) > 1.0e-10) // should be positive for a well posed least square problem
        computeLocalWeightsLeastSquares(dx, R[j], Wj);
      else{ // gradient is set to 0.0
        Wj[0] = 0.0;
        Wj[1] = 0.0;
        Wj[2] = 0.0;
      }
	  updatei = true;
	  updatej = true;
      for (int k=0; k<dim; ++k) {
		deltaVari[k] = var[j][k] - var[i][k];
		deltaVarj[k] = -deltaVari[k];
	  }
	}
	else {
	  if (LSS->isActive(0.0,i))
		if (countWstarij[l]!=0) {
		  LevelSetResult resij = LSS->getLevelSetDataAtEdgeCenter(0.0,l,true);
		  for (int k=0; k<3; ++k)
		    dx[k] = (X[j][k]-X[i][k]);//*(1.0-resij.alpha);
		  if (R[i][0]>0.0 && fabs(R[i][0]*R[i][3]*R[i][5]) > 1.0e-10)
		    computeLocalWeightsLeastSquares(dx, R[i], Wi);
		  else {
		    Wi[0] = 0.0;
		    Wi[0] = 0.0;
		    Wi[0] = 0.0;
		  }
		  updatei = true;
		  for (int k=0; k<dim; ++k) deltaVari[k] = Wstarij[l][k]-var[i][k];
		}
	  if (LSS->isActive(0.0,j))
		if (countWstarji[l]!=0) {
		  LevelSetResult resji = LSS->getLevelSetDataAtEdgeCenter(0.0,l,false);
		  for (int k=0; k<3; ++k)
		    dx[k] = (X[i][k]-X[j][k]);//*(1.0-resji.alpha);
		  if (R[j][0]>0.0 && fabs(R[j][0]*R[j][3]*R[j][5]) > 1.0e-10)
		    computeLocalWeightsLeastSquares(dx, R[j], Wj);
		  else {
		    Wj[0] = 0.0;
		    Wj[0] = 0.0;
		    Wj[0] = 0.0;
		  }
		  updatej = true;
		  for (int k=0; k<dim; ++k) deltaVarj[k] = Wstarji[l][k]-var[j][k];
		}
	}
	if (updatei)
	  for (int k=0; k<dim; ++k) {
        ddx[i][k] += Wi[0] * deltaVari[k];
        ddy[i][k] += Wi[1] * deltaVari[k];
        ddz[i][k] += Wi[2] * deltaVari[k];
	  }
	if (updatej)
	  for (int k=0; k<dim; ++k) {
        ddx[j][k] += Wj[0] * deltaVarj[k];
        ddy[j][k] += Wj[1] * deltaVarj[k];
        ddz[j][k] += Wj[2] * deltaVarj[k];
	  }
  }

//KW: set gradients = 0 for cells near interface.
  if(!linRecFSI)
    for (int l=0; l<edges.size(); ++l) {
      int i = edgePtr[l][0];
      int j = edgePtr[l][1];
      if (fluidId[i]!=fluidId[j] || (LSS && LSS->edgeIntersectsStructure(0.0,l))) {
        for (int k=0; k<dim; ++k)
          ddx[i][k] = ddy[i][k] = ddz[i][k] = ddx[j][k] = ddy[j][k] = ddz[j][k] = 0.0;
      }
    }

}

//------------------------------------------------------------------------------
// Included (YC)
template<int dim, class Scalar>
void SubDomain::computeDerivativeOfGradientsLeastSquares(
               RectangularSparseMat<double,3,dim> *dddxdX,
               RectangularSparseMat<double,3,dim> *dddydX,
               RectangularSparseMat<double,3,dim> *dddzdX,
               RectangularSparseMat<double,6,dim> *dddxdR,
               RectangularSparseMat<double,6,dim> *dddydR,
               RectangularSparseMat<double,6,dim> *dddzdR,
               RectangularSparseMat<double,dim,dim> *dddxdV,
               RectangularSparseMat<double,dim,dim> *dddydV,
               RectangularSparseMat<double,dim,dim> *dddzdV,
               SVec<double,3> &dX,
               SVec<double,6> &dR,
               SVec<double,dim> &dV,
               SVec<Scalar,dim> &dddx,
               SVec<Scalar,dim> &dddy, 
               SVec<Scalar,dim> &dddz)
{

  dddx = (Scalar) 0.0;
  dddy = (Scalar) 0.0;
  dddz = (Scalar) 0.0;

  SVec<Scalar,dim> xxx(dddx), xxx2(dddx); 
  SVec<Scalar,dim> yyy(dddy), yyy2(dddy); 
  SVec<Scalar,dim> zzz(dddz), zzz2(dddz); 

  xxx = (Scalar) 0;    yyy = (Scalar) 0;   zzz = (Scalar) 0;
  xxx2 = (Scalar) 0;   yyy2 = (Scalar) 0;  zzz2 = (Scalar) 0;

  dddxdX->apply(dX, xxx);
  dddydX->apply(dX, yyy);
  dddzdX->apply(dX, zzz);

  dddxdV->apply(dV, xxx2);
  dddydV->apply(dV, yyy2);
  dddzdV->apply(dV, zzz2);

  dddxdR->apply(dR, dddx);
  dddydR->apply(dR, dddy);
  dddzdR->apply(dR, dddz);

  dddx += xxx;    dddy += yyy;    dddz += zzz;
  dddx += xxx2;   dddy += yyy2;   dddz += zzz2;

}

//------------------------------------------------------------------------------
// Included (YC)
template<int dim, class Scalar>
void SubDomain::computeTransposeDerivativeOfGradientsLeastSquares(
               RectangularSparseMat<double,3,dim> *dddxdX,
               RectangularSparseMat<double,3,dim> *dddydX,
               RectangularSparseMat<double,3,dim> *dddzdX,
               RectangularSparseMat<double,6,dim> *dddxdR,
               RectangularSparseMat<double,6,dim> *dddydR,
               RectangularSparseMat<double,6,dim> *dddzdR,
               RectangularSparseMat<double,dim,dim> *dddxdV,
               RectangularSparseMat<double,dim,dim> *dddydV,
               RectangularSparseMat<double,dim,dim> *dddzdV,
               SVec<Scalar,dim> &dddx,
               SVec<Scalar,dim> &dddy, 
               SVec<Scalar,dim> &dddz,
               SVec<double,3> &dX2,
               SVec<double,6> &dR2,
               SVec<double,dim> &dV2)
{

  SVec<Scalar,3> dummyX(dX2);
  SVec<Scalar,6> dummyR(dR2);
  SVec<Scalar,dim> dummyV(dV2);

  dummyX = 0.0;
  dddxdX->applyTranspose(dddx, dummyX);
  dX2 += dummyX;
  dummyX = 0.0;
  dddydX->applyTranspose(dddy, dummyX);
  dX2 += dummyX;
  dummyX = 0.0;
  dddzdX->applyTranspose(dddz, dummyX);
  dX2 += dummyX;

  dddxdR->applyTranspose(dddx, dummyR);
  dR2 += dummyR;
  dddydR->applyTranspose(dddy, dummyR);
  dR2 += dummyR;
  dddzdR->applyTranspose(dddz, dummyR);
  dR2 += dummyR;

// TODO: uncomment below
  dddxdV->applyTranspose(dddx, dummyV);
  dV2 += dummyV;
  dddydV->applyTranspose(dddy, dummyV);
  dV2 += dummyV;
  dddzdV->applyTranspose(dddz, dummyV);
  dV2 += dummyV;

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar>
void SubDomain::computeDerivativeOfGradientsLeastSquares(
               SVec<double,3> &X, SVec<double,3> &dX,
							 SVec<double,6> &R, SVec<double,6> &dR,
							 SVec<Scalar,dim> &var, SVec<Scalar,dim> &dvar, SVec<Scalar,dim> &dddx,
							 SVec<Scalar,dim> &dddy, SVec<Scalar,dim> &dddz)
{

  dddx = (Scalar) 0.0;
  dddy = (Scalar) 0.0;
  dddz = (Scalar) 0.0;

  SVec<Scalar,dim> xxx(dddx), dddx2(dddx), uux(dddx);
  SVec<Scalar,dim> yyy(dddy), dddy2(dddy), uuy(dddy);
  SVec<Scalar,dim> zzz(dddz), dddz2(dddz), uuz(dddz);

  dddx = (Scalar) 0.0;
  dddy = (Scalar) 0.0;
  dddz = (Scalar) 0.0;

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    double Wi[3], Wj[3], deltaVar;
    double dWi[3]={0}, dWj[3]={0}, dDeltaVar;
    double dWi2[3]={0}, dWj2[3]={0};
    double dWidX[3][6] = {0}, dWidR[3][6] = {0};

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    double ddx[3] = {dX[j][0] - dX[i][0], dX[j][1] - dX[i][1], dX[j][2] - dX[i][2]};

    computeDerivativeOfLocalWeightsLeastSquares(dx, ddx, R[i], dR[i], Wi, dWi);

    dx[0] = -dx[0]; dx[1] = -dx[1]; dx[2] = -dx[2];
    ddx[0] = -ddx[0]; ddx[1] = -ddx[1]; ddx[2] = -ddx[2];

    computeDerivativeOfLocalWeightsLeastSquares(dx, ddx, R[j], dR[j], Wj, dWj);

    for (int k=0; k<dim; ++k) {
      deltaVar = var[j][k] - var[i][k];
      dDeltaVar = dvar[j][k] - dvar[i][k];
      dddx[i][k] += (dWi[0] * deltaVar + Wi[0] * dDeltaVar);
      dddy[i][k] += (dWi[1] * deltaVar + Wi[1] * dDeltaVar);
      dddz[i][k] += (dWi[2] * deltaVar + Wi[2] * dDeltaVar);
      dddx[j][k] -= (dWj[0] * deltaVar + Wj[0] * dDeltaVar);
      dddy[j][k] -= (dWj[1] * deltaVar + Wj[1] * dDeltaVar);
      dddz[j][k] -= (dWj[2] * deltaVar + Wj[2] * dDeltaVar);
    }

  }
}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim, class Scalar>
void SubDomain::computeDerivativeOperatorsOfGradientsLeastSquares(SVec<double,3> &X, SVec<double,6> &R, SVec<Scalar,dim> &var, 
                                                                  RectangularSparseMat<double,3,dim> &dddxdX,
                                                                  RectangularSparseMat<double,3,dim> &dddydX,
                                                                  RectangularSparseMat<double,3,dim> &dddzdX,
                                                                  RectangularSparseMat<double,6,dim> &dddxdR,
                                                                  RectangularSparseMat<double,6,dim> &dddydR,
                                                                  RectangularSparseMat<double,6,dim> &dddzdR,
                                                                  RectangularSparseMat<double,dim,dim> &dddxdV,
                                                                  RectangularSparseMat<double,dim,dim> &dddydV,
                                                                  RectangularSparseMat<double,dim,dim> &dddzdV)
   
{

  double dW0dWi[2][3] = {0}, dW0dWj[2][3] = {0}, dW1dWi[2][3] = {0}, dW1dWj[2][3] = {0}, dW2dWi[2][3] = {0}, dW2dWj[2][3] = {0};
  dW0dWi[0][0] = 1;  dW0dWj[1][0] = 1;
  dW1dWi[0][1] = 1;  dW1dWj[1][1] = 1;
  dW2dWi[0][2] = 1;  dW2dWj[1][2] = 1;

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    double Wi[3], Wj[3], deltaVar;
    double dWi[3]={0}, dWj[3]={0}, dWi2[3], dWj2[3];

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    double dWidX[3][6] = {0}, dWidR[3][6] = {0};
    compute_dWdXAnddWdR( 1, dx, R[i], Wi, dWidX, dWidR);
    dx[0] = -dx[0]; dx[1] = -dx[1]; dx[2] = -dx[2];
    double dWjdX[3][6] = {0}, dWjdR[3][6] = {0};
    compute_dWdXAnddWdR(-1, dx, R[j], Wj, dWjdX, dWjdR);

    double dddxdVarray[2*dim][2*dim] = {0}, dddydVarray[2*dim][2*dim] = {0}, dddzdVarray[2*dim][2*dim] = {0};

    for (int k=0; k<dim; ++k) {
      dddxdVarray[k][k] = -Wi[0];      dddxdVarray[k][dim+k] = Wi[0];
      dddydVarray[k][k] = -Wi[1];      dddydVarray[k][dim+k] = Wi[1];
      dddzdVarray[k][k] = -Wi[2];      dddzdVarray[k][dim+k] = Wi[2];
      dddxdVarray[dim+k][k] = Wj[0];   dddxdVarray[dim+k][dim+k] = -Wj[0];
      dddydVarray[dim+k][k] = Wj[1];   dddydVarray[dim+k][dim+k] = -Wj[1];
      dddzdVarray[dim+k][k] = Wj[2];   dddzdVarray[dim+k][dim+k] = -Wj[2];
    }


    double ddddW[2*dim][2] = {0};

    for (int k=0; k<dim; ++k) {
      deltaVar = var[j][k] - var[i][k];

      ddddW[k][0] = deltaVar;
      ddddW[dim+k][1] = -deltaVar;
    }
    
    double dW0dX[2][6] = {0}, dW0dR[2][12] = {0};
    double dW1dX[2][6] = {0}, dW1dR[2][12] = {0};
    double dW2dX[2][6] = {0}, dW2dR[2][12] = {0};
    double dQxdX[2*dim][6] = {0}, dQxdR[2*dim][12] = {0};
    double dQydX[2*dim][6] = {0}, dQydR[2*dim][12] = {0};
    double dQzdX[2*dim][6] = {0}, dQzdR[2*dim][12] = {0};

    for(int k=0; k<2; ++k) 
      for(int m=0; m<6; ++m) 
        for(int n=0; n<3; ++n) {
          dW0dX[k][m] += dW0dWi[k][n] * dWidX[n][m] + dW0dWj[k][n] * dWjdX[n][m];
          dW1dX[k][m] += dW1dWi[k][n] * dWidX[n][m] + dW1dWj[k][n] * dWjdX[n][m];
          dW2dX[k][m] += dW2dWi[k][n] * dWidX[n][m] + dW2dWj[k][n] * dWjdX[n][m];
          dW0dR[k][m] += dW0dWi[k][n] * dWidR[n][m]; 
          dW1dR[k][m] += dW1dWi[k][n] * dWidR[n][m]; 
          dW2dR[k][m] += dW2dWi[k][n] * dWidR[n][m]; 
          dW0dR[k][m+6] += dW0dWj[k][n] * dWjdR[n][m];
          dW1dR[k][m+6] += dW1dWj[k][n] * dWjdR[n][m];
          dW2dR[k][m+6] += dW2dWj[k][n] * dWjdR[n][m];
        }
      
    for(int k=0; k<2*dim; ++k) 
      for(int m=0; m<6; ++m)
        for(int n=0; n<2; ++n) {
          dQxdX[k][m] += ddddW[k][n] * dW0dX[n][m];
          dQydX[k][m] += ddddW[k][n] * dW1dX[n][m];
          dQzdX[k][m] += ddddW[k][n] * dW2dX[n][m];
          dQxdR[k][m] += ddddW[k][n] * dW0dR[n][m];
          dQydR[k][m] += ddddW[k][n] * dW1dR[n][m];
          dQzdR[k][m] += ddddW[k][n] * dW2dR[n][m];
          dQxdR[k][m+6] += ddddW[k][n] * dW0dR[n][m+6];
          dQydR[k][m+6] += ddddW[k][n] * dW1dR[n][m+6];
          dQzdR[k][m+6] += ddddW[k][n] * dW2dR[n][m+6];
        }

    int ndList[2] = {i, j};
    dddxdV.addContrib(2,ndList,dddxdVarray[0]);
    dddydV.addContrib(2,ndList,dddydVarray[0]);
    dddzdV.addContrib(2,ndList,dddzdVarray[0]);
    dddxdX.addContrib(2,ndList,dQxdX[0]);
    dddydX.addContrib(2,ndList,dQydX[0]);
    dddzdX.addContrib(2,ndList,dQzdX[0]);
    dddxdR.addContrib(2,ndList,dQxdR[0]); 
    dddydR.addContrib(2,ndList,dQydR[0]); 
    dddzdR.addContrib(2,ndList,dQzdR[0]); 
     
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
void SubDomain::computeGradientsGalerkin(Vec<double> &ctrlVol, 
					 SVec<double,3> &wii, SVec<double,3> &wij, SVec<double,3> &wji, 
					 SVec<Scalar,dim> &var,
					 SVec<Scalar,dim> &ddx, SVec<Scalar,dim> &ddy, SVec<Scalar,dim> &ddz)
{

  int i, j, k, l;

  for (i=0; i<var.size(); ++i)  {
    for (k=0; k<dim; ++k)  {
      ddx[i][k] = var[i][k] * wii[i][0];
      ddy[i][k] = var[i][k] * wii[i][1];
      ddz[i][k] = var[i][k] * wii[i][2];
    }
  }

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (l=0; l<edges.size(); ++l) {

    //if (!edgeFlag[l]) continue;
    i = edgePtr[l][0];
    j = edgePtr[l][1];

    for (k=0; k<dim; ++k)  {

      ddx[i][k] += wij[l][0] * var[j][k];
      ddx[j][k] += wji[l][0] * var[i][k];

      ddy[i][k] += wij[l][1] * var[j][k];
      ddy[j][k] += wji[l][1] * var[i][k];

      ddz[i][k] += wij[l][2] * var[j][k];
      ddz[j][k] += wji[l][2] * var[i][k];

    }

  }

  for (i=0; i<var.size(); ++i)  {

    double coef = 1.0 / (4.0*ctrlVol[i]);

    for (k=0; k<dim; ++k) {
      ddx[i][k] *= coef;
      ddy[i][k] *= coef;
      ddz[i][k] *= coef;
    }
  }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar>
void SubDomain::computeDerivativeOfGradientsGalerkin(Vec<double> &ctrlVol, Vec<double> &dCtrlVol,
					 SVec<double,3> &wii, SVec<double,3> &wij, SVec<double,3> &wji,
					 SVec<double,3> &dwii, SVec<double,3> &dwij, SVec<double,3> &dwji,
					 SVec<Scalar,dim> &var, SVec<Scalar,dim> &dvar, SVec<Scalar,dim> &ddx,
                     SVec<Scalar,dim> &ddy, SVec<Scalar,dim> &ddz, SVec<Scalar,dim> &dddx,
                     SVec<Scalar,dim> &dddy, SVec<Scalar,dim> &dddz)
{

  int i, j, k, l;

  ddx=0.0;
  ddy=0.0;
  ddz=0.0;

  for (i=0; i<var.size(); ++i)  {
    for (k=0; k<dim; ++k)  {
      ddx[i][k] = var[i][k] * wii[i][0];
      ddy[i][k] = var[i][k] * wii[i][1];
      ddz[i][k] = var[i][k] * wii[i][2];

      dddx[i][k] = var[i][k] * dwii[i][0] + dvar[i][k] * wii[i][0];
      dddy[i][k] = var[i][k] * dwii[i][1] + dvar[i][k] * wii[i][1];
      dddz[i][k] = var[i][k] * dwii[i][2] + dvar[i][k] * wii[i][2];
    }
  }

  int (*edgePtr)[2] = edges.getPtr();

  for (l=0; l<edges.size(); ++l) {

    i = edgePtr[l][0];
    j = edgePtr[l][1];

    for (k=0; k<dim; ++k)  {

      ddx[i][k] += wij[l][0] * var[j][k];
      ddx[j][k] += wji[l][0] * var[i][k];

      ddy[i][k] += wij[l][1] * var[j][k];
      ddy[j][k] += wji[l][1] * var[i][k];

      ddz[i][k] += wij[l][2] * var[j][k];
      ddz[j][k] += wji[l][2] * var[i][k];

      dddx[i][k] += dwij[l][0] * var[j][k] + wij[l][0] * dvar[j][k];
      dddx[j][k] += dwji[l][0] * var[i][k] + wji[l][0] * dvar[i][k];

      dddy[i][k] += dwij[l][1] * var[j][k] + wij[l][1] * dvar[j][k];
      dddy[j][k] += dwji[l][1] * var[i][k] + wji[l][1] * dvar[i][k];

      dddz[i][k] += dwij[l][2] * var[j][k] + wij[l][2] * dvar[j][k];
      dddz[j][k] += dwji[l][2] * var[i][k] + wji[l][2] * dvar[i][k];

    }

  }

  for (i=0; i<var.size(); ++i)  {

    double coef = 1.0 / (4.0*ctrlVol[i]);

    double dcoef = -1.0 / (4.0*ctrlVol[i]*ctrlVol[i]) * dCtrlVol[i];

    for (k=0; k<dim; ++k) {
      dddx[i][k] *= coef;
      dddx[i][k] += ddx[i][k]*dcoef;
      dddy[i][k] *= coef;
      dddy[i][k] += ddy[i][k]*dcoef;
      dddz[i][k] *= coef;
      dddz[i][k] += ddz[i][k]*dcoef;
    }

  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
void SubDomain::computeGradientsGalerkinT(Vec<double> &ctrlVol,
                SVec<double,3> &wii, SVec<double,3> &wij, SVec<double,3> &wji,
                SVec<Scalar,dim> &var, SVec<Scalar,dim> &var1,
                SVec<Scalar,dim> &var2,SVec<Scalar,dim> &ddx,
                SVec<Scalar,dim> &ddy, SVec<Scalar,dim> &ddz)  {

  int i, j, k, l;

  for (i=0; i<var.size(); ++i)  {
    for (k=0; k<dim; ++k)  {
      ddx[i][k] = var[i][k] * wii[i][0] / (4.0*ctrlVol[i]);
      ddy[i][k] = var1[i][k] * wii[i][1] / (4.0*ctrlVol[i]);
      ddz[i][k] = var2[i][k] * wii[i][2] / (4.0*ctrlVol[i]);
    }
  }

  //bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (l=0; l<edges.size(); ++l) {
    i = edgePtr[l][0];
    j = edgePtr[l][1];

    for (k=0; k<dim; ++k)  {

      ddx[i][k] += wji[l][0] * var[j][k] / (4.0*ctrlVol[j]);
      ddx[j][k] += wij[l][0] * var[i][k] / (4.0*ctrlVol[i]);

      ddy[i][k] += wji[l][1] * var1[j][k] / (4.0*ctrlVol[j]);
      ddy[j][k] += wij[l][1] * var1[i][k] / (4.0*ctrlVol[i]);

      ddz[i][k] += wji[l][2] * var2[j][k] / (4.0*ctrlVol[j]);
      ddz[j][k] += wij[l][2] * var2[i][k] / (4.0*ctrlVol[i]);


    }
  }

}


//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeMinMaxStencilValues(SVec<double,dim> &V, SVec<double,dim> &Vmin,
														 SVec<double,dim> &Vmax, LevelSetStructure *LSS)
{
  Vmin = V;
  Vmax = V;

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

	for(int l=0; l<edges.size(); ++l) 
	{
    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

		if(LSS)
		{
			if(!LSS->isActive(0.0, i) || !LSS->isActive(0.0, j)) continue;
		}

		for(int k=0; k<dim; ++k) 
		{
      Vmin[i][k] = min(Vmin[i][k], V[j][k]);
      Vmax[i][k] = max(Vmax[i][k], V[j][k]);
      Vmin[j][k] = min(Vmin[j][k], V[i][k]);
      Vmax[j][k] = max(Vmax[j][k], V[i][k]);
    }

  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeDerivativeOfMinMaxStencilValues(SVec<double,dim> &V, SVec<double,dim> &dV, SVec<double,dim> &Vmin, SVec<double,dim> &dVmin,
					   SVec<double,dim> &Vmax, SVec<double,dim> &dVmax)
{

  Vmin = V;
  Vmax = V;
  dVmin = dV;
  dVmax = dV;

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    for (int k=0; k<dim; ++k) {

      if (min(Vmin[i][k], V[j][k]) == V[j][k]) {
        Vmin[i][k] = V[j][k];
        dVmin[i][k] = dV[j][k];
      }

      if (max(Vmax[i][k], V[j][k]) == V[j][k]) {
        Vmax[i][k] = V[j][k];
        dVmax[i][k] = dV[j][k];
      }

      if (min(Vmin[j][k], V[i][k]) == V[i][k]) {
        Vmin[j][k] = V[i][k];
        dVmin[j][k] = dV[i][k];
      }

      if (max(Vmax[j][k], V[i][k]) == V[i][k]) {
        Vmax[j][k] = V[i][k];
        dVmax[j][k] = dV[i][k];
      }

    }

  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeMultiDimLimiter(RecLimiter *recFcn, SVec<double,3> &X,
				       Vec<double> &ctrlVol, SVec<double,dim> &V,
				       SVec<double,dim> &dVdx, SVec<double,dim> &dVdy,
				       SVec<double,dim> &dVdz, SVec<double,dim> &Vmin,
													SVec<double,dim> &Vmax, SVec<double,dim> &phi,
													LevelSetStructure *LSS)
{

  double ddVij[dim], ddVji[dim], Vi[dim], Vj[dim];

  phi = 1.0;
  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

	for(int l=0; l<edges.size(); ++l) 
	{
    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

		if(LSS)
		{
			if(!LSS->isActive(0.0, i) || !LSS->isActive(0.0, j)) continue;
		}

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};

		for(int k=0; k<dim; ++k) 
		{
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
    }

    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);

    recFcn->computeLimiter(Vmax[i], Vmin[i], V[i], Vi, ctrlVol[i],
			   Vmax[j], Vmin[j], V[j], Vj, ctrlVol[j], phi[i], phi[j]);

  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeDerivativeOfMultiDimLimiter(RecLimiter *recFcn, SVec<double,3> &X, SVec<double,3> &dX,
				       Vec<double> &ctrlVol, Vec<double> &dCtrlVol, SVec<double,dim> &V, SVec<double,dim> &dV,
				       SVec<double,dim> &dVdx, SVec<double,dim> &dVdy, SVec<double,dim> &dVdz,
				       SVec<double,dim> &ddVdx, SVec<double,dim> &ddVdy, SVec<double,dim> &ddVdz,
                       SVec<double,dim> &Vmin, SVec<double,dim> &dVmin, SVec<double,dim> &Vmax,
                       SVec<double,dim> &dVmax, SVec<double,dim> &phi, SVec<double,dim> &dphi)
{

  double ddVij[dim], ddVji[dim], Vi[dim], Vj[dim];

  double dddVij[dim], dddVji[dim], dVi[dim], dVj[dim];

  phi = 1.0;
  dphi = 0.0;
  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    double ddx[3] = {dX[j][0] - dX[i][0], dX[j][1] - dX[i][1], dX[j][2] - dX[i][2]};
    for (int k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
      dddVij[k] = ddx[0]*dVdx[i][k] + dx[0]*ddVdx[i][k] + ddx[1]*dVdy[i][k] + dx[1]*ddVdy[i][k] + ddx[2]*dVdz[i][k] + dx[2]*ddVdz[i][k];
      dddVji[k] = ddx[0]*dVdx[j][k] + dx[0]*ddVdx[j][k] + ddx[1]*dVdy[j][k] + dx[1]*ddVdy[j][k] + ddx[2]*dVdz[j][k] + dx[2]*ddVdz[j][k];
    }

    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);

    recFcn->computeDerivative(V[i], dV[i], ddVij, dddVij, V[j], dV[j], ddVji, dddVji, dVi, dVj);

    recFcn->computeDerivativeOfLimiter(Vmax[i], dVmax[i], Vmin[i], dVmin[i], V[i], dV[i], Vi, dVi, ctrlVol[i], dCtrlVol[i],
			   Vmax[j], dVmax[j], Vmin[j], dVmin[j], V[j], dV[j], Vj, dVj, ctrlVol[j], dCtrlVol[j], phi[i], dphi[i], phi[j], dphi[j]);

  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computePressureSensor(SVec<double,3>& X, SVec<double,dim>& V,
				      SVec<double,dim>& dVdx, SVec<double,dim>& dVdy,
				      SVec<double,dim>& dVdz, SVec<double,3>& sensor)
{

  sensor = 0.0;

  bool* edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {
    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    double dx[3];
    dx[0] = X[j][0] - X[i][0];
    dx[1] = X[j][1] - X[i][1];
    dx[2] = X[j][2] - X[i][2];

    double dpi = dx[0]*dVdx[i][4] + dx[1]*dVdy[i][4] + dx[2]*dVdz[i][4];
    double dpj = dx[0]*dVdx[j][4] + dx[1]*dVdy[j][4] + dx[2]*dVdz[j][4];

    double s0 = fabs(dpj - dpi);
    double s1 = fabs(dpj) + fabs(dpi);
    double s2 = fabs(V[i][4]) + fabs(V[j][4]);

    sensor[i][0] += s0;
    sensor[i][1] += s1;
    sensor[i][2] += s2;
    sensor[j][0] += s0;
    sensor[j][1] += s1;
    sensor[j][2] += s2;
  }

}

//------------------------------------------------------------------------------

/*
@INPROCEEDINGS{dervieux-85,
  author = "Dervieux, A.",
  title = "Steady {E}uler simulations using unstructured meshes",
  booktitle = "Proceedings of the VKI Lectures Series 1985-04",
  series = "16th Computational Fluid Dynamics",
  month = mar,
  year = 1985,
  address = "von Karman Institute, Brussels, Belgium",
}
@INBOOK{barth-92,
  author = "Barth, T. J.",
  title = "Aspects of Unstructured Grids and Finite-Volume Solvers
           for the {E}uler and {N}avier-{S}tokes Equations",
  series = "Special Course on Unstructured Grid Methods for Advection
            Dominated Flows",
  publisher = "AGARD R-787",
  month = may,
  year = 1992,
  chapter = "6",
}
@ARTICLE{vanleer-79,
  author = "van Leer, B.",
  title = "Towards the Ultimate Conservative Difference Scheme. {V}.
           {A} Second-Order Sequel to {G}odunov's Method",
  journal = jcp,
  year = 1979,
  volume = 32,
  pages = "101--136",
}
*/
template<int dim>
int SubDomain::computeFiniteVolumeTerm(Vec<double> &irey, FluxFcn** fluxFcn, RecFcn* recFcn,
                                       BcData<dim>& bcData, GeoState& geoState,
                                       SVec<double,3>& X, SVec<double,dim>& V,
                                       NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
                                       SVec<double,dim>& fluxes, SVec<int,2>& tag,
                                       int failsafe, int rshift)
{


  int ierr = edges.computeFiniteVolumeTerm(locToGlobNodeMap, irey, fluxFcn, recFcn, elems, geoState,
                                           X, V, ngrad, egrad, fluxes, tag, failsafe, rshift);

  faces.computeFiniteVolumeTerm(fluxFcn, bcData, geoState, V, fluxes);

  return(ierr);

}

//------------------------------------------------------------------------------

template<int dim>
int SubDomain::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann, 
                                       Vec<double> &irey, FluxFcn** fluxFcn, RecFcn* recFcn,
                                       BcData<dim>& bcData, GeoState& geoState,
				       SVec<double,3>& X, SVec<double,dim>& V,
				       NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
				       SVec<double,dim>& fluxes, SVec<int,2>& tag,
                                       int failsafe, int rshift)
{

  int ierr = edges.computeFiniteVolumeTerm(locToGlobNodeMap, irey, fluxFcn, recFcn, elems, geoState,
                                           X, V, ngrad, egrad, fluxes, tag, failsafe, rshift);

  faces.computeFiniteVolumeTerm(riemann, fluxFcn, bcData, geoState, V, fluxes);

  return(ierr);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeDerivativeOfFiniteVolumeTerm(Vec<double> &irey, Vec<double> &dIrey, FluxFcn** fluxFcn, RecFcn* recFcn,
                                                    BcData<dim>& bcData, GeoState& geoState,
                                                    SVec<double,3>& X, SVec<double,3>& dX, SVec<double,dim>& V, SVec<double,dim>& dV,
                                                    NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad, double dMach,
                                                    SVec<double,dim>& dFluxes)
{
  edges.computeDerivativeOfFiniteVolumeTerm(irey, dIrey, fluxFcn, recFcn, elems, geoState, X, dX, V, dV, ngrad, egrad, dMach, dFluxes);
  faces.computeDerivativeOfFiniteVolumeTerm(fluxFcn, bcData, geoState, V, dFluxes);
}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void SubDomain::computeDerivativeOfFiniteVolumeTerm(
                                        RectangularSparseMat<double,dim,dim> *dFluxdddx,
                                        RectangularSparseMat<double,dim,dim> *dFluxdddy,
                                        RectangularSparseMat<double,dim,dim> *dFluxdddz,
                                        RectangularSparseMat<double,3,dim> *dFluxdEdgeNorm,
                                        RectangularSparseMat<double,3,dim> *dFluxdX,
                                        RectangularSparseMat<double,3,dim> *dFluxdFaceNormal,
                                        RectangularSparseMat<double,1,dim> *dFluxdFaceNormalVel,
                                        RectangularSparseMat<double,dim,dim> *dFluxdUb,
                                        BcData<dim>& bcData, GeoState& geoState,
                                        SVec<double,3>& dX, 
                                        NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
                                        SVec<double,dim>& dddx,
                                        SVec<double,dim>& dddy,
                                        SVec<double,dim>& dddz,
                                        Vec<Vec3D>& dNormal,
                                        Vec<Vec3D>& dn,
                                        Vec<double>& dndot, 
                                        SVec<double,dim>& dFluxes)
{
  edges.computeDerivativeOfFiniteVolumeTerm(dFluxdddx, dFluxdddy, dFluxdddz, dFluxdX, dFluxdEdgeNorm,
                                            elems, geoState, dX, ngrad, egrad, dddx, dddy, dddz, dNormal, dFluxes);
  faces.computeDerivativeOfFiniteVolumeTerm(dFluxdFaceNormal, dFluxdFaceNormalVel, dFluxdUb, 
                                            bcData, geoState, dn, dndot, dFluxes);
}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void SubDomain::computeTransposeDerivativeOfFiniteVolumeTerm(
                                        RectangularSparseMat<double,dim,dim> *dFluxdddx,
                                        RectangularSparseMat<double,dim,dim> *dFluxdddy,
                                        RectangularSparseMat<double,dim,dim> *dFluxdddz,
                                        RectangularSparseMat<double,3,dim> *dFluxdEdgeNorm,
                                        RectangularSparseMat<double,3,dim> *dFluxdX,
                                        RectangularSparseMat<double,3,dim> *dFluxdFaceNormal,
                                        RectangularSparseMat<double,1,dim> *dFluxdFaceNormalVel,
                                        RectangularSparseMat<double,dim,dim> *dFluxdUb,
                                        BcData<dim>& bcData, GeoState& geoState,
                                        SVec<double,dim>& dFluxes, 
                                        NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad, 
                                        SVec<double,3>& dX2,
                                        SVec<double,dim>& dddx2,
                                        SVec<double,dim>& dddy2,
                                        SVec<double,dim>& dddz2,
                                        Vec<Vec3D>& dEdgeNormal2,
                                        Vec<Vec3D>& dFaceNormal2,
                                        Vec<double>& dFaceNormalVel2)
{

  faces.computeTransposeDerivativeOfFiniteVolumeTerm(dFluxdFaceNormal, dFluxdFaceNormalVel, dFluxdUb, bcData, geoState, dFluxes, dFaceNormal2, dFaceNormalVel2);

  edges.computeTransposeDerivativeOfFiniteVolumeTerm(dFluxdddx, dFluxdddy, dFluxdddz, dFluxdX, dFluxdEdgeNorm,
                                                     dFluxes, ngrad, egrad, elems, geoState, dX2, dddx2, dddy2, dddz2, dEdgeNormal2);

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void SubDomain::computeDerivativeOperatorsOfFiniteVolumeTerm(Vec<double> &irey, Vec<double> &dIrey, FluxFcn** fluxFcn, RecFcn* recFcn,
                                        BcData<dim>& bcData, GeoState& geoState,
                                        SVec<double,3>& X, SVec<double,dim>& V,
                                        NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad, double dMach,
                                        RectangularSparseMat<double,3,dim>& dFluxdEdgeNorm,
                                        RectangularSparseMat<double,3,dim> &dFluxdX,
                                        RectangularSparseMat<double,dim,dim> &dFluxdddx,
                                        RectangularSparseMat<double,dim,dim> &dFluxdddy,
                                        RectangularSparseMat<double,dim,dim> &dFluxdddz,
                                        RectangularSparseMat<double,3,dim> &dFluxdFaceNormal,
                                        RectangularSparseMat<double,1,dim> &dFluxdFaceNormalVel,
                                        RectangularSparseMat<double,dim,dim> &dFluxdUb)
{

  edges.computeDerivativeOperatorsOfFiniteVolumeTerm(irey, dIrey, fluxFcn, recFcn, elems, geoState, X, V, ngrad, egrad, dMach, 
                                                     dFluxdEdgeNorm, dFluxdX, dFluxdddx,dFluxdddy, dFluxdddz);

  faces.computeDerivativeOperatorsOfFiniteVolumeTerm(fluxFcn, bcData, geoState, V, dFluxdFaceNormal, dFluxdFaceNormalVel, dFluxdUb);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeDerivativeOfFiniteVolumeTerm(FluxFcn** fluxFcn, RecFcn* recFcn,
						    BcData<dim>& bcData, GeoState& geoState,
						    SVec<double,3>& X, LevelSetStructure &LSS,
						    bool linRecAtInterface, bool viscSecOrder,
						    Vec<int> &fluidId,
						    ExactRiemannSolver<dim>& riemann,
						    int Nriemann,
						    NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
						    double dMach,
						    SVec<double,dim>& V,
						    SVec<double,dim>& dFluxes)
{


  edges.computeDerivativeOfFiniteVolumeTerm(fluxFcn, recFcn, geoState, X, LSS, 
					    linRecAtInterface, fluidId, riemann,
					    Nriemann, ngrad, egrad, dMach, V, dFluxes);

  faces.computeDerivativeOfFiniteVolumeTerm(fluxFcn, bcData, geoState, V, dFluxes);

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
int SubDomain::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                       FluxFcn** fluxFcn, RecFcn* recFcn,
                                       BcData<dim>& bcData, GeoState& geoState,
                                       SVec<double,3>& X, SVec<double,dim>& V,
                                       Vec<int> &fluidId, FluidSelector &fluidSelector,
                                       NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
				       SVec<double,dimLS>& phi,
                                       NodalGrad<dimLS>& ngradLS,
				       EdgeGrad<dimLS>* egradLS,
                                       SVec<double,dim>& fluxes, int it,
                                       SVec<int,2>& tag, int failsafe, int rshift)
{

  int ierr = edges.computeFiniteVolumeTerm(riemann, locToGlobNodeMap, fluxFcn,
                                           recFcn, elems, geoState, X, V, fluidId, fluidSelector,
                                           ngrad, egrad, phi, ngradLS, egradLS, fluxes, it,
                                           tag, failsafe, rshift);

  faces.computeFiniteVolumeTerm(fluxFcn, bcData, geoState, V, fluidId, fluxes);

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
int SubDomain::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                       FluxFcn** fluxFcn, RecFcn* recFcn,
                                       BcData<dim>& bcData, GeoState& geoState,
                                       SVec<double,3>& X, SVec<double,dim>& V,
                                       SVec<double,dim>& Wstarij, SVec<double,dim>& Wstarji, 
                                       LevelSetStructure& LSS, bool linRecAtInterface, 
                                       Vec<int> &fluidId, int Nriemann,
                                       FluidSelector &fluidSelector,
                                       NodalGrad<dim>& ngrad, 
				       EdgeGrad<dim>* egrad,
				       SVec<double,dimLS>& phi,
                                       NodalGrad<dimLS>& ngradLS,
				       EdgeGrad<dimLS>* egradLS,
                                       SVec<double,dim>& fluxes, int it,
                                       SVec<int,2>& tag, int failsafe, int rshift)
{

  int ierr = edges.computeFiniteVolumeTerm(riemann, locToGlobNodeMap, fluxFcn,
                                           recFcn, elems, geoState, X, V, Wstarij, Wstarji, LSS, linRecAtInterface,
                                           fluidId, Nriemann, fluidSelector,
                                           ngrad, egrad, phi, ngradLS, egradLS, fluxes, it,
                                           tag, failsafe, rshift);

  faces.computeFiniteVolumeTerm(fluxFcn, bcData, geoState, V, fluidId, fluxes, &LSS);

  return ierr;

}

//------------------------------------------------------------------------------
//d2d embedded structure
template<int dim>
int SubDomain::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                       FluxFcn** fluxFcn, RecFcn* recFcn,
                                       BcData<dim>& bcData, GeoState& geoState,
                                       SVec<double,3>& X, SVec<double,dim>& V,
                                       SVec<double,dim>& Wstarij, SVec<double,dim>& Wstarji,
                                       SVec<double,dim>& Wext, LevelSetStructure &LSS, 
													bool linRecAtInterface, Vec<int> &fluidId,
                                       int Nriemann, 
													NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
                                       SVec<double,dim>& fluxes, 
													int it, SVec<int,2>& tag, 
													int failsafe, int rshift, bool externalSI)
{

	int ierr;

	if(externalSI) 
	{		
		ierr = edges.computeFiniteVolumeTerm(riemann, locToGlobNodeMap, fluxFcn,
														 recFcn, elems, geoState, X, V, Wstarij, Wstarji, Wext, LSS, 
														 fluidId, Nriemann, ngrad, egrad, fluxes, it,
														 tag, failsafe, rshift);
	}
	else
		ierr = edges.computeFiniteVolumeTerm(riemann, locToGlobNodeMap, fluxFcn,
                                           recFcn, elems, geoState, X, V, Wstarij, Wstarji, LSS, 
                                           linRecAtInterface, fluidId, Nriemann, ngrad, egrad, fluxes, it,
                                          tag, failsafe, rshift);

  faces.computeFiniteVolumeTerm(fluxFcn, bcData, geoState, V, fluidId, fluxes, &LSS);

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim>
int SubDomain::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                       FluxFcn** fluxFcn, RecFcn* recFcn,
                                       BcData<dim>& bcData, GeoState& geoState,
                                       SVec<double,3>& X, SVec<double,dim>& V,
                                       SVec<double,dim>& Wstarij, SVec<double,dim>& Wstarji,
									   Vec<int>& countWstarij, Vec<int>& countWstarji,
                                       LevelSetStructure &LSS, bool linRecAtInterface, 
									   Vec<int> &fluidId, int Nriemann,
									   double dt, double alpha, 
									   NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
                                       SVec<double,dim>& fluxes, int it,
                                       SVec<int,2>& tag, int failsafe, int rshift) 
{

  V6NodeData (*v6data)[2];
  v6data = 0;
  findEdgeTetrahedra(X, v6data);
  int ierr = edges.computeFiniteVolumeTerm(riemann, locToGlobNodeMap, fluxFcn,
                                           recFcn, elems, geoState, X, V, Wstarij, Wstarji, 
										   countWstarij, countWstarji, LSS, linRecAtInterface, 
										   fluidId, Nriemann, dt, alpha, ngrad,
  									       egrad, fluxes, it, tag, failsafe, rshift, v6data); 
  faces.computeFiniteVolumeTerm(fluxFcn, bcData, geoState, V, fluidId, fluxes, &LSS);

  delete [] v6data;

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void SubDomain::computeFiniteVolumeTermLS(FluxFcn** fluxFcn, RecFcn* recFcn, RecFcn* recFcnLS,
                                        BcData<dim>& bcData, GeoState& geoState,
                                        SVec<double,3>& X, SVec<double,dim>& V,
                                        Vec<int>& fluidId,
                                        NodalGrad<dim>& ngrad,     EdgeGrad<dim>* egrad,
					NodalGrad<dimLS> &ngradLS, EdgeGrad<dimLS>* egradLS,
                                        SVec<double,dimLS>& Phi, SVec<double,dimLS> &PhiF,
                                        LevelSetStructure* LSS, int ls_order)
{

  edges.computeFiniteVolumeTermLS(fluxFcn, recFcn, recFcnLS, elems, geoState, X, V,
                                  fluidId, ngrad, egrad,
				  ngradLS, egradLS, Phi, PhiF, LSS, ls_order);

  faces.computeFiniteVolumeTermLS(fluxFcn, bcData, geoState, V, Phi, PhiF);
    //Note: LSS is not needed because we assume that farfield nodes cannot be covered by structure.
}

//------------------------------------------------------------------------------

template<int dim>
int SubDomain::computeFiniteVolumeBar_Step1(Vec<double> &irey, FluxFcn** fluxFcn, RecFcn* recFcn,
                                            BcData<dim>& bcData, GeoState& geoState,
                                            SVec<double,3>& X, SVec<double,dim>& VBar,
                                            NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
                                            SVec<double,dim>& sigma, SVec<int,2> &tag,
                                            int failsafe, int rshift)
{

  int ierr = edges.computeFiniteVolumeTerm(locToGlobNodeMap, irey, fluxFcn, recFcn, elems, geoState,
                                           X, VBar, ngrad, egrad, sigma, tag, failsafe, rshift);

  faces.computeFiniteVolumeTerm(fluxFcn, bcData, geoState, VBar, sigma);

  return(ierr);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeFiniteVolumeBar_Step2(MacroCellSet **macroCells,
                                             SVec<double,1> &volRatio,
                                             SVec<double,dim>& sigma,
                                             SVec<double,dim>& fluxes,
                                             int scopeDepth)
{

  Connectivity &nToMN = *nodesToMCNodes[scopeDepth-1];

  for (int i=0; i<nToMN.csize(); ++i) {
    if (macroCells[scopeDepth-1]->containing(i) != -1) {
       for (int j=0; j<nToMN.num(i); ++j) {
         int idx = nToMN[i][j];
         for (int k=0; k<dim; ++k)
               fluxes[idx][k] += volRatio[idx][0] * sigma[i][k];
       }
     }
  }

}

//------------------------------------------------------------------------------
template<int dim, class Scalar, int neq>
void SubDomain::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
                                                GeoState &geoState, Vec<double> &irey,
                                                SVec<double,3> &X, Vec<double> &ctrlVol,
                                                SVec<double,dim> &V, GenMat<Scalar,neq> &A,
                                                CommPattern<double>* flag)
{
  if (!flag){
    edges.computeJacobianFiniteVolumeTerm(fluxFcn, geoState, irey, X, ctrlVol, V, A);

    faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A);
  }else{
    edges.computeJacobianFiniteVolumeTerm(fluxFcn, geoState, irey, X, ctrlVol, V, A, nodeType);

    faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A, nodeType);
  }

  for (int i=0; i<ctrlVol.size(); ++i) {
    double voli = 1.0 / ctrlVol[i];
    Scalar *Aii = A.getElem_ii(i);
    for (int k=0; k<neq*neq; ++k)
      Aii[k] *= voli;
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void SubDomain::computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                                FluxFcn **fluxFcn, BcData<dim> &bcData,
                                                GeoState &geoState, Vec<double> &irey,
                                                SVec<double,3> &X, Vec<double> &ctrlVol,
                                                SVec<double,dim> &V, GenMat<Scalar,neq> &A,
                                                CommPattern<double>* flag)
{

  if (!flag){
    edges.computeJacobianFiniteVolumeTerm(fluxFcn, geoState, irey, X, ctrlVol, V, A);

    faces.computeJacobianFiniteVolumeTerm(riemann, fluxFcn, bcData, geoState, V, A);
  }else{
    edges.computeJacobianFiniteVolumeTerm(fluxFcn, geoState, irey, X, ctrlVol, V, A, nodeType);

    fprintf(stderr,"WARNING: Exact Riemann solver at the wall has not been implemented in this case!\n");
    faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A, nodeType);
  }

  for (int i=0; i<ctrlVol.size(); ++i) {
    double voli = 1.0 / ctrlVol[i];
    Scalar *Aii = A.getElem_ii(i);
    for (int k=0; k<neq*neq; ++k)
      Aii[k] *= voli;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::recomputeRHS(VarFcn* vf, SVec<double,dim>& V, SVec<double,dim>& rhs,
                             Extrapolation<dim>* xpol, BcData<dim>& bcData,
                             GeoState& geoState, SVec<double,3> &X)
{
  inletNodes.recomputeRHS(vf, xpol, elems, V, bcData, geoState, rhs, X, locToGlobNodeMap);
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::recomputeRHS(VarFcn* vf, SVec<double,dim>& V, Vec<int> &fluidId,
                            SVec<double,dim>& rhs, Extrapolation<dim>* xpol,
                            BcData<dim>& bcData, GeoState& geoState, SVec<double,3> &X)
{
  inletNodes.recomputeRHS(vf, xpol, elems, V, fluidId, bcData, geoState, rhs, X, locToGlobNodeMap);
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::recomputeResidual(SVec<double,dim> &F, SVec<double,dim> &Finlet)
{
  Finlet = 0.0;
  inletNodes.recomputeResidual(F,Finlet);
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeRealFluidResidual(SVec<double, dim> &F, SVec<double,dim> &Freal, LevelSetStructure &lss)
{
  Freal = 0.0;
  for (int iNode=0; iNode<numNodes(); iNode++)
    if (lss.isActive(0,iNode))
        for (int j=0; j<dim; j++)  Freal[iNode][j] = F[iNode][j];
}




template<class Scalar, int dim>
void SubDomain::checkRHS(Scalar (*rhs)[dim])
{
  int node;
  for (int i = 0; i<inletNodes.size(); i++){
    node = inletNodes[i].getNodeNum();
  }
}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq, int dimLS>
void SubDomain::computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                                FluxFcn **fluxFcn, BcData<dim> &bcData,
                                                GeoState &geoState,
                                                NodalGrad<dim> &ngrad, NodalGrad<dimLS> &ngradLS,
                                                SVec<double,3> &X, Vec<double> &ctrlVol,
                                                SVec<double,dim> &V, GenMat<Scalar,neq> &A,
                                                FluidSelector &fluidSelector,
                                                Vec<int> &fluidId, CommPattern<double>* flag)
{
  if (!flag){
    edges.computeJacobianFiniteVolumeTerm(riemann, fluxFcn, geoState, ngrad, ngradLS, X, ctrlVol, V, A, fluidSelector, fluidId);
    faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A, fluidId);
  }else{
    edges.computeJacobianFiniteVolumeTerm(riemann, fluxFcn, geoState, ngrad, ngradLS, X, ctrlVol, V, A, fluidSelector, fluidId, nodeType);
    faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A, fluidId, nodeType);
  }

  for (int i=0; i<ctrlVol.size(); ++i) {
    double voli = 1.0 / ctrlVol[i];
    Scalar *Aii = A.getElem_ii(i);
    for (int k=0; k<neq*neq; ++k)
      Aii[k] *= voli;
  }
}
//-------------------------------------------------------------------------------
//d2d emebedded
template<class Scalar,int dim,int neq>
void SubDomain::computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                       FluxFcn** fluxFcn, 
                                       BcData<dim>& bcData, GeoState& geoState,
                                       SVec<double,3>& X, SVec<double,dim>& V,Vec<double>& ctrlVol,
                                       LevelSetStructure &LSS, Vec<int> &fluidId,
                                       int Nriemann,
																GenMat<Scalar,neq>& A,Vec<double>& irey,  bool externalSI)
{

	if(externalSI)
		edges.computeJacobianFiniteVolumeTerm(riemann, fluxFcn,
														  geoState, X, V, ctrlVol, LSS, 
														  fluidId, Nriemann, A);
	else
  edges.computeJacobianFiniteVolumeTerm(riemann,fluxFcn,
                                        geoState, X, V, ctrlVol, LSS, 
                                        fluidId, Nriemann, A, irey);

  faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A, fluidId, &LSS); 

	for (int i=0; i<ctrlVol.size(); ++i) 
	{
    double voli = 1.0 / ctrlVol[i];
    Scalar *Aii = A.getElem_ii(i);

		for (int k=0; k<neq*neq; ++k) Aii[k] *= voli;
  }
  
}

//-------------------------------------------------------------------------------

template<int dim, class Scalar, int neq, int dimLS>
void SubDomain::computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                       FluxFcn** fluxFcn, 
                                       BcData<dim>& bcData, GeoState& geoState,
                                       SVec<double,3>& X, SVec<double,dim>& V,Vec<double>& ctrlVol,
                                       NodalGrad<dimLS> &ngradLS,
                                       LevelSetStructure &LSS,Vec<int> &fluidId,
                                       int Nriemann,
                                       FluidSelector &fluidSelector,
                                       GenMat<Scalar,neq>& A) {

  edges.computeJacobianFiniteVolumeTerm(riemann,locToGlobNodeMap,
                                        fluxFcn, geoState, X, V, LSS, fluidId, Nriemann,
                                        fluidSelector, ngradLS, ctrlVol, A);

  faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A, fluidId, &LSS);

  for (int i=0; i<ctrlVol.size(); ++i) {
    double voli = 1.0 / ctrlVol[i];
    Scalar *Aii = A.getElem_ii(i);
    for (int k=0; k<neq*neq; ++k)
      Aii[k] *= voli;
  }
}

//-------------------------------------------------------------------------------

template<int dim, class Scalar, int dimLS>
void SubDomain::computeJacobianFiniteVolumeTermLS(RecFcn* recFcn, RecFcn* recFcnLS,
						  GeoState &geoState,SVec<double,3>& X,SVec<double,dim> &V,
						  NodalGrad<dim>& ngrad, NodalGrad<dimLS> &ngradLS,
						  EdgeGrad<dim>* egrad,Vec<double> &ctrlVol,SVec<double,dimLS>& Phi,
						  GenMat<Scalar,dimLS> &A, LevelSetStructure* LSS, CommPattern<double> * flag)
{

  if (!flag){
    edges.computeJacobianFiniteVolumeTermLS(recFcn,recFcnLS,geoState,X,V,ngrad ,ngradLS,
					    egrad, ctrlVol , Phi, A,LSS);
    faces.computeJacobianFiniteVolumeTermLS(geoState, V, A);
  }else{
    edges.computeJacobianFiniteVolumeTermLS(recFcn,recFcnLS,geoState,X,V,ngrad ,ngradLS,
					    egrad,ctrlVol , Phi, A,LSS);
    faces.computeJacobianFiniteVolumeTermLS(geoState, V, A );
  }

  for (int i=0; i<ctrlVol.size(); ++i) {
    double voli = 1.0 / ctrlVol[i];
    Scalar *Aii = A.getElem_ii(i);
    for (int k=0; k<dimLS*dimLS; ++k)
      Aii[k] *= voli;
  }
}
//-------------------------------------------------------------------------------

template<class Scalar, int neq>
void SubDomain::finishJacobianGalerkinTerm(Vec<double> &ctrlVol, GenMat<Scalar,neq> &A)  {

  for (int i=0; i<ctrlVol.size(); ++i) {
    double voli = 1.0 / ctrlVol[i];
    Scalar *Aii = A.getElem_ii(i);
    for (int k=0; k<neq*neq; ++k)
      Aii[k] *= voli;
  }
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeVolumicForceTerm(VolumicForceTerm *volForce, Vec<double> &ctrlVol,
                               SVec<double,dim> &V, SVec<double,dim> &fluxes)
{

  double r[5];
  for (int i=0; i<ctrlVol.size(); i++){
    volForce->computeVolumeTerm(ctrlVol[i], V[i], r);
    for (int k=0; k<5; k++)
      fluxes[i][k] += r[k];
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeDerivativeOfVolumicForceTerm(VolumicForceTerm *volForce, Vec<double> &ctrlVol, Vec<double> &dCtrlVol,
                               SVec<double,dim> &V, SVec<double,dim> &dV, SVec<double,dim> &dFluxes)
{

  double dr[5];
  for (int i=0; i<ctrlVol.size(); i++){
    volForce->computeDerivativeOfVolumeTerm(ctrlVol[i], dCtrlVol[i], V[i], dV[i], dr);
    for (int k=0; k<5; k++)
      dFluxes[i][k] += dr[k];
  }

}

//------------------------------------------------------------------------------
/*
@TECHREPORT{fezoui-lanteri-larrouturou-olivier-89,
  author = "Fezoui, F. and Lanteri, S. and Larrouturou, B. and Olivier, C.",
  title = "R\'esolution num\'erique des \'equations de {N}avier-{S}tokes pour un fluide
  compressible en maillage triangulaire",
  institution = "INRIA, France",
  number = 1033,
  year = 1989,
}
@UNPUBLISHED{barth-91a,
  author = "Barth, T. J.",
  title = "Numerical Aspects of Computing High {R}eynolds Number Flows
           on Unstructured Meshes",
  month = jan,
  year = 1991,
  note = "{AIAA} paper 91-0721",
}
*/
template<int dim>
void SubDomain::computeGalerkinTerm(FemEquationTerm *fet, BcData<dim> &bcData,
				    GeoState &geoState, SVec<double,3> &X,
				    SVec<double,dim> &V, SVec<double,dim> &R,
												Vec<GhostPoint<dim>*> *ghostPoints,
												LevelSetStructure *LSS, 
												bool externalSI)
{

	elems.computeGalerkinTerm(fet, geoState, X, V, R,ghostPoints, LSS, externalSI);

  faces.computeGalerkinTerm(elems, fet, bcData, geoState, X, V, R,LSS);

}


//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeDerivativeOfGalerkinTerm(FemEquationTerm *fet, BcData<dim> &bcData,
				    GeoState &geoState, SVec<double,3> &X, SVec<double,3> &dX,
				    SVec<double,dim> &V, SVec<double,dim> &dV, double dMach, SVec<double,dim> &dR)
{

  elems.computeDerivativeOfGalerkinTerm(fet, geoState, X, dX, V, dV, dMach, dR);
  //std::cout<<"\033[91mmelems.computeDerivativeOfGalerkinTerm fnished on\033[00m"<<this->globSubNum<<std::endl;//TODO delete line

  faces.computeDerivativeOfGalerkinTerm(elems, fet, bcData, geoState, X, dX, V, dV, dMach, dR);
  //std::cout<<"\033[91mfaces.computeDerivativeOfGalerkinTerm fnished on\033[00m"<<this->globSubNum<<std::endl;//TODO delete line

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeDerivativeOfGalerkinTerm(RectangularSparseMat<double,3,dim> *dViscousFluxdX,
            FemEquationTerm *fet, BcData<dim> &bcData,
				    GeoState &geoState, SVec<double,3> &X, SVec<double,3> &dX,
				    SVec<double,dim> &V, SVec<double,dim> &dV, double dMach, SVec<double,dim> &dR)
{ //YC

  dViscousFluxdX->apply(dX,dR);
//  faces.computeDerivativeOfGalerkinTerm(elems, fet, bcData, geoState, X, dX, V, dV, dMach, dR2); // for Turbulent flow
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeTransposeDerivativeOfGalerkinTerm(RectangularSparseMat<double,3,dim> *dViscousFluxdX,
				                                                 SVec<double,dim> &dR, SVec<double,3> &dX)
{ //YC

  SVec<double,3> dummy(dX);
  dViscousFluxdX->applyTranspose(dR,dummy);
  dX += dummy;

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void SubDomain::computeDerivativeOperatorsOfGalerkinTerm(FemEquationTerm *fet, BcData<dim> &bcData,
            GeoState &geoState, SVec<double,3> &X, SVec<double,dim> &V, RectangularSparseMat<double,3,dim> &dViscousFluxdX)
{

  elems.computeDerivativeOperatorsOfGalerkinTerm(fet, geoState, X, V, dViscousFluxdX);

}

//------------------------------------------------------------------------------


template<int dim>
void SubDomain::computeSmagorinskyLESTerm(SmagorinskyLESTerm *smag, SVec<double,3> &X,
					  SVec<double,dim> &V, SVec<double,dim> &R, 
				          Vec<GhostPoint<dim>*> *ghostPoints,
                                          LevelSetStructure *LSS, bool externalSI)
{
	
	elems.computeSmagorinskyLESTerm(smag, X, V, R, ghostPoints, LSS, externalSI);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeWaleLESTerm(WaleLESTerm *wale, SVec<double,3> &X,
			           SVec<double,dim> &V, SVec<double,dim> &R,
				   Vec<GhostPoint<dim>*> *ghostPoints,
                                   LevelSetStructure *LSS, bool externalSI)
{

	elems.computeWaleLESTerm(wale, X, V, R, ghostPoints, LSS, externalSI);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeTestFilterAvgs(SVec<double,dim> &VCap, SVec<double,16> &Mom_Test,
                                      SVec<double,6> &Sij_Test, Vec<double> &modS_Test,
                                      SVec<double,8> &Eng_Test, SVec<double,3> &X,
				      SVec<double,dim> &V, double gam, double R,
				      Vec<GhostPoint<dim>*> *ghostPoints,
                                      LevelSetStructure *LSS, bool externalSI)
{

	elems.computeTestFilterAvgs(VCap, Mom_Test, Sij_Test, modS_Test, Eng_Test, X, V, gam, R, ghostPoints, LSS, externalSI);

}

//------------------------------------------------------------------------------
// Computes Cs and Pt values in dynamic LES procedure
//
template<int dim>
void SubDomain::computeCsValues(SVec<double,dim> &VCap, SVec<double,16> &Mom_Test,
                                SVec<double,6> &Sij_Test, Vec<double> &modS_Test,
                                SVec<double,8> &Eng_Test, SVec<double,2> &Cs,
				Vec<int> &Ni, SVec<double,3> &X, double gam, 
                                double R, LevelSetStructure *LSS)
{

 for (int i=0; i<nodes.size(); ++i) {

   if (LSS && Mom_Test[i][0] == 0.0) {
      // Node i is a part of an inactive tet in embedded simulation.
      // Otherwise Mom_Test[i][0] cannot be zero (unless there is a bug).
      // Do nothing.
      continue;
   }

   double vc[5];
   double r_u[3];
   double r_u_u[6];
   double r_s_p[6];
   double ratdelta = pow(Ni[i],(2.0/3.0)); // should precompute this also

   double r_e;
   double r_e_plus_p;
   double r_s_dtdxj[3];
   double dtdxj[3];

   double num, denom;
   double sqrt2S2;
   double Pij[3][3], Bij[3][3], Lij[3][3];
   double Zi[3], Li[3];
   double oogam1 = 1.0/(gam - 1.0);
   double Cp = R*gam*oogam1; // specific heat at constant pressure

   // Compute Smagorinsky Coefficient at each node
   // ----------------------------------------------
   // Ref: Large Eddy Simulation of Bluff-Body flow on Unstructured Grids
   // International Journal of Numerical Methods in Fluids 2002
   // Vol : 40, pgs:1431-1460
   // Authors; Camarri, Salvetti, Koobus, Dervieux
   // ---------------------------------------------------------------------

    r_u[0] = Mom_Test[i][1]/Mom_Test[i][0];
    r_u[1] = Mom_Test[i][2]/Mom_Test[i][0];
    r_u[2] = Mom_Test[i][3]/Mom_Test[i][0];

    r_u_u[0] = Mom_Test[i][4]/Mom_Test[i][0];
    r_u_u[1] = Mom_Test[i][5]/Mom_Test[i][0];
    r_u_u[2] = Mom_Test[i][6]/Mom_Test[i][0];
    r_u_u[3] = Mom_Test[i][7]/Mom_Test[i][0];
    r_u_u[4] = Mom_Test[i][8]/Mom_Test[i][0];
    r_u_u[5] = Mom_Test[i][9]/Mom_Test[i][0];

    r_s_p[0] = Mom_Test[i][10]/Mom_Test[i][0];
    r_s_p[1] = Mom_Test[i][11]/Mom_Test[i][0];
    r_s_p[2] = Mom_Test[i][12]/Mom_Test[i][0];
    r_s_p[3] = Mom_Test[i][13]/Mom_Test[i][0];
    r_s_p[4] = Mom_Test[i][14]/Mom_Test[i][0];
    r_s_p[5] = Mom_Test[i][15]/Mom_Test[i][0];

    vc[0] = VCap[i][0]/Mom_Test[i][0];
    vc[1] = VCap[i][1]/Mom_Test[i][0];
    vc[2] = VCap[i][2]/Mom_Test[i][0];
    vc[3] = VCap[i][3]/Mom_Test[i][0];
    vc[4] = VCap[i][4]/Mom_Test[i][0];

    Pij[0][0] = Sij_Test[i][0]/Mom_Test[i][0];
    Pij[1][1] = Sij_Test[i][1]/Mom_Test[i][0];
    Pij[2][2] = Sij_Test[i][2]/Mom_Test[i][0];
    Pij[0][1] = Sij_Test[i][3]/Mom_Test[i][0];
    Pij[0][2] = Sij_Test[i][4]/Mom_Test[i][0];
    Pij[1][2] = Sij_Test[i][5]/Mom_Test[i][0];
    Pij[1][0] = Pij[0][1];
    Pij[2][0] = Pij[0][2];
    Pij[2][1] = Pij[1][2];

    sqrt2S2 = modS_Test[i]/Mom_Test[i][0];

    r_e = Eng_Test[i][0]/Mom_Test[i][0];

    r_e_plus_p = Eng_Test[i][1]/Mom_Test[i][0];

    r_s_dtdxj[0] = Eng_Test[i][2]/Mom_Test[i][0];
    r_s_dtdxj[1] = Eng_Test[i][3]/Mom_Test[i][0];
    r_s_dtdxj[2] = Eng_Test[i][4]/Mom_Test[i][0];

    dtdxj[0] = Eng_Test[i][5]/Mom_Test[i][0];
    dtdxj[1] = Eng_Test[i][6]/Mom_Test[i][0];
    dtdxj[2] = Eng_Test[i][7]/Mom_Test[i][0];

    // computing Smagorinsky constant //

    computeLij(Lij, r_u, r_u_u, vc);
    computeBij(Bij, r_s_p, sqrt2S2, Pij, ratdelta, vc);


    num = 0.0;
    denom = 0.0;

    // least squares procedure

    for (int j=0; j<3; ++j){
      for (int k=0; k<3; ++k){
         num += Bij[j][k]*Lij[j][k];
         denom += Bij[j][k]*Bij[j][k];
       }
    }

    if(denom < 1.0e-7) denom = 1.0e-7;
    Cs[i][0] = (num/denom);

    // computing turbulent Prandtl number //

    computeZi(Zi, ratdelta, sqrt2S2, dtdxj, r_s_dtdxj, vc, Cp);
    computeLi(Li, r_e, r_e_plus_p, r_u, vc);

    num = 0.0;
    denom = 0.0;

    // least squares procedure

    for (int j=0; j<3; ++j){
       num += Li[j]*Zi[j];
       denom += Zi[j]*Zi[j];
    }

    if(denom < 1.0e-7) denom = 1.0e-7;
    Cs[i][1] = (num/denom);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeDynamicLESTerm(DynamicLESTerm *dles, SVec<double,2> &Cs,
                                      SVec<double,3> &X, SVec<double,dim> &V, 
                                      SVec<double,dim> &R,
				      Vec<GhostPoint<dim>*> *ghostPoints,
                                      LevelSetStructure *LSS, bool externalSI)
{

	elems.computeDynamicLESTerm(dles, Cs, X, V, R, ghostPoints, LSS, externalSI);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeVMSLES_Step1(VMSLESTerm *vmst,
                                    SVec<double,dim> &VBar,
                                    SVec<double,3> &X,
                                    SVec<double,dim> &V,
                                    SVec<double,dim> &Sigma)
{
  elems.computeVMSLESTerm(vmst, VBar,X, V, Sigma);
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeVMSLES_Step2(SVec<double,1> &volRatio,
                                    MacroCellSet *macroCells,
                                    SVec<double,dim> &Sigma,
                                    SVec<double,dim> &R,
                                    int scopeDepth)
{

  Connectivity &nToMN = *nodesToMCNodes[scopeDepth-1];

  for (int i=0; i<nToMN.csize(); ++i) {
    if (macroCells->containing(i) != -1) {
      // excluding nodes on SD boundary
      // that are assigned to a MC in another SD
      for (int j=0; j<nToMN.num(i); ++j) {
	int idx = nToMN[i][j];
	if (i == idx) {
	  for (int k=0; k<dim; ++k)
	    R[idx][k] += (1.0 - volRatio[idx][0]) * Sigma[i][k];
	}
	else {
	  for (int k=0; k<dim; ++k)
	    R[idx][k] += -1.0 * volRatio[idx][0] * Sigma[i][k];
	}
      }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeGalerkinBar_Step1(FemEquationTerm *fet,
                                         BcData<dim> &bcData,
                                         GeoState &geoState,
                                         SVec<double,3> &X,
                                         SVec<double,dim> &VBar,
                                         SVec<double,dim> &Sigma)
{

  elems.computeGalerkinTerm(fet, geoState, X, VBar, Sigma);

  faces.computeGalerkinTerm(elems, fet, bcData, geoState, X, VBar, Sigma);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeGalerkinBar_Step2(MacroCellSet **macroCells,
                                         SVec<double,1> &volRatio,
                                         SVec<double,dim> &Sigma,
                                         SVec<double,dim> &RBar,
                                         int scopeDepth)
{

  Connectivity &nToMN = *nodesToMCNodes[scopeDepth-1];

  for (int i=0; i<nToMN.csize(); ++i) {
    if (macroCells[scopeDepth-1]->containing(i) != -1) {
       for (int j=0; j<nToMN.num(i); ++j) {
         int idx = nToMN[i][j];
         for (int k=0; k<dim; ++k)
               RBar[idx][k] += volRatio[idx][0] * Sigma[i][k];
       }
     }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeMBarAndM_Step1(DynamicVMSTerm *dvmst,
                                      SVec<double,dim> **VBar,
                                      SVec<double,1> **volRatio,
                                      SVec<double,3> &X,
                                      SVec<double,dim> &V,
                                      SVec<double,dim> &SigmaBar,
                                      SVec<double,dim> &Sigma)
{

  elems.computeMBarAndM(dvmst, VBar, volRatio, X, V, SigmaBar, Sigma);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeMBarAndM_Step2(MacroCellSet **macroCells,
                                      SVec<double,1> **volRatio,
                                      SVec<double,dim> &MBar,
                                      SVec<double,dim> &M,
                                      SVec<double,dim> &SigmaBar,
                                      SVec<double,dim> &Sigma,
                                      int scopeDepth1,
                                      int scopeDepth2)
{

  Connectivity &nToMN1 = *nodesToMCNodes[scopeDepth1-1];
  Connectivity &nToMN2 = *nodesToMCNodes[scopeDepth2-1];

  for (int i=0; i<nToMN1.csize(); ++i) {
    if (macroCells[scopeDepth1-1]->containing(i) != -1) {
       for (int j=0; j<nToMN1.num(i); ++j) {
         for (int k=0; k<dim; ++k)
            MBar[ nToMN1[i][j] ][k] += (*volRatio[0])[ nToMN1[i][j] ][0] * SigmaBar[i][k];
         if (i == nToMN1[i][j]) {
            for (int k=0; k<dim; ++k)
               M[ nToMN1[i][j] ][k] += (1.0 - (*volRatio[0])[ nToMN1[i][j] ][0]) * Sigma[i][k];
         }
         else {
            for (int k=0; k<dim; ++k)
               M[ nToMN1[i][j] ][k] += -1.0 * (*volRatio[0])[ nToMN1[i][j] ][0] * Sigma[i][k];
         }
       }
    }
  }

  for (int i=0; i<nToMN2.csize(); ++i) {
    if (macroCells[scopeDepth2-1]->containing(i) != -1) {
       for (int j=0; j<nToMN2.num(i); ++j) {
         for (int k=0; k<dim; ++k)
           MBar[ nToMN2[i][j] ][k] += -1.0 * (*volRatio[1])[ nToMN2[i][j] ][0] * SigmaBar[i][k];
       }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeDynamicVMSTerm_Step1(DynamicVMSTerm *dvmst,
                                            SVec<double,dim> **VBar,
                                            SVec<double,3> &X,
                                            SVec<double,dim> &V,
                                            SVec<double,dim> &Sigma,
                                            Vec<double> &CsDelSq,
                                            Vec<double> &PrT,
                                            Vec<double> *Cs,
                                            Vec<double> &Delta)
{

  elems.computeDynamicVMSTerm(dvmst, VBar, X, V, Sigma, CsDelSq, PrT, Cs, Delta);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeDynamicVMSTerm_Step2(MacroCellSet **macroCells,
                                            SVec<double,1> **volRatio,
                                            SVec<double,dim> &Sigma,
                                            SVec<double,dim> &R,
                                            int scopeDepth)
{

  Connectivity &nToMN = *nodesToMCNodes[scopeDepth-1];

  for (int i=0; i<nToMN.csize(); ++i) {
    if (macroCells[scopeDepth-1]->containing(i) != -1) {
       for (int j=0; j<nToMN.num(i); ++j) {
         if (i == nToMN[i][j]) {
            for (int k=0; k<dim; ++k)
               R[ nToMN[i][j] ][k] += (1.0 - (*volRatio[0])[ nToMN[i][j] ][0]) * Sigma[i][k];
         }
         else {
            for (int k=0; k<dim; ++k)
               R[ nToMN[i][j] ][k] += -1.0 * (*volRatio[0])[ nToMN[i][j] ][0] * Sigma[i][k];
         }
       }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computedWBar_dt(MacroCellSet **macroCells,
                                SVec<double,1> **volRatio,
                                SVec<double,dim> &Sigma,
                                SVec<double,dim> &dWBardt,
                                int scopeDepth)
{

  Connectivity &nToMN = *nodesToMCNodes[scopeDepth-1];

  for (int i=0; i<nToMN.csize(); ++i) {
    if (macroCells[scopeDepth-1]->containing(i) != -1) {
       for (int j=0; j<nToMN.num(i); ++j) {
         for (int k=0; k<dim; ++k)
               dWBardt[ nToMN[i][j] ][k] += (*volRatio[0])[ nToMN[i][j] ][0] * Sigma[i][k];
       }
     }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeCsDeltaSq(SVec<double,dim> &R,
                                 SVec<double,dim> &RBar,
                                 SVec<double,dim> &M,
                                 SVec<double,dim> &MBar,
                                 SVec<double,dim> &dWdt,
                                 SVec<double,dim> &dWBardt,
                                 Vec<double> &CsDeltaSq,
                                 Vec<double> &PrT,
                                 int method)
{

 for (int i=0; i<nodes.size(); ++i) {

    switch (method) {

    case (0):     // Variational Germano Identity
      {
      double num = 0.0;
      double denom = 0.0;
      for(int j=1; j<4; ++j) {
        num += (M[i][j] - MBar[i][j])*((RBar[i][j] + dWBardt[i][j]) - (R[i][j] + dWdt[i][j]));
        denom += (M[i][j] - MBar[i][j])*(M[i][j] - MBar[i][j]);
      }
      if (fabs(denom) < 0.000001) CsDeltaSq[i] = 0.0;
      else CsDeltaSq[i] = num / denom;

      num = CsDeltaSq[i] * (M[i][4] - MBar[i][4]);
      denom = (RBar[i][4] + dWBardt[i][4]) - (R[i][4] + dWdt[i][4]);
      if (fabs(denom) < 0.000001) PrT[i] = 0.9;
      else PrT[i] = num / denom;
      }
      break;

    case(1):     // Full Least Squares
      {
      double num = 0.0;
      double denom = 0.0;
      for(int j=1; j<4; ++j) {
        num += -(MBar[i][j] * (RBar[i][j] + dWBardt[i][j]))
               -(M[i][j] * (R[i][j] + dWdt[i][j]));
        denom += (MBar[i][j] * MBar[i][j]) + (M[i][j] * M[i][j]);
      }
      if (fabs(denom) < 0.000001) CsDeltaSq[i] = 0.0;
      else CsDeltaSq[i] = num / denom;

      num = -CsDeltaSq[i] * ((R[i][4] + dWdt[i][4])*M[i][4]
                           + (RBar[i][4] + dWBardt[i][4])*MBar[i][4]);
      denom = pow((RBar[i][4] + dWBardt[i][4]),2.0) + pow((R[i][4] + dWdt[i][4]),2.0);
      if (fabs(denom) < 0.000001) PrT[i] = 0.9;
      else PrT[i] = num / denom;
      }
      break;

    case(2):   // Special Clipping Procedure
      {
      double num = 0.0;
      double denom = 0.0;
      for(int j=1; j<4; ++j) {
        num += (M[i][j] - MBar[i][j])*((RBar[i][j] + dWBardt[i][j]) - (R[i][j] + dWdt[i][j]));
        denom += (M[i][j] - MBar[i][j])*(M[i][j] - MBar[i][j]);
      }
      if (fabs(denom) < 0.000001) CsDeltaSq[i] = 0.0;
      else CsDeltaSq[i] = num / denom;

      if (CsDeltaSq[i] < 0.0) {
        num = 0.0;
        denom = 0.0;
        for(int j=1; j<4; ++j) {
          num += pow((RBar[i][j] + dWBardt[i][j]) - (R[i][j] + dWdt[i][j]) , 2.0);
          denom += pow((M[i][j]) - MBar[i][j], 2.0);
        }
        if (fabs(denom) < 0.000001) CsDeltaSq[i] = 0.0;
        else CsDeltaSq[i] = sqrt(num / denom);
      }

      num = CsDeltaSq[i] * (M[i][4] - MBar[i][4]);
      denom = (RBar[i][4] + dWBardt[i][4]) - (R[i][4] + dWdt[i][4]);
      if (denom < 0.000001) PrT[i] = 0.9;
      else PrT[i] = num / denom;
      }
      break;

    default:
      fprintf(stderr,"Error :: Method to Solve the Residual Equation is Not Correct...Aborting !!\n");
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeMutOMuSmag(SmagorinskyLESTerm *smag, SVec<double,3> &X,
                                  SVec<double,dim> &V, Vec<double> &mutOmu)
{

  for (int tetNum=0; tetNum < elems.size(); ++tetNum) {
    double dp1dxj[4][3];
    double vol = elems[tetNum].computeGradientP1Function(X, dp1dxj);
    double *v[4] = {V[elems[tetNum][0]], V[elems[tetNum][1]],
                    V[elems[tetNum][2]], V[elems[tetNum][3]]};
    double mut = smag->computeMutOMu(vol, dp1dxj, v, X, elems[tetNum]);
    for (int i=0; i<4; ++i)
      mutOmu[elems[tetNum][i]] += mut * vol;
  }

}

//--------------------------------------------------------------------------

template<int dim>
void SubDomain::computeMutOMuVMS(VMSLESTerm *vmst, SVec<double,dim> &VBar, SVec<double,3> &X,
                                 SVec<double,dim> &V, Vec<double> &mutOmu)
{

   for (int tetNum=0; tetNum < elems.size(); ++tetNum) {
     double dp1dxj[4][3];
     double vol = elems[tetNum].computeGradientP1Function(X, dp1dxj);
     double *v[4] = {V[elems[tetNum][0]], V[elems[tetNum][1]],
                    V[elems[tetNum][2]], V[elems[tetNum][3]]};
     double *vbar[4] = {VBar[elems[tetNum][0]], VBar[elems[tetNum][1]],
                        VBar[elems[tetNum][2]], VBar[elems[tetNum][3]]};

     double mut = vmst->computeMutOMu(vol, dp1dxj, vbar, v, X, elems[tetNum]);
    for (int i=0; i<4; ++i)
      mutOmu[elems[tetNum][i]] += mut * vol;

  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeMutOMuDynamicVMS(DynamicVMSTerm *dvmst, SVec<double,dim> &VBar, SVec<double,3> &X,
                                        SVec<double,dim> &V, Vec<double> &Cs, Vec<double> &mutOmu)
{

   for (int tetNum=0; tetNum < elems.size(); ++tetNum) {
     double dp1dxj[4][3];
     double vol = elems[tetNum].computeGradientP1Function(X, dp1dxj);
     double *v[4] = {V[elems[tetNum][0]], V[elems[tetNum][1]],
                    V[elems[tetNum][2]], V[elems[tetNum][3]]};
     double *vbar[4] = {VBar[elems[tetNum][0]], VBar[elems[tetNum][1]],
                        VBar[elems[tetNum][2]], VBar[elems[tetNum][3]]};
     double cs[4] = {Cs[elems[tetNum][0]], Cs[elems[tetNum][1]],
                     Cs[elems[tetNum][2]], Cs[elems[tetNum][3]]};

     double mut = dvmst->computeMutOMu(vol, dp1dxj, vbar, v, cs, X, elems[tetNum]);
     for (int i=0; i<4; ++i)
       mutOmu[elems[tetNum][i]] += mut * vol;

  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeMutOMuWale(WaleLESTerm *wale, SVec<double,3> &X,
                                  SVec<double,dim> &V, Vec<double> &mutOmu)
{

  for (int tetNum=0; tetNum < elems.size(); ++tetNum) {
    double dp1dxj[4][3];
    double vol = elems[tetNum].computeGradientP1Function(X, dp1dxj);
    double *v[4] = {V[elems[tetNum][0]], V[elems[tetNum][1]],
                    V[elems[tetNum][2]], V[elems[tetNum][3]]};
    double mut = wale->computeMutOMu(vol, dp1dxj, v, X, elems[tetNum]);
    for (int i=0; i<4; ++i)
      mutOmu[elems[tetNum][i]] += mut * vol;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeMutOMuDynamicLES(DynamicLESTerm *dles, SVec<double,2> &Cs,
                                  SVec<double,3> &X, SVec<double,dim> &V,
				  Vec<double> &mutOmu)

{

  for (int tetNum=0; tetNum < elems.size(); ++tetNum) {
    double dp1dxj[4][3];
    double vol = elems[tetNum].computeGradientP1Function(X, dp1dxj);
    double *v[4] = {V[elems[tetNum][0]], V[elems[tetNum][1]],
                    V[elems[tetNum][2]], V[elems[tetNum][3]]};
    double cs[4] = {Cs[elems[tetNum][0]][0], Cs[elems[tetNum][1]][0],
                    Cs[elems[tetNum][2]][0], Cs[elems[tetNum][3]][0]};
    double pt[4] = {Cs[elems[tetNum][0]][1], Cs[elems[tetNum][1]][1],
                    Cs[elems[tetNum][2]][1], Cs[elems[tetNum][3]][1]};

    double mut = dles->computeMutOMu(vol, dp1dxj, v, cs, pt, X, elems[tetNum]);
    for (int i=0; i<4; ++i)
      mutOmu[elems[tetNum][i]] += mut * vol;
  }

}

//--------------------------------------------------------------------------
template<int dim, class Scalar, int neq>
void SubDomain::computeJacobianGalerkinTerm(FemEquationTerm *fet, BcData<dim> &bcData,
					    GeoState &geoState, SVec<double,3> &X,
					    Vec<double> &ctrlVol, SVec<double,dim> &V,
					    GenMat<Scalar,neq> &A,
														  Vec<GhostPoint<dim>*>* ghostPoints, 
														  LevelSetStructure *LSS, bool externalSI)
{

	elems.computeJacobianGalerkinTerm(fet, geoState, X, ctrlVol, V, A, ghostPoints, LSS, externalSI);

  faces.computeJacobianGalerkinTerm(elems, fet, bcData, geoState, X, ctrlVol, V, A);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeBCsJacobianWallValues(FemEquationTerm *fet, BcData<dim> &bcData,
					    GeoState &geoState, SVec<double,3> &X,
					    SVec<double,dim> &V)
{

  faces.computeBCsJacobianWallValues(elems, fet, bcData, geoState, X, V);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void SubDomain::computeJacobianVolumicForceTerm(VolumicForceTerm *volForce,
                                           Vec<double> &ctrlVol, SVec<double,dim> &V,
                                           GenMat<Scalar,neq> &A)
{
  /* computation of the jacobian part due to gravity
   * There is only a block diagonal term to compute for each node.
   */
  if (neq>=5){
    Scalar *Aii;
    double jac[neq*neq];
    for (int i=0; i<ctrlVol.size(); i++){
      Aii = A.getElem_ii(i);
      volForce->computeJacobianVolumeTerm(neq, ctrlVol[i], V[i], jac);
      for (int m=0; m<neq*neq; m++)
        Aii[m] += jac[m];

    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::getExtrapolationValue(Extrapolation<dim>* xpol,SVec<double,dim> &V, SVec<double,dim> &Ubc,
				      VarFcn *vf, BcData<dim>& bcData, GeoState& geoState, SVec<double,3>& X)
{
        inletNodes.getExtrapolationValue(xpol, V, Ubc, vf, bcData, geoState, elems, locToGlobNodeMap, X);
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::applyExtrapolationToSolutionVector(Extrapolation<dim>* xpol,SVec<double,dim> &U,
						   SVec<double,dim> &Ubc)
{
        inletNodes.applyExtrapolationToSolutionVector(xpol, U, Ubc, locToGlobNodeMap);
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::applyBCsToSolutionVector(BcFcn *bcFcn, BcData<dim> &bcData,
                                         SVec<double,dim> &U, LevelSetStructure *LSS)
{
  SVec<double,dim> &Vwall = bcData.getNodeStateVector();

  // In the case of an Embedded Simulation, we want inactive nodes to stay at initial state. Thus no BC should be apply 
  // to them. 
  bool isActive = true; 

  for (int i=0; i<nodes.size(); ++i) {
    if(LSS) isActive = LSS->isActive(0.0,i);
    if (nodeType[i] != BC_INTERNAL && isActive)
      bcFcn->applyToSolutionVector(nodeType[i], Vwall[i], U[i]);
  }
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::applyBCsToTurbSolutionVector(BcFcn *bcFcn, BcData<dim> &bcData,
                                         SVec<double,dim> &U, LevelSetStructure *LSS)
{
  SVec<double,dim> &Vwall = bcData.getNodeStateVector();

  if (offWallNode && dim>5) {
    for (int i=0; i<nodes.size(); ++i) {
      if (offWallNode[i]) 
        bcFcn->applyToTurbSolutionVector(nodeType[i], Vwall[i], U[i]);
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::applyBCsToResidual(BcFcn *bcFcn, BcData<dim> &bcData,
											  SVec<double,dim> &U, SVec<double,dim> &F, 
											  LevelSetStructure *LSS)
{

  SVec<double,dim> &Vwall = bcData.getNodeStateVector();

	// In the case of an Embedded Simulation, we want inactive nodes to stay at initial state. 
	// Thus no BC should be apply to them. 
  bool isActive = true; 

	if(sampleMesh) 
	{
		int i;
		for(int iNode=0; iNode<numSampledNodes; ++iNode) 
		{
			i = locSampleNodes[iNode];
			if (nodeType[i] != BC_INTERNAL && isActive)
				bcFcn->applyToResidualTerm(nodeType[i], Vwall[i], U[i], F[i]);
		}
	}
	else 
	{
		for(int i=0; i<nodes.size(); ++i) 
		{
			if (nodeType[i] != BC_INTERNAL)
				bcFcn->applyToResidualTerm(nodeType[i], Vwall[i], U[i], F[i]);
		}

		if(offWallNode && dim>5) 
		{
			for(int i=0; i<nodes.size(); ++i) 
			{
		    if (offWallNode[i]) 
			bcFcn->applyToTurbResidualTerm(nodeType[i], Vwall[i], U[i], F[i]);
		  }
		}
	}
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::applyBCsToDerivativeOfResidual(BcFcn *bcFcn, BcData<dim> &bcData,
				   SVec<double,dim> &U, SVec<double,dim> &dU, SVec<double,dim> &dF)
{

  SVec<double,dim> &Vwall = bcData.getNodeStateVector();
  SVec<double,dim> &dVwall = bcData.getdNodeStateVector();

  for (int i=0; i<nodes.size(); ++i)
    if (nodeType[i] != BC_INTERNAL)
      bcFcn->applyToDerivativeOfResidualTerm(nodeType[i], Vwall[i], dVwall[i], U[i], dU[i], dF[i]);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void SubDomain::applyBCsToH2Jacobian(BcFcn *bcFcn, BcData<dim> &bcs,
				   SVec<double,dim> &U, GenMat<Scalar,neq> &A)
{
  SVec<double,dim> &Vwall = bcs.getNodeStateVector();

  int (*edgePtr)[2] = edges.getPtr();

  int k;
  for (int l=0; l<edges.size(); ++l) {

    if (bcMap.find(l) != bcMap.end())  {
      int i = edgePtr[l][0];
      int j = edgePtr[l][1];

      if (nodeType[i] != BC_INTERNAL)  {
        Scalar *Aij = A.getBcElem_ij(bcMap[l]);
        Scalar *Aij_orig = A.getElem_ij(l);  // the Aij is the off-diagonal term for this equation
        if (Aij && Aij_orig)  {
          for (k = 0; k < neq*neq; k++)
            Aij[k] = Aij_orig[k];
        }

        Scalar *Aji = A.getBcElem_ji(bcMap[l]);
        Scalar *Aji_orig = A.getElem_ji(l);  // Aij is the diag term for the ith eq.
        if (Aji && Aji_orig)  {
          for (k = 0; k < neq*neq; k++)
            Aji[k] = Aji_orig[k];
        }

        bcFcn->applyToOffDiagonalTerm(nodeType[i], Aij);
        bcFcn->applyToOffDiagonalTerm(nodeType[i], Aji);

      }
      if (nodeType[j] != BC_INTERNAL) {
        Scalar *Aij = A.getBcElem_ij(bcMap[l]+numBcNodes[l]);
        Scalar *Aij_orig = A.getElem_ij(l);  // Aij is the diag term for the jth eq.
        if (Aij && Aij_orig)  {
          for (k = 0; k < neq*neq; k++)
            Aij[k] = Aij_orig[k];
        }

        Scalar *Aji = A.getBcElem_ji(bcMap[l]+numBcNodes[l]);
        Scalar *Aji_orig = A.getElem_ji(l);  // Aji is the off-diag term for the jth eq.
        if (Aji && Aji_orig)  {
          for (k = 0; k < neq*neq; k++)
            Aji[k] = Aji_orig[k];
        }

        bcFcn->applyToOffDiagonalTerm(nodeType[j], Aij);
        bcFcn->applyToOffDiagonalTerm(nodeType[j], Aji);

      }
    }
  }
  for (int i=0; i<nodes.size(); ++i) {
    if (nodeType[i] != BC_INTERNAL) {
      Scalar *Aii = A.getElem_ii(i);
      if (Aii)
        bcFcn->applyToOffDiagonalTerm(nodeType[i], Aii);
    }
  }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar>
void SubDomain::applyBCsToH2Jacobian(BcFcn *bcFcn, BcData<dim> &bcs,
				   SVec<double,dim> &U, GenMat<Scalar,dim> &A)
{

  SVec<double,dim> &Vwall = bcs.getNodeStateVector();

  int (*edgePtr)[2] = edges.getPtr();

  int k;
  for (int l=0; l<edges.size(); ++l) {

    if (bcMap.find(l) != bcMap.end())  {
      int i = edgePtr[l][0];
      int j = edgePtr[l][1];

      if (nodeType[i] != BC_INTERNAL)  {
        Scalar *Aij = A.getBcElem_ij(bcMap[l]);
        Scalar *Aij_orig = A.getElem_ij(l);  // the Aij is the off-diagonal term for this equation
        if (Aij && Aij_orig)  {
          for (k = 0; k < dim*dim; k++)
            Aij[k] = Aij_orig[k];
        }

        Scalar *Aji = A.getBcElem_ji(bcMap[l]);
        Scalar *Aji_orig = A.getElem_ji(l);  // Aij is the diag term for the ith eq.
        if (Aji && Aji_orig)  {
          for (k = 0; k < dim*dim; k++)
            Aji[k] = Aji_orig[k];
        }

        bcFcn->applyToOffDiagonalTerm(nodeType[i], Aij);
        bcFcn->applyToOffDiagonalTerm(nodeType[i], Aji);

      }
      if (nodeType[j] != BC_INTERNAL) {
        Scalar *Aij = A.getBcElem_ij(bcMap[l]+numBcNodes[l]);
        Scalar *Aij_orig = A.getElem_ij(l);  // Aij is the diag term for the jth eq.
        if (Aij && Aij_orig)  {
          for (k = 0; k < dim*dim; k++)
            Aij[k] = Aij_orig[k];
        }

        Scalar *Aji = A.getBcElem_ji(bcMap[l]+numBcNodes[l]);
        Scalar *Aji_orig = A.getElem_ji(l);  // Aji is the off-diag term for the jth eq.
        if (Aji && Aji_orig)  {
          for (k = 0; k < dim*dim; k++)
            Aji[k] = Aji_orig[k];
        }

        bcFcn->applyToOffDiagonalTerm(nodeType[j], Aij);
        bcFcn->applyToOffDiagonalTerm(nodeType[j], Aji);

      }
    }
  }
  for (int i=0; i<nodes.size(); ++i) {
    if (nodeType[i] != BC_INTERNAL) {
      Scalar *Aii = A.getElem_ii(i);
      if (Aii)
        bcFcn->applyToOffDiagonalTerm(nodeType[i], Aii);
    }
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar>
void SubDomain::applyBCsToProduct(BcFcn *bcFcn, BcData<dim> &bcs, SVec<double,dim> &U, SVec<Scalar,dim> &Prod)
{

  for (int i=0; i<nodes.size(); ++i)
    if (nodeType[i] != BC_INTERNAL)
        bcFcn->applyToProductTerm(nodeType[i], Prod[i]);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void SubDomain::applyBCsToJacobian(BcFcn *bcFcn, BcData<dim> &bcs,
                                   SVec<double,dim> &U, GenMat<Scalar,neq> &A, LevelSetStructure *LSS)
{
  SVec<double,dim> &Vwall = bcs.getNodeStateVector();

  int (*edgePtr)[2] = edges.getPtr();
  bool *edgeFlag = edges.getMasterFlag();

  for (int l=0; l<edges.size(); ++l) {
    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    if (nodeType[i] != BC_INTERNAL)  {
        Scalar *Aij = A.getElem_ij(l);
        if (Aij)
          bcFcn->applyToOffDiagonalTerm(nodeType[i], Aij);
    }

    if (nodeType[j] != BC_INTERNAL) {
      Scalar *Aji = A.getElem_ji(l);
      if (Aji)
        bcFcn->applyToOffDiagonalTerm(nodeType[j], Aji);
    }
  }

  for (int i=0; i<nodes.size(); ++i) {
    if (nodeType[i] != BC_INTERNAL) {
      Scalar *Aii = A.getElem_ii(i);
      if (Aii)
        bcFcn->applyToDiagonalTerm(nodeType[i], Vwall[i], U[i], Aii);
    }
  }

  if ( LSS && offWallNode && neq!=5 ) {
    for (int l=0; l<edges.size(); ++l) {
      int i = edgePtr[l][0];
      int j = edgePtr[l][1];

      if (offWallNode[i]) {
	Scalar *Aij = 0;
        if (LSS->edgeIntersectsStructure(0.0,l))
          Aij = A.getRealNodeElem_ij(i,j);
	else
          Aij = A.getElem_ij(l);

        if (Aij)
          bcFcn->applyToTurbOffDiagonalTerm(nodeType[i], Aij);


        Scalar *Aii = A.getElem_ii(i);
        if (Aii)
          bcFcn->applyToTurbDiagonalTerm(nodeType[i], Vwall[i], U[i], Aii);

      }

      if (offWallNode[j]) {
	Scalar *Aji = 0;
        if (LSS->edgeIntersectsStructure(0.0,l))
          Aji = A.getRealNodeElem_ij(j,i);
	else
          Aji = A.getElem_ji(l);

        if (Aji)
          bcFcn->applyToTurbOffDiagonalTerm(nodeType[j], Aji);

        Scalar *Ajj = A.getElem_ii(j);
        if (Ajj)
          bcFcn->applyToTurbDiagonalTerm(nodeType[j], Vwall[j], U[j], Ajj);

      }
    }
  }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim, class Scalar, int neq>
void SubDomain::applyBCsToJacobianWallValues(BcFcn *bcFcn, BcData<dim> &bcs,
                                   SVec<double,dim> &U, GenMat<Scalar,neq> &A)
{

  SVec<double,dim> &Vwall = bcs.getNodeStateVector();
  SVec<double,dim> &dVwall = bcs.getdNodeStateVectorSA();

  for (int i=0; i<nodes.size(); ++i) {
    if (nodeType[i] != BC_INTERNAL) {
      Scalar *Aii = A.getElem_ii(i);
      if (Aii)
        bcFcn->applyToDiagonalTerm(nodeType[i], Vwall[i], dVwall[i], U[i], Aii);
    }
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
SparseMat<Scalar,dim> *SubDomain::createMaskJacobian(int *ndType, MemoryPool *mp)
{

  nodeToNodeMaskJacobian = createElemBasedConnectivity();

  int *ia = (*nodeToNodeMaskJacobian).ptr();
  int *ja = (*nodeToNodeMaskJacobian)[0];
  int n = nodes.size();
  int nnz = ia[n];

  Scalar (*a)[dim*dim] = 0;

  if (mp)
    a = reinterpret_cast<Scalar (*)[dim*dim]>(mp->request(nnz * dim*dim * sizeof(Scalar)));

  SparseMat<Scalar,dim> *A = new SparseMat<Scalar,dim>(n, nnz, ia, ja, a, 0, ndType);

  return A;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
MvpMat<Scalar,dim> *SubDomain::createMaskMatVecProd(bool nsFlag)
{

// Original
//  int numOffDiagEntries = 0;

// Included (MB*)
  if ((nsFlag) && (!numOffDiagEntries))  {
// Original
//  if (nsFlag)  {
    numBcNodes = new int[edges.size()];
    int (*edgePtr)[2] = edges.getPtr();

    for (int iEdge = 0; iEdge < edges.size(); iEdge++)  {
      int i = edgePtr[iEdge][0];
      int j = edgePtr[iEdge][1];

      numBcNodes[iEdge] = 0;

      if (nodeType[i] != BC_INTERNAL)   {
        bcMap[iEdge] = numOffDiagEntries;
        numOffDiagEntries++;
      }

     if (nodeType[j] != BC_INTERNAL)  {
       if (bcMap.find(iEdge) == bcMap.end())
         bcMap[iEdge] = numOffDiagEntries;
       else
         numBcNodes[iEdge] = 1;
       numOffDiagEntries++;
     }
    }
  }

// Included (MB*)
  if (nsFlag) {
    MvpMat<Scalar,dim> *A = new MvpMat<Scalar,dim>(nodes.size(), edges.size(), numOffDiagEntries);

    return A;
  }
  else {
    MvpMat<Scalar,dim> *A = new MvpMat<Scalar,dim>(nodes.size(), edges.size(), 0);

    return A;
  }
// Original
//  MvpMat<Scalar,dim> *A = new MvpMat<Scalar,dim>(nodes.size(), edges.size(), numOffDiagEntries);
//
//  return A;

}

//------------------------------------------------------------------------------

template<int dim, int dim2>
RectangularSparseMat<double,dim,dim2> *SubDomain::create_NodeBaseddRdXoperators()
{
  Connectivity *nodeToNode = createEdgeBasedConnectivity();
  int numNodes = nodes.size();
  int *ia = (*nodeToNode).ptr();
  int *ja = (*nodeToNode)[0];
  double (*a)[dim*dim2] = 0;

  RectangularSparseMat<double, dim, dim2> *A = new RectangularSparseMat<double, dim, dim2>(numNodes, ia[numNodes], ia, ja, a, 0, 0);
  return A;
}

//------------------------------------------------------------------------------

template<int dim, int dim2>
RectangularSparseMat<double,dim,dim2> *SubDomain::create_NodeToConstantBaseddRdXoperators()
{
  Connectivity *nodeToConstant = createNodeToConstantConnectivity();
  int numNodes = nodes.size();
  int *ia = (*nodeToConstant).ptr();
  int *ja = (*nodeToConstant)[0];
  double (*a)[dim*dim2] = 0;

  RectangularSparseMat<double, dim, dim2> *A = new RectangularSparseMat<double, dim, dim2>(numNodes, ia[numNodes], ia, ja, a, 0, 0);
  return A;
}

//------------------------------------------------------------------------------

template<int dim, int dim2>
RectangularSparseMat<double,dim,dim2> *SubDomain::create_ConstantToNodeBaseddRdXoperators()
{
  Connectivity *constantToNode = createConstantToNodeConnectivity();
  int numNodes = nodes.size();
  int *ia = (*constantToNode).ptr();
  int *ja = (*constantToNode)[0];
  double (*a)[dim*dim2] = 0;

  RectangularSparseMat<double, dim, dim2> *A = new RectangularSparseMat<double, dim, dim2>(1, ia[1], ia, ja, a, 0, 0);
  return A;
}

//------------------------------------------------------------------------------

template<int dim, int dim2>
RectangularSparseMat<double,dim,dim2> *SubDomain::create_ConstantToConstantBaseddRdXoperators()
{
  Connectivity *constantToConstant = createConstantToConstantConnectivity();
  int numNodes = nodes.size();
  int *ia = (*constantToConstant).ptr();
  int *ja = (*constantToConstant)[0];
  double (*a)[dim*dim2] = 0;

  RectangularSparseMat<double, dim, dim2> *A = new RectangularSparseMat<double, dim, dim2>(1, ia[1], ia, ja, a, 0, 0);
  return A;
}

//------------------------------------------------------------------------------


template<int dim, int dim2>
RectangularSparseMat<double,dim,dim2> *SubDomain::create_EdgeBaseddRdXoperators()
{
  Connectivity *edgeToNode = createElementBasedEdgeToNodeConnectivity();
  int numEdges = edges.size();
  int *ia = (*edgeToNode).ptr();
  int *ja = (*edgeToNode)[0];
  double (*a)[dim*dim2] = 0;

  RectangularSparseMat<double, dim, dim2> *A = new RectangularSparseMat<double, dim, dim2>(numEdges, ia[numEdges], ia, ja, a, 0, 0);
  return A;
}

//------------------------------------------------------------------------------

template<int dim, int dim2>
RectangularSparseMat<double,dim,dim2> *SubDomain::create_NodeToEdgeBaseddRdXoperators()
{
  Connectivity *nodeToEdge = createElementBasedNodeToEdgeConnectivity();
  int numNodes = nodes.size();
  int *ia = (*nodeToEdge).ptr();
  int *ja = (*nodeToEdge)[0];
  double (*a)[dim*dim2] = 0;

  RectangularSparseMat<double, dim, dim2> *A = new RectangularSparseMat<double, dim, dim2>(numNodes, ia[numNodes], ia, ja, a, 0, 0);
  return A;
}

//------------------------------------------------------------------------------

template<int dim, int dim2>
RectangularSparseMat<double,dim,dim2> *SubDomain::create_NodeToFaceBaseddRdXoperators()
{
  Connectivity *nodeToFace = createNodeToFaceConnectivity();
  int numNodes = nodes.size();
  int *ia = (*nodeToFace).ptr();
  int *ja = (*nodeToFace)[0];
  double (*a)[dim*dim2] = 0;
  
  RectangularSparseMat<double, dim, dim2> *A = new RectangularSparseMat<double, dim, dim2>(numNodes, ia[numNodes], ia, ja, a, 0, 0);
  return A;
} 

//------------------------------------------------------------------------------

template<int dim, int dim2>
RectangularSparseMat<double,dim,dim2> *SubDomain::create_FaceBaseddRdXoperators()
{
  Connectivity *faceToNode = createFaceToNodeConnectivity();
  int numFaces = faces.size();
  int *ia = (*faceToNode).ptr();
  int *ja = (*faceToNode)[0];
  double (*a)[dim*dim2] = 0;

  RectangularSparseMat<double, dim, dim2> *A = new RectangularSparseMat<double, dim, dim2>(numFaces, ia[numFaces], ia, ja, a, 0, 0);
  return A;
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
DiagMat<Scalar,dim> *SubDomain::createMaskDiagonal(typename DiagMat<Scalar,dim>::Type type,
						   int *ndType)
{

  DiagMat<Scalar,dim> *A = new DiagMat<Scalar,dim>(type, nodes.size(), ndType);

  return A;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
SparseMat<Scalar,dim> *SubDomain::createMaskILU(int fill, int renum, int *ndType)
{

  nodeToNodeMaskILU = createEdgeBasedConnectivity();

  compStruct *nodeRenum = createRenumbering(nodeToNodeMaskILU, renum, 0);

  int *ia = (*nodeToNodeMaskILU).ptr();
  int *ja = (*nodeToNodeMaskILU)[0];
  int n = nodes.size();
  int nnz = ia[n];

  SparseMat<Scalar,dim> *A = new SparseMat<Scalar,dim>(n, nnz, ia, ja, 0, nodeRenum, ndType);

  A->symbolicILU(fill);

  A->createPointers(edges);

  return A;

}

//------------------------------------------------------------------------------
//
template<class Scalar, int dim>
void SubDomain::computeH1(FluxFcn **fluxFcn, BcData<dim> &bcData,
                          GeoState &geoState, Vec<double> &ctrlVol,
                          SVec<double,dim> &V, GenMat<Scalar,dim> &A)
{

  int k;

  double dfdUi[dim*dim], dfdUj[dim*dim];

  Scalar *Aii, *Ajj, *Aij, *Aji;

  // contribution of the edges

  Vec<Vec3D> &edgeNorm = geoState.getEdgeNormal();
  Vec<double> &edgeNormVel = geoState.getEdgeNormalVel();

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); ++l) {

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    fluxFcn[0]->computeJacobians(1.0, 0.0, edgeNorm[l], edgeNormVel[l],
                                 V[i], V[j], dfdUi, dfdUj);

    Aii = A.getElem_ii(i);
    Ajj = A.getElem_ii(j);
    Aij = A.getElem_ij(l);
    Aji = A.getElem_ji(l);

    if (edgeFlag[l])
      for (k=0; k<dim*dim; ++k) { Aii[k] += dfdUi[k]; Ajj[k] -= dfdUj[k]; }

    if (Aij && Aji) {

      double voli = 1.0 / ctrlVol[i];
      double volj = 1.0 / ctrlVol[j];
      for (k=0; k<dim*dim; ++k) {
        Aij[k] += dfdUj[k] * voli;
        Aji[k] -= dfdUi[k] * volj;
      }
    }

  }
  // contribution of the boundary faces
  faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A);

  for (int i=0; i<ctrlVol.size(); ++i) {

    double voli = 1.0 / ctrlVol[i];

    Aii = A.getElem_ii(i);

    for (k=0; k<dim*dim; ++k) Aii[k] *= voli;

  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void SubDomain::computeH2(FluxFcn **fluxFcn, RecFcn *recFcn, BcData<dim> &bcData,
			  GeoState &geoState, SVec<double,3> &X, SVec<double,dim> &V,
			  NodalGrad<dim> &ngrad, GenMat<Scalar,neq> &A)
{

  //std::cout << "$$$$$ IN SUBDOMAIN computeH2\n";
  
  double ddVij[dim], ddVji[dim], Vi[dim], Vj[dim], dfdVi[dim*dim], dfdVj[dim*dim];

  Scalar *Aij, *Aji;

  // contribution of the edges

  Vec<Vec3D> &edgeNorm = geoState.getEdgeNormal();
  Vec<double> &edgeNormVel = geoState.getEdgeNormalVel();

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  SVec<double,dim> &dVdx = ngrad.getX();
  SVec<double,dim> &dVdy = ngrad.getY();
  SVec<double,dim> &dVdz = ngrad.getZ();

  for (int l=0; l<edges.size(); ++l) {

    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    int k;
    for (k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
    }

    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);

    fluxFcn[BC_INTERNAL]->computeJacobians(1.0, 0.0, edgeNorm[l], edgeNormVel[l], Vi, Vj, dfdVi, dfdVj);

    Aij = A.getElem_ij(l);
    Aji = A.getElem_ji(l);

    if (Aij && Aji)  {
      for (k=0; k<neq*neq; ++k) {
        Aij[k] += dfdVj[k];
        Aji[k] += dfdVi[k];
      }
    }

  }

  // contribution of the boundary faces
  faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void SubDomain::computeH2transpose(FluxFcn **fluxFcn, RecFcn *recFcn, BcData<dim> &bcData,
                                   GeoState &geoState, SVec<double,3> &X, SVec<double,dim> &V,
                                   NodalGrad<dim> &ngrad, GenMat<Scalar,neq> &A)
{

  double ddVij[dim], ddVji[dim], Vi[dim], Vj[dim], dfdVi[dim*dim], dfdVj[dim*dim];

  Scalar *Aij, *Aji;

  // contribution of the edges

  Vec<Vec3D> &edgeNorm = geoState.getEdgeNormal();
  Vec<double> &edgeNormVel = geoState.getEdgeNormalVel();

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  SVec<double,dim> &dVdx = ngrad.getX();
  SVec<double,dim> &dVdy = ngrad.getY();
  SVec<double,dim> &dVdz = ngrad.getZ();

  for (int l=0; l<edges.size(); ++l) {

    if (!edgeFlag[l]) continue;

    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    int k;
    for (k=0; k<dim; ++k) {
      ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
      ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
    }

    recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);

    fluxFcn[BC_INTERNAL]->computeJacobians(1.0, 0.0, edgeNorm[l], edgeNormVel[l], Vi, Vj, dfdVi, dfdVj);

    Aij = A.getElem_ij(l);
    Aji = A.getElem_ji(l);

    if (Aij && Aji)  {
      for (k=0; k<neq*neq; ++k) {
        Aji[k] += dfdVj[k];
        Aij[k] += dfdVi[k];
      }
    }

  }

  // contribution of the boundary faces
  faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void SubDomain::computeH2(FluxFcn **fluxFcn, RecFcn *recFcn, BcData<dim> &bcData,
			  GeoState &geoState, SVec<double,3> &X, SVec<double,dim> &V,
			  NodalGrad<dim> &ngrad, 
			  ExactRiemannSolver<dim>& riemann,
			  LevelSetStructure &LSS, 
			  Vec<int> &fluidId, int Nriemann,
			  GenMat<Scalar,neq> &A,
			  SVec<double,dim> &aij, SVec<double,dim> &aji,
			  SVec<double,dim> &bij, SVec<double,dim> &bji,
			  SVec<double,dim> &betaij, SVec<double,dim> &betaji)
{

  //std::cout << "$$$$$ IN SUBDOMAIN EMB computeH2\n";

  double  Vi[2*dim],  Vj[2*dim], Vstar[2*dim];

  double ddVij[dim], ddVji[dim];

  double  dfdVi[dim*dim],  dfdVj[dim*dim], dVsdV[dim*dim];
  double dflux1[dim*dim], dflux2[dim*dim];
  double betai[dim], betaj[dim];

  Scalar *Aij, *Aji;

  Vec<Vec3D>     &edgeNorm = geoState.getEdgeNormal();
  Vec<double> &edgeNormVel = geoState.getEdgeNormalVel();

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  SVec<double,dim> &dVdx = ngrad.getX();
  SVec<double,dim> &dVdy = ngrad.getY();
  SVec<double,dim> &dVdz = ngrad.getZ();

  Vec3D normalDir;

  int k;

  int farfieldFluid = 0;
  double d_gradPhi;

   double alpha = 0.1;

   VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();

   for (int l=0; l<edges.size(); ++l) {

     if (!edgeFlag[l]) continue;

     int i = edgePtr[l][0];
     int j = edgePtr[l][1];

     bool intersect = LSS.edgeIntersectsStructure(0,l);

     bool iActive = LSS.isActive(0.0,i);
     bool jActive = LSS.isActive(0.0,j);

     double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
     double length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

     for (k=0; k<dim; ++k) {
       ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
       ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
     }

     if (iActive && jActive && !intersect){
       recFcn->compute(V[i], ddVij, V[j], ddVji, Vi, Vj);
     } else {
       for(k=0; k<dim; k++) {
	 Vi[k] = V[i][k];
	 Vj[k] = V[j][k];
       }
     }

     for(k=0; k<dim; k++) {
       Vi[k+dim] = V[i][k];
       Vj[k+dim] = V[j][k];
     }

     if (!intersect) {

       fluxFcn[BC_INTERNAL]->computeJacobians(1.0, 0.0, edgeNorm[l], edgeNormVel[l], Vi, Vj, dfdVi, dfdVj);

       Aij = A.getElem_ij(l);
       Aji = A.getElem_ji(l);

       if (Aij && Aji)  {
	 for (k=0; k<dim*dim; ++k) {
	   Aij[k] += dfdVj[k];
	   Aji[k] += dfdVi[k];
	 }
       }

     } else {

       if(iActive) {

	 LevelSetResult resij = LSS.getLevelSetDataAtEdgeCenter(0.0, l, true);
	 switch (Nriemann) {
	 case 0: //structure normal
	   d_gradPhi = dx[0]*resij.gradPhi[0]+dx[1]*resij.gradPhi[1]+dx[2]*resij.gradPhi[2];
	   normalDir = (d_gradPhi>=0.0) ? -1.0*resij.gradPhi : resij.gradPhi;
	   break;
	 case 1: //fluid normal
	   normalDir = -1.0/(edgeNorm[l].norm())*edgeNorm[l];
	   break;
	 default:
	   fprintf(stderr,"ERROR: Unknown RiemannNormal code!\n");
	   exit(-1);
	 }

	 for (k=0; k<dim; ++k) betai[k] = 0.0;

	 //*************************************
	 if (higherOrderFSI) {

	   double ri[dim];
	   higherOrderFSI->estimateR(l, 0, i, V, ngrad, X, fluidId, ri);
	  
	   for (k=0; k<dim; ++k) betai[k] = 1.0;

	   if (higherOrderFSI->limitExtrapolation()) {
	     if (V[i][1]*dx[0]+V[i][2]*dx[1]+V[i][3]*dx[2] < 0.0) {
	       for (int k = 0; k < dim; ++k) {
		 betai[k] = std::min<double>(betai[k],ri[k]);
	       }
	     }
	   }

	   for (k=0; k<dim; ++k){
	     Vi[k] = V[i][k] + (1.0 - resij.alpha)*ddVij[k]*betai[k];
	   }

	 }
	 //*************************************

	 riemann.computeFSIRiemannSolution(Vi, resij.normVel, normalDir, varFcn, Vstar, j,        fluidId[i]);
	 riemann.computeFSIRiemannJacobian(Vi, resij.normVel, normalDir, varFcn, Vstar, j, dVsdV, fluidId[i]);

	 //*************************************
	 if(higherOrderFSI) {

	   V6NodeData (*v6data)[2] = higherOrderFSI->getV6Data();

	   if (v6data == NULL) {
	     for (int k=0; k<dim; k++) {
	       Vstar[k] = V[i][k] + (0.5/max(1.0-resij.alpha, alpha)) * (Vstar[k] - V[i][k]);
	     }
	   } else {
	     higherOrderFSI->extrapolateV6(l, 0, i, V, Vi, Vstar, X, resij.alpha, length, fluidId, betai);
	   }

	 }
	 //*************************************

	 for (int k=0; k<dim; k++) {
	   aij[l][k] = Vi[k];
	   bij[l][k] = Vstar[k];
	   betaij[l][k] = betai[k];
	 }

	 Aji = A.getElem_ji(l);
	 //Aij = A.getElem_ij(l);
	 for (int k=0; k<dim*dim; ++k) {
	   Aji[k] += dVsdV[k];
	   //Aij[k] += 0.0;
	 }	

       }

       if(jActive){

	 LevelSetResult resji = LSS.getLevelSetDataAtEdgeCenter(0.0, l, false);
	 switch (Nriemann) {
	 case 0: //structure normal
	   d_gradPhi = dx[0]*resji.gradPhi[0]+dx[1]*resji.gradPhi[1]+dx[2]*resji.gradPhi[2];
	   normalDir = (d_gradPhi>=0.0) ? resji.gradPhi : -1.0*resji.gradPhi;
	   break;
	 case 1: //fluid normal
	   normalDir = 1.0/(edgeNorm[l].norm())*edgeNorm[l];
	   break;
	 default:
	   fprintf(stderr,"ERROR: Unknown RiemannNormal code!\n");
	   exit(-1);
	 }

	 for (k=0; k<dim; ++k) betaj[k] = 0.0;

	 //*************************************
	 if (higherOrderFSI) {

	   double rj[dim];

	   higherOrderFSI->estimateR(l, 1, j, V, ngrad, X, fluidId, rj);

	   for (k=0; k<dim; ++k) betaj[k] = 1.0;

	   if (higherOrderFSI->limitExtrapolation()) {
	     if (V[j][1]*dx[0]+V[j][2]*dx[1]+V[j][3]*dx[2] > 0.0) {
	       for (int k = 0; k < dim; ++k) {
		 betaj[k] = std::min<double>(betaj[k],rj[k]);
	       }
	     }
	   }
	   for (k=0; k<dim; ++k){
	     Vj[k] = V[j][k]-(1.0-resji.alpha)*ddVji[k]*betaj[k];
	   }

	 }
	 //*************************************

	 riemann.computeFSIRiemannSolution(Vj, resji.normVel, normalDir, varFcn, Vstar, i,        fluidId[j]);
	 riemann.computeFSIRiemannJacobian(Vj, resji.normVel, normalDir, varFcn, Vstar, i, dVsdV, fluidId[j]);

	 //*************************************
	 if (higherOrderFSI) {

	   V6NodeData (*v6data)[2] = higherOrderFSI->getV6Data();

	   if (v6data==NULL) {
	     for (int k=0; k<dim; k++){
	       Vstar[k] = V[j][k]+(0.5/max(1.0 - resji.alpha, alpha))*(Vstar[k] - V[j][k]);
	     }
	   } else {
	     higherOrderFSI->extrapolateV6(l, 1, j, V, Vj, Vstar, X, 1.0-resji.alpha, length, fluidId, betaj);
	   }

	 }
	 //*************************************

	 for (int k=0; k<dim; k++) {
	   aji[l][k] = Vstar[k];
	   bji[l][k] = Vj[k];
	   betaji[l][k] = betaj[k];
	 }

	 Aij = A.getElem_ij(l);
	 //Aji = A.getElem_ji(l);
	 for (int k=0; k<dim*dim; ++k) {
	   Aij[k] += dVsdV[k];
	   //Aji[k] += 0.0;
	 }

       }    

     }

   }

   // contribution of the boundary faces
   faces.computeJacobianFiniteVolumeTerm(fluxFcn, bcData, geoState, V, A);

 }

 //------------------------------------------------------------------------------

 template<class Scalar, int dim>
 void SubDomain::precomputeRec(RecFcn *recFcn, SVec<double,3> &X,
			       SVec<double,dim> &V, NodalGrad<dim> &ngrad,
			       SVec<Scalar,dim> &aij, SVec<Scalar,dim> &aji,
			       SVec<Scalar,dim> &bij, SVec<Scalar,dim> &bji)
 {

   double ddVij[dim], ddVji[dim];

   SVec<double,dim> &dVdx = ngrad.getX();
   SVec<double,dim> &dVdy = ngrad.getY();
   SVec<double,dim> &dVdz = ngrad.getZ();

   int (*edgePtr)[2] = edges.getPtr();

   for (int l=0; l<edges.size(); ++l) {

     int i = edgePtr[l][0];
     int j = edgePtr[l][1];

     double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
     for (int k=0; k<dim; ++k) {
       ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
       ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
     }

     recFcn->precompute(V[i], ddVij, V[j], ddVji, aij[l], aji[l], bij[l], bji[l]);

   }

 }

 //------------------------------------------------------------------------------

 template<class Scalar, int dim>
 void SubDomain::precomputeRec(RecFcn *recFcn, SVec<double,3> &X,
			       SVec<double,dim> &V, NodalGrad<dim> &ngrad,
			       LevelSetStructure &LSS, Vec<int> &fluidId,
			       SVec<Scalar,dim> &aij, SVec<Scalar,dim> &aji,
			       SVec<Scalar,dim> &bij, SVec<Scalar,dim> &bji)
 {

   //std::cout << "$$$$$ IN SUBDOMAIN EMB precompRec\n";

   double ddVij[dim], ddVji[dim];

   SVec<double,dim> &dVdx = ngrad.getX();
   SVec<double,dim> &dVdy = ngrad.getY();
   SVec<double,dim> &dVdz = ngrad.getZ();

   int (*edgePtr)[2] = edges.getPtr();

   for (int l=0; l<edges.size(); ++l) {

     int i = edgePtr[l][0];
     int j = edgePtr[l][1];

     bool intersect = LSS.edgeIntersectsStructure(0,l);

     bool iActive = LSS.isActive(0.0,i);
     bool jActive = LSS.isActive(0.0,j);

     double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};

     if (iActive && jActive && !intersect){
 
       for (int k=0; k<dim; ++k) {
	 ddVij[k] = dx[0]*dVdx[i][k] + dx[1]*dVdy[i][k] + dx[2]*dVdz[i][k];
	 ddVji[k] = dx[0]*dVdx[j][k] + dx[1]*dVdy[j][k] + dx[2]*dVdz[j][k];
       }
       recFcn->precompute(V[i], ddVij, V[j], ddVji, aij[l], aji[l], bij[l], bji[l]);

     }

   }

 }

 //------------------------------------------------------------------------------

 template<class Scalar, int dim>
 void SubDomain::computeMatVecProdH1(bool *nodeFlag, GenMat<Scalar,dim> &A,
				     SVec<double,dim> &p, SVec<double,dim> &prod)
 {

   int i, j, l;

   int numNodes = nodes.size();
   int numEdges = edges.size();

   bool *edgeFlag = edges.getMasterFlag();
   int (*edgePtr)[2] = edges.getPtr();

   Scalar (*a)[dim*dim] = A.data();

   prod = 0.0;

 #pragma ivdep
   for (i=0; i<numNodes; ++i)
     if (nodeFlag[i])
       DenseMatrixOp<Scalar,dim,dim*dim>::applyAndAddToVector(a, i, p.v, i, prod.v, i);

 #pragma ivdep
   for (l=0; l<numEdges; ++l) {

     if (edgeFlag[l]) {

       i = edgePtr[l][0];
       j = edgePtr[l][1];

       DenseMatrixOp<Scalar,dim,dim*dim>::applyAndAddToVector(a, numNodes + 2*l, p.v, j, prod.v, i);
       DenseMatrixOp<Scalar,dim,dim*dim>::applyAndAddToVector(a, numNodes + 2*l + 1, p.v, i, prod.v, j);

     }

   }

 }

 //------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::computeMatVecProdH1transpose(bool *nodeFlag, GenMat<Scalar,dim> &A,
                                             SVec<double,dim> &p, SVec<double,dim> &prod)
{

  int i, j, l;

  int numNodes = nodes.size();
  int numEdges = edges.size();

  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  Scalar (*a)[dim*dim] = A.data();

  prod = 0.0;

#pragma ivdep
  for (i=0; i<numNodes; ++i)
    if (nodeFlag[i])
      DenseMatrixOp<Scalar,dim,dim*dim>::applyTransAndAddToVector(a, i, p.v, i, prod.v, i);

#pragma ivdep
  for (l=0; l<numEdges; ++l) {

    if (edgeFlag[l]) {

      i = edgePtr[l][0];
      j = edgePtr[l][1];

      DenseMatrixOp<Scalar,dim,dim*dim>::applyTransAndAddToVector(a, numNodes + 2*l, p.v, i, prod.v, j);
      DenseMatrixOp<Scalar,dim,dim*dim>::applyTransAndAddToVector(a, numNodes + 2*l + 1, p.v, j, prod.v, i);

    }

  }

}

 //------------------------------------------------------------------------------

 template<class Scalar, int dim>
 void SubDomain::computeMatVecProdH1(bool *nodeFlag, GenMat<Scalar,dim> &A,
				     SVec<double,dim> &p, SVec<double,dim> &prod,
				     SVec<double,dim> &ghostP, SVec<double,dim>& ghostProd)
 {

   int i, j, l;

   int numNodes = nodes.size();
   int numEdges = edges.size();

   bool *edgeFlag = edges.getMasterFlag();
   int (*edgePtr)[2] = edges.getPtr();

   Scalar (*a)[dim*dim] = A.data();

   prod = 0.0;

 #pragma ivdep
   for (i=0; i<numNodes; ++i)
     if (nodeFlag[i])
       DenseMatrixOp<Scalar,dim,dim*dim>::applyAndAddToVector(a, i, p.v, i, prod.v, i);

 #pragma ivdep
   for (l=0; l<numEdges; ++l) {

     if (edgeFlag[l]) {

       i = edgePtr[l][0];
       j = edgePtr[l][1];

       DenseMatrixOp<Scalar,dim,dim*dim>::applyAndAddToVector(a, numNodes + 2*l,     p.v, j, prod.v, i);
       DenseMatrixOp<Scalar,dim,dim*dim>::applyAndAddToVector(a, numNodes + 2*l + 1, p.v, i, prod.v, j);

     }

   }

   typename GenMat<Scalar,dim>::AuxilliaryIterator* myItr = A.begin_realNodes();
   if (myItr) {
     do { 
       DenseMatrixOp<Scalar,dim,dim*dim>::applyAndAddToVector(reinterpret_cast<Scalar(*)[dim*dim]>(myItr->pData), 0, ghostP.v, myItr->col , prod.v, myItr->row);
     } while (A.next(myItr));
     A.free(myItr);
   }

   myItr = A.begin_ghostNodes();
   if (myItr) {
     do { 
       DenseMatrixOp<Scalar,dim,dim*dim>::applyAndAddToVector(reinterpret_cast<Scalar(*)[dim*dim]>(myItr->pData), 0, p.v, myItr->col , ghostProd.v, myItr->row);
     } while (A.next(myItr));
     A.free(myItr);
   }

   myItr = A.begin_ghostGhostNodes();
   if (myItr) {
     do { 
      DenseMatrixOp<Scalar,dim,dim*dim>::applyAndAddToVector(reinterpret_cast<Scalar(*)[dim*dim]>(myItr->pData), 0, ghostP.v, myItr->col , ghostProd.v, myItr->row);
     } while (A.next(myItr));
     A.free(myItr);
   }

 }


 template<class Scalar, int dim>
 void SubDomain::
 computeMatVecProdH1FarFieldHH(bool *nodeFlag, GenMat<Scalar,dim> &A, SVec<double,dim> &p_u,
			       SVec<double,dim> &prod_u,Vec<double>& p_hh, 
			       Vec<double>& prod_hh) {

   faces.computeMatVecProdH1FarFieldHH(A, p_u, prod_u,p_hh, prod_hh);

 }

 //------------------------------------------------------------------------------

 template<class Scalar1, class Scalar2, int dim>
 void SubDomain::computeMatVecProdH2(RecFcn *recFcn, SVec<double,3> &X,
				     Vec<double> &ctrlVol, GenMat<Scalar1,dim> &A,
				     SVec<double,dim> &aij, SVec<double,dim> &aji,
				     SVec<double,dim> &bij, SVec<double,dim> &bji,
				     SVec<Scalar2,dim> &p, NodalGrad<dim, Scalar2> &dpdxj,
				     SVec<Scalar2,dim> &prod) {

   int i, j, l;

   Scalar2 ddpij[dim], ddpji[dim], pij[1][dim], pji[1][dim];
   Scalar2 tmp[1][dim], tmpi[1][dim], tmpj[1][dim];

   Scalar1 (*a)[dim*dim] = A.data();

   SVec<Scalar2,dim> &dpdx = dpdxj.getX();
   SVec<Scalar2,dim> &dpdy = dpdxj.getY();
   SVec<Scalar2,dim> &dpdz = dpdxj.getZ();

   prod = (Scalar2) 0.0;

   int numNodes = nodes.size();
   int numEdges = edges.size();

   int (*edgePtr)[2] = edges.getPtr();

   bool *masterFlag = edges.getMasterFlag();

   if (bcMap.size() > 0)  {

     int index;
     for (l=0; l<numEdges; ++l) {
       if (!masterFlag[l]) continue;
       i = edgePtr[l][0];
       j = edgePtr[l][1];

       double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
       for (int k=0; k<dim; ++k) {
	 ddpij[k] = dx[0]*dpdx[i][k] + dx[1]*dpdy[i][k] + dx[2]*dpdz[i][k];
	 ddpji[k] = dx[0]*dpdx[j][k] + dx[1]*dpdy[j][k] + dx[2]*dpdz[j][k];

       }

       // result of the reconstructed-limited states are in pij, pji
       recFcn->template compute<Scalar2, dim>(p[i], ddpij, p[j], ddpji, aij[l], aji[l],
					      bij[l], bji[l], pij[0], pji[0]);

       // A is applied to reconstructed-limited states and stored in tmpi, tmpj
       // address of a is shifted by the number of diagonal entries (numnodes)

       if (bcMap.find(l) != bcMap.end())  {
	 if (nodeType[i] != BC_INTERNAL)
	   index = numNodes+2*numEdges+2*bcMap[l];
	 else
	   index = numNodes + 2*l;

	 DenseMatrixOp<Scalar1,dim,dim*dim>::applyAndAddToVector(a, index, pji, 0, prod.v, i);
	 DenseMatrixOp<Scalar1,dim,dim*dim>::applyAndAddToVector(a, index+1, pij, 0, prod.v, i);

	 if (nodeType[j] != BC_INTERNAL)
	   index = numNodes+2*numEdges+2*bcMap[l]+2*numBcNodes[l];
	 else
	   index = numNodes+2*l;

	 DenseMatrixOp<Scalar1,dim,dim*dim>::applyAndSubToVector(a, index, pji, 0, prod.v, j);
	 DenseMatrixOp<Scalar1,dim,dim*dim>::applyAndSubToVector(a, index+1, pij, 0, prod.v, j);

       }
       else  {
	 DenseMatrixOp<Scalar1,dim,dim*dim>::applyToVector(a, numNodes + 2*l, pji, 0, tmpi, 0);
	 DenseMatrixOp<Scalar1,dim,dim*dim>::applyToVector(a, numNodes + 2*l + 1, pij, 0, tmpj, 0);

	 VectorOp<Scalar2,dim>::sum(tmpi, 0, tmpj, 0, tmp, 0);
	 VectorOp<Scalar2,dim>::add(tmp, 0, prod.v, i);
	 VectorOp<Scalar2,dim>::sub(tmp, 0, prod.v, j);
       }
     }
   }
   else  {

     for (l=0; l<numEdges; ++l) {
       if (!masterFlag[l]) continue;
       i = edgePtr[l][0];
       j = edgePtr[l][1];

       double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
       for (int k=0; k<dim; ++k) {
	 ddpij[k] = dx[0]*dpdx[i][k] + dx[1]*dpdy[i][k] + dx[2]*dpdz[i][k];
	 ddpji[k] = dx[0]*dpdx[j][k] + dx[1]*dpdy[j][k] + dx[2]*dpdz[j][k];

       }

       // result of the reconstructed-limited states are in pij, pji
       recFcn->template compute<Scalar2, dim>(p[i], ddpij, p[j], ddpji, aij[l], aji[l],
					      bij[l], bji[l], pij[0], pji[0]);

       DenseMatrixOp<Scalar1,dim,dim*dim>::applyToVector(a, numNodes + 2*l,     pji, 0, tmpi, 0);
       DenseMatrixOp<Scalar1,dim,dim*dim>::applyToVector(a, numNodes + 2*l + 1, pij, 0, tmpj, 0);

       VectorOp<Scalar2,dim>::sum(tmpi, 0, tmpj, 0, tmp, 0);
       VectorOp<Scalar2,dim>::add(tmp, 0, prod.v, i); 
       VectorOp<Scalar2,dim>::sub(tmp, 0, prod.v, j); 
     }
   }

   // contribution from diagonal entries of A
   for (i=0; i<numNodes; ++i) {

     DenseMatrixOp<Scalar1,dim,dim*dim>::applyAndAddToVector(a, i, p.v, i, prod.v, i);

     double voli = 1.0 / ctrlVol[i];
     for (int k=0; k<dim; ++k) prod[i][k] *= voli;

   }
 }

 //------------------------------------------------------------------------------

template<class Scalar1, class Scalar2, int dim>
void SubDomain::computeMatVecProdH2transposeNew(IoData& iod, SVec<double,3> &X,
                                                Vec<double> &ctrlVol, GenMat<Scalar1,dim> &A,
                                                SVec<double,dim> &aij, SVec<double,dim> &aji,
                                                SVec<double,dim> &bij, SVec<double,dim> &bji,
                                                SVec<Scalar2,dim> &cij, SVec<Scalar2,dim> &cji,
                                                SVec<Scalar2,dim> &dij, SVec<Scalar2,dim> &dji,
                                                NodalGrad<dim, Scalar2> &dpdxj,
                                                SVec<Scalar2,dim> &p, SVec<Scalar2,dim> &prod) {

  int i, j, l, k;
  double voli, volj;

  Scalar2 tmpi[1][dim], tmpj[1][dim];
  Scalar1 (*a)[dim*dim] = A.data();

  prod = (Scalar2) 0.0;

  int numNodes = nodes.size();
  int numEdges = edges.size();

  int (*edgePtr)[2] = edges.getPtr();
  bool *masterFlag = edges.getMasterFlag();

  if(iod.schemes.ns.reconstruction == SchemeData::LINEAR) {
    double a, b, c, d;
    for (l=0; l<numEdges; ++l) {
      if (!masterFlag[l]) continue;
      i = edgePtr[l][0];
      j = edgePtr[l][1];
      voli = 1.0 / ctrlVol[i];
      volj = 1.0 / ctrlVol[j];

      for (k=0; k<dim; ++k) {
        a = aji[l][k]*cij[l][k];
        b = aij[l][k]*cji[l][k];
        c = aji[l][k]*dij[l][k];
        d = aij[l][k]*dji[l][k];
        prod[i][k] += voli*(a + cji[l][k] - b) - volj*(c + dji[l][k] - d);
        prod[j][k] += voli*(cij[l][k] - a + b) - volj*(dij[l][k] - c + d);
      }
    }
  } else {
    for (l=0; l<numEdges; ++l) {
      if (!masterFlag[l]) continue;
      i = edgePtr[l][0];
      j = edgePtr[l][1];
      voli = 1.0 / ctrlVol[i];
      volj = 1.0 / ctrlVol[j];

      DenseMatrixOp<Scalar1,dim,dim*dim>::applyTransToVector(a, numNodes + 2*l, p.v, i, tmpi, 0);
      DenseMatrixOp<Scalar1,dim,dim*dim>::applyTransToVector(a, numNodes + 2*l + 1, p.v, i, tmpj, 0);
      for (k=0; k<dim; ++k) {
        prod[i][k] += voli*tmpj[0][k];
        prod[j][k] += voli*tmpi[0][k];
      }

      DenseMatrixOp<Scalar1,dim,dim*dim>::applyTransToVector(a, numNodes + 2*l, p.v, j, tmpi, 0);
      DenseMatrixOp<Scalar1,dim,dim*dim>::applyTransToVector(a, numNodes + 2*l + 1, p.v, j, tmpj, 0);
      for (k=0; k<dim; ++k) {
        prod[i][k] += -volj*tmpj[0][k];
        prod[j][k] += -volj*tmpi[0][k];
      }
    }
  }
}

 //------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::addToMatVecProdH2transposeGalerkinNew(Vec<double> &ctrlVol, SVec<Scalar,dim> &ddxt,
                                                      SVec<Scalar,dim> &ddyt, SVec<Scalar,dim> &ddzt,
                                                      NodalGrad<dim, Scalar> &dpdxj, SVec<Scalar,dim> &prod) {
  int i, j;
  double voli, volj;

  SVec<double,3> wii = dpdxj.getWii();
  SVec<double,3> wij = dpdxj.getWij();
  SVec<double,3> wji = dpdxj.getWji();

  for (i=0; i<nodes.size(); ++i) {
    double coef = 1.0 / (4.0*ctrlVol[i]);
    for(int k=0; k<dim; ++k) {
      prod[i][k] += coef*(ddxt[i][k]*wii[i][0] + ddyt[i][k]*wii[i][1] + ddzt[i][k]*wii[i][2]);
    }
  }

  int numEdges = edges.size();

  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<numEdges; ++l) {
    i = edgePtr[l][0];
    j = edgePtr[l][1];
    voli = 1.0/(4.0*ctrlVol[i]);
    volj = 1.0/(4.0*ctrlVol[j]);

    for (int k=0; k<dim; ++k)  {
      prod[i][k] += volj*(ddxt[j][k]*wji[l][0] + ddyt[j][k]*wji[l][1] + ddzt[j][k]*wji[l][2]);
      prod[j][k] += voli*(ddxt[i][k]*wij[l][0] + ddyt[i][k]*wij[l][1] + ddzt[i][k]*wij[l][2]);
    }
  }
}

 //------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::addToMatVecProdH2transposeLeastSquareNew(SVec<double,3> &X, SVec<double,6> &R, SVec<Scalar,dim> &ddxt,
                                                         SVec<Scalar,dim> &ddyt, SVec<Scalar,dim> &ddzt,
                                                         NodalGrad<dim, Scalar> &dpdxj, SVec<Scalar,dim> &prod) {
  int i, j;
  double Wi[3], Wj[3], a, b;
  int numEdges = edges.size();
  int (*edgePtr)[2] = edges.getPtr();
  bool *edgeFlag = edges.getMasterFlag();

  for (int l=0; l<numEdges; ++l) {
    if (!edgeFlag[l]) continue;
    i = edgePtr[l][0];
    j = edgePtr[l][1];

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
    computeLocalWeightsLeastSquares(dx, R[i], Wi);
    dx[0] = -dx[0]; dx[1] = -dx[1]; dx[2] = -dx[2];
    computeLocalWeightsLeastSquares(dx, R[j], Wj);

    for (int k=0; k<dim; ++k)  {
      a = ddxt[j][k]*Wj[0] + ddyt[j][k]*Wj[1] + ddzt[j][k]*Wj[2];
      b = ddxt[i][k]*Wi[0] + ddyt[i][k]*Wi[1] + ddzt[i][k]*Wi[2];
      prod[i][k] += a-b;
      prod[j][k] += b-a;
    }
  }
}

 //------------------------------------------------------------------------------

template<class Scalar1, class Scalar2, int dim>
void SubDomain::computeGradientsTransposeNew(SVec<double,3> &X,
                                             Vec<double> &ctrlVol, GenMat<Scalar1,dim> &A,
                                             SVec<double,dim> &bij, SVec<double,dim> &bji,
                                             SVec<Scalar2,dim> &cij, SVec<Scalar2,dim> &cji,
                                             SVec<Scalar2,dim> &dij, SVec<Scalar2,dim> &dji,
                                             SVec<Scalar2,dim> &p, SVec<Scalar2, dim> &ddxt,
                                             SVec<Scalar2, dim> &ddyt, SVec<Scalar2, dim> &ddzt) {

  int i, j, l, k;
  double b, c, d, e;

  Scalar1 (*a)[dim*dim] = A.data();

  ddxt = (Scalar2) 0.0;
  ddyt = (Scalar2) 0.0;
  ddzt = (Scalar2) 0.0;
  cij = (Scalar2) 0.0;
  cji = (Scalar2) 0.0;
  dij = (Scalar2) 0.0;
  dji = (Scalar2) 0.0;

  int numNodes = nodes.size();
  int numEdges = edges.size();

  int (*edgePtr)[2] = edges.getPtr();
  bool *masterFlag = edges.getMasterFlag();

  for (l=0; l<numEdges; ++l) {
    if (!masterFlag[l]) continue;
    i = edgePtr[l][0];
    j = edgePtr[l][1];
    double voli = 1.0 / ctrlVol[i];
    double volj = 1.0 / ctrlVol[j];

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};

    DenseMatrixOp<Scalar1,dim,dim*dim>::applyTransToVector(a, numNodes + 2*l, p.v, i, cij.v, l);
    DenseMatrixOp<Scalar1,dim,dim*dim>::applyTransToVector(a, numNodes + 2*l + 1, p.v, i, cji.v, l);
    DenseMatrixOp<Scalar1,dim,dim*dim>::applyTransToVector(a, numNodes + 2*l, p.v, j, dij.v, l);
    DenseMatrixOp<Scalar1,dim,dim*dim>::applyTransToVector(a, numNodes + 2*l + 1, p.v, j, dji.v, l);
    for(k=0; k<dim; ++k) {
      b = bij[l][k]*(voli*cji[l][k] - volj*dji[l][k]);
      c = bji[l][k]*(-voli*cij[l][k] + volj*dij[l][k]);
      ddxt[i][k] += b*dx[0];  ddyt[i][k] += b*dx[1];   ddzt[i][k] += b*dx[2];
      ddxt[j][k] += c*dx[0];  ddyt[j][k] += c*dx[1];   ddzt[j][k] += c*dx[2];
    }
  }
}

 //------------------------------------------------------------------------------

template<class Scalar1, class Scalar2, int dim>
void SubDomain::addDiagonalInMatVecProdH2transpose(Vec<double> &ctrlVol, GenMat<Scalar1,dim> &A,
                                                   SVec<Scalar2,dim> &p, SVec<Scalar2,dim> &prod) {

  Scalar1 (*a)[dim*dim] = A.data();
  int numNodes = nodes.size();
  // contribution from diagonal entries of A
  for (int i=0; i<numNodes; ++i) {
    double voli = 1.0 / ctrlVol[i];
    for (int k=0; k<dim; ++k) p[i][k] *= voli;
    DenseMatrixOp<Scalar1,dim,dim*dim>::applyTransAndAddToVector(a, i, p.v, i, prod.v, i);
    for (int k=0; k<dim; ++k) p[i][k] /= voli;
  }

}

 //------------------------------------------------------------------------------

 template<class Scalar1, class Scalar2, int dim>
 void SubDomain::computeMatVecProdH2(FluxFcn **fluxFcn, RecFcn *recFcn, GeoState &geoState,
				     SVec<double,3> &X, Vec<double> &ctrlVol, 
				     ExactRiemannSolver<dim>& riemann,
				     LevelSetStructure &LSS,
				     Vec<int> &fluidId, int Nriemann,
				     GenMat<Scalar1,dim> &A,
				     SVec<double,dim> &aij, SVec<double,dim> &aji,
				     SVec<double,dim> &bij, SVec<double,dim> &bji,
				     SVec<double,dim> &betaij, SVec<double,dim> &betaji,
				     SVec<Scalar2,dim> &p, NodalGrad<dim, Scalar2> &dpdxj,
				     SVec<Scalar2,dim> &prod) {

   int i, j, l;

   Scalar2 ddpij[dim], ddpji[dim], pij[1][dim], pji[1][dim];
   Scalar2 tmp[1][dim], tmpi[1][dim], tmpj[1][dim];

   Scalar2 Vi[2*dim], Vj[2*dim], Vstar[2*dim];

   Scalar2 dfdVi[dim*dim],  dfdVj[dim*dim];
   Scalar2 dVsdV[dim*dim], dflux1[dim*dim], dflux2[dim*dim];

   Scalar1 (*a)[dim*dim] = A.data();

   Scalar1 Atmpij[1][dim*dim], Atmpji[1][dim*dim];

   SVec<Scalar2,dim> &dpdx = dpdxj.getX();
   SVec<Scalar2,dim> &dpdy = dpdxj.getY();
   SVec<Scalar2,dim> &dpdz = dpdxj.getZ();

   prod = (Scalar2) 0.0;

   int numNodes = nodes.size();
   int numEdges = edges.size();

   Vec<Vec3D>     &edgeNorm = geoState.getEdgeNormal();
   Vec<double> &edgeNormVel = geoState.getEdgeNormalVel();

   int (*edgePtr)[2] = edges.getPtr();

   bool *masterFlag = edges.getMasterFlag();

   double d_gradPhi;

   Vec3D normalDir;

   int k;

   int farfieldFluid = 0;

   double alpha_lim = 0.1;

   VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();

   if (bcMap.size() > 0)  {
     std::cout << "****Not sure about what to do *****\n";
     exit(-1);
   }

   for (l=0; l<numEdges; ++l) {

     if (!masterFlag[l]) continue;

     i = edgePtr[l][0];
     j = edgePtr[l][1];

     double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
     double length = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);

     bool intersect = LSS.edgeIntersectsStructure(0,l);

     bool iActive = LSS.isActive(0.0,i);
     bool jActive = LSS.isActive(0.0,j);    

     if( !iActive && !jActive ) {
       continue;
     }

     for (k=0; k<dim; ++k) {
       ddpij[k] = dx[0]*dpdx[i][k] + dx[1]*dpdy[i][k] + dx[2]*dpdz[i][k];
       ddpji[k] = dx[0]*dpdx[j][k] + dx[1]*dpdy[j][k] + dx[2]*dpdz[j][k];
     }

     if(iActive && jActive && !intersect) {

       recFcn->template compute<Scalar2, dim>(p[i], ddpij, p[j], ddpji, 
					      aij[l], aji[l], bij[l], bji[l], 
					      pij[0], pji[0]);

       DenseMatrixOp<Scalar1,dim,dim*dim>::applyToVector(a, numNodes + 2*l,     pji, 0, tmpi, 0);
       DenseMatrixOp<Scalar1,dim,dim*dim>::applyToVector(a, numNodes + 2*l + 1, pij, 0, tmpj, 0);

       VectorOp<Scalar2,dim>::sum(tmpi, 0, tmpj, 0, tmp, 0);
       VectorOp<Scalar2,dim>::add(tmp, 0, prod.v, i); 
       VectorOp<Scalar2,dim>::sub(tmp, 0, prod.v, j); 

     } else {

       if(iActive) {

	 LevelSetResult resij = LSS.getLevelSetDataAtEdgeCenter(0.0, l, true);
	 switch (Nriemann) {
	 case 0: //structure normal
	   d_gradPhi = dx[0]*resij.gradPhi[0]+dx[1]*resij.gradPhi[1]+dx[2]*resij.gradPhi[2];
	   normalDir = (d_gradPhi>=0.0) ? -1.0*resij.gradPhi : resij.gradPhi;
	   break;
	 case 1: //fluid normal
	   normalDir = -1.0/(edgeNorm[l].norm())*edgeNorm[l];
	   break;
	 default:
	   fprintf(stderr,"ERROR: Unknown RiemannNormal code!\n");
	   exit(-1);
	 }

	 for(k=0; k<dim; ++k) {	  
	   Vi[k]        = aij[l][k];
	   Vi[k+dim]    = Vi[k];
	   Vstar[k]     = bij[l][k];
	   Vstar[k+dim] = Vstar[k];
	 }

	 convert2(a[numNodes+2*l+1], dim, dVsdV);

	 fluxFcn[BC_INTERNAL]->computeJacobians(1.0, 0.0, edgeNorm[l], edgeNormVel[l], Vi, Vstar, dfdVi, dfdVj);

	 if (!higherOrderFSI) {

	   for(k=0; k<dim; ++k) pij[0][k] = p[i][k];

	   DenseMatrixOp<Scalar2, dim, dim*dim>::applyToDenseMatrix(&dfdVj, 0, &dVsdV, 0, &dflux1, 0);
	   for(k=0; k<dim*dim; ++k) Atmpij[0][k] = dflux1[k] + dfdVi[k];

	   DenseMatrixOp<Scalar2, dim, dim*dim>::applyToVector(Atmpij, 0, pij, 0, tmp, 0);

	 } else {

	   LevelSetResult resij = LSS.getLevelSetDataAtEdgeCenter(0.0, l, true);

	   higherOrderFSI->RcnExtrap(l, 0, i, length, resij.alpha, p, ddpij, X, fluidId, betaij[l], dVsdV, pij[0], pji[0]);

	   for(k=0; k<dim*dim; ++k) Atmpij[0][k] = dfdVi[k];
	   for(k=0; k<dim*dim; ++k) Atmpji[0][k] = dfdVj[k];

	   DenseMatrixOp<Scalar2, dim, dim*dim>::applyToVector(Atmpij, 0, pij, 0, tmpi, 0);
	   DenseMatrixOp<Scalar2, dim, dim*dim>::applyToVector(Atmpji, 0, pji, 0, tmpj, 0);
	   VectorOp<Scalar2,dim>::sum(tmpi, 0, tmpj, 0, tmp, 0);

	 }

	 VectorOp<Scalar2,dim>::add(tmp, 0, prod.v, i);

       }

       if(jActive) {

	 LevelSetResult resji = LSS.getLevelSetDataAtEdgeCenter(0.0, l, false);
	 switch (Nriemann) {
	 case 0: //structure normal
	   d_gradPhi = dx[0]*resji.gradPhi[0]+dx[1]*resji.gradPhi[1]+dx[2]*resji.gradPhi[2];
	   normalDir = (d_gradPhi>=0.0) ? resji.gradPhi : -1.0*resji.gradPhi;
	   break;
	 case 1: //fluid normal
	   normalDir = 1.0/(edgeNorm[l].norm())*edgeNorm[l];
	   break;
	 default:
	   fprintf(stderr,"ERROR: Unknown RiemannNormal code!\n");
	   exit(-1);
	 }

	 for(k=0; k<dim; ++k) {
	   Vj[k]        = bji[l][k];
	   Vj[k+dim]    = Vj[k];
	   Vstar[k]     = aji[l][k];
	   Vstar[k+dim] = Vstar[k];
	 }

	 convert2(a[numNodes+2*l], dim, dVsdV);

	 fluxFcn[BC_INTERNAL]->computeJacobians(1.0, 0.0, edgeNorm[l], edgeNormVel[l], Vstar, Vj, dfdVi, dfdVj);

	 if (!higherOrderFSI) {

	   for(k=0; k<dim; ++k) pji[0][k] = p[j][k];

	   DenseMatrixOp<Scalar2, dim, dim*dim>::applyToDenseMatrix(&dfdVi, 0, &dVsdV, 0, &dflux2, 0);
	   for(k=0; k<dim*dim; ++k) Atmpji[0][k] = dflux2[k] + dfdVj[k];

	   DenseMatrixOp<Scalar2, dim, dim*dim>::applyToVector(Atmpji, 0, pji, 0, tmp, 0);

	 }else{

	   LevelSetResult resji = LSS.getLevelSetDataAtEdgeCenter(0.0, l, false);

	   higherOrderFSI->RcnExtrap(l, 1, j, length, resji.alpha, p, ddpji, X, fluidId, betaji[l], dVsdV, pji[0], pij[0]);

	   for(k=0; k<dim*dim; ++k) Atmpij[0][k] = dfdVi[k];
	   for(k=0; k<dim*dim; ++k) Atmpji[0][k] = dfdVj[k];

	   DenseMatrixOp<Scalar2, dim, dim*dim>::applyToVector(Atmpij, 0, pij, 0, tmpi, 0);
	   DenseMatrixOp<Scalar2, dim, dim*dim>::applyToVector(Atmpji, 0, pji, 0, tmpj, 0);
	   VectorOp<Scalar2,dim>::sum(tmpi, 0, tmpj, 0, tmp, 0);	  

	 }

	 VectorOp<Scalar2,dim>::sub(tmp, 0, prod.v, j);
 
      }

    }

  }

  // contribution from diagonal entries of A
  for (i=0; i<numNodes; ++i) {

    DenseMatrixOp<Scalar2,dim,dim*dim>::applyAndAddToVector(a, i, p.v, i, prod.v, i);
    
    double voli = 1.0 / ctrlVol[i];
    for (int k=0; k<dim; ++k) prod[i][k] *= voli;

    if(!LSS.isActive(0.0,i)) for (int k=0; k<dim; ++k) prod[i][k] = 0.0;

  }

}
//------------------------------------------------------------------------------

template<class Scalar1, class Scalar2, int dim>
void SubDomain::computeMatVecProdH2T(RecFcn *recFcn, SVec<double,3> &X,
                Vec<double> &ctrlVol, GenMat<Scalar1,dim> &A,
                SVec<double,dim> &aij, SVec<double,dim> &aji,
                SVec<double,dim> &bij, SVec<double,dim> &bji,
                SVec<Scalar2,dim> &p, SVec<Scalar2,dim> &zu,
                SVec<Scalar2,dim> &zgx, SVec<Scalar2,dim> &zgy,
                SVec<Scalar2,dim> &zgz) {
  int i, j, l, k;

  // Same notations as in the Fortran code
  Scalar2 zu_is1[1][dim], zu_is2[1][dim], zgx_is1[1][dim], zgx_is2[1][dim];
  Scalar2 zgy_is1[1][dim], zgy_is2[1][dim], zgz_is1[1][dim], zgz_is2[1][dim];
  //Scalar2 tmp[1][dim];
  Scalar2 tmpi[1][dim], tmpj[1][dim];
  Scalar2 tmp1[1][dim];
  Scalar1 (*a)[dim*dim] = A.data();


  zu = (Scalar2)0.0;
  zgx = (Scalar2)0.0;
  zgy = (Scalar2)0.0;
  zgz = (Scalar2)0.0;

  int numNodes = nodes.size();
  int numEdges = edges.size();

  int (*edgePtr)[2] = edges.getPtr();

  bool *masterFlag = edges.getMasterFlag();
  for (l=0; l<numEdges; ++l) {

    if (!masterFlag[l]) continue;
    i = edgePtr[l][0];
    j = edgePtr[l][1];

    double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};

    for (k=0; k<dim; ++k)
      tmp1[0][k] = p[i][k] - p[j][k];

    DenseMatrixOp<Scalar1,dim, dim*dim>::applyTransToVector(a, numNodes + 2*l, tmp1, 0, tmpi, 0);
    DenseMatrixOp<Scalar1,dim, dim*dim>::applyTransToVector(a, numNodes + 2*l + 1, tmp1, 0, tmpj, 0);
    // These routines do not exist anymore
    //denseMatrixTransTimesVector(a, numNodes + 2*l, tmp1, 0, tmpi, 0);
    //denseMatrixTransTimesVector(a, numNodes + 2*l + 1, tmp1, 0, tmpj, 0);

    recFcn->template computeT<Scalar2, dim> (dx, tmpi[0], tmpj[0], aij[l], aji[l], bij[l], bji[l], i,
                j, zu_is1[0], zu_is2[0], zgx_is1[0], zgx_is2[0], zgy_is1[0], zgy_is2[0], zgz_is1[0], zgz_is2[0]);


    VectorOp<Scalar2,dim>::add(zu_is1, 0, zu.v, i);
    VectorOp<Scalar2,dim>::add(zu_is2, 0, zu.v, j);
    VectorOp<Scalar2,dim>::add(zgx_is1, 0, zgx.v, i);
    VectorOp<Scalar2,dim>::add(zgx_is2, 0, zgx.v, j);
    VectorOp<Scalar2,dim>::add(zgy_is1, 0, zgy.v, i);
    VectorOp<Scalar2,dim>::add(zgy_is2, 0, zgy.v, j);
    VectorOp<Scalar2,dim>::add(zgz_is1, 0, zgz.v, i);
    VectorOp<Scalar2,dim>::add(zgz_is2, 0, zgz.v, j);


/*
    // These routines do not exist anymore
    addVector(p1, 0, prod2.v, i);
    addVector(p2, 0, prod2.v, j);
    addVector(p3, 0, prod.v, i);
    addVector(p4, 0, prod.v, j);
    addVector(p5, 0, prod3.v, i);
    addVector(p6, 0, prod3.v, j);
    addVector(p7, 0, prod4.v, i);
    addVector(p8, 0, prod4.v, j);
*/
  }
  for (i=0; i<numNodes; ++i)
    DenseMatrixOp<Scalar1,dim, dim*dim>::applyTransAndAddToVector(a, i, p.v, i, zu.v, i);

    // This routine does not exist anymore
    //addDenseMatrixTransTimesVector(a, i, p.v, i, prod2.v, i);

}

//------------------------------------------------------------------------------

template<class Scalar1, class Scalar2, int dim>
void SubDomain::computeMatVecProdH2Tb(RecFcn *recFcn, SVec<double,3> &X,
                Vec<double> &ctrlVol, GenMat<Scalar1,dim> &A,
                NodalGrad<dim, Scalar2> &dpdxj, SVec<Scalar2,dim> &p,
                SVec<Scalar2,dim> &prod, SVec<Scalar2,dim> &zu)  {

  int i, l;
  Scalar2 p1[1][dim];

  SVec<Scalar2,dim> &dpdx = dpdxj.getX();
  SVec<Scalar2,dim> &dpdy = dpdxj.getY();
  SVec<Scalar2,dim> &dpdz = dpdxj.getZ();

  prod = (Scalar2) 0.0;

  int numNodes = nodes.size();
  int numEdges = edges.size();

  int (*edgePtr)[2] = edges.getPtr();

  for (l=0; l<numNodes; ++l) {

    i = l;
    recFcn->computeTb(zu, dpdx, dpdy, dpdz, i, p1[0]);
    VectorOp<Scalar2,dim>::add(p1, 0, prod.v, i);
    //addVector(p1, 0, prod.v, i);
  }
}

//------------------------------------------------------------------------------

template<class Scalar>
void SubDomain::setComLenNodes(int dim, CommPattern<Scalar> &cp)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub)
    cp.setLen(sndChannel[iSub], sharedNodes->num(iSub)*dim);

}

//------------------------------------------------------------------------------

template<class Scalar>
void SubDomain::setComLenInletNodes(int dim, CommPattern<Scalar> &cp)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub)
    cp.setLen(sndChannel[iSub], sharedInletNodes->num(iSub)*dim);

}

//------------------------------------------------------------------------------

template<class Scalar>
void SubDomain::setComLenEdges(int dim, CommPattern<Scalar> &cp)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub)
    cp.setLen(sndChannel[iSub], numSharedEdges[iSub]*dim);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::sndData(CommPattern<Scalar> &sp, Scalar (*w)[dim])
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode) {
      for (int j = 0; j < dim; ++j)
	buffer[iNode][j] = w[ (*sharedNodes)[iSub][iNode] ][j];
    }

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::addRcvData(CommPattern<Scalar> &sp, Scalar (*w)[dim])
{
  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode)
      for (int j = 0; j < dim; ++j)  {	// KTC: COULD DO ONLY FOR ONE LAYER OF NODES AWAY
	w[ (*sharedNodes)[iSub][iNode] ][j] += buffer[iNode][j];
      }
  }
}
//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::RcvData(CommPattern<Scalar> &sp, Scalar (*w)[dim])
{
  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode)
      for (int j = 0; j < dim; ++j)  {	// KTC: COULD DO ONLY FOR ONE LAYER OF NODES AWAY
			if(w[ (*sharedNodes)[iSub][iNode] ][j] == 0.0)
				w[ (*sharedNodes)[iSub][iNode] ][j] = buffer[iNode][j];
      }
  }
}
//------------------------------------------------------------------------------

template<int dim>
void SubDomain::sndGhostStates(CommPattern<double> &sp, Vec<GhostPoint<dim>*> &ghostPoints, int shift)
{

  double *v;

	for(int iSub = 0; iSub < numNeighb; ++iSub) 
	{
    SubRecInfo<double> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    double (*buffer)[dim] = reinterpret_cast<double (*)[dim]>(sInfo.data);

		for(int iNode=0; iNode<sharedNodes->num(iSub); ++iNode)
		{
      if(ghostPoints[ (*sharedNodes)[iSub][iNode] ])
	{
	  v = ghostPoints[ (*sharedNodes)[iSub][iNode] ]->getState();
				
				for (int j=0; j<dim; ++j) buffer[iNode][j] = v[j+shift];
	}
      else 
	{
				// It can happened that all the edges containing the 
				// ghostPoint are not master in this SubDomain.
	  // In which case we do not want to send garbage.
	  for (int j = 0; j < dim; ++j) buffer[iNode][j] = 0.0;
	}
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::sndNumGhostStates(CommPattern<int> &sp, Vec<GhostPoint<dim>*> &ghostPoints)
{

	for(int iSub = 0; iSub < numNeighb; ++iSub) 
	{
    SubRecInfo<int> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    int (*buffer)[1] = reinterpret_cast<int (*)[1]>(sInfo.data);

		for(int iNode=0; iNode<sharedNodes->num(iSub); ++iNode) 
		{
      if(ghostPoints[ (*sharedNodes)[iSub][iNode] ])
	  buffer[iNode][0] = ghostPoints[ (*sharedNodes)[iSub][iNode] ]->ng;
      else 
	  buffer[iNode][0] = 0;
	}
    }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::sndGhostWeights(CommPattern<double> &sp, Vec<GhostPoint<dim>*> &ghostPoints, int shift)
{

  double *v;

	for(int iSub=0; iSub<numNeighb; ++iSub) 
	{
    SubRecInfo<double> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    double (*buffer)[dim] = reinterpret_cast<double (*)[dim]>(sInfo.data);

		for(int iNode=0; iNode<sharedNodes->num(iSub); ++iNode) 
		{			
      if(ghostPoints[ (*sharedNodes)[iSub][iNode] ])
	{
	  v = ghostPoints[ (*sharedNodes)[iSub][iNode] ]->Ws;

				for(int j=0; j<dim; ++j) buffer[iNode][j] = v[j+shift];
	}
      else 
	{
				// It can happened that all the edges containing 
            // the ghostPoint are not master in this SubDomain.
	  // In which case we do not want to send garbage.
	  for (int j = 0; j < dim; ++j) buffer[iNode][j] = 0.0;
	}
    }
  }
}
//------------------------------------------------------------------------------

template<int dim>
void SubDomain::sndGhostTags(CommPattern<int> &sp, Vec<GhostPoint<dim>*> &ghostPoints)
{

	for(int iSub = 0; iSub < numNeighb; ++iSub) 
	{
    SubRecInfo<int> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    int (*buffer)[1] = reinterpret_cast<int (*)[1]>(sInfo.data);

		for(int iNode=0; iNode<sharedNodes->num(iSub); ++iNode) 
		{
      if(ghostPoints[ (*sharedNodes)[iSub][iNode] ])
	  buffer[iNode][0] = ghostPoints[ (*sharedNodes)[iSub][iNode] ]->ghostTag;
      else 
	  buffer[iNode][0] = -2;
	}
    }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::rcvGhostStates(CommPattern<double> &sp, Vec<GhostPoint<dim>*> &ghostPoints, int shift)
{

  GhostPoint<dim> *gp;

	for(int iSub = 0; iSub < numNeighb; ++iSub) 
	{	
    SubRecInfo<double> sInfo = sp.recData(rcvChannel[iSub]);
    double (*buffer)[dim] = reinterpret_cast<double (*)[dim]>(sInfo.data);

		for (int iNode=0; iNode<sharedNodes->num(iSub); ++iNode) 
		{
      if(ghostPoints[ (*sharedNodes)[iSub][iNode] ])
	{
	  gp = ghostPoints[ (*sharedNodes)[iSub][iNode] ];

				for(int j=0; j<dim; ++j) gp->Vg[j+shift] += buffer[iNode][j];
			}				
		}
	}

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::rcvNumGhostStates(CommPattern<int> &sp, Vec<GhostPoint<dim>*> &ghostPoints, VarFcn *varFcn)
{

	for(int iSub = 0; iSub < numNeighb; ++iSub) 
	{
    SubRecInfo<int> sInfo = sp.recData(rcvChannel[iSub]);
    int (*buffer)[1] = reinterpret_cast<int (*)[1]>(sInfo.data);

		for(int iNode=0; iNode < sharedNodes->num(iSub); ++iNode) 
		{
			// if the weight is zero, the whole buffered state is gonna be zero.
			if( buffer[iNode][0] == 0) continue;

      if(!ghostPoints[ (*sharedNodes)[iSub][iNode] ])
	  ghostPoints[ (*sharedNodes)[iSub][iNode] ] = new GhostPoint<dim>(varFcn);

      ghostPoints[ (*sharedNodes)[iSub][iNode] ]->ng += buffer[iNode][0];    
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::rcvGhostWeights(CommPattern<double> &sp, Vec<GhostPoint<dim>*> &ghostPoints, int shift)
{

  GhostPoint<dim> *gp;

	for(int iSub = 0; iSub < numNeighb; ++iSub) 
	{
    SubRecInfo<double> sInfo = sp.recData(rcvChannel[iSub]);
    double (*buffer)[dim] = reinterpret_cast<double (*)[dim]>(sInfo.data);

		for(int iNode=0; iNode<sharedNodes->num(iSub); ++iNode) 
		{
      if(ghostPoints[ (*sharedNodes)[iSub][iNode] ])
	{
	  gp = ghostPoints[ (*sharedNodes)[iSub][iNode] ];

				for (int j=0; j<dim; ++j) gp->Ws[j+shift] += buffer[iNode][j];
			}
		}
	}

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::rcvGhostTags(CommPattern<int> &sp, Vec<GhostPoint<dim>*> &ghostPoints)
{
  GhostPoint<dim> *gp;

	for(int iSub = 0; iSub < numNeighb; ++iSub) 
	{
    SubRecInfo<int> sInfo = sp.recData(rcvChannel[iSub]);
    int (*buffer)[1] = reinterpret_cast<int (*)[1]>(sInfo.data);

		for(int iNode=0; iNode<sharedNodes->num(iSub); ++iNode) 
		{			
      if( buffer[iNode][0] < 0) continue; // if the ghostTag is not set 
      if(ghostPoints[ (*sharedNodes)[iSub][iNode] ])
      {
        gp = ghostPoints[ (*sharedNodes)[iSub][iNode] ];
        if (gp->ghostTag<0) gp->ghostTag = buffer[iNode][0];    
        else if (gp->ghostTag != buffer[iNode][0])
        {
          fprintf(stderr,"The two ghost States refer to different Fluids in comm\n");
          fprintf(stderr,"ghostTag: %i, buffer: %i\n",gp->ghostTag,buffer[iNode][0]);
          exit(-1);
        }
      }  
    }
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class OpType>
void SubDomain::operateRcvData(CommPattern<Scalar> &sp, Scalar (*w)[dim], const OpType &oper)
{
  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode)
      for (int j = 0; j < dim; ++j)  {
        w[ (*sharedNodes)[iSub][iNode] ][j] = OpType::apply(w[ (*sharedNodes)[iSub][iNode]][j], buffer[iNode][j]);
      }
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::sndInletData(CommPattern<Scalar> &sp, Scalar (*w)[dim])
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<Scalar> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);
    for (int iNode = 0; iNode < sharedInletNodes->num(iSub); ++iNode)
      for (int j = 0; j < dim; ++j)
        buffer[iNode][j] = w[ (*sharedInletNodes)[iSub][iNode] ][j];
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::addRcvInletData(CommPattern<Scalar> &sp, Scalar (*w)[dim], bool ForExtrapolation)
{
  if(!ForExtrapolation){
    /* the values in w are accessed via the inlet node number
     * in the inletNode number chart.
     * This communication of values is used for
     * computation of normals at inlet nodes
     * and computation of number of faces surrounding
     * a given inlet node.
     */

    for (int iSub = 0; iSub < numNeighb; ++iSub) {
      SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
      Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);

      for (int iNode = 0; iNode < sharedInletNodes->num(iSub); ++iNode)
        for (int j = 0; j < dim; ++j)
          w[ (*sharedInletNodes)[iSub][iNode] ][j] += buffer[iNode][j];

    }
  }else{
    /* the values in w are accessed via the inlet node number
     * in the inletNode number chart.
     * This communication of values is used for
     * computation of extrapolation values
     * in the whole domain.
     */

    for (int iSub = 0; iSub < numNeighb; ++iSub) {
      SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
      Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);
      int inletnode;
      for (int iNode = 0; iNode < sharedInletNodes->num(iSub); ++iNode){
        inletnode = (*sharedInletNodes)[iSub][iNode] ;
        if (w[inletnode][0] == 0.0){
          for (int j = 0; j < dim; ++j)
            w[inletnode][j] += buffer[iNode][j];
        }
        else{
          //if(buffer[iNode][0]!=0.0)
          //      fprintf(stderr, "two shared inlet nodes have their own extrapolation and their densities are %f vs %f\n", w[ (*sharedInletNodes)[iSub][iNode] ][0], buffer[iNode][0]);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::sndInletRhsData(CommPattern<Scalar> &sp, Scalar (*w)[dim])
{
  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<Scalar> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);
    int inletnode, node;

    for (int iNode = 0; iNode < sharedInletNodes->num(iSub); ++iNode) {
      inletnode = (*sharedInletNodes)[iSub][iNode] ;
      node = inletNodes[inletnode].getNodeNum();
      //fprintf(stdout, "node in two domains %d for subd %d\n", node, locSubNum);
      for (int j = 0; j < dim; ++j)
        buffer[iNode][j] = w[ node ][j];
    }
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::addRcvInletRhsData(CommPattern<Scalar> &sp, Scalar (*w)[dim])
{
    /* the values in w are accessed via the inlet node number
     * in the subDomain number chart.
     * This communication of values is used for
     * the recomputation of the rhs term, which
     * contains the flux (computeFVT = computeResidual)
     */

        // we want to pass the data of a shared inlet node from one subdomain to a
        // neighbouring subdomain. Only one of the data has a physical value, while
        // the other ones are set to 0.0, except in one case where the intersection
        // between the normal of the inlet node and the opposite face of the tetrahedron
        // gives one of the edges of the tetrahedra.
  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);
    int inletnode, node;
    for (int iNode = 0; iNode < sharedInletNodes->num(iSub); ++iNode){
      inletnode = (*sharedInletNodes)[iSub][iNode] ;
      node = inletNodes[inletnode].getNodeNum();
      if (w[node][0] == 0.0){
        for (int j = 0; j < dim; ++j)
          w[node][j] += buffer[iNode][j];
      }
      else{
        //if(buffer[iNode][0]!=0.0){
          // fprintf(stderr, "two shared inlet nodes (inlet node %d and node %d) have their own extrapolation\n     and their densitiesRHS are %.14e vs %.14e\n", inletnode, locToGlobNodeMap[node]+1, w[ node ][0], buffer[iNode][0]);
        //}
      }
    }
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::minRcvData(CommPattern<Scalar> &sp, Scalar (*w)[dim])
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode)
      for (int j = 0; j < dim; ++j)
        if (buffer[iNode][j] < w[ (*sharedNodes)[iSub][iNode] ][j])
          w[ (*sharedNodes)[iSub][iNode] ][j] = buffer[iNode][j];

  }

}

//------------------------------------------------------------------------------
template<class Scalar, int dim>
void SubDomain::maxRcvData(CommPattern<Scalar> &sp, Scalar (*w)[dim])
{
  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);
    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode)
      for (int j = 0; j < dim; ++j)
        if (buffer[iNode][j] > w[ (*sharedNodes)[iSub][iNode] ][j])
          w[ (*sharedNodes)[iSub][iNode] ][j] = buffer[iNode][j];
  }
}

//------------------------------------------------------------------------------
// Adam 2011.09: Called from Domain::pseudoFastMarchingMethod. 
// Need to take into account the updated value after communication
template<class Scalar, int dim>
void SubDomain::maxRcvDataAndCountUpdates(CommPattern<Scalar> &sp, Scalar (*w)[dim],int &nSortedNodes, Vec<int> &sortedNodes)
{
  assert(dim == 1); // if you intend to use it in Vectorial mode, modify it your way

  int sharedNodeID;
  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);
    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode) {
      sharedNodeID = (*sharedNodes)[iSub][iNode];
      for (int j = 0; j < dim; ++j) {
        if (buffer[iNode][j] > w[sharedNodeID][j]) {
          w[sharedNodeID][j]        = buffer[iNode][j];
	  sortedNodes[nSortedNodes] = sharedNodeID;
	  nSortedNodes++;
	}
      }
    }
  }
}

//------------------------------------------------------------------------------

template<class Scalar1, class Scalar2, int dim1, int dim2>
void SubDomain::TagPsiExchangeData(CommPattern<Scalar1> &splevel, Scalar1 (*level)[dim1],
                                   CommPattern<Scalar2> &sppsi, Scalar2 (*psi)[dim2])
{
  /* it is assumed that dim1 = 1, dim2 = 1 */

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar1> sInfolevel = splevel.recData(rcvChannel[iSub]);
    Scalar1 (*blevel)[dim1] = reinterpret_cast<Scalar1 (*)[dim1]>(sInfolevel.data);
    SubRecInfo<Scalar2> sInfopsi = sppsi.recData(rcvChannel[iSub]);
    Scalar2 (*bpsi)[dim2] = reinterpret_cast<Scalar2 (*)[dim2]>(sInfopsi.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode)
      if (level[ (*sharedNodes)[iSub][iNode] ][0] == 0 &&
          blevel[iNode][0] > 0){
        level[ (*sharedNodes)[iSub][iNode] ][0] = blevel[iNode][0];
        psi[ (*sharedNodes)[iSub][iNode] ][0]   = bpsi[iNode][0];
      }

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::sndDiagBlocks(CommPattern<Scalar> &sp, GenMat<Scalar,dim> &A)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    Scalar (*buffer)[dim*dim] = reinterpret_cast<Scalar (*)[dim*dim]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode) {

      Scalar *a = A.getElem_ii((*sharedNodes)[iSub][iNode]);

      for (int j=0; j<dim*dim; ++j) buffer[iNode][j] = a[j];

    }

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::addRcvDiagBlocks(CommPattern<Scalar> &sp, GenMat<Scalar,dim> &A)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    Scalar (*buffer)[dim*dim] = reinterpret_cast<Scalar (*)[dim*dim]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode) {

      Scalar *a = A.getElem_ii((*sharedNodes)[iSub][iNode]);

      for (int j=0; j<dim*dim; ++j) a[j] += buffer[iNode][j];

    }

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::sndDiagInletBlocks(CommPattern<Scalar> &sp, GenMat<Scalar,dim> &A)
{
  int node;
  int type;
  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    Scalar (*buffer)[dim*dim] = reinterpret_cast<Scalar (*)[dim*dim]>(sInfo.data);

    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode) {
      node = (*sharedNodes)[iSub][iNode];
      type = nodeType[node];
      if (!(type == BC_INLET_MOVING || type == BC_OUTLET_MOVING ||
            type == BC_INLET_FIXED  || type == BC_OUTLET_FIXED) ){

        Scalar *a = A.getElem_ii((*sharedNodes)[iSub][iNode]);
        for (int j=0; j<dim*dim; ++j) buffer[iNode][j] = a[j];
      }
    }
  }


}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::addRcvDiagInletBlocks(CommPattern<Scalar> &sp, GenMat<Scalar,dim> &A)
{
  int node;
  int type;

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    Scalar (*buffer)[dim*dim] = reinterpret_cast<Scalar (*)[dim*dim]>(sInfo.data);
    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode) {
      node = (*sharedNodes)[iSub][iNode];
      type = nodeType[node];
      if (!(type == BC_INLET_MOVING || type == BC_OUTLET_MOVING ||
            type == BC_INLET_FIXED  || type == BC_OUTLET_FIXED) ){

        Scalar *a = A.getElem_ii((*sharedNodes)[iSub][iNode]);
        for (int j=0; j<dim*dim; ++j) a[j] += buffer[iNode][j];
      }else{
        Scalar *a = A.getElem_ii((*sharedNodes)[iSub][iNode]);
        if (a[0] == 0.0)
          for (int k=0; k<dim*dim; k++) a[k] += buffer[iNode][k];
      }

    }

  }

}
//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::sndEdgeData(CommPattern<Scalar> &sp, Scalar (*w)[dim])
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);

    for (int iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge)
      for (int k=0; k<dim; ++k)
	buffer[iEdge][k] = w[ sharedEdges[iSub][iEdge].edgeNum ][k];

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::addRcvEdgeData(CommPattern<Scalar> &sp, Scalar (*w)[dim])
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.getSendBuffer(rcvChannel[iSub]);
    Scalar (*buffer)[dim] = reinterpret_cast<Scalar (*)[dim]>(sInfo.data);

    for (int iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge)
      for (int k=0; k<dim; ++k)
	w[ sharedEdges[iSub][iEdge].edgeNum ][k] += buffer[iEdge][k];

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::sndOffDiagBlocks(CommPattern<Scalar> &sp, GenMat<Scalar,dim> &A)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    Scalar (*buffer)[2][dim*dim] = reinterpret_cast<Scalar (*)[2][dim*dim]>(sInfo.data);

    for (int iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge) {

      int edgeNum = sharedEdges[iSub][iEdge].edgeNum;

      Scalar *aij, *aji;

      if (sharedEdges[iSub][iEdge].sign > 0) {
	aij = A.getElem_ij(edgeNum);
	aji = A.getElem_ji(edgeNum);
      }
      else {
	aij = A.getElem_ji(edgeNum);
	aji = A.getElem_ij(edgeNum);
      }

      if (aij && aji) {
	for (int k=0; k<dim*dim; ++k) {
	  buffer[iEdge][0][k] = aij[k];
	  buffer[iEdge][1][k] = aji[k];
	}
      }
      else {
	for (int k=0; k<dim*dim; ++k) {
	  buffer[iEdge][0][k] = 0.0;
	  buffer[iEdge][1][k] = 0.0;
	}
      }

    }

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::addRcvOffDiagBlocks(CommPattern<Scalar> &sp, GenMat<Scalar,dim> &A)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    Scalar (*buffer)[2][dim*dim] = reinterpret_cast<Scalar (*)[2][dim*dim]>(sInfo.data);

    for (int iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge) {

      int edgeNum = sharedEdges[iSub][iEdge].edgeNum;

      Scalar *aij, *aji;

      if (sharedEdges[iSub][iEdge].sign > 0) {
	aij = A.getElem_ij(edgeNum);
	aji = A.getElem_ji(edgeNum);
      }
      else {
	aij = A.getElem_ji(edgeNum);
	aji = A.getElem_ij(edgeNum);
      }

      if (aij && aji) {
	for (int k=0; k<dim*dim; ++k) {
	  aij[k] += buffer[iEdge][0][k];
	  aji[k] += buffer[iEdge][1][k];
	}
      }

    }

  }

}

template<class Scalar, int dim>
void SubDomain::sndGhostOffDiagBlocks(CommPattern<Scalar> &sp, GenMat<Scalar,dim> &A)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    Scalar (*buffer)[2][dim*dim] = reinterpret_cast<Scalar (*)[2][dim*dim]>(sInfo.data);
    int (*ptr)[2] = edges.getPtr();

    for (int iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge) {

      int edgeNum = sharedEdges[iSub][iEdge].edgeNum;
      int i = ptr[edgeNum][0],j = ptr[edgeNum][1];

      Scalar *aij, *aji;

      if (sharedEdges[iSub][iEdge].sign > 0) {
	aij = A.queryRealNodeElem_ij(i,j);
	if (!aij)
	  aij = A.queryGhostNodeElem_ij(i,j);
	aji = A.queryRealNodeElem_ij(j,i);
	if (!aji)
	  aji = A.queryGhostNodeElem_ij(j,i);
      }
      else {
	aji = A.queryRealNodeElem_ij(i,j);
	if (!aji)
	  aji = A.queryGhostNodeElem_ij(i,j);
	aij = A.queryRealNodeElem_ij(j,i);
	if (!aij)
	  aij = A.queryGhostNodeElem_ij(j,i);
      }

      if (aij && aji) {
	for (int k=0; k<dim*dim; ++k) {
	  buffer[iEdge][0][k] = aij[k];
	  buffer[iEdge][1][k] = aji[k];
	}
      }
      else {
	for (int k=0; k<dim*dim; ++k) {
	  buffer[iEdge][0][k] = 0.0;
	  buffer[iEdge][1][k] = 0.0;
	}
      }

    }

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::addRcvGhostOffDiagBlocks(CommPattern<Scalar> &sp, GenMat<Scalar,dim> &A)
{

  for (int iSub = 0; iSub < numNeighb; ++iSub) {

    SubRecInfo<Scalar> sInfo = sp.recData(rcvChannel[iSub]);
    int (*ptr)[2] = edges.getPtr();
    Scalar (*buffer)[2][dim*dim] = reinterpret_cast<Scalar (*)[2][dim*dim]>(sInfo.data);

    for (int iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge) {
      
      int edgeNum = sharedEdges[iSub][iEdge].edgeNum;
      int i = ptr[edgeNum][0],j = ptr[edgeNum][1];

      Scalar *aij, *aji;

      if (sharedEdges[iSub][iEdge].sign > 0) {
	aij = A.queryRealNodeElem_ij(i,j);
	if (!aij)
	  aij = A.queryGhostNodeElem_ij(i,j);
	aji = A.queryRealNodeElem_ij(j,i);
	if (!aji)
	  aji = A.queryGhostNodeElem_ij(j,i);
      }
      else {
	aji = A.queryRealNodeElem_ij(i,j);
	if (!aji)
	  aji = A.queryGhostNodeElem_ij(i,j);
	aij = A.queryRealNodeElem_ij(j,i);
	if (!aij)
	  aij = A.queryGhostNodeElem_ij(j,i);
      }

      if (aij && aji) {
	for (int k=0; k<dim*dim; ++k) {
	  aij[k] += buffer[iEdge][0][k];
	  aji[k] += buffer[iEdge][1][k];
	}
      }
    }

  }

}

//------------------------------------------------------------------------------
template<class Scalar, int dim>
bool SubDomain::checkIfFileExists(const char *prefix)
{

  char name[MAXLINE];
  sprintf(name, "%s%s", prefix, suffix);

  struct stat buf;
  
  return (stat(name, &buf) != -1) ? true : false;

}


//------------------------------------------------------------------------------

template<class Scalar, int dim>
double SubDomain::readTagFromFile(const char *prefix, int no, int *neq, int *nsol)
{

  char name[MAXLINE];
  sprintf(name, "%s%s", prefix, suffix);

  BinFileHandler file(name, "rb");

  int info[3];
  file.read(info, 3);

  if (info[0] != numClusNodes) {
    fprintf(stderr, "*** Error: mismatch in size for \'%s\' (%d vs %d)\n",
	    name, info[0], numClusNodes);
    exit(1);
  }

  *neq = info[1];
  *nsol = info[2];
//	fprintf(stderr, "numStep is %d.\n", info[2]);
  double tag = 0.0;

  if (no < *nsol) {
    file.seek(3*sizeof(int) + no*(sizeof(double) + numClusNodes*(*neq)*sizeof(Scalar)));
    file.read(&tag, 1);
  }

  return tag;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::openFileForWriting(const char *prefix, int no)
{

  if (clusSubNum != 0) return;

  char name[MAXLINE];
  sprintf(name, "%s%s", prefix, suffix);

#ifdef SYNCHRO_WRITE
  const char *flag = "ws+";
  if (no == 0) flag = "ws";
#else
  const char *flag = "w+";
  if (no == 0) flag = "w";
#endif

  BinFileHandler file(name, flag);

  int info[3] = {numClusNodes, dim, no + 1};
  file.write(info, 3);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::writeTagToFile(const char *prefix, int no, double tag)
{

  if (clusSubNum != 0) return;

  BinFileHandler::OffType unit = dim * sizeof(Scalar);

  char name[MAXLINE];
  sprintf(name, "%s%s", prefix, suffix);

#ifdef SYNCHRO_WRITE
  BinFileHandler file(name, "ws+");
#else
  BinFileHandler file(name, "w+");
#endif

  file.seek(3*sizeof(int) + no*(sizeof(double) + numClusNodes*unit));
  file.write(&tag, 1);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::readVectorFromFile(const char *prefix, int no, int neq,
				   SVec<Scalar,dim> &U, Scalar* scale)
{

  char name[MAXLINE];
  sprintf(name, "%s%s", prefix, suffix);

  BinFileHandler file(name, "rb");

  BinFileHandler::OffType unit = neq * sizeof(Scalar);

  BinFileHandler::OffType pos = 3*sizeof(int) +
    no*(sizeof(double) + numClusNodes*unit) + sizeof(double);

  Scalar *data;

  if (neq == dim)
    data = reinterpret_cast<Scalar *>(U.data());
  else
    data = new Scalar[U.size()*neq];

  int i, count = 0;
  for (i=0; i<numNodeRanges; ++i) {
    file.seek(pos + nodeRanges[i][1]*unit);
    file.read(data + count*neq, nodeRanges[i][0]*neq);
    count += nodeRanges[i][0];
  }

  if (scale) {
    for (i=0; i<U.size()*neq; ++i)
      data[i] *= *scale;
  }

  if (neq != dim) {
    int minsize = min(neq, dim);
    for (i=0; i<U.size(); ++i) {
      for (int j=0; j<minsize; ++j)
        U[i][j] = data[neq*i + j];
    }
    delete [] data;
  }

}

//------------------------------------------------------------------------------

template<class Scalar>
void SubDomain::readVectorFromFile(const char *prefix, int no, Vec<Scalar> &U) {
    char name[MAXLINE];
    sprintf(name, "%s%s", prefix, suffix);

    BinFileHandler file(name, "rb");

    BinFileHandler::OffType unit = sizeof(Scalar);

    BinFileHandler::OffType pos = 3 * sizeof(int) + sizeof(double) +
                                  no * (sizeof(double) + numClusNodes * unit);

    Scalar *data = new Scalar[U.size()];

    int i, count = 0;
    for (i = 0; i < numNodeRanges; ++i) {
        file.seek(pos + nodeRanges[i][1] * unit);
        file.read(data + count, nodeRanges[i][0]);
        count += nodeRanges[i][0];
    }

    for (i = 0; i < U.size(); ++i) {
        U[i] = data[i];
    }

    delete[] data;
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void SubDomain::writeVectorToFile(const char *prefix, int no,
				  SVec<Scalar,dim> &U, Scalar* scale)
{
  char name[MAXLINE];
  sprintf(name, "%s%s", prefix, suffix);
#ifdef SYNCHRO_WRITE
  BinFileHandler file(name, "ws+");
#else
  BinFileHandler file(name, "w+");
#endif

  BinFileHandler::OffType unit = dim * sizeof(Scalar);
  BinFileHandler::OffType pos = 3*sizeof(int) +
    no*(sizeof(double) + numClusNodes*unit) + sizeof(double);

  Scalar (*data)[dim];
  if (scale) {
    data = new Scalar[U.size()][dim];
    Scalar* v = reinterpret_cast<Scalar*>(data);
    Scalar* u = reinterpret_cast<Scalar*>(U.data());
    for (int i=0; i<U.size()*dim; ++i)
      v[i] = (*scale) * u[i];
  } else {
    data = U.data();
  }

  int count = 0;
  for (int i=0; i<numNodeRanges; ++i) {
    if (nodeRanges[i][2]) {
      file.seek(pos + nodeRanges[i][1]*unit);
      file.write(reinterpret_cast<Scalar *>(data + count), nodeRanges[i][0]*dim);
    }
    count += nodeRanges[i][0];
  }

  if (scale)
    delete [] data;

}

//------------------------------------------------------------------------------

template<class Scalar>
void SubDomain::writeVectorToFile(const char *prefix, int no,
                                  Vec<Scalar> &U, Scalar* scale)
{
    char name[MAXLINE];
    sprintf(name, "%s%s", prefix, suffix);
#ifdef SYNCHRO_WRITE
    BinFileHandler file(name, "ws+");
#else
    BinFileHandler file(name, "w+");
#endif

    BinFileHandler::OffType unit =  sizeof(Scalar);
    BinFileHandler::OffType pos = 3*sizeof(int) +
                                  no*(sizeof(double) + numClusNodes*unit) + sizeof(double);

    Scalar *data;
    if (scale) {
        data = new Scalar[U.size()];
        Scalar* v = data;
        Scalar* u = U.data();
        for (int i=0; i<U.size(); ++i)
            v[i] = (*scale) * u[i];
    } else {
        data = U.data();
    }

    int count = 0;
    for (int i=0; i<numNodeRanges; ++i) {
        if (nodeRanges[i][2]) {
            file.seek(pos + nodeRanges[i][1]*unit);
            file.write(data + count, nodeRanges[i][0]);
        }
        count += nodeRanges[i][0];
    }

    if (scale)
        delete [] data;

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::assignFreeStreamValues2(SVec<double,dim> &Uin, SVec<double,dim> &Uout,
					SVec<double,dim> &U, SVec<double,dim> &Uinlet)
{
  int node;
  for (int j=0; j<inletNodes.size(); j++){
    node = inletNodes[j].getNodeNum();
    inletNodes[j].template assignFreeStreamValues<dim>(nodeType[node],
				      Uin[node], Uout[node], Uinlet[j]);
  }

  for (int i=0; i<faces.size(); ++i) 
    faces[i].template assignFreeStreamValues2<dim>(Uin, Uout, U[i]);
  
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::assignFreeStreamValues(double *Uin, double *Uout, SVec<double,dim> &U, SVec<double,dim> &Uinlet)
{

  for (int j=0; j<inletNodes.size(); j++)
    inletNodes[j].template assignFreeStreamValues<dim>(nodeType[inletNodes[j].getNodeNum()], Uin, Uout, Uinlet[j]);

  for (int i=0; i<faces.size(); ++i)
    faces[i].template assignFreeStreamValues<dim>(Uin, Uout, U[i]);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::assignPorousWallValues(SVec<double,dim> &Uin, SVec<double,dim> &U)
{

  for (int i=0; i<faces.size(); ++i)
    faces[i].template assignPorousWallValues<dim>(Uin, U[i]);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::setNodeBcValue(double* Vin, SVec<double,dim>& Unode)
{

  for (int i=0; i<nodes.size(); ++i) {
    if (nodeType[i] == BC_INLET_MOVING || nodeType[i] == BC_INLET_FIXED) {
      Unode[i][0] = Vin[0];
      Unode[i][1] = Vin[1];
      Unode[i][2] = Vin[2];
      Unode[i][3] = Vin[3];
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::setNodeBcValue2(double* Uin, SVec<double,dim>& Unode)
{

  for (int i=0; i<nodes.size(); ++i) {
    if (nodeType[i] == BC_INLET_MOVING || nodeType[i] == BC_INLET_FIXED ||
        nodeType[i] == BC_OUTLET_MOVING || nodeType[i] == BC_OUTLET_FIXED) {
      Unode[i][0] = Uin[0];
      Unode[i][1] = Uin[1];
      Unode[i][2] = Uin[2];
      Unode[i][3] = Uin[3];
      Unode[i][4] = Uin[4];
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeFaceBcValue(SVec<double,dim> &Unode, SVec<double,dim> &Uface)
{

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeFaceBcValue(Unode, Uface[i]);

}

//------------------------------------------------------------------------------
// compute node values of the last (dim2-1) face values of Uface

template<int dim1, int dim2>
void SubDomain::computeNodeBcValue(SVec<double,3> &X, SVec<double,dim1> &Uface,
				   SVec<double,dim2> &Unode)
{

  Unode = 0.0;

	if (sampleMesh) {
		int i;
		for (int iFace=0; iFace<faces.getNumSampledFaces(); ++iFace) {
			i = faces.facesConnectedToSampleNode[iFace];
			faces[i].template computeNodeBcValue<dim1,dim2>(X, Uface[i], Unode);
		}

	}
	else {
		for (int i=0; i<faces.size(); ++i)
			faces[i].template computeNodeBcValue<dim1,dim2>(X, Uface[i], Unode);
	}

}

//------------------------------------------------------------------------------
// compute node values of the last (dim2-1) face values of Uface

// Included (MB)
template<int dim1, int dim2>
void SubDomain::computeDerivativeOfNodeBcValue(SVec<double,3> &X, SVec<double,3> &dX, SVec<double,dim1> &Uface, SVec<double,dim1> &dUface,
				   SVec<double,dim2> &dUnode)
{

  dUnode = 0.0;

  for (int i=0; i<faces.size(); ++i)
    faces[i].template computeDerivativeOfNodeBcValue<dim1,dim2>(X, dX, Uface[i], dUface[i], dUnode);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeNodeBCsWallValues(SVec<double,3> &X, SVec<double,1> &dNormSA, SVec<double,dim> &dUfaceSA, SVec<double,dim> &dUnodeSA)
{

  dUnodeSA = 0.0;
  dNormSA = 0.0;

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeNodeBCsWallValues(X, dNormSA, dUfaceSA[i], dUnodeSA);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeNodalForce(PostFcn *postFcn, BcData<dim> &bcData,
				  GeoState &geoState, SVec<double,3> &X,
				  SVec<double,dim> &V, Vec<double> &Pin,
				  SVec<double,3> &F)
{

  F = 0.0;

  Vec<double> &d2wall = geoState.getDistanceToWall();
  SVec<double,dim> &Vwall = bcData.getFaceStateVector();

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeNodalForce(elems, postFcn, X, d2wall, Vwall[i], V, Pin[i], F, gradP);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeDerivativeOfNodalForce(PostFcn *postFcn, BcData<dim> &bcData,
				  GeoState &geoState, SVec<double,3> &X, SVec<double,3> &dX,
				  SVec<double,dim> &V, SVec<double,dim> &dV, Vec<double> &Pin,
				  double dS[3], SVec<double,3> &dF)
{

  dF = 0.0;

  Vec<double> &d2wall = geoState.getDistanceToWall();
  SVec<double,dim> &Vwall = bcData.getFaceStateVector();
  SVec<double,dim> &dVwall = bcData.getdFaceStateVector();

// Remark: Maybe it is necessary to transform conservative in primitive variables.

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeDerivativeOfNodalForce(elems, postFcn, X, dX, d2wall, Vwall[i], dVwall[i], V, dV, Pin[i], dS, dF, gradP, dGradP);

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void SubDomain::computeDerivativeOfNodalForce(RectangularSparseMat<double,3,3> *dForcedX,
                                              RectangularSparseMat<double,3,3> *dForcedGradP,
                                              RectangularSparseMat<double,dim,3> *dForcedV,
                                              RectangularSparseMat<double,3,3> *dForcedS,
                                              SVec<double,3> &dX, SVec<double,dim> &dV, 
                                              double dS[3], SVec<double,3> &dF,
                                              SVec<double,3> &dSSVec, SVec<double,3> &dGradPSVec)
{

  dF = 0.0;
  SVec<double,3> dummy(dF);
  dForcedX->apply(dX, dummy);
  dF += dummy;

  for (int i=0;i<nodes.size();i++) {
    dGradPSVec[i][0] = dGradP[0][i];
    dGradPSVec[i][1] = dGradP[1][i];
    dGradPSVec[i][2] = dGradP[2][i];
  }

  dForcedGradP->apply(dGradPSVec, dummy);
  dF += dummy;
  dForcedV->apply(dV, dummy);
  dF += dummy;

  dSSVec[0][0] = dS[0];
  dSSVec[0][1] = dS[1];
  dSSVec[0][2] = dS[2];

  dForcedS->apply(dSSVec, dummy);
  dF += dummy;

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void SubDomain::computeTransposeDerivativeOfNodalForce(RectangularSparseMat<double,3,3> *dForcedX,
                                              RectangularSparseMat<double,3,3> *dForcedGradP,
                                              RectangularSparseMat<double,dim,3> *dForcedV,
                                              RectangularSparseMat<double,3,3> *dForcedS,
                                              SVec<double,3> &dF, SVec<double,3> &dGradPSVec, 
                                              SVec<double,3> &dX, SVec<double,dim> &dV, SVec<double,3> dSSVec)
{
  SVec<double,3> dummydX(dX); 
  dForcedX->applyTranspose(dF, dummydX);
  dX += dummydX;

  SVec<double,3> dummydGradPSVec(dGradPSVec); 
  dForcedGradP->applyTranspose(dF, dummydGradPSVec);
  dGradPSVec += dummydGradPSVec;

  SVec<double,dim> dummydV(dV);
  dForcedV->applyTranspose(dF, dummydV);
  dV += dummydV;

  SVec<double,3> dummydSSVec(dSSVec);
  dForcedS->applyTranspose(dF, dummydSSVec);
  dSSVec += dummydSSVec;

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void SubDomain::computeDerivativeOperatorsOfNodalForce(PostFcn *postFcn, SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &Pin,
                                                       RectangularSparseMat<double,3,3> &dForcedX,
                                                       RectangularSparseMat<double,3,3> &dForcedGradP,
                                                       RectangularSparseMat<double,dim,3> &dForcedV,
                                                       RectangularSparseMat<double,3,3> &dForcedS)
{

  for (int i=0; i<faces.size(); ++i)
	  faces[i].computeDerivativeOperatorsOfNodalForce(elems, postFcn, X, V, Pin[i], gradP, dForcedX, dForcedGradP, dForcedV, dForcedS);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeNodalHeatPower(PostFcn* postFcn, BcData<dim>& bcData,
				      GeoState& geoState, SVec<double,3>& X,
				      SVec<double,dim>& V, Vec<double>& P)
{

  P = 0.0;

  Vec<double>& d2wall = geoState.getDistanceToWall();
  SVec<double,dim>& Vwall = bcData.getFaceStateVector();

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeNodalHeatPower(elems, postFcn, X, d2wall, Vwall[i], V, P);

}

//------------------------------------------------------------------------------
template<int dim>
void SubDomain::computeNodalHeatFluxRelatedValues(PostFcn* postFcn, BcData<dim>& bcData,
                                      GeoState& geoState, SVec<double,3>& X,
                                      SVec<double,dim>& V, Vec<double>& P, Vec<double>& N, bool includeKappa)
{

  P = 0.0;
  N = -1.0;

  Vec<double>& d2wall = geoState.getDistanceToWall();
  SVec<double,dim>& Vwall = bcData.getFaceStateVector();

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeNodalHeatFluxRelatedValues(elems, postFcn, X, d2wall, Vwall[i], V, P, N, includeKappa);

}

//------------------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeDerivativeOfNodalHeatPower(PostFcn* postFcn, BcData<dim>& bcData,
				      GeoState& geoState, SVec<double,3>& X, SVec<double,3>& dX,
				      SVec<double,dim>& V, SVec<double,dim>& dV, double dS[3], Vec<double>& dP)
{

  dP = 0.0;

  Vec<double>& d2wall = geoState.getDistanceToWall();
  SVec<double,dim>& Vwall = bcData.getFaceStateVector();
  SVec<double,dim>& dVwall = bcData.getdFaceStateVector();

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeDerivativeOfNodalHeatPower(elems, postFcn, X, dX, d2wall, Vwall[i], dVwall[i], V, dV, dS, dP);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeForceAndMoment(map<int,int> & surfOutMap, PostFcn *postFcn, BcData<dim> &bcData,
				      GeoState &geoState, SVec<double,3> &X,
				      SVec<double,dim> &V, Vec3D &x0, Vec3D *Fi,
				      Vec3D *Mi, Vec3D *Fv, Vec3D *Mv, int hydro,
				      SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX, Vec<double> *genCF)
{

  Vec<double> &d2wall = geoState.getDistanceToWall();
  SVec<double,dim> &Vwall = bcData.getFaceStateVector();

  for (int i=0; i<faces.size(); ++i) {
    int idx;
    map<int,int>::iterator it = surfOutMap.find(faces[i].getSurfaceID());
    if(it != surfOutMap.end() && it->second != -2)
      idx = it->second;
    else {
      if(faces[i].getCode() == BC_ISOTHERMAL_WALL_MOVING ||
         faces[i].getCode() == BC_ADIABATIC_WALL_MOVING  ||
         faces[i].getCode() == BC_SLIP_WALL_MOVING ||
         faces[i].getCode() == BC_POROUS_WALL_MOVING)
        idx = 0;
      else
        idx = -1;
    }

    if(idx >= 0)  {
      faces[i].computeForceAndMoment(elems, postFcn, X, d2wall, Vwall[i], V, x0,
                       Fi[idx], Mi[idx], Fv[idx], Mv[idx], gradP, hydro, mX, genCF);

    }
  }

}

//------------------------------------------------------------------------------
// KW: FS Riemann based force calculation.
template<int dim>
void SubDomain::computeForceAndMoment(ExactRiemannSolver<dim> &riemann, VarFcn *varFcn, 
                                      map<int,int> & surfOutMap, PostFcn *postFcn, BcData<dim> &bcData,
                                      GeoState &geoState, SVec<double,3> &X,
                                      SVec<double,dim> &V, Vec3D &x0, Vec3D *Fi,
                                      Vec3D *Mi, Vec3D *Fv, Vec3D *Mv, int hydro,
                                      SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX, Vec<double> *genCF)
{

  Vec<double> &d2wall = geoState.getDistanceToWall();
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();

  SVec<double,dim> &Vwall = bcData.getFaceStateVector();

  for (int i=0; i<faces.size(); ++i) {
    int idx;
    map<int,int>::iterator it = surfOutMap.find(faces[i].getSurfaceID());
    if(it != surfOutMap.end() && it->second != -2)
      idx = it->second;
    else {
      if(faces[i].getCode() == BC_ISOTHERMAL_WALL_MOVING ||
         faces[i].getCode() == BC_ADIABATIC_WALL_MOVING  ||
         faces[i].getCode() == BC_SLIP_WALL_MOVING ||
         faces[i].getCode() == BC_POROUS_WALL_MOVING)
        idx = 0;
      else
        idx = -1;
    }

    if(idx >= 0)  {
      faces[i].computeForceAndMoment(riemann, varFcn, n, ndot, elems, postFcn, X, d2wall, Vwall[i], V, x0,
                       Fi[idx], Mi[idx], Fv[idx], Mv[idx], gradP, hydro, mX, genCF);
    }
  }

}

/*
template<int dim>
void SubDomain::computeLiftSurfaces(map<int,int> & surfOutMap, PostFcn *postFcn, BcData<dim> &bcData,
				      GeoState &geoState, SVec<double,3> &X,
				      SVec<double,dim> &V, Vec3D &x0, Vec3D *Fi,
				      Vec3D *Mi, Vec3D *Fv, Vec3D *Mv, int hydro,
                                      SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX, Vec<double> *genCF)
{

  Vec<double> &d2wall = geoState.getDistanceToWall();
  SVec<double,dim> &Vwall = bcData.getFaceStateVector();

  for (int i=0; i<faces.size(); ++i) {
    int idx;
    map<int,int>::iterator it = surfOutMap.find(faces[i].getSurfaceID());
    if(it != surfOutMap.end() && it->second != -2)
      idx = it->second;
    else {
      if(faces[i].getCode() == BC_ISOTHERMAL_WALL_MOVING ||
         faces[i].getCode() == BC_ADIABATIC_WALL_MOVING  ||
         faces[i].getCode() == BC_SLIP_WALL_MOVING ||
         faces[i].getCode() == BC_POROUS_WALL_MOVING)
        idx = 0;
      else
        idx = -1;
    }

    if(idx >= 0)  {
      faces[i].computeForceAndMoment(elems, postFcn, X, d2wall, Vwall[i], V, x0,
                       Fi[idx], Mi[idx], Fv[idx], Mv[idx], gradP, hydro, mX, genCF);
    }
  }

}
*/
//------------------------------------------------------------------------------
template<int dim>
void SubDomain::computeHeatFluxes(map<int,int> & surfOutMapHF, PostFcn* postFcn, BcData<dim>& bcData,
                                      GeoState& geoState, SVec<double,3>& X,
                                      SVec<double,dim>& V, double* HF)
{
  Vec<double>& d2wall = geoState.getDistanceToWall();
  SVec<double,dim>& Vwall = bcData.getFaceStateVector();

  for (int i=0; i<faces.size(); ++i){
    int idx;
    map<int,int>::iterator it = surfOutMapHF.find(faces[i].getSurfaceID());
    if(it != surfOutMapHF.end() && it->second != -2)
      idx = it->second;
    else {
      if(faces[i].getCode() == BC_ISOTHERMAL_WALL_MOVING)  
        idx = 0;
      else
        idx = -1;
    }
    if(idx >= 0)  {
   double hp = faces[i].computeHeatFluxes(elems, postFcn, X, d2wall, Vwall[i], V);
    HF[idx] += hp;
    }
  }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::computeDerivativeOfForceAndMoment(map<int,int> & surfOutMap, PostFcn *postFcn, BcData<dim> &bcData,
		                                          GeoState &geoState, SVec<double,3> &X, SVec<double,3> &dX,
		                                          SVec<double,dim> &V, SVec<double,dim> &dV, double dS[3],
		                                          Vec3D &x0, Vec3D *dFi, Vec3D *dMi, Vec3D *dFv, Vec3D *dMv, int hydro)
{

  Vec<double> &d2wall = geoState.getDistanceToWall();
  SVec<double,dim> &Vwall = bcData.getFaceStateVector();
  SVec<double,dim> &dVwall = bcData.getdFaceStateVector();

  for (int i=0; i<faces.size(); ++i) {
    int idx;
    map<int,int>::iterator it = surfOutMap.find(faces[i].getSurfaceID());
    if(it != surfOutMap.end() && it->second != -2)
      idx = it->second;
    else {
      if(faces[i].getCode() == BC_ISOTHERMAL_WALL_MOVING ||
         faces[i].getCode() == BC_ADIABATIC_WALL_MOVING  ||
	 faces[i].getCode() == BC_SLIP_WALL_MOVING ||
	 faces[i].getCode() == BC_POROUS_WALL_MOVING)
        idx = 0;
      else
        idx = -1;
    }

    if(idx >= 0){
      faces[i].computeDerivativeOfForceAndMoment(elems, postFcn, X, dX, d2wall, Vwall[i], dVwall[i], V, dV, dS, x0, dFi[idx], dMi[idx], dFv[idx], dMv[idx], gradP, dGradP, hydro);
    }
  }

}




// Included (YC)
template<int dim>
void SubDomain::computeDerivativeOfForceAndMoment(RectangularSparseMat<double,3,3> *dFidGradP,
                                                  RectangularSparseMat<double,3,3> *dFidX,
                                                  RectangularSparseMat<double,dim,3> *dFidV,
                                                  RectangularSparseMat<double,3,3> *dFvdX,
                                                  RectangularSparseMat<double,dim,3> *dFvdV,
                                                  RectangularSparseMat<double,3,3> *dFidS,
                                                  RectangularSparseMat<double,3,3> *dMidGradP,
                                                  RectangularSparseMat<double,3,3> *dMidX,
                                                  RectangularSparseMat<double,dim,3> *dMidV,
                                                  RectangularSparseMat<double,3,3> *dMidS,
                                                  RectangularSparseMat<double,3,3> *dMvdX,
                                                  RectangularSparseMat<double,dim,3> *dMvdV,
                                                  SVec<double,3> &dX,
                                                  SVec<double,dim> &dV, double dS[3], SVec<double,3> &dGradPSVec,
                                                  Vec3D *dFi, Vec3D *dMi, Vec3D *dFv, Vec3D *dMv, int hydro)
{

  if(hydro != 0) { fprintf(stderr, " *** Error: hydro must be zero for sparse format\n");  exit(-1); }

  SVec<double,3> dFiSVec(1), dummy(1), dSSVec(1), dFvSVec(1);
  dFiSVec[0][0] = dFi[0][0];
  dFiSVec[0][1] = dFi[0][1];
  dFiSVec[0][2] = dFi[0][2];
  dFvSVec[0][0] = dFv[0][0];
  dFvSVec[0][1] = dFv[0][1];
  dFvSVec[0][2] = dFv[0][2];
  dSSVec[0][0]  = dS[0];
  dSSVec[0][1]  = dS[1];
  dSSVec[0][2]  = dS[2];

  dFidV->apply(dV, dummy);
  dFiSVec += dummy;
  dFidGradP->apply(dGradPSVec, dummy);
  dFiSVec += dummy;
  dFidX->apply(dX, dummy);
  dFiSVec += dummy;
  dFidS->apply(dSSVec, dummy);
  dFiSVec += dummy;
  dFvdV->apply(dV, dummy);
  dFvSVec += dummy;
  dFvdX->apply(dX, dummy);
  dFvSVec += dummy;

  dFi[0][0] = dFiSVec[0][0];
  dFi[0][1] = dFiSVec[0][1];
  dFi[0][2] = dFiSVec[0][2];
  dFv[0][0] = dFvSVec[0][0];
  dFv[0][1] = dFvSVec[0][1];
  dFv[0][2] = dFvSVec[0][2];

  SVec<double,3> dMiSVec(1);
  dMiSVec[0][0] = dMi[0][0];
  dMiSVec[0][1] = dMi[0][1];
  dMiSVec[0][2] = dMi[0][2];
  dMidV->apply(dV, dummy);
  dMiSVec += dummy;
  dMidGradP->apply(dGradPSVec, dummy);
  dMiSVec += dummy;
  dMidX->apply(dX, dummy);
  dMiSVec += dummy;
  dMidS->apply(dSSVec, dummy);
  dMiSVec += dummy;

  dMi[0][0] = dMiSVec[0][0];
  dMi[0][1] = dMiSVec[0][1];
  dMi[0][2] = dMiSVec[0][2];

  SVec<double,3> dMvSVec(1);
  dMvSVec[0][0] = dMv[0][0];
  dMvSVec[0][1] = dMv[0][1];
  dMvSVec[0][2] = dMv[0][2];
  dMvdX->apply(dX, dummy);
  dMvSVec += dummy;
  dMvdV->apply(dV, dummy);
  dMvSVec += dummy;
  dMv[0][0] = dMvSVec[0][0];
  dMv[0][1] = dMvSVec[0][1];
  dMv[0][2] = dMvSVec[0][2];

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void SubDomain::computeTransposeDerivativeOfForceAndMoment(RectangularSparseMat<double,3,3> *dFidGradP,
                                                           RectangularSparseMat<double,3,3> *dFidX,
                                                           RectangularSparseMat<double,dim,3> *dFidV,
                                                           RectangularSparseMat<double,3,3> *dFvdX,
                                                           RectangularSparseMat<double,dim,3> *dFvdV,
                                                           RectangularSparseMat<double,3,3> *dFidS,
                                                           RectangularSparseMat<double,3,3> *dMidGradP,
                                                           RectangularSparseMat<double,3,3> *dMidX,
                                                           RectangularSparseMat<double,dim,3> *dMidV,
                                                           RectangularSparseMat<double,3,3> *dMidS,
                                                           RectangularSparseMat<double,3,3> *dMvdX,
                                                           RectangularSparseMat<double,dim,3> *dMvdV,
                                                           SVec<double,3> &dFiSVec, SVec<double,3> &dFvSVec,
                                                           SVec<double,3> &dMiSVec, SVec<double,3> &dMvSVec, SVec<double,3> &dX,
                                                           SVec<double,dim> &dV, SVec<double,3> &dSSVec,
                                                           SVec<double,3> &dGradPSVec, int hydro)
{

  if(hydro != 0) { fprintf(stderr, " *** Error: hydro must be zero for sparse format\n");  exit(-1); }


  SVec<double,3> dummy3(dGradPSVec);
  SVec<double,dim> dVdummy(dV);
  dFidV->applyTranspose(dFiSVec, dVdummy);
  dV += dVdummy;
  dFvdV->applyTranspose(dFvSVec, dVdummy);
  dV += dVdummy;
  dFidGradP->applyTranspose(dFiSVec, dummy3);
  dGradPSVec += dummy3;
  dFidX->applyTranspose(dFiSVec, dummy3);
  dX += dummy3;
  dFvdX->applyTranspose(dFvSVec, dummy3);
  dX += dummy3;
  dFidS->applyTranspose(dFiSVec, dummy3);
  dSSVec += dummy3;

  dMidV->applyTranspose(dMiSVec, dVdummy);
  dV += dVdummy;
  dMvdV->applyTranspose(dMvSVec, dVdummy);
  dV += dVdummy;
  dMidGradP->applyTranspose(dMiSVec, dummy3);
  dGradPSVec += dummy3;
  dMidX->applyTranspose(dMiSVec, dummy3);
  dX += dummy3;
  dMidS->applyTranspose(dMiSVec, dummy3);
  dSSVec += dummy3;

  dMvdX->applyTranspose(dMvSVec, dummy3);
  dX += dummy3;

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeDerivativeOperatorsOfForceAndMoment(map<int,int> & surfOutMap, PostFcn *postFcn, BcData<dim> &bcData,
                                                           GeoState &geoState, SVec<double,3> &X,
                                                           SVec<double,dim> &V, Vec3D &x0, int hydro,
                                                           RectangularSparseMat<double,3,3> &dFidGradP,
                                                           RectangularSparseMat<double,3,3> &dFidX,
                                                           RectangularSparseMat<double,dim,3> &dFidV,
                                                           RectangularSparseMat<double,3,3> &dFvdX,
                                                           RectangularSparseMat<double,dim,3> &dFvdV,
                                                           RectangularSparseMat<double,3,3> &dFidS,
                                                           RectangularSparseMat<double,3,3> &dMidGradP,
                                                           RectangularSparseMat<double,3,3> &dMidX,
                                                           RectangularSparseMat<double,dim,3> &dMidV,
                                                           RectangularSparseMat<double,3,3> &dMidS,
                                                           RectangularSparseMat<double,3,3> &dMvdX,
                                                           RectangularSparseMat<double,dim,3> &dMvdV)
{

  Vec<double> &d2wall = geoState.getDistanceToWall();
  SVec<double,dim> &Vwall = bcData.getFaceStateVector();

  for (int i=0; i<faces.size(); ++i) {
    int idx;
    map<int,int>::iterator it = surfOutMap.find(faces[i].getSurfaceID());
    if(it != surfOutMap.end() && it->second != -2)
      idx = it->second;
    else {
      if(faces[i].getCode() == BC_ISOTHERMAL_WALL_MOVING ||
         faces[i].getCode() == BC_ADIABATIC_WALL_MOVING  ||
         faces[i].getCode() == BC_SLIP_WALL_MOVING ||
         faces[i].getCode() == BC_POROUS_WALL_MOVING)
        idx = 0;
      else
        idx = -1;
    }

    if(idx >= 0)
      faces[i].computeDerivativeOperatorsOfForceAndMoment(elems, postFcn, X, d2wall, Vwall[i], V, x0, gradP, hydro,
                                                          dFidGradP, dFidX, dFidV, dFvdX, dFvdV, dFidS, dMidGradP, dMidX, dMidV, dMidS, dMvdX, dMvdV);
  }

}

//------------------------------------------------------------------------------

template<int dim>
double SubDomain::computeInterfaceWork(PostFcn* postFcn, BcData<dim>& bcData,
				       GeoState& geoState, SVec<double,3>& X,
				       SVec<double,dim>& V, Vec<double>& Pin)
{

  Vec<double>& ndot = geoState.getFaceNormalVel();
  Vec<double>& d2wall = geoState.getDistanceToWall();
  SVec<double,dim>& Vwall = bcData.getFaceStateVector();

  double E = 0.0;
  for (int i=0; i<faces.size(); ++i)
    E += faces[i].computeInterfaceWork(elems, postFcn, X, d2wall, ndot[i], Vwall[i], V, Pin[i]);

  return E;

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeFaceScalarQuantity(PostFcn::ScalarType type, PostFcn *postFcn,
					  BcData<dim> &bcData, GeoState &geoState,
					  SVec<double,3> &X, SVec<double,dim> &V,
					  SVec<double,2> &Q)
{

  Q = 0.0;

  Vec<double> &d2wall = geoState.getDistanceToWall();
  SVec<double,dim> &Vwall = bcData.getFaceStateVector();

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeScalarQuantity(type, elems, postFcn, X, d2wall, Vwall[i], V, Q);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeNodeScalarQuantity(PostFcn::ScalarType type, PostFcn *postFcn,
					  SVec<double,dim> &V, SVec<double,3> &X,
					  Vec<double> &Q)
{
  double phi = 1.0;
  int fluidId = 0;
  for (int i=0; i<Q.size(); ++i)
    Q[i] = postFcn->computeNodeScalarQuantity(type, V[i], X[i],fluidId,&phi);
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void
SubDomain::computeDerivativeOfNodeScalarQuantity(PostFcn::ScalarDerivativeType type, PostFcn *postFcn, double dS[3], SVec<double,dim> &V, SVec<double,dim> &dV, SVec<double,3> &X, SVec<double,3> &dX, Vec<double> &dQ)
{

  for (int i=0; i<dQ.size(); ++i)
    dQ[i] = postFcn->computeDerivativeOfNodeScalarQuantity(type, dS, V[i], dV[i], X[i], dX[i], 1.0);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeXP(PostFcn *postFcn, SVec<double,dim> &V, SVec<double,3> &X, Vec<double> &Q, int dir)
{
  double phi = 1.0;
  int fluidId = 0;
  for (int i=0; i<Q.size(); ++i) {
    if (nodeType[i] == BC_ADIABATIC_WALL_MOVING || nodeType[i] == BC_ISOTHERMAL_WALL_MOVING)  {
      Q[i] = postFcn->computeNodeScalarQuantity(PostFcn::DIFFPRESSURE, V[i], X[i], &phi, fluidId);
      Q[i] *= X[i][dir];
    }
  }

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
inline void SubDomain::computeNodeScalarQuantity(PostFcn::ScalarType type, PostFcn *postFcn,
                                         SVec<double,dim> &V, SVec<double,3> &X,
                                          Vec<double> &Q, Vec<int> &fluidId,SVec<double,dimLS>* phi) 
{
  
  if (phi) {
    for (int i=0; i<Q.size(); ++i)
      Q[i] = postFcn->computeNodeScalarQuantity(type, V[i], X[i], fluidId[i],(*phi)[i]);
  } else { 
    for (int i=0; i<Q.size(); ++i)
      Q[i] = postFcn->computeNodeScalarQuantity(type, V[i], X[i], fluidId[i],NULL);
  }
}

template<int dim, int dimLS>
inline double SubDomain::computeNodeScalarQuantity(PostFcn::ScalarType type, PostFcn *postFcn,
						   SVec<double,dim> &V, SVec<double,3> &X,
						   Vec<int> &fluidId,int i,SVec<double,dimLS>* phi) 
{
  
  if (phi) {
    return postFcn->computeNodeScalarQuantity(type, V[i], X[i], fluidId[i],(*phi)[i]);
  } else { 
    return postFcn->computeNodeScalarQuantity(type, V[i], X[i], fluidId[i],NULL);
  }
}

//------------------------------------------------------------------------------

template<class S1, class S2>
void SubDomain::computeStiffAndForce(DefoMeshMotionData::Element typeElement,
				     SVec<double,3>& X, SVec<double,3>& F, GenMat<S1,3>& K, GenMat<S2,3>* P,
                                     double volStiff, int* ndType)
{
  const int MaxSize = (3*Elem::MaxNumNd);
  double kEl[MaxSize*MaxSize];
  double fEl[MaxSize];
		double minVolume = 10000.0;
  int minElemNum = 0;

  F = 0.0;
  K = 0.0;
  if (P) *P = 0.0;

  int i;
  for (i=0; i<elems.size(); i++)  {
    int j;
    double *fEl_loc;

    // Compute stiffness depending on type of structural analogy
    switch (typeElement) {

    case DefoMeshMotionData::LINEAR_FE : {
      elems[i].computeStiffAndForceLIN(kEl, X, nodes);
      break;
    }

    case DefoMeshMotionData::NON_LINEAR_FE : {
      double vol = elems[i].computeVolume(X); 
      if(minVolume > vol) {
								minElemNum = i+1;
        minVolume = vol;	
      }
      elems[i].computeStiffAndForce(fEl, kEl, X, nodes, volStiff);
      for (j=0, fEl_loc = fEl;
       	   j<elems[i].numNodes();
	          j++, fEl_loc+=3) {
	       F[ elems[i][j] ][0] -= fEl_loc[0];
	       F[ elems[i][j] ][1] -= fEl_loc[1];
	       F[ elems[i][j] ][2] -= fEl_loc[2];
      }
      break;
    }

    case DefoMeshMotionData::TORSIONAL_SPRINGS : {
      elems[i].computeStiffTorsionSpring(kEl, X, volStiff);
      break;
    }

    case DefoMeshMotionData::BALL_VERTEX : {
      elems[i].computeStiffBallVertex(kEl, X, nodes, volStiff);
      break;
    }
      
    case DefoMeshMotionData::NL_BALL_VERTEX : {
      elems[i].computeStiffAndForceBallVertex(fEl, kEl, X, nodes, volStiff);
      for (j=0, fEl_loc = fEl; j<elems[i].numNodes(); j++, fEl_loc+=3) {
	       F[ elems[i][j] ][0] -= fEl_loc[0];
       	F[ elems[i][j] ][1] -= fEl_loc[1];
       	F[ elems[i][j] ][2] -= fEl_loc[2];
      }
      break;
    }

    }

    // Add contribution to global matrix
    K.addContrib(elems[i].numNodes(), elems[i], kEl);
    if (P)
      P->addContrib(elems[i].numNodes(), elems[i], kEl);

  }
//		fprintf(stderr,"element %d with minimum volume of %6.3e.\n", minElemNum, minVolume);  

  if(ndType){
    for (i = 0; i < nodes.size(); ++i) {
      if (ndType[i] != BC_INTERNAL ) {

        F[i][0] = 0.0;
        F[i][1] = 0.0;
        F[i][2] = 0.0;
      }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
int SubDomain::checkSolution(VarFcn *varFcn, SVec<double,dim> &U, LevelSetStructure *LSS)
{

  bool vflag;
  int ierr = 0;
  int pclipping = 0;
  int rhoclipping = 0;
  int temp;
  double V[dim];
  double rho,p;

	for(int i=0; i<U.size(); i++) 
	{
		if(LSS && !LSS->isActive(0.0, i)) continue;

		if((vflag = varFcn->doVerification()))
		{
      temp = varFcn->conservativeToPrimitiveVerification(locToGlobNodeMap[i]+1,U[i],V);
      rhoclipping += temp % 2;
      pclipping += temp/2;
    }
		else 
		{
      varFcn->conservativeToPrimitive(U[i], V);
    }
    // Even if doVerification returns true, we still check for negative density and pressure in
    // case CheckSolution is used. However, by convention we don't print the error message to stderr.
    rho = varFcn->getDensity(V);
    p = varFcn->checkPressure(V);

		if(rho <= 0.0) 
		{
      if(!vflag) fprintf(stderr, "*** Error: negative density (%e) for node %d\n",
                         rho, locToGlobNodeMap[i] + 1);
      ++ierr;
    }

		if(p <= 0.0) 
		{
      if(!vflag) fprintf(stderr, "*** Error: negative pressure (%e) for node %d\n",
                         p, locToGlobNodeMap[i] + 1);
      ++ierr;
    }
    // Check for abnormally large velocities, which may be the result of an instability
		if(fabs(V[1]) > 1e6 || fabs(V[2]) > 1e6 || fabs(V[3]) > 1e6) 
		{
      errorHandler->localErrors[ErrorHandler::LARGE_VELOCITY] += 1;
      fprintf(stderr,"*** Warning: Abnormally large velocity: [%lf, %lf, %lf] detected at node %d."
                     " This may be a symptom of an instability\n",V[1],V[2],V[3],locToGlobNodeMap[i] + 1);
    }
  }

  errorHandler->localErrors[ErrorHandler::PRESSURE_CLIPPING] += pclipping;
  errorHandler->localErrors[ErrorHandler::DENSITY_CLIPPING] += rhoclipping;
  errorHandler->localErrors[ErrorHandler::UNPHYSICAL] += ierr;
  return ierr;
}

//------------------------------------------------------------------------------

template<int dim>
int SubDomain::checkSolution(VarFcn *varFcn, SVec<double,dim> &U, Vec<int> &fluidId, LevelSetStructure *LSS)
{

  bool vflag;
  int ierr = 0;
  int pclipping = 0;
  int rhoclipping = 0;
  int temp;
  double V[dim];
  double rho,p;

	for(int i=0; i<U.size(); i++) 
	{
		if(LSS && !LSS->isActive(0.0, i)) continue;

		if((vflag = varFcn->doVerification())) 
		{
      temp = varFcn->conservativeToPrimitiveVerification(locToGlobNodeMap[i]+1, U[i], V, fluidId[i]);
      rhoclipping += temp % 2;
      pclipping += temp/2;
    }   
		else 
		{
      varFcn->conservativeToPrimitive(U[i], V, fluidId[i]);
    }
    // Even if doVerification returns true, now we still check for negative density and pressure in
    // case CheckSolution is used. However, by convention we don't print the error message to stderr.
    rho = varFcn->getDensity(V, fluidId[i]);
    p = varFcn->checkPressure(V, fluidId[i]);

		if(rho <= 0.0) 
		{
      if(!vflag) fprintf(stderr, "*** Error: negative density (%e) for node %d with fluidId=%d\n",
                         rho, locToGlobNodeMap[i] + 1, fluidId[i]);
      ++ierr;
    }
		if(p <= 0.0) 
		{
      if(!vflag) fprintf(stderr, "*** Error: negative pressure (%e) for node %d with fluidId=%d\n",
                         p, locToGlobNodeMap[i] + 1, fluidId[i]);
      ++ierr;
    }
    // Check for abnormally large velocities, which may be the result of an instability
    if (fabs(V[1]) > 1e6 || fabs(V[2]) > 1e6 || fabs(V[3]) > 1e6) {
      errorHandler->localErrors[ErrorHandler::LARGE_VELOCITY] += 1;
      fprintf(stderr,"*** Warning: Abnormally large velocity: [%lf, %lf, %lf] detected at node %d with fluidId=%d."
                     " This may be a symptom of an instability\n",V[1],V[2],V[3],locToGlobNodeMap[i]+1,fluidId[i]);
    }
  }

  errorHandler->localErrors[ErrorHandler::PRESSURE_CLIPPING] += pclipping;
  errorHandler->localErrors[ErrorHandler::DENSITY_CLIPPING] += rhoclipping;
  errorHandler->localErrors[ErrorHandler::UNPHYSICAL] += ierr;

  return ierr;
}

//------------------------------------------------------------------------------

template<int dim>
int SubDomain::checkSolution(VarFcn *varFcn, Vec<double> &ctrlVol, SVec<double,dim> &U,
                             Vec<int> &fluidId, Vec<int> &fluidIdn)
{
  bool vflag;
  int ierr = 0;
  int pclipping = 0;
  int rhoclipping = 0;
  int temp;
  double V[dim];
  double rho,p;

  for(int i=0; i<U.size(); i++) {
    if((vflag = varFcn->doVerification())) {
      temp = varFcn->conservativeToPrimitiveVerification(locToGlobNodeMap[i]+1, U[i], V, fluidId[i]);
      rhoclipping += temp % 2;
      pclipping += temp/2;
    }   
    else {
      varFcn->conservativeToPrimitive(U[i], V, fluidId[i]);
    }
    // Even if doVerification returns true, now we still check for negative density and pressure in
    // case CheckSolution is used. However, by convention we don't print the error message to stderr.
    rho = varFcn->getDensity(V, fluidId[i]);
    p = varFcn->checkPressure(V, fluidId[i]);
    if (rho <= 0.0) {
      if(!vflag) fprintf(stderr, "*** Error: negative density (%e) for node %d with fluidID=%d (previously %d)\n",
                         rho, locToGlobNodeMap[i] + 1, fluidId[i], fluidIdn[i]);
      ++ierr;
    }
    if (p <= 0.0) {
      if(!vflag) fprintf(stderr, "*** Error: negative pressure (%e) for node %d with fluidID=%d (previously %d)\n",
                         p, locToGlobNodeMap[i] + 1, fluidId[i], fluidIdn[i]);
      ++ierr;
    }
    // Check for abnormally large velocities, which may be the result of an instability
    if (fabs(V[1]) > 1e6 || fabs(V[2]) > 1e6 || fabs(V[3]) > 1e6) {
      errorHandler->localErrors[ErrorHandler::LARGE_VELOCITY] += 1;
      fprintf(stderr,"*** Warning: Abnormally large velocity: [%lf, %lf, %lf] detected at node %d with fluidID=%d (previously %d)."
                     " This may be a symptom of an instability\n",V[1],V[2],V[3],locToGlobNodeMap[i]+1,fluidId[i],fluidIdn[i]);
    }
  }

  errorHandler->localErrors[ErrorHandler::PRESSURE_CLIPPING] += pclipping;
  errorHandler->localErrors[ErrorHandler::DENSITY_CLIPPING] += rhoclipping;
  errorHandler->localErrors[ErrorHandler::UNPHYSICAL] += ierr;

  return ierr;
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::restrictionOnPhi(SVec<double,dim> &initial, Vec<int> &fluidId,
                                 SVec<double,dim> &restriction, int fluidIdTarget){

  int idim;
  restriction = 0.0;
  for (int i=0; i<nodes.size(); i++)
    if(fluidId[i] == fluidIdTarget)
      for(idim=0; idim<dim; idim++) restriction[i][idim] = initial[i][idim];

}
//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
int SubDomain::fixSolution(VarFcn *varFcn, SVec<double,dim> &U, SVec<double,dim> &dU, Vec<int>* fluidId,int verboseFlag)
{

  int ierr = 0;

  for (int i=0; i<U.size(); ++i) {
    double V[dim];
    double Un[dim];

    for (int j=0; j<dim; ++j)
      Un[j] = U[i][j] + dU[i][j];

      int id = 0;
      if (fluidId)
	id = (*fluidId)[i];
      
      varFcn->conservativeToPrimitive(Un, V,id);
      double rho = varFcn->getDensity(V,id);
      double p = varFcn->checkPressure(V,id);
      
      if (rho <= 0.0) {
      if (verboseFlag == 4)
        fprintf(stderr, "*** Warning: negative density (%e) was fixed for node %d\n", rho, locToGlobNodeMap[i] + 1);

      for (int j=0; j<dim; ++j)
        dU[i][j] = 0.0;

      ++ierr;
    }
    if (p <= 0.0) {
      if (verboseFlag == 4)
        fprintf(stderr, "*** Warning: negative pressure (%e) was fixed for node %d (rho = %e)\n", p, locToGlobNodeMap[i] + 1, rho);

      for (int j=0; j<dim; ++j)
        dU[i][j] = 0.0;

      ++ierr;
    }

  }

  return ierr;

}

// Included (MB)
template<int dim>
int SubDomain::fixSolution2(VarFcn *varFcn, SVec<double,dim> &U, SVec<double,dim> &dU, Vec<int>* fluidId,int verboseFlag)
{

  int ierr = 0;

  for (int i=0; i<U.size(); ++i) {
    double V[dim];
    double Un[dim];

    for (int j=0; j<dim; ++j)
      Un[j] = U[i][j] + dU[i][j];

    int id = 0;
    if (fluidId)
      id = (*fluidId)[i];
      
    varFcn->conservativeToPrimitive(Un, V,id);
    double rho = varFcn->getDensity(V,id);
    double p = varFcn->checkPressure(V,id);
    
    varFcn->conservativeToPrimitive(U[i], V,id);
    double rho0 = varFcn->getDensity(V,id);
    double p0 = varFcn->checkPressure(V,id);

    double rhomin = varFcn->getVarFcnBase(id)->rhomin;  
    double pmin = varFcn->getVarFcnBase(id)->pmin;  
    if ((rhomin < 0.0 && pmin < 0.0) || 
        (rho > rhomin && p > pmin))
      continue;
  
    std::cout << "In fixSolution2 for node " << locToGlobNodeMap[i]+1 << std::endl; 
    double alpha = 1.0,alphamax = 1.0; 
    double alphamin = 0.0;
    while (fabs(alphamax-alphamin) > 1.0e-8) {

      alpha = 0.5*(alphamin+alphamax);
      for (int j=0; j<dim; ++j)
        Un[j] = U[i][j] + alpha*dU[i][j];

      varFcn->conservativeToPrimitive(Un, V,id);
      rho = varFcn->getDensity(V,id);
      p = varFcn->checkPressure(V,id);
      if (p < pmin || rho < rhomin)
        alphamax = alpha;
      else
        alphamin = alpha;
    }
   
    std::cout << "Alpha = " << alpha << std::endl;
    for (int j=0; j<dim; ++j)
      dU[i][j] *= alpha;
    
  }

  return ierr;

}

//------------------------------------------------------------------------------

template<int dim, int neq>
int SubDomain::clipSolution(TsData::Clipping ctype, BcsWallData::Integration wtype,
			    VarFcn* varFcn, double* Uin, bool* flag, SVec<double,dim>& U,
			    int* cmin, int* pmin, double* vmin)
{
  int ierr = 0;
  int pclipping = 0;
  int rhoclipping = 0;
  int temp;
  bool vflag;
  double V[dim];
  double rho,p;
  varFcn->conservativeToPrimitive(U[0], V);

  int k;
  for (k=0; k<neq; ++k) {
    cmin[k] = 0;
    pmin[k] = locToGlobNodeMap[0] + 1;
    vmin[k] = V[dim-neq+k];
  }

  for (int i=0; i<U.size(); ++i) {
    if((vflag = varFcn->doVerification())) {
      temp = varFcn->conservativeToPrimitiveVerification(locToGlobNodeMap[i]+1,U[i],V);
      rhoclipping += temp % 2;
      pclipping += temp/2;
    }
    else {
      varFcn->conservativeToPrimitive(U[i], V);
    }
    rho = varFcn->getDensity(V);
    p = varFcn->checkPressure(V);
    if (rho <= 0.0) {
      if(!vflag) fprintf(stderr, "*** Error: negative density (%e) for node %d\n",
	                 rho, locToGlobNodeMap[i] + 1);
      ++ierr;
    }
    if (p <= 0.0) {
      if(!vflag) fprintf(stderr, "*** Error: negative pressure (%e) for node %d\n",
	                 p, locToGlobNodeMap[i] + 1);
      ++ierr;
    }
    if (fabs(V[1]) > 1e6 || fabs(V[2]) > 1e6 || fabs(V[3]) > 1e6) {
      errorHandler->localErrors[ErrorHandler::LARGE_VELOCITY] += 1;
      fprintf(stderr,"*** Warning: Abnormally large velocity: [%lf, %lf, %lf] detected at node %d."
                     " This may be a symptom of an instability\n",V[1],V[2],V[3],locToGlobNodeMap[i] + 1);
    }

    if ((wtype == BcsWallData::WALL_FUNCTION) ||
	(wtype == BcsWallData::FULL &&
	 nodeType[i] != BC_ISOTHERMAL_WALL_MOVING &&
	 nodeType[i] != BC_ISOTHERMAL_WALL_FIXED &&
	 nodeType[i] != BC_ADIABATIC_WALL_MOVING &&
	 nodeType[i] != BC_ADIABATIC_WALL_FIXED)) {
      for (k=0; k<neq; ++k) {
	if (V[dim-neq+k] < 0.0) {
	  if (flag[i]) {
	    cmin[k]++;
	    if (V[dim-neq+k] < vmin[k]) {
	      pmin[k] = locToGlobNodeMap[i] + 1;
	      vmin[k] = V[dim-neq+k];
	    }
	  }
	  if (ctype == TsData::ABS_VALUE)
	    U[i][dim-neq+k] = fabs(U[i][dim-neq+k]);
	 else if (ctype == TsData::FREESTREAM)
	    U[i][dim-neq+k] = Uin[dim-neq+k];
	  else if (ctype == TsData::CUTOFF)
	    U[i][dim-neq+k] = 0.0;
	}
      }
    }
  }

  errorHandler->localErrors[ErrorHandler::PRESSURE_CLIPPING] += pclipping;
  errorHandler->localErrors[ErrorHandler::DENSITY_CLIPPING] += rhoclipping;
  errorHandler->localErrors[ErrorHandler::UNPHYSICAL] += ierr;
  return ierr;
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::checkFailSafe(VarFcn* varFcn, SVec<double,dim>& U,
                  SVec<bool,2>& tag, Vec<int> *fluidId)
{

  for (int i=0; i<U.size(); ++i) {
    tag[i][0] = false;
    tag[i][1] = false;
    double V[dim];
    double rho, p;
    if(!fluidId){
      varFcn->conservativeToPrimitive(U[i], V);
      rho = varFcn->getDensity(V);
      p = varFcn->checkPressure(V);
    }else{
      varFcn->conservativeToPrimitive(U[i], V, (*fluidId)[i]);
      rho = varFcn->getDensity(V,(*fluidId)[i]);
      p = varFcn->checkPressure(V, (*fluidId)[i]);
    }
    if (rho <= 0.0 || p <= 0.0)
      tag[i][0] = true;
  }

  int (*edgePtr)[2] = edges.getPtr();
  for (int l=0; l<edges.size(); ++l) {
    int i = edgePtr[l][0];
    int j = edgePtr[l][1];
    tag[i][1] = tag[i][1] || tag[j][0];
    tag[j][1] = tag[j][1] || tag[i][0];
  }

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::checkGradientsSetUp(SVec<double,3> &X, SVec<double,dim> &V)
{

  double a = 3.0, b = -2.0, c = 1.0, d = -10.0;

  for (int i=0; i<X.size(); ++i)
    for (int k=0; k<dim; ++k)
      V[i][k] = a*X[i][0] + b*X[i][1] + c*X[i][2] + d;

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::checkGradientsWrite(SVec<double,3> &X, NodalGrad<dim> &ngrad)
{

  double a = 3.0, b = -2.0, c = 1.0;

  char fileName[MAXLINE];
  sprintf(fileName, "gradients.%d", globSubNum+1);

  FILE *fp = fopen(fileName, "w");

  SVec<double,dim> &dVdx = ngrad.getX();
  SVec<double,dim> &dVdy = ngrad.getY();
  SVec<double,dim> &dVdz = ngrad.getZ();

  int j = 2;
  for (int i=0; i<dVdx.size(); ++i)
    fprintf(fp, "%d %e %e %e\n", locToGlobNodeMap[i]+1,
	    dVdx[i][j]-a, dVdy[i][j]-b, dVdz[i][j]-c);

  fclose(fp);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::checkMatVecProd(SVec<double,dim> &prod, const char *msg)
{

  Vec<int> nodeCodes(prod.size());
  nodeCodes = 0;

  int i, j;
  for (i=0; i<faces.size(); ++i)
    for (j=0; j<faces[i].numNodes(); ++j)
      nodeCodes[ faces[i][j] ] += 1;

  char fname1[MAXLINE], fname2[MAXLINE];
  sprintf(fname1, "%s.in.%d", msg, globSubNum+1);
  sprintf(fname2, "%s.bound.%d", msg, globSubNum+1);

  FILE *fp1 = fopen(fname1, "w");
  FILE *fp2 = fopen(fname2, "w");

 /*
  bool *edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int iSub = 0; iSub < numNeighb; ++iSub) {
    fprintf(fp1, "neighbor %d:\n", neighb[iSub]);
    for (int iEdge = 0; iEdge < numSharedEdges[iSub]; ++iEdge) {
      int glLeft  = sharedEdges[iSub][iEdge].glLeft;
      int glRight = sharedEdges[iSub][iEdge].glRight;
      int sign = sharedEdges[iSub][iEdge].sign;
      int num = sharedEdges[iSub][iEdge].edgeNum;
      int flag = edgeFlag[num];
      int i = locToGlobNodeMap[ edgePtr[num][0] ];
      int j = locToGlobNodeMap[ edgePtr[num][1] ];
      fprintf(fp1, "   %d (%d) %d (%d) %d %d\n", glLeft, i, glRight, j, sign, flag);
    }
  }
 */

  for (i=0; i<prod.size(); ++i) {
    FILE* fp;
    if (nodeCodes[i] == 0)
      fp = fp1;
    else
      fp = fp2;

    fprintf(fp, "%d", locToGlobNodeMap[i]+1);
    for (j=0; j<dim; ++j)
      fprintf(fp, " %e", prod[i][j]);
    fprintf(fp, "\n");
  }

  fclose(fp1);
  fclose(fp2);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeForceDerivs(VarFcn *varFcn, SVec<double,3> &X,
                                   SVec<double,dim> &V, SVec<double,dim> &deltaU,
                                   Vec<double> &modalF, SVec<double, 3> **locMX)  {

  modalF = 0.0;

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeForceDerivs(elems, varFcn, X, V, deltaU, modalF, locMX);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeForceCoefficients(PostFcn *postFcn, Vec3D &x0, GeoState &geoState,
                                         BcData<dim> &bcData, SVec<double,3> &X, SVec<double,dim> &V,
					 double pInfty, Vec3D &CFi, Vec3D &CMi, Vec3D &CFv, Vec3D &CMv,
                                         VecSet< SVec<double,3> > *mX , Vec<double> *genCF)
{

  CFi = 0.0;
  CMi = 0.0;
  CFv = 0.0;
  CMv = 0.0;

  Vec<double> &d2wall = geoState.getDistanceToWall();
  SVec<double,dim> &Vwall = bcData.getFaceStateVector();

  for (int i=0; i<faces.size(); ++i)
    faces[i].computeForceCoefficients(postFcn, x0, elems, X, V, d2wall, Vwall, pInfty,
                                      CFi, CMi, CFv, CMv, gradP, mX, genCF);

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::zeroInternalVals(SVec<double, dim> &v)  {

  for (int i = 0; i < nodes.size(); ++i)
    if (nodeType[i] == 0)
      for (int j = 0; j < dim; j++)
        v[i][j] = 0.0;
}

//------------------------------------------------------------------------------

// HB
template<int dim>
void SubDomain::zeroMeshMotionBCDofs(SVec<double,dim> &x, int* DofType)
{
  int (*dofType)[dim] = reinterpret_cast<int (*)[dim]>(DofType);
  for(int i=0;i<nodes.size(); i++)
    for(int l=0; l<dim; l++)
      if(dofType[i][l]!=BC_FREE) x[i][l] = 0.0;
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::setupUVolumesInitialConditions_Step1(const int volid, double UU[dim],
                                                     SVec<double, dim>& U,
                                                     CommPattern<double>& sp) {
  Vec<bool> flag(nodes.size());
  flag = false;

  for(int iElem = 0; iElem < elems.size(); iElem++) {
    if(elems[iElem].getVolumeID() == volid) {
      int *nodeNums = elems[iElem].nodeNum();
      for(int iNode = 0; iNode < elems[iElem].numNodes(); iNode++) {
        flag[nodeNums[iNode]] = true;
        for(int idim = 0; idim < dim; idim++) {
          U[nodeNums[iNode]][idim] = UU[idim];
        }
      }
    }
  }

  for(int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<double> sInfo = sp.getSendBuffer(sndChannel[iSub]);
    double(*buffer)[dim] = reinterpret_cast<double(*)[dim]>(sInfo.data);
    for(int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode) {
      for(int j = 0; j < dim; ++j) {
        buffer[iNode][j] = (flag[iNode]) ? U[(*sharedNodes)[iSub][iNode]][j] : -std::numeric_limits<double>::infinity();
      }
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::setupUVolumesInitialConditions_Step2(CommPattern<double>& sp,
                                                     SVec<double, dim>& U) {

  for(int iSub = 0; iSub < numNeighb; ++iSub) {
    SubRecInfo<double> sInfo = sp.recData(rcvChannel[iSub]);
    double(*buffer)[dim] = reinterpret_cast<double(*)[dim]>(sInfo.data);
    for(int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode) {
      for(int j = 0; j < dim; ++j) {
        if(buffer[iNode][j] > -std::numeric_limits<double>::infinity())
          U[(*sharedNodes)[iSub][iNode]][j] = buffer[iNode][j];
      }
    }
  }
}

//------------------------------------------------------------------------------
/*
template<int dim>
void SubDomain::setupUMultiFluidInitialConditionsSphere(FluidModelData &fm,
                       SphereData &ic, SVec<double,3> &X, SVec<double,dim> &U){

  double dist = 0.0;
  double x = ic.cen_x;
  double y = ic.cen_y;
  double z = ic.cen_z;
  double r = ic.radius;

  for (int i=0; i<U.size(); i++){
    dist = (X[i][0] - x)*(X[i][0] - x) + (X[i][1] - y)*(X[i][1] - y) + (X[i][2] - z)*(X[i][2] - z);
    if(sqrt(dist) < r) //it is inside the sphere
      for (int idim=0; idim<dim; idim++)
        U[i][idim] = UU[idim];
  }

}
*/
//------------------------------------------------------------------------------
/*
template<int dim>
void SubDomain::setupUMultiFluidInitialConditionsPlane(FluidModelData &fm,
                       PlaneData &ip, SVec<double,3> &X, SVec<double,dim> &U){


  double scalar = 0.0;
  double x = ip.cen_x;
  double y = ip.cen_y;
  double z = ip.cen_z;
  double nx = ip.nx;
  double ny = ip.ny;
  double nz = ip.nz;

  for (int i=0; i<U.size(); i++){
    scalar = nx*(X[i][0] - x)+ny*(X[i][1] - y)+nz*(X[i][2] - z);
    if(scalar > 0.0) //node is on the same side indicated by vector
      for (int idim=0; idim<dim; idim++)
        U[i][idim] = UU[idim];
  }

}
*/
//------------------------------------------------------------------------------
/*
template<int dim>
void SubDomain::setupUMultiFluidInitialConditionsPlane(FluidModelData &fm,
                       PlaneData &ip, SVec<double,3> &X, SVec<double,dim> &U, Vec<int> &nodeTag){

  double scalar = 0.0;
  double x = ip.cen_x;
  double y = ip.cen_y;
  double z = ip.cen_z;
  double nx = ip.nx;
  double ny = ip.ny;
  double nz = ip.nz;

  for (int i=0; i<U.size(); i++){
    scalar = nx*(X[i][0] - x)+ny*(X[i][1] - y)+nz*(X[i][2] - z);
    if(scalar > 0.0) {//node is on the same side indicated by vector
      nodeTag[i] = -1;
      for (int idim=0; idim<dim; idim++)
        U[i][idim] = UU[idim];
    }
  }

}
*/
//------------------------------------------------------------------------------
// TODO: should distinguish master nodes and non-master nodes
template<int dim>
void SubDomain::computeWeightsForEmbeddedStruct(SVec<double,dim> &V, SVec<double,dim> &VWeights,
																Vec<double> &Weights, LevelSetStructure &LSS, 
																SVec<double,3> &X, Vec<int> &init, Vec<int> &next_init,
	                                             bool externalSI)
{
	
  const Connectivity &nToN = *getNodeToNode();
  for(int currentNode=0;currentNode<numNodes();++currentNode)
	{		

		 // if(init[currentNode]!=1)
		 // {
		 // 	 std::cout << init[currentNode] << " " 
		 // 				  << X[currentNode][0] << " " << X[currentNode][1] << " " << X[currentNode][2] << " " 
		 // 				  << std::boolalpha << LSS.isActive(0.0,currentNode) << " " << LSS.isSwept(0.0,currentNode);

		 // 	 if(    !LSS.isActive(0.0,currentNode) && LSS.isSwept(0.0,currentNode)) std::cout << " ___ R2G";
		 // 	 else if(LSS.isActive(0.0,currentNode) && LSS.isSwept(0.0,currentNode)) std::cout << " ___ G2R";
			
		 // 	 std::cout<<"\n";
		 // }

		if(init[currentNode]<1 && LSS.isActive(0.0,currentNode))
		{		
			for(int j=0; j<nToN.num(currentNode); ++j)
			{				
        int neighborNode=nToN[currentNode][j];
        if(currentNode == neighborNode || init[neighborNode]<1) continue;
        int l = edges.findOnly(currentNode,neighborNode);

				bool intEdge;

				if(externalSI) 
					intEdge = LSS.edgeWithSI(l) || LSS.edgeIntersectsStructure(0.0, l);
				else
					intEdge = LSS.edgeIntersectsStructure(0.0, l);

				if(intEdge) 
					continue;
				else if(Weights[currentNode] < 1e-6)
				{
          Weights[currentNode]=1.0;
          next_init[currentNode]=1;

					for(int i=0; i<dim; ++i) VWeights[currentNode][i] = V[neighborNode][i];
				} 
				else 
				{
          Weights[currentNode] += 1.0;

					for(int i=0; i<dim; ++i) VWeights[currentNode][i] += V[neighborNode][i];
				}
			}
		}
	}

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeWeightsForFluidFluid(SVec<double,dim> &V, SVec<double,dim> &VWeights,
					    Vec<double> &Weights, LevelSetStructure *LSS, SVec<double,3> &X, Vec<int> &init, Vec<int> &next_init,
					    Vec<int>& fluidId)
{
  const Connectivity &nToN = *getNodeToNode();
  for(int currentNode=0;currentNode<numNodes();++currentNode)
    if(init[currentNode]<1/* && LSS.isActive(0.0,currentNode)*/){
      for(int j=0;j<nToN.num(currentNode);++j){
        int neighborNode=nToN[currentNode][j];
        if(currentNode == neighborNode || init[neighborNode]<1) continue;
        int l = edges.findOnly(currentNode,neighborNode);
        if(fluidId[currentNode] != fluidId[neighborNode]) continue;
        else if(Weights[currentNode] < 1e-6){
          Weights[currentNode]=1.0;
          next_init[currentNode]=1;
          for(int i=0;i<dim;++i)
            VWeights[currentNode][i] = V[neighborNode][i];
        } else {
          Weights[currentNode] += 1.0;
          for(int i=0;i<dim;++i)
            VWeights[currentNode][i] += V[neighborNode][i];
        }
      }
    }
}
//------------------------------------------------------------------------------
template<int dim>
void SubDomain::computeWeightsLeastSquaresForEmbeddedStruct(SVec<double,3> &X, SVec<double,10> &R, SVec<double,dim> &V, 
																				Vec<double> &Weights, SVec<double,dim> &VWeights, 
																				LevelSetStructure &LSS, Vec<int> &init, Vec<int> &next_init, 
																				NodalGrad<dim>& DX, bool limit, Vec<int>* fluidId, bool externalSI) 
{

  const Connectivity &nToN = *getNodeToNode();
  bool *masterFlag = edges.getMasterFlag();
  double lin_extrap[dim];
  for (int currentNode=0; currentNode<numNodes(); ++currentNode)
	{
		if(init[currentNode]<1 && LSS.isActive(0.0,currentNode)) 
		{
			for(int j=0; j<nToN.num(currentNode); ++j) 
			{
		int neighborNode = nToN[currentNode][j];
		if (currentNode==neighborNode || init[neighborNode]<1) continue;
		int l = edges.findOnly(currentNode,neighborNode);
		if (!masterFlag[l]) continue;

				bool intEdge;

				if(externalSI) 
					intEdge = LSS.edgeWithSI(l) || LSS.edgeIntersectsStructure(0.0, l);
				else
					intEdge = LSS.edgeIntersectsStructure(0.0, l);

				if(intEdge) continue;

		double dx[3] = {X[neighborNode][0]-X[currentNode][0],
				X[neighborNode][1]-X[currentNode][1],
				X[neighborNode][2]-X[currentNode][2]};

				/*
				  if (fabs(Weights[currentNode])<1e-6) {
		  next_init[currentNode] = 1;
		  if (R[currentNode][0]>0.0) {
		    double dx[3] = {X[neighborNode][0]-X[currentNode][0],
		    				X[neighborNode][1]-X[currentNode][1],
		    				X[neighborNode][2]-X[currentNode][2]};
		    double W[4];
		    Weights[currentNode] = -1.0;
		    computeLocalWeightsLeastSquaresForEmbeddedStruct(dx,R[currentNode],W);
		    for (int k=0; k<dim; ++k) VWeights[currentNode][k] = W[3]*V[neighborNode][k];
		  } else {
			Weights[currentNode] = 1.0;
		    for (int k=0; k<dim; ++k) VWeights[currentNode][k] = V[neighborNode][k];
		  }
		} else {
		  if (R[currentNode][0]>0.0) {
		    double dx[3] = {X[neighborNode][0]-X[currentNode][0],
		    				X[neighborNode][1]-X[currentNode][1],
		    				X[neighborNode][2]-X[currentNode][2]};
		    double W[4];
		    computeLocalWeightsLeastSquaresForEmbeddedStruct(dx,R[currentNode],W);
		    for (int k=0; k<dim; ++k) VWeights[currentNode][k] += W[3]*V[neighborNode][k];
		  }
		  else {
			Weights[currentNode] += 1.0;
		    for (int k=0; k<dim; ++k) VWeights[currentNode][k] += V[neighborNode][k];
		  }
		  }
		*/
		//if (fluidId && fluidId[currentNode] != fluidId[neighborNode]) continue;
				next_init[currentNode] = 1;
				double weight = 1.0;/*std::max<double>(1e-8,-V[neighborNode][1]*dx[0]-
								 V[neighborNode][2]*dx[1]-
								 V[neighborNode][3]*dx[2]);*/
		Weights[currentNode] += weight;

				for(int k=0; k<dim; ++k) 
				{
					lin_extrap[k] = V[neighborNode][k] 
						           - DX.getX()[neighborNode][k]*dx[0]
						           - DX.getY()[neighborNode][k]*dx[1]
						           - DX.getZ()[neighborNode][k]*dx[2];
		}

		double alpha = 1.0;
				if(limit) alpha = higherOrderFSI->computeAlpha<dim>(neighborNode, V[neighborNode], lin_extrap);
		
				for(int k=0; k<dim; ++k) 
				{
					VWeights[currentNode][k] += weight * lin_extrap[k]*(alpha) + V[neighborNode][k]*(1.0-alpha);
		  // std::cout << currentNode << " " << lin_extrap[k] << " " <<  V[neighborNode][k] << std::endl;
		}
			}
	  }
    }
}

template<int dim>
void SubDomain::computeWeightsLeastSquaresForFluidFluid(
		SVec<double,3> &X, SVec<double,10> &R, SVec<double,dim> &V, Vec<double> &Weights, 
		SVec<double,dim> &VWeights, LevelSetStructure *LSS, Vec<int> &init, Vec<int> &next_init,Vec<int>& fluidId,NodalGrad<dim>& DX,bool limit) 
{
  const Connectivity &nToN = *getNodeToNode();
  bool *masterFlag = edges.getMasterFlag();
  double lin_extrap[dim];
  for (int currentNode=0; currentNode<numNodes(); ++currentNode)
    if (init[currentNode]<1/* && LSS.isActive(0.0,currentNode)*/) {
	  for (int j=0; j<nToN.num(currentNode); ++j) {
		int neighborNode = nToN[currentNode][j];
		if (currentNode==neighborNode || init[neighborNode]<1) continue;
		int l = edges.findOnly(currentNode,neighborNode);
	        if (!masterFlag[l]) continue;
		if (fluidId[currentNode] != fluidId[neighborNode]) continue;
		double dx[3] = {X[neighborNode][0]-X[currentNode][0],
				X[neighborNode][1]-X[currentNode][1],
				X[neighborNode][2]-X[currentNode][2]};

		next_init[currentNode] = 1;
		/*else if (fabs(Weights[currentNode])<1e-6) {
		  next_init[currentNode] = 1;
		  if (R[currentNode][0]>0.0) {
		    double dx[3] = {X[neighborNode][0]-X[currentNode][0],
		    				X[neighborNode][1]-X[currentNode][1],
		    				X[neighborNode][2]-X[currentNode][2]};
		    double W[4];
		    Weights[currentNode] = -1.0;
		    computeLocalWeightsLeastSquaresForEmbeddedStruct(dx,R[currentNode],W);
		    //double alpha = higherOrderMF->computeAlpha<dim>(neighborNode,V[neighborNode],
								    
		    for (int k=0; k<dim; ++k) VWeights[currentNode][k] = W[3]*V[neighborNode][k];
		  } else {
			Weights[currentNode] = 1.0;
		    for (int k=0; k<dim; ++k) VWeights[currentNode][k] = V[neighborNode][k];
		  }
		} else {
		  if (R[currentNode][0]>0.0) {
		    double dx[3] = {X[neighborNode][0]-X[currentNode][0],
		    				X[neighborNode][1]-X[currentNode][1],
		    				X[neighborNode][2]-X[currentNode][2]};
		    double W[4];
		    computeLocalWeightsLeastSquaresForEmbeddedStruct(dx,R[currentNode],W);
		    for (int k=0; k<dim; ++k) VWeights[currentNode][k] += W[3]*V[neighborNode][k];
		  }
		  else {
			Weights[currentNode] += 1.0;
		    for (int k=0; k<dim; ++k) VWeights[currentNode][k] += V[neighborNode][k];
		  }
		  }*/

		/*std::cout << "Computing weights using linear extrapolation, limit = " << limit << std::endl;
		std::cout << "dx = [ " << DX.getX()[neighborNode][0]*dx[0] << " " << 
		  DX.getY()[neighborNode][0]*dx[1] << " "  << DX.getZ()[neighborNode][0]*dx[2] << "]\n";
		std::cout << "V[neighbor_node] = " << V[neighborNode][0] << std::endl;
		std::cout << "curr_node = " << currentNode << std::endl;*/
		Weights[currentNode] += 1.0;
		for (int k=0; k<dim; ++k) {

		  lin_extrap[k] = V[neighborNode][k]-DX.getX()[neighborNode][k]*dx[0]-
		    DX.getY()[neighborNode][k]*dx[1]-
		    DX.getZ()[neighborNode][k]*dx[2];
		}

		double alpha = 1.0;
		if (limit)
		  alpha = higherOrderMF->computeAlpha<dim>(neighborNode,V[neighborNode],
							   lin_extrap);
		
		//std::cout << "alpha = " << alpha << std::endl;
		for (int k=0; k<dim; ++k)
		  VWeights[currentNode][k] += lin_extrap[k]*(alpha)+
		    V[neighborNode][k]*(1.0-alpha);
		
	  }
    }
}


//------------------------------------------------------------------------------
// who: active and "swept" nodes
// what: phi and U
// how: pull from neighbours which are visible, not swept, and have the same fluid Id.
// TODO: more efficent in an "edge loop".
template<int dim, int dimLS>
void SubDomain::computeWeightsForEmbeddedStruct(SVec<double,dim> &V, SVec<double,dim> &VWeights,
                                                SVec<double,dimLS> &Phi, SVec<double,dimLS> &PhiWeights, 
                                                Vec<double> &Weights, LevelSetStructure &LSS, SVec<double,3> &X,
                                                Vec<int> &init, Vec<int> &next_init, Vec<int> &fluidId)
{
  const Connectivity &nToN = *getNodeToNode();
  for(int currentNode=0;currentNode<numNodes();++currentNode) {

    int caught = 0;
    if(locToGlobNodeMap[currentNode]+1 == 72844)
      caught = 1;

    if(init[currentNode]<1.0 && !LSS.isOccluded(0.0,currentNode)){
      int myId = fluidId[currentNode]; 
      for(int j=0;j<nToN.num(currentNode);++j){
        int neighborNode=nToN[currentNode][j];
        int yourId = fluidId[neighborNode];

      if(caught && currentNode!=neighborNode)
        fprintf(stderr,"SubDomain->computeWeights..., Nei of 72844: %d, init = %d, fluidId = %d, occluded = %d, swept = %d, X = %d\n",
            locToGlobNodeMap[neighborNode]+1, init[neighborNode], fluidId[neighborNode], LSS.isOccluded(0.0,neighborNode), 
            LSS.isSwept(0.0,neighborNode), (int)LSS.edgeIntersectsStructure(0.0,edges.findOnly(currentNode,neighborNode)));

        if(currentNode==neighborNode || init[neighborNode]<1 || myId!=yourId) continue;
        int l = edges.findOnly(currentNode,neighborNode);
        if(LSS.edgeIntersectsStructure(0.0,l)) continue;
        else if(Weights[currentNode] < 1e-6){
          Weights[currentNode]=1.0;
          next_init[currentNode]=1;
          for(int i=0;i<dim;++i)
            VWeights[currentNode][i] = V[neighborNode][i];
          for(int i=0;i<dimLS;++i)
            PhiWeights[currentNode][i] = Phi[neighborNode][i];
        } else {
          Weights[currentNode] += 1.0;
          for(int i=0;i<dim;++i)
            VWeights[currentNode][i] += V[neighborNode][i];
          for(int i=0;i<dimLS;++i)
            PhiWeights[currentNode][i] += Phi[neighborNode][i];
        }
      }
    } else if (init[currentNode]<1.0) {//occluded node, pull (1) velocity, (2) phi from every neighbor.
      int myId = fluidId[currentNode];
      for(int j=0;j<nToN.num(currentNode);++j){
        int neighborNode=nToN[currentNode][j];
        int yourId = fluidId[neighborNode];
        if(currentNode==neighborNode || init[neighborNode]<1) continue;
        int l = edges.findOnly(currentNode,neighborNode);
        if(Weights[currentNode] < 1e-6){
          Weights[currentNode]=1.0;
          next_init[currentNode]=1;
          for(int i=1;i<4;++i)  //only get velocity
            VWeights[currentNode][i] = V[neighborNode][i];
          for(int i=0;i<dimLS;++i)
            PhiWeights[currentNode][i] = Phi[neighborNode][i];
        } else {
          Weights[currentNode] += 1.0;
          for(int i=1;i<4;++i)  //only get velocity
            VWeights[currentNode][i] += V[neighborNode][i];
          for(int i=0;i<dimLS;++i)
            PhiWeights[currentNode][i] += Phi[neighborNode][i];
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

template<int dimLS>
void SubDomain::extrapolatePhiV(LevelSetStructure &LSS, SVec<double,dimLS> &PhiV)
{
  int (*edgePtr)[2] = edges.getPtr();
  bool *masterFlag = edges.getMasterFlag();
  
  for(int l=0; l<edges.size(); l++) {
    if(!masterFlag[l])
      continue;
    int i = edgePtr[l][0];
    int j = edgePtr[l][1];

    if(!LSS.withCracking()) {
      int Idi = LSS.fluidModel(0.0,i), Idj = LSS.fluidModel(0.0,j);
      if(Idi!=0 || Idj!=0) //meaning one of them is isolated by structure
        continue;
    }

    bool iSwept = LSS.isSwept(0.0,i), jSwept = LSS.isSwept(0.0,j);
 
    // pull data from j to i?
    if(iSwept && !jSwept)
      for(int k=0; k<dimLS; k++)
        PhiV[i][k] += PhiV[j][k];

    // pull data from i to j?
    if(jSwept && !iSwept)
      for(int k=0; k<dimLS; k++)
        PhiV[j][k] += PhiV[i][k];
  }
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::populateGhostPoints(Vec<GhostPoint<dim>*> &ghostPoints, SVec<double,3> &X, 
												SVec<double,dim> &U, NodalGrad<dim, double> &ngrad, 
												VarFcn *varFcn,LevelSetStructure &LSS,bool linRecFSI,Vec<int> &tag)
{

  int i, j, k;
  double alpha;

  bool* edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  double *Vi,*Vj,*weights;
  Vi = new double[dim];
  Vj = new double[dim];
  weights = new double[dim];

  for (int l=0; l<edges.size(); l++) {
    if(!edgeFlag[l]) continue; //not a master edge
    i = edgePtr[l][0];
    j = edgePtr[l][1];
    if(LSS.edgeIntersectsStructure(0.0,l)) { // at interface
      int tagI = tag[i];
      int tagJ = tag[j];
      bool iIsActive = LSS.isActive(0.0,i);
      bool jIsActive = LSS.isActive(0.0,j);

      if(iIsActive) {
        LevelSetResult resij = LSS.getLevelSetDataAtEdgeCenter(0.0, l, true);
        varFcn->conservativeToPrimitive(U[i],Vi,tagI);

// Initialize all variables and weights
        for (int k=0; k<dim; ++k) {
          Vj[k] = Vi[k];
	  weights[k] = 1.0;
        }

// Replace fourth variable with temperature
        double T = varFcn->computeTemperature(Vi,tagI);
        Vj[4] = T;

// Determine intersection alpha 
        if (!linRecFSI) { //first order
          alpha = 0.5;
        }
        else {
          alpha = resij.alpha;
          if (alpha > 0.5) alpha = 0.5; // Set limit for stability 
        }

// Update velocity
        for (int k=1; k<4; ++k) {
	  Vj[k] = ((resij.normVel)[k-1] - alpha*Vi[k])/(1.0-alpha);
          weights[k] = (1.0-alpha)*(1.0-alpha);
        }

	if (dim==6) {  // One Equation Turbulent Model
	  Vj[5] = 0.0;//-alpha*Vi[5]/(1.0-alpha);
          weights[5] = 1.0;//(1.0-alpha)*(1.0-alpha);
	}
	else if (dim==7) { // Two Equations Turbulent Model
	  Vj[5] = -alpha*Vi[5]/(1.0-alpha);
	  Vj[6] = -alpha*Vi[6]/(1.0-alpha);
          weights[5] = (1.0-alpha)*(1.0-alpha);
          weights[6] = (1.0-alpha)*(1.0-alpha);
	}

        if(!ghostPoints[j]) // GP has not been created
        {ghostPoints[j]=new GhostPoint<dim>(varFcn);}

        ghostPoints[j]->addNeighbour(Vj,weights,tagI);
      }
      if(jIsActive) {
        LevelSetResult resji = LSS.getLevelSetDataAtEdgeCenter(0.0, l, false);
        varFcn->conservativeToPrimitive(U[j],Vj,tagJ);

// Initialize all variables and weights
        for (int k=0; k<dim; ++k) {
          Vi[k] = Vj[k];
	  weights[k] = 1.0;
        }

// Replace fourth variable with temperature
        double T = varFcn->computeTemperature(Vj,tagJ);
        Vi[4] = T;

// Determine intersection alpha 
        if (!linRecFSI) { //first order
          alpha = 0.5;
        }
        else {
          alpha = resji.alpha;
          if (alpha > 0.5) alpha = 0.5; // Set limit for stability 
        }

// Update velocity
        for (int k=1; k<4; ++k) {
	  Vi[k] = ((resji.normVel)[k-1] - alpha*Vj[k])/(1.0-alpha);
          weights[k] = (1.0-alpha)*(1.0-alpha);
        }

	if (dim==6) {  // One Equation Turbulent Model
	  Vi[5] = 0.0;//-alpha*Vj[5]/(1.0-alpha);
          weights[5] = 1.0;//(1.0-alpha)*(1.0-alpha);
	}
	else if (dim==7) { // Two Equations Turbulent Model
	  Vi[5] = -alpha*Vj[5]/(1.0-alpha);
	  Vi[6] = -alpha*Vj[6]/(1.0-alpha);
          weights[5] = (1.0-alpha)*(1.0-alpha);
          weights[6] = (1.0-alpha)*(1.0-alpha);
	}

        if(!ghostPoints[i]) // GP has not been created
        {ghostPoints[i]=new GhostPoint<dim>(varFcn);}

        ghostPoints[i]->addNeighbour(Vi,weights,tagJ);
      }
    }
  }
  delete[] Vi;
  delete[] Vj;
  delete[] weights;

}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::reduceGhostPoints(Vec<GhostPoint<dim>*> &ghostPoints, SVec<double,3> &X)
{	
	bool isOK;
  for (int i=0; i<nodes.size(); i++)
      if(ghostPoints[i]) {
			//std::cout << X[i][0] << " "<< X[i][1] << " "<< X[i][2] << " ";
	  ghostPoints[i]->reduce();
	}
    } 

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void SubDomain::populateGhostJacobian(Vec<GhostPoint<dim>*> &ghostPoints,SVec<double,dim> &U,
												  FluxFcn** fluxFcn, VarFcn *varFcn, LevelSetStructure &LSS, 
												  Vec<int> &tag, GenMat<Scalar,neq>& A)
{

  int i, j, k;

  bool* edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  double B[neq*neq],tmp[neq*neq];
  double dUdV[neq*neq],dVdU[neq*neq];
  memset(tmp,0,sizeof(double)*neq*neq);
  memset(dUdV,0,sizeof(double)*neq*neq);
  memset(dVdU,0,sizeof(double)*neq*neq);
  memset(B,0,sizeof(double)*neq*neq);

	if(neq > 2)  
	{
    for (k = 1; k < 4; ++k) B[k*neq+k] = 1.0;
    B[0] = B[4*neq+4] = -1.0;
  }

	for(int l=0; l<edges.size(); l++) 
	{
    i = edgePtr[l][0];
    j = edgePtr[l][1];

		if(LSS.edgeIntersectsStructure(0.0, l)) 
		{ 
         //at interface
      int tagI = tag[i];
      int tagJ = tag[j];
      bool iIsActive = LSS.isActive(0.0,i);
      bool jIsActive = LSS.isActive(0.0,j);
	  
			if(iIsActive) 
			{
        double (*Aji)[neq*neq] = reinterpret_cast< double (*)[neq*neq]>(A.getGhostNodeElem_ij(j,i));

				if(Aji) 
				{
					if(neq > 2) 
					{
            Vec<double> Vi(dim);
            varFcn->conservativeToPrimitive(U[i],Vi.v,tagI);
            fluxFcn[BC_INTERNAL]->getFluxFcnBase(tagI)->getVarFcnBase()->computedVdU(Vi.v,dVdU);
            fluxFcn[BC_INTERNAL]->getFluxFcnBase(tagI)->getVarFcnBase()->computedUdV(ghostPoints[j]->getPrimitiveState(), dUdV);
            DenseMatrixOp<double,neq,neq*neq>::applyToDenseMatrix(&dUdV, 0, &B, 0, &tmp, 0);
            DenseMatrixOp<double,neq,neq*neq>::applyToDenseMatrix(&tmp, 0, &dVdU, 0, Aji, 0);
          }
					else 
					{
            for (k = 0; k < neq; ++k) tmp[k*neq+k] = 1.0;
            DenseMatrixOp<double,neq,neq*neq>::applyToDenseMatrix(&tmp, 0, &B, 0, Aji, 0);
          }
	}

        double (*Ajj) = reinterpret_cast< double (*)>(A.getGhostGhostElem_ij(j,j));

				if(Ajj) 
				{
          memset(Ajj,0,sizeof(double)*neq*neq);
          for (k = 0; k < neq; ++k) Ajj[k*neq+k] = ghostPoints[j]->ng;
	}
      }

			if(jIsActive) 
			{
        double (*Aij)[neq*neq] = reinterpret_cast< double (*)[neq*neq]>(A.getGhostNodeElem_ij(i,j));

				if(Aij) 
				{
					if(neq > 2) 
					{
            Vec<double> Vj(dim);
            varFcn->conservativeToPrimitive(U[j],Vj.v,tagJ);
            fluxFcn[BC_INTERNAL]->getFluxFcnBase(tagJ)->getVarFcnBase()->computedVdU(Vj.v,dVdU);
            fluxFcn[BC_INTERNAL]->getFluxFcnBase(tagJ)->getVarFcnBase()->computedUdV(ghostPoints[i]->getPrimitiveState(), dUdV);
            DenseMatrixOp<double,neq,neq*neq>::applyToDenseMatrix(&dUdV, 0, &B, 0, &tmp, 0);
            DenseMatrixOp<double,neq,neq*neq>::applyToDenseMatrix(&tmp, 0, &dVdU, 0, Aij, 0);
	    }
					else
					{
            for (k = 0; k < neq; ++k) tmp[k*neq+k] = 1.0;
            DenseMatrixOp<double,neq,neq*neq>::applyToDenseMatrix(&tmp, 0, &B, 0, Aij, 0);
          }
	}

        double (*Aii) = reinterpret_cast< double (*)>(A.getGhostGhostElem_ij(i,i));

				if(Aii) 
				{
          memset(Aii,0,sizeof(double)*neq*neq);
          for (k = 0; k < neq; ++k) Aii[k*neq+k] = ghostPoints[i]->ng;
	}
      }
    }
  }

}

//------------------------------------------------------------------------------
//d2d
template<int dim>
void SubDomain::populateGhostPoints(Vec<GhostPoint<dim>*> &ghostPoints, SVec<double,3> &X, 
												SVec<double,dim> &U, NodalGrad<dim, double> &ngrad, 
												VarFcn *varFcn, LevelSetStructure &LSS, Vec<int> &fluidId, 
												FemEquationTerm *fet)
{

	double *Vg1_i = new double[dim];
	double *Vg2_i = new double[dim];

	double *weights = new double[dim];

	for(int k=0; k<dim; ++k) weights[k] = 1.0;

	bool wRcn, isIsoTherm, gotIt, gotIt1, gotIt2;
	double TWall;

	int fId1, fId2;

	Vec3D xWall, nWall, vWall;

	for(int i=0; i<nodes.size(); ++i)
	{
		wRcn = LSS.xWallNode(i, xWall);

		if(wRcn) 
		{
			if(!ghostPoints[i]) ghostPoints[i] = new GhostPoint<dim>(varFcn);

			bool dummy = LSS.vWallNode(i, vWall);

			isIsoTherm = LSS.getTwall(TWall);

			LSS.nWallNode(i, nWall);

			gotIt1 = higherOrderFSI->setFEGhostPoint(1, i, varFcn, U, ngrad, X, fluidId,
																  xWall, vWall, nWall, isIsoTherm, TWall,
																  fet, Vg1_i, fId1);

			gotIt2 = higherOrderFSI->setFEGhostPoint(-1, i, varFcn, U, ngrad, X, fluidId,
																  xWall, vWall, nWall, isIsoTherm, TWall,
																  fet, Vg2_i, fId2);

			if( fabs(X[i][0]-0.974286)<1.0e-4 && fabs(X[i][1]-0.005102) < 1.0e-4 && fabs(X[i][2]) <1.0e-5) 

			if(gotIt1 || gotIt2) ghostPoints[i]->addNeighbour(gotIt1, Vg1_i, fluidId[i], 
																			  gotIt2, Vg2_i, fluidId[i], weights);			
		}
	}

	delete [] Vg1_i;
	delete [] Vg2_i;
	delete [] weights;

}
//--------------------------------------------------------------------------

template<int dim>
void SubDomain::checkGhostPoints(Vec<GhostPoint<dim>*> &ghostPoints, SVec<double,3> &X, 
											SVec<double,dim> &U, NodalGrad<dim, double> &ngrad, 
											VarFcn *varFcn, LevelSetStructure &LSS, Vec<int> &fluidId)
{

	Vec3D xWall, nWall, vWall;
	bool isIsoTherm, gValid;
	double TWall;

	double *Vf = new double[dim];
	int tag;

	double *V_;
	double *W_;

	bool dummy;

	for(int i=0; i<nodes.size(); ++i)
	{
		if(ghostPoints[i])
		{
			gValid = ghostPoints[i]->getStatus();

			V_ = ghostPoints[i]->V;
			W_ = ghostPoints[i]->Ws;

			if(gValid) continue;			

			dummy = LSS.xWallNode(i, xWall);
			dummy = LSS.vWallNode(i, vWall);

			isIsoTherm = LSS.getTwall(TWall);

			LSS.nWallNode(i, nWall);

			Vec3D Xi = X[i];
				
			int Nsize;
			int* Nlist;
			Nlist = getNeiNodeOfNode(i, Nsize);
			
			double minDist = FLT_MAX;
			
			bool gotIt = false;

			for(int j=0; j<Nsize; ++j) 
			{
				int Nj = Nlist[j];

				if(!LSS.isActive(0.0, Nj)) continue;

				Vec3D Xj = X[Nj];
					
				Vec3D vdist = Xi - Xj;
					
				double dist = sqrt(vdist*vdist);	

				if(dist < minDist)
				{
					minDist = dist;

					tag = fluidId[Nj];

					varFcn->conservativeToPrimitive(U[Nj], Vf, fluidId[Nj]);

					gotIt = true;
				}
			}

			if(!gotIt) 
			{
				fprintf(stderr, " *** error check ghost point \n");
				exit(-1);
			}

			for(int k=1; k<4; ++k) Vf[k] = vWall[k-1];
			
			if(isIsoTherm) Vf[4] = TWall;
			else  
			{
				double T = varFcn->computeTemperature(Vf, tag);
				Vf[4] = T;
			}

			if(dim>5) for(int k=5; k<dim; ++k) Vf[k] = 0.0;

			ghostPoints[i]->set(Vf, tag);
		}
	}	

	delete[] Vf;

}
//--------------------------------------------------------------------------

//d2d
template<int dim>
void SubDomain::setSIstencil(SVec<double,3> &X, LevelSetStructure &LSS, 
									  Vec<int> &fluidId, SVec<double,dim> &U)
{

	int i, j;

	bool* edgeFlag    = edges.getMasterFlag();
	int (*edgePtr)[2] = edges.getPtr();

	Vec3D xWall, normWall;

	V6NodeData (*SiStencilData);
   SiStencilData = 0;

	if(!SiStencilData) 
		SiStencilData = new V6NodeData[edges.size()]; 
	
	bool withSI = false;

	for(int l=0; l<edges.size(); l++)
	{
		//if( edgeFlag[l] ) continue;???

		if(!LSS.edgeWithSI(l)) continue;

		i = edgePtr[l][0];
		j = edgePtr[l][1];		

		Vec3D xwall, nwall;
	   LSS.xWallWithSI(l, xwall);
	   LSS.nWallWithSI(l, nwall);

		bool gotIt = getSIstencil(i, j, X, LSS, fluidId, nwall, xwall, SiStencilData[l]);
	  
		withSI = withSI || gotIt;
	}

	if(withSI) 
		higherOrderFSI->setSIstencil(SiStencilData, U);
	else
		delete [] SiStencilData;

}

//--------------------------------------------------------------------------
//d2d
template<int dim>
void SubDomain::setFEMstencil(SVec<double,3> &X, LevelSetStructure &LSS, 
										Vec<int> &fluidId, SVec<double,dim> &U)
{
	
	Vec3D xWall, normWall;

	V6NodeData (*NodeStencilData_p);
   NodeStencilData_p = 0;

	V6NodeData (*NodeStencilData_m);
   NodeStencilData_m = 0;
	
	if(!NodeStencilData_p) NodeStencilData_p = new V6NodeData[nodes.size()]; 
	if(!NodeStencilData_m) NodeStencilData_m = new V6NodeData[nodes.size()]; 
	
	bool withGhost = false;

	for(int i=0; i<nodes.size(); i++)
	{
		Vec3D xwall, nwall;

	   bool isIGhost = LSS.xWallNode(i, xwall);

		if(!isIGhost) continue;

	   LSS.nWallNode(i, nwall);

		//bool gotIt = getFEMstencil(i, X, LSS, fluidId, nwall, xwall, NodeStencilData[i]);
		bool gotIt = getFEMstencil2(i, X, LSS, fluidId, nwall, xwall, 
											 NodeStencilData_p[i], NodeStencilData_m[i]);

		withGhost = withGhost || gotIt;
	}

	if(withGhost)
	 	higherOrderFSI->setFEMstencil(NodeStencilData_p, NodeStencilData_m, U);
	else
	{
	 	delete [] NodeStencilData_p;
	 	delete [] NodeStencilData_m;
	}

}

//--------------------------------------------------------------------------
//TODO: should distinguish master edges and non-master edges.
template<int dim>
void SubDomain::computeRiemannWeightsForEmbeddedStruct(SVec<double,dim> &V, SVec<double,dim> &Wstarij,
                      SVec<double,dim> &Wstarji, SVec<double,dim> &VWeights, Vec<double> &Weights, 
                      LevelSetStructure &LSS, SVec<double,3> &X)
{

  int i, j, k;

  bool* edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); l++){
    i = edgePtr[l][0];
    j = edgePtr[l][1];

    if(LSS.edgeIntersectsStructure(0.0,l)){ //at interface
      if (LSS.isActive(0.0,i)) {// add Wstarij on node j.
        if(Weights[j]<1.e-6) {
          Weights[j] = 1.0;
          if(Wstarij[l][0]>1.0e-8) //use Wstarij.
            for(k=0; k<dim; k++)
              VWeights[j][k] = Wstarij[l][k];
          else { 
//            fprintf(stderr,"Storing weights for Node %d: have to use nodal state for edge (%d->%d) .\n", locToGlobNodeMap[j]+1, locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1);
            for(k=0; k<dim; k++)
              VWeights[j][k] = V[i][k];
          }
        }else{
          Weights[j] += 1.0;
          if(Wstarij[l][0]>1.0e-8)//use Wstarij
            for(k=0; k<dim; k++)
              VWeights[j][k] += Wstarij[l][k];
          else {
//            fprintf(stderr,"Storing weights for Node %d: have to use nodal state for edge (%d->%d) .\n", locToGlobNodeMap[j]+1, locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1);
            for(k=0; k<dim; k++)
              VWeights[j][k] += V[i][k];
          }
        }
      }
 
      if (LSS.isActive(0.0,j)) {// add Wstarji on node i
        if(Weights[i]<1.e-6) {
          Weights[i] = 1.0;
          if(Wstarji[l][0]>1.0e-8) //use Wstarji.
            for(k=0; k<dim; k++)
              VWeights[i][k] = Wstarji[l][k];
          else {
//            fprintf(stderr,"Storing weights for Node %d: have to use nodal state for edge (%d->%d) .\n", locToGlobNodeMap[i]+1, locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1);
            for(k=0; k<dim; k++)
              VWeights[i][k] = V[j][k];
          }
        }else{
          Weights[i] += 1.0;
          if(Wstarji[l][0]>1.0e-8)//use Wstarji.
            for(k=0; k<dim; k++)
              VWeights[i][k] += Wstarji[l][k];
          else {
//            fprintf(stderr,"Storing weights for Node %d: have to use nodal state for edge (%d->%d) .\n", locToGlobNodeMap[i]+1, locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1);
            for(k=0; k<dim; k++)
              VWeights[i][k] += V[j][k];
          }
        }
      }
 
    }

  }
}

//--------------------------------------------------------------------------
//TODO: should distinguish master edges and non-master edges.
template<int dim, int dimLS>
void SubDomain::computeRiemannWeightsForEmbeddedStruct(SVec<double,dim> &V, SVec<double,dim> &Wstarij,
                      SVec<double,dim> &Wstarji, SVec<double,dim> &VWeights, Vec<double> &Weights,
                      SVec<double,dimLS> &Phi, SVec<double,dimLS> &PhiWeights,  
                      LevelSetStructure &LSS, SVec<double,3> &X, Vec<int> &fluidId0, Vec<int> &fluidId)
{
  int i, j, k;
  bool* edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); l++){
    for (int enode = 0; enode<2; enode++) { //loop thru the two nodes on this edge
      SVec<double,dim> &Wstar = (enode==0) ? Wstarji : Wstarij;
      i = edgePtr[l][enode];
      j = edgePtr[l][(int)(!enode)];
      // for V and Phi
      if(LSS.isSwept(0.0,i) && !LSS.isOccluded(0.0,i)) { // phase change occurred && need an update
        if(Wstar[l][0]>1.0e-8 && fluidId0[j]==fluidId[i]) { //use Wstar 
          if(Weights[i]<1.0e-6) { // first touch of node i
            Weights[i] = 1.0;
            for(k=0; k<dim; k++) VWeights[i][k] = Wstar[l][k];
            for(k=0; k<dimLS; k++) PhiWeights[i][k] = Phi[j][k];
          } else {
            Weights[i] ++;
            for(k=0; k<dim; k++) VWeights[i][k] += Wstar[l][k];
            for(k=0; k<dimLS; k++) PhiWeights[i][k] += Phi[j][k];
          }
        } else if(!LSS.isSwept(0.0,j) && fluidId[i]==fluidId[j] && !LSS.edgeIntersectsStructure(0.0,l)) { // use V[j]
          if(Weights[i]<1.0e-6) { // first touch of node i
            Weights[i] = 1.0;
            for(k=0; k<dim; k++) VWeights[i][k] = V[j][k];
            for(k=0; k<dimLS; k++) PhiWeights[i][k] = Phi[j][k];
          } else {
            Weights[i] ++;
            for(k=0; k<dim; k++) VWeights[i][k] += V[j][k];
            for(k=0; k<dimLS; k++) PhiWeights[i][k] += Phi[j][k];
          }
        }
      }
    }
  }
}

//--------------------------------------------------------------------------
template<int dimLS>
void SubDomain::checkNodePhaseChange(SVec<double,dimLS> &PhiProduct)
{

  for(int i=0; i<nodes.size(); i++)
    for(int idim=0; idim<dimLS; idim++)
      if(PhiProduct[i][idim]<0.0)
        fprintf(stdout, "***Error: node %d (%d) has changed phase during reinitialization\n", i, locToGlobNodeMap[i]+1);

}
//------------------------------------------------------------------------------
template<int dim>
void SubDomain::storePreviousPrimitive(SVec<double,dim> &V, Vec<int> &fluidId, 
                                    SVec<double,3> &X, SVec<double,dim> &Vupdate, 
                                    Vec<double> &weight){

  int i, j, k;

  bool* edgeFlag = edges.getMasterFlag();
  int (*edgePtr)[2] = edges.getPtr();

  for (int l=0; l<edges.size(); l++){
    i = edgePtr[l][0];
    j = edgePtr[l][1];

// selective averaged extrapolation (based on direction of the flow)

    if(fluidId[i] != fluidId[j]){ //at interface
      double dx[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
      double normdx2 = dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
      double normUi2 = V[i][1]*V[i][1]+V[i][2]*V[i][2]+V[i][3]*V[i][3];
      double normUj2 = V[j][1]*V[j][1]+V[j][2]*V[j][2]+V[j][3]*V[j][3];
      double udotdx = 0.0;

      if(normdx2*normUj2 > 0.0)
        udotdx = -(dx[0]*V[j][1]+dx[1]*V[j][2]+dx[2]*V[j][3])/sqrt(normdx2*normUj2);
      if(udotdx > 0.0){
        weight[i] += udotdx;
        for(k=0; k<dim; k++)
          Vupdate[i][k] += udotdx*V[j][k];
      }

      if(normdx2*normUi2 > 0.0)
        udotdx = (dx[0]*V[i][1]+dx[1]*V[i][2]+dx[2]*V[i][3])/sqrt(normdx2*normUi2);
      if(udotdx > 0.0){
        weight[j] += udotdx;
        for(k=0; k<dim; k++)
          Vupdate[j][k] += udotdx*V[i][k];
      }
    }

  }

}
//------------------------------------------------------------------------------
template<int dim>
void SubDomain::IncreasePressure(double p, VarFcn *vf, SVec<double,dim> &U){
// only the pressure in the volumes that have an id=0 (outside part of a cylinder phi>0)
// are updated. It is assumed that the input and output states are uniform!

  bool found = false;
  double Uc[dim];

  for(int i=0; i<nodes.size(); i++) {
    if(!found) {//fill Uc[dim]
      double V[dim];
      vf->conservativeToPrimitive(U[i],V);
      vf->setPressure(p,V); //for Tait it actually modifies density.
      vf->primitiveToConservative(V,Uc);
      found = true;
    }
    for(int k=0; k<dim; k++)
      U[i][k] = Uc[k];
  }

/*  double V[dim];
  vf->conservativeToPrimitive(U[0],V);
  vf->setPressure(p,V);
  double rhoe = vf->computeRhoEnergy(V);


  for (int iElem = 0; iElem < elems.size(); iElem++)  {
    if (elems[iElem].getVolumeID() == 0)  {
      int *nodeNums = elems[iElem].nodeNum();
      for (int iNode = 0; iNode < elems[iElem].numNodes(); iNode++)
        U[nodeNums[iNode]][4] = rhoe;
    }
  }
*/
}
//------------------------------------------------------------------------------
template<int dim>
void SubDomain::IncreasePressure(double p, VarFcn *vf, SVec<double,dim> &U, Vec<int> &fluidId){
// only the pressure with fluidId = 0 (outside part)
// are updated. It is assumed that the input and output states are uniform!

  bool found = false;
  double Uc[dim];

  for(int i=0; i<nodes.size(); i++) {
    if(fluidId[i]!=0)
      continue;
    if(!found) {//fill Uc[dim]
      double V[dim];
      vf->conservativeToPrimitive(U[i],V,fluidId[i]);
      vf->setPressure(p,V,fluidId[i]); //for Tait it actually modifies density.
      vf->primitiveToConservative(V,Uc,fluidId[i]);
      found = true;
    }
    for(int k=0; k<dim; k++)
      U[i][k] = Uc[k];
  }
/*
  bool found = false;
  double rhoe = -1.0; 

// only the pressure with fluidId = 0 (outside part)
// are updated. It is assumed that the input and output states are uniform!
  for(int i=0; i<nodes.size(); i++) {
    if(fluidId[i]!=0) 
      continue;
    if(!found){ //rhoe is not computed yet.
      double V[dim];
      vf->conservativeToPrimitive(U[i],V,fluidId[i]);
      vf->setPressure(p,V,fluidId[i]);
      rhoe = vf->computeRhoEnergy(V,fluidId[i]);
      found = true;
    }    
    U[i][4] = rhoe;
  }
*/
}
//------------------------------------------------------------------------------
template<int dim>
void SubDomain::printVariable(SVec<double,dim> &V, VarFcn *vf)
{
  inletNodes.printVariable(V, sharedInletNodes, vf);
}

//------------------------------------------------------------------------------

template<int dim>
void SubDomain::printInletVariable(SVec<double,dim> &V)
{
  inletNodes.printInletVariable(V, sharedInletNodes);
}

//-----------------------------------------------------------------------------

template<int dim>
void SubDomain::printAllVariable(Vec<int> &X, SVec<double,dim> &U, int numSub, int it)
{

  int glob;
  for (int i=0; i<nodes.size(); i++){
    glob = locToGlobNodeMap[i]+1;
    fprintf(stdout, "Tag[%d,%d] = %d\n", i, glob, X[i]);
  }

}
//------------------------------------------------------------------------------

template<int dim>
void SubDomain::checkExtrapolationValue(SVec<double,dim> &V, VarFcn *vf,
                                        BcData<dim>& bcData, GeoState& geoState)
{
  inletNodes.checkExtrapolationValue(V, sharedInletNodes, nodeType, vf, bcData, geoState);
}

//------------------------------------------------------------------------------

template<class Scalar, int neq>
void SubDomain::printAllMatrix(GenMat<Scalar,neq> &A, int it)
{
    fprintf(stdout, "subDomain %d:\n", locSubNum);
    int glob;
    for (int i=0; i<nodes.size(); i++){
      glob = locToGlobNodeMap[i]+1;
      Scalar *Aii = A.getElem_ii(i);
      fprintf(stdout, "%d %d ", i, glob);
      for (int j = 0; j<neq; j++){
        if(j==0)
          fprintf(stdout, "Aii[%d] = ", i);
        else
          fprintf(stdout, "          ");
        for(int k=0; k<neq; k++){
          fprintf(stdout, "%.6e ", Aii[j*neq+k]);
          fprintf(stdout, "\n");
        }
      }
      fprintf(stdout, "\n");
    }

}
//--------------------------------------------------------------------------
template<int dim>
void SubDomain::padeReconstruction(SVec<double, dim> **dataCoarse, SVec<double, dim> **data, int *stepParam, double *freqCoarse, double deltaFreqFine, int nStrMode,int L, int M, int nPoints)
{


  int i, j;
  int numFreqCoarse = stepParam[1];
  int nSnapsCoarse = stepParam[3];
  int nSteps = stepParam[0];
  int numPadeDeriv = int(ceil((L+M+1.0)/((double)nPoints))-1.0);
  int size = L+M+1;


  bcomp *compMat = new bcomp[numFreqCoarse*(numPadeDeriv+1)];

  double tempMat[dim*nSnapsCoarse];
  bcomp *padeMat = new bcomp[size*size];
  bcomp *padeVec = new bcomp[size];

  bcomp *snaps = new bcomp[nSteps+1];
  int midFreq[1] ;
  double deltaFreqCoarse[numFreqCoarse];
  buildDeltaFreq(deltaFreqCoarse, numFreqCoarse, freqCoarse, midFreq);
  for (int iNode = 0; iNode < (*(data[0])).len; iNode++) {


    extractElementsRelativeToANode(dataCoarse, tempMat, iNode, nSnapsCoarse);
    for (int iDim = 0; iDim < dim; iDim++) {
      for (int iStrMode = 0; iStrMode < nStrMode; iStrMode++) {
        extractElementsRelativeToAComponentAndAMode(tempMat,compMat,iDim,iStrMode,numPadeDeriv,numFreqCoarse,nStrMode,nSnapsCoarse, freqCoarse[0]);
        multiPade(compMat, stepParam, deltaFreqCoarse, padeMat, padeVec, L, M, nPoints, deltaFreqFine, freqCoarse[midFreq[0]], snaps, freqCoarse);
        snapshotsConstruction(data,snaps,nSteps,iDim,iStrMode,iNode,nStrMode,freqCoarse[0]);

      }
    }
  }
}

//--------------------------------------------------------------------------
template<int dim>
void SubDomain::extractElementsRelativeToANode(SVec<double, dim> **dataCoarse, double *tempMat, int iNode, int nSnapsCoarse)
{

  for (int jSnapsCoarse = 0; jSnapsCoarse < nSnapsCoarse; jSnapsCoarse++) {

    SVec<double, dim> *X = dataCoarse[jSnapsCoarse];
    for (int iDim = 0; iDim < dim; iDim++) {
      *(tempMat + iDim*nSnapsCoarse + jSnapsCoarse)  = (*X)[iNode][iDim];  //tempMat[iDim][jSnapsCoarse]

    }

  }
}

//--------------------------------------------------------------------------
template<int dim>
void SubDomain::snapshotsConstruction(SVec<double, dim> **data, bcomp* snaps, int nSteps, int iDim, int iStrMode, int iNode, int nStrMode, double freq1)
{

  SVec<double,dim> *Y = 0;
  SVec<double,dim> *Z = 0;
  if (freq1 == 0.0) {

    for (int jFineFreq=0; jFineFreq<nSteps+1; jFineFreq++) {
      if (jFineFreq==0) {
        Y = data[iStrMode];
        (*Y)[iNode][iDim] = real(snaps[0]);
      }
      else {
        Y = data[nStrMode+2*nStrMode*(jFineFreq-1)+2*iStrMode];
        (*Y)[iNode][iDim] = real(snaps[jFineFreq]);
        Z = data[nStrMode+2*nStrMode*(jFineFreq-1)+2*iStrMode+1];
        (*Z)[iNode][iDim] = imag(snaps[jFineFreq]);
      }
    }
  }
  else {
    for (int jFineFreq=0; jFineFreq<nSteps+1; jFineFreq++) {
      Y = data[2*nStrMode*jFineFreq+2*iStrMode];
      (*Y)[iNode][iDim] = real(snaps[jFineFreq]);
      Z = data[2*nStrMode*jFineFreq+2*iStrMode+1];
      (*Z)[iNode][iDim] = imag(snaps[jFineFreq]);
    }
  }
  Y = 0;
  Z = 0;

}

//------------------------------------------------------------------------------
template<int dim>
void SubDomain::hardyInterpolationLogMap(SVec<double, dim> ***dataCoarse, SVec<double, dim> **dataInterp, int nData, int numPod, int iDataMin, FullM &B, FullM &b)
{

  double tempMat[dim*nData];
  double *hardyCoefs = new double[nData];
  //for each vector
  for (int iPod = 0; iPod < numPod; ++iPod) {
    //for each node
    for (int iNode = 0; iNode < (*(dataInterp[0])).len; ++iNode) {
      extractElementsRelativeToANodeAndAVector(dataCoarse,tempMat,iNode,nData,iDataMin,iPod);
      // for each component
      for (int iDim = 0; iDim < dim; ++iDim) {
        //hardyCoefs = B*data(component)
        for (int iData = 0; iData < nData; ++iData) {
          hardyCoefs[iData] = 0.0;
          for (int jData = 0; jData < nData; ++jData)
            hardyCoefs[iData] += B[iData][jData]*tempMat[iDim*nData+jData];
        }
        //compute the interpolated analog quantity
        (*(dataInterp[iPod]))[iNode][iDim] = 0.0;
        for (int iData = 0; iData < nData; ++iData)
          (*(dataInterp[iPod]))[iNode][iDim] += b[iData][0]*hardyCoefs[iData];
      }
    }
  }
  delete [] hardyCoefs;

}

//------------------------------------------------------------------------------
template<int dim>
void SubDomain::extractElementsRelativeToANodeAndAVector(SVec<double, dim> ***dataCoarse, double *tempMat, int iNode, int nData, int jDataMin, int iPod)
{


  for (int jData = 0; jData < nData; ++jData) {
    if (jData !=jDataMin) {
      SVec<double, dim> *X = dataCoarse[jData][iPod];
      for (int iDim = 0; iDim < dim; ++iDim)
        *(tempMat + iDim*nData + jData)  = (*X)[iNode][iDim];  //tempMat[iDim][jData]
    }
    else {
      for (int iDim = 0; iDim < dim; ++iDim)
        *(tempMat + iDim*nData + jData)  = 0.0;
    }

  }
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::getGradP(NodalGrad<dim>& ngrad)
{

  // gradents stored for future use in computeForce and computeForceTransmitted
  SVec<double,dim> &dVdx = ngrad.getX();
  SVec<double,dim> &dVdy = ngrad.getY();
  SVec<double,dim> &dVdz = ngrad.getZ();

  for (int i=0;i<nodes.size();i++) {
    gradP[0][i] = dVdx[i][4];
    gradP[1][i] = dVdy[i][4];
    gradP[2][i] = dVdz[i][4];
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void SubDomain::getDerivativeOfGradP(NodalGrad<dim>& ngrad)
{

  SVec<double,dim> &ddVdx = ngrad.getXderivative();
  SVec<double,dim> &ddVdy = ngrad.getYderivative();
  SVec<double,dim> &ddVdz = ngrad.getZderivative();

  for (int i=0;i<nodes.size();i++) {
    dGradP[0][i] = ddVdx[i][4];
    dGradP[1][i] = ddVdy[i][4];
    dGradP[2][i] = ddVdz[i][4];
    (*dGradPSVec)[i][0] = ddVdx[i][4];
    (*dGradPSVec)[i][1] = ddVdy[i][4];
    (*dGradPSVec)[i][2] = ddVdz[i][4];
  }

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void SubDomain::getDerivativeOfGradP(RectangularSparseMat<double,dim,3> &dGradPdddx, 
                                     RectangularSparseMat<double,dim,3> &dGradPdddy, 
                                     RectangularSparseMat<double,dim,3> &dGradPdddz, 
                                     SVec<double,dim>& ddVdx, 
                                     SVec<double,dim>& ddVdy, 
                                     SVec<double,dim>& ddVdz, 
                                     SVec<double,3>& dGradPSVec)
{

  SVec<double,3> dummy(dGradPSVec);
  dGradPdddx.apply(ddVdx, dGradPSVec);
  dGradPdddy.apply(ddVdy, dummy);
  dGradPSVec += dummy;
  dGradPdddz.apply(ddVdz, dummy);
  dGradPSVec += dummy;

  for(int i=0; i<3; ++i)
    for(int j=0; j<numNodes(); ++j)
      dGradP[i][j] = dGradPSVec[j][i];
}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void SubDomain::getTransposeDerivativeOfGradP(RectangularSparseMat<double,dim,3> &dGradPdddx, 
                                              RectangularSparseMat<double,dim,3> &dGradPdddy, 
                                              RectangularSparseMat<double,dim,3> &dGradPdddz, 
                                              SVec<double,3>& dGradPSVec,
                                              SVec<double,dim>& ddVdx, 
                                              SVec<double,dim>& ddVdy, 
                                              SVec<double,dim>& ddVdz)
{

  SVec<double,3> dummy(dGradPSVec);
  dGradPdddx.applyTranspose(dGradPSVec, ddVdx);
  dGradPdddy.applyTranspose(dGradPSVec, ddVdy);
  dGradPdddz.applyTranspose(dGradPSVec, ddVdz);

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void SubDomain::computeDerivativeOperatorOfGradP(RectangularSparseMat<double,dim,3> &dGradPdddx,
                                                  RectangularSparseMat<double,dim,3> &dGradPdddy,
                                                  RectangularSparseMat<double,dim,3> &dGradPdddz)
{

  for (int i=0;i<nodes.size();i++) {
    double dGradPddVdxarray[3][dim] = {0};
    double dGradPddVdyarray[3][dim] = {0};
    double dGradPddVdzarray[3][dim] = {0};
    dGradPddVdxarray[0][4] = 1.0;
    dGradPddVdyarray[1][4] = 1.0;
    dGradPddVdzarray[2][4] = 1.0;
    int ndList[1] = {i};
    dGradPdddx.addContrib(1, ndList, dGradPddVdxarray[0]);
    dGradPdddy.addContrib(1, ndList, dGradPddVdyarray[0]);
    dGradPdddz.addContrib(1, ndList, dGradPddVdzarray[0]);
  }

}

//-----------------------------------------------------------------------------

template<int dim>
void SubDomain::computePrdtWCtrlVolRatio(SVec<double,dim> &ratioTimesU, SVec<double,dim> &U, Vec<double> &ctrlVol, GeoState &geoState)
{
   Vec<double>& ctrlVol_n = geoState.getCtrlVol_n();

   for (int i=0; i<nodes.size(); ++i) {
     double ratio = ctrlVol_n[i]/ctrlVol[i];
     for (int j=0; j<dim; ++j) {
       ratioTimesU[i][j] = ratio * U[i][j];
     }
   }

}

//------------------------------------------------------------------------------
//             LEVEL SET SOLUTION AND REINITIALIZATION                       ---
//------------------------------------------------------------------------------
template<int dimLS>
void SubDomain::avoidNewPhaseCreation(SVec<double,dimLS> &Phi, SVec<double,dimLS> &Phin)
{
  if(!NodeToNode)
     NodeToNode = createEdgeBasedConnectivity();
  fprintf(stderr, "***Error: this routine is not valid for multiple subdomains\n");
  exit(1);

  for(int i=0; i<nodes.size(); i++){
    for(int j=0; j<dimLS; j++){
      if(Phi[i][j]*Phin[i][j]<0.0){
        // check if node i HAD a neighbour with a different levelset sign
        bool diffNeigh = false;
        for(int iNeigh=0; iNeigh<NodeToNode->num(i); iNeigh++)
          if(Phin[i][j]*Phin[(*NodeToNode)[i][iNeigh]][j]<0.0){
            diffNeigh = true;
            break;
          }
        if(!diffNeigh) {Phi[i][j] = Phin[i][j]; fprintf(stdout, "node %d has levelset %d clipped to avoid phase creation\n", locToGlobNodeMap[i]+1, j);}
      }
    }
  }

}
//------------------------------------------------------------------------------
template<int dimLS>
void SubDomain::avoidNewPhaseCreation(SVec<double,dimLS> &Phi, SVec<double,dimLS> &Phin, Vec<double> &weight, LevelSetStructure *LSS, Vec<int>* fluidIdToSet)
{

  for(int i=0; i<nodes.size(); i++){
    int fModel;//if fModel>0 (isolated), Phi is not used at all
    if(LSS && !LSS->withCracking()) 
      fModel = LSS->fluidModel(0.0, i); 
    else
      fModel = 0;
    //bool swept = LSS ? LSS->isSwept(0.0, i) : 0;
    for(int j=0; j<dimLS; j++){
      if(Phi[i][j]*Phin[i][j]<0.0/* && fModel==0 && !swept*/){
        // check if node i HAD a neighbour with a different levelset sign
        if(weight[i] <= 0.0 || (fluidIdToSet && (*fluidIdToSet)[i] != j)){
//          fprintf(stdout, "node %d (loc %d in %d) has weight = %f and has levelset %d"
//                          " moving from %e to %e\n", locToGlobNodeMap[i]+1,i,
//                         globSubNum,weight[i],j,Phin[i][j],Phi[i][j]);
          Phi[i][j] = Phin[i][j];
        }
      }
    }
  }

}

//------------------------------------------------------------------------------
template<int dimLS>
void SubDomain::TagInterfaceNodes(int lsdim, Vec<int> &Tag, SVec<double,dimLS> &Phi, int level,LevelSetStructure *LSS)
{
  if(!NodeToNode)
     NodeToNode = createEdgeBasedConnectivity();

  if(level==1){
  // tag nodes that are closest to interface by looking at phi[i]*phi[j]
    Tag = 0;
    edges.TagInterfaceNodes(lsdim,Tag,Phi,LSS);

  }else{
  // tag nodes that are neighbours of already tagged nodes.
    int nNeighs,nei,k,lowerLevel=level-1;
    for(int i=0; i<nodes.size(); i++){

      if(Tag[i]==lowerLevel){

        nNeighs = NodeToNode->num(i);
        for(k=0;k<nNeighs;k++){
          nei = (*NodeToNode)[i][k];
          if(nei==i) continue;
          if(Tag[nei]==0) Tag[nei] = level;
        }

      }
    }
  }
}

//------------------------------------------------------------------------------

template<int dimLS>
void SubDomain::pseudoFastMarchingMethod(Vec<int> &Tag, SVec<double,3> &X,
					 SVec<double,dimLS> &d2wall, int level, int iterativeLevel,
					 Vec<int> &sortedNodes, int &nSortedNodes, int &firstCheckedNode,
					 LevelSetStructure *LSS)
{
  if(!NodeToNode)
     NodeToNode = createEdgeBasedConnectivity();
  if(!NodeToElem) 
     NodeToElem = createNodeToElementConnectivity();
  if(level > 0 && level == iterativeLevel) {  
    nSortedNodes     = 0;
    firstCheckedNode = 0;
    for(int i=0;i<Tag.size();++i) {
      if(Tag[i] < iterativeLevel-1) {
	sortedNodes[nSortedNodes] = i;
	nSortedNodes++;
      }
    }
    for(int i=0;i<Tag.size();++i) {
      if(Tag[i] == iterativeLevel-1) {
	sortedNodes[nSortedNodes] = i;
	nSortedNodes++;
      }
      firstCheckedNode = nSortedNodes;
    }
    for(int i=0;i<Tag.size();++i) {
      if(Tag[i] == iterativeLevel) {
	sortedNodes[nSortedNodes] = i;
	nSortedNodes++;
      }
      if(Tag[i] > iterativeLevel) {
	Tag[i] = -1;
      }
    }
  }
  else if(level == 0) { // just get inactive nodes
    nSortedNodes     = 0;
    firstCheckedNode = 0;
    for(int i=0;i<Tag.size();++i) {
      if(!LSS->isActive(0.0,i))
      {
	d2wall[i][0] = 0.0;
	Tag[i]       = 0;
	sortedNodes[nSortedNodes] = i;
	nSortedNodes++;
      }
    }
  }
  else if(level==1){
//    Tag = -1;  // Tag is globally set to -1. 0 level are inactive nodes
    firstCheckedNode = nSortedNodes;
    edges.pseudoFastMarchingMethodInitialization(X,Tag,d2wall,sortedNodes,nSortedNodes,LSS);
  }
  else{
  // Tag nodes that are neighbours of already Tagged nodes and compute their distance
    int nNeighs,nTets,nei,tet,lowerLevel=level-1;
    int inter = nSortedNodes,fixedNode;
    for(int i=firstCheckedNode;i<inter; i++){
      fixedNode = sortedNodes[i];
//      if(Tag[fixedNode] == 0) continue; //structure nodes
      // Should be useless. Remove the following assert.Adam 2011.09
      // if(Tag[fixedNode]==lowerLevel){
      assert(Tag[fixedNode] == lowerLevel);
      nNeighs = NodeToNode->num(fixedNode);
      for(int k=0;k<nNeighs;k++){
        nei = (*NodeToNode)[fixedNode][k];
	 // Not necessary because of following test.
       // if(nei==i) continue;
        if(Tag[nei]<0) {
          Tag[nei] = level;
          sortedNodes[nSortedNodes] = nei;
          nTets = NodeToElem->num(nei);
          for(int j=0;j<nTets;j++) { 
            tet        = (*NodeToElem)[nei][j];
            elems[tet].FastMarchingDistanceUpdate(nei,Tag,lowerLevel,X,d2wall);
          }
          nSortedNodes++;
        }  
      }
    }
    firstCheckedNode = inter;
  }
}
//------------------------------------------------------------------------------

template<int dimLS>
void SubDomain::TagInterfaceNodes(int lsdim, SVec<bool,2> &Tag, SVec<double,dimLS> &Phi, LevelSetStructure *LSS)
{
  Tag = false;

  int i,j;
  int (*ptr)[2] = edges.getPtr(); 
  for(int l=0; l<edges.size(); l++) {
    i = ptr[l][0];
    j = ptr[l][1];

    if(Phi[i][lsdim]*Phi[j][lsdim]<=0.0) {
      if(LSS->edgeIntersectsStructure(0.0,l))
        Tag[i][0] = Tag[j][0] = true;
      else
        Tag[i][1] = Tag[j][1] = true;
    } else {
      if(LSS->edgeIntersectsStructure(0.0,l)) {
        Tag[i][0] = Tag[j][0] = true;
//        fprintf(stderr,"BUG: Sub %d: (%d,%d) intersects but phi[i]=%e, phi[j]=%e.\n", globSubNum, 
//                locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1, Phi[i][lsdim], Phi[j][lsdim]);
//        fprintf(stderr,"  %d: occluded(%d), swept(%d).\n", locToGlobNodeMap[i]+1, LSS->isOccluded(0,i),LSS->isSwept(0,i));
//        fprintf(stderr,"  %d: occluded(%d), swept(%d).\n", locToGlobNodeMap[j]+1, LSS->isOccluded(0,j),LSS->isSwept(0,j));
      }
    }
  }
}

//------------------------------------------------------------------------------
template<int dimLS>
void SubDomain::printPhi(SVec<double, 3> &X, SVec<double,dimLS> &Phi, int it)
{
  fprintf(stdout, "\nPhi - subDomain %d: \n", locSubNum);
  int glob, sh;
  /*for (int iSub = 0; iSub < numNeighb; ++iSub) {
    for (int iNode = 0; iNode < sharedNodes->num(iSub); ++iNode) {
      sh = (*sharedNodes)[iSub][iNode] ;
      glob = locToGlobNodeMap[sh]+1;
      fprintf(stderr, "%d %.14e %.14e %.14e %.14e %.14e\n", glob, V[ sh ][0],
                                                                  V[ sh ][1],
                                                                  V[ sh ][2],
                                                                  V[ sh ][3],
                                                                  V[ sh ][4]);
    }
  }*/

  for (int i=0; i<nodes.size(); i++){
    glob = locToGlobNodeMap[i]+1;
    fprintf(stdout, " Phi : %d %d   ", i,glob);
    for(int j=0; j<dimLS; j++)
      fprintf(stdout, "%e ", Phi[i][j]);
    fprintf(stdout, "\n");
  }

}
//------------------------------------------------------------------------------

template<int dimLS>
void SubDomain::setupPhiVolumesInitialConditions(const int volid, 
                    const int fluidId, SVec<double,dimLS> &Phi){
  for (int iElem = 0; iElem < elems.size(); iElem++)  {
    if (elems[iElem].getVolumeID() == volid)  {
      int *nodeNums = elems[iElem].nodeNum();
      for (int iNode = 0; iNode < elems[iElem].numNodes(); iNode++){
        for (int iDim = 0; iDim < dimLS; iDim++){
          if(iDim+1 == fluidId)
            Phi[nodeNums[iNode]][iDim] = 1.0;
        }
      }
    }
  }

}

//--------------------------------------------------------------------------
// for mesh motion (with RK2 time-integration)
template<int dimLS>
void SubDomain::computePrdtPhiCtrlVolRatio(SVec<double,dimLS> &ratioTimesPhi,
           SVec<double,dimLS> &Phi, Vec<double> &ctrlVol, GeoState &geoState)
{
   Vec<double>& ctrlVol_n = geoState.getCtrlVol_n();

   for (int i=0; i<nodes.size(); ++i) {
     double ratio = ctrlVol_n[i]/ctrlVol[i];
     for (int j=0; j<dimLS; j++)
       ratioTimesPhi[i][j] = ratio * Phi[i][j];
   }

}

//-----------------------------------------------------------------------------
/*
template<int dimLS>
void SubDomain::FinishReinitialization(Vec<int> &Tag, SVec<double,dimLS> &Psi, int level)
{

  if(!NodeToNode)
    NodeToNode = createEdgeBasedConnectivity();

  int nNeighs,nei,k;
  for (int i=0; i<nodes.size(); i++){
    if (Tag[i]==level){

      nNeighs = NodeToNode->num(i);
      for (k=0; k<nNeighs; k++){
        nei = (*NodeToNode)[i][k];
        if(Tag[nei]==0){
          Tag[nei] = level+1;
          Psi[nei][0] = Psi[i][0];
        }else if(Tag[nei]==level+1){
          if( (Psi[i][0] > 0.0 && Psi[i][0] > Psi[nei][0]) ||
              (Psi[i][0] < 0.0 && Psi[i][0] < Psi[nei][0])  )
            Psi[nei][0] = Psi[i][0];
        }

      }
    }
  }

}
*/
//-----------------------------------------------------------------------------
template<int dimLS>
void SubDomain::copyCloseNodes(int lsdim, int level, Vec<int> &Tag,SVec<double,dimLS> &Phi,SVec<double,1> &Psi)
{
  for(int i=0; i<nodes.size(); i++) {
    if(Tag[i]==level)
      Psi[i][0] = fabs(Phi[i][lsdim]);
    //DEBUG
//    if(locToGlobNodeMap[i]+1==115697)
//      fprintf(stderr,"Sub %d, level = %d, Tag[%d] = %d. Phi = %e, Psi = %e.\n", globSubNum, level, locToGlobNodeMap[i]+1,
//              Tag[i], Phi[i][0], Psi[i][0]);
  }

}
//------------------------------------------------------------------------------
template<int dimLS>
void SubDomain::computeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                       NodalGrad<dimLS> &grad,
                                       SVec<double,dimLS> &Phi,SVec<double,1> &Psi)
{
  for(int i=0; i<nodes.size(); i++){
    if(Tag[i]==1) Tag[i]=-1;
  }
  SVec<double,dimLS>& ddx  = grad.getX();
  SVec<double,dimLS>& ddy  = grad.getY();
  SVec<double,dimLS>& ddz  = grad.getZ();
  elems.computeDistanceCloseNodes(lsdim,Tag,X,ddx,ddy,ddz,Phi,Psi);
}
//-------------------------------------------------------------------------------
template<int dimLS>
void SubDomain::recomputeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                      NodalGrad<dimLS> &grad, SVec<double,dimLS> &Phi,
                                      SVec<double,1> &Psi)
{

  SVec<double,dimLS>& ddx  = grad.getX();
  SVec<double,dimLS>& ddy  = grad.getY();
  SVec<double,dimLS>& ddz  = grad.getZ();
  elems.recomputeDistanceCloseNodes(lsdim,Tag,X,ddx,ddy,ddz,Phi,Psi);

}
//-------------------------------------------------------------------------------
template<int dimLS>
double SubDomain::computeDistanceLevelNodes(int lsdim, Vec<int> &Tag, int level,
                                       SVec<double,3> &X, SVec<double,1> &Psi,SVec<double,dimLS> &Phi)
{

  if(level==2)
		for(int i=0; i<nodes.size(); i++) if(Tag[i]==-1) Tag[i]=1;

  elems.computeDistanceLevelNodes(lsdim,Tag,level,X,Psi,Phi);
  double res = 0.0;
  for(int i=0; i<nodes.size(); i++)
		if(Tag[i]==level)	res += Psi[i][0]*Psi[i][0];

  return res;
}
//-------------------------------------------------------------------------------

template<int dimLS>
void SubDomain::getSignedDistance(int lsdim, SVec<double,1> &Psi, SVec<double,dimLS> &Phi)
{
  for(int i=0; i<nodes.size(); i++){
    if(Phi[i][lsdim]<0.0)
      Psi[i][0] = -Psi[i][0];
    if(Phi[i][lsdim]<0.0 && Psi[i][0]>0.0)
      fprintf(stdout, "globnode %d (%d) has changed phase %e %e\n", locToGlobNodeMap[i]+1,i,Phi[i][lsdim],Psi[i][0]);
  }
}

//-----------------------------------------------------------------------------------------------

class ElemForceCalcValid {

  public:

    bool Valid(Elem*) { return true; }
};

//-----------------------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeEmbSurfBasedForceLoad(IoData &iod, int forceApp, int order, 
															SVec<double,3> &X, double (*Fs)[3], int sizeFs, 
															int numStructElems, int (*stElem)[3], Vec<Vec3D>& Xstruct, 
															LevelSetStructure &LSS, double pInfty, 
                                             SVec<double,dim> &Wstarij, SVec<double,dim> &Wstarji, 
															SVec<double,dim> &V, Vec<GhostPoint<dim>*> *ghostPoints, 
															PostFcn *postFcn, NodalGrad<dim, double> &ngrad, 
															VarFcn* vf, Vec<int>* fid)
{

  if (forceApp!=2) 
   {
		fprintf(stderr,"ERROR: force method (%d) not recognized! Abort..\n", forceApp); 
		exit(-1);
	}

   if (iod.problem.framework == ProblemData::EMBEDDEDALE)
     myTree->reconstruct<&Elem::computeBoundingBox>(X, elems.getPointer(), elems.size());

  int qOrder = iod.embed.qOrder; // default is 3
  Quadrature quadrature_formula(qOrder);
  int nqPoint = quadrature_formula.n_point;
  double (*qloc)[3](quadrature_formula.qloc);
  double *qweight(quadrature_formula.weight);

  int iElem = -1;
  int T[4];       //nodes in a tet.

  Vec3D vectorIJ, gradP, flocal;
  SVec<double,dim> gradX = ngrad.getX();
  SVec<double,dim> gradY = ngrad.getY();
  SVec<double,dim> gradZ = ngrad.getZ();

  CrackingSurface* cs = LSS.getCrackingSurface();
 
  double Vext[dim]; 
  int stNode[3];
  Vec3D Xst[3];
  Vec3D Xp;

	for(int nSt=0; nSt<numStructElems; ++nSt) 
	{	  
		for(int j=0; j<3; ++j) 
		{
      stNode[j] = stElem[nSt][j];
         Xst[j] = Xstruct[stNode[j]]; 
    }
    Vec3D normal = 0.5*(Xst[1]-Xst[0])^(Xst[2]-Xst[0]);

		for(int nq=0; nq<nqPoint; ++nq) 
		{
			for(int j=0; j<3; ++j) 
	Xp[j] = qloc[nq][0]*Xst[0][j] + qloc[nq][1]*Xst[1][j] + qloc[nq][2]*Xst[2][j];

      ElemForceCalcValid myObj;
      Elem* E = myTree->search<&Elem::isPointInside, ElemForceCalcValid,
	             &ElemForceCalcValid::Valid>(&myObj, X, Xp);

      if (!E) continue;      

      // Check to see if the structure has cracked, and this quadrature point
      // falls on a portion of the structure that is phantom (i.e., no longer
      // exists)
      if (cs && cs->getPhi(nSt, qloc[nq][0], qloc[nq][1]) < 0.0)
      	continue;
 

      for (int i=0; i<4; i++) T[i] = (*E)[i];

			Vec3D Xf[4]; 
			for (int i=0; i<4; i++)
				for(int j=0;j<3;++j) Xf[i][j] = X[T[i]][j];

			// Compute barycentric coordinates
			Vec3D bary;
			E->computeBarycentricCoordinates(X,Xp,bary); 
			if (bary[0] < 0.0 || bary[1] < 0.0 || bary[2] < 0.0 || bary[0]+bary[1]+bary[2] > 1.0) 
			{
				E = 0;
				continue;
			}

			Vec3D dbary[4];
			dbary[0] = Vec3D(1.0-bary[0],bary[1],bary[2]);
			dbary[1] = Vec3D(bary[0],1.0-bary[1],bary[2]);
			dbary[2] = Vec3D(bary[0],bary[1],1.0-bary[2]);
			dbary[3] = Vec3D(bary[0],bary[1],bary[2]);

			// Determine the side of the nodes of the tet on intersected edges
			int norm[4] = {0, 0, 0, 0};
			for (int e=0; e<6; ++e) 
			{
				int l = E->edgeNum(e);
				if(LSS.edgeIntersectsStructure(0,l)) 
				{
					int i = E->edgeEnd(e,0);
					int j = E->edgeEnd(e,1);

					LevelSetResult lsResij = LSS.getLevelSetDataAtEdgeCenter(0.0, l, (T[i]<T[j]));

					norm[i] = (lsResij.gradPhi*(Xstruct[lsResij.trNodes[0]]-Xf[i]) <= 0) ? -1 : 1;

					LevelSetResult lsResji = LSS.getLevelSetDataAtEdgeCenter(0.0, l, (T[i]>=T[j]));

					norm[j] = (lsResji.gradPhi*(Xstruct[lsResji.trNodes[0]]-Xf[j]) <= 0) ? -1 : 1;
				}
			}

			// Check for the dual volume using barycentric coordinates) of the tet on the either side of the surface element
			double mindist[2] = {FLT_MAX,FLT_MAX};

			int node[2] = {-1,-1};

			Vec3D nf[2] = {-normal, normal};

			for (int i=0; i<4; i++) 
			{
				double dist = dbary[i].norm();
				if (norm[i] < 0) 
				{
	  
					// Bug fix for cracking simulations (also below)
					// Note that when we are doing cracking, isActive() always
					// returns false.  In this case, the node is assumed to be active.
					// so we only need to check if the node is occluded.
					if( (LSS.isActive(0,T[i]) || (cs && !LSS.isOccluded(0,T[i]))) && dist < mindist[0] && normal*(Xp-Xf[i]) <= 0. ) 
					{
						mindist[0] = dist;
						node[0] = T[i];
					}
				}
				else if(norm[i] > 0) 
				{
					if( (LSS.isActive(0,T[i])|| (cs && !LSS.isOccluded(0,T[i]))) && dist < mindist[1] && normal*(Xp-Xf[i]) > 0. ) 
					{
						mindist[1] = dist;
						node[1] = T[i];
					}
				}
			}

			// For viscous simulation
			double dp1dxj[4][3]; // Gradient of the P1 basis functions

			for(int i=0;i<4;++i) for(int j=0;j<3;++j) dp1dxj[i][j] = 0.0;

			double d2w[3]; // not used, but required by postFcn->computeViscousForce(...)

			d2w[0] = d2w[1] = d2w[2] = 0.0;

			double *Vwall = 0;

			double *Vface[3] = {0,0,0};
			
			double *vtet[2][4];
			if(ghostPoints) 
			{
				E->computeGradientP1Function(X, dp1dxj);
				for(int i=0; i<4; ++i) 
				{
					vtet[0][i] = V[T[i]];
					vtet[1][i] = V[T[i]];
				}

				GhostPoint<dim> *gp;
				for(int i=0; i<4; ++i) 
				{
					gp = (*ghostPoints)[T[i]];
					if(gp) 
					{
						if(norm[i] <= 0.) vtet[1][i] = gp->getPrimitiveState();
						else              vtet[0][i] = gp->getPrimitiveState();
					}
				}
			}  

			flocal = 0.0;

			for(int n = 0; n < 2; ++n) 
			{
				int i = node[n];

				if (i < 0) continue;

				double *v = V[i];

				/*gradP[0] = gradX[i][4];
				  gradP[1] = gradY[i][4];
				  gradP[2] = gradZ[i][4];
				*/
				for(int m=0;m<3;++m) vectorIJ[m] = Xp[m] - X[i][m];

				for (int k = 0; k < dim; ++k) 
				{
					Vext[k] = v[k] + gradX[i][k]*vectorIJ[0] 
						            + gradY[i][k]*vectorIJ[1] 
                					+ gradZ[i][k]*vectorIJ[2];
				}

				// check for neg pressures/densities
				if(vf->doVerification()) 
				{
					double Udummy[dim];

					Vext[0] = std::max(Vext[0], 1e-10);

					vf->getVarFcnBase(fid?(*fid)[i]:0)->verification(0,Udummy, Vext);
				}

				double pp = vf->getPressure(Vext, fid?(*fid)[i]:0);

				flocal += (pp - pInfty)*nf[n];

				if(ghostPoints) {
					// Viscous Simulation
					flocal += postFcn->computeViscousForce(dp1dxj,nf[n],d2w,Vwall,Vface,vtet[n]);
				}
			}

			for(int j=0; j<3; ++j) 
			{
				Fs[stNode[0]][j] += qweight[nq]*flocal[j]*qloc[nq][0];
				Fs[stNode[1]][j] += qweight[nq]*flocal[j]*qloc[nq][1];
				Fs[stNode[2]][j] += qweight[nq]*flocal[j]*qloc[nq][2];
			}
		}
	}

}

//-----------------------------------------------------------------------------------------------
//d2d
template<int dim,int dimLS>
void SubDomain::computeEMBNodeScalarQuantity(SVec<double,3> &X, SVec<double,dim> &V, 
															PostFcn *postFcn, VarFcn *varFcn, 
															Vec<int> &fluidId, SVec<double,dimLS>* phi,
															double (*Qnty)[3], int sizeQnty, int numStructElems, int (*stElem)[3],
															Vec<Vec3D>& Xstruct, LevelSetStructure &LSS,
															double pInfty, 
															Vec<GhostPoint<dim>*> *ghostPoints,
															NodalGrad<dim, double> &ngrad)
{

  myTree->reconstruct<&Elem::computeBoundingBox>(X, elems.getPointer(), elems.size());

  int qOrder = 1;
  Quadrature quadrature_formula(qOrder);
  int nqPoint = quadrature_formula.n_point;
  double (*qloc)[3](quadrature_formula.qloc);
  double *qweight(quadrature_formula.weight);

  int iElem = -1;
  int T[4];       //nodes in a tet.


  Vec3D vectorIJ, gradP, flocal;
  SVec<double,dim> gradX = ngrad.getX();
  SVec<double,dim> gradY = ngrad.getY();
  SVec<double,dim> gradZ = ngrad.getZ();
  
  CrackingSurface* cs = LSS.getCrackingSurface();

  // if(cs) {
  //   fprintf(stderr, "Warning: wall quantities cannot be computed on a cracked embedded surface.\n");
  //   return;
  // }

  double Vext[dim]; 
  int stNode[3];
  Vec3D Xst[3];
  Vec3D Xp;

  int nq = 0;
  
  for(int nSt = 0; nSt < numStructElems; ++nSt) {

    for (int j=0; j<3; ++j) {
      stNode[j] = stElem[nSt][j];
         Xst[j] = Xstruct[stNode[j]]; 
    }
    Vec3D normal = 0.5*(Xst[1]-Xst[0])^(Xst[2]-Xst[0]);

    for (int j=0; j<3; ++j){
	Xp[j] = qloc[nq][0]*Xst[0][j] + qloc[nq][1]*Xst[1][j] + qloc[nq][2]*Xst[2][j];
    }
    
    ElemForceCalcValid myObj;
    Elem* E = myTree->search<&Elem::isPointInside, ElemForceCalcValid,
			     &ElemForceCalcValid::Valid>(&myObj, X, Xp);

    if (!E) 
      continue;

    if (cs && cs->getPhi(nSt, qloc[nq][0], qloc[nq][1]) < 0.0)
      continue;
 
    for (int i=0; i<4; i++) T[i] = (*E)[i];

      Vec3D Xf[4]; 
      for (int i=0; i<4; i++)
        for(int j=0;j<3;++j) Xf[i][j] = X[T[i]][j];

// Compute barycentric coordinates
      Vec3D bary;
      E->computeBarycentricCoordinates(X,Xp,bary); 
      if (bary[0] < 0.0 || bary[1] < 0.0 || bary[2] < 0.0 || bary[0]+bary[1]+bary[2] > 1.0) {
        E = 0;
	continue;
      }

      Vec3D dbary[4];
      dbary[0] = Vec3D(1.0-bary[0],bary[1],bary[2]);
      dbary[1] = Vec3D(bary[0],1.0-bary[1],bary[2]);
      dbary[2] = Vec3D(bary[0],bary[1],1.0-bary[2]);
      dbary[3] = Vec3D(bary[0],bary[1],bary[2]);

// Determine the side of the nodes of the tet on intersected edges
      int norm[4] = { 0, 0, 0, 0 };
      for (int e=0; e<6; ++e) {
	int l = E->edgeNum(e);
	if (LSS.edgeIntersectsStructure(0,l)) {
	  int i = E->edgeEnd(e,0);
	  int j = E->edgeEnd(e,1);
          LevelSetResult lsResij = LSS.getLevelSetDataAtEdgeCenter(0.0, l, (T[i]<T[j]));
          norm[i] = (lsResij.gradPhi*(Xstruct[lsResij.trNodes[0]]-Xf[i]) <= 0) ? -1 : 1;
          LevelSetResult lsResji = LSS.getLevelSetDataAtEdgeCenter(0.0, l, (T[i]>=T[j]));
          norm[j] = (lsResji.gradPhi*(Xstruct[lsResji.trNodes[0]]-Xf[j]) <= 0) ? -1 : 1;
	}
      }

// Check for the dual volume using barycentric coordinates) of the tet on the either side of the surface element
      double mindist[2] = {FLT_MAX,FLT_MAX};
      int node[2] = {-1,-1};
      Vec3D nf[2] = {-normal,normal};
      for (int i=0; i<4; i++) {
	double dist = dbary[i].norm();
        if (norm[i] < 0) {  
	  
	if( (LSS.isActive(0,T[i]) || (cs && !LSS.isOccluded(0,T[i]))) && dist < mindist[0] && normal*(Xp-Xf[i]) <= 0. ) {
	  mindist[0] = dist;
	  node[0] = T[i];
	}
      } else if(norm[i] > 0) {
	if( (LSS.isActive(0,T[i])|| (cs && !LSS.isOccluded(0,T[i]))) && dist < mindist[1] && normal*(Xp-Xf[i]) > 0. ) {
	  mindist[1] = dist;
	  node[1] = T[i];
	}
      }
    }

    
    // ---- For viscous simulation ----

    // Gradient of the P1 basis functions
    double dp1dxj[4][3]; 
    for(int i=0;i<4;++i) {
      for(int j=0;j<3;++j) { 
	dp1dxj[i][j] = 0.0;
      }
    }

    // not used, but required by postFcn->computeViscousForce(...)
    double d2w[3]; 
    d2w[0] = d2w[1] = d2w[2] = 0.0;
    double *Vwall = 0;
    double *Vface[3] = {0,0,0};
      
    double *vtet[2][4];
    if(ghostPoints) {
      E->computeGradientP1Function(X, dp1dxj);
      for(int i=0; i<4; ++i) {
	vtet[0][i] = V[T[i]];
	vtet[1][i] = V[T[i]];
      }
      
      GhostPoint<dim> *gp;
      for(int i=0; i<4; ++i) {
	gp = (*ghostPoints)[T[i]];
	if (gp) {
	  if (norm[i] <= 0.) vtet[1][i] = gp->getPrimitiveState();
	  else               vtet[0][i] = gp->getPrimitiveState();
	}
      }
    }
        
    double Cplocal = 0.0, Cflocal = 0.0;
    for (int n = 0; n < 2; ++n) {

      int i = node[n];
      if (i < 0) continue;
      double *v = V[i];

      for(int m=0;m<3;++m) {
	vectorIJ[m] = Xp[m] - X[i][m];
      }
      for (int k = 0; k < dim; ++k) {
	Vext[k] = v[k] + gradX[i][k]*vectorIJ[0]+
                         gradY[i][k]*vectorIJ[1]+
                         gradZ[i][k]*vectorIJ[2];
      }

      double S = sqrt(nf[n]*nf[n]);

      int fid(0);
      fid = fluidId[i]?fluidId[i]:0;	
      double pp = postFcn->computeNodeScalarQuantity(PostFcn::PRESSURECOEFFICIENT, Vext, Xp, fid, NULL);
      Cplocal += pp;

      // Viscous Simulation
      if(ghostPoints) {
	Vec3D t(1.0, 0.0, 0.0);
	Vec3D F = postFcn->computeViscousForce(dp1dxj, nf[n], d2w, Vwall, Vface, vtet[n]);
	Cflocal += 2.0 * t * F / S;
      }

      Qnty[stNode[0]][0] += S;
      Qnty[stNode[1]][0] += S;
      Qnty[stNode[2]][0] += S;

      Qnty[stNode[0]][1] += Cplocal*S;
      Qnty[stNode[1]][1] += Cplocal*S;
      Qnty[stNode[2]][1] += Cplocal*S;

      Qnty[stNode[0]][2] += Cflocal*S;
      Qnty[stNode[1]][2] += Cflocal*S;
      Qnty[stNode[2]][2] += Cflocal*S;

    }	
    
  }

}

//-----------------------------------------------------------------------------------------------


template<int dim>
void SubDomain::computederivativeEmbSurfBasedForceLoad(IoData &iod, int forceApp, int order, SVec<double,3> &X,
						       double (*dFs)[3], int sizeFs, int numStructElems, int (*stElem)[3], 
						       Vec<Vec3D>& Xstruct, Vec<Vec3D>& dXstruct, LevelSetStructure &LSS, 
						       double pInfty, double dpInfty, 
						       SVec<double,dim> &Wstarij, SVec<double,dim> &Wstarji, 
						       SVec<double,dim> &V, SVec<double,dim> &dV_, 
						       Vec<GhostPoint<dim>*> *ghostPoints, PostFcn *postFcn, 
						       NodalGrad<dim, double> &gradV, NodalGrad<dim, double> &graddV, VarFcn* vf, Vec<int>* fid){

  int qOrder = iod.embed.qOrder;
  Quadrature quadrature_formula(qOrder);
  int nqPoint = quadrature_formula.n_point;
  double (*qloc)[3](quadrature_formula.qloc);
  double *qweight(quadrature_formula.weight);

  int iElem = -1;
  int T[4];

  Vec3D vectorIJ, dvectorIJ, gradP, dflocal;

  SVec<double,dim> gradVX = gradV.getX();
  SVec<double,dim> gradVY = gradV.getY();
  SVec<double,dim> gradVZ = gradV.getZ();

  SVec<double,dim> graddVX = graddV.getX();
  SVec<double,dim> graddVY = graddV.getY();
  SVec<double,dim> graddVZ = graddV.getZ();

  int stNode[3];
  Vec3D Xst[3], dXst[3];
  Vec3D Xp, dXp;

  for(int nSt = 0; nSt < numStructElems; ++nSt) {

    for (int j=0; j<3; ++j) {
      stNode[j] = stElem[nSt][j];
         Xst[j] =  Xstruct[stNode[j]]; 
        dXst[j] = dXstruct[stNode[j]]; 
    }

    Vec3D normal  = 0.5*(Xst[1]-Xst[0])^(Xst[2]-Xst[0]);

    Vec3D dnormal = 0.5*(( (dXst[1] - dXst[0])^( Xst[2] -  Xst[0]) )
                       + ( ( Xst[1] -  Xst[0])^(dXst[2] - dXst[0]) ));

    for(int nq=0; nq<nqPoint; ++nq) {

      for (int j=0; j<3; ++j){
	 Xp[j] = qloc[nq][0]* Xst[0][j] + qloc[nq][1]* Xst[1][j] + qloc[nq][2]* Xst[2][j];
	dXp[j] = qloc[nq][0]*dXst[0][j] + qloc[nq][1]*dXst[1][j] + qloc[nq][2]*dXst[2][j];
      }

      ElemForceCalcValid myObj;
      Elem* E = myTree->search<&Elem::isPointInside, ElemForceCalcValid,
	                       &ElemForceCalcValid::Valid>(&myObj, X, Xp);

      if (!E) continue;

      for (int i=0; i<4; i++) T[i] = (*E)[i];

      Vec3D Xf[4]; 
      for (int i=0; i<4; i++)
        for(int j=0;j<3;++j) Xf[i][j] = X[T[i]][j];

      // Compute barycentric coordinates
      Vec3D bary;
      E->computeBarycentricCoordinates(X,Xp,bary); 

      if (bary[0] < 0.0 || bary[1] < 0.0 || bary[2] < 0.0 || bary[0]+bary[1]+bary[2] > 1.0) {
        E = 0;
	continue;
      }

      Vec3D dbary[4];
      dbary[0] = Vec3D(1.0-bary[0],bary[1],bary[2]);
      dbary[1] = Vec3D(bary[0],1.0-bary[1],bary[2]);
      dbary[2] = Vec3D(bary[0],bary[1],1.0-bary[2]);
      dbary[3] = Vec3D(bary[0],bary[1],bary[2]);

      // Determine the side of the nodes of the tet on intersected edges
      int norm[4] = { 0, 0, 0, 0 };

      for (int e=0; e<6; ++e) {

	int l = E->edgeNum(e);

	if (LSS.edgeIntersectsStructure(0,l)) {

	  int i = E->edgeEnd(e,0);
	  int j = E->edgeEnd(e,1);

          LevelSetResult lsResij = LSS.getLevelSetDataAtEdgeCenter(0.0, l, (T[i]<T[j]));
          norm[i] = (lsResij.gradPhi*(Xstruct[lsResij.trNodes[0]]-Xf[i]) <= 0) ? -1 : 1;

          LevelSetResult lsResji = LSS.getLevelSetDataAtEdgeCenter(0.0, l, (T[i]>=T[j]));
          norm[j] = (lsResji.gradPhi*(Xstruct[lsResji.trNodes[0]]-Xf[j]) <= 0) ? -1 : 1;
	
	}

      }

      double mindist[2] = {FLT_MAX, FLT_MAX};

      int node[2] = {-1,-1};

      Vec3D  nf[2] = { -normal,  normal};
      Vec3D dnf[2] = {-dnormal, dnormal};
      
      for (int i=0; i<4; i++) {

	double dist = dbary[i].norm();

        if (norm[i] < 0) {

	  if( LSS.isActive(0,T[i]) && dist < mindist[0] && normal*(Xp-Xf[i]) <= 0.0 ) {
	    mindist[0] = dist;
	    node[0] = T[i];
	  }

	} else if(norm[i] > 0) {

	  if( LSS.isActive(0,T[i]) && dist < mindist[1] && normal*(Xp-Xf[i]) > 0.0 ) {
	    mindist[1] = dist;
	    node[1] = T[i];
	  }

	}

      }

      dflocal = 0.0;
      for (int n = 0; n < 2; ++n) {

	int i = node[n];
	if (i < 0) continue;

	double *Vi = V[i];

        for(int m=0; m<3; ++m){
           vectorIJ[m] =  Xp[m] - X[i][m];
	  dvectorIJ[m] = dXp[m];
        }

	double Pe = Vi[4] + gradVX[i][4]*vectorIJ[0]+
                            gradVY[i][4]*vectorIJ[1]+
	                    gradVZ[i][4]*vectorIJ[2];

        double dPeS = gradVX[i][4]*dvectorIJ[0]+
	              gradVY[i][4]*dvectorIJ[1]+
 	              gradVZ[i][4]*dvectorIJ[2];

	double dPew = dV_[i][4] + graddVX[i][4]*vectorIJ[0]+
                                  graddVY[i][4]*vectorIJ[1]+
                                  graddVZ[i][4]*vectorIJ[2];

	double dPe = dPeS + dPew;
	//           ****

	dflocal += ( (dPe - dpInfty)*nf[n] + (Pe - pInfty)*dnf[n] );

      }	

      for (int j=0; j<3; ++j) {
        dFs[stNode[0]][j] += qweight[nq]*dflocal[j]*qloc[nq][0];
        dFs[stNode[1]][j] += qweight[nq]*dflocal[j]*qloc[nq][1];
        dFs[stNode[2]][j] += qweight[nq]*dflocal[j]*qloc[nq][2];
      }
    }
  }

}

//-----------------------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeEmbSurfBasedForceLoad_e(IoData &iod, int forceApp, int order, 
															  SVec<double,3> &X, double (*Fs)[3], int sizeFs, 
															  int numStructElems, int (*stElem)[3], Vec<Vec3D>& Xstruct, 
															  LevelSetStructure &LSS, double pInfty,
															  SVec<double,dim> &V, Vec<GhostPoint<dim>*> *ghostPoints, 
															  PostFcn *postFcn, NodalGrad<dim, double> &ngrad, 
															  VarFcn* vf, Vec<int>* fid,
															  int** stNodeDir, double** stX1, double** stX2)
{

	if (forceApp!=2)
   {
		fprintf(stderr,"ERROR: force method (%d) not recognized! Abort..\n", forceApp); 
		exit(-1);
	}
	
	int qOrder = iod.embed.qOrder; // default is 3
	Quadrature quadrature_formula(qOrder);

	int nqPoint = quadrature_formula.n_point;

	double (*qloc)[3](quadrature_formula.qloc);
	double *qweight(quadrature_formula.weight);

	int T[4];

	Vec3D DX, gradP, flocal;

	SVec<double,dim> gradX = ngrad.getX();
	SVec<double,dim> gradY = ngrad.getY();
	SVec<double,dim> gradZ = ngrad.getZ();

	CrackingSurface* cs = LSS.getCrackingSurface();

	double Vext[dim]; 
	double *Vtet[4];
	int stNode[3];
	Vec3D Xst[3];
	Vec3D Xp;

	double Surf;

	double dp1dxj[4][3]; 

	for(int i=0;i<4;++i) 
		for(int j=0;j<3;++j) dp1dxj[i][j] = 0.0;

	for(int nSt=0; nSt<numStructElems; ++nSt) 
	{
		for(int j=0; j<3; ++j) 
		{
			stNode[j] = stElem[nSt][j];
            Xst[j] = Xstruct[stNode[j]]; 
		}

		Vec3D stNormal = 0.5*(Xst[1]-Xst[0])^(Xst[2]-Xst[0]);

		Surf = sqrt(stNormal*stNormal);
		
		if(Surf != 0) stNormal *= (1.0/Surf);

		int node_e[2] = {-1, -1};

		Elem* elem_e1;
		Elem* elem_e2;

		double mindist[2] = {FLT_MAX, FLT_MAX};

		for(int dir=0; dir<2; ++dir)
		{			
			if(stNodeDir[nSt][dir] == 0) continue;
				
			Vec3D Xe;
			if(dir == 0) for(int i=0; i<3; ++i) Xe[i] = stX1[nSt][i];
			else         for(int i=0; i<3; ++i) Xe[i] = stX2[nSt][i];
	
			ElemForceCalcValid myObj;
			Elem* En = myTree->search<&Elem::isPointInside, ElemForceCalcValid,
											  &ElemForceCalcValid::Valid>(&myObj, X, Xe);
						
			if(!En) continue;

			if(dir == 0) elem_e1 = En;
			else         elem_e2 = En;

			for(int k=0; k<4; ++k) 
			{
				int Nk = (*En)[k];
				
				if(LSS.isActive(0.0, Nk))
				{
					Vec3D Xk;
					for(int i=0; i<3; ++i) Xk[i] = X[Nk][i];
					
					Vec3D vDist = Xp - Xk;
					double dist = sqrt(vDist*vDist);
						
					if(dist < mindist[dir])
					{
						mindist[dir] = dist;
						node_e[dir] = Nk;
					}
				}
			}
		}

		// -------------------------------------------------

		for(int nq=0; nq<nqPoint; ++nq) 
		{
			for(int j=0; j<3; ++j) Xp[j] = qloc[nq][0]*Xst[0][j] 
												  + qloc[nq][1]*Xst[1][j] 
 												  + qloc[nq][2]*Xst[2][j];
	
			flocal = 0.0;

			Vec3D stN[2] = {-stNormal*Surf, stNormal*Surf};

			for(int dir=0; dir<2; ++dir)
			{
				int Ni = node_e[dir];

				if(Ni < 0) continue;

            // ------------------ Inviscid part  ------------------ //
				for(int j=0; j<3; ++j) DX[j] = Xp[j] - X[Ni][j];
					
				for(int k = 0; k < dim; ++k) Vext[k] = V[Ni][k]
													          + gradX[Ni][k]*DX[0] 
													          + gradY[Ni][k]*DX[1] 
															    + gradZ[Ni][k]*DX[2];

				double pp = vf->getPressure(Vext, fid?(*fid)[Ni]:0);

				flocal += (pp - pInfty)*stN[dir];
				// ---------------------------------------------------- //

				// ------------------- Viscous part  ------------------ //
				if(ghostPoints)
				{
					GhostPoint<dim> *gp;

					Elem* elem_tmp;
					if(dir == 0) elem_tmp = elem_e1;
					else         elem_tmp = elem_e2;						

					elem_tmp->computeGradientP1Function(X, dp1dxj);

					for(int k=0; k<4; ++k)
					{
						int Nk = (*elem_tmp)[k];

						Vec3D Xwall;
						bool isNGhost = LSS.xWallNode(Nk, Xwall);

						if(isNGhost)
						{
							gp = (*ghostPoints)[Nk];						

							Vec3D Xiw, Xkw;
							for(int j=0; j<3; ++j) 
							{
								Xiw[j] = X[Ni][j] - Xwall[j];
								Xkw[j] = X[Nk][j] - Xwall[j];
							}
					
							double d1 = Xiw*stN[dir];
							double d2 = Xkw*stN[dir];
							
							int df = (d1*d2 > 0.0) ? 1 : -1;

							Vtet[k] = gp->getPrimitiveState(df);
						}
						else
							Vtet[k] = V[Nk];
					}

					// Dummy values (not used)
					double d2w[3] = {0.0, 0.0, 0.0}; 
					double *Vwall = 0;
					double *Vface[3] = {0,0,0};
					
					flocal += postFcn->computeViscousForce(dp1dxj, stN[dir], d2w, Vwall, Vface, Vtet);
				}
				// ---------------------------------------------------- //

				for(int j=0; j<3; ++j) 
				{
					Fs[stNode[0]][j] += qweight[nq]*flocal[j]*qloc[nq][0];
					Fs[stNode[1]][j] += qweight[nq]*flocal[j]*qloc[nq][1];
					Fs[stNode[2]][j] += qweight[nq]*flocal[j]*qloc[nq][2];
				}
			}
		}
	}
	
}


//-----------------------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeEMBNodeScalarQuantity_step1(SVec<double,3> &X, SVec<double,dim> &V,
																	int numStructElems, int (*stElem)[3],
																	Vec<Vec3D>& Xstruct, LevelSetStructure &LSS,
	                                                int** stNodeDir, double** stX1, double** stX2,
																	bool rebuildTree)
{

	if(rebuildTree) myTree->reconstruct<&Elem::computeBoundingBox>(X, elems.getPointer(), elems.size());

	const int qOrder = 1;
	Quadrature quadrature_formula(qOrder);

	int nqPoint = quadrature_formula.n_point;
	const int nq = 0;

	double (*qloc)[3](quadrature_formula.qloc);
	double *qweight(quadrature_formula.weight);

	int T[4];

	CrackingSurface* cs = LSS.getCrackingSurface();

	int stNode[3];
	Vec3D Xst[3];
	Vec3D Xp;

	double Surf;

	const double dtol_plus  = 1.1;
	const double dtol_minus = 0.9;

	for(int nSt=0; nSt<numStructElems; ++nSt) 
	{	
		for(int j=0; j<3; ++j) 
		{
			stNode[j] = stElem[nSt][j];
            Xst[j] = Xstruct[stNode[j]]; 
		}

		Vec3D stNormal = 0.5*(Xst[1]-Xst[0])^(Xst[2]-Xst[0]);
		
		Surf = sqrt(stNormal*stNormal);
		
		if(Surf != 0) stNormal *= (1.0/Surf);

		for(int j=0; j<3; ++j) Xp[j] = qloc[nq][0]*Xst[0][j] 
 										     + qloc[nq][1]*Xst[1][j] 
										     + qloc[nq][2]*Xst[2][j];

		if(cs && cs->getPhi(nSt, qloc[nq][0], qloc[nq][1]) < 0.0) continue; 
	
		ElemForceCalcValid myObj;
		Elem* E = myTree->search<&Elem::isPointInside, ElemForceCalcValid,
										 &ElemForceCalcValid::Valid>(&myObj, X, Xp);

		if(!E) continue;
		
		for(int i=0; i<4; i++) T[i] = (*E)[i];

		Vec3D XO[4]; 
		for(int i=0; i<4; i++)
			for(int j=0; j<3; ++j) XO[i][j] = X[T[i]][j];
		
		Vec3D bary;
		E->computeBarycentricCoordinates(X, Xp, bary); 
		if(bary[0] < 0.0 || bary[1] < 0.0 || bary[2] < 0.0 || bary[0]+bary[1]+bary[2] > 1.0) 
		{
			E = 0;
			continue;
		}

		Vec3D dbary[4];
		dbary[0] = Vec3D(1.0-bary[0],bary[1],bary[2]);
		dbary[1] = Vec3D(bary[0],1.0-bary[1],bary[2]);
		dbary[2] = Vec3D(bary[0],bary[1],1.0-bary[2]);
		dbary[3] = Vec3D(bary[0],bary[1],bary[2]);

		int norm[4] = {0, 0, 0, 0};

		for(int e=0; e<6; ++e)
		{
			int l = E->edgeNum(e);
			
			if(LSS.edgeIntersectsStructure(0,l)) 
			{
				int i = E->edgeEnd(e, 0);
				int j = E->edgeEnd(e, 1);

				LevelSetResult lsResij = LSS.getLevelSetDataAtEdgeCenter(0.0, l, (T[i]<T[j]));

				norm[i] = (lsResij.gradPhi*(Xstruct[lsResij.trNodes[0]] - XO[i]) <= 0) ? -1 : 1;

				LevelSetResult lsResji = LSS.getLevelSetDataAtEdgeCenter(0.0, l, (T[i]>=T[j]));
				
				norm[j] = (lsResji.gradPhi*(Xstruct[lsResji.trNodes[0]] - XO[j]) <= 0) ? -1 : 1;
			}				
		}

		double mindist[2] = {FLT_MAX, FLT_MAX};

		int node[2] = {-1, -1};

		Vec3D stN[2] = {-stNormal, stNormal};

		for (int i=0; i<4; i++) 
		{
			double dist = dbary[i].norm();
				
			if(norm[i] < 0)
			{
				if(!LSS.isOccluded(0.0, T[i]) && dist < mindist[0] && stNormal*(Xp - XO[i]) <= 0.0 ) 
				{
					mindist[0] = dist;
					node[0] = T[i];
				}
			}
			else if(norm[i] > 0) 
			{
				if(!LSS.isOccluded(0.0, T[i]) && dist < mindist[1] && stNormal*(Xp - XO[i]) > 0.0 ) 
				{
					mindist[1] = dist;
					node[1] = T[i];
				}
			}
		}
			
		for(int dir=0; dir<2; ++dir)
		{
			if(node[dir] < 0) continue;

			stNodeDir[nSt][dir] = 1;

			Vec3D Xn;
			for(int i=0; i<3; ++i) Xn[i] = X[node[dir]][i];

			double h = (Xn - Xp)*stNormal;

			Vec3D Xe;
			if(!LSS.isActive(0.0, node[dir]))
				Xe = Xn + stNormal*h*dtol_plus;
			else
				Xe = Xp;
			 
			if(dir == 0) for(int i=0; i<3; ++i) stX1[nSt][i] = Xe[i];
			else         for(int i=0; i<3; ++i) stX2[nSt][i] = Xe[i];
		}
	}
}

//-----------------------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeEMBNodeScalarQuantity_step2(SVec<double,3> &X, SVec<double,dim> &V, 
																	PostFcn *postFcn, VarFcn *varFcn, 
																	Vec<int> &fluidId,
																	double (*Qnty)[3], int sizeQnty, int numStructElems, int (*stElem)[3],
																	Vec<Vec3D>& Xstruct, LevelSetStructure &LSS,
																	double pInfty, 
																	Vec<GhostPoint<dim>*> *ghostPoints,
																	NodalGrad<dim, double> &ngrad, 
	                                                int** stNodeDir, double** stX1, double** stX2)
{

	const int qOrder = 1;
	Quadrature quadrature_formula(qOrder);

	int nqPoint = quadrature_formula.n_point;
	const int nq = 0;

	double (*qloc)[3](quadrature_formula.qloc);
	double *qweight(quadrature_formula.weight);

	int T[4];

	Vec3D DX, gradP, flocal;

	SVec<double,dim> gradX = ngrad.getX();
	SVec<double,dim> gradY = ngrad.getY();
	SVec<double,dim> gradZ = ngrad.getZ();

	CrackingSurface* cs = LSS.getCrackingSurface();

	double Vext[dim]; 
	double *Vtet[4];
	int stNode[3];
	Vec3D Xst[3];
	Vec3D Xp;

	double Surf;

	double dp1dxj[4][3]; 

	for(int i=0;i<4;++i) 
		for(int j=0;j<3;++j) dp1dxj[i][j] = 0.0;

	for(int nSt=0; nSt<numStructElems; ++nSt) 
	{
		for(int j=0; j<3; ++j) 
		{
			stNode[j] = stElem[nSt][j];
            Xst[j] = Xstruct[stNode[j]]; 
		}

		Vec3D stNormal = 0.5*(Xst[1]-Xst[0])^(Xst[2]-Xst[0]);
		
		Surf = sqrt(stNormal*stNormal);
		
		if(Surf != 0) stNormal *= (1.0/Surf);

		for(int j=0; j<3; ++j) Xp[j] = qloc[nq][0]*Xst[0][j] 
 										     + qloc[nq][1]*Xst[1][j] 
										     + qloc[nq][2]*Xst[2][j];

		if(cs && cs->getPhi(nSt, qloc[nq][0], qloc[nq][1]) < 0.0) continue;

		int node_e[2] = {-1, -1};
		//Elem* elem_e[2];
		Elem* elem_e1;
		Elem* elem_e2;

		double mindist[2] = {FLT_MAX, FLT_MAX};

		for(int dir=0; dir<2; ++dir)
		{
			if(stNodeDir[nSt][dir] == 0) continue;
				
			Vec3D Xe;
			if(dir == 0) for(int i=0; i<3; ++i) Xe[i] = stX1[nSt][i];
			else         for(int i=0; i<3; ++i) Xe[i] = stX2[nSt][i];

			ElemForceCalcValid myObj;
			Elem* En = myTree->search<&Elem::isPointInside, ElemForceCalcValid,
											  &ElemForceCalcValid::Valid>(&myObj, X, Xe);
						
			if(!En) continue;

			if(dir == 0) elem_e1 = En;
			else         elem_e2 = En;

			for(int k=0; k<4; ++k) 
			{
				int Nk = (*En)[k];
					
				if(LSS.isActive(0.0, Nk))
				{
					Vec3D Xk;
					for(int i=0; i<3; ++i) Xk[i] = X[Nk][i];
						
					Vec3D vDist = Xp - Xk;
					double dist = sqrt(vDist*vDist);
						
					if(dist < mindist[dir])
					{
						mindist[dir] = dist;
						node_e[dir] = Nk;
					}
				}
			}
		}

		// ------------------------------------------------------ 

		double Cplocal = 0.0;
		double Cflocal = 0.0;

		Vec3D stN[2] = {-stNormal*Surf, stNormal*Surf};

		for(int dir=0; dir<2; ++dir)
		{	
			int Ni = node_e[dir];			

			if(Ni < 0) continue;

			// ------------------ Inviscid part  ------------------ //
			for(int j=0; j<3; ++j) DX[j] = Xp[j] - X[Ni][j];
					
			for(int k = 0; k < dim; ++k) Vext[k] = V[Ni][k]
													       + gradX[Ni][k]*DX[0] 
													       + gradY[Ni][k]*DX[1] 
															 + gradZ[Ni][k]*DX[2];

			int fid;
			fid = fluidId[Ni] ? fluidId[Ni]:0;	

			double pp = postFcn->computeNodeScalarQuantity(PostFcn::PRESSURECOEFFICIENT, Vext, Xp, fid, NULL);
			
			Cplocal += pp;
			// ---------------------------------------------------- //

			// ------------------- Viscous part  ------------------ //
			if(ghostPoints)
			{
				GhostPoint<dim> *gp;				

				Elem* elem_tmp;
				if(dir == 0) elem_tmp = elem_e1;
				else         elem_tmp = elem_e2;

				elem_tmp->computeGradientP1Function(X, dp1dxj);

				double dist_tmp = 1.0e16;
				int N_e = -1, ne = -1;

				for(int k=0; k<4; ++k)
				{
					int Nk = (*elem_tmp)[k];

					Vec3D Xwall;
					bool isNGhost = LSS.xWallNode(Nk, Xwall);

					Vec3D X_da, X_db;
					for(int j=0; j<3; ++j)
					{
						X_da[j] = X[Nk][j] - Xp[j];
						X_db[j] = X[Ni][j] - Xp[j];
					}

					double dadb = X_da * X_db;
					double dist_x = sqrt(X_da*X_da);

					if(dist_x < dist_tmp && dadb > 0.0)  
					{
						dist_tmp = dist_x;
						N_e = Nk;
						ne  = k;
					}

					if(isNGhost)
					{						
						gp = (*ghostPoints)[Nk];						

						Vec3D Xiw, Xkw;
						for(int j=0; j<3; ++j) 
						{
							Xiw[j] = X[Ni][j] - Xwall[j];
							Xkw[j] = X[Nk][j] - Xwall[j];
						}

						double d1 = Xiw*stN[dir];
						double d2 = Xkw*stN[dir];
							
						int df = (d1*d2 > 0.0) ? 1 : -1;

						Vtet[k] = gp->getPrimitiveState(df);
					}
					else 
						Vtet[k] = V[Nk];
				}

				Vec3D X_ = X[N_e] - Xp;
				double dw = sqrt(X_ * X_);
				
				// Dummy values (not used)
				double d2w[3] = {0, 0, 0}; 
				//double *Vwall = 0;
				double *Vface[3] = {0,0,0};
				
				d2w[0] = dw; d2w[1] = dw; d2w[2] = dw;
				Vface[0] = Vtet[ne]; Vface[1] = Vtet[ne]; Vface[2] = Vtet[ne];

				double* Vwall = Vtet[ne];
				Vwall[1] = Vwall[2] = Vwall[3] = 0.0;

				Vec3D F = postFcn->computeViscousForce(dp1dxj, stN[dir], d2w, Vwall, Vface, Vtet);
				
				Vec3D tdir(1.0, 0.0, 0.0);
					
				Cflocal += 2.0 * tdir * F / Surf;
			}
			// ---------------------------------------------------- //

			Qnty[stNode[0]][0] += Surf;
			Qnty[stNode[1]][0] += Surf;
			Qnty[stNode[2]][0] += Surf;
 
			Qnty[stNode[0]][1] += Cplocal*Surf;
			Qnty[stNode[1]][1] += Cplocal*Surf;
			Qnty[stNode[2]][1] += Cplocal*Surf;
 
			Qnty[stNode[0]][2] += Cflocal*Surf;
			Qnty[stNode[1]][2] += Cflocal*Surf;
			Qnty[stNode[2]][2] += Cflocal*Surf;
		} 
	}

}
//-----------------------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeRecSurfBasedForceLoad(int forceApp, int order, SVec<double,3> &X,
                                             double (*Fs)[3], int sizeFs, LevelSetStructure &LSS, double pInfty, 
                                             SVec<double,dim> &Wstarij, SVec<double,dim> &Wstarji, SVec<double,dim> &V, 
                                             Vec<GhostPoint<dim>*> *ghostPoints, PostFcn *postFcn, VarFcn* vf, Vec<int>* fid)
{
  if (forceApp!=3) 
    {fprintf(stderr,"ERROR: force method (%d) not recognized! Abort..\n", forceApp); exit(-1);}

  // ---------------
  //   Preparation
  // ---------------
  int T[4];       //nodes in a tet.

  double d2w[3]; // not used, but required by postFcn->computeForces(...)..
  int fid_local[4];
  double dPdx[3][3];  // not used, but required by postFcn->computeForces(...)..
  double dp1dxj[4][3]; // Gradient of the P1 basis functions
  d2w[0] = d2w[1] = d2w[2] = 0.0;
  for(int i=0;i<3;++i) for(int j=0;j<3;++j) dPdx[i][j] = 0.0;
  for(int i=0;i<4;++i) for(int j=0;j<3;++j) dp1dxj[i][j] = 0.0;

  // -------------------------------------------------------------------------
  //   Loop through elements. Reconstruct the interface and compute the force 
  // -------------------------------------------------------------------------
  for (int iElem=0; iElem<elems.size(); iElem++) {
    map<int,int> global2local;
    for (int i=0; i<4; i++) {
      T[i] = elems[iElem][i];
      global2local[elems[iElem][i]] = i;
    }

    if(ghostPoints) //dp1dxj is needed only for the computation of viscous force.
      elems[iElem].computeGradientP1Function(X, dp1dxj);
    PolygonReconstructionData polygons[4];

    // find polygons.
    int numberOfPolygons = getPolygons(elems[iElem], LSS, polygons);
    assert(numberOfPolygons<=4);
    int index,vertex1,vertex2;

    // loop through each polygon
    for(int i=0; i < numberOfPolygons; ++i){

      PolygonReconstructionData& polygon=polygons[i];
      int nEdges = polygon.numberOfEdges;

      // get intersection information.
      std::vector<LevelSetResult> lsRes(nEdges);
      if(!nEdges) continue; //KW: why need this?
      for (int k=0; k<nEdges; ++k) {
        lsRes[k] = LSS.getLevelSetDataAtEdgeCenter(0,polygon.edge[k],polygon.edgeWithVertex[k][0]<polygon.edgeWithVertex[k][1]);
        if (lsRes[k].alpha<0) {fprintf(stderr,"Unable to get intersection results at edge center! Abort...\n"); exit(-1);}
      }

      // determine if we are applying a force from the "real fluid" or "ghost" side. From the "ghost" side, only a pressure
      // force (using pInfty) is applied, even for viscous flows.
      bool applyRealForce = true;
      for(int k=0; k<nEdges; k++)
        if(!LSS.isActive(0,polygon.edgeWithVertex[k][0])) {
          applyRealForce = false;
          break;
        }

      // get the correct states (mix of real and ghost) for this tet. This info is needed for viscous simulation only.
      double *v[4];
      if(ghostPoints && applyRealForce) {
        for(int k=0; k<4; k++)
          v[k] = V[T[k]];

        GhostPoint<dim> *gp;
        for(int k=0; k<nEdges; k++) {
          int vertex2 = polygon.edgeWithVertex[k][1];
          gp = ghostPoints->operator[](vertex2);
          v[global2local[vertex2]] = gp->getPrimitiveState(); 
        }
      }

      // get information at intersection points (3 for triangles, 4 for quadrangles).
      std::vector<Vec3D> Xinter(nEdges);
      double Vinter[nEdges][dim];
      Vec3D start_vertex;for(int m=0;m<3;++m) start_vertex[m]=X[polygon.nodeToLookFrom][m];
      for(int k=0; k<nEdges; ++k) {
        vertex1 = polygon.edgeWithVertex[k][0];
        vertex2 = polygon.edgeWithVertex[k][1];
        int l = polygon.edge[k];
        double alpha = lsRes[k].alpha;
        for(int m=0; m<3; m++)
          Xinter[k][m] = alpha*X[vertex1][m] + (1.0-alpha)*X[vertex2][m];

        fid_local[k] = fid?(*fid)[vertex1]:0;
        if(!applyRealForce)
          Vinter[k][4] = pInfty; //this is enough even for Tait (look at postFcn->computeEmbeddedForce(...))
        else {
          if(ghostPoints)  //viscous
            for(int m=0; m<dim; m++)
              Vinter[k][m] = alpha*v[global2local[vertex1]][m] + (1.0-alpha)*v[global2local[vertex2]][m];
          else { // inviscid
            if(vertex1<vertex2) // apply Wstarij
              for(int m=0; m<dim; m++)
                Vinter[k][m] = Wstarij[l][m];
            else //apply Wstarji
              for(int m=0; m<dim; m++)
                Vinter[k][m] = Wstarji[l][m];
          }
        }
      } 

      // compute force...
      double dist13,dist02;
      double fac1,fac2;
      double *Xface[3], *Vface[3];
      int fid_face[3];
      double oneThird = 1.0/3.0;
      Vec3D nf;
      Vec3D fi0,fi1,fi2,fv,FNodal; // forces storage

      switch(nEdges){
        case 3: //got a triangle
          nf = 0.5*(Xinter[1]-Xinter[0])^(Xinter[2]-Xinter[0]);
          if(nf*(Xinter[1]-start_vertex) <= 0) nf *=-1;
          for(int m=0; m<3; m++) {
            Xface[m] = &Xinter[m][0];
            Vface[m] = &Vinter[m][0];
            fid_face[m] = fid_local[m]; 
          }
          postFcn->computeForceEmbedded(order, dp1dxj, 
					Xface, nf, d2w,
					0/*Vwall*/, Vface, v, pInfty,
                                        fi0, fi1, fi2, fv, 
					dPdx, 0/*"hydro"*/, 
					fid_face, applyRealForce);

          if(ghostPoints && applyRealForce) {
            fv *= oneThird;
            fi0+=fv;
            fi1+=fv;
            fi2+=fv;
          }
          sendLocalForce(fi0,lsRes[0],Fs);
          sendLocalForce(fi1,lsRes[1],Fs);
          sendLocalForce(fi2,lsRes[2],Fs);
          break;
        case 4: //got a quadrangle. cut it into two triangles.
          //dist02 = (Xinter[2]-Xinter[0]).norm();
          //dist13 = (Xinter[3]-Xinter[1]).norm();
          fac1 = lsRes[0].gradPhi*lsRes[2].gradPhi;                                                    
          fac2 = lsRes[1].gradPhi*lsRes[3].gradPhi;
          if(fac1 > fac2){ // connect 0,2.
            //for triangle 012
            nf = 0.5*(Xinter[1]-Xinter[0])^(Xinter[2]-Xinter[0]);
            if(nf*(Xinter[1]-start_vertex) <= 0) nf *=-1;
            Xface[0] = &Xinter[0][0];  Vface[0] = &Vinter[0][0];  fid_face[0] = fid_local[0];
            Xface[1] = &Xinter[1][0];  Vface[1] = &Vinter[1][0];  fid_face[1] = fid_local[1];
            Xface[2] = &Xinter[2][0];  Vface[2] = &Vinter[2][0];  fid_face[2] = fid_local[2];

            postFcn->computeForceEmbedded(order,dp1dxj,Xface,nf,d2w,0/*Vwall*/,Vface,v, pInfty,
                                          fi0,fi1,fi2,fv,dPdx,0/*"hydro"*/, fid_face, applyRealForce);

            if(ghostPoints && applyRealForce) {
              fv *= oneThird;
              fi0+=fv;
              fi1+=fv;
              fi2+=fv;
            }
            sendLocalForce(fi0,lsRes[0],Fs);
            sendLocalForce(fi1,lsRes[1],Fs);
            sendLocalForce(fi2,lsRes[2],Fs);

            //for triangle 023
            nf = 0.5*(Xinter[2]-Xinter[0])^(Xinter[3]-Xinter[0]);
            if(nf*(Xinter[1]-start_vertex) <= 0) nf *=-1;
            Xface[0] = &Xinter[0][0];  Vface[0] = &Vinter[0][0];  fid_face[0] = fid_local[0];
            Xface[1] = &Xinter[2][0];  Vface[1] = &Vinter[2][0];  fid_face[1] = fid_local[2];
            Xface[2] = &Xinter[3][0];  Vface[2] = &Vinter[3][0];  fid_face[2] = fid_local[3];

            postFcn->computeForceEmbedded(order,dp1dxj,Xface,nf,d2w,0/*Vwall*/,Vface,v, pInfty,
                                          fi0,fi1,fi2,fv,dPdx,0/*"hydro"*/, fid_face, applyRealForce);

            if(ghostPoints && applyRealForce) {
              fv *= oneThird;
              fi0+=fv;
              fi1+=fv;
              fi2+=fv;
            }
            sendLocalForce(fi0,lsRes[0],Fs);
            sendLocalForce(fi1,lsRes[2],Fs);
            sendLocalForce(fi2,lsRes[3],Fs);

          }else{ // connect 1,3.
            //for triangle 123
            nf = 0.5*(Xinter[2]-Xinter[1])^(Xinter[3]-Xinter[1]);
            if(nf*(Xinter[1]-start_vertex) <= 0) nf *=-1;
            Xface[0] = &Xinter[1][0];  Vface[0] = &Vinter[1][0];  fid_face[0] = fid_local[1]; 
            Xface[1] = &Xinter[2][0];  Vface[1] = &Vinter[2][0];  fid_face[1] = fid_local[2];
            Xface[2] = &Xinter[3][0];  Vface[2] = &Vinter[3][0];  fid_face[2] = fid_local[3];

            postFcn->computeForceEmbedded(order,dp1dxj,Xface,nf,d2w,0/*Vwall*/,Vface,v, pInfty,
                                          fi0,fi1,fi2,fv,dPdx,0/*"hydro"*/, fid_face, applyRealForce);

            if(ghostPoints && applyRealForce) {
              fv *= oneThird;
              fi0+=fv;
              fi1+=fv;
              fi2+=fv;
            }
            sendLocalForce(fi0,lsRes[1],Fs);
            sendLocalForce(fi1,lsRes[2],Fs);
            sendLocalForce(fi2,lsRes[3],Fs);

            //for triangle 013
            nf = 0.5*(Xinter[1]-Xinter[0])^(Xinter[3]-Xinter[0]);
            if(nf*(Xinter[1]-start_vertex) <= 0) nf *=-1;
            Xface[0] = &Xinter[0][0];  Vface[0] = &Vinter[0][0];  fid_face[0] = fid_local[0];
            Xface[1] = &Xinter[1][0];  Vface[1] = &Vinter[1][0];  fid_face[1] = fid_local[1];
            Xface[2] = &Xinter[3][0];  Vface[2] = &Vinter[3][0];  fid_face[2] = fid_local[3];

            postFcn->computeForceEmbedded(order,dp1dxj,Xface,nf,d2w,0/*Vwall*/,Vface,v, pInfty,
                                          fi0,fi1,fi2,fv,dPdx,0/*"hydro"*/, fid_face, applyRealForce);

            if(ghostPoints && applyRealForce) {
              fv *= oneThird;
              fi0+=fv;
              fi1+=fv;
              fi2+=fv;
            }
            sendLocalForce(fi0,lsRes[0],Fs);
            sendLocalForce(fi1,lsRes[1],Fs);
            sendLocalForce(fi2,lsRes[3],Fs);
          }

          break;
#if 0
        default: fprintf(stderr,"ANOTHER PROBLEM!!!\n");break;
#endif
        default: break;
      }
    }
  }
}
//-----------------------------------------------------------------------------------------------

template<int dim>
void SubDomain::computeCVBasedForceLoad(int forceApp, int orderOfAccuracy, GeoState& geoState,
                                        SVec<double,3> &X, double (*Fs)[3], int sizeFs,
                                        LevelSetStructure &LSS, double pInfty, 
                                        SVec<double,dim> &Wstarij, SVec<double,dim> &Wstarji,
                                        SVec<double,dim> &V, Vec<GhostPoint<dim>*> *ghostPoints,
                                        PostFcn *postFcn, NodalGrad<dim, double> &ngrad, VarFcn *vf, Vec<int>* fid)
{
  if (forceApp!=1) 
	{
		fprintf(stderr,"ERROR: force method (%d) not recognized! Abort..\n", forceApp); 
		exit(-1);
	}

  Vec<Vec3D>& normal = geoState.getEdgeNormal();
  bool* masterFlag = edges.getMasterFlag();
  int (*ptr)[2];
	ptr = edges.getPtr();

  int i,j;

  bool iActive, jActive, intersect;

  Vec3D vectorIJ, gradP;
  SVec<double,dim> gradX = ngrad.getX();
  SVec<double,dim> gradY = ngrad.getY();
  SVec<double,dim> gradZ = ngrad.getZ();
  double dudxj[3][3];

	for(int l=0; l<edges.size(); l++) 
	{
    if (!masterFlag[l]) continue;
    i = ptr[l][0];
    j = ptr[l][1];
		
		if(LSS.withCracking()) 
		{
      iActive = !LSS.isOccluded(0,i);
      jActive = !LSS.isOccluded(0,j);
		} 
		else 
		{
      iActive = LSS.isActive(0,i);
      jActive = LSS.isActive(0,j);
    }
    intersect = LSS.edgeIntersectsStructure(0,l);

      //both inside structure
		if(!iActive && !jActive) continue; 

    if (!intersect) continue; 

    double *v[2] = {V[i],V[j]};  

		if(iActive) 
		{
      Vec3D flocal(0.0,0.0,0.0); 

      LevelSetResult lsRes = LSS.getLevelSetDataAtEdgeCenter(0.0, l, true);

			if(ghostPoints)
			{
				// Viscous Simulation
        // Replace the state of node j by its corresponding ghost state
        /* This is useless. We use ngrad to compute the velocity gradient. 
        GhostPoint<dim> *gp;
        gp = ghostPoints->operator[](j);
        v[1] = gp->getPrimitiveState();
        */
				for(int m=0;m<3;++m) vectorIJ[m] = X[j][m] - X[i][m];

        gradP[0] = gradX[i][4];
        gradP[1] = gradY[i][4];
        gradP[2] = gradZ[i][4];

        double pp = vf->getPressure(v[0], fid?(*fid)[i]:0);

        flocal = (pp - pInfty + 0.5*(gradP*vectorIJ))*normal[l];

        // Compute Velocity Gradient

        // ***************** Method 1 **********************
        dudxj[0][0] = gradX[i][1]; dudxj[0][1] = gradY[i][1]; dudxj[0][2] = gradZ[i][1];
        dudxj[1][0] = gradX[i][2]; dudxj[1][1] = gradY[i][2]; dudxj[1][2] = gradZ[i][2];
        dudxj[2][0] = gradX[i][3]; dudxj[2][1] = gradY[i][3]; dudxj[2][2] = gradZ[i][3];

        /*
        // ***************** Method 2 **********************
        double norm = vectorIJ.norm();
        double deltaU = (v[0][1]-v[1][1])/(norm*norm);
        dudxj[0][0] = deltaU*vectorIJ[0];dudxj[0][1] = deltaU*vectorIJ[1];dudxj[0][2] = deltaU*vectorIJ[2];
        double deltaU = (v[0][2]-v[1][2])/(norm*norm);
        dudxj[1][0] = deltaU*vectorIJ[0];dudxj[1][1] = deltaU*vectorIJ[1];dudxj[1][2] = deltaU*vectorIJ[2];
        double deltaU = (v[0][3]-v[1][3])/(norm*norm);
        dudxj[2][0] = deltaU*vectorIJ[0];dudxj[2][1] = deltaU*vectorIJ[1];dudxj[2][2] = deltaU*vectorIJ[2];
        */

        flocal += postFcn->computeViscousForceCVBoundary(normal[l],v[0],dudxj);

			} 
			else 
			{
        double pp = vf->getPressure(Wstarij[l],fid?(*fid)[i]:0);
        flocal = (pp - pInfty)*normal[l];
      }

      sendLocalForce(flocal, lsRes, Fs);
		}

		if(jActive) 
		{
      Vec3D flocal(0.0,0.0,0.0);
      LevelSetResult lsRes = LSS.getLevelSetDataAtEdgeCenter(0.0, l, false);

			if(ghostPoints) 
			{
            // Viscous Simulation
        // Replace the state of node i by its corresponding ghost state
        /* This is useless. We use ngrad to compute the velocity gradient.
        GhostPoint<dim> *gp;
        gp = ghostPoints->operator[](i);
        v[0] = gp->getPrimitiveState();
        */
				for(int m=0;m<3;++m) vectorIJ[m] = X[j][m] - X[i][m];

        gradP[0] = gradX[j][4];
        gradP[1] = gradY[j][4];
        gradP[2] = gradZ[j][4];

        double pp = vf->getPressure(v[1], fid?(*fid)[j]:0);

        // Minus, cause the normal points toward j
        flocal = -(pp - pInfty - 0.5*(gradP*vectorIJ))*normal[l];

        // Compute Velocity Gradient

        // ***************** Method 1 **********************
        dudxj[0][0] = gradX[j][1]; dudxj[0][1] = gradY[j][1]; dudxj[0][2] = gradZ[j][1];
        dudxj[1][0] = gradX[j][2]; dudxj[1][1] = gradY[j][2]; dudxj[1][2] = gradZ[j][2];
        dudxj[2][0] = gradX[j][3]; dudxj[2][1] = gradY[j][3]; dudxj[2][2] = gradZ[j][3];

        // Minus, cause the normal points toward j
        flocal -= postFcn->computeViscousForceCVBoundary(normal[l],v[1],dudxj);
			} 
			else 
			{
        double pp = vf->getPressure(Wstarji[l],fid?(*fid)[j]:0);
        flocal = -(pp - pInfty)*normal[l];
      }
      sendLocalForce(flocal, lsRes, Fs);
    }
  }
}


//------------------------------------------------------------------------------

template<int dim>
void SubDomain::blur(SVec<double,dim> &U,SVec<double,dim> &U0, Vec<double>& weight)
{
  const Connectivity &nToN = *getNodeToNode(); 
  for(int currentNode=0;currentNode<numNodes();++currentNode) {
        
    for (int k = 0; k < dim; ++k) {
      U0[currentNode][k] = 0.0;
    }

    for(int j=0;j<nToN.num(currentNode);++j){
      int neighborNode=nToN[currentNode][j];
      for (int k = 0; k < dim; ++k) {
	U0[currentNode][k] += U[neighborNode][k];
      }
      
    }
    weight[currentNode] = (double)nToN.num(currentNode);
  }
}

//------------------------------------------------------------------------------

template<int dimLS>
void SubDomain::updateFluidIdFS2(LevelSetStructure &LSS, SVec<double,dimLS> &PhiV, SVec<bool,3> &poll, 
                                 Vec<int> &fluidId, bool *masterFlag)
{
  fprintf(stderr,"ERROR: This function should not be called anymore. Moved into FluidSelector.\n");

  const Connectivity &Node2Node = *getNodeToNode();
  // ------- Determine status for grid-points swept by FS interface -------
  // Rule No.1: If this grid point is "occluded", set its status to "numPhases".
  // Rule No.2: If its "visible && !occluded && !swept" neighbors have the same status, use this one. 
  // Rule No.3: Otherwise, consider the sign of "PhiV". (PhiV should have been "blurred".)
  
  int rnk;
  MPI_Comm_rank(MPI_COMM_WORLD,&rnk);

  for(int i=0; i<PhiV.size(); i++) {
    bool swept = LSS.isSwept(0.0,i);
    bool occluded = LSS.isOccluded(0.0,i);

    if(rnk==6 && i==30903)
      fprintf(stderr,"Rank 6, i = %d, swept = %d, occluded = %d, dimLS = %d, poll = %d %d %d.\n", i, swept, occluded, dimLS, poll[i][0], poll[i][1], poll[i][2]);

    if(!swept) {//nothing to be done
      continue;

//KW: I DONT KNOW WHO ADDED THIS, NOR WHY
//      if(!occluded && fluidId[i]!=dimLS+1) //this "if" is false when the structural elment covering node i got deleted in Element Deletion.
//        continue;
    }

    if(occluded) { // Rule No.1
      fluidId[i] = LSS.numOfFluids();
      if(!poll[i][2/*LSS.numOfFluids()*/]) fprintf(stderr,"TOO BAD TOO!\n");
      continue;
    }

    int count = (int)poll[i][0] + (int)poll[i][1] + (int)poll[i][2];
    switch (count) {
      case 0: //no info
        fprintf(stderr,"More than one layer of nodes are swept in one step (near Node %d).\n", locToGlobNodeMap[i]+1);
        DebugTools::SpitRank();
        break;
      case 1: // Rule No.2
        if (poll[i][0]) fluidId[i] = 0;
        else if (poll[i][1]) fluidId[i] = dimLS;
        else if (poll[i][2]) fluidId[i] = dimLS+1;
        break;
    }
    
    // consider programmed burn

    if(count==1) //already applied Rule No.2
      continue;

    // Rule No.3
    bool done = false;
    for(int k=0; k<dimLS; k++)
      if(PhiV[i][k]>0.0) {
        fluidId[i] = k+1;
        done = true;
        break;
      }
    if(!done)
      fluidId[i] = 0;
  }
}

//------------------------------------------------------------------------------
/*
template<int dimLS> 
void SubDomain::updateFluidIdFS2(LevelSetStructure &LSS, SVec<double,dimLS> &PhiV, Vec<int> &fluidId, bool *masterFlag)
{
  // ------- Determine status for grid-points swept by FS interface -------
  // Rule No.1: If this grid point is "occluded", set its status to "numPhases".
  // Rule No.2: If its "visible && !occluded && !swept" neighbors have the same status, use this one. 
  // Rule No.3: Otherwise, consider the sign of "PhiV". (PhiV should have been "blurred".)

  const Connectivity &Node2Node = *getNodeToNode();
  
  for(int i=0; i<PhiV.size(); i++) {
    bool swept = LSS.isSwept(0.0,i);
    bool occluded = LSS.isOccluded(0.0,i);

    //DEBUG
    int myNode = 284121;
    if(locToGlobNodeMap[i]+1==myNode){
      int myrank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
      fprintf(stderr,"Node %d(%d), CPU = %d. master = %d, swept = %d, occluded = %d, id = %d, phi = %e.\n", myNode, i, myrank, masterFlag[i], swept, occluded, fluidId[i], PhiV[i][0]);
      for(int j=0; j<Node2Node.num(i); j++) { 
        if(Node2Node[i][j]==i) continue;
        fprintf(stderr,"  Nei(%d,%d) on CPU %d--> GlobId(%d), occluded(%d), swept(%d), intersect(%d), id(%d), phi(%e).\n", myNode,i,myrank,
                          locToGlobNodeMap[Node2Node[i][j]], LSS.isOccluded(0.0,Node2Node[i][j]), LSS.isSwept(0.0,Node2Node[i][j]),
                          LSS.edgeIntersectsStructure(0.0,i,Node2Node[i][j]), fluidId[Node2Node[i][j]], PhiV[Node2Node[i][j]][0]); 
      }

    }


    if(!swept) // nothing to be done.
      continue;
    if(!masterFlag[i]) {//will get Id from another subdomain (TODO: May need something more sophisticated... (KW).)
      fluidId[i] = 0;
      continue;
    }

    // Rule No.1
    if(occluded) {
      fluidId[i] = LSS.numOfFluids();
      continue;
    }

    // Rule No.2
    int myId = -1;
    int iNei, count = 0;
    bool consistent = false;
    for(int j=0; j<Node2Node.num(i); j++) {
      iNei = Node2Node[i][j];
      if(i==iNei)
        continue;
      if(LSS.isOccluded(0.0,iNei) || LSS.isSwept(0.0,iNei) || LSS.edgeIntersectsStructure(0.0,i,iNei))
        continue;
      count++;
      if(myId==-1) {
        myId = fluidId[iNei];
        consistent = true;
      } else if(myId!=fluidId[iNei]) {
        consistent = false;
        break;
      }
    }
    if(count==0)
      fprintf(stderr,"WARNING: More than one layer of nodes are swept in one step (near Node %d).\n", locToGlobNodeMap[i]+1);

    if(consistent) {
      fluidId[i] = myId;
      continue;
    }

    // Rule No.3
    bool done = false;
    for(int k=0; k<dimLS; k++)
      if(PhiV[i][k]>0.0) {
        fluidId[i] = k+1;
        done = true;
      }
    if(!done)
      fluidId[i] = 0;
  }
}
*/

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void SubDomain::debugMultiPhysics(LevelSetStructure &LSS, SVec<double,dimLS> &PhiV, Vec<int> &fluidId, SVec<double,dim> &U)
{
  for(int i=0; i<nodes.size(); i++) {
    switch (fluidId[i]) {
      case 0:
        if(PhiV[i][dimLS-1]>=0.0)
          fprintf(stderr,"BUGGY: Sub %d, Node %d: Id = %d but PhiV = %e.\n", globSubNum, locToGlobNodeMap[i]+1, fluidId[i], PhiV[i][dimLS-1]);
        if(LSS.isOccluded(0,i))
          fprintf(stderr,"BUGGY: Sub %d, Node %d: Id = %d, PhiV = %e,  but occluded!\n", globSubNum, locToGlobNodeMap[i]+1, 
                          fluidId[i], PhiV[i][dimLS-1]);
        break;
      case 1:
        if(PhiV[i][dimLS-1]<=0.0)
          fprintf(stderr,"BUGGY: Sub %d, Node %d: Id = %d but PhiV = %e.\n", globSubNum, locToGlobNodeMap[i]+1, fluidId[i], PhiV[i][dimLS-1]);
        if(LSS.isOccluded(0,i))
          fprintf(stderr,"BUGGY: Sub %d, Node %d: Id = %d, PhiV = %e,  but occluded!\n", globSubNum, locToGlobNodeMap[i]+1, 
                          fluidId[i], PhiV[i][dimLS-1]);
        break;
      case 2:
        if(fabs(PhiV[i][dimLS-1])>1.0e-8)
          fprintf(stderr,"BUGGY: Sub %d, Node %d: Id = %d but PhiV = %e.\n", globSubNum, locToGlobNodeMap[i]+1, fluidId[i], PhiV[i][dimLS-1]);
        if(!LSS.isOccluded(0,i))
          fprintf(stderr,"BUGGY: Sub %d, Node %d: Id = %d, PhiV = %e,  but NOT occluded!\n", globSubNum, locToGlobNodeMap[i]+1, 
                          fluidId[i], PhiV[i][dimLS-1]);
        break;
      default:
        fprintf(stderr,"BUGGY: Sub %d, Node %d: Id = %d!\n",  globSubNum, locToGlobNodeMap[i]+1, fluidId[i]);
    }
  }

/*
  int debugNode1 = 1033424, debugNode2 = -1;
  for(int i=0; i<nodes.size(); i++)
    if(locToGlobNodeMap[i]+1==debugNode1 || locToGlobNodeMap[i]+1==debugNode2) {
      fprintf(stderr,"%d: Node %d: Id = %d, PhiV = %e, occluded(%d), swept(%d), d2wall(%e), U(%e,%e,%e,%e,%e).\n",
              globSubNum, locToGlobNodeMap[i]+1, fluidId[i], PhiV[i][0], LSS.isOccluded(0,i), LSS.isSwept(0,i), LSS.distToInterface(0,i),
              U[i][0], U[i][1], U[i][2], U[i][3], U[i][4]);
      for(int j=0; j<NodeToNode->num(i); j++) {
        int yId = (*NodeToNode)[i][j];
        if(yId==i) continue;
        fprintf(stderr,"%d:    Nei(%d)=%d: Id = %d, PhiV = %e, occluded(%d), swept(%d), d2wall(%e), X(%d).\n", globSubNum,
                locToGlobNodeMap[i]+1, locToGlobNodeMap[yId]+1, fluidId[yId], PhiV[yId][0], LSS.isOccluded(0,yId), LSS.isSwept(0,yId),
                LSS.distToInterface(0,yId), LSS.edgeIntersectsStructure(0,edges.findOnly(i,yId)));
        if((locToGlobNodeMap[i]+1==debugNode1 && locToGlobNodeMap[yId]+1==debugNode2) ||
           (locToGlobNodeMap[i]+1==debugNode2 && locToGlobNodeMap[yId]+1==debugNode1)) {
          LevelSetResult resij = LSS.getLevelSetDataAtEdgeCenter(0.0,edges.findOnly(i,yId),i<yId);
          fprintf(stderr,"** (%d,%d): alpha(%e), Tr(%d,%d,%d), Normal(%e,%e,%e), xi(%e,%e,%e).\n",locToGlobNodeMap[i]+1,locToGlobNodeMap[yId]+1,
                  resij.alpha, resij.trNodes[0], resij.trNodes[1], resij.trNodes[2], resij.gradPhi[0], resij.gradPhi[1], resij.gradPhi[2],
                  resij.xi[0], resij.xi[1], resij.xi[2]);
        }
      }
    }
*/
}

//------------------------------------------------------------------------------

template<int dim, class Obj>
void SubDomain::integrateFunction(Obj* obj,SVec<double,3> &X,SVec<double,dim>& V, void (Obj::*F)(int node, const double* loc,double* f),
				  int npt) 
{
  elems.integrateFunction(obj,X,V,F,npt);
}

template<int dim> 
void SubDomain::interpolateSolution(SVec<double,3>& X, SVec<double,dim>& U, 
                                    const std::vector<Vec3D>& locs, double (*sol)[dim],
                                    int* status,int* last,int* nid,
                                    LevelSetStructure* LSS, Vec<GhostPoint<dim>*>* ghostPoints,
                                    VarFcn *varFcn, bool assumeCache, Vec<int> *fluidId) {

  elems.interpolateSolution(X,U,locs,sol,status,last,LSS,ghostPoints,varFcn,
			    assumeCache);
  for (int i = 0; i < locs.size(); ++i) {
    if(!status[i]) continue;
    int eid = last[i],nn;
    Elem& E = elems[eid];
    double mindist = std::numeric_limits<double>::max(),dst;
    for (int j = 0; j < E.numNodes(); ++j) {
      nn = E.nodeNum(j);
      dst = sqrt((X[nn][0]-locs[i][0])*(X[nn][0]-locs[i][0])+
                 (X[nn][1]-locs[i][1])*(X[nn][1]-locs[i][1])+
                 (X[nn][2]-locs[i][2])*(X[nn][2]-locs[i][2]));
      if (dst < mindist) {
        mindist = dst;
        nid[i] = nn;
      }
    }
    if(fluidId) {
      int fid0 = (*fluidId)[E.nodeNum(0)];
      for(int j = 1; j < E.numNodes(); ++j) {
        if((*fluidId)[E.nodeNum(j)] != fid0) {
          for(int k=0; k<5; ++k) sol[i][k] = 0;
          status[i] = 0;
          break;
        }
      }
    }
  }
}

template<int dim>
void SubDomain::interpolatePhiSolution(SVec<double,3>& X, SVec<double,dim>& U,
				       const std::vector<Vec3D>& locs, double (*sol)[dim],
				       int* status,int* last,int* nid,
				       bool assumeCache) {


  elems.interpolateSolution(X,U,locs,sol,status,last,NULL,
			    (Vec<GhostPoint<dim>*>*)0, NULL,
			    assumeCache);
  for (int i = 0; i < locs.size(); ++i) {
    if(!status[i]) continue;
    int eid = last[i],nn;
    Elem& E = elems[eid];
    double mindist = std::numeric_limits<double>::max(),dst;
    for (int j = 0; j < E.numNodes(); ++j) {
      nn = E.nodeNum(j);
      dst = sqrt((X[nn][0]-locs[i][0])*(X[nn][0]-locs[i][0])+
                 (X[nn][1]-locs[i][1])*(X[nn][1]-locs[i][1])+
                 (X[nn][2]-locs[i][2])*(X[nn][2]-locs[i][2]));
      if (dst < mindist) {
        mindist = dst;
        nid[i] = nn;
      }
    }
  }
}


//------------------------------------------------------------------------------

// Functions to compute the error (that is, the difference between two state vectors)
template <int dim>
void SubDomain::computeL1Error(bool* nodeFlag,SVec<double,dim>& U, SVec<double,dim>& Uexact, Vec<double>& vol,double error[dim], LevelSetStructure* LSS) {

  for (int k = 0; k < dim; ++k)
    error[k] = 0.0;

  for(int i=0; i<nodes.size(); i++) {

    if (nodeFlag[i] && (!LSS || LSS->isActive(0.0,i))) {
      
      for (int k = 0; k < dim; ++k) {
	
	error[k] += fabs(U[i][k]-Uexact[i][k])*vol[i];
      }
    }
  }
}

// Functions to compute the error (that is, the difference between two state vectors)
template <int dim>
void SubDomain::computeL2Error(bool* nodeFlag,SVec<double,dim>& U, SVec<double,dim>& Uexact, Vec<double>& vol,double error[dim], LevelSetStructure* LSS) {

  for (int k = 0; k < dim; ++k)
    error[k] = 0.0;

  for(int i=0; i<nodes.size(); i++) {

    if (nodeFlag[i] && (!LSS || LSS->isActive(0.0,i))) {
      
      for (int k = 0; k < dim; ++k) {
	
	error[k] += pow(U[i][k]-Uexact[i][k],2.0)*vol[i];
      }
    }
  }
}

template <int dim>
void SubDomain::computeLInfError(bool* nodeFlag,SVec<double,dim>& U, SVec<double,dim>& Uexact, double error[dim], LevelSetStructure* LSS) {

  for (int k = 0; k < dim; ++k)
    error[k] = 0.0;

  for(int i=0; i<nodes.size(); i++) {

    if (nodeFlag[i] && (!LSS || LSS->isActive(0.0,i))) {
      
      for (int k = 0; k < dim; ++k) {
	
        if (fabs(U[i][k]-Uexact[i][k]) > 0.05) {

//          std::cout << locToGlobNodeMap[i] << " " <<  U[i][k] << " " << Uexact[i][k] << std::endl;
        }
	
	error[k] = max(error[k],fabs(U[i][k]-Uexact[i][k]));
      }
    }
  }
}
 
template <int dim>
void SubDomain::computeHHBoundaryTermResidual(BcData<dim> &bcData,SVec<double,dim> &U,Vec<double>& res, VarFcn* vf) {

  faces.computeHHBoundaryTermResidual(bcData,U,res, vf);
}

template<int dim, class Scalar, int neq>
void SubDomain::computeJacobianFiniteVolumeTermHH(FluxFcn **fluxFcn, BcData<dim> &bcData,
						  GeoState& geoState,
						  Vec<double> &ctrlVol,
						  SVec<double,dim> &U, 
						  GenMat<Scalar,neq> &A, VarFcn* vf) {

  faces.computeHHBoundaryTermJacobian(fluxFcn, bcData, U, geoState, A,vf);

  for (int i=0; i<faces.size(); ++i) {
    Face& F = faces[i];
    Scalar *Auh = A.getElemUH(i);
    double voli;
    for (int l = 0; l < F.numNodes(); ++l) {
      voli = 1.0 / ctrlVol[F[l]];
      for (int k=0; k<neq; ++k)
        Auh[l*neq+k] *= voli;
    }
  }

}

