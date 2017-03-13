#include <Face.h>

#include <FluxFcn.h>
#include <FemEquationTerm.h>
#include <BcData.h>
#include <BcDef.h>
#include <Elem.h>
#include <ExactRiemannSolver.h>
#include <GeoState.h>
#include <Vector3D.h>
#include <Vector.h>
#include <GenMatrix.h>
#include <LowMachPrec.h>
#include <LevelSet/LevelSetStructure.h>
#include <HigherOrderMultiFluid.h>

#include <cmath>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
#endif

//------------------------------------------------------------------------------

template<class NodeMap>
void Face::renumberNodes(NodeMap &nodemap)
{
  for (int j=0; j<numNodes(); ++j)
    nodeNum(j) = nodemap[ nodeNum(j) ];
}

//------------------------------------------------------------------------------

template<int dim>
void Face::assignFreeStreamValues2(SVec<double,dim> &Uin, SVec<double,dim> &Uout, double *U)
{

  // UH (07/2012)
  // The next test is to handle the cases where code is set
  // to BC_KIRCHHOFF_SURFACE (== 9)
  if ((code < BC_MIN_CODE) | (code > BC_MAX_CODE))
    return;

  int k, j;

  NOT_CORRECTED("Divide by numNodes? Or take surface into account ?");
  if (code == BC_INLET_MOVING || code == BC_INLET_FIXED ||
      code == BC_DIRECTSTATE_INLET_MOVING || code == BC_DIRECTSTATE_INLET_FIXED ||
      code == BC_MASSFLOW_INLET_MOVING || code == BC_MASSFLOW_INLET_FIXED) 
    for (k=0; k<dim; ++k) {
      for (j=0, U[k] = 0.0; j<numNodes(); ++j) 
	U[k] += Uin[nodeNum(j)][k];
      U[k] /= numNodes();
    }
  else if (code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED ||
           code == BC_DIRECTSTATE_OUTLET_MOVING || code == BC_DIRECTSTATE_OUTLET_FIXED ||
           code == BC_MASSFLOW_OUTLET_MOVING || code == BC_MASSFLOW_OUTLET_FIXED)
    for (k=0; k<dim; ++k) {
      for (j=0, U[k] = 0.0; j<numNodes(); ++j) 
	U[k] += Uout[nodeNum(j)][k];
      U[k] /= numNodes();
    }
  else
    for (k=0; k<dim; ++k)
      U[k] = 0.0;

}

//------------------------------------------------------------------------------

template<int dim>
void Face::assignFreeStreamValues(double *Uin, double *Uout, double *U)
{

  // UH (07/2012)
  // The next test is to handle the cases where code is set
  // to BC_KIRCHHOFF_SURFACE (== 9)
  if ((code < BC_MIN_CODE) | (code > BC_MAX_CODE))
    return;

  int k;

  if (code == BC_INLET_MOVING || code == BC_INLET_FIXED)
    for (k=0; k<dim; ++k)
      U[k] = Uin[k];
  else if (code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED)
    for (k=0; k<dim; ++k)
      U[k] = Uout[k];
  else
    for (k=0; k<dim; ++k)
      U[k] = 0.0;

}

//------------------------------------------------------------------------------

template<int dim>
void Face::assignPorousWallValues(SVec<double,dim> &Uin, double *U)
{

  if ((code < BC_MIN_CODE) | (code > BC_MAX_CODE))
    return;

  int k, j;

  NOT_CORRECTED("Divide by numNodes? Or take surface into account ?");
  if (code == BC_POROUS_WALL_MOVING || code == BC_POROUS_WALL_FIXED) {
    for (k=0; k<dim; ++k) {
      for (j=0, U[k] = 0.0; j<numNodes(); ++j) 
	U[k] += Uin[nodeNum(j)][k];
      U[k] /= numNodes();
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void Face::computeFaceBcValue(SVec<double,dim> &Unode, double *Uface)
{
  int j, k;

  NOT_CORRECTED("Divide by numNodes? Or take surface into account ?");

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED ||
      code == BC_ADIABATIC_WALL_MOVING || code == BC_ADIABATIC_WALL_FIXED)
    for (k=0; k<dim; ++k) {
      for (j=0, Uface[k] = 0.0; j<numNodes(); ++j) 
	Uface[k] += Unode[nodeNum(j)][k];
      Uface[k] /= numNodes();
    }
}

//------------------------------------------------------------------------------

template<int dim1, int dim2>
void Face::computeNodeBcValue(SVec<double,3> &X, double *Uface, SVec<double,dim2> &Unode)
{

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED ||
      code == BC_ADIABATIC_WALL_MOVING || code == BC_ADIABATIC_WALL_FIXED) {
    Vec3D n;
    computeNormal(X, n);
    double S = sqrt(n*n);

    for (int j=0; j<numNodes(); ++j) {
      Unode[ nodeNum(j) ][0] += S;
      for (int k=1; k<dim2; ++k) 
	Unode[ nodeNum(j) ][k] += S * Uface[dim1-dim2+k];
    }
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim1, int dim2>
void Face::computeDerivativeOfNodeBcValue(SVec<double,3> &X, SVec<double,3> &dX, double *Uface, double *dUface, SVec<double,dim2> &dUnode)
{

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED ||
      code == BC_ADIABATIC_WALL_MOVING || code == BC_ADIABATIC_WALL_FIXED) {

    Vec3D n;
    Vec3D dn;

    computeNormalAndDerivative(X, dX, n, dn);

    double S = sqrt(n*n);
    double dS = 1.0/(2.0*S) * (dn*n + n*dn);

    for (int j=0; j<numNodes(); ++j) {
      dUnode[ nodeNum(j) ][0] += dS;
      for (int k=1; k<dim2; ++k)
        dUnode[ nodeNum(j) ][k] += dS * Uface[dim1-dim2+k] + S * dUface[dim1-dim2+k];
    }
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Face::computeNodeBCsWallValues(SVec<double,3> &X, SVec<double,1> &dNormSA, double *dUfaceSA, SVec<double,dim> &dUnodeSA)
{

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED ||
      code == BC_ADIABATIC_WALL_MOVING || code == BC_ADIABATIC_WALL_FIXED) {
    Vec3D n;
    computeNormal(X, n);
    double S = sqrt(n*n);

    for (int j=0; j<numNodes(); ++j) {
      dNormSA[ nodeNum(j) ][0] += S;
      for (int k=0; k<dim; ++k) 
        dUnodeSA[ nodeNum(j) ][k] += S * dUfaceSA[k];
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void Face::computeTimeStep(VarFcn *varFcn, Vec<Vec3D> &normals, Vec<double> &normalVel,
			   SVec<double,dim> &V, Vec<double> &dt, 
			   TimeLowMachPrec &tprec, Vec<int> &fluidId)
{
  Vec3D normal = getNormal(normals);
  double S = sqrt(normal * normal);
  double invS = 1.0 / S;

  Vec3D n = invS * normal;
  double ndot = invS * getNormalVel(normalVel);

  NOT_CORRECTED("Divide by numNodes? Or take surface into account ?");

  for (int l=0; l<numNodes(); ++l) {
    Vec3D u = varFcn->getVelocity(V[ nodeNum(l) ], fluidId[nodeNum(l)]);
    double a = varFcn->computeSoundSpeed(V[ nodeNum(l) ], fluidId[nodeNum(l)]);
    double un = u * n - ndot;
    double locMach = varFcn->computeMachNumber(V[ nodeNum(l) ], fluidId[nodeNum(l)]);
    double locbeta = tprec.getBeta(locMach,true);
    
    double beta2 = locbeta * locbeta;
    double coeff1 = (1.0+beta2)*un;
    double coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);

    dt[ nodeNum(l) ] += min(0.5*(coeff1-coeff2), 0.0)* S/numNodes();
  }   
}

//------------------------------------------------------------------------------

template<int dim>
 void Face::computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn, 
									 Vec<Vec3D> &normals, Vec<double> &normalVel, 
									 SVec<double,3> &X, SVec<double,dim> &V, 
			   Vec<double> &idti, Vec<double> &idtv,
									 TimeLowMachPrec &tprec, LevelSetStructure *LSS)
{
  Vec3D normal = getNormal(normals);
  double S = sqrt(normal * normal);
  double invS = 1.0 / S;

  Vec3D n = invS * normal;
  double ndot = invS * getNormalVel(normalVel);

  NOT_CORRECTED("Divide by numNodes? Or take surface into account ?");

	 bool isValid = true;

	 if(LSS)
	 { 		 
		 for(int k=0; k<numNodes(); ++k)
		 {
			 int e = edgeNum(k);

			 int i = nodeNum(edgeEnd(k,0));
			 int j = nodeNum(edgeEnd(k,1));

			 if(LSS->edgeIntersectsStructure(0, e)) isValid = false;

			 if(!LSS->isActive(0,i) || !LSS->isActive(0,j)) isValid = false;
		 }
	 }

	 if(!isValid) return;

	 for(int l=0; l<numNodes(); ++l) 
	 {
    Vec3D u = varFcn->getVelocity(V[ nodeNum(l) ]);
    double a = varFcn->computeSoundSpeed(V[ nodeNum(l) ]);
    double un = u * n - ndot;

    // Low-Mach Preconditioner
    double locMach = varFcn->computeMachNumber(V[ nodeNum(l) ]);
    double locbeta = tprec.getBeta(locMach,true);
    double beta2 = locbeta * locbeta;
    double coeff1 = (1.0+beta2)*un;
    double coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);

    idti[ nodeNum(l) ] += min(0.5*(coeff1-coeff2), 0.0)* S/numNodes();

  }
    
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void Face::computeDerivativeOfTimeStep(FemEquationTerm *fet, VarFcn *varFcn, Vec<Vec3D>  &normals, Vec<Vec3D>  &dNormals, Vec<double> normalVel, Vec<double> dNormalVel,
			   SVec<double,3> &X, SVec<double,3> &dX, SVec<double,dim> &V, SVec<double,dim> &dV, 
			   Vec<double> &dIdti, Vec<double> &dIdtv, double dMach,
                           TimeLowMachPrec &tprec)
{

  Vec3D normal = getNormal(normals);
  double S = sqrt(normal * normal);
  Vec3D dNormal = getdNormal(dNormals);
  double dS = (normal * dNormal) / sqrt(normal * normal);
  double invS = 1.0 / S;
  double dInvS = - dS / (S*S);

  Vec3D n = invS * normal;
  Vec3D dn = dInvS * normal + invS * dNormal;
  double ndot = invS * getNormalVel(normalVel);
  double dndot = dInvS * getNormalVel(normalVel) + invS *getdNormalVel(dNormalVel);

  for (int l=0; l<numNodes(); ++l) {
    Vec3D u = varFcn->getVelocity(V[ nodeNum(l) ]);
    Vec3D du = varFcn->getVelocity(dV[ nodeNum(l) ]);
    double a = varFcn->computeSoundSpeed(V[ nodeNum(l) ]);
    double da = varFcn->computeDerivativeOfSoundSpeed(V[ nodeNum(l) ], dV[ nodeNum(l) ], dMach);
    double un = u * n - ndot;
    double dun = du * n + u * dn - dndot;
    double locMach = varFcn->computeMachNumber(V[ nodeNum(l) ]);
    //double locMach = fabs(un/a); //local Preconditioning (ARL)
    double locbeta = tprec.getBeta(locMach,true);
    double dLocMach = varFcn->computeDerivativeOfMachNumber(V[ nodeNum(l) ], dV[ nodeNum(l) ], dMach);
    double dbeta = tprec.getdBeta(locMach,dLocMach,true);

    double beta2 = locbeta * locbeta;
    double dbeta2 = 2.0 * locbeta * dbeta;
    double coeff1 = (1.0+beta2)*un;
    double dCoeff1 = dbeta2*un + (1.0+beta2)*dun;
    double coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);
    double dCoeff2 = (((1.0-beta2)*un)*((-dbeta2*un) + ((1.0-beta2)*dun)) + (2.0*locbeta*a)*(2.0*dbeta*a+2.0*locbeta*da)) / pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);

    if (min(0.5*(coeff1-coeff2), 0.0) != 0.0)
      dIdti[ nodeNum(l) ] += 0.5*(dCoeff1-dCoeff2)* S/numNodes() + 0.5*(coeff1-coeff2)* dS/numNodes();
      
    double vis = 0.0;
    double dvis = 0.0;
    if (fet) {
      vis = fet->computeViscousTimeStep(X[nodeNum(l)],V[nodeNum(l)]);
      dvis = fet->computeDerivativeOfViscousTimeStep(X[nodeNum(l)],dX[nodeNum(l)],V[nodeNum(l)],dV[nodeNum(l)],dMach);
    }
    dIdtv[ nodeNum(l) ] += dvis*S*S + vis*2.0*S*dS;

  }
    
}

//------------------------------------------------------------------------------

template<int dim>
inline
void Face::computeFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normals, 
				   Vec<double> &normalVel, SVec<double,dim> &V,
				   double *Ub, SVec<double,dim> &fluxes, double UbHH)
{

  // UH (07/2012)
  // The next test is to handle the cases where code is set
  // to BC_KIRCHHOFF_SURFACE (== 9)
  if ((code < BC_MIN_CODE) | (code > BC_MAX_CODE))
    return;

  if(fluxFcn[code]){
    bool farfield = (code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED || code == BC_INLET_MOVING || code == BC_INLET_FIXED);
    const int dim0 = 7;
    double flux[2*dim0];
    double hh[6];

    fluxFcn[code]->setHHCoeffPointer(hh);
    for(int j=0; j<3; j++)
      hh[j] = faceCenter[j];
    hh[4] = UbHH;

    for (int l=0; l<numNodes(); ++l) {

      fluxFcn[code]->compute(0.0, 0.0, getNormal(normals, l), getNormalVel(normalVel, l),
                             V[nodeNum(l)], Ub, flux);

      for (int k=0; k<dim; ++k){ 
        fluxes[ nodeNum(l) ][k] += flux[k];
      }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void Face::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                   FluxFcn **fluxFcn, Vec<Vec3D> &normals,
                                   Vec<double> &normalVel, SVec<double,dim> &V,
                                   double *Ub, SVec<double,dim> &fluxes, double UbHH)
{

  // UH (07/2012)
  // The next test is to handle the cases where code is set
  // to BC_KIRCHHOFF_SURFACE (== 9)
  if ((code < BC_MIN_CODE) | (code > BC_MAX_CODE))
    return;

  if(code == BC_ADIABATIC_WALL_MOVING  || code == BC_ADIABATIC_WALL_FIXED ||
     code == BC_SLIP_WALL_MOVING       || code == BC_SLIP_WALL_FIXED      ||
     code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED ||
     code == BC_POROUS_WALL_MOVING     || code == BC_POROUS_WALL_FIXED) {
  // FS Riemann based flux calculation.
    double flux[dim], Wstar[2*dim], Vi[2*dim];
    int k;
    Vec3D wallVel, unitNormal;
    VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();

    for (int l=0; l<numNodes(); ++l) {
      k = nodeNum(l);
      unitNormal = getNormal(normals,l)/(getNormal(normals,l).norm());
      wallVel = getNormalVel(normalVel, l)/(getNormal(normals,l).norm())*unitNormal;
      for(int iDim=0; iDim<dim; iDim++)
        Vi[iDim] = Vi[iDim+dim] = V[k][iDim];

      riemann.computeFSIRiemannSolution(Vi, wallVel, -1.0*unitNormal, varFcn, Wstar, k/*dummy variable here*/);

//      fluxFcn[BC_INTERNAL]->compute(0.0, 0.0, getNormal(normals,l), getNormalVel(normalVel,l), V[k], Wstar, flux);
      fluxFcn[code]->compute(0.0, 0.0, getNormal(normals, l), getNormalVel(normalVel, l), Wstar, Ub, flux);

      for (int i=0; i<dim; ++i)
        fluxes[k][i] += flux[i];
    }
    return;
  }

  if(fluxFcn[code]) {
    bool farfield = (code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED || code == BC_INLET_MOVING || code == BC_INLET_FIXED);
    const int dim0 = 7;
//    double* flux;
    double flux[2*dim0];
    double hh[6];
    fluxFcn[code]->setHHCoeffPointer(hh);

    for(int j=0; j<3; j++)
      hh[j] = faceCenter[j];

    hh[4] = UbHH;

    for (int l=0; l<numNodes(); ++l) {

      fluxFcn[code]->compute(0.0, 0.0, getNormal(normals, l), getNormalVel(normalVel, l),
                           V[nodeNum(l)], Ub, flux);

      for (int k=0; k<dim; ++k) {
	fluxes[ nodeNum(l) ][k] += flux[k];
      }

    }
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
inline
void Face::computeDerivativeOfFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normals,
				      Vec<Vec3D> &dNormals, Vec<double> normalVel, Vec<double> dNormalVel,
				      SVec<double,dim> &V, double *Ub,
				      double *dUb, SVec<double,dim> &dFluxes)
{

  // UH (07/2012)
  // The next test is to handle the cases where code is set
  // to BC_KIRCHHOFF_SURFACE (== 9)
  if ((code < BC_MIN_CODE) | (code > BC_MAX_CODE))
    return;
    
  if(fluxFcn[code]){
    double flux[dim];
    double dFlux[dim];
    for (int l=0; l<numNodes(); ++l) {
//      double dFluxdNormal[7][3] = {0}, dFluxdNormalVel[7][1] = {0}, dFluxdub[7][7] = {0}, dFLUX[7] = {0}, diff[dim] = {0};
//      for(int i=0; i<dim; i++) dFLUX[i] = dFlux[i] = 0.0;
      fluxFcn[code]->computeDerivative(0.0, 0.0, getNormal(normals, l), getdNormal(dNormals, l), getNormalVel(normalVel, l), getdNormalVel(dNormalVel, l), V[nodeNum(l)], Ub, dUb, flux, dFlux);
/*
      fluxFcn[code]->computeDerivativeOperators(0.0, 0.0, getNormal(normals, l), getNormalVel(normalVel, l), V[nodeNum(l)], Ub, dFluxdNormal, dFluxdNormalVel, dFluxdub);
      for(int i=0; i<dim; ++i) {
        dFLUX[i] += 1.0/3.0*dFluxdNormalVel[i][0]*dNormalVel[normNum];
        for(int j=0; j<dim; ++j)
          dFLUX[i] += dFluxdub[i][j]*dUb[j];
        for(int j=0; j<3; ++j)
          dFLUX[i] += 1.0/3.0*dFluxdNormal[i][j]*dNormals[normNum][j];
      }
      double diffnorm(0), dFluxnorm(0), dFLUXnorm(0);
      for(int i=0; i<dim; ++i) {
        diff[i] = dFlux[i] - dFLUX[i];
        diffnorm += diff[i]*diff[i];
        dFluxnorm += dFlux[i]*dFlux[i];
        dFLUXnorm += dFLUX[i]*dFLUX[i];
      }
      diffnorm = sqrt(diffnorm);
      dFluxnorm = sqrt(dFluxnorm);
      dFLUXnorm = sqrt(dFLUXnorm);
      if(dFluxnorm != 0) fprintf(stderr, " .*.*. rel. diff = %e, dFluxnorm = %e, dFLUXnorm = %e\n", diffnorm/dFluxnorm, dFluxnorm, dFLUXnorm);
      else fprintf(stderr, " .*.*. abs. diff = %e, dFluxnorm = %e, dFLUXnorm = %e\n", diffnorm, dFluxnorm, dFLUXnorm);
*/
      for (int k=0; k<dim; ++k){
        dFluxes[ nodeNum(l) ][k] += dFlux[k];
      }
    }
  }

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
inline
void Face::computeDerivativeOperatorsOfFiniteVolumeTerm(int faceNum, FluxFcn **fluxFcn, Vec<Vec3D> &normals,
				      Vec<double> normalVel, SVec<double,dim> &V, double *Ub,
              RectangularSparseMat<double,3,dim> &dFluxdFaceNormal,
              RectangularSparseMat<double,1,dim> &dFluxdFaceNormalVel,
              RectangularSparseMat<double,dim,dim> &dFluxdUb)
{

  // UH (07/2012)
  // The next test is to handle the cases where code is set
  // to BC_KIRCHHOFF_SURFACE (== 9)
  if ((code < BC_MIN_CODE) | (code > BC_MAX_CODE))
    return;
    
  if(fluxFcn[code]){
    double flux[dim];
    double dFlux[dim];
    for (int l=0; l<numNodes(); ++l) {
      double dFluxdNormal[7][3] = {0}, dFluxdNormalVel[7][1] = {0}, dFluxdub[7][7] = {0};
      fluxFcn[code]->computeDerivativeOperators(0.0, 0.0, getNormal(normals, l), getNormalVel(normalVel, l), V[nodeNum(l)], Ub, dFluxdNormal, dFluxdNormalVel, dFluxdub);
      double dFluxdNormal2[dim][3] = {0}, dFluxdNormalVel2[dim][1] = {0}, dFluxdub2[dim][dim] = {0};
      for(int i=0; i<dim; ++i) {
        dFluxdNormalVel2[i][0] = dFluxdNormalVel[i][0];
        for(int j=0; j<dim; ++j)
          dFluxdub2[i][j] = dFluxdub[i][j];
        for(int j=0; j<3; ++j)
          dFluxdNormal2[i][j] = dFluxdNormal[i][j];
      }

      dFluxdFaceNormal.addContrib(nodeNum(l), faceNum, dFluxdNormal2[0]);
      dFluxdFaceNormalVel.addContrib(nodeNum(l), faceNum, dFluxdNormalVel2[0]);
      dFluxdUb.addContrib(nodeNum(l), faceNum, dFluxdub2[0]);
    }
  }

}

//------------------------------------------------------------------------------
 //d2d embedded
template<int dim>
inline
void Face::computeFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normals,
				   Vec<double> &normalVel, SVec<double,dim> &V,
				   double *Ub, Vec<int> &fluidId, 
				   SVec<double,dim> &fluxes, LevelSetStructure *LSS,
                                   double UbHH)
{

  // UH (07/2012)
  // The next test is to handle the cases where code is set
  // to BC_KIRCHHOFF_SURFACE (== 9)
  if ((code < BC_MIN_CODE) | (code > BC_MAX_CODE))
    return;
  
  Vec3D normal = getNormal(normals);
  const int dim0 = 7;
  //double* flux;
  double flux[2*dim0];
  bool cracking = LSS ? LSS->withCracking() : false;
  bool farfield = (code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED || code == BC_INLET_MOVING || code == BC_INLET_FIXED);
    
  double hh[6];
  // else
  //  flux = new double [dim];

  if(fluxFcn[code]){
    fluxFcn[code]->setHHCoeffPointer(hh);
    for(int j=0; j<3; j++)
      hh[j] = faceCenter[j];
    hh[4] = UbHH;

    for (int l=0; l<numNodes(); ++l) {
      if(cracking) {
        if(LSS->isOccluded(0.0,nodeNum(l))) continue;}
      else {
        if(LSS && !LSS->isActive(0.0, nodeNum(l))) continue;}

	fluxFcn[code]->compute(0.0, 0.0, getNormal(normals, l), getNormalVel(normalVel, l), 
			       V[nodeNum(l)], Ub, flux, fluidId[nodeNum(l)]);

	for (int k=0; k<dim; ++k)
	  fluxes[ nodeNum(l) ][k] += flux[k];
    }
  }

//  delete [] flux; 
}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
inline
void Face::computeFiniteVolumeTermLS(FluxFcn **fluxFcn, Vec<Vec3D> &normals,
				     Vec<double> &normalVel, SVec<double,dim> &V,
				     SVec<double,dimLS> &Phi, SVec<double,dimLS> &PhiF)
{

  // UH (07/2012)
  // The next test is to handle the cases where code is set
  // to BC_KIRCHHOFF_SURFACE (== 9)
  if ((code < BC_MIN_CODE) | (code > BC_MAX_CODE))
    return;
  
  double Uf = 0.0;
  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED ||
      code == BC_ADIABATIC_WALL_MOVING  || code == BC_ADIABATIC_WALL_FIXED  ||
      code == BC_SLIP_WALL_MOVING       || code == BC_SLIP_WALL_FIXED       ||
      code == BC_POROUS_WALL_MOVING     || code == BC_POROUS_WALL_FIXED     ||
      code == BC_SYMMETRY
      ) {
    //at wall either U.n = Uwall.n (Euler) or U = Uwall (Navier-Stokes)
    //and thus the flux is 0.0
  }else{
    for (int l=0; l<numNodes(); ++l) {
      Vec3D normal = getNormal(normals, l);
      Uf   = ( V[nodeNum(l)][1]*normal[0] +
               V[nodeNum(l)][2]*normal[1] +
               V[nodeNum(l)][3]*normal[2] ) -
             getNormalVel(normalVel, l);
      for (int k=0; k<dimLS; k++)
        PhiF[ nodeNum(l) ][k] += Uf*(-1.0);//Phi[nodeNum(l)][k];
    }

  }

}

template<int dim>
inline
void Face::computeHHBoundaryTermResidual(SVec<double,dim> &U,double* Ub,double& UbHH, double& res, VarFcn* vf) {

  if ((code < BC_MIN_CODE) | (code > BC_MAX_CODE))
    return;
  
  bool farfield = (code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED || code == BC_INLET_MOVING || code == BC_INLET_FIXED);
  if (!farfield)
    return;

  double radius = sqrt(faceCenter[0]*faceCenter[0]+
                       faceCenter[1]*faceCenter[1]+
                       faceCenter[2]*faceCenter[2]);

  double r[3] = {faceCenter[0]/radius, faceCenter[1]/radius, faceCenter[2]/radius};

  double hh[6];
  double locgam = 0;
  switch (vf->getType()) {

    case VarFcnBase::PERFECTGAS:
    case VarFcnBase::STIFFENEDGAS:
      locgam = vf->getGamma();
      break;
    case VarFcnBase::TAIT:
      locgam = vf->getBetaWater();
      break;
    default:
      std::cout << "Unsupported farfield fluid for Modified Ghidaglia!" << std::endl;
      exit(-1);
      break;
  }

  double Vb[dim];
  vf->conservativeToPrimitive(Ub,Vb);

  double c = 0.0;
  double V[dim];
  //double Rplus = 0.0;
  for (int l = 0; l < numNodes(); ++l) {
  
    vf->conservativeToPrimitive(U[nodeNum(l)],V);
    c += vf->computeSoundSpeed(V);
    //Rplus += V[1]*r[0]+V[2]*r[1]+V[3]*r[2];
  }
  c /= (numNodes());
  //Rplus /= (numNodes());

  double cinf = vf->computeSoundSpeed(Vb);
  
  //double Rinf = Vb[1]*r[0]+Vb[2]*r[1]+Vb[3]*r[2];
  
  double dSdt3 = 2.0*(c-cinf)*cinf/(radius*(locgam-1.0));
  if (aerof_isnan(dSdt3))
    std::cout << dSdt3 << " " << c << " " << cinf << " " << radius << " " << locgam << std::endl;
  
  res = -dSdt3;
}

template<class Scalar,int dim,int neq>
inline
void Face::computeHHBoundaryTermJacobian(int faceid,FluxFcn **fluxFcn, SVec<double,dim> &U,
				       double* Ub, GenMat<Scalar,neq> &A, VarFcn* vf,
                                       double& UbHH,
				       Vec<Vec3D> &normals, Vec<double> &normalVel) {

  if ((code < BC_MIN_CODE) | (code > BC_MAX_CODE))
    return;
  
  bool farfield = (code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED || code == BC_INLET_MOVING || code == BC_INLET_FIXED);
  if (!farfield)
    return;

  SVec<Scalar,neq*3>& jac = *A.getHU();

  double radius = sqrt(faceCenter[0]*faceCenter[0]+
                       faceCenter[1]*faceCenter[1]+
                       faceCenter[2]*faceCenter[2]);

  double r[3] = {faceCenter[0]/radius, faceCenter[1]/radius, faceCenter[2]/radius};
  double hh[6];
  fluxFcn[code]->setHHCoeffPointer(hh);
  for(int j=0; j<3; j++)
    hh[j] = faceCenter[j];
  hh[4] = UbHH;


  double locgam = 0;
  switch (vf->getType()) {

    case VarFcnBase::PERFECTGAS:
    case VarFcnBase::STIFFENEDGAS:
      locgam = vf->getGamma();
      break;
    case VarFcnBase::TAIT:
      locgam = vf->getBetaWater();
      break;
    default:
      std::cout << "Unsupported farfield fluid for Modified Ghidaglia!" << std::endl;
      exit(-1);
      break;
  }

  double Vb[dim];
  vf->conservativeToPrimitive(Ub,Vb);

  double c = 0.0;
  double V[dim];
  //double Rplus = 0.0;
  int nn = numNodes();
  double locjac[dim],jacu[dim];
  double jacuv[dim*dim];
  memset(locjac,0,sizeof(double)*dim);
  double cinf = vf->computeSoundSpeed(Vb);
  double gg = 2.0*cinf/(radius*(locgam-1.0));
  for (int l = 0; l < numNodes(); ++l) {
  
    vf->conservativeToPrimitive(U[nodeNum(l)],V);
    c = vf->computeSoundSpeed(V);
    switch (vf->getType()) {

    case VarFcnBase::PERFECTGAS:
    case VarFcnBase::STIFFENEDGAS:

      locjac[0] = -0.5*c / V[0];
      locjac[4] = 0.5/c*locgam/V[0];
      break;
    case VarFcnBase::TAIT:
      locjac[0] = 0.5*c*(locgam-1.0) / V[0];
      
      break;
    default:
      std::cout << "Unsupported farfield fluid for Modified Ghidaglia!" << std::endl;
      exit(-1);
      break;
    }

    vf->multiplyBydVdUT(V, locjac, jacu);
    for (int k = 0; k < neq; ++k) 
      jac[faceid][l*neq+k] += -gg*jacu[k]/nn;

    Vec3D  normal = getNormal(normals, l);
    double normVel= getNormalVel(normalVel, l);


    fluxFcn[code]->computeJacobianFarfield(1.0, 0.0, normal, normVel, V, Ub, jacu);
    Scalar *Auh = A.getElemUH(faceid);
    for (int k=0; k<neq; ++k) 
      Auh[neq*l+k] += jacu[k];

    //Rplus += V[1]*r[0]+V[2]*r[1]+V[3]*r[2];
  }
  //Rplus /= (numNodes());

  
  //double Rinf = Vb[1]*r[0]+Vb[2]*r[1]+Vb[3]*r[2];
  
  /*double dSdt3 = 2.0*(c-cinf)*cinf/(radius*(locgam-1.0));
  if (aerof_isnan(dSdt3))
    std::cout << dSdt3 << " " << c << " " << cinf << " " << radius << " " << locgam << std::endl;
  
  res = -dSdt3;
  */
}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
inline
void Face::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normals, 
					   Vec<double> &normalVel, SVec<double,dim> &V, 
					   double *Ub, GenMat<Scalar,neq> &A,double UbHH)
{

  // UH (07/2012)
  // The next test is to handle the cases where code is set
  // to BC_KIRCHHOFF_SURFACE (== 9)
  if ((code < BC_MIN_CODE) | (code > BC_MAX_CODE))
    return;
  
  double jac[neq*neq];
  bool farfield = (code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED || code == BC_INLET_MOVING || code == BC_INLET_FIXED);
  double hh[6];
  fluxFcn[code]->setHHCoeffPointer(hh);

  for(int j=0; j<3; j++)
    hh[j] = faceCenter[j];
  hh[4] = UbHH;
 
  for (int l=0; l<numNodes(); ++l) {
    Vec3D  normal = getNormal(normals, l);
    double normVel= getNormalVel(normalVel, l);


    fluxFcn[code]->computeJacobian(1.0, 0.0, normal, normVel, V[nodeNum(l)], Ub, jac);
    Scalar *Aii = A.getElem_ii(nodeNum(l));
    for (int k=0; k<neq*neq; ++k) 
      Aii[k] += jac[k];
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
inline
void Face::computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann, FluxFcn **fluxFcn, Vec<Vec3D> &normals,
                                           Vec<double> &normalVel, SVec<double,dim> &V,
                                           double *Ub, GenMat<Scalar,neq> &A,double UbHH)
{

  // UH (07/2012)
  // The next test is to handle the cases where code is set
  // to BC_KIRCHHOFF_SURFACE (== 9)
  if ((code < BC_MIN_CODE) | (code > BC_MAX_CODE))
    return;
  
  if(code == BC_ADIABATIC_WALL_MOVING  || code == BC_ADIABATIC_WALL_FIXED ||
     code == BC_SLIP_WALL_MOVING       || code == BC_SLIP_WALL_FIXED      ||
     code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED ||
     code == BC_POROUS_WALL_MOVING     || code == BC_POROUS_WALL_FIXED) {
  // FS Riemann based flux calculation.
    double flux[dim], Wstar[2*dim], Vi[2*dim];
    double dUdU[neq*neq],dfdUi[neq*neq],dkk[neq*neq],dWdW[dim*dim],dWdU[dim*dim];

    int k;
    Vec3D wallVel, unitNormal;
    VarFcn *varFcn = fluxFcn[BC_INTERNAL]->getVarFcn();

    for (int l=0; l<numNodes(); ++l) {
      Scalar *Aii = A.getElem_ii(nodeNum(l));
      k = nodeNum(l);
      unitNormal = getNormal(normals,l)/(getNormal(normals,l).norm());
      wallVel = getNormalVel(normalVel, l)/(getNormal(normals,l).norm())*unitNormal;
      for(int iDim=0; iDim<dim; iDim++)
        Vi[iDim] = Vi[iDim+dim] = V[k][iDim];

      riemann.computeFSIRiemannSolution(Vi,wallVel,-1.0*unitNormal,varFcn,Wstar,k);
      riemann.computeFSIRiemannJacobian(Vi,wallVel,-1.0*unitNormal,varFcn,Wstar,k,dWdW);

      varFcn->postMultiplyBydVdU(Wstar, dWdW, dWdU);
      varFcn->preMultiplyBydUdV(Vi, dWdU, dUdU);

      fluxFcn[code]->computeJacobian(1.0, 0.0, getNormal(normals,l), getNormalVel(normalVel,l), Wstar, Ub, dfdUi);
      DenseMatrixOp<double, neq, neq*neq>::applyToDenseMatrix(&dfdUi,0,&dUdU, 0, &dkk,0);
      Aii = A.getElem_ii(k);
      for (int s=0; s<neq*neq; ++s) {
        Aii[s] += dkk[s];
      }
    }
    return;
  }

  bool farfield = (code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED || code == BC_INLET_MOVING || code == BC_INLET_FIXED);
  double hh[6];
  fluxFcn[code]->setHHCoeffPointer(hh);

  for(int j=0; j<3; j++)
    hh[j] = faceCenter[j];
  hh[4] = UbHH;
 
  double jac[neq*neq];
  for (int l=0; l<numNodes(); ++l) {
    Vec3D  normal = getNormal(normals, l);
    double normVel= getNormalVel(normalVel, l);

    fluxFcn[code]->computeJacobian(1.0, 0.0, normal, normVel, V[nodeNum(l)], Ub, jac);
    Scalar *Aii = A.getElem_ii(nodeNum(l));
    for (int k=0; k<neq*neq; ++k)
      Aii[k] += jac[k];
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
inline
void Face::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normals,
					   Vec<double> &normalVel, SVec<double,dim> &V,
					   double *Ub, GenMat<Scalar,neq> &A, int* nodeType,double UbHH)
{

  // UH (07/2012)
  // The next test is to handle the cases where code is set
  // to BC_KIRCHHOFF_SURFACE (== 9)
  if ((code < BC_MIN_CODE) | (code > BC_MAX_CODE))
    return;
  
  bool farfield = (code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED || code == BC_INLET_MOVING || code == BC_INLET_FIXED);
    
  double hh[6];
  for(int j=0; j<3; j++)
    hh[j] = faceCenter[j];
  hh[4] = UbHH;

  double jac[neq*neq];
  for (int l=0; l<numNodes(); ++l) {
    if(!(code == BC_INLET_MOVING || code == BC_OUTLET_MOVING ||
         code == BC_INLET_FIXED  || code == BC_OUTLET_FIXED)) {
      Vec3D normal = getNormal(normals, l);
      double normVel= getNormalVel(normalVel, l);

      fluxFcn[code]->computeJacobian(1.0, 0.0, normal, normVel, V[nodeNum(l)], Ub, jac);
      Scalar *Aii = A.getElem_ii(nodeNum(l));
      for (int k=0; k<neq*neq; ++k)
        Aii[k] += jac[k];

    }
  }
}

//------------------------------------------------------------------------------

//d2d embedded
template<int dim, class Scalar, int neq>
inline
void Face::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normals,
					   Vec<double> &normalVel, SVec<double,dim> &V, 
					   double *Ub, GenMat<Scalar,neq> &A, Vec<int> &fluidId,
                                           LevelSetStructure* LSS,double UbHH)
{

  // UH (07/2012)
  // The next test is to handle the cases where code is set
  // to BC_KIRCHHOFF_SURFACE (== 9)
  if ((code < BC_MIN_CODE) | (code > BC_MAX_CODE))
    return;
  
  bool farfield = (code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED || code == BC_INLET_MOVING || code == BC_INLET_FIXED);
  double hh[6];
  fluxFcn[code]->setHHCoeffPointer(hh);
 
  for(int j=0; j<3; j++)
    hh[j] = faceCenter[j];
  hh[4] = UbHH;
 
  double jac[neq*neq];
  for (int l=0; l<numNodes(); ++l) {
    Vec3D  normal = getNormal(normals, l);
    double normVel= getNormalVel(normalVel, l);

    if(LSS && !LSS->isActive(0.0, nodeNum(l))) continue;

    fluxFcn[code]->computeJacobian(1.0, 0.0, normal, normVel, V[nodeNum(l)], Ub, jac, fluidId[nodeNum(l)]);
    Scalar *Aii = A.getElem_ii(nodeNum(l));

    for (int k=0; k<neq*neq; ++k) 
      Aii[k] += jac[k];
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
inline
void Face::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normals,
                                           Vec<double> &normalVel, SVec<double,dim> &V,
                                           double *Ub, GenMat<Scalar,neq> &A, 
                                           Vec<int> &fluidId, int* nodeType,double UbHH)
{

  // UH (07/2012)
  // The next test is to handle the cases where code is set
  // to BC_KIRCHHOFF_SURFACE (== 9)
  if ((code < BC_MIN_CODE) | (code > BC_MAX_CODE))
    return;
  
  bool farfield = (code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED || code == BC_INLET_MOVING || code == BC_INLET_FIXED);
  double hh[6];
  fluxFcn[code]->setHHCoeffPointer(hh);
  for(int j=0; j<3; j++)
    hh[j] = faceCenter[j];
  hh[4] = UbHH;

  double jac[neq*neq];
  for (int l=0; l<numNodes(); ++l) {
    if(!(code == BC_INLET_MOVING || code == BC_OUTLET_MOVING ||
         code == BC_INLET_FIXED  || code == BC_OUTLET_FIXED)) {
      Vec3D normal = getNormal(normals, l);
      double normVel= getNormalVel(normalVel, l);

      fluxFcn[code]->computeJacobian(1.0, 0.0, normal, normVel, V[nodeNum(l)], Ub, jac, fluidId[nodeNum(l)]);
      Scalar *Aii = A.getElem_ii(nodeNum(l));
      for (int k=0; k<neq*neq; ++k)
        Aii[k] += jac[k];

    }
  }

}

//------------------------------------------------------------------------------
template<int dim, class Scalar, int dimLS>
inline
void Face::computeJacobianFiniteVolumeTermLS(Vec<Vec3D> &normals,
					     Vec<double> &normalVel, SVec<double,dim> &V,
					     GenMat<Scalar,dimLS> &A)
{

  // UH (07/2012)
  // The next test is to handle the cases where code is set
  // to BC_KIRCHHOFF_SURFACE (== 9)
  if ((code < BC_MIN_CODE) | (code > BC_MAX_CODE))
    return;
  
  double Uf = 0.0;
  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED ||
      code == BC_ADIABATIC_WALL_MOVING  || code == BC_ADIABATIC_WALL_FIXED  ||
      code == BC_SLIP_WALL_MOVING       || code == BC_SLIP_WALL_FIXED       ||
      code == BC_POROUS_WALL_MOVING     || code == BC_POROUS_WALL_FIXED     ||
      code == BC_SYMMETRY
      ) {
    //at wall either U.n = Uwall.n (Euler) or U = Uwall (Navier-Stokes)
    //and thus the flux is 0.0
  }else{
    for (int l=0; l<numNodes(); ++l) {
      Scalar *Aii = A.getElem_ii(nodeNum(l));
      Vec3D normal = getNormal(normals, l);
      Uf   = ( V[nodeNum(l)][1]*normal[0] +
               V[nodeNum(l)][2]*normal[1] +
               V[nodeNum(l)][3]*normal[2] ) -
             getNormalVel(normalVel, l);
      *Aii /*+= PhiF[ nodeNum(l) ]*/ += 0.0;//Uf;
    }

  }
  /*  double jac;
  for (int l=0; l<numNodes(); ++l) {
    Vec3D normal = getNormal(normals, l);
    double normVel= getNormalVel(normalVel, l);

    jac  = V[nodeNum(l)][1]*normal[0]  + V[nodeNum(l)][2]*normal[1]  + V[nodeNum(l)][3]*normal[2];
    jac *= V[nodeNum(l)][0]; // why did sriram write this? not true!!!
    exit(1);
    Scalar *Aii = A.getElem_ii(nodeNum(l));
    for (int k=0; k<dimLS*dimLS; ++k)
      Aii[k] += jac;
      }*/
}

template <class T>
inline double mult_local_face_hh(T t, double a) {
  return t*a;
}

template <>
inline double mult_local_face_hh<bcomp>(bcomp t, double a) {
  fprintf(stderr,"Error, Face.C %d, incompatible types\n",__LINE__);
  return 0.0;
}

template <class Scalar,int dim>
void Face::computeMatVecProdH1FarFieldHH(int faceid,GenMat<Scalar,dim> &A, SVec<double,dim> &p_u,
	                                 SVec<double,dim> &prod_u,double& p_hh, double& prod_hh) {

  bool farfield = (code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED || code == BC_INLET_MOVING || code == BC_INLET_FIXED);
  if (!farfield)
    return;

  Scalar *Auh = A.getElemUH(faceid);
  Scalar *Ahu = A.getElemHU(faceid);
  Scalar *Ahh = A.getElemHH(faceid);
  
  for (int l=0; l<numNodes(); ++l) {

    double* ppu = p_u[nodeNum(l)];
    double* produ = prod_u[nodeNum(l)];
    for (int k = 0; k < dim; ++k) {
   
      produ[k] += mult_local_face_hh(Auh[l*dim+k],p_hh);
    }
  
    for (int k = 0; k < dim; ++k) {
      
      prod_hh += mult_local_face_hh(Ahu[l*dim+k],ppu[k]);
    }
  }

  prod_hh +=  mult_local_face_hh((*Ahh),p_hh);
}

//------------------------------------------------------------------------------

template<int dim>
void FaceSet::computeTimeStep(VarFcn *varFcn, GeoState &geoState, 
			      SVec<double,dim> &V, Vec<double> &dt,
			      TimeLowMachPrec &tprec,
                              Vec<int> &fluidId)
{

  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();

  for (int i=0; i<numFaces; ++i)
    faces[i]->computeTimeStep(varFcn, n, ndot, V, dt, tprec, fluidId);

}

//------------------------------------------------------------------------------

template<int dim>
void FaceSet::computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn, GeoState &geoState, 
			      SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &idti,
										Vec<double> &idtv, TimeLowMachPrec &tprec, LevelSetStructure *LSS)
{
  
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();

  for (int i=0; i<numFaces; ++i)
	  faces[i]->computeTimeStep(fet, varFcn, n, ndot, X, V, idti, idtv, tprec, LSS);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void FaceSet::computeDerivativeOfTimeStep(FemEquationTerm *fet, VarFcn *varFcn, GeoState &geoState, 
			      SVec<double,3> &X, SVec<double,3> &dX, SVec<double,dim> &V, SVec<double,dim> &dV, 
			      Vec<double> &dIdti,Vec<double> &dIdtv, double dMach, 
                              TimeLowMachPrec &tprec)
{

  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<Vec3D> &dn = geoState.getdFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  Vec<double> &dndot = geoState.getdFaceNormalVel();

  for (int i=0; i<numFaces; ++i)
    faces[i]->computeDerivativeOfTimeStep(fet, varFcn, n, dn, ndot, dndot, X, dX, V, dV, dIdti, dIdtv, dMach, tprec);

}

//------------------------------------------------------------------------------

template<int dim>
void FaceSet::computeFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
				      GeoState &geoState, SVec<double,dim> &V, 
				      SVec<double,dim> &fluxes)
{
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();
  Vec<double>* UbHH = bcData.getBoundaryStateHH();

	if (sampleMesh) {
		int i;
		for (int iFace=0; iFace<numSampledFaces; ++iFace) {
			i = facesConnectedToSampleNode[iFace];
                        if (!UbHH)
  			  faces[i]->computeFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], fluxes);
                        else
                          faces[i]->computeFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], fluxes,
                                                            (*UbHH)[i]);
		}
	}
	else {
		for (int i=0; i<numFaces; ++i) {
                  if (!UbHH)
  	            faces[i]->computeFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], fluxes);
                  else
                    faces[i]->computeFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], fluxes,
                                                      (*UbHH)[i]);
 	  }
	}
}


//------------------------------------------------------------------------------

template<int dim>
void FaceSet::computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann, 
                                      FluxFcn **fluxFcn, BcData<dim> &bcData,
                                      GeoState &geoState, SVec<double,dim> &V,
                                      SVec<double,dim> &fluxes)
{
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();

  Vec<double>* UbHH = bcData.getBoundaryStateHH();

  for (int i=0; i<numFaces; ++i) {
    if (!UbHH)
      faces[i]->computeFiniteVolumeTerm(riemann, fluxFcn, n, ndot, V, Ub[i], fluxes);
    else
      faces[i]->computeFiniteVolumeTerm(riemann, fluxFcn, n, ndot, V, Ub[i], fluxes,
                                        (*UbHH)[i]);
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void FaceSet::computeDerivativeOfFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
				      GeoState &geoState, SVec<double,dim> &V,
				      SVec<double,dim> &dFluxes)
{

  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<Vec3D> &dn = geoState.getdFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  Vec<double> &dndot = geoState.getdFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();
  SVec<double,dim> &dUb = bcData.getdFaceStateVector();
/*
//  SVec<double,dim> dFlux2(dFluxes), dummy(dFluxes), diff(dFluxes);
  SVec<double,dim> dummy(dFluxes);

  if(isSparse) {
    dn *= 1.0/3.0;
    dndot *= 1.0/3.0;
    dummy = 0.0;
    dFluxdFaceNormal->apply(dn, dummy); 
    dFluxes += dummy;
    dFluxdFaceNormalVel->apply(dndot, dummy);
    dFluxes += dummy;
    dummy = 0.0;
    dFluxdUb->apply(dUb, dummy);
    dFluxes += dummy; 
//    diff = 0.0;
    dn *= 3.0;
    dndot *= 3.0;
    return;
  }
*/
  for (int i=0; i<numFaces; ++i)
    faces[i]->computeDerivativeOfFiniteVolumeTerm(fluxFcn, n, dn, ndot, dndot, V, Ub[i], dUb[i], dFluxes);
/*
  diff = dFluxes - dFlux2;
  double diffnorm(0), dFluxesnorm(0), dFlux2norm(0);
  diffnorm = diff.norm();
  dFluxesnorm = dFluxes.norm();
  dFlux2norm = dFlux2.norm();
  if(dFluxesnorm != 0) fprintf(stderr, " ... rel. diff = %e, dFlux2norm = %e, dFluxesnorm = %e\n", diffnorm/dFluxesnorm, dFlux2norm, dFluxesnorm);
  else fprintf(stderr, " ... abs. diff = %e, dFlux2norm = %e, dFluxesnorm = %e\n", diffnorm, dFlux2norm, dFluxesnorm); 
*/
}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void FaceSet::computeDerivativeOfFiniteVolumeTerm(
              RectangularSparseMat<double,3,dim> *dFluxdFaceNormal,
              RectangularSparseMat<double,1,dim> *dFluxdFaceNormalVel,
              RectangularSparseMat<double,dim,dim> *dFluxdUb,
              BcData<dim> &bcData,
				      GeoState &geoState, 
              Vec<Vec3D>& dFaceNormal,
              Vec<double>& dFaceNormalVel,
				      SVec<double,dim> &dFluxes)
{

  SVec<double,dim> dummy(dFluxes);

  dFaceNormal *= 1.0/3.0;
  dFaceNormalVel *= 1.0/3.0;
  dummy = 0.0;
  dFluxdFaceNormal->apply(dFaceNormal, dummy); 
  dFluxes += dummy;
  dFluxdFaceNormalVel->apply(dFaceNormalVel, dummy);
  dFluxes += dummy;
/*  dummy = 0.0;
  dFluxdUb->apply(dUb, dummy);  //TODO: assumed dUb is zero
  dFluxes += dummy; 
*/
  dFaceNormal *= 3.0;
  dFaceNormalVel *= 3.0;

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void FaceSet::computeTransposeDerivativeOfFiniteVolumeTerm(RectangularSparseMat<double,3,dim> *dFluxdFaceNormal,
                                                           RectangularSparseMat<double,1,dim> *dFluxdFaceNormalVel,
                                                           RectangularSparseMat<double,dim,dim> *dFluxdUb,
                                                           BcData<dim> &bcData,
                                                           GeoState &geoState,
                                                           SVec<double,dim> &dFluxes,
                                                           Vec<Vec3D>& dFaceNormal2,
                                                           Vec<double>& dFaceNormalVel2)
{

//  Vec<Vec3D> &dn = geoState.getdFaceNormal();
//  Vec<double> &dndot = geoState.getdFaceNormalVel();
//  SVec<double,dim> &dUb = bcData.getdFaceStateVector();

  Vec<Vec3D> dFaceNormal2dummy(dFaceNormal2);
  Vec<double> dFaceNormalVel2dummy(dFaceNormalVel2); 

  dFluxes *= 1.0/3.0;
  dFluxdFaceNormal->applyTranspose(dFluxes, dFaceNormal2dummy); 
  dFaceNormal2 += dFaceNormal2dummy;
// TODO: uncomment dFaceNormalVel2 part
  dFluxdFaceNormalVel->applyTranspose(dFluxes, dFaceNormalVel2dummy);
  dFaceNormalVel2 += dFaceNormalVel2dummy;
  dFluxes *= 3.0;
// TODO: assumed dUb is zero
//  dFluxdUb->applyTranspose(dFluxes, dUb);
//  fprintf(stderr, " ... norm of dUb is %e\n", dUb.norm());

//  dFaceNormal2 = dn;
//  dFaceNormalVel2 = dndot;

}

//------------------------------------------------------------------------------

template<int dim>
void FaceSet::computeDerivativeOperatorsOfFiniteVolumeTerm(
              FluxFcn **fluxFcn, BcData<dim> &bcData,
              GeoState &geoState, SVec<double,dim> &V,
              RectangularSparseMat<double,3,dim> &dFluxdFaceNormal,
              RectangularSparseMat<double,1,dim> &dFluxdFaceNormalVel,
              RectangularSparseMat<double,dim,dim> &dFluxdUb)
{
  
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();

  for(int i=0; i<numFaces; ++i)
    faces[i]->computeDerivativeOperatorsOfFiniteVolumeTerm(i, fluxFcn, n, ndot, V, Ub[i], dFluxdFaceNormal, dFluxdFaceNormalVel, dFluxdUb);



}

//------------------------------------------------------------------------------

template<int dim>
void FaceSet::computeFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
				      GeoState &geoState, SVec<double,dim> &V,
				      Vec<int> &fluidId, SVec<double,dim> &fluxes,
                                      LevelSetStructure *LSS)
{
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();
                          
  Vec<double>* UbHH = bcData.getBoundaryStateHH(); 

  for (int i=0; i<numFaces; ++i)  {
    if (!UbHH)
      faces[i]->computeFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], fluidId, fluxes, LSS);
    else
      faces[i]->computeFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], fluidId, fluxes, LSS,
                                        (*UbHH)[i]);
  }

}

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void FaceSet::computeFiniteVolumeTermLS(FluxFcn **fluxFcn, BcData<dim> &bcData,
					GeoState &geoState, SVec<double,dim> &V,
					SVec<double,dimLS> &Phi, SVec<double,dimLS> &PhiF)
{
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();

  for (int i=0; i<numFaces; ++i)  {
    faces[i]->computeFiniteVolumeTermLS(fluxFcn, n, ndot, V, Phi, PhiF);
  }
}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void FaceSet::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
					      GeoState &geoState, SVec<double,dim> &V,
					      GenMat<Scalar,neq> &A)
{
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();

  Vec<double>* UbHH = bcData.getBoundaryStateHH(); 

  for (int i=0; i<numFaces; ++i)  {
   
    if (!UbHH)
      faces[i]->computeJacobianFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], A);
    else
      faces[i]->computeJacobianFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], A,
                                                (*UbHH)[i]);
  }
}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void FaceSet::computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim> &riemann,
                                              FluxFcn **fluxFcn, BcData<dim> &bcData,
                                              GeoState &geoState, SVec<double,dim> &V,
                                              GenMat<Scalar,neq> &A)
{
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();

  Vec<double>* UbHH = bcData.getBoundaryStateHH();
  for (int i=0; i<numFaces; ++i) {

    if (!UbHH)
      faces[i]->computeJacobianFiniteVolumeTerm(riemann, fluxFcn, n, ndot, V, Ub[i], A);
    else
      faces[i]->computeJacobianFiniteVolumeTerm(riemann, fluxFcn, n, ndot, V, Ub[i], A,
                                                (*UbHH)[i]);
  }
}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void FaceSet::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
					      GeoState &geoState, SVec<double,dim> &V,
					      GenMat<Scalar,neq> &A, int* nodeType)
{
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();
  Vec<double>* UbHH = bcData.getBoundaryStateHH();
                                                                                                  
  for (int i=0; i<numFaces; ++i) {
    if (!UbHH)
      faces[i]->computeJacobianFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], A, nodeType);
    else
      faces[i]->computeJacobianFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], A,
                                                nodeType, (*UbHH)[i]);
  }
                                                                                                  
}

//------------------------------------------------------------------------------
//d2d emebedded
template<int dim, class Scalar, int neq>
void FaceSet::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
					      GeoState &geoState, SVec<double,dim> &V,
					      GenMat<Scalar,neq> &A, Vec<int> &fluidId,
                                              LevelSetStructure* LSS)
{
  
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();
  Vec<double>* UbHH = bcData.getBoundaryStateHH();

  for (int i=0; i<numFaces; ++i) {

    if (!UbHH)
      faces[i]->computeJacobianFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], A, fluidId,LSS);
    else
      faces[i]->computeJacobianFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], A, fluidId,LSS,
                                                (*UbHH)[i]);

  }

}
//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void FaceSet::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
                                              GeoState &geoState, SVec<double,dim> &V,
                                              GenMat<Scalar,neq> &A, Vec<int> &fluidId, 
                                              int* nodeType)
{
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  SVec<double,dim> &Ub = bcData.getFaceStateVector();
  Vec<double>* UbHH = bcData.getBoundaryStateHH();
                                                                                                  
  for (int i=0; i<numFaces; ++i) {
    
    if (!UbHH)
      faces[i]->computeJacobianFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], A, fluidId, nodeType);
    else
      faces[i]->computeJacobianFiniteVolumeTerm(fluxFcn, n, ndot, V, Ub[i], A, fluidId, nodeType,
                                                (*UbHH)[i]);
  }
                                                                                                  
}


//------------------------------------------------------------------------------

template<int dim, class Scalar, int dimLS>
void FaceSet::computeJacobianFiniteVolumeTermLS(GeoState &geoState, SVec<double,dim> &V,
						GenMat<Scalar,dimLS> &A)
{
  Vec<Vec3D> &n = geoState.getFaceNormal();
  Vec<double> &ndot = geoState.getFaceNormalVel();
  for (int i=0; i<numFaces; ++i)
    faces[i]->computeJacobianFiniteVolumeTermLS(n, ndot, V, A);

}

//------------------------------------------------------------------------------

template<int dim>
void FaceSet::computeGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, 
				  BcData<dim> &bcData, GeoState &geoState, 
				  SVec<double,3> &X, SVec<double,dim> &V, 
				  SVec<double,dim> &R, LevelSetStructure *LSS)
{

  SVec<double,dim> &Vwall = bcData.getFaceStateVector();
  Vec<double> &d2wall = geoState.getDistanceToWall();

	if (sampleMesh) {
		int i;
		for (int iFace=0; iFace<numSampledFaces; ++iFace) {
			i = facesConnectedToSampleNode[iFace];
			faces[i]->computeGalerkinTerm(elems, fet, X, d2wall, Vwall[i], V, R);
		}
	}
	else {
		for (int i=0; i<numFaces; ++i) 
    faces[i]->computeGalerkinTerm(elems, fet, X, d2wall, Vwall[i], V, R, LSS);
	}

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void FaceSet::computeDerivativeOfGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, BcData<dim> &bcData,
				  GeoState &geoState, SVec<double,3> &X, SVec<double,3> &dX,
				  SVec<double,dim> &V, SVec<double,dim> &dV, double dMach, SVec<double,dim> &dR)
{

  SVec<double,dim> &Vwall = bcData.getFaceStateVector();
  SVec<double,dim> &dVwall = bcData.getdFaceStateVector();
  Vec<double> &d2wall = geoState.getDistanceToWall();

  for (int i=0; i<numFaces; ++i)
    faces[i]->computeDerivativeOfGalerkinTerm(elems, fet, X, dX, d2wall, Vwall[i], dVwall[i], V, dV, dMach, dR);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void FaceSet::computeJacobianGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, 
					  BcData<dim> &bcData, GeoState &geoState, 
					  SVec<double,3> &X, Vec<double> &ctrlVol, 
					  SVec<double,dim> &V, GenMat<Scalar,neq> &A)
{

  SVec<double,dim> &Vwall = bcData.getFaceStateVector();
  Vec<double> &d2wall = geoState.getDistanceToWall();

  for (int i=0; i<numFaces; ++i) 
    faces[i]->computeJacobianGalerkinTerm(elems, fet, X, ctrlVol, d2wall, Vwall[i], V, A);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void FaceSet::computeBCsJacobianWallValues(ElemSet &elems, FemEquationTerm *fet, 
					  BcData<dim> &bcData, GeoState &geoState, 
					  SVec<double,3> &X, SVec<double,dim> &V)
{

  SVec<double,dim> &dVwallface = bcData.getdFaceStateVectorSA();

  SVec<double,dim> &Vwall = bcData.getFaceStateVector();
  Vec<double> &d2wall = geoState.getDistanceToWall();

  for (int i=0; i<numFaces; ++i) 
    faces[i]->computeBCsJacobianWallValues(elems, fet, X, d2wall, Vwall[i], dVwallface[i], V);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void FaceSet::computeHHBoundaryTermResidual(BcData<dim> &bcData,SVec<double,dim> &V, Vec<double> &res,VarFcn* vf) {

  SVec<double,dim> &Ub = bcData.getFaceStateVector();
  Vec<double>& UbHH = *bcData.getBoundaryStateHH();
  for (int i=0; i<numFaces; ++i) 
    faces[i]->computeHHBoundaryTermResidual(V, Ub[i], UbHH[i],res[i], vf);
}

template<class Scalar,int dim,int neq>
inline
void FaceSet::computeHHBoundaryTermJacobian(FluxFcn **fluxFcn, BcData<dim> &bcData,SVec<double,dim> &U,
					    GeoState& geoState,
					    GenMat<Scalar,neq> &A, VarFcn* vf) {

  SVec<double,dim> &Ub = bcData.getFaceStateVector();
  Vec<double>& UbHH = *bcData.getBoundaryStateHH();

  Vec<Vec3D> &normals = geoState.getFaceNormal();
  Vec<double> &normalVel = geoState.getFaceNormalVel();
  
  for (int i=0; i<numFaces; ++i) 
    faces[i]->computeHHBoundaryTermJacobian(i,fluxFcn,U, Ub[i], A,vf,
                                            UbHH[i],normals, normalVel);
}

template<class Scalar, int dim>
void  FaceSet::computeMatVecProdH1FarFieldHH(GenMat<Scalar,dim> &A, SVec<double,dim> &p_u,
					     SVec<double,dim> &prod_u,Vec<double>& p_hh, 
					     Vec<double>& prod_hh) {

  
  for (int i=0; i<numFaces; ++i) 
    faces[i]->computeMatVecProdH1FarFieldHH(i,A,p_u,prod_u, p_hh[i],prod_hh[i]);
}
