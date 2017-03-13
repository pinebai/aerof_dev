#include <AgglomeratedFace.h>

#include "NavierStokesTerm.h"

#include "FemEquationTerm.h"

AgglomeratedFace::AgglomeratedFace() : node(0), code(0), normal(0.0), area(0.0),masterFlag(true) {

}

AgglomeratedFace::AgglomeratedFace(int node, int code) : node(node), code(code), normal(0.0),
                                    area(0.0),masterFlag(true) {

//  if (code == BC_INLET_MOVING || code == BC_INLET_FIXED ||
//      code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED)
//    normal = -1.0*normal;
}

AgglomeratedFace::AgglomeratedFace(const AgglomeratedFace& oth) :
  node(oth.node), code(oth.code), normal(oth.normal),masterFlag(oth.masterFlag) {
/*
  if (code == BC_INLET_MOVING || code == BC_INLET_FIXED ||
      code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED)
    normal = -1.0*normal;

*/
}

AgglomeratedFace::~AgglomeratedFace() {

}

template<int dim>
void AgglomeratedFace::assignFreeStreamValues2(SVec<double,dim> &Uin, SVec<double,dim> &Uout, double *U)
{
  int k, j;

  if (code == BC_INLET_MOVING || code == BC_INLET_FIXED ||
      code == BC_DIRECTSTATE_INLET_MOVING || code == BC_DIRECTSTATE_INLET_FIXED ||
      code == BC_MASSFLOW_INLET_MOVING || code == BC_MASSFLOW_INLET_FIXED)
    for (k=0; k<dim; ++k) {
      U[k] = Uin[node][k];
    }
  else if (code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED ||
           code == BC_DIRECTSTATE_OUTLET_MOVING || code == BC_DIRECTSTATE_OUTLET_FIXED ||
           code == BC_MASSFLOW_OUTLET_MOVING || code == BC_MASSFLOW_OUTLET_FIXED)
    for (k=0; k<dim; ++k) {
      U[k] = Uout[node][k];
    }
  else
    for (k=0; k<dim; ++k)
      U[k] = 0.0;
}


template<int dim>
void AgglomeratedFace::computeFiniteVolumeTerm(FluxFcn **fluxFcn, 
	                                       SVec<double,dim> &V, 
			                       double *Ub, SVec<double,dim> &fluxes) {
  if (!masterFlag)
    return;

  if(fluxFcn[code]){
    double flux[dim];
    fluxFcn[code]->compute(0.0, 0.0,  normal, 0.0, 
  		           V[node], Ub, flux);
    for (int k=0; k<dim; ++k){
      fluxes[ node ][k] += flux[k];
    } 
  }
}

template<int dim>
void AgglomeratedFace::computeFiniteVolumeTerm(FluxFcn **fluxFcn, 
	                                       SVec<double,dim> &V, 
			                       double *Ub, SVec<double,dim> &fluxes,
					       LevelSetStructure& LSS) {
  if (!masterFlag)
    return;
  
  if (!LSS.isActive(0.0, node))
    return;
  
  Vec3D n2 = normal;
   
  if(fluxFcn[code]){
    double flux[dim];
    fluxFcn[code]->compute(0.0, 0.0,  n2, 0.0, 
  		           V[node], Ub, flux);
//    if (code == BC_INLET_MOVING || code == BC_INLET_FIXED ||
//        code == BC_OUTLET_MOVING || code == BC_OUTLET_FIXED)
//      std::cout << Ub[0] << " " << Ub[1] << " " << Ub[3] <<  " " << V[node][0] << " "  << V[node][1] << " " << V[node][3] << std::endl;
 
    for (int k=0; k<dim; ++k){
      fluxes[ node ][k] += flux[k]; //= normal[0]+normal[1]+normal[2];//flux[k];
    } 
  }
}

template<int dim>
void AgglomeratedFace::computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn,
			   SVec<double,3> &X, SVec<double,dim> &V, 
			   Vec<double> &idti, Vec<double> &idtv,
			   TimeLowMachPrec &tprec)
{
  if (!masterFlag)
    return;

  double S = sqrt(normal * normal);
  double invS = 1.0 / S;

  Vec3D n = invS * normal;
  double ndot = 0.0;//invS * getNormalVel(normalVel);

  Vec3D u = varFcn->getVelocity(V[ node ]);
  double a = varFcn->computeSoundSpeed(V[ node]);
  double un = u * n - ndot;

  // Low-Mach Preconditioner
  double locMach = varFcn->computeMachNumber(V[ node ]);
  double locbeta = tprec.getBeta(locMach,true);
  double beta2 = locbeta * locbeta;
  double coeff1 = fabs((1.0+beta2)*un);
  double coeff2 = pow(pow((1.0-beta2)*un,2.0) + pow(2.0*locbeta*a,2.0),0.5);

  idti[ node ] += 0.5*(coeff1+coeff2)* area;
    
}


void AgglomeratedFace::setNodeType(int* priority, int* nodeType) {
 
  if (priority[code] > priority[ nodeType[ node ] ])
    nodeType[ node ] = code;
}

AgglomeratedFaceSet::AgglomeratedFaceSet(int size) : numFaces(size) {

  myFaces = new AgglomeratedFace[numFaces];
}

AgglomeratedFaceSet::~AgglomeratedFaceSet() {

  delete [] myFaces;
}

template<int dim, class Scalar, int neq>
void AgglomeratedFace::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, 
			 	                       SVec<double,dim> &V, 
 				                       double *Ub, GenMat<Scalar,neq> &A) {
  if (!masterFlag)
    return;

  double jac[neq*neq];

  double normVel= 0.0;

  fluxFcn[code]->computeJacobian(1.0, 0.0, normal, normVel, V[node], Ub, jac);
  Scalar *Aii = A.getElem_ii(node);
  for (int k=0; k<neq*neq; ++k) 
    Aii[k] += jac[k];

}
template<int dim, class Scalar, int neq>
void AgglomeratedFace::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, 
			 	                       SVec<double,dim> &V, 
 				                       double *Ub, GenMat<Scalar,neq> &A,
						       LevelSetStructure& LSS) {
  if (!masterFlag)
    return;

  if (!LSS.isActive(0.0, node))
    return;

  double jac[neq*neq];

  double normVel= 0.0;
  Vec3D n2 = normal;
  fluxFcn[code]->computeJacobian(1.0, 0.0, n2, normVel, V[node], Ub, jac);
  Scalar *Aii = A.getElem_ii(node);
  for (int k=0; k<neq*neq; ++k) 
    Aii[k] += jac[k];

}

template<int dim>
void AgglomeratedFace::
computeThinLayerViscousFiniteVolumeTerm(class FemEquationTerm* fet,
                                        VarFcn* varFcn,
                                        SVec<double,dim> &V, 
                                        SVec<double,dim> &dX,
                                        SVec<double,dim> &dY,
                                        SVec<double,dim> &dZ,
                                        Vec<double>& d2wall,
                                        double* Vwall,
			                SVec<double,dim> &fluxes) {
  if (!masterFlag)
    return;

  if (code == BC_SYMMETRY) return;

  if (!fet->doesFaceTermExist(code))
    return;

  if (fet->doesFaceNeedGradientP1Function())
    return;

  double d2w[3] = {d2wall[node], d2wall[node],d2wall[node]};
  double *v[3] = {V[node], V[node], V[node]};
  double r[dim];
  fet->computeSurfaceTerm(code, normal, d2w, Vwall, v, r);

  for (int k = 0; k < dim; ++k) 
    fluxes[node][k] -= r[k];  

/*
  double ooreynolds_mu = ns->get_ooreynolds_mu(); 
  double area = normal.norm();
  
  Vec3D xhat = normal;
  xhat /= area;

  double Tcg = varFcn->computeTemperature(V[node]);

  double Tg[5];
  varFcn->computeTemperatureGradient(V[node],Tg);

  double mu     = ns->getViscoFcn()->compute_mu(Tcg);
  double lambda = ns->getViscoFcn()->compute_lambda(Tcg,mu);
  double kappa  = ns->getThermalCondFcn()->compute(Tcg);

  mu     *= ooreynolds_mu;
  lambda *= ooreynolds_mu;
  kappa  *= ooreynolds_mu;

  double uhat = V[node][1]*xhat[0]+V[node][2]*xhat[1]+V[node][3]*xhat[2];

  // Not necessarily p
  double dpdxhat = dX[node][4]*xhat[0]+dY[node][4]*xhat[1]+dZ[node][4]*xhat[2];
  // Not necessarily rho
  double drhodxhat = dX[node][0]*xhat[0]+dY[node][0]*xhat[1]+dZ[node][0]*xhat[2];

  double dudxhat[3] = {dX[node][1]*xhat[0]+dY[node][1]*xhat[1]+dZ[node][1]*xhat[2],
                       dX[node][2]*xhat[0]+dY[node][2]*xhat[1]+dZ[node][2]*xhat[2],
                       dX[node][3]*xhat[0]+dY[node][3]*xhat[1]+dZ[node][3]*xhat[2] };

  double duhatdx = dX[node][1]*xhat[0]+dX[node][2]*xhat[1]+dX[node][3]*xhat[2];
  double duhatdy = dY[node][1]*xhat[0]+dY[node][2]*xhat[1]+dY[node][3]*xhat[2]; 
  double duhatdz = dZ[node][1]*xhat[0]+dZ[node][2]*xhat[1]+dZ[node][3]*xhat[2]; 
  
  double duhatdxhat = duhatdx*xhat[0]+duhatdy*xhat[1]+duhatdz*xhat[2];

  double dTdxhat = drhodxhat*Tg[0]+dudxhat[0]*Tg[1]+dudxhat[1]*Tg[2]+
                   dudxhat[2]*Tg[3]+dpdxhat*Tg[4];
 
  //double velnorm = sqrt(V[node][1]*V[node][1]+V[node][2]*V[node][2]+V[node][3]*V[node][3]);

  double flux[5];
  double f1 = (lambda+2.0*mu)*duhatdxhat;
  flux[0] = 0.0;
  for (int k = 0; k < 3; ++k) {

    flux[k+1] = xhat[k]*f1+mu*(dudxhat[k]-duhatdxhat*xhat[k]);
  }
  
  flux[4] = uhat*f1+mu*(V[node][1]*dudxhat[0]+V[node][2]*dudxhat[1]+V[node][3]*dudxhat[2] -
                        uhat*duhatdxhat)-kappa*dTdxhat ;

  for (int k = 0; k < 5; ++k) {

    fluxes[node][k] += flux[k]*area;
  }
*/
}

template<int dim,class Scalar,int neq>
void AgglomeratedFace::
computeJacobianThinLayerViscousFiniteVolumeTerm(class FemEquationTerm* fet,
                                        VarFcn* varFcn,
                                        SVec<double,dim> &V, 
                                        SVec<double,dim> &dX,
                                        SVec<double,dim> &dY,
                                        SVec<double,dim> &dZ,
/*			                SVec<double,neq*neq> &JacX,
                                        SVec<double,neq*neq> &JacY,
                                        SVec<double,neq*neq> &JacZ, 
*/
                                        Vec<double>& ctrlVol,
                                        Vec<double>& d2wall,
                                        double* Vwall,
                                        GenMat<Scalar,neq>& A) {
  if (!masterFlag)
    return;

  if (code == BC_SYMMETRY) return;

  if (fet->doesFaceNeedGradientP1Function())
    return;

  double d2w[3] = {d2wall[node], d2wall[node],d2wall[node]};
  double *v[3] = {V[node], V[node], V[node]};
  
  double dRdU[3][neq*neq];

  fet->computeJacobianSurfaceTerm(code, normal, d2w, Vwall, v, reinterpret_cast<double *>(dRdU));

  Scalar *Aii = A.getElem_ii(node);
  for (int m=0; m<neq*neq; ++m)
    Aii[m] -= dRdU[0][m];  

/*
  Scalar* jac = A.getElem_ii(node);
  
  Scalar* jacx = JacX[node];
  Scalar* jacy = JacY[node];
  Scalar* jacz = JacZ[node];

  //memset(jac,0,sizeof(double)*neq*neq);
  memset(jacx,0,sizeof(Scalar)*neq*neq);
  memset(jacy,0,sizeof(Scalar)*neq*neq);
  memset(jacz,0,sizeof(Scalar)*neq*neq);

  double ooreynolds_mu = ns->get_ooreynolds_mu(); 
  double area = normal.norm();
  
  Vec3D xhat = normal;
  xhat /= area;

  double Tcg = varFcn->computeTemperature(V[node]);

  double Trr,Trp,Tpp;

  double Tg[5];
  varFcn->computeTemperatureGradient(V[node],Tg);
  varFcn->computeTemperatureHessian(V[node],Trr,Trp,Tpp);

  double mu     = ns->getViscoFcn()->compute_mu(Tcg);
  double lambda = ns->getViscoFcn()->compute_lambda(Tcg,mu);
  double kappa  = ns->getThermalCondFcn()->compute(Tcg);

  mu     *= ooreynolds_mu;
  lambda *= ooreynolds_mu;
  kappa  *= ooreynolds_mu;

  double uhat = V[node][1]*xhat[0]+V[node][2]*xhat[1]+V[node][3]*xhat[2];

  // Not necessarily p
  double dpdxhat = dX[node][4]*xhat[0]+dY[node][4]*xhat[1]+dZ[node][4]*xhat[2];
  // Not necessarily rho
  double drhodxhat = dX[node][0]*xhat[0]+dY[node][0]*xhat[1]+dZ[node][0]*xhat[2];

  double dudxhat[3] = {dX[node][1]*xhat[0]+dY[node][1]*xhat[1]+dZ[node][1]*xhat[2],
                       dX[node][2]*xhat[0]+dY[node][2]*xhat[1]+dZ[node][2]*xhat[2],
                       dX[node][3]*xhat[0]+dY[node][3]*xhat[1]+dZ[node][3]*xhat[2] };

  double duhatdx = dX[node][1]*xhat[0]+dX[node][2]*xhat[1]+dX[node][3]*xhat[2];
  double duhatdy = dY[node][1]*xhat[0]+dY[node][2]*xhat[1]+dY[node][3]*xhat[2]; 
  double duhatdz = dZ[node][1]*xhat[0]+dZ[node][2]*xhat[1]+dZ[node][3]*xhat[2]; 
  
  double duhatdxhat = duhatdx*xhat[0]+duhatdy*xhat[1]+duhatdz*xhat[2];

  double dTdxhat = drhodxhat*Tg[0]+dudxhat[0]*Tg[1]+dudxhat[1]*Tg[2]+
                   dudxhat[2]*Tg[3]+dpdxhat*Tg[4];
 
  //double velnorm = sqrt(V[node][1]*V[node][1]+V[node][2]*V[node][2]+V[node][3]*V[node][3]);

  double flux[5];
  double f1 = (lambda+2.0*mu)*duhatdxhat;
  flux[0] = 0.0;
  double xhatsum = xhat[0]+xhat[1]+xhat[2];
  for (int k = 0; k < 3; ++k) {

    for (int l = 0; l < 3; ++l) {
      jacx[(k+1)*neq+l+1] = xhat[k]*(lambda+2.0*mu)*xhat[l]*xhat[0];
      jacy[(k+1)*neq+l+1] = xhat[k]*(lambda+2.0*mu)*xhat[l]*xhat[1];
      jacz[(k+1)*neq+l+1] = xhat[k]*(lambda+2.0*mu)*xhat[l]*xhat[2];
  
      jacx[(k+1)*neq+l+1] += mu*((k==l?xhat[0]:0.0)-(lambda+2.0*mu)*xhat[l]*xhat[0]*xhat[k]);
      jacy[(k+1)*neq+l+1] += mu*((k==l?xhat[1]:0.0)-(lambda+2.0*mu)*xhat[l]*xhat[1]*xhat[k]);
      jacz[(k+1)*neq+l+1] += mu*((k==l?xhat[2]:0.0)-(lambda+2.0*mu)*xhat[l]*xhat[2]*xhat[k]);
    } 
    //flux[k+1] = xhat[k]*f1+mu*(dudxhat[k]-duhatdxhat*xhat[k]);
  }

  for (int l = 0; l < 3; ++l) {
    jacx[4*neq+l+1] = ((lambda+2.0*mu)*xhat[l]*xhat[0])*uhat;
    jacy[4*neq+l+1] = ((lambda+2.0*mu)*xhat[l]*xhat[1])*uhat;
    jacz[4*neq+l+1] = ((lambda+2.0*mu)*xhat[l]*xhat[2])*uhat;

    jac[4*neq+l+1] += xhat[l]*f1*area;

    jac[4*neq+l+1] += mu*(dudxhat[l]-xhat[l]*duhatdxhat)*area;

    jacx[4*neq+l+1] += mu*(V[node][1]*(l==0?xhat[0]:0.0) - uhat*(lambda+2.0*mu)*xhat[l]*xhat[0]);
    jacy[4*neq+l+1] += mu*(V[node][2]*(l==1?xhat[1]:0.0) - uhat*(lambda+2.0*mu)*xhat[l]*xhat[1]);
    jacz[4*neq+l+1] += mu*(V[node][3]*(l==2?xhat[2]:0.0) - uhat*(lambda+2.0*mu)*xhat[l]*xhat[2]);

    jacx[4*neq] = -kappa*Tg[0]*xhat[0];
    jacx[4*neq+4] = -kappa*Tg[4]*xhat[0];
    jacy[4*neq] = -kappa*Tg[0]*xhat[1];
    jacy[4*neq+4] = -kappa*Tg[4]*xhat[1];
    jacz[4*neq] = -kappa*Tg[0]*xhat[2];
    jacz[4*neq+4] = -kappa*Tg[4]*xhat[2];

    jac[4*neq] += -kappa*(drhodxhat*Trr+dpdxhat*Trp)*area;
    jac[4*neq+4] += -kappa*(drhodxhat*Trp+dpdxhat*Tpp)*area;
  } 
   
  //flux[4] = uhat*f1+mu*(V[node][1]*dudxhat[0]+V[node][2]*dudxhat[1]+V[node][3]*dudxhat[2] -
  //                      uhat*duhatdxhat)-kappa*dTdxhat ;

  for (int k = 0; k < neq*neq; ++k) {

    jacx[k] *= area;
    jacy[k] *= area;
    jacz[k] *= area;
  }

  double tmpx[neq*neq],tmpy[neq*neq],tmpz[neq*neq];
  varFcn->postMultiplyBydVdU(V[node], jacx,tmpx);
  varFcn->postMultiplyBydVdU(V[node], jacy,tmpy);
  varFcn->postMultiplyBydVdU(V[node], jacz,tmpz);
  
  double invvol = 1.0/ctrlVol[node];
  for (int k = 0; k < neq*neq; ++k) {

    jac[k] += invvol*(normal[0]*jacx[k]+normal[1]*jacy[k]+normal[2]*jacz[k]);
  }
*/
}

template<int dim>
void AgglomeratedFaceSet::computeFiniteVolumeTerm(FluxFcn **fluxFcn, 
                                               SVec<double,dim> &V, 
			                       SVec<double,dim>& Ub, SVec<double,dim> &fluxes) {

  for (int i = 0; i < numFaces; ++i) {

    myFaces[i].computeFiniteVolumeTerm(fluxFcn, V, Ub[i], fluxes);
  }
}

template<int dim, class Scalar, int neq>
void AgglomeratedFaceSet::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn,
			                               SVec<double,dim> &V, 
				                       SVec<double,dim>& Ub, GenMat<Scalar,neq> &A) {

  for (int i = 0; i < numFaces; ++i) {
 
    myFaces[i].computeJacobianFiniteVolumeTerm(fluxFcn,V, Ub[i], A);
  }
}

template<int dim>
void AgglomeratedFaceSet::
computeFiniteVolumeTerm(FluxFcn **fluxFcn, 
			SVec<double,dim> &V, 
			SVec<double,dim>& Ub, SVec<double,dim> &fluxes,
			LevelSetStructure& LSS) {
  
  for (int i = 0; i < numFaces; ++i) {

    myFaces[i].computeFiniteVolumeTerm(fluxFcn, V, Ub[i], fluxes, LSS);
  }
}

template<int dim, class Scalar, int neq>
void AgglomeratedFaceSet::
computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn,
				SVec<double,dim> &V, 
				SVec<double,dim>& Ub, GenMat<Scalar,neq> &A,
				LevelSetStructure& LSS) {
  
  for (int i = 0; i < numFaces; ++i) {
 
    myFaces[i].computeJacobianFiniteVolumeTerm(fluxFcn,V, Ub[i], A, LSS);
  }
}

template<int dim>
void AgglomeratedFaceSet::
computeThinLayerViscousFiniteVolumeTerm(class FemEquationTerm* ns,
                                        VarFcn* varFcn,
                                        SVec<double,dim> &V, 
                                        SVec<double,dim> &dX,
                                        SVec<double,dim> &dY,
                                        SVec<double,dim> &dZ,
                                        Vec<double>& d2wall,
                                        SVec<double,dim>& Vwall, 
			                SVec<double,dim> &fluxes) {

  for (int i = 0; i < numFaces; ++i) {
 
    myFaces[i].computeThinLayerViscousFiniteVolumeTerm(ns,varFcn, V,dX,dY, dZ,
                                                       d2wall, Vwall[i],fluxes);
  }
}

template<int dim,class Scalar,int neq>
void AgglomeratedFaceSet::
computeJacobianThinLayerViscousFiniteVolumeTerm(class FemEquationTerm* fet,
                                        VarFcn* varFcn,
                                        SVec<double,dim> &V, 
                                        SVec<double,dim> &dX,
                                        SVec<double,dim> &dY,
                                        SVec<double,dim> &dZ,
/*			                SVec<double,neq*neq> &JacX,
                                        SVec<double,neq*neq> &JacY,
                                        SVec<double,neq*neq> &JacZ, 
*/
                                        Vec<double>& ctrlVol,
                                        Vec<double>& d2wall,
                                        SVec<double,dim>& Vwall,
                                        GenMat<Scalar,neq>& A) {

  for (int i = 0; i < numFaces; ++i) {

    myFaces[i].computeJacobianThinLayerViscousFiniteVolumeTerm(fet,
      varFcn,V,dX,dY,dZ/*,JacX,JacY,JacZ*/, ctrlVol, d2wall, Vwall[i], A);
  }
}

template<int dim>
void AgglomeratedFaceSet::computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn,
			   SVec<double,3> &X, SVec<double,dim> &V, 
			   Vec<double> &idti, Vec<double> &idtv,
			   TimeLowMachPrec &tprec)
{
  for (int i = 0; i < numFaces; ++i) {
   
    myFaces[i].computeTimeStep(fet,varFcn,X,V,idti,idtv,tprec);
  } 
}



#define INST_HELPER(dim) \
template void AgglomeratedFace::computeFiniteVolumeTerm(FluxFcn **fluxFcn, \
         SVec<double,dim> &V, double* Ub, SVec<double,dim> &fluxes); \
template void AgglomeratedFaceSet::computeFiniteVolumeTerm(FluxFcn **fluxFcn, \
         SVec<double,dim> &V, SVec<double,dim>& Ub, SVec<double,dim> &fluxes); \
template void AgglomeratedFace::computeFiniteVolumeTerm(FluxFcn **fluxFcn, \
							SVec<double,dim> &V, double* Ub, SVec<double,dim> &fluxes, LevelSetStructure&); \
template void AgglomeratedFaceSet::computeFiniteVolumeTerm(FluxFcn **fluxFcn, \
         SVec<double,dim> &V, SVec<double,dim>& Ub, SVec<double,dim> &fluxes, LevelSetStructure&); \
template void AgglomeratedFace::assignFreeStreamValues2(SVec<double,dim> &Uin, SVec<double,dim> &Uout, double *U); \
  template void AgglomeratedFace::computeThinLayerViscousFiniteVolumeTerm(class FemEquationTerm*, \
                                        VarFcn* varFcn, \
                                        SVec<double,dim> &V,  \
                                        SVec<double,dim> &dX, \
                                        SVec<double,dim> &dY, \
                                        SVec<double,dim> &dZ, \
                                        Vec<double>& d2wall,\
                                        double* Vwall, \
			                SVec<double,dim> &fluxes);\
  template void AgglomeratedFaceSet::computeThinLayerViscousFiniteVolumeTerm(class FemEquationTerm*, \
                                        VarFcn* varFcn, \
                                        SVec<double,dim> &V,  \
                                        SVec<double,dim> &dX, \
                                        SVec<double,dim> &dY, \
                                        SVec<double,dim> &dZ, \
                                        Vec<double>& d2wall,\
                                        SVec<double,dim>& Vwall, \
		                        SVec<double,dim> &fluxes); \
template void AgglomeratedFace::computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn, \
			   SVec<double,3> &X, SVec<double,dim> &V,  \
			   Vec<double> &idti, Vec<double> &idtv, \
			   TimeLowMachPrec &tprec); \
template void AgglomeratedFaceSet::computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn, \
			   SVec<double,3> &X, SVec<double,dim> &V,  \
			   Vec<double> &idti, Vec<double> &idtv, \
			   TimeLowMachPrec &tprec);

#define INST_HELPER2(dim,Scalar,neq) \
template void AgglomeratedFace::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, \
                         SVec<double,dim> &V, double *Ub, GenMat<Scalar,neq> &A); \
template void AgglomeratedFace::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, \
								SVec<double,dim> &V, double *Ub, GenMat<Scalar,neq> &A, LevelSetStructure&); \
template  void AgglomeratedFace:: \
  computeJacobianThinLayerViscousFiniteVolumeTerm(class FemEquationTerm* ns,\
                                        VarFcn* varFcn, \
                                        SVec<double,dim> &V,  \
                                        SVec<double,dim> &dX, \
                                        SVec<double,dim> &dY, \
                                        SVec<double,dim> &dZ,\
                                        Vec<double>& ctrlVol, \
                                        Vec<double>& d2wall, \
                                        double* Vwall, \
                                        GenMat<Scalar,neq>& A);\
template  void AgglomeratedFaceSet:: \
  computeJacobianThinLayerViscousFiniteVolumeTerm(class FemEquationTerm* ns,\
                                        VarFcn* varFcn, \
                                        SVec<double,dim> &V,  \
                                        SVec<double,dim> &dX, \
                                        SVec<double,dim> &dY, \
                                        SVec<double,dim> &dZ,\
                                        Vec<double>& ctrlVol, \
                                        Vec<double>& d2wall, \
                                        SVec<double,dim>& Vwall, \
                                        GenMat<Scalar,neq>& A);\
template void AgglomeratedFaceSet::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, \
                  SVec<double,dim> &V, SVec<double,dim>& Ub, GenMat<Scalar,neq> &A);\
template void AgglomeratedFaceSet::computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, \
								   SVec<double,dim> &V, SVec<double,dim>& Ub, GenMat<Scalar,neq> &A,LevelSetStructure&);


INST_HELPER(1);
INST_HELPER(2);
INST_HELPER(3);
INST_HELPER(5);
INST_HELPER(6);
INST_HELPER(7);

INST_HELPER2(1,double,1);
INST_HELPER2(2,double,2);
INST_HELPER2(5,double,5);
INST_HELPER2(6,double,6);
INST_HELPER2(7,double,5);
INST_HELPER2(6,double,5);
INST_HELPER2(7,double,7);

INST_HELPER2(6,double,1);
INST_HELPER2(7,double,2);
