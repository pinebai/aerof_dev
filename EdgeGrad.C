#include <EdgeGrad.h>
#include <IoData.h>
#include <Elem.h>
#include <SubDomain.h>
#include <Vector.h>
#include <Vector3D.h>
#include "LevelSet/LevelSetStructure.h"
#include "LevelSet/MultiGridLevelSetStructure.h"
#include <cmath>

//------------------------------------------------------------------------------

template<int dim>
EdgeGrad<dim>::EdgeGrad(IoData& iod)
{

  v6data = 0;

  beta = iod.schemes.ns.beta;

  tag = 0;
  
  if (iod.schemes.ns.gradient == SchemeData::NON_NODAL){
    xiu = 0.0;
    xic = 0.0;
  }
  
  if(iod.schemes.ns.dissipation == SchemeData::SIXTH_ORDER) {
    xiu = iod.schemes.ns.xiu;
    xic = iod.schemes.ns.xic;
    if (beta != 0.0) {
      xiu /= beta;
      xic /= beta;
    }   
  }    

}

//------------------------------------------------------------------------------

template<int dim>
EdgeGrad<dim>::~EdgeGrad()
{
  
  if (v6data) delete [] v6data;

}

//------------------------------------------------------------------------------

template<int dim>
void EdgeGrad<dim>::computeUpwindGradient(Elem& elem, double rij[3], SVec<double,3>& X,
					  SVec<double,dim>& V, double* grad)
{

  double dp1dxi[4][3];
  elem.computeGradientP1Function(X, dp1dxi);
  
  int i0 = elem.nodeNum(0);
  int i1 = elem.nodeNum(1);
  int i2 = elem.nodeNum(2);
  int i3 = elem.nodeNum(3);

  for (int k=0; k<dim; ++k) {
    grad[k] = ( (dp1dxi[0][0]*V[i0][k] + dp1dxi[1][0]*V[i1][k] + 
		 dp1dxi[2][0]*V[i2][k] + dp1dxi[3][0]*V[i3][k]) * rij[0] +
		(dp1dxi[0][1]*V[i0][k] + dp1dxi[1][1]*V[i1][k] + 
		 dp1dxi[2][1]*V[i2][k] + dp1dxi[3][1]*V[i3][k]) * rij[1] +
		(dp1dxi[0][2]*V[i0][k] + dp1dxi[1][2]*V[i1][k] + 
		 dp1dxi[2][2]*V[i2][k] + dp1dxi[3][2]*V[i3][k]) * rij[2] );
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void EdgeGrad<dim>::computeDerivativeOfUpwindGradient(Elem& elem, double rij[3], double drij[3], SVec<double,3>& X, SVec<double,3>& dX,
					  SVec<double,dim>& V, SVec<double,dim>& dV, double* dGrad)
{

  double dp1dxi[4][3];
  elem.computeGradientP1Function(X, dp1dxi);

  double ddp1dxi[4][3];
  elem.computeDerivativeOfGradientP1Function(X, dX, ddp1dxi);

  int i0 = elem.nodeNum(0);
  int i1 = elem.nodeNum(1);
  int i2 = elem.nodeNum(2);
  int i3 = elem.nodeNum(3);

  for (int k=0; k<dim; ++k) {
    dGrad[k] = ( (ddp1dxi[0][0]*V[i0][k] + dp1dxi[0][0]*dV[i0][k] + ddp1dxi[1][0]*V[i1][k] + dp1dxi[1][0]*dV[i1][k] +
		 ddp1dxi[2][0]*V[i2][k] + dp1dxi[2][0]*dV[i2][k] + ddp1dxi[3][0]*V[i3][k] + dp1dxi[3][0]*dV[i3][k]) * rij[0] +
		(dp1dxi[0][0]*V[i0][k] + dp1dxi[1][0]*V[i1][k] +
		 dp1dxi[2][0]*V[i2][k] + dp1dxi[3][0]*V[i3][k]) * drij[0] +
		(ddp1dxi[0][1]*V[i0][k] + dp1dxi[0][1]*dV[i0][k] + ddp1dxi[1][1]*V[i1][k] + dp1dxi[1][1]*dV[i1][k] +
		 ddp1dxi[2][1]*V[i2][k] + dp1dxi[2][1]*dV[i2][k] + ddp1dxi[3][1]*V[i3][k] + dp1dxi[3][1]*dV[i3][k]) * rij[1] +
		(dp1dxi[0][1]*V[i0][k] + dp1dxi[1][1]*V[i1][k] +
		 dp1dxi[2][1]*V[i2][k] + dp1dxi[3][1]*V[i3][k]) * drij[1] +
		(ddp1dxi[0][2]*V[i0][k] + dp1dxi[0][2]*dV[i0][k] + ddp1dxi[1][2]*V[i1][k] + dp1dxi[1][2]*dV[i1][k] +
		 ddp1dxi[2][2]*V[i2][k] + dp1dxi[2][2]*dV[i2][k] + ddp1dxi[3][2]*V[i3][k] + dp1dxi[3][2]*dV[i3][k]) * rij[2] +
		(dp1dxi[0][2]*V[i0][k] + dp1dxi[1][2]*V[i1][k] +
		 dp1dxi[2][2]*V[i2][k] + dp1dxi[3][2]*V[i3][k]) * drij[2] );
  }

}

//------------------------------------------------------------------------------


//------------------------------------------------------------------------------

template<int dim>
void EdgeGrad<dim>::computeUpwindGradient(Elem& elem, double rij[3], SVec<double,3>& X, SVec<double,dim>& V,
					  double* grad, bool* v6_flag,
					  Vec<int>& fluidId, LevelSetStructure &LSS)
{
  
  // Check if all the nodes of the element are "valid".
  // If not, flag the element and the gradient is set to 0
  //---------------------------------------------------
  int le, N1, N2;
  bool e_flag = false;

  for(int j=0; j<6; j++)
  {
    le = elem.edgeNum(j);

    N1 = elem.nodeNum( elem.edgeEnd(j,0) );
    N2 = elem.nodeNum( elem.edgeEnd(j,1) );
   
	 bool isValid = true;

	 if(fluidId[N1] != fluidId[N2]) isValid = false;

	 if(LSS.edgeIntersectsStructure(0.0, le)) isValid = false;
	 if(!LSS.isActive(0.0, N1) || !LSS.isActive(0.0, N2)) isValid = false;
	 
    if(!isValid)
	 {
    	 e_flag = true;
    	 break;
      }
            
  }
  //---------------------------------------------------

  if(!e_flag)
  {
    *v6_flag = false;

    double dp1dxi[4][3];
    elem.computeGradientP1Function(X, dp1dxi);
  
    int i0 = elem.nodeNum(0);
    int i1 = elem.nodeNum(1);
    int i2 = elem.nodeNum(2);
    int i3 = elem.nodeNum(3);

    for(int k=0; k<dim; ++k) 
	 {
      grad[k] = ( (dp1dxi[0][0]*V[i0][k] + dp1dxi[1][0]*V[i1][k] + 
		   dp1dxi[2][0]*V[i2][k] + dp1dxi[3][0]*V[i3][k]) * rij[0] +
		  (dp1dxi[0][1]*V[i0][k] + dp1dxi[1][1]*V[i1][k] + 
		   dp1dxi[2][1]*V[i2][k] + dp1dxi[3][1]*V[i3][k]) * rij[1] +
		  (dp1dxi[0][2]*V[i0][k] + dp1dxi[1][2]*V[i1][k] + 
		   dp1dxi[2][2]*V[i2][k] + dp1dxi[3][2]*V[i3][k]) * rij[2] );
    }
  }
  else
  {
    *v6_flag = true;

    for (int k=0; k<dim; ++k) grad[k] = 0.0;
  }

}


//------------------------------------------------------------------------------

template<int dim>
void EdgeGrad<dim>::computeUpwindGradient(Elem& elem, double rij[3], SVec<double,3>& X, SVec<double,dim>& V,
					  double* grad, bool* v6_flag,
					  Vec<int>& fluidId)
{
  
  // Check if all the nodes of the element are "valid".
  // If not, flag the element and the gradient is set to 0
  //---------------------------------------------------
  int le, N1, N2;
  bool e_flag = false;
  for (int j=0; j<6; j++){

    le = elem.edgeNum(j);

    N1 = elem.nodeNum( elem.edgeEnd(j,0) );
    N2 = elem.nodeNum( elem.edgeEnd(j,1) );
   
    if(fluidId[N1] != fluidId[N2]) 
	 {
    	 e_flag = true;
    	 break;
      }
            
  }
  //---------------------------------------------------

  if(!e_flag){

    *v6_flag = false;

    double dp1dxi[4][3];
    elem.computeGradientP1Function(X, dp1dxi);
  
    int i0 = elem.nodeNum(0);
    int i1 = elem.nodeNum(1);
    int i2 = elem.nodeNum(2);
    int i3 = elem.nodeNum(3);

    for (int k=0; k<dim; ++k) {
      grad[k] = ( (dp1dxi[0][0]*V[i0][k] + dp1dxi[1][0]*V[i1][k] + 
		   dp1dxi[2][0]*V[i2][k] + dp1dxi[3][0]*V[i3][k]) * rij[0] +
		  (dp1dxi[0][1]*V[i0][k] + dp1dxi[1][1]*V[i1][k] + 
		   dp1dxi[2][1]*V[i2][k] + dp1dxi[3][1]*V[i3][k]) * rij[1] +
		  (dp1dxi[0][2]*V[i0][k] + dp1dxi[1][2]*V[i1][k] + 
		   dp1dxi[2][2]*V[i2][k] + dp1dxi[3][2]*V[i3][k]) * rij[2] );
    }

  }else{

    *v6_flag = true;

    for (int k=0; k<dim; ++k) {
      grad[k] = 0.0;
    }

  }

}


//------------------------------------------------------------------------------


template<int dim>
void EdgeGrad<dim>::computeFaceGradient(ElemSet& elems, V6NodeData& data, double rij[3],
					SVec<double,dim>& dVdx, SVec<double,dim>& dVdy,
					SVec<double,dim>& dVdz, double* grad)
{
  Elem& elem = elems[data.tet];
  int n0_loc = elem.faceDef(data.face, 0);
  int n1_loc = elem.faceDef(data.face, 1);
  int n2_loc = elem.faceDef(data.face, 2);
  int n0 = elem.nodeNum(n0_loc);
  int n1 = elem.nodeNum(n1_loc);
  int n2 = elem.nodeNum(n2_loc);

  for (int k=0; k<dim; ++k) {
    grad[k] = ( (dVdx[n2][k] + data.r * (dVdx[n0][k]-dVdx[n2][k]) + 
		 data.t * (dVdx[n1][k]-dVdx[n2][k])) * rij[0] +
		(dVdy[n2][k] + data.r * (dVdy[n0][k]-dVdy[n2][k]) + 
		 data.t * (dVdy[n1][k]-dVdy[n2][k])) * rij[1] +
		(dVdz[n2][k] + data.r * (dVdz[n0][k]-dVdz[n2][k]) + 
		 data.t * (dVdz[n1][k]-dVdz[n2][k])) * rij[2] );
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void EdgeGrad<dim>::computeDerivativeOfFaceGradient(ElemSet& elems, V6NodeData& data, double rij[3], double drij[3],
					SVec<double,dim>& dVdx, SVec<double,dim>& dVdy,
					SVec<double,dim>& dVdz, SVec<double,dim>& ddVdx,
					SVec<double,dim>& ddVdy, SVec<double,dim>& ddVdz,
					double* dGrad)
{

  Elem& elem = elems[data.tet];
  int n0_loc = elem.faceDef(data.face, 0);
  int n1_loc = elem.faceDef(data.face, 1);
  int n2_loc = elem.faceDef(data.face, 2);
  int n0 = elem.nodeNum(n0_loc);
  int n1 = elem.nodeNum(n1_loc);
  int n2 = elem.nodeNum(n2_loc);

  for (int k=0; k<dim; ++k) {
    dGrad[k] = ( (ddVdx[n2][k] + data.r * (ddVdx[n0][k]-ddVdx[n2][k]) +
		 data.t * (ddVdx[n1][k]-ddVdx[n2][k])) * rij[0] +
		(dVdx[n2][k] + data.r * (dVdx[n0][k]-dVdx[n2][k]) +
		 data.t * (dVdx[n1][k]-dVdx[n2][k])) * drij[0] +     
		(ddVdy[n2][k] + data.r * (ddVdy[n0][k]-ddVdy[n2][k]) +
		 data.t * (ddVdy[n1][k]-ddVdy[n2][k])) * rij[1] +
		(dVdy[n2][k] + data.r * (dVdy[n0][k]-dVdy[n2][k]) +
		 data.t * (dVdy[n1][k]-dVdy[n2][k])) * drij[1] +
		(ddVdz[n2][k] + data.r * (ddVdz[n0][k]-ddVdz[n2][k]) +
		 data.t * (ddVdz[n1][k]-ddVdz[n2][k])) * rij[2] +
		(dVdz[n2][k] + data.r * (dVdz[n0][k]-dVdz[n2][k]) +
		 data.t * (dVdz[n1][k]-dVdz[n2][k])) * drij[2]);
  }

}  

//------------------------------------------------------------------------------

template<int dim>
void EdgeGrad<dim>::compute(int l, int i, int j, ElemSet& elems, 
			    SVec<double,3>& X, SVec<double,dim>& V, 
			    SVec<double,dim>& dVdx, SVec<double,dim>& dVdy, 
			    SVec<double,dim>& dVdz, double* ddVij, double* ddVji)
{

  double rij[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};

  if (beta == 0.0) {

    for (int k=0; k<dim; ++k) {
      ddVij[k] = 0.0;
      ddVji[k] = 0.0;
    }

  } else if (v6data[l][0].tet == -1 || v6data[l][1].tet == -1) {
    
    for (int k=0; k<dim; ++k) {
      ddVij[k] = dVdx[i][k]*rij[0] + dVdy[i][k]*rij[1] + dVdz[i][k]*rij[2];
      ddVji[k] = dVdx[j][k]*rij[0] + dVdy[j][k]*rij[1] + dVdz[j][k]*rij[2];
    }

  } else {

    double ddVij_u[dim];
    computeUpwindGradient(elems[v6data[l][0].tet], rij, X, V, ddVij_u);
    double ddVji_u[dim];
    computeUpwindGradient(elems[v6data[l][1].tet], rij, X, V, ddVji_u);
    double ddVij_f[dim];
    computeFaceGradient(elems, v6data[l][0], rij, dVdx, dVdy, dVdz, ddVij_f);
    double ddVji_f[dim];
    computeFaceGradient(elems, v6data[l][1], rij, dVdx, dVdy, dVdz, ddVji_f);

    for (int k=0; k<dim; ++k) {
      double ddVij_c = V[j][k] - V[i][k];
      double ddVi = dVdx[i][k]*rij[0] + dVdy[i][k]*rij[1] + dVdz[i][k]*rij[2];
      double ddVj = dVdx[j][k]*rij[0] + dVdy[j][k]*rij[1] + dVdz[j][k]*rij[2];
      ddVij[k] = 0.5 * (ddVij_c + ddVij_u[k] + xiu * (ddVij_f[k] - 2.0*ddVi + ddVj) +
                                               xic * (ddVij_u[k] - 2.0*ddVij_c + ddVji_u[k]));
      ddVji[k] = 0.5 * (ddVij_c + ddVji_u[k] + xiu * (ddVji_f[k] - 2.0*ddVj + ddVi) +
                                               xic * (ddVji_u[k] - 2.0*ddVij_c + ddVij_u[k]));
    }
  }

  if (tag[i])
    for (int k=0; k<dim; ++k)
      ddVij[k] = 0.0;

  if (tag[j])
    for (int k=0; k<dim; ++k)
      ddVji[k] = 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void EdgeGrad<dim>::computeDerivative(int l, int i, int j, ElemSet& elems,
			    SVec<double,3>& X, SVec<double,3>& dX,
			    SVec<double,dim>& V, SVec<double,dim>& dV,
			    SVec<double,dim>& dVdx, SVec<double,dim>& dVdy,
			    SVec<double,dim>& dVdz, SVec<double,dim>& ddVdx,
			    SVec<double,dim>& ddVdy, SVec<double,dim>& ddVdz,
			    double* dddVij, double* dddVji)          
{

  double rij[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};

  double drij[3] = {dX[j][0] - dX[i][0], dX[j][1] - dX[i][1], dX[j][2] - dX[i][2]};

  if (beta == 0.0) {
    for (int k=0; k<dim; ++k) {
      dddVij[k] = 0.0;
      dddVji[k] = 0.0;
    }
  }
  else if (v6data[l][0].tet == -1 || v6data[l][1].tet == -1) {
    for (int k=0; k<dim; ++k) {
      dddVij[k] = ddVdx[i][k]*rij[0] + dVdx[i][k]*drij[0] + ddVdy[i][k]*rij[1] + dVdy[i][k]*drij[1] + ddVdz[i][k]*rij[2] + dVdz[i][k]*drij[2];
      dddVji[k] = ddVdx[j][k]*rij[0] + dVdx[j][k]*drij[0] + ddVdy[j][k]*rij[1] + dVdy[j][k]*drij[1] + ddVdz[j][k]*rij[2] + dVdz[j][k]*drij[2];
    }
  }
  else {
    double ddVij_u[dim];
    computeUpwindGradient(elems[v6data[l][0].tet], rij, X, V, ddVij_u);

    double dddVij_u[dim];
    computeDerivativeOfUpwindGradient(elems[v6data[l][0].tet], rij, drij, X, dX, V, dV, dddVij_u);

    double ddVji_u[dim];
    computeUpwindGradient(elems[v6data[l][1].tet], rij, X, V, ddVji_u);

    double dddVji_u[dim];
    computeDerivativeOfUpwindGradient(elems[v6data[l][1].tet], rij, drij, X, dX, V, dV, dddVji_u);

    double ddVij_f[dim];
    computeFaceGradient(elems, v6data[l][0], rij, dVdx, dVdy, dVdz, ddVij_f);

    double dddVij_f[dim];
    computeDerivativeOfFaceGradient(elems, v6data[l][0], rij, drij, dVdx, dVdy, dVdz, ddVdx, ddVdy, ddVdz, dddVij_f);

    double ddVji_f[dim];
    computeFaceGradient(elems, v6data[l][1], rij, dVdx, dVdy, dVdz, ddVji_f);

    double dddVji_f[dim];
    computeDerivativeOfFaceGradient(elems, v6data[l][1], rij, drij, dVdx, dVdy, dVdz, ddVdx, ddVdy, ddVdz, dddVji_f);

    for (int k=0; k<dim; ++k) {
      double dddVij_c = dV[j][k] - dV[i][k];
      double dddVi = ddVdx[i][k]*rij[0] + dVdx[i][k]*drij[0] + ddVdy[i][k]*rij[1] + dVdy[i][k]*drij[1] + ddVdz[i][k]*rij[2] + dVdz[i][k]*drij[2];
      double dddVj = ddVdx[j][k]*rij[0] + dVdx[j][k]*drij[0] + ddVdy[j][k]*rij[1] + dVdy[j][k]*drij[1] + ddVdz[j][k]*rij[2] + dVdz[j][k]*drij[2];
      dddVij[k] = 0.5 * (dddVij_c + dddVij_u[k] + xiu * (dddVij_f[k] - 2.0*dddVi + dddVj) +
			xic * (dddVij_u[k] - 2.0*dddVij_c + dddVji_u[k]));
      dddVji[k] = 0.5 * (dddVij_c + dddVji_u[k] + xiu * (dddVji_f[k] - 2.0*dddVj + dddVi) +
			xic * (dddVji_u[k] - 2.0*dddVij_c + dddVij_u[k]));
    }
  }

  if (tag[i])
    for (int k=0; k<dim; ++k)
      dddVij[k] = 0.0;

  if (tag[j])
    for (int k=0; k<dim; ++k)
      dddVji[k] = 0.0;

}

//------------------------------------------------------------------------------

//d2d
template<int dim>
void EdgeGrad<dim>::compute(int l, int i, int j, ElemSet& elems, 
			    SVec<double,3>& X, SVec<double,dim>& V, 
			    SVec<double,dim>& dVdx, 
			    SVec<double,dim>& dVdy, 
			    SVec<double,dim>& dVdz, 
			    Vec<int> &fluidId, double* ddVij, double* ddVji,
			    LevelSetStructure &LSS)
{

	 bool isValid = true;

	 if(fluidId[i] != fluidId[j]) isValid = false;
	 if(LSS.edgeIntersectsStructure(0.0, l)) isValid = false;
	 if(!LSS.isActive(0.0, i) || !LSS.isActive(0.0, j)) isValid = false;
	 
	 if(!isValid)
	 {
		 for(int k=0; k<dim; ++k)
		 {
      ddVij[k] = 0.0;
      ddVji[k] = 0.0;
    }
    return;
  }

  // ~~~~~~

  double rij[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
  
  if (beta == 0.0) {

    for (int k=0; k<dim; ++k) {
      ddVij[k] = 0.0;
      ddVji[k] = 0.0;
    }

  }else if (v6data[l][0].tet == -1 || v6data[l][1].tet == -1) {

    for (int k=0; k<dim; ++k) {
      ddVij[k] = dVdx[i][k]*rij[0] + dVdy[i][k]*rij[1] + dVdz[i][k]*rij[2];
      ddVji[k] = dVdx[j][k]*rij[0] + dVdy[j][k]*rij[1] + dVdz[j][k]*rij[2];
    }

  } else {

    bool v6_flag1, v6_flag2;

    double ddVij_u[dim];
    computeUpwindGradient(elems[v6data[l][0].tet], rij, X, V, ddVij_u, &v6_flag1, fluidId, LSS);
    double ddVji_u[dim];
    computeUpwindGradient(elems[v6data[l][1].tet], rij, X, V, ddVji_u, &v6_flag2, fluidId, LSS);

    if(v6_flag1 || v6_flag2) {
      for (int k=0; k<dim; ++k){
	ddVij[k] = dVdx[i][k]*rij[0] + dVdy[i][k]*rij[1] + dVdz[i][k]*rij[2];
	ddVji[k] = dVdx[j][k]*rij[0] + dVdy[j][k]*rij[1] + dVdz[j][k]*rij[2];
      }

    }else{

      double ddVij_f[dim];
      computeFaceGradient(elems, v6data[l][0], rij, dVdx, dVdy, dVdz, ddVij_f);
      double ddVji_f[dim];
      computeFaceGradient(elems, v6data[l][1], rij, dVdx, dVdy, dVdz, ddVji_f);

      for (int k=0; k<dim; ++k) {

	double ddVij_c = V[j][k] - V[i][k];

	double ddVi = dVdx[i][k]*rij[0] + dVdy[i][k]*rij[1] + dVdz[i][k]*rij[2];
	double ddVj = dVdx[j][k]*rij[0] + dVdy[j][k]*rij[1] + dVdz[j][k]*rij[2];

	ddVij[k] = 0.5 * (ddVij_c + ddVij_u[k] + xiu * (ddVij_f[k] - 2.0*ddVi + ddVj) +
                                                 xic * (ddVij_u[k] - 2.0*ddVij_c + ddVji_u[k]));
	ddVji[k] = 0.5 * (ddVij_c + ddVji_u[k] + xiu * (ddVji_f[k] - 2.0*ddVj + ddVi) +
                                                 xic * (ddVji_u[k] - 2.0*ddVij_c + ddVij_u[k]));
      }

    }

  }

  if (tag[i])
    for (int k=0; k<dim; ++k)
      ddVij[k] = 0.0;

  if (tag[j])
    for (int k=0; k<dim; ++k)
      ddVji[k] = 0.0;

}


//------------------------------------------------------------------------------

//d2d
template<int dim>
void EdgeGrad<dim>::compute(int l, int i, int j, ElemSet& elems, 
			    SVec<double,3>& X, SVec<double,dim>& V, 
			    SVec<double,dim>& dVdx, 
			    SVec<double,dim>& dVdy, 
			    SVec<double,dim>& dVdz, 
			    Vec<int> &fluidId, double* ddVij, double* ddVji)
{

  if((fluidId[i] != fluidId[j])) {

    for (int k=0; k<dim; ++k) {
      ddVij[k] = 0.0;
      ddVji[k] = 0.0;
    }
    return;

  }

  // ~~~~~~

  double rij[3] = {X[j][0] - X[i][0], X[j][1] - X[i][1], X[j][2] - X[i][2]};
  
  if (beta == 0.0) {

    for (int k=0; k<dim; ++k) {
      ddVij[k] = 0.0;
      ddVji[k] = 0.0;
    }

  }else if (v6data[l][0].tet == -1 || v6data[l][1].tet == -1) {

    for (int k=0; k<dim; ++k) {
      ddVij[k] = dVdx[i][k]*rij[0] + dVdy[i][k]*rij[1] + dVdz[i][k]*rij[2];
      ddVji[k] = dVdx[j][k]*rij[0] + dVdy[j][k]*rij[1] + dVdz[j][k]*rij[2];
    }

  } else {

    bool v6_flag1, v6_flag2;

    double ddVij_u[dim];
    computeUpwindGradient(elems[v6data[l][0].tet], rij, X, V, ddVij_u, &v6_flag1, fluidId);
    double ddVji_u[dim];
    computeUpwindGradient(elems[v6data[l][1].tet], rij, X, V, ddVji_u, &v6_flag2, fluidId);

    if(v6_flag1 || v6_flag2) {
      for (int k=0; k<dim; ++k){
	ddVij[k] = dVdx[i][k]*rij[0] + dVdy[i][k]*rij[1] + dVdz[i][k]*rij[2];
	ddVji[k] = dVdx[j][k]*rij[0] + dVdy[j][k]*rij[1] + dVdz[j][k]*rij[2];
      }

    }else{

      double ddVij_f[dim];
      computeFaceGradient(elems, v6data[l][0], rij, dVdx, dVdy, dVdz, ddVij_f);
      double ddVji_f[dim];
      computeFaceGradient(elems, v6data[l][1], rij, dVdx, dVdy, dVdz, ddVji_f);

      for (int k=0; k<dim; ++k) {

	double ddVij_c = V[j][k] - V[i][k];

	double ddVi = dVdx[i][k]*rij[0] + dVdy[i][k]*rij[1] + dVdz[i][k]*rij[2];
	double ddVj = dVdx[j][k]*rij[0] + dVdy[j][k]*rij[1] + dVdz[j][k]*rij[2];

	ddVij[k] = 0.5 * (ddVij_c + ddVij_u[k] + xiu * (ddVij_f[k] - 2.0*ddVi + ddVj) +
                                                 xic * (ddVij_u[k] - 2.0*ddVij_c + ddVji_u[k]));
	ddVji[k] = 0.5 * (ddVij_c + ddVji_u[k] + xiu * (ddVji_f[k] - 2.0*ddVj + ddVi) +
                                                 xic * (ddVji_u[k] - 2.0*ddVij_c + ddVij_u[k]));
      }

    }

  }

  if (tag[i])
    for (int k=0; k<dim; ++k)
      ddVij[k] = 0.0;

  if (tag[j])
    for (int k=0; k<dim; ++k)
      ddVji[k] = 0.0;

}

