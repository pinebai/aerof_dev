#include <InletNode.h>

#include <BcData.h>
#include <Node.h>
#include <Elem.h>
#include <Vector.h>
#include <Vector3D.h>
#include <Extrapolation.h>
#include <GeoState.h>

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>

//------------------------------------------------------------------------------

template <int dim>
void InletNode::assignFreeStreamValues(int type, double *Uin, double *Uout, double *U)
{

  int k;

  if (type == BC_INLET_MOVING || type == BC_INLET_FIXED)
    for (k=0; k<dim; ++k)
      U[k] = Uin[k];
  else if (type == BC_OUTLET_MOVING || type == BC_OUTLET_FIXED)
    for (k=0; k<dim; ++k)
      U[k] = Uout[k];
  else
    fprintf(stderr, "the type of the inlet node should not be anything else than BC_INLET or BC_OUTLET");

}

//------------------------------------------------------------------------------

template<int dim>
void InletNode::computeZeroExtrapolation(VarFcn* vf, bool flag, Vec3D& normal,
				double* Ub, double* Ufar, double* Vinter1, double* Vinter2, 
				int* tet, int* mask, SVec<double,dim> &Ubc,
				SVec<double,3> &X, int i,
				int *locToGlobNodeMap)
// function called only for single-phase flow //
{

  double S = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
  double ooS = 1.0/S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};

  double Vfar[dim];
  double Vinter[dim];	//interpolated values from inside the domain
  double Vextra[dim];	//values that will be imposed
  double dV[dim];

  vf->conservativeToPrimitive(Ufar,Vfar);
  double cb = vf->computeSoundSpeed(Vfar);
  double unb = Vfar[1]*n[0] + Vfar[2]*n[1] + Vfar[3]*n[2];

  //for linear extrapolation:
  // if two inlet nodes get their extrapolated values from
  // one same interior node, then one value is extrapolated linearly
  // and the second is extrapolated constantly
  chooseExtrapolation(dim, unb, tet, mask, Vinter1, Vinter2, Vinter);
  for (int idim=0; idim<dim; idim++)
    dV[idim] = Vinter[idim]-Vfar[idim];
  master = flag;

  if(!flag){
    if(!locToGlobNodeMap)
      for(int j = 0;  j<dim; j++)
        Ubc[i][j]=0.0;
    else
      for(int j = 0;  j<dim; j++)
        Ubc[node][j]=0.0;

  }else{
    vf->extrapolatePrimitive(unb, cb, Vfar, Vinter, Vextra);
    //vf->extrapolateBoundaryCharacteristic(n,unb,cb,Vfar,dV);
    //for (int idim=0; idim<dim; idim++){
    //  Vextra[idim] = Vfar[idim]+dV[idim];
    //}
    if(Vextra[0]<=0.0 || Vextra[4]<=0.0){
      fprintf(stdout, "*** Error: negative density or pressure for inlet nodes\n");
      exit(1);
    }
    if(!locToGlobNodeMap){
      vf->primitiveToConservative(Vextra, Ubc[i]);
    }else{
       vf->primitiveToConservative(Vextra, Ubc[node]);
    }
  }
}

//------------------------------------------------------------------------------
template<int dim>
void InletNode::computeZeroExtrapolation(VarFcn* vf, bool flag, Vec3D& normal,
                                double* Ub, double* Ufar, int fluidId, double* Vinter1,
				double* Vinter2, int* tet, int* mask,
				SVec<double,dim> &Ubc, SVec<double,3> &X, int i,
                                int *locToGlobNodeMap)
{
  double S = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
  double ooS = 1.0/S;
  double n[3] = {normal[0]*ooS, normal[1]*ooS, normal[2]*ooS};

  double Vfar[dim];
  double Vinter[dim];
  double Vextra[dim];
  double dV[dim];

  vf->conservativeToPrimitive(Ufar,Vfar, fluidId);
  double cb = vf->computeSoundSpeed(Vfar, fluidId);
  double unb = Vfar[1]*n[0] + Vfar[2]*n[1] + Vfar[3]*n[2];

  chooseExtrapolation(dim, unb, tet, mask, Vinter1, Vinter2, Vinter);
  for (int idim=0; idim<dim; idim++)
    dV[idim] = Vinter[idim]-Vfar[idim];
  master = flag;

  if(!flag){
    if(!locToGlobNodeMap)
      for(int j = 0;  j<dim; j++)
        Ubc[i][j]=0.0;
    else
      for(int j = 0;  j<dim; j++)
        Ubc[node][j]=0.0;
  }else{
    vf->extrapolatePrimitive(unb, cb, Vfar, Vinter, Vextra, fluidId);
    //vf->extrapolateBoundaryCharacteristic(n,unb,cb,Vfar,dV,fluidId);
    //for (int idim=0; idim<dim; idim++)
    //  Vextra[idim] = Vfar[idim]+dV[idim];
    if(!locToGlobNodeMap){
      vf->primitiveToConservative(Vextra, Ubc[i], fluidId);
    }else{
       vf->primitiveToConservative(Vextra, Ubc[node], fluidId);
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
void InletNode::computeDifference(VarFcn* vf, SVec<double, dim> &V, SVec<double,dim> &rhs)
{

  if (master){
    double Utemp[dim];
    vf->primitiveToConservative(V[node], Utemp);
    for(int kk=0; kk<5; kk++)
      rhs[node][kk] -= Utemp[kk];
  }
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
template<int dim>
void InletNodeSet::recomputeRHS(VarFcn* vf, Extrapolation<dim>* xpol, ElemSet &elems, SVec<double,dim>&V,
  				BcData<dim>& bcData, GeoState& geoState, SVec<double,dim>& rhs,
				SVec<double,3>& X, int* locToGlobNodeMap)
{
//extrapolation routine for implicit time stepping
  double Vinter1[dim];
  double Vinter2[dim];
  bool master ;
  int node;
  int maskLinear[V.size()];  //mask of interior nodes used for linear extrapolation
                         // linear extrapolation uses two interior points
                         // the most interior point is an interpolation of three defined nodes(nodes from the mesh)
                         // at the beginning of the procedure, all the nodes are tagged to 0
                         // as the procedure goes along, the tags can stay 0, become 1, 2 or 3
                         //  0 - unused
                         //  1 - inlet extrapolation use
                         //  2 - outlet extrapolation use
                         //  3 - nor inlet nor outlet extrapolation use (u.n = 0)
                         // if a node is needed for the linear extrapolation:
                         //       -if its tag is 0, it can be used
                         //       -if its tag does not correspond to the way it will be used for this node, it cannot be used
                         //       -if its tag corresponds to the way it will be used for this node, it will be used
  for (int i=0; i<V.size(); i++)
    maskLinear[i] = 0;

  int LinTet[3];
  Vec<Vec3D> &n = geoState.getInletNodeNormal();
  SVec<double,dim> &Ub = bcData.getInletNodeStateVector();
  SVec<double,dim> &Ufar = bcData.getInletBoundaryVector();
  for ( int i=0; i<numInletNodes; i++){
    node = inletNodes[i].getNodeNum();
    xpol->computeFaceInterpolation(i, master, node, elems, Ufar, V, Vinter1, Vinter2, LinTet, locToGlobNodeMap, X);
    inletNodes[i].computeZeroExtrapolation(vf, master, n[i], Ub[i], Ufar[node], Vinter1, Vinter2, LinTet, maskLinear, rhs, X, i, locToGlobNodeMap);
    inletNodes[i].computeDifference(vf, V, rhs);
  }
}  

//------------------------------------------------------------------------------
template<int dim>
void InletNodeSet::recomputeRHS(VarFcn* vf, Extrapolation<dim>* xpol, ElemSet &elems, SVec<double,dim>&V,
                                Vec<int> &fluidId, BcData<dim>& bcData, GeoState& geoState,
                                SVec<double,dim>& rhs, SVec<double,3>& X, int* locToGlobNodeMap)
{
// This routine is for LevelSet method
  fprintf(stderr, "extrapolation recomputeRHS for implicitLS\n");
  double Vinter1[dim];
  double Vinter2[dim];
  bool master ;
  int node;
  int maskLinear[V.size()]; 

  for (int i=0; i<V.size(); i++)
    maskLinear[i] = 0;

  int LinTet[3];
  Vec<Vec3D> &n = geoState.getInletNodeNormal();
  SVec<double,dim> &Ub = bcData.getInletNodeStateVector();
  SVec<double,dim> &Ufar = bcData.getInletBoundaryVector();

  for ( int i=0; i<numInletNodes; i++){
    node = inletNodes[i].getNodeNum();
    xpol->computeFaceInterpolation(i, master, node, elems, Ufar, V, Vinter1, Vinter2, LinTet, locToGlobNodeMap, X);
    inletNodes[i].computeZeroExtrapolation(vf, master, n[i], Ub[i], Ufar[node], fluidId[node], Vinter1, 
					   Vinter2, LinTet, maskLinear, rhs, X, i, locToGlobNodeMap);
    inletNodes[i].computeDifference(vf, V, rhs);
  }

}
	
//------------------------------------------------------------------------------
template<int dim>
void InletNodeSet::recomputeResidual(SVec<double,dim> &F,SVec<double,dim> &Finlet)
{

  int node;
  for (int i = 0; i<numInletNodes; i++){
      node = inletNodes[i].getNodeNum();
      for (int j = 0; j<dim; j++)
        Finlet[node][j] = F[node][j];
  }

}
//------------------------------------------------------------------------------
template<int dim>
void InletNodeSet::getExtrapolationValue(Extrapolation<dim>* xpol,SVec<double,dim> &V, SVec<double,dim> &Ubc,
			 VarFcn* vf, BcData<dim>& bcData, GeoState& geoState, ElemSet &elems, int* locToGlobNodeMap, 
			 SVec<double,3>& X)
{
//extrapolation routine for explicit time stepping	
  Vec<Vec3D> &n = geoState.getInletNodeNormal();
  SVec<double,dim> &Ub = bcData.getInletNodeStateVector();
  SVec<double,dim> &Ufar = bcData.getInletBoundaryVector();

  double Vinter1[dim];
  double Vinter2[dim];
  bool master;
  int node;
  int LinTet[3];
  int maskLinear[V.size()];  //mask of interior nodes used for linear extrapolation
                         // linear extrapolation uses two interior points
                         // the most interior point is an interpolation of three defined nodes(nodes from the mesh)
                         // at the beginning of the procedure, all the nodes are tagged to 0
                         // as the procedure goes along, the tags can stay 0, become 1, 2 or 3
                         //  0 - unused
                         //  1 - inlet extrapolation use
                         //  2 - outlet extrapolation use
                         //  3 - nor inlet nor outlet extrapolation use (u.n = 0)
                         // if a node is needed for the linear extrapolation:
                         //       -if its tag is 0, it can be used
                         //       -if its tag does not correspond to the way it will be used for this node, it cannot be used
                         //       -if its tag corresponds to the way it will be used for this node, it will be used
  for (int i=0; i<V.size(); i++)
    maskLinear[i] = 0;

  for (int i=0; i<numInletNodes; i++){
    node = inletNodes[i].getNodeNum();
    xpol->computeFaceInterpolation(i, master, node, elems, Ufar, V, Vinter1, Vinter2, LinTet, locToGlobNodeMap, X);
    inletNodes[i].computeZeroExtrapolation(vf, master, n[i], Ub[i], Ufar[node], Vinter1, Vinter2, LinTet, maskLinear, Ubc, X, i);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void InletNodeSet::applyExtrapolationToSolutionVector(Extrapolation<dim>* xpol, SVec<double,dim> &U, SVec<double,dim> &Ubc, int* locToGlobNodeMap)
{

  for (int i=0; i<numInletNodes; i++){
    for (int k=0; k<dim; k++){
      U[inletNodes[i].getNodeNum()][k] = Ubc[i][k];
    }
  }
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template<int dim>
void InletNodeSet::printVariable(SVec<double,dim> &V, Connectivity* sharedInletNodes, VarFcn *vf)
{
	double U[5];
	int node;
	fprintf(stderr, "printing variables: \n");
	for (int i=0; i<numInletNodes; i++){
            node = inletNodes[i].getNodeNum();
	    vf->conservativeToPrimitive(V[node], U);// assumption : only one phase at infinity boundary
	    fprintf(stderr, "V[%d] = ", node);                                   
	    for (int k=0; k<dim; k++){
		fprintf(stderr, " %f  ", U[k]);
	    }
            fprintf(stderr,"\n");
	}

}

//-----------------------------------------------------------------------------

template<int dim>
void InletNodeSet::printInletVariable(SVec<double,dim> &V, Connectivity* sharedInletNodes)
{
                                                                                                                                                                                                     
        int node;
        fprintf(stderr, "printing inlet variables: \n");
        for (int i=0; i<numInletNodes; i++){
	    node = inletNodes[i].getNodeNum();
            fprintf(stderr, "V[%d] = ", node);
            for (int k=0; k<dim; k++){
                fprintf(stderr, " %f  ", V[i][k]);
            }
            fprintf(stderr,"\n");
        }
                                                                                                                                                                                                     
}
                                                                                                                                                                                                     
//-----------------------------------------------------------------------------
                                                                                                                                                                                                     
template<int dim>
void InletNodeSet::checkExtrapolationValue(SVec<double,dim> &U, Connectivity* sharedInletNodes, int* nodeType, VarFcn* vf, BcData<dim>& bcData, GeoState& geoState)
{

	Vec<Vec3D> &n = geoState.getInletNodeNormal();
        SVec<double,dim> &Ub = bcData.getInletNodeStateVector();	

	double Vb[dim];
        double machb ;
        double cb;
	double unb;

	fprintf(stderr, "checkExtrapolationValue constoprim\n");

	int node, i, k;
	for (i=0; i<numInletNodes; i++){
		node = inletNodes[i].getNodeNum();
 		vf->conservativeToPrimitive(Ub[i],Vb); // assumption : only one phase at infinity boundary
		machb = vf->computeMachNumber(Vb);
	        cb = vf->computeSoundSpeed(Vb);
		unb = Vb[1]*n[i][0] + Vb[2]*n[i][1] + Vb[3]*n[i][2];
		if (U[i][0]!=1.0 && unb+cb>0 && unb<0){
			fprintf(stderr, "the density should be imposed by the infinity points and should not be %f for node %d\n", U[i][0], node);
		}
		if (U[i][4]!=Ub[i][4] && -unb+cb>0 && unb>0){
			fprintf(stderr, "the energy should be imposed by the infinity points and should not be %f for node %d\n", U[i][4], node);
		}
	}
}

//-----------------------------------------------------------------------------





