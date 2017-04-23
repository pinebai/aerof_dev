#include <Elem.h>
#include <Face.h>

#include <FemEquationTerm.h>
#include <MacroCell.h>
#include <VMSLESTerm.h>
#include <DynamicVMSTerm.h>
#include <SmagorinskyLESTerm.h>
#include <WaleLESTerm.h>
#include <DynamicLESTerm.h>
#include <GenMatrix.h>
#include <cmath>
#include <GeoState.h>
#include <limits>
#include <VarFcn.h>
#include "LevelSet/LevelSetStructure.h"

#include <cstring> // for memset

//------------------------------------------------------------------------------
//--------------functions in ElemSet class
//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeTimeStep(
                FemEquationTerm *fet,
                SVec<double,3> &X,
                SVec<double,dim> &V,
                Vec<double> &idtv,
                LevelSetStructure *LSS)
{
  int *nodeNumber,numberOfNodes;
  double Vmid[dim],oneOnDofs;
  double Xmid[3];
  double localDt;
  double nGrad[4][3];
  int    nodeID;
  for (int i=0; i<numElems; ++i)
    {
      nodeNumber    = elems[i]->nodeNum();
      numberOfNodes = elems[i]->numNodes();

      oneOnDofs = 1.0/((double) numberOfNodes);

	  bool isValid = true;

	  if(LSS)
	  {		  
		  int numberOfEdges = elems[i]->numEdges();
		  
		  for(int j=0; j<numberOfEdges; ++j)
		  {
			  int l = elems[i]->edgeNum(j);

			  if(LSS->edgeIntersectsStructure(0, l)) isValid = false;

			  int e1 = nodeNumber[elems[i]->edgeEnd(j, 0)];
			  int e2 = nodeNumber[elems[i]->edgeEnd(j, 1)];

			  if(!LSS->isActive(0, e1) || !LSS->isActive(0, e2)) isValid = false;
		  }
	  }

	  if(!isValid) continue;

      for (int k=0; k<dim; ++k) 
	{
	  Vmid[k] = 0.0;
	  for(int j=0;j<numberOfNodes;++j)
	    { 
	      nodeID   = nodeNumber[j];
	      Vmid[k] += V[nodeID][k];	  
	    }
	  Vmid[k] *= oneOnDofs;
	}
      // localDt_i = mu/rho*||n_i||^2/|T|
      // This line gives: 9.0*mu/rho*|T|
      localDt = 9.0*fet->computeViscousTimeStep(Xmid, Vmid)*elems[i]->computeVolume(X);
      // gradPhi = n_i/(3.0*|T|)
      elems[i]->computeGradientP1Function(X,nGrad);
      for(int j=0;j<numberOfNodes;++j)
	{
	  nodeID = nodeNumber[j];

	  idtv[nodeID]  += localDt*(nGrad[j][0]*nGrad[j][0] + nGrad[j][1]*nGrad[j][1] + nGrad[j][2]*nGrad[j][2]);
	}
    }
}
//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeGalerkinTerm(FemEquationTerm *fet, GeoState &geoState, 
                SVec<double,3> &X, SVec<double,dim> &V,
                SVec<double,dim> &R,
                Vec<GhostPoint<dim>*> *ghostPoints,
                LevelSetStructure *LSS, bool externalSI)
{

  Vec<double> &d2wall = geoState.getDistanceToWall();

  if(!externalSI)
  {
    if (sampleMesh)
    {
      for (int iElem=0; iElem<numSampledElems; ++iElem)
      elems[ (elemsConnectedToSampleNode[iElem]) ]->computeGalerkinTerm(fet, X, d2wall, V, R, ghostPoints,LSS);
    }
    else
    {
      for (int iElem=0; iElem<numSampledElems; ++iElem)//TODO print the R norm or each element and see how it develops
      elems[ iElem ]->computeGalerkinTerm(fet, X, d2wall, V, R, ghostPoints,LSS);
    }
  }
  else
  {
    if (sampleMesh)
    {
      for (int iElem=0; iElem<numSampledElems; ++iElem)
        elems[ (elemsConnectedToSampleNode[iElem]) ]->computeGalerkinTerm_e(fet, X, d2wall, V, R, ghostPoints, LSS);
    }
    else
    {
      for (int iElem=0; iElem<numSampledElems; ++iElem)
        elems[ iElem ]->computeGalerkinTerm_e(fet, X, d2wall, V, R, ghostPoints, LSS);
    }
  }

}

//------------------------------------------------------------------------------

// Not used? KMW

template<int dim>
void ElemSet::computeGalerkinTermRestrict(FemEquationTerm *fet, GeoState &geoState, 
				  SVec<double,3> &X, SVec<double,dim> &V, 
				  SVec<double,dim> &R,const std::vector<int> &sampledLocElem,
				  Vec<GhostPoint<dim>*> *ghostPoints,LevelSetStructure *LSS)
{

  Vec<double> &d2wall = geoState.getDistanceToWall();

	int i;
	for (int iElem=0; iElem<numSampledElems; ++iElem) {
		i = sampleMesh ? elemsConnectedToSampleNode[iElem]: iElem;
    elems[i]->computeGalerkinTerm(fet, X, d2wall, V, R, ghostPoints,LSS);
	}

}
//------------------------------------------------------------------------------

/****************************************************************************************
 * Computes the derivative of the viscous term for non-embedded simulations.            *
 * This is the non-sparse implementation                                           (MB) *
 ****************************************************************************************/
template<int dim>
void ElemSet::computeDerivativeOfGalerkinTerm(FemEquationTerm *fet,
                GeoState &geoState,
                SVec<double,3> &X,
                SVec<double,3> &dX,
                SVec<double,dim> &V, SVec<double,dim> &dV,
                double dMach,
                SVec<double,dim> &dR)
{
  Vec<double> &d2wall = geoState.getDistanceToWall();

  for (int i=0; i<numElems; ++i)
    elems[i]->computeDerivativeOfGalerkinTerm(fet, X, dX, d2wall, V, dV, dMach, dR);
}




//TODO VISCOUSDERIV
/****************************************************************************************
 * Computes the derivative of the viscous term for non-embedded simulations.            *
 * This is the non-sparse implementation                                           (MB) *
 ****************************************************************************************/
template<int dim>
void ElemSet::computeDerivativeOfGalerkinTermEmb(
                FemEquationTerm *fet,
                GeoState &geoState,
                SVec<double,3> &X,
                SVec<double,3> &dX,
                SVec<double,dim> &V, SVec<double,dim> &dV,
                double dMach,
                SVec<double,dim> &dR,
                Vec<GhostPoint<dim>*> *ghostPoints,LevelSetStructure *LSS)
{
  Vec<double> &d2wall = geoState.getDistanceToWall();

  for (int i=0; i<numElems; ++i)
    elems[i]->computeDerivativeOfGalerkinTermEmb(fet, X, dX, d2wall, V, dV, dMach, dR, ghostPoints,LSS);
}




//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeDerivativeOperatorsOfGalerkinTerm(FemEquationTerm *fet, GeoState &geoState,
				 SVec<double,3> &X, SVec<double,dim> &V, RectangularSparseMat<double,3,dim> &dViscousFluxdX)
{

  Vec<double> &d2wall = geoState.getDistanceToWall();

  for (int i=0; i<numElems; ++i)
    elems[i]->computeDerivativeOperatorsOfGalerkinTerm(fet, X, d2wall, V, dViscousFluxdX);
}

//------------------------------------------------------------------------------
template<int dim>
void ElemSet::computeMBarAndM(DynamicVMSTerm *dvmst,
			      SVec<double,dim> **VBar,
			      SVec<double,1> **volRatio,
			      SVec<double,3> &X,
			      SVec<double,dim> &V,
			      SVec<double,dim> &MBar,
			      SVec<double,dim> &M)
{

  for (int i=0; i<numElems; ++i)
   elems[i]->computeMBarAndM(dvmst, VBar, volRatio, X, V, MBar, M);

}

//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeDynamicVMSTerm(DynamicVMSTerm *dvmst,
				    SVec<double,dim> **VBar,
				    SVec<double,3> &X,
				    SVec<double,dim> &V, SVec<double,dim> &S,
				    Vec<double> &CsDelSq, Vec<double> &PrT,
				    Vec<double> *Cs, Vec<double> &Delta)
{

  for (int i=0; i<numElems; ++i)
    elems[i]->computeDynamicVMSTerm(dvmst, VBar, X, V, S, CsDelSq, PrT, Cs, Delta);

}

//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeVMSLESTerm(VMSLESTerm *vmst,
				SVec<double,dim> &VBar,
				SVec<double,3> &X,
				SVec<double,dim> &V,
				SVec<double,dim> &Sigma)
                                                                                                                          
{
                                                                                                                          
  for (int i=0; i<numElems; ++i)
    elems[i]->computeVMSLESTerm(vmst, VBar, X, V, Sigma);
                                                                                                                          
}

//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeSmagorinskyLESTerm(SmagorinskyLESTerm *smag, SVec<double,3> &X,
					SVec<double,dim> &V, SVec<double,dim> &R,
													 Vec<GhostPoint<dim>*> *ghostPoints, 
													 LevelSetStructure *LSS, bool externalSI)

{

	if(!externalSI)
	{
  for (int i=0; i<numElems; ++i)
    elems[i]->computeSmagorinskyLESTerm(smag, X, V, R, ghostPoints, LSS);
}
	else
	{
		for(int i=0; i<numElems; ++i)	  
			elems[i]->computeSmagorinskyLESTerm_e(smag, X, V, R, ghostPoints, LSS);
	}

}

//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeWaleLESTerm(WaleLESTerm *wale, SVec<double,3> &X,
				SVec<double,dim> &V, SVec<double,dim> &R,
											Vec<GhostPoint<dim>*> *ghostPoints,
											LevelSetStructure *LSS, bool externalSI)

{

	if(!externalSI)
	{
  for (int i=0; i<numElems; ++i)
    elems[i]->computeWaleLESTerm(wale, X, V, R, ghostPoints, LSS);
}
	else
	{
		for(int i=0; i<numElems; ++i)
			elems[i]->computeWaleLESTerm_e(wale, X, V, R, ghostPoints, LSS);
	}

}

//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeDynamicLESTerm(DynamicLESTerm *dles, SVec<double,2> &Cs, 
                                    SVec<double,3> &X, SVec<double,dim> &V, 
												SVec<double,dim> &R,	Vec<GhostPoint<dim>*> *ghostPoints, 
												LevelSetStructure *LSS, bool externalSI)

{
	if(!externalSI)
	{
 for (int i=0; i<numElems; ++i)
    elems[i]->computeDynamicLESTerm(dles, Cs, X, V, R, ghostPoints, LSS);
	}
	else
	{
		for (int i=0; i<numElems; ++i)
			elems[i]->computeDynamicLESTerm_e(dles, Cs, X, V, R, ghostPoints, LSS);
	}
}

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
void ElemSet::computeJacobianGalerkinTerm(FemEquationTerm *fet, GeoState &geoState, 
					  SVec<double,3> &X, Vec<double> &ctrlVol,
					  SVec<double,dim> &V, GenMat<Scalar,neq> &A,
                                          Vec<GhostPoint<dim>*>* ghostPoints, 
														LevelSetStructure *LSS, bool externalSI)
{

  Vec<double> &d2wall = geoState.getDistanceToWall();

  if(!externalSI)
  {
  for (int i=0; i<numElems; ++i)
    elems[i]->computeJacobianGalerkinTerm(fet, X, ctrlVol, d2wall, V, A, ghostPoints,LSS);
  }
  else
  {
	  for (int i=0; i<numElems; ++i)
	  	  elems[i]->computeJacobianGalerkinTerm_e(fet, X, ctrlVol, d2wall, V, A, ghostPoints, LSS);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void ElemSet::computeTestFilterAvgs(SVec<double,dim> &VCap, SVec<double,16> &Mom_Test,
                                   SVec<double,6> &Sij_Test, Vec<double> &modS_Test, 
                                   SVec<double,8> &Eng_Test, SVec<double,3> &X, SVec<double,dim> &V, 
                                   double gam, double R,
												Vec<GhostPoint<dim>*>* ghostPoints, 
												LevelSetStructure *LSS, bool externalSI)
{

	if(!externalSI)
	{
 for (int i=0; i<numElems; ++i)
   elems[i]->computeP1Avg(VCap, Mom_Test, Sij_Test, modS_Test, Eng_Test, X, V, gam, R, ghostPoints, LSS);
	}
	else
	{
		for(int i=0; i<numElems; ++i)
			elems[i]->computeP1Avg_e(VCap, Mom_Test, Sij_Test, modS_Test, Eng_Test, X, V, gam, R, ghostPoints, LSS);
	}

}

//------------------------------------------------------------------------------
// Level Set Reinitialization functions

template<int dimLS>
void ElemSet::computeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                       SVec<double,dimLS> &ddx, SVec<double,dimLS> &ddy,
                                       SVec<double,dimLS> &ddz,
                                       SVec<double,dimLS> &Phi,SVec<double,1> &Psi)
{

  for (int i=0; i<numElems; i++)
    elems[i]->computeDistanceCloseNodes(lsdim,Tag,X,ddx,ddy,ddz,Phi,Psi);

}
//------------------------------------------------------------------------------
template<int dimLS>
void ElemSet::recomputeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                       SVec<double,dimLS> &ddx, SVec<double,dimLS> &ddy,
                                       SVec<double,dimLS> &ddz,
                                       SVec<double,dimLS> &Phi,SVec<double,1> &Psi)
{

  for (int i=0; i<numElems; i++)
    elems[i]->recomputeDistanceCloseNodes(lsdim,Tag,X,ddx,ddy,ddz,Phi,Psi);

}
//------------------------------------------------------------------------------
template<int dimLS>
void ElemSet::computeDistanceLevelNodes(int lsdim, Vec<int> &Tag, int level,
                                       SVec<double,3> &X, SVec<double,1> &Psi, SVec<double,dimLS> &Phi)
{

  for (int i=0; i<numElems; i++)
    elems[i]->computeDistanceLevelNodes(lsdim,Tag,level,X,Psi,Phi);

}
// End of Level Set Reinitialization functions
//-------------------------------------------------------------------------------

template<int dim, class Obj>
void ElemSet::integrateFunction(Obj* obj,SVec<double,3> &X,SVec<double,dim>& V, void (Obj::*F)(int node, const double* loc,double* f),int npt) {

  for (int i=0; i<numElems; ++i)
    elems[i]->integrateFunction(obj, X,V,F,npt);
}

template<int dim> 
void ElemSet::interpolateSolution(SVec<double,3>& X, SVec<double,dim>& U, 
                                  const std::vector<Vec3D>& locs, double (*sol)[dim],
                                  int* status, int* last, LevelSetStructure* LSS,
                                  Vec<GhostPoint<dim>*>* ghostPoints, VarFcn* varFcn,
				  bool assumeCache) {
  
  int nn;
  Vec3D bbox[2];
  std::memset(status, 0, sizeof(int)*locs.size());
 
  int found_all = 1; 
  for (int j = 0; j < locs.size(); ++j) {
    if(last[j] >= numElems) status[j] = 0;
    else {
      Elem& E = *elems[last[j]];
      status[j] = E.interpolateSolution(X, U, locs[j], sol[j], LSS, ghostPoints, varFcn);
    }
    if (!status[j]) found_all = 0;
  }

  if (found_all || assumeCache) 
    return;

  for (int i = 0; i < numElems; ++i)  {
    
    Elem& E = *elems[i];
    nn = E.numNodes(); 
    bbox[0] = Vec3D( std::numeric_limits<double>::max() ); 
    bbox[1] = Vec3D( -std::numeric_limits<double>::max() );
    for (int j = 0; j < nn; ++j) {
      const Vec3D& x = X[ E[j] ];
      bbox[0] = min( bbox[0], x );
      bbox[1] = max( bbox[1], x );
    }
    for (int j = 0; j < locs.size(); ++j) {
      if (!status[j]) { 
        if (bbox[0][0] <= locs[j][0] && bbox[1][0] >= locs[j][0] &&
            bbox[0][1] <= locs[j][1] && bbox[1][1] >= locs[j][1] &&
            bbox[0][2] <= locs[j][2] && bbox[1][2] >= locs[j][2]) {
          
          status[j] = E.interpolateSolution(X, U, locs[j], sol[j], LSS, ghostPoints, varFcn);
          if (status[j]) 
            last[j] = i;
        }
      }
    }
  } 
}

//------------------------------------------------------------------------------
