#ifndef _ELEM_H_
#define _ELEM_H_

#include <MacroCell.h>
#include <Vector.h>
#include <BlockAlloc.h>
#include <GhostPoint.h>
#include <VarFcn.h>

#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>

#ifdef OLD_STL
#include <map.h>
#include <vector.h>
#else
#include <map>
#include <vector>
using std::map;
#endif

/*
  The introduction of multiple element types in a code is traditionally 
  handled by creating a base class from which all other element types are 
  derived.
  
  In the base class, all of the required functions are declared as virtual 
  and then all the sub-classes can declare and implement their own variations
  of such functions.
  
  The problem we faced when wanting to generalize the code to any given set 
  of element (face) types is that the code heavily relies on templated 
  functions. C++ does not allow to have virtual templated functions. 
  Consequently we had to come up with an alternative approach.

  The approach we have used relies on the idea that a templated class can have 
  virtual functions. So if we have a templated function, that templated 
  function can construct a templated wrapper object that will handle the call 
  to the templated function on the actual element. The next question is how to 
  build that function when the type of the class is not specifically known at 
  compile time.  The way we do that is to call a virtual function to which we 
  pass a pointer to a Helper object. In fact that object is templated but each 
  of the templated classes is a subclass of a non-templated class that has a 
  virtual function to return the Wrapper we will need.

  To add an element of type E

  1) declare E as a subclass of A and give it all the required templated 
  functions.
  2) in E, add a function "void *getWrapper(int, char *memory)" That 
  function always has the exact same form. See in the example code.

  3) In GenHelper add the purely virtual function asClassE.
  4) in Helper: add the actual asClassE function. It is always the same 
  form and short

  To add a function:

  1) Add the function to the base class and in it create the helper, call 
  getWrapper on itself and call the function on the wrapper
  2) add the function to all subclasses.
*/

class FemEquationTerm;
class MacroCellSet;
class VMSLESTerm;
class DynamicVMSTerm;
class SmagorinskyLESTerm;
class WaleLESTerm;
class DynamicLESTerm;
class EdgeSet;
class BinFileHandler;
class GeoState;
class FaceSet;
class Face;
class PorousMedia;
class VolumeData;
class LevelSetStructure;

struct Vec3D;

template<class Scalar, int dim> class GenMat;
template<class Scalar, int dim, int dim2> class RectangularSparseMat;

class ElemTet;

/*
#ifdef OLD_STL
typedef map<Face, int, less<Face> > MapFaces;
#else
typedef map<Face, int> MapFaces;
#endif
*/
#include <MapFace.h>

//-------------- GENERAL HELPERS -----------------------------------------------
class GenElemHelper_dim {
public:
  virtual void *forClassTet(ElemTet *, int size, char *memorySpace) = 0;
};

class GenElemHelper_Scalar_dim_neq {
public:
  virtual void *forClassTet(ElemTet *, int size, char *memorySpace) = 0;
};

class GenElemHelper_dim_obj {
public:
  virtual void *forClassTet(ElemTet *, int size, char *memorySpace) = 0;
};

//-------------- GENERAL WRAPPERS ----------------------------------------------
template<int dim>
class GenElemWrapper_dim {
public:
  virtual 
  void computeGalerkinTerm(FemEquationTerm *, SVec<double,3> &, Vec<double> &,
									SVec<double,dim> &, SVec<double,dim> &, 
									Vec<GhostPoint<dim>*> * ghostPoints=0, 
									LevelSetStructure *LSS=0) = 0;

  virtual 
  void computeGalerkinTerm_e(FemEquationTerm *, SVec<double,3> &, Vec<double> &,
									  SVec<double,dim> &, SVec<double,dim> &, 
									  Vec<GhostPoint<dim>*> * ghostPoints=0, 
			   LevelSetStructure *LSS=0) = 0;
  
  virtual
  void computeVMSLESTerm(VMSLESTerm *, SVec<double,dim> &, SVec<double,3> &, 
			 SVec<double,dim> &, SVec<double,dim> &) = 0;
  
  virtual
  void computeMBarAndM(DynamicVMSTerm *, SVec<double,dim> **, SVec<double,1> **, 
		       SVec<double,3> &, SVec<double,dim> &,
		       SVec<double,dim> &, SVec<double,dim> &) = 0;  
  
  virtual
  void computeDynamicVMSTerm(DynamicVMSTerm *, SVec<double,dim> **, SVec<double,3> &,
			     SVec<double,dim> &, SVec<double,dim> &, Vec<double> &, Vec<double> &,
			     Vec<double> *, Vec<double> &) = 0;
  
  virtual
  void computeSmagorinskyLESTerm(SmagorinskyLESTerm *, SVec<double,3> &, SVec<double,dim> &V,
											SVec<double,dim> &R,  Vec<GhostPoint<dim>*> * ghostPoints=0, 
											LevelSetStructure *LSS=0) = 0;

  virtual
  void computeSmagorinskyLESTerm_e(SmagorinskyLESTerm *, SVec<double,3> &, SVec<double,dim> &V,
											SVec<double,dim> &R,  Vec<GhostPoint<dim>*> * ghostPoints=0, 
			         LevelSetStructure *LSS=0) = 0;

  virtual
  void computeWaleLESTerm(WaleLESTerm *, SVec<double,3> &, SVec<double,dim> &V, 
                          SVec<double,dim> &R, 
                          Vec<GhostPoint<dim>*> * ghostPoints=0, 
			  LevelSetStructure *LSS=0) = 0;
  
  virtual
  void computeWaleLESTerm_e(WaleLESTerm *, SVec<double,3> &, SVec<double,dim> &V, 
									 SVec<double,dim> &R, 
									 Vec<GhostPoint<dim>*> * ghostPoints=0, 
									 LevelSetStructure *LSS=0) = 0;
  
  virtual
  void computeDynamicLESTerm(DynamicLESTerm *, SVec<double,2> &, SVec<double,3> &, 
                             SVec<double,dim> &, SVec<double,dim> &, 
                             Vec<GhostPoint<dim>*> * ghostPoints=0, 
			     LevelSetStructure *LSS=0) = 0;    

  virtual
  void computeDynamicLESTerm_e(DynamicLESTerm *, SVec<double,2> &, SVec<double,3> &, 
										 SVec<double,dim> &, SVec<double,dim> &, 
										 Vec<GhostPoint<dim>*> * ghostPoints=0, 
										 LevelSetStructure *LSS=0) = 0;

  virtual
  void computeFaceGalerkinTerm(FemEquationTerm *, int [3], int, Vec3D &, 
			       SVec<double,3> &, Vec<double> &, double *, 
			       SVec<double,dim> &, SVec<double,dim> &) = 0;
  
  virtual
  void computeP1Avg(SVec<double,dim> &, SVec<double,16> &, SVec<double,6> &, Vec<double> &, 
                    SVec<double,8> &, SVec<double,3> &, SVec<double,dim> &, double, double,
                    Vec<GhostPoint<dim>*> * ghostPoints=0, 
		    LevelSetStructure *LSS=0) = 0;    
  
  virtual
  void computeP1Avg_e(SVec<double,dim> &, SVec<double,16> &, SVec<double,6> &, Vec<double> &, 
							 SVec<double,8> &, SVec<double,3> &, SVec<double,dim> &, double, double,
							 Vec<GhostPoint<dim>*> * ghostPoints=0, 
							 LevelSetStructure *LSS=0) = 0;  

// Included (MB)
  virtual
  void computeDerivativeOfGalerkinTerm(FemEquationTerm *, SVec<double,3> &, SVec<double,3> &, Vec<double> &,
			   SVec<double,dim> &, SVec<double,dim> &, double, SVec<double,dim> &) = 0;

  virtual
  void computeDerivativeOperatorsOfGalerkinTerm(FemEquationTerm *, SVec<double,3> &, Vec<double> &,
			   SVec<double,dim> &, RectangularSparseMat<double,3,dim> &) = 0; // YC

  virtual
  void computeDerivativeOfFaceGalerkinTerm(FemEquationTerm *, int [3], int, Vec3D &, Vec3D &,
			       SVec<double,3> &, SVec<double,3> &, Vec<double> &, double *, double *,
			       SVec<double,dim> &, SVec<double,dim> &, double, SVec<double,dim> &) = 0;

// Level Set Reinitialization

  virtual
  void computeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                 SVec<double,dim> &ddx, SVec<double,dim> &ddy,
                                 SVec<double,dim> &ddz,
                                 SVec<double,dim> &Phi,SVec<double,1> &Psi) = 0;

  virtual
  void recomputeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                 SVec<double,dim> &ddx, SVec<double,dim> &ddy,
                                 SVec<double,dim> &ddz,
                                 SVec<double,dim> &Phi,SVec<double,1> &Psi) = 0;

  virtual
  void computeDistanceLevelNodes(int lsdim, Vec<int> &Tag, int level,
                                 SVec<double,3> &X, SVec<double,1> &Psi, SVec<double,dim> &Phi) = 0;
  virtual
  void FastMarchingDistanceUpdate(int node, Vec<int> &Tag, int level,
                              SVec<double,3> &X,SVec<double,dim> &d2wall) = 0;


  // X is the deformed nodal location vector
  virtual
  int interpolateSolution(SVec<double,3>& X, SVec<double,dim>& U, const Vec3D& loc, double sol[dim], LevelSetStructure* LSS,
                          Vec<GhostPoint<dim>*>* ghostPoints, VarFcn* varFcn) = 0;

};

template<class Scalar, int dim, int neq>
class GenElemWrapper_Scalar_dim_neq {
public:
  virtual 
  void computeJacobianGalerkinTerm(FemEquationTerm *, SVec<double,3> &, 
				   Vec<double> &, Vec<double> &, 
				   SVec<double,dim> &, GenMat<Scalar,neq> &,
                                   Vec<GhostPoint<dim>*>*gp=0,LevelSetStructure *LSS=0) = 0;
  
  virtual 
  void computeJacobianGalerkinTerm_e(FemEquationTerm *, SVec<double,3> &, 
												 Vec<double> &, Vec<double> &, 
												 SVec<double,dim> &, GenMat<Scalar,neq> &,
												 Vec<GhostPoint<dim>*>*gp=0,LevelSetStructure *LSS=0) = 0;
  
  virtual 
  void computeFaceJacobianGalerkinTerm(FemEquationTerm *, int [3], int, Vec3D &, 
				       SVec<double,3> &, Vec<double> &, Vec<double> &, 
				       double *, SVec<double,dim> &, GenMat<Scalar,neq> &) = 0;

};

template<int dim, class Obj>
class GenElemWrapper_dim_obj {
public:
  virtual void integrateFunction(Obj* obj,SVec<double,3> &X,SVec<double,dim>& V, void (Obj::*F)(int node, const double* loc,double* f),int) = 0; 
  
};


//-------------- REAL WRAPPERS -------------------------------------------------
template<class Target, int dim>
class  ElemWrapper_dim : public GenElemWrapper_dim<dim> {
  Target *t;

public:
  ElemWrapper_dim(Target *tt) : t(tt) { };

  void computeGalerkinTerm(FemEquationTerm *fet, SVec<double,3> &X, 
			   Vec<double> &d2wall, SVec<double,dim> &V, 
									SVec<double,dim> &R, 
									Vec<GhostPoint<dim>*> *ghostPoints=0, 
									LevelSetStructure *LSS=0) {
    t->computeGalerkinTerm(fet, X, d2wall, V, R,ghostPoints,LSS);
  }
  
  void computeGalerkinTerm_e(FemEquationTerm *fet, SVec<double,3> &X, 
									  Vec<double> &d2wall, SVec<double,dim> &V, 
									  SVec<double,dim> &R, 
									  Vec<GhostPoint<dim>*> *ghostPoints=0, 
									  LevelSetStructure *LSS=0) {
	  t->computeGalerkinTerm_e(fet, X, d2wall, V, R, ghostPoints, LSS);
  }
  

  void computeVMSLESTerm(VMSLESTerm *vmst, SVec<double,dim> &VBar,
			 SVec<double,3> &X, SVec<double,dim> &V,
			 SVec<double,dim> &Sigma) {
    t->computeVMSLESTerm(vmst, VBar, X, V, Sigma);
  }
  
  void computeMBarAndM(DynamicVMSTerm *dvmst, SVec<double,dim> **VBar,
		       SVec<double,1> **volRatio, SVec<double,3> &X,
		       SVec<double,dim> &V, SVec<double,dim> &MBar,
		       SVec<double,dim> &M) {
    t->computeMBarAndM(dvmst,VBar, volRatio, X, V, MBar, M);
  }
  
  void computeDynamicVMSTerm(DynamicVMSTerm *dvmst, SVec<double,dim> **VBar,
			     SVec<double,3> &X, SVec<double,dim> &V,
			     SVec<double,dim> &S, Vec<double> &CsDelSq,
			     Vec<double> &PrT, Vec<double> *Cs,
			     Vec<double> &Delta) {
    t->computeDynamicVMSTerm(dvmst, VBar, X, V, S, CsDelSq, PrT, Cs, Delta);
  }
  
  void computeSmagorinskyLESTerm(SmagorinskyLESTerm *smag, SVec<double,3> &X,
				 SVec<double,dim> &V, SVec<double,dim> &R,
											Vec<GhostPoint<dim>*> *ghostPoints=0, 
											LevelSetStructure *LSS=0) {
    t->computeSmagorinskyLESTerm(smag, X, V, R, ghostPoints, LSS);
  }

  void computeSmagorinskyLESTerm_e(SmagorinskyLESTerm *smag, SVec<double,3> &X,
											  SVec<double,dim> &V, SVec<double,dim> &R,
											  Vec<GhostPoint<dim>*> *ghostPoints=0, 
											  LevelSetStructure *LSS=0) {
	  t->computeSmagorinskyLESTerm_e(smag, X, V, R, ghostPoints, LSS);
  }


  void computeWaleLESTerm(WaleLESTerm *wale, SVec<double,3> &X,
		          SVec<double,dim> &V, SVec<double,dim> &R,
			  Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0) {
    t->computeWaleLESTerm(wale, X, V, R, ghostPoints, LSS);
  }  

  void computeWaleLESTerm_e(WaleLESTerm *wale, SVec<double,3> &X,
									 SVec<double,dim> &V, SVec<double,dim> &R,
									 Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0) {
	  t->computeWaleLESTerm_e(wale, X, V, R, ghostPoints, LSS); 
  }  

  void computeDynamicLESTerm(DynamicLESTerm *dles, SVec<double,2> &Cs, 
                   SVec<double,3> &X, SVec<double,dim> &V, SVec<double,dim> &R,
	           Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0) {
    t->computeDynamicLESTerm(dles, Cs, X, V, R, ghostPoints, LSS);
  }

  void computeDynamicLESTerm_e(DynamicLESTerm *dles, SVec<double,2> &Cs, 
										 SVec<double,3> &X, SVec<double,dim> &V, SVec<double,dim> &R,
										 Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0) {
	  t->computeDynamicLESTerm_e(dles, Cs, X, V, R, ghostPoints, LSS);
  }
  
  void computeFaceGalerkinTerm(FemEquationTerm *fet, int face[3], int code, Vec3D &n, 
			       SVec<double,3> &X, Vec<double> &d2wall, double *Vwall, 
			       SVec<double,dim> &V, SVec<double,dim> &R) {
    t->computeFaceGalerkinTerm(fet, face, code, n, X, d2wall, Vwall, V, R);
  }


  void computeP1Avg(SVec<double,dim> &VCap, SVec<double,16> &Mom_Test, SVec<double,6> &Sij_Test, 
                    Vec<double> &modS_Test, SVec<double,8> &Eng_Test, SVec<double,3> &X, 
                    SVec<double,dim> &V, double gam, double R,
	            Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0) {
    t->computeP1Avg(VCap, Mom_Test, Sij_Test, modS_Test, Eng_Test, X, V, gam, R, ghostPoints, LSS);
  }
  

  void computeP1Avg_e(SVec<double,dim> &VCap, SVec<double,16> &Mom_Test, SVec<double,6> &Sij_Test, 
                    Vec<double> &modS_Test, SVec<double,8> &Eng_Test, SVec<double,3> &X, 
                    SVec<double,dim> &V, double gam, double R,
						  Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0) {
	  t->computeP1Avg_e(VCap, Mom_Test, Sij_Test, modS_Test, Eng_Test, X, V, gam, R, ghostPoints, LSS);
  }
  

// Included (MB)
  void computeDerivativeOfGalerkinTerm(FemEquationTerm *fet, SVec<double,3> &X, SVec<double,3> &dX,
			      Vec<double> &d2wall, SVec<double,dim> &V, SVec<double,dim> &dV, double dMach,
			      SVec<double,dim> &dR) {
    t->computeDerivativeOfGalerkinTerm(fet, X, dX, d2wall, V, dV, dMach, dR);
  }

  //YC
  void computeDerivativeOperatorsOfGalerkinTerm(FemEquationTerm *fet, SVec<double,3> &X,
             Vec<double> &d2wall, SVec<double,dim> &V, RectangularSparseMat<double,3,dim> &dViscousFluxdX) {
     t->computeDerivativeOperatorsOfGalerkinTerm(fet, X, d2wall, V, dViscousFluxdX);
  }

  void computeDerivativeOfFaceGalerkinTerm(FemEquationTerm *fet, int face[3], int code, Vec3D &n, Vec3D &dn,
				  SVec<double,3> &X, SVec<double,3> &dX, Vec<double> &d2wall, double *Vwall, double *dVwall,
				  SVec<double,dim> &V, SVec<double,dim> &dV, double dMach, SVec<double,dim> &dR) {
    t->computeDerivativeOfFaceGalerkinTerm(fet, face, code, n, dn, X, dX, d2wall, Vwall, dVwall, V, dV, dMach, dR);
  }

// Level Set Reinitialization

  void computeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                 SVec<double,dim> &ddx, SVec<double,dim> &ddy,
                                 SVec<double,dim> &ddz,
                                 SVec<double,dim> &Phi,SVec<double,1> &Psi){
    t->computeDistanceCloseNodes(lsdim,Tag,X,ddx,ddy,ddz,Phi,Psi);
  }

  void recomputeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                 SVec<double,dim> &ddx, SVec<double,dim> &ddy,
                                 SVec<double,dim> &ddz,
                                 SVec<double,dim> &Phi,SVec<double,1> &Psi){
    t->recomputeDistanceCloseNodes(lsdim,Tag,X,ddx,ddy,ddz,Phi,Psi);
  }
  void FastMarchingDistanceUpdate(int node, Vec<int> &Tag, int level,
                              SVec<double,3> &X,SVec<double,dim> &d2wall){
    t->FastMarchingDistanceUpdate(node,Tag,level,X,d2wall);
  }

  void computeDistanceLevelNodes(int lsdim, Vec<int> &Tag, int level,
                                 SVec<double,3> &X, SVec<double,1> &Psi, SVec<double,dim> &Phi){
    t->computeDistanceLevelNodes(lsdim,Tag,level,X,Psi,Phi);
  }

  // X is the deformed nodal location vector
  int interpolateSolution(SVec<double,3>& X, SVec<double,dim>& U, const Vec3D& loc,
			  double sol[dim], LevelSetStructure* LSS,
                          Vec<GhostPoint<dim>*>* ghostPoints, VarFcn* varFcn) {
    return t->interpolateSolution(X,U,loc,sol,LSS,ghostPoints,varFcn);
  }

};

template<class Target, class Scalar, int dim, int neq>
class  ElemWrapper_Scalar_dim_neq : public 
GenElemWrapper_Scalar_dim_neq<Scalar,dim,neq> {
  
  Target *t;
  
public:
  ElemWrapper_Scalar_dim_neq(Target *tt) : t(tt) { };
  
  void computeJacobianGalerkinTerm(FemEquationTerm *fet, SVec<double,3> &X, 
				   Vec<double> &ctrlVol, Vec<double> &d2wall, 
				   SVec<double,dim> &V, GenMat<Scalar,neq> &A,
				   Vec<GhostPoint<dim>*> *gp=0,LevelSetStructure *LSS=0) {
    t->computeJacobianGalerkinTerm(fet, X, ctrlVol, d2wall, V, A,gp,LSS);
  }
  
  void computeJacobianGalerkinTerm_e(FemEquationTerm *fet, SVec<double,3> &X, 
												 Vec<double> &ctrlVol, Vec<double> &d2wall, 
												 SVec<double,dim> &V, GenMat<Scalar,neq> &A,
												 Vec<GhostPoint<dim>*> *gp=0,LevelSetStructure *LSS=0) {
	  t->computeJacobianGalerkinTerm_e(fet, X, ctrlVol, d2wall, V, A,gp,LSS);
  }
  
  void computeFaceJacobianGalerkinTerm(FemEquationTerm *fet, int face[3], int code, 
				       Vec3D &n, SVec<double,3> &X, Vec<double> &ctrlVol,
				       Vec<double> &d2wall, double *Vwall, 
				       SVec<double,dim> &V, GenMat<Scalar,neq> &A) {
    t->computeFaceJacobianGalerkinTerm(fet, face, code, n, X, ctrlVol,
				       d2wall, Vwall, V, A);
  }
  
};

template<class Target,int dim, class Obj>
class  ElemWrapper_dim_obj : public 
GenElemWrapper_dim_obj<dim,Obj> {
  
  Target *t;
  
public:
  ElemWrapper_dim_obj(Target *tt) : t(tt) { };
  
    void integrateFunction(Obj* obj,SVec<double,3> &X,SVec<double,dim>& V, void (Obj::*F)(int node, const double* loc,double* f),int npt) {
      t->integrateFunction(obj,X,V,F,npt);
  }
};

//-------------- REAL HELPERS --------------------------------------------------
template<int dim>
class ElemHelper_dim : public GenElemHelper_dim {
public:

  void *forClassTet(ElemTet *tet, int size, char *memorySpace) {
    if(size < sizeof(ElemWrapper_dim<ElemTet, dim>) ) {
      fprintf(stderr, "Error: programming error in ElemHelper");
      exit(1);
    }
    return new (memorySpace) ElemWrapper_dim<ElemTet, dim>(tet);
  }

};

template<class Scalar, int dim, int neq>
class ElemHelper_Scalar_dim_neq : public GenElemHelper_Scalar_dim_neq {
public:

  void *forClassTet(ElemTet *tet, int size, char *memorySpace) {
    if(size < sizeof(ElemWrapper_Scalar_dim_neq<ElemTet, Scalar, dim, neq>) ) {
      fprintf(stderr, "Error: programming error in ElemHelper");
      exit(1);
    }
    return new (memorySpace) ElemWrapper_Scalar_dim_neq<ElemTet, Scalar, dim, neq>(tet);
  }

};

  template<int dim, class Obj>
class ElemHelper_dim_obj : public GenElemHelper_dim_obj {
public:

  void *forClassTet(ElemTet *tet, int size, char *memorySpace) {
    if(size < sizeof(ElemWrapper_dim_obj<ElemTet, dim, Obj>) ) {
      fprintf(stderr, "Error: programming error in ElemHelper");
      exit(1);
    }
    return new (memorySpace) ElemWrapper_dim_obj<ElemTet, dim, Obj>(tet);
  }
};


//--------------- BASE ELEM CLASS ----------------------------------------------
class Elem {
public:
  enum Type {TET=5};

  static const int MaxNumNd = 8;
  static const int MaxNumEd = 12;
  static const int MaxNumFc = 6;

protected:
  virtual void *getWrapper_dim(GenElemHelper_dim *, 
			       int size, char *memorySpace) = 0;
  virtual void *getWrapper_Scalar_dim_neq(GenElemHelper_Scalar_dim_neq *, 
					  int size, char *memorySpace) = 0;  
  virtual void *getWrapper_dim_obj(GenElemHelper_dim_obj *, 
				   int size, char *memorySpace) = 0;
  int volume_id;

public:

  virtual int* nodeNum() = 0;
  virtual int& nodeNum(int i) = 0;
  virtual int& edgeNum(int i) = 0;
  virtual int  edgeEnd(int i, int k) = 0;
  virtual int  edgeFace(int i, int k) = 0;
  virtual int  faceDef(int i, int k) = 0;
  virtual int  faceNnd(int i) = 0;
  virtual Type type() = 0;

  virtual bool isPointInside(SVec<double,3>& X, const Vec3D&) = 0;

  void computeBoundingBox(SVec<double,3>& X, double*);

  // Number of nodes
  virtual int numNodes() = 0;

  // Number of edges
  virtual int numEdges() = 0;

  // Number of faces
  virtual int numFaces() = 0;

  // Operator Elem[i]: get reference to node i
  int &operator[](int i) { return nodeNum(i); }  

  // Operator *(Elem): get pointer to node list
  operator int *() { return nodeNum(); }

  // fills the nd array with the list of nodes
  void nodes(int *nd) { for(int j=0; j<numNodes(); ++j) nd[j] = nodeNum(j); }
  
  // Renumber nodes according to nodemap
  template<class NodeMap>
  void renumberNodes(NodeMap &nodemap) {
    for (int j=0; j<numNodes(); ++j)
      nodeNum(j) = nodemap[ nodeNum(j) ];
  }

  // Set element volume ID
  void setVolumeID(int i) {volume_id = i;}

  // Get element volume ID
  int getVolumeID() { return volume_id; }


  //--------------Functions in ElemCore.C

  double computeLongestEdge(SVec<double,3> &X);
  void numberEdges(EdgeSet &);
  void renumberEdges(Vec<int> &);
  int  countNodesOnBoundaries(Vec<bool> &);
  int  setFaceToElementConnectivity(int i, Vec<bool> &, MapFaces &, int (*)[2]);


  //--------------Virtual functions in ElemTYPECore.C

  virtual double computeVolume(SVec<double,3> &) = 0;
  virtual double computeControlVolumes(SVec<double,3> &, Vec<double> &) = 0;
  virtual void printInvalidElement(int, double, int, int *, int *, 
				   SVec<double,3> &, SVec<double,3> &) = 0;
  virtual void computeEdgeNormalsConfig(SVec<double,3> &Xconfig, SVec<double,3> &Xdot,
                                        Vec<Vec3D> &edgeNorm, Vec<double> &edgeNormVel) = 0;
  virtual void computeEdgeNormalsGCL1(SVec<double,3> &, SVec<double,3> &, SVec<double,3> &, 
				      Vec<Vec3D> &, Vec<double> &) = 0;
  virtual void computeEdgeNormalsEZGCL1(double, SVec<double,3> &, SVec<double,3> &, 
					Vec<Vec3D> &, Vec<double> &) = 0;
  virtual void computeWeightsGalerkin(SVec<double,3> &, SVec<double,3> &, 
				      SVec<double,3> &, SVec<double,3> &) = 0;
  virtual void computeEdgeWeightsGalerkin(SVec<double,3> &, SVec<double,9> &) = 0;
  
  virtual double computeGradientP1Function(SVec<double,3> &, double (*)[3], 
					   double * = NULL) = 0;
  virtual void computeStiffAndForce(double *, double *, 
				    SVec<double, 3> &, SVec<double,3> &, 
				    double volStiff = 0.0) = 0;
  virtual void computeStiffAndForceBallVertex(double *f, double *K, 
				    SVec<double, 3> &X, SVec<double,3> &X0, 
				    double volStiff = 0.0) = 0;
  virtual void computeStiffAndForceLIN(double *, SVec<double,3> &, SVec<double,3> &) = 0;

  virtual void computeStiffBallVertex(double *, SVec<double,3> &X, SVec<double,3> &X0, double volStiff) = 0;
  virtual void computeStiffTorsionSpring(double *, SVec<double,3> &, double volStiff) = 0;
  
// Included (MB)
  virtual double computeDerivativeOfVolume(SVec<double,3> &, SVec<double,3> &) = 0;
  virtual void computeDerivativeOperatorsOfVolume(SVec<double,3> &, double [][3], double [][3], double [][3], double [][3]) = 0;

  virtual double computeDerivativeOfControlVolumes(SVec<double,3> &, SVec<double,3> &, Vec<double> &) = 0;
  virtual void computeDerivativeOperatorsOfControlVolumes(SVec<double,3> &, RectangularSparseMat<double,3,1> &) = 0;

  virtual void computeDerivativeOfEdgeNormals(SVec<double,3> &, SVec<double,3> &, Vec<Vec3D> &, Vec<Vec3D> &, Vec<double> &, Vec<double> &) = 0;
  virtual void computeDerivativeOperatorsOfEdgeNormals(SVec<double,3> &, RectangularSparseMat<double,3,3> &) = 0;

  virtual void computeDerivativeOfWeightsGalerkin(SVec<double,3> &, SVec<double,3> &, SVec<double,3> &, SVec<double,3> &, SVec<double,3> &) = 0;
  virtual void computeDerivativeTransposeOfWeightsGalerkin(SVec<double,3> &, SVec<double,3> &, SVec<double,3> &, SVec<double,3> &, SVec<double,3> &) = 0;

  virtual double computeDerivativeOfGradientP1Function(SVec<double,3> &, SVec<double,3> &, double [4][3]) = 0;
  virtual double computeDerivativeOfGradientP1Function2(SVec<double,3> &, SVec<double,3> &, double [4][3], double dX[4][3]) = 0;
  virtual void  computeDerivativeOperatorOfGradientP1Function(SVec<double,3> &, double [4][3], double[4][3][4][3], int [4]) = 0;
  virtual void computeDerivativeTransposeOfGradientP1Function(SVec<double,3> &, double, double [4][3], double [4][3], SVec<double,3> &) = 0;

  virtual void computeBarycentricCoordinates(SVec<double,3>&, const Vec3D& , double [3]) = 0;

  //-----Virtual template functions (handled through helpers classes, defined ElemTYPE.C)

  template<int dim>
  void computeGalerkinTerm(FemEquationTerm *fet, SVec<double,3> &X, 
			   Vec<double> &d2wall, SVec<double,dim> &V, 
									SVec<double,dim> &R, 
									Vec<GhostPoint<dim>*> *ghostPoints=0, 
									LevelSetStructure *LSS=0) {
    ElemHelper_dim<dim> h;
    char xx[64];
	  GenElemWrapper_dim<dim> *wrapper = (GenElemWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeGalerkinTerm(fet, X, d2wall, V, R,ghostPoints,LSS);
  }
  
  template<int dim>
	  void computeGalerkinTerm_e(FemEquationTerm *fet, SVec<double,3> &X,
										  Vec<double> &d2wall, SVec<double,dim> &V, 
										  SVec<double,dim> &R,
										  Vec<GhostPoint<dim>*> *ghostPoints=0, 
										  LevelSetStructure *LSS=0) {
	  ElemHelper_dim<dim> h;
	  char xx[64];
	  GenElemWrapper_dim<dim> *wrapper = (GenElemWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
	  wrapper->computeGalerkinTerm_e(fet, X, d2wall, V, R,ghostPoints, LSS);
  }
  
  template<int dim>
  void computeVMSLESTerm(VMSLESTerm *vmst, SVec<double,dim> &VBar,
			 SVec<double,3> &X, SVec<double,dim> &V,
			 SVec<double,dim> &Sigma) {
    ElemHelper_dim<dim> h;
    char xx[64];
    GenElemWrapper_dim<dim> *wrapper=
      (GenElemWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeVMSLESTerm(vmst, VBar, X, V, Sigma);
  }
  
  template<int dim>
  void computeMBarAndM(DynamicVMSTerm *dvmst, SVec<double,dim> **VBar,
		       SVec<double,1> **volRatio, SVec<double,3> &X,
		       SVec<double,dim> &V, SVec<double,dim> &MBar,
		       SVec<double,dim> &M) {
    ElemHelper_dim<dim> h;
    char xx[64];
    GenElemWrapper_dim<dim> *wrapper=
      (GenElemWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeMBarAndM(dvmst,VBar, volRatio, X, V, MBar, M);
  }
  
  template<int dim>
  void computeDynamicVMSTerm(DynamicVMSTerm *dvmst, SVec<double,dim> **VBar,
			     SVec<double,3> &X, SVec<double,dim> &V,
			     SVec<double,dim> &S, Vec<double> &CsDelSq,
			     Vec<double> &PrT, Vec<double> *Cs,
			     Vec<double> &Delta) {
    ElemHelper_dim<dim> h;
    char xx[64];
    GenElemWrapper_dim<dim> *wrapper=
      (GenElemWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDynamicVMSTerm(dvmst, VBar, X, V, S, CsDelSq, PrT, Cs, Delta);
  }
  
  template<int dim>
  void computeSmagorinskyLESTerm(SmagorinskyLESTerm *smag, SVec<double,3> &X,
				 SVec<double,dim> &V, SVec<double,dim> &R, 
                                 Vec<GhostPoint<dim>*> *ghostPoints=0, 
											LevelSetStructure *LSS=0) 
  {
    ElemHelper_dim<dim> h;
    char xx[64];
    GenElemWrapper_dim<dim> *wrapper=
      (GenElemWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeSmagorinskyLESTerm(smag, X, V, R, ghostPoints, LSS);
  }
 
  template<int dim>
  void computeSmagorinskyLESTerm_e(SmagorinskyLESTerm *smag, SVec<double,3> &X,
											  SVec<double,dim> &V, SVec<double,dim> &R, 
											  Vec<GhostPoint<dim>*> *ghostPoints=0, 
											  LevelSetStructure *LSS=0) 
  {
	  ElemHelper_dim<dim> h;
	  char xx[64];
	  GenElemWrapper_dim<dim> *wrapper=
		  (GenElemWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
	  wrapper->computeSmagorinskyLESTerm_e(smag, X, V, R, ghostPoints, LSS);
  }
 
 
   template<int dim>
   void computeWaleLESTerm(WaleLESTerm *wale, SVec<double,3> &X,
		           SVec<double,dim> &V, SVec<double,dim> &R,
                           Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0) {
    ElemHelper_dim<dim> h;
    char xx[64];
    GenElemWrapper_dim<dim> *wrapper=
      (GenElemWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeWaleLESTerm(wale, X, V, R, ghostPoints, LSS);
  }
 
	template<int dim>
		void computeWaleLESTerm_e(WaleLESTerm *wale, SVec<double,3> &X,
										SVec<double,dim> &V, SVec<double,dim> &R,
										Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0) {
		ElemHelper_dim<dim> h;
		char xx[64];
		GenElemWrapper_dim<dim> *wrapper=
			(GenElemWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
		wrapper->computeWaleLESTerm_e(wale, X, V, R, ghostPoints, LSS);
	}
 

  template<int dim>
  void computeDynamicLESTerm(DynamicLESTerm *dles, SVec<double,2> &Cs, 
                             SVec<double,3> &X, SVec<double,dim> &V, SVec<double,dim> &R,  
                             Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0) {
    ElemHelper_dim<dim> h;
    char xx[64];
    GenElemWrapper_dim<dim> *wrapper=
      (GenElemWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDynamicLESTerm(dles, Cs, X, V, R, ghostPoints, LSS);
  }
  
  template<int dim>
  void computeDynamicLESTerm_e(DynamicLESTerm *dles, SVec<double,2> &Cs, 
										 SVec<double,3> &X, SVec<double,dim> &V, SVec<double,dim> &R,  
										 Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0) {
    ElemHelper_dim<dim> h;
    char xx[64];
    GenElemWrapper_dim<dim> *wrapper=
		 (GenElemWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDynamicLESTerm_e(dles, Cs, X, V, R, ghostPoints, LSS);
  }
  
  template<int dim>
  void computeFaceGalerkinTerm(FemEquationTerm *fet, int face[3], int code, Vec3D &n, 
			       SVec<double,3> &X, Vec<double> &d2wall, double *Vwall, 
			       SVec<double,dim> &V, SVec<double,dim> &R) {
    ElemHelper_dim<dim> h;
    char xx[64];
    GenElemWrapper_dim<dim> *wrapper=
      (GenElemWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeFaceGalerkinTerm(fet, face, code, n, X, d2wall, Vwall, V, R);
  }


  template<int dim>
  void computeP1Avg(SVec<double,dim> &VCap, SVec<double,16> &Mom_Test, SVec<double,6> &Sij_Test, 
                       Vec<double> &modS_Test, SVec<double,8> &Eng_Test, SVec<double,3> &X, 
                       SVec<double,dim> &V, double gam, double R,
                       Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0) {
    ElemHelper_dim<dim> h;
    char xx[64];
    GenElemWrapper_dim<dim> *wrapper=
      (GenElemWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeP1Avg(VCap, Mom_Test, Sij_Test, modS_Test, Eng_Test, X, V, gam, R, ghostPoints, LSS);
  }  

  template<int dim>
  void computeP1Avg_e(SVec<double,dim> &VCap, SVec<double,16> &Mom_Test, SVec<double,6> &Sij_Test, 
							 Vec<double> &modS_Test, SVec<double,8> &Eng_Test, SVec<double,3> &X, 
							 SVec<double,dim> &V, double gam, double R,
							 Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0) {
	  ElemHelper_dim<dim> h;
	  char xx[64];
	  GenElemWrapper_dim<dim> *wrapper=
		  (GenElemWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
	  wrapper->computeP1Avg_e(VCap, Mom_Test, Sij_Test, modS_Test, Eng_Test, X, V, gam, R, ghostPoints, LSS);
  }  

  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(FemEquationTerm *fet, SVec<double,3> &X, 
				   Vec<double> &ctrlVol, Vec<double> &d2wall, 
				   SVec<double,dim> &V, GenMat<Scalar,neq> &A,
				   Vec<GhostPoint<dim>*> *gp=0,LevelSetStructure *LSS=0) {
    ElemHelper_Scalar_dim_neq<Scalar, dim, neq> h;
    char xx[64];
    GenElemWrapper_Scalar_dim_neq<Scalar, dim, neq> *wrapper=
      (GenElemWrapper_Scalar_dim_neq<Scalar, dim, neq> *)getWrapper_Scalar_dim_neq(&h, 64, xx);
    wrapper->computeJacobianGalerkinTerm(fet, X, ctrlVol, d2wall, V, A,gp,LSS);
  }
  
  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm_e(FemEquationTerm *fet, SVec<double,3> &X, 
												 Vec<double> &ctrlVol, Vec<double> &d2wall, 
												 SVec<double,dim> &V, GenMat<Scalar,neq> &A,
												 Vec<GhostPoint<dim>*> *gp=0,LevelSetStructure *LSS=0) {
	  ElemHelper_Scalar_dim_neq<Scalar, dim, neq> h;
	  char xx[64];
	  GenElemWrapper_Scalar_dim_neq<Scalar, dim, neq> *wrapper=
		  (GenElemWrapper_Scalar_dim_neq<Scalar, dim, neq> *)getWrapper_Scalar_dim_neq(&h, 64, xx);
	  wrapper->computeJacobianGalerkinTerm_e(fet, X, ctrlVol, d2wall, V, A,gp,LSS);
  }

  template<int dim, class Scalar, int neq>
  void computeFaceJacobianGalerkinTerm(FemEquationTerm *fet, int face[3], int code, 
				       Vec3D &n, SVec<double,3> &X, Vec<double> &ctrlVol,
				       Vec<double> &d2wall, double *Vwall, 
				       SVec<double,dim> &V, GenMat<Scalar,neq> &A) {
    ElemHelper_Scalar_dim_neq<Scalar, dim, neq> h;
    char xx[64];
    GenElemWrapper_Scalar_dim_neq<Scalar, dim, neq> *wrapper=
      (GenElemWrapper_Scalar_dim_neq<Scalar, dim, neq> *)getWrapper_Scalar_dim_neq(&h, 64, xx);
    wrapper->computeFaceJacobianGalerkinTerm(fet, face, code, n, X, ctrlVol,
					     d2wall, Vwall, V, A);
  }

// Included (MB)
  template<int dim>
  void computeDerivativeOfGalerkinTerm(FemEquationTerm *fet, SVec<double,3> &X, SVec<double,3> &dX,
			      Vec<double> &d2wall, SVec<double,dim> &V, SVec<double,dim> &dV, double dMach,
			      SVec<double,dim> &dR) {
    ElemHelper_dim<dim> h;
    char xx[64];
    GenElemWrapper_dim<dim> *wrapper=
      (GenElemWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDerivativeOfGalerkinTerm(fet, X, dX, d2wall, V, dV, dMach, dR);
  }

  template<int dim>
  void computeDerivativeOperatorsOfGalerkinTerm(FemEquationTerm *fet, SVec<double,3> &X,
			      Vec<double> &d2wall, SVec<double,dim> &V, RectangularSparseMat<double,3,dim> &dViscousFluxdX) {  //YC
    ElemHelper_dim<dim> h;
    char xx[64];
    GenElemWrapper_dim<dim> *wrapper=
      (GenElemWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDerivativeOperatorsOfGalerkinTerm(fet, X, d2wall, V, dViscousFluxdX);
  }

  template<int dim>
  void computeDerivativeOfFaceGalerkinTerm(FemEquationTerm *fet, int face[3], int code, Vec3D &n, Vec3D &dn,
				  SVec<double,3> &X, SVec<double,3> &dX, Vec<double> &d2wall, double *Vwall, double *dVwall,
				  SVec<double,dim> &V, SVec<double,dim> &dV, double dMach, SVec<double,dim> &dR) {
    ElemHelper_dim<dim> h;
    char xx[64];
    GenElemWrapper_dim<dim> *wrapper=
      (GenElemWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDerivativeOfFaceGalerkinTerm(fet, face, code, n, dn, X, dX, d2wall, Vwall, dVwall, V, dV, dMach, dR);
  }
  
  // X is the deformed nodal location vector
  template<int dim> 
  int interpolateSolution(SVec<double,3>& X, SVec<double,dim>& U, const Vec3D& loc, double sol[dim], LevelSetStructure* LSS,
                          Vec<GhostPoint<dim>*>* ghostPoints, VarFcn* varFcn) {

    ElemHelper_dim<dim> h;
    char xx[64];
    GenElemWrapper_dim<dim> *wrapper=
      (GenElemWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    return wrapper->interpolateSolution(X,U,loc,sol,LSS,ghostPoints,varFcn);
  }
  
// Level Set Reinitialization

  template<int dimLS>
  void computeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                 SVec<double,dimLS> &ddx, SVec<double,dimLS> &ddy,
                                 SVec<double,dimLS> &ddz,
                                 SVec<double,dimLS> &Phi,SVec<double,1> &Psi){
    ElemHelper_dim<dimLS> h;
    char xx[64];
    GenElemWrapper_dim<dimLS> *wrapper=
      (GenElemWrapper_dim<dimLS> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDistanceCloseNodes(lsdim,Tag,X,ddx,ddy,ddz,Phi,Psi);
  }

  template<int dimLS>
  void recomputeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                 SVec<double,dimLS> &ddx, SVec<double,dimLS> &ddy,
                                 SVec<double,dimLS> &ddz,
                                 SVec<double,dimLS> &Phi,SVec<double,1> &Psi){
    ElemHelper_dim<dimLS> h;
    char xx[64];
    GenElemWrapper_dim<dimLS> *wrapper=
      (GenElemWrapper_dim<dimLS> *)getWrapper_dim(&h, 64, xx);
    wrapper->recomputeDistanceCloseNodes(lsdim,Tag,X,ddx,ddy,ddz,Phi,Psi);
  }
  template<int dimLS>
  void FastMarchingDistanceUpdate(int node, Vec<int> &Tag, int level,
                               SVec<double,3> &X,SVec<double,dimLS> &d2wall)
  {
    ElemHelper_dim<dimLS> h;
    char xx[64];
    GenElemWrapper_dim<dimLS> *wrapper=
      (GenElemWrapper_dim<dimLS> *)getWrapper_dim(&h, 64, xx);
    wrapper->FastMarchingDistanceUpdate(node,Tag,level,X,d2wall);
  }

  template<int dimLS>
  void computeDistanceLevelNodes(int lsdim, Vec<int> &Tag, int level,
                                 SVec<double,3> &X, SVec<double,1> &Psi, SVec<double,dimLS> &Phi){
    ElemHelper_dim<dimLS> h;
    char xx[64];
    GenElemWrapper_dim<dimLS> *wrapper=
      (GenElemWrapper_dim<dimLS> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDistanceLevelNodes(lsdim,Tag,level,X,Psi,Phi);
  }

  template<int dim, class Obj>
    void integrateFunction(Obj* obj,SVec<double,3> &X,SVec<double,dim>& V, void (Obj::*F)(int node, const double* loc,double* f),int npt) { 
    ElemHelper_dim_obj<dim,Obj> h;
    char xx[64];
    GenElemWrapper_dim_obj<dim,Obj> *wrapper=
      (GenElemWrapper_dim_obj<dim,Obj> *)getWrapper_dim_obj(&h, 64, xx);
    wrapper->integrateFunction(obj,X,V,F,npt);
  }
};

//--------------- DUMMY ELEM CLASS ---------------------------------------------
class ElemDummy :  public Elem {

public:

  // Ensure that when "Virtual template" are not defined in derived classes
  // an error is thrown to avoid infinite loop in wrapper function
  
  template<int dim>
  void computeGalerkinTerm(FemEquationTerm *fet, SVec<double,3> &X, 
			   Vec<double> &d2wall, SVec<double,dim> &V, 
									SVec<double,dim> &R, 
									Vec<GhostPoint<dim>*> *ghostPoints=0, 
									LevelSetStructure *LSS=0) {
    fprintf(stderr, "Error: undefined function for this elem type\n"); exit(1);
  }

  template<int dim>
  void computeGalerkinTerm_e(FemEquationTerm *fet, SVec<double,3> &X, 
									  Vec<double> &d2wall, SVec<double,dim> &V, 
									  SVec<double,dim> &R, 
									  Vec<GhostPoint<dim>*> *ghostPoints=0, 
									  LevelSetStructure *LSS=0) {
    fprintf(stderr, "Error: undefined function for this elem type\n"); exit(1);
  }
  
  template<int dim>
  void computeVMSLESTerm(VMSLESTerm *vmst, SVec<double,dim> &VBar,
			 SVec<double,3> &X, SVec<double,dim> &V,
			 SVec<double,dim> &Sigma) {
    fprintf(stderr, "Error: undefined function for this elem type\n"); exit(1);
  }
  
  template<int dim>
  void computeMBarAndM(DynamicVMSTerm *dvmst, SVec<double,dim> **VBar,
		       SVec<double,1> **volRatio, SVec<double,3> &X,
		       SVec<double,dim> &V, SVec<double,dim> &MBar,
		       SVec<double,dim> &M) {
    fprintf(stderr, "Error: undefined function for this elem type\n"); exit(1);
  }
  
  template<int dim>
  void computeDynamicVMSTerm(DynamicVMSTerm *dvmst, SVec<double,dim> **VBar,
			     SVec<double,3> &X, SVec<double,dim> &V,
			     SVec<double,dim> &S, Vec<double> &CsDelSq,
			     Vec<double> &PrT, Vec<double> *Cs,
			     Vec<double> &Delta) {
    fprintf(stderr, "Error: undefined function for this elem type\n"); exit(1);
  }
  
  template<int dim>
  void computeSmagorinskyLESTerm(SmagorinskyLESTerm *smag, SVec<double,3> &X,
				 SVec<double,dim> &V, SVec<double,dim> &R,
                                 Vec<GhostPoint<dim>*> *ghostPoints=0, 
											LevelSetStructure *LSS=0) {
    fprintf(stderr, "Error: undefined function for this elem type\n"); exit(1);
  }
 
  template<int dim>
  void computeSmagorinskyLESTerm_e(SmagorinskyLESTerm *smag, SVec<double,3> &X,
											  SVec<double,dim> &V, SVec<double,dim> &R,
											  Vec<GhostPoint<dim>*> *ghostPoints=0, 
											  LevelSetStructure *LSS=0) {
    fprintf(stderr, "Error: undefined function for this elem type\n"); exit(1);
  }

  template<int dim>
  void computeWaleLESTerm(WaleLESTerm *wale, SVec<double,3> &X,
		          SVec<double,dim> &V, SVec<double,dim> &R,
                          Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0) {
    fprintf(stderr, "Error: undefined function for this elem type\n"); exit(1);
  }


  template<int dim>
  void computeWaleLESTerm_e(WaleLESTerm *wale, SVec<double,3> &X,
									 SVec<double,dim> &V, SVec<double,dim> &R,
									 Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0) {
    fprintf(stderr, "Error: undefined function for this elem type\n"); exit(1);
  }

  template<int dim>
  void computeDynamicLESTerm(DynamicLESTerm *dles, SVec<double,2> &Cs, 
			     SVec<double,3> &X, SVec<double,dim> &V, SVec<double,dim> &R,
                             Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0) {
    fprintf(stderr, "Error: undefined function for this elem type\n"); exit(1);
  }
  
  template<int dim>
  void computeDynamicLESTerm_e(DynamicLESTerm *dles, SVec<double,2> &Cs, 
										 SVec<double,3> &X, SVec<double,dim> &V, SVec<double,dim> &R,
										 Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0) {
	  fprintf(stderr, "Error: undefined function for this elem type\n"); exit(1);
  }
  
  template<int dim>
  void computeFaceGalerkinTerm(FemEquationTerm *fet, int face[3], int code, Vec3D &n, 
			       SVec<double,3> &X, Vec<double> &d2wall, double *Vwall, 
			       SVec<double,dim> &V, SVec<double,dim> &R) {
    fprintf(stderr, "Error: undefined function for this elem type\n"); exit(1);
  }
  

  template<int dim>
  void computeP1Avg(SVec<double,dim> &VCap, SVec<double,16> &Mom_Test, SVec<double,6> &Sij_Test, 
                    Vec<double> &modS_Test, SVec<double,8> &Eng_Test, SVec<double,3> &X, 
                    SVec<double,dim> &V, double gam, double R,
                    Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0) {
    fprintf(stderr, "Error: undefined function for this elem type\n"); exit(1);
  }

  template<int dim>
  void computeP1Avg_e(SVec<double,dim> &VCap, SVec<double,16> &Mom_Test, SVec<double,6> &Sij_Test, 
							 Vec<double> &modS_Test, SVec<double,8> &Eng_Test, SVec<double,3> &X, 
							 SVec<double,dim> &V, double gam, double R,
							 Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0) {
	  fprintf(stderr, "Error: undefined function for this elem type\n"); exit(1);
  }

  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(FemEquationTerm *, SVec<double,3> &, 
				   Vec<double> &, Vec<double> &, 
				   SVec<double,dim> &, GenMat<Scalar,neq> &,
                                   Vec<GhostPoint<dim>*> *gp=0,LevelSetStructure *LSS=0) {
    fprintf(stderr, "Error: undefined function for this elem type\n"); exit(1);
  }
  
  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm_e(FemEquationTerm *, SVec<double,3> &, 
												 Vec<double> &, Vec<double> &, 
												 SVec<double,dim> &, GenMat<Scalar,neq> &,
												 Vec<GhostPoint<dim>*> *gp=0,LevelSetStructure *LSS=0) {
	  fprintf(stderr, "Error: undefined function for this elem type\n"); exit(1);
  }
  
  template<int dim, class Scalar, int neq>
  void computeFaceJacobianGalerkinTerm(FemEquationTerm *fet, int face[3], int code, 
				       Vec3D &n, SVec<double,3> &X, Vec<double> &ctrlVol,
				       Vec<double> &d2wall, double *Vwall, 
				       SVec<double,dim> &V, GenMat<Scalar,neq> &A) {
    fprintf(stderr, "Error: undefined function for this elem type\n"); exit(1);
  }

// Included (MB)
  template<int dim>
  void computeDerivativeOfGalerkinTerm(FemEquationTerm *fet, SVec<double,3> &X, SVec<double,3> &dX,
			      Vec<double> &d2wall, SVec<double,dim> &V, SVec<double,dim> &dV, double dMach,
			      SVec<double,dim> &dR) {
    fprintf(stderr, "Error: undefined function (computeDerivativeOfGalerkinTerm) for this elem type\n"); exit(1);
  }

  template<int dim>
  void computeDerivativeOperatorsOfGalerkinTerm(FemEquationTerm *fet, SVec<double,3> &X,
			      Vec<double> &d2wall, SVec<double,dim> &V, RectangularSparseMat<double,3,dim> &dViscousFluxdX) { //YC
    fprintf(stderr, "Error: undefined function (computeDerivativeOperatorsOfGalerkinTerm) for this elem type\n"); exit(1);
  }


  template<int dim>
  void computeDerivativeOfFaceGalerkinTerm(FemEquationTerm *fet, int face[3], int code, Vec3D &n, Vec3D &dn,
				  SVec<double,3> &X, SVec<double,3> &dX, Vec<double> &d2wall, double *Vwall, double *dVwall,
				  SVec<double,dim> &V, SVec<double,dim> &dV, double dMach, SVec<double,dim> &dR) {
    fprintf(stderr, "Error: undefined function (computeDerivativeOfFaceGalerkinTerm) for this elem type\n"); exit(1);
  }
  
  // X is the deformed nodal location vector
  template<int dim> 
  int interpolateSolution(SVec<double,3>& X, SVec<double,dim>& U, const Vec3D& loc, double sol[dim], LevelSetStructure* LSS,
                          Vec<GhostPoint<dim>*>* ghostPoints, VarFcn* varFcn) {

    fprintf(stderr, "Error: undefined function (interpolateSolution0 for this elem type\n"); exit(1);
    return -1;
  }


//  // Included (YC)
//    template<int dim>
//    void computeDerivativeOperatorsOfGalerkinTerm(FemEquationTerm *, GeoState &, SVec<double,3> &,
//  			   SVec<double,dim> &, RectangularSparseMat<double,3,dim> &);

// Level Set Reinitialization

  template<int dimLS>
  void computeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                 SVec<double,dimLS> &ddx, SVec<double,dimLS> &ddy,
                                 SVec<double,dimLS> &ddz,
                                 SVec<double,dimLS> &Phi,SVec<double,1> &Psi){
    fprintf(stderr, "Error: undefined function (computeDistanceCloseNodes) for this elem type\n"); exit(1);
  }

  template<int dimLS>
  void recomputeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                 SVec<double,dimLS> &ddx, SVec<double,dimLS> &ddy,
                                 SVec<double,dimLS> &ddz,
                                 SVec<double,dimLS> &Phi,SVec<double,1> &Psi){
    fprintf(stderr, "Error: undefined function (recomputeDistanceCloseNodes) for this elem type\n"); exit(1);
  }

  template<int dimLS>
  void computeDistanceLevelNodes(int lsdim, Vec<int> &Tag, int level,
                                 SVec<double,3> &X, SVec<double,1> &Psi, SVec<double,dimLS> &Phi){
    fprintf(stderr, "Error: undefined function (computeDistanceLevelNodes) for this elem type\n"); exit(1);
  }

  template<int dim, class Obj>
    void integrateFunction(Obj* obj,SVec<double,3> &X,SVec<double,dim>& V, void (Obj::*F)(int node, const double* loc,double* f),int) {
    fprintf(stderr, "Error: undefined function (integrateFunction) for this elem type\n"); exit(1);
  }
  template<int dim>
  void FastMarchingDistanceUpdate(int node, Vec<int> &Tag, int level,
                              SVec<double,3> &X,SVec<double,dim> &d2wall)
  {
	fprintf(stderr, "Error: undefined function (FastMarchingDistanceUpdate) for this elem type\n"); exit(1);
  }
  
};

//------------------------------------------------------------------------------

class ElemSet {

  int numElems;
	int numSampledElems;

  Elem **elems;
  BlockAlloc memElems;
	bool sampleMesh;
	std::vector<int> elemsConnectedToSampleNode;	// for Gappy ROM

public:

  BlockAlloc& getBlockAllocator() { return memElems; }

  // Dummy constructor.  Sometimes we need to create a dummy elem set
  // object.
  ElemSet() : elems(NULL) { }

  ElemSet(int);
  ~ElemSet();

  Elem &operator[](int i) const { return *elems[i]; }

  Elem** getPointer() { return elems; }

  void addElem(int i, Elem *elem) { elems[i] = elem; }

  int size() const { return numElems; }

  int read(BinFileHandler&, int, int (*)[2], int *, map<int, VolumeData *> &);

  // compute Viscous Contribution Adam 2010.06.01
  template<int dim>
	  void computeTimeStep(FemEquationTerm *,SVec<double,3> &,SVec<double,dim> &,Vec<double> &, LevelSetStructure *LSS = 0);

  // Included (YC) //TODO UNSURE
  template<int dim>
  void computeDerivativeOperatorsOfGalerkinTerm(FemEquationTerm *, GeoState &, SVec<double,3> &,
                                      		   SVec<double,dim> &, RectangularSparseMat<double,3,dim> &);

  template<int dim>
  void computeGalerkinTerm(FemEquationTerm *, GeoState &, SVec<double,3> &, 
			   SVec<double,dim> &, SVec<double,dim> &,
									Vec<GhostPoint<dim>*> *ghostPoints=0, 
									LevelSetStructure *LSS=0, bool externalSI=false);

  template<int dim>
  void computeGalerkinTermRestrict(FemEquationTerm *, GeoState &, 
											  SVec<double,3> &, SVec<double,dim> &, 
											  SVec<double,dim> &, const std::vector<int> &,
											  Vec<GhostPoint<dim>*> *ghostPoints=0, 
											  LevelSetStructure *LSS=0);
  template<int dim>
  void computeVMSLESTerm(VMSLESTerm *, SVec<double,dim> &, SVec<double,3> &, 
			 SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeDynamicVMSTerm(DynamicVMSTerm *, SVec<double,dim> **, SVec<double,3> &,
                             SVec<double,dim> &, SVec<double,dim> &, Vec<double> &, Vec<double> &,
                             Vec<double> *, Vec<double> &);

  template<int dim>
  void computeMBarAndM(DynamicVMSTerm *, SVec<double,dim> **, SVec<double,1> **, SVec<double,3> &,
                       SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &);
  
  template<int dim>
  void computeSmagorinskyLESTerm(SmagorinskyLESTerm *, SVec<double,3> &, SVec<double,dim> &,
				 SVec<double,dim> &,
											Vec<GhostPoint<dim>*> *ghostPoints=0, 
											LevelSetStructure *LSS=0, bool externalSI = false);

  template<int dim>
  void computeWaleLESTerm(WaleLESTerm *, SVec<double,3> &, SVec<double,dim> &,
		          SVec<double,dim> &,
								  Vec<GhostPoint<dim>*> *ghostPoints=0, 
								  LevelSetStructure *LSS=0, bool externalSI = false);

  template<int dim>
  void computeDynamicLESTerm(DynamicLESTerm *, SVec<double,2> &,
                             SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &,
									  Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0,
									  bool externalSI = false);

  template<int dim>
  void computeTestFilterAvgs(SVec<double,dim> &, SVec<double,16> &, SVec<double,6> &,Vec<double> &,
                             SVec<double,8> &, SVec<double,3> &, SVec<double,dim> &, double, double,
									  Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0, bool externalSI=false);

  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(FemEquationTerm *fet, GeoState &geoState, 
				   SVec<double,3> &X, Vec<double> &ctrlVol,
				   SVec<double,dim> &V, GenMat<Scalar,neq> &A,
											  Vec<GhostPoint<dim>*>* ghostPoints=0, 
											  LevelSetStructure *LSS=0, bool externalSI=false);
    
// Included (MB)
  template<int dim>
  void computeDerivativeOfGalerkinTerm(FemEquationTerm *, GeoState &, SVec<double,3> &, SVec<double,3> &,
			   SVec<double,dim> &, SVec<double,dim> &, double, SVec<double,dim> &);

// Level Set Reinitialization
  template<int dimLS>
  void computeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                 SVec<double,dimLS> &ddx, SVec<double,dimLS> &ddy,
                                 SVec<double,dimLS> &ddz,
                                 SVec<double,dimLS> &Phi,SVec<double,1> &Psi);
  template<int dimLS>
  void recomputeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                 SVec<double,dimLS> &ddx, SVec<double,dimLS> &ddy,
                                 SVec<double,dimLS> &ddz,
                                 SVec<double,dimLS> &Phi,SVec<double,1> &Psi);
  template<int dimLS>
  void computeDistanceLevelNodes(int lsdim, Vec<int> &Tag, int level,
                                 SVec<double,3> &X, SVec<double,1> &Psi, SVec<double,dimLS> &Phi);

	void computeConnectedElems(const std::vector<int> &);

  template<int dim, class Obj>
  void integrateFunction(Obj* obj,SVec<double,3> &X,SVec<double,dim>& V, void (Obj::*F)(int node, const double* loc,double* f),int);

  // X is the deformed nodal location vector
  template<int dim> 
  void interpolateSolution(SVec<double,3>& X, SVec<double,dim>& U, const std::vector<Vec3D>& locs,
                           double (*sol)[dim], int* status,int* last, LevelSetStructure* LSS = 0,
                           Vec<GhostPoint<dim>*>* ghostPoints = 0, VarFcn* varFcn = 0,
			   bool assumeCache = false);
};

#ifdef TEMPLATE_FIX
#include <Elem.C>
#endif

#include <ElemTet.h>

#endif
