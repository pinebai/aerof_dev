#ifndef _FACE_H_
#define _FACE_H_

#include <PostFcn.h>
#include <MapFace.h>
#include <BlockAlloc.h>

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

class VarFcn;
class FluxFcn;
class FemEquationTerm;
class EdgeSet;
class ElemSet;
class GeoState;
class BinFileHandler;
class TimeLowMachPrec;
class LevelSetStructure;

struct Vec3D;

template<int dim> class BcData;
template<class Scalar> class Vec;
template<class Scalar, int dim> class GenMat;
template<class Scalar, int dim> class SVec;
template<class VecType> class VecSet;
template<class VecType, class VT2> class SubVecSet;
template<int dim> class ExactRiemannSolver;

class FaceTria;

#define NOT_CORRECTED(msg) {						        \
    static int first_time = 1;                                                  \
    if (type()!=Face::TRIA && first_time) {				        \
      fprintf(stderr, "---------------------------------------------------\n"); \
      fprintf(stderr, "  WARNING: Using non triangular faces\n");               \
      fprintf(stderr, "  %s:%d: Function %s not fully corrected yet\n\n",	\
	      __FILE__, __LINE__, __FUNCTION__);				\
      fprintf(stderr, "  %s\n", msg);                                           \
      fprintf(stderr, "---------------------------------------------------\n"); \
      first_time = 0;                                                           \
    }                                                                           \
  }


//-------------- GENERAL HELPERS -----------------------------------------------

struct HHCoeffs {
  double s0[3], s1[3];  //each face is split into three facets. 
  double currentDt;
};

class GenFaceHelper_dim {
public:
  virtual  void *forClassTria(FaceTria *, int size, char *memorySpace) = 0;
};

class GenFaceHelper_Scalar_dim_neq {
public:
  virtual  void *forClassTria(FaceTria *, int size, char *memorySpace) = 0;
};


//-------------- GENERAL WRAPPERS ----------------------------------------------
template<int dim>
class GenFaceWrapper_dim {
public:
  virtual void computeNodalForce(ElemSet &, PostFcn *, SVec<double,3> &, Vec<double> &, 
				 double *, SVec<double,dim> &, double, SVec<double,3> &, double* gradP[3]) = 0;
  virtual void computeNodalHeatPower(ElemSet &, PostFcn*, SVec<double,3>&, Vec<double>&, 
				     double*, SVec<double,dim>&, Vec<double>&) = 0;
  virtual double computeHeatFluxes(ElemSet &,PostFcn*, SVec<double,3>&, Vec<double>&,
                                    double*, SVec<double,dim>&) = 0;

  virtual void computeNodalHeatFluxRelatedValues(ElemSet &, PostFcn*, SVec<double,3>&, Vec<double>&, double*,
                                           SVec<double,dim>&, Vec<double>&, Vec<double>&, bool)=0;
  virtual void computeForceAndMoment(ElemSet &, PostFcn *, SVec<double,3> &, Vec<double> &, 
				     double *, SVec<double,dim> &, Vec3D &, Vec3D &, Vec3D &, 
				     Vec3D &, Vec3D &,  double* gradP[3], int, 
                                     SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX, Vec<double> *genCF) = 0;
  virtual void computeForceAndMoment(ExactRiemannSolver<dim>&, VarFcn*, Vec<Vec3D> &, Vec<double> &,
                                     ElemSet &, PostFcn *, SVec<double,3> &, Vec<double> &,
                                     double *, SVec<double,dim> &, Vec3D &, Vec3D &, Vec3D &,
                                     Vec3D &, Vec3D &,  double* gradP[3], int,
                                     SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX, Vec<double> *genCF) = 0;
  virtual double computeInterfaceWork(ElemSet &, PostFcn*, SVec<double,3>&, Vec<double>&, 
				      double, double*, SVec<double,dim>&, double) = 0;
  virtual void computeScalarQuantity(PostFcn::ScalarType, ElemSet &, PostFcn *, SVec<double,3> &, 
				     Vec<double> &, double *, SVec<double,dim> &, SVec<double,2> &) = 0;
  virtual void computeGalerkinTerm(ElemSet &, FemEquationTerm *, SVec<double,3> &, 
				   Vec<double> &, double *, SVec<double,dim> &, SVec<double,dim> &, LevelSetStructure *LSS=0) = 0;
  virtual void computeForceDerivs(ElemSet &, VarFcn *, SVec<double,3> &, 
				  SVec<double,dim> &,SVec<double,dim> &, 
				  Vec<double> &, SVec<double,3> **) = 0;
  virtual void computeForceCoefficients(PostFcn *, Vec3D &, ElemSet &, SVec<double,3> &, 
					SVec<double,dim> &, Vec<double> &, 
					SVec<double, dim> &,  double, Vec3D &, Vec3D &, 
					Vec3D &, Vec3D &, double* gradP[3], VecSet< SVec<double,3> > *mX = 0,
                                        Vec<double> *genCF = 0) = 0;
  virtual void computeFDerivs(ElemSet &, VarFcn *, SVec<double,3> &, 
			      SVec<double,dim> &, Vec3D (*)) = 0;


// Included (MB)
  virtual void computeDerivativeOfNodalForce(ElemSet &, PostFcn *, SVec<double,3> &, SVec<double,3> &,
                                             Vec<double> &, double *, double *,
                                             SVec<double,dim> &, SVec<double,dim> &,
                                             double, double [3], SVec<double,3> &, 
                                             double* gradP[3], double* dGradP[3]) = 0;

  virtual void computeDerivativeOperatorsOfNodalForce(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X, SVec<double,dim> &V, double Pin, double* gradP[3],
                                                      RectangularSparseMat<double,3,3> &dForcedX,
                                                      RectangularSparseMat<double,3,3> &dForcedGradP,
                                                      RectangularSparseMat<double,dim,3> &dForcedV,
                                                      RectangularSparseMat<double,3,3> &dForcedS) = 0;

  virtual void computeDerivativeOfNodalHeatPower(ElemSet&, PostFcn*, SVec<double,3>&, 
                                                 SVec<double,3>&, Vec<double>&, 
			                         double*, double*, SVec<double,dim>&, 
                                                 SVec<double,dim>&, double [3], Vec<double>&) = 0;

  virtual void computeDerivativeOfForceAndMoment(ElemSet &, PostFcn *, SVec<double,3> &, SVec<double,3> &,
                                                 Vec<double> &, double *, double *,
                                                 SVec<double,dim> &, SVec<double,dim> &, double [3], 
                                                 Vec3D &, Vec3D &, Vec3D &, Vec3D &, Vec3D &, 
                                                 double* gradP[3], double* dGradP[3], int = 0) = 0;
  /*
    virtual void computeDerivativeOfForceAndMoment2(ElemSet &, PostFcn *, SVec<double,3> &, SVec<double,3> &,
                                                    Vec<double> &, double *, double *,
                                                    SVec<double,dim> &, SVec<double,dim> &, double [3],
                                                    Vec3D &, Vec3D &, Vec3D &, Vec3D &, Vec3D &,
                                                    double* gradP[3], double* dGradP[3], int = 0) = 0;
  */
  virtual void computeDerivativeOperatorsOfForceAndMoment(ElemSet &, PostFcn *, SVec<double,3> &,
                                                          Vec<double> &, double *, SVec<double,dim> &, Vec3D &,
                                                          double* gradP[3], int ,
                                                          RectangularSparseMat<double,3,3> &,
                                                          RectangularSparseMat<double,3,3> &,
                                                          RectangularSparseMat<double,dim,3> &,
                                                          RectangularSparseMat<double,3,3> &,
                                                          RectangularSparseMat<double,dim,3> &,
                                                          RectangularSparseMat<double,3,3> &,
                                                          RectangularSparseMat<double,3,3> &,
                                                          RectangularSparseMat<double,3,3> &,
                                                          RectangularSparseMat<double,dim,3> &,
                                                          RectangularSparseMat<double,3,3> &,
                                                          RectangularSparseMat<double,3,3> &,
                                                          RectangularSparseMat<double,dim,3> &) = 0;

  virtual void computeDerivativeOfGalerkinTerm(ElemSet &, FemEquationTerm *, SVec<double,3> &, SVec<double,3> &,
                                               Vec<double> &, double *, double *, SVec<double,dim> &, 
                                               SVec<double,dim> &, double, SVec<double,dim> &) = 0;
  
  virtual void computeDerivativeOfGalerkinTermEmb(ElemSet &, FemEquationTerm *, SVec<double,3> &, SVec<double,3> &,
                                               Vec<double> &, double *, double *, SVec<double,dim> &,
                                               SVec<double,dim> &, double, SVec<double,dim> &,LevelSetStructure*) = 0;

  virtual void computeBCsJacobianWallValues(ElemSet &, FemEquationTerm *, SVec<double,3> &, 
                                            Vec<double> &, double *, double *, SVec<double,dim> &) = 0;

};

template<class Scalar, int dim, int neq>
class GenFaceWrapper_Scalar_dim_neq {
public:
  virtual void computeJacobianGalerkinTerm(ElemSet &, FemEquationTerm *, SVec<double,3> &, 
                 Vec<double> &, Vec<double> &, double *,
                 SVec<double,dim> &, GenMat<Scalar,neq> &) = 0;
};


//-------------- REAL WRAPPERS -------------------------------------------------
template<class Target, int dim>
class  FaceWrapper_dim : public GenFaceWrapper_dim<dim> {
  Target *t;

public:
  FaceWrapper_dim(Target *tt) : t(tt) { };

  void computeNodalForce(ElemSet &elems,
			 PostFcn *postFcn, SVec<double,3> &X, 
			 Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V,
			 double pin, SVec<double,3> &F, double* gradP[3]) {
    t->computeNodalForce(elems, postFcn, X, d2wall, Vwall, V,
			 pin, F, gradP);
  }

  void computeNodalHeatPower(ElemSet &elems,
			     PostFcn* postFcn, SVec<double,3>& X, 
			     Vec<double>& d2wall, double* Vwall, 
			     SVec<double,dim>& V, Vec<double>& P) {
    t->computeNodalHeatPower(elems, postFcn, X, 
			     d2wall, Vwall, V, P);
  }

  double computeHeatFluxes(ElemSet &elems,
                             PostFcn* postFcn, SVec<double,3>& X,
                             Vec<double>& d2wall, double* Vwall,
                             SVec<double,dim>& V) {
    return t->computeHeatFluxes(elems, postFcn, X,
                             d2wall, Vwall, V);
  }



  void computeNodalHeatFluxRelatedValues(ElemSet &elems, PostFcn* postFcn, SVec<double,3>& X, 
                                 Vec<double>& d2wall, double* Vwall, SVec<double,dim>& V, Vec<double>& P, Vec<double>& N, bool includeKappa){
     t->computeNodalHeatFluxRelatedValues(elems, postFcn, X,
                             d2wall, Vwall, V, P, N, includeKappa);
  }

  void computeForceAndMoment(ElemSet &elems,
			     PostFcn *postFcn, SVec<double,3> &X, 
			     Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V, 
			     Vec3D &x0, Vec3D &Fi, Vec3D &Mi, Vec3D &Fv, Vec3D &Mv, 
			      double* gradP[3], int hydro, SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX,
                                        Vec<double> *genCF) {
    t->computeForceAndMoment(elems, postFcn, X,  d2wall, Vwall, V, 
			     x0, Fi, Mi, Fv, Mv, gradP, hydro, mX, genCF);
  }

  void computeForceAndMoment(ExactRiemannSolver<dim> &riemann,
                             VarFcn *varFcn, Vec<Vec3D> &n, Vec<double> &nVel, ElemSet &elems,
                             PostFcn *postFcn, SVec<double,3> &X,
                             Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V,
                             Vec3D &x0, Vec3D &Fi, Vec3D &Mi, Vec3D &Fv, Vec3D &Mv,
                              double* gradP[3], int hydro, SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX,
                                        Vec<double> *genCF) {
    t->computeForceAndMoment(riemann, varFcn, n, nVel, elems, postFcn, X,  d2wall, Vwall, V,
                             x0, Fi, Mi, Fv, Mv, gradP, hydro, mX, genCF);
  }
  
  double computeInterfaceWork(ElemSet &elems, PostFcn* postFcn, 
			      SVec<double,3>& X, Vec<double>& d2wall, double ndot, 
			      double* Vwall, SVec<double,dim>& V, double pin) {
    return t->computeInterfaceWork(elems, postFcn, X, d2wall, ndot, Vwall, V, pin);
  }

  void computeScalarQuantity(PostFcn::ScalarType type, ElemSet &elems, PostFcn *postFcn, 
			     SVec<double,3> &X, Vec<double> &d2wall, double *Vwall, 
			     SVec<double,dim> &V, SVec<double,2> &Q) {
    t->computeScalarQuantity(type, elems, postFcn, X, d2wall, Vwall, V, Q);
  }

  void computeGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, 
			   Vec<double> &d2wall, double *Vwall,
			   SVec<double,dim> &V, SVec<double,dim> &R,LevelSetStructure *LSS=0) {
    t->computeGalerkinTerm(elems, fet, X, d2wall, Vwall, V, R, LSS);
  }
  
  void computeForceDerivs(ElemSet &elems, VarFcn *varFcn, SVec<double,3> &X, 
			  SVec<double,dim> &V, SVec<double,dim> &deltaU, Vec<double> &modalF, 
			  SVec<double,3> **localMX) {
    t->computeForceDerivs(elems, varFcn, X, V, deltaU, modalF, localMX);
  }

  void computeForceCoefficients(PostFcn *postFcn, Vec3D &x0, ElemSet &elems, 
				SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &d2wall, 
				SVec<double, dim> &Vwall, double pInfty, Vec3D &CFi, Vec3D &CMi, 
				Vec3D &CFv, Vec3D &CMv, double* gradP[3], VecSet< SVec<double,3> > *mX = 0,
                                        Vec<double> *genCF = 0) {
    t->computeForceCoefficients(postFcn, x0, elems, X, V, d2wall, Vwall, pInfty, 
				CFi, CMi, CFv, CMv, gradP, mX, genCF);
  }

  void computeFDerivs(ElemSet &elems, VarFcn *varFcn, SVec<double,3> &X, 
		      SVec<double,dim> &Vgl, Vec3D (*F)) {
    t->computeFDerivs(elems, varFcn, X, Vgl, F);
  }

// Included (MB)
  void computeDerivativeOfNodalForce(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X, SVec<double,3> &dX,
			     Vec<double> &d2wall, double *Vwall, double *dVwall, SVec<double,dim> &V, SVec<double,dim> &dV,
			     double pin, double dS[3], SVec<double,3> &dF, double* gradP[3], double* dGradP[3]) {
    t->computeDerivativeOfNodalForce(elems, postFcn, X, dX, d2wall, Vwall, dVwall, V, dV, pin, dS, dF, gradP, dGradP);
  }

  void computeDerivativeOperatorsOfNodalForce(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X, SVec<double,dim> &V, double Pin, double* gradP[3],
                                              RectangularSparseMat<double,3,3> &dForcedX,
                                              RectangularSparseMat<double,3,3> &dForcedGradP,
                                              RectangularSparseMat<double,dim,3> &dForcedV,
                                              RectangularSparseMat<double,3,3> &dForcedS) {
	t->computeDerivativeOperatorsOfNodalForce(elems, postFcn, X, V, Pin, gradP, dForcedX, dForcedGradP, dForcedV, dForcedS);
  }

  void computeDerivativeOfNodalHeatPower(ElemSet& elems, PostFcn* postFcn, SVec<double,3>& X, 
                                         SVec<double,3>& dX, Vec<double>& d2wall, double* Vwall, 
                                         double* dVwall, SVec<double,dim>& V, SVec<double,dim>& dV, 
                                         double dS[3], Vec<double>& dP) {
    t->computeDerivativeOfNodalHeatPower(elems, postFcn, X, dX, d2wall, Vwall, dVwall, V, dV, dS, dP);
  }

  void computeDerivativeOfForceAndMoment(ElemSet &elems, PostFcn *postFcn,
                                         SVec<double,3> &X, SVec<double,3> &dX,
                                         Vec<double> &d2wall, double *Vwall, double *dVwall,
                                         SVec<double,dim> &V, SVec<double,dim> &dV, double dS[3],
                                         Vec3D &x0, Vec3D &dFi, Vec3D &dMi, Vec3D &dFv, Vec3D &dMv, 
                                         double* gradP[3], double* dGradP[3], int hydro) {
    t->computeDerivativeOfForceAndMoment(elems, postFcn, X, dX, d2wall, Vwall, dVwall, V, dV, dS, x0, dFi, dMi, dFv, dMv, gradP, dGradP, hydro);
  }
  /*
    void computeDerivativeOfForceAndMoment2(ElemSet &elems, PostFcn *postFcn,
                                           SVec<double,3> &X, SVec<double,3> &dX,
                                           Vec<double> &d2wall, double *Vwall, double *dVwall,
                                           SVec<double,dim> &V, SVec<double,dim> &dV, double dS[3],
                                           Vec3D &x0, Vec3D &dFi, Vec3D &dMi, Vec3D &dFv, Vec3D &dMv,
                                           double* gradP[3], double* dGradP[3], int hydro) {
      t->computeDerivativeOfForceAndMoment2(elems, postFcn, X, dX, d2wall, Vwall, dVwall, V, dV, dS, x0, dFi, dMi, dFv, dMv, gradP, dGradP, hydro);
    }
  */
  void computeDerivativeOperatorsOfForceAndMoment(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X,
                                                  Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V, Vec3D &x0,
                                                  double* gradP[3], int hydro,
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
                                                  RectangularSparseMat<double,dim,3> &dMvdV) {
    t->computeDerivativeOperatorsOfForceAndMoment(elems, postFcn, X, d2wall, Vwall, V, x0, gradP, hydro,
                                                  dFidGradP, dFidX, dFidV, dFvdX, dFvdV, dFidS,
                                                  dMidGradP, dMidX, dMidV, dMidS, dMvdX, dMvdV);
  }

  void computeDerivativeOfGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, SVec<double,3> &dX,
                                       Vec<double> &d2wall, double *Vwall, double *dVwall, SVec<double,dim> &V, 
                                       SVec<double,dim> &dV, double dMach, SVec<double,dim> &dR) {
    t->computeDerivativeOfGalerkinTerm(elems, fet, X, dX, d2wall, Vwall, dVwall, V, dV, dMach, dR);
  }


  void computeDerivativeOfGalerkinTermEmb(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, SVec<double,3> &dX,
                                       Vec<double> &d2wall, double *Vwall, double *dVwall, SVec<double,dim> &V,
                                       SVec<double,dim> &dV, double dMach, SVec<double,dim> &dR,LevelSetStructure *LSS) {
    t->computeDerivativeOfGalerkinTermEmb(elems, fet, X, dX, d2wall, Vwall, dVwall, V, dV, dMach, dR, LSS);
  }

  void computeBCsJacobianWallValues(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, Vec<double> &d2wall, 
                                    double *Vwall, double *dVwall, SVec<double,dim> &V) {
    t->computeBCsJacobianWallValues(elems, fet, X, d2wall, Vwall, dVwall, V);
  }

};

template<class Target, class Scalar, int dim, int neq>
class  FaceWrapper_Scalar_dim_neq : public 
GenFaceWrapper_Scalar_dim_neq<Scalar,dim,neq> {

  Target *t;

public:
  FaceWrapper_Scalar_dim_neq(Target *tt) : t(tt) { };

  void computeJacobianGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, 
				   SVec<double,3> &X, Vec<double> &ctrlVol,
				   Vec<double> &d2wall, double *Vwall, 
				   SVec<double,dim> &V, GenMat<Scalar,neq> &A) {
    t->computeJacobianGalerkinTerm(elems, fet, X, ctrlVol, d2wall, Vwall, V, A);
  };
};


//-------------- REAL HELPERS --------------------------------------------------
template<int dim>
class FaceHelper_dim : public GenFaceHelper_dim {
public:

  void *forClassTria(FaceTria *tface, int size, char *memorySpace) {
    if(size < sizeof(FaceWrapper_dim<FaceTria, dim>) ) {
      fprintf(stderr, "Error: programming error in FaceHelper");
      exit(1);
    }
    return new (memorySpace) FaceWrapper_dim<FaceTria, dim>(tface);
  }

};


template<class Scalar, int dim, int neq>
class FaceHelper_Scalar_dim_neq : public GenFaceHelper_Scalar_dim_neq {
public:

  void *forClassTria(FaceTria *tface, int size, char *memorySpace) {
    if(size < sizeof(FaceWrapper_Scalar_dim_neq<FaceTria, Scalar, dim, neq>) ) {
      fprintf(stderr, "Error: programming error in FaceHelper");
      exit(1);
    }
    return new (memorySpace) FaceWrapper_Scalar_dim_neq<FaceTria, Scalar, dim, neq>(tface);
  }
  
};


//--------------- BASE FACE CLASS ----------------------------------------------
class Face {
public:
  enum Type {TRIA=4};

  static const int MaxNumNd = 4;

  Face();
  virtual int nodeNum(int i) const = 0;

protected:
  virtual void *getWrapper_dim(GenFaceHelper_dim *, 
			       int size, char *memorySpace) = 0;
  virtual void *getWrapper_Scalar_dim_neq(GenFaceHelper_Scalar_dim_neq *, 
					  int size, char *memorySpace) = 0;
  
  int code; //boundary condition code
  int elemNum;
  int surface_id;
  int normNum;

  //double ffWeight; // allows user to weight the farfield fluxes (for nonlinear ROM simulations)

  //for farfield flux
  Vec3D faceCenter;

  class HigherOrderMultiFluid* higherOrderMF;

  virtual int* nodeNum() = 0;
  virtual int& nodeNum(int i) = 0;
  virtual int& edgeNum(int i) = 0;
  virtual int  edgeEnd(int i, int k) = 0;  

public:

  int getEdgeNum(int i) { return edgeNum(i); }

  virtual void setEdgeNum(int edge_id, int l) = 0;

  void attachHigherOrderMF(class HigherOrderMultiFluid* mf) { higherOrderMF = mf; }

  // Number of nodes
  virtual int numNodes() = 0;

  // Number of normals to be stored
  virtual int numNorms() = 0;

  // Get element type
  virtual Type type() = 0;

  int operator[](int i) { return nodeNum(i); }  
  operator int *() { return nodeNum(); }
  operator MaxFace() { return MaxFace(numNodes(), nodeNum()); }

  /* WARNING : IS THIS THE RIGHT DEFINITION, WHEN NUMNODES()!=F.NUMNODES() ??? */
  /*           removed const ... function is actualy never used */
  bool operator<(Face &f)  {
    if(numNodes() != f.numNodes()) return numNodes() < f.numNodes();
    for (int i=0; i<numNodes(); i++) {
      if (nodeNum(i) < f.nodeNum(i)) return true;
      if (nodeNum(i) > f.nodeNum(i)) return false;
    }
    
    return false;
  }

  int* nodes(int* d = NULL) {  
    if (d)
      for (int i=0; i<numNodes(); i++) d[i] = nodeNum(i); 
    return nodeNum();
  }
  
  int getCode()  { return code; }
  int getSurfaceID() { return surface_id; }
  int getElementNumber() const { return elemNum; }

  //void setFarFieldBCWeight(double weight) {ffWeight = weight;};
  void setup(int, int *, int, int surface_id = 0);
  void setType(int *);
  void setType(int t) { code = t; }
  void setNodeType(int*, int*);
  void setNodeFaceType(int*);
  void setElementNumber(int elemNum, int rotDir);
  void tagNodesOnBoundaries(Vec<bool> &);
  void tagEdgesOnBoundaries(Vec<bool> &);
  void reorder();
  void numberEdges(EdgeSet &);
  void computeEdgeNormals(SVec<double,3>&, int*, SVec<double,6>&);

  // WARNING: THIS IS A FUNCTION FOR TFACES ONLY
  static double computeVolume(Vec3D &xa_n, Vec3D &xb_n, Vec3D &xc_n, 
			      Vec3D &xa_np1, Vec3D &xb_np1, Vec3D &xc_np1);
  
  // Compute total normal and return in Vec3D  
  virtual void computeNormal(SVec<double,3> &, Vec3D &) = 0;

  // Compute subface normals and save then in Vec<Vec3D>
  virtual void computeNormal(SVec<double,3> &, Vec<Vec3D> &) = 0;
  virtual void computeNormalConfig(SVec<double,3> &, SVec<double,3> &,
                                   Vec<Vec3D> &, Vec<double> &) = 0;
  virtual void computeNormalGCL1(SVec<double,3> &, SVec<double,3> &, 
				 SVec<double,3> &, Vec<Vec3D> &, Vec<double> &) = 0;
  virtual void computeNormalEZGCL1(double, SVec<double,3> &, SVec<double,3> &, 
				   Vec<Vec3D> &, Vec<double> &) = 0;

  // Get face total normal from Vec<Vec3D>
  virtual Vec3D getNormal(Vec<Vec3D> &) = 0;

  // Get subface i normal from Vec<Vec3D>
  virtual Vec3D getNormal(Vec<Vec3D> &, int) = 0;

  // Get face total normal velocity from Vec<double>
  virtual double getNormalVel(Vec<double> &) = 0;

  // Get subface i normal velocity from Vec<double>
  virtual double getNormalVel(Vec<double> &, int) = 0;

  template<class NodeMap>
  void renumberNodes(NodeMap &nodemap);

  template<int dim>
  void assignFreeStreamValues2(SVec<double,dim> &Uin, SVec<double,dim> &Uout, double *U);

  template<int dim>
  void assignFreeStreamValues(double *Uin, double *Uout, double *U);
  
  template<int dim>
  void assignPorousWallValues(SVec<double,dim> &Uin, double *U);
  
  template<int dim>
  void computeFaceBcValue(SVec<double,dim> &Unode, double *Uface);

  template<int dim1, int dim2>
  void computeNodeBcValue(SVec<double,3> &X, double *Uface, SVec<double,dim2> &Unode);
    
  // "Virtual template" functions (implemented through wrapper and helper functions)

  template<int dim>
  void computeNodalForce(ElemSet &elems,
			 PostFcn *postFcn, SVec<double,3> &X, 
			 Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V,
			 double pin, SVec<double,3> &F, double* gradP[3]) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeNodalForce(elems, postFcn, X, d2wall, Vwall, V,
			       pin, F, gradP);
  }

  template<int dim>
  void computeNodalHeatPower(ElemSet &elems,
			     PostFcn* postFcn, SVec<double,3>& X, 
			     Vec<double>& d2wall, double* Vwall, 
			     SVec<double,dim>& V, Vec<double>& P) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeNodalHeatPower(elems, postFcn, X, 
			     d2wall, Vwall, V, P);
  }

  template<int dim>
  double computeHeatFluxes(ElemSet &elems,
                             PostFcn* postFcn, SVec<double,3>& X,
                             Vec<double>& d2wall, double* Vwall,
                             SVec<double,dim>& V) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    return wrapper->computeHeatFluxes(elems, postFcn, X,
                             d2wall, Vwall, V);
  }

  template<int dim>
  void computeNodalHeatFluxRelatedValues(ElemSet &elems, PostFcn* postFcn, SVec<double,3>& X,
                                 Vec<double>& d2wall, double* Vwall, SVec<double,dim>& V, Vec<double>& P, Vec<double>& N, bool includeKappa){
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeNodalHeatFluxRelatedValues(elems, postFcn, X,
                             d2wall, Vwall, V, P, N, includeKappa);
  }

  template<int dim>
  void computeForceAndMoment(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X, 
			     Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V, 
			     Vec3D &x0, Vec3D &Fi, Vec3D &Mi, Vec3D &Fv, Vec3D &Mv, 
			     double* gradP[3], int hydro, SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX,
                                        Vec<double> *genCF) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeForceAndMoment(elems, postFcn, X,  d2wall, Vwall, V, 
			     x0, Fi, Mi, Fv, Mv, gradP, hydro, mX, genCF);
  }

  template<int dim>
  void computeForceAndMoment(ExactRiemannSolver<dim>& riemann, 
                             VarFcn *varFcn, Vec<Vec3D> &n, Vec<double> &nVel,
                             ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X,
                             Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V,
                             Vec3D &x0, Vec3D &Fi, Vec3D &Mi, Vec3D &Fv, Vec3D &Mv,
                             double* gradP[3], int hydro, SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX,
                                        Vec<double> *genCF) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeForceAndMoment(riemann, varFcn, n, nVel, elems, postFcn, X,  d2wall, Vwall, V,
                             x0, Fi, Mi, Fv, Mv, gradP, hydro, mX, genCF);
  }
 
  template<int dim>
  double computeInterfaceWork(ElemSet &elems, PostFcn* postFcn, 
			      SVec<double,3>& X, Vec<double>& d2wall, double ndot, 
			      double* Vwall, SVec<double,dim>& V, double pin) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    return wrapper->computeInterfaceWork(elems, postFcn, X, d2wall, ndot, Vwall, V, pin);
  }

  template<int dim>
  void computeScalarQuantity(PostFcn::ScalarType type, ElemSet &elems, PostFcn *postFcn, 
			     SVec<double,3> &X, Vec<double> &d2wall, double *Vwall, 
			     SVec<double,dim> &V, SVec<double,2> &Q) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeScalarQuantity(type, elems, postFcn, X, d2wall, Vwall, V, Q);
  }

  template<int dim>
  void computeTimeStep(VarFcn *varFcn, Vec<Vec3D> &normal, Vec<double> &normalVel,
		       SVec<double,dim> &V, Vec<double> &dt, 
		       TimeLowMachPrec &tprec, Vec<int> &fluidId);

  template<int dim>
  void computeTimeStep(FemEquationTerm *fet, VarFcn *varFcn, 
		       Vec<Vec3D> &normal, Vec<double> &normalVel,
		       SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &idti, 
							  Vec<double> &idtv, TimeLowMachPrec &tprec, LevelSetStructure *LSS=0);

  template<int dim>
  void computeFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normal, 
			       Vec<double> &normalVel, SVec<double,dim> &V, 
			       double *Ub, SVec<double,dim> &fluxes, double UbHH = 0);

  template<int dim>
  void computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                               FluxFcn **fluxFcn, Vec<Vec3D> &normal,
                               Vec<double> &normalVel, SVec<double,dim> &V,
                               double *Ub, SVec<double,dim> &fluxes,double UbHH = 0);

  template<int dim>
  void computeFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normal,
			       Vec<double> &normalVel, SVec<double,dim> &V,
			       double *Ub, Vec<int> &fluidId, 
			       SVec<double,dim> &fluxes,
                               LevelSetStructure* = 0,double UbHH = 0);

  template<int dim, int dimLS>
  void computeFiniteVolumeTermLS(FluxFcn **fluxFcn, Vec<Vec3D> &normal,
				 Vec<double> &normalVel, SVec<double,dim> &V,
				 SVec<double,dimLS> &Phi, SVec<double,dimLS> &PhiF);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normal, 
				       Vec<double> &normalVel, SVec<double,dim> &V, 
				       double *Ub, GenMat<Scalar,neq> &A,double UbHH = 0
                                    );

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim> &riemann,
                                       FluxFcn **fluxFcn, Vec<Vec3D> &normal, 
				       Vec<double> &normalVel, SVec<double,dim> &V, 
				       double *Ub, GenMat<Scalar,neq> &A,double UbHH = 0);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normal,
				       Vec<double> &normalVel, SVec<double,dim> &V,
				       double *Ub, GenMat<Scalar,neq> &A, int* nodeType,double UbHH = 0);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normal, 
				       Vec<double> &normalVel, SVec<double,dim> &V, 
				       double *Ub, GenMat<Scalar,neq> &A, Vec<int> &fluidId,
                                       LevelSetStructure* LSS,double UbHH = 0);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normal,
				       Vec<double> &normalVel, SVec<double,dim> &V,
				       double *Ub, GenMat<Scalar,neq> &A, 
                                       Vec<int> &fluidId, int* nodeType,double UbHH = 0);

  template<int dim, class Scalar, int dimLS>
  void computeJacobianFiniteVolumeTermLS(Vec<Vec3D> &normal,
					 Vec<double> &normalVel, SVec<double,dim> &V,
					 GenMat<Scalar,dimLS> &A);

  template<int dim>
  void computeGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, 
			   Vec<double> &d2wall, double *Vwall,
			   SVec<double,dim> &V, SVec<double,dim> &R, LevelSetStructure *LSS=0) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeGalerkinTerm(elems, fet, X, d2wall, Vwall, V, R, LSS);
  }
  
  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, 
				   SVec<double,3> &X, Vec<double> &ctrlVol,
				   Vec<double> &d2wall, double *Vwall, 
				   SVec<double,dim> &V, GenMat<Scalar,neq> &A) {
    FaceHelper_Scalar_dim_neq<Scalar,dim,neq> h;
    char xx[64];
    GenFaceWrapper_Scalar_dim_neq<Scalar,dim,neq> *wrapper=
      (GenFaceWrapper_Scalar_dim_neq<Scalar,dim,neq> *)getWrapper_Scalar_dim_neq(&h, 64, xx);
    wrapper->computeJacobianGalerkinTerm(elems, fet, X, ctrlVol, d2wall, Vwall, V, A);
  }

  template<int dim>
  void computeForceDerivs(ElemSet &elems, VarFcn *varFcn, SVec<double,3> &X, 
			  SVec<double,dim> &V, SVec<double,dim> &deltaU, Vec<double> &modalF, 
			  SVec<double,3> **localMX) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeForceDerivs(elems, varFcn, X, V, deltaU, modalF, localMX);
  }

  template<int dim>
  void computeForceCoefficients(PostFcn *postFcn, Vec3D &x0, ElemSet &elems, 
				SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &d2wall, 
				SVec<double, dim> &Vwall, double pInfty, Vec3D &CFi, Vec3D &CMi, 
				Vec3D &CFv, Vec3D &CMv, double* gradP[3]) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeForceCoefficients(postFcn, x0, elems, X, V, d2wall, Vwall, pInfty, 
				CFi, CMi, CFv, CMv, gradP);
  }

  template<int dim>
  void computeFDerivs(ElemSet &elems,
		      VarFcn *varFcn, SVec<double,3> &X, 
		      SVec<double,dim> &Vgl, Vec3D (*F)) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeFDerivs(elems, varFcn, X, Vgl, F);
  }
 
// Included (MB)
  virtual void computeNormalAndDerivative(SVec<double,3> &, SVec<double,3> &, Vec3D &, Vec3D&) = 0;

  virtual void computeDerivativeOfNormal(SVec<double,3> &, SVec<double,3> &, Vec3D &, Vec3D &, double &, double &) = 0;
  virtual void computeDerivativeOperatorsOfNormal(int, SVec<double,3> &, RectangularSparseMat<double,3,3> &) = 0;

  // Get face total normal derivative from Vec<Vec3D>
  virtual Vec3D getdNormal(Vec<Vec3D> &) = 0;

  // Get subface i normal derivative from Vec<Vec3D>
  virtual Vec3D getdNormal(Vec<Vec3D> &, int) = 0;

  // Get face total normal velocity derivative from Vec<double>
  virtual double getdNormalVel(Vec<double> &) = 0;

  // Get subface i normal velocity derivative from Vec<double>
  virtual double getdNormalVel(Vec<double> &, int) = 0;

  template<int dim>
  void computeDerivativeOfFiniteVolumeTerm(FluxFcn **fluxFcn, Vec<Vec3D> &normals,
				      Vec<Vec3D> &dNormals, Vec<double> normalVel, Vec<double> dNormalVel,
				      SVec<double,dim> &V, double *Ub,
				      double *dUb, SVec<double,dim> &dFluxes);

  template<int dim>
  void computeDerivativeOperatorsOfFiniteVolumeTerm(
              int, FluxFcn **fluxFcn, Vec<Vec3D> &normals,
				      Vec<double> normalVel, SVec<double,dim> &V, double *Ub,
				      RectangularSparseMat<double,3,dim> &dFluxdFaceNormal,
              RectangularSparseMat<double,1,dim> &dFluxdFaceNormalVel,
              RectangularSparseMat<double,dim,dim> &dFluxdUb);

  template<int dim1, int dim2>
  void computeDerivativeOfNodeBcValue(SVec<double,3> &X, SVec<double,3> &dX, double *Uface, double *dUface, SVec<double,dim2> &dUnode);

  template<int dim>
  void computeNodeBCsWallValues(SVec<double,3> &X, SVec<double,1> &dNormSA, double *dUfaceSA, SVec<double,dim> &dUnodeSA);

  template<int dim>
  void computeDerivativeOfTimeStep(FemEquationTerm *fet, VarFcn *varFcn, Vec<Vec3D>  &normals, Vec<Vec3D>  &dNormals, Vec<double> normalVel, Vec<double> dNormalVel,
			   SVec<double,3> &X, SVec<double,3> &dX, SVec<double,dim> &V, SVec<double,dim> &dV, 
			   Vec<double> &dIdti, Vec<double> &dIdtv, double dMach, 
                           TimeLowMachPrec &tprec);

  template<int dim>
  void computeDerivativeOfNodalForce(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X, SVec<double,3> &dX,
			     Vec<double> &d2wall, double *Vwall, double *dVwall, SVec<double,dim> &V, SVec<double,dim> &dV,
			     double pin, double dS[3], SVec<double,3> &dF, double* gradP[3],  double* dGradP[3]) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDerivativeOfNodalForce(elems, postFcn, X, dX, d2wall, Vwall, dVwall, V, dV, pin, dS, dF, gradP, dGradP);
  }

  template<int dim>
  void computeDerivativeOperatorsOfNodalForce(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X, SVec<double,dim> &V, double Pin, double *gradP[3],
                                              RectangularSparseMat<double,3,3> &dForcedX,
                                              RectangularSparseMat<double,3,3> &dForcedGradP,
                                              RectangularSparseMat<double,dim,3> &dForcedV,
                                              RectangularSparseMat<double,3,3> &dForcedS) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDerivativeOperatorsOfNodalForce(elems, postFcn, X, V, Pin, gradP, dForcedX, dForcedGradP, dForcedV, dForcedS);
  }

  template<int dim>
  void computeDerivativeOfNodalHeatPower(ElemSet& elems, PostFcn* postFcn, SVec<double,3>& X, SVec<double,3>& dX, 
				 Vec<double>& d2wall, double* Vwall, double* dVwall, SVec<double,dim>& V, SVec<double,dim>& dV, double dS[3], Vec<double>& dP) {

    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDerivativeOfNodalHeatPower(elems, postFcn, X, dX, d2wall, Vwall, dVwall, V, dV, dS, dP);
  }

  template<int dim>
  void computeDerivativeOfForceAndMoment(ElemSet &elems, PostFcn *postFcn,
                                         SVec<double,3> &X, SVec<double,3> &dX,
                                         Vec<double> &d2wall, double *Vwall, double *dVwall,
                                         SVec<double,dim> &V, SVec<double,dim> &dV, double dS[3],
                                         Vec3D &x0, Vec3D &dFi, Vec3D &dMi, Vec3D &dFv, Vec3D &dMv,
                                         double* gradP[3], double* dGradP[3], int hydro) {

    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDerivativeOfForceAndMoment(elems, postFcn, X, dX, d2wall, Vwall, dVwall, V, dV, dS, x0, dFi, dMi, dFv, dMv, gradP, dGradP, hydro);
  }
/*
  template<int dim>
  void computeDerivativeOfForceAndMoment2(ElemSet &elems, PostFcn *postFcn,
                                             SVec<double,3> &X, SVec<double,3> &dX,
                                             Vec<double> &d2wall, double *Vwall, double *dVwall,
                                             SVec<double,dim> &V, SVec<double,dim> &dV, double dS[3],
                                             Vec3D &x0, Vec3D &dFi, Vec3D &dMi, Vec3D &dFv, Vec3D &dMv, 
				             double* gradP[3], double* dGradP[3], int hydro) {

    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDerivativeOfForceAndMoment2(elems, postFcn, X, dX, d2wall, Vwall, dVwall, V, dV, dS, x0, dFi, dMi, dFv, dMv, gradP, dGradP, hydro);
  }
*/
  template<int dim>
  void computeDerivativeOperatorsOfForceAndMoment(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X,
                                                  Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V, Vec3D &x0,
                                                  double* gradP[3], int hydro,
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
                                                  RectangularSparseMat<double,dim,3> &dMvdV) {
    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDerivativeOperatorsOfForceAndMoment(elems, postFcn, X, d2wall, Vwall, V, x0, gradP, hydro,
                                                        dFidGradP, dFidX, dFidV, dFvdX, dFvdV, dFidS,
                                                        dMidGradP, dMidX, dMidV, dMidS, dMvdX, dMvdV);
  }

  template<int dim>
  void computeDerivativeOfGalerkinTerm(
         ElemSet &elems, FemEquationTerm *fet,
         SVec<double,3> &X, SVec<double,3> &dX,
         Vec<double> &d2wall, double *Vwall, double *dVwall,
         SVec<double,dim> &V, SVec<double,dim> &dV, double dMach,
         SVec<double,dim> &dR) {

    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDerivativeOfGalerkinTerm(elems, fet, X, dX, d2wall, Vwall, dVwall, V, dV, dMach, dR);
  }
  

  template<int dim>
  void computeDerivativeOfGalerkinTermEmb(
         ElemSet &elems, FemEquationTerm *fet,
         SVec<double,3> &X, SVec<double,3> &dX,
         Vec<double> &d2wall, double *Vwall, double *dVwall,
         SVec<double,dim> &V, SVec<double,dim> &dV, double dMach,
         SVec<double,dim> &dR,
         LevelSetStructure *LSS) {

    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeDerivativeOfGalerkinTermEmb(elems, fet, X, dX, d2wall, Vwall, dVwall, V, dV, dMach, dR, LSS);
  }



  template<int dim>
  void computeBCsJacobianWallValues(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, Vec<double> &d2wall, 
                                        double *Vwall, double *dVwall, SVec<double,dim> &V) {

    FaceHelper_dim<dim> h;
    char xx[64];
    GenFaceWrapper_dim<dim> *wrapper=
      (GenFaceWrapper_dim<dim> *)getWrapper_dim(&h, 64, xx);
    wrapper->computeBCsJacobianWallValues(elems, fet, X, d2wall, Vwall, dVwall, V);
  }

  template<int dim>
    void computeHHBoundaryTermResidual(SVec<double,dim> &U,double* Ub,double& UbHH,double& res, VarFcn* vf);

  template<class Scalar,int dim,int neq>
    inline
    void computeHHBoundaryTermJacobian(int faceid,FluxFcn **fluxFcn, SVec<double,dim> &U,
				       double* Ub, GenMat<Scalar,neq> &A, VarFcn* vf,
                                       double& UbHH,
				       Vec<Vec3D> &normals, Vec<double> &normalVel);

  template <class Scalar,int dim>
    void computeMatVecProdH1FarFieldHH(int, GenMat<Scalar,dim> &A, SVec<double,dim> &p_u,
	                      SVec<double,dim> &prod_u,double& p_hh, double& prod_hh);

};


//--------------- FACE CLASS ---------------------------------------------------
class FaceDummy :  public Face {

public:

  // Ensure that when "Virtual template" are not defined in derived classes
  // an error is thrown to avoid infinite loop in wrapper function
  
  template<int dim>
  void computeNodalForce(ElemSet &elems,
			 PostFcn *postFcn, SVec<double,3> &X, 
			 Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V,
			 double pin, SVec<double,3> &F, double* gradP[3]) {
    fprintf(stderr, "Error: undefined function for this face type\n"); exit(1);
  }

  template<int dim>
  void computeNodalHeatPower(ElemSet &elems,
			     PostFcn* postFcn, SVec<double,3>& X, 
			     Vec<double>& d2wall, double* Vwall, 
			     SVec<double,dim>& V, Vec<double>& P) {
    fprintf(stderr, "Error: undefined function for this face type\n"); exit(1);
  }

  template<int dim>
  double computeNodalHeatFluxes(ElemSet &elems,
                             PostFcn* postFcn, SVec<double,3>& X,
                             Vec<double>& d2wall, double* Vwall,
                             SVec<double,dim>& V) {
    fprintf(stderr, "Error: undefined function for this face type\n"); exit(1);
  }


  template<int dim>
  void computeNodalHeatFluxRelatedValues(ElemSet &elems, PostFcn* postFcn, SVec<double,3>& X,
                                 Vec<double>& d2wall, double* Vwall, SVec<double,dim>& V, Vec<double>& P, Vec<double>& N, bool includeKappa){
    fprintf(stderr, "Error: undefined function for this face type\n"); exit(1);
  }

  template<int dim>
  void computeForceAndMoment(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X, 
			     Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V, 
			     Vec3D &x0, Vec3D &Fi, Vec3D &Mi, Vec3D &Fv, Vec3D &Mv, 
			     double* gradP[3], int hydro, SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX,
                                        Vec<double> *genCF) {
    fprintf(stderr, "Error: undefined function for this face type\n"); exit(1);
  }

  template<int dim>
  void computeForceAndMoment(ExactRiemannSolver<dim>& riemann,
                             VarFcn *varFcn,  Vec<Vec3D> &n, Vec<double> &nVel,
                             ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X,
                             Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V,
                             Vec3D &x0, Vec3D &Fi, Vec3D &Mi, Vec3D &Fv, Vec3D &Mv,
                             double* gradP[3], int hydro, SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX,
                                        Vec<double> *genCF) {
    fprintf(stderr, "Error: undifined function for this face type\n"); exit(1);
  }
  
  /* WARNING : THIS FUNCTION IS RETURNING A DOUBLE ? ... IS THIS A PROBLEM ? */
  template<int dim>
  double computeInterfaceWork(ElemSet &elems, PostFcn* postFcn, 
			      SVec<double,3>& X, Vec<double>& d2wall, double ndot, 
			      double* Vwall, SVec<double,dim>& V, double pin) {
    fprintf(stderr, "Error: undefined function for this face type\n"); exit(1);
  }

  template<int dim>
  void computeScalarQuantity(PostFcn::ScalarType type, ElemSet &elems, PostFcn *postFcn, 
			     SVec<double,3> &X, Vec<double> &d2wall, double *Vwall, 
			     SVec<double,dim> &V, SVec<double,2> &Q) {
    fprintf(stderr, "Error: undefined function for this face type\n"); exit(1);
  }

  template<int dim>
  void computeGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, 
			   Vec<double> &d2wall, double *Vwall,
			   SVec<double,dim> &V, SVec<double,dim> &R, LevelSetStructure *LSS=0) {
    fprintf(stderr, "Error: undefined function for this face type\n"); exit(1);
  }
  
  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, 
				   SVec<double,3> &X, Vec<double> &ctrlVol,
				   Vec<double> &d2wall, double *Vwall, 
				   SVec<double,dim> &V, GenMat<Scalar,neq> &A) {
    fprintf(stderr, "Error: undefined function for this face type\n"); exit(1);
  }

  template<int dim>
  void computeForceDerivs(ElemSet &elems, VarFcn *varFcn, SVec<double,3> &X, 
			  SVec<double,dim> &V, SVec<double,dim> &deltaU, Vec<double> &modalF, 
			  SVec<double,3> **localMX) {
    fprintf(stderr, "Error: undefined function for this face type\n"); exit(1);
  }

  template<int dim>
  void computeForceCoefficients(PostFcn *postFcn, Vec3D &x0, ElemSet &elems, 
				SVec<double,3> &X, SVec<double,dim> &V, Vec<double> &d2wall, 
				SVec<double, dim> &Vwall, double pInfty, Vec3D &CFi, Vec3D &CMi, 
				Vec3D &CFv, Vec3D &CMv, double* gradP[3], VecSet< SVec<double,3> > *mX,
                                        Vec<double> *genCF) {
    fprintf(stderr, "Error: undefined function for this face type\n"); exit(1);
  }

  template<int dim>
  void computeFDerivs(ElemSet &elems, VarFcn *varFcn, SVec<double,3> &X, 
		      SVec<double,dim> &Vgl, Vec3D (*F)) {
    fprintf(stderr, "Error: undefined function for this face type\n"); exit(1);
  }

// Included (MB)
  template<int dim>
  void computeDerivativeOfNodalForce(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X, SVec<double,3> &dX,
			     Vec<double> &d2wall, double *Vwall, double *dVwall, SVec<double,dim> &V, SVec<double,dim> &dV,
			     double pin, double dS[3], SVec<double,3> &dF, double* gradP[3], double* dGradP[3]) {
    fprintf(stderr, "Error: undefined function (computeDerivativeOfNodalForce) for this face type\n"); exit(1);
  }

  template<int dim>
  void computeDerivativeOperatorsOfNodalForce(ElemSet *elems, PostFcn *postFcn, SVec<double,3> &X,
                          SVec<double,dim> &V, double pin, double* gradP[3],
                          RectangularSparseMat<double,3,3> &dForcedX,
                          RectangularSparseMat<double,3,3> &dForcedGradP,
                          RectangularSparseMat<double,dim,3> &dForcedV,
                          RectangularSparseMat<double,3,3> &dForcedS) {
    fprintf(stderr, "Error: undefined function (computeDerivativeOperatorsOfNodalForce) for this face type\n"); exit(1);
  }

  template<int dim>
  void computeDerivativeOfNodalHeatPower(ElemSet& elems, PostFcn* postFcn, SVec<double,3>& X, SVec<double,3>& dX, 
                                         Vec<double>& d2wall, double* Vwall, double* dVwall, SVec<double,dim>& V, 
                                         SVec<double,dim>& dV, double dS[3], Vec<double>& dP) {

    fprintf(stderr, "Error: undefined function (computeDerivativeOfNodalHeatPower) for this face type\n"); exit(1);
  }

  template<int dim>
  void computeDerivativeOfForceAndMoment(ElemSet &elems, PostFcn *postFcn,
                                             SVec<double,3> &X, SVec<double,3> &dX,
                                             Vec<double> &d2wall, double *Vwall, double *dVwall,
                                             SVec<double,dim> &V, SVec<double,dim> &dV, double dS[3],
                                             Vec3D &x0, Vec3D &dFi, Vec3D &dMi, Vec3D &dFv, Vec3D &dMv, 
				             double* gradP[3], double* dGradP[3], int hydro) {

    fprintf(stderr, "Error: undefined function (computeDerivativeOfForceAndMoment) for this face type\n"); exit(1);
  }
  /*
    template<int dim>
    void computeDerivativeOfForceAndMoment2(ElemSet &elems, PostFcn *postFcn,
                                               SVec<double,3> &X, SVec<double,3> &dX,
                                               Vec<double> &d2wall, double *Vwall, double *dVwall,
                                               SVec<double,dim> &V, SVec<double,dim> &dV, double dS[3],
                                               Vec3D &x0, Vec3D &dFi, Vec3D &dMi, Vec3D &dFv, Vec3D &dMv,
  				             double* gradP[3], double* dGradP[3], int hydro) {

      fprintf(stderr, "Error: undefined function (computeDerivativeOfForceAndMoment2) for this face type\n"); exit(1);
    }
  */
  template<int dim>
  void computeDerivativeOperatorsOfForceAndMoment(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X,
                                                  Vec<double> &d2wall, double *Vwall, SVec<double,dim> &V, Vec3D &x0,
                                                  double* gradP[3], int hydro,
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
                                                  RectangularSparseMat<double,dim,3> &dMvdV) {
    fprintf(stderr, "Error: undefined function (computeDerivativeOperatorsOfForceAndMoment) for this face type\n"); exit(1);
  }


  template<int dim>
  void computeDerivativeOfGalerkinTerm(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, SVec<double,3> &dX,
                                       Vec<double> &d2wall, double *Vwall, double *dVwall, SVec<double,dim> &V, 
                                       SVec<double,dim> &dV, double dMach, SVec<double,dim> &dR) {
    fprintf(stderr, "Error: undefined function (computeDerivativeOfGalerkinTerm) for this face type\n"); exit(1);
  }

  template<int dim>
  void computeDerivativeOfGalerkinTermEmb(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, SVec<double,3> &dX,
                                       Vec<double> &d2wall, double *Vwall, double *dVwall, SVec<double,dim> &V,
                                       SVec<double,dim> &dV, double dMach, SVec<double,dim> &dR, LevelSetStructure *LSS) {
    fprintf(stderr, "Error: undefined function (computeDerivativeOfGalerkinTermEmb) for this face type\n"); exit(1);
  }

  template<int dim>
  void computeBCsJacobianWallValues(ElemSet &elems, FemEquationTerm *fet, SVec<double,3> &X, Vec<double> &d2wall, 
                                    double *Vwall, double *dVwall, SVec<double,dim> &V) {
    fprintf(stderr, "Error: undefined function (computeBCsJacobianWallValues) for this face type\n"); exit(1);
  }

};


//------------------------------------------------------------------------------
class FaceSet {

  int numFaces;
  int numFaceNorms;
  int numSampledFaces;

  Face **faces;
  BlockAlloc memFaces;
  
	bool sampleMesh;

	class HigherOrderMultiFluid* higherOrderMF;

public:

  /* Need to define constructor and destructor! */
  FaceSet(int);
  ~FaceSet();

  BlockAlloc& getBlockAllocator() { return memFaces; }

  Face &operator[](int i) { return *faces[i]; }

  void addFace(int i, Face *face) { faces[i] = face; }

  int size() const { return numFaces; }

  // Number of face normals to be stored
  int sizeNorms() const { return numFaceNorms; }

  int read(BinFileHandler &, int, int (*)[2], int *);

  

  template<int dim>
  void computeTimeStep(VarFcn *, GeoState &, SVec<double,dim> &, Vec<double> &,
                       TimeLowMachPrec &, Vec<int> &);
  template<int dim>
  void computeTimeStep(FemEquationTerm *, VarFcn *, GeoState &, 
		       SVec<double,3> &, SVec<double,dim> &, Vec<double> &,
                       Vec<double> &, TimeLowMachPrec &, LevelSetStructure *LSS=0);

  template<int dim>
  void computeFiniteVolumeTerm(FluxFcn **, BcData<dim> &, GeoState &, 
			       SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann, 
                               FluxFcn **, BcData<dim> &, GeoState &, 
			       SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeFiniteVolumeTerm(FluxFcn **, BcData<dim> &, GeoState &, 
			       SVec<double,dim> &, Vec<int> &, 
                               SVec<double,dim> &, LevelSetStructure* =0);

  template<int dim, int dimLS>
  void computeFiniteVolumeTermLS(FluxFcn **, BcData<dim> &, GeoState &, 
				 SVec<double,dim> &, SVec<double,dimLS> &, SVec<double,dimLS> &);

  // DEBUG /* Not implemented */
  template<int dim>
  void computeInviscidFluxes(FluxFcn **, BcData<dim> &, GeoState &,
                             SVec<double,dim> &, SVec<double,dim> &);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, BcData<dim> &, GeoState &, 
				       SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim> &, 
                                       FluxFcn **, BcData<dim> &, GeoState &, 
				       SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, BcData<dim> &, GeoState &,
                                       SVec<double,dim> &, GenMat<Scalar,neq> &, int*);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, BcData<dim> &, GeoState &, 
				       SVec<double,dim> &, GenMat<Scalar,neq> &,
                                       Vec<int> &,LevelSetStructure* = 0);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, BcData<dim> &, GeoState &,
                                       SVec<double,dim> &, GenMat<Scalar,neq> &, 
                                       Vec<int> &, int*);

  template<int dim, class Scalar, int dimLS>
  void computeJacobianFiniteVolumeTermLS(GeoState &, 
					 SVec<double,dim> &, GenMat<Scalar,dimLS> &);

  template<int dim>
  void computeGalerkinTerm(ElemSet &, FemEquationTerm *, BcData<dim> &, GeoState &, 
			   SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &, LevelSetStructure *LSS=0);

  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(ElemSet &, FemEquationTerm *, BcData<dim> &, 
				   GeoState &, SVec<double,3> &, Vec<double> &, 
				   SVec<double,dim> &, GenMat<Scalar,neq> &);

// Included (MB)
  template<int dim>
  void computeDerivativeOfFiniteVolumeTerm(FluxFcn **, BcData<dim> &, GeoState &,
                                           SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeDerivativeOfFiniteVolumeTerm(RectangularSparseMat<double,3,dim> *dFluxdFaceNormal,
                                           RectangularSparseMat<double,1,dim> *dFluxdFaceNormalVel,
                                           RectangularSparseMat<double,dim,dim> *dFluxdUb,
                                           BcData<dim> &, GeoState &, 
                                           Vec<Vec3D>&, Vec<double>&, SVec<double,dim> &);

  template<int dim>
  void computeTransposeDerivativeOfFiniteVolumeTerm(RectangularSparseMat<double,3,dim> *dFluxdFaceNormal,
                                                    RectangularSparseMat<double,1,dim> *dFluxdFaceNormalVel,
                                                    RectangularSparseMat<double,dim,dim> *dFluxdUb,
                                                    BcData<dim> &, GeoState &, SVec<double,dim> &,
                                                    Vec<Vec3D>&, Vec<double>&);

  template<int dim>
  void computeDerivativeOperatorsOfFiniteVolumeTerm(FluxFcn **fluxFcn, BcData<dim> &bcData,
              GeoState &geoState, SVec<double,dim> &V,
              RectangularSparseMat<double,3,dim> &dFluxdFaceNormal,
              RectangularSparseMat<double,1,dim> &dFluxdFaceNormalVel,
              RectangularSparseMat<double,dim,dim> &dFluxdUb);

  template<int dim>
  void computeDerivativeOfGalerkinTerm(ElemSet &, FemEquationTerm *, BcData<dim> &, GeoState &,
                                       SVec<double,3> &, SVec<double,3> &, SVec<double,dim> &, 
                                       SVec<double,dim> &, double, SVec<double,dim> &);

  template<int dim>
  void computeDerivativeOfGalerkinTermEmb(ElemSet &, FemEquationTerm *, BcData<dim> &, GeoState &,
                                       SVec<double,3> &, SVec<double,3> &, SVec<double,dim> &,
                                       SVec<double,dim> &, double, SVec<double,dim> &, LevelSetStructure *LSS);



  template<int dim>
  void computeBCsJacobianWallValues(ElemSet &, FemEquationTerm *, BcData<dim> &, 
				   GeoState &, SVec<double,3> &, SVec<double,dim> &);
  template<int dim>
  void computeDerivativeOfTimeStep(FemEquationTerm *, VarFcn *, GeoState &, 
			      SVec<double,3> &, SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &, 
			      Vec<double> &, Vec<double> &, double, 
                              TimeLowMachPrec &);

  void computeConnectedFaces(const std::vector<int> &);
  std::vector<int> facesConnectedToSampleNode;	// for Gappy ROM

  int getNumSampledFaces() {return numSampledFaces;}

  void attachHigherOrderMF(class HigherOrderMultiFluid*);

  template <int dim>
  void updateHHState(SVec<double,dim>& V, VarFcn* vf, double dt) {

//    for (int i = 0; i < numFaces; ++i)
//      faces[i]->updateHHBoundaryTerm(V, vf, dt);
  }

  template<int dim>
    void computeHHBoundaryTermResidual(BcData<dim> &bcData,SVec<double,dim> &V,Vec<double>& res, VarFcn* vf);

  template<class Scalar,int dim,int neq>
  inline
    void computeHHBoundaryTermJacobian(FluxFcn **fluxFcn, BcData<dim> &bcData,SVec<double,dim> &U,
				       GeoState& geoState,
				       GenMat<Scalar,neq> &A, VarFcn* vf);

  
  template<class Scalar, int dim>
  void computeMatVecProdH1FarFieldHH(GenMat<Scalar,dim> &A, SVec<double,dim> &p_u,
 	                      SVec<double,dim> &prod_u,Vec<double>& p_hh, 
                              Vec<double>& prod_hh);
 //void initializeHHCoeffs(double cc) {for(int i=0; i<numFaces; i++) faces[i]->initializeHHCoeffs(cc);}
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <Face.C>
#endif

#include <FaceTria.h>

#endif
