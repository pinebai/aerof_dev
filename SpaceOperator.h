#ifndef _SPACE_OPERATOR_H_
#define _SPACE_OPERATOR_H_

#include <IoData.h>
#include <GhostPoint.h>
#include <RestrictionMapping.h>
#include <complex>

#include <TriangulatedInterface.h>

typedef std::complex<double> bcomp;

class VarFcn;
class FluidSelector;
class BcFcn;
class RecFcn;
class FluxFcn;
class FemEquationTerm;
class VolumicForceTerm;
class SmagorinskyLESTerm;
class WaleLESTerm;
class Domain;
class DistGeoState;
class Communicator;
class Timer;
class DistStructureLevelSet;

template<int dim> class DistVMSLESTerm;
template<int dim> class DistDynamicVMSTerm;
template<int dim> class DistDynamicLESTerm;
template<int dim> class DistEdgeGrad;
template<int dim> class DistExtrapolation;
template<int dim> class DistBcData;
template<int dim> class DistTimeState;
template<class Scalar, int dim> class DistSVec;
template<class Scalar, int dim> class DistMat;
template<int dim> class DistExactRiemannSolver;
template<int dim> struct dRdXoperators;

#ifndef _DNDGRAD_TMPL_
#define _DNDGRAD_TMPL_
template<int dim, class Scalar = double> class DistNodalGrad;
#endif

//------------------------------------------------------------------------------

template<int dim>
class SpaceOperator {

protected:

  VarFcn *varFcn;
  DistBcData<dim> *bcData;
  DistGeoState *geoState;
  DistSVec<double,dim> *V;

// Included (MB)
  DistSVec<double,dim> *dU;
  DistSVec<double,dim> *dV;
  DistSVec<double,dim> *dRm;
  int jacobianAA;
  int jacobianSA;
  int reconstructionAA;
  int reconstructionSA;
  FluxFcn **fluxFcnAA;
  RecFcn *recFcnAA;
  FluxFcn **fluxFcnSA;
  RecFcn *recFcnSA;

protected:

  bool locAlloc;

  BcFcn *bcFcn;
  FluxFcn **fluxFcn;
  RecFcn *recFcn;

  DistNodalGrad<dim, double> *ngrad;
  DistNodalGrad<dim, bcomp> *compNodalGrad;

  DistNodalGrad<dim, double> *ngraddV;

  DistEdgeGrad<dim> *egrad;
  DistExtrapolation<dim> *xpol;
  DistVMSLESTerm<dim> *vms;
  DistDynamicLESTerm<dim> *dles;
  DistDynamicVMSTerm<dim> *dvms;
  FemEquationTerm *fet;
  SmagorinskyLESTerm *smag;
  WaleLESTerm *wale;
  VolumicForceTerm *volForce;

  Domain *domain;

  Timer *timer;
  Communicator *com;

  bool use_modal;
  bool use_complex;
  int order;
  int failsafe;
  int rshift;
// Included (MB)
  IoData *iod;

  enum DescriptorCase {
    DESCRIPTOR, HYBRID, NONDESCRIPTOR
  };
  DescriptorCase descriptorCase;

private:
  bool externalSI; //d2d
  int ccc;
public:

  SpaceOperator(IoData &, VarFcn *, DistBcData<dim> *, DistGeoState *,
		Domain *, DistSVec<double,dim> * = 0);
  SpaceOperator(const SpaceOperator<dim> &, bool);
  ~SpaceOperator();

  DistNodalGrad<dim, double> *getDistNodalGrad(DistSVec<double,dim> &)  { return ngrad; }
  DistNodalGrad<dim, bcomp> *getDistNodalGrad(DistSVec<bcomp,dim> &)  { return compNodalGrad; }

  void updateFixes() { ngrad->updateFixes(); }

  int getSpaceOrder() {return order;}

  int **getNodeType() {return domain->getNodeType(); }

  FluxFcn** getFluxFcn() { return fluxFcn; }
  DistBcData<dim>* getDistBcData() { return bcData; }

  DistGeoState* getGeoState() { return geoState; }

  BcFcn *createBcFcn(IoData &);
  FluxFcn **createFluxFcn(IoData &);
  RecFcn *createRecFcn(IoData &);
  //RecFcn *createRecFcnLS(IoData &);
  FemEquationTerm *createFemEquationTerm(IoData &);
  VolumicForceTerm *createVolumicForceTerm(IoData &);
  VarFcn* getVarFcn() { return varFcn; }
  void setBcFcn(BcFcn *);
  void setFluxFcn(FluxFcn **);
  void setRecFcn(RecFcn *);
  void setFemEquationTerm(FemEquationTerm *);
  void fix(DistSVec<bool,2>&);
  void resetTag();

  Domain* getDomain() { return domain; }

  BcFcn* getBcFcn() { return bcFcn; }

  DistSVec<double,dim>* getCurrentPrimitiveVector() { return V; }

  FemEquationTerm *getFemEquationTerm() { return fet;}

  void conservativeToPrimitive(DistSVec<double,dim> &U, DistVec<int>* fid = 0) 
  {
	  varFcn->conservativeToPrimitive(U, *V, fid);
  }

  void conservativeToPrimitive(DistSVec<double,dim> &U, 
										 DistLevelSetStructure *distLSS,
										 DistVec<int>* fid = 0) 
  {
	  varFcn->conservativeToPrimitive(U, *V, distLSS, fid);
  }


// Included (MB)
  void computeResidual(DistSVec<double,3> &, DistVec<double> &,
		       DistSVec<double,dim> &, DistSVec<double,dim> &,
                       DistTimeState<dim> *, bool=true);
  void computeResidualRestrict(DistSVec<double,3> &, DistVec<double> &,
			DistSVec<double,dim> &, DistSVec<double,dim> &, DistTimeState<dim> *,
		        const std::vector< std::vector<int> > &, bool=true);
// Included (MB)
  void computeResidual(DistExactRiemannSolver<dim> *,
		       DistSVec<double,3> &, DistVec<double> &,
		       DistSVec<double,dim> &, DistSVec<double,dim> &,
		       DistTimeState<dim> *, bool=true);
  
  //d2d embedded
  void computeResidual(DistSVec<double,3> &, DistVec<double> &,
                       DistSVec<double,dim> &, DistSVec<double,dim> &,
                       DistSVec<double,dim> &, DistSVec<double,dim> &Wext,
							  DistLevelSetStructure *,
                       bool, bool, DistVec<int> &, 
		       DistSVec<double,dim> &,
                       DistExactRiemannSolver<dim> *, int, 
							  int it = 0, DistVec<GhostPoint<dim>*> *ghostPoints = 0, 
		       bool=true);

  void computeResidual(DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, 
                       DistSVec<double,dim> &, DistSVec<double,dim> &, 
                       DistVec<int> &, DistVec<int> &, DistLevelSetStructure *,
                       bool, bool, DistVec<int> &, DistSVec<double,dim> &,
					   DistExactRiemannSolver<dim> *, int,
					   double, double, int it = 0, DistVec<GhostPoint<dim>*> *ghostPoints = 0);

  void updateSweptNodes(DistSVec<double,3> &X,DistVec<double> &ctrlVol,
			int phaseChangeChoice, int phaseChangeAlg,
                        DistSVec<double,dim> &U, DistSVec<double,dim> &V,
                        DistVec<double> &Weights, DistSVec<double,dim> &VWeights,
                        DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji,
                        DistLevelSetStructure *distLSS, double *vfar, bool limit,
			DistVec<int> *fluidId);

  void populateGhostPoints(DistVec<GhostPoint<dim>*> *ghostPoints,DistSVec<double,3> &X,DistSVec<double,dim> &U,VarFcn *varFcn,DistLevelSetStructure *distLSS,bool linFSI,DistVec<int> &tag);
  
  template <int neq>
  void populateGhostPoints(DistVec<GhostPoint<dim>*> *ghostPoints,DistSVec<double,neq> &U,VarFcn *varFcn,DistLevelSetStructure *distLSS,DistVec<int> &tag) {
    fprintf(stderr,"PopulateGhostPoints<%d> not implemented!\n",neq);
    exit(-1);
  }


  // d2d
  void setSIstencil( DistSVec<double,3> &X, DistLevelSetStructure *distLSS, DistVec<int> &fluidId, DistSVec<double,dim> &U);
  void setFEMstencil(DistSVec<double,3> &X, DistLevelSetStructure *distLSS, DistVec<int> &fluidId, DistSVec<double,dim> &U);

  void computeRiemannWeightsForEmbeddedStruct(DistSVec<double,3> &X,
                           DistSVec<double,dim> &U, DistSVec<double,dim> &V,
                           DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji,
                           DistVec<double> &Weights, DistSVec<double,dim> &VWeights,
                           DistLevelSetStructure *distLSS, DistVec<int> *fluidId =  0);


  double recomputeResidual(DistSVec<double,dim> &, DistSVec<double,dim> &);
 
  double computeRealFluidResidual(DistSVec<double, dim> &, DistSVec<double,dim> &, DistLevelSetStructure &);

  void recomputeRHS(DistSVec<double,3> &, DistSVec<double,dim> &,
                                     DistSVec<double,dim> &);
  void recomputeRHS(DistSVec<double,3> &, DistSVec<double,dim> &,
                    DistVec<int> &, DistSVec<double,dim> &);


  void computePostOpDVMS(DistSVec<double,3> &, DistVec<double> &,
                         DistSVec<double,dim> &, DistVec<double> *,
                         DistTimeState<dim> *);

  template<class Scalar, int neq>
  void computeJacobian(DistSVec<double,3> &, DistVec<double> &,
		       DistSVec<double,dim> &, DistMat<Scalar,neq> &,
		       DistTimeState<dim> *);

  template<class Scalar, int neq>
  void computeJacobian(DistExactRiemannSolver<dim> *riemann, DistSVec<double,3> &, DistVec<double> &,
		       DistSVec<double,dim> &, DistMat<Scalar,neq> &,
		       DistTimeState<dim> *);

  template <class Scalar,int neq>
  void computeJacobian(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
                       DistSVec<double,dim> &U,
                       DistLevelSetStructure *distLSS,
                       DistVec<int> &fluidId, 
                       DistExactRiemannSolver<dim> *riemann, 
                       int Nriemann,
                       DistVec<GhostPoint<dim>*> *ghostPoints,
                       DistMat<Scalar,neq>& A,
                       DistTimeState<dim>*);
  
  void getExtrapolationValue(DistSVec<double,dim>&, DistSVec<double,dim>&, DistSVec<double,3>&);
  void applyExtrapolationToSolutionVector(DistSVec<double,dim>&, DistSVec<double,dim>&);

  template<class Scalar, int neq>
  void computeViscousJacobian(DistSVec<double,3> &, DistVec<double> &, DistMat<Scalar,neq> &);

  void applyBCsToSolutionVector(DistSVec<double,dim> &,DistLevelSetStructure *distLSS=0);

  void applyBCsToTurbSolutionVector(DistSVec<double,dim> &,DistLevelSetStructure *distLSS=0);

  void applyBCsToResidual(DistSVec<double,dim> &, DistSVec<double,dim> &, DistLevelSetStructure *distLSS=0);

  template<class Scalar, int neq>
  void applyBCsToJacobian(DistSVec<double,dim> &, DistMat<Scalar,neq> &, DistLevelSetStructure *distLSS=0);

  template<class Scalar, int neq>
  void applyBCsToH2Jacobian(DistSVec<double,dim> &, DistMat<Scalar,neq> &);

  template<class Scalar>
  void computeH1(DistSVec<double,3> &, DistVec<double> &,
                 DistSVec<double,dim> &, DistMat<Scalar,dim> &);

  template<class Scalar, int neq>
  void computeH2(DistSVec<double,3> &, DistVec<double> &,
		 DistSVec<double,dim> &, DistMat<Scalar,neq> &,
		 DistSVec<double,dim> &, DistSVec<double,dim> &,
		 DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<class Scalar, int neq>
  void computeH2transpose(DistSVec<double,3> &, DistVec<double> &,
                          DistSVec<double,dim> &, DistMat<Scalar,neq> &,
                          DistSVec<double,dim> &, DistSVec<double,dim> &,
                          DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<class Scalar, int neq>
  void computeH2(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
		 DistSVec<double,dim> &U, 
		 DistLevelSetStructure *distLSS,
		 DistVec<int> &fluidId, 
		 DistExactRiemannSolver<dim> *riemann, 
		 int Nriemann,
		 DistVec<GhostPoint<dim>*> *ghostPoints,
		 DistMat<Scalar,neq> &H2,
		 DistSVec<double,dim> &aij, DistSVec<double,dim> &aji,
		 DistSVec<double,dim> &bij, DistSVec<double,dim> &bji,
		 DistSVec<double,dim> &betaij, DistSVec<double,dim> &betaji);

  template<class Scalar1, class Scalar2>
  void applyH2(DistSVec<double,3> &, DistVec<double> &, 
	       DistSVec<double,dim> &,
               DistMat<Scalar1,dim> &, 
	       DistSVec<double,dim> &, DistSVec<double,dim> &,
               DistSVec<double,dim> &, DistSVec<double,dim> &,
               DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);
  
  template<class Scalar1, class Scalar2>
  void applyH2transposeNew(DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &,
                           DistMat<Scalar1,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
                           DistSVec<double,dim> &, DistSVec<double,dim> &,
                           DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);

  template<class Scalar1, class Scalar2>
  void applyH2(DistSVec<double,3> &X, DistVec<double> &ctrlVol,
	       DistSVec<double,dim> &U,   
	       DistLevelSetStructure *distLSS,
	       DistVec<int> &fluidId,
	       bool linRecAtInterface, bool viscSecOrder,
	       DistExactRiemannSolver<dim> *riemann, 
	       int Nriemann,
	       DistVec<GhostPoint<dim>*> *ghostPoints,
	       DistMat<Scalar1,dim> &H2,
	       DistSVec<double,dim> &aij, DistSVec<double,dim> &aji,
	       DistSVec<double,dim> &bij, DistSVec<double,dim> &bji,
	       DistSVec<double,dim> &betaij, DistSVec<double,dim> &betaji,
	       DistSVec<Scalar2,dim> &p, DistSVec<Scalar2,dim> &prod);
  
  template<class Scalar1, class Scalar2>
  void applyH2T(DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistMat<Scalar1,dim> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &,
                DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);

  template<class Scalar, int neq>
  void printAllMatrix(DistMat<Scalar,neq> &, int);

  void printAllVariable(DistSVec<double,3>&, DistSVec<double,dim> &, int);
  void printVariable(DistSVec<double,dim> &);


  // Included (MB)
  void rstFluxFcn(IoData &);

  // Included (YC)
  void computeDerivativeOperators(Vec3D &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, double, DistSVec<double,dim> &, DistVec<double> &,
                                  DistTimeState<dim> * = 0, PostOperator<dim> * = 0, dRdXoperators<dim> * =0);

  // Included (MB)
  void computeDerivativeOfResidual(
	   DistSVec<double,3> &,   //X->nodal position
	   DistSVec<double,3> &,   //dX->derivative of nodal postion
	   DistVec<double> &,      //ctrlVol->assiciated control volume
	   DistVec<double> &,      //dCtrlVol->derivative od associated control volume
       DistSVec<double,dim> &, //U->fluid state vector in conservative form
	   double,                 //dmach
	   DistSVec<double,dim> &, //R->Residual of the fluid equations
	   DistSVec<double,dim> &,   //dR->derivative of the fluid equation residual
	   DistTimeState<dim> * = 0);//timeState


  void computeDerivativeOfResidual(dRdXoperators<dim> *, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double> &, DistVec<double> &, DistVec<Vec3D>&, DistVec<Vec3D>&, DistVec<double>&,
                                   DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,6> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &, double, DistTimeState<dim> *);

  void computeTransposeDerivativeOfResidual(dRdXoperators<dim> *, DistSVec<double,dim> &, DistSVec<double,dim> &,
                                            DistVec<double> &, DistVec<double> &, DistSVec<double,3> &,
                                            DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistVec<Vec3D>&,
                                            DistVec<Vec3D>&, DistVec<double>&, DistSVec<double,6>& );

  void computeDerivativeOfResidual(DistSVec<double,3> &X,
				   DistVec<double> &ctrlVol,
				   DistSVec<double,dim> &U,
				   DistLevelSetStructure *distLSS,
				   bool linRecAtInterface, bool viscSecOrder, 
				   DistVec<int> &fluidId, 
				   DistExactRiemannSolver<dim> *riemann, 
				   int Nriemann,
				   DistVec<GhostPoint<dim>*> *ghostPoints,
				   double dMach,
				   DistSVec<double,dim> &R, DistSVec<double,dim> &dR,
				   DistTimeState<dim> *timeState);

  // Included (MB)
  void applyBCsToDerivativeOfResidual(DistSVec<double,dim> &, DistSVec<double,dim> &);


  // Included (MB)
  void rstVarFet(IoData &ioData) 
  {
    if (fet) fet->rstVar(ioData, com);
  }


  // Included (MB)
  bool useModal() 
  {return use_modal;}

  // Included (MB)
  /// \note This function is implemented.
  /// It is called only from MatVecProdFD::evaluateInviscid and
  /// MatVecProdFD::applyInviscid.
  void computeInviscidResidual(DistSVec<double,3> &, DistVec<double> &,
		       DistSVec<double,dim> &, DistSVec<double,dim> &, DistTimeState<dim> * = 0, bool=true);


  // Included (MB)
  /// \note This function is implemented.
  /// It is called only from MatVecProdFD::evaluateViscous and
  /// MatVecProdFD::applyViscous.
  void computeViscousResidual(DistSVec<double,3> &, DistVec<double> &,
		       DistSVec<double,dim> &, DistSVec<double,dim> &, DistTimeState<dim> * = 0, bool=true);


  // Included (MB)
  /// \note This routine is not implemented.
  //void computeInviscidResidual
  //(
  //  DistSVec<double,3> &, DistVec<double> &,
  //  DistSVec<double,dim> &, DistVec<double> &,
  //  DistSVec<double,dim> &, bool=true
  //);


  // Included (MB)
  /// \note This routine is not implemented.
  //void computeViscousResidual
  //(
  //  DistSVec<double,3> &, DistVec<double> &,
  //  DistSVec<double,dim> &, DistVec<double> &,
  //  DistSVec<double,dim> &, bool=true
  //);

  template<class Scalar>
  void applyBCsToH2Jacobian(DistSVec<double,dim> &, DistMat<Scalar,dim> &);

  void computeNodalGrad(DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, 
						  DistVec<int> *fluidId=0, DistLevelSetStructure *distLSS=0);



  // Included (MB)
  /// \note This routine is called from FluidSensitivityAnalysis.
  void computeGradP(DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, 
						  DistVec<int> *fluidId=0, DistLevelSetStructure *distLSS=0);


  // Included (MB)
  /// \note This routine is called from FluidSensitivityAnalysis.
  void computeDerivativeOfGradP
  (
    DistSVec<double,3> &X, DistSVec<double,3> &dX,
    DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol,
    DistSVec<double,dim> &U, DistSVec<double,dim> &dU
  );

  // Included (YC)
  /// \note This routine is called from FluidSensitivityAnalysis.
  void computeDerivativeOfGradP
  (
    dRdXoperators<dim> *dRdXop,
    DistSVec<double,3> &dX, DistVec<double> &dCtrlVol, DistSVec<double,dim> &dU,
    DistSVec<double,dim> &dddx, DistSVec<double,dim> &dddy,
    DistSVec<double,dim> &dddz, DistSVec<double,6> &dR,
    DistSVec<double,3> &dGradP
  );

  // Included (YC)
  /// \note This routine is called from FluidSensitivityAnalysis.
  void computeTransposeDerivativeOfGradP
  (
    dRdXoperators<dim> *dRdXop,
    DistSVec<double,3> &dGradP,
    DistSVec<double,dim> &dddx,
    DistSVec<double,dim> &dddy, 
    DistSVec<double,dim> &dddz,
    DistSVec<double,6> &dR,
    DistVec<double> &dCtrlVol,
    DistSVec<double,3> &dX,
	DistSVec<double,dim> &dU,
	bool assembleDX = true
  );

  // Included (MB)
  /// \note This routine is redundant. It is never called.
  //void computeDerivativeOfGradP
  //(
  //  DistSVec<double,3> &X, DistSVec<double,3> &dX,
  //  DistVec<double> &ctrlVol, DistVec<double> &dCtrlVol,
  //  DistSVec<double,dim> &U
  //);


  void computeForceLoad(int forceApp, int orderOfAccuracy, DistSVec<double,3> &X, DistVec<double> &ctrlVol, 
                        double (*Fs)[3], int sizeFs, DistLevelSetStructure *distLSS,
                        DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji, DistSVec<double,dim> *Wextij,
			DistVec<GhostPoint<dim>*> *ghostPoints = 0, PostFcn *postFcn = 0,DistVec<int>* fid = 0);

  void computederivativeOfForceLoad(int forceApp, int orderOfAccuracy, DistSVec<double,3> &X, DistVec<double> &ctrlVol,
   				    double (*dFs)[3], int sizeFs, DistLevelSetStructure *distLSS,
   				    DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji,
				    DistSVec<double,dim> &dV_, double dS[3],
   				    DistVec<GhostPoint<dim>*> *ghostPoints = 0, PostFcn *postFcn = 0,DistVec<int>* fid = 0);
};

//------------------------------------------------------------------------------
// derived class of SpaceOperator to run multi-phase flow problems
// the added functions allow to advance the Euler equations when
// different EOS are considered and to advance the level-set 
// advection equation(s).
template<int dim, int dimLS>
class MultiPhaseSpaceOperator : public SpaceOperator<dim> {

protected:

  RecFcn *recFcnLS;
  DistNodalGrad<dimLS, double> *ngradLS;
  DistEdgeGrad<dimLS>          *egradLS; //d2d

  DistSVec<double, dimLS>* normals[3]; 
  DistNodalGrad<dimLS, double> *ngradLS_second[3];
  DistSVec<double,dimLS>* curvature;

  RecFcn *createRecFcnLS(IoData &);

  struct {

    DistVec<int>* counts[2];
    DistSVec<double,dim>* vals[2];
    DistNodalGrad<dim,double>* cutgrads[2];    
  } higherOrderData;

public:

  MultiPhaseSpaceOperator(IoData &, VarFcn *, DistBcData<dim> *, DistGeoState *,
		Domain *, DistSVec<double,dim> * = 0);
  MultiPhaseSpaceOperator(const MultiPhaseSpaceOperator<dim,dimLS> &, bool);
  ~MultiPhaseSpaceOperator();

  void computeResidual(DistSVec<double,3> &, DistVec<double> &,
                       DistSVec<double,dim> &, DistSVec<double,dimLS> &,
                       FluidSelector &, DistSVec<double,dim> &,
                       DistExactRiemannSolver<dim> *,DistTimeState<dim>*, int it);
  void computeResidual(DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &, 
                       DistSVec<double,dim> &, DistLevelSetStructure *, bool, bool, DistExactRiemannSolver<dim> *, 
                       int, DistSVec<double,dimLS> &, FluidSelector &, 
                       DistSVec<double,dim> &, int, DistVec<GhostPoint<dim>*> *);

  void computeResidualLS(DistSVec<double,3> &, DistVec<double> &,
                         DistSVec<double,dimLS> &, DistVec<int> &, 
                         DistSVec<double,dim> &,DistSVec<double,dimLS> &, DistLevelSetStructure* =0, bool = true,
			 int method = 0, int ls_order = 1);

  template<class Scalar, int neq>
  void computeJacobian(DistSVec<double,3> &, DistVec<double> &,
                       DistSVec<double,dim> &, DistMat<Scalar,neq> &,
                       FluidSelector &, DistExactRiemannSolver<dim> *,DistTimeState<dim>*);

  template<class Scalar, int neq>
  void computeJacobian(DistExactRiemannSolver<dim>* riemann,
                       DistSVec<double,3>& X, DistSVec<double,dim>& U,DistVec<double>& ctrlVol,
                       DistLevelSetStructure *distLSS,
                       int Nriemann,
                       FluidSelector &fluidSelector,
                       DistMat<Scalar,neq>& A,DistTimeState<dim>* timeState);

  template<class Scalar>
  void computeJacobianLS(DistSVec<double,3> &X,DistSVec<double,dim> &V, DistVec<double> &ctrlVol,
			 DistSVec<double,dimLS> &Phi,DistMat<Scalar,dimLS> &A,DistVec<int> &fluidId,DistLevelSetStructure* distLSS,
			 int method);

  // for phase-change update
  void extrapolatePhiV(DistLevelSetStructure *distLSS, DistSVec<double,dimLS> &PhiV);
  void extrapolatePhiV2(DistLevelSetStructure *distLSS, DistSVec<double,dimLS> &PhiV);
  void updateSweptNodes(DistSVec<double,3> &X, int &phaseChangeChoice, DistSVec<double,dim> &U, DistSVec<double,dim> &V,
                        DistVec<double> &Weights, DistSVec<double,dim> &VWeights,
                        DistSVec<double,dimLS> &Phi, DistSVec<double,dimLS> &PhiWeights,
                        DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji,
                        DistLevelSetStructure *distLSS, double *vfar, bool updateWithCracking,
                        DistVec<int> *fluidId0, DistVec<int> *fluidId);
  void resetFirstLayerLevelSetFS(DistSVec<double,dimLS> &PhiV, DistLevelSetStructure *distLSS, DistVec<int> &fluidId, 
                                 DistSVec<bool,2> &Tag);

    void extrapolatePhaseChange(DistSVec<double,3> &X, DistVec<double> &ctrlVol,int phaseChangeAlg,
			   DistSVec<double,dim> &U, DistSVec<double,dim> &V,
			   DistVec<double> &Weights, DistSVec<double,dim> &VWeights,
			   DistLevelSetStructure *distLSS, DistVec<int> &fluidId,
				DistVec<int> &fluidIdn,bool limit = false);

  void attachTriangulatedInterface(TriangulatedInterface*);

  void setLastPhaseChangeValues(DistExactRiemannSolver<dim>*);

};
//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <SpaceOperator.C>
#endif

#endif
