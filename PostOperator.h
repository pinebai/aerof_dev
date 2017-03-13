#ifndef _POST_OPERATOR_H_
#define _POST_OPERATOR_H_

#include <PostFcn.h>
#include <VectorSet.h>
#include <GhostPoint.h>
#include <map>
#include <LevelSet/LevelSetStructure.h>

using std::map;


class IoData;
class VarFcn;
class SubDomain;
class Domain;
class DistGeoState;
class Communicator;
class SmagorinskyLESTerm;
class WaleLESTerm;
class DynamicLESTerm;

struct Vec3D;

template<int dim> class DistBcData;
template<class Scalar> class DistVec;
template<class Scalar, int dim> class DistSVec;
template<int dim> class DistVMSLESTerm;
template<int dim> class DistDynamicLESTerm;
template<int dim> class DistDynamicVMSTerm;
template<int dim> class SpaceOperator;
template<int dim> class DistTimeState;
template<int dim> class DistExactRiemannSolver;

template<int dim> struct dRdXoperators;


template<int dim>
class ForceGenerator {
  public:
    virtual void getForcesAndMoments(map<int,int> &surfMap, DistSVec<double,dim> &U, DistSVec<double,3> &X,
                                           Vec3D* Fi, Vec3D* Mi) = 0;

    virtual void getderivativeOfForcesAndMoments(map<int,int> & surfOutMap, 
						 DistSVec<double,dim> &U, DistSVec<double,dim> &dU, 
						 DistSVec<double,3> &X, double dS[3],
						 Vec3D *Fi, Vec3D *Mi) = 0;

};

//------------------------------------------------------------------------------

template<int dim>
class PostOperator {

  VarFcn *varFcn;
  DistBcData<dim> *bcData;
  DistGeoState *geoState;
  DistSVec<double,dim> *V;

  int numLocSub;
  SubDomain **subDomain;
  Domain *domain;
  Communicator *com;
  int numSurf;
  int numSurfHF;
  map<int,int> surfOutMap;
  map<int,int> surfOutMapHF;
  map<int,int> surfComputeMap; //AS far as I can figure out this map is never used
  ForceGenerator<dim> *forceGen;
  
 // Coefficients to Compute nodal force transfer
  double nodalForceWeights[2];

// Included (MB)
  DistSVec<double,dim> *dV;

private:

  double threshold;
  double refLengthSq;
  double pressInfty;
  DistSVec<double,2>* tmp2;
  SmagorinskyLESTerm *smag;
  WaleLESTerm *wale;  
  DistVMSLESTerm<dim> *vms;
  DistDynamicLESTerm<dim> *dles;
  DistDynamicVMSTerm<dim> *dvms;
  SpaceOperator<dim> *spaceOp;
  CommPattern<double>* vec2Pat;
  PostFcn *postFcn;
  DistVec<double> *mutOmu;
  DistVec<double> *Cs;
  DistVec<double> *CsDvms;
  DistVec<double> *CsDles;
  bool built_dVdU;

public:

  PostOperator(IoData &, VarFcn *, DistBcData<dim> *, DistGeoState *, 
	       Domain *, DistSVec<double,dim> * = 0);
  ~PostOperator();

  void computeNodalForce(DistSVec<double,3> &, DistSVec<double,dim> &, 
			 DistVec<double> &, DistSVec<double,3> &,
			 DistVec<int> * = 0);

  void computeNodalHeatPower(DistSVec<double,3> &, DistSVec<double,dim> &, 
			     DistVec<double> &);
  void computeNodalHeatFluxRelatedValues(DistSVec<double,3> &, DistSVec<double,dim> &,
                                               DistVec<double> &, bool includeKappa);
  void computeForceAndMoment(Vec3D &, DistSVec<double,3> &, DistSVec<double,dim> &,
                             DistVec<int> *,
			     Vec3D *, Vec3D *, Vec3D *, Vec3D *, int = 0, 
                             VecSet< DistSVec<double,3> > *mX = 0, Vec<double> *genCF = 0);
  void computeForceAndMoment(DistExactRiemannSolver<dim>&, 
                             Vec3D &, DistSVec<double,3> &, DistSVec<double,dim> &,
                             DistVec<int> *,
                             Vec3D *, Vec3D *, Vec3D *, Vec3D *, int = 0,
                             VecSet< DistSVec<double,3> > *mX = 0, Vec<double> *genCF = 0);

  void computeHeatFluxes(DistSVec<double,3> &,
                                          DistSVec<double,dim> &, double*);

  double computeInterfaceWork(DistSVec<double,3>&, DistSVec<double,dim>&, DistVec<double>&);


  void computeScalarQuantity(PostFcn::ScalarType, DistSVec<double,3> &,
			     DistSVec<double,dim> &, DistVec<double> &, 
                             DistVec<double> &, DistTimeState<dim> *);
  template<int dimLS>
  void computeScalarQuantity(PostFcn::ScalarType, DistSVec<double,3> &,
                             DistSVec<double,dim> &, DistVec<double> &,
                             DistVec<double> &, DistTimeState<dim> *,
                             DistVec<int> &,DistSVec<double,dimLS>* = NULL);

  template<int dimLS>
  void computeScalarQuantity(PostFcn::ScalarType, DistSVec<double,3> &,
			     DistSVec<double,dim> &, DistVec<double> &,
			     DistTimeState<dim> *,
			     DistVec<int> &,int* subId,int* locNodeId,
			     int* last,int count, double* result,
                             std::vector<Vec3D>& locations,
			     DistSVec<double,dimLS>* = NULL,
                             DistLevelSetStructure *distLSS = 0,
                             DistVec<GhostPoint<dim>*> *ghostPoints = 0);

   void computeCP(DistSVec<double,3>& X, DistSVec<double,dim>& U, Vec3D &cp);
  //void computeScalarQuantity(PostFcn::ScalarType, DistSVec<double,3> &,
  //                           DistSVec<double,dim> &, DistVec<double> &,
  //                           DistSVec<double,1> &, DistVec<int> &);

   template<int dimLS>
  void computeEMBScalarQuantity(DistSVec<double,3>& X,
				DistSVec<double,dim>& U,
				DistVec<double>& A,
				double** EmbQs,
				DistTimeState<dim> *timeState,
				DistVec<int>& fluidId, 
										  DistSVec<double,dim>* Wextij,
				DistSVec<double,dimLS> *Phi, 
				DistLevelSetStructure *distLSS,
										  DistVec<GhostPoint<dim>*> *ghostPoints, bool ExternalSI);

  void computeVectorQuantity(PostFcn::VectorType, DistSVec<double,3> &,
			     DistSVec<double,dim> &, DistSVec<double,3> &);
  void computeVectorQuantity(PostFcn::VectorType, DistSVec<double,3> &,
                             DistSVec<double,dim> &, DistSVec<double,3> &, DistVec<int> &);

  void computeVectorQuantity(PostFcn::VectorType, DistSVec<double,3> &,
                             DistSVec<double,dim> &,
			     int* subId,int* locNodeId,int* last,
			     int count, double* results, 
                             std::vector<Vec3D>& locations,
                             DistVec<int> &, DistLevelSetStructure *distLSS = 0,
                             DistVec<GhostPoint<dim>*> *ghostPoints = 0);

  void computeVectorQuantity(PostFcn::VectorType type,
									  DistSVec<double,3> &X,
									  DistSVec<double,dim> &U,
									  DistSVec<double,3> &Q,
									  DistLevelSetStructure *distLSS,
									  DistVec<int> &fluidId);
	  
  void computeForceDerivs(DistSVec<double,3> &, DistSVec<double,dim> &,
                          DistSVec<double,dim> &,Vec<double> &,VecSet< DistSVec<double, 3> > &);

  void computeForceCoefficients(Vec3D &, DistSVec<double,3> &, DistSVec<double,dim> &,
                                Vec3D &, Vec3D &, Vec3D &, Vec3D &,  
                                VecSet< DistSVec<double,3> > *mX = 0, DistVec<double> *genCF = 0);
  int getNumSurf() { return numSurf; }
  map<int, int> &getSurfMap() { return surfOutMap; }

  int getNumSurfHF() { return numSurfHF; }
  map<int, int> &getSurfMapHF() { return surfOutMapHF; }

  PostFcn* getPostFcn() {return postFcn;}
  
// Included (MB)
  void computeDerivativeOfScalarQuantity(PostFcn::ScalarDerivativeType, double [3], DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistVec<double> &, DistTimeState<dim> *);

  void computeDerivativeOfVectorQuantity(PostFcn::VectorDerivativeType, DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,3> &);

  void computeDerivativeOfForceAndMoment(Vec3D &, DistSVec<double,3> &, DistSVec<double,3> &,
                                                                           DistSVec<double,dim> &, DistSVec<double,dim> &, double [3],
                                                                           Vec3D *, Vec3D *, Vec3D *, Vec3D *, int = 0);

  void computeDerivativeOfForceAndMoment(Vec3D &x0, DistSVec<double,3> &X, 
					 DistSVec<double,dim> &U, 
					 DistSVec<double,dim> &dU, 
					 DistVec<int> *fluidId,
					 double dS[3],
					 Vec3D *dFi, Vec3D *dMi, 
					 Vec3D *dFv, Vec3D *dMv, 
					 int hydro=0, VecSet< DistSVec<double,3> > *mX=0, Vec<double> *genCF=0);

  void computeDerivativeOfForceAndMoment(dRdXoperators<dim> *, DistSVec<double,3> &,
                                         DistSVec<double,dim> &, double [3], DistSVec<double,3> &,
                                         Vec3D *, Vec3D *, Vec3D *, Vec3D *, int = 0);
  void computeTransposeDerivativeOfForceAndMoment(dRdXoperators<dim> *p, SVec<double,3> &,
                                                SVec<double,3> &, SVec<double,3> &, SVec<double,3> &,
                                                DistSVec<double,3> &, DistSVec<double,dim> &,
                                                SVec<double,3> &, DistSVec<double,3> &, int = 0);

  void computeDerivativeOperatorsOfForceAndMoment(dRdXoperators<dim> &dRdXop,
                                                  Vec3D &x0, DistSVec<double,3> &X,
                                                  DistSVec<double,dim> &U, int hydro);


  void computeDerivativeOfNodalForce(DistSVec<double,3> &, DistSVec<double,3> &,
                                     DistSVec<double,dim> &, DistSVec<double,dim> &,
                                     DistVec<double> &, double [3], DistSVec<double,3> &);

  void computeDerivativeOfNodalForce(RectangularSparseMat<double,3,3> **dForcedX,
                                     RectangularSparseMat<double,3,3> **dForcedGradP,
                                     RectangularSparseMat<double,dim,3> **dForcedV,
                                     RectangularSparseMat<double,3,3> **dForcedS,
                                     RectangularSparseMat<double,dim,dim> **dVdU,
                                     DistSVec<double,3> &dX, DistSVec<double,3> &dGradPSVec, DistSVec<double,dim> &dU,
                                     double dS[3], DistSVec<double,3> &dF);

  void computeTransposeDerivativeOfNodalForce(RectangularSparseMat<double,3,3> **dForcedX,
                                     RectangularSparseMat<double,3,3> **dForcedGradP,
                                     RectangularSparseMat<double,dim,3> **dForcedV,
                                     RectangularSparseMat<double,3,3> **dForcedS,
                                     RectangularSparseMat<double,dim,dim> **dVdU,
                                     DistSVec<double,3> &dF, DistSVec<double,3> &dX, 
                                     DistSVec<double,3> &dGradPSVec, DistSVec<double,dim> &dU,
                                     double dS[3]);

  void computeDerivativeOperatorsOfNodalForce(DistSVec<double,3> &X, DistSVec<double,dim> &U, DistVec<double> &Pin, 
                                              RectangularSparseMat<double,3,3> **dForcedX,
                                              RectangularSparseMat<double,3,3> **dForcedGradP,
                                              RectangularSparseMat<double,dim,3> **dForcedV,
                                              RectangularSparseMat<double,3,3> **dForcedS,
                                              RectangularSparseMat<double,dim,dim> **dVdU,
                                              RectangularSparseMat<double,1,dim> **dVdPstiff);
								
  void rstVar(IoData &iod) {pressInfty = iod.aero.pressure;}								

  void rstVarPostFcn(IoData &ioData) {postFcn->rstVar(ioData, com);}								

  void computeDerivativeOfNodalHeatPower(DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &, double [3], DistVec<double> &);

  void checkVec(DistSVec<double,3> &);

  void setForceGenerator(ForceGenerator<dim> *fg) { forceGen = fg; }
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <PostOperator.C>
#endif

#endif

