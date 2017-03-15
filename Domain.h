#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include <IoData.h>
#include <DistInfo.h>
#include <Timer.h>
#include <VectorSet.h>
#include <Vector.h>
#include <DenseMatrix.h>
#include <GhostPoint.h>
#include <complex>
#include <HigherOrderMultiFluid.h>
#include <LevelSet/LevelSetStructure.h>
typedef std::complex<double> bcomp;
#include <iostream>
#include <ErrorHandler.h>
#include <boost/unordered_map.hpp>
using std::cout;
using std::endl;

class IoData;
class VarFcn;
class BcFcn;
class RecFcn;
class FemEquationTerm;
class VolumicForceTerm;
class FluxFcn;
class TimeData;
class SubDTopo;
class SubDomain;
class GeoSource;
class DistGeoState;
class DistMacroCellSet;
class DynamicVMSTerm;
class VMSLESTerm;
class SmagorinskyLESTerm;
class WaleLESTerm;
class DynamicLESTerm;
class ViscoFcn;
class PostFcn;
class TimeLowMachPrec;
class SpatialLowMachPrec;
class DistLevelSetStructure;
class FluidSelector;

class BCApplier; //HB
class MatchNodeSet;

struct Vec3D;

template<class VecType> class VecSet;
template<typename Scalar> class GenFullM;
typedef GenFullM<double> FullM;

template<int dimLS> class LevelSet;
template<int dim> class RecFcnLtdMultiDim;
template<int dim> class DistEdgeGrad;
template<int dim> class DistExtrapolation;
template<int dim> class DistBcData;
template<class T> class CommPattern;
template<class Scalar> class DistVec;
template<class Scalar, int dim> class DistSVec;
template<class Scalar, int dim> class DistMat;
template<int dim> class DistExactRiemannSolver;
template<int dim> struct dRdXoperators;
template<class Scalar, int dim, int dim2> class RectangularSparseMat;

#ifndef _DNDGRAD_TMPL_
#define _DNDGRAD_TMPL_
template<int dim, class Scalar = double> class DistNodalGrad;
#endif

template<class Scalar>
class operAdd {
  public:
   static inline Scalar apply(Scalar a, Scalar b) { return a+b; }
};

template<class Scalar>
class operMin {
  public:
   static inline Scalar apply(Scalar a, Scalar b) { return std::min(a,b); }
};

template<class Scalar>
class operMax {
  public:
   static inline Scalar apply(Scalar a, Scalar b) { return std::max(a,b); }
};

//------------------------------------------------------------------------------
/** \brief Class of all data for this MPI process
 *  Class containing the geometry of all the subdomains used in the current MPI process
 *
 * The Domain class contains all the SubDomain s and offers entry points for all
 * parallel routines, distributing the work to each SubDomain.
 */
class Domain {

  int numLocSub;

  SubDomain **subDomain;
  SubDTopo *subTopo;

  int **nodeType;

  int **nodeFaceType;

  int **offWallNode;

  DistInfo *nodeDistInfo;
  DistInfo *edgeDistInfo;
  DistInfo *edgeDistInfoMF;
  DistInfo *faceDistInfo;
  DistInfo *faceNormDistInfo;
  DistInfo *inletNodeDistInfo;
  DistInfo *kirchhoffNodeDistInfo;
  DistInfo *sampledNodeDistInfo;
  DistInfo *oldSampledNodeDistInfo;

  CommPattern<double> *vecPat;
  CommPattern<double> *phiVecPat;
  CommPattern<bcomp> *compVecPat;
  CommPattern<double> *vec3DPat;
  CommPattern<double> *volPat;
  CommPattern<int> *levelPat;
  CommPattern<bool> *bool2Pat;
  CommPattern<bool> *bool3Pat;
  CommPattern<bool> *bool4Pat; //KW(02/2014): for FluidIdFS2

  CommPattern<double> *weightPat;
  CommPattern<double> *weightPhaseChangePat;
  CommPattern<double> *edgePat;
  CommPattern<double> *scalarEdgePat;
  CommPattern<double> *momPat;
  CommPattern<double> *csPat;
  CommPattern<double> *engPat;
  CommPattern<int> *fsPat;

  CommPattern<double> *inletVec3DPat;
  CommPattern<int> *inletCountPat;
  CommPattern<double> *inletRhsPat;

  Communicator *com;
  Communicator *strCom;
  Communicator *heatCom;
  Communicator *globCom;
  Communicator *embedCom;

  Timer *timer;
  Timer *strTimer;
  Timer *heatTimer;

  ErrorHandler *errorHandler;

  DistVec<double> *Delta;
  DistVec<double> *CsDelSq;
  DistVec<double> *PrT;
  DistVec<double> *WCsDelSq;
  DistVec<double> *WPrT;

  DistSVec<int,2> *tag;
  DistSVec<int,2> *tagBar;

  BCApplier* meshMotionBCs; //HB

// Included (MB)
  CommPattern<double> *weightDerivativePat;

	int numGlobNode;
	void computeNumGlobNode();
	

  Connectivity* mySubToSub;
  
  GeoSource* pGeoSource;

  // For Multifluid interface tracking
  class TriangulatedInterface* multiFluidInterface;
  
  
  // for outputting ROM snapshots from FOM
  // TODO: move these to a more appropriate location (KMW)
  int outputTimeIt;
  int outputNewtonIt;
  double outputNewtonTag;
  int outputNewtonStateStep;
  int outputNewtonResidualStep;
  int outputKrylovStep;  
  int numKrylovVecsOutputPrevNewtonIt;
  int numResidualsOutputCurrentNewtonIt;

public:

  Domain();
  Domain(Communicator *com);
  ~Domain();

  int getNumLocSub() const { return numLocSub; }
  int getNumGlobNode(); // compute if needed
  SubDomain **getSubDomain() const { return subDomain; }
  SubDTopo *getSubTopo() const { return subTopo; }
  int **getNodeType() const { return nodeType; }
  int **getNodeFaceType() const { return nodeFaceType; }
  int **getNodeTypeExtrapolation() const{
    if (inletRhsPat) return nodeType;
    int** empty = 0;
    return empty;
  }
 
  Connectivity* getSubToSub() { return mySubToSub; } 
  
  // for outputting ROM snapshots from FOM
  // TODO: move these to a more appropriate location (KMW)
  int *getTimeIt() { return &outputTimeIt; }
  int *getNewtonIt() { return &outputNewtonIt; }
  double *getNewtonTag() { return &outputNewtonTag; } 
  int *getNewtonStateStep() { return &outputNewtonStateStep; }
  int *getNewtonResidualStep() { return &outputNewtonResidualStep; }
  int *getKrylovStep() { return &outputKrylovStep; }
  int *getNumKrylovVecsOutputPrevNewtonIt() { return &numKrylovVecsOutputPrevNewtonIt; }
  int *getNumResidualsOutputCurrentNewtonIt() { return &numResidualsOutputCurrentNewtonIt; } 

  BCApplier* getMeshMotionBCs() const { return meshMotionBCs; } //HB
  CommPattern<double> *getVecPat() const { return vecPat; }
  CommPattern<bcomp> *getCompVecPat() const { return compVecPat; }
  CommPattern<double> *getVec3DPat() const { return vec3DPat; }
  CommPattern<double> *getVolPat() const { return volPat; }
  CommPattern<int> *getLevelPat() const { return levelPat; }
  CommPattern<double> *getWeightPat() const { return weightPat; }
  CommPattern<double> *getWeightPhaseChangePat() const {return weightPhaseChangePat; }
  CommPattern<double> *getMomPat() const { return momPat; }
  CommPattern<double> *getCsPat() const { return csPat; }
  CommPattern<double> *getEngPat() const { return engPat; }
  CommPattern<int> *getFsPat() const { return fsPat; }


  template<int dim>
  CommPattern<double> *getCommPat(DistSVec<double,dim> &vec) { return vecPat; }
  template<int dim>
  CommPattern<bcomp> *getCommPat(DistSVec<bcomp,dim> &vec) { return compVecPat; }
  CommPattern<double> *getCommPat(DistVec<double> &vec) { return volPat; }

  Communicator *getCommunicator() const { return com; }
  Communicator *getStrCommunicator() { return strCom; }
  Communicator *getHeatCommunicator() { return heatCom; }
  Communicator *getEmbedCommunicator() { return embedCom; }
  Timer *getTimer() const { return timer; }
  Timer *getStrTimer() const { return strTimer; }
  Timer *getHeatTimer() const { return heatTimer; }

  DistInfo &getNodeDistInfo() const { return *nodeDistInfo; }
  DistInfo &getEdgeDistInfo() const { return *edgeDistInfo; }
  DistInfo &getEdgeDistInfoMF() const { return *edgeDistInfoMF; }
  DistInfo &getFaceDistInfo() const { return *faceDistInfo; }
  DistInfo &getFaceNormDistInfo() const { return *faceNormDistInfo; }
  DistInfo &getInletNodeDistInfo() const { return *inletNodeDistInfo; }
  DistInfo &getKirchhoffNodeDistInfo() const { return *kirchhoffNodeDistInfo; }
  DistInfo &getSampledNodeDistInfo() const { return *sampledNodeDistInfo; }
  DistInfo &getOldSampledNodeDistInfo() const { return *oldSampledNodeDistInfo; }

  ErrorHandler *getErrorHandler() const {return errorHandler;}

  GeoSource& getGeoSource() { return *pGeoSource; }

  void attachTriangulatedInterface(TriangulatedInterface*);

  void getGeometry(GeoSource &, IoData&);
  void makeSampledNodeDistInfo(const std::vector<int> &cpuSample, const std::vector<int> &locSubSample);
  void makeSampledNodeDistInfo(const std::vector<int> &globalSampleNodesUnion, const boost::unordered_map<int, int> &globalNodeToCpuMap,
                               const boost::unordered_map<int, int> &globalNodeToLocSubDomainsMap);
  void makeOldSampledNodeDistInfo(const std::vector<int> &cpuSample, const std::vector<int> &locSubSample);
  void makeOldSampledNodeDistInfo(const std::vector<int> &globalSampleNodesUnion, const boost::unordered_map<int, int> &globalNodeToCpuMap,
                               const boost::unordered_map<int, int> &globalNodeToLocSubDomainsMap);
  void createRhsPat(int, IoData&);
  void createVecPat(int, IoData * = 0);
  void createPhiVecPat(int, IoData * = 0);
  void numberEdges();
  void setNodeType(IoData &);
  void setInletNodes(IoData &);
  void computeOffWallNode(DistLevelSetStructure *distLSS=0);
  void makeRotationOwnership(IoData &);
  void setFaceToElementConnectivity();
  void printElementStatistics();
  int computeControlVolumes(double, DistSVec<double,3> &, DistVec<double> &);
  void computeFaceNormals(DistSVec<double,3> &, DistVec<Vec3D> &);
  void computeNormalsConfig(DistSVec<double,3> &Xconfig, DistSVec<double,3> &Xdot,
                            DistVec<Vec3D> &edgeNorm, DistVec<double> &edgeNormVel,
                            DistVec<Vec3D> &faceNorm, DistVec<double> &faceNormVel,
                            bool communicate = true);
  void computeNormalsGCL1(DistSVec<double,3> &, DistSVec<double,3> &,
			  DistSVec<double,3> &, DistVec<Vec3D> &, DistVec<double> &,
			  DistVec<Vec3D> &, DistVec<double> &);
  void computeNormalsGCL2(TimeData &, DistVec<Vec3D> &, DistVec<Vec3D> &,
			  DistVec<double> &, DistVec<double> &, DistVec<Vec3D> &,
			  DistVec<Vec3D> &, DistVec<double> &, DistVec<double> &);
  void computeNormalsEZGCL1(double, DistSVec<double,3>&, DistSVec<double,3>&, DistVec<Vec3D>&,
			    DistVec<double>&, DistVec<Vec3D>&, DistVec<double>&);
  void computeNormalsEZGCL2(TimeData&, DistVec<double>&, DistVec<double>&,
			    DistVec<double>&, DistVec<double>&);
  void computeNormalsEZGCL3(TimeData&, DistVec<double>&, DistVec<double>&, DistVec<double>&,
			    DistVec<double>&, DistVec<double>&, DistVec<double>&);
  void computeInletNormals(DistVec<Vec3D>&, DistVec<Vec3D>&, DistVec<int> &);
  void computeVelocities(DGCLData::Velocities, TimeData &, DistSVec<double,3> &,
			 DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,3> &,
			 DistSVec<double,3> &);

  void computeWeightsLeastSquares(DistSVec<double,3> &, DistSVec<double,6> &);
  void computeWeightsLeastSquares(DistSVec<double,3> &, const DistVec<int>&, DistSVec<double,6> &, 
				  DistLevelSetStructure* =0, bool includeSweptNodes = true);
  void computeWeightsLeastSquares(DistSVec<double,3> &, const DistVec<int>&, DistSVec<double,6> &,
				  DistVec<int> &, DistVec<int> &, DistLevelSetStructure* =0);

  void computeWeightsGalerkin(DistSVec<double,3> &, 
			      DistSVec<double,3> &,
			      DistSVec<double,3> &, 
			      DistSVec<double,3> &);
  void computeWeightsGalerkin(DistSVec<double,3> &, const DistVec<int>&,
			      DistSVec<double,3> &,
			      DistSVec<double,3> &, 
			      DistSVec<double,3> &,
			      DistLevelSetStructure* =0, bool includeSweptNodes = true); //d2d


  void getReferenceMeshPosition(DistSVec<double,3> &);
  void getNdAeroLists(int *&, int **&, int *&, int **&, int *&, int **&, MatchNodeSet** matchNodes=0);

  void computeDelRatios(DistMacroCellSet *, DistVec<double> &, int, double *, double *, double *, double *, int *);
  void applySmoothing(DistVec<double> &, DistVec<double> &);
  void applySmoothing(DistVec<double> &, DistSVec<double,2> &);
  void computeTetsConnectedToNode(DistVec<int> &);
  void outputCsDynamicLES(DynamicLESTerm *, DistVec<double> &, DistSVec<double,2> &,
                          DistSVec<double,3> &, DistVec<double> &);
  void findNodeBoundingBoxes(DistSVec<double,3> &X, DistSVec<double,3> &Xmin, DistSVec<double,3> &Xmax);

  template<class MatScalar, class PrecScalar>
  void computeStiffAndForce(DefoMeshMotionData::Element, DistSVec<double,3>&,
			    DistSVec<double,3>&, DistMat<MatScalar,3>&,
			    DistMat<PrecScalar,3>*, double volStiff, int** ndType = 0);

  template<int dim>
  void computeTimeStep(double, double, double, FemEquationTerm *, VarFcn *, DistGeoState &,
		       DistSVec<double,3> &, DistVec<double> &,
		       DistSVec<double,dim> &, DistVec<double> &, DistVec<double> &,
		       DistVec<double> &, DistVec<double> &, DistVec<double> &,
                       TimeLowMachPrec &, SpatialLowMachPrec &, DistLevelSetStructure *distLSS=0);

  template<int dim>
  void computeTimeStep(double, double, double, FemEquationTerm *, VarFcn *, DistGeoState &, DistVec<double> &,
                       DistSVec<double,dim> &, DistVec<double> &, DistVec<double> &,
		       DistVec<double> &, DistVec<double> &, TimeLowMachPrec &,
		       DistVec<int> &, DistVec<double>* = NULL);


  template<int dim, class Scalar>
  void computeGradientsLeastSquares(DistSVec<double,3> &, DistSVec<double,6> &,
				    DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &,
				    DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &);

  template<int dim, class Scalar>
  void computeGradientsLeastSquares(DistSVec<double,3> &, DistVec<int> &,
                                    DistSVec<double,6> &,
                                    DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &,
                                    DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &, bool linFSI = true,
                                    DistLevelSetStructure* =0,bool includeSweptNodes = true);

  template<class Scalar>
  void computeGradientLeastSquares(DistSVec<double,3> &, DistVec<int> &,
                                    DistSVec<double,6> &,
                                    DistVec<Scalar> &, DistVec<Scalar> &,
                                    DistVec<Scalar> &, DistVec<Scalar> &,
                                    DistLevelSetStructure* =0);

  template<int dim, class Scalar>
  void computeGradientsLeastSquares(DistSVec<double,3> &, DistVec<int> &,
                                    DistSVec<double,6> &,
                                    DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &,
									DistSVec<Scalar,dim> &, DistVec<int> &, 
									DistVec<int> &, DistSVec<Scalar,dim> &,
                                    DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &,
									bool linFSI = true, DistLevelSetStructure* =0); 


  template<int dim, class Scalar>
  void computeGradientsGalerkin(DistVec<double> &, DistSVec<double,3> &,
				DistSVec<double,3> &, DistSVec<double,3> &,
				DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &,
				DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &);

  template<int dim, class Scalar>
  void computeGradientsGalerkinT(DistVec<double> &, DistSVec<double,3> &,
                DistSVec<double,3> &, DistSVec<double,3> &,
                DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &,
                DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &,
                DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &);

  template<int dim>
  void computePointWiseSourceTerm(DistGeoState &geoState, DistVec<double> &ctrlVol,
				  DistNodalGrad<dim> &dVdxj, DistSVec<double,dim> &VV,
				  DistSVec<double,dim> &RR);
  template<int dim>
  void computeMultiDimLimiter(RecFcnLtdMultiDim<dim> *, DistSVec<double,3> &,
			      DistVec<double> &, DistSVec<double,dim> &,
			      DistSVec<double,dim> &, DistSVec<double,dim> &,
			      DistSVec<double,dim> &, DistSVec<double,dim> &,
			      DistSVec<double,dim> &, DistSVec<double,dim> &);
  template<int dim>
  void computeMultiDimLimiter(RecFcnLtdMultiDim<dim> *, DistSVec<double,3> &,
			      DistVec<double> &, DistSVec<bcomp,dim> &,
			      DistSVec<bcomp,dim> &, DistSVec<bcomp,dim> &,
			      DistSVec<bcomp,dim> &, DistSVec<bcomp,dim> &,
			      DistSVec<bcomp,dim> &, DistSVec<bcomp,dim> &)  {
       std::cout << "computeMultiDimLimiter not implemented for complex data" << endl; }

  template<int dim>
	  void computeMultiDimLimiter(DistSVec<double,3> &X,
											DistVec<double> &A, DistSVec<double,dim> &V,
											DistSVec<double,dim> &dVdx, 
											DistSVec<double,dim> &dVdy,
											DistSVec<double,dim> &dVdz,
											DistLevelSetStructure *distLSS);

  template<int dim>
  void computePressureSensor(double, DistSVec<double,3>&,
                             DistSVec<double,dim>&, DistSVec<double,dim>&,
                             DistSVec<double,dim>&, DistSVec<double,dim>&,
                             DistSVec<double,3>&, DistVec<double>&);
  template<int dim>
  void computePressureSensor(double, DistSVec<double,3>&,
                             DistSVec<bcomp,dim>&, DistSVec<bcomp,dim>&,
                             DistSVec<bcomp,dim>&, DistSVec<bcomp,dim>&,
                             DistSVec<bcomp,3>&, DistVec<bcomp>&) {
	 std::cout << "computePressureSensor not implemented for complex operations" <<endl; }

  // ----- BEGIN LEVELSET - MULTIPHASE FLOW SPECIFIC FUNCTIONS ----- //
  template<int dimLS>
  void setupPhiVolumesInitialConditions(const int volid, const int fluidId, DistSVec<double,dimLS> &Phi);
  template<int dimLS>
  void TagInterfaceNodes(int lsdim, DistVec<int> &Tag, DistSVec<double,dimLS> &Phi, int level,DistLevelSetStructure *distLSS=0);
  template<int dimLS>
  void TagInterfaceNodes(int lsdim, DistSVec<bool,2> &Tag, DistSVec<double,dimLS> &Phi, DistLevelSetStructure *distLSS);
  template<int dimLS>
  void pseudoFastMarchingMethod(DistVec<int> &Tag, DistSVec<double,3> &X, 
				DistSVec<double,dimLS> &d2wall, int level, int iterativeLevel,
				DistVec<int> &sortedNodes, int *nSortedNodes,
				int *firstCheckedNode,DistLevelSetStructure *distLSS=0);
  //template<int dimLS>
  //void FinishReinitialization(DistVec<int> &Tag, DistSVec<double,dimLS> &Psi, int level);

  template<int dim>
  void setupUVolumesInitialConditions(const int volid, double UU[dim], DistSVec<double,dim> &U);

  void setupFluidIdVolumesInitialConditions(const int volid, const int myId, DistVec<int> &fluidId);

  //template<int dim>
  //void setupUMultiFluidInitialConditionsSphere(FluidModelData &fm,
  //           SphereData &ic, DistSVec<double,3> &X, DistSVec<double,dim> &U);
  //template<int dim>
  //void setupUMultiFluidInitialConditionsPlane(FluidModelData &fm,
  //           PlaneData &ip, DistSVec<double,3> &X, DistSVec<double,dim> &U);
  //template<int dim>
  //void setupUMultiFluidInitialConditionsPlane(FluidModelData &fm,
  //           PlaneData &ip, DistSVec<double,3> &X, DistSVec<double,dim> &U, DistVec<int> &nodeTag);

  template<int dim>
  void storeGhost(DistSVec<double,dim> &, DistSVec<double,dim> &, DistVec<double> &);

  template<int dim>
  void computeWeightsForEmbeddedStruct(DistSVec<double,3> &X, DistSVec<double,dim> &V, 
                          DistVec<double> &Weights, DistSVec<double,dim> &VWeights, 
                          DistVec<int> &init, DistVec<int> &next_init,
													DistLevelSetStructure *distLSS, bool externalSI=false);

  template<int dim>
  void computeWeightsForFluidFluid(DistSVec<double,3> &X, DistSVec<double,dim> &V, 
                          DistVec<double> &Weights, DistSVec<double,dim> &VWeights, 
                          DistVec<int> &init, DistVec<int> &next_init,
				   DistLevelSetStructure *distLSS,DistVec<int>& fid);

  template<int dim>
  void computeWeightsLeastSquaresForEmbeddedStruct(DistSVec<double,3> &X, DistSVec<double,dim> &V, 
                          DistVec<double> &Weights, DistSVec<double,dim> &VWeights, 
                          DistVec<int> &init, DistVec<int> &next_init,
                          DistLevelSetStructure *distLSS,
			   DistNodalGrad<dim>& dX, bool limit,
																	DistVec<int>* fid, bool externalSI=false);

  template<int dim>
  void computeWeightsLeastSquaresForFluidFluid(DistSVec<double,3> &X, DistSVec<double,dim> &V, 
                          DistVec<double> &Weights, DistSVec<double,dim> &VWeights, 
                          DistVec<int> &init, DistVec<int> &next_init,
					       DistLevelSetStructure *distLSS,DistVec<int>& fid,
					       DistNodalGrad<dim>&, bool limit);

  template<int dim, int dimLS>
  void computeWeightsForEmbeddedStruct(DistSVec<double,3> &X, DistSVec<double,dim> &V, 
                          DistVec<double> &Weights, DistSVec<double,dim> &VWeights, 
                          DistSVec<double,dimLS> &Phi, DistSVec<double,dimLS> &PhiWeights, 
                          DistVec<int> &init, DistVec<int> &next_init,
                          DistLevelSetStructure *distLSS, DistVec<int> *fluidId);
  template<int dimLS>
  void extrapolatePhiV(DistLevelSetStructure *distLSS, DistSVec<double,dimLS> &PhiV);

  template<int dim>
  void populateGhostPoints(DistVec<GhostPoint<dim>*> *ghostPoints, DistSVec<double,3> &X, 
									DistSVec<double,dim> &U, DistNodalGrad<dim, double> *ngrad, VarFcn *varFcn,
									DistLevelSetStructure *distLSS,bool linRecAtInterface, DistVec<int> &tag, 
									bool externalSI=false, FemEquationTerm *fet = 0);

  template<int dim, class Scalar, int neq>
  void populateGhostJacobian(DistVec<GhostPoint<dim>*> *ghostPoints,DistSVec<double,dim> &U,FluxFcn** fluxFcn, VarFcn *varFcn,DistLevelSetStructure *distLSS,DistVec<int> &tag, DistMat<Scalar,neq>& A);

  template<int dim>
  void setSIstencil(DistSVec<double,3> &X, DistLevelSetStructure *distLSS, DistVec<int> &fluidId, DistSVec<double,dim> &U);

  template<int dim>
  void setFEMstencil(DistSVec<double,3> &X, DistLevelSetStructure *distLSS, DistVec<int> &fluidId, DistSVec<double,dim> &U);

  template<int dim>
  void computeRiemannWeightsForEmbeddedStruct(DistSVec<double,3> &X, DistSVec<double,dim> &V,
               DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji,
               DistVec<double> &Weights, DistSVec<double,dim> &VWeights, DistLevelSetStructure *distLSS);

  template<int dim, int dimLS>
  void computeRiemannWeightsForEmbeddedStruct(DistSVec<double,3> &X, DistSVec<double,dim> &V,
               DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji,
               DistVec<double> &Weights, DistSVec<double,dim> &VWeights, 
               DistSVec<double,dimLS> &Phi, DistSVec<double,dimLS> &PhiWeights,
               DistLevelSetStructure *distLSS, DistVec<int> *fluidId0, DistVec<int> *fluidId);

  template<int dimLS>
  void computeDistanceCloseNodes(int lsdim, DistVec<int> &Tag, DistSVec<double,3> &X,
                                 DistNodalGrad<dimLS> &lsgrad,
                                 DistSVec<double,dimLS> &Phi,DistSVec<double,1> &Psi,
                                 MultiFluidData::CopyCloseNodes copy);
  template<int dimLS>
  void computeDistanceLevelNodes(int lsdim, DistVec<int> &Tag, int level,
                                 DistSVec<double,3> &X,DistSVec<double,1> &Psi,
                                 double &res, DistSVec<double,dimLS> &Phi,
                                 MultiFluidData::CopyCloseNodes copy);
  template<int dimLS>
  void checkNodePhaseChange(DistSVec<double,dimLS> &PhiProduct);
  template<int dim>
  void getSignedDistance(int lsdim, DistSVec<double,1> &Psi, DistSVec<double,dim> &Phi);
  template<int dimLS>
  void avoidNewPhaseCreation(DistSVec<double,dimLS> &Phi, DistSVec<double,dimLS> &Phin);
  template<int dimLS>
  void avoidNewPhaseCreation(DistSVec<double,dimLS> &Phi, DistSVec<double,dimLS> &Phin, DistVec<double> &weight, DistLevelSetStructure *distLSS = 0, DistVec<int>* fluidIdToSet = 0);
  template<int dim>
  void storePreviousPrimitive(DistSVec<double,dim> &V, DistVec<int> &fluidId, 
                              DistSVec<double,3> &X, DistSVec<double,dim> &Vupdate, 
                              DistVec<double> &weight);
  template<int dim>
  void IncreasePressure(double p, VarFcn *vf, DistSVec<double,dim> &U);
  template<int dim>
  void IncreasePressure(double p, VarFcn *vf, DistSVec<double,dim> &U, DistVec<int> &fluidId);

  // ----- END   LEVELSET - MULTIPHASE FLOW SPECIFIC FUNCTIONS ----- //

  template<int dim>
  void computeFiniteVolumeTerm(DistVec<double> &, DistVec<double> &, FluxFcn**, RecFcn*, DistBcData<dim>&, DistGeoState&,
			       DistSVec<double,3>&, DistSVec<double,dim>&,
			       DistNodalGrad<dim>&, DistEdgeGrad<dim>*,
			       DistSVec<double,dim>&, int, int);

  template<int dim>
  void computeFiniteVolumeTerm(DistExactRiemannSolver<dim>&,
                               DistVec<double> &, DistVec<double> &, FluxFcn**, RecFcn*, DistBcData<dim>&, DistGeoState&,
			       DistSVec<double,3>&, DistSVec<double,dim>&,
			       DistNodalGrad<dim>&, DistEdgeGrad<dim>*,
			       DistSVec<double,dim>&, int, int);

  template<int dim, int dimLS>
  void computeFiniteVolumeTerm(DistVec<double> &, 
										 DistExactRiemannSolver<dim> &,
                               FluxFcn**, RecFcn*, 
										 DistBcData<dim>&, 
										 DistGeoState&,
                               DistSVec<double,3>&, 
										 DistSVec<double,dim>&,
                               FluidSelector &,
                               DistNodalGrad<dim>&, DistEdgeGrad<dim>*,
			       DistSVec<double,dimLS>& phi,
                               DistNodalGrad<dimLS>&, 
			       DistEdgeGrad<dimLS>*,
			       DistSVec<double,dim>&,
                               int, int, int);

  // for multi-phase fluid-structure interaction under embedded framework 
  template<int dim, int dimLS>
  void computeFiniteVolumeTerm(DistVec<double> &ctrlVol, 
										 DistExactRiemannSolver<dim>&riemann, 
                               FluxFcn** fluxFcn, RecFcn* recFcn, 
										 DistBcData<dim>& bcData, 
										 DistGeoState& geoState,
                               DistSVec<double,3>& X, 
										 DistSVec<double,dim>& V, 
										 DistSVec<double,dim>& Wstarij, DistSVec<double,dim>& Wstarji, 
										 DistLevelSetStructure* LSS, bool linRecAtInterface, 
										 FluidSelector &fluidSelector, int Nriemann,
										 DistNodalGrad<dim>& ngrad,   DistEdgeGrad<dim>* egrad, 
										 DistSVec<double,dimLS>&  phi,
                               DistNodalGrad<dimLS>& ngradLS, DistEdgeGrad<dimLS>* egradLS,
										 DistSVec<double,dim>& R, int it, int failsave, int rshift);

  // embedded structure
  template<int dim>
  void computeFiniteVolumeTerm(DistVec<double> &ctrlVol, 
										 DistExactRiemannSolver<dim>&riemann,
                               FluxFcn** fluxFcn, RecFcn* recFcn, 
										 DistBcData<dim>& bcData, 
										 DistGeoState& geoState,
                               DistSVec<double,3>& X, 
										 DistSVec<double,dim>& V,
                               DistSVec<double,dim>& Wstarij, DistSVec<double,dim>& Wstarji,
										 DistSVec<double,dim>& Wext,
                               DistLevelSetStructure *LSS, bool linRecAtInterface, DistVec<int> &fluidId, 
										 int Nriemann,
                               DistNodalGrad<dim>& ngrad, DistEdgeGrad<dim>* egrad,
                               DistSVec<double,dim>& R,
										 int it, int failsafe, int rshift, bool externalSI=false);

  template<int dim>
  void computeFiniteVolumeTerm(DistVec<double> &, DistExactRiemannSolver<dim>&,
                               FluxFcn**, RecFcn*, DistBcData<dim>&, DistGeoState&,
                               DistSVec<double,3>&, DistSVec<double,dim>&,
                               DistSVec<double,dim>&, DistSVec<double,dim>&, 
                               DistVec<int> &, DistVec<int> &, DistLevelSetStructure *, 
                               bool, DistVec<int> &, int,
                               double, double, DistNodalGrad<dim>&, DistEdgeGrad<dim>*,
                               DistSVec<double,dim>&, int, int, int); 

  template<int dim, int dimLS>
  void computeFiniteVolumeTermLS(FluxFcn**, RecFcn*, RecFcn*, DistBcData<dim>&, DistGeoState&,
                               DistSVec<double,3>&, DistSVec<double,dim>&,DistVec<int>& fluidId,
                               DistNodalGrad<dim>&,   DistEdgeGrad<dim>*,
			       DistNodalGrad<dimLS>&, DistEdgeGrad<dimLS>*,
                               DistSVec<double,dimLS> &, DistSVec<double,dimLS> &,
                               DistLevelSetStructure * = 0, int lsorder = 1);

  template<int dim>
  void computeFiniteVolumeBarTerm(DistVec<double> &, DistVec<double> &, FluxFcn**,
                                  RecFcn*, DistBcData<dim>&, DistGeoState&, DistSVec<double,3>&,
                                  DistMacroCellSet *, DistSVec<double,dim> &, DistSVec<double,1> &,
                                  DistNodalGrad<dim>&, DistEdgeGrad<dim>*,
                                  DistSVec<double,dim>&, int, int, int, int);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, DistBcData<dim> &, DistGeoState &,
				       DistVec<double> &,
				       DistSVec<double,3> &,
				       DistVec<double> &, DistSVec<double,dim> &,
				       DistMat<Scalar,neq> &);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(DistExactRiemannSolver<dim> &,
                                       FluxFcn **, DistBcData<dim> &, DistGeoState &,
				       DistVec<double> &,
				       DistSVec<double,3> &,
				       DistVec<double> &, DistSVec<double,dim> &,
				       DistMat<Scalar,neq> &);

  template<int dim, class Scalar, int neq, int dimLS>
  void computeJacobianFiniteVolumeTerm(DistExactRiemannSolver<dim> &,
                                       FluxFcn **, DistBcData<dim> &, DistGeoState &,
                                       DistNodalGrad<dim>&, DistNodalGrad<dimLS>&,
                                       DistSVec<double,3> &,
                                       DistVec<double> &, DistSVec<double,dim> &,
                                       DistMat<Scalar,neq> &, FluidSelector &);
  template<class Scalar,int dim,int neq>
  void computeJacobianFiniteVolumeTerm(DistVec<double> &ctrlVol,
                                       DistExactRiemannSolver<dim> &riemann,
                                       FluxFcn** fluxFcn,
                                       DistBcData<dim>& bcData, DistGeoState& geoState,
                                       DistSVec<double,3>& X, DistSVec<double,dim>& V,
                                       DistLevelSetStructure *LSS, DistVec<int> &fluidId, 
                                       int Nriemann,
                                       DistMat<Scalar,neq>& A,DistVec<double>& irey, bool externalSI=false);
  
  template<int dim, class Scalar, int neq, int dimLS>
  void computeJacobianFiniteVolumeTerm(DistExactRiemannSolver<dim>& riemann,
                                       FluxFcn** fluxFcn, 
                                       DistBcData<dim>& bcData, DistGeoState& geoState,
                                       DistSVec<double,3>& X, DistSVec<double,dim>& V,DistVec<double>& ctrlVol,
                                       DistNodalGrad<dimLS> &ngradLS,
                                       DistLevelSetStructure *LSS,
                                       int Nriemann,
                                       FluidSelector &fluidSelector,
                                       DistMat<Scalar,neq>& A);

   template<int dim, class Scalar, int dimLS>
    void computeJacobianFiniteVolumeTermLS(RecFcn* recFcn, RecFcn* recFcnLS,
					   DistGeoState &geoState,DistSVec<double,3>& X,DistSVec<double,dim> &V,
					   DistNodalGrad<dim>& ngrad,DistNodalGrad<dimLS> &ngradLS,
					   DistEdgeGrad<dim>* egrad,
					   DistVec<double> &ctrlVol,DistSVec<double,dimLS>& Phi,
					   DistMat<Scalar,dimLS> &A,DistLevelSetStructure* distLSS);
 
  template<int dim>
  void recomputeRHS(VarFcn*, DistSVec<double,dim> &, DistSVec<double,dim> &, DistExtrapolation<dim>*,
                   DistBcData<dim>&, DistGeoState&, DistSVec<double,3> &);
  template<int dim>
  void recomputeRHS(VarFcn*, DistSVec<double,dim> &, DistVec<int> &,
                   DistSVec<double,dim> &, DistExtrapolation<dim>*,
                   DistBcData<dim>&, DistGeoState&, DistSVec<double,3> &);
  template<int dim>
  double recomputeResidual(DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<int dim>
  double computeRealFluidResidual(DistSVec<double, dim> &, DistSVec<double,dim> &, DistLevelSetStructure &);

  template<int dim>
  void computeGalerkinTerm(FemEquationTerm *, DistBcData<dim> &,
			   DistGeoState &, DistSVec<double,3> &,
			   DistSVec<double,dim> &, DistSVec<double,dim> &,
									DistVec<GhostPoint<dim>*> *ghostPoints=0, 
									DistLevelSetStructure *LSS=0, 	
									bool externalSI=false);

  template<int dim>
  void computeVolumicForceTerm(VolumicForceTerm *, DistVec<double> &,
                               DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<int dim>
  void computeSmagorinskyLESTerm(SmagorinskyLESTerm *, DistSVec<double,3> &,
				 DistSVec<double,dim> &, DistSVec<double,dim> &,
			         DistVec<GhostPoint<dim>*> *ghostPoints=0, 
                                 DistLevelSetStructure *LSS=0, bool externalSI=false);

  template<int dim>
  void computeWaleLESTerm(WaleLESTerm *, DistSVec<double,3> &, DistSVec<double,dim> &, 
                          DistSVec<double,dim> &, DistVec<GhostPoint<dim>*> *ghostPoints=0, 
                          DistLevelSetStructure *LSS=0, bool externalSI=false);


  //---start computation of MutOMu terms

  template<int dim>
  void computeMutOMuSmag(SmagorinskyLESTerm *, DistVec<double> &, DistSVec<double,3> &,
                               DistSVec<double,dim> &, DistVec<double> &);

  template<int dim>
  void computeMutOMuVMS(VMSLESTerm *, DistMacroCellSet *, DistVec<double> &, bool,
                        DistSVec<double,dim> &, DistSVec<double,1> &, DistSVec<double,3> &,
                        DistSVec<double,dim> &, int, DistVec<double> &);

  template<int dim>
  void computeMutOMuDynamicVMS(DynamicVMSTerm *, DistVec<double> &, DistSVec<double,dim> &,
                               DistSVec<double,3> &, DistSVec<double,dim> &, DistVec<double> &, DistVec<double> &);

  template<int dim>
  void computeMutOMuWale(WaleLESTerm *, DistVec<double> &, DistSVec<double,3> &,
                         DistSVec<double,dim> &, DistVec<double> &);

  template<int dim>
  void computeMutOMuDynamicLES(DynamicLESTerm *, DistVec<double> &, DistSVec<double,2> &,
                               DistSVec<double,3> &, DistSVec<double,dim> &, DistVec<double> &);

  //---complete computaton of MutOMu terms


  template<int dim>
  void computeGalerkinBarTerm(bool, FemEquationTerm *, DistBcData<dim> &, DistGeoState &, DistSVec<double,3> &,
                              DistMacroCellSet *, DistSVec<double,dim> &, DistSVec<double,1> &,
                              DistSVec<double,dim> &, int, int);

  template<int dim>
  void computeVBar(DistMacroCellSet *, bool, DistGeoState &, DistSVec<double,dim> &,
                   DistSVec<double,dim> &, int, int);

  template<int dim>
  void computeBarTerm(DistMacroCellSet *, bool, DistVec<double> &, DistSVec<double,dim> **,
                      DistSVec<double,1> **, DistSVec<double,3> &, DistSVec<double,dim> &, int, int);

  template<int dim>
  void computeDynamicVMSTerm(DynamicVMSTerm *, DistMacroCellSet *, bool, DistVec<double> &,
                             DistSVec<double,dim> **, DistSVec<double,1> **, DistSVec<double,3> &,
                             DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
                             DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
                             DistSVec<double,dim> &, DistSVec<double,dim> &,
                             int, int, int, DistVec<double> * = 0);

  template<int dim>
  void computeVMSLESTerm(VMSLESTerm *, DistMacroCellSet *, bool,
                           DistVec<double> &, DistSVec<double,dim> &,
                           DistSVec<double,1> &, DistSVec<double,3> &,
                           DistSVec<double,dim> &, DistSVec<double,dim> &, int);




  template<int dim>
  void computeTestFilterValues(DistVec<double> &, DistSVec<double,dim> &,
                               DistSVec<double,16> &, DistSVec<double,6> &,
                               DistVec<double> &, DistSVec<double,8> &,
                               DistSVec<double,2> &, DistVec<int> &, DistBcData<dim> &,
                               DistSVec<double,3> &, DistSVec<double,dim> &, double, double,
			       DistVec<GhostPoint<dim>*> *ghostPoints=0, 
                               DistLevelSetStructure *LSS=0,
                             	 bool externalSI=false);

  template<int dim>
  void computeDynamicLESTerm(DynamicLESTerm *, DistSVec<double,2> &,
                             DistSVec<double,3> &, DistSVec<double,dim> &, 
                             DistSVec<double,dim> &,
			     DistVec<GhostPoint<dim>*> *ghostPoints=0, 
                             DistLevelSetStructure *LSS=0, bool externalSI=false);


  template<int dim>
  void assemble_dWdt(DistSVec<double, dim> &, DistSVec<double, dim> &);

  template<int dim>
  void computedWBar_dt(DistSVec<double, dim> &, DistSVec<double, dim> &, DistMacroCellSet *, DistSVec<double,1> **, int);

  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(FemEquationTerm *fet, DistBcData<dim> &bcData,
					 DistGeoState &geoState, DistSVec<double,3> &X,
					 DistVec<double> &ctrlVol, DistSVec<double,dim> &V,
					 DistMat<Scalar,neq> &A,
											  DistVec<GhostPoint<dim>*> *ghostPoints=0,
											  DistLevelSetStructure *distLSS=0, bool externalSI=false);

  template<int dim, class Scalar, int neq>
  void computeJacobianVolumicForceTerm(VolumicForceTerm *, DistVec<double> &,
                                       DistSVec<double,dim> &, DistMat<Scalar,neq> &);

  template<class Scalar, int neq>
  void finishJacobianGalerkinTerm(DistVec<double> &, DistMat<Scalar,neq> &);

  template<int dim>
  void getExtrapolationValue(DistExtrapolation<dim> *, DistSVec<double,dim>&, DistSVec<double,dim>&,
                             VarFcn*, DistBcData<dim>&, DistGeoState&, DistSVec<double,3>&);

  template<int dim>
  void applyExtrapolationToSolutionVector(DistExtrapolation<dim> *, DistSVec<double,dim>&,
                                          DistSVec<double,dim>&);

  template<int dim>
    void applyBCsToSolutionVector(BcFcn *, DistBcData<dim> &, DistSVec<double,dim> &, DistLevelSetStructure *distLSS=0);

  template<int dim>
    void applyBCsToTurbSolutionVector(BcFcn *, DistBcData<dim> &, DistSVec<double,dim> &, DistLevelSetStructure *distLSS=0);

  template<int dim>
  void applyBCsToResidual(BcFcn *, DistBcData<dim> &,
			  DistSVec<double,dim> &, DistSVec<double,dim> &, DistLevelSetStructure *distLSS=0);

  template<int dim, class Scalar, int neq>
  void applyBCsToJacobian(BcFcn *, DistBcData<dim> &,
			  DistSVec<double,dim> &, DistMat<Scalar,neq> &, DistLevelSetStructure *distLSS=0);

  template<int dim, class Scalar, int neq>
  void applyBCsToH2Jacobian(BcFcn *, DistBcData<dim> &,
			  DistSVec<double,dim> &, DistMat<Scalar,neq> &);

  template<class Scalar, int dim>
  void computeH1(FluxFcn **, DistBcData<dim> &, DistGeoState &,
                 DistVec<double> &, DistSVec<double,dim> &, DistMat<Scalar,dim> &);

  template<int dim, class Scalar, int neq>
  void computeH2(FluxFcn **, RecFcn *, 
		 DistBcData<dim> &, DistGeoState &,
		 DistSVec<double,3> &, DistSVec<double,dim> &, 
		 DistNodalGrad<dim, double> &,
		 DistMat<Scalar,neq> &, 
		 DistSVec<double,dim> &, DistSVec<double,dim> &,
		 DistSVec<double,dim> &, DistSVec<double,dim> &);
 
  template<int dim, class Scalar, int neq>
  void computeH2transpose(FluxFcn **, RecFcn *, DistBcData<dim> &, DistGeoState &,
                          DistSVec<double,3> &, DistSVec<double,dim> &, DistNodalGrad<dim, double> &,
                          DistMat<Scalar,neq> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
                          DistSVec<double,dim> &, DistSVec<double,dim> &);
 
  template<int dim, class Scalar, int neq>
  void computeH2(FluxFcn **fluxFcn, RecFcn *recFcn,
		 DistBcData<dim> &bcData, DistGeoState &geoState,
		 DistSVec<double,3> &X, DistSVec<double,dim> &V,
		 DistNodalGrad<dim, double> &ngrad,
		 DistExactRiemannSolver<dim> &riemann,
		 DistLevelSetStructure *distLSS,
		 DistVec<int> &fluidId, 
		 int Nriemann,
		 DistMat<Scalar,neq> &H2,
		 DistSVec<double,dim> &aij, DistSVec<double,dim> &aji,
		 DistSVec<double,dim> &bij, DistSVec<double,dim> &bji,
		 DistSVec<double,dim> &betaij, DistSVec<double,dim> &betaji);
  
  template<class Scalar1, class Scalar2, int dim>
  void computeMatVecProdH2(RecFcn *, DistSVec<double,3> &, DistVec<double> &,
                           DistMat<Scalar1,dim> &, DistSVec<double,dim> &,
                           DistSVec<double,dim> &, DistSVec<double,dim> &,
                           DistSVec<double,dim> &, DistSVec<Scalar2,dim> &,
                           DistNodalGrad<dim, Scalar2> &, DistSVec<Scalar2,dim> &);

  template<class Scalar1, class Scalar2, int dim>
  void computeMatVecProdH2transposeNew(IoData& iod, DistSVec<double,3> &X,
                                       DistVec<double> &ctrlVol, DistMat<Scalar1,dim> &H2,
                                       DistSVec<double,dim> &aij, DistSVec<double,dim> &aji,
                                       DistSVec<double,dim> &bij, DistSVec<double,dim> &bji,
                                       DistNodalGrad<dim, Scalar2> &dpdxj,
                                       DistSVec<Scalar2,dim> &p, DistSVec<Scalar2,dim> &prod);

  template<class Scalar1, class Scalar2, int dim>
  void computeMatVecProdH2(FluxFcn **fluxFcn, RecFcn *recFcn, DistGeoState &geoState, 
			   DistSVec<double,3> &X, DistVec<double> &ctrlVol, 
			   DistExactRiemannSolver<dim> &riemann,
			   DistLevelSetStructure *distLSS,
			   DistVec<int> &fluidId,
			   int Nriemann,
			   DistMat<Scalar1,dim> &H2,
			   DistSVec<double,dim> &aij, DistSVec<double,dim> &aji,
			   DistSVec<double,dim> &bij, DistSVec<double,dim> &bji,
			   DistSVec<double,dim> &betaij, DistSVec<double,dim> &betaji,
			   DistSVec<Scalar2,dim> &p, DistNodalGrad<dim, Scalar2> &dpdxj,
			   DistSVec<Scalar2,dim> &prod);

  template<class Scalar1, class Scalar2, int dim>
  void computeMatVecProdH2T(RecFcn *, DistSVec<double,3> &,
                DistVec<double> &, DistMat<Scalar1,dim> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &,
                DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &,
                DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &,
                DistSVec<Scalar2,dim> &);

  template<class Scalar1, class Scalar2, int dim>
  void computeMatVecProdH2Tb(RecFcn *, DistSVec<double,3> &,
                DistVec<double> &, DistMat<Scalar1,dim> &,
                DistNodalGrad<dim, Scalar2> &, DistSVec<Scalar2,dim> &,
                DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);

  template<int dim, class Scalar>
  void assemble(CommPattern<Scalar> *, DistSVec<Scalar,dim> &);

  template<int dim, class Scalar, class OpType >
  void assemble(CommPattern<Scalar> *, DistSVec<Scalar,dim> &, const OpType &);

  template<class Scalar>
  void assemble(CommPattern<Scalar> *, DistVec<Scalar> &);

  template<class Scalar, class OpType >
  void assemble(CommPattern<Scalar> *, DistVec<Scalar> &, const OpType &);

  void assemble(DistVec<double> &v) {
    assemble(getVecPat(), v);
  }

  void assemble(DistSVec<double,3> &v) {
    assemble(getVec3DPat(), v);
  }

  template<class OpType>
  void assemble(DistSVec<double,3> &v, const OpType &oper) {
    assemble(getVec3DPat(), v, oper);
  }

  void assembleEdge(CommPattern<double> *commPat, DistVec<double> &W);
  
  void readEigenValuesAndVectors(const char *eigFile, double &realEigV, double &imagEigV, double &, int &, complex<double>*&, complex<double>*&);

  template<int dim>
  void assembleGhostPoints(DistVec<GhostPoint<dim>*> &ghostPoints, VarFcn *varFcn);

  template<class Scalar, int dim>
  bool readTagFromFile(const char *, int, double *, int *);

  template<class Scalar, int dim>
  bool readVectorFromFile(const char *, int, double *, DistSVec<Scalar,dim> &, Scalar* = 0);

	template<class Scalar>
	bool readVectorFromFile(const char *, int, double *, DistVec<Scalar> &, Scalar* = 0); //<! for non-state vector, Lei Lei, 20 Mar 2016
	//template<int dim>
	//void readMultiPodBasis(const char *, VecSet< DistSVec<double, dim> > **, int *, int nBasesNeeded = 0, int *whichFiles = NULL); 	//KTC

	void computeConnectedTopology(const std::vector<std::vector<int> > & locSampleNodes);
	void computeConnectedNodes(const std::vector<std::vector<int> > &,
			std::vector<int> &);

	template<typename Scalar> void communicateMesh( std::vector <Scalar> *nodeOrEle , int arraySize, int *alreadyCommunicated = NULL);

	template<typename Scalar> void makeUnique( std::vector <Scalar> *nodeOrEle, int length = 1);

  //template<int dim>
  //void readPodBasis(const char *, int &nPod, VecSet<DistSVec<double ,dim> > &, bool snaps = false);

  void readInterpNode(const char *, int &, int *&, int *&); // for Gappy Pod

  void readInterpMatrix(const char *, int &, FullM &); // for Gappy Pod

  //void readSampleNodes(std::vector<int> &, int &, const char *);

  template<class Scalar, int dim>
  void writeVectorToFile(const char *, int, double, DistSVec<Scalar,dim> &, Scalar* = 0);

	template<class Scalar>
	void writeVectorToFile(const char *, int, double, DistVec<Scalar> &, Scalar* = 0); //<! for non-state vector, Lei Lei, 03 Feb 2016

  template<class Scalar, int dim>
  void scaleSolution(DistSVec<Scalar,dim> &, RefVal*);

  template<int dim>
  void computeForceDerivs(VarFcn *, DistSVec<double,3> &, DistSVec<double,dim> &,
                          DistSVec<double ,dim> &, Vec<double> &, VecSet< DistSVec<double,3> > &);
/*
  template<int dim>
  void computeForceCoefficients(PostFcn *, Vec3D &, DistSVec<double,3> &,
                                DistSVec<double,dim> &, Vec3D &, Vec3D &, Vec3D &, Vec3D &);
*/

  void testNormals(DistVec<Vec3D> &, DistVec<double> &, DistVec<Vec3D> &, DistVec<double> &);

  template<int dim>
  int checkSolution(VarFcn *, DistSVec<double,dim> &, DistLevelSetStructure *distLSS=0);

  template<int dim>
  int checkSolution(VarFcn *, DistSVec<double,dim> &, DistVec<int> &fluidId, DistLevelSetStructure *distLSS=0);

  template<int dim>
  int checkSolution(VarFcn *, DistVec<double> &, DistSVec<double,dim> &, DistVec<int> &, DistVec<int> &);


  template<int dim>
  void restrictionOnPhi(DistSVec<double,dim> &initial, DistVec<int> &fluidId,
                        DistSVec<double,dim> &restriction, int fluidIdTarget);

  template<int dim>
  void checkFailSafe(VarFcn*, DistSVec<double,dim>&, DistSVec<bool,2>&, DistVec<int> * = 0);

  template<int dim, int neq>
  int clipSolution(TsData::Clipping, BcsWallData::Integration,
		   VarFcn*, double*, DistSVec<double,dim>&);

  template<int dim>
  void checkGradients(DistSVec<double,3> &, DistVec<double> &,
		      DistSVec<double,dim> &, DistNodalGrad<dim> &);

  template<int dim>
  void checkMatVecProd(DistSVec<double,dim> &, const char *);
    

  template<int dim>
  void zeroInternalVals(DistSVec<double, dim> &);

  template<int dim>
  void printVariable(DistSVec<double,dim>&, VarFcn *);

  template<int dim>
  void printInletVariable(DistSVec<double,dim>&);
  template<int dim>
  void printAllVariable(DistVec<int> &, DistSVec<double,dim>&, int );
	// for printing distLSS->is_active() to stdout, Lei Lei
	void printDistVecBool(DistVec<bool> &, bool );

  template<int dimLS>
  void printPhi(DistSVec<double,3> &, DistSVec<double,dimLS> &, int);

  template<int dim>
  void checkExtrapolationValue(DistSVec<double,dim>&, VarFcn*,
                               DistBcData<dim>&, DistGeoState&);
  template<class Scalar, int neq>
  void printAllMatrix( DistMat<Scalar,neq> &, int );

  // XML Debugging routine
  void checkNormalSums(IoData &);

  template<int dim>
  void padeReconstruction(VecSet<DistSVec<double, dim> >&, VecSet<DistSVec<double, dim> >&, int*, double*, double, int, int, int, int );

  template<int dim>
  void hardyInterpolationLogMap(VecSet<DistSVec<double, dim> >**, VecSet<DistSVec<double, dim> >&, int, int, int, FullM &, FullM &);

// Included (MB)
  int computeDerivativeOfControlVolumes(double, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double> &);
  int computeDerivativeOfControlVolumes(RectangularSparseMat<double,3,1> **, DistSVec<double,3> &, DistVec<double> &);
  int computeTransposeDerivativeOfControlVolumes(RectangularSparseMat<double,3,1> **, DistVec<double> &, DistSVec<double,3> &);
  int computeDerivativeOperatorsOfControlVolumes(double, DistSVec<double,3> &, RectangularSparseMat<double,3,1> **);

  void computeDerivativeOfNormals(DistSVec<double,3> &, DistSVec<double,3> &, DistVec<Vec3D> &, DistVec<Vec3D> &,
                                  DistVec<double> &, DistVec<double> &, DistVec<Vec3D> &, DistVec<Vec3D> &, DistVec<double> &, DistVec<double> &);
  void computeDerivativeOfNormals(RectangularSparseMat<double,3,3> **, RectangularSparseMat<double,3,3> **, DistSVec<double,3> &, DistVec<Vec3D> &, 
                                  DistVec<double> &, DistVec<Vec3D> &, DistVec<double> &);
  void computeTransposeDerivativeOfNormals(RectangularSparseMat<double,3,3> **, RectangularSparseMat<double,3,3> **, DistVec<Vec3D> &, DistVec<Vec3D> &, DistSVec<double,3> &);
  void computeDerivativeOperatorsOfNormals(DistSVec<double,3> &, RectangularSparseMat<double,3,3> **, RectangularSparseMat<double,3,3> **); 

  void computeDerivativeOfWeightsLeastSquares(DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,6> &);
  void computeDerivativeOfWeightsLeastSquares(RectangularSparseMat<double,3,6> **, RectangularSparseMat<double,6,6> **, DistSVec<double,3> &, DistSVec<double,6> &);
  void computeTransposeDerivativeOfWeightsLeastSquares(RectangularSparseMat<double,3,6> **, RectangularSparseMat<double,6,6> **, DistSVec<double,6> &, DistSVec<double,3> &);
  void computeDerivativeOperatorsOfWeightsLeastSquares(DistSVec<double,3> &, RectangularSparseMat<double,3,6> **, RectangularSparseMat<double,6,6> **);
  void computeDerivativeTransposeOfWeightsLeastSquares(DistSVec<double,3> &, DistSVec<double,6> &, DistSVec<double,3> &);

  void computeDerivativeOfWeightsGalerkin(DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,3> &,
			      DistSVec<double,3> &, DistSVec<double,3> &);
  void computeDerivativeOperatorsOfWeightsGalerkin(DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,3> &,
			      DistSVec<double,3> &, DistSVec<double,3> &);
  void computeDerivativeTransposeOfWeightsGalerkin(DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,3> &,
			      DistSVec<double,3> &, DistSVec<double,3> &);

  template<int dim, class Scalar>
  void computeDerivativeOfGradientsLeastSquares(
            DistSVec<double,3> &, DistSVec<double,3> &,
            DistSVec<double,6> &, DistSVec<double,6> &,
            DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &,
            DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &);

  template<int dim, class Scalar>
  void computeDerivativeOfGradientsLeastSquares(dRdXoperators<dim> &,
            DistSVec<double,3> &, DistSVec<double,6> &, DistSVec<double,dim> &, 
            DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &);

  template<int dim, class Scalar>
  void computeTransposeDerivativeOfGradientsLeastSquares(dRdXoperators<dim> &,
            DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &,
            DistSVec<double,3> &,
            DistSVec<double,6> &,
            DistSVec<double,dim> &);

  template<int dim, class Scalar>
  void computeDerivativeOperatorsOfGradientsLeastSquares(DistSVec<double,3> &, DistSVec<double,6> &, DistSVec<Scalar,dim> &, dRdXoperators<dim> &);

  template<int dim, class Scalar>
  void computeDerivativeOfGradientsGalerkin(DistVec<double> &, DistVec<double> &,
                DistSVec<double,3> &, DistSVec<double,3> &,
				DistSVec<double,3> &, DistSVec<double,3> &,
				DistSVec<double,3> &, DistSVec<double,3> &,
				DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &,
				DistSVec<Scalar,dim> &, DistSVec<Scalar,dim> &);

  template<int dim>
  void computeDerivativeOfMultiDimLimiter(RecFcnLtdMultiDim<dim> *, DistSVec<double,3> &, DistSVec<double,3> &,
			      DistVec<double> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
			      DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
			      DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
			      DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &,
                  DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<int dim>
    void computeDerivativeOfFiniteVolumeTerm(DistVec<double> &, DistVec<double> &, DistVec<double> &, 
					     DistVec<double> &, FluxFcn**, RecFcn*, DistBcData<dim>&, 
					     DistGeoState&, DistSVec<double,3>&, DistSVec<double,3>&, 
					     DistSVec<double,dim>&, DistSVec<double,dim>&,
					     DistNodalGrad<dim>&, DistEdgeGrad<dim>*, double,
					     DistSVec<double,dim>&);

  template<int dim>
    void computeDerivativeOfFiniteVolumeTerm(dRdXoperators<dim> &, DistBcData<dim>&, DistGeoState&, DistSVec<double,3>&, 
					                                   DistNodalGrad<dim>&, DistEdgeGrad<dim>*, 
                                             DistSVec<double,dim>& dddx,
                                             DistSVec<double,dim>& dddy,
                                             DistSVec<double,dim>& dddz,
                                             DistVec<Vec3D>& dNormal,
                                             DistVec<Vec3D>& dn,
                                             DistVec<double>& dndot, DistSVec<double,dim>&);

  template<int dim>
    void computeTransposeDerivativeOfFiniteVolumeTerm(dRdXoperators<dim> &, 
					     DistBcData<dim>&, DistGeoState&, DistSVec<double,dim>&, 
					     DistNodalGrad<dim>&, DistEdgeGrad<dim>*, DistSVec<double,3>&,
               DistSVec<double,dim>&, DistSVec<double,dim>&, DistSVec<double,dim>&, DistVec<Vec3D>&,
               DistVec<Vec3D>&, DistVec<double>&);

  template<int dim>
    void computeDerivativeOperatorsOfFiniteVolumeTerm(DistVec<double> &, 
                                                      DistVec<double> &, FluxFcn**, RecFcn*, DistBcData<dim>&, 
                                                      DistGeoState&, DistSVec<double,3>&, DistSVec<double,dim>&, 
                                                      DistNodalGrad<dim>&, DistEdgeGrad<dim>*, double,
                                                      dRdXoperators<dim> &);

  template<int dim>
    void computeDerivativeOfFiniteVolumeTerm(FluxFcn** fluxFcn,  RecFcn* recFcn,
					     DistBcData<dim>& bcData, DistGeoState& geoState,
					     DistSVec<double,3> &X,
					     DistLevelSetStructure *distLSS,
					     bool linRecAtInterface, bool viscSecOrder, 
					     DistVec<int> &fluidId, 
					     DistExactRiemannSolver<dim> &riemann, 
					     int Nriemann,
					     DistNodalGrad<dim>& ngrad, 
					     DistEdgeGrad<dim>* egrad,
					     double dMach,
					     DistSVec<double,dim>& V,
					     DistSVec<double,dim>& dF);

  template<int dim>
  void computeDerivativeOfGalerkinTerm(FemEquationTerm *, DistBcData<dim> &,
			   DistGeoState &, DistSVec<double,3> &, DistSVec<double,3> &,
			   DistSVec<double,dim> &, DistSVec<double,dim> &, double, DistSVec<double,dim> &);

  template<int dim>
  void computeDerivativeOfGalerkinTerm(dRdXoperators<dim> &, FemEquationTerm *, DistBcData<dim> &,
			   DistGeoState &, DistSVec<double,3> &, DistSVec<double,3> &,
			   DistSVec<double,dim> &, DistSVec<double,dim> &, double, DistSVec<double,dim> &);

  template<int dim>
  void computeTransposeDerivativeOfGalerkinTerm(dRdXoperators<dim> &, DistSVec<double,dim> &, DistSVec<double,3> &);

  template<int dim>
  void computeDerivativeOperatorsOfGalerkinTerm(FemEquationTerm *, DistBcData<dim> &,
			   DistGeoState &, DistSVec<double,3> &,
			   DistSVec<double,dim> &, RectangularSparseMat<double,3,dim> **);

  template<int dim>
  void applyBCsToDerivativeOfResidual(BcFcn *, DistBcData<dim> &,
			  DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<int dim>
  void computeDerivativeOfSmagorinskyLESTerm(SmagorinskyLESTerm *, DistSVec<double,3> &,
				 DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<int dim>
  void computeOnlyGalerkinTerm(FemEquationTerm *, DistBcData<dim> &,
			   DistGeoState &, DistSVec<double,3> &,
			   DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<int dim>
  void computeBCsJacobianWallValues(FemEquationTerm *, DistBcData<dim> &,
				   DistGeoState &, DistSVec<double,3> &,
				   DistSVec<double,dim> &);

  template<int dim, class Scalar, int neq>
  void applyBCsToJacobianWallValues(BcFcn *, DistBcData<dim> &,
			  DistSVec<double,dim> &, DistMat<Scalar,neq> &);

  template<int dim, class Scalar2>
  void applyBCsToProduct(BcFcn *, DistBcData<dim> &, DistSVec<double,dim> &, DistSVec<Scalar2,dim> &);

  template<int dim>
  void computeDerivativeOfVolumicForceTerm(VolumicForceTerm *, DistVec<double> &, DistVec<double> &,
                               DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &);

  template<int dim>
  void computeDerivativeOfInvReynolds(FemEquationTerm *, VarFcn *, DistGeoState &,
			     DistSVec<double,3> &, DistSVec<double,3> &,
                             DistVec<double> &, DistVec<double> &, DistSVec<double,dim> &,
                             DistSVec<double,dim> &, DistVec<double> &, DistVec<double> &,
                             DistVec<double> &, DistVec<double> &, DistVec<double> &,
                             double, TimeLowMachPrec &, SpatialLowMachPrec &);

  template<int dim>
  void fixSolution(VarFcn *, DistSVec<double,dim> &, DistSVec<double,dim> &,DistVec<int>* fluidId = NULL);
  template<int dim>
  void fixSolution2(VarFcn *, DistSVec<double,dim> &, DistSVec<double,dim> &,DistVec<int>* fluidId = NULL);

  template<int dim>
  void getGradP(DistNodalGrad<dim>&);

  template<int dim>
  void getDerivativeOfGradP(DistNodalGrad<dim>&);

  template<int dim>
  void computeDerivativeOperatorOfGradP(RectangularSparseMat<double,dim,3> **dGradPdddx,
                                        RectangularSparseMat<double,dim,3> **dGradPdddy,
                                        RectangularSparseMat<double,dim,3> **dGradPdddz);

  template<int dim>
  void getDerivativeOfGradP(RectangularSparseMat<double,dim,3> **dGradPdddx,
                            RectangularSparseMat<double,dim,3> **dGradPdddy,
                            RectangularSparseMat<double,dim,3> **dGradPdddz,
                            DistSVec<double,dim> &dddx,
                            DistSVec<double,dim> &dddy,
                            DistSVec<double,dim> &dddz,
                            DistSVec<double,3> &dGradP);

  template<int dim>
  void getTransposeDerivativeOfGradP(RectangularSparseMat<double,dim,3> **dGradPdddx,
                                     RectangularSparseMat<double,dim,3> **dGradPdddy,
                                     RectangularSparseMat<double,dim,3> **dGradPdddz,
                                     DistSVec<double,3> &dGradP,
                                     DistSVec<double,dim> &dddx,
                                     DistSVec<double,dim> &dddy,
                                     DistSVec<double,dim> &dddz);

  void updateNodeTag(DistSVec<double,3> &, DistLevelSetStructure *, DistVec<int> &, DistVec<int> &);

  void computeMaterialVolumes(double*, int, DistVec<double> &, DistVec<int> *);

  template<int dim>
  void computeMaterialConservationScalars(double*, double*, double*, double*, double*, int, DistSVec<double, dim> &, DistVec<double> &, DistVec<int> *);

  void computeCharacteristicEdgeLength(DistSVec<double,3> &, double&, double&, double&, int&, const double, const double, const double, const double, const double, const double);

  template<int dim>
  void computeCVBasedForceLoad(int, int, DistGeoState&, DistSVec<double,3>&, double (*)[3], int, DistLevelSetStructure*, 
                               double, DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &V, 
                               DistVec<GhostPoint<dim>*> *ghostPoints, PostFcn *postFcn,DistNodalGrad<dim, double> *ngrad,
                               VarFcn* vf, DistVec<int>* fid);

  template<int dim>
  void computeEmbSurfBasedForceLoad(IoData &iod, int, int, DistSVec<double,3>&, double (*)[3], int, 
												DistLevelSetStructure*, double, 
												DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji, 
												DistSVec<double,dim> *Wextij, DistSVec<double,dim> &V, 
                                    DistVec<GhostPoint<dim>*> *ghostPoints, PostFcn *postFcn, 
												DistNodalGrad<dim, double> *ngrad, VarFcn* vf, DistVec<int>* fid, bool externalSI=false);

  template<int dim>
  void computederivativeEmbSurfBasedForceLoad(IoData &iod, int forceApp, int orderOfAccuracy, DistSVec<double,3> &X, 
					      double (*dFs)[3], int sizedFs, DistLevelSetStructure *distLSS, 
					      double pInfty, double dpInfty, 
					      DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji, 
					      DistSVec<double,dim> &V, DistSVec<double,dim> &dV_,
					      DistVec<GhostPoint<dim>*> *ghostPoints, PostFcn *postFcn, 
					      DistNodalGrad<dim, double> *gradV, DistNodalGrad<dim, double> *graddV, 
					      VarFcn* vf, DistVec<int> *fid);
  template<int dim>
  void computeRecSurfBasedForceLoad(int, int, DistSVec<double,3>&, double (*)[3], int, DistLevelSetStructure*, double, 
                                    DistSVec<double,dim> &Wstarij, DistSVec<double,dim> &Wstarji, DistSVec<double,dim> &V, 
                                    DistVec<GhostPoint<dim>*> *ghostPoints, PostFcn *postFcn, VarFcn* vf, DistVec<int>* fid);
  template<int dim>
  void computePrdtWCtrlVolRatio(DistSVec<double,dim> &, DistSVec<double,dim> &, DistVec<double> &, DistGeoState &);

  template<int dimLS>
  void computePrdtPhiCtrlVolRatio(DistSVec<double,dimLS> &, DistSVec<double,dimLS> &, DistVec<double> &, DistGeoState &);

  template<int dim>
  void blur(DistSVec<double,dim> &U, DistSVec<double,dim> &U0);

  template<int dimLS>
  void updateFluidIdFS2Prep(DistLevelSetStructure &distLSS, DistSVec<double,dimLS> &PhiV, DistVec<int> &fluidId, DistSVec<bool,4> &poll);

  template<int dim, int dimLS>
  void debugMultiPhysics(DistLevelSetStructure &distLSS, DistSVec<double,dimLS> &PhiV, DistVec<int> &fluidId, DistSVec<double,dim> &U);

  template<int dim, class Obj>
  void integrateFunction(Obj* obj,DistSVec<double,3> &X,DistSVec<double,dim>& V, void (Obj::*F)(int node, const double* loc,double* f),
                         int npt);

  void createHigherOrderMultiFluid();

  void createHigherOrderFSI();

  // When a cell is omitted when doing higher order multi-fluid calculations, we can grab
  // a value of the state for the cut cell using an extrapolated state

  // Functions to compute the error (that is, the difference between two state vectors)
  template <int dim>
    void computeL1Error(DistSVec<double,dim>& U, DistSVec<double,dim>& Uexact, 
			DistVec<double>& vol,double error[dim],
                        DistLevelSetStructure* = NULL);


  template <int dim>
    void computeL2Error(DistSVec<double,dim>& U, DistSVec<double,dim>& Uexact, 
			DistVec<double>& vol,double error[dim],
                        DistLevelSetStructure* = NULL);

  template <int dim>
    void computeLInfError(DistSVec<double,dim>& U, DistSVec<double,dim>& Uexact, double error[dim],
                          DistLevelSetStructure* = NULL);


  // Assign ErrorHandler to subdomains
  void assignErrorHandler();

  void setFarFieldMask(DistVec<double> &ffMask, DistVec<double> &neighborMask); // ffMask is one for farfield nodes,
                                                                                // neighborMask is one for nodes that are one layer removed from the farfield
  void setWallMask(DistVec<double> &wallMask, DistVec<double> &neighborMask); // wallMask is one for wall nodes,
                                                                              // neighbor mask is one for nodes that are one layer removed from the wall

  template <int dim>
    void computeHHBoundaryTermResidual(DistBcData<dim> &bcData,DistSVec<double,dim> &U,DistVec<double>& res,
				       VarFcn* vf);
 
  void maskHHVector(DistVec<double>& hh);


  template <int dim>
    void setExactBoundaryValues(DistSVec<double,dim>& U, DistSVec<double,3>& X,
				IoData& iod,double t,VarFcn* vf);

  template <int dim>
    void setExactBoundaryResidual(DistSVec<double,dim>& F, DistSVec<double,3>& X,
				  IoData& iod,double t,VarFcn* vf);

  
  template <int dim,int neq, class Scalar>
    void setExactBoundaryJacobian(DistSVec<double,dim>& U, DistSVec<double,3>& X,
				  IoData& iod,double t, VarFcn* varFcn,
				  DistMat<Scalar,neq>& A);
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <Domain.C>
#endif

#endif
