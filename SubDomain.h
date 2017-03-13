#ifndef _SUBDOMAIN_H_
#define _SUBDOMAIN_H_

#include <IoData.h>
#include <PostFcn.h>
#include <Node.h>
#include <InletNode.h>
#include <Edge.h>
#include <Face.h>
#include <Elem.h>
#include <DiagMatrix.h>
#include <DistInfo.h>
#include <BCond.h>
#include <TriangulatedSurface.h>
#include <DenseMatrix.h>
#include <GhostPoint.h>
#include <PolygonReconstructionData.h>

#include <HigherOrderMultiFluid.h>
#include <HigherOrderFSI.h>

#include <ErrorHandler.h>

#include <Aerof_unordered_set.h>

#ifdef OLD_STL
#include <map.h>
#else
#include <map>
using std::map;
#endif

#include <set>

#include <complex>
typedef std::complex<double> bcomp;

class VarFcn;
class PostFcn;
class BcFcn;
class RecFcn;
class FluxFcn;
class FemEquationTerm;
class MacroCellSet;
class VMSLESTerm;
class DynamicVMSTerm;
class DynamicLESTerm;
class SmagorinskyLESTerm;
class WaleLESTerm;
class SubDTopo;
class TimeData;
class GeoSource;
class GeoState;
class Connectivity;
class Communicator;
class SubMat;
class SubDiagPrec;
class MemoryPool;
class BinFileHandler;
class BCApplier;
class BCData;
class MatchNodeSet; // HB
class VolumicForceTerm;
class TriangulatedSurface;
class TimeLowMachPrec;
class FluidSelector;
class LevelSetStructure;

struct V6NodeData;
struct Vec3D;
struct compStruct;
struct ExtrapolationNodeData;

template<int dimLS> class LevelSet;
template<int dim, class Scalar> class NodalGrad;
template<int dim> class EdgeGrad;
template<int dim> class Extrapolation;
template<int dim> class BcData;
template<int dim> class ExactRiemannSolver;
template<class Scalar> class Vec;
template<class Scalar> class CommPattern;
template<class Scalar, int dim> class SVec;
template<class Scalar, int dim> class MvpMat;
template<class Scalar, int dim> class SparseMat;
template<class Scalar, int dim> class GenMat;

template <class T> class RTree;

//------------------------------------------------------------------------------

struct EdgeDef {

  EdgeDef() { }

  EdgeDef(int glLeft, int glRight, int edgeNum, 
          int sign) : glLeft(glLeft), glRight(glRight),
                      edgeNum(edgeNum), sign(sign) {

  }

  EdgeDef(const EdgeDef& oth) : glLeft(oth.glLeft), glRight(oth.glRight),
                                edgeNum(oth.edgeNum), sign(oth.sign) { }

  int glLeft, glRight, edgeNum, sign;

  bool operator<(const EdgeDef &e) const
  {
    return ( (glLeft < e.glLeft) || (glLeft == e.glLeft && glRight < e.glRight) );
  }

  void order()
  {
    if (glLeft < glRight) sign = 1;
    else { int tmp = glRight; glRight = glLeft; glLeft = tmp; sign = -1; }
  }

};

//------------------------------------------------------------------------------

/** \brief SubDomain data
 *
 */
class SubDomain {
  int testEdge;

  int locSubNum;
  int clusSubNum;
  int globSubNum;
  int numClusNodes;
  int numNodeRanges;

  char suffix[100];

  NodeSet &nodes;
  EdgeSet  edges;
  FaceSet &faces;
  ElemSet &elems;

  InletNodeSet inletNodes;

  // for weighting the residual in Nonlinear ROM simulations
  std::vector<int> farFieldNodes;
  std::vector<int> wallNodes;

  int *locToGlobNodeMap;
  int *locToGlobFaceMap;
  int *locToGlobElemMap;

  int numNeighb;
  int *neighb;
  int *sndChannel;
  int *rcvChannel;
  Connectivity *sharedNodes;
  Connectivity *sharedInletNodes;
  Connectivity** nodesToMCNodes;
  Connectivity *NodeToSubD;
  int *numSharedEdges;
  EdgeDef **sharedEdges;
  int (*nodeRanges)[3];
  int *nodeType;
  int *nodeFaceType;
  int *offWallNode;

  map<int, int> bcMap;
  int *numBcNodes;

  BCondSet *mmsBCs;
  int *rotOwn;

  int **nodeFlag;
  Connectivity *NodeToNode;
  Connectivity *NodeToElem;
  Connectivity *ElemToElem;
  // Adam 2011.01.19: For Memory Leak Purpose
  Connectivity *nodeToNodeMaskJacobian, *nodeToNodeMaskILU;

  int **totalNeiData;
  double *gradP[3];

// Included (MB*)
  int numOffDiagEntries;
  double *dGradP[3];
  SVec<double,3> *dGradPSVec;

  bool sampleMesh;
  int numSampledNodes;
  std::vector<int> locSampleNodes;	// for Gappy ROM
  

  // UH (07/2012)
  // List of nodes on a surface for the Kirchhoff integral
  std::set<int> kirchhoffNodesList;

  RTree<Elem>* myTree;

  ErrorHandler* errorHandler;

  LevelSetStructure* triangulatedInterfaceLSS;

  double geomTol;
  double bigNum;
  
public:
  
  HigherOrderMultiFluid* higherOrderMF;
  HigherOrderFSI* higherOrderFSI;

  SubDomain(int, int, int, int, char *, NodeSet *, FaceSet *, ElemSet *,
	    int, int *, Connectivity *, int *, int *, int *, int, int (*)[3]);
  ~SubDomain();

  // topology
  int *getNodeMap()    { return locToGlobNodeMap; }
  int *getElemMap()  { return locToGlobElemMap; }
  int getGlobSubNum()  { return globSubNum; }
  int getLocSubNum()   { return locSubNum; }
  int getNumNeighb()   { return numNeighb; }
  int *getNeighb()     { return neighb; }
  int *getSndChannel() { return sndChannel; }
  int *getRcvChannel() { return rcvChannel; }
  Connectivity* getSharedNodes() {return sharedNodes;}
  int numberEdges();

  EdgeDef** getSharedEdges() { return sharedEdges; }
  const int* getNumSharedEdges() const { return numSharedEdges; }

  void computeConnectedTopology(const std::vector<int> &locSampleNodes, const std::vector<int> &globalNeighborNodes_);

  Connectivity *createElemBasedConnectivity();
  Connectivity *createNodeToElementConnectivity();
  Connectivity *createElementToElementConnectivity();
  Connectivity *createElementToNodeConnectivity();
  Connectivity *createEdgeBasedConnectivity();
  Connectivity *createNodeToConstantConnectivity();
  Connectivity *createConstantToNodeConnectivity();
  Connectivity *createConstantToConstantConnectivity();
  Connectivity *createElementBasedEdgeToNodeConnectivity();
  Connectivity *createElementBasedNodeToEdgeConnectivity();
  Connectivity *createFaceToNodeConnectivity();
  Connectivity *createNodeToFaceConnectivity();
  Connectivity *createNodeToSubDomainConnectivity();
  Connectivity *createNodeToMacroCellNodeConnectivity(MacroCellSet *);
  Connectivity *agglomerate(Connectivity &, int, bool *);
  void createSharedInletNodeConnectivity(int);

  compStruct *createRenumbering(Connectivity *, int, int);

  int numNodes() { return(nodes.size()); }
  int numFaces() { return(faces.size()); }
  int numElems() { return(elems.size()); }
  int numEdges() { return(edges.size()); }
	FaceSet& getFaces() {return faces;};
	ElemSet& getElems() {return elems;};

  int* getElemNodeNum(int i) {return(elems[i].nodeNum()); }

  // Get the local node number in the subdomain of the global node <id>
  // Returns -1 if it does not exist.  Warning: This method is O(N)
  int getLocalNodeNum(int globNodeNum) const;
  // geometry

  void getSurfaceNodes(Aerof_unordered_set<int>::type& boundaryNodes) const;
  void getSolidBoundaryNodes(Aerof_unordered_set<int>::type& boundaryNodes) const;
  void getFarFieldBoundaryNodes(Aerof_unordered_set<int>::type& boundaryNodes) const;
  void getSubDomainBoundaryNodes(Aerof_unordered_set<int>::type& boundaryNodes) const;

  void constructLines(std::vector<std::vector<int>*>& pLines, int& numLines);

  void setFaceType(int *);
  void setNodeType(int*, CommPattern<int> &);
  void setNodeFaceType(CommPattern<int> &);
  int* completeNodeType(int*, CommPattern<int> &);
  int* completeNodeFaceType(CommPattern<int> &);
  void computeOffWallNode(LevelSetStructure *, CommPattern<int> &);
  int* completeOffWallNode(CommPattern<int> &);
  int setFaceToElementConnectivity();
  void getElementStatistics(int &, int &, int &, int &);
  int computeControlVolumes(int, double, SVec<double,3> &, Vec<double> &);
  void computeFaceNormals(SVec<double,3> &, Vec<Vec3D> &);
  void computeFaceEdgeNormals(SVec<double,3>&, SVec<double,6>&);
  void computeEdgeDihedralAngle(double, SVec<double,6>&, Vec<double>&);
  void propagateInfoAlongEdges(Vec<double>&);
  void computeNormalsGCL1(SVec<double,3> &, SVec<double,3> &, SVec<double,3> &,
			  Vec<Vec3D> &, Vec<double> &, Vec<Vec3D> &, Vec<double> &);
  void computeNormalsConfig(SVec<double,3> &Xconfig, SVec<double,3> &Xdot,
                            Vec<Vec3D> &edgeNorm, Vec<double> &edgeNormVel,
                            Vec<Vec3D> &faceNorm, Vec<double> &faceNormVel);
  void computeNormalsEZGCL1(double, SVec<double,3>&, SVec<double,3>&, Vec<Vec3D>&,
			    Vec<double>&, Vec<Vec3D>&, Vec<double>&);
  void computeSmoothedSensor(SVec<double,3>&, Vec<double>&, SVec<double,3>&);
  void computeWeightsLeastSquaresEdgePart(SVec<double,3> &, SVec<double,6> &);
  void computeWeightsLeastSquaresNodePart(SVec<double,6> &);
  void computeWeightsLeastSquaresEdgePart(SVec<double,3> &, const Vec<int> &,
					  SVec<int,1> &, SVec<double,6> &, LevelSetStructure* =0,
					  bool includeSweptNodes = true);
  void computeWeightsLeastSquaresEdgePart(SVec<double,3> &, const Vec<int> &,
					  SVec<int,1> &, SVec<double,6> &, Vec<int> &, Vec<int> &, 
					  LevelSetStructure* =0);
  void computeWeightsLeastSquaresNodePart(SVec<int,1> &, SVec<double,6> &);

  void computeWeightsGalerkin(SVec<double,3> &, 
			      SVec<double,3> &,
			      SVec<double,3> &, 
			      SVec<double,3> &);

  void computeWeightsGalerkin(SVec<double,3> &, const Vec<int> &,
			      SVec<double,3> &,
			      SVec<double,3> &, 
			      SVec<double,3> &,
			      LevelSetStructure* =0, 
			      bool includeSweptNodes = true); //d2d

  void computeEdgeWeightsGalerkin(SVec<double,3> &, Vec<double> &, SVec<double,9> &);
#define EDGE_LENGTH
#ifdef EDGE_LENGTH
  bool findTetrahedron(int, int, Vec<int>&, int**, SVec<double,3>&, V6NodeData&, bool=true, double* refLength=0);
#else
  bool findTetrahedron(int, int, Vec<int>&, int**, SVec<double,3>&, V6NodeData&, bool=true);
#endif
  bool findNormalTet1(Vec3D, Vec3D, Vec3D, Vec3D, Vec3D, int, int, ExtrapolationNodeData &, bool);
  void findNormalTet2(int*, Vec3D, Vec3D, int, SVec<double,3> , ExtrapolationNodeData &);
  bool findNormalTetrahedron(int , Vec3D , int *, int *, int ,
                             SVec<double,3> , ExtrapolationNodeData *, bool =true);
  void findEdgeTetrahedra(SVec<double,3>&, V6NodeData (*&)[2]);
  void findNormalTetrahedra(SVec<double,3>& , Vec<Vec3D>& , ExtrapolationNodeData (*&)[2]);
  void checkNormalTetrahedra(ExtrapolationNodeData (*&)[2]);
  void setInletNodes(IoData& );
  void setInletNodes2(IoData& );
  void setFarFieldNodes();
  void setWallNodes();
  int findOppositeTet(int* , int );
  bool findNodeInTet(int, int);
  void checkInletNodes();
  void sumInletNormals(Vec<Vec3D>&,     Vec<Vec3D>&, Vec<int>&);
  void numDivideNormals(Vec<Vec3D>&, Vec<int>&);
  void getReferenceMeshPosition(SVec<double,3> &);
  void computeDisplacement(SVec<double,3> &, SVec<double,3> &);
  void computeDisplacement(SVec<double,3> &, double* res,int id);

  //d2d
  bool getPiercedPoint(Vec3D va, Vec3D vb, Vec3D vc, Vec3D vd, Vec3D ve, double &r, double &t);
  bool getSIstencil(int Ni, int Nj, SVec<double,3> &X, LevelSetStructure &LSS, Vec<int> &fluidId,
						  Vec3D &normWall, Vec3D &xWall, V6NodeData &SiStencilData);
  bool getFEMstencil(int Ni, SVec<double,3> &X, LevelSetStructure &LSS, Vec<int> &fluidId,
							Vec3D &normWall, Vec3D &xWall, V6NodeData &NodeStencilData);

  bool getFEMstencil2(int Ni, SVec<double,3> &X, LevelSetStructure &LSS, 
							 Vec<int> &fluidId, Vec3D &normWall, Vec3D &xWall, 
							 V6NodeData &NodeStencilData_p, V6NodeData &NodeStencilData_m);

  void printPoint(int Ni, SVec<double,3> &X, Vec3D &xWall, V6NodeData &NSdata, int dir, int mode);

  void assimilateCells(int, int, int*, int**, bool *, int, int *, bool *);
  MacroCellSet** findAgglomerateMesh(int, int, bool *, double);
  void createMacroCellConnectivities(MacroCellSet **, int);

  void setBCond(BCondSet *subBC) { mmsBCs = subBC; }
  void applySmoothing(Vec<double> &, Vec<double> &);
  void applySmoothing(Vec<double> &, SVec<double,2> &);
  void computeTetsConnectedToNode(Vec<int> &Ni);
  void computeLocalAvg(SVec<double,3> &, Vec<double> &, Vec<double> &);
  void computeLocalAvg(SVec<double,3> &, SVec<double,2> &, SVec<double,2> &);
  void computeFilterWidth(SVec<double,3> &, Vec<double> &);
  void finalizeTags(SVec<int,2> &);
  template<int dimLS>
  void setupPhiVolumesInitialConditions(const int volid, const int fluidId, SVec<double,dimLS> &Phi);
  void outputCsDynamicLES(DynamicLESTerm *, SVec<double,2> &, SVec<double,3> &, Vec<double> &);

  template<int dimLS>
  void avoidNewPhaseCreation(SVec<double,dimLS> &Phi, SVec<double,dimLS> &Phin);
  template<int dimLS>
  void avoidNewPhaseCreation(SVec<double,dimLS> &Phi, SVec<double,dimLS> &Phin, Vec<double> &weight, LevelSetStructure *LSS = 0, 
          Vec<int>* fluidIdToSet = 0);
  template<int dim>
  void setupUVolumesInitialConditions(const int volid, double UU[dim], SVec<double,dim> &U);

  void setupFluidIdVolumesInitialConditions(const int volid, const int myId, Vec<int> &fluidId);
  //template<int dim>
  //void setupUMultiFluidInitialConditionsSphere(FluidModelData &fm,
  //           SphereData &ic, SVec<double,3> &X, SVec<double,dim> &U);
  //template<int dim>
  //void setupUMultiFluidInitialConditionsPlane(FluidModelData &fm,
  //           PlaneData &ip, SVec<double,3> &X, SVec<double,dim> &U);
  //template<int dim>
  //void setupUMultiFluidInitialConditionsPlane(FluidModelData &fm,
  //           PlaneData &ip, SVec<double,3> &X, SVec<double,dim> &U, Vec<int> &nodeTag);

  void computeLij(double [3][3], double [3], double [6], double [5]);
  void computeBij(double [3][3], double [6], double, double [3][3], double, double [5]);
  void computeZi(double [3], double, double, double [3], double [3], double [5], double);
  void computeLi(double [3], double, double, double [3], double [5]);

  void findNodeBoundingBoxes(SVec<double,3>&X, SVec<double,3> &Xmin, SVec<double,3> &Xmax);
  void findSubDomainBoundingBoxes(SVec<double,3>&X, double *Xmin, double *Xmax);

  // moving mesh

  void getNdAeroLists(int &, int *&, int &, int *&, int &, int *&, MatchNodeSet* matchNodes=0);

  template<class MatScalar, class PrecScalar>
  void computeStiffAndForce(DefoMeshMotionData::Element, SVec<double,3>&, SVec<double,3>&, GenMat<MatScalar,3>&,
                            GenMat<PrecScalar,3>*, double volStiff, int* ndType = 0);

  // spatial discretization
  template<int dim>
  void computeTimeStep(FemEquationTerm *, VarFcn *, GeoState &, SVec<double,3> &, SVec<double,dim> &, Vec<double> &,
		       Vec<double> &, Vec<double> &, Vec<double> &,
                       TimeLowMachPrec &, LevelSetStructure *LSS=0);

  template<int dim>
  void computeTimeStep(FemEquationTerm *, VarFcn *, GeoState &, SVec<double,dim> &, Vec<double> &,
		       Vec<double> &, Vec<double> &, Vec<double> &,
                       TimeLowMachPrec &, Vec<int> &, Vec<double>* = NULL);


  template<int dim, class Scalar>
  void computeGradientsLeastSquares(SVec<double,3> &, SVec<double,6> &,
				    SVec<Scalar,dim> &, SVec<Scalar,dim> &,
				    SVec<Scalar,dim> &, SVec<Scalar,dim> &);

  template<int dim, class Scalar>
  void computeGradientsLeastSquares(SVec<double,3> &, const Vec<int> &,
                                    SVec<double,6> &,
                                    SVec<Scalar,dim> &, SVec<Scalar,dim> &,
                                    SVec<Scalar,dim> &, SVec<Scalar,dim> &,
                                    bool linRecFSI = true, LevelSetStructure* =0,
                                    bool includeSweptNodes = true);

  template<class Scalar>
  void computeGradientLeastSquares(SVec<double,3> &, const Vec<int> &,
                                    SVec<double,6> &,
                                    Vec<Scalar> &, Vec<Scalar> &,
                                    Vec<Scalar> &, Vec<Scalar> &,
                                    LevelSetStructure* =0);

  template<int dim, class Scalar>
  void computeGradientsLeastSquares(SVec<double,3> &, const Vec<int> &,
                                    SVec<double,6> &,
                                    SVec<Scalar,dim> &, SVec<Scalar,dim> &,
									SVec<Scalar,dim> &, Vec<int> &,
									Vec<int> &, SVec<Scalar,dim> &,
                                    SVec<Scalar,dim> &, SVec<Scalar,dim> &,
                                    bool linRecFSI = true, LevelSetStructure* =0);

  template<int dim, class Scalar>
  void computeGradientsGalerkin(Vec<double> &, SVec<double,3> &, SVec<double,3> &,
				SVec<double,3> &, SVec<Scalar,dim> &, SVec<Scalar,dim> &,
				SVec<Scalar,dim> &, SVec<Scalar,dim> &);

  template<class Scalar1, class Scalar2, int dim>
  void computeGradientsTransposeNew(SVec<double,3> &X,
                                    Vec<double> &ctrlVol, GenMat<Scalar1,dim> &A,
                                    SVec<double,dim> &bij, SVec<double,dim> &bji,
                                    SVec<Scalar2,dim> &cij, SVec<Scalar2,dim> &cji,
                                    SVec<Scalar2,dim> &dij, SVec<Scalar2,dim> &dji,
                                    SVec<Scalar2,dim> &p, SVec<Scalar2, dim> &ddxt,
                                    SVec<Scalar2, dim> &ddyt, SVec<Scalar2, dim> &ddzt);

  template<int dim, class Scalar>
  void computeGradientsGalerkinT(Vec<double> &, SVec<double,3> &,
                SVec<double,3> &, SVec<double,3> &, SVec<Scalar,dim> &,
                SVec<Scalar,dim> &, SVec<Scalar,dim> &, SVec<Scalar,dim> &,
                SVec<Scalar,dim> &, SVec<Scalar,dim> &);

  template<int dim>
	  void computeMinMaxStencilValues(SVec<double,dim> &, SVec<double,dim> &, 
												 SVec<double,dim> &, LevelSetStructure *LSS = 0);

  template<int dim>
//  void computeMultiDimLimiter(RecFcnLtdMultiDim<dim> *, SVec<double,3> &, Vec<double> &,
  void computeMultiDimLimiter(RecLimiter *, SVec<double,3> &, Vec<double> &,
			      SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &,
			      SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &,
										SVec<double,dim> &, LevelSetStructure *LSS = 0);

  template<int dim>
  void computePressureSensor(SVec<double,3>&, SVec<double,dim>&, SVec<double,dim>&,
			     SVec<double,dim>&, SVec<double,dim>&, SVec<double,3>&);

  template<int dim>
  int computeFiniteVolumeTerm(Vec<double> &, FluxFcn**, RecFcn*, BcData<dim>&, GeoState&,
                              SVec<double,3>&, SVec<double,dim>&, NodalGrad<dim>&,
                              EdgeGrad<dim>*, SVec<double,dim>&, SVec<int,2>&, int, int);

  template<int dim>
  int computeFiniteVolumeTerm(ExactRiemannSolver<dim>&,
                              Vec<double> &, FluxFcn**, RecFcn*, BcData<dim>&, GeoState&,
                              SVec<double,3>&, SVec<double,dim>&, NodalGrad<dim>&,
                              EdgeGrad<dim>*, SVec<double,dim>&, SVec<int,2>&, int, int);

  template<int dim, int dimLS>
  int computeFiniteVolumeTerm(ExactRiemannSolver<dim>&,
                              FluxFcn**, RecFcn*, 
										BcData<dim>&, GeoState&,
                              SVec<double,3>&, SVec<double,dim>&, 
                              Vec<int> &, FluidSelector &,
                              NodalGrad<dim>&, EdgeGrad<dim>*,
			      SVec<double,dimLS>& phi,
                              NodalGrad<dimLS>&,
			      EdgeGrad<dimLS>*,
                              SVec<double,dim>&, int, SVec<int,2>&, int, int);

  template<int dim, int dimLS>
  int computeFiniteVolumeTerm(ExactRiemannSolver<dim>&,
                              FluxFcn**, RecFcn*,BcData<dim>&, GeoState&,
                              SVec<double,3>&, SVec<double,dim>&,
                              SVec<double,dim>&, SVec<double,dim>&, LevelSetStructure&, bool,
                              Vec<int> &, int, FluidSelector &,
                              NodalGrad<dim>&, EdgeGrad<dim>*,
			      SVec<double,dimLS>& phi,
                              NodalGrad<dimLS>&,
			      EdgeGrad<dimLS>*,
                              SVec<double,dim>&, int, SVec<int,2>&, int, int);

  // embedded structure
  template<int dim>
  int computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                              FluxFcn** fluxFcn, RecFcn* recFcn, 
										BcData<dim>& bcData, GeoState& geoState,
                              SVec<double,3>& X, SVec<double,dim>& V,
                              SVec<double,dim>& Wstarij, SVec<double,dim>& Wstarji, 
										SVec<double,dim>& Wext, LevelSetStructure &LSS, 
										bool linRecAtInterface, Vec<int> &fluidId, 
										int Nriemann,
                              NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
                              SVec<double,dim>& fluxes, 
										int it, SVec<int,2>& tag, 
										int failsafe, int rshift, bool externalSI=false);

  template<int dim>
  int computeFiniteVolumeTerm(ExactRiemannSolver<dim>&,
                              FluxFcn**, RecFcn*, 
										BcData<dim>&, GeoState&,
                              SVec<double,3>&, SVec<double,dim>&,
                              SVec<double,dim>&, SVec<double,dim>&, 
                              Vec<int>&, Vec<int>&, LevelSetStructure &, bool, 
                              Vec<int> &, int, double, double, 
                              NodalGrad<dim>&, EdgeGrad<dim>*,
                              SVec<double,dim>&, int, SVec<int,2>&, int, int); 

  template<int dim, int dimLS>
  void computeFiniteVolumeTermLS(FluxFcn**, RecFcn*, RecFcn*, BcData<dim>&, GeoState&,
                               SVec<double,3>&, SVec<double,dim>&,Vec<int>& fluidId,
                               NodalGrad<dim>&,   EdgeGrad<dim>*,
			       NodalGrad<dimLS>&, EdgeGrad<dimLS>*,
			       SVec<double,dimLS>&,
                               SVec<double,dimLS>&, LevelSetStructure* =0, int ls_order = 1);
  template<int dim>
  int computeFiniteVolumeBar_Step1(Vec<double> &, FluxFcn**, RecFcn*, BcData<dim>&, GeoState&, SVec<double,3>& ,
                                    SVec<double,dim>&, NodalGrad<dim> &, EdgeGrad<dim>* , SVec<double,dim>&,
                                    SVec<int,2> &, int, int);

  template<int dim>
  void computeFiniteVolumeBar_Step2(MacroCellSet **,SVec<double,1> &, SVec<double,dim> &, SVec<double,dim> &, int);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, BcData<dim> &, GeoState &,
                                       Vec<double> &, SVec<double,3> &, Vec<double> &,
                                       SVec<double,dim> &, GenMat<Scalar,neq> &,
                                       CommPattern<double> *);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>&,
                                       FluxFcn **, BcData<dim> &, GeoState &,
                                       Vec<double> &, SVec<double,3> &, Vec<double> &,
                                       SVec<double,dim> &, GenMat<Scalar,neq> &,
                                       CommPattern<double> *);

  template<int dim, class Scalar, int neq, int dimLS>
  void computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>&,
                                       FluxFcn **, BcData<dim> &, GeoState &,
                                       NodalGrad<dim> &, NodalGrad<dimLS> &,
                                       SVec<double,3> &, Vec<double> &,
                                       SVec<double,dim> &, GenMat<Scalar,neq> &,
                                       FluidSelector &, 
                                       Vec<int> &, CommPattern<double> *);

  template<class Scalar,int dim,int neq>
  void computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                       FluxFcn** fluxFcn,
                                       BcData<dim>& bcData, GeoState& geoState,
                                       SVec<double,3>& X, SVec<double,dim>& V,Vec<double>& ctrlVol,
                                       LevelSetStructure &LSS,Vec<int> &fluidId,
                                       int Nriemann,
                                       GenMat<Scalar,neq>& A,Vec<double>& irey,  bool externalSI=false);

  template<int dim, class Scalar, int neq, int dimLS>
  void computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,
                                       FluxFcn** fluxFcn, 
                                       BcData<dim>& bcData, GeoState& geoState,
                                       SVec<double,3>& X, SVec<double,dim>& V,Vec<double>& ctrlVol,
                                       NodalGrad<dimLS> &ngradLS,
                                       LevelSetStructure &LSS,Vec<int> &fluidId,
                                       int Nriemann,
                                       FluidSelector &fluidSelector,
                                       GenMat<Scalar,neq>& A) ;

  template<int dim, class Scalar, int dimLS>
    void computeJacobianFiniteVolumeTermLS(RecFcn* recFcn, RecFcn* recFcnLS,
					   GeoState &geoState,SVec<double,3>& X,SVec<double,dim> &V,
					   NodalGrad<dim>& ngrad,
					   NodalGrad<dimLS> &ngradLS,
					   EdgeGrad<dim>* egrad,
					   Vec<double> &ctrlVol,SVec<double,dimLS>& Phi,
					   GenMat<Scalar,dimLS> &A, LevelSetStructure* LSS,CommPattern<double> * flag);

  template<int dim>
  void recomputeRHS(VarFcn*, SVec<double,dim>& ,SVec<double,dim>& , Extrapolation<dim>*,
                                        BcData<dim>&, GeoState&, SVec<double,3> &);
  template<int dim>
  void recomputeRHS(VarFcn*, SVec<double,dim>& ,Vec<int> &, SVec<double,dim>&,
                    Extrapolation<dim>*, BcData<dim>&, GeoState&, SVec<double,3> &);

  template<int dim>
  void recomputeResidual(SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeRealFluidResidual(SVec<double, dim> &, SVec<double,dim> &, LevelSetStructure &);


  template<class Scalar,int dim>
  void checkRHS(Scalar (*)[dim]);

  template<int dim>
  void computeGalerkinTerm(FemEquationTerm *, BcData<dim> &, GeoState &,
									SVec<double,3> &, SVec<double,dim> &, 
									SVec<double,dim> &,
									Vec<GhostPoint<dim>*> *ghostPoints=0,
									LevelSetStructure *LSS=0, 
									bool externalSI=false);

  template<int dim>
  void computeVolumicForceTerm(VolumicForceTerm *, Vec<double> &,
                               SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeSmagorinskyLESTerm(SmagorinskyLESTerm *, SVec<double,3> &, SVec<double,dim> &,
				 SVec<double,dim> &, Vec<GhostPoint<dim>*> *ghostPoints=0,
                                 LevelSetStructure *LSS=0, bool externalSI=false);

  template<int dim>
  void computeWaleLESTerm(WaleLESTerm *, SVec<double,3> &, SVec<double,dim> &, 
                          SVec<double,dim> &, Vec<GhostPoint<dim>*> *ghostPoints=0,
                          LevelSetStructure *LSS=0, bool externalSI=false);

  template<int dim>
  void computeMutOMuSmag(SmagorinskyLESTerm *, SVec<double,3> &, SVec<double,dim> &, Vec<double> &);

  template<int dim>
  void computeMutOMuVMS(VMSLESTerm *, SVec<double,dim> &, SVec<double,3> &, SVec<double,dim> &, Vec<double> &);

  template<int dim>
  void computeMutOMuDynamicVMS(DynamicVMSTerm *, SVec<double,dim> &, SVec<double,3> &,
                               SVec<double,dim> &, Vec<double> &, Vec<double> &);

  template<int dim>
  void computeMutOMuWale(WaleLESTerm *, SVec<double,3> &, SVec<double,dim> &, Vec<double> &);

  template<int dim>
  void computeMutOMuDynamicLES(DynamicLESTerm *, SVec<double,2> &, SVec<double,3> &,
                               SVec<double,dim> &, Vec<double> &);

  template<int dim>
  void computeTestFilterAvgs(SVec<double,dim> &,  SVec<double,16> &, SVec<double,6> &, Vec<double> &,
                             SVec<double,8> &, SVec<double,3> &, SVec<double,dim> &, double, double, 
                             Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0, bool externalSI=false);

  template<int dim>
  void computeCsValues(SVec<double,dim> &,  SVec<double,16> &, SVec<double,6> &, Vec<double> &,
                       SVec<double,8> &, SVec<double,2> &, Vec<int> &, SVec<double,3> &, double, double,LevelSetStructure *LSS=0);

  template<int dim>
  void computeDynamicLESTerm(DynamicLESTerm *, SVec<double,2> &, SVec<double,3> &, 
                             SVec<double,dim> &, SVec<double,dim> &, 
                             Vec<GhostPoint<dim>*> *ghostPoints=0,
									  LevelSetStructure *LSS=0, bool externalSI=false);

  template<int dim>
  void computeVMSLES_Step1(VMSLESTerm *, SVec<double,dim> &, SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeVMSLES_Step2(SVec<double,1> &, MacroCellSet *, SVec<double,dim> &, SVec<double,dim> &, int);

  template<int dim>
  void computeGalerkinBar_Step1(FemEquationTerm *, BcData<dim> &, GeoState &, SVec<double,3> &,
                                SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeGalerkinBar_Step2(MacroCellSet **, SVec<double,1> &, SVec<double,dim> &,
                                SVec<double,dim> &, int);

  template<int dim>
  void computeMBarAndM_Step1(DynamicVMSTerm *, SVec<double,dim> **, SVec<double,1> **, SVec<double,3> &,
                                      SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeMBarAndM_Step2(MacroCellSet **, SVec<double,1> **, SVec<double,dim> &, SVec<double,dim> &,
                             SVec<double,dim> &, SVec<double,dim> &, int, int);

  template<int dim>
  void computeCsDeltaSq(SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &,
                                   SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &,
                                   Vec<double> &, Vec<double> &, int);

  template<int dim>
  void computeDynamicVMSTerm_Step1(DynamicVMSTerm *, SVec<double,dim> **, SVec<double,3> &,
                                   SVec<double,dim> &, SVec<double,dim> &, Vec<double> &, Vec<double> &,
                                   Vec<double> *, Vec<double> &);
  template<int dim>
  void computeDynamicVMSTerm_Step2(MacroCellSet **, SVec<double,1> **, SVec<double,dim> &, SVec<double,dim> &, int);

  template<int dim>
  void computedWBar_dt(MacroCellSet **, SVec<double,1> **, SVec<double,dim> &, SVec<double,dim> &, int);
  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(FemEquationTerm *fet, BcData<dim> &bcData,
				   GeoState &geoState, SVec<double,3> &X,
				   Vec<double> &ctrlVol, SVec<double,dim> &V,
				   GenMat<Scalar,neq> &A,
											  Vec<GhostPoint<dim>*>* ghostPoints=0, 
											  LevelSetStructure *LSS=0, bool externalSI=false);

  template<int dim, class Scalar, int neq>
  void computeJacobianVolumicForceTerm(VolumicForceTerm *, Vec<double> &,
                                       SVec<double,dim> &, GenMat<Scalar,neq> &) ;

  template<class Scalar, int neq>
  void finishJacobianGalerkinTerm(Vec<double> &, GenMat<Scalar,neq> &) ;

  template<int dim>
  void getExtrapolationValue(Extrapolation<dim>*,SVec<double,dim> &, SVec<double,dim> &, VarFcn*,
                                         BcData<dim>&, GeoState&, SVec<double,3>&);

  template<int dim>
  void applyExtrapolationToSolutionVector(Extrapolation<dim>*, SVec<double,dim> &,
					  SVec<double,dim> &);
  template<int dim>
    void applyBCsToSolutionVector(BcFcn *, BcData<dim> &, SVec<double,dim> &, LevelSetStructure *LSS=0);

  template<int dim>
    void applyBCsToTurbSolutionVector(BcFcn *, BcData<dim> &, SVec<double,dim> &, LevelSetStructure *LSS=0);

  template<int dim>
  void applyBCsToResidual(BcFcn *, BcData<dim> &, SVec<double,dim> &, SVec<double,dim> &, LevelSetStructure *LSS=0);

  template<int dim, class Scalar, int neq>
  void applyBCsToJacobian(BcFcn *, BcData<dim> &, SVec<double,dim> &, GenMat<Scalar,neq> &, LevelSetStructure *LSS=0);

  template<int dim, class Scalar, int neq>
  void applyBCsToH2Jacobian(BcFcn *, BcData<dim> &, SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<int dim, int dim2>
  RectangularSparseMat<double,dim,dim2> *create_NodeBaseddRdXoperators();

  template<int dim, int dim2>
  RectangularSparseMat<double,dim,dim2> *create_ConstantToNodeBaseddRdXoperators();

  template<int dim, int dim2>
  RectangularSparseMat<double,dim,dim2> *create_NodeToConstantBaseddRdXoperators();

  template<int dim, int dim2>
  RectangularSparseMat<double,dim,dim2> *create_ConstantToConstantBaseddRdXoperators();

  template<int dim, int dim2>
  RectangularSparseMat<double,dim,dim2> *create_NodeToFaceBaseddRdXoperators();

  template<int dim, int dim2>
  RectangularSparseMat<double,dim,dim2> *create_NodeToEdgeBaseddRdXoperators();

  template<int dim, int dim2>
  RectangularSparseMat<double,dim,dim2> *create_EdgeBaseddRdXoperators();

  template<int dim, int dim2>
  RectangularSparseMat<double,dim,dim2> *create_FaceBaseddRdXoperators();

  template<class Scalar, int dim>
  SparseMat<Scalar,dim> *createMaskJacobian(int *, MemoryPool *);

  template<class Scalar, int dim>
  MvpMat<Scalar,dim> *createMaskMatVecProd(bool flag = false);

  template<class Scalar, int dim>
  DiagMat<Scalar,dim> *createMaskDiagonal(typename DiagMat<Scalar,dim>::Type, int *);

  template<class Scalar, int dim>
  SparseMat<Scalar,dim> *createMaskILU(int, int, int *);

  template<class Scalar, int dim>
  void computeH1(FluxFcn **, BcData<dim> &, GeoState &, Vec<double> &,
                 SVec<double,dim> &, GenMat<Scalar,dim> &);

  template<int dim, class Scalar, int neq>
  void computeH2(FluxFcn **, RecFcn *, BcData<dim> &, GeoState &, SVec<double,3> &,
		 SVec<double,dim> &, NodalGrad<dim> &, GenMat<Scalar,neq> &);

  template<int dim, class Scalar, int neq>
  void computeH2transpose(FluxFcn **, RecFcn *, BcData<dim> &, GeoState &, SVec<double,3> &,
                          SVec<double,dim> &, NodalGrad<dim> &, GenMat<Scalar,neq> &);

  template<int dim, class Scalar, int neq>
  void computeH2(FluxFcn **fluxFcn, RecFcn *recFcn, BcData<dim> &bcData,
		 GeoState &geoState, SVec<double,3> &X, SVec<double,dim> &V,
		 NodalGrad<dim> &ngrad, 
		 ExactRiemannSolver<dim>& riemann,
		 LevelSetStructure &LSS, 
		 Vec<int> &fluidId, int Nriemann,
		 GenMat<Scalar,neq> &A,
		 SVec<double,dim> &aij, SVec<double,dim> &aji,
		 SVec<double,dim> &bij, SVec<double,dim> &bji, 
		 SVec<double,dim> &betaij, SVec<double,dim> &betaji);

  template<class Scalar, int dim>
  void precomputeRec(RecFcn *, SVec<double,3> &, SVec<double,dim> &,
		     NodalGrad<dim> &, SVec<Scalar,dim> &, SVec<Scalar,dim> &,
		     SVec<Scalar,dim> &, SVec<Scalar,dim> &);

  template<class Scalar, int dim>
  void precomputeRec(RecFcn *recFcn, SVec<double,3> &X, SVec<double,dim> &V,
		     NodalGrad<dim> &ngrad, 
		     LevelSetStructure &LSS,  Vec<int> &fluidId,
		     SVec<Scalar,dim> &aij, SVec<Scalar,dim> &aji,
		     SVec<Scalar,dim> &bij, SVec<Scalar,dim> &bji);

  template<class Scalar, int dim>
  void computeMatVecProdH1(bool *, GenMat<Scalar,dim> &, SVec<double,dim> &,
			   SVec<double,dim> &);

  template<class Scalar, int dim>
  void computeMatVecProdH1transpose(bool *, GenMat<Scalar,dim> &, SVec<double,dim> &,
         SVec<double,dim> &);

  template<class Scalar, int dim>
  void computeMatVecProdH1(bool *, GenMat<Scalar,dim> &, SVec<double,dim> &,
			   SVec<double,dim> &,SVec<double,dim>&, SVec<double,dim>&);
 
  template<class Scalar, int dim>
  void computeMatVecProdH1FarFieldHH(bool *, GenMat<Scalar,dim> &, SVec<double,dim> &,
			             SVec<double,dim> &,Vec<double>&, Vec<double>&);
  
 
  template<class Scalar1, class Scalar2, int dim>
  void computeMatVecProdH2(RecFcn *, SVec<double,3> &, Vec<double> &,
			   GenMat<Scalar1,dim> &, SVec<double,dim> &, SVec<double,dim> &,
			   SVec<double,dim> &, SVec<double,dim> &, SVec<Scalar2,dim> &,
			   NodalGrad<dim, Scalar2> &, SVec<Scalar2,dim> &);

  template<class Scalar1, class Scalar2, int dim>
  void computeMatVecProdH2transposeNew(IoData& iod, SVec<double,3> &X,
                                       Vec<double> &ctrlVol, GenMat<Scalar1,dim> &A,
                                       SVec<double,dim> &aij, SVec<double,dim> &aji,
                                       SVec<double,dim> &bij, SVec<double,dim> &bji,
                                       SVec<Scalar2,dim> &cij, SVec<Scalar2,dim> &cji,
                                       SVec<Scalar2,dim> &dij, SVec<Scalar2,dim> &dji,
                                       NodalGrad<dim, Scalar2> &dpdxj,
                                       SVec<Scalar2,dim> &p, SVec<Scalar2,dim> &prod);

  template<class Scalar, int dim>
  void addToMatVecProdH2transposeLeastSquareNew(SVec<double,3> &X, SVec<double,6> &R,
                                                SVec<Scalar,dim> &ddxt, SVec<Scalar,dim> &ddyt,
                                                SVec<Scalar,dim> &ddzt, NodalGrad<dim, Scalar> &dpdxj,
                                                SVec<Scalar,dim> &prod);

  template<class Scalar, int dim>
  void addToMatVecProdH2transposeGalerkinNew(Vec<double> &ctrlVol,
                                             SVec<Scalar,dim> &ddxt, SVec<Scalar,dim> &ddyt,
                                             SVec<Scalar,dim> &ddzt, NodalGrad<dim, Scalar> &dpdxj,
                                             SVec<Scalar,dim> &prod);

  template<class Scalar1, class Scalar2, int dim>
  void addDiagonalInMatVecProdH2transpose(Vec<double> &, GenMat<Scalar1,dim> &,
                                          SVec<Scalar2,dim> &, SVec<Scalar2,dim> &);

  template<class Scalar1, class Scalar2, int dim>
  void computeMatVecProdH2(FluxFcn **fluxFcn, RecFcn *recFcn, GeoState &geoState,
			   SVec<double,3> &X, Vec<double> &ctrlVol, 
			   ExactRiemannSolver<dim>& riemann,
			   LevelSetStructure &LSS,
			   Vec<int> &fluidId, int Nriemann,
			   GenMat<Scalar1,dim> &A,
			   SVec<double,dim> &aij, SVec<double,dim> &aji,
			   SVec<double,dim> &bij, SVec<double,dim> &bji,
			   SVec<double,dim> &betaij, SVec<double,dim> &betaji,
			   SVec<Scalar2,dim> &p, NodalGrad<dim, Scalar2> &dpdxj,
			   SVec<Scalar2,dim> &prod);

  void convert2(double *Aij, int ln, double *Aij_c) {
    for(int k=0; k<ln*ln; ++k) {;
      Aij_c[k] = Aij[k];
    }
  }

  void convert2(bcomp *Aij, int ln, double *Aij_c) {
    fprintf(stderr, " *** ERROR: Incompatible Types\n");
    exit(-1);
  }

  template<class Scalar1, class Scalar2, int dim>
  void computeMatVecProdH2T(RecFcn *, SVec<double,3> &, Vec<double> &,
                GenMat<Scalar1,dim> &, SVec<double,dim> &, SVec<double,dim> &,
                SVec<double,dim> &, SVec<double,dim> &, SVec<Scalar2,dim> &,
                SVec<Scalar2,dim> &, SVec<Scalar2,dim> &,
                SVec<Scalar2,dim> &, SVec<Scalar2,dim> &);

  template<class Scalar1, class Scalar2, int dim>
  void computeMatVecProdH2Tb(RecFcn *, SVec<double,3> &, Vec<double> &,
                GenMat<Scalar1,dim> &, NodalGrad<dim, Scalar2> &,
                SVec<Scalar2,dim> &, SVec<Scalar2,dim> &, SVec<Scalar2,dim> &);


  // I/O
  template<class Scalar, int dim>
  bool checkIfFileExists(const char *); // without opening file

  template<class Scalar, int dim>
  double readTagFromFile(const char *, int, int *, int *);

  template<class Scalar, int dim>
  void openFileForWriting(const char *, int);

  template<class Scalar, int dim>
  void writeTagToFile(const char *, int, double);

  template<class Scalar, int dim>
  void readVectorFromFile(const char *, int, int, SVec<Scalar,dim> &, Scalar *);

    template<class Scalar>
    void readVectorFromFile(const char *, int, Vec<Scalar> &); //<! for non-state vector, Lei Lei, 22 March 2016

  template<class Scalar, int dim>
  void writeVectorToFile(const char *, int, SVec<Scalar,dim> &, Scalar *);

    template<class Scalar>
    void writeVectorToFile(const char *, int, Vec<Scalar> &, Scalar *); //<! for non-state vector, Lei Lei, 03 Feb 2016

  template<int dim>
  void assignFreeStreamValues2(SVec<double,dim> &, SVec<double,dim> &,
			       SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void assignFreeStreamValues(double *, double *, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void assignPorousWallValues(SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void setNodeBcValue(double*, SVec<double,dim>&);

  template<int dim>
  void setNodeBcValue2(double*, SVec<double,dim>&);

  template<int dim>
  void computeFaceBcValue(SVec<double,dim> &, SVec<double,dim> &);

  template<int dim1, int dim2>
  void computeNodeBcValue(SVec<double,3> &, SVec<double,dim1> &, SVec<double,dim2> &);

  template<int dim>
  void computeNodalForce(PostFcn *, BcData<dim> &, GeoState &, SVec<double,3> &,
		SVec<double,dim> &, Vec<double> &, SVec<double,3> &);

  template<int dim>
  void computeNodalHeatPower(PostFcn*, BcData<dim>&, GeoState&, SVec<double,3>&,
			     SVec<double,dim>&, Vec<double>&);

  template<int dim>
  void computeNodalHeatFluxRelatedValues(PostFcn*, BcData<dim>&,
                                            GeoState&, SVec<double,3>&,
                                            SVec<double,dim>&, Vec<double>&, Vec<double>&, bool);

  template<int dim>
  void computeForceAndMoment(map<int,int> &surfIndexMap, PostFcn *,
                             BcData<dim> &, GeoState &, SVec<double,3> &,
			     SVec<double,dim> &, Vec3D &, Vec3D *, Vec3D *,
			     Vec3D *, Vec3D *, int ,
                             SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX = 0,
                             Vec<double> *genCF = 0);
  template<int dim>
  void computeForceAndMoment(ExactRiemannSolver<dim>&, VarFcn*, map<int,int> &surfIndexMap, PostFcn *,
                             BcData<dim> &, GeoState &, SVec<double,3> &,
                             SVec<double,dim> &, Vec3D &, Vec3D *, Vec3D *,
                             Vec3D *, Vec3D *, int ,
                             SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX = 0,
                             Vec<double> *genCF = 0);
  template<int dim>
  void computeHeatFluxes(map<int,int> &surfIndexMap, PostFcn * , BcData<dim> &,
                                      GeoState &, SVec<double,3> &,
                                      SVec<double,dim> &, double *);
  template<int dim>
  double computeInterfaceWork(PostFcn*, BcData<dim>&, GeoState&, SVec<double,3>&,
			      SVec<double,dim>&, Vec<double>&);
  template<int dim>
  void computeFaceScalarQuantity(PostFcn::ScalarType, PostFcn *, BcData<dim> &, GeoState &,
				 SVec<double,3> &, SVec<double,dim> &, SVec<double,2> &);

  
  template<int dim>
  void computeNodeScalarQuantity(PostFcn::ScalarType, PostFcn *,
				 SVec<double,dim> &, SVec<double,3> &, Vec<double> &);
  template<int dim>
  void computeXP(PostFcn *, SVec<double,dim> &V, SVec<double,3> &X, Vec<double> &XP, int);

  template<int dim,int dimLS>
  void computeNodeScalarQuantity(PostFcn::ScalarType, PostFcn *,
                                 SVec<double,dim> &, SVec<double,3> &,
				 Vec<double> &,Vec<int> &,SVec<double,dimLS>*);

  template<int dim, int dimLS>
    double computeNodeScalarQuantity(PostFcn::ScalarType type, PostFcn *postFcn,
				     SVec<double,dim> &V, SVec<double,3> &X,
				     Vec<int> &fluidId,int i,SVec<double,dimLS>* phi = 0); 
  
  template<int dim>
  void computeForceDerivs(VarFcn *, SVec<double,3> &, SVec<double,dim> &,
                          SVec<double,dim> &, Vec<double> &, SVec<double, 3> **);

  template<int dim>
  void computeForceCoefficients(PostFcn *, Vec3D &, GeoState &, BcData<dim> &, SVec<double,3> &,
                                SVec<double,dim> &, double, Vec3D &, Vec3D &, Vec3D &, Vec3D &,
                                VecSet< SVec<double,3> > *mX = 0 , Vec<double> *genCF = 0);

  // communication

  void markLenNodes(DistInfo &distInfo) { distInfo.setLen(locSubNum, nodes.size()); }
  void markLenEdges(DistInfo &distInfo) { distInfo.setLen(locSubNum, edges.size()); }
  void markLenFaces(DistInfo &distInfo) { distInfo.setLen(locSubNum, faces.size()); }
  void markLenFaceNorms(DistInfo &distInfo) { distInfo.setLen(locSubNum, faces.sizeNorms()); }
  void markLenInletNodes(DistInfo &distInfo) { distInfo.setLen(locSubNum, inletNodes.size()); }

  //! Function to set the length for the object distInfo with the nodes on a Kirchhoff surface
  ///
  /// Function to set the length for the object distInfo with the nodes on a Kirchhoff surface
  ///
  void markLenKirchhoffNodes(IoData &iod, DistInfo &distInfo);


  //! Function to obtain the list of nodes on a surface for the Kirchhoff integral
  ///
  /// This function returns the list of nodes on a surface for the Kirchhoff integral.
  ///
  /// UH (07/2012)
  ///
  const std::set<int>& getKirchhoffNodesList() const { return kirchhoffNodesList; };


  void markLenNull(DistInfo &distInfo) {distInfo.setLen(locSubNum, 0); }
  void makeMasterFlag(DistInfo &);
  void makeEdgeMasterFlag(DistInfo &, CommPattern<int> &edgeNumPat);
  void setChannelNums(SubDTopo &);
  void identifyEdges(CommPattern<int> &);
  void sndNormals(CommPattern<double> &, Vec3D *, double *);
  void rcvNormals(CommPattern<double> &, Vec3D *, double *);
  void sndEdgeData(CommPattern<double> &, double *);
  void rcvEdgeData(CommPattern<double> &, double *);
  void sndEdgeInfo(CommPattern<int> &);
  void rcvEdgeInfo(CommPattern<int> &);

  template<class Scalar>
  void setComLenNodes(int dim, CommPattern<Scalar> &);

  template<class Scalar>
  void setComLenInletNodes(int dim, CommPattern<Scalar> &);

  template<class Scalar>
  void setComLenEdges(int dim, CommPattern<Scalar> &);

  template<class Scalar, int dim>
  void sndData(CommPattern<Scalar> &, Scalar (*)[dim]);

  template<class Scalar, int dim>
  void addRcvData(CommPattern<Scalar> &, Scalar (*)[dim]);

  template<class Scalar, int dim>
  void RcvData(CommPattern<Scalar> &, Scalar (*)[dim]);

  template<int dim>
	  void sndGhostStates(CommPattern<double> &, Vec<GhostPoint<dim>*> &, int);
  template<int dim>
  void sndNumGhostStates(CommPattern<int> &, Vec<GhostPoint<dim>*> &);
  template<int dim>
	  void sndGhostWeights(CommPattern<double> &, Vec<GhostPoint<dim>*> &, int);
  template<int dim>
  void sndGhostTags(CommPattern<int> &, Vec<GhostPoint<dim>*> &);

  template<int dim>
	  void rcvGhostStates(CommPattern<double> &, Vec<GhostPoint<dim>*> &, int);
  template<int dim>
  void rcvNumGhostStates(CommPattern<int> &, Vec<GhostPoint<dim>*> &, VarFcn *varFcn);
  template<int dim>
	  void rcvGhostWeights(CommPattern<double> &, Vec<GhostPoint<dim>*> &, int);
  template<int dim>
  void rcvGhostTags(CommPattern<int> &, Vec<GhostPoint<dim>*> &);

  template<class Scalar, int dim, class OpType>
  void operateRcvData(CommPattern<Scalar> &, Scalar (*)[dim], const OpType &oper);

  template<class Scalar, int dim>
  void sndInletData(CommPattern<Scalar> &, Scalar (*)[dim]);

  template<class Scalar, int dim>
  void addRcvInletData(CommPattern<Scalar> &, Scalar (*)[dim], bool = false);

  template<class Scalar, int dim>
  void sndInletRhsData(CommPattern<Scalar> &, Scalar (*)[dim]);

  template<class Scalar, int dim>
  void addRcvInletRhsData(CommPattern<Scalar> &, Scalar(*)[dim]);

  template<class Scalar, int dim>
  void minRcvData(CommPattern<Scalar> &, Scalar (*)[dim]);

  template<class Scalar, int dim>
  void maxRcvData(CommPattern<Scalar> &, Scalar (*)[dim]);

  template<class Scalar, int dim>
  void maxRcvDataAndCountUpdates(CommPattern<Scalar> &, Scalar (*)[dim],int &, Vec<int> &);

  template<class Scalar1, class Scalar2, int dim1, int dim2>
  void TagPsiExchangeData(CommPattern<Scalar1> &splevel, Scalar1 (*level)[dim1],
                          CommPattern<Scalar2> &sppsi, Scalar2 (*psi)[dim2]);

  template<class Scalar, int dim>
  void sndDiagBlocks(CommPattern<Scalar> &, GenMat<Scalar,dim> &);

  template<class Scalar, int dim>
  void addRcvDiagBlocks(CommPattern<Scalar> &, GenMat<Scalar,dim> &);

  template<class Scalar, int dim>
  void sndDiagInletBlocks(CommPattern<Scalar> &, GenMat<Scalar,dim> &);

  template<class Scalar, int dim>
  void addRcvDiagInletBlocks(CommPattern<Scalar> &, GenMat<Scalar,dim> &);

  template<class Scalar, int dim>
  void checkDiagInletBlocks(CommPattern<Scalar> &, GenMat<Scalar,dim> &);

  template<class Scalar, int dim>
  void sndEdgeData(CommPattern<Scalar> &, Scalar (*)[dim]);

  template<class Scalar, int dim>
  void addRcvEdgeData(CommPattern<Scalar> &, Scalar (*)[dim]);

  template<class Scalar, int dim>
  void sndOffDiagBlocks(CommPattern<Scalar> &, GenMat<Scalar,dim> &);

  template<class Scalar, int dim>
  void addRcvOffDiagBlocks(CommPattern<Scalar> &, GenMat<Scalar,dim> &);

  template<class Scalar, int dim>
  void sndGhostOffDiagBlocks(CommPattern<Scalar> &, GenMat<Scalar,dim> &);

  template<class Scalar, int dim>
  void addRcvGhostOffDiagBlocks(CommPattern<Scalar> &, GenMat<Scalar,dim> &);
  // test

  void testNormals(Vec<Vec3D> &, Vec<double> &, Vec<Vec3D> &, Vec<double> &);

  template<int dim>
  int checkSolution(VarFcn *, SVec<double,dim> &, LevelSetStructure *LSS = 0);

  template<int dim>
  int checkSolution(VarFcn *, SVec<double,dim> &, Vec<int> &, LevelSetStructure *LSS = 0);

  template<int dim>
  int checkSolution(VarFcn *, Vec<double> &, SVec<double,dim> &, Vec<int> &, Vec<int> &);

  template<int dim>
  void restrictionOnPhi(SVec<double,dim> &initial, Vec<int> &fluidId,
                        SVec<double,dim> &restriction, int fluidIdTarget);

  template<int dim>
  void checkFailSafe(VarFcn*, SVec<double,dim>&, SVec<bool,2>&, Vec<int> * = 0);

  template<int dim, int neq>
  int clipSolution(TsData::Clipping, BcsWallData::Integration, VarFcn*,
		   double*, bool*, SVec<double,dim>&, int*, int*, double*);

  template<int dim>
  void checkGradientsSetUp(SVec<double,3> &, SVec<double,dim> &);

  template<int dim>
  void checkGradientsWrite(SVec<double,3> &, NodalGrad<dim> &);

  template<int dim>
  void checkMatVecProd(SVec<double,dim> &, const char *);

#ifdef CXFS
  void printInfo(FILE *);
#endif

  template<int dim>
  void zeroInternalVals(SVec<double, dim> &);

  int *getRotSurfaceOwnership(CommPattern<int> &, map<int,SurfaceData *> &surfaceMap);
  void completeRotateSurfaceOwnership(CommPattern<int> &);

  int *getSlipSurfOwnership(CommPattern<int> &, map<int,SurfaceData *> &surfaceMap);

  void createSlipSurfProjection(int *surfaceOwn, CommPattern<int> &, BCApplier *,
                                SurfaceData *slipSurfaces[]);

  int* getMeshMotionDofType(map<int,SurfaceData*>& surfaceMap, CommPattern<int> &ntP, MatchNodeSet* matchNodes=0);

  int* getEmbeddedALEMeshMotionDofType(map<int,SurfaceData*>& surfaceMap, CommPattern<int> &ntP, MatchNodeSet* matchNodes=0,
                                       DefoMeshMotionData::SlidingSurfaceTreatment slidingSurfaceTreatment=DefoMeshMotionData::Default);

  void completeMeshMotionDofType(int* DofType, CommPattern<int> &ntP);

  void completeEmbeddedALEMeshMotionDofType(int* DofType, CommPattern<int> &ntP);

  void changeSurfaceType(map<int,SurfaceData*>& surfaceMap, map<int,BoundaryData*>& bcMap);
  void markFaceBelongsToSurface(Vec<int> &faceFlag, CommPattern<int> &ntP);
  void completeFaceBelongsToSurface(Vec<int> &faceFlag, Vec<double> &nodeTemp, map<int,SurfaceData*>& surfaceMap, CommPattern<int> &ntP);
 
  template<int dim>
  void zeroMeshMotionBCDofs(SVec<double,dim> &x, int* DofType);

  int *getRotOwn() { return rotOwn; }
  NodeSet &getNodes() { return nodes; }

  template<int dimLS>
  void TagInterfaceNodes(int lsdim, Vec<int> &Tag, SVec<double,dimLS> &Phi, int level, LevelSetStructure *LSS=0);
  template<int dimLS>
  void TagInterfaceNodes(int lsdim, SVec<bool,2> &Tag, SVec<double,dimLS> &Phi, LevelSetStructure *LSS);
  template<int dimLS>
  void pseudoFastMarchingMethod(Vec<int> &Tag, SVec<double,3> &X,
				SVec<double,dimLS> &d2wall, int level, int iterativeLevel,
			        Vec<int> &sortedNodes, int& nSortedNodes, int &firstCheckedNode,
				LevelSetStructure *LSS=0);

  template<int dimLS>
  void FinishReinitialization(Vec<int> &Tag, SVec<double,dimLS> &Psi, int level);

  template<int dim>
  void computeWeightsForEmbeddedStruct(SVec<double,dim> &V, SVec<double,dim> &VWeights,
													Vec<double> &Weights, LevelSetStructure &LSS, SVec<double,3> &X, 
													Vec<int> &init, Vec<int> &next_init, bool externalSI=false);

  template<int dim>
  void computeWeightsForFluidFluid(SVec<double,dim> &V, SVec<double,dim> &VWeights,
				       Vec<double> &Weights, LevelSetStructure *LSS, SVec<double,3> &X, Vec<int> &init, Vec<int> &next_init,
				   Vec<int>& fluidId);

  void computeWeightsLeastSquaresEdgePartForEmbeddedStruct(LevelSetStructure &LSS, 
																			  SVec<double,3> &X, SVec<int,1> &count, SVec<double,10> &R, Vec<int> &init, bool externalSI=false);

  void computeWeightsLeastSquaresNodePartForEmbeddedStruct(SVec<int,1> &count, SVec<double,10> &R);

  void computeWeightsLeastSquaresEdgePartForFF(LevelSetStructure *LSS, 
					       SVec<double,3> &X, SVec<int,1> &count, SVec<double,10> &R, Vec<int> &init,
					       Vec<int>& fluidId);

  void computeWeightsLeastSquaresNodePartForFF(SVec<int,1> &count, SVec<double,10> &R);

  template<int dim>
  void computeWeightsLeastSquaresForEmbeddedStruct(SVec<double,3> &X, SVec<double,10> &R, 
		  SVec<double,dim> &V, Vec<double> &Weights, SVec<double,dim> &VWeights, 
						   LevelSetStructure &LSS, Vec<int> &init, Vec<int> &next_init,
																	NodalGrad<dim>& DX,bool limit,Vec<int>* fluidId, bool externalSI=false);

  template<int dim>
  void computeWeightsLeastSquaresForFluidFluid(SVec<double,3> &X, SVec<double,10> &R, 
						   SVec<double,dim> &V, Vec<double> &Weights, SVec<double,dim> &VWeights, 
						   LevelSetStructure *LSS, Vec<int> &init, Vec<int> &next_init,
					       Vec<int>& fluidId,NodalGrad<dim>& DX, bool limit);

  template<int dim, int dimLS>
  void computeWeightsForEmbeddedStruct(SVec<double,dim> &V, SVec<double,dim> &VWeights, 
                      SVec<double,dimLS> &Phi, SVec<double,dimLS> &PhiWeights, Vec<double> &Weights, 
                      LevelSetStructure &LSS, SVec<double,3> &X, Vec<int> &init, Vec<int> &next_init, Vec<int> &fluidId);
  template<int dimLS>
  void extrapolatePhiV(LevelSetStructure &LSS, SVec<double,dimLS> &PhiV);

  template<int dim>
    void populateGhostPoints(Vec<GhostPoint<dim>*> &ghostPoints, SVec<double,3> &X, SVec<double,dim> &U,
									  NodalGrad<dim, double> &ngrad, VarFcn *varFcn, LevelSetStructure &LSS, bool linRecFSI, Vec<int> &tag);

  template<int dim>
    void populateGhostPoints(Vec<GhostPoint<dim>*> &ghostPoints, SVec<double,3> &X, SVec<double,dim> &U,
									  NodalGrad<dim, double> &ngrad, VarFcn *varFcn, LevelSetStructure &LSS, Vec<int> &fluidId,
									  FemEquationTerm *fet);

  template<int dim>
	  void checkGhostPoints(Vec<GhostPoint<dim>*> &ghostPoints, SVec<double,3> &X, 
							 SVec<double,dim> &U, NodalGrad<dim, double> &ngrad, 
							 VarFcn *varFcn, LevelSetStructure &LSS, Vec<int> &fluidId);
  
  template<int dim, class Scalar, int neq>
    void populateGhostJacobian(Vec<GhostPoint<dim>*> &ghostPoints,SVec<double,dim> &U,FluxFcn** fluxFcn,VarFcn *varFcn,LevelSetStructure &LSS,Vec<int> &tag,GenMat<Scalar,neq>& A);

  template<int dim>
	  void setSIstencil(SVec<double,3> &X, LevelSetStructure &LSS, Vec<int> &fluidId, SVec<double,dim> &U);

  template<int dim>
	  void setFEMstencil(SVec<double,3> &X, LevelSetStructure &LSS, Vec<int> &fluidId, SVec<double,dim> &U);

  template<int dim>
	  void reduceGhostPoints(Vec<GhostPoint<dim>*> &ghostPoints, SVec<double,3> &X);

  template<int dim>
  void computeRiemannWeightsForEmbeddedStruct(SVec<double,dim> &V, SVec<double,dim> &Wstarij,
                      SVec<double,dim> &Wstarji, SVec<double,dim> &VWeights, Vec<double> &Weights,
                      LevelSetStructure &LSS, SVec<double,3> &X);

  template<int dim, int dimLS>
  void computeRiemannWeightsForEmbeddedStruct(SVec<double,dim> &V, SVec<double,dim> &Wstarij,
                      SVec<double,dim> &Wstarji, SVec<double,dim> &VWeights, Vec<double> &Weights,
                      SVec<double,dimLS> &Phi, SVec<double,dimLS> &PhiWeights, 
                      LevelSetStructure &LSS, SVec<double,3> &X, Vec<int> &fluidId0, Vec<int> &fluidId);

  template<int dim>
  void storeGhost(SVec<double,dim> &, SVec<double,dim> &, Vec<double> &);

  template<int dimLS>
  void copyCloseNodes(int lsdim, int level, Vec<int> &Tag,SVec<double,dimLS> &Phi,SVec<double,1> &Psi);
  template<int dimLS>
  void computeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                 NodalGrad<dimLS> &grad,
                                 SVec<double,dimLS> &Phi,SVec<double,1> &Psi);
  template<int dimLS>
  void recomputeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                 NodalGrad<dimLS> &grad,
                                 SVec<double,dimLS> &Phi,SVec<double,1> &Psi);
  template<int dimLS>
  double computeDistanceLevelNodes(int lsdim, Vec<int> &Tag, int level,
                                 SVec<double,3> &X, SVec<double,1> &Psi, SVec<double,dimLS> &Phi);

  template<int dimLS>
  void checkNodePhaseChange(SVec<double,dimLS> &PhiProduct);
  template<int dimLS>
  void getSignedDistance(int lsdim, SVec<double,1> &Psi, SVec<double,dimLS> &Phi);
  template<int dim>
  void storePreviousPrimitive(SVec<double,dim> &V, Vec<int> &fluidId, 
                              SVec<double,3> &X, SVec<double,dim> &Vupdate, 
                              Vec<double> &weight);
  template<int dim>
  void IncreasePressure(double p, VarFcn *vf, SVec<double,dim> &U);
  template<int dim>
  void IncreasePressure(double p, VarFcn *vf, SVec<double,dim> &U, Vec<int> &fluidId);

  template<int dim>
  void checkExtrapolationValue(SVec<double,dim>&,  VarFcn*,
                               BcData<dim>&, GeoState&);

  template<int dim>
  void printVariable(SVec<double,dim>&, VarFcn *vf);

  template<int dim>
  void printInletVariable(SVec<double,dim>&);

  template<int dim>
  void printAllVariable(Vec<int> &, SVec<double,dim>&, int , int);

  template<int dimLS>
  void printPhi(SVec<double,3> &, SVec<double,dimLS>&, int );

  template<class Scalar, int neq>
  void printAllMatrix(GenMat<Scalar,neq> &, int );

    void printDistVecBool(Vec<bool> &, int, bool );

  template<int dim>
  void hardyInterpolationLogMap(SVec<double, dim> ***, SVec<double, dim> **, int, int, int, FullM &, FullM &);

  template<int dim>
   void padeReconstruction(SVec<double, dim> **, SVec<double, dim> **, int *, double *, double, int, int, int, int );
  void buildPadeMatrix(bcomp *, int *, int, double *, bcomp *, int, int, int );

  void buildPadeRhs(bcomp *, int *, int, bcomp *, int, int, int );
  void solveLinearSystem(bcomp *, bcomp *, int );

  void padeSolution(bcomp *, int *, double , int , int , int , double , int , bcomp *, int* , int , double , double );

  void buildDeltaFreq(double *, int , double *, int *);

  void extractElementsRelativeToAComponentAndAMode(double* , bcomp* , int , int , int , int , int , int , double );

  template<int dim>
  void extractElementsRelativeToANode(SVec<double, dim> **, double *, int , int );

  template<int dim>
  void extractElementsRelativeToANodeAndAVector(SVec<double, dim> ***, double *, int , int , int , int );

  template<int dim>
  void snapshotsConstruction(SVec<double, dim> **, bcomp* , int , int , int, int , int , double );

  void multiPointsDeltaFreq(int, double*, int , double* , int *);

  int multiPointsFreq(int , int , double *, int , int );

  void multiPade(bcomp *, int *, double *, bcomp *, bcomp *, int , int , int , double , double , bcomp *, double *);
// Included (MB)
  int computeDerivativeOfControlVolumes(int, double, SVec<double,3> &, SVec<double,3> &, Vec<double> &);
  int computeDerivativeOfControlVolumes(RectangularSparseMat<double,3,1> &, SVec<double,3> &, Vec<double> &);
  int computeTransposeDerivativeOfControlVolumes(RectangularSparseMat<double,3,1> &, Vec<double> &, SVec<double,3> &);
  int computeDerivativeOperatorsOfControlVolumes(SVec<double,3> &, RectangularSparseMat<double,3,1> &);

  void computeDerivativeOfNormals(SVec<double,3> &, SVec<double,3> &, Vec<Vec3D> &,
                                  Vec<Vec3D> &, Vec<double> &, Vec<double> &, Vec<Vec3D> &, Vec<Vec3D> &, Vec<double> &, Vec<double> &);
  void computeDerivativeOfNormals(RectangularSparseMat<double,3,3> &, RectangularSparseMat<double,3,3> &,
                                  SVec<double,3> &, Vec<Vec3D> &, Vec<Vec3D> &);
  void computeTransposeDerivativeOfNormals(RectangularSparseMat<double,3,3> &, RectangularSparseMat<double,3,3> &,
                                           Vec<Vec3D> &, Vec<Vec3D> &, SVec<double,3> &);

  void computeDerivativeOperatorsOfNormals(SVec<double,3> &, RectangularSparseMat<double,3,3> &, RectangularSparseMat<double,3,3> &); 

  void computeDerivativeOfWeightsLeastSquaresEdgePart(SVec<double,3> &, SVec<double,3> &, SVec<double,6> &, SVec<double,6> &);
  void computeDerivativeOfWeightsLeastSquaresEdgePart(RectangularSparseMat<double,3,6> &, SVec<double,3> &, SVec<double,6> &);
  void computeTransposeDerivativeOfWeightsLeastSquaresEdgePart(RectangularSparseMat<double,3,6> &, SVec<double,6> &, SVec<double,3> &);
  void computeDerivativeTransposeOfWeightsLeastSquaresEdgePart(SVec<double,3> &, SVec<double,6> &, SVec<double,6> &, SVec<double,3> &);

  void computeDerivativeOfWeightsLeastSquaresNodePart(SVec<double,6> &, SVec<double,6> &);
  void computeDerivativeOfWeightsLeastSquaresNodePart(RectangularSparseMat<double,6,6> &, SVec<double,6> &);
  void computeTransposeDerivativeOfWeightsLeastSquaresNodePart(RectangularSparseMat<double,6,6> &, SVec<double,6> &);
  void computeDerivativeTransposeOfWeightsLeastSquaresNodePart(SVec<double,6> &, SVec<double,6> &);

  void compute_dRdR(SVec<double,6> &, RectangularSparseMat<double,6,6> &);
  void compute_dRdX(SVec<double,3> &, RectangularSparseMat<double,3,6> &); 

  void computeDerivativeOfWeightsGalerkin(SVec<double,3> &, SVec<double,3> &, SVec<double,3> &,
			      SVec<double,3> &, SVec<double,3> &);
  void computeDerivativeTransposeOfWeightsGalerkin(SVec<double,3> &, SVec<double,3> &, SVec<double,3> &,
			      SVec<double,3> &, SVec<double,3> &);

  template<int dim, class Scalar>
  void computeDerivativeOperatorsOfGradientsLeastSquares(SVec<double,3> &, SVec<double,6> &, SVec<Scalar,dim> &, 
                                                         RectangularSparseMat<double,3,dim> &, 
                                                         RectangularSparseMat<double,3,dim> &, 
                                                         RectangularSparseMat<double,3,dim> &, 
                                                         RectangularSparseMat<double,6,dim> &, 
                                                         RectangularSparseMat<double,6,dim> &, 
                                                         RectangularSparseMat<double,6,dim> &,
                                                         RectangularSparseMat<double,dim,dim> &, 
                                                         RectangularSparseMat<double,dim,dim> &, 
                                                         RectangularSparseMat<double,dim,dim> &); 

  template<int dim, class Scalar>
  void computeDerivativeOfGradientsLeastSquares(
            SVec<double,3> &, SVec<double,3> &,
            SVec<double,6> &, SVec<double,6> &,
            SVec<Scalar,dim> &, SVec<Scalar,dim> &, SVec<Scalar,dim> &,
            SVec<Scalar,dim> &, SVec<Scalar,dim> &);

  template<int dim, class Scalar>
  void computeDerivativeOfGradientsLeastSquares( 
            RectangularSparseMat<double,3,dim> *dddxdX,
            RectangularSparseMat<double,3,dim> *dddydX,
            RectangularSparseMat<double,3,dim> *dddzdX,
            RectangularSparseMat<double,6,dim> *dddxdR,
            RectangularSparseMat<double,6,dim> *dddydR,
            RectangularSparseMat<double,6,dim> *dddzdR,
            RectangularSparseMat<double,dim,dim> *dddxdV,
            RectangularSparseMat<double,dim,dim> *dddydV,
            RectangularSparseMat<double,dim,dim> *dddzdV,
            SVec<double,3> &,
            SVec<double,6> &,
            SVec<double,dim> &,
            SVec<Scalar,dim> &,
            SVec<Scalar,dim> &, 
            SVec<Scalar,dim> &);

  template<int dim, class Scalar>
  void computeTransposeDerivativeOfGradientsLeastSquares( 
            RectangularSparseMat<double,3,dim> *dddxdX,
            RectangularSparseMat<double,3,dim> *dddydX,
            RectangularSparseMat<double,3,dim> *dddzdX,
            RectangularSparseMat<double,6,dim> *dddxdR,
            RectangularSparseMat<double,6,dim> *dddydR,
            RectangularSparseMat<double,6,dim> *dddzdR,
            RectangularSparseMat<double,dim,dim> *dddxdV,
            RectangularSparseMat<double,dim,dim> *dddydV,
            RectangularSparseMat<double,dim,dim> *dddzdV,
            SVec<Scalar,dim> &,
            SVec<Scalar,dim> &, 
            SVec<Scalar,dim> &,
            SVec<double,3> &,
            SVec<double,6> &,
            SVec<double,dim> &);

  template<int dim, class Scalar>
  void computeDerivativeOfGradientsGalerkin(Vec<double> &, Vec<double> &, SVec<double,3> &, SVec<double,3> &,
                SVec<double,3> &, SVec<double,3> &, SVec<double,3> &, SVec<double,3> &,
                SVec<Scalar,dim> &, SVec<Scalar,dim> &, SVec<Scalar,dim> &,
				SVec<Scalar,dim> &, SVec<Scalar,dim> &, SVec<Scalar,dim> &,
				SVec<Scalar,dim> &, SVec<Scalar,dim> &);

  template<int dim>
  void computeDerivativeOfMinMaxStencilValues(SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeDerivativeOfMultiDimLimiter(RecLimiter *, SVec<double,3> &, SVec<double,3> &, Vec<double> &, Vec<double> &,
			      SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &,
			      SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &,
			      SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &,
                  SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeDerivativeOfFiniteVolumeTerm(Vec<double> &, Vec<double> &, FluxFcn**, RecFcn*, BcData<dim>&, GeoState&,
                                           SVec<double,3>&, SVec<double,3>&, SVec<double,dim>&, SVec<double,dim>&,
                                           NodalGrad<dim>&, EdgeGrad<dim>*, double, SVec<double,dim>&);

  template<int dim>
  void computeDerivativeOfFiniteVolumeTerm(RectangularSparseMat<double,dim,dim> *dFluxdddx,
                                           RectangularSparseMat<double,dim,dim> *dFluxdddy,
                                           RectangularSparseMat<double,dim,dim> *dFluxdddz,
                                           RectangularSparseMat<double,3,dim> *dFluxdEdgeNorm,
                                           RectangularSparseMat<double,3,dim> *dFluxdX,
                                           RectangularSparseMat<double,3,dim> *dFluxdFaceNormal,
                                           RectangularSparseMat<double,1,dim> *dFluxdFaceNormalVelX,
                                           RectangularSparseMat<double,dim,dim> *dFluxdUb,
                                           BcData<dim>&, GeoState&, SVec<double,3>&, 
                                           NodalGrad<dim>&, EdgeGrad<dim>*, SVec<double,dim>&, 
                                           SVec<double,dim>&, SVec<double,dim>&, Vec<Vec3D>&, Vec<Vec3D>&, Vec<double>&, SVec<double,dim>&);

  template<int dim>
  void computeTransposeDerivativeOfFiniteVolumeTerm(
                                           RectangularSparseMat<double,dim,dim> *dFluxdddx,
                                           RectangularSparseMat<double,dim,dim> *dFluxdddy,
                                           RectangularSparseMat<double,dim,dim> *dFluxdddz,
                                           RectangularSparseMat<double,3,dim> *dFluxdEdgeNorm,
                                           RectangularSparseMat<double,3,dim> *dFluxdX,
                                           RectangularSparseMat<double,3,dim> *dFluxdFaceNormal,
                                           RectangularSparseMat<double,1,dim> *dFluxdFaceNormalVelX,
                                           RectangularSparseMat<double,dim,dim> *dFluxdUb,
                                           BcData<dim>&, GeoState&, SVec<double,dim>&, 
                                           NodalGrad<dim>&, EdgeGrad<dim>*, SVec<double,3>&,
                                           SVec<double,dim>&, SVec<double,dim>&, SVec<double,dim>&, Vec<Vec3D>&,
                                           Vec<Vec3D>&, Vec<double>&);

  template<int dim>
  void computeDerivativeOperatorsOfFiniteVolumeTerm(Vec<double> &, Vec<double> &, FluxFcn**, RecFcn*, BcData<dim>&, GeoState&,
                                                    SVec<double,3>&, SVec<double,dim>&, 
                                                    NodalGrad<dim>&, EdgeGrad<dim>*, double, 
                                                    RectangularSparseMat<double,3,dim>& dFluxdEdgeNorm,
                                                    RectangularSparseMat<double,3,dim> &dFluxdX,
                                                    RectangularSparseMat<double,dim,dim> &dFluxdddx,
                                                    RectangularSparseMat<double,dim,dim> &dFluxdddy,
                                                    RectangularSparseMat<double,dim,dim> &dFluxdddz,
                                                    RectangularSparseMat<double,3,dim> &dFluxdFaceNormal,
                                                    RectangularSparseMat<double,1,dim> &dFluxdFaceNormalVel,
                                                    RectangularSparseMat<double,dim,dim> &dFluxdUb);

  template<int dim>
  void computeDerivativeOfFiniteVolumeTerm(FluxFcn** fluxFcn, RecFcn* recFcn,
					   BcData<dim>& bcData, GeoState& geoState,
					   SVec<double,3>& X, LevelSetStructure &LSS,
					   bool linRecAtInterface, bool viscSecOrder,
					   Vec<int> &fluidId,
					   ExactRiemannSolver<dim>& riemann,
					   int Nriemann,
					   NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
					   double dMach,
					   SVec<double,dim>& V,
					   SVec<double,dim>& dFluxes);

  template<int dim1, int dim2>
  void computeDerivativeOfNodeBcValue(SVec<double,3> &, SVec<double,3> &, SVec<double,dim1> &, SVec<double,dim1> &, SVec<double,dim2> &);

  template<int dim>
  void computeDerivativeOfNodalForce(PostFcn *, BcData<dim> &, GeoState &, SVec<double,3> &, SVec<double,3> &,
                                     SVec<double,dim> &, SVec<double,dim> &, Vec<double> &, double [3], SVec<double,3> &);

  template<int dim>
  void computeDerivativeOfNodalForce(RectangularSparseMat<double,3,3> *dForcedX,
                                     RectangularSparseMat<double,3,3> *dForcedGradP,
                                     RectangularSparseMat<double,dim,3> *dForcedV,
                                     RectangularSparseMat<double,3,3> *dForcedS,
                                     SVec<double,3> &dX, SVec<double,dim> &dV,
                                     double dS[3], SVec<double,3> &dF,
                                     SVec<double,3> &, SVec<double,3> &); 

  template<int dim>
  void computeTransposeDerivativeOfNodalForce(RectangularSparseMat<double,3,3> *dForcedX,
                                              RectangularSparseMat<double,3,3> *dForcedGradP,
                                              RectangularSparseMat<double,dim,3> *dForcedV,
                                              RectangularSparseMat<double,3,3> *dForcedS,
                                              SVec<double,3> &dF, SVec<double,3> &dGradPSVec,
                                              SVec<double,3> &dX, SVec<double,dim> &dV,
                                              SVec<double,3> dSSVec);                                   

  template<int dim>
  void computeDerivativeOperatorsOfNodalForce(PostFcn *, SVec<double,3> &, SVec<double,dim> &, Vec<double> &, 
                                              RectangularSparseMat<double,3,3> &dForcedX,
                                              RectangularSparseMat<double,3,3> &dForcedGradP,
                                              RectangularSparseMat<double,dim,3> &dForcedV,
                                              RectangularSparseMat<double,3,3> &dForcedS);

  template<int dim>
  void computeDerivativeOfForceAndMoment(map<int,int> &surfIndexMap, PostFcn *, BcData<dim> &, GeoState &, SVec<double,3> &, SVec<double,3> &,
                                         SVec<double,dim> &, SVec<double,dim> &, double [3],
                                         Vec3D &, Vec3D *, Vec3D *, Vec3D *, Vec3D *, int = 0);

  template<int dim>
  void computeDerivativeOfForceAndMoment(RectangularSparseMat<double,3,3> *dFidGradP,
                                         RectangularSparseMat<double,3,3> *dFidX,
                                         RectangularSparseMat<double,dim,3> *dFidV,
                                         RectangularSparseMat<double,3,3> *dFvdX,
                                         RectangularSparseMat<double,dim,3> *dFvdV,
                                         RectangularSparseMat<double,3,3> *dFidS,
                                         RectangularSparseMat<double,3,3> *dMidGradP,
                                         RectangularSparseMat<double,3,3> *dMidX,
                                         RectangularSparseMat<double,dim,3> *dMidV,
                                         RectangularSparseMat<double,3,3> *dMidS,
                                         RectangularSparseMat<double,3,3> *dMvdX,
                                         RectangularSparseMat<double,dim,3> *dMvdV,
                                         SVec<double,3> &dX,
                                         SVec<double,dim> &dV, double dS[3], SVec<double,3> &,
                                         Vec3D *dFi, Vec3D *dMi, Vec3D *dFv, Vec3D *dMv, int hydro = 0);

  template<int dim>
  void computeTransposeDerivativeOfForceAndMoment(RectangularSparseMat<double,3,3> *dFidGradP,
                                                  RectangularSparseMat<double,3,3> *dFidX,
                                                  RectangularSparseMat<double,dim,3> *dFidV,
                                                  RectangularSparseMat<double,3,3> *dFvdX,
                                                  RectangularSparseMat<double,dim,3> *dFvdV,
                                                  RectangularSparseMat<double,3,3> *dFidS,
                                                  RectangularSparseMat<double,3,3> *dMidGradP,
                                                  RectangularSparseMat<double,3,3> *dMidX,
                                                  RectangularSparseMat<double,dim,3> *dMidV,
                                                  RectangularSparseMat<double,3,3> *dMidS,
                                                  RectangularSparseMat<double,3,3> *dMvdX,
                                                  RectangularSparseMat<double,dim,3> *dMvdV,
                                                  SVec<double,3> &dFiSVec, SVec<double,3> &dFvSVec,
                                                  SVec<double,3> &dMiSVec, SVec<double,3> &dMvSVec, SVec<double,3> &dX,
                                                  SVec<double,dim> &dV, SVec<double,3> &dSSVec,
                                                  SVec<double,3> &dGradPSVec, int hydro);

  template<int dim>
  void computeDerivativeOperatorsOfForceAndMoment(map<int,int> &surfIndexMap, PostFcn *, BcData<dim> &, GeoState &, SVec<double,3> &,
                                                  SVec<double,dim> &, Vec3D &, int,
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
                                                  RectangularSparseMat<double,dim,3> &dMvdV);

  template<int dim>
  void computeDerivativeOfGalerkinTerm(FemEquationTerm *, BcData<dim> &, GeoState &,
			   SVec<double,3> &, SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &, double, SVec<double,dim> &);

  template<int dim>
  void computeDerivativeOfGalerkinTerm(RectangularSparseMat<double,3,dim> *, FemEquationTerm *, BcData<dim> &, GeoState &,
			   SVec<double,3> &, SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &, double, SVec<double,dim> &);

  template<int dim>
  void computeTransposeDerivativeOfGalerkinTerm(RectangularSparseMat<double,3,dim> *, SVec<double,dim> &, SVec<double,3> &);

  template<int dim>
  void computeDerivativeOperatorsOfGalerkinTerm(FemEquationTerm *, BcData<dim> &, GeoState &,
			   SVec<double,3> &, SVec<double,dim> &, RectangularSparseMat<double,3,dim> &);


  template<int dim>
  void applyBCsToDerivativeOfResidual(BcFcn *, BcData<dim> &, SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeDerivativeOfNodeScalarQuantity(PostFcn::ScalarDerivativeType, PostFcn *, double [3], SVec<double,dim> &, SVec<double,dim> &, SVec<double,3> &, SVec<double,3> &, Vec<double> &);

  template<int dim, class Scalar>
  void applyBCsToH2Jacobian(BcFcn *, BcData<dim> &, SVec<double,dim> &, GenMat<Scalar,dim> &);

  template<int dim, class Scalar, int neq>
  void applyBCsToJacobianWallValues(BcFcn *, BcData<dim> &, SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<int dim>
  void computeBCsJacobianWallValues(FemEquationTerm *, BcData<dim> &, GeoState &,
				   SVec<double,3> &, SVec<double,dim> &) ;

  template<int dim>
  void computeNodeBCsWallValues(SVec<double,3> &, SVec<double,1> &, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim, class Scalar>
  void applyBCsToProduct(BcFcn *, BcData<dim> &, SVec<double,dim> &, SVec<Scalar,dim> &);

  template<int dim>
  void computeDerivativeOfNodalHeatPower(PostFcn*, BcData<dim>&, GeoState&, SVec<double,3>&, SVec<double,3>&,
			     SVec<double,dim>&, SVec<double,dim>&, double [3], Vec<double>&);

  template<int dim>
  void computeDerivativeOfVolumicForceTerm(VolumicForceTerm *, Vec<double> &, Vec<double> &,
                               SVec<double,dim> &, SVec<double,dim> &, SVec<double,dim> &);

  template<int dim>
  void computeDerivativeOfTimeStep(FemEquationTerm *, VarFcn *, GeoState &,
                                SVec<double,3> &, SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &,
                                Vec<double> &, Vec<double> &, double, TimeLowMachPrec &);

  void checkVec(SVec<double,3> &);

  template<int dim>
  int fixSolution(VarFcn *, SVec<double,dim> &, SVec<double,dim> &, Vec<int>*, int);
  
  template<int dim>
  int fixSolution2(VarFcn *, SVec<double,dim> &, SVec<double,dim> &, Vec<int>*, int);

  template<int dim>
  void getGradP(NodalGrad<dim>&);

  template<int dim>
  void getDerivativeOfGradP(NodalGrad<dim>&);

  template<int dim>
  void getDerivativeOfGradP(RectangularSparseMat<double,dim,3> &dGradPdddx,
                            RectangularSparseMat<double,dim,3> &dGradPdddy,
                            RectangularSparseMat<double,dim,3> &dGradPdddz,
                            SVec<double,dim>& ddVdx,
                            SVec<double,dim>& ddVdy,
                            SVec<double,dim>& ddVdz,
                            SVec<double,3>& dGradPSVec);

  template<int dim>
  void getTransposeDerivativeOfGradP(RectangularSparseMat<double,dim,3> &dGradPdddx,
                                     RectangularSparseMat<double,dim,3> &dGradPdddy,
                                     RectangularSparseMat<double,dim,3> &dGradPdddz,
                                     SVec<double,3>& dGradPSVec,
                                     SVec<double,dim>& ddVdx,
                                     SVec<double,dim>& ddVdy,
                                     SVec<double,dim>& ddVdz);

  template<int dim>
  void computeDerivativeOperatorOfGradP(RectangularSparseMat<double,dim,3> &dGradPdddx,
                                        RectangularSparseMat<double,dim,3> &dGradPdddy,
                                        RectangularSparseMat<double,dim,3> &dGradPdddz);

  template<int dim>
  void computePrdtWCtrlVolRatio(SVec<double,dim> &, SVec<double,dim> &, Vec<double> &, GeoState &);

  template<int dimLS>
  void computePrdtPhiCtrlVolRatio(SVec<double,dimLS> &, SVec<double,dimLS> &, Vec<double> &, GeoState &);

//  void getTriangulatedSurfaceFromFace( SVec<double,3> &);
  void getTriangulatedSurfaceFromFace( TriangulatedSurface* );
  void getTriangulatedSurfaceFromFace( SVec<double,3> &,  TriangulatedSurface* );

//  void printTriangulatedSurface();
  void getGhostNodes(double*, int*, int&);
  bool isINodeinITet(Vec3D, int, SVec<double,3>&);
  void  localCoord(Vec3D, int, SVec<double, 3>&, Vec3D&);
  int* getNeiElemOfNode(int, int, int&);
  int* getNeiNodeOfNode(int, int&);
  void getNodeCoords(int, SVec<double,3> &, double&, double&, double&);
  void updateNodeTag(SVec<double,3>&, LevelSetStructure &, Vec<int>&, Vec<int>&);

  template<int dim>
  void computeCVBasedForceLoad(int, int, GeoState &,SVec<double,3>&, double (*)[3], int, LevelSetStructure&, double pInfty, 
                               SVec<double,dim> &Wstarij,SVec<double,dim> &Wstarji,SVec<double,dim> &V, 
                               Vec<GhostPoint<dim>*> *ghostPoints, PostFcn *postFcn,NodalGrad<dim,double> &ngrad, VarFcn *vf, Vec<int> *fid);
  template<int dim>
  void computeEmbSurfBasedForceLoad(IoData &iod, int, int, SVec<double,3>&, double (*)[3], int, int, int (*)[3], Vec<Vec3D>&, 
				    LevelSetStructure&, double pInfty, 
				    SVec<double,dim> &Wstarij,SVec<double,dim> &Wstarji,SVec<double,dim> &V, 
				    Vec<GhostPoint<dim>*> *ghostPoints, PostFcn *postFcn, NodalGrad<dim,double> &ngrad, 
				    VarFcn *vf, Vec<int>* fid);

  template<int dim>
  void computeEmbSurfBasedForceLoad_e(IoData &iod, int, int, SVec<double,3>&, double (*)[3], int, int, int (*)[3], Vec<Vec3D>&, 
												  LevelSetStructure&, double pInfty, SVec<double,dim> &V, 
												  Vec<GhostPoint<dim>*> *ghostPoints, PostFcn *postFcn, NodalGrad<dim,double> &ngrad, 
												  VarFcn *vf, Vec<int>* fid,
												  int** stNodeDir, double** stX1, double** stX2);

  // d2d
  template<int dim,int dimLS>
  void computeEMBNodeScalarQuantity(SVec<double,3> &X, SVec<double,dim> &V,
                         	    PostFcn *postFcn, VarFcn *varFcn,
				    Vec<int> &fluidId, SVec<double,dimLS>* phi,
				    double (*Qnty)[3], int sizeQnty, int numStructElems, int (*stElem)[3],
				    Vec<Vec3D>& Xstruct, LevelSetStructure &LSS,
				    double pInfty,
				    Vec<GhostPoint<dim>*> *ghostPoints,
				    NodalGrad<dim, double> &ngrad);

  template<int dim,int dimLS>
	  void computeEMBNodeScalarQuantity_e(SVec<double,3> &X, SVec<double,dim> &V, 
													  PostFcn *postFcn, VarFcn *varFcn, 
													  Vec<int> &fluidId, SVec<double,dim> &Wextij, SVec<double,dimLS>* phi,
													  double (*Qnty)[3], int sizeQnty, int numStructElems, int (*stElem)[3],
													  Vec<Vec3D>& Xstruct, LevelSetStructure &LSS,
													  double pInfty, 
													  Vec<GhostPoint<dim>*> *ghostPoints,
													  NodalGrad<dim, double> &ngrad);

  template<int dim>
  void computeEMBNodeScalarQuantity_step1(SVec<double,3> &X, SVec<double,dim> &V,
														int numStructElems, int (*stElem)[3],
														Vec<Vec3D>& Xstruct, LevelSetStructure &LSS,
														int** stNodeDir, double** stX1, double** stX2,
														bool rebuildTree=true);


  template<int dim>
  void computeEMBNodeScalarQuantity_step2(SVec<double,3> &X, SVec<double,dim> &V, 
														PostFcn *postFcn, VarFcn *varFcn, 
														Vec<int> &fluidId, 
														double (*Qnty)[3], int sizeQnty, int numStructElems, int (*stElem)[3],
														Vec<Vec3D>& Xstruct, LevelSetStructure &LSS,
														double pInfty, 
														Vec<GhostPoint<dim>*> *ghostPoints,
														NodalGrad<dim, double> &ngrad,
														int** stNodeDir, double** stX1, double** stX2);
  
  template<int dim>
  void computederivativeEmbSurfBasedForceLoad(IoData &iod, int forceApp, int order, SVec<double,3> &X,
					      double (*dFs)[3], int sizeFs, int numStructElems, int (*stElem)[3], 
					      Vec<Vec3D>& Xstruct, Vec<Vec3D>& dXstruct, LevelSetStructure &LSS, 
					      double pInfty, double dpInfty, 
					      SVec<double,dim> &Wstarij, SVec<double,dim> &Wstarji, 
					      SVec<double,dim> &V, SVec<double,dim> &dV_, 
					      Vec<GhostPoint<dim>*> *ghostPoints, PostFcn *postFcn, 
					      NodalGrad<dim, double> &gradV, NodalGrad<dim, double> &graddV, VarFcn* vf, Vec<int>* fid);
    
  template<int dim>
  void computeRecSurfBasedForceLoad(int, int, SVec<double,3>&, double (*)[3], int, LevelSetStructure&, double pInfty, 
                                    SVec<double,dim> &Wstarij,SVec<double,dim> &Wstarji,SVec<double,dim> &V, 
			                              Vec<GhostPoint<dim>*> *ghostPoints, PostFcn *postFcn, VarFcn *vf, Vec<int>* fid);
  void addLocalForce(int, Vec3D, double, double, double, LevelSetResult&, LevelSetResult&,
                     LevelSetResult&, double(*)[3]); //not used.
  void sendLocalForce(Vec3D, LevelSetResult&, double(*)[3]);

  void computeCharacteristicEdgeLength(SVec<double,3>&, double&, double&, double&, int&, const double, const double, const double, const double, const double, const double);
  double scalarNormalExtrap(double*, Vec3D, Vec3D, int, SVec<double,3> &, bool);
  bool insideOutside(double*, const double, const double, const double, const double,
                     const double, const double);

  void getMeshInBoundingBox(SVec<double,3> &, const double, const double,
                              const double, const double,
                              const double, const double,
                              int*, int&, int*,
                              int&, int (*)[4]);

  EdgeSet &getEdges() { return edges; }
  Connectivity *getNodeToNode() { if(!NodeToNode) NodeToNode = createEdgeBasedConnectivity();  return NodeToNode; }
  Connectivity *getNodeToSubD() { if(!NodeToSubD) NodeToSubD = createNodeToSubDomainConnectivity();  return NodeToSubD; }
  int findFarfieldNode();


  template<int dim>
  void blur(SVec<double,dim> &U, SVec<double,dim> &U0,Vec<double>& weight);

  void solicitFluidIdFS(LevelSetStructure &LSS, Vec<int> &fluidId, SVec<bool,4> &poll);

  template<int dimLS>
  void updateFluidIdFS2(LevelSetStructure &LSS, SVec<double,dimLS> &PhiV, SVec<bool,3> &poll, Vec<int> &fluidId, bool *masterFlag); //KW:obsolete

  template<int dim, int dimLS>
  void debugMultiPhysics(LevelSetStructure &LSS, SVec<double,dimLS> &PhiV, Vec<int> &fluidId, SVec<double,dim> &U);

  template<int dim, class Obj>
  void integrateFunction(Obj* obj,SVec<double,3> &X,SVec<double,dim>& V, void (Obj::*F)(int node, const double* loc,double* f),
                         int npt);

  

  template<int dim> 
  void interpolateSolution(SVec<double,3>& X, SVec<double,dim>& U, 
                           const std::vector<Vec3D>& locs, double (*sol)[dim],
			   int* status,int* last,int* nid,
			   LevelSetStructure* LSS = 0, Vec<GhostPoint<dim>*>* ghostPoints = 0,
                           VarFcn *varFcn = 0, bool assumeCache = false, Vec<int> *fluidId = 0);
  
  template<int dim>
  void interpolatePhiSolution(SVec<double,3>& X, SVec<double,dim>& U,
                           const std::vector<Vec3D>& locs, double (*sol)[dim],
                           int* status,int* last,int* nid, bool assumeCache = false); 
  void createHigherOrderMultiFluid();
  void createHigherOrderFSI();

  void assignErrorHandler(ErrorHandler* in);

  // Functions to compute the error (that is, the difference between two state vectors)
  template <int dim>
    void computeL1Error(bool* nodeFlag,SVec<double,dim>& U, SVec<double,dim>& Uexact, Vec<double>& vol, double error[dim], LevelSetStructure* = NULL);

  template <int dim>
    void computeL2Error(bool* nodeFlag,SVec<double,dim>& U, SVec<double,dim>& Uexact, Vec<double>& vol, double error[dim], LevelSetStructure* = NULL);
  
  template <int dim>
    void computeLInfError(bool* nodeFlag,SVec<double,dim>& U, SVec<double,dim>& Uexact, double error[dim], LevelSetStructure* = NULL);

  HigherOrderMultiFluid* getHigherOrderMF() { return higherOrderMF; }
  HigherOrderFSI* getHigherOrderFSI() { return higherOrderFSI; }

  template <int dim>
    void computeHHBoundaryTermResidual(BcData<dim> &bcData,SVec<double,dim> &U,Vec<double>& res,
				       VarFcn* vf);

  void maskHHVector(Vec<double>& hh);
  template<int dim, class Scalar, int neq>
    void computeJacobianFiniteVolumeTermHH(FluxFcn **fluxFcn, BcData<dim> &bcData,
					   GeoState& geoState,
					   Vec<double> &ctrlVol,
					   SVec<double,dim> &U, 
					   GenMat<Scalar,neq> &A,
					   VarFcn* vf);

  void setFarFieldMask(Vec<double>& ffMask, Vec<double>& neighborMask);
  void setWallMask(Vec<double>& wallMask, Vec<double>& neighborMask); 
 
  void attachTriangulatedInterfaceLSS(LevelSetStructure*);
  
  Elem* searchPoint(Vec3D Xp, SVec<double,3>& X);

};
//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <SubDomain.C>
#endif

#endif
