#ifndef _GAPPY_H_
#define _GAPPY_H_

class IoData;
class Domain;
#include <Communicator.h>
#include <DistVector.h>
#include <VectorSet.h>
#include <ParallelRom.h>
#include <TsDesc.h>
#include <set>
#include <vector>
#include <algorithm>
#include <iterator>
#include <NonlinearRom.h>
#include <boost/unordered_map.hpp>

template <int dim>
class VecSetArray {
	public:
		// use a to access pointer to vecset, use [] to access the vecset itself
		VecSet< DistSVec<double,dim> > *(a[2]);	// B[0] = *(a[0])
		VecSet< DistSVec<double,dim> >& operator [](int i) { return *(a[i]); };
};

template <int dim>
struct SubDomainData { 
	// need this struct because we are creating vectors of them (unknown number of
	// entries
	double (*a) [dim];
	double * operator [](int i) {return a[i];};
	SubDomainData& operator=(double (*b)[dim] ) {
		a = b;
		return *this;
	};

	SubDomainData(double (*b)[dim] ) {
		a = b;
	};
	SubDomainData() { };
};

template <int dim>
struct VecSubDomainData { 
	// need this struct because the error vector is a 4-D array and complicated
	// 
	// locError[ nPodBasis ][ nRhs[iPodBasis] ][ locNodeNum ][ iDim ]

	std::vector< SubDomainData<dim> > a [2];	// array of two vectors, each vector is nRhs long

	std::vector< SubDomainData<dim> > &operator [](int i) {return a[i];};
	void resize (int size) {a[0].resize(size); a[1].resize(size); };	// within constructor, specify maximum size
	void resize (int *size) {a[0].resize(size[0]); a[1].resize(size[1]); };	// within constructor, specify maximum size
};

template <typename Scalar, int size>
class StaticArray {	//used for the value of a map
	private:
		Scalar a [size];
	public:
		// use a to access pointer to vecset, use [] to access the vecset itself
		StaticArray() { for (int i=0; i<size; ++i) a[i] = static_cast<Scalar>(0);}
		StaticArray(Scalar b [size]){ for (int i=0; i<size; ++i) a[i] = b[i];}
		StaticArray(const StaticArray &other){ for (int i=0; i<size; ++i) a[i] = other.a[i];}
		StaticArray& operator=(const StaticArray &other){ for (int i=0; i<size; ++i) a[i] = other.a[i]; return *this;}
		Scalar& operator[] (int i){ return a[i];}
		const Scalar &operator[] (int i) const{ return a[i];}
		bool operator<(const StaticArray &other) const{ 
			int index = -1;
			for (int i = 0; i < size; ++i) {
				if (a[i] != other[i]) {	// if different, return less than
					index = i;
					break;
				}
			}
			if (index == -1)	// two are equal
				return false;
			else
				return a[index] < other[index];
		}
};

template <int dim>
class GappyPreprocessing : public NonlinearRom<dim> {

protected:
        GappyConstructionData* gappyIO;
        ApproximatedMetricData* approxMetricData;

	typedef VecSet< DistSVec<double,dim> > SetOfVec;
	bool debugging; 	// debugging flag
	bool twoLayers; 	// debugging flag
	int includeLiftFaces; 	// if the reduced mesh should include lift faces

	std::vector< ParallelRom<dim> *> parallelRom;	// object for all parallel operations
	void setUpPodResJac(int);
  void setUpBasisBasisProducts();
	void setUpPseudoInverse();

	int nSampleNodes;	// number of parent sample globalNodes
	int errorBasis[2];	// basis indices for determining error
	void newNeighbors();	// add unique neighbors of a node

	Domain &domain;
	Communicator *com;
	IoData *ioData;	
	TsInput *input;	
	DistGeoState *geoState;
	DistSVec<double,3> &X;

  DistVec<double>* wallMask;
  DistVec<double>* wallNeighborsMask;
  DistVec<double>* targetRegionMask;

	GeoSource *geoSourceTmp;
	TsDesc<dim> *tsDescTmp;
  DistBcData<dim> *bcDataTmp;
  VarFcn *varFcnTmp;
  PostOperator<dim> *postOp;

	const int residual;	// refer to residual as 0
	const int jacobian;
	int nPod [2];	// nPod[0] = nPodRes, nPod[1] = nPodJac
	int nPodMax;
	int dimGreedy; // number of basis vectors used for greedy

	std::vector<int> cpuSample, locSubSample, locNodeSample;
	std::vector<int> globalSampleNodes, reducedSampleNodes;
	boost::unordered_map<int, int > globalSampleNodeRankMap, reducedSampleNodeRankMap;
		//key: node #, value: rank of sample node (position in _SampleNodeSet)

	int nPodBasis;	// # of unique pod bases for residual/jac (either 1 or 2)
	VecSetArray<dim> pod;	// pod bases for residual and jacobian
	SetOfVec podRes, podJac;
	VecSetArray<dim> podHat;	// restricted pod bases 
	SetOfVec podHatRes, podHatJac;
	VecSetArray<dim> error;	// error vectors
	SetOfVec errorRes, errorJac;
        VecSetArray<dim> errorHat;
        SetOfVec errorHatRes, errorHatJac;

	int nPodState;

	// greedy data
	int nRhsMax, nGreedyIt, handledNodes;
	int *(nRhsGreedy [2]);
	int handledVectors [2];
	int nRhs [2];	// nRhs at a given greedy iteration
	int *nodesToHandle;	// how many globalNodes are handled at each greedy iteration
	VecSubDomainData<dim> locError;

	void parallelLSMultiRHSGap(int iPodBasis, double **lsCoeff);
        void serialLSMultiRHSGap(int iPodBasis, double **lsCoeff);
        void probabilisticLSMultiRHSGap(int iPodBasis, double **lsCoeff);
	bool onlyInletOutletBC;

	// greedy functions

	virtual void determineSampleNodes();
	void greedyIteration(int greedyIt);
        void readGreedyData(int iCluster, bool& breakloop);
	virtual void setUpGreedy(int iCluster);
	void computeNodeError(bool *locMasterFlag, int locNodeNum, double &nodeError);
	void findMaxAndFillPodHat(const double myMaxNorm, const int locSub, const int locNode, const int globalNode);
	void makeNodeMaxIfUnique(double nodeError, double &myMaxNorm, int iSub, int locNodeNum, int &locSub, int &locNode, int &globalNode);
	void getSubDomainError(int iSub);	// computes locError for iSub subdomain
	void leastSquaresReconstruction();
	void subDFindMaxError(int iSub, bool onlyInletBC, double &myMaxNorm, int &locSub, int &locNode, int &globalNode);	// computes locError for iSub subdomain
        void maskError();
        void setUpSampledNodeTargetRegionMask();

	// least squares parameters
	
	void initializeLeastSquares();

	// distribution info

	int numLocSub, nTotCpus, thisCPU; 
	DistInfo &nodeDistInfo;
	SubDomain** subD; 
	
	// mesh construction
	// each of these arrays has nSampleNodes elements
  std::vector <int> *globalNodes;	// globalNodes[iSampleNode][iNode] is the global node number of the iNode in the iSampleNode island 
  std::vector <int> *cpus;	// globalNodes[iSampleNode][iNode] is the global node number of the iNode in the iSampleNode island 
  std::vector <int> *locSubDomains;	// globalNodes[iSampleNode][iNode] is the global node number of the iNode in the iSampleNode island 
  std::vector <int> *localNodes;	// globalNodes[iSampleNode][iNode] is the global node number of the iNode in the iSampleNode island 
  std::vector <double> *(nodesXYZ [3]);	// nodesXYZ[iXYZ][iSampleNode][iNode] is the iXYZ coordinate of the iNode in the iSampleNode island
  std::vector <int> *elements;		// elements[iSampleNode][iEle] is the global element number of the iEle element in the iSampleNode island 
  std::vector <int> *(elemToNode [4]);	// elemToNode[iNode][iSampleNode][iEle] is the global node number of the iNode attached to the iEle element of the iSampleNode island 
  int *nodeOffset;
  int *totalNodesCommunicated;
  int *totalEleCommunicated;
  std::vector< int > *(bcFaces [2][3]);	// boundary faces. bcfaces[iSign][whichNode][BCtype][iFace] returns the global node number of whichNode on the iFace face corresponding to iSign/BCtype. iSign = 0 if the BC definition is negative, and iSign = 1 if positive. BCtype can be found in BcDefs.h
  std::vector< int > *(bcFaceSurfID [2]);	// codes for the above boundary faces. bcFaceSurfID[iSign][BCtype][iFace] returns the surfaceID of the iFace face corresponding to iSign/BCtype
  std::set<StaticArray <int, 3> > bcFacesInfo;	// {iCPU,iSub,iFace}

  boost::unordered_map<int, StaticArray <double, 3> > nodesXYZmap;	// key: global node #, values: x, y, z
  boost::unordered_map<int, int > globalNodeToCpuMap;	// key: global node #, values: x, y, z
  boost::unordered_map<int, int > globalNodeToLocSubDomainsMap;	// key: global node #, values: x, y, z
  boost::unordered_map<int, int > globalNodeToLocalNodesMap;	// key: global node #, values: x, y, z
  boost::unordered_map <int, StaticArray <int, 4> > elemToNodeMap;	// key: global elem #, values: global node #s
  boost::unordered_map <int, std::string > boundaryConditionsMap;	// mapping between BC numbers in BcDef.h and Sower's identification
  //above maps have been defined for vector entries [iSampleNode][0:j]
  int numFullNodes, nReducedNodes;	// number of nodes in full and reduced meshes
  bool outputOnlineMatricesFull, outputOnlineMatricesSample;
  bool initializeLeastSquaresDone; 

  // KTC!!! then, when outputting the TOP file, need another key that maps global
  // node # to reduced mesh node #... this mapping will be different for
  // each island!

  // also create a pointer of vectors to handle the faces that might be
  // boundary conditions
  bool cleanTempFiles;

  void computeXYZ(int iSub, int iLocNode, double *xyz);
  virtual void buildRemainingMesh();
  void computeBCFaces(bool);
  bool checkFaceInMesh(FaceSet& currentFaces, const int iFace, const int iSub, const int *locToGlobNodeMap);
  bool checkFaceAlreadyAdded(const int cpuNum, const int
        iSub, const int iFace);
  bool checkFaceContributesToLift(FaceSet& currentFaces, const int iFace,
	const int iSub, const int *locToGlobNodeMap);
  void addFaceNodesElements(FaceSet& currentFaces, const int iFace,
	const int iSub, const int *locToGlobNodeMap);
  void addNodesOnFace(FaceSet& currentFaces, const int iFace,
	const int iSub, const int *locToGlobNodeMap, int *locNodeNums);
  void addElementOfFace(FaceSet& currentFaces, const int iFace,
	const int iSub, const int *locToGlobNodeMap, const int *globalNodeNums);
  virtual void addSampleNodesAndNeighbors();
  void addNeighbors(int iSampleNodes, int startingNodeWithNeigh = 0);
  void communicateAll();
  void initialize();
  void defineMaps();
  void communicateBCFaces();
  virtual void outputTopFile(int);
  virtual void outputSampleNodes(int);
  void outputSampleNodesGeneral(const std::vector<int> &sampleNodes, const char *sampleNodeFile);

  // A and B matrices functions

  // pseudo-inverse functions
  double **(podHatPseudoInv [2]);	// each dimension: (nSampleNode*dim) x nPod[i]
  virtual void computePseudoInverse();
  void parallelPseudoInverse(int iPodBasis);
  void serialPseudoInverse(int iPodBasis);
  //void computePseudoInverseRHS();
  void checkConsistency();
  SetOfVec pseudoInvRhs;

  double **podTpod;	// stores phiJ^T * phiR
  virtual void computePodTPod();  // compute phiJ^T * phiR
 
  double **(onlineMatrices [2]);	// dimension: (nSampleNode*dim) x nPod[1]
		// onlineMatrices[0] is related to the residual: 
		// 		pod[1]^Tpod[0] * podHatPseudoInv[0]^T
		// onlineMatrices[1] is related to the jacobian: 
		// 		podHatPseudoInv[1]^T

  virtual void assembleOnlineMatrices();
  void outputOnlineMatrices(int);
  virtual void outputOnlineMatricesGeneral( int iCluster, 
	int numNodes, const boost::unordered_map<int,int> &sampleNodeMap, const
	std::vector<int> &sampleNodeVec);
  virtual void outputReducedToFullNodes();
  //virtual void outputStateReduced();
  virtual void outputApproxSnapsReduced(int iCluster);
  virtual void outputInitialConditionReduced();
  virtual void outputMatchStateReduced();
  virtual void outputMultiSolutionsReduced();
  virtual void outputClusterCentersReduced();
  virtual void outputLocalStateBasisReduced(int);
  virtual void outputLocalReferenceStateReduced(int);
  virtual void outputWallDistanceReduced();
  virtual void outputDisplacementReduced();
  virtual void outputShapeDerivativeReduced();

  void outputReducedSVec(const DistSVec<double,dim> &distSVec, FILE* myOutFile, double tag);
  void outputReduced3DSVec(const DistSVec<double,3> &distSVec, FILE* myOutFile, double tag);
  void outputReducedVec(const DistVec<double> &distVec, FILE* outFile , int iVector);
	//void determineFileName(const char *fileNameInput, const char
	//		*currentExtension, const char *(&fileNameBase), const char
	//		*(&fileNameExtension));

  int unionOfSampleNodes; // = -1 (makes the code a bit easier to read)
  std::set<int> globalSampleNodesUnionSet; // union of sample nodes from each cluster
  std::set<int> globalSampleNodesUnionSetForApproxMetricState; // union of sample nodes from each cluster for approximated metric
  std::vector<int> globalSampleNodesUnion; // union of sample nodes from each cluster (as a vector)
  std::vector<int> globalSampleNodesUnionForApproxMetricState; // union of sample nodes from each cluster for approx metric (as a vector)
  std::vector<std::vector<int> > globalSampleNodesForCluster; // stores sampled nodes for each of the clusters   
  void constructApproximatedMetric(const char* type, int iCluster = -1, std::vector<std::vector<double> >* corrMat=NULL);
  void computeCorrelationMatrixEVD(std::vector<std::vector<double> >* corrMat=NULL);
  void computePseudoInverseMaskedSnapshots(const char* type, int iCluster = -1);
  void computeMaskedSnapshots(const char* type, int iCluster = -1);
  void computePseudoInverseTranspose();
  void computeApproximatedMetricLowRankFactor();
  void computeApproxMetricNonlinearNNLS(int iCluster);
  void outputApproxMetricLowRankFactorFullCoords(const char* type, int iCluster = -1);
  void outputApproxMetricLowRankFactorReducedCoords(const char* type, int iCluster = -1);
  void testInnerProduct(const char *);
  double** lowRankModes;
  double* lowRankApproxMetricEigenvalues;
  int numEigen;
  int numSnapsForApproxMetric;
  int nApproxMetricSampleNodes;
  SetOfVec pseudoInverseMaskedSnapsTrans;
  SetOfVec snapHatApproxMetric;
  void setSampleNodes(int);
  void formMaskedNonlinearROBs();
  void reinitializeMapsForSampleNodes();
  //void outputMaskedNonlinearROBs(int, const boost::unordered_map<int,int> &, const std::vector<int> &);
  //void readMaskedNonlinearROBs( );
  void initializeLeastSquaresPseudoInv(int);
  void formReducedSampleNodeMap();

  bool surfaceMeshConstruction;

public:
	GappyPreprocessing(Communicator *, IoData &, Domain &, DistGeoState *);
	~GappyPreprocessing();
	virtual void buildReducedModel();	// build all offline info (do everything)

};
#endif
