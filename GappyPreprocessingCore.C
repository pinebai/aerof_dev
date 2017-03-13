#include <GappyPreprocessing.C>

#define INSTANTIATION_HELPER(dim)\
template \
GappyPreprocessing<dim>::GappyPreprocessing(Communicator *_com, IoData &_ioData, Domain\
    &dom, DistGeoState *_geoState);\
\
template \
GappyPreprocessing<dim>::~GappyPreprocessing();\
\
template \
void GappyPreprocessing<dim>::initialize();\
\
template \
void GappyPreprocessing<dim>::buildReducedModel();\
\
template \
void GappyPreprocessing<dim>::constructApproximatedMetric(const char *, int, std::vector<std::vector<double> >* corrMat);\
\
template \
void GappyPreprocessing<dim>::computeCorrelationMatrixEVD(std::vector<std::vector<double> >* corrMat);\
\
template \
void GappyPreprocessing<dim>::computePseudoInverseMaskedSnapshots(const char *, int);\
\
template \
void GappyPreprocessing<dim>::computeMaskedSnapshots(const char *, int);\
\
template \
void GappyPreprocessing<dim>::computePseudoInverseTranspose();\
\
template \
void GappyPreprocessing<dim>::computeApproximatedMetricLowRankFactor();\
\
template \
void GappyPreprocessing<dim>::outputApproxMetricLowRankFactorFullCoords(const char*, int);\
\
template \
void GappyPreprocessing<dim>::outputApproxMetricLowRankFactorReducedCoords(const char*, int);\
\
template \
void GappyPreprocessing<dim>::testInnerProduct(const char *snapshotType);\
\
template \
void GappyPreprocessing<dim>::setSampleNodes(int iCluster);\
\
template \
void GappyPreprocessing<dim>::setUpPodResJac(int iCluster);\
\
template \
void GappyPreprocessing<dim>::setUpBasisBasisProducts();\
\
template \
void GappyPreprocessing<dim>::setUpGreedy(int iCluster);\
\
template \
void GappyPreprocessing<dim>::findMaxAndFillPodHat(const double myMaxNorm, const int\
    locSub, const int locNode, const int globalNode);\
\
template \
void GappyPreprocessing<dim>::determineSampleNodes();\
\
template \
void GappyPreprocessing<dim>::greedyIteration(int greedyIt);\
\
template \
void GappyPreprocessing<dim>::initializeLeastSquares();\
\
template \
void GappyPreprocessing<dim>::initializeLeastSquaresPseudoInv(int numRhs);\
\
template \
void GappyPreprocessing<dim>::makeNodeMaxIfUnique(double nodeError, double\
    &myMaxNorm, int iSub, int locNodeNum, int &locSub, int &locNode, int\
    &globalNode);\
\
template \
void GappyPreprocessing<dim>::computeNodeError(bool *locMasterFlag, int locNodeNum, double &nodeError);\
\
template \
void GappyPreprocessing<dim>::getSubDomainError(int iSub);\
\
template \
void GappyPreprocessing<dim>::leastSquaresReconstruction();\
\
template void GappyPreprocessing<dim>::subDFindMaxError(int iSub, bool\
    onlyInletOutletBC, double &myMaxNorm, int &locSub, int &locNode, int\
    &globalNode);\
\
template \
void GappyPreprocessing<dim>::parallelLSMultiRHSGap(int iPodBasis, double **lsCoeff);\
\
template \
void GappyPreprocessing<dim>::buildRemainingMesh();\
\
template \
void GappyPreprocessing<dim>::addSampleNodesAndNeighbors();\
\
template \
void GappyPreprocessing<dim>::addNeighbors(int iIslands, int startingNodeWithNeigh);\
\
template \
void GappyPreprocessing<dim>::computeBCFaces(bool liftContribution);\
\
template \
void GappyPreprocessing<dim>::communicateAll();\
\
template \
void GappyPreprocessing<dim>::defineMaps();\
\
template \
void GappyPreprocessing<dim>::communicateBCFaces();\
  \
template \
bool GappyPreprocessing<dim>::checkFaceInMesh(FaceSet& currentFaces, const int iFace, const int iSub, const int *locToGlobNodeMap);\
\
template \
bool GappyPreprocessing<dim>::checkFaceAlreadyAdded(const int cpuNum, const int\
    iSub, const int iFace);\
  \
template \
void GappyPreprocessing<dim>::addFaceNodesElements(FaceSet&\
    currentFaces, const int iFace, const int iSub, const int\
    *locToGlobNodeMap);\
\
template \
void GappyPreprocessing<dim>::addNodesOnFace(FaceSet&\
    currentFaces, const int iFace, const int iSub, const int\
		*locToGlobNodeMap, int *locNodeNums);\
\
template \
void GappyPreprocessing<dim>::addElementOfFace(FaceSet&\
    currentFaces, const int iFace, const int iSub, const int\
    *locToGlobNodeMap, const int *locNodeNums);\
\
template \
bool GappyPreprocessing<dim>::checkFaceContributesToLift(FaceSet& faces, const int iFace, const int iSub, const int *locToGlobNodeMap );\
\
template \
void GappyPreprocessing<dim>::outputTopFile(int iCluster);\
\
template \
void GappyPreprocessing<dim>::outputSampleNodes(int iCluster);\
\
template \
void GappyPreprocessing<dim>::outputSampleNodesGeneral(const std::vector<int> &sampleNodes, const char *outSampleNodeFile);\
\
template \
void GappyPreprocessing<dim>::computeXYZ(int iSub, int iLocNode, double *xyz);\
\
template \
void GappyPreprocessing<dim>::computePseudoInverse();\
\
template \
void GappyPreprocessing<dim>::computePodTPod();\
\
template \
void GappyPreprocessing<dim>::assembleOnlineMatrices();\
\
template \
void GappyPreprocessing<dim>::outputOnlineMatrices(int iCluster);\
\
template \
void GappyPreprocessing<dim>::outputOnlineMatricesGeneral(int iCluster, int numNodes,\
    const boost::unordered_map<int,int> &sampleNodeMap, const std::vector<int>\
    &sampleNodeVec);\
\
template \
void GappyPreprocessing<dim>::outputLocalReferenceStateReduced(int iCluster);\
\
template \
void GappyPreprocessing<dim>::outputInitialConditionReduced();\
\
template \
void GappyPreprocessing<dim>::outputClusterCentersReduced();\
\
template \
void GappyPreprocessing<dim>::outputLocalStateBasisReduced(int iCluster);\
\
template \
void GappyPreprocessing<dim>::outputWallDistanceReduced();\
\
template \
void GappyPreprocessing<dim>::outputReducedVec(const DistVec<double> &distVec, FILE* outFile , int iVector);\
\
template \
void GappyPreprocessing<dim>::outputReducedToFullNodes();\
\
template \
void GappyPreprocessing<dim>::checkConsistency();\
\
template \
void GappyPreprocessing<dim>::formMaskedNonlinearROBs();\
\
template \
void GappyPreprocessing<dim>::formReducedSampleNodeMap();\
\
template \
void GappyPreprocessing<dim>::serialLSMultiRHSGap(int iPodBasis, double **lsCoeff);\
\
template \
void GappyPreprocessing<dim>::probabilisticLSMultiRHSGap(int iPodBasis, double **lsCoeff);\
\
template \
void GappyPreprocessing<dim>::parallelPseudoInverse(int iPodBasis);\
\
template \
void GappyPreprocessing<dim>::serialPseudoInverse(int iPodBasis);\
\
template \
void GappyPreprocessing<dim>::outputMultiSolutionsReduced();\
\
template \
void GappyPreprocessing<dim>::outputReducedSVec(const DistSVec<double,dim> &distSVec, FILE* myOutFile , double tag);\
\
template \
void GappyPreprocessing<dim>::reinitializeMapsForSampleNodes();\
\
template \
void GappyPreprocessing<dim>::maskError();\
\
template \
void GappyPreprocessing<dim>::setUpSampledNodeTargetRegionMask();\
\
template \
void GappyPreprocessing<dim>::outputApproxSnapsReduced(int iCluster);\
\
template \
void GappyPreprocessing<dim>::computeApproxMetricNonlinearNNLS(int iCluster);

INSTANTIATION_HELPER(5);
INSTANTIATION_HELPER(6);
INSTANTIATION_HELPER(7);

#undef INSTANTIATION_HELPER
