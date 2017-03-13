#ifndef _NONLINEAR_ROM_H_
#define _NONLINEAR_ROM_H_

#include <RestrictionMapping.h>
#include <vector>
#include <string>
#include <map>

template <int dim>
class NonlinearRom {

// This class handles all input and output for the nonlinear ROM database.
// Currently, the ROM databse is implemented as a simple file system.  In
// the future, if it is necessary to implement a higher performance
// database, it would only be necessary to change the methods in this class.
// Note that every other nonlinear ROM class is derived from this base class.

  protected:
  Domain& domain;
  Communicator* com; 
  IoData* ioData;
  Timer *timer;

  // IO directory information
  const char* databasePrefix;
  const char* databaseName;
  const char* clusterName;
  const char* sensitivityClusterName;

  NonlinearRomFilesData* romFiles;

  // Map from snapshot file to corresponding reference vector
  std::map<std::string,std::string> mapSnapToRef;

  // State snapshot clusters
  char* stateSnapsName;
  char* mapName;  
  char* indexName;
  char* connName;
  char* centersName;
  char* nearestName;
  char* centerNormsName;
  char* distanceMatrixName;

  // State bases
  char* stateBasisPrefix;
  char* stateBasisName;
  char* stateSingValsName;
  char* simpleUpdateInfoName;
  char* exactUpdateInfoPrefix;      // only user-specified value for exact update file names
  char* basisBasisProductsName;     // for exact updates (exactUpdateInfoPrefix.exactUpdates_F)
  char* basisUrefProductsName;      // for exact updates (exactUpdateInfoPrefix.exactUpdates_e)
  char* basisUicProductsName;       // for exact updates (exactUpdateInfoPrefix.exactUpdates_d)
  char* basisMultiUicProductsName;
  char* urefUicProductsName;        // for exact updates (exactUpdateInfoPrefix.exactUpdates_c)
  char* urefMultiUicProductsName;
  char* urefUrefProductsName;       // for exact updates (exactUpdateInfoPrefix.exactUpdates_g)   
  char* urefComponentwiseSumsName;  // for exact updates with uniform IC (exactUpdateInfoPrefix.exactUpdates_UrefComponentwiseSums)
  char* basisComponentwiseSumsName; // for exact updates with uniform IC (exactUpdateInfoPrefix.exactUpdates_StateBasisComponentwiseSums)
  //char* approxUpdateInfoName; // treating this as a GNAT quantity
  char* stateDistanceComparisonInfoName;
  char* basisNormalizedCenterProductsName;   // for velocity-based local ROMs (incremental bases)
  char* centerComponentWiseSumsName; // for fast distance calcs when using either no updates or exact updates ()
  char* stateDistanceComparisonInfoExactUpdatesName;
  char* stateDistanceComparisonInfoExactUpdatesMultiICName;
  char* projErrorName;
  char* refStateName;
  char* multiUicMultiUicProductsName;

  // Krylov snaps
  char* krylovSnapsName;

  // Krylov bases
  char* krylovBasisPrefix;
  char* krylovBasisName;
  char* krylovSingValsName;
  char* krylovDistanceComparisonInfoName;

  // Sensitivities
  char* sensitivitySnapsName;

  // Sensitivity Basis
  char* sensitivityBasisPrefix;
  char* sensitivityBasisName;
  char* sensitivitySingValsName;
  char* sensitivityDistanceComparisonInfoName;

  // Residual snaps
  char* residualSnapsName;

  // Residual Bases
  char* residualBasisPrefix;
  char* residualBasisName;
  char* residualSingValsName;

  // Action-of-Jacobian snaps
  char* jacActionSnapsName;

  // Action-of-Jacobian Bases
  char* jacActionBasisPrefix;
  char* jacActionBasisName;
  char* jacActionSingValsName;

  // GNAT quantities
  char* sampledNodesName;
  char* sampledNodesFullCoordsName;
  char* sampledCentersName;
  char* sampledStateBasisName;
  char* sampledKrylovBasisName;
  char* sampledSensitivityBasisName;
  char* sampledResidualBasisName;
  char* sampledJacActionBasisName;
  char* sampledMeshName;
  char* sampledSolutionName;
  char* sampledMatchStateName;
  char* sampledMultiSolutionsName;
  char* sampledRefStateName;
  char* sampledWallDistName;
  char* sampledDisplacementName;
  char* sampledShapeDerivativeName;
  char* gappyJacActionName;
  char* gappyResidualName;
  char* approxMetricStateLowRankName;
  char* approxMetricNonlinearLowRankName;
  char* approxMetricStateLowRankFullCoordsName;
  char* approxMetricNonlinearLowRankFullCoordsName;
  char* approxMetricStateLowRankSurfaceCoordsName;
  char* approxMetricNonlinearName;
  char* correlationMatrixName;
  char* sampledApproxMetricNonlinearSnapsName;

  // Surface quantities 
  char* surfaceCentersName;
  char* surfaceStateBasisName;
  char* surfaceRefStateName;
  char* surfaceSolutionName;
  char* surfaceMatchStateName;
  char* surfaceMultiSolutionsName;
  char* surfaceWallDistName;
  char* surfaceDisplacementName;
  char* surfaceShapeDerivativeName;
  char* surfaceMeshName;

  // ROM database data
  VecSet< DistSVec<double, dim> >* snap; // snap(nTotSnaps, domain.getNodeDistInfo())
  VecSet< DistSVec<double, dim> >* clusterCenters; // average of all snapshots in a cluster 
  VecSet< DistSVec<double, dim> >* nearestSnapsToCenters; // closest snapshot to each cluster center 
  int* snapsInCluster; // number of snaps in each cluster
  int* clusterIndex; // stores original cluster association for each snapshot (before any overlap is introduced)
  int** clusterSnapshotMap;  // one vector per cluster, lists snapshots to include in the cluster (including overlapping snapshots) 
  int** clusterNeighbors;  // one vector per cluster, lists neighboring clusters
  int* clusterNeighborsCount; // stores number of neighbors for each cluster
  int* clusterNewtonCount;  // counts number of residuals/jacActions stored in each cluster
  int* clusterKrylovCount;  // counts number of krylov vectors stored in each cluster
  DistSVec<double, dim>* snapRefState; 


  // thin svd update quantities
  // 1: common to all update methods (simple updates, exact updates, approx updates)
  double rTol;
  std::vector<double>* columnSumsV;
  std::vector<double>* sVals;
  DistSVec<double, dim>* Uref;
  int nBuffer; // number of buffer vectors added for updates 
  std::vector<double> cumulativeSnapWeights;  // needed in order to ouput a correct columnSumsV that accounts for normalizing snapshots
  // 2: unique to exact updates
  double uicNorm;
  std::vector<std::vector<std::vector<std::vector<double> > > > basisBasisProducts;  // [iCluster][pCluster][:][:]
  std::vector<std::vector<std::vector<double> > > basisUrefProducts;  // [Cluster_Basis][Cluster_Uref][:]
  std::vector<std::vector<double> > basisUicProducts;  // [iCluster][1:nPod] only precomputed if Uic specified
  std::vector<std::vector<std::vector<double> > > basisMultiUicProducts; // [iCluster][Uj][1:nPod]
  std::vector<double> urefUicProducts; // [iCluster] only precomputed if Uic specified
  std::vector<std::vector<double> > urefMultiUicProducts; // [iCluster][Uj]
  std::vector<std::vector<double> > urefUrefProducts; //[iCluster][jCluster] symmetric (lower triangular)
  std::vector<std::vector<double> > multiUicMultiUicProducts; //[Ui][Uj] symmetric (lower triangular)
  std::vector<std::vector<double> > urefComponentwiseSums; //[iCluster][1:dim]
  std::vector<std::vector<std::vector<double> > > basisComponentwiseSums;  // [iCluster][iVec][1:dim]
  std::vector<double> exactUpdatesAlpha;  // [jVec]
  std::vector<std::vector<double> > exactUpdatesBeta; //[iCluster][jVec]  
  std::vector<std::vector<std::vector<double> > > exactUpdatesN;  // [iCluster][iVec][jVec]
  double exactUpdatesAlphaSwitch;  // scalar
  std::vector<double> exactUpdatesBetaSwitch; //[iCluster]
  std::vector<std::vector<double> > exactUpdatesNSwitch;  // [iCluster][iVec]
  DistSVec<double, dim>* Uic;
  // 3: unique to approximate updates
  VecSet<DistSVec<double, dim> >* lowRankFactor; // low rank factor for approx metric

  // fast distance calculation quantities
  // 1: common to all GNAT update methods (no updates, exact updates, approx updates)
  bool specifiedIC;
  bool interpolatedMultiIC;
  VecSet<DistSVec<double, dim> >* multiUic;

  SVec<double, dim>* uniformIC;  // value of uniform initial condition at node 0 (should be representative)
  std::vector<std::vector<double> > centerNorms; // note: actually norm squared
  std::vector<std::vector<std::vector<std::vector<double> > > > stateBasisCentersDifProduct;  //[iCluster][mCenter][pCenter][:]
  std::vector<std::vector<std::vector<std::vector<double> > > > krylovBasisCentersDifProduct; //[iCluster][mCenter][pCenter][:]
  std::vector<std::vector<std::vector<double> > > sensitivityBasisCentersDifProduct;          //[mCenter][pCenter][:]
  std::vector<std::vector<double> > distanceComparisons;  // this is "z_(m,p)" from Amsallem et al., INJME 2012, but with p<m
  void checkUniformInitialCondition(DistSVec<double, dim> &);
  void checkInitialConditionScenario();
  // 2: unique to exact updates
  std::vector<std::vector<double> > initialConditionCentersDifProduct; 
  std::vector<std::vector<std::vector<double> > > multiUicCentersDifProduct;
  std::vector<std::vector<std::vector<double> > > refStateCentersDifProduct;
  // 3: unique to approximate updates
  double ***hForFastDistComp;
  double ***cForFastDistComp;
  
  // for fast cluster selection in the case of incremental snapshots
  std::vector<std::vector<std::vector<double> > > basisNormalizedCenterProducts; // [Cluster_Basis][Cluster_Center][:]

  // non-database IO function
  int readSnapshotFiles(const char*, bool);
  std::vector<int> stateSnapsFromFile;   // stateSnapsFromFile[iFile] = number of snapshots taken from file iFile
  std::vector<std::vector<double> > stateSnapshotTags; // stateSnapshotInfo[iFile][iSnap] = tag associated with snapshot
                                                       // iSnap from file iFile

  typedef pair<std::string, int> snapID;
  std::vector<snapID> originalSnapshotLocation;

  // database IO functions
  void createDirectories();
  void outputClusteredSnapshots(const char*);
  void readClusteredSnapshots(int, bool, const char*, int first = 0, int last = 0, bool snapshotsAlreadyStored = false);
  void outputClusteredBasis(int, int, const char*);  // readClusteredBasis is public
  void outputClusteredReferenceState(int, DistSVec<double, dim> &);  // automatically stores snapshot reference state
  void readClusteredReferenceState(int, const char*);  // read the reference state that was automatically stored for each cluster.
  void readNearestSnapsToCenters();
  void readReferenceState();  // read a reference state specified by the user
  void readSampleNodes(int iCluster, const char* sampledOrFull, bool deleteExistingRestrictionMapping);  // wrapper for readClusteredSampleNodes, where iCluster<0 indicates union
  void readClusteredSampleNodes(int iCluster, const char* sampledOrFull, bool deleteExistingRestrictionMapping = true);
  void readClusteredGappyMatrix(int, const char*);

  void outputCenterNorms(std::vector<std::vector<double> > &); //note: norm squared
  void readCenterNorms();  // note: norm squared
  void outputClusteredInfoASCII(int, const char*, std::vector<double>* vec1 = NULL, 
                                std::vector<std::vector<double> >* vec2 = NULL,
                                std::vector<std::vector<std::vector<double> > >* vec3 = NULL,
                                std::vector<std::vector<std::vector<std::vector<double> > > >* vec4 = NULL);
  void readClusteredInfoASCII(int, const char*, std::vector<double>* vec1 = NULL,
                              std::vector<std::vector<double> >* vec2 = NULL,
                              std::vector<std::vector<std::vector<double> > >* vec3 = NULL,
                              std::vector<std::vector<std::vector<std::vector<double> > > >* vec4 = NULL);
  void writeMultiVecASCII(char*, std::vector<double>* vec1 = NULL, 
                          std::vector<std::vector<double> >* vec2 = NULL,
                          std::vector<std::vector<std::vector<double> > >* vec3 = NULL,
                          std::vector<std::vector<std::vector<std::vector<double> > > >* vec4 = NULL);
  void readMultiVecASCII(char*, std::vector<double>* vec1 = NULL,
                         std::vector<std::vector<double> >* vec2 = NULL,
                         std::vector<std::vector<std::vector<double> > >* vec3 = NULL,
                         std::vector<std::vector<std::vector<std::vector<double> > > >* vec4 = NULL);


  // for local gappy preprocessing
  void freeMemoryForGappyPrepro();
  void partitionAndSowerForGappy(bool);
  void callSowerSplit(std::string, std::string, char*);

  // for local GNAT online simulations
  int nSampleNodes;
  std::vector<int> sampleNodes;
  int numResJacMat;
  VecSet<DistSVec<double, dim> >* resMat;
  VecSet<DistSVec<double, dim> >* jacMat;
  RestrictionMapping<dim>* restrictionMapping;

  // for storing all clustered quantities in memory (optional)
  bool storedAllOnlineQuantities;
  bool storedAllOfflineQuantities;
  std::vector<int>** allSampleNodes;
  VecSet<DistSVec<double, dim> >** allResMat;
  VecSet<DistSVec<double, dim> >** allJacMat;
  VecSet<DistSVec<double, dim> >** allMetrics;
  VecSet<DistSVec<double, dim> >** allStateBases;
  VecSet<DistSVec<double, dim> >** allKrylovBases;
  VecSet<DistSVec<double, dim> >* sensitivityBasis;
  std::vector<double>** allStateSVals; 
  std::vector<double>** allKrylovSVals;
  std::vector<double>* sensitivitySVals;
  VecSet<DistSVec<double, dim> >* allRefStates;
  std::vector<double>** allColumnSumsV; 
  RestrictionMapping<dim>** allRestrictionMappings;
  std::vector<int> allNBuffer;

  // ASCII output files 
  FILE* clustUsageFile;
  FILE* reducedCoordsFile;

  int nState, nKrylov, nSens;

  public:

  NonlinearRom(Communicator *, IoData &, Domain &);
  ~NonlinearRom();

  int nClusters;
  int nFullMeshNodes;
  int nLowRankFactors;
  VecSet< DistSVec<double, dim> >* basis;
  VecSet< DistSVec<double, dim> >* metric;

  std::vector<double> interpWeightsForMultiIC;

  // When duplicateSnaps is set to true the clustered snapshots are written to the file system, which effectively
  // doubles the required storage.  When false, only a small text file is written.
  bool duplicateSnaps;
  bool euclideanDistances;
  bool incrementalStateSnaps;

  bool jacActionSnapsFileNameSpecified;

  // online selection of closest cluster center (calls either closestCenterFull or closestCenterFast)
  void closestCenter(DistSVec<double, dim> &, int* index1=NULL);

  // calculate closest center to current state using full vectors
  void closestCenterFull(DistSVec<double, dim> &, int* index1=NULL, int* index2=NULL, double* dist1=NULL, double* dist2=NULL);
  void distancesToCentersFull(DistSVec<double, dim> &, std::vector<double> &, int* closest=NULL);
  double distanceFull(DistSVec<double, dim> &, DistSVec<double, dim> &);
  double euclideanFull(DistSVec<double, dim> &, DistSVec<double, dim> &);
  double angleFull(DistSVec<double, dim> &, DistSVec<double, dim> &);

  // calculate closest center to current state without using full vectors (approach depends on ROB update method)
  void closestCenterFast(int* index1=NULL);
  void initializeDistanceComparisons(DistSVec<double, dim> &);
  void resetDistanceComparisonQuantitiesApproxUpdates();
  void advanceDistanceComparisons(int, Vec<double> , Vec<double> );  // calls one of the following five functions
  void distanceComparisonsForIncrements(Vec<double> , int);
  void distanceComparisonsForProjection(Vec<double> , int);
  void incrementDistanceComparisonsForNoUpdates(Vec<double> , int);
  void incrementDistanceComparisonsForExactUpdates(Vec<double> , int);
  void incrementDistanceComparisonsForApproxUpdates(Vec<double> , int);

  // public database IO functions
  void determineFileName(const char*, const char*, const char*, char*&);
  void determinePrefixName(const char*, const char*, char*&);
  void determinePath(char*, int, char*&); // top-level database directory is cluster "-1", sensitivity basis is cluster "-2"
  void readClusteredBasis(int, const char*, bool relProjError = false);
  void readClusteredColumnSumsV(int, const char*);
  void readClusteredUpdateInfo(int, const char*);
  void readNonClusteredUpdateInfo(const char*);
  void readExactUpdateInfo();
  void readProjectionInfo();
  void readClusterCenters(const char*);
  void readAllClusteredOnlineQuantities();
  void readAllClusteredOfflineQuantities();
  void readApproxMetricStateLowRankFactor(const char *); // not clustered
  void readClusteredNonlinearMetric(int iCluster);
  void readDistanceComparisonInfo(const char*); 
  void writeClusteredBinaryVectors(int iCluster, DistSVec<double,dim> *U1 = NULL, DistSVec<double,dim> *U2 = NULL,
                                   DistSVec<double,dim> *U3 = NULL, char* originalSnapshotFile = NULL, int originalSnapshotNumber = 0);
  void initializeClusteredOutputs(); 

  // for online ROMs (both with and without hyper-reduction)
  virtual void projectSwitchStateOntoAffineSubspace(int, int, DistSVec<double, dim> &, Vec<double> &) {};
  virtual bool updateBasis(int, DistSVec<double, dim> &, Vec<double>* coords = NULL) {return false;};
  virtual void appendNonStateDataToBasis(int, const char*, bool relProjError = false) {};
  virtual void readClusteredOnlineQuantities(int) {};
  void truncateBufferedBasis();
  void writeReducedCoords(const int, bool, bool, int, Vec<double>); 
  void initializeProjectionQuantities(DistSVec<double, dim> &);
  void initializeFastExactUpdatesQuantities(DistSVec<double, dim> &);

  // for online ROMs with hyper-reduction
  void determineNumResJacMat(); 
  void deleteRestrictedQuantities();
  int getNumResJacMat() {return numResJacMat;}
  VecSet<DistSVec<double,dim> >* getResMat() {return resMat;}
  VecSet<DistSVec<double,dim> >* getJacMat() {if (numResJacMat==2) { return jacMat; } else { return resMat;} }
  const DistInfo& getRestrictedDistInfo () const {return restrictionMapping->restrictedDistInfo();}
  RestrictionMapping<dim>* restrictMapping() { return restrictionMapping; } 
  void formInterpolatedInitialCondition(DistSVec<double,dim> *U, IoData &iod);

  virtual void appendVectorToBasis(DistSVec<double, dim>&, int numVec = 0) {};

  // general
  void qr(VecSet< DistSVec<double, dim> >* Q, std::vector<std::vector<double> >* RT=NULL, bool testQR=false);
  void rSVD(VecSet< DistSVec<double, dim> >& Utrue, std::vector<double>& singularValues, FullM& Vtrue, bool testSVD=false);
  void probabilisticSVD(VecSet< DistSVec<double, dim> >& Utrue, std::vector<double>& singularValues, FullM& Vtrue, 
                        int k, int nPowerIts, bool testSVD=false);
  void probabilisticLSMultiRHS(VecSet< DistSVec<double, dim> >& LHS, VecSet< DistSVec<double, dim> >& RHS,
                                 std::vector<std::vector<double> >& lsCoeffVec, int k, int nPowerIts, bool testSVD=false);
  void linpackSVD(double* tmpMat, int rows, int columns, double*& yVec, double*& singVals, double*& zVec);
  void printDebug(int iDebug);
};

#include "NonlinearRom.C"
#endif
