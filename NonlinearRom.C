#include <TsInput.h>
#include <cmath>
//#include <time.h>
#include <algorithm>
#include <sys/time.h>
#include <algorithm>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>


using std::stable_sort;

extern "C"      {
   void F77NAME(dsvdc)(double *, int &, int &, int&, double *,
                        double *, double *, int &, double *, int &,
                        double *, const int &, int &);
}


template<int dim> 
NonlinearRom<dim>::NonlinearRom(Communicator *_com, IoData &_ioData, Domain &_domain)  : 
com(_com), ioData(&_ioData), domain(_domain)
{ 
  //NOTE: ioData->example, com->example, domain.example
  timer = domain.getTimer();

  nClusters = ioData->romDatabase.nClusters;  // overwritten later if there are actually fewer clusters
  nFullMeshNodes = 0;  // read from centerNorms file if necessary
  nLowRankFactors = 0;

  // Directory information
  databasePrefix = ioData->romDatabase.directories.prefix;
  databaseName = ioData->romDatabase.directories.databaseName;
  clusterName = ioData->romDatabase.directories.clusterName;
  sensitivityClusterName  = ioData->romDatabase.directories.sensitivityClusterName;

  romFiles = &(ioData->romDatabase.files); 
 
  // When duplicateSnaps is set to true the clustered snapshots are written to the file system, which effectively
  // doubles the required storage.  When false, only a small text file is written.
  duplicateSnaps = (romFiles->duplicateSnapshots==NonlinearRomFilesData::DUPLICATE_SNAPSHOTS_TRUE) ? true : false;

  // flag indicating that we're operating with incremental snapshots of the state
  incrementalStateSnaps = (ioData->romDatabase.avgIncrementalStates==NonlinearRomFileSystemData::AVG_INCREMENTAL_STATES_TRUE) ? true : false;

  // use either euclidean distances or angle between vectors 
  euclideanDistances = (ioData->romDatabase.distanceMetric==NonlinearRomFileSystemData::DIST_EUCLIDEAN) ? true : false;

  // State snapshot clusters
  determineFileName(romFiles->stateSnapsName, "snaps", romFiles->statePrefix, stateSnapsName);
  determineFileName(romFiles->mapName, "map", romFiles->statePrefix, mapName);
  determineFileName(romFiles->indexName, "index", romFiles->statePrefix, indexName);
  determineFileName(romFiles->connName, "conn", romFiles->statePrefix, connName);
  determineFileName(romFiles->centersName, "centers", romFiles->statePrefix, centersName);
  determineFileName(romFiles->nearestName, "nearest", romFiles->statePrefix, nearestName);
  determineFileName(romFiles->centerNormsName, "centerNorms", romFiles->statePrefix, centerNormsName);
  determineFileName(romFiles->distanceMatrixName, "distanceMatrix", romFiles->statePrefix, distanceMatrixName);

  // State bases
  determinePrefixName(romFiles->stateBasisPrefix, romFiles->statePrefix, stateBasisPrefix);
  determineFileName(romFiles->stateBasisName, "rob", stateBasisPrefix, stateBasisName);
  determineFileName(romFiles->stateSingValsName, "svals", stateBasisPrefix, stateSingValsName);
  determineFileName(romFiles->projErrorName, "proj", stateBasisPrefix, projErrorName);
  determineFileName(romFiles->refStateName, "refState", stateBasisPrefix, refStateName); 

  // Update info for state bases (this is a bit tricky)
  determineFileName(romFiles->simpleUpdateInfoName, "allUpdates", stateBasisPrefix, simpleUpdateInfoName);
  determineFileName(romFiles->stateDistanceComparisonInfoName, "distanceInfo", stateBasisPrefix, stateDistanceComparisonInfoName);
  determineFileName(romFiles->basisNormalizedCenterProductsName, "basisCenterProducts", stateBasisPrefix, basisNormalizedCenterProductsName);
  determinePrefixName(romFiles->exactUpdateInfoPrefix, stateBasisPrefix, exactUpdateInfoPrefix);
  determineFileName("", "exactUpdates_F", exactUpdateInfoPrefix, basisBasisProductsName);
  determineFileName("", "exactUpdates_e", exactUpdateInfoPrefix, basisUrefProductsName);
  determineFileName("", "exactUpdates_d", exactUpdateInfoPrefix, basisUicProductsName);
  determineFileName("", "exactUpdates_d_multiIC", exactUpdateInfoPrefix, basisMultiUicProductsName);
  determineFileName("", "exactUpdates_c", exactUpdateInfoPrefix, urefUicProductsName);
  determineFileName("", "exactUpdates_c_multiIC", exactUpdateInfoPrefix, urefMultiUicProductsName);
  determineFileName("", "exactUpdates_g", exactUpdateInfoPrefix, urefUrefProductsName);
  determineFileName("", "exactUpdates_UrefComponentwiseSums", exactUpdateInfoPrefix, urefComponentwiseSumsName);
  determineFileName("", "exactUpdates_StateBasisComponentwiseSums", exactUpdateInfoPrefix, basisComponentwiseSumsName);
  determineFileName("", "exactUpdates_multiICMultiICProducts", exactUpdateInfoPrefix, multiUicMultiUicProductsName);
  determineFileName(romFiles->stateDistanceComparisonInfoExactUpdatesName, "exactUpdatesDistanceInfo", exactUpdateInfoPrefix, stateDistanceComparisonInfoExactUpdatesName); 
  determineFileName(romFiles->stateDistanceComparisonInfoExactUpdatesMultiICName, "exactUpdatesDistanceInfoMultiIC", exactUpdateInfoPrefix, stateDistanceComparisonInfoExactUpdatesMultiICName);

  // Krylov snaps
  determineFileName(romFiles->krylovSnapsName, "snaps", romFiles->krylovPrefix, krylovSnapsName);

  // Krylov bases: NOTE multiple bases is not (yet) implemented up for hyper-reduction
  determinePrefixName(romFiles->krylovBasisPrefix, romFiles->krylovPrefix, krylovBasisPrefix);
  determineFileName(romFiles->krylovBasisName, "rob", krylovBasisPrefix, krylovBasisName);
  determineFileName(romFiles->krylovSingValsName, "svals", krylovBasisPrefix, krylovSingValsName);
  determineFileName(romFiles->krylovDistanceComparisonInfoName, "distanceInfo", krylovBasisPrefix, krylovDistanceComparisonInfoName);

  // Sensitivity snaps
  determineFileName(romFiles->sensitivitySnapsName, "snaps", romFiles->sensitivityPrefix, sensitivitySnapsName);

  // Sensitivity basis: NOTE multiple bases is not (yet) implemented up for hyper-reduction
  determinePrefixName(romFiles->sensitivityBasisPrefix, romFiles->sensitivityPrefix, sensitivityBasisPrefix);
  determineFileName(romFiles->sensitivityBasisName, "rob", sensitivityBasisPrefix, sensitivityBasisName);
  determineFileName(romFiles->sensitivitySingValsName, "svals", sensitivityBasisPrefix, sensitivitySingValsName);
  determineFileName(romFiles->sensitivityDistanceComparisonInfoName, "distanceInfo", sensitivityBasisPrefix, sensitivityDistanceComparisonInfoName);

  // Residual snaps
  determineFileName(romFiles->residualSnapsName, "snaps", romFiles->residualPrefix, residualSnapsName);

  // Residual bases
  determinePrefixName(romFiles->residualBasisPrefix, romFiles->residualPrefix, residualBasisPrefix);
  determineFileName(romFiles->residualBasisName, "rob", residualBasisPrefix, residualBasisName);
  determineFileName(romFiles->residualSingValsName, "svals", residualBasisPrefix, residualSingValsName);

  // Action-of-Jacobian snaps
  determineFileName(romFiles->jacActionSnapsName, "snaps", romFiles->jacActionPrefix, jacActionSnapsName);

  // Action-of-Jacobian bases
  determinePrefixName(romFiles->jacActionBasisPrefix, romFiles->jacActionPrefix, jacActionBasisPrefix);
  determineFileName(romFiles->jacActionBasisName, "rob", jacActionBasisPrefix, jacActionBasisName);
  determineFileName(romFiles->jacActionSingValsName, "svals", jacActionBasisPrefix, jacActionSingValsName);

  // GNAT quantities
  determineFileName(romFiles->sampledNodesName, "sampledNodes", romFiles->gappyPrefix, sampledNodesName);
  determineFileName(romFiles->sampledNodesFullCoordsName, "sampledNodesFullCoords", romFiles->gappyPrefix, sampledNodesFullCoordsName);
  determineFileName(romFiles->sampledCentersName, "sampledCenters", romFiles->gappyPrefix, sampledCentersName);
  determineFileName(romFiles->sampledStateBasisName, "sampledStateROB", romFiles->gappyPrefix, sampledStateBasisName);
  determineFileName(romFiles->sampledSensitivityBasisName, "sampledSensitivityROB", romFiles->gappyPrefix, sampledSensitivityBasisName);
  determineFileName(romFiles->sampledKrylovBasisName, "sampledKrylovROB", romFiles->gappyPrefix, sampledKrylovBasisName);
  determineFileName(romFiles->sampledResidualBasisName, "sampledResROB", romFiles->gappyPrefix, sampledResidualBasisName);
  determineFileName(romFiles->sampledJacActionBasisName, "sampledJacROB", romFiles->gappyPrefix, sampledJacActionBasisName);
  determineFileName(romFiles->sampledMeshName, "top", romFiles->gappyPrefix, sampledMeshName);
  determineFileName(romFiles->sampledSolutionName, "sampledSolution", romFiles->gappyPrefix, sampledSolutionName);
  determineFileName(romFiles->sampledMatchStateName, "sampledMatchState", romFiles->gappyPrefix, sampledMatchStateName);
  determineFileName(romFiles->sampledMultiSolutionsName, "sampledMultiSolutions", romFiles->gappyPrefix, sampledMultiSolutionsName);
  determineFileName(romFiles->sampledRefStateName, "sampledRefState", romFiles->gappyPrefix, sampledRefStateName);
  determineFileName(romFiles->sampledWallDistName, "dwall", romFiles->gappyPrefix, sampledWallDistName);
  determineFileName(romFiles->sampledShapeDerivativeName, "sampledShapeDerivative", romFiles->gappyPrefix, sampledShapeDerivativeName);
  determineFileName(romFiles->sampledDisplacementName, "sampledDisplacement", romFiles->gappyPrefix, sampledDisplacementName);
  determineFileName(romFiles->gappyJacActionName, "gappyJac", romFiles->gappyPrefix, gappyJacActionName);
  determineFileName(romFiles->gappyResidualName, "gappyRes", romFiles->gappyPrefix, gappyResidualName);
  determineFileName(romFiles->approxMetricStateLowRankName, "approxMetricState", romFiles->gappyPrefix, approxMetricStateLowRankName);
  determineFileName(romFiles->approxMetricStateLowRankFullCoordsName, "approxMetricStateFullCoords", romFiles->gappyPrefix, approxMetricStateLowRankFullCoordsName);
  determineFileName(romFiles->approxMetricNonlinearLowRankName, "approxMetricNonlinear", romFiles->gappyPrefix, approxMetricNonlinearLowRankName);
  determineFileName(romFiles->approxMetricNonlinearLowRankFullCoordsName, "approxMetricNonlinearFullCoords", romFiles->gappyPrefix, approxMetricNonlinearLowRankFullCoordsName);
  determineFileName(romFiles->approxMetricNonlinearName, "approxMetricNonlinear", romFiles->gappyPrefix, approxMetricNonlinearName);
  determineFileName(romFiles->correlationMatrixName, "correlationMatrix", romFiles->gappyPrefix, correlationMatrixName);
  determineFileName(romFiles->sampledApproxMetricNonlinearSnapsName, "sampledApproxMetricNLSnaps", romFiles->gappyPrefix, sampledApproxMetricNonlinearSnapsName);

  // Surface quantities 
  determineFileName(romFiles->surfaceCentersName, "surfaceCenters", romFiles->surfacePrefix, surfaceCentersName);
  determineFileName(romFiles->surfaceStateBasisName, "surfaceStateROB", romFiles->surfacePrefix, surfaceStateBasisName);
  determineFileName(romFiles->surfaceRefStateName, "surfaceRefState", romFiles->surfacePrefix, surfaceRefStateName);
  determineFileName(romFiles->surfaceShapeDerivativeName, "surfaceShapeDerivative", romFiles->surfacePrefix, surfaceShapeDerivativeName);
  determineFileName(romFiles->surfaceSolutionName, "surfaceSolution", romFiles->surfacePrefix, surfaceSolutionName);
  determineFileName(romFiles->surfaceMatchStateName, "surfaceMatchState", romFiles->surfacePrefix, surfaceMatchStateName);
  determineFileName(romFiles->surfaceMultiSolutionsName, "surfaceMultiSolutions", romFiles->surfacePrefix, surfaceMultiSolutionsName);
  determineFileName(romFiles->surfaceWallDistName, "dwall", romFiles->surfacePrefix, surfaceWallDistName);
  determineFileName(romFiles->surfaceDisplacementName, "surfaceDisplacement", romFiles->surfacePrefix, surfaceDisplacementName);
  determineFileName(romFiles->surfaceMeshName, "top", romFiles->surfacePrefix, surfaceMeshName);
  determineFileName(romFiles->approxMetricStateLowRankSurfaceCoordsName, "approxMetricStateSurfaceCoords", romFiles->surfacePrefix, approxMetricStateLowRankSurfaceCoordsName);

  basis = NULL;
  snap = NULL; // snap(nTotSnaps, domain.getNodeDistInfo())
  clusterCenters = NULL; // average of all snapshots in a cluster 
  nearestSnapsToCenters = NULL; // closest snapshot to each cluster center 
  snapsInCluster = NULL; // number of snaps in each cluster
  clusterIndex = NULL; // stores original cluster association for each snapshot (before any overlap is introduced)
  clusterSnapshotMap = NULL;  // one vector per cluster, lists snapshots to include in the cluster (including overlapping snapshots) 
  clusterNeighbors = NULL;  // one vector per cluster, lists neighboring clusters
  clusterNeighborsCount = NULL; // stores number of neighbors for each cluster
  snapRefState = NULL;
  columnSumsV = NULL;
  sVals = NULL;
  Uref = NULL;
  clusterNewtonCount = NULL;
  clusterKrylovCount = NULL;
  lowRankFactor = NULL;  
  hForFastDistComp = NULL;
  cForFastDistComp = NULL;
  nSampleNodes = 0;
  sampleNodes.clear();
  numResJacMat = 0;
  resMat = NULL;
  jacMat = NULL;
  metric = NULL;
  restrictionMapping = NULL;
  cumulativeSnapWeights.clear();

  nBuffer = 0;

  clustUsageFile = NULL;
  reducedCoordsFile = NULL;

  // initialize ASCII outputs
  if (!(strcmp(ioData->output.rom.clusterUsage,"")==0) && nClusters>0) {
    char *fullClustUsageName = new char[strlen(ioData->output.rom.prefix) + 1 + strlen(ioData->output.rom.clusterUsage) + 1];
    sprintf(fullClustUsageName, "%s%s", ioData->output.rom.prefix, ioData->output.rom.clusterUsage);
    if (com->cpuNum() == 0)  clustUsageFile = fopen(fullClustUsageName, "wt");
    delete [] fullClustUsageName;
  }

  if (strcmp(ioData->output.rom.reducedCoords,"")==0) {
    if (ioData->romOnline.systemApproximation == NonlinearRomOnlineData::GNAT 
        || ioData->romOnline.systemApproximation == NonlinearRomOnlineData::COLLOCATION 
        || ioData->romOnline.systemApproximation == NonlinearRomOnlineData::APPROX_METRIC_NL)
      com->fprintf(stderr, "\n*** Warning: Reduced coordinates output file not specified\n\n");
  } else {
    char *fullReducedCoordsName = new char[strlen(ioData->output.rom.prefix) + 1 + strlen(ioData->output.rom.reducedCoords) + 1];
    sprintf(fullReducedCoordsName, "%s%s", ioData->output.rom.prefix, ioData->output.rom.reducedCoords);
    if (com->cpuNum() == 0)  reducedCoordsFile = fopen(fullReducedCoordsName, "wt");
    delete [] fullReducedCoordsName;
  }

  nState = 0;
  nKrylov = 0;
  nSens = 0;

  // quantities for storing all necessary online/offline info in RAM to limit IO.
  storedAllOnlineQuantities = false;
  storedAllOfflineQuantities = false;
  allSampleNodes = NULL;
  allResMat = NULL;
  allJacMat = NULL;
  allMetrics = NULL;
  allStateBases = NULL;
  allKrylovBases = NULL;
  sensitivityBasis = NULL;
  allStateSVals = NULL;
  allKrylovSVals = NULL;
  sensitivitySVals = NULL;
  allRefStates = NULL;
  allColumnSumsV = NULL; 
  allRestrictionMappings = NULL;
  allNBuffer.clear();

  // for fast distance calculation quantities / exact update quantitiess
  specifiedIC = false;
  interpolatedMultiIC = false;
  uniformIC = NULL;
  multiUic = NULL;

  rTol = ioData->romOnline.basisUpdateTolerance;
  jacActionSnapsFileNameSpecified = (strcmp(jacActionSnapsName,"")!=0);

}

//----------------------------------------------------------------------------------

template<int dim> 
NonlinearRom<dim>::~NonlinearRom() 
{

  delete [] stateSnapsName;
  delete [] mapName;  
  delete [] indexName;
  delete [] connName;
  delete [] centersName;
  delete [] nearestName;
  delete [] centerNormsName;
  delete [] stateBasisName;
  delete [] stateSingValsName;
  delete [] simpleUpdateInfoName;
  delete [] exactUpdateInfoPrefix;
  delete [] basisBasisProductsName;
  delete [] basisUrefProductsName;
  delete [] basisNormalizedCenterProductsName;
  delete [] basisUicProductsName;  
  delete [] basisMultiUicProductsName;
  delete [] urefUicProductsName;
  delete [] urefMultiUicProductsName;
  delete [] urefUrefProductsName;    
  delete [] multiUicMultiUicProductsName;    
  delete [] urefComponentwiseSumsName;
  delete [] basisComponentwiseSumsName;
  delete [] stateDistanceComparisonInfoName;
  delete [] stateDistanceComparisonInfoExactUpdatesName;
  delete [] stateDistanceComparisonInfoExactUpdatesMultiICName;
  delete [] refStateName;
  delete [] projErrorName;
  delete [] krylovSnapsName;
  delete [] krylovBasisName;
  delete [] krylovSingValsName;
  delete [] krylovDistanceComparisonInfoName;
  delete [] residualSnapsName;
  delete [] jacActionSnapsName;
  delete [] residualBasisName;
  delete [] residualSingValsName;
  delete [] jacActionBasisName;
  delete [] jacActionSingValsName;
  delete [] sampledNodesName;
  delete [] sampledCentersName;
  delete [] sampledStateBasisName;
  delete [] sampledKrylovBasisName;
  delete [] sampledSensitivityBasisName;
  delete [] sampledResidualBasisName;
  delete [] sampledJacActionBasisName;
  delete [] sampledMeshName;
  delete [] sampledSolutionName;
  delete [] sampledShapeDerivativeName;
  delete [] sampledMultiSolutionsName;
  delete [] sampledRefStateName;
  delete [] sampledWallDistName;
  delete [] sampledDisplacementName;
  delete [] gappyJacActionName;
  delete [] gappyResidualName;
  delete [] approxMetricNonlinearLowRankName;
  delete [] approxMetricStateLowRankName;
  delete [] approxMetricStateLowRankFullCoordsName;
  delete [] approxMetricNonlinearLowRankFullCoordsName;
  delete [] approxMetricStateLowRankSurfaceCoordsName;
  delete [] approxMetricNonlinearName;
  delete [] correlationMatrixName;
  delete [] sampledApproxMetricNonlinearSnapsName;
  delete [] surfaceCentersName;
  delete [] surfaceStateBasisName;
  delete [] surfaceRefStateName;
  delete [] surfaceSolutionName;
  delete [] surfaceShapeDerivativeName;
  delete [] surfaceMultiSolutionsName;
  delete [] surfaceWallDistName;
  delete [] surfaceDisplacementName;
  delete [] surfaceMeshName;

  delete [] stateBasisPrefix;
  delete [] krylovBasisPrefix;
  delete [] sensitivityBasisPrefix;
  delete [] residualBasisPrefix;
  delete [] jacActionBasisPrefix;

  if (lowRankFactor) delete lowRankFactor;
  if (hForFastDistComp) {  
    for (int iCluster = 0; iCluster < nClusters; ++iCluster) {
      for (int jCluster = 0; jCluster < nClusters; ++jCluster) { 
        delete [] hForFastDistComp[iCluster][jCluster];
        delete [] cForFastDistComp[iCluster][jCluster];
      }
      delete [] hForFastDistComp[iCluster];
      delete [] cForFastDistComp[iCluster];
    } 
    delete [] hForFastDistComp;
    delete [] cForFastDistComp;
  }

  if (basis) delete basis;
  if (snap) delete snap; 
  if (clusterCenters) delete clusterCenters;
  if (nearestSnapsToCenters) delete nearestSnapsToCenters;
  if (snapsInCluster) delete [] snapsInCluster;
  if (clusterIndex) delete [] clusterIndex;
  if (clusterNeighborsCount) delete [] clusterNeighborsCount;
  if (clusterNewtonCount) delete [] clusterNewtonCount;
  if (clusterKrylovCount) delete [] clusterKrylovCount;
  if (snapRefState) delete snapRefState;
  if (columnSumsV) delete columnSumsV; 
  if (sVals) delete sVals;
  if (Uref) delete Uref;
  cumulativeSnapWeights.clear();

  //TODO
  //clusterSnapshotMap
  //clusterNeighbors

  if (resMat) delete resMat;
  if (jacMat) delete jacMat;
  if (metric) delete metric;

  if (com->cpuNum() == 0) {
    if (clustUsageFile) fclose(clustUsageFile);
    if (reducedCoordsFile) fclose(reducedCoordsFile);
  }

  if (storedAllOnlineQuantities || storedAllOfflineQuantities) {
    for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      if (allSampleNodes) delete allSampleNodes[iCluster];
      if (allResMat) delete allResMat[iCluster];
      if (allJacMat) delete allJacMat[iCluster];
      if (allMetrics) delete allMetrics[iCluster];
      if (allStateBases) delete allStateBases[iCluster];
      if (allStateSVals) delete allStateSVals[iCluster];
      if (allKrylovBases) delete allKrylovBases[iCluster];
      if (allKrylovSVals) delete allKrylovSVals[iCluster];
      if (allColumnSumsV) delete allColumnSumsV[iCluster];
      if (allRestrictionMappings) delete allRestrictionMappings[iCluster];
    }
    if (allSampleNodes) delete [] allSampleNodes;
    if (allResMat) delete [] allResMat;
    if (allJacMat) delete [] allJacMat;
    if (allMetrics) delete [] allMetrics;
    if (allStateBases) delete [] allStateBases;
    if (allStateSVals) delete [] allStateSVals;
    if (allKrylovBases) delete [] allKrylovBases;
    if (allKrylovSVals) delete [] allKrylovSVals;
    if (allColumnSumsV) delete [] allColumnSumsV;
    if (allRestrictionMappings) delete [] allRestrictionMappings; 
    if (sensitivityBasis) delete sensitivityBasis;
    if (sensitivitySVals) delete sensitivitySVals;
    if (allRefStates) delete allRefStates;
    allNBuffer.clear();
  } else {
    if (restrictionMapping) delete restrictionMapping;
  }

  if (uniformIC) delete uniformIC;
  if (multiUic) delete multiUic;

}

//---------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::freeMemoryForGappyPrepro() {

  if (basis) {
    delete basis;
    basis=NULL;
  }
  if (snap) {
    delete snap; 
    snap=NULL;
  }
  if (nearestSnapsToCenters) {
    delete nearestSnapsToCenters;
    nearestSnapsToCenters=NULL; 
  }
  if (snapsInCluster) {
    delete [] snapsInCluster;
    snapsInCluster=NULL;
  }
  if (clusterIndex) {
    delete [] clusterIndex;
    clusterIndex=NULL;
  }
  if (clusterNeighborsCount) {
    delete [] clusterNeighborsCount;
    clusterNeighborsCount=NULL;
  }
  if (clusterNewtonCount) {
    delete [] clusterNewtonCount;
    clusterNewtonCount=NULL;
  }
  if (clusterKrylovCount) {
    delete [] clusterKrylovCount;
    clusterKrylovCount=NULL;
  }
  if (snapRefState) {
    delete snapRefState;
    snapRefState=NULL;
  }
  if (columnSumsV) {
    delete columnSumsV;
    columnSumsV=NULL; 
  }
  if (sVals) {
    delete sVals;
    sVals=NULL;
  }
  if (Uref) {
    delete Uref;
    Uref=NULL;
  }

  if (storedAllOfflineQuantities) {
    for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      if (allStateBases) delete allStateBases[iCluster];
      if (allStateSVals) delete allStateSVals[iCluster];
      if (allKrylovBases) delete allKrylovBases[iCluster];
      if (allKrylovSVals) delete allKrylovSVals[iCluster];
      if (allColumnSumsV) delete allColumnSumsV[iCluster];
    }
    if (allStateBases) {
      delete [] allStateBases;
      allStateBases=NULL;
    }
    if (allStateSVals) {
      delete [] allStateSVals;
      allStateSVals=NULL;
    }
    if (allKrylovBases) {
      delete [] allKrylovBases;
      allKrylovBases=NULL;
    }
    if (allKrylovSVals) {
      delete [] allKrylovSVals;
      allKrylovSVals=NULL;
    }
    if (allColumnSumsV) {
      delete [] allColumnSumsV; 
      allColumnSumsV=NULL;
    }
    if (sensitivityBasis) {
      delete sensitivityBasis;
      sensitivityBasis=NULL;
    }
    if (sensitivitySVals) {
      delete sensitivitySVals;
      sensitivitySVals=NULL;
    }
    if (allRefStates) {
      delete allRefStates;
      allRefStates=NULL;
    }
    storedAllOfflineQuantities=false;
  }

  if (lowRankFactor) {
    delete lowRankFactor;
    lowRankFactor=NULL;
  }

  allNBuffer.clear();
  cumulativeSnapWeights.clear();
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::closestCenter(DistSVec<double, dim> &vec, int *index1) {

  // for determining the closest cluster center during an online ROM simulation

  if (nClusters == 1) {
    *index1 = 0;
    return;
  }

  if (ioData->romOnline.distanceComparisons) {
    closestCenterFast(index1);
  } else {
    closestCenterFull(vec, index1);
  }

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::distancesToCentersFull(DistSVec<double, dim> &vec, std::vector<double> &distances, int* closest) {
// Computes the distance from vec to every cluster center.

  distances.resize(nClusters);  

  for (int iCluster=0; iCluster<nClusters; ++iCluster) {
    distances[iCluster] = distanceFull( vec, (*clusterCenters)[iCluster]);
  }

  if (closest) {
    *closest = 0;
    for (int iCluster=1; iCluster<nClusters; ++iCluster)
      *closest = (distances[*closest]<distances[iCluster]) ? *closest : iCluster;
  }

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::closestCenterFull(DistSVec<double, dim> &vec, int *index1, int *index2, double *dist1, double *dist2) {
// Computes the closest cluster center to vec using full scale vectors.

  double tmp;
  
  int i1 = 0; // index of closest cluster
  int i2 = 0; // index of second closest cluster
  double d1 = 0;  // distance to closest cluster center
  double d2 = 0;  // distance to second closest cluster center

  for (int iCluster=0; iCluster<nClusters; ++iCluster) {
    tmp = distanceFull( vec, (*clusterCenters)[iCluster]);
    if ((tmp<d1) || (iCluster==0)) {
      d2 = d1;
      i2 = i1;
      d1 = tmp;
      i1 = iCluster;
    } else if ((tmp<d2) || (iCluster==1)) {
      d2 = tmp;
      i2 = iCluster;
    }
  }
 
  if (index1) *index1 = i1;
  if (index2) *index2 = i2;
  if (dist1) *dist1 = d1;
  if (dist2) *dist2 = d2;

}


//----------------------------------------------------------------------------------

template<int dim>
double NonlinearRom<dim>::distanceFull(DistSVec<double, dim> &U1, DistSVec<double, dim> &U2) {
// Distance calculation using full-size state vectors

  double dist;

  if (euclideanDistances) {
    dist = euclideanFull(U1,U2);
  } else {
    dist = angleFull(U1,U2);
  }

  return dist;
}

//----------------------------------------------------------------------------------

template<int dim>
double NonlinearRom<dim>::euclideanFull(DistSVec<double, dim> &U1, DistSVec<double, dim> &U2) {
// Euclidean distance calculation using full-size state vectors

  DistSVec<double, dim> diff(domain.getNodeDistInfo());
  diff = (U1 - U2);
  double dist = diff.norm();

  return dist;
}

//----------------------------------------------------------------------------------

template<int dim>
double NonlinearRom<dim>::angleFull(DistSVec<double, dim> &U1, DistSVec<double, dim> &U2) {
// Angle calculation between two full-size state vectors.  
// Returns an angle between 0 and pi/2 unless one vector is too small (then return -1)

  double norm1 = U1.norm();
  double norm2 = U2.norm();

  if ((norm1 < 1e-8) || (norm2 < 1e-8)) {
    // at least one vector is tiny... probably a bad idea to compute the angle
    com->fprintf(stderr, "*** Warning: one vector passed into function angleFull(Vec1,Vec2) is basically zero. Norm1=%e, Norm2=%e\n",
                 norm1,norm2);
    return -1;
  }
  double arg = (U1*U2)/(norm1*norm2);
  if (arg>1.0) arg=1.0;
  if (arg<-1.0) arg=-1.0;
  double angle = acos(arg);

  return angle;
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::closestCenterFast(int *index1) {
// returns the index of the closest cluster center (using precomputed distance quantities)

  int closestCluster = 0;

  for (int iCluster=1; iCluster<nClusters; ++iCluster) {
    closestCluster = (distanceComparisons[iCluster][closestCluster] < 0) ? iCluster : closestCluster;
  }

  *index1 = closestCluster;

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::checkInitialConditionScenario() {

  if (specifiedIC || interpolatedMultiIC || uniformIC) return;

  if (strcmp(ioData->input.multiSolutionsParams,"")!=0) {
    com->fprintf(stdout, " ... Interpolated initial condition \n");
    interpolatedMultiIC = true;
  } else if (strcmp(ioData->input.solutions,"")!=0) {
    com->fprintf(stdout, " ... Initial condition file is specified \n");
    specifiedIC = true;
  } else {
    com->fprintf(stdout, " ... No initial solution specified -- checking for uniform initial condition\n");
  }

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::checkUniformInitialCondition(DistSVec<double, dim> &ic) {

    checkInitialConditionScenario();

    if (specifiedIC || interpolatedMultiIC || uniformIC) return;
    
    // double check that the IC is indeed uniform
    double minIC[dim], maxIC[dim];
    ic.min(minIC);
    ic.max(maxIC);

    for (int iDim=0; iDim<dim; ++iDim) {
      if (minIC[iDim]!=maxIC[iDim]) { // should be exact, no tolerance needed
       fprintf(stderr, " *** Error: expected a uniform initial condition (uniformIC test failed)!");
       exit(-1);
      }
    }
 
    uniformIC = new SVec<double, dim>( ic(0) ); // value of initial condition at node 0 (should be representative)

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::resetDistanceComparisonQuantitiesApproxUpdates() {
// resets some fast distance comparison quantities that are unique to approximate updates

  if (nClusters == 1) return;

  int robSize = basis->numVectors();

  for (int iCluster = 0; iCluster < nClusters; ++iCluster) {
    for (int jCluster = 0; jCluster < nClusters; ++jCluster) {
      delete [] (hForFastDistComp)[iCluster][jCluster];
      (hForFastDistComp)[iCluster][jCluster] = new double[robSize];
    }
  } 
      
  double *temp = new double[nLowRankFactors];
  for (int iVec = 0; iVec < robSize; ++iVec) {
    for (int iRank = 0; iRank < nLowRankFactors; ++iRank)
      temp[iRank] = (*lowRankFactor)[iRank] * (*(basis))[iVec];
    for (int iCluster = 0; iCluster < nClusters; ++iCluster) {
      for (int jCluster = 0; jCluster < nClusters; ++jCluster) {
        (this->hForFastDistComp)[iCluster][jCluster][iVec] = 0.0;
        for (int iRank = 0; iRank < nLowRankFactors; ++iRank)
          (hForFastDistComp)[iCluster][jCluster][iVec] += ( (cForFastDistComp[iCluster][jCluster][iRank]) * temp[iRank]);
      } 
    } 
  }
  delete [] temp;
  temp = NULL;

}



//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::initializeProjectionQuantities(DistSVec<double, dim> &ic) {

  // we always have basisBasisProducts basisUrefProducts precomputed (regardless of IC)

  checkInitialConditionScenario();

  if (specifiedIC) {
    // do nothing; basisUicProducts should already be stored
  } else if (interpolatedMultiIC) { 
    // do nothing; basisUicProducts should already have been formed
    // using basisMultiUicProducts and the interpWeightsForMultiIC
  } else {
    // we have urefComponentwiseSums and basisComponentwiseSums, but we need urefUicProducts and basisUicProducts 
    checkUniformInitialCondition(ic); // checks if IC is uniform, sets variable uniformIC
    this->basisUicProducts.resize(nClusters);

    for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      int nVecs = this->basisComponentwiseSums[iCluster].size();
      this->basisUicProducts[iCluster].clear();
      this->basisUicProducts[iCluster].resize(nVecs, 0.0);
      for (int iVec=0; iVec<nVecs; ++iVec) {
        for (int iDim=0; iDim<dim; ++iDim) {
          this->basisUicProducts[iCluster][iVec] += (*uniformIC)[0][iDim] * this->basisComponentwiseSums[iCluster][iVec][iDim];  
        }
      }
    }
  }

}



//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::initializeFastExactUpdatesQuantities(DistSVec<double, dim> &ic) {

  // we always have basisBasisProducts, basisUrefProducts, and basisUrefProducts precomputed (regardless of IC)

  checkInitialConditionScenario();

  if (specifiedIC) {
    // do nothing; basisUicProducts and urefUicProducts should already be stored
  } else if (interpolatedMultiIC) {
    // do nothing; basisUicProducts and urefUicProducts should already have been formed
    // using basisMultiUicProducts, urefMultiUicProducts, and the interpWeightsForMultiIC
  } else {
    // we have urefComponentwiseSums and basisComponentwiseSums, but we need urefUicProducts and basisUicProducts 
    checkUniformInitialCondition(ic); // checks if IC is uniform, sets variable uniformIC
    this->urefUicProducts.clear();
    this->urefUicProducts.resize(nClusters, 0.0);
    this->basisUicProducts.resize(nClusters);

    for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      for (int iDim=0; iDim<dim; ++iDim) {
        this->urefUicProducts[iCluster] += (*uniformIC)[0][iDim] * this->urefComponentwiseSums[iCluster][iDim];
      }
 
      int nVecs = this->basisComponentwiseSums[iCluster].size();
      this->basisUicProducts[iCluster].clear();
      this->basisUicProducts[iCluster].resize(nVecs, 0.0);
      for (int iVec=0; iVec<nVecs; ++iVec) {
        for (int iDim=0; iDim<dim; ++iDim) {
          this->basisUicProducts[iCluster][iVec] += (*uniformIC)[0][iDim] * this->basisComponentwiseSums[iCluster][iVec][iDim];  
        }
      }
    }
  }

  Uic = new DistSVec<double, dim>(domain.getNodeDistInfo());
  *Uic = ic;  // TODO: this is stored in two places for steady GNAT simulations with exact updates

  if (ioData->romOnline.systemApproximation == NonlinearRomOnlineData::SYSTEM_APPROXIMATION_NONE) {
    this->uicNorm = Uic->norm();
  } else if (ioData->romOnline.systemApproximation == NonlinearRomOnlineData::GNAT 
             || ioData->romOnline.systemApproximation == NonlinearRomOnlineData::COLLOCATION
             || ioData->romOnline.systemApproximation == NonlinearRomOnlineData::APPROX_METRIC_NL) {
    if (specifiedIC) {
      double tag = 0.0;
      int numVecs = 0;
      int step = 0;
      char *solutionPath = new char[strlen(ioData->input.prefix) + strlen(ioData->input.solutions) + 1];
      sprintf(solutionPath, "%s%s", ioData->input.prefix, ioData->input.solutions);
      bool status = domain.readTagFromFile<double, dim>(solutionPath, step, &tag, &numVecs);  // if file DNE, returns false, tag=0, and numSteps=0
      delete [] solutionPath;
      if (!status) {
        fprintf(stderr, "\nCould not open file \"%s\"\n", solutionPath);
        exit(-1);
      }
      this->uicNorm = tag;
    } else if (interpolatedMultiIC) {
      // form uicNorm (using interpWeightsForMultiIC and multiUicMultiUicProducts)
      readClusteredInfoASCII(-1, "multiUicMultiUicProducts", NULL, &multiUicMultiUicProducts);
      int nSol = multiUicMultiUicProducts.size();
      std::vector<double> result;
      result.resize(nSol,0.0);
      for (int iSol=0; iSol<nSol; ++iSol) {
        for (int jSol=0; jSol<nSol; ++jSol) {
          if (iSol>=jSol) {
            result[iSol] += multiUicMultiUicProducts[iSol][jSol]*interpWeightsForMultiIC[jSol];
          } else {
            result[iSol] += multiUicMultiUicProducts[jSol][iSol]*interpWeightsForMultiIC[jSol];
          }
        }
      } 
      multiUicMultiUicProducts.clear();

      this->uicNorm = 0.0;
      for (int iSol=0; iSol<nSol; ++iSol) {
        this->uicNorm += result[iSol]*interpWeightsForMultiIC[iSol];
      }
      result.clear();
      this->uicNorm = pow(this->uicNorm,0.5);

    } else {
      if (nFullMeshNodes==0) {
         readCenterNorms();
         if (nFullMeshNodes==0) {
           fprintf(stderr, "\n*** Error: Number of full mesh nodes = 0?  (Remember to preprocess for exact updates!)\n");
           exit(-1);
         }
      }
      this->uicNorm = 0.0;
      for (int iDim = 0; iDim<dim; ++iDim)
        this->uicNorm += pow((*uniformIC)[0][iDim],2);
      this->uicNorm *= double(nFullMeshNodes); 
      this->uicNorm = pow(this->uicNorm,0.5);
    }
  } else {
    fprintf(stderr, "*** Error: Unexpected system approximation method\n");
    exit(-1);
  }

  exactUpdatesAlpha.clear();          // [jVec]
  exactUpdatesBeta.resize(nClusters);   // [iCluster][jVec]  
  exactUpdatesN.resize(nClusters);      // [iCluster][iVec][jVec]
  exactUpdatesAlphaSwitch = 1.0;        // scalar
  exactUpdatesBetaSwitch.clear();
  exactUpdatesBetaSwitch.resize(nClusters, 0.0); //[iCluster]
  exactUpdatesNSwitch.resize(nClusters);// [iCluster][iVec]
 
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::initializeDistanceComparisons(DistSVec<double, dim> &ic) {

  if (nClusters == 1) return;

  // the diagonal elements are all zero, but element [0][0] is inserted to make the indexing more intuitive ([mCenter][pCenter])
  distanceComparisons.clear();
  distanceComparisons.resize(nClusters);  

  if (incrementalStateSnaps) {

    // not needed at first iteration (need to use tryAllClusters instead, since there's no increment to compare)
    for (int mCenter=1; mCenter<nClusters; ++mCenter) {
      distanceComparisons[mCenter].reserve(mCenter);
      for (int pCenter=0; pCenter<mCenter; ++pCenter) {
        distanceComparisons[mCenter].push_back(0.0);
      }
    }

  } else {
  
    checkInitialConditionScenario();
  
    if (specifiedIC) { // preprocessing was performed for a specified IC 
      // centerNorms is nClusters-by-1 with format: ||Center_0||^2;  ||Center_1||^2; ...
      // initialConditionCentersDifProduct has been precomputed

    } else if (interpolatedMultiIC) {
      // centerNorms is nClusters-by-1 with format: ||Center_0||^2;  ||Center_1||^2; ...
      // initialConditionCentersDifProduct has been formed and stored

    } else { // preprocesing was performed assuming a uniform initial condition for the online ROM
  
      // centerNorms is nClusters-by-(dim+1) with format:
      //    ||Center_0||^2, sum( Center_0(iDim=0) ), ... , sum( Center_0(iDim=dim-1) ); ||Center_1||^2, ...

      // need to form initialConditionCentersDifProduct 
 
      checkUniformInitialCondition(ic); // checks if IC is uniform, sets variable uniformIC
  
      if (ioData->romOnline.basisUpdates == NonlinearRomOnlineData::UPDATES_FAST_EXACT) {
        initialConditionCentersDifProduct.resize(nClusters);
        for (int mCenter=1; mCenter<nClusters; ++mCenter) {
          initialConditionCentersDifProduct[mCenter].resize(mCenter);
          for (int pCenter=0; pCenter<mCenter; ++pCenter) { 
            initialConditionCentersDifProduct[mCenter][pCenter] = 0.0;
            for (int iDim=0; iDim<dim; ++iDim) {
              initialConditionCentersDifProduct[mCenter][pCenter] += 2.0 * (*uniformIC)[0][iDim] * 
                                                                  (centerNorms[pCenter][iDim+1] - centerNorms[mCenter][iDim+1]);
            }
          }
        }
      }
    }

    for (int mCenter=1; mCenter<nClusters; ++mCenter) {
      distanceComparisons[mCenter].resize(mCenter,0.0);
      for (int pCenter=0; pCenter<mCenter; ++pCenter) {
        distanceComparisons[mCenter][pCenter]=(centerNorms[mCenter][0] - centerNorms[pCenter][0] + initialConditionCentersDifProduct[mCenter][pCenter]);
      }
    }

  }

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::advanceDistanceComparisons(int currentCluster, Vec<double> dUromTimeIt, Vec<double> UromCurrentROB) {

  if (nClusters == 1) return;

  if (incrementalStateSnaps) {
    distanceComparisonsForIncrements(dUromTimeIt, currentCluster);
  } else { 
    switch (ioData->romOnline.basisUpdates) {
      case (NonlinearRomOnlineData::UPDATES_OFF):
        if (ioData->romOnline.projectSwitchStateOntoAffineSubspace!=NonlinearRomOnlineData::PROJECT_OFF) {
          distanceComparisonsForProjection(UromCurrentROB, currentCluster); // no need to increment, easy to calculate directly
        } else {
          incrementDistanceComparisonsForNoUpdates(dUromTimeIt, currentCluster);
        }
        break;
      case (NonlinearRomOnlineData::UPDATES_SIMPLE):
        fprintf(stderr, "*** Error: fast distance comparisons are incompatible with simple ROB updates (use Exact)\n");
        exit(-1);
        break;
      case (NonlinearRomOnlineData::UPDATES_FAST_EXACT):
        incrementDistanceComparisonsForExactUpdates(dUromTimeIt, currentCluster);
        break;
      case (NonlinearRomOnlineData::UPDATES_FAST_APPROX):
        incrementDistanceComparisonsForApproxUpdates(dUromTimeIt, currentCluster);
        break;
      default:
        fprintf(stderr, "*** Error: Unexpected ROB updates method\n");
        exit(-1);
    }  
  }
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::distanceComparisonsForIncrements(Vec<double> dUromTimeIt, int currentCluster) {

  if (nClusters == 1) return;

  for (int mCenter=1; mCenter<nClusters; ++mCenter) {

    double mComponent = 0.0;
    for (int iState=0; iState<dUromTimeIt.size(); ++iState)
      mComponent += basisNormalizedCenterProducts[currentCluster][mCenter][iState] * dUromTimeIt[iState];

    for (int pCenter=0; pCenter<mCenter; ++pCenter) {

      double pComponent = 0.0;
      for (int iState=0; iState<dUromTimeIt.size(); ++iState)
        pComponent += basisNormalizedCenterProducts[currentCluster][pCenter][iState] * dUromTimeIt[iState];

      distanceComparisons[mCenter][pCenter] = pComponent - mComponent; //abs(pComponent) - abs(mComponent); 
                                                
      for (int iKrylov=0; iKrylov<nKrylov; ++iKrylov) {
        // account for Krylov bases
      }
      for (int iSens=0; iSens<nSens; ++iSens) {
        // account for sensitivity bases
      }
    }
  }

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::distanceComparisonsForProjection(Vec<double> UromCurrentROB, int currentCluster) {
// note:
// refStateCentersDifProduct[i][m][p] = 2 * Uref_i^T *(center_p - center_m)
// stateBasisCentersDifProduct[i][m][p] = 2 * StateROB_i^T *(center_p - center_m)
// centerNorms[i] = center_i^T * center_i
 
  if (nClusters == 1) return;

  for (int mCenter=1; mCenter<nClusters; ++mCenter) {
    for (int pCenter=0; pCenter<mCenter; ++pCenter) {
      distanceComparisons[mCenter][pCenter] = refStateCentersDifProduct[currentCluster][mCenter][pCenter]
                                              + centerNorms[mCenter][0] - centerNorms[pCenter][0];
      for (int iState=0; iState<UromCurrentROB.size(); ++iState) {
        distanceComparisons[mCenter][pCenter] += stateBasisCentersDifProduct[currentCluster][mCenter][pCenter][iState] 
                                                 * UromCurrentROB[iState];
      }
      for (int iKrylov=0; iKrylov<nKrylov; ++iKrylov) {
        // account for Krylov bases
      }
      for (int iSens=0; iSens<nSens; ++iSens) {
        // account for sensitivity bases
      }
    }
  }

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::incrementDistanceComparisonsForNoUpdates(Vec<double> dUromTimeIt, int currentCluster) {
 
  if (nClusters == 1) return;

  for (int mCenter=1; mCenter<nClusters; ++mCenter) {
    for (int pCenter=0; pCenter<mCenter; ++pCenter) {
      for (int iState=0; iState<dUromTimeIt.size(); ++iState) {
        distanceComparisons[mCenter][pCenter] += stateBasisCentersDifProduct[currentCluster][mCenter][pCenter][iState] 
                                                 * dUromTimeIt[iState];
      }
      for (int iKrylov=0; iKrylov<nKrylov; ++iKrylov) {
        // account for Krylov bases
      }
      for (int iSens=0; iSens<nSens; ++iSens) {
        // account for sensitivity bases
      }
    }
  }

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::incrementDistanceComparisonsForExactUpdates(Vec<double> dUromTimeIt, int currentCluster) {

  // e_(m,p) = 2 * Uic^T *(center_p - center_m) for 0<=p<m<nCluster

  // f_(m,p) = 2 * Uref_i^T *(center_p - center_m) for 0<=p<m<nClusters
  // note that f_(m,p) = -f_(p,m)
  
  // g_(m,p) = 2 * basis^T *(center_p - center_m) for 0<=p<m<nClusters
  // note that g_(m,p) = -w_(p,m)

  // d_(m,p) is computed during initializeDistanceComparisons and is not needed here
  // e_(m,p) = initialConditionCentersDifProduct[mCenter][pCenter]
  // f_(m,p) = refStateCentersDifProduct[iRefState][mCenter][pCenter]
  // g_(m,p) = stateBasisCentersDifProduct[iCluster][mCenter][pCenter][iState]

  if (nClusters == 1) return;
 
  std::vector<double> tmp;

  for (int mCenter=1; mCenter<nClusters; ++mCenter) {
    for (int pCenter=0; pCenter<mCenter; ++pCenter) {
      tmp.clear();
      tmp.resize(dUromTimeIt.size(),0.0);
 
      for (int jState=0; jState<exactUpdatesAlpha.size(); ++jState) {
        tmp[jState] += initialConditionCentersDifProduct[mCenter][pCenter] * exactUpdatesAlpha[jState];
      }

      for (int iCluster=0; iCluster<nClusters; ++iCluster) {
        for (int jState=0; jState<exactUpdatesBeta[iCluster].size(); ++jState) {
          tmp[jState] += refStateCentersDifProduct[iCluster][mCenter][pCenter] * exactUpdatesBeta[iCluster][jState];
        }
      }

      for (int iCluster=0; iCluster<nClusters; ++iCluster) {
        if (exactUpdatesN[iCluster].size()>0) {
          for (int jState=0; jState<dUromTimeIt.size(); ++jState) {
            for (int iState=0; iState<exactUpdatesN[iCluster].size(); ++iState) {
              tmp[jState] += stateBasisCentersDifProduct[iCluster][mCenter][pCenter][iState]*exactUpdatesN[iCluster][iState][jState];
            }
          }
        }
      }

      for (int iState=0; iState<dUromTimeIt.size(); ++iState) {        
        distanceComparisons[mCenter][pCenter] += tmp[iState] * dUromTimeIt[iState];
      }
    }
  }

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::incrementDistanceComparisonsForApproxUpdates(Vec<double> dUromTimeIt, int currentCluster) {

  if (nClusters == 1) return;

  for (int mCenter=1; mCenter<nClusters; ++mCenter) {
    for (int pCenter=0; pCenter<mCenter; ++pCenter) {
      for (int iState=0; iState<dUromTimeIt.size(); ++iState) {
        distanceComparisons[mCenter][pCenter] += hForFastDistComp[mCenter][pCenter][iState] * dUromTimeIt[iState];
      }
      for (int iKrylov=0; iKrylov<nKrylov; ++iKrylov) {
        // account for Krylov bases
      }
      for (int iSens=0; iSens<nSens; ++iSens) {
        // account for sensitivity bases
      }
    }
  }


}

//----------------------------------------------------------------------------------

template<int dim>
int NonlinearRom<dim>::readSnapshotFiles(const char* snapType, bool preprocess) {

  // Check for snapshot command file
  char *vecFile;
  bool typeIsState = false;
  if (strcmp(snapType, "state")==0) {
    vecFile = new char[strlen(ioData->input.prefix) + strlen(ioData->input.stateSnapFile) + 1];
    sprintf(vecFile, "%s%s", ioData->input.prefix, ioData->input.stateSnapFile);
    typeIsState = true;
  } else if (strcmp(snapType,"sensitivity")==0) {
    vecFile = new char[strlen(ioData->input.prefix) + strlen(ioData->input.sensitivitySnapFile) + 1];
    sprintf(vecFile, "%s%s", ioData->input.prefix, ioData->input.sensitivitySnapFile);
  } else if (strcmp(snapType,"approxMetricState")==0) {
    vecFile = new char[strlen(ioData->input.prefix) + strlen(ioData->input.approxMetricStateSnapFile) + 1];
    sprintf(vecFile, "%s%s", ioData->input.prefix, ioData->input.approxMetricStateSnapFile);
  } else if (strcmp(snapType,"approxMetricNonlinear")==0) {
    vecFile = new char[strlen(ioData->input.prefix) + strlen(ioData->input.approxMetricNonlinearSnapFile) + 1];
    sprintf(vecFile, "%s%s", ioData->input.prefix, ioData->input.approxMetricNonlinearSnapFile);
  } else if (strcmp(snapType,"greedyData")==0) {
    vecFile = new char[strlen(ioData->input.prefix) + strlen(ioData->input.greedyDataFile) + 1];
    sprintf(vecFile, "%s%s", ioData->input.prefix, ioData->input.greedyDataFile);
//  } else if (strcmp(snapType,"residual")==0) {
//    vecFile = new char[strlen(ioData->input.prefix) + strlen(ioData->input.residualSnapFile) + 1];
//    sprintf(vecFile, "%s%s", ioData->input.prefix, ioData->input.residualSnapFile);
  } else if (strcmp(snapType,"projError")==0) {
    vecFile = new char[strlen(ioData->input.prefix) + strlen(ioData->input.projErrorSnapFile) + 1];
    sprintf(vecFile, "%s%s", ioData->input.prefix, ioData->input.projErrorSnapFile);
  } else if (strcmp(snapType,"initialClusterCenters")==0) {
    vecFile = new char[strlen(ioData->input.prefix) + strlen(ioData->input.initialClusterCentersFile) + 1];
    sprintf(vecFile, "%s%s", ioData->input.prefix, ioData->input.initialClusterCentersFile);
  } else {
    fprintf(stderr, "*** Error: unexpected snapshot type %s\n", snapType);
    exit (-1);
  }

  FILE *inFP = fopen(vecFile, "r");
  if (!inFP)  {
    fprintf(stderr, "*** Error: No snapshots FILES in \"%s\"\n", vecFile);
    exit (-1);
  }

  int nData, _n;
  _n = fscanf(inFP, "%d",&nData);
  com->fprintf(stdout, "Reading snapshots from %d files \n",nData);

  FILE *inFPref; 
  char refFile1[500];
  char *refFile;
  if (typeIsState) {
    stateSnapsFromFile.clear();
    stateSnapsFromFile.resize(nData, 0);
    if (!(strcmp(ioData->input.multiStateSnapRefSolution,"")==0) && !duplicateSnaps) {
      refFile = new char[strlen(ioData->input.prefix) + strlen(ioData->input.multiStateSnapRefSolution)+1];
      sprintf(refFile, "%s%s", ioData->input.prefix, ioData->input.multiStateSnapRefSolution);

      inFPref = fopen(refFile, "r");
      if (!inFPref)  {
        fprintf(stderr, "*** Error: No references FILES in \"%s\"\n", vecFile);
        exit (-1);
      }
      int nRef;
      _n = fscanf(inFPref, "%d",&nRef);
      if (nRef != nData) {
        fprintf(stdout,"Number of references must be equal to number of snapshot files.  Exiting...\n");
        exit(-1);
      }
    }
  }

  char** snapFile = new char*[nData];
  for (int iData = 0; iData < nData; ++iData)
    snapFile[iData] = new char[500];
  char snapFile1[500];
  int* numSnaps = new int[nData];
  int* startSnaps = new int[nData];
  int* endSnaps = new int[nData];
  int* sampleFreq = new int[nData];
  double* snapWeight = new double[nData];
  int nSnap, iStart, iEnd, iFreq;
  double weight;

  if (typeIsState && (ioData->romOffline.rob.state.dataCompression.energyOnly == DataCompressionData::ENERGY_ONLY_TRUE))
    com->fprintf(stderr, "*** Warning: EnergyOnly is not supported for multiple bases\n");

  // read snapshot command file
  for (int iData = 0; iData < nData; ++iData){
    _n = fscanf(inFP, "%s %d %d %d %lf", snapFile1,&iStart,&iEnd,&iFreq,&weight);
    if (_n<5) {
      fprintf(stderr, "*** Error: snapshot file %s is not formatted properly (path startSnap endSnap freq weight)\n", vecFile);
      exit(-1);
    }
    if (iStart < 1) iStart = 1;
    if (iEnd <= 0) {
      iEnd = 0;
    } else if (iEnd < iStart) {
      fprintf(stderr, "*** Error in %s: endSnap (%d) < startSnap (%d)\n", vecFile, iEnd, iStart);
      exit(-1);
    }
    if (iFreq < 1) iFreq = 1;
    //numSnaps[iData] = nSnap;
    strcpy(snapFile[iData],snapFile1);
    startSnaps[iData] = iStart - 1;
    endSnaps[iData] = iEnd;
    sampleFreq[iData] = iFreq;
    snapWeight[iData] = weight;
    com->fprintf(stdout, " ... Reading snapshots from %s \n", snapFile[iData]);

    if (typeIsState) {
      if (!(strcmp(ioData->input.multiStateSnapRefSolution,"")==0) && !duplicateSnaps) {
        _n = fscanf(inFPref, "%s", refFile1);
        mapSnapToRef[std::string(snapFile1)]=std::string(refFile1);
      }
    }
  }

  delete [] vecFile;
  vecFile = NULL;

  // compute the total number of snapshots
  int nTotSnaps = 0;
  int dummyStep = 0;
  double dummyTag = 0.0;
  for (int iData = 0; iData < nData; ++iData) {
    bool status = domain.template readTagFromFile<double, dim>(snapFile[iData], dummyStep, &dummyTag, &(numSnaps[iData]));

    if (!status) {
      fprintf(stderr, "*** Error: could not read snapshots from \"%s\" \n", snapFile[iData]);
      exit(-1);
    }

    if ((endSnaps[iData]==0) || (endSnaps[iData]>numSnaps[iData]))
      endSnaps[iData]=numSnaps[iData];
    for (int iSnap = startSnaps[iData]; iSnap<endSnaps[iData]; ++iSnap) {
      if (iSnap % sampleFreq[iData] == 0) {
        ++nTotSnaps;
        if (typeIsState) ++stateSnapsFromFile[iData]; 
      }
    }
  }

  bool subtractRefSol = false;
  if (preprocess) {
    if (ioData->romOffline.rob.relativeProjectionError.subtractRefSol) {
      subtractRefSol = true;
      if (!(ioData->input.stateSnapRefSolution)) {
        fprintf(stderr, "*** Error: Reference solution not found \n");
        exit (-1);
      }
      this->readReferenceState();
    }
  }

  VecSet< DistSVec<double, dim> >* snapshots;
  if (strcmp(snapType,"initialClusterCenters")==0) {
    if (nTotSnaps!=nClusters) {
      fprintf(stderr, "*** Error: number of specified centers (%d) does not match number of clusters (%d)\n",nTotSnaps,nClusters);
      exit (-1);
    }
    if (clusterCenters) delete clusterCenters;
    clusterCenters = new VecSet< DistSVec<double, dim> >(nTotSnaps, domain.getNodeDistInfo()); 
    snapshots = clusterCenters;
  } else {
    if (snap) delete snap;
    snap = new VecSet< DistSVec<double, dim> >(nTotSnaps, domain.getNodeDistInfo());
    snapshots = snap;
  }

  DistSVec<double, dim>* snapBufOld = new DistSVec<double, dim>(domain.getNodeDistInfo());
  DistSVec<double, dim>* snapBufNew = new DistSVec<double, dim>(domain.getNodeDistInfo());

  originalSnapshotLocation.clear();

  *snapBufOld = 0.0;
  *snapBufNew = 0.0;

  std::vector<double> tags;
  tags.clear();

  double tagOld;
  double tagNew;

  int numCurrentSnapshots = 0;
  bool status;

  if (subtractRefSol) {
    *snapBufOld = *(this->snapRefState);
    delete (this->snapRefState);
    (this->snapRefState) = NULL;
  }

  for (int iData=0; iData < nData; ++iData){
    // read in Snapshot Vectors
    for (int iSnap = startSnaps[iData]; iSnap<endSnaps[iData]; ++iSnap) {
      if (iSnap % sampleFreq[iData] == 0) { //TODO ignore 
        // snapshot must be between startSnaps and endSnaps, and a multiple of sampleFreq. 
        status = domain.readVectorFromFile(snapFile[iData], iSnap, &tagNew, *snapBufNew);
        (*snapshots)[numCurrentSnapshots] = *snapBufNew - *snapBufOld;
        double snapNorm = (*snapshots)[numCurrentSnapshots].norm();
        if (snapNorm > ioData->romOffline.rob.clustering.snapshotNormTolerance) {
          if (snapWeight[iData]>0.0 && snapWeight[iData]!=1.0) {
            (*snapshots)[numCurrentSnapshots] *= snapWeight[iData];
            com->fprintf(stderr, "*** Warning: basis updates should not be used for bases built with weighted snapshots (normalizing the snapshot matrix is supported, however)\n");  // supporting updates for bases built with weighted snapshots would require storing the weights when clustering (and not just the snapshots).  Since normalization happens after clustering there's no need to store anything in that case.
            if (!preprocess && nClusters>1) {
              fprintf(stderr, "*** Error: Clustering weighted snapshots is probably a bad idea.  Exiting...\n");
              exit(-1);
            }
          }
          tags.push_back(tagNew);
          originalSnapshotLocation.push_back(std::make_pair(string(snapFile[iData]),iSnap));
          ++numCurrentSnapshots;
        } else {
          com->fprintf(stderr, "*** Warning: skipping a snapshot with a norm of %e (iData=%d, iSnap=%d)\n",snapNorm,iData,iSnap);
          --stateSnapsFromFile[iData];
        }    
      }
    }
  }

  if (numCurrentSnapshots <  nTotSnaps) {
    nTotSnaps = numCurrentSnapshots;
    VecSet< DistSVec<double, dim> > *newSnap = new VecSet< DistSVec<double, dim> >(nTotSnaps, domain.getNodeDistInfo());
    for (int iVec=0; iVec<nTotSnaps; ++iVec)
      (*newSnap)[iVec] = (*snapshots)[iVec];
    delete snapshots;
    snapshots = newSnap;
  }

  if (typeIsState) {
    stateSnapshotTags.resize(nData);
    int snapCount=0;
    for (int iFile=0;iFile<nData;++iFile) {
      stateSnapshotTags[iFile].clear();
      stateSnapshotTags[iFile].resize(stateSnapsFromFile[iFile], -1.0);
      for (int iSnap=0;iSnap<stateSnapsFromFile[iFile];++iSnap) {
        stateSnapshotTags[iFile][iSnap] = tags[snapCount];
        ++snapCount;
      }
    }
  }

  delete snapBufOld;
  snapBufOld = NULL;
  delete snapBufNew;
  snapBufNew = NULL;

  for (int iData=0; iData < nData; ++iData) {
    delete [] snapFile[iData];
  }
  delete [] snapFile;
  snapFile = NULL;
  delete [] numSnaps;
  numSnaps = NULL;
  delete [] startSnaps;
  startSnaps = NULL;
  delete [] endSnaps;
  endSnaps = NULL;
  delete [] sampleFreq;
  sampleFreq = NULL;
  delete [] snapWeight;
  snapWeight = NULL;

  if (typeIsState && !(strcmp(ioData->input.multiStateSnapRefSolution,"")==0) && !duplicateSnaps) {
    delete[] refFile;
  }

  return nData;

}
//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::outputClusteredSnapshots(const char* snapType)  { 
 
  int nTotSnaps = snap->numVectors();

  if (strcmp(snapType, "state")==0) { 

    createDirectories();
  
    if (com->cpuNum() == 0) {
  
      // output cluster index as ASCII file (original clusters before overlap)
      FILE *clustIndexFile;
      char *clustIndexPath = 0; 
      determinePath(indexName, -1, clustIndexPath); 
      com->fprintf(stdout, "\nWriting cluster index to disk\n");
  
      clustIndexFile = fopen(clustIndexPath, "wt");
      if (clustIndexFile) {
         com->fprintf(clustIndexFile,"Snapshot# OriginalCluster\n");
      } else {  
         fprintf(stdout,"***Error: Cannot open cluster index \"%s\"\n",clustIndexPath);
         exit(-1);
      }

      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        com->fprintf(clustIndexFile,"%d %d\n", iSnap, clusterIndex[iSnap]);
      }
      fclose (clustIndexFile);
      delete [] clustIndexPath;
      clustIndexPath = NULL;
 
      // output cluster snapshot map as ASCII file (clusters after overlap)
      if (clusterSnapshotMap) {
        FILE *clustMapFile;
        char *clustMapPath = 0;
        determinePath(mapName, -1, clustMapPath);
        com->fprintf(stdout, "\nWriting cluster-snapshot map to disk\n");
  
        clustMapFile = fopen(clustMapPath, "wt");
        com->fprintf(clustMapFile,"Snapshot# Cluster\n");
  
        for (int iCluster=0; iCluster<nClusters; ++iCluster) {
          for (int iSnap=0; iSnap<snapsInCluster[iCluster]; ++iSnap) {
            com->fprintf(clustMapFile,"%d %d\n", clusterSnapshotMap[iCluster][iSnap], iCluster);
          }
        }
        fclose (clustMapFile);
        delete [] clustMapPath;
        clustMapPath = NULL;
      }

      // output cluster connectivity as ASCII file
      if (clusterNeighbors) {
        FILE *clustConnFile;
        char *clustConnPath = 0;
        determinePath(connName, -1, clustConnPath);
        com->fprintf(stdout, "\nWriting cluster connectivity to disk\n");
  
        clustConnFile = fopen(clustConnPath, "wt");
        com->fprintf(clustConnFile,"Cluster#  ...Neighbors...\n");
        for (int iCluster=0; iCluster<nClusters; ++iCluster) {
          com->fprintf(clustConnFile,"%d", iCluster);
          for (int iNeighbor=0; iNeighbor<clusterNeighborsCount[iCluster]; ++iNeighbor) {
            com->fprintf(clustConnFile," %d", clusterNeighbors[iCluster][iNeighbor]);
          }
          com->fprintf(clustConnFile,"\n");
        }
        fclose (clustConnFile);
        delete [] clustConnPath;
        clustConnPath = NULL;
      }
    }

    com->barrier();

    // output cluster centers
    char *clustCentersPath = 0;
    determinePath(centersName, -1, clustCentersPath);
    com->fprintf(stdout, "\nWriting cluster centers to disk\n");
  
    for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      com->barrier();
      domain.writeVectorToFile(clustCentersPath, iCluster, double(snapsInCluster[iCluster]), (*clusterCenters)[iCluster]);
    }
    delete [] clustCentersPath;
    clustCentersPath = NULL;  

    // output nearest snap to each cluster
    char *nearestSnapsPath = 0;
    determinePath(nearestName, -1, nearestSnapsPath);
    com->fprintf(stdout, "\nWriting nearest snapshot to each center to disk\n");
  
    for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      com->barrier();
      domain.writeVectorToFile(nearestSnapsPath, iCluster, double(snapsInCluster[iCluster]), (*nearestSnapsToCenters)[iCluster]);
    }
    delete [] nearestSnapsPath;
    nearestSnapsPath = NULL;  

    // output clustered snapshots
    for (int iCluster=0; iCluster<nClusters; iCluster++) {
      char *snapshotsPath = 0;
      determinePath(stateSnapsName, iCluster, snapshotsPath);

      if (duplicateSnaps) {
        // write snapshots to a new binary file
        com->fprintf(stdout, "\nWriting %d snapshots to a new binary file for cluster %d\n", snapsInCluster[iCluster], iCluster);
  
        int numWritten = 0;
        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          for (int jSnap=0; jSnap<snapsInCluster[iCluster]; ++jSnap) {
            if (iSnap == clusterSnapshotMap[iCluster][jSnap]) {
              domain.writeVectorToFile(snapshotsPath, numWritten, double(numWritten), (*snap)[iSnap] );
              ++numWritten;
            }
          }
        } 
      } else {
        // output a small ASCII file with pointers to the original snapshot files (to avoid doubling storage)
        if (com->cpuNum() == 0) {
          FILE *asciiSnapshotsFile;
          com->fprintf(stdout, "\nWriting clustered snapshot info to %s\n", snapshotsPath);

          asciiSnapshotsFile = fopen(snapshotsPath, "wt");

          if (!asciiSnapshotsFile) {
             fprintf(stdout,"***Error: Cannot open snapshot file \"%s\"\n",snapshotsPath);
             exit(-1);
          }

          for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
            for (int jSnap=0; jSnap<snapsInCluster[iCluster]; ++jSnap) {
              if (iSnap == clusterSnapshotMap[iCluster][jSnap]) {
                com->fprintf(asciiSnapshotsFile,"%s %d\n", (char*)(originalSnapshotLocation[iSnap].first).c_str(), originalSnapshotLocation[iSnap].second);
              }
            }
          }
          fclose (asciiSnapshotsFile);

        }
      }
      delete [] snapshotsPath;
      snapshotsPath = NULL;
    }

    if ((nClusters>1) || (strcmp(sensitivitySnapsName, "")!=0) ||
        (ioData->romOffline.rob.state.dataCompression.computePOD==DataCompressionData::COMPUTE_POD_FALSE)) {
      com->fprintf(stdout, "\nFreeing memory for parallel SVD; read in snapshots as needed\n");
      delete snap;
      snap = NULL;
    }                 
    if (clusterCenters) delete clusterCenters;
    clusterCenters = NULL;
    if (nearestSnapsToCenters) delete nearestSnapsToCenters;
    nearestSnapsToCenters = NULL;
    if (clusterIndex) delete [] clusterIndex;
    clusterIndex = NULL;
  
    for (int iCluster=0; iCluster<nClusters; ++iCluster) 
    if (clusterSnapshotMap) delete [] clusterSnapshotMap[iCluster];
    if (clusterSnapshotMap) delete [] clusterSnapshotMap;
    clusterSnapshotMap = NULL;
    if (snapsInCluster) delete [] snapsInCluster;
    snapsInCluster = NULL;
  
    for (int iCluster=0; iCluster<nClusters; ++iCluster) 
      if (clusterNeighbors) delete [] clusterNeighbors[iCluster];
    if (clusterNeighbors) delete [] clusterNeighbors;
    clusterNeighbors = NULL;
    if (clusterNeighborsCount) delete [] clusterNeighborsCount;
    clusterNeighborsCount = NULL;

  } else if (strcmp(snapType, "sensitivity")==0) {

    // output sensitivities

    char *sensitivityPath = 0;
    determinePath(sensitivitySnapsName, -2, sensitivityPath);

    if (duplicateSnaps) {
      com->fprintf(stdout, "\nWriting %d snapshots to sensitivity cluster\n", nTotSnaps);

      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        domain.writeVectorToFile(sensitivityPath, iSnap, double(iSnap), (*snap)[iSnap] );
      }

    } else {
      if (com->cpuNum() == 0) {
        FILE *asciiSnapshotsFile;
        com->fprintf(stdout, "\nWriting clustered sensitivity snapshot info to disk\n");

        asciiSnapshotsFile = fopen(sensitivityPath, "wt");

        if (!asciiSnapshotsFile) {
           fprintf(stdout,"***Error: Cannot open snapshot file \"%s\"\n",sensitivityPath);
           exit(-1);
        }

        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          com->fprintf(asciiSnapshotsFile,"%s %d\n", (char*)(originalSnapshotLocation[iSnap].first).c_str(), originalSnapshotLocation[iSnap].second);
        }
        fclose (asciiSnapshotsFile);
      }
    }
 
    delete [] sensitivityPath;
    sensitivityPath = NULL;  

    com->fprintf(stdout, "\nFreeing memory for parallel SVD; read in snapshots as needed\n");

    if (snap) delete snap;
    snap = NULL;
  }
}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredSnapshots(int iCluster, bool preprocess, const char *basisType, int first, int last, bool snapshotsAlreadyStored) {

  int nTotSnaps;
  int normalizeSnaps;

  // reconstruct cluster name
  char *snapshotsPath = 0;
  if (strcmp(basisType,"state")==0) {
    if (strcmp(stateSnapsName,"")==0) return;
    determinePath(stateSnapsName, iCluster, snapshotsPath);
    nTotSnaps = (last>=0) ? (last-first+1) : snapsInCluster[iCluster];
    com->fprintf(stdout, "\nReading %d snapshots from cluster %d \n", nTotSnaps, iCluster);
    normalizeSnaps = ioData->romOffline.rob.state.snapshots.normalizeSnaps;
    if ((ioData->romOffline.rob.state.snapshots.subtractRefState) && 
        !(ioData->romOffline.rob.state.snapshots.subtractCenters || ioData->romOffline.rob.state.snapshots.subtractNearestSnapsToCenters)){
      if (!(ioData->input.stateSnapRefSolution)) {
        fprintf(stderr, "*** Error: reference solution file not specified\n");
        exit (-1);
      }
      readReferenceState();
    }
  } else if (strcmp(basisType,"residual")==0) {
    if (strcmp(this->residualSnapsName,"")==0) return;
    determinePath(residualSnapsName, iCluster, snapshotsPath);
    nTotSnaps = (last>=0) ? (last-first+1) : 1;// value not stored -- resize when reading file
    com->fprintf(stdout, "\nReading residual snapshots from cluster %d \n", iCluster);
    normalizeSnaps = ioData->romOffline.rob.residual.snapshots.normalizeSnaps;
  } else if (strcmp(basisType,"jacAction")==0) {
    if (strcmp(this->jacActionSnapsName,"")==0) return;
    determinePath(jacActionSnapsName, iCluster, snapshotsPath);
    nTotSnaps = (last>=0) ? (last-first+1) : 1;// value not stored -- resize when reading file
    com->fprintf(stdout, "\nReading action-of-Jacobian snapshots from cluster %d \n", iCluster);
    normalizeSnaps = ioData->romOffline.rob.jacAction.snapshots.normalizeSnaps;
  } else if (strcmp(basisType,"krylov")==0) {
    if (strcmp(this->krylovSnapsName,"")==0) return;
    determinePath(krylovSnapsName, iCluster, snapshotsPath);
    nTotSnaps = (last>=0) ? (last-first+1) : 1;// value not stored -- resize when reading file
    com->fprintf(stdout, "\nReading krylov snapshots from cluster %d \n", iCluster);
    normalizeSnaps = ioData->romOffline.rob.krylov.snapshots.normalizeSnaps;
  } else if (strcmp(basisType,"sensitivity")==0) {
    if (strcmp(this->sensitivitySnapsName,"")==0) return;
    determinePath(sensitivitySnapsName, iCluster, snapshotsPath);
    nTotSnaps = (last>=0) ? (last-first+1) : 1;// value not stored -- resize when reading file
    com->fprintf(stdout, "\nReading sensitivity from snapshots\n");
    normalizeSnaps = ioData->romOffline.rob.sensitivity.snapshots.normalizeSnaps;
  } else {
    fprintf(stderr, "\n*** Error: invalid basisType (%s)\n", basisType);
    exit (-1);
  }

  if (snapshotsAlreadyStored) {

    com->fprintf(stdout, " ... global ROB requested: using previously stored state snapshots \n");
    nTotSnaps = snap->numVectors();

  } else {

    if (snap) delete snap;
    snap = new VecSet< DistSVec<double, dim> >(nTotSnaps, domain.getNodeDistInfo());
  
    // read in Snapshot Vectors
    double tmp;
    bool status = true; 
    int snapCount = 0;
    int snapIndex = first; 
    DistSVec<double, dim>* tmpSnap = new DistSVec<double, dim>(domain.getNodeDistInfo());
  
    if (duplicateSnaps) { 
      while (true) {
        status = domain.readVectorFromFile(snapshotsPath, snapIndex, &tmp, *tmpSnap);
        if (!status) {
          if (snapCount<nTotSnaps) {
            nTotSnaps = snapCount;
            VecSet< DistSVec<double, dim> >* snapNew = new VecSet< DistSVec<double, dim> >(nTotSnaps, domain.getNodeDistInfo());
            for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
              (*snapNew)[iSnap] = (*snap)[iSnap];
            }
            delete snap;
            snap = snapNew;
            snapNew = NULL;
          }
          break;
        }
        if (snapCount==nTotSnaps) {
          if (last>=0) break;
          ++nTotSnaps;
          VecSet< DistSVec<double, dim> >* snapNew = new VecSet< DistSVec<double, dim> >(nTotSnaps, domain.getNodeDistInfo());
          for (int iSnap=0; iSnap<(nTotSnaps-1); ++iSnap) {
            (*snapNew)[iSnap] = (*snap)[iSnap];
          }
          delete snap;
          snap = snapNew;
          snapNew = NULL;
        }
        (*snap)[snapCount] = *tmpSnap;
        ++snapCount;
        ++snapIndex;
      }
    } else {
      FILE *asciiSnapshotsFile;
      char snapFile[500];
      int snapNum;
      int asciiStatus = 2;
      int binaryStatus;
  
      asciiSnapshotsFile = fopen(snapshotsPath, "rt");
      for (int iSnap=0; iSnap<first; ++iSnap) {
        int asciiStatus = fscanf(asciiSnapshotsFile,"%s %d",snapFile,&snapNum);
        if (asciiStatus < 2) {
          fprintf(stdout, "***Error: Encountered error while reading snapshot file \"%s\"\n", snapshotsPath);
          exit(-1);
        }
      }
  
      while (true) {
  
        int asciiStatus = fscanf(asciiSnapshotsFile,"%s %d",snapFile,&snapNum);
        if (asciiStatus<2) {
          if (snapCount<nTotSnaps) {
            nTotSnaps = snapCount;
            VecSet< DistSVec<double, dim> >* snapNew = new VecSet< DistSVec<double, dim> >(nTotSnaps, domain.getNodeDistInfo());
            for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
              (*snapNew)[iSnap] = (*snap)[iSnap];
            }
            delete snap;
            snap = snapNew;
            snapNew = NULL;
          }
          break;
        }
  
        delete [] snapshotsPath;
        snapshotsPath = new char[500];
        strcpy(snapshotsPath, snapFile);
        
        binaryStatus = domain.readVectorFromFile(snapshotsPath, snapNum, &tmp, *tmpSnap);
        if (!binaryStatus) {
          fprintf(stdout, "***Error: Encountered error while reading snapshot #%d from file \"%s\"\n", snapNum, snapshotsPath);
          exit(-1);
        }
  
        if (snapCount==nTotSnaps) {
          if (last>=0) break;
          ++nTotSnaps;
          VecSet< DistSVec<double, dim> >* snapNew = new VecSet< DistSVec<double, dim> >(nTotSnaps, domain.getNodeDistInfo());
          for (int iSnap=0; iSnap<(nTotSnaps-1); ++iSnap) {
            (*snapNew)[iSnap] = (*snap)[iSnap];
          }
          delete snap;
          snap = snapNew;
          snapNew = NULL;
        }
        (*snap)[snapCount] = *tmpSnap;
        if (strcmp(basisType,"state")==0 && !(strcmp(ioData->input.multiStateSnapRefSolution,"")==0)) {
          binaryStatus = domain.readVectorFromFile(mapSnapToRef[std::string(snapshotsPath)].c_str(), 0, &tmp, *tmpSnap);
          if (!binaryStatus) {
            fprintf(stdout, "***Error: Encountered error while reading snapshot reference from file \"%s\"\n", snapshotsPath);
            exit(-1);
          }
          (*snap)[snapCount] -= *tmpSnap;
        } 
        ++snapCount;
      }
  
      fclose(asciiSnapshotsFile);
  
  
    }
  
    delete tmpSnap;
    tmpSnap = NULL;
    delete [] snapshotsPath;
    snapshotsPath = NULL;
  } // end read snapshots

  if (preprocess) {
    if (strcmp(basisType,"state")==0) { 
      if (ioData->romOffline.rob.state.snapshots.subtractNearestSnapsToCenters) {
        com->fprintf(stdout, " ... subtracting nearest snapshot to center from snapshots \n");
        if (ioData->romOffline.rob.state.snapshots.subtractCenters)
          com->fprintf(stderr, "*** Warning: Incompatible commands -- ignoring subtractCenters command \n");
        if (ioData->romOffline.rob.state.snapshots.subtractRefState) 
          com->fprintf(stderr, "*** Warning: Incompatible commands -- ignoring reference solution \n");
        --nTotSnaps; //
        VecSet< DistSVec<double, dim> >* snapNew = new VecSet< DistSVec<double, dim> >(nTotSnaps, domain.getNodeDistInfo());
        DistSVec<double, dim>* tmpSnap = new DistSVec<double, dim>(domain.getNodeDistInfo());
        outputClusteredReferenceState(iCluster, (*nearestSnapsToCenters)[iCluster]);
        com->barrier();
        int snapCount = 0;
        for (int iSnap=0; iSnap<=nTotSnaps; ++iSnap) {
          *tmpSnap = (*snap)[iSnap] - (*nearestSnapsToCenters)[iCluster];
          if ((*tmpSnap).norm() > 1e-6) { // one will be exactly zero; don't include this one
            (*snapNew)[snapCount] = *tmpSnap;
            ++snapCount;
          }
        }
        delete tmpSnap;
        tmpSnap = NULL;
        delete snap;
        snap = snapNew;
        snapNew = NULL;
      } else if (ioData->romOffline.rob.state.snapshots.subtractCenters) {
        com->fprintf(stdout, " ... subtracting cluster center from snapshots \n");
        if (ioData->romOffline.rob.state.snapshots.subtractRefState) 
          com->fprintf(stderr, "*** Warning: Incompatible commands -- ignoring reference solution \n");
        outputClusteredReferenceState(iCluster,(*clusterCenters)[iCluster]);
        com->barrier();
        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          (*snap)[iSnap] -= (*clusterCenters)[iCluster];
        }
      } else if (ioData->romOffline.rob.state.snapshots.subtractRefState) {
        com->fprintf(stdout, " ... subtracting reference solution from snapshots \n");
        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          (*snap)[iSnap] -= *snapRefState;
        }
        outputClusteredReferenceState(iCluster,(*snapRefState));
        com->barrier();
        delete snapRefState;
        snapRefState = NULL;
      } else {
        com->fprintf(stdout, " ... using raw states as snapshots\n");
        DistSVec<double, dim> zeroVec(domain.getNodeDistInfo());
        zeroVec = 0.0;
        outputClusteredReferenceState(iCluster, zeroVec);
        com->barrier();
      }
    }

  }

  cumulativeSnapWeights.clear();
  cumulativeSnapWeights.resize(nTotSnaps, 1.0); // updated below if normalizing snapshots
  if (preprocess && normalizeSnaps) {
    com->fprintf(stdout, " ... normalizing snapshots \n");
    double tmpNorm;
    double weight;
    double maxWeight = 1.0;
    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
      tmpNorm = ((*snap)[iSnap]).norm();
      weight = (tmpNorm > 1e-14) ? 1.0/tmpNorm : -1.0;
      maxWeight = (weight > maxWeight) ? weight : maxWeight;
      cumulativeSnapWeights[iSnap] = weight;
    }
    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) { 
      cumulativeSnapWeights[iSnap] = (cumulativeSnapWeights[iSnap]>0) ? cumulativeSnapWeights[iSnap] : maxWeight;
      //com->fprintf(stdout, " ... debugging: normalization weight for snapshot %d = %e\n", iSnap, cumulativeSnapWeights[iSnap]);
      (*snap)[iSnap] *= cumulativeSnapWeights[iSnap];
    }
  }

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::outputClusteredReferenceState(int iCluster, DistSVec<double, dim>& ref) {
  // Automatically stores the reference snapshot for cluster iCluster.
  // This functionality is in place to reduce the user's workload.

  char *refStatePath = 0;
  determinePath(refStateName, iCluster, refStatePath);

  com->fprintf(stdout, " ... storing reference state to disk for later use\n", iCluster);
  domain.writeVectorToFile(refStatePath, 0, double(iCluster), ref);

  delete [] refStatePath;
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredReferenceState(int iCluster, const char* refType) {
  // This function reads in the automatically stored reference snapshot for a cluster.
  // By storing these reference snapshots and reading them automatically it reduces the user's workload.

  if (Uref) {
    delete Uref;
    Uref = NULL;
  }

  if (storedAllOnlineQuantities || storedAllOfflineQuantities) {
    com->fprintf(stdout, " ... loading snapshot reference state for cluster %d\n", iCluster);
    Uref = new DistSVec<double, dim>(domain.getNodeDistInfo());
    *Uref = (*allRefStates)[iCluster];
    return;
  }

  char *refStatePath = 0;
  if ((ioData->problem.alltype == ProblemData::_STEADY_NONLINEAR_ROM_POST_ 
           || ioData->problem.alltype == ProblemData::_UNSTEADY_NONLINEAR_ROM_POST_) 
           && strcmp(surfaceRefStateName,"")!=0) {
      determinePath(surfaceRefStateName, iCluster, refStatePath);
  } else if (strcmp(refType,"state")==0) {
      determinePath(refStateName, iCluster, refStatePath);
  } else if (strcmp(refType,"sampledState")==0) {
      determinePath(sampledRefStateName, iCluster, refStatePath);
  } else {
      fprintf(stderr, "\n*** Error: invalid refType (%s)\n", refType);
      exit (-1);
  }

  double tmp;
  bool status;
 
  //com->fprintf(stdout, "Reading reference snapshot %s\n", refStatePath);
  Uref = new DistSVec<double, dim>(domain.getNodeDistInfo());
  status = domain.readVectorFromFile(refStatePath, 0, &tmp, *Uref);

  delete [] refStatePath;
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusterCenters(const char* centersType) {

  char *clustCentersPath = 0;
  if ((ioData->problem.alltype == ProblemData::_STEADY_NONLINEAR_ROM_POST_
           || ioData->problem.alltype == ProblemData::_UNSTEADY_NONLINEAR_ROM_POST_)
           && strcmp(surfaceCentersName,"")!=0) {
      determinePath(surfaceCentersName, -1, clustCentersPath);
  } else if (strcmp(centersType,"centers")==0) {
      determinePath(centersName, -1, clustCentersPath);
  } else if (strcmp(centersType,"sampledCenters")==0) {
      determinePath(sampledCentersName, -1, clustCentersPath);
  } else {
      fprintf(stderr, "\n*** Error: invalid centersType (%s)\n", centersType);
      exit (-1);
  }

  if (nClusters <= 0) {
    fprintf(stderr, "\n*** Error: invalid value for NumClusters (%d)\n", nClusters);
    exit(-1);
  }

  // read centers

  if (clusterCenters) delete clusterCenters;  
  clusterCenters = new VecSet< DistSVec<double, dim> >(nClusters, domain.getNodeDistInfo());

  if (snapsInCluster) delete [] snapsInCluster;
  snapsInCluster = new int[nClusters];

  com->fprintf(stdout, "\nReading cluster centers\n");
  bool status;
  int expectedClusters = nClusters;
  nClusters = 0;
  double tmp;
  DistSVec<double, dim>* tmpVec = new DistSVec<double, dim>(domain.getNodeDistInfo());

  while (true) {

    status = domain.readVectorFromFile(clustCentersPath, nClusters, &tmp, *tmpVec);
    if (!status) break;
 
    ++nClusters;

    if (nClusters > expectedClusters) {
      fprintf(stderr, "\n*** Error: found more clusters than expected (NumClusters was specified as %d)\n", expectedClusters);
      exit(-1);
    }   

    VecSet< DistSVec<double, dim> >* clusterCentersNew =  new VecSet< DistSVec<double, dim> >(nClusters, domain.getNodeDistInfo());
    int* snapsInClusterNew = new int[nClusters];

    for (int iCluster=0; iCluster<(nClusters-1); ++iCluster) {
      (*clusterCentersNew)[iCluster] = (*clusterCenters)[iCluster];
      snapsInClusterNew[iCluster] = snapsInCluster[iCluster];
    }

    (*clusterCentersNew)[nClusters-1] = *tmpVec;
    snapsInClusterNew[nClusters-1] = int(tmp);

    delete clusterCenters;
    delete [] snapsInCluster;

    clusterCenters = clusterCentersNew;
    snapsInCluster = snapsInClusterNew;

    clusterCentersNew = NULL;
    snapsInClusterNew = NULL;
  }

  com->fprintf(stdout, "\n ... found %d clusters in this database\n", nClusters);
  com->barrier();
  delete tmpVec;
  tmpVec = NULL;
  delete [] clustCentersPath;
  clustCentersPath = NULL;

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readNearestSnapsToCenters() {
  // NOTE: this function must only be run after readClusterCenters (this is because the
  // cluster centers file defines snapsInCluster and nClusters)

  char *nearestSnapsPath = 0;
  determinePath(nearestName, -1, nearestSnapsPath);

  if (nearestSnapsToCenters) delete nearestSnapsToCenters; 
  nearestSnapsToCenters = new VecSet< DistSVec<double, dim> >(nClusters, domain.getNodeDistInfo());

  bool status;

  com->fprintf(stdout, "\nReading closest snapshot to each center\n");

  // read in Snapshot Vectors
  double tmp;
  for (int iCluster = 0; iCluster<nClusters; ++iCluster) {
    status = domain.readVectorFromFile(nearestSnapsPath, iCluster, &tmp, (*nearestSnapsToCenters)[iCluster]);
    if (!status) {
       fprintf(stderr,"*** Error: readNearestSnapsToCenters attempted to read %d vecs from a file (%s) with only %d.\n", 
               nClusters,nearestSnapsPath,iCluster-1);
       exit(-1);
    }
  }

  delete [] nearestSnapsPath;
  
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredBasis(int iCluster, const char* basisType, bool relProjError) {

  if (storedAllOnlineQuantities || storedAllOfflineQuantities) {
    if (basis) delete basis;
    if (sVals) delete sVals;
    if ((strcmp(basisType,"state")==0) || (strcmp(basisType,"sampledState")==0)) {
      com->fprintf(stdout, " ... loading state ROB for cluster %d\n", iCluster);
      nState = allStateBases[iCluster]->numVectors();
      basis = new VecSet< DistSVec<double, dim> >(nState, domain.getNodeDistInfo());
      for (int iVec=0; iVec<nState; ++iVec)
        (*basis)[iVec] = (*(allStateBases[iCluster]))[iVec]; 
      sVals = new std::vector<double>;
      *sVals = *(allStateSVals[iCluster]);
      nBuffer = (allNBuffer.size()>0) ? allNBuffer[iCluster] : 0;
    } else if ((strcmp(basisType,"sensitivity")==0) || (strcmp(basisType,"sampledSensitivity")==0) ) {
      com->fprintf(stdout, " ... loading sensitivity ROB\n");
      nSens = sensitivityBasis->numVectors();
      basis = new VecSet< DistSVec<double, dim> >(nSens, domain.getNodeDistInfo());
      for (int iVec=0; iVec<nSens; ++iVec)
        (*basis)[iVec] = (*sensitivityBasis)[iVec]; 
      sVals = new std::vector<double>;
      *sVals = *sensitivitySVals;
    } else if ((strcmp(basisType,"krylov")==0) || (strcmp(basisType,"sampledKrylov")==0)) {
      com->fprintf(stdout, " ... loading Krylov ROB for cluster %d\n", iCluster);
      nKrylov = allKrylovBases[iCluster]->numVectors();
      basis = new VecSet< DistSVec<double, dim> >(nKrylov, domain.getNodeDistInfo());
      for (int iVec=0; iVec<nKrylov; ++iVec)
        (*basis)[iVec] = (*(allKrylovBases[iCluster]))[iVec]; 
      sVals = new std::vector<double>;
      *sVals = *(allKrylovSVals[iCluster]);
    }
    return;
  }

  int maxDimension;
  int minDimension;
  double energyTol;
  double bufferEnergyTol = 0.0;

  if (ioData->problem.type[ProblemData::NLROMOFFLINE]) {
    if (relProjError) {
      if (strcmp(basisType,"state")==0) {
        maxDimension = ioData->romOffline.rob.relativeProjectionError.maxDimension;
        minDimension = ioData->romOffline.rob.relativeProjectionError.minDimension;
        energyTol = ioData->romOffline.rob.relativeProjectionError.energy;
      } else if (strcmp(basisType,"residual")==0) {
        maxDimension = ioData->romOffline.rob.relativeProjectionError.maxDimension;
        minDimension = ioData->romOffline.rob.relativeProjectionError.minDimension;
        energyTol = ioData->romOffline.rob.relativeProjectionError.energy;
      } else if (strcmp(basisType,"jacAction")==0) {
        maxDimension = ioData->romOffline.rob.relativeProjectionError.maxDimension;
        minDimension = ioData->romOffline.rob.relativeProjectionError.minDimension;
        energyTol = ioData->romOffline.rob.relativeProjectionError.energy;
      } else if (strcmp(basisType,"sensitivity")==0) {
        maxDimension = ioData->romOffline.rob.relativeProjectionError.sensitivity.maxDimension;
        minDimension = ioData->romOffline.rob.relativeProjectionError.sensitivity.minDimension;
        energyTol = ioData->romOffline.rob.relativeProjectionError.sensitivity.energy;
      } else if (strcmp(basisType,"krylov")==0) {
        maxDimension = ioData->romOffline.rob.relativeProjectionError.krylov.maxDimension;
        minDimension = ioData->romOffline.rob.relativeProjectionError.krylov.minDimension;
        energyTol = ioData->romOffline.rob.relativeProjectionError.krylov.energy;
      }
    } else { 
      if (strcmp(basisType,"state")==0) {
        maxDimension = ioData->romOffline.gappy.maxDimensionState;
        minDimension = ioData->romOffline.gappy.minDimensionState;
        energyTol = ioData->romOffline.gappy.energyState;
      } else if (strcmp(basisType,"sensitivity")==0) {
        maxDimension = ioData->romOffline.gappy.maxDimensionSensitivity;
        minDimension = ioData->romOffline.gappy.minDimensionSensitivity;
        energyTol = ioData->romOffline.gappy.energySensitivity;
      } else if (strcmp(basisType,"krylov")==0) {
        maxDimension = ioData->romOffline.gappy.maxDimensionKrylov;
        minDimension = ioData->romOffline.gappy.minDimensionKrylov;
        energyTol = ioData->romOffline.gappy.energyKrylov;
      } else if (strcmp(basisType,"residual")==0) {
        maxDimension = ioData->romOffline.gappy.maxDimensionResidual;
        minDimension = ioData->romOffline.gappy.minDimensionResidual;
        energyTol = ioData->romOffline.gappy.energyResidual;
      } else if (strcmp(basisType,"jacAction")==0) {
        if (ioData->romOffline.gappy.maxDimensionJacAction <= 0) {
          com->fprintf(stdout, "*** Warning: JacAction greedy parameters not specified; using Residual parameters\n");
          maxDimension = ioData->romOffline.gappy.maxDimensionResidual;
          minDimension = ioData->romOffline.gappy.minDimensionResidual;
          energyTol = ioData->romOffline.gappy.energyResidual;
        } else {
          maxDimension = ioData->romOffline.gappy.maxDimensionJacAction;
          minDimension = ioData->romOffline.gappy.minDimensionJacAction;
          energyTol = ioData->romOffline.gappy.energyJacAction;
        }
      }
    }
  } else { // ONLINE
    if ((strcmp(basisType,"state")==0) || (strcmp(basisType,"sampledState")==0)) {
      maxDimension = ioData->romOnline.maxDimension;
      minDimension = ioData->romOnline.minDimension;
      energyTol = ioData->romOnline.energy;
      bufferEnergyTol = ioData->romOnline.bufferEnergy;
    } else if ((strcmp(basisType,"sensitivity")==0) || (strcmp(basisType,"sampledSensitivity")==0) ) {
      maxDimension = ioData->romOnline.sensitivity.maxDimension;
      minDimension = ioData->romOnline.sensitivity.minDimension;
      energyTol = ioData->romOnline.sensitivity.energy;
    } else if ((strcmp(basisType,"krylov")==0) || (strcmp(basisType,"sampledKrylov")==0)) {
      maxDimension = ioData->romOnline.krylov.maxDimension;
      minDimension = ioData->romOnline.krylov.minDimension;
      energyTol = ioData->romOnline.krylov.energy;
    }
  }

  char* singValsPath = 0;
  char* basisPath = 0;
  
  if (strcmp(basisType,"state")==0) {
      determinePath(stateSingValsName, iCluster, singValsPath);
      if ((ioData->problem.alltype == ProblemData::_STEADY_NONLINEAR_ROM_POST_
           || ioData->problem.alltype == ProblemData::_UNSTEADY_NONLINEAR_ROM_POST_)
           && strcmp(surfaceStateBasisName,"")!=0) {
        determinePath(surfaceStateBasisName, iCluster, basisPath);
      } else {
        determinePath(stateBasisName, iCluster, basisPath);
      }
  } else if (strcmp(basisType,"sampledState")==0) {
      determinePath(stateSingValsName, iCluster, singValsPath);
      if ((ioData->problem.alltype == ProblemData::_STEADY_NONLINEAR_ROM_POST_
           || ioData->problem.alltype == ProblemData::_UNSTEADY_NONLINEAR_ROM_POST_)
           && strcmp(surfaceStateBasisName,"")!=0) {
        determinePath(surfaceStateBasisName, iCluster, basisPath);
      } else {
        determinePath(sampledStateBasisName, iCluster, basisPath);
      }
  } else if (strcmp(basisType,"residual")==0) { 
      determinePath(residualSingValsName, iCluster, singValsPath);
      determinePath(residualBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"jacAction")==0) {  
      determinePath(jacActionSingValsName, iCluster, singValsPath);
      determinePath(jacActionBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"krylov")==0) {
      determinePath(krylovSingValsName, iCluster, singValsPath);
      determinePath(krylovBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"sampledKrylov")==0) {
      determinePath(krylovSingValsName, iCluster, singValsPath);
      determinePath(sampledKrylovBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"sensitivity")==0) {
      iCluster = -2;
      determinePath(sensitivitySingValsName, iCluster, singValsPath);
      determinePath(sensitivityBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"sampledSensitivity")==0) {
      iCluster = -2;
      determinePath(sensitivitySingValsName, iCluster, singValsPath);
      determinePath(sampledSensitivityBasisName, iCluster, basisPath);
  } else {
      fprintf(stderr, "*** Error: unrecognized basis type \"%s\" in readClusteredBasis()\n", basisType);
      exit (-1);
  }

  int step = 0;
  double tag = 0.0;
  int numSteps = 0;
  int status = domain.readTagFromFile<double, dim>(basisPath, step, &tag, &numSteps);

  if (!status)  {
    fprintf(stderr, "*** Error: unable to open file \"%s\" in readClusteredBasis()\n", basisPath);
    exit (-1);
  }

  if (minDimension<=1) minDimension = 1;
  if (maxDimension<=0) maxDimension = numSteps;

  FILE *singValFile = fopen(singValsPath, "r");

  if (!singValFile)  {
    fprintf(stderr, "*** Error: unable to open file \"%s\" in readClusteredBasis()\n", singValsPath);
    exit (-1);
  }

  delete [] singValsPath;
  singValsPath = NULL;

  int _n;
  int vecNumber;
  double percentEnergy;
  double tmpSVal;

  if (sVals) delete sVals;
  sVals = new std::vector<double>;
 
  nBuffer=0;

  for (int iVec=0; iVec<maxDimension; ++iVec) {
    _n = fscanf(singValFile,"%d %le %le", &vecNumber, &tmpSVal, &percentEnergy);
    if (_n == 3) {
      sVals->push_back(tmpSVal);
      if ((percentEnergy>=energyTol)&&((iVec+1)>=minDimension)) {
        if (percentEnergy<bufferEnergyTol) {
          ++nBuffer;
        } else { 
          break;
        }
      }
    } else if (feof(singValFile)) {
      break;
    } else {
      fprintf(stderr, "*** Error: fscanf interrupted by non-EOF error in readClusteredBasis() (file: \"%s\")\n", singValsPath);
      exit(-1);
    }
  }

  int basisSize = sVals->size();
  fclose(singValFile);

  // read in basis vectors

  com->fprintf(stdout, "\nReading basis %d\n", iCluster);

  if (basis) delete basis;
  basis = new VecSet< DistSVec<double, dim> >(basisSize, domain.getNodeDistInfo());

  int numBasisVecs = 0;

  for (int iVec = 0; iVec<basisSize; ++iVec) {
    status = domain.readVectorFromFile(basisPath, iVec, &tmpSVal, (*basis)[iVec]);
    if (!status) break;
    ++numBasisVecs; 
  }

  delete [] basisPath;
  basisPath = NULL;

  if (numBasisVecs != basisSize) {
    VecSet< DistSVec<double, dim> >* basisNew =  new VecSet< DistSVec<double, dim> >(numBasisVecs, domain.getNodeDistInfo());

    for (int iVec=0; iVec<numBasisVecs; ++iVec) {
      (*basisNew)[iVec]=(*basis)[iVec];
    }

    delete basis;
    basis = basisNew;
    basisNew = NULL; 

    if (nBuffer>0) {
      nBuffer = (numBasisVecs-(basisSize-nBuffer)>0) ? numBasisVecs-(basisSize-nBuffer) : 0 ;   
    }
  }    

  if (nBuffer>0) com->fprintf(stderr, " ... using a buffer of size %d for this basis\n", nBuffer);

  // for ASCII output files (clusterUsage and reducedCoords)
  if ((strcmp(basisType,"state")==0) || (strcmp(basisType,"sampledState")==0)) {
      nState = numBasisVecs;
  } else if ((strcmp(basisType,"krylov")==0) || (strcmp(basisType,"sampledKrylov")==0)) {
      nKrylov = numBasisVecs;
  } else if ((strcmp(basisType,"sensitivity")==0) || (strcmp(basisType,"sampledSensitivity")==0)) {
      nSens = numBasisVecs;
  }

  com->fprintf(stdout, "\n");

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredUpdateInfo(int iCluster, const char* basisType) {

  if (ioData->romOnline.basisUpdates == NonlinearRomOnlineData::UPDATES_OFF &&
      ioData->romOnline.projectSwitchStateOntoAffineSubspace == NonlinearRomOnlineData::PROJECT_ON) {
    readClusteredReferenceState(iCluster, basisType);
  } else {
    readClusteredReferenceState(iCluster, basisType);
    readClusteredColumnSumsV(iCluster, basisType);
  }
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readNonClusteredUpdateInfo(const char* sampledOrFull) {

  switch (ioData->romOnline.basisUpdates) {
    case (NonlinearRomOnlineData::UPDATES_OFF):
      if ((strcmp(sampledOrFull, "sampled")==0) 
           && (ioData->romOnline.projectSwitchStateOntoAffineSubspace!=NonlinearRomOnlineData::PROJECT_OFF))
        readProjectionInfo();
      break;
    case (NonlinearRomOnlineData::UPDATES_SIMPLE):
      break;
    case (NonlinearRomOnlineData::UPDATES_FAST_EXACT):
      readExactUpdateInfo();
      break;
    case (NonlinearRomOnlineData::UPDATES_FAST_APPROX):
      readApproxMetricStateLowRankFactor(sampledOrFull); 
      break;
    default:
      fprintf(stderr, "*** Error: Unexpected ROB updates method encountered in readNonClusteredUpdateInfo()\n");
      exit(-1);
  }
}




//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readProjectionInfo() {


  // Basis Basis Products (rob_i^T * rob_p)
  // std::vector<std::vector<std::vector<std::vector<double> > > > basisBasisProducts;  // [iCluster][pCluster][:][:]
  readClusteredInfoASCII(-1, "basisBasisProducts", NULL, NULL, NULL, &basisBasisProducts);

  // Basis Uref Products (rob_i^T * Uref_p)
  // std::vector<std::vector<std::vector<double> > > basisUrefProducts;  // [Cluster_Basis][Cluster_Uref][:]
  readClusteredInfoASCII(-1, "basisUrefProducts", NULL, NULL, &basisUrefProducts);

  checkInitialConditionScenario();

  if (specifiedIC) {
    // Basis Uic Products
    // std::vector<std::vector<double> > basisUicProducts;  // [iCluster][1:nPod] only precomputed if Uic specified
    readClusteredInfoASCII(-1, "basisUicProducts", NULL, &basisUicProducts);
  } else if (interpolatedMultiIC){
    // Uic was interpolated from multiple initial conditions, need to form basisUicProducts
    // basisMultiUicProducts;  // [iCluster][iSolution][1:nPod]
    readClusteredInfoASCII(-1, "basisMultiUicProducts", NULL, NULL, &basisMultiUicProducts);
    basisUicProducts.resize(nClusters);
    for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      basisUicProducts[iCluster].clear();
      basisUicProducts[iCluster].resize(basisMultiUicProducts[iCluster][0].size(),0.0);
      if (basisMultiUicProducts[iCluster].size() != interpWeightsForMultiIC.size()) {
        fprintf(stderr, "*** Error: number of interpolated solutions used "
                            "for this simulation does not match precomputed quantities\n");
        exit(-1);
      }
      for (int iSol=0; iSol<basisMultiUicProducts[iCluster].size(); ++iSol) {
        for (int iPod=0; iPod<basisMultiUicProducts[iCluster][iSol].size(); ++iPod) {
          basisUicProducts[iCluster][iPod] += basisMultiUicProducts[iCluster][iSol][iPod]
                                              * interpWeightsForMultiIC[iSol];
        }
      }
    } 
    basisMultiUicProducts.clear();
   
  } else {
    // Basis Componentwise Sums
    // std::vector<std::vector<std::vector<double> > > basisComponentwiseSums;  // [iCluster][iVec][1:dim]
    readClusteredInfoASCII(-1, "basisComponentwiseSums", NULL, NULL, &basisComponentwiseSums);
  }

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readExactUpdateInfo() {


  // Basis Basis Products (rob_i^T * rob_p)
  // std::vector<std::vector<std::vector<std::vector<double> > > > basisBasisProducts;  // [iCluster][pCluster][:][:]
  readClusteredInfoASCII(-1, "basisBasisProducts", NULL, NULL, NULL, &basisBasisProducts);

  // Basis Uref Products (rob_i^T * Uref_p)
  // std::vector<std::vector<std::vector<double> > > basisUrefProducts;  // [Cluster_Basis][Cluster_Uref][:]
  readClusteredInfoASCII(-1, "basisUrefProducts", NULL, NULL, &basisUrefProducts);

  // Uref Uref Products
  // std::vector<std::vector<double> > urefUrefProducts; //[iCluster][jCluster] symmetric (lower triangular)
  readClusteredInfoASCII(-1, "urefUrefProducts", NULL, &urefUrefProducts);

  checkInitialConditionScenario();

  if (specifiedIC) {

    // Basis Uic Products
    // std::vector<std::vector<double> > basisUicProducts;  // [iCluster][1:nPod] only precomputed if Uic specified
    readClusteredInfoASCII(-1, "basisUicProducts", NULL, &basisUicProducts);

    // Uref Uic Products
    // std::vector<double> urefUicProducts; // [iCluster] only precomputed if Uic specified
    readClusteredInfoASCII(-1, "urefUicProducts", &urefUicProducts);

  } else if (interpolatedMultiIC) {

    readClusteredInfoASCII(-1, "basisMultiUicProducts", NULL, NULL, &basisMultiUicProducts);
    basisUicProducts.resize(nClusters);
    for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      basisUicProducts[iCluster].clear();
      basisUicProducts[iCluster].resize(basisMultiUicProducts[iCluster][0].size(),0.0);
      if (basisMultiUicProducts[iCluster].size() != interpWeightsForMultiIC.size()) {
        fprintf(stderr, "*** Error: number of interpolated solutions used "
                            "for this simulation does not match precomputed quantities\n");
        exit(-1);
      }
      for (int iSol=0; iSol<basisMultiUicProducts[iCluster].size(); ++iSol) {
        for (int iPod=0; iPod<basisMultiUicProducts[iCluster][iSol].size(); ++iPod) {
          basisUicProducts[iCluster][iPod] += basisMultiUicProducts[iCluster][iSol][iPod]
                                              * interpWeightsForMultiIC[iSol];
        }
      }
    } 
    basisMultiUicProducts.clear();

    readClusteredInfoASCII(-1, "urefMultiUicProducts", NULL, &urefMultiUicProducts);
    urefUicProducts.clear();
    urefUicProducts.resize(nClusters,0.0);
    for (int iCluster=0; iCluster<nClusters; ++iCluster) {
      for (int iSol=0; iSol<urefMultiUicProducts[iCluster].size(); ++iSol) {
          urefUicProducts[iCluster] += urefMultiUicProducts[iCluster][iSol]
                                              * interpWeightsForMultiIC[iSol];
      }
    } 
    urefMultiUicProducts.clear();

  } else {
    // Uref Componentwise Sums
    // std::vector<std::vector<double> > urefComponentwiseSums; //[iCluster][1:dim]
    readClusteredInfoASCII(-1, "urefComponentwiseSums", NULL, &urefComponentwiseSums);

    // Basis Componentwise Sums
    // std::vector<std::vector<std::vector<double> > > basisComponentwiseSums;  // [iCluster][iVec][1:dim]
    readClusteredInfoASCII(-1, "basisComponentwiseSums", NULL, NULL, &basisComponentwiseSums);
  }


}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredColumnSumsV(int iCluster, const char* basisType) {

  if (columnSumsV) {
    delete columnSumsV;
    columnSumsV = NULL;
  }

  if ((strcmp(basisType, "state") == 0) || (strcmp(basisType, "sampledState") == 0) ) { 

    if (storedAllOnlineQuantities || storedAllOfflineQuantities) {
      com->fprintf(stdout, " ... loading columnSumsV for cluster %d\n", iCluster);
      columnSumsV = new std::vector<double>;
      *columnSumsV = *(allColumnSumsV[iCluster]);
      return;
    }

    char *basisUpdatePath = NULL;
    determinePath(simpleUpdateInfoName, iCluster, basisUpdatePath);

    com->fprintf(stdout, "\nReading update info for basis %d \n", iCluster);

    FILE *basisUpdateFile = fopen(basisUpdatePath, "r");

    if (!basisUpdateFile)  {
      fprintf(stderr, "*** Error: unable to open file \"%s\" in readClusteredColumnSumsV()\n", basisUpdatePath);
      exit (-1);
    }

    columnSumsV = new std::vector<double>;
    double tmp;

    while (true) {
      int _n = fscanf(basisUpdateFile, "%le", &tmp);
      if (_n == 1) {
        columnSumsV->push_back(tmp);
      } else if (feof(basisUpdateFile)) {
        break;
      } else {
        fprintf(stderr, "*** Error: fscanf interrupted by non-EOF error in readClusteredColumnSumsV() (file: \"%s\")\n", basisUpdatePath);
        exit(-1);
      }
    }

    fclose(basisUpdateFile);
    delete [] basisUpdatePath;
    basisUpdatePath = NULL;

  } else {
    fprintf(stderr, "*** Error: unexpected ROB type (%s) in readClusteredColumnSumsV()\n", basisType);
    exit(-1);
  }

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::outputClusteredBasis(int iCluster, int nTotSnaps, const char* basisType) {

  int podSize = basis->numVectors();

  char* singValsPath = 0;
  char* basisPath = 0;

  if (strcmp(basisType,"state")==0) {
      determinePath(stateSingValsName, iCluster, singValsPath);
      determinePath(stateBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"residual")==0) {
      determinePath(residualSingValsName, iCluster, singValsPath);
      determinePath(residualBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"jacAction")==0) {
      determinePath(jacActionSingValsName, iCluster, singValsPath);
      determinePath(jacActionBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"krylov")==0) {
      determinePath(krylovSingValsName, iCluster, singValsPath);
      determinePath(krylovBasisName, iCluster, basisPath);
  } else if (strcmp(basisType,"sensitivity")==0) {
      determinePath(sensitivitySingValsName, iCluster, singValsPath);
      determinePath(sensitivityBasisName, iCluster, basisPath);
  } else {
      fprintf(stderr, "\n*** Error: invalid basisType (%s)\n", basisType);
      exit (-1);
  }

  // output vectors
  com->fprintf(stdout, "\nWriting basis %d to disk\n", iCluster);
  for (int iVec=0; iVec<podSize; ++iVec) {
      domain.writeVectorToFile(basisPath, iVec, double(iVec), (*basis)[iVec] );
  }

  delete basis;
  basis = NULL;
  delete [] basisPath;
  basisPath = NULL;

  // write singular values to ASCII file
  com->fprintf(stdout, "\nWriting singular values to disk\n");

  double sValSqrdSum = 0;
  for (int iSnap = 0; iSnap<sVals->size(); ++iSnap) 
    sValSqrdSum += ((*sVals)[iSnap] * (*sVals)[iSnap]);
  double invSValSqrdSum = 1/sValSqrdSum;

  double sValSqrdPartialSum = 0;

  if (com->cpuNum() == 0) {
    FILE *singValFile = fopen(singValsPath, "wt");
    //com->fprintf(singValFile,"Vector# SVal SValSum\n");

    for (int iVec=0; iVec<sVals->size(); ++iVec) {
      sValSqrdPartialSum += ((*sVals)[iVec]*(*sVals)[iVec]*invSValSqrdSum);
      com->fprintf(singValFile,"%d %1.12e %1.12e\n", iVec, (*sVals)[iVec], sValSqrdPartialSum);
     }

    fclose (singValFile);
  }

 
  delete [] singValsPath;
  singValsPath = NULL;
  delete sVals;
  sVals = NULL;

  if (strcmp(basisType, "state") == 0) { 
    // write thinSVD update info to ASCII file
    com->fprintf(stdout, "\nWriting sums of right singular vectors to disk (for basis updates)\n");

    char *basisUpdatePath = 0;
    determinePath(simpleUpdateInfoName, iCluster, basisUpdatePath);

    if (com->cpuNum() == 0) {
      FILE *basisUpdateFile = fopen(basisUpdatePath, "wt");

      for (int iVec=0; iVec<(columnSumsV->size()); ++iVec) {
        com->fprintf(basisUpdateFile,"%le\n", (*columnSumsV)[iVec]);
      }

      fclose (basisUpdateFile);
    }

    delete [] basisUpdatePath;
    basisUpdatePath = NULL;
    delete columnSumsV;
    columnSumsV = NULL;
  }


}

//----------------------------------------------------------------------------------

template<int dim> 
void NonlinearRom<dim>::determineFileName(const char* fileNameInput, const char* fileNameExtension, const char* prefix, char*& fileName) { 

  if (strcmp(fileNameInput,"") == 0) {
    if (strcmp(prefix,"") == 0) {
      fileName = new char[1];
      fileName[0] = 0;
    } else {
      fileName = new char[strlen(prefix) + 1 + strlen(fileNameExtension) + 1];
      sprintf(fileName, "%s.%s", prefix, fileNameExtension);
    }
  } 
  else {
    fileName = new char[strlen(fileNameInput) + 1];
    sprintf(fileName, "%s", fileNameInput); 
  }  
} 

//----------------------------------------------------------------------------------

template<int dim> 
void NonlinearRom<dim>::determinePrefixName(const char* prefixInput, const char* prefixDefault, char*& prefix) { 

  if (strcmp(prefixInput,"") == 0) {
    if (strcmp(prefixDefault,"") == 0) {
      prefix = new char[1];
      prefix[0] = 0;
    } else {
      prefix = new char[strlen(prefixDefault) + 1];
      sprintf(prefix, "%s", prefixDefault);
    }
  } 
  else {
    prefix = new char[strlen(prefixInput) + 1];
    sprintf(prefix, "%s", prefixInput); 
  }  
} 


//----------------------------------------------------------------------------------

template<int dim> 
void NonlinearRom<dim>::determinePath(char* fileName, int iCluster, char*& path) { 

  if (iCluster == -1) { // top level of ROM database
    path = new char[strlen(databasePrefix) + strlen(databaseName) + 1 + strlen(fileName) + 1];
    sprintf(path, "%s%s/%s", databasePrefix, databaseName, fileName);
  } else if (iCluster == -2) { // sensitivity cluster
    path = new char[strlen(databasePrefix) + strlen(databaseName) + 1 + strlen(sensitivityClusterName) + 1 + strlen(fileName) + 1];
    sprintf(path, "%s%s/%s/%s", databasePrefix, databaseName, sensitivityClusterName, fileName);
  } else { // path to appropriate cluster directory
    int addedDigits = 1;
    if (iCluster > 0)  addedDigits = int(ceil(log10(double(iCluster)*10)));
    path = new char[strlen(databasePrefix) + strlen(databaseName) + 1 + strlen(clusterName) + addedDigits + 1 + strlen(fileName) + 1];
    sprintf(path, "%s%s/%s%d/%s", databasePrefix, databaseName, clusterName, iCluster, fileName);
  }  
} 
//----------------------------------------------------------------------------------

template<int dim> 
void NonlinearRom<dim>::createDirectories() { 

  if (com->cpuNum() == 0) {
    char *fullDatabaseName = new char[strlen(databasePrefix) + 1 + strlen(databaseName)];
    sprintf(fullDatabaseName, "%s%s", databasePrefix, databaseName);
    int status;
    status = mkdir(fullDatabaseName, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    for (int iCluster=0; iCluster<nClusters; iCluster++) {
      int addedDigits = 1;
      if (iCluster > 0)  addedDigits = int(ceil(log10(double(iCluster)*10)));
      char *fullClusterName = new char[strlen(fullDatabaseName) + 1 + strlen(clusterName) + addedDigits + 1];
      sprintf(fullClusterName, "%s/%s%d", fullDatabaseName, clusterName, iCluster);
      int status;
      status = mkdir(fullClusterName, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      delete [] fullClusterName;
      fullClusterName = NULL;
    }

    if (strcmp(sensitivityClusterName,"")!=0) {
      char *fullSensitivitiesClusterName = new char[strlen(fullDatabaseName) + 1 + strlen(sensitivityClusterName)  + 1];
      sprintf(fullSensitivitiesClusterName, "%s/%s", fullDatabaseName, sensitivityClusterName);
      int status;
      status = mkdir(fullSensitivitiesClusterName, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      delete [] fullSensitivitiesClusterName;
      fullSensitivitiesClusterName = NULL;
    }

    delete [] fullDatabaseName;
    fullDatabaseName = NULL;
  }
  com->barrier();
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readReferenceState() {

  char* fullRefName;
  const char* refSolName = ioData->input.stateSnapRefSolution;
  fullRefName = new char[strlen(ioData->input.prefix) + strlen(refSolName) + 1];
  sprintf(fullRefName, "%s%s", ioData->input.prefix, refSolName);

  int numSteps = 0;
  double tag = 0.0;
  bool status = domain.readTagFromFile<double, dim>(fullRefName, 0, &tag, &numSteps);
  if (status) {
    if (snapRefState) delete snapRefState;
    snapRefState = new DistSVec<double, dim>(domain.getNodeDistInfo());
    com->fprintf(stdout, "\nReading reference solution for snapshots from %s\n", fullRefName);
  } else {
    fprintf(stderr, "\n*** Error: no snapshots found in \"%s\"\n", fullRefName);
    exit(-1);
  }

  // different stencils for different time integrators
  int refVecIndex = 0;
  if ((ioData->ts.implicit.type == ImplicitData::THREE_POINT_BDF) && (numSteps>=2)) {
    refVecIndex = 1;
  } else if ((ioData->ts.implicit.type == ImplicitData::FOUR_POINT_BDF) && (numSteps>=3)) {
    refVecIndex = 2;
  }

  status = domain.readVectorFromFile(fullRefName, refVecIndex, &tag, *snapRefState);
  
  delete [] fullRefName;
  fullRefName = NULL;

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::initializeClusteredOutputs()
{ // This function is called prior to clustering any non-state snapshots (by Model II (online) simulations
  // or, if using snapshot collection method 0, during GNAT preprocessing). It allows the code
  // to accumulate residual and jacaction snapshots over restarts / multiple simulations

  if (strcmp(residualSnapsName,"")!=0) {
    clusterNewtonCount = new int[nClusters];
 
    if (ioData->output.rom.overwriteNonlinearSnaps == ROMOutputData::OVERWRITE_OFF) {
      for (int iCluster = 0; iCluster < nClusters; ++iCluster) {
        double tag = 0.0;
        int numSteps = 0;
        int step = 0;
        char *snapshotsPath = 0;
        determinePath(residualSnapsName, iCluster, snapshotsPath); 
        bool status = domain.readTagFromFile<double, dim>(snapshotsPath, step, &tag, &numSteps);  // if file DNE, returns false, tag=0, and numSteps=0
        clusterNewtonCount[iCluster] = numSteps;
        delete [] snapshotsPath;
        snapshotsPath = NULL;

        if (strcmp(jacActionSnapsName,"")) { // if outputting jac-action snaps, check that num_residuals == num_jac_actions
          determinePath(jacActionSnapsName, iCluster, snapshotsPath);
          tag = 0;
          numSteps = 0;
          bool status = domain.readTagFromFile<double, dim>(snapshotsPath, step, &tag, &numSteps);
          delete [] snapshotsPath;
          snapshotsPath = NULL;
          if (clusterNewtonCount[iCluster] != numSteps) {
            fprintf(stderr, "*** Error: %d residual snapshots found in cluster %d, %d jac-action snapshots found (should match)",
              clusterNewtonCount[iCluster], iCluster, numSteps);
            exit(-1);
          }
        }
      }
    } else {
      for (int iCluster = 0; iCluster < nClusters; ++iCluster) {
        clusterNewtonCount[iCluster] = 0;
        if (!duplicateSnaps && (com->cpuNum()==0)) {
         char *resSnapshotsPath = 0; 
         char *jacActionSnapshotsPath = 0;
         determinePath(residualSnapsName, iCluster, resSnapshotsPath);
         determinePath(jacActionSnapsName, iCluster, jacActionSnapshotsPath);
         remove(resSnapshotsPath);
         remove(jacActionSnapshotsPath);
         delete [] resSnapshotsPath;
         delete [] jacActionSnapshotsPath;
        }
      }
    }
  } else if (strcmp(jacActionSnapsName,"")!=0) {
    fprintf(stderr, "*** Error: Aero-F assumes that if jacAction snapshots are being collected, then residual snapshots are also being collected");
    exit(-1);
  }

  if (strcmp(krylovSnapsName,"")!=0) {
    clusterKrylovCount = new int[nClusters];
    for (int iCluster = 0; iCluster < nClusters; ++iCluster) {
      double tag = 0.0;
      int numSteps = 0;
      int step = 0;
      char *snapshotsPath = 0;
      determinePath(krylovSnapsName, iCluster, snapshotsPath);
      bool status = domain.readTagFromFile<double, dim>(snapshotsPath, step, &tag, &numSteps);  // if file DNE, returns false, tag=0, and numSteps=0
      clusterKrylovCount[iCluster] = numSteps;
      delete [] snapshotsPath;
      snapshotsPath = NULL;
    }  
  } 
}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::writeClusteredBinaryVectors(int iCluster, DistSVec<double,dim> *U1, DistSVec<double,dim> *U2, 
                                                     DistSVec<double,dim> *U3, char* originalSnapshotFile, int originalSnapshotNumber)
{ 
  // For writing PG residual/jacaction/krylov snapshots during online simulation, 
  // and also for writing FOM residual/krylov snapshots during ROM preprocessing (snapshot collection method 0).
  // Note that residuals and krylov snapshots from FOM simulations are originally output in TsOutput.
  if (strcmp(residualSnapsName,"") && U1)  { 
    char *residualSnapsPath = 0;
    determinePath(residualSnapsName, iCluster, residualSnapsPath);
    if (originalSnapshotFile && (!duplicateSnaps)) {
      if (com->cpuNum()==0) {
        FILE* residualSnapsFile = fopen(residualSnapsPath, "at"); //append
        com->fprintf(residualSnapsFile, "%s %d\n", originalSnapshotFile, originalSnapshotNumber);
        fclose(residualSnapsFile);
      }
    } else {
      domain.writeVectorToFile(residualSnapsPath, clusterNewtonCount[iCluster], 0.0, *U1);
      if (!duplicateSnaps && (com->cpuNum()==0)) {
        FILE* residualSnapsFile = fopen(residualSnapsPath, "at"); //append
        com->fprintf(residualSnapsFile, "%s %d\n", residualSnapsPath, clusterNewtonCount[iCluster]);
        fclose(residualSnapsFile);
      }
    }
    delete [] residualSnapsPath;
    residualSnapsPath = NULL;

    if (strcmp(jacActionSnapsName,"") && U2)  {
      char *jacActionSnapsPath = 0;
      determinePath(jacActionSnapsName, iCluster, jacActionSnapsPath);
      domain.writeVectorToFile(jacActionSnapsPath, clusterNewtonCount[iCluster], 0.0, *U2);
      if (!duplicateSnaps && (com->cpuNum()==0)) {
        FILE* jacActionSnapsFile = fopen(jacActionSnapsPath, "at"); //append
        com->fprintf(jacActionSnapsFile, "%s %d\n", jacActionSnapsPath, clusterNewtonCount[iCluster]);
        fclose(jacActionSnapsFile);
      }
      delete [] jacActionSnapsPath;
      jacActionSnapsPath = NULL;
    }
    
    ++(clusterNewtonCount[iCluster]);
  }

  if (strcmp(krylovSnapsName,"") && U3) {
    char *krylovSnapsPath = 0;
    determinePath(krylovSnapsName, iCluster, krylovSnapsPath);
    if (originalSnapshotFile && (!duplicateSnaps)) {
      if (com->cpuNum()==0) {
        FILE* krylovSnapsFile = fopen(krylovSnapsPath, "at"); //append
        com->fprintf(krylovSnapsFile, "%s %d\n", originalSnapshotFile, originalSnapshotNumber); 
        fclose(krylovSnapsFile);
      }
    } else {
      domain.writeVectorToFile(krylovSnapsPath, clusterKrylovCount[iCluster], 0.0, *U3);
      if (!duplicateSnaps && (com->cpuNum()==0)) {
        FILE* krylovSnapsFile = fopen(krylovSnapsPath, "at"); //append
        com->fprintf(krylovSnapsFile, "%s %d\n", krylovSnapsPath, clusterKrylovCount[iCluster]);
        fclose(krylovSnapsFile);
      }
    }
    delete [] krylovSnapsPath;
    krylovSnapsPath = NULL;

    ++(clusterKrylovCount[iCluster]);
  }

  // storing the entire JPhi is no longer supported since it's infeasible in practice

}


//----------------------------------------------------------------------------------

template<int dim> 
void NonlinearRom<dim>::determineNumResJacMat() { 

  // tests for jacMat; sets numResJacMat
  double tag = 0.0;
  int numSteps = 0;
  int tmp = 0;
  char *jacMatPath = 0;
  determinePath(gappyJacActionName, tmp, jacMatPath);  // check cluster 0
  numResJacMat = (domain.template readTagFromFile<double, dim>(jacMatPath, tmp, &tag, &numSteps)) ? 2 : 1;
  delete [] jacMatPath;
  jacMatPath = NULL;

} 

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readSampleNodes(int iCluster, const char* sampledOrFull, bool deleteExistingRestrictionMapping) {
// this is generalized for the case where iCluster = -1, indicating that the union of all sampled nodes is required
  if (iCluster<0) { // union (which isn't stored in the ROM database)
    std::set<int> sampleNodesUnion;
    for (int jCluster=0; jCluster<nClusters; ++jCluster) {
      readClusteredSampleNodes(jCluster, sampledOrFull, true);
      for (std::vector<int>::iterator it = sampleNodes.begin(); it != sampleNodes.end(); ++it) {
        sampleNodesUnion.insert(*it);
      }   
    }
    sampleNodes.clear();
    nSampleNodes = sampleNodesUnion.size();
    sampleNodes.reserve(nSampleNodes);
    for (std::set<int>::iterator it = sampleNodesUnion.begin(); it != sampleNodesUnion.end(); ++it) {
      sampleNodes.push_back(*it);
    }

  } else { // sample nodes for iCluster, which is stored in the ROM database
    readClusteredSampleNodes(iCluster, sampledOrFull, deleteExistingRestrictionMapping);
  }
}


//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredSampleNodes(int iCluster, const char* sampledOrFull, bool deleteExistingRestrictionMapping) {


  if (storedAllOnlineQuantities) {
    com->fprintf(stdout, " ... loading sampled nodes for cluster %d\n", iCluster);
    sampleNodes.clear();
    sampleNodes = *(allSampleNodes[iCluster]);
    nSampleNodes = sampleNodes.size();
    restrictionMapping = allRestrictionMappings[iCluster];  //setting pointer
    restrictionMapping->recomputeConnectedTopology();       //need to reinitialize the domain for this restriction mapping
    return;
  }

  char *sampleNodesPath = 0;
  if (strcmp(sampledOrFull,"sampled")==0) {
    determinePath(sampledNodesName, iCluster, sampleNodesPath);
  } else if (strcmp(sampledOrFull,"full")==0) {
    determinePath(sampledNodesFullCoordsName, iCluster, sampleNodesPath);
  } else {
     com->fprintf(stderr, "*** Error: unexpected value for sampledOrFull (%s)\n", sampledOrFull);
     exit (-1);
  }
  FILE *sampleNodeFile = fopen(sampleNodesPath, "r");
  if (!sampleNodeFile)  {
     com->fprintf(stderr, "*** Error: unable to open file %s\n", sampleNodesPath);
     exit (-1);
  }

  int _n;

  _n = fscanf(sampleNodeFile, "%d", &nSampleNodes);  // first entry is the number of sample nodes

  sampleNodes.clear();
  sampleNodes.reserve(nSampleNodes);  // know it will be nSampleNodes long (efficiency)
  int index, currentSampleNode;
  for (int i = 0; i < nSampleNodes; ++i){
    _n = fscanf(sampleNodeFile, "%d", &index);
    _n = fscanf(sampleNodeFile, "%d", &currentSampleNode);
    sampleNodes.push_back(currentSampleNode-1); // reads in the sample node plus one (TODO - change this after code validation)
    if (_n != 1) {
      fprintf(stderr, "*** Error: unexpected file format encountered while reading \"%s\" in readClusteredSampleNodes()\n", sampleNodesPath);
      exit (-1);
    }
  }

  delete [] sampleNodesPath;

  fclose(sampleNodeFile);

  if (ioData->problem.type[ProblemData::NLROMOFFLINE]) {
    // continue
  } else {
    if (restrictionMapping && deleteExistingRestrictionMapping) delete restrictionMapping;
    restrictionMapping = new RestrictionMapping<dim>(&domain, sampleNodes.begin(), sampleNodes.end());
  }
}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::deleteRestrictedQuantities() {

  if (resMat) {
    delete resMat;
    resMat = NULL;
  }

  if (jacMat) {
    delete jacMat;
    jacMat = NULL;
  }

  if (metric) {
    delete metric;
    metric = NULL;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredNonlinearMetric(int iCluster) {

  if (storedAllOnlineQuantities || storedAllOfflineQuantities) {
      com->fprintf(stdout, " ... loading nonlinear approximate metric for cluster %d\n", iCluster);
      if (metric) delete metric;
      int numVecs = (allMetrics[iCluster])->numVectors();
      metric = new VecSet<DistSVec<double, dim> >(numVecs, restrictionMapping->restrictedDistInfo());
      for (int iVec=0; iVec<numVecs; ++iVec) 
        (*metric)[iVec] = (*(allMetrics[iCluster]))[iVec];
    return;
  }

  char* metricPath = 0;
  if (metric) {
    delete metric;
    metric = NULL;
  }
  determinePath(approxMetricNonlinearName, iCluster, metricPath);
  com->fprintf(stdout, "\nReading nonlinear approximate metric for cluster %d\n", iCluster);

  double tag = 0.0;
  int numVecs = 0;
  int step = 0;
  bool status = domain.readTagFromFile<double, dim>(metricPath, step, &tag, &numVecs);  // if file DNE, returns false, tag=0, and numSteps=0

  if (!status) {
    fprintf(stderr, "\nCould not open file \"%s\" in readClusteredNonlinearMetric()\n", metricPath);
    exit(-1);
  }

  VecSet<DistSVec<double, dim> > fullDistInfoMetric(numVecs, domain.getNodeDistInfo());

  for (int iVec=0; iVec<numVecs; ++iVec) {
    status = domain.readVectorFromFile(metricPath, iVec, &tag, fullDistInfoMetric[iVec]);
  }

  delete [] metricPath;

  // restrictionMapping is set during readClusteredSampleNodes
  metric = new VecSet<DistSVec<double, dim> >(numVecs, restrictionMapping->restrictedDistInfo());

  for (int i = 0; i < numVecs; ++i) {
    restrictionMapping->restriction(fullDistInfoMetric[i],(*metric)[i]);
  }

}



//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredGappyMatrix(int iCluster, const char* matrixType) {

  if (storedAllOnlineQuantities || storedAllOfflineQuantities) {
    if (strcmp(matrixType,"resMatrix")==0) { 
      com->fprintf(stdout, " ... loading gappy residual matrix for cluster %d\n", iCluster);
      if (resMat) delete resMat;
      int numVecs = (allResMat[iCluster])->numVectors();
      resMat = new VecSet<DistSVec<double, dim> >(numVecs, restrictionMapping->restrictedDistInfo());
      for (int iVec=0; iVec<numVecs; ++iVec) 
        (*resMat)[iVec] = (*(allResMat[iCluster]))[iVec];
    } else if (strcmp(matrixType,"jacMatrix")==0) {
      com->fprintf(stdout, " ... loading gappy jacAction matrix for cluster %d\n", iCluster);
      if (jacMat) delete jacMat;       
      int numVecs = (allJacMat[iCluster])->numVectors();
      jacMat = new VecSet<DistSVec<double, dim> >(numVecs, restrictionMapping->restrictedDistInfo());
      for (int iVec=0; iVec<numVecs; ++iVec) 
        (*jacMat)[iVec] = (*(allJacMat[iCluster]))[iVec]; 
    }
    return;
  }

  char* matrixPath = 0;

  if (strcmp(matrixType,"resMatrix")==0) {
      if (resMat) {
        delete resMat;
        resMat = NULL;
      }
      determinePath(gappyResidualName, iCluster, matrixPath);
      com->fprintf(stdout, "\nReading gappy residual matrix for cluster %d\n", iCluster);
  } else if (strcmp(matrixType,"jacMatrix")==0) {
      if (jacMat) {
        delete jacMat;
        jacMat = NULL;
      }
      determinePath(gappyJacActionName, iCluster, matrixPath);
      com->fprintf(stdout, "\nReading gappy jacAction matrix for cluster %d\n", iCluster);
  } else {
      fprintf(stderr, "*** Error: unrecognized matrix type (%s) in readClusteredGappyMatrix()\n", matrixType);
      exit (-1);
  }

  double tag = 0.0;
  int numVecs = 0;
  int step = 0;
  bool status = domain.readTagFromFile<double, dim>(matrixPath, step, &tag, &numVecs);  // if file DNE, returns false, tag=0, and numSteps=0

  if (!status) {
    fprintf(stderr, "\nCould not open file \"%s\" in readClusteredGappyMatrix()\n", matrixPath);
    exit(-1);
  }

  VecSet<DistSVec<double, dim> > gappyMatrix(numVecs, domain.getNodeDistInfo());

  for (int iVec=0; iVec<numVecs; ++iVec) {
    status = domain.readVectorFromFile(matrixPath, iVec, &tag, gappyMatrix[iVec]);
  }

  delete [] matrixPath;

  // restrictionMapping is set during readClusteredSampleNodes
  VecSet<DistSVec<double, dim> >* restrictedGappyMatrix = new VecSet<DistSVec<double, dim> >(numVecs, restrictionMapping->restrictedDistInfo());

  for (int i = 0; i < numVecs; ++i) {
    restrictionMapping->restriction(gappyMatrix[i],(*restrictedGappyMatrix)[i]);
  }

  if (strcmp(matrixType,"resMatrix")==0) {
      resMat = restrictedGappyMatrix;
  } else if (strcmp(matrixType,"jacMatrix")==0) {
      jacMat = restrictedGappyMatrix;
  }

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readApproxMetricStateLowRankFactor(const char* sampledOrFull) {

  char *approxMetricPath;

  if ((ioData->problem.alltype == ProblemData::_STEADY_NONLINEAR_ROM_POST_
       || ioData->problem.alltype == ProblemData::_UNSTEADY_NONLINEAR_ROM_POST_)
       && strcmp(approxMetricStateLowRankSurfaceCoordsName,"")!=0) {
      determinePath(approxMetricStateLowRankSurfaceCoordsName,-1,approxMetricPath);
    if (!clusterCenters) readClusterCenters(""); // surface centers
  } else if (strcmp(sampledOrFull,"sampled")==0) {
    determinePath(approxMetricStateLowRankName,-1,approxMetricPath);
    if (!clusterCenters) readClusterCenters("sampledCenters"); 
  } else if (strcmp(sampledOrFull,"full")==0) {
    determinePath(approxMetricStateLowRankFullCoordsName,-1,approxMetricPath);
    if (!clusterCenters) readClusterCenters("centers");
  } else {
    fprintf(stderr, "*** Error: sampledOrFull not specified correctly in readApproxMetricStateLowRankFactor()\n");
    exit (-1);
  }
  
  int nRank = 0; 
  int dummyStep = 0;
  double dummyTag = 0.0;
  bool status = domain.readTagFromFile<double, dim>(approxMetricPath,dummyStep,&dummyTag,&nRank); 

  if (!status) {
    fprintf(stderr, "\nCould not open file \"%s\" in readApproxMetricStateLowRankFactor()\n", approxMetricPath);
    exit(-1);
  }
  nLowRankFactors = nRank;

  if (lowRankFactor) delete lowRankFactor;
  if (cForFastDistComp) {
    for (int iCluster = 0; iCluster < nClusters; ++iCluster) {
      for (int jCluster = 0; jCluster < nClusters; ++jCluster) {
        delete [] hForFastDistComp[iCluster][jCluster];
        delete [] cForFastDistComp[iCluster][jCluster];
      }
      delete [] hForFastDistComp[iCluster];
      delete [] cForFastDistComp[iCluster];
    }
    delete [] hForFastDistComp;
    delete [] cForFastDistComp;
  }

  lowRankFactor = new VecSet< DistSVec<double, dim> >(nLowRankFactors, domain.getNodeDistInfo());
  cForFastDistComp = new double**[nClusters];
  hForFastDistComp = new double**[nClusters];
  for (int iCluster = 0; iCluster < nClusters; ++iCluster) {
    hForFastDistComp[iCluster] = new double*[nClusters];
    cForFastDistComp[iCluster] = new double*[nClusters];
    for (int jCluster = 0; jCluster < nClusters; ++jCluster) {
      hForFastDistComp[iCluster][jCluster] = new double[1];
      cForFastDistComp[iCluster][jCluster] = new double[nLowRankFactors];
    }
  } 
  com->fprintf(stdout, "\nReading approximated metric low rank factor\n");
  double tmp;

  for (int iRank = 0; iRank < nLowRankFactors; ++iRank) { 
    status = domain.readVectorFromFile(approxMetricPath, iRank, &tmp, (*lowRankFactor)[iRank]);
    if (!status) {
      fprintf(stderr, "\nError reading the low rank vector #%d from file \"%s\"\n", iRank,approxMetricPath);
      exit(-1);
    }
  }

  // build c
  double **approxMetricMaskCenters = new double*[nClusters];
  for (int iCluster = 0; iCluster < nClusters; ++iCluster)
    approxMetricMaskCenters[iCluster] = new double[nLowRankFactors];

  for (int iCluster = 0; iCluster < nClusters; ++iCluster) {
    for (int iRank = 0; iRank < nLowRankFactors; ++iRank) {
      approxMetricMaskCenters[iCluster][iRank] = (*lowRankFactor)[iRank] * (*clusterCenters)[iCluster];    
    }
  }

  for (int iCluster = 0; iCluster < nClusters; ++iCluster) {
    for (int jCluster = 0; jCluster < nClusters; ++jCluster) {
      for (int iRank = 0; iRank < nLowRankFactors; ++iRank) {
        cForFastDistComp[iCluster][jCluster][iRank] = 2.0*(approxMetricMaskCenters[jCluster][iRank]-approxMetricMaskCenters[iCluster][iRank]);
      }
    }
  }
 
  com->barrier();
  delete [] approxMetricPath;
  approxMetricPath = NULL;  

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readAllClusteredOnlineQuantities() {

// stores all online quantities at the beginning of the simulation
// (as opposed to reading information at each cluster switch)

// note: fast distance calculation info is handled separately

// initial allocation

  allStateBases = new VecSet< DistSVec<double, dim> >*[nClusters];
  allStateSVals = new std::vector<double>*[nClusters];

  allNBuffer.resize(nClusters);

  if (ioData->romOnline.krylov.include==NonlinearRomOnlineNonStateData::INCLUDE_ON) {
    allKrylovBases = new VecSet< DistSVec<double, dim> >*[nClusters];
    allKrylovSVals = new std::vector<double>*[nClusters];
  }

  if (ioData->romOnline.basisUpdates!=NonlinearRomOnlineData::UPDATES_OFF ||
      ioData->romOnline.projectSwitchStateOntoAffineSubspace!=NonlinearRomOnlineData::PROJECT_OFF) {
    allRefStates = new VecSet< DistSVec<double, dim> >(nClusters, domain.getNodeDistInfo());
    if (ioData->romOnline.basisUpdates!=NonlinearRomOnlineData::UPDATES_OFF) {
      allColumnSumsV = new std::vector<double>*[nClusters];
    }
    if (false) { 
      // additional update info
    }
  }

  if (ioData->romOnline.systemApproximation == NonlinearRomOnlineData::GNAT){
    allSampleNodes = new std::vector<int>*[nClusters];
    allRestrictionMappings = new RestrictionMapping<dim>*[nClusters];
    allResMat = new VecSet< DistSVec<double, dim> >*[nClusters];
    if (numResJacMat==2) allJacMat = new VecSet< DistSVec<double, dim> >*[nClusters];
    allMetrics = NULL;
  } else if (ioData->romOnline.systemApproximation == NonlinearRomOnlineData::COLLOCATION ) {
    allSampleNodes = new std::vector<int>*[nClusters];
    allRestrictionMappings = new RestrictionMapping<dim>*[nClusters];
    allResMat = NULL;
    allJacMat = NULL;
    allMetrics = NULL;
  } else if (ioData->romOnline.systemApproximation == NonlinearRomOnlineData::APPROX_METRIC_NL ) {
    allSampleNodes = new std::vector<int>*[nClusters];
    allRestrictionMappings = new RestrictionMapping<dim>*[nClusters];
    allResMat = NULL;
    allJacMat = NULL;
    allMetrics = new VecSet< DistSVec<double, dim> >*[nClusters];
  }

// read and store online info for each cluster

  if (ioData->romOnline.systemApproximation == NonlinearRomOnlineData::SYSTEM_APPROXIMATION_NONE) {
      for (int iCluster=0; iCluster<nClusters; ++iCluster) {

        // read state ROB and sVals
        readClusteredBasis(iCluster, "state");
        allStateBases[iCluster] = new VecSet< DistSVec<double, dim> >(basis->numVectors(), domain.getNodeDistInfo());
        for (int iVec=0; iVec<basis->numVectors(); ++iVec)
          (*(allStateBases[iCluster]))[iVec] = (*basis)[iVec];
        delete basis;
        basis = NULL;
        allStateSVals[iCluster] = new vector<double>;
        *(allStateSVals[iCluster]) = *sVals;
        delete sVals;
        sVals = NULL;
        allNBuffer[iCluster]=nBuffer;
  
        // read update info
        if (ioData->romOnline.basisUpdates!=NonlinearRomOnlineData::UPDATES_OFF || 
            ioData->romOnline.projectSwitchStateOntoAffineSubspace!=NonlinearRomOnlineData::PROJECT_OFF) {
          readClusteredUpdateInfo(iCluster, "state");
          (*allRefStates)[iCluster] = *Uref;
          delete Uref;
          Uref = NULL;
          if (ioData->romOnline.basisUpdates!=NonlinearRomOnlineData::UPDATES_OFF) {
            allColumnSumsV[iCluster] = new vector<double>;
            *(allColumnSumsV[iCluster]) = *columnSumsV;
            delete columnSumsV;
            columnSumsV = NULL;
            if (false) { 
              // read additional update information
            }
          }
        }

        // read Krylov ROB and sVals
        if (ioData->romOnline.krylov.include==NonlinearRomOnlineNonStateData::INCLUDE_ON) {
          readClusteredBasis(iCluster, "krylov");
          allKrylovBases[iCluster] = new VecSet< DistSVec<double, dim> >(basis->numVectors(), domain.getNodeDistInfo());
          for (int iVec=0; iVec<basis->numVectors(); ++iVec)
            (*(allKrylovBases[iCluster]))[iVec] = (*basis)[iVec];
          delete basis;
          basis = NULL;
          allKrylovSVals[iCluster] = new vector<double>;
          *(allKrylovSVals[iCluster]) = *sVals;
          delete sVals;
          sVals = NULL;
        }

      }

      // read sensitivity ROB and sVals
      if (ioData->romOnline.sensitivity.include==NonlinearRomOnlineNonStateData::INCLUDE_ON) {
        readClusteredBasis(-2, "sensitivity");
        sensitivityBasis = new VecSet< DistSVec<double, dim> >(basis->numVectors(), domain.getNodeDistInfo());
        for (int iVec=0; iVec<basis->numVectors(); ++iVec)
          (*sensitivityBasis)[iVec] = (*basis)[iVec];
        delete basis;
        basis = NULL;
        sensitivitySVals = new vector<double>;
        *sensitivitySVals = *sVals;
        delete sVals;
        sVals = NULL;
      }
            
    } else if (   ioData->romOnline.systemApproximation == NonlinearRomOnlineData::GNAT 
               || ioData->romOnline.systemApproximation == NonlinearRomOnlineData::COLLOCATION 
               || ioData->romOnline.systemApproximation == NonlinearRomOnlineData::APPROX_METRIC_NL) {
      for (int iCluster=0; iCluster<nClusters; ++iCluster) {
        // read sample nodes
        readClusteredSampleNodes(iCluster, "sampled", false); // resets restriction map
        allSampleNodes[iCluster] = new vector<int>;
        *(allSampleNodes[iCluster]) = sampleNodes;
        allRestrictionMappings[iCluster] = restrictionMapping; // sets pointer to dynamically allocated memory
         
        if (ioData->romOnline.systemApproximation == NonlinearRomOnlineData::GNAT) {
          // read gappy POD matrix for residual
          readClusteredGappyMatrix(iCluster, "resMatrix");
          allResMat[iCluster] = new VecSet< DistSVec<double, dim> >(resMat->numVectors(), restrictionMapping->restrictedDistInfo());
          for (int iVec=0; iVec<resMat->numVectors(); ++iVec)
            (*(allResMat[iCluster]))[iVec] = (*resMat)[iVec];
          delete resMat;
          resMat = NULL;

          // read gappy POD matrix for jacobian
          if (numResJacMat==2) {
            readClusteredGappyMatrix(iCluster, "jacMatrix");
            allJacMat[iCluster] = new VecSet< DistSVec<double, dim> >(jacMat->numVectors(), restrictionMapping->restrictedDistInfo());
            for (int iVec=0; iVec<jacMat->numVectors(); ++iVec)
              (*(allJacMat[iCluster]))[iVec] = (*jacMat)[iVec];
            delete jacMat;
            jacMat = NULL;
          }
        }

        if (ioData->romOnline.systemApproximation == NonlinearRomOnlineData::APPROX_METRIC_NL) {
          // read metric
          readClusteredNonlinearMetric(iCluster);
          allMetrics[iCluster] = new VecSet< DistSVec<double, dim> >(metric->numVectors(), restrictionMapping->restrictedDistInfo());
          for (int iVec=0; iVec<metric->numVectors(); ++iVec)
            (*(allMetrics[iCluster]))[iVec] = (*metric)[iVec];
          delete metric;
          metric = NULL;
        }

        // read sampled state ROB and sVals
        readClusteredBasis(iCluster, "sampledState");
        allStateBases[iCluster] = new VecSet< DistSVec<double, dim> >(basis->numVectors(), domain.getNodeDistInfo());
        for (int iVec=0; iVec<basis->numVectors(); ++iVec)
          (*(allStateBases[iCluster]))[iVec] = (*basis)[iVec];
        delete basis;
        basis = NULL;
        allStateSVals[iCluster] = new vector<double>;
        *(allStateSVals[iCluster]) = *sVals;
        delete sVals;
        sVals = NULL;
        allNBuffer[iCluster]=nBuffer;

        // read ROB update information
        if (ioData->romOnline.basisUpdates!=NonlinearRomOnlineData::UPDATES_OFF ||
            ioData->romOnline.projectSwitchStateOntoAffineSubspace!=NonlinearRomOnlineData::PROJECT_OFF) {
          readClusteredUpdateInfo(iCluster, "sampledState");
          (*allRefStates)[iCluster] = *Uref;
          delete Uref;
          Uref = NULL;
          if (ioData->romOnline.basisUpdates!=NonlinearRomOnlineData::UPDATES_OFF) {
            allColumnSumsV[iCluster] = new vector<double>;
            *(allColumnSumsV[iCluster]) = *columnSumsV;
            delete columnSumsV;
            columnSumsV = NULL;
            if (false) {
              // read additional update information
            }
          }
        }

        // read sampled Krylov ROB
        if (ioData->romOnline.krylov.include==NonlinearRomOnlineNonStateData::INCLUDE_ON) { 
          readClusteredBasis(iCluster, "sampledKrylov");
          allKrylovBases[iCluster] = new VecSet< DistSVec<double, dim> >(basis->numVectors(), domain.getNodeDistInfo());
          for (int iVec=0; iVec<basis->numVectors(); ++iVec)
            (*(allKrylovBases[iCluster]))[iVec] = (*basis)[iVec];
          delete basis;
          basis = NULL;
          allKrylovSVals[iCluster] = new vector<double>;
          *(allKrylovSVals[iCluster]) = *sVals;
          delete sVals;
          sVals = NULL;
        }
 
      }

      // read sampled sensitivity ROB
      if (ioData->romOnline.sensitivity.include==NonlinearRomOnlineNonStateData::INCLUDE_ON) {
        readClusteredBasis(-2, "sampledSensitivity");
        sensitivityBasis = new VecSet< DistSVec<double, dim> >(basis->numVectors(), domain.getNodeDistInfo());
        for (int iVec=0; iVec<basis->numVectors(); ++iVec)
          (*sensitivityBasis)[iVec] = (*basis)[iVec];
        delete basis;
        basis = NULL;
        sensitivitySVals = new vector<double>;
        *sensitivitySVals = *sVals;
        delete sVals;
        sVals = NULL;
      }
  } else {
      fprintf(stderr, "*** Error:  Unexpected system approximation type in readAllClusteredOnlineQuantities()\n");
      exit(-1);
  }

  storedAllOnlineQuantities = true;

}


//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readAllClusteredOfflineQuantities() {

// stores all offline quantities after finishing POD

// initial allocation

  allStateBases = new VecSet< DistSVec<double, dim> >*[nClusters];
  allStateSVals = new std::vector<double>*[nClusters];

  /*if (ioData->romOnline.krylov.include==NonlinearRomOnlineNonStateData::INCLUDE_ON) {
    allKrylovBases = new VecSet< DistSVec<double, dim> >*[nClusters];
    allKrylovSVals = new std::vector<double>*[nClusters];
  }*/

  allRefStates = new VecSet< DistSVec<double, dim> >(nClusters, domain.getNodeDistInfo());
  allColumnSumsV = new std::vector<double>*[nClusters];

// read and store online info for each cluster
  for (int iCluster=0; iCluster<nClusters; ++iCluster) {

    // read state ROB and sVals
    readClusteredBasis(iCluster, "state");
    allStateBases[iCluster] = new VecSet< DistSVec<double, dim> >(basis->numVectors(), domain.getNodeDistInfo());
    for (int iVec=0; iVec<basis->numVectors(); ++iVec)
      (*(allStateBases[iCluster]))[iVec] = (*basis)[iVec];
    delete basis;
    basis = NULL;
    allStateSVals[iCluster] = new vector<double>;
    *(allStateSVals[iCluster]) = *sVals;
    delete sVals;
    sVals = NULL;

    // read update info
    readClusteredReferenceState(iCluster, "state");
    (*allRefStates)[iCluster] = *Uref;
    delete Uref;
    Uref = NULL;
    readClusteredColumnSumsV(iCluster, "state");
    allColumnSumsV[iCluster] = new vector<double>;
    *(allColumnSumsV[iCluster]) = *columnSumsV;
    delete columnSumsV;
    columnSumsV = NULL;

    // read Krylov ROB and sVals
   /* if (ioData->romOnline.krylov.include==NonlinearRomOnlineNonStateData::INCLUDE_ON) {
      readClusteredBasis(iCluster, "krylov");
      allKrylovBases[iCluster] = new VecSet< DistSVec<double, dim> >(basis->numVectors(), domain.getNodeDistInfo());
      for (int iVec=0; iVec<basis->numVectors(); ++iVec)
        (*(allKrylovBases[iCluster]))[iVec] = (*basis)[iVec];
      delete basis;
      basis = NULL;
      allKrylovSVals[iCluster] = new vector<double>;
      *(allKrylovSVals[iCluster]) = *sVals;
      delete sVals;
      sVals = NULL;
    } */

  }
 
  /* // read sensitivity ROB and sVals
  if (ioData->romOnline.sensitivity.include==NonlinearRomOnlineNonStateData::INCLUDE_ON) {
    readClusteredBasis(-2, "sensitivity");
    sensitivityBasis = new VecSet< DistSVec<double, dim> >(basis->numVectors(), domain.getNodeDistInfo());
    for (int iVec=0; iVec<basis->numVectors(); ++iVec)
      (*sensitivityBasis)[iVec] = (*basis)[iVec];
    delete basis;
    basis = NULL;
    sensitivitySVals = new vector<double>;
    *sensitivitySVals = *sVals;
    delete sVals;
    sVals = NULL;
  }*/

  storedAllOfflineQuantities = true;

}


//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::writeReducedCoords(const int totalTimeSteps, bool clusterSwitch, bool update, int iCluster, Vec<double> dUromTimeIt) {

  int nPod = basis->numVectors();

  if (com->cpuNum() == 0) {
    if (clustUsageFile)  // for plotting only
      com->fprintf(clustUsageFile,"%d %d %d %d %d %d\n", totalTimeSteps, iCluster, nPod, nState, nKrylov, nSens);
    if (reducedCoordsFile) {
      if (clusterSwitch) {
        if (update) {
          com->fprintf(reducedCoordsFile,"%d switch update %d %d %d %d %d\n",
            totalTimeSteps, iCluster, nPod, nState, nKrylov, nSens);
        } else {
          com->fprintf(reducedCoordsFile,"%d switch noUpdate %d %d %d %d %d\n",
            totalTimeSteps, iCluster, nPod, nState, nKrylov, nSens);
        }
      } else if (update) {
        com->fprintf(reducedCoordsFile,"%d noSwitch update %d %d %d %d %d\n",
          totalTimeSteps, iCluster, nPod, nState, nKrylov, nSens);
      } else { 
        com->fprintf(reducedCoordsFile,"%d noSwitch noUpdate %d %d %d %d %d\n",
          totalTimeSteps, iCluster, nPod, nState, nKrylov, nSens);
      }
      for (int iPod=0; iPod<nPod; ++iPod) {
        com->fprintf(reducedCoordsFile, "%23.15e\n", dUromTimeIt[iPod]);
      }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::outputCenterNorms(std::vector<std::vector<double> > &vec) {
// note: actually norm squared

  char *infoPath = 0;
  determinePath(centerNormsName, -1, infoPath);   

  int nFluidMeshNodes = domain.getNumGlobNode();    

  if (com->cpuNum() == 0) {
    // format:  numFullMeshNodes numClusters vecSize; element1_1; element1_2; ...
    FILE *outputFile = fopen(infoPath, "wt");

    com->fprintf(outputFile, "%d %d %d\n", nFluidMeshNodes, nClusters, vec[0].size());

    for (int iCenter=0; iCenter<nClusters; ++iCenter) {
      for (int iVec=0; iVec<vec[iCenter].size(); ++iVec) {
        com->fprintf(outputFile,"%23.15e\n", vec[iCenter][iVec]);
      }
    }
 
    fclose (outputFile);
  }

  delete [] infoPath;
  infoPath = NULL;

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readCenterNorms() {
// note: actually norm squared

  centerNorms.clear();

  char *infoPath = 0;
  determinePath(centerNormsName, -1, infoPath);   

  // format:  numFullMeshNodes numClusters vecSize; element1_1; element1_2; ...
  FILE *inputFile = fopen(infoPath, "r");
  int _n;

  int expectedNClusters;
  int vecSize;

  _n = fscanf(inputFile, "%d %d %d\n", &nFullMeshNodes, &expectedNClusters, &vecSize);

  assert(expectedNClusters == nClusters);

  checkInitialConditionScenario();
  if (specifiedIC) {
    assert(vecSize==1);
  } else if (interpolatedMultiIC) {
    assert(vecSize==1);
  } else {
    assert(vecSize==(dim+1));
  }

  centerNorms.resize(nClusters);
  for (int iCenter=0; iCenter<nClusters; ++iCenter)
    centerNorms[iCenter].reserve(vecSize);

  double tmpVal;

  for (int iCenter=0; iCenter<nClusters; ++iCenter) {
    for (int iVec=0; iVec<vecSize; ++iVec) {
      _n = fscanf(inputFile,"%le", &tmpVal);
      if (_n == 1) {
        centerNorms[iCenter].push_back(tmpVal);
      } else if (feof(inputFile)) {
        break;
      } else {
        fprintf(stderr, "*** Error: fscanf of centerNorms file interrupted by non-EOF error\n");
        exit(-1);
      }
    }
  }

  fclose (inputFile);

  delete [] infoPath;
  infoPath = NULL;

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::outputClusteredInfoASCII(int iCluster, const char* type, std::vector<double>* vec1,
                                             std::vector<std::vector<double> >* vec2,
                                             std::vector<std::vector<std::vector<double> > >* vec3,
                                             std::vector<std::vector<std::vector<std::vector<double> > > >* vec4) {
                                                           
  // interface for outputting small ASCII quantities

  // sanity check inputs
  if ( (vec1&&(vec2||vec3||vec4)) || (vec2&&(vec3||vec4)) || (vec3&&vec4)) {
    fprintf(stderr, "*** Error: too many containers provided to outputClusteredInfoAscii()\n");
    exit(-1);
  } else if ((!vec1)&&(!vec2)&&(!vec3)&&(!vec4)) {
    fprintf(stderr, "*** Error: no container provided to outputClusteredInfoAscii()\n");
    exit(-1); 
  }

  //TODO INTERP

  char *infoPath = NULL;
  if (strcmp(type, "referenceState") == 0) { // 2*(U_center_p - U_center_m)^T U_ref
    determinePath(stateDistanceComparisonInfoExactUpdatesName, iCluster, infoPath);
    assert(vec2);
  } else if (strcmp(type, "initialCondition") == 0) { // 2*(U_center_p - U_center_m)^T U_ic
    determinePath(stateDistanceComparisonInfoExactUpdatesName, -1, infoPath);
    assert(vec2);
  } else if (strcmp(type, "multiInitialCondition") == 0) { // 2*(U_center_p - U_center_m)^T U_ic
    determinePath(stateDistanceComparisonInfoExactUpdatesMultiICName, -1, infoPath);
    assert(vec3);
  } else if (strcmp(type, "state") == 0) { // 2*(U_center_p - U_center_m)^T V_state_k
    determinePath(stateDistanceComparisonInfoName, iCluster, infoPath);
    assert(vec3);
  } else if (strcmp(type, "krylov") == 0) { // 2*(U_center_p - U_center_m)^T V_krylov_k
    determinePath(krylovDistanceComparisonInfoName, iCluster, infoPath);
    assert(vec3);
  } else if (strcmp(type, "sensitivity") == 0) { // 2*(U_center_p - U_center_m)^T V_sens
    determinePath(sensitivityDistanceComparisonInfoName, -2, infoPath);
    assert(vec3);
  } else if (strcmp(type, "distanceMatrix") == 0) { // A_ij = ||U_i - U_j||_2 
    determinePath(distanceMatrixName, -1, infoPath);
    assert(vec2);
  } else if (strcmp(type, "basisNormalizedCenterProducts") == 0) { // (rob_i^T * Center_p / || Center_p ||_2)
    determinePath(basisNormalizedCenterProductsName, -1, infoPath);
    assert(vec3);
  } else if (strcmp(type, "basisBasisProducts") == 0) {               
    determinePath(basisBasisProductsName, -1, infoPath);
    assert(vec4);
  } else if (strcmp(type, "basisUrefProducts") == 0) {               
    determinePath(basisUrefProductsName, -1, infoPath);
    assert(vec3);
  } else if (strcmp(type, "urefUrefProducts") == 0) {                
    determinePath(urefUrefProductsName, -1, infoPath);
    assert(vec2);
  } else if (strcmp(type, "urefUicProducts") == 0) {               
    determinePath(urefUicProductsName, -1, infoPath);
    assert(vec1);
  } else if (strcmp(type, "urefMultiUicProducts") == 0) {
    determinePath(urefMultiUicProductsName, -1, infoPath);
    assert(vec2);
  } else if (strcmp(type, "basisUicProducts") == 0) {               
    determinePath(basisUicProductsName, -1, infoPath);
    assert(vec2);
  } else if (strcmp(type, "basisMultiUicProducts") == 0) {
    determinePath(basisMultiUicProductsName, iCluster, infoPath);
    assert(vec3);
  } else if (strcmp(type, "urefComponentwiseSums") == 0) {               
    determinePath(urefComponentwiseSumsName, -1, infoPath);
    assert(vec2);
  } else if (strcmp(type, "basisComponentwiseSums") == 0) {               
    determinePath(basisComponentwiseSumsName, -1, infoPath);
    assert(vec3);
  } else if (strcmp(type, "correlationMatrix") == 0) {
    determinePath(correlationMatrixName, iCluster, infoPath);
    assert(vec2);
  } else if (strcmp(type, "multiUicMultiUicProducts") == 0) {                
    determinePath(multiUicMultiUicProductsName, -1, infoPath);
    assert(vec2);
  } else {
    fprintf(stderr, "*** Error: unexpected TYPE in outputClusteredInfoAscii()\n");
    exit(-1);
  }
  
  writeMultiVecASCII(infoPath, vec1, vec2, vec3, vec4);

  delete [] infoPath;
  infoPath = NULL;

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readClusteredInfoASCII(int iCluster, const char* type, std::vector<double>* vec1,
                                             std::vector<std::vector<double> >* vec2,
                                             std::vector<std::vector<std::vector<double> > >* vec3,
                                             std::vector<std::vector<std::vector<std::vector<double> > > >* vec4) {

  // interface for reading small ASCII quantities
                                    
  // sanity check inputs
  if ( (vec1&&(vec2||vec3||vec4)) || (vec2&&(vec3||vec4)) || (vec3&&vec4)) {
    fprintf(stderr, "*** Error: too many containers provided to readClusteredInfoAscii()\n");
    exit(-1);
  } else if ((!vec1)&&(!vec2)&&(!vec3)&&(!vec4)) {
    fprintf(stderr, "*** Error: no container provided to readClusteredInfoAscii()\n");
    exit(-1); 
  }

  char *infoPath = NULL;

  if (strcmp(type, "referenceState") == 0) { // 2*(U_center_p - U_center_m)^T U_ref
    determinePath(stateDistanceComparisonInfoExactUpdatesName, iCluster, infoPath);
    assert(vec2);
  } else if (strcmp(type, "initialCondition") == 0) { // 2*(U_center_p - U_center_m)^T U_ic
    determinePath(stateDistanceComparisonInfoExactUpdatesName, -1, infoPath);
    assert(vec2);
  } else if (strcmp(type, "multiInitialCondition") == 0) { // 2*(U_center_p - U_center_m)^T U_ic
    determinePath(stateDistanceComparisonInfoExactUpdatesMultiICName, -1, infoPath);
    assert(vec3);
  } else if (strcmp(type, "state") == 0) { // 2*(U_center_p - U_center_m)^T V_state_k
    determinePath(stateDistanceComparisonInfoName, iCluster, infoPath);
    assert(vec3);
  } else if (strcmp(type, "krylov") == 0) { // 2*(U_center_p - U_center_m)^T V_krylov_k
    determinePath(krylovDistanceComparisonInfoName, iCluster, infoPath);
    assert(vec3);
  } else if (strcmp(type, "sensitivity") == 0) { // 2*(U_center_p - U_center_m)^T V_sens
    determinePath(sensitivityDistanceComparisonInfoName, -2, infoPath);
    assert(vec3);
  } else if (strcmp(type, "distanceMatrix") == 0) { // A_ij = ||U_i - U_j||_2 
    determinePath(distanceMatrixName, -1, infoPath);
    assert(vec2);
  } else if (strcmp(type, "basisNormalizedCenterProducts") == 0) { // (rob_i^T * Center_p / || Center_p ||_2)
    determinePath(basisNormalizedCenterProductsName, -1, infoPath);
    assert(vec3);
  } else if (strcmp(type, "basisBasisProducts") == 0) {               
    determinePath(basisBasisProductsName, -1, infoPath);
    assert(vec4);
  } else if (strcmp(type, "basisUrefProducts") == 0) {               
    determinePath(basisUrefProductsName, -1, infoPath);
    assert(vec3);
  } else if (strcmp(type, "urefUrefProducts") == 0) {                
    determinePath(urefUrefProductsName, -1, infoPath);
    assert(vec2);
  } else if (strcmp(type, "urefUicProducts") == 0) {               
    determinePath(urefUicProductsName, -1, infoPath);
    assert(vec1);
  } else if (strcmp(type, "urefMultiUicProducts") == 0) {
    determinePath(urefMultiUicProductsName, -1, infoPath);
    assert(vec2);
  } else if (strcmp(type, "basisUicProducts") == 0) {               
    determinePath(basisUicProductsName, -1, infoPath);
    assert(vec2);
  } else if (strcmp(type, "basisMultiUicProducts") == 0) {
    determinePath(basisMultiUicProductsName, iCluster, infoPath);
    assert(vec3);
  } else if (strcmp(type, "urefComponentwiseSums") == 0) {               
    determinePath(urefComponentwiseSumsName, -1, infoPath);
    assert(vec2);
  } else if (strcmp(type, "basisComponentwiseSums") == 0) {               
    determinePath(basisComponentwiseSumsName, -1, infoPath);
    assert(vec3);
  } else if (strcmp(type, "correlationMatrix") == 0) {
    determinePath(correlationMatrixName, iCluster, infoPath);
    assert(vec2);
  } else if (strcmp(type, "multiUicMultiUicProducts") == 0) {                
    determinePath(multiUicMultiUicProductsName, -1, infoPath);
    assert(vec2);
  } else {
    fprintf(stderr, "*** Error: unexpected TYPE in readClusteredInfoAscii()\n");
    exit(-1);
  }

  readMultiVecASCII(infoPath, vec1, vec2, vec3, vec4);

  delete [] infoPath;
  infoPath = NULL;

} 

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::writeMultiVecASCII(char* path, std::vector<double>* vec1,
                                           std::vector<std::vector<double> >* vec2,
                                           std::vector<std::vector<std::vector<double> > >* vec3,
                                           std::vector<std::vector<std::vector<std::vector<double> > > >* vec4) {
                                                           
  // general IO function for multi-vector quantities with arbitrary dimensions
 
  // sanity check inputs
  if ( (vec1&&(vec2||vec3||vec4)) || (vec2&&(vec3||vec4)) || (vec3&&vec4)) {
    fprintf(stderr, "*** Error: only one multivec can be output at a time\n");
    exit(-1);
  } else if ((!vec1)&&(!vec2)&&(!vec3)&&(!vec4)) {
    fprintf(stderr, "*** Error: no multivec specified\n");
    exit(-1); 
  }

  if (com->cpuNum() == 0) {

    int dim1, dim2, dim3, dim4, multiVecType;  

    FILE *outputFile = fopen(path, "wt");

    if (vec1) {
      dim1 = vec1->size();
      multiVecType = 1;
    } else if (vec2) {
      dim1 = vec2->size();
      multiVecType = 2;
    } else if (vec3) {
      dim1 = vec3->size();
      multiVecType = 3;
    } else {
      dim1 = vec4->size();
      multiVecType = 4;
    }

    com->fprintf(outputFile, "MultiVecType: %d\n", multiVecType);
    com->fprintf(outputFile, "Dimension#1: %d\n", dim1);
    for (int i=0; i<dim1; ++i){
      if (vec1) {
        com->fprintf(outputFile,"%23.15e\n", (*vec1)[i]);
      } else {
        if (vec2) dim2 = (*vec2)[i].size();
        if (vec3) dim2 = (*vec3)[i].size();
        if (vec4) dim2 = (*vec4)[i].size();
        com->fprintf(outputFile, "Dimension#2: %d\n", dim2);
        for (int j=0; j<dim2; ++j) {
          if (vec2) {
            com->fprintf(outputFile,"%23.15e\n", (*vec2)[i][j]);
          } else {
            if (vec3) dim3 = (*vec3)[i][j].size();
            if (vec4) dim3 = (*vec4)[i][j].size();
            com->fprintf(outputFile, "Dimension#3: %d\n", dim3);
            for (int k=0; k<dim3; ++k) {
              if (vec3) {
                com->fprintf(outputFile,"%23.15e\n", (*vec3)[i][j][k]);
              } else {
                dim4 = (*vec4)[i][j][k].size();
                com->fprintf(outputFile, "Dimension#4: %d\n", dim4);
                for (int l=0; l<dim4; ++l) {
                  com->fprintf(outputFile,"%23.15e\n", (*vec4)[i][j][k][l]);
                }
              }
            }
          }
        }
      }
    }
  
    fclose (outputFile);

  }

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readMultiVecASCII(char* path, std::vector<double>* vec1,
                                          std::vector<std::vector<double> >* vec2,
                                          std::vector<std::vector<std::vector<double> > >* vec3,
                                          std::vector<std::vector<std::vector<std::vector<double> > > >* vec4) {
                                                           
  // general IO function for multi-vector quantities with arbitrary dimensions
                                    
  // sanity check inputs
  if ( (vec1&&(vec2||vec3||vec4)) || (vec2&&(vec3||vec4)) || (vec3&&vec4)) {
    fprintf(stderr, "*** Error: only one multivec can be read at a time\n");
    exit(-1);
  } else if ((!vec1)&&(!vec2)&&(!vec3)&&(!vec4)) {
    fprintf(stderr, "*** Error: no multivec specified\n");
    exit(-1); 
  }

  FILE *inputFile = fopen(path, "r");
  if (inputFile==NULL) {
    fprintf(stderr, "*** Error: could not open %s\n",path);
    exit(-1); 
  }

  int expectedMultiVecType;

  if (vec1) { 
    vec1->clear();
    expectedMultiVecType = 1;
  } else if (vec2) {
    vec2->clear();
    expectedMultiVecType = 2;
  } else if (vec3) {
    vec3->clear();
    expectedMultiVecType = 3;
  } else {
    vec4->clear();
    expectedMultiVecType = 4;
  }

  int _n, dim1, dim2, dim3, dim4, multiVecType;
  double tmpVal;

  fscanf(inputFile, "MultiVecType: %d\n", &multiVecType);
  assert(expectedMultiVecType == multiVecType);
  fscanf(inputFile, "Dimension#1: %d\n", &dim1);

  if (vec1) vec1->reserve(dim1);
  if (vec2) vec2->resize(dim1);
  if (vec3) vec3->resize(dim1);
  if (vec4) vec4->resize(dim1);

  for (int i=0; i<dim1; ++i){
    if (vec1) {
      _n = fscanf(inputFile,"%le\n", &tmpVal);
      vec1->push_back(tmpVal);
    } else {
      _n = fscanf(inputFile, "Dimension#2: %d\n", &dim2);
      if (vec2) (*vec2)[i].reserve(dim2);
      if (vec3) (*vec3)[i].resize(dim2);
      if (vec4) (*vec4)[i].resize(dim2);
      for (int j=0; j<dim2; ++j) {
        if (vec2) {
          _n = fscanf(inputFile,"%le\n", &tmpVal);
          (*vec2)[i].push_back(tmpVal);
        } else {
          _n = fscanf(inputFile, "Dimension#3: %d\n", &dim3);
          if (vec3) (*vec3)[i][j].reserve(dim3);
          if (vec4) (*vec4)[i][j].resize(dim3);
          for (int k=0; k<dim3; ++k) {
            if (vec3) {
              _n = fscanf(inputFile,"%le\n", &tmpVal);
              (*vec3)[i][j].push_back(tmpVal);
            } else {
              _n = fscanf(inputFile, "Dimension#4: %d\n", &dim4);
              (*vec4)[i][j][k].reserve(dim4);
              for (int l=0; l<dim4; ++l) {
                _n = fscanf(inputFile,"%le\n", &tmpVal);
                (*vec4)[i][j][k].push_back(tmpVal);
              }
            }
          }
        }
      }
    }
  }

  fclose (inputFile);

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::readDistanceComparisonInfo(const char* updateType) {

  if (nClusters == 1) return;

  checkInitialConditionScenario();
 
  if (incrementalStateSnaps) {

    readClusteredInfoASCII(-1, "basisNormalizedCenterProducts", NULL, NULL, &basisNormalizedCenterProducts);
 
  } else if ((strcmp(updateType, "noUpdates") == 0) ||
             (strcmp(updateType, "exactUpdates") == 0) ||
             (strcmp(updateType, "project") == 0) ||
             (strcmp(updateType, "approxUpdates") == 0)){

    readCenterNorms();

    if ((strcmp(updateType, "noUpdates") == 0) ||
             (strcmp(updateType, "exactUpdates") == 0) ||
             (strcmp(updateType, "project") == 0)) {

      stateBasisCentersDifProduct.resize(nClusters);

      if (ioData->romOnline.krylov.include==NonlinearRomOnlineNonStateData::INCLUDE_ON)
        krylovBasisCentersDifProduct.resize(nClusters); 

      for (int iCluster=0; iCluster<nClusters; ++iCluster) {
        readClusteredInfoASCII(iCluster, "state", NULL, NULL, &stateBasisCentersDifProduct[iCluster]);
    
        if (ioData->romOnline.krylov.include==NonlinearRomOnlineNonStateData::INCLUDE_ON)
          readClusteredInfoASCII(iCluster, "krylov", NULL, NULL, &krylovBasisCentersDifProduct[iCluster]); 
      }

      if (ioData->romOnline.sensitivity.include==NonlinearRomOnlineNonStateData::INCLUDE_ON)
        readClusteredInfoASCII(-2, "sensitivity", NULL, NULL, &sensitivityBasisCentersDifProduct);      
    }
   
    if (specifiedIC) {
      readClusteredInfoASCII(-1, "initialCondition", NULL, &initialConditionCentersDifProduct);
    } else if (interpolatedMultiIC) {

      readClusteredInfoASCII(-1, "multiInitialCondition", NULL, NULL, &multiUicCentersDifProduct);
      initialConditionCentersDifProduct.resize(nClusters);
      for (int iCluster=0; iCluster<nClusters; ++iCluster) {
        int rowLength = multiUicCentersDifProduct[iCluster].size();
        initialConditionCentersDifProduct[iCluster].clear();
        initialConditionCentersDifProduct[iCluster].resize(rowLength,0.0);
        for (int jCluster=0; jCluster<rowLength; ++jCluster) {
          for (int iSol=0; iSol<multiUicCentersDifProduct[iCluster][jCluster].size(); ++iSol) {
            initialConditionCentersDifProduct[iCluster][jCluster] +=
              multiUicCentersDifProduct[iCluster][jCluster][iSol]*interpWeightsForMultiIC[iSol];
          }
        }
      } 
      multiUicCentersDifProduct.clear();
   
    } else {
      // this will be constructed during initializeDistanceComparisons 
    }

    if ((strcmp(updateType, "exactUpdates") == 0) || (strcmp(updateType, "project") == 0)) {
      refStateCentersDifProduct.resize(nClusters);
      for (int iCluster=0; iCluster<nClusters; ++iCluster) {
        readClusteredInfoASCII(iCluster, "referenceState", NULL, &refStateCentersDifProduct[iCluster]);
      }
    }

    if (strcmp(updateType, "approxUpdates") == 0) {   
      if (ioData->romOnline.systemApproximation == NonlinearRomOnlineData::SYSTEM_APPROXIMATION_NONE) {
        this->readClusterCenters("centers");
      } else {
        this->readClusterCenters("sampledCenters");
      }
    }  

  } else {
    fprintf(stderr, "*** Error: unexpected update method\n");
    exit(-1);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::qr(VecSet< DistSVec<double, dim> >* Q, std::vector<std::vector<double> >* RT, bool testQR) {

  int nVec = Q->numVectors();

  std::vector<std::vector<double> >* RT_testing = NULL;
  VecSet< DistSVec<double, dim> >* testVecSet = NULL;
  if (testQR) {
    testVecSet = new VecSet< DistSVec<double, dim> >(nVec, domain.getNodeDistInfo());
    for (int iVec = 0; iVec<nVec; ++iVec) {
      (*testVecSet)[iVec] = (*Q)[iVec];
    }  
    if (!RT) {
        RT_testing = new std::vector<std::vector<double> >;
        RT = RT_testing;  // need R to test the QR...
    }
  }

  if (RT) {
    RT->clear();
    RT->resize(nVec);
    for (int iVec = 0; iVec<nVec; ++iVec)
      (*RT)[iVec].resize(iVec+1, 0.0);
  }

  for (int iVec = 0; iVec<nVec; ++iVec) {

    double norm = (*Q)[iVec].norm();
    if (RT) (*RT)[iVec][iVec] = norm;
    if (norm>=1e-13) {
      (*Q)[iVec] *= 1/norm;
    } else {
      com->fprintf(stderr, "*** Warning: QR encountered a rank defficient matrix in GNAT preprocessing (norm %e)\n",norm);
      // resize everything
      int newSize = nVec-1;
      for (int jVec = iVec+1; jVec<nVec; ++jVec)
        (*Q)[jVec-1] = (*Q)[jVec];
      Q->resize(newSize);

      if (testVecSet) {
        for (int jVec = iVec+1; jVec<nVec; ++jVec)
          (*testVecSet)[jVec-1] = (*testVecSet)[jVec];
        testVecSet->resize(newSize);
      }

      if (RT)
        RT->resize(newSize);

      nVec = newSize;
    }

    for (int jVec = iVec+1; jVec<nVec; ++jVec) {
      double r = (*Q)[iVec] * (*Q)[jVec];
      (*Q)[jVec] -= (*Q)[iVec] * r;
      if (RT) (*RT)[jVec][iVec] = r;
    }

  }

  if (testQR) {
    // Q orthogonal?
    com->fprintf(stdout, "\n ... testing whether Q is orthogonal\n");
    std::vector<int> errorLog;
    int nErrorBins = 11;
    errorLog.resize(nErrorBins,0);
    for (int iVec = 0; iVec<nVec; ++iVec) {
      for (int jVec = 0; jVec<=iVec; ++jVec) {
        double product = (*Q)[iVec] * (*Q)[jVec];
        double error = (iVec==jVec) ? abs(product - 1.0) : abs(product);
        if (error<=pow(10.0,-1.0*(double)(nErrorBins-1))) {
          ++errorLog[nErrorBins-1];
        } else {
          for (int iErr = 0; iErr<nErrorBins-1; ++iErr) {
            if ((log10(error)<=-1.0*(double)iErr) && (log10(error)>-1.0*(double)(iErr+1))) 
              ++errorLog[iErr];
          }
        }
      }
    }
    for (int iErr = 0; iErr<nErrorBins-1; ++iErr)
      com->fprintf(stdout, " ... ... %d errors between 1e-%d and 1e-%d\n", errorLog[iErr], iErr, iErr+1);
    com->fprintf(stdout, " ... ... %d errors less than 1e-%d\n", errorLog[nErrorBins-1],nErrorBins-1);

    // QR == original matrix?
    double maxError = 0.0;
    double avgError = 0.0;
    DistSVec<double, dim> testVec(domain.getNodeDistInfo());
    com->fprintf(stdout, " ... testing whether QR recovers the original matrix\n");
    for (int iVec = 0; iVec<nVec; ++iVec) {
      testVec = (*testVecSet)[iVec];
      for (int jVec = 0; jVec<=iVec; ++jVec) {
        testVec -= (*Q)[jVec]*(*RT)[iVec][jVec];
      }
      double error = testVec.norm();
      if (error>maxError) maxError=error;
      avgError += error/nVec;
    }

    com->fprintf(stdout, " ... ... QR accuracy test: maxError = %e\n", maxError);
    com->fprintf(stdout, " ... ... QR accuracy test: avgError = %e\n", avgError);
    delete testVecSet;
    if (RT_testing) delete RT_testing;
  }


}

//---------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::rSVD(VecSet< DistSVec<double, dim> >& Utrue, std::vector<double>& singularValues, FullM& Vtrue, bool testSVD) {
  // Standard R-SVD algorithm (using gram schmidt for QR):
  //   QR = origMatrix
  //   [Utmp S V] = svd(R) (in serial using lapack)
  //   Communicate svd result to all mpi ranks
  //   U = Q*Utmp
  //
  // Note: the matrix is passed in via Utrue, which is overwritten

  int nVecs = Utrue.numVectors();

  VecSet< DistSVec<double, dim> > testVecSet(0, domain.getNodeDistInfo());
  if (testSVD) {// store the original matrix for comparison later
    testVecSet.resize(nVecs);
    for (int iVec=0; iVec<nVecs; ++iVec) 
      testVecSet[iVec]=Utrue[iVec];
  }

  VecSet< DistSVec<double, dim> > Q(nVecs, domain.getNodeDistInfo());
  Q.resize(nVecs);
  for (int iVec=0; iVec<nVecs; ++iVec)
    Q[iVec]=Utrue[iVec];

  // orthogonalization
  com->fprintf(stdout, " ... orthogonalizing");
  double orthTime = timer->getTime();
  std::vector<std::vector<double> > RTranspose;
  qr(&Q, &RTranspose, false); 
  orthTime = timer->getTime() - orthTime;
  com->fprintf(stdout," (%e sec)\n", orthTime);
  double* tmpMat = NULL;
  if (com->cpuNum()==0) {
    tmpMat = new double[nVecs*nVecs];
    for (int iVec=0; iVec<nVecs; ++iVec) {
      for (int jVec=0; jVec<nVecs; ++jVec) {
        tmpMat[iVec*nVecs + jVec] = (iVec >= jVec) ? RTranspose[iVec][jVec] : 0.0; 
      // tmpMat[iVec*nVecs + jVec] =  tmpMat[jVec][iVec]
      }
    }
  }

  // SVD
  double lapackTime = timer->getTime();
  com->fprintf(stdout, " ... computing SVD using linpack");

  double *yVec;
  double *zVec;
  double *singVals;

  if (com->cpuNum()==0) {
    linpackSVD(tmpMat, nVecs, nVecs, yVec, singVals, zVec);
    delete [] tmpMat;
  }

  com->barrier();
  lapackTime = timer->getTime() - lapackTime;
  com->fprintf(stdout," (%e sec)\n", lapackTime);

  com->fprintf(stdout, " ... broadcasting result (and forming U)");
  double broadcastTime = timer->getTime();
  // broadcast singular values
  if (com->cpuNum()!=0) {
    singVals = new double[nVecs];
  }
  com->broadcast(nVecs, singVals, 0);
  com->barrier();

  singularValues.resize(nVecs);
  for (int iVec=0; iVec<nVecs; ++iVec) {
    singularValues[iVec]=singVals[iVec];
  }
  delete [] singVals;

  // broadcast left singular vectors
  if (com->cpuNum()!=0) {
    yVec = new double[nVecs*nVecs];
  }
  com->broadcast(nVecs*nVecs, yVec, 0);
  com->barrier();

  for (int iVec=0; iVec<nVecs; ++iVec) {
    Utrue[iVec] = 0.0;
    for (int jVec=0; jVec<nVecs; ++jVec) {
      Utrue[iVec] += Q[jVec]*yVec[iVec*nVecs + jVec];
    }
  }
  Q.resize(0);
  delete [] yVec;

  // broadcast right singular vectors
  if (com->cpuNum()!=0) {
    zVec = new double[nVecs*nVecs]; 
  }
  com->broadcast(nVecs*nVecs, zVec, 0);
  com->barrier();

  Vtrue.setNewSize(nVecs);
  for (int iVec=0; iVec<nVecs; ++iVec) {
    for (int jVec=0; jVec<nVecs; ++jVec) {
      Vtrue[jVec][iVec] = zVec[iVec*nVecs + jVec];
    }
  }
  delete [] zVec;

  broadcastTime = timer->getTime() - broadcastTime;
  com->fprintf(stdout," (%e sec)\n", broadcastTime);

  if (testSVD) {
    //check svd
    com->fprintf(stdout," ... checking the SVD (debugging)");
    double checkTime = timer->getTime();
    double errorNorm,maxErr,avgErr;
    DistSVec<double,dim> errorVec( domain.getNodeDistInfo() );
    maxErr = 0.0;
    avgErr = 0.0;
    for (int iVec = 0; iVec < nVecs; ++iVec) {
      errorVec = testVecSet[iVec];
      for (int jVec = 0; jVec < nVecs; ++jVec)
        errorVec = errorVec - ((singularValues[jVec] * Vtrue[iVec][jVec]) * Utrue[jVec]);
      errorNorm = (((testVecSet[iVec]).norm()) > 1e-15) ? errorVec.norm()/((testVecSet[iVec]).norm()) : 0.0;
      avgErr += errorNorm;
      if (errorNorm > maxErr)
        maxErr = errorNorm;
    }
    avgErr /= nVecs;
  
    checkTime = timer->getTime() - checkTime;
    com->fprintf(stdout," (%e sec)\n", checkTime);
  
    com->fprintf(stdout, " ... Average error on Snapshots after SVD = %e\n", avgErr);
    com->fprintf(stdout, " ... Maximum error on Snapshots after SVD = %e\n", maxErr);
    // end check svd  
    testVecSet.resize(0);
  }


}

//-------------------------------------------------------------------------------------------


template<int dim>
void NonlinearRom<dim>::probabilisticSVD(VecSet< DistSVec<double, dim> >& Utrue, std::vector<double>& singularValues, FullM& Vtrue,
                                           int k, int nPowerIts, bool testSVD) {
  // Strategy:
  //   CPU0 generates random martrix, randMat, which is nVecs x k with k <= nVecs
  //    (if k==nVecs, then just perform a standard R-SVD)
  //   Communicate randMat to all mpi ranks
  //   tmpMat = origMatrix * randMat
  //   for i = 1:nPowerIts
  //      tmpMat = podHat * (podHat' * tmpMat)
  //   end
  //   [Q,~] = qr(tmpMat)
  //   tmpMat2 = Q' * origMatrix
  //   [Utmp S V] = svd(tmpMat2) (in serial using lapack)
  //   Communicate svd result to all mpi ranks
  //   U = Q*Utmp
  //
  // Note: the matrix is passed in via Utrue, which is overwritten


  int nVecs = Utrue.numVectors();
  k = min(nVecs,k);

  VecSet< DistSVec<double, dim> > testVecSet(0, domain.getNodeDistInfo());
  if (testSVD) {// store the original matrix for comparison later
    testVecSet.resize(nVecs);
    for (int iVec=0; iVec<nVecs; ++iVec) 
      testVecSet[iVec]=Utrue[iVec];
  }

  VecSet< DistSVec<double, dim> > tmpVecSet(k, domain.getNodeDistInfo());
  bool solveLinearSystem = (nPowerIts<=0 && k==nVecs); // for Q^T * Snapshots in the projection step

  double* randMat = NULL;
  if ((com->cpuNum()==0) && solveLinearSystem) randMat = new double[nVecs*k]; // only store for cpu0

  int seed = 1;
  boost::mt19937 randSeed(seed);
  boost::normal_distribution<> normalDistribution(0.0, 1.0);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > normalGenerator(randSeed, normalDistribution);

  // tmp = snaps * randMat
  com->fprintf(stdout, " ... right multiplying by %dx%d random matrix", nVecs, k);
  double multTime = timer->getTime();
  double randVal;
  for (int iVec=0; iVec<k; ++iVec) {
    tmpVecSet[iVec] = 0.0;
    for (int jVec=0; jVec<nVecs; ++jVec) {
      randVal = normalGenerator();
      tmpVecSet[iVec] += Utrue[jVec]*randVal;
      if (randMat) randMat[iVec*nVecs + jVec] = randVal;
    }
  }
  multTime = timer->getTime() - multTime;
  com->fprintf(stdout," (%e sec)\n", multTime);

  // power iteration
  for (int iPowerIt=0; iPowerIt<nPowerIts; ++iPowerIt) {
    double powerItTime = timer->getTime();
    com->fprintf(stdout, " ... power iteration #%d of %d",iPowerIt,nPowerIts);
    double* tmpMat = new double[nVecs*k]; // store as vector for simplicity  
    for (int iVec=0; iVec<k; ++iVec) {
      for (int jVec=0; jVec<nVecs; ++jVec) {
        tmpMat[iVec*nVecs + jVec] = Utrue[jVec]*tmpVecSet[iVec]; // tmpMat[jSnap][iSnap]
      }
    }
    for (int iVec=0; iVec<k; ++iVec) {
      tmpVecSet[iVec] = 0.0;
      for (int jVec=0; jVec<nVecs; ++jVec) {
        tmpVecSet[iVec] += Utrue[jVec]*tmpMat[iVec*nVecs + jVec];
      }
    }
    delete [] tmpMat;
    powerItTime = timer->getTime() - powerItTime;
    com->fprintf(stdout," (%e sec)\n", powerItTime);
  }

  // orthogonalization
  com->fprintf(stdout, " ... orthogonalizing");
  double orthTime = timer->getTime();
  std::vector<std::vector<double> > RTranspose;
  if (solveLinearSystem) {
    qr(&tmpVecSet, &RTranspose, false); 
  } else {
    qr(&tmpVecSet, NULL, false);
  }
  if (tmpVecSet.numVectors() < k) { //QR might return fewer vectors if tmpVecSet wasn't full rank.
    k=tmpVecSet.numVectors();
    solveLinearSystem = false;
    if (randMat) {
      delete [] randMat;
      randMat = NULL;
    }
  }
  orthTime = timer->getTime() - orthTime;
  com->fprintf(stdout," (%e sec)\n", orthTime);

  // projection
  com->fprintf(stdout, " ... projecting");
  double projTime = timer->getTime();
  double* tmpMat = NULL;
  if (com->cpuNum()==0) {
    tmpMat = new double[k*nVecs];
  } 
  if (solveLinearSystem) { // solve a small system
    if (com->cpuNum()==0) {
      // set up the left hand side
      FullM randMatTranspose(k);
      for (int iEntry=0; iEntry<k; ++iEntry) {
        for (int jEntry=0; jEntry<k; ++jEntry) {
          randMatTranspose[iEntry][jEntry] = randMat[k*iEntry+jEntry]; 
        }
      }
      delete [] randMat;
      randMat = NULL;
      randMatTranspose.Factor();

      double *RHS = new double[k];
      for (int iRHS=0; iRHS<k; ++iRHS) {
        // Fill RHS from RTranspose
        for (int iEntry=0; iEntry<k; ++iEntry) {
          RHS[iEntry] = (iEntry>=iRHS) ? RTranspose[iEntry][iRHS] : 0.0; 
        }
        // Solve for i-th row of Q^T * Snapshots
        randMatTranspose.ReSolve(RHS);
        // Store result in tmpMat
        for (int iEntry=0; iEntry<k; ++iEntry) {
          tmpMat[iEntry*k + iRHS] = RHS[iEntry];
        }
      }
      delete [] RHS;
    }
    // testing
    /*projTime = timer->getTime() - projTime;
    com->fprintf(stdout," (%e sec)\n", projTime);
    projTime = timer->getTime();
    double tmpProduct;
    for (int iVec=0; iVec<k; ++iVec) {
      for (int jVec=0; jVec<k; ++jVec) {
        tmpProduct = tmpVecSet[jVec]*Utrue[iVec]; // tmpMat[jVec][iVec]
        if (com->cpuNum()==0) {
          if (abs(tmpMat[iVec*k + jVec] - tmpProduct)>1e-6)  
            com->fprintf(stderr,"... iVec=%d, jVec=%d, solve=%e, direct=%e)\n", iVec, jVec, tmpMat[iVec*k + jVec], tmpProduct);
        }
      }
     } */
  } else { // just compute the matrix matrix product directly.
    double tmpProduct = 0.0;
    for (int iVec=0; iVec<nVecs; ++iVec) {
      for (int jVec=0; jVec<k; ++jVec) {
        tmpProduct = tmpVecSet[jVec]*Utrue[iVec]; // tmpMat[jVec][iVec]
        if (com->cpuNum()==0) tmpMat[iVec*k + jVec] = tmpProduct; 
      }
    }
  }
  projTime = timer->getTime() - projTime;
  com->fprintf(stdout," (%e sec)\n", projTime);

  // SVD
  double lapackTime = timer->getTime();
  com->fprintf(stdout, " ... computing SVD using linpack");

  double *yVec;
  double *zVec;
  double *singVals;

  if (com->cpuNum()==0) {
    linpackSVD(tmpMat, k, nVecs, yVec, singVals, zVec);
    delete [] tmpMat;
  }

  com->barrier();
  lapackTime = timer->getTime() - lapackTime;
  com->fprintf(stdout," (%e sec)\n", lapackTime);

  com->fprintf(stdout, " ... broadcasting result (and forming U)");
  double broadcastTime = timer->getTime();
  // broadcast singular values
  // only the first k should be nonzero
  if (com->cpuNum()!=0) {
    singVals = new double[k];
  }
  com->broadcast(k, singVals, 0);
  com->barrier();

  singularValues.resize(k);
  for (int iVec=0; iVec<k; ++iVec) {
    singularValues[iVec]=singVals[iVec];
  }
  delete [] singVals;

  // broadcast left singular vectors
  if (com->cpuNum()!=0) {
    yVec = new double[k*k];
  }
  com->broadcast(k*k, yVec, 0);
  com->barrier();

  Utrue.resize(k);
  for (int iVec=0; iVec<k; ++iVec) {
    Utrue[iVec] = 0.0;
    for (int jVec=0; jVec<k; ++jVec) {
      Utrue[iVec] += tmpVecSet[jVec]*yVec[iVec*k + jVec];
    }
  }
  tmpVecSet.resize(0);
  delete [] yVec;

  // broadcast right singular vectors
  if (com->cpuNum()!=0) {
    zVec = new double[nVecs*nVecs]; 
  }
  com->broadcast(nVecs*nVecs, zVec, 0);
  com->barrier();

  Vtrue.setNewSize(nVecs);
  for (int iVec=0; iVec<nVecs; ++iVec) {
    for (int jVec=0; jVec<nVecs; ++jVec) {
      Vtrue[jVec][iVec] = zVec[iVec*nVecs + jVec];
    }
  }
  delete [] zVec;

  broadcastTime = timer->getTime() - broadcastTime;
  com->fprintf(stdout," (%e sec)\n", broadcastTime);

  if (testSVD) {
    //check svd
    com->fprintf(stdout," ... checking the SVD (debugging)");
    double checkTime = timer->getTime();
    double errorNorm,maxErr,avgErr;
    DistSVec<double,dim> errorVec( domain.getNodeDistInfo() );
    maxErr = 0.0;
    avgErr = 0.0;
    for (int iVec = 0; iVec < nVecs; ++iVec) {
      errorVec = testVecSet[iVec];
      for (int jVec = 0; jVec < k; ++jVec)
        errorVec = errorVec - ((singularValues[jVec] * Vtrue[iVec][jVec]) * Utrue[jVec]);
      errorNorm = (((testVecSet[iVec]).norm()) > 1e-15) ? errorVec.norm()/((testVecSet[iVec]).norm()) : 0.0;
      avgErr += errorNorm;
      if (errorNorm > maxErr)
        maxErr = errorNorm;
    }
    avgErr /= nVecs;
  
    checkTime = timer->getTime() - checkTime;
    com->fprintf(stdout," (%e sec)\n", checkTime);
  
    com->fprintf(stdout, " ... Average error on Snapshots after SVD = %e\n", avgErr);
    com->fprintf(stdout, " ... Maximum error on Snapshots after SVD = %e\n", maxErr);
    // end check svd  
    testVecSet.resize(0);
  }


}

//-------------------------------------------------------------------------------------------

template<int dim>
void NonlinearRom<dim>::linpackSVD(double* tmpMat, int rows, int columns, double*& yVec, double*& singVals, double*& zVec) {
/*subroutine dsvdc(x,ldx,n,p,s,e,u,ldu,v,ldv,work,job,info)
ON ENTRY

x    double precision(ldx,p), where ldx.ge.n.  x contains the matrix whose singular value decomposition is to be computed. x is destroyed by dsvdc.
ldx  integer. ldx is the leading dimension of the array x.
n    integer. n is the number of rows of the matrix x.
p    integer. p is the number of columns of the matrix x.
ldu  integer. ldu is the leading dimension of the array u. (see below).
ldv  integer. ldv is the leading dimension of the array v. (see below).
work double precision(n). work is a scratch array.
job  integer. job controls the computation of the singular vectors.  it has the decimal expansion ab with the following meaning
         a.eq.0    do not compute the left singular vectors.
         a.eq.1    return the n left singular vector in u.
         a.ge.2    return the first min(n,p) singular vectors in u.
         b.eq.0    do not compute the right singular vectors.
         b.eq.1    return the right singular vectors in v.

ON RETURN

s   double precision(mm), where mm=min(n+1,p). the first min(n,p) entries of s contain the singular values of x arranged in descending order of magnitude.
e   double precision(p), e ordinarily contains zeros.  however see the discussion of info for exceptions.
u   double precision(ldu,k), where ldu.ge.n.  if joba.eq.1 then k.eq.n, if joba.ge.2 then k.eq.min(n,p). u contains the matrix of left singular vectors. u is not referenced if joba.eq.0.  if n.le.p or if joba.eq.2, then u may be identified with x in the subroutine call.
v   double precision(ldv,p), where ldv.ge.p. v contains the matrix of right singular vectors. v is not referenced if job.eq.0.  if p.le.n, then v may be identified with x in the subroutine call.*/
    double *err = new double[columns];
    double *work = new double[rows];
    int info;
    int job;
    if (rows<=columns) { // fat or square
      job = 11;
      yVec = new double[rows*rows]; // left singular vectors
    } else { // thin
      job = 21;
      yVec = new double[rows*columns]; // left singular vectors
    }
    zVec = new double[columns*columns]; // right singular vectors
    singVals = new double[min(columns,rows+1)];

    F77NAME(dsvdc)(tmpMat, rows, rows, columns, singVals, err, yVec, rows, zVec, columns, work, job, info);

    if (info>0) {
      fprintf(stderr,"\n*** Error (mesg 1 of 2): linpack failed -- unable to find the first %d singular values\n", info);
      exit(-1);
    }

    delete [] err;
    delete [] work;


}

//------------------------------------------------------------------------------
template<int dim>
void NonlinearRom<dim>::probabilisticLSMultiRHS(VecSet< DistSVec<double, dim> >& LHS, VecSet< DistSVec<double, dim> >& RHS,
                                                  std::vector<std::vector<double> >& lsCoeffVec, int k, int nPowerIts, bool testSVD) {

  // Multi-RHS Least Squares via probabilistic SVD
  // LHS is destroyed by probabilisticSVD, which returns Utrue (size k) in its place.  LHS is resized to zero before this function returns.

  int nRhs = RHS.numVectors();
  int nLhs = LHS.numVectors();

  // first call probabilistic SVD
  std::vector<double> singularValues;
  FullM Vtrue;
  probabilisticSVD(LHS, singularValues, Vtrue, k, nPowerIts, testSVD);

  // now compute LS coefficients
  //   lsCoeff = V * inv(S) * U' * RHS
  com->fprintf(stdout," ... computing least squares coefficients using SVD");

  double lscoeffTime = timer->getTime();
  lsCoeffVec.clear();
  lsCoeffVec.resize(nRhs);
  for (int iRhs=0; iRhs<nRhs; ++iRhs) {
    lsCoeffVec[iRhs].resize(nLhs,0.0);
  }

  std::vector< double > result;
  result.resize(k, 0.0);
  for (int iRhs=0; iRhs<nRhs; ++iRhs) {
    for (int iVec=0; iVec<LHS.numVectors(); ++iVec) {
      result[iVec] = LHS[iVec] * RHS[iRhs];  // note: LHS was destroyed by probabilisticSVD and is now Utrue (size k)
      result[iVec] = (singularValues[iVec]>1e-15) ? result[iVec] / singularValues[iVec] : 0;
    }

    for (int iLhs=0; iLhs<nLhs; ++iLhs) {
      lsCoeffVec[iRhs][iLhs] = 0.0;
      for (int iVec=0; iVec<LHS.numVectors(); ++iVec) {
        lsCoeffVec[iRhs][iLhs] += Vtrue[iLhs][iVec]*result[iVec];
      }
    }
  }
  result.clear();
  LHS.resize(0);

  lscoeffTime = timer->getTime() - lscoeffTime;
  com->fprintf(stdout," (%e sec)\n", lscoeffTime);

}

//------------------------------------------------------------------------------

/*
  if (staterom) {
    if (it0 != 0)
      fpStateRom = backupAsciiFile(staterom);
    if (it0 == 0 || fpStateRom == 0) {
      fpStateRom = fopen(staterom, "w");
      if (!fpStateRom) {
  fprintf(stderr, "*** Error: could not open \'%s\'\n", staterom);
  exit(1);
      }
      fprintf(fpStateRom, "# TimeIteration ElapsedTime StatePodCoords \n");
    }
    fflush(fpStateRom);
  }


  if (iod.output.rom.staterom[0] != 0) {
    staterom = new char[sprom + strlen(iod.output.rom.staterom)];
    sprintf(staterom, "%s%s", iod.output.rom.prefix, iod.output.rom.staterom);
  }
  else
    staterom = 0;

  if (fpStateRom) fclose(fpStateRom);
*/

//------------------------------------------------------------------------------
template<int dim>
void NonlinearRom<dim>::truncateBufferedBasis() {

  int nPod = basis->numVectors();
  int nPodNew = nPod-nBuffer;

  if (nBuffer>0) {

    com->fprintf(stdout, " ... truncating buffered basis from %d vectors to %d vectors\n", nPod, nPodNew);

    VecSet< DistSVec<double, dim> >* basisNew =  new VecSet< DistSVec<double, dim> >(nPodNew, domain.getNodeDistInfo());

    for (int iVec=0; iVec<nPodNew; ++iVec) {
      (*basisNew)[iVec]=(*basis)[iVec];
    }

    delete basis;
    basis = basisNew;
    basisNew = NULL;

    sVals->resize(nPodNew);
  }

  nBuffer = 0;

}

//------------------------------------------------------------------------------
template<int dim>
void NonlinearRom<dim>::partitionAndSowerForGappy(bool surfaceMeshConstruction) {

  if (strcmp(ioData->input.metis,"")==0 || strcmp(ioData->input.sower,"")==0) {
    com->fprintf(stdout, " ... consider specifying the METIS and SOWER executables to save yourself some work\n");
    return; 
  }

  if (com->cpuNum() == 0) {
    if (surfaceMeshConstruction) {

      char *topFilePath = NULL;
      determinePath(surfaceMeshName, -1, topFilePath);

      // call metis
      FILE *shell;
      std::string metisCommandString(ioData->input.metis);
      metisCommandString += " ";
      metisCommandString += topFilePath;
      metisCommandString += " ";
      metisCommandString += boost::lexical_cast<std::string>(ioData->input.nParts);
      const char *metisCommandChar = metisCommandString.c_str();

      com->fprintf(stdout, "\n%s\n", metisCommandChar);
      if (!(shell = popen(metisCommandChar, "r"))) {
        fprintf(stderr, " *** Error: attempt to use external METIS executable (%s) failed!\n", ioData->input.metis);
        exit(-1);
      } else {
        com->fprintf(stdout, "\n ... Calling external METIS executable (%s) ...\n", ioData->input.metis);
      }

      char buff[512];
      while (fgets(buff, sizeof(buff), shell)!=NULL){
        com->fprintf(stdout, "%s", buff);
      }
      pclose(shell);

      std::string decompositionPathString(topFilePath);
      decompositionPathString += ".dec.";
      decompositionPathString += boost::lexical_cast<std::string>(ioData->input.nParts);

      // initial call to "sower -fluid"
      std::string surfacePrefixPathString(databasePrefix);
      surfacePrefixPathString += databaseName;
      surfacePrefixPathString += "/";
      surfacePrefixPathString += romFiles->surfacePrefix;
      
      std::vector<int> cpuMaps;
      int nCores = 1;
      while (nCores < ioData->input.nParts) {
        cpuMaps.push_back(nCores);
        nCores *= 2;
      }
      cpuMaps.push_back(ioData->input.nParts);

      std::string sowerCommandString = ioData->input.sower;
      sowerCommandString += " -fluid -mesh ";
      sowerCommandString += topFilePath;
      sowerCommandString += " -dec ";
      sowerCommandString += decompositionPathString;
      std::vector<int>::iterator it;
      for (it = cpuMaps.begin(); it != cpuMaps.end(); ++it) { 
        sowerCommandString += " -cpu ";
        sowerCommandString += boost::lexical_cast<std::string>(*it);
      }
      sowerCommandString += " -output ";
      sowerCommandString += surfacePrefixPathString;
      sowerCommandString += " -cluster ";
      sowerCommandString += boost::lexical_cast<std::string>(ioData->input.nParts);
      const char *sowerCommandChar = sowerCommandString.c_str();

      com->fprintf(stdout, "\n%s\n", sowerCommandChar);
      if (!(shell = popen(sowerCommandChar, "r"))) {
        fprintf(stderr, " *** Error: attempt to use external SOWER executable (%s) failed!\n", ioData->input.sower);
        exit(-1);
      } else {
        com->fprintf(stdout, "\n ... Calling external SOWER executable (%s) ...\n", ioData->input.sower);
      }

      while (fgets(buff, sizeof(buff), shell)!=NULL){
        com->fprintf(stdout, "%s", buff);
      }
      pclose(shell);

      std::string meshPathString(surfacePrefixPathString);
      meshPathString += ".msh";

      std::string connectivityPathString(surfacePrefixPathString);
      connectivityPathString += ".con";
    
      delete [] topFilePath;

      // now for all of the "sower -fluid -split" calls
    
      char *surfaceCentersPath = NULL;
      determinePath(surfaceCentersName, -1, surfaceCentersPath);
      callSowerSplit(meshPathString, connectivityPathString, surfaceCentersPath);
      delete [] surfaceCentersPath;
    
      char *surfaceSolutionPath = NULL;
      determinePath(surfaceSolutionName, -1, surfaceSolutionPath);
      callSowerSplit(meshPathString, connectivityPathString, surfaceSolutionPath);
      delete [] surfaceSolutionPath;

      char *surfaceMatchStatePath = NULL;
      determinePath(surfaceMatchStateName, -1, surfaceMatchStatePath);
      callSowerSplit(meshPathString, connectivityPathString, surfaceMatchStatePath);
      delete [] surfaceMatchStatePath;

      char *surfaceMultiSolutionsPath = NULL;
      determinePath(surfaceMultiSolutionsName, -1, surfaceMultiSolutionsPath);
      callSowerSplit(meshPathString, connectivityPathString, surfaceMultiSolutionsPath);
      delete [] surfaceMultiSolutionsPath;   

      char *surfaceWallDistPath = NULL;
      determinePath(surfaceWallDistName, -1, surfaceWallDistPath);
      callSowerSplit(meshPathString, connectivityPathString, surfaceWallDistPath);
      delete [] surfaceWallDistPath;
    
      char *surfaceDisplacementPath = NULL;
      determinePath(surfaceDisplacementName, -1, surfaceDisplacementPath);
      callSowerSplit(meshPathString, connectivityPathString, surfaceDisplacementPath);
      delete [] surfaceDisplacementPath;

      char *surfaceShapeDerivativePath = NULL;
      determinePath(surfaceShapeDerivativeName, -1, surfaceShapeDerivativePath);
      callSowerSplit(meshPathString, connectivityPathString, surfaceShapeDerivativePath);
      delete [] surfaceShapeDerivativePath;

      char *approxMetricStateLowRankSurfaceCoordsPath = NULL;
      determinePath(approxMetricStateLowRankSurfaceCoordsName, -1, approxMetricStateLowRankSurfaceCoordsPath);
      callSowerSplit(meshPathString, connectivityPathString, approxMetricStateLowRankSurfaceCoordsPath);
      delete [] approxMetricStateLowRankSurfaceCoordsPath;   
 
      for (int iCluster=0; iCluster<nClusters; ++iCluster) {
    
        char *surfaceStateBasisPath = NULL;
        determinePath(surfaceStateBasisName, iCluster, surfaceStateBasisPath);
        callSowerSplit(meshPathString, connectivityPathString, surfaceStateBasisPath);
        delete [] surfaceStateBasisPath;
    
        char *surfaceRefStatePath = NULL;
        determinePath(surfaceRefStateName, iCluster, surfaceRefStatePath);
        callSowerSplit(meshPathString, connectivityPathString, surfaceRefStatePath);
        delete [] surfaceRefStatePath;
      }

    } else {
   
      char *topFilePath = NULL;
      determinePath(sampledMeshName, -1, topFilePath);

      // call metis
      FILE *shell;
      std::string metisCommandString(ioData->input.metis);
      metisCommandString += " ";
      metisCommandString += topFilePath;
      metisCommandString += " ";
      metisCommandString += boost::lexical_cast<std::string>(ioData->input.nParts);
      const char *metisCommandChar = metisCommandString.c_str();

      com->fprintf(stdout, "\n%s\n", metisCommandChar);
      if (!(shell = popen(metisCommandChar, "r"))) {
        fprintf(stderr, " *** Error: attempt to use external METIS executable (%s) failed!\n", ioData->input.metis);
        exit(-1);
      } else {
        com->fprintf(stdout, "\n ... Calling external METIS executable (%s) ...\n", ioData->input.metis);
      }

      char buff[512];
      while (fgets(buff, sizeof(buff), shell)!=NULL){
        com->fprintf(stdout, "%s", buff);
      }
      pclose(shell);

      std::string decompositionPathString(topFilePath);
      decompositionPathString += ".dec.";
      decompositionPathString += boost::lexical_cast<std::string>(ioData->input.nParts);

      // initial call to "sower -fluid"
      std::string gappyPrefixPathString(databasePrefix);
      gappyPrefixPathString += databaseName;
      gappyPrefixPathString += "/";
      gappyPrefixPathString += romFiles->gappyPrefix;
      
      std::vector<int> cpuMaps;
      int nCores = 1;
      while (nCores < ioData->input.nParts) {
        cpuMaps.push_back(nCores);
        nCores *= 2;
      }
      cpuMaps.push_back(ioData->input.nParts);

      std::string sowerCommandString = ioData->input.sower;
      sowerCommandString += " -fluid -mesh ";
      sowerCommandString += topFilePath;
      sowerCommandString += " -dec ";
      sowerCommandString += decompositionPathString;
      std::vector<int>::iterator it;
      for (it = cpuMaps.begin(); it != cpuMaps.end(); ++it) { 
        sowerCommandString += " -cpu ";
        sowerCommandString += boost::lexical_cast<std::string>(*it);
      }
      sowerCommandString += " -output ";
      sowerCommandString += gappyPrefixPathString;
      sowerCommandString += " -cluster ";
      sowerCommandString += boost::lexical_cast<std::string>(ioData->input.nParts);
      const char *sowerCommandChar = sowerCommandString.c_str();

      com->fprintf(stdout, "\n%s\n", sowerCommandChar);
      if (!(shell = popen(sowerCommandChar, "r"))) {
        fprintf(stderr, " *** Error: attempt to use external SOWER executable (%s) failed!\n", ioData->input.sower);
        exit(-1);
      } else {
        com->fprintf(stdout, "\n ... Calling external SOWER executable (%s) ...\n", ioData->input.sower);
      }

      while (fgets(buff, sizeof(buff), shell)!=NULL){
        com->fprintf(stdout, "%s", buff);
      }
      pclose(shell);

      std::string meshPathString(gappyPrefixPathString);
      meshPathString += ".msh";

      std::string connectivityPathString(gappyPrefixPathString);
      connectivityPathString += ".con";
    
      delete [] topFilePath;

      // now for all of the "sower -fluid -split" calls
    
      char *sampledCentersPath = NULL;
      determinePath(sampledCentersName, -1, sampledCentersPath);
      callSowerSplit(meshPathString, connectivityPathString, sampledCentersPath);
      delete [] sampledCentersPath;
    
      char *sampledSolutionPath = NULL;
      determinePath(sampledSolutionName, -1, sampledSolutionPath);
      callSowerSplit(meshPathString, connectivityPathString, sampledSolutionPath);
      delete [] sampledSolutionPath;
     
      char *sampledMatchStatePath = NULL;
      determinePath(sampledMatchStateName, -1, sampledMatchStatePath);
      callSowerSplit(meshPathString, connectivityPathString, sampledMatchStatePath);
      delete [] sampledMatchStatePath;

      char *sampledMultiSolutionsPath = NULL;
      determinePath(sampledMultiSolutionsName, -1, sampledMultiSolutionsPath);
      callSowerSplit(meshPathString, connectivityPathString, sampledMultiSolutionsPath);
      delete [] sampledMultiSolutionsPath;
    
      char *sampledWallDistPath = NULL;
      determinePath(sampledWallDistName, -1, sampledWallDistPath);
      callSowerSplit(meshPathString, connectivityPathString, sampledWallDistPath);
      delete [] sampledWallDistPath;

      char *sampledDisplacementPath = NULL;
      determinePath(sampledDisplacementName, -1, sampledDisplacementPath);
      callSowerSplit(meshPathString, connectivityPathString, sampledDisplacementPath);
      delete [] sampledDisplacementPath;

      char *sampledShapeDerivativePath = NULL;
      determinePath(sampledShapeDerivativeName, -1, sampledShapeDerivativePath);
      callSowerSplit(meshPathString, connectivityPathString, sampledShapeDerivativePath);
      delete [] sampledShapeDerivativePath;
    
      char *approxMetricStateLowRankPath = NULL;
      determinePath(approxMetricStateLowRankName, -1, approxMetricStateLowRankPath);
      callSowerSplit(meshPathString, connectivityPathString, approxMetricStateLowRankPath);
      delete [] approxMetricStateLowRankPath;

      for (int iCluster=0; iCluster<nClusters; ++iCluster) {
        char *gappyResidualPath = NULL;
        determinePath(gappyResidualName, iCluster, gappyResidualPath);
        callSowerSplit(meshPathString, connectivityPathString, gappyResidualPath);
        delete [] gappyResidualPath;
    
        char *gappyJacActionPath = NULL;
        determinePath(gappyJacActionName, iCluster, gappyJacActionPath);
        callSowerSplit(meshPathString, connectivityPathString, gappyJacActionPath);
        delete [] gappyJacActionPath;
    
        char *sampledStateBasisPath = NULL;
        determinePath(sampledStateBasisName, iCluster, sampledStateBasisPath);
        callSowerSplit(meshPathString, connectivityPathString, sampledStateBasisPath);
        delete [] sampledStateBasisPath;
    
        char *sampledRefStatePath = NULL;
        determinePath(sampledRefStateName, iCluster, sampledRefStatePath);
        callSowerSplit(meshPathString, connectivityPathString, sampledRefStatePath);
        delete [] sampledRefStatePath;

        char *approxMetricNonlinearPath = NULL;
        determinePath(approxMetricNonlinearName, iCluster, approxMetricNonlinearPath);
        callSowerSplit(meshPathString, connectivityPathString, approxMetricNonlinearPath);
        delete [] approxMetricNonlinearPath;
   
      }

    }
  }
  
  com->barrier();

}

//------------------------------------------------------------------------------
template<int dim>
void NonlinearRom<dim>::callSowerSplit(std::string meshPath, std::string conPath, char* result) {

    if (FILE *test = fopen(result, "r")) {
        fclose(test);
    } else {
        return;
    } 

    std::string sowerCommandString = ioData->input.sower;
    sowerCommandString += " -fluid -split -mesh ";
    sowerCommandString += meshPath;
    sowerCommandString += " -con ";
    sowerCommandString += conPath;
    sowerCommandString += " -result ";
    sowerCommandString += result;
    sowerCommandString += " -ascii -output ";
    sowerCommandString += result;
    const char *sowerCommandChar = sowerCommandString.c_str();

    com->fprintf(stdout, "\n%s\n", sowerCommandChar);
 
   FILE *shell;
    if (!(shell = popen(sowerCommandChar, "r"))) {
      fprintf(stderr, " *** Error: attempt to use external SOWER executable (%s) failed!\n", ioData->input.sower);
      exit(-1);
    } else {
      com->fprintf(stdout, " ... Calling external SOWER executable (%s) ...\n", ioData->input.sower);
    }

    char buff[512];
    while (fgets(buff, sizeof(buff), shell)!=NULL){
      com->fprintf(stdout, "%s", buff);
    }
    pclose(shell);

}

//----------------------------------------------

template<int dim>
void NonlinearRom<dim>::printDebug(int iDebug) {

    com->fprintf(stderr," ... Debugging: %d ...\n",iDebug);
    sleep(1);

}

//----------------------------------------------

template<int dim>
void NonlinearRom<dim>::formInterpolatedInitialCondition(DistSVec<double,dim> *U, IoData &iod) {

  // form the path to multi solutions
  *U = 0.0;
  char *sampledMultiSolutionsPath = NULL;
  if ((ioData->problem.alltype == ProblemData::_STEADY_NONLINEAR_ROM_POST_ 
           || ioData->problem.alltype == ProblemData::_UNSTEADY_NONLINEAR_ROM_POST_) 
           && strcmp(surfaceRefStateName,"")!=0) {
    determinePath(surfaceMultiSolutionsName, -1, sampledMultiSolutionsPath);
  } else {
    determinePath(sampledMultiSolutionsName, -1, sampledMultiSolutionsPath);
  }

  DistSVec<double, dim> solVec(domain.getNodeDistInfo());

  for (int iData=0; iData<interpWeightsForMultiIC.size(); ++iData) {
    // add this vector's contribution to U
    domain.readVectorFromFile(sampledMultiSolutionsPath, iData, 0, solVec);
    *U += interpWeightsForMultiIC[iData]*solVec;
  }

  delete [] sampledMultiSolutionsPath; 

}


