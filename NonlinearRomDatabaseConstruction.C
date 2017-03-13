#include <NonlinearRomDatabaseConstruction.h>
#include <Modal.h>
#include <TsInput.h>
#include <cmath>
//#include <time.h>
#include <algorithm>
#include <sys/time.h>
#include <algorithm>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <map>
#include <set>
//#include <random>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#ifdef DO_MODAL
#include <arpack++/include/ardsmat.h>
#include <arpack++/include/ardssym.h>
#endif

using std::stable_sort;

template<int dim> 
NonlinearRomDatabaseConstruction<dim>::NonlinearRomDatabaseConstruction(Communicator* _com, IoData& _ioData, Domain& _domain, GeoSource& _geoSource)  : 
  NonlinearRomOnlineII<dim>(_com, _ioData, _domain), geoSource(_geoSource)
{ 
  // ioData->example, com->example, this->domain.example
  nSnapshotFiles = 0;

  projErrorLog = NULL;
  initialCondition = NULL; 

  preproForArbitraryUniformIC = false;
  preproForInterpolatedIC = false;

  robConstruction = &(this->ioData->romOffline.rob);
  projError = &(robConstruction->relativeProjectionError);

  
}

//----------------------------------------------------------------------------------

template<int dim> 
NonlinearRomDatabaseConstruction<dim>::~NonlinearRomDatabaseConstruction() 
{
  if (initialCondition) delete initialCondition;
  initialCondition = NULL;
}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::constructDatabase() {

// The functions kmeans(), localPod(), and localRelProjError() are called here, depending on inputs.
  double tOffline = this->timer->getTime();

  // clustering
  if (robConstruction->clustering.useExistingClusters == ClusteringData::USE_EXISTING_CLUSTERS_FALSE) {

    double tRead = this->timer->getTime();
    nSnapshotFiles = this->readSnapshotFiles("state", false);
    this->timer->addReadSnapshotFileTime(tRead);

    double tCluster = this->timer->getTime();
    if (robConstruction->clustering.clusteringAlgorithm == ClusteringData::K_MEANS_WITH_BOUNDS) {
      kmeansWithBounds();
    } else {
      kmeans();
    }
    this->timer->addClusteringTime(tCluster);

    // option to form 2D representation of state data using multi-dimensional scaling (for visualization)
    if (robConstruction->clustering.computeMDS == ClusteringData::COMPUTE_MDS_TRUE) {
      double tMDS = this->timer->getTime();
      computeClassicalMultiDimensionalScaling();
      this->timer->addMDSTime(tMDS);
    }

    this->outputClusteredSnapshots("state");

    // for snapshot collection method 0 (all snapshots from FOM)
    if (strcmp(this->ioData->input.residualSnapFile,"")!=0) placeNonStateSnapshotsInClusters("residual");
    if (strcmp(this->ioData->input.krylovSnapFile,"")!=0) placeNonStateSnapshotsInClusters("krylov");

    if (strcmp(this->ioData->input.sensitivitySnapFile,"")!=0 && 
         strcmp(this->sensitivityClusterName,"")!=0) {
       int nSensitivityFiles = this->readSnapshotFiles("sensitivity", false);
       this->outputClusteredSnapshots("sensitivity");
    } 
    
  }

  // local POD
  if (robConstruction->state.dataCompression.computePOD) localPod("state");
  if (robConstruction->krylov.dataCompression.computePOD) localPod("krylov");
  if (robConstruction->sensitivity.dataCompression.computePOD) localPod("sensitivity");
  if (robConstruction->residual.dataCompression.computePOD) localPod("residual");
  if (robConstruction->jacAction.dataCompression.computePOD) localPod("jacAction");

  // projection error
  if (projError->relProjError!=RelativeProjectionErrorData::REL_PROJ_ERROR_OFF) localRelProjError();
  if ((projError->relProjError!=RelativeProjectionErrorData::REL_PROJ_ERROR_OFF) &&
      (projError->sweepFreq > 0)) localRelProjErrorSweep();

  // store ROBs in memory to avoid excessive IO time
  if (robConstruction->storeAllClusters && (robConstruction->basisUpdates.preprocessForNoUpdates ||
      robConstruction->basisUpdates.preprocessForProjections ||
      robConstruction->basisUpdates.preprocessForExactUpdates ||
      robConstruction->basisUpdates.preprocessForApproxUpdates)) this->readAllClusteredOfflineQuantities();

  // preprocessing for fast distance calculations (not currently supported for simple updates)
  if (robConstruction->basisUpdates.preprocessForNoUpdates ||
      robConstruction->basisUpdates.preprocessForProjections || 
      robConstruction->basisUpdates.preprocessForExactUpdates || 
      robConstruction->basisUpdates.preprocessForApproxUpdates) preprocessForDistanceComparisons();

  // preprocessing for exact basis updates 
  // (data for simple updates is ouput automatically; data for approx updates is output in GNAT preprocessing)
  if (robConstruction->basisUpdates.preprocessForExactUpdates ||
      robConstruction->basisUpdates.preprocessForProjections) preprocessForExactBasisUpdates();

  this->timer->addTotalOfflineTime(tOffline);

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::initializeClusterCenters() {

    if (this->clusterCenters) delete (this->clusterCenters);
    this->clusterCenters = new VecSet< DistSVec<double, dim> >((this->nClusters), this->domain.getNodeDistInfo());

    if (!(strcmp(this->ioData->input.initialClusterCentersFile,"")==0)) { 
      // use user-specified snapshots as initial cluster centers
      int nFiles = this->readSnapshotFiles("initialClusterCenters", false);

    } else {
      // pick random initial centers (use shuffle algorithm to ensure no duplicates)

      int nTotSnaps = this->snap->numVectors();
      int* shuffle = new int[nTotSnaps];     
      int kMeansRandSeed = robConstruction->clustering.kMeansRandSeed;
   
      if (this->com->cpuNum()==0) { 
        int randSeed;
    
        if (kMeansRandSeed == -1) {
          randSeed = time(NULL);
        } else {
          randSeed = kMeansRandSeed;
        }
     
        srand(randSeed);
     
        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          shuffle[iSnap] = iSnap; 
        }
    
        for (int iSnap=0; iSnap<(this->nClusters); iSnap++) { // only need to shuffle first nClusters snapshots
          int randPosition = iSnap + (rand() % (nTotSnaps-iSnap));
          int temp = shuffle[iSnap];
          shuffle[iSnap] = shuffle[randPosition];
          shuffle[randPosition] = temp;
        }
      }
    
      this->com->broadcast(nTotSnaps, shuffle, 0);
    
      for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
        (*(this->clusterCenters))[iCluster]=(*(this->snap))[shuffle[iCluster]];
      }
    
      delete [] shuffle;  
   }

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::placeNonStateSnapshotsInClusters(const char* snapType) {

  char *vecFile;

  if (strcmp(snapType,"residual")==0) {
    vecFile = new char[strlen(this->ioData->input.prefix) + strlen(this->ioData->input.residualSnapFile) + 1];
    sprintf(vecFile, "%s%s", this->ioData->input.prefix, this->ioData->input.residualSnapFile);
  } else if (strcmp(snapType,"krylov")==0) {
    vecFile = new char[strlen(this->ioData->input.prefix) + strlen(this->ioData->input.krylovSnapFile) + 1];
    sprintf(vecFile, "%s%s", this->ioData->input.prefix, this->ioData->input.krylovSnapFile);
  } else {
    this->com->fprintf(stderr, "*** Error: unexpected snapshot type %s\n", snapType);
    exit (-1);
  }

  FILE *inFP = fopen(vecFile, "r");
  if (!inFP)  {
    this->com->fprintf(stderr, "*** Error: No snapshots FILES in %s\n", vecFile);
    exit (-1);
  }

  int nData, _n;
  _n = fscanf(inFP, "%d",&nData);
  this->com->fprintf(stdout, "Reading snapshots from %d files \n",nData);

  if (nData!=nSnapshotFiles) {
    this->com->fprintf(stderr, "*** Error: incorrect number of files listed in %s (these files must correspond to the state snapshot files)\n", vecFile);
    exit(-1);
  }

  delete [] vecFile;
  vecFile = NULL;

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

  // read snapshot command file
  for (int iData = 0; iData < nData; ++iData){
    _n = fscanf(inFP, "%s %d %d %d %lf", snapFile1,&iStart,&iEnd,&iFreq,&weight);
    if (iStart < 1) iStart = 1;
    if (iEnd < 0) iEnd = 0;
    if (iFreq < 1) iFreq = 1;
    if (weight==0.0) weight = 1.0;
    //numSnaps[iData] = nSnap;
    strcpy(snapFile[iData],snapFile1);
    startSnaps[iData] = iStart - 1;
    endSnaps[iData] = iEnd;
    sampleFreq[iData] = iFreq;
    snapWeight[iData] = weight;
    this->com->fprintf(stdout, " ... Reading snapshots from %s \n", snapFile[iData]);
  }

  // compute the total number of snapshots
  int nTotSnaps = 0;
  int dummyStep = 0;
  double dummyTag = 0.0;
  for (int iData = 0; iData < nData; ++iData) {
    bool status = this->domain.template readTagFromFile<double, dim>(snapFile[iData], dummyStep, &dummyTag, &(numSnaps[iData]));
    if (!status) {
      this->com->fprintf(stdout, "*** Error: could not read snapshot file %s \n", snapFile[iData]);
    }
    if ((endSnaps[iData]==0) || (endSnaps[iData]>numSnaps[iData]))
      endSnaps[iData]=numSnaps[iData];
    for (int iSnap = startSnaps[iData]; iSnap<endSnaps[iData]; ++iSnap) {
      if (iSnap % sampleFreq[iData] == 0) {
        ++nTotSnaps;
      }
    }
  }

  if (nTotSnaps==0) {
    this->com->fprintf(stdout, "*** Error: expecting to read zero snapshots (be sure to include start, finish, freq, and weight in snapshot file %s) \n", vecFile);
    exit(-1);
  }

  DistSVec<double, dim>* snapBuf = new DistSVec<double, dim>(this->domain.getNodeDistInfo());
  *snapBuf = 0.0;
  double tag;
  bool status;

  this->snapsInCluster = new int[this->nClusters]; 
  for (int iCluster=0;iCluster<(this->nClusters);++iCluster)
    this->snapsInCluster[iCluster] = 0;

  this->initializeClusteredOutputs();

  for (int iData=0; iData < nData; ++iData){
    // read in Snapshot Vectors
    for (int iSnap = startSnaps[iData]; iSnap<endSnaps[iData]; ++iSnap) {
      if (iSnap % sampleFreq[iData] == 0) { //TODO ignore 
        // snapshot must be between startSnaps and endSnaps, and a multiple of sampleFreq. 
        if (this->duplicateSnaps) {
          status = this->domain.readVectorFromFile(snapFile[iData], iSnap, &tag, *snapBuf);        
          if (snapWeight[iData]) *snapBuf *= snapWeight[iData]; //CBM--check
        } else { // only need the tag
          int dummyNumSteps;
          status = this->domain.template readTagFromFile<double, dim>(snapFile[iData], iSnap, &tag, &dummyNumSteps);
        }
        // find associated state
        int stateSnap = -1;
        for (int iStateSnap=0; iStateSnap<this->stateSnapsFromFile[iData]; ++iStateSnap) {
          if (iStateSnap==(this->stateSnapsFromFile[iData]-1)) {
            stateSnap = iStateSnap;
          } else if (tag<this->stateSnapshotTags[iData][iStateSnap]) {
            // do nothing
          } else if ((tag>=this->stateSnapshotTags[iData][iStateSnap]) && (tag<this->stateSnapshotTags[iData][iStateSnap+1])) {
            stateSnap = iStateSnap;
            break;
          }
        }
        // store snapshot in all applicaple clusters
        if (strcmp(snapType,"residual")==0) {
          for (int iCluster=0;iCluster<(this->nClusters);++iCluster) {
            if (stateSnapshotClustersAfterOverlap[iData][stateSnap][iCluster]) { 
              this->com->fprintf(stdout, " ... training simulation #%d: residual snapshot #%d matched with state snapshot #%d; writing to cluster %d\n",
                iData,iSnap,stateSnap,iCluster); 
              this->writeClusteredBinaryVectors(iCluster, snapBuf, NULL, NULL, snapFile[iData], iSnap);
              ++(this->snapsInCluster[iCluster]);
            }
          }
        } else if (strcmp(snapType,"krylov")==0) {
          for (int iCluster=0;iCluster<(this->nClusters);++iCluster) {
            if (stateSnapshotClustersAfterOverlap[iData][stateSnap][iCluster]) {
              this->com->fprintf(stdout, " ... training simulation #%d: krylov snapshot #%d matched with state snapshot #%d; writing to cluster %d\n",
                iData,iSnap,stateSnap,iCluster);
              this->writeClusteredBinaryVectors(iCluster, NULL, NULL, snapBuf, snapFile[iData], iSnap);
              ++(this->snapsInCluster[iCluster]);
            }
          }
        }
      }
    }
  }

  this->com->fprintf(stdout, "\n"); 
  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
    this->com->fprintf(stdout, " ... %d %s vectors placed in cluster %d\n",this->snapsInCluster[iCluster], snapType, iCluster);
  }
  this->com->fprintf(stdout, "\n");


  delete [] (this->snapsInCluster);
  this->snapsInCluster = NULL;

  delete snapBuf;
  snapBuf = NULL;

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

}


//----------------------------------------------------------------------------------

template<int dim>
double NonlinearRomDatabaseConstruction<dim>::calcResidual(VecSet< DistSVec<double, dim> > &centers, VecSet< DistSVec<double, dim> > &centersOld) {

// Calculates the residual for the kmeans clustering.  The clustering converges when the cluster centers stop moving.

  double norm;
  double maxNorm = 0.0;

  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
    norm = this->distanceFull( centers[iCluster], centersOld[iCluster]);
    if (norm > maxNorm) maxNorm = norm;
  }

  return maxNorm;

}


//----------------------------------------------------------------------------------

// this struct is used in the kmeans algorithm
struct sortStruct {
  int snapIndex; // snapshot #
  double dist; // distance to second closest cluster

  bool operator<(const sortStruct& a) const {
    return dist < a.dist;  
  }
};

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::kmeans() {

  this->com->fprintf(stdout, "\nUsing K-Means algorithm to cluster snapshots\n");

  // parameters that control the kmeans clustering
  int minClusterSize = robConstruction->clustering.minClusterSize;
  double percentOverlap = robConstruction->clustering.percentOverlap;
  double kMeansTol = robConstruction->clustering.kMeansTol;

  int nTotSnaps = this->snap->numVectors();

  (this->clusterIndex) = new int[nTotSnaps];
  VecSet< DistSVec<double, dim> > clusterCentersOld((this->nClusters), this->domain.getNodeDistInfo());
  initializeClusterCenters();

  int iterMax = robConstruction->clustering.maxIter;  // max number of kmeans iterations
  int iterMaxAggressive = robConstruction->clustering.maxIterAggressive;  // number of aggressive kmeans iterations to use before switching to a more robust single update scheme

  int index1 = 0;
  int index2 = 0;

  (this->snapsInCluster) = new int[(this->nClusters)];
  for (int iCluster=0; iCluster<this->nClusters; iCluster++)
    (this->snapsInCluster)[iCluster] = 0;

  // k-means algorithm
  int iter=0;
  double residual=1.0;
  while (((residual > kMeansTol) || (iter == 0)) && (iter<iterMax)) {
    
    this->com->fprintf(stdout, "Clustering iteration #%d \n", iter+1);

    double prevResidual = residual;
  
    if ((iter<iterMaxAggressive) || (iter==0)) {
      this->com->fprintf(stdout, " ... updating all snapshots simultaneously\n");
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        this->closestCenterFull( (*(this->snap))[iSnap], &index1, &index2);
        (this->clusterIndex)[iSnap]=index1;
      }
  
      for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
        (this->snapsInCluster)[iCluster]=0;
        clusterCentersOld[iCluster] = (*(this->clusterCenters))[iCluster];
        (*(this->clusterCenters))[iCluster] = 0.0;
      }

      //KYLE_HERE
      DistSVec<double, dim>* tmpDistVec = new DistSVec<double, dim>(this->domain.getNodeDistInfo());
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        *tmpDistVec = (*(this->snap))[iSnap];
        // if using angles, cluster center is average of normalized snapshots
        if (!this->euclideanDistances) {
          double normalize = tmpDistVec->norm();
          //if (*tmpDistVec * clusterCentersOld[(this->clusterIndex)[iSnap]] < 0) normalize *= -1.0;
          if (normalize > 0) normalize = 1.0/normalize;
          *tmpDistVec *= normalize;
        }
        (*(this->clusterCenters))[(this->clusterIndex)[iSnap]] += *tmpDistVec;
        ++((this->snapsInCluster)[(this->clusterIndex)[iSnap]]);  
      }
      delete tmpDistVec;
      tmpDistVec = NULL;

      for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
        this->com->fprintf(stdout, " ... cluster %d has %d snaps \n", iCluster, (this->snapsInCluster)[iCluster]);
        double snapsInClusterInv = 0;
        if ((this->snapsInCluster)[iCluster] != 0) snapsInClusterInv = 1.0/double((this->snapsInCluster)[iCluster]);
        (*(this->clusterCenters))[iCluster] *= snapsInClusterInv;
      }

    residual = calcResidual(*(this->clusterCenters), clusterCentersOld);
  
    } else { // algorithm is likely stuck in a cycle -- begin updating updating one at a time

      this->com->fprintf(stdout, " ... updating one snapshot at a time\n");

      for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
        clusterCentersOld[iCluster] = (*(this->clusterCenters))[iCluster];
      }
    
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        int oldIndex = (this->clusterIndex)[iSnap];
        int newIndex;

        this->closestCenterFull( (*(this->snap))[iSnap], &newIndex);

        // if using angles, cluster center is average of normalized snapshots
        double normalize = 1.0;
        if (!this->euclideanDistances) {
          normalize = (*(this->snap))[iSnap].norm();
          if (normalize > 0) normalize = 1.0/normalize;
        }

        //KYLE_HERE
        if (oldIndex != newIndex) {
          (this->clusterIndex)[iSnap] = newIndex;
          (*(this->clusterCenters))[oldIndex] *= (double((this->snapsInCluster)[oldIndex])/(double((this->snapsInCluster)[oldIndex])-1.0));
          (*(this->clusterCenters))[oldIndex] -= (*(this->snap))[iSnap]*(normalize)*(1.0/(double((this->snapsInCluster)[oldIndex])-1.0));
          (*(this->clusterCenters))[newIndex] *= (double((this->snapsInCluster)[newIndex])/(double((this->snapsInCluster)[newIndex])+1.0));
          (*(this->clusterCenters))[newIndex] += (*(this->snap))[iSnap]*(normalize)*(1.0/(double((this->snapsInCluster)[newIndex])+1.0));
          --((this->snapsInCluster)[oldIndex]);
          ++((this->snapsInCluster)[newIndex]);
        } 
      }

      for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
        this->com->fprintf(stdout, " ... cluster %d has %d snaps \n", iCluster, (this->snapsInCluster)[iCluster]);
      }

      residual = calcResidual(*(this->clusterCenters), clusterCentersOld);

    //  if (pow((prevResidual - residual),2)<pow(kMeansTol,2)) {
    //    this->com->fprintf(stderr, "*** Warning: Clustering algorithm is stuck. Exiting now.\n");
    //    break;
    //  }
    }
    this->com->fprintf(stdout, " ... absolute residual = %e (tolerance is set to %e)\n", residual, kMeansTol);
    ++iter;
  }


// after clustering, assimilate small clusters
  int emptyClusterCount = 0;
  if (minClusterSize>nTotSnaps) minClusterSize = nTotSnaps;

  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
    if ((this->snapsInCluster)[iCluster]==0) {
      // if no snapshots in this cluster, skip
      ++emptyClusterCount;
    } else {
    // if smaller than tolerance, add to nearest cluster
      if ((this->snapsInCluster)[iCluster]<minClusterSize) {
        this->com->fprintf(stderr, "*** Warning: combining small cluster with nearest neighbor\n");
        int index1 = 0; // closest center (should be itself)
        int index2 = 0; // second closest center (the one we want)

        this->closestCenterFull((*(this->clusterCenters))[iCluster], &index1, &index2);
        (*(this->clusterCenters))[index2] =   (*(this->clusterCenters))[index2]*(double((this->snapsInCluster)[index2])/
                                                double((this->snapsInCluster)[index2]+(this->snapsInCluster)[iCluster]))
                                            + (*(this->clusterCenters))[iCluster]*(double((this->snapsInCluster)[iCluster])/
                                                double((this->snapsInCluster)[index2]+(this->snapsInCluster)[iCluster]));
 
        (*(this->clusterCenters))[iCluster] = 0.0;
 
        (this->snapsInCluster)[index2] = (this->snapsInCluster)[index2] + (this->snapsInCluster)[iCluster];
        (this->snapsInCluster)[iCluster] = 0;  

        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          if ((this->clusterIndex)[iSnap]==iCluster)  (this->clusterIndex)[iSnap]=index2;
        }
        
        ++emptyClusterCount;
      }
    }
  }

// remove any empty clusters, renumber existing clusters

  if (emptyClusterCount>0) {

    this->com->fprintf(stderr, "*** Warning: Deleting %d empty clusters\n", emptyClusterCount);

    VecSet< DistSVec<double, dim> >* clusterCentersNew =  new VecSet< DistSVec<double, dim> >(((this->nClusters)-emptyClusterCount), this->domain.getNodeDistInfo()); 
    int* snapsInClusterNew = new int[((this->nClusters)-emptyClusterCount)];

    int renumberedIndices[(this->nClusters)];
    int clusterCount = 0;

    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
      if ((this->snapsInCluster)[iCluster]==0) {
        renumberedIndices[iCluster] = -1;
      } else {
        renumberedIndices[iCluster] = clusterCount; 
        (*clusterCentersNew)[clusterCount]=(*(this->clusterCenters))[iCluster];
        snapsInClusterNew[clusterCount]=(this->snapsInCluster)[iCluster];
        ++clusterCount;
      }
    }

    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
      (this->clusterIndex)[iSnap] = renumberedIndices[(this->clusterIndex)[iSnap]];
    }

    delete (this->clusterCenters);
    (this->clusterCenters) = clusterCentersNew;
    clusterCentersNew = NULL;

    delete [] (this->snapsInCluster);
    (this->snapsInCluster) = snapsInClusterNew; 
    snapsInClusterNew = NULL;  

    (this->nClusters) = (this->nClusters) - emptyClusterCount;
  }
  

// find snapshot closest to each center
  (this->nearestSnapsToCenters) = new VecSet< DistSVec<double, dim> >((this->nClusters), this->domain.getNodeDistInfo());
  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
      
    int sortSize = (this->snapsInCluster)[iCluster]; 
    sortStruct* snapDist = new sortStruct[sortSize]; // this struct was defined earlier in this file 

    int snapCount = 0;
    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
      if ((this->clusterIndex)[iSnap]==iCluster) { // only need to sort snapshots that are inside the current cluster
        snapDist[snapCount].snapIndex = iSnap;
        snapDist[snapCount].dist = this->distanceFull((*(this->snap))[iSnap],(*(this->clusterCenters))[iCluster]);
        ++snapCount;
      }
    }

    sort(snapDist, snapDist+sortSize);

    (*(this->nearestSnapsToCenters))[iCluster] = (*(this->snap))[snapDist[0].snapIndex];

    delete [] snapDist;
    snapDist = NULL;
  }

 
// add overlap if required
  if ((percentOverlap>0) && (this->nClusters>1)) {

    this->com->fprintf(stdout, "Adding additional vectors to clusters (PercentOverlap = %2.1f%%)\n", percentOverlap);

    int index1 = 0;
    int index2 = 0;
    double dist1 = 0;
    double dist2 = 0;

    //determine which clusters are neighbors
    //approach: if a cluster is second closest to a snapshot in another cluster, then those two clusters are neighbors
    (this->clusterNeighbors) = new int*[(this->nClusters)];
    (this->clusterNeighborsCount) = new int[(this->nClusters)];

    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
      (this->clusterNeighborsCount)[iCluster] = 0;
      (this->clusterNeighbors)[iCluster] = new int[1];
    }

    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        if ((this->clusterIndex)[iSnap]==iCluster) { // only consider snapshots that are inside of the current cluster
          this->closestCenterFull( (*(this->snap))[iSnap], &index1, &index2, &dist1, &dist2);
          //add the second closest cluster as a neighbor of current cluster (if this is the first time we've found it)
          bool unique = true;
          for (int iNeighbor=0; iNeighbor<(this->clusterNeighborsCount)[iCluster]; ++iNeighbor) {
            if (index2 == (this->clusterNeighbors)[iCluster][iNeighbor]) unique = false;
          }
          if (unique) {
            int* clusterNeighborsNew = new int[((this->clusterNeighborsCount)[iCluster])+1];
            for (int iNeighbor=0; iNeighbor<(this->clusterNeighborsCount)[iCluster]; ++iNeighbor) {
              clusterNeighborsNew[iNeighbor] = (this->clusterNeighbors)[iCluster][iNeighbor];
            }
 
           clusterNeighborsNew[((this->clusterNeighborsCount)[iCluster])] = index2;
            delete (this->clusterNeighbors)[iCluster];
            (this->clusterNeighbors)[iCluster] = clusterNeighborsNew;
            clusterNeighborsNew = NULL;
            ++(this->clusterNeighborsCount)[iCluster];
         }
          //also add the current cluster as a neighbor of second closest cluster
          unique = true;
          for (int iNeighbor=0; iNeighbor<(this->clusterNeighborsCount)[index2]; ++iNeighbor) {
            if (iCluster == (this->clusterNeighbors)[index2][iNeighbor]) unique = false;
          }
          if (unique) {
            int* clusterNeighborsNew = new int[((this->clusterNeighborsCount)[index2])+1];
            for (int iNeighbor=0; iNeighbor<(this->clusterNeighborsCount)[index2]; ++iNeighbor) {
              clusterNeighborsNew[iNeighbor] = (this->clusterNeighbors)[index2][iNeighbor];
            }
 
            clusterNeighborsNew[((this->clusterNeighborsCount)[index2])] = iCluster;
            delete (this->clusterNeighbors)[index2];
            (this->clusterNeighbors)[index2] = clusterNeighborsNew;
            clusterNeighborsNew = NULL;
            ++(this->clusterNeighborsCount)[index2];
         }
        }
      }
    }

    // if clustering increments, need to actually cluster increments by the preceding increment 
    // Fairly simply, just need to shift all cluster associations.
    // For the time being , I won't shift the cluster association of any IC... 
    // (This is still consistent, provided that all ROBs are tested for the IC of the online simulation)
    if (this->ioData->romDatabase.avgIncrementalStates==NonlinearRomFileSystemData::AVG_INCREMENTAL_STATES_TRUE) {

      int nSnapsHandled = 0;
      std::vector<int> newSnapsInCluster(this->nClusters,0);
      std::vector<int> newClusterIndex(nTotSnaps,0);

      for (int iFile=0; iFile<nSnapshotFiles; ++iFile) {
        newClusterIndex[nSnapsHandled]=(this->clusterIndex)[nSnapsHandled];
        ++newSnapsInCluster[newClusterIndex[nSnapsHandled]];
        ++nSnapsHandled;
        for (int iSnap=1; iSnap<this->stateSnapsFromFile[iFile]; ++iSnap) {
          newClusterIndex[nSnapsHandled] = (this->clusterIndex)[nSnapsHandled-1];
          ++newSnapsInCluster[newClusterIndex[nSnapsHandled]];
          ++nSnapsHandled;
        }
      }
      

      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) (this->clusterIndex)[iSnap] = newClusterIndex[iSnap];                         

      emptyClusterCount = 0;
      for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
        (this->snapsInCluster)[iCluster] = newSnapsInCluster[iCluster];
        if ((this->snapsInCluster)[iCluster] == 0) ++emptyClusterCount;
      }

      if (emptyClusterCount>0) {
        this->com->fprintf(stderr, "*** Warning: Deleting %d empty clusters\n", emptyClusterCount);

        VecSet< DistSVec<double, dim> >* clusterCentersNew =  new VecSet< DistSVec<double, dim> >(((this->nClusters)-emptyClusterCount), this->domain.getNodeDistInfo());
        int* snapsInClusterNew = new int[((this->nClusters)-emptyClusterCount)];
    
        int renumberedIndices[(this->nClusters)];
        int clusterCount = 0;
    
        for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
          if ((this->snapsInCluster)[iCluster]==0) {
            renumberedIndices[iCluster] = -1;
          } else {
            renumberedIndices[iCluster] = clusterCount;
            (*clusterCentersNew)[clusterCount]=(*(this->clusterCenters))[iCluster];
            snapsInClusterNew[clusterCount]=(this->snapsInCluster)[iCluster];
            ++clusterCount;
          }
        }
    
        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          (this->clusterIndex)[iSnap] = renumberedIndices[(this->clusterIndex)[iSnap]];
        }
    
        delete (this->clusterCenters);
        (this->clusterCenters) = clusterCentersNew;
        clusterCentersNew = NULL;
    
        delete [] (this->snapsInCluster);
        (this->snapsInCluster) = snapsInClusterNew;
        snapsInClusterNew = NULL;
    
        (this->nClusters) = (this->nClusters) - emptyClusterCount;
      }

    }

    index1 = 0;
    index2 = 0;
    dist1 = 0;
    dist2 = 0;

    int origSnapsInCluster[(this->nClusters)];
    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) origSnapsInCluster[iCluster] = (this->snapsInCluster)[iCluster];

    //share shapshots between neighboring clusters
    (this->clusterSnapshotMap) = new int*[(this->nClusters)];
    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {  

      int numToAdd = int(ceil(double((this->snapsInCluster)[iCluster])*percentOverlap*.01/double((this->clusterNeighborsCount)[iCluster])));
      if (((this->snapsInCluster)[iCluster]+(numToAdd*(this->clusterNeighborsCount)[iCluster]))>nTotSnaps) 
        numToAdd=int(floor(double((nTotSnaps-(this->snapsInCluster)[iCluster]))/double((this->clusterNeighborsCount)[iCluster])));
      (this->clusterSnapshotMap)[iCluster] = new int[((this->snapsInCluster)[iCluster]+(numToAdd*(this->clusterNeighborsCount)[iCluster]))];

      // first add all snapshots that were originally in the cluster
      int mappedCount = 0;
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        if ((this->clusterIndex)[iSnap]==iCluster) {
          (this->clusterSnapshotMap)[iCluster][mappedCount] = iSnap;
          ++mappedCount;
        }
      }

      for (int iNeighbor=0; iNeighbor<(this->clusterNeighborsCount)[iCluster]; ++iNeighbor) {

        int sortSize = origSnapsInCluster[(this->clusterNeighbors)[iCluster][iNeighbor]]; 
        sortStruct* snapDist = new sortStruct[sortSize]; // this struct was defined earlier in this file 

        int snapCount = 0;

        if (this->euclideanDistances) {

          DistSVec<double, dim>* tmpDistVec = new DistSVec<double, dim>(this->domain.getNodeDistInfo());
          *tmpDistVec = (*(this->clusterCenters))[iCluster] - (*(this->clusterCenters))[(this->clusterNeighbors)[iCluster][iNeighbor]];
          double offset = (pow(((*(this->clusterCenters))[iCluster]).norm(),2) - pow(((*(this->clusterCenters))[(this->clusterNeighbors)[iCluster][iNeighbor]]).norm(),2))*(-0.5);
          for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
            if ((this->clusterIndex)[iSnap]==(this->clusterNeighbors)[iCluster][iNeighbor]) { // don't need to sort snapshots belonging to the current cluster
              snapDist[snapCount].snapIndex = iSnap;
              snapDist[snapCount].dist = abs( (*tmpDistVec) * (*(this->snap))[iSnap] + offset );  //don't need the abs? KMW
              ++snapCount;
            }
          }
          delete tmpDistVec;
          tmpDistVec = NULL;

        } else { // angles

          DistSVec<double, dim> projVec1(this->domain.getNodeDistInfo());
          DistSVec<double, dim> projVec2(this->domain.getNodeDistInfo());

          projVec1 = (*(this->clusterCenters))[iCluster];
          projVec1 *= 1.0/projVec1.norm(); // norm can't be zero
          projVec2 = (*(this->clusterCenters))[(this->clusterNeighbors)[iCluster][iNeighbor]];
          projVec2 = projVec2 - (projVec1*projVec2)*projVec1;
          projVec2 *= 1.0/projVec2.norm(); // again, norm can't be zero

          DistSVec<double, dim>* tmpDistVec = new DistSVec<double, dim>(this->domain.getNodeDistInfo());

          for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
            if ((this->clusterIndex)[iSnap]==(this->clusterNeighbors)[iCluster][iNeighbor]) { // don't need to sort snapshots belonging to the current cluster
              snapDist[snapCount].snapIndex = iSnap;
              *tmpDistVec = projVec1 * (projVec1*((*(this->snap))[iSnap]));
              *tmpDistVec += projVec2 * (projVec2*((*(this->snap))[iSnap]));

              snapDist[snapCount].dist = this->distanceFull(*tmpDistVec, (*(this->clusterCenters))[iCluster]);
              ++snapCount;
            }
          }
          delete tmpDistVec;
          tmpDistVec = NULL;

        }
        sort(snapDist, snapDist+sortSize);

        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          for (int jSnap=0; jSnap<numToAdd; ++jSnap) {
            if (snapDist[jSnap].snapIndex==iSnap) {
              (this->clusterSnapshotMap)[iCluster][mappedCount] = iSnap;
              ++mappedCount;
              ++((this->snapsInCluster)[iCluster]);
            }
          }
        }

        delete [] snapDist;
        snapDist = NULL;
      }

      this->com->fprintf(stdout, " ... cluster %d has %d snaps \n", iCluster, (this->snapsInCluster)[iCluster]);

    }
  } else { // no overlap
    (this->clusterSnapshotMap) = new int*[(this->nClusters)];
    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {  
      (this->clusterSnapshotMap)[iCluster] = new int[(this->snapsInCluster)[iCluster]];

      int mappedCount = 0;
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        if ((this->clusterIndex)[iSnap]==iCluster) {
          (this->clusterSnapshotMap)[iCluster][mappedCount] = iSnap;
          ++mappedCount;
        }
      }
    }
  }

  // create data structure needed for clustering non-state FOM information (explained in header file)
  stateSnapshotClustersAfterOverlap.resize(nSnapshotFiles);
  for (int iFile=0; iFile<nSnapshotFiles; ++iFile) {
    stateSnapshotClustersAfterOverlap[iFile].resize(this->stateSnapsFromFile[iFile]);
    for (int iSnap=0; iSnap<this->stateSnapsFromFile[iFile]; ++iSnap) {
      stateSnapshotClustersAfterOverlap[iFile][iSnap].clear();
      stateSnapshotClustersAfterOverlap[iFile][iSnap].resize(this->nClusters,false);
    }
  }

  int** snapInfo = new int*[nTotSnaps];  // temporary - for constructing stateSnapshotClustersAfterOverlap
  for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
    snapInfo[iSnap] = new int[2];
  }

  int snapCount = 0;
  for (int iFile=0; iFile<nSnapshotFiles; ++iFile) {
    for (int iSnap=0; iSnap<this->stateSnapsFromFile[iFile]; ++iSnap) {
      snapInfo[snapCount][0] = iFile;// file
      snapInfo[snapCount][1] = iSnap;// snapshot number from file
      ++snapCount;
    }
  }

  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
    for (int iSnap=0; iSnap<(this->snapsInCluster[iCluster]); ++iSnap) {
      int currentSnap = this->clusterSnapshotMap[iCluster][iSnap];
      stateSnapshotClustersAfterOverlap[snapInfo[currentSnap][0]][snapInfo[currentSnap][1]][iCluster] = true;
    }
  }

  for (int iSnap=0; iSnap<nTotSnaps; ++iSnap)
    delete [] snapInfo[iSnap];

  delete [] snapInfo;

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::kmeansWithBounds() {
// Utilizes triangle inequality bounds to skip distance calculations for states that
// couldn't have changed clusters between kmeans iterations.  This can be done with
// either tight bounds (which requires storing the cluster centers at every iteration) 
// or with loose bounds (bounds become more conservative at each iteration).
 
  bool tightBounds = false;
  bool looseBounds = false;

  this->com->fprintf(stdout, "\nUsing K-Means algorithm to cluster snapshots\n");

  switch (robConstruction->clustering.kmeansBoundType) {
    case (ClusteringData::TIGHT_BOUNDS):
      this->com->fprintf(stdout, " ... using triangle inequality bounds (computed exactly) to improve K-Means efficiency\n");
      tightBounds = true;
      break;
    case (ClusteringData::LOOSE_BOUNDS):
      this->com->fprintf(stdout, " ... using triangle inequality bounds (conservative/approx) to improve K-Means efficiency\n");
      looseBounds = true;
      break;
    default:
      this->com->fprintf(stderr, "*** ERROR: unexpected inputs in kmeansWithBounds()\n");
      exit(-1);
  }

  // parameters that control the kmeans clustering
  int minClusterSize = robConstruction->clustering.minClusterSize;
  double percentOverlap = robConstruction->clustering.percentOverlap;
  double kMeansTol = robConstruction->clustering.kMeansTol;
  int iterMax = robConstruction->clustering.maxIter;  // max number of kmeans iterations
  int iterMaxAggressive = robConstruction->clustering.maxIterAggressive;  // number of aggressive kmeans iterations to use before switching to a more robust single update scheme

  // for clustering each file separately
  VecSet< DistSVec<double, dim> >** clusterCentersAll;
  std::vector<std::vector<int> > clusterIndexAll;
  std::vector<std::vector<int> > snapsInClusterAll;

  int nFiles;
  if (robConstruction->clustering.clusterFilesSeparately) {
    nFiles = nSnapshotFiles;
    clusterCentersAll = new VecSet< DistSVec<double, dim> >*[nFiles];
    clusterIndexAll.resize(nFiles);
    snapsInClusterAll.resize(nFiles);
  } else {
    nFiles = 1;
  } 

  VecSet< DistSVec<double, dim> >* clusterCentersOld = NULL;
  VecSet< DistSVec<double, dim> >** clusterCentersLog = NULL;
  std::vector<std::vector<double> > distToCenters; // [iSnap][iCluster] distance from iSnap to iCluster (Loose and Tight)
  std::vector<std::vector<double> > uncertainty;   // [iSnap][iCluster] uncertainty of distToCenters entries (Loose and Tight) 
  std::vector<std::vector<double> > distToCentersIteration;  // [iSnap][iCluster] = kmeans iteration corresponding to distToCenters entries (Tight only)

  int nDistSkipped = 0;
  int nDistComputed = 0;
  int nTotSnaps = 0;

  VecSet< DistSVec<double, dim> >* snapshots = NULL;
  int snapCount = 0;

  for (int iFile=0; iFile<nFiles; ++iFile) {
 
    if (robConstruction->clustering.clusterFilesSeparately) {
      this->com->fprintf(stdout, "\nClustering data from snapshot file #%d \n", iFile);
      nTotSnaps = this->stateSnapsFromFile[iFile]; //set in readSnapshotFiles
      if (snapshots) delete snapshots;
      snapshots = new VecSet< DistSVec<double, dim> >(nTotSnaps, this->domain.getNodeDistInfo());
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        (*snapshots)[iSnap] = (*(this->snap))[snapCount];
        ++snapCount;
      }
    } else {
      nTotSnaps = this->snap->numVectors();
      snapshots = this->snap;
    }

    this->clusterIndex = new int[nTotSnaps];
    this->snapsInCluster = new int[(this->nClusters)];

    distToCenters.resize(nTotSnaps);
    uncertainty.resize(nTotSnaps);
    for (int iSnap=0;iSnap<nTotSnaps;++iSnap) {
      distToCenters[iSnap].clear();
      uncertainty[iSnap].clear();
      distToCenters[iSnap].resize(this->nClusters,0.0);
      uncertainty[iSnap].resize(this->nClusters,0.0);
      (this->clusterIndex)[iSnap] = -1;
    }

    if (tightBounds) { // need to store all cluster centers
      clusterCentersLog = new VecSet< DistSVec<double, dim> >*[iterMax];
      for (int iCluster=0;iCluster<this->nClusters;++iCluster) 
        clusterCentersLog[iCluster] = NULL;

      distToCentersIteration.resize(nTotSnaps);
      for (int iSnap=0;iSnap<nTotSnaps;++iSnap) {
        distToCentersIteration[iSnap].clear();
        distToCentersIteration[iSnap].resize(this->nClusters,0.0);
      }
    } else if (looseBounds) { // only need centers from previous iteration
      clusterCentersOld = new VecSet< DistSVec<double, dim> >((this->nClusters), this->domain.getNodeDistInfo());
    }

    // initialize cluster centers either randomly or manualy
    initializeClusterCenters(); 

    // start k-means algorithm
    for (int iCluster=0; iCluster<this->nClusters; iCluster++)
      (this->snapsInCluster)[iCluster] = 0;

    int iter=0;
    double residual=1.0;
    int index1 = 0;
    int index2 = 0;

    while (((residual > kMeansTol) || (iter == 0)) && (iter<iterMax)) {
      
      this->com->fprintf(stdout, "Clustering iteration #%d \n", iter+1);
  
      // set centersOld = centers
      if (tightBounds) {
        clusterCentersLog[iter] = new VecSet< DistSVec<double, dim> >((this->nClusters), this->domain.getNodeDistInfo());
        for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
          (*(clusterCentersLog[iter]))[iCluster] = (*(this->clusterCenters))[iCluster];
        }
        clusterCentersOld = clusterCentersLog[iter];
      } else {
        for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
          (*clusterCentersOld)[iCluster] = (*(this->clusterCenters))[iCluster];
         }
      }
  
      std::vector<int> clusterIndexOld;
      clusterIndexOld.resize(nTotSnaps);
      for (int iSnap=0;iSnap<nTotSnaps;++iSnap)
        clusterIndexOld[iSnap] = (this->clusterIndex)[iSnap];
  
      double prevResidual = residual;
  
      if (iter==0) {
        this->com->fprintf(stdout, " ... initializing the clusters\n");
  
        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          this->distancesToCentersFull((*snapshots)[iSnap], distToCenters[iSnap], &((this->clusterIndex)[iSnap]));
          clusterIndexOld[iSnap] = (this->clusterIndex)[iSnap];
        }
  
        nDistComputed = nTotSnaps*(this->nClusters); 
        for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
          (*(this->clusterCenters))[iCluster] = 0.0;
        }
        
        //KYLE_HERE
        DistSVec<double, dim>* tmpDistVec = new DistSVec<double, dim>(this->domain.getNodeDistInfo());
        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          *tmpDistVec = (*snapshots)[iSnap];
          // if using angles, cluster center is average of normalized snapshots
          if (!this->euclideanDistances) {
            double normalize = tmpDistVec->norm();
            if (normalize > 0) normalize = 1.0/normalize;
            *tmpDistVec *= normalize;
          }
          (*(this->clusterCenters))[(this->clusterIndex)[iSnap]] += *tmpDistVec;
          ++((this->snapsInCluster)[(this->clusterIndex)[iSnap]]);
        }
        delete tmpDistVec;
        tmpDistVec = NULL;

        for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
          double snapsInClusterInv = 0;
          if ((this->snapsInCluster)[iCluster] != 0) snapsInClusterInv = 1.0/double((this->snapsInCluster)[iCluster]);
          (*(this->clusterCenters))[iCluster] *= snapsInClusterInv;
        }
      } else {
  
        if (iter<iterMaxAggressive) {
          this->com->fprintf(stdout, " ... updating all snapshots simultaneously\n");
        } else {
          this->com->fprintf(stdout, " ... updating one snapshot at a time\n");
        }
  
        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
            if (iCluster != (this->clusterIndex)[iSnap]) {
              // test bounds to see if iSnap could possibly belong to iCluster
              if (tightBounds && (iter<iterMaxAggressive) && 
                  (0 > distToCenters[iSnap][(this->clusterIndex)[iSnap]] + distToCenters[iSnap][iCluster]
                       + uncertainty[iSnap][(this->clusterIndex)[iSnap]] - uncertainty[iSnap][iCluster])) {
                // unlikely, but in this case iSnap cannot belong to iCluster -- no need to calculate distances
                ++nDistSkipped;
                this->com->fprintf(stdout, " ... case 1 is unlikely, but apparently not impossible!\n");
              } if (0 > distToCenters[iSnap][(this->clusterIndex)[iSnap]] - distToCenters[iSnap][iCluster]
                        + uncertainty[iSnap][(this->clusterIndex)[iSnap]] + uncertainty[iSnap][iCluster]) {
                // iSnap cannot belong to iCluster -- no need to calculate distances
                ++nDistSkipped;
              } else {
                // can't rule out that iSnap belongs to iCluster; eliminate uncertainty one at a time
                if (uncertainty[iSnap][(this->clusterIndex)[iSnap]]>0) {
                  distToCenters[iSnap][(this->clusterIndex)[iSnap]] = // ...
                    this->distanceFull((*snapshots)[iSnap], (*(this->clusterCenters))[(this->clusterIndex)[iSnap]]);
                  uncertainty[iSnap][(this->clusterIndex)[iSnap]] = 0.0; 
                  ++nDistComputed;
                  if (tightBounds) distToCentersIteration[iSnap][(this->clusterIndex)[iSnap]] = iter;
                }
                if ( 0 > distToCenters[iSnap][(this->clusterIndex)[iSnap]] - distToCenters[iSnap][iCluster] 
                         + uncertainty[iSnap][iCluster]) {
                  // iSnap cannot belong to iCluster -- no further distance calculations required
                  ++nDistSkipped;
                } else {
                  if (uncertainty[iSnap][iCluster]>0) {
                    // remove remaining uncertainty by computing distToCenters[iSnap][iCluster] exactly
                    distToCenters[iSnap][iCluster] = this->distanceFull((*snapshots)[iSnap], (*(this->clusterCenters))[iCluster]);
                    uncertainty[iSnap][iCluster] = 0.0;
                    ++nDistComputed;
                    if (tightBounds) distToCentersIteration[iSnap][iCluster] = iter;
                  }
                  if ( 0 <= distToCenters[iSnap][(this->clusterIndex)[iSnap]] - distToCenters[iSnap][iCluster]) {
                    // iSnap is closer to iCluster than current cluster
                    int oldIndex = (this->clusterIndex)[iSnap];
                    int newIndex = iCluster;
                    (this->clusterIndex)[iSnap] = newIndex;
  
                    if (iter<iterMaxAggressive) {
                      // don't update centers or uncertainties until after looping through all snapshots
                    } else {
                      // this case is tricky for the tight bound.  I've decided to take a hybrid approach here where 
                      // all bounds are updated with the (worst case scenario) loose bounds when updating a single datum
                      // at a time, but in the tight bounds case the bounds are evaluated exactly after iterating through
                      // all snapshots.  In this case "tight bounds" is a bit of a misnomer because the bounds are only tight
                      // at the beginning of each kmeans iteration (this is still tighter than the "loose bounds" case though)
  
                      // KYLE_HERE
                      // update centers

                      // if using angles, cluster center is average of normalized snapshots
                      double normalize = 1.0;
                      if (!this->euclideanDistances) {
                        normalize = (*snapshots)[iSnap].norm();
                        if (normalize > 0) normalize = 1.0/normalize;
                      }

                      DistSVec<double, dim> originalOldCenter(this->domain.getNodeDistInfo());
                      DistSVec<double, dim> originalNewCenter(this->domain.getNodeDistInfo());
                      originalOldCenter = (*(this->clusterCenters))[oldIndex];
                      originalNewCenter = (*(this->clusterCenters))[newIndex];
                      (*(this->clusterCenters))[oldIndex] *= (double((this->snapsInCluster)[oldIndex])/(double((this->snapsInCluster)[oldIndex])-1.0));
                      (*(this->clusterCenters))[oldIndex] -= (*snapshots)[iSnap]*normalize*(1.0/(double((this->snapsInCluster)[oldIndex])-1.0));
                      (*(this->clusterCenters))[newIndex] *= (double((this->snapsInCluster)[newIndex])/(double((this->snapsInCluster)[newIndex])+1.0));
                      (*(this->clusterCenters))[newIndex] += (*snapshots)[iSnap]*normalize*(1.0/(double((this->snapsInCluster)[newIndex])+1.0));
                      --((this->snapsInCluster)[oldIndex]);
                      ++((this->snapsInCluster)[newIndex]);
                      
                      // update uncertainties
                      DistSVec<double, dim> deltaOld(this->domain.getNodeDistInfo());
                      DistSVec<double, dim> deltaNew(this->domain.getNodeDistInfo());
                      deltaOld = originalOldCenter - (*(this->clusterCenters))[oldIndex];
                      deltaNew = originalNewCenter - (*(this->clusterCenters))[newIndex];
                      double normDeltaOld = deltaOld.norm();
                      double normDeltaNew = deltaNew.norm();
                      for (int jSnap=0;jSnap<nTotSnaps;++jSnap) {
                          uncertainty[jSnap][newIndex] += normDeltaNew;
                          uncertainty[jSnap][oldIndex] += normDeltaOld;  
                      }
                      nDistComputed = nDistComputed+2;
                    }
                  }
                }
              }
            }
          }
        }
      } // end distance comparison loop
      if (iter<iterMaxAggressive) {
        // update cluster centers now if doing aggressive iterations (otherwise this was done inside the distance comparison loop)
        int nSwitches = 0;
        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          if ((this->clusterIndex)[iSnap]!=clusterIndexOld[iSnap]) ++nSwitches;
        }
        // if (4*nSwitches > nTotSnaps + nClusters) it's cheaper to just recalculate
        if (4*nSwitches > nTotSnaps + this->nClusters) {
          for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
            (this->snapsInCluster)[iCluster]=0;
            (*(this->clusterCenters))[iCluster] = 0.0;
          }
          //KYLE_HERE
          DistSVec<double, dim>* tmpDistVec = new DistSVec<double, dim>(this->domain.getNodeDistInfo());
          for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
            *tmpDistVec = (*snapshots)[iSnap];
            // if using angles, cluster center is average of normalized snapshots
            if (!this->euclideanDistances) {
              double normalize = tmpDistVec->norm();
              if (normalize > 0) normalize = 1.0/normalize;
              *tmpDistVec *= normalize;
            }
            (*(this->clusterCenters))[(this->clusterIndex)[iSnap]] += *tmpDistVec;
            ++((this->snapsInCluster)[(this->clusterIndex)[iSnap]]);
          }
          delete tmpDistVec;
          tmpDistVec = NULL;

          for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
            double snapsInClusterInv = 0;
            if ((this->snapsInCluster)[iCluster] != 0) snapsInClusterInv = 1.0/double((this->snapsInCluster)[iCluster]);
            (*(this->clusterCenters))[iCluster] *= snapsInClusterInv;
          }
        } else { // otherwise it's cheaper to update incrementally
          for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
            int oldIndex = clusterIndexOld[iSnap];
            int newIndex = (this->clusterIndex)[iSnap];
            //KYLE_HERE
            double normalize = 1.0;
            if (!this->euclideanDistances) {
              normalize = (*snapshots)[iSnap].norm();
              if (normalize > 0) normalize = 1.0/normalize;
            }
            if (oldIndex != newIndex) {
                (this->clusterIndex)[iSnap] = newIndex;
                (*(this->clusterCenters))[oldIndex] *= (double((this->snapsInCluster)[oldIndex])/(double((this->snapsInCluster)[oldIndex])-1.0));
                (*(this->clusterCenters))[oldIndex] -= (*snapshots)[iSnap]*normalize*(1.0/(double((this->snapsInCluster)[oldIndex])-1.0));
                (*(this->clusterCenters))[newIndex] *= (double((this->snapsInCluster)[newIndex])/(double((this->snapsInCluster)[newIndex])+1.0));
                (*(this->clusterCenters))[newIndex] += (*snapshots)[iSnap]*normalize*(1.0/(double((this->snapsInCluster)[newIndex])+1.0));
                --((this->snapsInCluster)[oldIndex]);
                ++((this->snapsInCluster)[newIndex]);
            } 
          }
  
        }
      } 
  
      if (looseBounds && ((iter<iterMaxAggressive) || (iter==0))) {
        // increment loose bounds (unnecessary if iter>=iterMaxAggressive)
        std::vector<double> deltaCenters;
        deltaCenters.resize(this->nClusters,0.0);
        for (int iCluster=0;iCluster<this->nClusters;++iCluster) {
          deltaCenters[iCluster] = this->distanceFull((*(this->clusterCenters))[iCluster],(*clusterCentersOld)[iCluster]);
        }
        for (int iSnap=0;iSnap<nTotSnaps;++iSnap) {
          for (int iCluster=0;iCluster<this->nClusters;++iCluster) {
            uncertainty[iSnap][iCluster] = uncertainty[iSnap][iCluster] + deltaCenters[iCluster];
          }
        }
        nDistComputed = nDistComputed+this->nClusters;
      }
  
      clusterIndexOld.clear();
  
      if (tightBounds) {
        // recalculate tight bounds (regardless of whether iter<iterMaxAggressive)
        // calculate uncertainties (this gets a bit complicated for the Tight bounds case to avoid unnecessary work)
  
        // determine which distances need to be calculated for triangle inequalties
        std::vector<std::set<int> > deltaCentersToCalc; // [iCluster][(iteration set)] (need distance of these centers from current)
        deltaCentersToCalc.resize(this->nClusters);
        for (int iSnap=0;iSnap<nTotSnaps;++iSnap) {
          for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
            deltaCentersToCalc[iCluster].insert(distToCentersIteration[iSnap][iCluster]);
          }
        }
        // now that we've determined the set of distances that need to be calculated for the uncertainties, calculate them 
        std::vector<double> deltaCenters; // use vector instead of map to make random access faster (shouldn't be large)
        deltaCenters.resize((this->nClusters)*(iter+1),0.0);
        for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
          for (std::set<int>::iterator it=deltaCentersToCalc[iCluster].begin(); it!=deltaCentersToCalc[iCluster].end(); ++it) {
            int index = (this->nClusters)*(*it)+iCluster;
            double dist = this->distanceFull((*(this->clusterCenters))[iCluster],(*(clusterCentersLog[*it]))[iCluster]);
            deltaCenters[index] = dist;
            ++nDistComputed;
          }
        }
        // finally, store these distances in the uncertainty matrix (there should be many repeated values)
        for (int iSnap=0;iSnap<nTotSnaps;++iSnap) {
          for (int iCluster=0;iCluster<this->nClusters;++iCluster) {
            int index = (this->nClusters)*distToCentersIteration[iSnap][iCluster]+iCluster;
            uncertainty[iSnap][iCluster] = deltaCenters[index];
          }
        }
      }
  
      for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
        this->com->fprintf(stdout, " ... cluster %d has %d snaps \n", iCluster, (this->snapsInCluster)[iCluster]);
      }
      this->com->fprintf(stdout, " ... number of distances computed so far = %d, number skipped = %d\n", nDistComputed, nDistSkipped);
  
      residual = calcResidual(*(this->clusterCenters), *clusterCentersOld);
  
      this->com->fprintf(stdout, " ... absolute residual = %e (tolerance is set to %e)\n", residual, kMeansTol);
      ++iter;
    }  // end of kmeans clustering loop

    // delete any old cluster centers that were stored
    if (tightBounds) {
      for (int i=0; i<iter; ++i) {
         if (clusterCentersLog[i]) delete clusterCentersLog[i];
      }
      if (clusterCentersLog) delete [] clusterCentersLog;
      clusterCentersLog = NULL;

      for (int iSnap=0;iSnap<nTotSnaps;++iSnap)
        distToCentersIteration[iSnap].clear();
      distToCentersIteration.clear();

      clusterCentersOld = NULL;
    } else if (looseBounds) {
      if (clusterCentersOld) delete clusterCentersOld;
      clusterCentersOld = NULL;
    }

    for (int iSnap=0;iSnap<nTotSnaps;++iSnap) {
      distToCenters[iSnap].clear();
      uncertainty[iSnap].clear();
    }
    distToCenters.clear();
    uncertainty.clear();


    if (robConstruction->clustering.clusterFilesSeparately)  {
      // store results
      snapsInClusterAll[iFile].resize(this->nClusters);
      clusterCentersAll[iFile] = new VecSet< DistSVec<double, dim> >((this->nClusters), this->domain.getNodeDistInfo());
      for (int iCluster=0;iCluster<this->nClusters;++iCluster) {
        snapsInClusterAll[iFile][iCluster] = (this->snapsInCluster)[iCluster];
        (*clusterCentersAll[iFile])[iCluster] = (*this->clusterCenters)[iCluster];
      }
      clusterIndexAll[iFile].resize(nTotSnaps);
      for (int iSnap=0;iSnap<nTotSnaps;++iSnap) {
        clusterIndexAll[iFile][iSnap] = (this->clusterIndex)[iSnap];
      }
    
      // delete
      delete [] (this->clusterIndex);
      this->clusterIndex = NULL;
      delete [] (this->snapsInCluster);
      this->snapsInCluster = NULL;
      delete (this->clusterCenters);
      this->clusterCenters = NULL;
      delete snapshots;
      snapshots = NULL;
    }   
 
  }

  if (robConstruction->clustering.clusterFilesSeparately)  {
    this->com->fprintf(stdout, "\nCombining clusters from each data set\n");

    // form the final data structures
    nTotSnaps = this->snap->numVectors();
    this->clusterIndex = new int[nTotSnaps];
    this->clusterCenters = new VecSet< DistSVec<double, dim> >((this->nClusters), this->domain.getNodeDistInfo());
    this->snapsInCluster = new int[(this->nClusters)];

    // initialize to the values found for the first file
    for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
      (this->snapsInCluster)[iCluster] = snapsInClusterAll[0][iCluster];
      (*this->clusterCenters)[iCluster] = (*(clusterCentersAll[0]))[iCluster];
    }
    for (int iSnap=0; iSnap<clusterIndexAll[0].size(); ++iSnap) {
      (this->clusterIndex)[iSnap] = clusterIndexAll[0][iSnap];
    }

    // associate clusters to ensure that each cluster has snapshots from every file
 
    std::vector<std::vector<double> > distMat;
    distMat.resize(this->nClusters);

    snapCount = this->stateSnapsFromFile[0];
    for (int iFile=1; iFile<nSnapshotFiles; ++iFile) {
      
      // calc all distances
      for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
        distMat[iCluster].resize(this->nClusters, 0.0);
        for (int jCluster=0; jCluster<this->nClusters; ++jCluster) {
          distMat[iCluster][jCluster] = this->distanceFull((*this->clusterCenters)[iCluster], (*clusterCentersAll[iFile])[jCluster]);
        }
      }

      // perform small combinatorial search to determine best way to combine clusters from 
      // next file with clusterIndex, clusterCenters, snapsInCluster

      int* best = new int[this->nClusters];
      int* test = new int[this->nClusters];
      for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
        test[iCluster] = iCluster;
      }

      std::sort(test,test+this->nClusters);
      double testObjective = 0;
      double bestObjective = 0;
      for (int iCluster=0; iCluster<this->nClusters; ++iCluster)
        bestObjective+=pow(distMat[iCluster][test[iCluster]],2);

      for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
          best[iCluster] = test[iCluster];
      }

      while ( std::next_permutation(test,test+this->nClusters) ) {
        testObjective = 0;
        for (int iCluster=0; iCluster<this->nClusters; ++iCluster)
          testObjective+=pow(distMat[iCluster][test[iCluster]],2);

        if (testObjective<bestObjective) {
          for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
            best[iCluster] = test[iCluster];
          }
          bestObjective = testObjective;
        }
      }

      for (int iSnap=0; iSnap<this->stateSnapsFromFile[iFile]; ++iSnap){
        (this->clusterIndex)[snapCount+iSnap] = best[clusterIndexAll[iFile][iSnap]];
      }
      snapCount += this->stateSnapsFromFile[iFile];
 
      for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
        int snapsInNewCluster = (this->snapsInCluster)[iCluster] + snapsInClusterAll[iFile][best[iCluster]];
        (*this->clusterCenters)[iCluster] = (((double) (this->snapsInCluster)[iCluster])/((double)snapsInNewCluster)) * (*this->clusterCenters)[iCluster]
                                            + (((double)snapsInClusterAll[iFile][best[iCluster]])/((double)snapsInNewCluster)) * (*clusterCentersAll[iFile])[best[iCluster]];
        (this->snapsInCluster)[iCluster] = snapsInNewCluster;
      }  

      delete [] best;
      delete [] test;
    }
 
    // clean up
    for (int iFile=0;iFile<nSnapshotFiles;++iFile) {
      clusterIndexAll[iFile].clear();
      snapsInClusterAll[iFile].clear();
      delete clusterCentersAll[iFile];
    }

    clusterIndexAll.clear();
    snapsInClusterAll.clear();
    delete [] clusterCentersAll;
    clusterCentersAll = NULL;

    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
        this->com->fprintf(stdout, " ... cluster %d has %d snaps \n", iCluster, (this->snapsInCluster)[iCluster]);
    }

  }

  // after clustering assimilate small clusters
  int emptyClusterCount = 0;
  if (minClusterSize>nTotSnaps) minClusterSize = nTotSnaps;

  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
    if ((this->snapsInCluster)[iCluster]==0) {
      // if no snapshots in this cluster, skip
      ++emptyClusterCount;
    } else {
    // if smaller than tolerance, add to nearest cluster
      if ((this->snapsInCluster)[iCluster]<minClusterSize) {
        this->com->fprintf(stderr, "*** Warning: combining small cluster with nearest neighbor\n");
        int index1 = 0; // closest center (should be itself)
        int index2 = 0; // second closest center (the one we want)

        this->closestCenterFull((*(this->clusterCenters))[iCluster], &index1, &index2);
        (*(this->clusterCenters))[index2] =   (*(this->clusterCenters))[index2]*(double((this->snapsInCluster)[index2])/
                                                double((this->snapsInCluster)[index2]+(this->snapsInCluster)[iCluster]))
                                            + (*(this->clusterCenters))[iCluster]*(double((this->snapsInCluster)[iCluster])/
                                                double((this->snapsInCluster)[index2]+(this->snapsInCluster)[iCluster]));
 
        (*(this->clusterCenters))[iCluster] = 0.0;
 
        (this->snapsInCluster)[index2] = (this->snapsInCluster)[index2] + (this->snapsInCluster)[iCluster];
        (this->snapsInCluster)[iCluster] = 0;  

        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          if ((this->clusterIndex)[iSnap]==iCluster)  (this->clusterIndex)[iSnap]=index2;
        }
        
        ++emptyClusterCount;
      }
    }
  }

  if (this->nClusters == emptyClusterCount) {
    this->com->fprintf(stderr, "*** Error: no clusters remain\n");
    exit(-1);
  }

  // remove any empty clusters, renumber existing clusters
  if (emptyClusterCount>0) {

    this->com->fprintf(stderr, "*** Warning: Deleting %d empty clusters\n", emptyClusterCount);

    VecSet< DistSVec<double, dim> >* clusterCentersNew =  new VecSet< DistSVec<double, dim> >(((this->nClusters)-emptyClusterCount), this->domain.getNodeDistInfo()); 
    int* snapsInClusterNew = new int[((this->nClusters)-emptyClusterCount)];

    int renumberedIndices[(this->nClusters)];
    int clusterCount = 0;

    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
      if ((this->snapsInCluster)[iCluster]==0) {
        renumberedIndices[iCluster] = -1;
      } else {
        renumberedIndices[iCluster] = clusterCount; 
        (*clusterCentersNew)[clusterCount]=(*(this->clusterCenters))[iCluster];
        snapsInClusterNew[clusterCount]=(this->snapsInCluster)[iCluster];
        ++clusterCount;
      }
    }

    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
      (this->clusterIndex)[iSnap] = renumberedIndices[(this->clusterIndex)[iSnap]];
    }

    delete (this->clusterCenters);
    (this->clusterCenters) = clusterCentersNew;
    clusterCentersNew = NULL;

    delete [] (this->snapsInCluster);
    (this->snapsInCluster) = snapsInClusterNew; 
    snapsInClusterNew = NULL;  

    (this->nClusters) = (this->nClusters) - emptyClusterCount;
  }
  

  // find snapshot closest to each center
  (this->nearestSnapsToCenters) = new VecSet< DistSVec<double, dim> >((this->nClusters), this->domain.getNodeDistInfo());
  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
      
    int sortSize = (this->snapsInCluster)[iCluster]; 
    sortStruct* snapDist = new sortStruct[sortSize]; // this struct was defined earlier in this file 

    snapCount = 0;
    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
      if ((this->clusterIndex)[iSnap]==iCluster) { // only need to sort snapshots that are inside the current cluster
        snapDist[snapCount].snapIndex = iSnap;
        snapDist[snapCount].dist = this->distanceFull((*(this->snap))[iSnap],(*(this->clusterCenters))[iCluster]);
        ++snapCount;
      }
    }

    sort(snapDist, snapDist+sortSize);

    (*(this->nearestSnapsToCenters))[iCluster] = (*(this->snap))[snapDist[0].snapIndex];

    delete [] snapDist;
    snapDist = NULL;
  }

 
  // add overlap if required
  if ((percentOverlap>0) && (this->nClusters>1)) {

    this->com->fprintf(stdout, "Adding additional vectors to clusters (PercentOverlap = %2.1f%%)\n", percentOverlap);

    int index1 = 0;
    int index2 = 0;
    double dist1 = 0;
    double dist2 = 0;

    //determine which clusters are neighbors
    //approach: if a cluster is second closest to a snapshot in another cluster, then those two clusters are neighbors
    (this->clusterNeighbors) = new int*[(this->nClusters)];
    (this->clusterNeighborsCount) = new int[(this->nClusters)];

    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
      (this->clusterNeighborsCount)[iCluster] = 0;
      (this->clusterNeighbors)[iCluster] = new int[1];
    }

    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        if ((this->clusterIndex)[iSnap]==iCluster) { // only consider snapshots that are inside of the current cluster
          this->closestCenterFull( (*(this->snap))[iSnap], &index1, &index2, &dist1, &dist2);
          //add the second closest cluster as a neighbor of current cluster (if this is the first time we've found it)
          bool unique = true;
          for (int iNeighbor=0; iNeighbor<(this->clusterNeighborsCount)[iCluster]; ++iNeighbor) {
            if (index2 == (this->clusterNeighbors)[iCluster][iNeighbor]) unique = false;
          }
          if (unique) {
            int* clusterNeighborsNew = new int[((this->clusterNeighborsCount)[iCluster])+1];
            for (int iNeighbor=0; iNeighbor<(this->clusterNeighborsCount)[iCluster]; ++iNeighbor) {
              clusterNeighborsNew[iNeighbor] = (this->clusterNeighbors)[iCluster][iNeighbor];
            }
 
           clusterNeighborsNew[((this->clusterNeighborsCount)[iCluster])] = index2;
            delete [] (this->clusterNeighbors)[iCluster];
            (this->clusterNeighbors)[iCluster] = clusterNeighborsNew;
            clusterNeighborsNew = NULL;
            ++(this->clusterNeighborsCount)[iCluster];
         }
          //also add the current cluster as a neighbor of second closest cluster
          unique = true;
          for (int iNeighbor=0; iNeighbor<(this->clusterNeighborsCount)[index2]; ++iNeighbor) {
            if (iCluster == (this->clusterNeighbors)[index2][iNeighbor]) unique = false;
          }
          if (unique) {
            int* clusterNeighborsNew = new int[((this->clusterNeighborsCount)[index2])+1];
            for (int iNeighbor=0; iNeighbor<(this->clusterNeighborsCount)[index2]; ++iNeighbor) {
              clusterNeighborsNew[iNeighbor] = (this->clusterNeighbors)[index2][iNeighbor];
            }
 
            clusterNeighborsNew[((this->clusterNeighborsCount)[index2])] = iCluster;
            delete [] (this->clusterNeighbors)[index2];
            (this->clusterNeighbors)[index2] = clusterNeighborsNew;
            clusterNeighborsNew = NULL;
            ++(this->clusterNeighborsCount)[index2];
         }
        }
      }
    } 

    // if clustering increments, need to actually cluster increments by the preceding increment 
    // Fairly simply, just need to shift all cluster associations.
    // For the time being , I won't shift the cluster association of any IC... 
    // (This is still consistent, provided that all ROBs are tested for the IC of the online simulation)
    if (this->ioData->romDatabase.avgIncrementalStates==NonlinearRomFileSystemData::AVG_INCREMENTAL_STATES_TRUE) {

      int nSnapsHandled = 0;
      std::vector<int> newSnapsInCluster(this->nClusters,0);
      std::vector<int> newClusterIndex(nTotSnaps,0);

      for (int iFile=0; iFile<nSnapshotFiles; ++iFile) {
        newClusterIndex[nSnapsHandled]=(this->clusterIndex)[nSnapsHandled];
        ++newSnapsInCluster[newClusterIndex[nSnapsHandled]];
        ++nSnapsHandled;
        for (int iSnap=1; iSnap<this->stateSnapsFromFile[iFile]; ++iSnap) {
          newClusterIndex[nSnapsHandled] = (this->clusterIndex)[nSnapsHandled-1];
          ++newSnapsInCluster[newClusterIndex[nSnapsHandled]];
          ++nSnapsHandled;
        }
      }
      

      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) (this->clusterIndex)[iSnap] = newClusterIndex[iSnap];                         

      emptyClusterCount = 0;
      for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
        (this->snapsInCluster)[iCluster] = newSnapsInCluster[iCluster];
        if ((this->snapsInCluster)[iCluster] == 0) ++emptyClusterCount;
      }

      if (emptyClusterCount>0) {
        this->com->fprintf(stderr, "*** Warning: Deleting %d empty clusters\n", emptyClusterCount);

        VecSet< DistSVec<double, dim> >* clusterCentersNew =  new VecSet< DistSVec<double, dim> >(((this->nClusters)-emptyClusterCount), this->domain.getNodeDistInfo());
        int* snapsInClusterNew = new int[((this->nClusters)-emptyClusterCount)];
    
        int renumberedIndices[(this->nClusters)];
        int clusterCount = 0;
    
        for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
          if ((this->snapsInCluster)[iCluster]==0) {
            renumberedIndices[iCluster] = -1;
          } else {
            renumberedIndices[iCluster] = clusterCount;
            (*clusterCentersNew)[clusterCount]=(*(this->clusterCenters))[iCluster];
            snapsInClusterNew[clusterCount]=(this->snapsInCluster)[iCluster];
            ++clusterCount;
          }
        }
    
        for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
          (this->clusterIndex)[iSnap] = renumberedIndices[(this->clusterIndex)[iSnap]];
        }
    
        delete (this->clusterCenters);
        (this->clusterCenters) = clusterCentersNew;
        clusterCentersNew = NULL;
    
        delete [] (this->snapsInCluster);
        (this->snapsInCluster) = snapsInClusterNew;
        snapsInClusterNew = NULL;
    
        (this->nClusters) = (this->nClusters) - emptyClusterCount;
      }

    }


    index1 = 0;
    index2 = 0;
    dist1 = 0;
    dist2 = 0;

    int origSnapsInCluster[(this->nClusters)];
    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) origSnapsInCluster[iCluster] = (this->snapsInCluster)[iCluster];

    //share shapshots between neighboring clusters
    (this->clusterSnapshotMap) = new int*[(this->nClusters)];
    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {  

      int numToAdd = int(ceil(double((this->snapsInCluster)[iCluster])*percentOverlap*.01/double((this->clusterNeighborsCount)[iCluster])));
      if (((this->snapsInCluster)[iCluster]+(numToAdd*(this->clusterNeighborsCount)[iCluster]))>nTotSnaps) 
        numToAdd=int(floor(double((nTotSnaps-(this->snapsInCluster)[iCluster]))/double((this->clusterNeighborsCount)[iCluster])));
      (this->clusterSnapshotMap)[iCluster] = new int[((this->snapsInCluster)[iCluster]+(numToAdd*(this->clusterNeighborsCount)[iCluster]))];

      // first add all snapshots that were originally in the cluster
      int mappedCount = 0;
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        if ((this->clusterIndex)[iSnap]==iCluster) {
          (this->clusterSnapshotMap)[iCluster][mappedCount] = iSnap;
          ++mappedCount;
        }
      }

      for (int iNeighbor=0; iNeighbor<(this->clusterNeighborsCount)[iCluster]; ++iNeighbor) {

        if (iCluster == (this->clusterNeighbors)[iCluster][iNeighbor]) {
          this->com->fprintf(stderr, "*** Warning: Cluster  %d is its own neighbor (likely because KMeans didin't full converge). Ignoring...\n", iCluster);
        } else {

          int sortSize = origSnapsInCluster[(this->clusterNeighbors)[iCluster][iNeighbor]]; 
          sortStruct* snapDist = new sortStruct[sortSize]; // this struct was defined earlier
  
          snapCount = 0;
              //KYLE_HERE
          if (this->euclideanDistances) {
            DistSVec<double, dim>* tmpDistVec = new DistSVec<double, dim>(this->domain.getNodeDistInfo());
            *tmpDistVec = (*(this->clusterCenters))[iCluster] - (*(this->clusterCenters))[(this->clusterNeighbors)[iCluster][iNeighbor]];
            double offset = (pow(((*(this->clusterCenters))[iCluster]).norm(),2) - pow(((*(this->clusterCenters))[(this->clusterNeighbors)[iCluster][iNeighbor]]).norm(),2))*(-0.5);
    
            for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
              if ((this->clusterIndex)[iSnap]==(this->clusterNeighbors)[iCluster][iNeighbor]) { // don't need to sort snapshots belonging to the current cluster
                snapDist[snapCount].snapIndex = iSnap;
                snapDist[snapCount].dist = abs( (*tmpDistVec) * (*(this->snap))[iSnap] + offset );  //don't need the abs? KMW
                ++snapCount;
              }
            }
            delete tmpDistVec;
            tmpDistVec = NULL;
          } else { // angles
  
            DistSVec<double, dim> projVec1(this->domain.getNodeDistInfo());
            DistSVec<double, dim> projVec2(this->domain.getNodeDistInfo());
  
            projVec1 = (*(this->clusterCenters))[iCluster];
            projVec1 *= 1.0/projVec1.norm(); // norm can't be zero
            projVec2 = (*(this->clusterCenters))[(this->clusterNeighbors)[iCluster][iNeighbor]];
            projVec2 = projVec2 - (projVec1*projVec2)*projVec1;
            projVec2 *= 1.0/projVec2.norm(); // again, norm can't be zero

            DistSVec<double, dim>* tmpDistVec = new DistSVec<double, dim>(this->domain.getNodeDistInfo());
  
            for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
              if ((this->clusterIndex)[iSnap]==(this->clusterNeighbors)[iCluster][iNeighbor]) { // don't need to sort snapshots belonging to the current cluster
                snapDist[snapCount].snapIndex = iSnap;
                *tmpDistVec = projVec1 * (projVec1*((*(this->snap))[iSnap]));
                *tmpDistVec += projVec2 * (projVec2*((*(this->snap))[iSnap]));
                snapDist[snapCount].dist = this->distanceFull(*tmpDistVec, (*(this->clusterCenters))[iCluster]);
                ++snapCount;
              }
            }
            delete tmpDistVec;
            tmpDistVec = NULL;
                
          }
          sort(snapDist, snapDist+sortSize);
  
          for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
            for (int jSnap=0; jSnap<numToAdd; ++jSnap) {
              if (snapDist[jSnap].snapIndex==iSnap) {
                (this->clusterSnapshotMap)[iCluster][mappedCount] = iSnap;
                ++mappedCount;
                ++((this->snapsInCluster)[iCluster]);
              }
            }
          }
  
          delete [] snapDist;
          snapDist = NULL;
        }
      }

      this->com->fprintf(stdout, " ... cluster %d has %d snaps \n", iCluster, (this->snapsInCluster)[iCluster]);

    }
  } else { // no overlap
    (this->clusterSnapshotMap) = new int*[(this->nClusters)];
    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {  
      (this->clusterSnapshotMap)[iCluster] = new int[(this->snapsInCluster)[iCluster]];

      int mappedCount = 0;
      for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
        if ((this->clusterIndex)[iSnap]==iCluster) {
          (this->clusterSnapshotMap)[iCluster][mappedCount] = iSnap;
          ++mappedCount;
        }
      }
    }
  }

  // create data structure needed for clustering non-state FOM information (explained in header file)
  stateSnapshotClustersAfterOverlap.resize(nSnapshotFiles);
  for (int iFile=0; iFile<nSnapshotFiles; ++iFile) {
    stateSnapshotClustersAfterOverlap[iFile].resize(this->stateSnapsFromFile[iFile]);
    for (int iSnap=0; iSnap<this->stateSnapsFromFile[iFile]; ++iSnap) {
      stateSnapshotClustersAfterOverlap[iFile][iSnap].clear();
      stateSnapshotClustersAfterOverlap[iFile][iSnap].resize(this->nClusters,false);
    }
  }

  int** snapInfo = new int*[nTotSnaps];  // temporary - for constructing stateSnapshotClustersAfterOverlap
  for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
    snapInfo[iSnap] = new int[2];
  }

  snapCount = 0;
  for (int iFile=0; iFile<nSnapshotFiles; ++iFile) {
    for (int iSnap=0; iSnap<this->stateSnapsFromFile[iFile]; ++iSnap) {
      snapInfo[snapCount][0] = iFile;// file
      snapInfo[snapCount][1] = iSnap;// snapshot number from file
      ++snapCount;
    }
  }

  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
    for (int iSnap=0; iSnap<(this->snapsInCluster[iCluster]); ++iSnap) {
      int currentSnap = this->clusterSnapshotMap[iCluster][iSnap];
      stateSnapshotClustersAfterOverlap[snapInfo[currentSnap][0]][snapInfo[currentSnap][1]][iCluster] = true;
    }
  }

  for (int iSnap=0; iSnap<nTotSnaps; ++iSnap)
    delete [] snapInfo[iSnap];

  delete [] snapInfo;

}




//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::localPod(const char* basisType) {

  // for limited memory SVD
  int maxVecStorage;

  //these parameters control the size of the local bases
  double singValTolerance;
  double maxEnergyRetained;
  int maxBasisSize;
  int minBasisSize;
  int podMethod;
  int randMatDimension;
  int nPowerIts;
  int compareSVDMethods;
  int initialCluster;
  bool testProbSVD;

  if (strcmp(basisType, "state")==0) {
    if (strcmp(this->stateBasisName,"")==0) return;
    if (robConstruction->state.snapshots.subtractNearestSnapsToCenters) this->readNearestSnapsToCenters();
    maxVecStorage = robConstruction->state.dataCompression.maxVecStorage; 
    singValTolerance = robConstruction->state.dataCompression.singValTolerance;
    maxEnergyRetained = robConstruction->state.dataCompression.maxEnergyRetained;
    maxBasisSize = robConstruction->state.dataCompression.maxBasisSize;
    minBasisSize = robConstruction->state.dataCompression.minBasisSize;
    podMethod = robConstruction->state.dataCompression.podMethod;
    randMatDimension = robConstruction->state.dataCompression.randMatDimension;
    nPowerIts = robConstruction->state.dataCompression.nPowerIts;
    compareSVDMethods = robConstruction->state.dataCompression.compareSVDMethods;
    initialCluster = robConstruction->state.dataCompression.initialCluster;
    testProbSVD = (robConstruction->state.dataCompression.testProbabilisticSVD==1) ? true : false;
  } else if (strcmp(basisType,"residual")==0) {
    if (strcmp(this->residualBasisName,"")==0) return;
    maxVecStorage = robConstruction->residual.dataCompression.maxVecStorage;
    singValTolerance = robConstruction->residual.dataCompression.singValTolerance;
    maxEnergyRetained = robConstruction->residual.dataCompression.maxEnergyRetained;
    maxBasisSize = robConstruction->residual.dataCompression.maxBasisSize;
    minBasisSize = robConstruction->residual.dataCompression.minBasisSize;
    podMethod = robConstruction->residual.dataCompression.podMethod;
    randMatDimension = robConstruction->residual.dataCompression.randMatDimension;
    nPowerIts = robConstruction->residual.dataCompression.nPowerIts;
    compareSVDMethods = robConstruction->residual.dataCompression.compareSVDMethods;
    initialCluster = robConstruction->residual.dataCompression.initialCluster;
    testProbSVD = (robConstruction->residual.dataCompression.testProbabilisticSVD==1) ? true : false;
  } else if (strcmp(basisType,"jacAction")==0) {
    if (strcmp(this->jacActionBasisName,"")==0) return;
    maxVecStorage = robConstruction->jacAction.dataCompression.maxVecStorage;
    singValTolerance = robConstruction->jacAction.dataCompression.singValTolerance;
    maxEnergyRetained = robConstruction->jacAction.dataCompression.maxEnergyRetained;
    maxBasisSize = robConstruction->jacAction.dataCompression.maxBasisSize;
    minBasisSize = robConstruction->jacAction.dataCompression.minBasisSize;
    podMethod = robConstruction->jacAction.dataCompression.podMethod;
    randMatDimension = robConstruction->jacAction.dataCompression.randMatDimension;
    nPowerIts = robConstruction->jacAction.dataCompression.nPowerIts;
    compareSVDMethods = robConstruction->jacAction.dataCompression.compareSVDMethods;
    initialCluster = robConstruction->jacAction.dataCompression.initialCluster;
    testProbSVD = (robConstruction->jacAction.dataCompression.testProbabilisticSVD==1) ? true : false;
  } else if (strcmp(basisType,"krylov")==0) {
    if (strcmp(this->krylovBasisName,"")==0) return;
    maxVecStorage = robConstruction->krylov.dataCompression.maxVecStorage;
    singValTolerance = robConstruction->krylov.dataCompression.singValTolerance;
    maxEnergyRetained = robConstruction->krylov.dataCompression.maxEnergyRetained;
    maxBasisSize = robConstruction->krylov.dataCompression.maxBasisSize;
    minBasisSize = robConstruction->krylov.dataCompression.minBasisSize;
    randMatDimension = robConstruction->krylov.dataCompression.randMatDimension;
    nPowerIts = robConstruction->krylov.dataCompression.nPowerIts;
    compareSVDMethods = robConstruction->krylov.dataCompression.compareSVDMethods;
    initialCluster = robConstruction->krylov.dataCompression.initialCluster;
    testProbSVD = (robConstruction->krylov.dataCompression.testProbabilisticSVD==1) ? true : false;
  } else if (strcmp(basisType,"sensitivity")==0) {
    if (strcmp(this->sensitivityBasisName,"")==0) return;
    maxVecStorage = robConstruction->sensitivity.dataCompression.maxVecStorage;
    singValTolerance = robConstruction->sensitivity.dataCompression.singValTolerance;
    maxEnergyRetained = robConstruction->sensitivity.dataCompression.maxEnergyRetained;
    maxBasisSize = robConstruction->sensitivity.dataCompression.maxBasisSize;
    minBasisSize = robConstruction->sensitivity.dataCompression.minBasisSize;
    podMethod = robConstruction->sensitivity.dataCompression.podMethod;
    randMatDimension = robConstruction->sensitivity.dataCompression.randMatDimension;
    nPowerIts = robConstruction->sensitivity.dataCompression.nPowerIts;
    compareSVDMethods = robConstruction->sensitivity.dataCompression.compareSVDMethods;
    initialCluster = robConstruction->sensitivity.dataCompression.initialCluster;
    testProbSVD = (robConstruction->sensitivity.dataCompression.testProbabilisticSVD==1) ? true : false;
  } else {
    this->com->fprintf(stderr, "*** Error: unexpected snapshot type %s\n", basisType);
    exit (-1);
  }

  if (podMethod==DataCompressionData::PROBABILISTIC_SVD && (randMatDimension>0)) 
    maxBasisSize=min(randMatDimension, maxBasisSize);

  // read cluster centers
  this->readClusterCenters("centers");

  bool limitedMemorySVD = (maxVecStorage <= 0 ) ? false : true;

  for (int iCluster=initialCluster; iCluster<(this->nClusters); ++iCluster) {

    if (strcmp(basisType,"sensitivity")==0)
      iCluster = -2;

    int nTotSnaps = 0;
    VecSet< DistSVec<double, dim> >* Utrue = NULL;
    std::vector<double> singVals;
    FullM* Vtrue = NULL;

    // read snapshots and preprocess
    if ((strcmp(basisType, "state")==0) && (this->snap!=NULL) && (this->nClusters)==1) {
      this->readClusteredSnapshots(iCluster, true, basisType, 0, (maxVecStorage-1), true);
    } else {
      this->readClusteredSnapshots(iCluster, true, basisType, 0, (maxVecStorage-1));  
    }
    if (limitedMemorySVD && (this->snap->numVectors() == maxVecStorage)) {  // limited memory SVD algorithm

      this->com->fprintf(stdout, " ... beginning limited memory SVD algorithm \n");

      if (maxBasisSize <= 0) {
          this->com->fprintf(stderr, "*** Error: the limited memory SVD algorithm requires maxBasisSize to be specified \n");
          exit(-1);
      }

      if (maxBasisSize >= maxVecStorage) {
          this->com->fprintf(stderr, "*** Error: maxBasisSize >= maxVecStorage.  This doesn't make any sense.  The intended scenario is:  maxBasisSize < maxVecStorage < nTotSnaps.\n");
          exit(-1);
      }

      VecSet< DistSVec<double, dim> >* fullSnaps = NULL; // matrix for final SVD

      int nStoredSnaps = 0;

      while (this->snap->numVectors() > maxBasisSize) {
        int nSnaps = this->snap->numVectors();

        VecSet< DistSVec<double, dim> > UtrueTmp(nSnaps, this->domain.getNodeDistInfo());
        std::vector<double> singValsTmp(nSnaps,0.0);
        FullM VtrueTmp(nSnaps,0.0);

        int nKeep = maxBasisSize;
        this->com->fprintf(stdout, " ... performing SVD on matrix of size %d, storing first %d vectors to final SVD matrix \n", nSnaps, nKeep);

        int randMatDimensionTmp = (randMatDimension>0) ? randMatDimension : nSnaps; 

        if (compareSVDMethods) testProbabilisticSVD(this->snap, UtrueTmp, singValsTmp, VtrueTmp, podMethod, randMatDimensionTmp, nPowerIts, true);

        SVD(this->snap, UtrueTmp, singValsTmp, VtrueTmp, podMethod, randMatDimensionTmp, nPowerIts, true);
       
        if (this->snap) delete (this->snap);
        this->snap = NULL;         
 
        VecSet< DistSVec<double, dim> >* fullSnapsNew = new VecSet< DistSVec<double, dim> >(nStoredSnaps + nKeep, this->domain.getNodeDistInfo());
        for (int iVec=0;iVec<nStoredSnaps;++iVec) {
          (*fullSnapsNew)[iVec] = (*fullSnaps)[iVec];
        }
        for (int iVec=nStoredSnaps;iVec<(nStoredSnaps+nKeep);++iVec) {
          (*fullSnapsNew)[iVec] = UtrueTmp[iVec-nStoredSnaps] * singValsTmp[iVec-nStoredSnaps];
        }
        nStoredSnaps += nKeep;
        if (fullSnaps) delete fullSnaps;
        fullSnaps = fullSnapsNew;
        fullSnapsNew = NULL;

        nTotSnaps += nSnaps;

        // read next chunk of snapshots
        this->readClusteredSnapshots(iCluster, true, basisType, nTotSnaps, (nTotSnaps + maxVecStorage - 1));  
      }

      // append extra (<maxBasisSize) snapshots to the fullSnaps matrix
      int nSnaps = this->snap->numVectors();
      if (nSnaps > 0) {
        this->com->fprintf(stdout, " ... appending trailing %d vectors to final SVD matrix \n", nSnaps);
        VecSet< DistSVec<double, dim> >* fullSnapsNew = new VecSet< DistSVec<double, dim> >(nStoredSnaps + nSnaps, this->domain.getNodeDistInfo());
        for (int iVec=0;iVec<nStoredSnaps;++iVec) {
          (*fullSnapsNew)[iVec] = (*fullSnaps)[iVec];
        }
        for (int iVec=nStoredSnaps;iVec<(nStoredSnaps+nSnaps);++iVec) {
          (*fullSnapsNew)[iVec] = ((*this->snap))[iVec-nStoredSnaps];
        }
        nStoredSnaps += nSnaps;
        if (fullSnaps) delete fullSnaps;
        fullSnaps = fullSnapsNew;
        fullSnapsNew = NULL;
      }

      // SVD to compute final U
      nTotSnaps = fullSnaps->numVectors();
      this->com->fprintf(stdout, " ... performing final SVD on matrix of size %d\n", nTotSnaps);
      Utrue = new VecSet< DistSVec<double, dim> >(nTotSnaps, this->domain.getNodeDistInfo());
      Vtrue = new FullM(nTotSnaps);

      if (compareSVDMethods) testProbabilisticSVD(this->snap, *Utrue, singVals, *Vtrue, podMethod, randMatDimension, nPowerIts, true);

      SVD(fullSnaps, *Utrue, singVals, *Vtrue, podMethod, randMatDimension, nPowerIts, true, testProbSVD);

      if (fullSnaps) delete fullSnaps;
      fullSnaps = NULL;

    } else { // regular SVD
    
      nTotSnaps = this->snap->numVectors();
      Utrue = new VecSet< DistSVec<double, dim> >(nTotSnaps, this->domain.getNodeDistInfo());
      Vtrue = new FullM(nTotSnaps);

      int randMatDimensionTmp = (randMatDimension>0) ? randMatDimension : nTotSnaps; 

      if (compareSVDMethods) testProbabilisticSVD(this->snap, *Utrue, singVals, *Vtrue, podMethod, randMatDimensionTmp, nPowerIts, true);

      SVD(this->snap, *Utrue, singVals, *Vtrue, podMethod, randMatDimensionTmp, nPowerIts, true, testProbSVD);
      if (this->snap) delete (this->snap);
      (this->snap) = NULL;

      /*SVD(this->snap, *Utrue, singVals->data(), *Vtrue, nTotSnaps, podMethod, true);
  
      // check svd
      double errorNorm,maxErr,avgErr;
      DistSVec<double,dim> error( this->domain.getNodeDistInfo() );
      maxErr = 0.0;
      avgErr = 0.0;
      for (int iSnap = 0; iSnap < nTotSnaps; ++iSnap) {
        error = (*this->snap)[iSnap];
        for (int jSnap = 0; jSnap < nTotSnaps; ++jSnap)
          error = error - (((*singVals)[jSnap]*(*Vtrue)[iSnap][jSnap])*(*Utrue)[jSnap]);
        errorNorm = ((((*this->snap)[iSnap]).norm()) > 1e-15) ? error.norm()/(((*this->snap)[iSnap]).norm()) : 0.0;
        avgErr += errorNorm;
        if (errorNorm > maxErr)
          maxErr = errorNorm;   
      }
      avgErr /= nTotSnaps;
     
      this->com->fprintf(stderr, " ... Average error on Snapshots after SVD = %e\n", avgErr);  
      this->com->fprintf(stderr, " ... Maximum error on Snapshots after SVD = %e\n", maxErr);
      // end check svd
 
      delete (this->snap);
      (this->snap) = NULL; */

    }
  
    //check how many vectors to keep
    double singValSquaresTotal = 0;
    for(int i = 0; i < nTotSnaps; ++i){
      singValSquaresTotal += pow(singVals[i],2);
    }
  
    int podSize = 0; 
    double singValSquaresPartialSum = 0; 
    double target = maxEnergyRetained * singValSquaresTotal; 

    for(int iSnap=0; iSnap<nTotSnaps; ++iSnap){ 
      singValSquaresPartialSum += pow(singVals[iSnap],2);  
      podSize = iSnap+1; 
      if (podSize == maxBasisSize) { // setting maxBasisSize <= 0 guarantees that this is always false 
        this->com->barrier();
        this->com->fprintf(stdout, "Retaining the specified limit of %d vectors in basis %d \n", podSize, iCluster); 
        this->com->fprintf(stdout, "This basis retains %2.10f%% of the total energy\n", (singValSquaresPartialSum/singValSquaresTotal)*100); 
        break; 
      } else if ((singValSquaresPartialSum >= target) && (podSize>=minBasisSize)){ 
        this->com->barrier();
        this->com->fprintf(stdout, "Reached specified energy (%2.10f%%) \n", (singValSquaresPartialSum/singValSquaresTotal)*100); 
        this->com->fprintf(stdout, "Retaining %i vectors in basis %i \n", podSize, iCluster); 
        break; 
      } else if ((singVals[iSnap] <= singValTolerance) && (podSize>=minBasisSize)) { // singValTolerance is 1e-6 by default 
        this->com->barrier();
        this->com->fprintf(stderr, "*** Warning: Reached the singular value tolerance (%e)\n", singValTolerance); 
        this->com->fprintf(stderr, "Retaining %i vectors in basis %i \n", podSize, iCluster); 
        break; 
      } else if ((podMethod==DataCompressionData::PROBABILISTIC_SVD) 
                  && (randMatDimension>0) && (podSize>=randMatDimension)) {
        this->com->barrier();
        this->com->fprintf(stderr, "*** Warning: Keeping the full basis generated by the probabilistic SVD algorithm (Random_Matrix_Dim=%d)\n", randMatDimension);
        this->com->fprintf(stderr, "Retaining %i vectors in basis %i \n", podSize, iCluster);
        break;
      }
    } 

    if (this->basis) delete (this->basis);
    (this->basis) = new VecSet< DistSVec<double, dim> >(podSize, this->domain.getNodeDistInfo());
    for (int iSnap=0; iSnap<podSize; ++iSnap) { 
      (*(this->basis))[iSnap] = (*Utrue)[iSnap]; 
    }
    delete Utrue;
    Utrue = NULL;
  
    if (this->columnSumsV) delete (this->columnSumsV);
    (this->columnSumsV) = new std::vector<double>;
    this->columnSumsV->clear();
    this->columnSumsV->resize(nTotSnaps, 0.0);
    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
      for (int jSnap=0; jSnap<nTotSnaps; ++jSnap) {
        (*(this->columnSumsV))[iSnap] += (*Vtrue)[jSnap][iSnap] * (this->cumulativeSnapWeights[jSnap]);
      }
    }  
    delete Vtrue;
    Vtrue = NULL;
    this->cumulativeSnapWeights.clear();
  
    if (this->sVals) delete (this->sVals);
    (this->sVals) = new std::vector<double>;
    this->sVals->resize(nTotSnaps);
    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
      (*this->sVals)[iSnap] = singVals[iSnap];
    } 
    singVals.clear();

    // output the basis and the update quantities
    this->outputClusteredBasis(iCluster, nTotSnaps, basisType); 
  
    this->com->barrier();

    if (strcmp(basisType,"sensitivity")==0)  break;
  }  

  delete [] (this->snapsInCluster);
  (this->snapsInCluster) = NULL;
  delete (this->clusterCenters);
  (this->clusterCenters) = NULL;

  if (this->nearestSnapsToCenters) {
    delete (this->nearestSnapsToCenters);
    (this->nearestSnapsToCenters) = NULL;
  }

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::SVD(VecSet< DistSVec<double, dim> >*& snapshots, VecSet< DistSVec<double, dim> > &Utrue,
    std::vector<double>& singularValues, FullM &Vtrue, int podMethod, int randMatDimension, int nPowerIts, bool computeV, bool testProbSVD) {

  double podTime = this->timer->getTime();

  if (podMethod==DataCompressionData::SCALAPACK_SVD) {
    scalapackSVD(snapshots, Utrue, singularValues, Vtrue, computeV);
  } else if (podMethod==DataCompressionData::PROBABILISTIC_SVD) {
    probabilisticSVDWrapper(snapshots, Utrue, singularValues, Vtrue, randMatDimension, nPowerIts, testProbSVD);
  } else if (podMethod==DataCompressionData::R_SVD) {
    rSVDWrapper(snapshots, Utrue, singularValues, Vtrue, testProbSVD);
  } else {
    this->com->fprintf(stderr, "*** Error: Unexpected POD method (EVD is not supported for nonlinear ROM prepro)\n");
    exit(-1);
  }

  this->timer->addPODTime(podTime);

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::scalapackSVD(VecSet< DistSVec<double, dim> >*& snapshots, VecSet< DistSVec<double, dim> > &Utrue,
    std::vector<double>& singularValuesVec, FullM &Vtrue, bool computeV) {
 
  #ifndef DO_SCALAPACK
    this->com->fprintf(stderr, "*** Error: Aero-F was not compiled with ScaLAPACK \n");
    exit(-1);
  #endif

  int nVecs=snapshots->numVectors();
  double* singularValuesArray = new double[nVecs];

  ParallelRom<dim> parallelRom( this->domain, this->com, this->domain.getNodeDistInfo());
  parallelRom.parallelSVD(*snapshots, Utrue, singularValuesArray, Vtrue, nVecs, computeV);

  this->com->broadcast(nVecs, singularValuesArray, 0);

  singularValuesVec.resize(nVecs);
  for (int iVec=0; iVec<nVecs; ++iVec) {
    singularValuesVec[iVec] = singularValuesArray[iVec];
  }
  delete [] singularValuesArray;

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::rSVDWrapper(VecSet< DistSVec<double, dim> >*& snapshots, VecSet< DistSVec<double, dim> > &Utrue,
    std::vector<double>& singularValues, FullM &Vtrue, bool testSVD) {
//  wrapper for R-SVD function
//  Note: moved core functionality to NonlinearRom.C since it proved widely applicable

  int nVecs = snapshots->numVectors();

  Utrue.resize(nVecs);
  for (int iVec=0; iVec<nVecs; ++iVec) {
    Utrue[iVec]=(*snapshots)[iVec];
  } 

  this->rSVD(Utrue, singularValues, Vtrue, testSVD);
  
}



//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::probabilisticSVDWrapper(VecSet< DistSVec<double, dim> >*& snapshots, VecSet< DistSVec<double, dim> > &Utrue,
    std::vector<double>& singularValues, FullM &Vtrue, int randMatDimension, int nPowerIts, bool testSVD) {
//  Note: moved core functionality to NonlinearRom.C since it proved widely applicable

//  Matlab code:
//    function [U,S,V] = probabilistic_svd(A,k,q)
//        
//    % Gaussian Random Matrix
//    Omega = randn(size(A,2),k);
//    % Power iteration
//    Y = A*Omega;
//    for i = 1:q
//        Y = A*(A'*Y);
//    end
//    [Q,~] = qr(Y,0);   % Orthogonalization
//    B = Q'*A;          % Projections
//    [tU,S,V]=svd(B,0); % SVD
//    U = Q*tU;          % Reconstruction

  int nVecs = snapshots->numVectors();

  Utrue.resize(nVecs);
  for (int iVec=0; iVec<nVecs; ++iVec) {
    Utrue[iVec]=(*snapshots)[iVec];
  } 

  std::vector<double> truncatedSingularValues;

  this->probabilisticSVD(Utrue, truncatedSingularValues, Vtrue, randMatDimension, nPowerIts, testSVD);
  
  //  The SVD only returns k non-zero singular values, which makes it strange to size bases based on singular energy.
  //  Assume the worst case scenario where final nVecs-k singular values do not decay...
  int k = Utrue.numVectors();
  singularValues.clear();
  singularValues.resize(nVecs,truncatedSingularValues[k-1]);
  for (int iVec=0; iVec<k; ++iVec) {
    singularValues[iVec] = truncatedSingularValues[iVec];
  }
  
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::preprocessForExactBasisUpdates() {

this->com->fprintf(stdout, "\nPreprocessing for fast exact basis updates (GNAT compatible)\n");
this->com->fprintf(stdout, " ... Note: This preprocessing step currently assumes that V^T * V = I\n");
double exactUpdatesPreproTime = this->timer->getTime();

  VecSet<DistSVec<double, dim> >* rob_i = NULL;
  VecSet<DistSVec<double, dim> >* rob_p = NULL;

  this->basisBasisProducts.resize(this->nClusters);
  this->basisUrefProducts.resize(this->nClusters);
  this->urefUrefProducts.resize(this->nClusters);

  if (preproForArbitraryUniformIC) {
    this->urefComponentwiseSums.resize(this->nClusters);
    this->basisComponentwiseSums.resize(this->nClusters);
  } else if (preproForInterpolatedIC) {
    this->urefMultiUicProducts.resize(this->nClusters);
    this->basisMultiUicProducts.resize(this->nClusters);
  } else {
    this->urefUicProducts.clear();
    this->urefUicProducts.resize(this->nClusters, 0.0);
    this->basisUicProducts.resize(this->nClusters);
  }

  for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
    //load rob_i
    this->readClusteredBasis(iCluster, "state");
    rob_i = new VecSet<DistSVec<double, dim> >(this->basis->numVectors(), this->domain.getNodeDistInfo());
    for (int iVec=0; iVec<this->basis->numVectors(); ++iVec)
      (*rob_i)[iVec] = (*this->basis)[iVec];
    this->basisBasisProducts[iCluster].resize(iCluster);
    for (int pCluster=0; pCluster<iCluster; ++pCluster) {  // assume each basis is orthogonal, no need to compute diagonal
      //load rob_p
      if (iCluster == pCluster) {
        rob_p = rob_i;
      } else {
        this->readClusteredBasis(pCluster, "state");
        rob_p = new VecSet<DistSVec<double, dim> >(this->basis->numVectors(), this->domain.getNodeDistInfo());
        for (int pVec=0; pVec<this->basis->numVectors(); ++pVec)
          (*rob_p)[pVec] = (*this->basis)[pVec];
      }
      this->basisBasisProducts[iCluster][pCluster].resize(rob_i->numVectors());
      for (int iVec=0; iVec<rob_i->numVectors(); ++iVec) {
        this->basisBasisProducts[iCluster][pCluster][iVec].clear();
        this->basisBasisProducts[iCluster][pCluster][iVec].resize(rob_p->numVectors(), 0.0);
        for (int pVec=0; pVec<rob_p->numVectors(); ++pVec) {
          this->basisBasisProducts[iCluster][pCluster][iVec][pVec] = (*rob_i)[iVec] * (*rob_p)[pVec];
        }
      }
      if (iCluster == pCluster) {
        rob_p = NULL;
      } else {
        delete rob_p;
        rob_p = NULL;
      }
    } // end pCluster loop (#1)
    this->basisUrefProducts[iCluster].resize(this->nClusters);
    this->urefUrefProducts[iCluster].clear();
    this->urefUrefProducts[iCluster].resize(iCluster+1, 0.0);
    // read Uref_i
    DistSVec<double, dim> Uref_i(this->domain.getNodeDistInfo());
    DistSVec<double, dim> Uref_p(this->domain.getNodeDistInfo());   
    this->readClusteredReferenceState(iCluster, "state");
    Uref_i = *(this->Uref);
    if (preproForArbitraryUniformIC) {
      this->urefComponentwiseSums[iCluster].clear();
      this->urefComponentwiseSums[iCluster].resize(dim, 0.0);
      double componentSums[dim];
      Uref_i.sum(componentSums);
      for (int iDim=0; iDim<dim; ++iDim) {
        this->urefComponentwiseSums[iCluster][iDim] = componentSums[iDim];
      }
      this->basisComponentwiseSums[iCluster].resize(rob_i->numVectors());
      for (int iVec=0; iVec<rob_i->numVectors(); ++iVec) {
        this->basisComponentwiseSums[iCluster][iVec].clear();
        this->basisComponentwiseSums[iCluster][iVec].resize(dim, 0.0);
        (*rob_i)[iVec].sum(componentSums);
        for (int iDim=0; iDim<dim; ++iDim) {
          this->basisComponentwiseSums[iCluster][iVec][iDim] = componentSums[iDim];
        }
      }
    } else if (preproForInterpolatedIC) {
      int nIC = this->multiUic->numVectors();
      this->urefMultiUicProducts[iCluster].resize(nIC);
      this->basisMultiUicProducts[iCluster].resize(nIC);
      for (int iIC=0; iIC<nIC; ++iIC) {
        this->urefMultiUicProducts[iCluster][iIC] = (*this->multiUic)[iIC] * Uref_i;
        this->basisMultiUicProducts[iCluster][iIC].resize(rob_i->numVectors()); 
        for (int iVec=0; iVec<rob_i->numVectors(); ++iVec) {
          this->basisMultiUicProducts[iCluster][iIC][iVec] = (*rob_i)[iVec] * (*this->multiUic)[iIC];
        }
      }
    } else {
      this->urefUicProducts[iCluster]= (*initialCondition) * Uref_i;
      this->basisUicProducts[iCluster].resize(rob_i->numVectors()); 
      for (int iVec=0; iVec<rob_i->numVectors(); ++iVec) {
        this->basisUicProducts[iCluster][iVec] = (*rob_i)[iVec] * (*initialCondition);
      }
    }
    for (int pCluster=0; pCluster<this->nClusters; ++pCluster) {
      if (pCluster == iCluster) {
        Uref_p = Uref_i;
      } else {
        this->readClusteredReferenceState(pCluster, "state");
        Uref_p = *(this->Uref);
      }
      this->basisUrefProducts[iCluster][pCluster].clear();
      this->basisUrefProducts[iCluster][pCluster].resize(rob_i->numVectors(), 0.0);
      for (int iVec=0; iVec<rob_i->numVectors(); ++iVec)
        this->basisUrefProducts[iCluster][pCluster][iVec] = (*rob_i)[iVec] * Uref_p;
      if (pCluster <= iCluster)
        this->urefUrefProducts[iCluster][pCluster] = Uref_i * Uref_p;
    } // end pCluster loop (#2)
    delete rob_i;
    rob_i=NULL;
  } // end iCluster loop

  // Basis Basis Products (rob_i^T * rob_p)
  // std::vector<std::vector<std::vector<std::vector<double> > > > basisBasisProducts;  // [iCluster][pCluster][:][:]
  this->outputClusteredInfoASCII(-1, "basisBasisProducts", NULL, NULL, NULL, &this->basisBasisProducts);

  // Basis Uref Products (rob_i^T * Uref_p)
  // std::vector<std::vector<std::vector<double> > > basisUrefProducts;  // [Cluster_Basis][Cluster_Uref][:]
  this->outputClusteredInfoASCII(-1, "basisUrefProducts", NULL, NULL, &this->basisUrefProducts);

  // Uref Uref Products
  // std::vector<std::vector<double> > urefUrefProducts; //[iCluster][jCluster] symmetric (lower triangular)
  this->outputClusteredInfoASCII(-1, "urefUrefProducts", NULL, &this->urefUrefProducts);

  if (preproForArbitraryUniformIC) {
    // Uref Componentwise Sums
    // std::vector<std::vector<double> > urefComponentwiseSums; //[iCluster][1:dim]
    this->outputClusteredInfoASCII(-1, "urefComponentwiseSums", NULL, &this->urefComponentwiseSums);

    // Basis Componentwise Sums
    // std::vector<std::vector<std::vector<double> > > basisComponentwiseSums;  // [iCluster][iVec][1:dim]
    this->outputClusteredInfoASCII(-1, "basisComponentwiseSums", NULL, NULL, &this->basisComponentwiseSums);
  } else if (preproForInterpolatedIC) {
    // Basis MultiUic Products
    // std::vector<std::vector<double> > basisMultiUicProducts;  // [iCluster][jIC][1:nPod] only precomputed if multiUic specified
    this->outputClusteredInfoASCII(-1, "basisMultiUicProducts", NULL, NULL, &this->basisMultiUicProducts);

    // Uref MultiUic Products
    // std::vector<double> urefMultiUicProducts; // [iCluster][jIC] only precomputed if multiUic specified
    this->outputClusteredInfoASCII(-1, "urefMultiUicProducts", NULL, &this->urefMultiUicProducts);
  } else {
    // Basis Uic Products
    // std::vector<std::vector<double> > basisUicProducts;  // [iCluster][1:nPod] only precomputed if Uic specified
    this->outputClusteredInfoASCII(-1, "basisUicProducts", NULL, &this->basisUicProducts);

    // Uref Uic Products
    // std::vector<double> urefUicProducts; // [iCluster] only precomputed if Uic specified
    this->outputClusteredInfoASCII(-1, "urefUicProducts", &this->urefUicProducts);
  }

  this->timer->addExactUpdatesPreproTime(exactUpdatesPreproTime);

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::preprocessForDistanceComparisons() {

  if (this->incrementalStateSnaps) {
    preprocessForDistanceComparisonsIncrements();
  } else {
    preprocessForDistanceComparisonsStandard();
  }

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::preprocessForDistanceComparisonsIncrements() {

  // only need V_i^T * (Center_p / || Center_p ||_p)
  this->com->fprintf(stdout, "\nPreprocessing for fast online cluster selection\n");
  double distCalcsTime = this->timer->getTime();

  this->basisNormalizedCenterProducts.resize(this->nClusters);

  if (this->clusterCenters==NULL) {
    this->readClusterCenters("centers");
  }

  DistSVec<double,dim> normalizedCenter(this->domain.getNodeDistInfo());
  for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
    this->readClusteredBasis(iCluster, "state");
    this->basisNormalizedCenterProducts[iCluster].resize(this->nClusters);
    for (int pCluster=0; pCluster<this->nClusters; ++pCluster) {
      this->basisNormalizedCenterProducts[iCluster][pCluster].clear();
      this->basisNormalizedCenterProducts[iCluster][pCluster].resize(this->basis->numVectors(), 0.0);
      normalizedCenter = (*this->clusterCenters)[pCluster];
      double norm = normalizedCenter.norm();
      if (norm>1e-16) normalizedCenter *= 1/norm; 
      for (int iVec=0; iVec<this->basis->numVectors(); ++iVec)   
        this->basisNormalizedCenterProducts[iCluster][pCluster][iVec] = (*this->basis)[iVec] * normalizedCenter;
    }
  }

  // Basis Center Products (rob_i^T * Center_p / || Center_p ||)
  // std::vector<std::vector<std::vector<double> > > basisNormalizedCenterProducts;  // [Cluster_Basis][Cluster_Center][:]
  this->outputClusteredInfoASCII(-1, "basisNormalizedCenterProducts", NULL, NULL, &this->basisNormalizedCenterProducts);

  this->timer->addDistCalcsPreproTime(distCalcsTime);

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::preprocessForDistanceComparisonsStandard() {

  this->com->fprintf(stdout, "\nPreprocessing for fast online cluster selection\n");
  double distCalcsTime = this->timer->getTime();

  readInitialCondition();
  
  this->readClusterCenters("centers");
  std::vector<std::vector<double> > clusterCenterNorms; // nClusters-by-1, or nClusters-by-(dim+1)
  clusterCenterNorms.resize(this->nClusters);
  
  if (preproForArbitraryUniformIC) {  // preprocess for an arbitrary uniform initial condition

    // loop through all cluster centers, store squared norm of center and component-wise sums of center
    for (int iCenter=0; iCenter<(this->nClusters); ++iCenter) {
      clusterCenterNorms[iCenter].clear();
      clusterCenterNorms[iCenter].reserve(dim+1);
      DistSVec<double,dim> difference(this->domain.getNodeDistInfo());
      double norm = ((*(this->clusterCenters))[iCenter]).norm();
      norm *= norm;
      (clusterCenterNorms[iCenter]).push_back(norm);
      double componentSums[dim];
      ((*(this->clusterCenters))[iCenter]).sum(componentSums);
      for (int iDim=0; iDim<dim; ++iDim) {
        (clusterCenterNorms[iCenter]).push_back(componentSums[iDim]);
      }
    } 

  } else { // preprocess for specified initial condition / multiUic

    for (int iCenter=0; iCenter<(this->nClusters); ++iCenter) {
        clusterCenterNorms[iCenter].clear();
        clusterCenterNorms[iCenter].resize(1);
        double centerNorm = ((*(this->clusterCenters))[iCenter]).norm();
        centerNorm *= centerNorm;
        clusterCenterNorms[iCenter][0]=centerNorm;
    }

  }
 
  this->outputCenterNorms(clusterCenterNorms);

  
  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
    productOfVectorAndCenterDifferences(iCluster, "referenceState");
  }

  productOfVectorAndCenterDifferences(-1, "initialCondition");

  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
    productOfBasisAndCenterDifferences(iCluster, "state");
    if (strcmp(this->krylovBasisName,"")!=0) productOfBasisAndCenterDifferences(iCluster, "krylov");
  }

  if (strcmp(this->sensitivityBasisName,"")!=0) productOfBasisAndCenterDifferences(-2, "sensitivity");

  this->timer->addDistCalcsPreproTime(distCalcsTime);

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::productOfBasisAndCenterDifferences(int iCluster, const char* basisType) {
  // helper function for fast distance calculation preprocessing:
  // computes w_(m,p) = 2 * basis^T *(center_p - center_m) for 0<=p<m<nClusters
  // note that w_(m,p) = -w_(p,m)

  // note that in this case 0 <= p < m < nClusters, which differs from Amsallem et al., INJME 2012
  // (this allows for easy indexing [m][p] without any unnecessary storage)

  this->readClusteredBasis(iCluster, basisType);
  int nPodVecs = this->basis->numVectors();
  std::vector<std::vector<std::vector<double> > > w;  // collection of vectors: [m][p][iPodVec]
  w.resize(this->nClusters);

  // w[0][0] is a vector of zeros, no need to compute this
    
  for (int mCenter=1; mCenter<this->nClusters; ++mCenter) {
    for (int pCenter=0; pCenter<mCenter; ++pCenter) {
      DistSVec<double,dim> difference(this->domain.getNodeDistInfo());
      difference = (*(this->clusterCenters))[pCenter] - (*(this->clusterCenters))[mCenter]; 
      std::vector<double> product;
      product.reserve(nPodVecs);
      for (int iVec=0; iVec<nPodVecs; ++iVec) { 
        product.push_back( 2.0 * ((*(this->basis))[iVec] * difference));
      }
      w[mCenter].push_back(product); 
    }
  }
 
  this->outputClusteredInfoASCII(iCluster, basisType, NULL, NULL, &w);

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::productOfVectorAndCenterDifferences(int iCluster, const char* vecType) {
  // helper function for fast distance calculation preprocessing:
  // computes w_(m,p) = 2 * vector^T *(center_p - center_m) for 0<=p<m<nClusters
  // note that w_(m,p) = -w_(p,m)

  // note that in this case 0 <= p < m < nClusters, which differs from Amsallem et al., INJME 2012
  // (this allows for easy indexing [m][p] without any unnecessary storage)

  if (preproForArbitraryUniformIC) return;

  if (this->nClusters == 1) return;

  if (this->clusterCenters==NULL) this->readClusterCenters("state");

  if (preproForInterpolatedIC && strcmp(vecType,"initialCondition")==0) {

    std::vector<std::vector<std::vector<double> > > result;  // collection of vectors: [m][p][iIC]
    result.resize(this->nClusters);
    // result[0][0] is a vector of zeros, no need to compute this
    int nIC = this->multiUic->numVectors();
    for (int mCenter=1; mCenter<this->nClusters; ++mCenter) {
      result[mCenter].resize(mCenter);
      for (int pCenter=0; pCenter<mCenter; ++pCenter) {
        DistSVec<double,dim> difference(this->domain.getNodeDistInfo());
        difference = (*(this->clusterCenters))[pCenter] - (*(this->clusterCenters))[mCenter];
        result[mCenter][pCenter].resize(nIC);
        for (int iIC=0; iIC<nIC; ++iIC) {
          double tmp = 2.0 * ((*this->multiUic)[iIC] * difference);
          result[mCenter][pCenter][iIC]=tmp;
        }
      }
    }
    this->outputClusteredInfoASCII(iCluster, "multiInitialCondition", NULL, NULL, &result);

    // also need to output multiUicMultiUicProducts
    std::vector<std::vector<double> > multiUicMultiUicProducs;  // collection of vectors: [iIC][jIC]
    multiUicMultiUicProducs.resize(nIC);
    for (int iIC=0; iIC<nIC; ++iIC) {
      multiUicMultiUicProducs[iIC].resize(iIC+1);
      for (int jIC=0; jIC<=iIC; ++jIC) {
        multiUicMultiUicProducs[iIC][jIC] = (*this->multiUic)[iIC] * (*this->multiUic)[jIC];
      }
    }
    this->outputClusteredInfoASCII(iCluster, "multiUicMultiUicProducts", NULL, &multiUicMultiUicProducs);

  } else {

    DistSVec<double,dim>* vec;
    vec = NULL;
  
    if (strcmp(vecType,"referenceState")==0) {
      this->readClusteredReferenceState(iCluster, "state");
      vec = this->Uref;
    } else if (strcmp(vecType, "initialCondition")==0 && initialCondition!=NULL) {
      vec = initialCondition;
    } else {
      this->com->fprintf(stderr, "*** Error: unanticipated vector type '%s' encountered (or initialCondition=NULL)", vecType);
      exit(-1);
    }
  
    std::vector<std::vector<double> > result;  // collection of vectors: [m][p]
    result.resize(this->nClusters);
  
    // result[0][0] is a vector of zeros, no need to compute this
  
    for (int mCenter=1; mCenter<this->nClusters; ++mCenter) {
      for (int pCenter=0; pCenter<mCenter; ++pCenter) {
        DistSVec<double,dim> difference(this->domain.getNodeDistInfo());
        difference = (*(this->clusterCenters))[pCenter] - (*(this->clusterCenters))[mCenter]; 
        double tmp = 2.0 * ((*vec) * difference);
        result[mCenter].push_back(tmp); 
      }
    }
  
    this->outputClusteredInfoASCII(iCluster, vecType, NULL, &result);

  }

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::readInitialCondition() {

  // check if this function has already been called
  if (initialCondition || preproForArbitraryUniformIC || preproForInterpolatedIC) return;

  if (strcmp(this->ioData->input.solutions,"")!=0) {

    // read user-specified initial condition that will be used for the online ROM simulation
    initialCondition = new DistSVec<double,dim>( this->domain.getNodeDistInfo() );

    char* icFile = new char[strlen(this->ioData->input.prefix) + strlen(this->ioData->input.solutions) + 1];
    sprintf(icFile, "%s%s", this->ioData->input.prefix, this->ioData->input.solutions);

    double tmp;
    bool status = this->domain.readVectorFromFile(icFile, 0, &tmp, *initialCondition);

    if (!status) {
      this->com->fprintf(stderr, "*** Error: unable to read vector from file %s\n", icFile);
      exit(-1);
    }

    delete [] icFile;
  
  } else if (strcmp(this->ioData->input.multiSolutionsParams,"")!=0) {

    // read file, store solutions and parameters
    FILE *paramsFile = fopen(this->ioData->input.multiSolutionsParams,"r");
    if (!paramsFile)  {
      this->com->fprintf(stderr, "*** Error: No solution data FILES in %s\n", this->ioData->input.multiSolutionsParams);
      exit (-1);
    }
    int nData, _n, nParams;
    _n = fscanf(paramsFile, "%d",&nData);
    this->com->fprintf(stdout, " ... reading %d solutions for interpolation\n",nData);
  
    _n = fscanf(paramsFile, "%d",&nParams);
    this->com->fprintf(stdout, " ... %d-dimensional parameter space\n",nParams);
 
    this->multiUic = new VecSet<DistSVec<double,dim> >(nData, this->domain.getNodeDistInfo() );
 
    char solnFile[500];
    for (int iData=0; iData < nData; ++iData) {
      // read solution
      _n = fscanf(paramsFile, "%s", solnFile);
      double tmp;
      bool status = this->domain.readVectorFromFile(solnFile, 0, &tmp, (*this->multiUic)[iData]);

      if (!status) {
        this->com->fprintf(stderr, "*** Error: unable to read vector from file %s\n", solnFile);
        exit(-1);
      }

      // read parameters
      for (int iParam=0; iParam < nParams; ++iParam) {
        _n = fscanf(paramsFile, "%lf", &tmp);
      } 
    }
    fclose(paramsFile);

    preproForInterpolatedIC = true;

  } else {

    this->com->fprintf(stdout, "\n ... no initial condition specified; preprocessing for a uniform initial condition\n");
    preproForArbitraryUniformIC = true;

  }

}

////----------------------------------------------------------------------------------
//
//template<int dim>
//void NonlinearRomDatabaseConstruction<dim>::localRelProjError() {
//
//  double projErrorTime = this->timer->getTime();
//
//  this->readClusterCenters("centers");
//  delete (this->clusterCenters);
//  (this->clusterCenters) = NULL;
// 
//  nSnapshotFiles = this->readSnapshotFiles("projError", true);
//  int nTotSnaps = this->snap->numVectors();
// 
//  projErrorLog = new VecSet<Vec<double> >((this->nClusters),nTotSnaps);
//
//  // if computing residuals
//  //  if (projError->basisUpdates!=RelativeProjectionErrorData::UPDATES_OFF)
//  //  ImplicitPGTsDesc<dim> tsDesc(ioData, geoSource, &domain);
//
//  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
//
//    switch (projError->relProjError) {
//      case (RelativeProjectionErrorData::REL_PROJ_ERROR_STATE):
//        this->readClusteredBasis(iCluster, "state", true);
//        if (projError->basisUpdates!=RelativeProjectionErrorData::UPDATES_OFF &&
//            this->snapRefState!=NULL) 
//          this->updateBasis(iCluster, *(this->snapRefState));
//        if (projError->krylov.include) this->appendNonStateDataToBasis(iCluster,"krylov",true);
//        if (projError->sensitivity.include) this->appendNonStateDataToBasis(iCluster,"sensitivity",true);
//        break;
//      case (RelativeProjectionErrorData::REL_PROJ_ERROR_RESIDUAL):
//        this->readClusteredBasis(iCluster, "residual", true);
//        break;
//      case (RelativeProjectionErrorData::REL_PROJ_ERROR_JACACTION):
//        this->readClusteredBasis(iCluster, "jacAction", true);
//        break;
//      default:
//        exit (-1);
//    }
//
//    int nPodVecs = this->basis->numVectors();
//
//    this->com->fprintf(stdout, "\nCalculating relative projection error for basis %d\n", iCluster);
//
//    // basis^T * snapshots    
//    this->com->fprintf(stdout, " ... tmp = basis^T * snapshots\n");
//    Vec<double> tmpVec(nPodVecs);
//    VecSet<Vec<double> >* tmpVecSet = new VecSet<Vec<double> >(nTotSnaps, nPodVecs);
//    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
//      for (int iVec = 0; iVec < nPodVecs; iVec++){
//        tmpVec[iVec] = (*(this->basis))[iVec] * (*(this->snap))[iSnap];
//      }
//      (*tmpVecSet)[iSnap] = tmpVec;
//    }
//
//    // basis * result 
//    this->com->fprintf(stdout, " ... projected snaps = basis * tmp\n");
//    VecSet<DistSVec<double, dim> >* projectedSnaps = new VecSet<DistSVec<double, dim> >(nTotSnaps, this->domain.getNodeDistInfo());
//    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
//      (*projectedSnaps)[iSnap] = 0.0;
//      for (int iVec = 0; iVec < nPodVecs; iVec++)
//        (*projectedSnaps)[iSnap] += ((*tmpVecSet)[iSnap])[iVec] * (*(this->basis))[iVec];
//    }
//    delete tmpVecSet;
//    tmpVecSet = NULL;
//
//    // option to output postprocesed projected states (projectedSnaps + snapRefState if snapRefState is defined)
//    
//    if (projError->postProProjectedStates == RelativeProjectionErrorData::POST_PRO_ON ) { //||
//       // projError->outputResidualOfProjStates == RelativeProjectionErrorData::CALC_RESIDUALS_ON) {
//  
//      TsDesc<dim>* tsDesc = new TsDesc<dim>(*(this->ioData), geoSource, &(this->domain));    
//
//      if (projError->postProProjectedStates == RelativeProjectionErrorData::POST_PRO_ON) {
//
//        DistSVec<double,dim> outVec(this->domain.getNodeDistInfo());
//
//        for (int iSnap=0;iSnap<nTotSnaps;++iSnap) {
//          if (projError->subtractRefSol) {
//             this->readReferenceState(); 
//             outVec = *(this->snapRefState);
//             outVec += (*projectedSnaps)[iSnap];
//             delete (this->snapRefState);
//             (this->snapRefState) = NULL;
//          } else {
//             outVec = (*projectedSnaps)[iSnap];
//          }
//          tsDesc->performPostProForState(outVec);
//        }
//      }
//
//      // option to output spatial residual of projected states (projectedSnaps + snapRefState if snapRefState is defined)
//      //if (projError->outputResidualOfProjStates == RelativeProjectionErrorData::CALC_RESIDUALS_ON) {
//       // this->spaceOp->computeResidual(*this->X, *this->A, Q, *R, this->timeState);
//      //}
//      if (tsDesc) delete tsDesc;
//    }
//
//    // snaps - projSnaps
//    this->com->fprintf(stdout, " ... difference = originalSnaps - projectedSnaps\n");
//    VecSet<DistSVec<double, dim> >* snapDifference = new VecSet<DistSVec<double, dim> >(nTotSnaps, this->domain.getNodeDistInfo());
//    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
//      (*snapDifference)[iSnap] = (*(this->snap))[iSnap] - (*projectedSnaps)[iSnap];
//    }
//
//    delete projectedSnaps;
//    projectedSnaps = NULL;
//   
//    // compute l2 norm for each 
//    Vec<double> numeratorNorms(nTotSnaps);
//    Vec<double> denominatorNorms(nTotSnaps);
//
//    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
//      numeratorNorms[iSnap] = (*snapDifference)[iSnap].norm();
//      denominatorNorms[iSnap] = (*(this->snap))[iSnap].norm();
//    }
//
//    delete snapDifference;
//    snapDifference = NULL;
//
//    this->com->fprintf(stdout, " ... relProjError = difference.norm / originalSnap.norm\n");
//    Vec<double> relProjError(nTotSnaps);
//    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
//    if (denominatorNorms[iSnap] != 0) {
//        relProjError[iSnap] = numeratorNorms[iSnap]/denominatorNorms[iSnap];
//      } else {
//        relProjError[iSnap] = 1;
//      }
//    }
//
//    (*projErrorLog)[iCluster] = relProjError;
//
//    delete (this->basis);
//    (this->basis) = NULL;
//
//  }  
//
//  delete [] (this->snapsInCluster);
//  (this->snapsInCluster) = NULL;
//  delete (this->snap);
//  (this->snap) = NULL;
//
//  writeProjErrorToDisk();
//
//  this->timer->addProjErrorTime(projErrorTime);
//
//}
//
//
//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::localRelProjError() {

  double projErrorTime = this->timer->getTime();

  this->readClusterCenters("centers");
  delete (this->clusterCenters);
  (this->clusterCenters) = NULL;
 
  nSnapshotFiles = this->readSnapshotFiles("projError", false);
  int nTotSnaps = this->snap->numVectors();
 
  projErrorLog = new VecSet<Vec<double> >((this->nClusters),nTotSnaps);

  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
    switch (projError->relProjError) {
      case (RelativeProjectionErrorData::REL_PROJ_ERROR_STATE):
        this->readClusteredBasis(iCluster, "state", true);
        if (projError->basisUpdates!=RelativeProjectionErrorData::UPDATES_OFF) {
          this->updateBasis(iCluster, (*(this->snap))[0]);
        } else {
          this->readClusteredReferenceState(iCluster, "state");
          if (projError->useFirstStateAsRefStateForIncrBasis) {
            if (this->Uref->norm()>1e-15)  {
              this->com->fprintf(stderr, "*** Error: UseFirstStateAsRefStateForIncrementalBasis is set to True, but the stored reference state for this basis is nonzero.  This shouldn't be the case for a basis built with incremental snapshots.  Exiting.\n");
              exit(-1);
            }
            *(this->Uref) = (*(this->snap))[0];
          }
        }
        if (projError->krylov.include) this->com->fprintf(stdout, "*** Warning: Krylov vecs will not be appended for reconstruction error sweep\n");
        if (projError->sensitivity.include) this->com->fprintf(stdout, "*** Warning: Sensitivities will not be appended for reconstruction error sweep\n");
        break;
      case (RelativeProjectionErrorData::REL_PROJ_ERROR_RESIDUAL):
        this->readClusteredBasis(iCluster, "residual", true);
        if (this->Uref) delete this->Uref;
        this->Uref = new DistSVec<double, dim>(this->domain.getNodeDistInfo());
        *(this->Uref) = 0.0;
        break;
      case (RelativeProjectionErrorData::REL_PROJ_ERROR_JACACTION):
        this->readClusteredBasis(iCluster, "jacAction", true);
        if (this->Uref) delete this->Uref;
        this->Uref = new DistSVec<double, dim>(this->domain.getNodeDistInfo());
        *(this->Uref) = 0.0;
        break;
      default:
        exit (-1);
    }

    int nPodVecs = this->basis->numVectors();

    DistSVec<double,dim> referencedSnap(this->domain.getNodeDistInfo());
    Vec<double> denominatorNorms(nTotSnaps);
    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
      referencedSnap = (*(this->snap))[iSnap] - *(this->Uref);
      denominatorNorms[iSnap] = referencedSnap.norm();
    }

    this->com->fprintf(stdout, "\nCalculating affine approximation error using basis %d\n", iCluster);

    // basis^T * snapshots    
    this->com->fprintf(stdout, " ... tmp = basis^T * snapshots\n");
    Vec<double> tmpVec(nPodVecs);
    VecSet<Vec<double> >* tmpVecSet = new VecSet<Vec<double> >(nTotSnaps, nPodVecs);
    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
      for (int iVec = 0; iVec < nPodVecs; iVec++){
        tmpVec[iVec] = (*(this->basis))[iVec] * ((*(this->snap))[iSnap] - *(this->Uref));
      }
      (*tmpVecSet)[iSnap] = tmpVec;
    }

    // basis * result 
    this->com->fprintf(stdout, " ... projected snaps = basis * tmp\n");
    VecSet<DistSVec<double, dim> >* projectedSnaps = new VecSet<DistSVec<double, dim> >(nTotSnaps, this->domain.getNodeDistInfo());
    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
      (*projectedSnaps)[iSnap] = 0.0;
      for (int iVec = 0; iVec < nPodVecs; iVec++)
        (*projectedSnaps)[iSnap] += ((*tmpVecSet)[iSnap])[iVec] * (*(this->basis))[iVec];
    }
    delete tmpVecSet;
    tmpVecSet = NULL;

    // option to output postprocesed projected states (projectedSnaps + snapRefState if snapRefState is defined)
    
    if (projError->postProProjectedStates == RelativeProjectionErrorData::POST_PRO_ON ) { //||
       // projError->outputResidualOfProjStates == RelativeProjectionErrorData::CALC_RESIDUALS_ON) {
      TsDesc<dim>* tsDesc = new TsDesc<dim>(*(this->ioData), geoSource, &(this->domain));    
      if (projError->postProProjectedStates == RelativeProjectionErrorData::POST_PRO_ON) {
        DistSVec<double,dim> outVec(this->domain.getNodeDistInfo());
        for (int iSnap=0;iSnap<nTotSnaps;++iSnap) {
          if (projError->subtractRefSol) {
             this->readReferenceState(); 
             outVec = *(this->snapRefState);
             outVec += (*projectedSnaps)[iSnap];
             delete (this->snapRefState);
             (this->snapRefState) = NULL;
          } else {
             outVec = (*projectedSnaps)[iSnap] + *(this->Uref);
          }
          tsDesc->performPostProForState(outVec, iSnap);
        }
      }

      // option to output spatial residual of projected states (projectedSnaps + snapRefState if snapRefState is defined)
      //if (projError->outputResidualOfProjStates == RelativeProjectionErrorData::CALC_RESIDUALS_ON) {
       // this->spaceOp->computeResidual(*this->X, *this->A, Q, *R, this->timeState);
      //}
      if (tsDesc) delete tsDesc;
    }

    // snaps - projSnaps
    this->com->fprintf(stdout, " ... difference = originalSnaps - projectedSnaps\n");
    VecSet<DistSVec<double, dim> >* snapDifference = new VecSet<DistSVec<double, dim> >(nTotSnaps, this->domain.getNodeDistInfo());
    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
      (*snapDifference)[iSnap] = (*(this->snap))[iSnap] - *(this->Uref) - (*projectedSnaps)[iSnap];
    }

    delete projectedSnaps;
    projectedSnaps = NULL;
   
    // compute l2 norm for each 
    Vec<double> numeratorNorms(nTotSnaps);

    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
      numeratorNorms[iSnap] = (*snapDifference)[iSnap].norm();
    }

    delete snapDifference;
    snapDifference = NULL;

    this->com->fprintf(stdout, " ... relProjError = difference.norm / originalSnap.norm\n");
    Vec<double> relProjError(nTotSnaps);
    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
    if (denominatorNorms[iSnap] != 0) {
        relProjError[iSnap] = numeratorNorms[iSnap]/denominatorNorms[iSnap];
      } else {
        relProjError[iSnap] = 1;
      }
    }

    (*projErrorLog)[iCluster] = relProjError;

    delete (this->basis);
    (this->basis) = NULL;

  }  

  delete [] (this->snapsInCluster);
  (this->snapsInCluster) = NULL;
  delete (this->snap);
  (this->snap) = NULL;

  writeProjErrorToDisk();

  this->timer->addProjErrorTime(projErrorTime);

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::writeProjErrorToDisk()  {

  // output projError as ASCII file in top level of ROM database
  FILE *projErrorFile;
  char *fullProjErrorName = new char[strlen(this->databasePrefix) + strlen(this->databaseName) + 1 + strlen(this->projErrorName) + 1];
  sprintf(fullProjErrorName, "%s%s/%s", this->databasePrefix, this->databaseName, this->projErrorName);
  this->com->fprintf(stdout, "\nWriting projection error to disk: %s \n", fullProjErrorName);

  int nTotSnaps = (*projErrorLog)[0].size();

  if (this->com->cpuNum() == 0) {  
    projErrorFile = fopen(fullProjErrorName, "wt");
    this->com->fprintf(projErrorFile,"Snapshot# ErrorUsingEachBasis\n");

    for (int iSnap=0; iSnap<nTotSnaps; ++iSnap) {
      this->com->fprintf(projErrorFile,"%d ", iSnap);
      for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
      this->com->fprintf(projErrorFile,"%e ", ((*projErrorLog)[iCluster])[iSnap]);
      }
      this->com->fprintf(projErrorFile,"\n");
    }

    fclose (projErrorFile);
  }

  delete [] fullProjErrorName;
  fullProjErrorName = NULL;
  delete projErrorLog;
  projErrorLog = NULL;

  this->com->barrier();

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::writeProjErrorSweepToDisk(std::vector<std::vector<int> > robSizes,
                                                                      std::vector<std::vector<double> > projErrors)  {

  // output projError as ASCII file in top level of ROM database
  for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {

    FILE *projErrorSweepFile;
    int extraDigits = (iCluster==0) ? 1 : floor(log10((double)iCluster));
    char *fullProjErrorSweepName = new char[strlen(this->databasePrefix) + strlen(this->databaseName) + 1 + 
                                            strlen(this->projErrorName) + extraDigits + 1];
    sprintf(fullProjErrorSweepName, "%s%s/%s%d", this->databasePrefix, this->databaseName, this->projErrorName, iCluster);
    this->com->fprintf(stdout, "\nWriting projection error sweep info to disk: %s \n", fullProjErrorSweepName);
  
    int nSweeps = projErrors[iCluster].size();
  
    if (this->com->cpuNum() == 0) {  
      projErrorSweepFile = fopen(fullProjErrorSweepName, "wt");
      this->com->fprintf(projErrorSweepFile,"#BasisVecs FrobeniusNormOfError\n");
  
      for (int iSweep=0; iSweep<nSweeps; ++iSweep)
        this->com->fprintf(projErrorSweepFile,"%d %e\n", robSizes[iCluster][iSweep], projErrors[iCluster][iSweep]);
  
      fclose (projErrorSweepFile);
    }
  
    delete [] fullProjErrorSweepName;
    fullProjErrorSweepName = NULL;
  
    this->com->barrier();

  }

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::localRelProjErrorSweep() {

  double projErrorTime = this->timer->getTime();

  this->readClusterCenters("centers");
  delete (this->clusterCenters);
  (this->clusterCenters) = NULL;
 
  nSnapshotFiles = this->readSnapshotFiles("projError", false);

  int nTotSnaps = this->snap->numVectors();
  //normalize snapshots
  //for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
  //  double norm = ((*(this->snap))[iSnap].norm());
  //  if (norm>1e-16)  (*(this->snap))[iSnap] *= 1/norm;
  //}

  std::vector<std::vector<int> > robSizes; // [basis][sweep]
  std::vector<std::vector<double> > projErrors; // [basis][sweep]
  robSizes.resize(this->nClusters);
  projErrors.resize(this->nClusters);

  for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {

    robSizes[iCluster].clear();
    projErrors[iCluster].clear();

    switch (projError->relProjError) {
      case (RelativeProjectionErrorData::REL_PROJ_ERROR_STATE):
        this->readClusteredBasis(iCluster, "state", true);
        if (projError->basisUpdates!=RelativeProjectionErrorData::UPDATES_OFF) {
          this->updateBasis(iCluster, (*(this->snap))[0]);
        } else {
          this->readClusteredReferenceState(iCluster, "state");
          if (projError->useFirstStateAsRefStateForIncrBasis) {
            if (this->Uref->norm()>1e-15)  {
              this->com->fprintf(stderr, "*** Error: UseFirstStateAsRefStateForIncrementalBasis is set to True, but the stored reference state for this basis is nonzero.  This shouldn't be the case for a basis built with incremental snapshots.  Exiting.\n"); 
              exit(-1);
            }
            *(this->Uref) = (*(this->snap))[0];
          }
        }
        if (projError->krylov.include) this->com->fprintf(stdout, "*** Warning: Krylov vecs will not be appended for reconstruction error sweep\n");
        if (projError->sensitivity.include) this->com->fprintf(stdout, "*** Warning: Sensitivities will not be appended for reconstruction error sweep\n");
        break;
      case (RelativeProjectionErrorData::REL_PROJ_ERROR_RESIDUAL):
        this->readClusteredBasis(iCluster, "residual", true);
        if (this->Uref) delete this->Uref;
        this->Uref = new DistSVec<double, dim>(this->domain.getNodeDistInfo());
        *(this->Uref) = 0.0;
        break;
      case (RelativeProjectionErrorData::REL_PROJ_ERROR_JACACTION):
        this->readClusteredBasis(iCluster, "jacAction", true);
        if (this->Uref) delete this->Uref;
        this->Uref = new DistSVec<double, dim>(this->domain.getNodeDistInfo());
        *(this->Uref) = 0.0;
        break;
      default:
        exit (-1);
    }

    this->com->fprintf(stdout, "\nStarting snapshot reconstruction error sweep for basis %d\n", iCluster);

    int nPodVecs = this->basis->numVectors();
    double sweepFreq = (double)projError->sweepFreq;
    sweepFreq = (sweepFreq>nPodVecs) ? nPodVecs : sweepFreq;
    int nSweeps = ceil(((double)nPodVecs)/(sweepFreq));

    // calculate frobenius norm of training state matrix
    double stateMatrixFroNorm = 0.0;
    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
      stateMatrixFroNorm += pow((*(this->snap))[iSnap].norm(),2.0);
    }
    stateMatrixFroNorm = pow(stateMatrixFroNorm,0.5);

    // form snapshot matrix by subtracting reference state from each training state
    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
      (*(this->snap))[iSnap] -= *(this->Uref);
    }

    // basis^T * snapshots    
    this->com->fprintf(stdout, " ... tmp = basis^T * snapshots\n");
    Vec<double> tmpVec(nPodVecs);
    VecSet<Vec<double> >* tmpVecSet = new VecSet<Vec<double> >(nTotSnaps, nPodVecs);
    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
      for (int iVec = 0; iVec < nPodVecs; iVec++){
        tmpVec[iVec] = (*(this->basis))[iVec] * (*(this->snap))[iSnap];
      }
      (*tmpVecSet)[iSnap] = tmpVec;
    }

    VecSet<DistSVec<double, dim> >* snapDifference = new VecSet<DistSVec<double, dim> >(nTotSnaps, this->domain.getNodeDistInfo());
    for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
      (*snapDifference)[iSnap] = (*(this->snap))[iSnap];
    }

    int nPodVecsPrev = 0;
    int nPodVecsCurrent = sweepFreq;
 
    for (int iSweep = 0; iSweep < nSweeps; ++iSweep) {

      this->com->fprintf(stdout, "\nCalculating snapshot reconstruction error for basis %d with %d basis vectors\n", iCluster, nPodVecsCurrent);
   
      projErrors[iCluster].push_back(0.0);
      robSizes[iCluster].push_back(nPodVecsCurrent);
      for (int iSnap = 0; iSnap < nTotSnaps; iSnap++) {
        for (int iVec = nPodVecsPrev; iVec < nPodVecsCurrent; iVec++) {
          (*snapDifference)[iSnap] -= ((*tmpVecSet)[iSnap])[iVec] * (*(this->basis))[iVec];
        }
        projErrors[iCluster][iSweep] += (*snapDifference)[iSnap] * (*snapDifference)[iSnap];
      }
      projErrors[iCluster][iSweep] = pow(projErrors[iCluster][iSweep],0.5); // frobenius norm  
      projErrors[iCluster][iSweep] = (stateMatrixFroNorm>0) ? projErrors[iCluster][iSweep]/stateMatrixFroNorm : 1.0 ;

      nPodVecsPrev=nPodVecsCurrent;
      nPodVecsCurrent=nPodVecsPrev+sweepFreq;
      if (nPodVecsCurrent>nPodVecs) nPodVecsCurrent=nPodVecs;
    }

    delete snapDifference;
    snapDifference = NULL;
  
    delete (this->basis);
    (this->basis) = NULL;

    delete tmpVecSet;
    tmpVecSet = NULL;

  }  

  // output
  writeProjErrorSweepToDisk(robSizes, projErrors);

  delete [] (this->snapsInCluster);
  (this->snapsInCluster) = NULL;
  delete (this->snap);
  (this->snap) = NULL;

  this->timer->addProjErrorTime(projErrorTime);

}



//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::computeClassicalMultiDimensionalScaling()  {
#ifdef DO_MODAL

  this->com->fprintf(stdout, "\nEntering classical multi-dimensional scaling routine\n");

  //form distance matrix (squared distances)
  std::vector<std::vector<double> > distanceMat;
  std::vector<std::vector<double> > tmpMat;

  //DEBUGGING
 // int nSnaps = 10;
 // distanceMat.resize(nSnaps);
 // tmpMat.resize(nSnaps);
 // for (int iSnap=0;iSnap<nSnaps;++iSnap) {
 //     distanceMat[iSnap].resize(iSnap+1);
 //     tmpMat[iSnap].resize(iSnap+1);
 // }

 // for (int i=0;i<nSnaps;++i) {
 //  for (int j=0;j<nSnaps;++j) {
 //      if (abs(i-j)<((nSnaps/2)+1)) {
 //         distanceMat[i][j] = double(abs(j-i));
 //      } else {
  //         distanceMat[i][j] = double(nSnaps-abs(j-i));
 //      }
 //  }
 // }

  int nSnaps = this->snap->numVectors();
  this->com->fprintf(stdout, " ... calculating symmetric %dx%d distance matrix\n", nSnaps,nSnaps);
  this->com->fprintf(stdout, " ... (this might take a few minutes)\n");
  distanceMat.resize(nSnaps);
  tmpMat.resize(nSnaps);
  for (int iSnap=0;iSnap<nSnaps;++iSnap) {
      distanceMat[iSnap].resize(iSnap+1);
      tmpMat[iSnap].resize(iSnap+1);
  }

  for (int iSnap=0;iSnap<nSnaps;++iSnap) {
    for (int jSnap=0;jSnap<iSnap;++jSnap) {
      distanceMat[iSnap][jSnap] = this->distanceFull((*(this->snap))[iSnap], (*(this->snap))[jSnap]);
      distanceMat[iSnap][jSnap] *= distanceMat[iSnap][jSnap];
    }
    distanceMat[iSnap][iSnap] = 0.0; //diagonals are all zero 
  }

  //output distance matrix (matlab is useful for visualizing in 3D)
  this->outputClusteredInfoASCII(-1, "distanceMatrix", NULL, &distanceMat);

  //apply "double centering" to distance matrix (-.5*H*D*H)
  this->com->fprintf(stdout, " ... applying double centering transformation to distance matrix\n", nSnaps,nSnaps);
  for (int iSnap=0;iSnap<nSnaps;++iSnap) {
    for (int jSnap=0;jSnap<=iSnap;++jSnap) {
      tmpMat[iSnap][jSnap] = 0.0;
      for (int pSnap=0;pSnap<iSnap;++pSnap) {
        tmpMat[iSnap][jSnap] += distanceMat[iSnap][pSnap];
      }
      for (int qSnap=iSnap;qSnap<nSnaps;++qSnap) {
        tmpMat[iSnap][jSnap] += distanceMat[qSnap][iSnap];
      }
      tmpMat[iSnap][jSnap] /= double(nSnaps);
      tmpMat[iSnap][jSnap] = distanceMat[iSnap][jSnap] - tmpMat[iSnap][jSnap];
    }
  }
  for (int iSnap=0;iSnap<nSnaps;++iSnap) {
    for (int jSnap=0;jSnap<=iSnap;++jSnap) {
      distanceMat[iSnap][jSnap] = 0.0;
      for (int pSnap=0;pSnap<jSnap;++pSnap) {
        distanceMat[iSnap][jSnap] += tmpMat[iSnap][pSnap];
      }
      for (int qSnap=jSnap;qSnap<nSnaps;++qSnap) {
        distanceMat[iSnap][jSnap] += tmpMat[qSnap][jSnap];
      }
      distanceMat[iSnap][jSnap] /= double(nSnaps);
      distanceMat[iSnap][jSnap] = tmpMat[iSnap][jSnap] - distanceMat[iSnap][jSnap];
      distanceMat[iSnap][jSnap] *= -.5;
    }
  }

  if (this->com->cpuNum() == 0) {
  
    double* symMat = new double[((nSnaps*nSnaps)+nSnaps)/2]; //symmetric matrix
    int count=0;
    for (int jSnap=0;jSnap<nSnaps;++jSnap) {
      for (int iSnap=jSnap;iSnap<nSnaps;++iSnap) {
        symMat[count] = distanceMat[iSnap][jSnap];
        count++;
      }
    }
  
    for (int iSnap=0;iSnap<nSnaps;++iSnap) {
      tmpMat[iSnap].clear();
      distanceMat[iSnap].clear();
    }
    tmpMat.clear();
    distanceMat.clear();
   
    //extract eigenpairs corresponding to 2 largest positive eigenvalues (using ARPACK++)
    this->com->fprintf(stdout, " ... computing eigendecomposition of distance matrix using ARPACK\n");
    int numEigVals = nSnaps-1;
    double* eVals = new double[numEigVals];
    double* eVecs = new double[numEigVals*nSnaps];
    ARdsSymMatrix<double> arpackMat(nSnaps, symMat, 'L'); // real dense symmetric matrix
    arpackMat.FactorA();
    ARluSymStdEig<double> eigProb(numEigVals, arpackMat, "LM");
    int nconv = eigProb.EigenValVectors(eVecs, eVals);
  
    int eig1 = 0;
    int eig2 = 0;
  
    this->com->fprintf(stdout, "\n... eigenvalues\n");
    for (int iEig=0;iEig<numEigVals;++iEig) {
      this->com->fprintf(stdout, "%e\n", eVals[iEig]);
      if (eVals[iEig] > eVals[eig1]) {
        eig2=eig1;
        eig1=iEig;
      } else if (eVals[iEig] > eVals[eig2]) {
        eig2=iEig;
      }
    }
  
    this->com->fprintf(stdout, "\n2D visualization information: coord1 coord2 primaryCluster\n");
    for (int iSnap=0;iSnap<nSnaps;++iSnap) {
      this->com->fprintf(stdout, "%e %e %d\n", sqrt(eVals[eig1])*eVecs[nSnaps*eig1+iSnap], sqrt(eVals[eig2])*eVecs[nSnaps*eig2+iSnap], this->clusterIndex[iSnap]);
    }
  
    delete [] eVals;
    delete [] eVecs;
    delete [] symMat;
  }

#else
  this->com->fprintf(stderr, "*** Error: code must be compiled with ARPACK and DO_MODAL flag to compute MDS\n");
  return;
#endif

}


//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomDatabaseConstruction<dim>::testProbabilisticSVD(VecSet< DistSVec<double, dim> >*& snapshots, VecSet< DistSVec<double, dim> > &Utrue, std::vector<double>& singularValues, FullM &Vtrue, int podMethod, int randMatDimension, int nPowerIts, bool computeV)  {

  int nTotSnaps = snapshots->numVectors();
  int kStep = 100; 
  int numSteps = 1 + floor((double) nTotSnaps / ((double) kStep));

  // run scalapack
  this->com->fprintf(stdout,"\nRunning ScaLAPACK SVD\n");
  double scalapackTime = this->timer->getTime();
  scalapackSVD(snapshots, Utrue, singularValues, Vtrue, computeV);
  scalapackTime = this->timer->getTime() - scalapackTime;
  this->com->fprintf(stdout," ... elapsed time = %es\n",scalapackTime);

  // run R-SVD 
  this->com->fprintf(stdout,"\nRunning R-SVD\n");
  VecSet< DistSVec<double, dim> > URSVD(nTotSnaps, this->domain.getNodeDistInfo());
  FullM* VRSVD = new FullM(nTotSnaps);
  std::vector<double> singularValuesRSVD;
  double rTime = this->timer->getTime();
  rSVDWrapper(snapshots, URSVD, singularValuesRSVD, *VRSVD, false);
  rTime = this->timer->getTime() - rTime;
  this->com->fprintf(stdout," ... elapsed time = %es\n",rTime);

  // scalapack SVD errors
  std::vector<double> avgReconstructionErrorScalapack(numSteps,0.0);
  std::vector<double> maxReconstructionErrorScalapack(numSteps,0.0);

  // probabilistic SVD errors
  std::vector<double> sValMaxAbsError(numSteps,0.0);
  std::vector<double> sValMaxRelError(numSteps,0.0);
  std::vector<double> avgErrorU(numSteps,0.0);
  std::vector<double> maxErrorU(numSteps,0.0);
  std::vector<double> avgReconstructionErrorProb(numSteps,0.0);
  std::vector<double> maxReconstructionErrorProb(numSteps,0.0);
  std::vector<double> probabilisticTime(numSteps,0.0);

  // R-SVD errors
  std::vector<double> sValMaxAbsErrorRSVD(numSteps,0.0);
  std::vector<double> sValMaxRelErrorRSVD(numSteps,0.0);
  std::vector<double> avgErrorURSVD(numSteps,0.0);
  std::vector<double> maxErrorURSVD(numSteps,0.0);
  std::vector<double> avgReconstructionErrorRSVD(numSteps,0.0);
  std::vector<double> maxReconstructionErrorRSVD(numSteps,0.0);

  // first compute R-SVD errors and free up memory
  int k = 1;
  for (int step=0; step<numSteps; ++step) {
    k = min(k,nTotSnaps);
    
    // compare singular values with scalapack
    for (int i=0; i<k; ++i) {
      double dif = abs(singularValuesRSVD[i] -  singularValues[i]);
      sValMaxAbsErrorRSVD[step] = (dif>sValMaxAbsErrorRSVD[step]) ? dif : sValMaxAbsErrorRSVD[step];
      dif = (singularValues[i]>1e-14) ? dif / singularValues[i] : 1e6;
      sValMaxRelErrorRSVD[step] = (dif>sValMaxRelErrorRSVD[step]) ? dif : sValMaxRelErrorRSVD[step];
    }    

    // compare basis vectors with scalapack
    DistSVec<double,dim> errorVec( this->domain.getNodeDistInfo() );
    double errorNorm;
    for (int i=0; i<k; ++i) {
      errorVec = URSVD[i] -  Utrue[i];
      errorNorm = errorVec.norm();
      avgErrorURSVD[step] += errorNorm/((double)k);
      maxErrorURSVD[step] = (errorNorm>maxErrorURSVD[step]) ? errorNorm : maxErrorURSVD[step]; 
    }
 
    // compare reconstruction error to scalapack reconstruction error
    for (int iVec = 0; iVec < nTotSnaps; ++iVec) {
      errorVec = (*snapshots)[iVec];
      for (int jVec = 0; jVec < k; ++jVec)
        errorVec = errorVec - ((singularValues[jVec] * Vtrue[iVec][jVec]) * Utrue[jVec]);
      errorNorm = ((((*snapshots)[iVec]).norm()) > 1e-15) ? errorVec.norm()/(((*snapshots)[iVec]).norm()) : 0.0;
      avgReconstructionErrorScalapack[step] += errorNorm;
      if (errorNorm > maxReconstructionErrorScalapack[step])
        maxReconstructionErrorScalapack[step] = errorNorm;
    }
    avgReconstructionErrorScalapack[step] /= nTotSnaps;
   
    for (int iVec = 0; iVec < nTotSnaps; ++iVec) {
      errorVec = (*snapshots)[iVec];
      for (int jVec = 0; jVec < k; ++jVec)
        errorVec = errorVec - ((singularValuesRSVD[jVec] * (*VRSVD)[iVec][jVec]) * URSVD[jVec]);
      errorNorm = ((((*snapshots)[iVec]).norm()) > 1e-15) ? errorVec.norm()/(((*snapshots)[iVec]).norm()) : 0.0;
      avgReconstructionErrorRSVD[step] += errorNorm;
      if (errorNorm > maxReconstructionErrorRSVD[step])
        maxReconstructionErrorRSVD[step] = errorNorm;
    }
    avgReconstructionErrorRSVD[step] /= nTotSnaps;

    k += kStep;
  }

  URSVD.resize(0);
  delete VRSVD;
  VRSVD = NULL;

  // run probabilistic with variable k and compute errors
  this->com->fprintf(stdout,"\nRunning Probabilistic SVD\n");
  k = 1;
  for (int step=0; step<numSteps; ++step) {
    k = min(k,nTotSnaps);
   this->com->fprintf(stdout, " k=%d/%d\n", k, nTotSnaps);
    
    // run probabilistic SVD
    VecSet< DistSVec<double, dim> > Uprob(nTotSnaps, this->domain.getNodeDistInfo());
    FullM Vprob(nTotSnaps);
    std::vector<double> singularValuesProb;

    probabilisticTime[step] = this->timer->getTime();
    probabilisticSVDWrapper(snapshots, Uprob, singularValuesProb, Vprob, k, 0, false);
    probabilisticTime[step] = this->timer->getTime() - probabilisticTime[step];
    this->com->fprintf(stdout," elapsed time = %es\n\n",probabilisticTime[step]);

    // compare singular values with scalapack
    for (int i=0; i<k; ++i) {
      double dif = abs(singularValuesProb[i] -  singularValues[i]);
      sValMaxAbsError[step] = (dif>sValMaxAbsError[step]) ? dif : sValMaxAbsError[step];
      dif = (singularValues[i]>1e-14) ? dif / singularValues[i] : 1e6;
      sValMaxRelError[step] = (dif>sValMaxRelError[step]) ? dif : sValMaxRelError[step];
    }    

    // compare basis vectors with scalapack
    DistSVec<double,dim> errorVec( this->domain.getNodeDistInfo() );
    double errorNorm;
    for (int i=0; i<k; ++i) {
      errorVec = Uprob[i] -  Utrue[i];
      errorNorm = errorVec.norm();
      avgErrorU[step] += errorNorm/((double)k);
      maxErrorU[step] = (errorNorm>maxErrorU[step]) ? errorNorm : maxErrorU[step]; 
    } 
 
    // compare reconstruction error to scalapack reconstruction error
    for (int iVec = 0; iVec < nTotSnaps; ++iVec) {
      errorVec = (*snapshots)[iVec];
      for (int jVec = 0; jVec < k; ++jVec)
        errorVec = errorVec - ((singularValuesProb[jVec] * Vprob[iVec][jVec]) * Uprob[jVec]);
      errorNorm = ((((*snapshots)[iVec]).norm()) > 1e-15) ? errorVec.norm()/(((*snapshots)[iVec]).norm()) : 0.0;
      avgReconstructionErrorProb[step] += errorNorm;
      if (errorNorm > maxReconstructionErrorProb[step])
        maxReconstructionErrorProb[step] = errorNorm;
    }
    avgReconstructionErrorProb[step] /= nTotSnaps;

    k += kStep;
  }

  this->com->fprintf(stdout, " # dim scalapackTime avgReconstructionErrorScalapack maxReconstructionErrorScalapack rTime sValMaxAbsErrorRSVD sValMaxRelErrorRSVD avgErrorURSVD maxErrorURSVD avgReconstructionErrorRSVD maxReconstructionErrorRSVD probabilisticTime sValMaxAbsErrorProb sValMaxRelErrorProb avgErrorUProb maxErrorUProb  avgReconstructionErrorProb maxReconstructionErrorProb\n");
  k = 1;
  for (int step=0; step<numSteps; ++step) {
    this->com->fprintf(stdout, " %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", k, scalapackTime, avgReconstructionErrorScalapack[step], maxReconstructionErrorScalapack[step], rTime, sValMaxAbsErrorRSVD[step], sValMaxRelErrorRSVD[step], avgErrorURSVD[step], maxErrorURSVD[step], avgReconstructionErrorRSVD[step], maxReconstructionErrorRSVD[step], probabilisticTime[step], sValMaxAbsError[step], sValMaxRelError[step], avgErrorU[step], maxErrorU[step], avgReconstructionErrorProb[step], maxReconstructionErrorProb[step]);
    k += kStep;
  }
}

