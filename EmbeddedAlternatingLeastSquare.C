//
// Created by lei on 2/3/16.
//
#include <EmbeddedAlternatingLeastSquare.h>
#include <als_lapack.h>
#include <ParallelRomExtension.h>

template<int dim>
EmbeddedAlternatingLeastSquare<dim>::EmbeddedAlternatingLeastSquare(Communicator *_com, IoData &_ioData,
                                                                    Domain &_domain/*, DistGeoState *_geoSource*/) :
        NonlinearRom<dim>(_com, _ioData, _domain)/*, geoSource(_geoSource)*/ {
    this->numMasters = countMasters(this->domain.getNodeDistInfo());
    this->com->fprintf(stderr, "testing new I/O:\n");
    this->com->fprintf(stderr, "maxBasisSize = %d, relativeMinimumEnergy = %f, maxIteration = %d, leastSquareSolver = %d\n",
                        _ioData.romOffline.rob.embeddedALS.maxBasisSize,
                        _ioData.romOffline.rob.embeddedALS.relativeMinimumEnergy,
                        _ioData.romOffline.rob.embeddedALS.maxIteration,
                        _ioData.romOffline.rob.embeddedALS.leastSquareSolver);

    this->maxBasisSize = _ioData.romOffline.rob.embeddedALS.maxBasisSize;
    this->relativeMinimumEnergy = _ioData.romOffline.rob.embeddedALS.relativeMinimumEnergy;
    this->maxIteration = _ioData.romOffline.rob.embeddedALS.maxIteration;

    /*
    this->com->fprintf(stderr, "NormalizedSnaps = %d, SubtractClusterCenters = %d, SubtractNearestSnapsToCenters = %d, SubtractReferenceState = %d\n",
                       _ioData.romOffline.rob.embeddedALS.stateSnapshotsData.normalizeSnaps,
                       _ioData.romOffline.rob.embeddedALS.stateSnapshotsData.subtractCenters,
                       _ioData.romOffline.rob.embeddedALS.stateSnapshotsData.subtractNearestSnapsToCenters,
                       _ioData.romOffline.rob.embeddedALS.stateSnapshotsData.subtractRefState);

    this->normalizeSnaps = _ioData.romOffline.rob.embeddedALS.stateSnapshotsData.normalizeSnaps ? true : false;
    this->subtractCenters = _ioData.romOffline.rob.embeddedALS.stateSnapshotsData.subtractCenters ? true : false;
    this->subtractReferenceState = _ioData.romOffline.rob.embeddedALS.stateSnapshotsData.subtractRefState ? true : false;
     */
}

template<int dim>
EmbeddedAlternatingLeastSquare<dim>::~EmbeddedAlternatingLeastSquare() {
    //delete [] mask;
}

template<int dim>
void EmbeddedAlternatingLeastSquare<dim>::constructDatabase() {
    // read in all snapshots
    this->numSnapshots = this->readSnapshotsFilesHelper("mask", false);
    this->nSnapShotFiles = this->readSnapshotsFilesHelper("state", false);

    // dummy clustering step
    this->com->printf(10, "dummy clustering stage\n");
    this->dummyClustering();
    this->outputClusteredSnapshots("state");
    // temporary hotfix
    this->snapsInCluster = new int[1];
    this->snapsInCluster[0] = this->numSnapshots;

    // run ALS
    this->com->printf(10, "load cluster center\n");
    if(this->ioData->romOffline.rob.state.snapshots.subtractCenters)
        this->readClusterCenters("centers"); //todo: program exits here
    this->com->printf(10, "load nearest snaps\n");
    if(this->ioData->romOffline.rob.state.snapshots.subtractNearestSnapsToCenters)
        this->readNearestSnapsToCenters();
    this->com->printf(10, "reading clustered snapshots, with preprocessing\n");
    this->readClusteredSnapshots(0, true, "state", 0, -1);
    this->com->printf(10, "reduced basis construction\n");
    this->ReducedOrderBasisConstruction();
}

template<int dim>
int EmbeddedAlternatingLeastSquare<dim>::readSnapshotsFilesHelper(char *keyword, bool preprocess) {
    if (strcmp(keyword, "mask") == 0)
        return readStateMaskFile();
    else
        return readSnapshotFiles(keyword, preprocess);
}

/**
 * Assuming that
 * 1. all snapshots in one cluster;
 * 2. readSnapshotFiles("state") has been called or this->snap has been set;
 * 3. readSnapshotFiles("mask") has been called or this->numSnapshots has been set.
 * This function sets the following variables (see ROM Database data):
 * int* this->clusterIndex = all ones;
 * int** this->clusterSnapshotMap[0] = all ones;
 * int** this->clusterNeighbours[0] = all ones;
 * this->clusterCenters[0] = mean of all snaps;
 * this->nearestSnapsToCenters[0] = closest snap to mean;
 */
template<int dim>
void EmbeddedAlternatingLeastSquare<dim>::dummyClustering(){
    //only works for one cluster
    assert(this->nClusters == 1);
    assert(this->nClusters == 1);
    // set all ROM Database data
    this->clusterIndex = new int[this->numSnapshots];
    for(int i = 0; i < this->numSnapshots; i++){
        this->clusterIndex[i] = 0;
    }
    // cluster with overlap
    this->snapsInCluster = new int[1];
    this->snapsInCluster[0] = this->numSnapshots;
    this->clusterSnapshotMap = new int*[1];
    this->clusterSnapshotMap[0] = new int[this->numSnapshots];
    for(int i = 0; i < this->numSnapshots; i++){
        this->clusterSnapshotMap[0][i] = i;
    }
    // overlap data
    this->clusterNeighborsCount = new int[1];
    this->clusterNeighborsCount[0] = 0;
    this->clusterNeighbors = new int*[1];
    this->clusterNeighbors[0] = new int[1];
    this->clusterNeighbors[0][0] = 0;
    // compute mean of all snapshots
    this->clusterCenters = new VecSet<DistSVec<double, dim> >(1, this->domain.getNodeDistInfo());
    for(int i = 0; i < this->numSnapshots; i++){
        (*this->clusterCenters)[0] += (*this->snap)[i];
    }
    (*this->clusterCenters)[0] /= this->numSnapshots;
    // compute nearestSnapsToCenters;
    this->nearestSnapsToCenters = new VecSet<DistSVec<double, dim> >(1, this->domain.getNodeDistInfo());
    double distance = 1e8;
    DistSVec<double, dim> tempDistSVec = DistSVec<double, dim>(this->domain.getNodeDistInfo());
    int minIndex = -1;
    for(int i = 0; i < this->numSnapshots; i++){
        tempDistSVec = (*this->snap)[i] - (*this->clusterCenters)[0];
        if (tempDistSVec.norm() > distance + 1e-8)
            continue;
        distance = tempDistSVec.norm();
        minIndex = i;
    }
    (*this->nearestSnapsToCenters)[0] = (*this->snap)[minIndex];
}


//TODO: all snapshots stroed in one file. need to separate them
template<int dim>
int EmbeddedAlternatingLeastSquare<dim>::readStateMaskFile() {
    // get the name of state mask snapshot summary
    char *filename = new char[strlen(this->ioData->input.prefix) + strlen(this->ioData->input.stateMaskSnapFile) + 1];
    sprintf(filename, "%s%s", this->ioData->input.prefix, this->ioData->input.stateMaskSnapFile);

    // open summary file
    FILE *in = fopen(filename, "r");
    if(!in){
        this->com->fprintf(stderr, "***Error: NO snapshots FILES in %s\n", filename);
        exit(-1);
    }

    // read the first line from summary
    int nData, _n;
    _n = fscanf(in, "%d", &nData);
    this->com->fprintf(stdout, "Reading mask snapshots from %d files \n", nData);

    // scan all lines from summary
    char** maskFile = new char*[nData];
    for(int iData = 0; iData < nData; iData++)
        maskFile[iData] = new char[500];
    int* startIndex = new int[nData];
    int* endIndex = new int[nData];
    int* sampleFrequency = new int[nData];

    for(int iData = 0; iData < nData; iData++){
        char snapshotFilename[500];
        int start, end, freq,
        _n = fscanf(in, "%s %d %d %d", snapshotFilename, &start, &end, &freq);
        strcpy(maskFile[iData], snapshotFilename);
        startIndex[iData] = start < 1 ? 0 : start - 1;
        endIndex[iData] = end < 0 ? 0 : end;
        sampleFrequency[iData] = freq < 1 ? 1: freq;
        this->com->fprintf(stdout, " ... Reading snapshots from %s \n", maskFile[iData]);
    }

    // compute total number of snapshots from all files
    int numSnapshots = 0;
    int* numSnaps = new int[nData];
    for(int iData = 0; iData < nData; iData++){
        int dummyStep = 0;
        double dummyTag = 0.0;
        bool status = this->domain.template readTagFromFile<double, dim>(maskFile[iData], dummyStep, &dummyTag, &(numSnaps[iData]));
        if(!status){
            this->com->fprintf(stdout, "*** ERROR: could not read mask from %s \n", maskFile[iData]);
            exit(-1);
        }
        if(endIndex[iData] == 0 || endIndex[iData] > numSnaps[iData])
            endIndex[iData] = numSnaps[iData];
        for(int i = startIndex[iData]; i < endIndex[iData]; i++){
            if(i % sampleFrequency[iData] != 0) continue;
            numSnapshots++;
        }
    }
    this->com->fprintf(stderr, "%d state mask snapshots found\n", numSnapshots);
    this->numSnapshots = numSnapshots;
    //assert(this->numSnapshots == numSnapshots);

    // load the all snapshots into this->mask
    this->mask = new VecSet< DistSVec<char, dim> >(numSnapshots, this->domain.getNodeDistInfo());
    DistSVec<char, dim> *stateMask = new DistSVec<char, dim>(this->domain.getNodeDistInfo());
    int currentIndex = 0;
    for(int iData = 0; iData < nData; iData++){
        for(int i = startIndex[iData]; i < endIndex[iData]; i++){
            if(i % sampleFrequency[iData] != 0) continue;
            double tag;
            int status = this->domain.readVectorFromFile(maskFile[iData], i, &tag, *stateMask);
            (*this->mask)[currentIndex] = *stateMask;
            currentIndex++;
        }
    }
    assert(currentIndex == numSnapshots);

    //clean up all variables
    delete [] filename;
    delete [] maskFile;
    delete [] numSnaps;
    delete [] startIndex;
    delete [] endIndex;
    delete [] sampleFrequency;
    return numSnapshots;
}


/**
 * suppose M is the number of master nodes
 */
template<int dim>
void EmbeddedAlternatingLeastSquare<dim>::freeSlaves(double *&mem, const VecSet<DistSVec<double, dim> > &X, const int M,
                                                     const int N) {
    double* ptr = new double[dim * M * N];
    assert(ptr != NULL);
    for (int j = 0; j < N; j++) {
        int i = 0;
//#pragma omp parallel for
        for (int iSub = 0; iSub <X[j].numLocSub() ; iSub++) {
            bool *localMasterFlag = X[j].getMasterFlag(iSub);
            double (*temp)[dim] = X[j].subData(iSub);
            for (int iNode = 0; iNode < X[j].subSize(iSub); iNode++) {
                if (localMasterFlag[iNode]) {
                    for(int k = 0; k < dim; k++){
                        ptr[dim * M *j + i * dim + k] = temp[iNode][k];
                    }
                    //memcpy(&ptr[dim * M * j + i * dim], &X[j][iSub][iNode], dim * sizeof(double));
                    i++;
                }
            }
        }
        //assert(i * dim == M); // should be true
        if( i != M)
            this->com->fprintf(stderr, "double array, dim is %d, %d != %d\n", dim, i, M);
    }
    mem = ptr;
}


template<int dim>
void EmbeddedAlternatingLeastSquare<dim>::freeSlaves(char *&mem, const VecSet<DistSVec<char, dim> > &X, const int M,
                                                     const int N) {
    char* ptr = new char[dim * M * N];
    assert(ptr != NULL);
    for (int j = 0; j < N; j++) {
        int i = 0;
//#pragma omp parallel for
        for (int iSub = 0; iSub <X[j].numLocSub() ; iSub++) {
            bool *localMasterFlag = X[j].getMasterFlag(iSub);
            char (*temp)[dim] = X[j].subData(iSub);
            for (int iNode = 0; iNode < X[j].subSize(iSub); iNode++) {
                if (localMasterFlag[iNode]) {
                    for(int k = 0; k < dim; k++){
                        ptr[dim * M *j + i * dim + k] = temp[iNode][k];
                    }
                    //memcpy(&ptr[dim * M * j + i * dim], &X[j][iSub][iNode], dim * sizeof(double));
                    i++;
                }
            }
        }
        //assert(i  == M); // should be true
        if( i != M)
            this->com->fprintf(stderr, "char array, dim is %d, %d != %d\n", dim, i, M);
    }
    mem = ptr;
}

template<int dim>
void EmbeddedAlternatingLeastSquare<dim>::summonSlaves(double *&mem, VecSet<DistSVec<double, dim> > &X, const int M,
                                                       const int N) {
    SubDomain **subDomain = this->domain.getSubDomain();
    CommPattern<double> *vecPattern = this->domain.getVecPat();
    DistInfo& distInfo = this->domain.getNodeDistInfo();

    double *ptr = mem;
    for (int j = 0; j < N; j++) {
        int i = 0;
//#pragma omp parallel for
        for (int iSub = 0; iSub < X[j].numLocSub(); iSub++) {
            bool *localMasterFlag = X[j].getMasterFlag(iSub);
            double (*temp)[dim] = X[j].subData(iSub);
            for (int iNode = 0; iNode < X[j].subSize(iSub); iNode++) {
                if (localMasterFlag[iNode]) {
                    for(int k = 0; k < dim; k++){
                        temp[iNode][k] = ptr[dim * M * j + i * dim + k];
                    }
                    //memcpy(&X[j][iSub][iNode], &ptr[dim * M * j + i * dim], dim * sizeof(double));
                    i++;
                }
            }
            subDomain[iSub]->sndData(*vecPattern, X[j].subData(iSub));
        }
        //assert(i * dim == M); // should be true
        if(i != M) this->com->fprintf(stderr, "dim is %d, i(%d) != M(%d)", dim, i, M);
        vecPattern->exchange();

        for (int iSub = 0; iSub < X[j].numLocSub(); iSub++)
            subDomain[iSub]->addRcvData(*vecPattern, X[j].subData(iSub));
    }
}


template<int dim>
void EmbeddedAlternatingLeastSquare<dim>::summonZombies(double *&mem, VecSet<DistSVec<double, dim> > &X, const int M,
                                                       const int N) {
    SubDomain **subDomain = this->domain.getSubDomain();
    CommPattern<double> *vecPattern = this->domain.getVecPat();
    DistInfo& distInfo = this->domain.getNodeDistInfo();

    double *ptr = mem;
    for (int j = 0; j < N; j++) {
        int i = 0;
//#pragma omp parallel for
        for (int iSub = 0; iSub < X[j].numLocSub(); iSub++) {
            double (*temp)[dim] = X[j].subData(iSub);
            bool *localMasterFlag = X[j].getMasterFlag(iSub);
            for (int iNode = 0; iNode < X[j].subSize(iSub); iNode++) {
                if (localMasterFlag[iNode]) {
                    for(int k = 0; k < dim; k++){
                        temp[iNode][k] = ptr[dim * M * j + i * dim + k];
                    }
                   // memcpy(&X[j][iSub][iNode], &ptr[dim * M * j + i * dim], dim * sizeof(double));
                    i++;
                }
            }
        }
        assert(i * dim == M); // should be true
    }
}

//todo: save singular values for multiple runs
template<int dim>
int EmbeddedAlternatingLeastSquare<dim>::initialization(VecSet<DistSVec<double, dim> > &basisInit) {
    const int n = min(this->maxBasisSize, this->numSnapshots);
    VecSet<DistSVec<double, dim> >* U = new VecSet<DistSVec<double, dim> >(n, this->domain.getNodeDistInfo());
    int ncol = this->numSnapshots;
    FullM V(n, ncol);
    double* singularValues = new double[this->maxBasisSize];
    this->com->fprintf(stderr, " ... calling parallelROM SVD, n is %d\n", n);
    ParallelRom<dim> parallelRom(this->domain, this->com, this->domain.getNodeDistInfo());
    parallelRom.parallelSVD(*(this->snap), *U, singularValues, V, n, true);

    double singularValuesSum = 0;
    double remainingSumEstimate = 0;
    this->reducedDimension = -1;
    for(int i = 0; i < n; i++){
        double s = singularValues[i];
        remainingSumEstimate = s * s * (ncol - i);
        double percentage = singularValuesSum / (remainingSumEstimate + singularValuesSum);
        bool stopped = singularValuesSum > this->relativeMinimumEnergy * (remainingSumEstimate + singularValuesSum);
        this->com->fprintf(stderr, "s = %f, current sum = %f, remaining sum ~= %f, percentage = %f\n", s, singularValuesSum, remainingSumEstimate, percentage);
        if (stopped) {
            this->reducedDimension = i;
            break;
        }
        singularValuesSum += s * s;
    }
    if(this->reducedDimension < 0) this->reducedDimension = n;
    this->com->fprintf(stderr, "... reduced dimension is %d, initializing Basis accordingly\n", this->reducedDimension);
    basisInit.resize(this->reducedDimension);
    for(int i = 0; i < this->reducedDimension; i++){
        basisInit[i] = (*U)[i];
    }
    delete U;
    return this->reducedDimension;
}

/**
 * When no reduced dimension is provided, use relative minimum energy
 * to determine reduced dimension.
 */
template<int dim>
void EmbeddedAlternatingLeastSquare<dim>::ReducedOrderBasisConstruction() {
    int n = min(this->maxBasisSize, (*(this->snap)).numVectors());
    VecSet<DistSVec<double, dim> > *basisInitTemp = new VecSet<DistSVec<double, dim> >(n, this->domain.getNodeDistInfo());
    VecSet<DistSVec<double, dim> > basisInit = *basisInitTemp;
    int k = initialization(basisInit);
    delete basisInitTemp;
    ReducedOrderBasisConstruction(k);
}

template<int dim>
void EmbeddedAlternatingLeastSquare<dim>::ReducedOrderBasisConstructionTesting(int k) {
    VecSet< DistSVec<double, dim> >* basisInit = new VecSet< DistSVec<double, dim> >(k, this->domain.getNodeDistInfo());
    VecSet< DistSVec<double, dim> > UInit = *basisInit;
    ParallelRomExtension<dim> parallelRomExtension(this->domain, this->com, this->domain.getNodeDistInfo());
    parallelRomExtension.parallelALS(*(this->snap), *(this->mask), UInit, this->maxIteration);
    // do QR before writing to disk
    std::vector<std::vector<double> > RT(k, std::vector<double>(k));
    qr(&UInit, &RT, true);
    this->com->fprintf(stderr, "... using Kyle's QR with no pivoting\n");
    outputBasis(UInit);
}

//TODO: put this as part of NonlinearRomDatabaseConstruction<dim>::SVD()
template<int dim>
void EmbeddedAlternatingLeastSquare<dim>::ReducedOrderBasisConstruction(int _dim) {
    //read snapshots and mask
    //assert(numSnapshots == numStateSnapshots);
    // compute nrow, ncol, dim
    //todo: determine dim from initial SVD
    int nrow = 0;
    for(int i = 0; i < this->numMasters.size(); i++){
        nrow += this->numMasters[i];
    }
    int ncol = (*(this->snap)).numVectors();
    int k = _dim; //TODO: change it
    this->reducedDimension = _dim; //TODO: change it
    this->com->fprintf(stderr, "... (M, N) is [%d, %d], reduced dimension is %d\n", nrow, ncol, this->reducedDimension);
    // set up X, M
    VecSet< DistSVec<double, dim> > Snap = *(this->snap);
    VecSet< DistSVec<char, dim> > Mask = *(this->mask);
    double *X = NULL;
    char *M = NULL;
    freeSlaves(X, Snap, nrow, ncol); // freeSlaves allocates the memory for X
    freeSlaves(M, Mask, nrow, ncol); // freeSlaves allocates the memory for M
    // set up Uinit
    this->com->fprintf(stderr, "... initializing U\n");
    VecSet< DistSVec<double, dim> >* basisInit = new VecSet< DistSVec<double, dim> >(k, this->domain.getNodeDistInfo());
    VecSet< DistSVec<double, dim> > basis = *basisInit;
    ParallelRom<dim> parallelRom(this->domain, this->com, this->domain.getNodeDistInfo());
    double *singularValues = new double[k]; //todo: change definition.
    FullM VInitDummy(this->reducedDimension, ncol);
    this->com->fprintf(stderr, "... calling parallelRom.parallelSVD()\n");
    parallelRom.parallelSVD(Snap, basis, singularValues, VInitDummy, this->reducedDimension, true);
    this->com->fprintf(stderr, "... U and V initialized, V dimension is [%d, %d]\n", VInitDummy.numRow(), VInitDummy.numCol());
    double *U = NULL;
    double *UT = new double[k * nrow * dim];
    this->com->fprintf(stderr, "... allocating space for U, k * nrow * dim = %d\n", k * nrow * dim);
    freeSlaves(U, basis, nrow, k); // freeSlaves allocate the memory for U
    transpose(U, UT, nrow * dim, k); // does not allocate memory
    this->com->fprintf(stderr, "... U being initialized\n");
    this->com->barrier();
    // launch ALS external library
    AlternatingLeastSquare ALS((double *)X, (unsigned char *) M, (double *) UT, nrow * dim, ncol, this->reducedDimension, this->com->getMPIComm(), this->com->cpuNum());
    ALS.run(this->maxIteration);
    // write basis to file using parent class methods
    this->com->fprintf(stderr, "... ALS finished %d iterations\n", this->maxIteration);
    this->com->barrier();
    transpose(UT, U, k, nrow * dim);
    this->com->fprintf(stderr, "... transpose done\n");
    this->com->barrier();
    summonSlaves(U, basis, nrow, k);
    //todo: do QR before writing to disk
    std::vector<std::vector<double> > RT(k, std::vector<double>(k));
    qr(&basis, &RT, true);
    this->com->fprintf(stderr, "... using Kyle's QR with no pivoting\n");
    outputBasis(basis);
    //this->com->barrier();
    this->com->fprintf(stderr, "... cleaning up memories\n");
    // clean up all allocated memory
    if(X) delete[] X;
    if(M) delete[] M;
    if(U) delete[] U;
    if(UT) delete[] UT;
    if(singularValues) delete[] singularValues;
    //this->com->barrier();
    this->com->fprintf(stdout, "... all stuff cleaned\n");
}


template<int dim>
void EmbeddedAlternatingLeastSquare<dim>::outputBasis(const VecSet<DistSVec<double, dim> >& U) {
    char *basisPath = NULL;
    determinePath(this->stateBasisName, 0 /* only one cluster */, basisPath);
    int reducedDimension = U.numVectors();
    this->com->fprintf(stdout, "\n writing %d basis from cluster %d to disk\n", reducedDimension, 0);
    for(int i = 0; i < reducedDimension; i++)
        this->domain.writeVectorToFile(basisPath, i, double(i), U[i]);
}


template<int dim>
void EmbeddedAlternatingLeastSquare<dim>::readBasisFiles(VecSet<DistSVec<double, dim> >& U) {
    char *basisPath = NULL;
    determinePath(this->stateBasisName, 0 /* only one cluster */, basisPath);
    int reducedDimension;
    int dummyStep = 0;
    double dummyTag = 0.0;
    bool status = this->domain.template readTagFromFile<double, dim>(basisPath, dummyStep, &dummyTag, &reducedDimension);
    if(!status){
        this->com->fprintf(stdout, "*** ERROR: could not read basis from %s \n", basisPath);
        exit(-1);
    }
    U.resize(reducedDimension);
    this->com->fprintf(stdout, "\n reading %d basis from cluster %d to disk \n", reducedDimension, 0);
    for(int i = 0; i < reducedDimension; i++) {
        int status = this->domain.readVectorFromFile(basisPath, i, &dummyTag, U[i]);
    }
}


template<int dim>
void EmbeddedAlternatingLeastSquare<dim>::readReferenceStateFiles(DistSVec<double, dim>& U) {
    char *refStatePath = NULL;
    determinePath(this->refStateName, 0 /* only one cluster */, refStatePath);
    double dummyTag = 0.0;
    int status = this->domain.readVectorFromFile(refStatePath, 0, &dummyTag, U);
    if(!status){
        this->com->fprintf(stdout, "*** ERROR: could not read reference state from %s \n", refStatePath);
        exit(-1);
    }
}

/**
 *
 */
template<int dim>
void EmbeddedAlternatingLeastSquare<dim>::testingSnapshotIO() {
    // compute M and N
    int M = 0;
    for(int i = 0; i < this->numMasters.size(); i++){
        M += this->numMasters[i];
    }
    int N = this->numSnapshots;
    this->com->fprintf(stdout, "(M, N) is [%d, %d]\n", M, N);
    // set up X, Y, mem
    VecSet< DistSVec<double, dim> > X = *(this->snap);
    VecSet< DistSVec<double, dim> > *YY = new VecSet< DistSVec<double, dim> >(N, this->domain.getNodeDistInfo());
    VecSet< DistSVec<double, dim> > Y = *YY;
    VecSet< DistSVec<double, dim> > *ZZ = new VecSet< DistSVec<double, dim> >(N, this->domain.getNodeDistInfo());
    VecSet< DistSVec<double, dim> > Z = *ZZ;
    VecSet< DistSVec<char, dim> > W = *(this->mask);
    double* mem = NULL;
    // testing isEqualMatrices() implementation
    //this->com->fprintf(stdout, "testing equal matrix\n");
    //bool res = isEqualMatrices(X, X);
    bool res;
    this->com->fprintf(stdout, "dim is %d\n", dim);
    this->com->fprintf(stdout, "testing mask matrix\n");
    res = isEqualMaskMatrices(W, W);
    this->com->fprintf(stdout, "testing snapshot matrix\n");
    res = isEqualStateMatrices(X, X);
    //this->com->fprintf(stdout, "result of X == X is %d\n", res);
    freeSlaves(mem, X, M, N);
    summonSlaves(mem, Y, M, N);
    summonZombies(mem, Z, M, N);
    this->com->fprintf(stdout, "testing freeslave() and summonslave()\n");
    res = isEqualStateMatrices(X, Z);
    this->com->fprintf(stdout, "result of X == Z is %d\n", res);
    res = isEqualStateMatrices(X, Y);
    this->com->fprintf(stdout, "result of X == Y is %d\n", res);
}

template<int dim>
void EmbeddedAlternatingLeastSquare<dim>::testingALS(){
    AlternatingLeastSquare als(25000, 2000, 20, this->com->getMPIComm(), this->com->cpuNum());
    als.run(8);
    MPI_Barrier(this->com->getMPIComm());
    this->com->fprintf(stdout, "testing ALS completed\n");
}

template<int dim>
void EmbeddedAlternatingLeastSquare<dim>::testingInitialization() {
    int n = min(this->numSnapshots, this->maxBasisSize);
    VecSet<DistSVec<double, dim> > *basisInitTemp = new VecSet<DistSVec<double, dim> >(n, this->domain.getNodeDistInfo());
    VecSet<DistSVec<double, dim> > basisInit = *basisInitTemp;
    int k = initialization(basisInit);
    this->com->fprintf(stderr, "reduced dimension is %d\n", k);
}

/**
 * count the number of master nodes in each subdomain.
 * results stored in this->numMasters[i].
 */
template<int dim>
std::vector<int> EmbeddedAlternatingLeastSquare<dim>::countMasters(DistInfo &distinfo){
    vector<int> answer(distinfo.numLocSub, 0);
//#pragma omp parallel for
    for(int iSub = 0; iSub < distinfo.numLocSub; iSub++){
        bool *localMasterFlag = distinfo.getMasterFlag(iSub);
        for(int iNode = 0; iNode < distinfo.subLen[iSub]; iNode++){
            answer[iSub] += localMasterFlag[iNode] ? 1 : 0;
        }
    }
    return answer;
}

template<int dim>
bool EmbeddedAlternatingLeastSquare<dim>::isEqualStateMatrices(const VecSet<DistSVec<double, dim> > &X,
                                                          const VecSet<DistSVec<double, dim> > &Y){
    int N1 = X.numVectors();
    int N2 = Y.numVectors();
    int M1 = X[0].len;
    int M2 = Y[0].len;
    //dimension mismatch
    this->com->fprintf(stdout, "X = [%d, %d], Y = [%d, %d]\n", M1, N1, M2, N2);
    if(M1 != M2 || N1 != N2) return false;
    //distInfo mismatch
    this->com->fprintf(stdout, "numLocSub: %d %d\n", X[0].numLocSub(), Y[0].numLocSub());
    if(X[0].numLocSub() != Y[0].numLocSub()) return false;
    for(int i = 0; i < X[0].numLocSub(); i++){
        this->com->fprintf(stdout, "subdomain %d: %d, %d\n", i, X[0].subSize(i), Y[0].subSize(i));
        if(X[0].subSize(i) != Y[0].subSize(i)) return false;
    }
    DistInfo& distinfo = this->domain.getNodeDistInfo();
//#pragma omp parallel for
    for(int j = 0; j < N1; j++) {
        for(int iSub = 0; iSub < distinfo.numLocSub; iSub++){
            bool* XLocalMasterFlag = X[j].getMasterFlag(iSub);
            bool* YLocalMasterFlag = Y[j].getMasterFlag(iSub);
            double (*xtemp)[dim] = X[j].subData(iSub);
            double (*ytemp)[dim] = Y[j].subData(iSub);
            for(int iNode = 0; iNode < distinfo.subLen[iSub]; iNode++){
                /*this->com->fprintf(stdout, "[%d, %d, %d]: master(X, Y) = %d, %d \n",
                                   j, iSub, iNode,
                                   XLocalMasterFlag[iNode], YLocalMasterFlag[iNode]);*/
                for(int k = 0; k < dim; k++){
                   // this->com->fprintf(stdout, "[%d]: %f, %f", k, xtemp[iNode][k], ytemp[iNode][k]);
                    if(xtemp[iNode][k] != ytemp[iNode][k]){
                        this->com->fprintf(stdout, "[%d]: %f, %f\n", k, xtemp[iNode][k], ytemp[iNode][k]);
                        return false;
                    }
                }
               // this->com->fprintf(stdout, "\n");
            }
        }
    }
    this->com->fprintf(stdout, "all test passed\n");
    return true;
}

template<int dim>
bool EmbeddedAlternatingLeastSquare<dim>::isEqualMaskMatrices(const VecSet<DistSVec<char, dim> > &X,
                                                          const VecSet<DistSVec<char, dim> > &Y){
    int N1 = X.numVectors();
    int N2 = Y.numVectors();
    int M1 = X[0].len;
    int M2 = Y[0].len;
    //dimension mismatch
    this->com->fprintf(stdout, "X = [%d, %d], Y = [%d, %d]\n", M1, N1, M2, N2);
    if(M1 != M2 || N1 != N2) return false;
    //distInfo mismatch
    this->com->fprintf(stdout, "numLocSub: %d %d\n", X[0].numLocSub(), Y[0].numLocSub());
    if(X[0].numLocSub() != Y[0].numLocSub()) return false;
    for(int i = 0; i < X[0].numLocSub(); i++){
        this->com->fprintf(stdout, "subdomain %d: %d, %d\n", i, X[0].subSize(i), Y[0].subSize(i));
        if(X[0].subSize(i) != Y[0].subSize(i)) return false;
    }
    DistInfo& distinfo = this->domain.getNodeDistInfo();
//#pragma omp parallel for
    for(int j = 0; j < N1; j++) {
        for(int iSub = 0; iSub < distinfo.numLocSub; iSub++){
            bool* XLocalMasterFlag = X[j].getMasterFlag(iSub);
            bool* YLocalMasterFlag = Y[j].getMasterFlag(iSub);
            char (*xtemp)[dim] = X[j].subData(iSub);
            char (*ytemp)[dim] = Y[j].subData(iSub);
            for(int iNode = 0; iNode < distinfo.subLen[iSub]; iNode++){
                /*this->com->fprintf(stdout, "[%d, %d, %d]: master(X, Y) = %d, %d \n",
                                   j, iSub, iNode,
                                   XLocalMasterFlag[iNode], YLocalMasterFlag[iNode]);
                                   */
                for(int k = 0; k < dim; k++){
                    //this->com->fprintf(stdout, "[%d]: %c, %c", k, xtemp[iNode][k], ytemp[iNode][k]);
                    if(xtemp[iNode][k] != ytemp[iNode][k]){
                        this->com->fprintf(stdout, "[%d]: %d, %d", k, xtemp[iNode][k] == 0 ? 0 : 1 , ytemp[iNode][k] == 0 ? 0 : 1);
                        return false;
                    }
                }
                //this->com->fprintf(stdout, "\n");
            }
        }
    }
    this->com->fprintf(stdout, "all test passed\n");
    return true;
}


template<int dim>
void EmbeddedAlternatingLeastSquare<dim>::transpose(double* &buff1, double* &buff2, int nrow, int ncol){
    double *ptr1 = buff1;
    double *ptr2 = buff2;
    assert(ptr2 != NULL);
    for(int j = 0; j < ncol; j++){
        for(int i = 0; i < nrow; i++){
            ptr2[i * ncol + j] = ptr1[j * nrow + i];
        }
    }
}
/*
template<int dim>
GenFullM<double> EmbeddedAlternatingLeastSquare<dim>::getRow(VecSet<DistSVec<double, dim> >& X, int numCol, int subVecIndex, int nodeIndex){
    GenFullM<double> result(dim, numCol);
    for(int i = 0; i < numCol; i++){
        result[i] = X[i].subData(subVecIndex)[nodeIndex];
    }
    return result;
}

template<int dim>
void EmbeddedAlternatingLeastSquare<dim>::setRow(VecSet<DistSVec<double, dim> &X, GenFullM<double> &row, int numCol, int subVecIndex, int nodeIndex){
    for(int i = 0; i < numCol; ++i){
        X[i].subData(subVecIndex)[nodeIndex] = row[i]; // not sure if pointer copy or deep copy
    }
}

template<int dim>
GenFullM<double> EmbeddedAlternatingLeastSquare<dim>::kroneckerProduct(double *v, double *u, int m, int n){
    GenFullM<double> result(m, n);
    for(int i = 0; i < m; ++i){
        for(int j = 0; j < n; ++j){
            result[i][j] = v[i] * u[j];
        }
    }
    return result;
}
*/

/*
template <int dim>
int EmbeddedAlternatingLeastSquare<dim>::writeReducedOrderBasisToFile() {

}

template  <int dim>
int EmbeddedAlternatingLeastSquare<dim>::AlternatingLeastSquareMethodI(
        VecSet<DistSVec<double, dim> > &X,
        VecSet<DistVec<bool> > &M,
        VecSet<DistSVec<double, dim> > &U,
        GenFullM &V,
        int rank,
        int maxIteration,
        double lambda ) {
    int numCol = X.numVectors();
    int numRow = X[0].size();
    int it = 0;
    double error = 1e10;
    while(it < maxIteration && error > lambda){
#pragma omp parallel for
        for(int j = 0; j < numCol; ++j){ // loop over columns
            GenFullM<double> M(numCol, numCol);
            GenFullM<double> v(numCol, 1);
            M = 0;
            v = 0;
            DistInfo distInfo = mask[0].info();
            for(int iSub = 0; isub < numLocSub; ++iSub) { // loop over subdomains
                bool *localMask = mask[j].subData(iSub);
                bool *localMasterFlag = X[j].getMasterFlag(iSub);
                for(int iNode = 0; iNode < distInfo.subLen[iSub]; ++iNode) { // loop over rows
                    if(!localMasterFlag && localMasterFlag[iNode] && !localMask[iNode]){ // loop over unmasked and master nodes
                        GenFullM<double> A = getRow(U, numCol, iSub, iNode);
                        M += A^A; // or M.add(A^A, 0, 0);
                        v += A^ X[j][iSub][iNode];
                    }
                }
            }
            M.Factor(1e-9);
            M.ReSolve(v.data());
            V[j] = v.data();
        }
        // iterate over rows of X to update rows of U
#pragma omp parallel for
        for(int iSub = 0; isub < numLocSub; ++iSub){ //loop over subdomains
            for(int iNode = 0; iNode < distInfo.subLen[iSub]; ++iNode){ //loop over rows
                GenFullM<double> M(numCol, numCol);
                GenFullM<double> v(numCol, 3);
                M = 0;
                v = 0;
                for(int j = 0; j < numCol; j++){ //loop over cols
                    if(!X[j].getMasterFLag(iSub) && X[j].getMasterFlag(iSub)[iNode] && !mask[j].subData(iSub)[iNode]) { // loop over unmasked and master nodes}
                        M += kroneckerProduct(v, v, numCol, numCol);
                        X[j].subData(iSub)[iNode]
                        v += kronecherProduct(v, X[j].subData(iSub)[iNode], numCol, dim);

                    }
                }
            }
            // equivalent matlab code
            x = X(i, M(i, :));
            v = V(:, M(i, :));
            U(i, :) = ((v * v')\(v * x'))';


}
it++;
}

return 1;
}
*/

/* Lei Lei, 21 March 2016: does not work because NonlinearRom.h does not use this methods.
#define INSTANTIATION_HELPER(T, d) \
    template bool EmbeddedAlternatingLeastSquare<d>::isEqualMatrices(const VecSet<DistSVec<T, d> > &, const VecSet<DistSVec<T, d> > &); \
    template std::vector<int> EmbeddedAlternatingLeastSquare<d>::countMasters(DistInfo &); \
    template void EmbeddedAlternatingLeastSquare<d>::testingSnapshotIO(); \
    template void EmbeddedAlternatingLeastSquare<d>::ReducedOrderBasisConstruction(); \
    template void EmbeddedAlternatingLeastSquare<d>::summonSlaves(const void *&, VecSet<DistSVec<double, dim> > &, const int M, const int N); \
    template void EmbeddedAlternatingLeastSquare<d>::freeSlaves(void *&, const VecSet<DistSVec<double, dim> > &X, const int M, const int N); \
    template int EmbeddedAlternatingLeastSquare<d>::readStateMaskFile(); \
    template int EmbeddedAlternatingLeastSquare<d>::readSnapshotsFilesHelper(char *, bool); \
    template EmbeddedAlternatingLeastSquare<d>::EmbeddedAlternatingLeastSquare(Communicator *, IoData &, Domain &); \
    template EmbeddedAlternatingLeastSquare<d>::~EmbeddedAlternatingLeastSquare();

INSTANTIATION_HELPER(double, 1);
INSTANTIATION_HELPER(double, 2);
INSTANTIATION_HELPER(double, 5);
INSTANTIATION_HELPER(double, 6);
INSTANTIATION_HELPER(double, 7);
 */
