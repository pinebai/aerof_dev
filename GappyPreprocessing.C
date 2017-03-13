#include <GappyPreprocessing.h>
#ifdef DO_MODAL
#include <arpack++/include/ardsmat.h>
#include <arpack++/include/ardssym.h>
#endif

extern "C" {
  // Approximately solve the sparse non-negative least-squares problem
  //   min support(x) st ||A * x - b|| < reltol * ||b|| and x >= 0
  // Input: A is (mda x n), b is (m x 1), reltol is scalar
  // Output: A <- Q A, b <- Q b where Q is (m x m) orthogonal,
  //         rnorm <- ||b - Ax||_2,
  //         x is the (n x 1) primal solution,
  //         w <- A^T(b - Ax) is the (n x 1) dual solution
  // Work: zz is (m x 1), zz2 is (n x 1), index is (n x 1)
  // Info: mode: 1 => success, 2 => bad dim, 3 => too many iter
  void F77NAME(spnnls)(double *a, const long int *mda, const long int *m, const long int *n,
                        double *b, double *x, const double *reltol, double *rnorm, double *w,
                        double *zz, double *zz2, long int *index, long int *mode, long int *prtflg,
                        long int *sclflg, const double *maxsze, const double *maxite, double *dtime);
}


template<int dim>
GappyPreprocessing<dim>::GappyPreprocessing(Communicator *_com, IoData &_ioData, Domain
    &dom, DistGeoState *_geoState) :
  NonlinearRom<dim>(_com, _ioData, dom), 
  domain(dom), com(_com), ioData(&_ioData), 
  residual(0), jacobian(1), parallelRom(2),
  debugging(false), outputOnlineMatricesFull(false), outputOnlineMatricesSample(true),
  podRes(0, dom.getNodeDistInfo() ),
  podJac(0, dom.getNodeDistInfo() ),
  podHatRes(0, dom.getSampledNodeDistInfo() ),
  podHatJac(0, dom.getSampledNodeDistInfo() ),
  errorHatRes(0, dom.getSampledNodeDistInfo() ),
  errorHatJac(0, dom.getSampledNodeDistInfo() ),
  pseudoInvRhs(0, dom.getSampledNodeDistInfo() ),
  errorRes(0, dom.getNodeDistInfo() ),
  errorJac(0, dom.getNodeDistInfo() ),
  pseudoInverseMaskedSnapsTrans(0, dom.getNodeDistInfo() ),
  snapHatApproxMetric(0, dom.getNodeDistInfo() ),
  handledNodes(0), nPodBasis(0),
  // distribution info
  numLocSub(dom.getNumLocSub()), nTotCpus(_com->size()), thisCPU(_com->cpuNum()),
  nodeDistInfo(dom.getNodeDistInfo()), subD(dom.getSubDomain()),
  geoState(_geoState), X(_geoState->getXn())
{

  gappyIO = &(ioData->romOffline.gappy);
  approxMetricData = NULL;

  // create temporary objects to build postOp, which is needed for surfaces

  twoLayers = gappyIO->layers == 2;
  geoSourceTmp = new GeoSource(*ioData);
  tsDescTmp = new TsDesc<dim>(*ioData, *geoSourceTmp, &domain);
  bcDataTmp = tsDescTmp->createBcData(*ioData);
  varFcnTmp = new VarFcn(*ioData);
  postOp = new PostOperator<dim>(*ioData, varFcnTmp, bcDataTmp, geoState, &domain);
  input = new TsInput(_ioData);
  includeLiftFaces = gappyIO->includeLiftFaces;

  globalNodes = NULL;

  lowRankApproxMetricEigenvalues = NULL;

  numFullNodes = domain.getNumGlobNode();  // # globalNodes in full mesh
  com->fprintf(stdout," ... Number of full nodes in domain is %d ...\n",numFullNodes);

  if ((gappyIO->greedyLeastSquaresSolver == GappyConstructionData::GREEDY_LS_SCALAPACK) || 
      (gappyIO->pseudoInverseSolver == GappyConstructionData::PSEUDO_INVERSE_SCALAPACK)) {
    for(int i=0; i<2; ++i) parallelRom[i] = new ParallelRom<dim>(dom,_com,dom.getSampledNodeDistInfo());
  } else {
    for(int i=0; i<2; ++i) parallelRom[i] = NULL;
  }

  unionOfSampleNodes = -1;  // for readability
  globalSampleNodesForCluster.resize(0);
  globalSampleNodesUnion.resize(0);
  globalSampleNodesUnionForApproxMetricState.resize(0);
  globalSampleNodesUnionSet.clear();
  globalSampleNodesUnionSetForApproxMetricState.clear();

 /* d2wall = NULL;
  surfaceMask = NULL;
  targetRegionMask = NULL;
  if (gappyIO->minFractionOfSampledNodesOnSurfaceInTargetRegion > 0 ||
      gappyIO->minFractionOfSampledNodesInTargetRegion > 0) {
    if (ioData->input.d2wall[0] != 0) {
      d2wall = geoState->getd2wall();
      setUpSampledNodeTargetRegionMasks();    
    } else {
      com->fprintf(stdout,"*** ERROR: wall-biased sampled node selection has been requested, "
                          "but no WallDistance vector was supplied in the input file\n");
      exit(-1);     
    }
  }*/

  targetRegionMask = NULL;
  if (gappyIO->minFractionOfSampledNodesOnSurfaceInTargetRegion > 0 || gappyIO->minFractionOfSampledNodesInTargetRegion > 0) {
      setUpSampledNodeTargetRegionMask();    
  }

  wallMask = NULL;
  wallNeighborsMask = NULL;
  if (gappyIO->minFractionOfSampledNodesOnSurfaceInTargetRegion > 0) {
    wallMask = new DistVec<double>(domain.getNodeDistInfo());
    wallNeighborsMask = new DistVec<double>(domain.getNodeDistInfo());
    domain.setWallMask(*wallMask, *wallNeighborsMask);
    delete wallNeighborsMask;
    wallNeighborsMask = NULL;
    com->fprintf(stdout," ... wallMask norm = %e\n", wallMask->norm());
  }

  podTpod = NULL;

// bcFaces and bcFaceSurfID are deallocated in the write top file function
  int BC_CODE_EXTREME = max(BC_MAX_CODE, -BC_MIN_CODE)+1;  // there is a zero condition
  for (int i = 0; i < 2; ++i){
    for (int j = 0; j < 3; ++j)
      bcFaces[i][j] = new std::vector< int > [BC_CODE_EXTREME];
    bcFaceSurfID[i] = new std::vector< int > [BC_CODE_EXTREME];
  }

  surfaceMeshConstruction = false;

  cleanTempFiles=false;

 
}

//----------------------------------------------

template<int dim>
GappyPreprocessing<dim>::~GappyPreprocessing() 
{

  if (wallMask) delete wallMask;
  if (wallNeighborsMask) delete wallNeighborsMask;
  if (targetRegionMask) delete targetRegionMask;
  if (input) delete input;
  if (postOp) delete postOp;
  if (varFcnTmp) delete varFcnTmp;
  if (bcDataTmp) delete bcDataTmp;
  if (tsDescTmp) delete tsDescTmp;
  if (geoSourceTmp) delete geoSourceTmp;


  if (globalNodes) delete [] globalNodes;
  
  for(int i=0; i<2; ++i) {if (parallelRom[i]) delete parallelRom[i];}
}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::initialize() 
{
  // clear VecSets and reset counters/flags
  podRes.resize(0);
  podJac.resize(0);
  podHatRes.resize(0);
  podHatJac.resize(0);
  errorRes.resize(0);
  errorJac.resize(0);
  errorHatRes.resize(0);
  errorHatJac.resize(0);
  pseudoInvRhs.resize(0);
  pseudoInverseMaskedSnapsTrans.resize(0);
  snapHatApproxMetric.resize(0);
  nPodBasis = 0;
  initializeLeastSquaresDone = false;
  onlyInletOutletBC = true;  // first node should be on inlet/outlet
  handledNodes = 0;
  handledVectors[0] = 0;
  handledVectors[1] = 0;
  nRhs[0] = 0;
  nRhs[1] = 0;

  // reorient various pointers (just in case)
  pod.a[0] = &podRes;  // make pod point to res and jac
  pod.a[1] = &podJac;
  podHat.a[0] = &podHatRes;  // make pod point to res and jac
  podHat.a[1] = &podHatJac;
  error.a[0] = &errorRes;  // make pod point to res and jac
  error.a[1] = &errorJac;
  errorHat.a[0] = &errorHatRes;
  errorHat.a[1] = &errorHatJac;


  //if (globalNodes) delete [] globalNodes;
  //globalNodes = NULL;    // defined for union of sampled meshes (don't delete)

  // The following assignments should be unneccesary, since all of these pointers should already be NULL.  Better safe than sorry though -- wild pointers are very difficult to debug KMW

  nodesToHandle = NULL;  // allocated in setup, deallocated in determineSampleNodes
  
  for (int i = 0; i < 2; ++i) {
    podHatPseudoInv[i] = NULL;  //deallocated in outputOnlineMatrices
    nRhsGreedy[i] = NULL;  //deallocatd in determineSampleNodes
    onlineMatrices[i] = NULL;  //deallocated in outputOnlineMatrices (if needed)
  }

  // deallocated in buildMaps
  cpus = NULL;
  locSubDomains = NULL;
  localNodes = NULL;
  totalNodesCommunicated = NULL; 
  totalEleCommunicated = NULL;
  for (int i = 0; i < 3; ++i)
    nodesXYZ[i] = NULL;
  //elements = NULL;  // defined for union of sampled meshes (don't delete)
  for (int i = 0; i < 4; ++i)
    elemToNode[i] = NULL;

  cpuSample.clear();
  locSubSample.clear();
  locNodeSample.clear();
  globalSampleNodes.clear();
  reducedSampleNodes.clear();
  globalSampleNodeRankMap.clear();   // defined for a set of sample nodes in findMaxAndFillPodHat (in greedyIteration) 
                                     // used in defineMaps, outputTopFile, outputOnlineMats, (setSampleNodes)
  reducedSampleNodeRankMap.clear();  // defined for a set of sample nodes in outputTopFile
//  no need to delete following objects; defined for entire sample mesh and map only adds unique objects
//    globalNodeToCpuMap
//    globalNodeToLocSubDomainsMap
//    globalNodeToLocalNodesMap
//    nodesXYZmap
//    elemToNodeMap

}


//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::buildReducedModel() {

  this->freeMemoryForGappyPrepro();// clear unnecessary NonlinearRom.C objects

  double totalGappyOfflineTime = this->timer->getTime();

  if (gappyIO->doPrepro) {
    this->readClusterCenters("centers"); // double checks the value of nClusters that was given in the input file
    globalSampleNodesUnion.clear();
    globalSampleNodesForCluster.clear();
 
    double sampledMeshConstructionTime = this->timer->getTime(); // only used if sampled mesh is built
    double surfaceMeshConstructionTime = this->timer->getTime(); // only used if surface mesh is built
      
      //======================================
      // PURPOSE
      //   select sample globalNodes for each cluster.  Store sample nodes (mask) separately
      //   for each cluster, but also form union of all sample nodes for the construction
      //   of the sample mesh
      // INPUTS (for each cluster)
      //   nPod[0], nPod[1], nSampleNodes
      //   full domain decomposition: pod[0], pod[1]
      // OUTPUTS
      //   reduced domain decomposition: mesh, podHat[0], podHat[1],
      //======================================
  
      //======================================
      // NOTES
      // nRhsMax = number of POD vectors for each interpolation node that is
      // selected. It can be 1 to dim, where dim corresponds to interpolation.
      //
      // Fix the size of PhiHat using trick 1
      //
      // First compute a node, then compute the associated mask
      //
      //======================================
  
      //======================================
      // compute nRhsMax, nGreedyIt, nodesToHandle, possibly fix nSampleNodes
      //======================================

    if (gappyIO->selectSampledNodes || ioData->romOffline.rob.basisUpdates.preprocessForApproxUpdates || surfaceMeshConstruction) {

      globalSampleNodesForCluster.resize(this->nClusters);

      int nTargetSampleNodesForApproxMetricState;
      int nAddedSamplesLoc = 0;
      int nAddedSamplesGlob = 0;
      if (!surfaceMeshConstruction) {
        for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
          com->fprintf(stdout,"\n ... Selecting sample nodes for cluster %d ...\n", iCluster);
    
          initialize();
    
          bool breakloop=false;  // if using specified snapshots (not clustered)
          readGreedyData(iCluster,breakloop);
    
          setUpGreedy(iCluster);
    
          determineSampleNodes();  // use greedy algorithm to determine sample nodes
    
          globalSampleNodesForCluster[iCluster] = globalSampleNodes;  // store local sample nodes (the local mask)
        
          nTargetSampleNodesForApproxMetricState = floor(ioData->romOffline.rob.basisUpdates.approxMetricState.sampledMeshFraction*double(globalSampleNodes.size()));
          nAddedSamplesGlob = globalSampleNodesUnionSet.size();   
          for (int iNode=0; iNode<int(globalSampleNodes.size()); iNode++) { 
            globalSampleNodesUnionSet.insert(globalSampleNodes[iNode]);  // add local sample nodes union
            nAddedSamplesLoc = int(globalSampleNodesUnionSet.size())-nAddedSamplesGlob;
            if (nAddedSamplesLoc<=nTargetSampleNodesForApproxMetricState)
              globalSampleNodesUnionSetForApproxMetricState.insert(globalSampleNodes[iNode]);
          }

          if (breakloop) {
            for (int jCluster=iCluster; jCluster<(this->nClusters); ++jCluster) {
              globalSampleNodesForCluster[jCluster] = globalSampleNodes;
            }
            break; 
          }
        
        } 
    
        // place the union of the sample nodes into a vector
        for (std::set<int>::iterator it=globalSampleNodesUnionSet.begin(); it!=globalSampleNodesUnionSet.end(); ++it) {
          globalSampleNodesUnion.push_back(*it);
        }
        globalSampleNodesUnionSet.clear();
    
        for (std::set<int>::iterator it=globalSampleNodesUnionSetForApproxMetricState.begin(); it!=globalSampleNodesUnionSetForApproxMetricState.end(); ++it) {
          globalSampleNodesUnionForApproxMetricState.push_back(*it);
        }
        globalSampleNodesUnionSetForApproxMetricState.clear();
     
      }
    }

    domain.makeSampledNodeDistInfo(globalSampleNodesUnion, globalNodeToCpuMap, globalNodeToLocSubDomainsMap);

    initialize();
 
    setSampleNodes(unionOfSampleNodes); // also calls buildRemainingMesh()
  
    if (surfaceMeshConstruction) {
      this->timer->addSurfaceMeshConstructionTime(surfaceMeshConstructionTime);
    } else {
      this->timer->addSampledMeshConstructionTime(sampledMeshConstructionTime);
    }
  
    if (thisCPU == 0) outputTopFile(unionOfSampleNodes);

    outputMatchStateReduced();  
    outputInitialConditionReduced();
    outputMultiSolutionsReduced();
    outputWallDistanceReduced();  // distributed info (parallel)
    outputDisplacementReduced();
    outputShapeDerivativeReduced();
 
    for (int iCluster=0; iCluster<(this->nClusters); ++iCluster) {
  
      initialize();
  
      if (surfaceMeshConstruction) {
        com->fprintf(stdout,"\n-----------------------------------------------------\n");
        com->fprintf(stdout," ... Writing surface quantities for cluster %d\n\n", iCluster);
   
        setSampleNodes(unionOfSampleNodes);
        outputLocalStateBasisReduced(iCluster);
        outputLocalReferenceStateReduced(iCluster);
  
      } else {
        com->fprintf(stdout,"\n-----------------------------------------------------\n");
        com->fprintf(stdout," ... Writing sampled quantities for cluster %d\n\n", iCluster);
   
        int sampleNodes = (gappyIO->useUnionOfSampledNodes==GappyConstructionData::UNION_TRUE) ? unionOfSampleNodes : iCluster; 
        setSampleNodes(sampleNodes);

        if ((thisCPU == 0) && gappyIO->selectSampledNodes) outputSampleNodes(iCluster);
        outputLocalStateBasisReduced(iCluster);
        outputLocalReferenceStateReduced(iCluster);
 
      }
    }
  
    this->freeMemoryForGappyPrepro();// clear unnecessary NonlinearRom.C objects 
  
    if (surfaceMeshConstruction) {
  
      com->fprintf(stdout," \n... Finished with Surface Mesh Construction\n");
  
    } else { // build and output online matrices

      if (gappyIO->doPreproGNAT) {  
        for (int iCluster=gappyIO->initialCluster; iCluster<(this->nClusters); ++iCluster) {
      
          com->fprintf(stdout,"\n-----------------------------------------------------\n");
          com->fprintf(stdout," ... Computing online GNAT matrices for cluster %d\n\n", iCluster);
      
          initialize();                  // deallocate for everything but the sampled nodes (local and union)
      
          int sampleNodes = (gappyIO->useUnionOfSampledNodes==GappyConstructionData::UNION_TRUE) ? unionOfSampleNodes : iCluster; 
          setSampleNodes(sampleNodes);   // use local sampled nodes
      
          setUpPodResJac(iCluster);      // read in Res/Jac bases
       
          setUpBasisBasisProducts();     // computes ROB[1]T*ROB[0] and ROB[0]T*ROB[0] if necessary
      
          formMaskedNonlinearROBs();
      
          initializeLeastSquares();
       
          double pseudoInvTime = this->timer->getTime();
          computePseudoInverse();        // only requires sample nodes and masked nonlinear ROBs
          this->timer->addPseudoInvTime(pseudoInvTime);
      
          // STRATEGY:
          // - put zeros in the pod[0], pod[1] for any node that is not part of the local mask
          // - don't evaluate rhat and jphihat for nodes that are not part of the local mask
      
          // TRICK: adding zero rows to matrices has no effect on the qr decomposition
      
          if (thisCPU == 0) {
            assembleOnlineMatrices();   // handle online matrices so you can free up pseudo inverse matrix
          }
      
          outputOnlineMatrices(iCluster);
      
          if (podTpod) {
            for (int iVec=0; iVec<nPod[1]; ++iVec)
              delete [] podTpod[iVec];
            delete podTpod;
            podTpod = NULL;
          }
      
        }
      } // end GNAT online matrix preprocessing

      if (gappyIO->doPreproApproxMetricNonlinear || gappyIO->doPreproApproxMetricNonlinearNNLS ) {
        initialize();
        this->freeMemoryForGappyPrepro();

        // TODO if union == true and reading snaps from ascii file, only need to do this once
        std::vector<std::vector<double> >* corrMat = NULL;
        for (int iCluster=gappyIO->initialCluster; iCluster<(this->nClusters); ++iCluster) {
          setSampleNodes(iCluster); 
          com->fprintf(stdout,"\nPreprocessing for approximate gappy solver\n");
          double approxMetricTime = this->timer->getTime(); 
  
          constructApproximatedMetric("nonlinear",iCluster,corrMat);

          //if (gappyIO->testApproxMetric) {
          //  testInnerProduct("approxMetricNonlinear");
          //  testInnerProduct("nonlinear");
          //}
          this->timer->addApproxMetricPreproTime(approxMetricTime);
        }
        if (corrMat) delete corrMat;
      } // end Approx Metric Nonlinear preprocessing
     
    }

    setSampleNodes(unionOfSampleNodes);

    if (this->ioData->romOffline.rob.basisUpdates.preprocessForApproxUpdates) {
      this->freeMemoryForGappyPrepro();// clear unnecessary NonlinearRom.C objects
      initialize(); 
 
      if (!surfaceMeshConstruction) {
        com->fprintf(stdout,"\nPreprocessing for approximate basis updates\n");
        double approxUpdatesPreproTime = this->timer->getTime(); 
  
        outputClusterCentersReduced();
  
        constructApproximatedMetric("state");
  
        if (gappyIO->testApproxMetric) {
          testInnerProduct("approxMetricState");
          testInnerProduct("state");
        }
        this->timer->addApproxMetricPreproTime(approxUpdatesPreproTime);
      } else {
        com->fprintf(stdout,"\nOutputting approximated metric in surface coordinates\n");
        this->readApproxMetricStateLowRankFactor("full");
        outputApproxMetricLowRankFactorReducedCoords("state");
      }
    }

    if (surfaceMeshConstruction) outputClusterCentersReduced();

  } // end if(gappyIO->doPrepro)

  com->fprintf(stdout," \n ... Finished with preprocessing calculations ...\n");

  this->freeMemoryForGappyPrepro();

  if (gappyIO->sowerInputs) this->partitionAndSowerForGappy(surfaceMeshConstruction);

  this->timer->addTotalGappyOfflineTime(totalGappyOfflineTime);

} 

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::constructApproximatedMetric(const char* type, int iCluster, std::vector<std::vector<double> >* corrMat) {
  // the notations follow the paper "Fast Local Reduced Bais Updates Based on Approximated Metric"

  if (strcmp(type,"state")==0) {
    com->fprintf(stdout,"\nComputing approximate metric for inner products of gappy state vectors\n");
    approxMetricData = &(ioData->romOffline.rob.basisUpdates.approxMetricState);
    this->readSnapshotFiles("approxMetricState", 0);
  } else if (strcmp(type,"nonlinear")==0) {
    com->fprintf(stdout,"\nComputing approximate metric for inner products of gappy nonlinear quantities (cluster %d)\n", iCluster);
    approxMetricData = &(ioData->romOffline.gappy.approxMetricNonlinear);
    if (strcmp(this->ioData->input.approxMetricNonlinearSnapFile,"")!=0) {
      this->readSnapshotFiles("approxMetricNonlinear", 0);
      if (corrMat==NULL) {
        corrMat = new std::vector<std::vector<double> >;
        corrMat->resize(this->snap->numVectors());
        for (int iSnap=0; iSnap<this->snap->numVectors(); ++iSnap) {
          (*corrMat)[iSnap].resize(iSnap+1,0.0); //lower triangular
          for (int jSnap=0; jSnap<=iSnap; ++jSnap) {
            (*corrMat)[iSnap][jSnap] = -1; //(*this->snap)[iSnap] * (*this->snap)[jSnap];
          }
          (*corrMat)[iSnap][iSnap] = (*this->snap)[iSnap] * (*this->snap)[iSnap];
        }
      }
    } else {
      this->readClusteredSnapshots(iCluster, true, "residual", 0, ioData->romOffline.gappy.maxClusteredSnapshotsNonlinearApproxMetric);
      if (corrMat) delete corrMat;
      corrMat = new std::vector<std::vector<double> >;
      corrMat->resize(this->snap->numVectors());
      for (int iSnap=0; iSnap<this->snap->numVectors(); ++iSnap) {
        (*corrMat)[iSnap].resize(iSnap+1,0.0); //lower triangular
        for (int jSnap=0; jSnap<=iSnap; ++jSnap) {
          (*corrMat)[iSnap][jSnap] = -1; //(*this->snap)[iSnap] * (*this->snap)[jSnap];
        }
        (*corrMat)[iSnap][iSnap] = (*this->snap)[iSnap] * (*this->snap)[iSnap];
      }
    }
  }

  // Step 1: compute an EVD of the correlation matrix X'*X and truncate the EVD based on the decay of EV
  if (strcmp(type,"state")==0 || (gappyIO->doPreproApproxMetricNonlinear && strcmp(type,"nonlinear")==0))
    computeCorrelationMatrixEVD(corrMat);

  // Step 2: compute the pseudo-inverse of the masked snapshots Y
  computePseudoInverseMaskedSnapshots(type, iCluster);

  // Step 3: compute the low rank factor G of the approximated metric
  if (strcmp(type,"state")==0 || (gappyIO->doPreproApproxMetricNonlinear && strcmp(type,"nonlinear")==0))
    computeApproximatedMetricLowRankFactor();

  // Step 4: output low rank factor
  if (strcmp(type,"state")==0 || (gappyIO->doPreproApproxMetricNonlinear && strcmp(type,"nonlinear")==0)) {
    outputApproxMetricLowRankFactorReducedCoords(type, iCluster);
    outputApproxMetricLowRankFactorFullCoords(type, iCluster);
  }
}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::computeCorrelationMatrixEVD(std::vector<std::vector<double> >* corrMat) {

#ifdef DO_MODAL
  numSnapsForApproxMetric = this->snap->numVectors();
  com->fprintf(stderr, " ... Forming correlation matrix using %d snapshots\n",numSnapsForApproxMetric);
  // allocate for upper half of sym. eigprob
  double *rVals = new double[numSnapsForApproxMetric*(numSnapsForApproxMetric+1)/2];
  for (int i = 0; i < numSnapsForApproxMetric; i++)
    for (int j = 0; j <= i; j++)
      rVals[(i+1)*i/2 + j] = (corrMat) ? (*corrMat)[i][j] : ((*this->snap)[j]) * ((*this->snap)[i]);

  double tolerance = approxMetricData->tolerance;

  ARdsSymMatrix<double> corr(numSnapsForApproxMetric, rVals, 'U');

  com->fprintf(stderr, " ... Factoring correlation matrix\n");

  corr.FactorA();
  int ncv = numSnapsForApproxMetric-1;
  ARluSymStdEig<double> corrEigProb(ncv, corr, "LM", ncv+1, tolerance, 300*ncv);

  com->fprintf(stderr, " ... Solving eigenproblem\n");
  int nconv = corrEigProb.FindEigenvectors();

  com->fprintf(stderr, " ... Found %d converged eigenvectors out of %d snaps\n", nconv, numSnapsForApproxMetric);

  //Test eigenvalue decomposition
  double errrValsEV, maxErr, avgErr;
/*  maxErr = 0.0;
  avgErr = 0.0;
  for (int i = 0; i < numSnapsForApproxMetric; i++) {
    for (int j = 0; j <= i; j++) {
      errrValsEV = 0.0;
      for (int k = 0; k < nconv; ++k)      
        errrValsEV += (corrEigProb.Eigenvalue(nconv-k-1)*corrEigProb.Eigenvector(nconv-k-1, i)*corrEigProb.Eigenvector(nconv-k-1, j));
      errrValsEV = abs(rVals[(i+1)*i/2 + j] - errrValsEV) / abs(rVals[(i+1)*i/2 + j]);
      avgErr += errrValsEV;
      if (errrValsEV > maxErr)
        maxErr = errrValsEV;
    }
  }
  avgErr /= (numSnapsForApproxMetric*(numSnapsForApproxMetric+1)/2);
  com->fprintf(stderr, " ... Average error on EVD = %e\n", avgErr);
  com->fprintf(stderr, " ... Maximum error on EVD = %e\n", maxErr); */

  //Retain a low rank approximation
  double totalEnergy = 0.0;
  for (int i = 0; i < nconv; ++i) {
    //com->fprintf(stderr, "Eig[%d] = %f\n",i,corrEigProb.Eigenvalue(nconv-i-1));
    totalEnergy += corrEigProb.Eigenvalue(nconv-i-1);  
  }

  double energyRetained = totalEnergy*approxMetricData->lowRankEnergy;
  totalEnergy = 0.0;
  numEigen = 0;
  while (totalEnergy < energyRetained && numEigen< nconv-1) {
   totalEnergy += corrEigProb.Eigenvalue(nconv-numEigen-1);
   numEigen++;  
  }
  // retain first numEigen eigenmodes
  com->fprintf(stderr, " ... Retaining the first %d eigenvectors out of %d\n", numEigen, nconv);
  lowRankApproxMetricEigenvalues  = new double[numEigen];
  for (int j = 0; j < numEigen; ++j)
    lowRankApproxMetricEigenvalues[j] = corrEigProb.Eigenvalue(nconv-j-1);
  lowRankModes = new double*[numEigen];
  for (int j = 0; j < numEigen; ++j)
    lowRankModes[j] = new double[numSnapsForApproxMetric];
  for (int j = 0; j < numEigen; ++j)
    for (int k = 0; k < numSnapsForApproxMetric; ++k)
      lowRankModes[j][k] = corrEigProb.Eigenvector(nconv-j-1, k)*pow(lowRankApproxMetricEigenvalues[j],0.5);//scale by the square root of the eigenvalues

  //Test eigenvalue decomposition
  if (gappyIO->testApproxMetric) {
    maxErr = 0.0;
    avgErr = 0.0;
    for (int i = 0; i < numSnapsForApproxMetric; i++) {
      for (int j = 0; j <= i; j++) {
        errrValsEV = 0.0;
        for (int k = 0; k < numEigen; ++k)
          errrValsEV += (corrEigProb.Eigenvalue(nconv-k-1)*corrEigProb.Eigenvector(nconv-k-1, i)*corrEigProb.Eigenvector(nconv-k-1, j));
        errrValsEV = abs(rVals[(i+1)*i/2 + j] - errrValsEV) / abs(rVals[(i+1)*i/2 + j]);
        avgErr += errrValsEV;
        if (errrValsEV > maxErr)
          maxErr = errrValsEV;
      }
    }
    avgErr /= (numSnapsForApproxMetric*(numSnapsForApproxMetric+1)/2);
    com->fprintf(stderr, " ... Average error on EVD after truncation= %e\n", avgErr);
    com->fprintf(stderr, " ... Maximum error on EVD after truncation= %e\n", maxErr);
  }

  delete [] rVals;

#else
  com->fprintf(stderr, "*** ERROR: REQUIRES COMPILATION WITH ARPACK and DO_MODAL Flag\n");
  com->barrier();
  exit(-1);
#endif

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::computePseudoInverseMaskedSnapshots(const char* type, int iCluster) {

  if (strcmp(type,"nonlinear")==0 && gappyIO->doPreproApproxMetricNonlinearNNLS) {
    computeApproxMetricNonlinearNNLS(iCluster);
  } else if (strcmp(type,"state")==0 || (gappyIO->doPreproApproxMetricNonlinear && strcmp(type,"nonlinear")==0)) {
    computeMaskedSnapshots(type, iCluster);
    computePseudoInverseTranspose();
  } else {
    this->com->fprintf(stderr, "*** Error in computePseudoInverseMaskedSnapshots\n");
    sleep(1);
    exit(-1);
  }

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::computeApproxMetricNonlinearNNLS(int iCluster) {

  /* Adapted by Julien Cortial at Stanford University

    GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE A
    SPARSE N-VECTOR, X, THAT VERIFIES
  
        ||A * X - B||_2 <= RELTOL * ||B||_2  SUBJECT TO X .GE. 0   
    ------------------------------------------------------------------
                    Subroutine Arguments
           
    A(),MDA,M,N     MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE   
                    ARRAY, A().   ON ENTRY A() CONTAINS THE M BY N    
                    MATRIX, A.           ON EXIT A() CONTAINS 
                    THE PRODUCT MATRIX, Q*A , WHERE Q IS AN   
                    M BY M ORTHOGONAL MATRIX GENERATED IMPLICITLY BY  
                    THIS SUBROUTINE.  
    B()     ON ENTRY B() CONTAINS THE M-VECTOR, B.   ON EXIT B() CON- 
            TAINS Q*B. 
    X()     ON ENTRY X() NEED NOT BE INITIALIZED.  ON EXIT X() WILL   
            CONTAIN THE SOLUTION VECTOR. 
    RELTOL  RELATIVE TOLERANCE
            (STOPPING CRITERION: ||B - A*X||_2 < RELTOL * ||B||_2).
    RNORM   ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE  
            RESIDUAL VECTOR.  
    W()     AN N-ARRAY OF WORKING SPACE.  ON EXIT W() WILL CONTAIN    
            THE DUAL SOLUTION VECTOR.   W WILL SATISFY W(I) = 0.  
            FOR ALL I IN SET P  AND W(I) .LE. 0. FOR ALL I IN SET Z   
    ZZ()     AN M-ARRAY OF WORKING SPACE.     
    ZZ2()    AN N-ARRAY OF WORKING SPACE.     
    INDEX()     AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N.
                ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS    
                P AND Z AS FOLLOWS..  
  
                INDEX(1)   THRU INDEX(NSETP) = SET P.     
                INDEX(IZ1) THRU INDEX(IZ2)   = SET Z.     
                IZ1 = NSETP + 1 = NPP1
                IZ2 = N   
    MODE    THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING 
            MEANINGS. 
            1     THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY.
            2     THE DIMENSIONS OF THE PROBLEM ARE BAD.  
                  EITHER M .LE. 0 OR N .LE. 0.
            3     ITERATION COUNT EXCEEDED.
                  MORE THAN MAXITE*N ITERATIONS. 
    MAXSZE  ALTERNATIVE STOPPING CRITERION BASED ON ACTIVE SET SIZE
            NSETP <= MIN(M,MAXSZE*N)
  
    ------------------------------------------------------------------
     SUBROUTINE SPNNLS (A,MDA,M,N,B,X,RELTOL,RNORM,W,ZZ,ZZ2,INDEX,MODE,
    +                   PRTFLG,SCAFLG,MAXSZE,MAXITE,DTIME)  

  void F77NAME(spnnls)(double *a, const long int *mda, const long int *m, const long int *n,
                       double *b, double *x, const double *reltol, double *rnorm, double *w,
                       double *zz, double *zz2, long int *index, long int *mode, long int *prtflg,
                       long int *sclflg, const double *maxsze, const double *maxite, double *dtime);*/

  for (int iSnap=0; iSnap<this->snap->numVectors(); ++iSnap) {
    double norm = (*this->snap)[iSnap].norm();
    if (norm>0) (*this->snap)[iSnap] *= 1/norm;
  }

  long int nSnapsTrain = int(floor(double(this->snap->numVectors())/2.0));
  long int nSnapsTest = this->snap->numVectors();
  long int nVars = dim*nSampleNodes;
  long int nEqns = dim*nSampleNodes + nSnapsTrain;

  double *lhs = new double[nEqns*nVars];
  double *lhsUpper = new double[nSnapsTrain*nVars]; // for minimal communication
  double *lhsTest = new double[nSnapsTest*nVars];
  double *rhs = new double[nEqns];
  double *rhsTest = new double[nSnapsTest];
  double *sol = new double[nVars];
  double *dualSol = new double[nVars];
  double *work1 = new double[nEqns];
  double *work2 = new double[nVars];
  long int *index = new long int[nVars];

  // initialize everything 
  for (long int iEqn=0; iEqn<nEqns; ++iEqn) {
    work1[iEqn] = 0.0;
    rhs[iEqn] = 0.0;
    if (iEqn < nSnapsTrain) { 
      for (long int iVar=0; iVar<nVars; ++iVar) {
        lhsUpper[iEqn + iVar*nSnapsTrain] = 0.0;
      }
    }
    if (iEqn < nSnapsTest) { 
      for (long int iVar=0; iVar<nVars; ++iVar) {
        lhsTest[iEqn + iVar*nSnapsTest] = 0.0;
      }
    }
  }
  for (long int iVar=0; iVar<nVars; ++iVar) {
    sol[iVar] = 0.0;
    dualSol[iVar] = 0.0;
    work2[iVar] = 0.0;
    index[iVar] = 0;
  }
 
  // fill the top of lhs (masked snapshots)
  int nTotalDOF = 0;
  for (int iSampleNode=0; iSampleNode<nSampleNodes; ++iSampleNode) {
    int currentSampleNode = globalSampleNodes[iSampleNode];
    int locSub = globalNodeToLocSubDomainsMap.find(currentSampleNode)->second;
    int locNode = globalNodeToLocalNodesMap.find(currentSampleNode)->second;
    int currentCpu = globalNodeToCpuMap.find(currentSampleNode)->second;
    if (thisCPU == currentCpu) {
      for (long int iEqn=0; iEqn<nSnapsTrain; ++iEqn) {
        SubDomainData<dim> locSnap = (*this->snap)[iEqn*2].subData(locSub); 
        for (int iDim=0; iDim<dim; ++iDim) {
          long int iVar = dim*iSampleNode + iDim;
          lhsUpper[iEqn + iVar*nSnapsTrain] = pow(locSnap[locNode][iDim],2);
          ++nTotalDOF;
        }
      }
    }
  }

  // communicate
  com->globalSum(nSnapsTrain*nVars, lhsUpper); 
  com->globalSum(1,&nTotalDOF);
  this->com->fprintf(stdout, "expected %d DOF, found %d DOF\n",dim*nSampleNodes*nSnapsTrain,nTotalDOF); 
 
  // also fill test data (superset of masked snapshots)
  for (int iSampleNode=0; iSampleNode<nSampleNodes; ++iSampleNode) {
    int currentSampleNode = globalSampleNodes[iSampleNode];
    int locSub = globalNodeToLocSubDomainsMap.find(currentSampleNode)->second;
    int locNode = globalNodeToLocalNodesMap.find(currentSampleNode)->second;
    int currentCpu = globalNodeToCpuMap.find(currentSampleNode)->second;
    if (thisCPU == currentCpu) {
      for (long int iEqn=0; iEqn<nSnapsTest; ++iEqn) {
        SubDomainData<dim> locSnap = (*this->snap)[iEqn].subData(locSub);
        for (int iDim=0; iDim<dim; ++iDim) {
          long int iVar = dim*iSampleNode + iDim;
          lhsTest[iEqn + iVar*nSnapsTest] = pow(locSnap[locNode][iDim],2);
        }
      }
    }
  }

  for (long int iEqn=0; iEqn<nSnapsTest; ++iEqn) {
      rhsTest[iEqn] = 1.0;  // snapshots were normalized above
  }
 
  // communicate test data
  com->globalSum(nSnapsTest*nVars, lhsTest);  

  double nnlsTime = this->timer->getTime();

  // pick lambda
  double lambdaMin = -12;
  double lambdaMax = 4;
  double lambda = lambdaMin + (lambdaMax-lambdaMin)*double(thisCPU)/double(nTotCpus);

  // fill lhs and rhs
  for (long int iEqn=0; iEqn<nEqns; ++iEqn) {
    if (iEqn < nSnapsTrain) {
      rhs[iEqn] = 0.0; 
      for (long int iVar=0; iVar<nVars; ++iVar) {
        lhs[iEqn + iVar*nEqns] = lhsUpper[iEqn + iVar*nSnapsTrain];
        rhs[iEqn] -= lhsUpper[iEqn + iVar*nSnapsTrain];  // due to change of variables: q = x+1, x>=0
      }
    } else {
      rhs[iEqn] = 0.0;
      for (long int iVar=0; iVar<nVars; ++iVar) {
        lhs[iEqn + iVar*nEqns] = 0.0;
        if ((iEqn - nSnapsTrain) == iVar) lhs[iEqn + iVar*nEqns] = pow(10,lambda);
      }
    }
  }
  delete [] lhsUpper;

  for (long int iEqn=0; iEqn<nSnapsTrain; ++iEqn) {
    rhs[iEqn] += 1.0;  // snapshots were normalized above
  }

  /*
  for (long int iEqn=0; iEqn<nEqns; ++iEqn) {
    for (long int iVar=0; iVar<nVars; ++iVar) {
      this->com->fprintf(stdout, "lhs[%d,%d]=%e ",iEqn,iVar,lhs[iEqn + iVar*nEqns]); 
    }
    this->com->fprintf(stdout, "\n");
  }
  this->com->fprintf(stdout, "\n");

  for (long int iEqn=0; iEqn<nEqns; ++iEqn) {
    this->com->fprintf(stdout, "rhs[%d]=%e\n",iEqn,rhs[iEqn]); 
  }
  this->com->fprintf(stdout, "\n");
 
  for (long int iEqn=0; iEqn<nSnapsTest; ++iEqn) {
    for (long int iVar=0; iVar<nVars; ++iVar) {
      this->com->fprintf(stdout, "lhsTest[%d,%d]=%e ",iEqn,iVar,lhsTest[iEqn + iVar*nSnapsTest]); 
    }
    this->com->fprintf(stdout, "\n");
  }
  this->com->fprintf(stdout, "\n");*/ 

  // problem set up
  long int printFlag = 0; // zero for silent
  long int scaleFlag = 1; //zero for no scaling
  long int status = -1;   // success/failure flag
  double errorMag = -1; // norm of least squares error
  double relTol = 1e-16;
  double maxSizeRatio = 1.0;
  double maxItsRatio = 1.0;
  double dtime = 0.0;

  F77NAME(spnnls)(lhs, &nEqns, &nEqns, &nVars, rhs, sol, &relTol, &errorMag, dualSol,
                  work1, work2, index, &status, &printFlag, &scaleFlag, &maxSizeRatio, &maxItsRatio, &dtime);

  switch (status) {
    case (1):
      com->fprintf(stdout, "... success: NNLS converged\n");
      break;
    case (2):
      fprintf(stderr, "*** Error: Illegal dimensions passed to NNLS! (CPU %d, nEqns=%ld, nVars=%ld)\n", thisCPU, nEqns, nVars);
      sleep(1);
      exit(-1);
      break;
    case (3):
      fprintf(stderr, "*** Warning: NNLS hit max iterations! (CPU %d, errorMag= %e\n", thisCPU, errorMag);
      break;
    default:
      fprintf(stderr, "*** Error: Unexpected status reported by NNLS (CPU %d)\n", thisCPU);
      sleep(1); 
      exit(-1);
  }

  // find min error
  for (long int iVar=0; iVar<nVars; ++iVar) {
    sol[iVar] += 1.0;
    //this->com->fprintf(stdout, "... sol[%d]=%e\n",iVar,sol[iVar]); 
  }

  // product = lhsTest * sol ...  ( error = product - rhsTest)
  double *product = new double[nSnapsTest];
  for (long int iEqn=0; iEqn<nSnapsTest; ++iEqn) {
    product[iEqn] = 0.0;
    for (long int iVar=0; iVar<nVars; ++iVar) {
      product[iEqn] += lhsTest[iEqn + iVar*nSnapsTest] * sol[iVar];
    }
  }
 
  // normalize product and rhsTest to eliminate any scaling issues (which won't affect online performance) 
  double productNorm = 0.0;
  double rhsNorm = 0.0;
  for (long int iEqn=0; iEqn<nSnapsTest; ++iEqn) {
    productNorm += pow(product[iEqn],2);
    rhsNorm += pow(rhsTest[iEqn],2);
  }
  productNorm = pow(productNorm,0.5);
  rhsNorm = pow(rhsNorm,0.5);

  double* error = new double[nTotCpus];
  int* nonzeros = new int[nTotCpus];
  double* nnlsTimes = new double[nTotCpus];
  double* lambdas = new double[nTotCpus];
  for (int iCPU = 0; iCPU<nTotCpus; ++iCPU) {
    error[iCPU] = 0.0;
    nonzeros[iCPU] = 0.0;
    nnlsTimes[iCPU] = 0.0;
    lambdas[iCPU] = 0.0;
  }

  nnlsTimes[thisCPU] = this->timer->getTime() - nnlsTime;
  lambdas[thisCPU] = pow(10,lambda);

  for (long int iVar=0; iVar<nVars; ++iVar)
    nonzeros[thisCPU] += (sol[iVar]==1.0) ? 0 : 1;

  for (long int iEqn=0; iEqn<nSnapsTest; ++iEqn)
    error[thisCPU] += pow((product[iEqn]/productNorm) - (rhsTest[iEqn]/rhsNorm),2);
  error[thisCPU] = pow(error[thisCPU],0.5);

  com->globalSum(nTotCpus, error);  
  com->globalSum(nTotCpus, nonzeros);
  com->globalSum(nTotCpus, nnlsTimes);
  com->globalSum(nTotCpus, lambdas);
  com->barrier();

  for (int iCPU = 0; iCPU<nTotCpus; ++iCPU) {
    this->com->fprintf(stdout, "... CPU %d | lambda: %e | nonzero: %d/%d | error: %e | time %e\n", 
                              iCPU, lambdas[iCPU], nonzeros[iCPU], nVars, error[iCPU], nnlsTimes[iCPU]);
  }

  int minValCPU = nTotCpus-1;
  for (int iCPU = nTotCpus-1; iCPU >= 0; iCPU--){
    if ((error[iCPU] <= error[minValCPU]) && (nonzeros[iCPU] > (double(nVars)*0.25))) {
      minValCPU = iCPU;
    //} else {
    //  break;
    }            
  }
  double minError = error[minValCPU];
  this->com->fprintf(stdout, "... CPU %d selected (error %e, nonzero %d/%d)\n",
                            minValCPU, minError, nonzeros[minValCPU], nVars);

  // form metric
  com->broadcast(nVars,sol,minValCPU);

  DistSVec<double,dim> metric(this->domain.getNodeDistInfo());
  metric = 0.0;
  for (int iSampleNode=0; iSampleNode<nSampleNodes; ++iSampleNode) {
    int currentSampleNode = globalSampleNodes[iSampleNode];
    int locSub = globalNodeToLocSubDomainsMap.find(currentSampleNode)->second;
    int locNode = globalNodeToLocalNodesMap.find(currentSampleNode)->second;
    int currentCpu = globalNodeToCpuMap.find(currentSampleNode)->second;
    if (thisCPU == currentCpu) {
      SubDomainData<dim> locMetric = metric.subData(locSub);
      for (int iDim=0; iDim<dim; ++iDim) {
        locMetric[locNode][iDim] = sol[dim*iSampleNode + iDim];
      }
    }
  }

  // store metric
  com->fprintf(stdout,"\n ... Writing NNLS approximated metric ...\n");

  char *filePath = NULL;
  this->determinePath(this->approxMetricNonlinearName, iCluster, filePath);

  std::string header("Vector ApproxMetric under load for FluidNodesRed");

  FILE* myOutFile = NULL;
  if (thisCPU==0) {
    myOutFile = fopen(filePath, "wt");
    fprintf(myOutFile,"%s\n", header.c_str());
    fprintf(myOutFile,"%d\n", nReducedNodes);
  }

  outputReducedSVec(metric, myOutFile, 0.0);

  if (thisCPU==0) fclose(myOutFile);

  if (filePath) {
    delete [] filePath;
    filePath = NULL;
  }

  // clean up
  delete [] lhs;
  delete [] lhsTest;
  delete [] rhs;
  delete [] rhsTest;
  delete [] sol;
  delete [] product;
  delete [] dualSol;
  delete [] work1;
  delete [] work2;
  delete [] index;
  delete [] error;
  delete [] nonzeros;
  delete [] nnlsTimes;
  delete [] lambdas;
}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::computeMaskedSnapshots(const char* type, int iCluster) {

  //initialize
  int nCols = this->snap->numVectors();
  snapHatApproxMetric.resize(nCols);
  for (int iCol = 0; iCol < nCols; ++iCol) 
    snapHatApproxMetric[iCol] = 0.0;
  
  // fill
  std::vector<int> approxMask; 
  if (strcmp(type,"state")==0) {
    approxMask = globalSampleNodesUnionForApproxMetricState;
  } else if (strcmp(type,"nonlinear")==0) {
    approxMask = globalSampleNodes;
  } else { 
    exit(-1);
  }

  int nNodesApproxMetric = approxMask.size();
  for (int iNode=0; iNode < nNodesApproxMetric; ++iNode) {
    int currentNode = approxMask[iNode];
    int locSub = globalNodeToLocSubDomainsMap.find(currentNode)->second;
    int locNode = globalNodeToLocalNodesMap.find(currentNode)->second;
    int currentCpu = globalNodeToCpuMap.find(currentNode)->second;
    if (thisCPU == currentCpu) {
      // fill out sampled matrices (all columns for the current rows)
      SubDomainData<dim> locSnap, locSnapHat;
      for (int iCol = 0; iCol < nCols; ++iCol) {
        locSnap = (*this->snap)[iCol].subData(locSub);  // cannot access iDim entry
        locSnapHat = snapHatApproxMetric[iCol].subData(locSub);
        for (int iDim = 0; iDim < dim ; ++iDim) {
          locSnapHat[locNode][iDim] = locSnap[locNode][iDim];
          // zeros everywhere except at sample nodes
        }
      }
    }
    com->barrier();
  }

  this->freeMemoryForGappyPrepro();


}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::computePseudoInverseTranspose() {

  // svd of the masked snapHatApproxMetricshots USV'
  int nSnaps = snapHatApproxMetric.numVectors();
  SetOfVec UTrueTmp(nSnaps, this->domain.getNodeDistInfo());
  Vec<double> singValsTmp(nSnaps);
  FullM VtrueTmp(nSnaps);
  int nKeep = nSnaps;

  ParallelRom<dim> parallelRom( this->domain, this->com, this->domain.getNodeDistInfo());
  parallelRom.parallelSVD(snapHatApproxMetric, UTrueTmp, singValsTmp.data(), VtrueTmp, nSnaps, true);

  // check the svd
  // (looks good: errors less than 1e-15)  Note that Vtrue is not Vtranspose
  double errorNorm,maxErr,avgErr;
  DistSVec<double,dim> errorVec( domain.getNodeDistInfo() );
  /*maxErr = 0.0;
  avgErr = 0.0;
  for (int iSnap = 0; iSnap <nSnaps; ++iSnap) {
    errorVec = snapHatApproxMetric[iSnap];
    for (int jSnap = 0; jSnap < nSnaps; ++jSnap)
      errorVec = errorVec - ((singValsTmp[jSnap]*VtrueTmp[iSnap][jSnap])*UTrueTmp[jSnap]);
    errorNorm = ( ((snapHatApproxMetric[iSnap]).norm())> 1e-15) ? errorVec.norm()/((snapHatApproxMetric[iSnap]).norm()): 0.0;

    avgErr += errorNorm;
    if (errorNorm > maxErr)
      maxErr = errorNorm;   
  }
  avgErr /= (double) nSnaps;
 
  com->fprintf(stderr, " ... Average error on Snapshots after SVD = %e\n", avgErr);  
  com->fprintf(stderr, " ... Maximum error on Snapshots after SVD = %e\n", maxErr);*/
  
  double totalEnergy = 0.0;
  for (int i=0; i <nSnaps; ++i)
    totalEnergy += singValsTmp[i];

  double truncatedEnergy = totalEnergy*approxMetricData->lowRankEnergy;
  double energy = singValsTmp[0];
  int nSnapsRetained = 0;
  while (energy < truncatedEnergy && nSnapsRetained < nSnaps-1) {
    nSnapsRetained++;
    energy += singValsTmp[nSnapsRetained];
  }

  // compute transpose of the pseudo-inverse: US^\dagger V'
  for (int iSnap = 0; iSnap < nSnaps; ++iSnap) {
    for (int jSnap = 0; jSnap < nSnapsRetained; ++jSnap) 
        VtrueTmp[iSnap][jSnap] /= max(singValsTmp[jSnap],1e-16);
  }

  pseudoInverseMaskedSnapsTrans.resize(nSnaps);
  for (int iSnap = 0; iSnap < nSnaps; ++iSnap) {
    pseudoInverseMaskedSnapsTrans[iSnap] = 0.0;
    for (int jSnap = 0; jSnap < nSnapsRetained; ++jSnap)
      pseudoInverseMaskedSnapsTrans[iSnap] += UTrueTmp[jSnap] * VtrueTmp[iSnap][jSnap];
  }

  // check the pseudo inverse

  if (gappyIO->testApproxMetric) {
    double *rVals = new double[nSnaps*(nSnaps+1)/2];
    for (int i = 0; i < nSnaps; i++)
      for (int j = 0; j <= i; j++)
        rVals[(i+1)*i/2 + j] = (snapHatApproxMetric[j]) * (snapHatApproxMetric[i]);

    double temp;

    maxErr = 0.0;
    avgErr = 0.0;
    for (int iSnap = 0; iSnap <nSnaps; ++iSnap) {
      errorVec = snapHatApproxMetric[iSnap];
      for (int jSnap = 0; jSnap < nSnaps; ++jSnap) {
        if (jSnap <= iSnap)
          temp = rVals[(iSnap+1)*iSnap/2 + jSnap];
        else
           temp = rVals[(jSnap+1)*jSnap/2 + iSnap];
        errorVec = errorVec - (temp*pseudoInverseMaskedSnapsTrans[jSnap]);
      }
      errorNorm = errorVec.norm()/((snapHatApproxMetric[iSnap]).norm() + 1e-16);
      avgErr += (errorNorm / (double) nSnaps);
      if (errorNorm > maxErr)
        maxErr = errorNorm;
    }
  
    com->fprintf(stderr, " ... Average error on Snapshots after pseudo-inverse of transpose = %e\n", avgErr);
    com->fprintf(stderr, " ... Maximum error on Snapshots after pseudo-inverse of transpose = %e\n", maxErr);

    delete [] rVals;
  }

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::computeApproximatedMetricLowRankFactor() {

 // compute the SetOfVec lowRankFactor
  if (this->lowRankFactor) delete (this->lowRankFactor); 
  this->lowRankFactor = new VecSet< DistSVec<double, dim> >(numEigen, this->domain.getNodeDistInfo());
  for (int iEigen = 0; iEigen < numEigen; ++iEigen) {
    (*this->lowRankFactor)[iEigen] = 0.0;
    for (int iSnap = 0; iSnap < numSnapsForApproxMetric; ++iSnap)
      (*this->lowRankFactor)[iEigen] += (pseudoInverseMaskedSnapsTrans[iSnap]*lowRankModes[iEigen][iSnap]);
  }

  pseudoInverseMaskedSnapsTrans.resize(0);

  for (int i = 0; i < numEigen; ++i) {
    if (lowRankModes[i]) {
      delete [] lowRankModes[i];
      lowRankModes[i] = NULL;
    }
  }
  if (lowRankModes) {
    delete [] lowRankModes;
    lowRankModes = NULL;
  }
}
//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputApproxMetricLowRankFactorFullCoords(const char* type, int iCluster) {

 // output the SetOfVec lowRankFactor
   com->fprintf(stdout,"\n ... Writing Low Rank Factor of Approximated Metric in Full Mesh Coordinates...\n");

  char *filePath = NULL;
  if (strcmp(type,"state")==0) {
    this->determinePath(this->approxMetricStateLowRankFullCoordsName, -1, filePath);
  } else if (strcmp(type,"nonlinear")==0) {
    this->determinePath(this->approxMetricNonlinearLowRankFullCoordsName, iCluster, filePath);
  } else {
    com->fprintf(stderr,"*** Error: unexpected type in outputApproxMetricLowRankFactorFullCoords\n");
    com->barrier();
    exit(-1);
  }

  for (int iVec = 0; iVec < this->lowRankFactor->numVectors(); ++iVec) {  
     if (lowRankApproxMetricEigenvalues) {
       domain.writeVectorToFile(filePath,iVec,lowRankApproxMetricEigenvalues[iVec],(*this->lowRankFactor)[iVec]);
     } else {
       domain.writeVectorToFile(filePath,iVec,0.0,(*this->lowRankFactor)[iVec]);
     }
  } 
    
  if (filePath) {
    delete [] filePath;
    filePath = NULL;
  }
  delete this->lowRankFactor;
  this->lowRankFactor = NULL;
  if (lowRankApproxMetricEigenvalues) {
    delete [] lowRankApproxMetricEigenvalues;
    lowRankApproxMetricEigenvalues = NULL;
  }

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputApproxMetricLowRankFactorReducedCoords(const char* type, int iCluster) {

 // output the SetOfVec lowRankFactor
  com->fprintf(stdout,"\n ... Writing low rank factor of approximated metric ...\n");

  char *filePath = NULL;
  if (strcmp(type,"state")==0) {
    if (surfaceMeshConstruction) {
      this->determinePath(this->approxMetricStateLowRankSurfaceCoordsName, -1, filePath);
    } else { 
      this->determinePath(this->approxMetricStateLowRankName, -1, filePath);
    }
  } else if (strcmp(type,"nonlinear")==0) {
    this->determinePath(this->approxMetricNonlinearLowRankName, iCluster, filePath);
  } else {
    com->fprintf(stderr,"*** Error: unexpected type in outputApproxMetricLowRankFactorReducedCoords\n");
    com->barrier();
    exit(-1);
  } 

  std::string header("Vector ApproxMetric under load for FluidNodesRed");

  FILE* myOutFile = NULL;
  if (thisCPU==0) {
    myOutFile = fopen(filePath, "wt");
    fprintf(myOutFile,"%s\n", header.c_str());
    fprintf(myOutFile,"%d\n", nReducedNodes);
  }

  int percentComplete = 0;
  for (int iVec = 0; iVec < this->lowRankFactor->numVectors(); ++iVec) {  // # rows in A and B
    outputReducedSVec((*this->lowRankFactor)[iVec],myOutFile,double(iVec));
    if ((iVec+1)%((this->lowRankFactor->numVectors()/4)+1)==0) {
        percentComplete += 25;
        com->fprintf(stdout," ... %3d%% complete ...\n", percentComplete);
    }
  }

  if (thisCPU==0) fclose(myOutFile);

  if (filePath) {
    delete [] filePath;
    filePath = NULL;
  }

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::testInnerProduct(const char *snapshotType) {

  if (strcmp(snapshotType,"approxMetricState")==0) {
    com->fprintf(stdout," \n Computing state approximate metric errors on training data\n");
  } else if (strcmp(snapshotType,"state")==0) {
    com->fprintf(stdout," \n Computing state approximate metric errors on training+validation data\n");
  } else {
    com->fprintf(stderr,"*** Error: unexpected snapshotType in testInnerProduct\n");
    com->barrier();
    exit(-1);
  }

  this->readSnapshotFiles(snapshotType, 0);

  this->readApproxMetricStateLowRankFactor("full");

  std::vector<std::vector<double> > lowRankTimesSnaps;
  lowRankTimesSnaps.resize(numSnapsForApproxMetric);
  for (int i = 0; i < numSnapsForApproxMetric; i++)
    lowRankTimesSnaps[i].resize(this->lowRankFactor->numVectors(), 0.0);

  for (int iSnap = 0; iSnap < numSnapsForApproxMetric; iSnap++) {
    for (int iVec = 0; iVec < this->lowRankFactor->numVectors(); ++iVec) { 
      lowRankTimesSnaps[iSnap][iVec] = ((*this->lowRankFactor)[iVec]) * ((*this->snap)[iSnap]);
    }
  }
 
  double maxRelativeError = 0.0;
  double averageRelativeError = 0.0;
  double errVals = 0.0;
  double approxProd, trueProd;
  for (int i = 0; i < numSnapsForApproxMetric; i++) {
    for (int j = 0; j <= i; j++) {
      approxProd = 0.0;
      for (int iVec = 0; iVec < this->lowRankFactor->numVectors(); ++iVec)
        approxProd += (lowRankTimesSnaps[i][iVec]) *  (lowRankTimesSnaps[j][iVec]);
      trueProd = ((*this->snap)[j]) * ((*this->snap)[i]);
      errVals = abs(trueProd - approxProd) / abs(trueProd);
      averageRelativeError += errVals;
      if (errVals  > maxRelativeError)
        maxRelativeError = errVals;
    }
  }
  averageRelativeError /= (numSnapsForApproxMetric*(numSnapsForApproxMetric+1)/2);
  
  com->fprintf(stdout, "Average relative error for metric = %e\n",averageRelativeError);
  com->fprintf(stdout, "Maximum relative error for metric = %e\n",maxRelativeError);

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::setSampleNodes(int iCluster) {

  globalSampleNodes.clear();

  // for surface mesh construction
  if (surfaceMeshConstruction) {
    if (ioData->romOffline.rob.basisUpdates.preprocessForApproxUpdates==BasisUpdatesData::APPROX_UPDATES_FALSE) { // use a true surface mesh
      com->fprintf(stderr, "*** Note: This surface mesh is not appropriate for approximate basis upates.\n");
      nSampleNodes = 1; // should be zero previously
    } else { // need to include all sampled nodes in surface mesh
      if (int(globalSampleNodesUnion.size())>0) {
        globalSampleNodes = globalSampleNodesUnion;
      } else {
        this->readSampleNodes(-1, "full", false);
        globalSampleNodes = this->sampleNodes;
        this->sampleNodes.clear();
        reinitializeMapsForSampleNodes();
        globalSampleNodesUnion = globalSampleNodes;
      }
      nSampleNodes = globalSampleNodes.size();

      globalSampleNodeRankMap.clear();

      for (int iSampleNode = 0; iSampleNode < nSampleNodes; ++iSampleNode) {
        int globalSampleNode = globalSampleNodes[iSampleNode];
        globalSampleNodeRankMap.insert(pair<int, int > (globalSampleNode, iSampleNode));
      }
    }
  } else {
    if (iCluster<0) {
      if (int(globalSampleNodesUnion.size())>0) {
        globalSampleNodes = globalSampleNodesUnion;
      } else { 
        this->readSampleNodes(iCluster, "full", false);
        globalSampleNodes = this->sampleNodes;
        this->sampleNodes.clear();
        reinitializeMapsForSampleNodes();
        globalSampleNodesUnion = globalSampleNodes;
      }
    } else {
      if (int(globalSampleNodesForCluster.size())>0 && int(globalSampleNodesForCluster[iCluster].size())>0) {
        globalSampleNodes = globalSampleNodesForCluster[iCluster];
      } else {
        this->readSampleNodes(iCluster, "full", false);
        globalSampleNodes = this->sampleNodes;
        this->sampleNodes.clear();
      }
    }

    nSampleNodes = globalSampleNodes.size();

    globalSampleNodeRankMap.clear();

    for (int iSampleNode = 0; iSampleNode < nSampleNodes; ++iSampleNode) {
      int globalSampleNode = globalSampleNodes[iSampleNode];
      globalSampleNodeRankMap.insert(pair<int, int > (globalSampleNode, iSampleNode));
    }
  }
 
  if (!globalNodes && iCluster==unionOfSampleNodes) {
    buildRemainingMesh();
  }

  reducedSampleNodes.clear();
  reducedSampleNodes.resize(nSampleNodes);
  formReducedSampleNodeMap();

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::readGreedyData(int iCluster, bool& breakloop) {

  if (gappyIO->greedyData == GappyConstructionData::UNSPECIFIED_GREEDY  ||
      gappyIO->greedyData == GappyConstructionData::RESIDUAL_ROB_GREEDY ||
      gappyIO->greedyData == GappyConstructionData::JACOBIAN_ROB_GREEDY ||
      gappyIO->greedyData == GappyConstructionData::RESIDUAL_AND_JACOBIAN_ROBS_GREEDY) {
    setUpPodResJac(iCluster);
  } else {

    nPodBasis = 1;
    errorBasis[0] = 0;
    errorBasis[1] = 0;
    
    if (gappyIO->greedyData == GappyConstructionData::STATE_ROB_GREEDY) {
      com->fprintf(stdout, " ... Using state ROB for sample mesh construction ...\n");
      this->readClusteredBasis(iCluster, "state");
      nPod[0] = this->basis->numVectors(); 
      pod[0].resize(nPod[0]);
      for (int iVec=0; iVec<nPod[0]; ++iVec) {
        pod[0][iVec] = (*(this->basis))[iVec];
      }
      delete (this->basis);
      this->basis = NULL;
    }
    else if (gappyIO->greedyData == GappyConstructionData::SPECIFIED_SNAPS_GREEDY) {
      com->fprintf(stdout, " ... Reading specified snapshots for sample mesh construction ...\n");
      this->readSnapshotFiles("greedyData",true);
      nPod[0] = this->snap->numVectors();
      pod[0].resize(nPod[0]);
      for (int iVec=0; iVec<nPod[0]; ++iVec) {
        pod[0][iVec] = (*(this->snap))[iVec];
      }
      delete (this->snap);
      this->snap = NULL; 
      breakloop=true;
    }
 
    nPod[1] = nPod[0]; 
    pod.a[1] = &podRes;
    podHat.a[1] = &podHatRes;
    error.a[1] = &errorRes;  
    errorHat.a[1] = &errorHatRes;
    nPodMax = nPod[0];
     
 }

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::setUpPodResJac(int iCluster) {

  // use one basis if 1) same file name, or 2) jacAction basis unspecified

  if (strcmp(this->jacActionBasisName,this->residualBasisName)==0 ||
      strcmp(this->jacActionBasisName,"")==0) 
    nPodBasis = 1;
  else 
    nPodBasis = 2;
  
  if (gappyIO->greedyData == GappyConstructionData::UNSPECIFIED_GREEDY ||
      gappyIO->greedyData == GappyConstructionData::RESIDUAL_AND_JACOBIAN_ROBS_GREEDY) {
    errorBasis[0] = 0;
    errorBasis[1] = nPodBasis-1;
  }
  else if (gappyIO->greedyData == GappyConstructionData::RESIDUAL_ROB_GREEDY) {
    errorBasis[0] = 0;
    errorBasis[1] = 0;
  }
  else if (gappyIO->greedyData == GappyConstructionData::JACOBIAN_ROB_GREEDY) {
    errorBasis[0] = 1;
    errorBasis[1] = 1;
  }

  com->fprintf(stdout, " ... Reading POD bases for the residual and/or Jacobian ...\n");

  this->readClusteredBasis(iCluster, "residual");
  nPod[0] = this->basis->numVectors(); 
  pod[0].resize(nPod[0]);
  for (int iVec=0; iVec<nPod[0]; ++iVec) {
    pod[0][iVec] = (*(this->basis))[iVec];
  }
  delete (this->basis);
  this->basis = NULL;

  if (nPodBasis == 2) {
    this->readClusteredBasis(iCluster, "jacAction");
    nPod[1] = this->basis->numVectors();
    pod[1].resize(nPod[1]);
    for (int iVec=0; iVec<nPod[1]; ++iVec) {
      pod[1][iVec] = (*(this->basis))[iVec];
    }
    delete (this->basis); 
    this->basis = NULL;   
  } else {
    nPod[1] = nPod[0]; 
    pod.a[1] = &podRes;
    podHat.a[1] = &podHatRes;  // make pod point to res and jac
    error.a[1] = &errorRes;
    errorHat.a[1] = &errorHatRes; 
  }

  nPodMax = max(nPod[0],nPod[1]);

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::setUpBasisBasisProducts() {

// compute pod[1]^Tpod[i] (so you can delete these from memory sooner)

  if (nPodBasis == 2) {
    // pod[1]^T * pod[0]
    computePodTPod();
  }

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::setUpGreedy(int iCluster) {

  //==================================================================
  // PURPOSE: compute number of POD basis vectors (nRhsMax) and globalNodes and handled by each
  // iteration of the greedy algorithm
  // OUTPUTS
  //   nRhsMax, nGreedyIt, nodesToHandle
  // STRATEGY:
  // nothing special happens when nSampleNodes < dimGreedy < nSampleNodes * dim 
  //   1) require dimGreedy < nSampleNodes * dim to avoid underdetermined system
  //   2) if nSampleNodes > dimGreedy, need to treat more globalNodes per iteration
  //   (nRhsMax is 1)
  //==================================================================
  
  if (gappyIO->dimGreedyAlgorithmFactor > 1.0)
    com->fprintf(stderr,"*** Warning: Greedy Factor cannot be larger than 1\n");

  int maxDimGreedyAlgorithm = (gappyIO->maxDimGreedyAlgorithm == -1) ? nPodMax : min(nPodMax,gappyIO->maxDimGreedyAlgorithm);
  dimGreedy = min(maxDimGreedyAlgorithm,
                   max(static_cast<int>(ceil(nPodMax*(gappyIO->dimGreedyAlgorithmFactor))), gappyIO->minDimGreedyAlgorithm));

  double sampleNodeFactor =  gappyIO->sampledNodesFactor;

  if (sampleNodeFactor==-1.0 && gappyIO->maxSampledNodes<=0) {
    com->fprintf(stderr,"*** Error: at least one of sampleNodeFactor and maxSampledNodes must be specified\n");
    com->barrier();
    exit(-1);
  }

  nSampleNodes = 0;
  if (sampleNodeFactor != -1.0) {
    if ((sampleNodeFactor<1.0) &&
             (gappyIO->useUnionOfSampledNodes==GappyConstructionData::UNION_FALSE) &&
             (gappyIO->doPreproGNAT == GappyConstructionData::DO_PREPRO_GNAT_TRUE)) {
      com->fprintf(stderr,"*** Warning: Sampled Node Factor must be greater than or equal to 1\n");
      sampleNodeFactor = 1.0;
      nSampleNodes = static_cast<int>(ceil(double(nPodMax * sampleNodeFactor)/double(dim)));  // this will give interpolation or the smallest possible least squares
    }
  }

  nSampleNodes = (nSampleNodes <= 0) ? gappyIO->maxSampledNodes : min(gappyIO->maxSampledNodes,max(nSampleNodes,gappyIO->minSampledNodes));

  if (nSampleNodes * dim < dimGreedy) {  
    int nSampleNodesOld = nSampleNodes; 
    nSampleNodes = static_cast<int>(ceil(double(dimGreedy)/double(dim))); 
    com->fprintf(stderr,"Warning: not enough sample nodes! Increasing number of sample nodes from %d to %d",nSampleNodesOld,nSampleNodes);
  }

  nRhsMax = static_cast<int>(ceil(double(dimGreedy)/double(nSampleNodes))); // nSampleNodes * nRhsMax >= max(nPod[0],nPod[1])

  // the following should always hold because of the above fix (safeguard)
  
  if (nRhsMax > dim) {
    com->fprintf(stderr,"*** Error: nRhsMax > dim. More nodes should have been added.\n");
    com->barrier();
    exit(-1);
  }

  //==================================================================
  // 2) if nSampleNodes > dimGreedy, need to treat more nodes per iteration
  // strategy: fill more nodes at the earlier iterations because POD basis vectors are optimally ordered
  //==================================================================

  nGreedyIt = min(dimGreedy, nSampleNodes);  // number of greedy iterations (at most dimGreedy; if nSampleNodes > dimGreedy need to take care of more nodes per iteration)
  nodesToHandle = new int[nGreedyIt];  // number of nodes for each greedy iteration

  for (int iGreedyIt = 0; iGreedyIt < nGreedyIt; ++iGreedyIt)  {
    nodesToHandle[iGreedyIt] = (nSampleNodes * nRhsMax) / dimGreedy;
    if (iGreedyIt < nSampleNodes % dimGreedy && nRhsMax ==1)  // only in the dangerous case with nRhsMax = 1
      ++nodesToHandle[iGreedyIt];
  }

  for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis){
    nRhsGreedy[iPodBasis] = new int [nGreedyIt];
    int nRhsMin = min(nPod[iPodBasis],dimGreedy) / nGreedyIt;
    int nRhsExtra = min(nPod[iPodBasis],dimGreedy) % nGreedyIt;
    for (int iGreedyIt = 0; iGreedyIt < nGreedyIt; ++iGreedyIt) {
      nRhsGreedy[iPodBasis][iGreedyIt] = nRhsMin;
      if (iGreedyIt < nRhsExtra) ++nRhsGreedy[iPodBasis][iGreedyIt];
    }
  }

  // initialize sampled pod basis and error vectors

  for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis){
    podHat[iPodBasis].resize(nPod[iPodBasis]);
    for (int i = 0; i < nPod[iPodBasis]; ++i)
      podHat[iPodBasis][i] = 0.0;
    error[iPodBasis].resize(nRhsMax);
    for (int i = 0; i < nRhsMax; ++i) //first iteration, just pick out largest element
      error[iPodBasis][i] = pod[iPodBasis][i];
    errorHat[iPodBasis].resize(nRhsMax);
    for (int i = 0; i < nRhsMax; ++i)
      errorHat[iPodBasis][i] = 0.0;
   }

  //===============================================
  // initialize the least squares problems
  //===============================================
  
  initializeLeastSquares();  // no least squares the first greedy it

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::findMaxAndFillPodHat(const double myMaxNorm, const int
    locSub, const int locNode, const int globalNode) {

  //===============================================
  // PURPOSE: fill locPodHat
  // INPUTS
  //   Local: myMaxNorm, locSub, locNode, globalNode
  // OUTPUTS
  //   Global: locPodHat for maximum entry, globally summed cpuSample, locSubSample,
  //   locNodeSample, globalSampleNodes, xyz
  //===============================================
    
  double globalMaxNorm = myMaxNorm;
  com->barrier();
  com->globalMax(1, &globalMaxNorm);  // find the maximum value over all cpus

  // define variables to be summed globally
  //com->fprintf(stdout," ... max error val = %e\n", globalMaxNorm);

  int cpuTemp = 0;
  int locSubTemp = 0;
  int locNodeTemp = 0;
  int globalNodeTemp = 0;
  double xyz [3];
  for (int i=0; i<3; ++i)
    xyz[i]=0.0;

  // ensure only one cpu enters this loop

  int cpuHasMaxVal = 0;  // indicates if CPU has max value
  int cpuNumWithMaxVal = nTotCpus;
  if (myMaxNorm == globalMaxNorm) {
    cpuHasMaxVal = 1;
    cpuNumWithMaxVal = thisCPU;
  }
  com->globalSum(1, &cpuHasMaxVal);  // total CPUs with maximum value
  if (cpuHasMaxVal > 1) { 
    com->globalMin(1, &cpuNumWithMaxVal);  // take CPU with smallest number that has max val
  }

  //com->fprintf(stdout," ... cpu has max val = %d\n", cpuHasMaxVal);
  if (thisCPU == cpuNumWithMaxVal) {  // if this CPU has the maximum value

    // save the global subdomain and local node indices (sum at the very end of
    // algorithm)

    cpuTemp = thisCPU;
    locSubTemp = locSub;
    locNodeTemp = locNode;
    globalNodeTemp = globalNode;
    assert(locSubTemp!=-1 && locNodeTemp !=-1 && globalNodeTemp != -1);
    computeXYZ(locSub, locNode, xyz);

  }

  // make sure all cpus have the same copy
  
  com->barrier();
  com->globalSum(1, &cpuTemp);
  com->globalSum(1, &locSubTemp);
  com->globalSum(1, &locNodeTemp);
  com->globalSum(1, &globalNodeTemp);
  com->globalSum(3, xyz);

  if (debugging){
   com->fprintf(stdout, "CPU %d has sample node: globalNode = %d, locNode = %d, locSub = %d  \n",cpuTemp,globalNodeTemp,locNodeTemp,locSubTemp);
  }

  // add information to all cpus copies

  cpuSample.push_back(cpuTemp);
  locSubSample.push_back(locSubTemp);
  locNodeSample.push_back(locNodeTemp);
  globalSampleNodes.push_back(globalNodeTemp);

  // define maps for SAMPLE nodes
  globalSampleNodeRankMap.insert(pair<int, int > (globalNodeTemp,
        handledNodes)); globalNodeToCpuMap.insert(pair<int, int >
        (globalNodeTemp, cpuTemp));
  globalNodeToLocSubDomainsMap.insert(pair<int, int > (globalNodeTemp,
        locSubTemp)); globalNodeToLocalNodesMap.insert(pair<int, int >
        (globalNodeTemp, locNodeTemp)); StaticArray<double, 3> XYZ(xyz);
  nodesXYZmap.insert(pair<int, StaticArray <double, 3> > (globalNodeTemp,
        XYZ));

  ++handledNodes;

  if (thisCPU == cpuNumWithMaxVal) {  // if this CPU has the maximum value
    SetOfVec podHatRes_copy(podHatRes.numVectors(), domain.getOldSampledNodeDistInfo());
    SetOfVec podHatJac_copy((nPodBasis==2)?errorHatJac.numVectors():0, domain.getOldSampledNodeDistInfo());
    for (int iPod = 0 ; iPod < podHatRes.numVectors(); ++iPod) { 
      podHatRes_copy[iPod] = podHatRes[iPod];
    }
    if (nPodBasis==2) {
      for (int iPod = 0 ; iPod < podHatJac.numVectors(); ++iPod) {
        podHatJac_copy[iPod] = podHatJac[iPod];
      }
    }

    VecSetArray<dim> podHat_copy; 
    podHat_copy.a[0] = &podHatRes_copy;
    podHat_copy.a[1] = (nPodBasis==2) ? &podHatJac_copy : &podHatRes_copy;

    domain.makeSampledNodeDistInfo(cpuSample, locSubSample);

    podHatRes.resize(podHatRes.numVectors());
    errorHatRes.resize(errorHatRes.numVectors());
    if (nPodBasis==2) podHatJac.resize(podHatJac.numVectors());
    if (nPodBasis==2) errorHatJac.resize(errorHatJac.numVectors());

    SubDomainData<dim> locPod, locPodHat, locPodHat_copy;

    for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {
      for (int iPod = 0 ; iPod < nPod[iPodBasis]; ++iPod) {
        for (int iSub = 0; iSub<numLocSub; ++iSub) {
          int locSampledNode = (iSub==locSub) ? podHat[iPodBasis][iPod].subSize(iSub) - 1 : podHat[iPodBasis][iPod].subSize(iSub);
          locPod = pod[iPodBasis][iPod].subData(iSub); 
          locPodHat = podHat[iPodBasis][iPod].subData(iSub);
          locPodHat_copy = podHat_copy[iPodBasis][iPod].subData(iSub);
          for (int iDim = 0; iDim < dim ; ++iDim) {
            for(int j=0; j<locSampledNode; ++j) {
              locPodHat[j][iDim] = locPodHat_copy[j][iDim];
            }
            if (iSub==locSub) {
              locPodHat[locSampledNode][iDim] = locPod[locNode][iDim];
            }
          }
        }
      }
    }
    domain.makeOldSampledNodeDistInfo(cpuSample, locSubSample);
  }

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::determineSampleNodes() {

  for (int greedyIt = 0; greedyIt < nGreedyIt; ++greedyIt)  {
    com->fprintf(stdout," ... Greedy iteration %d ...\n", greedyIt);
    greedyIteration(greedyIt);
  }

  if (nodesToHandle) {
    delete [] nodesToHandle;
    nodesToHandle = NULL;
  }
  for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {
    if (nRhsGreedy[iPodBasis]) {
      delete [] nRhsGreedy[iPodBasis];
      nRhsGreedy[iPodBasis] = NULL;
    }
  }

  for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis){
    // pod no longer needed
    pod[iPodBasis].resize(0);
    error[iPodBasis].resize(0);
    errorHat[iPodBasis].resize(0);
  }

  if (debugging){
    com->fprintf(stdout,"globalSampleNodes are:");
    for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes)
      com->fprintf(stdout,"%d ",globalSampleNodes[iSampleNodes]);
    com->fprintf(stdout,"\n");
  }

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::greedyIteration(int greedyIt) {

  // Differences for 1st iteration compared with other greedy iterations:
  // 1) no least squares problem is solved (just take the maximum entry)
  // 2) look at the inlet face to ensure boundary condition is handled

  double myMaxNorm;
  int locSub, locNode, globalNode;  // temporary global subdomain and local node

  bool doLeastSquares = true;

  if (greedyIt == 0) {  
    doLeastSquares = false;  // don't do least squares if first iteration
  }

  // determine number of rhs for each
  // typically have nRhs = nRhsMax. exception: there are less than nRhsMax remaining vectors

  for (int iPodBasis = 0; iPodBasis  < nPodBasis; ++iPodBasis)  
    nRhs[iPodBasis] = nRhsGreedy[iPodBasis][greedyIt];
  assert(max(nRhs[0],nRhs[1]) > 0);    // must have at least one RHS

  double totalLSTime = this->timer->getTime();
  if (doLeastSquares) {
    leastSquaresReconstruction();    // solve the least-squares reconstruction
  }
  totalLSTime = this->timer->getTime() - totalLSTime;
  com->fprintf(stdout," ... Total Least Squares Time = %e(s)\n", totalLSTime);

  //double maskErrorTime = this->timer->getTime();
  maskError();
  //maskErrorTime = this->timer->getTime() - maskErrorTime;
  //com->fprintf(stdout," ... Total Masking Time = %e(s)\n", maskErrorTime);


  double totalFillTime = this->timer->getTime();
  for (int iFillNode = 0; iFillNode < nodesToHandle[greedyIt]; ++iFillNode) {  // fill up the appropriate number of nodes

    // initialize parameters
    myMaxNorm = 0.0;  // initial maximum is zero
    locSub = -1; locNode = -1; globalNode = -1;  // where the maximum is located

    // loop over nodes, and add to set if it is the maximum
    // subdomains -> nodes
    double loopTime = this->timer->getTime();
    for (int iSub = 0; iSub < numLocSub; ++iSub) {

      // get subdomain info for all RHS

      getSubDomainError(iSub);

      // find maximum error on the subdomain
      
      subDFindMaxError(iSub, onlyInletOutletBC, myMaxNorm, locSub, locNode, globalNode);
      
    }
    //com->fprintf(stdout," ... ... numLocSub loop = %e(s)\n", this->timer->getTime() - loopTime);
    //com->barrier(); //debugging
    // find global subdomain number and local node number for node with maximum norm
    //double fillTime = this->timer->getTime();
    findMaxAndFillPodHat(myMaxNorm, locSub, locNode, globalNode);
    //com->fprintf(stdout," ... ... fill  = %e(s)\n", this->timer->getTime() - fillTime);

    if (onlyInletOutletBC == true)  // only add one 
      onlyInletOutletBC = false;
  }
  totalFillTime = this->timer->getTime() - totalFillTime;
  com->fprintf(stdout," ... Total Fill Time = %e(s)\n", totalFillTime);


  for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis)  
    handledVectors[iPodBasis] += nRhs[iPodBasis];

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::initializeLeastSquares() {

  // initialize scalapack for least squares problems
  // TODO: only allocate memory for required columns!
  if (gappyIO->greedyLeastSquaresSolver != GappyConstructionData::GREEDY_LS_SCALAPACK) return;
  if (initializeLeastSquaresDone == true) return;

  for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {
    parallelRom[iPodBasis]->parallelLSMultiRHSClean();
    parallelRom[iPodBasis]->parallelLSMultiRHSInit(podHat[iPodBasis], errorHat[iPodBasis], dimGreedy);
  }

  initializeLeastSquaresDone = true;

}


//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::initializeLeastSquaresPseudoInv(int numRhs) {

  // initialize scalapack for least squares problems
  // TODO: only allocate memory for required columns!
  if (gappyIO->pseudoInverseSolver != GappyConstructionData::PSEUDO_INVERSE_SCALAPACK) return;

  for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {
    parallelRom[iPodBasis]->parallelLSMultiRHSClean();
    parallelRom[iPodBasis]->parallelLSMultiRHSInit(podHat[iPodBasis], pseudoInvRhs, numRhs);
  }

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::makeNodeMaxIfUnique(double nodeError, double
    &myMaxNorm, int iSub, int locNodeNum, int &locSub, int &locNode, int
    &globalNode) {

  // PURPOSE: make the node the current maximum on this cpu if it hasn't been added already
  // INPUT
  //   Local: nodeError, myMaxNorm (can change), iSub, locNodeNum (in subdomain node numbering
  // system)
  //   Global: handledNodes, globalSampleNodes,
  // OUTPUT: 
  //   Local: myMaxNorm (can change), locSub, locNode, globalNode
  
  // only do if the node could be the maximum

  if (nodeError >= myMaxNorm) {

    bool newNode = true;
    int *locToGlobNodeMap = subD[iSub]->getNodeMap();
    int thisGlobalNode = locToGlobNodeMap[locNodeNum];

    // check that the node hasn't already been added (look at globalSampleNodes)

    for (int iIslandCheck = 0; iIslandCheck < handledNodes; ++iIslandCheck) {
      if (thisGlobalNode == globalSampleNodes[iIslandCheck]) { 
        newNode = false;
        break;
      }
    }

    // make local maximum if it is maximum and not already in the set

    if (newNode) {
      myMaxNorm = nodeError;
      locSub = iSub;
      locNode = locNodeNum;
      globalNode = thisGlobalNode;
    }
  }
}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::computeNodeError(bool *locMasterFlag, int locNodeNum, double &nodeError) {

  // PURPOSE: compute the sum of squares of node error for all RHS, both bases
  //   at the locNodeNum node
  // INPUT
  //   Local: locMasterFlag, locNodeNum 
  //   Global: nRhs
  // OUTPUT
  //   Global: nodeError

   nodeError = 0.0; // initialize normed error to zero

   if (!locMasterFlag || locMasterFlag[locNodeNum]) {
     for (int iPodBasis = errorBasis[0]; iPodBasis  <= errorBasis[1] ; ++iPodBasis)
       for (int iRhs = 0; iRhs < nRhs[iPodBasis]; ++iRhs){  // add components of error from all vectors on this node
         for (int k = 0; k < dim ; ++k){  // add contributions from residual and jacobian reconstruction errors where possible
           double componentError = locError[iPodBasis][iRhs][locNodeNum][k];
           nodeError += componentError * componentError;
       }
     }
   }
}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::getSubDomainError(int iSub) {

  // INPUT 
  //   Passed: iSub
  //   Global: nPodBasis, error
  // OUTPUT
  //   Global: locError

  locError.resize(nRhs);  // create locError entity of maximal size nRhsMax

  for (int iPodBasis = 0; iPodBasis < nPodBasis ; ++iPodBasis) {
    for (int iRhs = 0; iRhs < nRhs[iPodBasis] ; ++iRhs){  // add components of error from all vectors on this node
      locError[iPodBasis][iRhs] = error[iPodBasis][iRhs].subData(iSub);  // first iteration, it is just the pod vectors themselves
    }
  }
}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::leastSquaresReconstruction() {

  // PURPOSE: compute least squares reconstruction error of the pod basis
  // vectors
  // INPUT
  //   Global: nPodBasis, nRhs, error, podHat, error

  double ** (lsCoeff[2]);
  for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {
    lsCoeff[iPodBasis] = new double * [ nRhs[iPodBasis] ];
    for (int iRhs = 0; iRhs < nRhs[iPodBasis]; ++iRhs)  {
      // temporarily fill the error vector with the RHS (solving nRhs[iPodBasis] problems)
      errorHat[iPodBasis][iRhs] = podHat[iPodBasis][handledVectors[iPodBasis] + iRhs];
      lsCoeff[iPodBasis][iRhs] = new double [handledVectors[iPodBasis]];
    }

    if (gappyIO->greedyLeastSquaresSolver == GappyConstructionData::GREEDY_LS_SCALAPACK) {
      double scalapackTime = this->timer->getTime();
      parallelLSMultiRHSGap(iPodBasis,lsCoeff[iPodBasis]);
      scalapackTime = this->timer->getTime() - scalapackTime;
      com->fprintf(stdout," ... ... Scalapack LS Time = %e(s)\n", scalapackTime);
    } else if (gappyIO->greedyLeastSquaresSolver == GappyConstructionData::GREEDY_LS_LINPACK) {
      double linpackTime = this->timer->getTime();
      serialLSMultiRHSGap(iPodBasis,lsCoeff[iPodBasis]);
      linpackTime = this->timer->getTime() - linpackTime;
      com->fprintf(stdout," ... ... linpack LS Time = %e(s)\n", linpackTime);
    } else {
      //double probTime = this->timer->getTime();
      //probabilisticLSMultiRHSGap(iPodBasis,lsCoeff[iPodBasis]);
      //probTime = this->timer->getTime() - probTime;
      //com->fprintf(stdout," ... ... probabilistic LS Time = %e(s)\n", probTime);
      com->fprintf(stderr," ... ... Probabilistic LS has not been checked with Phil's additions\n");
      sleep(1);
      exit(-1);
    }
 
    // compute reconstruction error

    double reconstructionTime = this->timer->getTime();
    for (int iRhs = 0; iRhs < nRhs[iPodBasis]; ++iRhs) { 
      error[iPodBasis][iRhs] = pod[iPodBasis][handledVectors[iPodBasis] + iRhs];  // NOTE: POD
      for (int jPod = 0; jPod < handledVectors[iPodBasis]; ++jPod) {
        error[iPodBasis][iRhs] -= pod[iPodBasis][jPod] * lsCoeff[iPodBasis][iRhs][jPod];
      }
    }
    for (int iRhs = 0; iRhs < nRhs[iPodBasis]; ++iRhs)  {
      if (lsCoeff[iPodBasis][iRhs]) delete [] lsCoeff[iPodBasis][iRhs];
      lsCoeff[iPodBasis][iRhs] = NULL;
    }
    if (lsCoeff[iPodBasis]) {
      delete [] lsCoeff[iPodBasis];
      lsCoeff[iPodBasis] = NULL;
    }
    reconstructionTime = this->timer->getTime() - reconstructionTime;
    com->fprintf(stdout," ... ... Reconstruction Time = %e(s)\n", reconstructionTime);
  }
}

//----------------------------------------------

template<int dim> void GappyPreprocessing<dim>::subDFindMaxError(int iSub, bool
    onlyInletOutletBC, double &myMaxNorm, int &locSub, int &locNode, int
    &globalNode) {

  // PURPOSE: Search the iSub subdomain for possible maximum error
  // Inputs:
  //   Passed: iSub, onlyInletOutletBC, myMaxNorm, locSub, locNode, globalNode
  // Outputs:
  //   Passed: (all can change) myMaxNorm, locSub, locNode, globalNode

  bool *locMasterFlag = nodeDistInfo.getMasterFlag(iSub); // master nodes on subdomain
  int nLocNodes = nodeDistInfo.subSize(iSub);  // number of nodes in this subdomain
  double nodeError; // the error at the node

  if (onlyInletOutletBC) {  // consider only nodes on inlet BC
    FaceSet& currentFaces = subD[iSub]->getFaces();
    for (int iFace = 0; iFace < subD[iSub]->numFaces(); ++iFace) {  
      // only consider inlet boundary conditions
      if (currentFaces[iFace].getCode() == BC_INLET_MOVING ||
          currentFaces[iFace].getCode() == BC_INLET_FIXED ||
          currentFaces[iFace].getCode() == BC_OUTLET_MOVING ||
          currentFaces[iFace].getCode() == BC_OUTLET_FIXED ||
          currentFaces[iFace].getCode() == BC_DIRECTSTATE_INLET_MOVING ||
          currentFaces[iFace].getCode() == BC_DIRECTSTATE_INLET_FIXED ||
          currentFaces[iFace].getCode() == BC_DIRECTSTATE_OUTLET_MOVING ||
          currentFaces[iFace].getCode() == BC_DIRECTSTATE_OUTLET_FIXED ||
          currentFaces[iFace].getCode() == BC_MASSFLOW_INLET_MOVING ||
          currentFaces[iFace].getCode() == BC_MASSFLOW_INLET_FIXED ||
          currentFaces[iFace].getCode() == BC_MASSFLOW_OUTLET_MOVING ||
          currentFaces[iFace].getCode() == BC_MASSFLOW_OUTLET_FIXED) {
       locMasterFlag = nodeDistInfo.getMasterFlag(iSub); // master nodes on subdomain
       nLocNodes = currentFaces[iFace].numNodes(); // number of nodes on this face
       for (int iLocNode = 0; iLocNode < nLocNodes ; ++iLocNode){
         int locNodeNum = currentFaces[iFace][iLocNode];  // local node number (in subdomain numbering) from the face
         computeNodeError(locMasterFlag, locNodeNum, nodeError);  // compute the error at the node
         // make the node max if it is the maximum so far and has not yet been added
         makeNodeMaxIfUnique(nodeError, myMaxNorm, iSub, locNodeNum, locSub, locNode, globalNode);
       }
      }
    }
  }

  else {  // consider all nodes
    for (int locNodeNum = 0; locNodeNum < nLocNodes ; ++locNodeNum){
      computeNodeError(locMasterFlag, locNodeNum, nodeError);  
      // make the node max if it is the maximum so far and has not yet been added
      makeNodeMaxIfUnique(nodeError, myMaxNorm, iSub, locNodeNum, locSub, locNode, globalNode);
    }
  }
}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::serialLSMultiRHSGap(int iPodBasis, double **lsCoeff) {

  // communicate full vector to all ranks
  int nLHS = handledVectors[iPodBasis];
  int nRHS = nRhs[iPodBasis];
  int sampledDOF = int(globalSampleNodes.size())*dim;
  int arraySizeLHS = nLHS*sampledDOF;
  int arraySizeRHS = nRHS*sampledDOF;

  double* sampledLHS = new double[arraySizeLHS];
  double* sampledRHS = new double[arraySizeRHS];

  for (int iLHS = 0; iLHS < nLHS; ++iLHS) {
    for (int iDOF = 0; iDOF < sampledDOF; ++iDOF) {
      sampledLHS[iLHS*sampledDOF + iDOF] =  0.0;
    }
  }
  for (int iRHS = 0; iRHS < nRHS; ++iRHS) {
    for (int iDOF = 0; iDOF < sampledDOF; ++iDOF) {
      sampledRHS[iRHS*sampledDOF + iDOF] =  0.0;
    }
  }

  for (int iLHS = 0; iLHS < nLHS; ++iLHS) {
    std::vector<int> locSampleNodeCount(numLocSub, 0);
    for (int iNode = 0; iNode < int(globalSampleNodes.size()); ++iNode) {
      int currentGlobalNode = globalSampleNodes[iNode];
      int iCpu = globalNodeToCpuMap.find(currentGlobalNode)->second;
      if (thisCPU == iCpu){
        int iSubDomain = globalNodeToLocSubDomainsMap.find(currentGlobalNode)->second;
        //int iLocalNode = globalNodeToLocalNodesMap.find(currentGlobalNode)->second;
        int locSampleNode = locSampleNodeCount[iSubDomain];
        ++locSampleNodeCount[iSubDomain];
        SubDomainData<dim> locValue = podHat[iPodBasis][iLHS].subData(iSubDomain);
        for (int iDim = 0; iDim < dim; ++iDim) {
          //sampledLHS[iLHS*sampledDOF + iNode*dim + iDim] = locValue[iLocalNode][iDim];
          sampledLHS[iLHS*sampledDOF + iNode*dim + iDim] = locValue[locSampleNode][iDim];
        }
      }
    }
  }
  com->globalSumRoot(arraySizeLHS, sampledLHS);
  if (thisCPU!=0) delete [] sampledLHS;

  for (int iRHS = 0; iRHS < nRHS; ++iRHS) {
    std::vector<int> locSampleNodeCount(numLocSub, 0);
    for (int iNode = 0; iNode < int(globalSampleNodes.size()); ++iNode) {
      int currentGlobalNode = globalSampleNodes[iNode];
      int iCpu = globalNodeToCpuMap.find(currentGlobalNode)->second;
      if (thisCPU == iCpu){
        int iSubDomain = globalNodeToLocSubDomainsMap.find(currentGlobalNode)->second;
        //int iLocalNode = globalNodeToLocalNodesMap.find(currentGlobalNode)->second;
        //SubDomainData<dim> locValue = error[iPodBasis][iRHS].subData(iSubDomain);
        int locSampleNode = locSampleNodeCount[iSubDomain];
        ++locSampleNodeCount[iSubDomain];
        SubDomainData<dim> locValue = errorHat[iPodBasis][iRHS].subData(iSubDomain);
        for (int iDim = 0; iDim < dim; ++iDim) {
          //sampledRHS[iRHS*sampledDOF + iNode*dim + iDim] = locValue[iLocalNode][iDim];
          sampledRHS[iRHS*sampledDOF + iNode*dim + iDim] = locValue[locSampleNode][iDim];
        }
      }
    }
  }
  com->globalSumRoot(arraySizeRHS, sampledRHS);
  if (thisCPU!=0) delete [] sampledRHS;

  // cpu0 calls linpack's SVD routine and then computes least squares coeffs
  if (thisCPU==0) {
    double *yVec;     // allocated in linpackSVD sampledDOF*min(nLHS,sampledDOF)
    double *singVals; // allocated in linpackSVD min(nLHS,sampledDOF)
    double *zVec;     // allocated in linpackSVD nLHS*nLHS
    this->linpackSVD(sampledLHS, sampledDOF, nLHS, yVec, singVals, zVec);
    delete [] sampledLHS;  
    int nVecsU = min(sampledDOF,nLHS);
    std::vector< double > result;
    result.resize(nVecsU, 0.0);
    for (int iRHS=0; iRHS<nRHS; ++iRHS) {
      for (int iVec=0; iVec<nVecsU; ++iVec) {
        result[iVec] = 0.0;
        for (int iDOF=0; iDOF<sampledDOF; ++iDOF) {
          result[iVec] += yVec[iVec*sampledDOF + iDOF] * sampledRHS[iRHS*sampledDOF + iDOF];
        }
        result[iVec] = (singVals[iVec]>1e-15) ? result[iVec] / singVals[iVec] : 0;
      }
      for (int iLHS=0; iLHS<nLHS; ++iLHS) {
        lsCoeff[iRHS][iLHS] = 0.0;
        for (int iVec=0; iVec<nVecsU; ++iVec) {
          lsCoeff[iRHS][iLHS] += zVec[iVec*nLHS + iLHS]*result[iVec];
        }
      }
    }
    delete [] yVec;
    delete [] singVals;
    delete [] zVec;
    delete [] sampledRHS;
  }

  // communicate lsCoeff to all ranks
  for (int iRHS=0; iRHS<nRHS; ++iRHS) {
    com->broadcast(nLHS, lsCoeff[iRHS], 0);
  }
  com->barrier();

}


//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::parallelLSMultiRHSGap(int iPodBasis, double **lsCoeff) {

  bool lsCoeffAllCPU = true; // all cpus need solution
  parallelRom[iPodBasis]->parallelLSMultiRHS(podHat[iPodBasis],errorHat[iPodBasis],
      handledVectors[iPodBasis], nRhs[iPodBasis], lsCoeff, lsCoeffAllCPU);

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::probabilisticLSMultiRHSGap(int iPodBasis, double **lsCoeff) {

  // Solve LS using probabilsitic SVD
  //   System is:  podHat[iPodBasis][1:handledVectors[iPodBasis]] * lsCoeff = error[iPodBasis][1:nRhs[iPodBasis]]
  //   Note that all CPUs need lsCoeff

  // Strategy:
  //   CPU0 generates random martrix, randMat, which is nVecs x k with k <= nVecs
  //   Communicate randMat to all mpi ranks
  //   tmpMat = podHat * randMat
  //   for i = 1:q
  //      tmpMat = podHat * (podHat' * tmpMat)
  //   end
  //   [Q,~] = qr(tmpMat)
  //   tmpMat2 = Q' * podHat
  //   [Utmp S V] = svd(tmpMat2) (in serial using lapack)
  //   Communicate svd result to all mpi ranks
  //   U = Q*Utmp
  //   lsCoeff[iPodBasis] = V * inv(S) * U' * error[iPodBasis]

  int nVecs = handledVectors[iPodBasis];
  bool testSVD=true;
  std::vector<std::vector<double> > lsCoeffVec;

  VecSet< DistSVec<double, dim> > LHS( handledVectors[iPodBasis], this->domain.getNodeDistInfo());
  for (int iVec=0; iVec<handledVectors[iPodBasis]; ++iVec) 
    LHS[iVec] = podHat[iPodBasis][iVec];

  VecSet< DistSVec<double, dim> > RHS( nRhs[iPodBasis], this->domain.getNodeDistInfo());
  for (int iRhs=0; iRhs<nRhs[iPodBasis]; ++iRhs) 
    RHS[iRhs] = error[iPodBasis][iRhs];

  this->probabilisticLSMultiRHS(LHS, RHS, lsCoeffVec, gappyIO->randMatDimension, gappyIO->nPowerIts, testSVD);

  for (int iRhs=0; iRhs<nRhs[iPodBasis]; ++iRhs) {
    for (int iVec=0; iVec<nVecs; ++iVec) {
      lsCoeff[iRhs][iVec] = lsCoeffVec[iRhs][iVec];
    }
  }

  lsCoeffVec.clear();

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::buildRemainingMesh() {

  // globalNodes[iIsland][0] is the sample node itself 
  // globalNodes[iIsland][iNode] is the iNode neighbor of iIsland.


  // quantities for entire sample mesh (not just sample nodes)

  globalNodes = new std::vector <int> [nSampleNodes];
  cpus = new std::vector <int> [nSampleNodes];
  locSubDomains = new std::vector <int> [nSampleNodes];
  localNodes = new std::vector <int> [nSampleNodes];
  for (int i = 0; i < 3; ++i)
    nodesXYZ[i] = new std::vector <double> [nSampleNodes];
  elements = new std::vector <int> [nSampleNodes];
  for (int i = 0; i < 4; ++i)
    elemToNode[i] = new std::vector <int> [nSampleNodes];
  totalNodesCommunicated = new int [nSampleNodes];
  totalEleCommunicated = new int [nSampleNodes];
  for (int i = 0; i < nSampleNodes; ++i) {
    totalNodesCommunicated[i] = 0;
    totalEleCommunicated[i] = 0;
  }

  // compute two layers of globalNodes
  
  nodeOffset = new int [nSampleNodes];
  for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes) 
    nodeOffset[iSampleNodes] = 0;

  com->fprintf(stdout," ... Adding sample nodes and neighbors ...\n");
  addSampleNodesAndNeighbors();  

  com->fprintf(stdout," ... Computing BC faces ...\n");
  computeBCFaces(true);  // compute BC faces and add nodes/elements for lift surfaces

  // add all neighbor globalNodes and elements of the sample node's neighbors
  // and neighbor nodes of the lift surface nodes (need 1 layer of nodes for lift)

  communicateAll();  // all cpus need all nodes/elements for adding neighbors

  if (twoLayers) {// && !surfaceMeshConstruction) {
    com->fprintf(stdout," ... Adding second layer of neighbors ...\n");
    for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes) {
      com->fprintf(stdout," ... Adding neighbors for sample node %d of %d...\n",iSampleNodes+1,nSampleNodes);
      addNeighbors(iSampleNodes, nodeOffset[iSampleNodes]);
    }
    communicateAll();  // all cpus need all nodes/elements to define maps
  }

  com->fprintf(stdout," ... Defining maps ...\n");
  defineMaps();  // define element and node maps
    // deletes cpus, locSubDomains, locNodes, totalNodesCommunicated,
    // totalEleCommunicated, nodesXYZ, elemToNode
    // (want to keep globalNodes and elements: see below)

  // from now on, use maps and not these vectors

  cpuSample.erase(cpuSample.begin(),cpuSample.end());
  locSubSample.erase(locSubSample.begin(),locSubSample.end());
  locNodeSample.erase(locNodeSample.begin(),locNodeSample.end());

  // remove redundant entries from the globalNodes and elements

  domain.makeUnique(globalNodes, nSampleNodes);
  domain.makeUnique(elements, nSampleNodes);

  com->fprintf(stdout," ... Building remaining BC faces ...\n");
  computeBCFaces(false);  // compute BC faces already in mesh
  communicateBCFaces();

  // after domain.makeUnique() globalNodes is now a single vector
  nReducedNodes = globalNodes[0].size();  // number of nodes in the sample mesh

  if (nodeOffset) {
    delete [] nodeOffset; 
    nodeOffset = NULL;
  }
  com->fprintf(stdout," ... Finished building mesh topology ...\n");

} 

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::addSampleNodesAndNeighbors() {

  if (surfaceMeshConstruction && ioData->romOffline.rob.basisUpdates.preprocessForApproxUpdates==BasisUpdatesData::APPROX_UPDATES_FALSE)
    return;
 
  for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes) {
    com->fprintf(stdout," ... Adding neighbors for sample node %d of %d...\n",iSampleNodes,nSampleNodes);

    // add the sample node itself to all processors

    globalNodes[iSampleNodes].push_back(globalSampleNodes[iSampleNodes]);
    cpus[iSampleNodes].push_back(cpuSample[iSampleNodes]);
    locSubDomains[iSampleNodes].push_back(locSubSample[iSampleNodes]);
    localNodes[iSampleNodes].push_back(locNodeSample[iSampleNodes]);

    StaticArray<double, 3> xyzVals = nodesXYZmap.find(globalSampleNodes[iSampleNodes])->second;
    for (int iXYZ = 0; iXYZ < 3; ++iXYZ) {
      nodesXYZ[iXYZ][iSampleNodes].push_back(xyzVals[iXYZ]);
    }

    // add all neighbor globalNodes and elements of the sample node itself
    nodeOffset[iSampleNodes] = globalNodes[iSampleNodes].size();
    addNeighbors(iSampleNodes);
  }

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::addNeighbors(int iIslands, int startingNodeWithNeigh) {

  // add all global neighbor globalNodes/elements in the iIslands row of and elements to the iIsland node set
  
  Connectivity *nodeToNode, *nodeToEle, *eleToNode;
  double xyz [3] = {0.0, 0.0, 0.0};
  int globalNodeNum;  // NOTE: must work with global node numbers because the greedy selection only operated on master globalNodes. Here, we care about the node even if it wasn't a master node.
  int locNodeNum, locEleNum;  // local node number of the node/element to be added
  int nNodesToAdd, nEleToAdd;  // number of globalNodes that should be added
  bool elemBasedConnectivityCreated;  // whether createElemBasedConnectivity has been created for the current subdomain
  int *nToNpointer, *nToEpointer;  // pointers to connectivity information

  int nNodeWithNeigh = globalNodes[iIslands].size();  // number of globalNodes whose neighbors will be added
  bool *locMasterFlag;

  for (int iSub = 0 ; iSub < numLocSub ; ++iSub) {  // all subdomains
    int *locToGlobNodeMap = subD[iSub]->getNodeMap();
    int *locToGlobElemMap = subD[iSub]->getElemMap();
    locMasterFlag = nodeDistInfo.getMasterFlag(iSub); // master nodes on subdomain

    elemBasedConnectivityCreated = false;
    for (int iLocNode = 0; iLocNode < subD[iSub]->numNodes(); ++iLocNode) {  // all local globalNodes in subdomain
      //only add master nodes
      globalNodeNum = locToGlobNodeMap[iLocNode];  // global node number
      for (int iNodeWithNeigh = startingNodeWithNeigh; iNodeWithNeigh < nNodeWithNeigh; ++iNodeWithNeigh) {  // check if this local node is in the current row
        if (globalNodeNum == globalNodes[iIslands][iNodeWithNeigh]) {// add all neighbors of this node
          if (!elemBasedConnectivityCreated) {
            // taken from createElemBasedConnectivity
            nodeToNode = subD[iSub]->createElemBasedConnectivity();
            nodeToEle = subD[iSub]->createNodeToElementConnectivity();
            eleToNode = subD[iSub]->createElementToNodeConnectivity();
            elemBasedConnectivityCreated = true;
          }

          // add neighbors of iLocNode

          nToNpointer = nodeToNode->ptr();
          nToEpointer = nodeToEle->ptr();
          nNodesToAdd = nToNpointer[iLocNode+1] - nToNpointer[iLocNode];  // number of neighbors this node has
          nEleToAdd = nToEpointer[iLocNode+1] - nToEpointer[iLocNode];  // number of neighbors this node has

          for (int iNodeToAdd = 0; iNodeToAdd < nNodesToAdd; ++iNodeToAdd) { 
            locNodeNum = *((*nodeToNode)[iLocNode]+iNodeToAdd);
            bool newNode = true;
            for (std::vector<int>::iterator it = globalNodes[iIslands].begin(); it != globalNodes[iIslands].end(); ++it) {
              if (*it == locToGlobNodeMap[locNodeNum]) {
                newNode=false;
                break;
              }
            }
            if (newNode) {
              globalNodes[iIslands].push_back(locToGlobNodeMap[locNodeNum]);
              cpus[iIslands].push_back(thisCPU);
              locSubDomains[iIslands].push_back(iSub);
              localNodes[iIslands].push_back(locNodeNum);
              for (int iXYZ = 0; iXYZ < 3 ; ++iXYZ) xyz[iXYZ] = 0.0; // KTCREMOVE
              computeXYZ(iSub, locNodeNum, xyz);
              for (int iXYZ = 0; iXYZ < 3; ++iXYZ) {
                nodesXYZ[iXYZ][iIslands].push_back(xyz[iXYZ]);
              }
            }
          }
          for (int iEleToAdd = 0; iEleToAdd < nEleToAdd; ++iEleToAdd) { 
            // add global element
            locEleNum = *((*nodeToEle)[iLocNode]+iEleToAdd);
            bool newElem = true;
            for (std::vector<int>::iterator it = elements[iIslands].begin(); it != elements[iIslands].end(); ++it) {
              if (*it == locToGlobElemMap[locEleNum]) {
                newElem=false;
                break;
              }
            }
            if (newElem) {
              elements[iIslands].push_back(locToGlobElemMap[locEleNum]);
              
              // determine globalNodes connected to the element
              for (int iNodesConn = 0; iNodesConn < 4; ++iNodesConn){
                int locNodeNumTmp = *((*eleToNode)[locEleNum] + iNodesConn);
                int globalNodeNumTmp = locToGlobNodeMap[locNodeNumTmp];
                elemToNode[iNodesConn][iIslands].push_back(globalNodeNumTmp);
              }
            }
          }
        }

      }
    }
  }
}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::computeBCFaces(bool liftContribution) {

  // PURPOSE: determine which faces of the sample mesh are on the boundary
  // INPUT: liftContribution (true: adding faces contributing to lift; false:
  //   adding other BC faces)
  // ASSUME: call with liftContribution = true first!
  // METHOD: loop over ALL FACES
  // NOTE: this is done in parallel; thus, a communication is needed to make
  // sure all cpus have the same copy
    // determine if any faces are boundary conditions

  int faceBCCode = 0;
  bool includeFace = false;
  int globalFaceNodes [3];  // global node numbers of current face
  for (int iSub = 0 ; iSub < numLocSub ; ++iSub) {  // all subdomains
    FaceSet& currentFaces = subD[iSub]->getFaces();  // faces on subdomain
    int *locToGlobNodeMap = subD[iSub]->getNodeMap();  // global node numbering
    for (int iFace = 0; iFace < subD[iSub]->numFaces(); ++iFace) {  // check all faces  
      faceBCCode = currentFaces[iFace].getCode();
      int codeIsPos = faceBCCode >0;
      if ((faceBCCode <= BC_MAX_CODE) && (faceBCCode != 0)) {
        
        // check if the face is in the sample mesh
         if (liftContribution && includeLiftFaces > 0) {  // including lift faces
           includeFace = checkFaceContributesToLift(currentFaces,
               iFace, iSub, locToGlobNodeMap);
           if (includeFace) {
             addFaceNodesElements(currentFaces, iFace, iSub, locToGlobNodeMap);// add nodes
           }
         }
         else if (!liftContribution) {// include faces already in the mesh
           includeFace = checkFaceInMesh(currentFaces, iFace, iSub, locToGlobNodeMap);
           if (includeFace)
             includeFace *= checkFaceAlreadyAdded(thisCPU, iSub, iFace);
         }

        if (includeFace) { // add face to extra set
          for (int iFaceNode = 0; iFaceNode < 3; ++iFaceNode) { 
            int localFaceNode = currentFaces[iFace][iFaceNode];
            globalFaceNodes[iFaceNode] = locToGlobNodeMap[localFaceNode];
          }
          for (int iFaceNode = 0; iFaceNode < 3; ++iFaceNode) {
            bcFaces[codeIsPos][iFaceNode][abs(faceBCCode)].push_back(globalFaceNodes[iFaceNode]);
          }
          int surfaceID = currentFaces[iFace].getSurfaceID();
          bcFaceSurfID[codeIsPos][abs(faceBCCode)].push_back(surfaceID);
          if (liftContribution)  {// only do when considering lift surfaces
            // NOTE: each face is on one subdomain, so no communication of
            // bcFacesInfo is required
            int faceLoc [3]  = {thisCPU, iSub, iFace};
            StaticArray <int,3> faceLocation(faceLoc);
            bcFacesInfo.insert(faceLocation);
          }
        }
      }
    }
  }
}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::communicateAll() {

  domain.communicateMesh(globalNodes, nSampleNodes, totalNodesCommunicated);
  domain.communicateMesh(cpus, nSampleNodes, totalNodesCommunicated);
  domain.communicateMesh(locSubDomains, nSampleNodes, totalNodesCommunicated);
  domain.communicateMesh(localNodes, nSampleNodes, totalNodesCommunicated);
  for (int i = 0; i < 3; ++i)
    domain.communicateMesh(nodesXYZ[i], nSampleNodes, totalNodesCommunicated);
  domain.communicateMesh(elements, nSampleNodes, totalEleCommunicated);
  for (int i = 0; i < 4; ++i)
    domain.communicateMesh(elemToNode[i], nSampleNodes, totalEleCommunicated);
  for (int i = 0; i < nSampleNodes; ++i) { 
    totalNodesCommunicated[i] = globalNodes[i].size();
    totalEleCommunicated[i] = elements[i].size();
  }
}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::defineMaps() {

  // defines nodesXYZmap and elemToNodeMap
  int globalNodeNumTmp;
  StaticArray<double, 3> nodesXYZTmp;
  int globalEleNumTmp;
  StaticArray<int, 4> elemToNodeTmp;
  boost::unordered_map<int,int>::const_iterator sampleNodeRank; 

  // first time, establish that no maps have been defined
  for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes) {
    // define nodesXYZmap    
    for (int iNeighbor = 0; iNeighbor < int(globalNodes[iSampleNodes].size()); ++iNeighbor) {
      // do not re-define map for sample nodes (already defined as master) 
      globalNodeNumTmp = globalNodes[iSampleNodes][iNeighbor];
      sampleNodeRank = globalSampleNodeRankMap.find(globalNodeNumTmp);
      if (sampleNodeRank != globalSampleNodeRankMap.end()) {continue;}  // skip if a sample node

      int cpuTemp = cpus[iSampleNodes][iNeighbor];
      int locSubDomTemp = locSubDomains[iSampleNodes][iNeighbor];
      int localNodesTemp = localNodes[iSampleNodes][iNeighbor];
      for (int iXYZ = 0 ; iXYZ < 3; ++iXYZ)
        nodesXYZTmp[iXYZ] = nodesXYZ[iXYZ][iSampleNodes][iNeighbor];

      globalNodeToCpuMap.insert(pair<int, int > (globalNodeNumTmp, cpuTemp));
      globalNodeToLocSubDomainsMap.insert(pair<int, int > (globalNodeNumTmp, locSubDomTemp));
      globalNodeToLocalNodesMap.insert(pair<int, int > (globalNodeNumTmp, localNodesTemp));
      nodesXYZmap.insert(pair<int, StaticArray <double, 3> > (globalNodeNumTmp, nodesXYZTmp));

    }

    // define elemToNodeMap

    for (int iEle = 0; iEle < int(elements[iSampleNodes].size()); ++iEle) {
      globalEleNumTmp = elements[iSampleNodes][iEle];
      for (int iNodesConn = 0 ; iNodesConn  < 4; ++iNodesConn)
        elemToNodeTmp[iNodesConn] = elemToNode[iNodesConn][iSampleNodes][iEle];
      elemToNodeMap.insert(pair<int, StaticArray <int, 4> > (globalEleNumTmp, elemToNodeTmp));
    }

  }

  // no longer need this information (use in the form of maps)

  if (cpus) {
    delete [] cpus;
    cpus = NULL;
  }
  if (locSubDomains) {
    delete [] locSubDomains;
    locSubDomains = NULL;
  }
  if (localNodes) {
    delete [] localNodes;
    localNodes = NULL;
  }
  if (totalNodesCommunicated) {
    delete [] totalNodesCommunicated;
    totalNodesCommunicated = NULL;
  }
  if (totalEleCommunicated) {
    delete [] totalEleCommunicated;
    totalEleCommunicated = NULL;
  }
  for (int i = 0; i < 3; ++i) {
    if (nodesXYZ[i]) { 
      delete [] nodesXYZ[i];
      nodesXYZ[i] = NULL;
    }
  }
  for (int i = 0; i < 4; ++i) {
    if (elemToNode[i]) {
      delete [] elemToNode[i];
      elemToNode[i] = NULL;
    }
  }
}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::communicateBCFaces(){
  
  int BC_CODE_EXTREME = max(BC_MAX_CODE, -BC_MIN_CODE)+1;

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      domain.communicateMesh(bcFaces[i][j], BC_CODE_EXTREME, NULL);
    }
    domain.communicateMesh(bcFaceSurfID[i], BC_CODE_EXTREME, NULL);
  }
  
}

//----------------------------------------------

template<int dim>
bool GappyPreprocessing<dim>::checkFaceInMesh(FaceSet& currentFaces, const int iFace, const int iSub, const int *locToGlobNodeMap){

  // PURPOSE: determine wheteher or not currentFace is in the sample mesh
  // OUTPUT: faceInMesh = true if the face is in the sample mesh

  bool faceInMesh = false;
  bool faceSomewhereInMesh;
  bool *nodeSomewhereInMesh = new bool [currentFaces[iFace].numNodes()];  // one bool for each face node
  int * globalNodeNum = new int [currentFaces[iFace].numNodes()];
  for (int iNodeFace = 0; iNodeFace < currentFaces[iFace].numNodes(); ++iNodeFace) { // globalNodes on face
    nodeSomewhereInMesh[iNodeFace] = false;
  }
  for (int iNodeFace = 0; iNodeFace < currentFaces[iFace].numNodes(); ++iNodeFace) { // globalNodes on face
    globalNodeNum[iNodeFace] = locToGlobNodeMap[currentFaces[iFace][iNodeFace]];
  }
  for (int iNodeFace = 0; iNodeFace < currentFaces[iFace].numNodes(); ++iNodeFace) { // globalNodes on face
    // check to see if the iNodeFace is in iIsland
    for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes) {
      for (int iReducedMeshNode = 0; iReducedMeshNode < int(globalNodes[iSampleNodes].size()); ++iReducedMeshNode) { // globalNodes on island
        if (globalNodeNum[iNodeFace] == globalNodes[iSampleNodes][iReducedMeshNode]){ 
          nodeSomewhereInMesh[iNodeFace] = true;
          break;
        }
      }
      if (nodeSomewhereInMesh[iNodeFace]) break;
    }
  }

  faceSomewhereInMesh = true;
  for (int iNodeFace = 0; iNodeFace < currentFaces[iFace].numNodes(); ++iNodeFace)
    faceSomewhereInMesh *= nodeSomewhereInMesh[iNodeFace];

  //check if a given element has all the globalNodes

  bool faceInElement = true;
  StaticArray <int, 4> globalNodesInElem;
  if (faceSomewhereInMesh) {  // if the face is in the island
    for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes) {
      for (int iEle = 0; iEle < int(elements[iSampleNodes].size());++iEle) {  // check all elements
        faceInElement = true;
        int globalEleNum = elements[iSampleNodes][iEle];
        globalNodesInElem = elemToNodeMap.find(globalEleNum)->second;
        for (int iNodeFace = 0; iNodeFace < currentFaces[iFace].numNodes(); ++iNodeFace){
          bool nodeInElement = false;
          for (int iNode = 0; iNode < 4; ++iNode) {
            if (globalNodeNum[iNodeFace] == globalNodesInElem[iNode]){
              nodeInElement = true;
              break;
            }
          }
          faceInElement *=nodeInElement;
        }
        if (faceInElement)
          break;
      }
      if (faceInElement)
        break;

    }
  }

  if (faceSomewhereInMesh && faceInElement) {
    faceInMesh = true;
  }
  
  if (globalNodeNum) {
    delete [] globalNodeNum;
    globalNodeNum = NULL;
  }
  if (nodeSomewhereInMesh) {
    delete [] nodeSomewhereInMesh;
    nodeSomewhereInMesh = NULL;
  }

  return faceInMesh;

}

//----------------------------------------------

template<int dim>
bool GappyPreprocessing<dim>::checkFaceAlreadyAdded(const int cpuNum, const int
    iSub, const int iFace){
  
  bool includeFace;
  std::set<StaticArray <int, 3> >::iterator bcFacesInfoIt;  // {iCPU,iSub,iFace}
  int faceLoc [3] = {cpuNum, iSub, iFace};
  StaticArray <int,3> faceLocation(faceLoc);
  bcFacesInfoIt = bcFacesInfo.find(faceLocation);
  if (bcFacesInfoIt == bcFacesInfo.end())  // new face
    includeFace = true;
  else   // already added face
    includeFace = false;
  return includeFace;

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::addFaceNodesElements(FaceSet&
    currentFaces, const int iFace, const int iSub, const int
    *locToGlobNodeMap){

  int *locNodeNums = new int [currentFaces[iFace].numNodes()];
  addNodesOnFace(currentFaces, iFace, iSub, locToGlobNodeMap, locNodeNums);// add nodes
  addElementOfFace(currentFaces, iFace, iSub, locToGlobNodeMap, locNodeNums);
  if (locNodeNums) {
    delete [] locNodeNums;
    locNodeNums = NULL;
  }

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::addNodesOnFace(FaceSet&
    currentFaces, const int iFace, const int iSub, const int
    *locToGlobNodeMap, int *locNodeNums){
  // output: locNodeNums

  int globalNodeNum;
  int localNodeNum;
  double xyz [3] = {0.0, 0.0, 0.0};

  for (int iNodeFace = 0; iNodeFace < currentFaces[iFace].numNodes(); ++iNodeFace) { // globalNodes on face
    localNodeNum = currentFaces[iFace][iNodeFace];
    globalNodeNum = locToGlobNodeMap[localNodeNum];
    for (int iXYZ = 0; iXYZ < 3 ; ++iXYZ) xyz[iXYZ] = 0.0;
    computeXYZ(iSub, localNodeNum, xyz);

    // add to the global set
    globalNodes[0].push_back(globalNodeNum);
    cpus[0].push_back(thisCPU);
    locSubDomains[0].push_back(iSub);
    localNodes[0].push_back(localNodeNum);
    for (int iXYZ = 0; iXYZ < 3 ; ++iXYZ) 
      nodesXYZ[iXYZ][0].push_back(xyz[iXYZ]);
    locNodeNums[iNodeFace] = localNodeNum;
  } 

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::addElementOfFace(FaceSet&
    currentFaces, const int iFace, const int iSub, const int
    *locToGlobNodeMap, const int *locNodeNums){


  bool faceInElement;
  int *locToGlobElemMap = subD[iSub]->getElemMap();

  Connectivity *nodeToEle = subD[iSub]->createNodeToElementConnectivity();
  int *nToEpointer = nodeToEle->ptr();
  int nEleConnected;
  int locEleNum, globEleNum;
  std::set<int> locEleNumsCurrent, locEleNums;
  std::vector<int> intersection;

  // loop on globalNodeNums and add element they all have in common
  for (int iLocNode = 0; iLocNode < currentFaces[iFace].numNodes(); ++iLocNode){
    locEleNumsCurrent.clear();
    int locNodeNum = locNodeNums[iLocNode];
    nEleConnected = nToEpointer[ locNodeNum + 1 ] - nToEpointer[ locNodeNum ];
      //number of elements this node belongs to

    // local element numbers connected to node
    for (int iEle = 0; iEle < nEleConnected; ++iEle) { 
      // add global element
      locEleNumsCurrent.insert(*((*nodeToEle)[locNodeNum]+iEle));
    }

    // which elements the nodes have in common
    if (iLocNode == 0) {
      locEleNums = locEleNumsCurrent;
    }
    else {
      intersection.clear();
      std::set_intersection(locEleNumsCurrent.begin(),locEleNumsCurrent.end(),
          locEleNums.begin(), locEleNums.end(),
          std::back_inserter(intersection));
      locEleNums.clear();
      for (int iIntersect = 0; iIntersect < int(intersection.size()); ++iIntersect)
        locEleNums.insert(intersection[iIntersect]);
    }
  }
  assert(int(locEleNums.size()) == 1);  // must only have one element in common
  // find most common element (shared by all nodes)
  locEleNum = *(locEleNums.begin());
  elements[0].push_back(locToGlobElemMap[locEleNum]);

  // add extra node

  // loop on nodes of element, and add if it is not in the set

  Connectivity *eleToNode;
  eleToNode = subD[iSub]->createElementToNodeConnectivity();
  int locNodeNum;
  for (int iNodesConn = 0; iNodesConn < 4; ++iNodesConn){
    locNodeNum = *((*eleToNode)[locEleNum] + iNodesConn);
    int globalNodeNumTmp = locToGlobNodeMap[locNodeNum];
    elemToNode[iNodesConn][0].push_back(globalNodeNumTmp);
    bool locNodeNumNew = true;
    for (int iLocNode = 0; iLocNode < currentFaces[iFace].numNodes(); ++iLocNode) {
      if (locNodeNum == locNodeNums[iLocNode] ) {
        locNodeNumNew = false;
        break;
      }
    }
    if (locNodeNumNew)  {// must add node if it is new
      double xyz [3] = {0.0, 0.0, 0.0};
      int *locToGlobNodeMap = subD[iSub]->getNodeMap();
      int globalNodeNum = locToGlobNodeMap[locNodeNum];
      computeXYZ(iSub, locNodeNum, xyz);

      globalNodes[0].push_back(globalNodeNum);
      cpus[0].push_back(thisCPU);
      locSubDomains[0].push_back(iSub);
      localNodes[0].push_back(locNodeNum);
      for (int iXYZ = 0; iXYZ < 3 ; ++iXYZ) 
        nodesXYZ[iXYZ][0].push_back(xyz[iXYZ]);
    }
  }

  if (nodeToEle) {
    delete nodeToEle;
    nodeToEle = NULL;
  }

}

//----------------------------------------------

template<int dim>
bool GappyPreprocessing<dim>::checkFaceContributesToLift(FaceSet& faces, const int iFace, const int iSub, const int *locToGlobNodeMap ){

  bool faceContributesToLift;
  map<int, int> surfOutMap = postOp->getSurfMap();
  int idx;
  map<int,int>::iterator it = surfOutMap.find(faces[iFace].getSurfaceID());
  if(it != surfOutMap.end() && it->second != -2)
    idx = it->second;
  else if (includeLiftFaces == 2){  // include face if it is on any moving wall (surfaceMeshConstruction)
    if(faces[iFace].getCode() == BC_ISOTHERMAL_WALL_MOVING ||
        faces[iFace].getCode() == BC_ADIABATIC_WALL_MOVING  ||
        faces[iFace].getCode() == BC_SLIP_WALL_MOVING ||
        faces[iFace].getCode() == BC_POROUS_WALL_MOVING)
      idx = 0;
    else
      idx = -1;
  }
  else
    idx = -1;

  if(idx >= 0) 
    faceContributesToLift = true;
  else
    faceContributesToLift = false;
  return faceContributesToLift; 

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputTopFile(int iCluster) {

  com->fprintf(stdout," ... Writing TOP file ...\n");

  // initialize file

  char *sampledMeshPath = 0;
  if (surfaceMeshConstruction) {
    this->determinePath(this->surfaceMeshName, iCluster, sampledMeshPath);
  } else {
    this->determinePath(this->sampledMeshName, iCluster, sampledMeshPath);
  }
  FILE *reducedMesh = 0;
  if (thisCPU == 0) reducedMesh = fopen(sampledMeshPath, "wt");

  // write out globalNodes

  com->fprintf(reducedMesh,"Nodes FluidNodesRed\n");

  int globalNodeNum, globalEleNum;
  StaticArray <double, 3> xyzVals;
  StaticArray <int, 4> globalNodesTmp, reducedNodes;
  boost::unordered_map<int, int > globalToReducedNodeNumbering;  // one mapping for each island

  reducedSampleNodes.resize(nSampleNodes);

  // save the reduced node number for the sample node
  double refLength = ioData->ref.rv.length;
  for (int j = 0; j < nReducedNodes; ++j) {
    // compute xyz position of the node
    globalNodeNum = globalNodes[0][j];  // global node numbers on sample mesh have been sorted in increasing order
    xyzVals = nodesXYZmap.find(globalNodeNum)->second;
    com->fprintf(reducedMesh, "%d %8.15e %8.15e %8.15e \n", j+1, xyzVals[0]*refLength, xyzVals[1]*refLength, xyzVals[2]*refLength);  

    // associate global node to the current reduced node
    globalToReducedNodeNumbering.insert(pair<int, int> (globalNodeNum, j));
  }

  formReducedSampleNodeMap();

  // write elements

  com->fprintf(reducedMesh,"Elements FluidMeshRed using FluidNodesRed\n");

  for (int iEle = 0; iEle < int(elements[0].size()); ++iEle) {
    globalEleNum = elements[0][iEle];
    globalNodesTmp = elemToNodeMap.find(globalEleNum)->second;
    for (int k = 0; k < 4; ++k) {
      reducedNodes[k] = globalToReducedNodeNumbering.find(globalNodesTmp[k])->second;
    }
    com->fprintf(reducedMesh, "%d 5 %d %d %d %d \n", iEle + 1,
        reducedNodes[0]+1, reducedNodes[1]+1, reducedNodes[2]+1,
        reducedNodes[3]+1);  
  }

  // write out boundary faces

  // from sower user manual
  boundaryConditionsMap.insert(pair<int, std::string > (BC_MASSFLOW_OUTLET_MOVING, "MassFlowOutletMoving"));    // guessing
  boundaryConditionsMap.insert(pair<int, std::string > (BC_MASSFLOW_INLET_MOVING, "MassFlowInletMoving"));     // guessing
  boundaryConditionsMap.insert(pair<int, std::string > (BC_DIRECTSTATE_OUTLET_MOVING, "DirectStateOutletMoving")); // guessing
  boundaryConditionsMap.insert(pair<int, std::string > (BC_DIRECTSTATE_INLET_MOVING, "DirectStateInletMoving"));  // guessing
  boundaryConditionsMap.insert(pair<int, std::string > (BC_POROUS_WALL_MOVING, "PorousMoving"));            // guessing
  boundaryConditionsMap.insert(pair<int, std::string > (BC_OUTLET_MOVING, "OutletMoving"));
  boundaryConditionsMap.insert(pair<int, std::string > (BC_INLET_MOVING, "InletMoving"));
  boundaryConditionsMap.insert(pair<int, std::string > (BC_ADIABATIC_WALL_MOVING, "StickMoving"));
  boundaryConditionsMap.insert(pair<int, std::string > (BC_SLIP_WALL_MOVING, "SlipMoving"));
  boundaryConditionsMap.insert(pair<int, std::string > (BC_ISOTHERMAL_WALL_MOVING, "IsothermalMoving"));  // not sure (not in manual)
  boundaryConditionsMap.insert(pair<int, std::string > (BC_INTERNAL, "Internal"));
  boundaryConditionsMap.insert(pair<int, std::string > (BC_ISOTHERMAL_WALL_FIXED, "IsothermalFixed"));  // not sure (not in manual)
  boundaryConditionsMap.insert(pair<int, std::string > (BC_SLIP_WALL_FIXED, "SlipFixed"));
  boundaryConditionsMap.insert(pair<int, std::string > (BC_ADIABATIC_WALL_FIXED, "StickFixed"));
  boundaryConditionsMap.insert(pair<int, std::string > (BC_INLET_FIXED, "InletFixed"));
  boundaryConditionsMap.insert(pair<int, std::string > (BC_OUTLET_FIXED, "OutletFixed"));
  boundaryConditionsMap.insert(pair<int, std::string > (BC_POROUS_WALL_FIXED, "PorousFixed"));            // guessing
  boundaryConditionsMap.insert(pair<int, std::string > (BC_DIRECTSTATE_INLET_FIXED, "DirectStateInletFixed"));  // guessing
  boundaryConditionsMap.insert(pair<int, std::string > (BC_DIRECTSTATE_OUTLET_FIXED, "DirectStateOutletFixed")); // guessing
  boundaryConditionsMap.insert(pair<int, std::string > (BC_MASSFLOW_INLET_FIXED, "MassFlowInletFixed"));     // guessing
  boundaryConditionsMap.insert(pair<int, std::string > (BC_MASSFLOW_OUTLET_FIXED, "MassFlowOutletFixed"));    // guessing
  boundaryConditionsMap.insert(pair<int, std::string > (BC_SYMMETRY, "Symmetry"));

  if (int(boundaryConditionsMap.size()) != (BC_MAX_CODE-BC_MIN_CODE+1)) {
    this->com->fprintf(stderr,"*** Error: new boundary conditions have been added -- the GNAT preprocessing code needs to be updated\n");
    com->barrier();
    exit(-1);  
  }

  for (int iSign = 0; iSign < 2; ++iSign) {
    for (int iBCtype = 0; iBCtype <= max(BC_MAX_CODE, -BC_MIN_CODE); ++iBCtype) {
      int maxID = 0;
      if (int(bcFaceSurfID[iSign][iBCtype].size()) > 0) {
        std::vector<int>::iterator maxIDit = max_element(bcFaceSurfID[iSign][iBCtype].begin(),bcFaceSurfID[iSign][iBCtype].end());
        maxID = *maxIDit;
        if (maxID == 0) {
          this->com->fprintf(stderr,"*** Error: The GNAT preprocessing code apparently assumes that mesh surfaces are labeled with unique, non-zero numerical identifiers (for example: StickMoving_1, StickMoving_2)\n\n\n");
          com->barrier();
          exit(-1);
        }
      }
      
      for (int iID = 0; iID < maxID; ++iID) {
        bool firstTime = true;
        int faceCounter = 0;
        for (int iFace = 0; iFace < int(bcFaces[iSign][0][iBCtype].size()) ; ++iFace) {
          int currentID = bcFaceSurfID[iSign][iBCtype][iFace];
          if (currentID == iID + 1) {  // face is in the current set
            if (firstTime) {  // only output the first time (and only if size > 0)
              int boundaryCondNumber = iBCtype * ((iSign > 0)*2 - 1) ;  // returns boundary condition number
              std::string boundaryCond = boundaryConditionsMap.find(boundaryCondNumber)->second;
              // note: multiplying bcFaceSurfID by 100 so you can visualize it in xpost along-side the full mesh (error if surfID is repeated)
              com->fprintf(reducedMesh,"Elements %sSurface_%d using FluidNodesRed\n",boundaryCond.c_str(), bcFaceSurfID[iSign][iBCtype][iFace]*100);
              firstTime = false;
            }
            for (int k = 0; k < 3; ++k) {
              int globalNode = bcFaces[iSign][k][iBCtype][iFace];
              reducedNodes[k] = globalToReducedNodeNumbering[globalNode];
            }
            com->fprintf(reducedMesh, "%d 4 %d %d %d \n", ++faceCounter, reducedNodes[0]+1, reducedNodes[1]+1, reducedNodes[2]+1);  
          }
        }
      }
    }
  }

  if (thisCPU == 0) fclose(reducedMesh);
  if (sampledMeshPath) {
    delete [] sampledMeshPath;
    sampledMeshPath = NULL;
  }

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (bcFaces[i][j]) {
        delete [] bcFaces[i][j];
        bcFaces[i][j] = NULL;
      }
    }
    if (bcFaceSurfID[i]) {
      delete [] bcFaceSurfID[i];
      bcFaceSurfID[i] = NULL;
    }
  }
  if (elements) {
    delete [] elements;
    elements = NULL;
  }

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputSampleNodes(int iCluster) {

  // output sample nodes (local mask) in sample mesh node numbering system
  char *sampledNodesPath = NULL;
  this->determinePath(this->sampledNodesName, iCluster, sampledNodesPath);   // memory for name is deallocated in outputSampleNodesGeneral
  com->fprintf(stdout," ... Writing sample node file with respect to sample mesh...\n");
  outputSampleNodesGeneral(reducedSampleNodes, sampledNodesPath);

  // output sample nodes (local mask) in full mesh node numbering system
  char *sampledNodesFullCoordsPath = NULL;
  this->determinePath(this->sampledNodesFullCoordsName, iCluster, sampledNodesFullCoordsPath);  // memory for name is deallocated in outputSampleNodesGeneral
  com->fprintf(stdout," ... Writing sample node file with respect to full mesh...\n");
  outputSampleNodesGeneral(globalSampleNodes, sampledNodesFullCoordsPath);

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputSampleNodesGeneral(const std::vector<int> &sampleNodes, const char *outSampleNodeFile) {

  FILE *writingFile;
  if (thisCPU ==0) writingFile = fopen(outSampleNodeFile, "wt");

  com->fprintf(writingFile, "%d", sampleNodes.size());  // first print number of sample nodes
  for (int i = 0; i < int(sampleNodes.size()); ++i) {
    com->fprintf(writingFile, "\n");
    com->fprintf(writingFile, "%d %d", i+1, sampleNodes[i]+1);  
  }

  delete [] outSampleNodeFile;
  outSampleNodeFile = NULL;
  if (thisCPU == 0) fclose(writingFile);

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::computeXYZ(int iSub, int iLocNode, double *xyz) {

  // input: iSub, iLocNode
  // output: x,y,z coordinates of the iLocNode located on iSub
  DistSVec<double,3>& X = geoState->getXn();
  SVec<double,3>& Xsub = X(iSub);  // X is of type DistSVec<double,3>
  for (int i = 0; i < 3; ++i) xyz[i] = Xsub[iLocNode][i];  

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::serialPseudoInverse(int iPodBasis) {
// just call linpack's SVD routine 
// inputs: vecset podHat[iPodBasis]
// outputs: double** podHatPseudoInv[iPodBasis] ([iDof<nSampledDOF][iPod<nPod[iPodBasis]])
 
  // parallelRom[iPodBasis]->parallelLSMultiRHS(podHat[iPodBasis], pseudoInvRhs, nPod[iPodBasis], numRhs, podHatPseudoInvTmp, false);
  int nVecs = nPod[iPodBasis];
  int sampledDOF = nSampleNodes*dim;
  int arraySize = nVecs*sampledDOF;
  double* podHatArray = new double[arraySize];

  for (int iVec = 0; iVec < nVecs; ++iVec) {
    for (int iDOF = 0; iDOF < sampledDOF; ++iDOF) {
      podHatArray[iVec*sampledDOF + iDOF] =  0.0;
    }
  }

  for (int iVec = 0; iVec < nVecs; ++iVec) {
    std::vector<int> locSampleNodeCount(numLocSub, 0);
    for (int iNode = 0; iNode < nSampleNodes; ++iNode) {
      int currentGlobalNode = globalSampleNodes[iNode];
      int iCpu = globalNodeToCpuMap.find(currentGlobalNode)->second;
      if (thisCPU == iCpu){
        int iSubDomain = globalNodeToLocSubDomainsMap.find(currentGlobalNode)->second;
        //int iLocalNode = globalNodeToLocalNodesMap.find(currentGlobalNode)->second;
        int locSampleNode = locSampleNodeCount[iSubDomain];
        ++locSampleNodeCount[iSubDomain];
        SubDomainData<dim> locValue = podHat[iPodBasis][iVec].subData(iSubDomain);
        for (int iDim = 0; iDim < dim; ++iDim) {
          //podHatArray[iVec*sampledDOF + iNode*dim + iDim] = locValue[iLocalNode][iDim];
          podHatArray[iVec*sampledDOF + iNode*dim + iDim] = locValue[locSampleNode][iDim];
        }
      }
    }
  }
  com->globalSumRoot(arraySize, podHatArray);
  if (thisCPU!=0) delete [] podHatArray;

  // cpu0 calls linpack's SVD routine and then computes least squares coeffs
  if (thisCPU==0) {

    double *yVec;     // allocated in linpackSVD sampledDOF*min(nLHS,sampledDOF)
    double *singVals; // allocated in linpackSVD min(nLHS,sampledDOF)
    double *zVec;     // allocated in linpackSVD nLHS*nLHS
    this->linpackSVD(podHatArray, sampledDOF, nVecs, yVec, singVals, zVec);
    delete [] podHatArray;  
    int nVecsU = min(sampledDOF,nVecs);
    assert(nVecsU==nVecs);

    //podHatPseudoInv[iPodBasis][iDOF][iVec]
    podHatPseudoInv[iPodBasis] = new double * [sampledDOF] ;
    for (int iDOF = 0; iDOF < sampledDOF; ++iDOF)  
      podHatPseudoInv[iPodBasis][iDOF] = new double [nVecs] ;

    double** tmpPodHatPseudoInv = new double * [sampledDOF] ; // temporary storage for computing podHatPseudoInv
    for (int iDOF = 0; iDOF < sampledDOF; ++iDOF)
      tmpPodHatPseudoInv[iDOF] = new double [nVecs] ;

    for (int iDOF=0; iDOF<sampledDOF; ++iDOF) {
      for (int iVec=0; iVec<nVecsU; ++iVec) {
        tmpPodHatPseudoInv[iDOF][iVec] = (singVals[iVec]>1e-15) ? yVec[iVec*sampledDOF+iDOF] / singVals[iVec] : 0;
      }
    }

    for (int iDOF=0; iDOF<sampledDOF; ++iDOF) {
      for (int iVec=0; iVec<nVecs; ++iVec) {
        podHatPseudoInv[iPodBasis][iDOF][iVec]=0.0;
        for (int jVec=0; jVec<nVecsU; ++jVec) {
          podHatPseudoInv[iPodBasis][iDOF][iVec] += zVec[jVec*nVecs + iVec]*tmpPodHatPseudoInv[iDOF][jVec];
        }
      }
    }

    delete [] yVec;
    delete [] singVals;
    delete [] zVec;
    for (int iVec = 0; iVec < nVecs; ++iVec) delete [] tmpPodHatPseudoInv[iVec];
    delete [] tmpPodHatPseudoInv; 
  }

  com->barrier();

  if (thisCPU == 0 && iPodBasis == 0) {
    podHatPseudoInv[1] = podHatPseudoInv[0];
  }

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::parallelPseudoInverse(int iPodBasis) {

//======================================
// Purpose
//   compute pseudo inverses of podHatRes or podHatJac
// Inputs
//   iPodBasis (0 or 1)
// Outputs
//   the pseudo inverse
// Approach
//   the ith column of the pseudo inverse solves min || podHat xi = ei ||_2
//   where ei is zero except 1 at the node and dim of interest
//======================================

  // generate pseudoInvRhs in chunks

  int nNodesAtATime = gappyIO->pseudoInverseNodes;  
  bool lastTime = false;
  nNodesAtATime = min(nSampleNodes,nNodesAtATime);  // fix if needed
  int numRhs = nNodesAtATime * dim;

  int nHandledVectors = 0;

  pseudoInvRhs.resize(numRhs);  // make correct size (possible memory problem!)
  for (int iRhs = 0; iRhs < numRhs; ++iRhs)
    pseudoInvRhs[iRhs] = 0.0;

  if (thisCPU == 0) {
    podHatPseudoInv[iPodBasis] = new double * [nSampleNodes * dim] ;
    for (int iRhs = 0; iRhs < nSampleNodes * dim; ++iRhs)  
      podHatPseudoInv[iPodBasis][iRhs] = new double [nPod[iPodBasis] ] ;
  }

  initializeLeastSquaresPseudoInv(numRhs);

  // all processors store the temporary solution
  double **podHatPseudoInvTmp;

  int iVector = 0;
  std::vector<int> locSampleNodeCount(numLocSub, 0);
  for (int iSampleNodes = 0; iSampleNodes < nSampleNodes; ++iSampleNodes) {
    // compute unit vectors
    int currentGlobalNode = globalSampleNodes[iSampleNodes];
    int currentCPU = globalNodeToCpuMap.find(currentGlobalNode)->second;
    if (thisCPU == currentCPU) {
      int iSubDomain = globalNodeToLocSubDomainsMap.find(currentGlobalNode)->second;
      int iLocalNode = locSampleNodeCount[iSubDomain];
      locSampleNodeCount[iSubDomain]++;
      //int iLocalNode = globalNodeToLocalNodesMap.find(currentGlobalNode)->second;
      for (int iDim = 0; iDim < dim; ++iDim) {
        SubDomainData<dim> locValue = pseudoInvRhs[iVector+iDim].subData(iSubDomain);
        locValue[iLocalNode][iDim] = 1.0;
      }
    }
    iVector += dim;
    com->barrier();  // temporary debugging
    // compute part of pseudo-inverse

    if (thisCPU == 0)  // directly store computed value on cpu 0
      podHatPseudoInvTmp = podHatPseudoInv[iPodBasis]+nHandledVectors;

    if (((iSampleNodes + 1)  % nNodesAtATime == 0 && !lastTime) || iSampleNodes == nSampleNodes - 1) {
      com->fprintf(stdout," computing pseudo inverse at sample node %d of %d \n", iSampleNodes+1, nSampleNodes);
      parallelRom[iPodBasis]->parallelLSMultiRHS(podHat[iPodBasis], pseudoInvRhs,
          nPod[iPodBasis], numRhs, podHatPseudoInvTmp, false);

      nHandledVectors+=numRhs;
      if (nNodesAtATime > nSampleNodes - iSampleNodes - 1) {
        nNodesAtATime = nSampleNodes - iSampleNodes - 1;
        lastTime = true;
      }
      
      int numRhsOld = numRhs;
      numRhs = nNodesAtATime * dim;
      if (numRhs != numRhsOld)
        pseudoInvRhs.resize(numRhs);
      for (int iRhs = 0; iRhs < numRhs; ++iRhs)
        pseudoInvRhs[iRhs] = 0.0;
      iVector = 0;  // re-filling rhs
    }
  }

  if (thisCPU == 0 && iPodBasis == 0) {
    podHatPseudoInv[1] = podHatPseudoInv[0];
  }
}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::computePseudoInverse() {

  for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {
    if (gappyIO->pseudoInverseSolver == GappyConstructionData::PSEUDO_INVERSE_SCALAPACK) {
      parallelPseudoInverse(iPodBasis);
    } else if (gappyIO->pseudoInverseSolver == GappyConstructionData::PSEUDO_INVERSE_LINPACK) {
      serialPseudoInverse(iPodBasis);
    }
  }
  //checkConsistency();  // debugging check

  // podHat no longer needed

  for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis)
    podHat[iPodBasis].resize(0);

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::computePodTPod() {

  // podTpod is nPod[1] x nPod[0] array
  //
  // since podTpod = I for nPodBasis = 1, only call this function for
  // nPodBasis = 2

  if (thisCPU == 0) {
    podTpod = new double * [nPod[1]];
    for (int i = 0; i < nPod[1]; ++i)
      podTpod[i] = new double [nPod[0]];
  }

  double podTpodTmp;
  for (int i = 0; i < nPod[1]; ++i) {
    com->fprintf(stdout," ... computing podJac^TpodRes row %d of %d ...\n",i,nPod[1]);
    for (int j = 0; j < nPod[0]; ++j) {
      podTpodTmp = pod[1][i]*pod[0][j];
      if (thisCPU == 0)
        podTpod[i][j] = podTpodTmp;
    }
  }
}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::assembleOnlineMatrices() {

  // Purpose: assemble matrices that are used online
  // Inputs: podHatPseudoInv, podTpod
  // Outputs: onlineMatrices

  int numCols = nSampleNodes * dim;

  if (nPodBasis == 1) { 
    onlineMatrices[0] = podHatPseudoInv[0];  // related to the jacobian
  }
  else {  // nPod[1] != 0 because nPodBasis == 2
    int nPodJac = (nPod[1] == 0) ? nPod[0] : nPod[1];
    onlineMatrices[1] = podHatPseudoInv[1];  // related to the jacobian

    onlineMatrices[0] = new double * [numCols] ;
    for (int i = 0; i < numCols; ++i) {
      onlineMatrices[0][i] = new double [nPod[1] ] ;
      for (int iPod = 0; iPod < nPod[1]; ++iPod) 
        onlineMatrices[0][i][iPod] = 0.0;
    }

    if (thisCPU == 0) {
      for (int iPod = 0; iPod < nPod[1]; ++iPod) {
        for (int jPod = 0; jPod < nPod[0]; ++jPod) {
          for (int k = 0; k < numCols; ++k) { 
            onlineMatrices[0][k][iPod] += podTpod[iPod][jPod] * podHatPseudoInv[0][k][jPod];
          }
        }
      }
    }

  }

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputOnlineMatrices(int iCluster) {

  // output matrices A and B in ASCII form as VecSet< DistSVec> with the
  // DistSVec defined on the sample mesh. Each column in this VecSet
  // corresponds to a row of A or B.

  outputReducedToFullNodes();

  com->fprintf(stdout," ... Writing online matrices ...\n");

  // sample mesh

  if (outputOnlineMatricesSample) 
  outputOnlineMatricesGeneral(iCluster, nReducedNodes, reducedSampleNodeRankMap, reducedSampleNodes);
  
  if (outputOnlineMatricesFull)  // for multi-step preprocessing where the full mesh *is* the sampled mesh 
    outputOnlineMatricesGeneral(iCluster, numFullNodes, globalSampleNodeRankMap, globalSampleNodes);

  for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {
    if (podHatPseudoInv[iPodBasis]) {
      // onlineMatrices might just be a pointer to podHatPseudoInv
      if (onlineMatrices[iPodBasis] == podHatPseudoInv[iPodBasis]) {
        onlineMatrices[iPodBasis] = NULL;
      } else {
        for (int iRhs = 0; iRhs < nSampleNodes * dim; ++iRhs) {
          if (onlineMatrices[iPodBasis][iRhs]) {
            delete [] onlineMatrices[iPodBasis][iRhs];
            onlineMatrices[iPodBasis][iRhs] = NULL;
          }
        }
        delete [] onlineMatrices[iPodBasis];
        onlineMatrices[iPodBasis] = NULL;
      }
      for (int iRhs = 0; iRhs < nSampleNodes * dim; ++iRhs) {
        if (podHatPseudoInv[iPodBasis][iRhs]) {
          delete [] podHatPseudoInv[iPodBasis][iRhs];
          podHatPseudoInv[iPodBasis][iRhs] = NULL;
        }
      }
      delete [] podHatPseudoInv[iPodBasis];
      podHatPseudoInv[iPodBasis] = NULL;
    }
  }

 /* for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {
    if (podHatPseudoInv[iPodBasis]) {
      for (int iRhs = 0; iRhs < nSampleNodes * dim; ++iRhs) {
        if (podHatPseudoInv[iPodBasis][iRhs]) {
          delete [] podHatPseudoInv[iPodBasis][iRhs];
          podHatPseudoInv[iPodBasis][iRhs] = NULL;
        }
      }
      delete [] podHatPseudoInv[iPodBasis];
      podHatPseudoInv[iPodBasis] = NULL;
    }
    if (ffWeight!=1 || (nPodBasis==2 && iPodBasis==0)) {
      if (onlineMatrices[iPodBasis]) {
        for (int i = 0; i < nSampleNodes * dim; ++i) {
          if (onlineMatrices[iPodBasis][i]) {
            delete [] onlineMatrices[iPodBasis][i];
            onlineMatrices[iPodBasis][i] = NULL;
          }
        }
        delete [] onlineMatrices[iPodBasis];
        onlineMatrices[iPodBasis] = NULL;
      }
    } else {
      onlineMatrices[iPodBasis] = NULL;
    }
  }*/

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputOnlineMatricesGeneral(int iCluster, int numNodes,
    const boost::unordered_map<int,int> &sampleNodeMap, const std::vector<int>
    &sampleNodeVec) {

  double sampledOutputTime = this->timer->getTime();  

  if (thisCPU ==0){
  
    // create isSampleNodeVec
    std::vector<bool> isSampleNodeVec(numNodes,false);
    std::vector<int> sampleNodeRankVec(numNodes, -1);
    for(boost::unordered_map<int,int>::const_iterator it = sampleNodeMap.begin(); it != sampleNodeMap.end(); ++it) {
      sampleNodeRankVec[it->first] = it->second;
      isSampleNodeVec[it->first] = true;
    }
  
    // determine file names
    char* gappyResPath = NULL;
    char* gappyJacPath = NULL;
    this->determinePath(this->gappyResidualName, iCluster, gappyResPath);
    this->determinePath(this->gappyJacActionName, iCluster, gappyJacPath);
  
    int nPodJac = (nPod[1] == 0) ? nPod[0] : nPod[1];
  
    for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {  
     
      char* outFilePath = (iPodBasis==0) ? gappyResPath : gappyJacPath;

      // clean up any temporary files that may be hanging around from a crashed simulation
      if (cleanTempFiles) {
        FILE *shell;
        std::string rmCommandString("rm ");
        rmCommandString += outFilePath;
        rmCommandString += "_step_* ";
        const char *rmCommandChar = rmCommandString.c_str();
   
        //fprintf(stdout, "\n%s\n", rmCommandChar);
        if (!(shell = popen(rmCommandChar, "r"))) {
          com->fprintf(stderr, " *** Error: attempt to use system call (rm) failed!\n");
          sleep(1);
          exit(-1);
        }
   
        char buff[512];
        while (fgets(buff, sizeof(buff), shell)!=NULL){}
        pclose(shell);
      }
    

      FILE* myOutFile = fopen(outFilePath, "wt");
      fprintf(myOutFile,"Vector abMatrix under load for FluidNodes\n"); // header
      fprintf(myOutFile,"%d\n",numNodes);
 
      for (int iPod = 0; iPod < nPodJac; ++iPod) {  // # rows in A and B

        fprintf(myOutFile,"%d\n", iPod); // tag
  
        for (int iNode = 0; iNode < numNodes; ++iNode) {
  
          bool isSampleNode = isSampleNodeVec[iNode];
          int sampleNodeRank = sampleNodeRankVec[iNode];
  
          if (dim==5) {
            if (isSampleNode) {
                fprintf(myOutFile,"%8.15e %8.15e %8.15e %8.15e %8.15e\n",
                onlineMatrices[iPodBasis][sampleNodeRank*dim][iPod],
                onlineMatrices[iPodBasis][sampleNodeRank*dim + 1][iPod],
                onlineMatrices[iPodBasis][sampleNodeRank*dim + 2][iPod],
                onlineMatrices[iPodBasis][sampleNodeRank*dim + 3][iPod],
                onlineMatrices[iPodBasis][sampleNodeRank*dim + 4][iPod]);
            }
            else // must output zeros if it is not a sample node
                fprintf(myOutFile,"%8.15e %8.15e %8.15e %8.15e %8.15e\n",0.0,0.0,0.0,0.0,0.0); 
          } else if (dim==6) {
            if (isSampleNode) {
                fprintf(myOutFile,"%8.15e %8.15e %8.15e %8.15e %8.15e %8.15e\n",
                onlineMatrices[iPodBasis][sampleNodeRank*dim][iPod],
                onlineMatrices[iPodBasis][sampleNodeRank*dim + 1][iPod],
                onlineMatrices[iPodBasis][sampleNodeRank*dim + 2][iPod],
                onlineMatrices[iPodBasis][sampleNodeRank*dim + 3][iPod],
                onlineMatrices[iPodBasis][sampleNodeRank*dim + 4][iPod],
                onlineMatrices[iPodBasis][sampleNodeRank*dim + 5][iPod]);
            }
            else // must output zeros if it is not a sample node
                fprintf(myOutFile,"%8.15e %8.15e %8.15e %8.15e %8.15e %8.15e\n",0.0,0.0,0.0,0.0,0.0,0.0);  
          } else if (dim==7) {
            if (isSampleNode) {
                fprintf(myOutFile,"%8.15e %8.15e %8.15e %8.15e %8.15e %8.15e %8.15e\n",
                onlineMatrices[iPodBasis][sampleNodeRank*dim][iPod],
                onlineMatrices[iPodBasis][sampleNodeRank*dim + 1][iPod],
                onlineMatrices[iPodBasis][sampleNodeRank*dim + 2][iPod],
                onlineMatrices[iPodBasis][sampleNodeRank*dim + 3][iPod],
                onlineMatrices[iPodBasis][sampleNodeRank*dim + 4][iPod],
                onlineMatrices[iPodBasis][sampleNodeRank*dim + 5][iPod],
                onlineMatrices[iPodBasis][sampleNodeRank*dim + 6][iPod]);
            }
            else // must output zeros if it is not a sample node
                fprintf(myOutFile,"%8.15e %8.15e %8.15e %8.15e %8.15e %8.15e %8.15e\n",0.0,0.0,0.0,0.0,0.0,0.0,0.0);
          } else { // general, but slow
            for (int iDim = 0; iDim < dim; ++iDim) { 
              if (isSampleNode) {
                fprintf(myOutFile,"%8.15e ",
                onlineMatrices[iPodBasis][sampleNodeRank*dim+iDim][iPod]);
              }
              else // must output zeros if it is not a sample node
                fprintf(myOutFile,"%8.15e ", 0.0);
            }
            fprintf(myOutFile,"\n");
          }
        } // end iNode loop (rows)
  
      } // end iPod loop (columns)
  
      fclose(myOutFile);

    } // end iPodBasis
  
    
    if (gappyResPath) {
       delete [] gappyResPath;
       gappyResPath = NULL;
    }
    if (gappyJacPath) {
       delete [] gappyJacPath;
       gappyJacPath = NULL;
    }

  }

  com->barrier(); 
  this->timer->addSampledOutputTime(sampledOutputTime);
}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputLocalReferenceStateReduced(int iCluster) {

  com->fprintf(stdout,"\n ... Writing reference state for local state basis %d in sample mesh coordinates ...\n", iCluster);

  char *sampledRefStatePath = NULL;
  if (surfaceMeshConstruction) {
    this->determinePath(this->surfaceRefStateName, iCluster, sampledRefStatePath);
  } else {
    this->determinePath(this->sampledRefStateName, iCluster, sampledRefStatePath);
  }

  std::string header("Vector ReferenceState under load for FluidNodesRed");

  // read in full reference state
  this->readClusteredReferenceState(iCluster,"state");

  FILE* myOutFile = NULL;
  if (thisCPU==0) {
    myOutFile = fopen(sampledRefStatePath, "wt");
    fprintf(myOutFile,"%s\n", header.c_str());
    fprintf(myOutFile,"%d\n", nReducedNodes);
  }

  // output
  outputReducedSVec(*(this->Uref),myOutFile,this->Uref->norm());

  if (thisCPU==0) fclose(myOutFile);

  if (this->Uref) {
    delete this->Uref;
    this->Uref=NULL;
  }

  if (sampledRefStatePath) {
    delete [] sampledRefStatePath;
    sampledRefStatePath = NULL;
  }

}
//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputMatchStateReduced() {

  if (strcmp(ioData->input.matchStateFile,"")!=0) {

    com->fprintf(stdout," ... Writing comparison state in sample mesh coordinates ...\n");
  
    char *sampledSolutionPath = NULL;
    if (surfaceMeshConstruction) {
      this->determinePath(this->surfaceMatchStateName, -1, sampledSolutionPath);
    } else {
      this->determinePath(this->sampledMatchStateName, -1, sampledSolutionPath);
    }
  
    // read solution

    char *solFile = new char[strlen(ioData->input.prefix) + strlen(ioData->input.matchStateFile) + 1];
    sprintf(solFile, "%s%s", ioData->input.prefix, ioData->input.matchStateFile);
    
    DistSVec<double,dim> *matchState = new DistSVec<double,dim>( domain.getNodeDistInfo() );
    double tmp;
    bool status = domain.readVectorFromFile(solFile, 0, &tmp, *matchState);

    if (!status) {
      com->fprintf(stderr, "*** Error: unable to read vector from file %s\n", solFile);
      com->barrier();
      exit(-1);
    }

    delete [] solFile;

    // output
    std::string header("Vector MatchState under load for FluidNodesRed");
    FILE* myOutFile = NULL;
    if (thisCPU==0) {
      myOutFile = fopen(sampledSolutionPath, "wt");
      fprintf(myOutFile,"%s\n", header.c_str());
      fprintf(myOutFile,"%d\n", nReducedNodes);
    }

    outputReducedSVec(*matchState,myOutFile,matchState->norm());

    if (thisCPU==0) fclose(myOutFile);

    if (matchState) {
      delete matchState;
      matchState = NULL;
    }

    if (sampledSolutionPath) {
      delete [] sampledSolutionPath;
      sampledSolutionPath = NULL;
    }

  }
}

//----------------------------------------------


template<int dim>
void GappyPreprocessing<dim>::outputInitialConditionReduced() {
  //INPUTS
  // ioData, nReducedNodes, domain
  // needed by outputReducedSVec: globalNodes, globalNodeToCpuMap,
  // globalNodeToLocSubDomainsMap, globalNodeToLocalNodesMap 

  if (strcmp(ioData->input.solutions,"")!=0) {

    com->fprintf(stdout," ... Writing initial condition in sample mesh coordinates ...\n");
  
    char *sampledSolutionPath = NULL;
    if (surfaceMeshConstruction) {
      this->determinePath(this->surfaceSolutionName, -1, sampledSolutionPath);
    } else {
      this->determinePath(this->sampledSolutionName, -1, sampledSolutionPath);
    }
  

    // read in initial condition

    char *icFile = new char[strlen(ioData->input.prefix) + strlen(ioData->input.solutions) + 1];
    sprintf(icFile, "%s%s", ioData->input.prefix, ioData->input.solutions);
    
    DistSVec<double,dim> *initialCondition = new DistSVec<double,dim>( domain.getNodeDistInfo() );
    double tmp;
    bool status = domain.readVectorFromFile(icFile, 0, &tmp, *initialCondition);

    if (!status) {
      com->fprintf(stderr, "*** Error: unable to read vector from file %s\n", icFile);
      com->barrier();
      exit(-1);
    }

    delete [] icFile;

    // output
    std::string header("Vector InitialCondition under load for FluidNodesRed");
    FILE* myOutFile = NULL;
    if (thisCPU==0) {
      myOutFile = fopen(sampledSolutionPath, "wt");
      fprintf(myOutFile,"%s\n", header.c_str());
      fprintf(myOutFile,"%d\n", nReducedNodes);
    }

    outputReducedSVec(*initialCondition,myOutFile,initialCondition->norm());

    if (thisCPU==0) fclose(myOutFile);


    if (initialCondition) {
      delete initialCondition;
      initialCondition = NULL;
    }

    if (sampledSolutionPath) {
      delete [] sampledSolutionPath;
      sampledSolutionPath = NULL;
    }

  }
}
//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputApproxSnapsReduced(int iCluster) {
  //INPUTS
  // nReducedNodes, snapHat
  // needed by outputReducedSVec: globalNodes, globalNodeToCpuMap,
  // globalNodeToLocSubDomainsMap, globalNodeToLocalNodesMap 

  com->fprintf(stdout," ... Writing approximate snapshots in sample mesh coordinates ...\n");


  // create mask
  DistVec<double> mask(domain.getNodeDistInfo());
  mask = 0.0;
  for (int iSampleNode=0; iSampleNode<nSampleNodes; ++iSampleNode) {
    int currentSampleNode = globalSampleNodes[iSampleNode];
    int locSub = globalNodeToLocSubDomainsMap.find(currentSampleNode)->second;
    int locNode = globalNodeToLocalNodesMap.find(currentSampleNode)->second;
    int currentCpu = globalNodeToCpuMap.find(currentSampleNode)->second;
    if (thisCPU == currentCpu) {
      double* locMask = mask.subData(locSub);
      locMask[locNode] = 1.0;
    }
  }
  com->barrier();
  com->fprintf(stdout," ... norm of mask squared = %e \n", mask*mask);

  char *sampledApproxSnapsPath = NULL;
  this->determinePath(this->sampledApproxMetricNonlinearSnapsName, iCluster, sampledApproxSnapsPath);

  std::string header("Vector SampledApproxMertricSnaps under load for FluidNodesRed");

  // output
  FILE* myOutFile = NULL;
  if (thisCPU==0) {
    myOutFile = fopen(sampledApproxSnapsPath, "wt");
    fprintf(myOutFile,"%s\n", header.c_str());
    fprintf(myOutFile,"%d\n", nReducedNodes);
  }

  DistSVec<double,dim> maskedVec(domain.getNodeDistInfo());
  for (int iSnap=0; iSnap<this->snap->numVectors(); ++iSnap) {  // # rows in A and B
    maskedVec = (*this->snap)[iSnap];
    maskedVec *= mask;
    outputReducedSVec(maskedVec,myOutFile,(*this->snap)[iSnap].norm());
  }

  if (thisCPU==0) fclose(myOutFile);
  if (sampledApproxSnapsPath) {
    delete [] sampledApproxSnapsPath;
    sampledApproxSnapsPath = NULL;
  }

}


//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputClusterCentersReduced() {
  //INPUTS
  // nReducedNodes
  // needed by outputReducedSVec: globalNodes, globalNodeToCpuMap,
  // globalNodeToLocSubDomainsMap, globalNodeToLocalNodesMap 

  com->fprintf(stdout," ... Writing cluster centers in sample mesh coordinates ...\n");

  if (!(this->clusterCenters)) this->readClusterCenters("centers");

  char *sampledCentersPath = NULL;
  if (surfaceMeshConstruction) {
    this->determinePath(this->surfaceCentersName, -1, sampledCentersPath);
  } else {
    this->determinePath(this->sampledCentersName, -1, sampledCentersPath);
  }

  std::string header("Vector ClusterCenters under load for FluidNodesRed");
  FILE* myOutFile = NULL;
  if (thisCPU==0) {
    myOutFile = fopen(sampledCentersPath, "wt");
    fprintf(myOutFile,"%s\n", header.c_str());
    fprintf(myOutFile,"%d\n", nReducedNodes);
  }

  // output
  for (int iCluster=0; iCluster < this->nClusters; ++iCluster) {  // # rows in A and B
    outputReducedSVec((*this->clusterCenters)[iCluster],myOutFile,(*this->clusterCenters)[iCluster].norm());
  }

  if (thisCPU==0) fclose(myOutFile);
  if (sampledCentersPath) {
    delete [] sampledCentersPath;
    sampledCentersPath = NULL;
  }

}


//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputLocalStateBasisReduced(int iCluster) {
  //INPUTS
  // ioData, nReducedNodes, domain
  // needed by outputReducedSVec: globalNodes, globalNodeToCpuMap,
  // globalNodeToLocSubDomainsMap, globalNodeToLocalNodesMap 


  if (gappyIO->outputReducedBases) {
  
    com->fprintf(stdout, " ... Reading local state ROB ...\n");
    this->readClusteredBasis(iCluster, "state");
 
    if (surfaceMeshConstruction) { 
      com->fprintf(stdout," ... Writing local state ROB in surface mesh coordinates ...\n");
    } else {
      com->fprintf(stdout," ... Writing local state ROB in sample mesh coordinates ...\n");
    }  

    double t0 = this->timer->getTime(); 

    char *filePath = NULL;
    if (surfaceMeshConstruction) {
      this->determinePath(this->surfaceStateBasisName, iCluster, filePath);
    } else {
      this->determinePath(this->sampledStateBasisName, iCluster, filePath);
    }
  
    std::string header("Vector PodState under load for FluidNodesRed");
    FILE* myOutFile = NULL;
    if (thisCPU==0) {
      myOutFile = fopen(filePath, "wt");
      fprintf(myOutFile,"%s\n", header.c_str());
      fprintf(myOutFile,"%d\n", nReducedNodes);
    }

    int percentComplete = 0;
    for (int iPod = 0; iPod < this->basis->numVectors(); ++iPod) {  // # rows in A and B
      outputReducedSVec((*(this->basis))[iPod],myOutFile,double(iPod));
      if ((iPod+1)%((this->basis->numVectors()/4)+1)==0) {
          percentComplete += 25;
          com->fprintf(stdout," ... %3d%% complete ...\n", percentComplete);
      }
    }
  
    if (thisCPU==0) fclose(myOutFile);
    if (filePath) {
      delete [] filePath;
      filePath = NULL;
    }
  
    com->fprintf(stdout, " ... time to write basis (as recorded by CPU 0) = %e\n", this->timer->getTime()-t0);

  } else {
    com->fprintf(stdout," ... Skipping output of reduced bases...\n");
  }

  if (this->basis) {
    delete (this->basis);
    (this->basis) = NULL;
  }

  if (this->sVals) {
    delete (this->sVals);
    (this->sVals) = NULL;
  }

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputMultiSolutionsReduced() {
  //INPUTS
  // nReducedNodes
  // needed by outputReducedSVec: globalNodes, globalNodeToCpuMap,
  // globalNodeToLocSubDomainsMap, globalNodeToLocalNodesMap 

  if (ioData->input.multiSolutions[0] == 0) return;

  com->fprintf(stdout," ... Writing multiple solutions in sample mesh coordinates ...\n");

  // open file for reading
  FILE *inFP = fopen(input->multiSolutions,"r");
  if (!inFP)  {
    com->fprintf(stderr, "*** Error: No solution data FILES in %s\n", input->multiSolutions);
    com->barrier();
    exit (-1);
  }
  int nData, _n;
  _n = fscanf(inFP, "%d",&nData);


  char *sampledMultiSolutionsPath = NULL;
  if (surfaceMeshConstruction) {
    this->determinePath(this->surfaceMultiSolutionsName, -1, sampledMultiSolutionsPath);
  } else {
    this->determinePath(this->sampledMultiSolutionsName, -1, sampledMultiSolutionsPath);
  }

  char solnFile[500];
  DistSVec<double,dim> solVec(domain.getNodeDistInfo());
  std::string header("Vector MultiSoutions under load for FluidNodesRed");

  FILE* myOutFile = NULL;
  if (thisCPU==0) {
    myOutFile = fopen(sampledMultiSolutionsPath, "wt");
    fprintf(myOutFile,"%s\n", header.c_str());
    fprintf(myOutFile,"%d\n", nReducedNodes);
  }

  if (this->multiUic) {
    for (int iIC=0; iIC < this->multiUic->numVectors(); ++iIC) {
    outputReducedSVec((*this->multiUic)[iIC],myOutFile,double(iIC));
    }
  } else {
    for (int iData=0; iData < nData; ++iData) {
      _n = fscanf(inFP, "%s", solnFile);
      domain.readVectorFromFile(solnFile, 0, 0, solVec);
      outputReducedSVec(solVec,myOutFile,double(iData));
    }
  }

  if (thisCPU==0) fclose(myOutFile);
  if (sampledMultiSolutionsPath) {
    delete [] sampledMultiSolutionsPath;
    sampledMultiSolutionsPath = NULL;
  }

}



//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputReducedSVec(const DistSVec<double,dim> &distSVec, FILE* myOutFile, double tag) {

  double sampledOutputTime = this->timer->getTime();  

  // communicate full vector to CPU0
  double* sampledVec;
  int sampledDOF = nReducedNodes*dim;
  sampledVec = new double[sampledDOF];
  for (int iDOF = 0; iDOF < sampledDOF; ++iDOF) {
    sampledVec[iDOF] =  0.0;
  }

  for (int iNode = 0; iNode < nReducedNodes; ++iNode) {
    int currentGlobalNode = globalNodes[0][iNode];
    int iCpu = globalNodeToCpuMap.find(currentGlobalNode)->second;
    if (thisCPU == iCpu){
      int iSubDomain = globalNodeToLocSubDomainsMap.find(currentGlobalNode)->second;
      int iLocalNode = globalNodeToLocalNodesMap.find(currentGlobalNode)->second;
      SubDomainData<dim> locValue = distSVec.subData(iSubDomain);
      for (int iDim = 0; iDim < dim; ++iDim) {
        sampledVec[iNode*dim + iDim] = locValue[iLocalNode][iDim];
      }
    }
  }

  com->globalSumRoot(sampledDOF, sampledVec); 
  if (thisCPU != 0) delete [] sampledVec;

  if (thisCPU == 0) { 
    fprintf(myOutFile,"%e\n", tag);
    for (int iNode = 0; iNode < nReducedNodes; ++iNode) {
      if (dim==5) {
        fprintf(myOutFile,"%8.15e %8.15e %8.15e %8.15e %8.15e\n",
            sampledVec[iNode*dim],
            sampledVec[iNode*dim + 1],
            sampledVec[iNode*dim + 2],
            sampledVec[iNode*dim + 3],
            sampledVec[iNode*dim + 4]);
      } else if (dim==6) {
        fprintf(myOutFile,"%8.15e %8.15e %8.15e %8.15e %8.15e %8.15e\n",
            sampledVec[iNode*dim],
            sampledVec[iNode*dim + 1],
            sampledVec[iNode*dim + 2],
            sampledVec[iNode*dim + 3],
            sampledVec[iNode*dim + 4],
            sampledVec[iNode*dim + 5]);
      } else if (dim==7) {
        fprintf(myOutFile,"%8.15e %8.15e %8.15e %8.15e %8.15e %8.15e %8.15e\n",
            sampledVec[iNode*dim],
            sampledVec[iNode*dim + 1],
            sampledVec[iNode*dim + 2],
            sampledVec[iNode*dim + 3],
            sampledVec[iNode*dim + 4],
            sampledVec[iNode*dim + 5],
            sampledVec[iNode*dim + 6]);
      } else { // general, but slow
        for (int iDim = 0; iDim < dim; ++iDim) {
          fprintf(myOutFile,"%8.15e ", sampledVec[iNode*dim + iDim]);
        }
        fprintf(myOutFile,"\n");
      }
    }
    delete [] sampledVec;
  }

  com->barrier();
  if (surfaceMeshConstruction) {
    this->timer->addSurfaceOutputTime(sampledOutputTime); 
  } else {     
    this->timer->addSampledOutputTime(sampledOutputTime);
  }
 
}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputReduced3DSVec(const DistSVec<double,3> &distSVec, FILE* myOutFile, double tag) {

  double sampledOutputTime = this->timer->getTime();

  // communicate full vector to CPU0
  double* sampledVec;
  int sampledDOF = nReducedNodes*3;
  sampledVec = new double[sampledDOF];
  for (int iDOF = 0; iDOF < sampledDOF; ++iDOF) {
    sampledVec[iDOF] =  0.0;
  }

  for (int iNode = 0; iNode < nReducedNodes; ++iNode) {
    int currentGlobalNode = globalNodes[0][iNode];
    int iCpu = globalNodeToCpuMap.find(currentGlobalNode)->second;
    if (thisCPU == iCpu){
      int iSubDomain = globalNodeToLocSubDomainsMap.find(currentGlobalNode)->second;
      int iLocalNode = globalNodeToLocalNodesMap.find(currentGlobalNode)->second;
      SubDomainData<3> locValue = distSVec.subData(iSubDomain);
      for (int iDim = 0; iDim < 3; ++iDim) {
        sampledVec[iNode*3 + iDim] = locValue[iLocalNode][iDim];
      }
    }
  }

  com->globalSumRoot(sampledDOF, sampledVec); 
  if (thisCPU != 0) delete [] sampledVec;

  if (thisCPU == 0) { 
    fprintf(myOutFile,"%e\n", tag);
    for (int iNode = 0; iNode < nReducedNodes; ++iNode) {
        fprintf(myOutFile,"%8.15e %8.15e %8.15e\n",
            sampledVec[iNode*3],
            sampledVec[iNode*3 + 1],
            sampledVec[iNode*3 + 2]);
    }
    delete [] sampledVec;
  }

  com->barrier();
  if (surfaceMeshConstruction) {
    this->timer->addSurfaceOutputTime(sampledOutputTime); 
  } else {     
    this->timer->addSampledOutputTime(sampledOutputTime);
  }
 
}
//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputWallDistanceReduced() {
  // INPUTS: geoState, ioData, nReducedNodes
  // needed by outputReducedVec: globalNodes, globalNodeToCpuMap, globalNodeToLocSubDomainsMap, globalNodeToLocalNodesMap

  com->fprintf(stdout," ... Writing wall distance for sample mesh ...\n");

  // load in wall distance

  DistVec<double> *d2wall = geoState->getd2wall();
  DistVec<double> d2wallOutput(d2wall->info());
  d2wallOutput = -ioData->bc.wall.delta;
    // must subtract off for input files (see DistGeoState.C constructor)
  d2wallOutput += *d2wall;

  char *wallDistPath = NULL;
  if (surfaceMeshConstruction) {
    this->determinePath(this->surfaceWallDistName, -1, wallDistPath);  //a bit silly, but easier than the alternative
  } else {
    this->determinePath(this->sampledWallDistName, -1, wallDistPath);  
  }


  FILE *outWallDist;
  if (thisCPU ==0) outWallDist = fopen(wallDistPath, "wt");

  if (thisCPU == 0) fprintf(outWallDist,"Scalar walldist under load for FluidNodesRed\n");
  if (thisCPU == 0) fprintf(outWallDist,"%d\n", nReducedNodes);
  outputReducedVec(d2wallOutput,outWallDist,0);

  if (wallDistPath) {
    delete [] wallDistPath;
    wallDistPath = NULL;
  }
  if (thisCPU == 0) fclose(outWallDist);

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputShapeDerivativeReduced() {


  if (strcmp(ioData->input.shapederivatives,"")!=0) {

    com->fprintf(stdout," ... Writing shape derivatives in sample mesh coordinates ...\n");

    char *sampledShapeDerivativePath = NULL;
    if (surfaceMeshConstruction) {
      this->determinePath(this->surfaceShapeDerivativeName, -1, sampledShapeDerivativePath);
    } else {
      this->determinePath(this->sampledShapeDerivativeName, -1, sampledShapeDerivativePath);
    }
    char *derFile = new char[strlen(ioData->input.prefix) + strlen(ioData->input.shapederivatives) + 1];
    sprintf(derFile, "%s%s", ioData->input.prefix, ioData->input.shapederivatives);

    DistSVec<double,3> *shapeDerivative = new DistSVec<double,3>( domain.getNodeDistInfo() );
    double tmp = 0.0;
    bool status;
    int nShapeDer=-1;
    while (true) {
      nShapeDer += 1;
      status = domain.readVectorFromFile(derFile, nShapeDer, &tmp, *shapeDerivative);

      if (!status && nShapeDer == 0) {
        com->fprintf(stderr, "*** Error: unable to read vector from file %s\n", derFile);
        exit(-1);
      } else if (!status && nShapeDer > 0) {
        break;
      }
    }

    std::string header("Vector ShapeDerivative under load for FluidNodesRed");
    FILE* myOutFile = NULL;
    if (thisCPU==0) {
      myOutFile = fopen(sampledShapeDerivativePath, "wt");
      fprintf(myOutFile,"%s\n", header.c_str());
      fprintf(myOutFile,"%d\n", nReducedNodes);
    }

    int shapeDerNum = -1;
    while (true) {
      shapeDerNum += 1;
      status = domain.readVectorFromFile(derFile, shapeDerNum, &tmp, *shapeDerivative);

      if (!status && shapeDerNum == 0) {
        com->fprintf(stderr, "*** Error: unable to read vector from file %s\n", derFile);
        exit(-1);
      } else if (!status && shapeDerNum > 0) {
        break;
      }
      //outputReducedSVec(*shapeDerivative,sampledShapeDerivativePath,double(shapeDerNum),shapeDerNum,nShapeDer,header.c_str());
      outputReduced3DSVec(*shapeDerivative,myOutFile,double(shapeDerNum));
    }

    if (thisCPU==0) fclose(myOutFile);
    if (shapeDerivative) {
      delete shapeDerivative;
      shapeDerivative = NULL;
    }

    if (sampledShapeDerivativePath) {
      delete [] sampledShapeDerivativePath;
      sampledShapeDerivativePath = NULL;
    }

    delete [] derFile;
  }
}


//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputDisplacementReduced() {
  //INPUTS
  // ioData, nReducedNodes, domain
  // needed by outputReducedSVec: globalNodes, globalNodeToCpuMap,
  // globalNodeToLocSubDomainsMap, globalNodeToLocalNodesMap 

  if (strcmp(ioData->input.displacements,"")!=0) {

    com->fprintf(stdout," ... Writing volume deformation vector in sample mesh coordinates ...\n");
  
    char *sampledDisplacementPath = NULL;
    if (surfaceMeshConstruction) {
      this->determinePath(this->surfaceDisplacementName, -1, sampledDisplacementPath);
    } else {
      this->determinePath(this->sampledDisplacementName, -1, sampledDisplacementPath);
    }
  
    // read in initial condition

    char *dispFile = new char[strlen(ioData->input.prefix) + strlen(ioData->input.displacements) + 1];
    sprintf(dispFile, "%s%s", ioData->input.prefix, ioData->input.displacements);
    
    DistSVec<double,3> *displacement = new DistSVec<double,3>(domain.getNodeDistInfo());

    double tmp;
    bool status = domain.readVectorFromFile(dispFile, 0, &tmp, *displacement);

    if (!status) {
      com->fprintf(stderr, "*** Error: unable to read vector from file %s\n", dispFile);
      com->barrier();
      exit(-1);
    }

    delete [] dispFile;

    // output
    std::string header("Vector InitialDisplacement under load for FluidNodesRed");
    FILE* myOutFile = NULL;
    if (thisCPU==0) {
      myOutFile = fopen(sampledDisplacementPath, "wt");
      fprintf(myOutFile,"%s\n", header.c_str());
      fprintf(myOutFile,"%d\n", nReducedNodes);
    }

    outputReduced3DSVec(*displacement,myOutFile,displacement->norm());

    if (thisCPU==0) fclose(myOutFile);

    if (displacement) {
      delete displacement;
      displacement = NULL;
    }

    if (sampledDisplacementPath) {
      delete [] sampledDisplacementPath;
      sampledDisplacementPath = NULL;
    }

  }

}




//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputReducedVec(const DistVec<double> &distVec, FILE* outFile , int iVector) {

  com->fprintf(outFile,"%d\n", iVector);

  if (gappyIO->useOldReducedSVecFunction) {

    for (int j = 0; j < nReducedNodes; ++j) {
      double value = 0.0; // initialize value to zero
      int currentGlobalNode = globalNodes[0][j];
      int iCpu = globalNodeToCpuMap.find(currentGlobalNode)->second;
  
      if (thisCPU == iCpu){
        int iSubDomain = globalNodeToLocSubDomainsMap.find(currentGlobalNode)->second;
        int iLocalNode = globalNodeToLocalNodesMap.find(currentGlobalNode)->second;
        double *locValue = distVec.subData(iSubDomain);
        value = locValue[iLocalNode];
      }
  
      com->globalSum(1, &value);
      com->fprintf(outFile,"%8.15e ", value);
      com->fprintf(outFile,"\n");
    }
  
  } else {

    double* sampledVec;
    sampledVec = new double[nReducedNodes];
    for (int iDOF = 0; iDOF < nReducedNodes; ++iDOF) {
      sampledVec[iDOF] =  0.0;
    }

    for (int iNode = 0; iNode < nReducedNodes; ++iNode) {
      int currentGlobalNode = globalNodes[0][iNode];
      int iCpu = globalNodeToCpuMap.find(currentGlobalNode)->second;
      if (thisCPU == iCpu){
        int iSubDomain = globalNodeToLocSubDomainsMap.find(currentGlobalNode)->second;
        int iLocalNode = globalNodeToLocalNodesMap.find(currentGlobalNode)->second;
        double* locValue = distVec.subData(iSubDomain);
        sampledVec[iNode] = locValue[iLocalNode];
      }
    }

    com->globalSumRoot(nReducedNodes, sampledVec);

    for (int iNode = 0; iNode < nReducedNodes; ++iNode) {
      if (thisCPU == 0) fprintf(outFile,"%8.15e\n", sampledVec[iNode]);
    }

    delete [] sampledVec;
  
  }

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::outputReducedToFullNodes() {
/*
  int sp = strlen(ioData->output.rom.prefix);

  const char *fileName;
  const char *fileNameExtension;
  determineFileName(ioData->output.rom.reducedfullnodemap, ".reducedFullNodeMap",fileName,fileNameExtension);
  char *outMeshFile = new char[sp + strlen(fileName)+strlen(fileNameExtension)+1];
  if (thisCPU ==0) sprintf(outMeshFile, "%s%s%s", ioData->output.rom.prefix, fileName, fileNameExtension);
  FILE *outMesh;
  if (thisCPU ==0) outMesh = fopen(outMeshFile, "wt");

  // save the reduced node number for the sample node
  for (int j = 0; j < globalNodes[0].size(); ++j) {
    com->fprintf(outMesh,"%d \n", globalNodes[0][j]+1);
  }
  if (thisCPU == 0) fclose(outMesh);
*/
}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::checkConsistency() {

  // PURPOSE: debugging

  if(com->cpuNum() == 0) std::cerr << " *** WARNING: GappyPreprocessing::checkConsistency is not up-to-date\n";

  int numRhs = nSampleNodes * dim;  // number of RHS treated
  double **consistency = new double * [numRhs];
  for (int i = 0; i < numRhs; ++i)
    consistency[i] = new double [numRhs];

  for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis){
    for (int i = 0; i < numRhs; ++i){
      for (int j = 0; j < numRhs; ++j) {
        consistency[i][j] = 0.0;
      }
    }

    for (int i = 0; i < nPod[iPodBasis]; ++i){
      for (int j = 0; j < nPod[iPodBasis]; ++j) {
        int counterk = 0;
        for (int k = 0; k < numRhs; ++k) {
          int sampleNodeIndex = k/dim;
          int iDim = k%dim;
          int currentSampleNode = globalSampleNodes[sampleNodeIndex];
          int currentCpu = globalNodeToCpuMap.find(currentSampleNode)->second; 
          if (thisCPU == currentCpu){
            int currentSub = globalNodeToLocSubDomainsMap.find(currentSampleNode)->second; 
            int currentLocNode = globalNodeToLocalNodesMap.find(currentSampleNode)->second; 
            assert(domain.getSubDomain()[currentSub]->getNodeMap()[currentLocNode] == currentSampleNode);
            SubDomainData<dim> locValue;
             locValue  = podHat[iPodBasis][j].subData(currentSub);
            consistency[i][j] += podHatPseudoInv[iPodBasis][k][i]*locValue[currentLocNode][iDim];
            ++counterk;
          }

        }
        com->globalSum(1, &counterk);
        assert(counterk==numRhs);
      }
      com->globalSum(nPod[iPodBasis], consistency[i]);
    }
    int asdf = 0;
  }

  for (int i = 0; i < numRhs; ++i) {
    if (consistency[i]) {
      delete [] consistency[i];
      consistency[i] = NULL;
    }
  }
  if (consistency) {
    delete [] consistency;
    consistency = NULL;
  }

}

//----------------------------------------------

template<int dim>
void GappyPreprocessing<dim>::formMaskedNonlinearROBs()
{
  // initialize (from setUpGreedy())
  //domain.makeSampledNodeDistInfo(cpuSample, locSubSample);
  for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis){
    podHat[iPodBasis].resize(nPod[iPodBasis]);
    for (int i = 0; i < nPod[iPodBasis]; ++i) podHat[iPodBasis][i] = 0.0;
  }

  // fill (from findMaxAndFillPodHat())
  std::vector<int> locSampleNodeCount(numLocSub, 0);
  for (int iSampleNode=0; iSampleNode<nSampleNodes; ++iSampleNode) {
    int currentSampleNode = globalSampleNodes[iSampleNode];
    int locSub = globalNodeToLocSubDomainsMap.find(currentSampleNode)->second;
    int locNode = globalNodeToLocalNodesMap.find(currentSampleNode)->second;
    int currentCpu = globalNodeToCpuMap.find(currentSampleNode)->second;
    if (thisCPU == currentCpu) {
      int locSampleNode = locSampleNodeCount[locSub];
      locSampleNodeCount[locSub]++;
      SubDomainData<dim> locPod, locPodHat;
      for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis) {
        for (int iPod = 0 ; iPod < nPod[iPodBasis]; ++iPod) {
          locPod = pod[iPodBasis][iPod].subData(locSub);  // cannot access iDim entry
          locPodHat = podHat[iPodBasis][iPod].subData(locSub);
          for (int iDim = 0; iDim < dim ; ++iDim) {
            locPodHat[locSampleNode][iDim] = locPod[locNode][iDim];
          }
        }
      }
    }
    com->barrier();
  }

  for (int iPodBasis = 0; iPodBasis < nPodBasis; ++iPodBasis){
    pod[iPodBasis].resize(0);  // pod no longer needed
  }

}

//----------------------------------------------


template<int dim>
void GappyPreprocessing<dim>::formReducedSampleNodeMap()
{
  reducedSampleNodeRankMap.clear();

  for (int j = 0; j < nReducedNodes; ++j) {
    
    // compute xyz position of the node
    int globalNodeNum = globalNodes[0][j];  // global node numbers on sample mesh have been sorted in increasing order

    // determine if a sample node
    boost::unordered_map<int, int>::const_iterator sampleNodeMapLoc = globalSampleNodeRankMap.find(globalNodeNum);
    if ( sampleNodeMapLoc != globalSampleNodeRankMap.end()) {
      int globalNodeRank = sampleNodeMapLoc->second;
      reducedSampleNodes[globalNodeRank] = j;  // the globalNodeRank sample node is node j in sample mesh
      reducedSampleNodeRankMap.insert(pair<int,int>( j , globalNodeRank));
    }
  }
}

//----------------------------------------------


template<int dim>
void GappyPreprocessing<dim>::reinitializeMapsForSampleNodes() {

  // Initializes maps when reading sample nodes from file (otherwise they should have been initialized during the greedy algorithm)
  // note: I'm assuming that globalSampleNodes is properly defined before this function is called

  assert( int(globalSampleNodes.size()) > 0);
  nSampleNodes = globalSampleNodes.size();

  cpuSample.clear();
  cpuSample.reserve(nSampleNodes);
  locSubSample.clear();
  locSubSample.reserve(nSampleNodes);
  locNodeSample.clear();
  locNodeSample.reserve(nSampleNodes);
  

  globalNodeToCpuMap.clear();
  globalNodeToLocSubDomainsMap.clear();
  globalNodeToLocalNodesMap.clear();
  nodesXYZmap.clear();

  com->barrier();
  for (int iSampleNode = 0; iSampleNode < nSampleNodes; ++iSampleNode) {
    int cpuTmp = 0;
    int subDTmp = 0;
    int locNodeTmp = 0;
    int globalSampleNode = globalSampleNodes[iSampleNode];
    double xyz [3];
    for (int i=0; i<3; ++i)
      xyz[i]=0.0;

    SubDomainData<dim> locPodHat, locPodHatTmp;
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      bool foundNode = false;
      int nLocNodes = nodeDistInfo.subSize(iSub);       // number of nodes in this subdomain
      int *locToGlobNodeMap = subD[iSub]->getNodeMap();
      bool *locMasterFlag = nodeDistInfo.getMasterFlag(iSub); // master nodes on subdomain
      for (int iLocNode = 0; iLocNode < subD[iSub]->numNodes(); ++iLocNode) {   // all local globalNodes in subdomain
        if (locToGlobNodeMap[iLocNode] == globalSampleNode && locMasterFlag[iLocNode]) {
          cpuTmp = thisCPU;
          subDTmp = iSub;
          locNodeTmp = iLocNode;
          computeXYZ(iSub, iLocNode, xyz);
          foundNode = true;
          break;
        }
      }
      if (foundNode == true) {
        break;
      }
    }

    com->barrier();
    com->globalSum(1,&cpuTmp);
    com->globalSum(1,&subDTmp);
    com->globalSum(1,&locNodeTmp);
    com->globalSum(3, xyz);

    cpuSample.push_back(cpuTmp);
    locSubSample.push_back(subDTmp);
    locNodeSample.push_back(locNodeTmp);

    globalNodeToCpuMap.insert(pair<int, int > (globalSampleNode, cpuTmp));
    globalNodeToLocSubDomainsMap.insert(pair<int, int > (globalSampleNode, subDTmp));
    globalNodeToLocalNodesMap.insert(pair<int, int > (globalSampleNode, locNodeTmp));
    StaticArray<double, 3> XYZ(xyz);
    nodesXYZmap.insert(pair<int, StaticArray <double, 3> > (globalSampleNode, XYZ));
  }

  domain.makeSampledNodeDistInfo(cpuSample, locSubSample);
  domain.makeOldSampledNodeDistInfo(cpuSample, locSubSample);
}


//----------------------------------------------


template<int dim>
void GappyPreprocessing<dim>::setUpSampledNodeTargetRegionMask() {

  // if there's a wall function then no nodes sit exactly on the wall
  //double d2WallMin = (ioData->bc.wall.delta>0.0) ? ioData->bc.wall.delta : 0.0;
  //  set up surface mask
  //surfaceMask = new DistVec<double>( domain.getNodeDistInfo() );
  //*surfaceMask = 0.0;
  //for (int iSub=0; iSub<numLocSub; ++iSub) {
  //  double *d2w = d2wall->subData(iSub);
  //  double *mask = surfaceMask->subData(iSub);
  //  for (int i=0; i<d2wall->subSize(iSub); ++i) {
  //    if (d2w[i]<=(d2WallMin+1e-10))  mask[i] = 1.0;
  //  }
  //}
  //com->fprintf(stdout," ... surfaceMask norm = %e\n", surfaceMask->norm());


  // form target region mask

  double spheres[20][4];
  double boxes[20][2][3];
  double cones[20][2][4];
  int nspheres = 0, nboxes = 0, ncones = 0;

  for (int j=0; j<gappyIO->sampledMeshTargetRegion.num; ++j) {
    if (gappyIO->sampledMeshTargetRegion.spheres[j]->r > 0.0) {
      spheres[nspheres][0] = gappyIO->sampledMeshTargetRegion.spheres[j]->x0;
      spheres[nspheres][1] = gappyIO->sampledMeshTargetRegion.spheres[j]->y0;
      spheres[nspheres][2] = gappyIO->sampledMeshTargetRegion.spheres[j]->z0;
      spheres[nspheres][3] = gappyIO->sampledMeshTargetRegion.spheres[j]->r;
      ++nspheres;
      if (gappyIO->sampledMeshTargetRegion.symmetry == SchemeFixData::X) {
        spheres[nspheres][0] = - gappyIO->sampledMeshTargetRegion.spheres[j]->x0;
        spheres[nspheres][1] = gappyIO->sampledMeshTargetRegion.spheres[j]->y0;
        spheres[nspheres][2] = gappyIO->sampledMeshTargetRegion.spheres[j]->z0;
        spheres[nspheres][3] = gappyIO->sampledMeshTargetRegion.spheres[j]->r;
        ++nspheres;
      }
      else if (gappyIO->sampledMeshTargetRegion.symmetry == SchemeFixData::Y) {
        spheres[nspheres][0] = gappyIO->sampledMeshTargetRegion.spheres[j]->x0;
        spheres[nspheres][1] = - gappyIO->sampledMeshTargetRegion.spheres[j]->y0;
        spheres[nspheres][2] = gappyIO->sampledMeshTargetRegion.spheres[j]->z0;
        spheres[nspheres][3] = gappyIO->sampledMeshTargetRegion.spheres[j]->r;
        ++nspheres;
      }
      else if (gappyIO->sampledMeshTargetRegion.symmetry == SchemeFixData::Z) {
        spheres[nspheres][0] = gappyIO->sampledMeshTargetRegion.spheres[j]->x0;
        spheres[nspheres][1] = gappyIO->sampledMeshTargetRegion.spheres[j]->y0;
        spheres[nspheres][2] = - gappyIO->sampledMeshTargetRegion.spheres[j]->z0;
        spheres[nspheres][3] = gappyIO->sampledMeshTargetRegion.spheres[j]->r;
        ++nspheres;
      }
    }
    if (gappyIO->sampledMeshTargetRegion.boxes[j]->x0 < gappyIO->sampledMeshTargetRegion.boxes[j]->x1) {
      boxes[nboxes][0][0] = gappyIO->sampledMeshTargetRegion.boxes[j]->x0;
      boxes[nboxes][0][1] = gappyIO->sampledMeshTargetRegion.boxes[j]->y0;
      boxes[nboxes][0][2] = gappyIO->sampledMeshTargetRegion.boxes[j]->z0;
      boxes[nboxes][1][0] = gappyIO->sampledMeshTargetRegion.boxes[j]->x1;
      boxes[nboxes][1][1] = gappyIO->sampledMeshTargetRegion.boxes[j]->y1;
      boxes[nboxes][1][2] = gappyIO->sampledMeshTargetRegion.boxes[j]->z1;
      ++nboxes;
      if (gappyIO->sampledMeshTargetRegion.symmetry == SchemeFixData::X) {
        boxes[nboxes][0][0] = -gappyIO->sampledMeshTargetRegion.boxes[j]->x1;
        boxes[nboxes][0][1] = gappyIO->sampledMeshTargetRegion.boxes[j]->y0;
        boxes[nboxes][0][2] = gappyIO->sampledMeshTargetRegion.boxes[j]->z0;
        boxes[nboxes][1][0] = -gappyIO->sampledMeshTargetRegion.boxes[j]->x0;
        boxes[nboxes][1][1] = gappyIO->sampledMeshTargetRegion.boxes[j]->y1;
        boxes[nboxes][1][2] = gappyIO->sampledMeshTargetRegion.boxes[j]->z1;
        ++nboxes;
      }
      if (gappyIO->sampledMeshTargetRegion.symmetry == SchemeFixData::Y) {
        boxes[nboxes][0][0] = gappyIO->sampledMeshTargetRegion.boxes[j]->x0;
        boxes[nboxes][0][1] = -gappyIO->sampledMeshTargetRegion.boxes[j]->y1;
        boxes[nboxes][0][2] = gappyIO->sampledMeshTargetRegion.boxes[j]->z0;
        boxes[nboxes][1][0] = gappyIO->sampledMeshTargetRegion.boxes[j]->x1;
        boxes[nboxes][1][1] = -gappyIO->sampledMeshTargetRegion.boxes[j]->y0;
        boxes[nboxes][1][2] = gappyIO->sampledMeshTargetRegion.boxes[j]->z1;
        ++nboxes;
      }
     if (gappyIO->sampledMeshTargetRegion.symmetry == SchemeFixData::Z) {
        boxes[nboxes][0][0] = gappyIO->sampledMeshTargetRegion.boxes[j]->x0;
        boxes[nboxes][0][1] = gappyIO->sampledMeshTargetRegion.boxes[j]->y0;
        boxes[nboxes][0][2] = -gappyIO->sampledMeshTargetRegion.boxes[j]->z1;
        boxes[nboxes][1][0] = gappyIO->sampledMeshTargetRegion.boxes[j]->x1;
        boxes[nboxes][1][1] = gappyIO->sampledMeshTargetRegion.boxes[j]->y1;
        boxes[nboxes][1][2] = -gappyIO->sampledMeshTargetRegion.boxes[j]->z0;
        ++nboxes;
      }
    }
    if (gappyIO->sampledMeshTargetRegion.cones[j]->r0 >= 0.0 && gappyIO->sampledMeshTargetRegion.cones[j]->r1 >= 0.0) {
      cones[ncones][0][0] = gappyIO->sampledMeshTargetRegion.cones[j]->x0;
      cones[ncones][0][1] = gappyIO->sampledMeshTargetRegion.cones[j]->y0;
      cones[ncones][0][2] = gappyIO->sampledMeshTargetRegion.cones[j]->z0;
      cones[ncones][0][3] = gappyIO->sampledMeshTargetRegion.cones[j]->r0;
      cones[ncones][1][0] = gappyIO->sampledMeshTargetRegion.cones[j]->x1;
      cones[ncones][1][1] = gappyIO->sampledMeshTargetRegion.cones[j]->y1;
      cones[ncones][1][2] = gappyIO->sampledMeshTargetRegion.cones[j]->z1;
      cones[ncones][1][3] = gappyIO->sampledMeshTargetRegion.cones[j]->r1;
      ++ncones;
      if (gappyIO->sampledMeshTargetRegion.symmetry == SchemeFixData::X) {
        cones[ncones][0][0] = -gappyIO->sampledMeshTargetRegion.cones[j]->x0;
        cones[ncones][0][1] = gappyIO->sampledMeshTargetRegion.cones[j]->y0;
        cones[ncones][0][2] = gappyIO->sampledMeshTargetRegion.cones[j]->z0;
        cones[ncones][0][3] = gappyIO->sampledMeshTargetRegion.cones[j]->r0;
        cones[ncones][1][0] = -gappyIO->sampledMeshTargetRegion.cones[j]->x1;
        cones[ncones][1][1] = gappyIO->sampledMeshTargetRegion.cones[j]->y1;
        cones[ncones][1][2] = gappyIO->sampledMeshTargetRegion.cones[j]->z1;
        cones[ncones][1][3] = gappyIO->sampledMeshTargetRegion.cones[j]->r1;
        ++ncones;
      }
      if (gappyIO->sampledMeshTargetRegion.symmetry == SchemeFixData::Y) {
        cones[ncones][0][0] = gappyIO->sampledMeshTargetRegion.cones[j]->x0;
        cones[ncones][0][1] = -gappyIO->sampledMeshTargetRegion.cones[j]->y0;
        cones[ncones][0][2] = gappyIO->sampledMeshTargetRegion.cones[j]->z0;
        cones[ncones][0][3] = gappyIO->sampledMeshTargetRegion.cones[j]->r0;
        cones[ncones][1][0] = gappyIO->sampledMeshTargetRegion.cones[j]->x1;
        cones[ncones][1][1] = -gappyIO->sampledMeshTargetRegion.cones[j]->y1;
        cones[ncones][1][2] = gappyIO->sampledMeshTargetRegion.cones[j]->z1;
        cones[ncones][1][3] = gappyIO->sampledMeshTargetRegion.cones[j]->r1;
        ++ncones;
      }
      if (gappyIO->sampledMeshTargetRegion.symmetry == SchemeFixData::Z) {
        cones[ncones][0][0] = gappyIO->sampledMeshTargetRegion.cones[j]->x0;
        cones[ncones][0][1] = gappyIO->sampledMeshTargetRegion.cones[j]->y0;
        cones[ncones][0][2] = -gappyIO->sampledMeshTargetRegion.cones[j]->z0;
        cones[ncones][0][3] = gappyIO->sampledMeshTargetRegion.cones[j]->r0;
        cones[ncones][1][0] = gappyIO->sampledMeshTargetRegion.cones[j]->x1;
        cones[ncones][1][1] = gappyIO->sampledMeshTargetRegion.cones[j]->y1;
        cones[ncones][1][2] = -gappyIO->sampledMeshTargetRegion.cones[j]->z1;
        cones[ncones][1][3] = gappyIO->sampledMeshTargetRegion.cones[j]->r1;
        ++ncones;
      }
    }
  }

  if (nspheres > 0 || nboxes > 0 || ncones > 0) {
    targetRegionMask = new DistVec<double>( domain.getNodeDistInfo() );
    *targetRegionMask = 1e-16;

    for (int j=0; j<nspheres; ++j)
      com->fprintf(stdout," ... target sampled mesh region includes the sphere [center: (%g, %g, %g), radius: %g]\n",
                  spheres[j][0], spheres[j][1], spheres[j][2], spheres[j][3]);
    for (int j=0; j<nboxes; ++j)
      com->fprintf(stdout," ... target sampled mesh region includes the box [corner1: (%g, %g, %g), corner2: (%g, %g, %g)]\n",
                  boxes[j][0][0], boxes[j][0][1], boxes[j][0][2],
                  boxes[j][1][0], boxes[j][1][1], boxes[j][1][2]);

    for (int j=0; j<ncones; ++j)
      com->fprintf(stdout," ... target sampled mesh region includes the cone [center1: (%g, %g, %g), radius1: %g; center2: (%g, %g, %g), radius2: %g]\n",
                  cones[j][0][0], cones[j][0][1], cones[j][0][2], cones[j][0][3],
                  cones[j][1][0], cones[j][1][1], cones[j][1][2], cones[j][1][3]);


    DistSVec<double,3> X0(domain.getNodeDistInfo());
    domain.getReferenceMeshPosition(X0);

#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      double* mask = targetRegionMask->subData(iSub);
      double (*x0)[3] = X0.subData(iSub);
      for (int i=0; i<targetRegionMask->subSize(iSub); ++i) {
        for (int j=0; j<nspheres; ++j) {
          double r = sqrt( (x0[i][0] - spheres[j][0])*(x0[i][0] - spheres[j][0]) +
                           (x0[i][1] - spheres[j][1])*(x0[i][1] - spheres[j][1]) +
                           (x0[i][2] - spheres[j][2])*(x0[i][2] - spheres[j][2]) );
          if (r <= spheres[j][3])
            mask[i] = 1.0;
        }
        for (int j=0; j<nboxes; ++j) {
          if ((x0[i][0] >= boxes[j][0][0]) && (x0[i][0] <= boxes[j][1][0]) &&
              (x0[i][1] >= boxes[j][0][1]) && (x0[i][1] <= boxes[j][1][1]) &&
              (x0[i][2] >= boxes[j][0][2]) && (x0[i][2] <= boxes[j][1][2]))
            mask[i] = true;
        }
        for (int j=0; j<ncones; ++j)  {
          Vec3D dr(cones[j][1][0]-cones[j][0][0], cones[j][1][1]-cones[j][0][1], cones[j][1][2]-cones[j][0][2]);
          double height = dr.norm();
          dr /= height;
          Vec3D xp;
          Vec3D pr0(x0[i][0]-cones[j][0][0], x0[i][1]-cones[j][0][1], x0[i][2]-cones[j][0][2]);
          double h = pr0*dr;
          if (h >= 0.0 && h <= height)  {
            xp = pr0 - (h*dr);
            double r = cones[j][0][3] + (cones[j][1][3]-cones[j][0][3]) * h / height;
            if (xp.norm() < r)
              mask[i] = true;
          }
        }
      }
    }
    com->fprintf(stdout," ... targetRegionMask = %e\n", targetRegionMask->norm());
  }


  
}

//----------------------------------------------


template<int dim>
void GappyPreprocessing<dim>::maskError() {

  if (onlyInletOutletBC) return; 

  int minLocSampledNodesOnSurfaceInTargetRegion = int(ceil(double(nSampleNodes)*gappyIO->minFractionOfSampledNodesOnSurfaceInTargetRegion));
  int minLocSampledNodesInTargetRegion = int(ceil(double(nSampleNodes)*gappyIO->minFractionOfSampledNodesInTargetRegion));

  if (int(globalSampleNodes.size()) < minLocSampledNodesOnSurfaceInTargetRegion) {
    com->fprintf(stdout," ... only considering surface nodes in target region ...\n");
    DistVec<double> tempWallMask(domain.getNodeDistInfo());
    tempWallMask = 1e-16;
    tempWallMask = tempWallMask + *wallMask;  // findAndFill is dramatically faster if this is nonzero
    if (targetRegionMask) {
      tempWallMask = tempWallMask + *targetRegionMask;
      tempWallMask *= 0.5;
    }
    for (int iPodBasis = 0; iPodBasis < nPodBasis ; ++iPodBasis) {
      for (int iRhs = 0; iRhs < nRhs[iPodBasis] ; ++iRhs){
        error[iPodBasis][iRhs] *= tempWallMask;
        //error[iPodBasis][iRhs] *= *wallMask;
        //if (targetRegionMask) error[iPodBasis][iRhs] *= *targetRegionMask;
      }
    }
  } else if (int(globalSampleNodes.size()) < minLocSampledNodesInTargetRegion) {
    com->fprintf(stdout," ... only considering target region ...\n");
    for (int iPodBasis = 0; iPodBasis < nPodBasis ; ++iPodBasis) {
      for (int iRhs = 0; iRhs < nRhs[iPodBasis] ; ++iRhs){
        if (targetRegionMask) error[iPodBasis][iRhs] *= *targetRegionMask;
      }
    }
  }

}
