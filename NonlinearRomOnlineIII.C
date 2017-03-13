#include <NonlinearRomOnlineIII.h>
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

extern "C"      {
   void F77NAME(dsvdc)(double *, int &, int &, int&, double *,
                        double *, double *, int &, double *, int &,
                        double *, const int &, int &);
}

template<int dim> 
NonlinearRomOnlineIII<dim>::NonlinearRomOnlineIII(Communicator* _com, IoData& _ioData, Domain& _domain, std::vector<double>* _weights)  : 
  NonlinearRom<dim>(_com, _ioData, _domain)
{ 

  if (_weights)
    this->interpWeightsForMultiIC = *_weights;

  // this->ioData->example, this->com->example, this->domain.example

  readClosestCenterInfoModelIII();
  
  if (this->ioData->romOnline.systemApproximation == NonlinearRomOnlineData::GNAT)
    this->determineNumResJacMat();

  if (this->ioData->romOnline.storeAllClusters==NonlinearRomOnlineData::STORE_ALL_CLUSTERS_TRUE)
    this->readAllClusteredOnlineQuantities(); // all ROBs, sVals, and update info

  this->readNonClusteredUpdateInfo("sampled");

}

//----------------------------------------------------------------------------------

template<int dim> 
NonlinearRomOnlineIII<dim>::~NonlinearRomOnlineIII() 
{
  
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomOnlineIII<dim>::readClosestCenterInfoModelIII() {


  if (this->nClusters == 1) return;

  if (this->ioData->romOnline.distanceComparisons) {
    switch (this->ioData->romOnline.basisUpdates) {
      case (NonlinearRomOnlineData::UPDATES_OFF):
        if (this->ioData->romOnline.projectSwitchStateOntoAffineSubspace!=NonlinearRomOnlineData::PROJECT_OFF) {
          this->readDistanceComparisonInfo("project");
        } else {
          this->readDistanceComparisonInfo("noUpdates");
        }
        break;
      case (NonlinearRomOnlineData::UPDATES_SIMPLE):
        this->com->fprintf(stderr, "*** Error: fast distance caluclations are incompatible with simple ROB updates (use Exact)\n");
        exit(-1);
        break;
      case (NonlinearRomOnlineData::UPDATES_FAST_EXACT):
        if (this->ioData->romOnline.sensitivity.include==NonlinearRomOnlineNonStateData::INCLUDE_ON ||
            this->ioData->romOnline.krylov.include==NonlinearRomOnlineNonStateData::INCLUDE_ON) {
          this->com->fprintf(stderr, "*** Error: exact updates not yet supported for heterogeneous ROBs\n");
          exit(-1);
        }
        this->readDistanceComparisonInfo("exactUpdates");
        break;
      case (NonlinearRomOnlineData::UPDATES_FAST_APPROX):
        this->readDistanceComparisonInfo("approxUpdates");
        break;
      default:
        this->com->fprintf(stderr, "*** Error: Unexpected ROB updates method\n");
        exit(-1);
    }
  } else {
    this->com->fprintf(stderr, "*** Error: fast cluster selection must be used for local hyper-reduction\n");
    exit(-1);
  }

}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomOnlineIII<dim>::readClusteredOnlineQuantities(int iCluster) {

    // read in sampled state ROB
    this->readClusteredBasis(iCluster, "sampledState");

    // read in sampled krylov ROB

    // read in sampled sensitivity ROB

    if ((this->ioData->problem.alltype == ProblemData::_STEADY_NONLINEAR_ROM_POST_
           || this->ioData->problem.alltype == ProblemData::_UNSTEADY_NONLINEAR_ROM_POST_)
           && strcmp(this->surfaceRefStateName,"")!=0) {
      return;
    }

    // read in sample nodes
    this->readClusteredSampleNodes(iCluster, "sampled");

    if (this->ioData->romOnline.systemApproximation == NonlinearRomOnlineData::GNAT) {
      this->readClusteredGappyMatrix(iCluster, "resMatrix");  // read in gappy POD matrix for residual
      if (this->numResJacMat==2) this->readClusteredGappyMatrix(iCluster, "jacMatrix"); // read in gappy POD matrix for Jacobian
    }

    if (this->ioData->romOnline.systemApproximation == NonlinearRomOnlineData::APPROX_METRIC_NL) {
      this->readClusteredNonlinearMetric(iCluster);
    }
 
}


//----------------------------------------------------------------------------------

template<int dim>
bool NonlinearRomOnlineIII<dim>::updateBasis(int iCluster, DistSVec<double, dim> &U, Vec<double>* coords) {

bool updatePerformed;

 switch (this->ioData->romOnline.basisUpdates) {
      case (NonlinearRomOnlineData::UPDATES_OFF):
        break;
      case (NonlinearRomOnlineData::UPDATES_SIMPLE):
        this->com->fprintf(stderr,"*** Error: simple basis updates are incompatible with hyper-reduction\n");
        exit(-1); 
        break;
      case (NonlinearRomOnlineData::UPDATES_FAST_EXACT):
        this->com->fprintf(stdout, " ... Applying exact rank one basis update using precomputed quantities\n");
        updatePerformed = updateBasisFastExact(iCluster, U, coords);
        break;
      case (NonlinearRomOnlineData::UPDATES_FAST_APPROX):
        this->com->fprintf(stdout, " ... Applying rank one basis update using approximate metric\n");
        updatePerformed = updateBasisFastApprox(iCluster, U);
        break;
      default:
        this->com->fprintf(stderr, "*** Error: Unexpected ROB updates method\n");
        exit(-1);
   }

   return updatePerformed;

}

//----------------------------------------------------------------------------------
    
template<int dim>
bool NonlinearRomOnlineIII<dim>::updateBasisFastExact(int currentCluster, DistSVec<double, dim> &U, Vec<double>* coords) {
    
/* This function will use notation from Amsallem et al 2013 

  Input:
    - Original bases, singlar values, and columnSumsV (stored if storeAllClusteredOnlineQuantities is on, otherwise read on demand)
    - Precomputed Quantities 
      > a = this->uicNorm              // double
      > b = this->urefNorms            // [iCluster]
      > c = this->urefUicProducts;     // [iCluster]
      > d = this->basisUicProducts;    // [iCluster][1:nPod]
      > e = this->basisUrefProducts;   // [Cluster_Basis][Cluster_Uref][:]
      > F = this->basisBasisProducts;  // [iCluster][pCluster][:][:] symmetric wrt clusters (strict lower triangular with identity assumed on the diagonal)
      > g = this->urefUrefProducts;    // [iCluster][jCluster] symmetric (lower triangular)
    - Masked initial condition of simulation (this->Uic)
    - Masked reference states (stored if storeAllClusteredOnlineQuantities is on, otherwise read on demand)
    - Current online state reduced coordinates (coords)
    - Current exact update quantities 
      > alpha = this->exactUpdatesAlpha; // [jVec]
      > beta  = this->exactUpdatesBeta;  // [iCluster][jVec]
      > N     = this->exactUpdatesN;     // [iCluster][iVec][jVec]
      > alpha_switch = this->exactUpdatesAlphaSwitch // (double) 
      > beta_switch  = this->exactUpdatesBetaSwitch  // [iCluster]
      > n_switch     = this->exactUpdatesNSwitch     // [iCluster][iVec]

  Output:
    - Updated basis
    - Updated exact update quantities

*/

  bool updatePerformed = false;

  this->readClusteredUpdateInfo(currentCluster, "sampledState");
  this->readClusteredReferenceState(currentCluster, "sampledState"); // reads Uref

  int robSize = this->basis->numVectors();
  int kSize = robSize+1;

  // use Kahan summation to minimize roundoff error.  Volatile should prevent intel compiler from using unsafe (associative) optimizations here.
  volatile double compensation;
  volatile double temp1;
  volatile double temp2;
  double totalCompensation = 0.0;

  // update alpha_switch, beta_switch, and n_switch
  // kVec is used for the dimension of the previous rob, jVec for current, iVec for cluster i
  assert((this->exactUpdatesAlpha.size() == 0) || (this->exactUpdatesAlpha.size() == coords->size()) );
  for (int kVec=0; kVec<this->exactUpdatesAlpha.size(); ++kVec) {// note that exactUpdatesAlpha might be empty, indicating that it is all zeros
    this->exactUpdatesAlphaSwitch += this->exactUpdatesAlpha[kVec]*(*coords)[kVec];
  }
  for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
    assert((this->exactUpdatesBeta[iCluster].size() == 0) || (this->exactUpdatesBeta[iCluster].size() == coords->size()));
    compensation = 0.0;
    for (int kVec=0; kVec<this->exactUpdatesBeta[iCluster].size(); ++kVec) {
        temp1 = this->exactUpdatesBeta[iCluster][kVec]*(*coords)[kVec] - compensation;
        temp2 = this->exactUpdatesBetaSwitch[iCluster] + temp1;
        compensation = (temp2 - this->exactUpdatesBetaSwitch[iCluster]) - temp1;
        totalCompensation += compensation;
        this->exactUpdatesBetaSwitch[iCluster] = temp2;
        //this->exactUpdatesBetaSwitch[iCluster] += this->exactUpdatesBeta[iCluster][kVec]*(*coords)[kVec];
    }

    assert((this->exactUpdatesN[iCluster].size() == 0) || (this->exactUpdatesN[iCluster][0].size() == coords->size()));
    compensation = 0.0;
    for (int iVec=0; iVec<this->exactUpdatesN[iCluster].size(); ++iVec) {
      if (this->exactUpdatesNSwitch[iCluster].size()==0) {
        this->exactUpdatesNSwitch[iCluster].resize(this->exactUpdatesN[iCluster].size(),0.0);
      }
      for (int kVec=0; kVec<this->exactUpdatesN[iCluster][iVec].size(); ++kVec) {
          temp1 = this->exactUpdatesN[iCluster][iVec][kVec]*(*coords)[kVec] - compensation;
          temp2 = this->exactUpdatesNSwitch[iCluster][iVec] + temp1;
          compensation = (temp2 - this->exactUpdatesNSwitch[iCluster][iVec]) - temp1;
          totalCompensation += compensation;
          this->exactUpdatesNSwitch[iCluster][iVec] = temp2;
      //    this->exactUpdatesNSwitch[iCluster][iVec] += this->exactUpdatesN[iCluster][iVec][kVec]*(*coords)[kVec];     
      }
    }
  }

  double alphaA = -1.0*this->exactUpdatesAlphaSwitch;
  std::vector<double> betaA = this->exactUpdatesBetaSwitch;
  std::vector<std::vector<double> > nA = this->exactUpdatesNSwitch;

  for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
    betaA[iCluster] = -1.0*this->exactUpdatesBetaSwitch[iCluster];
    for (int iVec=0; iVec<this->exactUpdatesNSwitch[iCluster].size(); ++iVec) {
      nA[iCluster][iVec] = -1.0*this->exactUpdatesNSwitch[iCluster][iVec];
    }
  }
  betaA[currentCluster] += 1.0;

  std::vector<double> sumVec;

  std::vector<double> m;
  m.resize(robSize,0.0);
  assert(this->basisUicProducts[currentCluster].size() >= robSize); // preprocessing might have used fewer vectors

  for (int jVec=0; jVec<robSize; ++jVec) {
    sumVec.push_back( alphaA*(this->basisUicProducts[currentCluster][jVec]) );
    for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
      assert(this->basisUrefProducts[currentCluster][iCluster].size() >= robSize);  // again, checking preprocessing
      sumVec.push_back( betaA[iCluster]*(this->basisUrefProducts[currentCluster][iCluster][jVec]) );
      if (iCluster < currentCluster) {
        assert(this->basisBasisProducts[currentCluster][iCluster].size() >= robSize);
        assert(this->basisBasisProducts[currentCluster][iCluster][0].size() >= nA[iCluster].size());
      } else if (iCluster > currentCluster) {
        assert(this->basisBasisProducts[iCluster][currentCluster].size() >= nA[iCluster].size());
        assert(this->basisBasisProducts[iCluster][currentCluster][0].size() >= robSize);
      }
      for (int iVec=0; iVec<nA[iCluster].size(); ++iVec) {
        if (iCluster < currentCluster) {
          sumVec.push_back( (this->basisBasisProducts[currentCluster][iCluster][jVec][iVec])*nA[iCluster][iVec] );
        } else if (iCluster > currentCluster) { // basisBasis products is lower triangular, so swap iCluster and currentCluster and take transpose
          sumVec.push_back( (this->basisBasisProducts[iCluster][currentCluster][iVec][jVec])*nA[iCluster][iVec] );
        } else {// iCluster == currentCluster, basisBasisProducts are assumed to be identity
          if (jVec == iVec) sumVec.push_back( nA[currentCluster][jVec] );
        }
      }
    }
    
    compensation = 0.0;
    for (int iEntry=0; iEntry<sumVec.size(); ++iEntry) {
       temp1 = sumVec[iEntry] - compensation;
       temp2 = m[jVec] + temp1;
       compensation = (temp2 - m[jVec]) - temp1;
       totalCompensation += compensation;
       m[jVec] = temp2;
    }
    sumVec.clear();
  }

  std::vector<std::vector<double> > nATmp = nA; //nATmp = nA[iCluser] - delta_(iCluster,currentCluster)*m
  if (nATmp[currentCluster].size()>0) {
    for (int iVec=0; iVec<nA[currentCluster].size(); ++iVec) {
      nATmp[currentCluster][iVec] -= m[iVec];
    }
  } else {
    nATmp[currentCluster].resize(m.size());
    for (int iVec=0; iVec<m.size(); ++iVec) {
      nATmp[currentCluster][iVec] = -1.0*m[iVec];
    }
  }

  // construct r
  sumVec.clear();
  sumVec.push_back( pow(alphaA * this->uicNorm, 2) );

  for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
    for (int pCluster=0; pCluster<this->nClusters; ++pCluster) {
      if (iCluster >= pCluster) {
        sumVec.push_back( betaA[iCluster]*betaA[pCluster]*this->urefUrefProducts[iCluster][pCluster] );
      } else {
        sumVec.push_back( betaA[iCluster]*betaA[pCluster]*this->urefUrefProducts[pCluster][iCluster] );
      }
    }
  }

  for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
    for (int pCluster=0; pCluster<this->nClusters; ++pCluster) {
      for (int iVec=0; iVec<nATmp[iCluster].size(); ++iVec) {
        sumVec.push_back( 2*betaA[pCluster]*this->basisUrefProducts[iCluster][pCluster][iVec]*nATmp[iCluster][iVec] );
      }
    }
  }

  for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
    for (int pCluster=0; pCluster<this->nClusters; ++pCluster) {
      if (iCluster > pCluster) {
        for (int pVec=0; pVec<nATmp[pCluster].size(); ++pVec) {
          for (int iVec=0; iVec<nATmp[iCluster].size(); ++iVec) {
            sumVec.push_back( nATmp[iCluster][iVec] * nATmp[pCluster][pVec] * this->basisBasisProducts[iCluster][pCluster][iVec][pVec] );
          }
        }
      } else if (iCluster < pCluster) {
        for (int pVec=0; pVec<nATmp[pCluster].size(); ++pVec) {
          for (int iVec=0; iVec<nATmp[iCluster].size(); ++iVec) {
            sumVec.push_back( nATmp[iCluster][iVec] * nATmp[pCluster][pVec] * this->basisBasisProducts[pCluster][iCluster][pVec][iVec] );
          }
        }
      } else { // iCluster = pCluster
        for (int iVec=0; iVec<nATmp[iCluster].size(); ++iVec) {
          sumVec.push_back( pow(nATmp[iCluster][iVec],2) );
        }
      }
    }
  }

  for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
    sumVec.push_back( 2.0*alphaA*betaA[iCluster]*this->urefUicProducts[iCluster] );
  }

  for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
    for (int iVec=0; iVec<nATmp[iCluster].size(); ++iVec) {
      sumVec.push_back( 2.0*alphaA*nATmp[iCluster][iVec]*this->basisUicProducts[iCluster][iVec] );
    }
  }

  double r = 0.0;
  compensation = 0.0;
  for (int iEntry=0; iEntry<sumVec.size(); ++iEntry) {
     temp1 = sumVec[iEntry] - compensation;
     temp2 = r + temp1;
     compensation = (temp2 - r) - temp1;
     totalCompensation += compensation;
     r = temp2;
  }
  sumVec.clear();
  //this->com->fprintf(stdout,"totalCompensation = %1.12e\n", totalCompensation);

  if ( r < pow(this->rTol,2.0) ) { // no need to update basis -- just need to update alpha, beta, and N
    this->com->fprintf(stdout, " ... r is less than the specified tolerance of %e -- skipping the update\n", this->rTol);
    this->exactUpdatesAlpha.clear();
    for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
      this->exactUpdatesBeta[iCluster].clear();
      if (iCluster == currentCluster) { // identity
         this->exactUpdatesN[iCluster].resize(robSize);
         for (int iVec=0; iVec<robSize; ++iVec) {
           this->exactUpdatesN[iCluster][iVec].clear();
           this->exactUpdatesN[iCluster][iVec].resize(robSize, 0.0);
           this->exactUpdatesN[iCluster][iVec][iVec] = 1.0;
         }
      } else { // empty
        this->exactUpdatesN[iCluster].clear();
      }
    }
  } else { // update basis 

    updatePerformed = true;

    r = pow(r, 0.5);

    // form alphaP, betaP, and nP
    double alphaP = alphaA/r;
    std::vector<double> betaP = betaA;
    std::vector<std::vector<double> > nP = nATmp;
    for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
      betaP[iCluster] = betaA[iCluster]/r;
      for (int jVec=0; jVec<nATmp[iCluster].size(); ++jVec) {
        nP[iCluster][jVec] = nATmp[iCluster][jVec]/r;
      }
    }

    double *K = new double[(kSize)*(kSize)];
    for (int iCol = 0; iCol < (kSize); ++iCol){
      for (int iRow = 0; iRow < (kSize); ++iRow) {
        if ((iCol == iRow) && (iCol < kSize-1)) {
          K[iCol*(kSize) + iRow] = (*this->sVals)[iCol];
        } else {
          K[iCol*(kSize) + iRow] = 0.0;
        }
      }
    }

    double q = 0;
    for (int iVec=robSize; iVec<(this->columnSumsV->size()); ++iVec) {
      q += pow((*(this->columnSumsV))[iVec], 2);
    }
    q = pow(q, 0.5);

    for (int iRow = 0; iRow < (kSize-1); ++iRow) {
      for (int iCol = 0; iCol < (kSize-1); ++iCol){
        K[iCol*(kSize) + iRow] += m[iRow] * ((*(this->columnSumsV))[iCol]);
      }
      K[(kSize-1)*(kSize) + iRow] = m[iRow] * q;
    }

    for (int iCol = 0; iCol < kSize-1; ++iCol){
      K[iCol*(kSize) + (kSize-1)] += r * ((*(this->columnSumsV))[iCol]);
    }
    K[(kSize-1)*(kSize) + kSize-1] += r * q;

    double *sigma = new double[kSize];
    double *error = new double[kSize];
    double *work = new double[kSize];
    int info;
    double *zVec = new double[kSize*kSize]; // right singular vectors
    double *yVec = new double[kSize*kSize]; // left singular vectors

    this->com->fprintf(stdout, " ... computing rank one update to basis using current state\n");
    F77NAME(dsvdc)(K, kSize, kSize, kSize, sigma, error, yVec, kSize, zVec, kSize, work, 11, info);

    // update alpha, beta, N
    std::vector<double> lastRowOfC; // last row of left singular vectors (excluding the bottom right element)
    lastRowOfC.resize(kSize-1, 0.0);
    for (int iCol = 0; iCol < (kSize-1); ++iCol){
      lastRowOfC[iCol] = yVec[iCol*(kSize) + (kSize-1)];
    }

    this->exactUpdatesAlpha.resize(robSize);
    for (int jVec=0; jVec<robSize; ++jVec) {  // note that robSize = kSize-1
      this->exactUpdatesAlpha[jVec] = alphaP*lastRowOfC[jVec];
    }

    for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
      this->exactUpdatesBeta[iCluster].resize(robSize);
      for (int jVec=0; jVec<robSize; ++jVec) {
        this->exactUpdatesBeta[iCluster][jVec] = betaP[iCluster]*lastRowOfC[jVec];
      }
      if (iCluster == currentCluster) {
        this->exactUpdatesN[iCluster].resize(robSize);
      } else {
        this->exactUpdatesN[iCluster].resize(nP[iCluster].size());
      }
      for (int iVec=0; iVec<this->exactUpdatesN[iCluster].size(); ++iVec) {
        this->exactUpdatesN[iCluster][iVec].clear();
        this->exactUpdatesN[iCluster][iVec].resize(robSize,0.0);
        for (int jVec=0; jVec<robSize; ++jVec) {
          if (iCluster == currentCluster) {
            this->exactUpdatesN[iCluster][iVec][jVec] = yVec[jVec*(kSize) + iVec];
          }
          this->exactUpdatesN[iCluster][iVec][jVec] += nP[iCluster][iVec] * lastRowOfC[jVec];
        }
      }
    }

    // update basis
    VecSet< DistSVec<double, dim> >* updatedBasis = new VecSet< DistSVec<double, dim> >(robSize, this->domain.getNodeDistInfo());
    for (int jVec=0; jVec<robSize; ++jVec) {
      (*updatedBasis)[jVec] =  this->exactUpdatesAlpha[jVec] * (*this->Uic);
    }

    for (int iCluster=0; iCluster<this->nClusters; ++iCluster) {
      this->readClusteredBasis(iCluster, "sampledState");
      this->readClusteredReferenceState(iCluster, "sampledState"); // reads Uref
      for (int jVec=0; jVec<robSize; ++jVec) {
        (*updatedBasis)[jVec] = (*updatedBasis)[jVec] + this->exactUpdatesBeta[iCluster][jVec] * (*this->Uref);
        for (int iVec=0; iVec<this->exactUpdatesN[iCluster].size(); ++iVec) {
          (*updatedBasis)[jVec] = (*updatedBasis)[jVec] + this->exactUpdatesN[iCluster][iVec][jVec] * (*(this->basis))[iVec];
        }
      }
    }

    delete (this->basis);
    delete (this->sVals);
    this->sVals = NULL;

    this->basis = updatedBasis;
    this->nState = this->basis->numVectors();

    delete [] K;
    delete [] sigma;
    delete [] error;
    delete [] work;
    delete [] zVec;
    delete [] yVec;
  }

  return updatePerformed;

}

//----------------------------------------------------------------------------------

template<int dim>
bool NonlinearRomOnlineIII<dim>::updateBasisFastApprox(int iCluster, DistSVec<double, dim> &U) {

/* 
  When updateBasis is called the following quantities are available:
    - nClusters: number of clusters
    - snapsInCluster:  a vector containing the number of snapshots in each cluster, size nClusters
    - columnSumsV: a vector containing the sums of the columns of V, size ny + buffer
    - basis: the local ROB, ny + buffer
    - sVals: the singular values, size ny + buffer
  (As well as all of the ioData values)

*/

  // Approximate updates case

  bool updatePerformed = false;

  this->readClusteredUpdateInfo(iCluster, "sampledState"); 
  this->readClusteredReferenceState(iCluster,"sampledState"); // reads Uref

  int robSize = this->basis->numVectors();
  int kSize = robSize+1;

  DistSVec<double, dim> a(this->domain.getNodeDistInfo());
  a = *(this->Uref) - U;

  double *m = new double[robSize];
  double temp1[this->nLowRankFactors];
  double temp2[this->nLowRankFactors];
  for (int iRank = 0; iRank<this->nLowRankFactors; ++iRank) 
    temp1[iRank] = (*this->lowRankFactor)[iRank] * a;

  for (int iVec=0; iVec<robSize; ++iVec) {
    m[iVec]  = 0.0;
    for (int iRank = 0; iRank<this->nLowRankFactors; ++iRank) {
      temp2[iRank] = (*this->lowRankFactor)[iRank] * (*(this->basis))[iVec];
      m[iVec] += (temp2[iRank]*temp1[iRank]);
    }
  }

  DistSVec<double, dim> p(this->domain.getNodeDistInfo());
  p = a;

  for (int iVec=0; iVec<robSize; ++iVec) 
    p -= (*(this->basis))[iVec] * m[iVec];
    

  double Ra = 0.0;
  for (int iRank = 0; iRank<this->nLowRankFactors; ++iRank) {
    temp1[iRank] = (*this->lowRankFactor)[iRank] * p;
    Ra += pow(temp1[iRank],2);
  }
 
  if (Ra >= pow(this->rTol,2.0)) {  // only update if Uref is different than U (this handles the case of time=0) 
   
    updatePerformed = true;
   
    Ra = pow(Ra,0.5);
    double RaInv = 1/Ra; 
    p *= RaInv;

    double *K = new double[(kSize)*(kSize)];

    for (int iCol = 0; iCol < (kSize); ++iCol){
      for (int iRow = 0; iRow < (kSize); ++iRow) {
        if ((iCol == iRow) && (iCol < kSize-1)) {
          K[iCol*(kSize) + iRow] = (*this->sVals)[iCol];
        } else {
          K[iCol*(kSize) + iRow] = 0.0;
        }
      }
    }

    double q = 0;
    for (int iVec=robSize; iVec<(this->columnSumsV->size()); ++iVec) {
      q += pow((*(this->columnSumsV))[iVec], 2);
    }
    q = pow(q, 0.5);

    for (int iRow = 0; iRow < (kSize-1); ++iRow) {
      for (int iCol = 0; iCol < (kSize-1); ++iCol){
        K[iCol*(kSize) + iRow] += m[iRow] * ((*(this->columnSumsV))[iCol]);
      }
      K[(kSize-1)*(kSize) + iRow] = m[iRow] * q;
    }

    for (int iCol = 0; iCol < kSize-1; ++iCol){
      K[iCol*(kSize) + (kSize-1)] += Ra * ((*(this->columnSumsV))[iCol]);
    }
    K[(kSize-1)*(kSize) + kSize-1] += Ra * q;

    double *sigma = new double[kSize];
    double *error = new double[kSize];
    double *work = new double[kSize];
    int info;
    double *zVec = new double[kSize*kSize]; // right singular vectors
    double *yVec = new double[kSize*kSize]; // left singular vectors

    this->com->fprintf(stdout, " ... computing rank one update to basis using current state\n");
    F77NAME(dsvdc)(K, kSize, kSize, kSize, sigma, error, yVec, kSize, zVec, kSize, work, 11, info);

    VecSet< DistSVec<double, dim> > basisOld(robSize, this->domain.getNodeDistInfo());

    for (int iVec=0; iVec<robSize; ++iVec)
        basisOld[iVec] = (*(this->basis))[iVec];

    for (int iVec=0; iVec<robSize; ++iVec) {
      (*(this->basis))[iVec] = p * yVec[(iVec*kSize) + kSize-1];
      for (int jVec=0; jVec<robSize; ++jVec) {
        (*(this->basis))[iVec] += yVec[(iVec*kSize) + jVec] * basisOld[jVec];
      }
    }

    delete [] K;
    delete [] m;
    delete [] sigma;
    delete [] error;
    delete [] work;
    delete [] zVec;
    delete [] yVec;

  } else {
    this->com->fprintf(stdout, " ... r is less than the specified tolerance of %e -- skipping the update\n", this->rTol);
  }

  delete (this->Uref);
  this->Uref = NULL;

  if (this->ioData->romOnline.distanceComparisons) 
    this->resetDistanceComparisonQuantitiesApproxUpdates();

  return updatePerformed;
}

//----------------------------------------------------------------------------------

template<int dim>
void NonlinearRomOnlineIII<dim>::appendNonStateDataToBasis(int cluster, const char* basisType, bool relProjError) {
/*
  int robSize = this->basis->numVectors();
  VecSet< DistSVec<double, dim> > basisOld(robSize, this->domain.getNodeDistInfo());

  for (int iVec=0; iVec<robSize; ++iVec)
    basisOld[iVec] = (*(this->basis))[iVec];
 
  this->readClusteredBasis(cluster, basisType);
  int nonStateSize = this->basis->numVectors();
  VecSet< DistSVec<double, dim> > nonStateBasis(nonStateSize, this->domain.getNodeDistInfo());

  for (int iVec=0; iVec<nonStateSize; ++iVec)
    nonStateBasis[iVec] = (*(this->basis))[iVec];

  this->basis->resize(robSize + nonStateSize);

  for (int iVec=0; iVec<robSize; ++iVec)
    (*(this->basis))[iVec] = basisOld[iVec];

  for (int iVec=robSize; iVec<(robSize+nonStateSize); iVec++)
    (*(this->basis))[iVec] = nonStateBasis[iVec-robSize];

  bool gramSchmidt;
  if (strcmp(basisType, "krylov")==0) {
    gramSchmidt = (this->ioData->romOnline.krylov.gramSchmidt==NonlinearRomOnlineNonStateData::GRAMSCHMIDT_ON) ? true : false;
  } else if (strcmp(basisType, "sensitivity")==0) {
    gramSchmidt = (this->ioData->romOnline.sensitivity.gramSchmidt==NonlinearRomOnlineNonStateData::GRAMSCHMIDT_ON) ? true : false;
  } else {
    this->com->fprintf(stderr, "*** Error: unexpected basis type passed to appendNonStateDataToBasis (%s)\n",basisType);
    exit(-1);
  }

  if (gramSchmidt) { 
    int uniqueVecs = robSize;
    for (int iVec = robSize; iVec<(robSize+nonStateSize); iVec++) {
      for (int jVec = 0; jVec<iVec; jVec++) {
        (*(this->basis))[iVec] -= (*(this->basis))[jVec] * ((*(this->basis))[iVec] * (*(this->basis))[jVec]);
        double norm = (*(this->basis))[iVec].norm();
        if (norm>=1e-14) {
          (*(this->basis))[iVec] *= 1/norm;
          ++uniqueVecs;
        } else {
          this->com->fprintf(stderr, "*** Warning: removing linearly dependent vector (norm=%e)\n",norm);
          (*(this->basis))[iVec] *= 0;    
        }
      }
    }
    if (uniqueVecs < (robSize+nonStateSize)) {
      basisOld.resize(robSize+nonStateSize);
      for (int iVec=0; iVec<(robSize+nonStateSize); ++iVec)
        basisOld[iVec] = (*(this->basis))[iVec];
      this->basis->resize(uniqueVecs);
      int added = 0;
      for (int iVec=0; iVec<(robSize+nonStateSize); ++iVec) {
        if (basisOld[iVec].norm()>=1e-14) {
          (*(this->basis))[added] = basisOld[iVec];
          ++added;
        }
      }
    }
  }
  */
}


//----------------------------------------------------------------------------------


//-------------------------------------------------------------------

template<int dim>
void NonlinearRomOnlineIII<dim>::projectSwitchStateOntoAffineSubspace(int currentCluster, int prevCluster, 
                                                                      DistSVec<double, dim> &U, Vec<double> &UromCurrentROB) {

  this->com->fprintf(stdout, " ... projecting switch state onto affine subspace\n");

  // if first iteration
  // currentReducedCoords = basisUicProducts[currentCluster] - basisUrefProducts[currentCluster,currentCluster];
  //
  // (if starting with an interpolated Uic, basisUicProducts[currentCluster] is formed in TsDesc.C)
  //
  // otherwise
  // currentReducedCoords = basisUrefProducts[currentCluster,prevCluster] - basisUrefProducts[currentCluster,currentCluster]
  //                      + basisBasisProducts[currentCluster,prevCluster,:,:] * prevReducedCoords;
  //                      
  // newState = currentReference + currentROB * currentReducedCoords;

  int nPodVecs = this->basis->numVectors();
  int nPodVecsPrev = UromCurrentROB.size();

  if (prevCluster==-1) {
    UromCurrentROB.resize(nPodVecs);
    for (int iVec = 0; iVec < nPodVecs; iVec++)
      UromCurrentROB[iVec] = this->basisUicProducts[currentCluster][iVec] - this->basisUrefProducts[currentCluster][currentCluster][iVec];
  } else if (currentCluster==prevCluster) {
    // do nothing -- (should never enter here, just being overly cautious because basisBasisProducts isn't defined for this case...)
  } else {
    Vec<double> UromPrevROB(UromCurrentROB);
    UromCurrentROB.resize(nPodVecs);
    for (int iVec = 0; iVec < nPodVecs; iVec++) {
      UromCurrentROB[iVec] = -1.0* this->basisUrefProducts[currentCluster][currentCluster][iVec];
      //if (prevCluster<currentCluster) {
        UromCurrentROB[iVec] += this->basisUrefProducts[currentCluster][prevCluster][iVec];
      //} else {
      //  UromCurrentROB[iVec] -= this->basisUrefProducts[prevCluster][currentCluster][iVec];
      //}
      if (prevCluster<currentCluster) {
        for (int jVec = 0; jVec < nPodVecsPrev; jVec++) {
          UromCurrentROB[iVec] += this->basisBasisProducts[currentCluster][prevCluster][iVec][jVec] * UromPrevROB[jVec];
        }                  
      } else {
        for (int jVec = 0; jVec < nPodVecsPrev; jVec++) {
          UromCurrentROB[iVec] += this->basisBasisProducts[prevCluster][currentCluster][jVec][iVec] * UromPrevROB[jVec];
        }
      }
    }
  }

  this->readClusteredUpdateInfo(currentCluster, "sampledState");
  U = *(this->Uref);
  for (int iVec = 0; iVec < nPodVecs; iVec++)
    U += UromCurrentROB[iVec] * (*(this->basis))[iVec];

  return;
}

