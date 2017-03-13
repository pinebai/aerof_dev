#include <TsInput.h>
#include <sys/time.h>
#include <cstdlib>
#include <cmath>

using std::stable_sort;

template<class VecType>  
KspBinaryOutput<VecType>::KspBinaryOutput(Communicator *_com, IoData *_ioData, Domain *_domain)  : 
com(_com), ioData(_ioData), domain(_domain)
{ 
  krylovFreqTime = ioData->output.rom.krylovOutputFreqTime;
  krylovFreqNewton = ioData->output.rom.krylovOutputFreqNewton;
  krylovEnergy = ioData->output.rom.krylovVectorEnergy;
  timeIt = domain->getTimeIt();
  newtonIt = domain->getNewtonIt();

  fileName = new char[strlen(ioData->output.rom.prefix) + strlen(ioData->output.rom.krylovVector) + 1];
  sprintf(fileName, "%s%s", ioData->output.rom.prefix, ioData->output.rom.krylovVector);

  U = new VecType(domain->getNodeDistInfo());

}

//----------------------------------------------------------------------------------

template<class VecType> 
KspBinaryOutput<VecType>::~KspBinaryOutput() 
{

  delete [] fileName;
  delete U;
}

//----------------------------------------------------------------------------------

// this struct is used in writeKrylovVectors
struct kspSortStruct {
  int kspIndex; // index of corresponding krylov Vector
  double energy; // distance to second closest cluster

  bool operator<(const kspSortStruct& a) const {
    return energy < a.energy;
  }
};

//----------------------------------------------------------------------------------

template<class VecType>
void KspBinaryOutput<VecType>::setCurrentState(VecType &state) {

  *U = state;

}

//----------------------------------------------------------------------------------

template<class VecType>
template<int dim>
void KspBinaryOutput<VecType>::writeKrylovVectors(VecSet<DistSVec<double, dim> >& kspVecs, Vec<double> kspCoords, int numVecs) {

  if ((*(timeIt)%krylovFreqTime==0) && (*(newtonIt)%krylovFreqNewton==0)) {

    // note: numVecs is the total number of krylov vectors for this newton step (usually less than the max)
    Vec<double> kspCoordsTrunc(numVecs);
    for (int iVec=0; iVec<numVecs; ++iVec) kspCoordsTrunc[iVec] = kspCoords[iVec];
    //kspCoordsTrunc *= (1/kspCoordsTrunc.norm());
    kspSortStruct* kspIndexedCoords = new kspSortStruct[numVecs]; 

    double totalEnergy = 0;

    for (int iVec=0; iVec<numVecs; ++iVec) {
      kspIndexedCoords[iVec].kspIndex = iVec;
      kspIndexedCoords[iVec].energy = pow(kspCoordsTrunc[iVec],2);
      totalEnergy += kspIndexedCoords[iVec].energy;
    }

    sort(kspIndexedCoords, kspIndexedCoords+numVecs); //ascending order
    double cumEnergy = 0;
    int vecsOutput = 0;

    DistSVec<double, dim> outVec(domain->getNodeDistInfo());
    if (ioData->output.rom.addStateToKrylov) outVec = *U;
    
    DistSVec<double, dim> currentKspVec(domain->getNodeDistInfo());
    double currentKspCoord;

    while (((cumEnergy/totalEnergy)<krylovEnergy)&&(vecsOutput<numVecs)) {
      currentKspVec = kspVecs[kspIndexedCoords[numVecs-vecsOutput-1].kspIndex];
      currentKspCoord = kspCoordsTrunc[kspIndexedCoords[numVecs-vecsOutput-1].kspIndex];
      if (ioData->output.rom.addStateToKrylov) { 
        outVec += (currentKspCoord*currentKspVec);
      } else {
        outVec = (currentKspCoord*currentKspVec);
      }
      domain->writeVectorToFile(fileName, *(domain->getKrylovStep()), *(domain->getNewtonTag()), outVec);
      cumEnergy += kspIndexedCoords[numVecs-vecsOutput-1].energy; 
      ++(*(domain->getKrylovStep()));
      com->fprintf(stdout, "vecsOutput %d, energy %e, cumEnergy %e, originalIndex %d\n",
                   vecsOutput+1, kspIndexedCoords[numVecs-vecsOutput-1].energy/totalEnergy, cumEnergy/totalEnergy,
                   kspIndexedCoords[numVecs-vecsOutput-1].kspIndex);
      ++vecsOutput;
    }

    *(domain->getNumKrylovVecsOutputPrevNewtonIt()) = vecsOutput;
    delete [] kspIndexedCoords;
  }
}

//----------------------------------------------------------------------------------

template<class VecType>
template<int dim>
void KspBinaryOutput<VecType>::writeKrylovVectors(VecSet<DistSVec<bcomp, dim> >& kspVecs, Vec<bcomp> kspCoords, int numVecs) {
// do nothing
}

//----------------------------------------------------------------------------------

template<class VecType>
template<class Scalar, int dim>
void KspBinaryOutput<VecType>::writeKrylovVectors(VecSet<DistEmbeddedVec<Scalar, dim> >& kspVecs, Vec<Scalar> kspCoords, int numVecs) {
// do nothing
}

