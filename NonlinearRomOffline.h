#ifndef _NLROM_H_
#define _NLROM_H_

class IoData;
class Domain;
class GeoSource;

#include <DistVector.h>
#include <VectorSet.h>
#include <SpaceOperator.h>
#include <Communicator.h>
#include <TsInput.h>
#include <TsOutput.h>
#include <TsRestart.h>
#include <GappyPreprocessing.h>
#include <SurfMeshGen.h>
#include <ReducedMeshShapeChanger.h>
#include <ParallelRom.h> 
#include <time.h> 

//-----------------------------------------------------------------------------------

template <int dim>
class NonlinearRomOfflineSolver {

    Communicator *com;
    Domain &domain;
    IoData *ioData;
    GeoSource &geoSource;
    DistBcData<dim> *bcData;
    DistGeoState *geoState;
    DistVec<double> controlVol;
    DistSVec<double,3> Xref;
    TsInput *tInput;

  public:
    NonlinearRomOfflineSolver(Communicator *, IoData &, Domain &, GeoSource &);
    void solve();
};

#include "NonlinearRomOffline.C"
#endif
