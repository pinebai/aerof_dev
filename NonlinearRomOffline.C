#include <cstdio>
#include <cmath>
#include <sys/time.h>
#include <algorithm>
#include <cstdlib>
#include <Domain.h>
#include <NonlinearRomOffline.h>
#include <IoData.h>
#include <NonlinearRom.h>
#include <NonlinearRomDatabaseConstruction.h>
#include <DistVector.h>
#include <VectorSet.h>
#include <MatVecProd.h>
#include <Timer.h>
#include <DistBcData.h>
#include <DistGeoState.h>
#include <DistTimeState.h>
#include <PostOperator.h>
#include <ParallelRom.h>
#include <EmbeddedAlternatingLeastSquare.h> //<! Lei Lei, 02/13/2016


template <int dim>
NonlinearRomOfflineSolver<dim>::NonlinearRomOfflineSolver(Communicator *_com, IoData &_ioData, Domain &dom, GeoSource &_geoSource) :
          domain(dom), Xref(dom.getNodeDistInfo()), controlVol(dom.getNodeDistInfo()), geoSource(_geoSource) {

 com = _com;
 ioData = &_ioData;
 geoState = 0;

}

//-------------------------------------------------------------------------------

template <int dim>
void NonlinearRomOfflineSolver<dim>::solve()  {

 // set up Timer
 Timer *modalTimer = domain.getTimer();
 double t0;


 if (ioData->problem.alltype == ProblemData::_NONLINEAR_ROM_OFFLINE_) {
     //Lei Lei, 04 July 2016, if (embedded) do ALS
     if(ioData->problem.framework == ProblemData::EMBEDDED) {
         com->fprintf(stdout, "reached  Embedded ALS ROM\n");
         geoState = new DistGeoState(*ioData, &domain);
         EmbeddedAlternatingLeastSquare<dim> embeddedAlternatingLeastSquare(com, *ioData, domain/*, geoState * not sure why it is needed */);
         // offline phase, construct the basis
         //embeddedAlternatingLeastSquare.readSnapshotsFilesHelper("mask");
         //embeddedAlternatingLeastSquare.readSnapshotsFilesHelper("state");
         //com->fprintf(stderr, "... read state file successfully\n");
         //testing
         //embeddedAlternatingLeastSquare.testingSnapshotIO();
         //embeddedAlternatingLeastSquare.testingInitialization();
         //embeddedAlternatingLeastSquare.testingALS();
         //embeddedAlternatingLeastSquare.ReducedOrderBasisConstruction();
         embeddedAlternatingLeastSquare.constructDatabase();
         // early termination
         modalTimer->print(this->domain.getStrTimer());
         if (geoState) delete geoState;
         return;
     }
   // create file system, cluster snapshots, construct ROBs
   if (ioData->romDatabase.nClusters > 0) {
     t0 = modalTimer->getTime();
     NonlinearRomDatabaseConstruction<dim> rom(com,*ioData,domain,geoSource);
     rom.constructDatabase();
     modalTimer->addPodConstrTime(t0);
   }
   // perform GNAT preprocessing (probably for snapshot collection method 0)
   const char *gappyPrefix = ioData->romDatabase.files.gappyPrefix;
   const char *sampledStateBasisName = ioData->romDatabase.files.sampledStateBasisName;
   if (strcmp(gappyPrefix,"")!=0 || strcmp(sampledStateBasisName,"")!=0) {
     geoState = new DistGeoState(*ioData, &domain);
     geoState->setup1(ioData->input.positions, &Xref, &controlVol);
     GappyPreprocessing<dim> gappy(com,*ioData,domain,geoState);
     gappy.buildReducedModel();
   }
   const char *surfacePrefix = ioData->romDatabase.files.surfacePrefix;
   const char *surfaceStateBasisName = ioData->romDatabase.files.surfaceStateBasisName;
   if (strcmp(surfacePrefix,"")!=0 || strcmp(surfaceStateBasisName,"")!=0) {
     geoState = new DistGeoState(*ioData, &domain);
     geoState->setup1(ioData->input.positions, &Xref, &controlVol);
     SurfMeshGen<dim> surfMeshGen(com,*ioData,domain,geoState);
     surfMeshGen.buildReducedModel();
   }
 }
 else if (ioData->problem.alltype == ProblemData::_NONLINEAR_ROM_PREPROCESSING_){
   geoState = new DistGeoState(*ioData, &domain);
   geoState->setup1(ioData->input.positions, &Xref, &controlVol);
   GappyPreprocessing<dim> gappy(com,*ioData,domain,geoState);
   gappy.buildReducedModel();
 }
 else if (ioData->problem.alltype == ProblemData::_SURFACE_MESH_CONSTRUCTION_){
   geoState = new DistGeoState(*ioData, &domain);
   geoState->setup1(tInput->positions, &Xref, &controlVol);
   SurfMeshGen<dim> surfMeshBuilder(com,*ioData,domain,geoState);
   surfMeshBuilder.buildReducedModel();
 }
 else if (ioData->problem.alltype == ProblemData::_SAMPLE_MESH_SHAPE_CHANGE_){
   //KTC: CHANGE!!!
   geoState = new DistGeoState(*ioData, &domain);
   geoState->setup1(tInput->positions, &Xref, &controlVol);
   ReducedMeshShapeChanger<dim> reducedMeshShapeChanger(com,*ioData,domain,geoState);
   reducedMeshShapeChanger.buildReducedModel();
 } else if (ioData->problem.alltype == ProblemData::_EMBEDDED_ALS_ROM_) {
     // Lei Lei, 02/13/2016
     com->fprintf(stdout, "reached  Embedded ALS ROM\n");
     com->fprintf(stdout, "skipped\n");
     geoState = new DistGeoState(*ioData, &domain);
     EmbeddedAlternatingLeastSquare<dim> embeddedAlternatingLeastSquare(com, *ioData, domain/*, geoState * not sure why it is needed */);
     // offline phase, construct the basis
     embeddedAlternatingLeastSquare.readSnapshotsFilesHelper("mask", false);
     embeddedAlternatingLeastSquare.readSnapshotsFilesHelper("state", false);
     com->fprintf(stderr, "... read state file successfully\n");
     //testing
     //embeddedAlternatingLeastSquare.testingSnapshotIO();
     //embeddedAlternatingLeastSquare.ReducedOrderBasisConstruction(40);
     embeddedAlternatingLeastSquare.ReducedOrderBasisConstructionTesting(40);

 }

    modalTimer->print(this->domain.getStrTimer());

 if (geoState) delete geoState;
}

//---------------------------------------------------------------------------------------
