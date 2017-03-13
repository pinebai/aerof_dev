#ifndef _FLUID_SELECTOR_H_
#define _FLUID_SELECTOR_H_

#include "DistVector.h"
#include "Vector.h"
#include "IoData.h"
#include "ProgrammedBurn.h"

//#include "Domain.h"
//#include "LevelSet/LevelSetStructure.h"
#include <map>
//#include <cassert>
//#define NDEBUG // if commented, assert statements are evaluated
using std::map;
using std::pair;

class Domain;
class DistLevelSetStructure;

//--------------------------------------------------------------------------
//
// This algorithm class allows to determine in which fluid a node is,
// by determining the fluid identification from the different 
// available level sets.
//
// The convention is as follows. For n different fluids, only
// n-1 (=dim) level-sets are necessary.
// Phi[i] is positive where fluid i+1 (=varFcn[i+1]) is found. 
// fluid 0 (=varFcn[0]) is found where all levelsets are negative.
//--------------------------------------------------------------------------

class FluidSelector {

  IoData* iodp;
  int numPhases;
  Domain *domain;
  map<int,pair<int,int> > ls2phases; //maps each levelset to the pair of fluid Ids for phi<0 and >0.
                                     //used only for multi-phase fluid-structure interactions right now.

  ProgrammedBurn* programmedBurn;

  // True if we own the fluid id variable.
  bool ownsData;

public:
  DistVec<int> *fluidId;
  DistVec<int> *fluidIdn;
  DistVec<int> *fluidIdnm1;
  DistVec<int> *fluidIdnm2;

protected:
  void setupFluidIdVolumesInitialConditions(IoData &ioData);
  void setupFluidIdOneDimensionalSolution(IoData &ioData, DistSVec<double,3> &X);
  void setupFluidIdMultiFluidInitialConditions(IoData &ioData, DistSVec<double,3> &X);

public:

  FluidSelector(const int nPhases, IoData &ioData, Domain *domain = 0);

  // Create a temporary fluid selector from the pointer to the node tag.
  FluidSelector(DistVec<int>& nodeTag,Domain* dom);

  ~FluidSelector();

  void attachProgrammedBurn(ProgrammedBurn* p) { programmedBurn = p; }

  ProgrammedBurn* getProgrammedBurn() const { return programmedBurn; }

//------------------------------------------------------------------------------

  int getNumOfPhases() {return numPhases;}
  void initializeFluidIds(DistVec<int> &fsId, DistSVec<double,3> &X, IoData &ioData);
  void getFluidId(DistVec<double> &Phi);
  void getFluidId(int &tag, double phi){ tag = (phi<0.0) ? 0 : 1; }
  void getFluidId(int &tag, double *phi);
  void getFluidId(Vec<int> &tag, Vec<double> &phi);
  void getFluidId(DistVec<int> &Tag, DistVec<double> &Phi);
  int getLevelSetDim(int fluidId1, int fluidId2, int node1 = 0, int node2 = 0);
  void printFluidId();
  void update(){
    if(fluidIdnm2) *fluidIdnm2 = *fluidIdnm1;
    if(fluidIdnm1) *fluidIdnm1 = *fluidIdn;
    *fluidIdn = *fluidId;
  }

//------------------------------------------------------------------------------

  template<int dim>
  void initializeFluidIds(DistSVec<double,dim> &Phin, DistSVec<double,dim> &Phinm1, DistSVec<double,dim> &Phinm2);

  template<int dim>
  void reinitializeFluidIds(DistVec<int> &fsId, DistSVec<double,dim> &Phin);
  
  template<int dim>
  void reinitializeFluidIdsWithCracking(DistVec<int> &fsId, DistSVec<double,dim> &Phin);

  template<int dim> /*this dim is actually dimLS*/
  void updateFluidIdFS(DistLevelSetStructure *distLSS, DistSVec<double,dim> &PhiV);

  template<int dim> /*this dim is actually dimLS*/
  void updateFluidIdFS2(DistLevelSetStructure *distLSS, DistSVec<double,dim> &PhiV, DistSVec<bool,4> &poll);

  template<int dim> /*this dim is actually dimLS*/
  void updateFluidIdFF(DistLevelSetStructure *distLSS, DistSVec<double,dim> &Phi);

  template<int dim> /*this dim is actually dimLS*/
  void updateFluidIdFF2(DistLevelSetStructure *distLSS, DistSVec<double,dim> &Phi);

  template<int dim>
  void getFluidId(DistVec<int> &Tag, DistSVec<double,dim> &Phi, DistVec<int>* fsId=0);

  template<int dim>
  void getFluidId(Vec<int> &tag, SVec<double,dim> &phi);

  template<int dim>
  void getFluidId(DistSVec<double,dim> &Phi);

  //template<int dim>
  void getFluidId(class TriangulatedInterface*);
  
  template<int dim>
  void getFluidId(int &tag, double *phi);

  template<int dim>
  void checkLSConsistency(DistSVec<double,dim> &Phi);

  void writeToDisk(const char* fn);
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <FluidSelector.C>
#endif

#endif
