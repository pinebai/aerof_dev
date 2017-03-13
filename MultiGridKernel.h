/* MultiGridKernel.h

 */

#pragma once

#include <DistMatrix.h>
#include <DistGeoState.h>
#include <MultiGridLevel.h>
#include <KspSolver.h>
#include <DistEmbeddedVector.h>
#include <MatVecProd.h>
#include <MultiGridOperator.h>
#include <DistTimeState.h>
#include <SpaceOperator.h>

class DistGeoState;
template<class Scalar,int dim> class DistSVec;

template<class Scalar>
class MultiGridKernel {

 public:

  MultiGridKernel(Domain *dom, DistGeoState& distGeoState, 
                  IoData&,int num_levels);

  ~MultiGridKernel();

  void setParameters(int v1, int v2, int
                     fine_sweeps, double relax, int do_out);

  void initialize(int dim, int neq1,int neq2);

  void setupAlgebraic();

  bool isInitialized() const { return initialized; }

  template<class Scalar2, int dim>
  void Restrict(int coarseLvl, DistSVec<Scalar2,dim>& fine, 
                DistSVec<Scalar2,dim>& coarse, bool average = true, 
		bool = false);

  template<class Scalar2, int dim>
  void ExtrapolateProlongation(int coarseLvl, DistSVec<Scalar2,dim>& coarseOld, 
			       DistSVec<Scalar2,dim>& coarse,	      
			       class DistLevelSetStructure* coarselss = NULL,
			       class DistLevelSetStructure* finelss = NULL);

  template<class Scalar2, int dim>
  void Prolong(int coarseLvl, DistSVec<Scalar2,dim>& coarseOld, 
               DistSVec<Scalar2,dim>& coarse, DistSVec<Scalar2,dim>& fine,
	       DistSVec<Scalar2,dim>& fine_ref,
               double relax,VarFcn* varFcn, 
	       class DistLevelSetStructure* coarselss = NULL,
	       class DistLevelSetStructure* finelss = NULL);

  int numLevels() const { return num_levels; }
 
  MultiGridLevel<Scalar> * getLevel(int i) { return multiGridLevels[i]; }

  void setUseVolumeWeightedAverage(bool b);

  template<class Scalar2, int dim>
  void fixNegativeValues(int,DistSVec<Scalar2,dim>& V, 
                         DistSVec<Scalar2,dim>& U, 
                         DistSVec<Scalar2,dim>& dx, 
                         DistSVec<Scalar2,dim>& f, 
                         DistSVec<Scalar2,dim>& forig, 
                         VarFcn*,
                         MultiGridOperator<Scalar2,dim>*);
                         
  template<class Scalar2, int dim>
  void applyFixes(int,DistSVec<Scalar2,dim>& f); 

 private:
 
  void
  setupFixes(IoData& ioData,int lvl,DistSVec<Scalar,3>& X0);

  double ref_length;

  bool isGeometric;
 
  bool coarsen4to1;

  int nSmooth1,nSmooth2;

  double relaxationFactor;
  
  const int num_levels, agglom_size, numLocSub;

  double beta; 
 
  double prolong_relax_factor;
  double restrict_relax_factor;

  bool smoothWithGMRES;

  MultiGridLevel<Scalar> ** multiGridLevels;
 
  DistGeoState& geoState;

  IoData& ioData;

  Domain* domain;

  int output;

  int fine_sweeps;

  bool initialized;

  std::set<int>** fixLocations;

  const char* agglomerationFile;

  double turbRelaxCutoff;
 
};
