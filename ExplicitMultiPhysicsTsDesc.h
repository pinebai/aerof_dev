#ifndef _EXPLICIT_MULTIPHYSICS_TS_DESC_H_
#define _EXPLICIT_MULTIPHYSICS_TS_DESC_H_

#include <MultiPhysicsTsDesc.h>
#include <DistVector.h>
#include <DebugTools.h>

class IoData;
class Domain;
class GeoSource;
//template<class Scalar, int dim> class DistVec;

//------------------------------------------------------------------------

template<int dim, int dimLS>
class ExplicitMultiPhysicsTsDesc : public MultiPhysicsTsDesc<dim,dimLS> {

 private:
  ExplicitData::Type timeType;

  DistSVec<double,dim> U0;
  DistSVec<double,dim> k1;
  DistSVec<double,dim> k2;
  DistSVec<double,dim> k3;
  DistSVec<double,dim> k4;
  DistSVec<double,dim> Ubc;

  DistSVec<double,dimLS> Phi0;
  DistVec<int> fluidId0; 
  DistSVec<double,dimLS> p1;
  DistSVec<double,dimLS> p2;
  DistSVec<double,dimLS> p3;
  DistSVec<double,dimLS> p4;

  DistVec<double>* hh1,*hh2,*hh3,*hh4,*hhorig;

 public:
  ExplicitMultiPhysicsTsDesc(IoData &, GeoSource &, Domain *);
  ~ExplicitMultiPhysicsTsDesc();

  int solveNonLinearSystem(DistSVec<double,dim> &U, int);

 private:
//  void solveNLSystemOneBlock(DistSVec<double,dim> &U);
  void solveNLSystemTwoBlocks(DistSVec<double,dim> &U);

  void recomputeIntersections();
  void updateFluidIdFS(DistSVec<double,dim> &U);
  void updatePhaseChangeFS(DistSVec<double,dim> &U);
  void populateGhostPointsForNavierStokes(DistSVec<double,dim> &U);
  void solveNLNavierStokes(DistSVec<double,dim> &U);
  void solveNLLevelSet(DistSVec<double,dim> &U);
  void updatePhaseChangeFF(DistSVec<double,dim> &U);

  void solveNLNavierStokesFE(DistSVec<double,dim> &U);
  void solveNLNavierStokesRK2(DistSVec<double,dim> &U);
  void solveNLNavierStokesRK4(DistSVec<double,dim> &U);
  void solveNLLevelSetFE(DistSVec<double,dim> &U);
  void solveNLLevelSetRK2(DistSVec<double,dim> &U);
  void solveNLLevelSetRK4(DistSVec<double,dim> &U);

  void computeRKUpdate(DistSVec<double,dim>& Ulocal, DistSVec<double,dim>& dU, int it);
  void computeRKUpdateLS(DistSVec<double,dimLS>& Philocal, DistVec<int> &localFluidId, 
                         DistSVec<double,dimLS>& dPhi, DistSVec<double,dim>& U);

  void computeRKUpdateHH(DistSVec<double,dim>& Ulocal,
			 DistVec<double>& dHH) ;
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ExplicitMultiPhysicsTsDesc.C>
#endif

#endif

