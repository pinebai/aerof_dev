/***********************************************************************************
*                                                                                  *
* Different explicit time-integration schemes are available:                       *
*  - Forward Euler where both fluid states and level set are advanced together     *
*  - RungeKutta2   where both fluid states and level set are advanced together     *
*                                                          (cf LeTallec's comment) *
*  - RungeKutta2   where fluid states and level set are advanced in a staggered    *
*                  fashion (original idea of Pr Farhat to increase stability,      *
*                  but then accuracy is lost)                                      *
*  - RungeKutta4   similar to staggered RungeKutta2                                *
*                                                                                  *
*                                                                                  *
*                                                                                  *
*                                                                                  *
***********************************************************************************/

#ifndef _EXPLICIT_LEVELSET_TS_DESC_H_
#define _EXPLICIT_LEVELSET_TS_DESC_H_

#include <LevelSetTsDesc.h>

struct DistInfo;

class IoData;
class Domain;
class GeoSource;
template<class Scalar, int dim> class DistSVec;

//------------------------------------------------------------------------

template<int dim, int dimLS>
class ExplicitLevelSetTsDesc : public LevelSetTsDesc<dim,dimLS> {

 private:
  ExplicitData::Type timeType;

  DistSVec<double,dim> U0;
  DistSVec<double,dim> k1;
  DistSVec<double,dim> k2;
  DistSVec<double,dim> k3;
  DistSVec<double,dim> k4;

  DistSVec<double,dimLS> Phi0;
  DistVec<int> fluidId0;
  DistSVec<double,dimLS> p1;
  DistSVec<double,dimLS> p2;
  DistSVec<double,dimLS> p3;
  DistSVec<double,dimLS> p4;

  DistVec<double>* hh1,*hh2,*hh3,*hh4,*hhorig;
// mesh motion modification for RK2
// otherwise equal to U and Phi respectively
  DistSVec<double,dim> ratioTimesU;
  DistSVec<double,dimLS> ratioTimesPhi;

  void computeRKUpdateHH(DistSVec<double,dim>& Ulocal,
                         DistVec<double>& dHH);
 
 public:
  ExplicitLevelSetTsDesc(IoData &, GeoSource &, Domain *);
  ~ExplicitLevelSetTsDesc();

  int solveNonLinearSystem(DistSVec<double,dim> &U, int);

 private:
  void solveNLSystemOneBlock(DistSVec<double,dim> &U);
  int solveNLSystemTwoBlocks(DistSVec<double,dim> &U);

// for solving the total system in one block (U and Phi at the same time)
  void solveNLAllFE(DistSVec<double,dim> &U);
  void solveNLAllRK2(DistSVec<double,dim> &U);
  void solveNLAllRK2bis(DistSVec<double,dim> &U);

// for solving the total system in two blocks (U first, Phi second)
  void solveNLEuler(DistSVec<double,dim> & U);
  void solveNLEulerRK2(DistSVec<double,dim> &U);
  void solveNLEulerRK4(DistSVec<double,dim> &U);
  void solveNLLevelSet(DistSVec<double,dim> &U);
  void solveNLLevelSetRK2(DistSVec<double,dim> &U);
  void solveNLLevelSetRK4(DistSVec<double,dim> &U);


  void computeRKUpdate(DistSVec<double,dim>& Ulocal, 
                       DistSVec<double,dim>& dU, int it);
  void computeRKUpdateLS(DistSVec<double,dimLS>& Philocal, DistVec<int> &localFluidId, DistSVec<double,dimLS>& dPhi, 
                         DistSVec<double,dim>& U);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ExplicitLevelSetTsDesc.C>
#endif

#endif
