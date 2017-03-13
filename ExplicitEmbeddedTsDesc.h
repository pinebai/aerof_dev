#ifndef _EXPLICIT_EMBEDDED_TS_DESC_H_
#define _EXPLICIT_EMBEDDED_TS_DESC_H_

#include <EmbeddedTsDesc.h>

#include <IoData.h>

struct DistInfo;

class GeoSource;
template<int dimLS> class LevelSet;
template<class Scalar, int dim> class DistSVec;

//------------------------------------------------------------------------

template<int dim>
class ExplicitEmbeddedTsDesc : public EmbeddedTsDesc<dim> {

 private:
  ExplicitData::Type timeType;

  DistSVec<double,dim> U0;
  DistSVec<double,dim> k1;
  DistSVec<double,dim> k2;
  DistSVec<double,dim> k3;
  DistSVec<double,dim> k4;

  DistVec<double> Phi0;
  DistVec<double> p1;
  DistVec<double> p2;
  DistVec<double> p3;
  DistVec<double> p4;

  DistVec<double>* hh1,*hh2,*hh3,*hh4,*hhorig;

  bool RK4;
  bool FE;

 public:

  ExplicitEmbeddedTsDesc(IoData &, GeoSource &, Domain *);
  ~ExplicitEmbeddedTsDesc();

  int solveNonLinearSystem(DistSVec<double,dim> &U, int);

 private:
  void solveNLSystemOneBlock(DistSVec<double,dim> &U);
//  void solveNLSystemTwoBlocks(DistSVec<double,dim> &U);

// for solving the total system in one block (U and Phi at the same time)
  void commonPart(DistSVec<double,dim> &U); // Common part to the two following functions.
  void solveNLAllFE(DistSVec<double,dim> &U, double t0, DistSVec<double,dim> &Ubc);
  void solveNLAllRK2(DistSVec<double,dim> &U, double t0, DistSVec<double,dim> &Ubc);
  void solveNLAllRK4(DistSVec<double,dim> &U, double t0, DistSVec<double,dim> &UBc);

// for solving the total system in two blocks (U first, Phi second)
//  void solveNLEuler(DistSVec<double,dim> & U);
//  void solveNLEulerRK2(DistSVec<double,dim> &U);
//  void solveNLEulerRK4(DistSVec<double,dim> &U);
//  void solveNLLevelSet(DistSVec<double,dim> &U);
//  void solveNLLevelSetRK2(DistSVec<double,dim> &U);
//  void solveNLLevelSetRK4(DistSVec<double,dim> &U);


  void computeRKUpdate(DistSVec<double,dim>& Ulocal,
                       DistSVec<double,dim>& dU, int it);

  void computeRKUpdateHH(DistSVec<double,dim>& Ulocal,
			 DistVec<double>& dHH) ;
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ExplicitEmbeddedTsDesc.C>
#endif

#endif

