#ifndef _EXPLICIT_TS_DESC_H_
#define _EXPLICIT_TS_DESC_H_

#include <IoData.h>
#include <TsDesc.h>

class GeoSource;
class Domain;

template<class Scalar> class DistVec;
template<class Scalar, int dim> class DistSVec;

//-------------------------------------------------------------------------

template<int dim>
class ExplicitTsDesc : public TsDesc<dim> {

private:
  LiquidModelData::YesNo check;
  DistSVec<double,dim> U0;
  DistSVec<double,dim> k1;
  DistSVec<double,dim> k2;
  DistSVec<double,dim> k3;
  DistSVec<double,dim> k4;
  DistSVec<double,dim> Ubc;
  DistSVec<double,dim> ratioTimesU;
  bool RK4;
  bool FE; //forward Euler

public:

  ExplicitTsDesc(IoData&, GeoSource&, Domain*);
  ~ExplicitTsDesc();

  virtual int solveNonLinearSystem(DistSVec<double,dim>&, int);

private:

  virtual void computeRKFourthOrder(DistSVec<double,dim>&);
  virtual void computeRKSecondOrder(DistSVec<double,dim>&);
  virtual void computeForwardEuler(DistSVec<double,dim>&);

  virtual void computeRKUpdate(DistSVec<double,dim>&, DistSVec<double,dim>&);

};

//-------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ExplicitTsDesc.C>
#endif

#endif
