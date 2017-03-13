#ifndef _EDGE_GALERKIN_H_
#define _EDGE_GALERKIN_H_

#include <IoData.h>
#include <DistVector.h>

class RefVal;
class VarFcn;
class ViscoFcn;
class ThermalCondFcn;
class Domain;
class SubDomain;
class Timer;

template<class T> class CommPattern;
template<int dim> class DistBcData;

//------------------------------------------------------------------------------

class EdgeGalerkin {

  int lastConfig;

  int numLocSub;

  SubDomain **subDomain;

  DistSVec<double,9> M;

  CommPattern<double> *cp;

  Timer *timer;

public:

  EdgeGalerkin(IoData &, Domain *);
  ~EdgeGalerkin();

  template<int dim>
  void compute(int, RefVal *, VarFcn *, ViscoFcn *, ThermalCondFcn *, DistBcData<dim> &, 
	       DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);

  void computeJacobian();

};

//------------------------------------------------------------------------------

#endif
