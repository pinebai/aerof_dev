#ifndef _HEAT_TRANSFER_HANDLER_H_
#define _HEAT_TRANSFER_HANDLER_H_

#include <DistVector.h>

class IoData;
class MatchNodeSet;
class Domain;
class StructExc;

template<int dim> class PostOperator;

//------------------------------------------------------------------------------

class HeatTransferHandler {

  bool steady;
  int it0;

  DistVec<double> P;

  StructExc* strExc;

  Domain* domain;
  Communicator* com;

public:

  HeatTransferHandler(IoData&, MatchNodeSet**, Domain*);
  ~HeatTransferHandler() {}

  void setup(int*, double*);
  double update(bool*, int, DistVec<double>&);
  double updateStep1(bool*, int, DistVec<double>&);
  double updateStep2(bool*, int, DistVec<double>&);

  template<int dim>
  void updateOutputToStructure(double, double, PostOperator<dim>*, 
			       DistSVec<double,3>&, DistSVec<double,dim>&);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <HeatTransferHandler.C>
#endif

#endif
