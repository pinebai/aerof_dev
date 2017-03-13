#ifndef _DIST_EXTRAPOLATION_H_
#define _DIST_EXTRAPOLATION_H_

class IoData;
class Domain;
class VarFcn;

template<int dim> class Extrapolation;

//------------------------------------------------------------------------------

template<int dim>
class DistExtrapolation {

  int lastConfig;
  int it0;
  int lastIt;
  int numLocSub;
  SubDomain** subDomain;
  Extrapolation<dim>** subExtrapolation;
  Communicator *myCom;

public:

  DistExtrapolation(IoData&, Domain*, VarFcn*);
  ~DistExtrapolation();

  Extrapolation<dim>& operator() (int i) const { return *subExtrapolation[i]; }

  void compute(int, DistVec<Vec3D> &, DistSVec<double,3>&);

// Included (MB)
   void computeDerivative(int, DistVec<Vec3D> &, DistSVec<double,3>&);
 
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <DistExtrapolation.C>
#endif


#endif
