#ifndef _DIST_EDGE_GRAD_H_
#define _DIST_EDGE_GRAD_H_

class IoData;
class Domain;

template<int dim> class EdgeGrad;
template<class Scalar, int dim> class DistSVec;

//------------------------------------------------------------------------------

template<int dim>
class DistEdgeGrad {

  int failSafeNewton;
  int lastConfig;
  int it0;
  int lastIt;
  int numLocSub;
  SubDomain** subDomain;
  EdgeGrad<dim>** subEdgeGrad;

  DistVec<bool> *tag;
  DistVec<bool> *backuptag;

public:

  DistEdgeGrad(IoData&, Domain*);
  ~DistEdgeGrad();

  EdgeGrad<dim>& operator() (int i) const { return *subEdgeGrad[i]; }

  void compute(int, DistSVec<double,3>&);

  void fix(DistSVec<bool,2> &fstag)
  {

#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subEdgeGrad[iSub]->fix(fstag.subData(iSub), tag->subSize(iSub));

  }

  void fix(DistSVec<int,2> &fstag)
  {

#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub)
      subEdgeGrad[iSub]->fix(fstag.subData(iSub), tag->subSize(iSub));

  }

  void resetTag()
  {
  
    *tag = *backuptag; // assigning backuptags to tag

  }

// Included (MB)
  void computeDerivative(int, DistSVec<double,3>&, DistSVec<double,3>&);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <DistEdgeGrad.C>
#endif

#endif
