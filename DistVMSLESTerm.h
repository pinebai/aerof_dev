#ifndef _DIST_VMS_LES_TERM_H_
#define _DIST_VMS_LES_TERM_H_

#include <DistMacroCell.h>
#include <DistVector.h>
#include <Domain.h>
#include <IoData.h>
#include <VMSLESTerm.h>

class VarFcn;


//------------------------------------------------------------------------

template<int dim>
class DistVMSLESTerm {

 private:

  VarFcn *varFcn;

  int lastConfig;
  int numLocSub;
  int it0;
  int lastIt;
  int scopeWidth;
  int scopeDepth;

  Domain            *domain;
  DistMacroCellSet  *macroCells;
  VMSLESTerm        *vmst;
  bool              **masterFlag;

  DistSVec<double,dim> *VtBar;
  DistSVec<double,1>   *volRatio;

 public:

  DistVMSLESTerm(VarFcn *, IoData &, Domain *);
  ~DistVMSLESTerm();

  void compute(int, DistVec<double> &, DistSVec<double,3> &,
	       DistSVec<double,dim> &, DistSVec<double,dim> &);

  void computeMutOMu(DistVec<double> &, DistSVec<double,3> &,
                     DistSVec<double,dim> &, DistVec<double> &);

// Included (MB)
  void computeDerivative(int, DistVec<double> &, DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &);

};

//------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <DistVMSLESTerm.C>
#endif

#endif

