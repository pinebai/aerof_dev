#ifndef _DIST_DYNAMIC_LES_TERM_H_
#define _DIST_DYNAMIC_LES_TERM_H_

#include <DistVector.h>
#include <IoData.h>
#include <DynamicLESTerm.h>

class VarFcn;

//------------------------------------------------------------------------

template<int dim>
class DistDynamicLESTerm {

 private:

  VarFcn *varFcn;

  int               numLocSub;
  Domain            *domain;
  double            gam, Rideal;

  DistSVec<double,dim> *VCap;
  DistSVec<double,16> *Mom_Test;
  DistSVec<double,6> *Sij_Test;
  DistVec<double> *modS_Test;
  DistSVec<double,8> *Eng_Test;
  DistSVec<double,2> *Cs;
  DistVec<int> *Ni;

  DynamicLESTerm *dlest;

 public:

  DistDynamicLESTerm(VarFcn *, IoData &, Domain *);
  ~DistDynamicLESTerm();

  void compute(DistVec<double> &, DistBcData<dim> &bcData, DistSVec<double,3> &, 
               DistSVec<double,dim> &, DistSVec<double,dim> &,
               DistVec<GhostPoint<dim>*> *ghostPoints=0, 
					DistLevelSetStructure *LSS=0, bool externalSI=false);

  void computeMutOMu(DistVec<double> &, DistBcData<dim> &bcData, DistSVec<double,3> &, 
                     DistSVec<double,dim> &, DistVec<double> &);
  void computeCsValue(DistVec<double> &, DistBcData<dim> &, DistSVec<double,3> &,
                      DistSVec<double,dim> &, DistVec<double> &);

};

//------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <DistDynamicLESTerm.C>
#endif

#endif
