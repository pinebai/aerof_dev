#ifndef _KSP_CONV_CRITERION_H_
#define _KSP_CONV_CRITERION_H_

#include <IoData.h>

//------------------------------------------------------------------------------

class KspConvCriterion {

  KspData::EpsFormula typeEpsFormula;

  int nlits;

  double eps0;
  double fnormPrev;
  double epsPrev;

public:

  KspConvCriterion(KspData &);
  ~KspConvCriterion() {}

  double compute(int, int, double);
  double computeEpsilonEisenstadt(double);

};

//------------------------------------------------------------------------------

#endif
