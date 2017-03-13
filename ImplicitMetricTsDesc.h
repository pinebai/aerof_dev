#ifndef IMPLICIT_METRIC_TS_DESC_H_
#define IMPLICIT_METRIC_TS_DESC_H_

#include <ImplicitGappyTsDesc.h>
#include <RestrictionMapping.h>
#include <DistLeastSquareSolver.h>
#include <NonlinearRomOnlineIII.h>

#include <memory>
#include <vector>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitMetricTsDesc : public ImplicitGappyTsDesc<dim> {

protected:
  Vec<double> From;
  Vec<double> rhs;

  void setProblemSize(DistSVec<double, dim> &);

  void solveNewtonSystem(const int &, double &, bool &, DistSVec<double, dim> &, const int & totalTimeSteps = 0);

  double meritFunction(int, DistSVec<double, dim> &, DistSVec<double, dim> &, DistSVec<double, dim> &, double);

  public:
  
  ImplicitMetricTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitMetricTsDesc();

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitMetricTsDesc.C>
#endif

#endif
