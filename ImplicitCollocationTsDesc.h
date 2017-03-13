#ifndef IMPLICIT_COLL_TS_DESC_H_
#define IMPLICIT_COLL_TS_DESC_H_

#include <ImplicitGappyTsDesc.h>
#include <RestrictionMapping.h>
#include <DistLeastSquareSolver.h>
#include <NonlinearRomOnlineIII.h>

#include <memory>
#include <vector>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitCollocationTsDesc : public ImplicitGappyTsDesc<dim> {

protected:
  Vec<double> From;
  Vec<double> rhs;

  void setProblemSize(DistSVec<double, dim> &);

  void solveNewtonSystem(const int &, double &, bool &, DistSVec<double, dim> &, const int & totalTimeSteps = 0);

  double meritFunction(int, DistSVec<double, dim> &, DistSVec<double, dim> &, DistSVec<double, dim> &, double);

  public:
  
  ImplicitCollocationTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitCollocationTsDesc();

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitCollocationTsDesc.C>
#endif

#endif
