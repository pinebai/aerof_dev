#ifndef IMPLICIT_GNAT_TS_DESC_H_
#define IMPLICIT_GNAT_TS_DESC_H_

#include <ImplicitGappyTsDesc.h>
#include <RestrictionMapping.h>
#include <DistLeastSquareSolver.h>
#include <NonlinearRomOnlineIII.h>

#include <memory>
#include <vector>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitGnatTsDesc : public ImplicitGappyTsDesc<dim> {

protected:

  int nPodJac;	//nPodJac specified under ROB{ NumROB2 }
  int numResJacMat ;	// number of matrices for A and B (1 if they use the same)

  DistLeastSquareSolver leastSquaresSolver;

  void setProblemSize(DistSVec<double, dim> &);

  void solveNewtonSystem(const int &, double &, bool &, DistSVec<double, dim> &, const int & totalTimeSteps = 0);

  double meritFunction(int, DistSVec<double, dim> &, DistSVec<double, dim> &, DistSVec<double, dim> &, double);

public:
  
  ImplicitGnatTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitGnatTsDesc();

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitGnatTsDesc.C>
#endif

#endif
