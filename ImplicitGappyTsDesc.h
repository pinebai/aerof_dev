#ifndef IMPLICIT_GAPPY_TS_DESC_H_
#define IMPLICIT_GAPPY_TS_DESC_H_

#include <ImplicitRomTsDesc.h>
#include <RestrictionMapping.h>
#include <DistLeastSquareSolver.h>
#include <NonlinearRomOnlineIII.h>

#include <memory>
#include <vector>

//------------------------------------------------------------------------------

template<int dim>
class ImplicitGappyTsDesc : public ImplicitRomTsDesc<dim> {

protected:

#if (__cplusplus >= 201103L)
  std::unique_ptr< VecSet < DistSVec<double, dim> > > AJRestrict;
  std::unique_ptr< DistSVec<double, dim> > ResRestrict;
#else
  std::auto_ptr< VecSet < DistSVec<double, dim> > > AJRestrict;
  std::auto_ptr< DistSVec<double, dim> > ResRestrict;
#endif

  void setProblemSize(DistSVec<double, dim> &) = 0;

  void solveNewtonSystem(const int &, double &, bool &, DistSVec<double, dim> &, const int & totalTimeSteps = 0) = 0;

  virtual void computeFullResidual(int, DistSVec<double, dim> &, bool applyWeighting = false, DistSVec<double, dim> *R = NULL, bool includeHomotopy = true);

  virtual void computeAJ(int, DistSVec<double, dim> &, bool applyWeighting = false, DistSVec<double, dim> *R = NULL);
  virtual bool breakloop1(const bool);
  virtual bool breakloop2(const bool);

  double computeGappyResidualNorm(DistSVec<double,dim> &);
  void setReferenceResidual();
  void deleteRestrictedQuantities();

  double *jactmp, *column;

public:
  
  bool checkForLastIteration(IoData &, int, double, double, DistSVec<double,dim> &);
  void monitorInitialState(int, DistSVec<double,dim> &);
  bool monitorConvergence(int, DistSVec<double,dim> &);
  void formInterpolatedInitialCondition(DistSVec<double,dim> *U, IoData &iod)  {
    this->rom->formInterpolatedInitialCondition(U, iod);
  }

  ImplicitGappyTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitGappyTsDesc();

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitGappyTsDesc.C>
#endif

#endif
