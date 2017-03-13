#ifndef _TIME_DATA_H_
#define _TIME_DATA_H_

#include <IoData.h>

template<class Scalar> class DistVec;
template<class Scalar, int dim> class DistSVec;

//------------------------------------------------------------------------------

class TimeData {

public:

  ImplicitData::Type typeIntegrator;
  ImplicitData::Startup typeStartup;
  TsData::TypeTimeStep typeTimeStep;

  double dt_imposed;
  double dt_n;
  double dt_nm1;
  double dt_nm2;

  int outputNewtonTag;
  int outputNewtonStateStep;
  int outputNewtonResidualStep;
  int outputKrylovStep;

  double dtau_switch;

  double tau_n;
  double tau_nm1;
  double alpha_np1;
  double alpha_n;
  double alpha_nm1;
  double alpha_nm2;

  double errorTol;

  bool exist_nm1;
  bool exist_nm2;
  bool use_nm1;
  bool use_nm2;

  bool use_freq;
  bool use_modal;

  int descriptor_form;

public:

  TimeData(IoData &);
  ~TimeData() {}

  void copy(TimeData& oth);

  void update();
  void computeCoefficients(DistVec<double> &, double);
  void computeVelocities(DGCLData::Velocities, DistSVec<double,3> &,
			 DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,3> &);

  double getTauN() const { return tau_n; }
  int getNewtonTag() const  { return outputNewtonTag; }
  int getNewtonStateStep() const  { return outputNewtonStateStep; }
  int getNewtonResidualStep() const  { return outputNewtonResidualStep; }
  int getKrylovStep() const  { return outputKrylovStep; }

// Included
  void rstVar(IoData &);

};

//------------------------------------------------------------------------------

#endif
