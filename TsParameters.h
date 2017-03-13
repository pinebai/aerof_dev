#ifndef _TS_PARAMETERS_H_
#define _TS_PARAMETERS_H_

#include <cstdio>
#include <ErrorHandler.h>

class IoData;

//------------------------------------------------------------------------------

class TsParameters {

  int cfllaw;

  double cfl0;
  double cflCoef1;
  double cflCoef2;
  double cflMax;
  double cflMin;
  double ser;
 
  double angle_growth;
  double angle_zero;
  int dft_history;
  int dft_freqcutoff;
  double dft_growth;
  int fixedunsteady_counter;
 
  double* reshistory;
  complex<double>* dft;

  ErrorHandler* errorHandler;

  FILE *cfl_output;
  bool nonlinearRomPostpro;

public:

  int maxIts;
  int resType;
  double eps;
  double epsabs;
  double maxTime;
  double dt_imposed;
  double cfl;
  double residual;
  double dualtimecfl;

  char *output;

  int checksol;
  int checklinsolve;
  int checkriemann;
  int checklargevelocity;
  int rapidpchangecutoff;
  int rapiddchangecutoff;
  int checkpclipping;
  int checkrhoclipping;

  bool allowstop;

public:

  TsParameters(IoData &);
  ~TsParameters();

  void computeCflNumber(int, double, double);
  void resolveErrors();
  double getCflMinOverCfl0(){return (cflMin/cfl0);}
  void assignErrorHandler(ErrorHandler* in){errorHandler = in; }

// Included (MB)
  void rstVar(IoData &);

};

//------------------------------------------------------------------------------

#endif
