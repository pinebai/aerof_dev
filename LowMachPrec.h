/****************************************************************************
 *                                                                          *
 *  These classes contain all the information for Low Mach Preconditioning  *
 *                                                                          *
 *  At low Mach number, steady-state solutions are difficult to obtain      *
 *  because of the ill-conditioning of the system due to the characteristic *
 *  velocities u, u-c and u+c.                                              *
 *  Also, for unsteady problems, it was shown that Godunov schemes do not   *
 *  respect the asymptotic behaviour of pressure as M->0. This can be       *
 *  remedied by slight modifications.                                       *
 *  The dissipation matrix of the Roe flux can be modified to obtain the    *
 *  the correct asymptotic behavior. And the wave speeds are more generally *
 *  modified to help convergence to the steady-state.                       *
 *                                                                          *
 *  References are:                                                         *
 *   Viozat ....                                                            *
 *       Implicit Upwind Schemes for Low Mach Number Compressible Flows     *
 *       INRIA Report, 1997                                                 *
 *   Turkel ....                                                            *
 *       Preconditioning Techniques in CFD                                  *
 *       Annual reviews in Fluid Mechanics, 1999, vol 31, pp 385-416        *
 *        ***********************************************************       *
 *                                                                          *
 *   A parameter beta is introduced which corresponds to a                  *
 *   characteristic mach number of the flow. More general formula is        *
 *     beta = min(max(beta_min, k*localMach), beta_max)                     *
 *                                                                          *
 *   Other refinements are possible for laminar or turbulent flows          *
 *   (cf inverse Reynolds number below)                                     *
 *                                                                          *
 *        ***********************************************************       *
 *   One virtual class LowMachPrec contains all the necessary constants     *
 *   for both the spatial and the time preconditioners.                     *
 *   A class SpatialLowMachPrec is specific to the spatial preconditioner   *
 *   A class TimeLowMachPrec    is specific to the time    preconditioner   *
 *                                                                          *
 *                                                                          *
 *                                                                          *
 ****************************************************************************
 */


#ifndef _LOW_MACH_PREC_H
#define _LOW_MACH_PREC_H

#include <cmath>

#include <IoData.h>
//------------------------------------------------------------------------------

class LowMachPrec {

//members
protected:
  int prec;
  bool timeac;
  double minMach;
  double maxMach;
  double slope;
  double betaviscous;
  double shockreducer;

//methods
protected:
  virtual void defineLowMachPrecType(IoData &iod){
    prec = 0;
    timeac = false;
    if(iod.problem.type[ProblemData::UNSTEADY]) timeac = true;

    if(iod.problem.prec == ProblemData::PRECONDITIONED){
      prec = 1; //at least spatial preconditioning
      if(iod.ts.prec == TsData::PREC){
        prec = 2;
        if(iod.problem.alltype != ProblemData::_SHAPE_OPTIMIZATION_ &&
           iod.problem.alltype != ProblemData::_SENSITIVITY_ANALYSIS_ && //TODO CHECK if neccessary
           iod.problem.alltype != ProblemData::_ROM_SHAPE_OPTIMIZATION_ &&
	       iod.problem.alltype != ProblemData::_STEADY_ &&
           iod.problem.alltype != ProblemData::_STEADY_AEROELASTIC_ &&
           iod.problem.alltype != ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_ &&
           iod.problem.alltype != ProblemData::_STEADY_THERMO_ &&
           iod.problem.alltype != ProblemData::_STEADY_AEROTHERMOELASTIC_){
          if(iod.ts.type != TsData::IMPLICIT)
            prec = 1; // no time precontioning if not steady and while using explicit time integration
        }
      }
    }
  }

  virtual void setupDefaultConstants(){
    minMach = 1.0; maxMach = 1.0; slope = 0.0;
    betaviscous = 0.0; shockreducer = 0.0;
  }

  virtual void setupIodConstants(IoData &iod){
    minMach      = iod.prec.mach;
    maxMach      = iod.prec.cmach;
    slope        = iod.prec.k;
    betaviscous  = iod.prec.betav;
    shockreducer = iod.prec.shockreducer;
  }

public:

  LowMachPrec() {
    prec = 0;
    setupDefaultConstants();
  }
  LowMachPrec(IoData &iod){
    defineLowMachPrecType(iod);
    setupDefaultConstants();
  }
  virtual ~LowMachPrec() {}; //destructor of base class should always be virtual

  virtual void setup(IoData &iod) {
    defineLowMachPrecType(iod);
    setupDefaultConstants();
  }

  virtual int    getPrecTag() const        { return prec;         }
  virtual double getMinMach() const        { return minMach;      }
  virtual double getCutOffMach() const     { return maxMach;      }
  virtual double getSlope() const          { return slope;        }
  virtual double getViscousRatio() const   { return betaviscous;  }
  virtual double getShockParameter() const { return shockreducer; }

  virtual double getBeta(double locMach,bool checktac = false) const = 0;


};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

class SpatialLowMachPrec : public LowMachPrec {

public:
  SpatialLowMachPrec() : LowMachPrec() {};
  SpatialLowMachPrec(IoData &iod) : LowMachPrec(iod) {if (prec>0) setupIodConstants(iod); };
  ~SpatialLowMachPrec() {};

  void setup(IoData &iod) {
    defineLowMachPrecType(iod);
    if(prec==0) setupDefaultConstants();
    else        setupIodConstants(iod);
  }

  double getBeta(double locMach, bool checktac = false) const
    {return fmin(fmax(slope*locMach, minMach),maxMach);}

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

class TimeLowMachPrec : public LowMachPrec {

public:
  TimeLowMachPrec() : LowMachPrec() {};
  TimeLowMachPrec(IoData &iod) : LowMachPrec(iod) {if (prec==2) setupIodConstants(iod); };
  ~TimeLowMachPrec() {};

  void setup(IoData &iod) {
    defineLowMachPrecType(iod);
    if(prec < 2) setupDefaultConstants();
    else         setupIodConstants(iod);
  }

  bool timePreconditioner() const { return prec==2; }

//  double getBeta(double locMach) const
//    {return fmin(fmax(slope*locMach, minMach),maxMach);}

  double getBeta(double locMach, bool checktac = false) const
    { double outbeta = (checktac && timeac) ? 1.0 : fmin(fmax(slope*locMach, minMach),maxMach);
      return outbeta; }

  double getBeta(double locMach, double irey) const
    { double outbeta = fmax(slope*locMach, minMach);
      return fmin((1.0+sqrt(irey))*outbeta,maxMach); }

  double getdBeta(double locMach, double dLocMach, bool checktac = false) const {
    double locbeta = getBeta(locMach,checktac);
    if (locbeta == maxMach)
      return 0.0;
    else
      if (fmax(slope*locMach, minMach) == minMach)
        return 0.0;
      else
        return slope*dLocMach;
  }

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#endif
