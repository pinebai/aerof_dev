#ifndef _FLUX_FCN_DESC_JWL_H_
#define _FLUX_FCN_DESC_JWL_H_

#include <FluxFcnDesc.h>

class IoData;
#include "VarFcnJwl.h"

//------------------------------------------------------------------------------

class FluxFcnJwlFDJacRoeEuler3D : public FluxFcnFDJacRoeEuler3D {

public:

  FluxFcnJwlFDJacRoeEuler3D(double gg, IoData &ioData, VarFcnJwl *varFcnJwl, Type tp = CONSERVATIVE) :
    FluxFcnFDJacRoeEuler3D(ioData, gg, varFcnJwl, tp) {}
  ~FluxFcnJwlFDJacRoeEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnJwlApprJacRoeEuler3D : public FluxFcnApprJacRoeEuler3D {

public:

  FluxFcnJwlApprJacRoeEuler3D(int rs, double gg, IoData &ioData, VarFcnJwl *varFcnJwl, Type tp = CONSERVATIVE) : 
    FluxFcnApprJacRoeEuler3D(ioData, rs, gg, varFcnJwl, tp) {}
  ~FluxFcnJwlApprJacRoeEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnJwlExactJacRoeEuler3D : public FluxFcnExactJacRoeEuler3D {

public:

  FluxFcnJwlExactJacRoeEuler3D(double gg, IoData &ioData, VarFcnJwl *varFcnJwl, Type tp = CONSERVATIVE) : 
    FluxFcnExactJacRoeEuler3D(ioData, gg, varFcnJwl, tp) {}
  ~FluxFcnJwlExactJacRoeEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnJwlWallEuler3D : public FluxFcnWallEuler3D {

public:
  FluxFcnJwlWallEuler3D(IoData &ioData, VarFcnJwl *varFcnJwl, Type tp=CONSERVATIVE) :
    FluxFcnWallEuler3D(varFcnJwl, tp) {}
  ~FluxFcnJwlWallEuler3D() { vf = 0; }
  
  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnJwlGhidagliaEuler3D : public FluxFcnGhidagliaEuler3D {

public:

  FluxFcnJwlGhidagliaEuler3D(IoData &ioData, VarFcnJwl *varFcnJwl, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaEuler3D(varFcnJwl, tp) {}
  ~FluxFcnJwlGhidagliaEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnJwlInflowEuler3D : public FluxFcnInflowEuler3D {

public:

  FluxFcnJwlInflowEuler3D(IoData &ioData, VarFcnJwl *varFcnJwl, Type tp = CONSERVATIVE) :
    FluxFcnInflowEuler3D(varFcnJwl, tp) {}
  ~FluxFcnJwlInflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnJwlInternalInflowEuler3D : public FluxFcnInternalInflowEuler3D {

public:

  FluxFcnJwlInternalInflowEuler3D(IoData &ioData, VarFcnJwl *varFcnJwl, Type tp = CONSERVATIVE) : 
    FluxFcnInternalInflowEuler3D(varFcnJwl, tp) {}
  ~FluxFcnJwlInternalInflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnJwlOutflowEuler3D : public FluxFcnOutflowEuler3D {

public:

  FluxFcnJwlOutflowEuler3D(IoData &ioData, VarFcnJwl *varFcnJwl, Type tp = CONSERVATIVE) :
    FluxFcnOutflowEuler3D(varFcnJwl, tp) {}
  ~FluxFcnJwlOutflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnJwlInternalOutflowEuler3D : public FluxFcnInternalOutflowEuler3D {

public:

  FluxFcnJwlInternalOutflowEuler3D(IoData &ioData, VarFcnJwl *varFcnJwl, Type tp = CONSERVATIVE) : 
    FluxFcnInternalOutflowEuler3D(varFcnJwl, tp) {}
  ~FluxFcnJwlInternalOutflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

#endif
