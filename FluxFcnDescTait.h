#ifndef _FLUX_FCN_DESC_TAIT_H_
#define _FLUX_FCN_DESC_TAIT_H_

#include <FluxFcnDesc.h>

class IoData;
#include "VarFcnBase.h"
#include "VarFcnTait.h"
#include "VarFcnTaitSA.h"
#include "VarFcnTaitKE.h"

//------------------------------------------------------------------------------

class FluxFcnTaitFDJacRoeEuler3D : public FluxFcnFDJacRoeEuler3D {

public:

  FluxFcnTaitFDJacRoeEuler3D(double gg, IoData &ioData, VarFcnTait *varFcnTait, Type tp = CONSERVATIVE) :
    FluxFcnFDJacRoeEuler3D(ioData, gg, varFcnTait, tp) {}
  ~FluxFcnTaitFDJacRoeEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnTaitApprJacRoeEuler3D : public FluxFcnApprJacRoeEuler3D {

public:

  FluxFcnTaitApprJacRoeEuler3D(int rs, double gg, IoData &ioData, VarFcnTait *varFcnTait, Type tp = CONSERVATIVE) : 
    FluxFcnApprJacRoeEuler3D(ioData, rs, gg, varFcnTait, tp) {}
  ~FluxFcnTaitApprJacRoeEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);


};

//------------------------------------------------------------------------------

class FluxFcnTaitExactJacRoeEuler3D : public FluxFcnExactJacRoeEuler3D {

public:

  FluxFcnTaitExactJacRoeEuler3D(double gg, IoData &ioData, VarFcnTait *varFcnTait, Type tp = CONSERVATIVE) : 
    FluxFcnExactJacRoeEuler3D(ioData, gg, varFcnTait, tp) {}
  ~FluxFcnTaitExactJacRoeEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);


};

//------------------------------------------------------------------------------

class FluxFcnTaitWallEuler3D : public FluxFcnWallEuler3D {

public:
  FluxFcnTaitWallEuler3D(IoData &ioData, VarFcnTait *varFcnTait, Type tp=CONSERVATIVE) :
    FluxFcnWallEuler3D(varFcnTait, tp) {}
  ~FluxFcnTaitWallEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnTaitGhidagliaEuler3D : public FluxFcnGhidagliaEuler3D {

public:

  FluxFcnTaitGhidagliaEuler3D(IoData &ioData, VarFcnTait *varFcnTait, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaEuler3D(varFcnTait, tp) {}
  ~FluxFcnTaitGhidagliaEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  //void computeJacobian(double length, double irey, double *normal, double normalVel, 
  //		       double *VL, double *Ub, double *jacL, bool useLimiter);
};

class FluxFcnTaitModifiedGhidagliaEuler3D : public FluxFcnModifiedGhidagliaEuler3D {

public:

  FluxFcnTaitModifiedGhidagliaEuler3D(IoData &ioData, VarFcnTait *varFcnTait, Type tp = CONSERVATIVE) :
    FluxFcnModifiedGhidagliaEuler3D(varFcnTait, tp) {}
  ~FluxFcnTaitModifiedGhidagliaEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

  void computeJacobianFarfield(double, double, double *, double, double *, double *, double *,
			       bool);
};


//------------------------------------------------------------------------------

class FluxFcnTaitInflowEuler3D : public FluxFcnInflowEuler3D {

public:

  FluxFcnTaitInflowEuler3D(IoData &ioData, VarFcnTait *varFcnTait, Type tp = CONSERVATIVE) :
    FluxFcnInflowEuler3D(varFcnTait, tp) {}
  ~FluxFcnTaitInflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnTaitOutflowEuler3D : public FluxFcnOutflowEuler3D {

public:

  FluxFcnTaitOutflowEuler3D(IoData &ioData, VarFcnTait *varFcnTait, Type tp = CONSERVATIVE) :
    FluxFcnOutflowEuler3D(varFcnTait, tp) {}
  ~FluxFcnTaitOutflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnTaitInternalInflowEuler3D : public FluxFcnInternalInflowEuler3D {

public:

  FluxFcnTaitInternalInflowEuler3D(IoData &ioData, VarFcnTait *varFcnTait, Type tp = CONSERVATIVE) : 
    FluxFcnInternalInflowEuler3D(varFcnTait, tp) {}
  ~FluxFcnTaitInternalInflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnTaitInternalOutflowEuler3D : public FluxFcnInternalOutflowEuler3D {

public:

  FluxFcnTaitInternalOutflowEuler3D(IoData &ioData, VarFcnTait *varFcnTait, Type tp = CONSERVATIVE) : 
    FluxFcnInternalOutflowEuler3D(varFcnTait, tp) {}
  ~FluxFcnTaitInternalOutflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

class FluxFcnTaitApprJacRoeSA3D : public FluxFcnApprJacRoeSA3D {

public:

  FluxFcnTaitApprJacRoeSA3D(int rs, double gg, IoData &ioData, VarFcnTaitSA *varFcnTait, Type tp = CONSERVATIVE) :
    FluxFcnApprJacRoeSA3D(ioData, rs, gg, varFcnTait, tp) {}
  ~FluxFcnTaitApprJacRoeSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);


};

//------------------------------------------------------------------------------

class FluxFcnTaitWallSA3D : public FluxFcnWallSA3D {

public:
  FluxFcnTaitWallSA3D(IoData &ioData, VarFcnTaitSA *varFcnTait, Type tp=CONSERVATIVE) :
    FluxFcnWallSA3D(varFcnTait, tp) {}
  ~FluxFcnTaitWallSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnTaitGhidagliaSA3D : public FluxFcnGhidagliaSA3D {

public:

  FluxFcnTaitGhidagliaSA3D(IoData &ioData, VarFcnTaitSA *varFcnTait, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaSA3D(varFcnTait, tp) {}
  ~FluxFcnTaitGhidagliaSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnTaitRoeSAturb3D : public FluxFcnRoeSAturb3D {

public:

  FluxFcnTaitRoeSAturb3D(double gg, IoData &ioData, VarFcnTaitSA *varFcnTait, Type tp = CONSERVATIVE) :
    FluxFcnRoeSAturb3D(ioData, gg, varFcnTait, tp) {}
  ~FluxFcnTaitRoeSAturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnTaitWallSAturb3D : public FluxFcnWallSAturb3D {

public:
  FluxFcnTaitWallSAturb3D(IoData &ioData, VarFcnTaitSA *varFcnTait, Type tp = CONSERVATIVE) :
    FluxFcnWallSAturb3D(varFcnTait, tp) {}
  ~FluxFcnTaitWallSAturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnTaitGhidagliaSAturb3D : public FluxFcnGhidagliaSAturb3D {

public:
  FluxFcnTaitGhidagliaSAturb3D(IoData &ioData, VarFcnTaitSA *varFcnTait, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaSAturb3D(varFcnTait, tp) {}
  ~FluxFcnTaitGhidagliaSAturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

class FluxFcnTaitApprJacRoeKE3D : public FluxFcnApprJacRoeKE3D {

public:

  FluxFcnTaitApprJacRoeKE3D(int rs, double gg, IoData &ioData, VarFcnTaitKE *varFcnTait, Type tp = CONSERVATIVE) :
    FluxFcnApprJacRoeKE3D(ioData, rs, gg, varFcnTait, tp) {}
  ~FluxFcnTaitApprJacRoeKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);


};

//------------------------------------------------------------------------------

class FluxFcnTaitWallKE3D : public FluxFcnWallKE3D {

public:
  FluxFcnTaitWallKE3D(IoData &ioData, VarFcnTaitKE *varFcnTait, Type tp=CONSERVATIVE) :
    FluxFcnWallKE3D(varFcnTait, tp) {}
  ~FluxFcnTaitWallKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnTaitGhidagliaKE3D : public FluxFcnGhidagliaKE3D {

public:

  FluxFcnTaitGhidagliaKE3D(IoData &ioData, VarFcnTaitKE *varFcnTait, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaKE3D(varFcnTait, tp) {}
  ~FluxFcnTaitGhidagliaKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnTaitRoeKEturb3D : public FluxFcnRoeKEturb3D {

public:

  FluxFcnTaitRoeKEturb3D(double gg, IoData &ioData, VarFcnTaitKE *varFcnTait, Type tp = CONSERVATIVE) :
    FluxFcnRoeKEturb3D(ioData, gg, varFcnTait, tp) {}
  ~FluxFcnTaitRoeKEturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnTaitWallKEturb3D : public FluxFcnWallKEturb3D {

public:
  FluxFcnTaitWallKEturb3D(IoData &ioData, VarFcnTaitKE *varFcnTait, Type tp = CONSERVATIVE) :
    FluxFcnWallKEturb3D(varFcnTait, tp) {}
  ~FluxFcnTaitWallKEturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnTaitGhidagliaKEturb3D : public FluxFcnGhidagliaKEturb3D {

public:
  FluxFcnTaitGhidagliaKEturb3D(IoData &ioData, VarFcnTaitKE *varFcnTait, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaKEturb3D(varFcnTait, tp) {}
  ~FluxFcnTaitGhidagliaKEturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};


//------------------------------------------------------------------------------


#endif
