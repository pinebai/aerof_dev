#ifndef _FLUX_FCN_DESC_H_
#define _FLUX_FCN_DESC_H_

#include <FluxFcnBase.h>

#include <LowMachPrec.h>

class VarFcnBase;
class IoData;

//------------------------------------------------------------------------------

class FluxFcnFDJacRoeEuler3D : public FluxFcnFD<5> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

 public:

  FluxFcnFDJacRoeEuler3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp):
    FluxFcnFD<5> (vf,tp) { sprec.setup(ioData), gamma = gg; } 

  ~FluxFcnFDJacRoeEuler3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacRoeEuler3D : public FluxFcnBase {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnApprJacRoeEuler3D(IoData &ioData, int rs, double gg, VarFcnBase *vf, Type tp) :
    FluxFcnBase(vf,tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }

  ~FluxFcnApprJacRoeEuler3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnExactJacRoeEuler3D : public FluxFcnBase {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnExactJacRoeEuler3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp) : 
    FluxFcnBase(vf, tp) { sprec.setup(ioData); gamma = gg; }

  ~FluxFcnExactJacRoeEuler3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnFDJacHLLEEuler3D : public FluxFcnFD<5> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

 public:

  FluxFcnFDJacHLLEEuler3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp) :
    FluxFcnFD<5> (vf,tp) { sprec.setup(ioData), gamma = gg; }

  ~FluxFcnFDJacHLLEEuler3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacHLLEEuler3D : public FluxFcnBase {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnApprJacHLLEEuler3D(IoData &ioData, int rs, double gg, VarFcnBase *vf, Type tp) :
    FluxFcnBase(vf,tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }

  ~FluxFcnApprJacHLLEEuler3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnFDJacHLLCEuler3D : public FluxFcnFD<5> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

 public:

  FluxFcnFDJacHLLCEuler3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp) :
    FluxFcnFD<5> (vf,tp) { sprec.setup(ioData), gamma = gg; }

  ~FluxFcnFDJacHLLCEuler3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacHLLCEuler3D : public FluxFcnBase {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnApprJacHLLCEuler3D(IoData &ioData, int rs, double gg, VarFcnBase *vf, Type tp) :
    FluxFcnBase(vf,tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }

  ~FluxFcnApprJacHLLCEuler3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnVanLeerEuler3D : public FluxFcnBase {

public:

  FluxFcnVanLeerEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnVanLeerEuler3D() {}
  
protected:
  void evalFlux(double, double, double *, double, double *, double *, int);
  void evalJac(double, double, double *, double, double *, double *, int);  

// Included (MB)
  void evalDerivativeOfFlux(double, double, double, double *, double *, double, double, double *, double *, double, double *, int);

};

//------------------------------------------------------------------------------

class FluxFcnWallEuler3D : public FluxFcnFD<5> {

public:

  FluxFcnWallEuler3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnWallEuler3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnPorousWallEuler3D : public FluxFcnFD<5> {

public:

  FluxFcnPorousWallEuler3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnPorousWallEuler3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnGhidagliaEuler3D : public FluxFcnFD<5> {

public:

  FluxFcnGhidagliaEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnGhidagliaEuler3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnModifiedGhidagliaEuler3D : public FluxFcnFD<5> {

public:

  FluxFcnModifiedGhidagliaEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnModifiedGhidagliaEuler3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnInflowEuler3D : public FluxFcnFD<5> {

public:

  FluxFcnInflowEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnInflowEuler3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnInternalInflowEuler3D : public FluxFcnBase {

public:

  FluxFcnInternalInflowEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnInternalInflowEuler3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnDirectStateInflowEuler3D : public FluxFcnFD<5> {

public:

  FluxFcnDirectStateInflowEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnDirectStateInflowEuler3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnMassFlowInflowEuler3D : public FluxFcnFD<5> {

public:

  FluxFcnMassFlowInflowEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnMassFlowInflowEuler3D() {}
  
};

//------------------------------------------------------------------------------
class FluxFcnOutflowEuler3D : public FluxFcnFD<5> {
  
public:

  FluxFcnOutflowEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnOutflowEuler3D() {} 
  
};

//------------------------------------------------------------------------------

class FluxFcnInternalOutflowEuler3D : public FluxFcnBase {

public:

  FluxFcnInternalOutflowEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnInternalOutflowEuler3D() {}
  
};

//------------------------------------------------------------------------------
//
class FluxFcnDirectStateOutflowEuler3D : public FluxFcnFD<5> {

public:

  FluxFcnDirectStateOutflowEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnDirectStateOutflowEuler3D() {}
  
};

//------------------------------------------------------------------------------
//
class FluxFcnMassFlowOutflowEuler3D : public FluxFcnFD<5> {

public:

  FluxFcnMassFlowOutflowEuler3D(VarFcnBase *vf, Type tp) :
    FluxFcnFD<5>(vf, tp) {}
  ~FluxFcnMassFlowOutflowEuler3D() {}
  
};

//------------------------------------------------------------------------------
//turbulence

class FluxFcnFDJacRoeSA3D : public FluxFcnFD<6> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;
  
 public:
  
  FluxFcnFDJacRoeSA3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) { sprec.setup(ioData), gamma = gg; }
  ~FluxFcnFDJacRoeSA3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnApprJacRoeSA3D : public FluxFcnBase {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;
  
public:

  FluxFcnApprJacRoeSA3D(IoData &ioData, int rs, double gg, VarFcnBase* vf, Type tp = CONSERVATIVE) : 
    FluxFcnBase(vf, tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }
  ~FluxFcnApprJacRoeSA3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnExactJacRoeSA3D : public FluxFcnBase {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnExactJacRoeSA3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) { sprec.setup(ioData); gamma = gg; }
  ~FluxFcnExactJacRoeSA3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnFDJacHLLESA3D : public FluxFcnFD<6> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

 public:

  FluxFcnFDJacHLLESA3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) { sprec.setup(ioData), gamma = gg; }
  ~FluxFcnFDJacHLLESA3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacHLLESA3D : public FluxFcnBase {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnApprJacHLLESA3D(IoData &ioData, int rs, double gg, VarFcnBase* vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }
  ~FluxFcnApprJacHLLESA3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnFDJacHLLCSA3D : public FluxFcnFD<6> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

 public:

  FluxFcnFDJacHLLCSA3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) { sprec.setup(ioData), gamma = gg; }
  ~FluxFcnFDJacHLLCSA3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacHLLCSA3D : public FluxFcnBase {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnApprJacHLLCSA3D(IoData &ioData, int rs, double gg, VarFcnBase* vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }
  ~FluxFcnApprJacHLLCSA3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnWallSA3D : public FluxFcnFD<6> {

public:

  FluxFcnWallSA3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) {}
  ~FluxFcnWallSA3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnPorousWallSA3D : public FluxFcnFD<6> {

public:

  FluxFcnPorousWallSA3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) {}
  ~FluxFcnPorousWallSA3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnGhidagliaSA3D : public FluxFcnFD<6>{

public:

  FluxFcnGhidagliaSA3D(VarFcnBase *vf, Type tp) :
    FluxFcnFD<6>(vf,tp) {}

  ~FluxFcnGhidagliaSA3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnOutflowSA3D : public FluxFcnFD<6> {

public:

  FluxFcnOutflowSA3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) {}
  ~FluxFcnOutflowSA3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnInternalInflowSA3D : public FluxFcnBase {

public:

  FluxFcnInternalInflowSA3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnInternalInflowSA3D() {}
  
};

//------------------------------------------------------------------------------
//
class FluxFcnDirectStateInflowSA3D : public FluxFcnFD<6> {

public:

  FluxFcnDirectStateInflowSA3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) {}
  ~FluxFcnDirectStateInflowSA3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnMassFlowInflowSA3D : public FluxFcnFD<6> {

public:

  FluxFcnMassFlowInflowSA3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) {}
  ~FluxFcnMassFlowInflowSA3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnInternalOutflowSA3D : public FluxFcnBase {

public:

  FluxFcnInternalOutflowSA3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnInternalOutflowSA3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnDirectStateOutflowSA3D : public FluxFcnFD<6> {

public:

  FluxFcnDirectStateOutflowSA3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) {}
  ~FluxFcnDirectStateOutflowSA3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnMassFlowOutflowSA3D : public FluxFcnFD<6> {

public:

  FluxFcnMassFlowOutflowSA3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<6>(vf, tp) {}
  ~FluxFcnMassFlowOutflowSA3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnRoeSAturb3D : public FluxFcnBase {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnRoeSAturb3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) { sprec.setup(ioData), gamma = gg; }
  ~FluxFcnRoeSAturb3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnGhidagliaSAturb3D : public FluxFcnBase{

public:
  FluxFcnGhidagliaSAturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf,tp) {}

  ~FluxFcnGhidagliaSAturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnWallSAturb3D : public FluxFcnBase {

public:
  FluxFcnWallSAturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnWallSAturb3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnPorousWallSAturb3D : public FluxFcnBase {

public:
  FluxFcnPorousWallSAturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnPorousWallSAturb3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnOutflowSAturb3D : public FluxFcnBase {

public:
  FluxFcnOutflowSAturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnOutflowSAturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnInternalInflowSAturb3D : public FluxFcnBase {

public:
  FluxFcnInternalInflowSAturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnInternalInflowSAturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnDirectStateInflowSAturb3D : public FluxFcnBase {

public:
  FluxFcnDirectStateInflowSAturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnDirectStateInflowSAturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnMassFlowInflowSAturb3D : public FluxFcnBase {

public:
  FluxFcnMassFlowInflowSAturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnMassFlowInflowSAturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnInternalOutflowSAturb3D : public FluxFcnBase {

public:
  FluxFcnInternalOutflowSAturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnInternalOutflowSAturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnDirectStateOutflowSAturb3D : public FluxFcnBase {

public:
  FluxFcnDirectStateOutflowSAturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnDirectStateOutflowSAturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnMassFlowOutflowSAturb3D : public FluxFcnBase {

public:
  FluxFcnMassFlowOutflowSAturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnMassFlowOutflowSAturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnFDJacRoeKE3D : public FluxFcnFD<7> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnFDJacRoeKE3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<7>(vf, tp) { sprec.setup(ioData), gamma = gg; }
  ~FluxFcnFDJacRoeKE3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacRoeKE3D : public FluxFcnBase {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:
  FluxFcnApprJacRoeKE3D(IoData &ioData, int rs, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }
  ~FluxFcnApprJacRoeKE3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnExactJacRoeKE3D : public FluxFcnBase {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

public:
  FluxFcnExactJacRoeKE3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) { sprec.setup(ioData); gamma = gg; }
  ~FluxFcnExactJacRoeKE3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnFDJacHLLEKE3D : public FluxFcnFD<7> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnFDJacHLLEKE3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<7>(vf, tp) { sprec.setup(ioData), gamma = gg; }
  ~FluxFcnFDJacHLLEKE3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacHLLEKE3D : public FluxFcnBase {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:
  FluxFcnApprJacHLLEKE3D(IoData &ioData, int rs, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }
  ~FluxFcnApprJacHLLEKE3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnFDJacHLLCKE3D : public FluxFcnFD<7> {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

public:

  FluxFcnFDJacHLLCKE3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<7>(vf, tp) { sprec.setup(ioData), gamma = gg; }
  ~FluxFcnFDJacHLLCKE3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnApprJacHLLCKE3D : public FluxFcnBase {

 protected:
  int rshift;
  double gamma;
  SpatialLowMachPrec sprec;

public:
  FluxFcnApprJacHLLCKE3D(IoData &ioData, int rs, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) { sprec.setup(ioData), rshift = rs; gamma = gg; }
  ~FluxFcnApprJacHLLCKE3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnWallKE3D : public FluxFcnFD<7> {

public:

  FluxFcnWallKE3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<7>(vf, tp) {}
  ~FluxFcnWallKE3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnPorousWallKE3D : public FluxFcnFD<7> {

public:

  FluxFcnPorousWallKE3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<7>(vf, tp) {}
  ~FluxFcnPorousWallKE3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnGhidagliaKE3D : public FluxFcnFD<7>{

public:

  FluxFcnGhidagliaKE3D(VarFcnBase *vf, Type tp) :
    FluxFcnFD<7>(vf,tp) {}

  ~FluxFcnGhidagliaKE3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnOutflowKE3D : public FluxFcnFD<7> {

public:

  FluxFcnOutflowKE3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<7>(vf, tp) {}
  ~FluxFcnOutflowKE3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnDirectStateInflowKE3D : public FluxFcnFD<7> {

public:

  FluxFcnDirectStateInflowKE3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<7>(vf, tp) {}
  ~FluxFcnDirectStateInflowKE3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnMassFlowInflowKE3D : public FluxFcnFD<7> {

public:

  FluxFcnMassFlowInflowKE3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<7>(vf, tp) {}
  ~FluxFcnMassFlowInflowKE3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnDirectStateOutflowKE3D : public FluxFcnFD<7> {

public:

  FluxFcnDirectStateOutflowKE3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<7>(vf, tp) {}
  ~FluxFcnDirectStateOutflowKE3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnMassFlowOutflowKE3D : public FluxFcnFD<7> {

public:

  FluxFcnMassFlowOutflowKE3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnFD<7>(vf, tp) {}
  ~FluxFcnMassFlowOutflowKE3D() {}
  
};

//------------------------------------------------------------------------------

class FluxFcnRoeKEturb3D : public FluxFcnBase {

 protected:
  double gamma;
  SpatialLowMachPrec sprec;

public:
  FluxFcnRoeKEturb3D(IoData &ioData, double gg, VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) { sprec.setup(ioData), gamma = gg; }

  ~FluxFcnRoeKEturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnWallKEturb3D : public FluxFcnBase {

public:

  FluxFcnWallKEturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnWallKEturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnPorousWallKEturb3D : public FluxFcnBase {

public:

  FluxFcnPorousWallKEturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnPorousWallKEturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnGhidagliaKEturb3D : public FluxFcnBase{

public:
  FluxFcnGhidagliaKEturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf,tp) {}

  ~FluxFcnGhidagliaKEturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnOutflowKEturb3D : public FluxFcnBase {

public:

  FluxFcnOutflowKEturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnOutflowKEturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnDirectStateInflowKEturb3D : public FluxFcnBase {

public:
  FluxFcnDirectStateInflowKEturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnDirectStateInflowKEturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnMassFlowInflowKEturb3D : public FluxFcnBase {

public:
  FluxFcnMassFlowInflowKEturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnMassFlowInflowKEturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnDirectStateOutflowKEturb3D : public FluxFcnBase {

public:
  FluxFcnDirectStateOutflowKEturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnDirectStateOutflowKEturb3D() {}

};

//------------------------------------------------------------------------------

class FluxFcnMassFlowOutflowKEturb3D : public FluxFcnBase {

public:
  FluxFcnMassFlowOutflowKEturb3D(VarFcnBase *vf, Type tp = CONSERVATIVE) :
    FluxFcnBase(vf, tp) {}
  ~FluxFcnMassFlowOutflowKEturb3D() {}

};

//------------------------------------------------------------------------------

#endif
