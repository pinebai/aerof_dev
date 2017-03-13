#ifndef _FLUX_FCN_DESC_SG_H_
#define _FLUX_FCN_DESC_SG_H_

#include <FluxFcnDesc.h>

class IoData;
#include "VarFcnBase.h"
#include "VarFcnSGEuler.h"
#include "VarFcnSGSA.h"
#include "VarFcnSGKE.h"

//------------------------------------------------------------------------------

class FluxFcnSGFDJacRoeEuler3D : public FluxFcnFDJacRoeEuler3D {

public:

  FluxFcnSGFDJacRoeEuler3D(double gg, IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) :
    FluxFcnFDJacRoeEuler3D(ioData, gg, varFcnSGEuler, tp) {}
  ~FluxFcnSGFDJacRoeEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  );

//  virtual void compute_dFluxdNormal_dFluxdNormalVel_dFluxdVL_dFluxdVR(double *n, double nv, double *vl, double *vr,
//                                                              double dmach, double *f, double *dfdn,
//                                                              double *fdnv, double *dfdvl, double *dfdvr);
 
  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  )
  {
    std::cout << "\n !!! FluxFcnSGFDJacRoeEuler3D::computeDerivative (11 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  virtual void computeDerivativeOperators(double *n, double nv, double *v,
                                          double *ub, double dfdn[7][3], double dfdv[7][1], double dfdub[7][7])
  {
    std::cout << "\n !!! FluxFcnSGFDJacRoeEuler3D::computeDerivativeOperators is not implemented !!!\n\n";
    exit(1);
  }

};

//------------------------------------------------------------------------------

class FluxFcnSGApprJacRoeEuler3D : public FluxFcnApprJacRoeEuler3D {

public:

  FluxFcnSGApprJacRoeEuler3D(int rs, double gg, IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) : 
    FluxFcnApprJacRoeEuler3D(ioData, rs, gg, varFcnSGEuler, tp) {}
  ~FluxFcnSGApprJacRoeEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  );

//  virtual void compute_dFluxdNormal_dFluxdNormalVel_dFluxdVL_dFluxdVR(double *n, double nv, double *vl, double *vr,
//                                                              double dmach, double *f, double *dfdn,
//                                                              double *fdnv, double *dfdvl, double *dfdvr);
 
  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  )
  {
    std::cout << "\n !!! FluxFcnSGApprJacRoeEuler3D::computeDerivative (11 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  virtual void computeDerivativeOperators(double *n, double nv, double *v,
                                          double *ub, double dfdn[7][3], double dfdv[7][1], double dfdub[7][7])
  {
    std::cout << "\n !!! FluxFcnSGApprJacRoeEuler3D::computeDerivativeOperators is not implemented !!!\n\n";
    exit(1);
  }

};

//------------------------------------------------------------------------------

class FluxFcnSGExactJacRoeEuler3D : public FluxFcnExactJacRoeEuler3D {

public:

  FluxFcnSGExactJacRoeEuler3D(double gg, IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) : 
    FluxFcnExactJacRoeEuler3D(ioData, gg, varFcnSGEuler, tp) {}
  ~FluxFcnSGExactJacRoeEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  );

  virtual void compute_dFluxdNormal_dFluxdNormalVel_dFluxdVL_dFluxdVR(double *n, double nv, double *vl, double *vr,
                                                              double dmach, double *f, double (*)[3],
                                                              double *, double (*)[7], double (*)[7]); 

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  )
  {
    std::cout << "\n !!! FluxFcnSGExactJacRoeEuler3D::computeDerivative (11 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  virtual void computeDerivativeOperators(double *n, double nv, double *v,
                                          double *ub, double dfdn[7][3], double dfdv[7][1], double dfdub[7][7])
  {
    std::cout << "\n !!! FluxFcnSGExactJacRoeEuler3D::computeDerivativeOperators is not implemented !!!\n\n";
    exit(1);
  }


};

//------------------------------------------------------------------------------

class FluxFcnSGFDJacHLLEEuler3D : public FluxFcnFDJacHLLEEuler3D {

public:

  FluxFcnSGFDJacHLLEEuler3D(double gg, IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) :
    FluxFcnFDJacHLLEEuler3D(ioData, gg, varFcnSGEuler, tp) {}
  ~FluxFcnSGFDJacHLLEEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGApprJacHLLEEuler3D : public FluxFcnApprJacHLLEEuler3D {

public:

  FluxFcnSGApprJacHLLEEuler3D(int rs, double gg, IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) :
    FluxFcnApprJacHLLEEuler3D(ioData, rs, gg, varFcnSGEuler, tp) {}
  ~FluxFcnSGApprJacHLLEEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGFDJacHLLCEuler3D : public FluxFcnFDJacHLLCEuler3D {

public:

  FluxFcnSGFDJacHLLCEuler3D(double gg, IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) :
    FluxFcnFDJacHLLCEuler3D(ioData, gg, varFcnSGEuler, tp) {}
  ~FluxFcnSGFDJacHLLCEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGApprJacHLLCEuler3D : public FluxFcnApprJacHLLCEuler3D {

public:

  FluxFcnSGApprJacHLLCEuler3D(int rs, double gg, IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) :
    FluxFcnApprJacHLLCEuler3D(ioData, rs, gg, varFcnSGEuler, tp) {}
  ~FluxFcnSGApprJacHLLCEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGVanLeerEuler3D : public FluxFcnVanLeerEuler3D {

public:

  FluxFcnSGVanLeerEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) : 
    FluxFcnVanLeerEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGVanLeerEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  );

 
  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  )
  {
    std::cout << "\n !!! FluxFcnSGVanLeerEuler3D::computeDerivative (11 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  virtual void computeDerivativeOperators(double *n, double nv, double *v,
                                          double *ub, double dfdn[7][3], double dfdv[7][1], double dfdub[7][7])
  {
    std::cout << "\n !!! FluxFcnSGVanLeerEuler3D::computeDerivativeOperators is not implemented !!!\n\n";
    exit(1);
  }

};

//------------------------------------------------------------------------------

class FluxFcnSGWallEuler3D : public FluxFcnWallEuler3D {

public:
  FluxFcnSGWallEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp=CONSERVATIVE) :
    FluxFcnWallEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGWallEuler3D() { vf = 0; }
  
  void compute(double, double, double *, double, double *, double *, double *, bool);

  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGWallEuler3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

  virtual void computeDerivativeOperators
  (
    double *normal, double normalVel, double *V,
    double *Ub, double dFluxdNormal[7][3], double dFluxdNormalVel[7][1], double dFluxdUb[7][7]
  );


};
//------------------------------------------------------------------------------

class FluxFcnSGPorousWallEuler3D : public FluxFcnPorousWallEuler3D {

public:
  FluxFcnSGPorousWallEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp=CONSERVATIVE) :
    FluxFcnPorousWallEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGPorousWallEuler3D() { vf = 0; }
  
  void compute(double, double, double *, double, double *, double *, double *, bool);

//  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGPorousWallEuler3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

};
//------------------------------------------------------------------------------

class FluxFcnSGGhidagliaEuler3D : public FluxFcnGhidagliaEuler3D {

public:

  FluxFcnSGGhidagliaEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGGhidagliaEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGModifiedGhidagliaEuler3D : public FluxFcnModifiedGhidagliaEuler3D {

public:

  FluxFcnSGModifiedGhidagliaEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) :
    FluxFcnModifiedGhidagliaEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGModifiedGhidagliaEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

  void computeJacobianFarfield(double, double, double *, double, double *, double *, double *,
			       bool);
};

//------------------------------------------------------------------------------

class FluxFcnSGInflowEuler3D : public FluxFcnInflowEuler3D {

public:

  FluxFcnSGInflowEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) :
    FluxFcnInflowEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGInflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGInflowEuler3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

  virtual void computeDerivativeOperators
  (
    double*, double, double*, double*,
    double [7][3], double[7][1], double [7][7]
  );

};

//------------------------------------------------------------------------------

class FluxFcnSGInternalInflowEuler3D : public FluxFcnInternalInflowEuler3D {

public:

  FluxFcnSGInternalInflowEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) : 
    FluxFcnInternalInflowEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGInternalInflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGInternalInflowEuler3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

  virtual void computeDerivativeOperators
  (
    double*, double, double*, double*,
    double [7][3], double[7][1], double [7][7]
  )
  {
    std::cout << "\n !!! FluxFcnSGInternalInflowEuler3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

};

//------------------------------------------------------------------------------

class FluxFcnSGDirectStateInflowEuler3D : public FluxFcnDirectStateInflowEuler3D {

public:

  FluxFcnSGDirectStateInflowEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) : 
    FluxFcnDirectStateInflowEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGDirectStateInflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
//  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGDirectStateInflowEuler3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

  virtual void computeDerivativeOperators
  (
    double*, double, double*, double*,
    double [7][3], double[7][1], double [7][7]
  )
  {
    std::cout << "\n !!! FluxFcnSGDirectStateInflowEuler3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

};

//------------------------------------------------------------------------------

class FluxFcnSGMassFlowInflowEuler3D : public FluxFcnMassFlowInflowEuler3D {

public:

  FluxFcnSGMassFlowInflowEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) : 
    FluxFcnMassFlowInflowEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGMassFlowInflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
//  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGMassFlowInflowEuler3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

  virtual void computeDerivativeOperators
  (
    double*, double, double*, double*,
    double [7][3], double[7][1], double [7][7]
  )
  {
    std::cout << "\n !!! FluxFcnSGMassFlowInflowEuler3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

};

//------------------------------------------------------------------------------

class FluxFcnSGOutflowEuler3D : public FluxFcnOutflowEuler3D {

public:

  FluxFcnSGOutflowEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) :
    FluxFcnOutflowEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGOutflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGOutflowEuler3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

  virtual void computeDerivativeOperators
  (
    double*, double, double*, double*,  
    double [7][3], double[7][1], double [7][7]
  );

};

//------------------------------------------------------------------------------

class FluxFcnSGInternalOutflowEuler3D : public FluxFcnInternalOutflowEuler3D {

public:

  FluxFcnSGInternalOutflowEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) : 
    FluxFcnInternalOutflowEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGInternalOutflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df,
    bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGInternalOutflowEuler3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

  virtual void computeDerivativeOperators
  (
    double*, double, double*, double*,
    double [7][3], double[7][1], double [7][7]
  )
  {
    std::cout << "\n !!! FluxFcnSGInternalOutflowEuler3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

};

//------------------------------------------------------------------------------
//
class FluxFcnSGDirectStateOutflowEuler3D : public FluxFcnDirectStateOutflowEuler3D {

public:

  FluxFcnSGDirectStateOutflowEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) : 
    FluxFcnDirectStateOutflowEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGDirectStateOutflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
//  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df,
    bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGDirectStateOutflowEuler3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

};

//------------------------------------------------------------------------------

class FluxFcnSGMassFlowOutflowEuler3D : public FluxFcnMassFlowOutflowEuler3D {

public:

  FluxFcnSGMassFlowOutflowEuler3D(IoData &ioData, VarFcnSGEuler *varFcnSGEuler, Type tp = CONSERVATIVE) : 
    FluxFcnMassFlowOutflowEuler3D(varFcnSGEuler, tp) {}
  ~FluxFcnSGMassFlowOutflowEuler3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
//  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df,
    bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGMassFlowOutflowEuler3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//turbulence
//NOTE : for the coupled solver, the varFcn correspond to the FluxFcn in the sense 
//       that if varFcn = VarFcnSGSA3D, then fluxFcn = FluxFcnSGApprJacRoeSA3D
//       HOWEVER, for the segregated solver, they do not necessarily correspond since we consider
//       two different fluxes ff1 and ff2. Then there must be two varFcn (vf1 and vf2) corresponding
//       to each fluxFcn, and not one varFcn (corresponding to the physical case we are considering)
//       for both ff1 and ff2
//
//
//NOTE2:  some FluxFcn do not really need varFcn, but they were given one varFcn still (easier to implement)
//        for instance all fluxFcn of type SAturb3D and KEturb3D do not need a real varFcn. Moreover,
//        there is no corresponding varFcn for these fluxes because the varFcn always consider at least
//        the Euler variables, but never just the turbulent variables.




class FluxFcnSGFDJacRoeSA3D : public FluxFcnFDJacRoeSA3D {

 public:
  
  FluxFcnSGFDJacRoeSA3D(double gg, IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnFDJacRoeSA3D(ioData, gg, varFcnSGSA) {}
  ~FluxFcnSGFDJacRoeSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  );

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  )
  {
    std::cout << "\n !!! FluxFcnSGFDJacRoeSA3D::computeDerivative (11 arg.) is not implemented !!!\n\n";
    exit(1);
  }

};

//------------------------------------------------------------------------------

class FluxFcnSGApprJacRoeSA3D : public FluxFcnApprJacRoeSA3D {

public:

  FluxFcnSGApprJacRoeSA3D(int rs, double gg, IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) : 
    FluxFcnApprJacRoeSA3D(ioData, rs,gg, varFcnSGSA, tp) {}
  ~FluxFcnSGApprJacRoeSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  );

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  )
  {
    std::cout << "\n !!! FluxFcnSGApprJacRoeSA3D::computeDerivative (11 arg.) is not implemented !!!\n\n";
    exit(1);
  }

};

//------------------------------------------------------------------------------

class FluxFcnSGExactJacRoeSA3D : public FluxFcnExactJacRoeSA3D {

public:

  FluxFcnSGExactJacRoeSA3D(double gg, IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnExactJacRoeSA3D(ioData, gg, varFcnSGSA, tp) {}
  ~FluxFcnSGExactJacRoeSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  );

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  )
  {
    std::cout << "\n !!! FluxFcnSGExactJacRoeSA3D::computeDerivative (11 arg.) is not implemented !!!\n\n";
    exit(1);
  }

};

//------------------------------------------------------------------------------

class FluxFcnSGFDJacHLLESA3D : public FluxFcnFDJacHLLESA3D {

 public:

  FluxFcnSGFDJacHLLESA3D(double gg, IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
  FluxFcnFDJacHLLESA3D(ioData, gg, varFcnSGSA, tp) {}
  ~FluxFcnSGFDJacHLLESA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGApprJacHLLESA3D : public FluxFcnApprJacHLLESA3D {

public:

  FluxFcnSGApprJacHLLESA3D(int rs, double gg, IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
  FluxFcnApprJacHLLESA3D(ioData, rs,gg, varFcnSGSA, tp) {}
  ~FluxFcnSGApprJacHLLESA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGFDJacHLLCSA3D : public FluxFcnFDJacHLLCSA3D {

 public:

  FluxFcnSGFDJacHLLCSA3D(double gg, IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
  FluxFcnFDJacHLLCSA3D(ioData, gg, varFcnSGSA, tp) {}
  ~FluxFcnSGFDJacHLLCSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGApprJacHLLCSA3D : public FluxFcnApprJacHLLCSA3D {

public:

  FluxFcnSGApprJacHLLCSA3D(int rs, double gg, IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
  FluxFcnApprJacHLLCSA3D(ioData, rs,gg, varFcnSGSA, tp) {}
  ~FluxFcnSGApprJacHLLCSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGWallSA3D : public FluxFcnWallSA3D {

public:

  FluxFcnSGWallSA3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnWallSA3D(varFcnSGSA, tp) {}
  ~FluxFcnSGWallSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGWallSA3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );
};

//------------------------------------------------------------------------------

class FluxFcnSGPorousWallSA3D : public FluxFcnPorousWallSA3D {

public:

  FluxFcnSGPorousWallSA3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnPorousWallSA3D(varFcnSGSA, tp) {}
  ~FluxFcnSGPorousWallSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

//  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGPorousWallSA3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );
};

//------------------------------------------------------------------------------

class FluxFcnSGGhidagliaSA3D : public FluxFcnGhidagliaSA3D {

public:

  FluxFcnSGGhidagliaSA3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaSA3D(varFcnSGSA, tp) {}
  ~FluxFcnSGGhidagliaSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGOutflowSA3D : public FluxFcnOutflowSA3D {

public:

  FluxFcnSGOutflowSA3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnOutflowSA3D(varFcnSGSA, tp) {}
  ~FluxFcnSGOutflowSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGOutflowSA3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

};

//------------------------------------------------------------------------------

class FluxFcnSGInternalInflowSA3D : public FluxFcnInternalInflowSA3D {

public:

  FluxFcnSGInternalInflowSA3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnInternalInflowSA3D(varFcnSGSA, tp) {}
  ~FluxFcnSGInternalInflowSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGInternalInflowSA3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

};

//------------------------------------------------------------------------------

class FluxFcnSGDirectStateInflowSA3D : public FluxFcnDirectStateInflowSA3D {

public:

  FluxFcnSGDirectStateInflowSA3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnDirectStateInflowSA3D(varFcnSGSA, tp) {}
  ~FluxFcnSGDirectStateInflowSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
//  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGDirectStateInflowSA3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

};

//------------------------------------------------------------------------------

class FluxFcnSGMassFlowInflowSA3D : public FluxFcnMassFlowInflowSA3D {

public:

  FluxFcnSGMassFlowInflowSA3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnMassFlowInflowSA3D(varFcnSGSA, tp) {}
  ~FluxFcnSGMassFlowInflowSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
//  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGMassFlowInflowSA3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

};

//------------------------------------------------------------------------------

class FluxFcnSGInternalOutflowSA3D : public FluxFcnInternalOutflowSA3D {

public:

  FluxFcnSGInternalOutflowSA3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnInternalOutflowSA3D(varFcnSGSA, tp) {}
  ~FluxFcnSGInternalOutflowSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGInternalOutflowSA3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

};

//------------------------------------------------------------------------------

class FluxFcnSGDirectStateOutflowSA3D : public FluxFcnDirectStateOutflowSA3D {

public:

  FluxFcnSGDirectStateOutflowSA3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnDirectStateOutflowSA3D(varFcnSGSA, tp) {}
  ~FluxFcnSGDirectStateOutflowSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
//  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGDirectStateOutflowSA3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

};

//------------------------------------------------------------------------------

class FluxFcnSGMassFlowOutflowSA3D : public FluxFcnMassFlowOutflowSA3D {

public:

  FluxFcnSGMassFlowOutflowSA3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnMassFlowOutflowSA3D(varFcnSGSA, tp) {}
  ~FluxFcnSGMassFlowOutflowSA3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
//  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGMassFlowOutflowSA3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

};

//------------------------------------------------------------------------------

class FluxFcnSGRoeSAturb3D : public FluxFcnRoeSAturb3D {

public:

  FluxFcnSGRoeSAturb3D(double gg, IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnRoeSAturb3D(ioData, gg, varFcnSGSA, tp) {}
  ~FluxFcnSGRoeSAturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGWallSAturb3D : public FluxFcnWallSAturb3D {

public:
  FluxFcnSGWallSAturb3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnWallSAturb3D(varFcnSGSA, tp) {}
  ~FluxFcnSGWallSAturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGPorousWallSAturb3D : public FluxFcnPorousWallSAturb3D {

public:
  FluxFcnSGPorousWallSAturb3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnPorousWallSAturb3D(varFcnSGSA, tp) {}
  ~FluxFcnSGPorousWallSAturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGGhidagliaSAturb3D : public FluxFcnGhidagliaSAturb3D {

public:
  FluxFcnSGGhidagliaSAturb3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaSAturb3D(varFcnSGSA, tp) {}
  ~FluxFcnSGGhidagliaSAturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGOutflowSAturb3D : public FluxFcnOutflowSAturb3D {

public:
  FluxFcnSGOutflowSAturb3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnOutflowSAturb3D(varFcnSGSA, tp) {}
  ~FluxFcnSGOutflowSAturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGInternalInflowSAturb3D : public FluxFcnInternalInflowSAturb3D {

public:
  FluxFcnSGInternalInflowSAturb3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnInternalInflowSAturb3D(varFcnSGSA, tp) {}
  ~FluxFcnSGInternalInflowSAturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGDirectStateInflowSAturb3D : public FluxFcnDirectStateInflowSAturb3D {

public:
  FluxFcnSGDirectStateInflowSAturb3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnDirectStateInflowSAturb3D(varFcnSGSA, tp) {}
  ~FluxFcnSGDirectStateInflowSAturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGMassFlowInflowSAturb3D : public FluxFcnMassFlowInflowSAturb3D {

public:
  FluxFcnSGMassFlowInflowSAturb3D(IoData &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnMassFlowInflowSAturb3D(varFcnSGSA, tp) {}
  ~FluxFcnSGMassFlowInflowSAturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGInternalOutflowSAturb3D : public FluxFcnInternalOutflowSAturb3D {

public:
  FluxFcnSGInternalOutflowSAturb3D(IoData  &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnInternalOutflowSAturb3D(varFcnSGSA, tp) {}
  ~FluxFcnSGInternalOutflowSAturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGDirectStateOutflowSAturb3D : public FluxFcnDirectStateOutflowSAturb3D {

public:
  FluxFcnSGDirectStateOutflowSAturb3D(IoData  &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnDirectStateOutflowSAturb3D(varFcnSGSA, tp) {}
  ~FluxFcnSGDirectStateOutflowSAturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGMassFlowOutflowSAturb3D : public FluxFcnMassFlowOutflowSAturb3D {

public:
  FluxFcnSGMassFlowOutflowSAturb3D(IoData  &ioData, VarFcnSGSA *varFcnSGSA, Type tp = CONSERVATIVE) :
    FluxFcnMassFlowOutflowSAturb3D(varFcnSGSA, tp) {}
  ~FluxFcnSGMassFlowOutflowSAturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGFDJacRoeKE3D : public FluxFcnFDJacRoeKE3D{

public:

  FluxFcnSGFDJacRoeKE3D(double gg, IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnFDJacRoeKE3D(ioData, gg, varFcnSGKE, tp) {}
  ~FluxFcnSGFDJacRoeKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  );

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  )
  {
    std::cout << "\n !!! FluxFcnSGFDJacRoeKE3D::computeDerivative (11 arg.) is not implemented !!!\n\n";
    exit(1);
  }

};

//------------------------------------------------------------------------------

class FluxFcnSGApprJacRoeKE3D : public FluxFcnApprJacRoeKE3D {

public:
  FluxFcnSGApprJacRoeKE3D(int rs, double gg, IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnApprJacRoeKE3D(ioData, rs, gg, varFcnSGKE, tp) {}
  ~FluxFcnSGApprJacRoeKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  );

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  )
  {
    std::cout << "\n !!! FluxFcnSGApprJacRoeKE3D::computeDerivative (11 arg.) is not implemented !!!\n\n";
    exit(1);
  }

};

//------------------------------------------------------------------------------

class FluxFcnSGFDJacHLLEKE3D : public FluxFcnFDJacHLLEKE3D{

public:

  FluxFcnSGFDJacHLLEKE3D(double gg, IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnFDJacHLLEKE3D(ioData, gg, varFcnSGKE, tp) {}
  ~FluxFcnSGFDJacHLLEKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);


};

//------------------------------------------------------------------------------

class FluxFcnSGApprJacHLLEKE3D : public FluxFcnApprJacHLLEKE3D {

public:
  FluxFcnSGApprJacHLLEKE3D(int rs, double gg, IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnApprJacHLLEKE3D(ioData, rs, gg, varFcnSGKE, tp) {}
  ~FluxFcnSGApprJacHLLEKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGFDJacHLLCKE3D : public FluxFcnFDJacHLLCKE3D{

public:

  FluxFcnSGFDJacHLLCKE3D(double gg, IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnFDJacHLLCKE3D(ioData, gg, varFcnSGKE, tp) {}
  ~FluxFcnSGFDJacHLLCKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);


};

//------------------------------------------------------------------------------

class FluxFcnSGApprJacHLLCKE3D : public FluxFcnApprJacHLLCKE3D {

public:
  FluxFcnSGApprJacHLLCKE3D(int rs, double gg, IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnApprJacHLLCKE3D(ioData, rs, gg, varFcnSGKE, tp) {}
  ~FluxFcnSGApprJacHLLCKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGExactJacRoeKE3D : public FluxFcnExactJacRoeKE3D {

public:
  FluxFcnSGExactJacRoeKE3D(double gg, IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnExactJacRoeKE3D(ioData, gg, varFcnSGKE, tp) {}
  ~FluxFcnSGExactJacRoeKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  );

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  )
  {
    std::cout << "\n !!! FluxFcnSGExactJacRoeKE3D::computeDerivative (11 arg.) is not implemented !!!\n\n";
    exit(1);
  }

};

//------------------------------------------------------------------------------

class FluxFcnSGWallKE3D : public FluxFcnWallKE3D {

public:

  FluxFcnSGWallKE3D(IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnWallKE3D(varFcnSGKE, tp) {}
  ~FluxFcnSGWallKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGWallKE3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

};

//------------------------------------------------------------------------------

class FluxFcnSGPorousWallKE3D : public FluxFcnPorousWallKE3D {

public:

  FluxFcnSGPorousWallKE3D(IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnPorousWallKE3D(varFcnSGKE, tp) {}
  ~FluxFcnSGPorousWallKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

//  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGPorousWallKE3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

};

//------------------------------------------------------------------------------

class FluxFcnSGGhidagliaKE3D : public FluxFcnGhidagliaKE3D {

public:

  FluxFcnSGGhidagliaKE3D(IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaKE3D(varFcnSGKE, tp) {}
  ~FluxFcnSGGhidagliaKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGOutflowKE3D : public FluxFcnOutflowKE3D {

public:

  FluxFcnSGOutflowKE3D(IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnOutflowKE3D(varFcnSGKE, tp) {}
  ~FluxFcnSGOutflowKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGOutflowKE3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

};

//------------------------------------------------------------------------------

class FluxFcnSGDirectStateInflowKE3D : public FluxFcnDirectStateInflowKE3D {

public:

  FluxFcnSGDirectStateInflowKE3D(IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnDirectStateInflowKE3D(varFcnSGKE, tp) {}
  ~FluxFcnSGDirectStateInflowKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
//  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGDirectStateInflowKE3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

};

//------------------------------------------------------------------------------

class FluxFcnSGMassFlowInflowKE3D : public FluxFcnMassFlowInflowKE3D {

public:

  FluxFcnSGMassFlowInflowKE3D(IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnMassFlowInflowKE3D(varFcnSGKE, tp) {}
  ~FluxFcnSGMassFlowInflowKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
//  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGMassFlowInflowKE3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

};

//------------------------------------------------------------------------------

class FluxFcnSGDirectStateOutflowKE3D : public FluxFcnDirectStateOutflowKE3D {

public:

  FluxFcnSGDirectStateOutflowKE3D(IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnDirectStateOutflowKE3D(varFcnSGKE, tp) {}
  ~FluxFcnSGDirectStateOutflowKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
//  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGDirectStateOutflowKE3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

};

//------------------------------------------------------------------------------

class FluxFcnSGMassFlowOutflowKE3D : public FluxFcnMassFlowOutflowKE3D {

public:

  FluxFcnSGMassFlowOutflowKE3D(IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnMassFlowOutflowKE3D(varFcnSGKE, tp) {}
  ~FluxFcnSGMassFlowOutflowKE3D() { vf = 0; }

  void compute(double, double, double *, double, double *, double *, double *, bool);
//  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *vl, double *dvl, double *vr, double *dvr,
    double dmach, double *f, double *df
    , bool useLimiter = false
  )
  {
    std::cout << "\n !!! FluxFcnSGMassFlowOutflowKE3D::computeDerivative (14 arg.) is not implemented !!!\n\n";
    exit(1);
  }

  //--- Sensitivity Analysis Function
  virtual void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv,
    double *v, double *ub, double *dub, double *f, double *df
  );

};

//------------------------------------------------------------------------------

class FluxFcnSGRoeKEturb3D : public FluxFcnRoeKEturb3D {

public:
  FluxFcnSGRoeKEturb3D(double gg, IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnRoeKEturb3D(ioData, gg, varFcnSGKE, tp) {}
  ~FluxFcnSGRoeKEturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobians(double, double, double *, double, double *, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGWallKEturb3D : public FluxFcnWallKEturb3D {

public:

  FluxFcnSGWallKEturb3D(IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnWallKEturb3D(varFcnSGKE, tp) {}
  ~FluxFcnSGWallKEturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGPorousWallKEturb3D : public FluxFcnPorousWallKEturb3D {

public:

  FluxFcnSGPorousWallKEturb3D(IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnPorousWallKEturb3D(varFcnSGKE, tp) {}
  ~FluxFcnSGPorousWallKEturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGGhidagliaKEturb3D : public FluxFcnGhidagliaKEturb3D {

public:
  FluxFcnSGGhidagliaKEturb3D(IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnGhidagliaKEturb3D(varFcnSGKE, tp) {}
  ~FluxFcnSGGhidagliaKEturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGOutflowKEturb3D : public FluxFcnOutflowKEturb3D{

public:

  FluxFcnSGOutflowKEturb3D(IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnOutflowKEturb3D(varFcnSGKE, tp) {}
  ~FluxFcnSGOutflowKEturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGDirectStateInflowKEturb3D : public FluxFcnDirectStateInflowKEturb3D {

public:
  FluxFcnSGDirectStateInflowKEturb3D(IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnDirectStateInflowKEturb3D(varFcnSGKE, tp) {}
  ~FluxFcnSGDirectStateInflowKEturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGMassFlowInflowKEturb3D : public FluxFcnMassFlowInflowKEturb3D {

public:
  FluxFcnSGMassFlowInflowKEturb3D(IoData &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnMassFlowInflowKEturb3D(varFcnSGKE, tp) {}
  ~FluxFcnSGMassFlowInflowKEturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGDirectStateOutflowKEturb3D : public FluxFcnDirectStateOutflowKEturb3D {

public:
  FluxFcnSGDirectStateOutflowKEturb3D(IoData  &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnDirectStateOutflowKEturb3D(varFcnSGKE, tp) {}
  ~FluxFcnSGDirectStateOutflowKEturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

class FluxFcnSGMassFlowOutflowKEturb3D : public FluxFcnMassFlowOutflowKEturb3D {

public:
  FluxFcnSGMassFlowOutflowKEturb3D(IoData  &ioData, VarFcnSGKE *varFcnSGKE, Type tp = CONSERVATIVE) :
    FluxFcnMassFlowOutflowKEturb3D(varFcnSGKE, tp) {}
  ~FluxFcnSGMassFlowOutflowKEturb3D() { vf = 0; }

  void compute(double, double, double *n, double nv, double *vl, double *vr, double *f, bool) {}
  void computeJacobian(double, double, double *, double, double *, double *, double *, bool);

};

//------------------------------------------------------------------------------

#endif

