#ifndef _FLUID_SHAPE_OPTIMIZATION_HANDLER_H_
#define _FLUID_SHAPE_OPTIMIZATION_HANDLER_H_

#include <ImplicitCoupledTsDesc.h>

#include <string>

#define Deg2Rad         0.01745329251994329576923
#define Rad2Deg        57.29577951308232087679815
#define perDeg2perRad  57.29577951308232087679815
#define perRad2perDeg   0.01745329251994329576923

class IoData;
class Domain;
class GeoSource;
class MeshMotionSolver;
//class StructExc;
class RigidMeshMotionHandler;

template<int dim, int neq> class MatVecProd;
template<int dim, class Scalar2> class KspPrec;
template<class Scalar, int dim> class DistSVec;
template<int dim> class ImplicitCoupledTsDesc;

#ifndef _KSPSLVR_TMPL_
#define _KSPSLVR_TMPL_
template<class VecType, class MvpOp, class PrecOp, class IoOp, class ScalarT = double> class KspSolver;
#endif

//------------------------------------------------------------------------------
template<int dim>
class FluidShapeOptimizationHandler : public ImplicitCoupledTsDesc<dim> {

private:

  Domain *domain;
//  TsSolver<ImplicitCoupledTsDesc<dim> > *tsSolver;
  MeshMotionSolver *mms;

  // UH (08/10) This pointer is never used.
  //StructExc *strExc;

  MatVecProd<dim,dim> *mvp;
  MatVecProd<dim,dim> *dRdX;
  KspPrec<dim> *pc;
  KspSolver<DistSVec<double,dim>, MatVecProd<dim,dim>, KspPrec<dim>, Communicator> *ksp;
  double steadyTol;

private:

  int step;
  int actvar;
  int numLocSub;

  double reynolds0;
  double kenergy0; 

  double length;
  double surface;
  double xmach;
  double alprad;
  double teta;
  double DFSPAR[3];


  struct{
	     bool mach;
		 bool alpha;
		 bool beta;
		 bool mesh;
  } senstype;

  DistVec<double> Pin;
  DistVec<double> dAdS;
  DistVec<double> *Ap;
  DistVec<double> *Am;
  
  DistSVec<double,3> p;
  DistSVec<double,3> dPdS;
  DistSVec<double,3> dXdS;//derivative of the mesh motion at all points
  DistSVec<double,3> dXdSb;//derivative of the mesh mesh motion at the interface
  DistSVec<double,3> lambdaSDisp;
  DistSVec<double,3> dfaX;
  DistSVec<double,3> dXb;
  DistSVec<double,3> Xc;
  DistSVec<double,3> *Xp;
  DistSVec<double,3> *Xm;
  DistSVec<double,3> *Lp;
  DistSVec<double,3> *Lm;
  DistSVec<double,3> *Z;
  DistSVec<double,3> *load;
  DistSVec<double,3> *dLoad;
  DistSVec<double,3> *dLoadref;

  DistSVec<double,dim> Flux;
  DistSVec<double,dim> FluxFD;
  DistSVec<double,dim> *Fp;
  DistSVec<double,dim> *Fm;
  DistSVec<double,dim> *Fp_inviscid;
  DistSVec<double,dim> *Fm_inviscid;
  DistSVec<double,dim> *Fp_viscous;
  DistSVec<double,dim> *Fm_viscous;
  DistSVec<double,dim> dFdS;
  DistSVec<double,dim> dFdS_viscous;
  DistSVec<double,dim> dFdS_inviscid;
  DistSVec<double,dim> dFdSref;
  DistSVec<double,dim> temp;
  DistSVec<double,dim> *Up;
  DistSVec<double,dim> *Um;
  DistSVec<double,dim> dUdS;
  DistSVec<double,dim> lambdaU;
  DistSVec<double,dim> dfaU;
  DistSVec<double,dim> Uc;

  DistSVec<double,3> Xplus;
  DistSVec<double,3> Xminus;
  DistSVec<double,3> dX;
  DistSVec<double,3> lambdaX;

  DistSVec<double,dim> dddx;  // nodal gradients or adjoint vectors
  DistSVec<double,dim> dddy;
  DistSVec<double,dim> dddz;

  DistVec<Vec3D> dEdgeNorm;      // derivative of edge normal or adjoint vectors
  DistVec<Vec3D> dFaceNorm;      // derivative of face normal or adjoint vectors
  DistVec<double> dFaceNormVel;  // derivative of face normal velocity or adjoint vectors

  DistSVec<double,6> dR;        // derivative of least square gradient coefficient or adjoint vectors
  DistSVec<double,3> dGradP;    // nodal pressure gradient

  FILE* outFile;
  void setDFSPAR(IoData &);  
  
public:

  /// Constructor
  /// \param[in] ioData  Reference to an 'IoData' object.
  /// \param[in] geoSource  Reference to a 'GeoSource' object.
  /// \param[in] dom  Pointer to the domain.
  /// \note The pointer 'dom' is passed to ImplicitCoupledTsDesc
  /// and stored in the member variable domain.
  /// It seems redundant (UH - 08/10)
  FluidShapeOptimizationHandler
  (
    IoData &ioData,
    GeoSource &geoSource,
    Domain *dom//,
  );

  ~FluidShapeOptimizationHandler();

  /// \note This function is implemented but never called.
  void fsoOutput1D(const char *, DistVec<double> &);

  /// \note This function is implemented but never called.
  void fsoOutput3D(const char *, DistSVec<double,3> &);

  /// \note This function is implemented but never called.
  void fsoOutputDimD(const char *, DistSVec<double,dim> &);

  void fsoPrintTextOnScreen(const char *);

  //TODO fix this does completely the same as the function above
  void fsaPrintTextOnScreen(const char *);

  void fsoRestartBcFluxs(IoData &);

  //TODO HACK
  void fsaRestartBcFluxs(IoData &);

  void fsoGetEfforts(IoData &, DistSVec<double,3> &, DistSVec<double,dim> &, Vec3D &, Vec3D &, Vec3D &);

  void fsoGetDerivativeOfEffortsFiniteDifference(IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double>&, DistSVec<double,dim> &, DistSVec<double,dim> &, Vec3D &, Vec3D &,Vec3D &);

  void fsoGetDerivativeOfEffortsAnalytical(bool, IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &, Vec3D &, Vec3D &, Vec3D &);

  void fsoGetDerivativeOfEffortsWRTStateAndMeshPositionAnalytical(IoData &,
		  Vec3D &,
		  Vec3D &,
		  Vec3D &,
		  DistSVec<double,3> &,
		  DistSVec<double,dim> &,
		  DistSVec<double,3> &,
		  DistSVec<double,dim> &);

  /// \note This function is implemented but never called.
  void fsoGetDerivativeOfLoadFiniteDifference(IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,3> &, DistSVec<double,3> &);

  /// \note This function is implemented but never called.
  void fsoGetDerivativeOfLoadAnalytical(bool, IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,3> &, DistSVec<double,3> &);
  void fsoGetTransposeDerivativeOfLoadAnalytical(IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,dim> &);
  void fsoSemiAnalytical(IoData &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void fsoAnalytical(bool, IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void fsoApply_dFdXtranspose(DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,3> &);
  void fsoSetUpLinearSolver(IoData &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void fsoSetUpAdjointLinearSolver(IoData &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void fsoLinearSolver(IoData &, DistSVec<double,dim> &, DistSVec<double,dim> &, bool=false);
  void fsoAdjointLinearSolver(IoData &, DistSVec<double,dim> &, DistSVec<double,dim> &, bool=false);
  int fsoHandler(IoData &, DistSVec<double,dim> &);
  int fsoAeroelasticHandler(IoData &, DistSVec<double,dim> &);
  void fsoComputeDerivativesOfFluxAndSolution(IoData &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, bool=false, bool=false);
  void fsoComputeSensitivities(bool, IoData &, const char *, const char *, DistSVec<double,3> &, DistSVec<double,dim> &);
  void fsoComputeAdjoint(IoData &, DistVec<double> &, DistSVec<double,3> &, DistSVec<double,dim> &, bool);
  void fsoComputeAndSendForceSensitivities(bool, IoData &, const char *, DistSVec<double,3> &, DistSVec<double,dim> &);
  void fsoInitialize(IoData &ioData, DistSVec<double,dim> &U);
  void fso_on_aeroelasticSensitivityFSI(bool, IoData &ioData, DistSVec<double,dim> &U);
  void fso_on_aeroelasticAdjointSensitivityFSI(IoData &ioData, DistSVec<double,dim> &U);
  void fso_on_sensitivityMesh(bool, IoData &ioData, DistSVec<double,dim> &U);
  void fso_on_sensitivityMach(bool, IoData &ioData, DistSVec<double,dim> &U);
  void fso_on_sensitivityAlpha(bool, IoData &ioData, DistSVec<double,dim> &U);
  void fso_on_sensitivityBeta(bool, IoData &ioData, DistSVec<double,dim> &U);

  void fso_on_AdjointSensitivityMesh(IoData &ioData, DistSVec<double,dim> &U);
//  void fso_on_AdjointSensitivityMach(IoData &ioData, DistSVec<double,dim> &U);
//  void fso_on_AdjointSensitivityAlpha(IoData &ioData, DistSVec<double,dim> &U);
//  void fso_on_AdjointSensitivityBeta(IoData &ioData, DistSVec<double,dim> &U);

  Communicator* getComm(){return this->com;};

  //TODO check if ever called
  void fsaOnlySolve(IoData &ioData);


  void Forces2Lifts(IoData &ioData, Vec3D &F, Vec3D &L)
  {
    double sin_a = sin(ioData.bc.inlet.alpha);
    double cos_a = cos(ioData.bc.inlet.alpha);
    double sin_b = sin(ioData.bc.inlet.beta);
    double cos_b = cos(ioData.bc.inlet.beta);
    L[0] =  F[0]*cos_a*cos_b + F[1]*cos_a*sin_b + F[2]*sin_a;
    L[1] = -F[0]*sin_b       + F[1]*cos_b;
    L[2] = -F[0]*sin_a*cos_b - F[1]*sin_a*sin_b + F[2]*cos_a;
  };

  void dForces2dLifts(IoData &ioData,Vec3D &F, Vec3D &dF,Vec3D &dL){
    Forces2Lifts(ioData,dF,dL);
    double sin_a = sin(ioData.bc.inlet.alpha);
    double cos_a = cos(ioData.bc.inlet.alpha);
    double sin_b = sin(ioData.bc.inlet.beta);
    double cos_b = cos(ioData.bc.inlet.beta);
    double dsin_a = cos_a*DFSPAR[1], dcos_a = -sin_a*DFSPAR[1];
    double dsin_b = cos_b*DFSPAR[2], dcos_b = -sin_b*DFSPAR[2];
    double convfac = ((ioData.sa.angleRad == ioData.sa.OFF_ANGLERAD) && (DFSPAR[1] || DFSPAR[2])) ? perRad2perDeg : 1.0;
    dL[0] += (F[0]*(dcos_a*cos_b + cos_a*dcos_b) +
              F[1]*(dcos_a*sin_b + cos_a*dsin_b) +
              F[2]*dsin_a                          )*convfac;
  
    dL[1] += (-F[0]*dsin_b + F[1]*dcos_b           )*convfac;
  
    dL[2] += (-F[0]*(dsin_a*cos_b + sin_a*dcos_b) -
              F[1]*(dsin_a*sin_b + sin_a*dsin_b) +
              F[2]*dcos_a                         )*convfac;
  };


};


//------------------------------------------------------------------------------


#ifdef TEMPLATE_FIX
#include <FluidShapeOptimizationHandler.C>
#endif

#endif


