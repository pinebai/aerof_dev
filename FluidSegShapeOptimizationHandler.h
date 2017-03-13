/*
 * FluidSegShapeOptimizationHandler.h
 *
 *  Created on: Dec 7, 2016
 *      Author: lscheuch
 */

#ifndef FLUIDSEGSHAPEOPTIMIZATIONHANDLER_H_
#define FLUIDSEGSHAPEOPTIMIZATIONHANDLER_H_



#ifndef _FLUID_SHAPE_OPTIMIZATION_HANDLER_H_
#define _FLUID_SHAPE_OPTIMIZATION_HANDLER_H_

#include <ImplicitSegTsDesc.h>

#include <string>

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
template<int dim, int neq1, int neq2>
class FluidSegShapeOptimizationHandler : public ImplicitSegTsDesc<dim,neq1,neq2> {

private:

  Domain *domain;
//  TsSolver<ImplicitCoupledTsDesc<dim> > *tsSolver;
  MeshMotionSolver *mms;

  // UH (08/10) This pointer is never used.
  //StructExc *strExc;

  MatVecProd<dim,dim> *mvp;//TODO this was the old version

  MatVecProd<dim,neq1> *mvp1;
  MatVecProd<dim,neq2> *mvp2;//TODO check if working

  MatVecProd<dim,dim> *dRdX;
  KspPrec<dim> *pc;
  KspSolver<DistSVec<double,dim>, MatVecProd<dim,dim>, KspPrec<dim>, Communicator> *ksp;
  double steadyTol;

  Communicator* getComm(){return this->com;};

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

  DistVec<double> Pin;
  DistVec<double> dAdS;
  DistVec<double> *Ap;
  DistVec<double> *Am;

  DistSVec<double,3> p;
  DistSVec<double,3> dPdS;
  DistSVec<double,3> dXdS;
  DistSVec<double,3> dXdSb;
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
  DistSVec<double,dim> dFdS;
  DistSVec<double,dim> dFdSref;
  DistSVec<double,dim> *Up;
  DistSVec<double,dim> *Um;
  DistSVec<double,dim> dUdS;
  DistSVec<double,dim> Uc;

  DistSVec<double,3> Xplus;
  DistSVec<double,3> Xminus;
  DistSVec<double,3> dX;

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
  FluidSegShapeOptimizationHandler
  (
    IoData &ioData,
    GeoSource &geoSource,
    Domain *dom//,
  );

  ~FluidSegShapeOptimizationHandler();

  /// \note This function is implemented but never called.
  void fsoOutput1D(const char *, DistVec<double> &);

  /// \note This function is implemented but never called.
  void fsoOutput3D(const char *, DistSVec<double,3> &);

  /// \note This function is implemented but never called.
  void fsoOutputDimD(const char *, DistSVec<double,dim> &);

  void fsoPrintTextOnScreen(const char *);

  //TODO fix: this is exactly the same as above
  void fsaPrintTextOnScreen(const char *);

  void fsoRestartBcFluxs(IoData &);

  //TODO this should be unneccesary
  void fsaRestartBcFluxs(IoData &);

  void fsoGetEfforts(IoData &, DistSVec<double,3> &, DistSVec<double,dim> &, Vec3D &, Vec3D &, Vec3D &);

  void fsoGetDerivativeOfEffortsFiniteDifference(IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double>&, DistSVec<double,dim> &, DistSVec<double,dim> &, Vec3D &, Vec3D &);

  void fsoGetDerivativeOfEffortsAnalytical(bool, IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &, Vec3D &, Vec3D &, Vec3D &);

  /// \note This function is implemented but never called.
  void fsoGetDerivativeOfLoadFiniteDifference(IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,3> &, DistSVec<double,3> &);

  /// \note This function is implemented but never called.
  void fsoGetDerivativeOfLoadAnalytical(bool, IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,3> &, DistSVec<double,3> &);
  void fsoGetTransposeDerivativeOfLoadAnalytical(IoData &, DistSVec<double,3> &, DistSVec<double,3> &, DistSVec<double,dim> &);
  void fsoSemiAnalytical(IoData &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void fsoAnalytical(bool, IoData &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void fsoAnalyticalTranspose(DistSVec<double,dim> &, DistSVec<double,3> &);
  void fsoSetUpLinearSolver(IoData &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void fsoLinearSolver(IoData &, DistSVec<double,dim> &, DistSVec<double,dim> &, bool=false);
  int fsoHandler(IoData &, DistSVec<double,dim> &);
  int fsoAeroelasticHandler(IoData &, DistSVec<double,dim> &);
  void fsoComputeDerivativesOfFluxAndSolution(IoData &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, bool=false, bool=false);
  void fsoComputeSensitivities(bool, IoData &, const char *, const char *, DistSVec<double,3> &, DistSVec<double,dim> &);
  void fsoComputeAndSendForceSensitivities(bool, IoData &, const char *, DistSVec<double,3> &, DistSVec<double,dim> &);
  void fsoInitialize(IoData &ioData, DistSVec<double,dim> &U);
  void fso_on_aeroelasticSensitivityFSI(bool, IoData &ioData, DistSVec<double,dim> &U);
  void fso_on_sensitivityMesh(bool, IoData &ioData, DistSVec<double,dim> &U);
  void fso_on_sensitivityMach(bool, IoData &ioData, DistSVec<double,dim> &U);
  void fso_on_sensitivityAlpha(bool, IoData &ioData, DistSVec<double,dim> &U);
  void fso_on_sensitivityBeta(bool, IoData &ioData, DistSVec<double,dim> &U);

};


//------------------------------------------------------------------------------


#ifdef TEMPLATE_FIX
#include <FluidSegShapeOptimizationHandler.C>
#endif

#endif










#endif /* FLUIDSEGSHAPEOPTIMIZATIONHANDLER_H_ */
