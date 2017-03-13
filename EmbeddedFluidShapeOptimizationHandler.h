#ifndef _EMB_FLUID_SHAPE_OPTIMIZATION_HANDLER_H_
#define _EMB_FLUID_SHAPE_OPTIMIZATION_HANDLER_H_

#include <ImplicitEmbeddedCoupledTsDesc.h>
#include <string>

class IoData;
class Domain;
class GeoSource;
//class MeshMotionSolver;
//class RigidMeshMotionHandler;

template<int dim, int neq>       class MatVecProd;
template<int dim, class Scalar2> class KspPrec;
template<class Scalar, int dim>  class DistSVec;
template<int dim>                class ImplicitEmbeddedCoupledTsDesc;

/*
 #ifndef _KSPSLVR_TMPL_
 #define _KSPSLVR_TMPL_
 template<class VecType, class MvpOp, class PrecOp, class IoOp, class ScalarT = double> class KspSolver;
 #endif
*/

//------------------------------------------------------------------------------
template<int dim>
class EmbeddedFluidShapeOptimizationHandler : public ImplicitEmbeddedCoupledTsDesc<dim> {

private:

  //  Domain *domain;
  MatVecProd<dim,dim> *mvp;
  KspPrec<dim> *pc;
  KspSolver<DistEmbeddedVec<double,dim>, MatVecProd<dim,dim>, KspPrec<dim>, Communicator> *ksp;
  double steadyTol;

private:

  double xmach;
  double alprad;
  double teta;
  double DFSPAR[3];

  char *dXdSb_file;

  //DistEmbeddedVec?????
  DistVec<double> Pin;

  DistSVec<double,3> *load;
  DistSVec<double,3> *dLoad;
  DistSVec<double,3> *dLoadref;

  DistSVec<double,dim> dFdS;

DistSVec<double,dim> dFdS_debug;   
DistSVec<double,dim> difference;

  DistSVec<double,3> dXdS;
  DistSVec<double,3> *X_;

  /* DistSVec<double,3>  Xc; */
  /* DistSVec<double,3> *Xm; */
  /* DistSVec<double,3> *Xp; */

  DistSVec<double,3> *Lp;
  DistSVec<double,3> *Lm;

  DistVec<double> *A_;
  /* DistVec<double> *Ap; */
  /* DistVec<double> *Am; */

  DistSVec<double,dim> Flux;
  DistSVec<double,dim> FluxFD;
  DistSVec<double,dim> *Fp;
  DistSVec<double,dim> *Fm;

  DistSVec<double,dim> *Up;
  DistSVec<double,dim> *Um;
  DistSVec<double,dim> dUdS;

  double reynolds0;
  double kenergy0;

  double length;
  double surface;

  int numLocSub;
  int actvar;
  int step;

  FILE* outFile;

public:

  /// Constructor
  /// \param[in] ioData  Reference to an 'IoData' object.
  /// \param[in] geoSource  Reference to a 'GeoSource' object.
  /// \param[in] dom  Pointer to the domain.
  /// \note The pointer 'dom' is passed to ImplicitCoupledTsDesc
  /// and stored in the member variable domain.
  /// It seems redundant (UH - 08/10)
  EmbeddedFluidShapeOptimizationHandler
  (
    IoData &ioData,
    GeoSource &geoSource,
    Domain *dom//,
  );

  ~EmbeddedFluidShapeOptimizationHandler();

  int fsoHandler(IoData &, 
		 DistSVec<double,dim> &);

  void fsoInitialize(IoData &ioData, 
		     DistSVec<double,dim> &U);
  
  void fsoSetUpLinearSolver(IoData &, 
			    DistSVec<double,3> &, 
			    DistVec<double> &, 			    
			    DistSVec<double,dim> &,
			    DistSVec<double,dim> &);

  void fsoRestartBcFluxs(IoData &);

  void fso_on_sensitivityMesh(IoData &,  DistSVec<double,dim> &);
  void fso_on_sensitivityMach(IoData &,  DistSVec<double,dim> &);
  void fso_on_sensitivityBeta(IoData &,  DistSVec<double,dim> &);
  void fso_on_sensitivityAlpha(IoData &, DistSVec<double,dim> &);

  void fsoComputeDerivativesOfFluxAndSolution(IoData &,
					      DistSVec<double,3> &, 
					      DistVec<double> &, 
					      DistSVec<double,dim> &, 
					      bool=false);

  void fsoSemiAnalytical(IoData &, 
			 DistSVec<double,3> &,
			 DistVec<double> &,
			 DistSVec<double,dim> &,
			 DistSVec<double,dim> &);

  void fsoAnalytical(IoData &, 
		     DistSVec<double,3> &,
		     DistVec<double> &,
		     DistSVec<double,dim> &,
		     DistSVec<double,dim> &);

  void fsoLinearSolver(IoData &, 
		       DistSVec<double,dim> &, 
		       DistSVec<double,dim> &,
		       bool);

  void fsoComputeSensitivities(IoData &,
			       const char *, 
			       const char *, 
			       DistSVec<double,3> &, 
			       DistSVec<double,dim> &);

  void fsoGetEfforts(IoData &, 
		     DistSVec<double,3> &, 
		     DistSVec<double,dim> &, 
		     Vec3D &, Vec3D &, Vec3D &);

  void fsoGetDerivativeOfEffortsAnalytical(IoData &ioData,  DistSVec<double,3>   &X,
					   DistSVec<double,dim> &U, DistSVec<double,dim> &dU,
					   Vec3D &dForces, Vec3D &dMoments, Vec3D &dL);

  void fsoGetDerivativeOfEffortsFiniteDifference(IoData &ioData, DistSVec<double,3> &X,
						 DistSVec<double,dim> &U, DistSVec<double,dim> &dU,
						 Vec3D &dForces, Vec3D &dMoments, Vec3D &dL);

  

  void fsoPrintTextOnScreen(const char *);

  void fsaPrintTextOnScreen(const char *);

  bool getdXdSb(int);

  Communicator* getComm(){return this->com;};//TODO HACK

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <EmbeddedFluidShapeOptimizationHandler.C>
#endif

#endif


