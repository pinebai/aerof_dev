#ifndef _EDGE_H_
#define _EDGE_H_

#ifdef OLD_STL
#include <map.h>
#else
#include <map>
#include <vector>
using std::map;
#endif

#if defined(__KCC)
#include <utility.h>
#else
#include <utility>
using std::pair;
#endif

// Included (MB)
#include <cstdio>

//#include <ProgrammedBurn.h>
//#include <HigherOrderMultiFluid.h>
//#include <HigherOrderFSI.h>
//#include <ErrorHandler.h>
//#include <LevelSet/LevelSetStructure.h>

class FluidSelector;
class VarFcn;
class RecFcn;
class FluxFcn;
class ElemSet;
class GeoState;
class FemEquationTerm;
class TimeLowMachPrec;
class LevelSetStructure;
class ProgrammedBurn;
class HigherOrderMultiFluid;
class HigherOrderFSI;
class ErrorHandler;

struct V6NodeData;
struct Vec3D;

#ifndef _NDGRAD_TMPL_
#define _NDGRAD_TMPL_
template<int dim, class Scalar = double> class NodalGrad;
#endif

template<int dim> class EdgeGrad;
template<int dim> class ExactRiemannSolver;
template<class Scalar> class Vec;
template<class Scalar, int dim> class SVec;
template<class Scalar, int dim> class GenMat;
template<class Scalar, int dim, int dim2> class RectangularSparseMat;

//------------------------------------------------------------------------------

class EdgeSet {

public:

  enum MultifluidRiemannNormal { MF_RIEMANN_NORMAL_MESH, MF_RIEMANN_NORMAL_REAL, MF_RIEMANN_NORMAL_LEGACYMESH };

private:

  typedef pair<int, int> Pair;

#ifdef OLD_STL
  typedef map<Pair, int, less<Pair> > MapPair;
#else
  typedef map<Pair, int> MapPair;
#endif

  MapPair *mp;

  int numEdges;
  bool sampleMesh;
  int numSampledEdges, numTwoLayerEdges;

  int (*ptr)[2];
  bool *masterFlag;

  double* edgeLength;

  ProgrammedBurn* programmedBurn;

  HigherOrderMultiFluid* higherOrderMF;
  HigherOrderFSI* higherOrderFSI;

  MultifluidRiemannNormal mfRiemannNormal;
  ErrorHandler* errorHandler;

  LevelSetStructure* triangulatedLSS;

public:


  EdgeSet();
  ~EdgeSet();

  int find(int, int);
  int findOnly(int, int) const;

  void createPointers(Vec<int> &);

  void setMultiFluidRiemannNormal(MultifluidRiemannNormal);

  template<int dim>
  void computeTimeStep(FemEquationTerm *, VarFcn *, GeoState &,
		       SVec<double,3> &, SVec<double,dim> &, Vec<double> &,
                       Vec<double> &, TimeLowMachPrec &, LevelSetStructure* lss=0);

  template<int dim>
  void computeTimeStep2(FemEquationTerm *, VarFcn *, GeoState &,
		       SVec<double,3> &, SVec<double,dim> &, Vec<double> &,
                       Vec<double> &, TimeLowMachPrec &, Vec<double>&);
 
  template<int dim>
  void computeTimeStep(VarFcn *, GeoState &, SVec<double,dim> &, Vec<double> &,
                       TimeLowMachPrec &, Vec<int> &, int,Vec<double>*);

  template<int dim>
  int computeFiniteVolumeTerm(int*, Vec<double> &, FluxFcn**, RecFcn*, ElemSet&, GeoState&,
                              SVec<double,3>&, SVec<double,dim>&, NodalGrad<dim>&, EdgeGrad<dim>*,
			      SVec<double,dim>&, SVec<int,2>&, int, int);
  
  template<int dim>
  int computeThinLayerViscousFiniteVolumeTerm(int* locToGlobNodeMap,
                                     VarFcn* varFcn,
                                     class FemEquationTerm*,
                                     GeoState& geoState, SVec<double,3>& X,
                                     SVec<double,dim>& V,
                                     SVec<double,dim>& fluxes);
  
  template<int dim>
    int computeViscousFiniteVolumeTerm(int* locToGlobNodeMap,
				   VarFcn* varFcn,
				   FemEquationTerm *fet,
				   GeoState& geoState, SVec<double,3>& X,
				   SVec<double,dim>& V,
				   SVec<double,dim> &dX,
				   SVec<double,dim> &dY,
				   SVec<double,dim> &dZ,
				   SVec<double,dim>& fluxes,
                                   LevelSetStructure* = NULL);

  template<int dim,class Scalar,int neq>
  int computeJacobianThinLayerViscousFiniteVolumeTerm(int* locToGlobNodeMap,
                                     VarFcn* varFcn,
                                     class FemEquationTerm*,
                                     GeoState& geoState, SVec<double,3>& X,
                                     SVec<double,dim>& V,
                                     Vec<double>& ctrlVol,
/*                                     SVec<double,3>& faceJacX,
                                     SVec<double,3>& faceJacY,
                                     SVec<double,3>& faceJacZ,
                                     bool* boundaryFlag,
*/
                                     GenMat<Scalar,neq>& A);

  template<int dim>
  int computeFiniteVolumeTermRestrict(int*, Vec<double> &, FluxFcn**, RecFcn*, ElemSet&, GeoState&,
                              SVec<double,3>&, SVec<double,dim>&, NodalGrad<dim>&, EdgeGrad<dim>*,
			      SVec<double,dim>&, SVec<int,2>&, int, int);

  template<int dim, int dimLS>
  int computeFiniteVolumeTerm(ExactRiemannSolver<dim>&, int*,
                              FluxFcn**, RecFcn*, ElemSet&, GeoState&, SVec<double,3>&,
                              SVec<double,dim>&, Vec<int> &, FluidSelector &,
                              NodalGrad<dim>&, EdgeGrad<dim>*,
			      SVec<double,dimLS>& phi,
                              NodalGrad<dimLS>&, EdgeGrad<dimLS>*,
                              SVec<double,dim>&, int,
                              SVec<int,2>&, int, int);

  template<int dim, int dimLS>
  int computeFiniteVolumeTerm(ExactRiemannSolver<dim>&, int*,
                              FluxFcn**, RecFcn*, ElemSet&, GeoState&, SVec<double,3>&,
                              SVec<double,dim>&, SVec<double,dim>&, SVec<double,dim>&,
                              LevelSetStructure&, bool, Vec<int> &, int, FluidSelector &,
                              NodalGrad<dim>&, EdgeGrad<dim>*,
			      SVec<double,dimLS>& phi,
                              NodalGrad<dimLS>&, EdgeGrad<dimLS>*,
                              SVec<double,dim>&, int,
                              SVec<int,2>&, int, int);

  /** compute flux for Riemann based FSI*/
  template<int dim>
  int computeFiniteVolumeTerm(ExactRiemannSolver<dim>&, int*,
                              FluxFcn**, RecFcn*, ElemSet&, GeoState&, SVec<double,3>&,
                              SVec<double,dim>&, SVec<double,dim>&, SVec<double,dim>&, LevelSetStructure &,
                              bool, Vec<int>&, int, NodalGrad<dim>&, EdgeGrad<dim>*,
                              SVec<double,dim>&, int,
                              SVec<int,2>&, int, int);

  template<int dim>
	  int computeFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann, int* locToGlobNodeMap,
											FluxFcn** fluxFcn, RecFcn* recFcn,
											ElemSet& elems, GeoState& geoState, SVec<double,3>& X,
											SVec<double,dim>& V, SVec<double,dim>& Vstarij,
											SVec<double,dim>& Vstarji, SVec<double,dim> &Vext,
											LevelSetStructure &LSS,
											Vec<int> &fluidId, int Nriemann,
											NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
											SVec<double,dim>& fluxes, int it,
											SVec<int,2>& tag, int failsafe, int rshift);
											
  template<int dim>
  int computeFiniteVolumeTerm(ExactRiemannSolver<dim>&, int*,
                              FluxFcn**, RecFcn*, ElemSet&, GeoState&, SVec<double,3>&,
                              SVec<double,dim>&, SVec<double,dim>&, SVec<double,dim>&,
                              Vec<int>&, Vec<int>&, LevelSetStructure &, bool, Vec<int>&, 
                              int, double, double, NodalGrad<dim>&, 
                              EdgeGrad<dim>*, SVec<double,dim>&, int,
                              SVec<int,2>&, int, int, V6NodeData (*v6Data)[2]=NULL); 

  template<int dim, int dimLS>
  void computeFiniteVolumeTermLS(FluxFcn**, RecFcn*, RecFcn*, ElemSet&, GeoState&, SVec<double,3>&,
                               SVec<double,dim>&,Vec<int>& fluidId, 
  			       NodalGrad<dim>&,   EdgeGrad<dim>*,
			       NodalGrad<dimLS>&, EdgeGrad<dimLS>*,
                               SVec<double,dimLS>&, SVec<double,dimLS>&, LevelSetStructure* =0, int order = 1);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, GeoState &,
                                       Vec<double> &, SVec<double,3> &, Vec<double> &,
                                       SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<int dim, class Scalar, int neq>
  void computeJacobianFiniteVolumeTerm(FluxFcn **, GeoState &,
                               Vec<double> &, SVec<double,3> &, Vec<double> &,
                               SVec<double,dim> &, GenMat<Scalar,neq> &,
                               int *);

  template<int dim, class Scalar, int neq, int dimLS>
  void computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>&, FluxFcn **, GeoState &,
                                       NodalGrad<dim> &, NodalGrad<dimLS> &,
                                       SVec<double,3> &, Vec<double> &,
                                       SVec<double,dim> &, GenMat<Scalar,neq> &,
                                       FluidSelector &, Vec<int> &);
  template<int dim, class Scalar, int neq, int dimLS>
  void computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>&, FluxFcn **, GeoState &,
                               NodalGrad<dim> &, NodalGrad<dimLS> &,
                               SVec<double,3> &, Vec<double> &,
                               SVec<double,dim> &, GenMat<Scalar,neq> &,
                               FluidSelector &, Vec<int> &, int * );
  
  template<class Scalar,int dim,int neq>
  void computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>&,
                              FluxFcn**, GeoState&, SVec<double,3>&,
                              SVec<double,dim>&, Vec<double>&, LevelSetStructure &,
                              Vec<int>&, int,
                              GenMat<Scalar,neq>&,Vec<double>& irey);

  template<class Scalar,int dim,int neq>
  void computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>&,
													FluxFcn**, GeoState&, SVec<double,3>&,
													SVec<double,dim>&, Vec<double>&, LevelSetStructure &,
													Vec<int>&, int,
													GenMat<Scalar,neq>&);

  template<class Scalar,int dim, int dimLS,int neq>
  void computeJacobianFiniteVolumeTerm(ExactRiemannSolver<dim>& riemann,int*,
                                     FluxFcn** fluxFcn, 
                                     GeoState& geoState, SVec<double,3>& X,
                                     SVec<double,dim>& V, 
                                     LevelSetStructure& LSS, Vec<int> &fluidId,
                                     int Nriemann, FluidSelector &fluidSelector,
                                     NodalGrad<dimLS>& ngradLS,Vec<double>&,GenMat<Scalar,neq>& A);

  template<class Scalar, int dim, int dimLS>
    void computeJacobianFiniteVolumeTermLS(RecFcn* recFcn, RecFcn* recFcnLS,
					   GeoState& geoState, SVec<double,3>& X,
					   SVec<double,dim>& V, NodalGrad<dim>& ngrad,
					   NodalGrad<dimLS> &ngradLS,
					   EdgeGrad<dim>* egrad,
					   Vec<double> &ctrlVol, SVec<double,dimLS>& Phi,
					   GenMat<Scalar,dimLS> &A,LevelSetStructure* LSS);
/*  template<int dim>
  void RiemannJacobianGasTait(int i, int j,
                              SVec<double,dim> &V, double Phii, double Phij,
		              double *nphi,
                              double *normal, double normalVel, VarFcn *varFcn,
                              FluxFcn** fluxFcn, double *dfdUi, double *dfdUj);
*/

  template<int dimLS>
  void TagInterfaceNodes(int lsdim, Vec<int> &Tag, SVec<double,dimLS> &Phi,LevelSetStructure *LSS=0);
  template<int dimLS>
  void pseudoFastMarchingMethodInitialization(SVec<double,3>& X,
		Vec<int> &Tag, SVec<double,dimLS> &d2wall, 
		Vec<int> &sortedNodes, int &nSortedNodes,
		LevelSetStructure *LSS=0);

  void setMasterFlag(bool *flag) { masterFlag = flag; }
  bool *getMasterFlag() const { return masterFlag; }
  int (*getPtr() const)[2] { return ptr; }
  int size() const { return numEdges; }

  void updateLength(SVec<double,3>& X);
  double length(int iedge) const {
    if(edgeLength && iedge<numEdges) return(edgeLength[iedge]);
    else return(0.0);
  }
  double* viewEdgeLength() { return(edgeLength); }

// Included (MB)
  template<int dim>
  void computeDerivativeOfFiniteVolumeTerm(Vec<double> &, Vec<double> &, FluxFcn**, RecFcn*, ElemSet&, GeoState&, SVec<double,3>&, SVec<double,3>&,
                                           SVec<double,dim>&, SVec<double,dim>&, NodalGrad<dim>&, EdgeGrad<dim>*, double,
                                           SVec<double,dim>&);
  
  //Direct senstitivity contribution of the viscous term in ALE simulation
  template<int dim>
  void computeDerivativeOfFiniteVolumeTerm(RectangularSparseMat<double,dim,dim> *,
                                           RectangularSparseMat<double,dim,dim> *dFluxdddy,
                                           RectangularSparseMat<double,dim,dim> *dFluxdddz,
                                           RectangularSparseMat<double,3,dim> *dFluxdX,
                                           RectangularSparseMat<double,3,dim> *dFluxdEdgeNorm,
                                           ElemSet&, GeoState&, SVec<double,3>&, 
                                           NodalGrad<dim>&, EdgeGrad<dim>*, 
                                           SVec<double,dim>& dddx,
                                           SVec<double,dim>& dddy,
                                           SVec<double,dim>& dddz,
                                           Vec<Vec3D>& dNormal,
                                           SVec<double,dim>&);
  
  // For Adjoint sensitivity contibution of the viscous term
  template<int dim>
  void computeTransposeDerivativeOfFiniteVolumeTerm(RectangularSparseMat<double,dim,dim> *,
                                           RectangularSparseMat<double,dim,dim> *dFluxdddy,
                                           RectangularSparseMat<double,dim,dim> *dFluxdddz,
                                           RectangularSparseMat<double,3,dim> *dFluxdX,
                                           RectangularSparseMat<double,3,dim> *dFluxdEdgeNorm,
                                           SVec<double,dim>& dFluxes,
                                           NodalGrad<dim>& ngrad,
                                           EdgeGrad<dim>* egrad,
                                           ElemSet& elems,
                                           GeoState& geoState,
                                           SVec<double,3>& dX,
                                           SVec<double,dim>& dddx,
                                           SVec<double,dim>& dddy,
                                           SVec<double,dim>& dddz,
                                           Vec<Vec3D>& dNormal);                          
  
  //Direct sensitivity contribution of the viscous term in Embedded simulation
  template<int dim>
  void computeDerivativeOfFiniteVolumeTerm(FluxFcn** fluxFcn, RecFcn* recFcn,
					   GeoState& geoState, SVec<double,3>& X, LevelSetStructure &LSS,
					   bool linRecAtInterface, Vec<int> &fluidId, 
					   ExactRiemannSolver<dim>& riemann, int Nriemann,
					   NodalGrad<dim>& ngrad, EdgeGrad<dim>* egrad,
					   double dMach, SVec<double,dim>& V, SVec<double,dim>& dFluxes);
   
  template<int dim>
  void computeDerivativeOperatorsOfFiniteVolumeTerm(Vec<double> &irey, Vec<double> &dIrey, FluxFcn** fluxFcn, RecFcn* recFcn,
              ElemSet& elems, GeoState& geoState, SVec<double,3>& X, SVec<double,dim>& V, NodalGrad<dim>& ngrad,
              EdgeGrad<dim>* egrad, double dMach, 
              RectangularSparseMat<double,3,dim> &dFluxdEdgeNorm,
              RectangularSparseMat<double,3,dim> &dFluxdX,
              RectangularSparseMat<double,dim,dim> &dFluxdddx,
              RectangularSparseMat<double,dim,dim> &dFluxdddy,
              RectangularSparseMat<double,dim,dim> &dFluxdddz);
 
  template<int dim>
  void computeDerivativeOfTimeStep(FemEquationTerm *, VarFcn *, GeoState &,
                              SVec<double,3> &, SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &,
			      Vec<double> &, Vec<double> &, double, TimeLowMachPrec &);
  int checkReconstructedValues(int i, int j, double *Vi, double *Vj, VarFcn *vf,
			       int *locToGlobNodeMap, int failsafe, SVec<int,2> &tag,
                               double *originalVi = 0, double *originalVj = 0,
                               int IDi = 0, int IDj = 0, bool iact = true, bool jact = true);

  void computeCharacteristicEdgeLength(SVec<double,3> &, double&, double&, double&, int&,
                                       const double, const double, const double,
                                       const double, const double, const double);

  void attachProgrammedBurn(ProgrammedBurn*);
  void attachHigherOrderMultiFluid(HigherOrderMultiFluid*);
  void attachHigherOrderFSI(HigherOrderFSI*);

  void assignErrorHandler(ErrorHandler* in){errorHandler = in;}

  void computeConnectedEdges(const std::vector<int> &);
  std::vector<int> edgesConnectedToSampleNode;	// for Gappy ROM
  std::vector<int> edgesTwoLayersSampleNode;	// for Gappy ROM
  int getNumSampledEdges() {return numSampledEdges;}
  int getNumTwoLayersEdges() {return numTwoLayerEdges;}
  void computeGlobalConnectedEdges(const std::vector<int> &globalNeighborNodes,
				   const int *locToGlobNodeMap) ;
  template<int dim>
  void ComputeAndAddActuatorDiskSourceTerm(Vec3D structureNormal, Vec3D fluidNormal,Vec3D normalDir,double controlVolumeArea,bool invertedEdge, double pressureJumpValue,double VelocityReconstructed[3] ,double (&fluxi)[dim],int i,bool isModifiedSourceTerm, double gamma);
  template<int dim>
  void ComputeActuatorVelocity(double (&VelocityReconstructed)[3],int reconstructionMethod,double (&Vi)[2*dim],double (&Vj)[2*dim], double alpha,double (&ddVij)[dim],double (&ddVji)[dim]);


  void attachTriangulatedInterfaceLSS(LevelSetStructure*);
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <Edge.C>
#endif

#endif
