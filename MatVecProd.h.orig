#ifndef _MAT_VEC_PROD_H_
#define _MAT_VEC_PROD_H_

#include <IoData.h>
#include <MvpMatrix.h>
#include <DistVector.h>
#include <DistMatrix.h>
#include <complex>
#include <RestrictionMapping.h>
//#include <TsDesc.h>
#include "Dev/devtools.h"

typedef std::complex<double> bcomp;

class VarFcn;
class FluxFcn;
class DistGeoState;
class MemoryPool;
//struct Vec3D;

template<int dimLS> class LevelSet;
template<int dim> class RecFcnConstant;
template<int dim> class DistTimeState;
template<int dim> class SpaceOperator;
template<int dim, int dimLS> class MultiPhaseSpaceOperator;
template<int dim> class DistExactRiemannSolver;
template<class Scalar, int dim> class DistEmbeddedVec;
template<int dim> class PostOperator;
template<int dim> class TsOutput;
template<int dim> class TsDesc;

// included for Rom, Lei Lei, 27 Sep 2016
template<class VecType> class VecSet;

//------------------------------------------------------------------------------

template<int dim, int neq>
class MatVecProd {

public:

  MatVecProd() : isFSI(false) {}
  virtual ~MatVecProd() {}

  virtual void constructOperators(Vec3D &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &,
                                  double, DistSVec<double,dim> &, DistVec<double> &, DistTimeState<dim> *, PostOperator<dim> *) = 0;
  virtual dRdXoperators<dim> *getdRdXop() = 0;
  virtual void exportMemory(MemoryPool *mp)
  {
//    Dev::Error(MPI_COMM_WORLD, "exportMemory not implemented",true); sleep(2000);
//    std::cout<<"*** Error: exportMemory not implemented for MatVecProd:"<<__LINE__<<std::endl; exit(-1);
  }//TODO SUSPICIOUS

  virtual void evaluate(int, DistSVec<double,3> &, DistVec<double> &, 
			DistSVec<double,dim> &, DistSVec<double,dim> &) = 0;

  virtual void evaluate(DistExactRiemannSolver<dim> &, int, DistSVec<double,3> &, DistVec<double> &, 
			DistSVec<double,dim> &, DistSVec<double,dim> &) = 0;

  virtual void evaluateHH(DistVec<double> &hhterm,
			  DistVec<double> &bcVal )
  {std::cout<<"*** Error: evaluateHH not implemented for this MatVecProd"<<std::endl; exit(-1);}

  virtual void apply(DistSVec<double,neq> &, DistSVec<double,neq> &) = 0;
  virtual void apply(DistSVec<bcomp,neq> &, DistSVec<bcomp,neq> &) = 0;
  virtual void apply(DistVec<double> &, DistVec<double> &)
  {std::cout<<"*** Error: apply not implemented for this MatVecProd"<<std::endl; exit(-1);}

  virtual void apply(DistEmbeddedVec<double,neq> &, DistEmbeddedVec<double,neq> &)
  {std::cout<<"*** Error: apply not implemented for this MatVecProd"<<std::endl; exit(-1);}
  
  virtual void applyT(DistSVec<double,neq> &, DistSVec<double,neq> &) = 0;
  virtual void applyT(DistSVec<bcomp,neq> &, DistSVec<bcomp,neq> &) = 0;

  virtual void applyTranspose(DistSVec<double,neq> &, DistSVec<double,neq> &) = 0;
  virtual void applyTranspose(DistEmbeddedVec<double,neq> &, DistEmbeddedVec<double,neq> &)
  {std::cout<<"*** Error: applyTranspose not implemented for this MatVecProd"<<std::endl; exit(-1);}

  virtual void applyTranspose(DistSVec<bcomp,neq> &, DistSVec<bcomp,neq> &)
  {std::cout<<"*** Error: applyTranspose not implemented for this MatVecProd"<<std::endl; exit(-1);}

  virtual void evaluateRestrict(int, DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &, RestrictionMapping<dim> &,
                TsDesc<dim>* probDesc=NULL,
                int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &)=NULL) {
    std::cout<<"*** Error: function evaluateRestrict not implemented"<<std::endl;
    sleep(1);
    exit(-1);
  }
  virtual void applyRestrict(DistSVec<double,neq> &, DistSVec<double,neq> &, RestrictionMapping<neq> &,
                TsDesc<dim>* probDesc=NULL,
                int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &)=NULL) { 
    std::cout<<"*** Error: function applyRestrict not implemented"<<std::endl;
    sleep(1);
    exit(-1);
  }


  // ROMs minimize the residual, so the weighting of the residual becomes very important.
  // These functions allow for Jacobians of residuals with non-constant weights.
  virtual void evaluateWeighted(int, DistSVec<double,3> &, DistVec<double> &, 
			DistSVec<double,dim> &, DistSVec<double,dim> &, VarFcn *) {
    std::cout<<"*** Error: function evaluateWeighted not implemented"<<std::endl;
    sleep(1);
    exit(-1);
  }

  virtual void applyWeighted(DistSVec<double,neq> &, DistSVec<double,neq> &, VarFcn *){
    std::cout<<"*** Error: function applyWeighted not implemented"<<std::endl;
    sleep(1);
    exit(-1);
  }

  virtual void evaluateWeightedRestrict(int, DistSVec<double,3> &, DistVec<double> &, 
                DistSVec<double,dim> &, DistSVec<double,dim> &, RestrictionMapping<dim> &, VarFcn *,
                TsDesc<dim>* probDesc=NULL, int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &)=NULL) {
    std::cout<<"*** Error: function evaluateWeighted not implemented"<<std::endl;
    sleep(1);
    exit(-1);
  }

  virtual void applyWeightedRestrict(DistSVec<double,neq> &, DistSVec<double,neq> &, RestrictionMapping<dim> &, VarFcn *,
               TsDesc<dim>* probDesc=NULL, int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &)=NULL) {
    std::cout<<"*** Error: function applyWeighted not implemented"<<std::endl;
    sleep(1);
    exit(-1);
  }


// Included (MB)
  virtual void evaluateInviscid(int , DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &){
    std::cout<<"*** Error: function evaluateInviscid not implemented"<<std::endl; exit(-1);}
  virtual void evaluateViscous(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &){
    std::cout<<"*** Error: function evaluateViscous not implemented"<<std::endl;exit(-1);}
  virtual void applyInviscid(DistSVec<double,neq> &, DistSVec<double,neq> &){
    std::cout<<"*** Error: function applyInviscid not implemented"<<std::endl;exit(-1);}
  virtual void applyViscous(DistSVec<double,neq> &, DistSVec<double,neq> &){
    std::cout<<"*** Error: function applyViscous not implemented"<<std::endl;exit(-1);}
  virtual void rstSpaceOp(IoData &, VarFcn *, SpaceOperator<dim> *, bool, SpaceOperator<dim> * = 0){
    std::cout<<"*** Error: function rstSpaceOp not implemented"<<std::endl;exit(-1);}

  virtual void attachHH(DistEmbeddedVec<double,dim>& v)
  {std::cout<<"*** Error: function attachHH not implemented"<<std::endl;exit(-1);}

  // Structure to enable fluid-structure interaction computations
  struct _fsi {
    DistLevelSetStructure* LSS;
    DistVec<int>* fluidId;
    DistExactRiemannSolver<dim>* riemann;
    bool linRecAtInterface, viscSecOrder;
    DistSVec<double,dim>* Wtemp;
    int Nriemann;
    DistVec<GhostPoint<dim>*>* ghostPoints;
  };

  void AttachStructure(const _fsi& f) {
    isFSI = true;
    fsi = f;
  }

  virtual void setTsOutput(TsOutput<dim>* outputPointer)
  {std::cout<<"*** Error: function setTsOutput not implemented"<<std::endl; exit(-1);}

protected:
  
  // Boolean; set to true if we are using a structure
  bool isFSI;
  _fsi fsi;
};

//------------------------------------------------------------------------------

template<int dim, int neq>
class MatVecProdFD : public MatVecProd<dim,neq> {

  DistSVec<double,dim> Qeps;
  DistSVec<double,dim> Feps;
  DistSVec<double,neq> Qepstmp;
  DistSVec<double,neq> Fepstmp;

  SpaceOperator<dim> *spaceOp;
  RecFcnConstant<dim> *recFcnCon;
  DistSVec<double,dim> *Rn;

  DistTimeState<dim> *timeState;
  DistGeoState *geoState;
  DistSVec<double,3> *X;
  DistVec<double> *ctrlVol;
  DistSVec<double,neq> Q;
  DistSVec<double,neq> F;

  Communicator *com;

  DistSVec<double, neq>* energyWeightVec;

// Included (MB)
  DistSVec<double,neq> Ftmp;
  IoData *iod;
  int fdOrder;
  double fdeps;

  DistVec<double>* hhRes,*hhEps,*hhVal;
  
  // for outputting residuals computed during FD
  TsOutput<dim>* output;
  void setTsOutput(TsOutput<dim>* outputPointer) {output = outputPointer;}
 
public:

  double computeEpsilon(DistSVec<double,neq> &, DistSVec<double,neq> &);

// Included (MB)
  MatVecProdFD(ImplicitData &, DistTimeState<dim> *, DistGeoState *, 
	       SpaceOperator<dim> *, Domain *, IoData &);
  
  ~MatVecProdFD();

  void attachHH(DistEmbeddedVec<double,dim>& v);

  void constructOperators(Vec3D &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &,
                          double, DistSVec<double,dim> &, DistVec<double> &, DistTimeState<dim> *, PostOperator<dim> *)
  {std::cout<<"*** Error: function constructOperators not implemented"<<std::endl; exit(-1);}

  dRdXoperators<dim> *getdRdXop() { return 0;}
  void evaluate(int, DistSVec<double,3> &, DistVec<double> &, 
		DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluateWeighted(int, DistSVec<double,3> &, DistVec<double> &, 
		DistSVec<double,dim> &, DistSVec<double,dim> &, VarFcn *);
  void evaluate(DistExactRiemannSolver<dim> &, int, DistSVec<double,3> &, DistVec<double> &, 
		DistSVec<double,dim> &, DistSVec<double,dim> &);

  void evaluateHH(DistVec<double> &hhterm,
		  DistVec<double> &bcVal );

  void evaluateRestrict(int, DistSVec<double,3> &, DistVec<double> &, 
		DistSVec<double,dim> &, DistSVec<double,dim> &, RestrictionMapping<dim> &,
                TsDesc<dim>* probDesc=NULL,
                int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &)=NULL);
  void evaluateWeightedRestrict(int, DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &, RestrictionMapping<dim> &, VarFcn *,
                TsDesc<dim>* probDesc=NULL,
                int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &)=NULL);

  void apply(DistSVec<double,neq> &, DistSVec<double,neq> &);
  void applyTranspose(DistSVec<double,neq> &, DistSVec<double,neq> &)
  {std::cout<<"*** Error: applyTranspose not implemented for this MatVecProdFD"<<std::endl; exit(-1);}

  void applyTranspose(DistEmbeddedVec<double,neq> &, DistEmbeddedVec<double,neq> &)
  {std::cout<<"*** Error: applyTranspose not implemented for this MatVecProdFD"<<std::endl; exit(-1);}

  void applyWeighted(DistSVec<double,neq> &, DistSVec<double,neq> &, VarFcn *);
  void applyRestrict(DistSVec<double,neq> &, DistSVec<double,neq> &, RestrictionMapping<neq> &,
                TsDesc<dim>* probDesc=NULL,
                int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &)=NULL);
  void applyWeightedRestrict(DistSVec<double,neq> &, DistSVec<double,neq> &, RestrictionMapping<neq> &, VarFcn*,
                TsDesc<dim>* probDesc=NULL,
                int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &)=NULL);

  void apply(DistEmbeddedVec<double,neq> &, DistEmbeddedVec<double,neq> &); //!!!

  void apply(DistSVec<bcomp,neq> &, DistSVec<bcomp,neq> &)  {
    std::cout << "... ERROR: ::apply function not implemented for class MatVecProdFD with complex arguments" << endl; }

  void applyT(DistSVec<double,neq> &, DistSVec<double,neq> &)  {
    std::cout << "... ERROR: ::applyT function not implemented for class MatVecProdFD" << endl; }
  void applyT(DistSVec<bcomp,neq> &x, DistSVec<bcomp,neq> &y)  {
    std::cout << "... ERROR: ::applyT function not implemented for class MatVecProdFD" << endl; }

// Included (MB)
  void evaluateInviscid(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluateViscous(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void applyInviscid(DistSVec<double,neq> &, DistSVec<double,neq> &);

  void applyViscous(DistSVec<double,neq> &, DistSVec<double,neq> &);

  void applyViscous(DistSVec<bcomp,neq> &, DistSVec<bcomp,neq> &)
  {
    std::cout << " ... ERROR: MatVecProdFD::applyViscous function not implemented";
    std::cout << " for complex arguments." << std::endl;
    exit(1);
  }

  void rstSpaceOp(IoData &, VarFcn *, SpaceOperator<dim> *, bool, SpaceOperator<dim> * = 0);


  
};

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
class MatVecProdH1 : public MatVecProd<dim,neq>, public DistMat<Scalar,neq> {

#ifdef _OPENMP
  int numLocSub; //BUG omp
#endif

  MvpMat<Scalar,neq> **A;

  DistTimeState<dim> *timeState;
  SpaceOperator<dim> *spaceOp;

  bool areHHTermsActive;

  DistVec<double>* hhVal;

public:

  /// Constructor.
  /// \note UH (09/10) This constructor is only called from MatVecProdH2.
  MatVecProdH1(DistTimeState<dim> *, SpaceOperator<dim> *, Domain *);

  /// Constructor.
  MatVecProdH1(DistTimeState<dim> *, SpaceOperator<dim> *, Domain *, IoData &);

  /// Destructor.
  ~MatVecProdH1();

  DistMat<Scalar,neq> &operator= (const Scalar);

  void constructOperators(Vec3D &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &,
                          double, DistSVec<double,dim> &, DistVec<double> &, DistTimeState<dim> *, PostOperator<dim> *)
  {std::cout<<"*** Error: function constructOperators not implemented"<<std::endl;exit(-1);}

  dRdXoperators<dim> *getdRdXop() { return 0; }

  GenMat<Scalar,neq> &operator() (int i) { return *A[i]; }

  void attachHH(DistEmbeddedVec<double,dim>& v);
  void evaluateHH(DistVec<double> &hhterm,
		  DistVec<double> &bcVal );

  void exportMemory(MemoryPool *);

  void evaluate(int, DistSVec<double,3> &, DistVec<double> &, 
		DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluate(DistExactRiemannSolver<dim>&, int, DistSVec<double,3> &, DistVec<double> &, 
		DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluateViscous(int, DistSVec<double,3> &, DistVec<double> &);

  void apply(DistSVec<double,neq> &, DistSVec<double,neq> &);
  void applyTranspose(DistSVec<double,neq> &, DistSVec<double,neq> &);
  void applyTranspose(DistEmbeddedVec<double,neq> &, DistEmbeddedVec<double,neq> &)
  {std::cout<<"*** Error: function applyTranspose not implemented"<<std::endl;exit(-1);}

  void apply(DistSVec<bcomp,neq> &, DistSVec<bcomp,neq> &)  {
    std::cout << "... ERROR: ::apply function not implemented for class MatVecProdH1 with complex arguments" << endl; exit(-1);}

  void evaluateRestrict(int, DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &, RestrictionMapping<dim> &,
                TsDesc<dim>* probDesc=NULL,
                int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &)=NULL);
  void applyRestrict(DistSVec<double,neq> &, DistSVec<double,neq> &, RestrictionMapping<neq> &,
                TsDesc<dim>* probDesc=NULL,
                int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &)=NULL);

  void apply(DistEmbeddedVec<double,neq> &, DistEmbeddedVec<double,neq> &);

  void applyT(DistSVec<double,neq> &, DistSVec<double,neq> &)  {
    std::cout << "... ERROR: ::applyT function not implemented for class MatVecProdH1" << endl; exit(-1);}
  void applyT(DistSVec<bcomp,neq> &x, DistSVec<bcomp,neq> &y) { 
    std::cout << "... ERROR: ::applyT function not implemented for class MatVecProdH1" << endl; exit(-1);}

  void rstSpaceOp(IoData &, VarFcn *, SpaceOperator<dim> *, bool, SpaceOperator<dim> * = 0);
 
  void clearGhost();
};

//------------------------------------------------------------------------------

///
/// The MatVecProdH2 class enables matrix vector product based on the exact
/// Jacobian matrix.
/// For turbulent problems with weak coupling, the parameters dim and neq will
/// be different.
///
/// \note (09/10) The implementation for weak turbulence-model coupling is
/// not optimal because it uses a vector of local length 'dim' for the
/// inviscid part.
/// For the viscous part (laminar + turbulent), the product is using
/// an object MatVecProdH1.
///
template<int dim, class Scalar, int neq>
class MatVecProdH2 : public MatVecProd<dim,neq>, public DistMat<Scalar,dim> {

#ifdef _OPENMP
  int numLocSub; //BUG omp
#endif

  // format of A is the diagonal terms and then the edge contributions
  MvpMat<Scalar,dim> **A;  

  // coefficients in the linearization of the reconstructed-limited primitive states
  /* see Lesoinne et al. A linearized method for the frequency analysis of three 
     dimensional fluid/structure interaction problems in all flow regimes, Comp. Meth. Appl. Mech. Eng.
     vol. 190 (2001) pp 3121-3146 */
  DistSVec<double,dim> aij; 
  DistSVec<double,dim> aji;
  DistSVec<double,dim> bij;
  DistSVec<double,dim> bji;
  DistSVec<double,dim> betaij;
  DistSVec<double,dim> betaji;

  DistTimeState<dim> *timeState;
  SpaceOperator<dim> *spaceOp;
  FluxFcn **fluxFcn;

  DistSVec<double,3> *X;
  DistVec<double> *ctrlVol;
  DistSVec<double,dim> *Q;
  DistSVec<double,dim> *F;

  // viscous flux jacobian and/or 1st order jac terms(ie face flux jac)
  MatVecProdH1<dim, Scalar, neq> *R;

  // Included (MB)
  MatVecProdFD<dim, neq> *RFD;
  DistSVec<double, neq> *vProd;

  //--------------------------------------
  /// \note (09/10) UH
  /// This nested class allows to specialize the matrix-vector product.
  /// In particular, we differentiate the cases where dim and neq are different
  /// for turbulent problems with weak coupling.
  template<int dd, int nn, class Scalar1, class Scalar2> struct Multiplier
  {
    void Apply
    (
      SpaceOperator<dd> *spaceOp
      , DistSVec<double,3> &X
      , DistVec<double> &ctrlVol
      , DistSVec<double,dd> &U
      , DistMat<Scalar1,dd> &H2
      , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
      , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
      , DistSVec<Scalar2,nn> &p, DistSVec<Scalar2,nn> &prod
      , MatVecProdH1<dd, Scalar1, nn> *R
      , MatVecProdFD<dd, nn> *RFD
      , DistSVec<Scalar2, nn> *vProd
    );
    void ApplyTranspose
    (
      SpaceOperator<dd> *spaceOp
      , DistSVec<double,3> &X
      , DistVec<double> &ctrlVol
      , DistSVec<double,dd> &U
      , DistMat<Scalar1,dd> &H2
      , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
      , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
      , DistSVec<Scalar2,nn> &p, DistSVec<Scalar2,nn> &prod
      , MatVecProdH1<dd, Scalar1, nn> *R
      , MatVecProdFD<dd, nn> *RFD
      , DistSVec<Scalar2, nn> *vProd
    );
    void ApplyT
    (
      SpaceOperator<dd> *spaceOp
      , DistSVec<double,3> &X
      , DistVec<double> &ctrlVol
      , DistSVec<double,dd> &U
      , DistMat<Scalar1,dd> &H2
      , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
      , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
      , DistSVec<Scalar2,nn> &p, DistSVec<Scalar2,nn> &prod
    );
  };

  template<int dd, class Scalar1, class Scalar2>
  struct Multiplier<dd,dd,Scalar1,Scalar2>
  {
    void Apply
    (
      SpaceOperator<dd> *spaceOp
      , DistSVec<double,3> &X
      , DistVec<double> &ctrlVol
      , DistSVec<double,dd> &U
      , DistMat<Scalar1,dd> &H2
      , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
      , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
      , DistSVec<Scalar2,dd> &p, DistSVec<Scalar2,dd> &prod
      , MatVecProdH1<dd, Scalar1, dd> *R
      , MatVecProdFD<dd, dd> *RFD
      , DistSVec<Scalar2, dd> *vProd
    );
    void ApplyTranspose
    (
      SpaceOperator<dd> *spaceOp
      , DistSVec<double,3> &X
      , DistVec<double> &ctrlVol
      , DistSVec<double,dd> &U
      , DistMat<Scalar1,dd> &H2
      , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
      , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
      , DistSVec<Scalar2,dd> &p, DistSVec<Scalar2,dd> &prod
      , MatVecProdH1<dd, Scalar1, dd> *R
      , MatVecProdFD<dd, dd> *RFD
      , DistSVec<Scalar2, dd> *vProd
    );
    void ApplyT
    (
      SpaceOperator<dd> *spaceOp
      , DistSVec<double,3> &X
      , DistVec<double> &ctrlVol
      , DistSVec<double,dd> &U
      , DistMat<Scalar1,dd> &H2
      , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
      , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
      , DistSVec<Scalar2,dd> &p, DistSVec<Scalar2,dd> &prod
    );
  };
  //--------------------------------------


public:

// Included (MB)
  MatVecProdH2(IoData &, VarFcn *, DistTimeState<dim> *, 
	       SpaceOperator<dim> *, Domain *, DistGeoState * = 0);

  ~MatVecProdH2();

  DistMat<Scalar,dim> &operator= (const Scalar);

  GenMat<Scalar,dim> &operator() (int i) { return *A[i]; }

  void constructOperators(Vec3D &, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &,
                          double, DistSVec<double,dim> &, DistVec<double> &, DistTimeState<dim> *, PostOperator<dim> *)
  {std::cout<<"*** Error: function constructOperators not implemented"<<std::endl; exit(-1);}

  dRdXoperators<dim> *getdRdXop() { return 0; }
/*
  void evaluate(int, DistSVec<double,3> &, DistVec<double> &, 
		DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluate(int , DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &, Scalar);
  void evaluate2(int, DistSVec<double,3> &, DistVec<double> &, 
		 DistSVec<double,dim> &, DistSVec<double,dim> &); 
*/

  void evaluate(int, DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluate(DistExactRiemannSolver<dim> &, int, DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &)
  {std::cout<<"*** Error: function constructOperators not implemented"<<std::endl; exit(-1);}

  void evaluate(int , DistSVec<double,3> &, DistVec<double> &, 
                DistSVec<double,dim> &, DistSVec<double,dim> &, Scalar);

  // UH (09/10)
  // The following function is never called and not implemented.
  //void evaluate(int , DistSVec<double,3> &, DistVec<double> &, DistVec<int> &,
  //              DistSVec<double,dim> &, DistSVec<double,dim> &, Scalar);

  void evaluate2(int, DistSVec<double,3> &, DistVec<double> &,
                 DistSVec<double,dim> &, DistSVec<double,dim> &);

  // UH (09/10)
  // The following functions are never called and not implemented.
  //
  //void evaluatestep1(int, DistSVec<double,3> &, DistVec<double> &,
  //               DistSVec<double,dim> &, DistSVec<double,dim> &);
  //void evaluate2step1(int, DistSVec<double,3> &, DistVec<double> &,
  //               DistSVec<double,dim> &, DistSVec<double,dim> &);

  void evalH(int , DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &);

  void evaluateRestrict(int, DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &, RestrictionMapping<dim> &,
                TsDesc<dim>* probDesc=NULL,
                int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &)=NULL);
  void applyRestrict(DistSVec<double,neq> &, DistSVec<double,neq> &, RestrictionMapping<neq> &,
                TsDesc<dim>* probDesc=NULL,
                int (TsDesc<dim>::*checkSolution)(DistSVec<double,dim> &)=NULL);

  void apply(DistSVec<double,neq>        &, DistSVec<double,neq> &);
  void apply(DistSVec<bcomp,neq>         &, DistSVec<bcomp,neq> &);
  void apply(DistEmbeddedVec<double,dim> &, DistEmbeddedVec<double,dim> &);

  void applyT(DistSVec<double,neq> &, DistSVec<double,neq> &);
  void applyT(DistSVec<bcomp,neq> &x, DistSVec<bcomp,neq> &y);
  void applyTranspose(DistSVec<double,neq> &, DistSVec<double,neq> &);
  void applyTranspose(DistEmbeddedVec<double,neq> &, DistEmbeddedVec<double,neq> &)
  {std::cout<<"*** Error: function MatVecProdH2::applyTranspose not implemented for Embedded"<<std::endl; exit(-1);}

// Included (MB)
  void evaluateInviscid(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);
  void evaluateViscous(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &);

  //void applyInviscid(DistSVec<double,neq> &, DistSVec<double,neq> &);
  //void applyViscous(DistSVec<double,neq> &, DistSVec<double,neq> &);

  void rstSpaceOp(IoData &, VarFcn *, SpaceOperator<dim> *, bool, SpaceOperator<dim> * = 0);

};

//----------------------------------------------------------------------------//
//                MatVecProd for Multiphase Euler equations                   //
//----------------------------------------------------------------------------//
// this class is not derived from the MatVecProd class!
// Its templates are different and additional members are different.
template<int dim, int dimLS>
class MatVecProdMultiPhase {

protected:
  DistTimeState<dim> *timeState;
  MultiPhaseSpaceOperator<dim,dimLS> *spaceOp;
  DistExactRiemannSolver<dim> *riemann;
  FluidSelector *fluidSelector;

public:

  MatVecProdMultiPhase(DistTimeState<dim> *ts, MultiPhaseSpaceOperator<dim,dimLS> *spo,
                       DistExactRiemannSolver<dim> *rsolver, FluidSelector *fs) : 
                       timeState(ts), spaceOp(spo), riemann(rsolver), fluidSelector(fs), isFSI(false) {}//TODO SUSPICIOUS
  virtual ~MatVecProdMultiPhase() { timeState=0; spaceOp = 0; riemann = 0; fluidSelector = 0; }

  virtual void exportMemory(MemoryPool *mp)
  {std::cout<<"*** Error: exportMemory not implemented"<<std::endl; exit(-1);}

  virtual void evaluate(int, DistSVec<double,3> &, DistVec<double> &, 
                        DistSVec<double,dim> &, DistSVec<double,dimLS> &,
                        DistSVec<double,dim> &) = 0;

  virtual void apply(DistSVec<double,dim> &, DistSVec<double,dim> &) = 0;

  virtual void apply(DistEmbeddedVec<double,dim> &, DistEmbeddedVec<double,dim> &)
  {std::cout<<"*** Error: apply not implemented for this MatVecProd:"<<__LINE__<<std::endl; exit(-1);}

  virtual void applyT(DistSVec<double,dim> &, DistSVec<double,dim> &) = 0;
  virtual void applyT(DistSVec<bcomp,dim> &x, DistSVec<bcomp,dim> &y) = 0;
  virtual void applyTranspose(DistSVec<double,dim> &, DistSVec<double,dim> &) = 0;
  
  virtual void applyTranspose(DistEmbeddedVec<double,dim> &, DistEmbeddedVec<double,dim> &)
  {std::cout<<"*** Error: applyTranspose not implemented for MultiPhase"<<std::endl; exit(-1);}

  virtual void evaluateHH(DistVec<double> &hhterm,
			  DistVec<double> &bcVal )
  {std::cout<<"*** Error: evaluateHH not implemented for MultiPhase"<<std::endl; exit(-1);}

  virtual void attachHH(DistEmbeddedVec<double,dim>& v)
  {std::cout<<"*** Error: attachHH not implemented for MultiPhase"<<std::endl; exit(-1);}
  
  // Structure to enable fluid-structure interaction computations
  struct _fsi {

    DistLevelSetStructure* LSS;
    DistVec<int>* fluidId;
    DistExactRiemannSolver<dim>* riemann;
    bool linRecAtInterface, viscSecOrder;
    DistSVec<double,dim>* Wtemp;
    int Nriemann;
    DistVec<GhostPoint<dim>*>* ghostPoints;
  };

  void AttachStructure(const _fsi& f) {
    isFSI = true;
    fsi = f;
  }
 
protected:
  
  // Boolean; set to true if we are using a structure
  bool isFSI;
  _fsi fsi;

};

//------------------------------------------------------------------------------
// Finite Difference method for Multiphase Euler equations
template<int dim, int dimLS>
class MatVecProdFDMultiPhase : public MatVecProdMultiPhase<dim,dimLS> {

  DistSVec<double,dim> Qeps;
  DistSVec<double,dim> Feps;
  DistSVec<double,dim> Q;
  DistSVec<double,dim> F;

  DistGeoState *geoState;
  DistSVec<double,3> *X;
  DistVec<double> *ctrlVol;
  DistSVec<double,dimLS> *Phi;

  Communicator *com;

  double computeEpsilon(DistSVec<double,dim> &, DistSVec<double,dim> &);

  IoData *iod;
  int fdOrder;

  DistVec<double>* hhRes,*hhEps,*hhVal;

public:

// Included (MB)
  MatVecProdFDMultiPhase(DistTimeState<dim> *, DistGeoState *, 
			 MultiPhaseSpaceOperator<dim,dimLS> *, DistExactRiemannSolver<dim> *,
			 FluidSelector *, Domain *,IoData& ioData);

  ~MatVecProdFDMultiPhase();

  void attachHH(DistEmbeddedVec<double,dim>& v);
  
  void evaluate(int, DistSVec<double,3> &, DistVec<double> &,
                     DistSVec<double,dim> &, DistSVec<double,dimLS> &,
                     DistSVec<double,dim> &);
  void apply(DistSVec<double,dim> &, DistSVec<double,dim> &);

  void apply(DistEmbeddedVec<double,dim> &, DistEmbeddedVec<double,dim> &);

  void applyT(DistSVec<double,dim> &, DistSVec<double,dim> &)
  {std::cout<<"*** Error: applyT not implemented for MultiPhase"<<std::endl; exit(-1);}

  void applyT(DistSVec<bcomp,dim> &x, DistSVec<bcomp,dim> &y)
  {std::cout<<"*** Error: applyT not implemented for MultiPhase"<<std::endl; exit(-1);}

  void applyTranspose(DistSVec<double,dim> &, DistSVec<double,dim> &)
  {std::cout<<"*** Error: applyTranspose not implemented for MultiPhase"<<std::endl; exit(-1);}

  void applyTranspose(DistEmbeddedVec<double,dim> &, DistEmbeddedVec<double,dim> &)
  {std::cout<<"*** Error: applyTranspose not implemented for MultiPhase"<<std::endl; exit(-1);}
  
  void evaluateHH(DistVec<double> &hhterm,
		  DistVec<double> &bcVal );
};

//------------------------------------------------------------------------------
// H1 MatVecProd for Multiphase Euler equations
template<int dim, int dimLS>
class MatVecProdH1MultiPhase : public MatVecProdMultiPhase<dim,dimLS>,
                               public DistMat<double,dim> {

  MvpMat<double,dim> **A;

  bool areHHTermsActive;

  DistVec<double>* hhVal;

public:

  MatVecProdH1MultiPhase(DistTimeState<dim> *, MultiPhaseSpaceOperator<dim,dimLS> *, 
                         DistExactRiemannSolver<dim> *, FluidSelector *, Domain *);
  ~MatVecProdH1MultiPhase();

  DistMat<double,dim> &operator= (const double);

  GenMat<double,dim> &operator() (int i) { return *A[i]; }

  void attachHH(DistEmbeddedVec<double,dim>& v);
  void evaluateHH(DistVec<double> &hhterm,
		  DistVec<double> &bcVal );
  
  void exportMemory(MemoryPool *);

  void evaluate(int, DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistSVec<double,dimLS> &,
                DistSVec<double,dim> &);

  void apply(DistSVec<double,dim> &, DistSVec<double,dim> &);

  void apply(DistEmbeddedVec<double,dim> &, DistEmbeddedVec<double,dim> &);
  void applyT(DistSVec<double,dim> &, DistSVec<double,dim> &)
  {std::cout<<"*** Error: applyT not implemented for MatVecProdH1MultiPhase"<<std::endl; exit(-1);}

  void applyT(DistSVec<bcomp,dim> &x, DistSVec<bcomp,dim> &y)
  {std::cout<<"*** Error: applyT not implemented for MatVecProdH1MultiPhase"<<std::endl; exit(-1);}

  void applyTranspose(DistSVec<double,dim> &, DistSVec<double,dim> &)
  {std::cout<<"*** Error: applyTranspose not implemented for MatVecProdH1MultiPhase"<<std::endl; exit(-1);}

  void applyTranspose(DistEmbeddedVec<double,dim> &, DistEmbeddedVec<double,dim> &)
  {std::cout<<"*** Error: applyTranspose not implemented for MatVecProdH1MultiPhase"<<std::endl; exit(-1);}
};

//----------------------------------------------------------------------------//
//                MatVecProd for Level Set equations                          //
//----------------------------------------------------------------------------//
// Finite Difference MatVecProd for LevelSet equation
template<int dim, int dimLS>
class MatVecProdLS : public DistMat<double,dimLS> {
                                                                                                                      
  DistTimeState<dim> *timeState;
  MultiPhaseSpaceOperator<dim,dimLS> *spaceOp;
  DistGeoState *geoState;
  LevelSet<dimLS> *levelSet;

  DistSVec<double,3> *X;
  DistVec<double> *ctrlVol;
  DistSVec<double,dim> *U;
  DistVec<int> *FluidId;
  DistSVec<double,dimLS> *Q;
  DistSVec<double,dim> *V;
  DistSVec<double,dimLS> *F;
  DistSVec<double,dimLS> Qeps;
  DistSVec<double,dimLS> Feps;

  Communicator *com;

  double computeEpsilon(DistSVec<double,dimLS> &, DistSVec<double,dimLS> &);

  MvpMat<double,dimLS> **A;
public:
  DistMat<double,dimLS> &operator= (const double);
  GenMat<double,dimLS> &operator() (int i) { return *A[i]; }

  MatVecProdLS(DistTimeState<dim> *, DistGeoState *,
               MultiPhaseSpaceOperator<dim,dimLS> *, Domain *, LevelSet<dimLS> *);
  ~MatVecProdLS();
                                                                                                                      
  void evaluate(int, DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dimLS> &, DistSVec<double,dim> &,
		DistSVec<double,dim> &,
                DistSVec<double,dimLS> &, DistVec<int> &,bool = false,DistLevelSetStructure* = NULL,
		int lsMethod = 0);

  void apply(DistSVec<double,dimLS> &, DistSVec<double,dimLS> &);

  void applyT(DistSVec<double,dimLS> &, DistSVec<double,dimLS> &)
  {std::cout<<"*** Error: function applyT not implemented for MatVecProdLS:"<<__LINE__<<std::endl; exit(-1);}

  void applyT(DistSVec<bcomp,dimLS> &x, DistSVec<bcomp,dimLS> &y)
  {std::cout<<"*** Error: function applyT not implemented for MatVecProdLS:"<<__LINE__<<std::endl; exit(-1);}

  void applyTranspose(DistSVec<double,dimLS> &, DistSVec<double,dimLS> &)
  {std::cout<<"*** Error: function applyTranspose not implemented for MatVecProdLS:"<<__LINE__<<std::endl; exit(-1);}

  void applyTranspose(DistEmbeddedVec<double,dimLS> &, DistEmbeddedVec<double,dimLS> &)
  {std::cout<<"*** Error: applyTranspose not implemented for MatVecProdLS:"<<__LINE__<<std::endl; exit(-1);}

};
                                                                                                                      
//------------------------------------------------------------------------------

template<int dim>
struct dRdXoperators
{
  int numLocSub;
  RectangularSparseMat<double,3,3> **dMidGradP;
  RectangularSparseMat<double,3,3> **dMidX;
  RectangularSparseMat<double,3,3> **dMvdX;
  RectangularSparseMat<double,dim,3> **dMvdV;
  RectangularSparseMat<double,dim,3> **dMidV;
  RectangularSparseMat<double,3,3> **dMidS;
  RectangularSparseMat<double,3,3> **dFidGradP;
  RectangularSparseMat<double,3,3> **dFidX;
  RectangularSparseMat<double,3,3> **dFvdX;
  RectangularSparseMat<double,dim,3> **dFidV;
  RectangularSparseMat<double,dim,3> **dFvdV;
  RectangularSparseMat<double,3,3> **dFidS;
  RectangularSparseMat<double,3,1> **dCtrlVoldX;
  RectangularSparseMat<double,3,3> **dEdgeNormdX;
  RectangularSparseMat<double,3,3> **dFaceNormdX;
  RectangularSparseMat<double,3,6> **dRdX;
  RectangularSparseMat<double,6,6> **dRdR;
  RectangularSparseMat<double,3,dim> **dddxdX;
  RectangularSparseMat<double,3,dim> **dddydX;
  RectangularSparseMat<double,3,dim> **dddzdX;
  RectangularSparseMat<double,dim,dim> **dddxdV;
  RectangularSparseMat<double,dim,dim> **dddydV;
  RectangularSparseMat<double,dim,dim> **dddzdV;
  RectangularSparseMat<double,6,dim> **dddxdR;
  RectangularSparseMat<double,6,dim> **dddydR;
  RectangularSparseMat<double,6,dim> **dddzdR;
  RectangularSparseMat<double,dim,dim> **dFluxdddx;
  RectangularSparseMat<double,dim,dim> **dFluxdddy;
  RectangularSparseMat<double,dim,dim> **dFluxdddz;
  RectangularSparseMat<double,3,dim> **dFluxdEdgeNorm;
  RectangularSparseMat<double,3,dim> **dFluxdFaceNormal;
  RectangularSparseMat<double,1,dim> **dFluxdFaceNormalVel;
  RectangularSparseMat<double,dim,dim> **dFluxdUb;
  RectangularSparseMat<double,3,dim> **dFluxdX;
  RectangularSparseMat<double,3,dim> **dViscousFluxdX;
  RectangularSparseMat<double,dim,3> **dGradPdddx;
  RectangularSparseMat<double,dim,3> **dGradPdddy;
  RectangularSparseMat<double,dim,3> **dGradPdddz;
  RectangularSparseMat<double,dim,3> **dForcedV;
  RectangularSparseMat<double,3,3> **dForcedX;
  RectangularSparseMat<double,3,3> **dForcedGradP;
  RectangularSparseMat<double,3,3> **dForcedS;
  RectangularSparseMat<double,1,dim> **dVdPstiff;
  RectangularSparseMat<double,dim,dim> **dVdU;

  // Constructor
  dRdXoperators() { dRdX = 0; dRdR = 0; dddxdX = 0; dddydX = 0; dddzdX = 0; dddxdR = 0; 
                    dddydR = 0; dddzdR = 0; dEdgeNormdX = 0; dFaceNormdX = 0; dCtrlVoldX = 0; dViscousFluxdX = 0;
                    dFidS = 0;  dFidX = 0;  dFvdX = 0;  dFidV = 0;  dFvdV = 0;     dFidGradP = 0;
                    dMidS = 0;  dMidX = 0;  dMidV = 0;       dMidGradP = 0;  dMvdX = 0;  dMvdV = 0;
                    dFluxdddx = 0;  dFluxdddy = 0;  dFluxdddz = 0;  dFluxdEdgeNorm = 0;  dFluxdX = 0; 
                    dFluxdFaceNormal = 0;  dFluxdFaceNormalVel = 0;  dFluxdUb = 0;
                    dGradPdddx = 0;  dGradPdddy = 0;  dGradPdddz = 0; 
                    dForcedV = 0;    dForcedX = 0;   dForcedGradP = 0;   dForcedS = 0; 
                    dVdPstiff = 0;   dVdU = 0;    
                    dddxdV = 0;  dddydV = 0;   dddzdV = 0;  }
  ~dRdXoperators() {

    if (dVdPstiff) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dVdPstiff[iSub]) delete dVdPstiff[iSub];

      delete [] dVdPstiff;
    }

    if (dVdU) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dVdU[iSub]) delete dVdU[iSub];

      delete [] dVdU;
    }

    if (dForcedS) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dForcedS[iSub]) delete dForcedS[iSub];

      delete [] dForcedS;
    }

    if (dForcedGradP) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dForcedGradP[iSub]) delete dForcedGradP[iSub];

      delete [] dForcedGradP;
    }

    if (dForcedX) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dForcedX[iSub]) delete dForcedX[iSub];

      delete [] dForcedX;
    }

    if (dForcedV) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dForcedV[iSub]) delete dForcedV[iSub];

      delete [] dForcedV;
    }

    if (dGradPdddx) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dGradPdddx[iSub]) delete dGradPdddx[iSub];

      delete [] dGradPdddx;
    }

    if (dGradPdddy) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dGradPdddy[iSub]) delete dGradPdddy[iSub];

      delete [] dGradPdddy;
    }

    if (dGradPdddz) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dGradPdddz[iSub]) delete dGradPdddz[iSub];

      delete [] dGradPdddz;
    }

    if (dRdX) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dRdX[iSub]) delete dRdX[iSub];

      delete [] dRdX;
    }

    if (dRdR) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dRdR[iSub]) delete dRdR[iSub];

      delete [] dRdR;
    }

    if (dddxdV) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dddxdV[iSub]) delete dddxdV[iSub];

      delete [] dddxdV;
    }

    if (dddydV) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dddydV[iSub]) delete dddydV[iSub];

      delete [] dddydV;
    }

    if (dddzdV) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dddzdV[iSub]) delete dddzdV[iSub];

      delete [] dddzdV;
    }

    if (dddxdX) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dddxdX[iSub]) delete dddxdX[iSub];

      delete [] dddxdX;
    }

    if (dddydX) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dddydX[iSub]) delete dddydX[iSub];

      delete [] dddydX;
    }

    if (dddzdX) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dddzdX[iSub]) delete dddzdX[iSub];

      delete [] dddzdX;
    }

    if (dddxdR) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dddxdR[iSub]) delete dddxdR[iSub];

      delete [] dddxdR;
    }

    if (dddydR) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dddydR[iSub]) delete dddydR[iSub];

      delete [] dddydR;
    }

    if (dddzdR) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dddzdR[iSub]) delete dddzdR[iSub];

      delete [] dddzdR;
    }

    if (dFaceNormdX) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dFaceNormdX[iSub]) delete dFaceNormdX[iSub];

      delete [] dFaceNormdX;
    }

    if (dEdgeNormdX) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dEdgeNormdX[iSub]) delete dEdgeNormdX[iSub];

      delete [] dEdgeNormdX;
    }

    if (dCtrlVoldX) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dCtrlVoldX[iSub]) delete dCtrlVoldX[iSub];

      delete [] dCtrlVoldX;
    }
    if (dViscousFluxdX) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub)
        if (dViscousFluxdX[iSub]) delete dViscousFluxdX[iSub];

      delete [] dViscousFluxdX;
    }

    if (dMidX) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub)
        if (dMidX[iSub]) delete dMidX[iSub];

      delete [] dMidX;
    }

    if (dMvdX) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub)
        if (dMvdX[iSub]) delete dMvdX[iSub];

      delete [] dMvdX;
    }

    if (dMvdV) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub)
        if (dMvdV[iSub]) delete dMvdV[iSub];

      delete [] dMvdV;
    }

    if (dMidV) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub)
        if (dMidV[iSub]) delete dMidV[iSub];

      delete [] dMidV;
    }

    if (dMidS) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub)
        if (dMidS[iSub]) delete dMidS[iSub];

      delete [] dMidS;
    }

    if (dMidGradP) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub)
        if (dMidGradP[iSub]) delete dMidGradP[iSub];

      delete [] dMidGradP;
    }


    if (dFvdX) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub)
        if (dFvdX[iSub]) delete dFvdX[iSub];

      delete [] dFvdX;
    }

    if (dFidX) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub)
        if (dFidX[iSub]) delete dFidX[iSub];

      delete [] dFidX;
    }

    if (dFidV) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub)
        if (dFidV[iSub]) delete dFidV[iSub];

      delete [] dFidV;
    }

    if (dFvdV) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub)
        if (dFvdV[iSub]) delete dFvdV[iSub];

      delete [] dFvdV;
    }

    if (dFidS) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub)
        if (dFidS[iSub]) delete dFidS[iSub];

      delete [] dFidS;
    }

    if (dFidGradP) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub)
        if (dFidGradP[iSub]) delete dFidGradP[iSub];

      delete [] dFidGradP;
    }

    if (dFluxdddx) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dFluxdddx[iSub]) delete dFluxdddx[iSub];

      delete [] dFluxdddx;
    }

    if (dFluxdddy) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dFluxdddy[iSub]) delete dFluxdddy[iSub];

      delete [] dFluxdddy;
    }

    if (dFluxdddz) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dFluxdddz[iSub]) delete dFluxdddz[iSub];

      delete [] dFluxdddz;
    }

    if (dFluxdX) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dFluxdX[iSub]) delete dFluxdX[iSub];

      delete [] dFluxdX;
    }

    if (dFluxdEdgeNorm) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dFluxdEdgeNorm[iSub]) delete dFluxdEdgeNorm[iSub];

      delete [] dFluxdEdgeNorm;
    }

    if (dFluxdFaceNormal) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dFluxdFaceNormal[iSub]) delete dFluxdFaceNormal[iSub];

      delete [] dFluxdFaceNormal;
    }

    if (dFluxdFaceNormalVel) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dFluxdFaceNormalVel[iSub]) delete dFluxdFaceNormalVel[iSub];

      delete [] dFluxdFaceNormalVel;
    }

    if (dFluxdUb) {
#pragma omp parallel for
      for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
        if (dFluxdUb[iSub]) delete dFluxdUb[iSub];

      delete [] dFluxdUb;
    }

  }


  

  void setNumLocSub(int nLocSub) { numLocSub = nLocSub; }

  void initialize() {

    if(numLocSub == 0) {   fprintf(stderr," ... Error: need to set numLocSub first in dRdXoperators\n");   exit(-1);   }
    dFaceNormdX = new RectangularSparseMat<double,3,3>*[numLocSub];
    dEdgeNormdX = new RectangularSparseMat<double,3,3>*[numLocSub];
    dCtrlVoldX = new RectangularSparseMat<double,3,1>*[numLocSub];
    dMidX = new RectangularSparseMat<double,3,3>*[numLocSub];
    dMvdX = new RectangularSparseMat<double,3,3>*[numLocSub];
    dMvdV = new RectangularSparseMat<double,dim,3>*[numLocSub];
    dMidV = new RectangularSparseMat<double,dim,3>*[numLocSub];
    dMidS = new RectangularSparseMat<double,3,3>*[numLocSub];
    dMidGradP = new RectangularSparseMat<double,3,3>*[numLocSub];
    dFidX = new RectangularSparseMat<double,3,3>*[numLocSub];
    dFvdX = new RectangularSparseMat<double,3,3>*[numLocSub];
    dFidV = new RectangularSparseMat<double,dim,3>*[numLocSub];
    dFvdV = new RectangularSparseMat<double,dim,3>*[numLocSub];
    dFidS = new RectangularSparseMat<double,3,3>*[numLocSub];
    dFidGradP = new RectangularSparseMat<double,3,3>*[numLocSub];
    dRdX = new RectangularSparseMat<double,3,6>*[numLocSub];
    dRdR = new RectangularSparseMat<double,6,6>*[numLocSub];
    dddxdX = new RectangularSparseMat<double,3,dim>*[numLocSub];
    dddydX = new RectangularSparseMat<double,3,dim>*[numLocSub];
    dddzdX = new RectangularSparseMat<double,3,dim>*[numLocSub];
    dddxdV = new RectangularSparseMat<double,dim,dim>*[numLocSub];
    dddydV = new RectangularSparseMat<double,dim,dim>*[numLocSub];
    dddzdV = new RectangularSparseMat<double,dim,dim>*[numLocSub];
    dddxdR = new RectangularSparseMat<double,6,dim>*[numLocSub];
    dddydR = new RectangularSparseMat<double,6,dim>*[numLocSub];
    dddzdR = new RectangularSparseMat<double,6,dim>*[numLocSub];
    dFluxdddx = new RectangularSparseMat<double,dim,dim>*[numLocSub];
    dFluxdddy = new RectangularSparseMat<double,dim,dim>*[numLocSub];
    dFluxdddz = new RectangularSparseMat<double,dim,dim>*[numLocSub];
    dFluxdEdgeNorm = new RectangularSparseMat<double,3,dim>*[numLocSub];
    dFluxdFaceNormal = new RectangularSparseMat<double,3,dim>*[numLocSub];
    dFluxdFaceNormalVel = new RectangularSparseMat<double,1,dim>*[numLocSub];
    dFluxdUb = new RectangularSparseMat<double,dim,dim>*[numLocSub];
    dFluxdX = new RectangularSparseMat<double,3,dim>*[numLocSub];
    dViscousFluxdX = new RectangularSparseMat<double,3,dim>*[numLocSub];//YC
    dGradPdddx = new RectangularSparseMat<double,dim,3>*[numLocSub];
    dGradPdddy = new RectangularSparseMat<double,dim,3>*[numLocSub];
    dGradPdddz = new RectangularSparseMat<double,dim,3>*[numLocSub];
    dForcedGradP = new RectangularSparseMat<double,3,3>*[numLocSub];
    dForcedX = new RectangularSparseMat<double,3,3>*[numLocSub];
    dForcedV = new RectangularSparseMat<double,dim,3>*[numLocSub];
    dForcedS = new RectangularSparseMat<double,3,3>*[numLocSub];
    dVdU = new RectangularSparseMat<double,dim,dim>*[numLocSub];
    dVdPstiff = new RectangularSparseMat<double,1,dim>*[numLocSub];
  }

};

//------------------------------------------------------------------------------

template<int dim, class Scalar, int neq>
class MatVecProd_dRdX : public MatVecProd<dim,neq> {

  int numLocSub;

  dRdXoperators<dim> *dRdXop;

  DistTimeState<dim> *timeState;
  SpaceOperator<dim> *spaceOp;
  FluxFcn **fluxFcn;
  SubDomain **subDomain;
  Communicator *com;
  IoData* iod;

  DistSVec<double,3> *X;
  DistVec<double> *ctrlVol;
  DistSVec<double,dim> *Q;
  DistSVec<double,dim> *F;

  //--------------------------------------
  /// \note (09/10) UH
  /// This nested class allows to specialize the matrix-vector product.
  /// In particular, we differentiate the cases where dim and neq are different
  /// for turbulent problems with weak coupling.
  template<int dd, int nn, class Scalar1, class Scalar2> struct Multiplier
  {
    void Apply
    (
      SpaceOperator<dd> *spaceOp
      , DistSVec<double,3> &X
      , DistVec<double> &ctrlVol
      , DistSVec<double,dd> &U
      , DistMat<Scalar1,dd> &H2
      , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
      , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
      , DistSVec<Scalar2,nn> &p, DistSVec<Scalar2,nn> &prod
      , MatVecProdH1<dd, Scalar1, nn> *R
      , MatVecProdFD<dd, nn> *RFD
      , DistSVec<Scalar2, nn> *vProd
    );
    void ApplyT
    (
      SpaceOperator<dd> *spaceOp
      , DistSVec<double,3> &X
      , DistVec<double> &ctrlVol
      , DistSVec<double,dd> &U
      , DistMat<Scalar1,dd> &H2
      , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
      , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
      , DistSVec<Scalar2,nn> &p, DistSVec<Scalar2,nn> &prod
    );
  };

  template<int dd, class Scalar1, class Scalar2>
  struct Multiplier<dd,dd,Scalar1,Scalar2>
  {
    void Apply
    (
      SpaceOperator<dd> *spaceOp
      , DistSVec<double,3> &X
      , DistVec<double> &ctrlVol
      , DistSVec<double,dd> &U
      , DistMat<Scalar1,dd> &H2
      , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
      , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
      , DistSVec<Scalar2,dd> &p, DistSVec<Scalar2,dd> &prod
      , MatVecProdH1<dd, Scalar1, dd> *R
      , MatVecProdFD<dd, dd> *RFD
      , DistSVec<Scalar2, dd> *vProd
    );

    void ApplyT
    (
      SpaceOperator<dd> *spaceOp
      , DistSVec<double,3> &X
      , DistVec<double> &ctrlVol
      , DistSVec<double,dd> &U
      , DistMat<Scalar1,dd> &H2
      , DistSVec<double,dd> &aij, DistSVec<double,dd> &aji
      , DistSVec<double,dd> &bij, DistSVec<double,dd> &bji
      , DistSVec<Scalar2,dd> &p, DistSVec<Scalar2,dd> &prod
    );
  };
  //--------------------------------------


public:

// Included (MB)
  MatVecProd_dRdX(IoData &, VarFcn *, DistTimeState<dim> *, 
	       SpaceOperator<dim> *, Domain *, DistGeoState * = 0);

  ~MatVecProd_dRdX();

//  GenMat<Scalar,dim> &operator() (int i) { return *A[i]; }

  void initializeOperators(double);
  dRdXoperators<dim> *getdRdXop() {return dRdXop;}
  void constructOperators(Vec3D &, DistSVec<double,3> &X, DistVec<double> &ctrlVol, DistSVec<double,dim> &U,
                          double dMach, DistSVec<double,dim> &R, DistVec<double> &, DistTimeState<dim> *timeState, PostOperator<dim> *);


  void evaluate(int, DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &)
  {std::cout<<__FILE__<<":"<<__LINE__<<" empty header declaration called"; exit(-1);}

  void evaluate(DistExactRiemannSolver<dim> &, int, DistSVec<double,3> &, DistVec<double> &,
                DistSVec<double,dim> &, DistSVec<double,dim> &)
  {std::cout<<__FILE__<<":"<<__LINE__<<" empty header declaration called"; exit(-1);}

  void evaluate(int , DistSVec<double,3> &, DistVec<double> &, 
                DistSVec<double,dim> &, DistSVec<double,dim> &, Scalar)
  {std::cout<<__FILE__<<":"<<__LINE__<<" empty header declaration called"; exit(-1);}

  void evaluate2(int, DistSVec<double,3> &, DistVec<double> &,
                 DistSVec<double,dim> &, DistSVec<double,dim> &)
  {std::cout<<__FILE__<<":"<<__LINE__<<" empty header declaration called"; exit(-1);}

  // UH (09/10)
  // The following functions are never called and not implemented.
  //
  //void evaluatestep1(int, DistSVec<double,3> &, DistVec<double> &,
  //               DistSVec<double,dim> &, DistSVec<double,dim> &);
  //void evaluate2step1(int, DistSVec<double,3> &, DistVec<double> &,
  //               DistSVec<double,dim> &, DistSVec<double,dim> &);

  void evalH(int , DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &);

  void apply(DistSVec<double,neq>        &, DistSVec<double,neq> &)
  {std::cout<<__FILE__<<":"<<__LINE__<<" empty header declaration called"; exit(-1);}

  void apply(DistSVec<bcomp,neq>         &, DistSVec<bcomp,neq> &)
  {std::cout<<__FILE__<<":"<<__LINE__<<" empty header declaration called"; exit(-1);}

  void apply(DistEmbeddedVec<double,dim> &, DistEmbeddedVec<double,dim> &)
  {std::cout<<__FILE__<<":"<<__LINE__<<" empty header declaration called"; exit(-1);}

  void applyTranspose(DistSVec<double,neq> &, DistSVec<double,neq> &)
  {std::cout<<__FILE__<<":"<<__LINE__<<" empty header declaration called"; exit(-1);}

  void applyTranspose(DistEmbeddedVec<double,neq> &, DistEmbeddedVec<double,neq> &)
  {std::cout<<__FILE__<<":"<<__LINE__<<" empty header declaration called"; exit(-1);}

  void applyT(DistSVec<double,neq> &, DistSVec<double,neq> &)
  {std::cout<<__FILE__<<":"<<__LINE__<<" empty header declaration called"; exit(-1);}

  void applyT(DistSVec<bcomp,neq> &x, DistSVec<bcomp,neq> &y)
  {std::cout<<__FILE__<<":"<<__LINE__<<" empty header declaration called"; exit(-1);}

// Included (MB)
  void evaluateInviscid(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &)
  {std::cout<<__FILE__<<":"<<__LINE__<<" empty header declaration called"; exit(-1);}

  void evaluateViscous(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,dim> &)
  {std::cout<<__FILE__<<":"<<__LINE__<<" empty header declaration called"; exit(-1);}

  //void applyInviscid(DistSVec<double,neq> &, DistSVec<double,neq> &);
  //void applyViscous(DistSVec<double,neq> &, DistSVec<double,neq> &);

  void rstSpaceOp(IoData &, VarFcn *, SpaceOperator<dim> *, bool, SpaceOperator<dim> * = 0)
  {std::cout<<__FILE__<<":"<<__LINE__<<" empty header declaration called"; exit(-1);}

};
/*
//----------------------------------------------------------------------------//
//                MatVecProd for H1 Reduced Order Model                       //
//----------------------------------------------------------------------------//
// Approximate Matrix Vector Product for Reduced Order model
template<int dim, class Scalar, int neq>
class MatVecProdRomH1: public MatVecProdH1<dim, Scalar, neq> {
protected:
    int reducedDimension;
    VecSet<DistSVec<Scalar, dim> > Phi; //<! reduced bases

public:
    MatVecProdRomH1(DistTimeState<dim> *, SpaceOperator<dim> *, Domain *, IoData &, VecSet<DistSVec<Scalar, dim> > &);
    ~MatVecProdRomH1();
    void apply(Vec<double> &b, Vec<double> &x); //<! matrix-vector-multiplication for reduced
};
*/
//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <MatVecProd.C>
#endif

#endif

