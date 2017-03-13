#ifndef _KSP_SOLVER_H_
#define _KSP_SOLVER_H_

#include <cstdio>
#include <Vector.h>
#include <DenseMatrix.h>
#include <VectorSet.h>
#include <KspBinaryOutput.h>

class KspData;
class KspConvCriterion;

#include <complex>
typedef std::complex<double> bcomp;

#ifndef _KSPSLVR_TMPL_
#define _KSPSLVR_TMPL_
template<class VecType, class MvpOp, class PrecOp, class IoOp, class ScalarT = double> class KspSolver;
#endif

//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class ScalarT>
class KspSolver {

protected:

  double eps;
  double absoluteEps;
  bool checkFinalRes;

  MatVecProdOp *mvpOp;
  PrecOp *pcOp;
  IoOp *ioOp;

  KspConvCriterion *kspConvCriterion;


  FILE *output;

public:
  int maxits;
  int typePrec;

  KspSolver() {}
  KspSolver(KspData &, MatVecProdOp *, PrecOp *, IoOp *);
  virtual ~KspSolver() { if (kspConvCriterion) delete kspConvCriterion; }

  void setup(int, int, VecType &);

  void setMaxIts(int i) { maxits = i; }
  void setEps(double e) { eps = e; }
  void disableOutput() { output = NULL; } 

  virtual int solve(VecType &, VecType &) = 0;
  virtual int solveNew(VecType &, VecType &) = 0;//TODO BUGHUNT
  virtual int solveT(VecType &, VecType &) = 0;
  virtual int solveLS(VecType &, VecType &) { std::cout<<"*** ERROR solveLS bas routine called"<<std::endl; exit(-1); return 0; };

  void printParam() { ioOp->fprintf(stderr, " solver params: %d maxits, %e eps\n", maxits, eps);  }

  void setKspBinaryOutput(KspBinaryOutput<VecType>*);
  KspBinaryOutput<VecType> *kspBinaryOutput;
  void setTypePrec(int tp) { typePrec = tp; }

};

//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp>
class RichardsonSolver : public KspSolver<VecType,MatVecProdOp,PrecOp,IoOp> {

  VecType dx, r;

public:

  RichardsonSolver(const typename VecType::InfoType &, KspData &, 
		   MatVecProdOp *, PrecOp *, IoOp *);
  ~RichardsonSolver() {}

  int solve(VecType &, VecType &);
  int solveNew(VecType &, VecType &){return 0;};//TODO BUGHUNT
  int solveLS(VecType &, VecType &);
  int solveT(VecType &, VecType &) { return 0; }

};

//------------------------------------------------------------------------------

template<class VecType, class MatVecProdOp, class PrecOp, class IoOp>
class CgSolver : public KspSolver<VecType,MatVecProdOp,PrecOp,IoOp> {

  VecType r, Ap, y, p;

public:

  CgSolver(const typename VecType::InfoType &, KspData &,
	   MatVecProdOp *, PrecOp *, IoOp *);
  ~CgSolver() {}

  int solve(VecType &, VecType &);
  int solveNew(VecType &, VecType &){return 0;};//TODO BUGHUNT
  int solveLS(VecType &, VecType &);
  int solveT(VecType &, VecType &) { return 0;}
  int solveMRhs(VecType &, VecType &);

};

//------------------------------------------------------------------------------
//IoOp refers to the communicator here
template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class
		ScalarT = double>
class GmresSolver : public KspSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT> {

  int numVec;

  GenFullM<ScalarT> cs, H;
  Vec<ScalarT> g, y;

  VecSet<VecType> V;
  VecType w, r;

  bool outputConvergenceInfo;

public:

  GmresSolver(const typename VecType::InfoType &, KspData &, 
	      MatVecProdOp *, PrecOp *, IoOp *);
  ~GmresSolver() {}

  void disableConvergenceInfo();

  int solve(VecType &, VecType &);
  int solveNew(VecType &, VecType &);//TODO BUGHUNT
  int solveLS(VecType &, VecType &);
  int solveT(VecType &, VecType &);
  int solve(VecSet<VecType> &, VecSet<VecType> &);

  void applyPreviousRotations(int, GenFullM<ScalarT> &, GenFullM<ScalarT> &);
  void applyNewRotation(int, GenFullM<ScalarT> &, GenFullM<ScalarT> &, Vec<ScalarT> &);
  void backwardSolve(int, GenFullM<ScalarT> &, Vec<ScalarT> &, Vec<ScalarT> &);
  
};

//------------------------------------------------------------------------------
                                                        
template<class VecType, class MatVecProdOp, class PrecOp, class IoOp, class
                ScalarT = double>
class GcrSolver : public KspSolver<VecType,MatVecProdOp,PrecOp,IoOp, ScalarT> {
                                                        
                                                        
                                                        
  int numVec;
                                                        
  ScalarT alpha, beta;
                                                        
  ScalarT *ApAp;
                                                        
  ScalarT *y;                                                                                                                  
  VecSet<VecType> p, Ap;
                                                        
  VecType w, r, R, AR, temp, w0, x0;
                                                        
                                                        
                                                        
public:
                                                        
  int numCalcVec;                                                                                                                  
  GcrSolver(const typename VecType::InfoType &, KspData &,
              MatVecProdOp *, PrecOp *, IoOp *);
  ~GcrSolver() {}
                      
  int solve(VecType &, VecType &);
  int solveNew(VecType &, VecType &){return 0;};//TODO BUGHUNT
  int solveT(VecType &, VecType &) { return 0; }
  int solveMRhs(VecType &, VecType &);
                                                        
};
                      
//------------------------------------------------------------------------------


#ifdef TEMPLATE_FIX
#include <KspSolver.C>
#endif

#endif
