#ifndef _IMPLICIT_TS_DESC_H_
#define _IMPLICIT_TS_DESC_H_

#include <IoData.h>
#include <TsDesc.h>
#include <KspPrec.h>
#include <KspBinaryOutput.h>

struct DistInfo;

class GeoSource;
class Domain;
class Communicator;

template<class Scalar, int dim> class DistSVec;
template<int dim, int neq> class MatVecProd;
template<class ProblemDescriptor> class NewtonSolver;
#ifndef _KSPSLVR_TMPL_
#define _KSPSLVR_TMPL_
template<class VecType, class MvpOp, class PrecOp, class IoOp, class ScalarT = double> class KspSolver;
#endif


//------------------------------------------------------------------------------

template<int dim>
class ImplicitTsDesc : public TsDesc<dim> {

protected:

  int failSafeNewton;
  int maxItsNewton;
  double epsNewton;
  double epsAbsResNewton, epsAbsIncNewton;
  FILE *outputNewton;
  int maxItsLS;
  double contractionLS, sufficDecreaseLS;

  DistSVec<bool,2> *tag;

  NewtonSolver<ImplicitTsDesc<dim> > *ns;

  IoData* myIoDataPtr;

protected:

  template<class Scalar, int neq>
  KspPrec<neq> *createPreconditioner(PcData &, Domain *);

  template<int neq>
  KspSolver<DistSVec<double,neq>, MatVecProd<dim,neq>, KspPrec<neq>, 
    Communicator> *createKrylovSolver(const DistInfo &, KspData &, MatVecProd<dim,neq> *, 
				      KspPrec<neq> *, Communicator *);

  KspBinaryOutput<DistSVec<double, dim> > *kspBinaryOutput; 
  DistSVec<double,dim> *prevOutputState;
  DistSVec<double,dim> *stateIncrement;

public:
  
  ImplicitTsDesc(IoData &, GeoSource &, Domain *);
  ~ImplicitTsDesc();

  virtual void computeJacobian(int, DistSVec<double,dim> &, DistSVec<double,dim> &) = 0;
  virtual void setOperators(DistSVec<double,dim> &) = 0;
  virtual int solveLinearSystem(int, DistSVec<double,dim> &, DistSVec<double,dim> &) = 0;
  int solveNonLinearSystem(DistSVec<double,dim> &, const int);
  void computeFunction(int, DistSVec<double,dim> &, DistSVec<double,dim> &);  
  void recomputeFunction(DistSVec<double,dim> &, DistSVec<double,dim> &);
  void doErrorEstimation(DistSVec<double,dim> &);

  int checkFailSafe(DistSVec<double,dim>&);
  void resetFixesTag();

  int getMaxItsNewton() const { return maxItsNewton; }
  double getEpsNewton() const { return epsNewton; }
  double getEpsAbsResNewton() const { return epsAbsResNewton; }
  double getEpsAbsIncNewton() const { return epsAbsIncNewton; }
  FILE* getOutputNewton() const { return outputNewton; }
  int getLineSearch() const { return (maxItsLS>0); }
  int getMaxItsLineSearch() const { return maxItsLS; }
  double getContractionLineSearch() const { return contractionLS; }
  double getSufficientDecreaseLineSearch() const { return sufficDecreaseLS; }

  void writeBinaryVectorsToDiskRom(bool lastNewtonIt, int timeStep, int newtonIt,
                                   DistSVec<double,dim> *state = NULL, DistSVec<double,dim> *residual = NULL);
  void calculateSpatialResidual(DistSVec<double,dim> &, DistSVec<double,dim> &); 
  void incrementNewtonOutputTag();
  void setCurrentStateForKspBinaryOutput(DistSVec<double,dim> &Q) {if (kspBinaryOutput) kspBinaryOutput->setCurrentState(Q);}

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <ImplicitTsDesc.C>
#endif

#endif
