//
// Created by lei on 5/16/16.
// for running with _EMBEDDED_ALS_ROM_ONLINE_

#ifndef PROJECT_IMPLICITEMBEDDEDROMTSDESC_H
#define PROJECT_IMPLICITEMBEDDEDROMTSDESC_H

#include <MatVecProd.h>                             // for Jacobian
#include <ImplicitEmbeddedCoupledTsDesc.h>          // base calss, time integrator
#include <VectorSet.h>                              // vector and matrix class
#include <ParallelRom.h>                            // for scalapack least square
#include <EmbeddedAlternatingLeastSquare.h>         // for reading basis

// TODO: inherit from ImplicitEmbeddedTsDesc instead
// TODO: inherit from ImplictRomTsDesc instead
// see ImplicitEmbeddedCoupledTsDesc.C
template<int dim>
class ImplicitEmbeddedRomTsDesc : public ImplicitEmbeddedCoupledTsDesc<dim> {
private:
    typedef ImplicitEmbeddedCoupledTsDesc<dim> super;
    ParallelRom<dim> *LeastSquareSolver;
    // two jacobians here, one for ROM and one for FOM
    // MatVecProd<dim, dim> *test_Jacobian; // not needed, use parent mvp object
    MatVecProd<dim, dim> *Jacobian;
    KspPrec<dim> *rom_pc;
    KspSolver<DistEmbeddedVec<double,dim>, MatVecProd<dim,dim>, KspPrec<dim>, Communicator> *rom_ksp;

    VecSet<DistSVec<double, dim> > reducedJacobian;
    Vec<double> reducedNewtonDirection;
    //TODO: two methods: project U into Qy and lift y to U

    int reducedDimension;
    VecSet<DistSVec<double, dim> > reducedBasis;
    DistSVec<double, dim> referenceState;
    double **result; //<! scratchpad for scalapack to store result, fixed size from initialization

    // internal methods
    void printGhostPoint();
    void projectStateOntoROB(DistSVec<double, dim> &U);
    void expandVector(Vec<double>& p, DistSVec<double, dim>& dQ);
    void projectVector(VecSet<DistSVec<double, dim> > &mat, DistSVec<double,dim> &vec, Vec<double> &buffer);

    // misc variables to pass c++ compiler hoops.
    EmbeddedAlternatingLeastSquare<dim> embeddedALS;

public:
    /** @name interface to NewtonSolver
     * REQUIRED functions to call NewtonSolver (backtracking linear search)
     * listed here explicitly for clarity
     */
    ///@{
    int getMaxItsNewton() { return super::getMaxItsNewtion(); };
    double getEpsNewton() { return super::getMaxItsNewtion(); };
    double getEpsAbsResNewton() { return super::getEpsAbsResNewton(); };
    double getEpsAbsIncNewton() {return super::getAbsIncNewton(); };
    FILE* getOutputNewton() { return super::getOutputNewton(); };
    bool getLineSearch() { return super::getLineSearch(); };
    double getContractionLinSearch() {return super::getContractionLinSearch(); };
    double getSufficientDecreaseLinSearch() {return super::getSufficientDecreaseLinSearch();};
    double getMaxItsLineSearch() {return super::getMaxItsLineSearch(); };
    double getNewtonIt(){return super::getNewtonIt(); };
    double getNumResidualOutputCurrentNewtonIt() {return super::getNumResidualOutputCurrentNewtonIt(); };
    // equivalent to computeFullResidual(), result in F ?
    void computeFunction(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> &F) { super::computeFunction(it, Q, F); };
    void recomputeResidual(DistSVec<double, dim> &F, DistSVec<double, dim> &Finlet) { super::recomputeResidual(F, Finlet); };
    // add inlet contribution, result in rhs ?
    void recomputeFunction(DistSVec<double, dim> &Q, DistSVec<double, dim> &rhs){ super::recomputeFunction(Q, rhs); };
    void setOperators(DistSVec<double, dim> &) {};
    // TODO: emulate coupledTsDesc (compute Hessian of function)
    void computeJacobian(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> &F);
    // TODO: emulate solveNewtonSystem() (get search direction in dQ)
    int solveLinearSystem(int it, DistSVec<double, dim> &rhs, DistSVec<double, dim> &dQ);
    ///@}
    int _solveLinearSystem(int it , DistSVec<double, dim> &rhs, Vec<double> &sol);

    /** @name interface to NewtonSolver
     * OPTIONAL functions to call NewtonSolver (backtracking linear search)
     */
    ///@{
    //void fprintf(FILE *, const char *, ...); // already in TsDesc.h
    void writeBinaryVectorsToDiskRom(bool, int, DistSVec<double, dim> &, DistSVec<double, dim> &) {};
    void setCurrentStateForKspBinaryOutput(DistSVec<double, dim> &Q) { super::setCurrentStateForKspBinaryOutput(Q); };
    int checkSolution(DistSVec<double, dim> &Q) {return super::checkSolution(Q); };
    void fixSolution(DistSVec<double, dim> &Q, DistSVec<double, dim> &dQ) { super::fixSolution(Q, dQ); };
    bool checkFailSafe(DistSVec<double, dim> &Q) {return super::checkFailSafe(Q); };
    void resetFixesTag() { super::resetFixesTag(); };
    ///@}

    /** @name interface to TsSolver
     * REQUIRED functions to run within TsSolver
     * override base methods
     */
    ///@{
    //TODO: overwrite this to use reduced coordinate ?
    int solveNonlinearSystem(DistSVec<double, dim> &Q, const int totalTimeSteps) { return super::solveNonLinearSystem(Q, totalTimeSteps); };
    ///@}
    // diagnostic/test function
    void test();
    void checkMVP(int it, DistSVec<double, dim> &Q, DistSVec<double, dim> &dx, DistSVec<double, dim> &orig);
    void outputMatrix(char *fn, VecSet<DistSVec<double, dim> > &A);
    void outputVector(char *fn, DistSVec<double, dim> &b);
    /*
     * rom krylov subspace solver
     * todo: make ksprec for vec
     * todo: make setup() work with vec<double>
     */
    //template<int dim, int dimLS, class MatVecProdOp>
    //KspSolver<Vec<double>, MatVecProdOp, KspPrec, Communicator> * createROMKrylovSolver(KspData &kspdata, MatVecProdOp *_mvp, KspPrec *_pc, Communicator *_com);

    /* must implement, or compilation error
    void solveNewtonSystem(const int &it, double &res, bool &breakloop,
                           DistSVec<double, dim> &U,
                           const int &totalTimeSteps = 0);
    */
    ImplicitEmbeddedRomTsDesc(IoData &, GeoSource &, Domain *);
    ~ImplicitEmbeddedRomTsDesc();
};

#ifdef TEMPLATE_FIX
#include <ImplicitEmbeddedRomTsDesc.cpp>
#endif

#endif //PROJECT_IMPLICITEMBEDDEDROMTSDESC_H
