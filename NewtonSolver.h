#ifndef _NEWTON_SOLVER_H_
#define _NEWTON_SOLVER_H_

#include <ErrorHandler.h>

#include <cstdlib>
#include <cmath>

/*@ARTICLE{brown-saad-90,
  author = "Brown, P. N. and Saad, Y.",
  title = "Hybrid {K}rylov Methods for
           Nonlinear Systems of Equations",
  journal = siamjscistat,
  year = 1990,
  volume = 11,
  number = 3,
  pages = "450--481",
}
@ARTICLE{keyes-venkatakrishnan-96,
  author = "Keyes, D. E. and Venkatakrishnan, V.",
  title = "{N}ewton-{K}rylov-{S}chwarz Methods: Interfacing Sparse Linear
  Solvers with Nonlinear Applications",
  journal = zamm,
  year = 1996,
  volume = 76,
  pages = "147--150",
} 
*/
//------------------------------------------------------------------------------

template<class ProblemDescriptor>
class NewtonSolver {

  ProblemDescriptor *probDesc;

  typename ProblemDescriptor::SolVecType F;  // nonlinear function
  typename ProblemDescriptor::SolVecType Finlet;  // nonlinear function at inlet nodes
  typename ProblemDescriptor::SolVecType dQ; // gradient of F
  typename ProblemDescriptor::SolVecType rhs; // right hand side
  typename ProblemDescriptor::PhiVecType dPhi; // 
  typename ProblemDescriptor::PhiVecType PhiF; // 
  typename ProblemDescriptor::PhiVecType rhsPhi; // 

public:

  NewtonSolver(ProblemDescriptor *);
  ~NewtonSolver() {}
  int solve(typename ProblemDescriptor::SolVecType &, const int timeStep = 0, const double time = 0.0);
  int solveLS(typename ProblemDescriptor::PhiVecType &, typename ProblemDescriptor::SolVecType &);

  typename ProblemDescriptor::SolVecType* GetResidual() { return &F; }

};

//------------------------------------------------------------------------------

template<class ProblemDescriptor>
NewtonSolver<ProblemDescriptor>::NewtonSolver(ProblemDescriptor *prbd) : 
  F(prbd->getVecInfo()), Finlet(prbd->getVecInfo()), dQ(prbd->getVecInfo()),
  rhs(prbd->getVecInfo()), dPhi(prbd->getVecInfo()), PhiF(prbd->getVecInfo()), 
  rhsPhi(prbd->getVecInfo())
{

  probDesc = prbd;

}

//------------------------------------------------------------------------------------------------------

template<class ProblemDescriptor>
int
NewtonSolver<ProblemDescriptor>::solve(typename ProblemDescriptor::SolVecType &Q , const int timeStep, const double t)
{

  double res, target;
  double qnorm;

  int fsIt = 0;
  int maxIts = probDesc->getMaxItsNewton();
  double eps = probDesc->getEpsNewton();
  double epsAbsRes = probDesc->getEpsAbsResNewton();
  double epsAbsInc = probDesc->getEpsAbsIncNewton();
  FILE *output = probDesc->getOutputNewton();
  bool finalRes = true; // if true, then a final residual evaluation will be done
                        // when maxIts is reached before terminating the loop.
  double res0, res2=0.0;

  double rho; //contraction factor for backtracking
  double c1; //sufficient decrease factor for backtracking
  double restrial, res2trial=0.0;
  double alpha;
  int maxItsLS;

	if(probDesc->getLineSearch()) 
	{
   rho = probDesc->getContractionLineSearch();
   c1 = probDesc->getSufficientDecreaseLineSearch();
   maxItsLS = probDesc->getMaxItsLineSearch(); 
  }
  int it, itLS;
  bool converged = false;

	for(it=0; finalRes || it<maxIts; ++it) 
	{

    *(probDesc->getNewtonIt()) = it;
    *(probDesc->getNumResidualsOutputCurrentNewtonIt()) = 0;


    // compute the nonlinear function value
    probDesc->computeFunction(it, Q, F);
    res2 = probDesc->recomputeResidual(F, Finlet);
    res = F*F-res2;
    probDesc->printf(9, " ... deubgging: it = %d, res2 = %f, res = %f\n", it, res2, res);
		if(res < 0.0)
		{
      probDesc->printf(1, "ERROR: negative residual captured in Newton Solver!\n");
      exit(1);
    }

    // UH (08/10) After the test, it is safe to take the square root.
    //res = sqrt(F*F-res2);
    res = sqrt(res);

		if(it == 0) 
		{
      target = eps*res; 
      res0 = res;
    }

    if(output) probDesc->fprintf(output,"Newton residual = %e, target = %e\n",res,target);
    if (res == 0.0 || res <= target) { converged = true; break; }
    if (it > 0 && res <= epsAbsRes && dQ.norm() <= epsAbsInc) { converged = true; break; } // alternative stopping criterion
    if (it == maxIts) break;

    rhs = -1.0 * F;

    if (probDesc->outputOnlySpatialResidual()) {
      typename ProblemDescriptor::SolVecType spatialRes(probDesc->getVecInfo());
      probDesc->calculateSpatialResidual(Q,spatialRes);
      probDesc->writeBinaryVectorsToDiskRom(false, timeStep, it, &Q, &spatialRes);
    } else {
      probDesc->writeBinaryVectorsToDiskRom(false, timeStep, it, &Q, &F);
    }

    probDesc->recomputeFunction(Q, rhs); // only for implicit level set methods; ignored in all other cases

    probDesc->computeJacobian(it, Q, F);

    // apply preconditioner if available
    probDesc->setOperators(Q);
    // set up krylov snapshots for ROM if applicable
    probDesc->setCurrentStateForKspBinaryOutput(Q);

    probDesc->solveLinearSystem(it, rhs, dQ);

    probDesc->printf(9, " debugging NewtonSolver::solve() : maxItsLs = %d, getLineSearch = %d\n", probDesc->getMaxItsLineSearch(), probDesc->getLineSearch());
   if (probDesc->getLineSearch()) { 
     for (itLS=0; itLS<maxItsLS; ++itLS) {
       if (itLS>0){
         alpha *= rho; 
         if (itLS==1)
           dQ *= (rho-1);
         else
           dQ *= rho;
       }
       else 
         alpha = 1.0;
       // increment or backtract from previous trial 
       probDesc->fixSolution(Q, dQ); // dQ[i] = negative_pressure_or_density(Q[i] + dQ[i]) ? 0.0 : dQ[i];
       // compute updated residual
       rhs = Q;
       Q += dQ;
       probDesc->computeFunction(it, Q, F);
       res2trial = probDesc->recomputeResidual(F, Finlet);
       restrial = F*F-res2trial;
       probDesc->printf(9, " debugging Newtonsolver::solve(): dQ.norm() = %f, res2trial = %f, restrial = %f\n", dQ.norm(), res2trial, restrial);
       if (restrial>=0.0) {
         if (sqrt(restrial) < sqrt(1-2.0*alpha*c1)*res || dQ.norm() <= epsAbsInc)
           break;
       }
       if (itLS == maxItsLS-1 && maxItsLS != 1) 
         probDesc->printf(1, "*** Warning: Line Search reached %d its ***\n", maxItsLS);
     }
   }
   else { 
// Included (MB)
      probDesc->fixSolution(Q, dQ);
      rhs = Q;
      Q += dQ;
   }
    probDesc->incrementNewtonOutputTag();

    // verify that the solution is physical
    if (probDesc->checkSolution(Q)) {
      if (probDesc->getTsParams()->checksol) {
        //probDesc->getErrorHandler()->localErrors[ErrorHandler::REDO_TIMESTEP] += 1;
        probDesc->checkFailSafe(Q);
        Q = rhs;
        --it;
        ++fsIt;
        return -10; // signal re-compute CFL number
      }
      else if (probDesc->checkFailSafe(Q) && fsIt < 5) {
        probDesc->printf(1, "*** Warning: Newton solver redoing iteration %d\n", it+1);
        Q = rhs;
        --it;
        ++fsIt;
      }
			else
			{
        probDesc->printf(1, "Newton solver failed\n");
        return -3;
      }
    }

  } // for (it=0; it<maxIts; ++it)

	if(fsIt > 0 && probDesc->checkFailSafe(Q) == 1)	probDesc->resetFixesTag();

	if(!converged && maxIts != 1) 
	{
    probDesc->printf(1, "*** Warning: Newton solver reached %d its", maxIts);
    probDesc->printf(1, " (Residual: initial=%.2e, reached=%.2e, target=%.2e)\n", res0, res, target);    
  }

  probDesc->writeBinaryVectorsToDiskRom(true, timeStep, it, &Q, NULL);  // save state after final iteration for ROM
  
  return it;

}

//------------------------------------------------------------------------------
template<class ProblemDescriptor>
int
NewtonSolver<ProblemDescriptor>::solveLS(typename ProblemDescriptor::PhiVecType &Phi,
					 typename ProblemDescriptor::SolVecType &U )
{

  double res, target, res1;

  int fsIt = 0;
  int maxIts = probDesc->getMaxItsNewton();
  double eps = probDesc->getEpsNewton();
  int it;
  for (it=0; it<maxIts; ++it) {

    // compute the nonlinear function value for the Level Set Equation
    probDesc->computeFunctionLS(it, U, PhiF);
    res = sqrt(PhiF*PhiF);
    //probDesc->printf(1,"Newton residual = %e,target = %e\n",res,target);
    if (it == 0){
      target = eps*res;
      res1  = res;
    }

    if (res == 0.0 || res <= target) break;

    rhsPhi = -1.0* PhiF;

    // compute the Jacobian for the Level Set Equation 
    probDesc->computeJacobianLS(it, U, PhiF);
 
    // Apply preconditioner
    probDesc->setOperatorsLS(Phi);
 
    // Solve the linearized system of equations
    probDesc->solveLinearSystemLS(it, rhsPhi, dPhi);

    rhsPhi  = Phi;
    Phi += dPhi;

  }
  if (it == maxIts && maxIts != 1) {
    probDesc->printf(1, "*** Warning: Newton solver for LS reached %d its", maxIts);
    probDesc->printf(1, " (initial=%.2e, res=%.2e, target=%.2e)\n", res1, res, target);
  }

  return it;
}
//------------------------------------------------------------------------------
#endif
