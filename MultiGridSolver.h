#pragma once

#include<IoData.h>

class IoData;

//------------------------------------------------------------------------------
/** Class which handles the algorithmic organization of the solution for all problems */
template<class ProblemDescriptor>
class MultiGridSolver {

  ProblemDescriptor *probDesc;

  int resolve(typename ProblemDescriptor::SolVecType &,
              IoData &);

  int num_presmooth,num_postsmooth;
  int num_finesweeps;

 
//  typename ProblemDescriptor::MultiGridKernelType* multiGridKernel; 

//  typedef typename ProblemDescriptor::MultiGridKernelType::MultiGridSmoother
//    MultiGridSmoother;

  typename ProblemDescriptor::SolVecType* U_vec_smooth;

public:

  MultiGridSolver(ProblemDescriptor *);
  ~MultiGridSolver() {}

  int solve(IoData &);

// Included (MB)
  int fsaSolve(IoData &);

};

//------------------------------------------------------------------------------

template<class ProblemDescriptor>
MultiGridSolver<ProblemDescriptor>::MultiGridSolver(ProblemDescriptor *prbd)
{

  probDesc = prbd;

}

//------------------------------------------------------------------------------

template<class ProblemDescriptor>
int MultiGridSolver<ProblemDescriptor>::solve(IoData &ioData)
{
  typename ProblemDescriptor::SolVecType U(probDesc->getVecInfo());

  // initialize solutions and geometry
  probDesc->setupTimeStepping(&U, ioData);
//  mySmoother = new MySmoother(probDesc->getMultiGridKernel(),probDesc);
  
  int status = resolve(U, ioData);
  return status;

}

//------------------------------------------------------------------------------

template<class ProblemDescriptor>
int MultiGridSolver<ProblemDescriptor>::resolve(typename ProblemDescriptor::SolVecType &U,
                                                IoData &ioData)
{
  bool lastIt = false;
  
  U_vec_smooth = &U;

  if (probDesc->getSmoothedVec()) {

    U_vec_smooth = probDesc->getSmoothedVec();
  }

  // dts is structural time step
  double dt, dts;
  int it = probDesc->getInitialIteration();
  double t = probDesc->getInitialTime();
  // setup solution output files
  probDesc->setupOutputToDisk(ioData, &lastIt, it, t, U);

  /** for embedded method: send force (if it>0) and receive disp (from Struct). */
  dts = probDesc->computePositionVector(&lastIt, it, t, U);

  if (lastIt)
    probDesc->outputPositionVectorToDisk(U);
    
  // compute control volumes and velocities
  probDesc->computeMeshMetrics();
    
  typename ProblemDescriptor::SolVecType zero(probDesc->getVecInfo());
  zero = 0.0;
/*
  multiGridKernel = probDesc->getMultiGridKernel();
  multiGridKernel->setParameters(1,0,2,1.0,1);
  multiGridKernel->initialize();
  multiGridKernel->setGeometric();
*/
  while (!lastIt) {
    probDesc->resetOutputToStructure(U);
    it++;
    
    bool solveOrNot = true;
    probDesc->setCurrentTime(t,U);
  
    double rf = 1.0;//std::min((double)it/10.0, 0.5);
//    multiGridKernel->setRestrictRelaxFactor( rf );
//    multiGridKernel->setProlongRelaxFactor( rf );

//    probDesc->printf( 0, "Current Restrict Relax factor = %lf\n", rf);

    probDesc->cycle(U);
//    multiGridKernel->cycleV(zero, U);
    //multiGridKernel->cycleW(0,zero, U);
    probDesc->updateStateVectors(U, 0);

    t += 1.0;

    // compute the current aerodynamic force
    if (it % 10 == 0) {
      probDesc->updateOutputToStructure(0.0, 0.0, U);
    }

// Modified (MB)
    lastIt = probDesc->checkForLastIteration(ioData, it, t, 0.00, U);

    if (it % 1 == 0) {
      probDesc->outputForces(ioData, &lastIt, it, 1, 1, t, 0.0, *U_vec_smooth);
      dts = probDesc->computePositionVector(&lastIt, it, t, U);

      probDesc->outputToDisk(ioData, &lastIt, it, 1, 1, t, 0.0, U);
    }
  }

  probDesc->outputForces(ioData, &lastIt, it, 1, 1, t, 0.0, *U_vec_smooth);
  dts = probDesc->computePositionVector(&lastIt, it, t, U);

  probDesc->outputToDisk(ioData, &lastIt, it, 1, 1, t, 0.0, U);
  return 0;

}

//------------------------------------------------------------------------------
