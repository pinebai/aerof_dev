#include <cstdlib>
#include <cmath>

#include <MeshMotionHandler.h>

#include <Domain.h>
#include <PostOperator.h>
#include <StructExc.h>
#include <MeshMotionSolver.h>

//------------------------------------------------------------------------------
// TODO variable time-step

template<int dim>
void MeshMotionHandler::computeInterfaceWork(double dt, PostOperator<dim>* postOp, 
					     DistSVec<double,3>& Xn, DistSVec<double,dim>& Un, 
					     DistSVec<double,3>& Xnp1, DistSVec<double,dim>& Unp1, 
					     double* Wnp1, DistVec<int> * fluidIdn, DistVec<int> *fluidIdnp1)
{
  // computing energy by method 1
  postOp->computeNodalForce(Xn, Un, Pin, dX, fluidIdn);   // These are the same forces that are transmitted to structure
  postOp->computeNodalForce(Xnp1, Unp1, Pin, F, fluidIdnp1);
  F = 0.5 * (dX + F);
  // F needs to be assembled in order to get the correct scalar product
  domain->assemble(domain->getVec3DPat(), F);
  Wnp1[0] = F * (Xn - Xnp1);

// This method is trying to shift the output to the same time step as the structure, but is not accurate
//  Wnp1[1] = 1./3.*Wn + 2./3.*dt* postOp->computeInterfaceWork(Xnp1, Unp1, Pin);
//  Wn = Wnp1[1];
    Wnp1[1] = 0.0; 
}

//------------------------------------------------------------------------------
// note: this initialization is not correct in case of a restart with an averaged force

template<int dim>
void AeroMeshMotionHandler::setup(int *rstrt, double *maxTime, PostOperator<dim>* postOp,
				  DistSVec<double,3> &X, DistSVec<double,dim> &U,
				  DistVec<int> *fluidId)
{

  postOp->computeNodalForce(X, U, Pin, F, fluidId);
  
  *rstrt = strExc->getRestartFrequency();
  *maxTime = strExc->getMaxTime();

  mms->setup(X);

}

//------------------------------------------------------------------------------

template<int dim>
void AeroMeshMotionHandler::resetOutputToStructure(PostOperator<dim>* postOp,
						   DistSVec<double,3> &X, 
						   DistSVec<double,dim> &U,
						   DistVec<int> *fluidId)
{

  if (forceComputation == AeroelasticData::AVERAGED) {
    *Favg = 0.0;
    if (timeIntegrator == IMPLICIT_SECOND_ORDER)
      postOp->computeNodalForce(X, U, Pin, *Fn, fluidId);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void AeroMeshMotionHandler::updateOutputToStructure(double dt, double dtLeft,
						    PostOperator<dim>* postOp,
						    DistSVec<double,3> &X, 
						    DistSVec<double,dim> &U,
						    DistVec<int> *fluidId)
{

  if (forceComputation == AeroelasticData::AVERAGED) {
    switch (timeIntegrator) {
      case IMPLICIT_FIRST_ORDER : {
        postOp->computeNodalForce(X, U, Pin, F, fluidId);
        break;
      }
      case IMPLICIT_SECOND_ORDER : {
        postOp->computeNodalForce(X, U, Pin, *Fnp1, fluidId);
        F = 0.5 * (*Fn + *Fnp1);
        *Fn = *Fnp1;
        break;
      }
    }
    *Favg += dt * F;
  }

  if (dtLeft == 0.0) {
    if (forceComputation == AeroelasticData::AVERAGED)
      F = (1.0 / strExc->getTimeStep()) * (*Favg);
    else if (forceComputation == AeroelasticData::LAST)
      postOp->computeNodalForce(X, U, Pin, F, fluidId);
    else if (forceComputation == AeroelasticData::LAST_KRIS) {
      postOp->computeNodalForce(X, U, Pin, F, fluidId);
      F = 1./3.*(*Fn) + 2./3.*F;
      *Fn = F;
    }
  }

}

