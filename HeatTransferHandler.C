#include <HeatTransferHandler.h>

#include <PostOperator.h>

//------------------------------------------------------------------------------

template<int dim>
void HeatTransferHandler::updateOutputToStructure(double dt, double dtLeft,
						  PostOperator<dim>* postOp,
						  DistSVec<double,3>& X, 
						  DistSVec<double,dim>& U)
{
  if (dtLeft == 0.0)
    postOp->computeNodalHeatPower(X, U, P);

}

//------------------------------------------------------------------------------
