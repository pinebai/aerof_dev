/*
 * FluidShapeOPtimizationHandler_temp.C
 *
 *  Created on: Mar 12, 2017
 *      Author: lscheuch
 */

//TODO this is probably wrong
template<int dim>
void FluidShapeOptimizationHandler<dim>::fsoGetDerivativeOfEffortsAnalytical(
                                     bool isSparse,
                     IoData &ioData,
                                         DistSVec<double,3> &X,   //mesh motion
                     DistSVec<double,3> &dX,  //derivative of mesh motion
                                         DistSVec<double,dim> &U, //state vector
                     DistSVec<double,dim> &dU,//derivative of state vector
                                         Vec3D &dForces,          //derivative of forces
                     Vec3D &dMoments,         //derivative of moments
                     Vec3D &dL)               //derivative of lift/drag
{


  double gamma     = ioData.eqs.fluidModel.gasModel.specificHeatRatio;

  //TODO I am not sure if that is correct, even if the refernce state is set to inflow, I dont think one has to differentiate here
  //!! The following formulas assume that the reference state is equal to the inflow parameters
  //formula for reference velocity(see Manual)
  // velocity
  double velocity  = ioData.ref.mach * sqrt(gamma * ioData.ref.pressure / ioData.ref.density);

  //dericative of the reference velocity with respect to the reference Mach number//TODO I think this is wrong if reference is not equal to inflow
  double dVelocity = sqrt(gamma * ioData.ref.pressure / ioData.ref.density)*DFSPAR[0];//0 exept for mach sensitivity

  double dForce =2.0*ioData.ref.density*ioData.ref.length*ioData.ref.length*velocity*dVelocity;//0 exept for mach sensitivity
  double dEnergy=2.0*ioData.ref.density*ioData.ref.length*
                  ioData.ref.length*ioData.ref.length*velocity*dVelocity;//0 exsptp for mach snesitivity

  int nSurfs = this->postOp->getNumSurf();
  std::cout<<"NUMSURF: "<<nSurfs<<std::endl;//TODO delete line

  Vec3D x0, F, dF, M, dM;
  Vec3D Fdim;//dimensionalized version of F

  Vec3D *Fi = new Vec3D[nSurfs]; //internal forces for all surfaces
  Vec3D *Mi = new Vec3D[nSurfs]; //internal moments for all surfaces
  Vec3D *Fv = new Vec3D[nSurfs]; //transmitted forces for all surfaces
  Vec3D *Mv = new Vec3D[nSurfs]; //transmitted moments for all surfaces

  Vec3D *dFi = new Vec3D[nSurfs]; //derivative of internal forces for all surfaces
  Vec3D *dMi = new Vec3D[nSurfs]; //derivative of internal moments for all surfaces
  Vec3D *dFv = new Vec3D[nSurfs]; //derivative of transmitted forces for all surfaces
  Vec3D *dMv = new Vec3D[nSurfs]; //derivative of transmitted moments for all surfaces

  x0[0] = ioData.output.transient.x0; //x-coordinate of the point around which the moments are computed
  x0[1] = ioData.output.transient.y0; //y-coordinate of the point around which the moments are computed
  x0[2] = ioData.output.transient.z0; //z-coordinate of the point around which the moments are computed

  //compute pressure gradient //TODO check if description is accurate
  //it seems, that the computed grad P is only stored at the subdomain itself and not returned lscheuch-12/2016
  this->spaceOp->computeGradP(X, *this->A, U);

  //computes Fi, Mi, Fv, Mv
  this->postOp->computeForceAndMoment(x0,  //reference point for moment calculation
                                  X,   //mesh motion
                    U,   //fluid state vector
                    0,   //fluid ID
                    Fi,  //internal forces
                    Mi,  //internal moments
                    Fv,  //transmitted forces
                    Mv); //transmitted moments

  F = 0.0;
  M = 0.0;

  F = Fi[0] + Fv[0]; //Add transmited and internal forces of first surface
  M = Mi[0] + Mv[0]; //Add transmitted and internam moments of first surface

  dRdXoperators<dim> *dRdXop = dRdX->getdRdXop();
//  this->spaceOp->computeDerivativeOfGradP(dRdXop, dX, dAdS, dU, dddx, dddy, dddz, dR, dGradP);
//  this->postOp->computeDerivativeOfForceAndMoment(x0, X, dX, U, dU, DFSPAR, dFi, dMi, dFv, dMv);


  if(isSparse) this->spaceOp->computeDerivativeOfGradP(dRdXop, dX, dAdS, dU, dddx, dddy, dddz, dR, dGradP);
  else this->spaceOp->computeDerivativeOfGradP(X, dX, *this->A, dAdS, U, dU);


  //computes dFi, dMi, dFv, dMv
  if(isSparse) {
    this->postOp->computeDerivativeOfForceAndMoment(dRdXop, dX, dU, DFSPAR, dGradP, dFi, dMi, dFv, dMv);
/* Verificaiton
    this->postOp->computeDerivativeOfForceAndMoment(x0, X, dX, U, dU, DFSPAR, dFi2, dMi2, dFv2, dMv2);
    for(int i=0; i<3; ++i) {
      double diff = abs(dFi2[0][i] - dFi[0][i]);
      if(dFi2[0][i] != 0) this->com->fprintf(stderr, "diff for dFi is %e\n", diff/abs(dFi2[0][i]));
      else this->com->fprintf(stderr, "diff for dFi is %e\n", diff);
      diff = abs(dMi2[0][i] - dMi[0][i]);
      if(dMi2[0][i] != 0) this->com->fprintf(stderr, "diff for dMi is %e\n", diff/abs(dMi2[0][i]));
      else this->com->fprintf(stderr, "diff for dMi is %e\n", diff);
      diff = abs(dFv2[0][i] - dFv[0][i]);
      if(dFv2[0][i] != 0) this->com->fprintf(stderr, "diff for dFv is %e\n", diff/abs(dFv2[0][i]));
      else this->com->fprintf(stderr, "diff for dFv is %e\n", diff);
      diff = abs(dMv2[0][i] - dMv[0][i]);
      if(dMv2[0][i] != 0) this->com->fprintf(stderr, "diff for dMv is %e\n", diff/abs(dMv2[0][i]));
      else this->com->fprintf(stderr, "diff for dMv is %e\n", diff);
    }
*/
  } else this->postOp->computeDerivativeOfForceAndMoment(x0, X, dX, U, dU, DFSPAR, dFi, dMi, dFv, dMv);


  dF = 0.0;               //derivative of force on first surface
  dM = 0.0;               //derivative of moment on first surface

  dF = dFi[0] + dFv[0];
  dM = dMi[0] + dMv[0];

  //forces and derivatives are non-dimensional at this point


  //Calculating the Lift out of the forces NEW version
  dL = 0.0;
  dForces2dLifts(ioData,dF,dL);
  dL *= this->refVal->force;
  std::cout<<"Force reference value: "<<this->refVal->force<<std::endl;
  ///////////////////////////////////////////////////////////////////////


  //Non-Dimensionalization
  if (this->refVal->mode == RefVal::NON_DIMENSIONAL) {
//  dF *= 2.0 * this->refVal->length*this->refVal->length / surface;
//  dM *= 2.0 * this->refVal->length*this->refVal->length*this->refVal->length / (surface * length);
//  dForces=dF;
//  dMoments=dM;
  std::cout<<"Non-dimensional case not really supported by sensitivity analysis"<<std::endl; exit(-1);//TOD
  }
  else {
  dF *= this->refVal->force; //redimenionalizing the derivative
  dM *= this->refVal->energy; //re-dimensionalize the energy
  F =F*dForce;//dimensionalization
  M *= dEnergy;//dimensionalization
  dForces = dF+F;//F should be zero exept for mach-sensitivity

  std::cout<<"F:  "<<F[0]<<" "<<F[1]<<" "<<F[2]<<std::endl;//TODO delete line
  std::cout<<"dF: "<<dF[0]<<" "<<dF[1]<<" "<<dF[2]<<std::endl;//TODO delete line
  dMoments = dM+M;
  }
  if(DFSPAR[1]==1.0 || DFSPAR[2]==1.0){//TODO check if this is dependent on whether it is done in degrees
    dF     *= Rad2Deg;
    dForces*= Rad2Deg;
  }






}

