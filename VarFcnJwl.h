#ifndef _VAR_FCN_JWL_H
#define _VAR_FCN_JWL_H

#include <VarFcnBase.h>

//--------------------------------------------------------------------------
// This class is the VarFcn class for the JWL EOS.
// Only elementary functions are declared and/or defined here.
// All arguments must be pertinent to only a single grid node or a single
// state, since it is assumed that the JWL EOS that must be used at this 
// point.
//
// lay-out of the base class is:
//  - 1 -  Transformation Operators
//  - 2 -  General Functions
//  - 3 -  Equations of State Parameters
//  - 4 -  EOS related functions
//
//--------------------------------------------------------------------------
//
// EOS: Pressure = omega*Density*Epsilon + f(Density)
// where: 
//   f(x) = A1*(1.0-omega*x/R1r)*exp(-R1r/x) 
//        + A2*(1.0-omega*x/R2r)*exp(-R2r/x)
//
//--------------------------------------------------------------------------
class VarFcnJwl : public VarFcnBase {

private:
  double omega;
  double invomega;
  double omegap1;
  double invomegap1;
  double A1, A2, R1, R2, rhoref;
  double R1r, R2r;

  void computedVdU(double *V, double *dVdU);
  void computedUdV(double *V, double *dUdV);
  int verification(int glob, double *U, double *V);

public:
  VarFcnJwl(FluidModelData &data);
  ~VarFcnJwl() { delete [] pname;}

  //----- Transformation Operators -----//
  void conservativeToPrimitive(double *U, double *V);
  void primitiveToConservative(double *V, double *U);

  void extrapolatePrimitive(double un, double c, double *Vb, double *Vinter, double *V);
  void extrapolateCharacteristic(double n[3], double un, double c, double *Vb, double *dV);
  void primitiveToCharacteristicVariations(double n[3], double *V, double *dV, double *dW);
  void characteristicToPrimitiveVariations(double n[3], double *V, double *dW, double *dV);

  //----- General Functions -----//
  //checks that the Euler equations are still hyperbolic
  double checkPressure(double *V) const{
    return V[4] - (computeFrho(V) - computeFrhop(V)*V[0])*invomegap1;
  }
  bool checkReconstructedValues(double *V, int nodeNum, int otherNodeNum, int phi, int otherPhi, int failsafe) const{
    bool error = false;
    if(V[0] <= 0.0){
      error = true;
      if (failsafe)
        fprintf(stdout, "*** Warning:  negative density (%e) for node %d after reconstruction on edge %d(%e) -> %d(%e)\n",
          V[0], nodeNum, nodeNum, double(phi), otherNodeNum, double(otherPhi));
      else
        fprintf(stderr, "*** Error:  negative density (%e) for node %d after reconstruction on edge %d(%e) -> %d(%e)\n",
          V[0], nodeNum, nodeNum, double(phi), otherNodeNum, double(otherPhi));
    }
    double pressureCheck = checkPressure(V);
    if(pressureCheck <= 0.0){
      error = true;
      if (failsafe)
        fprintf(stdout, "*** Warning:  negative pressure (%e) for node %d (rho = %e) after reconstruction on edge %d(%e) -> %d(%e)\n",
            pressureCheck, nodeNum, V[0], nodeNum, double(phi), otherNodeNum, double(otherPhi));
      else
        fprintf(stderr, "*** Error:  negative pressure (%e) for node %d (rho = %e) after reconstruction on edge %d(%e) -> %d(%e)\n",
            pressureCheck, nodeNum, V[0], nodeNum, double(phi), otherNodeNum, double(otherPhi));
    }
    return error;
  }

  double computeTemperature(double *V) const{
    return invomega * (V[4]-computeFrho(V)) / V[0]; 
  }
  void getV4FromTemperature(double *V, double T) const {
    V[4] = T*omega*V[0] + computeFrho(V);
  }
  double computeRhoEnergy(double *V)   const{
    return invomega * (V[4]-computeFrho(V)) + 0.5 * V[0] * getVelocitySquare(V);
  }
  //computes internal energy (=rho*e-0.5*rho*u^2)
  double computeRhoEpsilon(double *V)  const{ 
    return invomega * (V[4]-computeFrho(V));
  }
  double computeSoundSpeed(double *V)  const{ 
    return sqrt((omegap1*V[4] - computeFrho(V) + V[0]*computeFrhop(V))/V[0]);
  }
  double computeSoundSpeed(const double density, const double entropy) const{
    double pressure = entropy*pow(density,omegap1)+A1*exp(-R1r/density)+A2*exp(-R2r/density);
    double c2 = (omegap1*pressure - computeFrho(density) + density*computeFrhop(density))/density;
    if(c2>0) return sqrt(c2);
    return 0.0; 
  }
  double computeEntropy(const double density, const double pressure) const{
    return (pressure - A1*exp(-R1r/density) - A2*exp(-R2r/density))/pow(density,omegap1);
  }
  double computeIsentropicPressure(const double entropy, const double density) const{
    return entropy*pow(density,omegap1)+A1*exp(-R1r/density)+A2*exp(-R2r/density);
  }

  //----- Equation of State Parameters -----//
  double getOmega()      const{ return omega; }
  double getOmegap1()    const{ return omegap1; }
  double getInvOmega()   const{ return invomega; }
  double getInvOmegap1() const{ return invomegap1; }
  double getA1()         const{ return A1; }
  double getA2()         const{ return A2; }
  double getR1()         const{ return R1; }
  double getR2()         const{ return R2; }
  double getRhoref()     const{ return rhoref; }
  double getR1r()        const{ return R1r; }
  double getR2r()        const{ return R2r; }

  //----- EOS related functions -----//
  double computeExponentials(const double density) const{
    return A1*exp(-R1r/density) + A2*exp(-R2r/density);
  }
  double computeDerivativeOfExponentials(const double density) const{
    double invrho2 = 1.0/(density*density);
    return R1r*invrho2*A1*exp(-R1r/density) + R2r*invrho2*A2*exp(-R2r/density);
  }
  double computeExponentials2(const double density) const{
    return (A1*R1r*exp(-R1r/density) + A2*R2r*exp(-R2r/density))/(density*density);
  }
  double computeDerivativeOfExponentials2(const double density) const{
    double invrho2 = 1.0/(density*density);
    return (R1r*A1*(R1r-2.0*density)*exp(-R1r/density) + 
            R2r*A2*(R2r-2.0*density)*exp(-R2r/density)  )*invrho2*invrho2;
  }
  double computeFrho(double *V) const{
    return A1*(1.0-omega*V[0]/R1r)*exp(-R1r/V[0])
         + A2*(1.0-omega*V[0]/R2r)*exp(-R2r/V[0]); 
  }
  double computeFrho(const double rho) const{
    return A1*(1.0-omega*rho/R1r)*exp(-R1r/rho)
         + A2*(1.0-omega*rho/R2r)*exp(-R2r/rho);
  }
  double computeFrhop(double *V) const{
    double rho2 = V[0]*V[0];
    return A1*(-omega/R1r+(1.0-omega*V[0]/R1r)*R1r/rho2)*exp(-R1r/V[0])
         + A2*(-omega/R2r+(1.0-omega*V[0]/R2r)*R2r/rho2)*exp(-R2r/V[0]);
  }
  double computeFrhop(const double rho) const{
    double rho2 = rho*rho;
    return A1*(-omega/R1r+(1.0-omega*rho/R1r)*R1r/rho2)*exp(-R1r/rho)
         + A2*(-omega/R2r+(1.0-omega*rho/R2r)*R2r/rho2)*exp(-R2r/rho);
  }

};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
inline
VarFcnJwl::VarFcnJwl(FluidModelData &data) : VarFcnBase(data) {

  if(data.fluid != FluidModelData::JWL){
    fprintf(stderr, "*** Error: FluidModelData is not of type JWL\n");
    exit(1);
  }

  type = JWL;
  
  omega        = data.jwlModel.omega;
  omegap1      = data.jwlModel.omega + 1.0;
  invomega     = 1.0/data.jwlModel.omega;
  invomegap1   = 1.0/(data.jwlModel.omega+1.0);

  A1           = data.jwlModel.A1;
  A2           = data.jwlModel.A2;
  R1           = data.jwlModel.R1;
  R2           = data.jwlModel.R2;
  rhoref       = data.jwlModel.rhoref;
  R1r          = R1*rhoref;
  R2r          = R2*rhoref;

  pname = new const char*[5];
  pname[0] = "density";
  pname[1] = "x-velocity";
  pname[2] = "y-velocity";
  pname[3] = "z-velocity";
  pname[4] = "pressure";
}
//------------------------------------------------------------------------------
inline
void VarFcnJwl::conservativeToPrimitive(double *U, double *V){

  V[0] = U[0];

  double invRho = 1.0 / U[0];
   
  V[1] = U[1] * invRho;
  V[2] = U[2] * invRho;
  V[3] = U[3] * invRho;
      
  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];
    
  V[4] = omega*(U[4] - 0.5 * U[0] * vel2) + computeFrho(max(V[0],rhomin));
}
//------------------------------------------------------------------------------
inline
int VarFcnJwl::verification(int glob, double *U, double *V)
{
//verification of density and pressure value
//if pressure/density < pmin/rhomin, set pressure/density to pmin/rhomin
//and rewrite V and U!!
  int count = 0;

  if(V[0]<rhomin){
    if(verif_clipping)
      fprintf(stderr,"clip density[%d] in JWL from %e to %e\n", glob, V[0], rhomin);
    V[0] = rhomin;
    count += (count+1) % 2;
  }

  if(V[4]<pmin){
    if (verif_clipping)
      fprintf(stdout, "clip pressure[%d] in JWL from %e to %e\n", glob, V[4], pmin);
    V[4] = pmin;
    count += 2;
  }

  if(count) //also modify U
    primitiveToConservative(V,U);

  return count;
}
//------------------------------------------------------------------------------
inline
void VarFcnJwl::primitiveToConservative(double *V, double *U) {

  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];
                                            
  U[0] = V[0];
  U[1] = V[0] * V[1];
  U[2] = V[0] * V[2];
  U[3] = V[0] * V[3];
  U[4] = invomega*(V[4] - computeFrho(V[0])) + 0.5 * V[0] * vel2;

}
//------------------------------------------------------------------------------
inline
void VarFcnJwl::extrapolatePrimitive(double un, double c, double *Vb,
                                     double *Vinter, double *V)
{
  if (un == 0.0){
    V[0] = Vinter[0];
    V[1] = Vb[1];
    V[2] = Vb[2];
    V[3] = Vb[3];
    V[4] = Vb[4];
  }else{
    if (un<0.0){             // INLET
      if (-un-c > 0.0){      //    SUPERSONIC
        V[0] = Vb[0];
        V[1] = Vb[1];
        V[2] = Vb[2];
        V[3] = Vb[3];
        V[4] = Vb[4];
      }else{                 //    SUBSONIC
        V[0] = Vinter[0];
        V[1] = Vb[1];
        V[2] = Vb[2];
        V[3] = Vb[3];
        V[4] = Vb[4];
      }
    }else{                   // OUTLET
      if (un-c > 0.0){       //    SUPERSONIC
        V[0] = Vinter[0];
        V[1] = Vinter[1];
        V[2] = Vinter[2];
        V[3] = Vinter[3];
        V[4] = Vinter[4];
      }else{                //     SUBSONIC
        V[0] = Vb[0];
        V[1] = Vinter[1];
        V[2] = Vinter[2];
        V[3] = Vinter[3];
        V[4] = Vinter[4];
      }
    }
  }
}
//------------------------------------------------------------------------------
inline
void VarFcnJwl::extrapolateCharacteristic(double n[3], double un, double c,
                                          double *Vb, double *dV)
{
//cf Research/latex/notes/matrices, and look at L and L^{-1}
/* routine computes boundary conditions using characteristic methods
 * and assuming that values are small perturbations of values at infinity
 * initially dV contains perturbations of primitive variables to be extrapolated
 * at return, dV contains perturbations of primitive variables at the boundary
 */

  double dVn = dV[1]*n[0]+dV[2]*n[1]+dV[3]*n[2];
  double ooc2 = 1.0/(c*c);
  double oorhoc = sqrt(ooc2)/Vb[0];
  double coeff1 = dV[0] - ooc2*dV[4];

// step 1: primitive to characteristic variations
  double dW[5];
  dW[0] = n[0]*coeff1 + n[2]*dV[2] - n[1]*dV[3];
  dW[1] = n[1]*coeff1 + n[0]*dV[3] - n[2]*dV[1];
  dW[2] = n[2]*coeff1 + n[1]*dV[1] - n[0]*dV[2];
  dW[3] = dVn + oorhoc*dV[4];
  dW[4] =-dVn + oorhoc*dV[4];

// step 2: choose variations to be extrapolated
//         if incoming characteristic, then the
//         characteristic is set to 0.0
//         else, there is nothing to do
  if (un == 0.0){
    dW[0] = 0.0;
    dW[1] = 0.0;
    dW[2] = 0.0;
    dW[3] = 0.0;
  }else{
    if (un<0.0){             // INLET
      if (-un-c > 0.0){      //    SUPERSONIC
        dW[0] = 0.0;
        dW[1] = 0.0;
        dW[2] = 0.0;
        dW[3] = 0.0;
        dW[4] = 0.0;
      }else{                 //    SUBSONIC
        dW[0] = 0.0;
        dW[1] = 0.0;
        dW[2] = 0.0;
        dW[3] = 0.0;
      }
    }else{                   // OUTLET
      if (un-c > 0.0){       //    SUPERSONIC
      }else{                //     SUBSONIC
        dW[4] = 0.0;
      }
    }
  }

// step 3: characteristic to primitive variations
  double sum = dW[3]+dW[4];
  double diff = dW[3]-dW[4];

  dV[0] = dW[0]*n[0]+dW[1]*n[1]+dW[2]*n[2] + 0.5*Vb[0]*sum/c;
  dV[1] = n[1]*dW[2]-n[2]*dW[1] + 0.5*n[0]*diff;
  dV[2] = n[2]*dW[0]-n[0]*dW[2] + 0.5*n[1]*diff;
  dV[3] = n[0]*dW[1]-n[1]*dW[0] + 0.5*n[2]*diff;
  dV[4] = 0.5*Vb[0]*c*sum;

}
//---------------------------------------------------------------------------
inline
void VarFcnJwl::primitiveToCharacteristicVariations(double n[3], double *V, 
                                                    double *dV, double *dW)
{
  double dVn = dV[1]*n[0]+dV[2]*n[1]+dV[3]*n[2];
  double c = computeSoundSpeed(V);
  double ooc2 = 1.0/(c*c);
  double oorhoc = sqrt(ooc2)/V[0];
  double coeff1 = dV[0] - ooc2*dV[4];

  dW[0] = n[0]*coeff1 + n[2]*dV[2] - n[1]*dV[3];
  dW[1] = n[1]*coeff1 + n[0]*dV[3] - n[2]*dV[1];
  dW[2] = n[2]*coeff1 + n[1]*dV[1] - n[0]*dV[2];
  dW[3] = dVn + oorhoc*dV[4];
  dW[4] =-dVn + oorhoc*dV[4];
}
//---------------------------------------------------------------------------
inline
void VarFcnJwl::characteristicToPrimitiveVariations(double n[3], double *V, 
                                                    double *dW, double *dV)
{
  double sum = dW[3]+dW[4];
  double diff = dW[3]-dW[4];
  double c = computeSoundSpeed(V);

  dV[0] = dV[1]*n[0]+dV[2]*n[1]+dV[3]*n[2] + 0.5*V[0]*sum/c;
  dV[1] = n[1]*dV[2]-n[2]*dV[1] + 0.5*n[0]*diff;
  dV[2] = n[2]*dV[0]-n[0]*dV[2] + 0.5*n[1]*diff;
  dV[3] = n[0]*dV[1]-n[1]*dV[0] + 0.5*n[2]*diff;
  dV[4] = 0.5*V[0]*c*sum;
}
//---------------------------------------------------------------------------
inline
void VarFcnJwl::computedVdU(double *V, double *dVdU) {

  double invrho = 1.0 / V[0];
  dVdU[0]  = 1.0;
  dVdU[5]  = -invrho * V[1];
  dVdU[6]  = invrho;
  dVdU[10] = -invrho * V[2];
  dVdU[12] = invrho;
  dVdU[15] = -invrho * V[3];
  dVdU[18] = invrho;
  dVdU[20] = omega * 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) + computeFrhop(V);
  dVdU[21] = -omega * V[1];
  dVdU[22] = -omega * V[2];
  dVdU[23] = -omega * V[3];
  dVdU[24] = omega;

}
//------------------------------------------------------------------------------
inline
void VarFcnJwl::computedUdV(double *V, double *dUdV) {

  dUdV[0]  = 1.0;
  dUdV[5]  = V[1];
  dUdV[6]  = V[0];
  dUdV[10] = V[2];
  dUdV[12] = V[0];
  dUdV[15] = V[3];
  dUdV[18] = V[0];
  dUdV[20] = 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) - invomega*computeFrhop(V);
  dUdV[21] = V[0] * V[1];
  dUdV[22] = V[0] * V[2];
  dUdV[23] = V[0] * V[3];
  dUdV[24] = invomega;

}
//------------------------------------------------------------------------------

#endif
