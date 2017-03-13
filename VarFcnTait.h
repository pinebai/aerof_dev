#ifndef _VAR_FCN_TAIT_H
#define _VAR_FCN_TAIT_H

#include <VarFcnBase.h>

//--------------------------------------------------------------------------
// This class is the VarFcn class for the Tait EOS.
// Only elementary functions are declared and/or defined here.
// All arguments must be pertinent to only a single grid node or a single
// state, since it is assumed that the Tait EOS that must be used at this 
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
// EOS: Pressure = p + a*Density^b
//
// energy is given by
//      energy = enthalpy - Pressure/Density
//      enthalpy = cp * T
//
//--------------------------------------------------------------------------
class VarFcnTait : public VarFcnBase {

protected:
  double C_;
  double invC_;
  double a_;
  double b_;
  double p_;

  bool burnable;

  void computedVdU(double *V, double *dVdU);
  void computedUdV(double *V, double *dUdV);
  int verification(int glob, double *U, double *V);

  double* lookup_table;
  int lut_size;
  double lut_min,lut_max;

  inline double pow(double a, double b) const {

    if (!lookup_table || b != b_ || a >= lut_max)
      return ::pow(a,b);

    double lut_val = (a-lut_min)/(lut_max-lut_min)*lut_size;
    double intpart;
    double frac = modf(lut_val, &intpart);
    int ip = (int)intpart;

    return lookup_table[ip]*(1.0-frac)+lookup_table[ip+1]*frac;
  }

  void createLookupTable(int lsize) {

    lut_min = ::pow(-p_/a_,1.0/b_);
    lut_max = lut_min*1.8;
    lut_size = lsize;
    lookup_table = new double[lsize];
    for (int i = 0; i < lsize; ++i) {
      double loc = ((double)i)/(lsize-1)*(lut_max-lut_min)+lut_min;
      lookup_table[i] = ::pow(loc, b_);
    }
  }
 
public:
  // baseClass determines if VarFcnTait is used as a base class
  // for another class (like VarFcnTaitSA and VarFcnTaitKE)
  VarFcnTait(FluidModelData &data, bool baseClass = false);
  virtual ~VarFcnTait() {
  
    if (lookup_table)
      delete [] lookup_table;
  }

  //----- Transformation Operators -----//
  void conservativeToPrimitive(double *U, double *V);
  void primitiveToConservative(double *V, double *U);
  
  void extrapolatePrimitive(double un, double c, double *Vb, double *Vinter, double *V);
  void extrapolateCharacteristic(double n[3], double un, double c, double *Vb, double *dV);
  void primitiveToCharacteristicVariations(double n[3], double *V, double *dV, double *dW);
  void characteristicToPrimitiveVariations(double n[3], double *V, double *dW, double *dV);

  //----- General Functions -----//
  double getPressure(const double density) const{ return p_ + a_ * pow(density, b_); }
  double getPressure(double *V) const{ return p_ + a_*pow(V[0],b_); }

  void setPressure(const double p, double *V){V[0] = pow( (p-p_)/a_ , 1.0/b_); }
  void setPressure(double *V, double *Vorig) {V[0] = Vorig[0];}

  bool isBurnable() const { return burnable; }

  //checks that the Euler equations are still hyperbolic
  double checkPressure(double *V) const{ return getPressure(V); }
  bool checkReconstructedValues(double *V, int nodeNum, int otherNodeNum, int phi, int otherPhi, int failsafe) const{
    bool error = false;
    if(V[0] <= 0.0){
      error = true;
      if (failsafe)
        fprintf(stdout, "*** Warning: negative density (%e) for node %d after reconstruction on edge %d(%e) -> %d(%e)\n",
          V[0], nodeNum, nodeNum, double(phi), otherNodeNum, double(otherPhi));
      else
        fprintf(stderr, "*** Error: negative density (%e) for node %d after reconstruction on edge %d(%e) -> %d(%e)\n",
          V[0], nodeNum, nodeNum, double(phi), otherNodeNum, double(otherPhi));
    }
    // no check of pressure or temperature since hyperbolicity of the Euler equations with Tait EOS relies
    // only on the square of the speed of sound which is a_ * b_ * pow(V[0], b_ - 1.0)
    return error;
  }

  double computeTemperature(double *V) const{ return V[4]; }
  void getV4FromTemperature(double *V, double T) const {
    V[4] = T;
  }
  double computeRhoEnergy(double *V)   const{
    return computeRhoEpsilon(V) + 0.5 * V[0] * getVelocitySquare(V);
  }
  //computes internal energy (=rho*e-0.5*rho*u^2) knowing that h = cp * T
  double computeRhoEpsilon(double *V)  const{
    //return V[0] * C_ * V[4]; //version is epsilon = cv * T
    // In the case of a burnable fluid, override the correct rho*epsilon
    double pb = (!burnable?getPressure(V):0.0);
    return V[0] * C_ * V[4] - pb;
  }
  double computeSoundSpeed(double *V)  const{ 
    double c2 = a_ * b_ * pow(V[0], b_)/V[0];
    if (c2>0) return sqrt(c2);
    else {
      std::cout << "Negative c2 for Tait EOS: " << V[0] << "; c^2 = " << c2 << std::endl;
      return 0.0;
    }
  }
  double computeSoundSpeed(const double density, const double entropy) const{
    double c2 = a_ * b_ * pow(density, b_)/density;
    if (c2>0) return sqrt(c2);
    return 0.0;
  }
  double getPressureDerivative(double *V) const { return computeSoundSpeed(V); }
  double specificHeatCstPressure() const{ return C_; }
  double computePressureCoefficient(double *V, double pinfty, double mach, bool dimFlag) const {
    return 2.0 * (getPressure(V) - pinfty);
  }
  double computeTotalPressure(double machr, double* V) const {
    double mach = computeMachNumber(V);
    double opmach = 1.0 + 0.5*(b_-1.0)*mach*mach;
    double total_density = V[0] * pow(opmach, 1.0/(b_-1.0));
    return getPressure(total_density);
  }

  //----- Equation of State Parameters -----//
  double getCv()         const{ return C_; }
  double getAlphaWater() const{ return a_; }
  double getBetaWater()  const{ return b_; }
  double getPrefWater()  const{ return p_; }

  //----- EOS related functions -----//

};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
inline
VarFcnTait::VarFcnTait(FluidModelData &data, bool baseClass) : VarFcnBase(data) {

  if(data.fluid != FluidModelData::LIQUID){
    fprintf(stderr, "*** Error: FluidModelData is not of type Tait\n");
    exit(1);
  }

  burnable = (data.liquidModel.burnable == LiquidModelData::YES);

  type = TAIT;

  lookup_table = 0;
 
  C_    = data.liquidModel.specificHeat;
  if(C_ <= 0.0 ) fprintf(stderr, "*** Error: specific heat should not be %e in Tait EOS\n", C_);
  invC_ = 1.0/C_;
  a_     = data.liquidModel.alpha;
  b_     = data.liquidModel.beta;
  p_     = data.liquidModel.Pref;

  if(!baseClass){
    pname = new const char*[5];
    pname[0] = "density";
    pname[1] = "x-velocity";
    pname[2] = "y-velocity";
    pname[3] = "z-velocity";
    pname[4] = "temperature";
  }

  //createLookupTable(512);  
}
//------------------------------------------------------------------------------
inline
void VarFcnTait::conservativeToPrimitive(double *U, double *V){

  V[0] = U[0];

  double invRho;

  // Avoid divide by zero errors, in case this function is called with an 
  // invalid density
  if (U[0] > 0.0)
    invRho = 1.0 / U[0];
  else
    invRho = 0.0;
   
  V[1] = U[1] * invRho;
  V[2] = U[2] * invRho;
  V[3] = U[3] * invRho;
      
  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];
   
  double pb;
  // Avoid divide by zero errors, in case this function is called with an 
  // invalid density
  if (U[0] > 0.0)
    pb = (!burnable?getPressure(V[0]):0.0); 
  else
    pb = -1.0e20;

  V[4] = invRho * invC_ * (U[4] - 0.5 * U[0] * vel2 + pb);
  //note that h = cp * T and not epsilon = cv * T

}
//------------------------------------------------------------------------------
inline
int VarFcnTait::verification(int glob, double *U, double *V)
{
//verification of pressure value
//if pressure < pmin, set pressure to pmin
//and rewrite V and U!!
  double locPressure;
  if (V[0] > 0.0)
    locPressure = getPressure(V);
  else {
    return 1;
  }
  if(locPressure<pmin){
    if(verif_clipping)
      fprintf(stdout, "clip pressure[%d] in tait from %e to %e\n", glob, locPressure, pmin);
    V[0] = pow((pmin-p_)/a_, 1.0/b_);
    U[0] = V[0];
    U[1] = V[0]*V[1];
    U[2] = V[0]*V[2];
    U[3] = V[0]*V[3];
    double pb = (!burnable?getPressure(V[0]):0.0); 
    U[4] = V[0]*V[4]*C_+0.5*V[0]*(V[1]*V[1]+V[2]*V[2]+V[3]*V[3]) - pb;
    return 2;
  }
  return 0;

}
//------------------------------------------------------------------------------
inline
void VarFcnTait::primitiveToConservative(double *V, double *U) {
  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];
                                            
  U[0] = V[0];
  U[1] = V[0] * V[1];
  U[2] = V[0] * V[2];
  U[3] = V[0] * V[3];
  double pb = (!burnable?getPressure(V[0]):0.0); 
  U[4] = V[0] * C_ * V[4] + 0.5 * V[0] * vel2 - pb;

}
//------------------------------------------------------------------------------
inline
void VarFcnTait::extrapolatePrimitive(double un, double c, double *Vb,
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
void VarFcnTait::extrapolateCharacteristic(double n[3], double un, double c,
                                           double *Vb, double *dV)
{
//cf Research/latex/notes/matrices, and look at L and L^{-1}
/* routine computes boundary conditions using characteristic methods
 * and assuming that values are small perturbations of values at infinity
 * initially dV contains perturbations of primitive variables to be extrapolated
 * at return, dV contains perturbations of primitive variables at the boundary
 */

  double dVn = dV[1]*n[0]+dV[2]*n[1]+dV[3]*n[2];
  double rhooc = Vb[0]/c;
  double P = p_ +a_*pow(Vb[0],b_);
  double coeff1 = -P*dV[0]/Vb[0] + Vb[0]*C_*dV[4];

// step 1: primitive to characteristic variations
  double dW[5];
  dW[0] = coeff1*n[0] - Vb[0]*(n[2]*dV[2]-n[1]*dV[3]);
  dW[1] = coeff1*n[1] - Vb[0]*(n[0]*dV[3]-n[2]*dV[1]);
  dW[2] = coeff1*n[2] - Vb[0]*(n[1]*dV[0]-n[0]*dV[2]);
  dW[3] = 0.5*(dV[0] + rhooc*dVn);
  dW[4] = 0.5*(dV[0] - rhooc*dVn);

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
  double sum  = dW[3]+dW[4];
  double diff = dW[3]-dW[4];
  double corho = 1.0/rhooc;
  double oorho = 1.0/Vb[0];

  dV[0] = sum;
  dV[1] = oorho*(n[2]*dW[1]-n[1]*dW[2]) + corho*n[0]*diff;
  dV[2] = oorho*(n[0]*dW[2]-n[2]*dW[0]) + corho*n[1]*diff;
  dV[3] = oorho*(n[1]*dW[0]-n[0]*dW[1]) + corho*n[2]*diff;
  dV[4] = ( dW[0]*n[0]+dW[1]*n[1]+dW[2]*n[2] + P*sum*oorho)*oorho/C_;

}
//------------------------------------------------------------------------------
inline
void VarFcnTait::primitiveToCharacteristicVariations(double n[3], double *V, 
                                                     double *dV, double *dW)
{

  double Ptemp = -(p_+a_*pow(V[0],b_))/pow(V[0],2);
  double c = computeSoundSpeed(V)/V[0];
  double coeff1 = Ptemp*dV[0]+ C_*dV[4];
  double dVn = n[0]*dV[1]+n[1]*dV[2]+n[2]*dV[3];

  dW[0] = coeff1*n[0] - n[2]*dV[2] + n[1]*dV[3];
  dW[1] = coeff1*n[1] - n[0]*dV[3] + n[2]*dV[1];
  dW[2] = coeff1*n[2] - n[1]*dV[1] + n[0]*dV[2];
  dW[3] = 0.5*(c*dV[0]+dVn);
  dW[4] = 0.5*(c*dV[0]-dVn);
}
//------------------------------------------------------------------------------
inline
void VarFcnTait::characteristicToPrimitiveVariations(double n[3], double *V, 
                                                     double *dW, double *dV)
{
  double ooCv = 1.0/C_;
  double c = computeSoundSpeed(V);
  double Ptemp = (p_+a_*pow(V[0],b_))/(V[0]*C_*c);
  double rho = V[0]/c;

  dV[0] = rho*(dW[3]+dW[4]);
  dV[1] = n[2]*dV[1]-n[1]*dV[2] + n[0]*(dW[3]-dW[4]);
  dV[2] = n[0]*dV[2]-n[2]*dV[0] + n[1]*(dW[3]-dW[4]);
  dV[3] = n[1]*dV[0]-n[0]*dV[1] + n[2]*(dW[3]-dW[4]);
  dV[4] = ooCv*(n[0]*dW[0]+n[1]*dW[1]+n[2]*dW[2]) + Ptemp*(dW[3]+dW[4]);

}
//------------------------------------------------------------------------------
inline
void VarFcnTait::computedVdU(double *V, double *dVdU) {

  double invrho = 1.0 / V[0];
  double invrhoCp = invrho*invC_;
  dVdU[0]  = 1.0;
  dVdU[5]  = -invrho * V[1];
  dVdU[6]  = invrho;
  dVdU[10] = -invrho * V[2];
  dVdU[12] = invrho;
  dVdU[15] = -invrho * V[3];
  dVdU[18] = invrho;
  //dVdU[20] = invrhoCp * ( -C_*V[4]/b_ + 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) + a_*(b_-1.0)*pow(V[0], b_-1.0));
  double pb = (!burnable?getPressureDerivative(V) :0.0); 
  dVdU[20] = invrhoCp * (pb+ 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3])) - V[4]*invrho;
  dVdU[21] = -invrhoCp * V[1];
  dVdU[22] = -invrhoCp * V[2];
  dVdU[23] = -invrhoCp * V[3];
  dVdU[24] = invrhoCp;

}
//------------------------------------------------------------------------------
inline
void VarFcnTait::computedUdV(double *V, double *dUdV) {

  dUdV[0]  = 1.0;
  dUdV[5]  = V[1];
  dUdV[6]  = V[0];
  dUdV[10] = V[2];
  dUdV[12] = V[0];
  dUdV[15] = V[3];
  dUdV[18] = V[0];
  double pb = (!burnable?getPressureDerivative(V) :0.0); 
  dUdV[20] = C_ * V[4] + 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) - pb;
  dUdV[21] = V[0] * V[1];
  dUdV[22] = V[0] * V[2];
  dUdV[23] = V[0] * V[3];
  dUdV[24] = V[0] * C_;

}
//------------------------------------------------------------------------------

#endif
