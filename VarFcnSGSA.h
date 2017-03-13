#ifndef _VAR_FCN_SGSA_H
#define _VAR_FCN_SGSA_H

#include <VarFcnBase.h>

//--------------------------------------------------------------------------
// This class is the VarFcn class for the Stiffened Gas EOS in Spalat-Allmaras
// turbulent model. Only elementary functions are declared and/or defined here.
// All arguments must be pertinent to only a single grid node or a single
// state.
//
// lay-out of the base class is:
//  - 1 -  Transformation Operators
//  - 2 -  General Functions
//  - 3 -  Equations of State Parameters
//  - 4 -  EOS related functions
//
//--------------------------------------------------------------------------
//
// EOS: Pressure = (gam - 1)*Density*e - gam*Pc
// where
//   e  : internal energy per unit mass.
//   Pc : pressure constant.
//
//--------------------------------------------------------------------------

class VarFcnSGSA: public VarFcnBase {

private:
  double gam;
  double gam1;
  double invgam1;
  double Pstiff;
  double dPstiff;

  void computedVdU(double *V, double *dVdU);
  void computedUdV(double *V, double *dUdV);
  int verification(int glob, double *U, double *V);

public:
  VarFcnSGSA(FluidModelData &data);
  ~VarFcnSGSA() { delete [] pname; }

  //----- Transformation Operators -----//
  void conservativeToPrimitive(double *U, double *V);
  void primitiveToConservative(double *V, double *U);
  void conservativeToPrimitiveDerivative(double *, double *, double *, double *);
  void primitiveToConservativeDerivative(double *, double *, double *, double *);
  void multiplyBydVdU(double *, double *, double *);
  void multiplyBydVdU(double *, bcomp *, bcomp *) {fprintf(stderr,"ERROR: multiplyBydVdU needs to be implemented...\n");}
  void preMultiplyBydUdV(double *, double *, double *);
  void postMultiplyBydVdU(double *, double *, double *);
  void postMultiplyBydUdV(double *, double *, double *);
  void postMultiplyBydUdV(double *, bcomp *, bcomp *) {fprintf(stderr,"ERROR: postMultiplyBydUdV needs to be implemented...\n");}

  //----- General Functions -----//
  double checkPressure(double *V) const {
    return V[4]+Pstiff;
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

    if(V[4]+Pstiff <= 0.0){
      error = true;
      if (failsafe)
        fprintf(stdout, "*** Warning:  negative pressure (%e) for node %d (rho = %e) after reconstruction on edge %d(%e) -> %d(%e)\n",
            V[4]+Pstiff, nodeNum, V[0], nodeNum, double(phi), otherNodeNum, double(otherPhi));
      else
        fprintf(stderr, "*** Error:  negative pressure (%e) for node %d (rho = %e) after reconstruction on edge %d(%e) -> %d(%e)\n",
            V[4]+Pstiff, nodeNum, V[0], nodeNum, double(phi), otherNodeNum, double(otherPhi));
    }
    return error;
  }
  double computeTemperature(double *V) const {
    if (aerof_isnan(1.0/V[0])) {
      fprintf(stderr, "ERROR*** computeTemp\n");
      throw std::exception();
    }
    return invgam1 * (V[4]+Pstiff) / V[0];
  }
  void computeTemperatureGradient(double *V,double* Tg) const {
    if (aerof_isnan(1.0/V[0])) {
      fprintf(stderr, "ERROR*** computeTemp\n");
      throw std::exception();
    }
    Tg[0] =  -invgam1 * (V[4]+Pstiff) / (V[0]*V[0]);
    Tg[1] = Tg[2] = Tg[3] = 0.0;
    Tg[4] = invgam1 / V[0];
    Tg[5] = 0.0;
  }
  void getV4FromTemperature(double *V, double T) const {
    V[4] = T*V[0]*gam1 - Pstiff;
  }
  double computeRhoEnergy(double *V) const {
    return invgam1 * (V[4]+gam*Pstiff) + 0.5 * V[0] * (V[1]*V[1]+V[2]*V[2]+V[3]*V[3]);
  }
  double computeRhoEpsilon(double *V) const {
    return invgam1 * (V[4]+gam*Pstiff);
  }
  double computeSoundSpeed(double *V) const {
    return sqrt(gam * (V[4]+Pstiff) / V[0]);
  }
  double computeSoundSpeed(double density, double entropy) const {
    double c2 = gam * entropy*pow(density,gam-1.0);
    if(c2>0) return sqrt(c2);
    return 0.0;
  }
  double computeEntropy(double density, double pressure) const {
    return (pressure+Pstiff)/pow(density,gam);
  }
  double computeIsentropicPressure(double entropy, double density) const {
    return entropy*pow(density,gam)-Pstiff;
  }
  double computePressureCoefficient(double *V, double pinfty, double mach, bool dimFlag) const {
    if (dimFlag)
      return 2.0 * (V[4] - pinfty);
    else
      return 2.0 * (V[4] - 1.0/(gam*mach*mach));
  }
  double computeTotalPressure(double machr, double* V) const {
    double mach = computeMachNumber(V);
    double opmach = 1.0 + 0.5*gam1*mach*mach;
    return (V[4]+Pstiff)*pow(opmach, gam*invgam1) - Pstiff;
  }
  // specific heat at constant pressure is gamma for Perfect Gas
  //                                             and Stiffened Gas with h = cp * T
  double specificHeatCstPressure() const { return gam; }

  double computeDerivativeOfTemperature(double *V, double *dV) const {
    // Fix when Pstiff is non-zero.
    return ( invgam1 * dV[4] - computeTemperature(V) * dV[0] ) /V[0];
  }

  double computeDerivativeOfMachNumber(double *V, double *dV, double dMach) const
  {
    // Fix when the speed is 0
    double MyMach = computeMachNumber(V);
    if (MyMach < 100*std::numeric_limits<float>::min())
      return 0.0;
    //----
    return 1/(2.0*sqrt((V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) * V[0] / (gam * (V[4]+Pstiff)))) * ( ( (2.0*(V[1]*dV[1] + V[2]*dV[2] + V[3]*dV[3]) * V[0] + (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) * dV[0]) * (V[4]+Pstiff) - (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) * V[0] * (dV[4] + dPstiff*dMach) ) / ( (V[4]+Pstiff) * (V[4]+Pstiff) ) );
  }

  double computeDerivativeOfSoundSpeed(double *V, double *dV, double dMach) const {
    return 1.0/( 2.0*sqrt(gam * (V[4]+Pstiff) / V[0]) ) * gam * ( (dV[4]+dPstiff*dMach) * V[0] - (V[4]+Pstiff) * dV[0] ) / ( V[0] * V[0] );
  }
  double computeDerivativeOfTotalPressure(double machr, double dmachr, double* V, double* dV, double dMach) const {
    double mach = computeMachNumber(V);
    double dmach = computeDerivativeOfMachNumber(V, dV, dMach);
    double opmach = 1.0 + 0.5*gam1*mach*mach;
    double dopmach = gam1*mach*dmach;
    return dV[4]*pow(opmach, gam*invgam1) + (V[4]+Pstiff)*gam*invgam1*pow(opmach, (gam*invgam1-1))*dopmach;
  }
  void rstVar(IoData &iod) {
//    dPstiff = iod.eqs.fluidModel.gasModel.pressureConstant/iod.bc.inlet.pressure*(-2.0 / (gam * iod.bc.inlet.mach * iod.bc.inlet.mach * iod.bc.inlet.mach));
    dPstiff = iod.eqs.fluidModel.gasModel.pressureConstant*(-2.0 / (iod.bc.inlet.mach));
    rV(iod);
  }

  //----- Equation of State Parameters -----//
  double getGamma()                           const {return gam;}
  double getGamma1()                          const {return gam1;}
  double getPressureConstant()                const {return Pstiff;}
  double getPressure(double *V)               const {return V[4];}
  double getDerivativeOfPressureConstant()    const {return dPstiff;}
  double getTurbulentNuTilde(double *V)       const {return V[5];}

};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
inline
VarFcnSGSA::VarFcnSGSA(FluidModelData &data) : VarFcnBase(data) {

  if(data.fluid != FluidModelData::PERFECT_GAS && data.fluid != FluidModelData::STIFFENED_GAS){
    fprintf(stderr, "*** Error: FluidModelData is not of type GAS\n");
    exit(1);
  }

  if (data.gasModel.type == GasModelData::IDEAL)
    type = PERFECTGAS;
  else if(data.gasModel.type == GasModelData::STIFFENED)
    type = STIFFENEDGAS;

  gam = data.gasModel.specificHeatRatio;
  gam1 = gam -1.0;
  invgam1 = 1.0/gam1;
  Pstiff = data.gasModel.pressureConstant;
  dPstiff = 0.0;//iod.eqs.fluidModel.gasModel.pressureConstant/iod.bc.inlet.pressure*(-2.0 / (gam * iod.bc.inlet.mach * iod.bc.inlet.mach * iod.bc.inlet.mach));

  pname = new const char*[6];
  pname[0] = "density";
  pname[1] = "x-velocity";
  pname[2] = "y-velocity";
  pname[3] = "z-velocity";
  pname[4] = "pressure";
  pname[5] = "nut";
}
//------------------------------------------------------------------------------
inline
void VarFcnSGSA::conservativeToPrimitive(double *U, double *V)
{
  V[0] = U[0];

  double invRho = 1.0 / U[0];

  V[1] = U[1] * invRho;
  V[2] = U[2] * invRho;
  V[3] = U[3] * invRho;

  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];
  V[4] = (gam-1.0) * (U[4] - 0.5 * U[0] * vel2) - gam*Pstiff;
  V[5] = U[5] * invRho;

}
//------------------------------------------------------------------------------
inline
void VarFcnSGSA::conservativeToPrimitiveDerivative(double *U, double *dU, double *V, double *dV)
{
  dV[0] = dU[0];

  double invRho = 1.0 / V[0];

  dV[1] = ( dU[1]  - dV[0] * V[1] ) * invRho;
  dV[2] = ( dU[2]  - dV[0] * V[2] ) * invRho;
  dV[3] = ( dU[3]  - dV[0] * V[3] ) * invRho;

  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];

  double dvel2 = 2.0 * V[1] * dV[1] + 2.0 * V[2] * dV[2] + 2.0 * V[3] * dV[3];

  dV[4] = (gam-1.0) * (dU[4] - 0.5 * dU[0] * vel2  - 0.5 * U[0] * dvel2) - gam*dPstiff;

  dV[5] = ( dU[5]  - dV[0] * V[5] ) * invRho;
}
//------------------------------------------------------------------------------
inline
void VarFcnSGSA::primitiveToConservative(double *V, double *U)
{

  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];

  U[0] = V[0];
  U[1] = V[0] * V[1];
  U[2] = V[0] * V[2];
  U[3] = V[0] * V[3];
  U[4] = (V[4]+gam*Pstiff) * invgam1 + 0.5 * V[0] * vel2;
  U[5] = V[0] * V[5];

}
//------------------------------------------------------------------------------
inline
void VarFcnSGSA::primitiveToConservativeDerivative(double *V, double *dV, double *U, double *dU)
{

  double vel2 = V[1] * V[1] + V[2] * V[2] + V[3] * V[3];
  double dvel2 = 2.0 * V[1] * dV[1] + 2.0 * V[2] * dV[2] + 2.0 * V[3] * dV[3];

  dU[0] = dV[0];
  dU[1] = dV[0] * V[1] + V[0] * dV[1];
  dU[2] = dV[0] * V[2] + V[0] * dV[2];
  dU[3] = dV[0] * V[3] + V[0] * dV[3];
  dU[4] = (dV[4]+gam*dPstiff) * invgam1 + 0.5 * dV[0] * vel2 + 0.5 * V[0] * dvel2;
  dU[5] = dV[0] * V[5] + V[0] * dV[5];

}
//------------------------------------------------------------------------------
inline
int VarFcnSGSA::verification(int glob, double *U, double *V)
{
//verification of density and pressure value
//if pressure/density < pmin/rhomin, set pressure/density to pmin/rhomin
//and rewrite V and U!!
  int count = 0;

  if(V[0]<rhomin){
    if(verif_clipping)
      fprintf(stderr,"clip density[%d] in gas(SA) from %e to %e\n", glob, V[0], rhomin);
    V[0] = rhomin;
    count+= (count+1) % 2;
  }

  if(V[4]<pmin){
    if (verif_clipping)
      fprintf(stdout, "clip pressure[%d] in gas(SA) from %e to %e\n", glob, V[4], pmin);
    V[4] = pmin;
    count+=2;
  }

  if(count) //also modify U
    primitiveToConservative(V,U);

  return count;
}
//------------------------------------------------------------------------------
inline
void VarFcnSGSA::computedVdU(double *V, double *dVdU)
{
  double invrho = 1.0 / V[0];
  dVdU[0]  = 1.0;
  dVdU[6]  = -invrho * V[1];
  dVdU[7]  = invrho;
  dVdU[12] = -invrho * V[2];
  dVdU[14] = invrho;
  dVdU[18] = -invrho * V[3];
  dVdU[21] = invrho;
  dVdU[24] = gam1 * 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]);
  dVdU[25] = -gam1 * V[1];
  dVdU[26] = -gam1 * V[2];
  dVdU[27] = -gam1 * V[3];
  dVdU[28] = gam1;
  dVdU[30] = -invrho * V[5];
  dVdU[35] = invrho;
}
//------------------------------------------------------------------------------
inline
void VarFcnSGSA::computedUdV(double *V, double *dUdV)
{
  dUdV[0]  = 1.0;
  dUdV[6]  = V[1];
  dUdV[7]  = V[0];
  dUdV[12] = V[2];
  dUdV[14] = V[0];
  dUdV[18] = V[3];
  dUdV[21] = V[0];
  dUdV[24] = 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]);
  dUdV[25] = V[0] * V[1];
  dUdV[26] = V[0] * V[2];
  dUdV[27] = V[0] * V[3];
  dUdV[28] = invgam1;
  dUdV[30] = V[5];
  dUdV[35] = V[0];
}
//------------------------------------------------------------------------------
inline
void VarFcnSGSA::multiplyBydVdU(double *V, double *vec, double *res)
{
  double dVdU[36];
  computedVdU(V, dVdU);

  res[0] = dVdU[0]*vec[0];
  res[1] = dVdU[6]*vec[0]+dVdU[7]*vec[1];
  res[2] = dVdU[12]*vec[0]+dVdU[14]*vec[2];
  res[3] = dVdU[18]*vec[0]+dVdU[21]*vec[3];
  res[4] = dVdU[24]*vec[0]+dVdU[25]*vec[1]+dVdU[26]*vec[2]+dVdU[27]*vec[3]+dVdU[28]*vec[4];
  res[5] = dVdU[30]*vec[0]+dVdU[35]*vec[5];
}
//------------------------------------------------------------------------------
inline
void VarFcnSGSA::preMultiplyBydUdV(double *V, double *mat, double *res)
{
  double dUdV[36];
  computedUdV(V, dUdV);
  res[0] = dUdV[0]*mat[0];
  res[1] = dUdV[0]*mat[1];
  res[2] = dUdV[0]*mat[2];
  res[3] = dUdV[0]*mat[3];
  res[4] = dUdV[0]*mat[4];
  res[5] = dUdV[0]*mat[5];
  res[6] = dUdV[6]*mat[0]+dUdV[7]*mat[6];
  res[7] = dUdV[6]*mat[1]+dUdV[7]*mat[7];
  res[8] = dUdV[6]*mat[2]+dUdV[7]*mat[8];
  res[9] = dUdV[6]*mat[3]+dUdV[7]*mat[9];
  res[10] = dUdV[6]*mat[4]+dUdV[7]*mat[10];
  res[11] = dUdV[6]*mat[5]+dUdV[7]*mat[11];
  res[12] = dUdV[12]*mat[0]+dUdV[14]*mat[12];
  res[13] = dUdV[12]*mat[1]+dUdV[14]*mat[13];
  res[14] = dUdV[12]*mat[2]+dUdV[14]*mat[14];
  res[15] = dUdV[12]*mat[3]+dUdV[14]*mat[15];
  res[16] = dUdV[12]*mat[4]+dUdV[14]*mat[16];
  res[17] = dUdV[12]*mat[5]+dUdV[14]*mat[17];
  res[18] = dUdV[18]*mat[0]+dUdV[21]*mat[18];
  res[19] = dUdV[18]*mat[1]+dUdV[21]*mat[19];
  res[20] = dUdV[18]*mat[2]+dUdV[21]*mat[20];
  res[21] = dUdV[18]*mat[3]+dUdV[21]*mat[21];
  res[22] = dUdV[18]*mat[4]+dUdV[21]*mat[22];
  res[23] = dUdV[18]*mat[5]+dUdV[21]*mat[23];
  res[24] = dUdV[24]*mat[0]+dUdV[25]*mat[6]+dUdV[26]*mat[12]+dUdV[27]*mat[18]+dUdV[28]*mat[24];
  res[25] = dUdV[24]*mat[1]+dUdV[25]*mat[7]+dUdV[26]*mat[13]+dUdV[27]*mat[19]+dUdV[28]*mat[25];
  res[26] = dUdV[24]*mat[2]+dUdV[25]*mat[8]+dUdV[26]*mat[14]+dUdV[27]*mat[20]+dUdV[28]*mat[26];
  res[27] = dUdV[24]*mat[3]+dUdV[25]*mat[9]+dUdV[26]*mat[15]+dUdV[27]*mat[21]+dUdV[28]*mat[27];
  res[28] = dUdV[24]*mat[4]+dUdV[25]*mat[10]+dUdV[26]*mat[16]+dUdV[27]*mat[22]+dUdV[28]*mat[28];
  res[29] = dUdV[24]*mat[5]+dUdV[25]*mat[11]+dUdV[26]*mat[17]+dUdV[27]*mat[23]+dUdV[28]*mat[29];
  res[30] = dUdV[30]*mat[0]+dUdV[35]*mat[30];
  res[31] = dUdV[30]*mat[1]+dUdV[35]*mat[31];
  res[32] = dUdV[30]*mat[2]+dUdV[35]*mat[32];
  res[33] = dUdV[30]*mat[3]+dUdV[35]*mat[33];
  res[34] = dUdV[30]*mat[4]+dUdV[35]*mat[34];
  res[35] = dUdV[30]*mat[5]+dUdV[35]*mat[35];
}
//------------------------------------------------------------------------------
inline
void VarFcnSGSA::postMultiplyBydVdU(double *V, double *mat, double *res)
{
  double dVdU[36];
  computedVdU(V, dVdU);

  res[0] = mat[0]*dVdU[0]+mat[1]*dVdU[6]+mat[2]*dVdU[12]+
    mat[3]*dVdU[18]+mat[4]*dVdU[24]+mat[5]*dVdU[30];
  res[1] = mat[1]*dVdU[7]+mat[4]*dVdU[25];
  res[2] = mat[2]*dVdU[14]+mat[4]*dVdU[26];
  res[3] = mat[3]*dVdU[21]+mat[4]*dVdU[27];
  res[4] = mat[4]*dVdU[28];
  res[5] = mat[5]*dVdU[35];
  res[6] = mat[6]*dVdU[0]+mat[7]*dVdU[6]+mat[8]*dVdU[12]+
    mat[9]*dVdU[18]+mat[10]*dVdU[24]+mat[11]*dVdU[30];
  res[7] = mat[7]*dVdU[7]+mat[10]*dVdU[25];
  res[8] = mat[8]*dVdU[14]+mat[10]*dVdU[26];
  res[9] = mat[9]*dVdU[21]+mat[10]*dVdU[27];
  res[10] = mat[10]*dVdU[28];
  res[11] = mat[11]*dVdU[35];
  res[12] = mat[12]*dVdU[0]+mat[13]*dVdU[6]+mat[14]*dVdU[12]+
    mat[15]*dVdU[18]+mat[16]*dVdU[24]+mat[17]*dVdU[30];
  res[13] = mat[13]*dVdU[7]+mat[16]*dVdU[25];
  res[14] = mat[14]*dVdU[14]+mat[16]*dVdU[26];
  res[15] = mat[15]*dVdU[21]+mat[16]*dVdU[27];
  res[16] = mat[16]*dVdU[28];
  res[17] = mat[17]*dVdU[35];
  res[18] = mat[18]*dVdU[0]+mat[19]*dVdU[6]+mat[20]*dVdU[12]+
    mat[21]*dVdU[18]+mat[22]*dVdU[24]+mat[23]*dVdU[30];
  res[19] = mat[19]*dVdU[7]+mat[22]*dVdU[25];
  res[20] = mat[20]*dVdU[14]+mat[22]*dVdU[26];
  res[21] = mat[21]*dVdU[21]+mat[22]*dVdU[27];
  res[22] = mat[22]*dVdU[28];
  res[23] = mat[23]*dVdU[35];
  res[24] = mat[24]*dVdU[0]+mat[25]*dVdU[6]+mat[26]*dVdU[12]+
    mat[27]*dVdU[18]+mat[28]*dVdU[24]+mat[29]*dVdU[30];
  res[25] = mat[25]*dVdU[7]+mat[28]*dVdU[25];
  res[26] = mat[26]*dVdU[14]+mat[28]*dVdU[26];
  res[27] = mat[27]*dVdU[21]+mat[28]*dVdU[27];
  res[28] = mat[28]*dVdU[28];
  res[29] = mat[29]*dVdU[35];
  res[30] = mat[30]*dVdU[0]+mat[31]*dVdU[6]+mat[32]*dVdU[12]+
    mat[33]*dVdU[18]+mat[34]*dVdU[24]+mat[35]*dVdU[30];
  res[31] = mat[31]*dVdU[7]+mat[34]*dVdU[25];
  res[32] = mat[32]*dVdU[14]+mat[34]*dVdU[26];
  res[33] = mat[33]*dVdU[21]+mat[34]*dVdU[27];
  res[34] = mat[34]*dVdU[28];
  res[35] = mat[35]*dVdU[35];
}
//------------------------------------------------------------------------------
inline
void VarFcnSGSA::postMultiplyBydUdV(double *V, double *mat, double *res)
{
  double dUdV[36];
  computedUdV(V, dUdV);

  res[0] = mat[0]*dUdV[0]+mat[1]*dUdV[6]+mat[2]*dUdV[12]+
    mat[3]*dUdV[18]+mat[4]*dUdV[24]+mat[5]*dUdV[30];
  res[1] = mat[1]*dUdV[7]+mat[4]*dUdV[25];
  res[2] = mat[2]*dUdV[14]+mat[4]*dUdV[26];
  res[3] = mat[3]*dUdV[21]+mat[4]*dUdV[27];
  res[4] = mat[4]*dUdV[28];
  res[5] = mat[5]*dUdV[35];
  res[6] = mat[6]*dUdV[0]+mat[7]*dUdV[6]+mat[8]*dUdV[12]+
    mat[9]*dUdV[18]+mat[10]*dUdV[24]+mat[11]*dUdV[30];
  res[7] = mat[7]*dUdV[7]+mat[10]*dUdV[25];
  res[8] = mat[8]*dUdV[14]+mat[10]*dUdV[26];
  res[9] = mat[9]*dUdV[21]+mat[10]*dUdV[27];
  res[10] = mat[10]*dUdV[28];
  res[11] = mat[11]*dUdV[35];
  res[12] = mat[12]*dUdV[0]+mat[13]*dUdV[6]+mat[14]*dUdV[12]+
    mat[15]*dUdV[18]+mat[16]*dUdV[24]+mat[17]*dUdV[30];
  res[13] = mat[13]*dUdV[7]+mat[16]*dUdV[25];
  res[14] = mat[14]*dUdV[14]+mat[16]*dUdV[26];
  res[15] = mat[15]*dUdV[21]+mat[16]*dUdV[27];
  res[16] = mat[16]*dUdV[28];
  res[17] = mat[17]*dUdV[35];
  res[18] = mat[18]*dUdV[0]+mat[19]*dUdV[6]+mat[20]*dUdV[12]+
    mat[21]*dUdV[18]+mat[22]*dUdV[24]+mat[23]*dUdV[30];
  res[19] = mat[19]*dUdV[7]+mat[22]*dUdV[25];
  res[20] = mat[20]*dUdV[14]+mat[22]*dUdV[26];
  res[21] = mat[21]*dUdV[21]+mat[22]*dUdV[27];
  res[22] = mat[22]*dUdV[28];
  res[23] = mat[23]*dUdV[35];
  res[24] = mat[24]*dUdV[0]+mat[25]*dUdV[6]+mat[26]*dUdV[12]+
    mat[27]*dUdV[18]+mat[28]*dUdV[24]+mat[29]*dUdV[30];
  res[25] = mat[25]*dUdV[7]+mat[28]*dUdV[25];
  res[26] = mat[26]*dUdV[14]+mat[28]*dUdV[26];
  res[27] = mat[27]*dUdV[21]+mat[28]*dUdV[27];
  res[28] = mat[28]*dUdV[28];
  res[29] = mat[29]*dUdV[35];
  res[30] = mat[30]*dUdV[0]+mat[31]*dUdV[6]+mat[32]*dUdV[12]+
    mat[33]*dUdV[18]+mat[34]*dUdV[24]+mat[35]*dUdV[30];
  res[31] = mat[31]*dUdV[7]+mat[34]*dUdV[25];
  res[32] = mat[32]*dUdV[14]+mat[34]*dUdV[26];
  res[33] = mat[33]*dUdV[21]+mat[34]*dUdV[27];
  res[34] = mat[34]*dUdV[28];
  res[35] = mat[35]*dUdV[35];
}
//------------------------------------------------------------------------------

#endif
