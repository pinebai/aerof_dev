#ifndef _VAR_FCN_TAIT_KE_H
#define _VAR_FCN_TAIT_KE_H

#include "VarFcnBase.h"
#include "VarFcnTait.h"

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
//--------------------------------------------------------------------------
class VarFcnTaitKE : public VarFcnTait {

protected:

  void computedVdU(double *V, double *dVdU);
  void computedUdV(double *V, double *dUdV);

public:
  VarFcnTaitKE(FluidModelData &data);
  ~VarFcnTaitKE() {}

  //----- Transformation Operators -----//
  void conservativeToPrimitive(double *U, double *V);
  void primitiveToConservative(double *V, double *U);
  
  void multiplyBydVdU(double *, double *, double *);
  void preMultiplyBydUdV(double *, double *, double *);
  void postMultiplyBydVdU(double *, double *, double *);
  void postMultiplyBydUdV(double *, double *, double *);

};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
inline
VarFcnTaitKE::VarFcnTaitKE(FluidModelData &data) : VarFcnTait(data, true) {

  if(data.fluid != FluidModelData::LIQUID){
    fprintf(stderr, "*** Error: FluidModelData is not of type Tait in VarFcnTaitKE.\n");
    fflush(stderr);
    exit(1);
  }

  pname = new const char*[7];
  pname[0] = "density";
  pname[1] = "x-velocity";
  pname[2] = "y-velocity";
  pname[3] = "z-velocity";
  pname[4] = "temperature";
  pname[5] = "k";
  pname[6] = "epsilon";
}
//------------------------------------------------------------------------------
inline
void VarFcnTaitKE::conservativeToPrimitive(double *U, double *V){

  VarFcnTait::conservativeToPrimitive(U,V);
  V[5] = U[5] / U[0];
  V[6] = U[6] / U[0];

}
//------------------------------------------------------------------------------
inline
void VarFcnTaitKE::primitiveToConservative(double *V, double *U) {

  VarFcnTait::primitiveToConservative(V,U);
  U[5] = V[0] * V[5];
  U[6] = V[0] * V[6];

}
//------------------------------------------------------------------------------
inline
void VarFcnTaitKE::computedVdU(double *V, double *dVdU) {

  double invrho = 1.0 / V[0];
  double invrhoCp = invrho*invC_;
  dVdU[0]  = 1.0;
  dVdU[7]  = -invrho * V[1];
  dVdU[8]  = invrho;
  dVdU[14] = -invrho * V[2];
  dVdU[16] = invrho;
  dVdU[21] = -invrho * V[3];
  dVdU[24] = invrho;
  dVdU[28] = invrhoCp * (getPressureDerivative(V) + 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3])) - V[4]*invrho;
  dVdU[29] = -invrhoCp * V[1];
  dVdU[30] = -invrhoCp * V[2];
  dVdU[31] = -invrhoCp * V[3];
  dVdU[32] = invrhoCp;
  dVdU[35] = -invrho * V[5];
  dVdU[40] = invrho;
  dVdU[42] = -invrho * V[6];
  dVdU[48] = invrho;

}
//------------------------------------------------------------------------------
inline
void VarFcnTaitKE::computedUdV(double *V, double *dUdV) {

  dUdV[0]  = 1.0;
  dUdV[7]  = V[1];
  dUdV[8]  = V[0];
  dUdV[14] = V[2];
  dUdV[16] = V[0];
  dUdV[21] = V[3];
  dUdV[24] = V[0];
  dUdV[28] = C_ * V[4] + 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) - getPressureDerivative(V);
  dUdV[29] = V[0] * V[1];
  dUdV[30] = V[0] * V[2];
  dUdV[31] = V[0] * V[3];
  dUdV[32] = V[0] * C_;
  dUdV[35] = V[5];
  dUdV[40] = V[0];
  dUdV[42] = V[6];
  dUdV[48] = V[0];
}
//------------------------------------------------------------------------------
inline
void VarFcnTaitKE::multiplyBydVdU(double *V, double *vec, double *res)
{

  double dVdU[49];
  computedVdU(V, dVdU);

  res[0] = dVdU[0]*vec[0];
  res[1] = dVdU[7]*vec[0]+dVdU[8]*vec[1];
  res[2] = dVdU[14]*vec[0]+dVdU[16]*vec[2];
  res[3] = dVdU[21]*vec[0]+dVdU[24]*vec[3];
  res[4] = dVdU[28]*vec[0]+dVdU[29]*vec[1]+dVdU[30]*vec[2]+dVdU[31]*vec[3]+dVdU[32]*vec[4];
  res[5] = dVdU[35]*vec[0]+dVdU[40]*vec[5];
  res[6] = dVdU[42]*vec[0]+dVdU[48]*vec[6];

}
//------------------------------------------------------------------------------
inline
void VarFcnTaitKE::preMultiplyBydUdV(double *V, double *mat, double *res)
{
  double dUdV[49];
  computedUdV(V, dUdV);

  res[0] = dUdV[0]*mat[0];
  res[1] = dUdV[0]*mat[1];
  res[2] = dUdV[0]*mat[2];
  res[3] = dUdV[0]*mat[3];
  res[4] = dUdV[0]*mat[4];
  res[5] = dUdV[0]*mat[5];
  res[6] = dUdV[0]*mat[6];
  res[7] = dUdV[7]*mat[0]+dUdV[8]*mat[7];
  res[8] = dUdV[7]*mat[1]+dUdV[8]*mat[8];
  res[9] = dUdV[7]*mat[2]+dUdV[8]*mat[9];
  res[10] = dUdV[7]*mat[3]+dUdV[8]*mat[10];
  res[11] = dUdV[7]*mat[4]+dUdV[8]*mat[11];
  res[12] = dUdV[7]*mat[5]+dUdV[8]*mat[12];
  res[13] = dUdV[7]*mat[6]+dUdV[8]*mat[13];
  res[14] = dUdV[14]*mat[0]+dUdV[16]*mat[14];
  res[15] = dUdV[14]*mat[1]+dUdV[16]*mat[15];
  res[16] = dUdV[14]*mat[2]+dUdV[16]*mat[16];
  res[17] = dUdV[14]*mat[3]+dUdV[16]*mat[17];
  res[18] = dUdV[14]*mat[4]+dUdV[16]*mat[18];
  res[19] = dUdV[14]*mat[5]+dUdV[16]*mat[19];
  res[20] = dUdV[14]*mat[6]+dUdV[16]*mat[20];
  res[21] = dUdV[21]*mat[0]+dUdV[24]*mat[21];
  res[22] = dUdV[21]*mat[1]+dUdV[24]*mat[22];
  res[23] = dUdV[21]*mat[2]+dUdV[24]*mat[23];
  res[24] = dUdV[21]*mat[3]+dUdV[24]*mat[24];
  res[25] = dUdV[21]*mat[4]+dUdV[24]*mat[25];
  res[26] = dUdV[21]*mat[5]+dUdV[24]*mat[26];
  res[27] = dUdV[21]*mat[6]+dUdV[24]*mat[27];
  res[28] = dUdV[28]*mat[0]+dUdV[29]*mat[7]+dUdV[30]*mat[14]+dUdV[31]*mat[21]+dUdV[32]*mat[28];
  res[29] = dUdV[28]*mat[1]+dUdV[29]*mat[8]+dUdV[30]*mat[15]+dUdV[31]*mat[22]+dUdV[32]*mat[29];
  res[30] = dUdV[28]*mat[2]+dUdV[29]*mat[9]+dUdV[30]*mat[16]+dUdV[31]*mat[23]+dUdV[32]*mat[30];
  res[31] = dUdV[28]*mat[3]+dUdV[29]*mat[10]+dUdV[30]*mat[17]+dUdV[31]*mat[24]+dUdV[32]*mat[31];
  res[32] = dUdV[28]*mat[4]+dUdV[29]*mat[11]+dUdV[30]*mat[18]+dUdV[31]*mat[25]+dUdV[32]*mat[32];
  res[33] = dUdV[28]*mat[5]+dUdV[29]*mat[12]+dUdV[30]*mat[19]+dUdV[31]*mat[26]+dUdV[32]*mat[33];
  res[34] = dUdV[28]*mat[6]+dUdV[29]*mat[13]+dUdV[30]*mat[20]+dUdV[31]*mat[27]+dUdV[32]*mat[34];
  res[35] = dUdV[35]*mat[0]+dUdV[40]*mat[35];
  res[36] = dUdV[35]*mat[1]+dUdV[40]*mat[36];
  res[37] = dUdV[35]*mat[2]+dUdV[40]*mat[37];
  res[38] = dUdV[35]*mat[3]+dUdV[40]*mat[38];
  res[39] = dUdV[35]*mat[4]+dUdV[40]*mat[39];
  res[40] = dUdV[35]*mat[5]+dUdV[40]*mat[40];
  res[41] = dUdV[35]*mat[6]+dUdV[40]*mat[41];
  res[42] = dUdV[42]*mat[0]+dUdV[48]*mat[42];
  res[43] = dUdV[42]*mat[1]+dUdV[48]*mat[43];
  res[44] = dUdV[42]*mat[2]+dUdV[48]*mat[44];
  res[45] = dUdV[42]*mat[3]+dUdV[48]*mat[45];
  res[46] = dUdV[42]*mat[4]+dUdV[48]*mat[46];
  res[47] = dUdV[42]*mat[5]+dUdV[48]*mat[47];
  res[48] = dUdV[42]*mat[6]+dUdV[48]*mat[48];
}
//------------------------------------------------------------------------------
inline
void VarFcnTaitKE::postMultiplyBydVdU(double *V, double *mat, double *res)
{
  double dVdU[49];
  computedVdU(V, dVdU);

  res[0] = mat[0]*dVdU[0]+mat[1]*dVdU[7]+mat[2]*dVdU[14]+mat[3]*dVdU[21]+
    mat[4]*dVdU[28]+mat[5]*dVdU[35]+mat[6]*dVdU[42];
  res[1] = mat[1]*dVdU[8]+mat[4]*dVdU[29];
  res[2] = mat[2]*dVdU[16]+mat[4]*dVdU[30];
  res[3] = mat[3]*dVdU[24]+mat[4]*dVdU[31];
  res[4] = mat[4]*dVdU[32];
  res[5] = mat[5]*dVdU[40];
  res[6] = mat[6]*dVdU[48];
  res[7] = mat[7]*dVdU[0]+mat[8]*dVdU[7]+mat[9]*dVdU[14]+mat[10]*dVdU[21]+
    mat[11]*dVdU[28]+mat[12]*dVdU[35]+mat[13]*dVdU[42];
  res[8] = mat[8]*dVdU[8]+mat[11]*dVdU[29];
  res[9] = mat[9]*dVdU[16]+mat[11]*dVdU[30];
  res[10] = mat[10]*dVdU[24]+mat[11]*dVdU[31];
  res[11] = mat[11]*dVdU[32];
  res[12] = mat[12]*dVdU[40];
  res[13] = mat[13]*dVdU[48];
  res[14] = mat[14]*dVdU[0]+mat[15]*dVdU[7]+mat[16]*dVdU[14]+mat[17]*dVdU[21]+
    mat[18]*dVdU[28]+mat[19]*dVdU[35]+mat[20]*dVdU[42];
  res[15] = mat[15]*dVdU[8]+mat[18]*dVdU[29];
  res[16] = mat[16]*dVdU[16]+mat[18]*dVdU[30];
  res[17] = mat[17]*dVdU[24]+mat[18]*dVdU[31];
  res[18] = mat[18]*dVdU[32];
  res[19] = mat[19]*dVdU[40];
  res[20] = mat[20]*dVdU[48];
  res[21] = mat[21]*dVdU[0]+mat[22]*dVdU[7]+mat[23]*dVdU[14]+mat[24]*dVdU[21]+
    mat[25]*dVdU[28]+mat[26]*dVdU[35]+mat[27]*dVdU[42];
  res[22] = mat[22]*dVdU[8]+mat[25]*dVdU[29];
  res[23] = mat[23]*dVdU[16]+mat[25]*dVdU[30];
  res[24] = mat[24]*dVdU[24]+mat[25]*dVdU[31];
  res[25] = mat[25]*dVdU[32];
  res[26] = mat[26]*dVdU[40];
  res[27] = mat[27]*dVdU[48];
  res[28] = mat[28]*dVdU[0]+mat[29]*dVdU[7]+mat[30]*dVdU[14]+mat[31]*dVdU[21]+
    mat[32]*dVdU[28]+mat[33]*dVdU[35]+mat[34]*dVdU[42];
  res[29] = mat[29]*dVdU[8]+mat[32]*dVdU[29];
  res[30] = mat[30]*dVdU[16]+mat[32]*dVdU[30];
  res[31] = mat[31]*dVdU[24]+mat[32]*dVdU[31];
  res[32] = mat[32]*dVdU[32];
  res[33] = mat[33]*dVdU[40];
  res[34] = mat[34]*dVdU[48];
  res[35] = mat[35]*dVdU[0]+mat[36]*dVdU[7]+mat[37]*dVdU[14]+mat[38]*dVdU[21]+
    mat[39]*dVdU[28]+mat[40]*dVdU[35]+mat[41]*dVdU[42];
  res[36] = mat[36]*dVdU[8]+mat[39]*dVdU[29];
  res[37] = mat[37]*dVdU[16]+mat[39]*dVdU[30];
  res[38] = mat[38]*dVdU[24]+mat[39]*dVdU[31];
  res[39] = mat[39]*dVdU[32];
  res[40] = mat[40]*dVdU[40];
  res[41] = mat[41]*dVdU[48];
  res[42] = mat[42]*dVdU[0]+mat[43]*dVdU[7]+mat[44]*dVdU[14]+mat[45]*dVdU[21]+
    mat[46]*dVdU[28]+mat[47]*dVdU[35]+mat[48]*dVdU[42];
  res[43] = mat[43]*dVdU[8]+mat[46]*dVdU[29];
  res[44] = mat[44]*dVdU[16]+mat[46]*dVdU[30];
  res[45] = mat[45]*dVdU[24]+mat[46]*dVdU[31];
  res[46] = mat[46]*dVdU[32];
  res[47] = mat[47]*dVdU[40];
  res[48] = mat[48]*dVdU[48];
}
//------------------------------------------------------------------------------
inline
void VarFcnTaitKE::postMultiplyBydUdV(double *V, double *mat, double *res)
{
  double dUdV[49];
  computedUdV(V, dUdV);

  for (int i=0; i<7; i++){
    res[0+7*i] = mat[0+7*i]*dUdV[0]+mat[1+7*i]*dUdV[7]+mat[2+7*i]*dUdV[14]+
      mat[3+7*i]*dUdV[21]+mat[4+7*i]*dUdV[28]+mat[5+7*i]*dUdV[35]+mat[6+7*i]*dUdV[42];
    res[1+7*i] = mat[1+7*i]*dUdV[8]+mat[4+7*i]*dUdV[29];
    res[2+7*i] = mat[2+7*i]*dUdV[16]+mat[4+7*i]*dUdV[30];
    res[3+7*i] = mat[3+7*i]*dUdV[24]+mat[4+7*i]*dUdV[31];
    res[4+7*i] = mat[4+7*i]*dUdV[32];
    res[5+7*i] = mat[5+7*i]*dUdV[40];
    res[6+7*i] = mat[6+7*i]*dUdV[48];
  }
}
//------------------------------------------------------------------------------

#endif
