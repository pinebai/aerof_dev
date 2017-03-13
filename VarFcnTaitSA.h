#ifndef _VAR_FCN_TAIT_SA_H
#define _VAR_FCN_TAIT_SA_H

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
class VarFcnTaitSA : public VarFcnTait {

protected:

  void computedVdU(double *V, double *dVdU);
  void computedUdV(double *V, double *dUdV);

public:
  VarFcnTaitSA(FluidModelData &data);
  ~VarFcnTaitSA() {}

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
VarFcnTaitSA::VarFcnTaitSA(FluidModelData &data) : VarFcnTait(data, true) {

  if(data.fluid != FluidModelData::LIQUID){
    fprintf(stderr, "*** Error: FluidModelData is not of type Tait in VarFcnTaitSA.\n");
    fflush(stderr);
    exit(1);
  }

  pname = new const char*[6];
  pname[0] = "density";
  pname[1] = "x-velocity";
  pname[2] = "y-velocity";
  pname[3] = "z-velocity";
  pname[4] = "temperature";
  pname[5] = "nut";
}
//------------------------------------------------------------------------------
inline
void VarFcnTaitSA::conservativeToPrimitive(double *U, double *V){

  VarFcnTait::conservativeToPrimitive(U,V);
  V[5] = U[5] / U[0];

}
//------------------------------------------------------------------------------
inline
void VarFcnTaitSA::primitiveToConservative(double *V, double *U) {

  VarFcnTait::primitiveToConservative(V,U);
  U[5] = V[0] * V[5];

}
//------------------------------------------------------------------------------
inline
void VarFcnTaitSA::computedVdU(double *V, double *dVdU) {

  double invrho = 1.0 / V[0];
  double invrhoCp = invrho*invC_;
  dVdU[0]  = 1.0;
  dVdU[6]  = -invrho * V[1];
  dVdU[7]  = invrho;
  dVdU[12] = -invrho * V[2];
  dVdU[14] = invrho;
  dVdU[18] = -invrho * V[3];
  dVdU[21] = invrho;
  dVdU[24] = invrhoCp * (getPressureDerivative(V) + 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3])) - V[4]*invrho;
  dVdU[25] = -invrhoCp * V[1];
  dVdU[26] = -invrhoCp * V[2];
  dVdU[27] = -invrhoCp * V[3];
  dVdU[28] = invrhoCp;
  dVdU[30] = -invrho * V[5];
  dVdU[35] = invrho;

}
//------------------------------------------------------------------------------
inline
void VarFcnTaitSA::computedUdV(double *V, double *dUdV) {

  dUdV[0]  = 1.0;
  dUdV[6]  = V[1];
  dUdV[7]  = V[0];
  dUdV[12] = V[2];
  dUdV[14] = V[0];
  dUdV[18] = V[3];
  dUdV[21] = V[0];
  dUdV[24] = C_ * V[4] + 0.5 * (V[1]*V[1] + V[2]*V[2] + V[3]*V[3]) - getPressureDerivative(V);
  dUdV[25] = V[0] * V[1];
  dUdV[26] = V[0] * V[2];
  dUdV[27] = V[0] * V[3];
  dUdV[28] = V[0] * C_;
  dUdV[30] = V[5];
  dUdV[35] = V[0];

}
//------------------------------------------------------------------------------
inline
void VarFcnTaitSA::multiplyBydVdU(double *V, double *vec, double *res)
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
void VarFcnTaitSA::preMultiplyBydUdV(double *V, double *mat, double *res)
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
void VarFcnTaitSA::postMultiplyBydVdU(double *V, double *mat, double *res)
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
void VarFcnTaitSA::postMultiplyBydUdV(double *V, double *mat, double *res)
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
