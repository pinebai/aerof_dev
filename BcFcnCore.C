#include <cstdlib>
#include <cstdio>
#include <BcFcn.h>
#include <BcDef.h>
#include <IoData.h>


//------------------------------------------------------------------------------

void BcFcn::applyToSolutionVector(int t, double *v, double *u)
{

  fprintf(stderr, "*** Error: applyToSolutionVector function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToTurbSolutionVector(int t, double *v, double *u)
{

  fprintf(stderr, "*** Error: applyToTurbSolutionVector function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToResidualTerm(int t, double *v, double *u, double *f)
{

  fprintf(stderr, "*** Error: applyToResidualTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------
//
void BcFcn::applyToTurbResidualTerm(int t, double *v, double *u, double *f)
{

  fprintf(stderr, "*** Error: applyToTurbResidualTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

// Included (MB)
void BcFcn::applyToDerivativeOfResidualTerm(int t, double *v, double *dv, double *u, double *du, double *df)
{

  fprintf(stderr, "*** Error: applyToDerivativeOfResidualTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToDiagonalTerm(int t, double *v, double *u, float *a)
{

  fprintf(stderr, "*** Error: applyToDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToDiagonalTerm(int t, double *v, double *u, double *a)
{

  fprintf(stderr, "*** Error: applyToDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToDiagonalTerm(int t, double *v, double *u, bcomp *a)
{

  fprintf(stderr, "*** Error: applyToDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToTurbDiagonalTerm(int t, double *v, double *u, float *a)
{

  fprintf(stderr, "*** Error: applyToTurbDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToTurbDiagonalTerm(int t, double *v, double *u, double *a)
{

  fprintf(stderr, "*** Error: applyToTurbDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToTurbDiagonalTerm(int t, double *v, double *u, bcomp *a)
{

  fprintf(stderr, "*** Error: applyToTurbDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToOffDiagonalTerm(int t, float *a)
{

  fprintf(stderr, "*** Error: applyToOffDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToOffDiagonalTerm(int t, double *a)
{

  fprintf(stderr, "*** Error: applyToOffDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToOffDiagonalTerm(int t, bcomp *a)
{

  fprintf(stderr, "*** Error: applyToOffDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToTurbOffDiagonalTerm(int t, float *a)
{

  fprintf(stderr, "*** Error: applyToTurbOffDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToTurbOffDiagonalTerm(int t, double *a)
{

  fprintf(stderr, "*** Error: applyToTurbOffDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::applyToTurbOffDiagonalTerm(int t, bcomp *a)
{

  fprintf(stderr, "*** Error: applyToTurbOffDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::zeroDiagonalTerm(int t, float *a)
{

  fprintf(stderr, "*** Error: zeroDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::zeroDiagonalTerm(int t, double *a)
{

  fprintf(stderr, "*** Error: zeroDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

void BcFcn::zeroDiagonalTerm(int t, bcomp *a)
{

  fprintf(stderr, "*** Error: zeroDiagonalTerm function not implemented\n");
  exit(1);

}

//------------------------------------------------------------------------------

inline
void BcFcnNS::template_applyToSolutionVectorTerm(int type, double *Vwall, double *U)
{

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED) {
    U[1] = U[0] * Vwall[1];
    U[2] = U[0] * Vwall[2];
    U[3] = U[0] * Vwall[3];

    if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED)
      U[4] = U[0] * ( Vwall[4] + 0.5*(Vwall[1]*Vwall[1]+Vwall[2]*Vwall[2]+Vwall[3]*Vwall[3]) );
  }

}

//------------------------------------------------------------------------------

inline
void BcFcnNS::template_applyToResidualTerm(int type, double *Vwall, double *U, double *F)
{

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED) {
    F[1] = - Vwall[1] * U[0] + U[1];
    F[2] = - Vwall[2] * U[0] + U[2];
    F[3] = - Vwall[3] * U[0] + U[3];

    if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED)
      F[4] = (U[4] - U[0] * Vwall[4]) * U[0] - 0.5 * (U[1]*U[1] + U[2]*U[2] + U[3]*U[3]);
  }
  
#if defined(STRONG_INLET_BC)
  if (type == BC_INLET_MOVING || type == BC_INLET_FIXED) {
    F[0] = - Vwall[0] + U[0];
    F[1] = - Vwall[0]*Vwall[1] + U[1];
    F[2] = - Vwall[0]*Vwall[2] + U[2];
    F[3] = - Vwall[0]*Vwall[3] + U[3];
  }
#endif

#if defined(STRONG_FARFIELD_BC)
  if (type == BC_INLET_MOVING || type == BC_INLET_FIXED ||
      type == BC_OUTLET_MOVING || type == BC_OUTLET_FIXED) {
    F[0] = - Vwall[0] + U[0];
    F[1] = - Vwall[1] + U[1];
    F[2] = - Vwall[2] + U[2];
    F[3] = - Vwall[3] + U[3];
    F[4] = - Vwall[4] + U[4];
  }
#endif

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void BcFcnNS::template_applyToDerivativeOfResidualTerm(int type, double *Vwall, double *dVwall, double *U, double *dU, double *dF)
{

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED) {
    dF[1] = - dVwall[1] * U[0] - Vwall[1] * dU[0] + dU[1];
    dF[2] = - dVwall[2] * U[0] - Vwall[2] * dU[0] + dU[2];
    dF[3] = - dVwall[3] * U[0] - Vwall[3] * dU[0] + dU[3];

    if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED)
      dF[4] = (dU[4] - dU[0] * Vwall[4] - U[0] * dVwall[4]) * U[0] + (U[4] - U[0] * Vwall[4]) * dU[0] - (U[1]*dU[1] + U[2]*dU[2] + U[3]*dU[3]);
  }

#if defined(STRONG_INLET_BC)
  if (type == BC_INLET_MOVING || type == BC_INLET_FIXED) {
    dF[0] = - dVwall[0] + dU[0];
    dF[1] = - dVwall[0]*Vwall[1] - Vwall[0]*dVwall[1] + dU[1];
    dF[2] = - dVwall[0]*Vwall[2] - Vwall[0]*dVwall[2] + dU[2];
    dF[3] = - dVwall[0]*Vwall[3] - Vwall[0]*dVwall[3] + dU[3];
  }
#endif

#if defined(STRONG_FARFIELD_BC)
  if (type == BC_INLET_MOVING || type == BC_INLET_FIXED ||
      type == BC_OUTLET_MOVING || type == BC_OUTLET_FIXED) {
    dF[0] = - dVwall[0] + dU[0];
    dF[1] = - dVwall[1] + dU[1];
    dF[2] = - dVwall[2] + dU[2];
    dF[3] = - dVwall[3] + dU[3];
    dF[4] = - dVwall[4] + dU[4];
  }
#endif

}

//------------------------------------------------------------------------------

void BcFcnNS::applyToSolutionVector(int type, double *Vwall, double *U)
{

  template_applyToSolutionVectorTerm(type, Vwall, U);

}

//------------------------------------------------------------------------------

void BcFcnNS::applyToResidualTerm(int type, double *Vwall, double *U, double *F)
{

  template_applyToResidualTerm(type, Vwall, U, F);  

}

//------------------------------------------------------------------------------

// Included (MB)
void BcFcnNS::applyToDerivativeOfResidualTerm(int type, double *Vwall, double *dVwall, double *U, double *dU, double *dF)
{

  template_applyToDerivativeOfResidualTerm(type, Vwall, dVwall, U, dU, dF);

}

//------------------------------------------------------------------------------

void BcFcnNS::applyToDiagonalTerm(int type, double *Vwall, double *U, float *A)
{

  template_applyToDiagonalTerm<float,5>(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnNS::applyToDiagonalTerm(int type, double *Vwall, double *U, double *A)
{

  template_applyToDiagonalTerm<double,5>(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnNS::applyToDiagonalTerm(int type, double *Vwall, double *U, bcomp *A)
{

  template_applyToDiagonalTerm<bcomp,5>(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnNS::zeroDiagonalTerm(int type, float *A)
{

  template_zeroDiagonalTerm<float,5>(type, A);

}

//------------------------------------------------------------------------------

void BcFcnNS::zeroDiagonalTerm(int type, double *A)
{

  template_zeroDiagonalTerm<double,5>(type, A);

}

//------------------------------------------------------------------------------

void BcFcnNS::zeroDiagonalTerm(int type, bcomp *A)
{

  template_zeroDiagonalTerm<bcomp,5>(type, A);

}

//------------------------------------------------------------------------------

void BcFcnNS::applyToOffDiagonalTerm(int type, float *A)
{

  template_applyToOffDiagonalTerm<float,5>(type, A);

}

//------------------------------------------------------------------------------

void BcFcnNS::applyToOffDiagonalTerm(int type, double *A)
{

  template_applyToOffDiagonalTerm<double,5>(type, A);

}

//------------------------------------------------------------------------------

void BcFcnNS::applyToOffDiagonalTerm(int type, bcomp *A)
{

  template_applyToOffDiagonalTerm<bcomp,5>(type, A);

}

//------------------------------------------------------------------------------

BcFcnSA::BcFcnSA(IoData& iod)
{

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = true;
  else
    wallFcn = false;

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToSolutionVector(int type, double *Vwall, double *U)
{
  if (!wallFcn)
    BcFcnNS::template_applyToSolutionVectorTerm(type, Vwall, U);

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED)
    U[5] = U[0] * Vwall[5];

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToTurbSolutionVector(int type, double *Vwall, double *U)
{
    U[5] = U[0] * Vwall[5];
}

//------------------------------------------------------------------------------

void BcFcnSA::applyToResidualTerm(int type, double *Vwall, double *U, double *F)
{

  if (!wallFcn)
    BcFcnNS::template_applyToResidualTerm(type, Vwall, U, F);

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED)
    F[5] = - U[0] * Vwall[5] + U[5];

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToTurbResidualTerm(int type, double *Vwall, double *U, double *F)
{
    F[5] = - U[0] * Vwall[5] + U[5];
}

//------------------------------------------------------------------------------

// Included (MB)
void BcFcnSA::applyToDerivativeOfResidualTerm(int type, double *Vwall, double *dVwall, double *U, double *dU, double *dF)
{

  if (!wallFcn)
    BcFcnNS::template_applyToDerivativeOfResidualTerm(type, Vwall, dVwall, U, dU, dF);

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED)
    dF[5] = - dU[0] * Vwall[5] - U[0] * dVwall[5] + dU[5];

}

//------------------------------------------------------------------------------

// Included (MB)
void BcFcnSA::applyToProductTerm(int type, float *Prod)
{

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED)
    Prod[5] = 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
void BcFcnSA::applyToProductTerm(int type, double *Prod)
{

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED)
    Prod[5] = 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
void BcFcnSA::applyToProductTerm(int type, bcomp *Prod)
{

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED)
    Prod[5] = 0.0;

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToDiagonalTerm(int type, double *Vwall, double *U, float *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

// Included (MB)
void BcFcnSA::applyToDiagonalTerm(int type, double *Vwall, double *dVwall, double *U, float *A)
{

  template_applyToDiagonalTerm(type, Vwall, dVwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToDiagonalTerm(int type, double *Vwall, double *U, double *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToTurbDiagonalTerm(int type, double *Vwall, double *U, float *A)
{

  template_applyToTurbDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------


void BcFcnSA::applyToTurbDiagonalTerm(int type, double *Vwall, double *U, double *A)
{

  template_applyToTurbDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

// Included (MB)
void BcFcnSA::applyToDiagonalTerm(int type, double *Vwall, double *dVwall, double *U, double *A)
{

  template_applyToDiagonalTerm(type, Vwall, dVwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToDiagonalTerm(int type, double *Vwall, double *U, bcomp *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------
//
void BcFcnSA::applyToTurbDiagonalTerm(int type, double *Vwall, double *U, bcomp *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

// Included (MB)
void BcFcnSA::applyToDiagonalTerm(int type, double *Vwall, double *dVwall, double *U, bcomp *A)
{

  template_applyToDiagonalTerm(type, Vwall, dVwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToOffDiagonalTerm(int type, float *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToOffDiagonalTerm(int type, double *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToOffDiagonalTerm(int type, bcomp *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToTurbOffDiagonalTerm(int type, float *A)
{

  template_applyToTurbOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToTurbOffDiagonalTerm(int type, double *A)
{

  template_applyToTurbOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnSA::applyToTurbOffDiagonalTerm(int type, bcomp *A)
{

  template_applyToTurbOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToSolutionVector(int type, double *Vwall, double *U)
{

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED)
    U[5] = Vwall[5];


}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToTurbSolutionVector(int type, double *Vwall, double *U)
{
    U[5] = Vwall[5];
}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToDiagonalTerm(int type, double *Vwall, double *U, float *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToDiagonalTerm(int type, double *Vwall, double *U, double *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToDiagonalTerm(int type, double *Vwall, double *U, bcomp *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToTurbDiagonalTerm(int type, double *Vwall, double *U, float *A)
{

  template_applyToTurbDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToTurbDiagonalTerm(int type, double *Vwall, double *U, double *A)
{

  template_applyToTurbDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToTurbDiagonalTerm(int type, double *Vwall, double *U, bcomp *A)
{

  template_applyToTurbDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToOffDiagonalTerm(int type, float *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToOffDiagonalTerm(int type, double *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToOffDiagonalTerm(int type, bcomp *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToTurbOffDiagonalTerm(int type, float *A)
{

  template_applyToTurbOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToTurbOffDiagonalTerm(int type, double *A)
{

  template_applyToTurbOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnSAturb::applyToTurbOffDiagonalTerm(int type, bcomp *A)
{

  template_applyToTurbOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

BcFcnKE::BcFcnKE(IoData& iod)
{

  if (iod.bc.wall.integration == BcsWallData::WALL_FUNCTION)
    wallFcn = true;
  else
    wallFcn = false;

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToSolutionVector(int type, double *Vwall, double *U)
{
  if (!wallFcn)
    BcFcnNS::template_applyToSolutionVectorTerm(type, Vwall, U);

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED)
    U[5] = U[0] * Vwall[5];

}
//------------------------------------------------------------------------------

void BcFcnKE::applyToTurbSolutionVector(int type, double *Vwall, double *U)
{
  U[5] = U[0] * Vwall[5];
}
//------------------------------------------------------------------------------


void BcFcnKE::applyToResidualTerm(int type, double *Vwall, double *U, double *F)
{

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED) {
    F[5] = - Vwall[5] * U[0] + U[5];
    F[6] = - Vwall[6] * U[0] + U[6];
  }

}

//------------------------------------------------------------------------------
//
void BcFcnKE::applyToTurbResidualTerm(int type, double *Vwall, double *U, double *F)
{
  F[5] = - Vwall[5] * U[0] + U[5];
  F[6] = - Vwall[6] * U[0] + U[6];
}

//------------------------------------------------------------------------------

// Included (MB)
void BcFcnKE::applyToDerivativeOfResidualTerm(int type, double *Vwall, double *dVwall, double *U, double *dU, double *dF)
{

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED) {
    dF[5] = - dVwall[5] * U[0] - Vwall[5] * dU[0] + dU[5];
    dF[6] = - dVwall[6] * U[0] - Vwall[6] * dU[0] + dU[6];
  }

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToDiagonalTerm(int type, double *Vwall, double *U, float *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToDiagonalTerm(int type, double *Vwall, double *U, double *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToDiagonalTerm(int type, double *Vwall, double *U, bcomp *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToTurbDiagonalTerm(int type, double *Vwall, double *U, float *A)
{

  template_applyToTurbDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToTurbDiagonalTerm(int type, double *Vwall, double *U, double *A)
{

  template_applyToTurbDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToTurbDiagonalTerm(int type, double *Vwall, double *U, bcomp *A)
{

  template_applyToTurbDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToOffDiagonalTerm(int type, float *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToOffDiagonalTerm(int type, double *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToOffDiagonalTerm(int type, bcomp *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToTurbOffDiagonalTerm(int type, float *A)
{

  template_applyToTurbOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToTurbOffDiagonalTerm(int type, double *A)
{

  template_applyToTurbOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnKE::applyToTurbOffDiagonalTerm(int type, bcomp *A)
{

  template_applyToTurbOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnKEturb::applyToDiagonalTerm(int type, double *Vwall, double *U, float *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnKEturb::applyToDiagonalTerm(int type, double *Vwall, double *U, double *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnKEturb::applyToDiagonalTerm(int type, double *Vwall, double *U, bcomp *A)
{

  template_applyToDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnKEturb::applyToTurbDiagonalTerm(int type, double *Vwall, double *U, float *A)
{

  template_applyToTurbDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnKEturb::applyToTurbDiagonalTerm(int type, double *Vwall, double *U, double *A)
{

  template_applyToTurbDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnKEturb::applyToTurbDiagonalTerm(int type, double *Vwall, double *U, bcomp *A)
{

  template_applyToTurbDiagonalTerm(type, Vwall, U, A);

}

//------------------------------------------------------------------------------

void BcFcnKEturb::applyToOffDiagonalTerm(int type, float *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnKEturb::applyToOffDiagonalTerm(int type, double *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnKEturb::applyToOffDiagonalTerm(int type, bcomp *A)
{

  template_applyToOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnKEturb::applyToTurbOffDiagonalTerm(int type, float *A)
{

  template_applyToTurbOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnKEturb::applyToTurbOffDiagonalTerm(int type, double *A)
{

  template_applyToTurbOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------

void BcFcnKEturb::applyToTurbOffDiagonalTerm(int type, bcomp *A)
{

  template_applyToTurbOffDiagonalTerm(type, A);

}

//------------------------------------------------------------------------------
