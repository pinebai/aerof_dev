#include <BcFcn.h>
#include <BcDef.h>

//------------------------------------------------------------------------------

template<class Scalar, int neq>
inline
void BcFcnNS::template_applyToDiagonalTerm(int type, double *Vwall, double *U, Scalar *A)
{

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED) {
    for (int k=0; k<neq; ++k) {
      A[neq*1 + k] = 0.0;
      A[neq*2 + k] = 0.0;
      A[neq*3 + k] = 0.0;
      if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED)
	A[neq*4 + k] = 0.0;
    }
	 
    A[neq*1 + 0] = - Vwall[1];
    A[neq*1 + 1] = 1.0;

    A[neq*2 + 0] = - Vwall[2];
    A[neq*2 + 2] = 1.0;

    A[neq*3 + 0] = - Vwall[3];
    A[neq*3 + 3] = 1.0;

    if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED) {
      A[neq*4 + 0] = U[4] - 2.0*Vwall[4]*U[0];
      A[neq*4 + 1] = - U[1];
      A[neq*4 + 2] = - U[2];
      A[neq*4 + 3] = - U[3];
      A[neq*4 + 4] = U[0];
    }
  }
  
#if defined(STRONG_INLET_BC)
  if (type == BC_INLET_MOVING || type == BC_INLET_FIXED) {
    for (int k=0; k<neq; ++k) {
      A[neq*0 + k] = 0.0;
      A[neq*1 + k] = 0.0;
      A[neq*2 + k] = 0.0;
      A[neq*3 + k] = 0.0;
    }
	 
    A[neq*0 + 0] = 1.0;
    A[neq*1 + 1] = 1.0;
    A[neq*2 + 2] = 1.0;
    A[neq*3 + 3] = 1.0;
  }
#endif

#if defined(STRONG_FARFIELD_BC)
  if (type == BC_INLET_MOVING || type == BC_INLET_FIXED ||
      type == BC_OUTLET_MOVING || type == BC_OUTLET_FIXED) {
    for (int k=0; k<neq; ++k) {
      A[neq*0 + k] = 0.0;
      A[neq*1 + k] = 0.0;
      A[neq*2 + k] = 0.0;
      A[neq*3 + k] = 0.0;
      A[neq*4 + k] = 0.0;
    }
	 
    A[neq*0 + 0] = 1.0;
    A[neq*1 + 1] = 1.0;
    A[neq*2 + 2] = 1.0;
    A[neq*3 + 3] = 1.0;
    A[neq*4 + 4] = 1.0;
  }
#endif

}

//------------------------------------------------------------------------------

template<class Scalar, int neq>
inline
void BcFcnNS::template_applyToOffDiagonalTerm(int type, Scalar *A)
{

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED) {
    for (int k=0; k<neq; ++k) {
      A[neq*1 + k] = 0.0;
      A[neq*2 + k] = 0.0;
      A[neq*3 + k] = 0.0;
      if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED)
	A[neq*4 + k] = 0.0;
    }
  }

#if defined(STRONG_INLET_BC)
  if (type == BC_INLET_MOVING || type == BC_INLET_FIXED) {
    for (int k=0; k<neq; ++k) {
      A[neq*0 + k] = 0.0;
      A[neq*1 + k] = 0.0;
      A[neq*2 + k] = 0.0;
      A[neq*3 + k] = 0.0;
    }
  }
#endif

#if defined(STRONG_FARFIELD_BC)
  if (type == BC_INLET_MOVING || type == BC_INLET_FIXED ||
      type == BC_OUTLET_MOVING || type == BC_OUTLET_FIXED) {
    for (int k=0; k<neq; ++k) {
      A[neq*0 + k] = 0.0;
      A[neq*1 + k] = 0.0;
      A[neq*2 + k] = 0.0;
      A[neq*3 + k] = 0.0;
      A[neq*4 + k] = 0.0;
    }
  }
#endif

}

//------------------------------------------------------------------------------

template<class Scalar, int neq>
inline
void BcFcnNS::template_zeroDiagonalTerm(int type, Scalar *A)  {

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED) {
    for (int k=0; k<neq; ++k) {
      A[neq*1 + k] = 0.0;
      A[neq*2 + k] = 0.0;
      A[neq*3 + k] = 0.0;
      if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED)
        A[neq*4 + k] = 0.0;
    }

    A[neq*1 + 0] = 0.0;
    A[neq*1 + 1] = 0.0;

    A[neq*2 + 0] = 0.0;
    A[neq*2 + 2] = 0.0;

    A[neq*3 + 0] = 0.0;
    A[neq*3 + 3] = 0.0;

    if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED) {

      A[neq*4 + 0] = 0.0;
      A[neq*4 + 1] = 0.0;
      A[neq*4 + 2] = 0.0;
      A[neq*4 + 3] = 0.0;
      A[neq*4 + 4] = 0.0;

    }
  }

#if defined(STRONG_INLET_BC)
  if (type == BC_INLET_MOVING || type == BC_INLET_FIXED) {
    for (int k=0; k<neq; ++k) {
      A[neq*0 + k] = 0.0;
      A[neq*1 + k] = 0.0;
      A[neq*2 + k] = 0.0;
      A[neq*3 + k] = 0.0;

      A[neq*k+k] = 0.0;
    }
  }
#endif

#if defined(STRONG_FARFIELD_BC)
  if (type == BC_INLET_MOVING || type == BC_INLET_FIXED ||
      type == BC_OUTLET_MOVING || type == BC_OUTLET_FIXED) {
    for (int k=0; k<neq; ++k) {
      A[neq*0 + k] = 0.0;
      A[neq*1 + k] = 0.0;
      A[neq*2 + k] = 0.0;
      A[neq*3 + k] = 0.0;
      A[neq*4 + k] = 0.0;
    }
  }
#endif

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
void BcFcnSA::template_applyToDiagonalTerm(int type, double *Vwall, double *U, Scalar *A)
{

  const int neq = 6;

  if (!wallFcn)
    BcFcnNS::template_applyToDiagonalTerm<Scalar,neq>(type, Vwall, U, A);

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED) {
    for (int k=0; k<neq; ++k)
      A[neq*5 + k] = 0.0;
    A[neq*5 + 5] = 1.0;
  }

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
void BcFcnSA::template_applyToTurbDiagonalTerm(int type, double *Vwall, double *U, Scalar *A)
{

  const int neq = 6;

  for (int k=0; k<neq; ++k)
    A[neq*5 + k] = 0.0;
  A[neq*5 + 5] = 1.0;

}

//------------------------------------------------------------------------------

// Included (MB)
template<class Scalar>
inline
void BcFcnSA::template_applyToDiagonalTerm(int type, double *Vwall, double *dVwall, double *U, Scalar *A)
{

  const int neq = 6;

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED) {
    A[neq*5 + 0] = - Vwall[5] - U[0] * dVwall[0];
    A[neq*5 + 1] = - U[0] * dVwall[1];
    A[neq*5 + 2] = - U[0] * dVwall[2];
    A[neq*5 + 3] = - U[0] * dVwall[3];
    A[neq*5 + 4] = - U[0] * dVwall[4];
    A[neq*5 + 5] = - U[0] * dVwall[5] + 1.0;
  }

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
void BcFcnSA::template_applyToOffDiagonalTerm(int type, Scalar *A)
{

  const int neq = 6;

  if (!wallFcn)
    BcFcnNS::template_applyToOffDiagonalTerm<Scalar,neq>(type, A);

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED) {
    for (int k=0; k<neq; ++k)
      A[neq*5 + k] = 0.0;
  }

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
void BcFcnSA::template_applyToTurbOffDiagonalTerm(int type, Scalar *A)
{

  const int neq = 6;

  for (int k=0; k<neq; ++k)
    A[neq*5 + k] = 0.0;

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
void BcFcnSAturb::template_applyToDiagonalTerm(int type, double *Vwall, double *U, Scalar *A)
{

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED)
    A[0] = 1.0;

}

//------------------------------------------------------------------------------
//
template<class Scalar>
inline
void BcFcnSAturb::template_applyToTurbDiagonalTerm(int type, double *Vwall, double *U, Scalar *A)
{

  A[0] = 1.0;

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
void BcFcnSAturb::template_applyToOffDiagonalTerm(int type, Scalar *A)
{

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED)
    A[0] = 0.0;

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
void BcFcnSAturb::template_applyToTurbOffDiagonalTerm(int type, Scalar *A)
{

  A[0] = 0.0;

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
void BcFcnKE::template_applyToDiagonalTerm(int type, double *Vwall, double *U, Scalar *A)
{

  const int neq = 7;

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED) {
    for (int k=0; k<neq; ++k) {
      A[neq*5 + k] = 0.0;
      A[neq*6 + k] = 0.0;
    }
    A[neq*5 + 0] = - Vwall[5];
    A[neq*5 + 5] = 1.0;
    A[neq*5 + 0] = - Vwall[6];
    A[neq*6 + 6] = 1.0;
  }

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
void BcFcnKE::template_applyToTurbDiagonalTerm(int type, double *Vwall, double *U, Scalar *A)
{

  const int neq = 7;

  for (int k=0; k<neq; ++k) {
    A[neq*5 + k] = 0.0;
    A[neq*6 + k] = 0.0;
  }
  A[neq*5 + 0] = - Vwall[5];
  A[neq*5 + 5] = 1.0;
  A[neq*5 + 0] = - Vwall[6];
  A[neq*6 + 6] = 1.0;

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
void BcFcnKE::template_applyToOffDiagonalTerm(int type, Scalar *A)
{

  const int neq = 7;

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED) {
    for (int k=0; k<neq; ++k) {
      A[neq*5 + k] = 0.0;
      A[neq*6 + k] = 0.0;
    }
  }

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
void BcFcnKE::template_applyToTurbOffDiagonalTerm(int type, Scalar *A)
{

  const int neq = 7;

  for (int k=0; k<neq; ++k) {
    A[neq*5 + k] = 0.0;
    A[neq*6 + k] = 0.0;
  }

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
void BcFcnKEturb::template_applyToDiagonalTerm(int type, double *Vwall, double *U, Scalar *A)
{

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED) {
    A[0] = 1.0;
    A[1] = 0.0;
    A[2] = 0.0;
    A[3] = 1.0;
  }

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
void BcFcnKEturb::template_applyToTurbDiagonalTerm(int type, double *Vwall, double *U, Scalar *A)
{

  A[0] = 1.0;
  A[1] = 0.0;
  A[2] = 0.0;
  A[3] = 1.0;

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
void BcFcnKEturb::template_applyToOffDiagonalTerm(int type, Scalar *A)
{

  if (type == BC_ISOTHERMAL_WALL_MOVING || type == BC_ISOTHERMAL_WALL_FIXED ||
      type == BC_ADIABATIC_WALL_MOVING || type == BC_ADIABATIC_WALL_FIXED) {
    A[0] = 0.0;
    A[1] = 0.0;
    A[2] = 0.0;
    A[3] = 0.0;
  }

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
void BcFcnKEturb::template_applyToTurbOffDiagonalTerm(int type, Scalar *A)
{

  A[0] = 0.0;
  A[1] = 0.0;
  A[2] = 0.0;
  A[3] = 0.0;

}

//------------------------------------------------------------------------------
