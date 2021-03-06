/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE ARDSNSym.h.
   Arpack++ class ARluNonSymStdEig definition
   (dense matrix version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARDSNSYM_H
#define ARDSNSYM_H

#include <cstddef>
#include "arch.h"
#include "arsnsym.h"
#include "ardnsmat.h"


template<class ARFLOAT>
class ARluNonSymStdEig:
  public virtual ARNonSymStdEig<ARFLOAT, ARdsNonSymMatrix<ARFLOAT, ARFLOAT> > {

  using ARNonSymStdEig<ARFLOAT, ARdsNonSymMatrix<ARFLOAT, ARFLOAT> >::sigmaR;
  using ARNonSymStdEig<ARFLOAT, ARdsNonSymMatrix<ARFLOAT, ARFLOAT> >::sigmaI;
  using ARNonSymStdEig<ARFLOAT, ARdsNonSymMatrix<ARFLOAT, ARFLOAT> >::mode;
  using ARNonSymStdEig<ARFLOAT, ARdsNonSymMatrix<ARFLOAT, ARFLOAT> >::iparam;
  using ARNonSymStdEig<ARFLOAT, ARdsNonSymMatrix<ARFLOAT, ARFLOAT> >::objOP;
  using ARNonSymStdEig<ARFLOAT, ARdsNonSymMatrix<ARFLOAT, ARFLOAT> >::Restart;
  using ARStdEig<ARFLOAT, ARFLOAT, ARdsNonSymMatrix<ARFLOAT, ARFLOAT> >::NoShift;
  using ARStdEig<ARFLOAT, ARFLOAT, ARdsNonSymMatrix<ARFLOAT, ARFLOAT> >::ClearMem;
 public:

 // a) Public functions:

 // a.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(ARFLOAT sigmaRp);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(ARFLOAT sigmap);

 // a.2) Constructors and destructor.

  ARluNonSymStdEig() { }
  // Short constructor.

  ARluNonSymStdEig(int nevp, ARdsNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                   char* whichp = "LM", int ncvp = 0,
                   ARFLOAT tolp = 0.0, int maxitp = 0,
                   ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARluNonSymStdEig(int nevp, ARdsNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                   ARFLOAT sigma, char* whichp = "LM", int ncvp = 0,
                   ARFLOAT tolp = 0.0, int maxitp = 0,
                   ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARluNonSymStdEig(const ARluNonSymStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARluNonSymStdEig() { }
  // Destructor.

 // b) Operators.

  ARluNonSymStdEig& operator=(const ARluNonSymStdEig& other);
  // Assignment operator.

}; // class ARluNonSymStdEig.


// ------------------------------------------------------------------------ //
// ARluNonSymStdEig member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARluNonSymStdEig<ARFLOAT>::
ChangeShift(ARFLOAT sigmaRp)
{

  sigmaR    = sigmaRp;
  sigmaI    = 0.0;
  mode      = 3;
  iparam[7] = mode;

  objOP->FactorAsI(sigmaR);
  Restart();

} // ChangeShift.


template<class ARFLOAT>
inline void ARluNonSymStdEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, ARFLOAT, ARdsNonSymMatrix<ARFLOAT, ARFLOAT> >::
    SetRegularMode(objOP, &ARdsNonSymMatrix<ARFLOAT, ARFLOAT>::MultMv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARluNonSymStdEig<ARFLOAT>::SetShiftInvertMode(ARFLOAT sigmap)
{

  ARStdEig<ARFLOAT, ARFLOAT, ARdsNonSymMatrix<ARFLOAT, ARFLOAT> >::
    SetShiftInvertMode(sigmap, objOP, 
                       &ARdsNonSymMatrix<ARFLOAT, ARFLOAT>::MultInvv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline ARluNonSymStdEig<ARFLOAT>::
ARluNonSymStdEig(int nevp, ARdsNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                 char* whichp, int ncvp, ARFLOAT tolp,
                 int maxitp, ARFLOAT* residp, bool ishiftp)

{

  NoShift();
  this->DefineParameters(A.ncols(), nevp, &A, 
                         &ARdsNonSymMatrix<ARFLOAT, ARFLOAT>::MultMv,
                         whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARluNonSymStdEig<ARFLOAT>::
ARluNonSymStdEig(int nevp, ARdsNonSymMatrix<ARFLOAT, ARFLOAT>& A,
                 ARFLOAT sigmap, char* whichp, int ncvp, ARFLOAT tolp,
                 int maxitp, ARFLOAT* residp, bool ishiftp)

{

  DefineParameters(A.ncols(), nevp, &A, 
                   &ARdsNonSymMatrix<ARFLOAT, ARFLOAT>::MultInvv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);
  ChangeShift(sigmap);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARluNonSymStdEig<ARFLOAT>& ARluNonSymStdEig<ARFLOAT>::
operator=(const ARluNonSymStdEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARDSNSYM_H

