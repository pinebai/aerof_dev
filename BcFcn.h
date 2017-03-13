#ifndef _BC_FCN_H_
#define _BC_FCN_H_

#include<complex>
typedef std::complex<double> bcomp;
class IoData;

//------------------------------------------------------------------------------

class BcFcn {

public:

  BcFcn() {}
  virtual ~BcFcn() {}

  virtual void applyToSolutionVector(int, double *, double *);
  virtual void applyToTurbSolutionVector(int, double *, double *);
  virtual void applyToResidualTerm(int, double *, double *, double *);
  virtual void applyToTurbResidualTerm(int, double *, double *, double *);
  virtual void applyToDiagonalTerm(int, double *, double *, float *);
  virtual void applyToDiagonalTerm(int, double *, double *, double *);
  virtual void applyToDiagonalTerm(int, double *, double *, bcomp *);
  virtual void applyToOffDiagonalTerm(int, float *);
  virtual void applyToOffDiagonalTerm(int, double *);
  virtual void applyToOffDiagonalTerm(int, bcomp *);

  virtual void applyToTurbDiagonalTerm(int, double *, double *, float *);
  virtual void applyToTurbDiagonalTerm(int, double *, double *, double *);
  virtual void applyToTurbDiagonalTerm(int, double *, double *, bcomp *);
  virtual void applyToTurbOffDiagonalTerm(int, float *);
  virtual void applyToTurbOffDiagonalTerm(int, double *);
  virtual void applyToTurbOffDiagonalTerm(int, bcomp *);

  virtual void zeroDiagonalTerm(int, float *);
  virtual void zeroDiagonalTerm(int, double *);
  virtual void zeroDiagonalTerm(int, bcomp *);

// Included (MB)
  virtual void applyToDerivativeOfResidualTerm(int, double *, double *, double *, double *, double *);
  virtual void applyToDiagonalTerm(int, double *, double *, double *, float *) = 0;
  virtual void applyToDiagonalTerm(int, double *, double *, double *, double *) = 0;
  virtual void applyToDiagonalTerm(int, double *, double *, double *, bcomp *) = 0;
  virtual void applyToProductTerm(int , float *) = 0;
  virtual void applyToProductTerm(int , double *) = 0;
  virtual void applyToProductTerm(int , bcomp *) = 0;

};

//------------------------------------------------------------------------------

class BcFcnNS : public BcFcn {

public:

  BcFcnNS() {}
  ~BcFcnNS() {}

  static void template_applyToSolutionVectorTerm(int, double *, double *);

  static void template_applyToResidualTerm(int, double *, double *, double *);

  template<class Scalar, int neq>
  static void template_applyToDiagonalTerm(int, double *, double *, Scalar *);
  
  template<class Scalar, int neq>
  static void template_applyToOffDiagonalTerm(int, Scalar *);

  void applyToSolutionVector(int, double *, double *);
  void applyToResidualTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, float *);
  void applyToDiagonalTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, bcomp *);
  void applyToOffDiagonalTerm(int, float *);
  void applyToOffDiagonalTerm(int, double *);
  void applyToOffDiagonalTerm(int, bcomp *);

  void zeroDiagonalTerm(int, float *);
  void zeroDiagonalTerm(int, double *);
  void zeroDiagonalTerm(int, bcomp *);
  template<class Scalar, int neq>
  static void template_zeroDiagonalTerm(int, Scalar *);

// Included (MB)
  static void template_applyToDerivativeOfResidualTerm(int, double *, double *, double *, double *, double *);
  void applyToDerivativeOfResidualTerm(int, double *, double *, double *, double *, double *);
  void applyToDiagonalTerm(int c, double *vw, double *dvw, double *u, float *a) {}
  void applyToDiagonalTerm(int c, double *vw, double *dvw, double *u, double *a) {}
  void applyToDiagonalTerm(int c, double *vw, double *dvw, double *u, bcomp *a) {}
  void applyToProductTerm(int c, float *f) {}
  void applyToProductTerm(int c, double *f) {}
  void applyToProductTerm(int c, bcomp *f) {}

};

//------------------------------------------------------------------------------
/*
@BOOK{white-74,
  author = "White, F. M.",
  title = "Viscous fluid flow",
  publisher = "McGraw-Hill",
  year = 1974,
} 
*/

class BcFcnSA : public BcFcn {

  bool wallFcn;

  template<class Scalar>
  void template_applyToDiagonalTerm(int, double *, double *, Scalar *);
 
  template<class Scalar>
  void template_applyToOffDiagonalTerm(int, Scalar *);

  template<class Scalar>
  void template_applyToTurbDiagonalTerm(int, double *, double *, Scalar *);

  template<class Scalar>
  void template_applyToTurbOffDiagonalTerm(int, Scalar *);

// Included (MB)
  template<class Scalar>
  void template_applyToDiagonalTerm(int, double *, double *, double *, Scalar *);

  template<class Scalar>
  void template_applyToDiagonalTerm(int, double , double *, double *, double *, double *, Scalar *);

public:

  BcFcnSA(IoData&);
  ~BcFcnSA() {}

  void applyToSolutionVector(int, double *, double *);
  void applyToTurbSolutionVector(int, double *, double *);
  void applyToResidualTerm(int, double *, double *, double *);
  void applyToTurbResidualTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, float *);
  void applyToDiagonalTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, bcomp *);
  void applyToOffDiagonalTerm(int, float *);
  void applyToOffDiagonalTerm(int, double *);
  void applyToOffDiagonalTerm(int, bcomp *);

  void applyToTurbDiagonalTerm(int, double *, double *, float *);
  void applyToTurbDiagonalTerm(int, double *, double *, double *);
  void applyToTurbDiagonalTerm(int, double *, double *, bcomp *);
  void applyToTurbOffDiagonalTerm(int, float *);
  void applyToTurbOffDiagonalTerm(int, double *);
  void applyToTurbOffDiagonalTerm(int, bcomp *);

// Included (MB)
  void applyToDerivativeOfResidualTerm(int, double *, double *, double *, double *, double *);
  void applyToDiagonalTerm(int , double *, double *, double *, float *);
  void applyToDiagonalTerm(int , double *, double *, double *, double *);
  void applyToDiagonalTerm(int , double *, double *, double *, bcomp *);
  void applyToProductTerm(int , float *);
  void applyToProductTerm(int , double *);
  void applyToProductTerm(int , bcomp *);

};

//------------------------------------------------------------------------------

class BcFcnSAturb : public BcFcn {

  template<class Scalar>
  void template_applyToDiagonalTerm(int, double *, double *, Scalar *);
 
  template<class Scalar>
  void template_applyToOffDiagonalTerm(int, Scalar *);

  template<class Scalar>
  void template_applyToTurbDiagonalTerm(int, double *, double *, Scalar *);
 
  template<class Scalar>
  void template_applyToTurbOffDiagonalTerm(int, Scalar *);


public:

  BcFcnSAturb() {}
  ~BcFcnSAturb() {}

  void applyToSolutionVector(int, double *, double *);
  void applyToTurbSolutionVector(int, double *, double *);
  //void applyToResidualTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, float *);
  void applyToDiagonalTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, bcomp *);
  void applyToOffDiagonalTerm(int, float *);
  void applyToOffDiagonalTerm(int, double *);
  void applyToOffDiagonalTerm(int, bcomp *);

  void applyToTurbDiagonalTerm(int, double *, double *, float *);
  void applyToTurbDiagonalTerm(int, double *, double *, double *);
  void applyToTurbDiagonalTerm(int, double *, double *, bcomp *);
  void applyToTurbOffDiagonalTerm(int, float *);
  void applyToTurbOffDiagonalTerm(int, double *);
  void applyToTurbOffDiagonalTerm(int, bcomp *);

// Included (MB)
  void applyToDiagonalTerm(int c, double *vw, double *dvw, double *u, float *a) {}
  void applyToDiagonalTerm(int c, double *vw, double *dvw, double *u, double *a) {}
  void applyToDiagonalTerm(int c, double *vw, double *dvw, double *u, bcomp *a) {}
  void applyToProductTerm(int c, float *f) {}
  void applyToProductTerm(int c, double *f) {}
  void applyToProductTerm(int c, bcomp *f) {}

};

//------------------------------------------------------------------------------
/*
@ARTICLE{jaeger-dhatt-92,
  author = "Jaeger, M. and Dhatt, G.",
  title = "An extended k--$\epsilon$ finite element model",
  journal = ijnmf,
  year = 1992,
  volume = 14,
  pages = "1325--1345",
} 
*/

class BcFcnKE : public BcFcn {

  bool wallFcn;

  template<class Scalar>
  void template_applyToDiagonalTerm(int, double *, double *, Scalar *);
 
  template<class Scalar>
  void template_applyToOffDiagonalTerm(int, Scalar *);

  template<class Scalar>
  void template_applyToTurbDiagonalTerm(int, double *, double *, Scalar *);
 
  template<class Scalar>
  void template_applyToTurbOffDiagonalTerm(int, Scalar *);

public:

  BcFcnKE(IoData&);
  ~BcFcnKE() {}

  void applyToSolutionVector(int, double *, double *);
  void applyToTurbSolutionVector(int, double *, double *);
  void applyToResidualTerm(int, double *, double *, double *);
  void applyToTurbResidualTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, float *);
  void applyToDiagonalTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, bcomp *);
  void applyToOffDiagonalTerm(int, float *);
  void applyToOffDiagonalTerm(int, double *);
  void applyToOffDiagonalTerm(int, bcomp *);

  void applyToTurbDiagonalTerm(int, double *, double *, float *);
  void applyToTurbDiagonalTerm(int, double *, double *, double *);
  void applyToTurbDiagonalTerm(int, double *, double *, bcomp *);
  void applyToTurbOffDiagonalTerm(int, float *);
  void applyToTurbOffDiagonalTerm(int, double *);
  void applyToTurbOffDiagonalTerm(int, bcomp *);

// Included (MB)
  void applyToDerivativeOfResidualTerm(int, double *, double *, double *, double *, double *);
  void applyToDiagonalTerm(int c, double *vw, double *dvw, double *u, float *a) {}
  void applyToDiagonalTerm(int c, double *vw, double *dvw, double *u, double *a) {}
  void applyToDiagonalTerm(int c, double *vw, double *dvw, double *u, bcomp *a) {}
  void applyToProductTerm(int c, float *f) {}
  void applyToProductTerm(int c, double *f) {}
  void applyToProductTerm(int c, bcomp *f) {}

};

//------------------------------------------------------------------------------

class BcFcnKEturb : public BcFcn {

  template<class Scalar>
  void template_applyToDiagonalTerm(int, double *, double *, Scalar *);
 
  template<class Scalar>
  void template_applyToOffDiagonalTerm(int, Scalar *);

  template<class Scalar>
  void template_applyToTurbDiagonalTerm(int, double *, double *, Scalar *);
 
  template<class Scalar>
  void template_applyToTurbOffDiagonalTerm(int, Scalar *);

public:

  BcFcnKEturb() {}
  ~BcFcnKEturb() {}

  //void applyToResidualTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, float *);
  void applyToDiagonalTerm(int, double *, double *, double *);
  void applyToDiagonalTerm(int, double *, double *, bcomp *);
  void applyToOffDiagonalTerm(int, float *);
  void applyToOffDiagonalTerm(int, double *);
  void applyToOffDiagonalTerm(int, bcomp *);

  void applyToTurbDiagonalTerm(int, double *, double *, float *);
  void applyToTurbDiagonalTerm(int, double *, double *, double *);
  void applyToTurbDiagonalTerm(int, double *, double *, bcomp *);
  void applyToTurbOffDiagonalTerm(int, float *);
  void applyToTurbOffDiagonalTerm(int, double *);
  void applyToTurbOffDiagonalTerm(int, bcomp *);

// Included (MB)
  void applyToDiagonalTerm(int c, double *vw, double *dvw, double *u, float *a) {}
  void applyToDiagonalTerm(int c, double *vw, double *dvw, double *u, double *a) {}
  void applyToDiagonalTerm(int c, double *vw, double *dvw, double *u, bcomp *a) {}
  void applyToProductTerm(int c, float *f) {}
  void applyToProductTerm(int c, double *f) {}
  void applyToProductTerm(int c, bcomp *f) {}

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <BcFcn.C>
#endif

#endif
