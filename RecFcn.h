#ifndef _REC_FCN_H_
#define _REC_FCN_H_

#include <cstdio>
#include <cstdlib>
#include <cmath>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
using std::max;
#endif

//------------------------------------------------------------------------------

struct RecFcnBase {
  

  virtual void precompute(double*, double*, double*, double*, 
			  double*, double*, double*, double*) = 0;

  virtual void compute(double*, double*, double*, double*, double*, double*) = 0;
  virtual void compute(double*, double*, double*, double*, double*, double*, double, double) = 0;
  virtual void computeExtended(double*, double*, double*, double*, double*, double*,
                               double pl, double pr,int , int) = 0;
  virtual void computeDerivative(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*)=0;
  virtual void computeDerivativeOperators(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*)=0;
};

class RecFcn : virtual public RecFcnBase {

  double beta;
  double beta1;
  double beta2;
  double eps0;
  double eps3;

public:

  RecFcn(double b, double e) 
  { 
    beta = b;
    beta1 = 0.5 - beta;
    beta2 = 1.0 - 2.0*beta;
    eps0 = 1.e-16;
    eps3 = e*e*e;

  }
  virtual ~RecFcn() {}

  void preconstant(double&, double&, double&, double&);

  void constant(double, double, double&, double&);

  void prelinear(double&, double&, double&, double&);

  void linear(double, double, double, double, double&, double&);

  double computeVanAlbadaFcn(double, double, double);

  double computeDerivativeVanAlbadaFcn(double, double, double);

  void prevanalbada(double, double, double, double, 
		    double&, double&, double&, double&);

  void vanalbada(double, double, double, double, double&, double&);

  double computeBarthFcn(double, double, double, double, double);

  void barth(double, double, double, double, double,
	     double, double, double, double, double,
	     double&, double&);

  double computeVenkatFcn(double, double, double, double, double);

  void venkat(double, double, double, double, double,
	      double, double, double, double, double,
	      double&, double&);

  template<class Scalar, int dim>
  void compute(Scalar*, Scalar*, Scalar*, Scalar*, 
	       double*, double*, double*, double*, Scalar*, Scalar*);


  void precompute(double*, double*, double*, double*, 
		  double*, double*, double*, double*);

  void compute(double*, double*, double*, double*, double*, double*);
  void compute(double*, double*, double*, double*, double*, double*, double, double);

  void computeExtended(double*, double*, double*, double*, double*, double*,
                       double pl, double pr,int,int);
  
  template<int dim, class Scalar1, class Scalar2, class Scalar3>
  void compute(double *, SVec<Scalar1,dim> &, SVec<Scalar1,dim> &,
                SVec<Scalar1,dim> &, SVec<Scalar1,dim> &, Scalar2 *, Scalar2 *,
                Scalar2 *, Scalar2 *, int, int, Scalar3 *, Scalar3 *);

  template<class Scalar, int dim>
  void computeT(double *, Scalar *, Scalar *, double *, double *,
                double *, double *, int, int, Scalar *, Scalar *, Scalar *,
                Scalar *, Scalar *, Scalar *, Scalar *, Scalar *);

  template<int dim, class Scalar>
  void computeTb(SVec<Scalar,dim> &, SVec<Scalar,dim> &,
                SVec<Scalar,dim> &, SVec<Scalar,dim> &, int,
                Scalar *);

// Included (MB)
  void preconstantDerivative(double&, double&, double&, double&);

  void prelinearDerivative(double&, double&, double&, double&);

  void prevanalbadaDerivative(double, double, double, double, double, double, double, double,
		    double&, double&, double&, double&);

  double computeDerivativeOfDerivativeVanAlbadaFcn(double, double, double, double, double);

  void constantDerivative(double, double, double&, double&);

  void linearDerivative(double, double, double, double, double&, double&);

  void vanalbadaDerivative(double, double, double, double, double, double, double, double, double&, double&);
  void vanalbadaDerivativeOperators(double, double, double, double, double&, double&, double&, double&, double&, double&);
  double computeDerivativeOfBarthFcn(double, double, double, double, double, double, double, double, double, double);

  void barthDerivative(double, double, double, double, double, double, double, double, double, double,
	     double, double, double, double, double, double, double, double, double, double,
	     double&, double&, double&, double&);

  double computeDerivativeOfVanAlbadaFcn(double, double, double, double, double);
  void computeDerivativeOperatorsOfVanAlbadaFcn(double, double, double, double &, double &);

  double computeDerivativeOfVenkatFcn(double, double, double, double, double, double, double, double, double, double);

  void venkatDerivative(double, double, double, double, double, double, double, double, double, double,
	      double, double, double, double, double, double, double, double, double, double,
	      double&, double&, double&, double&);

  template<int dim>
  void computeDerivative(double*, double*, double*, double*,
	       double*, double*, double*, double*, double*, double*);

  template<int dim>
  void computeDerivativeOperators(double*, double*, double*, double*,
	       double*, double*, double*, double*, double*, double*) {}

  virtual void precomputeDerivative(double*, double*, double*, double*, double*, double*, double*, double*,
			  double*, double*, double*, double*);

  virtual void computeDerivative(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
  virtual void computeDerivativeOperators(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*) = 0;

  void interface(double, double, double, double, double&, double&, double, double);
};

class RecLimiter : virtual public RecFcnBase {
  public:
  
  virtual void computeLimiter(double *, double *, double *, double *, double,
			      double *, double *, double *, double *, double,
			      double *, double *) = 0;

// Included (MB)
  virtual void computeDerivativeOfLimiter(double *, double *, double *, double *, double *, double *, double *, double *, double, double,
		      double *, double *, double *, double *, double *, double *, double *, double *, double, double,
		      double *, double *, double *, double *) = 0;
};

//------------------------------------------------------------------------------

inline
void RecFcn::preconstant(double& aij, double& aji, double& bij, double& bji)
{

  aij = aji = bij = bji = 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void RecFcn::preconstantDerivative(double& daij, double& daji, double& dbij, double& dbji)
{

  daij = daji = dbij = dbji = 0.0;

}

//------------------------------------------------------------------------------

inline
void RecFcn::constant(double Vi, double Vj, double& Vij, double& Vji)
{
  Vij = Vi;
  Vji = Vj;

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void RecFcn::constantDerivative(double dVi, double dVj, double& dVij, double& dVji)
{

  dVij = dVi;
  dVji = dVj;

}

//------------------------------------------------------------------------------

inline
void RecFcn::prelinear(double& aij, double& aji, double& bij, double& bji)
{

  aij = aji = beta1;
  bij = bji = beta;

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void RecFcn::prelinearDerivative(double& daij, double& daji, double& dbij, double& dbji)
{

  daij = daji = 0.0;
  dbij = dbji = 0.0;

}

//------------------------------------------------------------------------------

inline
void RecFcn::linear(double Vi, double ddVij, double Vj, double ddVji, 
		    double& Vij, double& Vji)
{
 
  double dV = beta1 * (Vj - Vi);
  Vij = Vi + dV + beta * ddVij;
  Vji = Vj - dV - beta * ddVji;

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void RecFcn::linearDerivative(double dVi, double dddVij, double dVj, double dddVji,
		    double& dVij, double& dVji)
{

  double ddV = beta1 * (dVj - dVi);

  dVij = dVi + ddV + beta * dddVij;
  dVji = dVj - ddV - beta * dddVji;

}

//------------------------------------------------------------------------------

inline
double RecFcn::computeVanAlbadaFcn(double eps, double a, double b)
{
  double a2 = a * a;
  double b2 = b * b;

  if (a*b > 0.0)
    return (a*(b2+eps) + b*(a2+eps)) / (a2 + b2 + 2.0*eps);
  else if (a*b == 0.0)
    return 0.5 * eps * (a + b) / (a2 + b2 + 2.0*eps);
  else
    return 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
inline
double RecFcn::computeDerivativeOfVanAlbadaFcn(double eps, double a, double da, double b, double db)
{

  double a2 = a * a;
  double da2 = 2.0*a*da;
  double b2 = b * b;
  double db2 = 2.0*b*db;

  if (a*b > 0.0)
    return ( ( da*(b2+eps) + a*db2 +db*(a2+eps) +b*da2 ) * ( a2 + b2 + 2.0*eps ) -  ( a*(b2+eps) + b*(a2+eps) ) * ( da2 + db2 ) ) / ( ( a2 + b2 + 2.0*eps ) * ( a2 + b2 + 2.0*eps ) );
  else if (a*b == 0.0)
    return 0.5 * eps * ( ( da + db ) * ( a2 + b2 + 2.0*eps ) - ( a + b ) * (da2 + db2) ) / ( (a2 + b2 + 2.0*eps) * (a2 + b2 + 2.0*eps) );
  else
    return 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void RecFcn::computeDerivativeOperatorsOfVanAlbadaFcn(double eps, double a, double b, double &dda, double &ddb)
{

  double a2 = a * a;
  double b2 = b * b;

  if (a*b > 0.0) {
    double c1 = ( a2 + b2 + 2.0*eps );
    double c2 = ( a*(b2+eps) + b*(a2+eps) );
    double c3 = 1.0 / ( ( a2 + b2 + 2.0*eps ) * ( a2 + b2 + 2.0*eps ) );
    dda = c3*( c1*(b2+eps+b*2.0*a) - c2*2.0*a );
    ddb = c3*( c1*(a*2.0*b+a2+eps) - c2*2.0*b );
  } else if (a*b == 0.0) {
    double c4 = 0.5*eps/( (a2 + b2 + 2.0*eps) * (a2 + b2 + 2.0*eps) );
    double c5 = ( a2 + b2 + 2.0*eps );
    dda = (c4*c5-c4*(a+b)*2.0*a); 
    ddb = (c4*c5 - c4*(a+b)*2.0*b);
  } else {
    dda = 0.0;
    ddb = 0.0; 
  }

}

//------------------------------------------------------------------------------

inline
double RecFcn::computeDerivativeVanAlbadaFcn(double eps, double a, double b)
{

  double a2 = a * a;
  double b2 = b * b;

  if (a*b > 0.0)
    return (b2+eps) * (b2-a2+2.*eps+2.*a*b) / ((a2+b2+2.*eps)*(a2+b2+2.*eps));
  else if (a*b == 0.0)
    return 0.5 * (b2+eps) * (b2-a2+2.*eps) / ((a2+b2+2.*eps)*(a2+b2+2.*eps));
  else
    return 0.0;

}

//------------------------------------------------------------------------------

// Included (MB)
inline
double RecFcn::computeDerivativeOfDerivativeVanAlbadaFcn(double eps, double a, double da, double b, double db)
{

  double a2 = a * a;
  double da2 = 2.0 * a * da;
  double b2 = b * b;
  double db2 = 2.0 * b * db;

  if (a*b > 0.0)
    return ( ((db2) * (b2-a2+2.*eps+2.*a*b) + (b2+eps) * (db2-da2+2.*da*b+2.*a*db))*((a2+b2+2.*eps)*(a2+b2+2.*eps)) - ((b2+eps) * (b2-a2+2.*eps+2.*a*b))*( 2.0*(a2+b2+2.*eps)*(da2+db2)) ) / ( ((a2+b2+2.*eps)*(a2+b2+2.*eps))*((a2+b2+2.*eps)*(a2+b2+2.*eps)) ) ;
  else if (a*b == 0.0)
    return 0.5*( ((db2) * (b2-a2+2.*eps) + (b2+eps) * (db2-da2))*((a2+b2+2.*eps)*(a2+b2+2.*eps)) - ((b2+eps) * (b2-a2+2.*eps))*(2.0*(a2+b2+2.*eps)*(da2+db2)) ) / ( ((a2+b2+2.*eps)*(a2+b2+2.*eps))*((a2+b2+2.*eps)*(a2+b2+2.*eps)) );
  else
    return 0.0;

}

//------------------------------------------------------------------------------

inline
void RecFcn::prevanalbada(double Vi, double ddVij, double Vj, double ddVji,
			  double& aij, double& aji, double& bij, double& bji)
{

  double dVji = Vj - Vi;
  double dV   = beta2 * dVji;

  double dVi = dV + 2.0 * beta * ddVij;
  double dVj = dV + 2.0 * beta * ddVji;

  double lai = computeDerivativeVanAlbadaFcn(eps0, dVji, dVi);
  double lbi = computeDerivativeVanAlbadaFcn(eps0, dVi, dVji);

  double laj = computeDerivativeVanAlbadaFcn(eps0, dVji, dVj);
  double lbj = computeDerivativeVanAlbadaFcn(eps0, dVj, dVji);

  aij = 0.5 * (lai + beta2*lbi);
  bij = beta * lbi;
  aji = 0.5 * (laj + beta2*lbj);
  bji = beta * lbj;

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void RecFcn::prevanalbadaDerivative(double Vi, double dVi, double ddVij, double dddVij, double Vj, double dVj, double ddVji, double dddVji,
			  double& daij, double& daji, double& dbij, double& dbji)
{

  double d_Vji = Vj - Vi;
  double d_V   = beta2 * d_Vji;

  double dd_Vji = dVj - dVi;
  double dd_V   = beta2 * dd_Vji;

  double d_Vi = d_V + 2.0 * beta * ddVij;
  double d_Vj = d_V + 2.0 * beta * ddVji;

  double dd_Vi = dd_V + 2.0 * beta * dddVij;
  double dd_Vj = dd_V + 2.0 * beta * dddVji;

  double dlai = computeDerivativeOfDerivativeVanAlbadaFcn(eps0, d_Vji, dd_Vji, d_Vi, dd_Vi);
  double dlbi = computeDerivativeOfDerivativeVanAlbadaFcn(eps0, d_Vi, dd_Vi, d_Vji, dd_Vji);

  double dlaj = computeDerivativeOfDerivativeVanAlbadaFcn(eps0, d_Vji, dd_Vji, d_Vj, dd_Vj);
  double dlbj = computeDerivativeOfDerivativeVanAlbadaFcn(eps0, d_Vj, dd_Vj, d_Vji, dd_Vji);

  daij = 0.5 * (dlai + beta2*dlbi);
  dbij = beta * dlbi;
  daji = 0.5 * (dlaj + beta2*dlbj);
  dbji = beta * dlbj;

}

//------------------------------------------------------------------------------

inline
void RecFcn::vanalbada(double Vi, double ddVij, double Vj, double ddVji,
		       double& Vij, double& Vji)
{

  double dVji = Vj - Vi;
  double dV   = beta2 * dVji;

  double dVi = dV + 2.0 * beta * ddVij;
  double dVj = dV + 2.0 * beta * ddVji;

  Vij = Vi + 0.5 * computeVanAlbadaFcn(eps0, dVi, dVji);
  Vji = Vj - 0.5 * computeVanAlbadaFcn(eps0, dVj, dVji);

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void RecFcn::vanalbadaDerivative(double Vi, double dVi, double ddVij, double dddVij, double Vj, double dVj, double ddVji, double dddVji,
		       double& dVij, double& dVji)
{

  double d_Vji = Vj - Vi;
  double d_V   = beta2 * d_Vji;

  double dd_Vji = dVj - dVi;
  double dd_V   = beta2 * dd_Vji;

  double d_Vi = d_V + 2.0 * beta * ddVij;
  double d_Vj = d_V + 2.0 * beta * ddVji;

  double dd_Vi = dd_V + 2.0 * beta * dddVij;
  double dd_Vj = dd_V + 2.0 * beta * dddVji;

  dVij = dVi + 0.5 * computeDerivativeOfVanAlbadaFcn(eps0, d_Vi, dd_Vi, d_Vji, dd_Vji);
  dVji = dVj - 0.5 * computeDerivativeOfVanAlbadaFcn(eps0, d_Vj, dd_Vj, d_Vji, dd_Vji);

}

//------------------------------------------------------------------------------

// Included (YC)
inline
void RecFcn::vanalbadaDerivativeOperators(double Vi, double ddVij, double Vj, double ddVji, 
		       double &dVijdVi, double &dVijdVj, double &dVijdddVij, double &dVjidVi, double &dVjidVj, double &dVjidddVji)
{

  double d_Vji = Vj - Vi;
//  double dd_Vji = dVj - dVi;

  double d_Vi = beta2 * (Vj - Vi) + 2.0 * beta * ddVij;
  double d_Vj = beta2 * (Vj - Vi) + 2.0 * beta * ddVji;

//  double dd_Vi = beta2 * (dVj - dVi) + 2.0 * beta * dddVij;
//  double dd_Vj = beta2 * (dVj - dVi) + 2.0 * beta * dddVji;

  double dd_VjidVi = -1.0;
  double dd_VjidVj =  1.0;

  double dd_VidVi = -beta2;
  double dd_VidVj =  beta2;
  double dd_VidddVij = 2.0*beta;
  double dd_VjdVi = -beta2;
  double dd_VjdVj =  beta2;
  double dd_VjdddVji = 2.0*beta;

//  dVij = dVi;
//  dVij += 0.5 * computeDerivativeOfVanAlbadaFcn(eps0, d_Vi, dd_Vi, d_Vji, dd_Vji);
//  dVji = dVj;
//  dVji += -0.5 * computeDerivativeOfVanAlbadaFcn(eps0, d_Vj, dd_Vj, d_Vji, dd_Vji);

  dVijdVi = 1.0;
  dVjidVj = 1.0;

  double dVijdd_Vi, dVijdd_Vji, dVjidd_Vj, dVjidd_Vji;

  computeDerivativeOperatorsOfVanAlbadaFcn(eps0, d_Vi, d_Vji, dVijdd_Vi, dVijdd_Vji);
  dVijdVi += 0.5*(dVijdd_Vi*dd_VidVi + dVijdd_Vji*dd_VjidVi);
  dVijdVj =  0.5*(dVijdd_Vi*dd_VidVj + dVijdd_Vji*dd_VjidVj); 
  dVijdddVij = 0.5*dVijdd_Vi*dd_VidddVij;

  computeDerivativeOperatorsOfVanAlbadaFcn(eps0, d_Vj, d_Vji, dVjidd_Vj, dVjidd_Vji);
  dVjidVi = -0.5*(dVjidd_Vj*dd_VjdVi + dVjidd_Vji*dd_VjidVi);
  dVjidVj += -0.5*(dVjidd_Vj*dd_VjdVj + dVjidd_Vji*dd_VjidVj);
  dVjidddVji = -0.5*dVjidd_Vj*dd_VjdddVji;
  
//  dVijdddVij =  0.5 * computeDerivativeOperatorsOfVanAlbadaFcn(eps0, d_Vi, dd_Vi, d_Vji, dd_Vji, dVijdd_Vi, dVijdd_Vji);
//  dVjidddVji = -0.5 * computeDerivativeOperatorsOfVanAlbadaFcn(eps0, d_Vj, dd_Vj, d_Vji, dd_Vji, dVjidd_Vj, dVjidd_Vji);

}

//------------------------------------------------------------------------------

inline
double RecFcn::computeBarthFcn(double vmax, double vmin, double v,
			       double vrec, double eps2)
{

  double eps = 1.e-3;

  double num, den = vrec - v;

  if (den > eps)
    num = vmax - v;
  else if (den < -eps)
    num = vmin - v;
  else
    return 1.0;

  return min(1.0, num / den);

}

//------------------------------------------------------------------------------

// Included (MB)
inline
double RecFcn::computeDerivativeOfBarthFcn(double vmax, double dvmax, double vmin, double dvmin, double v, double dv,
				double vrec, double dvrec, double eps2, double deps2)
{

  double eps = 1.e-3;

  double num, den = vrec - v;

  double dnum, dden = dvrec - dv;

  if (den > eps) {
    num = vmax - v;
    dnum = dvmax - dv;
  }
  else if (den < -eps) {
    num = vmin - v;
    dnum = dvmin - dv;
  }
  else
    return 0.0;

  if (min(1.0, num / den) == 1.0  )
    return 0.0;
  else
    return ( dnum * den - num * dden ) / ( den * den );

}

//------------------------------------------------------------------------------

inline
void RecFcn::barth(double Vimax, double Vimin, double Vi, double Vij, double ctrlVoli,
		   double Vjmax, double Vjmin, double Vj, double Vji, double ctrlVolj,
		   double& phii, double& phij)
{

  phii = min(phii, computeBarthFcn(Vimax, Vimin, Vi, Vij, eps3 * ctrlVoli));

  phij = min(phij, computeBarthFcn(Vjmax, Vjmin, Vj, Vji, eps3 * ctrlVolj));

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void RecFcn::barthDerivative(double Vimax, double dVimax, double Vimin, double dVimin, double Vi, double dVi, double Vij, double dVij, double ctrlVoli, double dCtrlVoli,
		    double Vjmax, double dVjmax, double Vjmin, double dVjmin, double Vj, double dVj, double Vji, double dVji, double ctrlVolj, double dCtrlVolj,
		    double& phii, double& dphii, double& phij, double& dphij)
{

  if (min(phii, computeBarthFcn(Vimax, Vimin, Vi, Vij, eps3 * ctrlVoli)) != phii) {
    phii = computeBarthFcn(Vimax, Vimin, Vi, Vij, eps3 * ctrlVoli);
    dphii = computeDerivativeOfBarthFcn(Vimax, dVimax, Vimin, dVimin, Vi, dVi, Vij, dVij, eps3 * ctrlVoli, eps3 * dCtrlVoli);
  }

  if (min(phij, computeBarthFcn(Vjmax, Vjmin, Vj, Vji, eps3 * ctrlVolj)) != phij) {
    phij = computeBarthFcn(Vjmax, Vjmin, Vj, Vji, eps3 * ctrlVolj);
    dphij = computeDerivativeOfBarthFcn(Vjmax, dVjmax, Vjmin, dVjmin, Vj, dVj, Vji, dVji, eps3 * ctrlVolj, eps3 * dCtrlVolj);
  }

}

//------------------------------------------------------------------------------

inline
double RecFcn::computeVenkatFcn(double vmax, double vmin, double v,
				double vrec, double eps2)
{

  double eps = 1.e-3;

  double num, den = vrec - v;

  if (den > eps)
    num = vmax - v;
  else if (den < -eps)
    num = vmin - v;
  else
    return 1.0;

  return (num*num + 2.0*num*den + eps2) / (num*num + 2.0*den*den + num*den + eps2);

}

//------------------------------------------------------------------------------

// Included (MB)
inline
double RecFcn::computeDerivativeOfVenkatFcn(double vmax, double dvmax, double vmin, double dvmin, double v, double dv,
				double vrec, double dvrec, double eps2, double deps2)
{

  double eps = 1.e-3;

  double num, den = vrec - v;

  double dnum, dden = dvrec - dv;

  if (den > eps) {
    num = vmax - v;
    dnum = dvmax - dv;
  }
  else if (den < -eps) {
    num = vmin - v;
    dnum = dvmin - dv;
  }
  else
    return 0.0;

  return ( ( 2.0*num*dnum + 2.0*dnum*den + 2.0*num*dden + deps2 ) * ( num*num + 2.0*den*den + num*den + eps2 ) - (num*num + 2.0*num*den + eps2) * (2.0*num*dnum + 4.0*den*dden + dnum*den + num*dden + deps2) ) / ( (num*num + 2.0*den*den + num*den + eps2) * (num*num + 2.0*den*den + num*den + eps2) );

}

//------------------------------------------------------------------------------

inline
void RecFcn::venkat(double Vimax, double Vimin, double Vi, double Vij, double ctrlVoli,
		    double Vjmax, double Vjmin, double Vj, double Vji, double ctrlVolj,
		    double& phii, double& phij)
{

  phii = min(phii, computeVenkatFcn(Vimax, Vimin, Vi, Vij, eps3 * ctrlVoli));

  phij = min(phij, computeVenkatFcn(Vjmax, Vjmin, Vj, Vji, eps3 * ctrlVolj));

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void RecFcn::venkatDerivative(double Vimax, double dVimax, double Vimin, double dVimin, double Vi, double dVi, double Vij, double dVij, double ctrlVoli, double dCtrlVoli,
		    double Vjmax, double dVjmax, double Vjmin, double dVjmin, double Vj, double dVj, double Vji, double dVji, double ctrlVolj, double dCtrlVolj,
		    double& phii, double& dphii, double& phij, double& dphij)
{

  if (min(phii, computeVenkatFcn(Vimax, Vimin, Vi, Vij, eps3 * ctrlVoli)) != phii) {
    phii = computeVenkatFcn(Vimax, Vimin, Vi, Vij, eps3 * ctrlVoli);
    dphii = computeDerivativeOfVenkatFcn(Vimax, dVimax, Vimin, dVimin, Vi, dVi, Vij, dVij, eps3 * ctrlVoli, eps3 * dCtrlVoli);
  }

  if (min(phij, computeVenkatFcn(Vjmax, Vjmin, Vj, Vji, eps3 * ctrlVolj)) != phij) {
    phij = min(phij, computeVenkatFcn(Vjmax, Vjmin, Vj, Vji, eps3 * ctrlVolj));
    dphij = computeDerivativeOfVenkatFcn(Vjmax, dVjmax, Vjmin, dVjmin, Vj, dVj, Vji, dVji, eps3 * ctrlVolj, eps3 * dCtrlVolj);
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
void RecFcn::compute(Scalar* pi, Scalar* ddpij, Scalar* pj, Scalar* ddpji,
		     double* aij, double* aji, double* bij, double* bji,
		     Scalar* pij, Scalar* pji)
{

  for (int k=0; k<dim; ++k) {
    Scalar dpji = pj[k] - pi[k];
    pij[k] = pi[k] + aij[k] * dpji + bij[k] * ddpij[k];
    pji[k] = pj[k] - aji[k] * dpji - bji[k] * ddpji[k];
  }

}

//------------------------------------------------------------------------------

inline
void RecFcn::precompute(double* Vi, double* ddVij, double* Vj, double* ddVji,
			double* aij, double* aji, double* bij, double* bji)
{

  fprintf(stderr, "*** Error: RecFcn::precompute is not overloaded\n");
  exit(1);

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void RecFcn::precomputeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
			double* daij, double* daji, double* dbij, double* dbji)
{

  fprintf(stderr, "*** Error: RecFcn::precomputeDerivative is not overloaded\n");
  exit(1);

}

//------------------------------------------------------------------------------

inline
void RecFcn::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
		     double* Vij, double* Vji)
{

  fprintf(stderr, "*** Error: RecFcn::compute is not overloaded\n");
  exit(1);

}

//------------------------------------------------------------------------------

inline
void RecFcn::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
		     double* Vij, double* Vji, double phii, double phij)
{

  fprintf(stderr, "*** Error: RecFcn::compute is not overloaded\n");
  exit(1);

}

inline
void RecFcn::computeExtended(double* Vi, double* ddVij, double* Vj, double* ddVji,
                             double* Vij, double* Vji,
                             double pl, double pr,int i, int j) {

  this->compute(Vi,ddVij,Vj,ddVji, Vij,Vji);
}

//------------------------------------------------------------------------------

// Included (MB)
inline
void RecFcn::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
		     double* Vij, double* Vji)
{

  fprintf(stderr, "*** Error: RecFcn::computeDerivative is not overloaded\n");
  exit(1);

}

//------------------------------------------------------------------------------

template<int dim, class Scalar1, class Scalar2, class Scalar3>
inline
void RecFcn::compute(double *dx, SVec<Scalar1,dim> &p,
             SVec<Scalar1,dim> &dpdx, SVec<Scalar1,dim> &dpdy,
             SVec<Scalar1,dim> &dpdz, Scalar2 *aij, Scalar2 *aji, Scalar2 *bij,
             Scalar2 *bji, int i, int j, Scalar3 *pij, Scalar3 *pji) {

#pragma ivdep
  for (int k=0; k<dim; ++k) {
    pij[k] = p[i][k] + aij[k] * (p[j][k]-p[i][k]) +
      bij[k] * (dx[0]*dpdx[i][k] + dx[1]*dpdy[i][k] + dx[2]*dpdz[i][k]);
    pji[k] = p[j][k] - aji[k] * (p[j][k]-p[i][k]) -
      bji[k] * (dx[0]*dpdx[j][k] + dx[1]*dpdy[j][k] + dx[2]*dpdz[j][k]);

  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
void RecFcn::computeT(double *dx, Scalar *pi, Scalar *pj,
             double *aij, double *aji, double *bij, double *bji, int i,
             int j, Scalar *zu_is1, Scalar *zu_is2, Scalar *zgx_is1, Scalar *zgx_is2,
             Scalar *zgy_is1, Scalar *zgy_is2, Scalar *zgz_is1, Scalar *zgz_is2)
{

    // Same notations as in the Fortran code
#pragma ivdep
  for (int k=0; k<dim; ++k) {
    zu_is1[k] = pj[k] - aij[k]*pj[k] + aji[k]*pi[k];
    zu_is2[k] = pi[k] - aji[k]*pi[k] + aij[k]*pj[k];
    zgx_is1[k] =  bij[k] * pj[k] * dx[0];
    zgx_is2[k] = -bji[k] * pi[k] * dx[0];
    zgy_is1[k] =  bij[k] * pj[k] * dx[1];
    zgy_is2[k] = -bji[k] * pi[k] * dx[1];
    zgz_is1[k] =  bij[k] * pj[k] * dx[2];
    zgz_is2[k] = -bji[k] * pi[k] * dx[2];
  }

}

//------------------------------------------------------------------------------

template<int dim, class Scalar>
inline
void RecFcn::computeTb(SVec<Scalar,dim> &p,
             SVec<Scalar,dim> &dpdx, SVec<Scalar,dim> &dpdy,
             SVec<Scalar,dim> &dpdz, int i, Scalar *p1)
{

#pragma ivdep
  for (int k=0; k<dim; ++k){
    p1[k] = p[i][k] + (dpdx[i][k] + dpdy[i][k] + dpdz[i][k]);
 }
}

//------------------------------------------------------------------------------

inline
void RecFcn::interface(double Vi, double ddVij, double Vj, double ddVji,
                    double& Vij, double& Vji, double phii, double phij)
{

  if(phii*phij<0.0){

    double absphii = max(phii, -phii);
    double absphij = max(phij, -phij);
    double invTot = 1.0/(absphii+absphij);

    Vij = Vi + absphii*invTot * ddVij;
    Vji = Vj - absphij*invTot * ddVji;

  }else if(phii*phij==0){
    Vij = Vi;
    Vji = Vj;
  }else{
    fprintf(stderr, "***Error: RecFcn::interface should not be called! Exiting\n");
    exit(1);
  }

}

//------------------------------------------------------------------------------



#endif
