#ifndef _REC_FCN_DESC_H_
#define _REC_FCN_DESC_H_

#include <RecFcn.h>

//------------------------------------------------------------------------------

template<int dim>
class RecFcnConstant : public RecFcn {

public:

  RecFcnConstant() : RecFcn(0.0, 0.0) {}
  ~RecFcnConstant() {}

  void precompute(double *, double *, double *, double *, 
		  double *, double *, double *, double *);
  void compute(double *, double *, double *, double *, double *, double *);
  void compute(double *, double *, double *, double *, double *, double *, double, double);

// Included (MB)
  void precomputeDerivative(double *, double *, double *, double *, double *, double *, double *, double *,
                                               double *, double *, double *, double *);

  void computeDerivative(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
  void computeDerivativeOperators(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*) {}

};

//------------------------------------------------------------------------------

template<int dim>
class RecFcnLinear : public RecFcn {

public:

  RecFcnLinear(double b, double e) : RecFcn(b, e) {}
  ~RecFcnLinear() {}

  void precompute(double *, double *, double *, double *, 
		  double *, double *, double *, double *); 
  void compute(double *, double *, double *, double *, double *, double *);
  void compute(double *, double *, double *, double *, double *, double *, double, double);

// Included (MB)
  void precomputeDerivative(double *, double *, double *, double *, double *, double *, double *, double *,
                                               double *, double *, double *, double *);

  void computeDerivative(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
  void computeDerivativeOperators(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*) {}

};

//------------------------------------------------------------------------------

template<int dim>
class RecFcnVanAlbada : public RecFcn {

public:

  RecFcnVanAlbada(double b, double e) : RecFcn(b, e) {}
  ~RecFcnVanAlbada() {}

  void precompute(double *, double *, double *, double *, 
		  double *, double *, double *, double *); 
  void compute(double *, double *, double *, double *, double *, double *);
  void compute(double *, double *, double *, double *, double *, double *, double, double);

// Included (MB)
  void precomputeDerivative(double *, double *, double *, double *, double *, double *, double *, double *,
                                               double *, double *, double *, double *);

  void computeDerivative(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
  void computeDerivativeOperators(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);

};

template<int dim>
class RecFcnExtendedVanAlbada : public RecFcnVanAlbada<dim> {

public:

  RecFcnExtendedVanAlbada(double b, double e, double pc, double rhoc,
                          double xi_p, double vc,double xi_rho) : RecFcnVanAlbada<dim>(b, e),
    pc(pc), rhoc(rhoc), xi_p(xi_p),xi_rho(xi_rho), vc(vc) {}
  ~RecFcnExtendedVanAlbada() {}
  
  void computeExtended(double *, double *, double *, double *, double *, double *, double, double,int,int);

private:

  double pc,rhoc,vc;
  double xi_p, xi_rho;
};


//------------------------------------------------------------------------------

template<int dim>
class RecFcnLtdMultiDim : public RecFcnLinear<dim>, public RecLimiter {

public:

  RecFcnLtdMultiDim(double b, double e) : RecFcnLinear<dim>(b, e) {}
  ~RecFcnLtdMultiDim() {}


};

//------------------------------------------------------------------------------

template<int dim>
class RecFcnBarth : public RecFcnLtdMultiDim<dim> {

public:

  RecFcnBarth(double b, double e) : RecFcnLtdMultiDim<dim>(b, e) {}
  ~RecFcnBarth() {}

  void computeLimiter(double *, double *, double *, double *, double,
		      double *, double *, double *, double *, double,
		      double *, double *);

// Included (MB)
  void computeDerivativeOfLimiter(double *, double *, double *, double *, double *, double *, double *, double *, double, double,
		      double *, double *, double *, double *, double *, double *, double *, double *, double, double,
		      double *, double *, double *, double *);

};

//------------------------------------------------------------------------------

template<int dim>
class RecFcnVenkat : public RecFcnLtdMultiDim<dim> {

public:

  RecFcnVenkat(double b, double e) : RecFcnLtdMultiDim<dim>(b, e) {}
  ~RecFcnVenkat() {}

  void computeLimiter(double *, double *, double *, double *, double,
		      double *, double *, double *, double *, double,
		      double *, double *);

// Included (MB)
  void computeDerivativeOfLimiter(double *, double *, double *, double *, double *, double *, double *, double *, double, double,
		      double *, double *, double *, double *, double *, double *, double *, double *, double, double,
		      double *, double *, double *, double *);

};

//------------------------------------------------------------------------------

template<int dim>
class RecFcnLinearConstant : public RecFcn {

public:

  RecFcnLinearConstant(double b, double e) : RecFcn(b, e) {}
  ~RecFcnLinearConstant() {}

  void compute(double *, double *, double *, double *, double *, double *);
  void compute(double *, double *, double *, double *, double *, double *, double, double);

// Included (MB)
  void precompute(double *, double *, double *, double *, 
		  double *, double *, double *, double *); 
  void precomputeDerivative(double *, double *, double *, double *, double *, double *, double *, double *,
                                               double *, double *, double *, double *);
  void computeDerivative(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
  void computeDerivativeOperators(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*) {}

};

//------------------------------------------------------------------------------

template<int dim>
class RecFcnVanAlbadaConstant : public RecFcn {

public:

  RecFcnVanAlbadaConstant(double b, double e) : RecFcn(b, e) {}
  ~RecFcnVanAlbadaConstant() {}

  void compute(double *, double *, double *, double *, double *, double *);
  void compute(double *, double *, double *, double *, double *, double *, double, double);

// Included (MB)
  void precompute(double *, double *, double *, double *, 
		  double *, double *, double *, double *); 
  void precomputeDerivative(double *, double *, double *, double *, double *, double *, double *, double *,
                                               double *, double *, double *, double *);
  void computeDerivative(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
  void computeDerivativeOperators(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*) {}

};

//------------------------------------------------------------------------------

template<int dim>
class RecFcnLinearVanAlbada : public RecFcn {

public:

  RecFcnLinearVanAlbada(double b, double e) : RecFcn(b, e) {}
  ~RecFcnLinearVanAlbada() {}

  void compute(double *, double *, double *, double *, double *, double *);
  void compute(double *, double *, double *, double *, double *, double *, double, double);

// Included (MB)
  void precompute(double *, double *, double *, double *, 
		  double *, double *, double *, double *); 
  void precomputeDerivative(double *, double *, double *, double *, double *, double *, double *, double *,
                                               double *, double *, double *, double *);
  void computeDerivative(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
  void computeDerivativeOperators(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*) {}
};

//------------------------------------------------------------------------------

class RecFcnLtdSensor {

  double threshold;

public:

  RecFcnLtdSensor(double eps) { threshold = eps; }
  ~RecFcnLtdSensor() {}
  
  virtual double getThreshold() { return threshold; }

};

//------------------------------------------------------------------------------

template<int dim>
class RecFcnLtdLinear : public RecFcnLinear<dim>, public RecFcnLtdSensor {

public:

  RecFcnLtdLinear(double b, double e) : RecFcnLinear<dim>(b, e), RecFcnLtdSensor(e) {}
  ~RecFcnLtdLinear() {}

};

//------------------------------------------------------------------------------

template<int dim>
class RecFcnLtdLinearConstant : public RecFcnLinearConstant<dim>, public RecFcnLtdSensor {

public:

  RecFcnLtdLinearConstant(double b, double e) : RecFcnLinearConstant<dim>(b, e), RecFcnLtdSensor(e) {}
  ~RecFcnLtdLinearConstant() {}

};

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnConstant<dim>::precompute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				     double* aij, double* aji, double* bij, double* bji)
{

  for (int k=0; k<dim; ++k)
    preconstant(aij[k], aji[k], bij[k], bji[k]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
inline
void RecFcnConstant<dim>::precomputeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				     double* daij, double* daji, double* dbij, double* dbji)
{

  for (int k=0; k<dim; ++k)
    preconstantDerivative(daij[k], daji[k], dbij[k], dbji[k]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnConstant<5>::precompute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				   double* aij, double* aji, double* bij, double* bji)
{

  preconstant(aij[0], aji[0], bij[0], bji[0]);
  preconstant(aij[1], aji[1], bij[1], bji[1]);
  preconstant(aij[2], aji[2], bij[2], bji[2]);
  preconstant(aij[3], aji[3], bij[3], bji[3]);
  preconstant(aij[4], aji[4], bij[4], bji[4]);

}

//------------------------------------------------------------------------------
// Included (MB)
template<>
inline
void RecFcnConstant<5>::precomputeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				     double* daij, double* daji, double* dbij, double* dbji)
{

  preconstantDerivative(daij[0], daji[0], dbij[0], dbji[0]);
  preconstantDerivative(daij[1], daji[1], dbij[1], dbji[1]);
  preconstantDerivative(daij[2], daji[2], dbij[2], dbji[2]);
  preconstantDerivative(daij[3], daji[3], dbij[3], dbji[3]);
  preconstantDerivative(daij[4], daji[4], dbij[4], dbji[4]);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnConstant<dim>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				  double* Vij, double* Vji)
{
  for (int k=0; k<dim; ++k)
    constant(Vi[k], Vj[k], Vij[k], Vji[k]);

}

//------------------------------------------------------------------------------
template<int dim>
inline
void RecFcnConstant<dim>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
                                  double* Vij, double* Vji, double phii, double phij)
{

  for (int k=0; k<dim; ++k)
    interface(Vi[k], ddVij[k], Vj[k], ddVji[k], Vij[k], Vji[k], phii, phij);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
inline
void RecFcnConstant<dim>::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				double* dVij, double* dVji)
{

  for (int k=0; k<dim; ++k)
    constantDerivative(dVi[k], dVj[k], dVij[k], dVji[k]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnConstant<5>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				double* Vij, double* Vji)
{

  constant(Vi[0], Vj[0], Vij[0], Vji[0]);
  constant(Vi[1], Vj[1], Vij[1], Vji[1]);
  constant(Vi[2], Vj[2], Vij[2], Vji[2]);
  constant(Vi[3], Vj[3], Vij[3], Vji[3]);
  constant(Vi[4], Vj[4], Vij[4], Vji[4]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnConstant<5>::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				double* dVij, double* dVji)
{

  constantDerivative(dVi[0], dVj[0], dVij[0], dVji[0]);
  constantDerivative(dVi[1], dVj[1], dVij[1], dVji[1]);
  constantDerivative(dVi[2], dVj[2], dVij[2], dVji[2]);
  constantDerivative(dVi[3], dVj[3], dVij[3], dVji[3]);
  constantDerivative(dVi[4], dVj[4], dVij[4], dVji[4]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnConstant<6>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				double* Vij, double* Vji)
{

  constant(Vi[0], Vj[0], Vij[0], Vji[0]);
  constant(Vi[1], Vj[1], Vij[1], Vji[1]);
  constant(Vi[2], Vj[2], Vij[2], Vji[2]);
  constant(Vi[3], Vj[3], Vij[3], Vji[3]);
  constant(Vi[4], Vj[4], Vij[4], Vji[4]);
  constant(Vi[5], Vj[5], Vij[5], Vji[5]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnConstant<6>::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				double* dVij, double* dVji)
{

  constantDerivative(dVi[0], dVj[0], dVij[0], dVji[0]);
  constantDerivative(dVi[1], dVj[1], dVij[1], dVji[1]);
  constantDerivative(dVi[2], dVj[2], dVij[2], dVji[2]);
  constantDerivative(dVi[3], dVj[3], dVij[3], dVji[3]);
  constantDerivative(dVi[4], dVj[4], dVij[4], dVji[4]);
  constantDerivative(dVi[5], dVj[5], dVij[5], dVji[5]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnConstant<7>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				double* Vij, double* Vji)
{

  constant(Vi[0], Vj[0], Vij[0], Vji[0]);
  constant(Vi[1], Vj[1], Vij[1], Vji[1]);
  constant(Vi[2], Vj[2], Vij[2], Vji[2]);
  constant(Vi[3], Vj[3], Vij[3], Vji[3]);
  constant(Vi[4], Vj[4], Vij[4], Vji[4]);
  constant(Vi[5], Vj[5], Vij[5], Vji[5]);
  constant(Vi[6], Vj[6], Vij[6], Vji[6]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnConstant<7>::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				double* dVij, double* dVji)
{

  constantDerivative(dVi[0], dVj[0], dVij[0], dVji[0]);
  constantDerivative(dVi[1], dVj[1], dVij[1], dVji[1]);
  constantDerivative(dVi[2], dVj[2], dVij[2], dVji[2]);
  constantDerivative(dVi[3], dVj[3], dVij[3], dVji[3]);
  constantDerivative(dVi[4], dVj[4], dVij[4], dVji[4]);
  constantDerivative(dVi[5], dVj[5], dVij[5], dVji[5]);
  constantDerivative(dVi[6], dVj[6], dVij[6], dVji[6]);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnLinear<dim>::precompute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				   double* aij, double* aji, double* bij, double* bji)
{

  for (int k=0; k<dim; ++k)
    prelinear(aij[k], aji[k], bij[k], bji[k]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
inline
void RecFcnLinear<dim>::precomputeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				     double* daij, double* daji, double* dbij, double* dbji)
{

  for (int k=0; k<dim; ++k)
    prelinearDerivative(daij[k], daji[k], dbij[k], dbji[k]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnLinear<5>::precompute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				 double* aij, double* aji, double* bij, double* bji)
{

  prelinear(aij[0], aji[0], bij[0], bji[0]);
  prelinear(aij[1], aji[1], bij[1], bji[1]);
  prelinear(aij[2], aji[2], bij[2], bji[2]);
  prelinear(aij[3], aji[3], bij[3], bji[3]);
  prelinear(aij[4], aji[4], bij[4], bji[4]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnLinear<5>::precomputeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				     double* daij, double* daji, double* dbij, double* dbji)
{

  prelinearDerivative(daij[0], daji[0], dbij[0], dbji[0]);
  prelinearDerivative(daij[1], daji[1], dbij[1], dbji[1]);
  prelinearDerivative(daij[2], daji[2], dbij[2], dbji[2]);
  prelinearDerivative(daij[3], daji[3], dbij[3], dbji[3]);
  prelinearDerivative(daij[4], daji[4], dbij[4], dbji[4]);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnLinear<dim>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				double* Vij, double* Vji)
{
  for (int k=0; k<dim; ++k)
    linear(Vi[k], ddVij[k], Vj[k], ddVji[k], Vij[k], Vji[k]);

}

//------------------------------------------------------------------------------
template<int dim>
inline
void RecFcnLinear<dim>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
                                  double* Vij, double* Vji, double phii, double phij)
{

  for (int k=0; k<dim; ++k)
    interface(Vi[k], ddVij[k], Vj[k], ddVji[k], Vij[k], Vji[k], phii, phij);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
inline
void RecFcnLinear<dim>::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				double* dVij, double* dVji)
{

  for (int k=0; k<dim; ++k)
    linearDerivative(dVi[k], dddVij[k], dVj[k], dddVji[k], dVij[k], dVji[k]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnLinear<5>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
			      double* Vij, double* Vji)
{

  linear(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  linear(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  linear(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  linear(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  linear(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnLinear<5>::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
			      double* dVij, double* dVji)
{

  linearDerivative(dVi[0], dddVij[0], dVj[0], dddVji[0], dVij[0], dVji[0]);
  linearDerivative(dVi[1], dddVij[1], dVj[1], dddVji[1], dVij[1], dVji[1]);
  linearDerivative(dVi[2], dddVij[2], dVj[2], dddVji[2], dVij[2], dVji[2]);
  linearDerivative(dVi[3], dddVij[3], dVj[3], dddVji[3], dVij[3], dVji[3]);
  linearDerivative(dVi[4], dddVij[4], dVj[4], dddVji[4], dVij[4], dVji[4]);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnVanAlbada<dim>::precompute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				      double* aij, double* aji, double* bij, double* bji)
{

  for (int k=0; k<dim; ++k)
    prevanalbada(Vi[k], ddVij[k], Vj[k], ddVji[k], aij[k], aji[k], bij[k], bji[k]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
inline
void RecFcnVanAlbada<dim>::precomputeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				     double* daij, double* daji, double* dbij, double* dbji)
{

  for (int k=0; k<dim; ++k)
    prevanalbadaDerivative(Vi[k], dVi[k], ddVij[k], dddVij[k], Vj[k], dVj[k], ddVji[k], dddVji[k], daij[k], daji[k], dbij[k], dbji[k]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnVanAlbada<5>::precompute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				    double* aij, double* aji, double* bij, double* bji)
{

  prevanalbada(Vi[0], ddVij[0], Vj[0], ddVji[0], aij[0], aji[0], bij[0], bji[0]);
  prevanalbada(Vi[1], ddVij[1], Vj[1], ddVji[1], aij[1], aji[1], bij[1], bji[1]);
  prevanalbada(Vi[2], ddVij[2], Vj[2], ddVji[2], aij[2], aji[2], bij[2], bji[2]);
  prevanalbada(Vi[3], ddVij[3], Vj[3], ddVji[3], aij[3], aji[3], bij[3], bji[3]);
  prevanalbada(Vi[4], ddVij[4], Vj[4], ddVji[4], aij[4], aji[4], bij[4], bji[4]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnVanAlbada<5>::precomputeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				     double* daij, double* daji, double* dbij, double* dbji)
{

  prevanalbadaDerivative(Vi[0], dVi[0], ddVij[0], dddVij[0], Vj[0], dVj[0], ddVji[0], dddVji[0], daij[0], daji[0], dbij[0], dbji[0]);
  prevanalbadaDerivative(Vi[1], dVi[1], ddVij[1], dddVij[1], Vj[1], dVj[1], ddVji[1], dddVji[1], daij[1], daji[1], dbij[1], dbji[1]);
  prevanalbadaDerivative(Vi[2], dVi[2], ddVij[2], dddVij[2], Vj[2], dVj[2], ddVji[2], dddVji[2], daij[2], daji[2], dbij[2], dbji[2]);
  prevanalbadaDerivative(Vi[3], dVi[3], ddVij[3], dddVij[3], Vj[3], dVj[3], ddVji[3], dddVji[3], daij[3], daji[3], dbij[3], dbji[3]);
  prevanalbadaDerivative(Vi[4], dVi[4], ddVij[4], dddVij[4], Vj[4], dVj[4], ddVji[4], dddVji[4], daij[4], daji[4], dbij[4], dbji[4]);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnVanAlbada<dim>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				   double* Vij, double* Vji)
{
  for (int k=0; k<dim; ++k)
    vanalbada(Vi[k], ddVij[k], Vj[k], ddVji[k], Vij[k], Vji[k]);
}

//------------------------------------------------------------------------------
template<int dim>
inline
void RecFcnVanAlbada<dim>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
                                  double* Vij, double* Vji, double phii, double phij)
{
  for (int k=0; k<dim; ++k)
    interface(Vi[k], ddVij[k], Vj[k], ddVji[k], Vij[k], Vji[k], phii, phij);
}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
inline
void RecFcnVanAlbada<dim>::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				 double* dVij, double* dVji)
{

  for (int k=0; k<dim; ++k)
    vanalbadaDerivative(Vi[k], dVi[k], ddVij[k], dddVij[k], Vj[k], dVj[k], ddVji[k], dddVji[k], dVij[k], dVji[k]);

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
inline
void RecFcnVanAlbada<dim>::computeDerivativeOperators(double* Vi, double* ddVij, double* Vj, double* ddVji, 
                                                      double* dVijdVi, double* dVijdVj, double* dVijdddVij, double* dVjidVi, double* dVjidVj, double *dVjidddVji)
{

  for (int k=0; k<dim; ++k)
    vanalbadaDerivativeOperators(Vi[k], ddVij[k], Vj[k], ddVji[k], dVijdVi[k], dVijdVj[k], dVijdddVij[k], dVjidVi[k], dVjidVj[k], dVjidddVji[k]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnVanAlbada<5>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				 double* Vij, double* Vji)
{
  vanalbada(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  vanalbada(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  vanalbada(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  vanalbada(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  vanalbada(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);
}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnVanAlbada<5>::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				 double* dVij, double* dVji)
{

  vanalbadaDerivative(Vi[0], dVi[0], ddVij[0], dddVij[0], Vj[0], dVj[0], ddVji[0], dddVji[0], dVij[0], dVji[0]);
  vanalbadaDerivative(Vi[1], dVi[1], ddVij[1], dddVij[1], Vj[1], dVj[1], ddVji[1], dddVji[1], dVij[1], dVji[1]);
  vanalbadaDerivative(Vi[2], dVi[2], ddVij[2], dddVij[2], Vj[2], dVj[2], ddVji[2], dddVji[2], dVij[2], dVji[2]);
  vanalbadaDerivative(Vi[3], dVi[3], ddVij[3], dddVij[3], Vj[3], dVj[3], ddVji[3], dddVji[3], dVij[3], dVji[3]);
  vanalbadaDerivative(Vi[4], dVi[4], ddVij[4], dddVij[4], Vj[4], dVj[4], ddVji[4], dddVji[4], dVij[4], dVji[4]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnVanAlbada<6>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				 double* Vij, double* Vji)
{

  vanalbada(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  vanalbada(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  vanalbada(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  vanalbada(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  vanalbada(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);
  vanalbada(Vi[5], ddVij[5], Vj[5], ddVji[5], Vij[5], Vji[5]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnVanAlbada<6>::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				 double* dVij, double* dVji)
{

  vanalbadaDerivative(Vi[0], dVi[0], ddVij[0], dddVij[0], Vj[0], dVj[0], ddVji[0], dddVji[0], dVij[0], dVji[0]);
  vanalbadaDerivative(Vi[1], dVi[1], ddVij[1], dddVij[1], Vj[1], dVj[1], ddVji[1], dddVji[1], dVij[1], dVji[1]);
  vanalbadaDerivative(Vi[2], dVi[2], ddVij[2], dddVij[2], Vj[2], dVj[2], ddVji[2], dddVji[2], dVij[2], dVji[2]);
  vanalbadaDerivative(Vi[3], dVi[3], ddVij[3], dddVij[3], Vj[3], dVj[3], ddVji[3], dddVji[3], dVij[3], dVji[3]);
  vanalbadaDerivative(Vi[4], dVi[4], ddVij[4], dddVij[4], Vj[4], dVj[4], ddVji[4], dddVji[4], dVij[4], dVji[4]);
  vanalbadaDerivative(Vi[5], dVi[5], ddVij[5], dddVij[5], Vj[5], dVj[5], ddVji[5], dddVji[5], dVij[5], dVji[5]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnVanAlbada<7>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				 double* Vij, double* Vji)
{

  vanalbada(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  vanalbada(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  vanalbada(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  vanalbada(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  vanalbada(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);
  vanalbada(Vi[5], ddVij[5], Vj[5], ddVji[5], Vij[5], Vji[5]);
  vanalbada(Vi[6], ddVij[6], Vj[6], ddVji[6], Vij[6], Vji[6]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnVanAlbada<7>::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				 double* dVij, double* dVji)
{

  vanalbadaDerivative(Vi[0], dVi[0], ddVij[0], dddVij[0], Vj[0], dVj[0], ddVji[0], dddVji[0], dVij[0], dVji[0]);
  vanalbadaDerivative(Vi[1], dVi[1], ddVij[1], dddVij[1], Vj[1], dVj[1], ddVji[1], dddVji[1], dVij[1], dVji[1]);
  vanalbadaDerivative(Vi[2], dVi[2], ddVij[2], dddVij[2], Vj[2], dVj[2], ddVji[2], dddVji[2], dVij[2], dVji[2]);
  vanalbadaDerivative(Vi[3], dVi[3], ddVij[3], dddVij[3], Vj[3], dVj[3], ddVji[3], dddVji[3], dVij[3], dVji[3]);
  vanalbadaDerivative(Vi[4], dVi[4], ddVij[4], dddVij[4], Vj[4], dVj[4], ddVji[4], dddVji[4], dVij[4], dVji[4]);
  vanalbadaDerivative(Vi[5], dVi[5], ddVij[5], dddVij[5], Vj[5], dVj[5], ddVji[5], dddVji[5], dVij[5], dVji[5]);
  vanalbadaDerivative(Vi[6], dVi[6], ddVij[6], dddVij[6], Vj[6], dVj[6], ddVji[6], dddVji[6], dVij[6], dVji[6]);

}

template<int dim>
inline
void RecFcnExtendedVanAlbada<dim>::computeExtended(double* Vi, double* ddVij, double* Vj, double* ddVji,
		  	   	                   double* Vij, double* Vji,double pi,double pj,int i, int j) {

  double phi_pi = std::max<double>((pi-pc)/xi_p,0);
  double phi_pj = std::max<double>((pj-pc)/xi_p,0);
  double phi_rhoi = std::max<double>((Vi[0]-rhoc)/xi_rho,0);
  double phi_rhoj = std::max<double>((Vj[0]-rhoc)/xi_rho,0);

  double phi_vi = (vc*vc)/std::max<double>(0.01,Vi[1]*Vi[1]+Vi[2]*Vi[2]+Vi[3]*Vi[3]);
  double phi_vj = (vc*vc)/std::max<double>(0.01,Vj[1]*Vj[1]+Vj[2]*Vj[2]+Vj[3]*Vj[3]);

  double phi_i = std::min<double>(phi_pi,phi_rhoi);
  double phi_j = std::min<double>(phi_pj,phi_rhoj);
  
  phi_i = std::min<double>(phi_i,phi_vi);
  phi_j = std::min<double>(phi_j,phi_vj);

  for (int k = 0; k < dim; ++k) {

    this->vanalbada(Vi[k], ddVij[k], Vj[k], ddVji[k], Vij[k], Vji[k]);
    if (fabs(Vj[k]-Vi[k]) > 1.0e-10) {
 
      Vij[k] = fabs(2.0*(Vij[k]-Vi[k])/(Vj[k]-Vi[k]));
      Vij[k] = Vi[k] + 0.5*std::min<double>(phi_i, Vij[k])*(Vj[k]-Vi[k]); 
      Vji[k] = fabs(2.0*(Vji[k]-Vj[k])/(Vj[k]-Vi[k]));
      Vji[k] = Vj[k] - 0.5*std::min<double>(phi_j, Vji[k])*(Vj[k]-Vi[k]);
    }
  }
}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnBarth<dim>::computeLimiter(double *Vimax, double *Vimin, double *Vi, 
				      double *Vij, double ctrlVoli,
				      double *Vjmax, double *Vjmin, double *Vj, 
				      double *Vji, double ctrlVolj,
				      double *phii, double *phij)
{

  for (int k=0; k<dim; ++k)
    this->barth(Vimax[k], Vimin[k], Vi[k], Vij[k], ctrlVoli,
	  Vjmax[k], Vjmin[k], Vj[k], Vji[k], ctrlVolj, phii[k], phij[k]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
inline
void RecFcnBarth<dim>::computeDerivativeOfLimiter(double *Vimax, double *dVimax, double *Vimin, double *dVimin, double *Vi, double *dVi,
				      double *Vij, double *dVij, double ctrlVoli, double dCtrlVoli,
				      double *Vjmax, double *dVjmax, double *Vjmin, double *dVjmin, double *Vj, double *dVj,
				      double *Vji, double *dVji, double ctrlVolj, double dCtrlVolj,
				      double *phii, double *dphii, double *phij, double *dphij)
{

  for (int k=0; k<dim; ++k)
    this->barthDerivative(Vimax[k], dVimax[k], Vimin[k], dVimin[k], Vi[k], dVi[k], Vij[k], dVij[k], ctrlVoli, dCtrlVoli,
	   Vjmax[k], dVjmax[k], Vjmin[k], dVjmin[k], Vj[k], dVj[k], Vji[k], dVji[k], ctrlVolj, dCtrlVolj, phii[k], dphii[k], phij[k], dphij[k]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnBarth<5>::computeLimiter(double *Vimax, double *Vimin, double *Vi, 
				    double *Vij, double ctrlVoli,
				    double *Vjmax, double *Vjmin, double *Vj, 
				    double *Vji, double ctrlVolj,
				    double *phii, double *phij)
{

  this->barth(Vimax[0], Vimin[0], Vi[0], Vij[0], ctrlVoli,
	Vjmax[0], Vjmin[0], Vj[0], Vji[0], ctrlVolj, phii[0], phij[0]);
  this->barth(Vimax[1], Vimin[1], Vi[1], Vij[1], ctrlVoli,
	Vjmax[1], Vjmin[1], Vj[1], Vji[1], ctrlVolj, phii[1], phij[1]);
  this->barth(Vimax[2], Vimin[2], Vi[2], Vij[2], ctrlVoli,
	Vjmax[2], Vjmin[2], Vj[2], Vji[2], ctrlVolj, phii[2], phij[2]);
  this->barth(Vimax[3], Vimin[3], Vi[3], Vij[3], ctrlVoli,
	Vjmax[3], Vjmin[3], Vj[3], Vji[3], ctrlVolj, phii[3], phij[3]);
  this->barth(Vimax[4], Vimin[4], Vi[4], Vij[4], ctrlVoli,
	Vjmax[4], Vjmin[4], Vj[4], Vji[4], ctrlVolj, phii[4], phij[4]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnBarth<5>::computeDerivativeOfLimiter(double *Vimax, double *dVimax, double *Vimin, double *dVimin, double *Vi, double *dVi,
				    double *Vij, double *dVij, double ctrlVoli, double dCtrlVoli,
				    double *Vjmax, double *dVjmax, double *Vjmin, double *dVjmin, double *Vj, double *dVj,
				    double *Vji, double *dVji, double ctrlVolj, double dCtrlVolj,
				    double *phii, double *dphii, double *phij, double *dphij)
{

  this->barthDerivative(Vimax[0], dVimax[0], Vimin[0], dVimin[0], Vi[0], dVi[0], Vij[0], dVij[0], ctrlVoli, dCtrlVoli,
	 Vjmax[0], dVjmax[0], Vjmin[0], dVjmin[0], Vj[0], dVj[0], Vji[0], dVji[0], ctrlVolj, dCtrlVolj, phii[0], dphii[0], phij[0], dphij[0]);
  this->barthDerivative(Vimax[1], dVimax[1], Vimin[1], dVimin[1], Vi[1], dVi[1], Vij[1], dVij[1], ctrlVoli, dCtrlVoli,
	 Vjmax[1], dVjmax[1], Vjmin[1], dVjmin[1], Vj[1], dVj[1], Vji[1], dVji[1], ctrlVolj, dCtrlVolj, phii[1], dphii[1], phij[1], dphij[1]);
  barthDerivative(Vimax[2], dVimax[2], Vimin[2], dVimin[2], Vi[2], dVi[2], Vij[2], dVij[2], ctrlVoli, dCtrlVoli,
	 Vjmax[2], dVjmax[2], Vjmin[2], dVjmin[2], Vj[2], dVj[2], Vji[2], dVji[2], ctrlVolj, dCtrlVolj, phii[2], dphii[2], phij[2], dphij[2]);
  this->barthDerivative(Vimax[3], dVimax[3], Vimin[3], dVimin[3], Vi[3], dVi[3], Vij[3], dVij[3], ctrlVoli, dCtrlVoli,
	 Vjmax[3], dVjmax[3], Vjmin[3], dVjmin[3], Vj[3], dVj[3], Vji[3], dVji[3], ctrlVolj, dCtrlVolj, phii[3], dphii[3], phij[3], dphij[3]);
  this->barthDerivative(Vimax[4], dVimax[4], Vimin[4], dVimin[4], Vi[4], dVi[4], Vij[4], dVij[4], ctrlVoli, dCtrlVoli,
	 Vjmax[4], dVjmax[4], Vjmin[4], dVjmin[4], Vj[4], dVj[4], Vji[4], dVji[4], ctrlVolj, dCtrlVolj, phii[4], dphii[4], phij[4], dphij[4]);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnVenkat<dim>::computeLimiter(double *Vimax, double *Vimin, double *Vi, 
				       double *Vij, double ctrlVoli,
				       double *Vjmax, double *Vjmin, double *Vj, 
				       double *Vji, double ctrlVolj,
				       double *phii, double *phij)
{

  for (int k=0; k<dim; ++k)
    this->venkat(Vimax[k], Vimin[k], Vi[k], Vij[k], ctrlVoli,
	   Vjmax[k], Vjmin[k], Vj[k], Vji[k], ctrlVolj, phii[k], phij[k]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
inline
void RecFcnVenkat<dim>::computeDerivativeOfLimiter(double *Vimax, double *dVimax, double *Vimin, double *dVimin, double *Vi, double *dVi,
				       double *Vij, double *dVij, double ctrlVoli, double dCtrlVoli,
				       double *Vjmax, double *dVjmax, double *Vjmin, double *dVjmin, double *Vj, double *dVj,
				       double *Vji, double *dVji, double ctrlVolj, double dCtrlVolj,
				       double *phii, double *dphii, double *phij, double *dphij)
{

  for (int k=0; k<dim; ++k)
    this->venkatDerivative(Vimax[k], dVimax[k], Vimin[k], dVimin[k], Vi[k], dVi[k], Vij[k], dVij[k], ctrlVoli, dCtrlVoli,
	   Vjmax[k], dVjmax[k], Vjmin[k], dVjmin[k], Vj[k], dVj[k], Vji[k], dVji[k], ctrlVolj, dCtrlVolj, phii[k], dphii[k], phij[k], dphij[k]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnVenkat<5>::computeLimiter(double *Vimax, double *Vimin, double *Vi, 
				     double *Vij, double ctrlVoli,
				     double *Vjmax, double *Vjmin, double *Vj, 
				     double *Vji, double ctrlVolj,
				     double *phii, double *phij)
{

  this->venkat(Vimax[0], Vimin[0], Vi[0], Vij[0], ctrlVoli,
	 Vjmax[0], Vjmin[0], Vj[0], Vji[0], ctrlVolj, phii[0], phij[0]);
  this->venkat(Vimax[1], Vimin[1], Vi[1], Vij[1], ctrlVoli,
	 Vjmax[1], Vjmin[1], Vj[1], Vji[1], ctrlVolj, phii[1], phij[1]);
  this->venkat(Vimax[2], Vimin[2], Vi[2], Vij[2], ctrlVoli,
	 Vjmax[2], Vjmin[2], Vj[2], Vji[2], ctrlVolj, phii[2], phij[2]);
  this->venkat(Vimax[3], Vimin[3], Vi[3], Vij[3], ctrlVoli,
	 Vjmax[3], Vjmin[3], Vj[3], Vji[3], ctrlVolj, phii[3], phij[3]);
  this->venkat(Vimax[4], Vimin[4], Vi[4], Vij[4], ctrlVoli,
	 Vjmax[4], Vjmin[4], Vj[4], Vji[4], ctrlVolj, phii[4], phij[4]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnVenkat<5>::computeDerivativeOfLimiter(double *Vimax, double *dVimax, double *Vimin, double *dVimin, double *Vi, double *dVi,
				     double *Vij, double *dVij, double ctrlVoli, double dCtrlVoli,
				     double *Vjmax, double *dVjmax, double *Vjmin, double *dVjmin, double *Vj, double *dVj,
				     double *Vji, double *dVji, double ctrlVolj, double dCtrlVolj,
				     double *phii, double *dphii, double *phij, double *dphij)
{

  this->venkatDerivative(Vimax[0], dVimax[0], Vimin[0], dVimin[0], Vi[0], dVi[0], Vij[0], dVij[0], ctrlVoli, dCtrlVoli,
	 Vjmax[0], dVjmax[0], Vjmin[0], dVjmin[0], Vj[0], dVj[0], Vji[0], dVji[0], ctrlVolj, dCtrlVolj, phii[0], dphii[0], phij[0], dphij[0]);
  this->venkatDerivative(Vimax[1], dVimax[1], Vimin[1], dVimin[1], Vi[1], dVi[1], Vij[1], dVij[1], ctrlVoli, dCtrlVoli,
	 Vjmax[1], dVjmax[1], Vjmin[1], dVjmin[1], Vj[1], dVj[1], Vji[1], dVji[1], ctrlVolj, dCtrlVolj, phii[1], dphii[1], phij[1], dphij[1]);
  this->venkatDerivative(Vimax[2], dVimax[2], Vimin[2], dVimin[2], Vi[2], dVi[2], Vij[2], dVij[2], ctrlVoli, dCtrlVoli,
	 Vjmax[2], dVjmax[2], Vjmin[2], dVjmin[2], Vj[2], dVj[2], Vji[2], dVji[2], ctrlVolj, dCtrlVolj, phii[2], dphii[2], phij[2], dphij[2]);
  this->venkatDerivative(Vimax[3], dVimax[3], Vimin[3], dVimin[3], Vi[3], dVi[3], Vij[3], dVij[3], ctrlVoli, dCtrlVoli,
	 Vjmax[3], dVjmax[3], Vjmin[3], dVjmin[3], Vj[3], dVj[3], Vji[3], dVji[3], ctrlVolj, dCtrlVolj, phii[3], dphii[3], phij[3], dphij[3]);
  this->venkatDerivative(Vimax[4], dVimax[4], Vimin[4], dVimin[4], Vi[4], dVi[4], Vij[4], dVij[4], ctrlVoli, dCtrlVoli,
	 Vjmax[4], dVjmax[4], Vjmin[4], dVjmin[4], Vj[4], dVj[4], Vji[4], dVji[4], ctrlVolj, dCtrlVolj, phii[4], dphii[4], phij[4], dphij[4]);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnLinearConstant<dim>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
					double* Vij, double* Vji)
{

  fprintf(stderr, "*** Error: RecFcnLinearConstant<%d>::compute is not overloaded\n", dim);
  exit(1);

}

//------------------------------------------------------------------------------
template<int dim>
inline
void RecFcnLinearConstant<dim>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
                                  double* Vij, double* Vji, double phii, double phij)
{

  fprintf(stderr, "*** Error: RecFcnLinearConstant<%d>::compute is not overloaded\n", dim);
  exit(1);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
inline
void RecFcnLinearConstant<dim>::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				      double* dVij, double* dVji)
{

  fprintf(stderr, "*** Error: RecFcnLinearConstant<%d>::computeDerivative is not overloaded\n", dim);
  exit(1);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnLinearConstant<6>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				      double* Vij, double* Vji)
{

  linear(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  linear(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  linear(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  linear(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  linear(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);
  constant(Vi[5], Vj[5], Vij[5], Vji[5]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnLinearConstant<6>::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				      double* dVij, double* dVji)
{

  linearDerivative(dVi[0], dddVij[0], dVj[0], dddVji[0], dVij[0], dVji[0]);
  linearDerivative(dVi[1], dddVij[1], dVj[1], dddVji[1], dVij[1], dVji[1]);
  linearDerivative(dVi[2], dddVij[2], dVj[2], dddVji[2], dVij[2], dVji[2]);
  linearDerivative(dVi[3], dddVij[3], dVj[3], dddVji[3], dVij[3], dVji[3]);
  linearDerivative(dVi[4], dddVij[4], dVj[4], dddVji[4], dVij[4], dVji[4]);
  constantDerivative(dVi[5], dVj[5], dVij[5], dVji[5]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnLinearConstant<6>::precompute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				 double* aij, double* aji, double* bij, double* bji)
{

  prelinear(aij[0], aji[0], bij[0], bji[0]);
  prelinear(aij[1], aji[1], bij[1], bji[1]);
  prelinear(aij[2], aji[2], bij[2], bji[2]);
  prelinear(aij[3], aji[3], bij[3], bji[3]);
  prelinear(aij[4], aji[4], bij[4], bji[4]);
  preconstant(aij[5], aji[5], bij[5], bji[5]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnLinearConstant<6>::precomputeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				     double* daij, double* daji, double* dbij, double* dbji)
{

  prelinearDerivative(daij[0], daji[0], dbij[0], dbji[0]);
  prelinearDerivative(daij[1], daji[1], dbij[1], dbji[1]);
  prelinearDerivative(daij[2], daji[2], dbij[2], dbji[2]);
  prelinearDerivative(daij[3], daji[3], dbij[3], dbji[3]);
  prelinearDerivative(daij[4], daji[4], dbij[4], dbji[4]);
  preconstantDerivative(daij[5], daji[5], dbij[5], dbji[5]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnLinearConstant<7>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				      double* Vij, double* Vji)
{

  linear(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  linear(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  linear(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  linear(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  linear(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);
  constant(Vi[5], Vj[5], Vij[5], Vji[5]);
  constant(Vi[6], Vj[6], Vij[6], Vji[6]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnLinearConstant<7>::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				      double* dVij, double* dVji)
{

  linearDerivative(dVi[0], dddVij[0], dVj[0], dddVji[0], dVij[0], dVji[0]);
  linearDerivative(dVi[1], dddVij[1], dVj[1], dddVji[1], dVij[1], dVji[1]);
  linearDerivative(dVi[2], dddVij[2], dVj[2], dddVji[2], dVij[2], dVji[2]);
  linearDerivative(dVi[3], dddVij[3], dVj[3], dddVji[3], dVij[3], dVji[3]);
  linearDerivative(dVi[4], dddVij[4], dVj[4], dddVji[4], dVij[4], dVji[4]);
  constantDerivative(dVi[5], dVj[5], dVij[5], dVji[5]);
  constantDerivative(dVi[6], dVj[6], dVij[6], dVji[6]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnLinearConstant<7>::precompute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				 double* aij, double* aji, double* bij, double* bji)
{

  prelinear(aij[0], aji[0], bij[0], bji[0]);
  prelinear(aij[1], aji[1], bij[1], bji[1]);
  prelinear(aij[2], aji[2], bij[2], bji[2]);
  prelinear(aij[3], aji[3], bij[3], bji[3]);
  prelinear(aij[4], aji[4], bij[4], bji[4]);
  preconstant(aij[5], aji[5], bij[5], bji[5]);
  preconstant(aij[6], aji[6], bij[6], bji[6]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnLinearConstant<7>::precomputeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				     double* daij, double* daji, double* dbij, double* dbji)
{

  prelinearDerivative(daij[0], daji[0], dbij[0], dbji[0]);
  prelinearDerivative(daij[1], daji[1], dbij[1], dbji[1]);
  prelinearDerivative(daij[2], daji[2], dbij[2], dbji[2]);
  prelinearDerivative(daij[3], daji[3], dbij[3], dbji[3]);
  prelinearDerivative(daij[4], daji[4], dbij[4], dbji[4]);
  preconstantDerivative(daij[5], daji[5], dbij[5], dbji[5]);
  preconstantDerivative(daij[6], daji[6], dbij[6], dbji[6]);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnVanAlbadaConstant<dim>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
					   double* Vij, double* Vji)
{

  fprintf(stderr, "*** Error: RecFcnVanAlbadaConstant<%d>::compute is not overloaded\n", dim);
  exit(1);

}

//------------------------------------------------------------------------------
template<int dim>
inline
void RecFcnVanAlbadaConstant<dim>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
                                  double* Vij, double* Vji, double phii, double phij)
{

  fprintf(stderr, "*** Error: RecFcnVanAlbadaConstant<%d>::compute is not overloaded\n", dim);
  exit(1);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
inline
void RecFcnVanAlbadaConstant<dim>::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
					 double* dVij, double* dVji)
{

  fprintf(stderr, "*** Error: RecFcnVanAlbadaConstant<%d>::computeDerivative is not overloaded\n", dim);
  exit(1);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnVanAlbadaConstant<6>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
					 double* Vij, double* Vji)
{

  vanalbada(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  vanalbada(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  vanalbada(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  vanalbada(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  vanalbada(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);
  constant(Vi[5], Vj[5], Vij[5], Vji[5]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnVanAlbadaConstant<6>::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
					 double* dVij, double* dVji)
{

  vanalbadaDerivative(Vi[0], dVi[0], ddVij[0], dddVij[0], Vj[0], dVj[0], ddVji[0], dddVji[0], dVij[0], dVji[0]);
  vanalbadaDerivative(Vi[1], dVi[1], ddVij[1], dddVij[1], Vj[1], dVj[1], ddVji[1], dddVji[1], dVij[1], dVji[1]);
  vanalbadaDerivative(Vi[2], dVi[2], ddVij[2], dddVij[2], Vj[2], dVj[2], ddVji[2], dddVji[2], dVij[2], dVji[2]);
  vanalbadaDerivative(Vi[3], dVi[3], ddVij[3], dddVij[3], Vj[3], dVj[3], ddVji[3], dddVji[3], dVij[3], dVji[3]);
  vanalbadaDerivative(Vi[4], dVi[4], ddVij[4], dddVij[4], Vj[4], dVj[4], ddVji[4], dddVji[4], dVij[4], dVji[4]);
  constantDerivative(dVi[5], dVj[5], dVij[5], dVji[5]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnVanAlbadaConstant<6>::precompute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				 double* aij, double* aji, double* bij, double* bji)
{

  prevanalbada(Vi[0], ddVij[0], Vj[0], ddVji[0], aij[0], aji[0], bij[0], bji[0]);
  prevanalbada(Vi[1], ddVij[1], Vj[1], ddVji[1], aij[1], aji[1], bij[1], bji[1]);
  prevanalbada(Vi[2], ddVij[2], Vj[2], ddVji[2], aij[2], aji[2], bij[2], bji[2]);
  prevanalbada(Vi[3], ddVij[3], Vj[3], ddVji[3], aij[3], aji[3], bij[3], bji[3]);
  prevanalbada(Vi[4], ddVij[4], Vj[4], ddVji[4], aij[4], aji[4], bij[4], bji[4]);
  preconstant(aij[5], aji[5], bij[5], bji[5]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnVanAlbadaConstant<6>::precomputeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				     double* daij, double* daji, double* dbij, double* dbji)
{

  prevanalbadaDerivative(Vi[0], dVi[0], ddVij[0], dddVij[0], Vj[0], dVj[0], ddVji[0], dddVji[0], daij[0], daji[0], dbij[0], dbji[0]);
  prevanalbadaDerivative(Vi[1], dVi[1], ddVij[1], dddVij[1], Vj[1], dVj[1], ddVji[1], dddVji[1], daij[1], daji[1], dbij[1], dbji[1]);
  prevanalbadaDerivative(Vi[2], dVi[2], ddVij[2], dddVij[2], Vj[2], dVj[2], ddVji[2], dddVji[2], daij[2], daji[2], dbij[2], dbji[2]);
  prevanalbadaDerivative(Vi[3], dVi[3], ddVij[3], dddVij[3], Vj[3], dVj[3], ddVji[3], dddVji[3], daij[3], daji[3], dbij[3], dbji[3]);
  prevanalbadaDerivative(Vi[4], dVi[4], ddVij[4], dddVij[4], Vj[4], dVj[4], ddVji[4], dddVji[4], daij[4], daji[4], dbij[4], dbji[4]);
  preconstantDerivative(daij[5], daji[5], dbij[5], dbji[5]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnVanAlbadaConstant<7>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
					 double* Vij, double* Vji)
{

  vanalbada(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  vanalbada(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  vanalbada(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  vanalbada(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  vanalbada(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);
  constant(Vi[5], Vj[5], Vij[5], Vji[5]);
  constant(Vi[6], Vj[6], Vij[6], Vji[6]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnVanAlbadaConstant<7>::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
					 double* dVij, double* dVji)
{

  vanalbadaDerivative(Vi[0], dVi[0], ddVij[0], dddVij[0], Vj[0], dVj[0], ddVji[0], dddVji[0], dVij[0], dVji[0]);
  vanalbadaDerivative(Vi[1], dVi[1], ddVij[1], dddVij[1], Vj[1], dVj[1], ddVji[1], dddVji[1], dVij[1], dVji[1]);
  vanalbadaDerivative(Vi[2], dVi[2], ddVij[2], dddVij[2], Vj[2], dVj[2], ddVji[2], dddVji[2], dVij[2], dVji[2]);
  vanalbadaDerivative(Vi[3], dVi[3], ddVij[3], dddVij[3], Vj[3], dVj[3], ddVji[3], dddVji[3], dVij[3], dVji[3]);
  vanalbadaDerivative(Vi[4], dVi[4], ddVij[4], dddVij[4], Vj[4], dVj[4], ddVji[4], dddVji[4], dVij[4], dVji[4]);
  constantDerivative(dVi[5], dVj[5], dVij[5], dVji[5]);
  constantDerivative(dVi[6], dVj[6], dVij[6], dVji[6]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnVanAlbadaConstant<7>::precompute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				 double* aij, double* aji, double* bij, double* bji)
{

  prevanalbada(Vi[0], ddVij[0], Vj[0], ddVji[0], aij[0], aji[0], bij[0], bji[0]);
  prevanalbada(Vi[1], ddVij[1], Vj[1], ddVji[1], aij[1], aji[1], bij[1], bji[1]);
  prevanalbada(Vi[2], ddVij[2], Vj[2], ddVji[2], aij[2], aji[2], bij[2], bji[2]);
  prevanalbada(Vi[3], ddVij[3], Vj[3], ddVji[3], aij[3], aji[3], bij[3], bji[3]);
  prevanalbada(Vi[4], ddVij[4], Vj[4], ddVji[4], aij[4], aji[4], bij[4], bji[4]);
  preconstant(aij[5], aji[5], bij[5], bji[5]);
  preconstant(aij[6], aji[6], bij[6], bji[6]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnVanAlbadaConstant<7>::precomputeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				     double* daij, double* daji, double* dbij, double* dbji)
{

  prevanalbadaDerivative(Vi[0], dVi[0], ddVij[0], dddVij[0], Vj[0], dVj[0], ddVji[0], dddVji[0], daij[0], daji[0], dbij[0], dbji[0]);
  prevanalbadaDerivative(Vi[1], dVi[1], ddVij[1], dddVij[1], Vj[1], dVj[1], ddVji[1], dddVji[1], daij[1], daji[1], dbij[1], dbji[1]);
  prevanalbadaDerivative(Vi[2], dVi[2], ddVij[2], dddVij[2], Vj[2], dVj[2], ddVji[2], dddVji[2], daij[2], daji[2], dbij[2], dbji[2]);
  prevanalbadaDerivative(Vi[3], dVi[3], ddVij[3], dddVij[3], Vj[3], dVj[3], ddVji[3], dddVji[3], daij[3], daji[3], dbij[3], dbji[3]);
  prevanalbadaDerivative(Vi[4], dVi[4], ddVij[4], dddVij[4], Vj[4], dVj[4], ddVji[4], dddVji[4], daij[4], daji[4], dbij[4], dbji[4]);
  preconstantDerivative(daij[5], daji[5], dbij[5], dbji[5]);
  preconstantDerivative(daij[6], daji[6], dbij[6], dbji[6]);

}

//------------------------------------------------------------------------------

template<int dim>
inline
void RecFcnLinearVanAlbada<dim>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
					 double* Vij, double* Vji)
{

  fprintf(stderr, "*** Error: RecFcnLinearVanAlbada<%d>::compute is not overloaded\n", dim);
  exit(1);

}

//------------------------------------------------------------------------------
template<int dim>
inline
void RecFcnLinearVanAlbada<dim>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
                                  double* Vij, double* Vji, double phii, double phij)
{

  fprintf(stderr, "*** Error: RecFcnLinearVanAlbada<%d>::compute is not overloaded\n", dim);
  exit(1);

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
inline
void RecFcnLinearVanAlbada<dim>::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
					 double* dVij, double* dVji)
{

  fprintf(stderr, "*** Error: RecFcnLinearVanAlbada<%d>::computeDerivative is not overloaded\n", dim);
  exit(1);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnLinearVanAlbada<6>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				       double* Vij, double* Vji)
{

  linear(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  linear(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  linear(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  linear(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  linear(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);
  vanalbada(Vi[5], ddVij[5], Vj[5], ddVji[5], Vij[5], Vji[5]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnLinearVanAlbada<6>::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				       double* dVij, double* dVji)
{

  linearDerivative(dVi[0], dddVij[0], dVj[0], dddVji[0], dVij[0], dVji[0]);
  linearDerivative(dVi[1], dddVij[1], dVj[1], dddVji[1], dVij[1], dVji[1]);
  linearDerivative(dVi[2], dddVij[2], dVj[2], dddVji[2], dVij[2], dVji[2]);
  linearDerivative(dVi[3], dddVij[3], dVj[3], dddVji[3], dVij[3], dVji[3]);
  linearDerivative(dVi[4], dddVij[4], dVj[4], dddVji[4], dVij[4], dVji[4]);
  vanalbadaDerivative(Vi[5], dVi[5], ddVij[5], dddVij[5], Vj[5], dVj[5], ddVji[5], dddVji[5], dVij[5], dVji[5]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnLinearVanAlbada<6>::precompute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				 double* aij, double* aji, double* bij, double* bji)
{

  prelinear(aij[0], aji[0], bij[0], bji[0]);
  prelinear(aij[1], aji[1], bij[1], bji[1]);
  prelinear(aij[2], aji[2], bij[2], bji[2]);
  prelinear(aij[3], aji[3], bij[3], bji[3]);
  prelinear(aij[4], aji[4], bij[4], bji[4]);
  prevanalbada(Vi[5], ddVij[5], Vj[5], ddVji[5], aij[5], aji[5], bij[5], bji[5]);
  prevanalbada(Vi[6], ddVij[6], Vj[6], ddVji[6], aij[6], aji[6], bij[6], bji[6]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnLinearVanAlbada<6>::precomputeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				     double* daij, double* daji, double* dbij, double* dbji)
{

  prelinearDerivative(daij[0], daji[0], dbij[0], dbji[0]);
  prelinearDerivative(daij[1], daji[1], dbij[1], dbji[1]);
  prelinearDerivative(daij[2], daji[2], dbij[2], dbji[2]);
  prelinearDerivative(daij[3], daji[3], dbij[3], dbji[3]);
  prelinearDerivative(daij[4], daji[4], dbij[4], dbji[4]);
  prevanalbadaDerivative(Vi[5], dVi[5], ddVij[5], dddVij[5], Vj[5], dVj[5], ddVji[5], dddVji[5], daij[5], daji[5], dbij[5], dbji[5]);
  prevanalbadaDerivative(Vi[6], dVi[6], ddVij[6], dddVij[6], Vj[6], dVj[6], ddVji[6], dddVji[6], daij[6], daji[6], dbij[6], dbji[6]);

}

//------------------------------------------------------------------------------

template<>
inline
void RecFcnLinearVanAlbada<7>::compute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				       double* Vij, double* Vji)
{

  linear(Vi[0], ddVij[0], Vj[0], ddVji[0], Vij[0], Vji[0]);
  linear(Vi[1], ddVij[1], Vj[1], ddVji[1], Vij[1], Vji[1]);
  linear(Vi[2], ddVij[2], Vj[2], ddVji[2], Vij[2], Vji[2]);
  linear(Vi[3], ddVij[3], Vj[3], ddVji[3], Vij[3], Vji[3]);
  linear(Vi[4], ddVij[4], Vj[4], ddVji[4], Vij[4], Vji[4]);
  vanalbada(Vi[5], ddVij[5], Vj[5], ddVji[5], Vij[5], Vji[5]);
  vanalbada(Vi[6], ddVij[6], Vj[6], ddVji[6], Vij[6], Vji[6]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnLinearVanAlbada<7>::computeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				       double* dVij, double* dVji)
{

  linearDerivative(dVi[0], dddVij[0], dVj[0], dddVji[0], dVij[0], dVji[0]);
  linearDerivative(dVi[1], dddVij[1], dVj[1], dddVji[1], dVij[1], dVji[1]);
  linearDerivative(dVi[2], dddVij[2], dVj[2], dddVji[2], dVij[2], dVji[2]);
  linearDerivative(dVi[3], dddVij[3], dVj[3], dddVji[3], dVij[3], dVji[3]);
  linearDerivative(dVi[4], dddVij[4], dVj[4], dddVji[4], dVij[4], dVji[4]);
  vanalbadaDerivative(Vi[5], dVi[5], ddVij[5], dddVij[5], Vj[5], dVj[5], ddVji[5], dddVji[5], dVij[5], dVji[5]);
  vanalbadaDerivative(Vi[6], dVi[6], ddVij[6], dddVij[6], Vj[6], dVj[6], ddVji[6], dddVji[6], dVij[6], dVji[6]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnLinearVanAlbada<7>::precompute(double* Vi, double* ddVij, double* Vj, double* ddVji,
				 double* aij, double* aji, double* bij, double* bji)
{

  prelinear(aij[0], aji[0], bij[0], bji[0]);
  prelinear(aij[1], aji[1], bij[1], bji[1]);
  prelinear(aij[2], aji[2], bij[2], bji[2]);
  prelinear(aij[3], aji[3], bij[3], bji[3]);
  prelinear(aij[4], aji[4], bij[4], bji[4]);
  prevanalbada(Vi[5], ddVij[5], Vj[5], ddVji[5], aij[5], aji[5], bij[5], bji[5]);
  prevanalbada(Vi[6], ddVij[6], Vj[6], ddVji[6], aij[6], aji[6], bij[6], bji[6]);

}

//------------------------------------------------------------------------------

// Included (MB)
template<>
inline
void RecFcnLinearVanAlbada<7>::precomputeDerivative(double* Vi, double* dVi, double* ddVij, double* dddVij, double* Vj, double* dVj, double* ddVji, double* dddVji,
				     double* daij, double* daji, double* dbij, double* dbji)
{

  prelinearDerivative(daij[0], daji[0], dbij[0], dbji[0]);
  prelinearDerivative(daij[1], daji[1], dbij[1], dbji[1]);
  prelinearDerivative(daij[2], daji[2], dbij[2], dbji[2]);
  prelinearDerivative(daij[3], daji[3], dbij[3], dbji[3]);
  prelinearDerivative(daij[4], daji[4], dbij[4], dbji[4]);
  prevanalbadaDerivative(Vi[5], dVi[5], ddVij[5], dddVij[5], Vj[5], dVj[5], ddVji[5], dddVji[5], daij[5], daji[5], dbij[5], dbji[5]);
  prevanalbadaDerivative(Vi[6], dVi[6], ddVij[6], dddVij[6], Vj[6], dVj[6], ddVji[6], dddVji[6], daij[6], daji[6], dbij[6], dbji[6]);

}

//------------------------------------------------------------------------------

#endif
