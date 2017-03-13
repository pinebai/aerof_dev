#ifndef _DIST_NODAL_GRAD_H_
#define _DIST_NODAL_GRAD_H_

#include <IoData.h>

class RecFcn;
class Domain;
class DistLevelSetStructure;

#ifndef _NDGRAD_TMPL_
#define _NDGRAD_TMPL_
template<int dim, class Scalar = double> class NodalGrad;
#endif
template<class Scalar> class DistVec;
template<class Scalar, int dim> class DistSVec;
template<int dim> struct dRdXoperators;

#ifndef _DNDGRAD_TMPL_
#define _DNDGRAD_TMPL_
template<int dim, class Scalar = double> class DistNodalGrad;
#endif

//------------------------------------------------------------------------------

template<int dim, class Scalar>
class DistNodalGrad {

  SchemeData::Gradient typeGradient;

  int failSafeNewton;

  int lastConfig;

  int numLocSub;
  
  int iteration;

  IoData* myIoData;

  DistVec<bool> *tag;
  DistVec<bool> *backuptag;

  DistSVec<Scalar,dim> *Vmin;
  DistSVec<Scalar,dim> *Vmax;
  DistSVec<Scalar,dim> *phi;

  DistSVec<Scalar,3>* sensor;
  DistVec<Scalar>* sigma;

  DistSVec<double,6> *R;

  DistSVec<double,3> *wii;
  DistSVec<double,3> *wij;
  DistSVec<double,3> *wji;

  DistSVec<Scalar,dim> *ddx;
  DistSVec<Scalar,dim> *ddy;
  DistSVec<Scalar,dim> *ddz;

  DistVec<Scalar> *dTdx;
  DistVec<Scalar> *dTdy;
  DistVec<Scalar> *dTdz;

  Domain *domain;

  NodalGrad<dim, Scalar> **subNodalGrad;

// Included (MB)
  int lastConfigSA;
  DistSVec<double,dim> *dVmin;
  DistSVec<double,dim> *dVmax;
  DistSVec<double,dim> *dphi;
  DistSVec<double,6> *dR;
  DistSVec<double,3> *dwii;
  DistSVec<double,3> *dwij;
  DistSVec<double,3> *dwji;
  DistSVec<Scalar,dim> *dddx;
  DistSVec<Scalar,dim> *dddy;
  DistSVec<Scalar,dim> *dddz;

  double spheres[SchemeFixData::num * 2][4];
  double boxes[SchemeFixData::num * 2][2][3];
  double cones[SchemeFixData::num * 2][2][4];

public:

  DistNodalGrad(IoData &, Domain *);
  DistNodalGrad(IoData &, Domain *, int);
  ~DistNodalGrad();

  NodalGrad<dim, Scalar> &operator() (int i) const { return *subNodalGrad[i]; }

  DistSVec<Scalar,dim>& getX() { return *ddx; }
  DistSVec<Scalar,dim>& getY() { return *ddy; }
  DistSVec<Scalar,dim>& getZ() { return *ddz; }
  DistSVec<double,6>& getR() { return *R; }

  void updateFixes();

  void computeWeights(DistSVec<double,3> &);

  template<class Scalar2>
  void compute(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<Scalar2,dim> &);
  
  // $dd embedded
  template<class Scalar2>
  void compute(int, DistSVec<double,3> &, DistVec<double> &,
               DistVec<int> &, DistSVec<Scalar2,dim> &, bool linFSI = true, DistLevelSetStructure* =0,
               bool includeSweptNodes = true);

  template<class Scalar2>
  void computeTemperatureGradient(int, DistSVec<double,3> &, DistVec<double> &,
               DistVec<int> &, DistVec<Scalar2> &, DistLevelSetStructure* =0);

  template<class Scalar2>
  void compute(int, DistSVec<double,3> &, DistVec<double> &, DistVec<int> &, 
	       DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &,
	       DistVec<int> &, DistVec<int> &, bool linFSI = true, DistLevelSetStructure* =0);

  void compute(int config, DistSVec<double,3> &X, DistSVec<double,dim> &Psi);

  template<class Scalar2>
  void computeT(int, DistSVec<double,3> &, DistVec<double> &, DistSVec<Scalar2,dim> &,
		DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);

  template<class Scalar2>
	  void limit(RecFcn *, DistSVec<double,3> &, DistVec<double> &, 
					 DistSVec<Scalar2,dim> &,
             DistVec<int>* additionalFirstOrderNodes = 0);

  template<class Scalar2>
	  void limit(RecFcn *recFcn, DistSVec<double,3> &X,
					 DistVec<double> &ctrlVol,
					 DistLevelSetStructure *distLSS,  
					 DistSVec<Scalar2,dim> &V);

  void fix(DistSVec<bool,2>&);
  void fix(DistSVec<int,2>&);

  void resetTag();

// Included (MB)
  void computeDerivativeOfWeights(DistSVec<double,3> &, DistSVec<double,3> &);
  void computeDerivativeOfWeights(dRdXoperators<dim> &, DistSVec<double,3> &, DistSVec<double,6> &);
  void computeTransposeDerivativeOfWeights(dRdXoperators<dim> &, DistSVec<double,6> &, DistSVec<double,3> &);

  template<class Scalar2>
  void computeDerivative(int, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double> &, DistVec<double> &, DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);

  void computeDerivative(dRdXoperators<dim> *, DistSVec<double,3> &, DistVec<double> &, DistSVec<double,dim> &, 
                         DistSVec<double,6> &, DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &);

  void computeTransposeDerivative(dRdXoperators<dim> *,  
                                  DistSVec<double,dim> &, DistSVec<double,dim> &, DistSVec<double,dim> &, 
                                  DistSVec<double,6> &, DistVec<double> &, DistSVec<double,dim> &, DistSVec<double,3> &);

// Included (YC)
  void computeDerivativeOfWeightsOperators(DistSVec<double,3> &, dRdXoperators<dim> &);

  template<class Scalar2>
  void computeDerivativeOperators(DistSVec<double,3> &, DistVec<double> &, DistSVec<Scalar2,dim> &, dRdXoperators<dim> &);

  template<class Scalar2>
  void limitDerivative(RecFcn *, DistSVec<double,3> &, DistSVec<double,3> &, DistVec<double> &, DistVec<double> &, DistSVec<Scalar2,dim> &, DistSVec<Scalar2,dim> &);

  DistSVec<Scalar,dim> &getXderivative() const { return *dddx; }
  DistSVec<Scalar,dim> &getYderivative() const { return *dddy; }
  DistSVec<Scalar,dim> &getZderivative() const { return *dddz; }
  
};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <DistNodalGrad.C>
#endif

#endif
