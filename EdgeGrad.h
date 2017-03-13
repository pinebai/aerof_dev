#ifndef _EDGE_GRAD_H_
#define _EDGE_GRAD_H_

#include <Elem.h>
#include <Edge.h>

class IoData;
class Elem;
class ElemSet;
class SubDomain;
class Vec3D;

template<class Scalar, int dim> class SVec;
template<class Scalar> class Vec;

//------------------------------------------------------------------------------

#include "V6NodeData.h"
//------------------------------------------------------------------------------

template<int dim>
class EdgeGrad {

  double beta;
  double xiu;
  double xic;

  bool *tag;

  //V6NodeData (*v6data)[2];
  typedef V6NodeData (*V6NodeDataOf2)[2];
  V6NodeDataOf2 v6data;

private:

  void computeUpwindGradient(Elem&, double [3], SVec<double,3>&, SVec<double,dim>&, double*);

  void computeUpwindGradient(Elem&, double [3], SVec<double,3>&, SVec<double,dim>&, double*, 
			     bool*, Vec<int>&, LevelSetStructure&);
  void computeUpwindGradient(Elem&, double [3], SVec<double,3>&, SVec<double,dim>&, double*, 
			     bool*, Vec<int>&);

  void computeFaceGradient(ElemSet&, V6NodeData&, double [3], SVec<double,dim>&, SVec<double,dim>&, SVec<double,dim>&, double*);

// Included (MB)
  void computeDerivativeOfUpwindGradient(Elem&, double [3], double [3], SVec<double,3>&, SVec<double,3>&,
			     SVec<double,dim>&, SVec<double,dim>&, double*);
  void computeDerivativeOfFaceGradient(ElemSet&, V6NodeData&, double [3], double [3],
			   SVec<double,dim>&, SVec<double,dim>&, SVec<double,dim>&, SVec<double,dim>&,
			   SVec<double,dim>&, SVec<double,dim>&, double*);

public:

  EdgeGrad(IoData&);
  ~EdgeGrad();

  void settag(bool *ftag) { tag = ftag; }

  void compute(int, int, int, ElemSet&, SVec<double,3>&, SVec<double,dim>&, 
	       SVec<double,dim>&, SVec<double,dim>&, SVec<double,dim>&, double*, double*);

  //d2d
  void compute(int, int, int, ElemSet&, SVec<double,3>&, SVec<double,dim>&, 
	       SVec<double,dim>&, SVec<double,dim>&, SVec<double,dim>&, 
	       Vec<int>&, double*, double*, LevelSetStructure& LSS);
  void compute(int, int, int, ElemSet&, SVec<double,3>&, SVec<double,dim>&, 
	       SVec<double,dim>&, SVec<double,dim>&, SVec<double,dim>&, 
	       Vec<int>&, double*, double*);

  V6NodeDataOf2& getV6NodeData()  { return v6data; }


  void fix(bool (*fstag)[2], int size) {
    for (int i=0; i<size; ++i)
      tag[i] = tag[i] || fstag[i][1];
  }

  void fix(int (*fstag)[2], int size) {
    for (int i=0; i<size; ++i)
      tag[i] = tag[i] || bool(fstag[i][1]);
  }

// Included (MB)
  void computeDerivative(int, int, int, ElemSet&, SVec<double,3>&, SVec<double,3>&,
	       SVec<double,dim>&, SVec<double,dim>&, SVec<double,dim>&, SVec<double,dim>&,
	       SVec<double,dim>&, SVec<double,dim>&, SVec<double,dim>&, SVec<double,dim>&,
	       double*, double*);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <EdgeGrad.C>
#endif

#endif
