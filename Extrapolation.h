#ifndef _EXTRAPOLATION_H_
#define _EXTRAPOLATION_H_


class IoData;
class VarFcn;
class Elem;
class ElemSet;
class Vec3D;

template<class Scalar, int dim> class SVec;
//------------------------------------------------------------------------------
//	Extrapolation used: we determine which tetrahedra the line going through
//  the normal n to the boundary at node i intersects and what are the 
//	barycentric coordinates on the face opposite to i. The value of the variable
//  at that point is extrapolated
//------------------------------------------------------------------------------

struct ExtrapolationNodeData {
  int tet;
  int face;
  double r;
  double t;
  ExtrapolationNodeData() { tet = -1; face = -1; r = 0.0; t = 0.0; }
};

//------------------------------------------------------------------------------

template<int dim>
class Extrapolation {

  typedef ExtrapolationNodeData (*Extrapolationdata)[2];

  int type;

  FluidModelData::Fluid fluid;
  VarFcn *vf;

  double rho;
  double gravity;
  double ngravity[3];
  Extrapolationdata extrapolationdata;

			
public:

  Extrapolation(IoData&, VarFcn*);
  ~Extrapolation();
  
  void removeHydroStaticContribution(int i, int n[3], SVec<double,dim> &Ufar,SVec<double,dim> &V,
					 SVec<double,3> &X, double *Vinter, int node);
  void computeFaceInterpolation(int , bool &, int, ElemSet &, SVec<double,dim> &,  
				SVec<double,dim> &, double*,
				double*, int*, int*, SVec<double,3>&);
			
  Extrapolationdata &getExtrapolationData()  { return extrapolationdata; }

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <Extrapolation.C>
#endif

#endif
