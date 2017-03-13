#ifndef _NODAL_GRAD_H_
#define _NODAL_GRAD_H_

#include <Vector.h>


#ifndef _NDGRAD_TMPL_
#define _NDGRAD_TMPL_
template<int dim, class Scalar = double> class NodalGrad;
#endif

//------------------------------------------------------------------------------

template<int dim, class Scalar>
class NodalGrad {

  SVec<Scalar,dim> &ddx;
  SVec<Scalar,dim> &ddy;
  SVec<Scalar,dim> &ddz;

  SVec<double,3> *wii;
  SVec<double,3> *wij;
  SVec<double,3> *wji;

  Vec<Scalar> &dTdx;
  Vec<Scalar> &dTdy;
  Vec<Scalar> &dTdz;

// Included
  SVec<Scalar,dim> *dddx;
  SVec<Scalar,dim> *dddy;
  SVec<Scalar,dim> *dddz;

public:

// Included
  NodalGrad(SVec<Scalar,dim> &dx, SVec<Scalar,dim> &dy, SVec<Scalar,dim> &dz,
            Vec<Scalar> &DTdx, Vec<Scalar> &DTdy, Vec<Scalar> &DTdz,
            SVec<Scalar,dim> &dX, SVec<Scalar,dim> &dY, SVec<Scalar,dim> &dZ) :
            ddx(dx), ddy(dy), ddz(dz), dTdx(DTdx), dTdy(DTdy), dTdz(DTdz)
            { dddx = &dX; dddy = &dY; dddz = &dZ; }

  NodalGrad(SVec<Scalar,dim> &dx, SVec<Scalar,dim> &dy, SVec<Scalar,dim> &dz,
            Vec<Scalar> &DTdx, Vec<Scalar> &DTdy, Vec<Scalar> &DTdz,
            SVec<Scalar,dim> &dX, SVec<Scalar,dim> &dY, SVec<Scalar,dim> &dZ,
            SVec<double,3> &_wii, SVec<double,3> &_wij, SVec<double,3> &_wji) :
            ddx(dx), ddy(dy), ddz(dz), dTdx(DTdx), dTdy(DTdy), dTdz(DTdz), wii(&_wii), wij(&_wij), wji(&_wji)
            { dddx = &dX; dddy = &dY; dddz = &dZ; }

  NodalGrad(SVec<Scalar,dim> &dx, SVec<Scalar,dim> &dy, SVec<Scalar,dim> &dz, 
            Vec<Scalar> &DTdx, Vec<Scalar> &DTdy, Vec<Scalar> &DTdz) : 
    ddx(dx), ddy(dy), ddz(dz) , dTdx(DTdx), dTdy(DTdy), dTdz(DTdz) {}
  ~NodalGrad() {}

  SVec<Scalar,dim> &getX() const { return ddx; }
  SVec<Scalar,dim> &getY() const { return ddy; }
  SVec<Scalar,dim> &getZ() const { return ddz; }
  SVec<double,3> getWii() const { return *wii; }
  SVec<double,3> getWij() const { return *wij; }
  SVec<double,3> getWji() const { return *wji; }

  Vec<Scalar> &getTX() const { return dTdx; }
  Vec<Scalar> &getTY() const { return dTdy; }
  Vec<Scalar> &getTZ() const { return dTdz; }

// Included
  SVec<Scalar,dim> &getXderivative() const { return *dddx; }
  SVec<Scalar,dim> &getYderivative() const { return *dddy; }
  SVec<Scalar,dim> &getZderivative() const { return *dddz; }

};

//------------------------------------------------------------------------------

#endif
