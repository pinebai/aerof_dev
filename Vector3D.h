#ifndef _VECTOR3D_H_
#define _VECTOR3D_H_

#include <cstdio>
#include <cmath>
#include <algorithm>

//------------------------------------------------------------------------------

struct Vec3D {

  double v[3];

  Vec3D() {v[0] = v[1] = v[2] = 0.0; }
  Vec3D(double x[3]) { v[0] = x[0]; v[1] = x[1]; v[2] = x[2]; }
  Vec3D(double x, double y, double z) { v[0] = x; v[1] = y; v[2] = z; }
  Vec3D(const Vec3D &v2) { v[0] = v2.v[0]; v[1] = v2.v[1]; v[2] = v2.v[2]; }
  Vec3D(double x) { v[0] = v[1] = v[2] = x; }
  ~Vec3D() {}

  Vec3D &operator=(const double);
  Vec3D &operator=(const Vec3D &);
  Vec3D &operator+=(const Vec3D &);
  Vec3D &operator+=(const double &);
  Vec3D &operator-=(const Vec3D &);
  Vec3D &operator-=(const double &);
  Vec3D &operator*=(double);
  Vec3D &operator/=(double);
  Vec3D operator/(double) const;

  Vec3D operator+(const Vec3D &) const;
  Vec3D operator-(const Vec3D &) const;
  Vec3D operator-() const;
  Vec3D operator^(const Vec3D &) const;

  double operator*(const Vec3D &) const;

  operator double*() { return v; }

  double &operator[](int i) { return v[i]; }
  double operator[](int i) const { return v[i]; }

  void print(const char *msg = "") { fprintf(stdout, "%s(%e %e %e)\n", msg, v[0], v[1], v[2]); }

  double norm() { return(sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])); }
  double normsq() { return(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); }

  void crossProductSensitivityOperator(double dReturndv2[3][3]);
  void crossProductSensitivityOperator(const Vec3D &v2, double dReturndv1[3][3]);

};

//------------------------------------------------------------------------------

inline 
Vec3D &Vec3D::operator=(const double v2)
{

  v[0] = v2;
  v[1] = v2;
  v[2] = v2;

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec3D &Vec3D::operator=(const Vec3D &v2)
{

  v[0] = v2.v[0];
  v[1] = v2.v[1];
  v[2] = v2.v[2];

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec3D &Vec3D::operator+=(const Vec3D &v2)
{

  v[0] += v2.v[0];
  v[1] += v2.v[1];
  v[2] += v2.v[2];

  return *this;

}

//------------------------------------------------------------------------------

inline
Vec3D &Vec3D::operator+=(const double &c)
{

  v[0] += c;
  v[1] += c;
  v[2] += c;

  return *this;

}

//------------------------------------------------------------------------------

inline 
Vec3D &Vec3D::operator-=(const Vec3D &v2)
{

  v[0] -= v2.v[0];
  v[1] -= v2.v[1];
  v[2] -= v2.v[2];

  return *this;

}

//------------------------------------------------------------------------------

inline
Vec3D &Vec3D::operator-=(const double &c)
{

  v[0] -= c;
  v[1] -= c;
  v[2] -= c;

  return *this;

}

//------------------------------------------------------------------------------
// define vector cross product sensitivity operator with respect to v2

inline
void Vec3D::crossProductSensitivityOperator(double dReturndv2[3][3])
{

  dReturndv2[0][0] = 0.0;
  dReturndv2[0][1] = -v[2];
  dReturndv2[0][2] = v[1];
  dReturndv2[1][0] = v[2];
  dReturndv2[1][1] = 0.0;
  dReturndv2[1][2] = -v[0];
  dReturndv2[2][0] = -v[1];
  dReturndv2[2][1] = v[0];
  dReturndv2[2][2] = 0.0;

}

//------------------------------------------------------------------------------
// define vector cross product sensitivity operator with respect to v1

inline
void Vec3D::crossProductSensitivityOperator(const Vec3D &v2, double dReturndv1[3][3])
{

  dReturndv1[0][0] = 0.0;
  dReturndv1[0][1] = v2.v[2];
  dReturndv1[0][2] = -v2.v[1];
  dReturndv1[1][0] = -v2.v[2];
  dReturndv1[1][1] = 0.0;
  dReturndv1[1][2] = v2.v[0];
  dReturndv1[2][0] = v2.v[1];
  dReturndv1[2][1] = -v2.v[0];
  dReturndv1[2][2] = 0.0;

}

//------------------------------------------------------------------------------

inline 
Vec3D &Vec3D::operator*=(double cst)
{

  v[0] *= cst;
  v[1] *= cst;
  v[2] *= cst;

  return *this;

}

//------------------------------------------------------------------------------

inline
Vec3D &Vec3D::operator/=(double cst)
{
  cst = 1.0/cst;

  v[0] *= cst;
  v[1] *= cst;
  v[2] *= cst;

  return *this;
}

//------------------------------------------------------------------------------

inline
Vec3D Vec3D::operator/(double cst)  const
{
  Vec3D vnew;
  vnew[0] = v[0] / cst;
  vnew[1] = v[1] / cst;
  vnew[2] = v[2] / cst;

  return vnew;
}

//------------------------------------------------------------------------------

inline 
Vec3D Vec3D::operator+(const Vec3D &v2) const
{

  Vec3D res;

  res.v[0] = v[0]+v2.v[0];
  res.v[1] = v[1]+v2.v[1];
  res.v[2] = v[2]+v2.v[2];

  return res;

}

//------------------------------------------------------------------------------

inline 
Vec3D Vec3D::operator-(const Vec3D &v2) const
{

  Vec3D res;

  res.v[0] = v[0]-v2.v[0];
  res.v[1] = v[1]-v2.v[1];
  res.v[2] = v[2]-v2.v[2];

  return res;

}

//------------------------------------------------------------------------------

inline 
Vec3D Vec3D::operator-() const
{
  Vec3D res;

  res.v[0] = -v[0];
  res.v[1] = -v[1];
  res.v[2] = -v[2];

  return res;
}

//------------------------------------------------------------------------------
// define vector cross product

inline 
Vec3D Vec3D::operator^(const Vec3D &v2) const
{

  Vec3D res;

  res.v[0] = v[1]*v2.v[2]-v[2]*v2.v[1];
  res.v[1] = v[2]*v2.v[0]-v[0]*v2.v[2];
  res.v[2] = v[0]*v2.v[1]-v[1]*v2.v[0];

  return res;

}

//------------------------------------------------------------------------------

inline 
double Vec3D::operator*(const Vec3D &v2) const
{

  return v[0]*v2.v[0] + v[1]*v2.v[1] + v[2]*v2.v[2];

}

//------------------------------------------------------------------------------

inline 
Vec3D operator*(double c, const Vec3D &v)
{

  Vec3D res;

  res.v[0] = c*v.v[0];
  res.v[1] = c*v.v[1];
  res.v[2] = c*v.v[2];

  return res;

}

inline 
Vec3D operator*(const Vec3D &v,double c )
{

  Vec3D res;

  res.v[0] = c*v.v[0];
  res.v[1] = c*v.v[1];
  res.v[2] = c*v.v[2];

  return res;

}

//inline double min(double x,double y) { return (x<y?x:y); }
//inline double max(double x,double y) { return (x>y?x:y); }

inline
Vec3D min( const Vec3D& a, const Vec3D& b) {

  return Vec3D( std::min(a[0],b[0]),
                std::min(a[1],b[1]),
                std::min(a[2],b[2]));
}

inline
Vec3D max( const Vec3D& a, const Vec3D& b) {

  return Vec3D( std::max(a[0],b[0]),
                std::max(a[1],b[1]),
                std::max(a[2],b[2]));
}

//------------------------------------------------------------------------------

#endif
