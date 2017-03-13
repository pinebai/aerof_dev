#include <cstdio>
#include <cmath>
#include <complex>
typedef std::complex<double> bcomp;


//taken from numerical recipes in c++
//to find roots of a polynomial using Laguerre's method

//------------------------------------------------------------------------------
int laguer(bcomp *a, int degree, bcomp &x, int &its)
// Given the degree+1 complex coefficients a[0...degree] of the polynomial simu(i=0,...degree, a[i]x^i),
// and given a complex value x, this routine improves x by Laguerre's method until it converges,
// within the achievable roundoff limit, to a root of the given polynomial. The number of iterations
// taken is returned as its.
{

  const int MR=8, MT=10, MAXIT=MT*MR;
  const double EPS=1.0e-7;
  static const double frac[MR+1] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
  int iter, j;
  double abx,abp, abm, err;
  bcomp dx, x1, b, d, f, g, h, sq, gp, gm, g2;

  int m=degree;
  for (iter=1;iter<=MAXIT;iter++){
    its = iter;
    b = a[m];
    err = abs(b);
    d = f = 0.0;
    abx = abs(x);
    for (j=m-1;j>=0;j--){
      f = x*f+d;
      d = x*d+b;
      b = x*b+a[j];
      err = abs(b)+abx*err;
    }
    err *= EPS;
    if(abs(b) <= err) return 0;
    g = d/b;
    g2 = g*g;
    h = g2-2.0*f/b;
    sq = sqrt(double(m-1)*(double(m)*h-g2));
    gp = g+sq;
    gm = g-sq;
    abp = abs(gp);
    abm = abs(gm);
    if(abp<abm) gp = gm;
    dx = fmax(abp,abm) > 0.0 ? double(m)/gp : bcomp((1+abx)*cos(double(iter)),(1+abx)*sin(double(iter)));
    x1 = x-dx;
    if(x==x1) return 0;
    if(iter%MT != 0) x=x1;
    else x -= frac[iter/MT]*dx;
  }
  fprintf(stdout,"Too many iterations (%d) in laguer\n", iter);
  return 1;
}
//------------------------------------------------------------------------------
int zroots(bcomp *a, int degree, bcomp *roots, const bool &polish)
// Given the m+1 complex coefficients a[0...m] of the polynomial sum(i=0,...m, a[i]*x^i),
// this routine successively calls laguer and finds all m complex roots in roots[0...m-1].
// The boolean variable polish should be input as true if polishing (also by Laguerre's method)
// is desired, false if the roots will be subsequently polished by other means.
{

  int err = 0;
  const double EPS = 1.0e-14;
  int i,its,j,jj;
  bcomp x,b,c;

  int m=degree;
  bcomp *ad = new bcomp[m+1];
  for (j=0;j<=m;j++) ad[j]=a[j];
  for (j=m-1;j>=0;j--){
    x = 0.0;
    bcomp *ad_v = new bcomp[j+2];
    for (jj=0;jj<j+2;jj++) ad_v[jj] = ad[jj];
    err = laguer(ad_v,j+1,x,its);
    if(err>0) return 1;
    if(fabs(imag(x)) <= 2.0*EPS*fabs(real(x)))
      x = bcomp(real(x),0.0);
    roots[j] = x;
    b = ad[j+1];
    for (jj=j;jj>=0;jj--){
      c=ad[jj];
      ad[jj]=b;
      b = x*b+c;
    }
  }

  if(polish)
    for(j=0;j<m;j++){
      err = laguer(a,degree,roots[j],its);
      if(err>0) return 1;
    }
  for(j=1;j<m;j++){
    x=roots[j];
    for(i=j-1;i>=0;i--){
      if(real(roots[i]) <= real(x)) break;
      roots[i+1] = roots[i];
    }
    roots[i+1] = x;
  }
  return 0;
}
//------------------------------------------------------------------------------
