#ifndef _ROETURKELJAC5_H_
#define _ROETURKELJAC5_H_

#ifdef USE_EIGEN3
#include <algorithm>
#include <cmath>
#include <limits>
#include <Eigen/Core>
#include <AutoDiff/Function.h>
#include <AutoDiff/SpaceDerivatives.h>

template<typename Scalar>
class RoeTurkelFlux5Function : public VectorValuedFunction<5,5,Scalar,18,1,double>
{
    const Eigen::Array<double,18,1> &sconst;
    const Eigen::Array<int,1,1> &iconst;

  public:
    RoeTurkelFlux5Function(const Eigen::Array<double,18,1>& _sconst, const Eigen::Array<int,1,1>& _iconst)
      : sconst(_sconst), iconst(_iconst) {}

    Eigen::Matrix<Scalar,5,1> operator() (const Eigen::Matrix<Scalar,5,1>& ug, Scalar)
    {
      const double &gamma = sconst[0];
      const double &gam = sconst[1];
      const double &pstiff = sconst[2];
      const double *enormal = &sconst[3];
      const double &evitno = sconst[6];
      const double &mach = sconst[7];
      const double &k1 = sconst[8];
      const double &cmach = sconst[9];
      const double &shockreducer = sconst[10];
      const double &irey = sconst[11];
      const double &length = sconst[12];
      const double *ud = &sconst[13];

      const int &prec = iconst[0];

      Eigen::Matrix<Scalar,5,1> phi;

      // System generated locals
      Scalar d__1;

      // Builtin functions
      using std::sqrt;
      using std::abs;

      // Local variables 
      double gam1, invgam1, rnorm, invnorm, normal[3], vitno, vitd2, squsr2, ener2;
      Scalar h, invtbeta2, r, s, t, cr, cr2, vp1, vp4, vp5, qir,
         dif1, dif2, dif3, dif4, dif5, uar1, uar2, uar3, uar4, uar5, beta,
         usro, beta2, ener1, vitg2, flur1, flur2, flur3, flur4, flur5, shock, vdotn,
         cr2byt, squsr1, vdotnt, locmach, tet1, tet2, tet3;

      // -----------------------------------------------------------------------
      // This routine computes the Flux of Roe taken at the vectors Ug, Ud
      // normal is the normal of the boundary concerned by the flux.
      // phi stores the resulting flux.
      // gamma is the dissipation coefficient
      // gam is the ratio of cp/cv
      // The Roe-Turkel Preconditioning is applied for LowMach Simulations
      // -----------------------------------------------------------------------

      // Initialization

      gam1 = gam - 1.;

      rnorm = sqrt(enormal[0] * enormal[0] + enormal[1] * enormal[1] + enormal[2] * enormal[2]);
      invnorm = 1. / rnorm;

      normal[0] = enormal[0] * invnorm;
      normal[1] = enormal[1] * invnorm;
      normal[2] = enormal[2] * invnorm;
      vitno = evitno * invnorm;

      // Computation of the centered terms

      vdotn = ug[1] * normal[0] + ug[2] * normal[1] + ug[3] * normal[2];
      vitg2 = ug[1] * ug[1] + ug[2] * ug[2] + ug[3] * ug[3];
      h = gam * (ug[4] + pstiff) + gam1 * .5 * ug[0] * vitg2;
      h /= gam1 * ug[0];
      phi[0] = ug[0] * (vdotn - vitno);
      phi[1] = phi[0] * ug[1] + ug[4] * normal[0];
      phi[2] = phi[0] * ug[2] + ug[4] * normal[1];
      phi[3] = phi[0] * ug[3] + ug[4] * normal[2];
      phi[4] = phi[0] * h + ug[4] * vitno;

      vdotn = ud[1] * normal[0] + ud[2] * normal[1] + ud[3] * normal[2];
      vdotn -= vitno;
      vitd2 = ud[1] * ud[1] + ud[2] * ud[2] + ud[3] * ud[3];
      h = gam * (ud[4] + pstiff) + gam1 * .5 * ud[0] * vitd2;
      h /= gam1 * ud[0];
      phi[0] += ud[0] * vdotn;
      phi[1] = phi[1] + ud[0] * ud[1] * vdotn + ud[4] * normal[0];
      phi[2] = phi[2] + ud[0] * ud[2] * vdotn + ud[4] * normal[1];
      phi[3] = phi[3] + ud[0] * ud[3] * vdotn + ud[4] * normal[2];
      phi[4] = phi[4] + ud[0] * vdotn * h + ud[4] * vitno;

      // Computation of the Roe-averaged state

      squsr1 = sqrt(ug[0]);
      squsr2 = sqrt(ud[0]);

      ener1 = (ug[4] + gam * pstiff) / gam1 + ug[0] * .5 * vitg2;

      ener2 = (ud[4] + gam * pstiff) / gam1 + ud[0] * .5 * vitd2;

      usro = 1. / (squsr1 + squsr2);

      uar1 = (squsr1 * ug[0] + squsr2 * ud[0]) * usro;

      uar2 = (squsr1 * ug[1] + squsr2 * ud[1]) * usro;

      uar3 = (squsr1 * ug[2] + squsr2 * ud[2]) * usro;

      uar4 = (squsr1 * ug[3] + squsr2 * ud[3]) * usro;

      uar5 = ((ener1 + ug[4]) / squsr1 + (ener2 + ud[4]) / squsr2) * usro;

      // Computation of the dissipation term
      // if prec = 1 or 2 then the dissipation is preconditioned
      // else if prec = 0 then it is not preconditioned

      // Reference: Implicit Upwind Schemes for Lowmach number Compressible Flows
      //            By Cecile Viozat (INRIA Publication)
      vdotn = normal[0] * uar2 + normal[1] * uar3 + normal[2] * uar4;

      qir = (uar2 * uar2 + uar3 * uar3 + uar4 * uar4) * .5;

      tet1 = normal[2] * uar3 - normal[1] * uar4;
      tet2 = normal[0] * uar4 - normal[2] * uar2;
      tet3 = normal[1] * uar2 - normal[0] * uar3;

      cr2 = gam1 * (uar5 - qir);
      cr = sqrt(cr2);
      cr2 = 1. / cr2;

      if (prec == 0) {
          beta = 1.;
      } else {
          shock = abs(ug[4] - ud[4]) / (ug[4] + ud[4]) / length;
          locmach = sqrt(qir * 2. * cr2 + std::numeric_limits<double>::epsilon());
          beta = std::max(Scalar(k1 * locmach), Scalar(mach));
          beta = (sqrt(irey) + 1.) * beta + shockreducer * shock;
          beta = std::min(beta,Scalar(cmach));
      }

      beta2 = beta * beta;
      dif1 = -ug[0] + ud[0];
      dif2 = -ug[0] * ug[1] + ud[0] * ud[1];
      dif3 = -ug[0] * ug[2] + ud[0] * ud[2];
      dif4 = -ug[0] * ug[3] + ud[0] * ud[3];
      dif5 = -ener1 + ener2;
      vp1 = vdotn;
      d__1 = (1. - beta2) * vdotn;
      vp4 = ((beta2 + 1.) * vdotn + sqrt(d__1 * d__1 + beta2 * 4. * (cr * cr))) * .5;
      vp5 = ((beta2 + 1.) * vdotn - sqrt(d__1 * d__1 + beta2 * 4. * (cr * cr))) * .5;

      // Roe-Turkel coefficients

      r = vp4 - vp1 * beta2;
      s = vp5 - vp1 * beta2;
      t = (vp5 - vp4) * .5;

      // Dynamic mesh inclusion

      vp1 -= vitno;
      vp4 -= vitno;
      vp5 -= vitno;

      flur1 = abs(vp1) * ((normal[0] * (1. - gam1 * qir * cr2) - tet1) * dif1 + 
              normal[0] * gam1 * uar2 * cr2 * dif2 + (normal[2] + normal[0] * 
              gam1 * uar3 * cr2) * dif3 + (-normal[1] + normal[0] * gam1 * uar4 
              * cr2) * dif4 - normal[0] * gam1 * cr2 * dif5);

      flur2 = abs(vp1) * ((normal[1] * (1. - gam1 * qir * cr2) - tet2) * dif1 + 
              (-normal[2] + normal[1] * gam1 * uar2 * cr2) * dif2 + normal[1] * 
              gam1 * uar3 * cr2 * dif3 + (normal[0] + normal[1] * gam1 * uar4 * 
              cr2) * dif4 - normal[1] * gam1 * cr2 * dif5);

      flur3 = abs(vp1) * ((normal[2] * (1. - gam1 * qir * cr2) - tet3) * dif1 + 
              (normal[1] + normal[2] * gam1 * uar2 * cr2) * dif2 + (-normal[0] 
              + normal[2] * gam1 * uar3 * cr2) * dif3 + normal[2] * gam1 * uar4 
              * cr2 * dif4 - normal[2] * gam1 * cr2 * dif5);

      flur4 = abs(vp4) * ((cr * cr * vdotn / t + gam1 * qir * s / (t * 
              beta2)) * dif1 - (cr * cr * normal[0] / t + gam1 * uar2 * s / 
              (t * beta2)) * dif2 - (cr * cr * normal[1] / t + gam1 * uar3 *
               s / (t * beta2)) * dif3 - (cr * cr * normal[2] / t + gam1 * 
              uar4 * s / (t * beta2)) * dif4 + gam1 * s / (t * beta2) * dif5);

      flur5 = abs(vp5) * (-(cr * cr * vdotn / t + gam1 * qir * r / (beta2 
              * t)) * dif1 + (cr * cr * normal[0] / t + gam1 * uar2 * r / 
              (beta2 * t)) * dif2 + (cr * cr * normal[1] / t + gam1 * uar3 *
               r / (beta2 * t)) * dif3 + (cr * cr * normal[2] / t + gam1 *
               uar4 * r / (beta2 * t)) * dif4 - gam1 * r / (beta2 * t) * 
              dif5);

      // Final phi including the numerical viscosity parameter

      phi[0] -= gamma * (normal[0] * flur1 + normal[1] * flur2 + normal[2] * 
              flur3 + (flur4 + flur5) * .5 * cr2);

      phi[1] -= gamma * (uar2 * normal[0] * flur1 + (uar2 * normal[1] - normal[
              2]) * flur2 + (uar2 * normal[2] + normal[1]) * flur3 + (normal[0] 
              * (r * flur4 + s * flur5) + uar2 * (flur4 + flur5)) * .5 * cr2);

      phi[2] -= gamma * ((uar3 * normal[0] + normal[2]) * flur1 + uar3 * 
              normal[1] * flur2 + (uar3 * normal[2] - normal[0]) * flur3 + (
              normal[1] * (r * flur4 + s * flur5) + uar3 * (flur4 + flur5)) * 
              .5 * cr2);

      phi[3] -= gamma * ((uar4 * normal[0] - normal[1]) * flur1 + (uar4 * 
              normal[1] + normal[0]) * flur2 + uar4 * normal[2] * flur3 + (
              normal[2] * (r * flur4 + s * flur5) + uar4 * (flur4 + flur5)) * 
              .5 * cr2);

      phi[4] -= gamma * ((qir * normal[0] + tet1) * flur1 + (qir * normal[1] + 
              tet2) * flur2 + (qir * normal[2] + tet3) * flur3 + (vdotn * (r *
               flur4 + s * flur5) + uar5 * (flur4 + flur5)) * .5 * cr2);

      return (0.5*rnorm)*phi;
    }
};
#endif

void roeturkeljac5(int type, double gamma, double gam, double pstiff, double enormal[3],
                   double evitno, double *ug, double *ud, double *jac, double mach, double k1,
                   double cmach, double shockreducer, double irey, double length, int prec)
{
#ifdef USE_EIGEN3
  Eigen::Array<int,1,1> iconst;
  iconst << prec;

  Eigen::Array<double,18,1> sconst;
  sconst << gamma, gam, pstiff, enormal[0], enormal[1], enormal[2], evitno, mach, k1,
            cmach, shockreducer, irey, length, ud[0], ud[1], ud[2], ud[3], ud[4];

  Eigen::Matrix<double,5,1> v = Eigen::Map<Eigen::Matrix<double,5,1> >(ug);
  Jacobian<double,RoeTurkelFlux5Function> df5dv(sconst,iconst);
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> > J(jac,5+type,5+type);
  switch(type) {
    case 0 :
      J = df5dv(v,0.);
      break;
    case 1 :
      J.topLeftCorner(5,5) = df5dv(v,0.);
      J.bottomLeftCorner(1,5).setZero();
      J.topRightCorner(5,1).setZero();
      J.bottomRightCorner(1,1).setZero();
      break;
    case 2 :
      J.topLeftCorner(5,5) = df5dv(v,0.);
      J.bottomLeftCorner(2,5).setZero();
      J.topRightCorner(5,2).setZero();
      J.bottomRightCorner(2,2).setZero();
      break;
  }
#else
  std::cerr << "Error: roeturkeljac5 function requires Eigen 3 library.\n";
  exit(-1);
#endif
}

#endif 
