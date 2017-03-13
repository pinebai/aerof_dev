#include <WallFcn.h>

#include <IoData.h>
#include <BcDef.h>
#include <VarFcn.h>
#include <ViscoFcn.h>
#include <Vector3D.h>

#include <cmath>

//------------------------------------------------------------------------------

const double WallFcn::third = 1.0 / 3.0;
const double WallFcn::eleventh = 1.0 / 11.0;

//------------------------------------------------------------------------------

WallFcn::WallFcn(IoData &iod, VarFcn *varf, ViscoFcn *visf) : 
  varFcn(varf), viscoFcn(visf)
{

  gam = varFcn->specificHeatCstPressure();
  vkcst = 0.41;
  reynolds = iod.ref.reynolds_mu;

// Included (MB)
  prandtl = iod.eqs.thermalCondModel.prandtl;
  dRedMach = iod.ref.dRe_mudMach;

}

//------------------------------------------------------------------------------

// Included (MB)
void WallFcn::rstVar(IoData &iod, Communicator *com)
{

  reynolds = iod.ref.reynolds_mu;
  dRedMach = iod.ref.dRe_mudMach;
  viscoFcn->rstVar(iod);

}

//------------------------------------------------------------------------------

Vec3D WallFcn::computeTangentVector(Vec3D &n, Vec3D &u)
{

  Vec3D un = (u * n) * n;

  Vec3D t = u - un;

  double norm = sqrt(t*t);

  if (norm != 0.0)
    t *= 1.0 / norm;

  return t;

}

//------------------------------------------------------------------------------

// Included (MB)
Vec3D WallFcn::computeDerivativeOfTangentVector(Vec3D &n, Vec3D &dn, Vec3D &u, Vec3D &du)
{

  Vec3D un = (u * n) * n;

  Vec3D dun = (du * n + u * dn) * n + (u * n) * dn;

  Vec3D t = u - un;

  Vec3D dt = du - dun;

  double norm = sqrt(t*t);

  double dNorm;

  if (norm != 0.0) {
    dNorm = (1.0 / norm)*(t*dt);
    dt *= 1.0 / norm;
    dt += ( -1.0 / (norm * norm) * dNorm ) * t;
  }

  return dt;

}

//------------------------------------------------------------------------------

void WallFcn::computeFaceValues(double d2wall[3], double *Vwall, double *V[3],
				double &delta, Vec3D &du, double &dT, Vec3D &uw,
				double &rhow, double &muw)
{

  delta = third * (d2wall[0] + d2wall[1] + d2wall[2]);

  rhow = third * ( varFcn->getDensity(V[0]) + 
		   varFcn->getDensity(V[1]) +
		   varFcn->getDensity(V[2]) );

  uw[0] = Vwall[1];
  uw[1] = Vwall[2];
  uw[2] = Vwall[3];

  double Tw = Vwall[4];

  du = third * ( varFcn->getVelocity(V[0]) +
		 varFcn->getVelocity(V[1]) +
		 varFcn->getVelocity(V[2]) ) - uw;

  double T = third * ( varFcn->computeTemperature(V[0]) + 
		       varFcn->computeTemperature(V[1]) +
		       varFcn->computeTemperature(V[2]) );

  dT = T - Tw;

  muw = viscoFcn->compute_mu(T);

}

//------------------------------------------------------------------------------

// Included (MB)
void WallFcn::computeDerivativeOfFaceValues(double d2wall[3], double *Vwall, double *dVwall, double *V[3], double *dV[3],
				double &delta, Vec3D &ddu, double &ddT, Vec3D &duw,
				double &drhow, double &dmuw, double &dMach)
{

  delta = third * (d2wall[0] + d2wall[1] + d2wall[2]);

  drhow = third * ( varFcn->getDensity(dV[0]) +
		    varFcn->getDensity(dV[1]) +
		    varFcn->getDensity(dV[2]) );

  duw[0] = dVwall[1];
  duw[1] = dVwall[2];
  duw[2] = dVwall[3];

  double dTw = dVwall[4];

  ddu = third * ( varFcn->getDerivativeOfVelocity(dV[0]) +
		  varFcn->getDerivativeOfVelocity(dV[1]) +
		  varFcn->getDerivativeOfVelocity(dV[2]) ) - duw;

  double T = third * ( varFcn->computeTemperature(V[0]) +
		       varFcn->computeTemperature(V[1]) +
		       varFcn->computeTemperature(V[2]) );

  double dT = third * ( varFcn->computeDerivativeOfTemperature(V[0], dV[0]) +
		        varFcn->computeDerivativeOfTemperature(V[1], dV[1]) +
		        varFcn->computeDerivativeOfTemperature(V[2], dV[2]) );

  ddT = dT - dTw;

  dmuw = viscoFcn->compute_muDerivative(T, dT, dMach);

}

//------------------------------------------------------------------------------

double WallFcn::computeFrictionVelocity(double ut, double delta, double rho, double mu)
{

  int maxits = 20;
  double eps = 1.e-6;

  //double ut = u * t;

  if (ut < 0.0) ut = 0.0; 

  double utau = sqrt(mu * ut / (reynolds * rho * delta));

  double res;
  double target;
  int it;
  for (it=0; it<maxits; ++it) {

    double dplus = reynolds * utau * delta * rho / mu;

    if (dplus < 1.0) {
      dplus = 1.0;
      utau = mu / (reynolds * delta * rho);
      break;
    }

    double f = 2.5 * log(1.0 + vkcst*dplus) + 7.8 * (1.0 - exp(-eleventh*dplus) - eleventh*dplus*exp(-0.33*dplus));

    double F = utau * f - ut;

    res = F * F;

    if (it == 0) 
      target = eps * res;

    if (res == 0.0 || res <= target) 
      break;

    double dfdutau = 2.5*vkcst/(1.0 + vkcst*dplus) + 7.8*eleventh*(exp(-eleventh*dplus) - exp(-0.33*dplus) + 0.33*dplus*exp(-0.33*dplus));

    dfdutau *= reynolds * delta * rho / mu;

    utau -= F / (f + utau * dfdutau);

  }

  if (it == maxits) 
    fprintf(stderr, "*** Warning: Newton did not converge on utau; ut = %lf; delta = %lf rho = %lf; mu = %lf; res = %lf \n", ut, delta, rho, mu, res);

  if (utau <= 0.0)
    fprintf(stderr, "*** Warning: utau=%e\n", utau);

  return utau;

}

//------------------------------------------------------------------------------

// Included (MB)
double WallFcn::computeDerivativeOfFrictionVelocity(Vec3D &t, Vec3D &dt, double delta, double rho, double drho, Vec3D &u, Vec3D &du, double mu, double dmu, double dMach)
{

  int maxits = 20;
  double eps = 1.e-6;

  double ut = u * t;

  if(ut < 0.0) ut = 0.0;

  double dut = du * t + u * dt;

  double utau = sqrt( mu * ut / ( reynolds * rho * delta ) );

  double dutau;

  double target;

  if (ut != 0.0)
    dutau = 1.0 / (2.0 * utau) * ( ( dmu * ut + mu * dut ) * ( reynolds * rho ) - ( mu * ut ) * ( dRedMach * dMach * rho + reynolds * drho ) ) / ( reynolds * reynolds * rho * rho * delta );
  else
    dutau = 0.0;

  int it;
  for (it=0; it<maxits; ++it) {
    double dplus = reynolds * utau * delta * rho / mu;
    
    double ddplus = delta * ((dRedMach * dMach * utau * rho + reynolds * dutau * rho + reynolds * utau * drho) * mu - reynolds * utau * rho * dmu) / (mu * mu);    

    if (dplus < 1.0) {
      dplus = 1.0;
      ddplus = 0.0;
      utau = mu / (reynolds * delta * rho);
      dutau = (dmu * reynolds * rho - mu * (dRedMach * dMach * rho + reynolds * drho)) / (reynolds * reynolds * rho * rho * delta);
      break;
    }

    double f = 2.5 * log(1.0 + vkcst*dplus) + 7.8 * (1.0 - exp(-eleventh*dplus) - eleventh*dplus*exp(-0.33*dplus));

    double df = 2.5 * vkcst / ( 1.0 + vkcst * dplus ) * ddplus + 7.8 * eleventh * ( exp( -eleventh * dplus ) * ddplus - ddplus * exp( -0.33*dplus ) + dplus * 0.33 * exp( -0.33*dplus ) * ddplus );

    double F = utau * f - ut;

    double dF = dutau * f + utau * df - dut;

    double res = F * F;

    if (it == 0)
      target = eps * res;

    if (res == 0.0 || res <= target)
      break;

    double dfdutau = 2.5*vkcst/(1.0 + vkcst*dplus) + 7.8*eleventh*(exp(-eleventh*dplus) - exp(-0.33*dplus) + 0.33*dplus*exp(-0.33*dplus));

    double ddfdutau = ( -2.5 * vkcst * vkcst / ( ( 1.0 + vkcst * dplus ) * ( 1.0 + vkcst * dplus ) ) * ddplus + 7.8 * eleventh * ( - eleventh * exp( -eleventh * dplus ) * ddplus + 0.33 * exp( -0.33 * dplus ) * ddplus + 0.33 * exp( -0.33 * dplus ) * ddplus - 0.33 * 0.33 * dplus * exp( -0.33 * dplus ) * ddplus ) ); 

    ddfdutau *= reynolds * delta * rho / mu;

    ddfdutau += dfdutau * delta * ( ( dRedMach * dMach * rho + reynolds * drho ) * mu - reynolds * rho * dmu)/ ( mu * mu );

    dfdutau *= reynolds * delta * rho / mu;

    dutau -= ( dF * ( f + utau * dfdutau ) - F * ( df + dutau * dfdutau + utau * ddfdutau ) ) / ( ( f + utau * dfdutau ) * ( f + utau * dfdutau ) ) ;

    utau -= F / (f + utau * dfdutau);
  }

  if (it == maxits)
    fprintf(stderr, "*** Warning: Newton did not converge on utau\n");

  if (utau <= 0.0)
    fprintf(stderr, "*** Warning: utau=%e\n", utau);

  return dutau;

}

//------------------------------------------------------------------------------

double WallFcn::computeFrictionTemperature(double utau, double delta, double rho, 
					   double dT, double mu)
{

  double dplus = reynolds * utau * delta * rho / mu;

  double Tplus;

  if (dplus < 13.2)
    Tplus = prandtl * dplus;
  else
    Tplus = 2.195*log(dplus) + 13.2*prandtl - 5.66;

  return -dT / Tplus;

}

//------------------------------------------------------------------------------

// Included (MB)
double WallFcn::computeDerivativeOfFrictionTemperature(double utau, double dutau, double delta, double rho, double drho,
					   double dT, double ddT, double mu, double dmu, double dMach)
{

  double dplus = reynolds * utau * delta * rho / mu;

  double ddplus = delta * ( ( dRedMach * dMach * utau * rho + reynolds * dutau * rho + reynolds * utau * drho ) * mu - reynolds * utau * rho * dmu ) / ( mu * mu );

  double Tplus, dTplus;

  if (dplus < 13.2) {
    Tplus = prandtl * dplus;
    dTplus = prandtl * ddplus;
  }
  else {
    Tplus = 2.195*log(dplus) + 13.2*prandtl - 5.66;
    dTplus = 2.195 / dplus * ddplus;
  }

  return ( -ddT * Tplus + dT * dTplus ) / ( Tplus * Tplus );

}

//------------------------------------------------------------------------------

void WallFcn::computeSurfaceTerm(int code, Vec3D &normal, double d2wall[3],
				 double *Vwall, double *V[3], double *term)
{

  double delta, dT, rhow, muw;

  Vec3D du, uw;

  computeFaceValues(d2wall, Vwall, V, delta, du, dT, uw, rhow, muw);

  double norm = sqrt(normal*normal);

  Vec3D n = (1.0/norm) * normal;

  Vec3D t = computeTangentVector(n, du);

  double ut = du * t;

  double utau = computeFrictionVelocity(ut, delta, rhow, muw);

  double a = - rhow * utau*utau * norm;

  term[0] = 0.0;
  term[1] = a * t[0];
  term[2] = a * t[1];
  term[3] = a * t[2];
  term[4] = a * (uw * t);

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED) {
    double Ttau = computeFrictionTemperature(utau, delta, rhow, dT, muw);
    term[4] += gam * rhow * utau * Ttau * norm;
  }

  computeWallValues(utau, delta, rhow, du*t, muw, Vwall);

}

//------------------------------------------------------------------------------

// Included (MB)
void WallFcn::computeDerivativeOfSurfaceTerm(int code, Vec3D &normal, Vec3D &dNormal, double d2wall[3],
				 double *Vwall, double *dVwall, double *V[3], double *dV[3], double dMach, double *dTerm)
{

  double delta, dT, rhow, muw;

  double ddT, drhow, dmuw;

  Vec3D du, uw;

  Vec3D ddu, duw;

  computeFaceValues(d2wall, Vwall, V, delta, du, dT, uw, rhow, muw);

  computeDerivativeOfFaceValues(d2wall, Vwall, dVwall, V, dV, delta, ddu, ddT, duw, drhow, dmuw, dMach);

  double norm = sqrt(normal*normal);

  double dNorm = (1.0/norm) * (normal*dNormal);

  Vec3D n = (1.0/norm) * normal;

  Vec3D dn = (-1.0/(norm*norm)*dNorm) * normal + (1.0/norm) * dNormal;

  Vec3D t = computeTangentVector(n, du);

  Vec3D dt = computeDerivativeOfTangentVector(n, dn, du, ddu);

  double ut = du * t;

  double utau = computeFrictionVelocity(ut, delta, rhow, muw);

  double dutau = computeDerivativeOfFrictionVelocity(t, dt, delta, rhow, drhow, du, ddu, muw, dmuw, dMach);

  double a = - rhow * utau*utau * norm;

  double da = - drhow * utau*utau * norm - rhow * 2.0*utau*dutau * norm - rhow * utau*utau * dNorm;

  dTerm[0] = 0.0;
  dTerm[1] = da * t[0] + a * dt[0];
  dTerm[2] = da * t[1] + a * dt[1];
  dTerm[3] = da * t[2] + a * dt[2];
  dTerm[4] = da * ( uw * t ) + a * ( duw * t + uw * dt );

  if (code == BC_ISOTHERMAL_WALL_MOVING || code == BC_ISOTHERMAL_WALL_FIXED) {
    double Ttau = computeFrictionTemperature(utau, delta, rhow, dT, muw);
    double dTtau = computeDerivativeOfFrictionTemperature(utau, dutau, delta, rhow, drhow, dT, ddT, muw, dmuw, dMach);
    dTerm[4] += gam * drhow * utau * Ttau * norm + gam * rhow * dutau * Ttau * norm + gam * rhow * utau * dTtau * norm + gam * rhow * utau * Ttau * dNorm;
  }

  double dut = du * t;
  double ddut = ddu * t + du * dt;

  computeDerivativeOfWallValues(utau, dutau, delta, rhow, drhow, dut, ddut, muw, dmuw, dMach, Vwall, dVwall);

}

//------------------------------------------------------------------------------

double WallFcn::computedudT(double rho, double T, double du, double dT, 
								  double d2w, double &dudn, double &dTdn)
{

	double delta = d2w; // *********

	double mu = viscoFcn->compute_mu(T);

	double utau = computeFrictionVelocity(du, delta, rho, mu);

	dudn = reynolds * rho * utau*utau / mu;

	dTdn = 0.0; // To be done *******

	//double yp = reynolds * utau * delta * rho / mu;
	//fprintf(stdout, "  %f,%f,%f\n", delta, utau, yp);

	return utau;
}
//------------------------------------------------------------------------------

Vec3D WallFcn::computeForce(Vec3D &normal, double d2wall[3], double *Vwall, double *V[3])
{

  double delta, dT, rhow, muw;

  Vec3D du, uw;

  computeFaceValues(d2wall, Vwall, V, delta, du, dT, uw, rhow, muw);

  double norm = sqrt(normal*normal);

  Vec3D n = (1.0/norm) * normal;

  Vec3D t = computeTangentVector(n, du);

	double ut = du * t;

	double utau = computeFrictionVelocity(ut, delta, rhow, muw);

  Vec3D force = - rhow * utau*utau * norm * t;

  return force;

}

//------------------------------------------------------------------------------

// Included (MB)
Vec3D WallFcn::computeDerivativeOfForce(Vec3D &normal, Vec3D &dNormal, double d2wall[3], double *Vwall, double *dVwall, double *V[3], double *dV[3], double dMach)
{

  double delta, dT, rhow, muw;

  double ddT, drhow, dmuw;

  Vec3D du, uw;

  Vec3D ddu, duw;

  computeFaceValues(d2wall, Vwall, V, delta, du, dT, uw, rhow, muw);

  computeDerivativeOfFaceValues(d2wall, Vwall, dVwall, V, dV, delta, ddu, ddT, duw, drhow, dmuw, dMach);

  double norm = sqrt(normal*normal);

  double dNorm = (1.0/norm) * normal * dNormal;

  Vec3D n = (1.0/norm) * normal;

  Vec3D dn = (-1.0/(norm*norm)*dNorm) * normal + (1.0/norm) * dNormal;

  Vec3D t = computeTangentVector(n, du);

  Vec3D dt = computeDerivativeOfTangentVector(n, dn, du, ddu);

  double ut = du * t;

  double utau = computeFrictionVelocity(ut, delta, rhow, muw);

  double dutau = computeDerivativeOfFrictionVelocity(t, dt, delta, rhow, drhow, du, ddu, muw, dmuw, dMach);

  Vec3D dForce = - drhow * utau*utau * norm * t - rhow * 2.0*utau*dutau * norm * t - rhow * utau*utau * dNorm * t - rhow * utau*utau * norm * dt;

  return dForce;

}

//------------------------------------------------------------------------------

double WallFcn::computeHeatPower(Vec3D &normal, double d2wall[3], double *Vwall, double *V[3])
{

  double delta, dT, rhow, muw;

  Vec3D du, uw;

  computeFaceValues(d2wall, Vwall, V, delta, du, dT, uw, rhow, muw);

  double norm = sqrt(normal*normal);

  Vec3D n = (1.0/norm) * normal;

  Vec3D t = computeTangentVector(n, du);

  double ut = du * t;

  double utau = computeFrictionVelocity(ut, delta, rhow, muw);

  double Ttau = computeFrictionTemperature(utau, delta, rhow, dT, muw);

  double hp = - gam * rhow * utau * Ttau * norm;

  return hp;

}

//------------------------------------------------------------------------------

// Included (MB)
double WallFcn::computeDerivativeOfHeatPower(Vec3D &normal, Vec3D &dNormal, double d2wall[3], double *Vwall, double *dVwall, double *V[3], double *dV[3], double dMach)
{

  double delta, dT, rhow, muw;

  double ddT, drhow, dmuw;

  Vec3D du, uw;

  Vec3D ddu, duw;

  computeFaceValues(d2wall, Vwall, V, delta, du, dT, uw, rhow, muw);

  computeDerivativeOfFaceValues(d2wall, Vwall, dVwall, V, dV, delta, ddu, ddT, duw, drhow, dmuw, dMach);

  double norm = sqrt(normal*normal);

  double dNorm = (1.0/norm) * normal * dNormal;

  Vec3D n = (1.0/norm) * normal;

  Vec3D dn = (-1.0/(norm*norm)*dNorm) * normal + (1.0/norm) * dNormal;

  Vec3D t = computeTangentVector(n, du);

  Vec3D dt = computeDerivativeOfTangentVector(n, dn, du, ddu);

  double ut = du * t;

  double utau = computeFrictionVelocity(ut, delta, rhow, muw);

  double dutau = computeDerivativeOfFrictionVelocity(t, dt, delta, rhow, drhow, du, ddu, muw, dmuw, dMach);

  double Ttau = computeFrictionTemperature(utau, delta, rhow, dT, muw);

  double dTtau = computeDerivativeOfFrictionTemperature(utau, dutau, delta, rhow, drhow, dT, ddT, muw, dmuw, dMach);

  double dhp = - gam * drhow * utau * Ttau * norm - gam * rhow * dutau * Ttau * norm - gam * rhow * utau * dTtau * norm - gam * rhow * utau * Ttau * dNorm;

  return dhp;

}

//------------------------------------------------------------------------------

double WallFcn::computeInterfaceWork(Vec3D& normal, double d2wall[3], double* Vwall, double* V[3])
{

  double delta, dT, rhow, muw;

  Vec3D du, uw;

  computeFaceValues(d2wall, Vwall, V, delta, du, dT, uw, rhow, muw);

  double norm = sqrt(normal*normal);

  Vec3D n = (1.0/norm) * normal;

  Vec3D t = computeTangentVector(n, du);

  double ut = du * t;

  double utau = computeFrictionVelocity(ut, delta, rhow, muw);

  double W = - rhow * utau*utau * norm * (uw * t);

  return W;

}

//------------------------------------------------------------------------------

double WallFcn::computeDeltaPlus(Vec3D &normal, double d2wall[3], double *Vwall, double *V[3])
{

  double delta, dT, rhow, muw;

  Vec3D du, uw;

  computeFaceValues(d2wall, Vwall, V, delta, du, dT, uw, rhow, muw);

  double norm = sqrt(normal*normal);

  Vec3D n = (1.0/norm) * normal;

  Vec3D t = computeTangentVector(n, du);

  double ut = du * t;

  double utau = computeFrictionVelocity(ut, delta, rhow, muw);

  return reynolds * utau * delta * rhow / muw;

}

//------------------------------------------------------------------------------

WallFcnSA::WallFcnSA(IoData &iod, VarFcn *varf, ViscoFcn *visf) : 
  WallFcn(iod, varf, visf)
{

  double cv1 = iod.eqs.tc.tm.sa.cv1;
  cv1_pow3 = cv1*cv1*cv1;

}

//------------------------------------------------------------------------------

void WallFcnSA::computeWallValues(double utau, double delta, double rho, 
				  double ut, double mu, double *V)
{

  const int maxits = 20;
  const double eps = 1.e-6;
  const double coef0 = 0.4 * exp(-0.4*5.5);

  double kuplus = 0.4 * ut / utau;

  /*
  double dplus = reynolds * utau * delta * rho / mu;
  if (dplus > 200.0)
    kuplus = 0.4 * 75.0;
  */

  double mutomu = coef0 * (exp(kuplus) - 1.0 - kuplus - 0.5*kuplus*kuplus);

  if(mutomu<0.0) {
//     fprintf(stderr,"*** Warning: mutomu clipped at wall boundary (mutomu = %e)\n",mutomu);
     mutomu = 0.0;
  }

  double mutildeomu;
  if (mutomu < 3.5)
    mutildeomu = pow(mutomu*cv1_pow3, 0.25);
  else
    mutildeomu = mutomu;

  int it;
  double res, target;
  for (it=0; it<maxits; ++it) {
    double mutildeomu2 = mutildeomu*mutildeomu;
    double mutildeomu3 = mutildeomu2*mutildeomu;
    double f = mutildeomu2*mutildeomu2 - mutomu*(mutildeomu3 + cv1_pow3);
    res = f * f;
    if (it == 0)
      target = eps * res;
    if (res == 0.0 || res <= target)
      break;
    double df = 4.0*mutildeomu3 - 3.0*mutomu*mutildeomu2;
    mutildeomu -= f / df;
  }
  if (it == maxits && target > 1.e-16) {
    fprintf(stderr, "*** Warning: Newton reached %d its on mutilde", maxits);
    fprintf(stderr, " (res=%.2e, target=%.2e)\n", res, target);
  }

  V[5] = mu * mutildeomu;
 
}

//------------------------------------------------------------------------------

// Included (MB)
void WallFcnSA::computeDerivativeOfWallValues(double utau, double dutau, double delta, double rho, double drho,
				  double ut, double dut, double mu, double dmu, double dMach, double *V, double *dV)
{

  const int maxits = 20;
  const double eps = 1.e-6;
  const double coef0 = 0.4 * exp(-0.4*5.5);

  double kuplus = 0.4 * ut / utau;

  double dkuplus = 0.4 * ( dut * utau - ut * dutau ) / ( utau * utau );

  double mutomu = coef0 * (exp(kuplus) - 1.0 - kuplus - 0.5*kuplus*kuplus);
 
  double dmutomu = coef0 * (exp(kuplus) * dkuplus - dkuplus - kuplus*dkuplus);

  if(mutomu<0.0) {
//     fprintf(stderr,"*** Warning: mutomu clipped at wall boundary (mutomu = %e)\n",mutomu);
     mutomu = 0.0;
     dmutomu = 0.0;
  }

  double mutildeomu;
  double dmutildeomu;
  if (mutomu < 3.5)  {
    mutildeomu = pow(mutomu*cv1_pow3, 0.25);
    dmutildeomu = 0.25*pow(mutomu*cv1_pow3, -0.75)*dmutomu*cv1_pow3;
  }
  else {
    mutildeomu = mutomu;
    dmutildeomu = dmutomu;
  }

  int it;
  double res, target;
  for (it=0; it<maxits; ++it) {
    double mutildeomu2 = mutildeomu*mutildeomu;
    double mutildeomu3 = mutildeomu2*mutildeomu;
    double f = mutildeomu2*mutildeomu2 - mutomu*(mutildeomu3 + cv1_pow3);
    double dF = 4.0 * mutildeomu3 * dmutildeomu - dmutomu * (mutildeomu3 + cv1_pow3) - mutomu * 3.0 * mutildeomu2 * dmutildeomu;
    res = f * f;
    if (it == 0)
      target = eps * res;
    if (res == 0.0 || res <= target)
      break;
    double df = 4.0*mutildeomu3 - 3.0*mutomu*mutildeomu2;
    double ddf = 12.0 * mutildeomu2 * dmutildeomu - 3.0 * dmutomu * mutildeomu2 - 6.0 * mutomu * mutildeomu * dmutildeomu;
    mutildeomu -= f / df;
    dmutildeomu -= ( dF * df - f * ddf ) / ( df * df );
  }
  if (it == maxits) {
    fprintf(stderr, "*** Warning: Newton reached %d its on mutilde", maxits);
    fprintf(stderr, " (res=%.2e, target=%.2e)\n", res, target);
  }

  V[5] = mu * mutildeomu;

  dV[5] = dmu * mutildeomu + mu * dmutildeomu;

}

//------------------------------------------------------------------------------

WallFcnKE::WallFcnKE(IoData &iod, VarFcn *varf, ViscoFcn *visf) : 
  WallFcn(iod, varf, visf)
{

  orcmu = 1.0 / sqrt(iod.eqs.tc.tm.ke.c_mu);

}

//------------------------------------------------------------------------------

void WallFcnKE::computeWallValues(double utau, double delta, double rho, 
				  double ut, double mu, double *V)
{

  double nu = mu / rho;
  double dplus = reynolds * utau * delta / nu;
  double dplus2 = dplus * dplus;
  double utau2 = utau * utau;

  double k, eps;

  if (dplus < 10.0) {
    k = orcmu * utau2 * 0.01 * dplus2;
    eps = reynolds * utau2 * utau2 * 0.1 / (vkcst * nu);
    eps *= 0.01 * dplus2 + 0.2*vkcst*orcmu * (1.0 - 0.01 * dplus2);
  }
  else {
    k = orcmu * utau2;
    eps = utau2 * utau / (vkcst * delta);
  }

  V[5] = k;
  V[6] = eps;
 
}

//------------------------------------------------------------------------------

// Included (MB)
void WallFcnKE::computeDerivativeOfWallValues(double utau, double dutau, double delta, double rho, double drho,
				  double ut, double dut, double mu, double dmu, double dMach, double *V, double *dV)
{

  double nu = mu / rho;
  double dnu = ( dmu * rho - mu * drho ) / ( rho * rho );
  double dplus = reynolds * utau * delta / nu;
  double ddplus = delta * ( ( dRedMach * dMach * utau + reynolds * dutau ) * nu - reynolds * utau * dnu ) / ( nu * nu );

  double dplus2 = dplus * dplus;
  double utau2 = utau * utau;

  double k, eps;
  double dk, deps;

  if (dplus < 10.0) {
    k = orcmu * utau2 * 0.01 * dplus2;
    dk = orcmu * utau * dutau * 0.02 * dplus2 + orcmu * utau2 * 0.02 * dplus * ddplus;
    eps = reynolds * utau2 * utau2 * 0.1 / (vkcst * nu);
    deps = (0.1/vkcst) * ( ( dRedMach * dMach * utau2 * utau2 + reynolds * 4.0 * utau * utau2 * dutau ) * nu - reynolds * utau2 * utau2 * dnu ) / ( nu * nu );
    deps *= 0.01 * dplus2 + 0.2*vkcst*orcmu * (1.0 - 0.01 * dplus2);
    deps += eps * ( 0.02 * dplus * ddplus - 0.004 * vkcst * orcmu * dplus * ddplus );
    eps *= 0.01 * dplus2 + 0.2*vkcst*orcmu * (1.0 - 0.01 * dplus2);
  }
  else {
    k = orcmu * utau2;
    dk = orcmu * 2.0 * utau * dutau;
    eps = utau2 * utau / (vkcst * delta);
    deps = 3.0 * utau2 * dutau / (vkcst * delta);
  }

  V[5] = k;
  V[6] = eps;

  dV[5] = dk;
  dV[6] = deps;

}

//------------------------------------------------------------------------------
