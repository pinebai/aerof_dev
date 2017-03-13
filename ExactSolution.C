// Exact Solution.cpp

#include <IoData.h>

#include "ExactSolution.h"

#include <cmath>

#include <complex>

#ifdef AEROACOUSTIC

#include "gsl/gsl_sf.h"

#endif // AEROACOUSTIC


void ExactSolution::
AcousticBeam(IoData& iod,double x, double y, double z,
	     double t, double* V) {

 /* 
  double alpha = 3.810627283;//6.28267340353564;
  double omega0 = 2697.56348148;//26.975634814780758;
  double omegatilde = 0.60774519311;//0.9519666394290971;
  double H = 1.0;
  double omega = omega0*omegatilde;
  double what = 1.0e-6;
  double k = 2.0*3.14159265358979323846;
  double rhof = 1.3;
  */

  double alpha = 3.8116859813491843;//6.28267340353564;
  double omega0 = 2696.889343398955;//26.975634814780758;
  double omegatilde = 0.6077988242743707;//0.9519666394290971;
  double H = 1.0;
  double omega = omega0*omegatilde;
  double what = 1.0e-6;
  double k = 2.0*3.14159265358979323846;
  double rhof = 1.3;
 
  for (int i = 0; i < 5; ++i)
    V[i] = 0.0;
  
  double u = omega*k*what* 
    (cosh(alpha*y) /(alpha*sinh(alpha*H)))*(cos(k*x-omega*t*iod.ref.rv.time)+
					       cos(-k*x-omega*t*iod.ref.rv.time)) / iod.ref.rv.velocity;
  double v = omega*what*(sinh(alpha*y)/sinh(alpha*H))*
    (sin(k*x-omega*t*iod.ref.rv.time)-
     sin(-k*x-omega*t*iod.ref.rv.time)) / iod.ref.rv.velocity;
  double p = omega*omega*rhof*cosh(alpha*y)*what/(alpha*sinh(alpha*H))*
    (cos(k*x-omega*t*iod.ref.rv.time)-
     cos(-k*x-omega*t*iod.ref.rv.time)) / iod.ref.rv.pressure + iod.bc.inlet.pressure;

  V[0] = pow(p/(1.0e5/iod.ref.rv.pressure),(1.0/1.4))*(1.3 / iod.ref.rv.density);
  V[1] = u;
  V[2] = v;
  V[3] = 0.0;
  V[4] = p;

}

void ExactSolution::
AcousticBeamStructure(IoData& iod,double x, double y, double z,
  	              double t, double& uy, double& vy) {

  
  double alpha = 3.8116859813491843;//6.28267340353564;
  double omega0 = 2696.889343398955;//26.975634814780758;
  double omegatilde = 0.6077988242743707;//0.9519666394290971;
  double H = 1.0;
  double omega = omega0*omegatilde;
  double what = 1.0e-6;
  double k = 2.0*3.14159265358979323846;
  double rhof = 1.3;

  uy = what* (cos(k*x-omega*t*iod.ref.rv.time)-cos(-k*x-omega*t*iod.ref.rv.time));
  vy = omega*what* (sin(k*x-omega*t*iod.ref.rv.time)-
	             sin(-k*x-omega*t*iod.ref.rv.time)) / iod.ref.rv.velocity;
}

void ExactSolution::
AcousticViscousBeam(IoData& iod,double x, double y, double z,
	            double t, double* V) {

  t *= iod.ref.rv.time;

  double k = 2.0*3.14159265358979323846;
  double rhof = 1.3;
  double rhos = 10,h=0.02;
  double Ib = pow(h,3)/12.0;
  std::complex<double> I(0,1);
  double H = 1.0;
  std::complex<double> omega(59.0185,-0.824694);

  double c = sqrt(1.4*1e5/1.3);

  double nu = 0.01/rhof;

  std::complex<double> alpha = sqrt(k*k- omega*omega/(c*c-I*omega*nu));

  std::complex<double> albar = sqrt(k*k-I*omega/nu);

  std::complex<double> phi = 1.0/(nu*(alpha*alpha-k*k)+I*omega);

  std::complex<double>  q = I*omega/(rhof*c*c)-alpha*alpha*phi/rhof;

  std::complex<double> what(1e-6,0);
  
  std::complex<double> g2 = 1.4e6*pow(k,4)*Ib- rhos*h*pow(omega,2);

  std::complex<double> A = g2*what/(cosh(alpha*H)+sinh(alpha*H)*q*rhof/(alpha*albar*phi));
  
  std::complex<double> B = A*q*rhof/(alpha*albar*phi);
 
  std::complex<double> M = B*alpha*phi/rhof;
  std::complex<double> N = A*alpha*phi/rhof;

  std::complex<double> P = -I*omega*what-M*cosh(alpha*H)-N*sinh(alpha*H);

  std::complex<double> p = (A*cosh(alpha*y)+B*sinh(alpha*y))*(exp(I*k*x-I*omega*t)-exp(-I*k*x-I*omega*t)) / iod.ref.rv.pressure + iod.bc.inlet.pressure;

  std::complex<double> v = (M*cosh(alpha*y) + N*sinh(alpha*y) + P*exp(albar*(y-H))-M*exp(-albar*y))*(exp(I*k*x-I*omega*t)-exp(-I*k*x-I*omega*t)) / iod.ref.rv.velocity;
  std::complex<double> dv = (alpha*M*sinh(alpha*y) + alpha*N*cosh(alpha*y) + albar*P*exp(albar*(y-H))+M*albar*exp(-albar*y));

  std::complex<double> u = ((A*cosh(alpha*y)+B*sinh(alpha*y))*I*omega/(rhof*c*c)-dv)/(I*k)*(exp(I*k*x-I*omega*t)+exp(-I*k*x-I*omega*t)) / iod.ref.rv.velocity;

//  if (y > H)
//   std::cout <<  " " << x << " " << y << " " << z << " " << t << " " <<  (M*cosh(alpha*H) + N*sinh(alpha*H) + P*exp(albar*(H-H))-M*exp(-albar*H))*(exp(I*k*x-I*omega*t)-exp(-I*k*x-I*omega*t)) << " " << -I*omega*what*(exp(I*k*x-I*omega*t)-exp(-I*k*x-I*omega*t)) << std::endl;

  V[0] = pow(real(p)/(1.0e5/iod.ref.rv.pressure),(1.0/1.4))*(1.3 / iod.ref.rv.density);
  V[1] = real(u);
  V[2] = real(v);
  V[3] = 0.0;
  V[4] = real(p);


  
}

void ExactSolution::
AcousticViscousBeamStructure(IoData& iod,double x, double y, double z,
  	                     double t, double& uy, double& vy) {

  t *= iod.ref.rv.time;
  
  double k = 2.0*3.14159265358979323846;
  std::complex<double> omega(59.0185,-0.824694);

  std::complex<double> what(1e-6,0);

  std::complex<double> I(0,1);
 
  std::complex<double> w = what*(exp(I*k*x-I*omega*t)-exp(-I*k*x-I*omega*t));
  std::complex<double> v = -I*omega*what*(exp(I*k*x-I*omega*t)-exp(-I*k*x-I*omega*t)) / iod.ref.rv.velocity;


  uy = real(w);
  vy = real(v); 

}


#ifdef AEROACOUSTIC
double j0prime(double r) {

  return -gsl_sf_bessel_J1(r);
}

double y0prime(double r) {

  return -gsl_sf_bessel_Y1(r);
}

#endif

void ExactSolution::
CylindricalBubble(IoData& iod,double x, double y, double z,
		  double t, double* V, double* phi, int& fid) {

#ifdef AEROACOUSTIC
  double omega = 733.2541686104948;
  double Bn = 10.0;//1.0e-1;
  double r = sqrt(x*x+y*y);
  double rhoo = 10.0;
  double co = sqrt(4.4*1e5/rhoo);
  double rhoi = 1.0;
  double ci = sqrt(1.4*1e5/rhoi);
  double a = 0.5;
  double R = 1.0;
  double pinf = 1e5;
  
  t *= iod.ref.rv.time;

  double theta = atan2(y,x);

  double Cn = -Bn*j0prime(omega*R / co) / y0prime(omega*R / co);

  double An = (Bn*gsl_sf_bessel_J0(omega*a / co) + Cn*gsl_sf_bessel_Y0(omega*a / co))/gsl_sf_bessel_J0(omega*a / ci);

  double ur, utheta = 0;
  double p;


  if (r < a) {
    
    ur = -An/(rhoi*ci)*j0prime(omega*r/ci)*cos(omega*t);
    p = -An*gsl_sf_bessel_J0(omega*r / ci) *sin(omega*t);
  } else {
    ur = -1.0/(rhoo*co)*(Bn*j0prime(omega*r/co)+Cn*y0prime(omega*r/co))*cos(omega*t);
    p = -(Bn*gsl_sf_bessel_J0(omega*r / co) + Cn*gsl_sf_bessel_Y0(omega*r / co))*sin(omega*t);
  }

  double ux = ur*cos(theta), uy = ur*sin(theta);

  p += pinf;

  p /= iod.ref.rv.pressure;

  if (r < a)
    V[0] = pow(p/(1.0e5/iod.ref.rv.pressure),(1.0/1.4))*(1.0 / iod.ref.rv.density);
  else
    V[0] = pow(p/(1.0e5/iod.ref.rv.pressure),(1.0/4.4))*(10.0 / iod.ref.rv.density);

  V[1] = ux/ iod.ref.rv.velocity;
  V[2] = uy/ iod.ref.rv.velocity;
  
  V[3] = 0.0;
  
  V[4] = p;  

  fid = (r <= a);

  *phi = a-r;
  
#else
  std::cout << "Error: Cylindrical bubble exact solution requires the code to be compiled with the AEROACOUSTIC option" << std::endl;
  exit(-1);
#endif
}

void ExactSolution::
AcousticTwoFluid(IoData& iod,double x, double y, double z,
		 double t, double* V, double* phi, int& fid) {


  double alpha1 = 3.017303450887337;
  double alpha2 = 4.920327228498049;
  double H = 1.0;
  double omega = 1462.0641776748855;
  double k = 2.0*3.14159265358979323846;
  double rhof1 = 10,rhof2 = 1.0;
  double c1 = sqrt(4.4*1e5/rhof1);
  double c2 = sqrt(1.4*1e5/rhof2);
  double a = 0.5;
 
  for (int i = 0; i < 5; ++i)
    V[i] = 0.0;

  float f = cos(k*x - omega*t*iod.ref.rv.time);
  float fs = sin(k*x - omega*t*iod.ref.rv.time);
  //float A2overA1 = -rhof2/rhof1 *alpha1/alpha2* sinh(alpha1*a)/sinh(alpha2*(H-a));
  float A2overA1 = -rhof2/rhof1 *alpha1/alpha2* sin(alpha1*a)/sinh(alpha2*(H-a));

  float A1 = 0.1;
  float A2 = A1*A2overA1;

  std::cout.precision(15);

  //std::cout << (omega*omega)/(c1*c1) << " " << k*k+alpha1*alpha1 << std::endl;
  //std::cout << (omega*omega)/(c2*c2) << " " << k*k-alpha2*alpha2 << std::endl;

  if (y < a) {
    /*
    double u = k*f/(rhof1*omega)*A1*cosh(alpha1*y) / iod.ref.rv.velocity;
    double v = A1*fs*alpha1*sinh(alpha1*y)/(rhof1*omega) / iod.ref.rv.velocity;
    double p = A1*cosh(alpha1*y)*f / iod.ref.rv.pressure + iod.bc.inlet.pressure;
    */

    double u = -k*f/(rhof1*omega)*A1*cos(alpha1*y) / iod.ref.rv.velocity;
    double v = A1*fs*alpha1*sin(alpha1*y)/(rhof1*omega) / iod.ref.rv.velocity;
    double p = -A1*cos(alpha1*y)*f / iod.ref.rv.pressure + iod.bc.inlet.pressure;
    
    V[0] = pow(p/(1.0e5/iod.ref.rv.pressure),(1.0/4.4))*(10.0 / iod.ref.rv.density);
    V[1] = u;
    V[2] = v;
    V[3] = 0.0; 
    V[4] = p;
  } else {
    
    double u = k*f/(rhof2*omega)*A2*cosh(alpha2*(H-y)) / iod.ref.rv.velocity;
    double v = -A2*fs*alpha2*sinh(alpha2*(H-y))/(rhof2*omega) / iod.ref.rv.velocity;
    double p = A2*cosh(alpha2*(H-y))*f / iod.ref.rv.pressure + iod.bc.inlet.pressure;
    /*
    double u = k*f/(rhof2*omega)*A2*cosh(alpha2*(H-y)) / iod.ref.rv.velocity;
    double v = -A2*fs*alpha2*sinh(alpha2*(H-y))/(rhof2*omega) / iod.ref.rv.velocity;
    double p = A2*cosh(alpha2*(H-y))*f / iod.ref.rv.pressure + iod.bc.inlet.pressure;
    */
    V[0] = pow(p/(1.0e5/iod.ref.rv.pressure),(1.0/1.4))*(1.0 / iod.ref.rv.density);
    V[1] = u;
    V[2] = v;
    V[3] = 0.0;
    V[4] = p;
  }

  *phi = y-a;

  fid = (y > a);

}
