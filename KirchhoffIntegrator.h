#ifndef _KIRCHHOFF_INTEGRATOR_H
#define _KIRCHHOFF_INTEGRATOR_H


//------------------------------------------------------------------------------
#include <complex>
#include <vector>


class Domain;
class IoData;
class Timer;


template<class Scalar, int dim> class DistSVec;


class KirchhoffIntegrator {
  
private:
  
  // Do not define these functions
  KirchhoffIntegrator(const KirchhoffIntegrator &ref);
  KirchhoffIntegrator & operator=(const KirchhoffIntegrator &ref);

  
  ///////////////////////
  // PRIVATE Utility functions
  ///////////////////////
  
    
  // This function returns quadrature points on a face.
  //
  // \note UH (08/2012) It is implemented only for triangles
  // with a {3, 6, 12}-points quadrature rule.
  void getQuadrature
  (
   std::vector<double> &_xigauss,
   std::vector<double> &w,
   int rule = 3
   ) const;
  
  
  // This function converts cartesian coordinates to spherical coordinates.
  // theta is between 0 and pi, phi between 0 and 2*pi
  //
  // \note The cut is for {x > 0, y = 0}
  void cart2sph_x
  (
   double xo, double yo, double zo,
   double &ro, double &to, double &po
   ) const;
  
  
  // This function converts cartesian coordinates to spherical coordinates.
  // theta is between 0 and pi, phi between -pi/2 and 3*pi/2
  //
  // \note The cut is for {x = 0, y < 0}
  void cart2sph_y
  (
   double xo, double yo, double zo,
   double &ro, double &to, double &po
   ) const;
  
  
  // Function to evaluate the cylindrical Hankel H_n at a point x.
  std::complex<double> hankel
  (
   const int n, const double x
   ) const;
  
  
  // Function to evaluate the derivative of the Hankel H_n at a point x.
  std::complex<double> hankel_prime
  (
   const int n, const double x
   ) const;
  
  
  // Function to evaluate the spherical bessel h_n at a point x.
  std::complex<double> besselh
  (
   const int n, const double x
   ) const;
  
  
  // Function to evaluate the derivative of a spherical bessel h_n.
  std::complex<double> besselh_prime
  (
   const int n, const double x
   ) const;
  
  
  // Function to return the elemental area at a point
  // The integration point (xp, yp, zp) is placed on the surface.
  //
  // \note UH (09/2012)
  // This function works only with spherical and cylindrical surface.
  // The cylinder is assumed to be aligned with the y-axis.
  double dSigma
  (
   const double *x, const double *y, const double *z,
   double &xp, double &yp, double &zp
   )
  const;
  
  
  // This function computes the normalized spherical harmonics
  // at a point (theta, phi).
  // For a given value n, all values between -n and n are computed.
  // They are stored in the output array Yn.
  // Note: Yn may be resized.
  void sphericalHarmonic
  (
   const int n,
   const double theta,
   const double phi,
   std::vector< std::complex<double> > &Yn
   ) const;
  
  
protected:
  
  IoData &d_iod;
  Domain *d_domain_p;
  DistSVec<double,3> *d_X_p;
  
  std::vector<int> d_globalToNodeID;
  
  Timer *d_timer_p;
  
  enum TypeGamma {SPHERE = 0, CYLINDER = 1} d_SurfType;
  double d_R;
  double d_cyl_axis[2];
  
  std::vector<double> d_omega_t;
  
  
  ///////////////////////
  // Utility functions
  ///////////////////////
  
  
  // This function computes the values and the FFP for a spherical surface.
  void Sphere
  (
   std::vector< std::complex<double> > &nodal,
   const int *vecLen
   );
  
  
  // This function computes the values and the FFP for a cylindrical surface.
  void Cylinder_TensorGrid
  (
   std::vector< std::complex<double> > &nodal,
   const int *vecLen
   );
  
  
  // Function to check whether the grid is orthogonal.
  // Returns a flag whether the grid is orthogonal or not
  // and the numbers of terms to use in the series expansion.
  //
  // \note UH (09/2012)
  // This function sets also the radius of the cylinder (d_R)
  // and the axial interval (d_cyl_axis[0] and d_cyl_axis[1]).
  //
  void CylinderGrid
  (
   bool &tensorGrid,
   int &nKappa,
   int &nTheta
   );
  
  
  // This function extract the probe locations.
  //
  void getObservationPoints
  (
   std::vector<double> &observations
   ) const;
  
  
  // This function writes the pressure values in an ASCII file.
  void writePValues
  (
   std::vector< std::complex<double> > &pnoise
   ) const;
  
  
  // This function writes the FFP values in an ASCII file.
  void writeFFPValues
  (
   std::vector< std::complex<double> > &ffp
   ) const;
  
  
  // This function converts nodal values on a spherical surface
  // to the coefficients in a spherical harmonic series expansion.
  //
  // \note UH (09/2012)
  // This function sets the value for d_R (sphere radius).
  void convertToSHSeries
  (
   std::vector< std::complex<double> > &nodal,
   const int *vecLen,
   std::vector< std::complex<double> > &series,
   int &nmax
   );
  
  
  // This function computes the coefficients in a spherical harmonic series
  // expansion for the normal derivative of the pressure.
  // The pressure is also specified by the coefficients in a spherical
  // harmonic series.
  void getdpdnSHseries
  (
   const std::vector< std::complex<double> > &coeff,
   std::vector< std::complex<double> > &coeff_dudn,
   const int nMax
   ) const;
  
  
  // This function computes the Kirchhoff integral on a spherical
  // surface.
  //
  // \note UH (09/2012)
  // This function is not called but it is kept for comparison.
  void integrateOnSphere
  (
   const std::complex<double> *pvalues,
   const int *vecLen,
   const std::complex<double> *dpdn,
   int nmax,
   std::complex<double> *p_series
   );
  
  
  // This function evaluates a series in spherical harmonics
  // at probes locations.
  //
  void EvaluateSHSatProbes
  (
   std::vector< std::complex<double> > &p_coeff,
   int nmax
   );

  // This function computes the far-field pattern
  // when the surface is a sphere.
  void ffpDataOnSphere
  (
   const std::complex<double> *pvalues,
   const int *vecLen,
   const std::complex<double> *dpdn,
   const int nmax
   ) const;
  
  
  // This function evaluates a series in spherical harmonics
  // at a point (theta, phi) with the input coefficients.
  std::complex<double> evaluateSHS
  (
   const std::complex<double> *coeff,
   int nmax,
   double theta,
   double phi,
   double r = -1.0,
   double kappa = 0.0
   ) const;
  
  
public:
  
  
  //! Constructor
  ///
  ///
  KirchhoffIntegrator
  (
   IoData &iod,
   Domain *domain
   );
  
  
  //! Destructor
  ///
  ///
  ~KirchhoffIntegrator();
  
  
  //! Function to compute the pressure at different points
  //! and to compute the far-field pattern at specified directions.
  ///
  /// This function evaluates the Kirchhoff integral from snapshots
  /// of pressure nodal values.
  ///
  void Compute();
  
  
};


//------------------------------------------------------------------------------


#endif

