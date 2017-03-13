#ifndef _POST_FCN_H_
#define _POST_FCN_H_

#include <NavierStokesTerm.h>
#include <SpalartAllmarasTerm.h>
#include <DESTerm.h>
#include <KEpsilonTerm.h>

class VarFcn;
class WallFcn;

// Included (MB)
class Communicator;

struct Vec3D;

//------------------------------------------------------------------------------

class PostFcn {

public:

// Included (MB)
  enum ScalarType {DENSITY = 0, MACH = 1, PRESSURE = 2, TEMPERATURE = 3, TOTPRESSURE = 4,
		   VORTICITY = 5, CSDLES = 6, CSDVMS = 7, SKIN_FRICTION = 8, NUT_TURB = 9, 
                   K_TURB = 10, EPS_TURB = 11, EDDY_VISCOSITY = 12, DELTA_PLUS = 13, 
                   PSENSOR = 14, MUT_OVER_MU = 15, PHILEVEL = 16,PHILEVEL2 = 17, DIFFPRESSURE = 18, 
                   SPEED = 19, HYDROSTATICPRESSURE = 20, HYDRODYNAMICPRESSURE = 21, 
                   WTMACH = 22, WTSPEED = 23, VELOCITY_NORM = 24, TEMPERATURE_NORMAL_DERIVATIVE = 25, 
                   SURFACE_HEAT_FLUX = 26, PRESSURECOEFFICIENT = 27, CONTROL_VOLUME = 28, FLUIDID = 29,
                   D2WALL = 30, SPATIAL_RES = 31, SSIZE = 32};

  enum VectorType {VELOCITY = 0, DISPLACEMENT = 1, FLIGHTDISPLACEMENT = 2, LOCALFLIGHTDISPLACEMENT = 3, VSIZE = 4};
  enum ScalarAvgType {DENSITYAVG = 0, MACHAVG = 1, PRESSUREAVG = 2, TEMPERATUREAVG = 3,
                      TOTPRESSUREAVG = 4, VORTICITYAVG = 5, CSDLESAVG = 6, CSDVMSAVG = 7, 
                      SKIN_FRICTIONAVG =8, AVSSIZE = 9};
  enum VectorAvgType {VELOCITYAVG = 0, DISPLACEMENTAVG = 1, AVVSIZE = 2};

// Included (MB)
  //list of ScalarDerivatives
  enum ScalarDerivativeType {DERIVATIVE_DENSITY = 0,
                             DERIVATIVE_MACH = 1,
			                 DERIVATIVE_PRESSURE = 2,
                             DERIVATIVE_TEMPERATURE = 3,
			                 DERIVATIVE_TOTPRESSURE = 4,
			                 DERIVATIVE_NUT_TURB = 5,
                             DERIVATIVE_EDDY_VISCOSITY = 6,
			                 DERIVATIVE_VELOCITY_SCALAR = 7,
			                 DERIVATIVE_SPATIAL_RES = 8,
			                 DSSIZE = 9 //just an auxiliary parameter that holds the size of the list
                            };

  //list of vector derivatives
  enum VectorDerivativeType {DERIVATIVE_VELOCITY_VECTOR = 0,
                             DERIVATIVE_DISPLACEMENT = 1,
			     DVSIZE = 2 //just an auxiliary parameter that holds the size of the list
                            };

protected:

  VarFcn *varFcn;

public:

  PostFcn(VarFcn *);
  virtual ~PostFcn() {}

  virtual double computeNodeScalarQuantity(ScalarType, double *, double *, int = 0,double* = NULL);
  virtual double computeFaceScalarQuantity(ScalarType, double [4][3], Vec3D&, double [3], 
					   double*, double* [3], double* [4]);
  virtual void computeForce(double [4][3], double *[3], Vec3D &, double [3], double *, double *[3],
		double *[4], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], int = 0, int fid = 0) = 0;

  virtual void computeForceEmbedded(int, double [4][3], 
				    double *[3], Vec3D &, double [3], 
				    double *, double *[3], double *[4], double, 
				    Vec3D &, Vec3D &, Vec3D &, Vec3D &, 
				    double[3][3], int = 0, 
				    int *fid = 0, bool = true) = 0;

  virtual Vec3D computeViscousForceCVBoundary(Vec3D& n,  double* Vi, double dudxj[3][3])
  {
    fprintf(stderr,"Calling a PostFcn Function for Viscous Forces. Doesn't make sense!\n");
    exit(-1);
  }
  virtual Vec3D computeViscousForce(double [4][3], Vec3D&, double [3], double*, double* [3], double* [4])
  {
    fprintf(stderr,"Calling a PostFcn Function for Viscous Forces. Doesn't make sense!\n");
    exit(-1);
  }
  virtual void computeForceTransmitted(double [4][3], double *[3], Vec3D &, double [3], double *, double *[3], 
	double *[4], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], int = 0, int fid = 0) = 0;

  virtual double computeHeatPower(double [4][3], Vec3D&, double [3],
				  double*, double* [3], double* [4]) = 0;

  virtual double computeHeatFluxRelatedValues(double [4][3], Vec3D& , double [3],
                                   double* , double* [3], double* [4], bool) =0;

  virtual double computeInterfaceWork(double [4][3], Vec3D&, double, double [3], double*, 
				      double* [3], double* [4], double) = 0;
  virtual bool doesFaceNeedGradientP1Function() { return false; }

  virtual double* getMeshVel()  { return varFcn->getMeshVel(); }
  
// Included (MB)
  virtual double computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType, double [3], double *, double *, double *, double *, double = 0);
  virtual void computeDerivativeOfForce(double [4][3], double [4][3], double *[3], double *[3], Vec3D &, Vec3D &,
                                        double [3], double *, double *, double *[3], double *[3],
                                        double *[4], double *[4], double [3], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], double[3][3], int = 0) = 0;
  virtual void computeDerivativeOfForce2(double [4][3], double [4][3], double *[3], double *[3], Vec3D &, Vec3D &,
                                         double [3], double *, double *, double *[3], double *[3],
                                         double *[4], double *[4], double [3], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], double[3][3], int = 0) = 0;
  virtual void computeDerivativeOperatorsOfForce(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double *Vface[3], double *Vtet[4], double *pin,
                                                 double dPdx[3][3], int hydro, double dFi0dn[3], double dFi1dn[3], double dFi2dn[3],
                                                 double dFi0ddPdx[3][3], double dFi1ddPdx[3][3], double dFi2ddPdx[3][3],
                                                 double dFi0dXface0[3][3], double dFi0dXface1[3][3], double dFi0dXface2[3][3],
                                                 double dFi1dXface0[3][3], double dFi1dXface1[3][3], double dFi1dXface2[3][3],
                                                 double dFi2dXface0[3][3], double dFi2dXface1[3][3], double dFi2dXface2[3][3],
                                                 double dFi0dS[3][3], double dFi1dS[3][3], double dFi2dS[3][3],
                                                 double dFi0dVface[3][5], double dFi1dVface[3][5], double dFi2dVface[3][5],
                                                 double dFvddp1dxj[3][4][3], double dFvdn[3][3], double dFvdV[3][4][5]) = 0;
  virtual void computeDerivativeOfForceTransmitted(double [4][3], double [4][3], double *[3], double *[3], Vec3D &, Vec3D &,
                                        double [3], double *, double *, double *[3], double *[3],
                                        double *[4], double *[4], double [3], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], double[3][3], int = 0) = 0;
  virtual void computeDerivativeOperatorsOfForceTransmitted(double [4][3],
                                                            double *Xface[3], Vec3D &n, double *Vface[3], double *Vtet[4], double *pin, double dPdx[3][3], int hydro,
                                                            double dFi0dn[3], double dFi0dS[3][3], double dFi0dVface[3][3][5],
                                                            double dFi0ddPdx[3][3][3], double dFi0dXface[3][3][3],
                                                            double dFi1dn[3], double dFi1dS[3][3], double dFi1dVface[3][3][5],
                                                            double dFi1ddPdx[3][3][3], double dFi1dXface[3][3][3],
                                                            double dFi2dn[3], double dFi2dS[3][3], double dFi2dVface[3][3][5],
															double dFi2ddPdx[3][3][3], double dFi2dXface[3][3][3],
															double [3][4][3], double [3][3], double [3][4][5]) = 0;
  virtual void rstVar(IoData &, Communicator*) = 0;
  virtual double computeDerivativeOfHeatPower(double [4][3], double [4][3], Vec3D&, Vec3D&, double [3], double*, double*, double* [3], double* [3], double* [4], double* [4], double [3]) = 0;

};

//------------------------------------------------------------------------------

class PostFcnEuler : public PostFcn {

protected:

  static const double third;

  double mach;
  double pinfty;
  double depth;
  double gravity;
  double alpha;
  double beta;
  double nGravity[3];

// Included (MB)
  bool dimFlag;
  double dpinfty;
  double dPin;

public:

  PostFcnEuler(IoData &, VarFcn *);
  virtual ~PostFcnEuler() {}

  virtual double computeNodeScalarQuantity(ScalarType, double *, double *, int = 0,double* = NULL);
  virtual void computeForce(double [4][3], double *[3], Vec3D &, double [3], double *, double *[3],
                double *[4], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], int = 0, int fid = 0);

  virtual void computeForceEmbedded(int, double [4][3], double *[3], Vec3D &, 
				    double [3], double *, double *[3],  double *[4], double, 
				    Vec3D &, Vec3D &, Vec3D &, Vec3D &, 
				    double[3][3], int = 0, int* fid = 0, bool = true);

  virtual Vec3D computeViscousForceCVBoundary(Vec3D& n,  double* Vi, double dudxj[3][3])
  {
    fprintf(stderr,"Calling a PostFcnEuler Function for Viscous Forces. Doesn't make sense!\n");
    exit(-1);
  }
  virtual Vec3D computeViscousForce(double [4][3], Vec3D&, double [3], double*, double* [3], double* [4])
  {
    fprintf(stderr,"Calling a PostFcnEuler Function for Viscous Forces. Doesn't make sense!\n");
    exit(-1);
  }
  virtual void computeForceTransmitted(double [4][3], double *[3], Vec3D &, double [3], double *, double *[3],
				       double *[4], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], int = 0, int fid = 0);
  virtual double computeHeatPower(double [4][3], Vec3D&, double [3],
				  double*, double* [3], double* [4]);
  virtual double computeHeatFluxRelatedValues(double [4][3], Vec3D& , double [3],
                                   double* , double* [3], double* [4], bool);
  virtual double computeInterfaceWork(double [4][3], Vec3D&, double, double [3], double*, 
				      double* [3], double* [4], double);


// Included (MB)
  virtual double computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType, double [3], double *, double *, double *, double *, double = 0);
  virtual void computeDerivativeOfForce(double [4][3], double [4][3], double *[3], double *[3], Vec3D &, Vec3D &,
                                        double [3], double *, double *, double *[3], double *[3],
                                        double *[4], double *[4], double [3], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], double[3][3], int = 0);
  virtual void computeDerivativeOperatorsOfForce(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double *Vface[3], double *Vtet[4], double *pin,
                                                 double dPdx[3][3], int hydro, double dFi0dn[3], double dFi1dn[3], double dFi2dn[3],
                                                 double dFi0ddPdx[3][3], double dFi1ddPdx[3][3], double dFi2ddPdx[3][3],
                                                 double dFi0dXface0[3][3], double dFi0dXface1[3][3], double dFi0dXface2[3][3],
                                                 double dFi1dXface0[3][3], double dFi1dXface1[3][3], double dFi1dXface2[3][3],
                                                 double dFi2dXface0[3][3], double dFi2dXface1[3][3], double dFi2dXface2[3][3],
                                                 double dFi0dS[3][3], double dFi1dS[3][3], double dFi2dS[3][3],
                                                 double dFi0dVface[3][5], double dFi1dVface[3][5], double dFi2dVface[3][5],
                                                 double dFvddp1dxj[3][4][3], double dFvdn[3][3], double dFvdV[3][4][5]);
  virtual void computeDerivativeOfForce2(double [4][3], double [4][3], double *[3], double *[3], Vec3D &, Vec3D &,
                                         double [3], double *, double *, double *[3], double *[3],
                                         double *[4], double *[4], double [3], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], double[3][3], int = 0);
  virtual void computeDerivativeOfForceTransmitted(double [4][3], double [4][3], double *[3], double *[3], Vec3D &, Vec3D &,
                                        double [3], double *, double *, double *[3], double *[3],
                                        double *[4], double *[4], double [3], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], double[3][3], int = 0);
  virtual void computeDerivativeOperatorsOfForceTransmitted(double [4][3],
                                                            double *Xface[3], Vec3D &n, double *Vface[3], double *Vtet[4], double *pin, double dPdx[3][3], int hydro,
                                                            double dFi0dn[3], double dFi0dS[3][3], double dFi0dVface[3][3][5],
                                                            double dFi0ddPdx[3][3][3], double dFi0dXface[3][3][3],
                                                            double dFi1dn[3], double dFi1dS[3][3], double dFi1dVface[3][3][5],
                                                            double dFi1ddPdx[3][3][3], double dFi1dXface[3][3][3],
                                                            double dFi2dn[3], double dFi2dS[3][3], double dFi2dVface[3][3][5],
															double dFi2ddPdx[3][3][3], double dFi2dXface[3][3][3], double [3][4][3], double [3][3], double [3][4][5]);
  void rstVar(IoData &, Communicator*);
  virtual double computeDerivativeOfHeatPower(double [4][3], double [4][3], Vec3D&, Vec3D&, double [3], double*, double*, double* [3], double* [3], double* [4], double* [4], double [3]);

};

//------------------------------------------------------------------------------

class PostFcnNS : public PostFcnEuler, public NavierStokesTerm {

protected:

  WallFcn* wallFcn;

private:

// Included (MB)
  Vec3D computeDerivativeOfViscousForce(double [4][3], double [4][3], Vec3D&, Vec3D&, double [3], double*, double*, double* [3], double* [3], double* [4], double* [4], double [3]);
  void computeDerivativeOperatorsOfViscousForce(double dp1dxj[4][3], Vec3D& n, double* Vtet[4], double dFvddp1dxj[3][4][3], double dFvdn[3][3], double [3][4][5]); //YC

public:

  PostFcnNS(IoData &, VarFcn *);
  virtual ~PostFcnNS();

  double computeFaceScalarQuantity(ScalarType, double [4][3], Vec3D&, double [3], 
				   double*, double* [3], double* [4]);
  virtual void computeForce(double [4][3], double *[3], Vec3D &, double [3], double *, double *[3],
                double *[4], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], int = 0, int fid = 0);

  virtual void computeForceEmbedded(int, double [4][3], 
				    double *[3], Vec3D &, double [3], 
				    double *, double *[3], double *[4], double, 
				    Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], 
				    int = 0, int* fid = 0, bool = true);

  virtual Vec3D computeViscousForceCVBoundary(Vec3D& n,  double* Vi, double dudxj[3][3]);
  virtual Vec3D computeViscousForce(double [4][3], Vec3D&, double [3], double*, double* [3], double* [4]);

  virtual void computeForceTransmitted(double [4][3], double *[3], Vec3D &, double [3], double *, double *[3],
				       double *[4], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], int = 0, int fid = 0);

  double computeHeatPower(double [4][3], Vec3D&, double [3],
			  double*, double* [3], double* [4]);
  virtual double computeHeatFluxRelatedValues(double [4][3], Vec3D& , double [3],
                                   double* , double* [3], double* [4], bool);
  double computeInterfaceWork(double [4][3], Vec3D&, double, double [3], double*, 
			      double* [3], double* [4], double);
  bool doesFaceNeedGradientP1Function() { return ((wallFcn) ? false : true); }
  
// Included (MB)
  double computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType, double [3], double *, double *, double *, double *, double = 0);

  void computeDerivativeOfForce(double [4][3], double [4][3], double *[3], double *[3], Vec3D &, Vec3D &,
                                        double [3], double *, double *, double *[3], double *[3],
                                        double *[4], double *[4], double [3], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], double[3][3], int = 0);
  void computeDerivativeOfForceTransmitted(double [4][3], double [4][3], double *[3], double *[3], Vec3D &, Vec3D &,
                                        double [3], double *, double *, double *[3], double *[3],
                                        double *[4], double *[4], double [3], double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &, double[3][3], double[3][3], int = 0);

  virtual void computeDerivativeOperatorsOfForceTransmitted(double [4][3], double *Xface[3], Vec3D &n, double *Vface[3], double *Vtet[4], double *pin, double dPdx[3][3], int hydro,
                                                    double dFi0dn[3], double dFi0dS[3][3], double dFi0dVface[3][3][5],
                                                    double dFi0ddPdx[3][3][3], double dFi0dXface[3][3][3],
                                                    double dFi1dn[3], double dFi1dS[3][3], double dFi1dVface[3][3][5],
                                                    double dFi1ddPdx[3][3][3], double dFi1dXface[3][3][3],
                                                    double dFi2dn[3], double dFi2dS[3][3], double dFi2dVface[3][3][5],
                                                    double dFi2ddPdx[3][3][3], double dFi2dXface[3][3][3], double [3][4][3], double [3][3], double [3][4][5]);
  virtual void computeDerivativeOperatorsOfForce(double dp1dxj[4][3], double *Xface[3], Vec3D &n, double *Vface[3], double *Vtet[4], double *pin,
                                                 double dPdx[3][3], int hydro, double dFi0dn[3], double dFi1dn[3], double dFi2dn[3],
                                                 double dFi0ddPdx[3][3], double dFi1ddPdx[3][3], double dFi2ddPdx[3][3],
                                                 double dFi0dXface0[3][3], double dFi0dXface1[3][3], double dFi0dXface2[3][3],
                                                 double dFi1dXface0[3][3], double dFi1dXface1[3][3], double dFi1dXface2[3][3],
                                                 double dFi2dXface0[3][3], double dFi2dXface1[3][3], double dFi2dXface2[3][3],
                                                 double dFi0dS[3][3], double dFi1dS[3][3], double dFi2dS[3][3],
                                                 double dFi0dVface[3][5], double dFi1dVface[3][5], double dFi2dVface[3][5],
                                                 double dFvddp1dxj[3][4][3], double dFvdn[3][3], double dFvdV[3][4][5]);

  void rstVar(IoData &, Communicator*);
  double computeDerivativeOfHeatPower(double [4][3], double [4][3], Vec3D&, Vec3D&, double [3], double*, double*, double* [3], double* [3], double* [4], double* [4], double [3]);

};

//------------------------------------------------------------------------------

class PostFcnSA : public PostFcnNS, public SATerm {

public:

  PostFcnSA(IoData &, VarFcn *);
  virtual ~PostFcnSA();

  double computeNodeScalarQuantity(ScalarType, double *, double *, int = 0,double* = NULL);
  
// Included (MB)
  void rstVar(IoData &, Communicator*);
  double computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType, double [3], double *, double *, double *, double *, double = 0);

};

//------------------------------------------------------------------------------

class PostFcnDES : public PostFcnNS, public DESTerm {

public:

  PostFcnDES(IoData &, VarFcn *);
  virtual ~PostFcnDES() {}

  double computeNodeScalarQuantity(ScalarType, double *, double *, int = 0,double* = NULL);
  
// Included (MB)
  void rstVar(IoData &, Communicator*);
  double computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType, double [3], double *, double *, double *, double *, double = 0);

};

//------------------------------------------------------------------------------

class PostFcnKE : public PostFcnNS, public KEpsilonTerm {

public:

  PostFcnKE(IoData &, VarFcn *);
  virtual ~PostFcnKE() {}

  double computeNodeScalarQuantity(ScalarType, double *, double *, int = 0,double* = NULL);
  
// Included (MB)
  void rstVar(IoData &, Communicator*);
  double computeDerivativeOfNodeScalarQuantity(ScalarDerivativeType, double [3], double *, double *, double *, double *, double = 0);

};

//------------------------------------------------------------------------------

#endif
