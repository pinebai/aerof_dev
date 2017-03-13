#ifndef _FEM_EQUATION_TERM_DESC_H_
#define _FEM_EQUATION_TERM_DESC_H_

#include <FemEquationTerm.h>
#include <NavierStokesTerm.h>
#include <SpalartAllmarasTerm.h>
#include <KEpsilonTerm.h>
#include <DESTerm.h>
#include <Vector3D.h>

#include <cstdlib>
#include <cstdio>

class IoData;

//------------------------------------------------------------------------------

class FemEquationTermNS : public FemEquationTerm, public NavierStokesTerm {
  
public:
  
  //map<int, PorousMedia *> &volInfo;
  double velocity, density, length;
  
private:
  void computeTransportCoefficients(const double T, double &mu, double &lambda, double &kappa);

public:
  FemEquationTermNS(IoData &, VarFcn *);
  ~FemEquationTermNS() {}
  
  double computeViscousTimeStep(double *, double *);
  
  bool computeVolumeTerm(double [4][3], double [4], double *[4],
                         double *, double *, double *, double, 
                         SVec<double,3> &, int [4], int);
  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], 
                                 double *, double *, double *, double,
                                 SVec<double,3> &, int [4], int);
  void computeSurfaceTerm(int, Vec3D &, double [3], 
                          double *, double *[3], double *);

  double computeNormDerivWallFcn(double rho, double T, double Du1, double DT1, 
										 double d2w, double &dudn, double &dTdn);

  void computeJacobianSurfaceTerm(int, Vec3D &, double [3], double *, double *[3], double *);
  void computeSurfaceTerm(double [4][3], int, Vec3D &, double [4], 
                          double *, double *[4], double *);

  void computeDerivativeOperatorsOfVolumeTerm // YC
  (double [4][3], double *[4], double (*)[5][4][3], double (*)[5][4][5], double (*)[5]);

  void computeJacobianSurfaceTerm(double [4][3], int, Vec3D &, double [4], 
                                  double *, double *[4], double *);
  
  // Included (MB)
  bool computeDerivativeOfVolumeTerm
  (
   double [4][3], double [4][3], double [4], double *[4], double *[4],
   double, double *, double *, double *, double, SVec<double,3> &, int [4], int
   );
  
  
  void computeDerivativeOfSurfaceTerm
  (
   int, Vec3D &, Vec3D &, double [3], double *, double *, double *[3], 
   double *[3], double, double *
   );
  
  
  void computeDerivativeOfSurfaceTerm
  (
   double [4][3], double [4][3], int, Vec3D &, Vec3D &, double [4],
   double *, double *, double *[4], double *[4], double, double *
   );
 
  void rstVar
  (IoData &ioData, Communicator *com) 
  { 
    rstVarNS(ioData, com); 
    if (wallFcn) wallFcn->rstVar(ioData, com);
  };
  
  
  void computeBCsJacobianWallValues
  (
   int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]
   ) 
  {};
  
  
  double computeDerivativeOfViscousTimeStep
  (
   double *, double *, double *, double *, double
   );
  
  void computeTransportCoefficientsPublic(const double T, double &mu, 
                                          double &lambda, double &kappa);

  bool withWallFcn() 
  { 
	  if(wallFcn) return true;
	  else        return false;
  };

};

//------------------------------------------------------------------------------

class FemEquationTermSA : public FemEquationTerm, public NavierStokesTerm, 
public SATerm {
  
public:
  
  double x0,y0,z0,x1,y1,z1;
  bool trip, usefv3;
  //map<int, PorousMedia *> &volInfo;
  double velocity, density, length;
  ConstantPrandtlThermalCondFcn turbThermalCondFcn;
  
private:
  void computeLaminarTransportCoefficients(const double T, 
		  double &mu, double &lambda, double &kappa);

public:
  FemEquationTermSA(IoData &, VarFcn *);
  ~FemEquationTermSA() {}
  
  double computeViscousTimeStep(double *, double *);
  
  void computeTurbulentTransportCoefficients(double *V[], int nodeNum[], SVec<double,3> &X,
		  const double mul, const double lambdal, double kappal,
		  double &mutilde, double &mut, double &lambdat, double &kappat);
  
  bool computeVolumeTerm
  (
   double dp1dxj[4][3], double d2w[4], double *V[4],
   double *r, double *S, double *PR, 
   double tetVol, SVec<double,3> &X, 
   int nodeNum[4], int material_id
   );
  
  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], double *, double *, double *, double,
                                 SVec<double,3> &, int [4], int);
  void computeDerivativeOperatorsOfVolumeTerm // YC
  (
    double [4][3], double *[4], double (*)[5][4][3], double (*)[5][4][5], double (*)[5]) {
    fprintf(stderr, "*** Error: FemEquationTermSA::computeDerivativeOperatorsOfVolumeTerm is not implemented yet\n");
    exit(-1);
  }
  void computeSurfaceTerm(int, Vec3D &, double [3], 
                          double *, double *[3], double *);
  void computeJacobianSurfaceTerm(int, Vec3D &, double [3], double *, double *[3], double *);
  void computeSurfaceTerm(double [4][3], int, Vec3D &, double [4], 
                          double *, double *[4], double *);

  double computeNormDerivWallFcn(double rho, double T, double Du1, double DT1,
										 double d2w, double &dudn, double &dTdn);

  void computeJacobianSurfaceTerm(double [4][3], int, Vec3D &, double [4], 
                                  double *, double *[4], double *);
  bool doesSourceTermExist() { return true; }
  
  // Included (MB)
  bool computeDerivativeOfVolumeTerm
  (
   double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4],
   double *V[4], double *dV[4],
   double dMach, double *dr, double *dS, double *dPR, double dtetVol,
   SVec<double,3> &X, int nodeNum[4], int material_id
   );
  
  void computeDerivativeOfSurfaceTerm(int, Vec3D &, Vec3D &, double [3], double *, double *, double *[3], double *[3], double, double *);
  
  void computeDerivativeOfSurfaceTerm(double [4][3], double [4][3], int, Vec3D &, Vec3D &, double [4],
                                      double *, double *, double *[4], double *[4], double, double *);
  
  void computeSourceTerm(double dudxj[3][3],double dnudx[3],double d2wall, 
			 double *V,  double *S);
  
  void rstVar(IoData &ioData, Communicator *com) 
  { 
    rstVarNS(ioData, com); 
	  if (wallFcn)  wallFcn->rstVar(ioData, com);
    rstVarSA(ioData);
  }
  
  
  void computeBCsJacobianWallValues
  (
   int, Vec3D &, double [3], double *, double *, double *[3]
   );
  
  
  double computeDerivativeOfViscousTimeStep
  (
   double *, double *, double *, double *, double
   );
  
  bool withWallFcn() 
  { 
	  if(wallFcn) return true;
	  else        return false;
  };
  
};

//------------------------------------------------------------------------------

class FemEquationTermDES : public FemEquationTerm, public NavierStokesTerm, 
public DESTerm {
  
public:
  
  double x0,y0,z0,x1,y1,z1;
  double cdes;
  bool trip, usefv3;
  //map<int, PorousMedia *> &volInfo;
  double velocity, density, length;
  ConstantPrandtlThermalCondFcn turbThermalCondFcn;
  
private:
  void computeLaminarTransportCoefficients(const double T, 
		  double &mu, double &lambda, double &kappa);
  void computeTurbulentTransportCoefficients(double *V[], int nodeNum[], SVec<double,3> &X,
		  const double mul, const double lambdal, double kappal,
		  double &mutilde, double &mut, double &lambdat, double &kappat);

public:
  FemEquationTermDES(IoData &, VarFcn *);
  ~FemEquationTermDES() {}
  
  double computeViscousTimeStep(double *, double *);
  
  double max(double a, double b) { return (a>b) ? a : b; }
  double min(double a, double b) { return (a<b) ? a : b; }
  
  bool computeVolumeTerm(double [4][3], double [4], double *[4],
                         double *, double *, double *, double, SVec<double,3> &, int [4], int);
  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], double *, double *, double *, double,
                                 SVec<double,3> &, int [4], int);
  void computeSurfaceTerm(int, Vec3D &, double [3], 
                          double *, double *[3], double *);

  double computeNormDerivWallFcn(double rho, double T, double Du1, double DT1,
										 double d2w, double &dudn, double &dTdn);

  void computeJacobianSurfaceTerm(int, Vec3D &, double [3], double *, double *[3], double *);
  void computeSurfaceTerm(double [4][3], int, Vec3D &, double [4], 
                          double *, double *[4], double *);
  void computeJacobianSurfaceTerm(double [4][3], int, Vec3D &, double [4], 
                                  double *, double *[4], double *);
  bool doesSourceTermExist() { return true; }
  
  // Included (MB)
  /// UH (08/10) The following function results in exit (Not Implemented).
  bool computeDerivativeOfVolumeTerm
  (
   double [4][3], double [4][3], double [4], double *[4], double *[4],
   double, double *, double *, double *, double, SVec<double,3> &, int [4], int
   );
  
   void computeDerivativeOperatorsOfVolumeTerm // YC
   (double [4][3], double *[4], double (*)[5][4][3], double (*)[5][4][5], double (*)[5]) {
     fprintf(stderr, "*** Error: FemEquationTermDES::computeDerivativeOperatorsOfVolumeTerm is not implemented yet\n");
     exit(-1);
   }

  
  /// UH (08/10) The following function results in exit (Not Implemented).
  void computeDerivativeOfSurfaceTerm
  (
   int, Vec3D &, Vec3D &, double [3], double *, double *, 
   double *[3], double *[3], double, double *
   );
  
  
  /// UH (08/10) The following function results in exit (Not Implemented).
  void computeDerivativeOfSurfaceTerm
  (
   double [4][3], double [4][3], int, Vec3D &, Vec3D &, double [4],
   double *, double *, double *[4], double *[4], double, double *
   );
  
  
  void rstVar(IoData &ioData, Communicator *com)
  { 
    rstVarNS(ioData, com); 
    if (wallFcn) 
      wallFcn->rstVar(ioData, com);
    // Missing restart ?
    if (com->cpuNum() == 0)
      fprintf(stderr, "\n\n !!! WARNING DESTerm - Variables are not restarted !!!\n");
  }
  
  void computeBCsJacobianWallValues
  (
   int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]
   ) 
  {
    fprintf(stderr, "*** Error: computeBCsJacobianWallValues should not be called\n");
    exit(1);
  }
  
  
  double computeDerivativeOfViscousTimeStep
  (
   double *, double *, double *, double *, double
   );
  
  bool withWallFcn() 
  { 
	  if(wallFcn) return true;
	  else        return false;
  };
   
};

//------------------------------------------------------------------------------

class FemEquationTermSAmean : public FemEquationTerm, public NavierStokesTerm, 
public SATerm {
  
public:
  
  double x0,y0,z0,x1,y1,z1;
  bool trip;
  //map<int, PorousMedia *> &volInfo;
  double velocity, density, length;
  ConstantPrandtlThermalCondFcn turbThermalCondFcn;
  
private:
  void computeLaminarTransportCoefficients(const double T, 
		  double &mu, double &lambda, double &kappa);
  void computeTurbulentTransportCoefficients(double *V[], int nodeNum[], SVec<double,3> &X,
		  const double mul, const double lambdal, double kappal,
		  double &mutilde, double &mut, double &lambdat, double &kappat);

public:
  FemEquationTermSAmean(IoData &, VarFcn *);
  ~FemEquationTermSAmean() {}
  
  double computeViscousTimeStep(double *, double *);
  
  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], 
                                 double *, double *, double *, double, SVec<double,3> &, int [4], int);
  void computeJacobianSurfaceTerm(int, Vec3D &, double [3], 
                                  double *, double *[3], double *);
  void computeJacobianSurfaceTerm(double [4][3], int, Vec3D &, double [4], 
                                  double *, double *[4], double *);
  
  bool computeVolumeTerm(double dp1dxj[4][3], double d2w[4], double *v[4],
                         double *r, double *s, double *, double, 
                         SVec<double,3> &, int [4], int) {
    fprintf(stderr, "*** Error: computeVolumeTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(int c, Vec3D &n, double d2w[3], 
                          double *vw, double *v[3], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
                          double *vw, double *v[4], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  
  double computeNormDerivWallFcn(double rho, double T, double Du1, double DT1,
										 double d2w, double &dudn, double &dTdn)
  {
    fprintf(stderr, "*** Error: computeNormDerivWallFcn should not be called\n");
    exit(1);
  }
  
  // Included (MB)
  bool computeDerivativeOfVolumeTerm
  (
   double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4], 
   double *v[4], double *dv[4],
   double dMach, double *dr, double *ds, double *dpr, 
   double dtetvol, SVec<double,3> &x, int nodesnum[4], int volid
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermSAmean::computeDerivativeVolumeTerm should not be called\n");
    exit(1);
  }
  
  void computeDerivativeOperatorsOfVolumeTerm // YC
  (double [4][3], double *[4], double (*)[5][4][3], double (*)[5][4][5], double (*)[5]) {
    fprintf(stderr, "*** Error: FemEquationTermSAmean::computeDerivativeOperatorsOfVolumeTerm is not implemented yet\n");
    exit(-1);
  }
  
  void computeDerivativeOfSurfaceTerm
  (
   int c, Vec3D &n, Vec3D &dn, double d2w[3],
   double *vw, double *dvw, double *v[3], double *dv[3], double dMach, double *dr
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermSAmean::computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }
  
  
  void computeDerivativeOfSurfaceTerm
  (
   double dp1dxj[4][3], double ddp1dxj[4][3], int c, Vec3D &n, Vec3D &dn, 
   double d2w[4],
   double *vw, double *dvw, double *v[4], double *dv[4], double dMach, double *dr
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermSAmean::computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }
  
 
  /// \note (09/10)
  /// This function is only called with weak turbulence model coupling
  /// and with the exact matrix vector product for the Jacobian (H2). 
  /// The viscous part is contained in the H1 Jacobian component. 
  /// So no operation is performed.
  void computeBCsJacobianWallValues
  (
   int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]
   ) 
  { };


  void rstVar(IoData &ioData, Communicator *com) 
  {
    fprintf(stderr, "*** Error: FemEquationTermSAmean::rstVar should not be called\n");
    exit(1);
  }
  
  double computeDerivativeOfViscousTimeStep
  (
   double *, double *, double *, double *, double
   );
  
  bool withWallFcn() 
  { 
	  if(wallFcn) return true;
	  else        return false;
  };
  
};

//------------------------------------------------------------------------------
class FemEquationTermDESmean : public FemEquationTerm, public NavierStokesTerm, 
public DESTerm {
  
public:
  
  double x0,y0,z0,x1,y1,z1;
  bool trip;
  //map<int, PorousMedia *> &volInfo;
  double velocity, density, length;
  ConstantPrandtlThermalCondFcn turbThermalCondFcn;

private:
  void computeLaminarTransportCoefficients(const double T, 
		  double &mu, double &lambda, double &kappa);
  void computeTurbulentTransportCoefficients(double *V[], int nodeNum[], SVec<double,3> &X,
		  const double mul, const double lambdal, double kappal,
		  double &mutilde, double &mut, double &lambdat, double &kappat);

public:
  FemEquationTermDESmean(IoData &, VarFcn *);
  ~FemEquationTermDESmean() {}
  
  double computeViscousTimeStep(double *, double *);
  
  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], 
                                 double *, double *, double *, double, SVec<double,3> &, int [4], int);
  void computeJacobianSurfaceTerm(int, Vec3D &, double [3], 
                                  double *, double *[3], double *);
  void computeJacobianSurfaceTerm(double [4][3], int, Vec3D &, double [4], 
                                  double *, double *[4], double *);
  
  bool computeVolumeTerm(double dp1dxj[4][3], double d2w[4], double *v[4],
                         double *r, double *s, double *, double, 
                         SVec<double,3> &, int [4], int) {
    fprintf(stderr, "*** Error: computeVolumeTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(int c, Vec3D &n, double d2w[3], 
                          double *vw, double *v[3], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
                          double *vw, double *v[4], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  
  double computeNormDerivWallFcn(double rho, double T, double Du1, double DT1,
										 double d2w, double &dudn, double &dTdn)
  {
    fprintf(stderr, "*** Error: computeNormDerivWallFcn should not be called\n");
    exit(1);
  }
  
  // Included (MB)
  bool computeDerivativeOfVolumeTerm
  (
   double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4], 
   double *v[4], double *dv[4],
   double dMach, double *dr, double *ds, double *dpr, double dtetvol, 
   SVec<double,3> &x, int nodesnum[4], int volid
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermDESmean::computeDerivativeVolumeTerm should not be called\n");
    exit(1);
  }
  
  void computeDerivativeOperatorsOfVolumeTerm // YC
  (
   double [4][3], double *[4],
   double (*)[5][4][3], double (*)[5][4][5], double (*)[5]) {
    fprintf(stderr, "*** Error: FemEquationTermDESmean::computeDerivativeOperatorsOfVolumeTerm is not implemented yet\n");
   exit(-1);
  }

  
  void computeDerivativeOfSurfaceTerm
  (
   int c, Vec3D &n, Vec3D &dn, double d2w[3],
   double *vw, double *dvw, double *v[3], double *dv[3], double dMach, double *dr
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermDESmean::computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }
  
  
  void computeDerivativeOfSurfaceTerm
  (
   double dp1dxj[4][3], double ddp1dxj[4][3], int c, 
   Vec3D &n, Vec3D &dn, double d2w[4],
   double *vw, double *dvw, double *v[4], double *dv[4], double dMach, double *dr
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermDESmean::computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }
  
  
  void rstVar(IoData &ioData, Communicator *com) {
    fprintf(stderr, "*** Error: FemEquationTermDESmean::rstVar should not be called\n");
    exit(1);
  }
  
  void computeBCsJacobianWallValues
  (
   int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermDESmean::computeBCsJacobianWallValues should not be called\n");
    exit(1);
  }
  
  
  double computeDerivativeOfViscousTimeStep
  (
   double *, double *, double *, double *, double
   );
  
  bool withWallFcn() 
  { 
	  if(wallFcn) return true;
	  else        return false;
  };
  
};

//------------------------------------------------------------------------------

class FemEquationTermSAturb : public FemEquationTerm, public NavierStokesTerm,
public SATerm {
  
public:
  
  double x0,y0,z0,x1,y1,z1;
  bool trip;
  //map<int, PorousMedia *> &volInfo;
  
  FemEquationTermSAturb(IoData &, VarFcn *);
  ~FemEquationTermSAturb() {}
  
  double computeViscousTimeStep(double *, double *){
    fprintf(stderr, "*** Error: computeViscousTimeStep should not be called in FemSAturb\n");
    exit(1);
  }
  
  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], 
                                 double *, double *, double *, double, SVec<double,3> &, int [4], int);
  bool doesFaceTermExist(int code) { return false; }
  bool doesFaceNeedGradientP1Function() { return false; }
  bool doesSourceTermExist() { return true; }
  
  bool computeVolumeTerm(double dp1dxj[4][3], double d2w[4], double *v[4],
                         double *r, double *s, double *, double, 
                         SVec<double,3> &, int [4], int) {
    fprintf(stderr, "*** Error: computeVolumeTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(int c, Vec3D &n, double d2w[3], 
                          double *vw, double *v[3], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeJacobianSurfaceTerm(int c, Vec3D &n, double d2w[3], 
                                  double *vw, double *v[3], double *drdu) {
    fprintf(stderr, "*** Error: computeJacobianSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
                          double *vw, double *v[4], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }

  double computeNormDerivWallFcn(double rho, double T, double Du1, double DT1,
										 double d2w, double &dudn, double &dTdn)
  {
    fprintf(stderr, "*** Error: computeNormDerivWallFcn should not be called\n");
    exit(1);
  }

  void computeJacobianSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
                                  double *vw, double *v[4], double *drdu) {
    fprintf(stderr, "*** Error: computeJacobianSurfaceTerm should not be called\n");
    exit(1);
  }
  
  // Included (MB)
  bool computeDerivativeOfVolumeTerm
  (
   double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4], 
   double *v[4], double *dv[4],
   double dMach, double *dr, double *ds, double *dpr, double dtetvol, 
   SVec<double,3> &x, int nodesnum[4], int volid
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermSAturb::computeDerivativeOfVolumeTerm should not be called\n");
    exit(1);
  }
  
  void computeDerivativeOperatorsOfVolumeTerm // YC
  (
    double [4][3], double *[4], double (*)[5][4][3], double (*)[5][4][5], double (*)[5]) {
    fprintf(stderr, "*** Error: FemEquationTermSAturb::computeDerivativeOperatorsOfVolumeTerm is not implemented yet\n");
    exit(-1);
  }

  void computeDerivativeOfSurfaceTerm
  (
   int c, Vec3D &n, Vec3D &dn, double d2w[3],
   double *vw, double *dvw, double *v[3], double *dv[3], double dMach, double *dr
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermSAturb::computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }
  
  
  void computeDerivativeOfSurfaceTerm
  (
   double dp1dxj[4][3], double ddp1dxj[4][3], int c, Vec3D &n, Vec3D &dn, 
   double d2w[4],
   double *vw, double *dvw, double *v[4], double *dv[4], double dMach, double *dr
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermSAturb::computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }
  
  
  void rstVar(IoData &ioData, Communicator *com) {
    fprintf(stderr, "*** Error: FemEquationTermSAturb::rstVar should not be called\n");
    exit(1);
  }
  
  void computeBCsJacobianWallValues
  (
   int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermSAturb::computeBCsJacobianWallValues should not be called\n");
    exit(1);
  }
  
  
  double computeDerivativeOfViscousTimeStep
  (
   double *, double *, double *, double *, double
   )
  {
    fprintf(stderr, "*** Error: FemEquationTermSAturb::computeViscousDerivativeOfViscousTimeStep should not be called\n");
    exit(1);
  }
  
  bool withWallFcn() 
  { 
	  if(wallFcn) return true;
	  else        return false;
  };

};

//------------------------------------------------------------------------------
class FemEquationTermDESturb : public FemEquationTerm, public NavierStokesTerm,
public DESTerm {
  
public:
  
  double x0,y0,z0,x1,y1,z1;
  bool trip;
  //map<int, PorousMedia *> &volInfo;
  
  FemEquationTermDESturb(IoData &, VarFcn *);
  ~FemEquationTermDESturb() {}
  
  double computeViscousTimeStep(double *, double *){
    fprintf(stderr, "*** Error: computeViscousTimeStep should not be called in FemDESturb\n");
    exit(1);
  }
  
  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], 
                                 double *, double *, double *, double, SVec<double,3> &, int [4], int);
  bool doesFaceTermExist(int code) { return false; }
  bool doesFaceNeedGradientP1Function() { return false; }
  bool doesSourceTermExist() { return true; }
  
  bool computeVolumeTerm(double dp1dxj[4][3], double d2w[4], double *v[4],
                         double *r, double *s, double *, 
                         double, SVec<double,3> &, int [4], int) {
    fprintf(stderr, "*** Error: computeVolumeTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(int c, Vec3D &n, double d2w[3], 
                          double *vw, double *v[3], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeJacobianSurfaceTerm(int c, Vec3D &n, double d2w[3], 
                                  double *vw, double *v[3], double *drdu) {
    fprintf(stderr, "*** Error: computeJacobianSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
                          double *vw, double *v[4], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeJacobianSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
                                  double *vw, double *v[4], double *drdu) {
    fprintf(stderr, "*** Error: computeJacobianSurfaceTerm should not be called\n");
    exit(1);
  }
  
  double computeNormDerivWallFcn(double rho, double T, double Du1, double DT1,
										 double d2w, double &dudn, double &dTdn)
  {
    fprintf(stderr, "*** Error: computeNormDerivWallFcn should not be called\n");
    exit(1);
  }
  
  // Included (MB)
  bool computeDerivativeOfVolumeTerm
  (
   double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4], 
   double *v[4], double *dv[4],
   double dMach, double *dr, double *ds, double *dpr, double dtetvol, 
   SVec<double,3> &x, int nodesnum[4], int volid
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermDESturb::computeDerivativeOfVolumeTerm should not be called\n");
    exit(1);
  }

  void computeDerivativeOperatorsOfVolumeTerm // YC
  (
   double [4][3], double *[4], double (*)[5][4][3], double (*)[5][4][5], double (*)[5]) {
   fprintf(stderr, "*** Error: FemEquationTermDESturb::computeDerivativeOperatorsOfVolumeTerm is not implemented yet\n");
   exit(-1);
  }
  
  void computeDerivativeOfSurfaceTerm
  (
   int c, Vec3D &n, Vec3D &dn, double d2w[3],
   double *vw, double *dvw, double *v[3], double *dv[3], double dMach, double *dr
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermDESturb::computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }
  
  
  void computeDerivativeOfSurfaceTerm
  (
   double dp1dxj[4][3], double ddp1dxj[4][3], int c, Vec3D &n, Vec3D &dn, 
   double d2w[4],
   double *vw, double *dvw, double *v[4], double *dv[4], double dMach, double *dr
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermDESturb::computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }
  
  
  void rstVar(IoData &ioData, Communicator *com) {
    fprintf(stderr, "*** Error: FemEquationTermDESturb::rstVar should not be called\n");
    exit(1);
  }
  
  
  void computeBCsJacobianWallValues
  (
   int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermDESturb::computeBCsJacobianWallValues should not be called\n");
    exit(1);
  }
  
  
  double computeDerivativeOfViscousTimeStep
  (
   double *, double *, double *, double *, double
   )
  {
    fprintf(stderr, "*** Error: FemEquationTermDESturb::computeDerivativeOfViscousTimeStep should not be called\n");
    exit(1);
  }
  
  bool withWallFcn() 
  { 
	  if(wallFcn) return true;
	  else        return false;
  };
  
};

//------------------------------------------------------------------------------

class FemEquationTermKE : public FemEquationTerm, public NavierStokesTerm, 
public KEpsilonTerm {
  
public:
  
  double x0,y0,z0,x1,y1,z1;
  bool trip;
  //map<int, PorousMedia *> &volInfo;
  double velocity, density, length;
  ConstantPrandtlThermalCondFcn turbThermalCondFcn;
  
private:
  void computeLaminarTransportCoefficients(const double T, 
		  double &mu, double &lambda, double &kappa);
  void computeTurbulentTransportCoefficients(double *V[], int nodeNum[], SVec<double,3> &X,
		  const double mul, const double lambdal, double kappal,
		  double &rhok, double &rhoeps, double &mut, double &lambdat, double &kappat);

public:
  FemEquationTermKE(IoData &, VarFcn *);
  ~FemEquationTermKE() {}
  
  double computeViscousTimeStep(double *, double *);
  
  bool computeVolumeTerm(double [4][3], double [4], double *[4],
                         double *, double *, double *, double, 
                         SVec<double,3> &, int [4], int);
  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], double *, double *,
                                 double *, double, SVec<double,3> &, int [4], int);
  void computeSurfaceTerm(int, Vec3D &, double [3], 
                          double *, double *[3], double *);
  void computeJacobianSurfaceTerm(int, Vec3D &, double [3], double *, double *[3], double *);
  bool doesSourceTermExist() { return true; }
  
  void computeSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
                          double *vw, double *v[4], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeJacobianSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
                                  double *vw, double *v[4], double *drdu) {
    fprintf(stderr, "*** Error: computeJacobianSurfaceTerm should not be called\n");
    exit(1);
  }
  
  double computeNormDerivWallFcn(double rho, double T, double Du1, double DT1,
										 double d2w, double &dudn, double &dTdn);
  
  // Included (MB)
  bool computeDerivativeOfVolumeTerm
  (
   double [4][3], double [4][3], double [4], double *[4], double *[4],
   double, double *, double *, double *, double, SVec<double,3> &, int [4], int
   );
  
  void computeDerivativeOperatorsOfVolumeTerm // YC WTF is this code style?
  (
   double [4][3], double *[4], double (*)[5][4][3], double (*)[5][4][5], double (*)[5]) {
    fprintf(stderr, "*** Error: FemEquationTermKE::computeDerivativeOperatorsOfVolumeTerm is not implemented yet\n");
    exit(-1);
  }

  
  void computeDerivativeOfSurfaceTerm
  (
   int c, Vec3D &, Vec3D &, double [3],
   double *, double *, double *[3], double *[3], double, double *
   );
  
  
  void computeDerivativeOfSurfaceTerm
  (
   double dp1dxj[4][3], double ddp1dxj[4][3], int c, Vec3D &n, Vec3D &dn, 
   double d2w[4],
   double *vw, double *dvw, double *v[4], double *dv[4], double dMach, double *dr
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermKE::computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }
  
  void rstVar(IoData &ioData, Communicator *com) 
  { 
    rstVarNS(ioData, com); 
    if (wallFcn) wallFcn->rstVar(ioData, com);
    if (com->cpuNum() == 0)
      fprintf(stderr, "\n\n !!! WARNING: kEpsilon - Variables are not restarted !!!\n");
    //----------
  }
  
  void computeBCsJacobianWallValues
  (
   int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermKE::computeBCsJacobianWallValues should not be called\n");
    exit(1);
  }
  
  double computeDerivativeOfViscousTimeStep
  (
   double *, double *, double *, double *, double
   );
  
  bool withWallFcn() 
  { 
	  if(wallFcn) return true;
	  else        return false;
  };
  
};

//------------------------------------------------------------------------------

class FemEquationTermKEmean : public FemEquationTerm, public NavierStokesTerm, 
public KEpsilonTerm {
  
public:
  
  double x0,y0,z0,x1,y1,z1;
  bool trip;
  //map<int, PorousMedia *> &volInfo;
  double velocity, density, length;
  ConstantPrandtlThermalCondFcn turbThermalCondFcn;
  
private:
  void computeLaminarTransportCoefficients(const double T, 
		  double &mu, double &lambda, double &kappa);
  void computeTurbulentTransportCoefficients(double *V[], int nodeNum[], SVec<double,3> &X,
		  const double mul, const double lambdal, double kappal,
		  double &rhok, double &rhoeps, double &mut, double &lambdat, double &kappat);

public:
  FemEquationTermKEmean(IoData &, VarFcn *);
  ~FemEquationTermKEmean() {}
  
  double computeViscousTimeStep(double *, double *);
  
  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], 
                                 double *, double *, double *, double, SVec<double,3> &, int [4], int);
  void computeJacobianSurfaceTerm(int, Vec3D &, double [3], 
                                  double *, double *[3], double *);
  
  bool computeVolumeTerm(double dp1dxj[4][3], double d2w[4], double *v[4],
                         double *r, double *s, double *, double, 
                         SVec<double,3> &, int [4], int) {
    fprintf(stderr, "*** Error: computeVolumeTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(int c, Vec3D &n, double d2w[3], 
                          double *vw, double *v[3], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
                          double *vw, double *v[4], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeJacobianSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
                                  double *vw, double *v[4], double *drdu) {
    fprintf(stderr, "*** Error: computeJacobianSurfaceTerm should not be called\n");
    exit(1);
  }
  
  // Included (MB)
  bool computeDerivativeOfVolumeTerm
  (
   double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4], 
   double *v[4], double *dv[4],
   double dMach, double *dr, double *ds, double *dpr, double dtetvol, 
   SVec<double,3> &x, int nodesnum[4], int volid
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermKEmean::computeDerivativeOfVolumeTerm should not be called\n");
    exit(1);
  }
  
  void computeDerivativeOperatorsOfVolumeTerm // YC
   (
    double [4][3], double *[4], double (*)[5][4][3], double (*)[5][4][5], double (*)[5]) {
    fprintf(stderr, "*** Error: FemEquationTermKEmean::computeDerivativeOperatorsOfVolumeTerm is not implemented yet\n");
    exit(-1);
  }

  void computeDerivativeOfSurfaceTerm
  (
   int c, Vec3D &n, Vec3D &dn, double d2w[3],
   double *vw, double *dvw, double *v[3], double *dv[3], double dMach, double *dr
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermKEmean::computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }
  
  double computeNormDerivWallFcn(double rho, double T, double Du1, double DT1,
										 double d2w, double &dudn, double &dTdn)
  {
    fprintf(stderr, "*** Error: computeNormDerivWallFcn should not be called\n");
    exit(1);
  }
  
  void computeDerivativeOfSurfaceTerm
  (
   double dp1dxj[4][3], double ddp1dxj[4][3], int c, Vec3D &n, Vec3D &dn, 
   double d2w[4],
   double *vw, double *dvw, double *v[4], double *dv[4], double dMach, double *dr
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermKEmean::computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }
  
  
  void rstVar(IoData &ioData, Communicator *com) {
    fprintf(stderr, "*** Error: FemEquationTermKEmean::rstVar should not be called\n");
    exit(1);
  }
  
  
  void computeBCsJacobianWallValues
  (
   int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermKEmean::computeBCsJacobianWallValues should not be called\n");
    exit(1);
  }
  
  
  double computeDerivativeOfViscousTimeStep
  (
   double *, double *, double *, double *, double
   );
  
  bool withWallFcn() 
  { 
	  if(wallFcn) return true;
	  else        return false;
  };

};

//------------------------------------------------------------------------------

class FemEquationTermKEturb : public FemEquationTerm, public NavierStokesTerm,
public KEpsilonTerm {
  
public:
  
  double x0,y0,z0,x1,y1,z1;
  bool trip;
  //map<int, PorousMedia *> &volInfo;
  
  FemEquationTermKEturb(IoData &, VarFcn *);
  ~FemEquationTermKEturb() {}
  
  double computeViscousTimeStep(double *, double *){
    fprintf(stderr, "*** Error: computeViscousTimeStep should not be called in FemKEturb\n");
    exit(1);
  }
  
  bool computeJacobianVolumeTerm(double [4][3], double [4], double *[4], 
                                 double *, double *, double *, double, SVec<double,3> &, int [4], int);
  bool doesFaceTermExist(int code) { return false; }
  bool doesFaceNeedGradientP1Function() { return false; }
  bool doesSourceTermExist() { return true; }
  
  bool computeVolumeTerm(double dp1dxj[4][3], double d2w[4], double *v[4],
                         double *r, double *s, double *, double , 
                         SVec<double,3> &, int [4], int) {
    fprintf(stderr, "*** Error: computeVolumeTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(int c, Vec3D &n, double d2w[3], 
                          double *vw, double *v[3], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeJacobianSurfaceTerm(int c, Vec3D &n, double d2w[3], 
                                  double *vw, double *v[3], double *drdu) {
    fprintf(stderr, "*** Error: computeJacobianSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
                          double *vw, double *v[4], double *r) {
    fprintf(stderr, "*** Error: computeSurfaceTerm should not be called\n");
    exit(1);
  }
  void computeJacobianSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
                                  double *vw, double *v[4], double *drdu) {
    fprintf(stderr, "*** Error: computeJacobianSurfaceTerm should not be called\n");
    exit(1);
  }
  
  // Included (MB)
  bool computeDerivativeOfVolumeTerm
  (
   double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4], 
   double *v[4], double *dv[4],
   double dMach, double *dr, double *ds, double *dpr, double dtetvol, 
   SVec<double,3> &x, int nodesnum[4], int volid
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermKEturb::computeDerivativeOfVolumeTerm should not be called\n");
    exit(1);
  }
  
  void computeDerivativeOperatorsOfVolumeTerm // YC
  (
   double [4][3], double *[4], double (*)[5][4][3], double (*)[5][4][5], double (*)[5]) {
    fprintf(stderr, "*** Error: FemEquationTermKEturb::computeDerivativeOperatorsOfVolumeTerm is not implemented yet\n");
    exit(-1);
  }

  
  void computeDerivativeOfSurfaceTerm
  (
   int c, Vec3D &n, Vec3D &dn, double d2w[3],
   double *vw, double *dvw, double *v[3], double *dv[3], double dMach, double *dr
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermKEturb::computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }
  
  double computeNormDerivWallFcn(double rho, double T, double Du1, double DT1,
										 double d2w, double &dudn, double &dTdn)
  {
	  fprintf(stderr, "*** Error: computeNormDerivWallFcn should not be called\n");
	  exit(1);
  }
  
  void computeDerivativeOfSurfaceTerm
  (
   double dp1dxj[4][3], double ddp1dxj[4][3], int c, Vec3D &n, Vec3D &dn, 
   double d2w[4],
   double *vw, double *dvw, double *v[4], double *dv[4], double dMach, double *dr
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermKEturb::computeDerivativeOfSurfaceTerm should not be called\n");
    exit(1);
  }
  
  
  void rstVar(IoData &ioData, Communicator *com) {
    fprintf(stderr, "*** Error: FemEquationTermKEturb::rstVar should not be called\n");
    exit(1);
  }
  
  
  void computeBCsJacobianWallValues
  (
   int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]
   ) 
  {
    fprintf(stderr, "*** Error: FemEquationTermKEturb::computeBCsJacobianWallValues should not be called\n");
    exit(1);
  }
  
  
  double computeDerivativeOfViscousTimeStep
  (
   double *, double *, double *, double *, double
   )
  {
    fprintf(stderr, "*** Error: FemEquationTermKEturb::computeDerivativeOfViscousTimeStep should not be called\n");
    exit(1);
  }
  
  bool withWallFcn() 
  { 
	  if(wallFcn) return true;
	  else        return false;
  };

};

//------------------------------------------------------------------------------

#endif
