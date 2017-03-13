#ifndef _FEM_EQUATION_TERM_H_
#define _FEM_EQUATION_TERM_H_

#include <BcDef.h>
#include <WallFcn.h>
#include <Vector.h>

// Included
class Communicator;

struct Vec3D;

//------------------------------------------------------------------------------

class FemEquationTerm {

// Included (MB)
public:
  bool completeJac;

protected:

  WallFcn* wallFcn;
  map<int, PorousMedia *> volInfo;

public:

  FemEquationTerm(map<int, VolumeData *> &volData) {
    wallFcn = 0; 
    //construction of volInfo (map from id to porousMedia)
    //...loop on all the VolumeData and check which one is a PorousMedia...
    map<int, VolumeData *>::iterator it;
    if(!volData.empty()){
      for(it=volData.begin(); it!=volData.end(); it++){
        //...if it is a PorousMedia, add it to volInfo
        if(it->second->type == VolumeData::POROUS){
          //...check if it already exists...
          map<int, PorousMedia *>::iterator pmit = volInfo.find(it->first);
          //...otherwise add it...
          if(pmit == volInfo.end()) 
            volInfo[it->first] = &(it->second->porousMedia);
        }
      }
    }
  }

  virtual ~FemEquationTerm() { if (wallFcn) delete wallFcn; }

  virtual double computeViscousTimeStep(double *, double *) = 0;

  virtual bool computeVolumeTerm(double dp1dxj[4][3], double d2w[4], double *v[4],
				 double *r, double *s, double *, double, SVec<double,3> &, int [4], int) = 0;
  virtual bool computeJacobianVolumeTerm(double dp1dxj[4][3], double d2w[4], double *v[4], 
					 double *drdu, double *dsdu, double *dpdu, double, SVec<double,3> &, int [4], int) = 0;
  virtual void computeSurfaceTerm(int c, Vec3D &n, double d2w[3], 
				  double *vw, double *v[3], double *r) = 0;
  virtual void computeJacobianSurfaceTerm(int c, Vec3D &n, double d2w[3], 
					  double *vw, double *v[3], double *drdu) = 0;
  virtual void computeSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
				  double *vw, double *v[4], double *r) = 0;
  virtual void computeJacobianSurfaceTerm(double dp1dxj[4][3], int c, Vec3D &n, double d2w[4], 
					  double *vw, double *v[4], double *drdu) = 0;

  virtual double computeNormDerivWallFcn(double rho, double T, 
													double Du1, double DT1,	double d2w, 
													double &dudn, double &dTdn) = 0;

  virtual bool doesFaceTermExist(int code) {
    if (wallFcn)
      return (code == BC_ADIABATIC_WALL_MOVING || 
	      code == BC_ADIABATIC_WALL_FIXED ||
	      code == BC_ISOTHERMAL_WALL_MOVING || 
	      code == BC_ISOTHERMAL_WALL_FIXED) ? true : false; 
    else
      return (code == BC_ADIABATIC_WALL_MOVING || 
	      code == BC_ADIABATIC_WALL_FIXED) ? true : false; 
  }

  virtual bool doesFaceNeedGradientP1Function() {
    if (wallFcn) return false;
    else return true;
  }

  virtual bool doesSourceTermExist() { return false; }

  virtual bool withWallFcn() {return false; };

// Included (MB)
  virtual bool computeDerivativeOfVolumeTerm(double dp1dxj[4][3], double ddp1dxj[4][3], double d2w[4], double *v[4],
				 double *dv[4], double dMach, double *dr, double *ds, double *dpr, double dtetvol, SVec<double,3> &x, int nodenum[4], int volid) = 0;
  virtual void computeDerivativeOperatorsOfVolumeTerm(double dp1dxj[4][3], double *v[4],
  			 double (*drddp1dxj)[5][4][3], double (*drdV)[5][4][5], double (*drdMach)[5]) = 0;
  virtual void computeDerivativeOfSurfaceTerm(int c, Vec3D &n, Vec3D &dn, double d2w[3],
				  double *vw, double *dvw, double *v[3], double *dv[3], double dMach, double *dr) = 0;
  virtual void computeDerivativeOfSurfaceTerm(double dp1dxj[4][3], double ddp1dxj[4][3], int c, Vec3D &n, Vec3D &dn, double d2w[4],
				  double *vw, double *dvw, double *v[4], double *dv[4], double dMach, double *dr) = 0;
  virtual void rstVar(IoData &ioData, Communicator *com) = 0;
  virtual void computeBCsJacobianWallValues(int c, Vec3D &n, double d2w[3], double *vw, double *dvw, double *v[3]) = 0;

  virtual double computeDerivativeOfViscousTimeStep(double *, double *, double *, double *, double) = 0;


  void computePermittivityTensor(double alpha[3], double beta[3], double ucg[3], double R[9], double *K)
  {
   
    double onehalf = 1.0/2.0;
    double norm_u = pow(ucg[0]*ucg[0] + ucg[1]*ucg[1] + ucg[2]*ucg[2], onehalf); 
    
    double diag[3] = {alpha[0]*norm_u + beta[0], 
                      alpha[1]*norm_u + beta[1],
                      alpha[2]*norm_u + beta[2]};

    double Rt_diag[9] = { R[0]*diag[0], R[3]*diag[1], R[6]*diag[2],
                          R[1]*diag[0], R[4]*diag[1], R[7]*diag[2],
                          R[2]*diag[0], R[5]*diag[1], R[8]*diag[2] };   

    for(int i=0; i<9; ++i)  K[i] = 0.0;
 
    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
         for (int k=0; k<3; ++k)
             K[3*i+j] += Rt_diag[3*i+k]*R[3*k+j];

  }
  // for porous modelling
  void computeGradPermittivityTensor(double alpha[3], double ucg[3], double R[9] ,double *B)
  {

    double onehalf = 1.0/2.0;
    double onefourth = 1.0/4.0;

    double norm_u = pow(ucg[0]*ucg[0] + ucg[1]*ucg[1] + ucg[2]*ucg[2], onehalf);
    double inv_norm = 1.0/norm_u;
                                                                                                                                                         
    double diag[3] = {onefourth*inv_norm*alpha[0], 
                      onefourth*inv_norm*alpha[1], 
                      onefourth*inv_norm*alpha[2]};
                                                                                                                                                         
    double Rt_diag[9] = { R[0]*diag[0], R[3]*diag[1], R[6]*diag[2],
                          R[1]*diag[0], R[4]*diag[1], R[7]*diag[2],
                          R[2]*diag[0], R[5]*diag[1], R[8]*diag[2] };
                                                                                                                                                         
    for(int i=0; i<9; ++i)  B[i] = 0.0;
                                                                                                                                                         
    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
         for (int k=0; k<3; ++k)
             B[3*i+j] += Rt_diag[3*i+k]*R[3*k+j];
    
  }
  
  // for porous modelling
  double computePorousTurbulentViscosity(
      map<int,PorousMedia *>::iterator it, double up[3], double length)
  {
    double cmu = 0.09;
    double coeff = 1.2247*pow(cmu,0.25);
    double Idr = it->second->idr;                    // average turbulence intensity
    double Ldr = it->second->ldr/length;             // 0.1*characteristic passage dimension
    double vel = sqrt(up[0]*up[0] + up[1]*up[1] + up[2]*up[2]);

    return coeff*Idr*Ldr*vel;
  }

  // for porous modelling
  double computeSecondPorousTurbulentViscosity(double lambdal, double mul, double mut)
  {
    //simple model that remains true when the Stokes' hypothesis is assumed
    return -2.0*mut/3.0;
  }

  // for porous modelling
  bool computeVolumeTermPorousCore(double tetVol,
      map<int,PorousMedia *>::iterator it,
      double length, double density, double velocity,
      double up[3], double *V[], double *PR){

    double RR[9], K[9];
    double alpha[3], beta[3];
    double volten = tetVol *(1.0/10.0);

    // non-dimensionalization //

    alpha[0] = it->second->alphax*length/density;
    alpha[1] = it->second->alphay*length/density;
    alpha[2] = it->second->alphaz*length/density;

    beta[0]  = it->second->betax*length/(density*velocity);
    beta[1]  = it->second->betay*length/(density*velocity);
    beta[2]  = it->second->betaz*length/(density*velocity);

    // transformation matrix //

    RR[0] = it->second->iprimex; RR[1] = it->second->iprimey; RR[2] = it->second->iprimez;
    RR[3] = it->second->jprimex; RR[4] = it->second->jprimey; RR[5] = it->second->jprimez;
    RR[6] = it->second->kprimex; RR[7] = it->second->kprimey; RR[8] = it->second->kprimez;

    // permittivity matrix  //

    computePermittivityTensor(alpha, beta, up, RR, K);

    double SS[4][3];

    SS[0][0] = volten * (V[0][1] + 0.5 *(V[1][1] + V[2][1] + V[3][1]));
    SS[0][1] = volten * (V[0][2] + 0.5 *(V[1][2] + V[2][2] + V[3][2]));
    SS[0][2] = volten * (V[0][3] + 0.5 *(V[1][3] + V[2][3] + V[3][3]));

    SS[1][0] = volten * (V[1][1] + 0.5 *(V[0][1] + V[2][1] + V[3][1]));
    SS[1][1] = volten * (V[1][2] + 0.5 *(V[0][2] + V[2][2] + V[3][2]));
    SS[1][2] = volten * (V[1][3] + 0.5 *(V[0][3] + V[2][3] + V[3][3]));

    SS[2][0] = volten * (V[2][1] + 0.5 *(V[1][1] + V[0][1] + V[3][1]));
    SS[2][1] = volten * (V[2][2] + 0.5 *(V[1][2] + V[0][2] + V[3][2]));
    SS[2][2] = volten * (V[2][3] + 0.5 *(V[1][3] + V[0][3] + V[3][3]));

    SS[3][0] = volten * (V[3][1] + 0.5 *(V[1][1] + V[2][1] + V[0][1]));
    SS[3][1] = volten * (V[3][2] + 0.5 *(V[1][2] + V[2][2] + V[0][2]));
    SS[3][2] = volten * (V[3][3] + 0.5 *(V[1][3] + V[2][3] + V[0][3]));

    // FE flux for the porous sink term //

    for (int j=0; j<4; ++j) {
      for (int k=0; k<3; ++k)
        PR[3*j+k] += (K[3*k+0] * SS[j][0] + K[3*k+1] * SS[j][1] + K[3*k+2] * SS[j][2]);
    }

    return true;
  }

  // for porous modelling
  template<int dim>
  bool computeJacobianVolumeTermPorousCore(double tetVol,
      map<int,PorousMedia *>::iterator it,
      double length, double density, double velocity,
      double up[3], double *V[], double (*dPdU)[4][dim][dim]){

    double RR[9], K[9], B[9];
    double alpha[3], beta[3];

    double volten = tetVol *(1.0/10.0);

    // non-dimensionalization correction

    alpha[0] = it->second->alphax*length/density;
    alpha[1] = it->second->alphay*length/density;
    alpha[2] = it->second->alphaz*length/density;

    beta[0]  = it->second->betax*length/(density*velocity);
    beta[1]  = it->second->betay*length/(density*velocity);
    beta[2]  = it->second->betaz*length/(density*velocity);

    // transformation matrix
    RR[0] = it->second->iprimex; RR[1] = it->second->iprimey; RR[2] = it->second->iprimez;
    RR[3] = it->second->jprimex; RR[4] = it->second->jprimey; RR[5] = it->second->jprimez;
    RR[6] = it->second->kprimex; RR[7] = it->second->kprimey; RR[8] = it->second->kprimez;

    // permittivity matrix
    computePermittivityTensor(alpha, beta, up, RR, K);

    // gradient of permittivity matrix
    computeGradPermittivityTensor(alpha, up, RR, B);

    double SS[4][3];

    SS[0][0] = volten * (V[0][1] + 0.5 *(V[1][1] + V[2][1] + V[3][1]));
    SS[0][1] = volten * (V[0][2] + 0.5 *(V[1][2] + V[2][2] + V[3][2]));
    SS[0][2] = volten * (V[0][3] + 0.5 *(V[1][3] + V[2][3] + V[3][3]));

    SS[1][0] = volten * (V[1][1] + 0.5 *(V[0][1] + V[2][1] + V[3][1]));
    SS[1][1] = volten * (V[1][2] + 0.5 *(V[0][2] + V[2][2] + V[3][2]));
    SS[1][2] = volten * (V[1][3] + 0.5 *(V[0][3] + V[2][3] + V[3][3]));

    SS[2][0] = volten * (V[2][1] + 0.5 *(V[1][1] + V[0][1] + V[3][1]));
    SS[2][1] = volten * (V[2][2] + 0.5 *(V[1][2] + V[0][2] + V[3][2]));
    SS[2][2] = volten * (V[2][3] + 0.5 *(V[1][3] + V[0][3] + V[3][3]));

    SS[3][0] = volten * (V[3][1] + 0.5 *(V[1][1] + V[2][1] + V[0][1]));
    SS[3][1] = volten * (V[3][2] + 0.5 *(V[1][2] + V[2][2] + V[0][2]));
    SS[3][2] = volten * (V[3][3] + 0.5 *(V[1][3] + V[2][3] + V[0][3]));

    for (int k=0; k<4; ++k) {
      double BB[3];
      BB[0]  = B[0]*SS[k][0] +  B[1]*SS[k][1] +  B[2]*SS[k][2];
      BB[1]  = B[3]*SS[k][0] +  B[4]*SS[k][1] +  B[5]*SS[k][2];
      BB[2]  = B[6]*SS[k][0] +  B[7]*SS[k][1] +  B[8]*SS[k][2];

      double BV[9];
      BV[0] = up[0]*BB[0]; BV[1] = up[1]*BB[0]; BV[2] = up[2]*BB[0];
      BV[3] = up[0]*BB[1]; BV[4] = up[1]*BB[1]; BV[5] = up[2]*BB[1];
      BV[6] = up[0]*BB[2]; BV[7] = up[1]*BB[2]; BV[8] = up[2]*BB[2];

      for (int j=0; j<4; ++j) {
        double v[4] = {V[j][0], V[j][1], V[j][2], V[j][3]};

        double KU[25];
        multiplyBydVdU(v, K, KU, volten);

        double BU[25];
        multiplyBydVdU(v, BV, BU, 1.0);

        double cKU = (k==j) ? 1.0 : 0.5;

        for (int l=0; l<5; l++) {
          dPdU[k][j][0][l] = 0.0;
          dPdU[k][j][1][l] = cKU*KU[ 5+l] + BU[ 5+l];
          dPdU[k][j][2][l] = cKU*KU[10+l] + BU[10+l];
          dPdU[k][j][3][l] = cKU*KU[15+l] + BU[15+l];
          dPdU[k][j][4][l] = 0.0;
        }

        if(dim == 6 || dim == 7){
          for (int l=0; l<6; l++) {
            dPdU[k][j][5][l] = 0.0;
            dPdU[k][j][l][5] = 0.0;
          }
        }

      }
    }
    return true;

  }

  // for porous modelling
  void multiplyBydVdU(double V[4], double mat[9], double* res, double scaling) 
  {

   double rho = V[0];
   double invRho = 1.0/rho;
   double u = V[1];
   double v = V[2];
   double w = V[3];

   res[0] = 0.0;
   res[1] = 0.0;
   res[2] = 0.0;
   res[3] = 0.0;
   res[4] = 0.0;

   res[5] = -invRho*(mat[0]*u + mat[1]*v + mat[2]*w)*scaling;
   res[6] = mat[0]*invRho*scaling;
   res[7] = mat[1]*invRho*scaling;
   res[8] = mat[2]*invRho*scaling;
   res[9] = 0.0;

   res[10] = -invRho*(mat[3]*u + mat[4]*v + mat[5]*w)*scaling;
   res[11] = mat[3]*invRho*scaling;
   res[12] = mat[4]*invRho*scaling;
   res[13] = mat[5]*invRho*scaling;
   res[14] = 0.0;

   res[15] = -invRho*(mat[6]*u + mat[7]*v + mat[8]*w)*scaling;
   res[16] = mat[6]*invRho*scaling;
   res[17] = mat[7]*invRho*scaling;
   res[18] = mat[8]*invRho*scaling;
   res[19] = 0.0;

   res[20] = 0.0;
   res[21] = 0.0;
   res[22] = 0.0;
   res[23] = 0.0;
   res[24] = 0.0;

  }

/*
  void computeTranformationMatrix(Vec3D iprime, Vec3D jprime, Vec3D kprime, double *R)
  {
     Vec3D nx, ny, nz;

     nx.v[0] = 1.0; nx.v[1] = 0.0; nx.v[2] = 0.0;
     ny.v[0] = 0.0; ny.v[1] = 1.0; ny.v[2] = 0.0;
     nz.v[0] = 0.0; nz.v[1] = 0.0; nz.v[2] = 1.0;

     R[0] = nx*iprime; R[1] = ny*iprime; R[2] = nz*iprime;
     R[3] = nx*jprime; R[4] = ny*jprime; R[5] = nz*jprime;
     R[6] = nx*kprime; R[7] = ny*kprime; R[8] = nz*kprime;

    // to be used in pre-processing stage

     double onehalf = 1.0/2.0;
     Vec3D ny, nz, jprime, kprimei, temp1, temp2;
     ny.v[0] = 0.0; ny.v[1] = 1.0; ny.v[2] = 0.0;
     nz.v[0] = 0.0; nz.v[1] = 0.0; nz.v[2] = 1.0;
     
     temp1 = ny^iprime;
     double v1 = pow(temp1.v[0]*temp1.v[0] + temp1.v[1]*temp1.v[1] + temp3.v[2]*temp3.v[2], onehalf);

     temp2 = nz^iprime;
     double v2 = pow(temp2.v[0]*temp2.v[0] + temp2.v[1]*temp2.v[1] + temp3.v[2]*temp3.v[2], onehalf);
     
     if (v1 > v2) jprime = temp1;
     else jprime = temp2;

     kprime = iprime^jprime;
   
     R[0][0] = 

  }
*/      

};

//------------------------------------------------------------------------------

#endif
