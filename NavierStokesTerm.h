#ifndef _NAVIER_STOKES_TERM_H_
#define _NAVIER_STOKES_TERM_H_

#include <IoData.h>
#include <VarFcn.h>
#include <ViscoFcn.h>
#include <ThermalCondFcn.h>

// Included (MB)
class Communicator;

struct Vec3D;
//-----------------------------------------------------------------------------

class NavierStokesTerm {

protected:

  static const double third;
  static const double twothird;
  static const double fourth;

  double ooreynolds;
  double ooreynolds_mu;

  VarFcn *varFcn;
  ViscoFcn *viscoFcn;
  ThermalCondFcn *thermalCondFcn;

  void computeTemperature(double *, double &);
  void computeVelocity(double *[4], double [4][3], double [3]);
  void computeTemperature(double *[4], double [4], double &);
  void computeVelocityGradient(double [4][3], double [4][3], double [3][3]);
  void computeTemperatureGradient(double [4][3], double [4], double [3]);
  void computeStressTensor(double, double, double [3][3], double [3][3]);
  void computeHeatFluxVector(double, double [3], double [3]);

  template<int dim>
  void computeVolumeTermNS(double, double, double, double [3], double [3][3], 
			   double [3], double (*)[dim]);

  template<int neq>
  void computeJacobianVolumeTermNS(double [4][3], double, double,  double, double *[4], 
				   double [4], double (*)[3][neq][neq]);

  void computeSurfaceTermNS(double [4][3], Vec3D &, double *, double *[4], double *);

  template<int neq>
  void computeJacobianSurfaceTermNS(double [4][3], Vec3D &, double *, 
				    double *[4], double (*)[neq][neq]);
  
// Included (MB*)
  double reynoldsNS;
  double reynolds_muNS;
  double dRedMachNS;
  double dRe_mudMachNS;
  void computeDerivativeOfVelocity(double *[4], double [4][3], double [3]);
  void computeDerivativeOperatorsOfVelocity(double [4][3][4][4], double [3][4][4]);

  void computeDerivativeOfTemperature(double *[4], double *[4], double [4], double &);
  void computeDerivativeOperatorsOfTemperature(double *[4], double [4][5], double [4]);
  void computeDerivativeOperatorsOfTemperature2(double *V[4], double dTdV[4][4][5], double dTcgdV[4][5]);

  void computeDerivativeOfVelocityGradient(double [4][3], double [4][3], double [4][3], double [4][3], double [3][3]);
  void computeDerivativeOperatorsOfVelocityGradient(double [4][3], double [4][3], double [3][3][4][3], double [3][3][4][3]);

  void computeDerivativeOfTemperatureGradient(double [4][3], double [4][3], double [4], double [4], double [3]);
  void computeDerivativeOperatorsOfTemperatureGradient(double [4][3], double [4], double [3][4][3], double [3][4]);

  void computeDerivativeOfStressTensor(double, double, double, double, double [3][3], double [3][3], double [3][3]);
  void computeDerivativeOperatorsOfStressTensor(double, double, double[3][3], double [3][3][3][3], double [3][3], double [3][3]);

  void computeDerivativeOfHeatFluxVector(double, double, double [3], double [3], double [3]);
  void computeDerivativeOperatorsOfHeatFluxVector(double, double [3], double [3][3], double [3]);

  template<int dim>
  void computeDerivativeOfVolumeTermNS(double, double, double, double, double, double, double [3], double [3], double [3][3], double [3][3],
			   double [3], double [3], double (*)[dim]);
  template<int dim>
  void computeDerivativeOperatorsOfVolumeTermNS(double, double, double, double [3],
         double [3][3], double [3], double [3][dim], double [3][dim], double [3][dim],
         double [3][dim][3], double [3][dim][3][3], double [3][dim][3]);

  void computeDerivativeOfSurfaceTermNS(double [4][3], double [4][3], Vec3D &, Vec3D &, double *, double *, double *[4], double *[4], double, double *);
  void rstVarNS(IoData &, Communicator*);
  void rstVar(IoData &ioData, Communicator* com);
  template<int neq>
  void computeJacobianVolumeTermNS(double [4][3], double, double [4][5], double, double [4][5], double, double [4][5], double *[4], 
				   double [4], double (*)[3][neq][neq]);
  template<int neq>
  void computeJacobianVolumeTermNS(double [4][3], double, double [4][6], double, double [4][6], double, double [4][6], double *[4], 
				   double [4], double (*)[3][neq][neq]);
  void computeDerivativeOfTemperature(double *, double *, double &);

public:

  NavierStokesTerm(IoData &, VarFcn *);
  virtual ~NavierStokesTerm();
  
  ViscoFcn * getViscoFcn() const { return viscoFcn; }
  double get_ooreynolds_mu() const { return ooreynolds_mu; }

  ThermalCondFcn* getThermalCondFcn() const { return thermalCondFcn; }
};

//------------------------------------------------------------------------------

inline
NavierStokesTerm::NavierStokesTerm(IoData &iod, VarFcn *vf) : varFcn(vf)
{

// Included (MB)
  reynoldsNS = iod.ref.reynolds_mu;
  reynolds_muNS = iod.ref.reynolds_mu;
  dRedMachNS = iod.ref.dRe_mudMach;
  dRe_mudMachNS = iod.ref.dRe_mudMach;

  ooreynolds = 1.0/iod.ref.reynolds_mu;
  ooreynolds_mu = 1.0/iod.ref.reynolds_mu;

  viscoFcn = 0;
  thermalCondFcn = 0;

  if (iod.eqs.viscosityModel.type == ViscosityModelData::CONSTANT)
    viscoFcn = new ConstantViscoFcn(iod);
  else if (iod.eqs.viscosityModel.type == ViscosityModelData::SUTHERLAND)
    viscoFcn = new SutherlandViscoFcn(iod);
  else if (iod.eqs.viscosityModel.type == ViscosityModelData::PRANDTL)
    viscoFcn = new PrandtlViscoFcn(iod);

  if (iod.eqs.thermalCondModel.type == ThermalCondModelData::CONSTANT_PRANDTL)
    thermalCondFcn = new ConstantPrandtlThermalCondFcn(iod, viscoFcn, varFcn);
  else if (iod.eqs.thermalCondModel.type == ThermalCondModelData::CONSTANT)
    thermalCondFcn = new ConstantThermalCondFcn(iod);

}

//------------------------------------------------------------------------------

inline
NavierStokesTerm::~NavierStokesTerm()
{

  if (viscoFcn) delete viscoFcn;
  if (thermalCondFcn) delete thermalCondFcn;

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void NavierStokesTerm::rstVarNS(IoData &iod, Communicator *com)
{

  reynoldsNS = iod.ref.reynolds_mu;
  reynolds_muNS = iod.ref.reynolds_mu;
  dRedMachNS = iod.ref.dRe_mudMach;
  dRe_mudMachNS = iod.ref.dRe_mudMach;
  ooreynolds = 1.0 / iod.ref.reynolds_mu;
  ooreynolds_mu = 1.0 / iod.ref.reynolds_mu;
  viscoFcn->rstVar(iod);
  thermalCondFcn->rstVar(iod);

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void NavierStokesTerm::rstVar(IoData &iod, Communicator *com)
{

  reynoldsNS = iod.ref.reynolds_mu;
  reynolds_muNS = iod.ref.reynolds_mu;
  dRedMachNS = iod.ref.dRe_mudMach;
  dRe_mudMachNS = iod.ref.dRe_mudMach;
  ooreynolds = 1.0 / iod.ref.reynolds_mu;
  ooreynolds_mu = 1.0 / iod.ref.reynolds_mu;
  viscoFcn->rstVar(iod);
  thermalCondFcn->rstVar(iod);

}

//------------------------------------------------------------------------------

inline
void NavierStokesTerm::computeTemperature(double *V, double &T)
{
  T = varFcn->computeTemperature(V);
}

//------------------------------------------------------------------------------

// Included (MB)
inline
void NavierStokesTerm::computeDerivativeOfTemperature(double *V, double *dV, double &dT)
{
  dT = varFcn->computeDerivativeOfTemperature(V,dV);
}

//------------------------------------------------------------------------------

inline
void NavierStokesTerm::computeVelocity(double *V[4], double u[4][3], double ucg[3])
{

  u[0][0] = V[0][1];
  u[0][1] = V[0][2];
  u[0][2] = V[0][3];

  u[1][0] = V[1][1];
  u[1][1] = V[1][2];
  u[1][2] = V[1][3];

  u[2][0] = V[2][1];
  u[2][1] = V[2][2];
  u[2][2] = V[2][3];

  u[3][0] = V[3][1];
  u[3][1] = V[3][2];
  u[3][2] = V[3][3];

  ucg[0] = fourth * (u[0][0] + u[1][0] + u[2][0] + u[3][0]);
  ucg[1] = fourth * (u[0][1] + u[1][1] + u[2][1] + u[3][1]);
  ucg[2] = fourth * (u[0][2] + u[1][2] + u[2][2] + u[3][2]);

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void NavierStokesTerm::computeDerivativeOfVelocity(double *dV[4], double du[4][3], double ducg[3])
{

  du[0][0] = dV[0][1];
  du[0][1] = dV[0][2];
  du[0][2] = dV[0][3];

  du[1][0] = dV[1][1];
  du[1][1] = dV[1][2];
  du[1][2] = dV[1][3];

  du[2][0] = dV[2][1];
  du[2][1] = dV[2][2];
  du[2][2] = dV[2][3];

  du[3][0] = dV[3][1];
  du[3][1] = dV[3][2];
  du[3][2] = dV[3][3];

  ducg[0] = fourth * (du[0][0] + du[1][0] + du[2][0] + du[3][0]);
  ducg[1] = fourth * (du[0][1] + du[1][1] + du[2][1] + du[3][1]);
  ducg[2] = fourth * (du[0][2] + du[1][2] + du[2][2] + du[3][2]);

}

// Included (YC)
inline
void NavierStokesTerm::computeDerivativeOperatorsOfVelocity(double dudV[4][3][4][4], double ducgdV[3][4][4])
{

  dudV[0][2][0][3] = dudV[0][1][0][2] = dudV[0][0][0][1] = 1.0;
  dudV[1][2][1][3] = dudV[1][1][1][2] = dudV[1][0][1][1] = 1.0;
  dudV[2][2][2][3] = dudV[2][1][2][2] = dudV[2][0][2][1] = 1.0;
  dudV[3][2][3][3] = dudV[3][1][3][2] = dudV[3][0][3][1] = 1.0;
  double ducgdu[3][4][3] = {0};
  ducgdu[0][0][0] = ducgdu[0][1][0] = ducgdu[0][2][0] = ducgdu[0][3][0] = fourth;
  ducgdu[1][0][1] = ducgdu[1][1][1] = ducgdu[1][2][1] = ducgdu[1][3][1] = fourth;
  ducgdu[2][0][2] = ducgdu[2][1][2] = ducgdu[2][2][2] = ducgdu[2][3][2] = fourth;
  for(int i=0; i<3; ++i)
    for(int j=0; j<4; ++j)
      for(int k=0; k<3; ++k)
        for(int l=0; l<4; ++l)
          for(int m=0; m<4; ++m)
            ducgdV[i][l][m] += ducgdu[i][j][k]*dudV[j][k][l][m];
}

//------------------------------------------------------------------------------

inline
void NavierStokesTerm::computeTemperature(double *V[4], double T[4], double &Tcg)
{

  T[0] = varFcn->computeTemperature(V[0]);
  T[1] = varFcn->computeTemperature(V[1]);
  T[2] = varFcn->computeTemperature(V[2]);
  T[3] = varFcn->computeTemperature(V[3]);

  Tcg = fourth * (T[0] + T[1] + T[2] + T[3]);

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void NavierStokesTerm::computeDerivativeOfTemperature(double *V[4], double *dV[4], double dT[4], double &dTcg)
{

  dT[0] = varFcn->computeDerivativeOfTemperature(V[0], dV[0]);
  dT[1] = varFcn->computeDerivativeOfTemperature(V[1], dV[1]);
  dT[2] = varFcn->computeDerivativeOfTemperature(V[2], dV[2]);
  dT[3] = varFcn->computeDerivativeOfTemperature(V[3], dV[3]);

  dTcg = fourth * (dT[0] + dT[1] + dT[2] + dT[3]);

}

//------------------------------------------------------------------------------

// Included (YC)
inline
void NavierStokesTerm::computeDerivativeOperatorsOfTemperature2(double *V[4], double dTdV[4][4][5], double dTcgdV[4][5])
{
  double dTdV0[5] = {0};
  varFcn->computeDerivativeOperatorsOfTemperature(V[0], dTdV0);
  dTdV[0][0][0] = dTdV0[0];  dTdV[0][0][4] = dTdV0[4];
  double dTdV1[5] = {0};
  varFcn->computeDerivativeOperatorsOfTemperature(V[1], dTdV1);
  dTdV[1][1][0] = dTdV1[0];  dTdV[1][1][4] = dTdV1[4];
  double dTdV2[5] = {0};
  varFcn->computeDerivativeOperatorsOfTemperature(V[2], dTdV2);
  dTdV[2][2][0] = dTdV2[0];  dTdV[2][2][4] = dTdV2[4];
  double dTdV3[5] = {0};
  varFcn->computeDerivativeOperatorsOfTemperature(V[3], dTdV3);
  dTdV[3][3][0] = dTdV3[0];  dTdV[3][3][4] = dTdV3[4];

  for(int k=0; k<4; ++k)
    for(int j=0; j<4; ++j)
      for(int i=0; i<5; ++i)
        dTcgdV[j][i] += dTdV[k][j][i]*fourth;
}

//------------------------------------------------------------------------------

// Included (YC)
inline
void NavierStokesTerm::computeDerivativeOperatorsOfTemperature(double *V[4], double dTdV[4][5], double dTcgdV[4])
{

  varFcn->computeDerivativeOperatorsOfTemperature(V[0], dTdV[0]);
  varFcn->computeDerivativeOperatorsOfTemperature(V[1], dTdV[1]);
  varFcn->computeDerivativeOperatorsOfTemperature(V[2], dTdV[2]);
  varFcn->computeDerivativeOperatorsOfTemperature(V[3], dTdV[3]);

  for(int j=0; j<4; ++j)
    for(int i=0; i<5; ++i)
    dTcgdV[j] = dTdV[j][i]*fourth;
}

//------------------------------------------------------------------------------


inline
void NavierStokesTerm::computeVelocityGradient(double dp1dxj[4][3], double u[4][3],
					       double dudxj[3][3])
{

  dudxj[0][0] = dp1dxj[0][0]*u[0][0] + dp1dxj[1][0]*u[1][0] + 
    dp1dxj[2][0]*u[2][0] + dp1dxj[3][0]*u[3][0];

  dudxj[0][1] = dp1dxj[0][1]*u[0][0] + dp1dxj[1][1]*u[1][0] + 
    dp1dxj[2][1]*u[2][0] + dp1dxj[3][1]*u[3][0];

  dudxj[0][2] = dp1dxj[0][2]*u[0][0] + dp1dxj[1][2]*u[1][0] + 
    dp1dxj[2][2]*u[2][0] + dp1dxj[3][2]*u[3][0];

  dudxj[1][0] = dp1dxj[0][0]*u[0][1] + dp1dxj[1][0]*u[1][1] + 
    dp1dxj[2][0]*u[2][1] + dp1dxj[3][0]*u[3][1];

  dudxj[1][1] = dp1dxj[0][1]*u[0][1] + dp1dxj[1][1]*u[1][1] + 
    dp1dxj[2][1]*u[2][1] + dp1dxj[3][1]*u[3][1];

  dudxj[1][2] = dp1dxj[0][2]*u[0][1] + dp1dxj[1][2]*u[1][1] + 
    dp1dxj[2][2]*u[2][1] + dp1dxj[3][2]*u[3][1];

  dudxj[2][0] = dp1dxj[0][0]*u[0][2] + dp1dxj[1][0]*u[1][2] + 
    dp1dxj[2][0]*u[2][2] + dp1dxj[3][0]*u[3][2];

  dudxj[2][1] = dp1dxj[0][1]*u[0][2] + dp1dxj[1][1]*u[1][2] + 
    dp1dxj[2][1]*u[2][2] + dp1dxj[3][1]*u[3][2];

  dudxj[2][2] = dp1dxj[0][2]*u[0][2] + dp1dxj[1][2]*u[1][2] + 
    dp1dxj[2][2]*u[2][2] + dp1dxj[3][2]*u[3][2];

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void NavierStokesTerm::computeDerivativeOfVelocityGradient(double dp1dxj[4][3], double ddp1dxj[4][3], double u[4][3], double du[4][3], double ddudxj[3][3])
{

  ddudxj[0][0] = ddp1dxj[0][0]*u[0][0] + ddp1dxj[1][0]*u[1][0] +
                 ddp1dxj[2][0]*u[2][0] + ddp1dxj[3][0]*u[3][0] +
                 dp1dxj[0][0]*du[0][0] + dp1dxj[1][0]*du[1][0] +
                 dp1dxj[2][0]*du[2][0] + dp1dxj[3][0]*du[3][0];

  ddudxj[0][1] = ddp1dxj[0][1]*u[0][0] + ddp1dxj[1][1]*u[1][0] +
                 ddp1dxj[2][1]*u[2][0] + ddp1dxj[3][1]*u[3][0] +
                 dp1dxj[0][1]*du[0][0] + dp1dxj[1][1]*du[1][0] +
                 dp1dxj[2][1]*du[2][0] + dp1dxj[3][1]*du[3][0];

  ddudxj[0][2] = ddp1dxj[0][2]*u[0][0] + ddp1dxj[1][2]*u[1][0] +
                 ddp1dxj[2][2]*u[2][0] + ddp1dxj[3][2]*u[3][0] +
                 dp1dxj[0][2]*du[0][0] + dp1dxj[1][2]*du[1][0] +
                 dp1dxj[2][2]*du[2][0] + dp1dxj[3][2]*du[3][0];

  ddudxj[1][0] = ddp1dxj[0][0]*u[0][1] + ddp1dxj[1][0]*u[1][1] +
                 ddp1dxj[2][0]*u[2][1] + ddp1dxj[3][0]*u[3][1] +
                 dp1dxj[0][0]*du[0][1] + dp1dxj[1][0]*du[1][1] +
                 dp1dxj[2][0]*du[2][1] + dp1dxj[3][0]*du[3][1];

  ddudxj[1][1] = ddp1dxj[0][1]*u[0][1] + ddp1dxj[1][1]*u[1][1] +
                 ddp1dxj[2][1]*u[2][1] + ddp1dxj[3][1]*u[3][1] +
                 dp1dxj[0][1]*du[0][1] + dp1dxj[1][1]*du[1][1] +
                 dp1dxj[2][1]*du[2][1] + dp1dxj[3][1]*du[3][1];

  ddudxj[1][2] = ddp1dxj[0][2]*u[0][1] + ddp1dxj[1][2]*u[1][1] +
                 ddp1dxj[2][2]*u[2][1] + ddp1dxj[3][2]*u[3][1] +
                 dp1dxj[0][2]*du[0][1] + dp1dxj[1][2]*du[1][1] +
                 dp1dxj[2][2]*du[2][1] + dp1dxj[3][2]*du[3][1];

  ddudxj[2][0] = ddp1dxj[0][0]*u[0][2] + ddp1dxj[1][0]*u[1][2] +
                 ddp1dxj[2][0]*u[2][2] + ddp1dxj[3][0]*u[3][2] +
                 dp1dxj[0][0]*du[0][2] + dp1dxj[1][0]*du[1][2] +
                 dp1dxj[2][0]*du[2][2] + dp1dxj[3][0]*du[3][2];

  ddudxj[2][1] = ddp1dxj[0][1]*u[0][2] + ddp1dxj[1][1]*u[1][2] +
                 ddp1dxj[2][1]*u[2][2] + ddp1dxj[3][1]*u[3][2] +
                 dp1dxj[0][1]*du[0][2] + dp1dxj[1][1]*du[1][2] +
                 dp1dxj[2][1]*du[2][2] + dp1dxj[3][1]*du[3][2];

  ddudxj[2][2] = ddp1dxj[0][2]*u[0][2] + ddp1dxj[1][2]*u[1][2] +
                 ddp1dxj[2][2]*u[2][2] + ddp1dxj[3][2]*u[3][2] +
                 dp1dxj[0][2]*du[0][2] + dp1dxj[1][2]*du[1][2] +
                 dp1dxj[2][2]*du[2][2] + dp1dxj[3][2]*du[3][2];

}

//------------------------------------------------------------------------------



// Included (YC)
inline
void NavierStokesTerm::computeDerivativeOperatorsOfVelocityGradient(double dp1dxj[4][3], double u[4][3],
                                                                    double ddudxjddp1dxj[3][3][4][3], double ddudxjdu[3][3][4][3])
{

  if(ddudxjddp1dxj) {
    ddudxjddp1dxj[0][0][0][0] = u[0][0];  ddudxjddp1dxj[0][0][1][0] = u[1][0];  ddudxjddp1dxj[0][0][2][0] = u[2][0];  ddudxjddp1dxj[0][0][3][0] = u[3][0];
    ddudxjddp1dxj[0][1][0][1] = u[0][0];  ddudxjddp1dxj[0][1][1][1] = u[1][0];  ddudxjddp1dxj[0][1][2][1] = u[2][0];  ddudxjddp1dxj[0][1][3][1] = u[3][0];
    ddudxjddp1dxj[0][2][0][2] = u[0][0];  ddudxjddp1dxj[0][2][1][2] = u[1][0];  ddudxjddp1dxj[0][2][2][2] = u[2][0];  ddudxjddp1dxj[0][2][3][2] = u[3][0];
    ddudxjddp1dxj[1][0][0][0] = u[0][1];  ddudxjddp1dxj[1][0][1][0] = u[1][1];  ddudxjddp1dxj[1][0][2][0] = u[2][1];  ddudxjddp1dxj[1][0][3][0] = u[3][1];
    ddudxjddp1dxj[1][1][0][1] = u[0][1];  ddudxjddp1dxj[1][1][1][1] = u[1][1];  ddudxjddp1dxj[1][1][2][1] = u[2][1];  ddudxjddp1dxj[1][1][3][1] = u[3][1];
    ddudxjddp1dxj[1][2][0][2] = u[0][1];  ddudxjddp1dxj[1][2][1][2] = u[1][1];  ddudxjddp1dxj[1][2][2][2] = u[2][1];  ddudxjddp1dxj[1][2][3][2] = u[3][1];
    ddudxjddp1dxj[2][0][0][0] = u[0][2];  ddudxjddp1dxj[2][0][1][0] = u[1][2];  ddudxjddp1dxj[2][0][2][0] = u[2][2];  ddudxjddp1dxj[2][0][3][0] = u[3][2];
    ddudxjddp1dxj[2][1][0][1] = u[0][2];  ddudxjddp1dxj[2][1][1][1] = u[1][2];  ddudxjddp1dxj[2][1][2][1] = u[2][2];  ddudxjddp1dxj[2][1][3][1] = u[3][2];
    ddudxjddp1dxj[2][2][0][2] = u[0][2];  ddudxjddp1dxj[2][2][1][2] = u[1][2];  ddudxjddp1dxj[2][2][2][2] = u[2][2];  ddudxjddp1dxj[2][2][3][2] = u[3][2];
  }

  if(ddudxjdu) {
    ddudxjdu[0][0][0][0] = dp1dxj[0][0];  ddudxjdu[0][0][1][0] = dp1dxj[1][0];  ddudxjdu[0][0][2][0] = dp1dxj[2][0];  ddudxjdu[0][0][3][0] = dp1dxj[3][0];
    ddudxjdu[0][1][0][0] = dp1dxj[0][1];  ddudxjdu[0][1][1][0] = dp1dxj[1][1];  ddudxjdu[0][1][2][0] = dp1dxj[2][1];  ddudxjdu[0][1][3][0] = dp1dxj[3][1];
    ddudxjdu[0][2][0][0] = dp1dxj[0][2];  ddudxjdu[0][2][1][0] = dp1dxj[1][2];  ddudxjdu[0][2][2][0] = dp1dxj[2][2];  ddudxjdu[0][2][3][0] = dp1dxj[3][2];
    ddudxjdu[1][0][0][1] = dp1dxj[0][0];  ddudxjdu[1][0][1][1] = dp1dxj[1][0];  ddudxjdu[1][0][2][1] = dp1dxj[2][0];  ddudxjdu[1][0][3][1] = dp1dxj[3][0];
    ddudxjdu[1][1][0][1] = dp1dxj[0][1];  ddudxjdu[1][1][1][1] = dp1dxj[1][1];  ddudxjdu[1][1][2][1] = dp1dxj[2][1];  ddudxjdu[1][1][3][1] = dp1dxj[3][1];
    ddudxjdu[1][2][0][1] = dp1dxj[0][2];  ddudxjdu[1][2][1][1] = dp1dxj[1][2];  ddudxjdu[1][2][2][1] = dp1dxj[2][2];  ddudxjdu[1][2][3][1] = dp1dxj[3][2];
    ddudxjdu[2][0][0][2] = dp1dxj[0][0];  ddudxjdu[2][0][1][2] = dp1dxj[1][0];  ddudxjdu[2][0][2][2] = dp1dxj[2][0];  ddudxjdu[2][0][3][2] = dp1dxj[3][0];
    ddudxjdu[2][1][0][2] = dp1dxj[0][1];  ddudxjdu[2][1][1][2] = dp1dxj[1][1];  ddudxjdu[2][1][2][2] = dp1dxj[2][1];  ddudxjdu[2][1][3][2] = dp1dxj[3][1];
    ddudxjdu[2][2][0][2] = dp1dxj[0][2];  ddudxjdu[2][2][1][2] = dp1dxj[1][2];  ddudxjdu[2][2][2][2] = dp1dxj[2][2];  ddudxjdu[2][2][3][2] = dp1dxj[3][2];
  }

}

//------------------------------------------------------------------------------



inline
void NavierStokesTerm::computeTemperatureGradient(double dp1dxj[4][3], double T[4],
						  double dTdxj[3])
{

  dTdxj[0] = dp1dxj[0][0]*T[0] + dp1dxj[1][0]*T[1] + 
    dp1dxj[2][0]*T[2] + dp1dxj[3][0]*T[3];

  dTdxj[1] = dp1dxj[0][1]*T[0] + dp1dxj[1][1]*T[1] + 
    dp1dxj[2][1]*T[2] + dp1dxj[3][1]*T[3];

  dTdxj[2] = dp1dxj[0][2]*T[0] + dp1dxj[1][2]*T[1] + 
    dp1dxj[2][2]*T[2] + dp1dxj[3][2]*T[3];

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void NavierStokesTerm::computeDerivativeOfTemperatureGradient(double dp1dxj[4][3], double ddp1dxj[4][3], double T[4], double dT[4],
						  double ddTdxj[3])
{

  ddTdxj[0] = ddp1dxj[0][0]*T[0] + dp1dxj[0][0]*dT[0] + ddp1dxj[1][0]*T[1] + dp1dxj[1][0]*dT[1] +
    ddp1dxj[2][0]*T[2] + dp1dxj[2][0]*dT[2] + ddp1dxj[3][0]*T[3] + dp1dxj[3][0]*dT[3];

  ddTdxj[1] = ddp1dxj[0][1]*T[0] + dp1dxj[0][1]*dT[0] + ddp1dxj[1][1]*T[1] +  dp1dxj[1][1]*dT[1] +
    ddp1dxj[2][1]*T[2] + dp1dxj[2][1]*dT[2] + ddp1dxj[3][1]*T[3] + dp1dxj[3][1]*dT[3];

  ddTdxj[2] = ddp1dxj[0][2]*T[0] + dp1dxj[0][2]*dT[0] + ddp1dxj[1][2]*T[1] + dp1dxj[1][2]*dT[1] +
    ddp1dxj[2][2]*T[2] + dp1dxj[2][2]*dT[2] + ddp1dxj[3][2]*T[3]  + dp1dxj[3][2]*dT[3];

}

//------------------------------------------------------------------------------

// Included (YC)
inline
void NavierStokesTerm::computeDerivativeOperatorsOfTemperatureGradient(double dp1dxj[4][3], double T[4],
						  double ddTdxjddp1dxj[3][4][3], double ddTdxjdT[3][4])
{

  ddTdxjddp1dxj[0][0][0] = T[0];
  ddTdxjddp1dxj[0][1][0] = T[1];
  ddTdxjddp1dxj[0][2][0] = T[2];
  ddTdxjddp1dxj[0][3][0] = T[3];
  ddTdxjdT[0][0] = dp1dxj[0][0];
  ddTdxjdT[0][1] = dp1dxj[1][0];
  ddTdxjdT[0][2] = dp1dxj[2][0];
  ddTdxjdT[0][3] = dp1dxj[3][0];

  ddTdxjddp1dxj[1][0][1] = T[0];
  ddTdxjddp1dxj[1][1][1] = T[1];
  ddTdxjddp1dxj[1][2][1] = T[2];
  ddTdxjddp1dxj[1][3][1] = T[3];
  ddTdxjdT[1][0] = dp1dxj[0][1];
  ddTdxjdT[1][1] = dp1dxj[1][1];
  ddTdxjdT[1][2] = dp1dxj[2][1];
  ddTdxjdT[1][3] = dp1dxj[3][1];

  ddTdxjddp1dxj[2][0][2] = T[0];
  ddTdxjddp1dxj[2][1][2] = T[1];
  ddTdxjddp1dxj[2][2][2] = T[2];
  ddTdxjddp1dxj[2][3][2] = T[3];
  ddTdxjdT[2][0] = dp1dxj[0][2];
  ddTdxjdT[2][1] = dp1dxj[1][2];
  ddTdxjdT[2][2] = dp1dxj[2][2];
  ddTdxjdT[2][3] = dp1dxj[3][2];

}

//------------------------------------------------------------------------------



inline
void NavierStokesTerm::computeStressTensor(double mu, double lambda, double dudxj[3][3], double tij[3][3])
{

  double div = dudxj[0][0] + dudxj[1][1] + dudxj[2][2]; 

  tij[0][0] = lambda * div + 2.0 * mu *dudxj[0][0];
  tij[1][1] = lambda * div + 2.0 * mu *dudxj[1][1];
  tij[2][2] = lambda * div + 2.0 * mu *dudxj[2][2];
  tij[0][1] = mu * (dudxj[1][0] + dudxj[0][1]);
  tij[0][2] = mu * (dudxj[2][0] + dudxj[0][2]);
  tij[1][2] = mu * (dudxj[2][1] + dudxj[1][2]);
  tij[1][0] = tij[0][1];
  tij[2][0] = tij[0][2];
  tij[2][1] = tij[1][2];

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void NavierStokesTerm::computeDerivativeOfStressTensor(double mu, double dmu, double lambda, double dlambda, double dudxj[3][3], double ddudxj[3][3], double dtij[3][3])
{

  double div = dudxj[0][0] + dudxj[1][1] + dudxj[2][2]; 
  double ddiv = ddudxj[0][0] + ddudxj[1][1] + ddudxj[2][2]; 

  dtij[0][0] = dlambda * div + lambda * ddiv + 2.0 * dmu *dudxj[0][0] + 2.0 * mu *ddudxj[0][0];
  dtij[1][1] = dlambda * div + lambda * ddiv + 2.0 * dmu *dudxj[1][1] + 2.0 * mu *ddudxj[1][1];
  dtij[2][2] = dlambda * div + lambda * ddiv + 2.0 * dmu *dudxj[2][2] + 2.0 * mu *ddudxj[2][2];
  dtij[0][1] = dmu * (dudxj[1][0] + dudxj[0][1]) + mu * (ddudxj[1][0] + ddudxj[0][1]);
  dtij[0][2] = dmu * (dudxj[2][0] + dudxj[0][2]) + mu * (ddudxj[2][0] + ddudxj[0][2]);
  dtij[1][2] = dmu * (dudxj[2][1] + dudxj[1][2]) + mu * (ddudxj[2][1] + ddudxj[1][2]);
  dtij[1][0] = dtij[0][1];
  dtij[2][0] = dtij[0][2];
  dtij[2][1] = dtij[1][2];

}

//------------------------------------------------------------------------------


// Included (YC)
inline
void NavierStokesTerm::computeDerivativeOperatorsOfStressTensor(double mu, double lambda, double dudxj[3][3],
                                                                double dtijddudxj[3][3][3][3], double dtijdmu[3][3], double dtijdlambda[3][3])
{

  double div = dudxj[0][0] + dudxj[1][1] + dudxj[2][2];
  if(dtijdlambda) {
    dtijdlambda[0][0] = dtijdlambda[1][1] = dtijdlambda[2][2] = div;
  }
  if(dtijdmu) {
    dtijdmu[0][0] =  2.0*dudxj[0][0];
    dtijdmu[1][1] =  2.0*dudxj[1][1];
    dtijdmu[2][2] =  2.0*dudxj[2][2];
    dtijdmu[1][0] = dtijdmu[0][1] = (dudxj[1][0] + dudxj[0][1]);
    dtijdmu[2][0] = dtijdmu[0][2] = (dudxj[2][0] + dudxj[0][2]);
    dtijdmu[2][1] = dtijdmu[1][2] = (dudxj[2][1] + dudxj[1][2]);
  }
  if(dtijddudxj) {
    dtijddudxj[0][0][0][0] = dtijddudxj[1][1][0][0] = dtijddudxj[2][2][0][0] = lambda;
    dtijddudxj[0][0][1][1] = dtijddudxj[1][1][1][1] = dtijddudxj[2][2][1][1] = lambda;
    dtijddudxj[0][0][2][2] = dtijddudxj[1][1][2][2] = dtijddudxj[2][2][2][2] = lambda;
    dtijddudxj[0][0][0][0] += 2.0 * mu;
    dtijddudxj[1][1][1][1] += 2.0 * mu;
    dtijddudxj[2][2][2][2] += 2.0 * mu;
    dtijddudxj[1][0][0][1] = dtijddudxj[1][0][1][0] = dtijddudxj[0][1][0][1] = dtijddudxj[0][1][1][0] = mu;
    dtijddudxj[2][0][0][2] = dtijddudxj[2][0][2][0] = dtijddudxj[0][2][0][2] = dtijddudxj[0][2][2][0] = mu;
    dtijddudxj[2][1][1][2] = dtijddudxj[2][1][2][1] = dtijddudxj[1][2][1][2] = dtijddudxj[1][2][2][1] = mu;
  }

}

//------------------------------------------------------------------------------



inline
void NavierStokesTerm::computeHeatFluxVector(double kappa, double dTdxj[3], double qj[3])
{

  qj[0] = - kappa * dTdxj[0];
  qj[1] = - kappa * dTdxj[1];
  qj[2] = - kappa * dTdxj[2];

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void NavierStokesTerm::computeDerivativeOfHeatFluxVector(double kappa, double dkappa, double dTdxj[3], double ddTdxj[3], double dqj[3])
{

  dqj[0] = - dkappa * dTdxj[0] - kappa * ddTdxj[0];
  dqj[1] = - dkappa * dTdxj[1] - kappa * ddTdxj[1];
  dqj[2] = - dkappa * dTdxj[2] - kappa * ddTdxj[2];

}

//------------------------------------------------------------------------------


// Included (YC)
inline
void NavierStokesTerm::computeDerivativeOperatorsOfHeatFluxVector(double kappa, double dTdxj[3], double dqjddTdxj[3][3], double dqjdkappa[3])
{

  dqjddTdxj[2][2] = dqjddTdxj[1][1] = dqjddTdxj[0][0] = -kappa;
  dqjdkappa[0] = -dTdxj[0];
  dqjdkappa[1] = -dTdxj[1];
  dqjdkappa[2] = -dTdxj[2];

}

//------------------------------------------------------------------------------



template<int dim>
inline
void NavierStokesTerm::computeVolumeTermNS(double mu, double lambda, double kappa, double u[3], 
					   double dudxj[3][3], double dTdxj[3], 
					   double (*r)[dim])
{

  double tij[3][3];
  computeStressTensor(mu, lambda, dudxj, tij);

  double qj[3];
  computeHeatFluxVector(kappa, dTdxj, qj);

  r[0][0] = 0.0;
  r[0][1] = tij[0][0];
  r[0][2] = tij[1][0];
  r[0][3] = tij[2][0];
  r[0][4] = u[0] * tij[0][0] + u[1] * tij[1][0] + u[2] * tij[2][0] - qj[0];

  r[1][0] = 0.0;
  r[1][1] = tij[0][1];
  r[1][2] = tij[1][1];
  r[1][3] = tij[2][1];
  r[1][4] = u[0] * tij[0][1] + u[1] * tij[1][1] + u[2] * tij[2][1] - qj[1]; 

  r[2][0] = 0.0;
  r[2][1] = tij[0][2];
  r[2][2] = tij[1][2];
  r[2][3] = tij[2][2];
  r[2][4] = u[0] * tij[0][2] + u[1] * tij[1][2] + u[2] * tij[2][2] - qj[2];

}

//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
inline
void NavierStokesTerm::computeDerivativeOperatorsOfVolumeTermNS(double mu, double lambda, double kappa, double u[3],
					   double dudxj[3][3], double dTdxj[3],
					   double drdmu[3][dim], double drdlambda[3][dim], double drdkappa[3][dim],
             double drdu[3][dim][3], double drddudxj[3][dim][3][3], double drddTdxj[3][dim][3])
{
  double tij[3][3];
  computeStressTensor(mu, lambda, dudxj, tij);

  double dtijdlambda[3][3] = {0}, dtijddudxj[3][3][3][3] = {0}, dtijdmu[3][3] = {0};
  computeDerivativeOperatorsOfStressTensor(mu, lambda, dudxj, dtijddudxj, dtijdmu, dtijdlambda);

  double dqjdkappa[3] = {0}, dqjddTdxj[3][3] = {0};
  computeDerivativeOperatorsOfHeatFluxVector(kappa, dTdxj, dqjddTdxj, dqjdkappa);

  double drdtij[3][dim][3][3] = {0};
  drdtij[0][3][2][0] = drdtij[0][2][1][0] = drdtij[0][1][0][0] = 1.0;
  drdtij[1][3][2][1] = drdtij[1][2][1][1] = drdtij[1][1][0][1] = 1.0;
  drdtij[2][3][2][2] = drdtij[2][2][1][2] = drdtij[2][1][0][2] = 1.0;
  drdtij[0][4][0][0] = drdtij[1][4][0][1] = drdtij[2][4][0][2] = u[0];
  drdtij[0][4][1][0] = drdtij[1][4][1][1] = drdtij[2][4][1][2] = u[1];
  drdtij[0][4][2][0] = drdtij[1][4][2][1] = drdtij[2][4][2][2] = u[2];
  drdu[0][4][0] = tij[0][0];  drdu[0][4][1] = tij[1][0];  drdu[0][4][2] = tij[2][0];
  drdu[1][4][0] = tij[0][1];  drdu[1][4][1] = tij[1][1];  drdu[1][4][2] = tij[2][1];
  drdu[2][4][0] = tij[0][2];  drdu[2][4][1] = tij[1][2];  drdu[2][4][2] = tij[2][2];
  double drdqj[3][dim][3] = {0};
  drdqj[2][4][2] = drdqj[1][4][1] = drdqj[0][4][0] = -1.0;

  for(int i=0; i<3; ++i)
    for(int j=0; j<dim; ++j)
      for(int k=0; k<3; ++k) {
        drdkappa[i][j] += drdqj[i][j][k]*dqjdkappa[k];
        for(int l=0; l<3; ++l) {
          drddTdxj[i][j][l] += drdqj[i][j][k]*dqjddTdxj[k][l];
          drdmu[i][j] += drdtij[i][j][k][l]*dtijdmu[k][l];
          drdlambda[i][j] += drdtij[i][j][k][l]*dtijdlambda[k][l];
          for(int m=0; m<3; ++m)
            for(int n=0; n<3; ++n) {
              drddudxj[i][j][m][n] += drdtij[i][j][k][l]*dtijddudxj[k][l][m][n];
            }
        }
      }

}

//------------------------------------------------------------------------------


// Included (MB)
template<int dim>
inline
void NavierStokesTerm::computeDerivativeOfVolumeTermNS(double mu, double dmu, double lambda, double dlambda, double kappa, double dkappa, double u[3], double du[3],
					   double dudxj[3][3], double ddudxj[3][3], double dTdxj[3], double ddTdxj[3],
					   double (*dr)[dim])
{

  double tij[3][3];
  computeStressTensor(mu, lambda, dudxj, tij);

  double dtij[3][3];
  computeDerivativeOfStressTensor(mu, dmu, lambda, dlambda, dudxj, ddudxj, dtij);

  double dqj[3];
  computeDerivativeOfHeatFluxVector(kappa, dkappa, dTdxj, ddTdxj, dqj);

  dr[0][0] = 0.0;
  dr[0][1] = dtij[0][0];
  dr[0][2] = dtij[1][0];
  dr[0][3] = dtij[2][0];
  dr[0][4] = du[0] * tij[0][0] + u[0] * dtij[0][0] + du[1] * tij[1][0] + u[1] * dtij[1][0] + du[2] * tij[2][0]  + u[2] * dtij[2][0] - dqj[0];

  dr[1][0] = 0.0;
  dr[1][1] = dtij[0][1];
  dr[1][2] = dtij[1][1];
  dr[1][3] = dtij[2][1];
  dr[1][4] = du[0] * tij[0][1] + u[0] * dtij[0][1] + du[1] * tij[1][1] + u[1] * dtij[1][1] + du[2] * tij[2][1]  + u[2] * dtij[2][1] - dqj[1];

  dr[2][0] = 0.0;
  dr[2][1] = dtij[0][2];
  dr[2][2] = dtij[1][2];
  dr[2][3] = dtij[2][2];
  dr[2][4] = du[0] * tij[0][2] + u[0] * dtij[0][2] + du[1] * tij[1][2] + u[1] * dtij[1][2] + du[2] * tij[2][2]  + u[2] * dtij[2][2] - dqj[2];

}

//------------------------------------------------------------------------------

// Included (MB*)
template<int neq>
inline
void NavierStokesTerm::
computeJacobianVolumeTermNS(double dp1dxj[4][3], double mu, double dmu[4][5], double lambda, double dlambda[4][5], double kappa, double dkappa[4][5], 
			    double *V[4], double T[4], double (*dRdU)[3][neq][neq])
{

  double dTdxj[3];
  computeTemperatureGradient(dp1dxj, T, dTdxj);

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double tij[3][3];
  computeStressTensor(mu, lambda, dudxj, tij);

  double dtijdu0[3][3];
  double dtijdu1[3][3];
  double dtijdu2[3][3];
  double dtijdu3[3][3];
  double dtijdu4[3][3];

  double ducgdu0[3];
  double ducgdu1[3];
  double ducgdu2[3];
  double ducgdu3[3];
  double ducgdu4[3];

  double dTdu0[4];
  double dTdu1[4];
  double dTdu2[4];
  double dTdu3[4];
  double dTdu4[4];
  
  double ddTdxjdu0[3];
  double ddTdxjdu1[3];
  double ddTdxjdu2[3];
  double ddTdxjdu3[3];
  double ddTdxjdu4[3];

  double dqjdu0[3];
  double dqjdu1[3];
  double dqjdu2[3];
  double dqjdu3[3];
  double dqjdu4[3];

  for (int k=0; k<4; ++k) {

    ducgdu0[0] = - fourth * V[k][1] / V[k][0];
    ducgdu1[0] = fourth / V[k][0];
    ducgdu2[0] = 0.0;
    ducgdu3[0] = 0.0;
    ducgdu4[0] = 0.0;

    ducgdu0[1] = - fourth * V[k][2] / V[k][0];
    ducgdu1[1] = 0.0;               
    ducgdu2[1] = fourth / V[k][0];
    ducgdu3[1] = 0.0;
    ducgdu4[1] = 0.0;

    ducgdu0[2] = - fourth * V[k][3] / V[k][0];
    ducgdu1[2] = 0.0;
    ducgdu2[2] = 0.0;
    ducgdu3[2] = fourth / V[k][0];
    ducgdu4[2] = 0.0;

// Take a look on this piece of code and rederive the stress tensor
/*
    double nu = mu / rho;
    double lu = lambda / rho;

    double dtxxdu0 = -(2.0 * nu * dp1dx * u + lu * (dp1dx * u + dp1dy * v + dp1dz * w));
    double dtxxdu1 = (2.0 * nu + lu) * dp1dx;
    double dtxxdu2 = lu * dp1dy;
    double dtxxdu3 = lu * dp1dz;

    double dtyydu0 = -(2.0 * nu * dp1dy * v + lu * (dp1dx * u + dp1dy * v + dp1dz * w));
    double dtyydu1 = lu * dp1dx;
    double dtyydu2 = (2.0 * nu + lu) * dp1dy;
    double dtyydu3 = lu * dp1dz;

    double dtzzdu0 = -(2.0 * nu * dp1dz * w + lu * (dp1dx * u + dp1dy * v + dp1dz * w));
    double dtzzdu1 = lu * dp1dx;
    double dtzzdu2 = lu * dp1dy;
    double dtzzdu3 = (2.0 * nu + lu) * dp1dz;

    double dtxydu0 = - nu * (dp1dx * v + dp1dy * u);
    double dtxydu1 = nu * dp1dy;
    double dtxydu2 = nu * dp1dx;

    double dtxzdu0 = -nu * (dp1dx * w + dp1dz * u);
    double dtxzdu1 = nu * dp1dz;
    double dtxzdu3 = nu * dp1dx;

    double dtyzdu0 = -nu * (dp1dy * w + dp1dz * v);
    double dtyzdu2 = nu * dp1dz;
    double dtyzdu3 = nu * dp1dy;
*/
//    dtijdu0[0][0] = twothird * mu * (-2.0 * dp1dxj[k][0] * V[k][1] + dp1dxj[k][1] * V[k][2] + dp1dxj[k][2] * V[k][3]) / V[k][0] + dmu[k][0] * twothird * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);
//    dtijdu1[0][0] = 2.0 * twothird * mu * dp1dxj[k][0] / V[k][0] + dmu[k][1] * twothird * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);
//    dtijdu2[0][0] = - twothird * mu * dp1dxj[k][1] / V[k][0] + dmu[k][2] * twothird * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);
//    dtijdu3[0][0] = - twothird * mu * dp1dxj[k][2] / V[k][0] + dmu[k][3] * twothird * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);
//    dtijdu4[0][0] = dmu[k][4] * twothird * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);

//    dtijdu0[1][1] = twothird * mu * (-2.0 * dp1dxj[k][1] * V[k][2] + dp1dxj[k][0] * V[k][1] + dp1dxj[k][2] * V[k][3]) / V[k][0] + dmu[k][0] * twothird * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);
//    dtijdu1[1][1] = - twothird * mu * dp1dxj[k][0] / V[k][0] + dmu[k][1] * twothird * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);
//    dtijdu2[1][1] = 2.0 * twothird * mu * dp1dxj[k][1] / V[k][0] + dmu[k][2] * twothird * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);
//    dtijdu3[1][1] = - twothird * mu * dp1dxj[k][2] / V[k][0] + dmu[k][3] * twothird * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);
//    dtijdu4[1][1] = dmu[k][4] * twothird * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);

//    dtijdu0[2][2] = twothird * mu * (-2.0 * dp1dxj[k][2] * V[k][3] + dp1dxj[k][0] * V[k][1] + dp1dxj[k][1] * V[k][2]) / V[k][0] + dmu[k][0] * twothird * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);
//    dtijdu1[2][2] = - twothird * mu * dp1dxj[k][0] / V[k][0] + dmu[k][1] * twothird * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);
//    dtijdu2[2][2] = - twothird * mu * dp1dxj[k][1] / V[k][0] + dmu[k][2] * twothird * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);
//    dtijdu3[2][2] = 2.0 * twothird * mu * dp1dxj[k][2] / V[k][0] + dmu[k][3] * twothird * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);
//    dtijdu4[2][2] = dmu[k][4] * twothird * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);
//    dtijdu0[2][2] = 2.0 * dmu[k][0] * dudxj[2][2] - 2.0 * mu * dudxj[2][2] / V[k][0] + (-1.0)*twothird*dmu[k][0] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) - (-twothird*mu) * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) / V[k][0];

    dtijdu0[0][0] = 2.0 * dmu[k][0] * dudxj[0][0] - 2.0 * mu * dp1dxj[k][0] * V[k][1] / V[k][0] + dlambda[k][0] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) - lambda * (dp1dxj[k][0] * V[k][1] + dp1dxj[k][1] * V[k][2] + dp1dxj[k][2] * V[k][3]) / V[k][0];
    dtijdu1[0][0] = 2.0 * dmu[k][1] * dudxj[0][0] + 2.0 * mu * dp1dxj[k][0] / V[k][0] + dlambda[k][1] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][0] / V[k][0];
    dtijdu2[0][0] = 2.0 * dmu[k][2] * dudxj[0][0] + dlambda[k][2] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][1] / V[k][0];
    dtijdu3[0][0] = 2.0 * dmu[k][3] * dudxj[0][0] + dlambda[k][3] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][2] / V[k][0];
    dtijdu4[0][0] = 2.0 * dmu[k][4] * dudxj[0][0] + dlambda[k][4] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]);

    dtijdu0[1][1] = 2.0 * dmu[k][0] * dudxj[1][1] - 2.0 * mu * dp1dxj[k][1] * V[k][2] / V[k][0] + dlambda[k][0] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) - lambda * (dp1dxj[k][1] * V[k][2] + dp1dxj[k][0] * V[k][1] + dp1dxj[k][2] * V[k][3]) / V[k][0];
    dtijdu1[1][1] = 2.0 * dmu[k][1] * dudxj[1][1] + dlambda[k][1] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][0] / V[k][0];
    dtijdu2[1][1] = 2.0 * dmu[k][2] * dudxj[1][1] + 2.0 * mu * dp1dxj[k][1] / V[k][0] + dlambda[k][2] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][1] / V[k][0];
    dtijdu3[1][1] = 2.0 * dmu[k][3] * dudxj[1][1] + dlambda[k][3] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][2] / V[k][0];
    dtijdu4[1][1] = 2.0 * dmu[k][4] * dudxj[1][1] + dlambda[k][4] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]);

    dtijdu0[2][2] = 2.0 * dmu[k][0] * dudxj[2][2] - 2.0 * mu * dp1dxj[k][2] * V[k][3] / V[k][0] + dlambda[k][0] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) - lambda * (dp1dxj[k][2] * V[k][3] + dp1dxj[k][0] * V[k][1] + dp1dxj[k][1] * V[k][2]) / V[k][0];
    dtijdu1[2][2] = 2.0 * dmu[k][1] * dudxj[2][2] + dlambda[k][1] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][0] / V[k][0];
    dtijdu2[2][2] = 2.0 * dmu[k][2] * dudxj[2][2] + dlambda[k][2] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][1] / V[k][0];
    dtijdu3[2][2] = 2.0 * dmu[k][3] * dudxj[2][2] + 2.0 * mu * dp1dxj[k][2] / V[k][0] + dlambda[k][3] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][2] / V[k][0];
    dtijdu4[2][2] = 2.0 * dmu[k][4] * dudxj[2][2] + dlambda[k][4] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]);

    dtijdu0[0][1] = - mu * (dp1dxj[k][0] * V[k][2] + dp1dxj[k][1] * V[k][1]) / V[k][0] + dmu[k][0] * (dudxj[1][0] + dudxj[0][1]);
    dtijdu1[0][1] = mu * dp1dxj[k][1] / V[k][0] + dmu[k][1] * (dudxj[1][0] + dudxj[0][1]);
    dtijdu2[0][1] = mu * dp1dxj[k][0] / V[k][0] + dmu[k][2] * (dudxj[1][0] + dudxj[0][1]);
    dtijdu3[0][1] = dmu[k][3] * (dudxj[1][0] + dudxj[0][1]);
    dtijdu4[0][1] = dmu[k][4] * (dudxj[1][0] + dudxj[0][1]);

    dtijdu0[0][2] = - mu * (dp1dxj[k][0] * V[k][3] + dp1dxj[k][2] * V[k][1]) / V[k][0] + dmu[k][0] * (dudxj[2][0] + dudxj[0][2]);
    dtijdu1[0][2] = mu * dp1dxj[k][2] / V[k][0] + dmu[k][1] * (dudxj[2][0] + dudxj[0][2]);
    dtijdu2[0][2] = dmu[k][2] * (dudxj[2][0] + dudxj[0][2]);
    dtijdu3[0][2] = mu * dp1dxj[k][0] / V[k][0] + dmu[k][3] * (dudxj[2][0] + dudxj[0][2]);
    dtijdu4[0][2] = dmu[k][4] * (dudxj[2][0] + dudxj[0][2]);

    dtijdu0[1][2] = - mu * (dp1dxj[k][1] * V[k][3] + dp1dxj[k][2] * V[k][2]) / V[k][0] + dmu[k][0] * (dudxj[2][1] + dudxj[1][2]);
    dtijdu1[1][2] = dmu[k][1] * (dudxj[2][1] + dudxj[1][2]);
    dtijdu2[1][2] = mu * dp1dxj[k][2] / V[k][0] + dmu[k][2] * (dudxj[2][1] + dudxj[1][2]);
    dtijdu3[1][2] = mu * dp1dxj[k][1] / V[k][0] + dmu[k][3] * (dudxj[2][1] + dudxj[1][2]);
    dtijdu4[1][2] = dmu[k][4] * (dudxj[2][1] + dudxj[1][2]);

    dtijdu0[1][0] = dtijdu0[0][1];
    dtijdu1[1][0] = dtijdu1[0][1];
    dtijdu2[1][0] = dtijdu2[0][1];
    dtijdu3[1][0] = dtijdu3[0][1];
    dtijdu4[1][0] = dtijdu4[0][1];

    dtijdu0[2][0] = dtijdu0[0][2];
    dtijdu1[2][0] = dtijdu1[0][2];
    dtijdu2[2][0] = dtijdu2[0][2];
    dtijdu3[2][0] = dtijdu3[0][2];
    dtijdu4[2][0] = dtijdu4[0][2];

    dtijdu0[2][1] = dtijdu0[1][2];
    dtijdu1[2][1] = dtijdu1[1][2];
    dtijdu2[2][1] = dtijdu2[1][2];
    dtijdu3[2][1] = dtijdu3[1][2];
    dtijdu4[2][1] = dtijdu4[1][2];

    dTdu0[k] = - 1.0 / V[k][0] * (T[k] - 0.5 * (V[k][1]*V[k][1] + V[k][2]*V[k][2] + V[k][3]*V[k][3]));
    dTdu1[k] = - 1.0 / V[k][0] * V[k][1];
    dTdu2[k] = - 1.0 / V[k][0] * V[k][2];
    dTdu3[k] = - 1.0 / V[k][0] * V[k][3];
    dTdu4[k] = 1.0 / V[k][0];

    ddTdxjdu0[0] = dp1dxj[k][0]*dTdu0[k];
    ddTdxjdu1[0] = dp1dxj[k][0]*dTdu1[k];
    ddTdxjdu2[0] = dp1dxj[k][0]*dTdu2[k];
    ddTdxjdu3[0] = dp1dxj[k][0]*dTdu3[k];
    ddTdxjdu4[0] = dp1dxj[k][0]*dTdu4[k];

    ddTdxjdu0[1] = dp1dxj[k][1]*dTdu0[k];
    ddTdxjdu1[1] = dp1dxj[k][1]*dTdu1[k];
    ddTdxjdu2[1] = dp1dxj[k][1]*dTdu2[k];
    ddTdxjdu3[1] = dp1dxj[k][1]*dTdu3[k];
    ddTdxjdu4[1] = dp1dxj[k][1]*dTdu4[k];
    
    ddTdxjdu0[2] = dp1dxj[k][2]*dTdu0[k];
    ddTdxjdu1[2] = dp1dxj[k][2]*dTdu1[k];
    ddTdxjdu2[2] = dp1dxj[k][2]*dTdu2[k];
    ddTdxjdu3[2] = dp1dxj[k][2]*dTdu3[k];
    ddTdxjdu4[2] = dp1dxj[k][2]*dTdu4[k];
     
    dqjdu0[0] = - dkappa[k][0] * dTdxj[0] - kappa * ddTdxjdu0[0];
    dqjdu1[0] = - dkappa[k][1] * dTdxj[0] - kappa * ddTdxjdu1[0];
    dqjdu2[0] = - dkappa[k][2] * dTdxj[0] - kappa * ddTdxjdu2[0];
    dqjdu3[0] = - dkappa[k][3] * dTdxj[0] - kappa * ddTdxjdu3[0];
    dqjdu4[0] = - dkappa[k][4] * dTdxj[0] - kappa * ddTdxjdu4[0];
    
    dqjdu0[1] = - dkappa[k][0] * dTdxj[1] - kappa * ddTdxjdu0[1];
    dqjdu1[1] = - dkappa[k][1] * dTdxj[1] - kappa * ddTdxjdu1[1];
    dqjdu2[1] = - dkappa[k][2] * dTdxj[1] - kappa * ddTdxjdu2[1];
    dqjdu3[1] = - dkappa[k][3] * dTdxj[1] - kappa * ddTdxjdu3[1];
    dqjdu4[1] = - dkappa[k][4] * dTdxj[1] - kappa * ddTdxjdu4[1];

    dqjdu0[2] = - dkappa[k][0] * dTdxj[2] - kappa * ddTdxjdu0[2];
    dqjdu1[2] = - dkappa[k][1] * dTdxj[2] - kappa * ddTdxjdu1[2];
    dqjdu2[2] = - dkappa[k][2] * dTdxj[2] - kappa * ddTdxjdu2[2];
    dqjdu3[2] = - dkappa[k][3] * dTdxj[2] - kappa * ddTdxjdu3[2];
    dqjdu4[2] = - dkappa[k][4] * dTdxj[2] - kappa * ddTdxjdu4[2];

    // dRxdU

    dRdU[k][0][0][0] = 0.0;
    dRdU[k][0][0][1] = 0.0;
    dRdU[k][0][0][2] = 0.0;
    dRdU[k][0][0][3] = 0.0;
    dRdU[k][0][0][4] = 0.0;

    dRdU[k][0][1][0] = dtijdu0[0][0];
    dRdU[k][0][1][1] = dtijdu1[0][0];
    dRdU[k][0][1][2] = dtijdu2[0][0];
    dRdU[k][0][1][3] = dtijdu3[0][0];
    dRdU[k][0][1][4] = dtijdu4[0][0];

    dRdU[k][0][2][0] = dtijdu0[0][1];
    dRdU[k][0][2][1] = dtijdu1[0][1];
    dRdU[k][0][2][2] = dtijdu2[0][1];
    dRdU[k][0][2][3] = dtijdu3[0][1];
    dRdU[k][0][2][4] = dtijdu4[0][1];

    dRdU[k][0][3][0] = dtijdu0[0][2];
    dRdU[k][0][3][1] = dtijdu1[0][2];
    dRdU[k][0][3][2] = dtijdu2[0][2];
    dRdU[k][0][3][3] = dtijdu3[0][2];
    dRdU[k][0][3][4] = dtijdu4[0][2];

    dRdU[k][0][4][0] = ducgdu0[0] * tij[0][0] + ucg[0] * dtijdu0[0][0] + ducgdu0[1] * tij[0][1] + ucg[1] * dtijdu0[0][1] + ducgdu0[2] * tij[0][2] + ucg[2] * dtijdu0[0][2] - dqjdu0[0];
    dRdU[k][0][4][1] = ducgdu1[0] * tij[0][0] + ucg[0] * dtijdu1[0][0] + ducgdu1[1] * tij[0][1] + ucg[1] * dtijdu1[0][1] + ducgdu1[2] * tij[0][2] + ucg[2] * dtijdu1[0][2] - dqjdu1[0];
    dRdU[k][0][4][2] = ducgdu2[0] * tij[0][0] + ucg[0] * dtijdu2[0][0] + ducgdu2[1] * tij[0][1] + ucg[1] * dtijdu2[0][1] + ducgdu2[2] * tij[0][2] + ucg[2] * dtijdu2[0][2] - dqjdu2[0];
    dRdU[k][0][4][3] = ducgdu3[0] * tij[0][0] + ucg[0] * dtijdu3[0][0] + ducgdu3[1] * tij[0][1] + ucg[1] * dtijdu3[0][1] + ducgdu3[2] * tij[0][2] + ucg[2] * dtijdu3[0][2] - dqjdu3[0];
    dRdU[k][0][4][4] = ducgdu4[0] * tij[0][0] + ucg[0] * dtijdu4[0][0] + ducgdu4[1] * tij[0][1] + ucg[1] * dtijdu4[0][1] + ducgdu4[2] * tij[0][2] + ucg[2] * dtijdu4[0][2] - dqjdu4[0];

    // dRydU

    dRdU[k][1][0][0] = 0.0;
    dRdU[k][1][0][1] = 0.0;
    dRdU[k][1][0][2] = 0.0;
    dRdU[k][1][0][3] = 0.0;
    dRdU[k][1][0][4] = 0.0;
  
    dRdU[k][1][1][0] = dtijdu0[0][1];
    dRdU[k][1][1][1] = dtijdu1[0][1];
    dRdU[k][1][1][2] = dtijdu2[0][1];
    dRdU[k][1][1][3] = dtijdu3[0][1];
    dRdU[k][1][1][4] = dtijdu4[0][1];

    dRdU[k][1][2][0] = dtijdu0[1][1];
    dRdU[k][1][2][1] = dtijdu1[1][1];
    dRdU[k][1][2][2] = dtijdu2[1][1];
    dRdU[k][1][2][3] = dtijdu3[1][1];
    dRdU[k][1][2][4] = dtijdu4[1][1];

    dRdU[k][1][3][0] = dtijdu0[1][2];
    dRdU[k][1][3][1] = dtijdu1[1][2];
    dRdU[k][1][3][2] = dtijdu2[1][2];
    dRdU[k][1][3][3] = dtijdu3[1][2];
    dRdU[k][1][3][4] = dtijdu4[1][2];

    dRdU[k][1][4][0] = ducgdu0[0] * tij[0][1] + ucg[0] * dtijdu0[0][1] + ducgdu0[1] * tij[1][1] + ucg[1] * dtijdu0[1][1] + ducgdu0[2] * tij[1][2] + ucg[2] * dtijdu0[1][2] - dqjdu0[1];
    dRdU[k][1][4][1] = ducgdu1[0] * tij[0][1] + ucg[0] * dtijdu1[0][1] + ducgdu1[1] * tij[1][1] + ucg[1] * dtijdu1[1][1] + ducgdu1[2] * tij[1][2] + ucg[2] * dtijdu1[1][2] - dqjdu1[1];
    dRdU[k][1][4][2] = ducgdu2[0] * tij[0][1] + ucg[0] * dtijdu2[0][1] + ducgdu2[1] * tij[1][1] + ucg[1] * dtijdu2[1][1] + ducgdu2[2] * tij[1][2] + ucg[2] * dtijdu2[1][2] - dqjdu2[1];
    dRdU[k][1][4][3] = ducgdu3[0] * tij[0][1] + ucg[0] * dtijdu3[0][1] + ducgdu3[1] * tij[1][1] + ucg[1] * dtijdu3[1][1] + ducgdu3[2] * tij[1][2] + ucg[2] * dtijdu3[1][2] - dqjdu3[1];
    dRdU[k][1][4][4] = ducgdu4[0] * tij[0][1] + ucg[0] * dtijdu4[0][1] + ducgdu4[1] * tij[1][1] + ucg[1] * dtijdu4[1][1] + ducgdu4[2] * tij[1][2] + ucg[2] * dtijdu4[1][2] - dqjdu4[1];

    // dRzdU

    dRdU[k][2][0][0] = 0.0;
    dRdU[k][2][0][1] = 0.0;
    dRdU[k][2][0][2] = 0.0;
    dRdU[k][2][0][3] = 0.0;
    dRdU[k][2][0][4] = 0.0;

    dRdU[k][2][1][0] = dtijdu0[0][2];
    dRdU[k][2][1][1] = dtijdu1[0][2];
    dRdU[k][2][1][2] = dtijdu2[0][2];
    dRdU[k][2][1][3] = dtijdu3[0][2];
    dRdU[k][2][1][4] = dtijdu4[0][2];

    dRdU[k][2][2][0] = dtijdu0[1][2];
    dRdU[k][2][2][1] = dtijdu1[1][2];
    dRdU[k][2][2][2] = dtijdu2[1][2];
    dRdU[k][2][2][3] = dtijdu3[1][2];
    dRdU[k][2][2][4] = dtijdu4[1][2];

    dRdU[k][2][3][0] = dtijdu0[2][2];
    dRdU[k][2][3][1] = dtijdu1[2][2];
    dRdU[k][2][3][2] = dtijdu2[2][2];
    dRdU[k][2][3][3] = dtijdu3[2][2];
    dRdU[k][2][3][4] = dtijdu4[2][2];

    dRdU[k][2][4][0] = ducgdu0[0] * tij[0][2] + ucg[0] * dtijdu0[0][2] + ducgdu0[1] * tij[1][2] + ucg[1] * dtijdu0[1][2] + ducgdu0[2] * tij[2][2] + ucg[2] * dtijdu0[2][2] - dqjdu0[2];
    dRdU[k][2][4][1] = ducgdu1[0] * tij[0][2] + ucg[0] * dtijdu1[0][2] + ducgdu1[1] * tij[1][2] + ucg[1] * dtijdu1[1][2] + ducgdu1[2] * tij[2][2] + ucg[2] * dtijdu1[2][2] - dqjdu1[2];
    dRdU[k][2][4][2] = ducgdu2[0] * tij[0][2] + ucg[0] * dtijdu2[0][2] + ducgdu2[1] * tij[1][2] + ucg[1] * dtijdu2[1][2] + ducgdu2[2] * tij[2][2] + ucg[2] * dtijdu2[2][2] - dqjdu2[2];
    dRdU[k][2][4][3] = ducgdu3[0] * tij[0][2] + ucg[0] * dtijdu3[0][2] + ducgdu3[1] * tij[1][2] + ucg[1] * dtijdu3[1][2] + ducgdu3[2] * tij[2][2] + ucg[2] * dtijdu3[2][2] - dqjdu3[2];
    dRdU[k][2][4][4] = ducgdu4[0] * tij[0][2] + ucg[0] * dtijdu4[0][2] + ducgdu4[1] * tij[1][2] + ucg[1] * dtijdu4[1][2] + ducgdu4[2] * tij[2][2] + ucg[2] * dtijdu4[2][2] - dqjdu4[2];

  }

}

//------------------------------------------------------------------------------

// Included (MB*)
template<int neq>
inline
void NavierStokesTerm::
computeJacobianVolumeTermNS(double dp1dxj[4][3], double mu, double dmu[4][6], double lambda, double dlambda[4][6], double kappa, double dkappa[4][6], 
			    double *V[4], double T[4], double (*dRdU)[3][neq][neq])
{

  double dTdxj[3];
  computeTemperatureGradient(dp1dxj, T, dTdxj);

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double tij[3][3];
  computeStressTensor(mu, lambda, dudxj, tij);

  double dtijdu0[3][3];
  double dtijdu1[3][3];
  double dtijdu2[3][3];
  double dtijdu3[3][3];
  double dtijdu4[3][3];
  double dtijdu5[3][3];

  double ducgdu0[3];
  double ducgdu1[3];
  double ducgdu2[3];
  double ducgdu3[3];
  double ducgdu4[3];
  double ducgdu5[3];

  double dTdu0[4];
  double dTdu1[4];
  double dTdu2[4];
  double dTdu3[4];
  double dTdu4[4];
  double dTdu5[4];
  
  double ddTdxjdu0[3];
  double ddTdxjdu1[3];
  double ddTdxjdu2[3];
  double ddTdxjdu3[3];
  double ddTdxjdu4[3];
  double ddTdxjdu5[3];

  double dqjdu0[3];
  double dqjdu1[3];
  double dqjdu2[3];
  double dqjdu3[3];
  double dqjdu4[3];
  double dqjdu5[3];

  for (int k=0; k<4; ++k) {

    ducgdu0[0] = - fourth * V[k][1] / V[k][0];
    ducgdu1[0] = fourth / V[k][0];
    ducgdu2[0] = 0.0;
    ducgdu3[0] = 0.0;
    ducgdu4[0] = 0.0;
    ducgdu5[0] = 0.0;

    ducgdu0[1] = - fourth * V[k][2] / V[k][0];
    ducgdu1[1] = 0.0;               
    ducgdu2[1] = fourth / V[k][0];
    ducgdu3[1] = 0.0;
    ducgdu4[1] = 0.0;
    ducgdu5[1] = 0.0;

    ducgdu0[2] = - fourth * V[k][3] / V[k][0];
    ducgdu1[2] = 0.0;
    ducgdu2[2] = 0.0;
    ducgdu3[2] = fourth / V[k][0];
    ducgdu4[2] = 0.0;
    ducgdu5[2] = 0.0;

// Take a look on this piece of code and rederive the stress tensor
/*
    double nu = mu / rho;
    double lu = lambda / rho;

    double dtxxdu0 = -(2.0 * nu * dp1dx * u + lu * (dp1dx * u + dp1dy * v + dp1dz * w));
    double dtxxdu1 = (2.0 * nu + lu) * dp1dx;
    double dtxxdu2 = lu * dp1dy;
    double dtxxdu3 = lu * dp1dz;

    double dtyydu0 = -(2.0 * nu * dp1dy * v + lu * (dp1dx * u + dp1dy * v + dp1dz * w));
    double dtyydu1 = lu * dp1dx;
    double dtyydu2 = (2.0 * nu + lu) * dp1dy;
    double dtyydu3 = lu * dp1dz;

    double dtzzdu0 = -(2.0 * nu * dp1dz * w + lu * (dp1dx * u + dp1dy * v + dp1dz * w));
    double dtzzdu1 = lu * dp1dx;
    double dtzzdu2 = lu * dp1dy;
    double dtzzdu3 = (2.0 * nu + lu) * dp1dz;

    double dtxydu0 = - nu * (dp1dx * v + dp1dy * u);
    double dtxydu1 = nu * dp1dy;
    double dtxydu2 = nu * dp1dx;

    double dtxzdu0 = -nu * (dp1dx * w + dp1dz * u);
    double dtxzdu1 = nu * dp1dz;
    double dtxzdu3 = nu * dp1dx;

    double dtyzdu0 = -nu * (dp1dy * w + dp1dz * v);
    double dtyzdu2 = nu * dp1dz;
    double dtyzdu3 = nu * dp1dy;
*/
//    dtijdu0[0][0] = twothird * mu * (-2.0 * dp1dxj[k][0] * V[k][1] + dp1dxj[k][1] * V[k][2] + dp1dxj[k][2] * V[k][3]) / V[k][0] + dmu[k][0] * twothird * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);
//    dtijdu1[0][0] = 2.0 * twothird * mu * dp1dxj[k][0] / V[k][0] + dmu[k][1] * twothird * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);
//    dtijdu2[0][0] = - twothird * mu * dp1dxj[k][1] / V[k][0] + dmu[k][2] * twothird * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);
//    dtijdu3[0][0] = - twothird * mu * dp1dxj[k][2] / V[k][0] + dmu[k][3] * twothird * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);
//    dtijdu4[0][0] = dmu[k][4] * twothird * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);
//    dtijdu5[0][0] = dmu[k][5] * twothird * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);

//    dtijdu0[1][1] = twothird * mu * (-2.0 * dp1dxj[k][1] * V[k][2] + dp1dxj[k][0] * V[k][1] + dp1dxj[k][2] * V[k][3]) / V[k][0] + dmu[k][0] * twothird * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);
//    dtijdu1[1][1] = - twothird * mu * dp1dxj[k][0] / V[k][0] + dmu[k][1] * twothird * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);
//    dtijdu2[1][1] = 2.0 * twothird * mu * dp1dxj[k][1] / V[k][0] + dmu[k][2] * twothird * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);
//    dtijdu3[1][1] = - twothird * mu * dp1dxj[k][2] / V[k][0] + dmu[k][3] * twothird * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);
//    dtijdu4[1][1] = dmu[k][4] * twothird * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);
//    dtijdu5[1][1] = dmu[k][5] * twothird * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);

//    dtijdu0[2][2] = twothird * mu * (-2.0 * dp1dxj[k][2] * V[k][3] + dp1dxj[k][0] * V[k][1] + dp1dxj[k][1] * V[k][2]) / V[k][0] + dmu[k][0] * twothird * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);
//    dtijdu1[2][2] = - twothird * mu * dp1dxj[k][0] / V[k][0] + dmu[k][1] * twothird * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);
//    dtijdu2[2][2] = - twothird * mu * dp1dxj[k][1] / V[k][0] + dmu[k][2] * twothird * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);
//    dtijdu3[2][2] = 2.0 * twothird * mu * dp1dxj[k][2] / V[k][0] + dmu[k][3] * twothird * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);
//    dtijdu4[2][2] = dmu[k][4] * twothird * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);
//    dtijdu5[2][2] = dmu[k][5] * twothird * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);

    dtijdu0[0][0] = 2.0 * dmu[k][0] * dudxj[0][0] - 2.0 * mu * dp1dxj[k][0] * V[k][1] / V[k][0] + dlambda[k][0] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) - lambda * (dp1dxj[k][0] * V[k][1] + dp1dxj[k][1] * V[k][2] + dp1dxj[k][2] * V[k][3]) / V[k][0];
    dtijdu1[0][0] = 2.0 * dmu[k][1] * dudxj[0][0] + 2.0 * mu * dp1dxj[k][0] / V[k][0] + dlambda[k][1] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][0] / V[k][0];
    dtijdu2[0][0] = 2.0 * dmu[k][2] * dudxj[0][0] + dlambda[k][2] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][1] / V[k][0];
    dtijdu3[0][0] = 2.0 * dmu[k][3] * dudxj[0][0] + dlambda[k][3] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][2] / V[k][0];
    dtijdu4[0][0] = 2.0 * dmu[k][4] * dudxj[0][0] + dlambda[k][4] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]);
    dtijdu5[0][0] = 2.0 * dmu[k][5] * dudxj[0][0] + dlambda[k][5] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]);

    dtijdu0[1][1] = 2.0 * dmu[k][0] * dudxj[1][1] - 2.0 * mu * dp1dxj[k][1] * V[k][2] / V[k][0] + dlambda[k][0] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) - lambda * (dp1dxj[k][1] * V[k][2] + dp1dxj[k][0] * V[k][1] + dp1dxj[k][2] * V[k][3]) / V[k][0];
    dtijdu1[1][1] = 2.0 * dmu[k][1] * dudxj[1][1] + dlambda[k][1] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][0] / V[k][0];
    dtijdu2[1][1] = 2.0 * dmu[k][2] * dudxj[1][1] + 2.0 * mu * dp1dxj[k][1] / V[k][0] + dlambda[k][2] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][1] / V[k][0];
    dtijdu3[1][1] = 2.0 * dmu[k][3] * dudxj[1][1] + dlambda[k][3] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][2] / V[k][0];
    dtijdu4[1][1] = 2.0 * dmu[k][4] * dudxj[1][1] + dlambda[k][4] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]);
    dtijdu5[1][1] = 2.0 * dmu[k][5] * dudxj[1][1] + dlambda[k][5] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]);

    dtijdu0[2][2] = 2.0 * dmu[k][0] * dudxj[2][2] - 2.0 * mu * dp1dxj[k][2] * V[k][3] / V[k][0] + dlambda[k][0] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) - lambda * (dp1dxj[k][2] * V[k][3] + dp1dxj[k][0] * V[k][1] + dp1dxj[k][1] * V[k][2]) / V[k][0];
    dtijdu1[2][2] = 2.0 * dmu[k][1] * dudxj[2][2] + dlambda[k][1] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][0] / V[k][0];
    dtijdu2[2][2] = 2.0 * dmu[k][2] * dudxj[2][2] + dlambda[k][2] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][1] / V[k][0];
    dtijdu3[2][2] = 2.0 * dmu[k][3] * dudxj[2][2] + 2.0 * mu * dp1dxj[k][2] / V[k][0] + dlambda[k][3] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][2] / V[k][0];
    dtijdu4[2][2] = 2.0 * dmu[k][4] * dudxj[2][2] + dlambda[k][4] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]);
    dtijdu5[2][2] = 2.0 * dmu[k][5] * dudxj[2][2] + dlambda[k][5] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]);

    dtijdu0[0][1] = - mu * (dp1dxj[k][0] * V[k][2] + dp1dxj[k][1] * V[k][1]) / V[k][0] + dmu[k][0] * (dudxj[1][0] + dudxj[0][1]);
    dtijdu1[0][1] = mu * dp1dxj[k][1] / V[k][0] + dmu[k][1] * (dudxj[1][0] + dudxj[0][1]);
    dtijdu2[0][1] = mu * dp1dxj[k][0] / V[k][0] + dmu[k][2] * (dudxj[1][0] + dudxj[0][1]);
    dtijdu3[0][1] = dmu[k][3] * (dudxj[1][0] + dudxj[0][1]);
    dtijdu4[0][1] = dmu[k][4] * (dudxj[1][0] + dudxj[0][1]);
    dtijdu5[0][1] = dmu[k][5] * (dudxj[1][0] + dudxj[0][1]);

    dtijdu0[0][2] = - mu * (dp1dxj[k][0] * V[k][3] + dp1dxj[k][2] * V[k][1]) / V[k][0] + dmu[k][0] * (dudxj[2][0] + dudxj[0][2]);
    dtijdu1[0][2] = mu * dp1dxj[k][2] / V[k][0] + dmu[k][1] * (dudxj[2][0] + dudxj[0][2]);
    dtijdu2[0][2] = dmu[k][2] * (dudxj[2][0] + dudxj[0][2]);
    dtijdu3[0][2] = mu * dp1dxj[k][0] / V[k][0] + dmu[k][3] * (dudxj[2][0] + dudxj[0][2]);
    dtijdu4[0][2] = dmu[k][4] * (dudxj[2][0] + dudxj[0][2]);
    dtijdu5[0][2] = dmu[k][5] * (dudxj[2][0] + dudxj[0][2]);

    dtijdu0[1][2] = - mu * (dp1dxj[k][1] * V[k][3] + dp1dxj[k][2] * V[k][2]) / V[k][0] + dmu[k][0] * (dudxj[2][1] + dudxj[1][2]);
    dtijdu1[1][2] = dmu[k][1] * (dudxj[2][1] + dudxj[1][2]);
    dtijdu2[1][2] = mu * dp1dxj[k][2] / V[k][0] + dmu[k][2] * (dudxj[2][1] + dudxj[1][2]);
    dtijdu3[1][2] = mu * dp1dxj[k][1] / V[k][0] + dmu[k][3] * (dudxj[2][1] + dudxj[1][2]);
    dtijdu4[1][2] = dmu[k][4] * (dudxj[2][1] + dudxj[1][2]);
    dtijdu5[1][2] = dmu[k][5] * (dudxj[2][1] + dudxj[1][2]);

    dtijdu0[1][0] = dtijdu0[0][1];
    dtijdu1[1][0] = dtijdu1[0][1];
    dtijdu2[1][0] = dtijdu2[0][1];
    dtijdu3[1][0] = dtijdu3[0][1];
    dtijdu4[1][0] = dtijdu4[0][1];
    dtijdu5[1][0] = dtijdu5[0][1];

    dtijdu0[2][0] = dtijdu0[0][2];
    dtijdu1[2][0] = dtijdu1[0][2];
    dtijdu2[2][0] = dtijdu2[0][2];
    dtijdu3[2][0] = dtijdu3[0][2];
    dtijdu4[2][0] = dtijdu4[0][2];
    dtijdu5[2][0] = dtijdu5[0][2];

    dtijdu0[2][1] = dtijdu0[1][2];
    dtijdu1[2][1] = dtijdu1[1][2];
    dtijdu2[2][1] = dtijdu2[1][2];
    dtijdu3[2][1] = dtijdu3[1][2];
    dtijdu4[2][1] = dtijdu4[1][2];
    dtijdu5[2][1] = dtijdu5[1][2];

    dTdu0[k] = - 1.0 / V[k][0] * (T[k] - 0.5 * (V[k][1]*V[k][1] + V[k][2]*V[k][2] + V[k][3]*V[k][3]));
    dTdu1[k] = - 1.0 / V[k][0] * V[k][1];
    dTdu2[k] = - 1.0 / V[k][0] * V[k][2];
    dTdu3[k] = - 1.0 / V[k][0] * V[k][3];
    dTdu4[k] = 1.0 / V[k][0];
    dTdu5[k] = 0.0;

    ddTdxjdu0[0] = dp1dxj[k][0]*dTdu0[k];
    ddTdxjdu1[0] = dp1dxj[k][0]*dTdu1[k];
    ddTdxjdu2[0] = dp1dxj[k][0]*dTdu2[k];
    ddTdxjdu3[0] = dp1dxj[k][0]*dTdu3[k];
    ddTdxjdu4[0] = dp1dxj[k][0]*dTdu4[k];
    ddTdxjdu5[0] = dp1dxj[k][0]*dTdu5[k];

    ddTdxjdu0[1] = dp1dxj[k][1]*dTdu0[k];
    ddTdxjdu1[1] = dp1dxj[k][1]*dTdu1[k];
    ddTdxjdu2[1] = dp1dxj[k][1]*dTdu2[k];
    ddTdxjdu3[1] = dp1dxj[k][1]*dTdu3[k];
    ddTdxjdu4[1] = dp1dxj[k][1]*dTdu4[k];
    ddTdxjdu5[1] = dp1dxj[k][1]*dTdu5[k];
    
    ddTdxjdu0[2] = dp1dxj[k][2]*dTdu0[k];
    ddTdxjdu1[2] = dp1dxj[k][2]*dTdu1[k];
    ddTdxjdu2[2] = dp1dxj[k][2]*dTdu2[k];
    ddTdxjdu3[2] = dp1dxj[k][2]*dTdu3[k];
    ddTdxjdu4[2] = dp1dxj[k][2]*dTdu4[k];
    ddTdxjdu5[2] = dp1dxj[k][2]*dTdu5[k];
     
    dqjdu0[0] = - dkappa[k][0] * dTdxj[0] - kappa * ddTdxjdu0[0];
    dqjdu1[0] = - dkappa[k][1] * dTdxj[0] - kappa * ddTdxjdu1[0];
    dqjdu2[0] = - dkappa[k][2] * dTdxj[0] - kappa * ddTdxjdu2[0];
    dqjdu3[0] = - dkappa[k][3] * dTdxj[0] - kappa * ddTdxjdu3[0];
    dqjdu4[0] = - dkappa[k][4] * dTdxj[0] - kappa * ddTdxjdu4[0];
    dqjdu5[0] = - dkappa[k][5] * dTdxj[0] - kappa * ddTdxjdu5[0];
    
    dqjdu0[1] = - dkappa[k][0] * dTdxj[1] - kappa * ddTdxjdu0[1];
    dqjdu1[1] = - dkappa[k][1] * dTdxj[1] - kappa * ddTdxjdu1[1];
    dqjdu2[1] = - dkappa[k][2] * dTdxj[1] - kappa * ddTdxjdu2[1];
    dqjdu3[1] = - dkappa[k][3] * dTdxj[1] - kappa * ddTdxjdu3[1];
    dqjdu4[1] = - dkappa[k][4] * dTdxj[1] - kappa * ddTdxjdu4[1];
    dqjdu5[1] = - dkappa[k][5] * dTdxj[1] - kappa * ddTdxjdu5[1];

    dqjdu0[2] = - dkappa[k][0] * dTdxj[2] - kappa * ddTdxjdu0[2];
    dqjdu1[2] = - dkappa[k][1] * dTdxj[2] - kappa * ddTdxjdu1[2];
    dqjdu2[2] = - dkappa[k][2] * dTdxj[2] - kappa * ddTdxjdu2[2];
    dqjdu3[2] = - dkappa[k][3] * dTdxj[2] - kappa * ddTdxjdu3[2];
    dqjdu4[2] = - dkappa[k][4] * dTdxj[2] - kappa * ddTdxjdu4[2];
    dqjdu5[2] = - dkappa[k][5] * dTdxj[2] - kappa * ddTdxjdu5[2];

    // dRxdU

    dRdU[k][0][0][0] = 0.0;
    dRdU[k][0][0][1] = 0.0;
    dRdU[k][0][0][2] = 0.0;
    dRdU[k][0][0][3] = 0.0;
    dRdU[k][0][0][4] = 0.0;
    dRdU[k][0][0][5] = 0.0;

    dRdU[k][0][1][0] = dtijdu0[0][0];
    dRdU[k][0][1][1] = dtijdu1[0][0];
    dRdU[k][0][1][2] = dtijdu2[0][0];
    dRdU[k][0][1][3] = dtijdu3[0][0];
    dRdU[k][0][1][4] = dtijdu4[0][0];
    dRdU[k][0][1][5] = dtijdu5[0][0];

    dRdU[k][0][2][0] = dtijdu0[0][1];
    dRdU[k][0][2][1] = dtijdu1[0][1];
    dRdU[k][0][2][2] = dtijdu2[0][1];
    dRdU[k][0][2][3] = dtijdu3[0][1];
    dRdU[k][0][2][4] = dtijdu4[0][1];
    dRdU[k][0][2][5] = dtijdu5[0][1];

    dRdU[k][0][3][0] = dtijdu0[0][2];
    dRdU[k][0][3][1] = dtijdu1[0][2];
    dRdU[k][0][3][2] = dtijdu2[0][2];
    dRdU[k][0][3][3] = dtijdu3[0][2];
    dRdU[k][0][3][4] = dtijdu4[0][2];
    dRdU[k][0][3][5] = dtijdu5[0][2];

    dRdU[k][0][4][0] = ducgdu0[0] * tij[0][0] + ucg[0] * dtijdu0[0][0] + ducgdu0[1] * tij[0][1] + ucg[1] * dtijdu0[0][1] + ducgdu0[2] * tij[0][2] + ucg[2] * dtijdu0[0][2] - dqjdu0[0];
    dRdU[k][0][4][1] = ducgdu1[0] * tij[0][0] + ucg[0] * dtijdu1[0][0] + ducgdu1[1] * tij[0][1] + ucg[1] * dtijdu1[0][1] + ducgdu1[2] * tij[0][2] + ucg[2] * dtijdu1[0][2] - dqjdu1[0];
    dRdU[k][0][4][2] = ducgdu2[0] * tij[0][0] + ucg[0] * dtijdu2[0][0] + ducgdu2[1] * tij[0][1] + ucg[1] * dtijdu2[0][1] + ducgdu2[2] * tij[0][2] + ucg[2] * dtijdu2[0][2] - dqjdu2[0];
    dRdU[k][0][4][3] = ducgdu3[0] * tij[0][0] + ucg[0] * dtijdu3[0][0] + ducgdu3[1] * tij[0][1] + ucg[1] * dtijdu3[0][1] + ducgdu3[2] * tij[0][2] + ucg[2] * dtijdu3[0][2] - dqjdu3[0];
    dRdU[k][0][4][4] = ducgdu4[0] * tij[0][0] + ucg[0] * dtijdu4[0][0] + ducgdu4[1] * tij[0][1] + ucg[1] * dtijdu4[0][1] + ducgdu4[2] * tij[0][2] + ucg[2] * dtijdu4[0][2] - dqjdu4[0];
    dRdU[k][0][4][5] = ducgdu5[0] * tij[0][0] + ucg[0] * dtijdu5[0][0] + ducgdu5[1] * tij[0][1] + ucg[1] * dtijdu5[0][1] + ducgdu5[2] * tij[0][2] + ucg[2] * dtijdu5[0][2] - dqjdu5[0];

    // dRydU

    dRdU[k][1][0][0] = 0.0;
    dRdU[k][1][0][1] = 0.0;
    dRdU[k][1][0][2] = 0.0;
    dRdU[k][1][0][3] = 0.0;
    dRdU[k][1][0][4] = 0.0;
    dRdU[k][1][0][5] = 0.0;
  
    dRdU[k][1][1][0] = dtijdu0[0][1];
    dRdU[k][1][1][1] = dtijdu1[0][1];
    dRdU[k][1][1][2] = dtijdu2[0][1];
    dRdU[k][1][1][3] = dtijdu3[0][1];
    dRdU[k][1][1][4] = dtijdu4[0][1];
    dRdU[k][1][1][5] = dtijdu5[0][1];

    dRdU[k][1][2][0] = dtijdu0[1][1];
    dRdU[k][1][2][1] = dtijdu1[1][1];
    dRdU[k][1][2][2] = dtijdu2[1][1];
    dRdU[k][1][2][3] = dtijdu3[1][1];
    dRdU[k][1][2][4] = dtijdu4[1][1];
    dRdU[k][1][2][5] = dtijdu5[1][1];

    dRdU[k][1][3][0] = dtijdu0[1][2];
    dRdU[k][1][3][1] = dtijdu1[1][2];
    dRdU[k][1][3][2] = dtijdu2[1][2];
    dRdU[k][1][3][3] = dtijdu3[1][2];
    dRdU[k][1][3][4] = dtijdu4[1][2];
    dRdU[k][1][3][5] = dtijdu5[1][2];

    dRdU[k][1][4][0] = ducgdu0[0] * tij[0][1] + ucg[0] * dtijdu0[0][1] + ducgdu0[1] * tij[1][1] + ucg[1] * dtijdu0[1][1] + ducgdu0[2] * tij[1][2] + ucg[2] * dtijdu0[1][2] - dqjdu0[1];
    dRdU[k][1][4][1] = ducgdu1[0] * tij[0][1] + ucg[0] * dtijdu1[0][1] + ducgdu1[1] * tij[1][1] + ucg[1] * dtijdu1[1][1] + ducgdu1[2] * tij[1][2] + ucg[2] * dtijdu1[1][2] - dqjdu1[1];
    dRdU[k][1][4][2] = ducgdu2[0] * tij[0][1] + ucg[0] * dtijdu2[0][1] + ducgdu2[1] * tij[1][1] + ucg[1] * dtijdu2[1][1] + ducgdu2[2] * tij[1][2] + ucg[2] * dtijdu2[1][2] - dqjdu2[1];
    dRdU[k][1][4][3] = ducgdu3[0] * tij[0][1] + ucg[0] * dtijdu3[0][1] + ducgdu3[1] * tij[1][1] + ucg[1] * dtijdu3[1][1] + ducgdu3[2] * tij[1][2] + ucg[2] * dtijdu3[1][2] - dqjdu3[1];
    dRdU[k][1][4][4] = ducgdu4[0] * tij[0][1] + ucg[0] * dtijdu4[0][1] + ducgdu4[1] * tij[1][1] + ucg[1] * dtijdu4[1][1] + ducgdu4[2] * tij[1][2] + ucg[2] * dtijdu4[1][2] - dqjdu4[1];
    dRdU[k][1][4][5] = ducgdu5[0] * tij[0][1] + ucg[0] * dtijdu5[0][1] + ducgdu5[1] * tij[1][1] + ucg[1] * dtijdu5[1][1] + ducgdu5[2] * tij[1][2] + ucg[2] * dtijdu5[1][2] - dqjdu5[1];

    // dRzdU

    dRdU[k][2][0][0] = 0.0;
    dRdU[k][2][0][1] = 0.0;
    dRdU[k][2][0][2] = 0.0;
    dRdU[k][2][0][3] = 0.0;
    dRdU[k][2][0][4] = 0.0;
    dRdU[k][2][0][5] = 0.0;

    dRdU[k][2][1][0] = dtijdu0[0][2];
    dRdU[k][2][1][1] = dtijdu1[0][2];
    dRdU[k][2][1][2] = dtijdu2[0][2];
    dRdU[k][2][1][3] = dtijdu3[0][2];
    dRdU[k][2][1][4] = dtijdu4[0][2];
    dRdU[k][2][1][5] = dtijdu5[0][2];

    dRdU[k][2][2][0] = dtijdu0[1][2];
    dRdU[k][2][2][1] = dtijdu1[1][2];
    dRdU[k][2][2][2] = dtijdu2[1][2];
    dRdU[k][2][2][3] = dtijdu3[1][2];
    dRdU[k][2][2][4] = dtijdu4[1][2];
    dRdU[k][2][2][5] = dtijdu5[1][2];

    dRdU[k][2][3][0] = dtijdu0[2][2];
    dRdU[k][2][3][1] = dtijdu1[2][2];
    dRdU[k][2][3][2] = dtijdu2[2][2];
    dRdU[k][2][3][3] = dtijdu3[2][2];
    dRdU[k][2][3][4] = dtijdu4[2][2];
    dRdU[k][2][3][5] = dtijdu5[2][2];

    dRdU[k][2][4][0] = ducgdu0[0] * tij[0][2] + ucg[0] * dtijdu0[0][2] + ducgdu0[1] * tij[1][2] + ucg[1] * dtijdu0[1][2] + ducgdu0[2] * tij[2][2] + ucg[2] * dtijdu0[2][2] - dqjdu0[2];
    dRdU[k][2][4][1] = ducgdu1[0] * tij[0][2] + ucg[0] * dtijdu1[0][2] + ducgdu1[1] * tij[1][2] + ucg[1] * dtijdu1[1][2] + ducgdu1[2] * tij[2][2] + ucg[2] * dtijdu1[2][2] - dqjdu1[2];
    dRdU[k][2][4][2] = ducgdu2[0] * tij[0][2] + ucg[0] * dtijdu2[0][2] + ducgdu2[1] * tij[1][2] + ucg[1] * dtijdu2[1][2] + ducgdu2[2] * tij[2][2] + ucg[2] * dtijdu2[2][2] - dqjdu2[2];
    dRdU[k][2][4][3] = ducgdu3[0] * tij[0][2] + ucg[0] * dtijdu3[0][2] + ducgdu3[1] * tij[1][2] + ucg[1] * dtijdu3[1][2] + ducgdu3[2] * tij[2][2] + ucg[2] * dtijdu3[2][2] - dqjdu3[2];
    dRdU[k][2][4][4] = ducgdu4[0] * tij[0][2] + ucg[0] * dtijdu4[0][2] + ducgdu4[1] * tij[1][2] + ucg[1] * dtijdu4[1][2] + ducgdu4[2] * tij[2][2] + ucg[2] * dtijdu4[2][2] - dqjdu4[2];
    dRdU[k][2][4][5] = ducgdu5[0] * tij[0][2] + ucg[0] * dtijdu5[0][2] + ducgdu5[1] * tij[1][2] + ucg[1] * dtijdu5[1][2] + ducgdu5[2] * tij[2][2] + ucg[2] * dtijdu5[2][2] - dqjdu5[2];

  }

}

//------------------------------------------------------------------------------

template<int neq>
inline
void NavierStokesTerm::
computeJacobianVolumeTermNS(double dp1dxj[4][3], double mu, double lambda, double kappa, 
			    double *V[4], double T[4], double (*dRdU)[3][neq][neq])
{

  double uu[4][3], ucg[3];
  computeVelocity(V, uu, ucg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, uu, dudxj);

  double tij[3][3];
  computeStressTensor(mu, lambda, dudxj, tij);

  double txx = tij[0][0];
  double tyy = tij[1][1];
  double tzz = tij[2][2];
  double txy = tij[0][1];
  double txz = tij[0][2];
  double tyz = tij[1][2];

  for (int k=0; k<4; ++k) {
    
    double rho = V[k][0];
    double u = V[k][1];
    double v = V[k][2];
    double w = V[k][3];
    double Temp = T[k];

    double dp1dx = dp1dxj[k][0];
    double dp1dy = dp1dxj[k][1];
    double dp1dz = dp1dxj[k][2];

    double nu = mu / rho;
    double lu = lambda / rho;

/*//---------------------
    double adtxxdu0 = twothird * nu * (-2.0 * dp1dx * u + dp1dy * v + dp1dz * w);
    double adtxxdu1 = 2.0 * twothird * nu * dp1dx;
    double adtxxdu2 = - twothird * nu * dp1dy;
    double adtxxdu3 = - twothird * nu * dp1dz;

    double adtyydu0 = twothird * nu * (-2.0 * dp1dy * v + dp1dx * u + dp1dz * w);
    double adtyydu1 = - twothird * nu * dp1dx;
    double adtyydu2 = 2.0 * twothird * nu * dp1dy;
    double adtyydu3 = - twothird * nu * dp1dz;

    double adtzzdu0 = twothird * nu * (-2.0 * dp1dz * w + dp1dx * u + dp1dy * v);
    double adtzzdu1 = - twothird * nu * dp1dx;
    double adtzzdu2 = - twothird * nu * dp1dy;
    double adtzzdu3 = 2.0 * twothird * nu * dp1dz;

    double adtxydu0 = - nu * (dp1dx * v + dp1dy * u);
    double adtxydu1 = nu * dp1dy;
    double adtxydu2 = nu * dp1dx;

    double adtxzdu0 = -nu * (dp1dx * w + dp1dz * u);
    double adtxzdu1 = nu * dp1dz;
    double adtxzdu3 = nu * dp1dx;

    double adtyzdu0 = -nu * (dp1dy * w + dp1dz * v);
    double adtyzdu2 = nu * dp1dz;
    double adtyzdu3 = nu * dp1dy;

//--------------------*/

    double dtxxdu0 = -(2.0 * nu * dp1dx * u + lu * (dp1dx * u + dp1dy * v + dp1dz * w));
    double dtxxdu1 = (2.0 * nu + lu) * dp1dx;
    double dtxxdu2 = lu * dp1dy;
    double dtxxdu3 = lu * dp1dz;

    double dtyydu0 = -(2.0 * nu * dp1dy * v + lu * (dp1dx * u + dp1dy * v + dp1dz * w));
    double dtyydu1 = lu * dp1dx;
    double dtyydu2 = (2.0 * nu + lu) * dp1dy;
    double dtyydu3 = lu * dp1dz;

    double dtzzdu0 = -(2.0 * nu * dp1dz * w + lu * (dp1dx * u + dp1dy * v + dp1dz * w));
    double dtzzdu1 = lu * dp1dx;
    double dtzzdu2 = lu * dp1dy;
    double dtzzdu3 = (2.0 * nu + lu) * dp1dz;

    double dtxydu0 = - nu * (dp1dx * v + dp1dy * u);
    double dtxydu1 = nu * dp1dy;
    double dtxydu2 = nu * dp1dx;

    double dtxzdu0 = -nu * (dp1dx * w + dp1dz * u);
    double dtxzdu1 = nu * dp1dz;
    double dtxzdu3 = nu * dp1dx;

    double dtyzdu0 = -nu * (dp1dy * w + dp1dz * v);
    double dtyzdu2 = nu * dp1dz;
    double dtyzdu3 = nu * dp1dy;

    double alpha = kappa / rho;
    double c0 = alpha * (Temp - 0.5 * (u*u + v*v + w*w));
    double c1 = alpha * u;
    double c2 = alpha * v;
    double c3 = alpha * w;
    double c4 = - alpha;

    double dqxdu0 = c0 * dp1dx;
    double dqxdu1 = c1 * dp1dx;
    double dqxdu2 = c2 * dp1dx;
    double dqxdu3 = c3 * dp1dx;
    double dqxdu4 = c4 * dp1dx;

    double dqydu0 = c0 * dp1dy;
    double dqydu1 = c1 * dp1dy;
    double dqydu2 = c2 * dp1dy;
    double dqydu3 = c3 * dp1dy;
    double dqydu4 = c4 * dp1dy;

    double dqzdu0 = c0 * dp1dz;
    double dqzdu1 = c1 * dp1dz;
    double dqzdu2 = c2 * dp1dz;
    double dqzdu3 = c3 * dp1dz;
    double dqzdu4 = c4 * dp1dz;

    double beta = 0.25 / rho;

    double dudu0 = - beta * u;
    double dudu1 = beta;

    double dvdu0 = - beta * v;
    double dvdu2 = beta;

    double dwdu0 = - beta * w;
    double dwdu3 = beta;


    //computation of dRdU -> dRdU[k][a][b][c]
// k is the corner of the tetrahedra we are considering
// a is the coordinate direction: x y z
// b is the equation we are considering: 0 for mass conservation
//					 1 for x-momentum
//					 2 for y-momentum
//					 3 for z-momentum
//					 4 for energy
// c is the variable wrt which we derive

    // dRxdU

    dRdU[k][0][0][0] = 0.0;
    dRdU[k][0][0][1] = 0.0;
    dRdU[k][0][0][2] = 0.0;
    dRdU[k][0][0][3] = 0.0;
    dRdU[k][0][0][4] = 0.0;

    dRdU[k][0][1][0] = dtxxdu0;
    dRdU[k][0][1][1] = dtxxdu1;
    dRdU[k][0][1][2] = dtxxdu2;
    dRdU[k][0][1][3] = dtxxdu3;
    dRdU[k][0][1][4] = 0.0;

    dRdU[k][0][2][0] = dtxydu0;
    dRdU[k][0][2][1] = dtxydu1;
    dRdU[k][0][2][2] = dtxydu2;
    dRdU[k][0][2][3] = 0.0;
    dRdU[k][0][2][4] = 0.0;

    dRdU[k][0][3][0] = dtxzdu0;
    dRdU[k][0][3][1] = dtxzdu1;
    dRdU[k][0][3][2] = 0.0;
    dRdU[k][0][3][3] = dtxzdu3;
    dRdU[k][0][3][4] = 0.0;

    dRdU[k][0][4][0] = dudu0 * txx + ucg[0] * dtxxdu0 + dvdu0 * txy + ucg[1] * dtxydu0 + 
                       dwdu0 * txz + ucg[2] * dtxzdu0 - dqxdu0;
    dRdU[k][0][4][1] = dudu1 * txx + ucg[0] * dtxxdu1 +
                       ucg[1] * dtxydu1 + ucg[2] * dtxzdu1 - dqxdu1;
    dRdU[k][0][4][2] = ucg[0] * dtxxdu2 + dvdu2 * txy + ucg[1] * dtxydu2 - dqxdu2;
    dRdU[k][0][4][3] = ucg[0] * dtxxdu3 + dwdu3 * txz + ucg[2] * dtxzdu3 - dqxdu3;
    dRdU[k][0][4][4] = - dqxdu4;

    // dRydU

    dRdU[k][1][0][0] = 0.0;
    dRdU[k][1][0][1] = 0.0;
    dRdU[k][1][0][2] = 0.0;
    dRdU[k][1][0][3] = 0.0;
    dRdU[k][1][0][4] = 0.0;
  
    dRdU[k][1][1][0] = dtxydu0;
    dRdU[k][1][1][1] = dtxydu1;
    dRdU[k][1][1][2] = dtxydu2;
    dRdU[k][1][1][3] = 0.0;
    dRdU[k][1][1][4] = 0.0;

    dRdU[k][1][2][0] = dtyydu0;
    dRdU[k][1][2][1] = dtyydu1;
    dRdU[k][1][2][2] = dtyydu2;
    dRdU[k][1][2][3] = dtyydu3;
    dRdU[k][1][2][4] = 0.0;

    dRdU[k][1][3][0] = dtyzdu0;
    dRdU[k][1][3][1] = 0.0;
    dRdU[k][1][3][2] = dtyzdu2;
    dRdU[k][1][3][3] = dtyzdu3;
    dRdU[k][1][3][4] = 0.0;

    dRdU[k][1][4][0] = dudu0 * txy + ucg[0] * dtxydu0 + dvdu0 * tyy + ucg[1] * dtyydu0 + 
                       dwdu0 * tyz + ucg[2] * dtyzdu0 - dqydu0;
    dRdU[k][1][4][1] = dudu1 * txy + ucg[0] * dtxydu1 + ucg[1] * dtyydu1 - dqydu1;
    dRdU[k][1][4][2] = ucg[0] * dtxydu2 + dvdu2 * tyy + 
                       ucg[1] * dtyydu2 + ucg[2] * dtyzdu2 - dqydu2;
    dRdU[k][1][4][3] = ucg[1] * dtyydu3 + dwdu3 * tyz + ucg[2] * dtyzdu3 - dqydu3;
    dRdU[k][1][4][4] = - dqydu4;

    // dRzdU

    dRdU[k][2][0][0] = 0.0;
    dRdU[k][2][0][1] = 0.0;
    dRdU[k][2][0][2] = 0.0;
    dRdU[k][2][0][3] = 0.0;
    dRdU[k][2][0][4] = 0.0;

    dRdU[k][2][1][0] = dtxzdu0;
    dRdU[k][2][1][1] = dtxzdu1;
    dRdU[k][2][1][2] = 0.0;
    dRdU[k][2][1][3] = dtxzdu3;
    dRdU[k][2][1][4] = 0.0;

    dRdU[k][2][2][0] = dtyzdu0;
    dRdU[k][2][2][1] = 0.0;
    dRdU[k][2][2][2] = dtyzdu2;
    dRdU[k][2][2][3] = dtyzdu3;
    dRdU[k][2][2][4] = 0.0;

    dRdU[k][2][3][0] = dtzzdu0;
    dRdU[k][2][3][1] = dtzzdu1;
    dRdU[k][2][3][2] = dtzzdu2;
    dRdU[k][2][3][3] = dtzzdu3;
    dRdU[k][2][3][4] = 0.0;

    dRdU[k][2][4][0] = dudu0 * txz + ucg[0] * dtxzdu0 + dvdu0 * tyz + ucg[1] * dtyzdu0 + 
                       dwdu0 * tzz + ucg[2] * dtzzdu0 - dqzdu0;
    dRdU[k][2][4][1] = dudu1 * txz + ucg[0] * dtxzdu1 + ucg[2] * dtzzdu1 - dqzdu1;
    dRdU[k][2][4][2] = dvdu2 * tyz + ucg[1] * dtyzdu2 + ucg[2] * dtzzdu2 - dqzdu2;
    dRdU[k][2][4][3] = ucg[0] * dtxzdu3 + ucg[1] * dtyzdu3 + 
                       dwdu3 * tzz + ucg[2] * dtzzdu3 - dqzdu3;
    dRdU[k][2][4][4] = - dqzdu4;

  }

}

//------------------------------------------------------------------------------

inline
void NavierStokesTerm::computeSurfaceTermNS(double dp1dxj[4][3], Vec3D &n, 
					    double *Vwall, double *Vtet[4], double *R)
{

  double u[4][3], ucg[3];
  computeVelocity(Vtet, u, ucg);

  double T[4], Tcg;
  computeTemperature(Vtet, T, Tcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double mu     = viscoFcn->compute_mu(Tcg);
  double lambda = viscoFcn->compute_lambda(Tcg, mu);
  mu     *= ooreynolds_mu;
  lambda *= ooreynolds_mu;

  double tij[3][3];
  computeStressTensor(mu, lambda, dudxj, tij);

  R[0] = 0.0;
  R[1] = 0.0;
  R[2] = 0.0;
  R[3] = 0.0;
  R[4] = (Vwall[1] * tij[0][0] + Vwall[2] * tij[1][0] + Vwall[3] * tij[2][0]) * n[0] +
    (Vwall[1] * tij[0][1] + Vwall[2] * tij[1][1] + Vwall[3] * tij[2][1]) * n[1] + 
    (Vwall[1] * tij[0][2] + Vwall[2] * tij[1][2] + Vwall[3] * tij[2][2]) * n[2];

}

//------------------------------------------------------------------------------

// Included (MB)
inline
void NavierStokesTerm::computeDerivativeOfSurfaceTermNS(double dp1dxj[4][3], double ddp1dxj[4][3], Vec3D &n, Vec3D &dn,
					    double *Vwall, double *dVwall, double *Vtet[4], double *dVtet[4], double dMach, double *dR)
{

  double u[4][3], ucg[3];
  computeVelocity(Vtet, u, ucg);

  double du[4][3], ducg[3];
  computeDerivativeOfVelocity(dVtet, du, ducg);

  double T[4], Tcg;
  computeTemperature(Vtet, T, Tcg);

  double dT[4], dTcg;
  computeDerivativeOfTemperature(Vtet, dVtet, dT, dTcg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double ddudxj[3][3];
  computeDerivativeOfVelocityGradient(dp1dxj, ddp1dxj, u, du, ddudxj);

  double dooreynolds_mu = -1.0 / ( reynolds_muNS * reynolds_muNS ) * dRe_mudMachNS * dMach;
  
  double mu     = viscoFcn->compute_mu(Tcg);
  double lambda = viscoFcn->compute_lambda(Tcg, mu);

  double dmu     = dooreynolds_mu * mu + ooreynolds_mu * viscoFcn->compute_muDerivative(Tcg, dTcg, dMach);
  double dlambda = dooreynolds_mu * lambda + ooreynolds_mu * viscoFcn->compute_lambdaDerivative(mu, dmu, dMach);

  mu     *= ooreynolds_mu;
  lambda *= ooreynolds_mu;

  double tij[3][3];
  computeStressTensor(mu, lambda, dudxj, tij);

  double dtij[3][3];
  computeDerivativeOfStressTensor(mu, dmu, lambda, dlambda, dudxj, ddudxj, dtij);

  dR[0] = 0.0;
  dR[1] = 0.0;
  dR[2] = 0.0;
  dR[3] = 0.0;
  dR[4] = (Vwall[1] * dtij[0][0] + dVwall[1] * tij[0][0] + Vwall[2] * dtij[1][0] + dVwall[2] * tij[1][0] + Vwall[3] * dtij[2][0] + dVwall[3] * tij[2][0]) * n[0] + (Vwall[1] * tij[0][0] + Vwall[2] * tij[1][0] + Vwall[3] * tij[2][0]) * dn[0] +
          (Vwall[1] * dtij[0][1] + dVwall[1] * tij[0][1] + Vwall[2] * dtij[1][1] + dVwall[2] * tij[1][1] + Vwall[3] * dtij[2][1] + dVwall[3] * tij[2][1]) * n[1] + (Vwall[1] * tij[0][1] + Vwall[2] * tij[1][1] + Vwall[3] * tij[2][1]) * dn[1] +
          (Vwall[1] * dtij[0][2] + dVwall[1] * tij[0][2] + Vwall[2] * dtij[1][2] + dVwall[2] * tij[1][2] + Vwall[3] * dtij[2][2] + dVwall[3] * tij[2][2]) * n[2] + (Vwall[1] * tij[0][2] + Vwall[2] * tij[1][2] + Vwall[3] * tij[2][2]) * dn[2];

}

//------------------------------------------------------------------------------

// Included (MB*)
template<int neq>
inline
void NavierStokesTerm::computeJacobianSurfaceTermNS(double dp1dxj[4][3], Vec3D &n, 
						    double *Vwall, double *V[4], 
						    double (*dRdU)[neq][neq])
{

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double u[4][3], ucg[3];
  computeVelocity(V, u, ucg);

  double dudxj[3][3];
  computeVelocityGradient(dp1dxj, u, dudxj);

  double mu     = viscoFcn->compute_mu(Tcg);
  double lambda = viscoFcn->compute_lambda(Tcg, mu);
  mu     *= ooreynolds_mu;
  lambda *= ooreynolds_mu;

  double dmu[4][5];
  double dlambda[4][5];
  
  double dTcgdu0;
  double dTcgdu1;
  double dTcgdu2;
  double dTcgdu3;
  double dTcgdu4;
  
  double dMach = 0.0;

  double dtijdu0[3][3];
  double dtijdu1[3][3];
  double dtijdu2[3][3];
  double dtijdu3[3][3];
  double dtijdu4[3][3];

  for (int k=0; k<4; ++k) {
    
    dTcgdu0 = - 0.25 / V[k][0] * (T[k] - 0.5 * (V[k][1]*V[k][1] + V[k][2]*V[k][2] + V[k][3]*V[k][3]));
    dTcgdu1 = - 0.25 / V[k][0] * V[k][1];
    dTcgdu2 = - 0.25 / V[k][0] * V[k][2];
    dTcgdu3 = - 0.25 / V[k][0] * V[k][3];
    dTcgdu4 = 0.25 / V[k][0];

    dmu[k][0] = viscoFcn->compute_muDerivative(Tcg, dTcgdu0, dMach);
    dmu[k][1] = viscoFcn->compute_muDerivative(Tcg, dTcgdu1, dMach);
    dmu[k][2] = viscoFcn->compute_muDerivative(Tcg, dTcgdu2, dMach);
    dmu[k][3] = viscoFcn->compute_muDerivative(Tcg, dTcgdu3, dMach); 
    dmu[k][4] = viscoFcn->compute_muDerivative(Tcg, dTcgdu4, dMach); 

    for(int i=0; i<5; ++i){
      dlambda[k][i] = ooreynolds_mu * viscoFcn->compute_lambdaDerivative(mu/ooreynolds_mu, dmu[k][i], dMach);
      dmu[k][i] *= ooreynolds_mu;
    }

//    dtijdu0[0][0] = twothird * mu * (-2.0 * dp1dxj[k][0] * V[k][1] + dp1dxj[k][1] * V[k][2] + dp1dxj[k][2] * V[k][3]) / V[k][0] + dmu[k][0] * twothird * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);
//    dtijdu1[0][0] = 2.0 * twothird * mu * dp1dxj[k][0] / V[k][0] + dmu[k][1] * twothird * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);
//    dtijdu2[0][0] = - twothird * mu * dp1dxj[k][1] / V[k][0] + dmu[k][2] * twothird * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);
//    dtijdu3[0][0] = - twothird * mu * dp1dxj[k][2] / V[k][0] + dmu[k][3] * twothird * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);
//    dtijdu4[0][0] = dmu[k][4] * twothird * (2.0 * dudxj[0][0] - dudxj[1][1] - dudxj[2][2]);

//    dtijdu0[1][1] = twothird * mu * (-2.0 * dp1dxj[k][1] * V[k][2] + dp1dxj[k][0] * V[k][1] + dp1dxj[k][2] * V[k][3]) / V[k][0] + dmu[k][0] * twothird * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);
//    dtijdu1[1][1] = - twothird * mu * dp1dxj[k][0] / V[k][0] + dmu[k][1] * twothird * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);
//    dtijdu2[1][1] = 2.0 * twothird * mu * dp1dxj[k][1] / V[k][0] + dmu[k][2] * twothird * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);
//    dtijdu3[1][1] = - twothird * mu * dp1dxj[k][2] / V[k][0] + dmu[k][3] * twothird * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);
//    dtijdu4[1][1] = dmu[k][4] * twothird * (2.0 * dudxj[1][1] - dudxj[0][0] - dudxj[2][2]);

//    dtijdu0[2][2] = twothird * mu * (-2.0 * dp1dxj[k][2] * V[k][3] + dp1dxj[k][0] * V[k][1] + dp1dxj[k][1] * V[k][2]) / V[k][0] + dmu[k][0] * twothird * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);
//    dtijdu1[2][2] = - twothird * mu * dp1dxj[k][0] / V[k][0] + dmu[k][1] * twothird * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);
//    dtijdu2[2][2] = - twothird * mu * dp1dxj[k][1] / V[k][0] + dmu[k][2] * twothird * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);
//    dtijdu3[2][2] = 2.0 * twothird * mu * dp1dxj[k][2] / V[k][0] + dmu[k][3] * twothird * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);
//    dtijdu4[2][2] = dmu[k][4] * twothird * (2.0 * dudxj[2][2] - dudxj[0][0] - dudxj[1][1]);

    dtijdu0[0][0] = 2.0 * dmu[k][0] * dudxj[0][0] - 2.0 * mu * dp1dxj[k][0] * V[k][1] / V[k][0] + dlambda[k][0] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) - lambda * (dp1dxj[k][0] * V[k][1] + dp1dxj[k][1] * V[k][2] + dp1dxj[k][2] * V[k][3]) / V[k][0];
    dtijdu1[0][0] = 2.0 * dmu[k][1] * dudxj[0][0] + 2.0 * mu * dp1dxj[k][0] / V[k][0] + dlambda[k][1] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][0] / V[k][0];
    dtijdu2[0][0] = 2.0 * dmu[k][2] * dudxj[0][0] + dlambda[k][2] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][1] / V[k][0];
    dtijdu3[0][0] = 2.0 * dmu[k][3] * dudxj[0][0] + dlambda[k][3] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][2] / V[k][0];
    dtijdu4[0][0] = 2.0 * dmu[k][4] * dudxj[0][0] + dlambda[k][4] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]);

    dtijdu0[1][1] = 2.0 * dmu[k][0] * dudxj[1][1] - 2.0 * mu * dp1dxj[k][1] * V[k][2] / V[k][0] + dlambda[k][0] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) - lambda * (dp1dxj[k][1] * V[k][2] + dp1dxj[k][0] * V[k][1] + dp1dxj[k][2] * V[k][3]) / V[k][0];
    dtijdu1[1][1] = 2.0 * dmu[k][1] * dudxj[1][1] + dlambda[k][1] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][0] / V[k][0];
    dtijdu2[1][1] = 2.0 * dmu[k][2] * dudxj[1][1] + 2.0 * mu * dp1dxj[k][1] / V[k][0] + dlambda[k][2] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][1] / V[k][0];
    dtijdu3[1][1] = 2.0 * dmu[k][3] * dudxj[1][1] + dlambda[k][3] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][2] / V[k][0];
    dtijdu4[1][1] = 2.0 * dmu[k][4] * dudxj[1][1] + dlambda[k][4] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]);

    dtijdu0[2][2] = 2.0 * dmu[k][0] * dudxj[2][2] - 2.0 * mu * dp1dxj[k][2] * V[k][3] / V[k][0] + dlambda[k][0] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) - lambda * (dp1dxj[k][2] * V[k][3] + dp1dxj[k][0] * V[k][1] + dp1dxj[k][1] * V[k][2]) / V[k][0];
    dtijdu1[2][2] = 2.0 * dmu[k][1] * dudxj[2][2] + dlambda[k][1] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][0] / V[k][0];
    dtijdu2[2][2] = 2.0 * dmu[k][2] * dudxj[2][2] + dlambda[k][2] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][1] / V[k][0];
    dtijdu3[2][2] = 2.0 * dmu[k][3] * dudxj[2][2] + 2.0 * mu * dp1dxj[k][2] / V[k][0] + dlambda[k][3] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]) + lambda * dp1dxj[k][2] / V[k][0];
    dtijdu4[2][2] = 2.0 * dmu[k][4] * dudxj[2][2] + dlambda[k][4] * (dudxj[0][0] + dudxj[1][1] + dudxj[2][2]);

    dtijdu0[0][1] = - mu * (dp1dxj[k][0] * V[k][2] + dp1dxj[k][1] * V[k][1]) / V[k][0] + dmu[k][0] * (dudxj[1][0] + dudxj[0][1]);
    dtijdu1[0][1] = mu * dp1dxj[k][1] / V[k][0] + dmu[k][1] * (dudxj[1][0] + dudxj[0][1]);
    dtijdu2[0][1] = mu * dp1dxj[k][0] / V[k][0] + dmu[k][2] * (dudxj[1][0] + dudxj[0][1]);
    dtijdu3[0][1] = dmu[k][3] * (dudxj[1][0] + dudxj[0][1]);
    dtijdu4[0][1] = dmu[k][4] * (dudxj[1][0] + dudxj[0][1]);

    dtijdu0[0][2] = - mu * (dp1dxj[k][0] * V[k][3] + dp1dxj[k][2] * V[k][1]) / V[k][0] + dmu[k][0] * (dudxj[2][0] + dudxj[0][2]);
    dtijdu1[0][2] = mu * dp1dxj[k][2] / V[k][0] + dmu[k][1] * (dudxj[2][0] + dudxj[0][2]);
    dtijdu2[0][2] = dmu[k][2] * (dudxj[2][0] + dudxj[0][2]);
    dtijdu3[0][2] = mu * dp1dxj[k][0] / V[k][0] + dmu[k][3] * (dudxj[2][0] + dudxj[0][2]);
    dtijdu4[0][2] = dmu[k][4] * (dudxj[2][0] + dudxj[0][2]);

    dtijdu0[1][2] = - mu * (dp1dxj[k][1] * V[k][3] + dp1dxj[k][2] * V[k][2]) / V[k][0] + dmu[k][0] * (dudxj[2][1] + dudxj[1][2]);
    dtijdu1[1][2] = dmu[k][1] * (dudxj[2][1] + dudxj[1][2]);
    dtijdu2[1][2] = mu * dp1dxj[k][2] / V[k][0] + dmu[k][2] * (dudxj[2][1] + dudxj[1][2]);
    dtijdu3[1][2] = mu * dp1dxj[k][1] / V[k][0] + dmu[k][3] * (dudxj[2][1] + dudxj[1][2]);
    dtijdu4[1][2] = dmu[k][4] * (dudxj[2][1] + dudxj[1][2]);

    dtijdu0[1][0] = dtijdu0[0][1];
    dtijdu1[1][0] = dtijdu1[0][1];
    dtijdu2[1][0] = dtijdu2[0][1];
    dtijdu3[1][0] = dtijdu3[0][1];
    dtijdu4[1][0] = dtijdu4[0][1];

    dtijdu0[2][0] = dtijdu0[0][2];
    dtijdu1[2][0] = dtijdu1[0][2];
    dtijdu2[2][0] = dtijdu2[0][2];
    dtijdu3[2][0] = dtijdu3[0][2];
    dtijdu4[2][0] = dtijdu4[0][2];

    dtijdu0[2][1] = dtijdu0[1][2];
    dtijdu1[2][1] = dtijdu1[1][2];
    dtijdu2[2][1] = dtijdu2[1][2];
    dtijdu3[2][1] = dtijdu3[1][2];
    dtijdu4[2][1] = dtijdu4[1][2];

    dRdU[k][0][0] = 0.0;
    dRdU[k][0][1] = 0.0;
    dRdU[k][0][2] = 0.0;
    dRdU[k][0][3] = 0.0;
    dRdU[k][0][4] = 0.0;

    dRdU[k][1][0] = 0.0;
    dRdU[k][1][1] = 0.0;
    dRdU[k][1][2] = 0.0;
    dRdU[k][1][3] = 0.0;
    dRdU[k][1][4] = 0.0;

    dRdU[k][2][0] = 0.0;
    dRdU[k][2][1] = 0.0;
    dRdU[k][2][2] = 0.0;
    dRdU[k][2][3] = 0.0;
    dRdU[k][2][4] = 0.0;

    dRdU[k][3][0] = 0.0;
    dRdU[k][3][1] = 0.0;
    dRdU[k][3][2] = 0.0;
    dRdU[k][3][3] = 0.0;
    dRdU[k][3][4] = 0.0;

    dRdU[k][4][0] = (Vwall[1] * dtijdu0[0][0] + Vwall[2] * dtijdu0[1][0] + Vwall[3] * dtijdu0[2][0]) * n[0] +
                    (Vwall[1] * dtijdu0[0][1] + Vwall[2] * dtijdu0[1][1] + Vwall[3] * dtijdu0[2][1]) * n[1] + 
                    (Vwall[1] * dtijdu0[0][2] + Vwall[2] * dtijdu0[1][2] + Vwall[3] * dtijdu0[2][2]) * n[2];

    dRdU[k][4][1] = (Vwall[1] * dtijdu1[0][0] + Vwall[2] * dtijdu1[1][0] + Vwall[3] * dtijdu1[2][0]) * n[0] +
                    (Vwall[1] * dtijdu1[0][1] + Vwall[2] * dtijdu1[1][1] + Vwall[3] * dtijdu1[2][1]) * n[1] + 
                    (Vwall[1] * dtijdu1[0][2] + Vwall[2] * dtijdu1[1][2] + Vwall[3] * dtijdu1[2][2]) * n[2];

    dRdU[k][4][2] = (Vwall[1] * dtijdu2[0][0] + Vwall[2] * dtijdu2[1][0] + Vwall[3] * dtijdu2[2][0]) * n[0] +
                    (Vwall[1] * dtijdu2[0][1] + Vwall[2] * dtijdu2[1][1] + Vwall[3] * dtijdu2[2][1]) * n[1] + 
                    (Vwall[1] * dtijdu2[0][2] + Vwall[2] * dtijdu2[1][2] + Vwall[3] * dtijdu2[2][2]) * n[2];

    dRdU[k][4][3] = (Vwall[1] * dtijdu3[0][0] + Vwall[2] * dtijdu3[1][0] + Vwall[3] * dtijdu3[2][0]) * n[0] +
                    (Vwall[1] * dtijdu3[0][1] + Vwall[2] * dtijdu3[1][1] + Vwall[3] * dtijdu3[2][1]) * n[1] + 
                    (Vwall[1] * dtijdu3[0][2] + Vwall[2] * dtijdu3[1][2] + Vwall[3] * dtijdu3[2][2]) * n[2];

    dRdU[k][4][4] = (Vwall[1] * dtijdu4[0][0] + Vwall[2] * dtijdu4[1][0] + Vwall[3] * dtijdu4[2][0]) * n[0] +
                    (Vwall[1] * dtijdu4[0][1] + Vwall[2] * dtijdu4[1][1] + Vwall[3] * dtijdu4[2][1]) * n[1] +
                    (Vwall[1] * dtijdu4[0][2] + Vwall[2] * dtijdu4[1][2] + Vwall[3] * dtijdu4[2][2]) * n[2];

  }

}
//------------------------------------------------------------------------------

// Original
/*
template<int neq>
inline
void NavierStokesTerm::computeJacobianSurfaceTermNS(double dp1dxj[4][3], Vec3D &n, 
						    double *Vwall, double *V[4], 
						    double (*dRdU)[neq][neq])
{

  double T[4], Tcg;
  computeTemperature(V, T, Tcg);

  double mu = ooreynolds_mu * viscoFcn->compute_mu(Tcg);
  double lambda = ooreynolds_lambda * viscoFcn->compute_lambda(Tcg, mu);

  for (int k=0; k<4; ++k) {
    
    double rho = V[k][0];
    double u = V[k][1];
    double v = V[k][2];
    double w = V[k][3];

    double dp1dx = dp1dxj[k][0];
    double dp1dy = dp1dxj[k][1];
    double dp1dz = dp1dxj[k][2];

    double nu = mu / rho;
    double lu = lambda / rho;

//    double adtxxdu0 = twothird * nu * (-2.0 * dp1dx * u + dp1dy * v + dp1dz * w);
//    double adtxxdu1 = 2.0 * twothird * nu * dp1dx;
//    double adtxxdu2 = - twothird * nu * dp1dy;
//    double adtxxdu3 = - twothird * nu * dp1dz;

//    double adtyydu0 = twothird * nu * (-2.0 * dp1dy * v + dp1dx * u + dp1dz * w);
//    double adtyydu1 = - twothird * nu * dp1dx;
//    double adtyydu2 = 2.0 * twothird * nu * dp1dy;
//    double adtyydu3 = - twothird * nu * dp1dz;

//    double adtzzdu0 = twothird * nu * (-2.0 * dp1dz * w + dp1dx * u + dp1dy * v);
//    double adtzzdu1 = - twothird * nu * dp1dx;
//    double adtzzdu2 = - twothird * nu * dp1dy;
//    double adtzzdu3 = 2.0 * twothird * nu * dp1dz;

    double dtxxdu0 = -(2.0 * nu * dp1dx * u + lu * (dp1dx * u + dp1dy * v + dp1dz * w));
    double dtxxdu1 = (2.0 * nu + lu) * dp1dx;
    double dtxxdu2 = lu * dp1dy;
    double dtxxdu3 = lu * dp1dz;

    double dtyydu0 = -(2.0 * nu * dp1dy * v + lu * (dp1dx * u + dp1dy * v + dp1dz * w));
    double dtyydu1 = lu * dp1dx;
    double dtyydu2 = (2.0 * nu + lu) * dp1dy;
    double dtyydu3 = lu * dp1dz;

    double dtzzdu0 = -(2.0 * nu * dp1dz * w + lu * (dp1dx * u + dp1dy * v + dp1dz * w));
    double dtzzdu1 = lu * dp1dx;
    double dtzzdu2 = lu * dp1dy;
    double dtzzdu3 = (2.0 * nu + lu) * dp1dz;
    
//---check

//    if (adtxxdu0!=dtxxdu0) fprintf(stderr, "computeJacobianVolumeTermNS, dtxxdu0 = %e and adtxxdu0 = %f differ\n", dtxxdu0, adtxxdu0);
//    if (adtxxdu1!=dtxxdu1) fprintf(stderr, "computeJacobianVolumeTermNS, dtxxdu1 = %e and adtxxdu1 = %f differ\n", dtxxdu1, adtxxdu1);
//    if (adtxxdu2!=dtxxdu2) fprintf(stderr, "computeJacobianVolumeTermNS, dtxxdu2 = %e and adtxxdu2 = %f differ\n", dtxxdu2, adtxxdu2);
//    if (adtxxdu3!=dtxxdu3) fprintf(stderr, "computeJacobianVolumeTermNS, dtxxdu3 = %e and adtxxdu3 = %f differ\n", dtxxdu3, adtxxdu3);
                                                                                                                                                                                                     
//    if (adtyydu0!=dtyydu0) fprintf(stderr, "computeJacobianVolumeTermNS, dtyydu0 = %e and adtyydu0 = %f differ\n", dtyydu0, adtyydu0);
//    if (adtyydu1!=dtyydu1) fprintf(stderr, "computeJacobianVolumeTermNS, dtyydu1 = %e and adtyydu1 = %f differ\n", dtyydu1, adtyydu1);
//    if (adtyydu2!=dtyydu2) fprintf(stderr, "computeJacobianVolumeTermNS, dtyydu2 = %e and adtyydu2 = %f differ\n", dtyydu2, adtyydu2);
//    if (adtyydu3!=dtyydu3) fprintf(stderr, "computeJacobianVolumeTermNS, dtyydu3 = %e and adtyydu3 = %f differ\n", dtyydu3, adtyydu3);
                                                                                                                                                                                                     
//    if (adtzzdu0!=dtzzdu0) fprintf(stderr, "computeJacobianVolumeTermNS, dtzzdu0 = %e and adtzzdu0 = %f differ\n", dtzzdu0, adtzzdu0);
//    if (adtzzdu1!=dtzzdu1) fprintf(stderr, "computeJacobianVolumeTermNS, dtzzdu1 = %e and adtzzdu1 = %f differ\n", dtzzdu1, adtzzdu1);
//    if (adtzzdu2!=dtzzdu2) fprintf(stderr, "computeJacobianVolumeTermNS, dtzzdu2 = %e and adtzzdu2 = %f differ\n", dtzzdu2, adtzzdu2);
//    if (adtzzdu3!=dtzzdu3) fprintf(stderr, "computeJacobianVolumeTermNS, dtzzdu3 = %e and adtzzdu3 = %f differ\n", dtzzdu3, adtzzdu3);

    double dtxydu0 = - nu * (dp1dx * v + dp1dy * u);
    double dtxydu1 = nu * dp1dy;
    double dtxydu2 = nu * dp1dx;

    double dtxzdu0 = -nu * (dp1dx * w + dp1dz * u);
    double dtxzdu1 = nu * dp1dz;
    double dtxzdu3 = nu * dp1dx;

    double dtyzdu0 = -nu * (dp1dy * w + dp1dz * v);
    double dtyzdu2 = nu * dp1dz;
    double dtyzdu3 = nu * dp1dy;

    dRdU[k][0][0] = 0.0;
    dRdU[k][0][1] = 0.0;
    dRdU[k][0][2] = 0.0;
    dRdU[k][0][3] = 0.0;
    dRdU[k][0][4] = 0.0;

    dRdU[k][1][0] = 0.0;
    dRdU[k][1][1] = 0.0;
    dRdU[k][1][2] = 0.0;
    dRdU[k][1][3] = 0.0;
    dRdU[k][1][4] = 0.0;

    dRdU[k][2][0] = 0.0;
    dRdU[k][2][1] = 0.0;
    dRdU[k][2][2] = 0.0;
    dRdU[k][2][3] = 0.0;
    dRdU[k][2][4] = 0.0;

    dRdU[k][3][0] = 0.0;
    dRdU[k][3][1] = 0.0;
    dRdU[k][3][2] = 0.0;
    dRdU[k][3][3] = 0.0;
    dRdU[k][3][4] = 0.0;

    dRdU[k][4][0] = (Vwall[1]*dtxxdu0 + Vwall[2]*dtxydu0 + Vwall[3]*dtxzdu0) * n[0] +
      (Vwall[1]*dtxydu0 + Vwall[2]*dtyydu0 + Vwall[3]*dtyzdu0) * n[1] + 
      (Vwall[1]*dtxzdu0 + Vwall[2]*dtyzdu0 + Vwall[3]*dtzzdu0) * n[2];

    dRdU[k][4][1] = (Vwall[1]*dtxxdu1 + Vwall[2]*dtxydu1 + Vwall[3]*dtxzdu1) * n[0] +
      (Vwall[1]*dtxydu1 + Vwall[2]*dtyydu1) * n[1] + 
      (Vwall[1]*dtxzdu1 + Vwall[3]*dtzzdu1) * n[2];

    dRdU[k][4][2] = (Vwall[1]*dtxxdu2 + Vwall[2]*dtxydu2) * n[0] +
      (Vwall[1]*dtxydu2 + Vwall[2]*dtyydu2 + Vwall[3]*dtyzdu2) * n[1] + 
      (Vwall[2]*dtyzdu2 + Vwall[3]*dtzzdu2) * n[2];

    dRdU[k][4][3] = (Vwall[1]*dtxxdu3 + Vwall[3]*dtxzdu3) * n[0] +
      (Vwall[2]*dtyydu3 + Vwall[3]*dtyzdu3) * n[1] + 
      (Vwall[1]*dtxzdu3 + Vwall[2]*dtyzdu3 + Vwall[3]*dtzzdu3) * n[2];

    dRdU[k][4][4] = 0.0;

  }

}
*/

//------------------------------------------------------------------------------

#endif
