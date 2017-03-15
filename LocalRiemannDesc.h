#ifndef _LOCAL_RIEMANN_DESC_H
#define _LOCAL_RIEMANN_DESC_H

#include "LinkF77.h"
#include <LocalRiemann.h>
#include <VarFcn.h>
#include "IoData.h"
#include "SparseGridCluster.h"
#include <cmath>
#include "ImplicitRiemann.h"
#include "DenseMatrixOps.h"
#include "DebugTools.h"
//----------------------------------------------------------------------------
// First the derived classes of LocalRiemannGfmp (no exact Riemann problem)
// Second the derived classes of LocalRiemannGfmpar (with exact Riemann prob)
// Third the derived classes of LocalRiemannFluidStructure
//----------------------------------------------------------------------------

extern "C" {
  void F77NAME(eriemanngw) (const double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&, const double&, const double &,
                            const double &, const double&, const double&,
                            const double &, const double&, const double&,
                            const int&);
  void F77NAME(eriemanngg) (const double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&, const double&, const double &,
                            const double &, const double&, const double&,
                            const double &, const double&,int &,
                            const double&, const double&, const double&,
                            const double&, const double&, const int&);
  void F77NAME(eriemannww) (const double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&,  int &,
                            const double&, const double&, const double&,
                            const double&, const double&, const int&);
};

//----------------------------------------------------------------------------

class LocalRiemannGfmpGasGas : public LocalRiemannGfmp {

public:
  LocalRiemannGfmpGasGas(VarFcn *vf, int tag1, int tag2) : LocalRiemannGfmp(vf,tag1,tag2) {}
  ~LocalRiemannGfmpGasGas() { vf_ = 0; }

  int computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *initWi, double *initWj,
                              double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it,bool isHigherOrder);

  void computeRiemannJacobian(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *Wi, double *Wj,
                              double dx[3],int it,
                              double* dWidWi,double*  dWidWj,
                              double* dWjdWi,double*  dWjdWj) {
    fprintf(stderr,"ERROR: computeRiemannJacobian is not implemeted in LocalRiemannGfmpGasGas!\n");}

 // FS Riemann problem (implemented here just to stop compiler's complaining.)
  int computeRiemannSolution(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpGasGas!\n"); return 0; } 

  void computeRiemannJacobian(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, double* dWstardU,int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpGasGas!\n");}

  void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpGasGas!\n");}

private:
  LocalRiemannGfmpGasGas();
};

//----------------------------------------------------------------------------

inline
int LocalRiemannGfmpGasGas::computeRiemannSolution(double *Vi, double *Vj,
    int IDi, int IDj, double *nphi,
    double *initWi, double *initWj,
    double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj,
						    double dx[3], int it,bool isHigherOrder)
{

  for (int i=0; i<10; i++){
    Wi[i] = Vj[i];
    Wj[i] = Vi[i];
  }

  return 0;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmpTaitTait : public LocalRiemannGfmp {

public:
  LocalRiemannGfmpTaitTait(VarFcn *vf, int tag1, int tag2) : LocalRiemannGfmp(vf,tag1,tag2) {}
  ~LocalRiemannGfmpTaitTait() { vf_ = 0; }

int computeRiemannSolution(double *Vi, double *Vj,
                            int IDi, int IDj, double *nphi,
                            double *initWi, double *initWj,
                            double *Wi, double *Wj,
                            double *rupdatei, double *rupdatej, 
                            double &weighti, double &weightj,
                            double dx[3], int it,bool isHigherOrder);

  void computeRiemannJacobian(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *Wi, double *Wj,
                              double dx[3],int it,
                              double* dWidWi,double*  dWidWj,
                              double* dWjdWi,double*  dWjdWj) {
    fprintf(stderr,"ERROR: computeRiemannJacobian is not implemeted in LocalRiemannGfmpTaitTait!\n");}

 // FS Riemann problem (implemented here just to stop compiler's complaining.)
  int computeRiemannSolution(double *Vi, double *Vstar, 
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpTaitTait!\n"); return 0; }

  void computeRiemannJacobian(double *Vi, double *Vstar, 
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, double* dWstardU,int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpTaitTait!\n");}

  void computeRiemannSolution(int tag, double *Vi, double *Vstar, 
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpTaitTait!\n");}

private:
  LocalRiemannGfmpTaitTait();
};

//----------------------------------------------------------------------------

inline
int LocalRiemannGfmpTaitTait::computeRiemannSolution(double *Vi, double *Vj,
    int IDi, int IDj, double *nphi,
    double *initWi, double *initWj,
    double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj,
    double dx[3], int it,bool isHigherOrder)
{

  for (int i=0; i<10; i++){
    Wi[i] = Vj[i];
    Wj[i] = Vi[i];
  }

  double a1 = vf_->getAlphaWater(this->fluid1);
  double b1 = vf_->getBetaWater(this->fluid1);
  double p1 = vf_->getPrefWater(this->fluid1);
  double a2 = vf_->getAlphaWater(this->fluid2);
  double b2 = vf_->getBetaWater(this->fluid2);
  double p2 = vf_->getPrefWater(this->fluid2);

  double temp = 0;
	
  if(IDi==fluid1){
    temp = p2+a2*pow(Vj[0],b2);
    Wi[0] = pow((temp-p1)/a1,1.0/b1);
    temp = p2+a2*pow(Vj[5],b2);
    Wi[5] = pow((temp-p1)/a1,1.0/b1);
    temp = p1+a1*pow(Vi[0],b1);
    Wj[0] = pow((temp-p2)/a2,1.0/b2);
    temp = p1+a1*pow(Vi[5],b1);
    Wj[5] = pow((temp-p2)/a2,1.0/b2);
  }else{
    temp = p2+a2*pow(Vi[0],b2);
    Wj[0] = pow((temp-p1)/a1,1.0/b1);
    temp = p2+a2*pow(Vi[5],b2);
    Wj[5] = pow((temp-p1)/a1,1.0/b1);
    temp = p1+a1*pow(Vj[0],b1);
    Wi[0] = pow((temp-p2)/a2,1.0/b2);
    temp = p1+a1*pow(Vj[5],b1);
    Wi[5] = pow((temp-p2)/a2,1.0/b2);
  }

  return 0;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmpJWLJWL : public LocalRiemannGfmp {

public:
  LocalRiemannGfmpJWLJWL(VarFcn *vf, int tag1, int tag2) : LocalRiemannGfmp(vf,tag1,tag2) {}
  ~LocalRiemannGfmpJWLJWL() { vf_ = 0; }

  int computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *initWi, double *initWj,
                              double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it,bool isHigherOrder);

  void computeRiemannJacobian(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *Wi, double *Wj,
                              double dx[3],int it,
                              double* dWidWi,double*  dWidWj,
                              double* dWjdWi,double*  dWjdWj) {
    fprintf(stderr,"ERROR: computeRiemannJacobian is not implemeted in LocalRiemannGfmpJWLJWL!\n");}

 // FS Riemann problem (implemented here just to stop compiler's complaining.)
  int computeRiemannSolution(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpJWLJWL!\n"); return 0; }

  void computeRiemannJacobian(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, double* dWstardU,int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpJWLJWL!\n");}

  void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpJWLJWL!\n");}

private:
  LocalRiemannGfmpJWLJWL();
};

//----------------------------------------------------------------------------

inline
int LocalRiemannGfmpJWLJWL::computeRiemannSolution(double *Vi, double *Vj,
    int IDi, int IDj, double *nphi,
    double *initWi, double *initWj,
    double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj,
    double dx[3], int it,bool isHigherOrder)
{

  for (int i=0; i<10; i++){
    Wi[i] = Vj[i];
    Wj[i] = Vi[i];
  }

  return 0;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmpGasJWL : public LocalRiemannGfmp {

public:
  LocalRiemannGfmpGasJWL(VarFcn *vf, int tag1, int tag2) : LocalRiemannGfmp(vf,tag1,tag2) {}
  ~LocalRiemannGfmpGasJWL() { vf_ = 0; }

  int computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *initWi, double *initWj,
                              double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it,bool isHigherOrder);

  void computeRiemannJacobian(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *Wi, double *Wj,
                              double dx[3],int it,
                              double* dWidWi,double*  dWidWj,
                              double* dWjdWi,double*  dWjdWj) {
    fprintf(stderr,"ERROR: computeRiemannJacobian is not implemeted in LocalRiemannGfmpGasJWL!\n");}

 // FS Riemann problem (implemented here just to stop compiler's complaining.)
  int computeRiemannSolution(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpGasJWL!\n"); return 0; }

  void computeRiemannJacobian(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, double* dWstardU,int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpGasJWL!\n");}

  void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmpGasJWL!\n");}

private:
  LocalRiemannGfmpGasJWL();
};

//----------------------------------------------------------------------------

inline
int LocalRiemannGfmpGasJWL::computeRiemannSolution(double *Vi, double *Vj,
						    int IDi, int IDj, double *nphi,
						    double *initWi, double *initWj,
						    double *Wi, double *Wj,
						    double *rupdatei, double *rupdatej, double &weighti, double &weightj,
						    double dx[3], int it,bool isHigherOrder)
{

  for (int i=0; i<10; i++){
    Wi[i] = Vj[i];
    Wj[i] = Vi[i];
  }

  return 0;
  //Fedkiw's GFM -- extrapolation of density in the ghost cell
  //fprintf(stderr, "using density to extrapolate in ghost cell\n");
  //Wi[0] = Vi[0]; Wi[5] = Vi[5];
  //Wj[0] = Vj[0]; Wj[5] = Vj[5];

}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmparGasGas : public LocalRiemannGfmpar {

  double tolpre_;
  int nriter_;

public:
  LocalRiemannGfmparGasGas(VarFcn *vf, int tag1, int tag2, MultiFluidData::TypePhaseChange typePhaseChange,
                           double tolpre, int nriter) : LocalRiemannGfmpar(vf,tag1,tag2,typePhaseChange)
                           { tolpre_ = tolpre; nriter_ = nriter; }
  ~LocalRiemannGfmparGasGas() { vf_ = 0; }

  int computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *initWi, double *initWj,
                              double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it,bool isHigherOrder);

  void computeRiemannJacobian(double *Vi, double *Vj,
			      int IDi, int IDj, double *nphi,
			      double *Wi, double *Wj,
			      double dx[3],int it,
			      double* dWidUi,double*  dWidUj,double* dWjdUi,double*  dWjdUj);
/*
  int eriemann(double rhol, double ul, double pl, 
                double rhor, double ur, double pr, 
                double &pi, double &ui, double &rhoil, double &rhoir){

  int err;
  F77NAME(eriemanngg)(  rhol,ul,pl,rhor,ur,pr,pi,ui,rhoil,rhoir,
                        vf_->getGamma(fluid2), vf_->getPressureConstant(fluid2), 
                        vf_->getGamma(fluid1), vf_->getPressureConstant(fluid1), err);
  return err; 
  }
*/
 // FS Riemann problem (implemented here just to stop compiler's complaining.)
  int computeRiemannSolution(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparGasGas!\n"); return 0; }

  void computeRiemannJacobian(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, double* dWstardU,int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparGasGas!\n");}

  void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparGasGas!\n");}

private:
  LocalRiemannGfmparGasGas();
};

//----------------------------------------------------------------------------

inline
int LocalRiemannGfmparGasGas::computeRiemannSolution(double *Vi, double *Vj,
	 	int IDi, int IDj, double *nphi,
                double *initWi, double *initWj,
		double *Wi, double *Wj,
                double *rupdatei, double *rupdatej, 
                double &weighti, double &weightj,
                double dx[3], int it,bool isHigherOrder)
{
  int dim = 5;
	
  double P_1, P_2, U_1, U_2, R_1, R_2;
  double P_i, U_i, R_i1, R_i2;

  double gam1  = vf_->getGamma(fluid1);
  double pref1 = vf_->getPressureConstant(fluid1);
  double gam2  = vf_->getGamma(fluid2);
  double pref2 = vf_->getPressureConstant(fluid2);

  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  double pmin1 = vf_->getVarFcnBase(fluid1)->pmin;
  double pmin2 = vf_->getVarFcnBase(fluid2)->pmin;
  double rhomin1 = vf_->getVarFcnBase(fluid1)->rhomin;
  double rhomin2 = vf_->getVarFcnBase(fluid2)->rhomin;

  double vmid[3] = {0.5*(vti[0]+vtj[0]),
		    0.5*(vti[1]+vtj[1]),
		    0.5*(vti[2]+vtj[2])};
		    

  int err;
  if (IDi==fluid1) {

    // cell i is fluid1
    // cell j is fluid2
    R_2  = Vj[0];     R_1 = Vi[0];
    U_2  = vnj;       U_1 = vni;
    P_2  = vf_->getPressure(Vj, IDj);
    P_1  = vf_->getPressure(Vi, IDi);


    F77NAME(eriemanngg)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,gam2,pref2,gam1,pref1,err,
                        pmin2,pmin1,rhomin2,rhomin1,tolpre_,nriter_);

    Wi[0]  = R_i1;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = P_i;                     Wi[dim+4]  = Wi[4];

    Wj[0]  = R_i2;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    Wj[4]  = P_i;                     Wj[dim+4]  = Wj[4];

  }else{
    // cell i is fluid2
    // cell j is fluid1
    R_2  = Vi[0];     R_1  = Vj[0];
    U_2  = vni;       U_1  = vnj;
    P_2  = vf_->getPressure(Vi, IDi);
    P_1  = vf_->getPressure(Vj, IDj);

    F77NAME(eriemanngg)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,gam2,pref2,gam1,pref1,err,
                        pmin2,pmin1,rhomin2,rhomin1,tolpre_,nriter_);

    Wi[0]  = R_i2;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = P_i;                     Wi[dim+4]  = Wi[4];

    Wj[0]  = R_i1;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    Wj[4]  = P_i;                     Wj[dim+4]  = Wj[4];
  }

// to update the nodes when they change fluids
// METHOD1: naive approach of averaging the Riemann solution
  /*if(it==1){
    weighti += 1.0;
    weightj += 1.0;
    for (int k=0; k<5; k++){
      rupdatei[k] += Wj[k];
      rupdatej[k] += Wi[k];
    }
  }*/

// METHOD 2 : combine averaging and direction of flow
  if (it == 1 && !isHigherOrder)
    updatePhaseChangingNodeValues(dx, Wi, Wj, weighti, rupdatei, weightj, rupdatej);

  return err;
}

inline 
void LocalRiemannGfmparGasGas::computeRiemannJacobian(double *Vi, double *Vj,
				                      int IDi, int IDj, double *nphi,
				                      double *Wi, double *Wj,
				                      double dx[3],int it,
						      double* dWidWi,double*  dWidWj,
						      double* dWjdWi,double*  dWjdWj) {
  
  int dim = 5;
  int k,l;
  
  double P_1, P_2, U_1, U_2, R_1, R_2;
  double P_i, U_i, R_i1, R_i2;

  double gam1  = vf_->getGamma(fluid1);
  double pref1 = vf_->getPressureConstant(fluid1);
  double gam2  = vf_->getGamma(fluid2);
  double pref2 = vf_->getPressureConstant(fluid2);

  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  // 3x3 Jacobians, directly from implicit riemann jacobian
  double dWidWi3[9],  dWidWj3[9], dWjdWj3[9], dWjdWi3[9];

  if (IDi==fluid1) {

    // cell i is fluid1
    // cell j is fluid2
    R_2  = Vj[0];     R_1 = Vi[0];
    U_2  = vnj;       U_1 = vni;
    P_2  = vf_->getPressure(Vj, IDj);
    P_1  = vf_->getPressure(Vi, IDi);

    ImplicitRiemann::computeGasGasJacobian(Wi[4], gam2,pref2,P_2,R_2, gam1, pref1, P_1,R_1, dWjdWj3, dWjdWi3,  dWidWi3, dWidWj3 );
  }else{
    // cell i is fluid2
    // cell j is fluid1
    R_2  = Vi[0];     R_1  = Vj[0];
    U_2  = vni;       U_1  = vnj;
    P_2  = vf_->getPressure(Vi, IDi);
    P_1  = vf_->getPressure(Vj, IDj);

    ImplicitRiemann::computeGasGasJacobian(Wi[4], gam2, pref2, P_2,R_2,gam1,pref1,P_1,R_1, dWidWi3, dWidWj3,  dWjdWj3, dWjdWi3 );
  }

  this->oneDtoThreeD(dWidWi3, dWidWj3,
		     dWjdWi3,dWjdWj3,nphi,
		     dWidWi, dWidWj,
		     dWjdWi, dWjdWj);
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmparGasTait: public LocalRiemannGfmpar {

public:
  LocalRiemannGfmparGasTait(VarFcn *vf, int tag1, int tag2, MultiFluidData::TypePhaseChange typePhaseChange) : LocalRiemannGfmpar(vf,tag1,tag2,typePhaseChange) {}
  ~LocalRiemannGfmparGasTait() { vf_ = 0; }

  int computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *initWi, double *initWj,
                              double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj, 
                              double dx[3], int it,bool isHigherOrder);

  void computeRiemannJacobian(double *Vi, double *Vj,
			      int IDi, int IDj, double *nphi,
			      double *Wi, double *Wj,
			      double dx[3],int it,
			      double* dWidUi,double*  dWidUj,double* dWjdUi,double*  dWjdUj);

 // FS Riemann problem (implemented here just to stop compiler's complaining.)
  int computeRiemannSolution(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparGasTait!\n"); return 0; }

  void computeRiemannJacobian(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, double* dWstardU,int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparGasTait!\n");}

  void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparGasTait!\n");}

private:
  LocalRiemannGfmparGasTait();
};

//----------------------------------------------------------------------------

inline void solveSGTait(double Rg,double Ug,double Pg, 
			double Rw,double Uw,double Pw,
			double &Pi,double &Ui, 
			double &Rig, double &Riw,
			double alpha,double beta,
			double pref, double gamma,
			double Pinf, double Pcg,
                        double rhocg, double Pcw,
                        double rhocw, int &ierr) {

  double Q,f,m,n,dQ,df,g,dg,db;
  double ag = sqrt(gamma/Rg*(Pg+Pinf));
  double aw = sqrt(alpha*beta*pow(Rw,beta-1.0)),b;
  Pi = sqrt(Pw*Pg);
  Riw = pow((Pi-pref)/alpha,1.0/beta);
  double dpdrho = alpha*beta*pow(Riw,beta-1.0);

  double Pmin = std::max<double>(  Pcg,   Pcw);
  double Rmin = std::max<double>(rhocg, rhocw);

  const int max_ite = 1000;

  ierr = 0;

  int k = 0;
  while (++k < max_ite) {

    // Gas relations
    if (Pi > Pg) {
      // Shock
      m = 2.0/((gamma+1.0)*Rg);
      n = (Pg+Pinf)*(gamma-1.0)/(gamma+1.0);
      Q = sqrt( (Pi+Pinf+n)/m );
      dQ = 0.5/(Q*m);
      f = (Pi-Pg)/Q;
      df = 1.0/Q - (Pi-Pg)/(Q*Q)*dQ;
    } else {

      // Rarefaction
      f = 2.0*ag/(gamma-1.0)*( pow((Pi+Pinf)/(Pg+Pinf), (gamma-1.0)/(2.0*gamma) )-1.0 );
      df = 2.0*ag/(gamma-1.0)*((gamma-1.0)/(2.0*gamma)/(Pg+Pinf)*pow((Pi+Pinf)/(Pg+Pinf), (gamma-1.0)/(2.0*gamma)-1.0 ));
    }

    // Liquid relations
    if (Pi > Pw) {

      g = max(1.0e-8,sqrt( alpha*(pow(Riw,beta)-pow(Rw,beta))*(Riw-Rw)/(Riw*Rw) ));
      dg = 0.5/g*(alpha*( beta*pow(Riw,beta-1.0)*(Riw-Rw)/(Riw*Rw) + 
		  (pow(Riw,beta)-pow(Rw,beta))/(Riw*Rw) - 
		  (pow(Riw,beta)-pow(Rw,beta))*(Riw-Rw)/(Riw*Riw*Rw) )); // PJSA
    } else {

      g = 2.0*aw/(beta-1.0)*( pow(Riw/Rw, (beta-1.0)*0.5) - 1.0);
      dg = 2.0*aw/(beta-1.0)*(beta-1.0)*0.5/Rw*pow(Riw/Rw, (beta-1.0)*0.5-1.0);
    }

    b = f+g-Uw+Ug;
    db = df+dg/dpdrho;   
    if (fabs(b/db) < 1.0e-6*Pi || fabs(b) < 1.0e-8)
      break;
    Pi -= b/db;
    if (Pi < 1.0e-10)
      Pi = 1.0e-10;

    if (Pi < Pmin)
      Pi = Pmin;

    Riw = pow((Pi-pref)/alpha,1.0/beta);

    if (Riw < rhocw)
      Riw = rhocw;

    dpdrho = alpha*beta*pow(Riw,beta-1.0);    

  }

  if (k >= max_ite) {
    std::cout << "*** Warning " << std::endl;
    std::cout << "Newton for gas-tait ERS reached max num. iterations " << max_ite << std::endl;
    std::cout << "without converging to the desired tolerance " << 1.0e-6 << std::endl;
    std::cout << "Input gas:  Rg = " << Rg << " Ug = " << Ug << " Pg = " << Pg << std::endl;
    std::cout << "Input Tait: Rw = " << Rw << " Uw = " << Uw << " Pw = " << Pw << std::endl;
    std::cout << "Output    : Pi = " << Pi << " Riw = " << Riw << " b = " << b << " db = " << db << std::endl;    
    std::cout << "*** " << std::endl;
    ierr = 1;
  }

  Ui = (0.5*(Uw+Ug)+0.5*(f-g));
  if (Pi > Pg) {
    double h = (gamma-1.0)/(gamma+1.0);
    double j = (Pi+Pinf)/(Pg+Pinf);
    Rig = Rg*(j+h)/(j*h+1.0);
  } else {
    Rig = Rg*pow( (Pi+Pinf)/(Pg+Pinf), 1.0/gamma);
  }
  Rig = std::max<double>(rhocg, Rig);

  /*if( Pi <= Pmin || Rig <= Rmin || Riw <= Rmin ){
    std::cout << "*** ERROR ERS SGTait: detected too small density or pressure " << std::endl;
    std::cout << " Rig, Riw, Pi " << std::endl;
    std::cout << Rig << " " << Riw << " " << Pi << std::endl;
    ierr = 1;
  }*/

}

inline
int LocalRiemannGfmparGasTait::computeRiemannSolution(double *Vi, double *Vj,
	  int IDi, int IDj, double *nphi,
          double *initWi, double *initWj,
	  double *Wi, double *Wj,
          double *rupdatei, double *rupdatej, double &weighti, double &weightj,
          double dx[3], int it,bool isHigherOrder)
{

  int dim = 5;

  double alpha   = vf_->getAlphaWater(fluid2);
  double beta    = vf_->getBetaWater(fluid2);
  double pref    = vf_->getPrefWater(fluid2);
  double gam     = vf_->getGamma(fluid1);
  double Pinf    = vf_->getPressureConstant(fluid1);

  double cp = vf_->specificHeatCstPressure(fluid2);

  double T_w, P_g, P_w, U_w, U_g, R_w, R_g;
  double P_i, U_i, R_il, R_ir;

  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  int err = 0;

  if (IDi==fluid2) {
    if(IDj!=fluid1) {
      fprintf(stderr,"ERROR: IDi = %d, IDj = %d, fluid1 = %d, fluid2 = %d.\n", IDi, IDj, fluid1, fluid2);
      exit(-1);
    }
    // cell j is gas
    // cell i is tait
    R_g  = Vj[0];     R_w  = Vi[0];
    U_g  = vnj;       U_w  = vni;
    P_g  = vf_->getPressure(Vj, IDj);
    P_w  = vf_->getPressure(Vi, IDi);


    //F77NAME(eriemanngw)(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_il,R_ir,alpha,beta,pref,gam);
    solveSGTait(R_g,U_g,P_g, 
		R_w,U_w,P_w,
		P_i,U_i, 
		R_il,R_ir,
		alpha,beta,
		pref,gam,
		Pinf, vf_->getVarFcnBase(IDj)->pmin,
                vf_->getVarFcnBase(IDj)->rhomin, vf_->getVarFcnBase(IDi)->pmin,
		vf_->getVarFcnBase(IDi)->rhomin, err);

    Wi[0]  = R_ir;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    if (vf_->isBurnable(IDi))
      Wi[4]  = Vi[4] + 1.0/cp*(-0.5*(P_i+P_w)*(1.0/R_ir-1.0/R_w));
    else
      Wi[4]  = Vi[4] + 1.0/cp*(P_i/R_ir-P_w/R_w - 0.5*(P_i+P_w)*(1.0/R_ir-1.0/R_w));

    // CHF: believe this line is missing
    Wi[dim+4] = Wi[4];

    Wj[0]  = R_il;                      Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];        Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];        Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];        Wj[dim+3]  = Wj[3];
    Wj[4]  = P_i;                       Wj[dim+4]  = P_i;

  }else{
    if(IDi!=fluid1 || IDj!=fluid2) {
      fprintf(stderr,"ERROR: IDi = %d, IDj = %d, fluid1 = %d, fluid2 = %d.\n", IDi, IDj, fluid1, fluid2);
      exit(-1);
    }
    // cell j is tait
    // cell i is gas
    R_g  = Vi[0];     R_w  = Vj[0];
    U_g  = vni;       U_w  = vnj;
    P_g  = vf_->getPressure(Vi, IDi);
    P_w  = vf_->getPressure(Vj, IDj);

    //F77NAME(eriemanngw)(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_il,R_ir,alpha,beta,pref,gam);
    solveSGTait(R_g,U_g,P_g, 
		R_w,U_w,P_w,
		P_i,U_i, 
		R_il,R_ir,
		alpha,beta,
		pref,gam,
		Pinf , vf_->getVarFcnBase(IDi)->pmin,
                vf_->getVarFcnBase(IDi)->rhomin, vf_->getVarFcnBase(IDj)->pmin,
		vf_->getVarFcnBase(IDj)->rhomin, err);

    //std::cout << "P_i = " << P_i << " " << R_il << " " << R_ir << std::endl;

    Wi[0]  = R_il;                      Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];        Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];        Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];        Wi[dim+3]  = Wi[3];
    Wi[4]  = P_i;                       Wi[dim+4]  = P_i;

    Wj[0]  = R_ir;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    if (vf_->isBurnable(IDj))
      Wj[4]  = Vj[4] + 1.0/cp*(-0.5*(P_i+P_w)*(1.0/R_ir-1.0/R_w)); // PJSA
    else
      Wj[4]  = Vj[4] + 1.0/cp*(P_i/R_ir-P_w/R_w - 0.5*(P_i+P_w)*(1.0/R_ir-1.0/R_w)); // PJSA
    /*vf_->computeTemperature(Wj, IDi);*/      Wj[dim+4]  = Wj[4];

  }
  
// to update the nodes when they change fluids
// METHOD1: naive approach of averaging the Riemann solution
  /*if(it==1){
    weighti += 1.0;
    weightj += 1.0;
    for (int k=0; k<5; k++){
      rupdatei[k] += Wj[k];
      rupdatej[k] += Wi[k];
    }
  }*/

// METHOD 2 : combine averaging and direction of flow
  if (it == 1 && !isHigherOrder)
    updatePhaseChangingNodeValues(dx, Wi, Wj, weighti, rupdatei, weightj, rupdatej);

  return err;
}
//----------------------------------------------------------------------------

inline 
void LocalRiemannGfmparGasTait::computeRiemannJacobian(double *Vi, double *Vj,
						       int IDi, int IDj, double *nphi,
						       double *Wi, double *Wj,
						       double dx[3],int it,
						       double* dWidWi,double*  dWidWj,
						       double* dWjdWi,double*  dWjdWj) {
  int dim = 5;

  double alpha   = vf_->getAlphaWater(fluid2);
  double beta    = vf_->getBetaWater(fluid2);
  double pref    = vf_->getPrefWater(fluid2);
  double cp      = vf_->specificHeatCstPressure(fluid2);
  bool   isBurnable = vf_->isBurnable(fluid2);
  double gam     = vf_->getGamma(fluid1);
  double Pinf    = vf_->getPressureConstant(fluid1);

  double T_w, P_g, P_w, U_w, U_g, R_w, R_g;
  double P_i, U_i, R_il, R_ir,R_i1,R_i2;

  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  // 3x3 Jacobians, directly from implicit riemann jacobian
  double dWidWi3[9],  dWidWj3[9], dWjdWj3[9], dWjdWi3[9];
  double dTdrho,dTdp;

  if (IDi==fluid2) {
    // cell j is gas
    // cell i is tait
    R_g  = Vj[0];     R_w  = Vi[0];
    U_g  = vnj;       U_w  = vni;
    P_g  = vf_->getPressure(Vj, IDj);
    P_w  = vf_->getPressure(Vi, IDi);

    //F77NAME(eriemanngw)(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_il,R_ir,alpha,beta,pref,gam);

    ImplicitRiemann::computeGasTaitJacobian(Wj[4], gam, Pinf, P_g, R_g, alpha,
					    beta, pref, P_w, R_w, 
					    dWjdWj3, dWjdWi3, dWidWi3, dWidWj3,
                                            cp, isBurnable); // PJSA

  }else{
    // cell j is tait
    // cell i is gas
    R_g  = Vi[0];     R_w  = Vj[0];
    U_g  = vni;       U_w  = vnj;
    P_g  = vf_->getPressure(Vi, IDi);
    P_w  = vf_->getPressure(Vj, IDj);

    //F77NAME(eriemanngw)(R_g,U_g,P_g,R_w,U_w,P_w,P_i,U_i,R_il,R_ir,alpha,beta,pref,gam);

    ImplicitRiemann::computeGasTaitJacobian(Wi[4], gam, Pinf, P_g, R_g, alpha, // PJSA
					    beta, pref, P_w, R_w, 
					    dWidWi3, dWidWj3, dWjdWj3, dWjdWi3,
                                            cp, isBurnable); // PJSA
  }

  this->oneDtoThreeD(dWidWi3, dWidWj3,
		     dWjdWi3,dWjdWj3,nphi,
		     dWidWi, dWidWj,
		     dWjdWi, dWjdWj);
}

//----------------------------------------------------------------------------

class LocalRiemannGfmparTaitTait: public LocalRiemannGfmpar {

  double tolpre_;
  int nriter_;

public:
  LocalRiemannGfmparTaitTait(VarFcn *vf, int tag1, int tag2, MultiFluidData::TypePhaseChange typePhaseChange,
                             double tolpre, int nriter) : LocalRiemannGfmpar(vf,tag1,tag2,typePhaseChange)
                            { tolpre_ = tolpre; nriter_ = nriter; }
  ~LocalRiemannGfmparTaitTait() { vf_ = 0; }

  int computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *initWi, double *initWj,
                              double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it,bool isHigherOrder);

  void computeRiemannJacobian(double *Vi, double *Vj,
			      int IDi, int IDj, double *nphi,
			      double *Wi, double *Wj,
			      double dx[3],int it,
			      double* dWidUi,double*  dWidUj,double* dWjdUi,double*  dWjdUj);

 // FS Riemann problem (implemented here just to stop compiler's complaining.)
  int computeRiemannSolution(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparTaitTait!\n"); return 0; }

  void computeRiemannJacobian(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, double* dWstardU,int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparTaitTait!\n");}

  void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparTaitTait!\n");}

private:
  LocalRiemannGfmparTaitTait();
};

//----------------------------------------------------------------------------

inline
int LocalRiemannGfmparTaitTait::computeRiemannSolution(double *Vi, double *Vj,
    int IDi, int IDj, double *nphi,
    double *initWi, double *initWj,
    double *Wi, double *Wj,
    double *rupdatei, double *rupdatej, double &weighti, double &weightj,
    double dx[3], int it,bool isHigherOrder)
{

  int dim = 5;

  double alpha1   = vf_->getAlphaWater(fluid1);
  double beta1    = vf_->getBetaWater(fluid1);
  double pref1    = vf_->getPrefWater(fluid1);
  double alpha2   = vf_->getAlphaWater(fluid2);
  double beta2    = vf_->getBetaWater(fluid2);
  double pref2    = vf_->getPrefWater(fluid2);

  //double T_w, P_g, P_w, U_w, U_g, R_w, R_g;
  double P_1, P_2, R_1, R_2, U_1, U_2, T_1, T_2;
  double P_i, U_i, R_i1, R_i2;
  
  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  int err; 
 
  double pmin1 = vf_->getVarFcnBase(fluid1)->pmin;
  double pmin2 = vf_->getVarFcnBase(fluid2)->pmin;
  double rhomin1 = vf_->getVarFcnBase(fluid1)->rhomin;
  double rhomin2 = vf_->getVarFcnBase(fluid2)->rhomin;


  if (IDi==fluid1) {
    // cell j is tait2
    // cell i is tait1
    R_1  = Vi[0];     R_2  = Vj[0];
    U_1  = vni;       U_2  = vnj;
    P_1  = vf_->getPressure(Vi, IDi);
    P_2  = vf_->getPressure(Vj, IDj);
    T_1  = vf_->computeTemperature(Vi, IDi);
    T_2  = vf_->computeTemperature(Vj, IDj);
    
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,
                        alpha2,beta2,pref2,alpha1,beta1,pref1,err,
                        pmin2,pmin1,rhomin2,rhomin1,tolpre_,nriter_);
    
    Wi[0]  = R_i1;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = T_1;                     Wi[dim+4]  = Wi[4];
    //Wi[4]  = T_2;                     Wi[dim+4]  = Wi[4];
    
    Wj[0]  = R_i2;                      Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];        Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];        Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];        Wj[dim+3]  = Wj[3];
    Wj[4]  = T_2;                       Wj[dim+4]  = Wj[4];
    //Wj[4]  = T_1;                       Wj[dim+4]  = Wj[4];
    
  }else{
    // cell j is tait1
    // cell i is tait2
    R_1  = Vj[0];     R_2  = Vi[0];
    U_1  = vnj;       U_2  = vni;
    P_1  = vf_->getPressure(Vj, IDj);
    P_2  = vf_->getPressure(Vi, IDi);
    T_1  = vf_->computeTemperature(Vj, IDj);
    T_2  = vf_->computeTemperature(Vi, IDi);
    
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,
                        alpha2,beta2,pref2,alpha1,beta1,pref1,err,
                        pmin2,pmin1,rhomin2,rhomin1,tolpre_,nriter_);
    
    Wi[0]  = R_i2;                      Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];        Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];        Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];        Wi[dim+3]  = Wi[3];
    Wi[4]  = T_2;                       Wi[dim+4]  = Wi[4];
    //Wi[4]  = T_1;                       Wi[dim+4]  = Wi[4];
    
    Wj[0]  = R_i1;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    Wj[4]  = T_1;                     Wj[dim+4]  = Wj[4];
    //Wj[4]  = T_2;                     Wj[dim+4]  = Wj[4];
  }

// to update the nodes when they change fluids
// METHOD1: naive approach of averaging the Riemann solution
  /*if(it==1){
    weighti += 1.0;
    weightj += 1.0;
    for (int k=0; k<5; k++){
      rupdatei[k] += Wj[k];
      rupdatej[k] += Wi[k];
    }
  }*/

// METHOD 2 : combine averaging and direction of flow
  if (it == 1 && !isHigherOrder)
    updatePhaseChangingNodeValues(dx, Wi, Wj, weighti, rupdatei, weightj, rupdatej);

  return err;
}

inline 
void LocalRiemannGfmparTaitTait::computeRiemannJacobian(double *Vi, double *Vj,
							int IDi, int IDj, double *nphi,
							double *Wi, double *Wj,
							double dx[3],int it,
							double* dWidWi,double*  dWidWj,
							double* dWjdWi,double*  dWjdWj) {
  
  int dim = 5;
  int k,l;
  
  double alpha1   = vf_->getAlphaWater(fluid1);
  double beta1    = vf_->getBetaWater(fluid1);
  double pref1    = vf_->getPrefWater(fluid1);
  double alpha2   = vf_->getAlphaWater(fluid2);
  double beta2    = vf_->getBetaWater(fluid2);
  double pref2    = vf_->getPrefWater(fluid2);

  double P_1, P_2, U_1, U_2, R_1, R_2;
  double P_i, U_i, R_i1, R_i2;

  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  // 3x3 Jacobians, directly from implicit riemann jacobian
  double dWidWi3[9],  dWidWj3[9], dWjdWj3[9], dWjdWi3[9];

  int err; 
 
  double pmin1 = vf_->getVarFcnBase(fluid1)->pmin;
  double pmin2 = vf_->getVarFcnBase(fluid2)->pmin;
  double rhomin1 = vf_->getVarFcnBase(fluid1)->rhomin;
  double rhomin2 = vf_->getVarFcnBase(fluid2)->rhomin;

  if (IDi==fluid1) {

    //std::cout << "ij" << std::endl << std::endl ;
    // cell i is fluid1
    // cell j is fluid2
    R_2  = Vj[0];     R_1 = Vi[0];
    U_2  = vnj;       U_1 = vni;
    P_2  = vf_->getPressure(Vj, IDj);
    P_1  = vf_->getPressure(Vi, IDi);
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,
                        alpha2,beta2,pref2,alpha1,beta1,pref1,err,
                        pmin2,pmin1,rhomin2,rhomin1,tolpre_,nriter_);
;

    ImplicitRiemann::computeTaitTaitJacobian(P_i, alpha2,beta2,pref2,P_2,R_2, alpha1, beta1,pref1, P_1,R_1, dWjdWj3, dWjdWi3,  dWidWi3, dWidWj3 );
    
  }else{
    // cell i is fluid2
    // cell j is fluid1
    //std::cout << "ji" << std::endl << std::endl ;
    R_2  = Vi[0];     R_1  = Vj[0];
    U_2  = vni;       U_1  = vnj;
    P_2  = vf_->getPressure(Vi, IDi);
    P_1  = vf_->getPressure(Vj, IDj);
    F77NAME(eriemannww)(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,
                        alpha2,beta2,pref2,alpha1,beta1,pref1,err,
                        pmin2,pmin1,rhomin2,rhomin1,tolpre_,nriter_);

    ImplicitRiemann::computeTaitTaitJacobian(P_i, alpha2,beta2,pref2, P_2,R_2,alpha1, beta1,pref1,P_1,R_1, dWidWi3, dWidWj3,  dWjdWj3, dWjdWi3 );
  }

  this->oneDtoThreeD(dWidWi3, dWidWj3,
		     dWjdWi3,dWjdWj3,nphi,
		     dWidWi, dWidWj,
		     dWjdWi, dWjdWj);
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmparJWLJWL : public LocalRiemannGfmpar {

public:
  LocalRiemannGfmparJWLJWL(VarFcn *vf, int tag1, int tag2, MultiFluidData::TypePhaseChange typePhaseChange) : LocalRiemannGfmpar(vf,tag1,tag2,typePhaseChange) {}
  ~LocalRiemannGfmparJWLJWL() { vf_ = 0; }

  int computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *initWi, double *initWj,
                              double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it,bool isHigherOrder);

  void computeRiemannJacobian(double *Vi, double *Vj,
			      int IDi, int IDj, double *nphi,
			      double *Wi, double *Wj,
			      double dx[3],int it,
			      double* dWidWi,double*  dWidWj,
			      double* dWjdWi,double*  dWjdWj);
/*
  void eriemann(double rhol, double ul, double pl, 
                double rhor, double ur, double pr, 
                double &pi, double &ui, double &rhoil, double &rhoir){ 
    eriemannjj(rhol,ul,pl,rhor,ur,pr,pi,ui,rhoil,rhoir); }
*/

 // FS Riemann problem (implemented here just to stop compiler's complaining.)
  int computeRiemannSolution(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparJWLJWL!\n"); return 0; }

  void computeRiemannJacobian(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, double* dWstardU,int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparJWLJWL!\n");}

  void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparJWLJWL!\n");}

private:
  LocalRiemannGfmparJWLJWL();
  void eriemannjj(double rhol, double ul, double pl, 
                  double rhor, double ur, double pr, 
                  double &pi, double &ui, double &rhoil, double &rhoir,int& err,
                  double pcl,double pcr, double rhocl,double rhocr);
};

//----------------------------------------------------------------------------

inline
int LocalRiemannGfmparJWLJWL::computeRiemannSolution(double *Vi, double *Vj,
	 	int IDi, int IDj, double *nphi,
                double *initWi, double *initWj,
		double *Wi, double *Wj,
                double *rupdatei, double *rupdatej, double &weighti, double &weightj,
                double dx[3], int it,bool isHigherOrder)
{

  bool computeRiemannSolutionJWLJWLimplemented = false;
  int dim = 5;
	
  double P_1, P_2, U_1, U_2, R_1, R_2;
  double P_i, U_i, R_i1, R_i2;


  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  int err; 
 
  double pmin1 = vf_->getVarFcnBase(fluid1)->pmin;
  double pmin2 = vf_->getVarFcnBase(fluid2)->pmin;
  double rhomin1 = vf_->getVarFcnBase(fluid1)->rhomin;
  double rhomin2 = vf_->getVarFcnBase(fluid2)->rhomin;

  if (IDi==fluid1) {

    // cell i is fluid1
    // cell j is fluid2
    R_2  = Vj[0];     R_1 = Vi[0];
    U_2  = vnj;       U_1 = vni;
    P_2  = vf_->getPressure(Vj, IDj);
    P_1  = vf_->getPressure(Vi, IDi);

   eriemannjj(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,err,
              pmin2,pmin1,rhomin2,rhomin1);

    Wi[0]  = R_i1;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = P_i;                     Wi[dim+4]  = Wi[4];

    Wj[0]  = R_i2;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    Wj[4]  = P_i;                     Wj[dim+4]  = Wj[4];

  }else{
    // cell i is fluid2
    // cell j is fluid1
    R_2  = Vi[0];     R_1  = Vj[0];
    U_2  = vni;       U_1  = vnj;
    P_2  = vf_->getPressure(Vi, IDi);
    P_1  = vf_->getPressure(Vj, IDj);

    eriemannjj(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,err,
              pmin2,pmin1,rhomin2,rhomin1);

    Wi[0]  = R_i2;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = P_i;                     Wi[dim+4]  = Wi[4];

    Wj[0]  = R_i1;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    Wj[4]  = P_i;                     Wj[dim+4]  = Wj[4];
  }

// to update the nodes when they change fluids
// METHOD1: naive approach of averaging the Riemann solution
  /*if(it==1){
    weighti += 1.0;
    weightj += 1.0;
    for (int k=0; k<5; k++){
      rupdatei[k] += Wj[k];
      rupdatej[k] += Wi[k];
    }
  }*/

// METHOD 2 : combine averaging and direction of flow
  if (it == 1 && !isHigherOrder)
    updatePhaseChangingNodeValues(dx, Wi, Wj, weighti, rupdatei, weightj, rupdatej);

  return err;
}
//----------------------------------------------------------------------------

inline 
void LocalRiemannGfmparJWLJWL::computeRiemannJacobian(double *Vi, double *Vj,
							int IDi, int IDj, double *nphi,
							double *Wi, double *Wj,
							double dx[3],int it,
							double* dWidWi,double*  dWidWj,
							double* dWjdWi,double*  dWjdWj) {
  
  int dim = 5;
  int k,l;
  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  // 3x3 Jacobians, directly from implicit riemann jacobian
  double dWidWi3[9],  dWidWj3[9], dWjdWj3[9], dWjdWi3[9];

  if (IDi==fluid1) {

    //std::cout << "ij" << std::endl << std::endl ;
    // cell i is fluid1
    // cell j is fluid2

    ImplicitRiemann::computeJwlJwlJacobian(vf_, IDi, IDj, Vi, Vj, Wi,Wj, dWidWi3, dWidWj3,  dWjdWj3, dWjdWi3 );
    
  }else{
    // cell i is fluid2
    // cell j is fluid1

    ImplicitRiemann::computeJwlJwlJacobian(vf_, IDj, IDi, Vj, Vi, Wj,Wi, dWjdWj3, dWjdWi3,  dWidWi3, dWidWj3 );

  }

  this->oneDtoThreeD(dWidWi3, dWidWj3,
		     dWjdWi3,dWjdWj3,nphi,
		     dWidWi, dWidWj,
		     dWjdWi, dWjdWj);
}

//----------------------------------------------------------------------------
inline
void LocalRiemannGfmparJWLJWL::eriemannjj(double rhol, double ul, double pl, 
                                          double rhor, double ur, double pr, 
                                          double &pi, double &ui,  
                                          double &rhoil, double &rhoir,int& err,
                                          double pcl,double pcr, double rhocl,double rhocr){

//initialize
  double uil, uir, pil, pir, duil, duir, dpil, dpir;
  double jacobian[4];/* uil, uir, pil, pir*/
  double function[2];
  double increment[2];
  bool convergence = false;
  double eps = 1.e-6;
  int MaxIts = 100;
  int it = 0;

  double pcut = std::max<double>(  pcl,   pcr);
  double rcut = std::max<double>(rhocl, rhocr);

  double vl  = 1.0/rhol;
  double vr  = 1.0/rhor;
  double vil = vl;
  double vir = vr;
  double omegal = vf_->getOmega(fluid2);
  double omegar = vf_->getOmega(fluid1);
  double omp1ooml = (omegal+1.0)/omegal;
  double omp1oomr = (omegar+1.0)/omegar;
  double frhol = vf_->computeFrho(1.0/vl,fluid2);
  double frhor = vf_->computeFrho(1.0/vr,fluid1);
  double frhoil = frhol;
  double frhoir = frhor;
  double frhopil = vf_->computeFrhop(1.0/vl,fluid2);
  double frhopir = vf_->computeFrhop(1.0/vr,fluid1);
//check vacuum ?

  err = 0;

//start newton iteration loop
  while(!convergence){
    //fprintf(stdout, "%e %e %e %e\n", 1.0/vil, 1.0/vir, ui, pi);

  //compute left  term (shock or rarefaction)
    if( vil < vl){
      //fprintf(stdout, "leftshock\n");
      frhoil  = vf_->computeFrho(1.0/vil,fluid2);
      frhopil = vf_->computeFrhop(1.0/vil,fluid2);
      shockJWL(-1.0, omegal, omp1ooml, frhol, frhoil, frhopil, vl, ul, pl, vil, uil, pil, duil, dpil);
    }else{
      //fprintf(stdout, "leftraref\n");
      rarefactionJWL(-1.0, vl, ul, pl, vil, uil, pil, duil, dpil, 
                     MultiFluidData::RK2, 1);
    }
  //compute right term (shock or rarefaction)
    if( vir < vr){
      //fprintf(stdout, "rightshock\n");
      frhoir  = vf_->computeFrho(1.0/vir,fluid1);
      frhopir = vf_->computeFrhop(1.0/vir,fluid1);
      shockJWL(1.0, omegar, omp1oomr, frhor, frhoir, frhopir, vr, ur, pr, vir, uir, pir, duir, dpir);
    }
    else{
      //fprintf(stdout, "rightraref\n");
      rarefactionJWL(1.0, vr, ur, pr, vir, uir, pir, duir, dpir, 
                     MultiFluidData::RK2, 1);
    }

  //solve2x2System: function = jacobian*increment
    function[0] = uil-uir;
    function[1] = pil-pir;
    jacobian[0] = duil; jacobian[1] = -duir;
    jacobian[2] = dpil; jacobian[3] = -dpir;
    increment[0] = 0.0; increment[1] = 0.0;
    
    solve2x2System(jacobian,function,increment);

  //update values and check bounds
    //fprintf(stdout, "1 -- vil = %e and vir = %e\n", vil, vir);
    //fprintf(stdout, "11-- vil = %e and vir = %e\n", vil-increment[0], vir-increment[1]);
    if(vil - increment[0] < 0.0)
      increment[0] = 0.5*vil;
    if(vir - increment[1] < 0.0)
      increment[1] = 0.5*vir;
    
    vil -= increment[0];
    vir -= increment[1];
    //fprintf(stdout, "2 -- vil = %e and vir = %e\n", vil, vir);
    //
    if(vil < vl){ // at next iteration, leftrarefaction => ensures that some conditions are fulfilled
      double temp = omegal*vl/(omegal+2.0);
      if(vil<temp)
        vil = 0.5*(vil+increment[0]+temp);
    }
    if(vir < vr){ // at next iteration, rightrarefaction => ensures that some conditions are fulfilled
      double temp = omegar*vr/(omegar+2.0);
      if(vir<temp)
        vir = 0.5*(vir+increment[1]+temp);
    }

    if (1.0/vil < rhocl)
      vil = 1.0/rhocl;
    if (1.0/vir < rhocr)
      vir = 1.0/rhocr;

    //fprintf(stdout, "3 -- vil = %e and vir = %e\n", vil, vir);
    it++;

  //check convergence criterion
    if(fabs(increment[0])<eps*fabs(vil) &&
       fabs(increment[1])<eps*fabs(vir) )
      convergence = true;
    if(it>MaxIts) break;


  }//end newton iteration loop
  if(convergence){
    //fprintf(stdout, "riemann has converged to an approximate solution in %d iterations\n", it);
    rhoil = 1.0/vil;
    rhoir = 1.0/vir;
    ui    = 0.5*(uil+uir);
    pi    = 0.5*(pil+pir);
    pi = std::max<double>(pcut, pi);
  }else{
    fprintf(stdout, "***WARNING: ERS for JWL-JWL did not converge in %d iterations to the desired tolerance %e \n", it, eps);
    err = 1;
  }

  /*if(pi <= pcut || rhoil <= rcut || rhoir <= rcut){
    std::cout << "*** ERROR ERS JJ: detected too small density or pressure " << std::endl;
    std::cout << " rL, rR, Pi " << std::endl;
    std::cout << rhoil << " " << rhoir << " " << " " << pi << std::endl; 
    err = 1;
  }*/ 

}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

class LocalRiemannGfmparGasJWL : public LocalRiemannGfmpar {

private:
  MultiFluidData::RiemannComputation riemannComputationType_;
  SparseGridCluster *sgCluster_;

public:
  LocalRiemannGfmparGasJWL(VarFcn *vf, int tag1, int tag2, SparseGridCluster *sgCluster, 
                           MultiFluidData::RiemannComputation riemannComputation,double rfac,double refdensity,double refentropy,
                           MultiFluidData::TypePhaseChange typePhaseChange = MultiFluidData::RIEMANN_SOLUTION) : 
    LocalRiemannGfmpar(vf,tag1,tag2,typePhaseChange,rfac,refdensity,refentropy) {
    riemannComputationType_ = riemannComputation;
    sgCluster_ = sgCluster;
  }
  ~LocalRiemannGfmparGasJWL(){ vf_ = 0; sgCluster_ = 0; }

  int computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *initWi, double *initWj,
                              double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it,bool isHigherOrder);

  void computeRiemannJacobian(double *Vi, double *Vj,
			      int IDi, int IDj, double *nphi,
			      double *Wi, double *Wj,
			      double dx[3],int it,
			      double* dWidWi,double*  dWidWj,
			      double* dWjdWi,double*  dWjdWj);
/*
  void eriemann(double rhol, double ul, double pl, 
                double rhor, double ur, double pr, 
                double &pi, double &ui, double &rhoil, double &rhoir){
    eriemanngj(rhol,ul,pl,rhor,ur,pr,pi,ui,rhoil,rhoir,-1.0,-1.0); }
*/
  void eriemanngj_wrapper(double *in, double *res, double *para);
  void riemannInvariantGeneral1stOrder_wrapper(
                double *in, double *res, double *para);
  void riemannInvariantGeneral2ndOrder_wrapper(
                double *in, double *res, double *para);

 // FS Riemann problem (implemented here just to stop compiler's complaining.)
  int computeRiemannSolution(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparGasJWL!\n"); return 0; }

  void computeRiemannJacobian(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, double* dWstardU,int Id = 0) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparGasJWL!\n");}

  void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it) {
    fprintf(stderr,"ERROR: Should not call the FS Riemann solver in LocalRiemannGfmparGasJWL!\n");}


private:
  LocalRiemannGfmparGasJWL();

protected:
  void eriemanngj_selector(double rhol, double ul, double pl, 
                           double rhor, double ur, double pr, 
                           double &pi, double &ui, double &rhoil, double &rhoir,
                           double initrhol, double initrhor,int& err,
                           double pcl,double pcr, double rhocl,double rhocr);
  bool eriemanngj(double rhol, double ul, double pl, 
                  double rhor, double ur, double pr, 
                  double &pi, double &ui, double &rhoil, double &rhoir,
                  double initrhol, double initrhor,int& err,
                  double pcl,double pcr, double rhocl,double rhocr);
 
  int riemannInvariantGeneralTabulation(double *in, double *res);
  bool vacuum(const double rhol, const double ul, const double pl,
              const double rhor, const double ur, const double pr,
              double vacuumValues[6]);
  double sgZeroDensityPJwlDensity(const double density, 
                                  const double pressure,
                                  const double rho_c0);
  double pressureEqGasDensity(const double gasDensity, 
                              const double gasPressure,
                              const double jwlDensity, 
                              const double jwlPressure,
                              const double interfacialJwlDensity);
};

//----------------------------------------------------------------------------

inline
int LocalRiemannGfmparGasJWL::computeRiemannSolution(double *Vi, double *Vj,
	 	int IDi, int IDj, double *nphi,
                double *initWi, double *initWj,
		double *Wi, double *Wj,
                double *rupdatei, double *rupdatej, double &weighti, double &weightj,
                double dx[3], int it,bool isHigherOrder)
{

  int dim = 5;
	
  double P_1, P_2, U_1, U_2, R_1, R_2;
  double P_i, U_i, R_i1, R_i2;


  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  int err; 
 
  double pmin1 = vf_->getVarFcnBase(fluid1)->pmin;
  double pmin2 = vf_->getVarFcnBase(fluid2)->pmin;
  double rhomin1 = vf_->getVarFcnBase(fluid1)->rhomin;
  double rhomin2 = vf_->getVarFcnBase(fluid2)->rhomin;
 
  if (IDi==fluid1) {

    // cell i is fluid1
    // cell j is fluid2
    R_2  = Vj[0];     R_1 = Vi[0];
    U_2  = vnj;       U_1 = vni;
    P_2  = vf_->getPressure(Vj, IDj);
    P_1  = vf_->getPressure(Vi, IDi);

    //eriemanngj_selector(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,initWj[0],initWi[0]); 
    eriemanngj_selector(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,R_2,R_1,err,
              pmin2,pmin1,rhomin2,rhomin1);
 
    initWi[0] = R_i1;
    initWi[1] = U_i;
    initWi[2] = P_i;
    initWj[0] = R_i2;
    initWj[1] = U_i;
    initWj[2] = P_i;

    Wi[0]  = R_i1;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = P_i;                     Wi[dim+4]  = Wi[4];

    Wj[0]  = R_i2;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    Wj[4]  = P_i;                     Wj[dim+4]  = Wj[4];

  }else{
    // cell i is fluid2
    // cell j is fluid1
    R_2  = Vi[0];     R_1  = Vj[0];
    U_2  = vni;       U_1  = vnj;
    P_2  = vf_->getPressure(Vi, IDi);
    P_1  = vf_->getPressure(Vj, IDj);

    //eriemanngj_selector(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,initWi[0],initWj[0]); 
    eriemanngj_selector(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,R_2,R_1,err,
              pmin2,pmin1,rhomin2,rhomin1);
 
    initWi[0] = R_i2;
    initWi[1] = U_i;
    initWi[2] = P_i;
    initWj[0] = R_i1;
    initWj[1] = U_i;
    initWj[2] = P_i;

    Wi[0]  = R_i2;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = P_i;                     Wi[dim+4]  = Wi[4];

    Wj[0]  = R_i1;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    Wj[4]  = P_i;                     Wj[dim+4]  = Wj[4];
  }

// to update the nodes when they change fluids
// METHOD1: naive approach of averaging the Riemann solution
  /*if(it==1){
    weighti += 1.0;
    weightj += 1.0;
    for (int k=0; k<5; k++){
      rupdatei[k] += Wj[k];
      rupdatej[k] += Wi[k];
    }
  }*/

// METHOD 2 : combine averaging and direction of flow
  if (it == 1 && !isHigherOrder)
    updatePhaseChangingNodeValues(dx, Wi, Wj, weighti, rupdatei, weightj, rupdatej);

  return err;
}
inline 
void LocalRiemannGfmparGasJWL::computeRiemannJacobian(double *Vi, double *Vj,
							int IDi, int IDj, double *nphi,
							double *Wi, double *Wj,
							double dx[3],int it,
							double* dWidWi,double*  dWidWj,
							double* dWjdWi,double*  dWjdWj) {
  
	
  double P_1, P_2, U_1, U_2, R_1, R_2;
  double P_i, U_i, R_i1, R_i2;

  int dim = 5;
  int k,l;
  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  // 3x3 Jacobians, directly from implicit riemann jacobian
  double dWidWi3[9],  dWidWj3[9], dWjdWj3[9], dWjdWi3[9];

  if (IDi==fluid1) {

    //std::cout << "ij" << std::endl << std::endl ;
    // cell i is fluid1
    // cell j is fluid2

    int err = 1;
    if (riemannComputationType_==MultiFluidData::TABULATION2) {
      double dVdv[2];
      err = rarefactionJWLderivs(-1.0, 1.0/Vj[0], vnj, Vj[4], 1.0/Wj[0] , dVdv, sgCluster_ );
      if(!err) ImplicitRiemann::computeGasJwlJacobian(vf_, IDi, IDj, Vi, Vj, Wi,Wj, dWidWi3, dWidWj3,  dWjdWj3, dWjdWi3, &dVdv[0]  );
    }
    if(err)
      ImplicitRiemann::computeGasJwlJacobian(vf_, IDi, IDj, Vi, Vj, Wi,Wj, dWidWi3, dWidWj3,  dWjdWj3, dWjdWi3, NULL);
  
    dWidWi3[1] *= -1.0;
    dWidWi3[3] *= -1.0;
    dWidWi3[5] *= -1.0;
    dWidWi3[7] *= -1.0;
    dWidWj3[1] *= -1.0;
    dWidWj3[3] *= -1.0;
    dWidWj3[5] *= -1.0;
    dWidWj3[7] *= -1.0;
    dWjdWi3[1] *= -1.0;
    dWjdWi3[3] *= -1.0;
    dWjdWi3[5] *= -1.0;
    dWjdWi3[7] *= -1.0;
    dWjdWj3[1] *= -1.0;
    dWjdWj3[3] *= -1.0;
    dWjdWj3[5] *= -1.0;
    dWjdWj3[7] *= -1.0;

  }else{
    // cell i is fluid2
    // cell j is fluid1

    //std::cout << "ji" << std::endl << std::endl ;

    int err = 1;
    if (riemannComputationType_==MultiFluidData::TABULATION2) {
      double dVdv[2];
      err = rarefactionJWLderivs(-1.0, 1.0/Vi[0], vni, Vi[4], 1.0/Wi[0] , dVdv, sgCluster_ );
      if(!err) ImplicitRiemann::computeGasJwlJacobian(vf_, IDj, IDi, Vj, Vi, Wj,Wi, dWjdWj3, dWjdWi3,  dWidWi3, dWidWj3, &dVdv[0]  );
    }
    if(err)
      ImplicitRiemann::computeGasJwlJacobian(vf_, IDj, IDi, Vj, Vi, Wj,Wi, dWjdWj3, dWjdWi3,  dWidWi3, dWidWj3, NULL );
    
    dWidWi3[1] *= -1.0;
    dWidWi3[3] *= -1.0;
    dWidWi3[5] *= -1.0;
    dWidWi3[7] *= -1.0;
    dWidWj3[1] *= -1.0;
    dWidWj3[3] *= -1.0;
    dWidWj3[5] *= -1.0;
    dWidWj3[7] *= -1.0;
    dWjdWi3[1] *= -1.0;
    dWjdWi3[3] *= -1.0;
    dWjdWi3[5] *= -1.0;
    dWjdWi3[7] *= -1.0;
    dWjdWj3[1] *= -1.0;
    dWjdWj3[3] *= -1.0;
    dWjdWj3[5] *= -1.0;
    dWjdWj3[7] *= -1.0;

  }

  this->oneDtoThreeD(dWidWi3, dWidWj3,
		     dWjdWi3,dWjdWj3,nphi,
		     dWidWi, dWidWj,
		     dWjdWi, dWjdWj);
}

//----------------------------------------------------------------------------
inline
void LocalRiemannGfmparGasJWL::eriemanngj_wrapper(
                               double *in, double *res, double *para)
{

  double dummy1, dummy2;
  int err;
  eriemanngj(in[0], 0.0, in[1], in[2], in[4], in[3], dummy1, dummy2, res[0], res[1], -1.0, -1.0,
             err, -1.0,-1.0,-1.0,-1.0);

}

//----------------------------------------------------------------------------
inline
void LocalRiemannGfmparGasJWL::eriemanngj_selector(
                               double rhol, double ul, double pl, 
                               double rhor, double ur, double pr, 
                               double &pi, double &ui,  
                               double &rhoil, double &rhoir,
                               double initrhol, double initrhor,int& err,
                               double pcl,double pcr, double rhocl,double rhocr)
{

  err = 0;

  if(riemannComputationType_==MultiFluidData::TABULATION5){
    double *in  = new double[5]; in[0]=rhol;in[1]=pl;in[2]=rhor;in[3]=pr;in[4]=ur-ul;
    double *res = new double[2]; res[0]=0.0;res[1]=0.0;
    sgCluster_->interpolate(1,&in,&res);
    rhoil = fmax(res[0],0.0); rhoir = fmax(res[1],0.0);
    double d[2]; //dummy variable
    double uir, pir, uil, pil;

    double omegal = vf_->getOmega(fluid2);
    double omp1ooml = (omegal+1.0)/omegal;
    double frhol = vf_->computeFrho(rhol,fluid2);
    double gamr = vf_->getGamma(fluid1);
    double prefr = vf_->getPressureConstant(fluid1);
    double gam1r = vf_->getGamma(fluid1)-1.0;
    double gamogam1r = gamr/gam1r;
    double Vr[5] = { rhor, ur, 0.0, 0.0, pr };
    double cr = vf_->computeSoundSpeed(Vr,fluid1);

    if( rhoil > rhol){
      double frhoil  = vf_->computeFrho(rhoil,fluid2);
      double frhopil = vf_->computeFrhop(rhoil,fluid2);
      shockJWL(-1.0, omegal, omp1ooml, frhol, frhoil, frhopil, 1.0/rhol, ul, pl, 1.0/rhoil, uil, pil, d[0], d[1]);
    }else
      rarefactionJWL(-1.0, 1.0/rhol, ul, pl, 1.0/rhoil, uil, pil, d[0], d[1], riemannComputationType_,1);
    if( rhoir > rhor)
      shockGAS(1.0, gamogam1r, prefr, 1.0/rhor, ur, pr, 1.0/rhoir, uir, pir, d[0], d[1]);
    else
      rarefactionGAS(1.0, gamr, gam1r, prefr, cr, 1.0/rhor, ur, pr, 1.0/rhoir, uir, pir, d[0], d[1]);

    ui = 0.5*(uil+uir);
    pi = 0.5*(pil+pir);

    pi = std::max<double>(pi, std::max<double>(pcl,pcr));
    rhoil = std::max<double>(rhoil, rhocl);
    rhoir = std::max<double>(rhoir, rhocr);

  /*double pcut = std::max<double>(  pcl,   pcr);
    double rcut = std::max<double>(rhocl, rhocr);

    if(pi <= pcut || rhoil <= rcut || rhoir <= rcut){
      std::cout << "*** ERROR ERS GJ: detected too small density or pressure " << std::endl;
      std::cout << " rL, rR, Pi " << std::endl;
      std::cout << rhoil << " " << rhoir << " " << " " << pi << std::endl; 
      err = 1;
    }*/ 

  }
  else
    eriemanngj(rhol,ul,pl,rhor,ur,pr,pi,ui,rhoil,rhoir,initrhol,initrhor, err,pcl,pcr,rhocl,rhocr);

}

//----------------------------------------------------------------------------
inline
bool LocalRiemannGfmparGasJWL::eriemanngj(double rhol, double ul, double pl, 
                                          double rhor, double ur, double pr, 
                                          double &pi, double &ui,  
                                          double &rhoil, double &rhoir,
                                          double initrhol, double initrhor,int& err,
                                          double pcl,double pcr, double rhocl,double rhocr){
// left  -- JWL -- phi = -1.0
// right -- GAS -- phi = +1.0
  int verbose = -1;
  if(verbose>0){
    fprintf(stdout, "---- new Riemann ----\n");
    fprintf(stdout, "initial rhoil, rhoir = %e %e\n", rhol, rhor);
    fprintf(stdout, "initial vil,   vir   = %e %e\n", 1.0/rhol, 1.0/rhor);
  }

//initialize
  double uil, uir, pil, pir, duil, duir, dpil, dpir;
  double jacobian[4];/* uil, uir, pil, pir*/
  double function[2];
  double increment[2];
  bool convergence = false;
  double eps = 1.e-3;
  int MaxIts = 100;
  int it = 0;
  double relaxationFactorJwl = relaxFactorJwl;//1.0; //0.85; // must be between 0 and 1
  double relaxationFactorGas = relaxFactorJwl;//1.0; //0.85; // must be between 0 and 1
  int count = 0;

  double pcut = std::max<double>(  pcl,   pcr);
  double rcut = std::max<double>(rhocl, rhocr);

  double vl  = 1.0/rhol;
  double vr  = 1.0/rhor;
  double vil = vl;
  double vir = vr;
  vil = initrhol>0.0 ? 1.0/initrhol : vl;
  vir = initrhor>0.0 ? 1.0/initrhor : vr;

  double omegal = vf_->getOmega(fluid2);
  double omp1ooml = (omegal+1.0)/omegal;
  double frhol = vf_->computeFrho(1.0/vl,fluid2);
  double frhoil = frhol;
  double frhopil = vf_->computeFrhop(1.0/vl,fluid2);


  double gamr = vf_->getGamma(fluid1);
  double prefr = vf_->getPressureConstant(fluid1);
  double gam1r = vf_->getGamma(fluid1)-1.0;
  double gamogam1r = gamr/gam1r;
  double Vr[5] = { 1.0/vr, ur, 0.0, 0.0, pr };
  double cr = vf_->computeSoundSpeed(Vr,fluid1);
  double pastiterates[100][2];

  err = 0;

//check vacuum
  if(verbose>4) fprintf(stdout, "checking vacuum possibilities\n");
  double vacuumValues[6]; /* rhoil, uil, pil, rhoir, uir, pir */
  vacuumValues[0] = -1.0; // positive if proper vacuum values are computed
  bool checkVacuumValues = false;
  if(checkVacuumValues){
    if(vacuum(rhol,ul,pl,rhor,ur,pr,vacuumValues)){
      if(verbose>-1){
        fprintf(stdout, "rhoil_vac = %e and rhoir_vac = %e\n", vacuumValues[0], vacuumValues[3]);
        fprintf(stdout, "uil_vac   = %e and uir_vac   = %e\n", vacuumValues[1], vacuumValues[4]);
        fprintf(stdout, "pil_vac   = %e and pir_vac   = %e\n", vacuumValues[2], vacuumValues[5]);
      }
      rhoil = vacuumValues[0];
      rhoir = vacuumValues[3];
      ui    = 0.5*(vacuumValues[1]+vacuumValues[4]);
      pi    = 0.5*(vacuumValues[2]+vacuumValues[5]);
      return true;
    }
    if(verbose>4) fprintf(stdout, "checking vacuum possibilities -- DONE\n");
  }else{
    if(verbose>4) fprintf(stdout, "no checking of vacuum possibilities\n");
  }

  double res = 1.0e20;

//start newton iteration loop
  while(!convergence){
    if(verbose>0) fprintf(stdout, "\n");
    int status = 1;

  //compute left  JWL-term (shock or rarefaction)
    if( vil < vl){
      if(verbose>0) fprintf(stdout, "shockJWL\n");
      frhoil  = vf_->computeFrho(1.0/vil,fluid2);
      frhopil = vf_->computeFrhop(1.0/vil,fluid2);
      shockJWL(-1.0, omegal, omp1ooml, frhol, frhoil, frhopil, vl, ul, pl, vil, uil, pil, duil, dpil);
    }else{
      if(verbose>0) fprintf(stdout, "rarefactionJWL\n");
      status = rarefactionJWL(-1.0, vl, ul, pl, vil, uil, pil, duil, dpil, riemannComputationType_,1);
    }
  //compute right GAS-term (shock or rarefaction)
    if( vir < vr){
      if(verbose>0) fprintf(stdout, "shockGAS\n");
      shockGAS(1.0, gamogam1r, prefr, vr, ur, pr, vir, uir, pir, duir, dpir);
    }
    else{
      if(verbose>0) fprintf(stdout, "rarefactionGAS\n");
      rarefactionGAS(1.0, gamr, gam1r, prefr, cr, vr, ur, pr, vir, uir, pir, duir, dpir);
    }
    if(verbose>1){
      fprintf(stdout, "uil  = %e and uir  = %e\n", uil, uir);
      fprintf(stdout, "pil  = %e and pir  = %e\n", pil, pir);
      fprintf(stdout, "duil = %e and duir = %e\n", duil, duir);
      fprintf(stdout, "dpil = %e and dpir = %e\n", dpil, dpir);
    }

    pastiterates[it][0] = vil;
    pastiterates[it][1] = vir;

    if (!status) {
      fprintf(stdout, "$$$$ status = 0\n");
      fprintf(stdout, "rhol, ul, pl = %e %e %e\n", rhol, ul, pl);
      fprintf(stdout, "rhor, ur, pr = %e %e %e\n", rhor, ur, pr);
      fprintf(stdout, "rhoil  = %e and rhoir  = %e\n", 1/vil, 1/vir);
      fprintf(stdout, "initrhoil  = %e and initrhoir  = %e\n", initrhol, initrhor);
      fprintf(stdout, "uil  = %e and uir  = %e\n", uil, uir);
      fprintf(stdout, "pil  = %e and pir  = %e\n", pil, pir);
      fprintf(stdout, "duil = %e and duir = %e\n", duil, duir);
      fprintf(stdout, "dpil = %e and dpir = %e\n", dpil, dpir);
      for (int kk = 0; kk <= it; ++kk) {
	fprintf(stdout, "it = %i: vil = %e, vir = %e\n",kk,pastiterates[kk][0],pastiterates[kk][1]);
      }
      exit(1);
    }

    //solve2x2System: function = jacobian*increment
    function[0] = uil-uir;
    function[1] = pil-pir;
    jacobian[0] = duil; jacobian[1] = -duir;
    jacobian[2] = dpil; jacobian[3] = -dpir;
    increment[0] = 0.0; increment[1] = 0.0;
    
    /*double res2 = fabs(function[1])+fabs(function[0]);
    if (res2 < res)
      res = res2;
    else {
      vil += increment[0]*0.5;
      vir += increment[1]*0.5;
      increment[0]*=0.5;
      increment[1]*=0.5;
      continue;
    }
    */
    bool solved = solve2x2System(jacobian,function,increment);
    if(!solved){
      fprintf(stdout, "$$$$\n");
      fprintf(stdout, "rhol, ul, pl = %e %e %e\n", rhol, ul, pl);
      fprintf(stdout, "rhor, ur, pr = %e %e %e\n", rhor, ur, pr);
      fprintf(stdout, "rhoil  = %e and rhoir  = %e\n", 1/vil, 1/vir);
      fprintf(stdout, "initrhoil  = %e and initrhoir  = %e\n", initrhol, initrhor);
      fprintf(stdout, "uil  = %e and uir  = %e\n", uil, uir);
      fprintf(stdout, "pil  = %e and pir  = %e\n", pil, pir);
      fprintf(stdout, "duil = %e and duir = %e\n", duil, duir);
      fprintf(stdout, "dpil = %e and dpir = %e\n", dpil, dpir);
      for (int kk = 0; kk <= it; ++kk) {
	fprintf(stdout, "it = %i: vil = %e, vir = %e\n",kk,pastiterates[kk][0],pastiterates[kk][1]);
      }
      rarefactionGAS(1.0, gamr, gam1r, prefr, cr, vr, ur, pr, vir, uir, pir, duir, dpir,1);
    }

    if(verbose>2) fprintf(stdout, "dvil = %e and dvir = %e\n", -increment[0],-increment[1]);

    //update values and check bounds

    if(verbose>3) fprintf(stdout, "increment/v = %e %e\n", increment[0]/vil, increment[1]/vir);

    // prevent large increases
    if(-increment[0]>2.0*vil) increment[0] = -2.0*vil;
    if(-increment[1]>2.0*vir) increment[1] = -2.0*vir;
    // prevent large decreases
    if(increment[0]>0.5*vil) increment[0] = 0.5*vil;
    if(increment[1]>0.5*vir) increment[1] = 0.5*vir;

    increment[0] *= relaxationFactorJwl;
    increment[1] *= relaxationFactorGas;
    
    vil -= increment[0];
    vir -= increment[1];
    if(verbose>2) fprintf(stdout, "2 -- vil = %e and vir = %e\n", vil, vir);

    if(vil < vl){ // at next iteration, leftrarefaction => ensures that some conditions are fulfilled
      double temp = omegal*vl/(omegal+2.0);
      if(vil<temp){
        vil += increment[0];
        vir += increment[1];
        double alpha = -0.5*(temp-vil)/increment[0];
        increment[0] *= alpha;
        increment[1] *= alpha;
        vil -= increment[0];
        vir -= increment[1];
        count++;
      }
    }
    if(vir < vr){ // at next iteration, rightrarefaction => ensures that some conditions are fulfilled
      double temp = (gamr-1.0)/(gamr+1.0)*vr;
      if(vir<temp){
        vil += increment[0];
        vir += increment[1];
        double alpha = -0.5*(temp-vir)/increment[1];
        increment[0] *= alpha;
        increment[1] *= alpha;
        vil -= increment[0];
        vir -= increment[1];
        count++;
      }
    }
    if (1.0/vil < rhocl)
      vil = 1.0/rhocl;
    if (1.0/vir < rhocr)
      vir = 1.0/rhocr;


    if(verbose>2) fprintf(stdout, "3 -- vil = %e and vir = %e\n", vil, vir);
    //check - in case of rarefaction at next iteration, 1.0/rhoil may not be above a certain value
    if(vacuumValues[0] > 0.0 && vil>1.0/vacuumValues[0]){ // vacuumValues is negative if it does not contain any proper value(see declaration and definition above)
      vil += increment[0];
      increment[0] = -0.5*(1.0/vacuumValues[0] - vil);
      vil -= increment[0];
    }
    if(verbose>2) fprintf(stdout, "4 -- vil = %e and vir = %e\n", vil, vir);
    if(verbose>0) fprintf(stdout, "rhoil = %e and rhoir = %e\n", 1.0/vil, 1.0/vir);
    it++;

  //check convergence criterion
    if(fabs(increment[0])<eps*fabs(vil) &&
       fabs(increment[1])<eps*fabs(vir) )
      convergence = true;
    if(it>MaxIts) break;


  }//end newton iteration loop

  if( vil < vl){
    frhoil  = vf_->computeFrho(1.0/vil,fluid2);
    frhopil = vf_->computeFrhop(1.0/vil,fluid2);
    shockJWL(-1.0, omegal, omp1ooml, frhol, frhoil, frhopil, vl, ul, pl, vil, uil, pil, duil, dpil);
  }else rarefactionJWL(-1.0, vl, ul, pl, vil, uil, pil, duil, dpil, riemannComputationType_,0);
  
  if( vir < vr) shockGAS(1.0, gamogam1r, prefr, vr, ur, pr, vir, uir, pir, duir, dpir);
  else rarefactionGAS(1.0, gamr, gam1r, prefr, cr, vr, ur, pr, vir, uir, pir, duir, dpir);

  rhoil = 1.0/vil;
  rhoir = 1.0/vir;
  ui    = 0.5*(uil+uir);
  pi    = 0.5*(pil+pir);
  pi = std::max<double>(pcut, pi);

  //std::cout << "JWL results: rhoil = " << 1.0/vil << " rhoir = "  << 1.0/vir << " ui = " << ui << " pi = " << pi << std::endl;

  //double Vpp[5] = {rhoil,0,0,0,pi};
  //double cll = vf_->computeSoundSpeed(Vpp,fluid2);

  if(convergence){
    if(verbose>-1) fprintf(stdout, "riemann has converged to an approximate solution in %d iterations\n", it);
  }else{
    fprintf(stderr, "riemann solver did not converged\n");
    fprintf(stderr, "Warning: solution will be state given by vacuum\n");
    fprintf(stderr, "rhol, ul, pl = %e %e %e\n", rhol, ul, pl);
    fprintf(stderr, "rhor, ur, pr = %e %e %e\n", rhor, ur, pr);
    fprintf(stderr, "rhoil  = %e and rhoir  = %e\n", 1/vil, 1/vir);
    fprintf(stderr, "initrhoil  = %e and initrhoir  = %e\n", initrhol, initrhor);
    fprintf(stderr, "uil  = %e and uir  = %e\n", uil, uir);
    fprintf(stderr, "pil  = %e and pir  = %e\n", pil, pir);
    fprintf(stderr, "duil = %e and duir = %e\n", duil, duir);
    fprintf(stderr, "dpil = %e and dpir = %e\n", dpil, dpir);
    for (int kk = 0; kk <= it; ++kk) {
      fprintf(stderr, "it = %i: vil = %e, vir = %e\n",kk,pastiterates[kk][0],pastiterates[kk][1]);
    }

    fflush(stderr);

    rhoil = vacuumValues[0];
    rhoir = vacuumValues[3];
    uil   = vacuumValues[1];
    uir   = vacuumValues[4];
    pil   = vacuumValues[2];
    pir   = vacuumValues[5];
    ui    = 0.5*(uil+uir);
    pi    = 0.5*(pil+pir);
    err = 1;
    if(verbose>-1) fprintf(stdout, "Warning: uil = %e and uir = %e\n", uil, uir);
  }

  if(verbose>-1){
    fprintf(stdout, "rhol, ul, pl = %e %e %e\n", rhol, ul, pl);
    fprintf(stdout, "rhor, ur, pr = %e %e %e\n", rhor, ur, pr);
    fprintf(stdout, "initrhol, initrhor = %e %e\n", initrhol, initrhor);
    fprintf(stdout, "rhoil  = %e and rhoir  = %e\n", 1/vil, 1/vir);
    fprintf(stdout, "uil  = %e and uir  = %e\n", uil, uir);
    fprintf(stdout, "pil  = %e and pir  = %e\n", pil, pir);
    fprintf(stdout, "duil = %e and duir = %e\n", duil, duir);
    fprintf(stdout, "dpil = %e and dpir = %e\n", dpil, dpir);
  }

  /*if(pi <= pcut || rhoil <= rcut || rhoir <= rcut){
    std::cout << "*** ERROR ERS GJ: detected too small density or pressure " << std::endl;
    std::cout << " rL, rR, Pi " << std::endl;
    std::cout << rhoil << " " << rhoir << " " << " " << pi << std::endl; 
    err = 1;
  }*/

  if(convergence) {    
    return true;
  } else {
    return false;
  }


}

//----------------------------------------------------------------------------
inline
bool LocalRiemannGfmparGasJWL::vacuum(const double rhol, const double ul, const double pl,
                                      const double rhor, const double ur, const double pr,
                                      double vacuumValues[6]){
// notation: JWL on the left and SG on the right
// remember vacuum can occur only between two rarefaction waves, thus decrease of densities

// 1st step: find JWL-density for which there is loss of positivity of c^2 in JWL gas
  double min1 = jwlZeroSoundSpeedJwlDensity(rhol,pl); // returns -1 if none found

// 2nd step: find JWL-density for which the SG-density would become zero when
//           expressing equality of pressures on both sides of interface
//           equivalent to JWL-density for JWL-pressure is below the
//           lowest value of the SG-pressure
  double min2 = sgZeroDensityPJwlDensity(rhol,pl,min1);    // returns -1 if none found

// 3rd step: find max3, JWL-density for which the SG-density would become zero when
//           expressing equality of velocities on both sides of interface
// WARNING: not done because it would be too costly. Way things are computed
//          defines a zero sg-density for JWL-density above max3

  //fprintf(stdout, "min1 = %e - min2 = %e\n", min1, min2);

// compute SG-density corresponding to JWL-density-bound = max(0,min1,min2)
  double rhoil_vac, rhoir_vac;
  if(min1 < 0 && min2 < 0){
    rhoil_vac = 1.0e-14; // 0.0
    rhoir_vac = pressureEqGasDensity(rhor,pr,rhol,pl,0.0);
  }else if(min1 > 0 && min2 < 0){
    rhoil_vac = min1;
    rhoir_vac = pressureEqGasDensity(rhor,pr,rhol,pl,min1);
  }else if(min1 < 0 && min2 > 0){
    rhoil_vac = min2;
    rhoir_vac = 1.0e-14; // 0.0
  }else{
    rhoil_vac = min1>min2 ? min1 : min2;
    rhoir_vac = min1>min2 ? pressureEqGasDensity(rhor,pr,rhol,pl,rhoil_vac) : 1.0e-14;
  }
  //fprintf(stdout, "rhoil_vac = %e and rhoir_vac = %e\n", rhoil_vac, rhoir_vac);

  double uil,pil,duil,dpil,uir,pir,duir,dpir;
  rarefactionJWL(-1.0, 1.0/rhol, ul, pl, 1.0/rhoil_vac, uil, pil, duil, dpil);
  double Vr[5] = { rhor, ur, 0.0, 0.0, pr};
  double cr = vf_->computeSoundSpeed(Vr,fluid1);
  rarefactionGAS(1.0, vf_->getGamma(fluid1), vf_->getGamma(fluid1)-1.0, vf_->getPressureConstant(fluid1), cr, 1.0/rhor, ur, pr, 1.0/rhoir_vac, uir, pir, duir, dpir);
  //fprintf(stdout, "uil_vac = %e and uir_vac = %e\n", uil,uir);

  vacuumValues[0] = rhoil_vac;
  vacuumValues[1] = uil;
  vacuumValues[2] = pil;
  vacuumValues[3] = rhoir_vac;
  vacuumValues[4] = uir;
  vacuumValues[5] = pir;

  if(uil<uir) return true;
  return false; //vacuumValues[0] then contains the lower bound for JWL-density
}

//----------------------------------------------------------------------------
inline
double LocalRiemannGfmpar::jwlZeroSoundSpeedJwlDensity(const double density, const double pressure)
{

    bool convergence = false;
    double tol = 1.0e-4;
    int it = 0, maxIt = 30;
    double relaxation = 1.0;

    double entropy = vf_->computeEntropy(density,pressure,fluid2);
    if(entropy >= 0.0) return -1.0;

    double xn = density;
    double xnm1 = xn;
    double dx, fn, dfn;
    double lowerBound = 0.0; // this equation may have more than one root (xn = 0 is always root)
                             // but we want the larger positive root below density
                             // so we use this lowerBound to reduce the search domain.

    while(!convergence){
        fn = (vf_->getOmega(fluid2)+1.0)*entropy*pow(xn,vf_->getOmega(fluid2)) + vf_->computeExponentials2(xn,fluid2);
        if(fn<0.0) lowerBound = xn;
        dfn = entropy*(vf_->getOmega(fluid2)+1.0)*vf_->getOmega(fluid2)*pow(xn,vf_->getOmega(fluid2)-1) + vf_->computeDerivativeOfExponentials2(xn,fluid2);
        //fprintf(stdout, "xn = %e - fn = %e - dfn = %e\n", xn, fn, dfn);
        if(dfn!=0) dx = -relaxation*fn/dfn;
        else{ dx = -0.75*dx; continue; }

        //check lower and upper bounds
        if(xn+dx>density)    dx = 0.25*relaxation*(density-xn);
        if(xn+dx<lowerBound) dx = 0.25*relaxation*(lowerBound-xn);

        if(fn<0 && dfn<0){
          dx = 0.25*(xn-xnm1);
          xn = xnm1;
        }

        if(fabs(2.0*dx/(2*xn+dx))<tol) convergence = true;

        xnm1 = xn;
        xn += dx;
        it++;
        if(it>maxIt) break;
    }

    if(!convergence){
        //if(it>maxIt) fprintf(stdout, "vacuumJwlDensity::maximum number of iteration %d reached\n", it);
        //if(dfn==0)   fprintf(stdout, "vacuumJwlDensity::zero derivative after %d iterations\n", it);
        xn = -1.0; // non-convergence value
    }//else fprintf(stdout, "vacuumJwlDensity::convergence after %d iterations\n", it);
    //fprintf(stdout, "vacuumJwlDensity::fn = %e, dfn = %e, dx = %e, xn = %e\n", fn, dfn, dx, xn);
    return xn;
}

//----------------------------------------------------------------------------
inline
double LocalRiemannGfmparGasJWL::pressureEqGasDensity(const double gasDensity, const double gasPressure,
                                                      const double jwlDensity, const double jwlPressure,
                                                      const double interfacialJwlDensity)
{
  if(interfacialJwlDensity == 0)
    return gasDensity/pow(1.0+gasPressure/vf_->getPressureConstant(fluid1),1.0/vf_->getGamma(fluid1));

  double jwlEntropy = vf_->computeEntropy(jwlDensity, jwlPressure, fluid2);
  double gasEntropy = vf_->computeEntropy(gasDensity, gasPressure, fluid1);
  return pow((jwlEntropy*pow(interfacialJwlDensity,vf_->getOmega(fluid2)+1.0)+vf_->computeExponentials(interfacialJwlDensity,fluid2)+vf_->getPressureConstant(fluid1))/gasEntropy,1.0/vf_->getGamma(fluid1));

}

//----------------------------------------------------------------------------
inline
double LocalRiemannGfmparGasJWL::sgZeroDensityPJwlDensity(const double density, const double pressure,
                                                          const double rho_c0)
{

  int verbose=0;
  if(verbose>0) fprintf(stdout, "sgZeroDensityPJwlDensity - density=%e and pressure=%e and rho_c0=%e\n", density, pressure, rho_c0);
  double entropy = vf_->computeEntropy(density,pressure,fluid2);

  double fn = entropy*pow(rho_c0>0.0 ? rho_c0 : 1.e-14,vf_->getOmega(fluid2)+1.0) + vf_->computeExponentials(rho_c0>0.0 ? rho_c0 : 1.e-14,fluid2) + vf_->getPressureConstant(fluid1);
  if(verbose>0) fprintf(stdout, "sgZeroDensityPJwlDensity - fn(max(rho_c0,0)) = %e\n", fn);
  if(fn>0.0) return -1.0;

  bool convergence = false;
  double tol = 1.0e-4;
  int it =0, maxIt = 30;
  double relaxation = 1.0;


  double xn = density;
  double xnm1 = xn;
  double dx=0, dfn;

  while(!convergence){
    fn = entropy*pow(xn,vf_->getOmega(fluid2)+1.0) + vf_->computeExponentials(xn,fluid2) + vf_->getPressureConstant(fluid1);
    dfn = entropy*(vf_->getOmega(fluid2)+1.0)*pow(xn,vf_->getOmega(fluid2)) + vf_->computeDerivativeOfExponentials(xn,fluid2);
    //fprintf(stdout, "it = %d - xn = %e - fn = %e - dfn = %e\n", it, xn, fn, dfn);
    if(dfn>0) dx = -relaxation*fn/dfn;
    
    if(xn+dx<0)       dx = -0.25*relaxation*xn;
    if(xn+dx>density) dx = 0.25*relaxation*(density-xn);

    if(fn<0 && dfn<0){
      dx = 0.5*(xn-xnm1);
      xn = xnm1;
    }

    if(fabs(2.0*dx/(2*xn+dx))<tol) convergence = true;

    xnm1 = xn;
    xn += dx;
    it++;
    if(it>maxIt) break;
  }

  if(!convergence){
    xn = -1.0;
  }
  return xn;

}

//----------------------------------------------------------------------------
inline
int LocalRiemannGfmparGasJWL::riemannInvariantGeneralTabulation(double *in, 
                                                                 double *res){
  //fprintf(stdout, "in-value = %e %e\n", in[0],in[1]);
  return sgCluster_->interpolate(1,&in,&res);

}

//----------------------------------------------------------------------------
inline
void LocalRiemannGfmparGasJWL::riemannInvariantGeneral1stOrder_wrapper(
                double *in, double *res, double *para){

  fprintf(stdout, "in[0] = %e\n", in[0]);
  fprintf(stdout, "in[1] = %e\n", in[1]);
  fprintf(stdout, "para[1] = %e\n", para[1]);
  double locin[3] = {in[0], in[1], para[1]};
  riemannInvariantGeneral1stOrder(locin, res, &(para[0]));

}
//----------------------------------------------------------------------------
inline
void LocalRiemannGfmparGasJWL::riemannInvariantGeneral2ndOrder_wrapper(
                double *in, double *res, double *para){

  double locin[3] = {in[0], in[1], para[1]};
  riemannInvariantGeneral2ndOrder(locin, res, &(para[0]));

}

//---------------------------------------------------------------------------
// Tait-JWL
//---------------------------------------------------------------------------
class LocalRiemannGfmparTaitJWL : public LocalRiemannGfmpar {

private:
  MultiFluidData::RiemannComputation riemannComputationType_;
  SparseGridCluster *sgCluster_;

  double sign;

public:
  LocalRiemannGfmparTaitJWL(VarFcn *vf, int tag1, int tag2, SparseGridCluster *sgCluster, 
			    MultiFluidData::RiemannComputation riemannComputation,double rfac,double refdensity,double refentropy,
			    MultiFluidData::TypePhaseChange typePhaseChange = MultiFluidData::RIEMANN_SOLUTION, 
			    double sgn = 1.0) : 
    LocalRiemannGfmpar(vf,tag1,tag2,typePhaseChange,rfac,refdensity,refentropy) , sign(sgn) {
    riemannComputationType_ = riemannComputation;
    sgCluster_ = sgCluster;
  }
  ~LocalRiemannGfmparTaitJWL(){ vf_ = 0; sgCluster_ = 0; }

  int computeRiemannSolution(double *Vi, double *Vj,
                              int IDi, int IDj, double *nphi,
                              double *initWi, double *initWj,
                              double *Wi, double *Wj,
                              double *rupdatei, double *rupdatej, 
                              double &weighti, double &weightj,
                              double dx[3], int it,bool isHigherOrder);

  void computeRiemannJacobian(double *Vi, double *Vj,
			      int IDi, int IDj, double *nphi,
			      double *Wi, double *Wj,
			      double dx[3],int it,
			      double* dWidWi,double*  dWidWj,
			      double* dWjdWi,double*  dWjdWj);
/*
  void eriemann(double rhol, double ul, double pl, 
                double rhor, double ur, double pr, 
                double &pi, double &ui, double &rhoil, double &rhoir){
    eriemanntj(rhol,ul,pl,rhor,ur,pr,pi,ui,rhoil,rhoir,-1.0,-1.0); }
*/
  void eriemanntj_wrapper(double *in, double *res, double *para);
  void riemannInvariantGeneral1stOrder_wrapper(
                double *in, double *res, double *para);
  void riemannInvariantGeneral2ndOrder_wrapper(
                double *in, double *res, double *para);

private:
  LocalRiemannGfmparTaitJWL();

protected:
  void eriemanntj_selector(double rhol, double ul, double pl, 
                           double rhor, double ur, double pr, 
                           double &pi, double &ui, double &rhoil, double &rhoir,
                           double initrhol, double initrhor,int& err,
                           double pcl,double pcr, double rhocl,double rhocr);
  bool eriemanntj(double rhol, double ul, double pl, 
                  double rhor, double ur, double pr, 
                  double &pi, double &ui, double &rhoil, double &rhoir,
                  double initrhol, double initrhor,int& err,
                  double pcl,double pcr, double rhocl,double rhocr);

  int riemannInvariantGeneralTabulation(double *in, double *res);
  bool vacuum(const double rhol, const double ul, const double pl,
              const double rhor, const double ur, const double pr,
              double vacuumValues[6]);
 
};

//----------------------------------------------------------------------------

inline
int LocalRiemannGfmparTaitJWL::computeRiemannSolution(double *Vi, double *Vj,
	 	int IDi, int IDj, double *nphi,
                double *initWi, double *initWj,
		double *Wi, double *Wj,
                double *rupdatei, double *rupdatej, double &weighti, double &weightj,
						       double dx[3], int it,bool isHigherOrder)
{

  int dim = 5;
	
  double P_1, P_2, U_1, U_2, R_1, R_2;
  double P_i, U_i, R_i1, R_i2;


  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  int err; 
 
  double pmin1 = vf_->getVarFcnBase(fluid1)->pmin;
  double pmin2 = vf_->getVarFcnBase(fluid2)->pmin;
  double rhomin1 = vf_->getVarFcnBase(fluid1)->rhomin;
  double rhomin2 = vf_->getVarFcnBase(fluid2)->rhomin;
 
  if (IDi==fluid1) {

    // cell i is fluid1
    // cell j is fluid2
    R_2  = Vj[0];     R_1 = Vi[0];
    U_2  = vnj;       U_1 = vni;
    P_2  = vf_->getPressure(Vj, IDj);
    P_1  = vf_->getPressure(Vi, IDi);

    double cp = vf_->specificHeatCstPressure(IDi);

    eriemanntj_selector(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,R_2,R_1,err,
              pmin2,pmin1,rhomin2,rhomin1); 
    initWi[0] = R_i1;
    initWi[1] = U_i;
    initWi[2] = P_i;
    initWj[0] = R_i2;
    initWj[1] = U_i;
    initWj[2] = P_i;

    Wi[0]  = R_i1;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    //Wi[4] = P_i;
    if (vf_->isBurnable(IDi))
      Wi[4]  = Vi[4] + 1.0/cp*(-0.5*(P_i+P_1)*(1.0/R_i1-1.0/R_1));
    else
      Wi[4]  = Vi[4] + 1.0/cp*(P_i/R_i1-P_1/R_1 - 0.5*(P_i+P_1)*(1.0/R_i1-1.0/R_1));
   /*vf_->computeTemperature(Wi, IDj);*/      Wi[dim+4]  = Wi[4];

    Wj[0]  = R_i2;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    Wj[4]  = P_i;                     Wj[dim+4]  = Wj[4];

  }else{
    // cell i is fluid2
    // cell j is fluid1
    R_2  = Vi[0];     R_1  = Vj[0];
    U_2  = vni;       U_1  = vnj;
    P_2  = vf_->getPressure(Vi, IDi);
    P_1  = vf_->getPressure(Vj, IDj);
    double cp = vf_->specificHeatCstPressure(IDj);

    eriemanntj_selector(R_2,U_2,P_2,R_1,U_1,P_1,P_i,U_i,R_i2,R_i1,R_2,R_1,err,
              pmin2,pmin1,rhomin2,rhomin1); 
    initWi[0] = R_i2;
    initWi[1] = U_i;
    initWi[2] = P_i;
    initWj[0] = R_i1;
    initWj[1] = U_i;
    initWj[2] = P_i;

    Wi[0]  = R_i2;                    Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+U_i*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+U_i*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+U_i*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = P_i;                     Wi[dim+4]  = Wi[4];

    Wj[0]  = R_i1;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+U_i*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+U_i*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+U_i*nphi[2];      Wj[dim+3]  = Wj[3];
    //Wj[4]  = P_i;
    if (vf_->isBurnable(IDj))
      Wj[4]  = Vj[4] + 1.0/cp*(-0.5*(P_i+P_1)*(1.0/R_i1-1.0/R_1));
    else
      Wj[4]  = Vj[4] + 1.0/cp*(P_i/R_i1-P_1/R_1 - 0.5*(P_i+P_1)*(1.0/R_i1-1.0/R_1));
    /*vf_->computeTemperature(Wj, IDi);*/      Wj[dim+4]  = Wj[4];
  }

// METHOD 2 : combine averaging and direction of flow
  if (it == 1 && !isHigherOrder)
    updatePhaseChangingNodeValues(dx, Wi, Wj, weighti, rupdatei, weightj, rupdatej);

  return err;
}
inline 
void LocalRiemannGfmparTaitJWL::computeRiemannJacobian(double *Vi, double *Vj,
							int IDi, int IDj, double *nphi,
							double *Wi, double *Wj,
							double dx[3],int it,
							double* dWidWi,double*  dWidWj,
							double* dWjdWi,double*  dWjdWj) {
  
	
  double P_1, P_2, U_1, U_2, R_1, R_2;
  double P_i, U_i, R_i1, R_i2;

  int dim = 5;
  int k,l;
  double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  // 3x3 Jacobians, directly from implicit riemann jacobian
  double dWidWi3[9],  dWidWj3[9], dWjdWj3[9], dWjdWi3[9];

  if (IDi==fluid1) {

    //std::cout << "ij" << std::endl << std::endl ;
    // cell i is fluid1
    // cell j is fluid2

    int err = 1;
    if (riemannComputationType_==MultiFluidData::TABULATION2) {
      double dVdv[2];
      err = rarefactionJWLderivs(-1.0, 1.0/Vj[0], vnj, Vj[4], 1.0/Wj[0] , dVdv, sgCluster_ );
      if (!err)
        ImplicitRiemann::computeTaitJwlJacobian(vf_, IDi, IDj, Vi, Vj, Wi,Wj, dWidWi3, dWidWj3,  dWjdWj3, dWjdWi3, &dVdv[0] );
    }
    if (err)
      ImplicitRiemann::computeTaitJwlJacobian(vf_, IDi, IDj, Vi, Vj, Wi,Wj, dWidWi3, dWidWj3,  dWjdWj3, dWjdWi3,NULL);

    dWidWi3[1] *= -1.0*sign;
    dWidWi3[3] *= -1.0*sign;
    dWidWi3[5] *= -1.0*sign;
    dWidWi3[7] *= -1.0*sign;
    dWidWj3[1] *= -1.0*sign;
    dWidWj3[3] *= -1.0*sign;
    dWidWj3[5] *= -1.0*sign;
    dWidWj3[7] *= -1.0*sign;
    dWjdWi3[1] *= -1.0*sign;
    dWjdWi3[3] *= -1.0*sign;
    dWjdWi3[5] *= -1.0*sign;
    dWjdWi3[7] *= -1.0*sign;
    dWjdWj3[1] *= -1.0*sign;
    dWjdWj3[3] *= -1.0*sign;
    dWjdWj3[5] *= -1.0*sign;
    dWjdWj3[7] *= -1.0*sign;

  }else{
    // cell i is fluid2
    // cell j is fluid1

    //std::cout << "ji" << std::endl << std::endl ;
    int err = 1;
    if (riemannComputationType_==MultiFluidData::TABULATION2) {
      double dVdv[2];
      err = rarefactionJWLderivs(-1.0, 1.0/Vi[0], vni, Vi[4], 1.0/Wi[0] , dVdv, sgCluster_ );
      if (!err)
        ImplicitRiemann::computeTaitJwlJacobian(vf_, IDj, IDi, Vj, Vi, Wj,Wi, dWjdWj3, dWjdWi3,  dWidWi3, dWidWj3, &dVdv[0] );
    }
    if (err)
      ImplicitRiemann::computeTaitJwlJacobian(vf_, IDj, IDi, Vj, Vi, Wj,Wi, dWjdWj3, dWjdWi3,  dWidWi3, dWidWj3,NULL);

    dWidWi3[1] *= -1.0*sign;
    dWidWi3[3] *= -1.0*sign;
    dWidWi3[5] *= -1.0*sign;
    dWidWi3[7] *= -1.0*sign;
    dWidWj3[1] *= -1.0*sign;
    dWidWj3[3] *= -1.0*sign;
    dWidWj3[5] *= -1.0*sign;
    dWidWj3[7] *= -1.0*sign;
    dWjdWi3[1] *= -1.0*sign;
    dWjdWi3[3] *= -1.0*sign;
    dWjdWi3[5] *= -1.0*sign;
    dWjdWi3[7] *= -1.0*sign;
    dWjdWj3[1] *= -1.0*sign;
    dWjdWj3[3] *= -1.0*sign;
    dWjdWj3[5] *= -1.0*sign;
    dWjdWj3[7] *= -1.0*sign;
  }

  this->oneDtoThreeD(dWidWi3, dWidWj3,
		     dWjdWi3,dWjdWj3,nphi,
		     dWidWi, dWidWj,
		     dWjdWi, dWjdWj);
}

//----------------------------------------------------------------------------
inline
void LocalRiemannGfmparTaitJWL::eriemanntj_wrapper(
                               double *in, double *res, double *para)
{

  double dummy1, dummy2;
  int err;
  eriemanntj(in[0], 0.0, in[1], in[2], in[4], in[3], dummy1, dummy2, res[0], res[1], -1.0, -1.0,
             err, -1.0, -1.0, -1.0, -1.0);

}

//----------------------------------------------------------------------------
inline
void LocalRiemannGfmparTaitJWL::eriemanntj_selector(
                               double rhol, double ul, double pl, 
                               double rhor, double ur, double pr, 
                               double &pi, double &ui,  
                               double &rhoil, double &rhoir,
                               double initrhol, double initrhor,int& err,
                               double pcl,double pcr, double rhocl,double rhocr)
{

  err = 0;

  if(riemannComputationType_==MultiFluidData::TABULATION5){
    double *in  = new double[5]; in[0]=rhol;in[1]=pl;in[2]=rhor;in[3]=pr;in[4]=ur-ul;
    double *res = new double[2]; res[0]=0.0;res[1]=0.0;
    sgCluster_->interpolate(1,&in,&res);
    rhoil = fmax(res[0],0.0); rhoir = fmax(res[1],0.0);
    double d[2]; //dummy variable
    double uir, pir, uil, pil;

    double omegal = vf_->getOmega(fluid2);
    double omp1ooml = (omegal+1.0)/omegal;
    double frhol = vf_->computeFrho(rhol,fluid2);
    double alphar = vf_->getAlphaWater(fluid1);
    double betar = vf_->getBetaWater(fluid1);
    double pinfr = vf_->getPrefWater(fluid1);
    double Vr[5] = { rhor, ur, 0.0, 0.0, pr };
    double cr = vf_->computeSoundSpeed(Vr,fluid1);

    if( rhoil > rhol){
      double frhoil  = vf_->computeFrho(rhoil,fluid2);
      double frhopil = vf_->computeFrhop(rhoil,fluid2);
      shockJWL(-1.0, omegal, omp1ooml, frhol, frhoil, frhopil, 1.0/rhol, ul, pl, 1.0/rhoil, uil, pil, d[0], d[1]);
    }else
      rarefactionJWL(-1.0, 1.0/rhol, ul, pl, 1.0/rhoil, uil, pil, d[0], d[1], riemannComputationType_,1);
    if( rhoir > rhor)
      shockTAIT(1.0, alphar,betar,pinfr, 1.0/rhor, ur, pr, 1.0/rhoir, uir, pir, d[0], d[1],0);
    else
      rarefactionTAIT(1.0, alphar,betar,pinfr, 1.0/rhor, ur, pr, 1.0/rhoir, uir, pir, d[0], d[1],0);

    ui = 0.5*(uil+uir);
    pi = 0.5*(pil+pir);
    pi = std::max<double>(pi, std::max<double>(pcl,pcr));
    rhoil = std::max<double>(rhoil, rhocl);
    rhoir = std::max<double>(rhoir, rhocr);
    // not checking for vacuum!

  /*double pcut = std::max<double>(  pcl,   pcr);
    double rcut = std::max<double>(rhocl, rhocr);

    if(pi <= pcut || rhoil <= rcut || rhoir <= rcut){
      std::cout << "*** ERROR ERS GJ: detected too small density or pressure " << std::endl;
      std::cout << " rL, rR, Pi " << std::endl;
      std::cout << rhoil << " " << rhoir << " " << " " << pi << std::endl;
      err = 1;
    }*/

  }
  else
    eriemanntj(rhol,ul,pl,rhor,ur,pr,pi,ui,rhoil,rhoir,initrhol,initrhor, err,pcl,pcr,rhocl,rhocr);

}

//----------------------------------------------------------------------------
inline
bool LocalRiemannGfmparTaitJWL::eriemanntj(double rhol, double ul, double pl, 
                                          double rhor, double ur, double pr, 
                                          double &pi, double &ui,  
                                          double &rhoil, double &rhoir,
                                          double initrhol, double initrhor,int& err,
                                          double pcl,double pcr, double rhocl,double rhocr){
// left  -- JWL -- phi = -1.0
// right -- GAS -- phi = +1.0
  int verbose = -1;
  if(verbose>0){
    fprintf(stdout, "---- new Riemann ----\n");
    fprintf(stdout, "initial rhoil, rhoir = %e %e\n", rhol, rhor);
    fprintf(stdout, "initial vil,   vir   = %e %e\n", 1.0/rhol, 1.0/rhor);
  }

//initialize
  double uil, uir, pil, pir, duil, duir, dpil, dpir;
  double jacobian[4];/* uil, uir, pil, pir*/
  double function[2];
  double increment[2];
  bool convergence = false;
  double eps = 1.e-3;
  int MaxIts = 100;
  int it = 0;
  double relaxationFactorJwl = relaxFactorJwl; //0.85; // must be between 0 and 1
  double relaxationFactorGas = relaxFactorJwl; //0.85; // must be between 0 and 1
  int count = 0;

  double pcut = std::max<double>(  pcl,   pcr);
  double rcut = std::max<double>(rhocl, rhocr);

  double vl  = 1.0/rhol;
  double vr  = 1.0/rhor;
  double vil = vl;
  double vir = vr;
  double pmin = max( vf_->getPmin(fluid1), vf_->getPmin(fluid2) );
  vil = initrhol>0.0 ? 1.0/initrhol : vl;
  vir = initrhor>0.0 ? 1.0/initrhor : vr;

  double omegal = vf_->getOmega(fluid2);
  double omp1ooml = (omegal+1.0)/omegal;
  double frhol = vf_->computeFrho(1.0/vl,fluid2);
  double frhoil = frhol;
  double frhopil = vf_->computeFrhop(1.0/vl,fluid2);


  double alphar = vf_->getAlphaWater(fluid1);
  double betar = vf_->getBetaWater(fluid1);
  double pinfr = vf_->getPrefWater(fluid1);
  double Vr[5] = { 1.0/vr, ur, 0.0, 0.0, pr };
  double cr = vf_->computeSoundSpeed(Vr,fluid1);

  double res = 1.0e20;

  err = 0;

  double vacuumValues[6]; /* rhoil, uil, pil, rhoir, uir, pir */
  vacuumValues[0] = -1.0; // positive if proper vacuum values are computed
  bool checkVacuumValues = false;
  if(checkVacuumValues){
    if(vacuum(rhol,ul,pl,rhor,ur,pr,vacuumValues)){
      if(verbose>-1){
        fprintf(stdout, "rhoil_vac = %e and rhoir_vac = %e\n", vacuumValues[0], vacuumValues[3]);
        fprintf(stdout, "uil_vac   = %e and uir_vac   = %e\n", vacuumValues[1], vacuumValues[4]);
        fprintf(stdout, "pil_vac   = %e and pir_vac   = %e\n", vacuumValues[2], vacuumValues[5]);
      }
      rhoil = vacuumValues[0];
      rhoir = vacuumValues[3];
      ui    = 0.5*(vacuumValues[1]+vacuumValues[4]);
      pi    = 0.5*(vacuumValues[2]+vacuumValues[5]);
      return true;
    }
    if(verbose>4) fprintf(stdout, "checking vacuum possibilities -- DONE\n");
  }else{
    if(verbose>4) fprintf(stdout, "no checking of vacuum possibilities\n");
  }

//start newton iteration loop
  while(!convergence){
    if(verbose>0) fprintf(stdout, "\n");

  //compute left  JWL-term (shock or rarefaction)
    if( vil < vl){
      if(verbose>0) fprintf(stdout, "shockJWL\n");
      frhoil  = vf_->computeFrho(1.0/vil,fluid2);
      frhopil = vf_->computeFrhop(1.0/vil,fluid2);
      shockJWL(-1.0, omegal, omp1ooml, frhol, frhoil, frhopil, vl, ul, pl, vil, uil, pil, duil, dpil);
    }else{
      if(verbose>0) fprintf(stdout, "rarefactionJWL\n");
      rarefactionJWL(-1.0, vl, ul, pl, vil, uil, pil, duil, dpil, riemannComputationType_,1);
    }
  //compute right GAS-term (shock or rarefaction)
    if( vir < vr){
      if(verbose>0) fprintf(stdout, "shockTAIT\n");
      shockTAIT(1.0, alphar,betar,pinfr, vr, ur, pr, vir, uir, pir, duir, dpir,0);
    }
    else{
      if(verbose>0) fprintf(stdout, "rarefactionTAIT\n");
      rarefactionTAIT(1.0, alphar,betar,pinfr, vr, ur, pr, vir, uir, pir, duir, dpir,0);
    }

    if(verbose>1){
      fprintf(stdout, "uil  = %e and uir  = %e\n", uil, uir);
      fprintf(stdout, "pil  = %e and pir  = %e\n", pil, pir);
      fprintf(stdout, "duil = %e and duir = %e\n", duil, duir);
      fprintf(stdout, "dpil = %e and dpir = %e\n", dpil, dpir);
    }

    //solve2x2System: function = jacobian*increment
    function[0] = uil-uir;
    function[1] = pil-pir;
    jacobian[0] = duil; jacobian[1] = -duir;
    jacobian[2] = dpil; jacobian[3] = -dpir;
    increment[0] = 0.0; increment[1] = 0.0;

    /*double res2 = fabs(function[1])+fabs(function[0]);
    if (res2 < res)
      res = res2;
    else {
      vil += increment[0]*0.5;
      vir += increment[1]*0.5;
      increment[0]*=0.5;
      increment[1]*=0.5;
      continue;
      }*/
    
    bool solved = solve2x2System(jacobian,function,increment);
    if(!solved){
      fprintf(stdout, "$$$$\n");
      fprintf(stdout, "rhol, ul, pl = %e %e %e\n", rhol, ul, pl);
      fprintf(stdout, "rhor, ur, pr = %e %e %e\n", rhor, ur, pr);
      fprintf(stdout, "rhoil  = %e and rhoir  = %e\n", 1/vil, 1/vir);
      fprintf(stdout, "uil  = %e and uir  = %e\n", uil, uir);
      fprintf(stdout, "pil  = %e and pir  = %e\n", pil, pir);
      fprintf(stdout, "duil = %e and duir = %e\n", duil, duir);
      fprintf(stdout, "dpil = %e and dpir = %e\n", dpil, dpir);
      rarefactionTAIT(1.0, alphar,betar,pinfr, vr, ur, pr, vir, uir, pir, duir, dpir,1);
    }

    if(verbose>2) fprintf(stdout, "dvil = %e and dvir = %e\n", -increment[0],-increment[1]);

    //update values and check bounds

    if(verbose>3) fprintf(stdout, "increment/v = %e %e\n", increment[0]/vil, increment[1]/vir);

    // prevent large increases
    if(-increment[0]>2.0*vil) increment[0] = -2.0*vil;
    if(-increment[1]>2.0*vir) increment[1] = -2.0*vir;
    // prevent large decreases
    if(increment[0]>0.5*vil) increment[0] = 0.5*vil;
    if(increment[1]>0.5*vir) increment[1] = 0.5*vir;

    increment[0] *= relaxationFactorJwl;
    increment[1] *= relaxationFactorGas;
    
    vil -= increment[0];
    vir -= increment[1];
    if(verbose>2) fprintf(stdout, "2 -- vil = %e and vir = %e\n", vil, vir);

    if(vil < vl){ // at next iteration, leftrarefaction => ensures that some conditions are fulfilled
      double temp = omegal*vl/(omegal+2.0);
      if(vil<temp){
        vil += increment[0];
        vir += increment[1];
        double alpha = -0.5*(temp-vil)/increment[0];
        increment[0] *= alpha;
        increment[1] *= alpha;
        vil -= increment[0];
        vir -= increment[1];
        count++;
      }
    }
    if(vir < vr){ // at next iteration, rightrarefaction => ensures that some conditions are fulfilled
      /*fprintf(stdout,"??\n");
      fprintf(stdout, "rhol, ul, pl = %e %e %e\n", rhol, ul, pl);
      fprintf(stdout, "rhor, ur, pr = %e %e %e\n", rhor, ur, pr);
      fprintf(stdout, "initrhol, initrhor = %e %e\n", initrhol, initrhor);
      fprintf(stdout, "rhoil  = %e and rhoir  = %e\n", 1/vil, 1/vir);
      fprintf(stdout, "uil  = %e and uir  = %e\n", uil, uir);
      fprintf(stdout, "pil  = %e and pir  = %e\n", pil, pir);
      fprintf(stdout, "duil = %e and duir = %e\n", duil, duir);
      fprintf(stdout, "dpil = %e and dpir = %e\n", dpil, dpir);
      double pit = pinfr+alphar*pow(vir, -betar);
      fprintf(stdout,"pit = %e pmin = %e\n",pit,pmin);
      if (pit < pmin) {

	vir = pow( (pmin-pinfr) / alphar, betar);
      }

      double Pt,Ut,Rigt,Riwt;
      solveSGTait(rhol,ul,pl, 
		  rhor,ur,pr,
		  Pt,Ut, 
		  Rigt, Riwt,
		  alphar,betar,
		  pinfr, omegal+1.0,0.0);
      fprintf(stdout, "SGTAIT: rhoil  = %e and rhoir  = %e\n", Rigt, Riwt);
      fprintf(stdout, "SGTAIT: uil  = %e and uir  = %e\n", Ut, Ut);
      fprintf(stdout, "SGTAIT: pil  = %e and pir  = %e\n", Pt,Pt);*/
      /*double temp = (gamr-1.0)/(gamr+1.0)*vr;
      if(vir<temp){
        vil += increment[0];
        vir += increment[1];
        double alpha = -0.5*(temp-vir)/increment[1];
        increment[0] *= alpha;
        increment[1] *= alpha;
        vil -= increment[0];
        vir -= increment[1];
        count++;
	}*/
      }
    if(verbose>2) fprintf(stdout, "3 -- vil = %e and vir = %e\n", vil, vir);
    //check - in case of rarefaction at next iteration, 1.0/rhoil may not be above a certain value
    /*if(vacuumValues[0] > 0.0 && vil>1.0/vacuumValues[0]){ // vacuumValues is negative if it does not contain any proper value(see declaration and definition above)
      vil += increment[0];
      increment[0] = -0.5*(1.0/vacuumValues[0] - vil);
      vil -= increment[0];
      }*/
    if(verbose>2) fprintf(stdout, "4 -- vil = %e and vir = %e\n", vil, vir);
    if(verbose>0) fprintf(stdout, "rhoil = %e and rhoir = %e\n", 1.0/vil, 1.0/vir);
    it++;

    if (1.0/vil < rhocl)
      vil = 1.0/rhocl;
    if (1.0/vir < rhocr)
      vir = 1.0/rhocr;

    pi = alphar*pow(1.0/vir,betar)+pinfr;
    if (pi < pcut) {

      vir = pow((pi-pinfr)/alphar,1.0/betar);
    } 

  //check convergence criterion
    if(fabs(increment[0])<eps*fabs(vil) &&
       fabs(increment[1])<eps*fabs(vir) )
      convergence = true;
    if(it>MaxIts) break;


  }//end newton iteration loop

  if( vil < vl){
    frhoil  = vf_->computeFrho(1.0/vil,fluid2);
    frhopil = vf_->computeFrhop(1.0/vil,fluid2);
    shockJWL(-1.0, omegal, omp1ooml, frhol, frhoil, frhopil, vl, ul, pl, vil, uil, pil, duil, dpil);
  }else rarefactionJWL(-1.0, vl, ul, pl, vil, uil, pil, duil, dpil, riemannComputationType_,0);
  
  if( vir < vr) shockTAIT(1.0, alphar,betar,pinfr, vr, ur, pr, vir, uir, pir, duir, dpir,0);
  else rarefactionTAIT(1.0, alphar,betar,pinfr, vr, ur, pr, vir, uir, pir, duir, dpir,0);

  rhoil = 1.0/vil;
  rhoir = 1.0/vir;
  ui    = 0.5*(uil+uir);
  pi    = 0.5*(pil+pir);
  pi = std::max<double>(pcut, pi);

  if (pi == pcut) {

    rhoir = pow((pi-pinfr)/alphar,1.0/betar);
  }

  if(convergence){
    if(verbose>-1) fprintf(stdout, "riemann has converged to an approximate solution in %d iterations\n", it);
  }else{
    if(verbose>-1) fprintf(stdout, "riemann solver did not converged\n");
    if(verbose>-1) fprintf(stdout, "Warning: solution will be state given by vacuum\n");
    /*    rhoil = vacuumValues[0];
    rhoir = vacuumValues[3];
    uil   = vacuumValues[1];
    uir   = vacuumValues[4];
    pil   = vacuumValues[2];
    pir   = vacuumValues[5];
    ui    = 0.5*(uil+uir);
    pi    = 0.5*(pil+pir);*/
    err = 1;
    if(verbose>-1) fprintf(stdout, "Warning: uil = %e and uir = %e\n", uil, uir);
  }
  
  if(verbose>-1){
    fprintf(stdout, "rhol, ul, pl = %e %e %e\n", rhol, ul, pl);
    fprintf(stdout, "rhor, ur, pr = %e %e %e\n", rhor, ur, pr);
    fprintf(stdout, "initrhol, initrhor = %e %e\n", initrhol, initrhor);
    fprintf(stdout, "rhoil  = %e and rhoir  = %e\n", 1/vil, 1/vir);
    fprintf(stdout, "uil  = %e and uir  = %e\n", uil, uir);
    fprintf(stdout, "pil  = %e and pir  = %e\n", pil, pir);
    fprintf(stdout, "duil = %e and duir = %e\n", duil, duir);
    fprintf(stdout, "dpil = %e and dpir = %e\n", dpil, dpir);
  }

  /*if(pi <= pcut || rhoil <= rcut || rhoir <= rcut){
    std::cout << "*** ERROR ERS WJ: detected too small density or pressure " << std::endl;
    std::cout << " rL, rR, Pi " << std::endl;
    std::cout << rhoil << " " << rhoir << " " << " " << pi << std::endl; 
    err = 1;
  }*/

  if(convergence) return true;
  else            return false;


}

inline
bool LocalRiemannGfmparTaitJWL::vacuum(const double rhol, const double ul, const double pl,
                                      const double rhor, const double ur, const double pr,
                                      double vacuumValues[6]){
// notation: JWL on the left and Tait on the right
// remember vacuum can occur only between two rarefaction waves, thus decrease of densities

// 1st step: find JWL-density for which there is loss of positivity of c^2 in JWL gas
  double min1 = jwlZeroSoundSpeedJwlDensity(rhol,pl); // returns -1 if none found

  double rhoil_vac, rhoir_vac;
  if(min1 < 0){
    rhoil_vac = 1.0e-14; // 0.0
  } else {
    rhoil_vac = min1;

  }
  rhoir_vac = pow(-vf_->getPrefWater(fluid1)/vf_->getAlphaWater(fluid1), 
                  1.0/vf_->getBetaWater(fluid1))*(1.0+1.0e-8);

  double uil,pil,duil,dpil,uir,pir,duir,dpir;
  rarefactionJWL(-1.0, 1.0/rhol, ul, pl, 1.0/rhoil_vac, uil, pil, duil, dpil);
  double Vr[5] = { rhor, ur, 0.0, 0.0, pr};
  double cr = vf_->computeSoundSpeed(Vr,fluid1);
  rarefactionTAIT(1.0, vf_->getAlphaWater(fluid1), vf_->getBetaWater(fluid1), vf_->getPrefWater(fluid1),
                  1.0/rhor, ur, pr, 1.0/rhoir_vac, uir, pir, duir, dpir,0);
  //fprintf(stdout, "uil_vac = %e and uir_vac = %e\n", uil,uir);

  vacuumValues[0] = rhoil_vac;
  vacuumValues[1] = uil;
  vacuumValues[2] = pil;
  vacuumValues[3] = rhoir_vac;
  vacuumValues[4] = uir;
  vacuumValues[5] = pir;

  if(uil<uir) return true;
  return false; //vacuumValues[0] then contains the lower bound for JWL-density
}



//----------------------------------------------------------------------------
inline
int LocalRiemannGfmparTaitJWL::riemannInvariantGeneralTabulation(double *in, 
                                                                 double *res){
  //fprintf(stdout, "in-value = %e %e\n", in[0],in[1]);
  return sgCluster_->interpolate(1,&in,&res);

}

//----------------------------------------------------------------------------
inline
void LocalRiemannGfmparTaitJWL::riemannInvariantGeneral1stOrder_wrapper(
                double *in, double *res, double *para){

  fprintf(stdout, "in[0] = %e\n", in[0]);
  fprintf(stdout, "in[1] = %e\n", in[1]);
  fprintf(stdout, "para[1] = %e\n", para[1]);
  double locin[3] = {in[0], in[1], para[1]};
  riemannInvariantGeneral1stOrder(locin, res, &(para[0]));

}
//----------------------------------------------------------------------------
inline
void LocalRiemannGfmparTaitJWL::riemannInvariantGeneral2ndOrder_wrapper(
                double *in, double *res, double *para){

  double locin[3] = {in[0], in[1], para[1]};
  riemannInvariantGeneral2ndOrder(locin, res, &(para[0]));

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Adam 2010.08.17
// This Local Riemann is now templated by dim. It can now handle Turbulent simulations
template<int dim>
class LocalRiemannFluidStructure : public LocalRiemann {

public:
  LocalRiemannFluidStructure() : LocalRiemann(), stabil_alpha(0.0), viscous_switch(0.0), prec(false) {fluid1 = fluid2 = 0;}
  LocalRiemannFluidStructure(VarFcn *vf) : LocalRiemann(vf,0,0), stabil_alpha(0.0), viscous_switch(0.0), prec(false) {fluid1 = fluid2 = 0;}
  virtual ~LocalRiemannFluidStructure() { vf_ = 0; }

  void setStabilAlpha(double a) { stabil_alpha = a; }
  void setViscousSwitch(double v) { viscous_switch = v; }

  void setPreconditioner(double beta) { prec = true; mach = beta; }

  int computeRiemannSolution(double *Vi, double *Vstar,
                            double *nphi, VarFcn *vf,
                            double *Wstar, double *rupdatei,
                            double &weighti, int it, int Id = 0);
  void computeRiemannSolution(int tag, double *Vi, double *Vstar,
                            double *nphi, VarFcn *vf,
                            double *Wstar, double *rupdatei,
                            double &weighti, int it);//TODO:not needed!

  void computeRiemannJacobian(double *Vi, double *Vstar,
                            double *nphi, VarFcn *vf,
                            double *Wstar, double *rupdatei,
                            double &weighti, int it, double* WstardU,int Id = 0);

  void computeRiemannderivative(double *Vi, double *Vstar,
				double *nphi, VarFcn *vf,
				double *Wstar, double* dWstardn, int Id);

  // Multi-Phase Riemann solvers (implemented here just to stop compiler's complaining...)
  int computeRiemannSolution(double *Vi, double *Vj,
                            int IDi, int IDj, double *nphi,
                            double *initWi, double *initWj,
                            double *Wi, double *Wj,
                            double *rupdatei, double *rupdatej,
                            double &weighti, double &weightj,
                            double dx[3], int it) {
    fprintf(stderr,"ERROR: Should not call the two-phase Riemann solver in LocalRiemannFluidStructure!\n"); return 0; }

  void computeRiemannJacobian(double *Vi, double *Vj,
                                      int IDi, int IDj, double *nphi,
                                      double *Wi, double *Wj,
                                      double dx[3],int it,
                                      double* dWidWi,double*  dWidWj,
                                      double* dWjdWi,double*  dWjdWj) {
    fprintf(stderr,"ERROR: Should not call the two-phase Riemann solver in LocalRiemannFluidStructure!\n");}


private:

  bool prec;
  double mach;

  int eriemannfs(double rhol, double ul, double pl,
		 double &rhoi, double ui, double &pi,
		 VarFcn *vf, int Id,int& err,double pc, double rc,
 		 bool prec, double mach); //note: ui shouldn't be changed. so the value (instead of reference) is used.

  void eriemannfs_grad(double rho, double u, double p,
                       double &rhoi, double ui, double &pi,
                       VarFcn *vf, double* dWdWi,int Id,
                       bool prec, double mach, double& drdus); //Caution: "ui" will not be modified!
    
  void eriemannfs_tait(double rhol, double ul, double pl,
                  double &rhoi, double ui, double &pi,
                  VarFcn *vf, int Id,int& err,double pc, double rc); //note: ui shouldn't be changed. so the value (instead of reference) is used.

  void eriemannfs_tait_grad(double rho, double u, double p,
                       double &rhoi, double ui, double &pi,
                       VarFcn *vf, double* dWdWi,int Id); //Caution: "ui" will not be modified!

  template <int d>
  friend class FSJac;

  double stabil_alpha;
  double viscous_switch;
};

//------------------------------------------------------------------------------

template<int dim>
inline
int LocalRiemannFluidStructure<dim>::computeRiemannSolution(double *Vi, double *Vstar,
                            double *nphi, VarFcn *vf,
                            double *Wstar, double *rupdatej,
                            double &weightj, int it, int Id)
{
  // Commented by Adam on 2010.08.17 because it has to handle dim > 5
  //int dim = 5;

  double P_1, U_1, R_1; // pass to 1D-FSI Riemann solver
  double P_i, U_i, R_i; // solution given by 1D-FSI Riemann solver

  //---------------------------------------------------------------

  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};


  double rc = vf->getVarFcnBase(Id)->rhomin;
  double pc = vf->getVarFcnBase(Id)->pmin;

  R_1 = std::max(rc,Vi[0]);
  U_1 = vni;
  P_1  = std::max(pc,vf->getPressure(Vi,Id));
  U_i = Vstar[0]*nphi[0]+Vstar[1]*nphi[1]+Vstar[2]*nphi[2];
  double U_ti[3] = {Vstar[0] - U_i*nphi[0], Vstar[1] - U_i*nphi[1], Vstar[2] - U_i*nphi[2]};

/*
  double normv = Vi[1]*Vi[1]+Vi[2]*Vi[2]+Vi[3]*Vi[3];
  // Attempt at stabilization of structure normal.
  U_1 += stabil_alpha*sqrt(normv - U_1*U_1);
*/
  int err;

  switch (vf->getType(Id)) {
  case VarFcnBase::STIFFENEDGAS:
  case VarFcnBase::PERFECTGAS:
    eriemannfs(R_1,U_1,P_1,R_i,U_i,P_i,vf,Id,err,pc,rc, prec, mach); //caution: U_i will not be modified!
    break;
  case VarFcnBase::TAIT:
    eriemannfs_tait(R_1,U_1,P_1,R_i,U_i,P_i,vf,Id,err,pc,rc); //caution: U_i will not be modified!
    break;
  }

  if (err)
    return err;

  Wstar[0]  = R_i;
  Wstar[1]  = U_i*nphi[0] + viscous_switch*U_ti[0] + (1.0-viscous_switch)*(1.0-stabil_alpha)*vti[0];
  Wstar[2]  = U_i*nphi[1] + viscous_switch*U_ti[1] + (1.0-viscous_switch)*(1.0-stabil_alpha)*vti[1];
  Wstar[3]  = U_i*nphi[2] + viscous_switch*U_ti[2] + (1.0-viscous_switch)*(1.0-stabil_alpha)*vti[2];
  if (vf->getType(Id) == VarFcnBase::TAIT)
    Wstar[4] = vf->computeTemperature(Vi, Id);
  else
    Wstar[4]  = P_i;
  if(dim == 6)
    {
      Wstar[5]  = 0.0;// Boundary Condition: nuTilde = 0
    }
  else if(dim == 7)  // Boundary Condition for KE. To be improved with Wall Function...
    {
      Wstar[5] = 0.0;
      Wstar[6] = 0.0;
    }


  //---------------------------------------------------------------

  double Vi0[dim];
  for(int i=0; i<dim; i++)
    Vi0[i] = Vi[dim+i];

  vni = Vi0[1]*nphi[0]+Vi0[2]*nphi[1]+Vi0[3]*nphi[2];
  vti[0] = Vi0[1] - vni*nphi[0];
  vti[1] = Vi0[2] - vni*nphi[1];
  vti[2] = Vi0[3] - vni*nphi[2];

  R_1 = std::max(rc,Vi0[0]);
  U_1 = vni;
  P_1 = std::max(pc,vf->getPressure(Vi0,Id));

  // U_i is the same.
  switch (vf->getType(Id)) {
  case VarFcnBase::STIFFENEDGAS:
  case VarFcnBase::PERFECTGAS:
    eriemannfs(R_1,U_1,P_1,R_i,U_i,P_i,vf,Id,err,pc,rc, prec, mach); //caution: U_i will not be modified!
    break;
  case VarFcnBase::TAIT:
    eriemannfs_tait(R_1,U_1,P_1,R_i,U_i,P_i,vf,Id,err,pc,rc); //caution: U_i will not be modified!
    break;
  }

  Wstar[dim]    = R_i;
  Wstar[dim+1]  = U_i*nphi[0] + viscous_switch*U_ti[0] + (1.0-viscous_switch)*(1.0-stabil_alpha)*vti[0];
  Wstar[dim+2]  = U_i*nphi[1] + viscous_switch*U_ti[1] + (1.0-viscous_switch)*(1.0-stabil_alpha)*vti[1];
  Wstar[dim+3]  = U_i*nphi[2] + viscous_switch*U_ti[2] + (1.0-viscous_switch)*(1.0-stabil_alpha)*vti[2];
  if (vf->getType(Id) == VarFcnBase::TAIT) 
    Wstar[dim+4] = vf->computeTemperature(Vi0, Id);
  else
    Wstar[dim+4]  = P_i;
  if(dim == 6)
    {
      Wstar[dim+5]  = 0.0; // Boundary Condition: nuTilde = 0
    }
  else if(dim == 7) // Boundary Condition for KE. To be improved with Wall Function...
    {
      Wstar[dim+5]  = 0.0; 
      Wstar[dim+6]  = 0.0; 
    }

  //-----------------------------------------------------------------
/*
  if(it==1){
    weightj += 1.0;
    for (int k=0; k<dim; k++)
      rupdatej[k] += Wstar[k];  //TODO: rupdate is never used for FSI. (only used for MPF)
  }
*/

  return err;
}
/*
template<int dim>
class FSJac {
 
  int Id;
  VarFcn* vf;
  double U_i;
  LocalRiemannFluidStructure<dim>* ls; 
  public:
  FSJac(LocalRiemannFluidStructure<dim>* _ls,VarFcn* _vf, int _Id,double _ui) : vf(_vf), Id(_Id),ls(_ls),U_i(_ui) {}
  void Compute(const double u[3],double f[3]) const {

    ls->eriemannfs_tait(u[0],u[1],u[2],f[0],U_i,f[2],vf,Id); //caution: U_i will not be modified!
  }
};
*/

//------------------------------------------------------------------------------
template<int dim>
inline
void LocalRiemannFluidStructure<dim>::computeRiemannJacobian(double *Vi, double *Vstar,
                                      double *nphi, VarFcn *vf,
                                      double *Wstar, double *rupdatei,
                                      double &weighti, int it, double* dWstardU,int Id) {

  double P_1, U_1,R_1; // pass to 1D-FSI Riemann solver
  double P_i=Wstar[4], U_i, R_i=Wstar[0]; // solution given by 1D-FSI Riemann solver
  //---------------------------------------------------------------

  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

  double dWdW[9]={1,0,0,0,1,0,0,0,1};
  double dummy;

  R_1 = Vi[0];
  U_1 = vni;
  P_1  = vf->getPressure(Vi,Id);
  U_i = Vstar[0]*nphi[0]+Vstar[1]*nphi[1]+Vstar[2]*nphi[2];
  P_i = vf->getPressure(Wstar,Id);
/*
  double normv = Vi[1]*Vi[1]+Vi[2]*Vi[2]+Vi[3]*Vi[3];

  // Attempt at stabilization of structure normal.
  U_1 += stabil_alpha*sqrt(normv - U_1*U_1);
*/
  switch (vf->getType(Id)) {
  case VarFcnBase::STIFFENEDGAS:
  case VarFcnBase::PERFECTGAS:
    eriemannfs_grad(R_1,U_1,P_1,R_i,U_i,P_i,vf,dWdW,Id, prec, mach, dummy); //caution: U_i will not be modified!
    break;
  case VarFcnBase::TAIT:
    eriemannfs_tait_grad(R_1,U_1,P_1,R_i,U_i,P_i,vf,dWdW,Id);
    break;
  }

  // Checking jacobian
  double u[3] = {R_1,U_1,P_1};
  double f[3] = {R_i,U_i,P_i};
  //std::cout << R_1 << " " << U_1 << " " << P_1 << " " << R_i << " " << U_i << " " << P_i << std::endl;
  //if (Id == 0)
  //  DebugTools::CheckJacobian<3>(dWdW, u,f,FSJac<dim>(this,vf,Id,U_i), "FSJacobian\n-------------------------\n"); 

  if(dim == 6)
    {
      Wstar[5]  = 0.0;// Boundary Condition: nuTilde = 0
    }
  else if(dim == 7) // Boundary Condition for KE. To be improved with Wall Function...
    {
      Wstar[5] = 0.0;
      Wstar[6] = 0.0;
    }
 
  memset(dWstardU, 0, sizeof(double)*dim*dim);
  dWstardU[0] = dWdW[0];
  dWstardU[dim*4+4] = dWdW[8];

  for (int i = 0; i < 3; ++i) {
    dWstardU[(i+1)] = nphi[i]*dWdW[1];
    dWstardU[(i+1)+4*dim] = nphi[i]*dWdW[7];
  }

  for (int i = 0; i < 3; ++i) {
 
    for (int j = 0; j < 3; ++j) {

      dWstardU[(i+1)*dim+(j+1)] = ((i==j?1.0:0.0) - nphi[i]*nphi[j])*(1.0-viscous_switch)*(1.0-stabil_alpha);
    }
  }

  dWstardU[4] = dWdW[2];
  dWstardU[dim*4] = dWdW[6]; 
}
//------------------------------------------------------------------------------

template<int dim>
inline
void LocalRiemannFluidStructure<dim>::computeRiemannderivative(double *Vi, double *Vstar,
							       double *nphi, VarFcn *vf,
							       double *Wstar, double* dWstardn, int Id) {

  // Compute the derivative of the WStar w.r.t the normal vector (nphi) //
  // Works only for with perfect gas and no low-Mach preconditioner!!!!!!
  
  double P_1, U_1, R_1; // pass to 1D-FSI Riemann solver
  double P_i, U_i, R_i; // solution given by 1D-FSI Riemann solver
  double drdus, dpdus;

  double vni = Vi[1]*nphi[0] + Vi[2]*nphi[1] + Vi[3]*nphi[2];

  double dWdW[9]={1,0,0,0,1,0,0,0,1};

  R_i = Wstar[0];
  P_i = Wstar[4];

  R_1 = Vi[0];
  U_1 = vni;
  P_1 = vf->getPressure(Vi,Id);

  U_i = Vstar[0]*nphi[0] + Vstar[1]*nphi[1] + Vstar[2]*nphi[2];
  P_i = vf->getPressure(Wstar, Id);

  switch (vf->getType(Id)) {
  case VarFcnBase::STIFFENEDGAS:
  case VarFcnBase::PERFECTGAS:
    eriemannfs_grad(R_1, U_1, P_1, 
		    R_i, U_i, P_i, vf, dWdW, Id, prec, mach, drdus);
    break;
    case VarFcnBase::TAIT: 
      fprintf(stderr, " ***ERROR: Tait EOS -> dW*/dn not implementted\n");
      exit(-1);
    //eriemannfs_tait_grad(R_1,U_1,P_1,R_i,U_i,P_i,vf,dWdW,Id); 
    //break;
  }

  dpdus = -dWdW[7];
 
  memset(dWstardn, 0, sizeof(double)*dim*3);
  
  for (int i=0; i<3; ++i) {
    dWstardn[i]     = dWdW[1]*Vi[(i+1)] + drdus*Vstar[i];  // drho* / dnn_Wall
    dWstardn[i+4*3] = dWdW[7]*Vi[(i+1)] + dpdus*Vstar[i];  // dP*   / dnn_Wall
  } 

  double unn = (Vstar[0] - Vi[1])*nphi[0] 
             + (Vstar[1] - Vi[2])*nphi[1] 
             + (Vstar[2] - Vi[3])*nphi[2];

  dWstardn[3]  = (Vstar[0] - Vi[1])*nphi[0] + unn; // du* / dnx_Wall
  dWstardn[7]  = (Vstar[1] - Vi[2])*nphi[1] + unn; // dv* / dny_Wall
  dWstardn[11] = (Vstar[2] - Vi[3])*nphi[2] + unn; // dw* / dnz_Wall
  
  dWstardn[4]  = (Vstar[1] - Vi[2])*nphi[0]; // du* / dny_Wall
  dWstardn[5]  = (Vstar[2] - Vi[3])*nphi[0]; // du* / dnz_Wall

  dWstardn[6]  = (Vstar[0] - Vi[1])*nphi[1]; // dv* / dnx_Wall
  dWstardn[8]  = (Vstar[2] - Vi[3])*nphi[1]; // dv* / dnz_Wall

  dWstardn[9]  = (Vstar[0] - Vi[1])*nphi[2]; // dw* / dnx_Wall
  dWstardn[10] = (Vstar[1] - Vi[2])*nphi[2]; // dw* / dny_Wall

}

//------------------------------------------------------------------------------


template<int dim>
inline
void LocalRiemannFluidStructure<dim>::computeRiemannSolution(int tag, double *Vi, double *Vstar,
                            double *nphi, VarFcn *vf,
                            double *Wstar, double *rupdatej,
                            double &weightj, int it)
{
  // Adam 2010.08.18
  // This function doesn't seem to be used anymore.
  // To be removed in a couple of months
  fprintf(stderr,"Oh Sorry ! Please uncomment the function (LocalRiemmannDesc.h:1596). I thought it wasn't needed anymore\n");
  exit(-1);
  /*
  // Commented by Adam on 2010.08.17 because it has to handle dim > 5
  //  int dim = 5;

  double P_1, U_1, R_1; // pass to 1D-FSI Riemann solver
  double P_i, U_i, R_i; // solution given by 1D-FSI Riemann solver

  double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
  double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};


  R_1 = Vi[0];
  U_1 = vni;
  P_1  = vf->getPressure(Vi);

  U_i = Vstar[0]*nphi[0]+Vstar[1]*nphi[1]+Vstar[2]*nphi[2];
  eriemannfs(R_1,U_1,P_1,R_i,U_i,P_i,vf,tag); //caution: U_i will not be modified!

  Wstar[0]  = R_i;                     Wstar[dim]    = Wstar[0];
  Wstar[1]  = vti[0]+U_i*nphi[0];      Wstar[dim+1]  = Wstar[1];
  Wstar[2]  = vti[1]+U_i*nphi[1];      Wstar[dim+2]  = Wstar[2];
  Wstar[3]  = vti[2]+U_i*nphi[2];      Wstar[dim+3]  = Wstar[3];
  Wstar[4]  = P_i;                     Wstar[dim+4]  = Wstar[4];
  if(dim > 5)
    {
      Wstar[5]  = 0.0;                 Wstar[dim+5]  = 0.0; // Boundary Condition: nuTilde = 0
    }

  if(it==1){
    weightj += 1.0;
    for (int k=0; k<dim; k++)
      rupdatej[k] += Wstar[k];
  }
  */
}

//------------------------------------------------------------------------------
//Caution: "ui" will not be modified!
template<int dim>
inline
int LocalRiemannFluidStructure<dim>::eriemannfs(double rho, double u, double p,
                                            double &rhoi, double ui, double &pi,
                                            VarFcn *vf, int Id,int& err, 
                                            double pc, double rhoc, bool prec, double beta) 
{
  // assume structure on the left of the fluid
  // using the notation of Toro's paper

  err = 0;

  double gamma = vf->getGamma(Id);
  double pref  = vf->getPressureConstant(Id);

  if(u==ui){ // contact
    rhoi = rho;
    pi   = p;
    return 0;
  }

  if (!prec) {
    if(ui<u){ // rarefaction
      double power = gamma/(gamma-1.0);
      double a = sqrt(gamma*(p+pref)/rho);
      double pbar = p + pref;
      double dee = 0.5*(gamma-1.0)*(ui-u)/a + 1.0;
      if(dee < 0) { pi = pc; rhoi = rhoc; }
      else {
        pi = pbar*pow(dee,2*power)-pref;
        rhoi = rho*pow((pi+pref)/(p+pref), 1.0/gamma);
      }
    }
    else{ // shock
      double temp = ((gamma+1.0)*rho*(ui-u)*(ui-u))/2.0;
      pi = p + 0.5*temp + sqrt(0.25*temp*temp + 2.0*gamma*temp*(p+pref)/(gamma+1.0));
      temp = (gamma-1.0)/(gamma+1.0);
      double pstarbar = pi + pref;
      double pbar = p + pref;
      rhoi = rho*(pstarbar/pbar+temp)/(temp*pstarbar/pbar+1);
    }
  } else {

    double a = sqrt(gamma*(p+pref)/rho);
    double X = 4.0*a*a*beta*beta + ui*ui*pow(beta*beta-1.0,2.0);
    double sqrtX = sqrt(X);
    double lambda = 0.5*(ui*(1.0+beta*beta) - sqrtX);
    double dp = rho*a*a*(u-ui)*beta*beta/(lambda-beta*beta*ui);

    pi = p + dp;
    double power = 2.0*gamma/(gamma-1.0);
    double pbar = p + pref;
    rhoi = rho*pow((pi+pref)/(p+pref), 1.0/gamma);
  }

    pi = std::max<double>(  pi,   pc);
  rhoi = std::max<double>(rhoi, rhoc);

  /*if(pi <= pc || rhoi <= rhoc){
    std::cout << "*** ERROR FS-ERS: detected too small density or pressure " << std::endl;
    std::cout << " ri, Pi " << std::endl;
    std::cout << rhoi << " " << pi << std::endl; 
    err = 1;
  }*/ 

  return err;
}

template<int dim>
inline
void LocalRiemannFluidStructure<dim>::eriemannfs_grad(double rho, double u, double p,
						      double &rhoi, double ui, double &pi,
						      VarFcn *vf, double* dWidWi,int Id,
						      bool prec, double beta, double &drdus) //Caution: "ui" will not be modified!
{

  // assume structure on the left of the fluid
  // using the notation of Toro's paper

  double gamma = vf->getGamma(Id);
  double pref  = vf->getPressureConstant(Id);
  memset(dWidWi, 0,sizeof(double)*9);

  if(fabs(u-ui)<1e-14){ // contact
    dWidWi[0] = 1.0; 
    dWidWi[8] = 1.0;
    return;
  }

// CHF: check that the next line is correct
  double a = sqrt(gamma*(p+pref)/rho);
  double q = (gamma+1.0)/(gamma-1.0);
  if(ui<u){ // rarefaction
    double power = 2*gamma/(gamma-1.0);
    double pbar = p + pref;

    double s = 0.5*(gamma-1.0)*(ui-u)/a + 1.0;
//    if (s < 0.0) {
//      fprintf(stderr,"Warning: s (%lf) in fs_grad is < 0!\n",s);
//    }
    // CHF: debug think this line is correct? ... commented line could result in negative base in power operation
    double eta = pbar*power*pow(s*s,q*0.5);
    //double eta = pbar*power*pow(s,power-1.0);
    double xi = eta*(-0.5/(a*a)*(gamma-1.0)*(ui-u));
 
    double dadp = 0.5/a*(gamma/rho), dadrho = -0.5*a/rho;
     
    if (!prec) {
      // dpi/dp
      dWidWi[8] = (pi+pref)/pbar+xi*dadp;
      // dpidrho
      dWidWi[6] = xi*dadrho;
      // dpidu
      dWidWi[7] = eta*(-0.5/a*(gamma-1.0));
      
      double mu = rho/gamma*pow((pi+pref)/pbar, (1.0-gamma)/gamma); 
      dWidWi[2] = mu*(1.0/pbar*dWidWi[8]-(pi+pref)/(pbar*pbar));
      dWidWi[1] = mu*dWidWi[7]/pbar;
      dWidWi[0] = rhoi/rho+mu/pbar*dWidWi[6];

      drdus = -(mu/pbar)*dWidWi[7]; //**

    }
    else {

      double X = 4.0*a*a*beta*beta + ui*ui*pow(beta*beta-1.0,2.0);
      double sqrtX = sqrt(X);
      double lambda = 0.5*(ui*(1.0+beta*beta) - sqrtX);
      double q = lambda-beta*beta*ui;
  
      // dpi/dp
      dWidWi[8] = 1.0+beta*beta*(u-ui)/(q*q)*(gamma*q+gamma*a*a*beta*beta/sqrtX); // PJSA
      // dpidrho
      dWidWi[6] = beta*beta*(u-ui)*(a*a/q + (p+pref)*(-gamma/rho*q-a*a*beta*beta*gamma/(rho*sqrtX))/(q*q)); // PJSA
      // dpidu
      dWidWi[7] = beta*beta*rho*a*a/q;
      
      rhoi = rho*pow((pi+pref)/(p+pref), 1.0/gamma);
      double rhopp = rho*pow((pi+pref)/(p+pref), 1.0/gamma-1.0);
      double dpp = 1.0/gamma*rhopp/(p+pref);
      dWidWi[2] = dpp*dWidWi[8] - rhoi/gamma/(p+pref);
      dWidWi[1] = dpp*dWidWi[7];
      dWidWi[0] = dpp*dWidWi[6] + pow((pi+pref)/(p+pref), 1.0/gamma);

      //drdus = ?

    }
  }
  else{ // shock
    double power = 2*gamma/(gamma+1.0);
    double t = ((gamma+1)*rho*(ui-u)*(ui-u))/2.0;
    double pstarbar = pi + pref;
    double pbar = p + pref;
 
    double dtdrho = t/rho, dtdu = -(gamma+1.0)*rho*(ui-u);
    double xi = sqrt(0.25*t*t+power*t*pbar);
  
    double eta = 0.5+0.5/xi*(0.5*t+power*pbar);
    if (!prec) {
      dWidWi[8] = 1.0+0.5/xi*power*t;
      dWidWi[7] = eta*dtdu;
      dWidWi[6] = eta*dtdrho;
      
//    dWidWi[8] = 1.0+power*sqrt(t/(t/2.0+power*pbar))/2.0;
//    
//    double tmp = gamma*pbar*rho;
//    dWidWi[7] = dtdu/2.0-(dtdu*dtdu+8.0*tmp)/sqrt(8.0*dtdu*dtdu+64.0*tmp);
//
//    dWidWi[6] = dtdrho/2.0-(dtdu*dtdu*dtdu+8.0*tmp*dtdu)/sqrt(32.0*dtdu*dtdu+256.0*tmp)/(gamma+1.0)/(rho*rho);
      
      q = 1/q; // PJSA
      double s = q*pstarbar/pbar+1.0;
      double deriv = 1.0/(pbar*s)-(pstarbar/pbar+q)/(s*s)*(q/pbar);
      double deriv2 = -pstarbar/(pbar*pbar*s)+(pstarbar/pbar+q)*(q*pstarbar/(pbar*pbar))/(s*s);
      dWidWi[0] = rhoi/rho+rho*deriv*dWidWi[6];
      dWidWi[1] = rho*deriv*dWidWi[7];
      dWidWi[2] = rho*(deriv*dWidWi[8]+deriv2);

      drdus = -(rho*deriv)*dWidWi[7]; //**

    } else {

      double X = 4.0*a*a*beta*beta + ui*ui*pow(beta*beta-1.0,2.0);
      double sqrtX = sqrt(X);
      double lambda = 0.5*(ui*(1.0+beta*beta) - sqrtX);
      double q = lambda-beta*beta*ui;
  
      // dpi/dp
      dWidWi[8] = 1.0+beta*beta*(u-ui)/(q*q)*(gamma*q+gamma*a*a*beta*beta/sqrtX); // PJSA
      // dpidrho
      dWidWi[6] = beta*beta*(u-ui)*(a*a/q + (p+pref)*(-gamma/rho*q-a*a*beta*beta*gamma/(rho*sqrtX))/(q*q)); // PJSA
      // dpidu
      dWidWi[7] = beta*beta*rho*a*a/q;
      
      rhoi = rho*pow((pi+pref)/(p+pref), 1.0/gamma);
      double rhopp = rho*pow((pi+pref)/(p+pref), 1.0/gamma-1.0);
      double dpp = 1.0/gamma*rhopp/(p+pref);
      dWidWi[2] = dpp*dWidWi[8] - rhoi/gamma/(p+pref);
      dWidWi[1] = dpp*dWidWi[7];
      dWidWi[0] = dpp*dWidWi[6] + pow((pi+pref)/(p+pref), 1.0/gamma);

      //drdus = ?

    }
  }
}

template<int dim>
inline
void LocalRiemannFluidStructure<dim>::eriemannfs_tait(double rho, double u, double p,
						      double &rhoi, double ui, double &pi,
						      VarFcn *vf, int Id,int& err,
                                                      double pc, double rhoc) //Caution: "ui" will not be modified!
{
  // assume structure on the left of the fluid
  // using the notation of Toro's paper
 
  err = 0;

  int max_ite = 100;

  double a = vf->getAlphaWater(Id);
  double b = vf->getBetaWater(Id);
  double pref  = vf->getPrefWater(Id);
  double Udummy[5],Vdummy[5]={0,0,0,0,0};

  if(u==ui){ // contact
    rhoi = rho;
    pi   = p;
    return;
  }

  if(ui<u){ // rarefaction
    
    double ud = u-ui;
    double q = 2.0*sqrt(a*b)/(b-1.0);
    double qp = -ud/q+pow(rho,(b-1.0)*0.5);
    if (qp > 0.0)
      rhoi = pow(qp, 2.0/(b-1.0));
    else {
      // will be handled by verification below.
      rhoi = 0.0;
    }
     
    //pi = a*pow(rhoi,b)+pref;
    Vdummy[0] = rhoi;
    vf->getVarFcnBase(Id)->verification(0,Udummy,Vdummy);
    rhoi = Vdummy[0];
    pi = a*pow(rhoi,b)+pref;
    pi = std::max<double>(pi,pc);
    if (pi == pc)
    rhoi = pow((pi-pref)/a,1.0/b);
    rhoi = std::max<double>(rhoi,rhoc);
  }
  else{ // shock
    
    // Must solve a nonlinear equation.  I can't find an
    // analytical solution.
    rhoi = 1.001*rho;
    pi = a*pow(rhoi,b)+pref;
    double V,dV,dpdrho;
    int i = 0;
    do {
      //fprintf(stderr,"%d %lf %lf %lf %lf %lf %lf %lf\n",i,rho,u,p,rhoi,ui,pi,fabs(u-ui-V));
      dpdrho = a*b*pow(rhoi, b-1.0);
      V = sqrt(fabs( (p-pi)*(1.0/rhoi-1.0/rho) ) );
      dV = 0.5/V*(-dpdrho*(1.0/rhoi-1.0/rho) - (p-pi)*(1.0/(rhoi*rhoi) ) );
      assert(rhoi >= rho);

      if (fabs(V+u-ui) < 1.0e-6 || fabs(V+u-ui)/max(1e-8,fabs(ui)) < 1e-6)
	break;

      rhoi -= (V+u-ui) / dV;
      pi = a*pow(rhoi,b)+pref;
      pi = std::max<double>(pi,pc);
      if (pi == pc)
        rhoi = pow((pi-pref)/a,1.0/b);
      rhoi = std::max<double>(rhoi,rhoc);    

      ++i;
    } while (i < max_ite);

    if (i >= max_ite) {
      err = 1;
      // fprintf(stderr,"%d %lf %lf %lf %lf %lf %lf %lf\n",i,rho,u,p,rhoi,ui,pi,fabs(u-ui+V));
      std::cout << "*** Warning FS-ERS Tait: Newton reached max num. iterations " << max_ite  << std::endl;
      std::cout << "without converging to the desired tolerance " << 1.0e-6 << std::endl;
    }
  }

  /*if(pi <= pc || rhoi <= rhoc){
    std::cout << "*** ERROR FS-ERS Tait: detected too small density or pressure " << std::endl;
    std::cout << " ri, Pi " << std::endl;
    std::cout << rhoi << " " << pi << std::endl; 
    err = 1;
  }*/

}

template<int dim>
inline
void LocalRiemannFluidStructure<dim>::eriemannfs_tait_grad(double rho, double u, double p,
							   double &rhoi, double ui, double &pi,
							   VarFcn *vf, double* dWidWi,int Id) //Caution: "ui" will not be modified!
{

  // assume structure on the left of the fluid
  // using the notation of Toro's paper

  double a = vf->getAlphaWater(Id);
  double b = vf->getBetaWater(Id);
  double pref  = vf->getPrefWater(Id);
  memset(dWidWi, 0,sizeof(double)*9);
  if(u==ui){ // contact
    dWidWi[0] = 1.0; 
    dWidWi[8] = 1.0;
    return;
  }

  double dpdrho = a*b*pow(rho, b-1.0);
  double dpdrhos = a*b*pow(rhoi, b-1.0);
  double dVdrho,dVdrhos; 
  if(ui<u){ // rarefaction
     
    dVdrho = -sqrt(a*b*pow(rho,b-3.0));//sqrt(a*b)*(b-3.0)*(b-5.0)/4.0*pow(rho, (b-7.0)*0.5);
    dVdrhos = sqrt(a*b*pow(rhoi,b-3.0));//-sqrt(a*b)*(b-3.0)*(b-5.0)/4.0*pow(rhoi, (b-7.0)*0.5);    
  }
  else{ // shock
    
    double V;
    V = sqrt(fabs( (p-pi)*(1.0/rhoi-1.0/rho) ) );
    dVdrhos = 0.5/V*(-dpdrhos*(1.0/rhoi-1.0/rho) - (p-pi)*(1.0/(rhoi*rhoi) ) );
    dVdrho = 0.5/V*((dpdrho)*(1.0/rhoi-1.0/rho) + (p-pi)*(1.0/(rho*rho) ) );
  }

  dWidWi[0] = -dVdrho / dVdrhos;
  dWidWi[1] = -1.0 / dVdrhos;//*(ui<u?1.0:-1.0);
  dWidWi[2] = 0.0;

  dWidWi[8] = 1.0;
}







//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Daniel 2017.02.7
// This Local Riemann for actuator disk simulation
// The disk is immersed in single fluid.
// We solve the Riemann problem to force Dp pressure jump, but continuous velocity
// and continuous density
template<int dim>
class LocalRiemannActuatorDisk : public LocalRiemann {

public:
    LocalRiemannActuatorDisk() : LocalRiemann() { fluid1 = fluid2 = 0; }

    LocalRiemannActuatorDisk(VarFcn *vf) : LocalRiemann(vf, 0, 0){ fluid1 = fluid2 = 0; }

    virtual ~LocalRiemannActuatorDisk() { vf_ = 0; }


    int computeRiemannSolution(double *Vi, double *Vj, double *Vstar, double dp,
                               double *n_s, double *n_f, VarFcn *vf,
                               int it, double *Wi, double *Wj,int Id = 0);


    void computeSourceTerm(double *Vi, double *Vj,double dp,
                                                              double *n_s, double *n_f, VarFcn *vf,
                                                              double *flux, bool method = true, int Id = 0);

private:

    void riemannActuatorDisk(double rho_l, double v_l, double p_l,
                             double rho_r, double v_r, double p_r,
                             double vstar_n,
                             double dp,   VarFcn *vf, int Id,
                             double &rho_a, double &v_a, double &p_a,
                             int& err);

    void pressureFunction(double p, double rho_k, double p_k, double a_k, double A_k, double B_k,double gamma,
                          double &f_k, double &df_k);


    bool checkSolution(int diskCase, double rho_l, double v_l, double p_l,
                       double rho_r, double v_r, double p_r, double vstar_n,
                       double v_m, double p_m,
                       double dp_l,double dp_r,   double gamma,double&rho_mr, double &rho_ml,bool DEBUG=false);

    int solveContactDiscontinuity(double rho_l, double v_l, double p_l,
                                  double rho_r, double v_r, double p_r,
                                  double dp_l, double dp_r, double gamma,
                                  double &v_m, double &p_m);
};

template<int dim>
inline
int LocalRiemannActuatorDisk<dim>::computeRiemannSolution(double *Vi, double *Vj,double *Vstar,double dp,
                                                          double *n_s, double *n_f, VarFcn *vf,
                                                          int it, double *W_Ri, double *W_Rj,int Id) {
    /* Vi Vj are two fluid primitive state variables
     * dp pressure jump
     * n_s is unit structure normal, n_f is fluid normal(can be non unit)
     * vf state of equation function
     * W_Ri,WRj riemann primitive state variables on both side of actuator disk, near i and j
     * Id fluid id
     */
    //---------------------------------------------------------------

    double v_ni = Vi[1] * n_s[0] + Vi[2] * n_s[1] + Vi[3] * n_s[2]; //normal velocity of node i
    double v_ti[3] = {Vi[1] - v_ni * n_s[0], Vi[2] - v_ni * n_s[1], Vi[3] - v_ni * n_s[2]}; // tangential velocity of node i
    double v_nj = Vj[1] * n_s[0] + Vj[2] * n_s[1] + Vj[3] * n_s[2]; //normal velocity of node j
    double v_tj[3] = {Vj[1] - v_nj * n_s[0], Vi[2] - v_ni * n_s[1], Vi[3] - v_ni * n_s[2]}; // tangential velocity of node j
    double vstar_n = Vstar[0]*n_s[0] + Vstar[1]*n_s[1] + Vstar[2]*n_s[2];

    double rc = vf->getVarFcnBase(Id)->rhomin;
    double pc = vf->getVarFcnBase(Id)->pmin;

    bool iIsUpstream = (n_s[0] * n_f[0] + n_s[1] * n_f[1] + n_s[2] * n_f[2] >= 0);
    double rho_l,v_l,p_l,rho_r,v_r,p_r,rho_a,v_a,p_a;//project to 1D ,upstream fluid states *_l, downstream fluid states *_r, actuator disk upstream states *_a;
    if (iIsUpstream) { //  l  -> +dp r
        rho_l = std::max(rc, Vi[0]);
        v_l = v_ni;
        p_l = std::max(pc, vf->getPressure(Vi, Id));
        rho_r = std::max(rc, Vj[0]);
        v_r = v_nj;
        p_r = std::max(pc, vf->getPressure(Vj, Id));
    } else {
        rho_l = std::max(rc, Vj[0]);
        v_l = v_nj;
        p_l = std::max(pc, vf->getPressure(Vj, Id));
        rho_r = std::max(rc, Vi[0]);
        v_r = v_ni;
        p_r = std::max(pc, vf->getPressure(Vi, Id));
    }


    int err = 0;

    switch (vf->getType(Id)) {
        case VarFcnBase::STIFFENEDGAS:

            fprintf(stderr, "ERROR: NO IMPLEMENTATION FOR ACTUATOR DISK FOR STIFFENEDGAS use the one for PERFECTGAS!\n");
        case VarFcnBase::PERFECTGAS:
            riemannActuatorDisk(rho_l, v_l, p_l, rho_r, v_r, p_r, vstar_n, dp, vf, Id, rho_a, v_a, p_a,err);
            break;
        case VarFcnBase::TAIT:
            fprintf(stderr, "ERROR: NO IMPLEMENTATION FOR ACTUATOR DISK FOR TAIL GAS!\n");
            break;
    }

    if (err)
        return err;
    //The velocity is no-slip condition boundary
    W_Ri[0] = W_Rj[0] = rho_a;
    W_Ri[1] = v_a * n_s[0] + v_ti[0] ;
    W_Ri[2] = v_a * n_s[1] + v_ti[1] ;
    W_Ri[3] = v_a * n_s[2] + v_ti[2] ;
    W_Rj[1] = v_a * n_s[0] + v_tj[0] ;
    W_Rj[2] = v_a * n_s[1] + v_tj[1] ;
    W_Rj[3] = v_a * n_s[2] + v_tj[2] ;
    if (iIsUpstream) { //  l  -> +dp r
        W_Ri[4] = p_a;
        W_Rj[4] = p_a + dp;

    } else {
        W_Rj[4] = p_a;
        W_Ri[4] = p_a + dp;
    }
    if (dim == 6) {
        W_Ri[5] = W_Rj[5] = 0.0;// Boundary Condition: nuTilde = 0
    } else if (dim == 7)  // Boundary Condition for KE. To be improved with Wall Function...
    {
        W_Ri[5] = W_Rj[5] = 0.0;
        W_Ri[6] = W_Rj[6] = 0.0;
    }


    for (int i = 0; i < dim; i++){
        W_Ri[dim + i] = W_Ri[i];
        W_Rj[dim + i] = W_Rj[i];
    }

    return err;
}

template<int dim>
inline
void LocalRiemannActuatorDisk<dim>::computeSourceTerm(double *Vi, double *Vj,double dp,
                                                          double *n_s, double *n_f, VarFcn *vf,
                                                          double *flux, bool method, int Id) {
    /* Vi Vj are two fluid primitive state variables
     * dp pressure jump value
     * n_s is unit structure normal, n_f is fluid edge area normal(non unit)
     * vf state of equation function
     * flux source term
     * method, true for corrected one and false for traditional one
     * Id fluid id
     */


  double gamma = vf->getGamma(Id);
  double Vel[3] = {(Vi[1]+Vj[1])/2.0,(Vi[2]+Vj[2])/2.0,(Vi[3]+Vj[3])/2.0};//use average velocity
  double faceArea = abs(n_s[0]*n_f[0] +n_s[1]*n_f[1] +n_s[2]*n_f[2]);
  double normal[3] = {faceArea*n_s[0], faceArea*n_s[1],faceArea*n_s[2]};//use structure normal
  flux[0] = 0;
  flux[1] = dp*normal[0];
  flux[2] = dp*normal[1];
  flux[3] = dp*normal[2];

  double normalVelocity = Vel[0]*normal[0] + Vel[1]*normal[1] + Vel[2]*normal[2];
  flux[4] = method? gamma/(gamma-1)*dp*normalVelocity: dp*normalVelocity;
    /*
  fprintf(stderr, " ***ERROR: Actuator disk SourceTerm\n");
  fprintf(stderr, " ***ERROR: Actuator disk Vi %.10f,%.10f,%.10f,%.10f,%.10f, Vj %.10f,%.10f,%.10f,%.10f,%.10f\n",
                                               Vi[0],Vi[1],Vi[2],Vi[3],Vi[4],      Vj[0],Vj[1],Vj[2],Vj[3],Vj[4]);
  fprintf(stderr, " ***ERROR: n_s %.10f,%.10f,%.10f, n_f %.10f,%.10f,%.10f \n", n_s[0],n_s[1],n_s[2],n_f[0],n_f[1],n_f[2]);
  fprintf(stderr, " ***ERROR: normal  %.10f,%.10f,%.10f\n", normal[0],normal[1],normal[2]);
  fprintf(stderr, " ***ERROR: dp  %.10f\n", dp);
  fprintf(stderr, " ***ERROR: normalVelocity  %.10f\n", normalVelocity);
  fprintf(stderr, " ***ERROR: Actuator disk flux %.10f,%.10f,%.10f,%.10f,%.10f,\n\n",flux[0],flux[1],flux[2],flux[3],flux[4]);
*/

  }
//------------------------------------------------------------------------------
template<int dim>
inline
void LocalRiemannActuatorDisk<dim>::riemannActuatorDisk(double rho_l, double v_l, double p_l,
                                                    double rho_r, double v_r, double p_r, double vstar_n,
                                                    double dp,   VarFcn *vf, int Id,
                                                    double &rho_a, double &v_a, double &p_a,
                                                    int& err)
{
// assume rho_l, v_l and p_l are fluid states at upstream, and rho_r, v_r, p_r are fluid states at downstream,
// vstar_n is the actuator disk velocity
// using the notation of Toro's paper
// return the fluid states rho_a, v_a, p_a upstream the actuator disk
double gamma = vf->getGamma(Id);
double pref  = vf->getPressureConstant(Id);
double v_m, p_m, rho_ml, rho_mr;
double M_l = v_l/sqrt(gamma*p_l/rho_l);
const int LEFT = 0, CENTER_LEFT = 1, CENTER_RIGHT = 2, RIGHT = 3;
               // first wave, actuator disk, contact wave,  third wave
        // rho_l ,      rho_a,          rho_a=rho_ml   rho_mr         rho_r
        // v_l,         v_a             v_a  =v_ml     v_mr           v_r
        // p_l          p_a             p_a+dp =p_ml   p_mr           p_r
        err += solveContactDiscontinuity(rho_l, v_l, p_l, rho_r, v_r, p_r, -dp, 0, gamma,v_m,p_m);
        if(checkSolution(CENTER_LEFT, rho_l, v_l, p_l, rho_r, v_r, p_r, vstar_n,v_m, p_m,-dp,0,gamma, rho_mr, rho_ml)) {
            rho_a = rho_ml;
            v_a = v_m;
            p_a = p_m - dp;
            return;
        }

              //actuator disk, first wave, contact wave, third wave
        // rho_a=rho_l,   rho_a,      rho_ml   rho_mr         rho_r
        // v_a=v_l        v_a           v_ml     v_mr           v_r
        // p_a=p_l        p_a+dp        p_ml     p_mr           p_r
        err += solveContactDiscontinuity(rho_l, v_l, p_l+dp, rho_r,v_r,p_r, 0,0,gamma, v_m,p_m);
        if(checkSolution(LEFT, rho_l, v_l, p_l + dp, rho_r, v_r, p_r, vstar_n,v_m, p_m,0 ,0,gamma, rho_mr, rho_ml)) {
            rho_a = rho_l;
            v_a = v_l;
            p_a = p_l;
            return;
        }


             // first wave, contact wave, actuator disk, third wave
        // rho_l ,      rho_ml,      rho_mr=rho_a    rho_a         rho_r
        // v_l,         v_ml           v_mr=v_a      v_a=v_ml       v_r
        // p_l          p_ml           p_mr=p_a      p_a + dpl      p_r
        err += solveContactDiscontinuity(rho_l, v_l, p_l, rho_r, v_r, p_r, 0, dp, gamma, v_m, p_m);
        if(checkSolution(CENTER_RIGHT, rho_l, v_l, p_l, rho_r, v_r, p_r,vstar_n,v_m, p_m,0 ,dp,gamma, rho_mr, rho_ml)) {
            rho_a = rho_mr;
            v_a = v_m;
            p_a = p_m;
            return;
        }

            //  first wave, contact wave, third wave ,actuator disk
        // rho_l ,      rho_ml      rho_mr         rho_a=rho_r   rho_r
        // v_l,         v_ml        v_mr           v_a=v_r        v_r
        // p_l          p_ml        p_mr           p_a=p_r-dp     p_r
        err += solveContactDiscontinuity(rho_l, v_l, p_l, rho_r, v_r, p_r - dp, 0, 0, gamma, v_m, p_m);
        if(checkSolution(RIGHT, rho_l, v_l, p_l, rho_r, v_r, p_r-dp,vstar_n, v_m, p_m,0 ,dp,gamma, rho_mr, rho_ml)) {
            rho_a = rho_r;
            v_a = v_r;
            p_a = p_r - dp;
            return;
        }

    fprintf(stderr, " ***ERROR: Actuator disk Riemann solver has no solution, use approximate solution");
    if( M_l >= 0.0 ){
        rho_a = rho_l;
        v_a = v_l;
        p_a = p_l;
    }else{
        rho_a = rho_r;
        v_a = v_r;
        p_a = p_r - dp;
    }
}



template<int dim>
inline
void LocalRiemannActuatorDisk<dim>::pressureFunction(double p, double rho_k, double p_k, double a_k, double A_k, double B_k,double gamma,
                                                     double &f_k, double &df_k)  {

    if(p > p_k) { //left shock

        f_k = (p - p_k) * sqrt(A_k / (p + B_k));

        df_k = sqrt(A_k / (B_k + p)) * (1 - (p - p_k) / (2 * (B_k + p)));
    } else {  // rarefaction


        f_k = 2 * a_k / (gamma - 1) * ( pow(p / p_k, (gamma - 1) / (2 * gamma)) - 1);

        df_k = 1 / (rho_k * a_k) * pow(p / p_k, -(gamma + 1) / (2 * gamma));
    }

}

template<int dim>
inline
int LocalRiemannActuatorDisk<dim>::solveContactDiscontinuity(double rho_l, double v_l, double p_l,
                                                             double rho_r, double v_r, double p_r,
                                                             double dp_l, double dp_r, double gamma,
                                                             double &v_m, double &p_m) {
    // The pressure and density at contact discontinuity is p_m, v_m
    // depends on the position of the actuator disk we have
    // The pressure after  the first wave is p_m + dp_l
    // The pressure before the third wave is p_m + dp_r
    // p is the pressure at contact discontinuity

    int MAX_ITE = 100;
    double TOLERANCE = 1.0e-12;
    bool found = false;
    double f_l,df_l,f_r,df_r;


    double d_v = v_r - v_l;
    double a_l = sqrt(gamma * p_l / rho_l) , a_r = sqrt(gamma * p_r / rho_r);
    double A_l = 2 / ((gamma + 1) * rho_l), A_r = 2 / ((gamma + 1) * rho_r);
    double B_l = (gamma - 1) / (gamma + 1) * p_l, B_r = (gamma - 1) / (gamma + 1) * p_r;

    double p_old = (p_l + p_r)/2.0;

    for(int i = 0; i < MAX_ITE; i++) {

        pressureFunction(p_old + dp_l, rho_l, p_l, a_l, A_l, B_l, gamma, f_l, df_l);
        pressureFunction(p_old + dp_r, rho_r, p_r, a_r, A_r, B_r, gamma, f_r, df_r);

        p_m = p_old - (f_l + f_r + d_v) / (df_l + df_r);

        if (p_m < 0.0)        p_m = TOLERANCE;

        if (2 * fabs(p_m - p_old) / (p_m + p_old) < TOLERANCE) {
            found = true;
            break;
        }
        p_old = p_m;
    }
    if (!found)    {
        fprintf(stderr, " ***ERROR: Divergence in Newton-Raphason iteration in Actuator disk Riemann solver");
        return 1;
    }

    v_m = 0.5 * (v_l + v_r + f_r - f_l);
    return 0;

}

template<int dim>
inline
bool LocalRiemannActuatorDisk<dim>::checkSolution(int diskCase, double rho_l, double v_l, double p_l,
                                                  double rho_r, double v_r, double p_r, double vstar_n,
                                                  double v_m, double p_m,
                                                  double dp_l,double dp_r,   double gamma,double&rho_mr, double &rho_ml,bool DEBUG) {
    //  diskCase 0,1,2,3, i means disk is between i wave and i+1 wave;
    //  left fluid state variable rho_l, v_l, p_l;
    //right fluid state variable rho_r, v_r, p_r;
    //vstar_n: actuator disk velociy
    //velocity and pressure after the 1- wave are v_m, p_m + dp_l
    //velocity and pressure before the 3- wave are v_m, p_m + dp_r
    //density before and after contact discontinuity are rho_ml, rho_mr
    //return true of false if it is false means the pattern cannot match with the diskCase
    if(p_r <= 0.0 || p_l <= 0.0 || rho_l <= 0.0 || rho_r <= 0.0 || p_r + dp_r <= 0.0 || p_l+dp_l <= 0.0) {
        fprintf(stderr, "****ERROR, In LocalRiemannActuatorDisk, has negative pressure or negative density\n");
        return false;
    }
    bool result = true;
    double a_l = sqrt(gamma * p_l / rho_l), a_r = sqrt(gamma * p_r / rho_r); //sound speed

    //Left side
    double p_ml = p_m + dp_l;
    if (p_ml < p_l) { // left rarefaction wave
        double s_l = v_l - a_l;
        double a_ml = a_l * pow(p_ml / p_l, (gamma - 1) / (2 * gamma));
        double s_ml = v_m - a_ml;
        rho_ml = rho_l * pow(p_ml/ p_l, 1 / gamma);
        if(DEBUG) fprintf(stderr,"DEBUG: case %d, left rarefaction, velocity s_l s_ml are %f and %f\n",diskCase, s_l, s_ml);
        if((diskCase == 0 && s_l < vstar_n && s_ml < vstar_n)||(diskCase == 1 && s_l > vstar_n && s_ml > vstar_n)
           ||(diskCase == 2 && s_l > vstar_n && s_ml > vstar_n) ||(diskCase == 3 && s_l > vstar_n && s_ml > vstar_n))
            result = false;

    }
    else { // left shock wave
        double s_shock = v_l - a_l * sqrt((gamma + 1) *p_ml / (2 * gamma * p_l) + (gamma - 1) / (2 * gamma));
        rho_ml = rho_l * (p_ml/ p_l + (gamma - 1) / (gamma + 1)) / ((gamma - 1) * p_ml/ ((gamma + 1) * p_l) + 1);
        if(DEBUG) fprintf(stderr,"DEBUG: case %d, left shock, velocity s_shock_l is %f\n",diskCase, s_shock);
        if((diskCase == 0 && s_shock < vstar_n )||(diskCase == 1 && s_shock > vstar_n )
           ||(diskCase == 2 && s_shock > vstar_n) ||(diskCase == 3 && s_shock > vstar_n))
            result = false;
    }
    if(DEBUG) fprintf(stderr,"DEBUG: case %d, contact discontinuity velocity is %f\n",diskCase,  v_m);
    if((diskCase == 0 && v_m < vstar_n )||(diskCase == 1 && v_m < vstar_n )
       ||(diskCase == 2 && v_m > vstar_n) ||(diskCase == 3 && v_m > vstar_n))
        result = false;

    //Right side
    double p_mr = p_m + dp_r;
    if (p_mr < p_r) { // right rarefaction wave

        double s_r = v_r + a_r;
        double a_mr = a_r * pow(p_mr / p_r, (gamma - 1) / (2 * gamma));
        double s_mr = v_m + a_mr;
        rho_mr = rho_r * pow(p_mr / p_r, 1 / gamma);
        if(DEBUG) fprintf(stderr,"DEBUG: case %d,  right rarefaction, velocity s_mr s_r are %f and %f\n",diskCase, s_mr, s_r);
        if((diskCase == 0 && s_r < vstar_n && s_mr < vstar_n)||(diskCase == 1 && s_r < vstar_n && s_mr < vstar_n)
           ||(diskCase == 2 && s_r < vstar_n && s_mr < vstar_n) ||(diskCase == 3 && s_r > vstar_n && s_mr > vstar_n))
            result = false;
    }

    else {   //right shock wave
        double s_shock = v_r + a_r * sqrt((gamma + 1) * p_mr / (2 * gamma * p_r) + (gamma - 1) / (2 * gamma));
        rho_mr = rho_r * (p_mr / p_r + (gamma - 1) / (gamma + 1)) / ((gamma - 1) * p_mr / ((gamma + 1) * p_r) + 1);
        if(DEBUG) fprintf(stderr,"DEBUG: case %d, right shock, velocity s_shock_r is %f\n" ,diskCase, s_shock);
        if((diskCase == 0 && s_shock < vstar_n )||(diskCase == 1 && s_shock < vstar_n )
           ||(diskCase == 2 && s_shock < vstar_n) ||(diskCase == 3 && s_shock > vstar_n))
            result = false;
    }

    return result;
}
//------------------------------------------------------------------------------
#endif
