#ifndef _LOCAL_RIEMANN_H
#define _LOCAL_RIEMANN_H

#include <IoData.h>
#include "VarFcn.h"
#include "ODEIntegrator.h"
#include "SparseGrid.h"
#include "SparseGridCluster.h"

//----------------------------------------------------------------------------
// Virtual base class to provide nodal values to compute fluxes at the interface
// between two fluids
class LocalRiemann {

protected:
  VarFcn *vf_;
public:
  int fluid1, fluid2;  //             fluid1 ~ phi>0           fluid2 ~ phi<0
                       //          ~ "outside" ~ "right"    ~ "inside" ~ "left"
                       // GasGas            Gas1                    Gas2
                       // GasTait           Tait                    Gas
                       // TaitTait          Tait1                   Tait2	 
                       // JWLJWL            JWL1                    JWL2
                       // GasJWL            Gas                     JWL
                       // FluidStruct       Gas                     N/A

  LocalRiemann()           { vf_ = 0; fluid1 = 0; fluid2 = 1;}
  LocalRiemann(VarFcn *vf, int tag1, int tag2) { vf_ = vf; fluid1 = tag1; fluid2 = tag2;}
  virtual ~LocalRiemann()  { vf_ = 0; }

  // multiphase Riemann problem
  virtual int updatePhaseChange(double *V, int ID, int IDn, double *newV, double weight,bool isCellCut)
    {fprintf(stderr,"updatePhaseChange is not implemented here!\n");return 0;}
  virtual int computeRiemannSolution(double *Vi, double *Vj,
				      int IDi, int IDj, double *nphi,
				      double *initWi, double *initWj,
				      double *Wi, double *Wj,
				      double *rupdatei, double *rupdatej,
				      double &weighti, double &weightj, 
				      double dx[3], int it, bool isHigherOrder) { return 0; } 

  virtual void computeRiemannJacobian(double *Vi, double *Vj,
				      int IDi, int IDj, double *nphi,
				      double *Wi, double *Wj,
				      double dx[3],int it,
				      double* dWidWi,double*  dWidWj,
				      double* dWjdWi,double*  dWjdWj) {} 

  // FS Riemann problem
  virtual int computeRiemannSolution(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, int Id = 0) { return 0; } 

  virtual void computeRiemannJacobian(double *Vi, double *Vstar,
                              double *nphi, VarFcn *vf,
                              double *Wstar, double *rupdatei,
                              double &weighti, int it, double* dWstardU,int Id = 0) {}

  virtual void computeRiemannderivative(double *Vi, double *Vstar, 
					double *nphi, VarFcn *vf, double *Wstar, 
					double* dWstardn, int Id = 0) {}

  // Actuator disk problem
  virtual int computeRiemannSolution(double *Vi, double *Vj,double *Vstar,double dp,
				     double *n_s, double *n_f, VarFcn *vf,
				     int it, double *Wi, double *Wj,int Id = 0){return 0;}

  virtual  void computeSourceTerm(double *Vi, double *Vj,double dp,
                          double *n_s, double *n_f, VarFcn *vf,
                          double *flux, bool method = true, int Id = 0){};
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// class used when the original GFMP is used (no need to solve a two-phase
//                                            Riemann problem)
class LocalRiemannGfmp : public LocalRiemann {

public:
  LocalRiemannGfmp() : LocalRiemann() {}
  LocalRiemannGfmp(VarFcn *vf, int tag1, int tag2) : LocalRiemann(vf,tag1,tag2) {}
  virtual ~LocalRiemannGfmp() { vf_ = 0; }

  int updatePhaseChange(double *V, int ID, int IDn, double *newV, double weight,bool isCellCut){///*nothing to do for GFMP*/}
    //if(ID != IDn) fprintf(stdout, "node changes from phase %d to phase %d!\n", IDn, ID);
    return 0;
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// class used when the exact two-phase Riemann problem is solved
class LocalRiemannGfmpar : public LocalRiemann {

  MultiFluidData::TypePhaseChange phaseChangeType_;

public:
  LocalRiemannGfmpar() : LocalRiemann() {}
    LocalRiemannGfmpar(VarFcn *vf, int tag1, int tag2, MultiFluidData::TypePhaseChange phaseChangeType,
		       double relaxationFacJwl = 1.0,
		       double refd = 0.0,double refe = 0.0) : LocalRiemann(vf,tag1,tag2), phaseChangeType_(phaseChangeType), refdensity(refd), refentropy(refe),
      relaxFactorJwl(relaxationFacJwl) {}
  virtual ~LocalRiemannGfmpar() { vf_ = 0; }

  int updatePhaseChange(double *V, int ID, int IDn, double *newV, double weight,bool isCellCut){
    if(ID == IDn && !isCellCut) return 0; /* node does not change phase: nothing to do*/
    if(weight<=0.0)  { if (IDn >= 0 && ID >= 0) { fprintf(stdout, "*** Error: negative weight in LocalRiemannGfmpar::updatePhaseChange %d %d\n",ID, IDn);
      DebugTools::PrintBacktrace();
                       return -1;} }
    else
      for(int k=0; k<5; k++) V[k] = newV[k]/weight;

    return 0;
  }

protected:

  double refdensity,refentropy,relaxFactorJwl;

  // following functions used when Riemann problem formulated as a system of nonlinear equations.
  bool solve2x2System(double *mat, double *rhs, double *res);
  // valid for JWL phase
  int rarefactionJWL(double phi,
                   double v1, double u1, double p1, 
                   double v,  double &u, double &p, 
                   double &du, double &dp,
                   MultiFluidData::RiemannComputation type = MultiFluidData::RK2, int flag = 0);
  virtual int riemannInvariantGeneralTabulation(double *in, double *res);
  void riemannInvariantGeneral1stOrder(double *in, double *res, double *phi);
  void riemannInvariantGeneral2ndOrder(double *in, double *res, double *phi);

  int rarefactionJWLderivs(double phi,
		           double v1, double u1, double p1,
			   double v, double dVdv[2],SparseGridCluster *sgCluster_);

  struct RiemannInvParams {

    RiemannInvParams(int fluidId,double ent, VarFcn* v) : myFluidId(fluidId), entropy(ent), vf_(v) { }
    int myFluidId;
    double entropy;
    VarFcn* vf_;
  };
  double riemannInvariantKernel1(double density, const RiemannInvParams& J) ;

  double jwlZeroSoundSpeedJwlDensity(const double density, 
                                     const double pressure);
  
  void shockJWL(double phi, double omega,
                double omp1oom, double frho, double frhoi, 
                double frhopi,
                double v, double u, double p, double vi,
                double &ui, double &pi, double &dui, double &dpi);

  // valid for Gas (Perfect and Stiffened) phase
  void rarefactionGAS(double phi,
                   double gam, double gam1, double pref, double c1, 
                   double v1, double u1, double p1,
                   double v, double &u, double &p,
                   double &du, double &dp, int flag = 0);
  void shockGAS(double phi, double gamogam1,
                double pref,
                double v, double u, double p, 
                double vi, double &ui, double &pi,
                double &dui, double &dpi);

  void rarefactionTAIT(double phi,
					 double alpha, double beta, double Pinf,
                   double v1, double u1, double p1,
                   double v, double &u, double &p,
                   double &du, double &dp, int flag);
  void shockTAIT(double phi,
					 double alpha, double beta, double Pinf,
                   double v1, double u1, double p1,
                   double v, double &u, double &p,
                   double &du, double &dp, int flag);

  // function used for the multiphase flow algorithm to update
  // nodes that change phases
  void updatePhaseChangingNodeValues(double * const dx, 
                                     double * const Wi, double * const Wj,
                                     double &weighti, double *rupdatei, 
                                     double &weightj, double *rupdatej);

  // Helper function for Riemann problem jacobians
  void oneDtoThreeD(const double* dWidWi3, const double* dWidWj3,
		    const double* dWjdWi3, const double* dWjdWj3,
		    const double* nphi,
		    double* dWidWi,double*  dWidWj,
		    double* dWjdWi,double*  dWjdWj);

};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

inline
bool LocalRiemannGfmpar::solve2x2System(double *mat, double *rhs, double *res)
{
  double determinant = mat[0]*mat[3]-mat[1]*mat[2];
  double eps = 1.0e-15;
  double norm = fmax(fabs(mat[0])+fabs(mat[2]), fabs(mat[1])+fabs(mat[3]));
  if(fabs(determinant)>eps*norm){
    res[0] = ( mat[3]*rhs[0]-mat[1]*rhs[1])/determinant;
    res[1] = (-mat[2]*rhs[0]+mat[0]*rhs[1])/determinant;
    return true;
  }else{
    fprintf(stdout, "zero-determinant (mat = [%e , %e ; %e , %e] = %e)\n", mat[0],mat[3],mat[1],mat[2],determinant);
    res[0] = 0.0;
    res[1] = 0.0;
    return false;
  }

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

inline
int LocalRiemannGfmpar::rarefactionJWL(double phi,
                   double v1, double u1, double p1,
                   double v,  double &u, double &p,
                   double &du, double &dp, 
                   MultiFluidData::RiemannComputation type, int flag){

  int myFluidId = (phi>=0) ? fluid1 : fluid2;
  double entropy = vf_->computeEntropy(1.0/v1,p1, myFluidId);
  double *in = 0;
  double res1[1] = {0.0};
  int i1=0,i2=1;
  if(type == MultiFluidData::FE){
    in = new double[3];
    in[0] = 1.0/v1; in[1] = entropy; in[2] = 1.0/v;
    riemannInvariantGeneral1stOrder(in,res1,&phi);
  }else if(type == MultiFluidData::RK2){
    in = new double[3];
    in[0] = 1.0/v1; in[1] = entropy; in[2] = 1.0/v;
    riemannInvariantGeneral2ndOrder(in,res1,&phi);
  }else if(type == MultiFluidData::TABULATION2){
    in = new double[3];
    in[0] = 1.0/v1; in[1] = entropy;
    i1 = riemannInvariantGeneralTabulation(in,res1);
    if (!i1) { // Sparse grid failed
      in[2] = 1.0/v;
      fprintf(stderr,"*** Warning: Sparse grid failed on coordinate density = %lf, entropy = %lf; Reverting to 2nd order computation\n", in[0]*refdensity, in[1]*refentropy);
      riemannInvariantGeneral2ndOrder(in,res1,&phi);
    }
  }

  double res2[1] = {0.0};
  if(type == MultiFluidData::TABULATION2 && i1){
    in[0] = 1.0/v;
    i2 = riemannInvariantGeneralTabulation(in,res2);
    if (!i2) { // Sparse grid failed
      in[2] = 1.0/v;
      riemannInvariantGeneral2ndOrder(in,res2,&phi);
      i2 = 1;
    }
  }

  if (in)
    delete [] in;

  u = u1 - phi*(res2[0]-res1[0]);
  p = vf_->computeIsentropicPressure(entropy, 1.0/v, myFluidId);
  double c = vf_->computeSoundSpeed(1.0/v, entropy, myFluidId);
  //fprintf(stdout, "EOS = (%e %e %e %e %e), s=%e, v=%e, myFluidId=%d, c=%e\n", vf_->getOmega(myFluidId), vf_->getA1(myFluidId), vf_->getA2(myFluidId), vf_->getR1r(myFluidId), vf_->getR2r(myFluidId), entropy, v, myFluidId, c);
  du = -phi*c/v;
  dp = -c*c/(v*v);
  if (flag>0 && c<= 0.0) { 
    fprintf(stdout, "*** rarefaction: u1 = %e, v1 = %e, p1 = %e, v = %e\n",u1,v1,p1,v);
    fprintf(stdout, "*** rarefactionJWL returns c=%e, u=%e, p=%e, du=%e, dp=%e, s=%e\n", c,u,p,du,dp,entropy);
  }

  return 1;//(i1&&i2);
}

inline
int LocalRiemannGfmpar::rarefactionJWLderivs(double phi,
					      double v1, double u1, double p1,
					      double v, double dVdv[2],SparseGridCluster *sgCluster_) {

  // dVdv[0] = dV/dp;
  // dVdv[1] = dV/drho;
  int myFluidId = (phi>=0) ? fluid1 : fluid2;
  double entropy = vf_->computeEntropy(1.0/v1,p1, myFluidId);
  double V1[] = {1.0/v1,u1,0.0,0.0,p1};
  double c = vf_->computeSoundSpeed(V1,myFluidId);
  double in[2];
  double grad1[2];
  double* ip = &in[0];
  double* gp = &grad1[0];
  in[0] = 1.0/v1; in[1] = entropy;  
  int i1 = sgCluster_->interpolateGradient(1,&ip,&gp);

  double res2[1] = {0.0};
  double grad2[2];
  gp = &grad2[0];
  in[0] = 1.0/v;
  int i2 = sgCluster_->interpolateGradient(1,&ip,&gp);

  double dsdp = pow(v, vf_->getOmega(myFluidId) + 1.0 );
  double dsdrho = -c*c*dsdp;
  dVdv[0] = (grad2[0]-grad1[0])*dsdp;

  dVdv[1] = (grad2[1]-grad1[1]) + (grad2[0]-grad1[0])*dsdrho;

  return i1 || i2;

}
//----------------------------------------------------------------------------

inline
double LocalRiemannGfmpar::riemannInvariantKernel1(double density, const RiemannInvParams& J) {

  double c = J.vf_->computeSoundSpeed(density, J.entropy, J.myFluidId);
  return -c/density;
}

inline
void LocalRiemannGfmpar::riemannInvariantGeneral1stOrder(double *in, double *res,
                                                   double *phi){
// in contains density and pressure and density
// res is the output result and contains the variation of velocity
// integrates an ODE with first order integration

  ODEIntegrator myIntegrator(in[2],in[0],5000);

  int myFluidId = (*phi>=0) ? fluid1 : fluid2;
  res[0] = 0.0;
  myIntegrator.integrateFE(*this, res[0], &LocalRiemannGfmpar::riemannInvariantKernel1, 
  RiemannInvParams(myFluidId, in[1], vf_));
  /*int N  = 5000;
  double density = in[2]; double entropy = in[1];
  double ddensity = (in[0] - in[2])/N;
  double c = vf_->computeSoundSpeed(density,entropy,myFluidId);

  bool continueCondition = true;
  int it=0;
  while(continueCondition){
    res[0] -= c/density*ddensity;
    density  += ddensity;
    if(ddensity>0.0) density = density>in[0] ? in[0] : density;
    else             density = density<in[0] ? in[0] : density;
    c = vf_->computeSoundSpeed(density,entropy,myFluidId);
    if(ddensity>0.0)
      continueCondition = (density<in[0]-ddensity/2.0);
    else continueCondition = (density>in[0]-ddensity/2.0);
    it++;
    if(it==N) break;
  }
  */
}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpar::riemannInvariantGeneral2ndOrder(double *in, double *res,
                                                   double *phi){
// in contains density and entropy and density
// res is the output result and contains the variation of velocity
// integrates an ODE with second order integration (Mid-Point Rule)

  int myFluidId = (*phi>=0) ? fluid1 : fluid2;
  ODEIntegrator myIntegrator(in[2],in[0],2000);
  res[0] = 0.0;
  myIntegrator.integrateMidpoint(*this, res[0], &LocalRiemannGfmpar::riemannInvariantKernel1, 
  				 RiemannInvParams(myFluidId, in[1], vf_));

  /*int N  = 2000;
  double density = in[2]; double entropy = in[1];
  double ddensity = (in[0] - in[2])/N;
  double c = vf_->computeSoundSpeed(density,entropy,myFluidId);

  bool continueCondition = true;
  int it=0;
  while(continueCondition){
    res[0] -= c/density*ddensity;
    density  += ddensity;
    if(ddensity>0.0) density = density>in[0] ? in[0] : density;
    else             density = density<in[0] ? in[0] : density;
    //advance by second half density-step
    c = vf_->computeSoundSpeed(density,entropy,myFluidId);
    if(ddensity>0.0)
      continueCondition = (density<in[0]-ddensity/2.0);
    else continueCondition = (density>in[0]-ddensity/2.0);
    it++;
    if(it==N) break;
    }*/
  
}

//----------------------------------------------------------------------------

inline
int LocalRiemannGfmpar::riemannInvariantGeneralTabulation(double *in, double *res){
  fprintf(stderr, "*** Error: tabulation of the Riemann invariant is only available for JWL simulation\n");
  exit(1);
}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpar::shockJWL(double phi, double omega,
                   double omp1oom, double frho, double frhoi, 
                   double frhopi,
                   double v, double u, double p, double vi,
                   double &ui, double &pi, double &dui, double &dpi){
//phi = -1 => left
//phi = +1 => right
  double den=omp1oom*vi-0.5*(vi+v);

  pi = (omp1oom*v*p - 0.5*(vi+v)*p
     + (frhoi*vi-frho*v)/omega
       )/den;

  ui = u + phi * sqrt(-(pi-p)*(vi - v));

  dpi = ((frhoi - frhopi/vi)/omega
      - 0.5*p - (omp1oom - 0.5)*pi)/den;

  dui = -0.5*((vi-v)*dpi+pi-p)/(ui-u);

}

//---------------------------------------------------------------------------

inline
void LocalRiemannGfmpar::rarefactionGAS(double phi,
                   double gam, double gam1, double pref, double c1,
                   double v1, double u1, double p1,
                   double v, double &u, double &p,
                   double &du, double &dp, int flag){

  double ppref  = (p1+pref)*pow(v1/v,gam);
  p = ppref - pref;
  double c = sqrt(gam*ppref*v);

  u  = u1 - phi*2.0/gam1*(c1 - c);

  dp = -c*c/(v*v);
  du = -phi*c/v;

  if(flag>0){
    fprintf(stdout, "%e %e %e %e %e %e %e %e %e\n", phi,gam,gam1,pref,c1,1.0/v1,u1,p1,1.0/v);
    fprintf(stdout, "p = %e\n", p);
    fprintf(stdout, "c = %e\n", c);
    fprintf(stdout, "u = %e\n", u);
    fprintf(stdout, "dp = %e\n", dp);
    fprintf(stdout, "du = %e\n", du);
    exit(1);
  }
}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpar::shockGAS(double phi, double gamogam1,
                   double pref,
                   double v, double u, double p, 
                   double vi, double &ui, double &pi,
                   double &dui, double &dpi){
  double den=gamogam1*vi-0.5*(vi+v);

  pi = (gamogam1*v - 0.5*(vi+v))*(p+pref)/den - pref;

  ui = u + phi * sqrt(-(pi-p)*(vi - v));

  dpi = ( -0.5*(p+pref) - (gamogam1 - 0.5)*(pi+pref))/den;

  dui = -0.5*((vi-v)*dpi+pi-p)/(ui-u);
}

//----------------------------------------------------------------------------
inline
void LocalRiemannGfmpar::rarefactionTAIT(double phi,
					 double alpha, double beta, double Pinf,
                   double v1, double u1, double p1,
                   double v, double &u, double &p,
                   double &du, double &dp, int flag){

  double V;
  if (beta != 1.0)
    V = 2.0*sqrt(alpha*beta)/(beta-1.0)*(pow(v1, 0.5*(1.0-beta)) - pow(v, 0.5*(1.0-beta)));
  else
    V = sqrt(alpha*beta)*log(v/v1);

  p = Pinf+alpha*pow(v, -beta);

  u  = u1 - phi*V;

  dp = -alpha*beta*pow(v,-beta-1.0);
  du = -2.0*sqrt(alpha*beta)*0.5*(pow(v, 0.5*(-1.0-beta)));
  if(flag>0){
    //fprintf(stdout, "%e %e %e %e %e %e %e %e %e\n", phi,gam,gam1,pref,c1,1.0/v1,u1,p1,1.0/v);
    fprintf(stdout, "p = %e\n", p);
    //fprintf(stdout, "c = %e\n", c);
    fprintf(stdout, "u = %e\n", u);
    fprintf(stdout, "dp = %e\n", dp);
    fprintf(stdout, "du = %e\n", du);
    exit(1);
  }
}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpar::shockTAIT(double phi,
					 double alpha, double beta, double Pinf,
                   double v, double u, double p, 
                   double vi, double &ui, double &pi,
				   double &dui, double &dpi,int flag){

  //fprintf(stdout,"SHOCK!!\n");
  pi = Pinf+alpha*pow(vi, -beta);
  //fprintf(stdout,"SHOCK2!!\n");

  ui = u + phi * sqrt(max(0.0,-(pi-p)*(vi - v)));
  //fprintf(stdout,"SHOCK3!!\n");

  dpi = -alpha*beta*pow(vi,-beta-1.0);

  //fprintf(stdout,"SHOCK4!!\n");
  double div = ui-u;
  if (fabs(div) < 1.0e-8)
    div = 1.0e-8;
  //fprintf(stdout,"SHOCK5!!\n");
  dui = -0.5*((vi-v)*dpi+pi-p)/div;
  //fprintf(stdout,"SHOCK END!!\n");
}

//----------------------------------------------------------------------------

inline
void LocalRiemannGfmpar::updatePhaseChangingNodeValues( double * const dx, 
                      double * const Wi, double * const Wj,
                      double &weighti, double *rupdatei, 
                      double &weightj, double *rupdatej)
{
// In multiphase flow:
// From one iteration to the next, some nodes change phases.
// The state values at these nodes are not relevant to the thermodynamics
//     of the fluid they belong to at the end of the iteration and thus
//     must somehow be replaced or "updated".
// Appropriate state values are given by the interfacial states of the
//     solution of the two-phase Riemann problem.
// In three dimensions, there are several interfacial states to consider.
// The present routine uses all interfacial states that are upwind of
//    the node that may need to be "updated" and weighs them according
//    to their upwind position.

  int dim = 5;

  double temp = 0.0;
  double normdx2 = dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
  double normWi2 = Wi[1]*Wi[1]+Wi[2]*Wi[2]+Wi[3]*Wi[3];
  double normWj2 = Wj[1]*Wj[1]+Wj[2]*Wj[2]+Wj[3]*Wj[3];

  if(normdx2 > 0.0 && normWj2 > 0.0)
    temp = -(Wj[1]*dx[0]+Wj[2]*dx[1]+Wj[3]*dx[2])/sqrt(normdx2*normWj2);
  if (temp > 0.0){
    weighti += temp;
    for (int k=0; k<dim; k++)
      rupdatei[k] += temp*Wj[k];
  }
  temp = 0.0;
  if(normdx2 > 0.0 && normWi2 > 0.0)
    temp = (Wi[1]*dx[0]+Wi[2]*dx[1]+Wi[3]*dx[2])/sqrt(normdx2*normWi2);
  if(temp > 0.0){ // for update of node j
    weightj += temp;
    for (int k=0; k<dim; k++)
      rupdatej[k] += temp*Wi[k];
  }

}

inline
void LocalRiemannGfmpar::oneDtoThreeD(const double* dWidWi3, const double* dWidWj3,
				      const double* dWjdWi3, const double* dWjdWj3,
				      const double* nphi,
				      double* dWidWi,double*  dWidWj,
				      double* dWjdWi,double*  dWjdWj) {

  int k,l;
  dWidWi[0]  = dWidWi3[0];
  for (k = 0; k < 3; ++k)
    dWidWi[k+1] = dWidWi3[1]*nphi[k];
  dWidWi[4] = dWidWi3[2];
  for (k = 0; k < 3; ++k)
    dWidWi[k*5+5] = dWidWi3[3]*nphi[k];
  for (k = 0; k < 3; ++k) {
    for (l = 0; l < 3; ++l) {
      
      dWidWi[(k+1)*5+(l+1)] = (k==l?1.0:0.0)+(-1.0+dWidWi3[4])*nphi[k]*nphi[l];
    }
  }
  for (k = 0; k < 3; ++k)
    dWidWi[k*5+5+4] = dWidWi3[5]*nphi[k];
  
  dWidWi[20] = dWidWi3[6];
  for (k = 0; k < 3; ++k)
    dWidWi[k+21] = dWidWi3[7]*nphi[k];
  dWidWi[24] = dWidWi3[8];
  
  // Dwi/dwJ
  dWidWj[0]  = dWidWj3[0];
  for (k = 0; k < 3; ++k)
    dWidWj[k+1] = dWidWj3[1]*nphi[k];
  dWidWj[4] = dWidWj3[2];
  for (k = 0; k < 3; ++k)
    dWidWj[k*5+5] = dWidWj3[3]*nphi[k];
  for (k = 0; k < 3; ++k) {
    for (l = 0; l < 3; ++l) {
      
      dWidWj[(k+1)*5+(l+1)] = (dWidWj3[4])*nphi[k]*nphi[l];
    }
  }
  for (k = 0; k < 3; ++k)
    dWidWj[k*5+5+4] = dWidWj3[5]*nphi[k];
  
  dWidWj[20] = dWidWj3[6];
  for (k = 0; k < 3; ++k)
    dWidWj[k+21] = dWidWj3[7]*nphi[k];
  dWidWj[24] = dWidWj3[8];
  
  // DwJ/dWi
  dWjdWi[0]  = dWjdWi3[0];
  for (k = 0; k < 3; ++k)
    dWjdWi[k+1] = dWjdWi3[1]*nphi[k];
  dWjdWi[4] = dWjdWi3[2];
  for (k = 0; k < 3; ++k)
    dWjdWi[k*5+5] = dWjdWi3[3]*nphi[k];
  for (k = 0; k < 3; ++k) {
    for (l = 0; l < 3; ++l) {
      
      dWjdWi[(k+1)*5+(l+1)] = (dWjdWi3[4])*nphi[k]*nphi[l];
    }
  }
  for (k = 0; k < 3; ++k)
    dWjdWi[k*5+5+4] = dWjdWi3[5]*nphi[k];
  
  dWjdWi[20] = dWjdWi3[6];
  for (k = 0; k < 3; ++k)
    dWjdWi[k+21] = dWjdWi3[7]*nphi[k];
  dWjdWi[24] = dWjdWi3[8];
  
  // Dwj/Dwj
  dWjdWj[0]  = dWjdWj3[0];
  for (k = 0; k < 3; ++k)
    dWjdWj[k+1] = dWjdWj3[1]*nphi[k];
  dWjdWj[4] = dWjdWj3[2];
  for (k = 0; k < 3; ++k)
    dWjdWj[k*5+5] = dWjdWj3[3]*nphi[k];
  for (k = 0; k < 3; ++k) {
    for (l = 0; l < 3; ++l) {
      
      dWjdWj[(k+1)*5+(l+1)] = (k==l?1.0:0.0)+(-1.0+dWjdWj3[4])*nphi[k]*nphi[l];
    }
  }
  for (k = 0; k < 3; ++k)
    dWjdWj[k*5+5+4] = dWjdWj3[5]*nphi[k];
  
  dWjdWj[20] = dWjdWj3[6];
  for (k = 0; k < 3; ++k)
    dWjdWj[k+21] = dWjdWj3[7]*nphi[k];
  dWjdWj[24] = dWjdWj3[8];
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


class LocalRiemannLowMach : public LocalRiemann {

protected:

  double beta;

  int dim;

public:
  LocalRiemannLowMach()           { vf_ = 0; fluid1 = 0; fluid2 = 1;}
 LocalRiemannLowMach(VarFcn *vf, int tag1, int tag2, double beta, int dim) : LocalRiemann(vf, tag1, tag2), beta(beta), dim(dim)
  {}
  virtual ~LocalRiemannLowMach()  { vf_ = 0; }

  int updatePhaseChange(double *V, int ID, int IDn, double *newV, double weight,bool isCellCut){
    if(ID == IDn && !isCellCut) return 0; /* node does not change phase: nothing to do*/
    if(weight<=0.0)  { if (IDn >= 0 && ID >= 0) { fprintf(stdout, "*** Error: negative weight in LocalRiemannGfmpar::updatePhaseChange %d %d\n",ID, IDn);
                       return -1;} }
    else
      for(int k=0; k<5; k++) V[k] = newV[k]/weight;

    return 0;
  }
  
  // multiphase Riemann problem
  int computeRiemannSolution(double *Vi, double *Vj,
			     int IDi, int IDj, double *nphi,
			     double *initWi, double *initWj,
			     double *Wi, double *Wj,
			     double *rupdatei, double *rupdatej,
			     double &weighti, double &weightj, 
			     double dx[3], int it, bool isHigherOrder) {

    double aL = vf_->computeSoundSpeed(Vi, IDi);
    double aR = vf_->computeSoundSpeed(Vj, IDj);

    double beta2l = 1.0-beta*beta;
    double beta2r = beta2l;

    double vnj = Vj[1]*nphi[0]+Vj[2]*nphi[1]+Vj[3]*nphi[2];
    double vni = Vi[1]*nphi[0]+Vi[2]*nphi[1]+Vi[3]*nphi[2];
    double vtj[3] = {Vj[1] - vnj*nphi[0], Vj[2] - vnj*nphi[1], Vj[3] - vnj*nphi[2]};
    double vti[3] = {Vi[1] - vni*nphi[0], Vi[2] - vni*nphi[1], Vi[3] - vni*nphi[2]};

    double pl = vf_->getPressure(Vi, IDi);
    double pr = vf_->getPressure(Vj, IDj);
    
    
    double XL = beta2l*beta2l*vni*vni+4.0*beta*beta*aL*aL;
    double XR = beta2r*beta2r*vnj*vnj+4.0*beta*beta*aR*aR;
    
    double C_L = 0.5*Vi[0]*(sqrt(XL) + beta2l*vni);
    double C_R = 0.5*Vj[0]*(sqrt(XR) - beta2r*vnj);
    
    double Csum = C_L+C_R;
    double ustar = (C_L*vni+C_R*vnj)/Csum - (pr-pl)/Csum;
    double pstar = (C_L*pr+C_R*pl)/Csum - C_L*C_R*(vnj-vni)/Csum;
    
    // Use the secant method to find a root for the densities,
    // assuming the flow is isentropic across the two waves.
    double sL = vf_->computeEntropy(Vi[0], pl, IDi);
    double rhonm2 = Vi[0];
    double rhonm1 = Vi[0]*(1.0001);
    /*
    double fnm2 = vf_->computeIsentropicPressure(sL, rhonm2, IDi) - pl;
    double fnm1 = vf_->computeIsentropicPressure(sL, rhonm1, IDi) - pl;
    double rhon;

    while (fabs(rhonm1-rhonm2) > 1.0e-10) {

      rhon = rhonm1 - fnm1*(rhonm1-rhonm2)/(fnm1-fnm2);
      rhonm2 = rhonm1;
      rhonm1 = rhon;

      fnm2 = fnm1;
      fnm1 = vf_->computeIsentropicPressure(sL, rhonm1, IDi) - pl;
    }
    */
    double rhoistar = Vi[0];//rhon;
    /*
    double sR = vf_->computeEntropy(Vj[0], pr, IDj);
    rhonm2 = Vj[0];
    rhonm1 = Vj[0]*(1.0001);
    
    fnm2 = vf_->computeIsentropicPressure(sR, rhonm1, IDj) - pr;
    fnm1 = vf_->computeIsentropicPressure(sR, rhon, IDj) - pr;

    while (fabs(rhon-rhonm1) > 1.0e-10) {

      rhon = rhonm1 - fnm1*(rhonm1-rhonm2)/(fnm1-fnm2);
      rhonm2 = rhonm1;
      rhonm1 = rhon;

      fnm2 = fnm1;
      fnm1 = vf_->computeIsentropicPressure(sR, rhonm1, IDj) - pr;
    }
    */
    double rhojstar = Vj[0];//rhon;

    
    Wi[0]  = rhoistar;                Wi[dim]    = Wi[0];
    Wi[1]  = vti[0]+ustar*nphi[0];      Wi[dim+1]  = Wi[1];
    Wi[2]  = vti[1]+ustar*nphi[1];      Wi[dim+2]  = Wi[2];
    Wi[3]  = vti[2]+ustar*nphi[2];      Wi[dim+3]  = Wi[3];
    Wi[4]  = pstar;                     Wi[dim+4]  = Wi[4];

    Wj[0]  = rhojstar;                    Wj[dim]    = Wj[0];
    Wj[1]  = vtj[0]+ustar*nphi[0];      Wj[dim+1]  = Wj[1];
    Wj[2]  = vtj[1]+ustar*nphi[1];      Wj[dim+2]  = Wj[2];
    Wj[3]  = vtj[2]+ustar*nphi[2];      Wj[dim+3]  = Wj[3];
    Wj[4]  = pstar;                     Wj[dim+4]  = Wj[4];

    if (it == 1 && !isHigherOrder)
      updatePhaseChangingNodeValues(dx, Wi, Wj, weighti, rupdatei, weightj, rupdatej);

    return 0;
  }

  void computeRiemannJacobian(double *Vi, double *Vj,
			      int IDi, int IDj, double *nphi,
			      double *Wi, double *Wj,
			      double dx[3],int it,
			      double* dWidWi,double*  dWidWj,
			      double* dWjdWi,double*  dWjdWj)
  {
    fprintf(stderr,"Error: LocalRiemannLowMach::computeRiemannJacobian is not implemented\n");
    exit(-1);
  } 

  void updatePhaseChangingNodeValues( double * const dx, 
				      double * const Wi, double * const Wj,
				      double &weighti, double *rupdatei, 
				      double &weightj, double *rupdatej)
  {
    // In multiphase flow:
    // From one iteration to the next, some nodes change phases.
    // The state values at these nodes are not relevant to the thermodynamics
    //     of the fluid they belong to at the end of the iteration and thus
    //     must somehow be replaced or "updated".
    // Appropriate state values are given by the interfacial states of the
    //     solution of the two-phase Riemann problem.
    // In three dimensions, there are several interfacial states to consider.
    // The present routine uses all interfacial states that are upwind of
    //    the node that may need to be "updated" and weighs them according
    //    to their upwind position.
    
    int dim = 5;
    
    double temp = 0.0;
    double normdx2 = dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
    double normWi2 = Wi[1]*Wi[1]+Wi[2]*Wi[2]+Wi[3]*Wi[3];
    double normWj2 = Wj[1]*Wj[1]+Wj[2]*Wj[2]+Wj[3]*Wj[3];
    
    if(normdx2 > 0.0 && normWj2 > 0.0)
      temp = -(Wj[1]*dx[0]+Wj[2]*dx[1]+Wj[3]*dx[2])/sqrt(normdx2*normWj2);
    if (temp > 0.0){
      weighti += temp;
      for (int k=0; k<dim; k++)
	rupdatei[k] += temp*Wj[k];
    }
    temp = 0.0;
    if(normdx2 > 0.0 && normWi2 > 0.0)
      temp = (Wi[1]*dx[0]+Wi[2]*dx[1]+Wi[3]*dx[2])/sqrt(normdx2*normWi2);
    if(temp > 0.0){ // for update of node j
      weightj += temp;
      for (int k=0; k<dim; k++)
	rupdatej[k] += temp*Wi[k];
    }
    
  }

};

#endif

