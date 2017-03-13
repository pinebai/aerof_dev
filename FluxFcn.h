#ifndef _FLUX_FCN_H_
#define _FLUX_FCN_H_

#include "FluxFcnBase.h"
#include "FluxFcnDescSG.h"
#include "FluxFcnDescTait.h"
#include "FluxFcnDescJwl.h"

#include <IoData.h>
#include "BcDef.h"
#include <VarFcn.h>
#include <VarFcnJwl.h>
#include <VarFcnSGSA.h>
#include <VarFcnSGKE.h>
#include <VarFcnTait.h>
#include <VarFcnSGEuler.h>

#include <cassert>
#include <cmath>
#include <DebugTools.h>

//--------------------------------------------------------------------------
//
// This class is mostly a collection of FluxFcn used during the simulation
// and is created dynamically according to the needs of the simulation.
// The function of a specific FluxFcn can be accessed if a tag integer
// is provided for multiphase flows. A default tag value of zero
// always calls the member functions of the first FluxFcn.
//
//--------------------------------------------------------------------------
class FluxFcn {

private:
  int numPhases_;
  FluxFcnBase **ff_;

  VarFcn *vf_;
  
  void check(int tag=0) const{ 
    if(tag>=numPhases_){
		fprintf(stdout, "*** Error: An unknown fluid model with FluidID = %d is detected [%d]. Could be a software bug!\n", tag, numPhases_);
      DebugTools::PrintBacktrace();
      fflush(stdout);
      MPI_Barrier(MPI_COMM_WORLD);
      exit(1);
    }
    assert(tag<numPhases_);
  }

  FluxFcnBase *createFluxFcn(int rshift, int ffType, FluidModelData &fmodel, IoData &iod, VarFcnBase *vfb, FluxFcnBase::Type typeJac = FluxFcnBase::CONSERVATIVE, int fid = 0);

  //for Implicit Segregated Navier-Stokes solver ONLY!
  FluxFcnBase *createFluxFcnSeg1(int rshift, int ffType, FluidModelData &fmodel, IoData &iod, FluxFcnBase::Type typeJac = FluxFcnBase::CONSERVATIVE);
  FluxFcnBase *createFluxFcnSeg2(int rshift, int ffType, FluidModelData &fmodel, IoData &iod, VarFcnBase *vfb, FluxFcnBase::Type typeJac = FluxFcnBase::CONSERVATIVE);

public:
  FluxFcn(int rshift, int ffType, IoData &iod, VarFcn *vf, FluxFcnBase::Type typeJac = FluxFcnBase::CONSERVATIVE); 
  //for Implicit Segregated Navier-Stokes solver ONLY! (segPart = 1 or 2)
  FluxFcn(int rshift, int ffType, IoData &iod, VarFcn *vf, int segPart, FluxFcnBase::Type typeJac = FluxFcnBase::CONSERVATIVE); 

  FluxFcn() {}
  ~FluxFcn() {
    for(int i=0; i<numPhases_; i++) {
      if(ff_[i]) delete ff_[i];
    }
    delete [] ff_;
  }

  VarFcn *getVarFcn() { return vf_; }

  void setHHCoeffPointer(double* hh) { 

    for (int i = 0; i < numPhases_; i++)
       ff_[i]->setHHCoeffPointer(hh);
  }

  //----- General Functions -----//
  void compute(double length, double irey, double *normal, double normalVel, double *VL, double *VR, double *flux, int tag=0, bool useLimiter = true){
    check(tag);
    ff_[tag]->compute(length, irey, normal, normalVel, VL, VR, flux, useLimiter);
  }
  void computeJacobian(double length, double irey, double *normal, double normalVel, double *VL, double *VR, double *jacL, int tag=0, bool useLimiter = true){
    check(tag);
    ff_[tag]->computeJacobian(length, irey, normal, normalVel, VL, VR, jacL, useLimiter);
  }
  void computeJacobians(double length, double irey, double *normal, double normalVel, double *VL, double *VR, double *jacL, double *jacR, int tag=0, bool useLimiter = true){
    check(tag);
    ff_[tag]->computeJacobians(length, irey, normal, normalVel, VL, VR, jacL, jacR, useLimiter);
  }

  void computeJacobianFarfield(double length, double irey, double *normal, double normalVel, double * VL, double * Ub, double * jac,
			       int tag = 0,bool useLimiter = true) {

    check(tag);
    ff_[tag]->computeJacobianFarfield(length,irey,normal,normalVel,VL,Ub,jac,useLimiter);
  }

  //----- Sensitivity Analysis Functions -----//
  void computeDerivative
  (
     double ire, double dIre, double *n, double *dn, double nv, double dnv, 
     double *vl, double *dvl, double *vr, double *dvr, double dmach, 
     double *f, double *df, int tag=0
  )
  {
    assert(numPhases_==1);
    check(tag);
    ff_[tag]->computeDerivative(ire,dIre,n,dn,nv,dnv,vl,dvl,vr,dvr,dmach,f,df);
  }

  void compute_dFluxdNormal_dFluxdNormalVel_dFluxdVL_dFluxdVR(double *n, double nv, double *vl, double *vr,
                                                              double dmach, double *f, double dfdn[7][3],
                                                              double *fdnv, double dfdvl[7][7], double dfdvr[7][7], int tag=0) 
  {
    assert(numPhases_==1);
    check(tag);
    ff_[tag]->compute_dFluxdNormal_dFluxdNormalVel_dFluxdVL_dFluxdVR(n, nv, vl, vr, dmach, f, dfdn, fdnv, dfdvl, dfdvr);
  }

  void computeDerivative
  (
    double ire, double dIre, double *n, double *dn, double nv, double dnv, 
    double *v, double *ub, double *dub, double *f, double *df,
    int tag=0
  )
  {
    assert(numPhases_==1);
    check(tag);
    ff_[tag]->computeDerivative(ire,dIre,n,dn,nv,dnv,v,ub,dub,f,df);
  }

  void computeDerivativeOperators(double ire, double dIre, double *n, double nv, double *v, double *ub,
                                  double dfdn[7][3], double dfdsig[7][1], double dfdub[7][7], int tag=0)
  {
    assert(numPhases_==1);
    check(tag);
    ff_[tag]->computeDerivativeOperators(n, nv, v, ub, dfdn, dfdsig, dfdub);
  }


  FluxFcnBase* getFluxFcnBase(int tag = 0) const { check(tag); return ff_[tag]; }

};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

inline
FluxFcn::FluxFcn(int rshift, int ffType, IoData &iod, VarFcn *vf, FluxFcnBase::Type typeJac) : vf_(vf){

  numPhases_ = iod.eqs.numPhase;
  ff_ = new FluxFcnBase *[numPhases_];

  if(vf_==0 || vf_->varFcn==0){
    fprintf(stderr, "*** Error: VarFcn is NULL in FluxFcn constructor\n");
    exit(1);
  }
  if(numPhases_>0){
    if(vf_->varFcn[0] == 0){
      fprintf(stderr, "*** Error: member varFcn[0] of VarFcn has not been initialized correctly\n");
      exit(1);
    }
    ff_[0] = createFluxFcn(rshift, ffType, iod.eqs.fluidModel, iod, vf_->varFcn[0], typeJac,0);
  }

  for(int iPhase=1; iPhase<numPhases_; iPhase++){
    map<int, FluidModelData *>::iterator it = iod.eqs.fluidModelMap.dataMap.find(iPhase);
    if(it == iod.eqs.fluidModelMap.dataMap.end()){
      fprintf(stderr, "*** Error: no FluidModel[%d] was specified\n", iPhase);
      exit(1);
    }
    ff_[iPhase] = createFluxFcn(rshift, ffType, *it->second, iod, vf_->varFcn[iPhase], typeJac,
                                iPhase);
  }

}

//------------------------------------------------------------------------------

inline
FluxFcn::FluxFcn(int rshift, int ffType, IoData &iod, VarFcn *vf, int segPart, FluxFcnBase::Type typeJac) : vf_(vf){

  numPhases_ = iod.eqs.numPhase;
  ff_ = new FluxFcnBase *[numPhases_];

  if(vf_==0 || vf_->varFcn==0 || vf_->varFcn[0]==0 || numPhases_!=1 || ((segPart!=1)&&(segPart!=2)) ){
    fprintf(stderr, "*** Error: Unable to construct FluxFcn for Implicit Segregated NS solver!\n");
    exit(1);
  }
  if(segPart==1)  
    ff_[0] = createFluxFcnSeg1(rshift, ffType, iod.eqs.fluidModel, iod, typeJac); //an Euler (inviscid) VarFcn will be created
  if(segPart==2) 
    ff_[0] = createFluxFcnSeg2(rshift, ffType, iod.eqs.fluidModel, iod, vf_->varFcn[0], typeJac);
}

//------------------------------------------------------------------------------

inline
FluxFcnBase *FluxFcn::createFluxFcn(int rshift, int ffType, FluidModelData &fmodel, IoData &iod, 
                                    VarFcnBase *vfb, FluxFcnBase::Type typeJac,
                                    int iPhase){

  FluxFcnBase *localff = 0;
  double gamma = iod.schemes.ns.gamma;

  SchemeData::Flux ns_flux = iod.schemes.ns.flux; 

  // Check to see if the user has specified a particular flux for this material
  if (iod.schemes.ns.fluxMap.dataMap.find(iPhase) != iod.schemes.ns.fluxMap.dataMap.end()) {
    ns_flux = iod.schemes.ns.fluxMap.dataMap.find(iPhase)->second->flux;
    //fprintf(stderr,"*** using flux = %d for iPhase = %d\n", ns_flux, iPhase); 
  }

  if(fmodel.fluid == FluidModelData::PERFECT_GAS || fmodel.fluid == FluidModelData::STIFFENED_GAS){
    if (iod.eqs.type == EquationsData::NAVIER_STOKES &&
        iod.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
      if (iod.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
          iod.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES) {
        VarFcnSGSA *vfsgsa = dynamic_cast<VarFcnSGSA *>(vfb);
        if(vfsgsa == 0){
          fprintf(stderr, "*** Error: a VarFcnSGSA is expected to create the associated FluxFcn\n");
          exit(-1);
        }
// Spalart-Allmaras for Stiffened Gas
        switch(ffType){

          case BC_DIRECTSTATE_OUTLET_FIXED:
          case BC_DIRECTSTATE_OUTLET_MOVING:
            localff = new FluxFcnSGDirectStateOutflowSA3D(iod, vfsgsa, typeJac);
            break;

          case BC_DIRECTSTATE_INLET_FIXED:
          case BC_DIRECTSTATE_INLET_MOVING:
            localff = new FluxFcnSGDirectStateInflowSA3D(iod, vfsgsa, typeJac);
            break;

          case BC_MASSFLOW_OUTLET_FIXED:
          case BC_MASSFLOW_OUTLET_MOVING:
            localff = new FluxFcnSGMassFlowOutflowSA3D(iod, vfsgsa, typeJac);
            break;

          case BC_MASSFLOW_INLET_FIXED:
          case BC_MASSFLOW_INLET_MOVING:
            localff = new FluxFcnSGMassFlowInflowSA3D(iod, vfsgsa, typeJac);
            break;

          case BC_OUTLET_FIXED:
          case BC_OUTLET_MOVING:
            if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL &&
              iod.schemes.bc.type == BoundarySchemeData::STEGER_WARMING)
              localff = new FluxFcnSGOutflowSA3D(iod, vfsgsa, typeJac);
            else if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL &&
              iod.schemes.bc.type == BoundarySchemeData::GHIDAGLIA)
              localff = new FluxFcnSGGhidagliaSA3D(iod, vfsgsa, typeJac);
            else
              localff = new FluxFcnSGInternalOutflowSA3D(iod, vfsgsa, typeJac);
            break;
      
          case BC_INLET_FIXED:
          case BC_INLET_MOVING:
            if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL &&
              iod.schemes.bc.type == BoundarySchemeData::STEGER_WARMING)
              localff = new FluxFcnSGOutflowSA3D(iod, vfsgsa, typeJac);
            else if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL &&
              iod.schemes.bc.type == BoundarySchemeData::GHIDAGLIA)
              localff = new FluxFcnSGGhidagliaSA3D(iod, vfsgsa, typeJac);
            else
              localff = new FluxFcnSGInternalInflowSA3D(iod, vfsgsa, typeJac);
            break;

          case BC_ADIABATIC_WALL_MOVING:
          case BC_ADIABATIC_WALL_FIXED:
          case BC_SLIP_WALL_MOVING:
          case BC_SLIP_WALL_FIXED:
          case BC_SYMMETRY:
          case BC_ISOTHERMAL_WALL_MOVING:
          case BC_ISOTHERMAL_WALL_FIXED:
            localff = new FluxFcnSGWallSA3D(iod, vfsgsa, typeJac);
            break;

          case BC_POROUS_WALL_MOVING:
          case BC_POROUS_WALL_FIXED:
            localff = new FluxFcnSGPorousWallSA3D(iod, vfsgsa, typeJac);
            break;

          case BC_INTERNAL:
            if (ns_flux == SchemeData::ROE) {
              if (iod.ts.implicit.ffjacobian == ImplicitData::FINITE_DIFFERENCE)
                localff = new FluxFcnSGFDJacRoeSA3D(gamma, iod, vfsgsa, typeJac);
              else if (iod.ts.implicit.ffjacobian == ImplicitData::APPROXIMATE)
                localff = new FluxFcnSGApprJacRoeSA3D(rshift, gamma, iod, vfsgsa, typeJac);
              else if (iod.ts.implicit.ffjacobian == ImplicitData::EXACT)
                localff = new FluxFcnSGExactJacRoeSA3D(gamma, iod, vfsgsa, typeJac);
            }
            else if (ns_flux == SchemeData::HLLE) {
              if (iod.ts.implicit.ffjacobian == ImplicitData::FINITE_DIFFERENCE)
                localff = new FluxFcnSGFDJacHLLESA3D(gamma, iod, vfsgsa, typeJac);
              else if (iod.ts.implicit.ffjacobian == ImplicitData::APPROXIMATE) {
                localff = new FluxFcnSGApprJacHLLESA3D(rshift, gamma, iod, vfsgsa, typeJac);
              }
              else if (iod.ts.implicit.ffjacobian == ImplicitData::EXACT) {
                fprintf(stderr,"Error... HLLE with Exact Jacobian not Implemented.. Aborting !!");
                exit(1);
              }
            }
            else if (ns_flux == SchemeData::HLLC) {
              if (iod.ts.implicit.ffjacobian == ImplicitData::FINITE_DIFFERENCE)
                localff = new FluxFcnSGFDJacHLLCSA3D(gamma, iod, vfsgsa, typeJac);
              else if (iod.ts.implicit.ffjacobian == ImplicitData::APPROXIMATE) {
                localff = new FluxFcnSGApprJacHLLCSA3D(rshift, gamma, iod, vfsgsa, typeJac);
	      }
              else if (iod.ts.implicit.ffjacobian == ImplicitData::EXACT) {
                fprintf(stderr,"Error... HLLC with Exact Jacobian not Implemented.. Aborting !!");
                exit(1);
              }
            }
            break;

        }
      }
      else if (iod.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE) {
// k-epsilon turbulent model for Stiffened Gas
        VarFcnSGKE *vfsgke = dynamic_cast<VarFcnSGKE *>(vfb);
        if(vfsgke == 0){
          fprintf(stderr, "*** Error: a VarFcnSGKE is expected to create the associated FluxFcn\n");
          exit(-1);
        }
        switch(ffType){

          case BC_DIRECTSTATE_OUTLET_FIXED:
          case BC_DIRECTSTATE_OUTLET_MOVING:
            localff = new FluxFcnSGDirectStateOutflowKE3D(iod, vfsgke, typeJac);
            break;

          case BC_DIRECTSTATE_INLET_FIXED:
          case BC_DIRECTSTATE_INLET_MOVING:
            localff = new FluxFcnSGDirectStateInflowKE3D(iod, vfsgke, typeJac);
            break;

          case BC_MASSFLOW_OUTLET_FIXED:
          case BC_MASSFLOW_OUTLET_MOVING:
            localff = new FluxFcnSGMassFlowOutflowKE3D(iod, vfsgke, typeJac);
            break;

          case BC_MASSFLOW_INLET_FIXED:
          case BC_MASSFLOW_INLET_MOVING:
            localff = new FluxFcnSGMassFlowInflowKE3D(iod, vfsgke, typeJac);
            break;

          case BC_OUTLET_FIXED:
          case BC_OUTLET_MOVING:
          case BC_INLET_FIXED:
          case BC_INLET_MOVING:
            if(iod.schemes.bc.type == BoundarySchemeData::STEGER_WARMING)
              localff = new FluxFcnSGOutflowKE3D(iod, vfsgke, typeJac);
            else 
              localff = new FluxFcnSGGhidagliaKE3D(iod, vfsgke, typeJac);
            break;
  
          case BC_ADIABATIC_WALL_MOVING:
          case BC_ADIABATIC_WALL_FIXED:
          case BC_SLIP_WALL_MOVING:
          case BC_SLIP_WALL_FIXED:
          case BC_SYMMETRY:
          case BC_ISOTHERMAL_WALL_MOVING:
          case BC_ISOTHERMAL_WALL_FIXED:
            localff = new FluxFcnSGWallKE3D(iod, vfsgke, typeJac);
            break;
  
          case BC_POROUS_WALL_MOVING:
          case BC_POROUS_WALL_FIXED:
            localff = new FluxFcnSGPorousWallKE3D(iod, vfsgke, typeJac);
            break;

          case BC_INTERNAL:
            if (ns_flux == SchemeData::ROE) {
              if (iod.ts.implicit.ffjacobian == ImplicitData::FINITE_DIFFERENCE)
                localff = new FluxFcnSGFDJacRoeKE3D(gamma, iod, vfsgke, typeJac);
              else if (iod.ts.implicit.ffjacobian == ImplicitData::APPROXIMATE)
                localff = new FluxFcnSGApprJacRoeKE3D(rshift, gamma, iod, vfsgke, typeJac);
              else if (iod.ts.implicit.ffjacobian == ImplicitData::EXACT)
                localff = new FluxFcnSGExactJacRoeKE3D(gamma, iod, vfsgke, typeJac);
            }
            else if (ns_flux == SchemeData::HLLE) {
              if (iod.ts.implicit.ffjacobian == ImplicitData::FINITE_DIFFERENCE)
                localff = new FluxFcnSGFDJacHLLEKE3D(gamma, iod, vfsgke, typeJac);
              else if (iod.ts.implicit.ffjacobian == ImplicitData::APPROXIMATE) {
                localff = new FluxFcnSGApprJacHLLEKE3D(rshift, gamma, iod, vfsgke, typeJac);
              }
              else if (iod.ts.implicit.ffjacobian == ImplicitData::EXACT) {
                fprintf(stderr,"Error... HLLE with Exact Jacobian not Implemented.. Aborting !!");
                exit(1);
              }
            }
            else if (ns_flux == SchemeData::HLLC) {
              if (iod.ts.implicit.ffjacobian == ImplicitData::FINITE_DIFFERENCE)
                localff = new FluxFcnSGFDJacHLLCKE3D(gamma, iod, vfsgke, typeJac);
              else if (iod.ts.implicit.ffjacobian == ImplicitData::APPROXIMATE) {
                localff = new FluxFcnSGApprJacHLLCKE3D(rshift, gamma, iod, vfsgke, typeJac);
              }
              else if (iod.ts.implicit.ffjacobian == ImplicitData::EXACT) {
                fprintf(stderr,"Error... HLLC with Exact Jacobian not Implemented.. Aborting !!");
                exit(1);
              }
            }
            break;

        }
      } // end - k-epsilon turbulent model for Stiffened Gas
    } // end - turbulence
    else{
// Euler or Navier-Stokes for Stiffened Gas 
      VarFcnSGEuler *vfsgeuler = dynamic_cast<VarFcnSGEuler *>(vfb);
      if(vfsgeuler == 0){
        fprintf(stderr, "*** Error: a VarFcnSGEuler is expected to create the associated FluxFcn\n");
        exit(-1);
      }
      switch(ffType){

        case BC_DIRECTSTATE_OUTLET_FIXED:
        case BC_DIRECTSTATE_OUTLET_MOVING:
          localff = new FluxFcnSGDirectStateOutflowEuler3D(iod, vfsgeuler, typeJac);
          break;

        case BC_DIRECTSTATE_INLET_FIXED:
        case BC_DIRECTSTATE_INLET_MOVING:
          localff = new FluxFcnSGDirectStateInflowEuler3D(iod, vfsgeuler, typeJac);
          break;

        case BC_MASSFLOW_OUTLET_FIXED:
        case BC_MASSFLOW_OUTLET_MOVING:
          localff = new FluxFcnSGMassFlowOutflowEuler3D(iod, vfsgeuler, typeJac);
          break;

        case BC_MASSFLOW_INLET_FIXED:
        case BC_MASSFLOW_INLET_MOVING:
          localff = new FluxFcnSGMassFlowInflowEuler3D(iod, vfsgeuler, typeJac);
          break;

        case BC_OUTLET_FIXED:
        case BC_OUTLET_MOVING:
          if (iod.bc.outlet.type == BcsFreeStreamData::INTERNAL)
            localff = new FluxFcnSGInternalOutflowEuler3D(iod, vfsgeuler, typeJac);
          else if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL &&
                  iod.schemes.bc.type == BoundarySchemeData::STEGER_WARMING)
            localff = new FluxFcnSGOutflowEuler3D(iod, vfsgeuler, typeJac);
          else if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL &&
                  iod.schemes.bc.type == BoundarySchemeData::GHIDAGLIA)
            localff = new FluxFcnSGGhidagliaEuler3D(iod, vfsgeuler, typeJac);
          else if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL &&
                  iod.schemes.bc.type == BoundarySchemeData::MODIFIED_GHIDAGLIA)
            localff = new FluxFcnSGModifiedGhidagliaEuler3D(iod, vfsgeuler, typeJac);
          else{
            fprintf(stderr, "*** Error: no outlet boundary flux has been selected for Stiffened Gas\n");
            exit(-1);
          }
          break;
     
        case BC_INLET_FIXED:
        case BC_INLET_MOVING:
          if (iod.bc.inlet.type == BcsFreeStreamData::INTERNAL)
            localff = new FluxFcnSGInternalInflowEuler3D(iod, vfsgeuler, typeJac);
          else if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL &&
                  iod.schemes.bc.type == BoundarySchemeData::STEGER_WARMING)
            localff = new FluxFcnSGInflowEuler3D(iod, vfsgeuler, typeJac);
          else if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL &&
                  iod.schemes.bc.type == BoundarySchemeData::GHIDAGLIA)
            localff = new FluxFcnSGGhidagliaEuler3D(iod, vfsgeuler, typeJac);
          else if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL &&
                  iod.schemes.bc.type == BoundarySchemeData::MODIFIED_GHIDAGLIA)
            localff = new FluxFcnSGModifiedGhidagliaEuler3D(iod, vfsgeuler, typeJac);
          else{
            fprintf(stderr, "*** Error: no inlet boundary flux has been selected for Stiffened Gas\n");
            exit(-1);
          }
          break;

        case BC_ADIABATIC_WALL_MOVING:
        case BC_ADIABATIC_WALL_FIXED:
        case BC_SLIP_WALL_MOVING:
        case BC_SLIP_WALL_FIXED:
        case BC_SYMMETRY:
        case BC_ISOTHERMAL_WALL_MOVING:
        case BC_ISOTHERMAL_WALL_FIXED:
          localff = new FluxFcnSGWallEuler3D(iod, vfsgeuler, typeJac);
          break;

        case BC_POROUS_WALL_MOVING:
        case BC_POROUS_WALL_FIXED:
          localff = new FluxFcnSGPorousWallEuler3D(iod, vfsgeuler, typeJac);
          break;

        case BC_INTERNAL:
          if (ns_flux == SchemeData::VANLEER)
          {
            localff = new FluxFcnSGVanLeerEuler3D(iod, vfsgeuler, typeJac);
          }
          else if (ns_flux == SchemeData::ROE) 
          {
            if (iod.ts.implicit.ffjacobian == ImplicitData::FINITE_DIFFERENCE)
            {
              localff = new FluxFcnSGFDJacRoeEuler3D(gamma, iod, vfsgeuler, typeJac);
            }
            else if (iod.ts.implicit.ffjacobian == ImplicitData::APPROXIMATE)
            {
              localff = new FluxFcnSGApprJacRoeEuler3D(rshift, gamma, iod, vfsgeuler, typeJac);
            }
            else if (iod.ts.implicit.ffjacobian == ImplicitData::EXACT)
            {
              localff = new FluxFcnSGExactJacRoeEuler3D(gamma, iod, vfsgeuler, typeJac);
            }
          }
          else if (ns_flux == SchemeData::HLLE) 
          {
            if (iod.ts.implicit.ffjacobian == ImplicitData::FINITE_DIFFERENCE)
              localff = new FluxFcnSGFDJacHLLEEuler3D(gamma, iod, vfsgeuler, typeJac);
            else if (iod.ts.implicit.ffjacobian == ImplicitData::APPROXIMATE) {
              localff = new FluxFcnSGApprJacHLLEEuler3D(rshift, gamma, iod, vfsgeuler, typeJac);
            }
            else if (iod.ts.implicit.ffjacobian == ImplicitData::EXACT) {
              fprintf(stderr,"Error... HLLE with Exact Jacobian not Implemented.. Aborting !!");
              exit(1);
            }
          }
	  else if (ns_flux == SchemeData::HLLC) 
          {
            if (iod.ts.implicit.ffjacobian == ImplicitData::FINITE_DIFFERENCE)
              localff = new FluxFcnSGFDJacHLLCEuler3D(gamma, iod, vfsgeuler, typeJac);
            else if (iod.ts.implicit.ffjacobian == ImplicitData::APPROXIMATE) {
              localff = new FluxFcnSGApprJacHLLCEuler3D(rshift, gamma, iod, vfsgeuler, typeJac);
            }
            else if (iod.ts.implicit.ffjacobian == ImplicitData::EXACT) {
              fprintf(stderr,"Error... HLLC with Exact Jacobian not Implemented.. Aborting !!");
              exit(1);
            }
          }
          break;

      }
    } // end - Euler or Navier-Stokes for Stiffened Gas 
  } // end - All Stiffened Gas
  else if (fmodel.fluid == FluidModelData::LIQUID){
    if (iod.eqs.type == EquationsData::NAVIER_STOKES &&
        iod.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
      if (iod.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
          iod.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES) {
        VarFcnTaitSA *vftaitsa = dynamic_cast<VarFcnTaitSA *>(vfb);
        if(vftaitsa == 0){
          fprintf(stderr, "*** Error: a VarFcnTaitSA is expected to create the associated FluxFcn\n");
          exit(-1);
        }
// Spalart-Allmaras for Barotropic Liquids 
        switch(ffType){

          case BC_DIRECTSTATE_OUTLET_FIXED:
          case BC_DIRECTSTATE_OUTLET_MOVING:
          case BC_DIRECTSTATE_INLET_FIXED:
          case BC_DIRECTSTATE_INLET_MOVING:
          case BC_MASSFLOW_OUTLET_FIXED:
          case BC_MASSFLOW_OUTLET_MOVING:
          case BC_MASSFLOW_INLET_FIXED:
          case BC_MASSFLOW_INLET_MOVING:
            break;

          case BC_OUTLET_FIXED:
          case BC_OUTLET_MOVING:
            if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL)
              localff = new FluxFcnTaitGhidagliaSA3D(iod, vftaitsa, typeJac);
            else{
              fprintf(stderr,"Error... Internal flow boundary conditions are not available for a turbulent (SA) barotropic liquid. Aborting !!");
              exit(1);
	    }
            break;
      
          case BC_INLET_FIXED:
          case BC_INLET_MOVING:
            if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL)
              localff = new FluxFcnTaitGhidagliaSA3D(iod, vftaitsa, typeJac);
            else{
              fprintf(stderr,"Error... Internal flow boundary conditions are not available for a turbulent (SA) barotropic liquid. Aborting !!");
              exit(1);
	    }
            break;

          case BC_ADIABATIC_WALL_MOVING:
          case BC_ADIABATIC_WALL_FIXED:
          case BC_SLIP_WALL_MOVING:
          case BC_SLIP_WALL_FIXED:
          case BC_SYMMETRY:
          case BC_ISOTHERMAL_WALL_MOVING:
          case BC_ISOTHERMAL_WALL_FIXED:
          case BC_POROUS_WALL_MOVING:
          case BC_POROUS_WALL_FIXED:
            localff = new FluxFcnTaitWallSA3D(iod, vftaitsa, typeJac);
            break;

          case BC_INTERNAL:
            if (ns_flux == SchemeData::ROE) {
              localff = new FluxFcnTaitApprJacRoeSA3D(rshift, gamma, iod, vftaitsa, typeJac);
            }
            else{
              fprintf(stderr,"Error... only the Roe flux is available for a turbulent (SA) barotropic liquid. Aborting !!");
              exit(1);
            }
            break;

        }
      }
      else if (iod.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE) {
// k-epsilon turbulent model for Barotropic Liquids
        VarFcnTaitKE *vftaitke = dynamic_cast<VarFcnTaitKE *>(vfb);
        if(vftaitke == 0){
          fprintf(stderr, "*** Error: a VarFcnTaitKE is expected to create the associated FluxFcn\n");
          exit(-1);
        }
        switch(ffType){

          case BC_DIRECTSTATE_OUTLET_FIXED:
          case BC_DIRECTSTATE_OUTLET_MOVING:
          case BC_DIRECTSTATE_INLET_FIXED:
          case BC_DIRECTSTATE_INLET_MOVING:
          case BC_MASSFLOW_OUTLET_FIXED:
          case BC_MASSFLOW_OUTLET_MOVING:
          case BC_MASSFLOW_INLET_FIXED:
          case BC_MASSFLOW_INLET_MOVING:
            break;

          case BC_OUTLET_FIXED:
          case BC_OUTLET_MOVING:
            if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL)
              localff = new FluxFcnTaitGhidagliaKE3D(iod, vftaitke, typeJac);
            else{
              fprintf(stderr,"Error... Internal flow boundary conditions are not available for a turbulent (KE) barotropic liquid. Aborting !!");
              exit(1);
	    }
            break;
      
          case BC_INLET_FIXED:
          case BC_INLET_MOVING:
            if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL)
              localff = new FluxFcnTaitGhidagliaKE3D(iod, vftaitke, typeJac);
            else{
              fprintf(stderr,"Error... Internal flow boundary conditions are not available for a turbulent (KE) barotropic liquid. Aborting !!");
              exit(1);
	    }
            break;

          case BC_ADIABATIC_WALL_MOVING:
          case BC_ADIABATIC_WALL_FIXED:
          case BC_SLIP_WALL_MOVING:
          case BC_SLIP_WALL_FIXED:
          case BC_SYMMETRY:
          case BC_ISOTHERMAL_WALL_MOVING:
          case BC_ISOTHERMAL_WALL_FIXED:
          case BC_POROUS_WALL_MOVING:
          case BC_POROUS_WALL_FIXED:
            localff = new FluxFcnTaitWallKE3D(iod, vftaitke, typeJac);
            break;

          case BC_INTERNAL:
            if (ns_flux == SchemeData::ROE) {
              localff = new FluxFcnTaitApprJacRoeKE3D(rshift, gamma, iod, vftaitke, typeJac);
            }
            else{
              fprintf(stderr,"Error... only the Roe flux is available for a turbulent (KE) barotropic liquid. Aborting !!");
              exit(1);
            }
            break;

        }
      } // end - k-epsilon turbulent model for Stiffened Gas
    } // end - turbulence
    else{
// Euler or Navier-Stokes for Stiffened Gas 
      VarFcnTait *vftait = dynamic_cast<VarFcnTait *>(vfb);
      if(vftait == 0){
        fprintf(stderr, "*** Error: a VarFcnTait is expected to create the associated FluxFcn\n");
        exit(-1);
      }
      switch(ffType){

        case BC_DIRECTSTATE_OUTLET_FIXED:
        case BC_DIRECTSTATE_OUTLET_MOVING:
        case BC_DIRECTSTATE_INLET_FIXED:
        case BC_DIRECTSTATE_INLET_MOVING:
        case BC_MASSFLOW_OUTLET_FIXED:
        case BC_MASSFLOW_OUTLET_MOVING:
        case BC_MASSFLOW_INLET_FIXED:
        case BC_MASSFLOW_INLET_MOVING:
          break;

        case BC_OUTLET_FIXED:
        case BC_OUTLET_MOVING:
          if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL &&
             iod.schemes.bc.type == BoundarySchemeData::GHIDAGLIA)
            localff = new FluxFcnTaitGhidagliaEuler3D(iod, vftait, typeJac);
          else if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL &&
                  iod.schemes.bc.type == BoundarySchemeData::MODIFIED_GHIDAGLIA)
            localff = new FluxFcnTaitModifiedGhidagliaEuler3D(iod, vftait, typeJac);
          else if (iod.bc.outlet.type == BcsFreeStreamData::INTERNAL)
            localff = new FluxFcnTaitInternalOutflowEuler3D(iod, vftait, typeJac);
          else{
            fprintf(stderr, "*** Error: no outlet boundary flux has been selected for Tait\n");
            exit(-1);
          }
          break;
      
        case BC_INLET_FIXED:
        case BC_INLET_MOVING:
          if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL &&
             iod.schemes.bc.type == BoundarySchemeData::GHIDAGLIA)
            localff = new FluxFcnTaitGhidagliaEuler3D(iod, vftait, typeJac);
          else if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL &&
                  iod.schemes.bc.type == BoundarySchemeData::MODIFIED_GHIDAGLIA)
            localff = new FluxFcnTaitModifiedGhidagliaEuler3D(iod, vftait, typeJac);
          else if (iod.bc.inlet.type == BcsFreeStreamData::INTERNAL)
            localff = new FluxFcnTaitInternalInflowEuler3D(iod, vftait, typeJac);
          else{
            fprintf(stderr, "*** Error: no inlet boundary flux has been selected for Tait\n");
            exit(-1);
          }
          break;
  
        case BC_ADIABATIC_WALL_MOVING:
        case BC_ADIABATIC_WALL_FIXED:
        case BC_SLIP_WALL_MOVING:
        case BC_SLIP_WALL_FIXED:
        case BC_SYMMETRY:
        case BC_ISOTHERMAL_WALL_MOVING:
        case BC_ISOTHERMAL_WALL_FIXED:
        case BC_POROUS_WALL_MOVING:
        case BC_POROUS_WALL_FIXED:
          localff = new FluxFcnTaitWallEuler3D(iod, vftait, typeJac);
          break;
  
        case BC_INTERNAL:
          if (ns_flux == SchemeData::ROE &&
              iod.ts.implicit.ffjacobian == ImplicitData::APPROXIMATE)
            localff = new FluxFcnTaitApprJacRoeEuler3D(rshift, gamma, iod, vftait, typeJac);
          else{
            fprintf(stderr, "*** Error: only the Roe flux is available for Tait\n");
            exit(-1);
          }
          break;

      }
    } // end - Euler or Navier-Stokes for Tait EOS
  }
  else if (fmodel.fluid == FluidModelData::JWL){
// Euler or Navier-Stokes for JWL EOS
    VarFcnJwl *vfjwl = dynamic_cast<VarFcnJwl *>(vfb);
    if(vfjwl == 0){
      fprintf(stderr, "*** Error: a VarFcnJwl is expected to create the associated FluxFcn\n");
      exit(-1);
    }
    switch(ffType){

      case BC_DIRECTSTATE_OUTLET_FIXED:
      case BC_DIRECTSTATE_OUTLET_MOVING:
      case BC_DIRECTSTATE_INLET_FIXED:
      case BC_DIRECTSTATE_INLET_MOVING:
      case BC_MASSFLOW_OUTLET_FIXED:
      case BC_MASSFLOW_OUTLET_MOVING:
      case BC_MASSFLOW_INLET_FIXED:
      case BC_MASSFLOW_INLET_MOVING:
        break;

      case BC_OUTLET_FIXED:
      case BC_OUTLET_MOVING:
      case BC_INLET_FIXED:
      case BC_INLET_MOVING:
        localff = new FluxFcnJwlGhidagliaEuler3D(iod, vfjwl, typeJac);
        break;

      case BC_ADIABATIC_WALL_MOVING:
      case BC_ADIABATIC_WALL_FIXED:
      case BC_SLIP_WALL_MOVING:
      case BC_SLIP_WALL_FIXED:
      case BC_SYMMETRY:
      case BC_ISOTHERMAL_WALL_MOVING:
      case BC_ISOTHERMAL_WALL_FIXED:
      case BC_POROUS_WALL_MOVING:
      case BC_POROUS_WALL_FIXED:
        localff = new FluxFcnJwlWallEuler3D(iod, vfjwl, typeJac);
        break;

      case BC_INTERNAL:
        if (ns_flux == SchemeData::ROE &&
            iod.ts.implicit.ffjacobian == ImplicitData::APPROXIMATE){
          localff = new FluxFcnJwlApprJacRoeEuler3D(rshift, gamma, iod, vfjwl, typeJac);
        }else{
          fprintf(stderr, "*** Error: only the Roe flux is available for JWL\n");
          exit(-1);
        }
        break;

    }
    
  } // end - Euler or Navier-Stokes for JWL EOS

  return localff;

}

//------------------------------------------------------------------------------

inline
FluxFcnBase *FluxFcn::createFluxFcnSeg1(int rshift, int ffType, FluidModelData &fmodel, IoData &iod, 
                                        FluxFcnBase::Type typeJac){

  FluxFcnBase *localff;
  double gamma = iod.schemes.ns.gamma;

  if(fmodel.fluid == FluidModelData::PERFECT_GAS || fmodel.fluid == FluidModelData::STIFFENED_GAS){
// Euler or Navier-Stokes for Stiffened Gas 
    VarFcnSGEuler *vfsgeuler = new VarFcnSGEuler(fmodel);

    switch(ffType){

      case BC_OUTLET_FIXED:
      case BC_OUTLET_MOVING:
        if (iod.bc.outlet.type == BcsFreeStreamData::INTERNAL)
          localff = new FluxFcnSGInternalOutflowEuler3D(iod, vfsgeuler, typeJac);
        else if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL &&
                iod.schemes.bc.type == BoundarySchemeData::STEGER_WARMING)
          localff = new FluxFcnSGOutflowEuler3D(iod, vfsgeuler, typeJac);
        else if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL &&
                iod.schemes.bc.type == BoundarySchemeData::GHIDAGLIA)
          localff = new FluxFcnSGGhidagliaEuler3D(iod, vfsgeuler, typeJac);
        else if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL &&
                iod.schemes.bc.type == BoundarySchemeData::MODIFIED_GHIDAGLIA)
          localff = new FluxFcnSGModifiedGhidagliaEuler3D(iod, vfsgeuler, typeJac);
        else{
          fprintf(stderr, "*** Error: no outlet boundary flux has been selected for Stiffened Gas\n");
          exit(-1);
        }
        break;
     
      case BC_INLET_FIXED:
      case BC_INLET_MOVING:
        if (iod.bc.inlet.type == BcsFreeStreamData::INTERNAL)
          localff = new FluxFcnSGInternalInflowEuler3D(iod, vfsgeuler, typeJac);
        else if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL &&
                iod.schemes.bc.type == BoundarySchemeData::STEGER_WARMING)
          localff = new FluxFcnSGInflowEuler3D(iod, vfsgeuler, typeJac);
        else if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL &&
                iod.schemes.bc.type == BoundarySchemeData::GHIDAGLIA)
          localff = new FluxFcnSGGhidagliaEuler3D(iod, vfsgeuler, typeJac);
        else if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL &&
                iod.schemes.bc.type == BoundarySchemeData::MODIFIED_GHIDAGLIA)
          localff = new FluxFcnSGModifiedGhidagliaEuler3D(iod, vfsgeuler, typeJac);
        else{
          fprintf(stderr, "*** Error: no inlet boundary flux has been selected for Stiffened Gas\n");
          exit(-1);
        }
        break;

      case BC_ADIABATIC_WALL_MOVING:
      case BC_ADIABATIC_WALL_FIXED:
      case BC_SLIP_WALL_MOVING:
      case BC_SLIP_WALL_FIXED:
      case BC_SYMMETRY:
      case BC_ISOTHERMAL_WALL_MOVING:
      case BC_ISOTHERMAL_WALL_FIXED:
        localff = new FluxFcnSGWallEuler3D(iod, vfsgeuler, typeJac);
        break;

      case BC_POROUS_WALL_MOVING:
      case BC_POROUS_WALL_FIXED:
        localff = new FluxFcnSGPorousWallEuler3D(iod, vfsgeuler, typeJac);
        break;

      case BC_INTERNAL:
        localff = new FluxFcnSGApprJacRoeEuler3D(0, gamma, iod, vfsgeuler, typeJac);
        break;
    }
  } // end - All Stiffened Gas
  else if(fmodel.fluid == FluidModelData::LIQUID){
    // Euler or Navier-Stokes for Tait
    VarFcnTait *vftait = new VarFcnTait(fmodel);

    switch(ffType){

    case BC_OUTLET_FIXED:
    case BC_OUTLET_MOVING:
      if (iod.bc.outlet.type == BcsFreeStreamData::INTERNAL)
        localff = new FluxFcnTaitInternalOutflowEuler3D(iod, vftait, typeJac);
      else
        localff = new FluxFcnTaitGhidagliaEuler3D(iod, vftait, typeJac);
      break;

    case BC_INLET_FIXED:
    case BC_INLET_MOVING:
      if (iod.bc.inlet.type == BcsFreeStreamData::INTERNAL)
        localff = new FluxFcnTaitInternalInflowEuler3D(iod, vftait, typeJac);
      else
        localff = new FluxFcnTaitGhidagliaEuler3D(iod, vftait, typeJac);
      break;

    case BC_ADIABATIC_WALL_MOVING:
    case BC_ADIABATIC_WALL_FIXED:
    case BC_SLIP_WALL_MOVING:
    case BC_SLIP_WALL_FIXED:
    case BC_SYMMETRY:
    case BC_ISOTHERMAL_WALL_MOVING:
    case BC_ISOTHERMAL_WALL_FIXED:
    case BC_POROUS_WALL_MOVING:
    case BC_POROUS_WALL_FIXED:
      localff = new FluxFcnTaitWallEuler3D(iod, vftait, typeJac);
      break;

    case BC_INTERNAL:
      localff = new FluxFcnTaitApprJacRoeEuler3D(0, gamma, iod, vftait, typeJac);
      break;
    }
  } // end - All Tait
  else {
    fprintf(stderr, "Exiting: No turbulence model for JWL, or multiphase simulations\n");
    exit(-1);
  }

  return localff;

}

//------------------------------------------------------------------------------

inline
FluxFcnBase *FluxFcn::createFluxFcnSeg2(int rshift, int ffType, FluidModelData &fmodel, IoData &iod, 
                                        VarFcnBase *vfb, FluxFcnBase::Type typeJac){

  FluxFcnBase *localff;
  double gamma = iod.schemes.ns.gamma;

  if(fmodel.fluid == FluidModelData::PERFECT_GAS || fmodel.fluid == FluidModelData::STIFFENED_GAS){
    if (iod.eqs.type == EquationsData::NAVIER_STOKES &&
        iod.eqs.tc.type == TurbulenceClosureData::EDDY_VISCOSITY) {
      if (iod.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
          iod.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES) {
        VarFcnSGSA *vfsgsa = dynamic_cast<VarFcnSGSA *>(vfb);
        if(vfsgsa == 0){
          fprintf(stderr, "*** Error: a VarFcnSGSA is expected to create the associated FluxFcn\n");
          exit(-1);
        }
// Spalart-Allmaras for Stiffened Gas
        switch(ffType){

          case BC_DIRECTSTATE_OUTLET_FIXED:
          case BC_DIRECTSTATE_OUTLET_MOVING:
            localff = new FluxFcnSGDirectStateOutflowSAturb3D(iod, vfsgsa, typeJac);
            break;

          case BC_DIRECTSTATE_INLET_FIXED:
          case BC_DIRECTSTATE_INLET_MOVING:
            localff = new FluxFcnSGDirectStateInflowSAturb3D(iod, vfsgsa, typeJac);
            break;

          case BC_MASSFLOW_OUTLET_FIXED:
          case BC_MASSFLOW_OUTLET_MOVING:
            localff = new FluxFcnSGMassFlowOutflowSAturb3D(iod, vfsgsa, typeJac);
            break;

          case BC_MASSFLOW_INLET_FIXED:
          case BC_MASSFLOW_INLET_MOVING:
            localff = new FluxFcnSGMassFlowInflowSAturb3D(iod, vfsgsa, typeJac);
            break;

          case BC_OUTLET_FIXED:
          case BC_OUTLET_MOVING:
            if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL &&
              iod.schemes.bc.type == BoundarySchemeData::STEGER_WARMING)
              localff = new FluxFcnSGOutflowSAturb3D(iod, vfsgsa, typeJac);
            else if(iod.bc.outlet.type == BcsFreeStreamData::EXTERNAL &&
              iod.schemes.bc.type == BoundarySchemeData::GHIDAGLIA)
              localff = new FluxFcnSGGhidagliaSAturb3D(iod, vfsgsa, typeJac);
            else
              localff = new FluxFcnSGInternalOutflowSAturb3D(iod, vfsgsa, typeJac);
            break;
      
          case BC_INLET_FIXED:
          case BC_INLET_MOVING:
            if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL &&
              iod.schemes.bc.type == BoundarySchemeData::STEGER_WARMING)
              localff = new FluxFcnSGOutflowSAturb3D(iod, vfsgsa, typeJac);
            else if(iod.bc.inlet.type == BcsFreeStreamData::EXTERNAL &&
              iod.schemes.bc.type == BoundarySchemeData::GHIDAGLIA)
              localff = new FluxFcnSGGhidagliaSAturb3D(iod, vfsgsa, typeJac);
            else
              localff = new FluxFcnSGInternalInflowSAturb3D(iod, vfsgsa, typeJac);
            break;

          case BC_ADIABATIC_WALL_MOVING:
          case BC_ADIABATIC_WALL_FIXED:
          case BC_SLIP_WALL_MOVING:
          case BC_SLIP_WALL_FIXED:
          case BC_SYMMETRY:
          case BC_ISOTHERMAL_WALL_MOVING:
          case BC_ISOTHERMAL_WALL_FIXED:
            localff = new FluxFcnSGWallSAturb3D(iod, vfsgsa, typeJac);
            break;

          case BC_POROUS_WALL_MOVING:
          case BC_POROUS_WALL_FIXED:
            localff = new FluxFcnSGPorousWallSAturb3D(iod, vfsgsa, typeJac);
            break;

          case BC_INTERNAL:
            localff = new FluxFcnSGRoeSAturb3D(gamma, iod, vfsgsa, typeJac);
            break;
        }
      }
      else if (iod.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE) {
// k-epsilon turbulent model for Stiffened Gas
        VarFcnSGKE *vfsgke = dynamic_cast<VarFcnSGKE *>(vfb);
        if(vfsgke == 0){
          fprintf(stderr, "*** Error: a VarFcnSGKE is expected to create the associated FluxFcn\n");
          exit(-1);
        }
        switch(ffType){

          case BC_DIRECTSTATE_OUTLET_FIXED:
          case BC_DIRECTSTATE_OUTLET_MOVING:
            localff = new FluxFcnSGDirectStateOutflowKEturb3D(iod, vfsgke, typeJac);
            break;

          case BC_DIRECTSTATE_INLET_FIXED:
          case BC_DIRECTSTATE_INLET_MOVING:
            localff = new FluxFcnSGDirectStateInflowKEturb3D(iod, vfsgke, typeJac);
            break;

          case BC_MASSFLOW_OUTLET_FIXED:
          case BC_MASSFLOW_OUTLET_MOVING:
            localff = new FluxFcnSGMassFlowOutflowKEturb3D(iod, vfsgke, typeJac);
            break;

          case BC_MASSFLOW_INLET_FIXED:
          case BC_MASSFLOW_INLET_MOVING:
            localff = new FluxFcnSGMassFlowInflowKEturb3D(iod, vfsgke, typeJac);
            break;

          case BC_OUTLET_FIXED:
          case BC_OUTLET_MOVING:
          case BC_INLET_FIXED:
          case BC_INLET_MOVING:
            if (iod.schemes.bc.type == BoundarySchemeData::STEGER_WARMING)
              localff = new FluxFcnSGOutflowKEturb3D(iod, vfsgke, typeJac);
	    else
	      localff = new FluxFcnSGGhidagliaKEturb3D(iod, vfsgke, typeJac);
            break;
  
          case BC_ADIABATIC_WALL_MOVING:
          case BC_ADIABATIC_WALL_FIXED:
          case BC_SLIP_WALL_MOVING:
          case BC_SLIP_WALL_FIXED:
          case BC_SYMMETRY:
          case BC_ISOTHERMAL_WALL_MOVING:
          case BC_ISOTHERMAL_WALL_FIXED:
            localff = new FluxFcnSGWallKEturb3D(iod, vfsgke, typeJac);
            break;
  
          case BC_POROUS_WALL_MOVING:
          case BC_POROUS_WALL_FIXED:
            localff = new FluxFcnSGPorousWallKEturb3D(iod, vfsgke, typeJac);
            break;

          case BC_INTERNAL:
            localff = new FluxFcnSGRoeKEturb3D(gamma, iod, vfsgke, typeJac);
            break;

        }
      } // end - k-epsilon turbulent model for Stiffened Gas
      else {
        fprintf(stderr,"Error: Seg. solver is only implemented for SA and KE models.\n");
        exit(-1);
      }
    } // end - turbulence
    else {
      fprintf(stderr,"Error: Seg. solver is only implemented for Navier-Stokes Eqs.\n");
      exit(-1);
    }
  } // end - Stiffened Gas
  else if(fmodel.fluid == FluidModelData::LIQUID){
    if (iod.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_SPALART_ALLMARAS ||
        iod.eqs.tc.tm.type == TurbulenceModelData::ONE_EQUATION_DES) {
      VarFcnTaitSA *vftaitsa = dynamic_cast<VarFcnTaitSA *>(vfb);
      if(vftaitsa == 0){
        fprintf(stderr, "Error: a VarFcnTaitSA is expected to create the associated FluxFcn.\n");
	exit(-1);
      }

      // Spalart-Allmaras for Tait
      switch(ffType){

      case BC_DIRECTSTATE_OUTLET_FIXED:
      case BC_DIRECTSTATE_OUTLET_MOVING:
      case BC_DIRECTSTATE_INLET_FIXED:
      case BC_DIRECTSTATE_INLET_MOVING:
      case BC_MASSFLOW_OUTLET_FIXED:
      case BC_MASSFLOW_OUTLET_MOVING:
      case BC_MASSFLOW_INLET_FIXED:
      case BC_MASSFLOW_INLET_MOVING:
        break;

      case BC_OUTLET_FIXED:
      case BC_OUTLET_MOVING:
      case BC_INLET_FIXED:
      case BC_INLET_MOVING:
        localff = new FluxFcnTaitGhidagliaSAturb3D(iod, vftaitsa, typeJac);
        break;

      case BC_ADIABATIC_WALL_MOVING:
      case BC_ADIABATIC_WALL_FIXED:
      case BC_SLIP_WALL_MOVING:
      case BC_SLIP_WALL_FIXED:
      case BC_SYMMETRY:
      case BC_ISOTHERMAL_WALL_MOVING:
      case BC_ISOTHERMAL_WALL_FIXED:
      case BC_POROUS_WALL_MOVING:
      case BC_POROUS_WALL_FIXED:
        localff = new FluxFcnTaitWallSAturb3D(iod, vftaitsa, typeJac);
        break;

      case BC_INTERNAL:
        localff = new FluxFcnTaitRoeSAturb3D(gamma, iod, vftaitsa, typeJac);
        break;
      }
    }
    else if (iod.eqs.tc.tm.type == TurbulenceModelData::TWO_EQUATION_KE) {
      // k-epsilon turbulent model for Tait
      VarFcnTaitKE *vftaitke = dynamic_cast<VarFcnTaitKE *>(vfb);
      if(vftaitke == 0){
	fprintf(stderr, "Error: a VarFcnTaitKE is expected to create the associated FluxFcn.\n");
	exit(-1);
      }

      switch(ffType){

      case BC_DIRECTSTATE_OUTLET_FIXED:
      case BC_DIRECTSTATE_OUTLET_MOVING:
      case BC_DIRECTSTATE_INLET_FIXED:
      case BC_DIRECTSTATE_INLET_MOVING:
      case BC_MASSFLOW_OUTLET_FIXED:
      case BC_MASSFLOW_OUTLET_MOVING:
      case BC_MASSFLOW_INLET_FIXED:
      case BC_MASSFLOW_INLET_MOVING:
        break;

      case BC_OUTLET_FIXED:
      case BC_OUTLET_MOVING:
      case BC_INLET_FIXED:
      case BC_INLET_MOVING:
        localff = new FluxFcnTaitGhidagliaKEturb3D(iod, vftaitke, typeJac);
        break;

      case BC_ADIABATIC_WALL_MOVING:
      case BC_ADIABATIC_WALL_FIXED:
      case BC_SLIP_WALL_MOVING:
      case BC_SLIP_WALL_FIXED:
      case BC_SYMMETRY:
      case BC_ISOTHERMAL_WALL_MOVING:
      case BC_ISOTHERMAL_WALL_FIXED:
      case BC_POROUS_WALL_MOVING:
      case BC_POROUS_WALL_FIXED:
        localff = new FluxFcnTaitWallKEturb3D(iod, vftaitke, typeJac);
        break;

      case BC_INTERNAL:
        localff = new FluxFcnTaitRoeKEturb3D(gamma, iod, vftaitke, typeJac);
        break;

      }
    } // end - k-epsilon turbulent model for Tait
  } else {
    fprintf(stderr, "Exiting: No turbulence model for Tait, JWL, or multiphase simulations\n");
    exit(-1);
  }
 

  return localff;

}
//------------------------------------------------------------------------------

#endif
