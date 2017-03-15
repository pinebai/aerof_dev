#ifndef _DIST_BC_DATA_H_
#define _DIST_BC_DATA_H_

#include <DistVector.h>
#include <VarFcn.h>
#include <./Dev/devtools.h>//TODO delete line

class IoData;
class SubDomain;
class Communicator;
class Domain;

template<int dim> class BcData;

//------------------------------------------------------------------------------
/*
   angles[0] in the x-z plane (flow parallel to x-axis if alpha[0]=0)
   angles[1] in the x-y plane (flow parallel to x-axis if alpha[1]=0)
*/

/** Distributed Boundary Condition data.
 *
 */
template<int dim>
class DistBcData {

protected:
  
  DistInfo& nodeDistInfo,&inletNodeDistInfo,&faceDistInfo;
   
  enum BoundaryFluid { GAS=0, TAIT=1, JWL=2 } boundaryFluid;
  /** angles of the incoming flow
     *
     * - angles[0] is \f$\alpha\f$ in the x-z plane (flow parallel to x-axis if angles[0]=0).
     * - angles[1] is \f$\beta\f$ in the x'-y plane (flow parallel to x'-axis if angles[1]=0).
    */
  double angles[2];

  bool gravityOn; //!< Whether gravity is activated or not.
  double gravity; //!< intensity of the gravity
  double depth;
  double ngravity[3]; //!< direction of the gravity

  VarFcn *vf;

  DistSVec<double,3> Xdot;
  DistVec<double> Temp;

  DistSVec<double,dim> Ufarin;
  DistSVec<double,dim> Ufarout;
  DistSVec<double,dim> Uporouswall;
  double Uin[dim], Vin[dim];
  double Uout[dim], Vout[dim];
  double Ub[dim];

  DistSVec<double,dim> Uface;
  DistSVec<double,dim> Unode;
  DistSVec<double,dim> Uinletnode;

  int numLocSub;

  Communicator *com;
  SubDomain **subDomain;
  BcData<dim> **subBcData;

  map<int, RotationData *> &rotInfo;
  double tref;
  double vref;

// Included (MB)
  //for steady sensitivity analysis
  DistSVec<double,3> *dXdot;
  DistVec<double> *dTemp;
  DistSVec<double,dim> *dUface;
  DistSVec<double,dim> *dUnode;
  DistSVec<double,dim> *dUinletnode;
  DistSVec<double,dim> *dUfarin;
  DistSVec<double,dim> *dUfarout;
  DistSVec<double,dim> *dUporouswall;
  double dUin[dim], dVin[dim];
  double dUout[dim], dVout[dim];

  // only used for sensitivity analsysis ???
  DistSVec<double,dim> *dUfaceSA;
  DistSVec<double,dim> *dUnodeSA;

  double dtrefdMach;
  double dvrefdMach;
 
  DistVec<double> *boundaryStateHH; 

protected:

  void finalize(DistSVec<double,3> &);

// Included (MB)
  void finalizeSA(DistSVec<double,3> &, DistSVec<double,3> &, double &);

public:

  DistBcData(IoData &, VarFcn *, Domain *,DistInfo* nodeDistInfo = 0,
              DistInfo* inletNodeDistInfo = 0, DistInfo* faceDistInfo = 0);
  virtual ~DistBcData();

  BcData<dim> &operator() (int i) const { return *subBcData[i]; }

  void update(DistSVec<double,3> &);
  virtual void updateFarField(DistSVec<double,3> &) {}
  virtual void computeNodeValue(DistSVec<double,3> &) {}

  DistSVec<double,3> &getVelocityVector() { return Xdot; }
  DistVec<double> &getTemperatureVector() { return Temp; }
  DistSVec<double,dim> &getInletBoundaryVector()  { return Ufarin;  }
  DistSVec<double,dim> &getOutletBoundaryVector() { return Ufarout; }

  DistSVec<double,dim> &getFaceStateVector() { return Uface; }
  DistSVec<double,dim> &getNodeStateVector() { return Unode; }
  DistSVec<double,dim> &getUfarin() { return Ufarin; }
  DistSVec<double,dim> &getUfarout() { return Ufarout; }
  DistSVec<double,dim> &getUporouswall() { return Uporouswall; }

  DistVec<double>* getBoundaryStateHH() { return boundaryStateHH; }

  double *getInletAngles() { return angles; }
  double *getInletConservativeState() { return Uin; }
  double *getOutletConservativeState() { return Uout; }
  double *getInterface() { return Ub; }
  double *getInletPrimitiveState() { return Vin; }

// Included (MB)
  void updateSA(DistSVec<double,3> &, DistSVec<double,3> &, double &);
  void rstVar(IoData &ioData);

  virtual void computeDerivativeOfNodeValue(DistSVec<double,3> &, DistSVec<double,3> &)
  {
    std::cout<<"ERROR: header computeDerivativeOfNodeValue called"<<std::endl; exit(-1);
  }
  virtual void computeNodeWallValues(DistSVec<double,3> &)
  {std::cout<<"ERROR: header computeNodeWallValues called"<<std::endl; exit(-1);}//TODO BUUGHUNT

  DistVec<double> &getDerivativeOfTemperatureVector() { return (*dTemp); }

  virtual void updateFarFieldSA(DistSVec<double,3> &, DistSVec<double,3> &, double &)
  {std::cout<<"ERROR: header updateFarFieldSA called"<<std::endl; exit(-1);}//TODO BUUGHUNT

  virtual void initialize(IoData &, DistSVec<double,3> &)
  {std::cout<<"ERROR: header initialize called"<<std::endl; exit(-1);}//TODO BUUGHUNT

  virtual void initializeSA(IoData &, DistSVec<double,3> &, DistSVec<double,3> &, double &, double &, double &)
  {std::cout<<"ERROR: header initializeSA called"<<std::endl; exit(-1);}//TODO BUUGHUNT

  DistSVec<double,3> &getDerivativeOfVelocityVector() { return *dXdot; }

};

//------------------------------------------------------------------------------

template<int dim>
class DistBcDataEuler : public DistBcData<dim> {

private:
  void setBoundaryConditionsGas(IoData &, DistSVec<double,3> &);
  void setBoundaryConditionsLiquid(IoData &, DistSVec<double,3> &);
  void setBoundaryConditionsJWL(IoData &, DistSVec<double,3> &);
  void setBoundaryConditionsGasGas(IoData &, DistSVec<double,3> &);
  void setBoundaryConditionsLiquidLiquid(IoData &, DistSVec<double,3> &);
  void setBoundaryConditionsGasLiquid(IoData &, DistSVec<double,3> &);
  void setBoundaryConditionsLiquidGas(IoData &, DistSVec<double,3> &);
  void setBoundaryConditionsJWLGas(IoData &, DistSVec<double,3> &);
  void setBoundaryConditionsJWLLiquid(IoData &iod,
				      DistSVec<double,3> &X);

  void updateFarField(DistSVec<double,3> &);
  void updateFarFieldGas(DistSVec<double,3> &);
  void updateFarFieldLiquid(DistSVec<double,3> &);
  void updateFarFieldJWL(DistSVec<double,3> &);

// Included (MB)
  void initialize(IoData &, DistSVec<double,3> &);
  void initializeSA(
	     IoData &,  //
		 DistSVec<double,3> &, //X->nodal position
		 DistSVec<double,3> &, //dX->derivative of the nodal position
		 double &,  //dM->marks the derivative od mach number;          set at DFSPAR[0] in FluidShapeOptimizationHandler.C
		 double &,  //dA->marks the derivatice of angle of attck (AoA); set at DFSPAR[1] in FluidShapeOptimizationHandler.C
		 double &); //dB->marks the derivative of sidelipe angle;       set at DFSPAR[2] in FluidShapeOptimizationHandler.C

  void setDerivativeOfBoundaryConditionsGas(IoData &,DistSVec<double,3> &, DistSVec<double,3> &, double, double, double);
  void updateFarFieldSA(DistSVec<double,3> &, DistSVec<double,3> &, double &);
  void updateFarFieldGasSA(DistSVec<double,3> &, DistSVec<double,3> &, double &);

public:

  void computeDerivativeOfNodeValue(DistSVec<double,3> &, DistSVec<double,3> &);
  DistBcDataEuler(IoData &, VarFcn *, Domain *, DistSVec<double,3> &);
  ~DistBcDataEuler() {}

};

//------------------------------------------------------------------------------

template<int dim>
class DistBcDataSA : public DistBcDataEuler<dim> {

  DistSVec<double,2> *tmp;
  CommPattern<double> *vec2Pat;

// Included (MB)
  DistSVec<double,2> *dtmp;
  DistSVec<double,1> *dnormsa;
  DistSVec<double,dim> *dtmpsa;
  CommPattern<double> *vec1PatSA;
  CommPattern<double> *vec2PatSA;

public:

  DistBcDataSA(IoData &, VarFcn *, Domain *, DistSVec<double,3> &);
  ~DistBcDataSA();

  void computeNodeValue(DistSVec<double,3> &);

// Included (MB)
  void computeDerivativeOfNodeValue(DistSVec<double,3> &, DistSVec<double,3> &);
  void computeNodeWallValues(DistSVec<double,3> &);

};

//------------------------------------------------------------------------------

template<int dim>
class DistBcDataKE : public DistBcDataEuler<dim> {

  DistSVec<double,3> *tmp;
  CommPattern<double> *vec3Pat;

// Included (MB)
  DistSVec<double,3> *dtmp;

public:

  DistBcDataKE(IoData &, VarFcn *, Domain *, DistSVec<double,3> &);
  ~DistBcDataKE();

  void computeNodeValue(DistSVec<double,3> &);

// Included (MB)
  void computeDerivativeOfNodeValue(DistSVec<double,3> &, DistSVec<double,3> &);
  void computeNodeWallValues(DistSVec<double,3> &x)
  {std::cout<<"ERROR: header DistBcDataKE::computeNodeWallValues called"<<std::endl; exit(-1);}//TODO BUUGHUNT

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <DistBcData.C>
#endif

#endif
