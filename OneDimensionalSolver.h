#ifndef _ONE_DIMENSIONAL_SOLVER_H_
#define _ONE_DIMENSIONAL_SOLVER_H_


//------------------------------------------------------------------------------
// Created March 2010 by A. Rallu @Stanford University
//------------------------------------------------------------------------------
// This class is used to solve the one-dimensional Euler equations.
// They can be written in cartesian, cylindrical or spherical coordinates.
// All these cases are considered in this class.
// This class is geared toward solving two-phase one-dimensional Euler equations
// with a possible burn process of explosives.
//
// In order to allow some flexibility and ease of use, this class
// uses only some of the local-level classes of the rest of the AERO-F code
// (like SVec, FluxFcn, LocalRiemann, VarFcn). Hence, it does not use
// classes like DistTimeState or SpaceOperator. One reason is that these
// classes then always refer to Edges and Faces and Elems which are not
// available for the simulations intended with this class. (It would be possible
// to go through Edges, Faces, Elems and therefore use SpaceOperator and 
// DistTimeState and others, but that requires more work/time than I can 
// afford now).
//
//
// Note that state vectors have 5 coordinates. Only three of them are used.
// The coordinates [2] and [3] are unused but must be present in order to 
// be able to use the functions of VarFcn, FluxFcn, LocalRiemann (as are).
//
// Note: the spherical one-D Euler equations can be written as
// d(r^2 U)/dt + d(r^2 F(U))/dr = S2(U,r) = {0, 2*r*pressure, 0}
// or
// dU/dt + dF(U)/dx = S(U,r) = {-a*density*velocity/r, -a*density*velocity^2/r, -a*(density*energy+pressure)*velocity/r}
// where a = 2 for spherical coordinates (if a = 1, these are the cylindrical
// one-dimensional Euler equations)
// The first should require that the conserved quantities are r^2*U.
// The second requires that the source be integrated in a specific manner at r=0
//
//------------------------------------------------------------------------------
#include "SubDomain.h"
#include "Vector.h"
#include "FluidSelector.h"
#include "IoData.h"
#include "RefVal.h"
#include "RecFcn.h"

#include "RKIntegrator.h"
#include "ExactRiemannSolver.h"
#include "ProgrammedBurn.h"

#include "PostFcn.h"

#include <fstream>

class FluxFcn;
class VarFcn;
class LocalRiemann;
class OneDimensionalSourceTerm;
class Domain;

//------------------------------------------------------------------------------

class OneDimensional {
  const static int dim = 5, dimLS = 1;
  OneDimensionalInfo::CoordinateType coordType;
  OneDimensionalInfo::VolumeType volumeType;

  enum ProblemMode { MultiFluid, FSI};

  ProblemMode problemMode;

  int numPoints;
  double maxDistance;
  SVec<double,1> ctrlVol;
  // for cartesian, cylindrical and spherical one-D simulation with control surfaces
  SVec<double,1> ctrlSurf;
  SVec<double,1> X; // center of control volumes
  SVec<double,1> Y; // control surface locations

  double finalTime;
  double cfl;

  double BC[2][5];
  double BCphi[2];
  
  SVec<double,5> U;
  SVec<double,5> V;
  SVec<double,5> R;
  SVec<double,5> gradV;

  SVec<double,1> Phi;
  SVec<double,1> Phin;
  SVec<double,1> Rphi;
  SVec<double,1> gradPhi;
  Vec<int> fluidId;
  Vec<int> fluidIdn;

  FluxFcn **fluxFcn;
  VarFcn *varFcn;
  ExactRiemannSolver<5>* riemann;
  FluidSelector fluidSelector;
  // for cartesian, cylindrical and spherical one-D simulation with source term
  OneDimensionalSourceTerm *source; 

  int frequency; //postprocessing output frequency
  char *outfile;
  RefVal refVal;

  // Riemann solutions at the interface
  // Needed for appropriate handling of phase change updates.
  SVec<double,5> Wr;

  SVec<double,5> lastPhaseChangeValue;
  
  SVec<double,5> Vslope;
  SVec<double,1> Phislope;

  Vec<int> riemannStatus;
  Vec<int> cutCellStatus;

  RecFcn* recFcn, *recFcnLS;

  char* bubbleRadiusFile;

  RecFcn* createRecFcn(IoData &ioData);
  RecFcn* createRecFcnLS(IoData &ioData);

  RKIntegrator< SVec<double, 5> >* Vintegrator;
  RKIntegrator< SVec<double, 1> >* Phiintegrator;

  double time;

  SVec<double,dim> rupdate;
  Vec<double> weight;
  Vec<int> fidToSet;
  SVec<double,dim-2> interfacialWi;
  SVec<double,dim-2> interfacialWj;

  void loadSparseGrid(IoData&);

  SparseGridCluster* tabulationC;

  ProgrammedBurn* programmedBurn;

  double programmedBurnStopPercentDistance;
  bool programmedBurnIsUsed;

  double sscale[PostFcn::SSIZE];
  double vscale[PostFcn::SSIZE];

  char *solutions;
  char *scalars[PostFcn::SSIZE];
  char *vectors[PostFcn::VSIZE];

  struct {

    int numNodes;
    int step;
    std::vector<Vec3D> locations;
    std::vector<int> ids;
    std::vector<double> alpha;
  } nodal_output;
  
  char *nodal_scalars[PostFcn::SSIZE];
  char *nodal_vectors[PostFcn::VSIZE];

  void setupOutputFiles(IoData& iod);
  void setupFixes(IoData& iod);

  int* loctag;

  int typePhaseChange;

  bool isSixthOrder;

  bool isSinglePhase;

  double beta;
  double interfaceLocation;

  int interfaceTreatment;
  int interfaceExtrapolation;

  // 0 - constant reconstruction
  // 1 - linear reconstruction
  // 2 - limited
  int limiterLeft,limiterRight;

  int levelSetMethod;

  Timer* myTimer;

  void setupProbes(IoData& ioData);
  void outputProbes(double,int);

  class Veval {

  public:
    int findNode(const double* loc,double& localRadius) {

      if (spherical)
	localRadius = sqrt((loc[0]-bubble_x0)*(loc[0]-bubble_x0)+(loc[1]-bubble_y0)*(loc[1]-bubble_y0)+(loc[2]-bubble_z0)*(loc[2]-bubble_z0));
      else
	localRadius = loc[0];
	  
      // If the node is inside the sphere, set its values accordingly
      if (localRadius < max_distance) {
	
	// Do a binary search to find the closest point whose r
	// is less than or equal to localRadius
	int a = numPoints/2;
	int rmin = 0,rmax=numPoints-2;
	while (!(x_1D[a] <= localRadius && x_1D[a+1] > localRadius)  ) {
	  if (x_1D[a] < localRadius)
	    rmin = a;
	  else
	    rmax = a;
	  
	  if (a == rmin)
	    a = (rmax+rmin)/2 + 1;
	  else
	    a = (rmax+rmin)/2 ;
	}
	
	return a;
      } else
	return -1;
    }

    Veval(VarFcn* _vf,OneDimensionalInputData* _oned,FluidSelector* _fluidSelector,double* _x, double* _v,int* _fids,SVec<double,3>* _X,int _np,
	  double x0,double y0, double z0, bool sph = true)  : x_1D(_x), v_1D(_v), fids(_fids),X(_X), numPoints(_np),
      bubble_x0(x0), bubble_y0(y0), bubble_z0(z0),/* fluidSelector(_fluidSelector),*/ oned(_oned),varFcn(_vf) , spherical(sph)
      { 
      
      int i,fid_new = 0;
      memset(outletState,0,sizeof(double)*5);
      for (i = 0; i < numPoints; ++i) {
	if (fids[i] == 0) {

	  memset(boundaryStateL,0,sizeof(double)*5);
	  memset(boundaryStateR,0,sizeof(double)*5);
	  
	  boundaryStateL[0] = v_1D[(i-1)*5];
	  boundaryStateR[0] = v_1D[i*5];
	  boundaryStateL[1] = v_1D[(i-1)*5+1];
	  boundaryStateR[1] = v_1D[i*5+1];

	  //if (fluidSelector) {
	    fid_new = fids[i-1];
	    if (oned && oned->fluidRemap.dataMap.find(fids[i-1]) != oned->fluidRemap.dataMap.end())
	      fid_new = oned->fluidRemap.dataMap.find(fids[i-1])->second->newID;
	    //}
	  if (varFcn->getType(fid_new) == VarFcnBase::TAIT)
	    boundaryStateL[4] = v_1D[(i-1)*5+4];
	  else
	    boundaryStateL[4] = v_1D[(i-1)*5+2];

	  //if (fluidSelector) {
	    fid_new = fids[i];
	    if (oned && oned->fluidRemap.dataMap.find(fids[i]) != oned->fluidRemap.dataMap.end())
	      fid_new = oned->fluidRemap.dataMap.find(fids[i])->second->newID;
	    //}
	  if (varFcn->getType(fid_new) == VarFcnBase::TAIT)
	    boundaryStateR[4] = v_1D[i*5+4];
	  else
	    boundaryStateR[4] = v_1D[i*5+2];
	  
	  rad = (x_1D[i-1]*v_1D[i*5+3]-x_1D[i]*v_1D[(i-1)*5+3])/(v_1D[i*5+3]-v_1D[(i-1)*5+3]);
	  break;
	}
      }
      max_distance = x_1D[numPoints-1];
      outletState[0] = v_1D[(numPoints-1)*5];
      //if (fluidSelector) {
	fid_new = fids[numPoints-1];
	if (oned && oned->fluidRemap.dataMap.find(fids[numPoints-1]) != oned->fluidRemap.dataMap.end())
	  fid_new = oned->fluidRemap.dataMap.find(fids[numPoints-1])->second->newID;
	//}
      if (varFcn->getType(fid_new) == VarFcnBase::TAIT)
	outletState[4] = v_1D[(numPoints-1)*5+4];
      else
	outletState[4] = v_1D[(numPoints-1)*5+2];
    }
      
    ~Veval() { }

    void Eval(int node, const double* loc,double* f) {

      double localRadius;
      double ff[5];
      int lsdim;
      int i = node;
      int fid_new=0;
      double bub_x0[3] = { bubble_x0, bubble_y0, bubble_z0 };
      double xrad;
      if (spherical)
	xrad = sqrt(((*X)[i][0]-bubble_x0)*((*X)[i][0]-bubble_x0)+((*X)[i][1]-bubble_y0)*((*X)[i][1]-bubble_y0)+((*X)[i][2]-bubble_z0)*((*X)[i][2]-bubble_z0));	  
      else
	xrad = (*X)[i][0];

      int a = findNode(loc,localRadius);
      //if (xrad < 1.0)
      //	std::cout << xrad << " " << a<<std::endl;
      if (a >= 0) {

	double alpha = (localRadius-x_1D[a])/(x_1D[a+1]-x_1D[a]);
	if (fids[a] != fids[a+1]) {
	  if(rad-localRadius > 0.0) {
	    alpha = 0.0; 
	    fid_new = fids[a];
	    if (oned && oned->fluidRemap.dataMap.find(fids[a]) != oned->fluidRemap.dataMap.end())
	      fid_new = oned->fluidRemap.dataMap.find(fids[a])->second->newID;
          }
	  else {
	    alpha = 1.0;
	    fid_new = fids[a+1];
	    if (oned && oned->fluidRemap.dataMap.find(fids[a+1]) != oned->fluidRemap.dataMap.end())
	      fid_new = oned->fluidRemap.dataMap.find(fids[a+1])->second->newID;
          }
	} else {
	  fid_new = fids[a];
	  if (oned && oned->fluidRemap.dataMap.find(fids[a]) != oned->fluidRemap.dataMap.end())
	    fid_new = oned->fluidRemap.dataMap.find(fids[a])->second->newID;

        }

	
	//if (fluidSelector) {
	  //lsdim = fluidSelector->getLevelSetDim(0,fid_new);
	  //}

	  if (!spherical)
	    fid_new = 1-fid_new;
          fid_new = std::min<int>(fid_new, varFcn->size()-1);

	if ( (localRadius > rad && xrad > rad) || (localRadius <= rad && xrad <= rad)) {
	  //if (xrad < 1.0)
	  //  std::cout << "fid = " << fid_new <<std::endl;
	  ff[0] = v_1D[a*5]*(1.0-alpha)+v_1D[(a+1)*5]*(alpha);
	  if (spherical) {
	    ff[1] = (v_1D[a*5+1]*(1.0-alpha)+v_1D[(a+1)*5+1]*(alpha))*(loc[0]-bubble_x0)/max(localRadius,1.0e-8);
	    ff[2] = (v_1D[a*5+1]*(1.0-alpha)+v_1D[(a+1)*5+1]*(alpha))*(loc[1]-bubble_y0)/max(localRadius,1.0e-8);
	    ff[3] = (v_1D[a*5+1]*(1.0-alpha)+v_1D[(a+1)*5+1]*(alpha))*(loc[2]-bubble_z0)/max(localRadius,1.0e-8);
	  } else {
	    ff[1] = (v_1D[a*5+1]*(1.0-alpha)+v_1D[(a+1)*5+1]*(alpha));
	    ff[2] = 0.0;
	    ff[3] = 0.0;
	  }
	  if (varFcn->getType(fid_new) == VarFcnBase::TAIT)
	    ff[4] = v_1D[a*5+4]*(1.0-alpha)+v_1D[(a+1)*5+4]*(alpha);
	  else
	    ff[4] = v_1D[a*5+2]*(1.0-alpha)+v_1D[(a+1)*5+2]*(alpha);
	  varFcn->primitiveToConservative(ff,f,fid_new);
	} else if (xrad <= rad) {
	  fid_new = fids[0];
	  if (oned && oned->fluidRemap.dataMap.find(fids[0]) != oned->fluidRemap.dataMap.end())
	    fid_new = oned->fluidRemap.dataMap.find(fids[0])->second->newID;
	  varFcn->primitiveToConservative(boundaryStateL,f,fid_new);
	  if (spherical) {
	    for (int j = 1; j <= 3; ++j) {
	      f[j] = boundaryStateL[0]*boundaryStateL[1]*(loc[j-1]-bub_x0[j-1])/max(localRadius,1.0e-8);
	    }
	  } else {
	    f[1] = boundaryStateL[0]*boundaryStateL[1];
	    f[2] = f[3] = 0.0;
	  }
	    
	} else {
	  fid_new = fids[numPoints-1];
	  if (oned && oned->fluidRemap.dataMap.find(fids[numPoints-1]) != oned->fluidRemap.dataMap.end())
	    fid_new = oned->fluidRemap.dataMap.find(fids[numPoints-1])->second->newID;
	  varFcn->primitiveToConservative(boundaryStateR,f,fid_new);
	  if (spherical) {
	    for (int j = 1; j <= 3; ++j) {
	      f[j] = boundaryStateR[0]*boundaryStateR[1]*(loc[j-1]-bub_x0[j-1])/max(localRadius,1.0e-8);
	    }	  
	  } else {
	    f[1] = boundaryStateR[0]*boundaryStateR[1];
	    f[2] = f[3] = 0.0;
	  }
	}	  
      }
      else {
	//if (fluidSelector) {
	  fid_new = fids[numPoints-1];
	  if (oned && oned->fluidRemap.dataMap.find(fids[numPoints-1]) != oned->fluidRemap.dataMap.end())
	    fid_new = oned->fluidRemap.dataMap.find(fids[numPoints-1])->second->newID;
	  // lsdim = fluidSelector->getLevelSetDim(0,fid_new);
	  //}
	varFcn->primitiveToConservative(outletState,f,fid_new);	
      }
      assert(f[0] > 0.0);
      //if (f[4] <= 0.0)
      //	std::cout << "Error, neg energy! " << f[4] << " " << localRadius << std::endl;
      assert(f[4] > 0.0);
      double v[5];
      varFcn->conservativeToPrimitive(f,v,fid_new);
      for (int i = 0; i < dim; ++i)
	f[i] = v[i];
      //	  fid_new = fids[numPoints-1];
      //if (localRadius < 0.5)
      //	std::cout << "1: " << fid_new << " " <<v[4]/2.0e-6*100 << std::endl;
    }

  private:

    //FluidSelector* fluidSelector;
    bool spherical;
    double* x_1D,*v_1D;
    double boundaryStateL[5],boundaryStateR[5],outletState[5];
    SVec<double,3>* X;
    int* fids;
    double bubble_x0, bubble_y0,bubble_z0;
    int numPoints;
    double rad;
    OneDimensionalInputData* oned;
    VarFcn* varFcn;
    double max_distance;
  };

 public:

    void levelSetDerivative(double t0, Vec<double>& phi, Vec<double>& k);

  void EulerF(double t, SVec<double,5>& y,SVec<double,5>& k);
  void PhiF(double t, SVec<double,1>& y,SVec<double,1>& k);

  enum ReadMode { ModeU, ModePhi };
  template <int dimp,int dimLS>
  static void read1DSolution(IoData& iod, DistSVec<double,dimp>& Up, 
			     DistSVec<double,dimLS>* Phi,
			     FluidSelector* fluidSelector,
			     VarFcn* varFcn,
			     DistSVec<double,3>& X,
			     Domain& dom,
			     ReadMode mode,
			     bool spherical = true); 

  template <int dimp,int dimLS>
  static void read1DSolution(IoData& iod, const char* filename, DistSVec<double,dimp>& Up, 
			     DistSVec<double,dimLS>* Phi,
			     FluidSelector* fluidSelector,
			     VarFcn* varFcn,
			     DistSVec<double,3>& X,
			     Domain& dom,
			     ReadMode mode,
			     bool spherical = true) {

      // read 1D solution
      DistSVec<double,dimp> ut(Up);

      std::fstream input;
      input.open(filename, fstream::in);
      //cout << filename << endl;
      if (!input.is_open()) {
        cout<<"*** Error: could not open 1D solution file "<<filename<<endl;
        exit(1);
      }
      
      input.ignore(256,'\n');
      input.ignore(2,' ');
      int numPoints = 0;
      input >> numPoints;
      //cout << " Num 1d points = " << numPoints << endl;
      double* x_1D = new double[numPoints];
      double* v_1D = new double[numPoints*5];// rho, u, p, phi, T
      int* fids = new int[numPoints];
      
      double rad = 0;
      for(int i=0; i<numPoints; i++){
        input >> x_1D[i] >> v_1D[i*5] >> v_1D[i*5+1] >> v_1D[i*5+2] >> v_1D[i*5+3]>> fids[i] >> v_1D[i*5+4];
        fids[i] = std::min<int>(fids[i], varFcn->size()-1);
        x_1D[i]    /= iod.ref.rv.length;
        v_1D[i*5] /= iod.ref.rv.density;
        v_1D[i*5+1] /= iod.ref.rv.velocity;
        v_1D[i*5+2] /= iod.ref.rv.pressure;
        //v_1D[i][3] /= iod.ref.rv.length;
        v_1D[i*5+4] /= iod.ref.rv.temperature;
        //std::cout << v_1D[i*5+2] << std::endl;
        if (rad == 0 && fids[i] == 0) {
          rad = (x_1D[i-1]*v_1D[i*5+3]-x_1D[i]*v_1D[(i-1)*5+3])/(v_1D[i*5+3]-v_1D[(i-1)*5+3]);
        }
      }
      
      input.close();

      int lsdim = 0;
      
      // interpolation assuming 1D solution is centered on bubble_coord0
      double bubble_x0 = 0;//itr->second->x0;
      double bubble_y0 = 0;//itr->second->y0;
      double bubble_z0 = 0;//itr->second->z0;
      double max_distance = x_1D[numPoints-1]; 
      double localRadius;
      ut = 1.0;
#pragma omp parallel for
      for (int iSub=0; iSub<Up.numLocSub(); ++iSub) {
        SVec<double,dimp> &u(Up(iSub));
        SVec<double, 3> &x(X(iSub));
        for(int i=0; i<u.size(); i++) {
          localRadius = sqrt((x[i][0]-bubble_x0)*(x[i][0]-bubble_x0)+(x[i][1]-bubble_y0)*(x[i][1]-bubble_y0)+(x[i][2]-bubble_z0)*(x[i][2]-bubble_z0));
          for (int k = 0; k < dim; ++k)
            ut.subData(iSub)[i][k] = 0.0;
	      }

        Veval veval(varFcn,0,fluidSelector,x_1D,v_1D,fids,&x,numPoints,bubble_x0,bubble_y0,bubble_z0,spherical);

        double localRadius; int np;
        double localAlpha, velocity_r;
        double localV[5];
        if (mode == ModeU) {
          for(int i=0; i<u.size(); i++) {
            double v[5];
            localRadius = sqrt((x[i][0]-bubble_x0)*(x[i][0]-bubble_x0)+(x[i][1]-bubble_y0)*(x[i][1]-bubble_y0)+(x[i][2]-bubble_z0)*(x[i][2]-bubble_z0));
            veval.Eval(i, x[i], u[i]);
          }
        }
      }
}

  void resultsOutput(double time, int iteration);
  void restartOutput(double time, int iteration);

  void output2DVTK();

public:
  OneDimensional(int, double*, IoData &ioData, Domain *domain);
  ~OneDimensional();

  void spatialSetup();
  void temporalSetup();
  void stateInitialization(OneDimensionalInfo &data);
  void totalTimeIntegration();
  double computeMaxTimeStep();
  void singleTimeIntegration(double dt);
  void computeEulerFluxes(SVec<double,5>& y);
  void computeLevelSetFluxes(SVec<double,1>& y);

  static void load1DMesh(IoData& ioData,int& numPts,double* &meshPoints);

  template <int dim>
    void computeSlopes(SVec<double,dim>& VV, SVec<double,dim>& slopes,
		       Vec<int>& fid, bool crossInterface);

};

//#ifdef TEMPLATE_FIX
//#include <OneDimensionalSolver.C>
//#endif

#endif
