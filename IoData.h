#ifndef _IO_DATA_H_
#define _IO_DATA_H_

#include <RefVal.h>
#include <cstdio>
#include <map>
#include "parser/ParseTree.h"
#include "parser/Dictionary.h"
//#include "AlternatingLeastSquare/als_io.h" disable temporarily for testing

using std::map;

class Assigner;
class ClassAssigner;
class Communicator;

//------------------------------------------------------------------------------

template<class DataType>
class ObjectMap {

public:

  map<int, DataType *> dataMap;
  void setup(const char *name, ClassAssigner *);
  ~ObjectMap()
    {
      for(typename map<int, DataType *>::iterator it=dataMap.begin();it!=dataMap.end();++it)
      {
        delete it->second;
      }
    }
};

//------------------------------------------------------------------------------
struct FluidRemapData {

  FluidRemapData();
  ~FluidRemapData() {}

  int oldID,newID;

  void setup(const char *, ClassAssigner * = 0);

  Assigner *getAssigner();
};

struct OneDimensionalInputData {

  OneDimensionalInputData();
  ~OneDimensionalInputData() {}
  const char* file;
  double x0,y0,z0;

  ObjectMap<FluidRemapData> fluidRemap;

  Assigner *getAssigner();

  void setup(const char *, ClassAssigner * = 0);
};


struct InputData {

  enum OptimalPressureDimensionality {NON_DIMENSIONAL=0, DIMENSIONAL=1,NONE=2} optPressureDim;
  enum ShapeDerivativeType {WALL=0, VOLUME=1} shapederivativesType;
  enum UseMultiSolutionsGappy {MULTI_SOLUTIONS_GAPPY_FALSE=0, MULTI_SOLUTIONS_GAPPY_TRUE=1} useMultiSolutionsGappy;

  const char *prefix;
  const char *geometryprefix;
  const char *connectivity;
  const char *geometry;
  const char *decomposition;
  const char *cpumap;
  const char *match;
  const char *embmeshmatch;
  const char *embsurfmatch;
  const char *d2wall;
  const char *perturbed;
  const char *solutions;
  const char *referenceSolution;
  const char *multiSolutions;
  const char *multiSolutionsParams; // ROMs: path to a file listing file paths to solutions (not used) and their corresponding operating point (solutions are assumed to be in the same order as in the multiSolutions file)
  const char *parameters;  // ROMs: operating point for the current simulation (interpolates multiSolutionsParams to find an IC);
  double parametricDistanceExponent;
  int maxInterpolatedSolutions;
  const char *positions;
  const char *displacements;
  const char *embeddedpositions;
  const char *levelsets;
  const char *cracking;

  // We can now read the fluid ID from a file.
  // Added by Alex Main (September 2013)
  //
  const char *fluidId;
  const char *rstdata;

  // File Package for restart support.
  // Added by Alex Main (September 2013)
  //
  const char *restart_file_package;

  const char *podFile;
  const char *optimalPressureFile;
  const char *matchStateFile;
  const char *optimalPressureDim;
  const char *stateSnapFile;
  const char *stateSnapRefSolution;
    const char *stateMaskSnapFile; //<! for snapshot with embedded methods, Lei Lei, 02/03/2016
  const char *multiStateSnapRefSolution;
  const char *residualSnapFile;
  const char *krylovSnapFile;
  const char *sensitivitySnapFile;
  const char *approxMetricStateSnapFile;
  const char *approxMetricNonlinearSnapFile;
  const char *greedyDataFile;
  const char *projErrorSnapFile;
  const char *initialClusterCentersFile;
  const char *reducedCoords;
  const char *strModesFile;
  const char *embeddedSurface;
  const char *oneDimensionalSolution;

  // Paths to external executables.  Currently only hooked up for Gappy preprocessing, but could be used for general simulations.
  const char *sower;
  const char *metis;
  int nParts;

  //const char *stateVecFile;//CBM

  const char* convergence_file;

  const char* exactInterfaceLocation;

  //
  // String for the input files of pressure snapshots
  // This file is used for computing the Kirchhoff integral.
  // UH (08/2012)
  const char* strKPtraces;

  const char *wallsurfacedisplac; //YC
  const char *reducedEigState;
// Included (MB)
  const char *shapederivatives;
  const char *shapederivativestype;

  ObjectMap< OneDimensionalInputData > oneDimensionalInput;

  InputData();
  ~InputData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct Probes {

  const static int MAXNODES = 3176;
  struct Node {
    Node() { id = -1; locationX = locationY = locationZ = -1.0e20; subId = localNodeId = -1;
             isLocationBased = false; }
    int id;
    int subId;
    int localNodeId;
    double locationX,locationY,locationZ;
    bool isLocationBased;
    void setup(const char *, ClassAssigner * = 0);
  };

  Node myNodes[MAXNODES];

  const char *prefix;
  const char *density;
  const char *pressure;
  const char *diffpressure;
  const char *temperature;
  const char *velocity;
  const char *displacement;

  // UH >> for Aeroacoustic
  const char *farfieldpattern;

  Probes();
  ~Probes() {}

  void setup(const char *, ClassAssigner * = 0);
};

struct LinePlot {

  LinePlot();
  ~LinePlot() {}

  Assigner *getAssigner();
  double x0,y0,z0;
  double x1,y1,z1;

  int numPoints;

  const char *density;
  const char *pressure;
  const char *temperature;
  const char *velocity;
  //const char *displacement;
};

//------------------------------------------------------------------------------

struct TransientData {

  const char *prefix;
  const char *density;
  const char *tavdensity;
  const char *mach;
  const char *speed;
  const char *wtmach;
  const char *wtspeed;
  const char *absvelocity;
  const char *tavmach;
  const char *pressure;
  const char *diffpressure;
  const char *tavpressure;
  const char *hydrostaticpressure;
  const char *hydrodynamicpressure;
  const char *pressurecoefficient;
  const char *temperature;
  const char *tavtemperature;
  const char *totalpressure;
  const char *tavtotalpressure;
  const char *vorticity;
  const char *tavvorticity;
  const char *nutturb;
  const char *kturb;
  const char *epsturb;
  const char *eddyvis;
  const char *dplus;
  const char *sfric;
  const char *tavsfric;
  const char *psensor;
  const char *csdles;
  const char *tavcsdles;
  const char *csdvms;
  const char *tavcsdvms;
  const char *mutOmu;
  const char *velocity;
  const char *tavvelocity;
  const char *displacement;
  const char *tavdisplacement;
  const char *flightDisplacement;
  const char *localFlightDisplacement;
  const char *forces;
  const char *tavforces;
  const char *hydrostaticforces;
  const char *hydrodynamicforces;
  const char *generalizedforces;
  const char *lift;
  const char *matchpressure;
  const char *matchstate;
  const char *fluxnorm;
  const char *tavlift;
  const char *hydrostaticlift;
  const char *hydrodynamiclift;
  const char *residuals;
  const char *materialVolumes;
  const char *materialMassEnergy;
  const char *conservation;
  const char *podFile;
  const char *robProductFile;
  const char *rMatrixFile;
  const char *romFile;
  const char *gendispFile;
  const char *romInitialConditionFile;
  const char *philevel;
  const char *philevel2;
  const char *controlvolume;
  const char* fluidid;
  const char* d2wall;
  const char *embeddedsurface;
  const char *cputiming;
  const char *aeroelasticEigenvalues;
  const char *gamData;
  const char *gamFData;


// Included (MB)
  const char *velocitynorm;
  const char *dSpatialres;
  const char *dSpatialresnorm;
  const char *dSolutions;
  const char *dDensity;
  const char *dMach;
  const char *dPressure;
  const char *dMatchPressure;
  const char *dTemperature;
  const char *dTotalpressure;
  const char *dNutturb;
  const char *dEddyvis;
  const char *dVelocityScalar;
  const char *dVelocityVector;
  const char *dDisplacement;
  const char *dForces;
  const char *dLiftDrag;
  const char *dLiftx;
  const char *dLifty;
  const char *dLiftz;

  const char *tempnormalderivative;
  const char *surfaceheatflux;
  const char *heatfluxes;

  const char *sparseGrid;

  // For 1D solver
  const char* bubbleRadius;

  const char* multiSolnFluxNorm;

  int frequency;
  double x0, y0, z0;
  double length;
  double surface;
  double frequency_dt; //set to -1.0 by default. Used iff it is activated (>0.0) by user.

  Probes probes;

  ObjectMap<LinePlot> linePlots;

  TransientData();
  ~TransientData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct RestartData {

  enum Type {SINGLE = 0, DOUBLE = 1} type;

  const char *prefix;

  const char *solutions;
  const char *positions;
  const char *embeddedpositions;
  const char *levelsets;
  const char *cracking;
  const char *fluidId;
  const char *data;
  const char* filepackage;

  int frequency;
  double frequency_dt; //set to -1.0 by default. Used iff it is activated (>0.0) by user.

  /// UH (06/2012)
  ///
  /// The following member is used for computing the Kirchhoff integral.
  /// When active, this variable contains the prefix for saving the data.
  ///
  const char *strKPtraces;

  RestartData();
  ~RestartData() {}

  void setup(const char *, ClassAssigner * = 0);
};
//------------------------------------------------------------------------------

struct ROMOutputData {

  const char *prefix;

  const char *dFluxNorm;

  const char *stateVector;
    // lei lei, 02/01/2016: added mask vector for embedded Ts
    const char *stateMaskVector;
  int stateOutputFreqTime;
  int stateOutputFreqNewton;
  enum AvgStateIncrements {AVG_STATE_INCREMENTS_OFF = 0, AVG_STATE_INCREMENTS_ON = 1} avgStateIncrements;

  const char *residualVector;
  int residualOutputFreqTime;
  int residualOutputMaxNewton;
  enum FDResiduals {FD_RESIDUALS_OFF = 0, FD_RESIDUALS_ON = 1} fdResiduals;
  enum FDResidualsLimit {FD_RESIDUALS_LIMIT_OFF = 0, FD_RESIDUALS_LIMIT_ON = 1} fdResidualsLimit;
  enum OutputOnlySpatialResidual {OUTPUT_ONLY_SPATIAL_RES_OFF = 0, OUTPUT_ONLY_SPATIAL_RES_ON = 1} outputOnlySpatialResidual;

  const char *krylovVector;
  int krylovOutputFreqTime;
  int krylovOutputFreqNewton;
  double krylovVectorEnergy;
  enum AddStateToKrylov {ADD_STATE_TO_KRYLOV_OFF = 0, ADD_STATE_TO_KRYLOV_ON = 1} addStateToKrylov;

  const char *clusterUsage;
  const char *reducedCoords;  // generalized coords
  const char *dUnormAccum;

  const char *residualsForCoordRange;

  int resjacfrequency;

  enum OverwriteNonlinearSnaps {OVERWRITE_OFF = 0, OVERWRITE_ON = 1} overwriteNonlinearSnaps;

  ROMOutputData();
  ~ROMOutputData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct OutputData {

  TransientData transient;
  RestartData restart;
  ROMOutputData rom;

  OutputData();
  ~OutputData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct RestartParametersData {

  int iteration;

  double etime;
  double dt_nm1;
  double dt_nm2;
  double residual;
  double energy;
  double outputNewtonTag;
  int outputNewtonStateStep;
  int outputNewtonResidualStep;
  int outputKrylovStep;

  RestartParametersData();
  ~RestartParametersData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ProblemData {

  enum Type {UNSTEADY = 0, ACCELERATED = 1, AERO = 2, THERMO = 3, FORCED = 4,
             ROLL = 5, RBM = 6, LINEARIZED = 7, NLROMOFFLINE = 8, NLROMONLINE = 9, SIZE = 10};
  bool type[SIZE];

  enum AllType {_STEADY_                   = 0,
              _UNSTEADY_                   = 1,
              _ACC_UNSTEADY_               = 2,
              _STEADY_AEROELASTIC_         = 3,
              _UNSTEADY_AEROELASTIC_       = 4,
              _ACC_UNSTEADY_AEROELASTIC_   = 5,
              _STEADY_THERMO_              = 6,
              _UNSTEADY_THERMO_            = 7,
              _STEADY_AEROTHERMOELASTIC_   = 8,
              _UNSTEADY_AEROTHERMOELASTIC_ = 9,
              _FORCED_                     =10,
              _ACC_FORCED_                 =11,
              _ROLL_                       =12,
              _RBM_                        =13,
              _UNSTEADY_LINEARIZED_AEROELASTIC_ = 14,
              _UNSTEADY_LINEARIZED_        =15,
              _NONLINEAR_ROM_OFFLINE_      =16,
              _ROM_AEROELASTIC_            =17,
              _ROM_                        =18,
              _FORCED_LINEARIZED_          =19,
              _INTERPOLATION_              =20,
              _NONLINEAR_EIGEN_ERROR_INDICATOR_ = 21,
              _SPARSEGRIDGEN_              =22,
              _ONE_DIMENSIONAL_            =23,
              _UNSTEADY_NONLINEAR_ROM_     =24,
              _NONLINEAR_ROM_PREPROCESSING_ = 25,
              _SURFACE_MESH_CONSTRUCTION_  =26,
              _SAMPLE_MESH_SHAPE_CHANGE_   =27,
              _UNSTEADY_NONLINEAR_ROM_POST_ = 28,
              _POD_CONSTRUCTION_           =29,
              _ROB_INNER_PRODUCT_          =30,
              _AERO_ACOUSTIC_              =31,
              _SHAPE_OPTIMIZATION_         =32,
              _AEROELASTIC_SHAPE_OPTIMIZATION_ =  33,
              _AEROELASTIC_ANALYSIS_       =34,
              _GAM_CONSTRUCTION_           =35,
              _ACC_UNSTEADY_NONLINEAR_ROM_ =36,
              _STEADY_NONLINEAR_ROM_       =37,
              _FORCED_NONLINEAR_ROM_       =38,
              _ROM_SHAPE_OPTIMIZATION_     =39,
              _STEADY_NONLINEAR_ROM_POST_  =40,
              _EMBEDDED_ALS_ROM_           =41,
              _EMBEDDED_ALS_ROM_ONLINE_    =42,
              _SENSITIVITY_ANALYSIS_       =43} alltype;
  enum Mode {NON_DIMENSIONAL = 0, DIMENSIONAL = 1} mode;
  enum Test {REGULAR = 0} test;
  enum Prec {NON_PRECONDITIONED = 0, PRECONDITIONED = 1} prec;
  enum Framework {BODYFITTED = 0, EMBEDDED = 1, EMBEDDEDALE = 2} framework;
  enum SolveFluid {OFF = 0, ON = 1} solvefluid;
  enum SolutionMethod { TIMESTEPPING = 0, MULTIGRID = 1} solutionMethod;
  int verbose;

  enum SolveWithMultipleICs {_MULTIPLE_ICS_FALSE_ = 0, _MULTIPLE_ICS_TRUE_ = 1} solveWithMultipleICs;

  ProblemData();
  ~ProblemData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct PreconditionData {

  double mach;
  double cmach;
  double k;
  double betav;
  double shockreducer;

  PreconditionData();
  ~PreconditionData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ReferenceStateData {

  double mach;
  double velocity;
  double density;
  double pressure;
  double temperature;
  double reynolds_mu;
  double energy;
  double length;

// Included (MB)
  double dRe_mudMach;

  RefVal rv;

  ReferenceStateData();
  ~ReferenceStateData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct BcsFreeStreamData {

  enum Type {EXTERNAL = 0, INTERNAL = 1} type;

  double mach;
  double velocity;
  double density;
  double pressure;
  double temperature;
  double nutilde;
  double kenergy;
  double eps;
  double alpha;
  double beta;

  BcsFreeStreamData();
  ~BcsFreeStreamData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct BcsWallData {

  enum Type {ISOTHERMAL = 0, ADIABATIC = 1} type;
  enum Integration {AUTO = 0, WALL_FUNCTION = 1, FULL = 2} integration;
  enum Reconstruction {CONSTANT = 0, EXACT_RIEMANN = 1} reconstruction;

  double temperature;
  double delta;
	bool delta_given;

  BcsWallData();
  ~BcsWallData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct BcsHydroData {

  double depth;

  BcsHydroData();
  ~BcsHydroData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct BoundaryData  {

  static const int UNSPECIFIED = -1;
  enum Type {DIRECTSTATE = 1, MASSFLOW = 2, POROUSWALL = 3} type;

   enum vars {DENSITY = 0, VX = 1, VY = 2, VZ = 3, PRESSURE = 4, TEMPERATURE = 5, TOTALPRESSURE = 6, TOTALTEMPERATURE = 7, MDOT = 8, NUTILDE = 9, KENERGY = 10, EPSILON = 11, SIZE = 12};
  bool inVar[SIZE], outVar[SIZE];
  double density;
  double velocityX, velocityY, velocityZ;
  double pressure;
  double temperature;
  double totalPressure;
  double totalTemperature;
  double mdot;
  double nutilde;
  double kenergy;
  double epsilon;
  double porosity;

  BoundaryData();
  Assigner *getAssigner();

};

//-----------------------------------------------------------------------------

struct BcsData {

  BcsFreeStreamData inlet;
  BcsFreeStreamData outlet;
  BcsWallData wall;
  BcsHydroData hydro;
  ObjectMap<BoundaryData> bcMap;

  BcsData();
  ~BcsData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct GasModelData {

  enum Type {IDEAL = 0, STIFFENED = 1} type;

  double specificHeatRatio;
  double idealGasConstant;
  double pressureConstant;
  double specificHeatPressure;

  GasModelData();
  ~GasModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct JWLModelData {

  enum Type {IDEAL = 0, JWL = 1} type;

  double omega; // = specificHeatRatio-1.0
  double idealGasConstant;
  double A1,R1,rhoref,A2,R2;

  JWLModelData();
  ~JWLModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct LiquidModelData {

  enum Type { COMPRESSIBLE = 0 } type;

  enum YesNo {YES = 0, NO = 1 };
  YesNo check;

  YesNo burnable;

  // the state equation is derived from a linearization of the bulk modulus wrt
  // pressure: K = k1 + k2 * P
  // the integration constant of the ODE is given by the couple (RHOrefwater,Prefwater)
  double specificHeat;
  double k1water;
  double k2water;
  double Bwater;
  double Prefwater;
  double RHOrefwater;

  //the state equation can be put in the form P=Pref+alpha*rho^beta
  // with Pref, alpha and beta derived from k1, k2 and the 'initial' couple
  double Pref;
  double alpha;
  double beta;

  LiquidModelData();
  ~LiquidModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct FluidModelData {

  enum Fluid { PERFECT_GAS = 0, LIQUID = 1, JWL = 2, STIFFENED_GAS = 3, UNDEFINED = 4} fluid;
  double rhomin;
  double pmin;

  GasModelData gasModel;
  JWLModelData jwlModel;
  LiquidModelData liquidModel;

  FluidModelData();
  ~FluidModelData() {}

  Assigner *getAssigner();
  void setup(const char *, ClassAssigner * = 0);
  //FluidModelData &operator=(const FluidModelData &);

};

//------------------------------------------------------------------------------

struct ViscosityModelData {

  enum Type {CONSTANT = 0, SUTHERLAND = 1, PRANDTL = 2} type;

  double sutherlandReferenceTemperature;
  double sutherlandConstant;
  double dynamicViscosity;
  double bulkViscosity;

  ViscosityModelData();
  ~ViscosityModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ThermalCondModelData {

  enum Type {CONSTANT_PRANDTL = 0, CONSTANT = 1} type;

  double prandtl;
  double conductivity;

  ThermalCondModelData();
  ~ThermalCondModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct PorousMedia  {

  double iprimex, iprimey, iprimez;
  double jprimex, jprimey, jprimez;
  double kprimex, kprimey, kprimez;
  double alphax,alphay,alphaz;
  double betax,betay,betaz;
  double idr,ldr;

  PorousMedia();
  //Assigner *getAssigner();
  void setup(const char*, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct InitialConditions {

  double mach;
  double velocity;
  double alpha, beta;
  double pressure;
  double density;
  double temperature;

  InitialConditions();
  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct VolumeData  {

  enum Type {FLUID = 0, POROUS = 1} type;
  int fluidModelID;

  PorousMedia   porousMedia;
  InitialConditions initialConditions;

  VolumeData();
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct SAModelData {

  double cb1;
  double cb2;
  double cw2;
  double cw3;
  double cv1;
  double cv2;
  double sigma;
  double vkcst;
  enum Form {ORIGINAL = 0, FV3 = 1} form;

  SAModelData();
  ~SAModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct DESModelData {

  double cb1;
  double cb2;
  double cw2;
  double cw3;
  double cv1;
  double cv2;
  double cdes;
  double sigma;
  double vkcst;
  enum Form {ORIGINAL = 0, FV3 = 1} form;

  DESModelData();
  ~DESModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct KEModelData {

  double sigma_k;
  double sigma_eps;
  double sigma_eps1;
  double sigma_eps2;
  double c_mu;

  KEModelData();
  ~KEModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------
//
struct WallDistanceMethodData {

  enum Type {ITERATIVE = 0, NONITERATIVE = 1, HYBRID = 2} type;

  int maxIts;
  double eps;
  int iterativelvl;

  WallDistanceMethodData();
  ~WallDistanceMethodData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct TurbulenceModelData {

  enum Type {ONE_EQUATION_SPALART_ALLMARAS = 0, ONE_EQUATION_DES = 1, TWO_EQUATION_KE = 2} type;

  SAModelData sa;
  DESModelData des;
  KEModelData ke;
  WallDistanceMethodData d2wall;

  TurbulenceModelData();
  ~TurbulenceModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SmagorinskyLESData {

  double c_s;

  SmagorinskyLESData();
  ~SmagorinskyLESData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct WaleLESData {

  double c_w;

  WaleLESData();
  ~WaleLESData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ClippingData {


  double cs_max;
  double pt_min;
  double pt_max;

  ClippingData();
  ~ClippingData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct DynamicLESData {

  ClippingData clip;

  DynamicLESData();
  ~DynamicLESData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------


struct DynamicVMSData {

  enum Type {D1VMSLES = 0, D2VMSLES = 1, D3VMSLES = 2} type;

  double c_s_prime;
  int agglomeration_width;
  int agglomeration_depth1;
  int agglomeration_depth2;

  ClippingData clip;

  DynamicVMSData();
  ~DynamicVMSData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct VMSLESData {

  double c_s_prime;
  int agglomeration_width;
  int agglomeration_depth;

  VMSLESData();
  ~VMSLESData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct LESModelData {

  enum Type {SMAGORINSKY = 0, DYNAMIC = 1, VMS = 2, DYNAMICVMS = 3, WALE = 4} type;
  enum Delta {VOLUME = 0, SIDE = 1} delta;

  SmagorinskyLESData sma;
  DynamicLESData dles;
  VMSLESData vms;
  DynamicVMSData dvms;
  WaleLESData wale;

  LESModelData();
  ~LESModelData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct TBFixData {

  double x0;
  double y0;
  double z0;
  double x1;
  double y1;
  double z1;

  TBFixData();
  ~TBFixData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct TripDomainData {

  TBFixData bfix;

  TripDomainData();
  ~TripDomainData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct TurbulenceClosureData {

  enum Type {NONE = 0, EDDY_VISCOSITY = 1, LES = 2} type;

  double prandtlTurbulent;
  TurbulenceModelData tm;
  LESModelData les;
  TripDomainData tr;

  TurbulenceClosureData();
  ~TurbulenceClosureData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ProgrammedBurnData {

  int unburnedEOS,burnedEOS;
  double ignitionX0,ignitionY0,ignitionZ0;
  double e0;
  double cjDetonationVelocity;
  double cjPressure;
  double cjDensity;
  double cjEnergy;
  double ignitionTime;
  double factorB;
  double factorS;
  //double stopWhenShockReachesPercentDistance;
  int ignited;
  int limitPeak;

  ProgrammedBurnData();
  ~ProgrammedBurnData();

  void setup(const char*, ClassAssigner* = 0);

};

struct SphereData {

  double cen_x, cen_y, cen_z, radius;
  int fluidModelID;
  InitialConditions initialConditions;

  ProgrammedBurnData programmedBurn;

  SphereData();
  ~SphereData() {}
  Assigner *getAssigner();

};
//------------------------------------------------------------------------------
struct PrismData {

  double cen_x, cen_y, cen_z, w_x,w_y,w_z;
  double X0,Y0,Z0,X1,Y1,Z1;
  int fluidModelID;
  InitialConditions initialConditions;

  ProgrammedBurnData programmedBurn;

  bool inside(double x,double y,double z) const;

  PrismData();
  ~PrismData() {}
  Assigner *getAssigner();

};
//------------------------------------------------------------------------------

struct PlaneData {

  double cen_x, cen_y, cen_z, nx, ny, nz;
  int fluidModelID;
  InitialConditions initialConditions;

  PlaneData();
  ~PlaneData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct CylinderData {

  double cen_x, cen_y, cen_z, nx, ny, nz, r, L;
  int fluidModelID;
  InitialConditions initialConditions;

  CylinderData();
  ~CylinderData() {}
  Assigner *getAssigner();

};

//------------------------------------------------------------------------------

struct PointData {

  int fluidModelID;
  double x,y,z;
  InitialConditions initialConditions;

  ProgrammedBurnData programmedBurn;

  PointData();
  ~PointData() {}
  Assigner *getAssigner();

};

struct DummyPointData {

  int fluidModelID;

  DummyPointData();
  ~DummyPointData() {}
  Assigner *getAssigner();

};
//------------------------------------------------------------------------------

struct MultiInitialConditionsData {

  ObjectMap<SphereData> sphereMap;
  ObjectMap<PrismData>  prismMap;
  ObjectMap<PlaneData>  planeMap;
  ObjectMap<CylinderData>  cylinderMap;
  ObjectMap<PointData>  pointMap;
  ObjectMap<DummyPointData>  dummyPointMap;

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct SparseGridData {

  // to use already created sparse grids
  const char *tabulationFileName;
  int numberOfTabulations;

  // to generate sparse grids
  int verbose;

  int minPoints;
  int maxPoints;
  double relAccuracy;
  double absAccuracy;
  double dimAdaptDegree;

  double range1min, range1max, mapBaseValue1; int numDomainDim1;
  double range2min, range2max, mapBaseValue2; int numDomainDim2;
  double range3min, range3max, mapBaseValue3; int numDomainDim3;
  double range4min, range4max, mapBaseValue4; int numDomainDim4;
  double range5min, range5max, mapBaseValue5; int numDomainDim5;
  double range6min, range6max, mapBaseValue6; int numDomainDim6;
  typedef double Range[2];
  Range *range;
  double *mapBaseValue;
  int *numDomainDim;

  int numOutputs;
  int numInputs;

  SparseGridData();
  ~SparseGridData() {}
  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct MultiFluidData {
  enum Method {NONE = 0, GHOSTFLUID_FOR_POOR = 1, GHOSTFLUID_WITH_RIEMANN} method;
  enum InterfaceTracking {LINEAR = 0, GRADIENT = 1, HERMITE = 2};
  enum RiemannComputation {FE = 0, RK2 = 1, TABULATION2 = 2, TABULATION5 = 3} riemannComputation;
  int bandlevel;
  int frequency;
  double eps;
  int outputdiff;
  double jwlRelaxationFactor;
  enum TypePhaseChange {ASIS = 0, RIEMANN_SOLUTION = 1, EXTRAPOLATION = 2} typePhaseChange;
  enum CopyCloseNodes {FALSE = 0, TRUE = 1} copyCloseNodes;
  enum InterfaceType {FSF = 0, FF = 1, FSFandFF = 2} interfaceType;

  enum InterfaceTreatment {FIRSTORDER=0, SECONDORDER=1} interfaceTreatment;
  enum InterfaceExtrapolation {EXTRAPOLATIONFIRSTORDER=0, EXTRAPOLATIONSECONDORDER=1, AUTO=2} interfaceExtrapolation;
  enum InterfaceLimiter {LIMITERNONE = 0, LIMITERALEX1 = 1} interfaceLimiter;

  // TRIANGULATED refers to some tests Alex was doing where the multifluid interface is represented by a triangulated surface
  //    (akin to a massless embedded structure).  It is mostly implemented, but he left out a few parts he didn't need for his tests.
  enum LevelSetMethod { CONSERVATIVE = 0, HJWENO = 1, SCALAR=2, PRIMITIVE = 3,
                        TRIANGULATED = 4} levelSetMethod;

  enum RiemannNormal {REAL = 0, MESH = 1, LEGACYMESH = 2 } riemannNormal;
  double riemannEps;
  int riemannMaxIts;

  enum Prec {NON_PRECONDITIONED = 0, PRECONDITIONED = 1, SAME_AS_PROBLEM = 2} prec;

  MultiInitialConditionsData multiInitialConditions;

  SparseGridData sparseGrid;

  int testCase;

  MultiFluidData();
  ~MultiFluidData() {}
  void setup(const char *, ClassAssigner * = 0);
};


//------------------------------------------------------------------------------

struct Volumes  {

  ObjectMap<VolumeData> volumeMap;

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct EquationsData {

  int dimension;

  enum Type {EULER = 0, NAVIER_STOKES = 1} type;

  int numPhase;

// it is assumed that in a two-phase flow, fluidModel represents the surrounding fluid
// whereas fluidModel2 represents the fluid in the bubbles. This means that the
// uniform boundary conditions and the uniform initial conditions are computed
// with the values given to characterize the first fluid!

  double gravity_x, gravity_y, gravity_z;

  ObjectMap<FluidModelData> fluidModelMap;
  FluidModelData fluidModel;
  ViscosityModelData viscosityModel;
  ThermalCondModelData thermalCondModel;
  TurbulenceClosureData tc;

  EquationsData();
  ~EquationsData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SchemeData {

  enum AdvectiveOperator {FINITE_VOLUME = 0, FE_GALERKIN = 1} advectiveOperator;
  enum Flux {ROE = 0, VANLEER = 1, HLLE = 2, HLLC = 3} flux;

  enum Reconstruction {CONSTANT = 0, LINEAR = 1} reconstruction;

  enum Limiter {NONE = 0, VANALBADA = 1, BARTH = 2, VENKAT = 3, P_SENSOR = 4,
                EXTENDEDVANALBADA = 5} limiter;
  enum Gradient {LEAST_SQUARES = 0, GALERKIN = 1, NON_NODAL = 2} gradient;
  enum Dissipation {SECOND_ORDER = 0, SIXTH_ORDER = 1} dissipation;

  double beta;
  double gamma;
  double xiu;
  double xic;
  double eps;

  double xirho;
  double xip;
  double vel_fac;

  struct MaterialFluxData {

    Flux flux;

    Assigner *getAssigner();
  };

  // We now allow different flux functions to be used for different materials.
  // The behavior is that if the flux is specified for a fluid id in this map,
  // then it is used.  Otherwise, the default (schemedata.flux) is used for
  // that material.
  ObjectMap<MaterialFluxData> fluxMap;

  int allowsFlux;

  // allowsFlux = 0 for levelset equation (the choice of flux for the levelset is hardcoded)
  SchemeData(int allowsFlux = 1);
  ~SchemeData() {}

  void setup(const char *, ClassAssigner * = 0);

};
//------------------------------------------------------------------------------

struct CFixData {

  // nodal coordinates of first circle
  double x0;
  double y0;
  double z0;

  // nodal coordinates of second circle
  double x1;
  double y1;
  double z1;

  // radii of 1st and 2nd circle
  double r0;
  double r1;

  CFixData();
  ~CFixData() {}

  int failsafeN;
  enum {OFF=0, ON=1, ALWAYSON=2} failsafe;

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SFixData {

  double x0;
  double y0;
  double z0;
  double r;

  int failsafeN;
  enum {OFF=0, ON=1, ALWAYSON=2} failsafe;

  SFixData();
  ~SFixData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct BFixData {

  double x0;
  double y0;
  double z0;
  double x1;
  double y1;
  double z1;

  int failsafeN;
  enum {OFF=0, ON=1, ALWAYSON=2} failsafe;

  BFixData();
  ~BFixData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SchemeFixData {

  static const int num = 10;

  enum Symmetry {NONE = 0, X = 1, Y = 2, Z = 3} symmetry;
  double dihedralAngle;

  SFixData* spheres[num];
  BFixData* boxes[num];
  CFixData* cones[num];

  SFixData sfix1;
  SFixData sfix2;
  SFixData sfix3;
  SFixData sfix4;
  SFixData sfix5;
  SFixData sfix6;
  SFixData sfix7;
  SFixData sfix8;
  SFixData sfix9;
  SFixData sfix10;

  BFixData bfix1;
  BFixData bfix2;
  BFixData bfix3;
  BFixData bfix4;
  BFixData bfix5;
  BFixData bfix6;
  BFixData bfix7;
  BFixData bfix8;
  BFixData bfix9;
  BFixData bfix10;

  CFixData cfix1;
  CFixData cfix2;
  CFixData cfix3;
  CFixData cfix4;
  CFixData cfix5;
  CFixData cfix6;
  CFixData cfix7;
  CFixData cfix8;
  CFixData cfix9;
  CFixData cfix10;

  SchemeFixData();
  ~SchemeFixData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct BoundarySchemeData {

  enum Type { STEGER_WARMING = 0,
              CONSTANT_EXTRAPOLATION = 1,
              LINEAR_EXTRAPOLATION = 2,
              GHIDAGLIA = 3, MODIFIED_GHIDAGLIA = 4} type;

  BoundarySchemeData();
  ~BoundarySchemeData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SchemesData {

  SchemeData ns;
  SchemeData tm;
  SchemeData ls;
  SchemeFixData fixes;
  BoundarySchemeData bc;

  SchemesData();
  ~SchemesData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ExplicitData {

//time-integration scheme used
  enum Type {RUNGE_KUTTA_4 = 0, RUNGE_KUTTA_2 = 1, FORWARD_EULER = 2, ONE_BLOCK_RK2 = 3, ONE_BLOCK_RK2bis = 4} type;

  ExplicitData();
  ~ExplicitData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct PcData {

  enum Type {IDENTITY = 0, JACOBI = 1, AS = 2, RAS = 3, ASH = 4, AAS = 5, MG = 6} type;
  enum Renumbering {NATURAL = 0, RCM = 1} renumbering;

  enum MGSmoother { MGJACOBI = 0, MGLINEJACOBI = 1, MGRAS = 2 } mg_smoother;

  enum MGType { MGALGEBRAIC = 0, MGGEOMETRIC = 1} mg_type;

  int fill;

  int num_multigrid_smooth1,num_multigrid_smooth2;
  int num_multigrid_levels;

  int mg_output;

  double mg_smooth_relax;

  int num_fine_sweeps;

  PcData();
  ~PcData() {}

  void setup(const char *, ClassAssigner * = 0);

};

struct MultiGridData {

  enum MGSmoother { MGJACOBI = 0, MGLINEJACOBI = 1, MGRAS = 2, MGGMRES = 3 } mg_smoother;

  enum CycleScheme { VCYCLE = 0, WCYCLE = 1} cycle_scheme;

  enum RestrictMethod { VOLUME_WEIGHTED = 0, AVERAGE = 1 } restrictMethod;

  enum CoarseningRatio { TWOTOONE = 0, FOURTOONE = 1} coarseningRatio;

  int num_multigrid_smooth1,num_multigrid_smooth2;
  int num_multigrid_levels;

  int mg_output;

  int useGMRESAcceleration;

  double directional_coarsening_factor;

  double mg_smooth_relax;

  double prolong_relax_factor,restrict_relax_factor;

  int num_fine_sweeps;

  int addViscousTerms;

  int addTurbulenceTerms;

  SchemeFixData fixes;

  const char* agglomerationFile;

  double turbRelaxCutoff;

  double densityMin,densityMax;

  MultiGridData();
  ~MultiGridData() {}

  void setup(const char *, ClassAssigner * = 0);


};

//------------------------------------------------------------------------------

struct KspData {

  enum Type {RICHARDSON = 0, CG = 1, GMRES = 2, GCR = 3} type;
  enum EpsFormula {CONSTANT = 0, EISENSTADT = 1} epsFormula;
  enum CheckFinalRes {NO = 0, YES = 1} checkFinalRes;

  int maxIts;
  int numVectors;
  double eps;

  double absoluteEps;

  const char *output;

  PcData pc;

  KspData();
  ~KspData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct KspFluidData {

  KspData ns;
  KspData tm;
  KspData lsi;

  KspFluidData();
  ~KspFluidData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct LineSearchData {

  enum Type {NONE = 0, BACKTRACKING = 1} type;
  int maxIts;
  double rho;
  double c1;

  LineSearchData();
  ~LineSearchData() {}

  void setup(const char *, ClassAssigner * = 0);

};




//------------------------------------------------------------------------------

template<class GenericKrylov>
struct NewtonData {

  enum FailSafe {NO = 0, YES = 1, ALWAYS = 2} failsafe;
  int maxIts;
  double eps;
  int JacSkip;
  double epsAbsRes, epsAbsInc;
  GenericKrylov ksp;
  LineSearchData lineSearch;
  const char *output;

  NewtonData();
  ~NewtonData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ImplicitData {

  enum Type {BACKWARD_EULER = 0, CRANK_NICOLSON = 1, THREE_POINT_BDF = 2, FOUR_POINT_BDF = 3, SPATIAL_ONLY = 4} type;
  enum Startup {REGULAR = 0, MODIFIED = 1} startup;
  enum TurbulenceModelCoupling {WEAK = 0, STRONG = 1} tmcoupling;
  enum Mvp {FD = 0, H1 = 1, H2 = 2, H1FD = 3} mvp;
  enum FiniteDifferenceOrder {FIRST_ORDER = 1, SECOND_ORDER = 2} fdOrder;
  enum FVMERS3PBDFSchme { BDF_SCHEME1 = 1, BDF_SCHEME2 = 0 } fvmers_3pbdf;
  NewtonData<KspFluidData> newton;
  /// UH (09/10)
  /// This flag is not visible from the input file.
  /// It governs the computation of the Jacobian of the flux function,
  /// a component of the 'H' matrix (from the MatrixVectorProduct).
  enum FluxFcnJacobian {FINITE_DIFFERENCE = 0, APPROXIMATE = 1, EXACT = 2} ffjacobian;

  ImplicitData();
  ~ImplicitData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct CFLData {

  enum Strategy {DEFAULT = -1, RESIDUAL = 0, DIRECTION = 1, DFT = 2, HYBRID = 3, FIXEDUNSTEADY = 4, OLD = 5} strategy;

  // global cfl parameters
  double cfl0;
  double cflCoef1;
  double cflCoef2;
  double cflMax;
  double cflMin;
  double dualtimecfl;

  // residual based parameters
  double ser;

  // direction based parameters
  double angle_growth;
  double angle_zero;

  // dft based parameters
  int dft_history;
  int dft_freqcutoff;
  double dft_growth;

  const char *output;

  CFLData();
  ~CFLData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct TsData {

  enum Type {EXPLICIT = 0, IMPLICIT = 1} type;
  enum TypeTimeStep {AUTO = 0, LOCAL = 1, GLOBAL = 2} typeTimeStep;
  enum Clipping {NONE = 0, ABS_VALUE = 1, FREESTREAM = 2, CUTOFF = 3} typeClipping;
  enum TimeStepCalculation {CFL = 0, ERRORESTIMATION = 1} timeStepCalculation;
  enum DualTimeStepping {OFF = 0, ON = 1} dualtimestepping;

  enum Prec {NO_PREC = 0, PREC = 1} prec;
  enum Form {DESCRIPTOR = 1, NONDESCRIPTOR = 0, HYBRID = 2} form;
  double viscousCst;

  int maxIts;
  double eps;
  double epsabs;
  double timestep;
  double timestepinitial;
  double maxTime;

  int residual;
  double errorTol;

  // Kept for back compatibility
  double cfl0;
  double cflCoef1;
  double cflCoef2;
  double cflMax;
  double cflMin;
  double ser;
  double dualtimecfl;

  int checksol;
  int checkvelocity;
  int checkpressure;
  int checkdensity;
  int checklinsolve;
  int deltapressurethreshold;
  int deltadensitythreshold;

  double programmedBurnShockSensor;
  double rapidPressureThreshold;
  double rapidDensityThreshold;

  const char *output;

  ExplicitData expl;
  ImplicitData implicit;
  CFLData cfl;

  TsData();
  ~TsData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct DGCLData{

  enum Normals {AUTO = 0, IMPLICIT_FIRST_ORDER_GCL = 1, IMPLICIT_SECOND_ORDER_GCL = 2,
                IMPLICIT_FIRST_ORDER_EZGCL = 3, IMPLICIT_SECOND_ORDER_EZGCL = 4, IMPLICIT_THIRD_ORDER_EZGCL = 5,
                IMPLICIT_CURRENT_CFG = 6, IMPLICIT_LATEST_CFG = 7, EXPLICIT_RK2 = 8} normals;
  enum Velocities {AUTO_VEL = 0, IMPLICIT_BACKWARD_EULER_VEL = 1, IMPLICIT_THREE_POINT_BDF_VEL = 2,
                   IMPLICIT_IMPOSED_VEL = 3, IMPLICIT_IMPOSED_BACKWARD_EULER_VEL = 4,
                   IMPLICIT_IMPOSED_THREE_POINT_BDF_VEL = 5, IMPLICIT_ZERO = 6, EXPLICIT_RK2_VEL = 7} velocities;

  DGCLData();
  ~DGCLData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

// Included (MB)
struct SensitivityAnalysis {

  enum Method {DIRECT = 0, ADJOINT = 1} method;
  enum SensitivityComputation {ANALYTICAL = 0, SEMIANALYTICAL = 1,  FINITEDIFFERENCE = 2} scFlag;
  enum LsSolver {QR=0, NORMAL_EQUATIONS=1} lsSolver;
  enum Mvp {FD = 0, H1 = 1, H2 = 2, H1FD = 3} mvp;
  enum Compatible3D {OFF_COMPATIBLE3D = 0, ON_COMPATIBLE3D = 1} comp3d;
  enum AngleRadians {OFF_ANGLERAD = 0, ON_ANGLERAD = 1} angleRad;

  enum SensitivityFSI  {OFF_SENSITIVITYFSI  = 0, ON_SENSITIVITYFSI  = 1} sensFSI;
  enum SensitivityMesh {OFF_SENSITIVITYMESH = 0, ON_SENSITIVITYMESH = 1} sensMesh;
  enum SensitivityMach {OFF_SENSITIVITYMACH = 0, ON_SENSITIVITYMACH = 1} sensMach;
  enum SensitivityAOA {OFF_SENSITIVITYALPHA = 0, ON_SENSITIVITYALPHA = 1} sensAlpha;
  enum SensitivityYAW {OFF_SENSITIVITYBETA = 0, ON_SENSITIVITYBETA = 1} sensBeta;
  enum SensitivityLiftx {OFF_SENSITIVITYLIFTX = 0, ON_SENSITIVITYLIFTX = 1} sensLiftx;
  enum SensitivityLifty {OFF_SENSITIVITYLIFTY = 0, ON_SENSITIVITYLIFTY = 1} sensLifty;
  enum SensitivityLiftz {OFF_SENSITIVITYLIFTZ = 0, ON_SENSITIVITYLIFTZ = 1} sensLiftz;

  // This flag repeats the linear solves until the number of iterations
  // is smaller than the maximum allowed.
  // Default Value = OFF_EXACTSOLUTION
  enum ExactSolution {OFF_EXACTSOLUTION = 0, ON_EXACTSOLUTION = 1} excsol;

  enum HomotopyComputation {OFF_HOMOTOPY = 0, ON_HOMOTOPY = 1} homotopy;
  enum FixSolution {NONEFIX = 0, PREVIOUSVALEUSFIX = 1} fixsol;
  enum AdaptiveEpsFSI {OFF_ADAPTIVEEPSFSI = 0, ON_ADAPTIVEEPSFSI = 1} adaptiveEpsFSI;

  double machref;
  double alpharef;
  double betaref;

  const char* meshderiv;
  const char* sensoutput;

  //TODO temporary parameters for debugging
  const char* tempStateDeriv;
  const char* linsolverhs;
  const char* dFdS_final;

  bool densFlag;
  bool pressFlag;
  bool apressFlag;
  bool fsiFlag;

  int sparseFlag;
  int numShapeVariables;
  int avgsIt;

  double eps;
  double fres;

  KspData ksp;

  SensitivityAnalysis();
  ~SensitivityAnalysis() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SymmetryData {

  double nx, ny, nz;

  SymmetryData();
  ~SymmetryData() {};

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct BLMeshMotionData {
   int numLayers;
   int numIncrements;
   double power;
   double neiSelectionDistFactor;

   enum Type {PSEUDOSTRUCTURAL = 0, ALGEBRAIC = 1 } type;
   enum FractionalStrategy {Distance = 1, DotProduct = 2} fractionalStrategy;
   enum ClosestNeighbor {Fixed = 1, Variable = 0} bestNeiStrategy;
   int feedbackFrequency;

   BLMeshMotionData();
   ~BLMeshMotionData() {};

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct DefoMeshMotionData {

  enum Type {BASIC = 0, COROTATIONAL = 1} type;
  enum Element {LINEAR_FE = 0, NON_LINEAR_FE = 1, TORSIONAL_SPRINGS = 2, BALL_VERTEX = 3, NL_BALL_VERTEX = 4 } element;

  double volStiff;
  enum Mode {Recursive = 1, NonRecursive = 2} mode;
  int numIncrements;
  enum SlidingSurfaceTreatment {Default = 0, PrescribedAverage = 1} slidingSurfaceTreatment;

  BLMeshMotionData blmeshmotion;
  NewtonData<KspData> newton;
  SymmetryData symmetry;

  DefoMeshMotionData();
  ~DefoMeshMotionData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct VelocityPoints {

  double time;
  double velocityX;
  double velocityY;
  double velocityZ;

  VelocityPoints();
  ~VelocityPoints() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ForcePoints {

  double time;
  double force;

  ForcePoints();
  ~ForcePoints() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct RigidMeshMotionData {

  static const int num = 10;

  enum Tag {MACH = 0, TIME = 1, VELOCITY = 2} tag;
  enum LawType {VELOCITYPOINTS = 0, CONSTANTACCELERATION = 1} lawtype;

  double vx;
  double vy;
  double vz;

  double ax;
  double ay;
  double az;

  VelocityPoints* vpts[num];
  VelocityPoints vpts1;
  VelocityPoints vpts2;
  VelocityPoints vpts3;
  VelocityPoints vpts4;
  VelocityPoints vpts5;
  VelocityPoints vpts6;
  VelocityPoints vpts7;
  VelocityPoints vpts8;
  VelocityPoints vpts9;
  VelocityPoints vpts10;

  double timestep;

  RigidMeshMotionData();
  ~RigidMeshMotionData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct AeroelasticData {

  enum Force {LAST = 0, AVERAGED = 1, LAST_KRIS = 2} force;
  double pressure;
  double displacementScaling;
  double forceScaling;
  double powerScaling;

  AeroelasticData();
  ~AeroelasticData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct HeavingData {

  enum Domain {VOLUME = 0, SURFACE = 1} domain;

  double ax;
  double ay;
  double az;

  HeavingData();
  ~HeavingData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//----------------------------------------------------------

struct SpiralingData {

  enum Domain {VOLUME = 0, SURFACE = 1} domain;

  double xL;
  double x0;

  SpiralingData();
  ~SpiralingData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//----------------------------------------------------------
struct PitchingData {

  enum Domain {VOLUME = 0, SURFACE = 1} domain;

  double alpha_in;
  double alpha_max;
  double alpha_slope;
  double x11;
  double y11;
  double z11;
  double x21;
  double y21;
  double z21;

  double beta_in;
  double beta_max;
  double beta_slope;
  double x12;
  double y12;
  double z12;
  double x22;
  double y22;
  double z22;
  PitchingData();
  ~PitchingData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//----------------------------------------------------------

struct DeformingData {

  enum Domain {VOLUME = 0, SURFACE = 1} domain;

  const char *positions;

  double amplification;

  DeformingData();
  ~DeformingData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//-----------------------------------------------------------------------------

struct RotationData  {

  double nx, ny, nz;
  double x0, y0, z0;
  double omega;
  enum InfRadius {FALSE = 0, TRUE = 1} infRadius;

  RotationData();
  Assigner *getAssigner();
  void setup(const char *, ClassAssigner * = 0);

};

//-----------------------------------------------------------------------------

struct Velocity  {

  ObjectMap<RotationData> rotationMap;

  void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct ForcedData {

  enum Type {HEAVING = 0, PITCHING = 1, VELOCITY = 2, DEFORMING = 3, DEBUGDEFORMING=4,
             ACOUSTICBEAM=5, SPIRALING = 6, ACOUSTICVISCOUSBEAM=7} type;

  double frequency;
  double timestep;

  HeavingData hv;
  SpiralingData sp;
  PitchingData pt;
  Velocity vel;
  DeformingData df;

  double tsoffset;

  ForcedData();
  ~ForcedData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//----------------------------------------------------------

struct PadeData {

  static const int num = 11;

  double freq[num];
  double freq1;
  double freq2;
  double freq3;
  double freq4;
  double freq5;
  double freq6;
  double freq7;
  double freq8;
  double freq9;
  double freq10;
  double freq11;

  int nPoints;
  int degNum;
  int degDen;

  PadeData();
  ~PadeData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct NonlinearRomDirectoriesData {

  const char *prefix;
  const char *databaseName;
  const char *clusterName;
  const char *sensitivityClusterName;

  NonlinearRomDirectoriesData();
  ~NonlinearRomDirectoriesData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct NonlinearRomFilesData {

  // A nonlinear ROM database can consist of a very large number of files.
  // To make life easier for the user, this code introduces a "prefix" feature,
  // which tells aero-f to read and write database files using a built-in
  // naming convention.

  // If a prefix and a name are both given, the name overrides the prefix.

  enum DuplicateSnapshots {DUPLICATE_SNAPSHOTS_FALSE = 0, DUPLICATE_SNAPSHOTS_TRUE = 1} duplicateSnapshots;

  // State snapshot clusters
  const char *statePrefix;
  const char *stateSnapsName;
  const char *mapName;
  const char *indexName;
  const char *connName;
  const char *centersName;
  const char *nearestName;
  const char *centerNormsName;
  const char *distanceMatrixName;

  // State bases
  const char *stateBasisPrefix;
  const char *stateBasisName;
  const char *stateSingValsName;
  const char *simpleUpdateInfoName;
  const char *exactUpdateInfoPrefix;
  const char *stateDistanceComparisonInfoName;
  const char *stateDistanceComparisonInfoExactUpdatesName;
  const char *stateDistanceComparisonInfoExactUpdatesMultiICName;
  const char *basisNormalizedCenterProductsName;
  const char *refStateName;
  const char *projErrorName;

  // Krylov snaps
  const char *krylovPrefix;
  const char *krylovSnapsName;

  // Krylov bases
  const char *krylovBasisPrefix;
  const char *krylovBasisName;
  const char *krylovSingValsName;
  const char *krylovDistanceComparisonInfoName;

  // Sensitivities
  const char *sensitivityPrefix;
  const char *sensitivitySnapsName;

  // Sensitivity basis
  const char *sensitivityBasisPrefix;
  const char *sensitivityBasisName;
  const char *sensitivitySingValsName;
  const char *sensitivityDistanceComparisonInfoName;

  // Residual snaps
  const char *residualPrefix;
  const char *residualSnapsName;

  // Residual bases
  const char *residualBasisPrefix;
  const char *residualBasisName;
  const char *residualSingValsName;

  // Action-of-Jacobian snaps
  const char *jacActionPrefix;
  const char *jacActionSnapsName;

  // Action-of-Jacobian bases
  const char *jacActionBasisPrefix;
  const char *jacActionBasisName;
  const char *jacActionSingValsName;

  // Gappy quantities
  const char *gappyPrefix;
  const char *sampledNodesName;         //sampleNodes;
  const char *sampledNodesFullCoordsName; // sampled nodes in full mesh coordinates
  const char *sampledCentersName;
  const char *sampledStateBasisName;    //podStateRed;
  const char *sampledKrylovBasisName;
  const char *sampledSensitivityBasisName;
  const char *sampledResidualBasisName; //podFileResHat;
  const char *sampledJacActionBasisName; //podFileJacHat;
  const char *sampledMeshName;          //mesh;
  const char *sampledSolutionName;      //solution;
  const char *sampledMatchStateName;      //comparison state;
  const char *sampledShapeDerivativeName;
  const char *sampledMultiSolutionsName; // multiple solutions. Can start from one, or an arbitrary linear combination.
  const char *sampledRefStateName;
  const char *sampledWallDistName;      //wallDistanceRed;
  const char *sampledDisplacementName;  // sampled initial displacement vector
  const char *gappyJacActionName;             //jacMatrix in sampled coords;
  const char *gappyResidualName;             //resMatrix in sampled coords;
  const char *approxMetricStateLowRankName; // approximated metric in reduced mesh coordinates
  const char *approxMetricNonlinearLowRankName;
  const char *approxMetricStateLowRankFullCoordsName; // approximated metric in full mesh coordinates
  const char *approxMetricNonlinearLowRankFullCoordsName;
  const char *approxMetricStateLowRankSurfaceCoordsName;
  const char *approxMetricNonlinearName; // ascii file
  const char *correlationMatrixName;
  const char *sampledApproxMetricNonlinearSnapsName;

  // Surface quantities
  const char *surfacePrefix;
  const char *surfaceCentersName;
  const char *surfaceStateBasisName;
  const char *surfaceRefStateName;
  const char *surfaceSolutionName;
  const char *surfaceMatchStateName;
  const char *surfaceInitialDisplacementName;
  const char *surfaceShapeDerivativeName;
  const char *surfaceMultiSolutionsName;
  const char *surfaceWallDistName;
  const char *surfaceDisplacementName;
  const char *surfaceMeshName;

  NonlinearRomFilesData();
  ~NonlinearRomFilesData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct NonlinearRomFileSystemData {

  enum AvgIncrementalStates {AVG_INCREMENTAL_STATES_FALSE = 0, AVG_INCREMENTAL_STATES_TRUE = 1} avgIncrementalStates;
  enum DistanceMetric {DIST_EUCLIDEAN = 0, DIST_ANGLE = 1 } distanceMetric;

  int nClusters;

  NonlinearRomDirectoriesData directories;
  NonlinearRomFilesData files;

  NonlinearRomFileSystemData();
  ~NonlinearRomFileSystemData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct NonlinearRomOnlineNonStateData {

  enum Include {INCLUDE_OFF = 0, INCLUDE_ON = 1} include;
  enum GramSchmidt {GRAMSCHMIDT_OFF = 0, GRAMSCHMIDT_ON = 1} gramSchmidt;

	int minDimension;
  int maxDimension;
  double energy;

  int timeFreq;
  int newtonFreq;

  NonlinearRomOnlineNonStateData();
  ~NonlinearRomOnlineNonStateData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct NonlinearRomOnlineData {

  enum Projection {PETROV_GALERKIN = 0, GALERKIN = 1} projection;
  enum SystemApproximation {SYSTEM_APPROXIMATION_NONE = 0, GNAT = 1, COLLOCATION = 2, APPROX_METRIC_NL = 3} systemApproximation;
  enum LineSearch {LINE_SEARCH_FALSE = 0, LINE_SEARCH_BACKTRACKING = 1, LINE_SEARCH_WOLF = 2} lineSearch;
  enum LSSolver {QR = 0, NORMAL_EQUATIONS = 1, LEVENBERG_MARQUARDT_SVD = 2, PROBABILISTIC_SVD = 3} lsSolver;

  enum ResidualScaling {SCALING_OFF=0, SCALING_BALANCED=1, SCALING_ENERGY=2} residualScaling;
  double turbulenceWeight;
  double eddyLengthScale;

  double levenbergMarquardtWeight;

  int minDimension;
  int maxDimension;
  double energy;
  double bufferEnergy;

  int randMatDimension;

  double incrementCoordsTol;

  enum BasisUpdates {UPDATES_OFF = 0, UPDATES_SIMPLE = 1, UPDATES_FAST_EXACT = 2, UPDATES_FAST_APPROX = 3} basisUpdates;
  int basisUpdateFreq;
  int tryAllFreq;
  double basisUpdateTolerance;

  enum ProjectSwitchStateOntoAffineSubspace {PROJECT_OFF = 0, PROJECT_ON = 1} projectSwitchStateOntoAffineSubspace;

  enum DistanceComparisons {DISTANCE_COMPARISONS_OFF = 0, DISTANCE_COMPARISONS_ON = 1} distanceComparisons;
  enum StoreAllClusters {STORE_ALL_CLUSTERS_FALSE = 0, STORE_ALL_CLUSTERS_TRUE = 1} storeAllClusters;

  double romSpatialOnlyInitialHomotomyStep;
  double romSpatialOnlyMaxHomotomyStep;
  double romSpatialOnlyHomotomyStepExpGrowthRate;

  double newtonStepThreshold;

  enum MeritFunction {ROM_RESIDUAL=0, HDM_RESIDUAL=1} meritFunction;

  NonlinearRomOnlineNonStateData krylov;
  NonlinearRomOnlineNonStateData sensitivity;

  NonlinearRomOnlineData();
  ~NonlinearRomOnlineData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct RelativeProjectionErrorData {

  enum RelativeProjectionError {REL_PROJ_ERROR_OFF = 0, REL_PROJ_ERROR_STATE = 1, REL_PROJ_ERROR_RESIDUAL = 2, REL_PROJ_ERROR_JACACTION = 3} relProjError;
  enum ProjectIncrementalSnapshots {PROJECT_INCREMENTAL_SNAPS_FALSE = 0, PROJECT_INCREMENTAL_SNAPS_TRUE = 1} projectIncrementalSnaps;
  enum ProjectSnapshotsMinusRefSol {PROJECT_SNAPS_MINUS_REF_SOL_FALSE = 0, PROJECT_SNAPS_MINUS_REF_SOL_TRUE = 1} subtractRefSol;

  int minDimension;
  int maxDimension;
  double energy;

  int sweepFreq;
  enum UseFirstStateAsRefStateForIncrBasis {ASSUME_INCR_REFSTATE_FALSE = 0, ASSUME_INCR_REFSTATE_TRUE = 1} useFirstStateAsRefStateForIncrBasis;

  enum BasisUpdates {UPDATES_OFF = 0, UPDATES_SIMPLE = 1} basisUpdates;

  enum PostProProjectedStates {POST_PRO_OFF = 0, POST_PRO_ON = 1} postProProjectedStates;

  NonlinearRomOnlineNonStateData krylov;
  NonlinearRomOnlineNonStateData sensitivity;

  RelativeProjectionErrorData();
  ~RelativeProjectionErrorData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct StateSnapshotsData {

  enum NormalizeSnaps {NORMALIZE_FALSE = 0, NORMALIZE_TRUE = 1} normalizeSnaps;
  enum SubtractClusterCenters {SUBTRACT_CENTERS_FALSE = 0, SUBTRACT_CENTERS_TRUE = 1} subtractCenters;
  enum SubtractNearestSnapshotToCenter {SUBTRACT_NEAREST_FALSE = 0, SUBTRACT_NEAREST_TRUE = 1} subtractNearestSnapsToCenters;
  enum SubtractRefState {SUBTRACT_REF_STATE_FALSE = 0, SUBTRACT_REF_STATE_TRUE = 1} subtractRefState;

  StateSnapshotsData();
  ~StateSnapshotsData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SnapshotsData {

  enum NormalizeSnaps {NORMALIZE_FALSE = 0, NORMALIZE_TRUE = 1} normalizeSnaps;

  SnapshotsData();
  ~SnapshotsData() {}

  void setup(const char *, ClassAssigner * = 0);

};


//------------------------------------------------------------------------------

struct DataCompressionData {

  enum ComputePOD {COMPUTE_POD_FALSE = 0, COMPUTE_POD_TRUE = 1} computePOD;
  enum Type {POD = 0, BALANCED_POD = 1} type;
  enum PODMethod {SCALAPACK_SVD = 0, PROBABILISTIC_SVD = 1, R_SVD = 2,  Eig = 3} podMethod;
  int randMatDimension;
  int nPowerIts;
  enum CompareSVDMethods {COMPARE_SVD_FALSE = 0, COMPARE_SVD_TRUE = 1} compareSVDMethods;
  enum TestProbabilisticSVD {TEST_PROBABILISTIC_SVD_FALSE = 0, TEST_PROBABILISTIC_SVD_TRUE = 1} testProbabilisticSVD;
  int maxVecStorage;
  enum EnergyOnly {ENERGY_ONLY_FALSE = 0, ENERGY_ONLY_TRUE = 1} energyOnly;
  double tolerance;
  int maxBasisSize;
  int minBasisSize;
  double singValTolerance;
  double maxEnergyRetained;
  int initialCluster;

  DataCompressionData();
  ~DataCompressionData() {}

  void setup(const char *, ClassAssigner * = 0);

};


//------------------------------------------------------------------------------

struct StateData {

  StateSnapshotsData snapshots;
  DataCompressionData dataCompression;

  StateData();
  ~StateData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ResidualData {

  SnapshotsData snapshots;
  DataCompressionData dataCompression;

  ResidualData();
  ~ResidualData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct JacobianActionData {

  SnapshotsData snapshots;
  DataCompressionData dataCompression;

  JacobianActionData();
  ~JacobianActionData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SensitivityData {

  SnapshotsData snapshots;
  DataCompressionData dataCompression;

  SensitivityData();
  ~SensitivityData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct KrylovData {

  SnapshotsData snapshots;
  DataCompressionData dataCompression;

  KrylovData();
  ~KrylovData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct ClusteringData {

  enum ClusteringAlgorithm {K_MEANS_WITHOUT_BOUNDS = 0, K_MEANS_WITH_BOUNDS = 1} clusteringAlgorithm;
  enum KmeansBoundType {TIGHT_BOUNDS = 0, LOOSE_BOUNDS = 1} kmeansBoundType;
  double percentOverlap;
  int maxIter;
  int maxIterAggressive;
  int minClusterSize;
  double kMeansTol;
  int kMeansRandSeed;
  enum UseExistingClusters {USE_EXISTING_CLUSTERS_FALSE = 0, USE_EXISTING_CLUSTERS_TRUE = 1} useExistingClusters;
  enum ComputeMDS {COMPUTE_MDS_FALSE = 0, COMPUTE_MDS_TRUE = 1} computeMDS;
  enum ClusterFilesSeparately {CLUSTER_FILES_SEPARATELY_FALSE = 0, CLUSTER_FILES_SEPARATELY_TRUE = 1} clusterFilesSeparately;
  double snapshotNormTolerance;

  ClusteringData();
  ~ClusteringData() {}

  void setup(const char *, ClassAssigner * = 0);

};


//------------------------------------------------------------------------------

struct ApproximatedMetricData {

  double sampledMeshFraction;
  double lowRankEnergy;
  double tolerance;
  ApproximatedMetricData();
  ~ApproximatedMetricData() {}

  void setup(const char *, ClassAssigner * = 0);

};
//------------------------------------------------------------------------------

struct BasisUpdatesData {

  enum PreprocessForNoUpdates {NO_UPDATES_FALSE = 0, NO_UPDATES_TRUE = 1} preprocessForNoUpdates;
  enum PreprocessForProjections {PROJECTIONS_FALSE = 0, PROJECTIONS_TRUE = 1} preprocessForProjections;
  enum PreprocessForSimpleUpdates {SIMPLE_UPDATES_FALSE = 0, SIMPLE_UPDATES_TRUE = 1} preprocessForSimpleUpdates;
  enum PreprocessForExactUpdates {EXACT_UPDATES_FALSE = 0, EXACT_UPDATES_TRUE = 1} preprocessForExactUpdates;
  enum PreprocessForApproxUpdates {APPROX_UPDATES_FALSE = 0, APPROX_UPDATES_TRUE = 1} preprocessForApproxUpdates;
  ApproximatedMetricData approxMetricState;


  BasisUpdatesData();
  ~BasisUpdatesData() {}

  void setup(const char *, ClassAssigner * = 0);

};

struct EmbeddedAlternatingLeastSquareData {
    int maxBasisSize;
    double relativeMinimumEnergy;
    int maxIteration;
    enum LeastSquareSolver {QR = 0, SVD = 1} leastSquareSolver;
    StateSnapshotsData stateSnapshotsData;

    EmbeddedAlternatingLeastSquareData();
    ~EmbeddedAlternatingLeastSquareData() {}

    void setup(const char *, ClassAssigner * = 0);
};

//------------------------------------------------------------------------------

struct ROBConstructionData {

  StateData state;
  ResidualData residual;
  JacobianActionData jacAction;
  SensitivityData sensitivity;
  KrylovData krylov;

  ClusteringData clustering;
  BasisUpdatesData basisUpdates;
  RelativeProjectionErrorData relativeProjectionError;

    EmbeddedAlternatingLeastSquareData embeddedALS;

  enum StoreAllClusters {STORE_ALL_CLUSTERS_FALSE = 0, STORE_ALL_CLUSTERS_TRUE = 1} storeAllClusters;

  ROBConstructionData();
  ~ROBConstructionData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct GappyConstructionData {

  enum DoPrepro {DO_PREPRO_FALSE = 0, DO_PREPRO_TRUE = 1} doPrepro;
  enum SowerInputs {SOWER_INPUTS_FALSE = 0, SOWER_INPUTS_TRUE = 1} sowerInputs;
  enum DoPreproGNAT {DO_PREPRO_GNAT_FALSE = 0, DO_PREPRO_GNAT_TRUE = 1} doPreproGNAT;
  enum DoPreproApproxMetricNonlinear {DO_PREPRO_APPROX_METRIC_NL_FALSE = 0, DO_PREPRO_APPROX_METRIC_NL_TRUE = 1} doPreproApproxMetricNonlinear;
  enum DoPreproApproxMetricNonlinearNNLS {DO_PREPRO_NNLS_METRIC_NL_FALSE = 0, DO_PREPRO_NNLS_METRIC_NL_TRUE = 1} doPreproApproxMetricNonlinearNNLS;

  int maxDimensionState;
  int minDimensionState;
  double energyState;

  int maxDimensionSensitivity;
  int minDimensionSensitivity;
  double energySensitivity;

  int maxDimensionKrylov;
  int minDimensionKrylov;
  double energyKrylov;

  int maxDimensionResidual;
  int minDimensionResidual;
  double energyResidual;

  int maxDimensionJacAction;
  int minDimensionJacAction;
  double energyJacAction;

  enum SelectSampledNodes {SELECT_SAMPLED_NODES_FALSE = 0, SELECT_SAMPLED_NODES_TRUE = 1} selectSampledNodes;

  enum greedyData {UNSPECIFIED_GREEDY = -1, STATE_ROB_GREEDY = 0, RESIDUAL_ROB_GREEDY = 1,
                   JACOBIAN_ROB_GREEDY = 2, RESIDUAL_AND_JACOBIAN_ROBS_GREEDY = 3, SPECIFIED_SNAPS_GREEDY = 4} greedyData;
  enum GreedyLeastSquaresSolver {GREEDY_LS_PROBABILISTIC = 0, GREEDY_LS_SCALAPACK = 1, GREEDY_LS_LINPACK = 2} greedyLeastSquaresSolver;

  enum PseudoInverseSolver {PSEUDO_INVERSE_SCALAPACK = 0, PSEUDO_INVERSE_LINPACK = 1} pseudoInverseSolver;
  int pseudoInverseNodes;

  double initialCluster; // for online matrix computations (restart)

  int maxClusteredSnapshotsNonlinearApproxMetric;

  int randMatDimension;
  int nPowerIts;

  int maxDimGreedyAlgorithm;
  int minDimGreedyAlgorithm;
  double dimGreedyAlgorithmFactor;

  int maxSampledNodes;
  int minSampledNodes;
  double sampledNodesFactor;
  int layers;

  enum IncludeLiftFaces {NONE_LIFTFACE = 0,
		SPECIFIED_LIFTFACE  = 1, ALL_LIFTFACE = 2} includeLiftFaces;

  double minFractionOfSampledNodesOnSurfaceInTargetRegion;
  double minFractionOfSampledNodesInTargetRegion;
  SchemeFixData sampledMeshTargetRegion;  // use fix regions to specify target areas for sampled mesh construction

  ApproximatedMetricData approxMetricNonlinear;

  enum ComputeGappyRes {NO_GAPPYRES = 0, YES_GAPPYRES  = 1} computeGappyRes;

  enum UseUnionOfSampledNodes {UNION_FALSE = 0, UNION_TRUE = 1} useUnionOfSampledNodes;

  enum UseOldReducedSVecFunction {USE_OLD_FALSE = 0, USE_OLD_TRUE = 1} useOldReducedSVecFunction;

  enum SampledMeshUsed {SAMPLED_MESH_NOT_USED = 0, SAMPLED_MESH_USED = 1} sampledMeshUsed;

  enum OutputReducedBases {OUTPUT_REDUCED_BASES_FALSE = 0, OUTPUT_REDUCED_BASES_TRUE = 1} outputReducedBases;
  enum TestApproxMetric {TEST_APPROX_METRIC_FALSE = 0, TEST_APPROX_METRIC_TRUE = 1} testApproxMetric;

  GappyConstructionData();
  ~GappyConstructionData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct NonlinearRomOfflineData {

  ROBConstructionData rob;
  GappyConstructionData gappy;

  NonlinearRomOfflineData();
  ~NonlinearRomOfflineData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct LinearizedData {

  enum PadeReconstruction {TRUE = 1, FALSE = 0} padeReconst;
  // Type is a defunct variable - it is not documented and only used for development and testing
  enum Type {DEFAULT = 0, ROM = 1, FORCED = 2} type;

  enum Domain {TIME = 0, FREQUENCY = 1} domain;
  enum InitialCondition {DISPLACEMENT = 0, VELOCITY = 1} initCond;
  enum GramSchmidt {TRUE_GS = 1, FALSE_GS = 0} doGramSchmidt;
  enum ErrorIndicator {OIBEI = 0, RBEI1 = 1, RBEI2 = 2, RBEI3 = 3, RBEI4 = 4} errorIndicator;
  double amplification;
  double frequency;
  double stepsize;
  double stepsizeinitial;
  double freqStep;
  double eps;
  double eps2;
  double epsEV;
  double tolerance;
  double refLength;
  const char *strModesFile;
  int modeNumber;
  int numSteps;
  int numPOD;
  int numStrModes;
  int maxItEV;
  const char *romFile;
  static const int numFreq = 20;
  double gamFreq[numFreq];
  double gamFreq1;
  double gamFreq2;
  double gamFreq3;
  double gamFreq4;
  double gamFreq5;
  double gamFreq6;
  double gamFreq7;
  double gamFreq8;
  double gamFreq9;
  double gamFreq10;
  double gamFreq11;
  double gamFreq12;
  double gamFreq13;
  double gamFreq14;
  double gamFreq15;
  double gamFreq16;
  double gamFreq17;
  double gamFreq18;
  double gamFreq19;
  double gamFreq20;



  PadeData pade;

  LinearizedData();
  ~LinearizedData() {}

  void setup(const char *, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

struct SurfaceData  {

  double nx, ny, nz;
  int sBit;
  static const int UNSPECIFIED = -1;
  enum ComputeForces {FALSE = 0, TRUE = 1 } computeForces;
  enum ForceResults {NO = 0, YES = 1} forceResults;
  int rotationID;
  int forceID;
  int bcID;
  double velocity;

  enum Type { ADIABATIC = 1, ISOTHERMAL = 2 } type;
  double temp;

  enum ComputeHeatPower {FALSE_HF = 0, TRUE_HF = 1 } computeHeatFluxes;
  enum HeatFluxResults {UNSPECIFIED_HF = -1, NO_HF = 0, YES_HF = 1} heatFluxResults;
  //the HF (Heat Flux) index ensures that there is no confusion with the force related data.

  SurfaceData();
  Assigner *getAssigner();
  void setBit(int b) { sBit = b; }
};

//------------------------------------------------------------------------------

struct Surfaces  {

  ObjectMap<SurfaceData> surfaceMap;
  void setup(const char *);
};

//------------------------------------------------------------------------------

struct PadeFreq  {

  ObjectMap<double> freqMap;
  void setup(const char *);
};

//------------------------------------------------------------------------------

struct EmbeddedFramework {

  enum IntersectorName {PHYSBAM = 0, FRG = 1} intersectorName;
  enum StructureNormal {ELEMENT_BASED = 0, NODE_BASED = 1} structNormal;
  enum EOSChange {NODAL_STATE = 0, RIEMANN_SOLUTION = 1} eosChange;
  enum ForceAlgorithm {RECONSTRUCTED_SURFACE = 0, CONTROL_VOLUME_BOUNDARY = 1, EMBEDDED_SURFACE = 2} forceAlg;
  enum RiemannNormal {STRUCTURE = 0, FLUID = 1} riemannNormal;
  enum PhaseChangeAlgorithm {AVERAGE = 0, LEAST_SQUARES = 1, AUTO = 2} phaseChangeAlg;
  enum InterfaceAlgorithm {MID_EDGE = 0, INTERSECTION = 1} interfaceAlg;

  enum InterfaceLimiter {LIMITERNONE = 0, LIMITERALEX1 = 1} interfaceLimiter;
  // Low mach preconditioning of the exact Riemann problem.
  // Added by Alex Main (December 2013)
  //
  enum Prec {NON_PRECONDITIONED = 0, PRECONDITIONED = 1, SAME_AS_PROBLEM = 3 } prec;

  double alpha;   // In the case of solve Riemann problem at intersection, this parameter
                  // controls whether to switch to a first order method to avoid divided-by-zero

  // stabilizing alpha (attempt at stabilizing the structure normal)
  // Tries to add some dissipation.  should be small.
  double stabil_alpha;

  double interfaceThickness;

  MultiInitialConditionsData embedIC;

  int nLevelset; //number of level-sets. Currently only consider bubbles.

  int qOrder; // order of quadrature rule used for EMBEDDED_SURFACE forceAlg

  enum CrackingWithLevelSet {OFF = 0, ON = 1} crackingWithLevelset;
  enum Reconstruction {CONSTANT = 0, LINEAR = 1} reconstruct;
  enum ViscousInterfaceOrder {FIRST = 0, SECOND = 1} viscousinterfaceorder;
  enum ViscousBoundaryCondition {WEAK = 0, STRONG = 1} viscousboundarycondition;
  enum SurrogateInterface{HYBRID = 0, EXTERNAL = 1} surrogateinterface;

  int testCase;

  EmbeddedFramework();
  ~EmbeddedFramework() {}

  void setup(const char *);
};

//------------------------------------------------------------------------------

struct OneDimensionalInfo {
  enum CoordinateType {CARTESIAN = 0, CYLINDRICAL = 1, SPHERICAL = 2} coordType;
  enum VolumeType { CONSTANT_VOLUME = 0, REAL_VOLUME = 1} volumeType;
  double maxDistance; //mesh goes from 0 to maxDistance

  int numPoints; //mesh has numPoints elements
  int fluidId2;

  int sourceTermOrder;

  double interfacePosition;

  double density1, velocity1, pressure1,temperature1;
  double density2, velocity2, pressure2,temperature2;

  ProgrammedBurnData programmedBurn;

  enum Mode { NORMAL=0, CONVTEST1 = 1, CONVTEST2=2 } mode;

  enum ProblemMode { MULTIFLUID=0, FSI=1} problemMode;

  OneDimensionalInfo();
  ~OneDimensionalInfo() {}

  void setup(const char *);
};
//------------------------------------------------------------------------------

struct ImplosionSetup {
  enum Type{LINEAR=0, SMOOTHSTEP=1} type;
  double Prate, Pinit, tmax;
  int intersector_freq;

  ImplosionSetup();
  ~ImplosionSetup() {}
  void setup(const char *);
};


//------------------------------------------------------------------------------

struct KirchhoffData {

  /// UH (08/2012)
  ///
  /// This structure stores information for computing the Kirchhoff integral.
  /// Information is used with the problem type "Aeroacoustic".
  ///

  enum Type {CYLINDRICAL = 0, SPHERICAL = 1} d_surfaceType;
  double d_energyFraction;
  int d_angularIncrement;
  int d_nyquist;

  KirchhoffData();
  ~KirchhoffData() {}

  void setup(Communicator *communicator, const char *name, ClassAssigner * = 0);

};

//------------------------------------------------------------------------------

class IoData {

  char *cmdFileName;
  FILE *cmdFilePtr;

  Communicator *com;

public:

  InputData input;
  OutputData output;
  PreconditionData prec;
  RestartParametersData restart;
  ProblemData problem;
  ReferenceStateData ref;
  BcsData bc;
  EquationsData eqs;
  MultiFluidData mf;
  SchemesData schemes;

// Included (MB)
  SensitivityAnalysis sa;

  TsData ts;
  DGCLData dgcl;
  DefoMeshMotionData dmesh;
  RigidMeshMotionData rmesh;
  AeroelasticData aero;
  ForcedData forced;
  NonlinearRomOfflineData romOffline;
  NonlinearRomOnlineData romOnline;
  NonlinearRomFileSystemData romDatabase;
  LinearizedData linearizedData;
  Surfaces surfaces;
  Velocity rotations;
  Volumes volumes;
  EmbeddedFramework embed;
  OneDimensionalInfo oneDimensionalInfo;
  ImplosionSetup implosion;

  MultiGridData mg;

  // UH (08/2012)
  // The next member is used for the Kirchhoff integral.
  KirchhoffData surfKI;

public:

  IoData(Communicator *);
  ~IoData() {}

  void readCmdLine(int, char**);
  void setupCmdFileVariables();
  void readCmdFile();
  void resetInputValues();
  int checkFileNames();
  int checkInputValues();
  int checkInputValuesAeroAcoustic();
  int checkInputValuesAllEquationsOfState();
  int checkInputValuesProgrammedBurn();
  int checkProgrammedBurnLocal(ProgrammedBurnData& programmedBurn,
                               InitialConditions& IC);
  int checkCFLBackwardsCompatibility();
  int checkInputValuesAllInitialConditions();
  void nonDimensionalizeAllEquationsOfState();
  void nonDimensionalizeAllInitialConditions();
  void nonDimensionalizeForcedMotion();
  void nonDimensionalizeOneDimensionalProblem();
  int checkInputValuesNonDimensional();
  int checkInputValuesDimensional(map<int,SurfaceData*>& surfaceMap);
  int checkInputValuesEssentialBC();
  void checkInputValuesTurbulence();
  void checkInputValuesDefaultOutlet();
  int checkBoundaryValues();
  int checkSolverValues(map<int,SurfaceData*>& surfaceMap);
  int checkInputValuesInitialConditions(InitialConditions &initialConditions,
                                        int fluidModelID);
  int checkInputValuesEquationOfState(FluidModelData &fluidModel, int fluidModelID);
  void nonDimensionalizeInitialConditions(InitialConditions &initialConditions);
  void nonDimensionalizeFluidModel(FluidModelData &fluidModel);
  void nonDimensionalizeViscosityModel(ViscosityModelData &vm);
  void nonDimensionalizeThermalCondModel(ThermalCondModelData &tm);
  int checkInputValuesSparseGrid(SparseGridData &sparseGrid);
  int checkInputValuesEmbeddedFramework();
  void printDebug();

  void setupOneDimensional();

  int checkInputValuesNonlinearRomPreprocessing();
  int checkInputValuesNonlinearRomOnline();
  int checkInputValuesNonlinearRomPostprocessing();

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <IoData.C>
#endif

#endif
