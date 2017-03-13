#include <cstdlib>
#include <cstring>

#include <TsOutput.h>

#include <IoData.h>
#include <RefVal.h>
#include <Domain.h>
#include <PostOperator.h>
#include <MeshMotionHandler.h>
#include <DistVector.h>
#include <DistExactRiemannSolver.h>
#include <BinFileHandler.h>
#include <VectorSet.h>
#include <GhostPoint.h>
#include <SubDomain.h>

//------------------------------------------------------------------------------

template<int dim>
TsOutput<dim>::TsOutput(IoData &iod, RefVal *rv, Domain *dom, PostOperator<dim> *po) :
refVal(rv),
domain(dom),
postOp(po),
rmmh(0),
hmmh(NULL),
smmh(NULL),
pmmh(NULL),
dmmh(NULL),
com(NULL),
steady(false),
it0(0),
frequency(-1),
frequency_dt(-1.0),
prtout(-1.0),
numFluidPhases(-1),
length(-1.0),
surface(-1.0),
stateOutputFreqTime(-1),
stateOutputFreqNewton(-1),
residualOutputFreqTime(-1),
residualOutputMaxNewton(-1),
fdResiduals(false),
fdResidualsLimit(false),
externalSI(false),
fullOptPressureName(NULL),
optPressureDimensional(false),
forces(NULL),
tavforces(NULL),
hydrostaticforces(NULL),
hydrodynamicforces(NULL),
lift(NULL),
matchpressure(NULL),
matchstate(NULL),
fluxnorm(NULL),
tavlift(NULL),
hydrostaticlift(NULL),
hydrodynamiclift(NULL),
generalizedforces(NULL),
residuals(NULL),
material_volumes(NULL),
material_mass_energy(NULL),
conservation(NULL),
modeFile(NULL),
embeddedsurface(NULL),
embeddedsurfaceCp(NULL),
embeddedsurfaceCf(NULL),
cputiming(NULL),
stateVectors(NULL),
stateMaskVectors(NULL),
residualVectors(NULL),
tscale(NULL),
xscale(NULL),
TavF(NULL),
TavM(NULL),
TavL(NULL),
mX(NULL),
tprevf(-1.0),
tprevl(-1.0),
tinit(-1.0),
tener(-1.0),
tenerold(-1.0),
fpForces(NULL),
fpLift(NULL),
fpTavForces(NULL),
fpTavLift(NULL),
fpHydroStaticForces(NULL),
fpHydroDynamicForces(NULL),
fpHydroStaticLift(NULL),
fpHydroDynamicLift(NULL),
fpResiduals(NULL),
fpMatchPressure(NULL),
fpMatchState(NULL),
fpFluxNorm(NULL),
fpMatVolumes(NULL),
fpMaterialMassEnergy(NULL),
fpConservationErr(NULL),
fpGnForces(NULL),
fpStateRom(NULL),
fpError(NULL),
fpEmbeddedSurface(NULL),
fpCpuTiming(NULL),
fpEmbeddedSurfaceCp(NULL),
fpEmbeddedSurfaceCf(NULL),
Qs(NULL),
Qv(NULL),
Qs_match(NULL),
Qs_match_opt(NULL),
Uref(NULL),
Uref_norm(-1.0),
switchOpt(false),
dMatchPressure(NULL),
dForces(NULL),
dLiftDrag(NULL),
dFluxNorm(NULL),
fpdMatchPressure(NULL),
fpdForces(NULL),
fpdLiftDrag(NULL),
fpdFluxNorm(NULL),
heatfluxes(NULL),
fpHeatFluxes(NULL)
{

//  //temporary initialziaze string for the writing of TempStateDeriv
//  if (iod.sa.tempStateDeriv[0]!=0) {
//    const char* tempStateDerivName = iod.sa.tempStateDeriv;
//    this->tempStateDeriv = new char[strlen(iod.input.prefix) + strlen(tempStateDerivName) + 1];
//    sprintf(fullOptPressureName, "%s%s", iod.input.prefix, this->tempStateDeriv);
//  }

  int i;

  modeFile = 0;
  TavF = 0;
  TavM = 0;
  TavL = 0;
  Qs = 0;
  Qv = 0;
  Qs_match     = 0;
  Qs_match_opt = 0;
  Uref = 0;
  Uref_norm = -1.0;

  for (i=0; i<PostFcn::AVSSIZE; ++i) AvQs[i] = 0;
  for (i=0; i<PostFcn::AVVSIZE; ++i) AvQv[i] = 0;

  steady = !iod.problem.type[ProblemData::UNSTEADY];
  com = domain->getCommunicator();

  stateOutputFreqTime = iod.output.rom.stateOutputFreqTime;
  stateOutputFreqNewton = iod.output.rom.stateOutputFreqNewton;
  residualOutputFreqTime = iod.output.rom.residualOutputFreqTime;
  residualOutputMaxNewton = iod.output.rom.residualOutputMaxNewton;
  fdResiduals = (iod.output.rom.fdResiduals == ROMOutputData::FD_RESIDUALS_ON) ? true : false;
  fdResidualsLimit = (iod.output.rom.fdResidualsLimit == ROMOutputData::FD_RESIDUALS_LIMIT_ON) ? true : false;

  externalSI = (iod.embed.surrogateinterface == EmbeddedFramework::EXTERNAL) ? true : false;

  int sp = strlen(iod.output.transient.prefix) + 1;
  int spn = strlen(iod.output.transient.probes.prefix) + 1;
  int sprom = strlen(iod.output.rom.prefix) + 1;

  if (iod.input.optimalPressureFile[0]!=0) {
    const char* optPressureName = iod.input.optimalPressureFile;
    fullOptPressureName = new char[strlen(iod.input.prefix) + strlen(optPressureName) + 1];
    sprintf(fullOptPressureName, "%s%s", iod.input.prefix, optPressureName);

    if ((iod.input.optPressureDim == InputData::NON_DIMENSIONAL) || (iod.input.optPressureDim == InputData::NONE && iod.problem.mode==ProblemData::NON_DIMENSIONAL))
      optPressureDimensional=false;
    else
      optPressureDimensional=true;

    if (!Qs_match_opt) Qs_match_opt = new DistSVec<double,1>(domain->getNodeDistInfo());
    com->fprintf(stdout, "\nReading optimal pressure distribution from %s\n", fullOptPressureName);
    domain->readVectorFromFile(fullOptPressureName,0,0,*Qs_match_opt);
  }

  if (iod.input.matchStateFile[0]!=0) {
    const char* matchStateName = iod.input.matchStateFile;
    char *matchStatePath = new char[strlen(iod.input.prefix) + strlen(matchStateName) + 1];
    sprintf(matchStatePath, "%s%s", iod.input.prefix, matchStateName);

    if (!Uref) Uref = new DistSVec<double,dim>(domain->getNodeDistInfo());
    com->fprintf(stdout, "\nReading comparison state %s\n", matchStatePath);
    domain->readVectorFromFile(matchStatePath,0,0,*Uref);

    delete [] matchStatePath;
  }

  for (i=0; i<PostFcn::SSIZE; ++i) {
    sscale[i] = 1.0;
    scalars[i] = 0;
    nodal_scalars[i] = 0;
  }
  for (i=0; i<PostFcn::AVSSIZE; ++i) {
    avsscale[i] = 1.0;
    avscalars[i] = 0;
  }
  for (i=0; i<PostFcn::VSIZE; ++i) {
    vscale[i] = 1.0;
    vectors[i] = 0;
    nodal_vectors[i] = 0;
  }
  for (i=0; i<PostFcn::AVVSIZE; ++i) {
    avvscale[i] = 1.0;
    avvectors[i] = 0;
  }

  sscale[PostFcn::DENSITY] = iod.ref.rv.density;
  sscale[PostFcn::PRESSURE] = iod.ref.rv.pressure;
  sscale[PostFcn::TEMPERATURE] = iod.ref.rv.temperature;
  vscale[PostFcn::VELOCITY] = iod.ref.rv.velocity;
  vscale[PostFcn::DISPLACEMENT] = iod.ref.rv.tlength;

  if (iod.output.transient.density[0] != 0) {
    scalars[PostFcn::DENSITY] = new char[sp + strlen(iod.output.transient.density)];
    sprintf(scalars[PostFcn::DENSITY], "%s%s",
            iod.output.transient.prefix, iod.output.transient.density);
  }
  if (iod.output.transient.tavdensity[0] != 0) {
    avsscale[PostFcn::DENSITYAVG] = iod.ref.rv.density;
    avscalars[PostFcn::DENSITYAVG] = new char[sp + strlen(iod.output.transient.tavdensity)];
    sprintf(avscalars[PostFcn::DENSITYAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavdensity);
  }
  if (iod.output.transient.mach[0] != 0) {
    scalars[PostFcn::MACH] = new char[sp + strlen(iod.output.transient.mach)];
    sprintf(scalars[PostFcn::MACH], "%s%s",
            iod.output.transient.prefix, iod.output.transient.mach);
  }
  if (iod.output.transient.wtmach[0] != 0) {
    scalars[PostFcn::WTMACH] = new char[sp + strlen(iod.output.transient.wtmach)];
    sprintf(scalars[PostFcn::WTMACH], "%s%s", iod.output.transient.prefix, iod.output.transient.wtmach);
  }
  if (iod.output.transient.wtspeed[0] != 0) {
    scalars[PostFcn::WTSPEED] = new char[sp + strlen(iod.output.transient.wtmach)];
    sprintf(scalars[PostFcn::WTSPEED], "%s%s", iod.output.transient.prefix, iod.output.transient.wtspeed);
    sscale[PostFcn::WTSPEED] = iod.ref.rv.velocity;
  }
  if (iod.output.transient.speed[0] != 0) {
    sscale[PostFcn::SPEED] = iod.ref.rv.velocity;
    scalars[PostFcn::SPEED] = new char[sp + strlen(iod.output.transient.speed)];
    sprintf(scalars[PostFcn::SPEED], "%s%s",
            iod.output.transient.prefix, iod.output.transient.speed);
  }
  if (iod.output.transient.tavmach[0] != 0) {
    avscalars[PostFcn::MACHAVG] = new char[sp + strlen(iod.output.transient.tavmach)];
    sprintf(avscalars[PostFcn::MACHAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavmach);
  }
  if (iod.output.transient.pressure[0] != 0) {
    scalars[PostFcn::PRESSURE] = new char[sp + strlen(iod.output.transient.pressure)];
    sprintf(scalars[PostFcn::PRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.pressure);
  }
  if (iod.output.transient.diffpressure[0] != 0) {
    sscale[PostFcn::DIFFPRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::DIFFPRESSURE] = new char[sp + strlen(iod.output.transient.diffpressure)];
    sprintf(scalars[PostFcn::DIFFPRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.diffpressure);
  }
  if (iod.output.transient.tavpressure[0] != 0) {
    avsscale[PostFcn::PRESSUREAVG] = iod.ref.rv.pressure;
    avscalars[PostFcn::PRESSUREAVG] = new char[sp + strlen(iod.output.transient.tavpressure)];
    sprintf(avscalars[PostFcn::PRESSUREAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavpressure);
  }
  if (iod.output.transient.hydrostaticpressure[0] != 0) {
    sscale[PostFcn::HYDROSTATICPRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::HYDROSTATICPRESSURE] = new char[sp + strlen(iod.output.transient.hydrostaticpressure)];
    sprintf(scalars[PostFcn::HYDROSTATICPRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.hydrostaticpressure);
  }
  if (iod.output.transient.hydrodynamicpressure[0] != 0) {
    sscale[PostFcn::HYDRODYNAMICPRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::HYDRODYNAMICPRESSURE] = new char[sp + strlen(iod.output.transient.hydrodynamicpressure)];
    sprintf(scalars[PostFcn::HYDRODYNAMICPRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.hydrodynamicpressure);
  }
  if (iod.output.transient.pressurecoefficient[0] != 0) {
    sscale[PostFcn::PRESSURECOEFFICIENT] = 1.0;
    scalars[PostFcn::PRESSURECOEFFICIENT] = new char[sp + strlen(iod.output.transient.pressurecoefficient)];
    sprintf(scalars[PostFcn::PRESSURECOEFFICIENT], "%s%s",
            iod.output.transient.prefix, iod.output.transient.pressurecoefficient);
  }
  if (iod.output.transient.temperature[0] != 0) {
//    sscale[PostFcn::TEMPERATURE] = 1;
    scalars[PostFcn::TEMPERATURE] = new char[sp + strlen(iod.output.transient.temperature)];
    sprintf(scalars[PostFcn::TEMPERATURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.temperature);
  }
  if (iod.output.transient.tavtemperature[0] != 0) {
    avsscale[PostFcn::TEMPERATUREAVG] = iod.ref.rv.temperature;
    avscalars[PostFcn::TEMPERATUREAVG] = new char[sp + strlen(iod.output.transient.tavtemperature)];
    sprintf(avscalars[PostFcn::TEMPERATUREAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavtemperature);
  }
  if (iod.output.transient.totalpressure[0] != 0) {
    sscale[PostFcn::TOTPRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::TOTPRESSURE] = new char[sp + strlen(iod.output.transient.totalpressure)];
    sprintf(scalars[PostFcn::TOTPRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.totalpressure);
  }
  if (iod.output.transient.tavtotalpressure[0] != 0) {
    avsscale[PostFcn::TOTPRESSUREAVG] = iod.ref.rv.pressure;
    avscalars[PostFcn::TOTPRESSUREAVG] = new char[sp + strlen(iod.output.transient.tavtotalpressure)];
    sprintf(avscalars[PostFcn::TOTPRESSUREAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavtotalpressure);
  }
  if (iod.output.transient.vorticity[0] != 0) {
    sscale[PostFcn::VORTICITY] = iod.ref.rv.velocity/iod.ref.rv.tlength;
    scalars[PostFcn::VORTICITY] = new char[sp + strlen(iod.output.transient.vorticity)];
    sprintf(scalars[PostFcn::VORTICITY], "%s%s",
            iod.output.transient.prefix, iod.output.transient.vorticity);
  }
  if (iod.output.transient.tavvorticity[0] != 0) {
    avsscale[PostFcn::VORTICITYAVG] = iod.ref.rv.velocity/iod.ref.rv.tlength;
    avscalars[PostFcn::VORTICITYAVG] = new char[sp + strlen(iod.output.transient.tavvorticity)];
    sprintf(avscalars[PostFcn::VORTICITYAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavvorticity);
  }
  if (iod.output.transient.surfaceheatflux[0] != 0) {
    sscale[PostFcn::SURFACE_HEAT_FLUX] = iod.ref.rv.power /(iod.ref.rv.length * iod.ref.rv.length);
    scalars[PostFcn::SURFACE_HEAT_FLUX] = new char[sp + strlen(iod.output.transient.surfaceheatflux)];
    sprintf(scalars[PostFcn::SURFACE_HEAT_FLUX], "%s%s",
            iod.output.transient.prefix, iod.output.transient.surfaceheatflux);
  }
  if (iod.output.transient.tempnormalderivative[0] != 0) {
    sscale[PostFcn::TEMPERATURE_NORMAL_DERIVATIVE] = iod.ref.rv.temperature/iod.ref.rv.length;
    scalars[PostFcn::TEMPERATURE_NORMAL_DERIVATIVE] = new char[sp + strlen(iod.output.transient.tempnormalderivative)];
    sprintf(scalars[PostFcn::TEMPERATURE_NORMAL_DERIVATIVE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tempnormalderivative);
  }
  if (iod.output.transient.nutturb[0] != 0) {
    sscale[PostFcn::NUT_TURB] = iod.ref.rv.viscosity_mu/iod.ref.rv.density;
    scalars[PostFcn::NUT_TURB] = new char[sp + strlen(iod.output.transient.nutturb)];
    sprintf(scalars[PostFcn::NUT_TURB], "%s%s",
            iod.output.transient.prefix, iod.output.transient.nutturb);
  }
  if (iod.output.transient.kturb[0] != 0) {
    sscale[PostFcn::K_TURB] = iod.ref.rv.kenergy;
    scalars[PostFcn::K_TURB] = new char[sp + strlen(iod.output.transient.kturb)];
    sprintf(scalars[PostFcn::K_TURB], "%s%s",
            iod.output.transient.prefix, iod.output.transient.kturb);
  }
  if (iod.output.transient.epsturb[0] != 0) {
    sscale[PostFcn::EPS_TURB] = iod.ref.rv.epsilon;
    scalars[PostFcn::EPS_TURB] = new char[sp + strlen(iod.output.transient.epsturb)];
    sprintf(scalars[PostFcn::EPS_TURB], "%s%s",
            iod.output.transient.prefix, iod.output.transient.epsturb);
  }
  if (iod.output.transient.eddyvis[0] != 0) {
    sscale[PostFcn::EDDY_VISCOSITY] = iod.ref.rv.viscosity_mu;
    scalars[PostFcn::EDDY_VISCOSITY] = new char[sp + strlen(iod.output.transient.eddyvis)];
    sprintf(scalars[PostFcn::EDDY_VISCOSITY], "%s%s",
            iod.output.transient.prefix, iod.output.transient.eddyvis);
  }
  if (iod.output.transient.dplus[0] != 0) {
#if defined(HEAT_FLUX)
    /*
    double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
    double dT = iod.bc.wall.temperature - 1.0 / (gam*(gam-1.0)*iod.bc.inlet.mach*iod.bc.inlet.mach);
    sscale[PostFcn::DELTA_PLUS] = iod.ref.reynolds * iod.eqs.thermalCondModel.prandtl / (gam * dT);
    */
    sscale[PostFcn::DELTA_PLUS] = iod.ref.rv.tpower / (iod.ref.length*iod.ref.length);
#endif
    scalars[PostFcn::DELTA_PLUS] = new char[sp + strlen(iod.output.transient.dplus)];
    sprintf(scalars[PostFcn::DELTA_PLUS], "%s%s",
            iod.output.transient.prefix, iod.output.transient.dplus);
  }
  if (iod.output.transient.sfric[0] != 0) {
    scalars[PostFcn::SKIN_FRICTION] = new char[sp + strlen(iod.output.transient.sfric)];
    sprintf(scalars[PostFcn::SKIN_FRICTION], "%s%s",
            iod.output.transient.prefix, iod.output.transient.sfric);
  }
  if (iod.output.transient.tavsfric[0] != 0) {
    avscalars[PostFcn::SKIN_FRICTIONAVG] = new char[sp + strlen(iod.output.transient.tavsfric)];
    sprintf(avscalars[PostFcn::SKIN_FRICTIONAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavsfric);
  }
  if (iod.output.transient.psensor[0] != 0) {
    scalars[PostFcn::PSENSOR] = new char[sp + strlen(iod.output.transient.psensor)];
    sprintf(scalars[PostFcn::PSENSOR], "%s%s",
            iod.output.transient.prefix, iod.output.transient.psensor);
  }
  if (iod.output.transient.csdles[0] != 0) {
    scalars[PostFcn::CSDLES] = new char[sp + strlen(iod.output.transient.csdles)];
    sprintf(scalars[PostFcn::CSDLES], "%s%s",
            iod.output.transient.prefix, iod.output.transient.csdles);
  }
  if (iod.output.transient.tavcsdles[0] != 0) {
    avscalars[PostFcn::CSDLESAVG] = new char[sp + strlen(iod.output.transient.tavcsdles)];
    sprintf(avscalars[PostFcn::CSDLESAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavcsdles);
  }
  if (iod.output.transient.csdvms[0] != 0) {
    scalars[PostFcn::CSDVMS] = new char[sp + strlen(iod.output.transient.csdvms)];
    sprintf(scalars[PostFcn::CSDVMS], "%s%s",
            iod.output.transient.prefix, iod.output.transient.csdvms);
  }
  if (iod.output.transient.tavcsdvms[0] != 0) {
    avscalars[PostFcn::CSDVMSAVG] = new char[sp + strlen(iod.output.transient.tavcsdvms)];
    sprintf(avscalars[PostFcn::CSDVMSAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavcsdvms);
  }
  if (iod.output.transient.mutOmu[0] != 0) {
    scalars[PostFcn::MUT_OVER_MU] = new char[sp + strlen(iod.output.transient.mutOmu)];
    sprintf(scalars[PostFcn::MUT_OVER_MU], "%s%s",
            iod.output.transient.prefix, iod.output.transient.mutOmu);
  }
  if (iod.output.transient.philevel[0] != 0) {
    sscale[PostFcn::PHILEVEL] = 1.0;
    scalars[PostFcn::PHILEVEL] = new char[sp + strlen(iod.output.transient.philevel)];
    sprintf(scalars[PostFcn::PHILEVEL], "%s%s",
            iod.output.transient.prefix, iod.output.transient.philevel);
  }
  if (iod.output.transient.philevel2[0] != 0) {
    sscale[PostFcn::PHILEVEL2] = 1.0;
    scalars[PostFcn::PHILEVEL2] = new char[sp + strlen(iod.output.transient.philevel2)];
    sprintf(scalars[PostFcn::PHILEVEL2], "%s%s",
            iod.output.transient.prefix, iod.output.transient.philevel2);
  }
  if (iod.output.transient.fluidid[0] != 0) {
    sscale[PostFcn::FLUIDID] = 1.0;
    scalars[PostFcn::FLUIDID] = new char[sp + strlen(iod.output.transient.fluidid)];
    sprintf(scalars[PostFcn::FLUIDID], "%s%s",
            iod.output.transient.prefix, iod.output.transient.fluidid);
  }
  if (iod.output.transient.d2wall[0] != 0) {
    sscale[PostFcn::D2WALL] = 1.0;
    scalars[PostFcn::D2WALL] = new char[sp + strlen(iod.output.transient.d2wall)];
    sprintf(scalars[PostFcn::D2WALL], "%s%s",
            iod.output.transient.prefix, iod.output.transient.d2wall);
  }
  if (iod.output.transient.controlvolume[0] != 0) {
    sscale[PostFcn::CONTROL_VOLUME] = iod.ref.rv.length * iod.ref.rv.length * iod.ref.rv.length;
    scalars[PostFcn::CONTROL_VOLUME] = new char[sp + strlen(iod.output.transient.controlvolume)];
    sprintf(scalars[PostFcn::CONTROL_VOLUME], "%s%s",
            iod.output.transient.prefix, iod.output.transient.controlvolume);
  }
  if (iod.output.transient.velocity[0] != 0) {
    vectors[PostFcn::VELOCITY] = new char[sp + strlen(iod.output.transient.velocity)];
    sprintf(vectors[PostFcn::VELOCITY], "%s%s",
            iod.output.transient.prefix, iod.output.transient.velocity);
  }
  if (iod.output.transient.tavvelocity[0] != 0) {
    avvscale[PostFcn::VELOCITYAVG] = iod.ref.rv.velocity;
    avvectors[PostFcn::VELOCITYAVG] = new char[sp + strlen(iod.output.transient.tavvelocity)];
    sprintf(avvectors[PostFcn::VELOCITYAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavvelocity);
  }
  if (iod.output.transient.displacement[0] != 0) {
    vectors[PostFcn::DISPLACEMENT] = new char[sp + strlen(iod.output.transient.displacement)];
    sprintf(vectors[PostFcn::DISPLACEMENT], "%s%s",
            iod.output.transient.prefix, iod.output.transient.displacement);
  }
  if (iod.output.transient.tavdisplacement[0] != 0) {
    avvscale[PostFcn::DISPLACEMENTAVG] = iod.ref.rv.tlength;
    avvectors[PostFcn::DISPLACEMENTAVG] = new char[sp + strlen(iod.output.transient.tavdisplacement)];
    sprintf(avvectors[PostFcn::DISPLACEMENTAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavdisplacement);
  }

  if (iod.output.transient.flightDisplacement[0] != 0) {
    vscale[PostFcn::FLIGHTDISPLACEMENT] = iod.ref.rv.tlength;
    vectors[PostFcn::FLIGHTDISPLACEMENT] = new char[sp + strlen(iod.output.transient.flightDisplacement)];
    sprintf(vectors[PostFcn::FLIGHTDISPLACEMENT], "%s%s",
            iod.output.transient.prefix, iod.output.transient.flightDisplacement);
  }

  if (iod.output.transient.localFlightDisplacement[0] != 0) {
    vscale[PostFcn::LOCALFLIGHTDISPLACEMENT] = iod.ref.rv.tlength;
    vectors[PostFcn::LOCALFLIGHTDISPLACEMENT] = new char[sp + strlen(iod.output.transient.localFlightDisplacement)];
    sprintf(vectors[PostFcn::LOCALFLIGHTDISPLACEMENT], "%s%s",
            iod.output.transient.prefix, iod.output.transient.localFlightDisplacement);
  }

  if (iod.output.transient.forces[0] != 0) {
    forces = new char[sp + strlen(iod.output.transient.forces)];
    sprintf(forces, "%s%s", iod.output.transient.prefix, iod.output.transient.forces);
  }
  else
    forces = 0;

  if (iod.output.transient.tavforces[0] != 0) {
    tavforces = new char[sp + strlen(iod.output.transient.tavforces)];
    sprintf(tavforces, "%s%s", iod.output.transient.prefix, iod.output.transient.tavforces);
  }
  else
    tavforces = 0;

  if (iod.output.transient.hydrostaticforces[0] != 0) {
    hydrostaticforces = new char[sp + strlen(iod.output.transient.hydrostaticforces)];
    sprintf(hydrostaticforces, "%s%s", iod.output.transient.prefix, iod.output.transient.hydrostaticforces);
  }
  else
    hydrostaticforces = 0;

  if (iod.output.transient.hydrodynamicforces[0] != 0) {
    hydrodynamicforces = new char[sp + strlen(iod.output.transient.hydrodynamicforces)];
    sprintf(hydrodynamicforces, "%s%s", iod.output.transient.prefix, iod.output.transient.hydrodynamicforces);
  }
  else
    hydrodynamicforces = 0;

  mX = 0;
  if (iod.output.transient.generalizedforces[0] != 0 && strcmp(iod.input.strModesFile, "") != 0) {
    generalizedforces = new char[sp + strlen(iod.output.transient.generalizedforces)];
    sprintf(generalizedforces, "%s%s", iod.output.transient.prefix, iod.output.transient.generalizedforces);
    modeFile = new char [MAXLINE];
    sprintf(modeFile, "%s%s", iod.input.prefix, iod.input.strModesFile);
    DistSVec<double, dim> tmpVec(domain->getNodeDistInfo());
    DistSVec<double, 3> Xtmp(domain->getNodeDistInfo());
    double f;

    if (strcmp(modeFile, "") != 0)  {
      com->fprintf(stderr, " ... Reading Modefile %s\n", modeFile);
      domain->readVectorFromFile(modeFile, 0, &f, Xtmp);
      int nStrMode = int(f);

      // We read the modal deformations
      mX = new VecSet<DistSVec<double,3> >(nStrMode, domain->getNodeDistInfo());

      for(int iMode = 0; iMode < nStrMode; ++iMode)
        domain->readVectorFromFile(modeFile, iMode+1, &f, (*mX)[iMode]);
     }
     else
      mX = 0;
  }
  else
    {
      generalizedforces = 0;
      modeFile =0;
    }

  if (iod.output.transient.generalizedforces[0] != 0 && !(iod.problem.type[ProblemData::FORCED]) && strcmp(iod.input.strModesFile, "") == 0)
    fprintf(stderr, "Error : StrModes file for Generalized Forces not given.. Aborting !! \n");


  if (iod.output.transient.lift[0] != 0){
    lift = new char[sp + strlen(iod.output.transient.lift)];
    sprintf(lift, "%s%s", iod.output.transient.prefix, iod.output.transient.lift);
  }
  else
    lift = 0;

  if (iod.output.transient.matchpressure[0] != 0){
    matchpressure = new char[sp + strlen(iod.output.transient.matchpressure)];
    sprintf(matchpressure, "%s%s", iod.output.transient.prefix, iod.output.transient.matchpressure);
  }
  else
    matchpressure = 0;

  if (iod.output.transient.matchstate[0] != 0){
    matchstate = new char[sp + strlen(iod.output.transient.matchstate)];
    sprintf(matchstate, "%s%s", iod.output.transient.prefix, iod.output.transient.matchstate);
  }
  else
    matchstate = 0;

  if (iod.output.transient.fluxnorm[0] != 0){
    fluxnorm = new char[sp + strlen(iod.output.transient.fluxnorm)];
    sprintf(fluxnorm, "%s%s", iod.output.transient.prefix, iod.output.transient.fluxnorm);
  }
  else
    fluxnorm = 0;

  if (iod.output.transient.tavlift[0] != 0){
    tavlift = new char[sp + strlen(iod.output.transient.tavlift)];
    sprintf(tavlift, "%s%s", iod.output.transient.prefix, iod.output.transient.tavlift);
  }
  else
    tavlift = 0;

  if (iod.output.transient.hydrostaticlift[0] != 0) {
    hydrostaticlift = new char[sp + strlen(iod.output.transient.hydrostaticlift)];
    sprintf(hydrostaticlift, "%s%s", iod.output.transient.prefix, iod.output.transient.hydrostaticlift);
  }
  else
    hydrostaticlift = 0;

  if (iod.output.transient.hydrodynamiclift[0] != 0) {
    hydrodynamiclift = new char[sp + strlen(iod.output.transient.hydrodynamiclift)];
    sprintf(hydrodynamiclift, "%s%s", iod.output.transient.prefix, iod.output.transient.hydrodynamiclift);
  }
  else
    hydrodynamiclift = 0;

  if (iod.output.transient.heatfluxes[0] != 0) {
    heatfluxes = new char[sp + strlen(iod.output.transient.heatfluxes)];
    sprintf(heatfluxes, "%s%s", iod.output.transient.prefix, iod.output.transient.heatfluxes);
  }
  else
    heatfluxes = 0;

  if (iod.output.transient.residuals[0] != 0) {
    residuals = new char[sp + strlen(iod.output.transient.residuals)];
    sprintf(residuals, "%s%s", iod.output.transient.prefix, iod.output.transient.residuals);
  }
  else
    residuals = 0;

  if (iod.output.transient.materialVolumes[0] != 0) {
    material_volumes = new char[sp + strlen(iod.output.transient.materialVolumes)];
    sprintf(material_volumes, "%s%s", iod.output.transient.prefix, iod.output.transient.materialVolumes);
  }
  else
    material_volumes = 0;

  if (iod.output.transient.materialMassEnergy[0] != 0) {
    material_mass_energy = new char[sp + strlen(iod.output.transient.materialMassEnergy)];
    sprintf(material_mass_energy, "%s%s", iod.output.transient.prefix, iod.output.transient.materialMassEnergy);
  }
  else
    material_mass_energy = 0;

  if (iod.output.transient.embeddedsurface[0] != 0) {
    embeddedsurface = new char[sp + strlen(iod.output.transient.embeddedsurface)];
    sprintf(embeddedsurface, "%s%s", iod.output.transient.prefix, iod.output.transient.embeddedsurface);
  }
  else
    embeddedsurface = 0;

  if (iod.output.transient.cputiming[0] != 0) {
    cputiming = new char[sp + strlen(iod.output.transient.cputiming)];
    sprintf(cputiming, "%s%s", iod.output.transient.prefix, iod.output.transient.cputiming);
  }
  else
    cputiming = 0;

  if (iod.output.transient.conservation[0] != 0) {
    conservation = new char[sp + strlen(iod.output.transient.conservation)];
    sprintf(conservation, "%s%s", iod.output.transient.prefix, iod.output.transient.conservation);
  }
  else
    conservation = 0;

  // for embedded ROMs (only active in embeddedTsDesc), Lei Lei, 02/01/2016
  if (iod.output.rom.stateMaskVector[0] != 0) {
    stateMaskVectors = new char[sprom + strlen(iod.output.rom.stateMaskVector)];
    sprintf(stateMaskVectors, "%s%s", iod.output.rom.prefix, iod.output.rom.stateMaskVector);
  }
  else
    stateMaskVectors = 0;

  // for ROMs (only active for implicit time stepping)
  if (iod.output.rom.stateVector[0] != 0) {
    stateVectors = new char[sprom + strlen(iod.output.rom.stateVector)];
    sprintf(stateVectors, "%s%s", iod.output.rom.prefix, iod.output.rom.stateVector);
  }
  else
    stateVectors = 0;

  // for ROMs (only active for implicit time stepping)
  if (iod.output.rom.residualVector[0] != 0) {
    residualVectors = new char[sprom + strlen(iod.output.rom.residualVector)];
    sprintf(residualVectors, "%s%s", iod.output.rom.prefix, iod.output.rom.residualVector);
  }
  else
    residualVectors = 0;

  // If we want to output the pressure coeff. in the embedded framework d2d
  if ( (iod.problem.framework==ProblemData::EMBEDDED ||
	iod.problem.framework==ProblemData::EMBEDDEDALE) &&
       iod.output.transient.pressurecoefficient[0] != 0) {
    embeddedsurfaceCp = new char[sp + strlen(iod.output.transient.pressurecoefficient)];
    sprintf(embeddedsurfaceCp, "%s%s%s", iod.output.transient.prefix, "emb_", iod.output.transient.pressurecoefficient);
  }  else
    embeddedsurfaceCp = 0;

  // If we want to output the skin friction coeff. in the embedded framework d2d
  if ( (iod.problem.framework==ProblemData::EMBEDDED ||
	iod.problem.framework==ProblemData::EMBEDDEDALE) &&
       iod.output.transient.sfric[0] != 0) {
    embeddedsurfaceCf = new char[sp + strlen(iod.output.transient.sfric)];
    sprintf(embeddedsurfaceCf, "%s%s%s", iod.output.transient.prefix, "emb_", iod.output.transient.sfric);
  }  else
    embeddedsurfaceCf = 0;

  it0 = iod.restart.iteration;
  //std::cout << "it0 = " << it0 << std::endl;
  numFluidPhases = iod.eqs.numPhase;
  frequency = iod.output.transient.frequency;
  frequency_dt = iod.output.transient.frequency_dt;
  prtout = iod.restart.iteration ? frequency_dt*(floor(iod.restart.etime/frequency_dt)+1.0) : iod.restart.etime;
  length = iod.output.transient.length;
  surface = iod.output.transient.surface;
  x0[0] = iod.output.transient.x0;
  x0[1] = iod.output.transient.y0;
  x0[2] = iod.output.transient.z0;

  fpCpuTiming = 0;
  fpResiduals = 0;
  fpMatchState = 0;
  fpMatchPressure = 0;
  fpFluxNorm      = 0;
  fpMatVolumes = 0;
  fpConservationErr = 0;
  fpGnForces  = 0;
  fpEmbeddedSurface = 0;
  fpCpuTiming = 0;

  fpEmbeddedSurface = 0;
  fpCpuTiming = 0;

  //
  fpEmbeddedSurfaceCp = 0;
  fpEmbeddedSurfaceCf = 0;

  int nSurf = postOp->getNumSurf();
  int nSurfHF = postOp->getNumSurfHF();
  fpForces             = new FILE *[nSurf];
  fpLift               = new FILE *[nSurf];
  fpTavForces          = new FILE *[nSurf];
  fpTavLift            = new FILE *[nSurf];
  fpHydroStaticForces  = new FILE *[nSurf];
  fpHydroDynamicForces = new FILE *[nSurf];
  fpHydroStaticLift    = new FILE *[nSurf];
  fpHydroDynamicLift   = new FILE *[nSurf];
  fpHeatFluxes         = new FILE *[nSurfHF];

  for (int iSurf = 0; iSurf < nSurf; iSurf++)  {
    fpForces[iSurf]          = 0;
    fpLift[iSurf]            = 0;
    fpTavForces[iSurf]       = 0;
    fpTavLift[iSurf]         = 0;
    fpHydroStaticForces[iSurf]  = 0;
    fpHydroDynamicForces[iSurf] = 0;
    fpHydroStaticLift[iSurf]    = 0;
    fpHydroDynamicLift[iSurf]   = 0;
  }

  for (int iSurf = 0; iSurf < nSurfHF; iSurf++)  {
    fpHeatFluxes[iSurf]         = 0;
  }


// Included (MB)
  if (iod.output.transient.velocitynorm[0] != 0) {
    sscale[PostFcn::VELOCITY_NORM] = iod.ref.rv.velocity;
    scalars[PostFcn::VELOCITY_NORM] = new char[sp + strlen(iod.output.transient.velocitynorm)];
    sprintf(scalars[PostFcn::VELOCITY_NORM], "%s%s",
            iod.output.transient.prefix, iod.output.transient.velocitynorm);
  }

  if (iod.problem.alltype == ProblemData::_SHAPE_OPTIMIZATION_ ||
      iod.problem.alltype == ProblemData::_AEROELASTIC_SHAPE_OPTIMIZATION_ ||
      iod.problem.alltype == ProblemData::_ROM_SHAPE_OPTIMIZATION_ ||
      iod.problem.alltype == ProblemData::_SENSITIVITY_ANALYSIS_ ) { //TODO CHECK if needed

  int dsp = strlen(iod.output.transient.prefix) + 1;

  if (iod.output.transient.dSolutions[0] != 0) {
    dSolutions = new char[dsp + strlen(iod.output.transient.dSolutions)];
    sprintf(dSolutions, "%s%s", iod.output.transient.prefix, iod.output.transient.dSolutions);
  }
  else
    dSolutions = 0;

  int i;
  for (i=0; i<PostFcn::DSSIZE; ++i) {
    dSscale[i] = 1.0;
    dScalars[i] = 0;
  }
  for (i=0; i<PostFcn::DVSIZE; ++i) {
    dVscale[i] = 1.0;
    dVectors[i] = 0;
  }

  if (iod.output.transient.dDensity[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_DENSITY] = iod.ref.rv.density;
    dScalars[PostFcn::DERIVATIVE_DENSITY] = new char[dsp + strlen(iod.output.transient.dDensity)];
    sprintf(dScalars[PostFcn::DERIVATIVE_DENSITY], "%s%s",
            iod.output.transient.prefix, iod.output.transient.dDensity);
  }
  if (iod.output.transient.dMach[0] != 0) {
    dScalars[PostFcn::DERIVATIVE_MACH] = new char[dsp + strlen(iod.output.transient.dMach)];
    sprintf(dScalars[PostFcn::DERIVATIVE_MACH], "%s%s",
            iod.output.transient.prefix, iod.output.transient.dMach);
  }
  if (iod.output.transient.dPressure[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_PRESSURE] = iod.ref.rv.pressure;
    dScalars[PostFcn::DERIVATIVE_PRESSURE] = new char[dsp + strlen(iod.output.transient.dPressure)];
    sprintf(dScalars[PostFcn::DERIVATIVE_PRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.dPressure);
  }
  if (iod.output.transient.dTemperature[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_TEMPERATURE] = iod.ref.rv.temperature;
    dScalars[PostFcn::DERIVATIVE_TEMPERATURE] = new char[dsp + strlen(iod.output.transient.dTemperature)];
    sprintf(dScalars[PostFcn::DERIVATIVE_TEMPERATURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.dTemperature);
  }
  if (iod.output.transient.dTotalpressure[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_TOTPRESSURE] = iod.ref.rv.pressure;
    dScalars[PostFcn::DERIVATIVE_TOTPRESSURE] = new char[dsp + strlen(iod.output.transient.dTotalpressure)];
    sprintf(dScalars[PostFcn::DERIVATIVE_TOTPRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.dTotalpressure);
  }
  if (iod.output.transient.dNutturb[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_NUT_TURB] = iod.ref.rv.viscosity_mu/iod.ref.rv.density;
    dScalars[PostFcn::DERIVATIVE_NUT_TURB] = new char[dsp + strlen(iod.output.transient.dNutturb)];
    sprintf(dScalars[PostFcn::DERIVATIVE_NUT_TURB], "%s%s",
            iod.output.transient.prefix, iod.output.transient.dNutturb);
  }
  if (iod.output.transient.dEddyvis[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_EDDY_VISCOSITY] = iod.ref.rv.viscosity_mu;
    dScalars[PostFcn::DERIVATIVE_EDDY_VISCOSITY] = new char[dsp + strlen(iod.output.transient.dEddyvis)];
    sprintf(dScalars[PostFcn::DERIVATIVE_EDDY_VISCOSITY], "%s%s",
            iod.output.transient.prefix, iod.output.transient.dEddyvis);
  }
  if (iod.output.transient.dVelocityScalar[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_VELOCITY_SCALAR] = iod.ref.rv.velocity;
    dScalars[PostFcn::DERIVATIVE_VELOCITY_SCALAR] = new char[dsp + strlen(iod.output.transient.dVelocityScalar)];
    sprintf(dScalars[PostFcn::DERIVATIVE_VELOCITY_SCALAR], "%s%s",
            iod.output.transient.prefix, iod.output.transient.dVelocityScalar);
  }
  if (iod.output.transient.dVelocityVector[0] != 0) {
    dVscale[PostFcn::DERIVATIVE_VELOCITY_VECTOR] = iod.ref.rv.velocity;
    dVectors[PostFcn::DERIVATIVE_VELOCITY_VECTOR] = new char[dsp + strlen(iod.output.transient.dVelocityVector)];
    sprintf(dVectors[PostFcn::DERIVATIVE_VELOCITY_VECTOR], "%s%s",
            iod.output.transient.prefix, iod.output.transient.dVelocityVector);
  }
  if (iod.output.transient.dDisplacement[0] != 0) {
    dVscale[PostFcn::DERIVATIVE_DISPLACEMENT] = iod.ref.rv.tlength;
    dVectors[PostFcn::DERIVATIVE_DISPLACEMENT] = new char[dsp + strlen(iod.output.transient.dDisplacement)];
    sprintf(dVectors[PostFcn::DERIVATIVE_DISPLACEMENT], "%s%s",
            iod.output.transient.prefix, iod.output.transient.dDisplacement);
  }
  if (iod.output.transient.dMatchPressure[0] != 0 && iod.input.optimalPressureFile) {
    dMatchPressure = new char[dsp + strlen(iod.output.transient.dMatchPressure)];
    sprintf(dMatchPressure, "%s%s", iod.output.transient.prefix, iod.output.transient.dMatchPressure);
  }
  else
    dMatchPressure = 0;

  fpdMatchPressure = 0;


  if (iod.output.transient.dLiftDrag[0] != 0) {
    dLiftDrag = new char[dsp + strlen(iod.output.transient.dLiftDrag)];
    sprintf(dLiftDrag, "%s%s", iod.output.transient.prefix, iod.output.transient.dLiftDrag);
  }
  else
    dLiftDrag = 0;

  if (iod.output.transient.dLiftx[0] != 0) {
    dLiftx = new char[dsp + strlen(iod.output.transient.dLiftx)];
    sprintf(dLiftx, "%s%s", iod.output.transient.prefix, iod.output.transient.dLiftx);
  }
  else dLiftx = 0;

  if (iod.output.transient.dLifty[0] != 0) {
    dLifty = new char[dsp + strlen(iod.output.transient.dLifty)];
    sprintf(dLifty, "%s%s", iod.output.transient.prefix, iod.output.transient.dLifty);
  }
  else dLifty = 0;

  if (iod.output.transient.dLiftz[0] != 0) {
    dLiftz = new char[dsp + strlen(iod.output.transient.dLiftz)];
    sprintf(dLiftz, "%s%s", iod.output.transient.prefix, iod.output.transient.dLiftz);
  }
  else dLiftz = 0;


  fpdLiftDrag = 0;
  fpdLiftx = 0;
  fpdLifty = 0;
  fpdLiftz = 0;

  if (iod.output.rom.dFluxNorm[0] != 0) {
    dFluxNorm = new char[sprom+ strlen(iod.output.rom.dFluxNorm)];
    sprintf(dFluxNorm, "%s%s", iod.output.rom.prefix, iod.output.rom.dFluxNorm);
  }
  else
    dFluxNorm = 0;

  fpdFluxNorm = 0;

  if (iod.output.transient.dForces[0] != 0) {
    dForces = new char[dsp + strlen(iod.output.transient.dForces)];
    sprintf(dForces, "%s%s", iod.output.transient.prefix, iod.output.transient.dForces);
  }
  else
    dForces = 0;

  fpdForces = 0;

  switchOpt = true;
  }
  else {
    switchOpt = false;
    dSolutions = 0;
    dMatchPressure = 0;
    fpdMatchPressure = 0;
    dLiftDrag = 0;
    fpdLiftDrag = 0;
    dFluxNorm = 0;
    fpdFluxNorm = 0;
    dForces = 0;
    fpdForces = 0;

    int i;
    for (i=0; i<PostFcn::DSSIZE; ++i) {
      dScalars[i] = 0;
    }
    for (i=0; i<PostFcn::DVSIZE; ++i) {
      dVectors[i] = 0;
    }
  }

  // Initialize nodal output structures
  Probes& myProbes = iod.output.transient.probes;
  nodal_output.step = 0;
  nodal_output.results = new double[Probes::MAXNODES*3];
  nodal_output.subId = new int[Probes::MAXNODES];
  nodal_output.locNodeId = new int[Probes::MAXNODES];
  nodal_output.last = new int[Probes::MAXNODES];
  nodal_output.locations.resize(Probes::MAXNODES);

  nodal_output.step = 0;


  for (i = 0; i < Probes::MAXNODES; ++i) {
    nodal_output.locations[i] = Vec3D(myProbes.myNodes[i].locationX/iod.ref.rv.length,
                                      myProbes.myNodes[i].locationY/iod.ref.rv.length,
                                      myProbes.myNodes[i].locationZ/iod.ref.rv.length);

    nodal_output.last[i] = 0;
    if (myProbes.myNodes[i].id >= 0) {
      com->fprintf(stdout,"[Probe] Node %d: NodeId = %d.\n", i+1, myProbes.myNodes[i].id);

      int flag = -1;
      int locid = -1;
      int lis = -1;
//      bool abortOmp = false;
#pragma omp parallel for
      for (int iSub = 0; iSub < dom->getNumLocSub(); ++iSub) {
        //#pragma omp flush (abortOmp)
        //if (!abortOmp) {
          locid = dom->getSubDomain()[iSub]->getLocalNodeNum( myProbes.myNodes[i].id-1 );
//          fprintf(stdout,"locid = %i\n",locid);
          if (locid >= 0) {
            lis = iSub;
            flag = com->cpuNum();
           // abortOmp = true;
           // #pragma omp flush (abortOmp)
            break;
          }
        //}
      }

      com->globalMax(1,&flag);

      if ( flag == com->cpuNum() ) {
        nodal_output.locNodeId[i] = locid;
        nodal_output.subId[i] = lis;
      } else {
        nodal_output.locNodeId[i] = -1;
        nodal_output.subId[i] = -1;
      }
    } else if (myProbes.myNodes[i].locationX < -1.0e10)
      break;
    else
      com->fprintf(stdout,"[Probe] Node %d: Coords = (%e, %e, %e).\n", i+1,
                   myProbes.myNodes[i].locationX, myProbes.myNodes[i].locationY, myProbes.myNodes[i].locationZ);
  }

  com->fprintf(stdout,"[Probe] Number of probing nodes is %d\n",i);

  nodal_output.numNodes = i;


  if (iod.output.transient.probes.density[0] != 0) {
    nodal_scalars[PostFcn::DENSITY] = new char[spn + strlen(iod.output.transient.probes.density)];
    sprintf(nodal_scalars[PostFcn::DENSITY], "%s%s",
            iod.output.transient.probes.prefix, iod.output.transient.probes.density);
  }
  if (iod.output.transient.probes.pressure[0] != 0) {
    nodal_scalars[PostFcn::PRESSURE] = new char[spn + strlen(iod.output.transient.probes.pressure)];
    sprintf(nodal_scalars[PostFcn::PRESSURE], "%s%s",
            iod.output.transient.probes.prefix, iod.output.transient.probes.pressure);
  }
  if (iod.output.transient.probes.diffpressure[0] != 0) {
    nodal_scalars[PostFcn::DIFFPRESSURE] = new char[spn + strlen(iod.output.transient.probes.diffpressure)];
    sprintf(nodal_scalars[PostFcn::DIFFPRESSURE], "%s%s",
            iod.output.transient.probes.prefix, iod.output.transient.probes.diffpressure);
  }
  if (iod.output.transient.probes.temperature[0] != 0) {
    nodal_scalars[PostFcn::TEMPERATURE] = new char[spn + strlen(iod.output.transient.probes.temperature)];
    sprintf(nodal_scalars[PostFcn::TEMPERATURE], "%s%s",
            iod.output.transient.probes.prefix, iod.output.transient.probes.temperature);
  }
  if (iod.output.transient.probes.velocity[0] != 0) {
    nodal_vectors[PostFcn::VELOCITY] = new char[spn + strlen(iod.output.transient.probes.velocity)];
    sprintf(nodal_vectors[PostFcn::VELOCITY], "%s%s",
            iod.output.transient.probes.prefix, iod.output.transient.probes.velocity);
  }
  if (iod.output.transient.probes.displacement[0] != 0) {
    nodal_vectors[PostFcn::DISPLACEMENT] = new char[spn + strlen(iod.output.transient.probes.displacement)];
    sprintf(nodal_vectors[PostFcn::DISPLACEMENT], "%s%s",
            iod.output.transient.probes.prefix, iod.output.transient.probes.displacement);
  }

  map<int, LinePlot *>::iterator it;
  for (it = iod.output.transient.linePlots.dataMap.begin();
       it !=  iod.output.transient.linePlots.dataMap.end();
       ++it) {

    LinePlot& L = *it->second;
    line_output* Lp = new line_output;

    memset(Lp, 0, sizeof(line_output));

    Lp->numPoints = L.numPoints;
    Lp->x0 = L.x0;
    Lp->y0 = L.y0;
    Lp->z0 = L.z0;
    Lp->x1 = L.x1;
    Lp->y1 = L.y1;
    Lp->z1 = L.z1;

    if (L.density[0] != 0) {
      Lp->scalars[PostFcn::DENSITY] = new char[sp + strlen(L.density)];
      sprintf(Lp->scalars[PostFcn::DENSITY], "%s%s",
	      iod.output.transient.prefix, L.density);
    }
    if (L.pressure[0] != 0) {
      Lp->scalars[PostFcn::PRESSURE] = new char[sp + strlen(L.pressure)];
      sprintf(Lp->scalars[PostFcn::PRESSURE], "%s%s",
	      iod.output.transient.prefix, L.pressure);
    }
    if (L.temperature[0] != 0) {
      Lp->scalars[PostFcn::TEMPERATURE] = new char[sp + strlen(L.temperature)];
      sprintf(Lp->scalars[PostFcn::TEMPERATURE], "%s%s",
	      iod.output.transient.prefix, L.temperature);
    }
    if (L.velocity[0] != 0) {
      Lp->vectors[PostFcn::VELOCITY] = new char[sp + strlen(L.velocity)];
      sprintf(Lp->vectors[PostFcn::VELOCITY], "%s%s",
	      iod.output.transient.prefix, L.velocity);
    }
    line_outputs.push_back(Lp);
  }


  tscale = iod.ref.rv.time;
  xscale = iod.ref.rv.length;
}

//------------------------------------------------------------------------------

template<int dim>
TsOutput<dim>::~TsOutput()
{

  for (int i=0; i<PostFcn::AVSSIZE; ++i) {
    if(AvQs[i]) delete AvQs[i];
  }
  for (int i=0; i<PostFcn::AVVSIZE; ++i) {
    if(AvQv[i]) delete AvQv[i];
  }
  if (Qs_match) delete Qs_match;
  if (Qs_match_opt) delete Qs_match_opt;
  if (Qs) delete Qs;
  if (Qv) delete Qv;
  if (Uref) delete Uref;
  if(TavF) delete [] TavF;
  if(TavM) delete [] TavM;
  if(TavL) delete [] TavL;
  if(modeFile) delete [] modeFile;
  if(mX) delete mX;

  if (switchOpt)
    {
      delete[] dMatchPressure;
      delete[] dForces;
      delete[] dLiftDrag;
      delete[] dLiftx;
      delete[] dLifty;
      delete[] dLiftz;
      delete[] dFluxNorm;

      int i;
      for (i=0; i<PostFcn::DSSIZE; ++i) {
        delete[]  dScalars[i];
      }
      for (i=0; i<PostFcn::DVSIZE; ++i) {
        delete[] dVectors[i];
      }
      delete[] dSolutions;
    }

  delete[]  fpForces            ;
  delete[]  fpLift              ;
  delete[]  fpTavForces         ;
  delete[]  fpTavLift           ;
  delete[]  fpHydroStaticForces ;
  delete[]  fpHydroDynamicForces;
  delete[]  fpHydroStaticLift   ;
  delete[]  fpHydroDynamicLift  ;
  delete[]  fpHeatFluxes        ;

  delete[] heatfluxes;
  delete[] residuals;
  delete[] matchpressure;
  if (matchstate) delete[] matchstate;
  delete[] fluxnorm;
  delete[] material_volumes;
  delete[] embeddedsurface;
  delete[] embeddedsurfaceCp;
  delete[] embeddedsurfaceCf;
  delete[] cputiming;
  delete[] conservation;

  delete[] lift;
  delete[] tavlift;
  delete[] hydrostaticlift;
  delete[] hydrodynamiclift;

  delete mX;

  delete[] generalizedforces;
  delete[] modeFile;

  delete[] forces;
  delete[] tavforces;
  delete[] hydrostaticforces;
  delete[] hydrodynamicforces;

  int i;
  for (i=0; i<PostFcn::SSIZE; ++i) {
    delete[] scalars[i];
  }
  for (i=0; i<PostFcn::AVSSIZE; ++i) {
     delete[]  avscalars[i];
  }
  for (i=0; i<PostFcn::VSIZE; ++i) {
     delete[]  vectors[i];
  }
  for (i=0; i<PostFcn::AVVSIZE; ++i) {
     delete[]  avvectors[i];
  }

}

//------------------------------------------------------------------------------

template<int dim>
bool TsOutput<dim>::toWrite(int it, bool lastIt, double t)
{
  if(frequency_dt<=0.0)
    return (((frequency > 0) && (it % frequency == 0)) || lastIt);

  return (it==0 || t>=prtout || lastIt);

}

//------------------------------------------------------------------------------

template<int dim>
int TsOutput<dim>::getStep(int it, bool lastIt, double t)
{
  int step = 0;
  if(frequency_dt<=0.0) {
    if (frequency > 0) {
      step = it / frequency;
      if (lastIt && (it % frequency != 0))
        ++step;
    }
  } else {
    step = (int)(t / frequency_dt);
    if (lastIt && (t<prtout))
      ++step;
  }

  return step;
}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::updatePrtout(double t)
{
  if(frequency_dt<=0.0)
    return;
  if(t>=prtout)
    prtout += frequency_dt;
}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::setMeshMotionHandler(IoData &ioData, MeshMotionHandler *mmh)
{

  rmmh = dynamic_cast<RigidMeshMotionHandler *>(mmh);

  if (ioData.problem.type[ProblemData::FORCED]) {
     if (ioData.forced.type == ForcedData::HEAVING)
       hmmh = dynamic_cast<HeavingMeshMotionHandler *>(mmh);
     else if (ioData.forced.type  == ForcedData::SPIRALING)
       smmh = dynamic_cast<SpiralingMeshMotionHandler *>(mmh);
     else if (ioData.forced.type  == ForcedData::PITCHING)
       pmmh = dynamic_cast<PitchingMeshMotionHandler *>(mmh);
     else if (ioData.forced.type  == ForcedData::DEFORMING)
       dmmh = dynamic_cast<DeformingMeshMotionHandler *>(mmh);
  }

  if (ioData.output.transient.generalizedforces[0] != 0 && ioData.problem.type[ProblemData::FORCED] && strcmp(ioData.input.strModesFile, "") == 0) {
    int sp = strlen(ioData.output.transient.prefix) + 1;
    generalizedforces = new char[sp + strlen(ioData.output.transient.generalizedforces)];
    sprintf(generalizedforces, "%s%s", ioData.output.transient.prefix, ioData.output.transient.generalizedforces);

    mX = new VecSet<DistSVec<double,3> >(1, domain->getNodeDistInfo());

    if (hmmh)
       (*mX)[0] = hmmh->getModes();
    else if(smmh)
       (*mX)[0] = smmh->getModes();
    else if(pmmh)
       (*mX)[0] = pmmh->getModes();
    else if(dmmh)
       (*mX)[0] = dmmh->getModes();

  }
//  else
//    generalizedforces = 0;

}

//------------------------------------------------------------------------------

template<int dim>
FILE *TsOutput<dim>::backupAsciiFile(char *name)
{

  FILE *fp = fopen(name, "r");
  if (!fp)
    return 0;

  char nameback[MAXLINE], line[MAXLINE];
  sprintf(nameback, "%s.back", name);
  FILE *fpback = fopen(nameback, "w");
  if (!fpback) {
    fprintf(stderr, "*** Error: could not open backup file \'%s\'\n", nameback);
    return 0;
  }
  while (fgets(line, MAXLINE, fp) != 0)
    fprintf(fpback, "%s", line);
  fclose(fp);
  fclose(fpback);

  fpback = fopen(nameback, "r");
  fp = fopen(name, "w");
  char* toto = fgets(line, MAXLINE, fpback);
  fprintf(fp, "%s", line);
  while (fgets(line, MAXLINE, fpback) != 0) {
    int iter;
    sscanf(line, "%d", &iter);
    if (iter <= it0)
      fprintf(fp, "%s", line);
  }
  fclose(fpback);

  return fp;

}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::openAsciiFiles()
{

  if (com->cpuNum() != 0) return;

  int nSurf = postOp->getNumSurf();
  int *surfNums = 0;
  int iSurf;

  if (nSurf > 0)  {
    surfNums = new int[nSurf];
    map<int, int> surfMap = postOp->getSurfMap();
    map<int, int>::iterator it;
    iSurf = 1;
    surfNums[0] = 0;
    for (it = surfMap.begin(); it != surfMap.end(); it++)
      if (it->second > 0)
        surfNums[iSurf++] = it->first;
  }

  if (forces) {
    if (it0 != 0)
      fpForces[0] = backupAsciiFile(forces);
    if (it0 == 0 || fpForces[0] == 0) {
      fpForces[0] = fopen(forces, "w");
      if (!fpForces[0]) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", forces);
        exit(1);
      }
      const char *addvar = "";
      if (rmmh) addvar = rmmh->getTagName();
      fprintf(fpForces[0], "#TimeIteration Time SubCycles NewtonSteps ");
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        fprintf(fpForces[0], "Cfx Cfy Cfz Cmx Cmy Cmz Energy %s\n", addvar);
      else
        fprintf(fpForces[0], "Fx Fy Fz Mx My Mz Energy %s\n", addvar);
    }
    fflush(fpForces[0]);
  }

  if (forces) {
    for (iSurf = 1; iSurf < nSurf; iSurf++) {
      char filename[256];
      sprintf(filename,"%s%d", forces, surfNums[iSurf]);
      if (it0 != 0)
        fpForces[iSurf] = backupAsciiFile(filename);
      if (it0 == 0 || fpForces[iSurf] == 0) {
        fpForces[iSurf] = fopen(filename, "w");
        if (!fpForces[iSurf]) {
           fprintf(stderr, "*** Error: could not open \'%s\'\n", filename);
           exit(1);
        }
        const char *addvar = "";
        if (rmmh) addvar = rmmh->getTagName();
        fprintf(fpForces[iSurf], "#TimeIteration Time SubCycles NewtonSteps ");
        if (refVal->mode == RefVal::NON_DIMENSIONAL)
          fprintf(fpForces[iSurf], "Cfx Cfy Cfz Cmx Cmy Cmz Energy %s\n", addvar);
        else
          fprintf(fpForces[iSurf], "Fx Fy Fz Mx My Mz Energy %s\n", addvar);
      }
      fflush(fpForces[iSurf]);
    }
  }

// Included (MB)
  if (switchOpt) {
    if (dForces) {
      fpdForces = fopen(dForces, "w");
      if (!fpdForces) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", dForces);
        exit(1);
      }
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        fprintf(fpdForces, "Step Variable Cfx Cfy Cfz Cmx Cmy Cmz sboom dCfx dCfy dCfz dCmx dCmy dCmz dSboom \n");
      else
        fprintf(fpdForces, "Step Variable Fx Fy Fz Mx My Mz sboom dFx dFy dFz dMx dMy dMz dSboom \n");

      fflush(fpdForces);
    }
    if (dMatchPressure) {
      fpdMatchPressure = fopen(dMatchPressure, "w");
      if (!fpdMatchPressure) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", dMatchPressure);
        exit(1);
      }
      fprintf(fpdMatchPressure, "Step Variable 0.5*||P-P_opt||_2^2 d[0.5*||P-P_opt||_2^2]\n");
      fflush(fpdMatchPressure);
    }
    if (dFluxNorm) {
      fpdFluxNorm = fopen(dFluxNorm, "w");
      if (!fpdFluxNorm) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", dFluxNorm);
        exit(1);
      }
      fprintf(fpdFluxNorm, "Step Variable normF dnormF\n");
      fflush(fpdFluxNorm);
    }
    if (dLiftDrag) {
      fpdLiftDrag = fopen(dLiftDrag, "w");
      if (!fpdLiftDrag) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", dLiftDrag);
        exit(1);
      }
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        fprintf(fpdLiftDrag, "Step Variable CLx CLy CLz dCLx dCLy dCLz \n");
      else
        fprintf(fpdLiftDrag, "Step Variable Lx Ly Lz dLx dLy dLz \n");

      fflush(fpdLiftDrag);
    }
    if (dLiftx) {
      fpdLiftx = fopen(dLiftx, "w");
      if (!fpdLiftx) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", dLiftx);
        exit(1);
      }
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        fprintf(fpdLiftx, "dCLx\n");
      else
        fprintf(fpdLiftx, "dLx\n");

      fflush(fpdLiftx);
    }
    if (dLifty) {
      fpdLifty = fopen(dLifty, "w");
      if (!fpdLifty) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", dLifty);
        exit(1);
      }
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        fprintf(fpdLifty, "dCLy\n");
      else
        fprintf(fpdLifty, "dLy\n");

      fflush(fpdLifty);
    }
    if (dLiftz) {
      fpdLiftz = fopen(dLiftz, "w");
      if (!fpdLiftz) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", dLiftz);
        exit(1);
      }
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        fprintf(fpdLiftz, "dCLz\n");
      else
        fprintf(fpdLiftz, "dLz\n");

      fflush(fpdLiftz);
    }
  }
  if (hydrostaticforces) {
    if (it0 != 0){
      fpHydroStaticForces[0] = backupAsciiFile(hydrostaticforces);
    }
    if (it0 == 0 || fpHydroStaticForces[0] == 0) {
      fpHydroStaticForces[0] = fopen(hydrostaticforces, "w");
      if (!fpHydroStaticForces[0]) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", hydrostaticforces);
        exit(1);
      }
      const char *addvar = "";
      if (rmmh) addvar = rmmh->getTagName();
      fprintf(fpHydroStaticForces[0], "#TimeIteration Time SubCycles NewtonSteps ");
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        fprintf(fpHydroStaticForces[0], "Cfx Cfy Cfz Cmx Cmy Cmz Energy %s\n", addvar);
      else
        fprintf(fpHydroStaticForces[0], "Fx Fy Fz Mx My Mz Energy %s\n", addvar);
    }
    fflush(fpHydroStaticForces[0]);
  }

  if (hydrostaticforces) {
    for (iSurf = 1; iSurf < nSurf; iSurf++) {
      char filename[256];
      sprintf(filename,"%s%d", hydrostaticforces, surfNums[iSurf]);
      if (it0 != 0)
        fpHydroStaticForces[iSurf] = backupAsciiFile(filename);
      if (it0 == 0 || fpHydroStaticForces[iSurf] == 0) {
        fpHydroStaticForces[iSurf] = fopen(filename, "w");
        if (!fpHydroStaticForces[iSurf]) {
           fprintf(stderr, "*** Error: could not open \'%s\'\n", filename);
           exit(1);
        }
        const char *addvar = "";
        if (rmmh) addvar = rmmh->getTagName();
        fprintf(fpHydroStaticForces[iSurf], "#TimeIteration Time SubCycles NewtonSteps ");
        if (refVal->mode == RefVal::NON_DIMENSIONAL)
          fprintf(fpHydroStaticForces[iSurf], "Cfx Cfy Cfz Cmx Cmy Cmz Energy %s\n", addvar);
        else
          fprintf(fpHydroStaticForces[iSurf], "Fx Fy Fz Mx My Mz Energy %s\n", addvar);
      }
      fflush(fpHydroStaticForces[iSurf]);
    }
  }

  if (hydrodynamicforces) {
    if (it0 != 0)
      fpHydroDynamicForces[0] = backupAsciiFile(hydrodynamicforces);
    if (it0 == 0 || fpHydroDynamicForces[0] == 0) {
      fpHydroDynamicForces[0] = fopen(hydrodynamicforces, "w");
      if (!fpHydroDynamicForces[0]) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", hydrodynamicforces);
        exit(1);
      }
      const char *addvar = "";
      if (rmmh) addvar = rmmh->getTagName();
      fprintf(fpHydroDynamicForces[0], "#TimeIteration Time SubCycles NewtonSteps ");
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        fprintf(fpHydroDynamicForces[0], "Cfx Cfy Cfz Cmx Cmy Cmz Energy %s\n", addvar);
      else
        fprintf(fpHydroDynamicForces[0], "Fx Fy Fz Mx My Mz Energy %s\n", addvar);
    }
    fflush(fpHydroDynamicForces[0]);
  }

  if (hydrodynamicforces) {
    for (iSurf = 1; iSurf < nSurf; iSurf++) {
      char filename[256];
      sprintf(filename,"%s%d", hydrodynamicforces, surfNums[iSurf]);
      if (it0 != 0)
        fpHydroDynamicForces[iSurf] = backupAsciiFile(filename);
      if (it0 == 0 || fpHydroDynamicForces[iSurf] == 0) {
        fpHydroDynamicForces[iSurf] = fopen(filename, "w");
        if (!fpHydroDynamicForces[iSurf]) {
           fprintf(stderr, "*** Error: could not open \'%s\'\n", filename);
           exit(1);
        }
        const char *addvar = "";
        if (rmmh) addvar = rmmh->getTagName();
        fprintf(fpHydroDynamicForces[iSurf], "#TimeIteration Time SubCycles NewtonSteps ");
        if (refVal->mode == RefVal::NON_DIMENSIONAL)
          fprintf(fpHydroDynamicForces[iSurf], "Cfx Cfy Cfz Cmx Cmy Cmz Energy %s\n", addvar);
        else
          fprintf(fpHydroDynamicForces[iSurf], "Fx Fy Fz Mx My Mz Energy %s\n", addvar);
      }
      fflush(fpHydroDynamicForces[iSurf]);
    }
  }

  if (tavforces) {
    if (it0 != 0)
      fpTavForces[0] = backupAsciiFile(tavforces);
    if (it0 == 0 || fpTavForces[0] == 0) {
      fpTavForces[0] = fopen(tavforces, "w");
      if (!fpTavForces[0]) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", tavforces);
        exit(1);
      }
      const char *addvar = "";
      if (rmmh) addvar = rmmh->getTagName();
      fprintf(fpTavForces[0], "#TimeIteration Time SubCycles NewtonSteps ");
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        fprintf(fpTavForces[0], "Cfx Cfy Cfz Cmx Cmy Cmz Energy %s\n", addvar);
      else
        fprintf(fpTavForces[0], "Fx Fy Fz Mx My Mz Energy %s\n", addvar);
    }
    fflush(fpTavForces[0]);
  }

  if (tavforces) {
    for (iSurf = 1; iSurf < nSurf; iSurf++)  {
      char filename[256];
      sprintf(filename,"%s%d", tavforces, surfNums[iSurf]);
      if (it0 != 0)
        fpTavForces[iSurf] = backupAsciiFile(filename);
      if (it0 == 0 || fpTavForces[iSurf] == 0) {
        fpTavForces[iSurf] = fopen(filename, "w");
        if (!fpTavForces[iSurf]) {
          fprintf(stderr, "*** Error: could not open \'%s\'\n", filename);
          exit(1);
        }
        const char *addvar = "";
        if (rmmh) addvar = rmmh->getTagName();
        fprintf(fpTavForces[iSurf], "#TimeIteration Time SubCycles NewtonSteps ");
        if (refVal->mode == RefVal::NON_DIMENSIONAL)
          fprintf(fpTavForces[iSurf], "Cfx Cfy Cfz Cmx Cmy Cmz Energy %s\n", addvar);
        else
          fprintf(fpTavForces[iSurf], "Fx Fy Fz Mx My Mz Energy %s\n", addvar);
      }
      fflush(fpTavForces[iSurf]);
    }
  }


  if (generalizedforces) {
    if (it0 != 0)
      fpGnForces = backupAsciiFile(generalizedforces);
    if (it0 == 0 || fpGnForces == 0) {
      fpGnForces = fopen(generalizedforces, "w");
      if (!fpGnForces) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", generalizedforces);
        exit(1);
      }
      const char *addvar = "";
      if (rmmh) addvar = rmmh->getTagName();
      fprintf(fpGnForces, "#TimeIteration Time SubCycles NewtonSteps ");
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        fprintf(fpGnForces, "Coefficient of Generalized Forces %s\n", addvar);
      else
        fprintf(fpGnForces, "Generalized Forces %s\n", addvar);
    }
    fflush(fpGnForces);
  }

  if (generalizedforces) {
    for (iSurf = 1; iSurf < nSurf; iSurf++) {
      char filename[256];
      sprintf(filename,"%s%d", generalizedforces, surfNums[iSurf]);
      if (it0 != 0)
        fpGnForces = backupAsciiFile(filename);
      if (it0 == 0 || fpGnForces == 0) {
        fpGnForces = fopen(filename, "w");
        if (!fpGnForces) {
           fprintf(stderr, "*** Error: could not open \'%s\'\n", filename);
           exit(1);
        }
        const char *addvar = "";
        if (rmmh) addvar = rmmh->getTagName();
        fprintf(fpGnForces, "#TimeIteration Time SubCycles NewtonSteps ");
        if (refVal->mode == RefVal::NON_DIMENSIONAL)
          fprintf(fpGnForces, "Coefficient of Generalized Forces %s\n", addvar);
        else
          fprintf(fpGnForces, "Generalized Forces %s\n", addvar);
      }
      fflush(fpGnForces);
    }
  }

  if (lift) {
   if(it0 !=0)
      fpLift[0] = backupAsciiFile(lift);
   if (it0 == 0 || fpLift[0] == 0){
     fpLift[0] = fopen(lift,"w");
     if(!fpLift[0]){
       fprintf(stderr, "*** Error: could not open \'%s\'\n", lift);
       exit(1);
     }
     const char *addvar = "";
     if(rmmh) addvar = rmmh->getTagName();
     fprintf(fpLift[0], "#TimeIteration Time SubCycles NewtonSteps Lx Ly Lz %s\n", addvar);
     }
   fflush(fpLift[0]);
  }

  if (lift) {
    for (iSurf = 1; iSurf < nSurf; iSurf++)  {
      char filename[256];
      sprintf(filename,"%s%d", lift, surfNums[iSurf]);
      if(it0 !=0)
        fpLift[iSurf] = backupAsciiFile(filename);
      if (it0 == 0 || fpLift[iSurf] == 0){
        fpLift[iSurf] = fopen(filename,"w");
        if (!fpLift[iSurf]) {
          fprintf(stderr, "*** Error: could not open \'%s\'\n", filename);
          exit(1);
        }
        const char *addvar = "";
        if(rmmh) addvar = rmmh->getTagName();
        fprintf(fpLift[iSurf], "#TimeIteration Time SubCycles NewtonSteps Lx Ly Lz %s\n", addvar);
      }
      fflush(fpLift[iSurf]);
    }
  }

  if (hydrostaticlift) {
   if(it0 !=0)
      fpHydroStaticLift[0] = backupAsciiFile(hydrostaticlift);
   if (it0 == 0 || fpHydroStaticLift[0] == 0){
     fpHydroStaticLift[0] = fopen(hydrostaticlift,"w");
     if(!fpHydroStaticLift[0]){
       fprintf(stderr, "*** Error: could not open \'%s\'\n", hydrostaticlift);
       exit(1);
     }
     const char *addvar = "";
     if(rmmh) addvar = rmmh->getTagName();
     fprintf(fpHydroStaticLift[0], "#TimeIteration Time SubCycles NewtonSteps Lx Ly Lz %s\n", addvar);
     }
   fflush(fpHydroStaticLift[0]);
  }

  if (hydrostaticlift) {
    for (iSurf = 1; iSurf < nSurf; iSurf++)  {
      char filename[256];
      sprintf(filename,"%s%d", hydrostaticlift, surfNums[iSurf]);
      if(it0 !=0)
        fpHydroStaticLift[iSurf] = backupAsciiFile(filename);
      if (it0 == 0 || fpHydroStaticLift[iSurf] == 0){
        fpHydroStaticLift[iSurf] = fopen(filename,"w");
        if (!fpHydroStaticLift[iSurf]) {
          fprintf(stderr, "*** Error: could not open \'%s\'\n", filename);
          exit(1);
        }
        const char *addvar = "";
        if(rmmh) addvar = rmmh->getTagName();
        fprintf(fpHydroStaticLift[iSurf], "#TimeIteration Time SubCycles NewtonSteps Lx Ly Lz %s\n", addvar);
      }
      fflush(fpHydroStaticLift[iSurf]);
    }
  }

  if (hydrodynamiclift) {
   if(it0 !=0)
      fpHydroDynamicLift[0] = backupAsciiFile(hydrodynamiclift);
   if (it0 == 0 || fpHydroDynamicLift[0] == 0){
     fpHydroDynamicLift[0] = fopen(hydrodynamiclift,"w");
     if(!fpHydroDynamicLift[0]){
       fprintf(stderr, "*** Error: could not open \'%s\'\n", hydrodynamiclift);
       exit(1);
     }
     const char *addvar = "";
     if(rmmh) addvar = rmmh->getTagName();
     fprintf(fpHydroDynamicLift[0], "#TimeIteration Time SubCycles NewtonSteps Lx Ly Lz %s\n", addvar);
     }
   fflush(fpHydroDynamicLift[0]);
  }

  if (hydrodynamiclift) {
    for (iSurf = 1; iSurf < nSurf; iSurf++)  {
      char filename[256];
      sprintf(filename,"%s%d", hydrodynamiclift, surfNums[iSurf]);
      if(it0 !=0)
        fpHydroDynamicLift[iSurf] = backupAsciiFile(filename);
      if (it0 == 0 || fpHydroDynamicLift[iSurf] == 0){
        fpHydroDynamicLift[iSurf] = fopen(filename,"w");
        if (!fpHydroDynamicLift[iSurf]) {
          fprintf(stderr, "*** Error: could not open \'%s\'\n", filename);
          exit(1);
        }
        const char *addvar = "";
        if(rmmh) addvar = rmmh->getTagName();
        fprintf(fpHydroDynamicLift[iSurf], "#TimeIteration Time SubCycles NewtonSteps Lx Ly Lz %s\n", addvar);
      }
      fflush(fpHydroDynamicLift[iSurf]);
    }
  }

  if(tavlift) {
    if (it0 != 0)
      fpTavLift[0] = backupAsciiFile(tavlift);
   if (it0 == 0 || fpTavLift[0] == 0){
     fpTavLift[0] = fopen(tavlift,"w");
     if(!fpTavLift[0]){
       fprintf(stderr, "*** Error: could not open \'%s\'\n", tavlift);
       exit(1);
     }
     const char *addvar = "";
     if(rmmh) addvar = rmmh->getTagName();
     fprintf(fpTavLift[0], "#TimeIteration Time SubCycles NewtonSteps Lx Ly Lz %s\n", addvar);
   }
   fflush(fpTavLift[0]);
  }

  if(tavlift) {
    for (iSurf = 1; iSurf < nSurf; iSurf++)  {
      char filename[256];
      sprintf(filename,"%s%d", tavlift, surfNums[iSurf]);
      if(it0 !=0)
        fpTavLift[iSurf] = backupAsciiFile(filename);
      if (it0 == 0 || fpTavLift[iSurf] == 0){
        fpTavLift[iSurf] = fopen(filename,"w");
        if (!fpTavLift[iSurf]) {
          fprintf(stderr, "*** Error: could not open \'%s\'\n", filename);
          exit(1);
        }
        const char *addvar = "";
        if(rmmh) addvar = rmmh->getTagName();
        fprintf(fpTavLift[iSurf], "#TimeIteration Time SubCycles NewtonSteps Lx Ly Lz %s\n", addvar);
      }
      fflush(fpTavLift[iSurf]);
    }
  }

//For Heat Fluxes the flags are different

  int nSurfHF = postOp->getNumSurfHF();
   int *surfNumsHF = 0;
    if (nSurfHF > 0)  {
    surfNumsHF = new int[nSurfHF];
    map<int,int> surfMapHF = postOp->getSurfMapHF();
    map<int, int>::iterator it;
    iSurf = 1;
    surfNumsHF[0] = 0;
    for (it = surfMapHF.begin(); it != surfMapHF.end(); it++)
      if (it->second > 0)
        surfNumsHF[iSurf++] = it->first;
  }

  if (heatfluxes) {
    if (it0 != 0)
      fpHeatFluxes[0] = backupAsciiFile(heatfluxes);
    if (it0 == 0 || fpHeatFluxes[0] == 0) {
      fpHeatFluxes[0] = fopen(heatfluxes, "w");
      if (!fpHeatFluxes[0]) {
        fprintf(stderr, "*** HF Error: could not open \'%s\'\n", heatfluxes);
        exit(1);
      }
      const char *addvar = "";
 //     if (rmmh) addvar = rmmh->getTagName();
      fprintf(fpHeatFluxes[0], "#TimeIteration Time SubCycles NewtonSteps ");
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        fprintf(fpHeatFluxes[0], "Nondimensional HeatFlux %s\n", addvar);
      else
        fprintf(fpHeatFluxes[0], "HeatFlux %s\n", addvar);
    }
    fflush(fpHeatFluxes[0]);
  }

  if (heatfluxes) {
    for (iSurf = 1; iSurf < nSurfHF; iSurf++) {
      char filename[256];
      sprintf(filename,"%s%d", heatfluxes, surfNumsHF[iSurf]);
      if (it0 != 0)
        fpHeatFluxes[iSurf] = backupAsciiFile(filename);
      if (it0 == 0 || fpHeatFluxes[iSurf] == 0) {
        fpHeatFluxes[iSurf] = fopen(filename, "w");
        if (!fpHeatFluxes[iSurf]) {
           fprintf(stderr, "*** HF Error: could not open \'%s\'\n", filename);
           exit(1);
        }
        const char *addvar = "";
        if (rmmh) addvar = rmmh->getTagName();
        fprintf(fpHeatFluxes[iSurf], "#TimeIteration Time SubCycles NewtonSteps ");
        if (refVal->mode == RefVal::NON_DIMENSIONAL)
          fprintf(fpHeatFluxes[iSurf], "Nondimensional HeatFlux  %s\n", addvar);
        else
          fprintf(fpHeatFluxes[iSurf], "HeatFlux %s\n", addvar);
      }
      fflush(fpHeatFluxes[iSurf]);
    }
  }

  if (residuals) {
    if (it0 != 0)
      fpResiduals = backupAsciiFile(residuals);
    if (it0 == 0 || fpResiduals == 0) {
      fpResiduals = fopen(residuals, "w");
      if (!fpResiduals) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", residuals);
        exit(1);
      }
      fprintf(fpResiduals, "#TimeIteration ElapsedTime RelativeResidual CflNumber\n");
    }
    fflush(fpResiduals);
  }

  if (matchpressure) {
    if (it0 != 0)
      fpMatchPressure = backupAsciiFile(matchpressure);
    if (it0 == 0 || fpMatchPressure == 0) {
      fpMatchPressure = fopen(matchpressure, "w");
      if (!fpMatchPressure) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", matchpressure);
        exit(1);
      }
      fprintf(fpMatchPressure, "#TimeIteration Time SubCycles NewtonSteps 0.5||P-P^*||_2^2/||P^*||_2^2\n");
    }
    fflush(fpMatchPressure);
  }

  if (matchstate) {
    if (it0 != 0)
      fpMatchState = backupAsciiFile(matchstate);
    if (it0 == 0 || fpMatchState == 0) {
      fpMatchState = fopen(matchstate, "w");
      if (!fpMatchState) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", matchstate);
        exit(1);
      }
      fprintf(fpMatchState, "#TimeIteration Time ElapsedTime ||U-U_ref||_2/||U_ref||_2\n");
    }
    fflush(fpMatchState);
  }

   if (fluxnorm) {
    if (it0 != 0)
      fpFluxNorm = backupAsciiFile(fluxnorm);
    if (it0 == 0 || fpFluxNorm == 0) {
      fpFluxNorm = fopen(fluxnorm, "w");
      if (!fpFluxNorm) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", fluxnorm);
        exit(1);
      }
      fprintf(fpFluxNorm, "#TimeIteration Time SubCycles NewtonSteps 0.5||R||_2^2\n");
    }
    fflush(fpFluxNorm);
  }

  if (material_volumes) {
    if (it0 != 0)
      fpMatVolumes = backupAsciiFile(material_volumes);
    if (it0 == 0 || fpMatVolumes == 0) {
      fpMatVolumes = fopen(material_volumes, "w");
      if (!fpMatVolumes) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", material_volumes);
        exit(1);
      }
      fprintf(fpMatVolumes, "#TimeIteration ElapsedTime ");
      for(int i=0; i<numFluidPhases; i++)
        fprintf(fpMatVolumes, "Volume[FluidID==%d] ", i);
      fprintf(fpMatVolumes, "Volume[FluidID==%d(GhostSolid)] TotalVolume\n", numFluidPhases);
    }
    fflush(fpMatVolumes);
  }

  if (material_mass_energy) {
    if (it0 != 0)
      fpMaterialMassEnergy = backupAsciiFile(material_mass_energy);
    if (it0 == 0 || fpMaterialMassEnergy == 0) {
      fpMaterialMassEnergy = fopen(material_mass_energy, "w");
      if (!fpMaterialMassEnergy) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", material_mass_energy);
        exit(1);
      }
      fprintf(fpMaterialMassEnergy, "#TimeIteration ElapsedTime ");
      for(int i=0; i<numFluidPhases; i++)
        fprintf(fpMaterialMassEnergy, "Mass[FluidID==%d] Energy[FluidId==%d]", i,i);
      fprintf(fpMaterialMassEnergy, "Mass[FluidID==%d(GhostSolid)] Energy[FluidID==%d(GhostSolid)] TotalVolume\n", numFluidPhases, numFluidPhases);
    }
    fflush(fpMaterialMassEnergy);
  }







  if (embeddedsurface) {
    if (it0 != 0)
      fpEmbeddedSurface = backupAsciiFile(embeddedsurface);
    if (it0 == 0 || fpEmbeddedSurface == 0) {
      fpEmbeddedSurface = fopen(embeddedsurface, "w");
      if (!fpEmbeddedSurface) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", embeddedsurface);
        exit(1);
      }
      fprintf(fpEmbeddedSurface, "Vector Displacement under NLDynamic for EmbeddedNodes\n");
    }
    fflush(fpEmbeddedSurface);
  }

  if (cputiming) {
    if (it0 != 0)
      fpCpuTiming = backupAsciiFile(cputiming);
    if (it0 == 0 || fpCpuTiming== 0) {
      fpCpuTiming= fopen(cputiming, "w");
      if (!fpCpuTiming) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", cputiming);
        exit(1);
      }
    }
    fflush(fpCpuTiming);
  }

  if (conservation) {
    if (it0 != 0)
      fpConservationErr = backupAsciiFile(conservation);
    if (it0 == 0 || fpConservationErr == 0) {
      fpConservationErr = fopen(conservation, "w");
      if (!fpConservationErr) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", conservation);
        exit(1);
      }
      fprintf(fpConservationErr, "# 1TimeIteration 2Time ");
      fprintf(fpConservationErr, "3TotalExpectedMass 4TotalExpectedMomentumx 5TotalExpectedMomentumy 6TotalExpectedMomentumz 7TotalExpectedEnergy ");
      fprintf(fpConservationErr, "8Fluid1ExpectedMass 9Fluid1ExpectedMomentumx 10Fluid1ExpectedMomentumx 11Fluid1ExpectedMomentumz 12Fluid1ExpectedEnergy ");
      fprintf(fpConservationErr, "13Fluid2ExpectedMass 14Fluid2ExpectedMomentumx 15Fluid2ExpectedMomentumx 16Fluid2ExpectedMomentumz 17Fluid2ExpectedEnergy ");
      fprintf(fpConservationErr, "18TotalComputedMass 19TotalComputedMomentumx 20TotalComputedMomentumy 21TotalComputedMomentumz 22TotalComputedEnergy ");
      fprintf(fpConservationErr, "23Fluid1ComputedMass 24Fluid1ComputedMomentumx 25Fluid1ComputedMomentumx 26Fluid1ComputedMomentumz 27Fluid1ComputedEnergy ");
      fprintf(fpConservationErr, "28Fluid2ComputedMass 29Fluid2ComputedMomentumx 30Fluid2ComputedMomentumx 31Fluid2ComputedMomentumz 32Fluid2ComputedEnergy\n");
    }
    fflush(fpConservationErr);
 }


  if (embeddedsurfaceCp) {

    if (it0 != 0)
      fpEmbeddedSurfaceCp = backupAsciiFile(embeddedsurfaceCp);
    if (it0 == 0 || fpEmbeddedSurfaceCp == 0) {
      fpEmbeddedSurfaceCp = fopen(embeddedsurfaceCp, "w");
      if (!fpEmbeddedSurfaceCp) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", embeddedsurfaceCp);
        exit(1);
      }
      fprintf(fpEmbeddedSurfaceCp, "Scalar %s under load for FluidNodes\n", embeddedsurfaceCp);
    }
    fflush(fpEmbeddedSurfaceCp);

  }

  if (embeddedsurfaceCf) {

    if (it0 != 0)
      fpEmbeddedSurfaceCf = backupAsciiFile(embeddedsurfaceCf);
    if (it0 == 0 || fpEmbeddedSurfaceCf == 0) {
      fpEmbeddedSurfaceCf = fopen(embeddedsurfaceCf, "w");
      if (!fpEmbeddedSurfaceCf) {
        fprintf(stderr, "*** Error: could not open \'%s\'\n", embeddedsurfaceCf);
        exit(1);
      }
      fprintf(fpEmbeddedSurfaceCf, "Scalar %s under load for FluidNodes\n", embeddedsurfaceCf);
    }
    fflush(fpEmbeddedSurfaceCf);

  }


 delete [] surfNums;
 delete [] surfNumsHF;
}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::closeAsciiFiles()
{

  for (int iSurf = 0; iSurf < postOp->getNumSurf(); iSurf++)  {

    if (fpForces[iSurf]) fclose(fpForces[iSurf]);
    if (fpHydroDynamicForces[iSurf]) fclose(fpHydroDynamicForces[iSurf]);
    if (fpHydroStaticForces[iSurf]) fclose(fpHydroStaticForces[iSurf]);
    if (fpLift[iSurf]) { fclose(fpLift[iSurf]); fpLift[iSurf] = 0; }
    if (fpTavForces[iSurf]) fclose(fpTavForces[iSurf]);
    if (fpTavLift[iSurf]) fclose(fpTavLift[iSurf]);

  }

  for (int iSurf = 0; iSurf < postOp->getNumSurfHF(); iSurf++)  {
     if (fpHeatFluxes[iSurf]) fclose(fpHeatFluxes[iSurf]);
  }


  if (fpResiduals) fclose(fpResiduals);
  if (fpMatchPressure) fclose(fpMatchPressure);
  if (fpMatchState) fclose(fpMatchState);
  if (fpFluxNorm) fclose(fpFluxNorm);
  if (fpMatVolumes) fclose(fpMatVolumes);
  if (fpEmbeddedSurface) fclose(fpEmbeddedSurface);
  if (fpCpuTiming) fclose(fpCpuTiming);
  if (fpGnForces) fclose(fpGnForces);
  if (fpConservationErr) fclose(fpConservationErr);

  if (fpEmbeddedSurfaceCp) fclose(fpEmbeddedSurfaceCp);
  if (fpEmbeddedSurfaceCf) fclose(fpEmbeddedSurfaceCf);

}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeForcesToDisk(DistExactRiemannSolver<dim> &riemann,
                                      bool lastIt, int it, int itSc, int itNl, double t, double cpu,
                                      double* e, DistSVec<double,3> &X, DistSVec<double,dim> &U,
                                      DistVec<int> *fluidId)
{

  int nSurfs = postOp->getNumSurf();

  Vec3D *Fi = new Vec3D[nSurfs];
  Vec3D *Mi = new Vec3D[nSurfs];
  Vec3D *Fv = new Vec3D[nSurfs];
  Vec3D *Mv = new Vec3D[nSurfs];

  Vec<double> GF( mX ? mX->numVectors(): 0);
  if(mX) GF = 0.0;

  double del_t;

  if (TavF ==0 && TavM == 0) {
    TavF = new Vec3D[nSurfs];
    TavM = new Vec3D[nSurfs];
    for(int i = 0; i < nSurfs; ++i) {
      TavF[i] = 0.0;
      TavM[i] = 0.0;
    }
  }

  if(counter == 0) {
    tinit = refVal->time*t;
    tprevf = refVal->time*t;
    tener = 0.0; tenerold = 0.0;
    del_t = 0.0;
  }

  double time = refVal->time * t;

  if (forces || tavforces)  {
    Vec3D rVec = x0;
    // if single-phase flow -- fluidId is a null pointer
    // if multi-phase flow  -- fluidId points to a DistVec<int>
    postOp->computeForceAndMoment(riemann, rVec, X, U, fluidId, Fi, Mi, Fv, Mv, 0, mX, mX ? &GF: 0);
    if(mX != 0)
      com->globalSum(GF.size(), GF.data());

    Vec3D F = Fi[0] + Fv[0];
    Vec3D moment = Mi[0] + Mv[0];
    F *= refVal->force;
    moment *= refVal->energy;
  }

  int iSurf;

  if (fpForces[0] || fpTavForces[0]) {

    for (iSurf = 0; iSurf < nSurfs; iSurf++)  {
      Vec3D F = Fi[iSurf] + Fv[iSurf];
      Vec3D M = Mi[iSurf] + Mv[iSurf];
      if (refVal->mode == RefVal::NON_DIMENSIONAL) {
        F *= 2.0 * refVal->length*refVal->length / surface;
        M *= 2.0 * refVal->length*refVal->length*refVal->length / (surface * length);
      }
      else {
        F *= refVal->force;
        M *= refVal->energy;
      }
      double energy = refVal->energy * e[0];


      if (rmmh) {
        double tag = rmmh->getTagValue(t);
        fprintf(fpForces[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e \n",
                it, time, itSc, itNl, F[0], F[1], F[2], M[0], M[1], M[2], energy, tag);
      }
      else
        fprintf(fpForces[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e \n",
                it, time, itSc, itNl, F[0], F[1], F[2], M[0], M[1], M[2], energy);

      fflush(fpForces[iSurf]);

      if(fpTavForces[iSurf] && counter == 0){
        fprintf(fpTavForces[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e \n",
          it, time, itSc, itNl, F[0], F[1], F[2], M[0], M[1], M[2], energy);
      }

      del_t = time - tprevf;
      F *= del_t; TavF[iSurf] += F;
      M *= del_t; TavM[iSurf] += M;
      tener += energy*del_t;

      if(fpTavForces[iSurf] && counter > 0){
        F = TavF[iSurf]; M = TavM[iSurf]; energy = tener;
        F /= (time - tinit); M /= (time - tinit);
        energy /= (time - tinit);
        fprintf(fpTavForces[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e \n",
           it, time, itSc, itNl, F[0], F[1], F[2], M[0], M[1], M[2], energy);
      }

      fflush(fpTavForces[iSurf]);
    }
  }
  if (fpGnForces) {

     fprintf(fpGnForces,"%d %e %d %d", it, time, itSc, itNl);
     for (int i=0; i < mX->numVectors(); i++) {
       if (refVal->mode == RefVal::NON_DIMENSIONAL)
         fprintf(fpGnForces," %e", 2.0*refVal->length*refVal->length*GF[i]/surface);
       else
         fprintf(fpGnForces," %e", refVal->force*GF[i]);
     }

      fprintf(fpGnForces, "\n");

      fflush(fpGnForces);
  }

  if (!steady) {
    if (rmmh) {
      double tag = rmmh->getTagValue(t);
      com->printf(0, "It %d (%d,%d): Time=%g, Mach=%g, Elapsed Time=%.2e s\n",
                  it, itSc, itNl, t*refVal->time, tag, cpu);
    }
    else
      com->printf(0, "It %d (%d,%d): Time=%g, Elapsed Time=%.2e s\n",
                  it, itSc, itNl, t*refVal->time, cpu);
  }
  tprevf = time;

  delete[] Fi;
  delete[] Fv;
  delete[] Mi;
  delete[] Mv;
}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeForcesToDisk(bool lastIt, int it, int itSc, int itNl, double t, double cpu,
                                      double* e, DistSVec<double,3> &X, DistSVec<double,dim> &U,
                                      DistVec<int> *fluidId)
{
  int nSurfs = postOp->getNumSurf();

  Vec3D *Fi = new Vec3D[nSurfs];
  Vec3D *Mi = new Vec3D[nSurfs];
  Vec3D *Fv = new Vec3D[nSurfs];
  Vec3D *Mv = new Vec3D[nSurfs];

  Vec<double> GF( mX ? mX->numVectors(): 0);
  if(mX) GF = 0.0;

  double del_t;

  if (TavF ==0 && TavM == 0) {
    TavF = new Vec3D[nSurfs];
    TavM = new Vec3D[nSurfs];
    for(int i = 0; i < nSurfs; ++i) {
      TavF[i] = 0.0;
      TavM[i] = 0.0;
    }
  }

  if(counter == 0) {
    tinit = refVal->time*t;
    tprevf = refVal->time*t;
    tener = 0.0; tenerold = 0.0;
    del_t = 0.0;
  }

  double time = refVal->time * t;

  if (forces || tavforces)  {
    Vec3D rVec = x0;
    // if single-phase flow -- fluidId is a null pointer
    // if multi-phase flow  -- fluidId points to a DistVec<int>
    postOp->computeForceAndMoment(rVec, X, U, fluidId, Fi, Mi, Fv, Mv, 0, mX, mX ? &GF: 0);
    if(mX != 0)
      com->globalSum(GF.size(), GF.data());

    Vec3D F = Fi[0] + Fv[0];
    Vec3D moment = Mi[0] + Mv[0];
    F *= refVal->force;
    moment *= refVal->energy;
  }

  int iSurf;

  if (fpForces[0] || fpTavForces[0]) {

    for (iSurf = 0; iSurf < nSurfs; iSurf++)  {
      Vec3D F = Fi[iSurf] + Fv[iSurf];
      Vec3D M = Mi[iSurf] + Mv[iSurf];
      if (refVal->mode == RefVal::NON_DIMENSIONAL) {
        F *= 2.0 * refVal->length*refVal->length / surface;
        M *= 2.0 * refVal->length*refVal->length*refVal->length / (surface * length);
      }
      else {
        F *= refVal->force;
        M *= refVal->energy;
      }
      double energy = refVal->energy * e[0];


      if (rmmh) {
        double tag = rmmh->getTagValue(t);
        fprintf(fpForces[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e \n",
                it, time, itSc, itNl, F[0], F[1], F[2], M[0], M[1], M[2], energy, tag);
      }
      else
        fprintf(fpForces[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e \n",
                it, time, itSc, itNl, F[0], F[1], F[2], M[0], M[1], M[2], energy);

      fflush(fpForces[iSurf]);

      if(fpTavForces[iSurf] && counter == 0){
        fprintf(fpTavForces[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e \n",
          it, time, itSc, itNl, F[0], F[1], F[2], M[0], M[1], M[2], energy);
      }

      del_t = time - tprevf;
      F *= del_t; TavF[iSurf] += F;
      M *= del_t; TavM[iSurf] += M;
      tener += energy*del_t;

      if(fpTavForces[iSurf] && counter > 0){
        F = TavF[iSurf]; M = TavM[iSurf]; energy = tener;
        F /= (time - tinit); M /= (time - tinit);
        energy /= (time - tinit);
        fprintf(fpTavForces[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e \n",
           it, time, itSc, itNl, F[0], F[1], F[2], M[0], M[1], M[2], energy);
      }

      fflush(fpTavForces[iSurf]);
    }
  }
  if (fpGnForces) {

     fprintf(fpGnForces,"%d %e %d %d", it, time, itSc, itNl);
     for (int i=0; i < mX->numVectors(); i++) {
       if (refVal->mode == RefVal::NON_DIMENSIONAL)
         fprintf(fpGnForces," %e", 2.0*refVal->length*refVal->length*GF[i]/surface);
       else
         fprintf(fpGnForces," %e", refVal->force*GF[i]);
     }

      fprintf(fpGnForces, "\n");

      fflush(fpGnForces);
  }

  if (!steady) {
    if (rmmh) {
      double tag = rmmh->getTagValue(t);
      com->printf(0, "It %d (%d,%d): Time=%g, Mach=%g, Elapsed Time=%.2e s\n",
                  it, itSc, itNl, t*refVal->time, tag, cpu);
    }
    else
      com->printf(0, "It %d (%d,%d): Time=%g, Elapsed Time=%.2e s\n",
                  it, itSc, itNl, t*refVal->time, cpu);
  }
  tprevf = time;

  delete[] Fi;
  delete[] Fv;
  delete[] Mi;
  delete[] Mv;
}

//------------------------------------------------------------------------------

// Included (MZ)
template<int dim>
void TsOutput<dim>::writeDerivativeOfFluxNormToDisk(int it, int actvar, double normF, double dnormF)
{

  if (fpdFluxNorm) {
    fprintf(fpdFluxNorm, "%d %d %16.13e %16.13e \n", it, actvar, normF, dnormF);
    fflush(fpdFluxNorm);
  }

}


//------------------------------------------------------------------------------

// Included (YC)
template<int dim>
void TsOutput<dim>::writeDerivativeOfLiftDragToDisk(int it, int actvar, Vec3D & L, Vec3D & dL)
{

  if (fpdLiftDrag) {
    fprintf(fpdLiftDrag, "%d %d %16.13e %16.13e %16.13e %16.13e %16.13e %16.13e \n", it, actvar, L[0], L[1], L[2], dL[0], dL[1], dL[2]);
    fflush(fpdLiftDrag);
  }

}



//------------------------------------------------------------------------------
template<int dim>
void TsOutput<dim>::writeDerivativeOfLiftxToDisk(double& dLx)
{

  if (fpdLiftx) {
    fprintf(fpdLiftx, "%16.13e\n", dLx);
    fflush(fpdLiftx);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeDerivativeOfLiftyToDisk(double& dLy)
{

  if (fpdLifty) {
    fprintf(fpdLifty, "%16.13e\n", dLy);
    fflush(fpdLifty);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeDerivativeOfLiftzToDisk(double& dLz)
{

  if (fpdLiftz) {
    fprintf(fpdLiftz, "%16.13e\n", dLz);
    fflush(fpdLiftz);
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void TsOutput<dim>::writeDerivativeOfForcesToDisk(int it, int actvar, Vec3D & F, Vec3D & dF, Vec3D & M, Vec3D & dM, double &sboom, double &dSboom)
{

  if (fpdForces) {
    fprintf(fpdForces, "%d %d %16.13e %16.13e %16.13e %16.13e %16.13e %16.13e %16.13e %16.13e %16.13e %16.13e %16.13e %16.13e %16.13e %16.13e \n", it, actvar, F[0], F[1], F[2], M[0], M[1], M[2], sboom, dF[0], dF[1], dF[2], dM[0], dM[1], dM[2], dSboom);
    fflush(fpdForces);
  }

}

//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void TsOutput<dim>::writeDerivativeOfMatchPressureToDisk(int it, int actvar, DistSVec<double,1> &dPds, DistSVec<double,3> &X, DistSVec<double,dim> &U, DistVec<double> &A, DistTimeState<dim> *timeState)
{

  if (!Qs_match)     Qs_match     = new DistVec<double>(domain->getNodeDistInfo());
  if (!Qs_match_opt) Qs_match_opt = new DistSVec<double,1>(domain->getNodeDistInfo());

  com->fprintf(stdout, "\nReading optimal pressure distribution from %s\n", fullOptPressureName);
  //domain->readVectorFromFile(fullPressureName,0,0,*Qs_match);
  domain->readVectorFromFile(fullOptPressureName,0,0,*Qs_match_opt);
  postOp->computeScalarQuantity(PostFcn::PRESSURE, X, U, A, *Qs_match, timeState);
  DistSVec<double,1> Qs1(Qs_match->info(), reinterpret_cast<double (*)[1]>(Qs_match->data()));
  //DistSVec<double,1> Qs2(Qs_match_opt->info(), reinterpret_cast<double (*)[1]>(Qs_match_opt->data()));
  if (optPressureDimensional)
    Qs1-=((*Qs_match_opt)*(1/sscale[2]));
  else
    Qs1-=(*Qs_match_opt);
  //Qs1-=*Qs_match_opt;

  double normPdiff, dnormPdiff;
  normPdiff  = 0.5*(Qs1*Qs1);
  dnormPdiff = Qs1*dPds;

  normPdiff*=(sscale[2]*sscale[2]);
  dnormPdiff*=(sscale[2]*sscale[2]);
  //normPdiff  = 0.5*((*Qs_match)*(*Qs_match));
  //dnormPdiff = (*Qs_match)*dPds;

  com->fprintf(stderr,"normPdiff = %20.16e\n",normPdiff);
  com->fprintf(stderr,"dnormPdiff = %20.16e\n",dnormPdiff);

  if (fpdMatchPressure) {
    fprintf(fpdMatchPressure, "%d %d %16.13e %16.13e\n", it, actvar,normPdiff,dnormPdiff);
    fflush(fpdMatchPressure);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeHydroForcesToDisk(bool lastIt, int it, int itSc, int itNl, double t, double cpu,
                                           double* e, DistSVec<double,3> &X, DistSVec<double,dim> &U,
                                           DistVec<int> *fluidId)
{

  int nSurfs = postOp->getNumSurf();

  Vec3D *FiS = new Vec3D[nSurfs];
  Vec3D *MiS = new Vec3D[nSurfs];
  Vec3D *FvS = new Vec3D[nSurfs];
  Vec3D *MvS = new Vec3D[nSurfs];
  Vec3D *FiD = new Vec3D[nSurfs];
  Vec3D *MiD = new Vec3D[nSurfs];
  Vec3D *FvD = new Vec3D[nSurfs];
  Vec3D *MvD = new Vec3D[nSurfs];

  double time = refVal->time * t;

  // if single-phase flow -- fluidId is a null pointer
  // if multi-phase flow  -- fluidId points to a DistVec<int>
  if (hydrostaticforces)
    postOp->computeForceAndMoment(x0, X, U, fluidId, FiS, MiS, FvS, MvS, 1);
  if (hydrodynamicforces)
    postOp->computeForceAndMoment(x0, X, U, fluidId, FiD, MiD, FvD, MvD, 2);

  int iSurf;
  if (fpHydroStaticForces[0]) {

    for (iSurf = 0; iSurf < nSurfs; iSurf++)  {
      Vec3D FS = FiS[iSurf] + FvS[iSurf];
      Vec3D MS = MiS[iSurf] + MvS[iSurf];
      if (refVal->mode == RefVal::NON_DIMENSIONAL) {
        FS *= 2.0 * refVal->length*refVal->length / surface;
        MS *= 2.0 * refVal->length*refVal->length*refVal->length / (surface * length);
      }
      else {
        FS *= refVal->force;
        MS *= refVal->energy;
      }
      double energy = refVal->energy * e[0];

      if (rmmh) {
        double tag = rmmh->getTagValue(t);
        fprintf(fpHydroStaticForces[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e \n",
                it, time, itSc, itNl, FS[0], FS[1], FS[2], MS[0], MS[1], MS[2], energy, tag);
      }
      else
        fprintf(fpHydroStaticForces[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e \n",
                it, time, itSc, itNl, FS[0], FS[1], FS[2], MS[0], MS[1], MS[2], energy);

      fflush(fpHydroStaticForces[iSurf]);
    }
  }

  if (fpHydroDynamicForces[0]) {

    for (iSurf = 0; iSurf < nSurfs; iSurf++)  {
      Vec3D FD = FiD[iSurf] + FvD[iSurf];
      Vec3D MD = MiD[iSurf] + MvD[iSurf];
      if (refVal->mode == RefVal::NON_DIMENSIONAL) {
        FD *= 2.0 * refVal->length*refVal->length / surface;
        MD *= 2.0 * refVal->length*refVal->length*refVal->length / (surface * length);
      }
      else {
        FD *= refVal->force;
        MD *= refVal->energy;
      }
      double energy = refVal->energy * e[0];

      if (rmmh) {
        double tag = rmmh->getTagValue(t);
        fprintf(fpHydroDynamicForces[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e\n",
                it, time, itSc, itNl, FD[0], FD[1], FD[2], MD[0], MD[1], MD[2], energy, tag);
      }
      else
        fprintf(fpHydroDynamicForces[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e %13.16e \n",
                it, time, itSc, itNl, FD[0], FD[1], FD[2], MD[0], MD[1], MD[2], energy);

      fflush(fpHydroDynamicForces[iSurf]);
    }
  }

  delete []FiS;
  delete []MiS;
  delete []FvS;
  delete []MvS;
  delete []FiD;
  delete []MiD;
  delete []FvD;
  delete []MvD;

}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeLiftsToDisk(IoData &iod, bool lastIt, int it, int itSc, int itNl, double t, double cpu,
                                     double* e, DistSVec<double,3> &X, DistSVec<double,dim> &U,
                                     DistVec<int> *fluidId)
{

// This routine outputs both the non-averaged and time-averaged values of the lift and drag
// in ascii format

  int nSurfs = postOp->getNumSurf();

  Vec3D *Fi = new Vec3D[nSurfs];
  Vec3D *Mi = new Vec3D[nSurfs];
  Vec3D *Fv = new Vec3D[nSurfs];
  Vec3D *Mv = new Vec3D[nSurfs];

  double del_t;

  if (TavL == 0) {
  TavL = new Vec3D[nSurfs];
  for(int i = 0; i < nSurfs; ++i) {
      TavL[i] = 0.0;
  }
  }

  int iSurf;
  if(counter == 0) {
    tinit = refVal->time*t;
    tprevl = refVal->time*t;
    del_t = 0.0;
  }

  // if single-phase flow -- fluidId is a null pointer
  // if multi-phase flow  -- fluidId points to a DistVec<int>
  if (lift || tavlift)
    postOp->computeForceAndMoment(x0, X, U, fluidId, Fi, Mi, Fv, Mv);

  double time = refVal->time * t;

  if (fpLift[0] || fpTavLift[0]) {
    for (iSurf = 0; iSurf < nSurfs; iSurf++)  {
      Vec3D L;
      Vec3D F = Fi[iSurf] + Fv[iSurf];
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        F *= 2.0 * refVal->length*refVal->length / surface;
      else
       F *= refVal->force;

      L[0] = F[0]*cos(iod.bc.inlet.alpha)*cos(iod.bc.inlet.beta) +
             F[1]*cos(iod.bc.inlet.alpha)*sin(iod.bc.inlet.beta) +
             F[2]*sin(iod.bc.inlet.alpha);

      L[1] = -F[0]*sin(iod.bc.inlet.beta) + F[1]*cos(iod.bc.inlet.beta);

      L[2] = -F[0]*sin(iod.bc.inlet.alpha)*cos(iod.bc.inlet.beta) -
              F[1]*sin(iod.bc.inlet.alpha)*sin(iod.bc.inlet.beta) +
              F[2]*cos(iod.bc.inlet.alpha);

      if (rmmh) {
        double tag = rmmh->getTagValue(t);
        fprintf(fpLift[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e %13.16e \n",
                it, time, itSc, itNl, L[0], L[1], L[2], tag);
      }
      else
        fprintf(fpLift[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e \n",
                it, time, itSc, itNl, L[0], L[1], L[2]);
      fflush(fpLift[iSurf]);

      if(fpTavLift[iSurf] && counter == 0){
        fprintf(fpTavLift[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e \n",
             it, time, itSc, itNl, L[0], L[1], L[2]);
      }

      del_t = time - tprevl;
      L *= del_t;  TavL[iSurf] += L;

      if(fpTavLift[iSurf] && counter > 0){
        L = TavL[iSurf];
        L /= (time - tinit);
        fprintf(fpTavLift[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e \n",
                it, time, itSc, itNl, L[0], L[1], L[2]);
      }
      fflush(fpTavLift[iSurf]);
    }
  }

  tprevl = time;

  delete [] Fi;
  delete [] Fv;
  delete [] Mi;
  delete [] Mv;

}
//------------------------------------------------------------------------------


template<int dim>
void TsOutput<dim>::writeHydroLiftsToDisk(IoData &iod, bool lastIt, int it, int itSc, int itNl, double t, double cpu,
                                      double* e, DistSVec<double,3> &X, DistSVec<double,dim> &U,
                                      DistVec<int> *fluidId)
{

  int nSurfs = postOp->getNumSurf();

  Vec3D *FiS = new Vec3D[nSurfs];
  Vec3D *MiS = new Vec3D[nSurfs];
  Vec3D *FvS = new Vec3D[nSurfs];
  Vec3D *MvS = new Vec3D[nSurfs];
  Vec3D *FiD = new Vec3D[nSurfs];
  Vec3D *MiD = new Vec3D[nSurfs];
  Vec3D *FvD = new Vec3D[nSurfs];
  Vec3D *MvD = new Vec3D[nSurfs];

  double time = refVal->time * t;

  // if single-phase flow -- fluidId is a null pointer
  // if multi-phase flow  -- fluidId  points to a DistVec<int>
  if (hydrostaticlift)
    postOp->computeForceAndMoment(x0, X, U, fluidId, FiS, MiS, FvS, MvS, 1);
  if (hydrodynamiclift)
    postOp->computeForceAndMoment(x0, X, U, fluidId, FiD, MiD, FvD, MvD, 2);

  int iSurf;
  if (fpHydroStaticLift[0]) {
    for (iSurf = 0; iSurf < nSurfs; iSurf++)  {
      Vec3D LS;
      Vec3D FS = FiS[iSurf] + FvS[iSurf];
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        FS *= 2.0 * refVal->length*refVal->length / surface;
      else
       FS *= refVal->force;

      LS[0] = FS[0]*cos(iod.bc.inlet.alpha)*cos(iod.bc.inlet.beta) +
              FS[1]*cos(iod.bc.inlet.alpha)*sin(iod.bc.inlet.beta) +
              FS[2]*sin(iod.bc.inlet.alpha);

      LS[1] = -FS[0]*sin(iod.bc.inlet.beta) + FS[1]*cos(iod.bc.inlet.beta);

      LS[2] = -FS[0]*sin(iod.bc.inlet.alpha)*cos(iod.bc.inlet.beta) -
               FS[1]*sin(iod.bc.inlet.alpha)*sin(iod.bc.inlet.beta) +
               FS[2]*cos(iod.bc.inlet.alpha);

      if (rmmh) {
        double tag = rmmh->getTagValue(t);
        fprintf(fpHydroStaticLift[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e %13.16e \n",
                it, time, itSc, itNl, LS[0], LS[1], LS[2], tag);
      }
      else
        fprintf(fpHydroStaticLift[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e \n",
                it, time, itSc, itNl, LS[0], LS[1], LS[2]);
      fflush(fpHydroStaticLift[iSurf]);
    }
  }

  if (fpHydroDynamicLift[0]) {
    for (iSurf = 0; iSurf < nSurfs; iSurf++)  {
      Vec3D LD;
      Vec3D FD = FiD[iSurf] + FvD[iSurf];
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        FD *= 2.0 * refVal->length*refVal->length / surface;
      else
        FD *= refVal->force;

      LD[0] = FD[0]*cos(iod.bc.inlet.alpha)*cos(iod.bc.inlet.beta) +
              FD[1]*cos(iod.bc.inlet.alpha)*sin(iod.bc.inlet.beta) +
              FD[2]*sin(iod.bc.inlet.alpha);

      LD[1] = -FD[0]*sin(iod.bc.inlet.beta) + FD[1]*cos(iod.bc.inlet.beta);

      LD[2] = -FD[0]*sin(iod.bc.inlet.alpha)*cos(iod.bc.inlet.beta) -
               FD[1]*sin(iod.bc.inlet.alpha)*sin(iod.bc.inlet.beta) +
               FD[2]*cos(iod.bc.inlet.alpha);

      if (rmmh) {
        double tag = rmmh->getTagValue(t);
        fprintf(fpHydroDynamicLift[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e %13.16e \n",
                it, time, itSc, itNl, LD[0], LD[1], LD[2], tag);
      }
      else
        fprintf(fpHydroDynamicLift[iSurf], "%d %13.16e %d %d %13.16e %13.16e %13.16e \n",
                it, time, itSc, itNl, LD[0], LD[1], LD[2]);
      fflush(fpHydroDynamicLift[iSurf]);
    }
  }

  delete [] FiS;
  delete [] FvS;
  delete [] MiS;
  delete [] MvS;
  delete [] FiD;
  delete [] FvD;
  delete [] MiD;
  delete [] MvD;
}

//------------------------------------------------------------------------------
template<int dim>
  void TsOutput<dim>::writeHeatFluxesToDisk(bool lastIt, int it, int itSc, int itNl, double t, double cpu,
                                      double* e, DistSVec<double,3> &X, DistSVec<double,dim> &U,
                                      DistVec<int> *fluidId)
{
  int nSurfs = postOp->getNumSurfHF();

  double *HF = new double[nSurfs];
  for(int index =0; index < nSurfs; index++){
    HF[index] = 0.0;
  }

  double time = refVal->time * t;

  if (heatfluxes)
    postOp->computeHeatFluxes(X,U,HF);

  int iSurf;
  if (fpHeatFluxes[0]) {
    for (iSurf = 0; iSurf < nSurfs; iSurf++)  {
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
        HF[iSurf] *= 2.0 * refVal->length*refVal->length / surface;
      else
       HF[iSurf] *= refVal->power;
        fprintf(fpHeatFluxes[iSurf], "%d %e %d %d %e \n",
                it, time, itSc, itNl, HF[iSurf]);
      fflush(fpHeatFluxes[iSurf]);
    }
  }
  delete[] HF;
}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeResidualsToDisk(int it, double cpu, double res, double cfl)
{

  if (com->cpuNum() != 0) return;

  if (fpResiduals) {
    fprintf(fpResiduals, "%d %e %e %e\n", it, cpu, res, cfl);
    fflush(fpResiduals);
  }

  if (steady)
    com->printf(0, "It %5d: Res = %e, Cfl = %e, Elapsed Time = %.2e s\n", it, res, cfl, cpu);

}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeMatchStateToDisk(IoData &iod,  int it, double t, double cpu, DistSVec<double,dim> &U, DistVec<double> &A)
{

  // monitoring convergence to a reference solution (error = |Uref - U|/|Uref|)

  if (iod.output.transient.matchstate[0] == 0) return;

  double time = refVal->time * t;

  if (Uref_norm<0) {
    *Uref *= A;
    Uref_norm = Uref->norm();
  }

  DistSVec<double, dim> Udif(U);
  Udif *= A;
  Udif -= *Uref;

  double relDif;
  relDif = Udif.norm()/Uref_norm;

  if (com->cpuNum() != 0) return;
  if (fpMatchState) {
    fprintf(fpMatchState, "%d %e %e %e \n",it, time, cpu, relDif);
    fflush(fpMatchState);
  }

}


//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeMatchPressureToDisk(IoData &iod, bool lastIt, int it, int itSc, int itNl, double t, double cpu,
                                     double* e, DistSVec<double,3> &X, DistVec<double> &A, DistSVec<double,dim> &U,
                                     DistTimeState<dim> * timeState, DistVec<int> *fluidId)
{

  if (iod.output.transient.matchpressure[0] == 0) return;

  double time = refVal->time * t;

  if (!Qs_match)     Qs_match     = new DistVec<double>(domain->getNodeDistInfo());

  postOp->computeScalarQuantity(PostFcn::PRESSURE, X, U, A, *Qs_match, timeState);
  DistSVec<double,1> Qs1(Qs_match->info(), reinterpret_cast<double (*)[1]>(Qs_match->data()));

  if (optPressureDimensional)
    Qs1-=((*Qs_match_opt)*(1/sscale[2]));
  else
    Qs1-=(*Qs_match_opt);

  double normPdiff;
  normPdiff  = 0.5*(Qs1*Qs1);
  normPdiff*=(sscale[2]*sscale[2]);

  if (com->cpuNum() != 0) return;
  if (fpMatchPressure) {
    fprintf(fpMatchPressure, "%d %e %d %d %e \n",it, time, itSc, itNl, normPdiff);
    fflush(fpMatchPressure);
  }

}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeFluxNormToDisk(int it, int itSc, int itNl, double t, double normFlux)
{

 double time = refVal->time * t;

  if (com->cpuNum() != 0) return;

  if (fpFluxNorm) {
    fprintf(fpFluxNorm, "%d %e %d %d %e \n",it, time, itSc, itNl, normFlux);
    fflush(fpFluxNorm);
  }

}
//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeMaterialVolumesToDisk(int it, double t, DistVec<double> &A, DistVec<int> *fluidId)
{
 if(!material_volumes)
    return;

  int myLength = numFluidPhases + 1/*ghost*/;
  double Vol[myLength];
  for(int i=0; i<myLength; i++)
    Vol[i] = 0.0;

  domain->computeMaterialVolumes(Vol,myLength,A,fluidId); //computes Vol

  if (com->cpuNum() !=0 ) return;

  double length3 = length*length*length;
  for(int i=0; i<myLength; i++)
    Vol[i] *= length3; //dimensionalize

  fprintf(fpMatVolumes, "%d %e ", it, (refVal->time)*t);
  for(int i=0; i<numFluidPhases+1; i++)
    fprintf(fpMatVolumes, "%e ", Vol[i]);

  double totVol = 0.0;
  for(int i=0; i<myLength; i++)
    totVol += Vol[i];

  fprintf(fpMatVolumes, "%e\n", totVol);

  fflush(fpMatVolumes);

}



//------------------------------------------------------------------------------
/** Function writeMaterialMassEnergyToDisk
   *  U is the conservative state variables
   *  A is the volume of fluid control volume
   *  fluidId is the fluid Id vector for mutiphase problem, and NULL for single phase problem
   */
template<int dim>
void TsOutput<dim>::writeMaterialMassEnergyToDisk(int it, double t,DistSVec<double,dim> & U,  DistVec<double> &A, DistVec<int> *fluidId)
{
  if(!material_mass_energy)
    return;

  int myLength = numFluidPhases + 1/*ghost*/;
  double Mass[myLength];
  double Energy[myLength];

  for(int i=0; i<myLength; i++) {
    Mass[i] = 0.0;
    Energy[i] = 0.0;
  }

  domain->computeMaterialMassEnergy(Mass,Energy, myLength,U, A,fluidId); //computes Mass

  if (com->cpuNum() !=0 ) return;


  double massScale   = length*length*length*refVal->density;
  double energyScale = length*length*length*refVal->energy;

  for(int i=0; i<myLength; i++) {
    Mass[i] *= massScale; //dimensionalize
    Energy[i] *= energyScale;
  }



  fprintf(fpMaterialMassEnergy, "%d %e ", it, (refVal->time)*t);

  for(int i=0; i<numFluidPhases+1; i++) {
    fprintf(fpMaterialMassEnergy, "%15.10E %15.10E ", Mass[i], Energy[i]);
  }



  double totMass = 0.0;
  double totEnergy = 0.0;
  for(int i=0; i<myLength; i++) {
    totMass += Mass[i];
    totEnergy += Energy[i];
  }
  fprintf(fpMaterialMassEnergy, "%15.10E %15.10E\n", totMass,totEnergy);

  fflush(fpMaterialMassEnergy);

}


//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeEmbeddedSurfaceToDisk(bool lastIt, int it, double t, Vec<Vec3D>& solidX, Vec<Vec3D>& solidX0)
{
  if(!embeddedsurface)
    return;

  if(toWrite(it,lastIt,t)) {
    if(com->cpuNum() != 0) return;
    if(it == 0)  fprintf(fpEmbeddedSurface, "%d\n", solidX.size());
    fprintf(fpEmbeddedSurface, "%e\n", t*tscale);
    for(int i=0; i<solidX.size(); i++)
      fprintf(fpEmbeddedSurface, "%e %e %e\n", xscale*(solidX[i][0]-solidX0[i][0]), xscale*(solidX[i][1]-solidX0[i][1]),
              xscale*(solidX[i][2]-solidX0[i][2]));
    fflush(fpEmbeddedSurface);
    fprintf(stdout, "Wrote solution %d to \'%s\'\n", getStep(it, lastIt, t), embeddedsurface);
  }
}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeCPUTimingToDisk(bool lastIt, int it, double t, Timer *timer)
{
  if(!cputiming)
    return;

  if(toWrite(it,lastIt,t)) {
    com->fprintf(fpCpuTiming, "It %d, Time = %e.", it, t*tscale);
    timer->setRunTime();
    timer->print(NULL, fpCpuTiming);
    if(com->cpuNum()==0)
      fflush(fpCpuTiming);
    com->fprintf(stdout, "Wrote CPUTiming %d to \'%s\'\n", getStep(it, lastIt, t), cputiming);
  }
}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeConservationErrors(IoData &iod, int it, double t,
                    int numPhases, double **expected, double **computed)
{

  if (com->cpuNum() != 0) return;

  if (!steady)
    if (fpConservationErr) {
      fprintf(fpConservationErr, "%d %e ", it, t);
      for(int iPhase=0; iPhase<numPhases; iPhase++)
        for(int i=0; i<dim; i++)
          fprintf(fpConservationErr, "%e ", expected[iPhase][i]);
      for(int iPhase=0; iPhase<numPhases; iPhase++)
        for(int i=0; i<dim; i++)
          fprintf(fpConservationErr, "%e ", computed[iPhase][i]);
      fprintf(fpConservationErr, "\n");
      fflush(fpConservationErr);

    }

}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeDisplacementVectorToDisk(int step, double tag,
                      DistSVec<double,3> &X, DistSVec<double,dim> &U){

  if(vectors[PostFcn:: DISPLACEMENT]){
    if (!Qv) Qv = new DistSVec<double,3>(domain->getNodeDistInfo());
    postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(PostFcn::DISPLACEMENT), X, U, *Qv);
    domain->writeVectorToFile(vectors[PostFcn::DISPLACEMENT], 1, 1.0, *Qv, &(vscale[PostFcn::DISPLACEMENT]));
  }


}

//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writePositionSensitivityVectorToDisk(int step, double tag,
                      DistSVec<double,3> &dXsds){

  domain->writeVectorToFile(dVectors[PostFcn::DERIVATIVE_DISPLACEMENT], step, tag, dXsds, &(refVal->tlength));

}




//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeBinaryVectorsToDisk(bool lastIt,
					     int it,
					     double t,
					     DistSVec<double,3> &X,
                                             DistVec<double> &A,
					     DistSVec<double,dim> &U,
					     DistTimeState<dim> *timeState)
{

  if (toWrite(it,lastIt,t)) {
    int step = getStep(it,lastIt,t);
    double tag;
    if (rmmh)
      tag = rmmh->getTagValue(t);
    else
      tag = t * refVal->time;


    int i;
    for (i=0; i<PostFcn::SSIZE; ++i) {
      if (scalars[i]) {
        if (!Qs) Qs = new DistVec<double>(domain->getNodeDistInfo());
        postOp->computeScalarQuantity(static_cast<PostFcn::ScalarType>(i), X, U, A, *Qs, timeState);
        DistSVec<double,1> Qs1(Qs->info(), reinterpret_cast<double (*)[1]>(Qs->data()));
        domain->writeVectorToFile(scalars[i], step, tag, Qs1, &(sscale[i]));
      }
    }
    for (i=0; i<PostFcn::VSIZE; ++i) {
      if (vectors[i]) {
        if (!Qv) Qv = new DistSVec<double,3>(domain->getNodeDistInfo());

        if (static_cast<PostFcn::VectorType>(i) == PostFcn::FLIGHTDISPLACEMENT)  {

          if (rmmh) {
            DistSVec<double,3> &Xr = rmmh->getFlightPositionVector(t, X);
            postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), Xr, U, *Qv);
          }
          else
            com->fprintf(stderr, "WARNING: Flight Displacement Output not available\n");
        }
        else if (static_cast<PostFcn::VectorType>(i) == PostFcn::LOCALFLIGHTDISPLACEMENT)  {
          if (rmmh) {
            DistSVec<double,3> &Xr = rmmh->getRelativePositionVector(t, X);
            postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), Xr, U, *Qv);
          }
          else
            com->fprintf(stderr, "WARNING: Local Flight Displacement Output not available\n");

        }
        else
          postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), X, U, *Qv);
        domain->writeVectorToFile(vectors[i], step, tag, *Qv, &(vscale[i]));
      }
    }
  }

}

//------------------------------------------------------------------------------

template<int dim>
template<int dimLS>
void TsOutput<dim>::writeBinaryVectorsToDisk(bool lastIt, int it, double t, DistSVec<double,3> &X,
                                             DistVec<double> &A, DistSVec<double,dim> &U,
                                             DistTimeState<dim> *timeState,
                                             DistVec<int> &fluidId,DistSVec<double,dimLS>* Phi)
{

  if (toWrite(it,lastIt,t)) {
    int step = getStep(it,lastIt,t);
    double tag;
    if (rmmh)
      tag = rmmh->getTagValue(t);
    else
      tag = t * refVal->time;

    //if (stateVector) // MOR not supported for embedded framework
    //  domain->writeVectorToFile(stateVector, step, tag, U);

    int i;
    for (i=0; i<PostFcn::SSIZE; ++i) {
      if (scalars[i]) {
        if (!Qs) Qs = new DistVec<double>(domain->getNodeDistInfo());
        postOp->computeScalarQuantity(static_cast<PostFcn::ScalarType>(i), X, U, A, *Qs, timeState,fluidId,Phi);
        DistSVec<double,1> Qs1(Qs->info(), reinterpret_cast<double (*)[1]>(Qs->data()));
        domain->writeVectorToFile(scalars[i], step, tag, Qs1, &(sscale[i]));
      }
    }
    for (i=0; i<PostFcn::VSIZE; ++i) {
      if (vectors[i]) {
        if (!Qv) Qv = new DistSVec<double,3>(domain->getNodeDistInfo());

        if (static_cast<PostFcn::VectorType>(i) == PostFcn::FLIGHTDISPLACEMENT)  {

          if (rmmh) {
            DistSVec<double,3> &Xr = rmmh->getFlightPositionVector(t, X);
            postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), Xr, U, *Qv);
          }
          else
            com->fprintf(stderr, "WARNING: Flight Displacement Output not available\n");
        }
        else if (static_cast<PostFcn::VectorType>(i) == PostFcn::LOCALFLIGHTDISPLACEMENT)  {
          if (rmmh) {
            DistSVec<double,3> &Xr = rmmh->getRelativePositionVector(t, X);
            postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), Xr, U, *Qv);
          }
          else
            com->fprintf(stderr, "WARNING: Local Flight Displacement Output not available\n");

        }
        else
          postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), X, U, *Qv, fluidId);
        domain->writeVectorToFile(vectors[i], step, tag, *Qv, &(vscale[i]));
      }
    }
  }

}

static void copyFile(const char* fname) {

  FILE* f = fopen(fname,"rb");
  fseek (f , 0 , SEEK_END);
  size_t lSize = ftell (f);
  rewind (f);
  char* buffer = new char[lSize];
  size_t err = fread (buffer,1,lSize,f);
  fclose(f);

  char nn[256];
  sprintf(nn,"%s.back",fname);
  f = fopen(nn,"wb");
  fwrite(buffer,1,lSize,f);
  fclose(f);

  delete [] buffer;
}

template<int dim>
void TsOutput<dim>::cleanProbesFile() {

  char nn[256];
  int iter,i,n;
  double time,res;
  if (it0 == 0)
    return;

  for (i=0; i<PostFcn::SSIZE; ++i) {
    if (nodal_scalars[i]) {
      if (com->cpuNum() == 0) {
        copyFile(nodal_scalars[i]);
        sprintf(nn,"%s.back",nodal_scalars[i]);
        FILE* scalar_file = fopen(nodal_scalars[i],"w");
        FILE* scalar_file_old = fopen(nn,"r");
        while (!feof(scalar_file_old)) {
          n = fscanf(scalar_file_old,"%d",&iter);
          n = fscanf(scalar_file_old,"%lf",&time);
          if (iter > it0)
            break;
          fprintf(scalar_file,"%d %e ",iter,time);
          for (int k =0 ; k < nodal_output.numNodes; ++k) {
            n = fscanf(scalar_file_old,"%lf",&res);
            fprintf(scalar_file,"%e ",res);
          }
          fprintf(scalar_file,"\n");

        }
        fclose(scalar_file);
        fclose(scalar_file_old);
      }
    }
  }

  for (i=0; i<PostFcn::VSIZE; ++i) {
    if (nodal_vectors[i]) {

      if (com->cpuNum() == 0) {
        copyFile(nodal_vectors[i]);
        sprintf(nn,"%s.back",nodal_vectors[i]);
        FILE* vector_file = fopen(nodal_vectors[i],"w");
        FILE* vector_file_old = fopen(nn,"r");
        while (!feof(vector_file_old)) {
          n = fscanf(vector_file_old,"%d",&iter);
          n = fscanf(vector_file_old,"%lf",&time);
          if (iter > it0)
            break;
          fprintf(vector_file,"%d %e ",iter,time);
          for (int k =0 ; k < nodal_output.numNodes; ++k) {
            for (int l = 0; l < 3; ++l) {
              n = fscanf(vector_file_old,"%lf",&res);
              fprintf(vector_file,"%e ",res);
            }
          }
          fprintf(vector_file,"\n");
        }
        fclose(vector_file);
        fclose(vector_file_old);
      }
    }
  }
}

template<int dim>
template<int dimLS>
void TsOutput<dim>::writeProbesToDisk(bool lastIt, int it, double t, DistSVec<double,3> &X,
                                      DistVec<double> &A, DistSVec<double,dim> &U,
                                      DistTimeState<dim> *timeState, DistVec<int> &fluidId,
                                      DistSVec<double,dimLS>* Phi, DistLevelSetStructure *distLSS,
                                      DistVec<GhostPoint<dim>*> *ghostPoints)
{
  //if (toWrite(it,lastIt)) {
  if (nodal_output.numNodes == 0)
    return;
    double tag;
    if (rmmh)
      tag = rmmh->getTagValue(t);
    else
      tag = t * refVal->time;


    int i;
    const char* mode = nodal_output.step ? "a" : "w";
    if (it0 > 0)
      mode = "a";

    for (i=0; i<PostFcn::SSIZE; ++i) {
      if (nodal_scalars[i]) {

        postOp->computeScalarQuantity(static_cast<PostFcn::ScalarType>(i), X, U, A,
                                      timeState,fluidId,
                                      nodal_output.subId, nodal_output.locNodeId,
                                      nodal_output.last,nodal_output.numNodes,nodal_output.results,
                                      nodal_output.locations,
                                      Phi, distLSS, ghostPoints);
        if (com->cpuNum() == 0) {
          FILE* scalar_file = fopen(nodal_scalars[i],mode);
          if (scalar_file != 0) {
            fprintf(scalar_file,"%d %e ",nodal_output.step+it0, tag);
            for (int k =0 ; k < nodal_output.numNodes; ++k)
              fprintf(scalar_file,"%e ",nodal_output.results[k]*sscale[i]);
            fprintf(scalar_file,"\n");
            fclose(scalar_file);
          } else {

            this->com->fprintf(stderr,"Warning: Cannot open probe file %s",nodal_scalars[i]);
          }
        }
      }
    }

    for (i=0; i<PostFcn::VSIZE; ++i) {
      if (nodal_vectors[i]) {

        postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), X, U,
                                      nodal_output.subId, nodal_output.locNodeId,
                                      nodal_output.last,nodal_output.numNodes,nodal_output.results,
                                      nodal_output.locations,
                                      fluidId, distLSS, ghostPoints);

        if (com->cpuNum() == 0) {
          FILE* vector_file = fopen(nodal_vectors[i],mode);
          if (vector_file != 0) {
            fprintf(vector_file,"%d %e ",nodal_output.step+it0, tag);
            for (int k =0 ; k < nodal_output.numNodes; ++k)
            fprintf(vector_file,"%e %e %e ",
            nodal_output.results[k*3]*vscale[i],
            nodal_output.results[k*3+1]*vscale[i],
            nodal_output.results[k*3+2]*vscale[i]);
            fprintf(vector_file,"\n");
            fclose(vector_file);
          } else {

            this->com->fprintf(stderr,"Warning: Cannot open probe file %s",nodal_vectors[i]);
          }
        }
      }
    }
    // }
    ++nodal_output.step;
}

template<int dim>
template<int dimLS>
void TsOutput<dim>::writeLinePlotsToDisk(bool lastIt, int it, double t, DistSVec<double,3> &X,
					DistVec<double> &A, DistSVec<double,dim> &U,
					DistTimeState<dim> *timeState, DistVec<int> &fluidId,
					DistSVec<double,dimLS>* Phi, DistLevelSetStructure *distLSS,
					DistVec<GhostPoint<dim>*> *ghostPoints)
{
  //if (toWrite(it,lastIt,t)) {
  if (line_outputs.size() == 0)
    return;

    double tag;
    if (rmmh)
      tag = rmmh->getTagValue(t);
    else
      tag = t * refVal->time;

    // if (solutions)
    //  domain->writeVectorToFile(solutions, step, tag, U);

    int i;
    const char* mode = "w";//nodal_output.step ? "a" : "w";
    //if (it0 > 0)
    //  mode = "a";


    for (int j = 0; j < line_outputs.size(); ++j) {

      FILE* scalar_files[PostFcn::SSIZE],*vector_files[PostFcn::VSIZE];
      line_output* lout = line_outputs[j];

      memset(scalar_files,0, sizeof(scalar_files));
      memset(vector_files,0, sizeof(vector_files));


      for (i=0; i<PostFcn::SSIZE; ++i) {

	if (lout->scalars[i]) {
	  if (com->cpuNum() == 0) {
	    scalar_files[i] = fopen(lout->scalars[i],mode);
	    if (scalar_files[i] != 0) {
	    } else {

	      this->com->fprintf(stderr,"Warning: Cannot open line output file %s",lout->scalars[i]);
	    }
	  }
	}
      }

      for (i=0; i<PostFcn::VSIZE; ++i) {

	if (lout->vectors[i]) {
	  if (com->cpuNum() == 0) {
	    vector_files[i] = fopen(lout->vectors[i],mode);
	    if (vector_files[i] != 0) {
	    } else {

	      this->com->fprintf(stderr,"Warning: Cannot open line output file %s",lout->vectors[i]);
	    }
	  }
	}
      }

      int last = 0;
      int subId = -1, locNodeId = -1;
      double result[3];
      std::vector<Vec3D> location(1);
      for (int k = 0; k < lout->numPoints; ++k) {
	location[0] = Vec3D(((lout->x1-lout->x0)/(lout->numPoints-1)*k+lout->x0)/xscale,
			    ((lout->y1-lout->y0)/(lout->numPoints-1)*k+lout->y0)/xscale,
			    ((lout->z1-lout->z0)/(lout->numPoints-1)*k+lout->z0)/xscale);

	for (i=0; i<PostFcn::SSIZE; ++i) {

	  if (lout->scalars[i]) {

	    postOp->computeScalarQuantity(static_cast<PostFcn::ScalarType>(i), X, U, A,
					  timeState,fluidId,
					  &subId, &locNodeId,
					  &last,1,result,
					  location,
					  Phi, distLSS, ghostPoints);

	    if (com->cpuNum() == 0) {
	      fprintf(scalar_files[i],"%e %e %e %e %e\n",(double)k/(lout->numPoints-1),
		      location[0][0]*xscale,location[0][1]*xscale,
		      location[0][2]*xscale,result[0]*sscale[i]);
	    }

	  }
	}

	for (i=0; i<PostFcn::VSIZE; ++i) {
	  if (lout->vectors[i]) {

	    postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), X, U,
					  &subId, &locNodeId,
					  &last,1,result,
					  location,
					  fluidId, distLSS, ghostPoints);

	    if (com->cpuNum() == 0) {
	      fprintf(vector_files[i],"%e %e %e %e %e %e %e\n",(double)k/(lout->numPoints-1),
		      location[0][0]*xscale,location[0][1]*xscale,
		      location[0][2]*xscale,result[0]*vscale[i], result[1]*vscale[i],result[2]*vscale[i]);
	    }
	  }
	}
      }

      for (i=0; i<PostFcn::SSIZE; ++i) {

	if (scalar_files[i])
	  fclose(scalar_files[i]);
      }

      for (i=0; i<PostFcn::VSIZE; ++i) {

	if (vector_files[i])
	  fclose(vector_files[i]);
      }
    }
}

//----------------------------------------------------------------------------------------
// d2d
template<int dim>
void TsOutput<dim>::writeBinaryVectorsToDisk(bool lastIt, int it, double t, DistSVec<double,3> &X,
                                             DistVec<double> &A, DistSVec<double,dim> &U,
                                             DistTimeState<dim> *timeState,
                                             DistVec<int> &fluidId, DistSVec<double,dim> *Wextij,
															DistLevelSetStructure *distLSS,
					     DistVec<GhostPoint<dim>*> *ghostPoints)
{

  if (toWrite(it,lastIt,t)) {

    int step = getStep(it,lastIt,t);
    double tag;
    if (rmmh)
      tag = rmmh->getTagValue(t);
    else
      tag = t * refVal->time;

    if (dSolutions)
      domain->writeVectorToFile(dSolutions, step, tag, U);

    int i;
    for (i=0; i<PostFcn::SSIZE; ++i) {
      if (scalars[i]) {

        if (!Qs) Qs = new DistVec<double>(domain->getNodeDistInfo());

        postOp->computeScalarQuantity(static_cast<PostFcn::ScalarType>(i), X, U, A, *Qs, timeState,fluidId, (DistSVec<double,1>*)0);
        DistSVec<double,1> Qs1(Qs->info(), reinterpret_cast<double (*)[1]>(Qs->data()));

        domain->writeVectorToFile(scalars[i], step, tag, Qs1, &(sscale[i]));

      }

    }

    /////////////////////////////////////////////
    if(embeddedsurfaceCp || embeddedsurfaceCf) {

      int ns = distLSS->getNumStructNodes();

      double** EmbQs;
      EmbQs = new double* [ns];
      for(int i=0; i<ns; ++i) {
	EmbQs[i] = new double[3];
	EmbQs[i][0] = EmbQs[i][1] = EmbQs[i][2] = 0.0;
      }

      postOp->computeEMBScalarQuantity(X, U, A, EmbQs, timeState, fluidId, Wextij,
													(DistSVec<double,1>*)0, distLSS, ghostPoints, externalSI);


      double * cnt = new double[ns];
      for(int i=0; i<ns; ++i) cnt[i] = EmbQs[i][0] ? 1.0 : 0.0;
      com->globalSum(ns, cnt);

      double * Cp_ = new double[ns];
      double * Cf_ = new double[ns];

      if(embeddedsurfaceCp) {
	for(int i=0; i<ns; ++i) Cp_[i] = EmbQs[i][1];
	com->globalSum(ns, Cp_);
      }

      if(embeddedsurfaceCf) {
	for(int i=0; i<ns; ++i) Cf_[i] = EmbQs[i][2];
	com->globalSum(ns, Cf_);
      }

      if(embeddedsurfaceCp) {
	if(toWrite(it,lastIt,t)) {
	  if(com->cpuNum() == 0) {

	    if(it == 0) fprintf(fpEmbeddedSurfaceCp, "%i \n", ns);
	    fprintf(fpEmbeddedSurfaceCp, " %f \n", t*tscale);

	    for(int i=0; i<ns; i++) {
	      double val = cnt[i] ? Cp_[i] /= cnt[i] : 0.0;
	      fprintf(fpEmbeddedSurfaceCp, "%e \n", val);
	    }

	    fflush(fpEmbeddedSurfaceCp);
	    fprintf(stdout, "Wrote solution %d to \'%s\'\n", getStep(it, lastIt, t), embeddedsurfaceCp);

	  }
	}
      }

      // ~~~~

      if(embeddedsurfaceCf) {
	if(toWrite(it,lastIt,t)) {
	  if(com->cpuNum() == 0) {

	    if(it == 0) fprintf(fpEmbeddedSurfaceCf, "%i \n", ns);
	    fprintf(fpEmbeddedSurfaceCf, " %f \n", t*tscale);

	    for(int i=0; i<ns; i++) {
	      double val = cnt[i] ? Cf_[i] /= cnt[i] : 0.0;
	      fprintf(fpEmbeddedSurfaceCf, "%e \n", val);
	    }

	    fflush(fpEmbeddedSurfaceCf);
	    fprintf(stdout, "Wrote solution %d to \'%s\'\n", getStep(it, lastIt, t), embeddedsurfaceCf);

	  }
	}
      }

      for(int i=0; i<ns; ++i) delete [] EmbQs[i];
      delete [] EmbQs;
      delete [] cnt;
      delete [] Cp_;
      delete [] Cf_;
    }
    ///////////////////////////////////////////////

    for (i=0; i<PostFcn::VSIZE; ++i) {
      if (vectors[i]) {
        if (!Qv) Qv = new DistSVec<double,3>(domain->getNodeDistInfo());

        if (static_cast<PostFcn::VectorType>(i) == PostFcn::FLIGHTDISPLACEMENT)  {

          if (rmmh) {
            DistSVec<double,3> &Xr = rmmh->getFlightPositionVector(t, X);
            postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), Xr, U, *Qv);
          }
          else
            com->fprintf(stderr, "WARNING: Flight Displacement Output not available\n");
        }
        else if (static_cast<PostFcn::VectorType>(i) == PostFcn::LOCALFLIGHTDISPLACEMENT)  {
          if (rmmh) {
            DistSVec<double,3> &Xr = rmmh->getRelativePositionVector(t, X);
            postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), Xr, U, *Qv);
          }
          else
            com->fprintf(stderr, "WARNING: Local Flight Displacement Output not available\n");

        }
        else
		  {
			  //postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), X, U, *Qv, fluidId);
			  postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), X, U, *Qv, distLSS, fluidId);
		  }
        domain->writeVectorToFile(vectors[i], step, tag, *Qv, &(vscale[i]));
      }
    }
  }


}

template<int dim>
void TsOutput<dim>::writeProbesToDisk(bool lastIt, int it, double t, DistSVec<double,3> &X,
                                      DistVec<double> &A, DistSVec<double,dim> &U,
                                      DistTimeState<dim> *timeState,
                                      DistVec<int> &fluidId, DistLevelSetStructure *distLSS,
                                      DistVec<GhostPoint<dim>*> *ghostPoints)
{
  writeProbesToDisk(lastIt,it,t,X,A,U,timeState,fluidId, (DistSVec<double,1>*)0, distLSS, ghostPoints);
}

//----------------------------------------------------------------------------------------



//TODO BUGHUNT


// Included (MB)
template<int dim>
void TsOutput<dim>::writeAnyVectorToDisk(
		              const char* filename,
		              int it,
					  int tag,
					  DistSVec<double,dim> &vec)
{
  int    step = it-1;
  domain->writeVectorToFile(filename, step, tag, vec);
}


// Included (MB)
template<int dim>
void TsOutput<dim>::writeBinaryDerivativeOfVectorsToDisk(
		                     int it,
							 int actvar,
							 double dS[3],
							 DistSVec<double,3> &X,
							 DistSVec<double,3> &dX,
							 DistSVec<double,dim> &U,
							 DistSVec<double,dim> &dU,
							 DistTimeState<dim> *timeState,
							 DistVec<double>* A)
{
  int    step = it-1;
  double tag  = (double)actvar;

  if (dSolutions)
    domain->writeVectorToFile(dSolutions, step, tag, dU);

  int i;
  for (i=0; i<PostFcn::DSSIZE; ++i) {
    if (dScalars[i]) {
      if (!Qs) Qs = new DistVec<double>(domain->getNodeDistInfo());
      postOp->computeDerivativeOfScalarQuantity(static_cast<PostFcn::ScalarDerivativeType>(i), dS, X, dX, U, dU, *Qs, timeState);
      DistSVec<double,1> Qs1(Qs->info(), reinterpret_cast<double (*)[1]>(Qs->data()));
      domain->writeVectorToFile(dScalars[i], step, tag, Qs1, &(dSscale[i]));

      //Match Properties
      if (static_cast<PostFcn::ScalarDerivativeType>(i) == PostFcn::DERIVATIVE_PRESSURE  &&
    fullOptPressureName!=NULL)
        writeDerivativeOfMatchPressureToDisk(it,actvar,Qs1,X,U,*A,timeState);
    }
  }
  for (i=0; i<PostFcn::DVSIZE; ++i) {
    if (dVectors[i]) {
      if (!Qv) Qv = new DistSVec<double,3>(domain->getNodeDistInfo());
      if (static_cast<PostFcn::VectorType>(i) != PostFcn::FLIGHTDISPLACEMENT && static_cast<PostFcn::VectorType>(i) != PostFcn::LOCALFLIGHTDISPLACEMENT)
        postOp->computeDerivativeOfVectorQuantity(static_cast<PostFcn::VectorDerivativeType>(i), X, dX, U, dU, *Qv);
      domain->writeVectorToFile(dVectors[i], step, tag, *Qv, &(dVscale[i]));
    }
  }

}

//------------------------------------------------------------------------------
/*
template<int dim>
void TsOutput<dim>::writeBinaryVectorsToDisk(bool lastIt, int it, double t,
                                             DistSVec<double,3> &X,
                                             DistVec<double> &A, DistSVec<double,dim> &U,
                                             DistSVec<double,1> &Phi, DistVec<int> &fluidId)
{

  if (toWrite(it,lastIt,t)) {
    int step = getStep(it,lastIt,t);
    double tag;
    if (rmmh)
      tag = rmmh->getTagValue(t);
    else
      tag = t * refVal->time;

    if (solutions)
      domain->writeVectorToFile(solutions, step, tag, U);

    int i;
    for (i=0; i<PostFcn::SSIZE; ++i) {
      if (scalars[i]) {
        if (!Qs) Qs = new DistVec<double>(domain->getNodeDistInfo());
        postOp->computeScalarQuantity(static_cast<PostFcn::ScalarType>(i), X, U, *Qs, Phi, fluidId);
        DistSVec<double,1> Qs1(Qs->info(), reinterpret_cast<double (*)[1]>(Qs->data()));
        domain->writeVectorToFile(scalars[i], step, tag, Qs1, &(sscale[i]));
      }
    }
    for (i=0; i<PostFcn::VSIZE; ++i) {
      if (vectors[i]) {
        if (!Qv) Qv = new DistSVec<double,3>(domain->getNodeDistInfo());
        if (rmmh) {
          DistSVec<double,3> &Xr = rmmh->getRelativePositionVector(t, X);
          postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), Xr, U, *Qv);
        }
        else
          postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), X, U, *Qv, fluidId);
        domain->writeVectorToFile(vectors[i], step, tag, *Qv, &(vscale[i]));
      }
    }
  }

}
*/
//------------------------------------------------------------------------------

template<int dim>
void TsOutput<dim>::writeAvgVectorsToDisk(bool lastIt, int it, double t, DistSVec<double,3> &X,
                                             DistVec<double> &A, DistSVec<double,dim> &U, DistTimeState<dim> *timeState)
{

// This routine outputs the time averaged values of the scalar and vector output files
// in binary format
  int i;
  static double tprev,tinit;
  double time = refVal->time*t;
  double del_t;

  int step = getStep(it,lastIt,t);

  double tag;
  if (rmmh)
    tag = rmmh->getTagValue(t);
  else
    tag = t * refVal->time;


  if (counter == 0){
    for (i=0; i<PostFcn::AVSSIZE; ++i) {
      if(!AvQs[i]) AvQs[i] = new DistVec<double>(domain->getNodeDistInfo());
      *AvQs[i] = 0.0;
    }
    for (i=0; i<PostFcn::AVSSIZE; ++i) {
      if(avscalars[i]) {
        if (!Qs) Qs = new DistVec<double>(domain->getNodeDistInfo());
        postOp->computeScalarQuantity(static_cast<PostFcn::ScalarType>(i), X, U, A,*Qs, timeState);
        DistSVec<double,1> Qs1(Qs->info(), reinterpret_cast<double (*)[1]>(Qs->data()));
        domain->writeVectorToFile(avscalars[i], step, tag, Qs1, &(avsscale[i]));
      }
    }

    for (i=0; i<PostFcn::AVVSIZE; ++i) {
      if(!AvQv[i]) AvQv[i] = new DistSVec<double,3>(domain->getNodeDistInfo());
      *AvQv[i] = 0.0;
    }

    for (i=0; i<PostFcn::AVVSIZE; ++i) {
      if(avvectors[i]) {
        if (!Qv) Qv = new DistSVec<double,3>(domain->getNodeDistInfo());
        if (rmmh) {
          DistSVec<double,3> &Xr = rmmh->getRelativePositionVector(t, X);
          postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), Xr, U, *Qv);
        }
        else
        postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), X, U, *Qv);
        domain->writeVectorToFile(avvectors[i], step, tag, *Qv, &(avvscale[i]));
      }
    }

    tinit = time;
    tprev = time;
  }


  if (counter > 0){
    del_t = time - tprev;
    for (i=0; i<PostFcn::AVSSIZE; ++i) {
      if (avscalars[i]) {
        if (!Qs) Qs = new DistVec<double>(domain->getNodeDistInfo());
        postOp->computeScalarQuantity(static_cast<PostFcn::ScalarType>(i), X, U, A, *Qs, timeState);
        *Qs *= del_t;
        *AvQs[i] += *Qs;
      }
    }
    for (i=0; i<PostFcn::AVVSIZE; ++i) {
      if (avvectors[i]) {
        if (!Qv) Qv = new DistSVec<double,3>(domain->getNodeDistInfo());
        if (rmmh) {
          DistSVec<double,3> &Xr = rmmh->getRelativePositionVector(t, X);
          postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), Xr, U, *Qv);
        }
        else
          postOp->computeVectorQuantity(static_cast<PostFcn::VectorType>(i), X, U, *Qv);
        *Qv *= del_t;
        *AvQv[i] += *Qv;
      }
    }
    tprev = time;
  }

  if (toWrite(it,lastIt,t) && counter>0) {
    for (i=0; i<PostFcn::AVSSIZE; ++i) {
      if(avscalars[i]) {
        if (!Qs) Qs = new DistVec<double>(domain->getNodeDistInfo());
        *Qs = *AvQs[i];
        *Qs *= (1.0/(time-tinit));
        DistSVec<double,1> Qs1(Qs->info(), reinterpret_cast<double (*)[1]>(Qs->data()));
        domain->writeVectorToFile(avscalars[i], step, tag, Qs1, &(avsscale[i]));
      }
    }

    for (i=0; i<PostFcn::AVVSIZE; ++i) {
      if(avvectors[i]) {
        if (!Qv) Qv = new DistSVec<double,3>(domain->getNodeDistInfo());
        *Qv = *AvQv[i];
        *Qv *= (1.0/(time-tinit));
        domain->writeVectorToFile(avvectors[i], step, tag, *Qv, &(avvscale[i]));
      }
    }
  }

  if (lastIt) {
    // Before deletion, check that pointers have been allocated
    for (i=0; i<PostFcn::AVSSIZE; ++i)
      if ((avscalars[i]) && (AvQs[i])) {
        delete AvQs[i];
        AvQs[i] = 0;
      }
    for (i=0; i<PostFcn::AVVSIZE; ++i)
      if ((avvectors[i]) && (AvQv[i])) {
        delete AvQv[i];
        AvQv[i] = 0;
      }
  }

  counter += 1; // increment the counter for keeping track of the averaging

}


//------------------------------------------------------------------------------

// Included (MB)
template<int dim>
void TsOutput<dim>::rstVar(IoData &iod) {

  int sp = strlen(iod.output.transient.prefix) + 1;

  if (iod.output.transient.density[0] != 0) {
    sscale[PostFcn::DENSITY] = iod.ref.rv.density;
    scalars[PostFcn::DENSITY] = new char[sp + strlen(iod.output.transient.density)];
    sprintf(scalars[PostFcn::DENSITY], "%s%s",
            iod.output.transient.prefix, iod.output.transient.density);
  }
  if (iod.output.transient.tavdensity[0] != 0) {
    avsscale[PostFcn::DENSITYAVG] = iod.ref.rv.density;
    avscalars[PostFcn::DENSITYAVG] = new char[sp + strlen(iod.output.transient.tavdensity)];
    sprintf(avscalars[PostFcn::DENSITYAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavdensity);
  }
  if (iod.output.transient.speed[0] != 0) {
    sscale[PostFcn::SPEED] = iod.ref.rv.velocity;
    scalars[PostFcn::SPEED] = new char[sp + strlen(iod.output.transient.speed)];
    sprintf(scalars[PostFcn::SPEED], "%s%s",
            iod.output.transient.prefix, iod.output.transient.speed);
  }
  if (iod.output.transient.pressure[0] != 0) {
    sscale[PostFcn::PRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::PRESSURE] = new char[sp + strlen(iod.output.transient.pressure)];
    sprintf(scalars[PostFcn::PRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.pressure);
  }
  if (iod.output.transient.diffpressure[0] != 0) {
    sscale[PostFcn::DIFFPRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::DIFFPRESSURE] = new char[sp + strlen(iod.output.transient.diffpressure)];
    sprintf(scalars[PostFcn::DIFFPRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.diffpressure);
  }
  if (iod.output.transient.tavpressure[0] != 0) {
    avsscale[PostFcn::PRESSUREAVG] = iod.ref.rv.pressure;
    avscalars[PostFcn::PRESSUREAVG] = new char[sp + strlen(iod.output.transient.tavpressure)];
    sprintf(avscalars[PostFcn::PRESSUREAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavpressure);
  }
  if (iod.output.transient.hydrostaticpressure[0] != 0) {
    sscale[PostFcn::HYDROSTATICPRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::HYDROSTATICPRESSURE] = new char[sp + strlen(iod.output.transient.hydrostaticpressure)];
    sprintf(scalars[PostFcn::HYDROSTATICPRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.hydrostaticpressure);
  }
  if (iod.output.transient.hydrodynamicpressure[0] != 0) {
    sscale[PostFcn::HYDRODYNAMICPRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::HYDRODYNAMICPRESSURE] = new char[sp + strlen(iod.output.transient.hydrodynamicpressure)];
    sprintf(scalars[PostFcn::HYDRODYNAMICPRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.hydrodynamicpressure);
  }
  if (iod.output.transient.pressurecoefficient[0] != 0) {
    sscale[PostFcn::PRESSURECOEFFICIENT] = 1.0;
    scalars[PostFcn::PRESSURECOEFFICIENT] = new char[sp + strlen(iod.output.transient.pressurecoefficient)];
    sprintf(scalars[PostFcn::PRESSURECOEFFICIENT], "%s%s",
            iod.output.transient.prefix, iod.output.transient.pressurecoefficient);
  }
  if (iod.output.transient.temperature[0] != 0) {
    sscale[PostFcn::TEMPERATURE] = iod.ref.rv.temperature;
    scalars[PostFcn::TEMPERATURE] = new char[sp + strlen(iod.output.transient.temperature)];
    sprintf(scalars[PostFcn::TEMPERATURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.temperature);
  }
  if (iod.output.transient.tavtemperature[0] != 0) {
    avsscale[PostFcn::TEMPERATUREAVG] = iod.ref.rv.temperature;
    avscalars[PostFcn::TEMPERATUREAVG] = new char[sp + strlen(iod.output.transient.tavtemperature)];
    sprintf(avscalars[PostFcn::TEMPERATUREAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavtemperature);
  }
  if (iod.output.transient.totalpressure[0] != 0) {
    sscale[PostFcn::TOTPRESSURE] = iod.ref.rv.pressure;
    scalars[PostFcn::TOTPRESSURE] = new char[sp + strlen(iod.output.transient.totalpressure)];
    sprintf(scalars[PostFcn::TOTPRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.totalpressure);
  }
  if (iod.output.transient.tavtotalpressure[0] != 0) {
    avsscale[PostFcn::TOTPRESSUREAVG] = iod.ref.rv.pressure;
    avscalars[PostFcn::TOTPRESSUREAVG] = new char[sp + strlen(iod.output.transient.tavtotalpressure)];
    sprintf(avscalars[PostFcn::TOTPRESSUREAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavtotalpressure);
  }
  if (iod.output.transient.vorticity[0] != 0) {
    sscale[PostFcn::VORTICITY] = iod.ref.rv.velocity/iod.ref.rv.tlength;
    scalars[PostFcn::VORTICITY] = new char[sp + strlen(iod.output.transient.vorticity)];
    sprintf(scalars[PostFcn::VORTICITY], "%s%s",
            iod.output.transient.prefix, iod.output.transient.vorticity);
  }
  if (iod.output.transient.tavvorticity[0] != 0) {
    avsscale[PostFcn::VORTICITYAVG] = iod.ref.rv.velocity/iod.ref.rv.tlength;
    avscalars[PostFcn::VORTICITYAVG] = new char[sp + strlen(iod.output.transient.tavvorticity)];
    sprintf(avscalars[PostFcn::VORTICITYAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavvorticity);
  }
  if (iod.output.transient.surfaceheatflux[0] != 0) {
    scalars[PostFcn::SURFACE_HEAT_FLUX] = new char[sp + strlen(iod.output.transient.surfaceheatflux)];
    sprintf(scalars[PostFcn::SURFACE_HEAT_FLUX], "%s%s",
            iod.output.transient.prefix, iod.output.transient.surfaceheatflux);
  }
  if (iod.output.transient.tempnormalderivative[0] != 0) {
    scalars[PostFcn::TEMPERATURE_NORMAL_DERIVATIVE] = new char[sp + strlen(iod.output.transient.tempnormalderivative)];
    sprintf(scalars[PostFcn::TEMPERATURE_NORMAL_DERIVATIVE], "%s%s",
             iod.output.transient.prefix, iod.output.transient.tempnormalderivative);
  }
  if (iod.output.transient.nutturb[0] != 0) {
    sscale[PostFcn::NUT_TURB] = iod.ref.rv.viscosity_mu/iod.ref.rv.density;
    scalars[PostFcn::NUT_TURB] = new char[sp + strlen(iod.output.transient.nutturb)];
    sprintf(scalars[PostFcn::NUT_TURB], "%s%s",
            iod.output.transient.prefix, iod.output.transient.nutturb);
  }
  if (iod.output.transient.eddyvis[0] != 0) {
    sscale[PostFcn::EDDY_VISCOSITY] = iod.ref.rv.viscosity_mu;
    scalars[PostFcn::EDDY_VISCOSITY] = new char[sp + strlen(iod.output.transient.eddyvis)];
    sprintf(scalars[PostFcn::EDDY_VISCOSITY], "%s%s",
            iod.output.transient.prefix, iod.output.transient.eddyvis);
  }
  if (iod.output.transient.dplus[0] != 0) {
#if defined(HEAT_FLUX)
    /*
    double gam = iod.eqs.fluidModel.gasModel.specificHeatRatio;
    double dT = iod.bc.wall.temperature - 1.0 / (gam*(gam-1.0)*iod.bc.inlet.mach*iod.bc.inlet.mach);
    sscale[PostFcn::DELTA_PLUS] = iod.ref.reynolds * iod.eqs.thermalCondModel.prandtl / (gam * dT);
    */
    sscale[PostFcn::DELTA_PLUS] = iod.ref.rv.tpower / (iod.ref.length*iod.ref.length);
#endif
    scalars[PostFcn::DELTA_PLUS] = new char[sp + strlen(iod.output.transient.dplus)];
    sprintf(scalars[PostFcn::DELTA_PLUS], "%s%s",
            iod.output.transient.prefix, iod.output.transient.dplus);
  }


  if (iod.output.transient.philevel[0] != 0) {
    sscale[PostFcn::PHILEVEL] = 1.0;
    scalars[PostFcn::PHILEVEL] = new char[sp + strlen(iod.output.transient.philevel)];
    sprintf(scalars[PostFcn::PHILEVEL], "%s%s",
            iod.output.transient.prefix, iod.output.transient.philevel);
  }
  if (iod.output.transient.philevel2[0] != 0) {
    sscale[PostFcn::PHILEVEL2] = 1.0;
    scalars[PostFcn::PHILEVEL2] = new char[sp + strlen(iod.output.transient.philevel2)];
    sprintf(scalars[PostFcn::PHILEVEL2], "%s%s",
            iod.output.transient.prefix, iod.output.transient.philevel2);
  }
  if (iod.output.transient.fluidid[0] != 0) {
    sscale[PostFcn::FLUIDID] = 1.0;
    scalars[PostFcn::FLUIDID] = new char[sp + strlen(iod.output.transient.fluidid)];
    sprintf(scalars[PostFcn::FLUIDID], "%s%s",
            iod.output.transient.prefix, iod.output.transient.fluidid);
  }
  if (iod.output.transient.controlvolume[0] != 0) {
    sscale[PostFcn::CONTROL_VOLUME] = iod.ref.rv.length * iod.ref.rv.length * iod.ref.rv.length;
    scalars[PostFcn::CONTROL_VOLUME] = new char[sp + strlen(iod.output.transient.controlvolume)];
    sprintf(scalars[PostFcn::CONTROL_VOLUME], "%s%s",
            iod.output.transient.prefix, iod.output.transient.controlvolume);
  }
  if (iod.output.transient.velocity[0] != 0) {
    vscale[PostFcn::VELOCITY] = iod.ref.rv.velocity;
    vectors[PostFcn::VELOCITY] = new char[sp + strlen(iod.output.transient.velocity)];
    sprintf(vectors[PostFcn::VELOCITY], "%s%s",
            iod.output.transient.prefix, iod.output.transient.velocity);
  }
  if (iod.output.transient.tavvelocity[0] != 0) {
    avvscale[PostFcn::VELOCITYAVG] = iod.ref.rv.velocity;
    avvectors[PostFcn::VELOCITYAVG] = new char[sp + strlen(iod.output.transient.tavvelocity)];
    sprintf(avvectors[PostFcn::VELOCITYAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavvelocity);
  }
  if (iod.output.transient.displacement[0] != 0) {
    vscale[PostFcn::DISPLACEMENT] = iod.ref.rv.tlength;
    vectors[PostFcn::DISPLACEMENT] = new char[sp + strlen(iod.output.transient.displacement)];
    sprintf(vectors[PostFcn::DISPLACEMENT], "%s%s",
            iod.output.transient.prefix, iod.output.transient.displacement);
  }
  if (iod.output.transient.tavdisplacement[0] != 0) {
    avvscale[PostFcn::DISPLACEMENTAVG] = iod.ref.rv.tlength;
    avvectors[PostFcn::DISPLACEMENTAVG] = new char[sp + strlen(iod.output.transient.tavdisplacement)];
    sprintf(avvectors[PostFcn::DISPLACEMENTAVG], "%s%s",
            iod.output.transient.prefix, iod.output.transient.tavdisplacement);
  }

  if (iod.output.transient.flightDisplacement[0] != 0) {
    vscale[PostFcn::FLIGHTDISPLACEMENT] = iod.ref.rv.tlength;
    vectors[PostFcn::FLIGHTDISPLACEMENT] = new char[sp + strlen(iod.output.transient.flightDisplacement)];
    sprintf(vectors[PostFcn::FLIGHTDISPLACEMENT], "%s%s",
            iod.output.transient.prefix, iod.output.transient.flightDisplacement);
  }

  if (iod.output.transient.localFlightDisplacement[0] != 0) {
    vscale[PostFcn::LOCALFLIGHTDISPLACEMENT] = iod.ref.rv.tlength;
    vectors[PostFcn::LOCALFLIGHTDISPLACEMENT] = new char[sp + strlen(iod.output.transient.localFlightDisplacement)];
    sprintf(vectors[PostFcn::LOCALFLIGHTDISPLACEMENT], "%s%s",
            iod.output.transient.prefix, iod.output.transient.localFlightDisplacement);
  }

  if (iod.output.transient.velocitynorm[0] != 0) {
    sscale[PostFcn::VELOCITY_NORM] = iod.ref.rv.velocity;
    scalars[PostFcn::VELOCITY_NORM] = new char[sp + strlen(iod.output.transient.velocitynorm)];
    sprintf(scalars[PostFcn::VELOCITY_NORM], "%s%s",
            iod.output.transient.prefix, iod.output.transient.velocitynorm);
    sprintf(dScalars[PostFcn::DERIVATIVE_TOTPRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.dTotalpressure);
    sprintf(dScalars[PostFcn::DERIVATIVE_TOTPRESSURE], "%s%s",
            iod.output.transient.prefix, iod.output.transient.dTotalpressure);
  }

  int dsp = strlen(iod.output.transient.prefix) + 1;

  if (iod.output.transient.dNutturb[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_NUT_TURB] = iod.ref.rv.viscosity_mu/iod.ref.rv.density;
    dScalars[PostFcn::DERIVATIVE_NUT_TURB] = new char[dsp + strlen(iod.output.transient.dNutturb)];
    sprintf(dScalars[PostFcn::DERIVATIVE_NUT_TURB], "%s%s",
            iod.output.transient.prefix, iod.output.transient.dNutturb);
  }

  if (iod.output.transient.dEddyvis[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_EDDY_VISCOSITY] = iod.ref.rv.viscosity_mu;
    dScalars[PostFcn::DERIVATIVE_EDDY_VISCOSITY] = new char[dsp + strlen(iod.output.transient.dEddyvis)];
    sprintf(dScalars[PostFcn::DERIVATIVE_EDDY_VISCOSITY], "%s%s",
            iod.output.transient.prefix, iod.output.transient.dEddyvis);
  }

  if (iod.output.transient.dVelocityScalar[0] != 0) {
    dSscale[PostFcn::DERIVATIVE_VELOCITY_SCALAR] = iod.ref.rv.velocity;
    dScalars[PostFcn::DERIVATIVE_VELOCITY_SCALAR] = new char[dsp + strlen(iod.output.transient.dVelocityScalar)];
    sprintf(dScalars[PostFcn::DERIVATIVE_VELOCITY_SCALAR], "%s%s",
            iod.output.transient.prefix, iod.output.transient.dVelocityScalar);
  }

  if (iod.output.transient.dVelocityVector[0] != 0) {
    dVscale[PostFcn::DERIVATIVE_VELOCITY_VECTOR] = iod.ref.rv.velocity;
    dVectors[PostFcn::DERIVATIVE_VELOCITY_VECTOR] = new char[dsp + strlen(iod.output.transient.dVelocityVector)];
    sprintf(dVectors[PostFcn::DERIVATIVE_VELOCITY_VECTOR], "%s%s",
            iod.output.transient.prefix, iod.output.transient.dVelocityVector);
  }

  if (iod.output.transient.dDisplacement[0] != 0) {
    dVscale[PostFcn::DERIVATIVE_DISPLACEMENT] = iod.ref.rv.tlength;
    dVectors[PostFcn::DERIVATIVE_DISPLACEMENT] = new char[dsp + strlen(iod.output.transient.dDisplacement)];
    sprintf(dVectors[PostFcn::DERIVATIVE_DISPLACEMENT], "%s%s",
            iod.output.transient.prefix, iod.output.transient.dDisplacement);
  }

}

//------------------------------------------------------------------------------


template<int dim>
int TsOutput<dim>::writeBinaryVectorsToDiskRom(bool lastNewtonIt, int timeStep, int newtonIt,
                                                  DistSVec<double,dim> *state, DistSVec<double,dim> *residual)
{ // Outputs state and residual snapshots (from the FOM newton solver) for building a nonlinear ROM.
  // The logic tests ensure that every state from a given timestep is ouput (if requested), but that
  // no state snapshot is stored twice. (Need to be careful because the initial state of any given
  // timestep is the converged state from the previous timestep)

  int status = 0;

  timeStep = timeStep - 1; //need the counting to start at 0, not 1

  if (state && stateVectors) {
    double tag = *(domain->getNewtonTag());
    double prevTag = -1;
    int step = *(domain->getNewtonStateStep());
    int prevStep = step-1;
    int dummy;
    if (prevStep>=0) domain->readTagFromFile<double,dim>(stateVectors, prevStep, &prevTag, &dummy);
    if ((newtonIt==0) && (step!=0) && (prevTag==tag)) { //do nothing
    } else if (timeStep%stateOutputFreqTime==0) {
        if (((stateOutputFreqNewton==0)&&lastNewtonIt) ||  // special case: last newton iteration
            ((timeStep==0) && (newtonIt==0)) ||  // special case: initial condition
            ((stateOutputFreqNewton>0)&&(newtonIt%stateOutputFreqNewton==0))) {
          // output FOM state
          domain->writeVectorToFile(stateVectors, step, tag, *state);
          ++(*(domain->getNewtonStateStep()));
          ++status;
        }
    }
  }

  if (residual && residualVectors) {
    if (timeStep%residualOutputFreqTime==0) {
      if (residualOutputMaxNewton>=newtonIt) {
        // for FOM residuals only (residuals from PG are clustered during the online simulations)

        // if outputting krylov vects, limit number of residuals output per newton iteration to
        // number of krylov vecs output at previous it
        if ((fdResiduals && fdResidualsLimit) && (*(domain->getNumKrylovVecsOutputPrevNewtonIt())>0) &&
            (*(domain->getNumResidualsOutputCurrentNewtonIt()) >= *(domain->getNumKrylovVecsOutputPrevNewtonIt())))  return status;

        domain->writeVectorToFile(residualVectors, *(domain->getNewtonResidualStep()), *(domain->getNewtonTag()), *residual);
        ++(*(domain->getNewtonResidualStep()));
        ++(*(domain->getNumResidualsOutputCurrentNewtonIt()));
        status = (status==1) ? 3 : 2;
      }
    }
  }

  return status;
}

/**
 * output state and mask from a embedded full order model for training a ROM.
 * This is called from EmbeddedTsDesc::outputToDisk(); the frequency of saving is set in
 * IoData::RomOutputData.
 */
template<int dim>
void TsOutput<dim>::writeStateMaskVectorsToDiskRom(int timestep, DistSVec<double, dim> &state, DistSVec<char, dim> &mask){
  if (timestep == 0 || timestep % stateOutputFreqTime == 0) {
  // use a global variable to increment it
    int step = *(domain->getTimeIt());
  if(stateVectors)
      domain->writeVectorToFile(stateVectors, step, *(domain->getNewtonTag()), state);
  // saving mask, use overloaded writeVectorToFile function
  if(stateMaskVectors)
      domain->writeVectorToFile(stateMaskVectors, step, *(domain->getNewtonTag()), mask);
    ++(*(domain->getTimeIt()));
  }
}

//------------------------------------------------------------------------------
