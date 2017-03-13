#include <cstring>
#include <complex>
using std::complex;
#include <TsParameters.h>
#include <IoData.h>

#include <cmath>

#ifdef OLD_STL
#include <algo.h>
#else
#include <algorithm>
using std::min;
using std::max;
#endif
#include <cstring>

#ifndef PI
#define PI 3.14159
#endif

//------------------------------------------------------------------------------

TsParameters::TsParameters(IoData &ioData)
{

  resType = ioData.ts.residual;

  cfllaw = ioData.ts.cfl.strategy;

  cfl0 = ioData.ts.cfl.cfl0;
  cflCoef1 = ioData.ts.cfl.cflCoef1;
  cflCoef2 = ioData.ts.cfl.cflCoef2;
  cflMax = ioData.ts.cfl.cflMax;
  cflMin = ioData.ts.cfl.cflMin;  // not used... cflCoef1 is the real cflMin
  dualtimecfl = ioData.ts.cfl.dualtimecfl;

  checksol = ioData.ts.checksol;
  checklinsolve = ioData.ts.checklinsolve;
  checkriemann = checksol;
  checklargevelocity = ioData.ts.checkvelocity;
  checkpclipping = ioData.ts.checkpressure;
  checkrhoclipping = ioData.ts.checkdensity;
  rapidpchangecutoff = max(0,ioData.ts.deltapressurethreshold);
  rapiddchangecutoff = max(0,ioData.ts.deltadensitythreshold);

  ser = ioData.ts.cfl.ser;
  angle_growth = ioData.ts.cfl.angle_growth;
  angle_zero = ioData.ts.cfl.angle_zero;
  dft_history = ioData.ts.cfl.dft_history;
  dft_freqcutoff = ioData.ts.cfl.dft_freqcutoff;
  dft_growth = ioData.ts.cfl.dft_growth;
  fixedunsteady_counter = 0;

  maxIts = ioData.ts.maxIts;
  eps = ioData.ts.eps;
  epsabs = ioData.ts.epsabs; 
  maxTime = ioData.ts.maxTime;
  dt_imposed = ioData.ts.timestep;

  output = new char[strlen(ioData.ts.output) + 1];
  sprintf(output, "%s", ioData.ts.output);

  cfl = cfl0;
  residual = 1.0;

  reshistory = new double[dft_history];
  for (int i=0; i<dft_history; i++) reshistory[i] = 1.0;
  dft = new complex<double>[dft_history];
  for (int i=0; i<dft_history; i++) dft[i]=0.0;

  allowstop = true;

  errorHandler = NULL;

  if (strcmp(ioData.ts.cfl.output, "") == 0)
    cfl_output = 0;
  else if (strcmp(ioData.ts.cfl.output, "stdout") == 0)
    cfl_output = stdout;
  else if (strcmp(ioData.ts.cfl.output, "stderr") == 0)
    cfl_output = stderr;
  else {
    cfl_output = fopen(ioData.ts.cfl.output, "w");
    if (!cfl_output) {
      fprintf(stderr, "*** Error: could not open \'%s\'\n", ioData.ts.cfl.output);
      exit(1);
    }
  }

  nonlinearRomPostpro = (ioData.problem.alltype == ProblemData::_STEADY_NONLINEAR_ROM_POST_
                         || ioData.problem.alltype == ProblemData::_UNSTEADY_NONLINEAR_ROM_POST_);

}

//------------------------------------------------------------------------------

TsParameters::~TsParameters()
{

  if (output) delete [] output;
  if (reshistory) delete [] reshistory;
  if (dft) delete [] dft;
  if (cfl_output) fclose(cfl_output);
}

//------------------------------------------------------------------------------

//Figure out what to do with errors (related to time step)
void TsParameters::resolveErrors()
{
  if(cfllaw == CFLData::OLD && dt_imposed < 0) return; // for backwards compatilibity
  if(nonlinearRomPostpro) return;

  if (checklinsolve && errorHandler->globalErrors[ErrorHandler::SATURATED_LS]) {
    errorHandler->com->printf(1,"Detected saturated linear solver. Reducing time-step.\n");
    errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP] += 1;
    errorHandler->globalErrors[ErrorHandler::SATURATED_LS] = 0;
    errorHandler->localErrors[ErrorHandler::SATURATED_LS] = 0;
  }

  if (checksol && errorHandler->globalErrors[ErrorHandler::UNPHYSICAL]) {
    errorHandler->com->printf(1,"Detected unphysical solution. Reducing time-step.\n");
    errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP] += 1;
    errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP] += 1;
    errorHandler->globalErrors[ErrorHandler::UNPHYSICAL] = 0;
    errorHandler->localErrors[ErrorHandler::UNPHYSICAL] = 0;
  }
  
  if (checkriemann && errorHandler->globalErrors[ErrorHandler::BAD_RIEMANN]) {
    errorHandler->com->printf(1,"Detected error in Riemann solver. Reducing time-step.\n");
    errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP] += 1;
    errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP] += 1;
    errorHandler->globalErrors[ErrorHandler::BAD_RIEMANN] = 0;
    errorHandler->localErrors[ErrorHandler::BAD_RIEMANN] = 0;
  }

  if (checkpclipping && errorHandler->globalErrors[ErrorHandler::PRESSURE_CLIPPING]) {
    errorHandler->com->printf(1,"Detected pressure clipping. Reducing time-step. "
    "If the problem is expected to produce cavitation, set CheckPressure to Off.\n");
    errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP] += 1;
    errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP] += 1;
    errorHandler->globalErrors[ErrorHandler::PRESSURE_CLIPPING] = 0;
    errorHandler->localErrors[ErrorHandler::PRESSURE_CLIPPING] = 0;
  }

  if (checkrhoclipping && errorHandler->globalErrors[ErrorHandler::DENSITY_CLIPPING]) {
    errorHandler->com->printf(1,"Detected density clipping. Reducing time-step. "
    "If the problem is expected to produce cavitation, set CheckDensity to Off.\n");
    errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP] += 1;
    errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP] += 1;
    errorHandler->globalErrors[ErrorHandler::DENSITY_CLIPPING] = 0;
    errorHandler->localErrors[ErrorHandler::DENSITY_CLIPPING] = 0;
  }

  if (checklargevelocity && errorHandler->globalErrors[ErrorHandler::LARGE_VELOCITY]) {
    errorHandler->com->printf(1,"Detected abnormally large velocities. Reducing time-step.\n");
    errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP] += 1;
    errorHandler->globalErrors[ErrorHandler::REDO_TIMESTEP] += 1;
    errorHandler->globalErrors[ErrorHandler::LARGE_VELOCITY] = 0;
    errorHandler->localErrors[ErrorHandler::LARGE_VELOCITY] = 0;
  }

  if (rapidpchangecutoff && errorHandler->globalErrors[ErrorHandler::RAPIDLY_CHANGING_PRESSURE] >= rapidpchangecutoff) {
    errorHandler->com->printf(1,"Detected multiple rapidly changing pressures. Reducing time-step.\n");
    errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP] += 1;
    errorHandler->globalErrors[ErrorHandler::RAPIDLY_CHANGING_PRESSURE] = 0;
    errorHandler->localErrors[ErrorHandler::RAPIDLY_CHANGING_PRESSURE] = 0;
  }

  if (rapidpchangecutoff && errorHandler->globalErrors[ErrorHandler::RAPIDLY_CHANGING_DENSITY] >= rapiddchangecutoff) {
    errorHandler->com->printf(1,"Detected multiple rapidly changing density. Reducing time-step.\n");
    errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP] += 1;
    errorHandler->globalErrors[ErrorHandler::RAPIDLY_CHANGING_DENSITY] = 0;
    errorHandler->localErrors[ErrorHandler::RAPIDLY_CHANGING_DENSITY] = 0;
  }

  errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP_TIME] = errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP];
}

//------------------------------------------------------------------------------

void TsParameters::computeCflNumber(int its, double res, double angle)
{
  // First run automatic CFL checks


  if(errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP]) {
    errorHandler->globalErrors[ErrorHandler::REDUCE_TIMESTEP] = 0;
    double cflold = cfl;
    cfl *= 0.5;
    fixedunsteady_counter = 1;
    errorHandler->com->printf(1,"Reducing CFL. Previous cfl=%e, new cfl=%e.\n",cflold,cfl);
    if(cfl < cfl0/1000. && allowstop) {std::printf("Cannot further reduce CFL number. Aborting.\n"); exit(-1);}
    return;
  }

  // for backwards compatibility
  if (cfllaw == CFLData::OLD){
    cfl = min( max( max(cflCoef1, min(its*cflCoef2,cfl+cflCoef2)), cfl0/pow(res,ser) ), cflMax );
    return;
  }

  // Now run main code
  double cfl_res, cfl_dir, cfl_dft, hf_ratio;

  double cfl_prev = cfl;

  for (int i=dft_history-1; i>0; i--) reshistory[i]=reshistory[i-1];
  reshistory[0]=res;

  if (cfllaw == CFLData::RESIDUAL || cfllaw == CFLData::HYBRID){
    // compute residual strategy proposal
    cfl = cfl_prev;
    cfl *= pow(reshistory[1]/reshistory[0],ser);
    cfl_res = cfl;
    if(cfl_output) errorHandler->com->fprintf(cfl_output,"CFL residual strategy: old residual: %e, new residual: %e, CFL proposal: %e\n",reshistory[1],reshistory[0],cfl_res);
  }
  if (cfllaw == CFLData::DIRECTION || cfllaw == CFLData::HYBRID){
    // compute direction strategy proposal
    cfl = cfl_prev;
    if (angle != -2.0) cfl *= pow(angle_growth,angle-angle_zero);
    cfl_dir = cfl;
    if(cfl_output) errorHandler->com->fprintf(cfl_output,"CFL direction strategy: angle: %e, CFL proposal: %e\n",angle,cfl_dir);
  }
  if (cfllaw == CFLData::DFT || cfllaw == CFLData::HYBRID){
    // compute dft strategy proposal
    cfl = cfl_prev;

    if (reshistory[0] != 0.0){
      double e_hf, e_ac, e_dc, e_total;
      int cutofflow = dft_history/2 - (dft_freqcutoff-1)/2; 
      int cutoffhigh = (dft_history+1)/2 + (dft_freqcutoff-1)/2;

      for (int i=cutofflow; i<cutoffhigh+1; i++){
	dft[i]=0.0;
	for (int j=0; j<dft_history; j++) 
	  dft[i] += 1/sqrt((double)dft_history)*complex<double>(reshistory[j],0.0)*exp(-2.0*complex<double>(0.0,1.0)*PI*(double)i*(double)j/(double)dft_history);
      }
      e_hf = 0.0;
      for (int i=cutofflow; i<cutoffhigh+1; i++)
	e_hf += norm(dft[i]);

      e_total = 0.0;
      for (int i=0; i<dft_history; i++) e_total += pow(reshistory[i],2);
      e_dc = 0.0;
      for (int i=0; i<dft_history; i++) e_dc += reshistory[i]/sqrt((double)dft_history);
      e_dc *= e_dc;
      e_ac = e_total-e_dc;

      if (e_ac<e_dc*1e-5) hf_ratio = 1.0;
      else hf_ratio = e_hf/e_ac;

      if (hf_ratio<0 || hf_ratio>1) {
	std::printf("Found invalid hf_ratio: %e, check for bugs in CFL Law\n",hf_ratio); 
        char buffer[2500];
        sprintf(buffer,"Res history:");
        for (int i=0; i<dft_history; i++) sprintf(buffer,"%s %e",buffer,reshistory[i]);
        sprintf(buffer,"%s \n",buffer);
        sprintf(buffer,"%s DFT:",buffer);
        for (int i=0; i<dft_history; i++) sprintf(buffer,"%s %e+%ei", buffer,dft[i].real(),dft[i].imag());
        sprintf(buffer,"%s \n",buffer);
        errorHandler->com ->fprintf(stderr,"%s e_total=%e, e_dc=%e, e_ac=%e, e_hf=%e, hf_ratio=%e\n",buffer,e_total,e_dc,e_ac,e_hf,hf_ratio);
	exit(-1);
      }

      if(its != 0) cfl *= pow(dft_growth, 1-2*hf_ratio);
      if(cfl_output) errorHandler->com->fprintf(cfl_output,"CFL DFT strategy: e_ac: %e, e_hf: %e, CFL proposal: %e\n",e_ac,e_hf,cfl);
    }
    cfl_dft = cfl;
  }
  if (cfllaw == CFLData::HYBRID){
    // compute hybrid strategy cfl
    if (hf_ratio > 0.66) cfl = cfl_dft;
    else if (angle < angle_zero) cfl = cfl_dir;
    else cfl = max(cfl_dir, cfl_res);
    if(cfl_output) errorHandler->com->fprintf(cfl_output,"CFL Hybrid strategy: dft_proposal: %e, direction proposal: %e, residual proposal: %e, chosen cfl: %e\n",cfl_dft,cfl_dir,cfl_res,cfl);
  }
  if (cfllaw == CFLData::FIXEDUNSTEADY){
    // compute fixed unsteady cfl law
    // attempt to keep cfl fixed, if automatic reductions lower the CFL, re-increase it after a few iterations
    cfl = cfl_prev;
    if(fixedunsteady_counter >= 4){
      fixedunsteady_counter = 2;
      cfl *= 2.0;
      cfl = min(cfl,cfl0);
    }
    fixedunsteady_counter++;
  }

  cfl = (min(max(cfl, cflCoef1),cflMax));
  if(cfl_output) errorHandler->com->fprintf(cfl_output,"CFL number chosen: %e. cfl_prev=%e\n",cfl,cfl_prev);
  return;
}

//------------------------------------------------------------------------------

// Included (MB)
void TsParameters::rstVar(IoData &ioData) {

  maxTime = ioData.ts.maxTime;

}

//------------------------------------------------------------------------------

