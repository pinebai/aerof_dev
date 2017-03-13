
#ifdef USE_MPI
#include <mpi.h>
#endif

#include <cstring>

#include <TsRestart.h>

#include <IoData.h>
#include <RefVal.h>

//------------------------------------------------------------------------------

TsRestart::TsRestart(IoData &iod, RefVal *rv) : refVal(rv)
  , writeKPtraces(false)
{

  int sp = strlen(iod.output.restart.prefix) + 1;

  solutions[0] = new char[sp + strlen(iod.output.restart.solutions)];
  if (iod.output.restart.solutions[0] != 0)
    sprintf(solutions[0], "%s%s", iod.output.restart.prefix, iod.output.restart.solutions);
  else
    sprintf(solutions[0], "");

  positions[0] = new char[sp + strlen(iod.output.restart.positions)];
  if (iod.output.restart.positions[0] != 0) // && (iod.problem.framework == ProblemData::BODYFITTED || iod.problem.framework == ProblemData::EMBEDDEDALE))
    sprintf(positions[0], "%s%s", iod.output.restart.prefix, iod.output.restart.positions);
  else
    sprintf(positions[0], "");

  levelsets[0] = new char[sp + strlen(iod.output.restart.levelsets)];
  if (iod.output.restart.levelsets[0] != 0)
    sprintf(levelsets[0], "%s%s", iod.output.restart.prefix, iod.output.restart.levelsets);
  else
    sprintf(levelsets[0], "");

  fluidId[0] = new char[sp + strlen(iod.output.restart.fluidId)];
  if (iod.output.restart.fluidId[0] != 0)
    sprintf(fluidId[0], "%s%s", iod.output.restart.prefix, iod.output.restart.fluidId);
  else
    sprintf(fluidId[0], "");


  cracking[0] = new char[sp + strlen(iod.output.restart.cracking)];
  if (iod.output.restart.cracking[0] != 0)
    sprintf(cracking[0], "%s%s", iod.output.restart.prefix, iod.output.restart.cracking);
  else
    sprintf(cracking[0], "");


  data[0] = new char[sp + strlen(iod.output.restart.data)];
  if (iod.output.restart.data[0] != 0)
    sprintf(data[0], "%s%s", iod.output.restart.prefix, iod.output.restart.data);
  else
    sprintf(data[0], "");

  structPos = new char[sp + strlen(iod.output.restart.embeddedpositions)];
  if (iod.output.restart.embeddedpositions[0] != 0 && (iod.problem.framework == ProblemData::EMBEDDED || iod.problem.framework == ProblemData::EMBEDDEDALE) )
    sprintf(structPos, "%s%s", iod.output.restart.prefix, iod.output.restart.embeddedpositions);
  else
    sprintf(structPos, "");

  char restart_file_package[256];
  if (iod.output.restart.filepackage[0] != 0) {
    sprintf(restart_file_package, "%s%s", iod.output.restart.prefix, iod.output.restart.filepackage);
 
    writeRestartFileNames(restart_file_package);
  }

  if (iod.output.restart.type == RestartData::SINGLE) {
    deleteCharStar = false;
    for (int i=1; i<3; ++i) {
      solutions[i] = solutions[0];
      positions[i] = positions[0];
      levelsets[i] = levelsets[0];
      cracking[i] = cracking[0];
      fluidId[i] = fluidId[0];
      data[i] = data[0];

    }

    index = 0;
  }
  else {
    deleteCharStar = true;
    for (int i=1; i<3; ++i) {
      solutions[i] = new char[strlen(solutions[0]) + 5 + 1];
      if (solutions[0][0] != 0)
	sprintf(solutions[i], "%s.%drst", solutions[0], i);
      else
	sprintf(solutions[i], "");

      positions[i] = new char[strlen(positions[0]) + 5 + 1];
      if (positions[0][0] != 0) 
	sprintf(positions[i], "%s.%drst", positions[0], i);
      else
	sprintf(positions[i], "");

      levelsets[i] = new char[strlen(levelsets[0]) + 5 + 1];
      if (levelsets[0][0] != 0)
	sprintf(levelsets[i], "%s.%drst", levelsets[0], i);
      else
	sprintf(levelsets[i], "");

      cracking[i] = new char[strlen(cracking[0]) + 5 + 1];
      if (cracking[0][0] != 0)
	sprintf(cracking[i], "%s.%drst", cracking[0], i);
      else
	sprintf(cracking[i], "");

      fluidId[i] = new char[strlen(fluidId[0]) + 5 + 1];
      if (fluidId[0][0] != 0)
	sprintf(fluidId[i], "%s.%drst", fluidId[0], i);
      else
	sprintf(fluidId[i], "");


      data[i] = new char[strlen(data[0]) + 5 + 1];
      if (data[0][0] != 0)
	sprintf(data[i], "%s.%drst", data[0], i);
      else
	sprintf(data[i], "");

    }

    index = 1;
  }

  iteration = iod.restart.iteration;
  etime = iod.restart.etime;
  residual = iod.restart.residual;
  energy[0] = iod.restart.energy;
  energy[1] = iod.restart.energy;
  frequency = iod.output.restart.frequency;
  frequency_dt = iod.output.restart.frequency_dt;
  prtout = 0.0;

  //
  // Check whether pressure snapshots are needed.
  // UH (07/2012)
  //
  if (strlen(iod.output.restart.strKPtraces) > 0)
  {
    writeKPtraces = true;
  }

}

void TsRestart::writeRestartFileNames(const char* fn) {

  FILE* file = fopen(fn, "w");
  if (!file) { 
    fprintf(stderr, "*** Error: could not open \'%s\'\n", fn);
    exit(1);
  }
  fprintf(file,"%s\n",solutions[0]);
  fprintf(file,"%s\n",positions[0]);
  fprintf(file,"%s\n",levelsets[0]);
  fprintf(file,"%s\n",cracking[0]);
  fprintf(file,"%s\n",fluidId[0]);
  fprintf(file,"%s\n",data[0]);

  fprintf(file,"%s\n",structPos); 

  fclose(file); 
  
}

void TsRestart::readRestartFileNames(const char* fn,
				     char* sols,
				     char* posit,
				     char* ls,
				     char* crk,
				     char* fid,
				     char* dat,
				     char* spos,
                                     Communicator* com) {

  /*std::ifstream infile(fn);

  if (!infile.good()) {

    std::cout << "Error: Cannot read package file  " << fn << std::endl;
    exit(-1);
  }

  infile.getline(sols,256);
  infile.getline(posit,256);
  infile.getline(ls,256);
  infile.getline(crk,256);
  infile.getline(fid,256);
  infile.getline(dat,256);
  infile.getline(spos,256);
 */

  char tmp[7][256];
  FILE* fin;

  if (com == NULL || com->cpuNum() == 0) {
    fin = fopen(fn,"r");

    if (!fin) {

      std::cout << "Could not read restart package file " << fn << std::endl;
      exit(-1);
    }
  }
   
  for (int i = 0; i < 7; ++i) {

    if (com == NULL || com->cpuNum() == 0)
      fscanf(fin,"%s",tmp[i]); 

    if (com)
      com->broadcast(256, tmp[i]);
  }
    


  if (com == NULL || com->cpuNum() == 0) {
    fclose(fin); 
  }

  strcpy(sols, tmp[0]);
  strcpy(posit, tmp[1]);
  strcpy(ls, tmp[2]);
  strcpy(crk, tmp[3]);
  strcpy(fid, tmp[4]);
  strcpy(dat, tmp[5]);
  strcpy(spos, tmp[6]);
}


//------------------------------------------------------------------------------

bool TsRestart::toWrite(int it, bool lastIt, double t)
{
  if(frequency_dt<=0.0)
    return (((frequency > 0) && (it % frequency == 0)) || lastIt);

  return (t>=prtout || lastIt);
}

//------------------------------------------------------------------------------

void TsRestart::updatePrtout(double t)
{
  if(frequency_dt<=0.0)
    return;
  if(t>=prtout)
    prtout += frequency_dt;
}

//------------------------------------------------------------------------------

void TsRestart::writeStructPosToDisk(int cpuNum, bool lastIt, Vec<Vec3D>& Xs)
{
  
  if(cpuNum>0) return; //only Proc.#1 will work.

  if (structPos[0]!=0 && toWrite(iteration, lastIt, etime)) {
  //if ((lastIt || (frequency > 0 && iteration % frequency == 0)) && structPos[0]!=0) {
    FILE *fp = fopen(structPos,"w");
    if (!fp) {
      fprintf(stderr, "*** Error: could not open \'%s\'\n", structPos);
      exit(1);
    }

    for (int i=0; i<Xs.size(); i++)
      fprintf(fp,"%d %e %e %e\n", i+1, Xs[i][0], Xs[i][1], Xs[i][2]);

   fclose(fp);
  }
}

//------------------------------------------------------------------------------

// Included (MB)
void TsRestart::rstVar(IoData &ioData) {

  etime = ioData.restart.etime;
  energy[0] = ioData.restart.energy;
  energy[1] = ioData.restart.energy;

}

//------------------------------------------------------------------------------

void TsRestart::writeCrackingDataToDisk(int cpuNum, bool lastIt, int it, double t,
                                        DynamicNodalTransfer* dyn // to output cracking information
                                       )
{
  if (toWrite(it, lastIt, t)) {

    if (cpuNum == 0 && dyn) {

      std::ofstream ofile(cracking[index],std::ios::binary);
      dyn->writeCrackingData(ofile);
    }
  }
}
