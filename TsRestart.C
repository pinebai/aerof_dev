#include <cstdio>
#include <fstream>
#include <iostream>

#include <TsRestart.h>

#include <IoData.h>
#include <RefVal.h>
#include <DistGeoState.h>
#include <DistTimeState.h>
#include "DistVector.h"
#include "Domain.h"
#include <LevelSet.h>
#include <MeshMotionHandler.h>
#include "PostOperator.h"
#include "SubDomain.h"
#include <FSI/DynamicNodalTransfer.h>

//------------------------------------------------------------------------------

template<int dim, int dimLS>
void TsRestart::writeToDisk(int cpuNum, bool lastIt, int it, double t, double dt,
			    DistTimeState<dim> &timeState, DistGeoState &geoState,
			    LevelSet<dimLS> *levelSet, 
                            DynamicNodalTransfer* dyn, // to output cracking information
			    FluidSelector* fluidSelector
                            )
{
  iteration = it;
  etime = t;
  double dt_nm1 = timeState.getData().dt_nm1;
  double dt_nm2 = timeState.getData().dt_nm2;

  if (toWrite(iteration, lastIt, etime)) {

    if (lastIt) 
      index = 0;

    if (solutions[index][0] != 0)
      timeState.writeToDisk(solutions[index]);
    if (positions[index][0] != 0)
      geoState.writeToDisk(positions[index]);
    if (levelsets[index][0] != 0 && levelSet)
      levelSet->writeToDisk(levelsets[index]);
    
    if (fluidId[index][0] != 0 && fluidSelector)
      fluidSelector->writeToDisk(fluidId[index]);
      
/* PJSA moved to separate function writeCrackingDataToDisk

    if (cpuNum == 0 && dyn) {

      std::ofstream ofile(cracking[index],std::ios::binary);
      dyn->writeCrackingData(ofile);
    }
*/  
    if (cpuNum == 0 && data[index][0] != 0) {
      FILE *fp = fopen(data[index], "w");
      if (!fp) {
	fprintf(stderr, "*** Error: could not open \'%s\'\n", data[index]);
	exit(1);
      }
      
      if (refVal->mode == RefVal::NON_DIMENSIONAL)
	fprintf(fp, "// Restart file (values are non dimensional)\n\n");
      else
	fprintf(fp, "// Restart file (values are dimensional)\n\n");
      fprintf(fp, "under RestartParameters {\n");
      fprintf(fp, "  Iteration = %d;\n", iteration);
      fprintf(fp, "  Time = %e;\n", etime * refVal->time);
      fprintf(fp, "  TimeStep1 = %e;\n", dt_nm1 * refVal->time);
      fprintf(fp, "  TimeStep2 = %e;\n", dt_nm2 * refVal->time);
      fprintf(fp, "  Residual = %e;\n", residual);
      fprintf(fp, "  Energy = %e;\n", energy[0] * refVal->energy);
      fprintf(fp, "  NewtonOutputTag = %e;\n", timeState.getNewtonTag());
      fprintf(fp, "  NewtonStateOutputStep = %d;\n", timeState.getNewtonStateStep());
      fprintf(fp, "  NewtonResidualOutputStep = %d;\n", timeState.getNewtonResidualStep());
      fprintf(fp, "  KrylovOutputStep = %d;\n", timeState.getKrylovStep());
      fprintf(fp, "}\n");
      fclose(fp);
    }

    if (index == 1)
      index = 2;
    else if (index == 2)
      index = 1;

  }

}

//------------------------------------------------------------------------------




template<int dim>
void TsRestart::writeKPtracesToDisk
(
  IoData &iod,
  bool lastIt, int it, double t, 
  DistSVec<double,3> &X,
  DistVec<double> &A,
  DistSVec<double,dim> &U,
  DistTimeState<dim> *timeState,
  Domain *domain,
  PostOperator<dim> *postOp
)
{

  if (writeKPtraces == false)
    return;

  double time = t * refVal->time;
  int step = 0;
  if (iod.output.restart.frequency > 0)
  {
    step = it / iod.output.restart.frequency;      
    if (it % iod.output.restart.frequency != 0)
    {
      if (lastIt)
      {
        step += 1;
      }
      else
      {
        return;
      }
    }
  }
  else
  {
    step = it;
  }

  DistVec<double> Qs(domain->getNodeDistInfo());

  char *prefix = new char[strlen(iod.output.restart.prefix) + 1 + strlen(iod.output.restart.strKPtraces)];

  SubDomain **subDomain = domain->getSubDomain();

  struct Data {
    int nid;
    double value;
  };

  sprintf(prefix, "%s%s", iod.output.restart.prefix, iod.output.restart.strKPtraces);

  Qs = (double) 0.0;
  postOp->computeScalarQuantity(PostFcn::PRESSURE, X, U, A, Qs, timeState);

#pragma omp parallel for
  for (int iSub = 0; iSub < domain->getNumLocSub(); ++iSub) 
  {
    char filename[strlen(prefix)+2];
    sprintf(&filename[0], "%s_sub%d", prefix, subDomain[iSub]->getGlobSubNum());
    std::set<int> nList = subDomain[iSub]->getKirchhoffNodesList();
    ofstream myFile;
    if (step == 0)
    {
      myFile.open(&filename[0], ios::binary);
      Data x; x.nid = 0; x.value = (double) nList.size();
      myFile.write((char*) &x, sizeof(Data));
    }
    else 
    {
      myFile.open(&filename[0], ios::binary | ios::app);
    }
    myFile.write((char*) &time, sizeof(double));
    int nlSize = nList.size();
    if (nlSize > 0)
    {
      Data *myData = new Data[nlSize];
      int count = 0;
      for (std::set<int>::iterator it = nList.begin(); it != nList.end(); ++it)
      {
        myData[count].value = (Qs(iSub))[*it] * iod.ref.rv.pressure;
        myData[count].nid = *it;
        count += 1;
      } 
      myFile.write((char*) myData, sizeof(Data)*nlSize);
      myFile.flush();
      free(myData);
    }
    //
    myFile.flush();
    myFile.close(); 
    //
    if (lastIt == true)
    {
      Data y; y.nid = step + 1; y.value = (double) nList.size();
      myFile.open(&filename[0], ios::binary | ios::in | ios::out);
      myFile.seekp((long) 0, ios::beg);
      myFile.write((char*) &y, sizeof(Data));
      myFile.flush();
      myFile.close();
    }
    //
  } // for (int iSub = 0; iSub < domain->getNumLocSub(); ++iSub)

  Communicator *MyCom_p = domain->getCommunicator();
  if ((iod.output.restart.frequency > 0) && (MyCom_p->cpuNum() == 0))
  {
    // std::cout << "Wrote Kirchhoff trace " << step << " to '" << prefix << "'\n";
    if (lastIt == true)
    {
      std::cout << "Kirchhoff traces have been written to " << prefix << "'\n";
      std::cout << " --> Total Number of Snapshots N = " << step + 1 << std::endl;
      double TT = time;
      if (iod.output.restart.frequency > 0)
        TT += (time/it) * iod.output.restart.frequency;
      else
        TT += (time/it);
      std::cout << " --> Pseudo-period T = " << TT << std::endl;
      std::cout << " --> The discrete Fourier frequencies will be of the form 2*PI*j/T with j in [0, N-1].\n";
    }
  }

  delete[] prefix;

}
