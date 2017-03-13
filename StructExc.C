#include <StructExc.h>

#include <IoData.h>
#include <MatchNode.h>
#include <Domain.h>
#include <DistVector.h>

#include <cstdlib>
#include <cmath>

#define FORCE_TAG 1000
#define DISP_TAG 2000
#define INFO_TAG 3000
#define HEATPOWER_TAG 5000
#define TEMP_TAG 6000
#define STRUC_NUMPA_TAG 7500
#define FLUID_NUMPA_TAG 7550
#define STRUC_RELRES_TAG 7600
#define STRUC_CMD_TAG 8000
#define FLUID_CMD_TAG 9000
#define NEGO_NUM_TAG 10000
#define NEGO_BUF_TAG 10001

#define WET_SURF_TAG1 555
#define WET_SURF_TAG2 666
#define WET_SURF_TAG3 888
#define WET_SURF_TAG4 999
#define SUBCYCLING_TAG 777
#define SUGGEST_DT_TAG 444

#define CRACK_TAG1 22
#define CRACK_TAG2 33
#define CRACK_TAG3 44
#define CRACK_TAG4 55

//------------------------------------------------------------------------------

StructExc::StructExc(IoData& iod, MatchNodeSet** mns, int bs, Communicator* sc, Communicator* fluidCom, int nSub)
{

  bufsize = bs;
  com = fluidCom;
  strCom = sc;

  if (!strCom) {
    com->fprintf(stderr, "*** Error: structure communicator is null\n");
    exit(1);
  }

  oolscale = 1.0 / iod.ref.rv.tlength;
  ootscale = 1.0 / iod.ref.rv.time;
  oovscale = 1.0 / iod.ref.rv.tvelocity;
  ootempscale = 1.0 / iod.ref.rv.temperature;
  fscale = iod.ref.rv.tforce;
  pscale = iod.ref.rv.tpower;

  numLocSub = nSub;
  numStrCPU = strCom->remoteSize();
  numStrNodes = 0;

  sndParity = 0;
  recParity = 0;
  buffer = 0;

  matchNodes = mns;

}

//------------------------------------------------------------------------------

StructExc::~StructExc()
{

  if (numStrNodes) delete [] numStrNodes;
  if (buffer) delete [] buffer;

}

//------------------------------------------------------------------------------

void StructExc::updateNumStrNodes(int nn) 
{
  numStrNodes[0][0] = nn; 
  matchNodes[0]->updateNumStNodes(nn);
}

//------------------------------------------------------------------------------

void StructExc::negotiate()
{

  // compute the total number of matched nodes for this fluid CPU

  int (*ptr)[2] = new int[numLocSub][2];

  int iSub;
#pragma omp parallel for
  for (iSub=0; iSub<numLocSub; ++iSub)
    ptr[iSub][0] = matchNodes[iSub]->size();

  ptr[0][1] = 0;
  for (iSub=1; iSub<numLocSub; ++iSub)
    ptr[iSub][1] = ptr[iSub - 1][1] + ptr[iSub - 1][0];

  int numCPUMatchedNodes = ptr[numLocSub - 1][1] + ptr[numLocSub - 1][0];

  // gather the global matched point nodes for this fluid CPU

  int *ibuffer = 0;
  int (*cpuMatchNodes)[3] = 0;

  if (numCPUMatchedNodes > 0) {
    ibuffer = new int[numCPUMatchedNodes];
    cpuMatchNodes = new int[numCPUMatchedNodes][3];
    buffer = new double[bufsize * numCPUMatchedNodes];

#pragma omp parallel for
    for (iSub=0; iSub<numLocSub; ++iSub)
      matchNodes[iSub]->exportInfo(iSub, cpuMatchNodes + ptr[iSub][1]);
  }

  // send the matched node numbers of this fluid CPU to all the structure CPUs

  int i;
  for (i=0; i<numCPUMatchedNodes; ++i)
    ibuffer[i] = cpuMatchNodes[i][1];

  int iCpu;
  for (iCpu=0; iCpu<numStrCPU; ++iCpu) {
    strCom->sendTo(iCpu, NEGO_NUM_TAG, &numCPUMatchedNodes, 1);
    if (numCPUMatchedNodes > 0)
      strCom->sendTo(iCpu, NEGO_BUF_TAG, ibuffer, numCPUMatchedNodes);
  }
  strCom->waitForAllReq();

  // receive the list of matched nodes that each structure CPU contains

  if (numCPUMatchedNodes > 0) {

    numStrNodes = new int[numStrCPU][2];

    int pos = 0;

    for (iCpu=0; iCpu<numStrCPU; ++iCpu) {
      strCom->recFrom(iCpu, NEGO_NUM_TAG, &numStrNodes[iCpu][0], 1);

      if (numStrNodes[iCpu][0] > 0) {
	strCom->recFrom(iCpu, NEGO_BUF_TAG, ibuffer, numStrNodes[iCpu][0]);
	for (i=0; i<numStrNodes[iCpu][0]; ++i) {
	  int idx = ibuffer[i];
	  matchNodes[cpuMatchNodes[idx][2]]->setBufferPosition(cpuMatchNodes[idx][0], pos);
	  pos++;
	}
      }

    }

    if (pos != numCPUMatchedNodes) {
      fprintf(stderr, "*** Error: wrong number of matched nodes (%d instead of %d)\n",
	      pos, numCPUMatchedNodes);
      exit(1);
    }

    numStrNodes[0][1] = 0;
    for (iCpu=1; iCpu<numStrCPU; ++iCpu)
      numStrNodes[iCpu][1] = numStrNodes[iCpu-1][1] + numStrNodes[iCpu-1][0];

  }

  if (com->getMaxVerbose() >= 8)
    fprintf(stdout, "CPU %d has %d matched node%s\n", 
	    com->cpuNum(), numCPUMatchedNodes, numCPUMatchedNodes>1? "s":"");

  if (ptr) delete [] ptr;
  if (ibuffer) delete [] ibuffer;
  if (cpuMatchNodes) delete [] cpuMatchNodes;

}

//------------------------------------------------------------------------------

double StructExc::getInfo() 
{
  double info[5];

  if (strCom->cpuNum() == 0)
    strCom->recFrom(INFO_TAG, info, 5);

  com->broadcast(5, info);

  algNum = int(info[0]);
  dt = info[1] * ootscale;
  tmax = info[2] * ootscale;
  rstrt = int(info[3]);
  smode = int(info[4]);

  if (algNum == 6) tmax -= 0.5 * dt;
  if (algNum == 20) tmax -= 0.5 * dt;
  if (algNum == 21) tmax += 0.5 * dt;
  if (algNum == 22) tmax += 0.5 * dt;

  double mppFactor = 1.0;
  if (algNum == 8)
    mppFactor = 1.0/info[2];

  return mppFactor;

}

//------------------------------------------------------------------------------

void StructExc::getEmbeddedWetSurfaceInfo(int& elemType, bool& crack, int& nNodes, int& nElems)
{
  int info[4];
  if(strCom->cpuNum()==0) 
    strCom->recFrom(WET_SURF_TAG1, info, 4);
  com->broadcast(4, info);
  elemType = info[0]; //3~triangles, 4~quadrangles (triangles can be represented by degenerated quads)
  crack    = info[1] ? true : false;
  nNodes   = info[2];
  nElems   = info[3];
}

//------------------------------------------------------------------------------

void StructExc::getEmbeddedWetSurface(int nNodes, double *nodes, int nElems, int *elems, int eType)
{
  if(strCom->cpuNum()==0)
    strCom->recFrom(WET_SURF_TAG2, nodes, nNodes*3);
  com->broadcast(nNodes*3, nodes);

  if(strCom->cpuNum()==0)
    strCom->recFrom(WET_SURF_TAG3, elems, nElems*eType);
  com->broadcast(nElems*eType, elems);
}

//------------------------------------------------------------------------------

void StructExc::getInitialCrackingSetup(int& totalNodes, int& totalElems)
{
  int info[2];
  if(strCom->cpuNum()==0) 
    strCom->recFrom(WET_SURF_TAG4, info, 2);
  com->broadcast(2, info);
  totalNodes = info[0];
  totalElems = info[1];
}

//------------------------------------------------------------------------------

bool StructExc::getNewCrackingStats(int &numConnUpdate, int &numLSUpdate, int &newNodes)
{
  int size = 4;
  int nNew[size];
  if(strCom->cpuNum()==0) 
    strCom->recFrom(CRACK_TAG1, nNew, size);
  com->broadcast(size, nNew);
  numConnUpdate = nNew[1];
  numLSUpdate = nNew[2];
  newNodes = nNew[3];
  return nNew[0]>0;
}

//------------------------------------------------------------------------------

void StructExc::getInitialPhantomNodes(int newNodes, double(*xyz)[3], int nNodes)
{
  double coords[newNodes*3];
  if(strCom->cpuNum()==0)
    strCom->recFrom(CRACK_TAG4, coords, newNodes*3); //assume the correct ordering: nNodes, nNodes+1, ..., nNodes+newNodes-1
  com->broadcast(newNodes*3,coords);

  for(int i=0; i<newNodes; i++)
    for(int j=0; j<3; j++)
      xyz[nNodes+i][j] = coords[i*3+j];
}

//------------------------------------------------------------------------------

void StructExc::getNewCracking(int numConnUpdate, int numLSUpdate, int* phantoms, double* phi, int* phiIndex, int* new2old, int newNodes)
{
  int integer_pack_size = 5*numConnUpdate + numLSUpdate + 2*newNodes;
  int integer_pack[integer_pack_size]; //KW: This should be a short array since there will not be many new cracked elements in 
                                       //    one time-step. Therefore it should be OK to create and destroy it repeatedly.
  if(strCom->cpuNum()==0)
    strCom->recFrom(CRACK_TAG2, integer_pack, integer_pack_size);
  com->broadcast(integer_pack_size, integer_pack);
  for(int i=0; i<5*numConnUpdate; i++) {
    phantoms[i] = integer_pack[i];
    assert(integer_pack[i] >= 0);
  }
  for(int i=0; i<numLSUpdate; i++)
    phiIndex[i] = integer_pack[5*numConnUpdate+i];
  for(int i=0; i<newNodes; i++) {
    new2old[2*i] = integer_pack[5*numConnUpdate+numLSUpdate+2*i];
    new2old[2*i+1] = integer_pack[5*numConnUpdate+numLSUpdate+2*i+1];
  }

  if(strCom->cpuNum()==0)
    strCom->recFrom(CRACK_TAG3, phi, 4*numLSUpdate);
  com->broadcast(4*numLSUpdate, phi); 
}

//------------------------------------------------------------------------------
/*
     command communication between structure "0" and fluid "0"
     ! structure allways starts communication !

         structure sends            fluid has             fluid sends
     ------------------------------------------------------------------------
     small - talk

  1.     0 : continue               0 : continue          0: continue
  2.     0 : continue               1 : converged         0: continue
  3.     1 : converged              0 : continue          0: continue
  4.     1 : converged              1 : converged         1: stop analysis
                                                             next task
  5.     2 : converged              0/1                   1: stop analysis
                                                             next task
  5.    -1 : exit                   0/1                  -1: exit
  6.     0/1                       -1 : exit             -1: exit
     ------------------------------------------------------------------------
     command - talk

  7.   100 : analysis               0/1                 100: analysis
  8.   200 : senstivity             0/1                 200: senstivitiy
  9.   100/200                     -1 : exit            - 1: exit

*/

void StructExc::negotiateStopping(bool* lastIt)
{
  double xbuf;
  
  if (strCom->cpuNum() == 0) {
    strCom->recFrom(STRUC_CMD_TAG, &xbuf, 1);
    if (xbuf == 2.0)
      xbuf = 1.0;
    strCom->sendTo(0, FLUID_CMD_TAG, &xbuf, 1);
    strCom->waitForAllReq();
  }

  com->broadcast(1, &xbuf);  
  if (xbuf == 0.0)
    *lastIt = false;
  else
    *lastIt = true;

}

//------------------------------------------------------------------------------

void StructExc::sendNumParam(int numParam)
{
  double xbuf = double(numParam);

  if (strCom->cpuNum() == 0) {
    strCom->sendTo(0, FLUID_NUMPA_TAG, &xbuf, 1);
    strCom->waitForAllReq();
  }

  com->broadcast(1, &xbuf);
}

//------------------------------------------------------------------------------

void StructExc::getNumParam(int &numParam, int &actvar, double &steadyTol)
{
  double xbuf[3];

  if (strCom->cpuNum() == 0) {
    strCom->recFrom(STRUC_NUMPA_TAG, xbuf, 3);
    strCom->waitForAllReq();
  }

  com->broadcast(3, xbuf);
  numParam = int(xbuf[0]);
  actvar = int(xbuf[1]);
  steadyTol = xbuf[2];
}

//------------------------------------------------------------------------------

void StructExc::getRelResidual(double &relres)
{
  double xbuf;

  if (strCom->cpuNum() == 0) {
    strCom->recFrom(STRUC_RELRES_TAG, &xbuf, 1);
    strCom->waitForAllReq();
  }

  com->broadcast(1, &xbuf);
  relres = double(xbuf);
}

//------------------------------------------------------------------------------
/* 
   dX contains the displacement of the boundaries with respect to the CURRENT 
   configuration of the structure
*/

void StructExc::getDisplacement(DistSVec<double,3> &X0, DistSVec<double,3> &X, 
				DistSVec<double,3> &Xdot, DistSVec<double,3> &dX,
                                bool isEmbedded) 
{  
  double norms[2] = {0.0, 0.0};

  dX = 0.0;
  Xdot = 0.0;

  if (algNum == 4 || algNum == 5) recParity = 1 - recParity;

  if (numStrNodes) {

    for (int iCpu=0; iCpu<numStrCPU; ++iCpu) {
      if (numStrNodes[iCpu][0] > 0) {
        int size = bufsize * numStrNodes[iCpu][0];
        double *localBuffer = buffer + bufsize * numStrNodes[iCpu][1];
        strCom->recFrom(iCpu, DISP_TAG + recParity, localBuffer, size);
        com->printf(7, "[F] received displacement from structure\n");
      }
    }
    double (*disp)[2][3] = reinterpret_cast<double (*)[2][3]>(buffer);

#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      double (*x0)[3] = X0.subData(iSub);
      double (*x)[3] = X.subData(iSub);
      double (*xdot)[3] = Xdot.subData(iSub);
      double (*dx)[3] = dX.subData(iSub);

      double locNorms[2];
      matchNodes[iSub]->getDisplacement(algNum, dt, oolscale, oovscale, X.getMasterFlag(iSub), 
					disp, x0, x, xdot, dx, locNorms, isEmbedded);

#pragma omp critical
      norms[0] += locNorms[0];
#pragma omp critical
      norms[1] += locNorms[1];
    }

  }
  com->barrier(); //added for timing purposes (otherwise global comm timing can be very long
                  //if waiting for another cpu.
  com->globalSum(2, norms);

  com->printf(7, "Received total disp=%e and vel=%e from the structure\n", sqrt(norms[0]), sqrt(norms[1]));
}

//------------------------------------------------------------------------------
/* 
   dX contains the displacement of the boundaries with respect to the CURRENT 
   configuration of the structure
*/

void StructExc::getDisplacementSensitivity(DistSVec<double,3> &X, DistSVec<double,3> &dX, bool applyScale)
{  
  double norms[2] = {0.0, 0.0};

  dX = 0.0;

  if (algNum == 4 || algNum == 5) recParity = 1 - recParity;

  if (numStrNodes) {

    for (int iCpu=0; iCpu<numStrCPU; ++iCpu) {
      if (numStrNodes[iCpu][0] > 0) {
        int size = bufsize * numStrNodes[iCpu][0];
        double *localBuffer = buffer + bufsize * numStrNodes[iCpu][1];
        strCom->recFrom(iCpu, DISP_TAG + recParity, localBuffer, size);
        com->printf(7, "[F] received displacement sensitivity from structure\n"); 
      }
    }
    double (*disp)[2][3] = reinterpret_cast<double (*)[2][3]>(buffer);

#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      double (*dx)[3] = dX.subData(iSub);

      double locNorms[2];
      if(!applyScale) matchNodes[iSub]->getDisplacementSensitivity(X.getMasterFlag(iSub), 1.0, disp, dx, locNorms);  // for adjoint sensitivity analysis
      else matchNodes[iSub]->getDisplacementSensitivity(X.getMasterFlag(iSub), fscale, disp, dx, locNorms);  // for direct sensitivity analysis

#pragma omp critical
      norms[0] += locNorms[0];
#pragma omp critical
      norms[1] += locNorms[1];
    }

  }
  com->barrier(); //added for timing purposes (otherwise global comm timing can be very long
                  //if waiting for another cpu.
  com->globalSum(2, norms);

  com->printf(7, "Received total dispSensitivity=%e and velSensitivity=%e from the structure\n", sqrt(norms[0]), sqrt(norms[1]));
  //com->fprintf(stderr, "[StExc] Received total disp=%e and vel=%e from the structure\n", sqrt(norms[0]), sqrt(norms[1]));
}

//------------------------------------------------------------------------------

int StructExc::getSubcyclingInfo()
{
  double info;
  if (strCom->cpuNum() == 0)
    strCom->recFrom(SUBCYCLING_TAG, &info, 1);

  com->broadcast(1, &info);

  return (int)info;
}

//------------------------------------------------------------------------------

void StructExc::getTemperature(DistVec<double>& Temp)
{  

  double norm = 0.0;

  if (numStrNodes) {
    for (int iCpu=0; iCpu<numStrCPU; ++iCpu) {
      if (numStrNodes[iCpu][0] > 0) {
	int size = bufsize * numStrNodes[iCpu][0];
	double* localBuffer = buffer + bufsize * numStrNodes[iCpu][1];
	strCom->recFrom(iCpu, TEMP_TAG, localBuffer, size);
      }
    }

#pragma omp parallel for reduction(+: norm)
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      double* temp = Temp.subData(iSub);
      bool* flag = Temp.getMasterFlag(iSub);
      norm += matchNodes[iSub]->getTemperature(algNum, dt, ootempscale, flag, buffer, temp);
    }
  }

  com->barrier(); //added for timing purposes (otherwise global comm timing can be very long
                  //if waiting for another cpu.
  com->globalSum(1, &norm);

  com->printf(7, "Received temp=%e from the structure\n", norm);

}

//------------------------------------------------------------------------------
// note: the force vector is *** NOT *** assembled

void StructExc::sendForce(DistSVec<double,3> &F, bool applyScale)
{
  if (algNum == 4 || algNum == 5) sndParity = 1 - sndParity;

  double norm = 0.0;

  if (numStrNodes) {

    double (*forces)[3] = reinterpret_cast<double (*)[3]>(buffer);

//TODO BUGHUNT YUNGSOOS VERSION
//    if(F.info().masterFlag) {
//#pragma omp parallel for reduction (+: norm)
//      for (int iSub = 0; iSub < numLocSub; ++iSub) {
//        SVec<double,3> &f = F(iSub);
//        if(applyScale) {
//          norm += f*f * fscale*fscale;
//          matchNodes[iSub]->send(fscale, f.data(), forces); // for the direct sensitivity analysis
//        } else {
//          int locOffset = F.info().subOffset[iSub];
//          int locLen = F.info().subLen[iSub];
//          norm += f*f;
//          matchNodes[iSub]->sendWithMasterFlag(1, f.data(), forces, F.info().masterFlag, locOffset);  // for the adjoint sensitivity analysis
//
//        }
//      }
//    }

    //if(F.info().masterFlag) {
    if(true) {
#pragma omp parallel for reduction (+: norm)
      for (int iSub = 0; iSub < numLocSub; ++iSub) {
        SVec<double,3> &f = F(iSub);
        if(applyScale) {
          norm += f*f * fscale*fscale;
          matchNodes[iSub]->send(fscale, f.data(), forces); // for the direct sensitivity analysis
        }
        else if(F.info().masterFlag)
	{
          int locOffset = F.info().subOffset[iSub];
          int locLen = F.info().subLen[iSub];
          norm += f*f;
          matchNodes[iSub]->sendWithMasterFlag(1, f.data(), forces, F.info().masterFlag, locOffset);  // for the adjoint sensitivity analysis

        }
        else
        {
          std::cout<<"Cannot send when masterflag is NULL"<<std::endl;//TODO delete line
          exit(-1);
        }
      }
    }

    for (int iCpu=0; iCpu<numStrCPU; ++iCpu) {
      if (numStrNodes[iCpu][0] > 0) {
        int size = 3 * numStrNodes[iCpu][0];
        double *localBuffer = buffer + 3 * numStrNodes[iCpu][1];
        strCom->sendTo(iCpu, FORCE_TAG + sndParity, localBuffer, size);
      }
    }
    strCom->waitForAllReq();
  }

  com->barrier();
  com->globalSum(1, &norm);
  norm = sqrt(norm);

  com->printf(7, "Sent fluid force=%e to the structure\n", norm);
}

//------------------------------------------------------------------------------
// note: the heat power vector is *** NOT *** assembled

void StructExc::sendHeatPower(DistVec<double>& P) 
{

  double norm = 0.0;

  if (numStrNodes) {
    double (*buf)[1] = reinterpret_cast<double (*)[1]>(buffer);
#pragma omp parallel for reduction (+: norm)
    for (int iSub = 0; iSub < numLocSub; ++iSub) {
      Vec<double>& p = P(iSub);
      norm += p*p * pscale*pscale;
      matchNodes[iSub]->send(pscale, reinterpret_cast<double (*)[1]>(p.data()), buf);
    }

    for (int iCpu=0; iCpu<numStrCPU; ++iCpu) {
      if (numStrNodes[iCpu][0] > 0) {
	int size = numStrNodes[iCpu][0];
	double* localBuffer = buffer + numStrNodes[iCpu][1];
	strCom->sendTo(iCpu, HEATPOWER_TAG, localBuffer, size);
      }
    }
    strCom->waitForAllReq();
  }

  com->barrier();
  com->globalSum(1, &norm);

  com->printf(7, "Sent fluid heat power=%e to the structure\n", norm);

}

//------------------------------------------------------------------------------

void StructExc::sendFluidSuggestedTimestep(double dtf0)
{
  if(com->cpuNum()==0)
    for (int iCpu=0; iCpu<numStrCPU; ++iCpu) {
//      fprintf(stderr,"*** sending dtf = %e to DYNA!\n", dtf0);
      strCom->sendTo(iCpu, SUGGEST_DT_TAG, &dtf0, 1);
    }
}

//------------------------------------------------------------------------------

void StructExc::getMdFreq(int &nf, double *&f)
{

  double buffer[2000];
  if (strCom->cpuNum() == 0)
    strCom->recFrom(1200, buffer, 1000);

  com->broadcast(1000, buffer);

  nf = int(buffer[0]);
  f = new double[nf];

  for(int i=0; i<nf; ++i)
    f[i] = buffer[i+1];
}

//------------------------------------------------------------------------------

void StructExc::getMdStrDisp(int id, DistSVec<double,3> &X0,
                             DistSVec<double,3> &X, DistSVec<double,3> &dX)
{
  double norms[2] = {0.0, 0.0};

  DistSVec<double,3> Xdot( X.info() );

  dX = 0.0;
  Xdot = 0.0;

  if (numStrNodes) {

    for (int iCpu=0; iCpu<numStrCPU; ++iCpu) {
      if (numStrNodes[iCpu][0] > 0) {
        int size = bufsize * numStrNodes[iCpu][0];
        double *localBuffer = buffer + bufsize * numStrNodes[iCpu][1];
        strCom->recFrom(iCpu, 1201 + id, localBuffer, size);
      }
    }

    double (*disp)[2][3] = reinterpret_cast<double (*)[2][3]>(buffer);

#pragma omp parallel for
    for (int iSub = 0; iSub < numLocSub; ++iSub) {

      double (*x0)[3] = X0.subData(iSub);
      double (*x)[3] = X.subData(iSub);
      double (*xdot)[3] = Xdot.subData(iSub);
      double (*dx)[3] = dX.subData(iSub);

      double locNorms[2];

      matchNodes[iSub]->getDisplacement(algNum, dt, oolscale, oovscale, X.getMasterFlag(iSub),
                                        disp, x0, x, xdot, dx, locNorms, false);

#pragma omp critical
      norms[0] += locNorms[0];
#pragma omp critical
      norms[1] += locNorms[1];

    }

  }

  com->barrier();
  com->globalSum(2, norms);

  com->fprintf(stdout, "Received disp=%e and vel=%e from the structure \n", norms[0], norms[1]);
}
