#ifndef _DIST_INFO_H_
#define _DIST_INFO_H_

#include <Communicator.h>

//------------------------------------------------------------------------------

struct DistInfo {

  int numLocThreads;
  int numLocSub;
  int totLen;
  int *subLen;
  int *subOffset;
  int *subLenReg;
  int *subOffsetReg;
  bool *masterFlag;
  double *invNdWeight; // inverse of the number of subs touching a node

  int numGlobSub;
  int *locSubToGlobSub;

  Communicator *com;

  DistInfo(int _numLocThreads, int _numLocSub, int _numGlobSub, 
	   int *_locSubToGlobSub, Communicator *_com)
  { 
    numLocThreads   = _numLocThreads;
    numLocSub       = _numLocSub;
    numGlobSub      = _numGlobSub;
    locSubToGlobSub = _locSubToGlobSub;
    com             = _com;
    subLen          = new int[numLocSub]; 
    subOffset       = new int[numLocSub];

    subLenReg      = 0;
    subOffsetReg   = 0;

    masterFlag     = 0;
    invNdWeight    = 0;
  }

  void setLen(int sub, int len) { subLen[sub] = len; }

  int subSize(int sub) const { return (subLen) ? subLen[sub] : 0; } //HB

  bool* getMasterFlag(int iSub) const
  {
    return (masterFlag) ? masterFlag+subOffset[iSub] : 0; 
  }

  double* getInvWeight(int iSub) const
  {
    return (invNdWeight) ? invNdWeight+subOffset[iSub] : 0;
  }

  void finalize(bool makeFlag) 
  {
    subOffset[0] = 0;

    for (int iSub = 1; iSub < numLocSub; ++iSub)
      subOffset[iSub] = subOffset[iSub-1] + subLen[iSub-1];

    totLen = subOffset[numLocSub-1] + subLen[numLocSub-1];

    if (makeFlag) {
       masterFlag = new bool[totLen];
       invNdWeight= new double[totLen];
    } else {
       masterFlag = 0;
       invNdWeight= 0;
    }

    numLocThreads = numLocSub;
    subLenReg     = subLen;
    subOffsetReg  = subOffset;    

    /* 
    int load = totLen / numLocThreads;
    int remainder = totLen % numLocThreads;

    for (int iSub=0; iSub<numLocThreads; ++iSub)
      subLenReg[iSub] = load;

    for (int iSub=0; iSub<remainder; ++iSub)
      subLenReg[iSub] += 1;

    subOffsetReg[0] = 0;

    for (int iSub = 1; iSub < numLocThreads; ++iSub)
      subOffsetReg[iSub] = subOffsetReg[iSub-1] + subLenReg[iSub-1];
    */
  }

  void print(char* mssg=0) const 
  {
    com->sync();
    if(mssg) com->fprintf(stderr," DistInfo::print of %s\n",mssg);
    for(int iCPU=0; iCPU<com->size(); iCPU++) { 
      if(iCPU==com->cpuNum()) {
        fprintf(stderr," DistInfo::print on CPU %d:\n",com->cpuNum()); fflush(stderr);
        for(int iSub=0; iSub<numLocSub; ++iSub) {
          fprintf(stderr," * subd %3d, subOffset[%3d] = %6d, subLen[%3d] = %6d\n",
                locSubToGlobSub[iSub],iSub,subOffset[iSub],iSub,subLen[iSub]); fflush(stderr);
        }
      }
      com->sync();
    }
  }

  ~DistInfo() 
  {
    delete [] invNdWeight;
    delete [] masterFlag;
    delete [] subOffset;
    delete [] subLen; 
  }

  // Copy Constructor and Assignement Operator are declared Private
  // Compiler will return an error if used.
private:
  DistInfo(const DistInfo &);
  DistInfo& operator=(const DistInfo &);

};

//------------------------------------------------------------------------------
#endif
