#ifndef _COMMUNICATOR_H_
#define _COMMUNICATOR_H_

#include <cstdlib>
#include <cstdio>
#include <cstdarg>

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <Timer.h>
#include <ResizeArray.h>

class Connectivity;

extern void initCommunication(int &, char **&);
extern void closeCommunication();

//------------------------------------------------------------------------------

#ifdef USE_MPI
template <class T>
class CommTrace {

public:

  static MPI_Datatype MPIType;
  static int multiplicity;

};
#endif

//------------------------------------------------------------------------------

struct RecInfo {

  int cpu, len;

};

namespace Communication {
  template <typename Scalar>
    class Window;
}

//------------------------------------------------------------------------------
/** class for communicator with other processes */
class Communicator {

  int thisCPU;
  int numCPU;

public: //Needed by IntersectorPhysBAM
#ifdef USE_MPI
  MPI_Comm comm;
  int nPendReq;
  ResizeArray<MPI_Request> pendReq;
  ResizeArray<MPI_Status> reqStatus;
#endif

private:
  Timer *timer;
  int maxverbose;

public:

  Communicator();
#ifdef USE_MPI
  Communicator(MPI_Comm);
#endif
  ~Communicator(){};

  void split(int, int, Communicator**);
  Communicator *merge(bool high); //<! returns the intra-communicator for the two groups of an inter-communicator
  int remoteSize();
  int barrier();
  int sync() { return(barrier()); } ;
  void waitForAllReq();
  void printf(int, const char *, ...);
  void fprintf(FILE *, const char *, ...);

  void system(const char* message);

  template<class Scalar>
  void sendTo(int, int, Scalar *, int);

  template<class Scalar>
  RecInfo recFrom(int, Scalar *, int);

  template<class Scalar>
  RecInfo recFrom(int, int, Scalar *, int);

  template<class Scalar>
  void exchange(int, int, int *, int *, Scalar **, int *, Scalar **);

  template<class Scalar>
  void broadcast(int, Scalar *, int = 0);

#ifdef USE_MPI
  template<class Scalar>
  void globalOp(int, Scalar *, MPI_Op);
#endif

  template<class Scalar>
  void globalMin(int, Scalar *);

  template<class Scalar>
  void globalMax(int, Scalar *);

  template<class Scalar>
  void globalSum(int, Scalar *);

#ifdef USE_MPI
  template<class Scalar>
  void globalOpRoot(int, Scalar *, MPI_Op); // (MPI_Reduce vs. MPI_Reduceall)
#endif

  template<class Scalar>
  void globalMinRoot(int, Scalar *); // only cpu 0 gets result (MPI_Reduce vs. MPI_Reduceall)

  template<class Scalar>
  void globalMaxRoot(int, Scalar *); // only cpu 0 gets result (MPI_Reduce vs. MPI_Reduceall)

  template<class Scalar>
  void globalSumRoot(int, Scalar *); // only cpu 0 gets result (MPI_Reduce vs. MPI_Reduceall)

  int size() const { return numCPU; }
  int cpuNum() const { return thisCPU; }
  void setTimer(Timer *t) { timer = t; }
  void setMaxVerbose(int v) { maxverbose = v; }
  int getMaxVerbose() { return maxverbose; }
    MPI_Comm getMPIComm() {return comm;} //<! Lei Lei, 24 March 2016, needed for ALS

  template <typename Scalar>
    friend class Communication::Window;
};

namespace Communication {

  template<typename Scalar>
  class Window {
#ifdef USE_MPI
    MPI_Win win;
#else
    Scalar *data;
#endif
    Communicator &com;
    int length;
  public:
    static const int Add=0, Min=1, Max=2;
    Window(Communicator &c, int size, Scalar *s);
    ~Window();
    //void get(int locOff, int size, int prNum, int remOff);
//    void put(int locOff, int size, int prNum, int remOff);
    void put(Scalar *s, int locOff, int size, int prNum, int remOff);
    void accumulate(Scalar *s, int locOff, int size, int prNum, int remOff, int op);
    void fence(bool startOrEnd);
    int size() {return length;}
  };

}
/** allocate memory that can be used for one-sided communication */
 void* operator new(size_t, Communicator &c);
 void operator delete(void *p, Communicator &c);

 void *operator new[](size_t size, Communicator &c);

 void operator delete[](void *p, Communicator &c);

//------------------------------------------------------------------------------

class SubDTopo {

  struct CPair {

    int from, to, cpuID;

    CPair() {}
    CPair(int f, int t, int c) { from =f; to = t; cpuID =c; }

    // the following operator is required by the STL sort algorithm
    bool operator < (const CPair &x) const
    {
      // we want to order first by cpuID (with local coms first)
      // then by origin and then by destination
      return cpuID < x.cpuID || (cpuID  == x.cpuID &&
				 (from < x.from || (from == x.from && to < x.to)));
    }
    bool operator == (const CPair &x) const
    {
      return cpuID == x.cpuID && from == x.from && to == x.to;
    }

  };

  int numCPU;// number of CPUs
  int cpuNum;// this CPU's number
  int localNumSub;
  int *glSubToLocal;
  int *localSubToGl;
  int *subNumNeighb;
  int *glSubToCPU;

  int *neighbCPU;
  Connectivity *cpuSC;
  Connectivity *cpuRC;
  int numPairs;
  CPair *allPairs;

  void makeLocalPairIndex(Connectivity *subToSub);
  void makeCrossConnect(); // builds neighbCPU cpuSC and cpuRC

public:
  // The constructor only needs the connectivity of the subdomains of
  // this CPU to be correct, the connectivity of other subdomains can be
  // (and is in slave domain) ommited.
  SubDTopo (int CPU, Connectivity *subToSub, Connectivity *CPUToSub);
  ~SubDTopo();

  int getChannelID(int glFrom, int glTo);
  int numChannels() { return numPairs; }
  int numNeighbCPUs();
  int crossNeighb(int iCpu) { return neighbCPU[iCpu]; }
  int *crossNeighb() { return neighbCPU; }
  Connectivity *cpuSndChannels() { return cpuSC; }
  Connectivity *cpuRcvChannels() { return cpuRC; }
  int reverseChannel(int);
  int sourceIsLocal(int channel) { return glSubToCPU[allPairs[channel].from] == cpuNum; }
  bool isSubLocal(int sub) { return glSubToCPU[sub] == cpuNum; }
  int locSubNum(int sub) { return glSubToLocal[sub]; }
};

//------------------------------------------------------------------------------

template <class T>
struct SubRecInfo {

    T *data;
    int len;
    int leadDim;
    int nvec;

};

//------------------------------------------------------------------------------
/*
   CommPattern represent a communication pattern.
     Communication is based on the model that a message from one
     subdomain to another subdomain is made of a number of vectors.
     Such vectors are stored in a matrix form. That matrix may have
     larger columns than the actual message length. This allows a subdomain
     to send different subparts of a given matrix to different subdomains.

   A communication patterns contains the following information:
     - Length of a vector from one subdomain to its neighbors
     - Number of vectors in each message
     - Leading dimension of the matrix containing the vectors
*/
template <class T>
class CommPattern {

public:
  enum Mode { Share, CopyOnSend };
  enum Symmetry { Sym, NonSym };

  // Buffers for local copy-communication
protected:
  Mode mode;
  Symmetry sym;
  T *localDBuffer;

  // for non-local communication we need a few more info
  int numCPUs;
  int numNeighbCPUs;
  int *neighbCPUs;
  Connectivity *sndConnect;
  Connectivity *rcvConnect;
  // in cpuForChannel, we use the convention that this CPU is marked as -1
  int *cpuForChannel;
  int *isSend;  // wether this channel is a send or receive channel
  int *reverseChannel; // corresponding reverse channel
  // Buffers for cross-memory communication on a cpu per cpu basis
  int *crossSendLen, *crossRcvLen;
  T **crossSendBuffer;
  T **crossRcvBuffer;

  int numChannels;
  SubRecInfo<T> *sRecInfo;

  Communicator *communicator;

public:

  CommPattern(SubDTopo *, Communicator *, Mode = Share, Symmetry = Sym);
  ~CommPattern()
    {
      delete[] cpuForChannel;
      delete[] reverseChannel;
      delete[] isSend;
      delete[] sRecInfo;
      if(crossSendBuffer) {
        if(crossSendBuffer[0]) delete[] crossSendBuffer[0];
        delete[] crossSendBuffer;
      }
      if(crossRcvBuffer) delete[] crossRcvBuffer;
      if(crossSendLen) delete[] crossSendLen;
      if(crossRcvLen) delete[] crossRcvLen;
      if(localDBuffer) delete[] localDBuffer;
    }

#ifdef MEM_TMPL_FUNC
  template <class TB>
  CommPattern(CommPattern<TB> &, Communicator *, Mode = Share, Symmetry = Sym);
#endif

  void setLen(int ID, int size, int lDim=-1, int nv = 1);
  int getLen(int ID) const { return sRecInfo[ID].len; }
  void finalize(); // complete the internal setup
  SubRecInfo<T> recData(int channel);
  void sendData(int channel, T *);
  SubRecInfo<T> getSendBuffer(int channel);
  void exchange();
  int numCh(){ return numChannels; }

};

//------------------------------------------------------------------------------

#include <alloca.h>
#ifdef TEMPLATE_FIX
#include "Communicator.C"
#endif

#endif
