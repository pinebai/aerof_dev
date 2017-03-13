#include <Communicator.h>
#include <Connectivity.h>

#include <algorithm>

template <class T>
CommPattern<T>::CommPattern(SubDTopo *topo, Communicator *_com,
			    Mode _mode, Symmetry _sym)
{
 int i, j;
 communicator = _com;
 mode = _mode;
 sym = _sym;
 numChannels = topo->numChannels();
 numCPUs = communicator->size();
 sRecInfo = new SubRecInfo<T>[numChannels];
 isSend = new int[numChannels];
 reverseChannel = new int[numChannels];
 cpuForChannel = new int[numChannels]; // -1 if local, cpu number otherwise
 for(i = 0; i < numChannels; ++i) {
    sRecInfo[i].len = 1;
    sRecInfo[i].leadDim = 1;
    sRecInfo[i].nvec = 1;
    reverseChannel[i] = topo->reverseChannel(i);
    isSend[i] = topo->sourceIsLocal(i);
    cpuForChannel[i] = -1;
 }
 if(numCPUs > 1) {
   sndConnect = topo->cpuSndChannels();
   rcvConnect = topo->cpuRcvChannels();
   numNeighbCPUs = sndConnect->csize();
   neighbCPUs = topo->crossNeighb();
   for(i = 0; i < numNeighbCPUs; ++i) {
     for(j = 0; j < sndConnect->num(i); ++j)
       cpuForChannel[(*sndConnect)[i][j] ] = neighbCPUs[i];
     for(j = 0; j < rcvConnect->num(i); ++j)
       cpuForChannel[(*rcvConnect)[i][j] ] = neighbCPUs[i];
   }
 }
 crossSendLen = 0;
 crossRcvLen = 0;
 crossSendBuffer = 0;
 crossRcvBuffer = 0;
 localDBuffer = 0;
}

template <class T>
void
CommPattern<T>::setLen(int channel, int len, int ldim, int nvec)
{
 if(ldim < 0) ldim = len;
 sRecInfo[channel].len  = len;
 sRecInfo[channel].leadDim = ldim;
 sRecInfo[channel].nvec = nvec;
}

template <class T>
void
CommPattern<T>::finalize()
{
 // If the number of cpus is > 1 then we need to get the length of messages
 // for the other CPUs.
 if(numCPUs > 1) {
   if(sym == Sym) {
      int i;
      for(i = 0; i < numChannels; ++i) {
        if(cpuForChannel[i] >= 0 && isSend[i]) {
          sRecInfo[reverseChannel[i]].len  = sRecInfo[i].len;
          sRecInfo[reverseChannel[i]].leadDim = sRecInfo[i].leadDim;
          sRecInfo[reverseChannel[i]].nvec = sRecInfo[i].nvec;
        }
      }
   } else {
      // send and receive the message length and number of vectors (leadDim is
      // implicitely set to the message length)
      int nSndCh = sndConnect->numConnect();
      int nRcvCh = rcvConnect->numConnect();
      //fprintf(stderr, "Number of neighb: %d with %d and %d\n", numNeighbCPUs,nSndCh,nRcvCh);
      int *sndMsg = new int[2*nSndCh];
      int *rcvMsg = new int[2*nRcvCh];
      int **sndPtr = new int*[numNeighbCPUs];
      int **rcvPtr = new int*[numNeighbCPUs];
      int *sndLen = new int[numNeighbCPUs];
      int *rcvLen = new int[numNeighbCPUs];
      int offset = 0;
      int iCPU, i;
      int totLen = 0;
      for(iCPU = 0; iCPU < numNeighbCPUs; ++iCPU) {
        sndPtr[iCPU] = sndMsg + offset;
        rcvPtr[iCPU] = rcvMsg + offset;
        sndLen[iCPU] = 2* sndConnect->num(iCPU);
        rcvLen[iCPU] = 2* sndConnect->num(iCPU);
        for(i = 0; i < sndConnect->num(iCPU); ++i) {
          sndPtr[iCPU][2*i]   = sRecInfo[ (*sndConnect)[iCPU][i] ].len;
          sndPtr[iCPU][2*i+1] = sRecInfo[ (*sndConnect)[iCPU][i] ].nvec;
          totLen += sndMsg[offset]*sndMsg[offset+1];
          offset += 2;
        }
      }
      // Do the actual exchange MLX
      //fprintf(stderr, "Tot len is %d  for %d\n", totLen, communicator->cpuNum());
      communicator->exchange(101, numNeighbCPUs, neighbCPUs, sndLen, sndPtr,
			     rcvLen, rcvPtr);
      offset = 0;
      totLen = 0;
      for(iCPU = 0; iCPU < numNeighbCPUs; ++iCPU)
        for(i = 0; i < sndConnect->num(iCPU); ++i) {
          sRecInfo[ (*rcvConnect)[iCPU][i] ].len = rcvMsg[offset];
          sRecInfo[ (*rcvConnect)[iCPU][i] ].leadDim = rcvMsg[offset];
          sRecInfo[ (*rcvConnect)[iCPU][i] ].nvec = rcvMsg[offset+1];
          totLen += rcvMsg[offset]*rcvMsg[offset+1];
          offset += 2;
        }
      //fprintf(stderr, "Tot=len is %d  for %d\n", totLen, communicator->cpuNum());
      delete [] sndMsg;
      delete [] rcvMsg;
      delete [] sndPtr;
      delete [] rcvPtr;
      delete [] sndLen;
      delete [] rcvLen;
   }
 }
 // in the case of single CPU communication, Share mode does not require any
 // work, independently of whether we are symmetric or not
 if(mode == Share &&numCPUs == 1) return;

 if(numCPUs != 1 && numNeighbCPUs > 0){
    // When there are several other CPUs, we need to allocate
    // the cross CPU memory buffers
    int i, j, totLen=0;
    crossSendBuffer = new T *[numNeighbCPUs];
    crossRcvBuffer  = new T *[numNeighbCPUs];
    crossSendLen = new int[numNeighbCPUs];
    crossRcvLen  = new int[numNeighbCPUs];
    for(i = 0; i < numNeighbCPUs; ++i)
      crossSendLen[i] = crossRcvLen[i] = 0;
    for(i = 0; i < numNeighbCPUs; ++i) {
      for(j = 0; j < sndConnect->num(i); ++j) {
    	int channel = (*sndConnect)[i][j];
    	int msgLen = sRecInfo[channel].len*sRecInfo[channel].nvec;
    	crossSendLen[i] += msgLen;
    	totLen += msgLen;
       }
       for(j = 0; j < rcvConnect->num(i); ++j) {
    	int channel = (*rcvConnect)[i][j];
    	int msgLen = sRecInfo[channel].len*sRecInfo[channel].nvec;
    	crossRcvLen[i] += msgLen;
    	totLen += msgLen;
       }
    }
    T *crBuff = new T[totLen];
    for(i = 0; i < numNeighbCPUs; ++i) {
       crossSendBuffer[i] = crBuff;
       for(j = 0; j < sndConnect->num(i); ++j) {
    	 int channel = (*sndConnect)[i][j];
    	 int msgLen = sRecInfo[channel].len*sRecInfo[channel].nvec;
         sRecInfo[channel].data = crBuff;
         crBuff += msgLen;
       }
       crossRcvBuffer[i] = crBuff;
       for(j = 0; j < rcvConnect->num(i); ++j) {
    	 int channel = (*rcvConnect)[i][j];
    	 int msgLen = sRecInfo[channel].len*sRecInfo[channel].nvec;
         sRecInfo[channel].data = crBuff;
         crBuff += msgLen;
       }
    }
 }
 if(mode == Share) return;
 // In Copy On Send mode, we need to allocate a T array of appropriate length
 // Compute the length:
 int i, len=0;
 for(i = 0; i < numChannels; ++i)
   if(numCPUs == 1 || cpuForChannel[i] < 0)
        len += sRecInfo[i].len*sRecInfo[i].nvec;
 localDBuffer = new T[len];

 // Now update the pointers.
 T *cBuf = localDBuffer;
 for(i = 0; i < numChannels; ++i)
   if(numCPUs == 1 || cpuForChannel[i] < 0)
     {
      sRecInfo[i].data = cBuf;
      cBuf += sRecInfo[i].len*sRecInfo[i].nvec;
     }
}

template <class T>
void
CommPattern<T>::exchange()
{

  if (numCPUs == 1) return;
  // Do the actual exchange MLX

  communicator->exchange(101, numNeighbCPUs, neighbCPUs, crossSendLen,
			 crossSendBuffer, crossRcvLen, crossRcvBuffer);

}

template <class T>
SubRecInfo<T>
CommPattern<T>::recData(int channel)
{
 if(mode == Share) return sRecInfo[channel];
 SubRecInfo<T> ret = sRecInfo[channel];
 // in Copy On Send mode, for sending we have packed the vectors and therefore
 // the leading dimension of received data is the length of a vector
 ret.leadDim = ret.len;
 return ret;
}

template <class T>
void
CommPattern<T>::sendData(int channel, T *data)
{
 if(mode == Share && (numCPUs == 1 || cpuForChannel[channel] < 0)) {
   sRecInfo[channel].data = data;
   return;
 }
 // For Copy On Send, or non local communication,
 // copy the data into the right buffer
 int i, iVec;
 SubRecInfo<T> chObj = sRecInfo[channel];

 for(iVec = 0; iVec < chObj.nvec; ++iVec)
   for(i = 0; i < chObj.len; ++ i)
      chObj.data[iVec*chObj.len + i] = data[iVec*chObj.leadDim + i];
}

template <class T>
SubRecInfo<T>
CommPattern<T>::getSendBuffer(int channel)
{
 SubRecInfo<T> chObj = sRecInfo[channel];
 if(mode == Share) // Check if we have a send buffer. If not, return NULL
   chObj.data = NULL;
 return chObj;
}

//------------------------------------------------------------------------------

template <class Scalar>
void Communicator::sendTo(int cpu, int tag, Scalar *buffer, int len)
{

#ifdef USE_MPI
  double t0;
  if (timer) t0 = timer->getTime();

  int thisReq = nPendReq++;
  MPI_Request *req = pendReq+thisReq;

  MPI_Isend(buffer, CommTrace<Scalar>::multiplicity*len, CommTrace<Scalar>::MPIType, cpu, tag, comm, req);

  if (timer) timer->addInterComTime(t0);
#endif

}

//------------------------------------------------------------------------------

template<class Scalar>
RecInfo Communicator::recFrom(int tag, Scalar *buffer, int len)
{

  RecInfo rInfo;

#ifdef USE_MPI
  double t0;
  if (timer) t0 = timer->getTime();

  MPI_Status status;
  MPI_Recv(buffer, CommTrace<Scalar>::multiplicity*len, CommTrace<Scalar>::MPIType, MPI_ANY_SOURCE, tag, comm, &status);
  MPI_Get_count(&status, CommTrace<Scalar>::MPIType, &rInfo.len);
  rInfo.len /= CommTrace<Scalar>::multiplicity;
  rInfo.cpu = status.MPI_SOURCE;
  if (timer) timer->addInterComTime(t0);
#endif

  return rInfo;

}

//------------------------------------------------------------------------------

template<class Scalar>
RecInfo Communicator::recFrom(int cpu, int tag, Scalar *buffer, int len)
{

  RecInfo rInfo;

#ifdef USE_MPI
  double t0;
  if (timer) t0 = timer->getTime();

  MPI_Status status;
  MPI_Recv(buffer, CommTrace<Scalar>::multiplicity*len, CommTrace<Scalar>::MPIType, cpu, tag, comm, &status);
  MPI_Get_count(&status, CommTrace<Scalar>::MPIType, &rInfo.len);
  rInfo.len /= CommTrace<Scalar>::multiplicity;
  rInfo.cpu = status.MPI_SOURCE;

  if (timer) timer->addInterComTime(t0);
#endif

  return rInfo;

}

//------------------------------------------------------------------------------

template<class Scalar>
void Communicator::exchange(int tag, int numNeighb, int *cpus, int *sndLen,
			    Scalar **sndData, int *rcvLen, Scalar **rcvData)
{

  //fprintf(stderr, " ... Exchanging w/ %d Neigh CPUS, %d crossSendLen %d, crossRcvLen\n", numNeighb, sndLen, rcvLen);
#ifdef USE_MPI

  double t0;
  if (timer) t0 = timer->getTime();

  MPI_Request *sndId = reinterpret_cast<MPI_Request *>(alloca(sizeof(MPI_Request)*numNeighb));
  MPI_Request *rcvId = reinterpret_cast<MPI_Request *>(alloca(sizeof(MPI_Request)*numNeighb));

  int iCpu;
  int rcvReq = 0;
  for (iCpu = 0; iCpu < numNeighb; ++iCpu) {
    int len = rcvLen[iCpu];
    if (len == 0) continue;
    MPI_Irecv(rcvData[iCpu], CommTrace<Scalar>::multiplicity*len, CommTrace<Scalar>::MPIType,
	      cpus[iCpu], tag, comm, rcvId+rcvReq);
    rcvReq += 1;
  }

  int sendReq = 0;
  for (iCpu = 0; iCpu < numNeighb; ++iCpu) {
    int len = sndLen[iCpu];
    if (len == 0) continue;
    MPI_Isend(sndData[iCpu], CommTrace<Scalar>::multiplicity*len, CommTrace<Scalar>::MPIType,
	      cpus[iCpu], tag, comm, sndId+sendReq);
    sendReq += 1;
  }

  MPI_Status *status = reinterpret_cast<MPI_Status *>(alloca(sizeof(MPI_Status)*numNeighb));
  MPI_Waitall(rcvReq, rcvId, status);
  MPI_Waitall(sendReq, sndId, status);

  if (timer) timer->addLocalComTime(t0);
#endif

}

//------------------------------------------------------------------------------

template<class Scalar>
void Communicator::broadcast(int n, Scalar *b, int root)
{

#ifdef USE_MPI
  MPI_Bcast(b, CommTrace<Scalar>::multiplicity*n, CommTrace<Scalar>::MPIType, root, comm);
#endif

}

//------------------------------------------------------------------------------

#ifdef USE_MPI
template<class Scalar>
void Communicator::globalOp(int len, Scalar *x, MPI_Op op)
{

  if (numCPU < 2) return;

  double t0;
  if (timer) t0 = timer->getTime();

  int  segSize = (len > 65536) ? 65536 : len;

  Scalar *work;
  bool delete_work = false;

  if (segSize > 5000) {
    delete_work =  true;
    work = new Scalar[segSize];
  } else
    work = reinterpret_cast<Scalar *>(alloca(sizeof(Scalar)*segSize));

  for (int offset = 0; offset < len; offset += segSize) {
    if (offset+segSize > len)
      segSize = len-offset;
    MPI_Allreduce(x+offset, work, CommTrace<Scalar>::multiplicity*segSize, CommTrace<Scalar>::MPIType, op, comm);
    // Should there be a barrier here? MPI_Allreduce does not guarantee synchronization!
    for (int j = 0; j < segSize; ++j)
      x[offset+j] = work[j];
  }

  if (delete_work) delete [] work;

  if (timer) timer->addGlobalComTime(t0);

}
#endif

//------------------------------------------------------------------------------

template<class Scalar>
void Communicator::globalMin(int len, Scalar *x)
{

#ifdef USE_MPI
  globalOp(len, x, MPI_MIN);
#endif

}

//------------------------------------------------------------------------------

template<class Scalar>
void Communicator::globalMax(int len, Scalar *x)
{

#ifdef USE_MPI
  globalOp(len, x, MPI_MAX);
#endif

}

//------------------------------------------------------------------------------

template<class Scalar>
void Communicator::globalSum(int len, Scalar *x)
{

#ifdef USE_MPI
  globalOp(len, x, MPI_SUM);
#endif

}

//------------------------------------------------------------------------------

#ifdef USE_MPI
template<class Scalar>
void Communicator::globalOpRoot(int len, Scalar *x, MPI_Op op)
{

  if (numCPU < 2) return;

  double t0;
  if (timer) t0 = timer->getTime();

  int root = 0;
  int segSize = (len > 65536) ? 65536 : len;

  Scalar *work;
  bool delete_work = false;

  if (segSize > 5000) {
    delete_work =  true;
    work = new Scalar[segSize];
  } else
    work = reinterpret_cast<Scalar *>(alloca(sizeof(Scalar)*segSize));

  for (int offset = 0; offset < len; offset += segSize) {
    if (offset+segSize > len)
      segSize = len-offset;
    MPI_Reduce(x+offset, work, CommTrace<Scalar>::multiplicity*segSize, CommTrace<Scalar>::MPIType, op, root, comm);
    // Should there be a barrier here?
    for (int j = 0; j < segSize; ++j)
      x[offset+j] = work[j];
  }

  if (delete_work) delete [] work;

  if (timer) timer->addGlobalComTime(t0);

}
#endif


//------------------------------------------------------------------------------

template<class Scalar>
void Communicator::globalMinRoot(int len, Scalar *x)
{

#ifdef USE_MPI
  globalOpRoot(len, x, MPI_MIN);
#endif

}

//------------------------------------------------------------------------------

template<class Scalar>
void Communicator::globalMaxRoot(int len, Scalar *x)
{

#ifdef USE_MPI
  globalOpRoot(len, x, MPI_MAX);
#endif

}

//------------------------------------------------------------------------------

template<class Scalar>
void Communicator::globalSumRoot(int len, Scalar *x)
{

#ifdef USE_MPI
  globalOpRoot(len, x, MPI_SUM);
#endif

}

//------------------------------------------------------------------------------

#ifdef MEM_TMPL_FUNC
template <class T>
template <class TB>
CommPattern<T>::CommPattern(CommPattern<TB> &pb, Communicator *_com,
			    Mode _mode, Symmetry _sym)
{
  communicator = _com;
 mode = _mode;
 sym = _sym;
 numChannels = pb.numCh();
 sRecInfo = new SubRecInfo<T>[numChannels];
 for(int i = 0; i < numChannels; ++i) {
    SubRecInfo<TB> rInfo = pb.getSendBuffer(i);
    sRecInfo[i].len = rInfo.len;
    sRecInfo[i].leadDim = rInfo.leadDim;
    sRecInfo[i].nvec = rInfo.nvec;
 }
 finalize();
}
#endif

namespace Communication {

  template <typename Scalar>
  Window<Scalar>::Window(Communicator &c, int size, Scalar *s) : com(c) {
#ifdef USE_MPI
    MPI_Win_create(s, size, sizeof(Scalar), MPI_INFO_NULL,
        com.comm, &win);
#else
    data = s;
#endif
    length = size;
  }

  template <typename Scalar>
    Window<Scalar>::~Window() {
  #ifdef USE_MPI
      MPI_Win_free(&win);
  #endif
  }

  template <typename Scalar>
  void Window<Scalar>::put(Scalar *a, int locOff, int size, int prNum, int remOff) {
#ifdef USE_MPI
    double t0;
    if (com.timer) t0 = com.timer->getTime();
 
   static const MPI_Op mpiOp[] = { MPI_SUM, MPI_MIN, MPI_MAX };
    MPI_Put(a+locOff, size*CommTrace<Scalar>::multiplicity, CommTrace<Scalar>::MPIType,
        prNum,
        remOff*CommTrace<Scalar>::multiplicity, size*CommTrace<Scalar>::multiplicity,
        CommTrace<Scalar>::MPIType,
        win);

    if (com.timer) com.timer->addRMAComTime(t0);
#else
    const double * bufferBegin = a + locOff;
    std::copy(bufferBegin, bufferBegin + size, data + remOff);
#endif 
  }

  template <typename Scalar>
  void Window<Scalar>::accumulate(Scalar *a, int locOff, int size, int prNum, int remOff, int op) {
#ifdef USE_MPI
    double t0;
    if (com.timer) t0 = com.timer->getTime();

    static const MPI_Op mpiOp[] = { MPI_SUM, MPI_MIN, MPI_MAX };
    MPI_Accumulate(a+locOff, size*CommTrace<Scalar>::multiplicity, CommTrace<Scalar>::MPIType,
        prNum,
        remOff*CommTrace<Scalar>::multiplicity, size*CommTrace<Scalar>::multiplicity,
        CommTrace<Scalar>::MPIType, mpiOp[op],
        win);

    if (com.timer) com.timer->addRMAComTime(t0);
#else
    switch(op)
      {
      case Add:
	{
	  for(int i=0;i<size;++i)
	    {
	      data[remOff+i] += a[locOff+i];
	    }
	  break;
	}
      case Min:
	{
	  for(int i=0;i<size;++i)
	    {
	      data[remOff+i] = std::min(data[remOff+i],a[locOff+i]);
	    }
	  break;
	}
      case Max:
	{
	  for(int i=0;i<size;++i)
	    {
	      data[remOff+i] = std::max(data[remOff+i],a[locOff+i]);
	    }
	  break;
	}
      default:
	{
	fprintf(stderr,"WOW! Stop\n");exit(-1);
	}
      }
#endif
  }

  template <typename Scalar>
    void Window<Scalar>::fence(bool isBeginning) {
#ifdef USE_MPI
      double t0;
      if (com.timer) t0 = com.timer->getTime();

      if(isBeginning)
        MPI_Win_fence((MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE), win);
      else
        MPI_Win_fence(MPI_MODE_NOSUCCEED, win);

      if (com.timer) com.timer->addRMAComTime(t0);
#endif
  }
}
