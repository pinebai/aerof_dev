#include <alloca.h>
#include <cstdio>
#include <cstdlib>

#include <Communicator.h>
#include <complex>
#include <new>

#ifndef MPI_INTEGER
#define MPI_INTEGER MPI_INT
#endif

#include <complex>
using std::complex;

#ifdef USE_MPI
#ifdef MPI_MISSING_BOOL_WORKAROUND
template<>
MPI_Datatype CommTrace<bool>::MPIType = MPI_INTEGER;
#else
template<>
MPI_Datatype CommTrace<bool>::MPIType = MPI::BOOL;
#endif
template<>
MPI_Datatype CommTrace<int>::MPIType = MPI_INTEGER;
template<>
MPI_Datatype CommTrace<float>::MPIType = MPI_FLOAT;
template<>
MPI_Datatype CommTrace<double>::MPIType = MPI_DOUBLE;
template<>
MPI_Datatype CommTrace<char>::MPIType = MPI_CHAR;
//template<>
//MPI_Datatype CommTrace<std::complex<double> >::MPIType = MPI_DOUBLE_COMPLEX;
template<>
MPI_Datatype CommTrace<std::complex<double> >::MPIType = MPI_DOUBLE;

template<>
int CommTrace<bool>::multiplicity = 1;
template<>
int CommTrace<int>::multiplicity = 1;
template<>
int CommTrace<float>::multiplicity = 1;
template<>
int CommTrace<double>::multiplicity = 1;
template<>
int CommTrace<std::complex<double> >::multiplicity = 2;
template<>
int CommTrace<char>::multiplicity = 1;

static MPI_Request nullReq;
static MPI_Status nullStat;
#endif

//------------------------------------------------------------------------------

void initCommunication(int &argc, char **&argv)
{

#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif

}

//------------------------------------------------------------------------------

void closeCommunication()
{

#ifdef USE_MPI
  MPI_Finalize();
#endif

}

//------------------------------------------------------------------------------

/*
#ifdef USE_MPI
void exit(int status)
{

  MPI_Abort(MPI_COMM_WORLD, status);

}
#endif
*/

//------------------------------------------------------------------------------

Communicator::Communicator()
#ifdef USE_MPI
  : pendReq(nullReq), reqStatus(nullStat)
#endif
{

#ifdef USE_MPI
  comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &numCPU);
  MPI_Comm_rank(comm, &thisCPU);
  nPendReq = 0;
#else
  numCPU = 1;
  thisCPU = 0;
#endif

  timer = 0;
  maxverbose = 0;

}

//------------------------------------------------------------------------------

#ifdef USE_MPI
Communicator::Communicator(MPI_Comm c1)
  : pendReq(nullReq), reqStatus(nullStat)
{

  comm = c1;
  MPI_Comm_size(comm, &numCPU);
  MPI_Comm_rank(comm, &thisCPU);
  nPendReq = 0;

  timer = 0;
  maxverbose = 0;

}
#endif

//------------------------------------------------------------------------------
// note: color+1 is used in order to make the routine compatible with the
// fortran communication library (which is restricted to 4 codes)

void Communicator::split(int color, int maxcolor, Communicator** c)
{

  int i;
  for (i=0; i<maxcolor; ++i)
    c[i] = 0;

#ifdef USE_MPI
  int rank;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm comm1;
#ifdef __INTEL_COMPILER
  MPI_Comm_split(comm, color + 1, 0, &comm1); //possible bug in openmpi1.4.3 compiled with intel compilerpro-12.0.2.137
#else
  MPI_Comm_split(comm, color+1, rank, &comm1);
#endif
  c[color] = new Communicator(comm1);

  int* leaders = new int[maxcolor];
  int* newleaders = new int[maxcolor];
  for (i=0; i<maxcolor; ++i)
    leaders[i] = -1;

  int localRank;
  MPI_Comm_rank(comm1, &localRank);
  if (localRank == 0)
    leaders[color] = rank;

  MPI_Allreduce(leaders, newleaders, maxcolor, MPI_INTEGER, MPI_MAX, comm); // May be a prob

  for (i=0; i<maxcolor; ++i) {
    if (i != color && newleaders[i] >= 0) {
      int tag;
      if (color < i)
	tag = maxcolor*(color+1)+i+1;
      else
	tag = maxcolor*(i+1)+color+1;
      MPI_Comm comm2;
      MPI_Intercomm_create(comm1, 0, comm, newleaders[i], tag, &comm2);
      c[i] = new Communicator(comm2);
    }
  }

  if (leaders)
    delete [] leaders;
  if (newleaders)
    delete [] newleaders;
#else
  c[color] = this;
#endif

}

Communicator *Communicator::merge(bool high)
{
#ifdef USE_MPI
  MPI_Comm newComm;
  MPI_Intercomm_merge(comm, high ? 1 : 0, &newComm);
  return new Communicator(newComm);
#endif
  return 0;
}

//------------------------------------------------------------------------------

int Communicator::remoteSize()
{

  int numRemote = 0;

#ifdef USE_MPI
  MPI_Comm_remote_size(comm, &numRemote);
#endif

  return numRemote;

}

//------------------------------------------------------------------------------

int Communicator::barrier()
{

  int ierr = 0;

#ifdef USE_MPI
  ierr = MPI_Barrier(comm);
#endif

  return ierr;

}

//------------------------------------------------------------------------------


void Communicator::printf(int verbose, const char *format, ...)
{

  if (thisCPU == 0 && verbose <= maxverbose) {
    va_list args;
    va_start(args, format);
    vfprintf(stdout, format, args);
    ::fflush(stdout);
    va_end(args);
  }

}

//------------------------------------------------------------------------------

void Communicator::fprintf(FILE *file, const char *format, ...)
{

  if (thisCPU == 0) {
    va_list args;
    va_start(args, format);
    vfprintf(file, format, args);
    ::fflush(file);
    va_end(args);
  }

}

//------------------------------------------------------------------------------

void Communicator::system(const char* command)
{
  if (thisCPU == 0)
    std::system(command);
}

//------------------------------------------------------------------------------

void Communicator::waitForAllReq()
{

#ifdef USE_MPI
  // Just making sure that reqStatus has an appropriate length
  MPI_Status *safe = reqStatus+nPendReq;

  if (safe == 0)
    exit(1);

  int nSuccess = MPI_Waitall(nPendReq, pendReq+0, reqStatus+0);

  if (nSuccess) {
    fprintf(stderr, "*** Error: unexpected success number %d\n", nSuccess);
    exit(1);
  }

  nPendReq = 0;
#endif

}

//------------------------------------------------------------------------------

void* operator new(size_t size, Communicator &) {
  void *a;
#ifdef USE_MPI
	MPI_Alloc_mem(size, MPI_INFO_NULL, &a);

#else // USE_MPI
	a = malloc(size);
#endif // USE_MPI
	if( !a ) {

	    std::bad_alloc ba;

	    throw ba;
	  }
	return a;
}

void operator delete(void *p, Communicator &) {
#ifdef USE_MPI
	MPI_Free_mem(p);
#else // USE_MPI
	free(p);
#endif // USE_MPI
}

void* operator new[](size_t size, Communicator &) {
	void *a;
#ifdef USE_MPI
        MPI_Alloc_mem(size, MPI_INFO_NULL, &a);
#else // USE_MPI
	// Changed by Adam and Julien 2010.07.30
	//	a = malloc(size);
	a = operator new[](size);
#endif // USE_MPI
	if( !a ) {

	    std::bad_alloc ba;

	    throw ba;
	  }
	return a;
}

void operator delete[](void *p, Communicator &) {
#ifdef USE_MPI
	MPI_Free_mem(p);
#else // USE_MPI
	// Changed by Adam and Julien 2010.07.30
	//	free(p);
	operator delete[] (p);
#endif // USE_MPI
}
