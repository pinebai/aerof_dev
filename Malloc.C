#include <sys/types.h>
#include <sys/mman.h>
#include <malloc.h>
#include <fcntl.h>

#ifdef OLD_STL
#include <map.h>
#else
#include <map>
using namespace std;
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

/*
#define _OPENMP
int omp_get_thread_num() { return 0; }
*/

#ifdef sgi
int maxThread = 16;
void **arenas = 0;
//map<void *, int> memMap;
int initSize=2*1024*1024;
int growth = 32*1024*1024;
int zeroFd;

ulock_t mmapLock = 0;
usptr_t *usPtr;

//------------------------------------------------------------------------------

void *arenaGrow(size_t size, void *)
{

  void *gr = mmap(0, size, PROT_READ|PROT_WRITE , MAP_PRIVATE, zeroFd, 0);

  if (gr == 0) { fprintf(stderr,"Out of memory\n"); exit(-1); }

#ifdef _OPENMP
  int tNum = omp_get_thread_num();
  //fprintf(stderr, "Grew  at %p for %d\n", gr, tNum);
  //if(tNum > maxThread)
  //tNum %= maxThread;
  //ussetlock(mmapLock);
  //memMap[gr] = tNum;
  //usunsetlock(mmapLock);
#endif

  return gr;

}

//------------------------------------------------------------------------------

#ifdef MYMALLOC
void * operator new(size_t size)
{

#ifdef _OPENMP
  if (arenas == 0) {
    usPtr = usinit("/dev/zero");
    mmapLock = usnewlock(usPtr);
    zeroFd = open("/dev/zero",O_RDWR);
    arenas = (void **)malloc(maxThread*sizeof(void *));
    int iThread;
    for (iThread = 0; iThread < maxThread; ++iThread) {
      arenas[iThread] = malloc(initSize);
      acreate(arenas[iThread], initSize, 0, 0, &arenaGrow);
      amallopt(M_BLKSZ, growth, arenas[iThread]);
    }
  }

  int tNum = omp_get_thread_num();

  if (tNum > maxThread)
    tNum %= maxThread;

  if (size > 1000000) {
    ussetlock(mmapLock);
    //fprintf(stderr, "Large allocation in %d %d\n", tNum, size);
  }

  void *res = amalloc(size,arenas[tNum]);

  if(size > 1000000) {
    //fprintf(stderr, "Done with large allocation\n");
    usunsetlock(mmapLock);
  }
  return res;
#else
  return malloc(size);
#endif

}

//------------------------------------------------------------------------------

void operator delete(void *p)
{

#ifdef _OPENMP
  int tNum = omp_get_thread_num();
  if (tNum > maxThread)
    tNum %= maxThread;
  afree(p,arenas[tNum]);
#else
  free(p);
#endif

}
#endif

//------------------------------------------------------------------------------

void setNumArenas(int na)
{

  if (na < maxThread) return;

  void **newA = (void **) malloc(sizeof(void *)*na);

  int i;
  for (i = 0; i < maxThread; ++i)
    newA[i] = arenas[i];

  for (i = maxThread; i < na; ++i) {
    newA[i] = malloc(initSize);
    acreate(newA[i], initSize, 0, 0, &arenaGrow);
    amallopt(M_BLKSZ, growth, newA[i]);
  }

  arenas = newA;
  maxThread = na;

}

//------------------------------------------------------------------------------

#endif
