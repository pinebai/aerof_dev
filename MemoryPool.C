#include <cstdlib>
#include <cstdio>

#include <MemoryPool.h>

//------------------------------------------------------------------------------

MemoryPool::MemoryPool() : status(false), size(0), ptr(0)
{ 

  count = 0;

}

//------------------------------------------------------------------------------

void MemoryPool::set(int nbytes, void *p)
{

  bool exist = false;

  for (int i=0; i<count; ++i) {
    if (size[i] == 0) {
      status[i] = false;
      size[i] = nbytes;
      ptr[i] = p;
      exist = true;
      break;
    }
  }

  if (!exist) {
    status[count] = false;
    size[count] = nbytes;
    ptr[count] = p;
    ++count;
  }

}

//------------------------------------------------------------------------------

void *MemoryPool::request(int nbytes)
{

  void *p = 0;

  for (int i=0; i<count; ++i) {
    if (!status[i] && (nbytes <= size[i])) {
      p = ptr[i];
      status[i] = true;
      break;
    }
  }

  if (!p) {
    double s = nbytes / (1024.0 * 1024.0);
    fprintf(stderr, "*** Warning: could not find %3.2f MB of memory\n", s);
  }
  
  return p;

}

//------------------------------------------------------------------------------
