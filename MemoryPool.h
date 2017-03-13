#ifndef _MEMORY_POOL_H_
#define _MEMORY_POOL_H_

#define MAX_ACTIVE_PTR 10

#include <ResizeArray.h>

//------------------------------------------------------------------------------

class MemoryPool {

  int count;

  ResizeArray<bool> status;
  ResizeArray<int> size;
  ResizeArray<void *> ptr;

public:

  MemoryPool();
  ~MemoryPool() {}

  void set(int, void *);
  void *request(int);

};

//------------------------------------------------------------------------------

#endif
