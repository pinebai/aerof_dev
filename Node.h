#ifndef _NODE_H_
#define _NODE_H_

#include <Vector.h>

class BinFileHandler;

//------------------------------------------------------------------------------

class NodeSet : public SVec<double,3> {
  
public:

  NodeSet(int n) : SVec<double,3>(n) {}
  ~NodeSet() {}

  int read(BinFileHandler &, int, int (*)[2], int *, int *);

};

//------------------------------------------------------------------------------

#endif
