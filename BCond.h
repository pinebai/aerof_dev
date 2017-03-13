#ifndef _BCOND_H_
#define _BCOND_H_

#include <BinFileHandler.h>

//------------------------------------------------------------------------------

struct BCond 
{
  int nnum;
  int dofnum;
  double val;
};

//------------------------------------------------------------------------------

class BCondSet {

  int numBC;
  BCond *bcs;
	
public:

  BCondSet(int _numBC=0) { numBC = _numBC; bcs = (numBC) ? new BCond[numBC] : 0; }
  ~BCondSet() { if (bcs) delete[] bcs; }

  BCond &operator[](int i) const { return bcs[i]; }

  void read(BinFileHandler &);
  int size() const { return numBC; }
  void print();

};

//------------------------------------------------------------------------------

#endif
