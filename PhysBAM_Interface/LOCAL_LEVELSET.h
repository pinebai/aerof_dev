#ifndef __LOCAL_LEVELSET__
#define __LOCAL_LEVELSET__
#include<cstdio>

class LocalLevelSet {
public:
  //NOTE: trId starts from 0.
  virtual bool purelyPhantomPhysBAM(int trId) {
    fprintf(stderr,"ERROR: function purelyPhantom is not implemented in LoclLevelSet!\n");
    return false;
  }
  //NOTE: trId starts from 0.
  virtual double getPhiPhysBAM(int trId, double xi1, double xi2, bool* hasCracked=0, bool debug = false) {
    fprintf(stderr,"ERROR: function getPhi is not implemented in LocalLevelSet!\n");
    return 0; 
  }

}; 

#endif
