#ifndef _NONLINEAR_ROM_ONLINE_III_H_
#define _NONLINEAR_ROM_ONLINE_III_H_

#include <NonlinearRom.h>

template <int dim>
class NonlinearRomOnlineIII : public NonlinearRom<dim> {

  protected:


  public:

  NonlinearRomOnlineIII(Communicator *_com, IoData &_ioData, Domain &_domain,
                        std::vector<double>* _weights=NULL);
  ~NonlinearRomOnlineIII();

  bool updateBasis(int, DistSVec<double, dim> &, Vec<double>* coords = NULL);
  bool updateBasisFastExact(int, DistSVec<double, dim> &, Vec<double>*);
  bool updateBasisFastApprox(int, DistSVec<double, dim> &);
  void projectSwitchStateOntoAffineSubspace(int, int, DistSVec<double, dim> &, Vec<double> &);

  void appendNonStateDataToBasis(int, const char *, bool relProjError = false); 
  void readClusteredOnlineQuantities(int);
  void readClosestCenterInfoModelIII();

};

#include "NonlinearRomOnlineIII.C"
#endif
