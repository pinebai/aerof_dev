#ifndef _NONLINEAR_ROM_ONLINE_II_H_
#define _NONLINEAR_ROM_ONLINE_II_H_

#include <NonlinearRom.h>

template <int dim>
class NonlinearRomOnlineII : public NonlinearRom<dim> {

  protected:


  public:

  NonlinearRomOnlineII(Communicator *_com, IoData &_ioData, Domain &_domain,
                       std::vector<double>* _weights=NULL);
  ~NonlinearRomOnlineII();

  bool updateBasis(int, DistSVec<double, dim> &, Vec<double>* coords = NULL);
  bool updateBasisSimple(int, DistSVec<double, dim> &);
  bool updateBasisFastExact(int, DistSVec<double, dim> &, Vec<double>*);
  bool updateBasisFastApprox(int, DistSVec<double, dim> &);
  void projectSwitchStateOntoAffineSubspace(int, int, DistSVec<double, dim> &, Vec<double> &);

  void appendNonStateDataToBasis(int, const char *, bool relProjError = false); 
  void readClusteredOnlineQuantities(int);

  void readClosestCenterInfoModelII();

  void appendVectorToBasis(DistSVec<double, dim>&, int numVec = 0);
 
};

#include "NonlinearRomOnlineII.C"
#endif
