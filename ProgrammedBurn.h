/* ProgrammedBurn.h

*/

#pragma once

#include <DistVector.h>
#include <VarFcn.h>

class ProgrammedBurn {

 public:
  
  ProgrammedBurn(IoData& ioData, DistSVec<double,3>* _nodeSet);
  ProgrammedBurn(IoData& ioData, SVec<double,1>* _nodeSet);

  ~ProgrammedBurn();

  //template<int dimLS>
  //  void setLevelSet(double t, DistSVec<double,dimLS>& LS, int lsid);

  static int countBurnableFluids(IoData& ioData);

  void setFluidIds(double t, DistVec<int>& fluidIds,DistSVec<double,5>& U);
  void setFluidIds(double t, Vec<int>& fluidIds,SVec<double,5>& U);

  bool isDetonationInterface(int i, int j,int& tag) const;

  void getDetonationNormal(int tag,int i,int j, double xmid[3], double gradphi[3]);

  template <int dim>
    void setCurrentTime(double t,VarFcn* vf,DistSVec<double,dim>& U,DistVec<int>& fid,DistVec<int>& fidn);

  // 1D version!!
  template <int dim>
    void setCurrentTime(double t,VarFcn* vf,SVec<double,dim>& U,Vec<int>& fid,Vec<int>& fidn);

  bool isBurnedEOS(int eos,int& tag) const;

  bool isUnburnedEOS(int eos,int& tag) const;

  int getBurnedEOS(int tag) const;
  int getUnburnedEOS(int tag) const;

  bool nodeInside(int tag,int sub, int i);
  bool nodeInside(int tag,int i);

  bool numberOfBurns() { return myBurns.size(); }
  bool isIgnited(int tag) const { return myBurns[tag].ignited; }
  bool isFinished(int tag) const { return myBurns[tag].finished; }

  //bool isBurned(int iSub,int i) { return bBurned->subData(iSub)[i]; }

  // Static methods.  Compute assorted values associated with a programmed burn
  static void computeChapmanJouguetStateJWL(double A1,double A2,double R1,double R2,double omega, // JWL parameters
					    double p0, double rho0, double e0,
					    double& rho_cj,double& p_cj, double& e_cj, double& s);

  static void computeChapmanJouguetStatePG(double gamma, // Perfect Gas gamma
					   double p0, double rho0, double e0,
					   double& rho_cj,double& p_cj, double& e_cj, double& s);
  
  
 private:

  struct Burn {
    Burn() : x0subdom(-1), x0id(-1), ignited(false) {}
    ProgrammedBurnData* pgData;
    double x0[3];
    int x0subdom, x0id;
    bool ignited;
    bool finished;
  };

  template <int dim>
    void setCJInitialState(Burn& B,VarFcn* vf,DistSVec<double,dim>& U,DistVec<int>&,DistVec<int>&);
  template <int dim>
    void setCJInitialState(Burn& B,VarFcn* vf,SVec<double,dim>& U,Vec<int>&,Vec<int>&);
  
  std::vector<Burn> myBurns;

  void computeNearestNode(const double x0n[3], double x0[3], int& x0subdom, int& x0id);

  double lastTime;

  const DistInfo* distInfo;

  DistSVec<double,3>* nodeSet;
  SVec<double,1>* nodeSet0;

};

#ifdef TEMPLATE_FIX
#include "ProgrammedBurn.C"
#endif
