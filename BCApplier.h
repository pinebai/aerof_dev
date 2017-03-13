#ifndef _BC_APPLIER_H_
#define _BC_APPLIER_H_

#include <DistVector.h>
#include <IoData.h>

#ifdef OLD_STL
#include <list.h>
#include <map.h>
#else
#include <list>
#include <map>
using std::list;
using std::map;
#endif


class MatchNodeSet;
class Domain;
class SurfaceData;

struct ProjData {
  int node;
  double n[3];
  ProjData(int nd, double nrm[3]) {
    node = nd;
    n[0] = nrm[0];
    n[1] = nrm[1];
    n[2] = nrm[2];
  }
  void apply(double u[3]) {
    double dd = n[0]*u[0]+n[1]*u[1]+n[2]*u[2];
    u[0] -= dd*n[0];
    u[1] -= dd*n[1];
    u[2] -= dd*n[2];
  }
  void rotateNormal(double dRot[3][3]) {
    double v[3] = {n[0],n[1],n[2]};
    for (int j = 0; j < 3; j++)
      n[j] = dRot[j][0]*v[0] + dRot[j][1]*v[1] + dRot[j][2]*v[2];
                                                                                                                   
  }
};

class BCApplier {
  private:
    Domain* domain;
    int numLocSub;
    int **dofType;
    int **dofTypeStep1; 
    int **dofTypeStep2;
    map<int,SurfaceData*>& surfaceMap;
    list<ProjData> *proj;
    const DefoMeshMotionData::SlidingSurfaceTreatment &slidingSurfaceTreatment;

  public:
    BCApplier(int nLocSub, Domain *dom, IoData& ioData);
    ~BCApplier();

    int** getDofType() { return(dofType); } 
 
    // Adds a projection on a node. 
    // WARNING: Projections should be on orthogonal
    // planes to guarantee commutativity.
    void addProj(int iSub, int node, double n[3]);

    void rotateProjNormal(double dRot[3][3]);
       
    template<int dim> void applyP(DistSVec<double,dim> &X);
    template<int dim> void applyPt(DistSVec<double,dim> &X) { applyP(X); }

    template<int dim> void applyD(DistSVec<double,dim> &X);
    template<int dim> void applyDt(DistSVec<double,dim> &X) { applyD(X); }
    template<int dim> void applyD2(DistSVec<double,dim> &X, double dX[dim]);

    template<int dim> void applyPD(DistSVec<double,dim> &X) { applyD(X); applyP(X); }
    template<int dim> void applyPDt(DistSVec<double,dim> &X) { applyPt(X); applyDt(X); }


    void setDofType(MatchNodeSet** matchNodes=0);
    void setEmbeddedALEDofType(MatchNodeSet** matchNodes=0);
    void print();

};

#ifdef TEMPLATE_FIX
#include <BCApplier.C>
#endif

#endif
