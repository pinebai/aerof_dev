#ifndef _MULTIGRID_LEVEL_SET_STRUCTURE_H_
#define _MULTIGRID_LEVEL_SET_STRUCTURE_H_
 
#include "LevelSetStructure.h"

template <typename Scalar> class MultiGridLevel;
class SubDomain;
class EdgeSet;

/** Abstract class for finding levelset information */
class MultiGridLevelSetStructure : public LevelSetStructure {
  protected:
    Vec<int> &status;
    Vec<double> &distance;
    Vec<bool> &is_swept;
    Vec<bool> &is_active;
    Vec<bool> &is_occluded;
    Vec<bool> &edge_intersects;

    Vec<Vec3D> &surfaceNormals;

    class DistMultiGridLevelSetStructure& distLSS;

    LevelSetStructure* parent;

    SubDomain &subD;
    EdgeSet &edges;

    int mySub;

    MultiGridLevel<double>* myLevel;

  public:
    MultiGridLevelSetStructure(class DistMultiGridLevelSetStructure& lss,
			       SubDomain& sub,
			       Vec<int>& status,Vec<double>& distance,Vec<bool>& is_swept,
			       Vec<bool>& is_active,Vec<bool>& is_occluded,
			       Vec<bool>& edge_intersects,
			       Vec<Vec3D>& surfaceNormals,
			       LevelSetStructure* parent,
			       int mySub,MultiGridLevel<double>* myLevel);

    ~MultiGridLevelSetStructure()
    {}

    /** returns the normal and normal velocity at intersection between edge ni, nj and structure
     *
     * If ni, nj is not an edge of the fluid mesh, result is undefined.
     * */
    LevelSetResult
      getLevelSetDataAtEdgeCenter(double t, int l, bool i_less_j, double *Xr=0, double *Xg=0);

    int numOfFluids();

    void recompute();

    void computeEdgeCrossing(SVec<double,3>& nodeNormals);

    bool withCracking() const { return false; }

	 void xWallWithSI(int n, Vec3D &xWall) {exit(-1);}
	 void vWallWithSI(int n, Vec3D &vWall) {exit(-1);}

	 bool xWallNode(int i, Vec3D &xWall) {exit(-1);}
	 bool vWallNode(int i, Vec3D &vWall) {exit(-1);}
	 bool getTwall(double &Tw) {exit(-1);}

    bool isNearInterface(double, int) const { return false; }

    double isPointOnSurface(Vec3D, int, int, int) { return 0.0; }

    void findNodesNearInterface(SVec<double, 3>&, SVec<double, 3>&, SVec<double, 3>&) { }

    void setdXdSb(int, double*, double*, double*) {}

};

class DistMultiGridLevelSetStructure : public DistLevelSetStructure {

  friend class MultiGridLevelSetStructure;

  // LSS for the parent level
  DistLevelSetStructure* parent;

  Communicator *com;

  MultiGridLevel<double>* myLevel;

  int numLocSub;

  MultiGridLevelSetStructure** subLSS;

  Domain* domain;

  DistVec<ClosestPoint>* dummycp;

  DistVec<Vec3D>* surfaceNormals;

  public:
    DistMultiGridLevelSetStructure(IoData &iod, Communicator *comm,
				   DistLevelSetStructure* parent,
				   MultiGridLevel<double>*);

    virtual ~DistMultiGridLevelSetStructure()
      {}

    void initialize(Domain *, DistSVec<double,3> &X, DistSVec<double,3> &Xn, IoData &iod, DistVec<int> *point_based_id = 0, DistVec<int>* oldStatus = 0);
    LevelSetStructure & operator()(int subNum) const;

    void setStatus(DistVec<int> nodeTag) { }

    void updateStructure(double *Xs, double *Vs, int nNodes, int(*abc)[3]=0) {

      parent->updateStructure(Xs, Vs, nNodes, abc);
    }

    int recompute(double dtf, double dtfLeft, double dts, bool findStatus, bool retry = false);


    Vec<Vec3D> &getStructPosition()  { return parent->getStructPosition(); }
    Vec<Vec3D> &getStructPosition_0()  { return parent->getStructPosition_0(); }
    Vec<Vec3D> &getStructPosition_n() { return parent->getStructPosition_n(); }
    Vec<Vec3D> &getStructPosition_np1() { return parent->getStructPosition_np1(); }
    Vec<Vec3D> &getStructDerivative() { return parent->getStructDerivative(); }

    int getNumStructNodes() { return parent->getNumStructNodes(); }
    int getNumStructElems() { return parent->getNumStructElems(); }
    int (*getStructElems())[3]  { return parent->getStructElems(); }

    int getSurfaceID(int k) { return parent->getSurfaceID(k); }

    virtual DistVec<ClosestPoint> &getClosestPoints() { return *dummycp; }
    virtual DistVec<ClosestPoint> *getClosestPointsPointer() { return dummycp; }

    void setdXdSb(int, double*, double*, double*){}
    void updateXb(double) {} //  
    void testAlpha(Vec3D, Vec3D, int, int, int, double*, double){}
};

#endif
