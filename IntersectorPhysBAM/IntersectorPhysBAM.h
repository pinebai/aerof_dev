#ifndef _INTERSECTORPHYSBAM_H_
#define _INTERSECTORPHYSBAM_H_

#include <list>
#include <map>
#include <set>
#include <string>

#include <Edge.h>
#include <LevelSet/LevelSetStructure.h>
#include <Vector.h>

#include "PHYSBAM_INTERFACE.h"
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>

#define MAXLINE 500

using std::pair;
using std::map;
using std::list;
using PhysBAM::PhysBAMInterface;
using PhysBAM::ARRAY;
using PhysBAM::PAIR;
using PhysBAM::VECTOR;
using PhysBAM::IntersectionResult;

class Vec3D;
class Communicator;
class IntersectorPhysBAM;
class FloodFill;
class SubDomain;
class EdgeSet;
class Timer;
class IoData;
class CrackingSurface;

template<class Scalar, int dim> class SVec;

class DistIntersectorPhysBAM : public DistLevelSetStructure {

  friend class IntersectorPhysBAM;

  typedef pair<int, int> iipair;
  typedef pair<int, bool> ibpair;
  typedef pair<iipair, ibpair> EdgePair;

  using DistLevelSetStructure::status;
  using DistLevelSetStructure::distance;
  using DistLevelSetStructure::is_swept;
  using DistLevelSetStructure::is_active;
  using DistLevelSetStructure::is_occluded;
  using DistLevelSetStructure::edge_intersects;

  using DistLevelSetStructure::edge_SI;
  using DistLevelSetStructure::xi_SI;
  using DistLevelSetStructure::eta_SI;
  using DistLevelSetStructure::nWall_SI;
  using DistLevelSetStructure::TriID_SI;

  using DistLevelSetStructure::xi_node;
  using DistLevelSetStructure::eta_node;
  using DistLevelSetStructure::nWall_node;
  using DistLevelSetStructure::TriID_node;
  protected:

    IoData &iod;

    int numStNodes, numStElems;
    int totStNodes, totStElems;
    bool gotNewCracking;
    double xMin, xMax, yMin, yMax, zMin, zMax; //a bounding box over the struct body
    DistSVec<double,3> *boxMax, *boxMin; //fluid node bounding boxes

    FloodFill* floodFill;

    // closest point to CFD grid-points
    DistVec<ClosestPoint> *closest;

    // struct node coords
    Vec3D *Xs;
    Vec3D *Xs0;
    Vec3D *Xs_n;
    Vec3D *Xs_np1;
    Vec3D *dXdSb;

    Vec<Vec3D> *solidX;   //pointer to Xs
    Vec<Vec3D> *solidXn;  //pointer to Xs_n
    Vec<Vec3D> *solidXnp1;//pointer to Xs_np1
    Vec<Vec3D> *solidX0;  //pointer to Xs0
    Vec<Vec3D> *solidXdS; //pointer to dXdSb

    int *faceID;
    double *porosity;

    // surface rotation
    int *surfaceID;
    int *rotOwn;

    int (*stElem)[3]; //structure elements (topology)

    Vec3D *Xsdot; //velocity

    CrackingSurface *cracking; //only a pointer.

    double interface_thickness;

	 bool externalSI;

	 bool isoThermalwall;	 
	 double Twall;

	 bool withViscousTerms;

  public:
    DistVec<int> *status0;  //previous node status
    DistVec<bool> *occluded_node0;// previous occluded node status
    DistVec<int> *is_swept_helper;

    double *triSize;
    Vec3D *triNorms;
    Vec3D *nodalNormal; //memory allocated only if interpolatedNormal == true
  
    DistSVec<double,3> *X; //pointer to fluid node coords
    DistSVec<double,3> *Xn; //pointer to fluid node coords at previous time

    // parameters from input
    bool interpolatedNormal;

    // sophisticated stuff...
    Communicator *com;
    IntersectorPhysBAM **intersector;
    PhysBAMInterface<double> *physInterface;
    Domain *domain;

    void buildSolidNormals();
    void getBoundingBox();

    void findActiveNodesUsingFloodFill(const DistVec<bool>& tId,const list<pair<Vec3D,int> >& points);
    void findActiveNodes(const DistVec<bool>& tId);

  public: //TODO: a lot of them can be moved to "protected".
    DistIntersectorPhysBAM(IoData &iod, Communicator *comm, int nNodes = 0, double *xyz = 0, int nElems = 0, int (*abc)[3] = 0, CrackingSurface *cs = 0);
    ~DistIntersectorPhysBAM();

    void init(char *meshfile, char *restartfile, double XScale);
    void init(int nNodes, double *xyz, int nElems, int (*abc)[3], char *restartSolidSurface);
    void setPorosity();
    void makerotationownership();
    void updatebc();

    void setdXdSb(int, double*, double*, double*); //
    void updateXb(double); //  

    EdgePair makeEdgePair(int,int,int);
    bool checkTriangulatedSurface();
    void initializePhysBAM();

    void initialize(Domain *, DistSVec<double,3> &X, DistSVec<double,3> &Xn, IoData &iod, DistVec<int>* point_based_id = 0,
                    DistVec<int>* oldStatus = 0);
    void updateStructure(double *xs, double *Vs, int nNodes, int (*abc)[3]=0);
   
    void updateCracking(int (*abc)[3]);
    void expandScope();
    void updatePhysBAMInterface(Vec3D *particles, int size,const DistSVec<double,3>& fluid_nodes,const bool fill_scope,const bool retry);
    int recompute(double dtf, double dtfLeft, double dts, bool findStatus, bool retry = false);

    LevelSetStructure & operator()(int subNum) const;

    PhysBAMInterface<double> &getInterface() { return *physInterface; }
    const Vec3D &getSurfaceNorm(int i) const {return triNorms[i]; }
    const Vec3D &getNodalNorm(int i) const {if (!nodalNormal) {fprintf(stderr,"ERROR: nodal normal not initialized!\n");exit(-1);} return nodalNormal[i];}
    Vec<Vec3D> &getStructPosition() { return *solidX; }
    Vec<Vec3D> &getStructPosition_0() { return *solidX0; }
    Vec<Vec3D> &getStructPosition_n() { return *solidXn; }
    Vec<Vec3D> &getStructPosition_np1() { return *solidXnp1; }
    DistVec<ClosestPoint> * getClosestPointsPointer() {return closest;}
    DistVec<ClosestPoint> & getClosestPoints() {return *closest;}
    void setStatus(DistVec<int> nodeTag) { *status = nodeTag; }
    int getNumStructNodes () { return numStNodes; }
    int getNumStructElems () { return numStElems; }
    int (*getStructElems())[3] { return stElem; }

    Vec<Vec3D> &getStructDerivative() { return *solidXdS; }

    int getSurfaceID(int k) {
      if (!surfaceID)
        return 0;
 
      if (k >=0 && k < numStNodes) {
	return surfaceID[k]; 
      }
      else {
        fprintf(stderr,"Error:: SurfaceID requested for invalid point.\n");
        exit(-1);
      }
    }
};

class IntersectorPhysBAM : public LevelSetStructure {
  friend class DistIntersectorPhysBAM;

  using LevelSetStructure::status;
  using LevelSetStructure::distance;
  using LevelSetStructure::is_swept;
  using LevelSetStructure::is_active;
  using LevelSetStructure::is_occluded;
  using LevelSetStructure::edge_intersects;

  using LevelSetStructure::edge_SI;
  using LevelSetStructure::xi_SI;
  using LevelSetStructure::eta_SI;
  using LevelSetStructure::nWall_SI;
  using LevelSetStructure::TriID_SI;

  using LevelSetStructure::xi_node;
  using LevelSetStructure::eta_node;
  using LevelSetStructure::nWall_node;
  using LevelSetStructure::TriID_node;

  public:
    static const int OUTSIDE = -2, UNDECIDED = -1, INSIDE = 0; //INSIDE: inside real fluid, OUTSIDE: not a fluid
    static int OUTSIDECOLOR;

    int locIndex,globIndex;
    int *locToGlobNodeMap;
 
    std::map<int,IntersectionResult<double> > CrossingEdgeRes;
    std::map<int,IntersectionResult<double> > ReverseCrossingEdgeRes;

    virtual CrackingSurface* getCrackingSurface() {
      return distIntersector.cracking; 
    }

    SubDomain &subD;
    EdgeSet &edges;
    DistIntersectorPhysBAM &distIntersector;

    std::set<int> *package;
    map<int,int> sub2pack;
    Connectivity &nodeToSubD;

    Vec<int> &status0; //<! status at the previous time-step.
    Vec<ClosestPoint> &closest;
    Vec<bool> &occluded_node0;//<! occluded node status at the previous time-step.
    int nFirstLayer;
    ARRAY<int> reverse_mapping,forward_mapping;
    ARRAY<VECTOR<double,3> > xyz;
    ARRAY<VECTOR<double,3> > xyz_n;

    bool testIsActive(double t, int n) const {return (status[n] >= 0 && status[n]!=OUTSIDECOLOR);}
    int hasCloseTriangle(SVec<double,3>& X,SVec<double,3>& Xn,SVec<double,3> &boxMin, SVec<double,3> &boxMax, Vec<bool> &tId);
    int findIntersections(SVec<double,3>& X,Vec<bool>& tId,Communicator&);
    int computeSweptNodes(SVec<double,3>& X, Vec<bool>& tId,Communicator&,const double dt);

	 void reFlagRealNodes(SVec<double,3>& X, Vec<bool> *bk_ISactive = 0); //d2d
	 void ComputeSIbasedIntersections(int iSub, SVec<double,3>& X, SVec<double,3> &boxMin, SVec<double,3> &boxMax, bool withViscousTerms); //d2d

  public:
    IntersectorPhysBAM(SubDomain &, SVec<double, 3> &X, Vec<int> &status0, Vec<ClosestPoint> &closest,
                       Vec<bool>& occluded_node0,DistIntersectorPhysBAM &);
    virtual ~IntersectorPhysBAM();
    int numOfFluids() {return distIntersector.numOfFluids();}

    void reset(const bool findStatus,const bool retry); //<! set status0=status and reset status and nFirstLayer.
    /** find intersections for each edge that has nodes with different statuses */
    double isPointOnSurface(Vec3D pt, int N1, int N2, int N3);
    /** check the distance of apoint to a surface defined by a triangle. (used for debug only) */ 
    void printFirstLayer(SubDomain& sub, SVec<double,3>& X, int TYPE = 1);

    LevelSetResult getLevelSetDataAtEdgeCenter(double t, int l, bool i_less_j, double *Xr=0, double *Xg=0);
    void findNodesNearInterface(SVec<double, 3>&, SVec<double, 3>&, SVec<double, 3>&) {}

    bool isNearInterface(double t, int n) const {return closest[n].nearInterface();}
    bool withCracking() const {return distIntersector.cracking ? true : false;}

	 void xWallWithSI(int n, Vec3D &xWall);
	 void vWallWithSI(int n, Vec3D &vWall);
	 bool xWallNode(int i, Vec3D &xWall);
	 bool vWallNode(int i, Vec3D &vWall);

	 bool getTwall(double &Tw)
	 { 
		 Tw = distIntersector.Twall; 
		 return distIntersector.isoThermalwall;
	 }

  private:
    void addToPackage(const int i,const int candidate);
    void findNodeClosestPoint(const int nodeId, Vec3D &x0, ARRAY<int> &cand);
    double project(Vec3D x0, int tria, double& xi1, double& xi2) const;
    double edgeProject(Vec3D x0, Vec3D &xA, Vec3D &xB, double &alpha) const;
    double edgeProject(Vec3D x0, int n1, int n2, double &alpha) const;

	 double piercing(Vec3D x0, int tria, double xi[3]);

    void derivativeOFnormal(Vec3D  xA, Vec3D  xB, Vec3D  xC, 
			    Vec3D dxA, Vec3D dxB, Vec3D dxC, 
			    Vec3D &dnds);

    double derivativeOFalpha(Vec3D  xA, Vec3D  xB, Vec3D  xC, 
	  		     Vec3D dxA, Vec3D dxB, Vec3D dxC,
			     Vec3D X1, Vec3D X2);

    double testAlpha(Vec3D, Vec3D, Vec3D, Vec3D, Vec3D);

};

#endif
