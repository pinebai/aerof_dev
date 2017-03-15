#ifndef _LEVEL_SET_STRUCTURE_H_
#define _LEVEL_SET_STRUCTURE_H_
 
#include "Vector3D.h"
#include <DistVector.h>
#include <Vector.h>

class Domain;
class IoData;
template<class Scalar, int dim>
class DistSVec;
template <class Scalar> class Vec;
template <class Scalar, int dim> class SVec;
template <class Scalar> class DistVec;

/** Structure used to return levelset information */
struct LevelSetResult {
  double alpha;
  double xi[3];
  int trNodes[3];
  Vec3D gradPhi;
  Vec3D normVel; //NOTE: this is the velocity, NOT normal velocity.
  int structureType;
  double porosity;
  double actuatorDiskPressureJump;
	int actuatorDiskMethod;
  double gamma;
  bool isCorrectedMethod;
  int actuatorDiskReconstructionMethod;

  Vec3D dnds;
  double dads;

  LevelSetResult() {
    alpha = xi[0] = xi[1] = -1.0;
    trNodes[0] = trNodes[1] = trNodes[2] = -1;
    gradPhi = normVel = dnds = 0.0;
      structureType = 1;
    porosity = 0.0;
    actuatorDiskPressureJump = 0.0;
	  actuatorDiskMethod = 1;
    actuatorDiskReconstructionMethod = -1;
    gamma = 0.0;//error value
    isCorrectedMethod = false;
    dads = 0.0;
   }
  
   LevelSetResult(double gpx, double gpy, double gpz,
                  double nvx, double nvy, double nvz) :
		  gradPhi(gpx, gpy, gpz), normVel(nvx, nvy, nvz) {
		    gradPhi *= 1.0/gradPhi.norm();
		    alpha = -1.0;
		    xi[0] = xi[1] = xi[2] = -1.0;
		    trNodes[0] = trNodes[1] = trNodes[2] = -1;
		    porosity = 0.0;
	        structureType = 1;
	   actuatorDiskMethod = 1;
                    actuatorDiskPressureJump = 0.0;
                    actuatorDiskReconstructionMethod = -1;
                    gamma = 0.0;//error value
                    isCorrectedMethod = false;
		    dnds = 0.0;
		    dads = 0.0;
		   }

   class iterator {
	   double *xip;
	   int *nodep;
   public:
	   iterator(double *x, int *n) : xip(x), nodep(n) { }
	   iterator & operator++() { xip++; nodep++; return (*this); }
	   bool operator !=(const iterator &it) const {
		   return xip != it.xip || nodep != it.nodep;
	   }
	   double Ni() const { return *xip; }
	   int nodeNum() const { return *nodep; }
   };

   iterator begin() { return iterator(xi, trNodes); }
   iterator end() { return iterator(xi+3, trNodes+3); }
};

/** Storing the closest point on the interface */
struct ClosestPoint {
  int mode; //-2: unknown, -1: far from interface, 0: face, 1: edge, 2: vertex.
  double dist; //this is the unsigned distance. always >= 0.
  int tracker[2]; // for mode=0: tracker[0] = tria Id; for mode=1: the two vertices; for mode=2: the vertex
  double xi1,xi2; // local coordinates.
  ClosestPoint() {mode=-2;}
  bool known() {return mode!=-2;}
  bool nearInterface() {return mode!=-1;}
};

/** Abstract class for finding levelset information */
class LevelSetStructure {
  protected:
    Vec<int> &status;
    Vec<double> &distance;
    Vec<bool> &is_swept;
    Vec<bool> &is_active;
    Vec<bool> &is_occluded;
    Vec<bool> &edge_intersects;

	 // Surrogate based quantities
    Vec<bool>   &edge_SI; // ID of the edge that contains the SI 
	 Vec<double>   &xi_SI; // baricentric cooardinates of the the 
	 Vec<double>  &eta_SI; // intersection point on the wall 
	 Vec<Vec3D> &nWall_SI; // wall normal at the Wall
	 Vec<int>   &TriID_SI; // ID of the element of the embedded 
	                       //    surface intersected

	 // Node based quantities
	 Vec<double>   &xi_node; // baricentric cooardinates of the the 
	 Vec<double>  &eta_node; // intersection point on the wall 
	 Vec<Vec3D> &nWall_node; // wall normal at the Wall
	 Vec<int>   &TriID_node; // ID of the element of the embedded 
	                         //    surface intersected

  public:
LevelSetStructure(Vec<int>& status,
						Vec<double>& distance,
						Vec<bool>& is_swept,
						Vec<bool>& is_active,
						Vec<bool>& is_occluded,
						Vec<bool>& edge_intersects,
						Vec<bool>&  edge_SI,
						Vec<double>& xi_SI,
						Vec<double>& eta_SI,
						Vec<Vec3D>& nWall_SI,
						Vec<int>&   TriID_SI,
						Vec<double>& xi_node,
						Vec<double>& eta_node,
						Vec<Vec3D>& nWall_node,
						Vec<int>& TriID_node): status(status),
		                                   distance(distance),
		                                   is_swept(is_swept),
		                                   is_active(is_active),
		                                   is_occluded(is_occluded),
                                    	  edge_intersects(edge_intersects),
                                    	  edge_SI(edge_SI),
		                                   xi_SI(xi_SI),
		                                   eta_SI(eta_SI),
		                                   nWall_SI(nWall_SI),
                                         TriID_SI(TriID_SI),
		                                   xi_node(xi_node),
		                                   eta_node(eta_node),
		                                   nWall_node(nWall_node),
		                                   TriID_node(TriID_node)

    {}
    virtual ~LevelSetStructure()
    {}

    /** returns the normal and normal velocity at intersection between edge ni, nj and structure
     *
     * If ni, nj is not an edge of the fluid mesh, result is undefined.
     * */
    virtual LevelSetResult
      getLevelSetDataAtEdgeCenter(double t, int l, bool i_less_j, double *Xr=0, double *Xg=0) = 0;
    virtual bool withCracking() const = 0;
    virtual bool isNearInterface(double t, int n) const = 0;

    void forceOccluded(double t, int n) const                { is_swept[n] = true; is_occluded[n] = true; }
    int fluidModel(double t, int n) const                 { return status[n]; }
    double distToInterface(double t, int n) const         { return distance[n]; } 
    bool isSwept(double t, int n) const                   { return is_swept[n]; }
    bool isActive(double t, int n) const                  { return is_active[n]; }
    bool isOccluded(double t, int n) const                { return is_occluded[n]; }
    bool edgeIntersectsStructure(double t, int eij) const { return edge_intersects[eij]; }
    void computeSwept(Vec<int> &swept){
        for(int i = 0; i < swept.size(); ++i)
            swept[i] = is_swept[i] ? 1 : 0;
    }

	 /* Retrun true if the edge 'n' contains a SI */
	 bool edgeWithSI( int n) const { return edge_SI[n];  } 
	 int  eWallWithSI(int n) const { return TriID_SI[n]; }
 
	 void nWallWithSI(int n, Vec3D &nWall) { nWall = nWall_SI[n]; }
	 void nWallNode(  int i, Vec3D &nWall) { nWall = nWall_node[i]; }

	 virtual void xWallWithSI(int n, Vec3D &xWall) = 0;
	 virtual void vWallWithSI(int n, Vec3D &vWall) = 0;
	 virtual bool xWallNode(  int i, Vec3D &xWall) = 0;
	 virtual bool vWallNode(  int i, Vec3D &vWall) = 0;

	 virtual bool getTwall(double &Tw) = 0;

    Vec<int> & getStatus() { return status; }

    virtual class CrackingSurface* 
      getCrackingSurface() { return NULL; }

    virtual double isPointOnSurface(Vec3D, int, int, int) = 0;

    virtual int numOfFluids() = 0; 

    virtual void findNodesNearInterface(SVec<double,3>&, SVec<double,3>&, SVec<double,3>&) = 0;

};

class DistLevelSetStructure {
  protected:
    DistVec<int> *status;
    DistVec<double> *distance;
    DistVec<bool> *is_swept;
    DistVec<bool> *is_active;
    DistVec<bool> *is_occluded;
    DistVec<bool> *edge_intersects;
    DistVec<bool> *edge_SI;   // d2d
	 DistVec<int>  *TriID_SI;
	 DistVec<double> *xi_SI;
	 DistVec<double> *eta_SI;
	 DistVec<Vec3D> *nWall_SI;
	 DistVec<int>  *TriID_node;
	 DistVec<double> *xi_node;
	 DistVec<double> *eta_node;
	 DistVec<Vec3D> *nWall_node;

  protected:
    int numLocSub;
    int numFluid;
    bool with_sensitivity;

  public:
    DistLevelSetStructure()
		 : status(0), distance(0), is_swept(0), is_active(0), is_occluded(0), edge_intersects(0), 
		 edge_SI(0), xi_SI(0), eta_SI(0), nWall_SI(0), TriID_SI(0), xi_node(0), eta_node(0), nWall_node(0), TriID_node(0)
    {}
    virtual ~DistLevelSetStructure()
    {delete status;delete distance;delete is_swept;delete is_active;delete is_occluded;delete edge_intersects; 
	  delete edge_SI; delete nWall_SI; delete xi_SI; delete eta_SI; delete TriID_SI; 
	  delete nWall_node; delete xi_node; delete eta_node; delete TriID_node;}

    int numOfFluids() {return numFluid;}
    void setNumOfFluids(int nf) {numFluid = nf;}
    virtual void initialize(Domain *, DistSVec<double,3> &X, DistSVec<double,3> &Xn, IoData &iod, DistVec<int> *point_based_id = 0, DistVec<int>* oldStatus = 0) = 0;
    virtual LevelSetStructure & operator()(int subNum) const = 0;

    void getSwept(DistVec<int>& swept){
        for (int iSub=0; iSub<numLocSub; iSub++)
            (*this)(iSub).computeSwept((swept)(iSub));
    }

    DistVec<int> & getStatus()            const { return *status; }
    DistVec<double> & getDistance()       const { return *distance; }
    DistVec<bool> & getIsSwept()          const { return *is_swept; }
    DistVec<bool> & getIsActive()         const { return *is_active; }
    DistVec<bool> & getIsOccluded()       const { return *is_occluded; }
    DistVec<bool> & getIntersectedEdges() const { return *edge_intersects; }

    virtual DistVec<ClosestPoint> &getClosestPoints() = 0;
    virtual DistVec<ClosestPoint> *getClosestPointsPointer() = 0;
    virtual void setStatus(DistVec<int> nodeTag) = 0;                                

    virtual void updateStructure(double *Xs, double *Vs, int nNodes, int(*abc)[3]=0) = 0;
    virtual int recompute(double dtf, double dtfLeft, double dts, bool findStatus, bool retry = false) = 0;
    virtual Vec<Vec3D> &getStructPosition() = 0;
    virtual Vec<Vec3D> &getStructPosition_0() = 0;
    virtual Vec<Vec3D> &getStructPosition_n() = 0;
    virtual Vec<Vec3D> &getStructPosition_np1() = 0;
    virtual Vec<Vec3D> &getStructDerivative() = 0;

    virtual int getNumStructNodes() = 0;
    virtual int getNumStructElems() = 0;
    virtual int (*getStructElems())[3] = 0;
    virtual int getSurfaceID(int) = 0;

    virtual void updateXb(double) = 0;
    virtual void setdXdSb(int, double*, double*, double*) = 0;
};

#endif
