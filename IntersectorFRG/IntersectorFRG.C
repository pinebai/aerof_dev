#include <cstdio>
#include <iostream>
#include "IntersectorFRG.h"
#include "Vector3D.h"
#include "Communicator.h"
#include <Domain.h>
#include <IoData.h>
#include <Vector.h>
#include <DistVector.h>
#include <Timer.h>
#include "parser/Assigner.h"
#include "Geometry/KDTree.h"
#include <Connectivity.h>
#include <queue>

#include "TsRestart.h"

using std::pair;
using std::map;
using std::list;
using std::set;

typedef pair<int, int> iipair;
typedef pair<int, bool> ibpair;
typedef pair<iipair, ibpair> EdgePair;

const int IntersectorFRG::UNDECIDED, IntersectorFRG::INSIDE, IntersectorFRG::OUTSIDE;
int IntersectorFRG::OUTSIDECOLOR;

//int debug_FRG_count = 0;
//----------------------------------------------------------------------------

/** Utility class to find and store bounding boxes for triangles */
class MyTriangle {
  int id;
  double x[3], w[3];
public:
  MyTriangle() {}
  MyTriangle(int i, Vec<Vec3D> &coord, int *nd) {
    id = i;
//    const double expansion_factor = 2.0;

    for(int j = 0; j <3; ++j) {
      x[j] = std::min(std::min(coord[nd[0]][j], coord[nd[1]][j]), coord[nd[2]][j]);
      w[j] = std::max(std::max(coord[nd[0]][j], coord[nd[1]][j]), coord[nd[2]][j])-x[j];
    }
/*
    double extra = 0.5*expansion_factor*std::max(w[0], std::max(w[1],w[2]));
    for(int j=0; j<3; j++) {
      x[j] -= extra;
      w[j] += extra;
    }
*/
  }
  double val(int i) const { return x[i]; }
  double width(int i) const { return w[i]; }
  int trId() const { return id; }
};

//----------------------------------------------------------------------------

/** Utility class to find the signed distance to the surface of the structure */
class ClosestTriangle {
  int (*triNodes)[3];
  Vec3D *structX;
  Vec3D *structNorm;
  set<int> *node2node;
  set<int> *node2elem;

public: //for debug only
//protected:
  bool fail;
  bool isFirst;
  int bestTrId;
  int n1, n2; //!< if both ns are non negative, the best point is on an edge
  Vec3D n;
  double minDist; //!< Signed distance to the surface
  Vec3D x;
  map<int,int> nd2tri; //needed for the vertex case (only if isConsistent = false)
  map<int,int> nd2tester;
  double dist2othernode;

  static const int maxNPairs = 100;

  bool isConsistent, isPositive, hasBeenFixed;
  int nPairs;
  int periTri[maxNPairs];

  bool checkEdge(int trId, int p1, int p2, int p3, double trDist);
  void checkVertex(int vn, int trId, double trDist);
  int registerNodes(int ip1, int trId, int& repeated1, int& repeated2);
  double project(Vec3D x0, int tria, double& xi1, double& xi2) const;
  double edgeProject(Vec3D x0, int n1, int n2, double &alpha) const;
  double getSignedVertexDistance() const;
  double findSignedVertexDistance();
  int findTesterStatus(Vec3D) const;
  void checkEdgeForTester(Vec3D xt, int trId, int ip1, int ip2, int p3, double trDist, int &nn1, int &nn2,
                          double &mindist, int &myMode, int &bestTriangle, const double eps) const;

	double piercing(Vec3D x0, int tria, double xi[3]);

public:
  ClosestTriangle(int (*triNodes)[3], Vec3D *structX, Vec3D *sN, set<int> *n2n, set<int> *n2e);
  void start(Vec3D x);
  void checkTriangle(int trId);
  double signedDistance() {
    if(n1 < 0 || n2 >= 0) return minDist;
    else
      return findSignedVertexDistance();
      //return getSignedVertexDistance();
  }

  int bestTriangle() const { return bestTrId; }

  int mode;
};

//----------------------------------------------------------------------------

ClosestTriangle::ClosestTriangle(int (*nd)[3], Vec3D *sX, Vec3D *sN, set<int> *n2n, set<int> *n2e) {
  fail = false;
  triNodes = nd;
  structX = sX;
  structNorm = sN;
  node2node = n2n;
  node2elem = n2e;
}

//----------------------------------------------------------------------------

void
ClosestTriangle::start(Vec3D xp) {
  x = xp;
  isFirst = true;
  bestTrId = -1;
  n1 = n2 = -1;
  mode = -1;
  nPairs = 0;
  minDist = 1.0e10;
  nd2tri.clear();
  nd2tester.clear();
  dist2othernode = 1.0e10;
}
//----------------------------------------------------------------------------
double ClosestTriangle::edgeProject(Vec3D x0, int n1, int n2, double &alpha) const
{
  Vec3D xA =   structX[n1];
  Vec3D xB =   structX[n2];
  Vec3D AB= xB-xA;
  Vec3D AX = x0-xA;
  alpha = AB*AX/(AB*AB);
  Vec3D P = xA + alpha*AB;
  return (P-x0).norm();
}
//----------------------------------------------------------------------------

double ClosestTriangle::project(Vec3D x0, int tria, double& xi1, double& xi2) const
{
  int iA = triNodes[tria][0];
  int iB = triNodes[tria][1];
  int iC = triNodes[tria][2];
  Vec3D xA = structX[iA];
  Vec3D xB = structX[iB];
  Vec3D xC = structX[iC];

  Vec3D ABC = (xB-xA)^(xC-xA);
  double areaABC = ABC.norm();
  Vec3D dir = 1.0/areaABC*ABC;

  //calculate the projection.
  double dist = (x0-xA)*dir;
  Vec3D xp = x0 - dist*dir;

  //calculate barycentric coords.
  double areaPBC = (((xB-xp)^(xC-xp))*dir);
  double areaPCA = (((xC-xp)^(xA-xp))*dir);
  xi1 = areaPBC/areaABC;
  xi2 = areaPCA/areaABC;

  return dist;
}

//----------------------------------------------------------------------------

void
ClosestTriangle::checkTriangle(int trId) {
  double dist;
  int *nd = triNodes[trId];
  double xi1, xi2;
  dist = project(x, trId, xi1, xi2);
  double xi3 = 1-xi1-xi2;
  const double eps = 0;

  if(xi1 >= -eps && xi2 >= -eps && xi3 >= -eps) {
    if(isFirst || std::abs(minDist) > std::abs(dist)) {
      isFirst = false;
      minDist = dist;
      bestTrId = trId;
      n1 = n2 = -1;
      mode = 0;
    }
  }else {
    if(xi1 < -eps)
      checkEdge(trId, nd[1], nd[2], nd[0], dist);
    if(xi2 < -eps)
      checkEdge(trId, nd[0], nd[2], nd[1], dist);
    if(xi3 < -eps)
      checkEdge(trId, nd[0], nd[1], nd[2], dist);
  }
}

//----------------------------------------------------------------------------

void ClosestTriangle::checkVertex(int ip1, int trId, double trDist) {

  // If this node is already our best solution
  if(n1 == ip1 && n2 < 0) {
    if(nPairs == maxNPairs) {
      std::cerr << "WARNING: On the embedded structure surface, detected too many (>" << maxNPairs << ") peripheral triangles to a node!" << std::endl;
      throw "Too many peripheral triangles";
    }
    periTri[nPairs++] = trId;
    int repeated1, repeated2;
    int repeated_nodes = registerNodes(ip1,trId, repeated1, repeated2);
    bool thisSign = (trDist >= 0);

    if((thisSign != isPositive || !isConsistent)/* && !hasBeenFixed*/) { // need to figure out the correct sign...
      isConsistent = false;
      if(repeated1 != -1) {// this triangle and another traversed triangle share an edge that has this damn node.
        //Step 1. Determine if it is a "hill" or a "dent". 
        int repeated = repeated1;
        double dist2node = (structX[triNodes[trId][repeated]]-x).norm();
        if(dist2node<dist2othernode) {
          dist2othernode = dist2node;
          double xi1_temp, xi2_temp;
          double dist2 = project(structX[nd2tester[triNodes[trId][repeated]]], trId, xi1_temp, xi2_temp);
          int type = (dist2>0) ? 1 : 2; //1: dent, 2: hill
          //Step 2. Check the angles between (ip1->x) and the normal of the two triangles
          Vec3D ip1_x = x - structX[ip1];
          double alpha = ip1_x*structNorm[trId];
          double beta  = ip1_x*structNorm[nd2tri[triNodes[trId][repeated]]];
          if(type==1) 
            isPositive = (alpha>=0&&beta>=0) ? true : false;
          else
            isPositive = (alpha<0&&beta<0) ? false : true;
         
//          hasBeenFixed = true; //the sign is determined.
        }
      }
      if(repeated2 != -1) {// this triangle and another traversed triangle share an edge that has this damn node.
        //Step 1. Determine if it is a "hill" or a "dent". 
        int repeated = repeated2;
        double dist2node = (structX[triNodes[trId][repeated]]-x).norm();
        if(dist2node<dist2othernode) {
          dist2othernode = dist2node;
          double xi1_temp, xi2_temp;
          double dist2 = project(structX[nd2tester[triNodes[trId][repeated]]], trId, xi1_temp, xi2_temp);
          int type = (dist2>0) ? 1 : 2; //1: dent, 2: hill
          //Step 2. Check the angles between (ip1->x) and the normal of the two triangles
          Vec3D ip1_x = x - structX[ip1];
          double alpha = ip1_x*structNorm[trId];
          double beta  = ip1_x*structNorm[nd2tri[triNodes[trId][repeated]]];
          if(type==1) 
            isPositive = (alpha>=0&&beta>=0) ? true : false;
          else
            isPositive = (alpha<0&&beta<0) ? false : true;
         
//          hasBeenFixed = true; //the sign is determined.
        }
      }
    }
    return;
  }
  // Compute the distance to the node. Determining the sign of the distance
  // is delayed until the very end
  double dist = (x-structX[ip1]).norm();
  if(isFirst || dist < std::abs(minDist)) {
    isFirst = false;
    n1 = ip1;
    n2 = -1;
    nPairs = 1;
    minDist = dist;
    dist2othernode = 1.0e10; //distance to the other node on the edge that is shared with another triangle
    periTri[0] = trId;
    bestTrId = trId;
    nd2tri.clear();
    nd2tester.clear();
    int repeated1, repeated2;
    registerNodes(ip1,trId, repeated1, repeated2);
    mode = 2;
    hasBeenFixed = false;
    isConsistent = true;
    isPositive = trDist >= 0;
  }
}

//----------------------------------------------------------------------------

int ClosestTriangle::registerNodes(int ip1, int trId, int& repeated1, int& repeated2)
{
    map<int,int>::iterator it;
      
    unsigned int repeated_node = 0;
    repeated1 = repeated2 = -1;

    if(ip1==triNodes[trId][0]) {
        it = nd2tri.find(triNodes[trId][1]);
        if(it==nd2tri.end()) {
          nd2tri[triNodes[trId][1]] = trId;
          nd2tester[triNodes[trId][1]] = triNodes[trId][2];
        } else if(it->second!=trId) {
          repeated_node |= 1 << 1;
          repeated1 = 1;
        }
        it = nd2tri.find(triNodes[trId][2]);
        if(it==nd2tri.end()) {
           nd2tri[triNodes[trId][2]] = trId;
           nd2tester[triNodes[trId][2]] = triNodes[trId][1];
        } else if(it->second!=trId) {
          repeated_node |= 1 << 2;
          repeated2 = 2;
        }
    }
    else if(ip1==triNodes[trId][1]) {
        it = nd2tri.find(triNodes[trId][2]);
        if(it==nd2tri.end()) {
          nd2tri[triNodes[trId][2]] = trId;
          nd2tester[triNodes[trId][2]] = triNodes[trId][0];
        } else if(it->second!=trId) {
          repeated_node |= 1 << 2;
          repeated1 = 2;
        }
        it = nd2tri.find(triNodes[trId][0]);
        if(it==nd2tri.end()) {
           nd2tri[triNodes[trId][0]] = trId;
           nd2tester[triNodes[trId][0]] = triNodes[trId][2];
        } else if(it->second!=trId) {
          repeated_node |= 1 << 0;
          repeated2 = 0;
        }
    }
    else if(ip1==triNodes[trId][2]) {
        it = nd2tri.find(triNodes[trId][0]);
        if(it==nd2tri.end()) {
          nd2tri[triNodes[trId][0]] = trId;
          nd2tester[triNodes[trId][0]] = triNodes[trId][1];
        } else if(it->second!=trId) {
          repeated_node |= 1 << 0;
          repeated1 = 0;
        }
        it = nd2tri.find(triNodes[trId][1]);
        if(it==nd2tri.end()) {
           nd2tri[triNodes[trId][1]] = trId;
           nd2tester[triNodes[trId][1]] = triNodes[trId][0];
        } else if(it->second!=trId) {
          repeated_node |= 1 << 1;
          repeated2 = 1;
        }
    } else
      fprintf(stderr,"ERROR (for debug): node %d doesn't belong to triangle %d!\n", ip1+1, trId+1);

    return repeated_node;
} 

//----------------------------------------------------------------------------

bool
ClosestTriangle::checkEdge(int trId, int ip1, int ip2, int p3, double trDist) {
  int p1, p2;
  if(ip1 < ip2) {
    p1 = ip1;
    p2 = ip2;
  } else {
    p1 = ip2;
    p2 = ip1;
  }

  if(n1 == p1 && n2 == p2) { //<! If we've already looked at this edge before
    // If we have already examined this edge and the sign opinions disagree, we need to make the right decision
    // Check if the p3 is in the positive or negative side of the other triangle
    // When signs disagree, the true distance has the opposite sign. because the triangle surface is the
    // end of point with the same sign as p3;
    if(trDist*minDist >= 0)
      return true;
    double xi1, xi2;
    double d2 = project(structX[p3], bestTrId, xi1, xi2);
    if(d2*minDist>0)
      minDist = -minDist;
    return true;
  } else {
    double dist, alpha;
    double sign = trDist >= 0 ? 1 : -1;
    dist = sign*edgeProject(x, p1, p2, alpha);

    int cn1 = (alpha > 1) ? p2 : p1;
    int cn2 = p2;
    if(alpha < 0 || alpha > 1) {
      checkVertex(cn1, trId, trDist);
      return true;
    }

    if(isFirst || std::abs(dist) < std::abs(minDist)) {
      isFirst = false;
      bestTrId = trId;
      n = structNorm[bestTrId];
      minDist = dist;
      n1 = cn1;
      n2 = cn2;
      mode = 1;
      return true;
    }
  }
  return false;
}

//----------------------------------------------------------------------------

double ClosestTriangle::getSignedVertexDistance() const {

  return isPositive ? minDist : -minDist;

/*  if(isConsistent)
    return isPositive ? minDist : -minDist;
  else if(hasBeenFixed)
    return isPositive ? minDist : -minDist;
  return minDist;
*/
}

//----------------------------------------------------------------------------

double ClosestTriangle::findSignedVertexDistance()
{
  set<int> &vertices = node2node[n1]; //vertices in the direct neighborhood of n1
  Vec3D    &xp       = structX[n1];   //coordinate of n1

  //step 1: find a test direction.
  Vec3D dir(0,0,0); //the test direction (normalized).
  double rr = 0.0; //distance in the test direction.

  double ry = 0.0;
  Vec3D npx = x - xp;
  npx = 1.0/npx.norm()*npx; 
  for(set<int>::iterator it=vertices.begin(); it!=vertices.end(); it++) {
    Vec3D xr = structX[*it]; 

    double d_xp_xr = (xp-xr).norm();
    if(d_xp_xr>=ry)
      ry = d_xp_xr;

    if(npx*(xp-xr)<1.0e-14) continue;
    double t = npx*(x-xr)/(npx*(xp-xr)); 
    Vec3D xt = xr + t*(xp-xr);
    double d_x_xt = (xt - x).norm();
    if(d_x_xt>=rr) {
      dir = xt - x;
      rr  = d_x_xt;
    }
  }
  if(rr<1.0e-14) {fprintf(stderr,"ERROR: (in IntersectorFRG) distance = %e.\n",rr);exit(-1);}
  dir *= 1.0/rr; //normalize dir
  rr = std::min(1.5*rr, std::max((x-xp).norm(), ry)); 
  rr *= 0.05;

  int nTrial = 50;
  for(int iTrial=0; iTrial<nTrial; iTrial++) {
    Vec3D x_trial;
    if(iTrial%2==1)
      x_trial = x - (double(iTrial+1.0)*rr)*dir;
    else
      x_trial = x + (double(iTrial+1.0)*rr)*dir;

    int sign = findTesterStatus(x_trial); 
    if(sign!=0) {
      minDist = (double)sign*std::abs(minDist);
//      fprintf(stderr,"NOTE: x = (%e, %e, %e), node = %d, x_trial = (%e, %e, %e). sign = %d, minDist = %e\n", x[0], x[1], x[2], n1+1, x_trial[0], x_trial[1], x_trial[2], sign, minDist);
      return minDist;
    }// else 
 //     fprintf(stderr,"NOTE: x = (%e, %e, %e), node = %d, x_trial = (%e, %e, %e). I = %d, sign = %d\n", x[0], x[1], x[2], n1+1, x_trial[0], x_trial[1], x_trial[2], iTrial, sign);
  }
  
  fprintf(stderr,"ERROR: (in IntersectorFRG) failed in determining node status! nTrial = %d.\n", nTrial);
  fail = true;
//  exit(-1);

  return 0.0;
}

//----------------------------------------------------------------------------

int ClosestTriangle::findTesterStatus(Vec3D xt) const
{
//  Vec3D xdebug(-4.348610e+00, -5.030928e+00, -6.636450e-02);
//  bool debug = (xt-xdebug).norm()<1.0e-5 ? true : false;

  set<int> &vertices = node2node[n1]; //vertices in the direct neighborhood of n1
  set<int> &elements = node2elem[n1]; //elements in the direct neighborhood of n1
  double mindist = 1.0e14;
  const double eps = 1.0e-14;
  int nn1,nn2;
  int myMode = -1;
  int bestTriangle = -1;
  nn1 = nn2 = -1;

  for(set<int>::iterator it=elements.begin(); it!=elements.end(); it++) {
    double xi[3];
    double dist = project(xt, *it, xi[0], xi[1]); 
    xi[2] = 1.0 - xi[0] - xi[1]; // project onto the plane determined by this triangle

    if((n1==triNodes[*it][0] || xi[0] >= -eps) && 
       (n1==triNodes[*it][1] || xi[1] >= -eps) &&
       (n1==triNodes[*it][2] || xi[2] >= -eps)) { 
      if(std::abs(mindist) >= std::abs(dist)) {
        mindist = dist;
        myMode = 0;
        nn1 = nn2 = -1;
        bestTriangle = *it; 
      }
    }
    else { 
      if(!(n1==triNodes[*it][0] || xi[0] >= -eps))
        checkEdgeForTester(xt, *it, triNodes[*it][1], triNodes[*it][2], triNodes[*it][0], dist, nn1, nn2, mindist, myMode, bestTriangle, eps);
      if(!(n1==triNodes[*it][1] || xi[1] >= -eps))
        checkEdgeForTester(xt, *it, triNodes[*it][2], triNodes[*it][0], triNodes[*it][1], dist, nn1, nn2, mindist, myMode, bestTriangle, eps);
      if(!(n1==triNodes[*it][2] || xi[2] >= -eps))
        checkEdgeForTester(xt, *it, triNodes[*it][0], triNodes[*it][1], triNodes[*it][2], dist, nn1, nn2, mindist, myMode, bestTriangle, eps);
    }
//    if(debug)
//      fprintf(stderr,"NOTE (debug): node = %d, trId = %d. | mode = %d, nn1 = %d, nn2 = %d, bestTriangle = %d, mindist = %e.\n",
//              n1+1, *it+1, myMode, nn1+1, nn2+1, bestTriangle+1, mindist);
  }  
//  if(debug)
//    fprintf(stderr,"NOTE (final): node = %d, mode = %d, nn1 = %d, nn2 = %d, bestTriangle = %d, mindist = %e.\n",
//            n1+1, myMode, nn1+1, nn2+1, bestTriangle+1, mindist);


  int sign;
  if(myMode<0) {
//    fprintf(stderr,"WARNING: (in IntersectorFRG) Mode = %d!\n", myMode);
//    fprintf(stderr,"-- x = (%e %e %e), tester = (%e, %e, %e).\n", x[0], x[1], x[2], xt[0], xt[1], xt[2]);
    sign = 0;
  } else
    sign = mindist>=0.0 ? 1 : -1;

  return sign;
}

//----------------------------------------------------------------------------

void ClosestTriangle::checkEdgeForTester(Vec3D xt, int trId, int ip1, int ip2, int p3, double trDist, int &nn1, int &nn2, 
                                         double &mindist, int &myMode, int &bestTriangle, const double eps) const
{
  int p1, p2;
  if(ip1 < ip2) {
    p1 = ip1;  p2 = ip2;
  } else {
    p1 = ip2;  p2 = ip1;
  }

  if(nn1 == p1 && nn2 == p2) { //this edge has been traversed.
    if(trDist*mindist >= 0)
      return;
    double xi1,xi2;
    double d2 = project(structX[p3], bestTriangle, xi1, xi2);
    if(d2*mindist>0)
      mindist = -mindist;
    return; 
  } else {
    double dist, alpha;
    double sign = trDist >= 0 ? 1 : -1;
    dist = sign*edgeProject(xt, p1, p2, alpha);
    if((alpha<-eps && p1==n1) || (alpha>1.0+eps && p2==n1))
      return;
    
    if(std::abs(dist) < std::abs(mindist)) {
      bestTriangle = trId;
      nn1 = p1;
      nn2 = p2;
      myMode = 1;
      mindist = dist;
    }
  } 
}

//----------------------------------------------------------------------------

double ClosestTriangle::piercing(Vec3D x0, int tria, double xi[3])
{

   const double tol = 1.0e-12;
	
	double dist;

	dist = std::abs( project(x0, tria, xi[0], xi[1]) );

	xi[2] = 1.0 - xi[0] - xi[1];

	if(xi[0] >= tol && xi[1] >= tol && xi[2] >= tol)
	{		
		return dist;
	} 
	else
	{
		dist = 1.0e16;

		for(int i=0; i<3; ++i)
		{
			if(xi[i] < tol)
			{
				int p1 = triNodes[tria][(i+1)%3];
				int p2 = triNodes[tria][(i+2)%3];
				
				double alpha;
				double d2l = std::abs(edgeProject(x0, p1, p2, alpha));

				if(alpha >= tol && alpha <= 1.0 + tol) 
				{
					if(dist > d2l) 
					{
						dist        = d2l; 
						xi[i]       = 0.0; 
						xi[(i+1)%3] = 1.0-alpha; 
						xi[(i+2)%3] = alpha;
					}
				} 
				else 
				{
					if(alpha < tol) 
					{
						double d2p = (x0 - structX[p1]).norm();

						if(dist > d2p) 
						{
							dist        = d2p; 
							xi[i]       = 0.0;
							xi[(i+2)%3] = 0.0; 
							xi[(i+1)%3] = 1.0;
						}
					}

					if(alpha > 1.0+tol) 
					{
						double d2p = (x0 - structX[p2]).norm();

						if(dist > d2p) 
						{
							dist        = d2p; 
							xi[i]       = 0.0;
							xi[(i+1)%3] = 0.0; 
							xi[(i+2)%3] = 1.0;
						}
					}
				}
			}
		}

		return dist;
	}

}

//----------------------------------------------------------------------------

DistIntersectorFRG::DistIntersectorFRG(IoData &iodata, 
													Communicator *comm, 
													int nNodes,	double *xyz, 
													int nElems, int (*abc)[3]) : DistLevelSetStructure(),iod(iodata)
{

  externalSI = (iod.embed.surrogateinterface == EmbeddedFramework::EXTERNAL) ? true : false;

  if(externalSI) comm->fprintf(stdout, " +++ Using external-based surrogate interface +++ \n");

  withViscousTerms = false;
  if(iod.eqs.type == EquationsData::NAVIER_STOKES)
  {
	  withViscousTerms = true;

	  if(iod.bc.wall.type == BcsWallData::ISOTHERMAL) 
	  {
		  isoThermalwall = true;
		  Twall = iod.bc.wall.temperature;
	  }
	  else
	  {
		  isoThermalwall = false;
		  Twall = -1.0; // dummy
	  }
  }

  this->numFluid = iod.eqs.numPhase;
  twoPhase = false;

  com = comm;
  com->fprintf(stderr,"- Using Intersector: FRG\n");

  //get embedded structure surface mesh and restart pos
  char *struct_mesh, *struct_restart_pos;
  int sp = strlen(iod.input.prefix) + 1;

  struct_mesh        = new char[sp + strlen(iod.input.embeddedSurface)];
  sprintf(struct_mesh,"%s%s", iod.input.prefix, iod.input.embeddedSurface);

  if (iod.input.restart_file_package[0] == 0) 
  {
    struct_restart_pos = new char[sp + strlen(iod.input.embeddedpositions)];
    if(iod.input.embeddedpositions[0] != 0)
      sprintf(struct_restart_pos,"%s%s", iod.input.prefix, iod.input.embeddedpositions);
    else //no restart position file provided
      strcpy(struct_restart_pos,""); 
  } 
  else 
  {
    struct_restart_pos = new char[256];
    char dummy[256], tmp[256];
    sprintf(tmp,"%s%s", iod.input.prefix, iod.input.restart_file_package);
    TsRestart::readRestartFileNames(tmp,
				    dummy,
				    dummy,
				    dummy,
				    dummy,
				    dummy,
				    dummy,
				    struct_restart_pos, comm);
				    
				    
  }
  
  //nodal or facet normal?
  interpolatedNormal = (iod.embed.structNormal==EmbeddedFramework::NODE_BASED) ?  true : false;

  with_sensitivity = false;

  //initialize the following to 0(NULL)
  globPhysInterface = 0;
  triNorms = 0;
  nodalNormal = 0;
  status = 0;
  status0 = 0;
  boxMin = 0;
  boxMax = 0;
  distance = 0;
  tId = 0;
  poly = 0;

  faceID = NULL;
  surfaceID = NULL;

  //Load files. Compute structure normals. Initialize PhysBAM Interface
  if(nNodes && xyz && nElems && abc)
    init(nNodes, xyz, nElems, abc, struct_restart_pos);
  else {
    double XScale = (iod.problem.mode==ProblemData::NON_DIMENSIONAL) ? 1.0 : iod.ref.rv.length;
    init(struct_mesh, struct_restart_pos, XScale);
  }
  setPorosity();
  makerotationownership();
  updatebc();

  if(numStElems<=0||numStNodes<=0) {
    fprintf(stderr,"ERROR: Found %d nodes and %d elements in the embedded surface!\n", numStNodes, numStElems);
    exit(-1);
  }

  node2node = new set<int>[numStNodes];
  node2elem = new set<int>[numStNodes];
  buildConnectivity();

  delete[] struct_mesh;
  delete[] struct_restart_pos;
}

//----------------------------------------------------------------------------

DistIntersectorFRG::~DistIntersectorFRG() 
{
  delete [] stElem;
  delete [] node2node;
  delete [] node2elem;
  if(Xs)          delete[] Xs;
  if(Xs0)         delete[] Xs0;
  if(Xs_n)        delete[] Xs_n;
  if(Xs_np1)      delete[] Xs_np1;
  if(Xsdot)       delete[] Xsdot;
  delete solidX;
  delete solidX0;
  delete solidXn;
  delete solidXnp1;
  if(faceID) delete[] faceID;
  if(surfaceID) delete[] surfaceID;
  if(porosity) delete[] porosity;
  if(rotOwn) delete[] rotOwn;
  if(status0)     delete   status0;
  if(triNorms)    delete[] triNorms;
  if(nodalNormal) delete[] nodalNormal;
  if(boxMax)      delete   boxMax;
  if(boxMin)      delete   boxMin;
  if(tId)         delete   tId;
  if(poly)        delete   poly;
  delete globPhysInterface;
  for(int i=0;i<domain->getNumLocSub();++i)
    {
      delete intersector[i];
    }
  delete[] intersector;

  if(dXdSb)    delete[] dXdSb;
  if(solidXdS) delete solidXdS;

}

//----------------------------------------------------------------------------

LevelSetStructure &
DistIntersectorFRG::operator()(int subNum) const {
  return *intersector[subNum];
}

//----------------------------------------------------------------------------

/** Intersector initialization method
*
* \param dataTree the data read from the input file for this intersector.
*/
void DistIntersectorFRG::init(char *solidSurface, char *restartSolidSurface, double XScale) {

  // Read data from the solid surface input file.
  FILE *topFile;
  topFile = fopen(solidSurface, "r");
  if (topFile == NULL) {com->fprintf(stderr, "Embedded structure surface mesh doesn't exist :(\n"); exit(1); }

  char line[MAXLINE],key1[MAXLINE],key2[MAXLINE];

  // load the nodes and initialize all node-based variables.

  // load solid nodes at t=0
  int num0 = 0; 
  int num1 = 0;

  double x1,x2,x3;
  int node1, node2, node3;

  int type_read = 0;
  int surfaceid = 0;

  std::list<int> indexList;
  std::list<Vec3D> nodeList;
  std::list<int> elemIdList;    
  std::list<int> elemList1;
  std::list<int> elemList2;
  std::list<int> elemList3;
  std::list<int> surfaceIDList;

  int maxIndex = 0, maxElem = -1;

  while (fgets(line,MAXLINE,topFile) != 0) {
    sscanf(line, "%s", key1);
    bool skip = false;
    if (strcmp(key1, "Nodes") == 0) {
      sscanf(line, "%*s %s", key2);
      skip = true;
      type_read = 1;
    }
    if (strcmp(key1, "Elements") == 0) {
      sscanf(line, "%*s %s", key2);
      skip = true;
      type_read = 2;
      int underscore_pos = -1;
      int k = 0;
      while((key2[k] != '\0') && (k<MAXLINE)) {
        if(key2[k] == '_') underscore_pos = k;
        k++;
      }
      if(underscore_pos > -1)
        sscanf(key2+(underscore_pos+1),"%d",&surfaceid); 
    }
    if (!skip) {
      if (type_read == 1) {
	sscanf(line, "%d %lf %lf %lf", &num1, &x1, &x2, &x3);
	if(num1<1) {com->fprintf(stderr,"ERROR: detected a node with index %d in the embedded surface file!\n",num1); exit(-1);}
        x1 /= XScale;
        x2 /= XScale;
        x3 /= XScale;
        indexList.push_back(num1);
        if(num1>maxIndex) maxIndex = num1;
        nodeList.push_back(Vec3D(x1,x2,x3));
      }
      if (type_read == 2) {
	sscanf(line,"%d %d %d %d %d", &num0, &num1, &node1, &node2, &node3);
        elemList1.push_back(node1-1);
        elemList2.push_back(node2-1);
        elemList3.push_back(node3-1);
	surfaceIDList.push_back(surfaceid);
      }
    }
  }

  numStNodes = nodeList.size();
  if(numStNodes != maxIndex) {
    com->fprintf(stderr,"ERROR: The node set of the embedded surface have gap(s). \n");
    com->fprintf(stderr,"       Detected max index = %d, number of nodes = %d\n", maxIndex, numStNodes);
    com->fprintf(stderr,"NOTE: Currently the node set of the embedded surface cannot have gaps. Moreover, the index must start from 1.\n");
    exit(-1);
  }

  numStElems = elemList1.size();

  // feed data to Xss. 
  Xs      = new Vec3D[numStNodes];
  Xs0     = new Vec3D[numStNodes];
  Xs_n    = new Vec3D[numStNodes];
  Xs_np1  = new Vec3D[numStNodes];
  Xsdot   = new Vec3D[numStNodes];
  solidX  = new Vec<Vec3D>(numStNodes, Xs);
  solidX0 = new Vec<Vec3D>(numStNodes, Xs0);
  solidXn = new Vec<Vec3D>(numStNodes, Xs_n);
  solidXnp1 = new Vec<Vec3D>(numStNodes, Xs_np1);

     dXdSb = new Vec3D[numStNodes];
  solidXdS = new Vec<Vec3D>(numStNodes, dXdSb);

  faceID = new int[numStElems];

  surfaceID = new int[numStNodes];

  std::list<Vec3D>::iterator it1;
  std::list<int>::iterator it2;
  std::list<int>::iterator it_1;
  std::list<int>::iterator it_2;
  std::list<int>::iterator it_3;
  std::list<int>::iterator it_4;

  it2 = indexList.begin();
  for (it1=nodeList.begin(); it1!=nodeList.end(); it1++) {
    Xs[(*it2)-1] = *it1;
    it2++;
  }

  for (int k=0; k<numStNodes; k++) {
    Xs0[k]    = Xs[k];
    Xs_n[k]   = Xs[k];
    Xs_np1[k] = Xs[k];
    Xsdot[k]  = Vec3D(0.0, 0.0, 0.0);
    surfaceID[k] = 0;
    dXdSb[k]     = 0;
  }

  stElem = new int[numStElems][3];
  
  it_1 = elemList1.begin();
  it_2 = elemList2.begin();
  it_3 = elemList3.begin();
  it_4 = surfaceIDList.begin();
  for (int i=0; i<numStElems; i++) {
    stElem[i][0] = *it_1;
    stElem[i][1] = *it_2;
    stElem[i][2] = *it_3;
    faceID[i] = *it_4;
    surfaceID[*it_1] = *it_4;
    surfaceID[*it_2] = *it_4;
    surfaceID[*it_3] = *it_4;
    it_1++;
    it_2++;
    it_3++;
    it_4++;
  } 

  fclose(topFile);

  int nInputs;
  char c1[200], c2[200], c3[200];
  // load solid nodes at restart time.
  if (restartSolidSurface[0] != 0) {
    FILE* resTopFile = fopen(restartSolidSurface, "r");
    if(resTopFile==NULL) {com->fprintf(stderr, "restart topFile doesn't exist.\n"); exit(1);}
    int ndMax = 0, ndMax2 = 0;
    std::list<std::pair<int,Vec3D> > nodeList2;
    std::list<std::pair<int,Vec3D> >::iterator it2;

    while(1) {
      nInputs = fscanf(resTopFile,"%s", c1);
      if(nInputs!=1) break;    
      char *endptr;
      num1 = strtol(c1, &endptr, 10);
      if(endptr == c1) break;

      int toto = fscanf(resTopFile,"%lf %lf %lf\n", &x1, &x2, &x3);
      nodeList2.push_back(std::pair<int,Vec3D>(num1,Vec3D(x1,x2,x3)));
      ndMax = std::max(num1, ndMax);
    }
    if (ndMax!=numStNodes) {
      com->fprintf(stderr,"ERROR: number of nodes in restart top-file is wrong.\n");
      exit(1);
    }

    for (int k=0; k<numStNodes; k++)
      Xs[k] = Vec3D(0,0,0);
    
    for (it2=nodeList2.begin(); it2!=nodeList2.end(); it2++)
      Xs[it2->first-1] = it2->second;

    for (int k=0; k<numStNodes; k++) {
      Xs_n[k]         = Xs[k];
      Xs_np1[k]    = Xs[k];
    }
    fclose(resTopFile);
  }

  // Verify (1)triangulated surface is closed (2) normal's of all triangles point outward.
  com->fprintf(stderr,"- IntersectorFRG: Checking the embedded structure surface...   ");
  if (checkTriangulatedSurface()) 
    com->fprintf(stderr,"Ok.\n");
  else {
    com->fprintf(stderr,"\n");
    exit(-1); 
  }

//  getBoundingBox();
  initializePhysBAM();
}

//----------------------------------------------------------------------------
/** Intersector initialization method
*
* \param dataTree the data read from the input file for this intersector.
*/
void DistIntersectorFRG::init(int nNodes, double *xyz, int nElems, int (*abc)[3], char *restartSolidSurface) {

  // node set
  numStNodes = nNodes;
  // feed data to Xss. 
  Xs      = new Vec3D[numStNodes];
  Xs0     = new Vec3D[numStNodes];
  Xs_n    = new Vec3D[numStNodes];
  Xs_np1  = new Vec3D[numStNodes];
  Xsdot   = new Vec3D[numStNodes];
  solidX  = new Vec<Vec3D>(numStNodes, Xs);
  solidX0 = new Vec<Vec3D>(numStNodes, Xs0);
  solidXn = new Vec<Vec3D>(numStNodes, Xs_n);
  solidXnp1 = new Vec<Vec3D>(numStNodes, Xs_np1);

     dXdSb = new Vec3D[numStNodes];
  solidXdS = new Vec<Vec3D>(numStNodes, dXdSb);

  for (int k=0; k<numStNodes; k++) {
    Xs[k]     = Vec3D(xyz[3*k], xyz[3*k+1], xyz[3*k+2]);
    Xs0[k]    = Xs[k];
    Xs_n[k]   = Xs[k];
    Xs_np1[k] = Xs[k];
    Xsdot[k]  = Vec3D(0.0, 0.0, 0.0);
    dXdSb[k]   = 0;
  }

  // elem set
  numStElems = nElems;
  stElem = new int[numStElems][3];
  for (int i=0; i<numStElems; i++)
    for (int k=0; k<3; k++)
      stElem[i][k] = abc[i][k];

  // load solid nodes at restart time.
  if (restartSolidSurface[0] != 0) {
    FILE* resTopFile = fopen(restartSolidSurface, "r");
    if(resTopFile==NULL) {com->fprintf(stderr, "restart topFile doesn't exist.\n"); exit(1);}
    int ndMax2 = 0;
    std::list<std::pair<int,Vec3D> > nodeList2;
    std::list<std::pair<int,Vec3D> >::iterator it2;

    int ndMax = 0;
    while(1) {
      int nInputs, num1;
      double x1, x2, x3;
      char *endptr, c1[200];

      nInputs = fscanf(resTopFile,"%s", c1);
      if(nInputs!=1) break;
      num1 = strtol(c1, &endptr, 10);
      if(endptr == c1) break;

      int toto = fscanf(resTopFile,"%lf %lf %lf\n", &x1, &x2, &x3);
      nodeList2.push_back(std::pair<int,Vec3D>(num1,Vec3D(x1,x2,x3)));
      ndMax = std::max(num1, ndMax);
    }

    if (ndMax!=numStNodes) {
      com->fprintf(stderr,"ERROR: number of nodes in restart top-file is wrong.\n");
      com->fprintf(stderr,"ndMax = %d; numStNodes = %d\n", ndMax, numStNodes);
      exit(1);
    }

    for (int k=0; k<numStNodes; k++)
      Xs[k] = Vec3D(0,0,0);

    for (it2=nodeList2.begin(); it2!=nodeList2.end(); it2++)
      Xs[it2->first-1] = it2->second;

    for (int k=0; k<numStNodes; k++) {
      Xs_n[k]   = Xs[k];
      Xs_np1[k] = Xs[k];
    }
    fclose(resTopFile);
  }

  // Verify (1)triangulated surface is closed (2) normal's of all triangles point outward.
  com->fprintf(stderr,"- IntersectorFRG: Checking the embedded structure surface...   ");
  if (checkTriangulatedSurface())
    com->fprintf(stderr,"Ok.\n");
  else {
    com->fprintf(stderr,"\n");
    exit(-1);
  }

//  getBoundingBox();
  initializePhysBAM();
}

//----------------------------------------------------------------------------
void DistIntersectorFRG::setPorosity() {
  map<int,SurfaceData *> &surfaceMap = iod.surfaces.surfaceMap.dataMap;
  map<int,BoundaryData *> &bcMap = iod.bc.bcMap.dataMap;

  porosity = new double[numStElems];
  for(int i=0; i<numStElems; i++) {
    porosity[i] = 0.;
  }

  if(faceID) {
    for(int i=0; i<numStElems; i++) {
      map<int,SurfaceData*>::iterator it = surfaceMap.find(faceID[i]);
      if (it != surfaceMap.end()) {
        map<int,BoundaryData *>::iterator it2 = bcMap.find(it->second->bcID);
        if(it2 != bcMap.end()) { // the bc data have been defined
          if(it2->second->type == BoundaryData::POROUSWALL ) {
            porosity[i] = it2->second->porosity;
          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------------
void DistIntersectorFRG::setdXdSb(int N, double* dxdS, double* dydS, double* dzdS){
  
  if(N != numStNodes) {
    fprintf(stderr, "Error num nodes: %d %d\n", N,  numStNodes);
    exit(-1);
  }

  with_sensitivity = true;
  //dXdSb = new Vec3D[numStNodes];

  double Normx=0.0, Normy=0.0, Normz=0.0;

  for(int i=0; i<N; ++i){
    dXdSb[i][0] = dxdS[i];
    dXdSb[i][1] = dydS[i];
    dXdSb[i][2] = dzdS[i];

    Normx += fabs(dxdS[i]);
    Normy += fabs(dydS[i]);
    Normz += fabs(dzdS[i]);
  }

  if(Normx==0 && Normy==0 && Normz==0)
    com->fprintf(stderr, "!!! WARNING: No surface Sensitivity Perturbation\n\n");

}

//----------------------------------------------------------------------------
void DistIntersectorFRG::makerotationownership() {
  map<int,SurfaceData *> &surfaceMap = iod.surfaces.surfaceMap.dataMap;
  map<int,RotationData*> &rotationMap= iod.rotations.rotationMap.dataMap;
  map<int,SurfaceData *>::iterator it = surfaceMap.begin();
  rotOwn = 0;

  int numRotSurfs = 0;
  int numTransWalls = 0;
  while(it != surfaceMap.end()) {
    map<int,RotationData*>::iterator it1 = rotationMap.find(it->second->rotationID);
    if(it1!=rotationMap.end()) {
      if(it1->second->infRadius) {
        numTransWalls++;
        com->fprintf(stderr," ... surface %2d is ``translating''\n",it->first, it1->first);
        com->fprintf(stderr,"     -> uniform velocity V = %3.2e in direction %3.2e %3.2e %3.2e\n",
                     it1->second->omega, it1->second->nx,it1->second->ny,it1->second->nz);
      } else {
        numRotSurfs++;
        com->fprintf(stderr," ... surface %2d is ``rotating'' using rotation data %2d\n",it->first, it1->first);
        com->fprintf(stderr,"     -> omega = %3.2e, rotation axis = %3.2e %3.2e %3.2e\n",
                     it1->second->omega, it1->second->nx,it1->second->ny,it1->second->nz);
      }
    }
    it++;
  }

  if(numRotSurfs || numTransWalls) {
    rotOwn = new int[numStNodes];

    for (int k=0; k<numStNodes; k++) {
      rotOwn[k] = -1;

      map<int,SurfaceData *>::iterator it = surfaceMap.find(surfaceID[k]);
      if(it != surfaceMap.end()) {
         int rotID = it->second->rotationID; // = -1 if not defined in input file
         if ((rotOwn[k] != -1 && rotOwn[k] != rotID) ||
             (rotOwn[k] != -1 && rotOwn[k] != rotID) ||
             (rotOwn[k] != -1 && rotOwn[k] != rotID))  {
           com->fprintf(stderr, " ... WARNING: Embedded Node %d associated to more than 1 Rotation ID\n",k);
         }
         rotOwn[k] = rotID;
         rotOwn[k] = rotID;
         rotOwn[k] = rotID;
      }
    }
  }
}
//----------------------------------------------------------------------------
void DistIntersectorFRG::updatebc() {
  map<int,RotationData*> &rotationMap= iod.rotations.rotationMap.dataMap;
  if (rotOwn)  { 
    for (int k=0; k<numStNodes; k++) {
      if (rotOwn[k]>=0) {    // node belongs to a (potential) "rotating" surface
        map<int,RotationData *>::iterator it =  rotationMap.find(rotOwn[k]);
        if(it != rotationMap.end()) { // the rotation data have been defined
	  if(it->second->infRadius == RotationData::TRUE) {
            double vel = it->second->omega / iod.ref.rv.velocity;
	    Xsdot[k][0] = vel*it->second->nx + Xsdot[k][0];
	    Xsdot[k][1] = vel*it->second->ny + Xsdot[k][1];
	    Xsdot[k][2] = vel*it->second->nz + Xsdot[k][2];
	  } 
          else {
            double XScale = (iod.problem.mode==ProblemData::NON_DIMENSIONAL) ? 1.0 : iod.ref.rv.length;
	    double tref = iod.ref.rv.time;
	    double ox = tref*it->second->omega*it->second->nx;
	    double oy = tref*it->second->omega*it->second->ny;
	    double oz = tref*it->second->omega*it->second->nz;
	    double xd = Xs[k][0] - it->second->x0 / XScale;
	    double yd = Xs[k][1] - it->second->y0 / XScale;
	    double zd = Xs[k][2] - it->second->z0 / XScale;
	    Xsdot[k][0] = oy*zd-oz*yd + Xsdot[k][0];
	    Xsdot[k][1] = oz*xd-ox*zd + Xsdot[k][1];
	    Xsdot[k][2] = ox*yd-oy*xd + Xsdot[k][2];
	  }
	} 
        else  { // no rotation data -> use velocity from mesh motion if any
          Xsdot[k][0] = Xsdot[k][0];
          Xsdot[k][1] = Xsdot[k][1];
          Xsdot[k][2] = Xsdot[k][2];
	}
      }
      else  { // no rotation data -> use velocity from mesh motion if any
        Xsdot[k][0] = Xsdot[k][0];
        Xsdot[k][1] = Xsdot[k][1];
        Xsdot[k][2] = Xsdot[k][2];
      }
    }
  }
}
//----------------------------------------------------------------------------
/*
void DistIntersectorFRG::getBoundingBox() {
  xMin = xMax = Xs[0][0];
  yMin = yMax = Xs[0][1];
  zMin = zMax = Xs[0][2];
  for(int i = 1; i < numStNodes; ++i) {
    xMin = std::min(xMin, Xs[i][0]);
    xMax = std::max(xMax, Xs[i][0]);
    yMin = std::min(yMin, Xs[i][1]);
    yMax = std::max(yMax, Xs[i][1]);
    zMin = std::min(zMin, Xs[i][2]);
    zMax = std::max(zMax, Xs[i][2]);
  }
}
*/
//----------------------------------------------------------------------------

void
DistIntersectorFRG::initializePhysBAM() { //NOTE: In PhysBAM array index starts from 1 instead of 0
  // Initialize the Particles list
  PhysBAM::GEOMETRY_PARTICLES<PhysBAM::VECTOR<double,3> > *physbam_solids_particle=new PhysBAM::GEOMETRY_PARTICLES<PhysBAM::VECTOR<double,3> >();
  physbam_solids_particle->array_collection.Resize(numStNodes);
  for (int i=0; i<numStNodes; i++) physbam_solids_particle->X(i+1) = PhysBAM::VECTOR<double,3>(Xs[i][0],Xs[i][1], Xs[i][2]);
  
  // Initialize the Triangle list.
  PhysBAM::ARRAY<PhysBAM::VECTOR<int,3> > physbam_stElem(numStElems);
  for (int i=0; i<numStElems; i++) physbam_stElem(i+1) = PhysBAM::VECTOR<int,3>(stElem[i][0]+1,stElem[i][1]+1,stElem[i][2]+1);

  // Initialize the mesh.
//  com->fprintf(stderr,"Initializing the Mesh with %d particles and %d triangles\n",physbam_solids_particle->array_collection.Size(),physbam_stElem.Size());
  PhysBAM::TRIANGLE_MESH *mesh = new PhysBAM::TRIANGLE_MESH(numStNodes,physbam_stElem);
  mesh->Initialize_Adjacent_Elements();mesh->Set_Number_Nodes(numStNodes);

  // Construct TRIANGULATED_SURFACE.
  if(globPhysInterface) delete globPhysInterface;
  globPhysInterface = new PhysBAMInterface<double>(*mesh,*physbam_solids_particle);

  with_sensitivity = false;
  solidXdS = 0;

}

//----------------------------------------------------------------------------

EdgePair DistIntersectorFRG::makeEdgePair(int node1, int node2, int triangleNumber) {
if(node1 < node2)
 return EdgePair(iipair(node1, node2), ibpair(triangleNumber, true));
else
 return EdgePair(iipair(node2, node1), ibpair(triangleNumber, false));
}

//----------------------------------------------------------------------------

bool DistIntersectorFRG::checkTriangulatedSurface()
{
  map<iipair, ibpair> edgeMap;

  for (int iTriangle=0; iTriangle<numStElems; iTriangle++) {
    int from1, to1, from2, to2, from3, to3;
    bool found1, found2, found3;
    from1 = stElem[iTriangle][0];  to1 = stElem[iTriangle][1];  found1 = false;
    from2 = stElem[iTriangle][1];  to2 = stElem[iTriangle][2];  found2 = false;
    from3 = stElem[iTriangle][2];  to3 = stElem[iTriangle][0];  found3 = false;

    EdgePair ep[3];
    ep[0] = makeEdgePair(from1, to1, iTriangle);
    ep[1] = makeEdgePair(from2, to2, iTriangle);
    ep[2] = makeEdgePair(from3, to3, iTriangle);

    for(int i=0; i < 3; ++i) {
      map<iipair, ibpair>::iterator it = edgeMap.find(ep[i].first);
      if(it != edgeMap.end()) { // we found this edge
         if(it->second.second == ep[i].second.second)
           {com->fprintf(stderr,"ERROR: surface is not closed or a triangle orientation problem. exit.\n"); return false;}
      } else
        edgeMap[ep[i].first] = ep[i].second;
    }
  }
  return true;
}

//----------------------------------------------------------------------------

void DistIntersectorFRG::buildConnectivity() {
  for(int i=0; i<numStElems; i++) 
    for(int j=0; j<3; j++) {// assume triangle elements
      node2node[stElem[i][j]].insert(stElem[i][(j+1)%3]);
      node2node[stElem[i][j]].insert(stElem[i][(j+2)%3]);
      node2elem[stElem[i][j]].insert(i);
    }
}

//----------------------------------------------------------------------------

void DistIntersectorFRG::buildSolidNormals() {
  if(!triNorms) triNorms = new Vec3D[numStElems];
  if(interpolatedNormal) {
    if(!nodalNormal)
      nodalNormal = new Vec3D[numStNodes];
    for(int i=0; i<numStNodes; i++)
      nodalNormal[i] = 0.0;
  }

  // Also look to determine a point inside the solid but away from the structure.
  for (int iTriangle=0; iTriangle<numStElems; iTriangle++) {
    int n1 = stElem[iTriangle][0];
    int n2 = stElem[iTriangle][1];
    int n3 = stElem[iTriangle][2];
    double x1 = Xs[n1][0];
    double y1 = Xs[n1][1];
    double z1 = Xs[n1][2];
    double dx2 = Xs[n2][0]-x1;
    double dy2 = Xs[n2][1]-y1;
    double dz2 = Xs[n2][2]-z1;
    double dx3 = Xs[n3][0]-x1;
    double dy3 = Xs[n3][1]-y1;
    double dz3 = Xs[n3][2]-z1;

    // calculate the normal.
    triNorms[iTriangle] = Vec3D(dx2, dy2, dz2)^Vec3D(dx3,dy3,dz3);
    
    if(interpolatedNormal){ // compute nodal normal (weighted by area)
      nodalNormal[n1] += triNorms[iTriangle];
      nodalNormal[n2] += triNorms[iTriangle];
      nodalNormal[n3] += triNorms[iTriangle];
    }

    double nrm = triNorms[iTriangle].norm();
    // normalize the normal.
    if(nrm > 0.0)
      triNorms[iTriangle] /= nrm;
    else
      fprintf(stderr,"ERROR: area of Triangle %d is %e.\n", iTriangle+1, nrm);
  }

  if(interpolatedNormal) //normalize nodal normals.
    for(int i=0; i<numStNodes; i++)
      nodalNormal[i] /= nodalNormal[i].norm();
}

//----------------------------------------------------------------------------

void DistIntersectorFRG::expandScope()
{
  SubDomain **subs = domain->getSubDomain();

  // 1. setup communication pattern
  int numLocSub = domain->getNumLocSub();
  SubDTopo *subTopo = domain->getSubTopo();
  CommPattern<int> trader(subTopo, com, CommPattern<int>::CopyOnSend, CommPattern<int>::NonSym);
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++) {
    int *sndChannel = subs[iSub]->getSndChannel();
    for(int iNei=0; iNei<subs[iSub]->getNumNeighb(); iNei++)
      trader.setLen(sndChannel[iNei], 1+intersector[iSub]->package[iNei].size());     
  }    
  trader.finalize();

  // 2. send packages to neighbour subdomains.
  set<int>::iterator it;
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++) {
    int *sndChannel = subs[iSub]->getSndChannel();
    for(int iNei=0; iNei<subs[iSub]->getNumNeighb(); iNei++) {
      SubRecInfo<int> sInfo = trader.getSendBuffer(sndChannel[iNei]);
      int *buffer = reinterpret_cast<int*>(sInfo.data);
      buffer[0] = intersector[iSub]->package[iNei].size(); 
      int count = 0;
      for(it=intersector[iSub]->package[iNei].begin(); it!= intersector[iSub]->package[iNei].end(); it++)
        buffer[++count] = *it;
    }
  }

  // 3. exchange
  trader.exchange();

  // 4. expand scope
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++) {

/*    char ch[20] = "debug_cpuX";
    ch[9] = 'A' + subs[iSub]->getGlobSubNum();
    FILE* myDebug = fopen(ch,"w");*/

    int *rcvChannel = subs[iSub]->getRcvChannel();
    for(int iNei=0; iNei<subs[iSub]->getNumNeighb(); iNei++) {
      SubRecInfo<int> sInfo = trader.recData(rcvChannel[iNei]);
      int *buffer = reinterpret_cast<int*>(sInfo.data);
      for(int j=0; j<buffer[0]; j++)
        intersector[iSub]->scope.insert(buffer[j+1]);
    }
/*
    fprintf(myDebug,"Scope (%d) of Sub %d is:\n", intersector[iSub]->scope.size(), subs[iSub]->getGlobSubNum());
    for(set<int>::iterator itt=intersector[iSub]->scope.begin(); itt!=intersector[iSub]->scope.end(); itt++)
      fprintf(myDebug,"%d\n", *itt+1);
    fprintf(myDebug,"\n");
    fclose(myDebug);*/
  }
}

//----------------------------------------------------------------------------

/** compute the intersections, node statuses and normals for the initial geometry */
void
DistIntersectorFRG::initialize(Domain *d, DistSVec<double,3> &X, DistSVec<double,3> &Xn, IoData &iod, DistVec<int> *point_based_id, DistVec<int>* oldStatus) {
  if(this->numFluid<1) {
    fprintf(stderr,"ERROR: number of fluid = %d!\n", this->numFluid);
    exit(-1);
  }
  this->X = &X;
  this->Xn = &Xn;
  domain = d;
  numLocSub = d->getNumLocSub();
  intersector = new IntersectorFRG*[numLocSub];

  status0 = new DistVec<int>(domain->getNodeDistInfo());
  boxMin = new DistSVec<double,3>(domain->getNodeDistInfo());
  boxMax = new DistSVec<double,3>(domain->getNodeDistInfo());
  tId = new DistVec<int>(domain->getNodeDistInfo());
  status = new DistVec<int>(domain->getNodeDistInfo());
  distance = new DistVec<double>(domain->getNodeDistInfo());
  is_swept = new DistVec<bool>(domain->getNodeDistInfo());
  is_active = new DistVec<bool>(domain->getNodeDistInfo());
  is_occluded = new DistVec<bool>(domain->getNodeDistInfo());
  edge_intersects = new DistVec<bool>(domain->getEdgeDistInfo());

  edge_SI  = new DistVec<bool>(domain->getEdgeDistInfo());

     xi_SI = new DistVec<double>(domain->getEdgeDistInfo());
    eta_SI = new DistVec<double>(domain->getEdgeDistInfo());
  TriID_SI = new DistVec<int>(domain->getEdgeDistInfo());
  nWall_SI = new DistVec<Vec3D>(domain->getEdgeDistInfo());

     xi_node = new DistVec<double>(domain->getNodeDistInfo());
    eta_node = new DistVec<double>(domain->getNodeDistInfo());
  TriID_node = new DistVec<int>(domain->getNodeDistInfo());
  nWall_node = new DistVec<Vec3D>(domain->getNodeDistInfo());

  poly = new DistVec<bool>(domain->getNodeDistInfo());
  findPoly();

  buildSolidNormals();
  d->findNodeBoundingBoxes(X,*boxMin,*boxMax);

#pragma omp parallel for
  for(int i = 0; i < numLocSub; ++i) {
    intersector[i] = new IntersectorFRG(*(d->getSubDomain()[i]), X(i), (*status0)(i), *this);
    intersector[i]->getClosestTriangles(X(i), (*boxMin)(i), (*boxMax)(i), (*tId)(i), (*distance)(i), false);
    intersector[i]->computeFirstLayerNodeStatus((*tId)(i), (*distance)(i));
  }
  findInAndOut();
  finishStatusByPoints(iod, point_based_id);   

  if (oldStatus)
    *status = *oldStatus;
 
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) 
  {
    intersector[iSub]->findIntersections(X(iSub), false);
    for(int i = 0; i < (*is_active)(iSub).size(); ++i)
	  {
      (*is_active)(iSub)[i] = intersector[iSub]->testIsActive(0.0, i);
  }
  }

  

  ///////////////////////////////////////////////////////
  if(externalSI)
  {
	  DistVec<int>  intISactive(domain->getNodeDistInfo());
	  DistVec<bool> ISactive_bk(domain->getNodeDistInfo());

#pragma omp parallel for
	  for(int iSub = 0; iSub < numLocSub; ++iSub)
	  {
		  for(int i=0; i<X(iSub).size(); ++i) 
			  ISactive_bk(iSub)[i] = (*is_active)(iSub)[i];
		  
		  intersector[iSub]->reFlagRealNodes(X(iSub), ISactive_bk(iSub));

		  for(int i=0; i<X(iSub).size(); ++i)
		  {
			  intISactive(iSub)[i] = (*is_active)(iSub)[i] ? 1 : 0;
			  (*TriID_node)(iSub)[i] = -1;
		  }
	  }

	  operMin<int> minOp;
	  domain->assemble(domain->getLevelPat(), intISactive, minOp);

#pragma omp parallel for
	  for(int iSub = 0; iSub < numLocSub; ++iSub)
		   for(int i=0; i<X(iSub).size(); ++i)
				(*is_active)(iSub)[i] = (intISactive(iSub)[i] == 1) ? true : false;

#pragma omp parallel for
 	  for(int iSub = 0; iSub < numLocSub; ++iSub)
 		  intersector[iSub]->ComputeSIbasedIntersections(iSub, X(iSub), (*boxMin)(iSub), (*boxMax)(iSub), withViscousTerms);
	  
  }
  ///////////////////////////////////////////////////////


}

//----------------------------------------------------------------------------

void 
DistIntersectorFRG::findPoly() {
  if(!poly) {
//    com->fprintf(stderr,"ERROR: poly not initialized.\n"); 
    exit(-1);
  }

  (*poly) = false;
  DistVec<int> tester(domain->getNodeDistInfo());
  tester = 1;
  domain->assemble(domain->getLevelPat(), tester);

#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++) 
    for(int i=0; i<tester(iSub).size(); i++)
      if(tester(iSub)[i]>2)
        (*poly)(iSub)[i] = true;
}

//----------------------------------------------------------------------------

void DistIntersectorFRG::updateStructure(double *xs, double *Vs, int nNodes, int (*abc)[3]) 
{
//  com->fprintf(stderr,"DistIntersectorFRG::updateStructure called!\n");
  if(nNodes!=numStNodes) {
    com->fprintf(stderr,"Number of structure nodes has changed!\n");
    exit(-1);
  }

  for (int i=0; i<nNodes; i++) 
    for(int j=0; j<3; j++) {
      Xs_n[i][j] = Xs_np1[i][j];
      Xs_np1[i][j] = xs[3*i+j];
      Xs[i][j] = Xs_np1[i][j];
      Xsdot[i][j] = Vs[3*i+j];}

// add surface velocity contribution
  updatebc();
}

//----------------------------------------------------------------------------

void DistIntersectorFRG::updatePhysBAMInterface() 
{
  for (int i=0; i<numStNodes; i++)
    globPhysInterface->particles.X(i+1) = PhysBAM::VECTOR<double,3>(Xs[i][0], Xs[i][1], Xs[i][2]);
  globPhysInterface->Update(false); //also rebuild the topology (not really needed for now).
}

//----------------------------------------------------------------------------

/** compute the intersections, node statuses and normals for the initial geometry */
int DistIntersectorFRG::recompute(double dtf, double dtfLeft, double dts, bool findStatus, bool retry) 
{
  if(!findStatus) {
    fprintf(stderr,"ERROR: findStatus = %d. IntersectorFRG always needs to determine status.\n", findStatus);
    exit(-1);
  }

  //updateStructCoords(0.0, 1.0);
  updateStructCoords(( (dtfLeft-dtf)/dts ), ( 1.0 - (dtfLeft-dtf)/dts ));
  buildSolidNormals();
  domain->findNodeBoundingBoxes(*X,*boxMin,*boxMax);
  expandScope();

  //Debug only!
/*  com->fprintf(stderr,"count = %d.\n", ++debug_FRG_count);
  if(debug_FRG_count==50 && com->cpuNum()==0) {
    for(int i=0; i<numStNodes; i++)
      fprintf(stderr,"%d %e %e %e\n", i+1, Xs[i][0], Xs[i][1], Xs[i][2]);
    for(int i=0; i<numStElems; i++)
      fprintf(stderr,"%d 4 %d %d %d\n", i+1, stElem[i][0]+1, stElem[i][1]+1, stElem[i][2]+1);
  }*/

  int subdXing = 0;
  //bool abortOmp = false;
#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub) {
    //#pragma omp flush (abortOmp)
    //if (!abortOmp) {
      intersector[iSub]->reset(retry); 
      intersector[iSub]->rebuildPhysBAMInterface(Xs, numStNodes, stElem, numStElems);
      intersector[iSub]->getClosestTriangles((*X)(iSub), (*boxMin)(iSub), (*boxMax)(iSub), (*tId)(iSub), (*distance)(iSub), true);
      intersector[iSub]->computeFirstLayerNodeStatus((*tId)(iSub), (*distance)(iSub));
      subdXing += !(intersector[iSub]->finishStatusByHistory(*(domain->getSubDomain()[iSub])));
      if(twoPhase) {
//      intersector[iSub]->printFirstLayer(*(domain->getSubDomain()[iSub]), (*X)(iSub), 1);
        int error = intersector[iSub]->findIntersections((*X)(iSub), true);
        int nCalls = 0;
        while(error && nCalls<100) {
          nCalls++;
          com->fprintf(stderr,"Recomputing intersections (%d) ...\n", error);
          intersector[iSub]->CrossingEdgeRes.clear();
          intersector[iSub]->ReverseCrossingEdgeRes.clear();
          error = intersector[iSub]->findIntersections((*X)(iSub), true);
        }
        if(error)
          return -1;
          //abortOmp = true;
          //#pragma omp flush (abortOmp)
      }
   // }
  }

  //if (abortOmp) return -1;

  if(!twoPhase) {
    com->globalMax(1,&subdXing);
    if(subdXing) { //this happens if the structure is entering/leaving a subdomain. Rarely happens but needs to be delt with...
//      com->fprintf(stderr,"FS Interface entering/leaving a subdomain...\n");
      finalizeStatus(); //contains a local communication
    }
#pragma omp parallel for
    for(int iSub = 0; iSub < numLocSub; ++iSub) {
      int error = intersector[iSub]->findIntersections((*X)(iSub), true);
      int nCalls = 0;
      while(error && nCalls<100) {
        nCalls++;
        com->fprintf(stderr,"Recomputing intersections (%d) ...\n", error);
        intersector[iSub]->CrossingEdgeRes.clear();
        intersector[iSub]->ReverseCrossingEdgeRes.clear();
        error = intersector[iSub]->findIntersections((*X)(iSub), true);
      }
      if(error)
        return -1;
    }
  }

#pragma omp parallel for
  for(int iSub = 0; iSub < numLocSub; ++iSub)
    for(int i = 0; i < (*is_active)(iSub).size(); ++i) {
      (*is_swept)(iSub)[i] = (*status)(iSub)[i] != (*status0)(iSub)[i];
      (*is_active)(iSub)[i] = intersector[iSub]->testIsActive(0.0, i);
    }
  return 0;
}

//----------------------------------------------------------------------------

void DistIntersectorFRG::updateStructCoords(double c_n, double c_np1)
{
  for (int i=0; i<numStNodes; i++)
    Xs[i] = c_n*Xs_n[i] + c_np1*Xs_np1[i];
}

//----------------------------------------------------------------------------

void DistIntersectorFRG::findInAndOut()
{
  int nUndecided[numLocSub], total;
  DistSVec<int,2> status_and_weight(domain->getNodeDistInfo());
//  DistVec<int> one(domain->getNodeDistInfo());
//  one = 1;

#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++) 
    intersector[iSub]->floodFill(*(domain->getSubDomain()[iSub]),nUndecided[iSub]);

  while(1) { //get out only when all nodes are decided

    //1. check if all the nodes (globally) are determined
    total = 0;
    for(int iSub=0; iSub<numLocSub; iSub++)
      total += nUndecided[iSub];
    com->globalMax(1,&total);
//    com->fprintf(stderr,"total = %d\n",total);
    if(total==0) //done!
      break; 

    //2. try to get a seed from neighbor subdomains, then floodFill
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++)
      for(int i=0; i<status_and_weight(iSub).size(); i++) {
        status_and_weight(iSub)[i][0] = (*status)(iSub)[i] + 1; 
          //status(temp) = 0(UNDECIDED),1(INSIDE),or -1(OUTSIDE).
        status_and_weight(iSub)[i][1] = ((*status)(iSub)[i]==IntersectorFRG::UNDECIDED) ? 0 : 1;
      } 
         
    domain->assemble(domain->getFsPat(),status_and_weight);

#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++) {
      int nNewSeed = intersector[iSub]->findNewSeedsAfterMerging(status_and_weight(iSub), nUndecided[iSub]);
      if(nNewSeed)
        intersector[iSub]->floodFill(*(domain->getSubDomain()[iSub]),nUndecided[iSub]);      
    }
  }
}

//----------------------------------------------------------------------------
/*
void DistIntersectorFRG::findInAndOut()
{
  int nUndecided[numLocSub], total;
  DistVec<int> status_temp(domain->getNodeDistInfo());
  DistVec<int> one(domain->getNodeDistInfo());
  one = 1;

#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++)
    intersector[iSub]->floodFill(*(domain->getSubDomain()[iSub]),nUndecided[iSub]);

  while(1) { //get out only when all nodes are decided

    //1. check if all the nodes (globally) are determined
    total = 0;
    for(int iSub=0; iSub<numLocSub; iSub++)
      total += nUndecided[iSub];
    com->globalMax(1,&total);
//    com->fprintf(stderr,"total = %d\n",total);
    if(total==0) //done!
      break;

    //2. try to get a seed from neighbor subdomains, then floodFill
    status_temp = *status + one; //status_temp = 0(UNDECIDED),1(INSIDE),or -1(OUTSIDE).
    domain->assemble(domain->getLevelPat(),status_temp);
    status_temp -= one;

#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++) {
      int nNewSeed = intersector[iSub]->findNewSeedsAfterMerging(status_temp(iSub), (*poly)(iSub), nUndecided[iSub]);
      if(nNewSeed)
        intersector[iSub]->floodFill(*(domain->getSubDomain()[iSub]),nUndecided[iSub]);
    }
  }
}
*/
//----------------------------------------------------------------------------

void IntersectorFRG::reFlagRealNodes(SVec<double,3>& X, Vec<bool> &ISactive_bk)
{

	int (*ptr)[2] = edges.getPtr();

	bool i_less_j;

	for(int l=0; l<edges.size(); l++) 
	{
		if(edge_intersects[l])
		{
			int i = ptr[l][0];
			int j = ptr[l][1];

			bool isOK1 = false;
			if( (fabs(X[i][0] - 0.2127) < 1.0e-3 && fabs(X[i][1] + 0.01017) < 1.0e-4 /* && fabs(X[i][2] - 0.05) < 1.0e-4*/)  ||
				 (fabs(X[j][0] - 0.2127) < 1.0e-3 && fabs(X[j][1] + 0.01017) < 1.0e-4 /* && fabs(X[j][2] - 0.05) < 1.0e-4*/) ) isOK1 = true;
			bool isOK2 = false;
			if( (fabs(X[i][0] - 0.2127) < 1.0e-3 && fabs(X[i][1] + 0.02034) < 1.0e-4 /* && fabs(X[i][2] - 0.05) < 1.0e-4*/)  ||
				 (fabs(X[j][0] - 0.2127) < 1.0e-3 && fabs(X[j][1] + 0.02034) < 1.0e-4 /* && fabs(X[j][2] - 0.05) < 1.0e-4*/) ) isOK2 = true;

		   IntersectionResult<double>& Res_ij = CrossingEdgeRes[l];
			IntersectionResult<double>& Res_ji = ReverseCrossingEdgeRes[l];

			IntersectionResult<double> Res;
			
			if(ISactive_bk[i]) 
				i_less_j = true;
			else
				i_less_j = false;

			double alpha_1, alpha_2;
			
			if(Res_ij.triangleID > 0 && Res_ji.triangleID > 0 && Res_ij.triangleID == Res_ji.triangleID) 
			{
				Res = Res_ij;
				alpha_1 = i_less_j ? Res.alpha : 1.0 - Res.alpha;
			}  
			else if(Res_ij.triangleID > 0 && Res_ji.triangleID > 0 && Res_ij.triangleID != Res_ji.triangleID) 
			{
				Res = i_less_j ? Res_ij : Res_ji;
				alpha_1 = Res.alpha;
			}
			else if(Res_ij.triangleID > 0 && Res_ji.triangleID < 0) 
			{
				Res = Res_ij;
				alpha_1 = i_less_j ? Res.alpha : 1.0 - Res.alpha;
			}
			else if(Res_ij.triangleID<0 && Res_ji.triangleID>0) 
			{
				Res = Res_ji;
				alpha_1 = i_less_j ? 1.0 - Res.alpha : Res.alpha;
			}

			alpha_2 = 1.0 - alpha_1;
		
			if(alpha_1 > 0.5 && is_active[i]) is_active[i] = false;
			if(alpha_2 > 0.5 && is_active[j]) is_active[j] = false;
		}

	}

}

//----------------------------------------------------------------------------
void IntersectorFRG::ComputeSIbasedIntersections(int iSub, SVec<double,3>& X,  
																 SVec<double,3> &boxMin, SVec<double,3> &boxMax, 
																 bool withViscousTerms)
{	

	double l1, l2, l1_i, l1_j, l2_i, l2_j, lambda[3];

	int Tri_, Tri_i, Tri_j;

	Vec3D *structX = (distIntersector.getStructPosition()).data();

	int (*triNodes)[3] = distIntersector.stElem;

	int ntri = distIntersector.getNumStructElems();

	MyTriangle *myTris = new MyTriangle[ntri];

	for(int i = 0; i < ntri; ++i)
      myTris[i] = MyTriangle(i, distIntersector.getStructPosition(), triNodes[i]);

	KDTree<MyTriangle> structureTree(ntri, myTris);

	ClosestTriangle closestTriangle(triNodes, structX, distIntersector.triNorms, distIntersector.node2node, distIntersector.node2elem);

	int nMaxCand_i = 500, nMaxCand_j = 500;

	MyTriangle *candidates_i = new MyTriangle[nMaxCand_i];
	MyTriangle *candidates_j = new MyTriangle[nMaxCand_j];

	int (*ptr)[2] = edges.getPtr();

	for(int l=0; l<edges.size(); l++)
	{
		int i = ptr[l][0];
		int j = ptr[l][1];

		edge_SI[l] = false;
	  
		if(!is_active[i] || !is_active[j])
		{
			if(is_active[i] == is_active[j]) continue;

			Vec3D X_ij;
			for(int k=0; k<3; ++k) X_ij[k] = 0.5*(X[j][k] + X[i][k]);

			edge_SI[l] = true;

			int nFound_i = structureTree.findCandidatesInBox(boxMin[i], boxMax[i], candidates_i, nMaxCand_i);

			if(nFound_i > nMaxCand_i) 
			{
				nMaxCand_i = nFound_i;

				delete [] candidates_i;

				candidates_i = new MyTriangle[nMaxCand_i];

				structureTree.findCandidatesInBox(boxMin[i], boxMax[i], candidates_i, nMaxCand_i);
			}

			int nFound_j = structureTree.findCandidatesInBox(boxMin[j], boxMax[j], candidates_j, nMaxCand_j);

			if(nFound_j > nMaxCand_j) 
			{
				nMaxCand_j = nFound_j;

				delete [] candidates_j;

				candidates_j = new MyTriangle[nMaxCand_j];

				structureTree.findCandidatesInBox(boxMin[j], boxMax[j], candidates_j, nMaxCand_j);
			}

			int trId;
			double dist_si, dist_i, dist_j, xi[3];
			
			double min_dist_si = 1.0e16;
			double min_dist_i  = 1.0e16;
			double min_dist_j  = 1.0e16;

			Tri_ =-1;
			l1 = l2 = -1.0;

			Tri_i = Tri_j =-1;
			l1_i = l1_j = l2_i = l2_j = -1.0;

			if(nFound_i > 0)
			{
				for(int l = 0; l < nFound_i; ++l) 
				{
					int trId = candidates_i[l].trId();

					dist_si = closestTriangle.piercing(X_ij, trId, lambda);

					if(dist_si < min_dist_si) 
					{
						min_dist_si = dist_si;

						Tri_ = trId;
						l1   = lambda[0];
						l2   = lambda[1];
					}

				}
			}

			if(nFound_j > 0)
			{
				for(int l = 0; l < nFound_j; ++l) 
				{
					int trId = candidates_j[l].trId();

					dist_si = closestTriangle.piercing(X_ij, trId, lambda);

					if(dist_si < min_dist_si)
					{
						min_dist_si = dist_si;

						Tri_ = trId;
						l1   = lambda[0];
						l2   = lambda[1];
					}

				}
			}

			Vec3D d2, d3, nn;

			for(int id=0; id<3; ++id)
			{
				d2[id] = structX[triNodes[Tri_][1]][id] - structX[triNodes[Tri_][0]][id];
				d3[id] = structX[triNodes[Tri_][2]][id] - structX[triNodes[Tri_][0]][id];
			}

			nn = d2^d3;

			double norm = sqrt(nn * nn);
			if(norm != 0.0) nn *= (1.0 / norm);
			
			   xi_SI[l] = l1;
		     eta_SI[l] = l2;
			nWall_SI[l] = nn;
			TriID_SI[l] = Tri_;

			/*
			PhysBAMInterface<double>& physbam_interface=*distIntersector.physInterface;

			ARRAY<int> candidates_i;
			VECTOR<double,3> min_corner_i(boxMin[i][0],boxMin[i][1],boxMin[i][2]);
			VECTOR<double,3> max_corner_i(boxMax[i][0],boxMax[i][1],boxMax[i][2]);		
			bool t_i = physbam_interface.HasCloseTriangle(locIndex+1, VECTOR<double,3>(X[i][0],X[i][1],X[i][2]), 
																		 min_corner_i, max_corner_i, &shrunk_index, &occluded, &candidates_i);
																		 


			ARRAY<int> candidates_j;
			VECTOR<double,3>min_corner_j(boxMin[j][0],boxMin[j][1],boxMin[j][2]);
			VECTOR<double,3>max_corner_j(boxMax[j][0],boxMax[j][1],boxMax[j][2]);
			bool t_j= physbam_interface.HasCloseTriangle(locIndex+1, VECTOR<double,3>(X[i][0],X[i][1],X[i][2]), 
																		min_corner_j, max_corner_j, &shrunk_index, &occluded, &candidates_j);

			int trId;		
			double dist_si, dist_i, dist_j, xi[3];
			
			double min_dist_si = 1.0e16;
			double min_dist_i  = 1.0e16;
			double min_dist_j  = 1.0e16;

			Tri_ =-1;
			l1 = l2 = -1.0;

			Tri_i = Tri_j =-1;
			l1_i = l1_j = l2_i = l2_j = -1.0;
			
			if(t_i) 
			{
				for(int iArray=1; iArray <= candidates_i.Size(); iArray++) 
				{					
					trId = candidates_i(iArray) - 1;
			
					dist_si = piercing(X_ij, trId, lambda);

					if(dist_si < min_dist_si) 
					{
						min_dist_si = dist_si;

						Tri_ = trId;
						l1   = lambda[0];
						l2   = lambda[1];
					}

					if(withViscousTerms && !is_active[i])
					{
						dist_i = piercing(X[i], trId, lambda);

						if(dist_i < min_dist_i)
						{
							min_dist_i = dist_i;

							Tri_i = trId;
							l1_i  = lambda[0];
							l2_i  = lambda[1];
						}
					}
					
				}
			}

			if(t_j) 
			{		  
				for(int iArray=1; iArray <= candidates_j.Size(); iArray++) 
				{
					trId = candidates_j(iArray) - 1;

					dist_si = piercing(X_ij, trId, lambda);

					if(dist_si < min_dist_si) 
					{
						min_dist_si = dist_si;

						Tri_ = trId;
						l1   = lambda[0];
						l2   = lambda[1];
					}

					if(withViscousTerms && !is_active[j])
					{
						dist_j = piercing(X[j], trId, lambda);

						if(dist_j < min_dist_j)
						{
							min_dist_j = dist_j;

							Tri_j = trId;
							l1_j  = lambda[0];
							l2_j  = lambda[1];
						}
					}

				}
			}
			
			Vec3D d2, d3, nn;

			for(int id=0; id<3; ++id)
			{
				d2[id] = structX[triNodes[Tri_][1]][id] - structX[triNodes[Tri_][0]][id];
				d3[id] = structX[triNodes[Tri_][2]][id] - structX[triNodes[Tri_][0]][id];
			}

			nn = d2^d3;

			double norm = sqrt(nn * nn);
			if(norm != 0.0) nn *= (1.0 / norm);
			
			   xi_SI[l] = l1;
		     eta_SI[l] = l2;
			nWall_SI[l] = nn;
			TriID_SI[l] = Tri_;

			if(withViscousTerms && !is_active[i])
			{
				for(int id=0; id<3; ++id)
				{
					d2[id] = structX[triNodes[Tri_i][1]][id] - structX[triNodes[Tri_i][0]][id];
					d3[id] = structX[triNodes[Tri_i][2]][id] - structX[triNodes[Tri_i][0]][id];
				}

				nn = d2^d3;
				
				double norm = sqrt(nn * nn);
				if(norm != 0.0) nn *= 1.0 / norm;

				   xi_node[i] = l1_i;
 				  eta_node[i] = l2_i;
 				nWall_node[i] = nn;
				TriID_node[i] = Tri_i;
			}

			if(withViscousTerms && !is_active[j])
			{
				for(int id=0; id<3; ++id)
				{
					d2[id] = structX[triNodes[Tri_j][1]][id] - structX[triNodes[Tri_j][0]][id];
					d3[id] = structX[triNodes[Tri_j][2]][id] - structX[triNodes[Tri_j][0]][id];
				}

				nn = d2^d3;
				
				double norm = sqrt(nn * nn);
				if(norm != 0.0) nn *= 1.0 / norm;

				   xi_node[j] = l1_j;
 				  eta_node[j] = l2_j;
 				nWall_node[j] = nn;
				TriID_node[j] = Tri_j;
			}
		*/		
		}

	}

/* --------- debug --------- */
#if 0
	for (int l=0; l<edges.size(); l++)
	{
		int i = ptr[l][0]; int j = ptr[l][1];

		if(edge_SI[l])
		{			
			if( (!is_active[i] && !is_active[j]) ||  
				 ( is_active[i] &&  is_active[j]) )
			{		
				std::cout << "error SI edge and " << std::boolalpha << is_active[i] << " " << is_active[j] << " " << edge_intersects[l] << endl;
				std::cout << i << " Xi = " << X[i][0]<<" "<<X[i][1]<<" "<<X[i][2] 
							 << " " << j << " Xj = " << X[j][0]<<" "<<X[j][1]<<" "<<X[j][2] <<endl;
				exit(-1);
			}
		}
		else
		{
			if( is_active[i] != is_active[j] )
			{	
				std::cout << "error  edge without SI and " << std::boolalpha << is_active[i] << " " << is_active[j] 
							 << " " << edge_intersects[l] << endl;
				std::cout << i << " Xi = " << X[i][0]<<" "<<X[i][1]<<" "<<X[i][2] << " " 
							 << j << " Xj = " << X[j][0]<<" "<<X[j][1]<<" "<<X[j][2] <<endl;
				exit(-1);
			}
		}

	}
#endif

}

//----------------------------------------------------------------------------

void IntersectorFRG::printFirstLayer(SubDomain& sub, SVec<double,3>&X, int TYPE)
{
  int mySub = sub.getGlobSubNum();
  int myLocSub = sub.getLocSubNum();
  int (*ptr)[2] = edges.getPtr();
  Connectivity &nToN = *(sub.getNodeToNode());
  char fileName[50] = "firstLayera.top";
  char nodesName[2] = "a";

  fileName[10] += mySub;
  nodesName[0] += mySub;

  FILE* firstLayer = fopen(fileName,"w");
  fprintf(firstLayer, "Nodes InsideNodes%s\n", nodesName);
  for (int i=0; i<sub.numNodes(); i++)
    if (status[i]==TYPE) fprintf(firstLayer,"%d %e %e %e\n", i+1, X[i][0], X[i][1], X[i][2]);
  fprintf(firstLayer, "Elements FirstLayer%s using InsideNodes%s\n", nodesName, nodesName);
  for (int l=0; l<edges.size(); l++){
    int x1 = ptr[l][0], x2 = ptr[l][1];
    if (status[x1]!=TYPE || status[x2]!=TYPE) continue;
    int crit = 0;
    for (int i=0; i<nToN.num(x1); i++)
      if (status[nToN[x1][i]]!=TYPE) {crit++; break;}
    for (int i=0; i<nToN.num(x2); i++)
      if (status[nToN[x2][i]]!=TYPE) {crit++; break;}
    if (crit==2)
      fprintf(firstLayer,"%d %d %d %d\n", l+1, (int)1, x1+1, x2+1);
  }
  fclose(firstLayer);

}

//----------------------------------------------------------------------------

void DistIntersectorFRG::finishStatusByPoints(IoData &iod, DistVec<int> *point_based_id)
{
  if((numFluid==1 || numFluid==2) && iod.embed.embedIC.pointMap.dataMap.empty()) {
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++)
      for(int i=0; i<(*status)(iSub).size(); i++)
        if((*status)(iSub)[i]==IntersectorFRG::OUTSIDE)
          (*status)(iSub)[i]=1; 

    twoPhase = true; 
    return;
  }

  list< pair<Vec3D,int> > Points; //pair points with fluid model ID.
  map<int,int> pid2id;

  if(!iod.embed.embedIC.pointMap.dataMap.empty()){
    map<int, PointData *>::iterator pointIt;
    int count = 0;
    for(pointIt  = iod.embed.embedIC.pointMap.dataMap.begin();
        pointIt != iod.embed.embedIC.pointMap.dataMap.end();
        pointIt ++){
      int myID = pointIt->second->fluidModelID;
      int myPID = ++count;
      Vec3D xyz(pointIt->second->x, pointIt->second->y,pointIt->second->z);

      if(point_based_id) {
        Points.push_back(pair<Vec3D,int>(xyz, myPID));
        pid2id[myPID] = myID;
      } else
        Points.push_back(pair<Vec3D,int>(xyz, myID));

      if(myID>=numFluid) { //myID should start from 0
        com->fprintf(stderr,"ERROR:FluidModel %d doesn't exist! NumPhase = %d\n", myID, numFluid);
        exit(-1);
      } 
    }
  } else {
    com->fprintf(stderr, "ERROR: (INTERSECTOR) Point-based initial conditions are required for multi-phase FSI.\n");
    exit(-1);
  }

  list< pair<Vec3D,int> >::iterator iter;

  int nUndecided[numLocSub], total;
  DistSVec<int,2> status_and_weight(domain->getNodeDistInfo());

  // first round
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++) {
    nUndecided[iSub] = 0;

    // 1. move "OUTSIDE" nodes to "UNDECIDED".
    for(int i=0; i<(*status)(iSub).size(); i++)
      if((*status)(iSub)[i]==IntersectorFRG::OUTSIDE) {
        (*status)(iSub)[i] = IntersectorFRG::UNDECIDED;
        nUndecided[iSub]++;
      }
    // 2. find seeds by points
    int nSeeds = intersector[iSub]->findSeedsByPoints(*(domain->getSubDomain()[iSub]), (*X)(iSub), Points, nUndecided[iSub]);
    // 3. flood fill if seeds are found
    if(nSeeds>0)
      intersector[iSub]->noCheckFloodFill(*(domain->getSubDomain()[iSub]),nUndecided[iSub]);
  }   

  int total0 = 0;
  while(1) { //get out only when all nodes are decided
    //1. check if all the nodes (globally) are determined
    total = 0;
#pragma omp parallel for reduction(+: total)
    for(int iSub=0; iSub<numLocSub; iSub++)
      total += nUndecided[iSub];
    com->globalSum(1,&total);
    if(total==0 || total==total0)
      break; //done
    total0 = total;

    //2. try to get a seed from neighbor subdomains
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++)
      for(int i=0; i<status_and_weight(iSub).size(); i++) {
        status_and_weight(iSub)[i][0] = (*status)(iSub)[i] + 1;
          //status(temp) = 0(UNDECIDED),1(INSIDE),or 2,3, ...
        status_and_weight(iSub)[i][1] = ((*status)(iSub)[i]==IntersectorFRG::UNDECIDED) ? 0 : 1;
      }

    domain->assemble(domain->getFsPat(),status_and_weight);

#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++) {
      int nNewSeed = intersector[iSub]->findNewSeedsAfterMerging(status_and_weight(iSub), nUndecided[iSub]);
      if(nNewSeed)
        intersector[iSub]->noCheckFloodFill(*(domain->getSubDomain()[iSub]),nUndecided[iSub]);
    }
  }

  if(point_based_id) {
    *point_based_id = -999;
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++)
      for(int i=0; i<(*status)(iSub).size(); i++) {
        int myPID = (*status)(iSub)[i];
        map<int,int>::iterator pit = pid2id.find(myPID);
        if(pit!=pid2id.end()) {
          (*point_based_id)(iSub)[i] = myPID;
          (*status)(iSub)[i] = pit->second;
        }
      }
  }

  if(total) {//still have undecided nodes. They must be ghost nodes (i.e. covered by solid).
    com->fprintf(stderr,"- IntersectorFRG: Ghost node(s) detected...\n");
    twoPhase = false;
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++)
      for(int i=0; i<(*status)(iSub).size(); i++) {
        if((*status)(iSub)[i]==IntersectorFRG::UNDECIDED) {
          (*status)(iSub)[i] = IntersectorFRG::OUTSIDECOLOR;
        }
      }
  } else 
    twoPhase = (numFluid<3) ? true : false;
}

//----------------------------------------------------------------------------
/*
void DistIntersectorFRG::finishStatusByPoints(IoData &iod)
{
  if((numFluid==1 || numFluid==2) && iod.embed.embedIC.pointMap.dataMap.empty()) {
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++)
      for(int i=0; i<(*status)(iSub).size(); i++)
        if((*status)(iSub)[i]==IntersectorFRG::OUTSIDE)
          (*status)(iSub)[i]=1; 

    twoPhase = true; 
    return;
  }

  list< pair<Vec3D,int> > Points; //pair points with fluid model ID.
  if(!iod.embed.embedIC.pointMap.dataMap.empty()){
    map<int, PointData *>::iterator pointIt;
    for(pointIt  = iod.embed.embedIC.pointMap.dataMap.begin();
        pointIt != iod.embed.embedIC.pointMap.dataMap.end();
        pointIt ++){
      int myID = pointIt->second->fluidModelID;
      Vec3D xyz(pointIt->second->x, pointIt->second->y,pointIt->second->z);
      Points.push_back(pair<Vec3D,int>(xyz, myID));

      if(myID>=numFluid) { //myID should start from 0
        com->fprintf(stderr,"ERROR:FluidModel %d doesn't exist! NumPhase = %d\n", myID, numFluid);
        exit(-1);
      } 
    }
  } else {
    com->fprintf(stderr, "ERROR: (INTERSECTOR) Point-based initial conditions are required for multi-phase FSI.\n");
    exit(-1);
  }

  list< pair<Vec3D,int> >::iterator iter;
  for(iter = Points.begin(); iter!=Points.end(); iter++)
    com->fprintf(stderr,"  - Detected point (%e %e %e) with FluidModel %d\n", (iter->first)[0], (iter->first)[1], (iter->first)[2], iter->second);
  

  int nUndecided[numLocSub], total;
  DistVec<int> status_temp(domain->getNodeDistInfo());
  DistVec<int> one(domain->getNodeDistInfo());
  one = 1;

  // first round
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++) {
    nUndecided[iSub] = 0;

    // 1. move "OUTSIDE" nodes to "UNDECIDED".
    for(int i=0; i<(*status)(iSub).size(); i++)
      if((*status)(iSub)[i]==IntersectorFRG::OUTSIDE) {
        (*status)(iSub)[i] = IntersectorFRG::UNDECIDED;
        nUndecided[iSub]++;
      }
    // 2. find seeds by points
    int nSeeds = intersector[iSub]->findSeedsByPoints(*(domain->getSubDomain()[iSub]), (*X)(iSub), Points, nUndecided[iSub]);
    // 3. flood fill if seeds are found
    if(nSeeds>0)
      intersector[iSub]->noCheckFloodFill(*(domain->getSubDomain()[iSub]),nUndecided[iSub]);
  }   

  int total0 = 0;
  while(1) { //get out only when all nodes are decided
    //1. check if all the nodes (globally) are determined
    total = 0;
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++)
      total += nUndecided[iSub];
    com->globalSum(1,&total);
    com->fprintf(stderr,"Total number of undecided nodes = %d\n", total);    
    if(total==0 || total==total0)
      break; //done

    total0 = total;
    //2. try to get a seed from neighbor subdomains
    status_temp = *status + one; //status_temp = 0,1,2,3,....
    domain->assemble(domain->getLevelPat(),status_temp);
    status_temp -= one;
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++) {
      int nNewSeed = intersector[iSub]->findNewSeedsAfterMerging(status_temp(iSub), (*poly)(iSub), nUndecided[iSub]);
      if(nNewSeed)
        intersector[iSub]->noCheckFloodFill(*(domain->getSubDomain()[iSub]),nUndecided[iSub]);
    }
  }

  if(total) {//still have undecided nodes. They must be ghost nodes (i.e. covered by solid).
    twoPhase = false;
#pragma omp parallel for
    for(int iSub=0; iSub<numLocSub; iSub++)
      for(int i=0; i<(*status)(iSub).size(); i++) {

        if(com->cpuNum()==154&&intersector[iSub]->locToGlobNodeMap[i]==700111-1) {
          fprintf(stderr,"Node 700111 belongs to %d subdomains.\n", intersector[iSub]->nodeToSubD.num(i));
          Connectivity &nToN = *(intersector[iSub]->subD->getNodeToNode());
          for(int j=0; j<nToN.num(i); j++)
            fprintf(stderr,"  -- Neighbour: %d\n", intersector[iSub]->locToGlobNodeMap[nToN[i][j]]+1);
        }

        if(intersector[iSub]->locToGlobNodeMap[i]==700111-1)
          fprintf(stderr,"CPU %d has ndoe 700111 with status %d.\n", com->cpuNum(), (*status)(iSub)[i]);


        if((*status)(iSub)[i]==IntersectorFRG::UNDECIDED) {
          fprintf(stderr,"CPU %d: Node %d is undecided!!!\n", com->cpuNum(), intersector[iSub]->locToGlobNodeMap[i]+1);
          (*status)(iSub)[i] = IntersectorFRG::OUTSIDECOLOR;
        }
      }
  } else 
    twoPhase = (numFluid<3) ? true : false;

  exit(-1);
}
*/
//----------------------------------------------------------------------------

void DistIntersectorFRG::finalizeStatus()
{
  DistSVec<int,2> status_and_weight(domain->getNodeDistInfo());
#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++)
    for(int i=0; i<status_and_weight(iSub).size(); i++) {
      if((*status)(iSub)[i]==IntersectorFRG::OUTSIDE) {//this means its real status hasn't been decided.
        status_and_weight(iSub)[i][0] = 0;
        status_and_weight(iSub)[i][1] = 0;
      } else {
        status_and_weight(iSub)[i][0] = (*status)(iSub)[i];
        status_and_weight(iSub)[i][1] = 1;
      }
    
//      if(intersector[iSub]->locToGlobNodeMap[i]+1==151870)
//        fprintf(stderr,"Before Assemble: CPU%d: status/weight of 151870 is %d/%d...\n", com->cpuNum(), status_and_weight(iSub)[i][0],status_and_weight(iSub)[i][1]);
    }  


  domain->assemble(domain->getFsPat(),status_and_weight);

#pragma omp parallel for
  for(int iSub=0; iSub<numLocSub; iSub++)
    for(int i=0; i<(*status)(iSub).size(); i++) {
//      if(intersector[iSub]->locToGlobNodeMap[i]+1==151870)
//        fprintf(stderr,"After Assemble: CPU%d: status/weight of 151870 is %d/%d...\n", com->cpuNum(), status_and_weight(iSub)[i][0],status_and_weight(iSub)[i][1]);
      if((*status)(iSub)[i]==IntersectorFRG::OUTSIDE) {
        if(status_and_weight(iSub)[i][1]!=0)
          (*status)(iSub)[i] = status_and_weight(iSub)[i][0] / status_and_weight(iSub)[i][1];

        if((*status)(iSub)[i]<=0 || (*status)(iSub)[i]>numFluid)
          fprintf(stderr,"ERROR: failed at determining status for node %d. Structure may have crossed two or more layers of nodes in one time-step!\n", 
                  (intersector[iSub]->locToGlobNodeMap)[i]+1);
      }
//      if((intersector[iSub]->locToGlobNodeMap)[i]+1==81680)
//        fprintf(stderr,"Node 81680 has status %d.....\n", (*status)(iSub)[i]);
    }
}

//----------------------------------------------------------------------------

bool IntersectorFRG::finishStatusByHistory(SubDomain& sub)
{
  bool good = true;
  bool two_phase = distIntersector.twoPhase;

  int numNodes = status.size();
  Connectivity &nToN = *(sub.getNodeToNode());
  for(int i=0; i<numNodes; i++) {
    if(status[i]==OUTSIDE){
      if(two_phase) {
        status[i] = 1;
        continue;
      }

      if(status0[i]!=INSIDE)
        status[i] = status0[i];
      else {//status0[i]==INSIDE
        bool done = false;    
        for(int iNei=0; iNei<nToN.num(i); iNei++)
          if(status0[nToN[i][iNei]]!=INSIDE) {
            status[i] = status0[nToN[i][iNei]];
            done = true;
            break;
          }
        if(!done) {
          good = false; 
//          fprintf(stderr,"ERROR(CPU %d): F-S Interface moved over 2 (or more) layers of nodes in 1 time-step!\n", distIntersector.com->cpuNum()); 
//          fprintf(stderr,"Node %d: status = %d, status0 = %d\n",locToGlobNodeMap[i]+1,status[i],status0[i]);
//          for(int iNei=0; iNei<nToN.num(i); iNei++) 
//            fprintf(stderr,"Neighbor %d: status = %d, status0 = %d\n", locToGlobNodeMap[nToN[i][iNei]]+1, status[nToN[i][iNei]], status0[nToN[i][iNei]]);
        }
      }
    } 
    else if(status[i]==UNDECIDED) {
      if(status0[i]!=UNDECIDED)
        status[i] = status0[i];
      else  
        fprintf(stderr,"ERROR: Unable to determine node status for Node %d\n", locToGlobNodeMap[i]+1);
    }

//    if(locToGlobNodeMap[i]+1==151870)
//      fprintf(stderr,"CPU%d: status of 151870 is %d...\n", distIntersector.com->cpuNum(), status[i]);

  }

  return good;
}

//----------------------------------------------------------------------------

void IntersectorFRG::floodFill(SubDomain& sub, int& nUndecided)
{
  nUndecided = 0;
  int numNodes = status.size();
  Connectivity &nToN = *(sub.getNodeToNode());

  // Propogate status to the entire subdomain.
  // List is used as a FIFO queue
  Vec<int> seed(numNodes); //list of decided nodes  
  Vec<int> level(numNodes);
  // Look for a start point
  // lead: In pointer
  // next: out pointer (of FIFO)
  int next = 0, lead = 0;
  for(int i = 0; i < numNodes; ++i) {
    if(status[i] != UNDECIDED) {
      seed[lead++] = i;
      level[i] = 0;
    } else
      nUndecided++;
  }

  while(next < lead) { //still have seeds not used
    int cur = seed[next++];
    int curStatus = status[cur];
    int curLevel = level[cur];
    for(int i = 0; i < nToN.num(cur); ++i) {
      if(status[nToN[cur][i]] == UNDECIDED) {
        status[nToN[cur][i]] = curStatus;
        level[nToN[cur][i]] = curLevel+1;
        seed[lead++] = nToN[cur][i]; 
        nUndecided--;
      } else 
        if(status[nToN[cur][i]] != curStatus && ( curLevel != 0 || level[nToN[cur][i]] != 0))
          std::cerr << "Incompatible nodes have met: " << locToGlobNodeMap[cur]+1 << "("<< status[cur]
                    << ") and " << locToGlobNodeMap[nToN[cur][i]]+1 << "(" << status[nToN[cur][i]] << ") "
                    << " " << curLevel << " " << level[nToN[cur][i]] << std::endl;
    }
  }
}

//----------------------------------------------------------------------------

IntersectorFRG::IntersectorFRG(SubDomain &sub, SVec<double,3> &X,
                               Vec<int> &stat0, DistIntersectorFRG &distInt) :
    LevelSetStructure((*distInt.status)(sub.getLocSubNum()),(*distInt.distance)(sub.getLocSubNum()),
                      (*distInt.is_swept)(sub.getLocSubNum()),(*distInt.is_active)(sub.getLocSubNum()),
                      (*distInt.is_occluded)(sub.getLocSubNum()),(*distInt.edge_intersects)(sub.getLocSubNum()),
							 (*distInt.edge_SI)(sub.getLocSubNum()),
							 (*distInt.xi_SI)(sub.getLocSubNum()),(*distInt.eta_SI)(sub.getLocSubNum()),
		                (*distInt.nWall_SI)(sub.getLocSubNum()),
							 (*distInt.TriID_SI)(sub.getLocSubNum()),
							 (*distInt.xi_node)(sub.getLocSubNum()),(*distInt.eta_node)(sub.getLocSubNum()),
		                (*distInt.nWall_node)(sub.getLocSubNum()),
							 (*distInt.TriID_node)(sub.getLocSubNum()) ),
    subD(&sub), distIntersector(distInt), status0(stat0),
    edges(sub.getEdges()), globIndex(sub.getGlobSubNum()), nodeToSubD(*sub.getNodeToSubD())
{
  physInterface = 0;

  status0 = UNDECIDED;

  OUTSIDECOLOR = distInt.numOfFluids();

  locToGlobNodeMap = sub.getNodeMap();

  iscope = 0;
  package = new set<int>[sub.getNumNeighb()];
  int *neighb = sub.getNeighb();
  for(int i=0; i<sub.getNumNeighb(); i++)
    sub2pack[neighb[i]] = i;

  nFirstLayer = 0;
  status = UNDECIDED;
  distance = 0;
  is_swept = false;
  is_active = false;
  is_occluded = false;
  edge_intersects = false;
  edge_SI = false;

}

//----------------------------------------------------------------------------

IntersectorFRG::~IntersectorFRG()
{
  delete[] package;
  if(iscope) delete[] iscope;

  if(physInterface) delete physInterface;
}

//----------------------------------------------------------------------------

void IntersectorFRG::reset(const bool retry)
{
  for(int i=0; i<subD->getNumNeighb(); i++)
    package[i].clear();

  if(iscope) {delete[] iscope; iscope = 0;}
  n2p.clear();
  particle.clear();

  CrossingEdgeRes.clear();
  ReverseCrossingEdgeRes.clear();

  if(!retry) status0 = status;
  status = UNDECIDED;
  if(physInterface) {delete physInterface; physInterface = 0;}

  nFirstLayer = 0;
  distance = 0.0;
  is_swept = false;
  is_active = false;
  edge_intersects = false;
  edge_SI = false;
}

//----------------------------------------------------------------------------
//TODO: discuss with Jon
void IntersectorFRG::rebuildPhysBAMInterface(Vec3D *Xs, int nsNodes, int (*sElem)[3], int nsElems)
{ //IMPORTANT: In PhysBAM array index starts fr*om 1 instead of 0
  int nPar, nTri = scope.size();
  if(nTri==0)
    return;
  nPar = buildScopeTopology(sElem, nsElems);

  int count = 0;
  // Initialize the Particles list
  PhysBAM::GEOMETRY_PARTICLES<PhysBAM::VECTOR<double,3> > *physbam_solids_particle=new PhysBAM::GEOMETRY_PARTICLES<PhysBAM::VECTOR<double,3> >();
  physbam_solids_particle->array_collection.Resize(nPar);
  for (list<int>::iterator lit=particle.begin(); lit!=particle.end(); lit++)
    physbam_solids_particle->X(++count) = PhysBAM::VECTOR<double,3>(Xs[*lit][0],Xs[*lit][1], Xs[*lit][2]);

  // Initialize the Triangle list
  PhysBAM::ARRAY<PhysBAM::VECTOR<int,3> > physbam_stElem;
  physbam_stElem.Preallocate(nTri);
  for (set<int>::iterator it=scope.begin(); it!=scope.end(); it++)
    physbam_stElem.Append(PhysBAM::VECTOR<int,3>(n2p[sElem[*it][0]]+1, n2p[sElem[*it][1]]+1, n2p[sElem[*it][2]]+1));

  // Initialize the mesh.
  // fprintf(stderr,"Initializing the Mesh with %d particles and %d triangles\n",physbam_solids_particle->array_collection.Size(),physbam_stElem.Size());
  PhysBAM::TRIANGLE_MESH *mesh = new PhysBAM::TRIANGLE_MESH(count,physbam_stElem);
  mesh->Initialize_Adjacent_Elements();mesh->Set_Number_Nodes(count);

  // Construct TRIANGULATED_SURFACE.
  if(physInterface) delete physInterface;
  physInterface = new PhysBAMInterface<double>(*mesh,*physbam_solids_particle);





/*







  // Initialize the Particles list 
  int count = 0;
  physbam_solids_particle = new PhysBAM::GEOMETRY_PARTICLES<PhysBAM::VECTOR<double,3> >();
  physbam_solids_particle->array_collection.Resize(nPar);
  for (list<int>::iterator lit=particle.begin(); lit!=particle.end(); lit++) 
    physbam_solids_particle->X(++count) = PhysBAM::VECTOR<double,3>(Xs[*lit][0],Xs[*lit][1], Xs[*lit][2]);

  // Initialize the Triangle list
  physbam_stElem = new PhysBAM::ARRAY<PhysBAM::VECTOR<int,3> >();
  for (set<int>::iterator it=scope.begin(); it!=scope.end(); it++)
    physbam_stElem->Append(PhysBAM::VECTOR<int,3>(n2p[sElem[*it][0]]+1, n2p[sElem[*it][1]]+1, n2p[sElem[*it][2]]+1));

  // Construct TRIANGLE_MESH triangle_mesh (it stores its own copy of physbam_stElem).
  physbam_triangle_mesh = new PhysBAM::TRIANGLE_MESH(physbam_solids_particle->array_collection.Size(), *physbam_stElem);
  physbam_triangle_mesh->Initialize_Adjacent_Elements(); //need this?

  // Construct TRIANGULATED_SURFACE. (only stores the reference to physbam_triangle_mesh and physbam_solids_particle)
  physbam_triangulated_surface = new PhysBAM::TRIANGULATED_SURFACE<double>(*physbam_triangle_mesh, *physbam_solids_particle);
  physbam_triangulated_surface->Update_Triangle_List();

  // Construct PhysBAMInterface. (only stores a reference to physbam_triangulated_surface)
  physInterface = new PhysBAMInterface<double>(*physbam_triangulated_surface);
*/
/*  char ch[20] = "physbam_cpuX";
  ch[11] = '0' + subD->getGlobSubNum();
  FILE* myDebug = fopen(ch,"w");
  fprintf(myDebug,"Particles: \n");
  for(int i=0; i<nPar; i++)
    fprintf(myDebug, "%d %e %e %e\n", i+1, physbam_solids_particle.X(i+1)[1], physbam_solids_particle.X(i+1)[2], physbam_solids_particle.X(i+1)[3]);
  fprintf(myDebug,"Node to Particle map:\n");
  for(map<int,int>::iterator it=n2p.begin(); it!=n2p.end(); it++)
    fprintf(myDebug,"%d -> %d\n", it->first, it->second);
  fprintf(myDebug,"Elements:\n");
  for(int i=0; i<nTri; i++)
    fprintf(myDebug,"%d: %d %d %d\n", i+1, physbam_stElem(i+1)[1], physbam_stElem(i+1)[2], physbam_stElem(i+1)[3]);
  fclose(myDebug);
*/
}

//----------------------------------------------------------------------------
/*
void IntersectorFRG::buildSolidNormals(Vec3D *Xs, int nsNodes, int (*sElem)[3], int nsElem)
{
  int nTr = scope.size();
  if(nTr==0)
    return;
  int nPt = particle.size();
  bool interp = distIntersector.interpolatedNormal;

  // allocate memory.
  if(trNormal) delete[] trNormal;
  trNormal = new Vec3D[nTr];
  if(interp) {
    if(ndNormal) delete[] ndNormal;
    ndNormal = new Vec3D[nPt];
    for(int i=0; i<nPt; i++)
      ndNormal[i] = 0.0;
  }

  for(int iTr=0; iTr<nTr; iTr++) {
    int n1 = sElem[iscope[iTr]][0];
    int n2 = sElem[iscope[iTr]][1];
    int n3 = sElem[iscope[iTr]][2];
    double x1 = Xs[n1][0];
    double y1 = Xs[n1][1];
    double z1 = Xs[n1][2];
    double dx2 = Xs[n2][0]-x1;
    double dy2 = Xs[n2][1]-y1;
    double dz2 = Xs[n2][2]-z1;
    double dx3 = Xs[n3][0]-x1;
    double dy3 = Xs[n3][1]-y1;
    double dz3 = Xs[n3][2]-z1;

    // calculate the normal.
    trNormal[iTr] = Vec3D(dx2,dy2,dz2)^Vec3D(dx3,dy3,dz3);

    if(interp){ // compute nodal normal (weighted by 2*area)
      ndNormal[n2p[n1]] += trNormal[iTr];
      ndNormal[n2p[n2]] += trNormal[iTr];
      ndNormal[n2p[n3]] += trNormal[iTr];
    } 

    // normalize the normal.
    double nrm = trNormal[iTr].norm();
    if(nrm > 0.0)
       trNormal[iTr] /= nrm;
    else
      fprintf(stderr,"ERROR: Area of Triangle %d is %e.\n", iscope[iTr]+1, 0.5*nrm); 
  }

  if(interp) //normalize nodal normals.
    for(int i=0; i<nPt; i++) 
      ndNormal[i] /= ndNormal[i].norm(); //TODO: need a local communication
}
*/
//----------------------------------------------------------------------------

int IntersectorFRG::buildScopeTopology(int (*sElem)[3], int nsElem)
{ // construct iscope, n2p, and particles.
  iscope = new int[scope.size()]; //it's memory has been released in "reset".
  int nd, newID = 0, newTrID = 0;
  for(set<int>::iterator it=scope.begin(); it!=scope.end(); it++){
    iscope[newTrID++] = *it;
    for(int k=0; k<3; k++) {
      nd = sElem[*it][k];
      if(n2p.find(nd)==n2p.end()) {
        n2p[nd] = newID++;
        particle.push_back(nd);
      }
    }
  }
  return newID;
}

//----------------------------------------------------------------------------

/** Find the closest structural triangle for each node. If no triangle intersect the bounding box of the node,
* no closest triangle exists
*/
void IntersectorFRG::getClosestTriangles(SVec<double,3> &X, SVec<double,3> &boxMin, SVec<double,3> &boxMax, Vec<int> &tId, Vec<double> &dist, bool useScope) 
{
  int (*triNodes)[3] = distIntersector.stElem;
  Vec3D *structX = (distIntersector.getStructPosition()).data();
  int ntri;
  MyTriangle *myTris;

  // build the KDTree
  if(useScope) {
    ntri = scope.size();
    if(ntri==0) {
      tId = -1;
      return;
    }
    myTris = new MyTriangle[ntri];
    int count = 0;
    set<int>::iterator it;
    for(it=scope.begin(); it!=scope.end(); it++)
      myTris[count++] = MyTriangle(*it, distIntersector.getStructPosition(), triNodes[*it]);
  } else {
    ntri = distIntersector.getNumStructElems();
    myTris = new MyTriangle[ntri];
    for(int i = 0; i < ntri; ++i)
      myTris[i] = MyTriangle(i, distIntersector.getStructPosition(), triNodes[i]);
  }
  KDTree<MyTriangle> structureTree(ntri, myTris);

  // clear scope for refill 
  scope.clear();

  // find candidates
  ClosestTriangle closestTriangle(triNodes, structX, distIntersector.triNorms, distIntersector.node2node, distIntersector.node2elem);
  int nMaxCand = 500;
  MyTriangle *candidates = new MyTriangle[nMaxCand];
  for(int i = 0; i < X.size(); ++i) {
/*    double bMin[3], bMax[3], expansion_factor = 1.0; //bounding boxes
    double extra = std::max(std::max(boxMax[i][0]-boxMin[i][0], boxMax[i][1]-boxMin[i][1]), boxMax[i][2]-boxMin[i][2]);
    extra *= (0.5*expansion_factor);
    for(int k=0; k<3; k++) {
      bMin[k] = boxMin[i][k] - extra;
      bMax[k] = boxMax[i][k] + extra;
    } 
*/    
    int nFound = structureTree.findCandidatesInBox(/*bMin*/boxMin[i], /*bMax*/boxMax[i], candidates, nMaxCand);
    if(nFound > nMaxCand) {
//      std::cerr << "For Fluid node " << locToGlobNodeMap[i]+1 << ", number of candidates: " << nFound << std::endl;
      nMaxCand = nFound;
      delete [] candidates;
      candidates = new MyTriangle[nMaxCand];
      structureTree.findCandidatesInBox(/*bMin*/boxMin[i], /*bMax*/boxMax[i], candidates, nMaxCand);
    }
    closestTriangle.start(X[i]);
    for(int j = 0; j < nFound; ++j) {
      int myId = candidates[j].trId();
      scope.insert(myId);
      addToPackage(i, myId);    
      closestTriangle.checkTriangle(myId);
/*      //debug
      Vec3D Xi(X[i][0], X[i][1], X[i][2]);
      Vec3D Tester(-1.333232e+00, 6.207416e-01, 1.327290e-01);
      if((Xi-Tester).norm()<1e-5)
        fprintf(stderr,"debug: fluid node = %d, TrId %d, n1/n2 = %d/%d, mode = %d, minDist = %e.\n", locToGlobNodeMap[i]+1, myId+1, 
                closestTriangle.n1+1, closestTriangle.n2+1, closestTriangle.mode, closestTriangle.minDist);
      //end of debug
*/    }
    
    if(nFound <= 0) {
      tId[i] = -1;
    }
    else {
      if(!closestTriangle.fail) {
        tId[i] = closestTriangle.bestTriangle();
        dist[i] = closestTriangle.signedDistance();
      } else {
        tId[i] = -1;
      }
    }
  }

  delete [] candidates;
  delete [] myTris;
}

//----------------------------------------------------------------------------

void IntersectorFRG::computeFirstLayerNodeStatus(Vec<int> tId, Vec<double> dist)
{
  for (int i=0; i<tId.size(); i++) {
    if (tId[i]<0)
      continue;
    status[i] = (dist[i]>=0.0) ? INSIDE : OUTSIDE;
    nFirstLayer++;
  }
}

//----------------------------------------------------------------------------

int IntersectorFRG::findNewSeedsAfterMerging(SVec<int,2>& status_and_weight, int& nUndecided)
{
  int numNewSeeds = 0;
  for(int i=0; i<status.size(); i++){
    if(status[i]!=UNDECIDED || status_and_weight[i][1]==0)
      continue; //already decided || didn't get anything from neighbours
    status[i] = status_and_weight[i][0] / status_and_weight[i][1] - 1;//back to normal conventional
    nUndecided--;
    numNewSeeds++;
  }
  return numNewSeeds;
}

//----------------------------------------------------------------------------

int IntersectorFRG::findNewSeedsAfterMerging(Vec<int>& status_temp, Vec<bool>& poly, int& nUndecided)
{
  int numNodes = status.size();
  int myStatus, numNewSeeds = 0;

  for(int i=0; i<numNodes; i++){
    if(status_temp[i]==status[i] || poly[i])
      continue; // inside or on the boundary but the neighbor hasn't decided it, or lies on n>2 subdomains.

    myStatus = status_temp[i];
//    if(myStatus!=INSIDE && myStatus!=OUTSIDE){ 
//      fprintf(stderr,"horror!\n"); exit(-1);}
    if (status[i]==UNDECIDED) {
      status[i] = myStatus;
      nUndecided--;
      numNewSeeds++;
    }
//    else if (status[i]!=myStatus)
//      fprintf(stderr,"ERROR: Node %d got different statuses from different subdomains.\n", locToGlobNodeMap[i]+1);
  }
  return numNewSeeds;
}

//----------------------------------------------------------------------------

int IntersectorFRG::findSeedsByPoints(SubDomain& sub, SVec<double,3>& X, list<pair<Vec3D,int> > P, int& nUndecided)
{
  int *myNodes, nSeeds = 0;
  list<pair<Vec3D,int> >::iterator iP;

  for(int iElem=0; iElem<sub.numElems(); iElem++) {
    myNodes = sub.getElemNodeNum(iElem);
    if(status[myNodes[0]]==INSIDE && status[myNodes[1]]==INSIDE &&
       status[myNodes[2]]==INSIDE && status[myNodes[3]]==INSIDE)//this tet is inside farfield fluid
      continue;

    for(iP=P.begin(); iP!=P.end(); iP++)
      if(sub.isINodeinITet(iP->first, iElem, X)) 
        for(int i=0; i<4; i++)
          if(status[myNodes[i]]==UNDECIDED) {
            status[myNodes[i]] = iP->second;
            nUndecided--;
            nSeeds++;
          }
  }

  return nSeeds;
}

//----------------------------------------------------------------------------

void IntersectorFRG::noCheckFloodFill(SubDomain& sub, int& nUndecided)
{
  int numNodes = status.size();
  Connectivity &nToN = *(sub.getNodeToNode());

  // Propogate status to the entire subdomain.
  // List is used as a FIFO queue
  Vec<int> seed(numNodes); //list of decided nodes
  // Look for a start point
  // lead: In pointer
  // next: out pointer (of FIFO)
  int next = 0, lead = 0;
  for(int i = 0; i < numNodes; ++i)
    if(status[i] != UNDECIDED && status[i] != INSIDE) {
      seed[lead++] = i;
    }

  while(next < lead) { //still have seeds not used
    int cur = seed[next++];
    int curStatus = status[cur];
    for(int i = 0; i < nToN.num(cur); ++i)
      if(status[nToN[cur][i]] == UNDECIDED) {
        status[nToN[cur][i]] = curStatus;
        seed[lead++] = nToN[cur][i];
        nUndecided--;
      } 
  }
}

//----------------------------------------------------------------------------

int IntersectorFRG::findIntersections(SVec<double,3>&X, bool useScope)
{
  int error = 0;

  int (*ptr)[2] = edges.getPtr();
  const double TOL = 1.0e-4;
  int MAX_ITER = 20;
  int max_iter = 0;

  bool GlobPhysBAMUpdated = false;

  for (int l=0; l<edges.size(); l++) {
    int p = ptr[l][0], q = ptr[l][1];
    if(status[p]==status[q]) continue;
    edge_intersects[l] = true;

    //now need to find an intersection for this edge .
    Vec3D xp(X[p]), xq(X[q]);
    Vec3D dir = xq - xp;
    Vec3D xpPrime, xqPrime;

    IntersectionResult<double> res1;

    for (int iter=0; iter<MAX_ITER; iter++) {
      double coeff = iter*iter*TOL;
      Vec3D xpPrime = xp - coeff*dir;
      Vec3D xqPrime = xq + coeff*dir;
      VECTOR<double,3> xyz1(xpPrime[0],xpPrime[1],xpPrime[2]), 
                       xyz2(xqPrime[0],xqPrime[1],xqPrime[2]);

      if(useScope) {
        res1 = physInterface->Intersect(xyz1, xyz2, coeff*dir.norm());
        if(res1.triangleID>0)
          res1.triangleID = iscope[res1.triangleID-1] + 1;
      } else 
        res1 = distIntersector.getInterface().Intersect(xyz1,xyz2, coeff*dir.norm());
      // the triangle Id stored in edgeRes starts from 1, i.e. using PhysBAM convention.

      if (res1.triangleID>0 /*|| edgeRes(2).y.triangleID>0*/) {
        CrossingEdgeRes[l] = res1;
        ReverseCrossingEdgeRes[l] = res1; //TODO: redundent!
        if (iter>max_iter) max_iter = iter;
        break;
      }
    }

    if (res1.triangleID<0 && useScope) { // try global intersector as a 'fail-safe'.
      if(!GlobPhysBAMUpdated) {
        distIntersector.updatePhysBAMInterface();
        GlobPhysBAMUpdated = true;
//        fprintf(stderr,"Now using the global intersector...\n");
      }

      for (int iter=0; iter<MAX_ITER; iter++) {
        double coeff = iter*iter*TOL;
        Vec3D xpPrime = xp - coeff*dir;
        Vec3D xqPrime = xq + coeff*dir;
        VECTOR<double,3> xyz1(xpPrime[0],xpPrime[1],xpPrime[2]), 
                         xyz2(xqPrime[0],xqPrime[1],xqPrime[2]);
        res1 = distIntersector.getInterface().Intersect(xyz1,xyz2, coeff*dir.norm()); 

        if (res1.triangleID>0 /*|| edgeRes(2).y.triangleID>0*/) {
          CrossingEdgeRes[l] = res1;
          ReverseCrossingEdgeRes[l] = res1; //TODO: redundent!
          if (iter>max_iter) max_iter = iter;
          break;
        }
      }
 
      if(res1.triangleID<0) {
         error++;
//         fprintf(stderr,"WARNING: No intersection between node %d(status = %d, status0 = %d) and %d(status = %d, status0 = %d). \n",
//                        locToGlobNodeMap[p]+1,status[p], status0[p], locToGlobNodeMap[q]+1,status[q], status0[q]);
         if(status[p]!=status0[p] || status[q]!=status0[q]) {
           status[p] = status0[p];
           status[q] = status0[q];
         } else
           status[p] = status[q] = 0; //TODO: this is not a fix!
      }
    }
  }

  return error;
}

//----------------------------------------------------------------------------

void IntersectorFRG::addToPackage(int iNode, int trID)
{
  int nSub = nodeToSubD.num(iNode);
  for(int iSub=0; iSub<nSub; iSub++) {
    if(nodeToSubD[iNode][iSub]==globIndex) 
      continue;
    package[sub2pack[nodeToSubD[iNode][iSub]]].insert(trID);  
  }
}

//----------------------------------------------------------------------------

LevelSetResult
IntersectorFRG::getLevelSetDataAtEdgeCenter(double t, int l, bool i_less_j, double *Xr, double *Xg) {
  if (!edge_intersects[l]) {
    int (*ptr)[2] = edges.getPtr();
    int i=i_less_j ? ptr[l][0] : ptr[l][1],
        j=i_less_j ? ptr[l][1] : ptr[l][0];
    fprintf(stderr,"There is no intersection between node %d(status:%d) and %d(status:%d)! Abort...\n",
                   locToGlobNodeMap[i]+1, status[i], locToGlobNodeMap[j]+1, status[j]);
    exit(-1);
  }

  IntersectionResult<double> result; //need to determine which result to choose.
  IntersectionResult<double> rij = CrossingEdgeRes[l];
  IntersectionResult<double> rji = ReverseCrossingEdgeRes[l];
  double alpha0 = 0.0;

  if (rij.triangleID>0 && rji.triangleID>0 && rij.triangleID==rji.triangleID) {
    result = rij;
    alpha0 = i_less_j ? result.alpha : 1.0-result.alpha;
  }  
  else if (rij.triangleID>0 && rji.triangleID>0 && rij.triangleID!=rji.triangleID) {
    result = i_less_j ? rij : rji;
    alpha0 = result.alpha;
  }
  else if (rij.triangleID>0 && rji.triangleID<0) {
    result = rij;
    alpha0 = i_less_j ? result.alpha : 1.0-result.alpha;
  }
  else if (rij.triangleID<0 && rji.triangleID>0) {
    result = rji;
    alpha0 = i_less_j ? 1.0-result.alpha : result.alpha;
  }
  else { //we really have no intersection for this edge!
    int (*ptr)[2] = edges.getPtr();
    int i=i_less_j ? ptr[l][0] : ptr[l][1],
        j=i_less_j ? ptr[l][1] : ptr[l][0];
    fprintf(stderr,"ERROR: intersection between node %d and node %d can not be detected.\n", locToGlobNodeMap[i]+1, locToGlobNodeMap[j]+1);
  }

  int trueTriangleID = result.triangleID-1;

  LevelSetResult lsRes;
  lsRes.alpha = alpha0;
  lsRes.xi[0] = result.zeta[0];
  lsRes.xi[1] = result.zeta[1];
  lsRes.xi[2] = 1-result.zeta[0]-result.zeta[1];
  lsRes.trNodes[0] = distIntersector.stElem[trueTriangleID][0];
  lsRes.trNodes[1] = distIntersector.stElem[trueTriangleID][1];
  lsRes.trNodes[2] = distIntersector.stElem[trueTriangleID][2];
  lsRes.normVel = lsRes.xi[0]*distIntersector.Xsdot[lsRes.trNodes[0]]
                + lsRes.xi[1]*distIntersector.Xsdot[lsRes.trNodes[1]]
                + lsRes.xi[2]*distIntersector.Xsdot[lsRes.trNodes[2]]; 
  lsRes.porosity = distIntersector.porosity[trueTriangleID];

  if(!distIntersector.interpolatedNormal) {

    lsRes.gradPhi = distIntersector.getSurfaceNorm(trueTriangleID);

    if(distIntersector.with_sensitivity && Xr != 0 && Xg != 0){
      Vec3D *XX  = (distIntersector.getStructPosition()).data();

      Vec3D dNdS;  
      derivativeOFnormal(XX[lsRes.trNodes[0]], 
			 XX[lsRes.trNodes[1]], 
			 XX[lsRes.trNodes[2]],
			 distIntersector.dXdSb[lsRes.trNodes[0]],
			 distIntersector.dXdSb[lsRes.trNodes[1]],
			 distIntersector.dXdSb[lsRes.trNodes[2]], dNdS);
      lsRes.dnds = dNdS;

      double dads;
      dads = derivativeOFalpha(XX[lsRes.trNodes[0]], 
                               XX[lsRes.trNodes[1]], 
       			       XX[lsRes.trNodes[2]],
       			       distIntersector.dXdSb[lsRes.trNodes[0]],
			       distIntersector.dXdSb[lsRes.trNodes[1]],
			       distIntersector.dXdSb[lsRes.trNodes[2]],
			       Xr, Xg);
      lsRes.dads = dads;
    }

  } else { //use nodal normals.
    Vec3D ns0 = distIntersector.getNodalNorm(lsRes.trNodes[0]);
    Vec3D ns1 = distIntersector.getNodalNorm(lsRes.trNodes[1]);
    Vec3D ns2 = distIntersector.getNodalNorm(lsRes.trNodes[2]);
    lsRes.gradPhi = lsRes.xi[0]*ns0 + lsRes.xi[1]*ns1 + lsRes.xi[2]*ns2;
    lsRes.gradPhi /= lsRes.gradPhi.norm();
  }

  return lsRes;
}

//----------------------------------------------------------------------------

void IntersectorFRG::projection(Vec3D x0, int tria, double& xi1, double& xi2, double& dist)
{
  int iA = distIntersector.stElem[tria][0];
  int iB = distIntersector.stElem[tria][1];
  int iC = distIntersector.stElem[tria][2];
  Vec3D xA = distIntersector.Xs[iA];
  Vec3D xB = distIntersector.Xs[iB];
  Vec3D xC = distIntersector.Xs[iC];

  Vec3D ABC = 0.5*(xB-xA)^(xC-xA);
  double areaABC = ABC.norm();
  Vec3D dir = 1.0/areaABC*ABC;

  //calculate the projection.
  dist = (x0-xA)*dir;
  Vec3D xp = x0 - dist*dir;

  //calculate barycentric coords.
  double areaPBC = (0.5*(xB-xp)^(xC-xp))*dir;
  double areaPCA = (0.5*(xC-xp)^(xA-xp))*dir;
  xi1 = areaPBC/areaABC;
  xi2 = areaPCA/areaABC;

}

//----------------------------------------------------------------------------

double IntersectorFRG::isPointOnSurface(Vec3D pt, int N1, int N2, int N3) 
{
  Vec<Vec3D> &solidX = distIntersector.getStructPosition();
  Vec3D X1 = solidX[N1];
  Vec3D X2 = solidX[N2];
  Vec3D X3 = solidX[N3];

  Vec3D normal = (X2-X1)^(X3-X1);
  normal /=  normal.norm();

  return fabs((pt-X1)*normal);
}

//----------------------------------------------------------------------------
void IntersectorFRG::derivativeOFnormal(Vec3D  xA, Vec3D  xB, Vec3D  xC, 
		   		        Vec3D dxA, Vec3D dxB, Vec3D dxC, 
				        Vec3D &dnds){

  Vec3D V1, V2, Vp, dV1, dV2, dVp; 
  double sp, dsp;

  V1 = (xC - xA);
  V2 = (xB - xA);
    
  dV1 = dxC - dxA;
  dV2 = dxB - dxA;

  Vp = -(V1^V2);

  dVp = -(dV1^V2) - (V1^dV2);

  double d = Vp*Vp;

  double tmp = dVp*Vp;

  dnds = (1.0/sqrt(d)) * (dVp) - (Vp*tmp)/pow(d, 3.0/2.0);

}

//----------------------------------------------------------------------------

double IntersectorFRG::derivativeOFalpha(Vec3D  xA, Vec3D  xB, Vec3D  xC, 
		  		         Vec3D dxA, Vec3D dxB, Vec3D dxC, 
					 Vec3D  X1, Vec3D  X2){

  Vec3D  u =  xB -  xA;
  Vec3D du = dxB - dxA;

  Vec3D  v =  xC -  xA;
  Vec3D dv = dxC - dxA;

  Vec3D  n = u ^ v;
  Vec3D dn = (du ^ v) + (u ^ dv);
 
  Vec3D d = X1 - X2;

  Vec3D  w = X2 - xA;
  Vec3D dw =    -dxA;

  double tmp1 = (dn*w) + (n*dw);
  double tmp2 = n*d;
  double tmp3 = n*w;
  double tmp4 = dn*d;
  
  double dalpha_ds = -(tmp1/tmp2) + tmp3*tmp4/(tmp2*tmp2);

  return dalpha_ds;

}

//----------------------------------------------------------------------------

void 
DistIntersectorFRG::updateXb(double epsilon){

  for(int i=0; i<numStNodes; ++i){

    for(int j=0; j<3; ++j){
      //Xs0[i][j]
      Xs[i][j]     = Xs0[i][j] + dXdSb[i][j]*epsilon;
      Xs_n[i][j]   = Xs[i][j];
      Xs_np1[i][j] = Xs_n[i][j];
      Xsdot[i][j]  = 0.0;
    }
  }

  initializePhysBAM();
  updatebc();

}

//----------------------------------------------------------------------------

void IntersectorFRG::xWallWithSI(int n, Vec3D &xWall)
{

	double l1 =  xi_SI[n];
	double l2 = eta_SI[n];

	int tria  = TriID_SI[n];

	Vec3D X0 = distIntersector.Xs[distIntersector.stElem[tria][0]];
	Vec3D X1 = distIntersector.Xs[distIntersector.stElem[tria][1]];
	Vec3D X2 = distIntersector.Xs[distIntersector.stElem[tria][2]];

	for(int k=0; k<3; ++k) 
		xWall[k] = X2[k] + l1*(X0[k] - X2[k]) + l2*(X1[k] - X2[k]);

}

//----------------------------------------------------------------------------

void IntersectorFRG::vWallWithSI(int n, Vec3D &vWall)
{
  
	double l1 =  xi_SI[n];
	double l2 = eta_SI[n];

	int tria  = TriID_SI[n];

	Vec3D v0 = distIntersector.Xsdot[distIntersector.stElem[tria][0]];
	Vec3D v1 = distIntersector.Xsdot[distIntersector.stElem[tria][1]];
	Vec3D v2 = distIntersector.Xsdot[distIntersector.stElem[tria][2]];

	for(int k=0; k<3; ++k) 
		vWall[k] = v2[k] + l1*(v0[k] - v2[k]) + l2*(v1[k] - v2[k]);

}

//----------------------------------------------------------------------------

bool IntersectorFRG::xWallNode(int i, Vec3D &xWall)
{

	double l1 =  xi_node[i];
	double l2 = eta_node[i];

	int tria  = TriID_node[i];

	if(tria < 0) return false;
	
	Vec3D X0 = distIntersector.Xs[distIntersector.stElem[tria][0]];
	Vec3D X1 = distIntersector.Xs[distIntersector.stElem[tria][1]];
	Vec3D X2 = distIntersector.Xs[distIntersector.stElem[tria][2]];

	for(int k=0; k<3; ++k) 
		xWall[k] = X2[k] + l1*(X0[k] - X2[k]) + l2*(X1[k] - X2[k]);

	return true;

}

//----------------------------------------------------------------------------

bool IntersectorFRG::vWallNode(int i, Vec3D &vWall)
{
  
	double l1 =  xi_node[i];
	double l2 = eta_node[i];

	int tria  = TriID_node[i];

	if(tria < 0) return false;

	Vec3D v0 = distIntersector.Xsdot[distIntersector.stElem[tria][0]];
	Vec3D v1 = distIntersector.Xsdot[distIntersector.stElem[tria][1]];
	Vec3D v2 = distIntersector.Xsdot[distIntersector.stElem[tria][2]];

	for(int k=0; k<3; ++k) 
		vWall[k] = v2[k] + l1*(v0[k] - v2[k]) + l2*(v1[k] - v2[k]);

	return true;
}
