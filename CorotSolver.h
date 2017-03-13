#ifndef _COROT_SOLVER_H_
#define _COROT_SOLVER_H_

#include <IoData.h>
#include <DistVector.h>

class MatchNodeSet;
class Domain;
class Communicator;
class BCApplier;

//------------------------------------------------------------------------------

class CorotSolver {

  int numLocSub;

  int *nInterfNd;  // int nInerfNd[3] ==> number of nodes for sub 3
  int **interfNd;  // int *interfNd[3] ==> array of nodes for sub 3
  int *nInfNd;     // number of nodes at infinity per sub
  int **infNd;
  int *nInternalNd;
  int **internalNd;

  bool locAllocGap;
  
  double (**gapVec)[3]; // gap vectors on structure interface

  double cg0[3];
  double cgN[3];

  double R[3][3];

  DistSVec<double,3> X0;
  DistSVec<double,3> Xtilde;

  Domain *domain;
  Communicator *com;

  enum SymmetryAxis {NONE, AXIS_X, AXIS_Y, AXIS_Z} SymAxis;

private:

  void computeCG(DistSVec<double,3> &, double [3]);
  void computeRotGradAndJac(DistSVec<double,3> &, double [3][3], 
			    double [3], double [3], double [3][3]);
  void computeRotMat(double *, double [3][3]);
  double computeSubGradAndJac(SVec<double,3> &, double [3][3], double [3],
			      double [3], double [3][3], int);
  void rotLocVec(double [3][3], double [3]);
  void solveDeltaRot(DistSVec<double,3> &, double [3][3], double [3]);
  void solveRotMat(double [3][3], double [3]);
  void computeInfNodeRot(double [3][3], DistSVec<double, 3> &, double [3], double [3]);
  void computeNodeRot(double [3][3], DistSVec<double,3> &, 
		      DistSVec<double,3> &, double [3]);

  void printRotMat(double mat[3][3]);

public:

  CorotSolver(DefoMeshMotionData &,  MatchNodeSet **, Domain *);
  ~CorotSolver();

  void solve(DistSVec<double,3> &dX, DistSVec<double,3> &X, BCApplier* meshMotionBCs=0);
  // setup computes the rotation and CG that will fit best for restarting
  void setup(DistSVec<double,3> &);
		
};

//------------------------------------------------------------------------------

#endif
