#ifndef _ELEMTET_H_
#define _ELEMTET_H_

#include <Elem.h>

template<class Scalar, int dim, int dim2> class RectangularSparseMat;

//------------------------------------------------------------------------------
/** Tetrahedral element */
class ElemTet : public ElemDummy {

  static const double third;
  static const double fourth;
  static const double sixth;
  static const int edgeEndTet[6][2];
  static const int edgeFaceTet[6][2];
  static const int faceDefTet[4][3];

  int nodeNumTet[4];
  int edgeNumTet[6];

protected:
  void *getWrapper_dim(GenElemHelper_dim *h, 
		       int size, char *memorySpace) {
    return h->forClassTet(this,size,memorySpace);
  }
  void *getWrapper_Scalar_dim_neq(GenElemHelper_Scalar_dim_neq *h, 
				  int size, char *memorySpace) {
    return h->forClassTet(this,size,memorySpace);
  }

  void *getWrapper_dim_obj(GenElemHelper_dim_obj *h, 
				  int size, char *memorySpace) {

    return h->forClassTet(this,size,memorySpace);
  }

public:

  ElemTet() { volume_id = 0; }
  ~ElemTet() {}

  int* nodeNum() { return nodeNumTet; }
  int& nodeNum(int i) { return nodeNumTet[i]; }
  int& edgeNum(int i) { return edgeNumTet[i]; }
  int  edgeEnd(int i, int k) { return edgeEndTet[i][k]; }
  int  edgeFace(int i, int k) { return edgeFaceTet[i][k]; }
  int  faceDef(int i, int k) { return faceDefTet[i][k]; }
  int  faceNnd(int i) { return 3; }
  Type type() { return Elem::TET; }

  bool isPointInside(SVec<double,3> &,const Vec3D&);

  // Number of nodes
  int numNodes() { return 4; }

  // Number of edges
  int numEdges() { return 6; }

  // Number of faces
  int numFaces() { return 4; }

  //--------------functions in TetCore.C

  double computeLongestEdge(SVec<double,3> &);
  double computeVolume(SVec<double,3> &);
  double computeControlVolumes(SVec<double,3> &, Vec<double> &);
  void printInvalidElement(int, double, int, int *, int *, SVec<double,3> &, SVec<double,3> &);
  void computeEdgeNormalsConfig(SVec<double,3> &Xconfig, SVec<double,3> &Xdot,
                                Vec<Vec3D> &edgeNorm, Vec<double> &edgeNormVel);
  void computeEdgeNormalsGCL1(SVec<double,3> &, SVec<double,3> &, SVec<double,3> &,
			      Vec<Vec3D> &, Vec<double> &);
  void computeEdgeNormalsEZGCL1(double, SVec<double,3> &, SVec<double,3> &, 
				Vec<Vec3D> &, Vec<double> &);
  void computeWeightsGalerkin(SVec<double,3> &, SVec<double,3> &, 
			      SVec<double,3> &, SVec<double,3> &);
  void computeEdgeWeightsGalerkin(SVec<double,3> &, SVec<double,9> &);
  double computeGradientP1Function(SVec<double,3> &nodes, double ngrad[4][3], double * = NULL);
  double computeGradientP1Function(Vec3D &A, Vec3D &B, Vec3D &C, Vec3D &D, 
                                   double nGrad[4][3]);
  void computeStiffAndForce(double *, double *, 
			    SVec<double, 3> &, SVec<double,3> &, double volStiff = 0.0);
  void computeStiffAndForceBallVertex(double *f, double *K, 
				    SVec<double, 3> &X, SVec<double,3> &X0, 
				    double volStiff = 0.0);
  void computeStiffAndForceLIN(double *, SVec<double,3> &, SVec<double,3> &);
  void computeStiffBallVertex(double *, SVec<double, 3> &X, SVec<double,3> &X0, double volStiff);
  void computeStiffTorsionSpring(double *, SVec<double, 3> &X, double volStiff);

  void computeTemp(double *[4], double [4], double);

  void computeTempGradient(double [4][3], double [4], double [3]);

  
  //-----functions in Tet.C

  template<class NodeMap>
  void renumberNodes(NodeMap &);

  template<int dim>
  void computeGalerkinTerm(FemEquationTerm *, SVec<double,3> &, Vec<double> &,
			   SVec<double,dim> &, SVec<double,dim> &,
			   Vec<GhostPoint<dim>*> *ghostPoints=0,LevelSetStructure *LSS=0);

  template<int dim>
  void computeGalerkinTerm_e(FemEquationTerm *, SVec<double,3> &, Vec<double> &,
									  SVec<double,dim> &, SVec<double,dim> &,
									  Vec<GhostPoint<dim>*> *ghostPoints=0,LevelSetStructure *LSS=0);

  template<int dim>
  double* setGhostOccludedValue(int i, SVec<double,3> &X, SVec<double,dim> &V, 
										  LevelSetStructure *LSS);

  template<int dim>
  void computeVMSLESTerm(VMSLESTerm *, SVec<double,dim> &, SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &);
                                                                                                                          
  template<int dim>
  void computeMBarAndM(DynamicVMSTerm *, SVec<double,dim> **, SVec<double,1> **, SVec<double,3> &, SVec<double,dim> &,
                       SVec<double,dim> &, SVec<double,dim> &);
                                                                                                                          
  template<int dim>
  void computeDynamicVMSTerm(DynamicVMSTerm *, SVec<double,dim> **, SVec<double,3> &,
                             SVec<double,dim> &, SVec<double,dim> &, Vec<double> &, Vec<double> &,
                             Vec<double> *, Vec<double> &);
                                                                                                                          
  template<int dim>
  void computeSmagorinskyLESTerm(SmagorinskyLESTerm *, SVec<double,3> &, SVec<double,dim> &V,
				 SVec<double,dim> &R,
			         Vec<GhostPoint<dim>*> *ghostPoints=0,LevelSetStructure *LSS=0);

  template<int dim>
  void computeSmagorinskyLESTerm_e(SmagorinskyLESTerm *, SVec<double,3> &, 
											  SVec<double,dim> &V,	SVec<double,dim> &R,
											  Vec<GhostPoint<dim>*> *ghostPoints=0, LevelSetStructure *LSS=0);

  template<int dim>
  void computeDynamicLESTerm(DynamicLESTerm *, SVec<double,2> &,
                             SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &,
			     Vec<GhostPoint<dim>*> *ghostPoints=0,LevelSetStructure *LSS=0);

  template<int dim>
  void computeDynamicLESTerm_e(DynamicLESTerm *, SVec<double,2> &,
										 SVec<double,3> &, SVec<double,dim> &, SVec<double,dim> &,
										 Vec<GhostPoint<dim>*> *ghostPoints=0,LevelSetStructure *LSS=0);

  template<int dim>
  void computeWaleLESTerm(WaleLESTerm *, SVec<double,3> &, SVec<double,dim> &V,
		          SVec<double,dim> &R,
			  Vec<GhostPoint<dim>*> *ghostPoints=0,LevelSetStructure *LSS=0);

  template<int dim>
  void computeWaleLESTerm_e(WaleLESTerm *, SVec<double,3> &, SVec<double,dim> &V,
									 SVec<double,dim> &R,
									 Vec<GhostPoint<dim>*> *ghostPoints=0,LevelSetStructure *LSS=0);

  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(FemEquationTerm *, SVec<double,3> &, 
				   Vec<double> &, Vec<double> &, 
				   SVec<double,dim> &, GenMat<Scalar,neq> &,
     				   Vec<GhostPoint<dim>*> *gp=0,LevelSetStructure *LSS=0);

  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm_e(FemEquationTerm *, SVec<double,3> &, 
												 Vec<double> &, Vec<double> &, 
												 SVec<double,dim> &, GenMat<Scalar,neq> &,
												 Vec<GhostPoint<dim>*> *gp=0, LevelSetStructure *LSS=0);

  template<int dim>
  void computeFaceGalerkinTerm(FemEquationTerm *, int [3], int, Vec3D &, 
			       SVec<double,3> &, Vec<double> &, double *, 
			       SVec<double,dim> &, SVec<double,dim> &);

  template<int dim, class Scalar, int neq>
  void computeFaceJacobianGalerkinTerm(FemEquationTerm *, int [3], int, Vec3D &, 
				       SVec<double,3> &, Vec<double> &, Vec<double> &, 
				       double *, SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<int dim>
  void computeP1Avg(SVec<double,dim> &, SVec<double,16> &, SVec<double,6> &, Vec<double> &,
                    SVec<double,8> &, SVec<double,3> &, SVec<double,dim> &, double, double,
     		    Vec<GhostPoint<dim>*> *gp=0,LevelSetStructure *LSS=0);

  template<int dim>
  void computeP1Avg_e(SVec<double,dim> &, SVec<double,16> &, SVec<double,6> &, Vec<double> &,
							 SVec<double,8> &, SVec<double,3> &, SVec<double,dim> &, double, double,
							 Vec<GhostPoint<dim>*> *gp=0,LevelSetStructure *LSS=0);

// Included (MB)
  double computeDerivativeOfVolume(SVec<double,3> &, SVec<double,3> &);
  void computeDerivativeOperatorsOfVolume(SVec<double,3> &, double [][3], double [][3], double [][3], double [][3]);

  double computeDerivativeOfControlVolumes(SVec<double,3> &, SVec<double,3> &, Vec<double> &);
  void computeDerivativeOperatorsOfControlVolumes(SVec<double,3> &, RectangularSparseMat<double,3,1> &);

  void computeDerivativeOfEdgeNormals(SVec<double,3> &, SVec<double,3> &, Vec<Vec3D> &, Vec<Vec3D> &, Vec<double> &, Vec<double> &);
  void computeDerivativeOperatorsOfEdgeNormals(SVec<double,3> &, RectangularSparseMat<double,3,3> &);

  void computeDerivativeOfWeightsGalerkin(SVec<double,3> &, SVec<double,3> &, SVec<double,3> &,
			      SVec<double,3> &, SVec<double,3> &);
  void computeDerivativeTransposeOfWeightsGalerkin(SVec<double,3> &, SVec<double,3> &, SVec<double,3> &,
			      SVec<double,3> &, SVec<double,3> &);

  double computeDerivativeOfGradientP1Function(SVec<double,3> &, SVec<double,3> &, double [4][3]);
  double computeDerivativeOfGradientP1Function2(SVec<double,3> &, SVec<double,3> &, double [4][3], double [4][3]);
  void computeDerivativeOperatorOfGradientP1Function(SVec<double,3> &, double [4][3], double [4][3][4][3], int [4] = NULL);
  void computeDerivativeTransposeOfGradientP1Function(SVec<double,3> &, double, double [4][3], double [4][3], SVec<double,3> &);

  template<int dim>
  void computeDerivativeOfGalerkinTerm(FemEquationTerm *, SVec<double,3> &, SVec<double,3> &, Vec<double> &,
			   SVec<double,dim> &, SVec<double,dim> &, double, SVec<double,dim> &);

  template<int dim>
  void computeDerivativeOperatorsOfGalerkinTerm(FemEquationTerm *, SVec<double,3> &, Vec<double> &,
                                                SVec<double,dim> &, RectangularSparseMat<double,3,dim> &);
  template<int dim>
  void computeDerivativeOfFaceGalerkinTerm(FemEquationTerm *, int [3], int, Vec3D &, Vec3D &,
			       SVec<double,3> &, SVec<double,3> &, Vec<double> &, double *, double *,
			       SVec<double,dim> &, SVec<double,dim> &, double, SVec<double,dim> &);

// Level Set Reinitialization

  template<int dim>
  void computeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                 SVec<double,dim> &ddx, SVec<double,dim> &ddy,
                                 SVec<double,dim> &ddz,
                                 SVec<double,dim> &Phi,SVec<double,1> &Psi);

  template<int dim>
  void recomputeDistanceCloseNodes(int lsdim, Vec<int> &Tag, SVec<double,3> &X,
                                 SVec<double,dim> &ddx, SVec<double,dim> &ddy,
                                 SVec<double,dim> &ddz,
                                 SVec<double,dim> &Phi,SVec<double,1> &Psi);

  template<int dim>
  void computeDistanceLevelNodes(int lsdim, Vec<int> &Tag, int level,
                                 SVec<double,3> &X, SVec<double,1> &Psi, SVec<double,dim> &Phi);
  template<int dim>
  void FastMarchingDistanceUpdate(int node,Vec<int> &Tag,int level,
                                 SVec<double,3> &X,SVec<double,dim> &d2wall);


  template<int dim, class Obj>
    void integrateFunction(Obj* obj,SVec<double,3> &X,SVec<double,dim>& V, void (Obj::*F)(int node, const double* loc,double* f),int npt);

  // X is the deformed nodal location vector
  template<int dim> 
  int interpolateSolution(SVec<double,3>& X, SVec<double,dim>& U, const Vec3D& loc, double sol[dim], LevelSetStructure* LSS,
                          Vec<GhostPoint<dim>*>* ghostPoints, VarFcn* varFcn);

  //--------------functions in ElemTetCore.C

  void computeBarycentricCoordinates(SVec<double,3>&X, const Vec3D& loc, double bary[3]);

private:

//Level Set Reinitialization
  double findRootPolynomialNewtonRaphson(double f1, double f2, double fp1, double fp2);
  int findRootPolynomialLaguerre(double f1, double f2, double fp1, double fp2, double &root);
  bool computeDistancePlusPhiToOppFace(double phi[3], Vec3D Y0,
                                       Vec3D Y1, Vec3D Y2, double &mini, bool show = false);
  bool computeDistancePlusPhiToEdges(double phi[3], Vec3D Y0,
                                     Vec3D Y1, Vec3D Y2, double &mini, bool show = false);
  bool computeDistancePlusPhiToVertices(double phi[3], Vec3D Y0,
                                        Vec3D Y1, Vec3D Y2, double &mini, bool show = false);
  bool computeDistancePlusPhiToEdge(double phi0, double phi1,
                                    Vec3D Y0, Vec3D Y1, double &mini, bool show = false);
  int computeDistanceToAll(double phi[3],Vec3D Y0,Vec3D Y1,Vec3D Y2, double &psi);


  void getVelocityAndGradient(double *v[4], double dp1dxj[4][3],
										double u[4][3], double dudxj[3][3]);

  void getTemperatureAndGradient(double *v[4], double dp1dxj[4][3], double R,
											double T[4],  double dtdxj[3]);

  void ComputeStrainAndStressTensor(double dudxj[3][3], 
												double S[3][3], double Pij[6]);
  

  //--------------functions in ElemTet.C

//Level Set Reinitialization
  template<int dim>
  int findLSIntersectionPoint(int lsdim, SVec<double,dim> &Phi, SVec<double,dim> &ddx,
                              SVec<double,dim> &ddy, SVec<double,dim> &ddz,
 			      SVec<double,3> &X,
                              int reorder[4], Vec3D P[4], int typeTracking);
  template<int dim>
  void findLSIntersectionPointLinear(int lsdim, SVec<double,dim> &Phi, SVec<double,dim> &ddx,
                                 SVec<double,dim> &ddy, SVec<double,dim> &ddz,
				 SVec<double,3> &X,
                                 int reorder[4], Vec3D P[4], int scenario);
  template<int dim>
  void findLSIntersectionPointGradient(int lsdim, SVec<double,dim> &Phi,  SVec<double,dim> &ddx,
                                 SVec<double,dim> &ddy, SVec<double,dim> &ddz,
				 SVec<double,3> &X,
                                 int reorder[4], Vec3D P[4], int scenario);
  template<int dim>
  int findLSIntersectionPointHermite(int lsdim, SVec<double,dim> &Phi,  SVec<double,dim> &ddx,
                                 SVec<double,dim> &ddy, SVec<double,dim> &ddz,
                                 SVec<double,3> &X,
                                 int reorder[4], Vec3D P[4], int scenario);

  template<int dim>
  void computeDistanceToInterface(int type, SVec<double,3> &X, int reorder[4],
                                  Vec3D P[4], SVec<double,dim> &Psi, Vec<int> &Tag);
  template<int dim>
  void recomputeDistanceToInterface(int type, SVec<double,3> &X, int reorder[4],
                                Vec3D P[4], SVec<double,dim> &Psi, Vec<int> &Tag);

  template<int dim>
  double computeDistancePlusPhi(int i, SVec<double,3> &X, SVec<double,dim> &Psi);

};

#ifdef TEMPLATE_FIX
#include <ElemTet.C>
#endif

#endif
