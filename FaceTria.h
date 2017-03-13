#ifndef _FACETRIA_H_
#define _FACETRIA_H_

#include <Face.h>

template<class Scalar, int dim, int dim2> class RectangularSparseMat;

//------------------------------------------------------------------------------

class FaceTria : public FaceDummy {

  static const double third;
  static const int edgeEndT[Face::MaxNumNd][2];

  int nodeNumT[3];
  int edgeNumT[3];
  int nodeNumTet[4];

  // This is for triangles only (6 Vec3D's = 3 positions at t_n and t_n+1... 
  // make general function (also need to modify function calls)
  template<int dim>
  void computeForce(ElemSet &,
		    PostFcn *, SVec<double,3> &, Vec<double> &, double *, 
		    SVec<double,dim> &, double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &,
		    double*gradP[3], int = 0);

  template<int dim>
  void computeForce(ExactRiemannSolver<dim>&, VarFcn*, Vec<Vec3D> &, Vec<double> &, ElemSet &,
                    PostFcn *, SVec<double,3> &, Vec<double> &, double *,
                    SVec<double,dim> &, double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &,
                    double*gradP[3], int = 0);


  template<int dim>
  void computeForceTransmitted(ElemSet &,
                    PostFcn *, SVec<double,3> &, Vec<double> &, double *,
                    SVec<double,dim> &, double *, Vec3D &, Vec3D &, Vec3D &, Vec3D &,
                    double*gradP[3], int = 0);

protected:
  void *getWrapper_dim(GenFaceHelper_dim *h, 
		       int size, char *memorySpace) {
    return h->forClassTria(this,size,memorySpace);
  }
  void *getWrapper_Scalar_dim_neq(GenFaceHelper_Scalar_dim_neq *h, 
				  int size, char *memorySpace) {
    return h->forClassTria(this,size,memorySpace);
  }

  int* nodeNum() { return nodeNumT; }
  int& nodeNum(int i) { return nodeNumT[i]; }
  int& edgeNum(int i) { return edgeNumT[i]; }
  int  edgeEnd(int i, int k) { return edgeEndT[i][k]; }

public:

  // Number of nodes
  int numNodes() { return 3; }

  int nodeNum(int i) const  { return nodeNumT[i]; }
  
  void setEdgeNum(int edge_id, int l) { edgeNumT[edge_id] = l; }

  // Get element type
  Type type() { return Face::TRIA; }

  // Number of normals to be stored for a triangle 
  // (the normals for all nodes equal, thus only 1)
  int numNorms() { return 1; }

  // Compute total normal and return in Vec3D
  void computeNormal(SVec<double,3> &, Vec3D &);

  // Compute subface normals and save then in Vec<Vec3D>
  void computeNormal(SVec<double,3> &, Vec<Vec3D> &);
  void computeNormalConfig(SVec<double,3> &Xconfig, SVec<double,3> &Xdot,
                           Vec<Vec3D> &faceNorm, Vec<double> &faceNormVel);
  void computeNormalGCL1(SVec<double,3> &, SVec<double,3> &, 
			 SVec<double,3> &, Vec<Vec3D> &, Vec<double> &);
  void computeNormalEZGCL1(double, SVec<double,3> &, SVec<double,3> &, 
			   Vec<Vec3D> &, Vec<double> &);

  // Get face total normal from Vec<Vec3D>
  Vec3D  getNormal(Vec<Vec3D> &);
  double getNormalVel(Vec<double> &);

  // Get subface i normal from Vec<Vec3D>
  Vec3D  getNormal(Vec<Vec3D> &, int);
  double getNormalVel(Vec<double> &, int);

  template<int dim>
  void computeNodalForce(ElemSet &, PostFcn *, SVec<double,3> &, 
			 Vec<double> &, double *, SVec<double,dim> &, double, 
			 SVec<double,3> &, double* gradP[3]);

  template<int dim>
  void computeNodalHeatPower(ElemSet &,PostFcn*, SVec<double,3>&, Vec<double>&, 
			     double*, SVec<double,dim>&, Vec<double>&);

  template<int dim>
  double computeHeatFluxes(ElemSet &,PostFcn*, SVec<double,3>&, Vec<double>&,
                             double*, SVec<double,dim>&);

  template<int dim>
  void computeNodalHeatFluxRelatedValues(ElemSet &, PostFcn*, SVec<double,3>&, Vec<double>&, double*,
                                           SVec<double,dim>&, Vec<double>&, Vec<double>&, bool);

  template<int dim>
  void computeForceAndMoment(ElemSet &, PostFcn *, SVec<double,3> &, Vec<double> &, 
			     double *, SVec<double,dim> &, Vec3D &, Vec3D &, Vec3D &, 
			     Vec3D &, Vec3D &, double* gradP[3], int, 
                             SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX,
                             Vec<double> *genCF);
  template<int dim>
  void computeForceAndMoment(ExactRiemannSolver<dim>&, VarFcn *, Vec<Vec3D> &, Vec<double> &,
                             ElemSet &, PostFcn *, SVec<double,3> &, Vec<double> &,
                             double *, SVec<double,dim> &, Vec3D &, Vec3D &, Vec3D &,
                             Vec3D &, Vec3D &, double* gradP[3], int,
                             SubVecSet< DistSVec<double,3>, SVec<double,3> > *mX,
                             Vec<double> *genCF);

  template<int dim>
  double computeInterfaceWork(ElemSet &, PostFcn*, SVec<double,3>&, Vec<double>&, 
			      double, double*, SVec<double,dim>&, double);

  template<int dim>
  void computeScalarQuantity(PostFcn::ScalarType, ElemSet &, PostFcn *, SVec<double,3> &, 
			     Vec<double> &, double *, SVec<double,dim> &, SVec<double,2> &);

  template<int dim>
  void computeGalerkinTerm(ElemSet &, FemEquationTerm *, SVec<double,3> &, 
			   Vec<double> &, double *, SVec<double,dim> &, SVec<double,dim> &,LevelSetStructure *LSS=0);
  
  template<int dim, class Scalar, int neq>
  void computeJacobianGalerkinTerm(ElemSet &, FemEquationTerm *, SVec<double,3> &, 
				   Vec<double> &, Vec<double> &, double *,
				   SVec<double,dim> &, GenMat<Scalar,neq> &);

  template<int dim>
  void computeForceDerivs(ElemSet &, VarFcn *, SVec<double,3> &, 
                          SVec<double,dim> &,SVec<double,dim> &, Vec<double> &, SVec<double,3> **);

  template<int dim>
  void computeForceCoefficients(PostFcn *, Vec3D &, ElemSet &, SVec<double,3> &, 
				SVec<double,dim> &,  Vec<double> &, SVec<double, dim> &, 
				double, Vec3D &, Vec3D &, Vec3D &,
				Vec3D &, double* gradP[3], VecSet< SVec<double,3> > *mX = 0, 
                                Vec<double> *genCF = 0);

  template<int dim>
  void computeFDerivs(ElemSet &, VarFcn *, SVec<double,3> &, SVec<double,dim> &, Vec3D (*));

// Included (MB)
  // Get face total normal derivative from Vec<Vec3D>
  Vec3D  getdNormal(Vec<Vec3D> &);
  double getdNormalVel(Vec<double> &);

  // Get subface i normal derivative from Vec<Vec3D>
  Vec3D  getdNormal(Vec<Vec3D> &, int);
  double getdNormalVel(Vec<double> &, int);

  template<int dim>
  void computeDerivativeOfForce(ElemSet &, PostFcn *, SVec<double,3> &, SVec<double,3> &,
                                Vec<double> &, double *, double *, SVec<double,dim> &,
                                SVec<double,dim> &, double [3], double *, Vec3D &, Vec3D &, 
                                Vec3D &, Vec3D &, double * gradP[3], double* dGradP[3], int = 0);

  template<int dim>
  void computeDerivativeOperatorsOfForce(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X, SVec<double,dim> &V, double *pin, double* gradP[3], int hydro,
                                         double dFi0dGradP[3][3], double dFi1dGradP[3][3], double dFi2dGradP[3][3],
                                         double dFi0dX0[3][3], double dFi0dX1[3][3], double dFi0dX2[3][3],
                                         double dFi1dX0[3][3], double dFi1dX1[3][3], double dFi1dX2[3][3],
                                         double dFi2dX0[3][3], double dFi2dX1[3][3], double dFi2dX2[3][3],
                                         double dFidS[3][3],
                                         double dFi0dV[3][5], double dFi1dV[3][5], double dFi2dV[3][5],
                                         double dFvdX[3][3][3], double dFvdXtet[3][4][3], double dFvdV[3][4][dim]);

  template<int dim>
  void computeDerivativeOfForceTransmitted(ElemSet &, PostFcn *, SVec<double,3> &, SVec<double,3> &,
                                           Vec<double> &, double *, double *, SVec<double,dim> &,
                                           SVec<double,dim> &, double [3], double *, Vec3D &, Vec3D &, 
                                           Vec3D &, Vec3D &, double* gradP[3], double* dGradP[3], int = 0);

  template<int dim>
  void computeDerivativeOperatorsOfForceTransmitted(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X, SVec<double,dim> &V,
                                                    double *pin, double* gradP[3], int hydro,
                                                    double dFi0dV[3][3][dim], double dFi0dGradP[3][3][3], double dFi0dX[3][3][3],
                                                    double dFi0dn[3], double dFi0dS[3][3],
                                                    double dFi1dV[3][3][dim], double dFi1dGradP[3][3][3], double dFi1dX[3][3][3],
                                                    double dFi1dn[3], double dFi1dS[3][3],
                                                    double dFi2dV[3][3][dim], double dFi2dGradP[3][3][3], double dFi2dX[3][3][3],
													double dFi2dn[3], double dFi2dS[3][3],
													double dFvdX[3][3][3], double dFvdXtet[3][4][3], double dFvdV[3][4][dim]);

  template<int dim>
  void computeDerivativeOfNodalForce(ElemSet &, PostFcn *, SVec<double,3> &, SVec<double,3> &,
                                     Vec<double> &, double *, double *,
                                     SVec<double,dim> &, SVec<double,dim> &,
                                     double, double [3], SVec<double,3> &, 
                                     double * gradP[3], double* dGradP[3]);

  template<int dim>
  void computeDerivativeOperatorsOfNodalForce(ElemSet &elems, PostFcn *postFcn, SVec<double,3> &X,
                          SVec<double,dim> &V, double pin, double* gradP[3], 
                          RectangularSparseMat<double,3,3> &dForcedX,
                          RectangularSparseMat<double,3,3> &dForcedGradP,
                          RectangularSparseMat<double,dim,3> &dForcedV,
                          RectangularSparseMat<double,3,3> &dForcedS);

  template<int dim>
  void computeDerivativeOfNodalHeatPower(ElemSet&, PostFcn*, SVec<double,3>&, SVec<double,3>&, Vec<double>&, 
                                         double*, double*, SVec<double,dim>&, SVec<double,dim>&, double [3], Vec<double>&);

  template<int dim>
  void computeDerivativeOfForceAndMoment(ElemSet &, PostFcn *, SVec<double,3> &, SVec<double,3> &,
                                         Vec<double> &, double *, double *,
                                         SVec<double,dim> &, SVec<double,dim> &, double [3],
                                         Vec3D &, Vec3D &, Vec3D &, Vec3D &, Vec3D &, 
                                         double * gradP[3], double* dGradP[3], int = 0);

  //  template<int dim>
  //  void computeDerivativeOfForceAndMoment2(ElemSet &, PostFcn *, SVec<double,3> &, SVec<double,3> &,
  //                                          Vec<double> &, double *, double *,
  //                                          SVec<double,dim> &, SVec<double,dim> &, double [3],
  //                                          Vec3D &, Vec3D &, Vec3D &, Vec3D &, Vec3D &,
  //                                          double * gradP[3], double* dGradP[3], int = 0);

  template<int dim>
  void computeDerivativeOperatorsOfForceAndMoment(ElemSet &, PostFcn *, SVec<double,3> &,
                                                  Vec<double> &, double *, SVec<double,dim> &,
                                                  Vec3D &, double * gradP[3], int,
                                                  RectangularSparseMat<double,3,3> &dFidGradP,
                                                  RectangularSparseMat<double,3,3> &dFidX,
                                                  RectangularSparseMat<double,dim,3> &dFidV,
                                                  RectangularSparseMat<double,3,3> &dFvdX,
                                                  RectangularSparseMat<double,dim,3> &dFvdV,
                                                  RectangularSparseMat<double,3,3> &dFidS,
                                                  RectangularSparseMat<double,3,3> &dMidGradP,
                                                  RectangularSparseMat<double,3,3> &dMidX,
                                                  RectangularSparseMat<double,dim,3> &dMidV,
                                                  RectangularSparseMat<double,3,3> &dMidS,
                                                  RectangularSparseMat<double,3,3> &dMvdX,
                                                  RectangularSparseMat<double,dim,3> &dMvdV);

  template<int dim>
  void computeDerivativeOfGalerkinTerm(ElemSet &, FemEquationTerm *, SVec<double,3> &, SVec<double,3> &,
                                       Vec<double> &, double *, double *, SVec<double,dim> &, 
                                       SVec<double,dim> &, double, SVec<double,dim> &);
  
  template<int dim>
  void computeBCsJacobianWallValues(ElemSet &, FemEquationTerm *, SVec<double,3> &, 
				    Vec<double> &, double *, double *,
				    SVec<double,dim> &);
				   
  void computeNormalAndDerivative(SVec<double,3> &, SVec<double,3> &, Vec3D &, Vec3D &);

  void computeDerivativeOfNormal(SVec<double,3> &, SVec<double,3> &, Vec3D &, Vec3D &, double &, double &);
  void computeDerivativeOperatorsOfNormal(int, SVec<double,3> &, RectangularSparseMat<double,3,3> &);
  void compute_dndX(SVec<double,3> &, double dndX[3][3][3]);

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <FaceTria.C>
#endif

#endif
