#include <FaceTria.h>

#include <RefVal.h>
#include <BcDef.h>
#include <Edge.h>
#include <Vector3D.h>
#include <Vector.h>
#include <BinFileHandler.h>

#include <cstdlib>
#include <cstdio>
#include <cmath>

#ifdef OLD_STL
#include <defalloc.h>
#include <algo.h>
#else
#include <algorithm>
using std::stable_sort;
using std::min;
using std::max;
using std::swap;
#endif

const double FaceTria::third = 1.0/3.0;
const int FaceTria::edgeEndT[Face::MaxNumNd][2] = { {0,1}, {1,2}, {2,0}, {-1,-1} };


//------------------------------------------------------------------------------
// Computation of the OUTWARD subface normals (only one, they are all equal)
void FaceTria::computeNormal(SVec<double,3> &X, Vec<Vec3D> &faceNorm)
{
  Vec3D x[3] = {X[ nodeNum(0) ], X[ nodeNum(1) ], X[ nodeNum(2) ]};

  faceNorm[normNum] = 0.5 * ((x[2] - x[0]) ^ (x[1] - x[0]));

}

//------------------------------------------------------------------------------
// Computation of the OUTWARD face normal
void FaceTria::computeNormal(SVec<double,3> &X, Vec3D &faceNorm)
{

  Vec3D x[3] = {X[ nodeNum(0) ], X[ nodeNum(1) ], X[ nodeNum(2) ]};

  faceNorm = 0.5 * ((x[2] - x[0]) ^ (x[1] - x[0]));

}

//------------------------------------------------------------------------------

void FaceTria::computeNormalConfig(SVec<double,3> &Xconfig, SVec<double,3> &Xdot,
				   Vec<Vec3D> &faceNorm, Vec<double> &faceNormVel)
{

  Vec3D x[3] = {Xconfig[ nodeNum(0) ], Xconfig[ nodeNum(1) ], Xconfig[ nodeNum(2) ]};
  Vec3D xdot[3] = {Xdot[ nodeNum(0) ], Xdot[ nodeNum(1) ], Xdot[ nodeNum(2) ]};

  Vec3D configFaceNorm = 0.5 * ((x[2] - x[0]) ^ (x[1] - x[0]));
  faceNorm[normNum] += configFaceNorm;
  faceNormVel[normNum] += third * (xdot[0] + xdot[1] + xdot[2]) * configFaceNorm;
  
}

//------------------------------------------------------------------------------
// Computation of the OUTWARD face normal
void FaceTria::computeNormalGCL1(SVec<double,3> &Xn, SVec<double,3> &Xnp1, 
				 SVec<double,3> &Xdot, Vec<Vec3D> &faceNorm, 
				 Vec<double> &faceNormVel)
{

  static double twelfth = 1.0/12.0;

  Vec3D x_n[3] = {Xn[ nodeNum(0) ], Xn[ nodeNum(1) ], Xn[ nodeNum(2) ]};
  Vec3D x_np1[3] = {Xnp1[ nodeNum(0) ], Xnp1[ nodeNum(1) ], Xnp1[ nodeNum(2) ]};
  Vec3D xdot[3] = {Xdot[ nodeNum(0) ], Xdot[ nodeNum(1) ], Xdot[ nodeNum(2) ]};

  Vec3D x01_n = x_n[1] - x_n[0];
  Vec3D x02_n = x_n[2] - x_n[0];

  Vec3D x01_np1 = x_np1[1] - x_np1[0];
  Vec3D x02_np1 = x_np1[2] - x_np1[0];

  faceNorm[normNum] = twelfth * (((2.0*x02_np1 + x02_n) ^ x01_np1) + 
				 ((2.0*x02_n + x02_np1) ^ x01_n));

  faceNormVel[normNum] = third * (xdot[0] + xdot[1] + xdot[2]) * faceNorm[normNum];

  //TODO: for testing only. If the HH farfield flux works, this should be implemented in a better way!
  faceCenter = third * (x_n[0] + x_n[1] + x_n[2]);
}

//------------------------------------------------------------------------------
// computation of the OUTWARD face normal

// Included (MB)
void FaceTria::computeDerivativeOfNormal(SVec<double,3> &X, SVec<double,3> &dX, Vec3D &faceNorm,
                            Vec3D &dFaceNorm, double &faceNormVel, double &dFaceNormVel)
{

  Vec3D x[3] = {X[ nodeNum(0) ], X[ nodeNum(1) ], X[ nodeNum(2) ]};

  Vec3D dx[3] = {dX[ nodeNum(0) ], dX[ nodeNum(1) ], dX[ nodeNum(2) ]};

//  faceNorm = 0.5 * ((x[2] - x[0]) ^ (x[1] - x[0]));

  dFaceNorm = 0.5 * (((dx[2] - dx[0]) ^ (x[1] - x[0])) + ((x[2] - x[0]) ^ (dx[1] - dx[0])));

//  faceNormVel = 0.0;
  dFaceNormVel = 0.0;

}

//------------------------------------------------------------------------------
// computation of the OUTWARD face normal

void FaceTria::compute_dndX(SVec<double,3> &X, double dFaceNormdX[3][3][3]) 
{

  Vec3D x[3] = {X[ nodeNum(0) ], X[ nodeNum(1) ], X[ nodeNum(2) ]};

  double dx100 = 0.5 * (x[1][0] - x[0][0]);
  double dx101 = 0.5 * (x[1][1] - x[0][1]);
  double dx102 = 0.5 * (x[1][2] - x[0][2]);
  double dx200 = 0.5 * (x[2][0] - x[0][0]);
  double dx201 = 0.5 * (x[2][1] - x[0][1]);
  double dx202 = 0.5 * (x[2][2] - x[0][2]);

  dFaceNormdX[0][0][1] = dx202-dx102;
  dFaceNormdX[0][0][2] = dx101-dx201;
  dFaceNormdX[0][1][0] = dx102-dx202;
  dFaceNormdX[0][1][2] = dx200-dx100;
  dFaceNormdX[0][2][0] = dx201-dx101;
  dFaceNormdX[0][2][1] = dx100-dx200;

  dFaceNormdX[1][0][1] = -dx202;
  dFaceNormdX[1][0][2] =  dx201;
  dFaceNormdX[1][1][0] =  dx202;
  dFaceNormdX[1][1][2] = -dx200;
  dFaceNormdX[1][2][0] = -dx201;
  dFaceNormdX[1][2][1] =  dx200;

  dFaceNormdX[2][0][1] =  dx102;
  dFaceNormdX[2][0][2] = -dx101;
  dFaceNormdX[2][1][0] = -dx102;
  dFaceNormdX[2][1][2] =  dx100;
  dFaceNormdX[2][2][0] =  dx101;
  dFaceNormdX[2][2][1] = -dx100;

}

// Included (YC)
void FaceTria::computeDerivativeOperatorsOfNormal(int faceNum, SVec<double,3> &X, RectangularSparseMat<double,3,3> &dFaceNormdX)
{

  Vec3D x[3] = {X[ nodeNum(0) ], X[ nodeNum(1) ], X[ nodeNum(2) ]};
/*
  dFaceNorm = 0.5 * (((dx[2] - dx[0]) ^ (x[1] - x[0])) + ((x[2] - x[0]) ^ (dx[1] - dx[0])));

  dFaceNorm[0] = 0.5 * (  (dx[2][1] - dx[0][1])*(x[1][2] - x[0][2]) - (dx[2][2] - dx[0][2])*(x[1][1] - x[0][1]) 
                        + (x[2][1] - x[0][1])*(dx[1][2] - dx[0][2]) - (x[2][2] - x[0][2])*(dx[1][1] - dx[0][1]) )
  dFaceNorm[1] = 0.5 * (  (dx[2][2] - dx[0][2])*(x[1][0] - x[0][0]) - (dx[2][0] - dx[0][0])*(x[1][2] - x[0][2]) 
                        + (x[2][2] - x[0][2])*(dx[1][0] - dx[0][0]) - (x[2][0] - x[0][0])*(dx[1][2] - dx[0][2]) )
  dFaceNorm[2] = 0.5 * (  (dx[2][0] - dx[0][0])*(x[1][1] - x[0][1]) - (dx[2][1] - dx[0][1])*(x[1][0] - x[0][0]) 
                        + (x[2][0] - x[0][0])*(dx[1][1] - dx[0][1]) - (x[2][1] - x[0][1])*(dx[1][0] - dx[0][0]) )
*/
  double dx100 = 0.5 * (x[1][0] - x[0][0]);
  double dx101 = 0.5 * (x[1][1] - x[0][1]);
  double dx102 = 0.5 * (x[1][2] - x[0][2]);
  double dx200 = 0.5 * (x[2][0] - x[0][0]);
  double dx201 = 0.5 * (x[2][1] - x[0][1]);
  double dx202 = 0.5 * (x[2][2] - x[0][2]);
/*
  dFaceNormdX = 0.5 * [  (dx202-dx102)*dx[0][1] + (dx101-dx201)*dx[0][2] - dx202*dx[1][1] + dx201*dx[1][2] + dx102*dx[2][1] - dx101*dx[2][2] ]
                      [  (dx102-dx202)*dx[0][0] + (dx200-dx100)*dx[0][2] + dx202*dx[1][0] - dx200*dx[1][2] - dx102*dx[2][0] + dx100*dx[2][2] ]
                      [  (dx201-dx101)*dx[0][0] + (dx100-dx200)*dx[0][1] - dx201*dx[1][0] + dx200*dx[1][1] + dx101*dx[2][0] - dx100*dx[2][1] ]

  dFaceNormdX = 0.5 * [  0               (dx202-dx102)  (dx101-dx201)      0  -dx202   dx201       0     dx102     -dx101 ] [ dx[0][0] ]
                      [  (dx102-dx202)               0  (dx200-dx100)  dx202       0  -dx200  -dx102         0      dx100 ] [ dx[0][1] ]
                      [  (dx201-dx101)   (dx100-dx200)              0 -dx201   dx200       0   dx101    -dx100          0 ] [ dx[0][2] ]
                                                                                                                            [ dx[1][0] ]
                                                                                                                            [ dx[1][1] ]
                                                                                                                            [ dx[1][2] ]
                                                                                                                            [ dx[2][0] ]
                                                                                                                            [ dx[2][1] ]
                                                                                                                            [ dx[2][2] ] 
*/
  double dFaceNormdX0[3][3] = {0}, dFaceNormdX1[3][3] = {0}, dFaceNormdX2[3][3] = {0};
  dFaceNormdX0[0][1] = dx202-dx102;
  dFaceNormdX0[0][2] = dx101-dx201;
  dFaceNormdX0[1][0] = dx102-dx202;
  dFaceNormdX0[1][2] = dx200-dx100;
  dFaceNormdX0[2][0] = dx201-dx101;
  dFaceNormdX0[2][1] = dx100-dx200;

  dFaceNormdX1[0][1] = -dx202;
  dFaceNormdX1[0][2] =  dx201;
  dFaceNormdX1[1][0] =  dx202;
  dFaceNormdX1[1][2] = -dx200;
  dFaceNormdX1[2][0] = -dx201;
  dFaceNormdX1[2][1] =  dx200;

  dFaceNormdX2[0][1] =  dx102;
  dFaceNormdX2[0][2] = -dx101;
  dFaceNormdX2[1][0] = -dx102;
  dFaceNormdX2[1][2] =  dx100;
  dFaceNormdX2[2][0] =  dx101;
  dFaceNormdX2[2][1] = -dx100;
  
//  fprintf(stderr," ... in dFaceNormdX, faceNum = %d, nodeNum(0) = %d\n", faceNum, nodeNum(0));
//  fprintf(stderr," ... in dFaceNormdX, faceNum = %d, nodeNum(1) = %d\n", faceNum, nodeNum(1));
//  fprintf(stderr," ... in dFaceNormdX, faceNum = %d, nodeNum(2) = %d\n", faceNum, nodeNum(2));

  dFaceNormdX.addContrib(faceNum, nodeNum(0), dFaceNormdX0[0]);
  dFaceNormdX.addContrib(faceNum, nodeNum(1), dFaceNormdX1[0]);
  dFaceNormdX.addContrib(faceNum, nodeNum(2), dFaceNormdX2[0]);

}

//------------------------------------------------------------------------------
// Computation of the OUTWARD face normal
void FaceTria::computeNormalEZGCL1(double oodt, SVec<double,3> &Xn, SVec<double,3> &Xnp1, 
				   Vec<Vec3D> &faceNorm, Vec<double> &faceNormVel)
{

  Vec3D x_n[3] = {Xn[ nodeNum(0) ], Xn[ nodeNum(1) ], Xn[ nodeNum(2) ]};
  Vec3D x_np1[3] = {Xnp1[ nodeNum(0) ], Xnp1[ nodeNum(1) ], Xnp1[ nodeNum(2) ]};

  faceNorm[normNum] = 0.5 * ((x_np1[2] - x_np1[0]) ^ (x_np1[1] - x_np1[0]));
  //EZ1 faceNorm = 0.5 * ((x_n[2] - x_n[0]) ^ (x_n[1] - x_n[0]));

  double vol = Face::computeVolume(x_n[0], x_n[2], x_n[1], x_np1[0], x_np1[2], x_np1[1]);

  faceNormVel[normNum] = oodt * vol;

}

//------------------------------------------------------------------------------
// computation of the OUTWARD face normal

// Included (MB)
void FaceTria::computeNormalAndDerivative(SVec<double,3> &X, SVec<double,3> &dX, Vec3D &faceNorm, Vec3D &dFaceNorm)
{

  Vec3D x[3] = {X[ nodeNum(0) ], X[ nodeNum(1) ], X[ nodeNum(2) ]};

  Vec3D dx[3] = {dX[ nodeNum(0) ], dX[ nodeNum(1) ], dX[ nodeNum(2) ]};

  faceNorm = 0.5 * ((x[2] - x[0]) ^ (x[1] - x[0]));

  dFaceNorm = 0.5 * (((dx[2] - dx[0]) ^ (x[1] - x[0])) + ((x[2] - x[0]) ^ (dx[1] - dx[0])));

}

//------------------------------------------------------------------------------
// Get OUTWARD face normal
Vec3D FaceTria::getNormal(Vec<Vec3D> &faceNorm) {

  return faceNorm[normNum];

}

//------------------------------------------------------------------------------
// Get OUTWARD face normal

// Included (MB)
Vec3D FaceTria::getdNormal(Vec<Vec3D> &facedNorm) {

  return facedNorm[normNum];

}

//------------------------------------------------------------------------------

double FaceTria::getNormalVel(Vec<double> &faceNormVel) {
  return faceNormVel[normNum];
}

//------------------------------------------------------------------------------

// Included (MB)
double FaceTria::getdNormalVel(Vec<double> &facedNormVel) {
  return facedNormVel[normNum];
}

//------------------------------------------------------------------------------
// Get i-th OUTWARD subface normal
Vec3D FaceTria::getNormal(Vec<Vec3D> &faceNorm, int i) {

  return third*faceNorm[normNum];

}


//------------------------------------------------------------------------------
// Get i-th OUTWARD subface normal
// Included (MB)
Vec3D FaceTria::getdNormal(Vec<Vec3D> &facedNorm, int i) {

  return third*facedNorm[normNum];

}

//------------------------------------------------------------------------------

double FaceTria::getNormalVel(Vec<double> &faceNormVel, int i) {
  return third*faceNormVel[normNum];
}

//------------------------------------------------------------------------------

// Included (MB)
double FaceTria::getdNormalVel(Vec<double> &facedNormVel, int i) {
  return third*facedNormVel[normNum];
}

//------------------------------------------------------------------------------
