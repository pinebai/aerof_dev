#ifndef _DENSE_MATRIX_OPS_H_
#define _DENSE_MATRIX_OPS_H_

#include <Vector3D.h>

#include <cstdio>
#include <cstdlib>

#include <cmath>

#include <complex>

typedef std::complex<double> bcomp;

//------------------------------------------------------------------------------
// warning: b and c must be different

template <class Scalar, int dim, int dim2>
class DenseMatrixOp {

public:

  static void applyToDenseMatrix(Scalar (*a)[dim2], int k, Scalar (*b)[dim2], int i, 
				 Scalar (*c)[dim2], int j) {
    fprintf(stderr, "*** Warning: applyToDenseMatrix routine is not optimized for %d\n", dim);

    for (int jj=0; jj<dim; ++jj) {
      for (int ii=0; ii<dim; ++ii) {
        c[j][dim*ii+jj] = 0.0;
        for (int kk=0; kk<dim; ++kk) 
          c[j][dim*ii+jj] += a[k][dim*ii+kk] * b[i][dim*kk+jj];
      }
    }
  }

  template<class Scalar2>
  static void applyToVector(Scalar (*a)[dim2], int k, Scalar2 (*b)[dim], int i, 
			    Scalar2 (*c)[dim], int j) {
    fprintf(stderr, "*** Warning: applyToVector routine is not optimized for %d\n", dim);

    for (int ii=0; ii<dim; ++ii) {
      c[j][ii] = 0.0;
      for (int kk=0; kk<dim; ++kk)
        c[j][ii] += a[k][dim*ii+kk] * b[i][kk];
    }
  }

  static void applyToVector(bcomp (*a)[dim2], int k, double (*b)[dim], int i,
                            double (*c)[dim], int j) {

    fprintf(stderr, " *** WARNING: Incompatible Types\n");
    exit(-1);
  }


  static void applyTransToVector(Scalar (*a)[dim2], int k, double (*b)[dim], int i,
                                 double (*c)[dim], int j) {
    fprintf(stderr, "*** Warning: applyTransToVector routine is not optimized for %d\n", dim);
    fprintf(stderr, " *** WARNING: Bill's non-optimized routine\n");
    exit(-1);
  }

  static void applyTransToVector(Scalar (*a)[dim2], int k, bcomp (*b)[dim], int i,
                                 bcomp (*c)[dim], int j) {
    fprintf(stderr, "*** Warning: applyTransToVector routine is not optimized for %d\n", dim);
    fprintf(stderr, " *** WARNING: Bill's non-optimized routine for bcomp\n");
    exit(-1);
  }

  template<class Scalar2>
  static void applyAndAddToVector(Scalar (*a)[dim2], int k, Scalar2 (*b)[dim], int i, 
				  Scalar2 (*c)[dim], int j) {
    fprintf(stderr, "*** Warning: applyAndAddToVector routine is not optimized for %d\n", dim);

    for (int ii=0; ii<dim; ++ii)
      for (int kk=0; kk<dim; ++kk)
        c[j][ii] += a[k][dim*ii+kk] * b[i][kk];
  }

  template<class Scalar2>
  static void applyAndSubToVector(Scalar (*a)[dim2], int k, Scalar2 (*b)[dim], int i, 
				  Scalar2 (*c)[dim], int j) {
    fprintf(stderr, "*** Warning: applyAndSubToVector routine is not optimized for %d\n", dim);

    for (int ii=0; ii<dim; ++ii)
      for (int kk=0; kk<dim; ++kk)
        c[j][ii] -= a[k][dim*ii+kk] * b[i][kk];
  }

  template<class Scalar2>
  static void applyTransAndAddToVector(Scalar (*a)[dim2], int k, Scalar2 (*b)[dim], int i,
                                       Scalar2 (*c)[dim], int j) {
    fprintf(stderr, "*** Warning: applyTransAndAddToVector routine is not optimized for %d\n", dim);

    fprintf(stderr, " *** WARNING: Bill's non optimized routine\n");
    exit(-1);
  }

  static void applyTransAndSubToVector(Scalar (*a)[dim2], int k, double (*b)[dim], int i,
                                       double (*c)[dim], int j) {
    fprintf(stderr, "*** Warning: applyTransAndSubToVector routine is not optimized for %d\n", dim);

    for (int ii=0; ii<dim; ++ii)
      for (int kk=0; kk<dim; ++kk)
        c[j][ii] -= a[k][dim*kk+ii] * b[i][kk];
  }

  static void applyTransAndSubToVector(Scalar (*a)[dim2], int k, bcomp (*b)[dim], int i,
                                       bcomp (*c)[dim], int j) {
    fprintf(stderr, "*** Warning: applyTransAndSubToVector routine is not optimized for %d\n",
dim);

    for (int ii=0; ii<dim; ++ii)
      for (int kk=0; kk<dim; ++kk)
        c[j][ii] -= a[k][dim*kk+ii] * b[i][kk];
  }

  static void transpose(Scalar *a, Scalar *b) {
    fprintf(stderr, "*** Warning: transpose routine is not optimized for %d\n", dim);

    for (int i=0; i<dim; ++i)
      for (int j=0; j<dim; ++j)
        b[dim*i + j] = a[i + dim*j];
  }


  static void lu(Scalar *a, Scalar *b, int n) {
                                                                                                                 
    int *index;
    index = new int[n];
    double d;
    d  = 1.0;
    ludec(a,index,d,n);
    ludfdbksb(a,index,b,n);
    delete [] index;
  }
  
  inline static double abs(double x) { return (x > 0.0 ? x : (-x)); }
  inline static double abs(bcomp x) { return std::abs(x); }
                                                        
  static void ludec(Scalar *a, int *index, double d, int n)  {
                                                        
    int i, imax, j, k;
    double big, temp;
    Scalar sum, dum;
    double *vv = new double[n];
    double epsilon = 1.0e-20;
                                                        
                                                        
    for (i=0; i<n; i++) {
      big = 0.0;
      for (j=0; j<n; j++){
        if ((temp = abs(*(a+i*n+j))) > big) big = temp;
      }
      if (big == 0.0) {
        fprintf(stderr,"Singular Matrix in ludec\n");
        return;
      }
      vv[i] = 1.0 / big;
    }
   for (j=0; j<n; j++) {
      for (i=0; i<j; i++) {
        sum = (*(a+n*i+j));
        for (k=0; k<i; k++)
          sum -= *(a+n*i+k)*(*(a+k*n+j));
        *(a+i*n+j) = sum;
      }                                                 
      big = 0.0;
      for (i=j; i<n; i++) {
        sum = (*(a+i*n+j));
        for (k=0; k<j;k++)
          sum -= *(a+i*n+k)*(*(a+k*n+j));
        *(a+i*n+j) = sum;
        dum = vv[i]*abs(sum);
        if ( abs(dum)  >= big) {
          big = abs(dum);
          imax = i;
        }
      }
      if ( j != imax) {
        for (k=0; k<n; k++) {
          dum = (*(a+imax*n+k));
          *(a+imax*n+k) =(*(a+j*n+k));
          *(a+j*n+k) = dum;
        }
        d = -d;
        vv[imax] = vv[j];
      }
      index[j] = imax;
      if (*(a+j*n+j) == 0.0)  *(a+j*n+j) = epsilon;
                                                        
                                                        
      if (j != n) {
        dum = 1.0 / (*(a+j*n+j));
        for (i=j+1; i<n; i++)  *(a+i*n+j) *= dum;
      }
    }
    delete[] vv;
  }
                                                        
                                                        
  static void ludfdbksb(Scalar *a, int *index, Scalar *b, int n) {
                                                        
    int i, ip, j;
    int ii = -1;
    Scalar sum;
    for (i=0; i<n; i++) {
      ip = index[i];
      sum = b[ip];
      b[ip] = b[i];
      if (ii>= 0)
        for(j=ii; j<=i-1; j++)   sum -= (*(a+i*n+j))*b[j];
      else if (abs(sum)!=0.0)  ii = i;
      b[i] =sum;
    }
    for (i=n-1; i >=0; i--) {
      sum = b[i];
      for (j=i+1; j<n; j++) sum -= (*(a+i*n+j))*b[j];
      b[i] = sum / (*(a+i*n+i));
    }
                                                        
  }
};

//------------------------------------------------------------------------------

template <class Scalar, int dim>
class VectorOp {

public:

  static void sum(Scalar (*a)[dim], int k, Scalar (*b)[dim], int i, 
		  Scalar (*c)[dim], int j) {

// Modified (MB)
//    fprintf(stderr, "*** Warning: sum routine is not optimized for %d\n", dim);

    for (int ii=0; ii<dim; ++ii)
      c[j][ii] = a[k][ii] + b[i][ii];
  }

  static void add(Scalar (*b)[dim], int i, Scalar (*c)[dim], int j) {

// Modified (MB)
//    fprintf(stderr, "*** Warning: add routine is not optimized for %d\n", dim);

    for (int ii=0; ii<dim; ++ii)
      c[j][ii] += b[i][ii];
  }

  static void sub(Scalar (*b)[dim], int i, Scalar (*c)[dim], int j) {

// Modified (MB)
//    fprintf(stderr, "*** Warning: routine sub is not optimized for %d\n", dim);

    for (int ii=0; ii<dim; ++ii)
      c[j][ii] -= b[i][ii];
  }

};

//------------------------------------------------------------------------------

#define DENSEMATRIXTIMESDENSEMATRIX1(op, a, b, c) \
c[0] op a[0]*b[0];

#define DENSEMATRIXTIMESDENSEMATRIX2(op, a, b, c) \
c[0] op a[0]*b[0]+a[1]*b[2]; \
c[1] op a[0]*b[1]+a[1]*b[3]; \
c[2] op a[2]*b[0]+a[3]*b[2]; \
c[3] op a[2]*b[1]+a[3]*b[3];

#define DENSEMATRIXTIMESDENSEMATRIX3(op, a, b, c) \
c[0] op a[0]*b[0]+a[1]*b[3]+a[2]*b[6]; \
c[1] op a[0]*b[1]+a[1]*b[4]+a[2]*b[7]; \
c[2] op a[0]*b[2]+a[1]*b[5]+a[2]*b[8]; \
c[3] op a[3]*b[0]+a[4]*b[3]+a[5]*b[6]; \
c[4] op a[3]*b[1]+a[4]*b[4]+a[5]*b[7]; \
c[5] op a[3]*b[2]+a[4]*b[5]+a[5]*b[8]; \
c[6] op a[6]*b[0]+a[7]*b[3]+a[8]*b[6]; \
c[7] op a[6]*b[1]+a[7]*b[4]+a[8]*b[7]; \
c[8] op a[6]*b[2]+a[7]*b[5]+a[8]*b[8];

#define DENSEMATRIXTIMESDENSEMATRIX4(op, a, b, c) \
c[0] op a[0]*b[0]+a[1]*b[4]+a[2]*b[8]+a[3]*b[12]; \
c[1] op a[0]*b[1]+a[1]*b[5]+a[2]*b[9]+a[3]*b[13]; \
c[2] op a[0]*b[2]+a[1]*b[6]+a[2]*b[10]+a[3]*b[14]; \
c[3] op a[0]*b[3]+a[1]*b[7]+a[2]*b[11]+a[3]*b[15]; \
c[4] op a[4]*b[0]+a[5]*b[4]+a[6]*b[8]+a[7]*b[12]; \
c[5] op a[4]*b[1]+a[5]*b[5]+a[6]*b[9]+a[7]*b[13]; \
c[6] op a[4]*b[2]+a[5]*b[6]+a[6]*b[10]+a[7]*b[14]; \
c[7] op a[4]*b[3]+a[5]*b[7]+a[6]*b[11]+a[7]*b[15]; \
c[8] op a[8]*b[0]+a[9]*b[4]+a[10]*b[8]+a[11]*b[12]; \
c[9] op a[8]*b[1]+a[9]*b[5]+a[10]*b[9]+a[11]*b[13]; \
c[10] op a[8]*b[2]+a[9]*b[6]+a[10]*b[10]+a[11]*b[14]; \
c[11] op a[8]*b[3]+a[9]*b[7]+a[10]*b[11]+a[11]*b[15]; \
c[12] op a[12]*b[0]+a[13]*b[4]+a[14]*b[8]+a[15]*b[12]; \
c[13] op a[12]*b[1]+a[13]*b[5]+a[14]*b[9]+a[15]*b[13]; \
c[14] op a[12]*b[2]+a[13]*b[6]+a[14]*b[10]+a[15]*b[14]; \
c[15] op a[12]*b[3]+a[13]*b[7]+a[14]*b[11]+a[15]*b[15];

#define DENSEMATRIXTIMESDENSEMATRIX5(op, a, b, c) \
c[0] op a[0]*b[0]+a[1]*b[5]+a[2]*b[10]+a[3]*b[15]+a[4]*b[20]; \
c[1] op a[0]*b[1]+a[1]*b[6]+a[2]*b[11]+a[3]*b[16]+a[4]*b[21]; \
c[2] op a[0]*b[2]+a[1]*b[7]+a[2]*b[12]+a[3]*b[17]+a[4]*b[22]; \
c[3] op a[0]*b[3]+a[1]*b[8]+a[2]*b[13]+a[3]*b[18]+a[4]*b[23]; \
c[4] op a[0]*b[4]+a[1]*b[9]+a[2]*b[14]+a[3]*b[19]+a[4]*b[24]; \
c[5] op a[5]*b[0]+a[6]*b[5]+a[7]*b[10]+a[8]*b[15]+a[9]*b[20]; \
c[6] op a[5]*b[1]+a[6]*b[6]+a[7]*b[11]+a[8]*b[16]+a[9]*b[21]; \
c[7] op a[5]*b[2]+a[6]*b[7]+a[7]*b[12]+a[8]*b[17]+a[9]*b[22]; \
c[8] op a[5]*b[3]+a[6]*b[8]+a[7]*b[13]+a[8]*b[18]+a[9]*b[23]; \
c[9] op a[5]*b[4]+a[6]*b[9]+a[7]*b[14]+a[8]*b[19]+a[9]*b[24]; \
c[10] op a[10]*b[0]+a[11]*b[5]+a[12]*b[10]+a[13]*b[15]+a[14]*b[20]; \
c[11] op a[10]*b[1]+a[11]*b[6]+a[12]*b[11]+a[13]*b[16]+a[14]*b[21]; \
c[12] op a[10]*b[2]+a[11]*b[7]+a[12]*b[12]+a[13]*b[17]+a[14]*b[22]; \
c[13] op a[10]*b[3]+a[11]*b[8]+a[12]*b[13]+a[13]*b[18]+a[14]*b[23]; \
c[14] op a[10]*b[4]+a[11]*b[9]+a[12]*b[14]+a[13]*b[19]+a[14]*b[24]; \
c[15] op a[15]*b[0]+a[16]*b[5]+a[17]*b[10]+a[18]*b[15]+a[19]*b[20]; \
c[16] op a[15]*b[1]+a[16]*b[6]+a[17]*b[11]+a[18]*b[16]+a[19]*b[21]; \
c[17] op a[15]*b[2]+a[16]*b[7]+a[17]*b[12]+a[18]*b[17]+a[19]*b[22]; \
c[18] op a[15]*b[3]+a[16]*b[8]+a[17]*b[13]+a[18]*b[18]+a[19]*b[23]; \
c[19] op a[15]*b[4]+a[16]*b[9]+a[17]*b[14]+a[18]*b[19]+a[19]*b[24]; \
c[20] op a[20]*b[0]+a[21]*b[5]+a[22]*b[10]+a[23]*b[15]+a[24]*b[20]; \
c[21] op a[20]*b[1]+a[21]*b[6]+a[22]*b[11]+a[23]*b[16]+a[24]*b[21]; \
c[22] op a[20]*b[2]+a[21]*b[7]+a[22]*b[12]+a[23]*b[17]+a[24]*b[22]; \
c[23] op a[20]*b[3]+a[21]*b[8]+a[22]*b[13]+a[23]*b[18]+a[24]*b[23]; \
c[24] op a[20]*b[4]+a[21]*b[9]+a[22]*b[14]+a[23]*b[19]+a[24]*b[24];

#define DENSEMATRIXTIMESDENSEMATRIX6(op, a, b, c) \
c[0] op a[0]*b[0]+a[1]*b[6]+a[2]*b[12]+a[3]*b[18]+a[4]*b[24]+a[5]*b[30]; \
c[1] op a[0]*b[1]+a[1]*b[7]+a[2]*b[13]+a[3]*b[19]+a[4]*b[25]+a[5]*b[31]; \
c[2] op a[0]*b[2]+a[1]*b[8]+a[2]*b[14]+a[3]*b[20]+a[4]*b[26]+a[5]*b[32]; \
c[3] op a[0]*b[3]+a[1]*b[9]+a[2]*b[15]+a[3]*b[21]+a[4]*b[27]+a[5]*b[33]; \
c[4] op a[0]*b[4]+a[1]*b[10]+a[2]*b[16]+a[3]*b[22]+a[4]*b[28]+a[5]*b[34]; \
c[5] op a[0]*b[5]+a[1]*b[11]+a[2]*b[17]+a[3]*b[23]+a[4]*b[29]+a[5]*b[35]; \
c[6] op a[6]*b[0]+a[7]*b[6]+a[8]*b[12]+a[9]*b[18]+a[10]*b[24]+a[11]*b[30]; \
c[7] op a[6]*b[1]+a[7]*b[7]+a[8]*b[13]+a[9]*b[19]+a[10]*b[25]+a[11]*b[31]; \
c[8] op a[6]*b[2]+a[7]*b[8]+a[8]*b[14]+a[9]*b[20]+a[10]*b[26]+a[11]*b[32]; \
c[9] op a[6]*b[3]+a[7]*b[9]+a[8]*b[15]+a[9]*b[21]+a[10]*b[27]+a[11]*b[33]; \
c[10] op a[6]*b[4]+a[7]*b[10]+a[8]*b[16]+a[9]*b[22]+a[10]*b[28]+a[11]*b[34]; \
c[11] op a[6]*b[5]+a[7]*b[11]+a[8]*b[17]+a[9]*b[23]+a[10]*b[29]+a[11]*b[35]; \
c[12] op a[12]*b[0]+a[13]*b[6]+a[14]*b[12]+a[15]*b[18]+a[16]*b[24]+a[17]*b[30]; \
c[13] op a[12]*b[1]+a[13]*b[7]+a[14]*b[13]+a[15]*b[19]+a[16]*b[25]+a[17]*b[31]; \
c[14] op a[12]*b[2]+a[13]*b[8]+a[14]*b[14]+a[15]*b[20]+a[16]*b[26]+a[17]*b[32]; \
c[15] op a[12]*b[3]+a[13]*b[9]+a[14]*b[15]+a[15]*b[21]+a[16]*b[27]+a[17]*b[33]; \
c[16] op a[12]*b[4]+a[13]*b[10]+a[14]*b[16]+a[15]*b[22]+a[16]*b[28]+a[17]*b[34]; \
c[17] op a[12]*b[5]+a[13]*b[11]+a[14]*b[17]+a[15]*b[23]+a[16]*b[29]+a[17]*b[35]; \
c[18] op a[18]*b[0]+a[19]*b[6]+a[20]*b[12]+a[21]*b[18]+a[22]*b[24]+a[23]*b[30]; \
c[19] op a[18]*b[1]+a[19]*b[7]+a[20]*b[13]+a[21]*b[19]+a[22]*b[25]+a[23]*b[31]; \
c[20] op a[18]*b[2]+a[19]*b[8]+a[20]*b[14]+a[21]*b[20]+a[22]*b[26]+a[23]*b[32]; \
c[21] op a[18]*b[3]+a[19]*b[9]+a[20]*b[15]+a[21]*b[21]+a[22]*b[27]+a[23]*b[33]; \
c[22] op a[18]*b[4]+a[19]*b[10]+a[20]*b[16]+a[21]*b[22]+a[22]*b[28]+a[23]*b[34]; \
c[23] op a[18]*b[5]+a[19]*b[11]+a[20]*b[17]+a[21]*b[23]+a[22]*b[29]+a[23]*b[35]; \
c[24] op a[24]*b[0]+a[25]*b[6]+a[26]*b[12]+a[27]*b[18]+a[28]*b[24]+a[29]*b[30]; \
c[25] op a[24]*b[1]+a[25]*b[7]+a[26]*b[13]+a[27]*b[19]+a[28]*b[25]+a[29]*b[31]; \
c[26] op a[24]*b[2]+a[25]*b[8]+a[26]*b[14]+a[27]*b[20]+a[28]*b[26]+a[29]*b[32]; \
c[27] op a[24]*b[3]+a[25]*b[9]+a[26]*b[15]+a[27]*b[21]+a[28]*b[27]+a[29]*b[33]; \
c[28] op a[24]*b[4]+a[25]*b[10]+a[26]*b[16]+a[27]*b[22]+a[28]*b[28]+a[29]*b[34]; \
c[29] op a[24]*b[5]+a[25]*b[11]+a[26]*b[17]+a[27]*b[23]+a[28]*b[29]+a[29]*b[35]; \
c[30] op a[30]*b[0]+a[31]*b[6]+a[32]*b[12]+a[33]*b[18]+a[34]*b[24]+a[35]*b[30]; \
c[31] op a[30]*b[1]+a[31]*b[7]+a[32]*b[13]+a[33]*b[19]+a[34]*b[25]+a[35]*b[31]; \
c[32] op a[30]*b[2]+a[31]*b[8]+a[32]*b[14]+a[33]*b[20]+a[34]*b[26]+a[35]*b[32]; \
c[33] op a[30]*b[3]+a[31]*b[9]+a[32]*b[15]+a[33]*b[21]+a[34]*b[27]+a[35]*b[33]; \
c[34] op a[30]*b[4]+a[31]*b[10]+a[32]*b[16]+a[33]*b[22]+a[34]*b[28]+a[35]*b[34]; \
c[35] op a[30]*b[5]+a[31]*b[11]+a[32]*b[17]+a[33]*b[23]+a[34]*b[29]+a[35]*b[35];

#define DENSEMATRIXTIMESDENSEMATRIX7(op, a, b, c) \
c[0] op a[0]*b[0]+a[1]*b[7]+a[2]*b[14]+a[3]*b[21]+a[4]*b[28]+a[5]*b[35]+a[6]*b[42]; \
c[1] op a[0]*b[1]+a[1]*b[8]+a[2]*b[15]+a[3]*b[22]+a[4]*b[29]+a[5]*b[36]+a[6]*b[43]; \
c[2] op a[0]*b[2]+a[1]*b[9]+a[2]*b[16]+a[3]*b[23]+a[4]*b[30]+a[5]*b[37]+a[6]*b[44]; \
c[3] op a[0]*b[3]+a[1]*b[10]+a[2]*b[17]+a[3]*b[24]+a[4]*b[31]+a[5]*b[38]+a[6]*b[45]; \
c[4] op a[0]*b[4]+a[1]*b[11]+a[2]*b[18]+a[3]*b[25]+a[4]*b[32]+a[5]*b[39]+a[6]*b[46]; \
c[5] op a[0]*b[5]+a[1]*b[12]+a[2]*b[19]+a[3]*b[26]+a[4]*b[33]+a[5]*b[40]+a[6]*b[47]; \
c[6] op a[0]*b[6]+a[1]*b[13]+a[2]*b[20]+a[3]*b[27]+a[4]*b[34]+a[5]*b[41]+a[6]*b[48]; \
c[7] op a[7]*b[0]+a[8]*b[7]+a[9]*b[14]+a[10]*b[21]+a[11]*b[28]+a[12]*b[35]+a[13]*b[42]; \
c[8] op a[7]*b[1]+a[8]*b[8]+a[9]*b[15]+a[10]*b[22]+a[11]*b[29]+a[12]*b[36]+a[13]*b[43]; \
c[9] op a[7]*b[2]+a[8]*b[9]+a[9]*b[16]+a[10]*b[23]+a[11]*b[30]+a[12]*b[37]+a[13]*b[44]; \
c[10] op a[7]*b[3]+a[8]*b[10]+a[9]*b[17]+a[10]*b[24]+a[11]*b[31]+a[12]*b[38]+a[13]*b[45]; \
c[11] op a[7]*b[4]+a[8]*b[11]+a[9]*b[18]+a[10]*b[25]+a[11]*b[32]+a[12]*b[39]+a[13]*b[46]; \
c[12] op a[7]*b[5]+a[8]*b[12]+a[9]*b[19]+a[10]*b[26]+a[11]*b[33]+a[12]*b[40]+a[13]*b[47]; \
c[13] op a[7]*b[6]+a[8]*b[13]+a[9]*b[20]+a[10]*b[27]+a[11]*b[34]+a[12]*b[41]+a[13]*b[48]; \
c[14] op a[14]*b[0]+a[15]*b[7]+a[16]*b[14]+a[17]*b[21]+a[18]*b[28]+a[19]*b[35]+a[20]*b[42]; \
c[15] op a[14]*b[1]+a[15]*b[8]+a[16]*b[15]+a[17]*b[22]+a[18]*b[29]+a[19]*b[36]+a[20]*b[43]; \
c[16] op a[14]*b[2]+a[15]*b[9]+a[16]*b[16]+a[17]*b[23]+a[18]*b[30]+a[19]*b[37]+a[20]*b[44]; \
c[17] op a[14]*b[3]+a[15]*b[10]+a[16]*b[17]+a[17]*b[24]+a[18]*b[31]+a[19]*b[38]+a[20]*b[45]; \
c[18] op a[14]*b[4]+a[15]*b[11]+a[16]*b[18]+a[17]*b[25]+a[18]*b[32]+a[19]*b[39]+a[20]*b[46]; \
c[19] op a[14]*b[5]+a[15]*b[12]+a[16]*b[19]+a[17]*b[26]+a[18]*b[33]+a[19]*b[40]+a[20]*b[47]; \
c[20] op a[14]*b[6]+a[15]*b[13]+a[16]*b[20]+a[17]*b[27]+a[18]*b[34]+a[19]*b[41]+a[20]*b[48]; \
c[21] op a[21]*b[0]+a[22]*b[7]+a[23]*b[14]+a[24]*b[21]+a[25]*b[28]+a[26]*b[35]+a[27]*b[42]; \
c[22] op a[21]*b[1]+a[22]*b[8]+a[23]*b[15]+a[24]*b[22]+a[25]*b[29]+a[26]*b[36]+a[27]*b[43]; \
c[23] op a[21]*b[2]+a[22]*b[9]+a[23]*b[16]+a[24]*b[23]+a[25]*b[30]+a[26]*b[37]+a[27]*b[44]; \
c[24] op a[21]*b[3]+a[22]*b[10]+a[23]*b[17]+a[24]*b[24]+a[25]*b[31]+a[26]*b[38]+a[27]*b[45]; \
c[25] op a[21]*b[4]+a[22]*b[11]+a[23]*b[18]+a[24]*b[25]+a[25]*b[32]+a[26]*b[39]+a[27]*b[46]; \
c[26] op a[21]*b[5]+a[22]*b[12]+a[23]*b[19]+a[24]*b[26]+a[25]*b[33]+a[26]*b[40]+a[27]*b[47]; \
c[27] op a[21]*b[6]+a[22]*b[13]+a[23]*b[20]+a[24]*b[27]+a[25]*b[34]+a[26]*b[41]+a[27]*b[48]; \
c[28] op a[28]*b[0]+a[29]*b[7]+a[30]*b[14]+a[31]*b[21]+a[32]*b[28]+a[33]*b[35]+a[34]*b[42]; \
c[29] op a[28]*b[1]+a[29]*b[8]+a[30]*b[15]+a[31]*b[22]+a[32]*b[29]+a[33]*b[36]+a[34]*b[43]; \
c[30] op a[28]*b[2]+a[29]*b[9]+a[30]*b[16]+a[31]*b[23]+a[32]*b[30]+a[33]*b[37]+a[34]*b[44]; \
c[31] op a[28]*b[3]+a[29]*b[10]+a[30]*b[17]+a[31]*b[24]+a[32]*b[31]+a[33]*b[38]+a[34]*b[45]; \
c[32] op a[28]*b[4]+a[29]*b[11]+a[30]*b[18]+a[31]*b[25]+a[32]*b[32]+a[33]*b[39]+a[34]*b[46]; \
c[33] op a[28]*b[5]+a[29]*b[12]+a[30]*b[19]+a[31]*b[26]+a[32]*b[33]+a[33]*b[40]+a[34]*b[47]; \
c[34] op a[28]*b[6]+a[29]*b[13]+a[30]*b[20]+a[31]*b[27]+a[32]*b[34]+a[33]*b[41]+a[34]*b[48]; \
c[35] op a[35]*b[0]+a[36]*b[7]+a[37]*b[14]+a[38]*b[21]+a[39]*b[28]+a[40]*b[35]+a[41]*b[42]; \
c[36] op a[35]*b[1]+a[36]*b[8]+a[37]*b[15]+a[38]*b[22]+a[39]*b[29]+a[40]*b[36]+a[41]*b[43]; \
c[37] op a[35]*b[2]+a[36]*b[9]+a[37]*b[16]+a[38]*b[23]+a[39]*b[30]+a[40]*b[37]+a[41]*b[44]; \
c[38] op a[35]*b[3]+a[36]*b[10]+a[37]*b[17]+a[38]*b[24]+a[39]*b[31]+a[40]*b[38]+a[41]*b[45]; \
c[39] op a[35]*b[4]+a[36]*b[11]+a[37]*b[18]+a[38]*b[25]+a[39]*b[32]+a[40]*b[39]+a[41]*b[46]; \
c[40] op a[35]*b[5]+a[36]*b[12]+a[37]*b[19]+a[38]*b[26]+a[39]*b[33]+a[40]*b[40]+a[41]*b[47]; \
c[41] op a[35]*b[6]+a[36]*b[13]+a[37]*b[20]+a[38]*b[27]+a[39]*b[34]+a[40]*b[41]+a[41]*b[48]; \
c[42] op a[42]*b[0]+a[43]*b[7]+a[44]*b[14]+a[45]*b[21]+a[46]*b[28]+a[47]*b[35]+a[48]*b[42]; \
c[43] op a[42]*b[1]+a[43]*b[8]+a[44]*b[15]+a[45]*b[22]+a[46]*b[29]+a[47]*b[36]+a[48]*b[43]; \
c[44] op a[42]*b[2]+a[43]*b[9]+a[44]*b[16]+a[45]*b[23]+a[46]*b[30]+a[47]*b[37]+a[48]*b[44]; \
c[45] op a[42]*b[3]+a[43]*b[10]+a[44]*b[17]+a[45]*b[24]+a[46]*b[31]+a[47]*b[38]+a[48]*b[45]; \
c[46] op a[42]*b[4]+a[43]*b[11]+a[44]*b[18]+a[45]*b[25]+a[46]*b[32]+a[47]*b[39]+a[48]*b[46]; \
c[47] op a[42]*b[5]+a[43]*b[12]+a[44]*b[19]+a[45]*b[26]+a[46]*b[33]+a[47]*b[40]+a[48]*b[47]; \
c[48] op a[42]*b[6]+a[43]*b[13]+a[44]*b[20]+a[45]*b[27]+a[46]*b[34]+a[47]*b[41]+a[48]*b[48];

//------------------------------------------------------------------------------

#define DENSEMATRIXTIMESVECTOR1(op, a, b, c) \
c[0] op a[0]*b[0];

#define DENSEMATRIXTIMESVECTORLS(op, a, b, c) \
c op a[0]*b;

#define DENSEMATRIXTIMESVECTOR2(op, a, b, c) \
c[0] op a[0]*b[0]+a[1]*b[1]; \
c[1] op a[2]*b[0]+a[3]*b[1];

#define DENSEMATRIXTIMESVECTOR3(op, a, b, c) \
c[0] op a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; \
c[1] op a[3]*b[0]+a[4]*b[1]+a[5]*b[2]; \
c[2] op a[6]*b[0]+a[7]*b[1]+a[8]*b[2];

#define DENSEMATRIXTIMESVECTOR4(op, a, b, c) \
c[0] op a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3]; \
c[1] op a[4]*b[0]+a[5]*b[1]+a[6]*b[2]+a[7]*b[3]; \
c[2] op a[8]*b[0]+a[9]*b[1]+a[10]*b[2]+a[11]*b[3]; \
c[3] op a[12]*b[0]+a[13]*b[1]+a[14]*b[2]+a[15]*b[3];

#define DENSEMATRIXTIMESVECTOR5(op, a, b, c) \
c[0] op a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3]+a[4]*b[4]; \
c[1] op a[5]*b[0]+a[6]*b[1]+a[7]*b[2]+a[8]*b[3]+a[9]*b[4]; \
c[2] op a[10]*b[0]+a[11]*b[1]+a[12]*b[2]+a[13]*b[3]+a[14]*b[4]; \
c[3] op a[15]*b[0]+a[16]*b[1]+a[17]*b[2]+a[18]*b[3]+a[19]*b[4]; \
c[4] op a[20]*b[0]+a[21]*b[1]+a[22]*b[2]+a[23]*b[3]+a[24]*b[4];

#define DENSEMATRIXTRANSTIMESVECTOR5(op, a, b, c) \
c[0] op a[0]*b[0] + a[5]*b[1] + a[10]*b[2] + a[15]*b[3] + a[20]*b[4]; \
c[1] op a[1]*b[0] + a[6]*b[1] + a[11]*b[2] + a[16]*b[3] + a[21]*b[4]; \
c[2] op a[2]*b[0] + a[7]*b[1] + a[12]*b[2] + a[17]*b[3] + a[22]*b[4]; \
c[3] op a[3]*b[0] + a[8]*b[1] + a[13]*b[2] + a[18]*b[3] + a[23]*b[4]; \
c[4] op a[4]*b[0] + a[9]*b[1] + a[14]*b[2] + a[19]*b[3] + a[24]*b[4];

#define DENSEMATRIXTRANSTIMESVECTOR6(op, a, b, c) \
c[0] op a[0]*b[0] + a[6]*b[1] + a[12]*b[2] + a[18]*b[3] + a[24]*b[4] + a[30]*b[5]; \
c[1] op a[1]*b[0] + a[7]*b[1] + a[13]*b[2] + a[19]*b[3] + a[25]*b[4] + a[31]*b[5]; \
c[2] op a[2]*b[0] + a[8]*b[1] + a[14]*b[2] + a[20]*b[3] + a[26]*b[4] + a[32]*b[5]; \
c[3] op a[3]*b[0] + a[9]*b[1] + a[15]*b[2] + a[21]*b[3] + a[27]*b[4] + a[33]*b[5]; \
c[4] op a[4]*b[0] + a[10]*b[1] + a[16]*b[2] + a[22]*b[3] + a[28]*b[4] + a[34]*b[5]; \
c[5] op a[5]*b[0] + a[11]*b[1] + a[17]*b[2] + a[23]*b[3] + a[29]*b[4] + a[35]*b[5];

#define DENSEMATRIXTIMESVECTOR6(op, a, b, c) \
c[0] op a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3]+a[4]*b[4]+a[5]*b[5]; \
c[1] op a[6]*b[0]+a[7]*b[1]+a[8]*b[2]+a[9]*b[3]+a[10]*b[4]+a[11]*b[5]; \
c[2] op a[12]*b[0]+a[13]*b[1]+a[14]*b[2]+a[15]*b[3]+a[16]*b[4]+a[17]*b[5]; \
c[3] op a[18]*b[0]+a[19]*b[1]+a[20]*b[2]+a[21]*b[3]+a[22]*b[4]+a[23]*b[5]; \
c[4] op a[24]*b[0]+a[25]*b[1]+a[26]*b[2]+a[27]*b[3]+a[28]*b[4]+a[29]*b[5]; \
c[5] op a[30]*b[0]+a[31]*b[1]+a[32]*b[2]+a[33]*b[3]+a[34]*b[4]+a[35]*b[5];

#define DENSEMATRIXTIMESVECTOR7(op, a, b, c) \
c[0] op a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3]+a[4]*b[4]+a[5]*b[5]+a[6]*b[6]; \
c[1] op a[7]*b[0]+a[8]*b[1]+a[9]*b[2]+a[10]*b[3]+a[11]*b[4]+a[12]*b[5]+a[13]*b[6]; \
c[2] op a[14]*b[0]+a[15]*b[1]+a[16]*b[2]+a[17]*b[3]+a[18]*b[4]+a[19]*b[5]+a[20]*b[6]; \
c[3] op a[21]*b[0]+a[22]*b[1]+a[23]*b[2]+a[24]*b[3]+a[25]*b[4]+a[26]*b[5]+a[27]*b[6]; \
c[4] op a[28]*b[0]+a[29]*b[1]+a[30]*b[2]+a[31]*b[3]+a[32]*b[4]+a[33]*b[5]+a[34]*b[6]; \
c[5] op a[35]*b[0]+a[36]*b[1]+a[37]*b[2]+a[38]*b[3]+a[39]*b[4]+a[40]*b[5]+a[41]*b[6]; \
c[6] op a[42]*b[0]+a[43]*b[1]+a[44]*b[2]+a[45]*b[3]+a[46]*b[4]+a[47]*b[5]+a[48]*b[6];

//------------------------------------------------------------------------------

template <class Scalar>
class DenseMatrixOp<Scalar,1,1> {

public:

  static void applyToDenseMatrix(Scalar (*a)[1], int k, Scalar (*b)[1], int i, 
				 Scalar (*c)[1], int j) {
    DENSEMATRIXTIMESDENSEMATRIX1(=, a[k], b[i], c[j]);
  }

  template<class Scalar2>
  static void applyToVector(Scalar (*a)[1], int k, Scalar2 (*b)[1], int i, 
			    Scalar2 (*c)[1], int j) {
    DENSEMATRIXTIMESVECTOR1(=, a[k], b[i], c[j]);
  }

  template<class Scalar2>
  static void applyToVector(Scalar (*a)[1], int k, Scalar2 (*b), int i,
                            Scalar2 (*c), int j) {
    DENSEMATRIXTIMESVECTORLS(=, a[k], b[i], c[j]);
  }

  static void applyTransToVector(Scalar (*a)[1], int k, double (*b)[1], int i,
                                 double (*c)[1], int j) {
    fprintf(stderr, "*** Warning: applyTransToVector routine is not optimized for dim of 1\n");
    fprintf(stderr, " *** WARNING: Bill's non-optimized routine\n");
    exit(-1);
  }

  template<class Scalar2>
  static void applyAndAddToVector(Scalar (*a)[1], int k, Scalar2 (*b)[1], int i, 
				  Scalar2 (*c)[1], int j) {
    DENSEMATRIXTIMESVECTOR1(+=, a[k], b[i], c[j]);
  }

  template<class Scalar2>
  static void applyAndAddToVector(Scalar (*a)[1], int k, Scalar2 (*b), int i,
                                  Scalar2 (*c), int j) {
    DENSEMATRIXTIMESVECTORLS(+=, a[k], b[i], c[j]);
  }

  template<class Scalar2>
  static void applyAndSubToVector(Scalar (*a)[1], int k, Scalar2 (*b)[1], int i, 
				  Scalar2 (*c)[1], int j) {
    DENSEMATRIXTIMESVECTOR1(-=, a[k], b[i], c[j]); 
  }

  template<class Scalar2>
  static void applyTransAndAddToVector(Scalar (*a)[1], int k, Scalar2 (*b)[1], int i,
                                       Scalar2 (*c)[1], int j) {
    fprintf(stderr, "*** Warning: applyTransAndAddToVector routine not implemented for dim of 1 \n");
    exit(-1);
  }

  static void applyTransAndSubToVector(Scalar (*a)[1], int k, double (*b)[1], int i,
                                       double (*c)[1], int j) {
    fprintf(stderr, "*** Not Implemented applyTransAndSubToVector for dim of 1\n");
  }

  static void transpose(Scalar *a, Scalar *b) {
    b[0] = a[0]; 
  }

};

//------------------------------------------------------------------------------

template <class Scalar>
class DenseMatrixOp<Scalar,2,4> {

public:

  static void applyToDenseMatrix(Scalar (*a)[4], int k, Scalar (*b)[4], int i, 
				 Scalar (*c)[4], int j) {
    DENSEMATRIXTIMESDENSEMATRIX2(=, a[k], b[i], c[j]);
  }

  template<class Scalar2>
  static void applyToVector(Scalar (*a)[4], int k, Scalar2 (*b)[2], int i, 
			    Scalar2 (*c)[2], int j) {
    DENSEMATRIXTIMESVECTOR2(=, a[k], b[i], c[j]);
  }

  static void applyTransToVector(Scalar (*a)[4], int k, double (*b)[2], int i,
                                 double (*c)[2], int j) {
    fprintf(stderr, "*** Warning: applyTransToVector routine is not implemented for dim of 2\n");
    fprintf(stderr, " *** WARNING: Bill's non-optimized routine\n");
    exit(-1);
  }


  template<class Scalar2>
  static void applyAndAddToVector(Scalar (*a)[4], int k, Scalar2 (*b)[2], int i, 
				  Scalar2 (*c)[2], int j) {
    DENSEMATRIXTIMESVECTOR2(+=, a[k], b[i], c[j]);
  }

  template<class Scalar2>
  static void applyTransAndAddToVector(Scalar (*a)[4], int k, Scalar2 (*b)[2], int i,
                                       Scalar2 (*c)[2], int j) {
    fprintf(stderr, "*** Warning: applyTransAndAddToVector routine is not optimized for dim of 2\n");

    fprintf(stderr, " *** WARNING: Bill's non optimized routine\n");
    exit(-1);
  }


  template<class Scalar2>
  static void applyAndSubToVector(Scalar (*a)[4], int k, Scalar2 (*b)[2], int i, 
				  Scalar2 (*c)[2], int j) {
    DENSEMATRIXTIMESVECTOR2(-=, a[k], b[i], c[j]); 
  }

  static void applyTransAndSubToVector(Scalar (*a)[4], int k, double (*b)[2], int i,
                                       double (*c)[2], int j) {
    fprintf(stderr, "*** Not Implemented applyTransAndSubToVector for dim of 2\n");
  }

  static void transpose(Scalar *a, Scalar *b) {
    b[0] = a[0]; 
    b[1] = a[2];
    b[2] = a[1]; 
    b[3] = a[3];
  }

};

//------------------------------------------------------------------------------

template <class Scalar>
class DenseMatrixOp<Scalar,3,9> {

public:

  static void applyToDenseMatrix(Scalar (*a)[9], int k, Scalar (*b)[9], int i, 
				 Scalar (*c)[9], int j) {
    DENSEMATRIXTIMESDENSEMATRIX3(=, a[k], b[i], c[j]);    
  }

  template<class Scalar2>
  static void applyToVector(Scalar (*a)[9], int k, Scalar2 (*b)[3], int i, 
			    Scalar2 (*c)[3], int j) {
    DENSEMATRIXTIMESVECTOR3(=, a[k], b[i], c[j]);
  }

  static void applyTransToVector(Scalar (*a)[9], int k, double (*b)[3], int i,
                                 double (*c)[3], int j) {
    fprintf(stderr, " applyTransToVector not implemented for DenseMatrixOp<Scalar,3,9>\n");
    exit(-1);
  }

  template<class Scalar2>
  static void applyAndAddToVector(Scalar (*a)[9], int k, Scalar2 (*b)[3], int i, 
				  Scalar2 (*c)[3], int j) {
    DENSEMATRIXTIMESVECTOR3(+=, a[k], b[i], c[j]);
  }

  template<class Scalar2>
  static void applyTransAndAddToVector(Scalar (*a)[9], int k, Scalar2 (*b)[3], int i,
                                       Scalar2 (*c)[3], int j) {
    fprintf(stderr, "*** Warning: applyTransAndAddToVector routine is not implemented for dim of 3 \n");
    exit(-1);
  }

  template<class Scalar2>
  static void applyAndSubToVector(Scalar (*a)[9], int k, Scalar2 (*b)[3], int i, 
				  Scalar2 (*c)[3], int j) {
    DENSEMATRIXTIMESVECTOR3(-=, a[k], b[i], c[j]); 
  }

  template<class Scalar2>
  static void applyTransAndSubToVector(Scalar (*a)[9], int k, Scalar2 (*b)[3], int i,
                                       Scalar2 (*c)[3], int j) {
    DENSEMATRIXTIMESVECTOR3(-=, a[k], b[i], c[j]);
  }


  static void transpose(Scalar *a, Scalar *b) {
    b[0]=a[0]; b[1]=a[3]; b[2]=a[6];
    b[3]=a[1]; b[4]=a[4]; b[5]=a[7];
    b[6]=a[2]; b[7]=a[5]; b[8]=a[8];
  }

};

//------------------------------------------------------------------------------

template <class Scalar>
class DenseMatrixOp<Scalar,5,25> {

public:

  static void applyToDenseMatrix(Scalar (*a)[25], int k, Scalar (*b)[25], int i, 
				 Scalar (*c)[25], int j) {
    DENSEMATRIXTIMESDENSEMATRIX5(=, a[k], b[i], c[j]);    
  }

  template<class Scalar2>
  static void applyToVector(Scalar (*a)[25], int k, Scalar2 (*b)[5], int i, 
			    Scalar2 (*c)[5], int j) {
    DENSEMATRIXTIMESVECTOR5(=, a[k], b[i], c[j]);
  }

  static void applyToVector(bcomp (*a)[25], int k, double (*b)[5], int i,
                            double (*c)[5], int j) {
 
    fprintf(stderr, " ... Incompatible types in DenseMatrixOps<bcomp, 5,25>::applyToVector\n");
    exit(-1);
  }

  template<class Scalar2>
  static void applyTransToVector(Scalar (*a)[25], int k, Scalar2 (*b)[5], int i,
                            Scalar2 (*c)[5], int j) {
    DENSEMATRIXTRANSTIMESVECTOR5(=, a[k], b[i], c[j]);
  }
  static void applyTransToVector(bcomp (*a)[25], int k, double (*b)[5], int i,
                            double (*c)[5], int j) {
    fprintf(stderr, " ... Incompatible types in DenseMatrixOps<bcomp, 5,25>::applyTransToVector\n");
    exit(-1);
  }

  template<class Scalar2>
  static void applyAndAddToVector(Scalar (*a)[25], int k, Scalar2 (*b)[5], int i, 
				  Scalar2 (*c)[5], int j) {

    DENSEMATRIXTIMESVECTOR5(+=, a[k], b[i], c[j]);
  }

    static void applyAndAddToVector(bcomp (*a)[25], int k, double (*b)[5], int i,
                                  double (*c)[5], int j) {

    fprintf(stderr, " ... Incompatible types in DenseMatrixOps<bcomp, 5,25>::applyAndAddToVector\n");
    exit(-1);
  }


  template<class Scalar2>
  static void applyTransAndAddToVector(Scalar (*a)[25], int k, Scalar2 (*b)[5], int i,
                                       Scalar2 (*c)[5], int j) {

    DENSEMATRIXTRANSTIMESVECTOR5(+=, a[k], b[i], c[j]);
  }

  static void applyTransAndAddToVector(bcomp (*a)[25], int k, double (*b)[5], int i,
                                       double (*c)[5], int j) {

    fprintf(stderr, " ... Incompatible types in DenseMatrixOps<bcomp, 5,25>::applyTransAndAddToVector\n");
    exit(-1);

  }

  template<class Scalar2>
  static void applyAndSubToVector(Scalar (*a)[25], int k, Scalar2 (*b)[5], int i,
                                  Scalar2 (*c)[5], int j) {
    DENSEMATRIXTIMESVECTOR5(-=, a[k], b[i], c[j]);
  }

    static void applyAndSubToVector(bcomp (*a)[25], int k, double (*b)[5], int i,
                                  double (*c)[5], int j) {

    fprintf(stderr, " ... Incompatible types in DenseMatrixOps<bcomp, 5,25>::applyAndSubToVector\n");
    exit(-1);
  }

  template<class Scalar2>
  static void applyTransAndSubToVector(Scalar (*a)[25], int k, Scalar2 (*b)[5], int i,
                                       Scalar2 (*c)[5], int j) {
    DENSEMATRIXTRANSTIMESVECTOR5(-=, a[k], b[i], c[j]);
  }

  static void applyTransAndSubToVector(bcomp (*a)[25], int k, double (*b)[5], int i,
                                       double (*c)[5], int j) {

    fprintf(stderr, " ... Incompatible types in DenseMatrixOps<bcomp, 5,25>::applyTransAndSubToVector\n");
    exit(-1);
  }


  static void transpose(Scalar *a, Scalar *b) {
    b[0]=a[0]; b[1]=a[5]; b[2]=a[10]; b[3]=a[15]; b[4]=a[20];
    b[5]=a[1]; b[6]=a[6]; b[7]=a[11]; b[8]=a[16]; b[9]=a[21];
    b[10]=a[2]; b[11]=a[7]; b[12]=a[12]; b[13]=a[17]; b[14]=a[22];
    b[15]=a[3]; b[16]=a[8]; b[17]=a[13]; b[18]=a[18]; b[19]=a[23];
    b[20]=a[4]; b[21]=a[9]; b[22]=a[14]; b[23]=a[19]; b[24]=a[24];
  }

};


//------------------------------------------------------------------------------

template <class Scalar>
class DenseMatrixOp<Scalar,6,36> {

public:

  static void applyToDenseMatrix(Scalar (*a)[36], int k, Scalar (*b)[36], int i, 
				 Scalar (*c)[36], int j) {
    DENSEMATRIXTIMESDENSEMATRIX6(=, a[k], b[i], c[j]);    
  }

  template<class Scalar2>
  static void applyToVector(Scalar (*a)[36], int k, Scalar2 (*b)[6], int i, 
			    Scalar2 (*c)[6], int j) {
    DENSEMATRIXTIMESVECTOR6(=, a[k], b[i], c[j]);
  }

  static void applyToVector(bcomp (*a)[36], int k, double (*b)[6], int i,
			    double (*c)[6], int j) {
    fprintf(stderr, " ... Incompatible types in DenseMatrixOps<bcomp,6,36>::applyToVector\n");
	  exit(-1);
	} 


  template<class Scalar2>
  static void applyTransToVector(Scalar (*a)[36], int k, Scalar2 (*b)[6], int i,
                            Scalar2 (*c)[6], int j) {
    fprintf(stderr, " *** Warning: applyTransToVector routine is not optimized for dim of 6\n");
    exit(-1);
  }
  


  static void applyTransToVector(Scalar (*a)[36], int k, double (*b)[6], int i,
                                 double (*c)[6], int j) {
    fprintf(stderr, "*** Warning: applyTransToVector routine is not optimized for dim of 6\n");
    fprintf(stderr, " *** WARNING: Bill's non-optimized routine\n");
    exit(-1);
  }


  template<class Scalar2>
  static void applyAndAddToVector(Scalar (*a)[36], int k, Scalar2 (*b)[6], int i, 
				  Scalar2 (*c)[6], int j) {
    DENSEMATRIXTIMESVECTOR6(+=, a[k], b[i], c[j]);
  }

  static void applyAndAddToVector(bcomp (*a)[36], int k, double (*b)[6], int i,
			                                  double (*c)[6], int j) {
	  fprintf(stderr, " ... Incompatible types in DenseMatrixOps<bcomp, 6,36>::applyAndAddToVector\n");
	  exit(-1);
	} 

  template<class Scalar2>
  static void applyTransAndAddToVector(Scalar (*a)[36], int k, Scalar2 (*b)[6], int i,
                                       Scalar2 (*c)[6], int j) {
    fprintf(stderr, "*** Warning: applyTransAndAddToVector routine is not optimized for dim of 6\n");
    exit(-1);
  }

	
  static void applyTransAndAddToVector(bcomp (*a)[36], int k, double (*b)[6], int i,
			                                 double (*c)[6], int j) {
		fprintf(stderr, " ... Incompatible types in DenseMatrixOps<bcomp, 6,36>::applyTransAndAddToVector\n");
	  exit(-1);
	} 

  template<class Scalar2>
  static void applyAndSubToVector(Scalar (*a)[36], int k, Scalar2 (*b)[6], int i, 
				  Scalar2 (*c)[6], int j) {
    DENSEMATRIXTIMESVECTOR6(-=, a[k], b[i], c[j]); 
  }

	
  static void applyAndSubToVector(bcomp (*a)[36], int k, double (*b)[6], int i,                                   double (*c)[6], int j) {

		fprintf(stderr, " ... Incompatible types in DenseMatrixOps<bcomp, 6,36>::applyAndSubToVector\n");
	  exit(-1);
	}

	template<class Scalar2> 
  static void applyTransAndSubToVector(Scalar (*a)[36], int k, Scalar2 (*b)[6], int i,
                                  Scalar2 (*c)[6], int j) { 
    fprintf(stderr, "*** Not Implemented applyTransAndSubToVector for dim of 6\n");
  }

	
  static void applyTransAndSubToVector(bcomp (*a)[36], int k, double (*b)[6], int i,
			                                 double (*c)[6], int j) {

		fprintf(stderr, " ... Incompatible types in DenseMatrixOps<bcomp, 6,36>::applyTransAndSubToVector\n");
	  exit(-1);
	}


  static void transpose(Scalar *a, Scalar *b) {
    b[0]=a[0]; b[1]=a[6]; b[2]=a[12]; b[3]=a[18]; b[4]=a[24]; b[5]=a[30];
    b[6]=a[1]; b[7]=a[7]; b[8]=a[13]; b[9]=a[19]; b[10]=a[25]; b[11]=a[31];
    b[12]=a[2]; b[13]=a[8]; b[14]=a[14]; b[15]=a[20]; b[16]=a[26]; b[17]=a[32];
    b[18]=a[3]; b[19]=a[9]; b[20]=a[15]; b[21]=a[21]; b[22]=a[27]; b[23]=a[33];
    b[24]=a[4]; b[25]=a[10]; b[26]=a[16]; b[27]=a[22]; b[28]=a[28]; b[29]=a[34];
    b[30]=a[5]; b[31]=a[11]; b[32]=a[17]; b[33]=a[23]; b[34]=a[29]; b[35]=a[35];
  }

};

//------------------------------------------------------------------------------

template <class Scalar>
class DenseMatrixOp<Scalar,7,49> {

public:

  static void applyToDenseMatrix(Scalar (*a)[49], int k, Scalar (*b)[49], int i, 
				 Scalar (*c)[49], int j) {
    DENSEMATRIXTIMESDENSEMATRIX7(=, a[k], b[i], c[j]);    
  }

  template<class Scalar2>
  static void applyToVector(Scalar (*a)[49], int k, Scalar2 (*b)[7], int i, 
			    Scalar2 (*c)[7], int j) {
    DENSEMATRIXTIMESVECTOR7(=, a[k], b[i], c[j]);
  }

  template<class Scalar2>
  static void applyTransToVector(Scalar (*a)[49], int k, Scalar2 (*b)[7], int i,
                            Scalar2 (*c)[7], int j) {
    fprintf(stderr, " *** Warning: applyTransToVector routine is not optimized for dim of 7\n");
    exit(-1);
  }


  static void applyTransToVector(Scalar (*a)[49], int k, double (*b)[7], int i,
                                 double (*c)[7], int j) {
    fprintf(stderr, "*** Warning: applyTransToVector routine is not optimized for dim of 7\n");
    fprintf(stderr, " *** WARNING: Bill's non-optimized routine\n");
    exit(-1);
  }

  template<class Scalar2>
  static void applyAndAddToVector(Scalar (*a)[49], int k, Scalar2 (*b)[7], int i, 
				  Scalar2 (*c)[7], int j) {
    DENSEMATRIXTIMESVECTOR7(+=, a[k], b[i], c[j]);
  }

  template<class Scalar2>
  static void applyTransAndAddToVector(Scalar (*a)[49], int k, Scalar2 (*b)[7], int i,
                                       Scalar2 (*c)[7], int j) {
    fprintf(stderr, "*** Warning: applyTransAndAddToVector routine is not optimized for dim of 7\n");
    exit(-1);
  }

  template<class Scalar2>
  static void applyAndSubToVector(Scalar (*a)[49], int k, Scalar2 (*b)[7], int i, 
				  Scalar2 (*c)[7], int j) {
    DENSEMATRIXTIMESVECTOR7(-=, a[k], b[i], c[j]); 
  }

  static void applyTransAndSubToVector(Scalar (*a)[49], int k, double (*b)[7], int i,
                                       double (*c)[7], int j) {
    fprintf(stderr, "*** Not Implemented applyTransAndSubToVector for dim of 7\n");
  }

  static void transpose(Scalar *a, Scalar *b) {
    b[0]=a[0]; b[1]=a[7]; b[2]=a[14]; b[3]=a[21]; b[4]=a[28]; b[5]=a[35]; b[6]=a[42];
    b[7]=a[1]; b[8]=a[8]; b[9]=a[15]; b[10]=a[22]; b[11]=a[29]; b[12]=a[36]; b[13]=a[43];
    b[14]=a[2]; b[15]=a[9]; b[16]=a[16]; b[17]=a[23]; b[18]=a[30]; b[19]=a[37]; b[20]=a[44];
    b[21]=a[3]; b[22]=a[10]; b[23]=a[17]; b[24]=a[24]; b[25]=a[31]; b[26]=a[38]; b[27]=a[45];
    b[28]=a[4]; b[29]=a[11]; b[30]=a[18]; b[31]=a[25]; b[32]=a[32]; b[33]=a[39]; b[34]=a[46];
    b[35]=a[5]; b[36]=a[12]; b[37]=a[19]; b[38]=a[26]; b[39]=a[33]; b[40]=a[40]; b[41]=a[47];
    b[42]=a[6]; b[43]=a[13]; b[44]=a[20]; b[45]=a[27]; b[46]=a[34]; b[47]=a[41]; b[48]=a[48];
  }

};

//------------------------------------------------------------------------------

template <class Scalar>
class VectorOp<Scalar,5> {

public:

  static void sum(Scalar (*a)[5], int k, Scalar (*b)[5], int i, Scalar (*c)[5], int j) {
    c[j][0] = a[k][0] + b[i][0];
    c[j][1] = a[k][1] + b[i][1];
    c[j][2] = a[k][2] + b[i][2];
    c[j][3] = a[k][3] + b[i][3];
    c[j][4] = a[k][4] + b[i][4];
  }

  static void add(Scalar (*b)[5], int i, Scalar (*c)[5], int j) {
    c[j][0] += b[i][0];
    c[j][1] += b[i][1];
    c[j][2] += b[i][2];
    c[j][3] += b[i][3];
    c[j][4] += b[i][4];
  }

  static void sub(Scalar (*b)[5], int i, Scalar (*c)[5], int j) {
    c[j][0] -= b[i][0];
    c[j][1] -= b[i][1];
    c[j][2] -= b[i][2];
    c[j][3] -= b[i][3];
    c[j][4] -= b[i][4];
  }

};

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
void denseMatrixTimesDenseMatrixInPlace(Scalar m1[dim][dim], Scalar m2[dim][dim])
{

  int i, j;
  Scalar res[dim][dim];

  for (i = 0; i < dim; ++i)
    for (j = 0; j < dim; ++j) {
      res[i][j] = 0.0;
      for (int k = 0; k < dim; ++k)
        res[i][j] += m1[i][k] * m2[k][j];
    }

  for (i = 0; i < dim; ++i)
    for (j = 0; j < dim; ++j)
      m2[i][j] = res[i][j];

}

//------------------------------------------------------------------------------

template<class Scalar, int n>
int LINPACKdgefa(Scalar *a, int *ipvt)
{

  int     i__2, i__3, kp1, nm1, j, k, l,ll,kn,knp1,jn1;
  Scalar  t,*ax,*ay,*aa;
  Scalar  tmp,max;

  // gaussian elimination with partial pivoting

  // parameter adjustments
  --ipvt;
  a       -= n + 1;

  // function Body
  nm1 = n - 1;
  for (k = 1; k <= nm1; ++k) {
    kp1  = k + 1;
    kn   = k*n;
    knp1 = k*n + k;

    // find l = pivot index

    i__2 = n - k + 1;
    aa = &a[knp1];
    max = sqrt(sqNorm(aa[0]));
    l = 1;
    for ( ll=1; ll<i__2; ll++ ) {
      tmp = sqrt(sqNorm(aa[ll]));
      if (sqNorm(tmp) > sqNorm(max)) { max = tmp; l = ll+1;}
    }
    l += k - 1;
    ipvt[k] = l;

    if (a[l + kn] == 0.) return k;

    // interchange if necessary

    if (l != k) {
      t = a[l + kn];
      a[l + kn] = a[knp1];
      a[knp1] = t;
    }

    // compute multipliers

    t = -1. / a[knp1];
    i__2 = n - k;
    aa = &a[1 + knp1]; 
    for ( ll=0; ll<i__2; ll++ ) {
      aa[ll] *= t;
    }

    // row elimination with column indexing

    ax = aa;
    for (j = kp1; j <= n; ++j) {
      jn1 = j*n;
      t = a[l + jn1];
      if (l != k) {
	a[l + jn1] = a[k + jn1];
	a[k + jn1] = t;
      }

      i__3 = n - k;
      ay = &a[1+k+jn1];
      for ( ll=0; ll<i__3; ll++ ) {
	ay[ll] += t*ax[ll];
      }
    }
  }

  ipvt[n] = n;

  if (a[n + n * n] == 0.) return n;

  return 0;

} 

//------------------------------------------------------------------------------

template<class Scalar, int n>
int LINPACKdgedi(Scalar *a, int *ipvt, Scalar *work)
{

  int     i__2,kb, kp1, nm1,i, j, k, l, ll,kn,knp1,jn1;
  Scalar  *aa,*ax,*ay,tmp;
  Scalar  t;

  --work;
  --ipvt;
  a       -= n + 1;

  // compute inverse(u)

  for (k = 1; k <= n; ++k) {
    kn           = k*n;
    knp1         = kn + k;
    a[knp1]      = 1.0 / a[knp1];
    t            = -a[knp1];
    i__2         = k - 1;
    aa           = &a[1 + kn]; 
    for ( ll=0; ll<i__2; ll++ ) aa[ll] *= t;
    kp1 = k + 1;
    if (n < kp1) continue;
    ax = aa;
    for (j = kp1; j <= n; ++j) {
      jn1 = j*n;
      t = a[k + jn1];
      a[k + jn1] = 0.;
      ay = &a[1 + jn1];
      for ( ll=0; ll<k; ll++ ) {
	ay[ll] += t*ax[ll];
      }
    }
  }

  // form inverse(u)*inverse(l)

  nm1 = n - 1;
  if (nm1 < 1) return 0;

  for (kb = 1; kb <= nm1; ++kb) {
    k   = n - kb;
    kn  = k*n;
    kp1 = k + 1;
    aa  = a + kn;
    for (i = kp1; i <= n; ++i) {
      work[i] = aa[i];
      aa[i]   = 0.;
    }
    for (j = kp1; j <= n; ++j) {
      t = work[j];
      ax = &a[j * n + 1];
      ay = &a[kn + 1];
      for ( ll=0; ll<n; ll++ ) {
	ay[ll] += t*ax[ll];
      }
    }
    l = ipvt[k];
    if (l != k) {
      ax = &a[kn + 1]; 
      ay = &a[l * n + 1];
      for ( ll=0; ll<n; ll++ ) {
	tmp    = ax[ll];
	ax[ll] = ay[ll];
	ay[ll] = tmp;
      }
    }
  }

  return 0;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
int invertDenseMatrix(Scalar *a)
{

  int ipvt[dim];
  Scalar work[dim];
  Scalar tmp[dim*dim];

  DenseMatrixOp<Scalar,dim,dim*dim>::transpose(a, tmp);

  int ierr = LINPACKdgefa<Scalar,dim>(tmp, ipvt);

  if (ierr) return ierr;

  ierr = LINPACKdgedi<Scalar,dim>(tmp, ipvt, work);

  DenseMatrixOp<Scalar,dim,dim*dim>::transpose(tmp, a);

  return ierr;

}

//------------------------------------------------------------------------------
// warning: b and c must be different

template <class Scalar, int dim, int dim2, int dim3>
class RectangularDenseMatrixOp {

public:

  template<class Scalar2>
  static void applyAndAddToVector(Scalar (*a)[dim3], int k, Scalar2 (*b)[dim], int i, 
				  Scalar2 (*c)[dim2], int j) {
//    fprintf(stderr, "*** Warning: applyAndAddToVector routine is not optimized for %d\n", dim);

    for (int ii=0; ii<dim2; ++ii)
      for (int kk=0; kk<dim; ++kk)
        c[j][ii] += a[k][dim*ii+kk] * b[i][kk];
  }

  template<class Scalar2>
  static void applyTransposeAndAddToVector(Scalar (*a)[dim3], int k, Scalar2 (*b)[dim2], int i, 
				  Scalar2 (*c)[dim], int j) {
//    fprintf(stderr, "*** Warning: applyTransposeAndAddToVector routine is not optimized for %d\n", dim2);

    for (int ii=0; ii<dim2; ++ii)
      for (int kk=0; kk<dim; ++kk)
        c[j][kk] += a[k][dim*ii+kk] * b[i][ii];
  }

  template<class Scalar2>
  static void applyTransposeAndAddToVector(Scalar (*a)[dim3], int k, Scalar2 (*b)[dim2], int i, 
				  Scalar2 *c, int j) {
//    fprintf(stderr, "*** Warning: applyTransposeAndAddToVector routine is not optimized for %d\n", dim2);

    for (int ii=0; ii<dim2; ++ii)
      c[j] += a[k][dim*ii] * b[i][ii];
  }

  template<class Scalar2>
  static void applyAndAddToVector(Scalar (*a)[dim3], int k, Scalar2 (*b)[dim], int i, 
				  Vec3D *c, int j) {
//    fprintf(stderr, "*** Warning: applyAndAddToVector routine is not optimized for %d\n", dim);

    for (int ii=0; ii<dim2; ++ii)
      for (int kk=0; kk<dim; ++kk)
        c[j][ii] += a[k][dim*ii+kk] * b[i][kk];
  }

  template<class Scalar2>
  static void applyAndAddToVector(Scalar (*a)[dim3], int k, Vec3D *b, int i, Scalar2 (*c)[dim2], int j) { 
//    fprintf(stderr, "*** Warning: applyAndAddToVector routine is not optimized for %d\n", dim);

    for (int ii=0; ii<dim2; ++ii)
      for (int kk=0; kk<dim; ++kk)
        c[j][ii] += a[k][dim*ii+kk] * b[i][kk];
  }

  template<class Scalar2>
  static void applyAndAddToVector(Scalar (*a)[dim3], int k, Scalar2 *b, int i, Scalar2 (*c)[dim2], int j) { 
//    fprintf(stderr, "*** Warning: applyAndAddToVector routine is not optimized for %d\n", dim);

    for (int ii=0; ii<dim2; ++ii)
      c[j][ii] += a[k][dim*ii] * b[i];
  }

  template<class Scalar2>
  static void applyTransposeAndAddToVector(Scalar (*a)[dim3], int k, Vec3D *b, int i, 
				  Scalar2 (*c)[dim], int j) {
//    fprintf(stderr, "*** Warning: applyTransposeAndAddToVector routine is not optimized for %d\n", dim2);

    for (int ii=0; ii<dim2; ++ii)
      for (int kk=0; kk<dim; ++kk)
        c[j][kk] += a[k][dim*ii+kk] * b[i][ii];
//    fprintf(stderr, " ... in applyTransposeAndAddToVector, dim(5)=%d, dim2(3)=%d \n",dim, dim2);
//    exit(-1);
  }

  template<class Scalar2>
  static void applyTransposeAndAddToVector(Scalar (*a)[dim3], int k, Scalar2 (*b)[dim2], int i, 
				  Vec3D *c, int j) {
//    fprintf(stderr, "*** Warning: applyTransposeAndAddToVector routine is not optimized for %d\n", dim2);

    for (int ii=0; ii<dim2; ++ii)
      for (int kk=0; kk<dim; ++kk)
        c[j][kk] += a[k][dim*ii+kk] * b[i][ii];
//    fprintf(stderr, " ... in applyTransposeAndAddToVector, dim(5)=%d, dim2(3)=%d \n",dim, dim2);
//    exit(-1);
  }

  template<class Scalar2>
  static void applyAndAddToVector(Scalar (*a)[dim3], int k, Scalar2 (*b)[dim], int i, 
				  Scalar2 *c, int j) {
//    fprintf(stderr, "*** Warning: applyAndAddToVector routine is not optimized for %d\n", dim);

      for (int kk=0; kk<dim; ++kk)
        c[j] += a[k][kk] * b[i][kk];
  }

  template<class Scalar2>
  static void applyTransposeAndAddToVector(Scalar (*a)[dim3], int k, Scalar2 *b, int i, 
				  Scalar2 (*c)[dim], int j) {
//    fprintf(stderr, "*** Warning: applyTransposeAndAddToVector routine is not optimized for %d\n", dim2);

      for (int kk=0; kk<dim; ++kk)
        c[j][kk] += a[k][kk] * b[i];
  }

};

#endif
