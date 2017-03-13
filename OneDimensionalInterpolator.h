#pragma once

#include <DenseMatrixOps.h>

class OneDimensionalInterpolator {

  static int factorial(int i) { 
    return (i == 0 ? 1 : i*factorial(i-1)); 
  }

 public:

  template <int dim>
    static void Interpolate(SVec<double,dim> &V, double& Vi, int id,
			    SVec<double,1> &X, double r,int order) {

    int j,i;
    for (i = 0; i < X.size(); ++i) {
      if (X[i][0] < r && X[i+1][0] >= r)
	break;
    }

    j = i-order/2;
    if (j < 0)
      j = 0;
    else if (j+order > X.size())
      j = X.size()-order;
    
    double* a = new double[order*order];
    double* b = new double[order];
    for (int k = 0; k < order; ++k) {
      for (int l = 0; l < order; ++l) {
	if (l == 0)
	  a[k*order+l] = 1.0 / factorial(l);
	else
	  a[k*order+l] = pow( X[j+k][0]-X[i][0] , l) / factorial(l);
      }
      
      b[k] = V[j+k][id];
    }

    DenseMatrixOp<double,5,5>::lu(a, b, order);
    
    Vi = 0.0;
    for (int l = 0; l < order; ++l)
      Vi += (l>0?pow( r-X[i][0] , l):1.0) / factorial(l)*b[l];

    delete [] a;
    delete [] b;
  }
  
};
