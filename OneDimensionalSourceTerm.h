#ifndef _ONE_DIMENSIONAL_SOURCE_TERM_H_
#define _ONE_DIMENSIONAL_SOURCE_TERM_H_

#include <cmath>

#include <DenseMatrixOps.h>

#include <sstream>

const static int factorial_oned[] = {1,1,2,6,24,120,720,5040,40320};
const static int permutation_oned[] = {1,-1,1,-1,1,-1,1,-1};

class OneDimensionalSourceTerm {

  struct MyLU {
    
    double* a;
    int* index;
  };

  double* pow_ry_k,*pow_ri_k,*logr1r0;

 public:
 
  OneDimensionalSourceTerm() { logr1r0 = 0; pow_ry_k = pow_ri_k = 0; }

  ~OneDimensionalSourceTerm() { }

  inline static int factorial(int i) { 
    return factorial_oned[i];//(i == 0 ? 1 : i*factorial(i-1)); 
  }

  void initialize(double _alpha, int _order, SVec<double,1>& X) {
    
    alpha = _alpha;
    order = _order;
    ludata = new MyLU[X.size()];
    int j;

    pow_ry_k = new double[(X.size()+1)*order];
    pow_ri_k = new double[(X.size())*order];
    for (int i = 0; i < X.size(); ++i) {
      
      j = i-order/2;
      if (j < 0)
	j = 0;
      else if (j+order > X.size())
	j = X.size()-order;

      MyLU lu = {new double[order*order], new int[order]};
      for (int k = 0; k < order; ++k) {
	for (int l = 0; l < order; ++l) {
	  if (l == 0)
	    lu.a[k*order+l] = 1.0 / factorial(l);
	  else
	    lu.a[k*order+l] = pow( X[j+k][0]-X[i][0] , l) / factorial(l);
	}
      }
      DenseMatrixOp<double,5,5>::ludec(lu.a, lu.index,1.0, order);
      ludata[i] = lu;
  
    }
  }

  int binomialTerm(int n, int k) {
    return factorial(n)/(factorial(n-k)*factorial(k));
  }

  double computeIntegralTerm(double r1, double r0,double ri, int n,int i) {
    
    double I = 0.0;
    I = permutation_oned[n]*pow_ri_k[i*(order+1)+n]*logr1r0[i];//log(r1/r0);
    for (int k = 1; k <= n; ++k) {
      I += 1.0/k*binomialTerm(n,k)*(k==n?1.0 : permutation_oned[n]*pow_ri_k[i*(order+1)+n-k]*(pow_ry_k[(i+1)*(order+1)+k]-pow_ry_k[(i)*(order+1)+k]));
    }
    return I*alpha;
  }

  template <class FluxF,int dim>
    void compute(FluxF& f, SVec<double,dim>& V, SVec<double,dim>& F, SVec<double,1>& X, SVec<double,1>& Y, Vec<int>& fluidId) {
    
    double* local = new double[dim*V.size()];
    double* derivs = new double[order];
    
    if (!logr1r0) {
      logr1r0 = new double[X.size()];
      pow_ri_k = new double[X.size()*(order+1)];
      pow_ry_k = new double[(X.size()+1)*(order+1)];
      for (int i = 0; i < X.size(); ++i) {
        logr1r0[i] = (i > 0 ? log(Y[i+1][0]/Y[i][0]) : 0.0);
        for (int k = 0; k <= order; ++k) {   
          pow_ri_k[i*(order+1)+k] = pow(X[i][0],k);
        } 
      }
      for (int i = 0; i < X.size()+1; ++i) {
        for (int k = 0; k <= order; ++k) {   
          pow_ry_k[i*(order+1)+k] = pow(Y[i][0],k);
        } 
      }
    }

    int j;
    for (int i = 0; i < V.size(); ++i) {
      f.compute(V[i], local+i*dim,i);
    }

    for (int i = 0; i < V.size(); ++i) {
     
      double Yi = Y[i][0],Yip1 = Y[i+1][0],Xi = X[i][0];
      j = i-order/2;
      int fid = fluidId[i];
      if (j < 0)
	j = 0;
      else if (j+order > X.size())
	j = X.size()-order;

      MyLU& lu = ludata[i];
      double* Fiptr = F[i];
      for (int k = 0; k < dim; ++k) {
	for (int l = 0; l < order; ++l) {
          if (fluidId[j+l] == fid)
	    derivs[l] = local[(j+l)*dim+k];
          else
            derivs[l] = local[i*dim+k];
	}
 	
	DenseMatrixOp<double,5,5>::ludfdbksb(lu.a, lu.index, derivs, order);
	
	double term = 0.0;
	for (int l = 0; l < order; ++l) {
	  Fiptr[k] += computeIntegralTerm(Yip1,Yi,Xi,l,i)*derivs[l]/factorial(l);
	}
	//Fiptr[k] += term;
      }
    } 

    delete [] local;
    delete [] derivs;
  }
  
 private:

  double alpha; // 1 for cylindrical, 2 for spherical

  int order;

  MyLU* ludata;

};

#endif
