#include <cstdio>
#include <iostream>

#include "DenseMatrix.h"

// CBM: for Factor and ReSolve
// ---------------------------
#include <LinkF77.h>

extern "C" {

  // triangular factorization of a real general matrix using Gaussian elimination with complete pivoting
  void F77NAME(dgecp)(const int &, double *, const int &, const double &, const int &, int &, int *, int *, int &);

  // solves a real general linear system using the triangular factorization computed by dgecp
  void F77NAME(dgers)(const int &, const int &, double *, const int &, const int &, const int *,
                       const int *, double *, const int &, int &);

}
// ---------------------------
template<class Scalar>
GenFullM<Scalar>::GenFullM() 
{
 nrow    = 0;
 ncolumn = 0;
 //v = new double[1];
 v = 0;

 iprow = NULL;
 ipcol = NULL;

}

template<class Scalar>
GenFullM<Scalar>::GenFullM(int nr) 
{
 nrow    = nr;
 ncolumn = nr;
 v = new Scalar[nrow*ncolumn];

 iprow = NULL;
 ipcol = NULL;

}

template<class Scalar>
GenFullM<Scalar>::GenFullM(int nr, int nc) 
{
 nrow    = nr;
 ncolumn = nc;
 if(nrow*ncolumn == 0)
   v = new Scalar[1];
 else
   v = new Scalar[nrow*ncolumn];

 iprow = NULL;
 ipcol = NULL;

}

template<class Scalar>
GenFullM<Scalar>::GenFullM(const GenFullM &m) 
{
 nrow    = m.nrow;
 ncolumn = m.ncolumn;
 v = new Scalar[nrow*ncolumn];
 int i;
 for(i=0; i < nrow*ncolumn; ++i)
   v[i] = m.v[i];

 iprow = NULL;
 ipcol = NULL;

}

template<class Scalar>
GenFullM<Scalar>::GenFullM(const GenFullM &m, int nr, int sr, int nc, int sc)
{
 nrow    = nr;
 ncolumn = nc;
 v = new Scalar[nrow*ncolumn] ;

 int i,j;
 for(i=0; i < nrow; ++i)
  for(j=0; j < ncolumn; ++j)
    (*this)[i][j] = m[i+sr][j+sc] ;

 iprow = NULL;
 ipcol = NULL;

}

template<class Scalar>
GenFullM<Scalar>::~GenFullM()
{
 if(v) delete [] v;
 if(iprow) delete [] iprow;
 if(ipcol) delete [] ipcol;
}

template<class Scalar>
void
GenFullM<Scalar>::setNewSize(int nr, int nc, double initVal)
{
 nrow    = nr;
 ncolumn = nc;
 delete [] v;
 v = new Scalar[nrow*ncolumn];
 int i;
 for(i=0; i < nrow*ncolumn; ++i)
   v[i] = initVal;
}

template<class Scalar>
void
GenFullM<Scalar>::setNewSize(int nr, double initVal)
{
 nrow = nr;
 ncolumn = nr;
 delete [] v;
 v = new Scalar[nrow*ncolumn];
 int i;
 for(i=0; i<nrow*ncolumn; ++i)
   v[i] = initVal;
}

template<class Scalar>
void
GenFullM<Scalar>::operator=(const GenFullM<Scalar> &m)
{
 if(m.v == v) return ;

 if(nrow != m.nrow || ncolumn != m.ncolumn) {
   if(v) delete [] v ;

   nrow = m.nrow ; ncolumn = m.ncolumn ;

 
   v = new Scalar[nrow*ncolumn] ;
 }

 // copy data
 int i;
 for(i=0; i < nrow*ncolumn; ++i)
   v[i] = m.v[i] ;

}

template<class Scalar>
void
GenFullM<Scalar>::operator=(const Scalar c)
{
 int i;
 for(i=0; i<nrow*ncolumn; ++i)
   v[i] = c;
}

template<class Scalar>
void
GenFullM<Scalar>::operator*=(const Scalar c)
{
 int i;
 for(i=0; i<nrow*ncolumn; ++i)
   v[i] *= c;
}

template<class Scalar>
void GenFullM<Scalar>::operator+=(const GenFullM<Scalar> &m)
{
    if(this->ncolumn != m.ncolumn || this->nrow != m.nrow) return;
    for(int i = 0; i < this->ncolumn * this->nrow; ++i)
        this->v[i] += m[i];
}

template<class Scalar>
GenFullM<Scalar>
GenFullM<Scalar>::operator*(GenFullM &m)
{
 if(ncolumn != m.nrow) return GenFullM(1,1) ; //error
 GenFullM res(nrow,m.ncolumn) ;
 int i,j,k;
 for(i = 0 ; i < nrow ; ++i)
  for(j=0; j < m.ncolumn ; ++j) {
    res[i][j] = 0.0;
    for(k = 0;  k < ncolumn; ++k)
      res[i][j] += (*this)[i][k] * m[k][j];
   }
 return res;
}


//*************************************************
// I added the matrix times vector operation
//*************************************************
/*template<class Scalar>
Vector 
GenFullM::operator*(const Vector &v)
{
  if(ncolumn != v.size()) return Vector(1) ; //error
  Vector res(nrow) ;
  int i,j;
  for (i = 0 ; i < nrow ; ++i) {
     res[i] = 0.0;
     for(j=0; j < ncolumn ; ++j)
       res[i] += (*this)[i][j] * v[j];
   }
  return res;
}
*/
//*************************************************


template<class Scalar>
GenFullM<Scalar>
GenFullM<Scalar>::operator%(GenFullM &m)
{
 if(ncolumn != m.ncolumn) return GenFullM(1,1) ; //error
 GenFullM res(nrow,m.nrow) ;
 int i,j,k;
 for(i=0; i<nrow; ++i)
  for(j=0; j<m.nrow; ++j) {
    res[i][j] = 0.0;
    for(k=0;  k<ncolumn; ++k)
      res[i][j] += (*this)[i][k] * m[j][k];
   }
 return res;
}

template<class Scalar>
GenFullM<Scalar>
GenFullM<Scalar>::operator^(GenFullM<Scalar> &m)
{
 if(nrow != m.nrow) return GenFullM(1,1) ; //error

 GenFullM res(ncolumn,m.ncolumn) ;
 int i,j,k;
 for(i=0; i < ncolumn; ++i)
  for(j=0; j < m.ncolumn; ++j)
   {
    res[i][j] = 0.0 ;
    for(k = 0 ;  k < nrow ; ++k)
      res[i][j] += (*this)[k][i] * m[k][j];
   }
 return res;
}

template<class Scalar>
void
GenFullM<Scalar>::invert()
{

  if (nrow != ncolumn)
    std::cout << "*** Matrix to be inverted is not square ***" << std::endl;

  GenFullM tmp(*this);
  tmp.factor();
  double *rhs = new double[ncolumn];
  for (int i=0; i < ncolumn; ++i) { //for each i solve Axi = ei
    for (int j=0; j < ncolumn; ++j) {
      if (j!=i)
        rhs[j] = 0.0;
      else
        rhs[j] = 1.0;
    }
    tmp.reSolve(rhs);
    for (int j=0; j < ncolumn; ++j)
      (*this)[j][i] = rhs[j];
  }

}
/*
template<class Scalar>
GenFullM
GenFullM::invert()
{
 if(nrow != ncolumn) return GenFullM(1,1) ; //error
 GenFullM res(*this) ;
  int i,j,k;
  for(i = 0 ; i < nrow ; ++i)
    for(j=i+1; j < nrow ; ++j) {
       Scalar p = res[j][i] = -res[j][i]/res[i][i] ;
       for(k = i+1; k < ncolumn; ++k)
         res[j][k] += p*res[i][k] ;
    }
 GenFullM inv(nrow,nrow) ;
 for(i = nrow-1 ; i >=0 ; --i)
  for(j=0; j < ncolumn; ++j)
   {
    Scalar piv = 1/res[i][i] ;
    if(j < i) inv[i][j] = res[i][j] ;
    if(j == i) inv[i][j] = 1.0 ;
    if(j > i) inv[i][j] = 0.0 ;
    for(k = i+1; k <nrow; ++k)
       inv[i][j] -= res[i][k]*inv[k][j] ;
    inv[i][j] *= piv ;
   }
 return inv ;
}

template<class Scalar>
GenFullM
GenFullM::transpose()
{
 GenFullM res(ncolumn,nrow);

 int i,j;
 for(i=0; i<nrow; ++i)
   for(j=0; j<ncolumn; ++j)
      res[j][i] = (*this)[i][j]; 

 return res; 
}
*/
template<class Scalar>
void
GenFullM<Scalar>::factor()
{
 if(this->nrow != this->ncolumn) return; // Error actually
 int i,j,k;
 for(i=0; i<nrow; ++i) {
   double invD = (*this)[i][i] = 1.0/(*this)[i][i];
   for(j=i+1; j < nrow; ++j) {
      double c = ( (*this)[j][i] *=  invD );
      for(k = i+1; k < nrow; ++k)
        (*this)[j][k] -= c*(*this)[i][k];
   }
 }
}

//CBM: Factor and ReSolve are copied from the FEM code

template<class Scalar>
void
GenFullM<Scalar>::Factor(double tol)
{ // triangular factorization of a real general matrix using Gaussian elimination with complete pivoting
  if(nrow != ncolumn) cerr << " *** WARNING: GenFullM<Scalar>::Factor(double tol), nrow != ncolumn \n";

  if(iprow) { delete [] iprow; iprow = 0; }
  if(ipcol) { delete [] ipcol; ipcol = 0; }
  iprow = new int[nrow]; // output: row permutation
  ipcol = new int[ncolumn]; // output: column permutation
  int info; // output: error flag

  // change ordering of v to be [column][row] instead of [row][column]
  Scalar *v_copy = new Scalar[nrow*ncolumn];
  for(int i=0; i<nrow*ncolumn; ++i) v_copy[i] = v[i];
  for(int i=0; i<nrow; ++i)
    for(int j=0; j<ncolumn; ++j)
      v[i*nrow+j] = v_copy[j*ncolumn+i];
  delete [] v_copy;

  F77NAME(dgecp)(nrow, v, ncolumn, tol, 0, ndef, iprow, ipcol, info);

  if(info != 0) cerr << " *** WARNING: error in dgecp, info = " << info << endl;
  else if(ndef > 0) cerr << " GenFullM is rank deficient, ndef = " << ndef << endl;
}

template<class Scalar>
void
GenFullM<Scalar>::reSolve(double *x)
{
 int i,j;
 // Forward elimination
 for(i=0; i<nrow; ++i){
   for(j=i+1; j < nrow; ++j) 
      x[j] -= (*this)[j][i]*x[i];
 }

 // Backward substitution
 for(i=nrow; i--; ) {
   for(j = i+1; j < nrow; ++j)
     x[i] -= (*this)[i][j]*x[j];
   x[i] *= (*this)[i][i];
 }
}

template<class Scalar>
void
GenFullM<Scalar>::ReSolve(double *b)
{
 // solves a real general linear system using the triangular factorization computed by Factor(int)
  int info;
  F77NAME(dgers)(nrow, 1, v, nrow, ndef, iprow, ipcol, b, nrow, info);


  if(info != 0) cerr << " *** WARNING: error in dgers, info = " << info << endl;
}

template<class Scalar>
void
GenFullM<Scalar>::print(const char *msg)
{
 if(*msg) fprintf(stderr,"%s\n",msg);
 int i,j;
 for(i = 0 ; i < nrow ; ++i) {
   for(j=0; j < ncolumn ; ++j)
     fprintf(stderr,"%e ",(*this)[i][j]) ;
   fprintf(stderr,"\n") ;
 }
}
/*
template<class Scalar>
double
GenFullM::max()
{
 double max = (*this)[0][0];
 int i,j;
 for(i = 0; i<nrow; ++i)
   for(j = 0; j<ncolumn; ++j)
     if((*this)[i][j] > max) max = (*this)[i][j];

 return max;
}
*/
template<class Scalar>
void
GenFullM<Scalar>::zero()
{
 int i;
 for(i=0; i<nrow*ncolumn; ++i)
   v[i] = 0.0;
}

template<class Scalar>
double
GenFullM<Scalar>::norm()
{
 double frobNorm = 0;
 for(int i=0; i<nrow*ncolumn; ++i)
   frobNorm += pow(v[i],2);

 frobNorm = sqrt(frobNorm);

 return frobNorm;
}

template<class Scalar>
void
GenFullM<Scalar>::add(GenFullM &mat, int fRow, int fCol)
{
  int mrow = mat.numRow();
  int mcol = mat.numCol();

  int icol,irow;
  for(icol = 0; icol < mcol; ++icol) {
    for(irow = 0; irow < mrow; ++irow) {
      (*this)[fRow+irow][fCol+icol] += mat[irow][icol];
      }
  }
}

/*
template<class Scalar>
void
GenFullM::transposeAssign(GenFullM&source)
{
 if(source.nrow != ncolumn || source.ncolumn != nrow) return; //error

 // unroll to avoid cache trashing

 int i,j;
 for(i = 0; i + 7 < ncolumn; i+= 8) {
   for(j = 0; j < nrow; ++j)  {
     double *vv = v+j*ncolumn;
     double *ww = source.v+j;
     vv[i] = ww[i*nrow]; 
     vv[i+1] = ww[(i+1)*nrow]; 
     vv[i+2] = ww[(i+2)*nrow]; 
     vv[i+3] = ww[(i+3)*nrow]; 
     vv[i+4] = ww[(i+4)*nrow]; 
     vv[i+5] = ww[(i+5)*nrow]; 
     vv[i+6] = ww[(i+6)*nrow]; 
     vv[i+7] = ww[(i+7)*nrow]; 
   }
 }
}

template<class Scalar>
void
GenFullM::transposeMult(GenFullM& m, GenFullM& res)
{

 if(nrow != m.nrow) {
   fprintf(stderr," *** ERROR: incompatible dimensions"
                  " GenFullM::tranposeMult" );
 }

 int i,j,k;
 for(i=0; i < ncolumn; ++i)
  for(j=0; j < m.ncolumn; ++j) {
    res[i][j] = 0.0 ;
    for(k = 0 ;  k < nrow ; ++k)
      res[i][j] += (*this)[k][i] * m[k][j];
   }
 
}
*/

//--------------------------------------------------------------------------------
template<class Scalar>
SymFullM<Scalar>::SymFullM()
{
  n = 0;
  v = 0;
}

template<class Scalar>
SymFullM<Scalar>::SymFullM(int nr)
{
  n = nr;
  v = new Scalar[n*(n+1)/2];
}

template<class Scalar>
SymFullM<Scalar>::SymFullM(const SymFullM &m)
{
  n = m.n;
  int vSize = n*(n+1)/2;
  v = new Scalar[vSize];
  for (int i=0; i < vSize; ++i)
    v[i] = m.v[i];
}

template<class Scalar>
SymFullM<Scalar>::~SymFullM()
{
  if (v) delete [] v;
}

template<class Scalar>
void SymFullM<Scalar>::setNewSize(int nr, double initVal)
{
  n = nr;
  delete [] v;
  int vSize = n*(n+1)/2;
  v = new Scalar[vSize];
  for (int i=0; i<vSize; ++i)
    v[i] = initVal;
}

template<class Scalar>
void SymFullM<Scalar>::operator=(const SymFullM<Scalar> &m)
{
  if (m.v == v) return ;

  if (n != m.n) {
    if (v) delete [] v ;
    n = m.n;
    v = new Scalar[n*(n+1)/2] ;
  }

  // copy data
  for (int i=0; i < n*(n+1)/2; ++i)
    v[i] = m.v[i] ;

}

template<class Scalar>
void SymFullM<Scalar>::operator=(const Scalar c)
{
  for (int i=0; i<n*(n+1)/2; ++i)
    v[i] = c;
}

template<class Scalar>
void SymFullM<Scalar>::operator*=(const Scalar c)
{

  for (int i=0; i<n*(n+1)/2; ++i)
    v[i] *= c;
}

template<class Scalar>
SymFullM<Scalar> SymFullM<Scalar>::operator*(SymFullM &m)
{
  if (n != m.n) return SymFullM(1) ; //error
  SymFullM res(n) ;
  for (int i = 0 ; i < n ; ++i) {
    for (int j=0; j <= i ; ++j) {
      res[i][j] = 0.0;
      for (int k = 0;  k <= j ; ++k)
        res[i][j] += (*this)[i][k] * m[j][k];
      for (int k = j+1;  k <= i ; ++k)
        res[i][j] += (*this)[i][k] * m[k][j];
      for (int k = i+1;  k < n ; ++k)
        res[i][j] += (*this)[k][i] * m[k][j];
    }
  }
  return res;
}

template<class Scalar>
GenFullM<Scalar> SymFullM<Scalar>::operator*(GenFullM<Scalar> &m)
{
  if (n != m.nrow) return GenFullM<Scalar>(1) ; //error
  SymFullM res(n, m.ncolumn) ;
  for (int i = 0 ; i < n ; ++i) {
    for (int j = 0; j < m.ncolumn ; ++j) {
      res[i][j] = 0.0;
      for (int k = 0;  k <= i; ++k)
        res[i][j] += (*this)[i][k] * m[k][j];
      for (int k = i+1;  k < n; ++k)
        res[i][j] += (*this)[k][i] * m[k][j];
    }
  }
  return res;
}

template<class Scalar>
GenFullM<Scalar> SymFullM<Scalar>::operator%(GenFullM<Scalar> &m)
{
  if (n != m.ncolumn) return GenFullM<Scalar>(1,1) ; //error
  GenFullM<Scalar> res(n,m.nrow) ;

  for (int i=0; i< n; ++i) {
    for (int j=0; j< m.nrow; ++j) {
      res[i][j] = 0.0;
      for (int k=0;  k<= i; ++k)
        res[i][j] += (*this)[i][k] * m[j][k];
        for (int k=0;  k< n; ++k)
          res[i][j] += (*this)[k][i] * m[j][k];
    }
  }
  return res;
}


template<class Scalar>
void SymFullM<Scalar>::invert()
{

  SymFullM tmp(*this);
  tmp.factor();
  double *rhs = new double[n];
  for (int j=0; j < n; ++j) { //for each j solve Axj = ej
    for (int i=0; i < n; ++i) {
      if (j!=i)
        rhs[i] = 0.0;
      else
        rhs[i] = 1.0;
    }
    tmp.reSolve(rhs);
    for (int i=j; i < n; ++i)
      (*this)[i][j] = rhs[i];
  }

}

template<class Scalar>
void SymFullM<Scalar>::print(const char *msg)
{
  if (*msg) fprintf(stderr,"%s\n",msg);
  for (int i=0 ; i < n ; ++i) {
    for (int j=0; j <= i  ; ++j)
      fprintf(stderr,"%e ",(*this)[i][j]) ;
    for (int j=i+1; j < n  ; ++j)
      fprintf(stderr,"%e ",(*this)[j][i]) ;
    fprintf(stderr,"\n") ;
  }
}

template<class Scalar>
void SymFullM<Scalar>::zero()
{
 for (int i=0; i< n*(n+1)/2 ; ++i)
   v[i] = 0.0;
}

template<class Scalar>
double SymFullM<Scalar>::norm()
{
 double frobNorm = 0;
 for (int i = 0; i < n; ++i) {
   frobNorm += pow((*this)[i][i],2);
   for (int j = i+1; j < n; ++j)
     frobNorm += 2.0*pow((*this)[j][i],2);
 }
 frobNorm = sqrt(frobNorm);

 return frobNorm;
}


template<class Scalar>
void SymFullM<Scalar>::factor()
{
  for (int i=0; i<n; ++i) {
    if ((*this)[i][i] <= 0.0) {
       fprintf(stderr,"Error in Cholesky factorization:negative elements on the diagonal \n");
       return;
    }
    double invD = 1.0/(*this)[i][i];
    for (int j=i+1; j < n; ++j) {
      double c = (*this)[j][i] *invD ;
      for (int k=j; k < n; ++k)
        (*this)[k][j] -= c*(*this)[k][i];
    }
    for (int k=i; k<n; ++k)
      (*this)[k][i] *= sqrt(invD);
  }  
}

template<class Scalar>
void SymFullM<Scalar>::reSolve(double *x)
{
 // Forward elimination
 for (int i=0; i<n; ++i){
   x[i] = x[i]/(*this)[i][i];
   for (int j=i+1; j < n; ++j)
     x[j] -= (*this)[j][i]*x[i];
 }
 // Backward substitution
 for (int i=n; i--; ) {
   for (int j=i+1; j < n; ++j)
     x[i] -= (*this)[j][i]*x[j];
   x[i] /= (*this)[i][i];
 }
}


