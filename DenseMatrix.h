#ifndef _DENSE_MATRIX_H_
#define _DENSE_MATRIX_H_

// GenFullM = Full Matrix class
//         stores an mxn matrix
//         certain member functions only work for
//         the square matrix case (nxn)

#include <Vector.h>

typedef Vec<double> Vector;

template<class Scalar>
class GenFullM {
 protected:
   int nrow;	// number of rows
   int ncolumn; // number of columns
   Scalar *v;   // pointer to matrix data
   //CBM-pivot
   int *iprow; 
   int *ipcol;
   int ndef;

 public:

   // constructors
   GenFullM<Scalar>(); // Creates an empty matrix
   GenFullM(int _nr);
   GenFullM(int _nr, int _nc);
   GenFullM(const GenFullM &,int _nr, int sr, int _nc, int sc);
   GenFullM(const GenFullM &);

   // destructor
   ~GenFullM();

   void setNewSize(int _nr, int _nc, double d=0.0);
   void setNewSize(int _nr, double d=0.0);

   // OPERATORS
   void  operator = (const GenFullM &);
   void  operator = (const Scalar c);
   void  operator *= (const Scalar c);
    void operator += (const GenFullM &);
   GenFullM<Scalar> operator *(GenFullM<Scalar>&);
   GenFullM<Scalar> operator ^(GenFullM<Scalar>&); // product A^T*B
   GenFullM<Scalar> operator %(GenFullM<Scalar>&); // product A*B^T

//   Vector operator *(const Vector &v);
   void invert();
//   GenFullM invert();
//   GenFullM transpose();

   int dim()    { return nrow;    }
   int numRow() { return nrow;    }
   int numCol() { return ncolumn; }

   Scalar *operator[](int i) const;
   Scalar* data() const { return v; }

//   double max();
   void print(const char *msg = "");
   void factor();
   void Factor(double tol=1.0e-6); //CBM-pivot
   void reSolve(double *d);
   void ReSolve(double *d); //CBM-pivot
   void zero();
   void add(GenFullM&, int, int);
   double norm(); // Frobenius norm

//   void transposeAssign(GenFullM&);

//   void transposeMult(GenFullM&, GenFullM&);
};

template<class Scalar>
inline
Scalar *
GenFullM<Scalar>::operator[](int i) const
 { return v+i*ncolumn; }

/*
class StackFullM : public GenFullM {
 public:
   StackFullM(int nr, int nc, double *data);
   ~StackFullM() { v = 0; }
};

template<class Scalar>
inline
StackFullM::StackFullM(int nr, int nc, double *data)
{
 nrow    = nr;
 ncolumn = nc;
 v       = data;
}

*/
typedef GenFullM<double> FullM;

// SymFullM = Symmetric Full Matrix class
//         stores the lower part of an nxn matrix

template<class Scalar>
class SymFullM : public GenFullM<Scalar> {

 protected:

  //v stores the bottom half of the symmetric matrix
  int n;
  Scalar *v;


 public:

  SymFullM<Scalar>(); // Creates an empty matrix
  SymFullM(int _nr);
  SymFullM(const SymFullM &);

  // destructor
  ~SymFullM();

  void setNewSize(int _nr, double d=0.0);

  // OPERATORS
  void  operator = (const SymFullM &);
  void  operator = (const Scalar c);
  void  operator *= (const Scalar c);
  SymFullM<Scalar> operator *(SymFullM<Scalar>&);
  GenFullM<Scalar> operator *(GenFullM<Scalar>&); // product A*B
  GenFullM<Scalar> operator %(GenFullM<Scalar>&); // product A*B^T

  void invert();

  int dim()    { return n;    }

  // write this operator correctly to access symmetric values
  Scalar *operator[](int i) const;


  Scalar* data() const { return v; }

//double max();
  void print(const char *msg = "");

 // write a cholesky
 void factor();
 void reSolve(double *d);

 void zero();
 double norm();

};

template<class Scalar>
inline
Scalar *
SymFullM<Scalar>::operator[](int i) const
 { return v+i*(i+1)/2; }

typedef SymFullM<double> SFullM;


#ifdef TEMPLATE_FIX
#include <DenseMatrix.C>
#endif
#endif
