//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MATRIX_1X2__
#define __MATRIX_1X2__

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX_BASE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T>
class MATRIX<T,1,2>:public MATRIX_BASE<T,MATRIX<T,1,2> >
{
    struct UNUSABLE{};
public:
    typedef T SCALAR;
    enum WORKAROUND1 {m=1,n=2};
    typedef MATRIX_BASE<T,MATRIX<T,1,2> > BASE;using BASE::operator*;using BASE::Transpose_Times;using BASE::Times_Transpose;

    T x11,x12;

    explicit MATRIX(INITIAL_SIZE mm=INITIAL_SIZE(1),INITIAL_SIZE nn=INITIAL_SIZE(2))
        :x11(0),x12(0)
    {
        assert(mm==INITIAL_SIZE(1) && nn==INITIAL_SIZE(2));
    }

    MATRIX(const MATRIX& matrix)
        :x11(matrix.x11),x12(matrix.x12)
    {}

    template<class T2>
    explicit MATRIX(const MATRIX<T2,1,2>& matrix)
        :x11(matrix.x11),x12(matrix.x12)
    {}

    template<class T_MATRIX>
    explicit MATRIX(const MATRIX_BASE<T,T_MATRIX>& A)
        :x11(A(1,1)),x12(A(1,2))
    {
        assert(A.Rows()==1 && A.Columns()==2);
    }

    explicit MATRIX(const T x11,const T x12)
        :x11(x11),x12(x12)
    {}

    explicit MATRIX(const VECTOR<T,2>& v)
        :x11(v.x),x12(v.y)
    {}

    MATRIX& operator=(const MATRIX& matrix)
    {x11=matrix.x11;x12=matrix.x12;return *this;}

    int Rows() const
    {return 1;}

    int Columns() const
    {return 2;}

    T& operator()(const int i,const int j=1)
    {assert(i==1 && (j==1 || j==2));return (j==1?x11:x12);}

    const T& operator()(const int i,const int j=1) const
    {assert(i==1 && (j==1 || j==2));return (j==1?x11:x12);}

    bool Valid_Index(const int i,const int j) const
    {return i==1 && (j==1 || j==2);}

    VECTOR<T,1>& Column(const int j)
    {assert((j==1 || j==2));return *(VECTOR<T,1>*)(j==1?&x11:&x12);}

    const VECTOR<T,1>& Column(const int j) const
    {assert((j==1 || j==2));return  *(const VECTOR<T,1>*)(j==1?&x11:&x12);}

    bool operator==(const MATRIX& A) const
    {return x11==A.x11 && x12==A.x12;}

    bool operator!=(const MATRIX& A) const
    {return !(*this==A);}

    VECTOR<T,2> Column_Sum() const
    {return VECTOR<T,2>(x11,x12);}

    VECTOR<T,2> Column_Magnitudes() const
    {return VECTOR<T,2>(abs(x11),abs(x12));}

    MATRIX operator-() const
    {return MATRIX(-x11,-x12);}

    MATRIX& operator+=(const MATRIX& A)
    {x11+=A.x11;x12+=A.x12;return *this;}

    MATRIX& operator+=(const T& a)
    {x11+=a;x12+=a;return *this;}

    MATRIX& operator-=(const MATRIX& A)
    {x11-=A.x11;x12-=A.x12;return *this;}

    MATRIX& operator-=(const T& a)
    {x11-=a;x12-=a;return *this;}

    MATRIX& operator*=(const T a)
    {x11*=a;x12*=a;return *this;}

    MATRIX& operator/=(const T a)
    {x11/=a;x12/=a;return *this;}

    MATRIX operator+(const MATRIX& A) const
    {return MATRIX(x11+A.x11,x12+A.x12);}

    MATRIX operator+(const T a) const
    {return MATRIX(x11+a,x12+a);}

    MATRIX operator-(const MATRIX& A) const
    {return MATRIX(x11-A.x11,x12-A.x12);}

    MATRIX operator-(const T a) const
    {return MATRIX(x11-a,x12-a);}

    MATRIX operator*(const MATRIX<T,2>& A) const
    {return MATRIX(x11*A(1,1)+x12*A(2,1),x11*A(1,2)+x12*A(2,2));}

    MATRIX_MXN<T> operator*(const MATRIX_MXN<T>& A) const
    {assert(A.m==1);MATRIX_MXN<T> matrix(1,A.n);for(int i=1;i<=A.n;++i) matrix(1,i)=x11*A(1,i)+x12*A(2,i);return matrix;}

    MATRIX operator*(const T a) const
    {return MATRIX(a*x11,a*x12);}

    MATRIX operator/(const T a) const
    {return MATRIX(x11/a,x12/a);}

    VECTOR<T,1> operator*(const VECTOR<T,2>& v) const
    {return VECTOR<T,1>(x11*v.x+x12*v.y);}

    VECTOR<T,2> Transpose_Times(const VECTOR<T,1>& v) const
    {return VECTOR<T,2>(x11*v.x,x12*v.x);}

    static MATRIX Outer_Product(const VECTOR<T,1>& u,const VECTOR<T,2>& v)
    {return MATRIX(u.x*v.x,u.x*v.y);}

    MATRIX<T,2,1> Transposed() const
    {MATRIX<T,2,1> matrix;matrix(1,1)=x11;matrix(2,1)=x12;return matrix;}

    T Max() const
    {return max(x11,x12);}

    T Min() const
    {return min(x11,x12);}

    MATRIX<T,1> Times_Cross_Product_Matrix_Transpose(VECTOR<T,2> v) const
    {return MATRIX<T,1>(-v.y*x11+v.x*x12);}

    MATRIX<T,1> Times_Cross_Product_Matrix_Transpose_With_Symmetric_Result(VECTOR<T,2> v) const
    {return MATRIX<T,1>(-v.y*x11+v.x*x12);}

    MATRIX<T,2> Normal_Equations_Matrix() const
    {MATRIX<T,2> result;result(1,1)=x11*x11;result(1,2)=x12*x11;result(2,1)=x12*x11;result(2,2)=x12*x12;return result;}

    static MATRIX<T,1,2> Cross_Product_Matrix(const VECTOR<T,2>& v)
    {return  MATRIX<T,1,2>(-v.y,v.x);}
//#####################################################################
};

template<class T>
inline MATRIX<T,1,2> operator*(const T a,const MATRIX<T,1,2>& A)
{return A*a;}

}
#endif
