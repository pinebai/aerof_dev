//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUADRATIC
//#####################################################################
#ifndef __QUADRATIC__
#define __QUADRATIC__

#include <cmath>
#include <PhysBAM_Tools/Math_Tools/exchange.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
namespace PhysBAM{

template<class T>
class QUADRATIC:public NONLINEAR_FUNCTION<T(T)>
{
public:
    T a,b,c; // coefficients
    int roots; // number of roots, -1 indicates a=b=c=0 - always a root!
    T root1,root2; // root1 < root2

    QUADRATIC()
        :a(1),b(0),c(0),roots(1),root1(0),root2(0)
    {}

    QUADRATIC(const T a_input,const T b_input,const T c_input)
        :a(a_input),b(b_input),c(c_input),roots(0),root1(0),root2(0)
    {}

    T operator()(const T x) const PHYSBAM_OVERRIDE
    {return (a*x+b)*x+c;}

    void Compute_Roots()
    {if(a==0){
        if(b==0){
            if(c==0){roots=-1;return;} // function is identically zero - a=b=c=0 - always a root!
            else{roots=0;return;}} // when a=b=0 and c != 0, there are no roots
        else{roots=1;root1=-c/b;return;}} // when a=0 and b != 0, there is one root
    else{ // a != 0
        T d=Discriminant();
        if(d<0){roots=0;return;} // no real roots
        else if(d==0){roots=1;root1=-b/(2*a);return;} // one root
        else{ // d > 0 - two real roots
            using std::sqrt;
            T radical;if(b>0) radical=-b-sqrt(d);else radical=-b+sqrt(d);
            roots=2;root1=radical/(2*a);root2=2*c/radical;if(root1>root2) exchange(root1,root2);return;}}}

    void Compute_Roots_In_Interval(const T xmin,const T xmax)
    {Compute_Roots();
    if(roots==1){
        if(root1<xmin || root1>xmax) roots=0;
        else{roots=2;root2=root1;}}
    else if(roots==2){
        if(root2<xmin || root2>xmax) roots--;
        if(root1<xmin || root1>xmax){
            root1=root2;
            roots--;}}}

    T Discriminant() const
    {return sqr(b)-4*a*c;}

//#####################################################################
};
}
#endif
