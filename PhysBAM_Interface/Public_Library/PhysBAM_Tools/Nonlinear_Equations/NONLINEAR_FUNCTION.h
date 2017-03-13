//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NONLINEAR_FUNCTION
//#####################################################################
#ifndef __NONLINEAR_FUNCTION__
#define __NONLINEAR_FUNCTION__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
namespace PhysBAM{

template<class F> class NONLINEAR_FUNCTION; // F must be a function type (e.g., T(T) or T(T1,T2))
template<class T,class F> class PARAMETRIC_LINE; // F must be a two argument function type

template<class R,class T1>
class NONLINEAR_FUNCTION<R(T1)>
{
public:
    virtual ~NONLINEAR_FUNCTION(){}
//#####################################################################
    virtual R operator()(const T1 x) const=0;
    virtual R Prime(const T1 x) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual R Prime_Prime(const T1 x) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};

template<class R,class T1,class T2>
class NONLINEAR_FUNCTION<R(T1,T2)>
{
public:
    virtual ~NONLINEAR_FUNCTION() {}
//#####################################################################
    virtual R operator()(const T1 x,const T2 y) const=0;
    virtual R Partial_X(const T1 x,const T2 y) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual R Partial_Y(const T1 x,const T2 y) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};

template<class T,class R,class T1,class T2>
class PARAMETRIC_LINE<T,R(T1,T2)>:public NONLINEAR_FUNCTION<R(T)>
{
public:
    const NONLINEAR_FUNCTION<R(T1,T2)>& f;
    T1 x_not,direction_x;
    T2 y_not,direction_y;

    PARAMETRIC_LINE(const NONLINEAR_FUNCTION<R(T1,T2)>& f,const T1 x,const T2 y,const T1 a,const T2 b)
        :f(f),x_not(x),direction_x(a),y_not(y),direction_y(b)
    {}

    R operator()(const T t) const PHYSBAM_OVERRIDE
    {return f(x_not+t*direction_x,y_not+t*direction_y);}

//#####################################################################
};
}
#endif
