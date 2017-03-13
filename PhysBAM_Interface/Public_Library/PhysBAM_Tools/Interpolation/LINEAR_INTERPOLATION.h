//#####################################################################
// Copyright 2003-2005, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INTERPOLATION 
//#####################################################################
#ifndef __LINEAR_INTERPOLATION__
#define __LINEAR_INTERPOLATION__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T,class T2>
class LINEAR_INTERPOLATION
{
public:
    static T2 Linear(const T2& u_left,const T2& u_right,const T x)
    {return (1-x)*u_left+x*u_right;}

    static T2 Linear(const T2& u_left,const T2& u_right,const VECTOR<T,1>& X)
    {return Linear(u_left,u_right,X.x);}

    static T2 Linear(const T x_left,const T x_right,const T2& u_left,const T2& u_right,const T x)
    {return u_left+(x-x_left)*(u_right-u_left)/(x_right-x_left);}

    static T2 Linear(const T x_left,const T2& u_left,const T2& u_slope,const T x)
    {return u_left+(x-x_left)*u_slope;}

    static T2 Linear(const T x_left,const T2& u_left,const T2& u_slope,const VECTOR<T,1> X)
    {return u_left+(X.x-x_left)*u_slope;}

    static T2 Linear_Predivided(const T x_left,const T one_over_x_right_minus_x_left,const T2& u_left,const T2& u_right,const T x)
    {return u_left+(x-x_left)*one_over_x_right_minus_x_left*(u_right-u_left);}

    static T2 Linear_Normalized(const T2& u_left,const T2& u_slope,const T x)
    {return u_left+x*u_slope;}

    static T2 Bilinear(const T2& u1,const T2& u2,const T2& u3,const T2& u4,const VECTOR<T,2>& minimum_corner,const VECTOR<T,2>& maximum_corner,const VECTOR<T,2>& X)
    {T one_over_x_right_minus_x_left=1/(maximum_corner.x-minimum_corner.x);
        T2 u_bottom=Linear_Predivided(minimum_corner.x,one_over_x_right_minus_x_left,u1,u2,X.x),
            u_top=Linear_Predivided(minimum_corner.x,one_over_x_right_minus_x_left,u3,u4,X.x);
    return Linear(minimum_corner.y,maximum_corner.y,u_bottom,u_top,X.y);}

    // X in [0,1]x[0,1]
    static T2 Bilinear(const T2& u1,const T2& u2,const T2& u3,const T2& u4,const VECTOR<T,2>& X)
    {T2 u_bottom=Linear_Normalized(u1,u2-u1,X.x),u_top=Linear_Normalized(u3,u4-u3,X.x);
    return Linear_Normalized(u_bottom,u_top-u_bottom,X.y);}

    static T2 Bilinear(const T2& u1,const T2& u3,T one_over_y_top_minus_y_bottom,const T x_left,const T y_bottom,const T2& slope12,const T2& slope34,const VECTOR<T,2>& X)
    {T2 u_bottom=Linear(x_left,u1,slope12,X.x),u_top=Linear(x_left,u3,slope34,X.x);
    return Linear_Predivided(y_bottom,one_over_y_top_minus_y_bottom,u_bottom,u_top,X.y);}

    static T2 Trilinear(const T2& u1,const T2& u2,const T2& u3,const T2& u4,const T2& u5,const T2& u6,const T2& u7,const T2& u8,
                        const VECTOR<T,3>& minimum_corner,const VECTOR<T,3>& maximum_corner,const VECTOR<T,3>& X)
    {T one_over_x_right_minus_x_left=1/(maximum_corner.x-minimum_corner.x),one_over_y_right_minus_y_left=1/(maximum_corner.y-minimum_corner.y);
    T2 u_bottom=Linear_Predivided(minimum_corner.x,one_over_x_right_minus_x_left,u1,u2,X.x),
        u_top=Linear_Predivided(minimum_corner.x,one_over_x_right_minus_x_left,u3,u4,X.x);
    T2 u_front=Linear_Predivided(minimum_corner.y,one_over_y_right_minus_y_left,u_bottom,u_top,X.y);
    u_bottom=Linear_Predivided(minimum_corner.x,one_over_x_right_minus_x_left,u5,u6,X.x);
    u_top=Linear_Predivided(minimum_corner.x,one_over_x_right_minus_x_left,u7,u8,X.x);    
    T2 u_back=Linear_Predivided(minimum_corner.y,one_over_y_right_minus_y_left,u_bottom,u_top,X.y);
    return Linear(minimum_corner.z,maximum_corner.z,u_front,u_back,X.z);}

    // X in [0,1]x[0,1]x[0,1]
    static T2 Trilinear(const T2& u1,const T2& u2,const T2& u3,const T2& u4,const T2& u5,const T2& u6,const T2& u7,const T2& u8,const VECTOR<T,3>& X)
    {T2 u_bottom=Linear_Normalized(u1,u2-u1,X.x),u_top=Linear_Normalized(u3,u4-u3,X.x),u_front=Linear_Normalized(u_bottom,u_top-u_bottom,X.y);
    u_bottom=Linear_Normalized(u5,u6-u5,X.x);u_top=Linear_Normalized(u7,u8-u7,X.x);T2 u_back=Linear_Normalized(u_bottom,u_top-u_bottom,X.y);
    return Linear_Normalized(u_front,u_back-u_front,X.z);}

    static T2 Trilinear(const T2& u1,const T2& u3,const T2& u5,const T2& u7,T one_over_y_top_minus_y_bottom,T one_over_z_back_minus_z_front,const T x_left,const T y_bottom,const T z_front,
        const T2& slope12,const T2& slope34,const T2& slope56,const T2& slope78,const VECTOR<T,3>& X)
    {T2 u_bottom=Linear(x_left,u1,slope12,X.x),u_top=Linear(x_left,u3,slope34,X.x);
    T2 u_front=Linear_Predivided(y_bottom,one_over_y_top_minus_y_bottom,u_bottom,u_top,X.y);
    u_bottom=Linear(x_left,u5,slope56,X.x);u_top=Linear(x_left,u7,slope78,X.x);    
    T2 u_back=Linear_Predivided(y_bottom,one_over_y_top_minus_y_bottom,u_bottom,u_top,X.y);
    return Linear_Predivided(z_front,one_over_z_back_minus_z_front,u_front,u_back,X.z);}

    static T2 Linear(const T2 nodes[2],const VECTOR<T,1>& X)
    {return Linear(nodes[0],nodes[1],X.x);}

    static T2 Linear(const T2 nodes[4],const VECTOR<T,2>& X)
    {return Bilinear(nodes[0],nodes[1],nodes[2],nodes[3],X);}

    static T2 Linear(const T2 nodes[8],const VECTOR<T,3>& X)
    {return Trilinear(nodes[0],nodes[1],nodes[2],nodes[3],nodes[4],nodes[5],nodes[6],nodes[7],X);}

//#####################################################################
};
}
#endif
