//#####################################################################
// Copyright 2002-2005, Robert Bridson, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POLYGON  
//##################################################################### 
#include <PhysBAM_Geometry/Basic_Geometry/POLYGON.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> POLYGON<T>::
POLYGON() 
    :closed_polygon(true)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> POLYGON<T>::
POLYGON(const int number_of_vertices) 
    :closed_polygon(true),X(number_of_vertices)
{
}
//#####################################################################
// Function Area
//#####################################################################
// doesn't work if the polygon is self intersecting
template<class T> T POLYGON<T>::
Area()
{       
    if(!closed_polygon) return 0;
    T area=0;
    for(int k=1;k<X.m;k++) area+=(X(k).y+X(k+1).y)*(X(k+1).x-X(k).x);
    area+=(X(X.m).y+X(1).y)*(X(1).x-X(X.m).x); // last edge
    return abs(area)/2; // should have done ((*y)(k)+(*y)(k+1))/2 above
}
//#####################################################################
// Function Find_Closest_Point_On_Polygon
//#####################################################################
template<class T> VECTOR<T,2> POLYGON<T>::
Find_Closest_Point_On_Polygon(const VECTOR<T,2>& X_point,int& side)
{                  
    T distance_squared=FLT_MAX;side=0;VECTOR<T,2> result;

    // all sides except for the last one
    for(int k=1;k<X.m;k++){
        VECTOR<T,2> closest=SEGMENT_2D<T>(X(k),X(k+1)).Closest_Point_On_Segment(X_point);T d=(closest-X_point).Magnitude_Squared();
        if(distance_squared > d){distance_squared=d;result=closest;side=k;}}

    // last side - if the polygon is closed
    if(closed_polygon){
        VECTOR<T,2> closest=SEGMENT_2D<T>(X(X.m),X(1)).Closest_Point_On_Segment(X_point);T d=(closest-X_point).Magnitude_Squared();
        if(distance_squared > d){result=closest;side=X.m;}}
    return result;
}
//#####################################################################
// Function Distance_From_Polygon_To_Point
//#####################################################################
template<class T> T POLYGON<T>::
Distance_From_Polygon_To_Point(const VECTOR<T,2>& X_point)
{
    int side;return (X_point-Find_Closest_Point_On_Polygon(X_point,side)).Magnitude();
}
//#####################################################################
// Function Inside_Polygon
//#####################################################################
template<class T> bool POLYGON<T>::
Inside_Polygon(const VECTOR<T,2>& X_point)
{                  
    T theta_total=0;

    // all sides except for the last one
    for(int k=1;k<X.m;k++){
        VECTOR<T,2> X1=X(k)-X_point,X2=X(k+1)-X_point;
        if(X1==VECTOR<T,2>() || X2==VECTOR<T,2>()) return true; // (x,y) lies on the polygon
        T theta1=atan2(X1.y,X1.x),theta2=atan2(X2.y,X2.x); // atan2 returns values between -pi and pi, if (x,y) != (0,0)
        T theta=theta2-theta1;
        if(theta == (T)pi || theta == -(T)pi) return true; // (x,y) lies on the polygon
        if(theta > (T)pi) theta-=2*(T)pi;          // make sure the smaller angle is swept out
        else if(theta < -(T)pi) theta+=2*(T)pi; // make sure the smaller angle is swept out
        theta_total+=theta;}

    // last side
    VECTOR<T,2> X1=X(X.m)-X_point,X2=X(1)-X_point;
    if(X1==VECTOR<T,2>() || X2==VECTOR<T,2>()) return true; // (x,y) lies on the polygon
    T theta1=atan2(X1.y,X1.x),theta2=atan2(X2.y,X2.x); // atan2 returns values between -pi and pi, if (x,y) != (0,0)
    T theta=theta2-theta1;
    if(theta == (T)pi || theta == -(T)pi) return true; // (x,y) lies on the polygon
    if(theta > (T)pi) theta-=2*(T)pi;          // make sure the smaller angle is swept out
    else if(theta < -(T)pi) theta+=2*(T)pi; // make sure the smaller angle is swept out
    theta_total+=theta;

    // decide on inside or outside
    if(abs(theta_total) >= (T)pi)  return true; // theta_total = +2*pi or -2*pi for a point inside the polygon
    else return false;                          // theta_total = 0 for a point outside the polygon
}
//#####################################################################
template class POLYGON<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class POLYGON<double>;
#endif
