//#####################################################################
// Copyright 2002-2005, Robert Bridson, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POLYGON  
//##################################################################### 
#ifndef __POLYGON__
#define __POLYGON__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <cfloat>
namespace PhysBAM{

template<class T>
class POLYGON 
{
public:
    bool closed_polygon; // (true) closed with first and last vertex connected, (false) open - one less side
    ARRAY<VECTOR<T,2> > X; // vertex coordinates

    POLYGON();
    POLYGON(const int number_of_vertices);

    void Set_Closed_Polygon()
    {closed_polygon=true;}

    void Set_Open_Polygon()
    {closed_polygon=false;}

    void Set_Number_Of_Vertices(const int number_of_vertices)
    {X.Resize(number_of_vertices);}

//#####################################################################
    T Area();
    VECTOR<T,2> Find_Closest_Point_On_Polygon(const VECTOR<T,2>& X_point,int& side);                             
    T Distance_From_Polygon_To_Point(const VECTOR<T,2>& X_point);
    bool Inside_Polygon(const VECTOR<T,2>& X_point);
//#####################################################################
};   
//#####################################################################
}
#endif
