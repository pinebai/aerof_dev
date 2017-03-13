//#####################################################################
// Copyright 2008-2009, Geoffrey Irving, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POINT_CLOUD
//#####################################################################
#ifndef __POINT_CLOUD__
#define __POINT_CLOUD__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays_Computations/RANGE_COMPUTATIONS.h>
#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND_BASE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UNIFORM_FORWARD.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD_FORWARD.h>
#include <PhysBAM_Tools/Point_Clouds_Computations/CENTER.h>
#include <PhysBAM_Tools/Point_Clouds_Computations/EULER_STEP.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <memory>
namespace PhysBAM{

template<class TV>
class POINT_CLOUD:public CLONEABLE<POINT_CLOUD<TV> >
{
    typedef typename TV::SCALAR T;
public:
    typedef T SCALAR;
    typedef TV VECTOR_T;
    typedef ARRAY_COLLECTION<int> T_ARRAY_COLLECTION;

    T_ARRAY_COLLECTION array_collection;
    ARRAY<TV> X;
    ARRAY<int> id; //TODO: Replace dynamic one with this
    bool store_id;

    POINT_CLOUD()
        :store_id(false)
    {array_collection.Add_Array("X",X);}

    virtual ~POINT_CLOUD(){}
    
    void Store_Id(bool store=true)
    {store_id=store;if(store) array_collection.Add_Array("id",id);else array_collection.Remove_Array("id");}

    template<class T_INDICES>
    void Euler_Step(const T_INDICES& indices,const ARRAY<TV>& V,const T dt)
    {POINT_CLOUDS_COMPUTATIONS::Euler_Step(indices,X.array,V,dt);}

    TV Weighted_Center(const ARRAY<T>& weights)
    {return POINT_CLOUDS_COMPUTATIONS::Weighted_Center(X.array,weights);}

    RANGE<TV> Compute_Bounding_Box()
    {return ARRAYS_COMPUTATIONS::Compute_Range(X.array);}

    TV Centroid()
    {return ARRAYS_COMPUTATIONS::Average(X.array);};

    void Clone_Helper(const POINT_CLOUD<TV>& particles)
    {array_collection.Initialize(particles.array_collection);}

    bool operator==(const POINT_CLOUD<TV>& particles) const
    {return (this==&particles || array_collection==particles.array_collection);}

    bool operator!=(const POINT_CLOUD<TV>& particles) const
    {return !(*this==particles);}
};
}
#endif
