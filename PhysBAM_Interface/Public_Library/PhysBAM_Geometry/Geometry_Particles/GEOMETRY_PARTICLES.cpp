//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GEOMETRY_PARTICLES<TV>::
GEOMETRY_PARTICLES()
    :store_velocity(false)
{}
//#####################################################################
// Destructor 
//#####################################################################
template<class TV> GEOMETRY_PARTICLES<TV>::
~GEOMETRY_PARTICLES()
{}
//#####################################################################
template class GEOMETRY_PARTICLES<VECTOR<float,1> >;
template class GEOMETRY_PARTICLES<VECTOR<float,2> >;
template class GEOMETRY_PARTICLES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GEOMETRY_PARTICLES<VECTOR<double,1> >;
template class GEOMETRY_PARTICLES<VECTOR<double,2> >;
template class GEOMETRY_PARTICLES<VECTOR<double,3> >;
#endif 
}
