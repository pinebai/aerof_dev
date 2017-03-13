//#####################################################################
// Copyright 2007, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __POINT_SIMPLEX_MESH__
#define __POINT_SIMPLEX_MESH__

#include <PhysBAM_Geometry/Topology/SIMPLEX_MESH.h>
namespace PhysBAM{

class POINT_SIMPLEX_MESH:public SIMPLEX_MESH<0>
{
public:
    ARRAY<bool> directions; // false for left and true for right

//#####################################################################
};
}
#endif
