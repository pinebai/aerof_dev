//#####################################################################
// Copyright 2007-2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Interface for AERO-F
//##################################################################### 
#ifndef __PHYSBAM_INTERFACE__
#define __PHYSBAM_INTERFACE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include "LOCAL_LEVELSET.h"

#include <set>

namespace PhysBAM {

template<class T> class TRIANGLE_HIERARCHY;
template<class TV> class GEOMETRY_PARTICLES;
class TRIANGLE_MESH;

template<class T>
struct IntersectionResult {
  int triangleID; // -> -1 if no intersection
  T alpha; // Intersection is at alpha*edgeNode1 + (1-alpha)*edgeNode2
  T zeta[3];  // Intersection is at zeta[0]*triNode1+zeta[1]*triNode2+zeta[2]*triNode3
};

template<class T>
struct SubDInterface {
    ARRAY<int> scope;
    std::set<int> next_scope;
    ARRAY<TRIANGLE_3D<T> > triangle_list;
    TRIANGLE_MESH* scoped_triangle_mesh;
    TRIANGLE_HIERARCHY<T>* triangle_hierarchy;

    ARRAY<ARRAY<int> > candidates;

    SubDInterface() : scope(0),triangle_list(0),scoped_triangle_mesh(0),triangle_hierarchy(0),candidates(0) {}
    ~SubDInterface() {} // Leave Deletion responsibilities of scoped_triangle_mesh and triangle_hierarchy to DistPhysBAMInterface
};

template<class T>
class PhysBAMInterface {
    typedef VECTOR<T,3> TV;
    T thickness_parameter,thickness_over_two;
    LocalLevelSet *surface_levelset;

public:
    TRIANGLE_MESH& triangle_mesh;
    ARRAY<TRIANGLE_3D<T> > triangle_list;
    TRIANGLE_HIERARCHY<T>* triangle_hierarchy;

    ARRAY<SubDInterface<T> > SubD;

    GEOMETRY_PARTICLES<TV>& particles;
    PAIR<T,GEOMETRY_PARTICLES<TV>*> saved_state;

    PhysBAMInterface(TRIANGLE_MESH& triangle_mesh,GEOMETRY_PARTICLES<TV>& particles, LocalLevelSet* cs=0);

    ~PhysBAMInterface();

    void SetThickness(const T input_thickness)
    {thickness_parameter=input_thickness;thickness_over_two=(T).5*input_thickness;}

    // When particles have moved, then call the following
    // method to update the triangle hierarchy
    // pass rebuild_hierarchy=true when topology changes
    void Update(const int numLocSub,const bool rebuild_hierarchy=false);
    void UpdateScope(const int subD,const bool use_global_scope=false);
    std::set<int>& getScope(const int subD){return SubD(subD).next_scope;}

    void SaveOldState();

    bool HasCloseTriangle(const int subD,const TV position,const TV min_corner,const TV max_corner,int* index,bool* is_occluded,ARRAY<int>* cand=0);

    IntersectionResult<T> Intersect(const TV& start,const TV& end,const T thickness) const; // Rather than traversing the hierarchy here, could do union of triangles in the two boxes of the start / end nodes?
    void Intersect(const int subD,const ARRAY<TV>& node_positions,const ARRAY<bool>& occluded_node,ARRAY<TRIPLE<VECTOR<int,3>,IntersectionResult<T>,IntersectionResult<T> > >& edges_and_results) const;

    void computeSweptNodes(const int subD,const ARRAY<TV>& node_positions,const ARRAY<TV>& node_positions_initial,ARRAY<bool>& swept_node,const T dt) const;
};

}
#endif
