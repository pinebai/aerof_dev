//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Neil Molino, Eftychios Sifakis, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_MESH
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
SEGMENT_MESH::
SEGMENT_MESH()
    :connected_segments(0),ordered_loop_nodes(0)
{}
//#####################################################################
// Constructor
//#####################################################################
SEGMENT_MESH::
SEGMENT_MESH(const int number_nodes_input,const ARRAY<VECTOR<int,2> >& segment_list)
    :SIMPLEX_MESH<1>(number_nodes_input,segment_list),connected_segments(0),ordered_loop_nodes(0)
{}
//#####################################################################
// Constructor
//#####################################################################
SEGMENT_MESH::
SEGMENT_MESH(const SEGMENT_MESH& segment_mesh)
    :SIMPLEX_MESH<1>(segment_mesh),connected_segments(0),ordered_loop_nodes(0)
{}
//#####################################################################
// Destructor
//#####################################################################
SEGMENT_MESH::
~SEGMENT_MESH()
{
    delete connected_segments;delete ordered_loop_nodes;
}
//#####################################################################
// Function Delete_Auxiliary_Structures
//#####################################################################
void SEGMENT_MESH::
Delete_Auxiliary_Structures()
{
    SIMPLEX_MESH<1>::Delete_Auxiliary_Structures();
    delete connected_segments;connected_segments=0;delete ordered_loop_nodes;ordered_loop_nodes=0;
}
//#####################################################################
// Function Refresh_Auxiliary_Structures
//#####################################################################
void SEGMENT_MESH::
Refresh_Auxiliary_Structures()
{
    SIMPLEX_MESH<1>::Refresh_Auxiliary_Structures();
    if(connected_segments) Initialize_Connected_Segments();if(ordered_loop_nodes) Initialize_Ordered_Loop_Nodes();
}
//#####################################################################
// Function Initialize_Connected_Segments
//#####################################################################
void SEGMENT_MESH::
Initialize_Connected_Segments()
{
    delete connected_segments;connected_segments=new ARRAY<ARRAY<VECTOR<int,2> > >;
    UNION_FIND<> union_find(number_nodes);Add_Connectivity(union_find);
    ARRAY<int> buckets;
    for(int s=1;s<=elements.m;s++){
        int bucket=union_find.Find(elements(s)(1));
        int index=buckets.Find(bucket);
        if(!index){
            connected_segments->Append(ARRAY<VECTOR<int,2> >());
            buckets.Append(bucket);
            index=connected_segments->m;}
        (*connected_segments)(index).Append(elements(s));}
}
//#####################################################################
// Function Initialize_Ordered_Loop_Nodes
//#####################################################################
void SEGMENT_MESH::
Initialize_Ordered_Loop_Nodes()
{
    delete ordered_loop_nodes;ordered_loop_nodes=new ARRAY<ARRAY<int> >;

    bool created_neighbor_nodes=false,created_connected_segments=false;
    if(!neighbor_nodes){Initialize_Neighbor_Nodes();created_neighbor_nodes=true;} 
    if(!connected_segments){Initialize_Connected_Segments();created_connected_segments=true;}

    for(int i=1;i<=connected_segments->m;i++){
        for(int j=1;j<=(*connected_segments)(i).m;j++)for(int k=1;k<=2;k++)if((*neighbor_nodes)((*connected_segments)(i)(j)(k)).m!=2) goto not_closed;
        ordered_loop_nodes->Append(ARRAY<int>());
        {int start_node=(*connected_segments)(i)(1)(1),curr_node=(*connected_segments)(i)(1)(2);
        ordered_loop_nodes->Last().Append(start_node);
        do{
            int previous_node=ordered_loop_nodes->Last().Last(),neighbor1=(*neighbor_nodes)(curr_node)(1),neighbor2=(*neighbor_nodes)(curr_node)(2);
            ordered_loop_nodes->Last().Append(curr_node);
            curr_node=previous_node==neighbor1?neighbor2:neighbor1;}
        while(curr_node!=start_node);}
      not_closed: continue;}

    if(created_neighbor_nodes){delete neighbor_nodes;neighbor_nodes=0;}
    if(created_connected_segments){delete connected_segments;connected_segments=0;}
}
//#####################################################################
// Function Initialize_Straight_Mesh
//#####################################################################
void SEGMENT_MESH::
Initialize_Straight_Mesh(const int number_of_points,bool loop)
{
    Clean_Memory();
    number_nodes=number_of_points;
    elements.Exact_Resize(number_of_points-!loop);
    for(int i=1;i<number_of_points;i++)elements(i).Set(i,i+1);
    if(loop)elements(number_of_points).Set(number_of_points,1);
}
//#####################################################################
// Function Assert_Consistent
//#####################################################################
bool SEGMENT_MESH::
Assert_Consistent() const
{
    if(connected_segments && ordered_loop_nodes) assert(connected_segments->m==ordered_loop_nodes->m);
    return SIMPLEX_MESH<1>::Assert_Consistent();
}
//#####################################################################
