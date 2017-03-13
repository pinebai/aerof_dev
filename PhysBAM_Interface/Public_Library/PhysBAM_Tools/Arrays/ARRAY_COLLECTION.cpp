//#####################################################################
// Copyright 2008-2009, Geoffrey Irving, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_COLLECTION
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/SORT.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <sstream>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class ID> ARRAY_COLLECTION<ID>::
ARRAY_COLLECTION()
    :number(0),buffer_size(0)
{}
//#####################################################################
// Constructor
//#####################################################################
template<class ID> ARRAY_COLLECTION<ID>::
~ARRAY_COLLECTION()
{
    for(int i=1;i<=owned_elements.m;i++) delete arrays(owned_elements(i));
    for(int i=1;i<=arrays.m;++i) delete arrays(i);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class ID> void ARRAY_COLLECTION<ID>::
Initialize(const ARRAY_COLLECTION& elements)
{
    Clean_Memory();
    Add_Arrays(elements);
    Append(elements);
    name_to_id_map=elements.name_to_id_map;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class ID> void ARRAY_COLLECTION<ID>::
Initialize(ARRAY_VIEW<const ARRAY_COLLECTION* const> elements_per_cell)
{
    Clean_Memory();
    int total_number=0;for(int c=1;c<=elements_per_cell.Size();c++) if(elements_per_cell(c)){
        total_number+=elements_per_cell(c)->number;
        Add_Arrays(*elements_per_cell(c));} // include arrays that occur on any of the cell elements
    Preallocate(total_number);
    for(int c=1;c<=elements_per_cell.Size();c++) if(elements_per_cell(c)) Append(*elements_per_cell(c));
}
//#####################################################################
// Function Add_Arrays
//#####################################################################
template<class ID> void ARRAY_COLLECTION<ID>::
Add_Arrays(const ARRAY_COLLECTION& elements)
{
    for(ID i(1);i<=elements.arrays.m;i++) if(elements.arrays(i) && (arrays.m<i || !arrays(i))) Add_Array_Helper(i,elements.arrays(i)->Clone_Default());
}
//#####################################################################
// Function Add_Elements_From_Deletion_List
//#####################################################################
template<class ID> void ARRAY_COLLECTION<ID>::
Add_Elements_From_Deletion_List(const int count,ARRAY<int>& added_indices)
{
    added_indices.Preallocate(added_indices.Size()+count);
    int added=min(deletion_list.m,count);
    added_indices.Append_Elements(deletion_list.Pop_Elements(added));
    added_indices.Append_Elements(Add_Elements(count-added));
}
//#####################################################################
// Function Delete_Elements_On_Deletion_List
//#####################################################################
template<class ID> void ARRAY_COLLECTION<ID>::
Delete_Elements_On_Deletion_List(const bool preserve_order)
{
    Sort(deletion_list);
    if(preserve_order){
        for(int k=1;k<=deletion_list.m;k++){
            int next=k<deletion_list.m?deletion_list(k+1):number+1;
            for(int i=deletion_list(k)+1;i<next;i++) Copy_Element_Helper(i,i-k);}}
    else{
        int last=number;
        for(int k=deletion_list.m;k>=1;k--)
            Copy_Element_Helper(last--,deletion_list(k));}
    number-=deletion_list.m;
    Resize_Within_Buffer();
    deletion_list.Remove_All();
}
//#####################################################################
// Function Add_Array
//#####################################################################
template<class ID> void ARRAY_COLLECTION<ID>::
Add_Array_Helper(const ID array_id,ARRAY_COLLECTION_ELEMENT_BASE* array)
{
    array->Resize(number,buffer_size);
    arrays.Resize(max(array_id,arrays.Size()));
    if(arrays(array_id)) delete arrays(array_id);
    arrays(array_id)=array;
}
//#####################################################################
// Function Add_Array
//#####################################################################
template<class ID> ID ARRAY_COLLECTION<ID>::
Add_Array_Helper(ARRAY_COLLECTION_ELEMENT_BASE* array)
{
    ID index=Number_Of_Arrays()+1;
    Add_Array_Helper(index,array);
    return index;
}
//#####################################################################
// Function Remove_Array
//#####################################################################
template<class ID> void ARRAY_COLLECTION<ID>::
Remove_Array(const ID array_id)
{
    delete arrays(array_id);arrays(array_id)=0;
}
//#####################################################################
// Function Get_Array_Helper
//#####################################################################
template<class ID> ARRAY_COLLECTION_ELEMENT_BASE* ARRAY_COLLECTION<ID>::
Get_Array_Helper(const ID array_id)
{
    return arrays(array_id);
}
//#####################################################################
// Function Copy_Element_Helper
//#####################################################################
template<class ID> void ARRAY_COLLECTION<ID>::
Copy_Element_Helper(const ARRAY_COLLECTION& from_elements,const int from,const int to)
{
    for(ID i(1);i<=arrays.m;i++) if(arrays(i)){
        if(i<=from_elements.arrays.Size() && from_elements.arrays(i)) arrays(i)->Copy_Element(*from_elements.arrays(i),from,to);
        else arrays(i)->Clear(to);}
}
//#####################################################################
// Function Copy_All_Elements_Helper
//#####################################################################
template<class ID> void ARRAY_COLLECTION<ID>::
Copy_All_Elements_Helper(const ARRAY_COLLECTION& from_elements,const int offset)
{
    for(ID i(1);i<=arrays.m;i++) if(arrays(i)){
        if(i<=from_elements.arrays.Size() && from_elements.arrays(i)) arrays(i)->Copy_With_Offset(*from_elements.arrays(i),offset);
        else arrays(i)->Clear_Range(offset+1,offset+from_elements.number);}
}
//#####################################################################
template class ARRAY_COLLECTION<int>;
}
