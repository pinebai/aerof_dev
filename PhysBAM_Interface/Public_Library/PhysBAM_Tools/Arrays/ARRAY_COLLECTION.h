//#####################################################################
// Copyright 2008-2009, Geoffrey Irving, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_COLLECTION
//#####################################################################
#ifndef __ARRAY_COLLECTION__
#define __ARRAY_COLLECTION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION_ELEMENT.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <memory>
namespace PhysBAM{

template <class ID>
class ARRAY_COLLECTION:public CLONEABLE<ARRAY_COLLECTION<ID> >,public NONCOPYABLE
{
    template<class T_CLASS,class T_BASE> friend class CLONEABLE;
public:
    int number; // total length
    int buffer_size; // includes dead space
    ARRAY<ARRAY_COLLECTION_ELEMENT_BASE*,ID> arrays;
    HASHTABLE<std::string,ID> name_to_id_map;
    ARRAY<ID> owned_elements;
    ARRAY<int> deletion_list;

    ARRAY_COLLECTION();
    ~ARRAY_COLLECTION();

private:
    void Clone_Helper(const ARRAY_COLLECTION& collection)
    {Initialize(collection);}
public:

    int Size() const
    {return number;}

    void Clean_Memory()
    {number=buffer_size=0;name_to_id_map.Delete_All_Entries();name_to_id_map.Clean_Memory();for(ID i(1);i<=arrays.m;i++) if(arrays(i)){arrays(i)->Clean_Memory();}}

    bool operator==(const ARRAY_COLLECTION& collection) const
    {if(this==&collection) return true;
    if(arrays.m!=collection.arrays.m) return false;
    for(ID i(1);i<=arrays.m;i++) if(typeid(arrays(i))!=typeid(collection.arrays(i)) || arrays(i)!=collection.arrays(i)) return false;
    return true;}

    bool operator!=(const ARRAY_COLLECTION& collection) const
    {return !(*this==collection);}

    bool Valid_Index(const int p) const
    {return 1<=p && p<=number;}

    void Resize(const int new_size)
    {if(new_size>buffer_size){number=new_size;buffer_size=4*new_size/3+2;Reallocate_Buffer();}
    else if(new_size!=number){number=new_size;Resize_Within_Buffer();}}

    void Compact()
    {if(number!=buffer_size){buffer_size=number;Reallocate_Buffer();}}

    void Preallocate(const int max_size)
    {if(buffer_size<max_size){buffer_size=max_size;Reallocate_Buffer();}}

    int Add_Element()
    {Resize(number+1);return number;}

    ARRAY_PLUS_SCALAR<int,IDENTITY_ARRAY<> > Add_Elements(const int new_element)
    {int old_number=number;Resize(number+new_element);
    return IDENTITY_ARRAY<>(new_element)+old_number;}

    int Add_Element_From_Deletion_List()
    {return deletion_list.m?deletion_list.Pop():Add_Element();}

    void Delete_Element(const int p) 
    {Copy_Element_Helper(number--,p);Resize_Within_Buffer();}

    void Delete_All_Elements() 
    {if(number){number=0;Resize_Within_Buffer();}}

    void Add_To_Deletion_List(const int p)
    {assert(1<=p && p<=number);deletion_list.Append(p);}

    template<class T_OTHER_ARRAY> T_OTHER_ARRAY* Get_Array(const std::string array_name)
    {return const_cast<ARRAY_COLLECTION&>(*this).template Get_Array<T_OTHER_ARRAY>(name_to_id_map.Get(array_name));}

    template<class T_OTHER_ARRAY> const T_OTHER_ARRAY* Get_Array(const std::string array_name) const
    {return const_cast<ARRAY_COLLECTION&>(*this).template Get_Array<T_OTHER_ARRAY>(array_name);}

    template<class T_OTHER_ARRAY> T_OTHER_ARRAY* Get_Array(const ID array_id)
    {return dynamic_cast<ARRAY_COLLECTION_ELEMENT<T_OTHER_ARRAY>*>(Get_Array_Helper(array_id))->array;}

    template<class T_OTHER_ARRAY> const T_OTHER_ARRAY* Get_Array(const ID array_id) const
    {return const_cast<ARRAY_COLLECTION&>(*this).template Get_Array<T_OTHER_ARRAY>(array_id);}
    
    template<class T_OTHER_ARRAY> ID Add_Array(const std::string& array_name)
    {ID id=Add_Array(array_name,*new T_OTHER_ARRAY());owned_elements.Append(id);return id;}
    
    template<class T_OTHER_ARRAY> void Add_Array(const ID array_id)
    {Add_Array(array_id,new T_OTHER_ARRAY());owned_elements.Append(array_id);}

    template<class T_OTHER_ARRAY> ID Add_Array(const std::string& array_name,T_OTHER_ARRAY& array)
    {ID id=Add_Array(&array);name_to_id_map.Get_Or_Insert(array_name,id);return id;}
    
    template<class T_OTHER_ARRAY> void Add_Array(const ID array_id,T_OTHER_ARRAY* array)
    {Add_Array_Helper(array_id,new ARRAY_COLLECTION_ELEMENT<T_OTHER_ARRAY>(array));}
    
    template<class T_OTHER_ARRAY> ID Add_Array(T_OTHER_ARRAY* array)
    {return Add_Array_Helper(new ARRAY_COLLECTION_ELEMENT<T_OTHER_ARRAY>(array));}
    
    void Remove_Array(const std::string array_name)
    {ID array_id=name_to_id_map.Get_Default(array_name,ID());if(array_id!=ID()) Remove_Array(array_id);}

    int Pack_Size() const
    {int pack_size=0;for(ID i(1);i<=arrays.m;i++) pack_size+=arrays(i)->Pack_Size();return pack_size;}

    void Pack(ARRAY_VIEW<char> buffer,int& position,const int p) const
    {assert(Valid_Index(p));for(ID i(1);i<=arrays.m;i++) arrays(i)->Pack(buffer,position,p);}

    void Unpack(ARRAY_VIEW<const char> buffer,int& position,const int p)
    {assert(Valid_Index(p));for(ID i(1);i<=arrays.m;i++) arrays(i)->Unpack(buffer,position,p);}

    ID Number_Of_Arrays() const
    {return arrays.m;}

protected:
    void Resize_Within_Buffer()
    {for(ID i(1);i<=arrays.m;i++) if(arrays(i)) arrays(i)->Resize(number);}

    void Copy_Element_Helper(const int from,const int to)
    {for(ID i(1);i<=arrays.m;i++) if(arrays(i)) arrays(i)->Copy_Element(from,to);}

public:
    void Reallocate_Buffer()
    {for(ID i(1);i<=arrays.m;i++) if(arrays(i)) arrays(i)->Resize(number,buffer_size);}

    template<class T_ARRAY_COLLECTION> void
    Initialize(const ARRAY_VIEW<T_ARRAY_COLLECTION*>& elements_per_cell)
    {PHYSBAM_ASSERT(static_cast<void*>(static_cast<T_ARRAY_COLLECTION*>(0))==static_cast<ARRAY_COLLECTION*>(0)); // make sure the following cast is valid
    Initialize(ARRAY_VIEW<const ARRAY_COLLECTION* const>(elements_per_cell.Size(),reinterpret_cast<const ARRAY_COLLECTION* const*>(elements_per_cell.Get_Array_Pointer())));}

    template<class T_ARRAY_COLLECTION> void
    Initialize(const ARRAY<T_ARRAY_COLLECTION*>& elements_per_cell)
    {Initialize(static_cast<const ARRAY_VIEW<T_ARRAY_COLLECTION*>&>(elements_per_cell));}

    template<class T_ARRAY_COLLECTION,class T_ARRAYS> void
    Initialize(const ARRAY_BASE<T_ARRAY_COLLECTION*,T_ARRAYS,typename T_ARRAYS::INDEX>& elements_per_cell)
    {Initialize(elements_per_cell.array);}

    void Copy_Element(const ARRAY_COLLECTION& from_elements,const int from,const int to)
    {Copy_Element_Helper(from_elements,from,to);}

    int Append(const ARRAY_COLLECTION& source,int from)
    {Copy_Element(source,from,Add_Element());return number;}

    void Append(const ARRAY_COLLECTION& from_elements)
    {int offset=number;Add_Elements(from_elements.number);Copy_All_Elements_Helper(from_elements,offset);}

    int Take(ARRAY_COLLECTION& source,int from)
    {Append(source,from);source.Delete_Element(from);return number;}

    void Take(ARRAY_COLLECTION& source)
    {Append(source);source.Delete_All_Elements();}

//#####################################################################
    void Initialize(const ARRAY_COLLECTION& elements);
    void Initialize(ARRAY_VIEW<const ARRAY_COLLECTION* const> elements_per_cell);
    void Add_Arrays(const ARRAY_COLLECTION& elements);
    void Add_Elements_From_Deletion_List(const int count,ARRAY<int>& added_indices);
    void Delete_Elements_On_Deletion_List(const bool preserve_order=false);
private:
    void Add_Array_Helper(const ID array_id,ARRAY_COLLECTION_ELEMENT_BASE* array);
    ID Add_Array_Helper(ARRAY_COLLECTION_ELEMENT_BASE* array);
public:
    void Remove_Array(const ID array_id);
    ARRAY_COLLECTION_ELEMENT_BASE* Get_Array_Helper(const ID array_id);
    void Copy_Element_Helper(const ARRAY_COLLECTION& from_elements,const int from,const int to);
    void Copy_All_Elements_Helper(const ARRAY_COLLECTION& from_elements,const int offset);
//#####################################################################
};
}
#endif
