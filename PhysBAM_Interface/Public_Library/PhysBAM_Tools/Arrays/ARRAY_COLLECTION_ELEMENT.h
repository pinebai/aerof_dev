//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_COLLECTION_ELEMENT
//#####################################################################
#ifndef __ARRAY_COLLECTION_ELEMENT__
#define __ARRAY_COLLECTION_ELEMENT__

#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION_ELEMENT_BASE.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
namespace PhysBAM{

template<class T_ARRAY>
class ARRAY_COLLECTION_ELEMENT:public CLONEABLE<ARRAY_COLLECTION_ELEMENT<T_ARRAY>,ARRAY_COLLECTION_ELEMENT_BASE>
{
    typedef typename T_ARRAY::ELEMENT T;
public:
    T_ARRAY* array;

    ARRAY_COLLECTION_ELEMENT(T_ARRAY* input_array)
        :array(input_array)
    {}

    ARRAY_COLLECTION_ELEMENT()
    {array=new T_ARRAY();}

    ~ARRAY_COLLECTION_ELEMENT()
    {}

    void Clone_Helper(const ARRAY_COLLECTION_ELEMENT& element)
    {*array=*element.array;}

    void Clean_Memory()
    {array->Clean_Memory();}
    
    void Clear(const int p)
    {(*array)(p)=T();}

    void Clear_Range(const int start,const int end)
    {for(int i=start;i<=end;i++) (*array)(i)=T();}

    void Resize(const int new_size)
    {array->Resize(new_size);}
    
    void Resize(const int new_size,const int new_buffer_size)
    {array->Resize(new_size,new_buffer_size);}
    
    void Copy_Element(const int from,const int to)
    {(*array)(to)=(*array)(from);}

    void Copy_Element(const ARRAY_COLLECTION_ELEMENT_BASE& from_attribute,const int from,const int to)
    {(*array)(to)=(*dynamic_cast<const ARRAY_COLLECTION_ELEMENT<T_ARRAY>&>(from_attribute).array)(from);}

    void Copy_With_Offset(const ARRAY_COLLECTION_ELEMENT_BASE& from_attribute,const int offset)
    {array->Copy_With_Offset(*dynamic_cast<const ARRAY_COLLECTION_ELEMENT<T_ARRAY>&>(from_attribute).array,offset);}
    
    int Pack_Size() const
    {return array->Pack_Size();}
    
    void Pack(ARRAY_VIEW<char> buffer,int& position,const int p) const
    {array->Pack(buffer,position,p);}
    
    void Unpack(ARRAY_VIEW<const char> buffer,int& position,const int p)
    {array->Unpack(buffer,position,p);}
};
}
#endif
