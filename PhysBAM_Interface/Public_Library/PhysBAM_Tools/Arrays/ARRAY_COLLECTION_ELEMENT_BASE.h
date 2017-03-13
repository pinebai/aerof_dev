//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_COLLECTION_ELEMENT_BASE
//#####################################################################
#ifndef __ARRAY_COLLECTION_ELEMENT_BASE__
#define __ARRAY_COLLECTION_ELEMENT_BASE__

#include <PhysBAM_Tools/Clone/CLONEABLE.h>
namespace PhysBAM{

class ARRAY_COLLECTION_ELEMENT_BASE:public CLONEABLE_ABSTRACT<ARRAY_COLLECTION_ELEMENT_BASE>
{
public:
    virtual void Clean_Memory()=0;
    virtual void Resize(const int new_size)=0;
    virtual void Resize(const int new_size,const int new_buffer_size)=0;
    virtual void Clear(const int p)=0;
    virtual void Clear_Range(const int start,const int end)=0;
    virtual void Copy_Element(const int from,const int to)=0;
    virtual void Copy_Element(const ARRAY_COLLECTION_ELEMENT_BASE& from_attribute,const int from,const int to)=0;
    virtual void Copy_With_Offset(const ARRAY_COLLECTION_ELEMENT_BASE& from_attribute,const int offset)=0;
    virtual int Pack_Size() const=0;
    virtual void Pack(ARRAY_VIEW<char> buffer,int& position,const int p) const=0;
    virtual void Unpack(ARRAY_VIEW<const char> buffer,int& position,const int p)=0;
};
}
#endif
