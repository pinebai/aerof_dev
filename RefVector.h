#ifndef _REF_VECTOR_H_
#define _REF_VECTOR_H_

//------------------------------------------------------------------------------

template<class VecType>
class RefVec {

  VecType& refVec;

public:

  RefVec(VecType &);
  ~RefVec();

  const VecType &operator[] (int i) const {return refVec;}
  VecType &operator[] (int i) {return refVec;}
  
  int numVectors() const {return 1;}	// always 1 vector

};


//------------------------------------------------------------------------------

template<class VecType>
RefVec<VecType>::RefVec(VecType &vec) : refVec(vec)
{
}

template<class VecType>
RefVec<VecType>::~RefVec() 
{
}
#endif
