// DistVectorOp.h
// Templated vector operations (element by element)

#pragma once

#include "DistVector.h"

class DistVectorOp {
public:
  template <class Scalar1, int dim1,
            class Scalar2,int dim2,
            class Scalar3,class Scalar4,class Scalar5,class Operation>
  static void Op(DistSVec<Scalar1,dim1>& v1,DistSVec<Scalar2,dim2>& v2,
                DistVec<Scalar3>& v3, DistVec<Scalar4>& v4,DistVec<Scalar5>& v5, const Operation& oper) {

    const DistInfo& distInfo = v1.info();

#pragma omp parallel for
    for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

      int locOffset = distInfo.subOffsetReg[iSub];
      int locLen = distInfo.subLenReg[iSub];

      for (int i = 0; i < locLen; ++i) 
      //  this->v[locOffset+i] = v2.v[locOffset+i];
        oper.Perform(v1.v[locOffset+i],v2.v[locOffset+i],
                     v3.v[locOffset+i],v4.v[locOffset+i],
                     v5.v[locOffset+i]);
    }
  }

  template <class Scalar1, int dim1,
            class Scalar2,int dim2,
            class Scalar3,class Scalar4,class Operation>
  static void Op(DistSVec<Scalar1,dim1>& v1,DistSVec<Scalar2,dim2>& v2,
                DistVec<Scalar3>& v3, DistVec<Scalar4>& v4, const Operation& oper) {

    const DistInfo& distInfo = v1.info();

#pragma omp parallel for
    for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

      int locOffset = distInfo.subOffsetReg[iSub];
      int locLen = distInfo.subLenReg[iSub];

      for (int i = 0; i < locLen; ++i) 
      //  this->v[locOffset+i] = v2.v[locOffset+i];
        oper.Perform(v1.v[locOffset+i],v2.v[locOffset+i],
                     v3.v[locOffset+i],v4.v[locOffset+i]);
    }
  }


  template <class Scalar1, int dim1,
            class Scalar2,int dim2,
            class Scalar3,class Operation>
  static void Op(DistSVec<Scalar1,dim1>& v1,DistSVec<Scalar2,dim2>& v2,
                DistVec<Scalar3>& v3, const Operation& oper) {

    const DistInfo& distInfo = v1.info();

#pragma omp parallel for
    for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

      int locOffset = distInfo.subOffsetReg[iSub];
      int locLen = distInfo.subLenReg[iSub];

      for (int i = 0; i < locLen; ++i) 
      //  this->v[locOffset+i] = v2.v[locOffset+i];
        oper.Perform(v1.v[locOffset+i],v2.v[locOffset+i],
                     v3.v[locOffset+i]);
    }
  }

};
