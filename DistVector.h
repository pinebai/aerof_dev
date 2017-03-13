#ifndef _DIST_VECTOR_H_
#define _DIST_VECTOR_H_

#include <alloca.h>
#include <Vector.h>
#include <DistInfo.h>
#include <complex>
#include <cmath>

typedef std::complex<double> bcomp;
//------------------------------------------------------------------------------

template<class Scalar>
class DistVec : public Vec<Scalar> {

public:

  typedef DistInfo InfoType;

private:

  const DistInfo &distInfo;

  Vec<Scalar> **subVec;

public:
  
  DistVec(const DistInfo &);
  DistVec(const DistVec<Scalar> &);
  DistVec(const DistInfo &, Scalar *);
  ~DistVec();

  template<class T>
  DistVec<Scalar> &operator=(const T);

  template<class T>
  DistVec<Scalar> &operator*=(const T);

  DistVec<Scalar> &operator=(const Scalar&);
  DistVec<Scalar> &operator=(const DistVec<Scalar> &);
  DistVec<Scalar> &operator+=(const DistVec<Scalar> &);
  DistVec<Scalar> &operator-=(const DistVec<Scalar> &);

  Scalar operator*(const DistVec<Scalar> &);

  DistVec<Scalar> &conjugate(const DistVec<Scalar> &);
  DistVec<Scalar> &pow(const DistVec<Scalar> &, double);

  template<class T>
  DistVec<Scalar> &operator=(const Expr<T, Scalar> &);

  template<class T>
  DistVec<Scalar> &operator+=(const Expr<T, Scalar> &);

  template<class T>
  DistVec<Scalar> &operator-=(const Expr<T, Scalar> &);
  
  template<class T>
  Scalar operator*(const Expr<T, Scalar> &);

  Vec<Scalar> &operator() (int i) const { return *subVec[i]; }

  Vec<Scalar>* operator[] (int i) const { return subVec[i]; }

  void createSubVec();

  DistVec<Scalar> *alias() const;

  Scalar sum() const;

  Scalar min() const;

  Scalar max() const;

  double norm()  { return sqrt(*this * *this); }

  const InfoType &info() const { return distInfo; }

  int numLocSub() const { return distInfo.numLocSub; }

  int numGlobSub() const { return distInfo.numGlobSub; }

  int subSize(int isub) const { return distInfo.subLen[isub]; }

  Scalar *subData(int isub) const { return this->v+distInfo.subOffset[isub]; }

  bool *getMasterFlag(int i) const { return distInfo.getMasterFlag(i); }

  void zeroNonMaster(); 

  // Adam April 2010 : 
  // nullifyPointers has to be called is the object is a DistVec<Scalar *>
  // otherwise call nullify and Scalar::nullify must be defined.
  void nullifyPointers() 
  {
    for(int i=0;i<this->len;++i) 
      {
        this->v[i]=0;
      }
  }
  void nullify() 
  {
    for(int i=0;i<this->len;++i) this->v[i].nullify();
  }
  void deletePointers() 
  {
    for(int i=0;i<this->len;++i) 
      {
	delete this->v[i];
	this->v[i] = 0;
      }
  }
};

//------------------------------------------------------------------------------

template <class Scalar>
DistVec<Scalar>::DistVec(const DistInfo &dI) : 
  Vec<Scalar>(dI.totLen), distInfo(dI), subVec(0)
{
  createSubVec();
}

//------------------------------------------------------------------------------

template <class Scalar>
DistVec<Scalar>::DistVec(const DistVec<Scalar> &v2) :
  Vec<Scalar>(v2.size()), distInfo(v2.distInfo), subVec(0)
{
#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

     int locOffset = distInfo.subOffsetReg[iSub];
     int locLen = distInfo.subLenReg[iSub];

     for (int i = 0; i < locLen; ++i) 
       this->v[locOffset+i] = v2.v[locOffset+i];
  }
  createSubVec();
}

//------------------------------------------------------------------------------

template<class Scalar>
DistVec<Scalar>::DistVec(const DistInfo &dI, Scalar *vv) : 
  Vec<Scalar>(dI.totLen, vv), distInfo(dI), subVec(0)
{
  createSubVec();
}

//------------------------------------------------------------------------------

template <class Scalar>
DistVec<Scalar>::~DistVec() 
{ 
  if (subVec) {
    //#pragma omp parallel for BUG alloc
    for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub)
      if (subVec[iSub]) delete subVec[iSub];
    delete [] subVec;
  }
}

//------------------------------------------------------------------------------

template <class Scalar>
inline
void 
DistVec<Scalar>::createSubVec()
{
  subVec = new Vec<Scalar>*[distInfo.numLocSub];

  //#pragma omp parallel for BUG alloc
  for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub)
    subVec[iSub] = new Vec<Scalar>(distInfo.subLen[iSub], this->v+distInfo.subOffset[iSub]);
}

//------------------------------------------------------------------------------

template<class Scalar>
inline
DistVec<Scalar> *
DistVec<Scalar>::alias() const
{
  return new DistVec<Scalar>(distInfo, this->v); 
}

//------------------------------------------------------------------------------

template<class Scalar>
template<class T>
inline
DistVec<Scalar> &
DistVec<Scalar>::operator=(const T y)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = distInfo.subOffsetReg[iSub];
    int locLen = distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) 
      this->v[locOffset+i] = y;

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
DistVec<Scalar> &
DistVec<Scalar>::operator=(const Scalar& y)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = distInfo.subOffsetReg[iSub];
    int locLen = distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) 
      this->v[locOffset+i] = y;

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar>
template<class T>
inline
DistVec<Scalar> &
DistVec<Scalar>::operator*=(const T y)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = distInfo.subOffsetReg[iSub];
    int locLen = distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) 
      this->v[locOffset+i] *= y;

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
DistVec<Scalar> &
DistVec<Scalar>::operator=(const DistVec<Scalar> &v2)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = distInfo.subOffsetReg[iSub];
    int locLen = distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) 
      this->v[locOffset+i] = v2.v[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
DistVec<Scalar> &
DistVec<Scalar>::operator+=(const DistVec<Scalar> &v2)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = distInfo.subOffsetReg[iSub];
    int locLen = distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) 
      this->v[locOffset+i] += v2.v[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
DistVec<Scalar> &
DistVec<Scalar>::operator-=(const DistVec<Scalar> &v2)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub)  {

    int locOffset = distInfo.subOffsetReg[iSub];
    int locLen = distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) 
      this->v[locOffset+i] -= v2.v[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
DistVec<Scalar> &
DistVec<Scalar>::conjugate(const DistVec<Scalar> &v2)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub)  {

    int locOffset = distInfo.subOffsetReg[iSub];
    int locLen = distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i)
      this->v[locOffset+i] = conj(v2.v[locOffset+i]);

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
DistVec<Scalar> & DistVec<Scalar>::pow(const DistVec<Scalar> &v2, double exp)
{
  using std::pow;
#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub)  {

    int locOffset = distInfo.subOffsetReg[iSub];
    int locLen = distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i)
      this->v[locOffset+i] = pow(v2.v[locOffset+i],exp);

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
Scalar 
DistVec<Scalar>::operator*(const DistVec<Scalar> &x)
{

  int iSub;

  Scalar res = 0.0;

#ifndef MPI_OMP_REDUCTION
  Scalar *allres = reinterpret_cast<Scalar *>(alloca(sizeof(Scalar) * distInfo.numGlobSub));

  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) allres[iSub] = 0;
#endif

  if (distInfo.masterFlag) {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
    for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

      int locOffset = distInfo.subOffset[iSub];
      int locLen = distInfo.subLen[iSub];

      Scalar locres = 0.0;

      for (int i = 0; i < locLen; ++i)
        if (distInfo.masterFlag[locOffset+i])
          locres += this->v[locOffset+i] * x.v[locOffset+i];

#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif

    }

  } 
  else {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
    for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

      int locOffset = distInfo.subOffset[iSub];
      int locLen = distInfo.subLen[iSub];

      Scalar locres = 0.0;

      for (int i = 0; i < locLen; ++i)
        locres += this->v[locOffset+i] * x.v[locOffset+i];

#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif
      
    }

  }

#ifdef MPI_OMP_REDUCTION
  distInfo.com->globalSum(1, &res);
#else
  distInfo.com->globalSum(distInfo.numGlobSub, allres);

  res = 0;
  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) res += allres[iSub];
#endif

  return res;

}

//------------------------------------------------------------------------------

template<class Scalar>
template<class T>
inline
DistVec<Scalar> &
DistVec<Scalar>::operator=(const Expr<T, Scalar> &expr)
{

  const T &x = expr.x;

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = distInfo.subOffsetReg[iSub];
    int locLen = distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) 
      this->v[locOffset+i] = x[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar>
template<class T>
DistVec<Scalar> &
DistVec<Scalar>::operator+=(const Expr<T, Scalar> &expr)
{

  const T &x = expr.x;

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = distInfo.subOffsetReg[iSub];
    int locLen = distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) 
      this->v[locOffset+i] += x[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar>
template<class T>
inline
DistVec<Scalar> &
DistVec<Scalar>::operator-=(const Expr<T, Scalar> &expr)
{

  const T &x = expr.x;

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = distInfo.subOffsetReg[iSub];
    int locLen = distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) 
      this->v[locOffset+i] -= x[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------
    
template<class Scalar>
template<class T>
inline
Scalar 
DistVec<Scalar>::operator*(const Expr<T, Scalar> &expr)
{

  int iSub;

  Scalar res = 0;

  const T &x = expr.x;

#ifndef MPI_OMP_REDUCTION
  Scalar *allres = reinterpret_cast<Scalar *>(alloca(sizeof(Scalar) * distInfo.numGlobSub));

  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) allres[iSub] = 0;
#endif

  if (distInfo.masterFlag) {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
    for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

      int locOffset = distInfo.subOffset[iSub];
      int locLen = distInfo.subLen[iSub];

      Scalar locres = 0;

      for (int i = 0; i < locLen; ++i)
	if (distInfo.masterFlag[locOffset+i])
	  locres += this->v[locOffset+i] * x[locOffset+i];

#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif

    }

  } 
  else {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
    for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

      int locOffset = distInfo.subOffset[iSub];
      int locLen = distInfo.subLen[iSub];

      Scalar locres = 0;

      for (int i = 0; i < locLen; ++i)
	locres += this->v[locOffset+i] * x[locOffset+i];

#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif
      
    }

  }

#ifdef MPI_OMP_REDUCTION
  distInfo.com->globalSum(1, &res);
#else
  distInfo.com->globalSum(distInfo.numGlobSub, allres);

  res = 0;
  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) res += allres[iSub];
#endif

  return res;

}

//------------------------------------------------------------------------------
    
template<class Scalar>
inline
Scalar 
DistVec<Scalar>::sum() const
{
  
  int iSub;

  Scalar res = 0;

#ifndef MPI_OMP_REDUCTION
  Scalar *allres = reinterpret_cast<Scalar *>(alloca(sizeof(Scalar) * distInfo.numGlobSub));

  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) allres[iSub] = 0;
#endif

  if (distInfo.masterFlag) {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
    for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

      int locOffset = distInfo.subOffset[iSub];
      int locLen = distInfo.subLen[iSub];

      Scalar locres = 0;

      for (int i = 0; i < locLen; ++i)
	if (distInfo.masterFlag[locOffset+i])
	  locres += this->v[locOffset+i];

#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif

    }

  } 
  else {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
    for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

      int locOffset = distInfo.subOffset[iSub];
      int locLen = distInfo.subLen[iSub];

      Scalar locres = 0;

      for (int i = 0; i < locLen; ++i)
	locres += this->v[locOffset+i];

#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif

    }

  }

#ifdef MPI_OMP_REDUCTION
  distInfo.com->globalSum(1, &res);
#else
  distInfo.com->globalSum(distInfo.numGlobSub, allres);

  res = 0;
  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) res += allres[iSub];
#endif

  return res;

}

//------------------------------------------------------------------------------
    
template<class Scalar>
inline
Scalar 
DistVec<Scalar>::min() const
{
  
  int iSub;

  Scalar *allmin = reinterpret_cast<Scalar *>(alloca(sizeof(Scalar) * distInfo.numLocSub));

#pragma omp parallel for
  for (iSub = 0; iSub < distInfo.numLocSub; ++iSub)
    allmin[iSub] = subVec[iSub]->min();

  Scalar vmin = allmin[0];
  for (iSub = 1; iSub < distInfo.numLocSub; ++iSub)
    vmin = vmin < allmin[iSub] ? vmin : allmin[iSub];

  distInfo.com->globalMin(1, &vmin);

  return vmin;

}

//------------------------------------------------------------------------------
    
template<class Scalar>
inline
Scalar 
DistVec<Scalar>::max() const
{
  
  int iSub;

  Scalar *allmax = reinterpret_cast<Scalar *>(alloca(sizeof(Scalar) * distInfo.numLocSub));

#pragma omp parallel for
  for (iSub = 0; iSub < distInfo.numLocSub; ++iSub)
    allmax[iSub] = subVec[iSub]->max();

  Scalar vmax = allmax[0];
  for (iSub = 1; iSub < distInfo.numLocSub; ++iSub)
    vmax = vmax > allmax[iSub] ? vmax : allmax[iSub];

  distInfo.com->globalMax(1, &vmax);

  return vmax;

}

//------------------------------------------------------------------------------

template<class Scalar>
inline
void
DistVec<Scalar>::zeroNonMaster() 
{

  int iSub;

#pragma omp parallel for
    for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

      int locOffset = distInfo.subOffset[iSub];
      int locLen = distInfo.subLen[iSub];

      for (int i = 0; i < locLen; ++i)
        if (!distInfo.masterFlag[locOffset+i])
           this->v[locOffset+i] = 0;
   }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
class DistSVec : public SVec<Scalar,dim> {

public:

  typedef DistInfo InfoType;

private:

  const DistInfo &distInfo;

  SVec<Scalar,dim> **subVec;

public:

  DistSVec(const DistInfo &);
  DistSVec(const DistSVec<Scalar,dim> &);
  DistSVec(const DistInfo &, Scalar (*)[dim]);
  ~DistSVec();

  DistSVec<Scalar,dim> &operator=(const Scalar);
  DistSVec<Scalar,dim> &operator*=(const Scalar);
    DistSVec<Scalar, dim> &operator/=(const Scalar);
  DistSVec<Scalar,dim> &operator/=(const DistVec<Scalar> &);
  DistSVec<Scalar,dim> &operator*=(const DistVec<Scalar> &);
  DistSVec<Scalar,dim> operator* (const Scalar);
  DistSVec<Scalar,dim> &operator+=(const Scalar);


  DistSVec<Scalar,dim> &operator=(const DistSVec<Scalar,dim> &);
  DistSVec<Scalar,dim> &operator+=(const DistSVec<Scalar,dim> &);
  DistSVec<Scalar,dim> &operator-=(const DistSVec<Scalar,dim> &);
  DistSVec<Scalar,dim> &operator*=(const DistSVec<Scalar,dim> &);
  DistSVec<Scalar,dim> &operator=(const DistVec<Scalar> &);

  DistSVec<Scalar,dim> &conjugate();
  DistSVec<Scalar,dim> &randVec();
  DistSVec<Scalar,dim> &getReal(const DistSVec<bcomp,dim> &);
  DistSVec<Scalar,dim> &getImag(const DistSVec<bcomp,dim> &);
  DistSVec<Scalar,dim> &getAbs(const DistSVec<bcomp,dim> &);

  Scalar operator*(const DistSVec<Scalar,dim> &);

  void min(Scalar vmin[dim]) const;
  
  void max(Scalar vmax[dim]) const;

  template<class T>
  DistSVec<Scalar,dim> &operator=(const Expr<T, Scalar> &);

  template<class T>
  DistSVec<Scalar,dim> &operator+=(const Expr<T, Scalar> &);

  template<class T>
  DistSVec<Scalar,dim> &operator-=(const Expr<T, Scalar> &);

  template<class T>
  Scalar operator*(const Expr<T, Scalar> &);

  SVec<Scalar,dim> &operator() (int i) const { return *subVec[i]; }

  template<int dim1, int dim2>
  void split(DistSVec<Scalar,dim1> &, DistSVec<Scalar,dim2> &);

  template<int dim1, int dim2>
  void merge(DistSVec<Scalar,dim1> &, DistSVec<Scalar,dim2> &);

  template<int dim1>
  void strip(DistSVec<Scalar,dim1> &);

  template<int dim1>
  void pad(DistSVec<Scalar,dim1> &);

  void createSubVec();

  void set(const Scalar *);
  
  void restrict();

  void average();

  int nonOverlapSize() const;

  DistSVec<Scalar,dim> *alias() const;

  const InfoType &info() const { return distInfo; }

  int numLocSub() const { return distInfo.numLocSub; }

  int numGlobSub() const { return distInfo.numGlobSub; }

  int subSize(int i) const { return distInfo.subLen[i]; }

  Scalar (*subData(int i) const)[dim] { return this->v+distInfo.subOffset[i]; }

  bool *getMasterFlag(int i) const { return distInfo.getMasterFlag(i); }

  double norm();

  void sum(Scalar sumres[dim]) const;
};

//------------------------------------------------------------------------------

template<class Scalar, int dim>
DistSVec<Scalar,dim>::DistSVec(const DistInfo &dI) : 
  SVec<Scalar,dim>(dI.totLen), distInfo(dI), subVec(0)
{

  createSubVec();

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
DistSVec<Scalar,dim>::DistSVec(const DistSVec<Scalar,dim> &y) :
  SVec<Scalar,dim>(y.size()), distInfo(y.distInfo), subVec(0)
{

  const Scalar *yy = reinterpret_cast<Scalar *>(y.v);
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] = yy[locOffset+i];

  }

  createSubVec();

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
DistSVec<Scalar,dim>::DistSVec(const DistInfo &dI, Scalar (*vv)[dim]) :
  SVec<Scalar,dim>(dI.totLen, vv), distInfo(dI), subVec(0)
{

  createSubVec();

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
DistSVec<Scalar,dim>::~DistSVec() 
{ 

  if (subVec) {
    //#pragma omp parallel for BUG alloc
    for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub)
      if (subVec[iSub]) delete subVec[iSub];
    delete [] subVec;
  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
void 
DistSVec<Scalar,dim>::createSubVec()
{

  subVec = new SVec<Scalar,dim>*[distInfo.numLocSub];

  //#pragma omp parallel for BUG alloc
  for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub)
    subVec[iSub] = new SVec<Scalar,dim>(distInfo.subLen[iSub], this->v+distInfo.subOffset[iSub]);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
void
DistSVec<Scalar,dim>::set(const Scalar *y)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = distInfo.subOffsetReg[iSub];
    int locLen = distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) 
      for (int j = 0; j < dim; ++j)
	this->v[locOffset+i][j] = y[j];
    
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
double
DistSVec<Scalar,dim>::norm()
{
  int iSub;
  double res = 0;

#ifndef MPI_OMP_REDUCTION
  double *allres = reinterpret_cast<double *>(alloca(sizeof(double) * distInfo.numGlobSub));

  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) allres[iSub] = 0;
#endif

  if (distInfo.masterFlag) {
#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
    for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

      int locOffset = distInfo.subOffset[iSub];
      int locLen = distInfo.subLen[iSub];
      double locres = 0;

      for (int i = 0; i < locLen; ++i)
        if (distInfo.masterFlag[locOffset+i])
          for (int j = 0; j < dim; ++j)
            locres += sqNorm(this->v[locOffset+i][j]);

#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif
    }
  }
  else {
#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
    for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub) {
      int locOffset = distInfo.subOffset[iSub];
      int locLen = distInfo.subLen[iSub];
      double locres = 0;

      for (int i = 0; i < locLen; ++i)
        for (int j = 0; j < dim; ++j)
          locres += sqNorm(this->v[locOffset+i][j]);

#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif
    }
  }

#ifdef MPI_OMP_REDUCTION
  distInfo.com->globalSum(1, &res);
#else
  distInfo.com->globalSum(distInfo.numGlobSub, allres);
  res = 0;
  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) res += allres[iSub];
#endif

  return sqrt(res);
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> *
DistSVec<Scalar,dim>::alias() const
{

  return new DistSVec<Scalar,dim>(distInfo, this->v); 

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim>
DistSVec<Scalar,dim>::operator*(const Scalar y)
{

  DistSVec<Scalar, dim> result(*this);
  result *= y;

  return result;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator+=(const Scalar y)
{

  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] += y;

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::getImag(const DistSVec<bcomp,dim> &y)
{

  const bcomp *yy = reinterpret_cast<bcomp *>(y.v);
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] = imag(yy[locOffset+i]);

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator=(const Scalar y)
{

  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] = y;

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator*=(const Scalar y)
{

  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] *= y;

  }

  return *this;

}

template<class Scalar, int dim>
inline
DistSVec<Scalar, dim> &
DistSVec<Scalar, dim>::operator/=(const Scalar y)
{
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#pragma omp parallel for
  for(int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) vv[locOffset + i] /= y;
  }

  return *this;

};

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator/=(const DistVec<Scalar> &y)
{

  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = distInfo.subLenReg[iSub];
    double *ySubData = y.subData(iSub);

    if (y.subSize(iSub) != locLen)
      fprintf(stderr, "WARNING %d vs %d\n", y.subSize(iSub), locLen);

    // Suppressed by Adam 2010.07.30
    //    InfoType yInfo = y.info();

    for (int i = 0; i < locLen; ++i)
      for (int j = 0; j < dim; j++)
        vv[locOffset+i*dim+j] /= ySubData[i];
  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator=(const DistSVec<Scalar,dim> &y)
{

  const Scalar *yy = reinterpret_cast<Scalar *>(y.v);
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] = yy[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator=(const DistVec<Scalar> &y)
{

  const Scalar *yy = y.data();

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = distInfo.subOffsetReg[iSub];
    int locLen = distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) 
      for (int j = 0; j < dim; ++j)
	this->v[locOffset+i][j] = yy[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator*=(const DistVec<Scalar> &y)
{

  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = distInfo.subLenReg[iSub];
    Scalar *ySubData = reinterpret_cast<Scalar *>(y.subData(iSub));

    for (int i = 0; i < locLen; ++i)
      for (int j = 0; j < dim; j++)
        vv[locOffset+i*dim+j] *= ySubData[i];
  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator+=(const DistSVec<Scalar,dim> &y)
{

  const Scalar *yy = reinterpret_cast<Scalar *>(y.v);
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] += yy[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator*=(const DistSVec<Scalar,dim> &y)
{

  const Scalar *yy = reinterpret_cast<Scalar *>(y.v);
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] *= yy[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator-=(const DistSVec<Scalar,dim> &y)
{

  const Scalar *yy = reinterpret_cast<Scalar *>(y.v);
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] -= yy[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------
    
template<class Scalar, int dim>
inline
Scalar 
DistSVec<Scalar,dim>::operator*(const DistSVec<Scalar,dim> &x)
{

  int iSub;

  Scalar res = 0;

#ifndef MPI_OMP_REDUCTION
  Scalar *allres = reinterpret_cast<Scalar *>(alloca(sizeof(Scalar) * distInfo.numGlobSub));

  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) allres[iSub] = 0;
#endif

  if (distInfo.masterFlag) {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
    for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

      int locOffset = distInfo.subOffset[iSub];
      int locLen = distInfo.subLen[iSub];

      Scalar locres = 0;

      for (int i = 0; i < locLen; ++i)
        if (distInfo.masterFlag[locOffset+i])
          for (int j = 0; j < dim; ++j)
            locres += DotTerm(this->v[locOffset+i][j],x.v[locOffset+i][j]);

#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif

    }

  } 

  else {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
    for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

      int locOffset = distInfo.subOffset[iSub];
      int locLen = distInfo.subLen[iSub];

      Scalar locres = 0;

      for (int i = 0; i < locLen; ++i)
        for (int j = 0; j < dim; ++j)
          locres += DotTerm(this->v[locOffset+i][j], x.v[locOffset+i][j]);
      
#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif

    }

  }

#ifdef MPI_OMP_REDUCTION
  distInfo.com->globalSum(1, &res);
#else
  distInfo.com->globalSum(distInfo.numGlobSub, allres);

  res = 0;
  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) res += allres[iSub];
#endif

  return res;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::conjugate()
{

  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] = conj(vv[locOffset+i]);

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::getReal(const DistSVec<bcomp,dim> &y)
{

  const bcomp *yy = reinterpret_cast<bcomp *>(y.v);
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] = real(yy[locOffset+i]);

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::getAbs(const DistSVec<bcomp,dim> &y)
{

  const bcomp *yy = reinterpret_cast<bcomp *>(y.v);
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] = std::abs(yy[locOffset+i]);

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::randVec()
{

  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] = rand();

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<class T>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator=(const Expr<T, Scalar> &expr)
{

  const T &x = expr.x;
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] = x[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<class T>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator+=(const Expr<T, Scalar> &expr)
{

  const T &x = expr.x;
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] += x[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<class T>
inline
DistSVec<Scalar,dim> &
DistSVec<Scalar,dim>::operator-=(const Expr<T, Scalar> &expr)
{

  const T &x = expr.x;
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) vv[locOffset+i] -= x[locOffset+i];

  }

  return *this;

}

//------------------------------------------------------------------------------
    
template<class Scalar, int dim>
template<class T>
inline
Scalar 
DistSVec<Scalar,dim>::operator*(const Expr<T, Scalar> &expr)
{

  int iSub;

  Scalar res = 0;

  const T &x = expr.x;
  Scalar *vv = reinterpret_cast<Scalar *>(this->v);

#ifndef MPI_OMP_REDUCTION
  Scalar *allres = reinterpret_cast<Scalar *>(alloca(sizeof(Scalar) * distInfo.numGlobSub));

  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) allres[iSub] = 0;
#endif

  if (distInfo.masterFlag) {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
    for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

      int locOffset = distInfo.subOffset[iSub];
      int locLen = distInfo.subLen[iSub];

      Scalar locres = 0;

      for (int i = 0; i < locLen; ++i)
	if (distInfo.masterFlag[locOffset+i])
	  for (int j = 0; j < dim; ++j)
	    locres += vv[(locOffset+i)*dim + j] * x[(locOffset+i)*dim + j];

#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif

    }

  } 
  else {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: res)
#else
#pragma omp parallel for
#endif
    for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

      int locOffset = distInfo.subOffset[iSub];
      int locLen = distInfo.subLen[iSub];

      Scalar locres = 0;

      for (int i = 0; i < locLen; ++i)
	for (int j = 0; j < dim; ++j)
	  locres += vv[(locOffset+i)*dim + j] * x[(locOffset+i)*dim + j];
      
#ifdef MPI_OMP_REDUCTION
      res += locres;
#else
      allres[distInfo.locSubToGlobSub[iSub]] = locres;
#endif

    }

  }

#ifdef MPI_OMP_REDUCTION
  distInfo.com->globalSum(1, &res);
#else
  distInfo.com->globalSum(distInfo.numGlobSub, allres);

  res = 0;
  for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) res += allres[iSub];
#endif

  return res;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
int 
DistSVec<Scalar,dim>::nonOverlapSize() const
{

  int tsize = 0;

  if (distInfo.masterFlag) {

#pragma omp parallel for reduction(+: tsize)
    for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

      for (int i=0; i<distInfo.subLen[iSub]; ++i)
	if (distInfo.masterFlag[ distInfo.subOffset[iSub]+i ]) ++tsize;

    }
    distInfo.com->globalSum(1, &tsize);

  } 
  else
    tsize = distInfo.totLen;

  return tsize;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
void
DistSVec<Scalar,dim>::restrict()
{

  if (!distInfo.masterFlag) return;

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub) {
 
    int locOffset = distInfo.subOffset[iSub];
    
    for (int i=0; i<distInfo.subLen[iSub]; ++i)
      if (!distInfo.masterFlag[locOffset+i])
	for (int j=0; j<dim; ++j) this->v[locOffset+i][j] = 0;
      
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
void
DistSVec<Scalar,dim>::average()
{
   if (!distInfo.invNdWeight) return;
#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

    int locOffset = distInfo.subOffset[iSub];
   
    for (int i=0; i<distInfo.subLen[iSub]; ++i) {
       for (int j=0; j<dim; ++j) 
          this->v[locOffset+i][j] *= distInfo.invNdWeight[locOffset+i];
    }
   }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<int dim1, int dim2>
inline
void
DistSVec<Scalar,dim>::split(DistSVec<Scalar,dim1> &y, DistSVec<Scalar,dim2> &z)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub) {
 
    int locOffset = distInfo.subOffset[iSub];
    
    for (int i=0; i<distInfo.subLen[iSub]; ++i) {
      int j;
      for (j=0; j<dim1; ++j) 
	y.v[locOffset+i][j] = this->v[locOffset+i][j];
      for (j=0; j<dim2; ++j) 
	z.v[locOffset+i][j] = this->v[locOffset+i][dim1 + j];
    }
      
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<int dim1, int dim2>
inline
void
DistSVec<Scalar,dim>::merge(DistSVec<Scalar,dim1> &y, DistSVec<Scalar,dim2> &z)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub) {
 
    int locOffset = distInfo.subOffset[iSub];
    
    for (int i=0; i<distInfo.subLen[iSub]; ++i) {
      int j;
      for (j=0; j<dim1; ++j) 
	this->v[locOffset+i][j] = y.v[locOffset+i][j];
      for (j=0; j<dim2; ++j) 
	this->v[locOffset+i][dim1 + j] = z.v[locOffset+i][j];
    }
      
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<int dim1>
inline
void
DistSVec<Scalar,dim>::strip(DistSVec<Scalar,dim1> &y)
{

  int offset = 0;
  if (dim > dim1)
    if (dim1 < 5)
      offset = 5;

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

    int locOffset = distInfo.subOffset[iSub];

    for (int i=0; i<distInfo.subLen[iSub]; ++i) {
      int j;
      for (j=0; j<dim1; ++j)  {
        y.v[locOffset+i][j] = this->v[locOffset+i][offset+j];
      }
    }

  }
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<int dim1>
inline
void
DistSVec<Scalar,dim>::pad(DistSVec<Scalar,dim1> &y)
{

  int offset = 0;
  if (dim < dim1)
    if (dim < 5)
      offset = 5;

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

    int locOffset = distInfo.subOffset[iSub];

    for (int i=0; i<distInfo.subLen[iSub]; ++i) {
      int j;
      for (j=0; j<dim; ++j)  {
        y.v[locOffset+i][j+offset] = this->v[locOffset+i][j];
      }
    }

  }
}

//------------------------------------------------------------------------------
    
template<class Scalar, int dim>
inline
void DistSVec<Scalar,dim>::sum(Scalar sumres[dim]) const
{
  
  int iSub, idim;

  Scalar sum = 0;

  for(idim=0; idim<dim; idim++){

    sum = 0;

#ifndef MPI_OMP_REDUCTION
    Scalar *allsum = reinterpret_cast<Scalar *>(alloca(sizeof(Scalar) * distInfo.numGlobSub));

    for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) allsum[iSub] = 0;
#endif

    if (distInfo.masterFlag) {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: sum)
#else
#pragma omp parallel for
#endif
      for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

        int locOffset = distInfo.subOffset[iSub];
        int locLen = distInfo.subLen[iSub];

        Scalar locsum = 0;

        for (int i = 0; i < locLen; ++i)
	  if (distInfo.masterFlag[locOffset+i])
	    locsum += this->v[locOffset+i][idim];

#ifdef MPI_OMP_REDUCTION
        sum += locsum;
#else
        allsum[distInfo.locSubToGlobSub[iSub]] = locsum;
#endif

      }

    } 
    else {

#ifdef MPI_OMP_REDUCTION
#pragma omp parallel for reduction(+: sum)
#else
#pragma omp parallel for
#endif
      for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

        int locOffset = distInfo.subOffset[iSub];
        int locLen = distInfo.subLen[iSub];

        Scalar locsum = 0;

        for (int i = 0; i < locLen; ++i)
	  locsum += this->v[locOffset+i][idim];

#ifdef MPI_OMP_REDUCTION
        sum += locsum;
#else
        allsum[distInfo.locSubToGlobSub[iSub]] = locsum;
#endif

      }

    }

#ifdef MPI_OMP_REDUCTION
    distInfo.com->globalSum(1, &sum);
#else
    distInfo.com->globalSum(distInfo.numGlobSub, allsum);

    sum = 0;
    for (iSub=0; iSub<distInfo.numGlobSub; ++iSub) sum += allsum[iSub];
#endif

    sumres[idim] = sum;

  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
void DistSVec<Scalar,dim>::min(Scalar vmin[dim]) const
{

  int iSub;

  Scalar *allmin = reinterpret_cast<Scalar *>(alloca(sizeof(Scalar) * distInfo.numLocSub * dim));

#pragma omp parallel for
  for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {
    subVec[iSub]->min(vmin);
    for (int idim=0; idim<dim; idim++)
      allmin[iSub*dim+idim] = vmin[idim];
  }

  for (int idim=0; idim<dim; idim++)
    vmin[idim] = allmin[idim];
  for (iSub = 1; iSub < distInfo.numLocSub; ++iSub)
    for (int idim=0; idim<dim; idim++)
      vmin[idim] = vmin[idim] < allmin[iSub*dim+idim] ? vmin[idim] : allmin[iSub*dim+idim];

  distInfo.com->globalMin(dim, vmin);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
void DistSVec<Scalar,dim>::max(Scalar vmax[dim]) const
{

  int iSub;

  Scalar *allmax = reinterpret_cast<Scalar *>(alloca(sizeof(Scalar) * distInfo.numLocSub * dim));

#pragma omp parallel for
  for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {
    subVec[iSub]->max(vmax);
    for (int idim=0; idim<dim; idim++)
      allmax[iSub*dim+idim] = vmax[idim];
  }

  for (int idim=0; idim<dim; idim++)
    vmax[idim] = allmax[idim];
  for (iSub = 1; iSub < distInfo.numLocSub; ++iSub)
    for (int idim=0; idim<dim; idim++)
      vmax[idim] = vmax[idim] > allmax[iSub*dim+idim] ? vmax[idim] : allmax[iSub*dim+idim];

  distInfo.com->globalMax(dim, vmax);

}

//------------------------------------------------------------------------------
#endif
