#ifndef _DIST_EMBEDDED_VECTOR_H_
#define _DIST_EMBEDDED_VECTOR_H_

#include <DistVector.h>
#include <set>
#include <cmath>
//------------------------------------------------------------------------------

template<class T, class Scalar>
class EmbeddedExpr {
  
  const T &x;
 
  public:

    const std::set<int>& stencil(int iSub) const { return x.stencil(iSub); }

    EmbeddedExpr(const T& _x) : x(_x) { }

    Scalar real(int iSub,int i,int k) const { return x.real(iSub,i,k); }
    Scalar ghost(int iSub,int i,int k) const { return x.ghost(iSub,i,k); }

    bool hasHHBoundaryTerm() const { return x.hasHHBoundaryTerm(); }

    const DistInfo& hhInfo() const { return x.hhInfo(); }

    Scalar hh(int iSub, int i) const { return x.hh(iSub,i); }
};

template<class T1,class T2, class Scalar>
class EmbeddedSumExpr : public EmbeddedExpr< EmbeddedSumExpr<T1,T2,Scalar> , Scalar> {
  
  const T1 &x;
  const T2 &y;
 
  public:

    const std::set<int>& stencil(int iSub) const { return x.stencil(iSub); }
    
    EmbeddedSumExpr(const T1& _x, const T2& _y) : x(_x), y(_y), 
      EmbeddedExpr< EmbeddedSumExpr<T1,T2,Scalar> , Scalar>(*this) { }

    Scalar real(int iSub,int i,int k) const { return x.real(iSub,i,k) + y.real(iSub,i,k) ; }
    Scalar ghost(int iSub,int i,int k) const { return x.ghost(iSub,i,k) + y.ghost(iSub,i,k); }

    bool hasHHBoundaryTerm() const { return x.hasHHBoundaryTerm(); }
    const DistInfo& hhInfo() const { return x.hhInfo(); }

    Scalar hh(int iSub, int i) const { return x.hh(iSub,i)+y.hh(iSub,i); }
};

template<class T1,class T2, class Scalar>
class EmbeddedDiffExpr : public EmbeddedExpr< EmbeddedDiffExpr<T1,T2,Scalar> , Scalar> {
  
  const T1 &x;
  const T2 &y;
 
  public:
    
    const std::set<int>& stencil(int iSub) const { return x.stencil(iSub); }

    EmbeddedDiffExpr(const T1& _x, const T2& _y) : x(_x), y(_y),
      EmbeddedExpr< EmbeddedDiffExpr<T1,T2,Scalar> , Scalar>(*this) { }

    Scalar real(int iSub,int i,int k) const { return x.real(iSub,i,k) - y.real(iSub,i,k) ; }
    Scalar ghost(int iSub,int i,int k) const { return x.ghost(iSub,i,k) - y.ghost(iSub,i,k); }

    bool hasHHBoundaryTerm() const { return x.hasHHBoundaryTerm(); }
    const DistInfo& hhInfo() const { return x.hhInfo(); }

    Scalar hh(int iSub, int i) const { return x.hh(iSub,i)-y.hh(iSub,i); }
};

template<class T, class Scalar>
class EmbeddedScaleExpr : public EmbeddedExpr< EmbeddedScaleExpr<T,Scalar> , Scalar> {
  
  const T &x;
  
  Scalar a; 

  public:

  const std::set<int>& stencil(int iSub) const { return x.stencil(iSub); }
    
    EmbeddedScaleExpr(const T& _x, Scalar _a) : x(_x), a(_a),
      EmbeddedExpr< EmbeddedScaleExpr<T,Scalar> , Scalar>(*this) { }

    Scalar real(int iSub,int i,int k) const { return a*x.real(iSub,i,k); }
    Scalar ghost(int iSub,int i,int k) const { return a*x.ghost(iSub,i,k); }

    bool hasHHBoundaryTerm() const { return x.hasHHBoundaryTerm(); }
    const DistInfo& hhInfo() const { return x.hhInfo(); }

    Scalar hh(int iSub, int i) const { return a*x.hh(iSub,i); }
};

template <class T1, class T2, class Scalar>
EmbeddedSumExpr<EmbeddedExpr<T1,Scalar>, EmbeddedExpr<T2,Scalar>,Scalar> operator + (const EmbeddedExpr<T1,Scalar>& x,
                                          const EmbeddedExpr<T2,Scalar>& y) {
 
  return EmbeddedSumExpr<EmbeddedExpr<T1,Scalar>, EmbeddedExpr<T2,Scalar>,Scalar>(x,y);
}

template <class T1, class T2, class Scalar>
EmbeddedDiffExpr<EmbeddedExpr<T1,Scalar>,EmbeddedExpr<T2,Scalar>,Scalar> operator - (const EmbeddedExpr<T1,Scalar>& x,
                                           const EmbeddedExpr<T2,Scalar>& y) {
 
  return EmbeddedDiffExpr<EmbeddedExpr<T1,Scalar>,EmbeddedExpr<T2,Scalar>,Scalar>(x,y);
}

template <class T, class Scalar>
EmbeddedScaleExpr<EmbeddedExpr<T,Scalar>,Scalar> operator * (const EmbeddedExpr<T,Scalar>& x,
                                        const Scalar& a) {
 
  return EmbeddedScaleExpr<EmbeddedExpr<T,Scalar>,Scalar>(x,a);
}

template <class T, class Scalar>
EmbeddedScaleExpr<EmbeddedExpr<T,Scalar>,Scalar> operator * (const Scalar& a,
                                        const EmbeddedExpr<T,Scalar>& x) {
 
  return EmbeddedScaleExpr<EmbeddedExpr<T,Scalar>,Scalar>(x,a);
}

//------------------------------------------------------------------------------
template<class Scalar, int dim>
class DistEmbeddedVec : public EmbeddedExpr< DistEmbeddedVec<Scalar,dim>, Scalar> {

public:

  typedef DistInfo InfoType;

private:

  DistSVec<Scalar,dim> realVec,ghostVec;
  std::set<int>* ghostNodes;

  DistVec<Scalar>* hhBoundaryTerm;

public:

  DistEmbeddedVec(const DistInfo &);
  DistEmbeddedVec(const DistEmbeddedVec<Scalar,dim> &);
  //DistEmbeddedVec(const DistInfo &, Scalar (*)[dim]);
  ~DistEmbeddedVec();

  void addHHBoundaryTerm(const DistInfo & I) { hhBoundaryTerm = new DistVec<Scalar>(I); }

  DistEmbeddedVec<Scalar,dim> &operator=(const Scalar);
  DistEmbeddedVec<Scalar,dim> &operator*=(const Scalar);
  DistEmbeddedVec<Scalar,dim> &operator/=(const DistVec<Scalar> &);
  DistEmbeddedVec<Scalar,dim> &operator*=(const DistVec<Scalar> &);
  DistEmbeddedVec<Scalar,dim> operator* (const Scalar);
  DistEmbeddedVec<Scalar,dim> &operator+=(const Scalar);


  DistEmbeddedVec<Scalar,dim> &operator=(const DistEmbeddedVec<Scalar,dim> &);
  DistEmbeddedVec<Scalar,dim> &operator+=(const DistEmbeddedVec<Scalar,dim> &);
  DistEmbeddedVec<Scalar,dim> &operator-=(const DistEmbeddedVec<Scalar,dim> &);
  DistEmbeddedVec<Scalar,dim> &operator=(const DistVec<Scalar> &);

  DistEmbeddedVec<Scalar,dim> &conjugate();
  DistEmbeddedVec<Scalar,dim> &randVec();
  DistEmbeddedVec<Scalar,dim> &getReal(const DistEmbeddedVec<bcomp,dim> &);
  DistEmbeddedVec<Scalar,dim> &getImag(const DistEmbeddedVec<bcomp,dim> &);
  DistEmbeddedVec<Scalar,dim> &getAbs(const DistEmbeddedVec<bcomp,dim> &);

  Scalar operator*(const DistEmbeddedVec<Scalar,dim> &);

  bool hasHHBoundaryTerm() const { return hhBoundaryTerm; }

  void setHH(DistVec<Scalar>& hh);

  void min(Scalar vmin[dim]) const;
  
  void max(Scalar vmax[dim]) const;

  template<class T>
  DistEmbeddedVec<Scalar,dim> &operator=(const EmbeddedExpr<T, Scalar> &);

  template<class T>
  DistEmbeddedVec<Scalar,dim> &operator+=(const EmbeddedExpr<T, Scalar> &);

  template<class T>
  DistEmbeddedVec<Scalar,dim> &operator-=(const EmbeddedExpr<T, Scalar> &);

  template<class T>
  Scalar operator*(const EmbeddedExpr<T, Scalar> &);
  
  Scalar real(int iSub,int i,int k) const { return realVec.subData(iSub)[i][k]; }
  Scalar ghost(int iSub,int i,int k) const { return ghostVec.subData(iSub)[i][k]; }
  const std::set<int>& stencil(int iSub) const { return ghostNodes[iSub]; }

  Scalar hh(int iSub, int i) const { return hhBoundaryTerm->subData(iSub)[i]; }

  const DistInfo& hhInfo() const { return hhBoundaryTerm->info(); }

  void setGhostStencil(DistVec<GhostPoint<dim>*>& gp);
  void setGhost(DistVec<GhostPoint<dim>*>& gp,VarFcn* vf);
  void getGhost(DistVec<GhostPoint<dim>*>& gp,VarFcn* vf);

  DistSVec<Scalar,dim>& real() { return realVec; }
  DistSVec<Scalar,dim>& ghost() { return ghostVec; }

  DistVec<Scalar>& hh() { return *hhBoundaryTerm; }

  template<int dim1, int dim2>
  void split(DistEmbeddedVec<Scalar,dim1> &, DistEmbeddedVec<Scalar,dim2> &);

  template<int dim1, int dim2>
  void merge(DistEmbeddedVec<Scalar,dim1> &, DistEmbeddedVec<Scalar,dim2> &);

  template<int dim1>
  void strip(DistEmbeddedVec<Scalar,dim1> &);

  template<int dim1>
  void pad(DistEmbeddedVec<Scalar,dim1> &);

  void createSubVec();

  void set(const Scalar *);
  
  void restrict();

  void average();

  int nonOverlapSize() const;

  DistEmbeddedVec<Scalar,dim> *alias() const;

  const DistInfo &info() const { return real().info(); }
 
  double sizeMB() { return realVec.sizeMB()*2; }

  //int numLocSub() const { return distInfo.numLocSub; }

  //int numGlobSub() const { return distInfo.numGlobSub; }

  //int subSize(int i) const { return distInfo.subLen[i]; }

  //Scalar (*subData(int i) const)[dim] { return this->v+real().inf().subOffset[i]; }

  //bool *getMasterFlag(int i) const { return real().info().getMasterFlag(i); }

  double norm();

  void sum(Scalar sumres[dim]) const;
};

//------------------------------------------------------------------------------

template<class Scalar, int dim>
DistEmbeddedVec<Scalar,dim>::DistEmbeddedVec(const DistInfo &dI) : 
  realVec(dI), ghostVec(dI), EmbeddedExpr< DistEmbeddedVec<Scalar,dim>,Scalar>(*this)
{
  ghostNodes = new std::set<int>[ dI.numLocSub ];
  hhBoundaryTerm = NULL;
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
DistEmbeddedVec<Scalar,dim>::DistEmbeddedVec(const DistEmbeddedVec<Scalar,dim> &y) :
  realVec(y.realVec), ghostVec(y.ghostVec), EmbeddedExpr< DistEmbeddedVec<Scalar,dim>,Scalar>(*this)
  
{
  int nls = y.realVec.info().numLocSub ; 
  ghostNodes = new std::set<int>[ nls ];
#pragma omp parallel for
  for (int i = 0; i < nls; ++i) {
    ghostNodes[i] = y.ghostNodes[i];
  }
  
  if (y.hhBoundaryTerm)
    hhBoundaryTerm = new DistVec<Scalar>(*y.hhBoundaryTerm);
}

template<class Scalar, int dim>
void DistEmbeddedVec<Scalar,dim>::setHH(DistVec<Scalar>& hh) {

  if (!hhBoundaryTerm)
    hhBoundaryTerm = new DistVec<Scalar>(hh);
  else
    *hhBoundaryTerm = hh;
}

//------------------------------------------------------------------------------

/*template<class Scalar, int dim>
DistEmbeddedVec<Scalar,dim>::DistEmbeddedVec(const DistInfo &dI, Scalar (*vv)[dim]) :
  EmbeddedVec<Scalar,dim>(dI.totLen, vv), distInfo(dI), subVec(0)
{

  createSubVec();

}
*/
//------------------------------------------------------------------------------

template<class Scalar, int dim>
DistEmbeddedVec<Scalar,dim>::~DistEmbeddedVec() 
{ 
  delete [] ghostNodes;
  if (hhBoundaryTerm)
    delete hhBoundaryTerm;
}

//------------------------------------------------------------------------------
/*template<class Scalar, int dim>
inline
void
DistEmbeddedVec<Scalar,dim>::set(const Scalar *y)
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
*/
//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
double
DistEmbeddedVec<Scalar,dim>::norm()
{
  int iSub,k;
  double res = 0;

  const DistInfo& distInfo = realVec.info();
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
            locres += sqNorm(this->realVec.v[locOffset+i][j]);
    
      for (std::set<int>::iterator itr = ghostNodes[iSub].begin(); itr != ghostNodes[iSub].end(); ++itr) {
        k = *itr;
        if (distInfo.masterFlag[locOffset+k])
          for (int j = 0; j < dim; ++j)
            locres += sqNorm(this->ghostVec.v[locOffset+k][j]);
      }

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
          locres += sqNorm(this->realVec.v[locOffset+i][j]);

      for (std::set<int>::iterator itr = ghostNodes[iSub].begin(); itr != ghostNodes[iSub].end(); ++itr) {
        k = *itr;
        if (distInfo.masterFlag[locOffset+k])
          for (int j = 0; j < dim; ++j)
            locres += sqNorm(this->ghostVec.v[locOffset+k][j]);
      }

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

  if (hhBoundaryTerm)
    res += (*hhBoundaryTerm)*(*hhBoundaryTerm);

  return sqrt(res);
}

//------------------------------------------------------------------------------
/*
template<class Scalar, int dim>
inline
DistEmbeddedVec<Scalar,dim> *
DistEmbeddedVec<Scalar,dim>::alias() const
{

  return new DistEmbeddedVec<Scalar,dim>(distInfo, this->v); 

}
*/
//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistEmbeddedVec<Scalar,dim>
DistEmbeddedVec<Scalar,dim>::operator*(const Scalar y)
{

  DistEmbeddedVec<Scalar, dim> result(*this);
  result *= y;

  return result;

}

//------------------------------------------------------------------------------
/*
template<class Scalar, int dim>
inline
DistEmbeddedVec<Scalar,dim> &
DistEmbeddedVec<Scalar,dim>::operator+=(const Scalar y)
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
*/
//------------------------------------------------------------------------------
/*
template<class Scalar, int dim>
inline
DistEmbeddedVec<Scalar,dim> &
DistEmbeddedVec<Scalar,dim>::getImag(const DistEmbeddedVec<bcomp,dim> &y)
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
*/
//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistEmbeddedVec<Scalar,dim> &
DistEmbeddedVec<Scalar,dim>::operator=(const Scalar y)
{
  
  ghostVec = realVec = y;
  if (hhBoundaryTerm)
    *hhBoundaryTerm = y;
  return *this;
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistEmbeddedVec<Scalar,dim> &
DistEmbeddedVec<Scalar,dim>::operator*=(const Scalar y)
{

  const DistInfo& distInfo = realVec.info();
  realVec *= y;

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];
    for (std::set<int>::iterator itr = ghostNodes[iSub].begin(); 
         itr != ghostNodes[iSub].end(); ++itr) {
      ghostVec.v[locOffset+*itr]*=y;
    }
  }

  if (hhBoundaryTerm)
    *hhBoundaryTerm *= y;
  

  return *this;

}

//------------------------------------------------------------------------------

/*template<class Scalar, int dim>
inline
DistEmbeddedVec<Scalar,dim> &
DistEmbeddedVec<Scalar,dim>::operator/=(const DistVec<Scalar> &y)
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
*/
//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistEmbeddedVec<Scalar,dim> &
DistEmbeddedVec<Scalar,dim>::operator=(const DistEmbeddedVec<Scalar,dim> &y)
{
  ghostVec = y.ghostVec;
  realVec = y.realVec;
  
  int nls = y.realVec.info().numLocSub ; 
  delete [] ghostNodes; // fixes memory leak
  ghostNodes = new std::set<int>[ nls ];
#pragma omp parallel for
  for (int i = 0; i < nls; ++i) {
    ghostNodes[i] = y.ghostNodes[i];
  }

  if (hhBoundaryTerm)
    delete hhBoundaryTerm;

  hhBoundaryTerm = NULL;
  if (y.hhBoundaryTerm)
    hhBoundaryTerm = new DistVec<double>(*y.hhBoundaryTerm);

  return *this;

}

//------------------------------------------------------------------------------

/*template<class Scalar, int dim>
inline
DistEmbeddedVec<Scalar,dim> &
DistEmbeddedVec<Scalar,dim>::operator=(const DistVec<Scalar> &y)
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
*/

//------------------------------------------------------------------------------
/*
template<class Scalar, int dim>
inline
DistEmbeddedVec<Scalar,dim> &
DistEmbeddedVec<Scalar,dim>::operator*=(const DistVec<Scalar> &y)
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
*/
//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistEmbeddedVec<Scalar,dim> &
DistEmbeddedVec<Scalar,dim>::operator+=(const DistEmbeddedVec<Scalar,dim> &y)
{

  const DistInfo& distInfo = realVec.info();
  const Scalar *yy = reinterpret_cast<Scalar *>(y.realVec.v);
  Scalar *vv = reinterpret_cast<Scalar *>(this->realVec.v);

  const Scalar *g_yy = reinterpret_cast<Scalar *>(y.ghostVec.v);
  Scalar *g_vv = reinterpret_cast<Scalar *>(this->ghostVec.v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) { 
      vv[locOffset+i] += yy[locOffset+i];
    }

    for (std::set<int>::iterator itr = ghostNodes[iSub].begin(); 
         itr != ghostNodes[iSub].end(); ++itr) {
      g_vv[locOffset+*itr] += g_yy[locOffset+*itr];
    }
  }

  if (hhBoundaryTerm)
    *hhBoundaryTerm += *y.hhBoundaryTerm;

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
inline
DistEmbeddedVec<Scalar,dim> &
DistEmbeddedVec<Scalar,dim>::operator-=(const DistEmbeddedVec<Scalar,dim> &y)
{

  const DistInfo& distInfo = realVec.info();
  const Scalar *yy = reinterpret_cast<Scalar *>(y.realVec.v);
  Scalar *vv = reinterpret_cast<Scalar *>(this->realVec.v);

  const Scalar *g_yy = reinterpret_cast<Scalar *>(y.ghostVec.v);
  Scalar *g_vv = reinterpret_cast<Scalar *>(this->ghostVec.v);

#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocThreads; ++iSub) {

    int locOffset = dim * distInfo.subOffsetReg[iSub];
    int locLen = dim * distInfo.subLenReg[iSub];

    for (int i = 0; i < locLen; ++i) { 
      vv[locOffset+i] += yy[locOffset+i];
    }

    for (std::set<int>::iterator itr = ghostNodes[iSub].begin(); 
         itr != ghostNodes[iSub].end(); ++itr) {
      g_vv[locOffset+*itr] -= g_yy[locOffset+*itr];
    }
  }

  if (hhBoundaryTerm)
    *hhBoundaryTerm -= *y.hhBoundaryTerm;

  return *this;
}

//------------------------------------------------------------------------------
    
template<class Scalar, int dim>
inline
Scalar 
DistEmbeddedVec<Scalar,dim>::operator*(const DistEmbeddedVec<Scalar,dim> &x)
{

  int iSub,k;
  const DistInfo& distInfo = realVec.info();

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
	    locres += DotTerm(this->realVec.v[locOffset+i][j],x.realVec.v[locOffset+i][j]);
      
      for (std::set<int>::iterator itr = ghostNodes[iSub].begin(); itr != ghostNodes[iSub].end(); ++itr) {
        k = *itr;
        if (distInfo.masterFlag[locOffset+k])
          for (int j = 0; j < dim; ++j)
            locres += DotTerm(this->ghostVec.v[locOffset+k][j], x.ghostVec.v[locOffset+k][j]);
      }

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
	  locres += DotTerm(this->realVec.v[locOffset+i][j], x.realVec.v[locOffset+i][j]);
      
      for (std::set<int>::iterator itr = ghostNodes[iSub].begin(); itr != ghostNodes[iSub].end(); ++itr) {
        k = *itr;
        for (int j = 0; j < dim; ++j)
          locres += DotTerm(this->ghostVec.v[locOffset+k][j], x.ghostVec.v[locOffset+k][j]);
      }
      
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

  if (hhBoundaryTerm)
    res += (*hhBoundaryTerm)*(*x.hhBoundaryTerm);

  return res;

}

//------------------------------------------------------------------------------
/*
template<class Scalar, int dim>
inline
DistEmbeddedVec<Scalar,dim> &
DistEmbeddedVec<Scalar,dim>::conjugate()
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
DistEmbeddedVec<Scalar,dim> &
DistEmbeddedVec<Scalar,dim>::getReal(const DistEmbeddedVec<bcomp,dim> &y)
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
DistEmbeddedVec<Scalar,dim> &
DistEmbeddedVec<Scalar,dim>::getAbs(const DistEmbeddedVec<bcomp,dim> &y)
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
DistEmbeddedVec<Scalar,dim> &
DistEmbeddedVec<Scalar,dim>::randVec()
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
*/
//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<class T>
inline
DistEmbeddedVec<Scalar,dim> &
DistEmbeddedVec<Scalar,dim>::operator=(const EmbeddedExpr<T, Scalar> &expr)
{

  const DistInfo& distInfo = realVec.info();
#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

    int locLen = distInfo.subLen[iSub];
    for (int i = 0; i < locLen; ++i) {
      for (int j = 0; j < dim; ++j)
        realVec.subData(iSub)[i][j] = expr.real(iSub,i,j);
    }

    ghostNodes[iSub] = expr.stencil(iSub);
    for (std::set<int>::iterator itr = ghostNodes[iSub].begin(); itr != ghostNodes[iSub].end(); ++itr) {

      for (int j = 0; j < dim; ++j)
        ghostVec.subData(iSub)[*itr][j] = expr.ghost(iSub,*itr,j);
      
    }
  }

  if (expr.hasHHBoundaryTerm()) {

    if (!hhBoundaryTerm) {
      hhBoundaryTerm = new DistVec<Scalar>(expr.hhInfo());
    }

    const DistInfo& distInfoHH = expr.hhInfo();
    //#pragma omp parallel for
    for (int iSub = 0; iSub < distInfoHH.numLocSub; ++iSub) {

      int locLen = distInfoHH.subLen[iSub];
      for (int i = 0; i < locLen; ++i) {
	Scalar res =  expr.hh(iSub,i);
	hhBoundaryTerm->subData(iSub)[i] = res;
      }
    }
  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<class T>
inline
DistEmbeddedVec<Scalar,dim> &
DistEmbeddedVec<Scalar,dim>::operator+=(const EmbeddedExpr<T, Scalar> &expr)
{

  const DistInfo& distInfo = realVec.info();
#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

    int locLen = distInfo.subLen[iSub];
    for (int i = 0; i < locLen; ++i) {
      for (int j = 0; j < dim; ++j)
        realVec.subData(iSub)[i][j] += expr.real(iSub,i,j);
    }

    ghostNodes[iSub] = expr.stencil(iSub);
    for (std::set<int>::iterator itr = ghostNodes[iSub].begin(); itr != ghostNodes[iSub].end(); ++itr) {

      for (int j = 0; j < dim; ++j)
        ghostVec.subData(iSub)[*itr][j] += expr.ghost(iSub,*itr,j);
      
    }
  }

  if (expr.hasHHBoundaryTerm()) {

    if (!hhBoundaryTerm) {
      std::cout << "Internal error: DistEmbeddedVector.h:" << __LINE__ << std::endl;
      exit(-1);
    }

    const DistInfo& distInfoHH = expr.hhInfo();
#pragma omp parallel for
    for (int iSub = 0; iSub < distInfoHH.numLocSub; ++iSub) {

      int locLen = distInfoHH.subLen[iSub];
      for (int i = 0; i < locLen; ++i) {
	hhBoundaryTerm->subData(iSub)[i] += expr.hh(iSub,i);
      }
    }
  }

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<class T>
inline
DistEmbeddedVec<Scalar,dim> &
DistEmbeddedVec<Scalar,dim>::operator-=(const EmbeddedExpr<T, Scalar> &expr)
{

  const DistInfo& distInfo = realVec.info();
#pragma omp parallel for
  for (int iSub = 0; iSub < distInfo.numLocSub; ++iSub) {

    int locLen = distInfo.subLen[iSub];
    for (int i = 0; i < locLen; ++i) {
      for (int j = 0; j < dim; ++j)
        realVec.subData(iSub)[i][j] -= expr.real(iSub,i,j);
    }

    ghostNodes[iSub] = expr.stencil(iSub);
    for (std::set<int>::iterator itr = ghostNodes[iSub].begin(); itr != ghostNodes[iSub].end(); ++itr) {

      for (int j = 0; j < dim; ++j)
        ghostVec.subData(iSub)[*itr][j] -= expr.ghost(iSub,*itr,j);
      
    }
  }

  if (expr.hasHHBoundaryTerm()) {

    if (!hhBoundaryTerm) {
      std::cout << "Internal error: DistEmbeddedVector.h:" << __LINE__ << std::endl;
      exit(-1);
    }

    const DistInfo& distInfoHH = expr.hhInfo();
#pragma omp parallel for
    for (int iSub = 0; iSub < distInfoHH.numLocSub; ++iSub) {

      int locLen = distInfoHH.subLen[iSub];
      for (int i = 0; i < locLen; ++i) {
	hhBoundaryTerm->subData(iSub)[i] -= expr.hh(iSub,i);
      }
    }
  }

  return *this;


}

//------------------------------------------------------------------------------
    
template<class Scalar, int dim>
template<class T>
inline
Scalar 
DistEmbeddedVec<Scalar,dim>::operator*(const EmbeddedExpr<T, Scalar> &expr)
{

  int iSub,k;

  Scalar res = 0;

  const DistInfo& distInfo = realVec.info();
  const DistInfo& distInfoHH = expr.hhInfo();
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
	    locres += DotTerm(this->realVec.v[locOffset+i][j],expr.real(iSub,i,j));
      
      for (std::set<int>::iterator itr = ghostNodes[iSub].begin(); itr != ghostNodes[iSub].end(); ++itr) {
        k = *itr;
        if (distInfo.masterFlag[locOffset+k])
          for (int j = 0; j < dim; ++j)
            locres += DotTerm(this->ghostVec.v[locOffset+k][j], expr.ghost(iSub,k,j));
      }

      if (expr.hasHHBoundaryTerm()) {

	if (!hhBoundaryTerm) {
	  std::cout << "Internal error: DistEmbeddedVector.h:" << __LINE__ << std::endl;
	  exit(-1);
	}
	locOffset = distInfoHH.subOffset[iSub];
	locLen = distInfoHH.subLen[iSub];
	for (int i = 0; i < locLen; ++i) {
	  if (distInfoHH.masterFlag[locOffset+i])
	    locres += DotTerm(hhBoundaryTerm->v[locOffset+i], expr.hh(iSub,i));
	}
      }
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
	  locres += DotTerm(this->realVec.v[locOffset+i][j], expr.real(iSub,i,j));
      
      for (std::set<int>::iterator itr = ghostNodes[iSub].begin(); itr != ghostNodes[iSub].end(); ++itr) {
        k = *itr;
        for (int j = 0; j < dim; ++j)
          locres += DotTerm(this->ghostVec.v[locOffset+k][j], expr.ghost(iSub,k,j));
      }
      
      if (expr.hasHHBoundaryTerm()) {

	if (!hhBoundaryTerm) {
	  std::cout << "Internal error: DistEmbeddedVector.h:" << __LINE__ << std::endl;
	  exit(-1);
	}
	locOffset = distInfoHH.subOffset[iSub];
	locLen = distInfoHH.subLen[iSub];
	for (int i = 0; i < locLen; ++i) {
	  locres += DotTerm(hhBoundaryTerm->v[locOffset+i], expr.hh(iSub,i));
	}
      }

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
template<int dim1, int dim2>
inline
void
DistEmbeddedVec<Scalar,dim>::split(DistEmbeddedVec<Scalar,dim1> &y, DistEmbeddedVec<Scalar,dim2> &z)
{
  this->realVec.split(y.real(),z.real());
  this->ghostVec.split(y.ghost(),z.ghost());
}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
template<int dim1, int dim2>
inline
void
DistEmbeddedVec<Scalar,dim>::merge(DistEmbeddedVec<Scalar,dim1> &y, DistEmbeddedVec<Scalar,dim2> &z)
{
  this->realVec.merge(y.real(),z.real());
  this->ghostVec.merge(y.ghost(),z.ghost());
}

//------------------------------------------------------------------------------
/*
template<class Scalar, int dim>
inline
int 
DistEmbeddedVec<Scalar,dim>::nonOverlapSize() const
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
*/
//------------------------------------------------------------------------------
/*
template<class Scalar, int dim>
inline
void
DistEmbeddedVec<Scalar,dim>::restrict()
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
DistEmbeddedVec<Scalar,dim>::average()
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
*/
//------------------------------------------------------------------------------
/*
template<class Scalar, int dim>
template<int dim1, int dim2>
inline
void
DistEmbeddedVec<Scalar,dim>::split(DistEmbeddedVec<Scalar,dim1> &y, DistEmbeddedVec<Scalar,dim2> &z)
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
*/
//------------------------------------------------------------------------------
/*
template<class Scalar, int dim>
template<int dim1, int dim2>
inline
void
DistEmbeddedVec<Scalar,dim>::merge(DistEmbeddedVec<Scalar,dim1> &y, DistEmbeddedVec<Scalar,dim2> &z)
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
*/
//------------------------------------------------------------------------------
/*
template<class Scalar, int dim>
template<int dim1>
inline
void
DistEmbeddedVec<Scalar,dim>::strip(DistEmbeddedVec<Scalar,dim1> &y)
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
DistEmbeddedVec<Scalar,dim>::pad(DistEmbeddedVec<Scalar,dim1> &y)
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
*/
//------------------------------------------------------------------------------
/*    
template<class Scalar, int dim>
inline
void DistEmbeddedVec<Scalar,dim>::sum(Scalar sumres[dim]) const
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
*/
//------------------------------------------------------------------------------
/*
template<class Scalar, int dim>
inline
void DistEmbeddedVec<Scalar,dim>::min(Scalar vmin[dim]) const
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
void DistEmbeddedVec<Scalar,dim>::max(Scalar vmax[dim]) const
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
*/
//------------------------------------------------------------------------------
template <class Scalar, int dim>
void DistEmbeddedVec<Scalar,dim>::setGhostStencil(DistVec<GhostPoint<dim>*>& gp) {

  int iSub;
  const DistInfo& distInfo = realVec.info();
#pragma omp parallel for
  for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {
    ghostNodes[iSub].clear();
    for (int i = 0; i < distInfo.subLen[iSub]; ++i) {
      if (gp[i])
        ghostNodes[iSub].insert(i);
    }
  }     
}

template <class Scalar, int dim>
void DistEmbeddedVec<Scalar,dim>::setGhost(DistVec<GhostPoint<dim>*>& gp,VarFcn* vf) {

  int iSub;
  const DistInfo& distInfo = realVec.info();
  GhostPoint<dim>* gpi;
#pragma omp parallel for
  for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {
    ghostNodes[iSub].clear();
    for (int i = 0; i < distInfo.subLen[iSub]; ++i) {
      gpi = gp.subData(iSub)[i];
      if (gpi) {
        ghostNodes[iSub].insert(i);
        vf->primitiveToConservative(gpi->getPrimitiveState(), ghostVec.subData(iSub)[i]);
      } 
    }
  }     
}

template <class Scalar, int dim>
void DistEmbeddedVec<Scalar,dim>::getGhost(DistVec<GhostPoint<dim>*>& gp,VarFcn* vf) {

  int iSub;
  const DistInfo& distInfo = realVec.info();
  GhostPoint<dim>* gpi;
#pragma omp parallel for
  for (iSub = 0; iSub < distInfo.numLocSub; ++iSub) {
    //ghostNodes[iSub].clear();
    for (int i = 0; i < distInfo.subLen[iSub]; ++i) {
      gpi = gp.subData(iSub)[i];
      if (gpi) {
        //ghostNodes[iSub].insert(i);
        vf->conservativeToPrimitive(ghostVec.subData(iSub)[i],gpi->getPrimitiveState());
      } 
    }
  }     
}

#endif
