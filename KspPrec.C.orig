#include <cstdlib>
#include <cstdio>
#include <KspPrec.h>
#include <MemoryPool.h>
#include <BCApplier.h>

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
JacobiPrec<Scalar,dim, Scalar2>::JacobiPrec(typename DiagMat<Scalar,dim>::Type type, 
				   Domain *domain, int **nodeType, BCApplier* bcs) 
  : DistMat<Scalar,dim>(domain)
{ 

#ifdef _OPENMP
  this->numLocSub = DistMat<Scalar,dim>::numLocSub; //BUG omp
#endif

  A = new DiagMat<Scalar,dim>*[this->numLocSub];

  BCs = bcs; //HB

  double size = 0.0;

#pragma omp parallel for reduction (+: size)
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {

    int *ndType = (nodeType) ? nodeType[iSub] : 0;

    A[iSub] = this->subDomain[iSub]->template createMaskDiagonal<Scalar,dim>(type, ndType);

    size += double(A[iSub]->numNonZeroBlocks()*dim*dim*sizeof(Scalar)) / (1024.*1024.);

  }

  this->com->globalSum(1, &size);

  this->com->printf(2, "Memory for Jacobi preconditioner (dim=%d): %3.2f MB\n", dim, size);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
JacobiPrec<Scalar,dim, Scalar2>::~JacobiPrec()
{

  if (A) {
#pragma omp parallel for
    for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
      if (A[iSub]) delete A[iSub];

    delete [] A;
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
DistMat<Scalar,dim> &JacobiPrec<Scalar,dim, Scalar2>::operator= (const Scalar x)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
    *A[iSub] = x;

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void JacobiPrec<Scalar,dim, Scalar2>::setup()
{

#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
    A[iSub]->invert();

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void JacobiPrec<Scalar,dim, Scalar2>::apply(DistSVec<Scalar2,dim> &y, DistSVec<Scalar2,dim> &x) 
{

  int iSub;

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    A[iSub]->apply(y(iSub), x(iSub));

  x.restrict();

  CommPattern<Scalar2> *vPat = this->getCommPat(x);
#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    this->subDomain[iSub]->sndData(*vPat, x.subData(iSub));

  vPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    this->subDomain[iSub]->addRcvData(*vPat, x.subData(iSub));

  //PJSA (moved this to operator, see StiffMatrix.C)  
  //if(BCs) BCs->applyPD(x); //HB: x <- P.x (need to take care of x<-D.x ...)
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
template<class MatScal>
void JacobiPrec<Scalar,dim, Scalar2>::getData(DistMat<MatScal,dim> &B)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub)
    A[iSub]->getData(B(iSub));

}

//------------------------------------------------------------------------------
/*
@ARTICLE{cai-farhat-sarkis-98,
  author = "Cai, X. C. and Farhat, C. and Sarkis, M.",
  title = "A minimum overlap restricted additive {S}chwarz preconditioner and
  applications in 3{D} flow simulations",
  journal = "Contemporary Mathematics",
  year = 1998,
  volume = 218,
  pages = "478--484",
}
*/
template<class Scalar, int dim, class Scalar2>
IluPrec<Scalar,dim, Scalar2>::IluPrec(PcData &pcData, Domain *domain, int **nodeType) : 
  DistMat<Scalar,dim>(domain), jw(domain->getNodeDistInfo()) , domain(domain) {

#ifdef _OPENMP
  this->numLocSub = DistMat<Scalar,dim>::numLocSub; //BUG omp
#endif

  type = pcData.type;

  int renum;
  fill = pcData.fill;

  if (pcData.renumbering == PcData::RCM) renum = 2;
  else if (pcData.renumbering == PcData::NATURAL) renum = 0;

  A = new SparseMat<Scalar,dim>*[this->numLocSub];

  double size = 0.0;

#pragma omp parallel for reduction (+: size)
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {

    int *ndType = (nodeType) ? nodeType[iSub] : 0;

    A[iSub] = this->subDomain[iSub]->template createMaskILU<Scalar,dim>(fill, renum, ndType);
    
    size += double(A[iSub]->numNonZeroBlocks()*dim*dim*sizeof(Scalar)) / (1024.*1024.);

  }

  this->com->globalSum(1, &size);

  this->com->printf(2, "Memory for ILU(%d) preconditioner (dim=%d): %3.2f MB\n", fill, dim, size);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
IluPrec<Scalar,dim, Scalar2>::~IluPrec()
{

  if (A) {
#pragma omp parallel for
    for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
      if (A[iSub]) delete A[iSub];

    delete [] A;
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
DistMat<Scalar,dim> &IluPrec<Scalar,dim, Scalar2>::operator= (const Scalar x)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
    *A[iSub] = x;

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void IluPrec<Scalar,dim, Scalar2>::exportMemory(MemoryPool *mp)
{

  if (!mp) return;

  for (int iSub = 0; iSub < this->numLocSub; ++iSub)
    mp->set(A[iSub]->numNonZeroBlocks() * dim*dim * sizeof(Scalar), 
	    reinterpret_cast<void *>(A[iSub]->data()));

}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void IluPrec<Scalar,dim, Scalar2>::setup()  {

#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub)
    A[iSub]->numericILU(jw(iSub).data());
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void IluPrec<Scalar,dim, Scalar2>::setupTR()  {

#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub)
    A[iSub]->ILUTR();
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void IluPrec<Scalar,dim, Scalar2>::apply(DistSVec<Scalar2,dim> &y, DistSVec<Scalar2,dim> &x) 
{

  int iSub;

  static int cnt = 0;
  static double prec_time = 0.0;
  ++cnt; 
  double t = domain->getTimer()->getTime();
  
  DistSVec<Scalar2, dim> tmp(y);

  if (type == PcData::ASH) tmp.restrict();
  if (type == PcData::AAS) tmp.average();

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    A[iSub]->lusol(tmp(iSub), x(iSub));

  if (type == PcData::RAS) x.restrict();
  if (type == PcData::AAS) x.average();

  CommPattern<Scalar2> *vPat = this->getCommPat(x);
#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    this->subDomain[iSub]->sndData(*vPat, x.subData(iSub));
  vPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    this->subDomain[iSub]->addRcvData(*vPat, x.subData(iSub));

  prec_time += domain->getTimer()->getTime() - t;
//  if (cnt % 50 == 0)
//    domain->getCommunicator()->fprintf(stdout, "Prec time = %e\n",prec_time);
}

//------------------------------------------------------------------------------

template<class Scalar, int dim, class Scalar2>
void IluPrec<Scalar,dim, Scalar2>::applyT(DistSVec<Scalar2,dim> &y, DistSVec<Scalar2,dim> &x)
{

  int iSub;

  DistSVec<Scalar2, dim> tmp(y);

  // switched RAS to ASH

  if (type == PcData::ASH) tmp.restrict();
  if (type == PcData::AAS) tmp.average();

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    A[iSub]->lusolTR(tmp(iSub), x(iSub));

  // switched ASH to RAS
  if (type == PcData::RAS) x.restrict();
  if (type == PcData::AAS) x.average();

  CommPattern<Scalar2> *vPat = this->getCommPat(x);
#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    this->subDomain[iSub]->sndData(*vPat, x.subData(iSub));
  vPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    this->subDomain[iSub]->addRcvData(*vPat, x.subData(iSub));

}

////------------------------------------------------------------------------------
//
//template<class Scalar, int dim, class Scalar2>
//void IluPrec<Scalar,dim, Scalar2>::applyTranspose(DistSVec<Scalar2,dim> &y, DistSVec<Scalar2,dim> &x)
//{
//
//  int iSub;
//
//  DistSVec<Scalar2, dim> tmp(y);
//
//  // switched RAS to ASH
//
//  if (type == PcData::ASH) tmp.restrict();
//  if (type == PcData::AAS) tmp.average();
////  if (type == PcData::RAS) tmp.restrict();
//
//#pragma omp parallel for
//  for (iSub = 0; iSub < this->numLocSub; ++iSub)
//    A[iSub]->lusolTR(tmp(iSub), x(iSub));
//
//  // switched ASH to RAS
//  if (type == PcData::RAS) x.restrict();
//  if (type == PcData::AAS) x.average();
////  if (type == PcData::ASH) x.restrict();
//
//  CommPattern<Scalar2> *vPat = this->getCommPat(x);
//#pragma omp parallel for
//  for (iSub = 0; iSub < this->numLocSub; ++iSub)
//    this->subDomain[iSub]->sndData(*vPat, x.subData(iSub));
//  vPat->exchange();
//
//#pragma omp parallel for
//  for (iSub = 0; iSub < this->numLocSub; ++iSub)
//    this->subDomain[iSub]->addRcvData(*vPat, x.subData(iSub));
//
//}

//------------------------------------------------------------------------------


template<class Scalar, int dim, class Scalar2>
template<class MatScal>
void IluPrec<Scalar,dim, Scalar2>::getData(DistMat<MatScal,dim> &B)
{

  int iSub;
  if (fill > 0)  {
#pragma omp parallel for
    for (iSub = 0; iSub < this->numLocSub; ++iSub)
      A[iSub]->zeroData();
}

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    A[iSub]->getData(B(iSub));

}

//------------------------------------------------------------------------------
template<int dim, class Scalar2>
void
IdentityPrec<dim,Scalar2>::apply(DistSVec<Scalar2,dim> &x, DistSVec<Scalar2,dim> &Ix) 
{ 

  Ix = x; 

  //PJSA (moved this to operator, see StiffMatrix.C)
  //if(BCs) {
  //  BCs->applyPD(Ix); 
  //}

}

