#include <StiffMatrix.h>
#include <MemoryPool.h>
#include <BCApplier.h>

//------------------------------------------------------------------------------

template<class Scalar, int dim>
StiffMat<Scalar,dim>::StiffMat(Domain *domain, int **nodeType, MemoryPool *mp, BCApplier* bcs) 
  : DistMat<Scalar,dim>(domain)
{

#ifdef _OPENMP
  this->numLocSub = DistMat<Scalar,dim>::numLocSub; //BUG omp
#endif

  ndType = nodeType;

  BCs = bcs; // HB

  A = new SparseMat<Scalar,dim>*[this->numLocSub];

  double size = 0.0;

#pragma omp parallel for reduction (+: size)
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) {

    A[iSub] = this->subDomain[iSub]->template createMaskJacobian<Scalar,dim>(0, mp);

    size += double(A[iSub]->numNonZeroBlocks()*dim*dim*sizeof(Scalar)) / (1024.*1024.);

  }

  this->com->globalSum(1, &size);

  this->com->printf(2, "Memory for stiffness matrix (dim=%d): %3.2f MB\n", dim, size);

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
StiffMat<Scalar,dim>::~StiffMat()
{

  if (A) {
#pragma omp parallel for
    for (int iSub = 0; iSub < this->numLocSub; ++iSub)
      if (A[iSub]) delete A[iSub];

    delete [] A;
  }

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
DistMat<Scalar,dim> &StiffMat<Scalar,dim>::operator= (const Scalar x)
{

#pragma omp parallel for
  for (int iSub = 0; iSub < this->numLocSub; ++iSub) 
    *A[iSub] = x;

  return *this;

}

//------------------------------------------------------------------------------

template<class Scalar, int dim>
void StiffMat<Scalar,dim>::apply(DistSVec<double,dim> &x, DistSVec<double,dim> &Ax)
{
#ifdef HB_MVP_DEBUG
  fprintf(stderr," -> in StiffMat::apply, CPU %d, GET IN\n",this->com->cpuNum());
#endif
  int iSub;

  // PJSA (moved from preconditioner, see KspPrec.C)
  DistSVec<double,dim> x_copy(x);
  if(BCs) BCs->applyPD(x_copy);
  
#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub) {
    int* subNdType = (ndType) ? ndType[iSub] : 0;
    A[iSub]->apply(x_copy(iSub), Ax(iSub), subNdType);
    this->subDomain[iSub]->sndData(*this->vecPat, Ax.subData(iSub));
  }

  this->vecPat->exchange();

#pragma omp parallel for
  for (iSub = 0; iSub < this->numLocSub; ++iSub)
    this->subDomain[iSub]->addRcvData(*this->vecPat, Ax.subData(iSub));

#ifdef HB_MVP_DEBUG
  fprintf(stderr," -> in StiffMat::apply, CPU %d, HERE 0\n",this->com->cpuNum());
#endif
  if(BCs) BCs->applyPDt(Ax); //HB: apply (transpose) of projector onto sliding BCs
#ifdef HB_MVP_DEBUG
  fprintf(stderr," -> in StiffMat::apply, CPU %d, GET OUT\n",this->com->cpuNum());
#endif
}

//------------------------------------------------------------------------------
