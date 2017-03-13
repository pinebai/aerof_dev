#include <ParallelRom.C>
#include <ParallelRomExtension.cpp>
#include <RefVector.h>

#define INSTANTIATION_HELPER_1(dim)\
template \
ParallelRom<dim>::ParallelRom(Domain & _domain, Communicator *_com, const DistInfo& dI);\
\
template \
ParallelRom<dim>::~ParallelRom();\
\
template \
void ParallelRom<dim>::scalapackCpuDecomp(const int nCol);\
\
template \
void ParallelRom<dim>::transferDataBackLS (double *subMatB, int n, double \
    **lsSol, int nRhs, int subMatLLD, bool lsCoeffAllCPU);\
\
template \
void ParallelRom<dim>::setTransfer();\
\
template \
void ParallelRom<dim>::parallelLSMultiRHSClean();\
\
template \
ParallelRomExtension<dim>::ParallelRomExtension(Domain &, Communicator *, const DistInfo&);\
\
template \
ParallelRomExtension<dim>:: ~ParallelRomExtension();\
\
template \
void ParallelRomExtension<dim>::summonSlaves(double *&mem, VecSet<DistSVec<double, dim> > &X, const int M, const int N);\
\
template \
void  ParallelRomExtension<dim>::summonZombies(double *&mem, VecSet<DistSVec<double, dim> > &X, const int M, const int N);\
\
template \
void  ParallelRomExtension<dim>::transpose(double* &buff1, double* &buff2, int nrow, int ncol);\
\
template \
std::vector<int>  ParallelRomExtension<dim>::countMasters(const DistInfo &distinfo);\

#define INSTANTIATION_HELPER_2(dim, VecContainer)\
template \
void ParallelRom<dim>::transferData(VecContainer &snaps, double* subMat, int nSnaps);\
\
template \
void ParallelRom<dim>::transferDataBack(double *U, VecContainer &Utrue , int nSnaps);\
\

#define INSTANTIATION_HELPER_3(dim, VecContainer1, VecContainer2)\
template \
void ParallelRom<dim>::parallelSVD(VecContainer1 &snaps, VecContainer2 &Utrue,\
    double *S, FullM &Vtrue, int nSnaps, bool computeV);\
\
template \
void ParallelRom<dim>::parallelLSMultiRHSInit(const VecContainer1 &A, const VecContainer2 &B, int _nA);\
\
template \
void ParallelRom<dim>::parallelLSMultiRHS(const VecContainer1 &A,\
    const VecContainer2 &B, int n, int nRhs, double **lsSol, bool lsCoeffAllCPU);\

#define INSTANTIATION_HELPER_4(dim, DataType, MatType)\
template \
void ParallelRomExtension<dim>::freeSlaves(DataType *&mem, const MatType &X, const int M, const int N);

#define INSTANTIATION_HELPER_5(dim, MatType1, MatType2)\
template \
void ParallelRomExtension<dim>::parallelALS(const MatType1 &X, const MatType2 &M, MatType1 &UT, int maxIts);

INSTANTIATION_HELPER_1(5);
INSTANTIATION_HELPER_1(6);
INSTANTIATION_HELPER_1(7);

typedef VecSet< DistSVec<double,5> > VecSet5;
typedef VecSet< DistSVec<double,6> > VecSet6;
typedef VecSet< DistSVec<double,7> > VecSet7;

typedef VecSet< DistSVec<char,5> > CharMat5;
typedef VecSet< DistSVec<char,6> > CharMat6;
typedef VecSet< DistSVec<char,7> > CharMat7;

typedef RefVec<DistSVec<double, 5> > RefVec5;
typedef RefVec<DistSVec<double, 6> > RefVec6;
typedef RefVec<DistSVec<double, 7> > RefVec7;

INSTANTIATION_HELPER_2(5,VecSet5);
INSTANTIATION_HELPER_2(6,VecSet6);
INSTANTIATION_HELPER_2(7,VecSet7);

INSTANTIATION_HELPER_3(5,VecSet5,VecSet5);
INSTANTIATION_HELPER_3(6,VecSet6,VecSet6);
INSTANTIATION_HELPER_3(7,VecSet7,VecSet7);

INSTANTIATION_HELPER_3(5,VecSet5,RefVec5);
INSTANTIATION_HELPER_3(6,VecSet6,RefVec6);
INSTANTIATION_HELPER_3(7,VecSet7,RefVec7);

INSTANTIATION_HELPER_4(5,double,VecSet5);
INSTANTIATION_HELPER_4(6,double,VecSet6);
INSTANTIATION_HELPER_4(7,double,VecSet7);

INSTANTIATION_HELPER_4(5,char,CharMat5);
INSTANTIATION_HELPER_4(6,char,CharMat6);
INSTANTIATION_HELPER_4(7,char,CharMat7);

INSTANTIATION_HELPER_5(5,VecSet5,CharMat5);
INSTANTIATION_HELPER_5(6,VecSet6,CharMat6);
INSTANTIATION_HELPER_5(7,VecSet7,CharMat7);

#undef INSTANTIATION_HELPER_1
#undef INSTANTIATION_HELPER_2
#undef INSTANTIATION_HELPER_3
#undef INSTANTIATION_HELPER_4
