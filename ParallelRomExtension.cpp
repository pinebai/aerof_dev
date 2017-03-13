//
// Created by lei on 1/5/17.
//

#include <ParallelRomExtension.h>
#include <als_lapack.h>

/**
 * default constructor
 */
template <int dim>
ParallelRomExtension<dim>::ParallelRomExtension(Domain &domain, Communicator *comm, const DistInfo &distInfo) :
ParallelRom<dim>(domain, comm, distInfo) {

}

/**
 * default destructor
 */
template <int dim>
ParallelRomExtension<dim>::~ParallelRomExtension() {

}

/**
 * Alternating Least Square Methods for Embedded Framework
 *
 */
template <int dim>
template<class Mat1, class Mat2>
void ParallelRomExtension<dim>::parallelALS(const Mat1 &snap, const Mat2 &mask, Mat1 &UT, int maxIts) {
		//read snapshots and mask
		//assert(numSnapshots == numStateSnapshots);
		// step 1: comute dimension of X, M, UT from given data
		vector<int> masters = countMasters(this->distInfo);
		int nrow = 0;
		for(int i = 0; i < masters.size(); i++){
				nrow += masters[i];
		}
		int ncol = snap.numVectors();
		int k = UT.numVectors();
		this->com->fprintf(stderr, "... (M, N) is [%d, %d], reduced dimension is %d\n", nrow, ncol, k);
		// step 2: convert the data into contagious memory region
		double *X = NULL;
		char *M = NULL;
		freeSlaves<double, Mat1>(X, snap, nrow, ncol); // freeSlaves allocates the memory for X
		freeSlaves<char, Mat2>(M, mask, nrow, ncol); // freeSlaves allocates the memory for M
		// step 3: initialise UT
		this->com->fprintf(stderr, "... initializing U\n");
		//Mat1* basisInit = new Mat1(k, this->domain.getNodeDistInfo());
		//Mat1 basis = *basisInit;
		//ParallelRom<dim> parallelRom(this->domain, this->com, this->domain.getNodeDistInfo());
		double *singularValues = new double[k]; //todo: change definition.
		FullM VInitDummy(k, ncol);
		this->com->fprintf(stderr, "... calling parallelRom.parallelSVD()\n");
		this->parallelSVD(snap, UT, singularValues, VInitDummy, k, true);
		this->com->fprintf(stderr, "... U and V initialized, V dimension is [%d, %d]\n", VInitDummy.numRow(), VInitDummy.numCol());
		double *U = NULL;
		double *UT_ptr = new double[k * nrow * dim];
		this->com->fprintf(stderr, "... allocating space for U, k * nrow * dim = %d\n", k * nrow * dim);
		freeSlaves<double, Mat1>(U, UT, nrow, k); // freeSlaves allocates the memory for U
		transpose(U, UT_ptr, nrow * dim, k);
		this->com->fprintf(stderr, "... U being initialized\n");
		this->com->barrier();
		// step 4: launch external als library
		AlternatingLeastSquare ALS((double *)X, (unsigned char *) M, (double *) UT_ptr, nrow * dim, ncol, k, this->com->getMPIComm(), this->com->cpuNum());
		ALS.run(maxIts);
		// step 5: write result in UT
		this->com->fprintf(stderr, "... ALS finished %d iterations\n", maxIts);
		this->com->barrier();
		transpose(UT_ptr, U, k, nrow * dim);
		this->com->fprintf(stderr, "... transpose done\n");
		this->com->barrier();
		summonSlaves(U, UT, nrow, k);
		//UT = basis;
		// step 6: clean up memory
		this->com->fprintf(stderr, "... cleaning up memories\n");
		if(X) delete[] X;
		if(M) delete[] M;
		if(U) delete[] U;
		if(UT_ptr) delete[] UT_ptr;
		if(singularValues) delete[] singularValues;
		this->com->fprintf(stdout, "... all stuff cleaned\n");
}

/**
 * count the number of master nodes in each subdomain.
 * results stored in this->numMasters[i].
 */
template<int dim>
std::vector<int> ParallelRomExtension<dim>::countMasters(const DistInfo &distinfo){
		vector<int> answer(distinfo.numLocSub, 0);
//#pragma omp parallel for
		for(int iSub = 0; iSub < distinfo.numLocSub; iSub++){
				bool *localMasterFlag = distinfo.getMasterFlag(iSub);
				for(int iNode = 0; iNode < distinfo.subLen[iSub]; iNode++){
						answer[iSub] += localMasterFlag[iNode] ? 1 : 0;
				}
		}
		return answer;
}

/**
 * template version for both data type
 */
template<int dim>
template<class D, class Mat>
void ParallelRomExtension<dim>::freeSlaves(D *&mem, const Mat &X, const int M,
                                           const int N) {
		D* ptr = new D[dim * M * N];
		assert(ptr != NULL);
		for (int j = 0; j < N; j++) {
				int i = 0;
//#pragma omp parallel for
				for (int iSub = 0; iSub <X[j].numLocSub() ; iSub++) {
						bool *localMasterFlag = X[j].getMasterFlag(iSub);
						D (*temp)[dim] = X[j].subData(iSub);
						for (int iNode = 0; iNode < X[j].subSize(iSub); iNode++) {
								if (localMasterFlag[iNode]) {
										for(int k = 0; k < dim; k++){
												ptr[dim * M *j + i * dim + k] = temp[iNode][k];
										}
										//memcpy(&ptr[dim * M * j + i * dim], &X[j][iSub][iNode], dim * sizeof(double));
										i++;
								}
						}
				}
				//assert(i * dim == M); // should be true
				if( i != M)
						this->com->fprintf(stderr, "double array, dim is %d, %d != %d\n", dim, i, M);
		}
		mem = ptr;
}

/**
 * communication between subdomains to
 */
template<int dim>
template<class Mat>
void ParallelRomExtension<dim>::summonSlaves(double *&mem, Mat &X, const int M,
                                                       const int N) {
		SubDomain **subDomain = this->domain.getSubDomain();
		CommPattern<double> *vecPattern = this->domain.getVecPat();
		DistInfo& distInfo = this->domain.getNodeDistInfo();

		double *ptr = mem;
		for (int j = 0; j < N; j++) {
				int i = 0;
//#pragma omp parallel for
				for (int iSub = 0; iSub < X[j].numLocSub(); iSub++) {
						bool *localMasterFlag = X[j].getMasterFlag(iSub);
						double (*temp)[dim] = X[j].subData(iSub);
						for (int iNode = 0; iNode < X[j].subSize(iSub); iNode++) {
								if (localMasterFlag[iNode]) {
										for(int k = 0; k < dim; k++){
												temp[iNode][k] = ptr[dim * M * j + i * dim + k];
										}
										//memcpy(&X[j][iSub][iNode], &ptr[dim * M * j + i * dim], dim * sizeof(double));
										i++;
								}
						}
						subDomain[iSub]->sndData(*vecPattern, X[j].subData(iSub));
				}
				//assert(i * dim == M); // should be true
				if(i != M) this->com->fprintf(stderr, "dim is %d, i(%d) != M(%d)", dim, i, M);
				vecPattern->exchange();

				for (int iSub = 0; iSub < X[j].numLocSub(); iSub++)
						subDomain[iSub]->addRcvData(*vecPattern, X[j].subData(iSub));
		}
}


template<int dim>
template<class Mat>
void ParallelRomExtension<dim>::summonZombies(double *&mem, Mat &X, const int M,
                                                        const int N) {
		SubDomain **subDomain = this->domain.getSubDomain();
		CommPattern<double> *vecPattern = this->domain.getVecPat();
		DistInfo& distInfo = this->domain.getNodeDistInfo();

		double *ptr = mem;
		for (int j = 0; j < N; j++) {
				int i = 0;
//#pragma omp parallel for
				for (int iSub = 0; iSub < X[j].numLocSub(); iSub++) {
						double (*temp)[dim] = X[j].subData(iSub);
						bool *localMasterFlag = X[j].getMasterFlag(iSub);
						for (int iNode = 0; iNode < X[j].subSize(iSub); iNode++) {
								if (localMasterFlag[iNode]) {
										for(int k = 0; k < dim; k++){
												temp[iNode][k] = ptr[dim * M * j + i * dim + k];
										}
										// memcpy(&X[j][iSub][iNode], &ptr[dim * M * j + i * dim], dim * sizeof(double));
										i++;
								}
						}
				}
				assert(i * dim == M); // should be true
		}
}

template<int dim>
void ParallelRomExtension<dim>::transpose(double* &buff1, double* &buff2, int nrow, int ncol){
		double *ptr1 = buff1;
		double *ptr2 = buff2;
		assert(ptr2 != NULL);
		for(int j = 0; j < ncol; j++){
				for(int i = 0; i < nrow; i++){
						ptr2[i * ncol + j] = ptr1[j * nrow + i];
				}
		}
}