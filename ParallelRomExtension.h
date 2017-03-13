//
// Created by lei on 1/5/17.
//

#ifndef PROJECT_PARALLELROMEXTENSION_H
#define PROJECT_PARALLELROMEXTENSION_H

#include <ParallelRom.h>
#include <vector>

template <int dim>
class ParallelRomExtension : public ParallelRom<dim> {
	//<! extends parallelROM to have alternating least square to calculate basis
private:
	/**
	 * strips distSVec of slave nodes, put it into a continuous region of memory, pointed by mem
	 */
		template<class D, class Mat>
		void freeSlaves(D *&mem, const Mat &X, const int M, const int N);

		/**
		 * Refills distSVec with mem, a continuous region of memory, put it in X
		 */
		template<class Mat>
		void summonSlaves(double *&mem, Mat &X, const int M, const int N);
		template<class Mat>
		void summonZombies(double *&mem, Mat &X, const int M, const int N);

		/**
		 * convert a matrix stored in a continuous region of memory from row-major to column-major
		 */
		void transpose(double* &buff1, double* &buff2, int nrow, int ncol);


		/**
			* returns the number of master nodes in a vector
			* res[i] = # of master nodes in subdomain i
			*/
		std::vector<int> countMasters(const DistInfo &distinfo);

public:
	//<! constructor and destructor
	ParallelRomExtension(Domain &, Communicator *, const DistInfo&);
	~ParallelRomExtension();
	//<! only method available
		template<class Mat1, class Mat2>
		void parallelALS(const Mat1 &X, const Mat2 &M, Mat1 &UT, int maxIts);
};


#endif //PROJECT_PARALLELROMEXTENSION_H
