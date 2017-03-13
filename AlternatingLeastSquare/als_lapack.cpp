#include "als_lapack.h"
#include "als_util.h"	// various armadillo matrix function replacement
//#include <string>		// strncpy()
#include <stdlib.h>		// srand(), rand()
#include <iostream>   // std::cout
#include <algorithm>	// std::fill_n, std::copy
#include <cmath>		// sqrt()

/* user supplied all initial conditions and MPI process ranks */
AlternatingLeastSquare::AlternatingLeastSquare(double *_X, unsigned char *_M, double *UT_init,
											   int _nrow, int _ncol, int _dim,
											   MPI_Comm _comm, int _rank/*, char *name*/) :
		X(_X), M(_M), UT(UT_init),
		nrow(_nrow), ncol(_ncol), dim(_dim),
		comm(_comm), rank(_rank) {
	//strncpy(processor_name, name, MPI_MAX_PROCESSOR_NAME);
	this->MM = new double[nrow * ncol];
	char_2_double(M, MM, nrow, ncol);
	this->V = new double[dim * ncol];
	this->error_mat = new double[nrow * ncol];
	this->global_mat_data = new double[dim * dim];
	this->local_mat_data = new double[dim * dim];
	this->global_vec_data = new double[dim];
	this->local_vec_data = new double[dim];
}

/* user supplied nothing; purely for internal testing without MPI */
AlternatingLeastSquare::AlternatingLeastSquare(int _nrow, int _ncol, int _dim,
											   MPI_Comm _comm, int _rank/*, char *name*/) :
		nrow(_nrow), ncol(_ncol), dim(_dim), comm(_comm), rank(_rank) {
	srand(rank);
	// initialize X, M, UT and V
	this->X = new double[nrow * ncol];
	for (int i = 0; i < nrow * ncol; i++) X[i] = rand() % 100;
	this->M = new unsigned char[nrow * ncol];
	for (int i = 0; i < nrow * ncol; i++) M[i] = rand() % 2;
	this->MM = new double[nrow * ncol];
	char_2_double(M, MM, nrow, ncol);
	this->UT = new double[dim * nrow];
	this->V = new double[dim * ncol];
	for (int i = 0; i < dim * nrow; i++) UT[i] = rand() % 70;
	for (int i = 0; i < dim * ncol; i++) V[i] = 0;
	//strncpy(processor_name, name, MPI_MAX_PROCESSOR_NAME);
	this->error_mat = new double[nrow * ncol];
	this->global_mat_data = new double[dim * dim];
	this->local_mat_data = new double[dim * dim];
	this->global_vec_data = new double[dim];
	this->local_vec_data = new double[dim];
	std::cout << "initialization done \n";
}

AlternatingLeastSquare::~AlternatingLeastSquare() {
	if (global_mat_data) delete[] global_mat_data;
	if (local_mat_data) delete[] local_mat_data;
	if(global_vec_data) delete [] global_vec_data;
	if(local_vec_data) delete [] local_vec_data;
}

void AlternatingLeastSquare::display(){
	//std::cout << processor_name << "with rank " << rank;
	if(this->rank == 0){
		std::cout << "\n X is \n";
		print(this->X, this->nrow, this->ncol);
		std::cout << "M is \n";
		print(this->M, this->nrow, this->ncol);
		std::cout << "UT is \n";
		print(this->UT, this->dim, this->nrow);
		std::cout << "V is \n";
		print(this->V, this->dim, this->ncol);
	}
}

// TODO: implement run() with mpi
void AlternatingLeastSquare::run(int maxIterations) {
	this->error_trajectory = new double[maxIterations];
	//X_norm = error();
	//std::cout << "entering run(), initial norm is " << X_norm << std::endl;
	for (int p = 0; p < maxIterations; p++) {
		for (int j = 0; j < ncol; ++j){
			updateColOfV(this->X, this->M, this->UT, this->V, j);
			//std::cout << j << "-th col of V updated\n";
		}
		MPI_Barrier(comm);
		for (int i = 0; i < nrow; ++i) {
			updateRowOfU(this->X, this->M, this->UT, this->V, i);
			//std::cout << i << "-th row of U updated\n";
		}
		MPI_Barrier(comm);
		this->error_trajectory[p] = objectiveFunctionValue();
		if(this->rank == 0) std::cout << "iteration " << p << ", error now is " << this->error_trajectory[p] << std::endl;
		if (p > 0 && this->error_trajectory[p] > this->error_trajectory[p - 1]) break;
	}
	// calculate the projection error at the end.
	relativeProjectionError();
}

/*
 * solves rows of U in paralell by solving normal equation using SVD least square
 */
void AlternatingLeastSquare::updateRowOfU(double *X, unsigned char *M, double *UT, double *V, const int i) {
	std::vector<unsigned int> indices = find_col_index(M, nrow, ncol, i, 0, ncol);
	double *A = join_cols(V, dim, ncol, indices);
	double *b = get_row_slice(X, nrow, ncol, i, indices);
	int n = (int) indices.size();
	//TODO: no need to transpose it actually
	double *res = NULL;
	// only call the following when dim <= n
	int result;
	result = transpose_multiply(A, dim, n, local_mat_data);
	result = matrix_vector_multiply(A, b, dim, n, local_vec_data);
	res = linear_solve_QR(local_mat_data, local_vec_data, dim, dim);
	for(int j = 0; j < dim; j++)
		UT[dim * i + j] = res[j]; // set(UT, dim, nrow, j, i) = res[j];
	delete [] res;
	delete [] A;
	delete [] b;
}

/*
 * solves columns of V in paralell by solving normal equation using SVD least square
 */
void AlternatingLeastSquare::updateColOfV(double *X, unsigned char *M, double *UT, double *V, const int j){
	int result;
	std::vector<unsigned int> indices = find_row_index(M, nrow, ncol, 0, nrow, j);
	double *A = join_cols(UT, dim, nrow, indices);
	double *b = get_col_slice(X, nrow, ncol, j, indices);
	int n = (int) indices.size();
	// solve a normal equation
	//std::cout << " n is " << n << std::endl;
	result = transpose_multiply(A, dim, n, local_mat_data);
	result = matrix_vector_multiply(A, b, dim, n, local_vec_data);
	//std::cout << "transpose and matrix_vector computed\n";

	// TODO: MPI_Allreduce to combine all normal equation into one
	std::fill_n(global_mat_data, dim * dim, 0.0);
	std::fill_n(global_vec_data, dim, 0.0);
	MPI_Barrier(comm);
	MPI_Allreduce(local_mat_data, global_mat_data, dim * dim, MPI_DOUBLE, MPI_SUM, comm);
	MPI_Allreduce(local_vec_data, global_vec_data, dim, MPI_DOUBLE, MPI_SUM, comm);
	MPI_Barrier(comm);
	std::copy(global_vec_data, global_vec_data + dim, local_vec_data);
	std::copy(global_mat_data, global_mat_data + dim * dim, local_mat_data);

	//std::cout << "MPI_Allreduce done \n";
	double * res = linear_solve_QR(local_mat_data, local_vec_data, dim, dim);
	//std::cout << "linear_solve done \n";
	for(int i = 0; i < dim; i++)
		V[dim * j + i] = res[i];//set(V, dim, ncol, i, j) = res[i];
	//std::cout << "assigning element of V done \n";
	delete [] res;
	delete [] A;
	delete [] b;
	//std::cout << "existing updatingColOfV\n";
}

// done: implement mask in error
double AlternatingLeastSquare::objectiveFunctionValue() {
	double *temp_error_mat = new double[nrow * ncol];
	std::copy(X, X + nrow * ncol, temp_error_mat);
	matrix_multiply(UT, V, dim, nrow, ncol, temp_error_mat);
	matrix_matrix_elementwise_multiply(temp_error_mat, MM, nrow, ncol, error_mat);
	local_error = frobenius_norm(error_mat, nrow, ncol);
	local_error = local_error * local_error;
	global_error = 0.0;
	MPI_Barrier(comm);
	MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, comm);
	MPI_Barrier(comm);
	local_error = sqrt(global_error);
	delete [] temp_error_mat;
	return local_error;
}

//todo : seg fault in the function:
void AlternatingLeastSquare::relativeProjectionError() {
	double *temp_matrix = new (std::nothrow) double[nrow * ncol];
	std::cout << " ... rank is " << this->rank << " , nrow is " << nrow << " , ncol is " << ncol << std::endl;
	std::copy(X, X + nrow * ncol, temp_matrix);
	// todo: calculate original col 2-norm first in another function
	matrix_matrix_elementwise_multiply(temp_matrix, MM, nrow, ncol, error_mat);
	double *local_col_norm = new double[ncol];
	double *global_col_norm = new double[ncol];
	// scatter
	std::cout << " ... rank is " << this->rank << " , entering for loop" << std::endl;
	for(int i = 0; i < ncol; i++){
		local_col_norm[i] = frobenius_norm(error_mat + i * nrow, nrow, 1);
		local_col_norm[i] = local_col_norm[i] * local_col_norm[i];
	}
	std::cout << " ... rank is " << this->rank << " , calling std::fill_n" << std::endl;
	std::fill_n(global_col_norm, ncol, 0.0);
	std::cout << " ... rank is " << this->rank << " , entering MPI_Allreduce" << std::endl;
	MPI_Barrier(comm);
	MPI_Allreduce(local_col_norm, global_col_norm, ncol, MPI_DOUBLE, MPI_SUM, comm);
	MPI_Barrier(comm);
	// gather
	for(int i = 0; i < ncol; i++)
		local_col_norm[i] = sqrt(global_col_norm[i]);
	if(this->rank == 0) std::cout << " ... each col sum is computed." << std::endl;
	// todo: calculate difference of col 2-norm next
	std::copy(X, X + nrow * ncol, temp_matrix);
	std::fill_n(error_mat, nrow * ncol, 0.0);
	matrix_multiply(UT, V, dim, nrow, ncol, temp_matrix);
	matrix_matrix_elementwise_multiply(temp_matrix, MM, nrow, ncol, error_mat);
	double *relLocalProjError = new double[ncol];
	double *relGlobalProjError = new double[ncol];
	// scatter
	for(int i = 0; i < ncol; i++){
		relLocalProjError[i] = frobenius_norm(error_mat + i * nrow, nrow, 1);
		relLocalProjError[i] = relLocalProjError[i] * relLocalProjError[i];
	}
	std::fill_n(relGlobalProjError, ncol, 0.0);
	MPI_Barrier(comm);
	MPI_Allreduce(relLocalProjError, relGlobalProjError, ncol, MPI_DOUBLE, MPI_SUM, comm);
	MPI_Barrier(comm);
	// gather
	for(int i = 0; i < ncol; i++)
		relLocalProjError[i] = sqrt(relGlobalProjError[i]);
	if(this->rank == 0) std::cout << "each projected col is computed " << ncol << std::endl;
	delete [] temp_matrix;
	// todo: return this vector somehow and let it be written
	// todo:temporary solution: write it into terminal:
	if(this->rank == 0){
		std::cout << "relative projective errors are [" << std::endl;
		for(int i = 0; i < ncol; i++){
			std::cout << relLocalProjError[i]/local_col_norm[i] << ", ";
		}
		std::cout << "] " << std::endl;
	}
	delete [] local_col_norm;
	delete [] global_col_norm;
	delete [] relGlobalProjError;
	delete [] relLocalProjError;
}