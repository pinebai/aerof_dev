//
// Created by lei on 4/14/16.
// replacement for armadillo functionality
//

#include "als_util.h"
#include "lapack.h"
#include "blas.h"
#include <iostream>
#include <algorithm>    // std::max, std::min
#include <math.h>       //log()

/**
 * assumes matrix are stored in col-major format
 */
inline
unsigned char get(unsigned char *M, int m, int n, int r, int c) {
    return M[c * m + r];
}

inline
double get(double *M, int m, int n, int r, int c) {
    return M[c * m + r];
}

inline
double& set(double *M, int m, int n, int r, int c) {
    return M[c * m + r];
}

/**
 * give a matrix M pointed by mem, of size M-by-N
 * return the number of nonzero elements in the range
 * M[a:b, c:d] (matlab notation)
 */
std::vector<unsigned int> find_col_index(const unsigned char *mem, int M, int N,
                                         int RowIndex, int startColIndex, int endColIndex) {
    std::vector<unsigned int> result;
    unsigned char *Mat = (unsigned char *) mem;
    for (int j = startColIndex; j < endColIndex; j++) {
        if (get(Mat, M, N, RowIndex, j))
            result.push_back(j);
    }
    return result;
}

std::vector<unsigned int> find_row_index(const unsigned char* mem, int M, int N,
                                         int startRowIndex, int endRowIndex, int ColIndex) {
    std::vector<unsigned int> result;
    unsigned char *Mat = (unsigned char *) mem;
    for (int i = startRowIndex; i < endRowIndex; i++) {
        if (get(Mat, M, N, i, ColIndex))
            result.push_back(i);
    }
    return result;
}

/**
 * given a matrix X pointed by mem, of size M-by-N
 * and a vector of indices, from smallest to largest, say v
 * return X[:, v] (matlab notation)
 */
double *join_cols(double *mem, int M, int N, std::vector<unsigned int> colIndices) {
    int n = colIndices.size();
    double *result = new double[M * n];
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < M; i++) {
            set(result, M, n, i, j) = get(mem, M, N, i, colIndices[j]);
        }
    }
    return result;
}

/**
 * given a matrix X pointed by mem, of size M-by-N
 * and an integer r, and a set of indices v,
 * return X[r, v]' (matlab notation)
 */

double *get_row_slice(double *mem, int M, int N, int RowIndex,
                      std::vector<unsigned int> colIndices) {
    int n = colIndices.size();
    double *result = new double[n];
    for (int j = 0; j < n; j++) {
        result[j] = get(mem, M, N, RowIndex, colIndices[j]);
    }
    return result;
}

double *get_col_slice(double *mem, int M, int N, int ColIndex,
                      std::vector<unsigned int> rowIndices) {
    int m = rowIndices.size();
    double *result = new double[m];
    for (int i = 0; i < m; i++) {
        result[i] = get(mem, M, N, rowIndices[i], ColIndex);
    }
    return result;
}

/**
 * transpose a matrix A to A^T, A^T in col-major order
 * essentially store A in row-major order
 */
double *transpose(double *A, int M, int N){
    double *result = new double[M * N];
    for(int j = 0; j < N; j++){
        for( int i = 0; i < M; i++){
            set(result, N, M, j, i) = get(A, M, N, i, j);
        }
    }
    return result;
}


/**
 * solve Ax = b via
 * external call to lapack's least square solver
 * double precision least square solver using SVD: dgelsd()
 * result is stored in first ncol entries of b
 */
double* linear_solve_LS(double *A, double *b, int nrow, int ncol) {
    int m = nrow;
    int n = ncol;
    int nrhs = 1;
    int lda = m;
    int min_mn = std::min(m, n);
    int max_mn = std::max(m, n);
    int ldb = max_mn;
    double *res = new double [ldb];
    std::fill_n(res, ldb, 0.0);
    std::copy(b, b + nrow, res);
    char *transA = "N";
    int info;
    // work is to be constructed later
    /*
    double work_dummy;
    int info;     // query work size
    char *transA = "N";
    F77_NAME(dgels)(transA, &m, &n, &nrhs,
                    A, &lda, res, &ldb,
                    &work_dummy, &lwork, &info);
    lwork = (int) work_dummy;
    double *work = new double [lwork + 1];
    */
    // alternative work size computation
    int lwork = 3 * std::max(1, min_mn + std::max(min_mn, nrhs));
    double *work = new double [lwork + 1];
    // real computation
    F77_NAME(dgels)(transA, &m, &n, &nrhs,
                    A, &lda, res, &ldb,
                    work, &lwork, &info);
    delete [] work;
    return res;
}

double *linear_solve_QR(double *A, double *b, int nrow, int ncol) {
    int m = nrow;
    int n = ncol;
    int nrhs = 1;
    //double *A;
    int lda = nrow;
    //double *b;
    int ldb = m > n ? m : n;
    double *res = new double [ldb];
    std::fill_n(res, ldb, 0.0);
    std::copy(b, b + nrow, res);
    int *jpvt = new int [n];
    double rcond = 1e-8;
    int rank;
    int lwork = -1;
    double work_dummy;
    int info;
    // query work size
    F77_NAME(dgelsy)(&m, &n, &nrhs,
                     A, &lda, res, &ldb,
                     jpvt, &rcond, &rank,
                     &work_dummy, &lwork, &info);
    lwork = (int) work_dummy;
    double *work = new double [lwork + 1];
    // real computation
    F77_NAME(dgelsy)(&m, &n, &nrhs,
                     A, &lda, res, &ldb,
                     jpvt, &rcond, &rank,
                     work, &lwork, &info);
    delete [] work;
    return res;
}

double *linear_solve_SVD(double *A, double *b, int nrow, int ncol) {
    int m = nrow;
    int n = ncol;
    int nrhs = 1;
    int l = (m < n ? m : n);
    //double *A;
    int lda = nrow;
    //double *b;
    int minmn = std::min(m, n);
    int maxmn = std::max(m, n);
    int ldb = maxmn;
    double *res = new double [ldb];
    std::fill_n(res, ldb, 0.0);
    std::copy(b, b + nrow, res);
    double* s = new double [l];
    double rcond = 1e-8;
    int rank = 0;
    // see
    int ispec = 9;
    char* name = "DGELSD";
    char* opt = "";
    int smlsiz = std::max(25, F77_NAME(ilaenv)(&ispec, name, opt, &m, &n, &nrhs, &lda) );
    int nlvl = std::max( 1, 1 + (int)(log((double)minmn/(double)smlsiz)/(double) 0.69314718055994530942) );
    int liwork = std::max(1, 3 * minmn * nlvl + 11 * minmn);
    int iwork[liwork];
    //
    int lwork_query = -1;
    double work_query;
    int info = 0;
    // query work size
    F77_NAME(dgelsd)(&m, &n, &nrhs,
                     A, &lda, res, &ldb,
                     s, &rcond, &rank,
                     &work_query, &lwork_query, iwork, &info);
    int lwork = (int) work_query;
    double *work = new double [lwork];
    // real computation
    F77_NAME(dgelsd)(&m, &n, &nrhs,
                     A, &lda, res, &ldb,
                     s, &rcond, &rank,
                     work, &lwork, iwork, &info);
    delete [] work;
    return res;
}

/**
 * given matrix A, return
 * A A^T
 */
int transpose_multiply(double *A, int nrow, int ncol, double *C) {
    double *B = new double[nrow * ncol];
    std::copy(A, A + (nrow * ncol), B);
    int m = nrow;
    int n = nrow;
    int k = ncol;
    double alpha = 1.0;
    double beta = 0.0;
    int lda = m;
    int ldb = n;
    std::fill_n(C, nrow * nrow, 0.0);
    int ldc = m;

    const char* transA = "N";
    const char* transB = "T";
    F77_NAME(dgemm)(transA, transB,
                    &m, &n, &k,
                    &alpha, A, &lda,
                    B, &ldb, &beta,
                    C, &ldc);
    delete [] B;
    return 0;
}

/**
 * given matrix A and vector b, return
 * v = Ab
 */
int matrix_vector_multiply(double *A, double *b, int nrow, int ncol, double *v) {
    const char* trans = "N";
    int m = nrow;
    int n = ncol;
    double alpha = 1.0;
    int lda = m;
    //double *x = b;
    int incb = 1;
    double beta = 0.0;
    std::fill_n(v, nrow, 0.0);
    int incv = 1;

    F77_NAME(dgemv)(trans,
                    &m, &n,
                    &alpha, A, &lda,
                    b, &incb,
                    &beta, v, &incv);
    return 0;
}

/**
 * A is a-by-b, B is a-by-c
 * result <- A^T B - result, a b-by-c matrix
 */
int matrix_multiply(double *A, double *B, int a, int b, int c, double *result){
    char *transa = "T";
    char *transb = "N";
    int m = b;
    int n = c;
    int k = a;
    double alpha = 1.0;
    int lda = k;
    int ldb = k;
    double beta = -1.0;
    int ldc = m;

    F77_NAME(dgemm)(transa, transb,
                    &m, &n, &k,
                    &alpha, A, &lda,
                    B, &ldb,
                    &beta, result, &ldc);
    return 0;
}


/**
 * A is a-by-b, B is a-by-b
 * result <= A % B, where % is elementwise multiplication
 */
int matrix_matrix_elementwise_multiply(double *A, double *B, int M, int N, double *result){
    char *uplo = "L";
    int n = M * N;
    int k = 0; // only diagonal is provided
    double alpha = 1.0;
    int lda = 1;
    int incx = 1;
    double beta = 0.0;
    int incy = 1;
    F77_NAME(dsbmv)(uplo,
                    &n, &k, &alpha, A,
                    &lda, B, &incx,
                    &beta, result, &incy);
    return 0;
}

void char_2_double(unsigned char *X, double *Y, int M, int N){
    for(int i = 0; i < M * N; i++) Y[i] = (double) X[i];
}

/**
 * frobenius norm of matrix
 */
double frobenius_norm(double *X, int M, int N){
    char *trans = "F";
    int m = M;
    int n = N;
    int lda = m;
    double work;

    double result = F77_NAME(dlange)(trans, &m, &n, X, &lda, &work);
    return result;
}
/**
 * print out matrix, assuming col-major order
 */
void print(double *mat, int nrow, int ncol){
    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            std::cout << get(mat, nrow, ncol, i, j) << ", ";
        }
        std::cout << std::endl;
    }
}

void print(unsigned char *mat, int nrow, int ncol){
    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            char c = get(mat, nrow, ncol, i, j);
            int symbol = c == 1? 1 : 0;
            std::cout << symbol  << ", ";
        }
        std::cout << std::endl;
    }
}