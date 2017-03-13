#ifndef ALS_UTIL_H
#define ALS_UTIL_H

#include <vector>

inline unsigned char get(unsigned char *M, int m, int n, int r, int c);
inline double get(double *M, int m, int n, int r, int c);
inline double& set(double *M, int m, int n, int r, int c);
std::vector<unsigned int> find_col_index(const unsigned char* mem, int M, int N, int RowIndex, int startColIndex, int endColIndex);
std::vector<unsigned int> find_row_index(const unsigned char* mem, int M, int N, int startRowIndex, int endRowIndex, int ColIndex);
double *join_cols(double *mem, int M, int N, std::vector<unsigned int> colIndices);
double *get_row_slice(double *mem, int M, int N, int RowIndex, std::vector<unsigned int> colIndices);
double *get_col_slice(double *mem, int M, int N, int ColIndex, std::vector<unsigned int> rowIndices);
double* transpose(double* A, int M, int N);

double *linear_solve_LS(double *A, double *b, int nrow, int ncol);
double *linear_solve_QR(double *A, double *b, int nrow, int nocl);
double *linear_solve_SVD(double *A, double *b, int nrow, int ncol);

void char_2_double(unsigned char *X, double *Y, int M, int N);
int transpose_multiply(double *A, int M, int N, double *result);
int matrix_multiply(double *A, double *B, int M, int N, int L, double *result);
int matrix_vector_multiply(double *A, double *x, int M, int N, double *result);
int matrix_matrix_elementwise_multiply(double *A, double *B, int M, int N, double *result);
double frobenius_norm(double *X, int M, int N);

void print(double *mat, int nrow, int ncol);
void print(unsigned char *mat, int nrow, int col);

#endif
