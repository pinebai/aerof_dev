#ifndef ALS_LAPACK_H
#define ALS_LAPACK_H

#include <mpi.h>		// MPI_Comm, MPI_MAX_PROCESSOR_NAME

class AlternatingLeastSquare {
    /** MPI process related variable */
    MPI_Comm comm;
    //char processor_name[MPI_MAX_PROCESSOR_NAME];
    int rank;
    /** all related variables */
    int nrow, ncol, dim;
    double *X;
    unsigned char *M;
    double *MM; // double version of M, for calculating error
    double *UT;
    double *V;
    double X_norm;
    double *error_trajectory;
    /** MPI communication memory */
    double *global_vec_data;
    double *local_vec_data;
    double *global_mat_data;
    double *local_mat_data;
    double global_error;
    double local_error;
    double *error_mat;
    /** helper routine for performing ALS */
    void updateRowOfU(double *X, unsigned char *M, double *U, double *V, const int i);

    void updateColOfV(double *X, unsigned char *M, double *U, double *V, const int j);

public:
    /* user supplied all initial conditions and MPI process ranks */
    AlternatingLeastSquare(double *_X, unsigned char *_M, double *UT_init,
                           int _nrow, int _ncol, int _dim,
                           MPI_Comm _comm, int _rank/*, char *name*/) ;

    /* user supplied nothing; purely for internal testing without MPI */
    AlternatingLeastSquare(int _nrow, int _ncol, int _dim,
                           MPI_Comm _comm, int _rank/*, char *name*/);

    ~AlternatingLeastSquare();

    void run(int maxIterations);

    /** debugging output */
    void updateU(const int i) { updateRowOfU(X, M, UT, V, i); }

    void updateV(const int j) { updateColOfV(X, M, UT, V, j); }

    void display();

    /** objective function value */
    double objectiveFunctionValue();

    void relativeProjectionError();
};

#endif