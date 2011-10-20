/*! \file linear algebra.h : linear algebra functions
    \ingroup optking
*/

#ifndef _opt_linear_algebra_h_
#define _opt_linear_algebra_h_

// C functions called by opt which use BLAS/LAPACK routines
extern "C" {

// matrix multiplications
void opt_matrix_mult(double **A, bool tA, double **B, bool tB, double **C, bool tC,
  int nr, int nl, int nc, bool add);

// eigenvector/eigenvalues
bool opt_symm_matrix_eig(double **A, int dim, double *evals);

}

namespace opt {

// invert a symmetric matrix
double ** symm_matrix_inv(double **A, int dim, bool redundant=true);

// allocate memory and return a copy of a matrix
double ** matrix_return_copy(double **A, int nr, int nc);
bool ** matrix_return_copy(bool **A, int nr, int nc);

void matrix_copy(double **from, double **to, int nr, int nc);

void array_copy(double *from, double *to, int n);
double array_dot(double *v1, double *v2, int n);
double array_norm(double *v1, int n);
void array_normalize(double *v1, int n);
void array_scm(double *v1, double a, int n);
double array_abs_max(double *v1, int n);
double array_rms(double *v1, int n);
// Compute matrix ^1/2 or ^-1/2 if inverse=true
void matrix_root(double **A, int dim, bool inverse);

}

#endif

