/*! \file    linear_algebra.cc
    \ingroup optking
    \brief   linear algebra functions which call lapack and blas are in
             global namespace and preceded with opt
*/

#include "linear_algebra.h"

#include <cstdlib>
#include <cmath>
#include "mem.h"

#define EXTERN
#include "globals.h"

extern "C" {

// defines correct format of BLAS/LAPACK functions
#define F_DGEMM dgemm_
#define F_DSYEV dsyev_

// declations of BLAS/LAPACK routines
extern void F_DGEMM(char *transa, char *transb, int *m, int *n, int *k,
  double *alpha, double *A, int *lda, double *B, int *ldb,
  double *beta, double *C, int *ldc);

extern int F_DSYEV(char *, char *, int *, double *, int *, double *,
  double *, int *, int *);

// C functions called from opt that use the BLAS/LAPACK routines

/*
  matrix multiplication using DGEMM :  A * B += C 
  tA (tB) indicate if A (B) is transposed;
  nr = num of rows
  nl = num of links
  nc = num of cols
  add = indicates whether to add to C or overwrite
*/
void opt_matrix_mult(double **A, bool tA, double **B, bool tB, double **C, bool tC,
    int nr, int nl, int nc, bool add) {
  double alpha = 1.0;
  double beta = (add ? 1.0 : 0.0);
  char FtA, FtB;
  int nca, ncb;

  if (!nr || !nl || !nc) return;

  if (!tC)  { // reverse A and B to account for different C/Fortran stride
    FtA = (tA ? 'T' : 'N');
    FtB = (tB ? 'T' : 'N');
    nca = (tA ? nr : nl );
    ncb = (tB ? nl : nc );

    //F_DGEMM(&FtB, &FtA, &nc, &nr, &nl, &alpha, B[0], &ncb, A[0], &nca, &beta, C[0], &nr);
    F_DGEMM(&FtB, &FtA, &nc, &nr, &nl, &alpha, B[0], &ncb, A[0], &nca, &beta, C[0], &nc);
  }
  else { // C is transposed, so compute B^t * A^t = C^t
    FtA = (tA ? 'N' : 'T');
    FtB = (tB ? 'N' : 'T');
    nca = (tA ? nr : nl );
    ncb = (tB ? nl : nc );

    F_DGEMM(&FtA, &FtB, &nr, &nc, &nl, &alpha, A[0], &nca, B[0], &ncb, &beta, C[0], &nr);
  }
  return;
}

/*
  Compute eigenvectors and eigenvalues of a real, symmetric matrix.
  Eigenvectors are stored as rows of A.
  eigenvalues returned in ascending order.
*/
bool opt_symm_matrix_eig(double **A, int dim, double *evals) {
  char cv = 'V';  // compute evals and evects
//  char cl = 'L';  // lower triangle (upper in C) is necessary
  char cl = 'U';  // upper triangle (lower triangle in C) is necessary
  int lwork = 3 * dim;
  int rval, i, j;
  double tval;
  double * work = opt_init_array(lwork);

  F_DSYEV(&cv, &cl, &dim, A[0], &dim, evals, work, &lwork, &rval);

  opt_free_array(work);
  if (rval != 0) return false;
  return true;
}

}

namespace opt {

/* 
  Invert a real, symmetric matrix.  If "redundant" == true, then a 
  generalized inverse is permitted.  A is preserved.  Requires lower
  triangle (in C) of matrix.
*/
double ** symm_matrix_inv(double **A, int dim, bool redundant) {
  int i, j;
  double det = 1.0;
  double * evals = init_array(dim);
  double ** A_evects = matrix_return_copy(A, dim, dim);

  if (dim <= 0) return ( (double **) NULL);

  if (! opt_symm_matrix_eig(A_evects, dim, evals) )
    throw(INTCO_EXCEPT("symm_matrix_inv : opt_symm_matrix_eig could not diagonalize"));
  //for (i=0; i<dim; ++i) fprintf(outfile, "evals[%d] = %15.10lf\n", i, evals[i]);

  for (i=0;i<dim;++i)
    det *= evals[i];

  if (!redundant && fabs(det) < 1E-10)
    throw(INTCO_EXCEPT("symm_matrix_inv : Non-generalized inverse of matrix failed"));

  double ** A_inv = init_matrix(dim, dim);

  if (redundant) {
    for (i=0;i<dim;++i)
      if (fabs(evals[i]) > Opt_params.redundant_eval_tol)
        A_inv[i][i] = 1.0/evals[i];
  }
  else {
    for (i=0;i<dim;++i)
      A_inv[i][i] = 1.0/evals[i];
  }

  double ** A_temp = init_matrix(dim, dim);

  // A^-1 = P^t D^-1 P
  opt_matrix_mult(A_inv, 0, A_evects, 0, A_temp, 0, dim, dim, dim,0);
  opt_matrix_mult(A_evects, 1, A_temp, 0, A_inv, 0, dim, dim, dim, 0);

  free_matrix(A_temp);
  free_array(evals);
  free_matrix(A_evects);

  return A_inv;
}

// allocate memory and return a copy of a matrix
double ** matrix_return_copy(double **A, int nr, int nc) {
  double ** A_new = init_matrix(nr, nc);
  int i, j;
  for(i=0; i<nr; ++i)
    for(j=0; j<nc; ++j)
      A_new[i][j] = A[i][j];
  return A_new;
}

// allocate memory and return a copy of a matrix
bool ** matrix_return_copy(bool **A, int nr, int nc) {
  bool ** A_new = init_bool_matrix(nr, nc);
  int i, j;
  for(i=0; i<nr; ++i)
    for(j=0; j<nc; ++j)
      A_new[i][j] = A[i][j];
  return A_new;
}

void array_copy(double *v_from, double *v_to, int n) {
  for (int i=0; i<n; ++i)
    v_to[i] = v_from[i];
}

void matrix_copy(double **from, double **to, int nr, int nc) {
  double *from1 = from[0];
  double *to1 = to[0];
  for (int i=0; i<nr*nc; ++i) 
    to1[i] = from1[i];
}

double array_dot(double *v1, double *v2, int n) {
  double dot = 0;
  for (int i=0; i<n; ++i)
    dot += v1[i] * v2[i];
  return dot;
}

void array_normalize(double *v1, int n) {
  double norm = sqrt(array_dot(v1, v1, n));
  array_scm(v1, 1/norm, n);
}

double array_norm(double *v1, int n) {
  return sqrt(array_dot(v1, v1, n));
}

void array_scm(double *v1, double a, int n) {
  for (int i=0; i<n; ++i)
    v1[i] *= a;
}

double array_abs_max(double *v1, int n) {
  double max = 0.0;
  for (int i=0; i<n; ++i)
    if (fabs(v1[i]) > max) max = fabs(v1[i]);
  return max;
}

double array_rms(double *v1, int n) {
  double rms = array_dot(v1, v1, n);
  rms = rms / ((double) n);
  return sqrt(rms);
}

// Compute matrix ^1/2 or ^-1/2 if inverse=true
void matrix_root(double **A, int dim, bool inverse) {
  double **V = matrix_return_copy(A, dim, dim);
  double *A_evals = init_array(dim);

  opt_symm_matrix_eig(V, dim, A_evals);

  if (inverse) {
    for(int k=0; k<dim; k++)
      if(fabs(A_evals[k]) > Opt_params.redundant_eval_tol)
        A_evals[k] = 1.0/A_evals[k];
  }

  for(int k=0; k<dim; k++) {
    if(A_evals[k] > 0)
      A_evals[k] = sqrt(A_evals[k]);
    else
    {
      A_evals[k] = 0;
    }
  }

  zero_matrix(A, dim, dim);

  for(int i=0; i<dim; i++)
    for(int j=0; j<dim; j++)
      for(int k=0; k<dim; k++)
        A[i][j] += V[k][i] * A_evals[k] * V[k][j];

  free_matrix(V);
  free_array(A_evals);
}

} // namespace:: opt

