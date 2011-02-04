/*! \file
    \ingroup OPT10
    \brief memory allocation
*/

#include <cstdlib>
#include "mem.h"

namespace opt {

void zero_array(double *A, long int n) {
  for (int i=0; i<n; ++i)
    A[i]= 0.0;
}

void zero_matrix(double **A, long int m, long int n) {
  double *A1 = A[0];
  for(int i=0; i<m*n; ++i)
    A1[i] = 0.0;
}

void zero_int_array(int *A, long int n) {
  for (int i=0; i<n; ++i)
    A[i]= 0.0;
}

void zero_bool_array(bool *A, long int n) {
  for (int i=0; i<n; ++i)
    A[i]= false;
}

double *init_array(long int m) {
  if (!m) return NULL;
  double *A = (double *) malloc(m * sizeof(double));
  zero_array(A,m);
  return A;
}

int *init_int_array(int m) {
  int *A = (int *) malloc(m * sizeof(int));
  zero_int_array(A,m);
  return A;
}

bool *init_bool_array(int m) {
  bool *A = (bool *) malloc(m * sizeof(bool));
  for (int i=0; i<m; ++i)
     A[i] = false;
  return A;
}

void free_array(double *f) {
  if (f == NULL) return;
  free(f);
}

void free_int_array(int *f) {
  free(f);
}

void free_bool_array(bool *f) {
  free(f);
}

// allocate memory block of doubles for 2D matrix
double **init_matrix(long int m, long int n) {
  double **A = NULL;
  double *B = NULL;
  long int i;

  if(m<=0 || n<=0) return((double **) NULL);

  A = (double **) malloc(m * (long int)sizeof(double *));
  B = (double *) malloc(m*n * (long int)sizeof(double));

  if ((A == NULL) || (B == NULL)) throw("init_matrix : allocation error.");

  zero_array(B, m*n);

  for (i = 0; i < m; i++)
    A[i] = &(B[i*n]);

  return(A);
}

double **unit_matrix(long int m) {
  double **A = init_matrix(m, m);
  for (int i=0; i<m; ++i)
    A[i][i] = 1.0;
  return A;
}

void unit_matrix(double **A, long int m) {
  zero_matrix(A, m, m);
  for (int i=0; i<m; ++i)
    A[i][i] = 1.0;
}

// free memory block of doubles
void free_matrix(double **A) {
  if(A == NULL) return;
  free(A[0]);
  free(A);
}

int **init_int_matrix(long int m, long int n) {
  int **A = NULL;
  int *B = NULL;
  long int i;

  if(m<=0 || n<=0) return((int **) NULL);

  A = (int **) malloc(m * (long int)sizeof(int *));
  B = (int *) malloc(m*n * (long int)sizeof(int));

  if ((A == NULL) || (B == NULL)) throw("init_matrix : allocation error.");

  zero_int_array(B, m*n);

  for (i = 0; i < m; i++)
    A[i] = &(B[i*n]);

  return(A);
}

void free_int_matrix(int **A) {
  if(A == NULL) return;
  free(A[0]);
  free(A);
}

bool **init_bool_matrix(long int m, long int n) {
  bool **A = NULL;
  bool *B = NULL;
  long int i;

  if(m<=0 || n<=0) return((bool **) NULL);

  A = (bool **) malloc(m * (long int)sizeof(bool *));
  B = (bool *) malloc(m*n * (long int)sizeof(bool));

  if ((A == NULL) || (B == NULL)) throw("init_bool_matrix : allocation error.");

  zero_bool_array(B, m*n);

  for (i = 0; i < m; i++)
    A[i] = &(B[i*n]);

  return(A);
}

void free_bool_matrix(bool **A) {
  if(A == NULL) return;
  free(A[0]);
  free(A);
}

}

extern "C" {

double *opt_init_array(long int m) {
  double *A = (double *) malloc(m * sizeof(double));
  for (int i=0; i<m; ++i)
    A[i]= 0.0;
  return A;
}

void opt_free_array(double *f) {
  if (f == NULL) return;
  free(f);
}

double **opt_init_matrix(long int m, long int n) {
  double **A = NULL;
  double *B = NULL;
  long int i;

  if(m<=0 || n<=0) return((double **) NULL);

  A = (double **) malloc(m * (long int)sizeof(double *));
  B = (double *) malloc(m*n * (long int)sizeof(double));

  if ((A == NULL) || (B == NULL)) throw("init_matrix : allocation error.");

  for (i=0; i<m*n; ++i)
    B[i] = 0.0;

  for (i = 0; i < m; i++)
    A[i] = &(B[i*n]);

  return(A);
}

// free memory block of doubles
void opt_free_matrix(double **A) {
  if(A == NULL) return;
  free(A[0]);
  free(A);
}

}
