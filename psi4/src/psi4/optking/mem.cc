/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file mem.cc : memory allocation
    \ingroup optking
*/

#include <cstdlib>
#include "mem.h"
#include "opt_except.h"

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
  if(A == NULL)
    throw(opt::INTCO_EXCEPT("init_array : allocation error."));

  zero_array(A,m);
  return A;
}

int *init_int_array(int m) {
  int *A = (int *) malloc(m * sizeof(int));
  if(A == NULL)
    throw(opt::INTCO_EXCEPT("init_int_array : allocation error."));

  zero_int_array(A,m);
  return A;
}

bool *init_bool_array(int m) {
  bool *A = (bool *) malloc(m * sizeof(bool));
  if(A == NULL)
    throw(opt::INTCO_EXCEPT("init_bool_array : allocation error."));

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

  if ((A == NULL) || (B == NULL))
    throw(opt::INTCO_EXCEPT("init_matrix : allocation error."));

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

  if ((A == NULL) || (B == NULL)) throw(opt::INTCO_EXCEPT("init_int_matrix : allocation error."));

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

  if ((A == NULL) || (B == NULL)) throw(opt::INTCO_EXCEPT("init_bool_matrix : allocation error."));

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
  if(A == NULL)
    throw(opt::INTCO_EXCEPT("opt_init_array : allocation error."));

  for (int i=0; i<m; ++i)
    A[i]= 0.0;
  return A;
}

void opt_free_array(double *f) {
  if (f == NULL) return;
  free(f);
}

void opt_matrix_copy(double **from, double **to, long int nr, long int nc) {
  double *from1 = from[0];
  double *to1 = to[0];
  for (long int i=0; i<nr*nc; ++i)
    to1[i] = from1[i];
}

double **opt_init_matrix(long int m, long int n) {
  double **A = NULL;
  double *B = NULL;
  long int i;

  if(m<=0 || n<=0) return((double **) NULL);

  A = (double **) malloc(m * (long int)sizeof(double *));
  B = (double *) malloc(m*n * (long int)sizeof(double));

  if ((A == NULL) || (B == NULL))
    throw(opt::INTCO_EXCEPT("opt_init_matrix : allocation error."));

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
