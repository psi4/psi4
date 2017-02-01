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

bool opt_asymm_matrix_eig(double **A, int dim, double *evals);

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
double array_max(double *v1, int n);
double array_rms(double *v1, int n);
// Compute matrix ^1/2 or ^-1/2 if inverse=true
void matrix_root(double **A, int dim, bool inverse);

}

#endif