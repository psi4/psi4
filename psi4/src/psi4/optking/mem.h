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

/*! \file mem.h : header for memory allocation functions
    \ingroup optking
*/

#ifndef _opt_mem_h_
#define _opt_mem_h_

namespace opt {

void zero_array(double *A, long int n);
void zero_int_array(int *A, long int n);

double *init_array(long int m);
void free_array(double *f);

int *init_int_array(int m);
void free_int_array(int *f);

bool *init_bool_array(int m);
void free_bool_array(bool *f);

double **init_matrix(long int m, long int n);
void free_matrix(double **A);
void zero_matrix(double **A, long int m, long int n);

int **init_int_matrix(long int m, long int n);
void free_int_matrix(int **A);

bool **init_bool_matrix(long int m, long int n);
void free_bool_matrix(bool **A);

double **unit_matrix(long int m);
void unit_matrix(double **A, long int m);

}

// C versions for BLAS/LAPACK C functions in global namespace
extern "C" {

double *opt_init_array(long int m);
void opt_free_array(double *f);

double **opt_init_matrix(long int m, long int n);
void opt_free_matrix(double **A);

void opt_matrix_copy(double **A, double **B, long int m, long int n);

}

#endif