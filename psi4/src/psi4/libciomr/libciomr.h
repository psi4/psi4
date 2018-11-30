/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*
** Declarations for functions found in libciomr.a
**
** C. David Sherrill and T. Daniel Crawford
**
**
*/

#ifndef _psi_src_lib_libciomr_libciomr_h_
#define _psi_src_lib_libciomr_libciomr_h_

#include <cstdio>
#include <string>

#include "psi4/pragma.h"

namespace psi {

int psi_start(FILE **infile, FILE **outfile, char **psi_file_prefix, int argc, char *argv[], int overwrite_output);
int psi_stop(FILE *infile, FILE *outfile, char *psi_file_prefix);
char *psi_ifname();
char *psi_ofname();
char *psi_fprefix();

void balance(double **a, int n);

void eigsort(double *d, double **v, int n);
void eivout(double **a, double *b, int m, int n, std::string out);
void mosort(double *d, double **v, int *sym, int nso, int nmo);

void flin(double **a, double *b, int in, int im, double *det);
void free_matrix(double **array, size_t size);
double *init_array(size_t size);
double **init_matrix(size_t rows, size_t cols);

void lubksb(double **a, int n, int *indx, double *b);
void ludcmp(double **a, int n, int *indx, double *d);

/* Functions under mat_to_arr.c */
void mat_to_arr(double **a, double *b, int m, int n);
void arr_to_mat(double **a, double *b, int m, int n);

// void mmult(double **AF, int ta, double **BF, int tb, double **CF, int tc,
//            int nr, int nl, int nc, int add) ;
void print_array(double *a, int m, std::string out);
void print_mat(double **a, int rows, int cols, std::string out);

void rsp(int nm, int n, int nv, double *array, double *evals, int matz, double **evecs, double toler);
void sq_rsp(int nm, int n, double **array, double *evals, int matz, double **evecs, double toler);
void sq_to_tri(double **bmat, double *amat, int size);

/* Functions under tri_to_block.c */
void tri_to_block(double *a, double **b, int num_ir, int *num_so, int *ioff);
void block_to_tri(double *a, double **b, int num_ir, int *num_so, int *ioff);

void tri_to_sq(double *amat, double **bmat, int size);

/* Functions under tstart.c */
PSI_API
void tstart();
PSI_API
void tstop();

/* Functions in zero.c */
void zero_arr(double *a, int size);
void zero_mat(double **a, int rows, int cols);

/* Functions in int_array.c */
PSI_API int *init_int_array(int size);
void zero_int_array(int *a, int size);
PSI_API int **init_int_matrix(int rows, int cols);
void free_int_matrix(int **array);
void zero_int_matrix(int **array, int rows, int cols);
void print_int_mat(int **a, int m, int n, std::string out);

/* Functions in long_int_array.c */
long int *init_long_int_array(int size);
void zero_long_int_array(long int *a, int size);
long int **init_long_int_matrix(int rows, int cols);
void free_long_int_matrix(long int **array);
void zero_long_int_matrix(long int **array, int rows, int cols);
void print_long_int_mat(long int **a, int m, int n, std::string out);

/* Functions in block_matrix.c */
PSI_API double **block_matrix(size_t n, size_t m, bool mlock = false);
void free_block(double **array);

/* Functions in fndcor */
void fndcor(long int *maxcrb, std::string out_fname);
}  // namespace psi

#endif /* header guard */
