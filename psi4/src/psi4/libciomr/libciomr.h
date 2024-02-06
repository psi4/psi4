/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

PSI_API int psi_start(FILE **infile, FILE **outfile, char **psi_file_prefix, int argc, char *argv[],
                      int overwrite_output);
PSI_API int psi_stop(FILE *infile, FILE *outfile, char *psi_file_prefix);
PSI_API char *psi_ifname();
PSI_API char *psi_ofname();
PSI_API char *psi_fprefix();

PSI_API void balance(double **a, int n);

PSI_API void eigsort(double *d, double **v, int n);
PSI_API void eivout(double **a, const double *b, int m, int n, std::string out);
PSI_API void mosort(double *d, double **v, int *sym, int nso, int nmo);

PSI_API void flin(double **a, double *b, int in, int im, double *det);
PSI_API void free_matrix(double **array, size_t size);
PSI_API double *init_array(size_t size);
PSI_API double **init_matrix(size_t rows, size_t cols);

PSI_API void lubksb(double **a, int n, int *indx, double *b);
PSI_API void ludcmp(double **a, int n, int *indx, double *d);

/* Functions under mat_to_arr.c */
PSI_API void mat_to_arr(double **a, double *b, int m, int n);
PSI_API void arr_to_mat(double **a, double *b, int m, int n);

// void mmult(double **AF, int ta, double **BF, int tb, double **CF, int tc,
//            int nr, int nl, int nc, int add) ;
PSI_API void print_array(double *a, int m, std::string out);
PSI_API void print_mat(double **a, int rows, int cols, std::string out);
PSI_API void sq_to_tri(double **bmat, double *amat, int size);

[[nodiscard]] int DSYEV_ascending(const int N, const double *const *const array, double *e_vals,
                                  double *const *const e_vecs = nullptr);
[[nodiscard]] int DSYEV_descending(const int N, const double *const *const array, double *e_vals,
                                   double *const *const e_vecs = nullptr);

/* Functions under tri_to_block.c */
PSI_API void tri_to_block(double *a, double **b, int num_ir, int *num_so, int *ioff);
PSI_API void block_to_tri(double *a, double **b, int num_ir, int *num_so, int *ioff);

PSI_API void tri_to_sq(double *amat, double **bmat, int size);

/* Functions under tstart.c */
void tstart();
void tstop();

/* Functions in zero.c */
PSI_API void zero_arr(double *a, int size);
PSI_API void zero_mat(double **a, int rows, int cols);

/* Functions in int_array.c */
PSI_API int *init_int_array(int size);
PSI_API void zero_int_array(int *a, int size);
PSI_API int **init_int_matrix(int rows, int cols);
PSI_API void free_int_matrix(int **array);
PSI_API void zero_int_matrix(int **array, int rows, int cols);
PSI_API void print_int_mat(int **a, int m, int n, std::string out);

/* Functions in long_int_array.c */
PSI_API long int *init_long_int_array(int size);
PSI_API void zero_long_int_array(long int *a, int size);
PSI_API long int **init_long_int_matrix(int rows, int cols);
PSI_API void free_long_int_matrix(long int **array);
PSI_API void zero_long_int_matrix(long int **array, int rows, int cols);
PSI_API void print_long_int_mat(long int **a, int m, int n, std::string out);
PSI_API size_t *init_size_t_array(int size);

/* Functions in block_matrix.c */
[[nodiscard]] PSI_API double **block_matrix(size_t n, size_t m, bool mlock = false);
PSI_API void free_block(double **array);

/* Functions in fndcor */
PSI_API void fndcor(long int *maxcrb, std::string out_fname);
}  // namespace psi

#endif /* header guard */
