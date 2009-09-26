/*! \file
    \ingroup INPUT
    \brief Enter brief description of file here 
*/

#ifndef _psi3_bin_input_linalg_h_
#define _psi3_bin_input_linalg_h_

#include "float.h"
#include <cstdio>

namespace psi { namespace input {

FLOAT** create_matrix(int a, int b);
void delete_matrix(FLOAT** M);
void print_matrix(FLOAT** a, int m, int n, FILE* out);
FLOAT** convert_matrix(double **M, int a, int b, int transpose);
int matrix_mult(FLOAT** A, int arow, int acol, FLOAT** B, int brow, int bcol, FLOAT** C);
void lu_decom(FLOAT** a, int n, int* indx, FLOAT* d);

}} // namespace psi::input

#endif
