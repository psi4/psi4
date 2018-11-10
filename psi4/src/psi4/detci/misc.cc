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

/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here
*/
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <cstdio>

namespace psi {
namespace detci {

/*
** TRANSP_SIGMA(): Function adds the transpose (times a phase factor) of
**    a matrix to itself.
**
*/
void transp_sigma(double **a, int rows, int cols, int phase) {
    int i, j;

    if (rows != cols) {
        outfile->Printf("(transp_sigma): Error, rows != cols\n");
        outfile->Printf("\trows = %d, cols = %d\n", rows, cols);
        return;
    }

    /* do lower triangle */
    if (phase == 1) {
        for (i = 0; i < rows; i++) {
            for (j = 0; j <= i; j++) {
                a[i][j] += a[j][i];
            }
        }
    } else if (phase == -1) {
        for (i = 0; i < rows; i++) {
            for (j = 0; j <= i; j++) {
                a[i][j] -= a[j][i];
            }
        }
    }

    /* fix upper triangle (could remove me later if don't use upper tri) */
    if (phase == 1) {
        for (i = 0; i < rows; i++) {
            for (j = i; j < cols; j++) {
                a[i][j] = a[j][i];
            }
        }
    } else {
        for (i = 0; i < rows; i++) {
            for (j = i; j < cols; j++) {
                a[i][j] = -a[j][i];
            }
        }
    }
}

/*
** SET_ROW_PTRS()
**
** This function sets the row pointers for a 2D matrix allocated as one
** contiguous block of memory
**
*/
void set_row_ptrs(int rows, int cols, double **matrix) {
    int i;
    double *ptr;

    ptr = matrix[0];

    for (i = 1; i < rows; i++) {
        matrix[i] = matrix[0] + i * cols;
    }
}

}  // namespace detci
}  // namespace psi
