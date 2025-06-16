/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

/*!
** \file
** \brief Converts lower triangle to square matrix
** \ingroup CIOMR
*/

namespace psi {

/*!
** tri_to_sq(): converts lower triangle to square matrix
**
** \param amat = lower triangle matrix
** \param bmat = square matrix
** \param size = number of rows/cols of matrix
**
** Returns: none
**
** \ingroup CIOMR
*/
void tri_to_sq(double *amat, double **bmat, int size) {
    int i, j, ij;

    ij = 0;
    for (i = 0; i < size; i++) {
        for (j = 0; j <= i; j++) {
            bmat[i][j] = amat[ij];
            bmat[j][i] = amat[ij++];
        }
    }
}
}
