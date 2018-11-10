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

/*!
  \file
  \brief Fill a symmetric matrix from a lower triangle
  \ingroup QT
*/

namespace psi {

/*!
** fill_sym_matrix(): Fills a symmetric matrix by placing the elements of
** the lower triangle into the upper triangle.
**
** \param  A    = matrix to symmetrize
** \param  size = number of rows or columns (assume square)
**
** Returns: none
** \ingroup QT
*/
void fill_sym_matrix(double **A, int size) {
    double **row, *col;
    int rc, cc;

    row = A;
    for (rc = 0; rc < (size - 1); rc++) {
        col = *row;
        for (cc = 0; cc < size; cc++) {
            if (cc > rc) {
                *col = A[cc][rc];
            }
            col++;
        }
        row++;
    }
}

}  // namespace psi
