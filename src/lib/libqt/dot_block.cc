/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

/*!
  \file
  \brief Take dot product of two block matrices
  \ingroup QT
*/

#include <libqt/qt.h>

namespace psi {
	
/*!
** dot_block(): Find dot product of two block matrices
**
** \param A     = block matrix A
** \param B     = block matrix B 
** \param nrows = number of rows of A and B
** \param ncols = number of columns of A and B
** \param alpha = scale factor by which the dot product is multiplied
**
** Returns: dot product
** \ingroup QT
*/
double dot_block(double **A, double **B, int rows, int cols, double alpha)
{
    double value;
    long int size;
        
    size = ((long) rows) * ((long) cols);
        
    if(!size) return 0.0;
        
    C_DGEMM('T', 'N', 1, 1, size, alpha, A[0], 1, B[0], 1, 0.0, &value, 1);
        
    return value;
}

}