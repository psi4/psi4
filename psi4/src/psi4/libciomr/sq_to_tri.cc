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

/*!
** \file
** \brief Convert square matrix to lower triangle packing
** \ingroup CIOMR
*/

namespace psi {

/*!
** sq_to_tri(): converts square matrix to lower triangle
**
** \param bmat = matrix to convert
** \param amat = array to put lower triangle of bmat into
** \param size = number of rows/columns of bmat
**
** Returns: none
**
** \ingroup CIOMR
*/ 
void sq_to_tri(double **bmat, double *amat, int size)
{
  int i, j, ij;

  ij=0;
  for(i=0; i < size; i++) {
    for(j=0 ; j <= i; j++) {
      amat[ij++] = bmat[i][j];
    }
  }

}

}
