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

#include <cstdio>
#include <cstdlib>

/*!
  \file
  \brief Routines for 3d arrays
  \ingroup QT
*/

namespace psi {
	
/*!
** 3d_array(): Initialize a (non-contiguous) 3D array
**
** \param p = size of first dimension
** \param q = size of second dimension
** \param r = size of third dimension
**
** Returns: triple-star pointer to 3D array
**
** \ingroup QT
*/

double ***init_3d_array(int p, int q, int r)
{
  double ***A;
  int i,j,k;

  A = (double ***) malloc(p * sizeof(double **));
  for(i=0; i < p; i++) {
    A[i] = (double **) malloc(q * sizeof(double *));
    for(j=0; j < q; j++) {
      A[i][j] = (double *) malloc(r * sizeof(double));
      for(k=0; k < r; k++) {
        A[i][j][k] = 0.0;
      }
    }
  }

  return A;
}


/*!
** free_3d_array(): Free a (non-contiguous) 3D array
**
** \param A = triple-star pointer to the 3D array
** \param p = size of first dimension
** \param q = size of second dimension
**
** Returns: none
**
** \ingroup QT
*/
void free_3d_array(double ***A, int p, int q)
{
  int i,j;

  for(i=0; i < p; i++)
    for(j=0; j < q; j++)
      free(A[i][j]);

  for(i=0; i < p; i++) free(A[i]);

  free(A);
}

}
