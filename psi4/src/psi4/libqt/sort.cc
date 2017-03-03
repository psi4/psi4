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
  \file
  \brief Sort eigenvalues and eigenvectors into ascending order
  \ingroup QT
*/

namespace psi {
	
/* sort(): Sorts eigenvalues and corresponding eigenvectors into
** ascending order.  Based on eigsort.c in libciomr.
**
** TDC, July 2002
**
** \param A  = array of eigenvalues
** \param B  = matrix of eigenvectors
** \param n  = length of A
**
** Returns: none
** \ingroup QT
*/

void sort(double *A, double **B, int n)
{
  int i, j, k;
  double val;

  for (i=0; i < n-1; i++) {
    val = A[k=i];

    for (j=i+1; j < n; j++) 
      if(A[j] <= val) val = A[k=j];

    if (k != i) {
      A[k] = A[i];
      A[i] = val;
      for (j=0;j < n; j++) {
	val = B[j][i];
	B[j][i] = B[j][k];
	B[j][k] = val;
      }

    }
  }
}

/*!
** sort_vector(): Sort the elements of a vector into increasing order
**
** \param A = array to sort
** \param n = length of array
**
** Returns: none
** \ingroup QT
*/
void sort_vector(double *A, int n)
{
  int i, j, k;
  double val;

  for(i=0; i < n-1; i++) {
    val = A[k=i];

    for(j=i+1; j < n; j++)
      if(A[j] <= val) val = A[k=j];

    if(k != i) {
      A[k] = A[i];
      A[i] = val;
    }
  }
}

}
