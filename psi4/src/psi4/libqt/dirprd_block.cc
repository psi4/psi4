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
  \brief Take a direct product of two matrices
  \ingroup QT
*/

namespace psi {
	
/*!
 
   dirprd_block()
 
   This function takes two block matrices A and B and multiplies
   each element of B by the corresponding element of A
 
   \param A     = block matrix A
   \param B     = block matrix B 
   \param nrows = number of rows of A and B
   \param ncols = number of columns of A and B
 
   Returns: none

   \ingroup QT
*/
void dirprd_block(double **A, double **B, int rows, int cols)
{
  long int i;
  double *a, *b;
  long size;

  size = ((long) rows) * ((long) cols);

  if(!size) return;

  a = A[0]; b= B[0];

  for(i=0; i < size; i++, a++, b++) (*b) = (*a) * (*b);
}

}
