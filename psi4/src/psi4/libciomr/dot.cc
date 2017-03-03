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
** \brief Take dot product between two matrices   
** \ingroup CIOMR
*/

namespace psi {

/*!
** dot_mat(): Takes the dot product between two 2D matrices a and b 
** with dimensions n x n and returns the value
**
** \param a     = first matrix for dot product
** \param b     = second matrix for dot product
** \param n     = number of rows/columns for matrices a and b
**
** Returns: value of dot product
**
** \ingroup CIOMR
*/
double dot_mat(double **a, double **b, int n)
{
  int i,j;
  double *ta, *tb, tval;

  tval = 0.0;
  for (i=0; i < n; i++) {
    ta = a[i];
    tb = b[i];
    for (j=0; j < n; j++,ta++,tb++) {
      tval += (*ta) * (*tb);
    }
  }
  return(tval);
}


/*!
** dot_arr(): Take the dot product of the first n elements of two arrays 
** a and b and return the result
**
** \param a = first array to take dot product of
** \param b = second array to take dot product of
** \param n = number of elements in array
** \param value = pointer to hold dot product result
**
** Returns: none
**
** \ingroup CIOMR
*/
void dot_arr(double *a, double *b, int n, double *value)
{
  int i;
  double tval;

  tval = 0.0;
  for (i=0; i < n; i++) {
    tval += a[i]*b[i];
  }
  *value = tval;
}

}
