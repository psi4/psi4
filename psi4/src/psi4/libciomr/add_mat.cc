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
** \brief Add two matrices
** \ingroup CIOMR
*/

namespace psi {

/*!
** add_mat(): Add matrices a and b into c for n rows and m columns
**
** \param a = double star pointer to first matrix to add
** \param b = double star pointer to second matrix to add
** \param c = double star pointer to matrix to hold the result of a+b
** \param n = number of rows in a,b,c
** \param m = number of columns in a,b,c
**
** \ingroup CIOMR
*/
void add_mat(double **a, double **b, double **c, int n, int m)
{
  int i,j;

  if (n != m) {
    for (i=0; i < n ; i++) {
      for (j=0; j < m ; j++) {
        c[i][j] = a[i][j]+b[i][j];
      }
    }
  }
  else {
    for (i=0; i < n; i++) {
      for (j=0; j < i; j++) {
        c[i][j] = a[i][j]+b[i][j];
        c[j][i] = a[j][i]+b[j][i];
      }
      c[i][i] = a[i][i]+b[i][i];
    }
  }
}

}
