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
** \brief Zero arrays or matrices of doubles
** \ingroup CIOMR
*/

#include <strings.h>

namespace psi {

/*!
** zero_arr(): zero out an array of length 'size'
**
** \param a    = array to zero out
** \param size = how many elements of a to zero
**
** Returns: none
**
** \ingroup CIOMR
*/
void zero_arr(double *a, int size)
{
  bzero(a,sizeof(double)*size);
}

/*!
** zero_mat(): zero out a matrix 'a' with n rows and m columns 
** 
** \param a = matrix of doubles to zero out
** \param n = number of rows in a
** \param m = number of columns in a
**
** \ingroup CIOMR
*/
void zero_mat(double **a, int n, int m)
{
  int i;

  for (i=0; i < n; i++) {
    bzero(a[i],sizeof(double)*m);
  }
}

}
