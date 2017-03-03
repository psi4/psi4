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

/*! \defgroup CIOMR libciomr: The PSI I/O and Math Library */

/*!
  \file
  \brief Add two arrays
  \ingroup CIOMR
*/

namespace psi {
/*!
** add_arr(): Add arrays a and b and put the result in array c.  Adds
** the first n elements
**
** \param a = first array to add
** \param b = second array to add
** \param c = array to hold the result of a+b
** \param n = number of elements to add
**
** Returns: none
**
** \ingroup CIOMR
*/
void add_arr(double *a, double *b, double *c, int n)
{
  int i;

  for (i=0; i < n; i++) {
    c[i] = a[i]+b[i];
  }
}

}
