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
  \brief Gram-Schmidt orthogonalize vector and add to set
  \ingroup QT
*/

#include <cstdio>
#include <cmath>
#include "psi4/libciomr/libciomr.h"

namespace psi {

#define NORM_TOL 1.0E-5

/*!
** SCHMIDT_ADD(): Assume A is a orthogonal matrix.  This function Gram-Schmidt
** orthogonalizes a new vector v and adds it to matrix A.  A must contain
** a free row pointer for a new row.  Don't add orthogonalized v' if
** norm(v') < NORM_TOL.
**
** David Sherrill, Feb 1994
**
** \param A    = matrix to add new vector to
** \param rows = current number of rows in A
**               (A must have ptr for 'rows+1' row.)
** \param cols = columns in A
** \parm v     = vector to add to A after it has been made orthogonal
**               to rest of A
**
** Returns: 1 if a vector is added to A, 0 otherwise
** \ingroup QT
*/
int schmidt_add(double **A, int rows, int cols, double *v)
{
   double dotval, normval ;
   int i, I ;

   for (i=0; i<rows; i++) {
      dot_arr(A[i], v, cols, &dotval) ;
      for (I=0; I<cols; I++) v[I] -= dotval * A[i][I] ;
      }

   dot_arr(v, v, cols, &normval) ;
   normval = sqrt(normval) ;

   if (normval < NORM_TOL)
      return(0) ;
   else {
      if (A[rows] == NULL) A[rows] = init_array(cols) ;
      for (I=0; I<cols; I++) A[rows][I] = v[I] / normval ;
      return(1) ;
      }
}

}
