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
** \brief Normalize a set of vectors
** \ingroup QT
*/

#include <cstdio>
#include <cmath>
#include "psi4/libciomr/libciomr.h"

namespace psi {

/* #define STANDALONE */

/*!
** normalize(): Normalize a set of vectors
**
** Assume we're normalizing the ROWS
**
** \param A    = matrix holding vectors to normalize
** \param rows = number of rows in A
** \param cols = number of columns in A
**
** Returns: none
**
** David Sherrill, Feb 1994
** \ingroup QT
*/

void normalize(double **A, int rows, int cols)
{
  double normval;
  int i, j;

  /* divide each row by the square root of its norm */
  for (i=0; i<rows; i++) {
    dot_arr(A[i], A[i], cols, &normval);
    normval = sqrt(normval);
    for (j=0; j<cols; j++) A[i][j] /= normval;
  }

}

#ifdef STANDALONE
main()
{
std::string OutFileRMR ;
double **mat ;
void normalize(double **A, int rows, int cols) ;

   mat = init_matrix(3, 3) ;
   mat[0][0] = 1.0 ; mat[0][1] = 0.0 ; mat[0][2] = 1.0 ;
   mat[1][0] = 1.0 ; mat[1][1] = -1.0 ; mat[1][2] = 0.0 ;
   mat[2][0] = 0.0 ; mat[2][1] = 0.0 ; mat[2][2] = 0.5 ;

   ffile(&outfile, "output.dat", 0) ;
   outfile->Printf( "Matrix before normalization process\n") ;
   print_mat(mat,3,3,outfile) ;
   normalize(mat,3,3) ;
   outfile->Printf( "\nMatrix after normalization process\n") ;
   print_mat(mat,3,3,outfile) ;

   free_matrix(mat,3) ;
}
#endif

}
