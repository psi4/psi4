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
  \brief Gram-Schmidt orthogonalize a set of vectors
  \ingroup QT
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libciomr/libciomr.h"

namespace psi {

/* #define STANDALONE */

/*!
** SCHMIDT(): Gram-Schmidt orthogonalize a set of vectors
**
** Assume we're orthogonalizing the ROWS, since in C
** a vector is usually a row more often than a column.
**
** David Sherrill, Feb 1994
**
** \param A    = matrix to orthogonalize (matrix of doubles)
** \param rows = rows of A
** \param cols = columns of A
**
** Returns: none
** \ingroup QT
*/
void schmidt(double **A, int rows, int cols, std::string)
{
   double RValue;
   for(size_t i=0;i<cols;++i){
      dot_arr(A[i],A[i],cols,&RValue);
      RValue=sqrt(RValue);
      for(size_t I=0;I<cols;++I)A[i][I]/=RValue;
      for(size_t j=i+1;j<cols;++j){
         dot_arr(A[i],A[j],cols,&RValue);
         for(size_t I=0;I<cols;++I)A[j][I]-=RValue*A[i][I];
      }
   }

}



#ifdef STANDALONE
main()
{
   std::string OutFileRMR ;
   double **mat, **mat_copy, **mat_x_mat ;
   void schmidt(double **A, int rows, int cols) ;

   mat = init_matrix(3, 3) ;
   mat_copy = init_matrix(3, 3) ;
   mat_x_mat = init_matrix(3, 3) ;
   mat[0][0] = 1.0 ; mat[0][1] = 0.0 ; mat[0][2] = 1.0 ;
   mat[1][0] = 1.0 ; mat[1][1] = -1.0 ; mat[1][2] = 0.0 ;
   mat[2][0] = 0.0 ; mat[2][1] = 0.0 ; mat[2][2] = 0.5 ;

   ffile(&outfile, "output.dat", 0) ;
   outfile->Printf( "Matrix before Gram-Schmidt process\n") ;
   print_mat(mat,3,3,outfile) ;
   schmidt(mat,3,3) ;
   outfile->Printf( "\nMatrix after Gram-Schmidt process\n") ;
   print_mat(mat,3,3,outfile) ;

   outfile->Printf( "\nTest A * A = \n") ;

   mmult(mat, 0, mat, 1, mat_x_mat, 0, 3, 3, 3, 0) ;
   print_mat(mat_x_mat,3,3,outfile) ;

   free_matrix(mat,3) ;
   free_matrix(mat_copy,3) ;
   free_matrix(mat_x_mat,3) ;
}
#endif

}
