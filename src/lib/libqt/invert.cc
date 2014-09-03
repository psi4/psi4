/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*!
  \file
  \brief Invert a small matrix
  \ingroup QT
*/

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include "qt.h"
#include "libparallel/ParallelPrinter.h"
namespace psi {

#define SMALL_DET 1.0E-10

/*!
** INVERT_MATRIX(): The function takes the inverse of a matrix using the
**    C routines in Numerical Recipes in C.
**
** Matt Leininger, Summer 1994
**
** Parameters:
**    \param a       = matrix to take the inverse of
**                     (is modified by invert_matrix())
**    \param y       = the inverse matrix
**    \param N       = the size of the matrices
**    \param outfile = file for error messages
**
** Other variables:
**    col and indx are temporary arrays
**    d is 1 or -1
**
** Returns: double (determinant)
** Note: The original matrix is modified by invert_matrix()
** \ingroup QT
*/
double invert_matrix(double **a, double **y, int N, std::string out)
{
   boost::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            boost::shared_ptr<OutFile>(new OutFile(out)));
   double  d, *col, *colptr;
   register int i, j;
   int *indx ;

   col = init_array(N) ;
   indx = init_int_array(N) ;

   ludcmp(a,N,indx,&d) ;
   for (j=0; j<N; j++) d *= a[j][j];

   /* outfile->Printf("detH0 in invert = %lf\n", fabs(d));
    */

    if (fabs(d) < SMALL_DET) {
      printer->Printf("Warning (invert_matrix): Determinant is %g\n", d);
      printf("Warning (invert_matrix): Determinant is %g\n", d);
      
      }

   for (j=0; j<N; j++) {
       bzero(col,sizeof(double)*N);
       col[j] = 1.0 ;
       lubksb(a,N,indx,col) ;
       colptr = col;
       for (i=0; i<N; i++) y[i][j] = *colptr++;
       }

   free(col);
   free(indx);

   d = fabs(d);
   return(d);
}

}

