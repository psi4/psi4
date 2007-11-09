/*!
  \file invert.c
  \ingroup (QT)
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include "qt.h"

extern "C" {
	
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
** \ingroup (QT)
*/
double invert_matrix(double **a, double **y, int N, FILE *outfile)
{
   double  d, *col, *colptr;
   register int i, j;
   int *indx ;

   col = init_array(N) ;
   indx = init_int_array(N) ;

   ludcmp(a,N,indx,&d) ;
   for (j=0; j<N; j++) d *= a[j][j];
  
   /* fprintf(outfile,"detH0 in invert = %lf\n", fabs(d));
   fflush(outfile); */

    if (fabs(d) < SMALL_DET) {
      fprintf(outfile,"Warning (invert_matrix): Determinant is %g\n", d);
      printf("Warning (invert_matrix): Determinant is %g\n", d);
      fflush(outfile);
      }

   for (j=0; j<N; j++) {
       bzero(col,sizeof(double)*N);
       col[j] = 1.0 ;
       lubksb(a,N,indx,col) ;
       colptr = col;
       for (i=0; i<N; i++) y[i][j] = *colptr++;
       }

   d = fabs(d);
   return(d);
}

} /* extern "C" */