/*!
** \file
** \brief Normalize a set of vectors
** \ingroup QT
*/

#include <cstdio>
#include <cmath>
#include <libciomr/libciomr.h>

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
  register int i, j;

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
FILE *outfile ;
double **mat ;
void normalize(double **A, int rows, int cols) ;

   mat = init_matrix(3, 3) ;
   mat[0][0] = 1.0 ; mat[0][1] = 0.0 ; mat[0][2] = 1.0 ;
   mat[1][0] = 1.0 ; mat[1][1] = -1.0 ; mat[1][2] = 0.0 ;
   mat[2][0] = 0.0 ; mat[2][1] = 0.0 ; mat[2][2] = 0.5 ;

   ffile(&outfile, "output.dat", 0) ;
   fprintf(outfile, "Matrix before normalization process\n") ;
   print_mat(mat,3,3,outfile) ;
   normalize(mat,3,3) ;
   fprintf(outfile, "\nMatrix after normalization process\n") ;
   print_mat(mat,3,3,outfile) ;

   free_matrix(mat,3) ;
} 
#endif

}

