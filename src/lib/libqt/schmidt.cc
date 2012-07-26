/*!
  \file
  \brief Gram-Schmidt orthogonalize a set of vectors
  \ingroup QT
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>

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
void schmidt(double **A, int rows, int cols, FILE * /*outfile*/)
{
   double *tmp;
   double normval, dotval;
   int i, j;
   register int I;

   /* initialize working array */
   tmp = init_array(cols);

   /* always take the first vector (normalized) as given */
   dot_arr(A[0], A[0], cols, &normval) ; /* normval = dot (A0 * A0) */
   normval = sqrt(normval) ;

   for (i=0; i<cols; i++) A[0][i] /= normval;

   /* now, one at a time, get the new rows */
   for (i=1; i<rows; i++) {
      for (I=0; I<cols; I++) tmp[I] = A[i][I] ;
      for (j=0; j<i; j++) {
         dot_arr(A[i], A[j], cols, &dotval) ;
         for (I=0; I<cols; I++) tmp[I] -= dotval * A[j][I];
         }
      dot_arr(tmp, tmp, cols, &normval);
      normval = sqrt(normval);
      /* fprintf(outfile,"\n norm[%d] = %20.15f\n",i, (1.0/normval));
      fflush(outfile); */
      for (I=0; I<cols; I++) A[i][I] = tmp[I] / normval;
      }

   free(tmp);

}



#ifdef STANDALONE
main()
{
   FILE *outfile ;
   double **mat, **mat_copy, **mat_x_mat ;
   void schmidt(double **A, int rows, int cols) ;

   mat = init_matrix(3, 3) ;
   mat_copy = init_matrix(3, 3) ;
   mat_x_mat = init_matrix(3, 3) ;
   mat[0][0] = 1.0 ; mat[0][1] = 0.0 ; mat[0][2] = 1.0 ;
   mat[1][0] = 1.0 ; mat[1][1] = -1.0 ; mat[1][2] = 0.0 ;
   mat[2][0] = 0.0 ; mat[2][1] = 0.0 ; mat[2][2] = 0.5 ;

   ffile(&outfile, "output.dat", 0) ;
   fprintf(outfile, "Matrix before Gram-Schmidt process\n") ;
   print_mat(mat,3,3,outfile) ;
   schmidt(mat,3,3) ;
   fprintf(outfile, "\nMatrix after Gram-Schmidt process\n") ;
   print_mat(mat,3,3,outfile) ;

   fprintf(outfile, "\nTest A * A = \n") ;

   mmult(mat, 0, mat, 1, mat_x_mat, 0, 3, 3, 3, 0) ;
   print_mat(mat_x_mat,3,3,outfile) ;

   free_matrix(mat,3) ;
   free_matrix(mat_copy,3) ;
   free_matrix(mat_x_mat,3) ;
}
#endif

}

