/*!
  \file
  \brief Gram-Schmidt orthogonalize vector and add to set
  \ingroup QT
*/

#include <cstdio>
#include <cmath>
#include <libciomr/libciomr.h>

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

