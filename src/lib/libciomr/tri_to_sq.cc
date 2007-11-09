/*!
   \file tri_to_sq.cc
   \ingroup (CIOMR)
*/

#include "includes.h"
#include <libciomr/libciomr.h>

extern "C" {

/*!
** tri_to_sq: converts lower triangle to square matrix.
**
** \param amat = lower triangle matrix
** \param bmat = square matrix
** \param size = number of rows/cols of matrix
** \ingroup (CIOMR)
*/
void tri_to_sq(double *amat, double **bmat, int size)
{
  int i, j, ij;

  ij=0;
  for(i = 0 ; i < size ; i++) {
    for(j = 0 ; j <= i ; j++) {
      bmat[i][j] = amat[ij];
      bmat[j][i] = amat[ij++];
    }
  }
}

} /* extern "C" */
