/*!
** \file
** \brief Convert square matrix to lower triangle packing
** \ingroup CIOMR
*/

namespace psi {

/*!
** sq_to_tri(): converts square matrix to lower triangle
**
** \param bmat = matrix to convert
** \param amat = array to put lower triangle of bmat into
** \param size = number of rows/columns of bmat
**
** Returns: none
**
** \ingroup CIOMR
*/ 
void sq_to_tri(double **bmat, double *amat, int size)
{
  int i, j, ij;

  ij=0;
  for(i=0; i < size; i++) {
    for(j=0 ; j <= i; j++) {
      amat[ij++] = bmat[i][j];
    }
  }

}

}

