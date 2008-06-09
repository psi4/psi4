/*!
** \file
** \brief Zero arrays or matrices of doubles
** \ingroup CIOMR
*/

#include <strings.h>

namespace psi {

/*!
** zero_arr(): zero out an array of length 'size'
**
** \param a    = array to zero out
** \param size = how many elements of a to zero
**
** Returns: none
**
** \ingroup CIOMR
*/
void zero_arr(double *a, int size)
{
  bzero(a,sizeof(double)*size);
}

/*!
** zero_mat(): zero out a matrix 'a' with n rows and m columns 
** 
** \param a = matrix of doubles to zero out
** \param n = number of rows in a
** \param m = number of columns in a
**
** \ingroup CIOMR
*/
void zero_mat(double **a, int n, int m)
{
  register int i;

  for (i=0; i < n; i++) {
    bzero(a[i],sizeof(double)*m);
  }
}

}

