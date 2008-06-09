/*!
  \file
  \brief Sort eigenvalues and eigenvectors into ascending order
  \ingroup QT
*/

namespace psi {
	
/* sort(): Sorts eigenvalues and corresponding eigenvectors into
** ascending order.  Based on eigsort.c in libciomr.
**
** TDC, July 2002
**
** \param A  = array of eigenvalues
** \param B  = matrix of eigenvectors
** \param n  = length of A
**
** Returns: none
** \ingroup QT
*/

void sort(double *A, double **B, int n)
{
  int i, j, k;
  double val;

  for (i=0; i < n-1; i++) {
    val = A[k=i];

    for (j=i+1; j < n; j++) 
      if(A[j] <= val) val = A[k=j];

    if (k != i) {
      A[k] = A[i];
      A[i] = val;
      for (j=0;j < n; j++) {
	val = B[j][i];
	B[j][i] = B[j][k];
	B[j][k] = val;
      }

    }
  }
}

/*!
** sort_vector(): Sort the elements of a vector into increasing order
**
** \param A = array to sort
** \param n = length of array
**
** Returns: none
** \ingroup QT
*/
void sort_vector(double *A, int n)
{
  int i, j, k;
  double val;

  for(i=0; i < n-1; i++) {
    val = A[k=i];

    for(j=i+1; j < n; j++)
      if(A[j] <= val) val = A[k=j];

    if(k != i) {
      A[k] = A[i];
      A[i] = val;
    }
  }
}

}

