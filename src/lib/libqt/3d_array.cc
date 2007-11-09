#include <stdio.h>
#include <stdlib.h>

/*!
  \file 3d_array.c
  \ingroup (QT)
*/

extern "C" {
	
/*!
** 3d_array()
** Initialize and free 3d arrays  
** \ingroup (QT)
*/

double ***init_3d_array(int p, int q, int r)
{
  double ***A;
  int i,j,k;

  A = (double ***) malloc(p * sizeof(double **));
  for(i=0; i < p; i++) {
    A[i] = (double **) malloc(q * sizeof(double *));
    for(j=0; j < q; j++) {
      A[i][j] = (double *) malloc(r * sizeof(double));
      for(k=0; k < r; k++) {
        A[i][j][k] = 0.0;
      }
    }
  }

  return A;
}

void free_3d_array(double ***A, int p, int q)
{
  int i,j;

  for(i=0; i < p; i++)
    for(j=0; j < q; j++)
      free(A[i][j]);

  for(i=0; i < p; i++) free(A[i]);

  free(A);
}

} /* extern "C" */