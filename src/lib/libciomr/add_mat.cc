/*!
** \file add_mat.cc
** \ingroup (CIOMR)
*/

#include "includes.h"

extern "C" {

/*!
** add_mat(): Add matrices a and b into c for n rows and m columns
**
** \ingroup (CIOMR)
*/
void add_mat(double **a, double **b, double **c, int n, int m)
   {
      register int i,j;

      if (n != m) {
         for (i=0; i < n ; i++) {
            for (j=0; j < m ; j++) {
               c[i][j] = a[i][j]+b[i][j];
               }
            }
         }
      else {
         for (i=0; i < n; i++) {
            for (j=0; j < i; j++) {
               c[i][j] = a[i][j]+b[i][j];
               c[j][i] = a[j][i]+b[j][i];
               }
            c[i][i] = a[i][i]+b[i][i];
            }
         }
      }

} /* extern "C" */
