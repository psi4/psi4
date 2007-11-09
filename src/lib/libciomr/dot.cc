/*!
 \file dot.cc
 \ingroup (CIOMR)
*/

#define ALLOC_GLOBALS
#include "includes.h"
#undef ALLOC_GLOBALS
#include "common.h"

extern "C" {

/*!
** dot_mat():
** Takes the dot product between two 2D matrices a and b with dimensions
** n x n and returns the value
** \ingroup (CIOMR)
*/
void dot_mat(double **a, double **b, int n, double *value)
   {
      register int i,j;
      double *ta, *tb, tval;

      tval = 0.0;
      for (i=0; i < n; i++) {
         ta = a[i];
         tb = b[i];
         for (j=0; j < n; j++,ta++,tb++) {
            tval += (*ta) * (*tb);
            }
         }
      *value = tval;
      }

/*!
** dot_arr():
** Take the dot product of the first n elements of two arrays a and b
** and put the result in variable value.
** \ingroup (CIOMR)
*/
void dot_arr(double *a, double *b, int n, double *value)
   {
      register int i;
      double tval;

      tval = 0.0;
      for (i=0; i < n; i++) {
         tval += a[i]*b[i];
         }
      *value = tval;
      }

} /* extern "C" */
