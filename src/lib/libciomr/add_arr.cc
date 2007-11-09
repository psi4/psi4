/*! \defgroup CIOMR libciomr: The Old PSI I/O and Math Library */

/*!
  \file add_arr.cc
  \ingroup (CIOMR)
*/

extern "C" {
/*!
** add_arr: Add arrays a and b and put the result in array c.  Adds
** the first n elements
** \ingroup (CIOMR)
*/
void add_arr(double *a, double *b, double *c, int n)
   {
      register int i;

      for (i=0; i < n ; i++) {
         c[i] = a[i]+b[i];
         }
      }
} /* extern "C" */