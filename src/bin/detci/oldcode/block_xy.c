#include <cstdio>
#include <cmath>

#define HD_MIN 1.0E-4                     /* for H(diag) calculations */

/*
** block_xy()
** This function evaluates a block's contribution to the quantities x and y
**    x = C^(i) * (Hd - E)^-1 * C^(i)
**    y = C^(i) * (Hd - E)^-1 * sigma^(i)
**
** Parameters:
**   C   =  reference of current iteration's ci vector
**   S   =  reference of current iteration's sigma vector
**   Hd  =  reference of vector of diagonal elements of H
**   E   =  current iteration's energy
**   x   =  pointer to double to hold x
**   y   =  pointer to doulbe to hold y 
**
** Returns: none
*/
void block_xy(double *C, double *S, double *Hd, int len, double E, 
      double *x, double *y)
{
   int i;
   double ci, tval;
   double tx = 0.0, ty = 0.0;

   for (i=0; i<len; i++) {
      ci = C[i];
      tval = Hd[i] - E;
      if (fabs(tval) < HD_MIN) tval = HD_MIN; /* prevent /0 */
      tval = 1.0 / tval;
      tx += ci * ci * tval;
      ty += ci * S[i] * tval;
      }

   *x = tx;
   *y = ty;
}
      



