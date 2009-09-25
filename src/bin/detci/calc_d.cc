/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here 
*/

#define EXTERN
#include <cstdio>
#include <cmath>
#include "ci_tol.h"
#include "structs.h"
#include "globals.h"

namespace psi { namespace detci {

/*
** calc_d()
**
** Function calculates a block of the numerators for the Davidson 
** algorithm correction vector d.
**
** Parameters:
**    target  = array to store result
**    alpha   = coefficient
**    sigma   = sigma block
**    lambda  = energy coefficient
**    c       = c vector block
**    size    = size of block
**
** Returns: none
*/
void calc_d(double *target, double alpha, double *sigma, double lambda,
      double *c, int size)
{
   register int i;
   double tval;

   for (i=0; i<size; i++) {
      tval = alpha * (sigma[i] - lambda * c[i]);
      target[i] += tval;
      }
}


/*
** calc_d2()
**
** Function calculates a block of the denominators for the Davidson 
** algorithm correction vector d.
**
** Parameters:
**    target  = array to store result
**    lambda  = coefficient
**    Hd      = Diagonal Hamiltonian block
**    size    = size of block
**    precon  = type of preconditioner (Parameters.precon)
**
** Returns: sum of squares of coefficients
*/
double calc_d2(double *target, double lambda, double *Hd, int size, int precon)
{
   register int i;
   double norm = 0.0, tval, tval2;

   for (i=0; i<size; i++) {
      tval = lambda - Hd[i];
      if (precon == PRECON_LANCZOS) tval = 1.0;
      if (fabs(tval) > HD_MIN) { 
         tval2 = (target[i] /= tval);
         norm += tval2 * tval2;
         }
      else target[i] = 0.0;
      }

   return(norm);
}

/*
** calc_mpn_vec()
**
** Function calculates a block of the denominators for the kth order
** wavefunction in a perturbation series
**
** Parameters:
**    target  = array to store result
**    energy  = energy
**    Hd      = Diagonal Hamiltonian block
**    size    = size of block
**    sign1   = sign1*E + sign2*Hd 
**    sign2   = sign1*E + sign2*Hd 
**
** Returns: sum of squares of coefficients
*/
double calc_mpn_vec(double *target, double energy, double *Hd, int size, double
        sign1, double sign2, int precon)
{
   register int i;
   double norm = 0.0, tval, tval2;

   for (i=0; i<size; i++) {
      tval = sign1*energy + sign2*Hd[i];
      if (precon==1) 
        tval2 = (target[i] /= tval);
      else if (precon==0) 
        tval2 = (target[i] *= tval);
      norm += tval2 * tval2;
      }
   return(norm);

}

}} // namespace psi::detci

