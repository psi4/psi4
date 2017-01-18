/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here 
*/

#include <cstdio>
#include <cmath>
#include "psi4/detci/ci_tol.h"
#include "psi4/detci/structs.h"

namespace psi { namespace detci {

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
   int i;
   double norm = 0.0, tval, tval2;

   for (i=0; i<size; i++) {
      tval = lambda - Hd[i];
      if (precon == PRECON_LANCZOS) tval = 1.0;

      if (fabs(tval) > HD_MIN){ 
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
   int i;
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
