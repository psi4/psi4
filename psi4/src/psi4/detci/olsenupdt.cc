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

/*
** OLSENUPDT.C
** 
** Contains some code related to the Olsen iterative scheme
**
** David Sherrill
** Center for Computational Quantum Chemistry
** February 1996
*/

#include <cstdio>
#include <cmath>
#include "psi4/detci/structs.h"
#include "psi4/detci/ci_tol.h"

namespace psi { namespace detci {

/*
** buf_xy1()
** 
** Do some of the work to get x and y
**
**    x = C^(i) * (Hd - E)^-1 * C^(i)
**    y = C^(i) * (Hd - E)^-1 * sigma^(i)
**
** Specifically this function replaces the buffer 
** (not the vector on disk) of Hd vector
** with the vector C^(i) * (Hd - E)^-1 for use in 
** the evaluation of x and y as described above 
** while tx returns the result of C^(i) * the new
** Hd vector (i.e. x from above)
**        
*/
double buf_xy1(double *c, double *hd, double E, int len)
{
   int i;
   double ci, tval1, tval2;
   double tx = 0.0;

   for (i=0; i<len; i++) {
      ci = c[i];
      tval1 = hd[i] - E;
      if (fabs(tval1) < HD_MIN) tval1 = HD_MIN; /* prevent /0 */
      tval2 = ci / tval1;
      hd[i] = tval2;
      tx += ci * tval2;
      }

   return(tx);
}


/*
** buf_ols_denom()
**
** Get the denominator for the Olsen update
**
*/
void buf_ols_denom(double *a, double *hd, double E, int len)
{
   int i;
   double tval;

   for (i=0; i<len; i++) {
      tval = hd[i] - E;
      if (fabs(tval) < HD_MIN) tval = HD_MIN; /* prevent /0 */
      a[i] /= tval;
      } 
}


/*
** buf_ols_updt()
**
** Do the Olsen update for a buffer
**
*/
void buf_ols_updt(double *a, double *c, double *norm, double *ovrlap,  
      double *tmpnorm, int len)
{
   int i;
   double tval1, tval2, nx = 0.0, ox = 0.0, c1norm = 0.0;

   for (i=0; i<len; i++) {
      tval1 = c[i];
      tval2 = tval1 + a[i];
     /*
      outfile->Printf("C_0[%d] = %14.12lf " \
       "C_1[%d] = %14.12lf C_new[%d] = %14.12lf\n",i,c[i],i,a[i],i,tval2); 
     */
      c[i] = tval2;
      nx += tval2 * tval2;
      ox += tval2 * tval1;
      c1norm += a[i] * a[i];
      }      

   *norm = nx;
   *ovrlap = ox;
   *tmpnorm = c1norm;
}

}} // namespace psi::detci
