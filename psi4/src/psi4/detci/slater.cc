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
    \brief Some of the Slater determinant routines
*/

#include <cstdio>

namespace psi { namespace detci {

/*
** calc_orb_diff()
** 
** Function calculates the differences in the inputed alpha(beta) strings
** and returns whether this difference is in 0, 1, 2, or more spin orbitals
**
** extended  = 0 if strings which differ by > 2 orbitals are not considered
**             and 1 if the I and J diff array should be setup for > 2
**             orbital differences.
**
*/
int calc_orb_diff(int cnt, unsigned char *I, unsigned char *J, 
      int *I_alpha_diff, int *J_alpha_diff, int *sign, int *same, 
      int extended)
{
   int i = 0; 
   int j = 0; 
   int k = 0; 
   int icnt = 0;  /* number of counts in I array */
   int jcnt = 0;  /* number of counts in J array */ 
   int flipI = 0; /* number of flips in I string to have max coincidence */ 
   int flipJ = 0; /* number of flips in J string to have max coincidence */ 
   int ndoI = 0; /* number of different orbitals in I string */ 
   int ndoJ = 0; /* number of different orbitals in J string */

   while ((i < cnt) && (j < cnt)) {

      if (I[i] == J[j]) {
         same[k] = (int) I[i]; 
         i++;
         j++; 
         k++;
         }

      else if (I[i] < J[j]) {
         I_alpha_diff[icnt++] = (int) I[i];
         flipI += i - ndoI; 
         ndoI++; 
         i++;
         if ( ((icnt + jcnt) > 4) && !extended) return(-1);
         } 

      else if (I[i] > J[j]) {
         J_alpha_diff[jcnt++] = (int) J[j];
         flipJ += j - ndoJ; 
         ndoJ++; 
         j++;
         if ( ((icnt + jcnt) > 4) && !extended) return(-1);
         }

      } /* end while loop */

   /* Matt: this used to be a bit different...is this version
    * actually faster?  Seemed better the other way but I dunno ..CDS */
   if (i != j) {
      if (i<j) {
        if ( ((j-i+icnt)>=3) && !extended) return(-1);
        while (i<cnt) {
          I_alpha_diff[icnt++] = (int) I[i];
          flipI += i - ndoI;
          ndoI++;
          i++;
          } 
       }
      else {
          if ( ((i-j+jcnt)>=3) && !extended) return(-1);
          while (j<cnt) {
            J_alpha_diff[jcnt++] = (int) J[j];
            flipJ += j - ndoJ;
            ndoJ++;
            j++;
            } 
        } 
      } /* end if (i != j) loop */

   *sign += flipI + flipJ; 

   return icnt;
/*
   if (icnt == 2)
      return 2; 
   else if (icnt == 1)
      return 1; 
   else
      return 0; 
*/ 

}



/*
** common_orbs()
**
** This function creates the arrays common_docc, common_alpha_socc,
** and common_beta_socc used in the Slater Determinant Class.
**
*/
void common_orbs(int *same_alpha, int *same_beta, int cnt_alpha, int cnt_beta,
            int *common_docc, int *common_alpha_socc, int *common_beta_socc,
            int *cnt_docc, int *cnt_alpha_socc, int *cnt_beta_socc)  
{
   int i = 0; 
   int j = 0; 
   int k = 0; 
 
   while ((i<cnt_alpha) && (j<cnt_beta)) { 

      if (same_alpha[i] == same_beta[j]) {
         common_docc[(*cnt_docc)++] = same_alpha[i]; 
         i++; 
         j++; 
         }
   
      else if (same_alpha[i] < same_beta[j]) {
         common_alpha_socc[(*cnt_alpha_socc)++] = same_alpha[i]; 
         i++; 
         } 

      else if (same_alpha[i] > same_beta[j]) {
         common_beta_socc[(*cnt_beta_socc)++] = same_beta[j]; 
         j++; 
         } 

      } /* end top while loop */

   while (i<cnt_alpha) {
      common_alpha_socc[(*cnt_alpha_socc)++] = same_alpha[i];
      i++; 
      } 

   while (j<cnt_beta) {
      common_beta_socc[(*cnt_beta_socc)++] = same_beta[j]; 
      j++; 
      }

} 

}} // namespace psi::detci
