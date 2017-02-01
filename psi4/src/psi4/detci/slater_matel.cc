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
#include <cstdlib>
#include "psi4/psi4-dec.h"

#include "psi4/detci/slaterd.h"
#include "psi4/detci/ciwave.h"


namespace psi { namespace detci {

extern int calc_orb_diff(int cnt, unsigned char *I, unsigned char *J,
   int *I_alpha_diff, int *J_alpha_diff, int *sign, int *same,
   int extended);
extern void common_orbs(int *same_alpha_, int *same_beta_, int cnt_alpha,
   int cnt_beta, int *common_docc_, int *common_alpha_socc,
   int *common_beta_socc, int *cnt_docc, int *cnt_alpha_socc,
   int *cnt_beta_socc);


double CIWavefunction::matrix_element(SlaterDeterminant* I, SlaterDeterminant* J)
{
   int *common_socc1, *common_socc2;
   int sign=0;
   int nalp, nbet, nalp_same, nbet_same;
   int cnt_docc=0, cnt_alp_socc=0, cnt_bet_socc=0;
   int total_diff=-4, alpha_diff=-2, beta_diff=-2;
   int diffspin, samespin, common1, common2;
   int i, j, k, l, a, b, m;
   double val = 0.0;

   if (I->nalp_ != J->nalp_ || I->nbet_ != J->nbet_) {
      outfile->Printf("(matrix_element): unequal length alp/bet strings!\n");
      return(0.0);
    }

   nalp = I->nalp_;
   nbet = I->nbet_;

   if (sme_first_call_) {
      same_alpha_ = (int *) malloc(sizeof(int) * nalp);
      same_beta_ = (int *) malloc(sizeof(int) * nbet);
      common_alp_socc_ = (int *) malloc(sizeof(int) * nalp);
      common_bet_socc_ = (int *) malloc(sizeof(int) * nbet);
      common_docc_ = (int *) malloc(sizeof(int) * nalp);
      I_diff_ = (int **) malloc(sizeof(int *) * 2);
      I_diff_[0] = (int *) malloc(sizeof(int) * nalp);
      I_diff_[1] = (int *) malloc(sizeof(int) * nalp);
      J_diff_ = (int **) malloc(sizeof(int *) * 2);
      J_diff_[0] = (int *) malloc(sizeof(int) * nalp);
      J_diff_[1] = (int *) malloc(sizeof(int) * nalp);
      init_nalp_ = nalp;
      init_nbet_ = nbet;
      sme_first_call_ = 0;
   }

   // make sure that subsequent calls don't change the number of alpha or
   // beta electrons without re-initializing the static arrays
   else {
      if ((nalp != init_nalp_) || (nbet != init_nbet_)) {
         throw PsiException("(matrix_element): nalp/nbet != init_nalp_/nbet",__FILE__,__LINE__);
      }
   }

   alpha_diff = calc_orb_diff(nalp, I->Occs_[0], J->Occs_[0], I_diff_[0],
      J_diff_[0], &sign, same_alpha_, 0);
   beta_diff = calc_orb_diff(nbet, I->Occs_[1], J->Occs_[1], I_diff_[1],
      J_diff_[1], &sign, same_beta_, 0);

   total_diff = alpha_diff + beta_diff ;
   nalp_same = nalp - alpha_diff ;
   nbet_same = nbet - beta_diff ;

   /*
  outfile->Printf(" alpha_diff = %d\n", alpha_diff) ;
  outfile->Printf(" beta_diff = %d\n", beta_diff) ;
  outfile->Printf(" total_diff = %d\n", total_diff) ;
  outfile->Printf(" sign = %d\n", sign) ;

   for (i=0; i<alpha_diff; i++)
     outfile->Printf(" I_diff_[0][%d] = %d ", i, I_diff_[0][i]+1) ;
     outfile->Printf("\n") ;
   for (i=0; i<alpha_diff; i++)
     outfile->Printf(" J_diff_[0][%d] = %d ", i, J_diff_[0][i]+1) ;
     outfile->Printf("\n\n") ;

   for (i=0; i<beta_diff; i++)
     outfile->Printf(" I_diff_[1][%d] = %d ", i, I_diff_[1][i]+1) ;
     outfile->Printf("\n") ;
   for (i=0; i<beta_diff; i++)
     outfile->Printf(" J_diff_[1][%d] = %d ", i, J_diff_[1][i]+1) ;
     outfile->Printf("\n\n") ;
   for (i=0; i<nalp_same; i++)
     outfile->Printf(" same_alpha_[%d] = %d", i, same_alpha_[i]+1) ;
     outfile->Printf("\n\n") ;
   for (i=0; i<nbet_same; i++)
     outfile->Printf(" same_beta_[%d] = %d", i, same_beta_[i]+1) ;
     outfile->Printf("\n\n") ;
   */

   if ((alpha_diff == -2) || (beta_diff == -2)) {
      outfile->Printf( "(matrix_element): Problem with calc_orb_diff.\n");
      outfile->Printf( "  Returns -2 value. \n");
      }

   else if ((alpha_diff == -1) || (beta_diff == -1) || total_diff > 2) {
      return(0.0);
      }

   else if (total_diff == 2) {

      if (alpha_diff == 1) {     /* Case 1: 1 in alpha and 1 in beta */

         /* assign i,j,k,l for <ij||kl> */
         i = I_diff_[0][0];
         j = I_diff_[1][0];
         k = J_diff_[0][0];
         l = J_diff_[1][0];

         #ifdef PRINT_INTS
         if (sign % 2)outfile->Printf("-");
         #endif

         val = get_twoel(i,k,j,l);
         if (sign % 2) val = -val;

         return(val);
         } /* end Case 1 */

      else if ((alpha_diff == 2) || (beta_diff == 2)) { /* Case 2: 2 in alpha */

         if (alpha_diff == 2)
            diffspin = 0;
         else diffspin = 1;

         /* assign <ij||kl> */
         i = I_diff_[diffspin][0];
         j = I_diff_[diffspin][1];
         k = J_diff_[diffspin][0];
         l = J_diff_[diffspin][1];

         #ifdef PRINT_INTS
         if (sign % 2)outfile->Printf("- [ ");
         #endif

         val = get_twoel(i,k,j,l);
         val -= get_twoel(i,l,j,k);

         #ifdef PRINT_INTS
         if (sign % 2)outfile->Printf(" ] ");
         #endif

         if (sign % 2) val = -val;

         return(val);
         } /* end else if for differ by 2 in alpha or beta */

      else {
         throw PsiException("Error (matrix_element): total_diff != alpha_diff + beta_diff",__FILE__,__LINE__);
         }

      } /* end else if for differ by 2 spin orbitals */


   /* Differ by 1 spin orbital */
   else if (total_diff == 1) {

      common_orbs(same_alpha_, same_beta_, nalp_same, nbet_same, common_docc_,
         common_alp_socc_, common_bet_socc_, &cnt_docc, &cnt_alp_socc,
         &cnt_bet_socc);

      /*
     outfile->Printf("cnt_docc = %d\n", cnt_docc);
     outfile->Printf("cnt_alp_socc = %d\n", cnt_alp_socc);
     outfile->Printf("cnt_bet_socc = %d\n", cnt_bet_socc);

      for (i=0; i<cnt_docc; i++)
        outfile->Printf("common_docc_[%d] = %d\n", i, common_docc_[i]+1) ;
        outfile->Printf("\n") ;
      for (i=0; i<cnt_alp_socc; i++)
        outfile->Printf("common_alp_socc_[%d] = %d\n", i, common_alp_socc_[i]+1) ;
        outfile->Printf("\n") ;
      for (i=0; i<cnt_bet_socc; i++)
        outfile->Printf("common_bet_socc_[%d] = %d\n", i, common_bet_socc_[i]+1) ;
        outfile->Printf("\n") ;
      */

      if (alpha_diff == 1) { /* Case 1: 1 in alpha */
         diffspin = 0;
         samespin = 1;
         common1 = cnt_alp_socc;
         common2 = cnt_bet_socc;
         common_socc1 = common_alp_socc_;
         common_socc2 = common_bet_socc_;
         } /* end Case 1: 1 in alpha */

      else if (beta_diff == 1) { /* Case 1: 1 in beta */
         diffspin = 1;
         samespin = 0;
         common1 = cnt_bet_socc;
         common2 = cnt_alp_socc;
         common_socc1 = common_bet_socc_;
         common_socc2 = common_alp_socc_;
         } /* end Case 1: 1 in beta */

      else {
         throw PsiException("(matrix_element): impossible case abdiff",__FILE__,__LINE__);
         }

      i = I_diff_[diffspin][0];
      j = J_diff_[diffspin][0];

      val = get_onel(i,j);

      #ifdef PRINT_INTS
      if (sign % 2)
        outfile->Printf(" - [ \n");
      #endif

      for (b=0; b<cnt_docc; b++) { /* looping over all common docc electrons. */
         m = common_docc_[b];

         #ifdef PRINT_INTS
        outfile->Printf("2 * ");
         #endif

         val += 2.0 * get_twoel(i,j,m,m);

         #ifdef PRINT_INTS
        outfile->Printf("-");
         #endif

         val -= get_twoel(i,m,m,j);
         } /* end loop over all common docc electrons */

      for (b=0; b<common2; b++) { /* looping over all beta or alpha elec.;
                                    beta if differs in one by alpha */
         m = common_socc2[b];
         val += get_twoel(i,j,m,m);
         }

      for (b=0; b<common1; b++) { /* looping over all alpha or beta elec.;
                                    alpha if differs in one by alpha */
         m = common_socc1[b];

         val += get_twoel(i,j,m,m);

         #ifdef PRINT_INTS
        outfile->Printf("-");
         #endif

         val -= get_twoel(i,m,m,j);
         }

      if (sign % 2) val = -val;

      #ifdef PRINT_INTS
      if (sign % 2)
        outfile->Printf(" ] \n");
      #endif

      return(val);

      } /* end if (total_diff == 1) for Case 1: Differ by 1 spin orbital */


   else if (total_diff == 0) {

      common_orbs(same_alpha_, same_beta_, nalp_same, nbet_same, common_docc_,
         common_alp_socc_, common_bet_socc_, &cnt_docc, &cnt_alp_socc,
         &cnt_bet_socc);

      val = 0.0;

      #ifdef PRINT_INTS
      if (sign % 2)outfile->Printf("- [ \n");
     outfile->Printf(" 2 * [ ");
      #endif

      /* get one-electron integrals */

      for (b=0; b<cnt_docc; b++) {
         m = common_docc_[b];
         val += 2.0 * get_onel(m,m);
         }

      #ifdef PRINT_INTS
     outfile->Printf(" ] ");
      #endif

      for (b=0; b<cnt_alp_socc; b++) {
         m = common_alp_socc_[b];
         val += get_onel(m,m);
         }

      for (b=0; b<cnt_bet_socc; b++) {
         m = common_bet_socc_[b];
         val += get_onel(m,m);
         }

      #ifdef PRINT_INTS
     outfile->Printf("\n");
      #endif

      /* get two-electron integrals */

      for (a=0; a<cnt_docc; a++) {
         i = common_docc_[a];

         for (b=0; b<cnt_alp_socc; b++) {
            j = common_alp_socc_[b];

            #ifdef PRINT_INTS
           outfile->Printf("2 * ");
            #endif

            val += 2.0 * get_twoel(i,i,j,j);

            #ifdef PRINT_INTS
           outfile->Printf("- ");
            #endif

            val -= get_twoel(i,j,i,j);
            }

         for (b=0; b<cnt_bet_socc; b++) {
            j = common_bet_socc_[b];

            #ifdef PRINT_INTS
           outfile->Printf("2 * ");
            #endif

            val += 2.0 * get_twoel(i,i,j,j);

            #ifdef PRINT_INTS
           outfile->Printf("- ");
            #endif

            val -= get_twoel(i,j,i,j);
            }

         val += get_twoel(i,i,i,i);

         for (b=a+1; b<cnt_docc; b++) {

            j = common_docc_[b];

            #ifdef PRINT_INTS
           outfile->Printf("4 * ");
            #endif

            val += 4.0 * get_twoel(i,i,j,j);

            #ifdef PRINT_INTS
           outfile->Printf("- 2 * ");
            #endif

            val -= 2.0 * get_twoel(i,j,j,i);

            } /* end loop over pairs of doccs */

         } /* end loop a over doccs */

      /* loop over unique pairs of alpha electrons and alpha with beta */
      for (a=0; a<cnt_alp_socc; a++) {

         i = common_alp_socc_[a];

         for (b=a+1; b<cnt_alp_socc; b++) {

            j = common_alp_socc_[b];

            val += get_twoel(i,i,j,j);

            #ifdef PRINT_INTS
           outfile->Printf("- ");
            #endif

            val -= get_twoel(i,j,j,i);
            }

         for (b=0; b<cnt_bet_socc; b++) {
            j = common_bet_socc_[b];
            val += get_twoel(i,i,j,j);
            }
         }


      /* loop over unique pairs of beta electrons */
      for (a=0; a<cnt_bet_socc; a++) {
         for (b=a+1; b<cnt_bet_socc; b++) {

            i = common_bet_socc_[a];
            j = common_bet_socc_[b];

            val += get_twoel(i,i,j,j);

            #ifdef PRINT_INTS
           outfile->Printf("- ");
            #endif

            val -= get_twoel(i,j,j,i);
            }
         }

      if (sign % 2) val = -val;

      #ifdef PRINT_INTS
      if (sign % 2)
        outfile->Printf(" ] \n");
      #endif

      return(val);

      } /* end total_diff == 0 case */

   else {
      throw PsiException("(matrix_element): Impossible case for total_diff!",__FILE__,__LINE__);
      }

   return(0.0);
}

//#ifdef STANDALONE
//main() {
//   unsigned char string1[2], string2[2], string3[2];
//
//   SlaterDeterminant A;
//   SlaterDeterminant B;
//   SlaterDeterminant C;
//   double value = 0.0 ;
//
//   string1[0] = (unsigned char) 0;
//   string1[1] = (unsigned char) 1;
//   string2[0] = (unsigned char) 0;
//   string2[1] = (unsigned char) 2;
//   string3[0] = (unsigned char) 0;
//   string3[1] = (unsigned char) 3;
//
//   A.set(2, string1, 2, string1);
//   B.set(2, string1, 2, string3);
//
//  outfile->Printf("Slater determinant A \n");
//   A.print() ;
//  outfile->Printf("Slater determinant B \n");
//   B.print() ;
//  outfile->Printf("Matrix element = %lf\n", matrix_element(&A,&B));
//}

//#ifdef PRINT_INTS
//double get_twoel(int i, int j, int k, int l)
//{
//  outfile->Printf("(%d %d | %d %d) ", i, j, k, l);
//   return(0.0);
//}
//
//double get_onel(int i, int j)
//{
//  outfile->Printf("(%d|h|%d) ", i, j);
//   return(0.0) ;
//}
//#endif

//#endif


}} // namespace psi::detci
