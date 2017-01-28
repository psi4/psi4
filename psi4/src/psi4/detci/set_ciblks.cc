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
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/detci/structs.h"
#include "psi4/detci/ciwave.h"

namespace psi { namespace detci {

extern int og_lex_addr(struct olsen_graph *Graph, int *occs, int nel,
      int *listnum);
extern void print_ciblk_summary(std::string out);


void CIWavefunction::set_ciblks()
{
   int *occs, *occs2;
   int i, j, k, l, m, irrep, cnt, cnt2, betirrep;
   int betsocc;    /* number of beta singly occupied orbitals accounted for */
   int alp_sym, bet_sym;
   int nalp1, nalp3, alpcode, nbet1, nbet3, betcode;
   int nalp4, nbet4, maxblk, xlvl;
   int nas, nbs;
   int nblocks=0;
   double orbsum = 0.0;
   int set = 0;

   CalcInfo_->num_alp_str = AlphaG_->num_str;
   CalcInfo_->num_bet_str = BetaG_->num_str;
   xlvl = Parameters_->ex_lvl;

   if (print_) {
      outfile->Printf( "    There are %d alpha and %d beta strings\n", CalcInfo_->num_alp_str, CalcInfo_->num_bet_str);
      }

   /* Get the occupations for the reference alpha and beta strings.
    * This used to be one line (occupy orbitals 0 to CalcInfo_->num_alp_expl)
    * but now it's got to be more complicated due to orbital renumbering
    * (now by RAS then by symm then by energy) and due to open-shell
    * cases.  Assume that occupied orbs in reference belong only to RAS I
    * or RAS II.  Should add a check for this sometime.  The following
    * routine is by MLL and CDS 11/96.  Modified to put RAS II occs in
    * separate array because otherwise you get strings like |0 1 3 2>.
    */
   occs = init_int_array(CalcInfo_->num_alp_expl);
   occs2 = init_int_array(CalcInfo_->num_alp_expl);

   betsocc = 0; cnt = 0; cnt2 = 0;
   for (i=0; i<CalcInfo_->nirreps; i++) {
      j = CalcInfo_->docc[i] - CalcInfo_->dropped_docc[i];
      k = CalcInfo_->ras_opi[0][i];
      l = (j<k) ? j : k;
      for (m=0; m<l; m++) /* RAS I part */
         occs[cnt++] = CalcInfo_->ras_orbs[0][i][m];
      for (; m<j; m++)    /* RAS II part */
         occs2[cnt2++] = CalcInfo_->ras_orbs[1][i][m-k];
      if (CalcInfo_->socc[i] != 0) {
         if (Parameters_->opentype == PARM_OPENTYPE_NONE)
            outfile->Printf("Warning: ignoring socc since opentype=none\n");
         else if (Parameters_->opentype == PARM_OPENTYPE_HIGHSPIN) {
            j += CalcInfo_->socc[i];
            l = (j<k) ? j : k;
            for (; m<l; m++) /* RAS I part */
               occs[cnt++] = CalcInfo_->ras_orbs[0][i][m];
            for (; m<j; m++) /* RAS II part */
               occs2[cnt2++] = CalcInfo_->ras_orbs[1][i][m-k];
            }
         else if (Parameters_->opentype == PARM_OPENTYPE_SINGLET) {
            if (betsocc + CalcInfo_->socc[i] <= CalcInfo_->spab)
               betsocc += CalcInfo_->socc[i];
            else {
               j += CalcInfo_->socc[i];
               m += CalcInfo_->spab - betsocc;
               betsocc = CalcInfo_->spab;
               l = (j<k) ? j : k;
               for (; m<l; m++) /* RAS I part */
                  occs[cnt++] = CalcInfo_->ras_orbs[0][i][m];
               for (; m<j; m++) /* RAS II part */
                  occs2[cnt2++] = CalcInfo_->ras_orbs[1][i][m-k];
               }
            }
         } /* end socc[i] != 0 */
      }

   for (i=0; i<cnt2; i++) occs[cnt++] = occs2[i];
   CalcInfo_->ref_alp_rel = og_lex_addr(AlphaG_, occs, CalcInfo_->num_alp_expl,
      &(CalcInfo_->ref_alp_list));
   CalcInfo_->ref_alp = CalcInfo_->ref_alp_rel +
      AlphaG_->list_offset[CalcInfo_->ref_alp_list];
   alp_sym = CalcInfo_->ref_alp_list / AlphaG_->subgr_per_irrep;


   if (CalcInfo_->iopen) {
      betsocc = 0; cnt = 0; cnt2 = 0;
      zero_int_array(occs, CalcInfo_->num_alp_expl);
      for (i=0; i<CalcInfo_->nirreps; i++) {
         j = CalcInfo_->docc[i] - CalcInfo_->dropped_docc[i];
         k = CalcInfo_->ras_opi[0][i];
         l = (j<k) ? j : k;
         for (m=0; m<l; m++) /* RAS I part */
            occs[cnt++] = CalcInfo_->ras_orbs[0][i][m];
         for (; m<j; m++)    /* RAS II part */
            occs2[cnt2++] = CalcInfo_->ras_orbs[1][i][m-k];
         if (CalcInfo_->socc[i] != 0) {
            if (Parameters_->opentype == PARM_OPENTYPE_NONE)
               outfile->Printf("Warning: ignoring socc since opentype=none\n");
            else if (Parameters_->opentype == PARM_OPENTYPE_SINGLET &&
                     betsocc < CalcInfo_->spab) {
               if (betsocc + CalcInfo_->socc[i] <= CalcInfo_->spab) {
                  j += CalcInfo_->socc[i];
                  betsocc += CalcInfo_->socc[i];
                  l = (j<k) ? j : k;
                  for (; m<l; m++) /* RAS I part */
                     occs[cnt++] = CalcInfo_->ras_orbs[0][i][m];
                  for (; m<j; m++) /* RAS II part */
                     occs2[cnt2++] = CalcInfo_->ras_orbs[1][i][m-k];
                  }
               else {
                  j += CalcInfo_->spab - betsocc;
                  betsocc = CalcInfo_->spab;
                  l = (j<k) ? j : k;
                  for (; m<l; m++) /* RAS I part */
                     occs[cnt++] = CalcInfo_->ras_orbs[0][i][m];
                  for (; m<j; m++) /* RAS II part */
                     occs2[cnt2++] = CalcInfo_->ras_orbs[1][i][m-k];
                  }
               }
            }
         } /* end loop over irreps */

      for (i=0; i<cnt2; i++) occs[cnt++] = occs2[i];
      CalcInfo_->ref_bet_rel = og_lex_addr(BetaG_, occs, CalcInfo_->num_bet_expl,
         &(CalcInfo_->ref_bet_list));
      CalcInfo_->ref_bet = CalcInfo_->ref_bet_rel +
         BetaG_->list_offset[CalcInfo_->ref_bet_list];
      bet_sym = CalcInfo_->ref_bet_list / BetaG_->subgr_per_irrep;
      }

   else { /* closed-shell case */
      CalcInfo_->ref_bet = CalcInfo_->ref_alp;
      bet_sym = alp_sym;
      CalcInfo_->ref_bet_list = CalcInfo_->ref_alp_list;
      CalcInfo_->ref_bet_rel = CalcInfo_->ref_alp_rel;
      }

   if (Parameters_->ref_sym == -1)
      CalcInfo_->ref_sym = alp_sym ^ bet_sym;
   else
      CalcInfo_->ref_sym = Parameters_->ref_sym;


   /* form the new CIvect structure...watch for codex's with no strings...*/
   CIblks_->num_blocks = 0;
   CIblks_->num_alp_codes = AlphaG_->nirreps * AlphaG_->subgr_per_irrep;
   CIblks_->num_bet_codes = BetaG_->nirreps * BetaG_->subgr_per_irrep;
   CIblks_->decode = init_int_matrix(CIblks_->num_alp_codes, CIblks_->num_bet_codes);

   /* figure out the possible alpha/beta combinations */

   if (Parameters_->fci_strings) {
      for (irrep=0; irrep<AlphaG_->nirreps; irrep++) {
         betirrep = irrep ^ CalcInfo_->ref_sym;
         if (AlphaG_->sg[irrep][0].num_strings &&
             BetaG_->sg[betirrep][0].num_strings) {
            CIblks_->Ia_code[nblocks] = irrep;
            CIblks_->Ib_code[nblocks] = betirrep;
            CIblks_->Ia_size[nblocks] = AlphaG_->sg[irrep][0].num_strings;
            CIblks_->Ib_size[nblocks] = BetaG_->sg[betirrep][0].num_strings;
            nblocks++;
            }
         }
      }
   else {
      for (irrep=0; irrep<AlphaG_->nirreps; irrep++) {
         for (nalp1=AlphaG_->ras1_max; nalp1>=AlphaG_->ras1_min; nalp1--) {
            for (nalp3=0; nalp3<=AlphaG_->ras3_max; nalp3++) {
               for (nalp4=0; nalp4<=AlphaG_->ras4_max; nalp4++) {

                  alpcode=AlphaG_->decode[nalp1-AlphaG_->ras1_min][nalp3][nalp4];
                  if (alpcode == -1) continue;
                  nas = AlphaG_->sg[irrep][alpcode].num_strings;
                  if (!nas) continue;

                  for (nbet1=BetaG_->ras1_max; (nbet1>=Parameters_->ras1_min-nalp1
                     && nbet1>=BetaG_->ras1_min); nbet1--) {
                     for (nbet3=0; (nbet3<=Parameters_->ras3_max-nalp3 &&
                           nbet3<=BetaG_->ras3_max); nbet3++) {

                        if (!Parameters_->mixed && (nalp3 || nbet3) &&
                             (AlphaG_->ras1_max - nalp1 + BetaG_->ras1_max -
                              nbet1 > xlvl)) continue;

                        for (nbet4=0; (nbet4<=Parameters_->ras4_max-nalp4 &&
                           nbet4<=BetaG_->ras4_max); nbet4++) {

                           if (nalp3 + nalp4 + nbet3 + nbet4 >
                               Parameters_->ras34_max) continue;

                           if (!Parameters_->mixed4 && (nalp4 || nbet4) &&
                              (AlphaG_->ras1_max - nalp1 + BetaG_->ras1_max -
                               nbet1 > xlvl)) continue;

                           if (!Parameters_->cc_mixed &&
                               nalp3+nalp4+nbet3+nbet4 > xlvl &&
                               (nalp4>2 || nbet4>2 ||
                                nalp1-AlphaG_->ras1_min > 2 ||
                                nbet1-BetaG_->ras1_min > 2)) continue;

                           /* add special constraint if we want to kick
                              out any determinants which would not be
                              included in (spin-complete) DETCI translations
                              of Anna Krylov's SF-CI stuff
                              CDS 3/19/02
                            */
                           if (Parameters_->sf_restrict && (nalp4 || nbet4)
                               && (nalp1<AlphaG_->ras1_max ||
                                   nbet1<BetaG_->ras1_max) &&
                                  (nalp3+nbet3==0)) continue;

                           betcode =
                             BetaG_->decode[nbet1-BetaG_->ras1_min][nbet3][nbet4];
                           if (betcode == -1) continue;
                           betirrep = irrep ^ CalcInfo_->ref_sym;
                           nbs = BetaG_->sg[betirrep][betcode].num_strings;
                           if (!nbs) continue;

                           /* add nonstandard excitation types, such as
			      CID, CIST, CIDTQ, etc.
			      Parameters_->ex_allow[0] = Single excitations
			      Parameters_->ex_allow[1] = Double excitations
			      etc.
			      MLA 12/16/03
			   */
			   set = 0;
			   for (i=0; i<Parameters_->ex_lvl; i++) {
			     if (Parameters_->ex_allow[i] == 0) {
			       if (nalp3+nbet3 == i+1) set = 1;
			     }
			   }
                           if (set == 1) continue;

			   /* Remove all Singles: Testing
			   if (nalp3+nbet3==1)
			     continue;*/

			   /* Remove all Doubles: Testing
			   if (nalp3+nbet3==2)
			     continue;*/

                           CIblks_->Ia_code[nblocks] = AlphaG_->subgr_per_irrep *
                              irrep + alpcode;
                           CIblks_->Ia_size[nblocks] = nas;
                           CIblks_->Ib_code[nblocks] = BetaG_->subgr_per_irrep *
                              betirrep + betcode;
                           CIblks_->Ib_size[nblocks] = nbs;
                           nblocks++;
                           } /* end loop over nbet4 */
                        } /* end loop over nbet3 */
                     } /* end loop over nbet1 */
                  } /* end loop over nalp4 */
               } /* end loop over nalp3 */
            } /* end loop over nalp1 */
         } /* end loop over irrep */
      } /* end RAS case */

   /* get the first_iablk[] and last_iablk[] arrays */
   CIblks_->first_iablk = init_int_array(AlphaG_->nirreps);
   CIblks_->last_iablk = init_int_array(AlphaG_->nirreps);

   for (irrep=0; irrep < AlphaG_->nirreps; irrep++) {
      for (i=0; i<nblocks; i++) {
         j = CIblks_->Ia_code[i] / AlphaG_->subgr_per_irrep;
         if (j == irrep) break;
         }
      if (j == irrep) CIblks_->first_iablk[irrep] = i;
      else CIblks_->first_iablk[irrep] = -1;
      }
   for (irrep=0; irrep < AlphaG_->nirreps; irrep++) {
      maxblk = -2;
      for (i=0; i<nblocks; i++) {
         j = CIblks_->Ia_code[i] / AlphaG_->subgr_per_irrep;
         if (j == irrep && i > maxblk) maxblk = i;
         }
      CIblks_->last_iablk[irrep] = maxblk;
      }


   /* calculate the offsets */
   CIblks_->num_blocks = nblocks;

   if (nblocks > CI_BLK_MAX) {
      std::string str = "nblocks = ";
      str += std::to_string( nblocks) ;
      str += " > CI_BLK_MAX = ";
      str += std::to_string( CI_BLK_MAX) ;
      throw PsiException(str,__FILE__,__LINE__);
      }

   CIblks_->offset[0] = 0;
   for (i=1; i<nblocks; i++) {
      CIblks_->offset[i] = CIblks_->offset[i-1] +
         (BIGINT) CIblks_->Ia_size[i-1] *
         (BIGINT) CIblks_->Ib_size[i-1];
      }
   CIblks_->vectlen = CIblks_->offset[nblocks-1] +
                    (BIGINT) CIblks_->Ia_size[nblocks-1] *
                    (BIGINT) CIblks_->Ib_size[nblocks-1];

   if (print_) {
     outfile->Printf(
       "    The CI space requires %.0lf (%1.2E) determinants and %d blocks\n",
       (double) CIblks_->vectlen, (double) CIblks_->vectlen, nblocks);
     }

    // If we only have less than two dets we die
   if (CIblks_->vectlen < 2){
     throw PSIEXCEPTION("DETCI requires at least two determinants! Quitting...");
   }

   //if (print_)
   //   outfile->Printf( "\n   CI space contains %4d blocks\n", nblocks);

   /* set up the decode array */
   for (i=0; i<CIblks_->num_alp_codes; i++) {
     for (j=0; j<CIblks_->num_bet_codes; j++) {
        CIblks_->decode[i][j] = -1;
        for (k=0; k<nblocks; k++) {
           if (CIblks_->Ia_code[k] == i && CIblks_->Ib_code[k] == j)
              CIblks_->decode[i][j] = k;
           }
        }
     }
   free(occs);  free(occs2);

   if (Parameters_->print_ciblks){
     int blk;

     outfile->Printf( "\nCI Block Summary:\n");
     for (blk=0; blk<CIblks_->num_blocks; blk++) {
        outfile->Printf("Block %3d: Alp=%3d, Bet=%3d  Size = %4d x %4d = %ld\n",
                blk, CIblks_->Ia_code[blk], CIblks_->Ib_code[blk],
                CIblks_->Ia_size[blk], CIblks_->Ib_size[blk],
                (unsigned long) CIblks_->Ia_size[blk] *
                (unsigned long) CIblks_->Ib_size[blk]);

     };
   };

   // Copy a few things to make CIblks self contained
   CIblks_->subgr_per_irrep = AlphaG_->subgr_per_irrep;
   CIblks_->nirreps = CalcInfo_->nirreps;
   CIblks_->Ms0 = Parameters_->Ms0;
}

}} // namespace psi::detci
