/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace detci {

extern int og_lex_addr(struct olsen_graph *Graph, int *occs, int nel,
      int *listnum);
extern void print_ciblk_summary(FILE *outfile);


void set_ciblks(struct olsen_graph *AlphaG, struct olsen_graph *BetaG)
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

   CalcInfo.num_alp_str = AlphaG->num_str;
   CalcInfo.num_bet_str = BetaG->num_str;
   xlvl = Parameters.ex_lvl;

   if (Parameters.print_lvl) {
      fprintf(outfile, "\nThere are %d alpha strings\n", CalcInfo.num_alp_str);
      fprintf(outfile, "There are %d beta strings\n", CalcInfo.num_bet_str);
      }

   /* Get the occupations for the reference alpha and beta strings.
    * This used to be one line (occupy orbitals 0 to CalcInfo.num_alp_expl)
    * but now it's got to be more complicated due to orbital renumbering
    * (now by RAS then by symm then by energy) and due to open-shell
    * cases.  Assume that occupied orbs in reference belong only to RAS I 
    * or RAS II.  Should add a check for this sometime.  The following 
    * routine is by MLL and CDS 11/96.  Modified to put RAS II occs in
    * separate array because otherwise you get strings like |0 1 3 2>.
    */ 
   occs = init_int_array(CalcInfo.num_alp_expl);
   occs2 = init_int_array(CalcInfo.num_alp_expl);

   betsocc = 0; cnt = 0; cnt2 = 0;
   for (i=0; i<CalcInfo.nirreps; i++) {
      j = CalcInfo.docc[i];
      if (Parameters.fzc) j = j - CalcInfo.frozen_docc[i];
      k = CalcInfo.ras_opi[0][i];
      l = (j<k) ? j : k;
      for (m=0; m<l; m++) /* RAS I part */
         occs[cnt++] = CalcInfo.ras_orbs[0][i][m];
      for (; m<j; m++)    /* RAS II part */
         occs2[cnt2++] = CalcInfo.ras_orbs[1][i][m-k];
      if (CalcInfo.socc[i] != 0) {
         if (Parameters.opentype == PARM_OPENTYPE_NONE)
            fprintf(outfile,"Warning: ignoring socc since opentype=none\n");
         else if (Parameters.opentype == PARM_OPENTYPE_HIGHSPIN) {
            j += CalcInfo.socc[i];
            l = (j<k) ? j : k;
            for (; m<l; m++) /* RAS I part */
               occs[cnt++] = CalcInfo.ras_orbs[0][i][m];
            for (; m<j; m++) /* RAS II part */
               occs2[cnt2++] = CalcInfo.ras_orbs[1][i][m-k];
            }
         else if (Parameters.opentype == PARM_OPENTYPE_SINGLET) {
            if (betsocc + CalcInfo.socc[i] <= CalcInfo.spab)
               betsocc += CalcInfo.socc[i];
            else {
               j += CalcInfo.socc[i];
               m += CalcInfo.spab - betsocc;
               betsocc = CalcInfo.spab;
               l = (j<k) ? j : k;
               for (; m<l; m++) /* RAS I part */
                  occs[cnt++] = CalcInfo.ras_orbs[0][i][m];
               for (; m<j; m++) /* RAS II part */
                  occs2[cnt2++] = CalcInfo.ras_orbs[1][i][m-k];
               }                                
            }
         } /* end socc[i] != 0 */
      }

   for (i=0; i<cnt2; i++) occs[cnt++] = occs2[i];
   CalcInfo.ref_alp_rel = og_lex_addr(AlphaG, occs, CalcInfo.num_alp_expl,
      &(CalcInfo.ref_alp_list));
   CalcInfo.ref_alp = CalcInfo.ref_alp_rel + 
      AlphaG->list_offset[CalcInfo.ref_alp_list];
   alp_sym = CalcInfo.ref_alp_list / AlphaG->subgr_per_irrep;


   if (CalcInfo.iopen) {
      betsocc = 0; cnt = 0; cnt2 = 0;
      zero_int_array(occs, CalcInfo.num_alp_expl);
      for (i=0; i<CalcInfo.nirreps; i++) {
         j = CalcInfo.docc[i];
         if (Parameters.fzc) j = j - CalcInfo.frozen_docc[i];
         k = CalcInfo.ras_opi[0][i];
         l = (j<k) ? j : k;
         for (m=0; m<l; m++) /* RAS I part */
            occs[cnt++] = CalcInfo.ras_orbs[0][i][m];
         for (; m<j; m++)    /* RAS II part */
            occs2[cnt2++] = CalcInfo.ras_orbs[1][i][m-k];
         if (CalcInfo.socc[i] != 0) {
            if (Parameters.opentype == PARM_OPENTYPE_NONE)
               fprintf(outfile,"Warning: ignoring socc since opentype=none\n");
            else if (Parameters.opentype == PARM_OPENTYPE_SINGLET &&
                     betsocc < CalcInfo.spab) {
               if (betsocc + CalcInfo.socc[i] <= CalcInfo.spab) {
                  j += CalcInfo.socc[i];
                  betsocc += CalcInfo.socc[i];
                  l = (j<k) ? j : k;
                  for (; m<l; m++) /* RAS I part */
                     occs[cnt++] = CalcInfo.ras_orbs[0][i][m];
                  for (; m<j; m++) /* RAS II part */
                     occs2[cnt2++] = CalcInfo.ras_orbs[1][i][m-k];
                  }
               else {
                  j += CalcInfo.spab - betsocc;
                  betsocc = CalcInfo.spab;
                  l = (j<k) ? j : k;
                  for (; m<l; m++) /* RAS I part */
                     occs[cnt++] = CalcInfo.ras_orbs[0][i][m];
                  for (; m<j; m++) /* RAS II part */
                     occs2[cnt2++] = CalcInfo.ras_orbs[1][i][m-k];
                  }
               }
            }
         } /* end loop over irreps */

      for (i=0; i<cnt2; i++) occs[cnt++] = occs2[i];
      CalcInfo.ref_bet_rel = og_lex_addr(BetaG, occs, CalcInfo.num_bet_expl,
         &(CalcInfo.ref_bet_list));
      CalcInfo.ref_bet = CalcInfo.ref_bet_rel + 
         BetaG->list_offset[CalcInfo.ref_bet_list];
      bet_sym = CalcInfo.ref_bet_list / BetaG->subgr_per_irrep;
      }

   else { /* closed-shell case */
      CalcInfo.ref_bet = CalcInfo.ref_alp;
      bet_sym = alp_sym;
      CalcInfo.ref_bet_list = CalcInfo.ref_alp_list;
      CalcInfo.ref_bet_rel = CalcInfo.ref_alp_rel;
      }

   if (Parameters.ref_sym == -1)
      CalcInfo.ref_sym = alp_sym ^ bet_sym;
   else 
      CalcInfo.ref_sym = Parameters.ref_sym;


   /* form the new CIvect structure...watch for codex's with no strings...*/
   CIblks.num_blocks = 0;
   CIblks.num_alp_codes = AlphaG->nirreps * AlphaG->subgr_per_irrep;
   CIblks.num_bet_codes = BetaG->nirreps * BetaG->subgr_per_irrep;
   CIblks.decode = init_int_matrix(CIblks.num_alp_codes, CIblks.num_bet_codes);

   /* figure out the possible alpha/beta combinations */
 
   if (Parameters.fci_strings) {
      for (irrep=0; irrep<AlphaG->nirreps; irrep++) {
         betirrep = irrep ^ CalcInfo.ref_sym;
         if (AlphaG->sg[irrep][0].num_strings && 
             BetaG->sg[betirrep][0].num_strings) {
            CIblks.Ia_code[nblocks] = irrep;
            CIblks.Ib_code[nblocks] = betirrep;
            CIblks.Ia_size[nblocks] = AlphaG->sg[irrep][0].num_strings;
            CIblks.Ib_size[nblocks] = BetaG->sg[betirrep][0].num_strings;
            nblocks++;
            }
         }
      }
   else {
      for (irrep=0; irrep<AlphaG->nirreps; irrep++) {
         for (nalp1=AlphaG->ras1_max; nalp1>=AlphaG->ras1_min; nalp1--) {
            for (nalp3=0; nalp3<=AlphaG->ras3_max; nalp3++) {
               for (nalp4=0; nalp4<=AlphaG->ras4_max; nalp4++) {

                  alpcode=AlphaG->decode[nalp1-AlphaG->ras1_min][nalp3][nalp4];
                  if (alpcode == -1) continue;
                  nas = AlphaG->sg[irrep][alpcode].num_strings;
                  if (!nas) continue;

                  for (nbet1=BetaG->ras1_max; (nbet1>=Parameters.ras1_min-nalp1
                     && nbet1>=BetaG->ras1_min); nbet1--) {
                     for (nbet3=0; (nbet3<=Parameters.ras3_max-nalp3 &&
                           nbet3<=BetaG->ras3_max); nbet3++) {

                        if (!Parameters.mixed && (nalp3 || nbet3) &&
                             (AlphaG->ras1_max - nalp1 + BetaG->ras1_max -
                              nbet1 > xlvl)) continue;

                        for (nbet4=0; (nbet4<=Parameters.ras4_max-nalp4 &&
                           nbet4<=BetaG->ras4_max); nbet4++) {

                           if (nalp3 + nalp4 + nbet3 + nbet4 > 
                               Parameters.ras34_max) continue;

                           if (!Parameters.mixed4 && (nalp4 || nbet4) &&
                              (AlphaG->ras1_max - nalp1 + BetaG->ras1_max -
                               nbet1 > xlvl)) continue;

                           if (!Parameters.cc_mixed && 
                               nalp3+nalp4+nbet3+nbet4 > xlvl &&
                               (nalp4>2 || nbet4>2 || 
                                nalp1-AlphaG->ras1_min > 2 ||
                                nbet1-BetaG->ras1_min > 2)) continue;

                           /* add special constraint if we want to kick 
                              out any determinants which would not be
                              included in (spin-complete) DETCI translations
                              of Anna Krylov's SF-CI stuff
                              CDS 3/19/02
                            */
                           if (Parameters.sf_restrict && (nalp4 || nbet4)
                               && (nalp1<AlphaG->ras1_max || 
                                   nbet1<BetaG->ras1_max) && 
                                  (nalp3+nbet3==0)) continue;

                           betcode = 
                             BetaG->decode[nbet1-BetaG->ras1_min][nbet3][nbet4];
                           if (betcode == -1) continue;
                           betirrep = irrep ^ CalcInfo.ref_sym;
                           nbs = BetaG->sg[betirrep][betcode].num_strings;
                           if (!nbs) continue;

                           /* add nonstandard excitation types, such as
			      CID, CIST, CIDTQ, etc.
			      Parameters.ex_allow[0] = Single excitations
			      Parameters.ex_allow[1] = Double excitations
			      etc.
			      MLA 12/16/03
			   */
			   set = 0;
			   for (i=0; i<Parameters.ex_lvl; i++) {
			     if (Parameters.ex_allow[i] == 0) {
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
			   
                           CIblks.Ia_code[nblocks] = AlphaG->subgr_per_irrep * 
                              irrep + alpcode;
                           CIblks.Ia_size[nblocks] = nas;
                           CIblks.Ib_code[nblocks] = BetaG->subgr_per_irrep * 
                              betirrep + betcode;
                           CIblks.Ib_size[nblocks] = nbs;
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
   CIblks.first_iablk = init_int_array(AlphaG->nirreps);
   CIblks.last_iablk = init_int_array(AlphaG->nirreps);
 
   for (irrep=0; irrep < AlphaG->nirreps; irrep++) {
      for (i=0; i<nblocks; i++) {
         j = CIblks.Ia_code[i] / AlphaG->subgr_per_irrep;
         if (j == irrep) break;
         }
      if (j == irrep) CIblks.first_iablk[irrep] = i;
      else CIblks.first_iablk[irrep] = -1;
      } 
   for (irrep=0; irrep < AlphaG->nirreps; irrep++) {
      maxblk = -2;
      for (i=0; i<nblocks; i++) {
         j = CIblks.Ia_code[i] / AlphaG->subgr_per_irrep;
         if (j == irrep && i > maxblk) maxblk = i;
         }
      CIblks.last_iablk[irrep] = maxblk;
      }


   /* calculate the offsets */
   CIblks.num_blocks = nblocks;

   if (Parameters.print_lvl)
      fprintf(outfile, "CI space contains %4d blocks\n", nblocks);

   if (nblocks > CI_BLK_MAX) {
      std::string str = "nblocks = ";
      str += static_cast<std::ostringstream*>( &(std::ostringstream() << nblocks) )->str();
      str += " > CI_BLK_MAX = ";
      str += static_cast<std::ostringstream*>( &(std::ostringstream() << CI_BLK_MAX) )->str();
      throw PsiException(str,__FILE__,__LINE__);
      }

   CIblks.offset[0] = 0;
   for (i=1; i<nblocks; i++) {
      CIblks.offset[i] = CIblks.offset[i-1] +
         (BIGINT) CIblks.Ia_size[i-1] * 
         (BIGINT) CIblks.Ib_size[i-1];
      }
   CIblks.vectlen = CIblks.offset[nblocks-1] + 
                    (BIGINT) CIblks.Ia_size[nblocks-1] *
                    (BIGINT) CIblks.Ib_size[nblocks-1];

   if (Parameters.print_lvl) {
     fprintf(outfile,
       "\nCI space requires %.0lf determinants\n", (double) CIblks.vectlen);
     fflush(outfile);
     }

   /* set up the decode array */
   for (i=0; i<CIblks.num_alp_codes; i++) {
     for (j=0; j<CIblks.num_bet_codes; j++) {
        CIblks.decode[i][j] = -1;
        for (k=0; k<nblocks; k++) {
           if (CIblks.Ia_code[k] == i && CIblks.Ib_code[k] == j)
              CIblks.decode[i][j] = k;
           }
        } 
     }
   free(occs);  free(occs2);

   if (Parameters.print_ciblks) print_ciblk_summary(outfile);
}

}} // namespace psi::detci

