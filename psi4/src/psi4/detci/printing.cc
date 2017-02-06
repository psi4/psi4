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
** PRINTING.C
**
** File contains routines associated with printing CI space, vectors, etc.
**
** C. David Sherrill
** Center for Computational Quantum Chemistry
** University of Georgia
*/

#include <cstdio>
#include <cmath>
#include <cstring>
#include <sstream>
#include <cctype> // for toupper()
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"

#include "psi4/detci/structs.h"
#include "psi4/detci/ciwave.h"

namespace psi { namespace detci {

#define FLAG_NONBLOCKS
#define MIN_COEFF 1.0E-13

std::string orb2lbl(int orbnum, struct calcinfo *Cinfo, int* orbs_per_irr);
extern int str_rel2abs(int relidx, int listnum, struct olsen_graph *Graph);


/*
** PRINT_VEC()
**
** Print the Most Important Determinants in the CI vector
** David Sherrill, February 1995
*/
void CIWavefunction::print_vec(unsigned int nprint, int *Ialist, int *Iblist,
      int *Iaidx, int *Ibidx, double *coeff)
{
   int Ia_abs, Ib_abs;

   /* print out the list of most important determinants */
   outfile->Printf("\n   The %d most important determinants:\n\n", nprint) ;
   for (unsigned int i=0; i<nprint; i++) {
      if (fabs(coeff[i]) < MIN_COEFF) continue;

      Ia_abs = str_rel2abs(Iaidx[i], Ialist[i], AlphaG_);
      Ib_abs = str_rel2abs(Ibidx[i], Iblist[i], BetaG_);

      #ifdef FLAG_NONBLOCKS
      int found_inblock=0;
      for (unsigned int j=0, found_inblock=0; j<H0block_->size; j++) {
         if (Iaidx[i] == H0block_->alpidx[j] &&
             Ibidx[i] == H0block_->betidx[j] &&
             Ialist[i] == H0block_->alplist[j] &&
             Iblist[i] == H0block_->betlist[j]) {
            found_inblock = 1;
            break;
            }
         }
      outfile->Printf("    %c", found_inblock ? ' ' : '*');
      #endif

      outfile->Printf("%4d  %10.6lf  (%5d,%5d)  ", i+1, coeff[i],
         Ia_abs, Ib_abs);

      std::string configstring(print_config(AlphaG_->num_orb, AlphaG_->num_el_expl, BetaG_->num_el_expl,
         alplist_[Ialist[i]] + Iaidx[i], betlist_[Iblist[i]] + Ibidx[i],
         AlphaG_->num_drc_orbs));

      outfile->Printf("%s\n", configstring.c_str());

      } /* end loop over important determinants */

   outfile->Printf("\n");

}



/*
** PRINT_CONFIG()
**
** Function prints a configuration, given a list of
** alpha and beta string occupancies.
**
** David Sherrill, February 1995
**
*/
std::string CIWavefunction::print_config(int nbf, int num_alp_el, int num_bet_el,
	     struct stringwr *stralp, struct stringwr *strbet, int num_drc_orbs)
{
   int j,k;
   int afound, bfound;

   std::ostringstream oss;

   /* loop over orbitals */
   for (j=0; j<nbf; j++) {

     std::string olabel(orb2lbl(j+num_drc_orbs, CalcInfo_, nmopi_)); /* get label for orbital j */

      for (k=0,afound=0; k<num_alp_el; k++) {
         if ((stralp->occs)[k] > j) break;
         else if ((stralp->occs)[k] == j) {
            afound = 1;
            break;
            }
         }
      for (k=0, bfound=0; k<num_bet_el; k++) {
         if ((strbet->occs)[k] > j) break;
         else if ((strbet->occs)[k] == j) {
            bfound = 1;
            break;
            }
         }
      if (afound || bfound) oss << olabel;

      if (afound && bfound) oss << "X ";
      else if (afound) oss << "A ";
      else if (bfound) oss << "B ";
      } /* end loop over orbitals */

   return oss.str();
}

/*
** orb2lbl(): Function converts an absolute orbital number into a
**    label such as 4A1, 2B2, etc.
**
** Parameters:
**    orbnum = orbital number in CI order (add frozen core!)
**    label  = place to put constructed label
**
** Needs Global (CalcInfo):
**    orbs_per_irrep = number of orbitals per irrep
**    order          = ordering array which maps a CI orbital to a
**                     Pitzer orbital (the opposite mapping from the
**                     "reorder" array)
**    irreps         = number of irreducible reps
**    nmo            = num of molecular orbitals
**    labels         = labels for all the irreps
**
** Notes:
**    If there are frozen core (FZC) orbitals, they are not included in the
**       CI numbering (unless they're "restricted" or COR orbitals).  This
**       is bothersome because some of the arrays constructed in the CI program
**       do start numbering from FZC orbitals.  Thus, pass orbnum as the CI
**       orbital PLUS any frozen core orbitals.
**
** Updated 8/16/95 by CDS
**    Allow it to handle more complex spaces...don't assume QT orbital order.
**    It was getting labels all mixed up for RAS's.
*/
std::string orb2lbl(int orbnum, struct calcinfo *Cinfo, int* orbs_per_irr)
{

   int ir, j, pitzer_orb, rel_orb;

   /* get Pitzer ordering */
   pitzer_orb = Cinfo->order[orbnum];

   if (pitzer_orb > Cinfo->nmo) {
      outfile->Printf( "(orb2lbl): pitzer_orb > nmo!\n");
      }

   for (ir=0,j=0; ir<Cinfo->nirreps; ir++) {
      if (orbs_per_irr[ir] == 0) continue;
      if (j + orbs_per_irr[ir] > pitzer_orb) break;
      else j += orbs_per_irr[ir];
      }
   rel_orb = pitzer_orb - j;

   if (rel_orb < 0) {
      outfile->Printf( "(orb2lbl): rel_orb < 0\n");
      }
   else if (rel_orb > orbs_per_irr[ir]) {
      outfile->Printf( "(orb2lbl): rel_orb > orbs_per_irrep[ir]\n");
      }

   std::ostringstream oss;
   oss << rel_orb+1 << Cinfo->labels[ir];
   return oss.str();
}


/*
** lbl2orb(): Function converts a label such as 4A1, 2B2, etc., to
**   an absolute orbital number.  The reverse of the above function
**   orb2lbl().
**
** Parameters:
**    orbnum = orbital number in CI order (add frozen core!)
**    label  = place to put constructed label
**
** Returns:
**    absolute orbital number for the correlated calc (less frozen)
**
*/
//int lbl2orb(char *orbstring)
//{
//
//   int ir, i, j, pitzer_orb, rel_orb, corr_orb;
//   char *s, *t;
//   char orblbl[10];
//
//   sscanf(orbstring, "%d%s", &rel_orb, orblbl);
//
//   /* get the irrep */
//   for (i=0,ir=-1; i<CalcInfo.nirreps; i++) {
//     s = orblbl;
//     t = CalcInfo.labels[i];
//     j = 0;
//     while ((toupper(*s) == toupper(*t)) && (j < strlen(orblbl))) {
//       s++;
//       t++;
//       j++;
//     }
//     if (j == strlen(orblbl)) {
//       ir = i;
//       break;
//     }
//   }
//
//   if (ir == -1) {
//     outfile->Printf( "lbl2orb: can't find label %s!\n", orblbl);
//     return(0);
//   }
//
//   /* get Pitzer ordering */
//   for (i=0,pitzer_orb=0; i<ir; i++) {
//     pitzer_orb += CalcInfo.orbs_per_irr[i];
//   }
//   pitzer_orb += rel_orb - 1; /* 1A1 is orbital 0 in A1 stack ... */
//
//   /* get correlated ordering */
//   corr_orb = CalcInfo.reorder[pitzer_orb];
//   corr_orb -= CalcInfo.num_drc_orbs;
//
//   if (corr_orb < 0 || corr_orb > CalcInfo.num_ci_orbs) {
//     outfile->Printf( "lbl2orb: error corr_orb out of bounds, %d\n",
//       corr_orb);
//     return(0);
//   }
//
//   return(corr_orb);
//
//}


//void eivout_t(double **a, double *b, int m, int n)
//   {
//      int ii,jj,kk,nn,ll;
//      int i,j,k;
//
//      ii=0;jj=0;
//L200:
//      ii++;
//      jj++;
//      kk=10*jj;
//      nn=n;
//      if (nn > kk) nn=kk;
//      ll = 2*(nn-ii+1)+1;
//      outfile->Printf("\n");
//      for (i=ii; i <= nn; i++) outfile->Printf("       %5d",i);
//      outfile->Printf("\n");
//      for (i=0; i < m; i++) {
//         outfile->Printf("\n%5d",i+1);
//         for (j=ii-1; j < nn; j++) {
//            outfile->Printf("%12.7f",a[j][i]);
//            }
//         }
//      outfile->Printf("\n");
//      outfile->Printf("\n     ");
//      for (j=ii-1; j < nn; j++) {
//         outfile->Printf("%12.7f",b[j]);
//         }
//      outfile->Printf("\n");
//      if (n <= kk) {
//         return;
//         }
//      ii=kk; goto L200;
//}


/*
** PRINT_CIBLK_SUMMARY()
**
** C. David Sherrill
** April 1996
**
*/
//void print_ciblk_summary(std::string out)
//{
//   int blk;
//
//   outfile->Printf( "\nCI Block Summary:\n");
//   for (blk=0; blk<CIblks.num_blocks; blk++) {
//      outfile->Printf("Block %3d: Alp=%3d, Bet=%3d  Size = %4d x %4d = %ld\n",
//              blk, CIblks.Ia_code[blk], CIblks.Ib_code[blk],
//              CIblks.Ia_size[blk], CIblks.Ib_size[blk],
//              (unsigned long) CIblks.Ia_size[blk] *
//              (unsigned long) CIblks.Ib_size[blk]);
//      }
//}

}} // namespace psi::detci
