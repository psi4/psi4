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
#include <cctype> // for toupper()
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace detci {

#define CONFIG_STRING_MAX 200
#define FLAG_NONBLOCKS
#define MIN_COEFF 1.0E-13

void orb2lbl(int orbnum, char *label);
void print_config(int nbf, int num_alp_el, int num_bet_el, 
   struct stringwr *stralp, struct stringwr *strbet, 
   int num_fzc_orbs, char *outstring);
extern int str_rel2abs(int relidx, int listnum, struct olsen_graph *Graph);


/*
** PRINT_VEC()
**
** Print the Most Important Determinants in the CI vector
** David Sherrill, February 1995
*/
void print_vec(unsigned int nprint, int *Ialist, int *Iblist, 
      int *Iaidx, int *Ibidx, double *coeff,
      struct olsen_graph *AlphaG, struct olsen_graph *BetaG, 
      struct stringwr **alplist, struct stringwr **betlist,
      FILE *outfile)
{
   int i,j,k;
   unsigned long *index;
   double value, abs_value, minval=0.0;
   unsigned int alp_idx, bet_idx;
   char configstring[CONFIG_STRING_MAX];
   int Ia_abs, Ib_abs;

   #ifdef FLAG_NONBLOCKS 
   int found_inblock=0;
   #endif

   /* print out the list of most important determinants */
   fprintf(outfile, "\n\nThe %d most important determinants\n\n", nprint) ;
   for (i=0; i<nprint; i++) {

      if (fabs(coeff[i]) < MIN_COEFF) continue;

      Ia_abs = str_rel2abs(Iaidx[i], Ialist[i], AlphaG);
      Ib_abs = str_rel2abs(Ibidx[i], Iblist[i], BetaG);

      #ifdef FLAG_NONBLOCKS
      for (j=0, found_inblock=0; j<H0block.size; j++) {
         if (Iaidx[i] == H0block.alpidx[j] && 
             Ibidx[i] == H0block.betidx[j] &&
             Ialist[i] == H0block.alplist[j] && 
             Iblist[i] == H0block.betlist[j]) {
            found_inblock = 1;
            break; 
            }
         }
      fprintf(outfile, "%c", found_inblock ? ' ' : '*');
      #endif

      fprintf(outfile, "%4d  %10.6lf  (%5d,%5d)  ", i+1, coeff[i], 
         Ia_abs, Ib_abs);

      print_config(AlphaG->num_orb, AlphaG->num_el_expl, BetaG->num_el_expl,
         alplist[Ialist[i]] + Iaidx[i], betlist[Iblist[i]] + Ibidx[i],
         AlphaG->num_fzc_orbs, configstring);

      fprintf(outfile, "%s\n", configstring);

      } /* end loop over important determinants */

   fprintf(outfile, "\n\n");
 
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
void print_config(int nbf, int num_alp_el, int num_bet_el, 
   struct stringwr *stralp, struct stringwr *strbet, int num_fzc_orbs,
   char *outstring)
{
   int j,k;
   int afound, bfound;
   char olabel[10];

   sprintf(outstring, "");

   /* loop over orbitals */
   for (j=0; j<nbf; j++) {

      orb2lbl(j+num_fzc_orbs, olabel); /* get label for orbital j */

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
      if (afound || bfound) strcat(outstring, olabel);

      if (afound && bfound) strcat(outstring, "X  ");
      else if (afound) strcat(outstring, "A  ");
      else if (bfound) strcat(outstring, "B  ");
      } /* end loop over orbitals */

} 


/*
** PRINT_CI_SPACE()
** 
** This function is for debugging purposes.  It prints the 
** CI space and the associated single-replacement lists.
**
** Arguments:
**    strlist     = list of alpha/beta strings
**    num_strings = number of strings in list
**    nirreps     = number of irreducible representations in molecular pt grp
**    strtypes    = number of possible string types (nirreps * ncodes)
**    nel         = number of electrons explicitly included
**    outfile     = file to print to
*/
void print_ci_space(struct stringwr *strlist, int num_strings, 
      int nirreps, int strtypes, int nel, FILE *outfile) 
{
   int i, j, strsym, cnt=0 ;

   while (cnt != num_strings) {
      fprintf(outfile, "\nString %4d (", cnt++);
      for (i=0; i<nel; i++)
         fprintf(outfile, "%2d ", (int) (strlist->occs)[i]) ;
      fprintf(outfile, ")\n");

      if (!Parameters.repl_otf) {
         fprintf(outfile, "   Links:\n") ;
         for (strsym=0; strsym < strtypes; strsym++) {
            for (j=0; j<strlist->cnt[strsym]; j++) {
               fprintf(outfile, "   %3d [%3d] %c (%2d %3d)   %d\n",
                  strlist->ij[strsym][j], 
                  strlist->oij[strsym][j],
                  (strlist->sgn[strsym][j] == 1) ? '+' : '-', 
                  strsym, strlist->ridx[strsym][j], 
                  (int) strlist->sgn[strsym][j]);
               }
            } /* end loop over strsym */
         }
      strlist++;
      }
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
void orb2lbl(int orbnum, char *label)
{

   int ir, i, j, pitzer_orb, rel_orb;

   /* get Pitzer ordering */
   pitzer_orb = CalcInfo.order[orbnum];
   
   if (pitzer_orb > CalcInfo.nmo) {
      fprintf(outfile, "(orb2lbl): pitzer_orb > nmo!\n");
      }

   for (ir=0,j=0; ir<CalcInfo.nirreps; ir++) {
      if (CalcInfo.orbs_per_irr[ir] == 0) continue;
      if (j + CalcInfo.orbs_per_irr[ir] > pitzer_orb) break;
      else j += CalcInfo.orbs_per_irr[ir];
      }
   rel_orb = pitzer_orb - j;

   if (rel_orb < 0) {
      fprintf(outfile, "(orb2lbl): rel_orb < 0\n");
      }
   else if (rel_orb > CalcInfo.orbs_per_irr[ir]) {
      fprintf(outfile, "(orb2lbl): rel_orb > orbs_per_irrep[ir]\n");
      }
 
   sprintf(label, "%d%s", rel_orb+1, CalcInfo.labels[ir]);

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
int lbl2orb(char *orbstring)
{

   int ir, i, j, pitzer_orb, rel_orb, corr_orb;
   char *s, *t;
   char orblbl[10];

   sscanf(orbstring, "%d%s", &rel_orb, orblbl);

   /* get the irrep */
   for (i=0,ir=-1; i<CalcInfo.nirreps; i++) {
     s = orblbl;
     t = CalcInfo.labels[i];
     j = 0;
     while ((toupper(*s) == toupper(*t)) && (j < strlen(orblbl))) {
       s++; 
       t++;
       j++;
     }
     if (j == strlen(orblbl)) {
       ir = i;
       break;
     }
   }

   if (ir == -1) {
     fprintf(outfile, "lbl2orb: can't find label %s!\n", orblbl);
     return(0);
   }

   /* get Pitzer ordering */
   for (i=0,pitzer_orb=0; i<ir; i++) {
     pitzer_orb += CalcInfo.orbs_per_irr[i];
   }
   pitzer_orb += rel_orb - 1; /* 1A1 is orbital 0 in A1 stack ... */

   /* get correlated ordering */
   corr_orb = CalcInfo.reorder[pitzer_orb];

   /* probably need to subtract frozen here */
   corr_orb -= CalcInfo.num_fzc_orbs;

   if (corr_orb < 0 || corr_orb > CalcInfo.num_ci_orbs) {
     fprintf(outfile, "lbl2orb: error corr_orb out of bounds, %d\n", 
       corr_orb);
     return(0);
   }

   return(corr_orb);

}


void eivout_t(double **a, double *b, int m, int n, FILE *out)
   {
      int ii,jj,kk,nn,ll;
      int i,j,k;

      ii=0;jj=0;
L200:
      ii++;
      jj++;
      kk=10*jj;
      nn=n;
      if (nn > kk) nn=kk;
      ll = 2*(nn-ii+1)+1;
      fprintf (out,"\n");
      for (i=ii; i <= nn; i++) fprintf(out,"       %5d",i);
      fprintf (out,"\n");
      for (i=0; i < m; i++) {
         fprintf (out,"\n%5d",i+1);
         for (j=ii-1; j < nn; j++) {
            fprintf (out,"%12.7f",a[j][i]);
            }
         }
      fprintf (out,"\n");
      fprintf (out,"\n     ");
      for (j=ii-1; j < nn; j++) {
         fprintf(out,"%12.7f",b[j]);
         }
      fprintf (out,"\n");
      if (n <= kk) {
         fflush(out);
         return;
         }
      ii=kk; goto L200;
}


/*
** PRINT_CIBLK_SUMMARY()
**
** C. David Sherrill
** April 1996
**
*/
void print_ciblk_summary(FILE *outfile)
{
   int blk;

   fprintf(outfile, "\nCI Block Summary:\n");
   for (blk=0; blk<CIblks.num_blocks; blk++) {
      fprintf(outfile,"Block %3d: Alp=%3d, Bet=%3d  Size = %4d x %4d = %ld\n", 
              blk, CIblks.Ia_code[blk], CIblks.Ib_code[blk], 
              CIblks.Ia_size[blk], CIblks.Ib_size[blk],
              (unsigned long) CIblks.Ia_size[blk] * 
              (unsigned long) CIblks.Ib_size[blk]);
      }
}

/*
** WRITE_ENERGY
**
** This routine writes out the energies to an ASCII file
*/
void write_energy(int nroots, double *evals, double offset)
{
  FILE *efile;
  int i;

  ffile(&efile,"detci_energies.dat",1);
  for (i=0; i<nroots; i++) { 
    fprintf(efile, "%8.6lf ", evals[i]+offset);
  }
  fprintf(efile, "\n");
  fclose(efile);
}

}} // namespace psi::detci

