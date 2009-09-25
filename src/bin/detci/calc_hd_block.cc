/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here 
*/

#define EXTERN
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#include "globals.h"

namespace psi { namespace detci {

extern int calc_orb_diff(int cnt, unsigned char *I, unsigned char *J,
   int *I_alpha_diff, int *J_alpha_diff, int *sign, int *same, int extended);
extern int *ioff;

/* C "GLOBAL" VARIABLES FOR THIS MODULE */
extern struct stringwr **alplist;
extern struct stringwr **betlist;

 
#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))


/*
** calc_hd_block(): Function calculates a block of H0, the diagonal elements of
**    the Hamiltonian matrix.
**
** Parameters:
**    alplist_local = list of alpha strings with replacements (used to get occs)
**    betlist_local = list of beta strings with replacements
**    nas     = number of alpha strings in list
**    nbs     = number of beta strings in list
**    H0      = matrix to hold results (stored as H0[alpidx][betidx])
**    oei     = one-electron integrals
**    tei     = two-electron integrals
**    na      = number of explicit alpha electrons 
**    nb      = number of explicit beta electrons
**    nbf     = number of orbitals in CI
**    efzc    = frozen core energy
**
*/ 
void calc_hd_block(struct stringwr *alplist_local, struct stringwr *betlist_local,
      double **H0, double *oei, double *tei, double efzc,
      int nas, int nbs, int na, int nb, int nbf)
{
   int acnt, bcnt;
   int a1, a2, b1, b2;
   int i,j, ii, iii, jj, ij, iijj, ijij;
   double value;
   struct stringwr *betlist0;

   betlist0 = betlist_local;

   for (acnt=0; acnt<nas; acnt++) {
      
      for (bcnt=0, betlist_local=betlist0; bcnt<nbs; bcnt++) {

         /* add frozen core energy first */
/************************************************/
         value = efzc; 

         /* loop over alpha occs */
         for (a1=0; a1<na; a1++) {
            i = (int) alplist_local->occs[a1];
            ii = ioff[i] + i;
            value += oei[ii];  
            /* fprintf(outfile,"oei[%d] = %lf\n",ii,oei[ii]); */ 
            iii = ioff[ii];

            for (a2=0; a2<a1; a2++) {
               j = (int) alplist_local->occs[a2];
               jj = ioff[j] + j;
               iijj = iii + jj;
               ij = ioff[i] + j;
               ijij = ioff[ij] + ij;
               value += tei[iijj] - tei[ijij]; 
               }

            for (b1=0; b1<nb; b1++) {
               j = (int) betlist_local->occs[b1];
               jj = ioff[j] + j;
               iijj = ioff[MAX0(ii,jj)] + MIN0(ii,jj);
               value += tei[iijj]; 
               }
           }

         for (b1=0; b1<nb; b1++) {
            i = (int) betlist_local->occs[b1];
            ii = ioff[i] + i;
            value += oei[ii]; 
            iii = ioff[ii];

            for (b2=0; b2<b1; b2++) {
               j = (int) betlist_local->occs[b2];
               jj = ioff[j] + j;
               iijj = iii + jj;
               ij = ioff[i] + j;
               ijij = ioff[ij] + ij;
               value += tei[iijj] - tei[ijij]; 
               }
            } 

         H0[acnt][bcnt] = value;
      /*   
         fprintf(outfile,"H0[%d][%d] = %lf\n",acnt,bcnt,value); 
      */ 
         betlist_local++;
         } /* end loop over bcnt */

      alplist_local++;
      }

}
/*
** calc_hd_block_ave(): Function calculates a block of H0 and the diagonal elements
** of the Hamiltonian matrix averaged over spin-coupling sets to correct any
** spin contamination of the c and sigma vectors. 
**
** Parameters:
**    alplist_local = list of alpha strings with replacements (used to get occs)
**    betlist_local = list of beta strings with replacements
**    nas     = number of alpha strings in list
**    nbs     = number of beta strings in list
**    H0      = matrix to hold results (stored as H0[alpidx][betidx])
**    oei     = one-electron integrals
**    tei     = two-electron integrals
**    na      = number of explicit alpha electrons 
**    nb      = number of explicit beta electrons
**    nbf     = number of orbitals in CI
**    efzc    = frozen core energy
**
*/ 
void calc_hd_block_ave(struct stringwr *alplist_local, struct stringwr *betlist_local,
      double **H0, double *tf_oei, double *tei, double efzc,
      int nas, int nbs, int na, int nb, int nbf)
{
   int acnt, bcnt;
   int a1, a2, a3, b1, b2, b3;
   int i,j, ii, iii, jj, ij, iijj, ijij;
   double value, tval, tval2, Kave;
   struct stringwr *betlist0;
   double k_total; /* total number of K ints in energy expression */
   int k_combo; /* total combination of unique K ints over spin-coupling set */
   int *unique_occs; /* the uniquely occupied orbitals for a given determinant */
   int num_el;  /* total number of electrons explicitly treated */
   int num_unique; /* number of unique orbitals */
   betlist0 = betlist_local;
   

   k_total = combinations(na,2) + combinations(nb,2); 

   num_el = na + nb;
   unique_occs = init_int_array(num_el);

   for (acnt=0; acnt<nas; acnt++) {
      
      for (bcnt=0, betlist_local=betlist0; bcnt<nbs; bcnt++) {

         /* add frozen core energy first */
         value = efzc; 

         /* loop over alpha occs */
         for (a1=0; a1<na; a1++) {
            i = (int) alplist_local->occs[a1];
            ii = ioff[i] + i;
            /* h_ii bar alpha alpha */
            value += tf_oei[ii];
            /* fprintf(outfile,"tf_oei[%d] = %lf\n",ii,tf_oei[ii]); */
            iii = ioff[ii];

            /* loop over alpha occs */
            for (a2=0; a2<a1; a2++) {
               j = (int) alplist_local->occs[a2];
               jj = ioff[j] + j;
               iijj = iii + jj;
               /* J alpha alpha */
               value += tei[iijj]; 
               }

            /* loop over beta occs */
            for (b1=0; b1<nb; b1++) {
               j = (int) betlist_local->occs[b1];
               jj = ioff[j] + j;
               iijj = ioff[MAX0(ii,jj)] + MIN0(ii,jj);
               value += tei[iijj]; 
               }
           }

         /* loop over beta occs */
         for (b1=0; b1<nb; b1++) {
            i = (int) betlist_local->occs[b1];
            ii = ioff[i] + i;
            value += tf_oei[ii];
            /* fprintf(outfile,"tf_oei[%d] = %lf\n",ii,tf_oei[ii]); */
            iii = ioff[ii];

            /* loop over beta occs */
            for (b2=0; b2<b1; b2++) {
               j = (int) betlist_local->occs[b2];
               jj = ioff[j] + j;
               iijj = iii + jj;
               ij = ioff[i] + j; 
               ijij = ioff[ij] + ij;
               value += tei[iijj];  
               }
            }

         /* determine average K over spin-coupling set */
         num_unique = 0;
         for (a1=0; a1<na; a1++) unique_occs[num_unique++] = (int) alplist_local->occs[a1];
         /* for (j=0; j<num_unique; j++) 
            fprintf(outfile,"unique_occs[%d] = %d\n",j,unique_occs[j]); */
            for (b1=0; b1<nb; b1++) {
               j = (int) betlist_local->occs[b1];
               for (a1=0; a1<na; a1++) {
                  if (j==unique_occs[a1]) break;
                  if (a1==(na-1)) unique_occs[num_unique++] = j;
                  }
               }
         /* fprintf(outfile,"num_unique = %d\n",num_unique);
         fprintf(outfile,"num_el = %d\n",num_el);
         */
         if (num_unique>num_el) fprintf(outfile,"WARNING: The number of explicit electrons" \
                             "!= num_el\n");
               
       /*   
         for (j=0; j<na; j++) 
            fprintf(outfile,"alp_occs[%d] = %d\n",j,(int)alplist_local->occs[j]);
         for (j=0; j<nb; j++) 
            fprintf(outfile,"bet_occs[%d] = %d\n",j,(int)betlist_local->occs[j]);
         for (j=0; j<num_unique; j++) 
            fprintf(outfile,"unique_occs[%d] = %d\n",j,unique_occs[j]);
       */
 
         Kave = 0.0;
         for (a1=0; a1<num_unique; a1++) {
            i = unique_occs[a1];
            for (b1=0; b1<a1; b1++) {
               j = unique_occs[b1];
               ij = ioff[MAX0(i,j)] + MIN0(i,j);
               ijij = ioff[ij] + ij;
               Kave += tei[ijij];
               /* fprintf(outfile,"tei[%d] = %lf\n",ijij,tei[ijij]); */ 
               }
            }
         
         /* fprintf(outfile,"num_unique = %d\n",num_unique);
         fprintf(outfile,"ioff[num_unique-1] = %d\n",ioff[num_unique]);
         fprintf(outfile,"k_total = %d\n",k_total);
         */

         if (num_unique > 1) Kave /= ioff[num_unique-1];
         value -= 0.5 * Kave * k_total; 
         /* fprintf(outfile,"Kave = %lf\n",Kave); */

         if (Parameters.print_lvl > 5) {
           fprintf(outfile,"acnt = %d\t bcnt = %d\n",acnt,bcnt); 
           fprintf(outfile,"tval = %lf\n",tval);
           for(a1=0; a1<na; a1++)
             fprintf(outfile," %d",alplist_local->occs[a1]);
           fprintf(outfile," \n");
           for(b1=0; b1<nb; b1++)
             fprintf(outfile," %d",betlist_local->occs[b1]);
           fprintf(outfile," \n");
           } 

         H0[acnt][bcnt] = value;
         /* fprintf(outfile,"H0[%d][%d] = %lf\n",acnt,bcnt,value); */
         betlist_local++;
         } /* end loop over bcnt */

      alplist_local++;
      }

}

/*
** calc_hd_block_orbenergy(): Function calculates a block of H0 and the diagonal elements
** of the Hamiltonian matrix as the sum of orbital energies. 
**
** Parameters:
**    alplist_local = list of alpha strings with replacements (used to get occs)
**    betlist_local = list of beta strings with replacements
**    nas     = number of alpha strings in list
**    nbs     = number of beta strings in list
**    H0      = matrix to hold results (stored as H0[alpidx][betidx])
**    oei     = one-electron integrals
**    tei     = two-electron integrals
**    na      = number of explicit alpha electrons 
**    nb      = number of explicit beta electrons
**    nbf     = number of orbitals in CI
**    efzc    = frozen core energy
**
*/ 
void calc_hd_block_orbenergy(struct stringwr *alplist_local, 
      struct stringwr *betlist_local, double **H0, double *oei, 
      double *tei, double efzc, int nas, int nbs, int na, int nb, int nbf)
{
   int acnt, bcnt;
   int a1, b1, i,j; 
   double value, tval;
   struct stringwr *betlist0, *alplist0;
   double *orb_e_diff_alp, *orb_e_diff_bet;
   double sum_orb_energies = 0.0;

   betlist0 = betlist_local;
   alplist0 = alplist_local; 

   orb_e_diff_alp = init_array(nas);
   orb_e_diff_bet = init_array(nbs);
  /* if (Parameters.Ms0) orb_e_diff_bet = &orb_e_diff_alp;
   else orb_e_diff_bet = init_array(CalcInfo.num_bet_str);
  */

   for (acnt=0; acnt<nas; acnt++) {
      orb_e_diff_alp[acnt] = 0.0;
      for (a1=0; a1<na; a1++) {
         i = (int) alplist_local->occs[a1];
         i += CalcInfo.num_fzc_orbs;
         if(Parameters.zaptn) 
           orb_e_diff_alp[acnt] += CalcInfo.scfeigvala[i];
         else
           orb_e_diff_alp[acnt] += CalcInfo.scfeigval[i];
         }
      alplist_local++;
      }

   for (bcnt=0; bcnt<nbs; bcnt++) {
      orb_e_diff_bet[bcnt] = 0.0;
      for (b1=0; b1<nb; b1++) {
         j = (int) betlist_local->occs[b1];
         j += CalcInfo.num_fzc_orbs;
         if(Parameters.zaptn) 
           orb_e_diff_bet[bcnt] += CalcInfo.scfeigvalb[j];
         else
           orb_e_diff_bet[bcnt] += CalcInfo.scfeigval[j];
         }
      betlist_local++;
      }

   alplist_local = alplist0;
   betlist_local = betlist0;

   for (acnt=0; acnt<nas; acnt++) {
         tval = efzc + orb_e_diff_alp[acnt]; 
      for (bcnt=0; bcnt<nbs; bcnt++) {
         value = orb_e_diff_bet[bcnt] + tval; 
         H0[acnt][bcnt] = value;
        /* 
         fprintf(outfile,"H0[%d][%d] = %lf\n",acnt,bcnt,value); 
        */ 
         betlist_local++;
         } /* end loop over bcnt */
      alplist_local++;
      }

/* Free up memory */
free(orb_e_diff_alp);
free(orb_e_diff_bet);


}

/*
** calc_hd_block_evangelisti(): Function calculates a block of H0 and the diagonal elements
** of the Hamiltonian matrix averaged over spin-coupling sets to correct any
** spin contamination of the c and sigma vectors. 
**
** Parameters:
**    alplist_local = list of alpha strings with replacements (used to get occs)
**    betlist_local = list of beta strings with replacements
**    nas     = number of alpha strings in list
**    nbs     = number of beta strings in list
**    H0      = matrix to hold results (stored as H0[alpidx][betidx])
**    oei     = one-electron integrals
**    tei     = two-electron integrals
**    na      = number of explicit alpha electrons 
**    nb      = number of explicit beta electrons
**    nbf     = number of orbitals in CI
**    efzc    = frozen core energy
**
*/ 
void calc_hd_block_evangelisti(struct stringwr *alplist_local, struct stringwr *betlist_local,
      double **H0, double *tf_oei, double *tei, double efzc,
      int nas, int nbs, int na, int nb, int nbf)
{
   int acnt, bcnt;
   int a1, b1, i,j; 
   double value, tval;
   struct stringwr *betlist0, *alplist0;
   double *orb_e_diff_alp, *orb_e_diff_bet;
   int num_alp_diff, num_bet_diff;
   int **orb_diff, *jnk;
   int sign;

   betlist0 = betlist_local;
   alplist0 = alplist_local; 

   orb_diff = init_int_matrix(2,na);
   jnk = init_int_array(na);
   orb_e_diff_alp = init_array(nas);
   orb_e_diff_bet = init_array(nbs);

   for (acnt=0; acnt<nas; acnt++) {
      orb_e_diff_alp[acnt] = 0.0;
      num_alp_diff = calc_orb_diff(na,
                     alplist[CalcInfo.ref_alp_list][CalcInfo.ref_alp_rel].occs, 
                     alplist_local->occs, orb_diff[0], orb_diff[1], &sign,
                     jnk, 1);
      for (a1=0; a1<num_alp_diff; a1++) {
         i = orb_diff[0][a1]; 
         j = orb_diff[1][a1]; 
         i += CalcInfo.num_fzc_orbs;
         j += CalcInfo.num_fzc_orbs;
         orb_e_diff_alp[acnt] += CalcInfo.scfeigval[j] 
                                 - CalcInfo.scfeigval[i]; 
         }
      alplist_local++;
      }

   for (bcnt=0; bcnt<nbs; bcnt++) {
      orb_e_diff_bet[bcnt] = 0.0;
      num_bet_diff = calc_orb_diff(nb, 
                     betlist[CalcInfo.ref_bet_list][CalcInfo.ref_bet_rel].occs,
                     betlist_local->occs, orb_diff[0], orb_diff[1], &sign, 
                     jnk, 1);
      for (b1=0; b1<num_bet_diff; b1++) {
         i = orb_diff[0][b1];
         j = orb_diff[1][b1];  
         i += CalcInfo.num_fzc_orbs;
         j += CalcInfo.num_fzc_orbs;
         orb_e_diff_bet[bcnt] += CalcInfo.scfeigval[j]
                                 - CalcInfo.scfeigval[i];
         }
      betlist_local++;
      } 

   alplist_local = alplist0;
   betlist_local = betlist0;

   for (acnt=0; acnt<nas; acnt++) {
         /* add frozen core energy first */
         tval = CalcInfo.escf - CalcInfo.enuc; 
         tval += orb_e_diff_alp[acnt]; 
      for (bcnt=0; bcnt<nbs; bcnt++) {
         value = 0.0;
         value = orb_e_diff_bet[bcnt] + tval; 
         H0[acnt][bcnt] = value;
         /* fprintf(outfile,"H0[%d][%d] = %lf\n",acnt,bcnt,value); */ 
         betlist_local++;
         } /* end loop over bcnt */
      alplist_local++;
      }

/* Free memory */
/*
free(jnk);
free(orb_e_diff_alp);
free(orb_e_diff_bet);
free(orb_diff);
*/


}


/*
** calc_hd_block_mll(): Function calculates a block of H0 and the diagonal elements
** of the Hamiltonian matrix as the sum of orbital energies. 
**
** Parameters:
**    alplist_local = list of alpha strings with replacements (used to get occs)
**    betlist_local = list of beta strings with replacements
**    nas     = number of alpha strings in list
**    nbs     = number of beta strings in list
**    H0      = matrix to hold results (stored as H0[alpidx][betidx])
**    oei     = one-electron integrals
**    tei     = two-electron integrals
**    na      = number of explicit alpha electrons 
**    nb      = number of explicit beta electrons
**    nbf     = number of orbitals in CI
**    efzc    = frozen core energy
**
*/ 
void calc_hd_block_mll(struct stringwr *alplist_local, 
      struct stringwr *betlist_local, double **H0, double *oei, 
      double *tei, double efzc, int nas, int nbs, int na, int nb, int nbf)
{
   int acnt, bcnt;
   int a1, b1, i,j, i_offset, j_offset, ii, jj; 
   double value, tval;
   struct stringwr *betlist0, *alplist0;
   double *orb_e_diff_alp, *orb_e_diff_bet;
   double *oei_alp, *oei_bet, *eigval;

   betlist0 = betlist_local;
   alplist0 = alplist_local; 

   oei_alp = init_array(nas);
   oei_bet = init_array(nbs);
   orb_e_diff_alp = init_array(nas);
   orb_e_diff_bet = init_array(nbs);
  /* if (Parameters.Ms0) orb_e_diff_bet = &orb_e_diff_alp;
   else orb_e_diff_bet = init_array(nbs);
  */

   for (acnt=0; acnt<nas; acnt++) {
      orb_e_diff_alp[acnt] = oei_alp[acnt] = 0.0;
      for (a1=0; a1<na; a1++) {
         i = (int) alplist_local->occs[a1];
         ii = ioff[i] + i;
         i_offset = i + CalcInfo.num_fzc_orbs;
         oei_alp[acnt] += oei[ii]; 
         orb_e_diff_alp[acnt] += CalcInfo.scfeigval[i_offset] - oei[ii];
         }
      alplist_local++;
      }

   for (bcnt=0; bcnt<nbs; bcnt++) {
      orb_e_diff_bet[bcnt] = oei_bet[bcnt] = 0.0;
      for (b1=0; b1<nb; b1++) {
         j = (int) betlist_local->occs[b1];
         jj = ioff[j] + j;
         j_offset = j + CalcInfo.num_fzc_orbs;
         oei_bet[bcnt] += oei[jj];
         orb_e_diff_bet[bcnt] += CalcInfo.scfeigval[j_offset] - oei[jj];
         }
      betlist_local++;
      }

   alplist_local = alplist0;
   betlist_local = betlist0;

   for (acnt=0; acnt<nas; acnt++) {
         tval = efzc + 0.5 * orb_e_diff_alp[acnt] + oei_alp[acnt]; 
      for (bcnt=0; bcnt<nbs; bcnt++) {
         value = 0.5 * orb_e_diff_bet[bcnt] + oei_bet[bcnt] + tval; 
         H0[acnt][bcnt] = value;
         betlist_local++;
         } /* end loop over bcnt */
      alplist_local++;
      }

 free(oei_alp);
 free(oei_bet);
 free(orb_e_diff_alp);
 free(orb_e_diff_bet);
}
/*
** calc_hd_block_z_ave(): Function calculates a block of H0 and the diagonal elements
** of the Hamiltonian matrix averaged over spin-coupling sets to correct any
** spin contamination of the c and sigma vectors.
**
** Parameters:
**    alplist_local = list of alpha strings with replacements (used to get occs)
**    betlist_local = list of beta strings with replacements
**    nas     = number of alpha strings in list
**    nbs     = number of beta strings in list
**    H0      = matrix to hold results (stored as H0[alpidx][betidx])
**    oei     = one-electron integrals
**    tei     = two-electron integrals
**    na      = number of explicit alpha electrons
**    nb      = number of explicit beta electrons
**    nbf     = number of orbitals in CI
**    efzc    = frozen core energy
**
*/
void calc_hd_block_z_ave(struct stringwr *alplist_local, struct stringwr *betlist_local, 
double **H0, double pert_param, double *tei, double efzc, int nas, int nbs, int na, 
int nb, int nbf)
{
   int acnt, bcnt;
   int a1, a2, a3, b1, b2, b3;
   int i,j, ii, iii, jj, ij, iijj, ijij;
   double value, tval, tval2, Kave;
   struct stringwr *betlist0;
   double k_total; /* total number of K ints in energy expression */
   int k_combo; /* total combination of unique K ints over spin-coupling set */
   int *unique_occs; /* the uniquely occupied orbitals for a given determinant */
   int num_el;  /* total number of electrons explicitly treated */
   int num_unique; /* number of unique orbitals */
   betlist0 = betlist_local;


   k_total = combinations(na,2) + combinations(nb,2);
   num_el = na + nb;
   unique_occs = init_int_array(num_el);

   for (acnt=0; acnt<nas; acnt++) {

      for (bcnt=0, betlist_local=betlist0; bcnt<nbs; bcnt++) {

         /* add frozen core energy first */
         value = efzc;

         /* loop over alpha occs */
         for (a1=0; a1<na; a1++) {
            i = (int) alplist_local->occs[a1];
            value += CalcInfo.scfeigval[i+CalcInfo.num_fzc_orbs];
            ii = ioff[i] + i;
            /* h_ii bar alpha alpha */
            iii = ioff[ii];

            /* loop over alpha occs */
            for (a2=0; a2<a1; a2++) {
               j = (int) alplist_local->occs[a2];
               jj = ioff[j] + j;
               iijj = iii + jj;
               /* J alpha alpha */
               value -= pert_param * tei[iijj];
               }

            /* loop over beta occs */
            for (b1=0; b1<nb; b1++) {
               j = (int) betlist_local->occs[b1];
               jj = ioff[j] + j;
               iijj = ioff[MAX0(ii,jj)] + MIN0(ii,jj);
               value -= pert_param * tei[iijj];
               }
           }

         /* loop over beta occs */
         for (b1=0; b1<nb; b1++) {
            i = (int) betlist_local->occs[b1];
            value += CalcInfo.scfeigval[i+CalcInfo.num_fzc_orbs];
            ii = ioff[i] + i;
            iii = ioff[ii];

            /* loop over beta occs */
            for (b2=0; b2<b1; b2++) {
               j = (int) betlist_local->occs[b2];
               jj = ioff[j] + j;
               iijj = iii + jj;
               ij = ioff[i] + j;
               ijij = ioff[ij] + ij;
               value -= pert_param * tei[iijj];
               }
            }

         /* determine average K over spin-coupling set */
         num_unique = 0;
         for (a1=0; a1<na; a1++) unique_occs[num_unique++] = (int) alplist_local->occs[a1];
         /* for (j=0; j<num_unique; j++)
            fprintf(outfile,"unique_occs[%d] = %d\n",j,unique_occs[j]); */
            for (b1=0; b1<nb; b1++) {
               j = (int) betlist_local->occs[b1];
               for (a1=0; a1<na; a1++) {
                  if (j==unique_occs[a1]) break;
                  if (a1==(na-1)) unique_occs[num_unique++] = j;
                  }
               }
         /* fprintf(outfile,"num_unique = %d\n",num_unique);
         fprintf(outfile,"num_el = %d\n",num_el);
         */
         if (num_unique>num_el) fprintf(outfile,"WARNING: The number of explicit electrons" \
                             "!= num_el\n");

       /*
         for (j=0; j<na; j++)
            fprintf(outfile,"alp_occs[%d] = %d\n",j,(int)alplist_local->occs[j]);
         for (j=0; j<nb; j++)
            fprintf(outfile,"bet_occs[%d] = %d\n",j,(int)betlist_local->occs[j]);
         for (j=0; j<num_unique; j++)
            fprintf(outfile,"unique_occs[%d] = %d\n",j,unique_occs[j]);
       */

         Kave = 0.0;
         for (a1=0; a1<num_unique; a1++) {
            i = unique_occs[a1];
            for (b1=0; b1<a1; b1++) {
               j = unique_occs[b1];
               ij = ioff[MAX0(i,j)] + MIN0(i,j);
               ijij = ioff[ij] + ij;
               Kave += tei[ijij];
               /* fprintf(outfile,"tei[%d] = %lf\n",ijij,tei[ijij]); */
               }
            }

         /* fprintf(outfile,"num_unique = %d\n",num_unique);
         fprintf(outfile,"ioff[num_unique-1] = %d\n",ioff[num_unique]);
         fprintf(outfile,"k_total = %d\n",k_total);
         */
         if (num_unique > 1) Kave /= ioff[num_unique-1];
         value += 0.5 * Kave * k_total * pert_param;
         /* fprintf(outfile,"Kave = %lf\n",Kave); */

         if (Parameters.print_lvl > 5) {
           fprintf(outfile,"acnt = %d\t bcnt = %d\n",acnt,bcnt);
           fprintf(outfile,"tval = %lf\n",tval);
           for(a1=0; a1<na; a1++)
             fprintf(outfile," %d",alplist_local->occs[a1]);
           fprintf(outfile," \n");
           for(b1=0; b1<nb; b1++)
             fprintf(outfile," %d",betlist_local->occs[b1]);
           fprintf(outfile," \n");
           }
         H0[acnt][bcnt] = value;
       /*
         fprintf(outfile,"H0[%d][%d] = %lf\n",acnt,bcnt,value);
       */
         betlist_local++;
         } /* end loop over bcnt */

      alplist_local++;
      }

}

}} // namespace psi::detci

