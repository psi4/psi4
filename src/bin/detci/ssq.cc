/*! \file
**  \ingroup DETCI
**  \brief Compute expectation value of S^2
**
** Routine for computing the expectation value of S^2.
** Useful for determining if spin-contamination (due to the davidson
** procedure) is a problem.
**
** 24 June 1997
**
*/

/* #define DEBUG */
#include <cstdio>
#include <libciomr/libciomr.h>
#include "structs.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace detci {

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))


/*
** SSQ()
**
** Calculates the expectation value of S^2.
**
*/
double ssq(struct stringwr *alplist, struct stringwr *betlist,
     double **CL, double **CR, int nas, int nbs,
     int Ja_list, int Jb_list) 
{
   struct stringwr *Ia, *Ib ;
   unsigned int Ia_ex, Ib_ex;
   int Ia_idx, Ib_idx;
   int Ja_idx, Jb_idx;
   int Ja_sgn, Jb_sgn;
   int ij, ji, i1, j1, i2, j2;
   double tval, Ms, S2, smin_spls = 0.0;

   int Iacnt, Jbcnt, *Iaij, *Ibij;
   unsigned int *Iaridx, *Ibridx;
   signed char *Iasgn, *Ibsgn;

   /* <S^2> = <S_z> + <S_z>^2 + <S_S+> */
   /* First determine the expection value of <S_S+> */

   /* loop over Ia */
   #ifdef DEBUG
   fprintf(outfile,"number of alpha strings = %d\n",nas);
   #endif
   for (Ia=alplist,Ia_idx=0; Ia_idx < nas; Ia_idx++,Ia++) {

      /* loop over excitations E^a_{ji} from |A(I_a)> */
      Iacnt = Ia->cnt[Ja_list];
      Iaridx = Ia->ridx[Ja_list];
      Iasgn = Ia->sgn[Ja_list];
      Iaij = Ia->oij[Ja_list];
      for (Ia_ex=0; Ia_ex < Iacnt; Ia_ex++) {
         ji = *Iaij++;
         Ja_idx = *Iaridx++;
         Ja_sgn = *Iasgn++;
         i1 = ji/CalcInfo.num_ci_orbs;
         j1 = ji%CalcInfo.num_ci_orbs;

         /* loop over Ib */
         #ifdef DEBUG
         fprintf(outfile,"number of beta strings = %d\n",nbs);
         #endif
         for (Ib=betlist, Ib_idx=0; Ib_idx < nbs; Ib_idx++, Ib++) {

            /* loop over excitations E^b_{ij} from |B(I_b)> */
            Jbcnt = Ib->cnt[Jb_list];
            Ibridx = Ib->ridx[Jb_list];
            Ibsgn = Ib->sgn[Jb_list];
            Ibij = Ib->oij[Jb_list];
 
            tval = 0.0;
            for (Ib_ex=0; Ib_ex < Jbcnt; Ib_ex++) {
               ij = *Ibij++;
               Jb_idx = *Ibridx++;
               Jb_sgn = *Ibsgn++;
               i2 = ij/CalcInfo.num_ci_orbs;
               j2 = ij%CalcInfo.num_ci_orbs; 
               if (i1!=j2 || i2!=j1) continue;
               tval += CR[Ia_idx][Ib_idx] * CL[Ja_idx][Jb_idx] *
                   (double) Ja_sgn * (double) Jb_sgn;
               #ifdef DEBUG
               fprintf(outfile,"\n\nIa_idx = %d\n",Ia_idx);
               fprintf(outfile,"Ib_idx = %d\n",Ib_idx);
               fprintf(outfile,"Ja_idx = %d\n",Ja_idx);
               fprintf(outfile,"Jb_idx = %d\n",Jb_idx);
               fprintf(outfile,"tval_ssq = %lf\n",-tval);
               fprintf(outfile,"CR = %lf\n",CR[Ia_idx][Ib_idx]);
               fprintf(outfile,"LR = %lf\n",CL[Ja_idx][Jb_idx]);
               fprintf(outfile,"Ja_sgn = %lf\n",Ja_sgn);
               fprintf(outfile,"Jb_sgn = %lf\n",Jb_sgn);
               #endif
               }
            smin_spls += tval;
      
            } /* end loop over Ib */
         } /* end loop over Ia excitations */ 
     } /* end loop over Ia */ 

   S2 = -smin_spls;

   return(S2);
}

}} // namespace psi::detci

