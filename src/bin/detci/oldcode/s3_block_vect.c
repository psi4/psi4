#include <cstdio>
#include <libciomr/libciomr.h>
#include "structs.h"
#define EXTERN
#include "globals.h"

int form_ilist(struct stringwr *alplist, int Ja_list, int nas, int kl,
   int *L, int *R, double *Sgn);

/*
** S3_BLOCK_VECT()
**
** Calculate a block of the sigma3 vector in equation (9c) of
** Olsen, Roos, et al.  For diagonal blocks of sigma.
**
** currently assumes that (ij|ij)'s have not been halved
** Try to get the Olsen vector version working....again!!!!
*/
void s3_block_vect(struct stringwr *alplist, struct stringwr *betlist,
      double **C, double **S, double *tei, int nas, int nbs, int cnbs,
      int Ia_list, int Ja_list, int Jb_list, double **Cprime, double *F, 
      double *V, double *Sgn, int *L, int *R)
{
   struct stringwr *Ib;
   unsigned int Ib_ex;
   int ij, k, l, kl, I, J, RI;
   double tval, *CprimeI0, *CI0;
   int Ja_sym, Ia_sym;
   int ilen, Ib_idx, Jbcnt, *Iaij, *Ibij, *orbsym, norbs;
   unsigned int *Ibridx;
   signed char *Ibsgn;
   double *Tptr;

   norbs = CalcInfo.num_ci_orbs;
   orbsym = CalcInfo.orbsym + CalcInfo.num_fzc_orbs;

   /* assume fci for now */
   Ia_sym = Ia_list;
   Ja_sym = Ja_list;

   /* loop over k, l */
   for (k=0; k<norbs; k++) {
       for (l=0; l<=k; l++) {
           if ((orbsym[k] ^ orbsym[l] ^ Ja_sym ^ Ia_sym) != 0) continue;
           kl = ioff[k] + l;
           ilen = form_ilist(alplist, Ja_list, nas, kl, L, R, Sgn);
           
           if (!ilen) continue;

           Tptr = tei + ioff[kl];

           /* gather operation */
           for (I=0; I<ilen; I++) {
               CprimeI0 = Cprime[I];
               CI0 = C[L[I]];
               tval = Sgn[I];
               for (J=0; J<cnbs; J++) {
                   CprimeI0[J] = CI0[J] * tval;
               }
           }

           /* loop over Ib */
           for (Ib=betlist, Ib_idx=0; Ib_idx<nbs; Ib_idx++, Ib++) {

              zero_arr(F, cnbs);
           
               /* loop over excitations E^b_{ij} from |B(I_b)> */
               Jbcnt = Ib->cnt[Jb_list];
               Ibridx = Ib->ridx[Jb_list];
               Ibsgn = Ib->sgn[Jb_list];
               Ibij = Ib->ij[Jb_list];
               
               for (Ib_ex=0; Ib_ex < Jbcnt && (ij = *Ibij++)<=kl; Ib_ex++) {
                   J = *Ibridx++;
                   tval = *Ibsgn++;
                   if (ij == kl) tval *= 0.5;
                   F[J] += tval * Tptr[ij];
               }

               mmult(Cprime, 0, &F, 1, &V, 1, ilen, cnbs, 1, 0);

               for (I=0; I<ilen; I++) {
                   RI = R[I];
                   S[RI][Ib_idx] += V[I];
               }

           } /* end loop over Ib */

       } /* end loop over l */
   } /* end loop over k */
}              


int form_ilist(struct stringwr *alplist, int Ja_list, int nas, int kl,
   int *L, int *R, double *Sgn)
{

   int inum=0, Ia_idx, Ia_ex, Iacnt, ij;
   int *Iaij;
   struct stringwr *Ia;
   unsigned int *Iaridx;
   signed char *Iasgn;

   /* loop over Ia */
   for (Ia=alplist, Ia_idx=0; Ia_idx < nas; Ia_idx++,Ia++) {

      /* loop over excitations E^a_{kl} from |A(I_a)> */

      Iacnt = Ia->cnt[Ja_list];
      if (!Iacnt) continue;
      Iaridx = Ia->ridx[Ja_list];
      Iasgn = Ia->sgn[Ja_list];
      Iaij = Ia->ij[Ja_list];
      Ia_ex=0;
      while (Ia_ex < Iacnt && (ij = *Iaij++)<kl) Ia_ex++;
      if (ij == kl) {
         *R++ = Ia_idx;
         *L++ = Iaridx[Ia_ex];
         *Sgn++ = (double) Iasgn[Ia_ex];
         inum++;
         }
      }  /* end loop over Ia */

   return(inum);
}

