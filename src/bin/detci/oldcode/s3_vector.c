/*
** CALC_SIGMA3.C: This file contains the function needed to calculate
**   the sigma1 contribution to sigma = H c, as given by Olsen, Roos,
**   Jorgensen, and Aa. Jensen in equation (9c).  
** 
** This version uses a scatter/gather idea in a clever way that avoids
**   copying the C vector into the C' vector.  The C I(alpha) pointers
**   are merely rearranged to make C' ("virtual" scatter/gather?).  
**   This idea is similar to those alluded to by Mitrushenkov.
**   
** David Sherrill, 18 November 1994
**
*/

#include <cstdio>
#include <libciomr/libciomr.h>
#include "structs.h"
#include "CalcInfo.h"
#include "Parameters.h"
#include "globals.h"


#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))

/* prototypes within this module */
int form_ilist(int Iasym, struct stringwr *slist, unsigned int *L, 
   unsigned int *R, int *S, int Jasym, int klsym, int kl);


/*
** calc_sigma3(): Calculate the sigma3 vector in equation (9c) of
**    Olsen, Roos, et. al.
**
** Modified 4/8/94 to make C and s one-dimensional
** Modified 11/18/94 for virtual scatter/gather method
** Warning to C neophytes: C is case-sensitive! (e.g. s is not S)
**
*/
void calc_sigma3(struct stringwr *slist, double **C, double **s,
      double *tei, int nas, int nbs, 
      unsigned int *bfora, unsigned int *bfirst)
{

   struct stringwr *Ib ;
   int Ia_sym, Ja_sym, Ib_sym, Jb_sym, Ib_idx, Jb_idx;
   int Ib_offset, Ib_end, Jb_end, *Ibij;
   int k, l, ij, kl, klsym;
   int ioffk, Inum, I;
   unsigned int Jbcnt, *Ibridx, Ib_ex;
   signed char *Ibsgn;
   double Jb_sgn;

   static unsigned int *L, LI, *L0 = NULL;
   static unsigned int *R, RI, *R0 = NULL;
   static int *S, *S0 = NULL;  /* hmmm... */
   double tsgn;
   static double *F0 = NULL;
   static double *V = NULL;
   static double **Cprime = NULL;
   double *Tptr;
   int *orbsym;

   orbsym = CalcInfo.orbsym + CalcInfo.num_fzc_orbs;

   if (F0 == NULL) F0 = init_array(nbs);
   if (Cprime == NULL) {
      Cprime = init_matrix(nas, nbs);
      }

   if (V == NULL) V = init_array(nas);
   if (L0 == NULL) L0 = (unsigned int *) malloc(nas * sizeof(unsigned int));
   if (R0 == NULL) R0 = (unsigned int *) malloc(nas * sizeof(unsigned int));
   if (S0 == NULL) S0 = (int *) malloc (nas * sizeof(int));

   /* set up list L(I), R(I), and S(I) */
   for (Ia_sym=0; Ia_sym < CalcInfo.nirreps; Ia_sym++) {
      for (Ja_sym=0; Ja_sym < CalcInfo.nirreps; Ja_sym++) {
         Jb_sym = CalcInfo.ref_sym ^ Ja_sym;
         Jb_end = CalcInfo.bsymnum[Jb_sym];
 
         for (k=0; k<CalcInfo.num_ci_orbs; k++) {
            ioffk = ioff[k];
            for (l=0,kl=ioffk; l<=k; l++,kl++) {
            klsym = orbsym[k] ^ orbsym[l]; 
            Tptr = tei + ioff[kl];
            Inum = form_ilist(Ia_sym, slist, L0, R0, S0, Ja_sym, klsym, kl);
            if (!Inum) continue;

            /* gathering operation to form C' */
            for (I=0,L=L0,S=S0; I<Inum; I++) {
               LI = *L++;
               tsgn = *S++;               
               for (Jb_idx=0; Jb_idx < Jb_end; Jb_idx++) {
                  Cprime[I][Jb_idx] = C[LI][Jb_idx] * tsgn;
                  }
               }

            /* loop over Ib */
            Ib_sym = CalcInfo.ref_sym ^ Ia_sym;
            Ib_offset = CalcInfo.bsymst[Ib_sym];
            Ib_end = CalcInfo.bsymnum[Ib_sym];
            for (Ib_idx=0,Ib=slist+Ib_offset; Ib_idx<Ib_end; Ib_idx++,Ib++) {
               zero_arr(F0, nbs);
                
               /* loop over excitations E^{b}_{ij} from |B(I_b)> */
               Jbcnt = Ib->cnt[Jb_sym][klsym];
               Ibridx = Ib->ridx[Jb_sym][klsym];
               Ibsgn = Ib->sgn[Jb_sym][klsym];
               Ibij = Ib->ij[Jb_sym][klsym];
 
               for (Ib_ex=0; Ib_ex < Jbcnt && (ij = *Ibij++)<=kl; Ib_ex++) {
                  Jb_idx = *Ibridx++;
                  Jb_sgn = (double) *Ibsgn++;
                  F0[Jb_idx] += Jb_sgn * Tptr[ij];
                  }

               /* V(I) = \Sum{J_b} F(J_b) * C'(I, J_b) */
               mmult(Cprime, 0, &F0, 1, &V, 1, Inum, Jb_end, 1, 0);

               /* vectorized scattering */
               R = R0;
               for (I=0,R=R0; I<Inum; I++) {
                  RI = *R++;
                  s[RI][Ib_idx] += V[I];
                  }

               } /* end loop over Ib */

            } /* end loop over l */
         } /* end loop over k */
      } /* end loop over Ja_sym */
   } /* end loop over Ia_sym */

}


/*
** form_ilist()
**
** David Sherrill, 22 November 1994
**
*/
int form_ilist(int Iasym, struct stringwr *slist, unsigned int *L, 
   unsigned int *R, int *S, int Jasym, int klsym, int kl)
{
   struct stringwr *Ia;
   unsigned int Ia_idx, offset, *Iaridx, Ia_ex;
   int Iacnt, *Iaij;
   int Inum=0;
   signed char *Iasgn; 

      offset = CalcInfo.asymst[Iasym];
      for (Ia_idx = offset, Ia = slist + offset; Ia_idx < 
           offset + CalcInfo.asymnum[Iasym]; Ia_idx++, Ia++) {
         Iacnt = Ia->cnt[Jasym][klsym];
         Iaridx = Ia->ridx[Jasym][klsym];
         Iasgn = Ia->sgn[Jasym][klsym];
         Iaij = Ia->ij[Jasym][klsym];
         for (Ia_ex=0; Ia_ex < Iacnt; Ia_ex++) {
            if (Iaij[Ia_ex] == kl) {
               *R++ = Ia_idx;
               *L++ = Iaridx[Ia_ex] + CalcInfo.asymst[Jasym]; 
               *S++ = Iasgn[Ia_ex];
               Inum++;
               } 
            } 
         } /* end construction of I lists */

   return(Inum);

}


