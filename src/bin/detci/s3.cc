/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here 
*/

/*
** S3.C
** 
** Routines for calculating a block of the sigma3 (alpha-beta) vector
**
** C. David Sherrill
** Center for Computational Quantum Chemistry
** University of Georgia
** 21 June 1995
**
*/

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
** S3_BLOCK()
**
** Calculate a block of the sigma3 vector in equation (9c) of
** Olsen, Roos, et al.
**
*/
void s3_block(struct stringwr *alplist, struct stringwr *betlist,
      double **C, double **S, double *tei, int nas, int nbs,
      int Ja_list, int Jb_list)
{
   struct stringwr *Ia, *Ib ;
   unsigned int Ia_ex, Ib_ex;
   int Ia_idx, Ib_idx ;
   int Ja_idx, Jb_idx ;
   int Ja_sgn, Jb_sgn ;
   int ij, kl, ijkl ;
   double tval ;

   int Iacnt, Jbcnt, *Iaij, *Ibij;
   unsigned int *Iaridx, *Ibridx;
   signed char *Iasgn, *Ibsgn;
   int nirreps;
   double *Tptr;

   nirreps = CalcInfo.nirreps;

   /* loop over Ia */
   for (Ia=alplist,Ia_idx=0; Ia_idx < nas; Ia_idx++,Ia++) {

      /* loop over excitations E^a_{kl} from |A(I_a)> */
      Iacnt = Ia->cnt[Ja_list];
      Iaridx = Ia->ridx[Ja_list];
      Iasgn = Ia->sgn[Ja_list];
      Iaij = Ia->ij[Ja_list];
      for (Ia_ex=0; Ia_ex < Iacnt; Ia_ex++) {
         kl = *Iaij++;
         Tptr = tei + ioff[kl];
         Ja_idx = *Iaridx++;
         Ja_sgn = *Iasgn++;

         /* loop over Ib */
         for (Ib=betlist, Ib_idx=0; Ib_idx < nbs; Ib_idx++, Ib++) {

            /* loop over excitations E^b_{ij} from |B(I_b)> */
            Jbcnt = Ib->cnt[Jb_list];
            Ibridx = Ib->ridx[Jb_list];
            Ibsgn = Ib->sgn[Jb_list];
            Ibij = Ib->ij[Jb_list];
                  
            tval = 0.0;
            for (Ib_ex=0; Ib_ex < Jbcnt; Ib_ex++) {
               ij = *Ibij++;
               Jb_idx = *Ibridx++;
               Jb_sgn = *Ibsgn++;
               ijkl = ioff[MAX0(ij,kl)] + MIN0(ij,kl);
               tval += tei[ijkl] * C[Ja_idx][Jb_idx] * 
                  (double) Ja_sgn * (double) Jb_sgn;
               }
            S[Ia_idx][Ib_idx] += tval;

            } /* end loop over Ib */
         } /* end loop over Ia excitations */
   } /* end loop over Ia */
}



/*
** S3_BLOCK_DIAG()
**
** Calculate a block of the sigma3 vector in equation (9c) of
** Olsen, Roos, et al.  For diagonal blocks of sigma.
**
** currently assumes that (ij|ij)'s have not been halved
**
*/
void s3_block_diag(struct stringwr *alplist, struct stringwr *betlist,
      double **C, double **S, double *tei, int nas, int nbs,
      int Ja_list, int Jb_list)
{
   struct stringwr *Ia, *Ib;
   unsigned int Ia_ex, Ib_ex;
   int Ia_idx, Ib_idx;
   int Ja_idx, Jb_idx;
   int Ja_sgn, Jb_sgn;
   int ij, kl;
   double tval,tval2;

   int Iacnt, Jbcnt, *Iaij, *Ibij;
   unsigned int *Iaridx, *Ibridx;
   signed char *Iasgn, *Ibsgn;
   double *Tptr, *Cptr, *Sptr;

   /* loop over Ia */
   for (Ia=alplist,Ia_idx=0; Ia_idx < nas; Ia_idx++,Ia++) {

      /* loop over excitations E^a_{kl} from |A(I_a)> */

      Iacnt = Ia->cnt[Ja_list];
      Iaridx = Ia->ridx[Ja_list];
      Iasgn = Ia->sgn[Ja_list];
      Iaij = Ia->ij[Ja_list];
      for (Ia_ex=0; Ia_ex < Iacnt; Ia_ex++) {

         Sptr = S[Ia_idx];

         kl = *Iaij++;
         Tptr = tei + ioff[kl];
         Ja_idx = *Iaridx++;
         Cptr = C[Ja_idx];
         Ja_sgn = *Iasgn++;

         /* loop over Ib */
         for (Ib=betlist, Ib_idx=0; Ib_idx < nbs; Ib_idx++, Ib++) {

            /* loop over excitations E^b_{ij} from |B(I_b)> */
            Jbcnt = Ib->cnt[Jb_list];
            Ibridx = Ib->ridx[Jb_list];
            Ibsgn = Ib->sgn[Jb_list];
            Ibij = Ib->ij[Jb_list];
               
            tval = 0.0;
            for (Ib_ex=0; Ib_ex < Jbcnt && (ij = *Ibij++)<=kl; Ib_ex++) {
               Jb_idx = *Ibridx++;
               Jb_sgn = *Ibsgn++ * Ja_sgn;
               tval2 = Tptr[ij] * Cptr[Jb_idx] * Jb_sgn;
               if (ij == kl) tval2 *= 0.5; 
               tval += tval2;
               }
            *Sptr++ += tval;
            } /* end loop over Ib */
         } /* end loop over Ia excitations */
      } /* end loop over Ia */
}




/*
** S3_BLOCK_ROTF()
**
** Calculate a block of the sigma3 (alpha-beta) vector for a diagonal
**   block of sigma.  The string replacements are fed in through arrays.
**   Assumes that (ij|ij)'s have not been halved
**
** David Sherrill, 13 August 1995
**
*/
void s3_block_rotf(int *Cnt[2], int **Ij[2], 
      int **Ridx[2], signed char **Sgn[2], double **C, double **S,
      double *tei, int nas, int nbs)
{
   int Ia_ex, Ib_ex;
   int Ia_idx, Ib_idx ;
   int Ja_idx, Jb_idx ;
   int Ja_sgn, Jb_sgn ;
   int ij, kl, ijkl ;
   double tval;

   int Iacnt, Jbcnt, *Iaij, *Ibij;
   int *Iaridx, *Ibridx;
   signed char *Iasgn, *Ibsgn;
   double *Tptr;


   /* loop over Ia */
   for (Ia_idx=0; Ia_idx < nas; Ia_idx++) {

      /* loop over excitations E^a_{kl} from |A(I_a)> */
      Iacnt = Cnt[0][Ia_idx];
      Iaridx = Ridx[0][Ia_idx];
      Iasgn = Sgn[0][Ia_idx];
      Iaij = Ij[0][Ia_idx];
      for (Ia_ex=0; Ia_ex < Iacnt; Ia_ex++) {
         kl = *Iaij++;
         Tptr = tei + ioff[kl];
         Ja_idx = *Iaridx++;
         Ja_sgn = *Iasgn++;

         /* loop over Ib */
         for (Ib_idx=0; Ib_idx < nbs; Ib_idx++) {

            /* loop over excitations E^b_{ij} from |B(I_b)> */
            Jbcnt = Cnt[1][Ib_idx];
            Ibridx = Ridx[1][Ib_idx];
            Ibsgn = Sgn[1][Ib_idx];
            Ibij = Ij[1][Ib_idx];
               
            tval = 0.0;
            for (Ib_ex=0; Ib_ex < Jbcnt; Ib_ex++) {
               ij = *Ibij++;
               Jb_idx = *Ibridx++;
               Jb_sgn = *Ibsgn++;
               ijkl = ioff[MAX0(ij,kl)] + MIN0(ij,kl);
               tval += tei[ijkl] * C[Ja_idx][Jb_idx] * 
                  (double) Ja_sgn * (double) Jb_sgn;
               }
            S[Ia_idx][Ib_idx] += tval;
            } /* end loop over Ib */
         } /* end loop over Ia excitations */
      } /* end loop over Ia */
}


/*
** S3_BLOCK_DIAG_ROTF()
**
** Calculate a block of the sigma3 vector in equation 
** (9c) of Olsen, Roos, et al.  For diagonal blocks of sigma.
** The string replacements are fed in through arrays.  
** Assumes that (ij|ij)'s have not been halved.
**
*/
void s3_block_diag_rotf(int *Cnt[2], int **Ij[2], 
      int **Ridx[2], signed char **Sgn[2], double **C, double **S,
      double *tei, int nas, int nbs)
{
   int Ia_ex, Ib_ex;
   int Ia_idx, Ib_idx ;
   int Ja_idx, Jb_idx ;
   int Ja_sgn, Jb_sgn ;
   int ij, kl, ijkl ;
   double tval,tval2 ;

   int Iacnt, Jbcnt, *Iaij, *Ibij;
   int *Iaridx, *Ibridx;
   signed char *Iasgn, *Ibsgn;
   double *Tptr;


   /* loop over Ia */
   for (Ia_idx=0; Ia_idx < nas; Ia_idx++) {

      /* loop over excitations E^a_{kl} from |A(I_a)> */
      Iacnt = Cnt[0][Ia_idx];
      Iaridx = Ridx[0][Ia_idx];
      Iasgn = Sgn[0][Ia_idx];
      Iaij = Ij[0][Ia_idx];
      for (Ia_ex=0; Ia_ex < Iacnt; Ia_ex++) {
         kl = *Iaij++;
         Tptr = tei + ioff[kl];
         Ja_idx = *Iaridx++;
         Ja_sgn = *Iasgn++;

         /* loop over Ib */
         for (Ib_idx=0; Ib_idx < nbs; Ib_idx++) {

            /* loop over excitations E^b_{ij} from |B(I_b)> */
            Jbcnt = Cnt[1][Ib_idx];
            Ibridx = Ridx[1][Ib_idx];
            Ibsgn = Sgn[1][Ib_idx];
            Ibij = Ij[1][Ib_idx];
               
            tval = 0.0;
            for (Ib_ex=0; Ib_ex < Jbcnt; Ib_ex++) {
               ij = *Ibij++;
               Jb_idx = *Ibridx++;
               Jb_sgn = *Ibsgn++;
               if (ij > kl) continue;
               ijkl = ioff[MAX0(ij,kl)] + MIN0(ij,kl);
               tval2 = tei[ijkl] * C[Ja_idx][Jb_idx] * 
                  (double) Ja_sgn * (double) Jb_sgn;
               if (ij == kl) tval2 *= 0.5;
               tval += tval2;
               }
            S[Ia_idx][Ib_idx] += tval;
            } /* end loop over Ib */
         } /* end loop over Ia excitations */
      } /* end loop over Ia */
}


}} // namespace psi::detci

