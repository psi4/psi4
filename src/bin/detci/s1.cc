/*!
  \file
  \ingroup DETCI
  \brief Enter brief description of file here 
*/

/*
** S1.C
** 
** File contains code to calculate sigma1 in various ways, all
** block-at-a-time now.
**
** C. David Sherrill
** Center for Computational Quantum Chemistry
** University of Georgia
** 
*/

#include <cstdio>
#include <libciomr/libciomr.h>
#include "structs.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace detci {

extern unsigned char ***Occs;

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))

extern void b2brepl(unsigned char **occs, int *Jcnt, int **Jij, int **Joij, 
      int **Jridx, signed char **Jsgn, struct olsen_graph *Graph,
      int Ilist, int Jlist, int len);


/*!
** S1_BLOCK_FCI(): 
**
** Calculate the sigma_1 vector as described by
** equation (20) of RAS Paper (Olsen, Roos, Jorgensen, Aa. Jensen JCP 1988)
**
** This sigma1 routine is for Full CI's only.
** currently assumes that (ij|ij)'s have not been halved!! 
** 
** David Sherrill, 21 June 1995
** Based on previous code by David Sherrill, 1994
**
** Updated 3/27/94 to include g matrix for RAS
** Modified 4/8/94 to make C and s one-dimensional
** Modified 4/10/94 to make FCI-only (for now) and use new string structs
** Modified 6/21/95 for use in new RAS program
**
** \ingroup DETCI
*/
void s1_block_fci(struct stringwr **alplist, struct stringwr **betlist, 
      double **C, double **S, double *oei, double *tei, double *F,
      int nlists, int nas, int nbs, int Ib_list, int Jb_list, 
      int Jb_list_nbs)
{
   struct stringwr *Ib, *Kb;
   unsigned int Ia_idx, Ib_idx, Kb_idx, Jb_idx;
   unsigned int Ibcnt, Kbcnt, Kb_list, Ib_ex, Kb_ex;
   unsigned int *Ibridx, *Kbridx;
   int nirreps, *Ibij, *Kbij;
   signed char *Ibsgn, *Kbsgn;
   int ij,kl,ijkl;
   double Kb_sgn, Jb_sgn;
   double tval;

   nirreps = CalcInfo.nirreps;

   /* loop over I_b */
   for (Ib=betlist[Ib_list], Ib_idx=0; Ib_idx < nbs; Ib_idx++, Ib++) {
      zero_arr(F, Jb_list_nbs);

      /* loop over excitations E^b_{kl} from |B(I_b)> */
      for (Kb_list=0; Kb_list < nlists; Kb_list++) {
         Ibcnt = Ib->cnt[Kb_list];
         Ibridx = Ib->ridx[Kb_list];
         Ibsgn = Ib->sgn[Kb_list];
         Ibij = Ib->ij[Kb_list];
         for (Ib_ex=0; Ib_ex < Ibcnt; Ib_ex++) {
            kl = *Ibij++;
            Kb_idx = *Ibridx++;
            Kb_sgn = (double) *Ibsgn++;

            /* B(K_b) = sgn(kl) * E^b_{kl} |B(I_b)> */
            Kb = betlist[Kb_list] + Kb_idx;
            if (Kb_list == Jb_list) F[Kb_idx] += Kb_sgn * oei[kl];

            /* loop over excitations E^b_{ij} from |B(K_b)> */
            /* Jb_list pre-determined because of C blocking */
            Kbcnt = Kb->cnt[Jb_list];
            Kbridx = Kb->ridx[Jb_list];
            Kbsgn = Kb->sgn[Jb_list];
            Kbij = Kb->ij[Jb_list];
            for (Kb_ex=0; Kb_ex < Kbcnt; Kb_ex++) {
               Jb_idx = *Kbridx++;
               Jb_sgn = (double) *Kbsgn++;
               ij = *Kbij++;
               ijkl = ioff[MAX0(ij,kl)] + MIN0(ij,kl) ;
               F[Jb_idx] += 0.5 * Kb_sgn * Jb_sgn * tei[ijkl] ;
               }
            } /* end loop over Ib excitations */
         } /* end loop over Kb_list */

      
      for (Ia_idx=0; Ia_idx < nas; Ia_idx++) {
         tval = 0.0;
         for (Jb_idx=0; Jb_idx < Jb_list_nbs; Jb_idx++) {
            tval += C[Ia_idx][Jb_idx] * F[Jb_idx];
            }
         S[Ia_idx][Ib_idx] += tval;
         }
      } /* end loop over Ib */

}


/*
** S1_BLOCK_RAS.C: 
** 
** Calculate the sigma_1 vector as described by
** equation (20) of RAS Paper (Olsen, Roos, Jorgensen, Aa. Jensen JCP 1988)
**
** This sigma1 routine is for RAS CI's.
** currently assumes that (ij|ij)'s have not been halved!!
** 
** David Sherrill, 2 August 1995
** Based on previous code by David Sherrill, 1994
**
** Updated 3/27/94 to include g matrix for RAS
** Modified 4/8/94 to make C and s one-dimensional
** Modified 4/10/94 to make FCI-only (for now) and use new string structs
** Modified 6/21/95 for use in new RAS program (C, s now 2D again!)
** Modified 8/2/95 to make RAS again
**
*/
void s1_block_ras(struct stringwr **alplist, struct stringwr **betlist, 
      double **C, double **S, double *oei, double *tei, double *F,
      int nlists, int nas, int nbs, int Ib_list, int Jb_list, 
      int Jb_list_nbs)
{
   struct stringwr *Ib, *Kb;
   unsigned int Ia_idx, Ib_idx, Kb_idx, Jb_idx;
   unsigned int Ibcnt, Kbcnt, Kb_list, Ib_ex, Kb_ex;
   unsigned int *Ibridx, *Kbridx;
   int nirreps,  *Ibij, *Kbij, *Iboij, *Kboij;
   signed char *Ibsgn, *Kbsgn;
   int ij,kl,ijkl,oij,okl;
   double Kb_sgn, Jb_sgn;
   double tval;

   nirreps = CalcInfo.nirreps;

   /* loop over I_b */
   for (Ib=betlist[Ib_list], Ib_idx=0; Ib_idx < nbs; Ib_idx++, Ib++) {
      zero_arr(F, Jb_list_nbs);

      /* loop over excitations E^b_{kl} from |B(I_b)> */
      for (Kb_list=0; Kb_list < nlists; Kb_list++) {
         Ibcnt = Ib->cnt[Kb_list];
         Ibridx = Ib->ridx[Kb_list];
         Ibsgn = Ib->sgn[Kb_list];
         Ibij = Ib->ij[Kb_list];
         Iboij = Ib->oij[Kb_list];
         for (Ib_ex=0; Ib_ex < Ibcnt; Ib_ex++) {
            kl = *Ibij++;
            okl = *Iboij++;
            Kb_idx = *Ibridx++;
            Kb_sgn = (double) *Ibsgn++;

            /* B(K_b) = sgn(kl) * E^b_{kl} |B(I_b)> */
            Kb = betlist[Kb_list] + Kb_idx;
            /* note okl on next line, not kl */
            if (Kb_list == Jb_list) F[Kb_idx] += Kb_sgn * oei[okl];

            /* loop over excitations E^b_{ij} from |B(K_b)> */
            /* Jb_list pre-determined because of C blocking */
            Kbcnt = Kb->cnt[Jb_list];
            Kbridx = Kb->ridx[Jb_list];
            Kbsgn = Kb->sgn[Jb_list];
            Kbij = Kb->ij[Jb_list];
            Kboij = Kb->oij[Jb_list];
            for (Kb_ex=0; Kb_ex < Kbcnt; Kb_ex++) {
               Jb_idx = *Kbridx++;
               Jb_sgn = (double) *Kbsgn++;
               ij = *Kbij++;
               oij = *Kboij++;
               ijkl = ioff[MAX0(ij,kl)] + MIN0(ij,kl) ;
               if (oij > okl) 
                  F[Jb_idx] += Kb_sgn * Jb_sgn * tei[ijkl] ;
               else if (oij == okl) 
                  F[Jb_idx] += 0.5 * Kb_sgn * Jb_sgn * tei[ijkl] ;
               }
            } /* end loop over Ib excitations */
      } /* end loop over Kb_list */

      
   for (Ia_idx=0; Ia_idx < nas; Ia_idx++) {
      tval = 0.0;
      for (Jb_idx=0; Jb_idx < Jb_list_nbs; Jb_idx++) {
         tval += C[Ia_idx][Jb_idx] * F[Jb_idx];
         }
      S[Ia_idx][Ib_idx] += tval;
      }
   } /* end loop over Ib */

}



/*
** S1_BLOCK_RAS_ROTF
** 
** String replacements on-the-fly version
**
** This sigma1 routine is for RAS CI's.
** currently assumes that (ij|ij)'s have not been halved!! 
** 
** David Sherrill, 13 August 1995
** Based on previous code by David Sherrill, 1994
**
** Updated 3/27/94 to include g matrix for RAS
** Modified 4/8/94 to make C and s one-dimensional
** Modified 4/10/94 to make FCI-only (for now) and use new string structs
** Modified 6/21/95 for use in new RAS program (C, s now 2D again!)
** Modified 8/2/95 to make RAS again
**
*/
void s1_block_ras_rotf(int *Cnt[2], int **Ij[2], int **Oij[2],
      int **Ridx[2], signed char **Sgn[2], unsigned char **Toccs,
      double **C, double **S,
      double *oei, double *tei, double *F, int nlists, int nas, int nbs,
      int Ib_list, int Jb_list, int Jb_list_nbs)
{
   int Ia_idx, Ib_idx, Kb_idx, Jb_idx;
   int Ibcnt, Kbcnt, Kb_list, Ib_ex, Kb_ex;
   int *Ibridx, *Kbridx;
   int nirreps,  *Ibij, *Kbij, *Iboij, *Kboij;
   signed char *Ibsgn, *Kbsgn;
   int i,ij,kl,ijkl,oij,okl;
   double Kb_sgn, Jb_sgn;
   double tval;

   nirreps = CalcInfo.nirreps;

   for (Kb_list=0; Kb_list < nlists; Kb_list++) {
      b2brepl(Occs[Ib_list], Cnt[0], Ij[0], Oij[0], Ridx[0],
         Sgn[0], BetaG, Ib_list, Kb_list, nbs);

      /* loop over I_b */
      for (Ib_idx=0; Ib_idx < nbs; Ib_idx++) {

         if ((Ibcnt = Cnt[0][Ib_idx]) < 0) continue;
         zero_arr(F, Jb_list_nbs);

         /* loop over excitations E^b_{kl} from |B(I_b)> */
         Ibridx = Ridx[0][Ib_idx];
         Ibsgn = Sgn[0][Ib_idx];
         Ibij = Ij[0][Ib_idx];
         Iboij = Oij[0][Ib_idx];

         for (i=0; i<Ibcnt; i++) 
            Toccs[i] = Occs[Kb_list][Ibridx[i]];

         b2brepl(Toccs, Cnt[1], Ij[1], Oij[1], Ridx[1], Sgn[1],
            BetaG, Kb_list, Jb_list, Ibcnt);

         for (Ib_ex=0; Ib_ex < Ibcnt; Ib_ex++) {
            kl = *Ibij++;
            okl = *Iboij++;
            Kb_idx = *Ibridx++;
            Kb_sgn = (double) *Ibsgn++;

            /* B(K_b) = sgn(kl) * E^b_{kl} |B(I_b)> */
            /* note okl on next line, not kl */
            if (Kb_list == Jb_list) F[Kb_idx] += Kb_sgn * oei[okl];

            /* loop over excitations E^b_{ij} from |B(K_b)> */
            /* Jb_list pre-determined because of C blocking */
            Kbcnt = Cnt[1][Ib_ex];
            Kbridx = Ridx[1][Ib_ex];
            Kbsgn = Sgn[1][Ib_ex];
            Kbij = Ij[1][Ib_ex];
            Kboij = Oij[1][Ib_ex];
            for (Kb_ex=0; Kb_ex < Kbcnt; Kb_ex++) {
               Jb_idx = *Kbridx++;
               Jb_sgn = (double) *Kbsgn++;
               ij = *Kbij++;
               oij = *Kboij++;
               ijkl = ioff[MAX0(ij,kl)] + MIN0(ij,kl) ;
               if (oij > okl) 
                  F[Jb_idx] += Kb_sgn * Jb_sgn * tei[ijkl] ;
               else if (oij == okl) 
                  F[Jb_idx] += 0.5 * Kb_sgn * Jb_sgn * tei[ijkl] ;
               }
            } /* end loop over Ib excitations */

      for (Ia_idx=0; Ia_idx < nas; Ia_idx++) {
         tval = 0.0;
         for (Jb_idx=0; Jb_idx < Jb_list_nbs; Jb_idx++) {
            tval += C[Ia_idx][Jb_idx] * F[Jb_idx];
            }
         S[Ia_idx][Ib_idx] += tval;
         }
      } /* end loop over Ib */
   } /* end loop over Kb_list */

}

}} // namespace psi::detci


