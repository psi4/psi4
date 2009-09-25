/*! \file
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
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "structs.h"
#define EXTERN
#include "globals.h"
#include <pthread.h>
#include "tpool.h"

namespace psi { namespace detci {

extern unsigned char ***Occs;

extern void b2brepl(unsigned char **occs, int *Jcnt, int **Jij, int **Joij, 
      int **Jridx, signed char **Jsgn, struct olsen_graph *Graph,
      int Ilist, int Jlist, int len);
void s1_block_vfci_pthread(void *threadarg);
void s1_block_vras_pthread(void *threadarg);

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))


/*
** S1_BLOCK_VFCI(): 
**
** Calculate the sigma_1 vector as described by
** equation (20) of RAS Paper (Olsen, Roos, Jorgensen, Aa. Jensen JCP 1988)
**
** This sigma1 routine is for Full CI's only, assumes (ij|ij)'s have not
** been halved, and attempts to follow Olsen's vectorized algorithm more
** closely than previous versions, using sparsity of F.
** 
** David Sherrill, 10 May 1996
** Based on previous code by David Sherrill, 1994-5
**
** Updated 3/27/94 to include g matrix for RAS
** Modified 4/8/94 to make C and s one-dimensional
** Modified 4/10/94 to make FCI-only (for now) and use new string structs
** Modified 6/21/95 for use in new RAS program
** Modified 5/10/96 for new sparse-F method
*/
void s1_block_vfci(struct stringwr **alplist, struct stringwr **betlist, 
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
               ijkl = INDEX(ij,kl);
               F[Jb_idx] += 0.5 * Kb_sgn * Jb_sgn * tei[ijkl] ;
               }
            } /* end loop over Ib excitations */
         } /* end loop over Kb_list */

      
      /*
      for (Ia_idx=0; Ia_idx < nas; Ia_idx++) {
         tval = 0.0;
         for (Jb_idx=0; Jb_idx < Jb_list_nbs; Jb_idx++) {
            tval += C[Ia_idx][Jb_idx] * F[Jb_idx];
            }
         S[Ia_idx][Ib_idx] += tval;
         }
      */

      /* need to improve mem access pattern here! Above vers may be better! */
      /* min op cnt may also be better */
      for (Jb_idx=0; Jb_idx < Jb_list_nbs; Jb_idx++) {
         if ((tval=F[Jb_idx]) == 0.0) continue;

      #ifdef USE_BLAS
         C_DAXPY(nas,tval,(C[0]+Jb_idx),Jb_list_nbs,(S[0]+Ib_idx),nbs);
      #else
         for (Ia_idx=0; Ia_idx < nas; Ia_idx++) {
            S[Ia_idx][Ib_idx] += tval * C[Ia_idx][Jb_idx];
            }
      #endif
         }

      } /* end loop over Ib */

}

/*
** S1_BLOCK_VFCI_THREAD(): 
**
** Calculate the sigma_1 vector as described by
** equation (20) of RAS Paper (Olsen, Roos, Jorgensen, Aa. Jensen JCP 1988)
**
** This sigma1 routine is for Full CI's only, assumes (ij|ij)'s have not
** been halved, and attempts to follow Olsen's vectorized algorithm more
** closely than previous versions, using sparsity of F.
** 
** David Sherrill, 10 May 1996
** Based on previous code by David Sherrill, 1994-5
**
** Updated 3/27/94 to include g matrix for RAS
** Modified 4/8/94 to make C and s one-dimensional
** Modified 4/10/94 to make FCI-only (for now) and use new string structs
** Modified 6/21/95 for use in new RAS program
** Modified 5/10/96 for new sparse-F method
*/
void s1_block_vfci_thread(struct stringwr **alplist, struct stringwr **betlist, 
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
  struct pthreads_s1vfci **thread_info;
  int i;
   
  thread_info = (struct pthreads_s1vfci **)
                malloc(sizeof(struct pthreads_s1vfci *) * nbs);
  for (i=0; i<nbs; i++) {
      thread_info[i] = (struct pthreads_s1vfci *)
                       malloc(sizeof(struct pthreads_s1vfci));
    }

  tpool_queue_open(thread_pool);
   

  detci_time.s1_mt_before_time = wall_time_new();

  /* loop over I_b */
  for (Ib=betlist[Ib_list], Ib_idx=0; Ib_idx < nbs; Ib_idx++, Ib++) {
      thread_info[Ib_idx]->alplist=alplist;
      thread_info[Ib_idx]->betlist=betlist;
      thread_info[Ib_idx]->C=C;
      thread_info[Ib_idx]->S=S;
      thread_info[Ib_idx]->oei=oei;
      thread_info[Ib_idx]->tei=tei;
      thread_info[Ib_idx]->nlists=nlists;
      thread_info[Ib_idx]->nas=nas;
      thread_info[Ib_idx]->nbs=nbs;
      thread_info[Ib_idx]->Ib_list=Ib_list;
      thread_info[Ib_idx]->Jb_list=Jb_list;
      thread_info[Ib_idx]->Jb_list_nbs=Jb_list_nbs;
      thread_info[Ib_idx]->Ib=Ib;
      thread_info[Ib_idx]->Ib_idx=Ib_idx;
      tpool_add_work(thread_pool, s1_block_vfci_pthread, (void *) thread_info[Ib_idx]);
    } /* end loop over Ib */
  tpool_queue_close(thread_pool, 1);

  detci_time.s1_mt_after_time = wall_time_new();
  detci_time.s1_mt_total_time += detci_time.s1_mt_after_time - detci_time.s1_mt_before_time;

  for (i=0; i<nbs; i++) free(thread_info[i]);

   
}

/*
** S1_BLOCK_VFCI_PTHREAD(): 
**
** Calculate the sigma_1 vector as described by
** equation (20) of RAS Paper (Olsen, Roos, Jorgensen, Aa. Jensen JCP 1988)
**
** This sigma1 routine is for Full CI's only, assumes (ij|ij)'s have not
** been halved, and attempts to follow Olsen's vectorized algorithm more
** closely than previous versions, using sparsity of F.
** 
** David Sherrill, 10 May 1996
** Based on previous code by David Sherrill, 1994-5
**
** Updated 3/27/94 to include g matrix for RAS
** Modified 4/8/94 to make C and s one-dimensional
** Modified 4/10/94 to make FCI-only (for now) and use new string structs
** Modified 6/21/95 for use in new RAS program
** Modified 5/10/96 for new sparse-F method
*/
void s1_block_vfci_pthread(void *threadarg)
{
  struct stringwr *Ib, *Kb, **alplist, **betlist;
  unsigned int Ia_idx, Ib_idx, Kb_idx, Jb_idx;
  unsigned int Ibcnt, Kbcnt, Kb_list, Ib_ex, Kb_ex;
  unsigned int *Ibridx, *Kbridx;
  int nirreps, *Ibij, *Kbij;
  signed char *Ibsgn, *Kbsgn;
  int ij,kl,ijkl;
  double Kb_sgn, Jb_sgn;
  double tval;
  double *oei, *tei, **C, **S, *F;
  struct pthreads_s1vfci *thread_info;
  int nlists, nas, nbs, Ib_list, Jb_list, Jb_list_nbs;
   
  thread_info = (struct pthreads_s1vfci *) threadarg;
  alplist = thread_info->alplist;
  betlist = thread_info->betlist;
  C = thread_info->C;
  S = thread_info->S;
  oei = thread_info->oei;
  tei = thread_info->tei;
  nlists = thread_info->nlists;
  nas = thread_info->nas;
  nbs = thread_info->nbs;
  Ib_list = thread_info->Ib_list;
  Jb_list = thread_info->Jb_list;
  Jb_list_nbs = thread_info->Jb_list_nbs;
  Ib = thread_info->Ib;
  Ib_idx = thread_info->Ib_idx;
  
  nirreps = CalcInfo.nirreps;
  F = init_array(Jb_list_nbs);
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
              ijkl = INDEX(ij,kl);
              F[Jb_idx] += 0.5 * Kb_sgn * Jb_sgn * tei[ijkl] ;
            }
        } /* end loop over Ib excitations */
    } /* end loop over Kb_list */

      
  /*
    for (Ia_idx=0; Ia_idx < nas; Ia_idx++) {
    tval = 0.0;
    for (Jb_idx=0; Jb_idx < Jb_list_nbs; Jb_idx++) {
    tval += C[Ia_idx][Jb_idx] * F[Jb_idx];
    }
    S[Ia_idx][Ib_idx] += tval;
    }
  */

  /* need to improve mem access pattern here! Above vers may be better! */
  /* min op cnt may also be better */
  for (Jb_idx=0; Jb_idx < Jb_list_nbs; Jb_idx++) {
      if ((tval=F[Jb_idx]) == 0.0) continue;

#ifdef USE_BLAS
      C_DAXPY(nas,tval,(C[0]+Jb_idx),Jb_list_nbs,(S[0]+Ib_idx),nbs);
#else
      for (Ia_idx=0; Ia_idx < nas; Ia_idx++) {
          S[Ia_idx][Ib_idx] += tval * C[Ia_idx][Jb_idx];
        }
#endif
    }
  free(F);

}


/*
** S1_BLOCK_VRAS.C: 
** 
** Calculate the sigma_1 vector as described by
** equation (20) of RAS Paper (Olsen, Roos, Jorgensen, Aa. Jensen JCP 1988)
**
** This sigma1 routine is for Full CI's only, assumes (ij|ij)'s have not
** been halved, and attempts to follow Olsen's vectorized algorithm more
** closely than previous versions, using sparsity of F.
** 
** David Sherrill, 10 May 1996
** Based on previous code by David Sherrill, 1994-5
**
** Updated 3/27/94 to include g matrix for RAS
** Modified 4/8/94 to make C and s one-dimensional
** Modified 4/10/94 to make FCI-only (for now) and use new string structs
** Modified 6/21/95 for use in new RAS program (C, s now 2D again!)
** Modified 8/2/95 to make RAS again
** Modified 5/10/96 for new sparse-F method
*/
void s1_block_vras(struct stringwr **alplist, struct stringwr **betlist, 
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
               ijkl = INDEX(ij,kl);
               if (oij > okl) 
                  F[Jb_idx] += Kb_sgn * Jb_sgn * tei[ijkl] ;
               else if (oij == okl) 
                  F[Jb_idx] += 0.5 * Kb_sgn * Jb_sgn * tei[ijkl] ;
               }
            } /* end loop over Ib excitations */
      } /* end loop over Kb_list */

      
   /* 
   for (Ia_idx=0; Ia_idx < nas; Ia_idx++) {
      tval = 0.0;
      for (Jb_idx=0; Jb_idx < Jb_list_nbs; Jb_idx++) {
         tval += C[Ia_idx][Jb_idx] * F[Jb_idx];
         }
      S[Ia_idx][Ib_idx] += tval;
      }
   */

      /* need to improve mem access pattern here! Above vers may be better!  */
      /* min op cnt may also be better */
      for (Jb_idx=0; Jb_idx < Jb_list_nbs; Jb_idx++) {
         if ((tval=F[Jb_idx]) == 0.0) continue;

      #ifdef USE_BLAS
         C_DAXPY(nas,tval, (C[0]+Jb_idx), Jb_list_nbs, (S[0]+Ib_idx), nbs);
      #else
         for (Ia_idx=0; Ia_idx < nas; Ia_idx++) {
            S[Ia_idx][Ib_idx] += tval * C[Ia_idx][Jb_idx];
            }
      #endif
         }


   } /* end loop over Ib */

}


/*
** S1_BLOCK_VRAS_THREAD.C: 
** 
** Calculate the sigma_1 vector as described by
** equation (20) of RAS Paper (Olsen, Roos, Jorgensen, Aa. Jensen JCP 1988)
**
** This sigma1 routine is for Full CI's only, assumes (ij|ij)'s have not
** been halved, and attempts to follow Olsen's vectorized algorithm more
** closely than previous versions, using sparsity of F.
** 
** David Sherrill, 10 May 1996
** Based on previous code by David Sherrill, 1994-5
**
** Updated 3/27/94 to include g matrix for RAS
** Modified 4/8/94 to make C and s one-dimensional
** Modified 4/10/94 to make FCI-only (for now) and use new string structs
** Modified 6/21/95 for use in new RAS program (C, s now 2D again!)
** Modified 8/2/95 to make RAS again
** Modified 5/10/96 for new sparse-F method
*/
void s1_block_vras_thread(struct stringwr **alplist, struct stringwr **betlist, 
      double **C, double **S, double *oei, double *tei, double *F,
      int nlists, int nas, int nbs, int Ib_list, int Jb_list, 
      int Jb_list_nbs)
{
  struct stringwr *Ib;
  unsigned int Ib_idx;
  struct pthreads_s1vfci **thread_info;
  int i;

  thread_info = (struct pthreads_s1vfci **)
                malloc(sizeof(struct pthreads_s1vfci *) * nbs);
  for (i=0; i<nbs; i++) {
      thread_info[i] = (struct pthreads_s1vfci *)
                       malloc(sizeof(struct pthreads_s1vfci));
    }

  tpool_queue_open(thread_pool);
  
  detci_time.s1_mt_before_time = wall_time_new();
  /* loop over I_b */
  for (Ib=betlist[Ib_list], Ib_idx=0; Ib_idx < nbs; Ib_idx++, Ib++) {
      thread_info[Ib_idx]->alplist=alplist;
      thread_info[Ib_idx]->betlist=betlist;
      thread_info[Ib_idx]->C=C;
      thread_info[Ib_idx]->S=S;
      thread_info[Ib_idx]->oei=oei;
      thread_info[Ib_idx]->tei=tei;
      thread_info[Ib_idx]->nlists=nlists;
      thread_info[Ib_idx]->nas=nas;
      thread_info[Ib_idx]->nbs=nbs;
      thread_info[Ib_idx]->Ib_list=Ib_list;
      thread_info[Ib_idx]->Jb_list=Jb_list;
      thread_info[Ib_idx]->Jb_list_nbs=Jb_list_nbs;
      thread_info[Ib_idx]->Ib=Ib;
      thread_info[Ib_idx]->Ib_idx=Ib_idx;
      tpool_add_work(thread_pool, s1_block_vras_pthread, (void *) thread_info[Ib_idx]);
    } /* end loop over Ib */
  tpool_queue_close(thread_pool, 1);

  detci_time.s1_mt_after_time = wall_time_new();
  detci_time.s1_mt_total_time += detci_time.s1_mt_after_time - detci_time.s1_mt_before_time;

  for (i=0; i<nbs; i++) free(thread_info[i]);

}

/*
** S1_BLOCK_VRAS_PTHREAD.C: 
** 
** Calculate the sigma_1 vector as described by
** equation (20) of RAS Paper (Olsen, Roos, Jorgensen, Aa. Jensen JCP 1988)
**
** This sigma1 routine is for Full CI's only, assumes (ij|ij)'s have not
** been halved, and attempts to follow Olsen's vectorized algorithm more
** closely than previous versions, using sparsity of F.
** 
** David Sherrill, 10 May 1996
** Based on previous code by David Sherrill, 1994-5
**
** Updated 3/27/94 to include g matrix for RAS
** Modified 4/8/94 to make C and s one-dimensional
** Modified 4/10/94 to make FCI-only (for now) and use new string structs
** Modified 6/21/95 for use in new RAS program (C, s now 2D again!)
** Modified 8/2/95 to make RAS again
** Modified 5/10/96 for new sparse-F method
*/
void s1_block_vras_pthread(void *threadarg)
{
  struct stringwr *Ib, *Kb, **alplist, **betlist;
  unsigned int Ia_idx, Ib_idx, Kb_idx, Jb_idx;
  unsigned int Ibcnt, Kbcnt, Kb_list, Ib_ex, Kb_ex;
  unsigned int *Ibridx, *Kbridx;
  int nirreps,  *Ibij, *Kbij, *Iboij, *Kboij;
  signed char *Ibsgn, *Kbsgn;
  int ij,kl,ijkl,oij,okl;
  double Kb_sgn, Jb_sgn;
  double tval;
  double *oei, *tei, **C, **S, *F;
  struct pthreads_s1vfci *thread_info;
  int nlists, nas, nbs, Ib_list, Jb_list, Jb_list_nbs;
  
  thread_info = (struct pthreads_s1vfci *) threadarg;
  alplist = thread_info->alplist;
  betlist = thread_info->betlist;
  C = thread_info->C;
  S = thread_info->S;
  oei = thread_info->oei;
  tei = thread_info->tei;
  nlists = thread_info->nlists;
  nas = thread_info->nas;
  nbs = thread_info->nbs;
  Ib_list = thread_info->Ib_list;
  Jb_list = thread_info->Jb_list;
  Jb_list_nbs = thread_info->Jb_list_nbs;
  Ib = thread_info->Ib;
  Ib_idx = thread_info->Ib_idx;

  nirreps = CalcInfo.nirreps;
  F = init_array(Jb_list_nbs);
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
              ijkl = INDEX(ij,kl);
              if (oij > okl) 
                  F[Jb_idx] += Kb_sgn * Jb_sgn * tei[ijkl] ;
              else if (oij == okl) 
                  F[Jb_idx] += 0.5 * Kb_sgn * Jb_sgn * tei[ijkl] ;
            }
        } /* end loop over Ib excitations */
    } /* end loop over Kb_list */

      
  /* 
     for (Ia_idx=0; Ia_idx < nas; Ia_idx++) {
     tval = 0.0;
     for (Jb_idx=0; Jb_idx < Jb_list_nbs; Jb_idx++) {
     tval += C[Ia_idx][Jb_idx] * F[Jb_idx];
     }
     S[Ia_idx][Ib_idx] += tval;
     }
  */

  /* need to improve mem access pattern here! Above vers may be better!  */
  /* min op cnt may also be better */
  for (Jb_idx=0; Jb_idx < Jb_list_nbs; Jb_idx++) {
      if ((tval=F[Jb_idx]) == 0.0) continue;

#ifdef USE_BLAS
      C_DAXPY(nas,tval, (C[0]+Jb_idx), Jb_list_nbs, (S[0]+Ib_idx), nbs);
#else
      for (Ia_idx=0; Ia_idx < nas; Ia_idx++) {
          S[Ia_idx][Ib_idx] += tval * C[Ia_idx][Jb_idx];
        }
#endif
    }
  free(F);
}

/*
** S1_BLOCK_VRAS_ROTF
** 
** String replacements on-the-fly version
**
** This sigma1 routine is for RAS CI's.
** currently assumes that (ij|ij)'s have not been halved!! 
** 
** David Sherrill, 13 May 1996
** Based on previous code by David Sherrill, 1994-5
**
** Updated 3/27/94 to include g matrix for RAS
** Modified 4/8/94 to make C and s one-dimensional
** Modified 4/10/94 to make FCI-only (for now) and use new string structs
** Modified 6/21/95 for use in new RAS program (C, s now 2D again!)
** Modified 8/2/95 to make RAS again
** Modified 5/13/96 for new sparse-F method
**
*/
void s1_block_vras_rotf(int *Cnt[2], int **Ij[2], int **Oij[2],
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
               ijkl = INDEX(ij,kl);
               if (oij > okl) 
                  F[Jb_idx] += Kb_sgn * Jb_sgn * tei[ijkl] ;
               else if (oij == okl) 
                  F[Jb_idx] += 0.5 * Kb_sgn * Jb_sgn * tei[ijkl] ;
               }
            } /* end loop over Ib excitations */

      /*
      for (Ia_idx=0; Ia_idx < nas; Ia_idx++) {
         tval = 0.0;
         for (Jb_idx=0; Jb_idx < Jb_list_nbs; Jb_idx++) {
            tval += C[Ia_idx][Jb_idx] * F[Jb_idx];
            }
         S[Ia_idx][Ib_idx] += tval;
         }
      */

      /* need to improve mem access pattern here! Above vers may be better!  */
      /* min op cnt may also be better */
      for (Jb_idx=0; Jb_idx < Jb_list_nbs; Jb_idx++) {
         if ((tval=F[Jb_idx]) == 0.0) continue;
         for (Ia_idx=0; Ia_idx < nas; Ia_idx++) {
            S[Ia_idx][Ib_idx] += tval * C[Ia_idx][Jb_idx];
            }
         }

      } /* end loop over Ib */
   } /* end loop over Kb_list */

}

}} // namespace psi::detci

