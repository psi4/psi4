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
#include <pthread.h>
#include "tpool.h"

namespace psi { namespace detci {

extern unsigned char ***Occs;

extern void b2brepl(unsigned char **occs, int *Jcnt, int **Jij, int **Joij, 
      int **Jridx, signed char **Jsgn, struct olsen_graph *Graph,
      int Ilist, int Jlist, int len);
void s2_block_vfci_pthread(void *threadarg);
void s2_block_vras_pthread(void *threadarg);

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

/*
** S2_BLOCK_VFCI()
** 
** Calculate the sigma_2 vector as described by
** equation (20) of RAS Paper (Olsen, Roos, Jorgensen, Aa. Jensen JCP 1988)
**
** This sigma1 routine is for Full CI's only, assumes (ij|ij)'s have not
** been halved, and attempts to follow Olsen's vectorized algorithm more
** closely than previous versions, using sparsity of F.
** 
** David Sherrill, 18 April 1996
** Based on many previous versions by David Sherrill 1994-5
*/
void s2_block_vfci(struct stringwr **alplist, struct stringwr **betlist, 
      double **C, double **S, double *oei, double *tei, double *F,
      int nlists, int nas, int nbs, int Ia_list, int Ja_list, 
      int Ja_list_nas)
{
   struct stringwr *Ia, *Ka;
   unsigned int Ia_idx, Ib_idx, Ka_idx, Ja_idx;
   unsigned int Iacnt, Kacnt, Ka_list, Ia_ex, Ka_ex;
   unsigned int *Iaridx, *Karidx;
   int *Iaij, *Kaij;
   signed char *Iasgn, *Kasgn;
   int ij,kl,ijkl;
   double Ka_sgn, Ja_sgn;
   double tval;
   double *Sptr, *Cptr;


   /* loop over I_a */
   for (Ia=alplist[Ia_list], Ia_idx=0; Ia_idx < nas; Ia_idx++, Ia++) {

      Sptr = S[Ia_idx];
      zero_arr(F, Ja_list_nas);

      /* loop over excitations E^a_{kl} from |A(I_a)> */
      for (Ka_list=0; Ka_list < nlists; Ka_list++) {
            Iacnt = Ia->cnt[Ka_list];
            Iaridx = Ia->ridx[Ka_list];
            Iasgn = Ia->sgn[Ka_list];
            Iaij = Ia->ij[Ka_list];
            for (Ia_ex=0; Ia_ex < Iacnt; Ia_ex++) {
               kl = *Iaij++;
               Ka_idx = *Iaridx++;
               Ka_sgn = (double) *Iasgn++;

               /* A(K_a) = sgn(kl) * E^a_{kl} |A(I_a)> */
               Ka = alplist[Ka_list] + Ka_idx;
               if (Ka_list == Ja_list) F[Ka_idx] += Ka_sgn * oei[kl];

               /* loop over excitations E^a_{ij} from |A(K_a)> */
               /* Ja_list pre-determined because of C blocking */
               Kacnt = Ka->cnt[Ja_list];
               Karidx = Ka->ridx[Ja_list];
               Kasgn = Ka->sgn[Ja_list];
               Kaij = Ka->ij[Ja_list];
               for (Ka_ex=0; Ka_ex < Kacnt; Ka_ex++) {
                  Ja_idx = *Karidx++;
                  Ja_sgn = (double) *Kasgn++;
                  ij = *Kaij++;
                  ijkl = INDEX(ij,kl);
                  F[Ja_idx] += 0.5 * Ka_sgn * Ja_sgn * tei[ijkl] ;
                  }
               } /* end loop over Ia excitations */
         } /* end loop over Ka_list */

      
      /*
      for (Ib_idx=0; Ib_idx < nbs; Ib_idx++) {
         tval = 0.0;
         for (Ja_idx=0; Ja_idx < Ja_list_nas; Ja_idx++) {
            tval += C[Ja_idx][Ib_idx] * F[Ja_idx];
            }
         S[Ia_idx][Ib_idx] += tval;
         }
      */

      for (Ja_idx=0; Ja_idx < Ja_list_nas; Ja_idx++) {
         if ((tval=F[Ja_idx]) == 0.0) continue;
         Cptr = C[Ja_idx];

      #ifdef USE_BLAS
         C_DAXPY(nbs, tval, Cptr, 1, Sptr, 1);
      #else   
         for (Ib_idx=0; Ib_idx < nbs; Ib_idx++) {
            Sptr[Ib_idx] += tval * Cptr[Ib_idx];
            }
      #endif
         }

      } /* end loop over Ia */

}


/*
** S2_BLOCK_VFCI_THREAD()
** 
** Calculate the sigma_2 vector as described by
** equation (20) of RAS Paper (Olsen, Roos, Jorgensen, Aa. Jensen JCP 1988)
**
** This sigma1 routine is for Full CI's only, assumes (ij|ij)'s have not
** been halved, and attempts to follow Olsen's vectorized algorithm more
** closely than previous versions, using sparsity of F.
** 
** David Sherrill, 18 April 1996
** Based on many previous versions by David Sherrill 1994-5
*/
void s2_block_vfci_thread(struct stringwr **alplist, struct stringwr **betlist, 
      double **C, double **S, double *oei, double *tei, double *F,
      int nlists, int nas, int nbs, int Ia_list, int Ja_list, 
      int Ja_list_nas)
{
  struct stringwr *Ia, *Ka;
  unsigned int Ia_idx, Ib_idx, Ka_idx, Ja_idx;
  unsigned int Iacnt, Kacnt, Ka_list, Ia_ex, Ka_ex;
  unsigned int *Iaridx, *Karidx;
  int *Iaij, *Kaij;
  signed char *Iasgn, *Kasgn;
  int ij,kl,ijkl;
  double Ka_sgn, Ja_sgn;
  double tval;
  double *Sptr, *Cptr;
  struct pthreads_s2vfci **thread_info;
  int i;
   
  thread_info = (struct pthreads_s2vfci **)
                malloc(sizeof(struct pthreads_s2vfci *) * nas);
  for (i=0; i<nas; i++) {
      thread_info[i] = (struct pthreads_s2vfci *)
                       malloc(sizeof(struct pthreads_s2vfci));
    }

  tpool_queue_open(thread_pool);

  detci_time.s2_mt_before_time = wall_time_new();
  /* loop over I_a */
  for (Ia=alplist[Ia_list], Ia_idx=0; Ia_idx < nas; Ia_idx++, Ia++) {
      thread_info[Ia_idx]->alplist=alplist;
      thread_info[Ia_idx]->betlist=betlist;
      thread_info[Ia_idx]->C=C;
      thread_info[Ia_idx]->S=S;
      thread_info[Ia_idx]->oei=oei;
      thread_info[Ia_idx]->tei=tei;
      thread_info[Ia_idx]->nlists=nlists;
      thread_info[Ia_idx]->nas=nas;
      thread_info[Ia_idx]->nbs=nbs;
      thread_info[Ia_idx]->Ia_list=Ia_list;
      thread_info[Ia_idx]->Ja_list=Ja_list;
      thread_info[Ia_idx]->Ja_list_nas=Ja_list_nas;
      thread_info[Ia_idx]->Ia=Ia;
      thread_info[Ia_idx]->Ia_idx=Ia_idx;
      tpool_add_work(thread_pool, s2_block_vfci_pthread, (void *) thread_info[Ia_idx]);
    } /* end loop over Ia */
  tpool_queue_close(thread_pool, 1);

  detci_time.s2_mt_after_time = wall_time_new();
  detci_time.s2_mt_total_time += detci_time.s2_mt_after_time - detci_time.s2_mt_before_time;

  for (i=0; i<nas; i++) free(thread_info[i]);
}


/*
** S2_BLOCK_VFCI_PTHREAD()
** 
** Calculate the sigma_2 vector as described by
** equation (20) of RAS Paper (Olsen, Roos, Jorgensen, Aa. Jensen JCP 1988)
**
** This sigma1 routine is for Full CI's only, assumes (ij|ij)'s have not
** been halved, and attempts to follow Olsen's vectorized algorithm more
** closely than previous versions, using sparsity of F.
** 
** David Sherrill, 18 April 1996
** Based on many previous versions by David Sherrill 1994-5
*/
void s2_block_vfci_pthread(void *threadarg)
{
  struct stringwr *Ia, *Ka, **alplist, **betlist;
  unsigned int Ia_idx, Ib_idx, Ka_idx, Ja_idx;
  unsigned int Iacnt, Kacnt, Ka_list, Ia_ex, Ka_ex;
  unsigned int *Iaridx, *Karidx;
  int *Iaij, *Kaij;
  signed char *Iasgn, *Kasgn;
  int ij,kl,ijkl;
  double Ka_sgn, Ja_sgn;
  double tval;
  double *Sptr, *Cptr, *oei, *tei;
  double *F, **C, **S;
  struct pthreads_s2vfci *thread_info;
  int nlists, nas, nbs, Ia_list, Ja_list, Ja_list_nas;
  
  thread_info = (struct pthreads_s2vfci *) threadarg;
  alplist = thread_info->alplist;
  betlist = thread_info->betlist;
  C = thread_info->C;
  S = thread_info->S;
  oei = thread_info->oei;
  tei = thread_info->tei;
  nlists = thread_info->nlists;
  nas = thread_info->nas;
  nbs = thread_info->nbs;
  Ia_list = thread_info->Ia_list;
  Ja_list = thread_info->Ja_list;
  Ja_list_nas = thread_info->Ja_list_nas;
  Ia = thread_info->Ia;
  Ia_idx = thread_info->Ia_idx;

  F = init_array(Ja_list_nas);
  Sptr = S[Ia_idx];
  zero_arr(F, Ja_list_nas);

  /* loop over excitations E^a_{kl} from |A(I_a)> */
  for (Ka_list=0; Ka_list < nlists; Ka_list++) {
      Iacnt = Ia->cnt[Ka_list];
      Iaridx = Ia->ridx[Ka_list];
      Iasgn = Ia->sgn[Ka_list];
      Iaij = Ia->ij[Ka_list];
      for (Ia_ex=0; Ia_ex < Iacnt; Ia_ex++) {
          kl = *Iaij++;
          Ka_idx = *Iaridx++;
          Ka_sgn = (double) *Iasgn++;

          /* A(K_a) = sgn(kl) * E^a_{kl} |A(I_a)> */
          Ka = alplist[Ka_list] + Ka_idx;
          if (Ka_list == Ja_list) F[Ka_idx] += Ka_sgn * oei[kl];

          /* loop over excitations E^a_{ij} from |A(K_a)> */
          /* Ja_list pre-determined because of C blocking */
          Kacnt = Ka->cnt[Ja_list];
          Karidx = Ka->ridx[Ja_list];
          Kasgn = Ka->sgn[Ja_list];
          Kaij = Ka->ij[Ja_list];
          for (Ka_ex=0; Ka_ex < Kacnt; Ka_ex++) {
              Ja_idx = *Karidx++;
              Ja_sgn = (double) *Kasgn++;
              ij = *Kaij++;
              ijkl = INDEX(ij,kl);
              F[Ja_idx] += 0.5 * Ka_sgn * Ja_sgn * tei[ijkl] ;
            }
        } /* end loop over Ia excitations */
    } /* end loop over Ka_list */

      
  /*
    for (Ib_idx=0; Ib_idx < nbs; Ib_idx++) {
    tval = 0.0;
    for (Ja_idx=0; Ja_idx < Ja_list_nas; Ja_idx++) {
    tval += C[Ja_idx][Ib_idx] * F[Ja_idx];
    }
    S[Ia_idx][Ib_idx] += tval;
    }
  */

  for (Ja_idx=0; Ja_idx < Ja_list_nas; Ja_idx++) {
      if ((tval=F[Ja_idx]) == 0.0) continue;
      Cptr = C[Ja_idx];

#ifdef USE_BLAS
      C_DAXPY(nbs, tval, Cptr, 1, Sptr, 1);
#else   
      for (Ib_idx=0; Ib_idx < nbs; Ib_idx++) {
          Sptr[Ib_idx] += tval * Cptr[Ib_idx];
        }
#endif
    }
  free(F);
}

/*
** S2_BLOCK_VRAS()
** 
** Calculate the sigma_2 vector as described by
** equation (20) of RAS Paper (Olsen, Roos, Jorgensen, Aa. Jensen JCP 1988)
**
** This sigma2 routine is for RAS CI's.
** currently assumes that (ij|ij)'s have not been halved!! 
** 
** David Sherrill, 10 May 1996
** Based on previous code by David Sherrill, 1994-5
**
** Updated 3/27/94 to include g matrix for RAS
** Modified 4/8/94 to make C and s one-dimensional
** Modified 4/10/94 to make FCI-only (for now) and use new string structs
** Modified 6/21/95 for use in new RAS program
** Obtained 7/22/95 from s1 routine by changing a's to b's and vice versa
** Modified 5/10/96 for more vectorized approach
*/
void s2_block_vras(struct stringwr **alplist, struct stringwr **betlist, 
      double **C, double **S, double *oei, double *tei, double *F,
      int nlists, int nas, int nbs, int Ia_list, int Ja_list, 
      int Ja_list_nas)
{
   struct stringwr *Ia, *Ka;
   unsigned int Ia_idx, Ib_idx, Ka_idx, Ja_idx;
   unsigned int Iacnt, Kacnt, Ka_list, Ia_ex, Ka_ex;
   unsigned int *Iaridx, *Karidx;
   int nirreps, *Iaij, *Kaij, *Iaoij, *Kaoij;
   signed char *Iasgn, *Kasgn;
   int ij,kl,ijkl,oij,okl;
   double Ka_sgn, Ja_sgn;
   double tval;
   double *Sptr, *Cptr;

   nirreps = CalcInfo.nirreps;

   /* loop over I_a */
   for (Ia=alplist[Ia_list], Ia_idx=0; Ia_idx < nas; Ia_idx++, Ia++) {

      Sptr = S[Ia_idx];
      zero_arr(F, Ja_list_nas);

      /* loop over excitations E^a_{kl} from |A(I_a)> */
      for (Ka_list=0; Ka_list < nlists; Ka_list++) {
            Iacnt = Ia->cnt[Ka_list];
            Iaridx = Ia->ridx[Ka_list];
            Iasgn = Ia->sgn[Ka_list];
            Iaij = Ia->ij[Ka_list];
            Iaoij = Ia->oij[Ka_list];
            for (Ia_ex=0; Ia_ex < Iacnt; Ia_ex++) {
               kl = *Iaij++;
               okl = *Iaoij++;
               Ka_idx = *Iaridx++;
               Ka_sgn = (double) *Iasgn++;

               /* A(K_a) = sgn(kl) * E^a_{kl} |A(I_a)> */
               Ka = alplist[Ka_list] + Ka_idx;
               /* note okl on next line, not kl */
               if (Ka_list == Ja_list) F[Ka_idx] += Ka_sgn * oei[okl];

               /* loop over excitations E^a_{ij} from |A(K_a)> */
               /* Ja_list pre-determined because of C blocking */
               Kacnt = Ka->cnt[Ja_list];
               Karidx = Ka->ridx[Ja_list];
               Kasgn = Ka->sgn[Ja_list];
               Kaij = Ka->ij[Ja_list];
               Kaoij = Ka->oij[Ja_list];
               for (Ka_ex=0; Ka_ex < Kacnt; Ka_ex++) {
                  Ja_idx = *Karidx++;
                  Ja_sgn = (double) *Kasgn++;
                  ij = *Kaij++;
                  oij = *Kaoij++;
                  ijkl = INDEX(ij,kl);
                  if (oij > okl) 
                     F[Ja_idx] += Ka_sgn * Ja_sgn * tei[ijkl] ;
                  else if (oij == okl) 
                     F[Ja_idx] += 0.5 * Ka_sgn * Ja_sgn * tei[ijkl] ;
                  }
               } /* end loop over Ia excitations */
         } /* end loop over Ka_list */
      
      /*
      for (Ib_idx=0; Ib_idx < nbs; Ib_idx++) {
         tval = 0.0;
         for (Ja_idx=0; Ja_idx < Ja_list_nas; Ja_idx++) {
            tval += C[Ja_idx][Ib_idx] * F[Ja_idx];
            }
         S[Ia_idx][Ib_idx] += tval;
         }
      */

      for (Ja_idx=0; Ja_idx < Ja_list_nas; Ja_idx++) {
         if ((tval=F[Ja_idx]) == 0.0) continue;
         Cptr = C[Ja_idx];
      #ifdef USE_BLAS
         C_DAXPY(nbs, tval, Cptr, 1, Sptr, 1);
      #else
         for (Ib_idx=0; Ib_idx < nbs; Ib_idx++) {
            Sptr[Ib_idx] += tval * Cptr[Ib_idx];
            }
      #endif
         }

      } /* end loop over Ia */

}


/*
** S2_BLOCK_VRAS_THREAD()
** 
** Calculate the sigma_2 vector as described by
** equation (20) of RAS Paper (Olsen, Roos, Jorgensen, Aa. Jensen JCP 1988)
**
** This sigma2 routine is for RAS CI's.
** currently assumes that (ij|ij)'s have not been halved!! 
** 
** David Sherrill, 10 May 1996
** Based on previous code by David Sherrill, 1994-5
**
** Updated 3/27/94 to include g matrix for RAS
** Modified 4/8/94 to make C and s one-dimensional
** Modified 4/10/94 to make FCI-only (for now) and use new string structs
** Modified 6/21/95 for use in new RAS program
** Obtained 7/22/95 from s1 routine by changing a's to b's and vice versa
** Modified 5/10/96 for more vectorized approach
*/
void s2_block_vras_thread(struct stringwr **alplist, struct stringwr **betlist, 
      double **C, double **S, double *oei, double *tei, double *F,
      int nlists, int nas, int nbs, int Ia_list, int Ja_list, 
      int Ja_list_nas)
{
  struct stringwr *Ia;
  unsigned int Ia_idx;
  struct pthreads_s2vfci **thread_info;
  int i;

  thread_info = (struct pthreads_s2vfci **)
                malloc(sizeof(struct pthreads_s2vfci *) * nas);
  for (i=0; i<nas; i++) {
      thread_info[i] = (struct pthreads_s2vfci *)
                       malloc(sizeof(struct pthreads_s2vfci));
    }

  tpool_queue_open(thread_pool);
  detci_time.s2_mt_before_time = wall_time_new();

 
  /* loop over I_a */
  for (Ia=alplist[Ia_list], Ia_idx=0; Ia_idx < nas; Ia_idx++, Ia++) {
      thread_info[Ia_idx]->alplist=alplist;
      thread_info[Ia_idx]->betlist=betlist;
      thread_info[Ia_idx]->C=C;
      thread_info[Ia_idx]->S=S;
      thread_info[Ia_idx]->oei=oei;
      thread_info[Ia_idx]->tei=tei;
      thread_info[Ia_idx]->nlists=nlists;
      thread_info[Ia_idx]->nas=nas;
      thread_info[Ia_idx]->nbs=nbs;
      thread_info[Ia_idx]->Ia_list=Ia_list;
      thread_info[Ia_idx]->Ja_list=Ja_list;
      thread_info[Ia_idx]->Ja_list_nas=Ja_list_nas;
      thread_info[Ia_idx]->Ia=Ia;
      thread_info[Ia_idx]->Ia_idx=Ia_idx;
      tpool_add_work(thread_pool, s2_block_vras_pthread, (void *) thread_info[Ia_idx]);
    } /* end loop over Ia */

  tpool_queue_close(thread_pool, 1);
  detci_time.s2_mt_after_time = wall_time_new();
  detci_time.s2_mt_total_time += detci_time.s2_mt_after_time - detci_time.s2_mt_before_time;


  for (i=0; i<nas; i++) free(thread_info[i]);
   
}


/*
** S2_BLOCK_VRAS_PTHREAD()
** 
** Calculate the sigma_2 vector as described by
** equation (20) of RAS Paper (Olsen, Roos, Jorgensen, Aa. Jensen JCP 1988)
**
** This sigma2 routine is for RAS CI's.
** currently assumes that (ij|ij)'s have not been halved!! 
** 
** David Sherrill, 10 May 1996
** Based on previous code by David Sherrill, 1994-5
**
** Updated 3/27/94 to include g matrix for RAS
** Modified 4/8/94 to make C and s one-dimensional
** Modified 4/10/94 to make FCI-only (for now) and use new string structs
** Modified 6/21/95 for use in new RAS program
** Obtained 7/22/95 from s1 routine by changing a's to b's and vice versa
** Modified 5/10/96 for more vectorized approach
*/
void s2_block_vras_pthread(void *threadarg)
{
  struct stringwr *Ia, *Ka, **alplist, **betlist;
  unsigned int Ia_idx, Ib_idx, Ka_idx, Ja_idx;
  unsigned int Iacnt, Kacnt, Ka_list, Ia_ex, Ka_ex;
  unsigned int *Iaridx, *Karidx;
  int nirreps, *Iaij, *Kaij, *Iaoij, *Kaoij;
  signed char *Iasgn, *Kasgn;
  int ij,kl,ijkl,oij,okl;
  double Ka_sgn, Ja_sgn;
  double tval;
  double *Sptr, *Cptr, *oei, *tei, **C, **S, *F;
  struct pthreads_s2vfci *thread_info;
  int nlists, nas, nbs, Ia_list, Ja_list, Ja_list_nas;

  thread_info = (struct pthreads_s2vfci *) threadarg;
  alplist = thread_info->alplist;
  betlist = thread_info->betlist;
  C = thread_info->C;
  S = thread_info->S;
  oei = thread_info->oei;
  tei = thread_info->tei;
  nlists = thread_info->nlists;
  nas = thread_info->nas;
  nbs = thread_info->nbs;
  Ia_list = thread_info->Ia_list;
  Ja_list = thread_info->Ja_list;
  Ja_list_nas = thread_info->Ja_list_nas;
  Ia = thread_info->Ia;
  Ia_idx = thread_info->Ia_idx;

  F = init_array(Ja_list_nas);
  Sptr = S[Ia_idx];
  zero_arr(F, Ja_list_nas);

  nirreps = CalcInfo.nirreps;

  /* loop over excitations E^a_{kl} from |A(I_a)> */
  for (Ka_list=0; Ka_list < nlists; Ka_list++) {
      Iacnt = Ia->cnt[Ka_list];
      Iaridx = Ia->ridx[Ka_list];
      Iasgn = Ia->sgn[Ka_list];
      Iaij = Ia->ij[Ka_list];
      Iaoij = Ia->oij[Ka_list];
      for (Ia_ex=0; Ia_ex < Iacnt; Ia_ex++) {
          kl = *Iaij++;
          okl = *Iaoij++;
          Ka_idx = *Iaridx++;
          Ka_sgn = (double) *Iasgn++;

          /* A(K_a) = sgn(kl) * E^a_{kl} |A(I_a)> */
          Ka = alplist[Ka_list] + Ka_idx;
          /* note okl on next line, not kl */
          if (Ka_list == Ja_list) F[Ka_idx] += Ka_sgn * oei[okl];

          /* loop over excitations E^a_{ij} from |A(K_a)> */
          /* Ja_list pre-determined because of C blocking */
          Kacnt = Ka->cnt[Ja_list];
          Karidx = Ka->ridx[Ja_list];
          Kasgn = Ka->sgn[Ja_list];
          Kaij = Ka->ij[Ja_list];
          Kaoij = Ka->oij[Ja_list];
          for (Ka_ex=0; Ka_ex < Kacnt; Ka_ex++) {
              Ja_idx = *Karidx++;
              Ja_sgn = (double) *Kasgn++;
              ij = *Kaij++;
              oij = *Kaoij++;
              ijkl = INDEX(ij,kl);
              if (oij > okl) 
                  F[Ja_idx] += Ka_sgn * Ja_sgn * tei[ijkl] ;
              else if (oij == okl) 
                  F[Ja_idx] += 0.5 * Ka_sgn * Ja_sgn * tei[ijkl] ;
            }
        } /* end loop over Ia excitations */
    } /* end loop over Ka_list */
      
  /*
    for (Ib_idx=0; Ib_idx < nbs; Ib_idx++) {
    tval = 0.0;
    for (Ja_idx=0; Ja_idx < Ja_list_nas; Ja_idx++) {
    tval += C[Ja_idx][Ib_idx] * F[Ja_idx];
    }
    S[Ia_idx][Ib_idx] += tval;
    }
  */

  for (Ja_idx=0; Ja_idx < Ja_list_nas; Ja_idx++) {
      if ((tval=F[Ja_idx]) == 0.0) continue;
      Cptr = C[Ja_idx];
#ifdef USE_BLAS
      C_DAXPY(nbs, tval, Cptr, 1, Sptr, 1);
#else
      for (Ib_idx=0; Ib_idx < nbs; Ib_idx++) {
          Sptr[Ib_idx] += tval * Cptr[Ib_idx];
        }
#endif
    }
  free(F);
}


/*
** S2_BLOCK_VRAS_ROTF()
** 
** s2_block_vras_rotf(): Calculate the sigma_2 vector as described by
** equation (20) of RAS Paper (Olsen, Roos, Jorgensen, Aa. Jensen JCP 1988)
**
** String replacements on-the-fly version
** currently assumes that (ij|ij)'s have not been halved!! 
** 
** This sigma2 routine is for RAS CI's.
** 
** David Sherrill, 13 May 1996
** Based on previous code by David Sherrill, 1994-5
**
** Updated 3/27/94 to include g matrix for RAS
** Modified 4/8/94 to make C and s one-dimensional
** Modified 4/10/94 to make FCI-only (for now) and use new string structs
** Modified 6/21/95 for use in new RAS program
** Obtained 7/22/95 from s1 routine by changing a's to b's and vice versa
** Modified 5/13/96 for new sparse-F vectorized version
**
*/
void s2_block_vras_rotf(int *Cnt[2], int **Ij[2], int **Oij[2],
      int **Ridx[2], signed char **Sgn[2], unsigned char **Toccs,
      double **C, double **S,
      double *oei, double *tei, double *F, int nlists, int nas, int nbs,
      int Ia_list, int Ja_list, int Ja_list_nas)
{
   int Ia_idx, Ib_idx, Ka_idx, Ja_idx;
   int Iacnt, Kacnt, Ka_list, Ia_ex, Ka_ex;
   int *Iaridx, *Karidx;
   int nirreps, *Iaij, *Kaij, *Iaoij, *Kaoij;
   signed char *Iasgn, *Kasgn;
   int i,ij,kl,ijkl,oij,okl;
   double Ka_sgn, Ja_sgn;
   double tval, *Cptr, *Sptr;

   nirreps = CalcInfo.nirreps;

   for (Ka_list=0; Ka_list < nlists; Ka_list++) {
      b2brepl(Occs[Ia_list], Cnt[0], Ij[0], Oij[0], Ridx[0],
         Sgn[0], BetaG, Ia_list, Ka_list, nas);

      /* loop over I_a */
      for (Ia_idx=0; Ia_idx < nas; Ia_idx++) {

         if ((Iacnt = Cnt[0][Ia_idx]) < 0) continue;
         Sptr = S[Ia_idx];
         zero_arr(F, Ja_list_nas);

         /* loop over excitations E^a_{kl} from |A(I_a)> */
         Iaridx = Ridx[0][Ia_idx];
         Iasgn = Sgn[0][Ia_idx];
         Iaij = Ij[0][Ia_idx];
         Iaoij = Oij[0][Ia_idx];

         for (i=0; i<Iacnt; i++) 
            Toccs[i] = Occs[Ka_list][Iaridx[i]];

         b2brepl(Toccs, Cnt[1], Ij[1], Oij[1], Ridx[1], Sgn[1],
            AlphaG, Ka_list, Ja_list, Iacnt);

         for (Ia_ex=0; Ia_ex < Iacnt; Ia_ex++) {
            kl = *Iaij++;
            okl = *Iaoij++;
            Ka_idx = *Iaridx++;
            Ka_sgn = (double) *Iasgn++;

            /* A(K_a) = sgn(kl) * E^a_{kl} |A(I_a)> */
            /* note okl on next line, not kl */
            if (Ka_list == Ja_list) F[Ka_idx] += Ka_sgn * oei[okl];

            /* loop over excitations E^a_{ij} from |A(K_a)> */
            /* Ja_list pre-determined because of C blocking */
            Kacnt = Cnt[1][Ia_ex];
            Karidx = Ridx[1][Ia_ex];
            Kasgn = Sgn[1][Ia_ex];
            Kaij = Ij[1][Ia_ex];
            Kaoij = Oij[1][Ia_ex];
            for (Ka_ex=0; Ka_ex < Kacnt; Ka_ex++) {
               Ja_idx = *Karidx++;
               Ja_sgn = (double) *Kasgn++;
               ij = *Kaij++;
               oij = *Kaoij++;
               ijkl = INDEX(ij,kl);
               if (oij > okl) 
                  F[Ja_idx] += Ka_sgn * Ja_sgn * tei[ijkl] ;
               else if (oij == okl) 
                  F[Ja_idx] += 0.5 * Ka_sgn * Ja_sgn * tei[ijkl] ;
               }
            } /* end loop over Ia excitations */

      /*
      for (Ib_idx=0; Ib_idx < nbs; Ib_idx++) {
         tval = 0.0;
         for (Ja_idx=0; Ja_idx < Ja_list_nas; Ja_idx++) {
            tval += C[Ja_idx][Ib_idx] * F[Ja_idx];
            }
         S[Ia_idx][Ib_idx] += tval;
         }
      */

      for (Ja_idx=0; Ja_idx < Ja_list_nas; Ja_idx++) {
         if ((tval=F[Ja_idx]) == 0.0) continue;
         Cptr = C[Ja_idx];
         for (Ib_idx=0; Ib_idx < nbs; Ib_idx++) {
            Sptr[Ib_idx] += tval * Cptr[Ib_idx];
            }
         }

      } /* end loop over Ia */
   } /* end loop over Ka_list */

}

}} // namespace psi::detci

