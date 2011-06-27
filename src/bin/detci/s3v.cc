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

/* Global variables necessary for pthreads */
pthread_mutex_t inc_Ia_mutex;
struct stringwr *Ia_global;
int Ia_idx_global;


int form_ilist(struct stringwr *alplist, int Ja_list, int nas, int kl,
   int *L, int *R, double *Sgn);
int form_ilist_rotf(int *Cnt, int **Ridx, signed char **Sn, int **Ij,
   int nas, int kl, int *L, int *R, double *Sgn);
void s3_block_vdiag_pthread(void *threadarg);
void s3_block_v_pthread(void *threadarg);


#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))


/*
** S3_BLOCK_VDIAG()
**
** Calculate a block of the sigma3 vector in equation (9c) of
** Olsen, Roos, et al.  For diagonal blocks of sigma.
**
** currently assumes that (ij|ij)'s have not been halved
** Try to get the Olsen vector version working....again!!!!
*/
void s3_block_vdiag(struct stringwr *alplist, struct stringwr *betlist,
      double **C, double **S, double *tei, int nas, int nbs, int cnas,
      int Ib_list, int Ja_list, int Jb_list, int Ib_sym, int Jb_sym,
      double **Cprime, double *F, double *V, double *Sgn, int *L, int *R)
{
  struct stringwr *Ia;
  unsigned int Ia_ex;
  int ij, i, j, t, kl, I, J, RJ;
  double tval, VS, *CprimeI0, *CI0;
  int jlen, Jacnt, *Iaij, *orbsym, norbs, Ia_idx;
  unsigned int *Iaridx;
  signed char *Iasgn;
  double *Tptr;
  struct pthreads_s3diag **thread_info;
  int npthreads, rc, status;
  pthread_t *thread;
   
  npthreads = Parameters.nthreads-1;  /* subtract out the main thread */

  thread = (pthread_t *) malloc(sizeof(pthread_t)*Parameters.nthreads);

  thread_info = (struct pthreads_s3diag **)
                malloc(sizeof(struct pthreads_s3diag *) * nas);

  for (i=0; i<nas; i++) {
      thread_info[i] = (struct pthreads_s3diag *)
                       malloc(sizeof(struct pthreads_s3diag));
    }
  
  norbs = CalcInfo.num_ci_orbs;
  orbsym = CalcInfo.orbsym + CalcInfo.num_fzc_orbs;

  /* loop over i, j */
  for (i=0; i<norbs; i++) {
      for (j=0; j<=i; j++) {
          if ((orbsym[i] ^ orbsym[j] ^ Jb_sym ^ Ib_sym) != 0) continue;
          ij = ioff[i] + j;
          jlen = form_ilist(betlist, Jb_list, nbs, ij, L, R, Sgn);
           
          if (!jlen) continue;
          /*  fprintf(outfile,"S3_BLOCK_VDIAG: ij = %d\t jlen = %d\n", ij, jlen);
           */
          Tptr = tei + ioff[ij];
       
          /* gather operation */
          for (I=0; I<cnas; I++) {
              CprimeI0 = Cprime[I];
              CI0 = C[I];
              for (J=0; J<jlen; J++) {
                  tval = Sgn[J];
                  CprimeI0[J] = CI0[L[J]] * tval;
                }
            }


          /* loop over Ia */
          if (Parameters.nthreads > 1) {
              detci_time.s3_mt_before_time = wall_time_new();
              tpool_queue_open(thread_pool);
              for (Ia=alplist, Ia_idx=0; Ia_idx<nas; Ia_idx++, Ia++) {
                  thread_info[Ia_idx]->nas = nas;
                  thread_info[Ia_idx]->jlen = jlen;
                  thread_info[Ia_idx]->ij = ij;
                  thread_info[Ia_idx]->Cprime = Cprime;    
                  thread_info[Ia_idx]->Ja_list = Ja_list;
                  thread_info[Ia_idx]->Tptr = Tptr;
                  thread_info[Ia_idx]->S = S;
                  thread_info[Ia_idx]->R = R;
                  thread_info[Ia_idx]->Ia_local = Ia;
                  thread_info[Ia_idx]->Ia_idx_local = Ia_idx;
                  tpool_add_work(thread_pool , s3_block_vdiag_pthread, (void *) thread_info[Ia_idx]); 
                }
              tpool_queue_close(thread_pool, 1);
              detci_time.s3_mt_after_time = wall_time_new();
              detci_time.s3_mt_total_time += detci_time.s3_mt_after_time - detci_time.s3_mt_before_time;

            }
          else {
              for (Ia=alplist, Ia_idx=0; Ia_idx<nas; Ia_idx++, Ia++) {

                  /* loop over excitations E^a_{kl} from |A(I_a)> */
                  Jacnt = Ia->cnt[Ja_list];
                  Iaridx = Ia->ridx[Ja_list];
                  Iasgn = Ia->sgn[Ja_list];
                  Iaij = Ia->ij[Ja_list];

                  zero_arr(V, jlen);
/*                  fprintf(outfile,"Ia = %x\t Ia_idx = %d\n", Ia, Ia_idx);
                    fflush(outfile);

                    if (Jacnt) {
                    fprintf(outfile,"S3_BLOCK_VDIAG: Jacnt = %d\n", Jacnt);
                    fflush(outfile);
                    }
*/
                  for (Ia_ex=0; Ia_ex < Jacnt && (kl = *Iaij++)<=ij; Ia_ex++) {
                      I = *Iaridx++;
                      tval = *Iasgn++;
                      if (ij == kl) tval *= 0.5;
                      VS = Tptr[kl] * tval;
                      CprimeI0 = Cprime[I];
           
#ifdef USE_BLAS
                      C_DAXPY(jlen, VS, CprimeI0, 1, V, 1);
#else
                      for (J=0; J<jlen; J++) {
                          V[J] += VS * CprimeI0[J];
                        }
#endif
           
                    }

         
                  /* scatter */
                  for (J=0; J<jlen; J++) {
                      RJ = R[J];
                      S[Ia_idx][RJ] += V[J];
                    }

                } /* end loop over Ia */
            }
       
        } /* end loop over j */
    } /* end loop over i */

  for (i=0; i<nas; i++) free(thread_info[i]);
  free(thread);
  
}              


/*
** S3_BLOCK_VDIAG_PTHREAD()
**
** Multithreaded component of the S3_BLOCK_VDIAG routine
**
** Calculate a block of the sigma3 vector in equation (9c) of
** Olsen, Roos, et al.  For diagonal blocks of sigma.
**
** currently assumes that (ij|ij)'s have not been halved
** Try to get the Olsen vector version working....again!!!!
*/
void s3_block_vdiag_pthread(void *threadarg)
{

  struct stringwr *Ia_local;
  unsigned int Ia_ex;
  int I, J, RJ, kl;
  double tval, VS, *CprimeI0, *CI0;
  int Jacnt, *Iaij;
  unsigned int *Iaridx;
  signed char *Iasgn;
  struct pthreads_s3diag *thread_info_i;

  int nas, jlen, ij, Ja_list;
  double **S, **Cprime, *Tptr, *V;
  int *R, Ia_idx_local;
  int loop_var = 0;
  int thread_id;
  
  thread_info_i = (struct pthreads_s3diag *) threadarg;
  nas = thread_info_i->nas;
  jlen = thread_info_i->jlen;
  ij = thread_info_i->ij;
  Ja_list = thread_info_i->Ja_list;
  S = thread_info_i->S;
  Cprime = thread_info_i->Cprime;
  Tptr = thread_info_i->Tptr;
  R = thread_info_i->R;
  Ia_idx_local = thread_info_i->Ia_idx_local;
  Ia_local = thread_info_i->Ia_local;
  
  V = init_array(jlen);
  /* loop over excitations E^a_{kl} from |A(I_a)> */
  Jacnt = Ia_local->cnt[Ja_list];
  Iaridx = Ia_local->ridx[Ja_list];
  Iasgn = Ia_local->sgn[Ja_list];
  Iaij = Ia_local->ij[Ja_list];

/*
  fprintf(outfile,"Ia_local = %x\t Ia_idx_local = %d\n", Ia_local, Ia_idx_local);
  fflush(outfile);
  
  fprintf(outfile,"Jacnt = %d\n", Jacnt);
*/
      
  zero_arr(V, jlen);
         
  for (Ia_ex=0; Ia_ex < Jacnt && (kl = *Iaij++)<=ij; Ia_ex++) {
      I = *Iaridx++;
      tval = *Iasgn++;
      if (ij == kl) tval *= 0.5;
      VS = Tptr[kl] * tval;
      CprimeI0 = Cprime[I];
           
#ifdef USE_BLAS
      C_DAXPY(jlen, VS, CprimeI0, 1, V, 1);
#else
      for (J=0; J<jlen; J++) {
          V[J] += VS * CprimeI0[J];
        }
#endif
    }

  /* scatter */
  for (J=0; J<jlen; J++) {
      RJ = R[J];
      S[Ia_idx_local][RJ] += V[J];
    }

  free(V);
  /* return 0; */
  
}              


/*
** S3_BLOCK_V()
**
** Calculate a block of the sigma3 vector in equation (9c) of
** Olsen, Roos, et al.  For non-diagonal blocks of s3
**
*/
void s3_block_v(struct stringwr *alplist, struct stringwr *betlist,
      double **C, double **S, double *tei, int nas, int nbs, int cnas,
      int Ib_list, int Ja_list, int Jb_list, int Ib_sym, int Jb_sym,
      double **Cprime, double *F, double *V, double *Sgn, int *L, int *R)
{
   struct stringwr *Ia;
   unsigned int Ia_ex;
   int ij, i, j, kl, ijkl, I, J, RJ;
   double tval, VS, *CprimeI0, *CI0;
   int jlen, Ia_idx, Jacnt, *Iaij, *orbsym, norbs;
   unsigned int *Iaridx;
   signed char *Iasgn;
   double *Tptr;
   struct pthreads_s3diag **thread_info;
   int t, npthreads, rc, status;
   pthread_t *thread;
   
   norbs = CalcInfo.num_ci_orbs;
   orbsym = CalcInfo.orbsym + CalcInfo.num_fzc_orbs;

   thread = (pthread_t *) malloc(sizeof(pthread_t)*Parameters.nthreads);
   thread_info = (struct pthreads_s3diag **)
                  malloc(sizeof(struct pthreads_s3diag *) * nas);
   for (i=0; i<nas; i++) {
       thread_info[i] = (struct pthreads_s3diag *)
                        malloc(sizeof(struct pthreads_s3diag));
     }
   
   /* loop over i, j */
   for (i=0; i<norbs; i++) {
     for (j=0; j<=i; j++) {
       if ((orbsym[i] ^ orbsym[j] ^ Jb_sym ^ Ib_sym) != 0) continue;
       ij = ioff[i] + j;
       jlen = form_ilist(betlist, Jb_list, nbs, ij, L, R, Sgn);

       if (!jlen) continue;

       Tptr = tei + ioff[ij];

       /* gather operation */
       for (I=0; I<cnas; I++) {
         CprimeI0 = Cprime[I];
         CI0 = C[I];
         for (J=0; J<jlen; J++) {
           tval = Sgn[J];
           CprimeI0[J] = CI0[L[J]] * tval;
         }
       }


       /* loop over Ia */
       if (Parameters.nthreads > 1) {
           detci_time.s3_mt_before_time = wall_time_new();
           tpool_queue_open(thread_pool);
           for (Ia=alplist, Ia_idx=0; Ia_idx<nas; Ia_idx++, Ia++) {
               thread_info[Ia_idx]->nas = nas;
               thread_info[Ia_idx]->jlen = jlen;
               thread_info[Ia_idx]->ij = ij;
               thread_info[Ia_idx]->Cprime = Cprime;    
               thread_info[Ia_idx]->Ja_list = Ja_list;
               thread_info[Ia_idx]->Tptr = tei;
               thread_info[Ia_idx]->S = S;
               thread_info[Ia_idx]->R = R;
               thread_info[Ia_idx]->Ia_local = Ia;
               thread_info[Ia_idx]->Ia_idx_local = Ia_idx;
               tpool_add_work(thread_pool, s3_block_v_pthread, (void *) thread_info[Ia_idx]);
             }
           tpool_queue_close(thread_pool, 1);
           detci_time.s3_mt_after_time = wall_time_new();
           detci_time.s3_mt_total_time += detci_time.s3_mt_after_time - detci_time.s3_mt_before_time;

         }
       else {
           for (Ia=alplist, Ia_idx=0; Ia_idx<nas; Ia_idx++, Ia++) {

               /* loop over excitations E^a_{kl} from |A(I_a)> */
               Jacnt = Ia->cnt[Ja_list];
               Iaridx = Ia->ridx[Ja_list];
               Iasgn = Ia->sgn[Ja_list];
               Iaij = Ia->ij[Ja_list];

               zero_arr(V, jlen);

               for (Ia_ex=0; Ia_ex < Jacnt; Ia_ex++) {
                   kl = *Iaij++;
                   I = *Iaridx++;
                   tval = *Iasgn++;
                   ijkl = INDEX(ij,kl);
                   VS = tval * tei[ijkl];
                   CprimeI0 = Cprime[I];

#ifdef USE_BLAS
                   C_DAXPY(jlen, VS, CprimeI0, 1, V, 1);
#else
                   for (J=0; J<jlen; J++) {
                       V[J] += VS * CprimeI0[J];
                     }
#endif

                 }


               /* scatter */
               for (J=0; J<jlen; J++) {
                   RJ = R[J];
                   S[Ia_idx][RJ] += V[J];
                 }

             } /* end loop over Ia */
         }
       
     } /* end loop over j */
   } /* end loop over i */
  for (i=0; i<nas; i++) free(thread_info[i]);
  free(thread);
   
}

/*
** S3_BLOCK_V_PTHREAD()
**
** Multithreaded component of the S3_BLOCK_V routine
**
** Calculate a block of the sigma3 vector in equation (9c) of
** Olsen, Roos, et al.  For non-diagonal blocks of s3
**
*/
void s3_block_v_pthread(void *threadarg)
{

  struct stringwr *Ia_local;
  unsigned int Ia_ex;
  int I, J, RJ, kl, ijkl;
  double tval, VS, *CprimeI0, *CI0;
  int Jacnt, *Iaij;
  unsigned int *Iaridx;
  signed char *Iasgn;
  struct pthreads_s3diag *thread_info_i;

  int nas, jlen, ij, Ja_list;
  double **S, **Cprime, *tei, *V, *Tptr;
  int *R, Ia_idx_local;
  int loop_var = 0;
  int work_units = 0;
  int thread_id;

  thread_info_i = (struct pthreads_s3diag *) threadarg;
  nas = thread_info_i->nas;
  jlen = thread_info_i->jlen;
  ij = thread_info_i->ij;
  Ja_list = thread_info_i->Ja_list;
  S = thread_info_i->S;
  Cprime = thread_info_i->Cprime;
  tei = thread_info_i->Tptr;
  R = thread_info_i->R;
  thread_id = thread_info_i->thread_id;
  Ia_idx_local = thread_info_i->Ia_idx_local;
  Ia_local = thread_info_i->Ia_local;
  
  V = init_array(jlen);

  /* loop over excitations E^a_{kl} from |A(I_a)> */
  Jacnt = Ia_local->cnt[Ja_list];
  Iaridx = Ia_local->ridx[Ja_list];
  Iasgn = Ia_local->sgn[Ja_list];
  Iaij = Ia_local->ij[Ja_list];

  work_units++;
      
  zero_arr(V, jlen);
         
  for (Ia_ex=0; Ia_ex < Jacnt; Ia_ex++) {
      kl = *Iaij++;
      I = *Iaridx++;
      tval = *Iasgn++;
      ijkl = INDEX(ij,kl);
      VS = tval * tei[ijkl];
      CprimeI0 = Cprime[I];

#ifdef USE_BLAS
      C_DAXPY(jlen, VS, CprimeI0, 1, V, 1);
#else
      for (J=0; J<jlen; J++) {
          V[J] += VS * CprimeI0[J];
        }
#endif
    }

  /* scatter */
  for (J=0; J<jlen; J++) {
      RJ = R[J];
      S[Ia_idx_local][RJ] += V[J];
    }

  free(V);
  /* return 0; */
  
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
/*  if(inum) {
      fprintf(outfile,"form_ilist: nas = %d\n", nas);
      fprintf(outfile,"form_ilist: jlen = %d\n", inum);
      fflush(outfile);
    }
*/
    return(inum);
}


/*
** S3_BLOCK_VDIAG_ROTF()
**
** Calculate a block of the sigma3 vector in equation (9c) of
** Olsen, Roos, et al.  For diagonal blocks of sigma.
**
** currently assumes that (ij|ij)'s have not been halved
** Try to get the Olsen vector version working....again!!!!
*/

void s3_block_vdiag_rotf(int *Cnt[2], int **Ij[2], int **Ridx[2],
      signed char **Sn[2], double **C, double **S, 
      double *tei, int nas, int nbs, int cnas,
      int Ib_list, int Ja_list, int Jb_list, int Ib_sym, int Jb_sym,
      double **Cprime, double *F, double *V, double *Sgn, int *L, int *R)
{
   int Ia_ex;
   int ij, i, j, kl, I, J, RJ;
   double tval, VS, *CprimeI0, *CI0;
   int jlen, Ia_idx, Jacnt, *Iaij, *orbsym, norbs;
   int *Iaridx;
   signed char *Iasgn;
   double *Tptr;
   
   norbs = CalcInfo.num_ci_orbs;
   orbsym = CalcInfo.orbsym + CalcInfo.num_fzc_orbs;

   /* loop over i, j */
   for (i=0; i<norbs; i++) {
     for (j=0; j<=i; j++) {
       if ((orbsym[i] ^ orbsym[j] ^ Jb_sym ^ Ib_sym) != 0) continue;
       ij = ioff[i] + j;
       jlen = form_ilist_rotf(Cnt[1], Ridx[1], Sn[1], Ij[1], 
          nbs, ij, L, R, Sgn);
       
       if (!jlen) continue;
       
       Tptr = tei + ioff[ij];
       
       /* gather operation */
       for (I=0; I<cnas; I++) {
         CprimeI0 = Cprime[I];
         CI0 = C[I];
         for (J=0; J<jlen; J++) {
           tval = Sgn[J];
           CprimeI0[J] = CI0[L[J]] * tval;
         }
       }


       /* loop over Ia */
       for (Ia_idx=0; Ia_idx<nas; Ia_idx++) {

         /* loop over excitations E^a_{kl} from |A(I_a)> */
         Jacnt = Cnt[0][Ia_idx];
         Iaridx = Ridx[0][Ia_idx];
         Iasgn = Sn[0][Ia_idx];
         Iaij = Ij[0][Ia_idx];
         
         zero_arr(V, jlen);
         
         /* rotf doesn't yet ensure kl's in order */
         for (Ia_ex=0; Ia_ex < Jacnt; Ia_ex++) {
           kl = *Iaij++;
           I = *Iaridx++;
           tval = *Iasgn++;
           if (kl > ij) continue;
           if (ij == kl) tval *= 0.5;
           VS = Tptr[kl] * tval;
           CprimeI0 = Cprime[I];
           
         #ifdef USE_BLAS
           C_DAXPY(jlen, VS, CprimeI0, 1, V, 1);
         #else
           for (J=0; J<jlen; J++) {
             V[J] += VS * CprimeI0[J];
           }
         #endif
         }

         
         /* scatter */
         for (J=0; J<jlen; J++) {
           RJ = R[J];
           S[Ia_idx][RJ] += V[J];
         }

       } /* end loop over Ia */

     } /* end loop over j */
   } /* end loop over i */
   
}              


/*
** S3_BLOCK_VROTF()
**
** Calculate a block of the sigma3 vector in equation (9c) of
** Olsen, Roos, et al.  For non-diagonal blocks of s3
**
*/
void s3_block_vrotf(int *Cnt[2], int **Ij[2], int **Ridx[2],
      signed char **Sn[2], double **C, double **S, 
      double *tei, int nas, int nbs, int cnas,
      int Ib_list, int Ja_list, int Jb_list, int Ib_sym, int Jb_sym,
      double **Cprime, double *F, double *V, double *Sgn, int *L, int *R)
{
   int Ia_ex;
   int ij, i, j, kl, ijkl, I, J, RJ;
   double tval, VS, *CprimeI0, *CI0;
   int jlen, Ia_idx, Jacnt, *Iaij, *orbsym, norbs;
   int *Iaridx;
   signed char *Iasgn;
   double *Tptr;
   
   norbs = CalcInfo.num_ci_orbs;
   orbsym = CalcInfo.orbsym + CalcInfo.num_fzc_orbs;

   /* loop over i, j */
   for (i=0; i<norbs; i++) {
     for (j=0; j<=i; j++) {
       if ((orbsym[i] ^ orbsym[j] ^ Jb_sym ^ Ib_sym) != 0) continue;
       ij = ioff[i] + j;
       jlen = form_ilist_rotf(Cnt[1], Ridx[1], Sn[1], Ij[1], 
          nbs, ij, L, R, Sgn);
       
       if (!jlen) continue;
       
       Tptr = tei + ioff[ij];
       
       /* gather operation */
       for (I=0; I<cnas; I++) {
         CprimeI0 = Cprime[I];
         CI0 = C[I];
         for (J=0; J<jlen; J++) {
           tval = Sgn[J];
           CprimeI0[J] = CI0[L[J]] * tval;
         }
       }


       /* loop over Ia */
       for (Ia_idx=0; Ia_idx<nas; Ia_idx++) {

         /* loop over excitations E^a_{kl} from |A(I_a)> */
         Jacnt = Cnt[0][Ia_idx];
         Iaridx = Ridx[0][Ia_idx];
         Iasgn = Sn[0][Ia_idx];
         Iaij = Ij[0][Ia_idx];
         
         zero_arr(V, jlen);
         
         for (Ia_ex=0; Ia_ex < Jacnt; Ia_ex++) {
           kl = *Iaij++;
           I = *Iaridx++;
           tval = *Iasgn++;
           ijkl = INDEX(ij,kl);
           VS = tval * tei[ijkl];
           CprimeI0 = Cprime[I];
           
         #ifdef USE_BLAS
           C_DAXPY(jlen, VS, CprimeI0, 1, V, 1);
         #else
           for (J=0; J<jlen; J++) {
             V[J] += VS * CprimeI0[J];
           }
         #endif
           
         }

         
         /* scatter */
         for (J=0; J<jlen; J++) {
           RJ = R[J];
           S[Ia_idx][RJ] += V[J];
         }

       } /* end loop over Ia */

     } /* end loop over j */
   } /* end loop over i */
   
}              


int form_ilist_rotf(int *Cnt, int **Ridx, signed char **Sn, int **Ij,
   int nas, int kl, int *L, int *R, double *Sgn)
{

   int inum=0, Ia_idx, Ia_ex, Iacnt, ij;
   int *Iaij;
   int *Iaridx;
   signed char *Iasgn;

   /* loop over Ia */
   for (Ia_idx=0; Ia_idx < nas; Ia_idx++) {

      /* loop over excitations E^a_{kl} from |A(I_a)> */

      Iacnt = Cnt[Ia_idx];
      if (!Iacnt) continue;
      Iaridx = Ridx[Ia_idx];
      Iasgn = Sn[Ia_idx];
      Iaij = Ij[Ia_idx];
      Ia_ex=0;
      for (Ia_ex=0; Ia_ex<Iacnt; Ia_ex++) {
         ij = *Iaij++;
         if (ij == kl) {
            *R++ = Ia_idx;
            *L++ = Iaridx[Ia_ex];
            *Sgn++ = (double) Iasgn[Ia_ex];
            inum++;
            }
         }
      }  /* end loop over Ia */

   return(inum);
}

}} // namespace psi::detci

