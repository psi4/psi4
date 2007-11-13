/*
** IWL.H
** Header file for Integrals With Labels Library
**
** David Sherrill
** Center for Computational Quantum Chemistry, UGA
**
*/

#ifndef _psi_src_lib_libiwl_iwl_h_
#define _psi_src_lib_libiwl_iwl_h_

#include <stdio.h>
#include <libpsio/psio.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "config.h"

void iwl_buf_fetch(struct iwlbuf *Buf);
void iwl_buf_put(struct iwlbuf *Buf);

int iwl_rdone(int itap, char *label, double *ints, int ntri, int erase,
              int printflg, FILE *outfile);

void iwl_wrtone(int itap, char *label, int ntri, double *onel_ints);

void iwl_rdtwo(int itap, double *ints, int *ioff, int norbs,
                      int nfzc, int nfzv, int printflg, FILE *outfile);
void iwl_wrttwo(int itap, int nbfso, double *ints, int *ioff, 
      double toler, int printflg, FILE *outfile);
void sortbuf(struct iwlbuf *inbuf, struct iwlbuf *outbuf,
      double *ints, int fpq, int lpq, int *ioff, int *ioff2, 
      int nbfso, int elbert, int intermediate, int no_pq_perm, 
      int qdim, int add, int printflg, FILE *outfile); 
void iwl_buf_init(struct iwlbuf *Buf, int intape, double cutoff,
      int oldfile, int readflag);
int iwl_buf_rd(struct iwlbuf *Buf, int target_pq, double *ints,
      int *ioff_lt, int *ioff_rt, int mp2, int printflg, FILE *outfile);
int iwl_buf_rd_all(struct iwlbuf *Buf, double *ints,
      int *ioff_lt, int *ioff_rt, int mp2, int *ioff,
      int printflg, FILE *outfile);
int iwl_buf_rd_all2(struct iwlbuf *Buf, double **ints,
		   int *ioff_lt, int *ioff_rt, int no_pq_perm, int *ioff,
		    int printflg, FILE *outfile);
int iwl_buf_rd_all_act(struct iwlbuf *Buf, double *ints,
      int *ioff_lt, int *ioff_rt, int mp2, int *ioff,
      int fstact, int lstact, int printflg, FILE *outfile);
int iwl_buf_rd_all_mp2r12a(struct iwlbuf *Buf, double *ints,
      int *ioff_lt, int *ioff_rt, int bra_ket_symm, int *ioff,
      int printflg, FILE *outfile);
void iwl_buf_wrt_all(struct iwlbuf *Buf, int nbfso, double *ints, 
      int *ioff, int printflg, FILE *outfile);
void iwl_buf_wrt(struct iwlbuf *Buf, int p, int q, int pq, int pqsym,
      double *arr, int rmax, int *active, int *ioff, int *orbsym, int *firsti, 
      int *lasti, int sortby_rs, int printflag, FILE *outfile);
void iwl_buf_wrt_mp2(struct iwlbuf *Buf, int p, int q, int pq,
      int pqsym, double **arr, int rsym, int *firstr, int *lastr, 
      int *firsts, int *lasts, int *occ, int *vir, int *ioff, 
      int printflag, FILE *outfile);
void iwl_buf_wrt_mp2r12a(struct iwlbuf *Buf, int p, int q, int pq,
      int pqsym, double **arr, int rsym, int *firstr, int *lastr, 
      int *firsts, int *lasts, int *occ, int bra_ket_symm, int *ioff, 
      int printflag, FILE *outfile);
void iwl_buf_flush(struct iwlbuf *Buf, int lastbuf);
void iwl_buf_close(struct iwlbuf *Buf, int keep);
void iwl_buf_toend(struct iwlbuf *Buf);
void iwl_buf_wrt_arr(struct iwlbuf *Buf, double *arr, int *p, int *q,
      int *r, int *s, long int size);
void iwl_buf_wrt_arr_SI(struct iwlbuf *Buf, double *arr, 
      short int *p, short int *q, short int *r, short int *s, int size);
void iwl_buf_wrt_arr_SI_nocut(struct iwlbuf *Buf, double *arr,
      short int *p, short int *q, short int *r, short int *s, int size);
int iwl_buf_rd_arr(struct iwlbuf *Buf, int target_pq, double *ints,
      int *rlist, int *slist, int *size, int *ioff,
      int printflg, FILE *outfile);
int iwl_buf_rd_arr2(struct iwlbuf *Buf, double *ints, int *plist,
      int *qlist, int *rlist, int *slist, int *size, int *ioff,
      int printflg, FILE *outfile);
void iwl_buf_wrt_arr2(struct iwlbuf *Buf, double *arr, int p, int q, 
      int *rlist, int *slist, int size, int printflag, FILE *outfile);
void iwl_buf_wrt_mat(struct iwlbuf *Buf, int ptr, int qtr,
      double **mat, int rfirst, int rlast, int sfirst, int slast,
      int *reorder, int reorder_offset, int printflag, int *ioff,
      FILE *outfile);
void iwl_buf_wrt_mat2(struct iwlbuf *Buf, int ptr, int qtr,
      double **mat, int rfirst, int rlast, int sfirst, int slast,
      int *reorder, int reorder_offset, int printflag, int *ioff,
      FILE *outfile);
void iwl_buf_wrt_val(struct iwlbuf *Buf, int p, int q, int r, int s,
                     double value, int printflag, FILE *outfile, int dirac);
void iwl_buf_wrt_val_SI(struct iwlbuf *Buf, short int p, short int q,
                     short int r, short int s, double value, int printflag,
                     FILE *outfile, int dirac);
#ifdef __cplusplus
}
#endif

#endif /* end _psi_src_lib_libiwl_iwl_h */
