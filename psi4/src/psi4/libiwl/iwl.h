/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

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

#include <cstdio>
#include "psi4/libpsio/psio.h"
#include "config.h"
#include "psi4/psi4-dec.h"
namespace psi {

struct iwlbuf {
  int itap;                   /* tape number for input file */
  psio_address bufpos;        /* current page/offset */
  int ints_per_buf;           /* integrals per buffer */
  int bufszc;                 /* buffer size in characters (bytes) */
  double cutoff;              /* cutoff value for writing */
  int lastbuf;                /* is this the last IWL buffer? 1=yes,0=no */
  int inbuf;                  /* how many ints in current buffer? */
  int idx;                    /* index of integral in current buffer */
  Label *labels;              /* pointer to where integral values begin */
  Value *values;              /* integral values */
};


void iwl_buf_fetch(struct iwlbuf *Buf);
void iwl_buf_put(struct iwlbuf *Buf);

int iwl_rdone(int itap, const char *label, double *ints, int ntri, int erase,
              int printflg, std::string out);

void iwl_wrtone(int itap, const char *label, int ntri, double *onel_ints);

void iwl_rdtwo(int itap, double *ints, int *ioff, int norbs,
                      int nfzc, int nfzv, int printflg, std::string out);
void iwl_wrttwo(int itap, int nbfso, double *ints, int *ioff,
      double toler, int printflg, std::string out);
void sortbuf(struct iwlbuf *inbuf, struct iwlbuf *outbuf,
      double *ints, int fpq, int lpq, int *ioff, int *ioff2,
      int nbfso, int elbert, int intermediate, int no_pq_perm,
      int qdim, int add, int printflg, std::string out);
void sortbuf_pk(struct iwlbuf *Inbuf, int out_tape, int is_exch,
      double *ints, unsigned int fpq, unsigned int lpq, int *so2ind, int *so2sym, int *pksymoff,
      int printflg, std::string out);
void iwl_buf_init(struct iwlbuf *Buf, int intape, double cutoff,
      int oldfile, int readflag);
int iwl_buf_rd(struct iwlbuf *Buf, int target_pq, double *ints,
      int *ioff_lt, int *ioff_rt, int mp2, int printflg, std::string out);
int iwl_buf_rd_all(struct iwlbuf *Buf, double *ints,
      int *ioff_lt, int *ioff_rt, int mp2, int *ioff,
      int printflg, std::string out);
int iwl_buf_rd_all2(struct iwlbuf *Buf, double **ints,
		   int *ioff_lt, int *ioff_rt, int no_pq_perm, int *ioff,
		    int printflg, std::string out);
int iwl_buf_rd_all_act(struct iwlbuf *Buf, double *ints,
      int *ioff_lt, int *ioff_rt, int mp2, int *ioff,
      int fstact, int lstact, int printflg, std::string out);
int iwl_buf_rd_all_mp2r12a(struct iwlbuf *Buf, double *ints,
      int *ioff_lt, int *ioff_rt, int bra_ket_symm, int *ioff,
      int printflg, std::string out);
void iwl_buf_wrt_all(struct iwlbuf *Buf, int nbfso, double *ints,
      int *ioff, int printflg, std::string out);
void iwl_buf_wrt(struct iwlbuf *Buf, int p, int q, int pq, int pqsym,
      double *arr, int rmax, int *ioff, int *orbsym, int *firsti,
      int *lasti, int printflag, std::string out);
void iwl_buf_wrt_mp2(struct iwlbuf *Buf, int p, int q, int pq,
      int pqsym, double **arr, int rsym, int *firstr, int *lastr,
      int *firsts, int *lasts, int *occ, int *vir, int *ioff,
      int printflag, std::string out);
void iwl_buf_wrt_mp2r12a(struct iwlbuf *Buf, int p, int q, int pq,
      int pqsym, double **arr, int rsym, int *firstr, int *lastr,
      int *firsts, int *lasts, int *occ, int bra_ket_symm, int *ioff,
      int printflag, std::string out);
void iwl_buf_flush(struct iwlbuf *Buf, int lastbuf);
void iwl_buf_close(struct iwlbuf *Buf, int keep);
//void iwl_buf_toend(struct iwlbuf *Buf);
void iwl_buf_wrt_arr(struct iwlbuf *Buf, double *arr, int *p, int *q,
      int *r, int *s, long int size);
void iwl_buf_wrt_arr_SI(struct iwlbuf *Buf, double *arr,
      short int *p, short int *q, short int *r, short int *s, int size);
void iwl_buf_wrt_arr_SI_nocut(struct iwlbuf *Buf, double *arr,
      short int *p, short int *q, short int *r, short int *s, int size);
int iwl_buf_rd_arr(struct iwlbuf *Buf, int target_pq, double *ints,
      int *rlist, int *slist, int *size, int *ioff,
      int printflg, std::string out);
int iwl_buf_rd_arr2(struct iwlbuf *Buf, double *ints, int *plist,
      int *qlist, int *rlist, int *slist, int *size, int *ioff,
      int printflg, std::string out);
void iwl_buf_wrt_arr2(struct iwlbuf *Buf, double *arr, int p, int q,
      int *rlist, int *slist, int size, int printflag, std::string out);
void iwl_buf_wrt_mat(struct iwlbuf *Buf, int ptr, int qtr,
      double **mat, int rfirst, int rlast, int sfirst, int slast,
      int *reorder, int reorder_offset, int printflag, int *ioff,
      std::string out);
void iwl_buf_wrt_mat2(struct iwlbuf *Buf, int ptr, int qtr,
      double **mat, int rfirst, int rlast, int sfirst, int slast,
      int *reorder, int reorder_offset, int printflag, int *ioff,
      std::string out);
void iwl_buf_wrt_val(struct iwlbuf *Buf, int p, int q, int r, int s,
                     double value, int printflag, std::string out, int dirac);
void iwl_buf_wrt_val_SI(struct iwlbuf *Buf, short int p, short int q,
                     short int r, short int s, double value, int printflag,
                     std::string out, int dirac);

}

#endif /* end _psi_src_lib_libiwl_iwl_h */
