/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*! \file
    \ingroup TRANSQT
    \brief Enter brief description of file here 
*/
#ifndef _psi3_bin_transqt_yoshimine_h_
#define _psi3_bin_transqt_yoshimine_h_

namespace psi { namespace transqt {

/*
** YOSHIMINE.H 
** Function prototypes for Yoshimine sort object
**
** David Sherrill
** Center for Comptational Quantum Chemistry, UGA
** February 1995
*/

/* need to include iwl.h before including this file */

struct bucket {
   long int in_bucket;
   struct iwlbuf IWLBuf;
   int *p;
   int *q;
   int *r;
   int *s;
   double *val;
   int hi;
   int lo;
   };

struct yoshimine {
   int core_loads;
   int nbuckets;
   int *bucket_for_pq;
   unsigned long int bucketsize;
   struct bucket *buckets;
   int first_tmp_file;
   int pq_per_bucket;
   int bra_indices;
   int ket_indices;
   double cutoff;
   };

union psi_buffer {
  double *pki ;
  int *lbli ;
  unsigned char *lbl;
  double *val;
} ;

#ifdef YEXTERN
#undef YEXTERN
#define YEXTERN  
#else
#define YEXTERN extern
#endif

YEXTERN void yosh_init(struct yoshimine *YBuff, unsigned bra_indices, 
      unsigned ket_indices, long maxcor, long maxcord,
const int max_buckets,
      unsigned int first_tmp_file, double cutoff, std::string OutFileRMR);
YEXTERN void yosh_print(struct yoshimine *YBuff, std::string OutFileRMR); 
YEXTERN void yosh_init_buckets(struct yoshimine *YBuff);
YEXTERN void yosh_close_buckets(struct yoshimine *YBuff, int erase);
YEXTERN void yosh_rdtwo(struct yoshimine *YBuff, int itapERI, int del_tei_file, int *num_so,
      int nirreps, int *ioff, int elbert, int fzcflag, double *P,
      double *Hc, int matrix, int printflag, std::string OutFileRMR);
YEXTERN void yosh_rdtwo_uhf(struct yoshimine *YBuff, int itapERI, int del_tei_file, int *num_so,
      int nirreps, int *ioff, int elbert, int fzcflag, double *Pa, double *Pb,
      double *Hca, double *Hcb, int matrix, int printflag, std::string OutFileRMR);
YEXTERN void yosh_rdtwo_backtr(struct yoshimine *YBuff, int tei_file, 
      int *ioff, int symmetrize, int add_ref_pt, int del_tei_file, int prtflg, 
      std::string OutFileRMR);
YEXTERN void yosh_rdtwo_backtr_uhf(std::string, struct yoshimine *YBuff, int tei_file, 
      int *ioff, int swap_bk, int symm_pq, int del_tei_file, int prtflg, std::string OutFileRMR);
YEXTERN void flush_bucket(struct bucket *bptr, int lastbuf);
YEXTERN void yosh_sort(struct yoshimine *YBuff, int out_tape, int keep_bins,
      int *ioff, int *ioff2, int nbfso, int nbstri, 
      int elbert, int intermediate, int no_pq_perm, int qdim,
      int add, int print_lvl, std::string OutFileRMR);
YEXTERN void yosh_done(struct yoshimine *YBuff);
YEXTERN void yosh_flush(struct yoshimine *YBuff);
YEXTERN void yosh_wrt_arr(struct yoshimine *YBuff, int p, int q, int pq, 
   int pqsym, double *arr, int rmax, int *ioff, 
   int *orbsym, int *firsti, int *lasti, int sortby_rs, int printflag, 
   std::string OutFileRMR);
YEXTERN void yosh_wrt_arr2(struct yoshimine *YBuff, int size, double *arr,
   int p, int q, int *rlist, int *slist, int *ioff, int printflag,
   std::string OutFileRMR);
YEXTERN void yosh_wrt_arr_mp2(struct yoshimine *YBuff, int p, int q, int pq,
                      int pqsym, double **arr, int rsym, int *firstr,
                      int *lastr, int *firsts, int *lasts, int sortby_rs,
                      int ndocc, int nvirt, int *occ, int *vir, int *ioff3,
                      int printflag, std::string OutFileRMR);
YEXTERN void add_2pdm_ref_pt(struct yoshimine *YBuff,int *ioff,int prtflg,
                             std::string OutFileRMR);
YEXTERN void yosh_buff_put_val(struct yoshimine *YBuff, int *ioff, int pq,
                       int p, int q, int r, int s, double value, int prtflg,
                       std::string OutFileRMR);
YEXTERN void yosh_wrt_arr_mp2r12a(struct yoshimine *YBuff, int p, int q, int pq,
                          int pqsym, double **arr, int rsym, int *firstr,
                          int *lastr, int *firsts, int *lasts, int sortby_rs,
                          int *occ, int *ioff3,
                          int printflag, std::string OutFileRMR);
}} // end namespace psi::transqt
#endif // header guard
