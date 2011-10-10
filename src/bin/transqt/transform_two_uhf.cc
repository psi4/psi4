/*! \file
    \ingroup TRANSQT
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>

#include "MOInfo.h"
#include "Params.h"
#include "globals.h"
#include "yoshimine.h"

namespace psi {
extern FILE* outfile;
namespace transqt {

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))

double fzc_energy_uhf(int nbfso, int *sosym, double *Pa, double *Pb,
                      double *Hca, double *Hcb,double *H,
                      int *first_so, int *ioff);

/* transform_two_uhf(): Carry out the transformation of the SO-basis
** two-electron integrals into the AA, BB, and AB MO-basis classes.
** In these comments, I use p,q,r,s to denote SO-basis indices,
** I,J,K,L to denote alpha-MO indices, and i,j,k,l to denote beta-MO
** indices:
**
** (1) Yoshimine presort the SO-basis integrals into a supermatrix,
** with p>=q, r>=s, but with all pq and rs (just like the usual
** presort for integral transformations in transform_two.c).
**
** (2) Half-transform the integrals into the alpha-MO basis to
** form a (pq,IJ) matrix and Yoshimine sort the result into the
** transposed form, (IJ,pq).
**
** (3) Complete the second half-transformation for the AA and AB
** integral lists: (IJ,KL) and (IJ,kl).
**
** (4) Half-transform the SO-basis integrals into the beta-MO basis to
** form a (pq,ij) matrix and Yoshimine sort the result into the
** transposed form, (ij,pq).
**
** (5) Complete the second half-transformation for the BB integral
** list: (ij,kl).
**
** TDC
** June 2001
*/

void transform_two_uhf(void)
{
  int nirreps;
  long int maxcor, maxcord;
  int max_buckets, first_tmp_file;
  double tolerance;
  int print_lvl;
  int fzc_offset;
  int *src_first, *src_last, *src_orbspi, *src_orbsym, src_orbs, src_ntri;
  int *dst_first, *dst_last, *dst_orbspi, *dst_orbsym, dst_orbs, dst_ntri;
  int i, j, ifirst, ilast, jfirst, jlast, ij;
  int k, l, ksym, lsym, kfirst, klast, lfirst, llast, kl, klsym;
  int p, q, psym, qsym, pfirst, plast, qfirst, qlast, pq, pqsym;
  int r, s, rsym, ssym, rfirst, rlast, sfirst, slast, rs;
  int P, Q, R, S, K, L, ktr, ltr;
  struct yoshimine YBuffP;
  struct yoshimine YBuffJ;
  struct iwlbuf PBuff;
  struct iwlbuf JBuff;
  struct iwlbuf MBuff;
  double *P_block, *J_block;
  double **A, **B, ***C;
  int A_cols, B_cols, *C_colspi;
  int *reorder_alpha, *reorder_beta;

  maxcor = params.maxcor;
  maxcord = params.maxcord;
  max_buckets = params.max_buckets;
  first_tmp_file = params.first_tmp_file;
  tolerance = params.tolerance;
  print_lvl = params.print_lvl;

  nirreps = moinfo.nirreps;
  reorder_alpha = moinfo.order_alpha;
  reorder_beta = moinfo.order_beta;

  src_first = moinfo.first_so;
  src_last = moinfo.last_so;
  src_orbspi = moinfo.sopi;
  src_orbsym = moinfo.sosym;
  src_orbs = moinfo.nso;
  src_ntri = src_orbs * (src_orbs+1)/2;

  dst_first = (params.fzc && !params.do_all_tei) ? moinfo.fstact : moinfo.first;
  dst_last = (!params.do_all_tei) ? moinfo.lstact : moinfo.last;
  dst_orbspi = (!params.do_all_tei) ? moinfo.active : moinfo.orbspi;
  dst_orbsym = moinfo.orbsym;
  dst_orbs = (!params.do_all_tei) ? (moinfo.nmo - moinfo.nfzv - moinfo.nfzc) : moinfo.nmo;
  dst_ntri = moinfo.nmo * (moinfo.nmo+1)/2;
  C_colspi = dst_orbspi;

  fzc_offset = (params.fzc && !params.do_all_tei) ? moinfo.nfzc : 0;

  /** Presort the two-electron integrals **/

  if (params.print_lvl) {
    fprintf(outfile, "\n\tPre-sorting two-electron integrals...\n\n");
    fflush(outfile);
  }

  yosh_init(&YBuffP, src_ntri, src_ntri, maxcor, maxcord,
            max_buckets, first_tmp_file, tolerance, outfile);

  if (print_lvl > 1) {
    fprintf(outfile, "\tPresort Yoshimine Parameters");
    yosh_print(&YBuffP, outfile);
    fprintf(outfile, "\n");
    fflush(outfile);
  }

  yosh_init_buckets(&YBuffP);

  yosh_rdtwo_uhf(&YBuffP, params.src_tei_file,
                 params.delete_src_tei, src_orbspi, nirreps, ioff, 0,
                 params.fzc && moinfo.nfzc, moinfo.fzc_density_alpha,
                 moinfo.fzc_density_beta, moinfo.fzc_operator_alpha,
                 moinfo.fzc_operator_beta, 1, (print_lvl > 5), outfile);

  yosh_close_buckets(&YBuffP, 0);

  yosh_sort(&YBuffP, params.presort_file, 0, ioff, NULL, src_orbs, src_ntri,
            0, 1, 0, 0, 1, (print_lvl > 5), outfile);

  yosh_done(&YBuffP);

  /** Pre-sort complete **/

  /** Compute the frozen core energy **/
  moinfo.efzc = 0.0;
  if (params.fzc && moinfo.nfzc) {
    moinfo.efzc =
      fzc_energy_uhf(moinfo.nso, moinfo.sosym, moinfo.fzc_density_alpha,
                     moinfo.fzc_density_beta, moinfo.fzc_operator_alpha,
                     moinfo.fzc_operator_beta, moinfo.oe_ints,
                     moinfo.first_so, ioff);
    free(moinfo.fzc_density_alpha);
    free(moinfo.fzc_density_beta);
  }

  /* Write efzc to chkpt file */
  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_efzc(moinfo.efzc);
  chkpt_close();

  if (params.print_lvl)
    fprintf(outfile, "\n\tFrozen core energy = %20.15lf\n", moinfo.efzc);

  P_block = init_array(src_ntri);
  J_block = init_array(MAX0(src_ntri,dst_ntri));

  A_cols = MAX0(src_orbs, dst_orbs);
  B_cols = dst_orbs;
  A = block_matrix(A_cols,A_cols);
  B = block_matrix(src_orbs, dst_orbs);

  /** AA/AB First-half transformation **/

  fprintf(outfile, "\n\tBeginning AA/AB two-electron transform...\n");

  C = moinfo.evects_alpha;

  iwl_buf_init(&PBuff, params.presort_file, tolerance, 1, 1);
  yosh_init(&YBuffJ, dst_ntri, src_ntri, maxcor, maxcord, max_buckets,
            first_tmp_file, tolerance, outfile);

  yosh_init_buckets(&YBuffJ);

  if (print_lvl > 1) {
    fprintf(outfile, "\tHalf-transform Yoshimine Parameters");
    yosh_print(&YBuffJ, outfile);
    fprintf(outfile, "\n");
    fflush(outfile);
  }

  for (psym=0; psym < nirreps; psym++) {
    pfirst = src_first[psym];
    plast = src_last[psym];
    for (p=pfirst; p <= plast; p++) {
      for (qsym=0; qsym <= psym; qsym++) {
        qfirst = src_first[qsym];
        qlast = src_last[qsym];
        pqsym = psym^qsym;
        for (q=qfirst; (q<=qlast) && (q <= p); q++) {
          pq = ioff[p] + q;

          zero_arr(P_block,src_ntri);
          iwl_buf_rd(&PBuff, pq, P_block, ioff, ioff, 0, (print_lvl>4),
                     outfile);

          for (rsym=0; rsym < nirreps; rsym++) {
            rfirst = src_first[rsym];
            rlast = src_last[rsym];
            kfirst = dst_first[rsym];
            klast = dst_last[rsym];
            ssym = pqsym^rsym;
            if (ssym <= rsym) {
              sfirst = src_first[ssym];
              slast = src_last[ssym];
              lfirst = dst_first[ssym];
              llast = dst_last[ssym];
              if(ssym == rsym) {
                for(r=rfirst,R=0; r <= rlast; r++,R++) {
                  for(s=sfirst,S=0; (s <= slast) && (s <= r);
                      s++,S++) {
                    rs = INDEX(r,s);
                    A[R][S] = P_block[rs];
                    A[S][R] = P_block[rs];
                  }
                }
              }
              else {
                for (r=rfirst,R=0; r <= rlast; r++,R++) {
                  for (s=sfirst,S=0; s <= slast; s++,S++) {
                    rs = INDEX(r,s);
                    A[R][S] = P_block[rs];
                  }
                }
              }

              if (C_colspi[ssym] > 0)
                C_DGEMM('n', 'n', src_orbspi[rsym],
                        dst_orbspi[ssym], src_orbspi[ssym], 1.0,
                        A[0], A_cols, C[ssym][0], C_colspi[ssym], 0.0,
                        B[0], B_cols);
              if (C_colspi[rsym] > 0)
                C_DGEMM('t', 'n', dst_orbspi[rsym],
                        dst_orbspi[ssym], src_orbspi[rsym], 1.0,
                        C[rsym][0], C_colspi[rsym], B[0], B_cols, 0.0,
                        A[0], A_cols);

              zero_arr(J_block, dst_ntri);
              for (k=kfirst,K=0; k <= klast; k++,K++) {
                for (l=lfirst,L=0; (l <= llast) && (l <= k); l++,L++) {
                  kl = ioff[k] + l;
                  J_block[kl] = A[K][L];
                }
              }

              yosh_wrt_arr(&YBuffJ, p, q, pq, pqsym, J_block,
                           moinfo.nmo, ioff, dst_orbsym, dst_first, dst_last,
                           1, (print_lvl > 4), outfile);

            }
          }

        }
      }
    }
  }

  iwl_buf_close(&PBuff, 1);

  if (params.print_lvl)
    fprintf(outfile, "\tSorting AA/AB half-transformed integrals...\n");

  yosh_flush(&YBuffJ);
  yosh_close_buckets(&YBuffJ, 0);
  yosh_sort(&YBuffJ, params.jfile, 0, ioff, NULL, src_orbs, src_ntri, 0, 1, 0, 0,
            1, (print_lvl > 5), outfile);
  yosh_done(&YBuffJ);

  if (print_lvl) {
    fprintf(outfile, "\tFinished AA/AB half-transformation...\n");
    fflush(outfile);
  }

  /** AA Second-half transformation **/

  C = moinfo.evects_alpha;

  iwl_buf_init(&JBuff, params.jfile, tolerance, 1, 1);
  iwl_buf_init(&MBuff, params.aa_mfile, tolerance, 0, 0);

  for (ksym=0; ksym < nirreps; ksym++) {
    kfirst = dst_first[ksym];
    klast = dst_last[ksym];
    for (k=kfirst; k <= klast; k++) {
      for (lsym=0; lsym <= ksym; lsym++) {
        lfirst = dst_first[lsym];
        llast = dst_last[lsym];
        klsym = ksym^lsym;
        for (l=lfirst; (l <= llast) && (l <= k); l++) {
          kl = ioff[k] + l;
          if (!params.backtr) {
            ktr = reorder_alpha[k] - fzc_offset;
            ltr = reorder_alpha[l] - fzc_offset;
          }

          zero_arr(J_block, src_ntri);
          iwl_buf_rd(&JBuff, kl, J_block, ioff, ioff, 0, 0, outfile);

          for (psym=0; psym < nirreps; psym++) {
            pfirst = src_first[psym];
            plast = src_last[psym];
            ifirst = dst_first[psym];
            ilast = dst_last[psym];
            qsym = klsym^psym;
            if (qsym <= psym) {
              qfirst = src_first[qsym];
              qlast = src_last[qsym];
              jfirst = dst_first[qsym];
              jlast = dst_last[qsym];
              for (p=pfirst,P=0; p <= plast; p++,P++) {
                for (q=qfirst,Q=0; q <= qlast; q++,Q++) {
                  pq = INDEX(p,q);
                  A[P][Q] = J_block[pq];
                }
              }


              if (C_colspi[qsym] > 0)
                C_DGEMM('n', 'n', src_orbspi[psym],
                        dst_orbspi[qsym], src_orbspi[qsym], 1.0,
                        A[0], A_cols, C[qsym][0], C_colspi[qsym], 0.0,
                        B[0], B_cols);
              if (C_colspi[psym] > 0)
                C_DGEMM('t', 'n', dst_orbspi[psym],
                        dst_orbspi[qsym], src_orbspi[psym], 1.0,
                        C[psym][0], C_colspi[psym], B[0], B_cols, 0.0,
                        A[0], A_cols);

              iwl_buf_wrt_mat(&MBuff, ktr, ltr, A,
                              ifirst, ilast, jfirst, jlast,
                              reorder_alpha, fzc_offset, params.print_te_ints,
                              ioff, outfile);
            }
          }
        }
      }
    }
  }

  iwl_buf_close(&JBuff, 1);

  iwl_buf_flush(&MBuff, 1);
  iwl_buf_close(&MBuff, 1);

  if (print_lvl) {
    fprintf(outfile, "\n\tAA Transformation finished.\n");
    fprintf(outfile, "\tTwo-electron AA integrals written to file%d.\n",params.aa_mfile);
    fflush(outfile);
  }

  /** AB Second-half transformation **/

  C = moinfo.evects_beta;

  iwl_buf_init(&JBuff, params.jfile, tolerance, 1, 1);
  iwl_buf_init(&MBuff, params.ab_mfile, tolerance, 0, 0);

  for (ksym=0; ksym < nirreps; ksym++) {
    kfirst = dst_first[ksym];
    klast = dst_last[ksym];
    for (k=kfirst; k <= klast; k++) {
      for (lsym=0; lsym <= ksym; lsym++) {
        lfirst = dst_first[lsym];
        llast = dst_last[lsym];
        klsym = ksym^lsym;
        for (l=lfirst; (l <= llast) && (l <= k); l++) {
          kl = ioff[k] + l;
          if (!params.backtr) {
            ktr = reorder_alpha[k] - fzc_offset;
            ltr = reorder_alpha[l] - fzc_offset;
          }

          zero_arr(J_block, src_ntri);
          iwl_buf_rd(&JBuff, kl, J_block, ioff, ioff, 0, 0, outfile);

          for (psym=0; psym < nirreps; psym++) {
            pfirst = src_first[psym];
            plast = src_last[psym];
            ifirst = dst_first[psym];
            ilast = dst_last[psym];
            qsym = klsym^psym;
            if (qsym <= psym) {
              qfirst = src_first[qsym];
              qlast = src_last[qsym];
              jfirst = dst_first[qsym];
              jlast = dst_last[qsym];
              for (p=pfirst,P=0; p <= plast; p++,P++) {
                for (q=qfirst,Q=0; q <= qlast; q++,Q++) {
                  pq = INDEX(p,q);
                  A[P][Q] = J_block[pq];
                }
              }


              if (C_colspi[qsym] > 0)
                C_DGEMM('n', 'n', src_orbspi[psym],
                        dst_orbspi[qsym], src_orbspi[qsym], 1.0,
                        A[0], A_cols, C[qsym][0], C_colspi[qsym], 0.0,
                        B[0], B_cols);
              if (C_colspi[psym] > 0)
                C_DGEMM('t', 'n', dst_orbspi[psym],
                        dst_orbspi[qsym], src_orbspi[psym], 1.0,
                        C[psym][0], C_colspi[psym], B[0], B_cols, 0.0,
                        A[0], A_cols);

              iwl_buf_wrt_mat2(&MBuff, ktr, ltr, A,
                               ifirst, ilast, jfirst, jlast,
                               reorder_beta, fzc_offset, params.print_te_ints,
                               ioff, outfile);
            }
          }
        }
      }
    }
  }

  iwl_buf_close(&JBuff, 0);

  iwl_buf_flush(&MBuff, 1);
  iwl_buf_close(&MBuff, 1);

  if (print_lvl) {
    fprintf(outfile, "\n\tAB Transformation finished.\n");
    fprintf(outfile, "\tTwo-electron AB integrals written to file%d.\n",params.ab_mfile);
    fflush(outfile);
  }

  /** BB First-half transformation **/

  fprintf(outfile, "\n\tBeginning BB two-electron transform...\n");

  C = moinfo.evects_beta;

  iwl_buf_init(&PBuff, params.presort_file, tolerance, 1, 1);
  yosh_init(&YBuffJ, dst_ntri, src_ntri, maxcor, maxcord, max_buckets,
            first_tmp_file, tolerance, outfile);

  yosh_init_buckets(&YBuffJ);

  if (print_lvl > 1) {
    fprintf(outfile, "\tHalf-transform Yoshimine Parameters");
    yosh_print(&YBuffJ, outfile);
    fprintf(outfile, "\n");
    fflush(outfile);
  }

  for (psym=0; psym < nirreps; psym++) {
    pfirst = src_first[psym];
    plast = src_last[psym];
    for (p=pfirst; p <= plast; p++) {
      for (qsym=0; qsym <= psym; qsym++) {
        qfirst = src_first[qsym];
        qlast = src_last[qsym];
        pqsym = psym^qsym;
        for (q=qfirst; (q<=qlast) && (q <= p); q++) {
          pq = ioff[p] + q;

          zero_arr(P_block,src_ntri);
          iwl_buf_rd(&PBuff, pq, P_block, ioff, ioff, 0, (print_lvl>4),
                     outfile);

          for (rsym=0; rsym < nirreps; rsym++) {
            rfirst = src_first[rsym];
            rlast = src_last[rsym];
            kfirst = dst_first[rsym];
            klast = dst_last[rsym];
            ssym = pqsym^rsym;
            if (ssym <= rsym) {
              sfirst = src_first[ssym];
              slast = src_last[ssym];
              lfirst = dst_first[ssym];
              llast = dst_last[ssym];
              if(ssym == rsym) {
                for(r=rfirst,R=0; r <= rlast; r++,R++) {
                  for(s=sfirst,S=0; (s <= slast) && (s <= r);
                      s++,S++) {
                    rs = INDEX(r,s);
                    A[R][S] = P_block[rs];
                    A[S][R] = P_block[rs];
                  }
                }
              }
              else {
                for (r=rfirst,R=0; r <= rlast; r++,R++) {
                  for (s=sfirst,S=0; s <= slast; s++,S++) {
                    rs = INDEX(r,s);
                    A[R][S] = P_block[rs];
                  }
                }
              }

              if (C_colspi[ssym] > 0)
                C_DGEMM('n', 'n', src_orbspi[rsym],
                        dst_orbspi[ssym], src_orbspi[ssym], 1.0,
                        A[0], A_cols, C[ssym][0], C_colspi[ssym], 0.0,
                        B[0], B_cols);
              if (C_colspi[rsym] > 0)
                C_DGEMM('t', 'n', dst_orbspi[rsym],
                        dst_orbspi[ssym], src_orbspi[rsym], 1.0,
                        C[rsym][0], C_colspi[rsym], B[0], B_cols, 0.0,
                        A[0], A_cols);

              zero_arr(J_block, dst_ntri);
              for (k=kfirst,K=0; k <= klast; k++,K++) {
                for (l=lfirst,L=0; (l <= llast) && (l <= k); l++,L++) {
                  kl = ioff[k] + l;
                  J_block[kl] = A[K][L];
                }
              }

              yosh_wrt_arr(&YBuffJ, p, q, pq, pqsym, J_block,
                           moinfo.nmo, ioff, dst_orbsym, dst_first, dst_last,
                           1, (print_lvl > 4), outfile);

            }
          }

        }
      }
    }
  }

  iwl_buf_close(&PBuff, params.keep_presort);

  if (params.print_lvl)
    fprintf(outfile, "\tSorting BB half-transformed integrals...\n");

  yosh_flush(&YBuffJ);
  yosh_close_buckets(&YBuffJ, 0);
  yosh_sort(&YBuffJ, params.jfile, 0, ioff, NULL, src_orbs, src_ntri, 0, 1, 0, 0,
            1, (print_lvl > 5), outfile);
  yosh_done(&YBuffJ);

  if (print_lvl) {
    fprintf(outfile, "\tFinished BB half-transformation...\n");
    fflush(outfile);
  }

  /** BB Second-half transformation **/

  C = moinfo.evects_beta;

  iwl_buf_init(&JBuff, params.jfile, tolerance, 1, 1);
  iwl_buf_init(&MBuff, params.bb_mfile, tolerance, 0, 0);

  for (ksym=0; ksym < nirreps; ksym++) {
    kfirst = dst_first[ksym];
    klast = dst_last[ksym];
    for (k=kfirst; k <= klast; k++) {
      for (lsym=0; lsym <= ksym; lsym++) {
        lfirst = dst_first[lsym];
        llast = dst_last[lsym];
        klsym = ksym^lsym;
        for (l=lfirst; (l <= llast) && (l <= k); l++) {
          kl = ioff[k] + l;
          if (!params.backtr) {
            ktr = reorder_beta[k] - fzc_offset;
            ltr = reorder_beta[l] - fzc_offset;
          }

          zero_arr(J_block, src_ntri);
          iwl_buf_rd(&JBuff, kl, J_block, ioff, ioff, 0, 0, outfile);

          for (psym=0; psym < nirreps; psym++) {
            pfirst = src_first[psym];
            plast = src_last[psym];
            ifirst = dst_first[psym];
            ilast = dst_last[psym];
            qsym = klsym^psym;
            if (qsym <= psym) {
              qfirst = src_first[qsym];
              qlast = src_last[qsym];
              jfirst = dst_first[qsym];
              jlast = dst_last[qsym];
              for (p=pfirst,P=0; p <= plast; p++,P++) {
                for (q=qfirst,Q=0; q <= qlast; q++,Q++) {
                  pq = INDEX(p,q);
                  A[P][Q] = J_block[pq];
                }
              }


              if (C_colspi[qsym] > 0)
                C_DGEMM('n', 'n', src_orbspi[psym],
                        dst_orbspi[qsym], src_orbspi[qsym], 1.0,
                        A[0], A_cols, C[qsym][0], C_colspi[qsym], 0.0,
                        B[0], B_cols);
              if (C_colspi[psym] > 0)
                C_DGEMM('t', 'n', dst_orbspi[psym],
                        dst_orbspi[qsym], src_orbspi[psym], 1.0,
                        C[psym][0], C_colspi[psym], B[0], B_cols, 0.0,
                        A[0], A_cols);

              iwl_buf_wrt_mat(&MBuff, ktr, ltr, A,
                              ifirst, ilast, jfirst, jlast,
                              reorder_beta, fzc_offset, params.print_te_ints,
                              ioff, outfile);
            }
          }
        }
      }
    }
  }

  iwl_buf_close(&JBuff, 0);

  iwl_buf_flush(&MBuff, 1);
  iwl_buf_close(&MBuff, 1);

  if (print_lvl) {
    fprintf(outfile, "\n\tBB Transformation finished.\n");
    fprintf(outfile, "\tTwo-electron BB integrals written to file%d.\n",params.bb_mfile);
    fflush(outfile);
  }

  free(P_block);
  free(J_block);
  free_block(A);
  free_block(B);
}

}} // end namespace psi::transqt

