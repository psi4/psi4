/*! \file
    \ingroup TRANSQT
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include "MOInfo.h"
#include "Params.h"
#include "globals.h"
#include "yoshimine.h"
#include "backsort.h"

namespace psi {
extern FILE* outfile;
namespace transqt {

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))


void transform_two_mp2(int maxcor, int maxcord);

double fzc_energy(int nbfso, int *sosym, double *P, double *Hc, double *H,
      int *first_so, int *ioff);
void transform_two_uhf(void);
void transform_two_backtr_uhf(void);

void transform_two(void)
{

  int p,q,pq;
  int r,s,rs,pqrs;
  int k,l,kl,ktr,ltr;
  int i,j,ij,ijkl;
  int pfirst, qfirst, rfirst, sfirst;
  int plast, qlast, rlast, slast;
  int kfirst,lfirst,ifirst,jfirst;
  int klast,llast,ilast,jlast;
  int psym,qsym,rsym,ssym,pqsym,pqrsym;
  int ksym,lsym,klsym,klpsym;
  int P,Q,I,J,R,S,K,L;
  int nmo,nirreps,src_ntri,dst_ntri,src_orbs,dst_orbs;
  int *src_first, *dst_first, *src_last, *dst_last;
  int *src_orbspi, *dst_orbspi;
  int fzc_offset, *reorder, *dst_orbsym, *src_orbsym;
  int errcod,print_integrals;
  int lmax,count;
  char **labels;
  double **A, **B, *J_block, *P_block, *mbuf;
  double ***C;
  struct yoshimine YBuffP;
  struct yoshimine YBuffJ;
  struct iwlbuf PBuff;
  struct iwlbuf JBuff;
  struct iwlbuf MBuff;
  long int maxcor,maxcord;
  int max_buckets, print_lvl, first_tmp_file, presort_file;
  int keep_presort, jfile, keep_half_tf, mfile;
  double tolerance;
  int A_cols, B_cols, *C_cols;
  struct iwlbuf *twopdm_out;

  /* Special code for full UHF transformations */
  if(params.ref == "UHF") {
    if(params.backtr) transform_two_backtr_uhf();
    else transform_two_uhf();
    return;
  }

  nmo = moinfo.nmo;
  nirreps = params.backtr ? moinfo.backtr_nirreps : moinfo.nirreps;

  if (params.backtr) {  /* if we do a so-called "backtransform" */
    src_first = moinfo.backtr_mo_first;
    dst_first = moinfo.backtr_ao_first;
    src_last = moinfo.backtr_mo_lstact;
    dst_last = moinfo.backtr_ao_last;
    src_orbspi = moinfo.backtr_mo_active;
    dst_orbspi = moinfo.backtr_ao_orbspi;
    dst_orbsym = moinfo.backtr_ao_orbsym;
    fzc_offset = 0;
    src_orbs = moinfo.nmo - moinfo.nfzv;
    dst_orbs = moinfo.nao;
    C_cols = src_orbspi;
  }
  else { /* if we do the normal "forwards" transform */
    src_first = moinfo.first_so;
    src_last = moinfo.last_so;
    src_orbspi = moinfo.sopi;
    src_orbsym = moinfo.sosym;
    dst_first = (params.fzc && !params.do_all_tei)?
                moinfo.fstact : moinfo.first;
    dst_last = (params.do_all_tei) ? moinfo.last : moinfo.lstact;
    dst_orbspi = (params.do_all_tei) ? moinfo.orbspi : moinfo.active;
    dst_orbsym = moinfo.orbsym;
    fzc_offset = (params.fzc && !params.do_all_tei) ? moinfo.nfzc : 0;
    src_orbs = moinfo.nso;
    if (params.do_all_tei) dst_orbs = moinfo.nmo;
    else dst_orbs = moinfo.nmo - moinfo.nfzv - moinfo.nfzc;
    C_cols = dst_orbspi;
  }

  print_integrals = params.print_te_ints;
  labels = moinfo.labels;
  C = moinfo.evects;
  max_buckets = params.max_buckets;
  print_lvl = params.print_lvl;
  first_tmp_file = params.first_tmp_file;
  presort_file = params.presort_file;
  keep_presort = params.keep_presort;
  jfile = params.jfile;
  keep_half_tf = params.keep_half_tf;
  mfile = params.mfile;
  tolerance = params.tolerance;
  maxcor = params.maxcor;
  maxcord = params.maxcord;
  reorder = moinfo.order;

  /* compute the triangle sizes needed by the sorting algorithms...note
   *  that we can't use dst_orbs (at least for forward transforms) because
   *  we are in the Pitzer order, where not all frozen virts are up top
   */
  if (!params.backtr) {
    src_ntri = src_orbs * (src_orbs+1) / 2;
    dst_ntri = moinfo.nmo*(moinfo.nmo+1)/2;
  }
  else {
    src_ntri = src_orbs * (src_orbs+1) / 2;
    dst_ntri = moinfo.nao*(moinfo.nao+1)/2;
  }

  /** Pre-sort the two-electron Integrals **/

  if (params.print_lvl) {
    fprintf(outfile, "\n\tPre-sorting two-electron ints...\n\n");
    fflush(outfile);
  }

  yosh_init(&YBuffP, src_ntri, src_ntri, maxcor, maxcord, max_buckets,
            first_tmp_file, tolerance, outfile);

  if (print_lvl > 1) {
    fprintf(outfile, "\tPresort");
    yosh_print(&YBuffP, outfile);
    fprintf(outfile, "\n");
    fflush(outfile);
  }

  yosh_init_buckets(&YBuffP);

  if (!params.backtr) {
    yosh_rdtwo(&YBuffP, params.src_tei_file,
               params.delete_src_tei, src_orbspi, nirreps, ioff, 0,
               params.fzc && moinfo.nfzc, moinfo.fzc_density,
               moinfo.fzc_operator, 1, (print_lvl > 5), outfile);

  }
  else {
    yosh_rdtwo_backtr(&YBuffP, params.src_tei_file, ioff, 1,
                 params.tpdm_add_ref, params.delete_src_tei,
                 (print_lvl > 4), outfile);
  }

  yosh_close_buckets(&YBuffP, 0);

  yosh_sort(&YBuffP, presort_file, 0, ioff, NULL, src_orbs, src_ntri,
            0, 1, 0, 0, 1, (print_lvl > 5), outfile);

  yosh_done(&YBuffP);  /* Pre-transform complete */

  /* get the frozen core energy */
  moinfo.efzc = 0.0;
  if (params.fzc && moinfo.nfzc) {
    moinfo.efzc = fzc_energy(moinfo.nso, moinfo.sosym, moinfo.fzc_density,
                             moinfo.fzc_operator, moinfo.oe_ints,
                             moinfo.first_so, ioff);
    free(moinfo.fzc_density);
  }

  /* Write efzc to chkpt file */
  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_efzc(moinfo.efzc);
  chkpt_close();

  if (params.print_lvl)
    fprintf(outfile, "\n\tFrozen core energy = %20.15lf\n", moinfo.efzc);

  if (params.print_lvl) {
    fprintf(outfile, "\tTransforming two-electron ints...\n\n");
    fflush(outfile);
  }

  iwl_buf_init(&PBuff, presort_file, tolerance, 1, 1);
  yosh_init(&YBuffJ, dst_ntri, src_ntri, maxcor, maxcord, max_buckets,
            first_tmp_file, tolerance, outfile);

  yosh_init_buckets(&YBuffJ);

  if (print_lvl > 1) {
    fprintf(outfile, "\tHalf-transform");
    yosh_print(&YBuffJ, outfile);
    fprintf(outfile, "\n");
    fflush(outfile);
  }

  A_cols = MAX0(src_orbs,dst_orbs);
  A = block_matrix(A_cols, A_cols);
  B = block_matrix(src_orbs,dst_orbs);
  B_cols = dst_orbs;

  P_block = init_array(src_ntri);
  J_block = init_array(MAX0(src_ntri,dst_ntri));

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
          iwl_buf_rd(&PBuff, pq, P_block, ioff, ioff, 0, (params.print_lvl > 4),
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

#ifdef USE_BLAS
              if (C_cols[ssym] > 0)
                C_DGEMM('n', params.backtr ? 't' : 'n', src_orbspi[rsym],
                        dst_orbspi[ssym], src_orbspi[ssym], 1.0,
                        A[0], A_cols, C[ssym][0], C_cols[ssym], 0.0,
                        B[0], B_cols);
              if (C_cols[rsym] > 0)
                C_DGEMM(params.backtr ? 'n' : 't', 'n', dst_orbspi[rsym],
                        dst_orbspi[ssym], src_orbspi[rsym], 1.0,
                        C[rsym][0], C_cols[rsym], B[0], B_cols, 0.0,
                        A[0], A_cols);
#else
              mmult(A, 0, C[ssym], params.backtr ? 1 : 0, B, 0,
                    src_orbspi[rsym], src_orbspi[ssym], dst_orbspi[ssym], 0);
              mmult(C[rsym], params.backtr ? 0 : 1, B, 0, A, 0,
                    dst_orbspi[rsym], src_orbspi[rsym], dst_orbspi[ssym], 0);
#endif

              zero_arr(J_block, dst_ntri);
              for (k=kfirst,K=0; k <= klast; k++,K++) {
                for (l=lfirst,L=0; (l <= llast) && (l <= k); l++,L++) {
                  kl = ioff[k] + l;
                  J_block[kl] = A[K][L];
                }
              }
              yosh_wrt_arr(&YBuffJ, p, q, pq, pqsym, J_block,
                           params.backtr ? moinfo.nao : nmo, ioff, dst_orbsym,
                           dst_first, dst_last, 1, (print_lvl > 4), outfile);
            }
          }
        }
      }
    }
  }

  if (params.print_lvl)
    fprintf(outfile, "\tSorting half-transformed integrals...\n");

  iwl_buf_close(&PBuff, keep_presort);
  free(P_block);
  yosh_flush(&YBuffJ);
  yosh_close_buckets(&YBuffJ, 0);
  yosh_sort(&YBuffJ, jfile, 0, ioff, NULL, src_orbs, src_ntri, 0, 1, 0, 0,
            1, (print_lvl > 5), outfile);
  yosh_done(&YBuffJ);

  if (params.print_lvl) {
    fprintf(outfile, "\tFinished half-transform...\n");
    fprintf(outfile, "\tWorking on second half...\n");
    fflush(outfile);
  }

  /** Second half of transformation **/

  iwl_buf_init(&JBuff, jfile, tolerance, 1, 1);
  if (!params.backtr) iwl_buf_init(&MBuff, mfile, tolerance, 0, 0);
  else {
    backsort_prep(0);
    twopdm_out = (struct iwlbuf *) malloc(nbuckets * sizeof(struct iwlbuf));
    for(p=0; p < nbuckets; p++)
        iwl_buf_init(&twopdm_out[p], first_tmp_file+p, tolerance, 0, 0);
  }

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
            ktr = reorder[k] - fzc_offset;
            ltr = reorder[l] - fzc_offset;
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

#ifdef USE_BLAS
              if (C_cols[qsym] > 0)
                C_DGEMM('n', params.backtr ? 't' : 'n', src_orbspi[psym],
                        dst_orbspi[qsym], src_orbspi[qsym], 1.0,
                        A[0], A_cols, C[qsym][0], C_cols[qsym], 0.0,
                        B[0], B_cols);
              if (C_cols[psym] > 0)
                C_DGEMM(params.backtr ? 'n' : 't', 'n', dst_orbspi[psym],
                        dst_orbspi[qsym], src_orbspi[psym], 1.0,
                        C[psym][0], C_cols[psym], B[0], B_cols, 0.0,
                        A[0], A_cols);
#else
              mmult(A, 0, C[qsym], params.backtr ? 1 : 0, B, 0,
                    src_orbspi[psym], src_orbspi[qsym], dst_orbspi[qsym], 0);
              mmult(C[psym], params.backtr ? 0 : 1, B, 0, A, 0,
                    dst_orbspi[psym], src_orbspi[psym], dst_orbspi[qsym], 0);
#endif

              if (!params.backtr)
                iwl_buf_wrt_mat(&MBuff, ktr, ltr, A,
                                ifirst, ilast, jfirst, jlast,
                                reorder, fzc_offset, print_integrals,
                                ioff, outfile);

              else
                backsort_write(k, l, A, ifirst, ilast, jfirst, jlast,
                               print_integrals, outfile, twopdm_out, 0);
            }
          }
        }
      }
    }
  }

  free(J_block);
  free_block(A);
  free_block(B);
  iwl_buf_close(&JBuff, keep_half_tf);

  /* Need to flush last buffer, else it's not written to disk */
  if (!params.backtr) {
    iwl_buf_flush(&MBuff, 1);
    iwl_buf_close(&MBuff, 1);
  }
  else {
    for(p=0; p < nbuckets; p++) {
        iwl_buf_flush(&twopdm_out[p], 1);
        iwl_buf_close(&twopdm_out[p], 1);
      }
    free(twopdm_out);
    backsort(first_tmp_file, tolerance, 0);
  }

  if (params.print_lvl) {
    fprintf(outfile, "\n\tTransformation finished.\n");
    fprintf(outfile, "\tTwo-electron integrals written to file%d.\n",mfile);
    fflush(outfile);
  }

}

}} // end namespace psi::transqt

