/*! \file
    \ingroup TRANSQT
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <psifiles.h>

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

/* transform_two_backtr_uhf(): Carry out the backtransformation of the
** AA, BB, and AB two-particle density matrices (twopdms) for
** UHF-based correlated wave functions:
**
** (1) Yoshimine presort the AA and AB G(pqrs) twopdms into two
** supermatrices, with p>=q, r>=s, but with all pq and rs (just like
** the usual presort for integral transformations in transform_two.c
** and transform_two_uhf.c).
**
** (2) With a single set of p,q,r,s loops, half-transform the AA and
** AB twopdms to the AO basis, forming G(pq,ij), and Yoshimine sort
** the results into the transposed form, G(ij,pq).  Note that since
** both twopdms are now ordered (alpha-MO pq x AO ij), they are
** combined during the Yoshimine transposition.
**
** (3) Yoshimine presort the BB G(pqrs) as in step (1).
**
** (4) Half-transform the BB twopdm into G(pq,ij) and Yoshimine sort
** the result into G(ij,pq).
**
** (5) With a single set of i,j,p,q loops, complete the second
** half-transform for both the AA+AB and BB twopdms.  Since the final
** targets are entirely in the AO basis, they are combined during the
** final dump to disk.
**
** (6) The total AO-basis twopdm is Yoshimine sorted into the
** "shell-quartet" ordering needed by the derivative integral code.
** (See backsort.c.)
**
** NB that this structure is (necessarily) very different from the
** forward UHF-based transform of the two-electron integrals in
** transform_two_uhf.c, so I could not use the same code for forward
** and backwards transformations as was possible for the
** RHF/ROHF-based transforms in transform_two.c.
**
** TDC
** January 2003
*/

void transform_two_backtr_uhf(void)
{
  int nirreps;
  long int maxcor, maxcord;
  int max_buckets, first_tmp_file;
  double tolerance;
  int print_lvl;
  int src_orbs, src_ntri, *src_first, *src_last, *src_orbspi;
  int dst_orbs, dst_ntri, *dst_first, *dst_last, *dst_orbspi, *dst_orbsym;
  int p, q, psym, qsym, pfirst, plast, qfirst, qlast, pq, pqsym;
  int r, s, rsym, ssym, rfirst, rlast, sfirst, slast, rs;
  int i, j, isym, jsym, ifirst, ilast, jfirst, jlast, ij;
  int k, l, ksym, lsym, kfirst, klast, lfirst, llast, kl, klsym;
  int P, Q, R, S, K, L, I, J;
  struct yoshimine YBuffP;
  struct yoshimine YBuffJ;
  struct iwlbuf PAA_Buff, PAB_Buff, PBB_Buff;
  struct iwlbuf JA_Buff, JB_Buff;
  struct iwlbuf *twopdm_out;
  double *PAA_block, *PAB_block, *PBB_block, *J_block, *JA_block, *JB_block;
  double **A_AA, **A_AB, **A_BB, **B, ***CA, ***CB;
  int A_cols, B_cols, *C_colspi;
  std::string AA = "AA";
  std::string AB = "AB";
  std::string BB = "BB";

  int ktr, ltr, *reorder, ntei;
  struct iwlbuf MBuff, JBuff;
  double *ints, **dens, energy;

  int dim, pqrs;
  double *ab_block, **ab_dens1, **ab_dens2, **ab_ints;

  double AA_norm=0.0, AB_norm=0.0, BB_norm=0.0, AA_norm1=0.0, BB_norm1=0.0;

  nirreps = moinfo.backtr_nirreps;
  maxcor = params.maxcor;
  maxcord = params.maxcord;
  max_buckets = params.max_buckets;
  first_tmp_file = params.first_tmp_file;
  tolerance = params.tolerance;
  print_lvl = params.print_lvl;

  src_first = moinfo.backtr_mo_first;
  src_last = moinfo.backtr_mo_lstact;
  src_orbspi = moinfo.backtr_mo_active;
  src_orbs = moinfo.nmo - moinfo.nfzv;
  src_ntri = src_orbs * (src_orbs+1)/2;

  dst_first = moinfo.backtr_ao_first;
  dst_last = moinfo.backtr_ao_last;
  dst_orbspi = moinfo.backtr_ao_orbspi;
  dst_orbsym = moinfo.backtr_ao_orbsym;
  dst_orbs = moinfo.nao;
  dst_ntri = dst_orbs * (dst_orbs+1)/2;

  C_colspi = src_orbspi;

  /** Presort the AA twopdm **/

  if(params.print_lvl) {
    fprintf(outfile, "\n\tPre-sorting AA two-particle density...\n\n"); fflush(outfile);
  }

  yosh_init(&YBuffP, src_ntri, src_ntri, maxcor, maxcord, max_buckets,
            first_tmp_file, tolerance, outfile);
  if(print_lvl > 1) { yosh_print(&YBuffP, outfile); fflush(outfile); }
  yosh_init_buckets(&YBuffP);
  yosh_rdtwo_backtr_uhf(AA, &YBuffP, PSIF_MO_AA_TPDM, ioff, 1, 1, 1, 0, outfile);
  yosh_close_buckets(&YBuffP, 0);
  yosh_sort(&YBuffP, PSIF_AA_PRESORT, 0, ioff, NULL, src_orbs, src_ntri, 0, 1, 0, 0, 1, 0, outfile);
  yosh_done(&YBuffP);

  /** Presort the AB twopdm **/

  if(params.print_lvl) {
    fprintf(outfile, "\n\tPre-sorting AB two-particle density...\n\n"); fflush(outfile);
  }

  yosh_init(&YBuffP, src_ntri, src_ntri, maxcor, maxcord, max_buckets,
            first_tmp_file, tolerance, outfile);
  if(print_lvl > 1) { yosh_print(&YBuffP, outfile); fflush(outfile); }
  yosh_init_buckets(&YBuffP);
  yosh_rdtwo_backtr_uhf(AB, &YBuffP, PSIF_MO_AB_TPDM, ioff, 0, 1, 1, 0, outfile);
  yosh_close_buckets(&YBuffP, 0);
  yosh_sort(&YBuffP, PSIF_AB_PRESORT, 0, ioff, NULL, src_orbs, src_ntri, 0, 1, 0, 0, 1, 0, outfile);
  yosh_done(&YBuffP);

  A_cols = MAX0(src_orbs, dst_orbs);
  A_AA = block_matrix(A_cols,A_cols);
  A_AB = block_matrix(A_cols,A_cols);
  A_BB = block_matrix(A_cols,A_cols);
  B_cols = dst_orbs;
  B = block_matrix(src_orbs, dst_orbs);

  CA = moinfo.evects_alpha;
  CB = moinfo.evects_beta;

  /* Try an in-core backtr of the AA density */
  /*
  dim = MAX0(src_ntri, dst_ntri);

  ab_dens1 = block_matrix(dim, dim);
  ab_dens2 = block_matrix(dim, dim);

  ab_block = init_array(dim);
  iwl_buf_init(&JBuff, PSIF_AA_PRESORT, tolerance, 1, 1);
  for(p=0, pq=0; p < src_orbs; p++) {
    for(q=0; q <= p; q++, pq++) {

      zero_arr(ab_block, src_ntri);
      iwl_buf_rd(&JBuff, pq, ab_block, ioff, ioff, 0, 0, outfile);

      for(r=0,rs=0; r < src_orbs; r++) {
        for(s=0; s <= r; s++,rs++) {
          ab_dens1[pq][rs] = ab_block[rs];
        }
      }
    }
  }
  iwl_buf_close(&JBuff, 1);
  free(ab_block);
  */

  /*
    fprintf(outfile, "\n\tMO-basis AA SCF Twopdm:\n");
    print_mat(ab_dens1, dim, dim, outfile);
  */

  /* Check energy in the MO basis */
  /*
  ntei = src_ntri * (src_ntri + 1)/2;
  ints = init_array(ntei);
  iwl_buf_init(&JBuff, PSIF_MO_AA_TEI, tolerance, 1, 0);
  iwl_buf_rd_all(&JBuff, ints, ioff, ioff, 0, ioff, 0, outfile);
  iwl_buf_close(&JBuff, 1);
  energy = 0.0;
  for(p=0; p < src_orbs; p++)
    for(q=0; q < src_orbs ; q++) {
      pq = INDEX(p,q);
      for(r=0; r < src_orbs; r++) {
        for(s=0; s < src_orbs; s++) {
          rs = INDEX(r,s);
          pqrs = INDEX(pq,rs);
          energy += ab_dens1[pq][rs] * ints[pqrs];
        }
      }
    }
  fprintf(outfile, "\n\tAA energy from MO-twopdm: %20.14f\n", energy);
  free(ints);


  for(p=0, pq=0; p < src_orbs; p++) {
    for(q=0; q <= p; q++, pq++) {

      for(r=0,rs=0; r < src_orbs; r++) {
        for(s=0; s <= r; s++,rs++) {
          A_AB[r][s] = A_AB[s][r] = ab_dens1[pq][rs];
        }
      }

      C_DGEMM('n','t', src_orbs, dst_orbs, src_orbs, 1.0, A_AB[0], A_cols,
              CA[0][0], src_orbs, 0.0, B[0], B_cols);
      zero_mat(A_AB, A_cols, A_cols);
      C_DGEMM('n','n', dst_orbs, dst_orbs, src_orbs, 1.0, CA[0][0], src_orbs,
              B[0], B_cols, 0.0, A_AB[0], A_cols);

      for(r=0, rs=0; r < dst_orbs; r++) {
        for(s=0; s <= r; s++,rs++) {
          ab_dens1[pq][rs] = A_AB[r][s];
          AA_norm += A_AB[r][s] * A_AB[r][s];
        }
      }
    }
  }

  fprintf(outfile, "\n\tAA_norm (1st half) = %20.15f\n", AA_norm);
  AA_norm = 0.0;

  for(p=0, pq=0; p < src_orbs; p++) {
    for(q=0; q <= p; q++, pq++) {
      for(r=0,rs=0; r < dst_orbs; r++) {
        for(s=0; s <= r; s++,rs++) {
          ab_dens2[rs][pq] = ab_dens1[pq][rs];
        }
      }
    }
  }


  for(r=0,rs=0; r < dst_orbs; r++) {
    for(s=0; s <= r; s++,rs++) {

      for(p=0, pq=0; p < src_orbs; p++) {
        for(q=0; q <= p; q++,pq++) {
          A_AB[p][q] = A_AB[q][p] = ab_dens2[rs][pq];
        }
      }

      C_DGEMM('n','t',src_orbs,dst_orbs,src_orbs,1.0,A_AB[0], A_cols,
              CA[0][0], src_orbs, 0.0, B[0], B_cols);
      zero_mat(A_AB, A_cols, A_cols);
      C_DGEMM('n','n',dst_orbs,dst_orbs,src_orbs,1.0, CA[0][0], src_orbs,
              B[0], B_cols, 0.0, A_AB[0], A_cols);

      for(p=0,pq=0; p < dst_orbs; p++) {
        for(q=0; q <= p; q++,pq++) {
          ab_dens2[rs][pq] = A_AB[p][q];
          AA_norm += A_AB[p][q] * A_AB[p][q];
        }
      }

    }
  }

  fprintf(outfile, "\n\tAA_norm (2nd half) = %20.15f\n", AA_norm);
  AA_norm = 0.0;
  */

  /* compute the AA energy to test above backtr */
  /*
  ntei = dst_ntri * (dst_ntri + 1)/2;
  ints = init_array(ntei);
  iwl_buf_init(&JBuff, PSIF_SO_TEI, tolerance, 1, 0);
  iwl_buf_rd_all(&JBuff, ints, ioff, ioff, 0, ioff, 0, outfile);
  iwl_buf_close(&JBuff, 1);

  energy = 0.0;
  for(p=0; p < dst_orbs; p++) {
    for(q=0; q < dst_orbs; q++) {
      pq = INDEX(p,q);
      for(r=0; r < dst_orbs; r++) {
        for(s=0; s < dst_orbs; s++) {
          rs = INDEX(r,s);
          pqrs = INDEX(pq,rs);
          energy += ab_dens2[pq][rs] * ints[pqrs];
        }
      }
    }
  }
  fprintf(outfile, "\n\tAA energy from AO-twopdm: %20.14f\n", energy);
  free(ints);

  free_block(ab_dens1);
  free_block(ab_dens2);
  */


  /* Try an in-core backtr of the AB density */
  /*
  dim = MAX0(src_ntri, dst_ntri);

  ab_dens1 = block_matrix(dim, dim);
  ab_dens2 = block_matrix(dim, dim);

  ab_block = init_array(dim);
  iwl_buf_init(&JBuff, PSIF_AB_PRESORT, tolerance, 1, 1);
  for(p=0, pq=0; p < src_orbs; p++) {
    for(q=0; q <= p; q++, pq++) {

      zero_arr(ab_block, src_ntri);
      iwl_buf_rd(&JBuff, pq, ab_block, ioff, ioff, 0, 0, outfile);

      for(r=0,rs=0; r < src_orbs; r++) {
        for(s=0; s <= r; s++,rs++) {
          ab_dens1[pq][rs] = ab_block[rs];
        }
      }
    }
  }
  iwl_buf_close(&JBuff, 1);
  free(ab_block);
  */

  /*
    fprintf(outfile, "\n\tMO-basis AB SCF Twopdm:\n");
    print_mat(ab_dens1, dim, dim, outfile);
  */

  /* Check energy in the MO basis */
  /*
  ntei = src_ntri * (src_ntri + 1)/2;
  ab_ints = block_matrix(dim,dim);
  iwl_buf_init(&JBuff, PSIF_MO_AB_TEI, tolerance, 1, 0);
  iwl_buf_rd_all2(&JBuff, ab_ints, ioff, ioff, 0, ioff, 0, outfile);
  iwl_buf_close(&JBuff, 1);
  energy = 0.0;
  for(p=0; p < src_orbs; p++) {
    for(q=0; q < src_orbs; q++) {
      pq = INDEX(p,q);
      for(r=0; r < src_orbs; r++) {
        for(s=0; s < src_orbs; s++) {
          rs = INDEX(r,s);

          pqrs = INDEX(pq,rs);

          energy += ab_dens1[pq][rs] * ab_ints[pq][rs];
        }
      }
    }
  }
  fprintf(outfile, "\n\tAB energy from MO-twopdm: %20.14f\n", energy);
  free_block(ab_ints);


  for(p=0, pq=0; p < src_orbs; p++) {
    for(q=0; q <= p; q++, pq++) {

      for(r=0,rs=0; r < src_orbs; r++) {
        for(s=0; s <= r; s++,rs++) {
          A_AB[r][s] = A_AB[s][r] = ab_dens1[pq][rs];
        }
      }

      C_DGEMM('n','t', src_orbs, dst_orbs, src_orbs, 1.0, A_AB[0], A_cols,
              CB[0][0], src_orbs, 0.0, B[0], B_cols);
      zero_mat(A_AB, A_cols, A_cols);
      C_DGEMM('n','n', dst_orbs, dst_orbs, src_orbs, 1.0, CB[0][0], src_orbs,
              B[0], B_cols, 0.0, A_AB[0], A_cols);

      for(r=0, rs=0; r < dst_orbs; r++) {
        for(s=0; s <= r; s++,rs++) {
          ab_dens1[pq][rs] = A_AB[r][s];
        }
      }
    }
  }

  for(p=0, pq=0; p < src_orbs; p++) {
    for(q=0; q <= p; q++, pq++) {
      for(r=0,rs=0; r < dst_orbs; r++) {
        for(s=0; s <= r; s++,rs++) {
          ab_dens2[rs][pq] = ab_dens1[pq][rs];
        }
      }
    }
  }


  for(r=0,rs=0; r < dst_orbs; r++) {
    for(s=0; s <= r; s++,rs++) {

      for(p=0, pq=0; p < src_orbs; p++) {
        for(q=0; q <= p; q++,pq++) {
          A_AB[p][q] = A_AB[q][p] = ab_dens2[rs][pq];
        }
      }

      C_DGEMM('n','t',src_orbs,dst_orbs,src_orbs,1.0,A_AB[0], A_cols,
              CA[0][0], src_orbs, 0.0, B[0], B_cols);
      zero_mat(A_AB, A_cols, A_cols);
      C_DGEMM('n','n',dst_orbs,dst_orbs,src_orbs,1.0, CA[0][0], src_orbs,
              B[0], B_cols, 0.0, A_AB[0], A_cols);

      for(p=0,pq=0; p < dst_orbs; p++) {
        for(q=0; q <= p; q++,pq++) {
          ab_dens2[rs][pq] = A_AB[p][q];
        }
      }

    }
  }
  */

  /* compute to the AB energy to test above backtr */
  /*
  ntei = dst_ntri * (dst_ntri + 1)/2;
  ints = init_array(ntei);
  iwl_buf_init(&JBuff, PSIF_SO_TEI, tolerance, 1, 0);
  iwl_buf_rd_all(&JBuff, ints, ioff, ioff, 0, ioff, 0, outfile);
  iwl_buf_close(&JBuff, 1);

  energy = 0.0;
  for(p=0; p < dst_orbs; p++) {
    for(q=0; q < dst_orbs; q++) {
      pq = INDEX(p,q);
      for(r=0; r < dst_orbs; r++) {
        for(s=0; s < dst_orbs; s++) {
          rs = INDEX(r,s);

          pqrs = INDEX(pq,rs);

          energy += ab_dens2[pq][rs] * ints[pqrs];
        }
      }
    }
  }
  fprintf(outfile, "\n\tAB energy from AO-twopdm: %20.14f\n", energy);
  free(ints);

  free_block(ab_dens1);
  free_block(ab_dens2);
  */

  /** AA/AB First-half transformation **/
  fprintf(outfile, "\n\tBeginning AA/AB twopdm transform...\n");

  PAA_block = init_array(src_ntri);
  PAB_block = init_array(src_ntri);
  iwl_buf_init(&PAA_Buff, PSIF_AA_PRESORT, tolerance, 1, 1);
  iwl_buf_init(&PAB_Buff, PSIF_AB_PRESORT, tolerance, 1, 1);

  J_block = init_array(MAX0(src_ntri,dst_ntri));
  yosh_init(&YBuffJ, dst_ntri, src_ntri, maxcor, maxcord, max_buckets,
            first_tmp_file, tolerance, outfile);
  if(print_lvl > 1) { yosh_print(&YBuffJ, outfile); fflush(outfile); }
  yosh_init_buckets(&YBuffJ);

  for(psym=0; psym < nirreps; psym++) {
    pfirst = src_first[psym];
    plast = src_last[psym];
    for (p=pfirst; p <= plast; p++) {
      for (qsym=0; qsym <= psym; qsym++) {
        qfirst = src_first[qsym];
        qlast = src_last[qsym];
        pqsym = psym^qsym;
        for (q=qfirst; (q<=qlast) && (q <= p); q++) {
          pq = ioff[p] + q;

          zero_arr(PAA_block, src_ntri);
          zero_arr(PAB_block, src_ntri);

          iwl_buf_rd(&PAA_Buff, pq, PAA_block, ioff, ioff, 0, 0, outfile);
          iwl_buf_rd(&PAB_Buff, pq, PAB_block, ioff, ioff, 0, 0, outfile);

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
                    A_AA[R][S] = PAA_block[rs];
                    A_AA[S][R] = PAA_block[rs];
                    A_AB[R][S] = PAB_block[rs];
                    A_AB[S][R] = PAB_block[rs];
                  }
                }
              }
              else {
                for (r=rfirst,R=0; r <= rlast; r++,R++) {
                  for (s=sfirst,S=0; s <= slast; s++,S++) {
                    rs = INDEX(r,s);
                    A_AA[R][S] = PAA_block[rs];
                    A_AB[R][S] = PAB_block[rs];
                  }
                }
              }

              /** AA half-transform for current pq **/
              if(C_colspi[ssym] > 0)
                C_DGEMM('n','t',src_orbspi[rsym],dst_orbspi[ssym],src_orbspi[ssym],1.0,
                        A_AA[0], A_cols,CA[ssym][0],C_colspi[ssym],0.0,B[0],B_cols);
              zero_mat(A_AA, A_cols, A_cols);
              if(C_colspi[rsym] > 0)
                C_DGEMM('n','n',dst_orbspi[rsym],dst_orbspi[ssym],src_orbspi[rsym],1.0,
                        CA[rsym][0],C_colspi[rsym],B[0],B_cols,0.0,A_AA[0],A_cols);

              /** AB half-transform for current pq **/
              if(C_colspi[ssym] > 0)
                C_DGEMM('n','t',src_orbspi[rsym],dst_orbspi[ssym],src_orbspi[ssym],1.0,
                        A_AB[0], A_cols,CB[ssym][0],C_colspi[ssym],0.0,B[0],B_cols);
              zero_mat(A_AB, A_cols, A_cols);
              if(C_colspi[rsym] > 0)
                C_DGEMM('n','n',dst_orbspi[rsym],dst_orbspi[ssym],src_orbspi[rsym],1.0,
                        CB[rsym][0],C_colspi[rsym],B[0],B_cols,0.0,A_AB[0],A_cols);

              /* collect the results and sum AA and AB contributions into J_block */
              zero_arr(J_block, dst_ntri);
              for(k=kfirst,K=0; k <= klast; k++,K++) {
                for (l=lfirst,L=0; (l <= llast) && (l <= k); l++,L++) {
                  kl = ioff[k] + l;
                  J_block[kl] = A_AA[K][L] + A_AB[K][L];
                }
              }

              yosh_wrt_arr(&YBuffJ, p, q, pq, pqsym, J_block,
                           moinfo.nao, ioff, dst_orbsym, dst_first, dst_last, 1, 0, outfile);

            }
          }
        }
      }
    }
  }

  iwl_buf_close(&PAA_Buff, 0);
  iwl_buf_close(&PAB_Buff, 0);
  free(PAA_block);
  free(PAB_block);

  if (params.print_lvl) {
    fprintf(outfile, "\tSorting AA/AB half-transformed twopdm...\n"); fflush(outfile);
  }
  yosh_flush(&YBuffJ);
  yosh_close_buckets(&YBuffJ, 0);
  yosh_sort(&YBuffJ, params.jfile, 0, ioff, NULL, src_orbs, src_ntri, 0, 1, 0, 0, 1, 0, outfile);
  yosh_done(&YBuffJ);
  free(J_block);
  if (print_lvl) {
    fprintf(outfile, "\tFinished AA/AB half-transformation...\n"); fflush(outfile);
  }

  /** Presort the BB twopdm **/

  if(params.print_lvl) {
    fprintf(outfile, "\n\tPre-sorting BB two-particle density...\n\n"); fflush(outfile);
  }

  yosh_init(&YBuffP, src_ntri, src_ntri, maxcor, maxcord, max_buckets,
            first_tmp_file, tolerance, outfile);
  if(print_lvl > 1) { yosh_print(&YBuffP, outfile); fflush(outfile); }
  yosh_init_buckets(&YBuffP);
  yosh_rdtwo_backtr_uhf(BB, &YBuffP, PSIF_MO_BB_TPDM, ioff, 1, 1, 1, 0, outfile);
  yosh_close_buckets(&YBuffP, 0);
  yosh_sort(&YBuffP, PSIF_BB_PRESORT, 0, ioff, NULL, src_orbs, src_ntri, 0, 1, 0, 0, 1, 0, outfile);
  yosh_done(&YBuffP);

  /* Try an in-core backtr of the BB density */
  /*
  dim = MAX0(src_ntri, dst_ntri);

  ab_dens1 = block_matrix(dim, dim);
  ab_dens2 = block_matrix(dim, dim);

  ab_block = init_array(dim);
  iwl_buf_init(&JBuff, PSIF_BB_PRESORT, tolerance, 1, 1);
  for(p=0, pq=0; p < src_orbs; p++) {
    for(q=0; q <= p; q++, pq++) {

      zero_arr(ab_block, src_ntri);
      iwl_buf_rd(&JBuff, pq, ab_block, ioff, ioff, 0, 0, outfile);

      for(r=0,rs=0; r < src_orbs; r++) {
        for(s=0; s <= r; s++,rs++) {
          ab_dens1[pq][rs] = ab_block[rs];
        }
      }
    }
  }
  iwl_buf_close(&JBuff, 1);
  free(ab_block);
  */

  /*
    fprintf(outfile, "\n\tMO-basis BB SCF Twopdm:\n");
    print_mat(ab_dens1, dim, dim, outfile);
  */

  /* Check energy in the MO basis */
  /*
  ntei = src_ntri * (src_ntri + 1)/2;
  ints = init_array(ntei);
  iwl_buf_init(&JBuff, PSIF_MO_BB_TEI, tolerance, 1, 0);
  iwl_buf_rd_all(&JBuff, ints, ioff, ioff, 0, ioff, 0, outfile);
  iwl_buf_close(&JBuff, 1);
  energy = 0.0;
  for(p=0; p < src_orbs; p++) {
    for(q=0; q < src_orbs; q++) {
      pq = INDEX(p,q);
      for(r=0; r < src_orbs; r++) {
        for(s=0; s < src_orbs; s++) {
          rs = INDEX(r,s);

          pqrs = INDEX(pq,rs);

          energy += ab_dens1[pq][rs] * ints[pqrs];
        }
      }
    }
  }
  fprintf(outfile, "\n\tBB energy from MO-twopdm: %20.14f\n", energy);
  free(ints);


  for(p=0, pq=0; p < src_orbs; p++) {
    for(q=0; q <= p; q++, pq++) {

      for(r=0,rs=0; r < src_orbs; r++) {
        for(s=0; s <= r; s++,rs++) {
          A_AB[r][s] = A_AB[s][r] = ab_dens1[pq][rs];
        }
      }

      C_DGEMM('n','t', src_orbs, dst_orbs, src_orbs, 1.0, A_AB[0], A_cols,
              CB[0][0], src_orbs, 0.0, B[0], B_cols);
      zero_mat(A_AB, A_cols, A_cols);
      C_DGEMM('n','n', dst_orbs, dst_orbs, src_orbs, 1.0, CB[0][0], src_orbs,
              B[0], B_cols, 0.0, A_AB[0], A_cols);

      for(r=0, rs=0; r < dst_orbs; r++) {
        for(s=0; s <= r; s++,rs++) {
          ab_dens1[pq][rs] = A_AB[r][s];
          BB_norm += A_AB[r][s] * A_AB[r][s];
        }
      }
    }
  }

  fprintf(outfile, "\n\tBB_norm (1st half) = %20.15f\n", BB_norm);
  BB_norm = 0.0;

  for(p=0, pq=0; p < src_orbs; p++) {
    for(q=0; q <= p; q++, pq++) {
      for(r=0,rs=0; r < dst_orbs; r++) {
        for(s=0; s <= r; s++,rs++) {
          ab_dens2[rs][pq] = ab_dens1[pq][rs];
        }
      }
    }
  }


  for(r=0,rs=0; r < dst_orbs; r++) {
    for(s=0; s <= r; s++,rs++) {

      for(p=0, pq=0; p < src_orbs; p++) {
        for(q=0; q <= p; q++,pq++) {
          A_AB[p][q] = A_AB[q][p] = ab_dens2[rs][pq];
        }
      }

      C_DGEMM('n','t',src_orbs,dst_orbs,src_orbs,1.0,A_AB[0], A_cols,
              CB[0][0], src_orbs, 0.0, B[0], B_cols);
      zero_mat(A_AB, A_cols, A_cols);
      C_DGEMM('n','n',dst_orbs,dst_orbs,src_orbs,1.0, CB[0][0], src_orbs,
              B[0], B_cols, 0.0, A_AB[0], A_cols);

      for(p=0,pq=0; p < dst_orbs; p++) {
        for(q=0; q <= p; q++,pq++) {
          ab_dens2[rs][pq] = A_AB[p][q];
          BB_norm += A_AB[p][q] * A_AB[p][q];
        }
      }

    }
  }

  fprintf(outfile, "\n\tBB_norm (2nd half) = %20.15f\n", BB_norm);
  BB_norm = 0.0;
  */

  /* compute to the BB energy to test above backtr */
  /*
  ntei = dst_ntri * (dst_ntri + 1)/2;
  ints = init_array(ntei);
  iwl_buf_init(&JBuff, PSIF_SO_TEI, tolerance, 1, 0);
  iwl_buf_rd_all(&JBuff, ints, ioff, ioff, 0, ioff, 0, outfile);
  iwl_buf_close(&JBuff, 1);

  energy = 0.0;
  for(p=0; p < dst_orbs; p++) {
    for(q=0; q < dst_orbs; q++) {
      pq = INDEX(p,q);
      for(r=0; r < dst_orbs; r++) {
        for(s=0; s < dst_orbs; s++) {
          rs = INDEX(r,s);

          pqrs = INDEX(pq,rs);

          energy += ab_dens2[pq][rs] * ints[pqrs];
        }
      }
    }
  }
  fprintf(outfile, "\n\tBB energy from AO-twopdm: %20.14f\n", energy);
  free(ints);

  free_block(ab_dens1);
  free_block(ab_dens2);
  */

  /** BB First-half transformation **/
  fprintf(outfile, "\n\tBeginning BB twopdm transform...\n");

  PBB_block = init_array(src_ntri);
  iwl_buf_init(&PBB_Buff, PSIF_BB_PRESORT, tolerance, 1, 1);

  J_block = init_array(MAX0(src_ntri,dst_ntri));
  yosh_init(&YBuffJ, dst_ntri, src_ntri, maxcor, maxcord, max_buckets,
            first_tmp_file, tolerance, outfile);
  if(print_lvl > 1) { yosh_print(&YBuffJ, outfile); fflush(outfile); }
  yosh_init_buckets(&YBuffJ);

  for(psym=0; psym < nirreps; psym++) {
    pfirst = src_first[psym];
    plast = src_last[psym];
    for (p=pfirst; p <= plast; p++) {
      for (qsym=0; qsym <= psym; qsym++) {
        qfirst = src_first[qsym];
        qlast = src_last[qsym];
        pqsym = psym^qsym;
        for (q=qfirst; (q<=qlast) && (q <= p); q++) {
          pq = ioff[p] + q;

          zero_arr(PBB_block, src_ntri);

          iwl_buf_rd(&PBB_Buff, pq, PBB_block, ioff, ioff, 0, 0, outfile);

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
                    A_BB[R][S] = PBB_block[rs];
                    A_BB[S][R] = PBB_block[rs];
                  }
                }
              }
              else {
                for (r=rfirst,R=0; r <= rlast; r++,R++) {
                  for (s=sfirst,S=0; s <= slast; s++,S++) {
                    rs = INDEX(r,s);
                    A_BB[R][S] = PBB_block[rs];
                  }
                }
              }

              /** BB half-transform for current pq **/
              if(C_colspi[ssym] > 0)
                C_DGEMM('n','t',src_orbspi[rsym],dst_orbspi[ssym],src_orbspi[ssym],1.0,
                        A_BB[0], A_cols,CB[ssym][0],C_colspi[ssym],0.0,B[0],B_cols);
              zero_mat(A_BB, A_cols, A_cols);
              if(C_colspi[rsym] > 0)
                C_DGEMM('n','n',dst_orbspi[rsym],dst_orbspi[ssym],src_orbspi[rsym],1.0,
                        CB[rsym][0],C_colspi[rsym],B[0],B_cols,0.0,A_BB[0],A_cols);

              /* collect the results J_block */
              zero_arr(J_block, dst_ntri);
              for(k=kfirst,K=0; k <= klast; k++,K++) {
                for (l=lfirst,L=0; (l <= llast) && (l <= k); l++,L++) {
                  kl = ioff[k] + l;
                  J_block[kl] = A_BB[K][L];
                }
              }

              yosh_wrt_arr(&YBuffJ, p, q, pq, pqsym, J_block,
                           moinfo.nao, ioff, dst_orbsym, dst_first, dst_last, 1, 0, outfile);

            }
          }
        }
      }
    }
  }

  iwl_buf_close(&PBB_Buff, 0);
  free(PBB_block);

  if (params.print_lvl) {
    fprintf(outfile, "\tSorting BB half-transformed twopdm...\n"); fflush(outfile);
  }
  yosh_flush(&YBuffJ);
  yosh_close_buckets(&YBuffJ, 0);
  yosh_sort(&YBuffJ, params.jfile+1, 0, ioff, NULL, src_orbs, src_ntri, 0, 1, 0, 0, 1, 0, outfile);
  yosh_done(&YBuffJ);
  free(J_block);
  if (print_lvl) {
    fprintf(outfile, "\tFinished BB half-transformation...\n"); fflush(outfile);
  }

  if (print_lvl) {
    fprintf(outfile, "\tStarting final half-transformation...\n"); fflush(outfile);
  }

  JA_block = init_array(MAX0(src_ntri,dst_ntri));
  JB_block = init_array(MAX0(src_ntri,dst_ntri));
  iwl_buf_init(&JA_Buff, params.jfile, tolerance, 1, 1);
  iwl_buf_init(&JB_Buff, params.jfile+1, tolerance, 1, 1);

  backsort_prep(1);
  twopdm_out = (struct iwlbuf *) malloc(nbuckets * sizeof(struct iwlbuf));
  for(p=0; p < nbuckets; p++)
    iwl_buf_init(&twopdm_out[p], first_tmp_file+p, tolerance, 0, 0);

  /*  iwl_buf_init(&MBuff, 78, tolerance, 0, 0); */

  /*
  reorder = init_int_array(dst_orbs);
  for(p=0; p < dst_orbs; p++) reorder[p] = p;
  */

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

          zero_arr(JA_block, dst_ntri);
          zero_arr(JB_block, dst_ntri);
          iwl_buf_rd(&JA_Buff, kl, JA_block, ioff, ioff, 0, 0, outfile);
          iwl_buf_rd(&JB_Buff, kl, JB_block, ioff, ioff, 0, 0, outfile);

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
                  A_AA[P][Q] = JA_block[pq];
                  A_BB[P][Q] = JB_block[pq];
                }
              }

              /** AA/AB second-half-transform for current pq **/
              if (C_colspi[qsym] > 0)
                C_DGEMM('n', 't', src_orbspi[psym],dst_orbspi[qsym], src_orbspi[qsym], 1.0,
                        A_AA[0], A_cols, CA[qsym][0], C_colspi[qsym], 0.0, B[0], B_cols);
              zero_mat(A_AA, A_cols, A_cols);
              if (C_colspi[psym] > 0)
                C_DGEMM('n', 'n', dst_orbspi[psym], dst_orbspi[qsym], src_orbspi[psym], 1.0,
                        CA[psym][0], C_colspi[psym], B[0], B_cols, 0.0, A_AA[0], A_cols);

              /** BB second-half-transform for current pq **/
              if (C_colspi[qsym] > 0)
                C_DGEMM('n', 't', src_orbspi[psym],dst_orbspi[qsym], src_orbspi[qsym], 1.0,
                        A_BB[0], A_cols, CB[qsym][0], C_colspi[qsym], 0.0, B[0], B_cols);
              zero_mat(A_BB, A_cols, A_cols);
              if (C_colspi[psym] > 0)
                C_DGEMM('n', 'n', dst_orbspi[psym], dst_orbspi[qsym], src_orbspi[psym], 1.0,
                        CB[psym][0], C_colspi[psym], B[0], B_cols, 0.0, A_BB[0], A_cols);

              /** combine AA/AB and BB transformed twopdm's for final sort **/
              for (i=ifirst,I=0; i <= ilast; i++,I++) {
                for (j=jfirst,J=0; j <= jlast; j++,J++) {
                  A_AA[I][J] += A_BB[I][J];
                }
              }

              /*
              iwl_buf_wrt_mat2(&MBuff, k, l, A_AA, ifirst, ilast, jfirst, jlast, reorder, 0, 0,
                               ioff, outfile);
              */

              backsort_write(k, l, A_AA, ifirst, ilast, jfirst, jlast, 0, outfile, twopdm_out, 1);

            }
          }
        }
      }
    }
  }

  free(JA_block);
  free(JB_block);
  iwl_buf_close(&JA_Buff, 0);
  iwl_buf_close(&JB_Buff, 0);

  /*
  iwl_buf_flush(&MBuff, 1);
  iwl_buf_close(&MBuff, 1);
  free(reorder);
  */

  free_block(A_AA);
  free_block(A_AB);
  free_block(A_BB);
  free_block(B);

  for(p=0; p < nbuckets; p++) {
    iwl_buf_flush(&twopdm_out[p], 1);
    iwl_buf_close(&twopdm_out[p], 1);
  }
  free(twopdm_out);
  if(print_lvl) {
    fprintf(outfile, "\n\tSorting AO-basis twopdm...");
    fflush(outfile);
  }
  backsort(first_tmp_file, tolerance, 1);
  if(print_lvl) {
    fprintf(outfile, "\n\tdone.");
    fflush(outfile);
  }

  if (print_lvl) {
    fprintf(outfile, "\n\tAA/AB/BB twopdm transformation finished.\n");
    fprintf(outfile, "\tAO-basis twopdm written to file%d.\n",params.mfile);
    fflush(outfile);
  }

  /*
  ntei = dst_ntri * (dst_ntri + 1)/2;
  ints = init_array(ntei);
  dens = block_matrix(dst_ntri, dst_ntri);

  iwl_buf_init(&MBuff, 78, tolerance, 1, 0);
  iwl_buf_rd_all2(&MBuff, dens, ioff, ioff, 0, ioff, 0, outfile);
  iwl_buf_init(&JBuff, PSIF_SO_TEI, tolerance, 1, 0);
  iwl_buf_rd_all(&JBuff, ints, ioff, ioff, 0, ioff, 0, outfile);
  energy = 0.0;
  for(p=0; p < dst_orbs; p++) {
    for(q=0; q < dst_orbs; q++) {
      pq = INDEX(p,q);
      for(r=0; r < dst_orbs; r++) {
        for(s=0; s < dst_orbs; s++) {
          rs = INDEX(r,s);
          pqrs = INDEX(pq,rs);
          energy += dens[pq][rs] * ints[pqrs];
        }
      }
    }
  }
  fprintf(outfile, "\n\tTotal energy from AO-twopdm: %20.14f\n", energy);
  iwl_buf_close(&MBuff, 1);
  iwl_buf_close(&JBuff, 1);

  free(ints);
  free_block(dens);
  */
}

}} // end namespace psi::transqt
