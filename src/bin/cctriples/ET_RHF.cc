/*! \file \ingroup CCTRIPLES
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <exception.h>
#include <libqt/qt.h>
#include <pthread.h>
#include <psiconfig.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

//MKL Header
#ifdef HAVE_MKL
#include <mkl.h>
#endif


namespace psi { namespace cctriples {

pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;

struct thread_data {
 dpdfile2 *fIJ; dpdfile2 *fAB; dpdfile2 *fIA; dpdfile2 *T1;
 dpdbuf4 *T2; dpdbuf4 *Eints; dpdbuf4 *Dints; dpdbuf4 *Fints_local;
 double *ET_local; int Gi; int Gj; int Gk; int first_ijk; int last_ijk;
};

void *ET_RHF_thread(void *thread_data);

double ET_RHF(void)
{
  int i,j,k,I,J,K,Gi,Gj,Gk, h, nirreps, cnt;
  int nijk, nthreads, thread, *ijk_part, errcod;
  int *occpi, *virtpi, *occ_off, *vir_off;
  double ET, *ET_array;
  dpdfile2 fIJ, fAB, fIA, T1;
  dpdbuf4 T2, Eints, Dints, *Fints_array;
  FILE *ijkfile;
  pthread_t  *p_thread;
  struct thread_data *thread_data_array;

  timer_on("ET_RHF");

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off;
  vir_off = moinfo.vir_off;

  nthreads = params.nthreads;
  thread_data_array = (struct thread_data *) malloc(nthreads*sizeof(struct thread_data));
  p_thread = (pthread_t *) malloc(nthreads*sizeof(pthread_t));

#ifdef HAVE_MKL
  int old_threads = mkl_get_max_threads();
  mkl_set_num_threads(1);
#endif

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_init(&fIA);
  dpd_file2_mat_rd(&fIJ);
  dpd_file2_mat_rd(&fAB);
  dpd_file2_mat_rd(&fIA);

  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_mat_init(&T1);
  dpd_file2_mat_rd(&T1);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2, h);
    dpd_buf4_mat_irrep_rd(&T2, h);

    dpd_buf4_mat_irrep_init(&Eints, h);
    dpd_buf4_mat_irrep_rd(&Eints, h);

    dpd_buf4_mat_irrep_init(&Dints, h);
    dpd_buf4_mat_irrep_rd(&Dints, h);
  }
  ffile(&ijkfile,"ijk.dat", 0);

  /* each thread gets its own F buffer to assign memory and read blocks
     into and its own energy double - all else shared */
  Fints_array = (dpdbuf4 *) malloc(nthreads*sizeof(dpdbuf4));
  for (thread=0; thread<nthreads;++thread)
    dpd_buf4_init(&(Fints_array[thread]), CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  ET_array = (double *) malloc(nthreads*sizeof(double));
  ijk_part = new int [nthreads];

  for (thread=0;thread<nthreads;++thread) {
    thread_data_array[thread].fIJ = &fIJ;
    thread_data_array[thread].fAB = &fAB;
    thread_data_array[thread].fIA = &fIA;
    thread_data_array[thread].T1 = &T1;
    thread_data_array[thread].T2 = &T2;
    thread_data_array[thread].Eints = &Eints;
    thread_data_array[thread].Dints = &Dints;
    thread_data_array[thread].Fints_local = &(Fints_array[thread]);
    thread_data_array[thread].ET_local = &(ET_array[thread]);
  }

  /* Compute total number of IJK combinations */
  nijk = 0;
  for(Gi=0; Gi < nirreps; Gi++)
    for(Gj=0; Gj < nirreps; Gj++)
      for(Gk=0; Gk < nirreps; Gk++)
        for(i=0; i < occpi[Gi]; i++) {
          I = occ_off[Gi] + i;
          for(j=0; j < occpi[Gj]; j++) {
            J = occ_off[Gj] + j;
            for(k=0; k < occpi[Gk]; k++) {
              K = occ_off[Gk] + k;
              if(I >= J && J >= K) nijk++;
            }
          }
        }
  fprintf(ijkfile, "Total number of IJK combinations =: %d\n", nijk);
  fflush(ijkfile);

  ET = 0.0;
  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      for(Gk=0; Gk < nirreps; Gk++) {

        nijk = 0; // total number of ijk for these irreps
        for(i=0; i < occpi[Gi]; i++) {
          I = occ_off[Gi] + i;
          for(j=0; j < occpi[Gj]; j++) {
            J = occ_off[Gj] + j;
            for(k=0; k < occpi[Gk]; k++) {
              K = occ_off[Gk] + k;
              if(I >= J && J >= K) nijk++;
            }
          }
        }
        fprintf(ijkfile, "Num. of IJK with (Gi,Gj,Gk)=(%d,%d,%d) =: %d\n",
          Gi, Gj, Gk, nijk);
        fflush(ijkfile);
        if (nijk == 0) continue;

        for (thread=0; thread<nthreads;++thread) {
          thread_data_array[thread].Gi = Gi;
          thread_data_array[thread].Gj = Gj;
          thread_data_array[thread].Gk = Gk;
          ET_array[thread] = 0.0;
          ijk_part[thread] = nijk / nthreads; // number of ijk for each thread
          if (thread < (nijk % nthreads)) ++ijk_part[thread];
        }

        cnt = 0;
        for (thread=0; thread<nthreads; ++thread) {
          if (!ijk_part[thread]) continue;  // there are more threads than nijk
          thread_data_array[thread].first_ijk = cnt;
          cnt += ijk_part[thread];
          thread_data_array[thread].last_ijk = cnt-1;
        }

        /* execute threads */
        for (thread=0; thread<nthreads;++thread) {
          if (!ijk_part[thread]) continue;
          fprintf(ijkfile,"\tthread %d: first_ijk=%d,  last_ijk=%d\n", thread,
            thread_data_array[thread].first_ijk, thread_data_array[thread].last_ijk);
          fflush(ijkfile);
        }

        for (thread=0;thread<nthreads;++thread) {
          if (!ijk_part[thread]) continue;
          errcod = pthread_create(&(p_thread[thread]), NULL, ET_RHF_thread,
                   (void *) &thread_data_array[thread]);
          if (errcod) {
            throw PsiException("pthread_create in ET_RHF()",__FILE__,__LINE__);
          }
        }

        for (thread=0; thread<nthreads;++thread) {
          if (!ijk_part[thread]) continue;
          errcod = pthread_join(p_thread[thread], NULL);
          if (errcod) {
            throw PsiException("pthread_join in ET_RHF() failed",__FILE__,__LINE__);
          }
        }

        for (thread=0;thread<nthreads;++thread)
          ET += ET_array[thread];

      } /* Gk */
    } /* Gj */
  } /* Gi */

  fclose(ijkfile);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&T2, h);
    dpd_buf4_mat_irrep_close(&Eints, h);
    dpd_buf4_mat_irrep_close(&Dints, h);
  }
  dpd_buf4_close(&T2);
  dpd_buf4_close(&Eints);
  dpd_buf4_close(&Dints);

  dpd_file2_mat_close(&T1);
  dpd_file2_close(&T1);

  dpd_file2_mat_close(&fIJ);
  dpd_file2_mat_close(&fAB);
  dpd_file2_mat_close(&fIA);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fAB);
  dpd_file2_close(&fIA);

  for (thread=0; thread<nthreads; ++thread)
    dpd_buf4_close(&(Fints_array[thread]));

  free(Fints_array);
  free(ET_array);
  delete [] ijk_part;

  free(thread_data_array);
  free(p_thread);

  timer_off("ET_RHF");

#ifdef HAVE_MKL
  mkl_set_num_threads(old_threads);
#endif

  return ET;
}


void* ET_RHF_thread(void* thread_data_in)
{
  int h, nirreps, cnt_ijk;
  int Gp, p, nump;
  int nrows, ncols, nlinks;
  int Gijk, Gid, Gkd, Gjd, Gil, Gkl, Gjl;
  int Gab, Gba, Gbc, Gcb, Gac, Gca;
  int ab, ba, bc, cb, ac, ca;
  int cd, bd, ad, lc, lb, la;
  int il, jl, kl;
  int Gi, Gj, Gk, Ga, Gb, Gc, Gd, Gl;
  int Gij, Gji, Gjk, Gkj, Gik, Gki;
  int I, J, K, A, B, C, D, L;
  int i, j, k, a, b, c, d, l;
  int ij, ji, ik, ki, jk, kj;
  int *occpi, *virtpi, *occ_off, *vir_off;
  double t_ia, t_jb, t_kc, D_jkbc, D_ikac, D_ijab;
  double f_ia, f_jb, f_kc, t_jkbc, t_ikac, t_ijab;
  double dijk, value1, value2, value3, value4, value5, value6, denom, *ET_local;
  double ***W0, ***W1, ***V, ***X, ***Y, ***Z;
  dpdbuf4 *T2, *Eints, *Dints, *Fints;
  dpdfile2 *fIJ, *fAB, *fIA, *T1;
  int nijk, nthreads, first_ijk, last_ijk, thr_id;
  struct thread_data *data;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off;
  vir_off = moinfo.vir_off;

  data = (struct thread_data *) thread_data_in;
  fIJ   = data->fIJ;
  fAB   = data->fAB;
  fIA   = data->fIA;
  T1    = data->T1;
  T2    = data->T2;
  Eints = data->Eints;
  Dints = data->Dints;
  Fints = data->Fints_local;
  ET_local  = data->ET_local; // pointer to where thread E goes
  Gi        =  data->Gi;
  Gj        =  data->Gj;
  Gk        =  data->Gk;
  first_ijk = data->first_ijk;
  last_ijk  = data->last_ijk;

  W0 = (double ***) malloc(nirreps * sizeof(double **));
  W1 = (double ***) malloc(nirreps * sizeof(double **));
  V = (double ***) malloc(nirreps * sizeof(double **));
  X = (double ***) malloc(nirreps * sizeof(double **));
  Y = (double ***) malloc(nirreps * sizeof(double **));
  Z = (double ***) malloc(nirreps * sizeof(double **));

  Gkj = Gjk = Gk ^ Gj;
  Gji = Gij = Gi ^ Gj;
  Gik = Gki = Gi ^ Gk;
  Gijk = Gi ^ Gj ^ Gk;

  cnt_ijk = -1;
  for(i=0; i < occpi[Gi]; i++) {
    I = occ_off[Gi] + i;
    for(j=0; j < occpi[Gj]; j++) {
      J = occ_off[Gj] + j;
      for(k=0; k < occpi[Gk]; k++) {
        K = occ_off[Gk] + k;
        if(I >= J && J >= K) {

          ++cnt_ijk;
          /* check to see if this ijk is for this thread */
          if ( (cnt_ijk < first_ijk) || (cnt_ijk > last_ijk))
             continue;

          ij = T2->params->rowidx[I][J];
          ji = T2->params->rowidx[J][I];
          ik = T2->params->rowidx[I][K];
          ki = T2->params->rowidx[K][I];
          jk = T2->params->rowidx[J][K];
          kj = T2->params->rowidx[K][J];

          dijk = 0.0;
          if(fIJ->params->rowtot[Gi])
            dijk += fIJ->matrix[Gi][i][i];
          if(fIJ->params->rowtot[Gj])
            dijk += fIJ->matrix[Gj][j][j];
          if(fIJ->params->rowtot[Gk])
            dijk += fIJ->matrix[Gk][k][k];

                /* Malloc space for the W intermediate */

                // timer_on("malloc");
                for(Gab=0; Gab < nirreps; Gab++) {
                  Gc = Gab ^ Gijk;

                  W0[Gab] = dpd_block_matrix(Fints->params->coltot[Gab],virtpi[Gc]);
                  W1[Gab] = dpd_block_matrix(Fints->params->coltot[Gab],virtpi[Gc]);
                }
                // timer_off("malloc");

                // timer_on("N7 Terms");

                /* +F_idab * t_kjcd */
                for(Gd=0; Gd < nirreps; Gd++) {

                  Gab = Gid = Gi ^ Gd;
                  Gc = Gkj ^ Gd;

                  /* Set up F integrals */
                  Fints->matrix[Gid] = dpd_block_matrix(virtpi[Gd], Fints->params->coltot[Gid]);
          pthread_mutex_lock(&mut);
                  dpd_buf4_mat_irrep_rd_block(Fints, Gid, Fints->row_offset[Gid][I], virtpi[Gd]);
          pthread_mutex_unlock(&mut);

                  /* Set up T2 amplitudes */
                  cd = T2->col_offset[Gkj][Gc];

                  /* Set up multiplication parameters */
                  nrows = Fints->params->coltot[Gid];
                  ncols = virtpi[Gc];
                  nlinks = virtpi[Gd];

                  if(nrows && ncols && nlinks)
                    C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                            &(Fints->matrix[Gid][0][0]), nrows,
                            &(T2->matrix[Gkj][kj][cd]), nlinks, 0.0,
                            &(W0[Gab][0][0]), ncols);

                  dpd_free_block(Fints->matrix[Gid], virtpi[Gd], Fints->params->coltot[Gid]);
                }

                /* -E_jklc * t_ilab */
                for(Gl=0; Gl < nirreps; Gl++) {

                  Gab = Gil = Gi ^ Gl;
                  Gc = Gjk ^ Gl;

                  /* Set up E integrals */
                  lc = Eints->col_offset[Gjk][Gl];

                  /* Set up T2 amplitudes */
                  il = T2->row_offset[Gil][I];

                  /* Set up multiplication parameters */
                  nrows = T2->params->coltot[Gil];
                  ncols = virtpi[Gc];
                  nlinks = occpi[Gl];

                  if(nrows && ncols && nlinks)
                    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                            &(T2->matrix[Gil][il][0]), nrows,
                            &(Eints->matrix[Gjk][jk][lc]), ncols, 1.0,
                            &(W0[Gab][0][0]), ncols);
                }

                /* Sort W[ab][c] --> W[ac][b] */
                dpd_3d_sort(W0, W1, nirreps, Gijk, Fints->params->coltot, Fints->params->colidx,
                       Fints->params->colorb, Fints->params->rsym, Fints->params->ssym,
                       vir_off, vir_off, virtpi, vir_off, Fints->params->colidx, acb, 0);

                /* +F_idac * t_jkbd */
                for(Gd=0; Gd < nirreps; Gd++) {

                  Gac = Gid = Gi ^ Gd;
                  Gb = Gjk ^ Gd;

                  Fints->matrix[Gid] = dpd_block_matrix(virtpi[Gd], Fints->params->coltot[Gid]);
          pthread_mutex_lock(&mut);
                  dpd_buf4_mat_irrep_rd_block(Fints, Gid, Fints->row_offset[Gid][I], virtpi[Gd]);
          pthread_mutex_unlock(&mut);

                  bd = T2->col_offset[Gjk][Gb];

                  nrows = Fints->params->coltot[Gid];
                  ncols = virtpi[Gb];
                  nlinks = virtpi[Gd];

                  if(nrows && ncols && nlinks)
                    C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                            &(Fints->matrix[Gid][0][0]), nrows,
                            &(T2->matrix[Gjk][jk][bd]), nlinks, 1.0,
                            &(W1[Gac][0][0]), ncols);

                  dpd_free_block(Fints->matrix[Gid], virtpi[Gd], Fints->params->coltot[Gid]);
                }

                /* -E_kjlb * t_ilac */
                for(Gl=0; Gl < nirreps; Gl++) {

                  Gac = Gil = Gi ^ Gl;
                  Gb = Gkj ^ Gl;

                  lb = Eints->col_offset[Gkj][Gl];

                  il = T2->row_offset[Gil][I];

                  nrows = T2->params->coltot[Gil];
                  ncols = virtpi[Gb];
                  nlinks = occpi[Gl];

                  if(nrows && ncols && nlinks)
                    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                            &(T2->matrix[Gil][il][0]), nrows,
                            &(Eints->matrix[Gkj][kj][lb]), ncols, 1.0,
                            &(W1[Gac][0][0]), ncols);
                }

                /* Sort W[ac][b] --> W[ca][b] */
                dpd_3d_sort(W1, W0, nirreps, Gijk, Fints->params->coltot, Fints->params->colidx,
                       Fints->params->colorb, Fints->params->rsym, Fints->params->ssym,
                       vir_off, vir_off, virtpi, vir_off, Fints->params->colidx, bac, 0);

                /* +F_kdca * t_jibd */
                for(Gd=0; Gd < nirreps; Gd++) {

                  Gca = Gkd = Gk ^ Gd;
                  Gb = Gji ^ Gd;

                  Fints->matrix[Gkd] = dpd_block_matrix(virtpi[Gd], Fints->params->coltot[Gkd]);
          pthread_mutex_lock(&mut);
                  dpd_buf4_mat_irrep_rd_block(Fints, Gkd, Fints->row_offset[Gkd][K], virtpi[Gd]);
          pthread_mutex_unlock(&mut);

                  bd = T2->col_offset[Gji][Gb];

                  nrows = Fints->params->coltot[Gkd];
                  ncols = virtpi[Gb];
                  nlinks = virtpi[Gd];

                  if(nrows && ncols && nlinks)
                    C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                            &(Fints->matrix[Gkd][0][0]), nrows,
                            &(T2->matrix[Gji][ji][bd]), nlinks, 1.0,
                            &(W0[Gca][0][0]), ncols);

                  dpd_free_block(Fints->matrix[Gkd], virtpi[Gd], Fints->params->coltot[Gkd]);
                }

                /* -E_ijlb * t_klca */
                for(Gl=0; Gl < nirreps; Gl++) {

                  Gca = Gkl = Gk ^ Gl;
                  Gb = Gij ^ Gl;

                  lb = Eints->col_offset[Gij][Gl];

                  kl = T2->row_offset[Gkl][K];

                  nrows = T2->params->coltot[Gkl];
                  ncols = virtpi[Gb];
                  nlinks = occpi[Gl];

                  if(nrows && ncols && nlinks)
                    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                            &(T2->matrix[Gkl][kl][0]), nrows,
                            &(Eints->matrix[Gij][ij][lb]), ncols, 1.0,
                            &(W0[Gca][0][0]), ncols);
                }

                /* Sort W[ca][b] --> W[cb][a] */
                dpd_3d_sort(W0, W1, nirreps, Gijk, Fints->params->coltot, Fints->params->colidx,
                       Fints->params->colorb, Fints->params->rsym, Fints->params->ssym,
                       vir_off, vir_off, virtpi, vir_off, Fints->params->colidx, acb, 0);

                /* +F_kdcb * t_ijad */
                for(Gd=0; Gd < nirreps; Gd++) {

                  Gcb = Gkd = Gk ^ Gd;
                  Ga = Gij ^ Gd;

                  Fints->matrix[Gkd] = dpd_block_matrix(virtpi[Gd], Fints->params->coltot[Gkd]);
          pthread_mutex_lock(&mut);
                  dpd_buf4_mat_irrep_rd_block(Fints, Gkd, Fints->row_offset[Gkd][K], virtpi[Gd]);
          pthread_mutex_unlock(&mut);

                  ad = T2->col_offset[Gij][Ga];

                  nrows = Fints->params->coltot[Gkd];
                  ncols = virtpi[Ga];
                  nlinks = virtpi[Gd];

                  if(nrows && ncols && nlinks)
                    C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                            &(Fints->matrix[Gkd][0][0]), nrows,
                            &(T2->matrix[Gij][ij][ad]), nlinks, 1.0,
                            &(W1[Gcb][0][0]), ncols);

                  dpd_free_block(Fints->matrix[Gkd], virtpi[Gd], Fints->params->coltot[Gkd]);
                }

                /* -E_jila * t_klcb */
                for(Gl=0; Gl < nirreps; Gl++) {

                  Gcb = Gkl = Gk ^ Gl;
                  Ga = Gji ^ Gl;

                  la = Eints->col_offset[Gji][Gl];

                  kl = T2->row_offset[Gkl][K];

                  nrows = T2->params->coltot[Gkl];
                  ncols = virtpi[Ga];
                  nlinks = occpi[Gl];

                  if(nrows && ncols && nlinks)
                    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                            &(T2->matrix[Gkl][kl][0]), nrows,
                            &(Eints->matrix[Gji][ji][la]), ncols, 1.0,
                            &(W1[Gcb][0][0]), ncols);
                }

                /* Sort W[cb][a] --> W[bc][a] */
                dpd_3d_sort(W1, W0, nirreps, Gijk, Fints->params->coltot, Fints->params->colidx,
                       Fints->params->colorb, Fints->params->rsym, Fints->params->ssym,
                       vir_off, vir_off, virtpi, vir_off, Fints->params->colidx, bac, 0);

                /* +F_jdbc * t_ikad */
                for(Gd=0; Gd < nirreps; Gd++) {

                  Gbc = Gjd = Gj ^ Gd;
                  Ga = Gik ^ Gd;

                  Fints->matrix[Gjd] = dpd_block_matrix(virtpi[Gd], Fints->params->coltot[Gjd]);
          pthread_mutex_lock(&mut);
                  dpd_buf4_mat_irrep_rd_block(Fints, Gjd, Fints->row_offset[Gjd][J], virtpi[Gd]);
          pthread_mutex_unlock(&mut);

                  ad = T2->col_offset[Gik][Ga];

                  nrows = Fints->params->coltot[Gjd];
                  ncols = virtpi[Ga];
                  nlinks = virtpi[Gd];

                  if(nrows && ncols && nlinks)
                    C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                            &(Fints->matrix[Gjd][0][0]), nrows,
                            &(T2->matrix[Gik][ik][ad]), nlinks, 1.0,
                            &(W0[Gbc][0][0]), ncols);

                  dpd_free_block(Fints->matrix[Gjd], virtpi[Gd], Fints->params->coltot[Gjd]);
                }

                /* -E_kila * t_jlbc */
                for(Gl=0; Gl < nirreps; Gl++) {

                  Gbc = Gjl = Gj ^ Gl;
                  Ga = Gki ^ Gl;

                  la = Eints->col_offset[Gki][Gl];

                  jl = T2->row_offset[Gjl][J];

                  nrows = T2->params->coltot[Gjl];
                  ncols = virtpi[Ga];
                  nlinks = occpi[Gl];

                  if(nrows && ncols && nlinks)
                    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                            &(T2->matrix[Gjl][jl][0]), nrows,
                            &(Eints->matrix[Gki][ki][la]), ncols, 1.0,
                            &(W0[Gbc][0][0]), ncols);
                }

                /* Sort W[bc][a] --> W[ba][c] */
                dpd_3d_sort(W0, W1, nirreps, Gijk, Fints->params->coltot, Fints->params->colidx,
                       Fints->params->colorb, Fints->params->rsym, Fints->params->ssym,
                       vir_off, vir_off, virtpi, vir_off, Fints->params->colidx, acb, 0);

                /* +F_jdba * t_kicd */
                for(Gd=0; Gd < nirreps; Gd++) {

                  Gba = Gjd = Gj ^ Gd;
                  Gc = Gki ^ Gd;

                  Fints->matrix[Gjd] = dpd_block_matrix(virtpi[Gd], Fints->params->coltot[Gjd]);
          pthread_mutex_lock(&mut);
                  dpd_buf4_mat_irrep_rd_block(Fints, Gjd, Fints->row_offset[Gjd][J], virtpi[Gd]);
          pthread_mutex_unlock(&mut);

                  cd = T2->col_offset[Gki][Gc];

                  nrows = Fints->params->coltot[Gjd];
                  ncols = virtpi[Gc];
                  nlinks = virtpi[Gd];

                  if(nrows && ncols && nlinks)
                    C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
                            &(Fints->matrix[Gjd][0][0]), nrows,
                            &(T2->matrix[Gki][ki][cd]), nlinks, 1.0,
                            &(W1[Gba][0][0]), ncols);

                  dpd_free_block(Fints->matrix[Gjd], virtpi[Gd], Fints->params->coltot[Gjd]);
                }

                /* -E_iklc * t_jlba */
                for(Gl=0; Gl < nirreps; Gl++) {

                  Gba = Gjl = Gj ^ Gl;
                  Gc = Gik ^ Gl;

                  lc = Eints->col_offset[Gik][Gl];

                  jl = T2->row_offset[Gjl][J];

                  nrows = T2->params->coltot[Gjl];
                  ncols = virtpi[Gc];
                  nlinks = occpi[Gl];

                  if(nrows && ncols && nlinks)
                    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
                            &(T2->matrix[Gjl][jl][0]), nrows,
                            &(Eints->matrix[Gik][ik][lc]), ncols, 1.0,
                            &(W1[Gba][0][0]), ncols);
                }

                /* Sort W[ba][c] --> W[ab][c] */
                dpd_3d_sort(W1, W0, nirreps, Gijk, Fints->params->coltot, Fints->params->colidx,
                       Fints->params->colorb, Fints->params->rsym, Fints->params->ssym,
                       vir_off, vir_off, virtpi, vir_off, Fints->params->colidx, bac, 0);

                // timer_off("N7 Terms");

                // timer_on("malloc");
                for(Gab=0; Gab < nirreps; Gab++) {
                  Gc = Gab ^ Gijk;
                  dpd_free_block(W1[Gab],Fints->params->coltot[Gab],virtpi[Gc]);

                  V[Gab] = dpd_block_matrix(Fints->params->coltot[Gab],virtpi[Gc]);
                }
                // timer_off("malloc");

                /* Copy W intermediate into V */
                for(Gab=0; Gab < nirreps; Gab++) {
                  Gc = Gab ^ Gijk;

                  for(ab=0; ab < Fints->params->coltot[Gab]; ab++) {
                    for(c=0; c < virtpi[Gc]; c++) {

                      V[Gab][ab][c] = W0[Gab][ab][c];
                    }
                  }
                }

                // timer_on("EST Terms");

                /* Add EST terms to V */

                for(Gab=0; Gab < nirreps; Gab++) {

                  Gc = Gab ^ Gijk;

                  for(ab=0; ab < Fints->params->coltot[Gab]; ab++) {

                    A = Fints->params->colorb[Gab][ab][0];
                    Ga = Fints->params->rsym[A];
                    a = A - vir_off[Ga];
                    B = Fints->params->colorb[Gab][ab][1];
                    Gb = Fints->params->ssym[B];
                    b = B - vir_off[Gb];

                    Gbc = Gb ^ Gc;
                    Gac = Ga ^ Gc;

                    for(c=0; c < virtpi[Gc]; c++) {
                      C = vir_off[Gc] + c;

                      bc = Dints->params->colidx[B][C];
                      ac = Dints->params->colidx[A][C];

                      /* +t_ia * D_jkbc + f_ia * t_jkbc */
                      if(Gi == Ga && Gjk == Gbc) {
                        t_ia = D_jkbc = 0.0;

                        if(T1->params->rowtot[Gi] && T1->params->coltot[Gi]) {
                          t_ia = T1->matrix[Gi][i][a];
                          f_ia = fIA->matrix[Gi][i][a];
                        }

                        if(Dints->params->rowtot[Gjk] && Dints->params->coltot[Gjk]) {
                          D_jkbc = Dints->matrix[Gjk][jk][bc];
                          t_jkbc = T2->matrix[Gjk][jk][bc];
                        }

                        V[Gab][ab][c] += t_ia * D_jkbc + f_ia * t_jkbc;

                      }

                      /* +t_jb * D_ikac */
                      if(Gj == Gb && Gik == Gac) {
                        t_jb = D_ikac = 0.0;

                        if(T1->params->rowtot[Gj] && T1->params->coltot[Gj]) {
                          t_jb = T1->matrix[Gj][j][b];
                          f_jb = fIA->matrix[Gj][j][b];
                        }

                        if(Dints->params->rowtot[Gik] && Dints->params->coltot[Gik]) {
                          D_ikac = Dints->matrix[Gik][ik][ac];
                          t_ikac = T2->matrix[Gik][ik][ac];
                        }

                        V[Gab][ab][c] += t_jb * D_ikac + f_jb * t_ikac;
                      }

                      /* +t_kc * D_ijab */
                      if(Gk == Gc && Gij == Gab) {
                        t_kc = D_ijab = 0.0;

                        if(T1->params->rowtot[Gk] && T1->params->coltot[Gk]) {
                          t_kc = T1->matrix[Gk][k][c];
                          f_kc = fIA->matrix[Gk][k][c];
                        }

                        if(Dints->params->rowtot[Gij] && Dints->params->coltot[Gij]) {
                          D_ijab = Dints->matrix[Gij][ij][ab];
                          t_ijab = T2->matrix[Gij][ij][ab];
                        }

                        V[Gab][ab][c] += t_kc * D_ijab + f_kc * t_ijab;
                      }

                      V[Gab][ab][c] /= (1 + (A==B) + (B==C) + (A==C));
                    }
                  }
                }

                // timer_off("EST Terms");

                // timer_on("malloc");
                for(Gab=0; Gab < nirreps; Gab++) {
                  Gc = Gab ^ Gijk;

                  X[Gab] = dpd_block_matrix(Fints->params->coltot[Gab],virtpi[Gc]);
                  Y[Gab] = dpd_block_matrix(Fints->params->coltot[Gab],virtpi[Gc]);
                  Z[Gab] = dpd_block_matrix(Fints->params->coltot[Gab],virtpi[Gc]);
                }
                // timer_off("malloc");

                // timer_on("XYZ");
                /* Build X, Y, and Z intermediates */

                for(Gab=0; Gab < nirreps; Gab++) {

                  Gc = Gab ^ Gijk;

                  Gba = Gab;

                  for(ab=0; ab < Fints->params->coltot[Gab]; ab++) {

                    A = Fints->params->colorb[Gab][ab][0];
                    Ga = Fints->params->rsym[A];
                    a = A - vir_off[Ga];
                    B = Fints->params->colorb[Gab][ab][1];
                    Gb = Fints->params->ssym[B];
                    b = B - vir_off[Gb];

                    Gac = Gca = Ga ^ Gc;
                    Gbc = Gcb = Gb ^ Gc;

                    ba = Dints->params->colidx[B][A];

                    for(c=0; c < virtpi[Gc]; c++) {
                      C = vir_off[Gc] + c;

                      ac = Dints->params->colidx[A][C];
                      ca = Dints->params->colidx[C][A];
                      bc = Dints->params->colidx[B][C];
                      cb = Dints->params->colidx[C][B];

                      X[Gab][ab][c] =
                        W0[Gab][ab][c] * V[Gab][ab][c] + W0[Gac][ac][b] * V[Gac][ac][b] +
                        W0[Gba][ba][c] * V[Gba][ba][c] + W0[Gbc][bc][a] * V[Gbc][bc][a] +
                        W0[Gca][ca][b] * V[Gca][ca][b] + W0[Gcb][cb][a] * V[Gcb][cb][a];

                      Y[Gab][ab][c] = V[Gab][ab][c] + V[Gbc][bc][a] + V[Gca][ca][b];

                      Z[Gab][ab][c] = V[Gac][ac][b] + V[Gba][ba][c] + V[Gcb][cb][a];

                    }
                  }
                }
                // timer_off("XYZ");

                // timer_on("malloc");
                for(Gab=0; Gab < nirreps; Gab++) {
                  Gc = Gab ^ Gijk;

                  dpd_free_block(V[Gab], Fints->params->coltot[Gab], virtpi[Gc]);
                }
                // timer_off("malloc");

                // timer_on("Energy");
                for(Gab=0; Gab < nirreps; Gab++) {

                  Gc = Gab ^ Gijk;
                  Gba = Gab;

                  for(ab=0; ab < Fints->params->coltot[Gab]; ab++) {

                    A = Fints->params->colorb[Gab][ab][0];
                    Ga = Fints->params->rsym[A];
                    a = A - vir_off[Ga];
                    B = Fints->params->colorb[Gab][ab][1];
                    Gb = Fints->params->ssym[B];
                    b = B - vir_off[Gb];

                    if(A >= B) {

                      Gac = Gca = Ga ^ Gc;
                      Gbc = Gcb = Gb ^ Gc;

                      ba = Dints->params->colidx[B][A];

                      for(c=0; c < virtpi[Gc]; c++) {
                        C = vir_off[Gc] + c;

                        if(B >= C) {

                          ac = Dints->params->colidx[A][C];
                          ca = Dints->params->colidx[C][A];
                          bc = Dints->params->colidx[B][C];
                          cb = Dints->params->colidx[C][B];

                          value1 = Y[Gab][ab][c] - 2.0 * Z[Gab][ab][c];
                          value2 = Z[Gab][ab][c] - 2.0 * Y[Gab][ab][c];
                          value3 = W0[Gab][ab][c] + W0[Gbc][bc][a] + W0[Gca][ca][b];
                          value4 = W0[Gac][ac][b] + W0[Gba][ba][c] + W0[Gcb][cb][a];
                          value5 = 3.0 * X[Gab][ab][c];
                          value6 = 2 - ((I==J) + (J==K) + (I==K));

                          denom = dijk;
                          if(fAB->params->rowtot[Ga])
                            denom -= fAB->matrix[Ga][a][a];
                          if(fAB->params->rowtot[Gb])
                            denom -= fAB->matrix[Gb][b][b];
                          if(fAB->params->rowtot[Gc])
                            denom -= fAB->matrix[Gc][c][c];

                          *ET_local += (value1 * value3 + value2 * value4 + value5) * value6/denom;

                        }

                      }

                    }
                  }
                }
                // timer_off("Energy");

                /* Free the W and V intermediates */
                // timer_on("malloc");
                for(Gab=0; Gab < nirreps; Gab++) {
                  Gc = Gab ^ Gijk;

                  dpd_free_block(W0[Gab],Fints->params->coltot[Gab],virtpi[Gc]);
                  dpd_free_block(X[Gab],Fints->params->coltot[Gab],virtpi[Gc]);
                  dpd_free_block(Y[Gab],Fints->params->coltot[Gab],virtpi[Gc]);
                  dpd_free_block(Z[Gab],Fints->params->coltot[Gab],virtpi[Gc]);
                }
                // timer_off("malloc");

              }
            } /* k */
          } /* j */
        } /* i */

  pthread_exit(NULL);
}

}} // namespace psi::CCTRIPLES

