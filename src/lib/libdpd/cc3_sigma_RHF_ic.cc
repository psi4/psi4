/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <ccfiles.h>
#include <pthread.h>
#include <psiconfig.h>

//MKL Header
#ifdef HAVE_MKL
#include <mkl.h>
#endif


namespace psi {

/*
  This function computes contributions to singles and doubles of
  matrix elements of triples:
    SIA   <-- <S|(Dints)           <T|(Wmbij,Wabei) CMNEF |0> |T> / (w-wt)
    SIjAb <-- <D|(FME,WmAEf,WMnIe) <T|(Wmbij,Wabei) CMNEF |0> |T> / (w-wt)
  Irrep variables are:
    GS <--            GW                  GWX3    ^ GC
                                  (           GX3            )
  These are used to make X3 quantity in T3_RHF:
    CIjAb, WAbEi, WMbIj, fIJ, fAB, omega
*/

struct thread_data {
 dpdbuf4 *CIjAb; dpdbuf4 *WAbEi; dpdbuf4 *WMbIj; int do_singles; dpdbuf4 *Dints;
 dpdfile2 *SIA; int do_doubles; dpdfile2 *FME; dpdbuf4 *WmAEf; dpdbuf4 *WMnIe;
 dpdbuf4 *SIjAb; int *occpi; int *occ_off; int *virtpi; int *vir_off;
 double omega; dpdfile2 *fIJ; dpdfile2 *fAB; int Gi; int Gj; int Gk;
 int first_ijk; int last_ijk;
 FILE *outfile; int thr_id; dpdfile2 SIA_local; dpdbuf4 SIjAb_local;
 int newtrips;
};

void *cc3_sigma_RHF_ic_thread(void *thread_data);

void cc3_sigma_RHF_ic(dpdbuf4 *CIjAb, dpdbuf4 *WAbEi, dpdbuf4 *WMbIj,
    int do_singles, dpdbuf4 *Dints, dpdfile2 *SIA,
    int do_doubles, dpdfile2 *FME, dpdbuf4 *WmAEf, dpdbuf4 *WMnIe,
    dpdbuf4 *SIjAb, int *occpi, int *occ_off, int *virtpi, int *vir_off,
    double omega, FILE *outfile, int nthreads, int newtrips)
{
  int h, nirreps, thread, nijk, *ijk_part, errcod=0;
  int Gi, Gj, Gk, Gl, Ga, Gb, Gc, Gd;
  int i, j, k, l, a, b, c, d, row, col;
  int I, J, K, L, A, B, C, D;
  int kj, jk, ji, ij, ik, ki;
  int Gkj, Gjk, Gji, Gij, Gik, Gki, Gkd;
  int Gijk, GS, GC, GWX3, GW, GX3, nrows,ncols,nlinks;
  int ab, ba, ac, ca, bc, cb;
  int Gab, Gba, Gac, Gca, Gbc, Gcb, Gid, Gjd;
  int id, jd, kd, ad, bd, cd;
  int il, jl, kl, la, lb, lc, li, lk;
  int da, di, dj, dk, thr_id;
  int Gad, Gdi, Gdj, Gdk, Glc, Gli, Glk,cnt,cnt2;
  long int length;
  double value, F_val, t_val, E_val;
  double dijk, denom, *tvect, **Z;
  double value_ia, value_ka, denom_ia, denom_ka;
  dpdfile2 fIJ, fAB, *SIA_local;
  dpdbuf4 buf4_tmp, *SIjAb_local;
  pthread_t  *p_thread;
  struct thread_data *thread_data_array;
  char lbl[32];

  thread_data_array = (struct thread_data *) malloc(nthreads*sizeof(struct thread_data));
  p_thread = (pthread_t *) malloc(nthreads*sizeof(pthread_t));

#ifdef HAVE_MKL
  int old_threads = mkl_get_max_threads();
  mkl_set_num_threads(1);
#endif

  nirreps = CIjAb->params->nirreps;
  /* these are sent to T3 function */
  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");

  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_rd(&fIJ);
  dpd_file2_mat_rd(&fAB);
  dpd_file2_mat_init(FME);
  dpd_file2_mat_rd(FME);

  GC = CIjAb->file.my_irrep;
  GWX3 = WAbEi->file.my_irrep;
  GX3 = GWX3^GC;
  GW = WmAEf->file.my_irrep;
  GS = SIjAb->file.my_irrep;
  if (GS != (GX3^GW)) {
    fprintf(outfile,"problem with irreps in cc3_sigma_RHF()\n"); 
    exit(1);
  }

  if (do_singles) {
    dpd_file2_mat_init(SIA);
    dpd_file2_mat_rd(SIA);
  }

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(CIjAb, h);
    dpd_buf4_mat_irrep_rd(CIjAb, h);
    dpd_buf4_mat_irrep_init(WMbIj, h);
    dpd_buf4_mat_irrep_rd(WMbIj, h);
    dpd_buf4_mat_irrep_init(WAbEi, h);
    dpd_buf4_mat_irrep_rd(WAbEi, h);
    dpd_buf4_mat_irrep_init(WmAEf, h);
    dpd_buf4_mat_irrep_rd(WmAEf, h);

    if (do_singles) {
      dpd_buf4_mat_irrep_init(Dints, h);
      dpd_buf4_mat_irrep_rd(Dints, h);
    }
    if (do_doubles) {
      dpd_buf4_mat_irrep_init(WMnIe, h);
      dpd_buf4_mat_irrep_rd(WMnIe, h);
      dpd_buf4_mat_irrep_init(SIjAb, h);
      dpd_buf4_mat_irrep_rd(SIjAb, h);
    }
  }

  SIA_local = (dpdfile2 *) malloc(nthreads*sizeof(dpdfile2));
  SIjAb_local = (dpdbuf4 *) malloc(nthreads*sizeof(dpdbuf4));

  for (i=0;i<nthreads;++i) {
    if (do_singles) {
      sprintf(lbl, "%s %d", "CC3 SIA", i);
      dpd_file2_init(&(SIA_local[i]), CC_TMP1, GS, 0, 1, lbl);
      dpd_file2_mat_init(&(SIA_local[i]));
    }
    if (do_doubles) {
      sprintf(lbl, "%s %d", "CC3 SIjAb", i);
      dpd_buf4_init(&(SIjAb_local[i]), CC_TMP1, GS, 0, 5, 0, 5, 0, lbl);
      for (h=0;h<nirreps;++h) {
        dpd_buf4_mat_irrep_init(&(SIjAb_local[i]),h);
      }
    }
  }

  for (thread=0;thread<nthreads;++thread) {
    thread_data_array[thread].CIjAb = CIjAb;
    thread_data_array[thread].WAbEi = WAbEi;
    thread_data_array[thread].WMbIj = WMbIj;
    thread_data_array[thread].do_singles = do_singles;
    thread_data_array[thread].Dints = Dints;
    thread_data_array[thread].SIA= SIA;
    thread_data_array[thread].do_doubles = do_doubles;
    thread_data_array[thread].FME = FME;
    thread_data_array[thread].WmAEf = WmAEf;
    thread_data_array[thread].WMnIe = WMnIe;
    thread_data_array[thread].SIjAb = SIjAb;
    thread_data_array[thread].occpi = occpi;
    thread_data_array[thread].occ_off = occ_off;
    thread_data_array[thread].virtpi = virtpi;
    thread_data_array[thread].vir_off = vir_off;
    thread_data_array[thread].omega = omega;
    thread_data_array[thread].fIJ = &fIJ;
    thread_data_array[thread].fAB = &fAB;
    thread_data_array[thread].outfile = outfile;
    thread_data_array[thread].thr_id = thread;
    thread_data_array[thread].SIA_local = SIA_local[thread];
    thread_data_array[thread].SIjAb_local = SIjAb_local[thread];
    thread_data_array[thread].newtrips = newtrips;
  }

  ijk_part = new int [nthreads];

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      Gij = Gji = Gi ^ Gj; 
      for(Gk=0; Gk < nirreps; Gk++) {
        Gkj = Gjk = Gk ^ Gj;
        Gik = Gki = Gi ^ Gk;
        Gijk = Gi ^ Gj ^ Gk;

        nijk = occpi[Gi] * occpi[Gj] * occpi[Gk];
        if (nijk == 0) continue;

        for (thread=0; thread<nthreads;++thread) {
          thread_data_array[thread].Gi = Gi;
          thread_data_array[thread].Gj = Gj;
          thread_data_array[thread].Gk = Gk;
          ijk_part[thread] = nijk / nthreads;
          if (thread < (nijk % nthreads)) ++ijk_part[thread];
          if (do_singles) { /* zero S's */
            for (h=0; h<nirreps;++h) 
              zero_mat(SIA_local[thread].matrix[h], SIA_local[thread].params->rowtot[h],
                SIA_local[thread].params->coltot[h^GS]);
          }
          if (do_doubles) {
            for (h=0;h<nirreps;++h)
              zero_mat( SIjAb_local[thread].matrix[h], SIjAb_local[thread].params->rowtot[h],
                SIjAb_local[thread].params->coltot[h^GS]);
          }
        }

        cnt = 0;
        for (thread=0; thread<nthreads; ++thread) {
          if (!ijk_part[thread]) continue;  // there are more threads than nijk
          thread_data_array[thread].first_ijk = cnt;
          cnt += ijk_part[thread];
          thread_data_array[thread].last_ijk = cnt-1;
        }

        /* execute threads */
        for (thread=0;thread<nthreads;++thread) {
          if (!ijk_part[thread]) continue;
          errcod = pthread_create(&(p_thread[thread]), NULL, cc3_sigma_RHF_ic_thread,
                 (void *) &thread_data_array[thread]);
          if (errcod) {
            fprintf(stderr,"pthread_create in cc3_sigma_RHF_ic failed\n");
            exit(PSI_RETURN_FAILURE);
          }
        }

        for (thread=0; thread<nthreads;++thread) {
          if (!ijk_part[thread]) continue;
          errcod = pthread_join(p_thread[thread], NULL);
          if (errcod) {
            fprintf(stderr,"pthread_join in cc3_sigma_RHF_ic failed\n");
            exit(PSI_RETURN_FAILURE);
          }
        }

        for (thread=0;thread<nthreads;++thread) {
          if (do_singles) {
            for(h=0;h<nirreps;++h)
              for(row=0; row < SIA->params->rowtot[h]; row++)
                for(col=0; col < SIA->params->coltot[h^GS]; col++)
                  SIA->matrix[h][row][col] += SIA_local[thread].matrix[h][row][col];
          }
          if (do_doubles) {
            for (h=0;h<nirreps;++h) {
              length = ((long) SIjAb->params->rowtot[h])*((long) SIjAb->params->coltot[h^GS]);
              if(length)
                C_DAXPY(length, 1.0, &(SIjAb_local[thread].matrix[h][0][0]), 1,
                  &(SIjAb->matrix[h][0][0]), 1);
            }
          }
        } /* end adding up S's */
      } /* Gk */
    } /* Gj */
  } /* Gi */

  /* close up files and update sigma vectors */
  dpd_file2_mat_close(&fIJ);
  dpd_file2_mat_close(&fAB);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fAB);
  dpd_file2_mat_close(FME);

  for (i=0;i<nthreads;++i) {
    if (do_singles) {
      sprintf(lbl, "%s %d", "CC3 SIA", i);
      dpd_file2_mat_close(&(SIA_local[i]));
      dpd_file2_close(&(SIA_local[i]));
    }
    if (do_doubles) {
      sprintf(lbl, "%s %d", "CC3 SIjAb", i);
      for (h=0;h<nirreps;++h)
        dpd_buf4_mat_irrep_close(&(SIjAb_local[i]),h);
      dpd_buf4_close(&(SIjAb_local[i]));
    }
  }
  delete [] ijk_part;

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(WAbEi, h);
    dpd_buf4_mat_irrep_close(WmAEf, h);
  }

  if (do_singles) {
    dpd_file2_mat_wrt(SIA);
    dpd_file2_mat_close(SIA);
    for(h=0; h < nirreps; h++)
      dpd_buf4_mat_irrep_close(Dints, h);
  }

  if (do_doubles) {
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_close(WMnIe, h);
      dpd_buf4_mat_irrep_wrt(SIjAb, h);
      dpd_buf4_mat_irrep_close(SIjAb, h);
    }
  }
  free(thread_data_array);
  free(p_thread);

#ifdef HAVE_MKL
  mkl_set_num_threads(old_threads);
#endif

}

void* cc3_sigma_RHF_ic_thread(void* thread_data_in)
{
  int Ga, Gb, Gc, Gd, Gab, Gbc, Gl, Gik, Gki, Gkj, Glk, Gjk, Gli, Gijk, nirreps;
  int GC, GWX3, GX3, GW, GS, Gij, Gji, Gid, Gjd, Gkd;
  int i, j, k, l, I, J, K, L, nrows, ncols, nlinks,h,cnt_ijk;
  int ij, ji, ik, ki, jk, kj, lc, li, il, lk, kl, id, jd, kd;
  int a, b, c, d, A, B, C, D, ab, ba, bc, ad, da, cnt;
  double **Z, *tvect, ***W3, ***W3a, ***W, ***V, ***Wa, ***Va;
  dpdbuf4 *SIjAb, SIjAb_local;
  dpdfile2 *SIA, SIA_local;
  char lbl[32];
  struct thread_data data;

  int do_singles, do_doubles, *occpi, *occ_off, *virtpi, *vir_off;
  int Gi, Gj, Gk, thr_id, first_ijk, last_ijk;
  double omega;
  dpdfile2 *FME, *fIJ, *fAB;
  dpdbuf4 *CIjAb, *WAbEi, *WMbIj, *Dints, *WmAEf, *WMnIe;
  FILE *outfile;
  int newtrips;
  int Gcb, Gac, cb, ac;

  data = *(struct thread_data *) thread_data_in;

  SIjAb  = data.SIjAb;
  SIA    = data.SIA;
  CIjAb  = data.CIjAb;
  WAbEi  = data.WAbEi;
  WMbIj  = data.WMbIj;
  do_singles = data.do_singles;
  Dints  = data.Dints;
  do_doubles = data.do_doubles;
  FME    = data.FME;
  WmAEf  = data.WmAEf;
  WMnIe  = data.WMnIe;
  occpi  = data.occpi;
  occ_off= data.occ_off;
  virtpi = data.virtpi;
  vir_off= data.vir_off;
  omega  = data.omega;
  fIJ  =  data.fIJ;
  fAB  =  data.fAB;
  Gi    =  data.Gi;
  Gj    =  data.Gj;
  Gk    =  data.Gk;
  first_ijk = data.first_ijk;
  last_ijk = data.last_ijk;
  outfile = data.outfile;
  thr_id = data.thr_id;
  SIA_local = data.SIA_local;
  SIjAb_local = data.SIjAb_local;
  newtrips = data.newtrips;

  nirreps = CIjAb->params->nirreps;
  Gij = Gji = Gi ^ Gj; 
  Gkj = Gjk = Gk ^ Gj;
  Gik = Gki = Gi ^ Gk;
  Gijk = Gi ^ Gj ^ Gk;
  nirreps = CIjAb->params->nirreps;
  GC = CIjAb->file.my_irrep;
  GWX3 = WAbEi->file.my_irrep;
  GX3 = GWX3^GC;
  GW = WmAEf->file.my_irrep;
  GS = SIjAb->file.my_irrep;

  W3 = (double ***) malloc(nirreps * sizeof(double **));
  W3a = (double ***) malloc(nirreps * sizeof(double **));
  W = (double ***) malloc(nirreps * sizeof(double **));
  V = (double ***) malloc(nirreps * sizeof(double **));
  Wa = (double ***) malloc(nirreps * sizeof(double **));
  Va = (double ***) malloc(nirreps * sizeof(double **));
  /* allocate memory for all irrep blocks of (ab,c) */
  for(Gab=0; Gab < nirreps; Gab++) {
    Gc = Gab ^ Gijk ^ GX3;
    W3[Gab] = dpd_block_matrix(WAbEi->params->coltot[Gab], virtpi[Gc]);
    if(newtrips) {
      W[Gab] = dpd_block_matrix(WAbEi->params->coltot[Gab], virtpi[Gc]);
      V[Gab] = dpd_block_matrix(WAbEi->params->coltot[Gab], virtpi[Gc]);
    }
  }
  for(Ga=0; Ga < nirreps; Ga++) {
    Gbc = Ga ^ Gijk ^ GX3;
    if(newtrips) {
      Wa[Ga] = dpd_block_matrix(virtpi[Ga], WAbEi->params->coltot[Gbc]);
      Va[Ga] = dpd_block_matrix(virtpi[Ga], WAbEi->params->coltot[Gbc]);
    }
    else
      W3a[Ga] = dpd_block_matrix(virtpi[Ga], WAbEi->params->coltot[Gbc]);
  }

  cnt_ijk = -1;
  for(i=0; i < occpi[Gi]; i++) {
    I = occ_off[Gi] + i;
    for(j=0; j < occpi[Gj]; j++) {
      J = occ_off[Gj] + j;
      for(k=0; k < occpi[Gk]; k++) {
        K = occ_off[Gk] + k;

        ++cnt_ijk;
        /* check to see if this ijk is for this thread */
        if ( (cnt_ijk < first_ijk) || (cnt_ijk > last_ijk))
           continue;

        ij = CIjAb->params->rowidx[I][J];
        ji = CIjAb->params->rowidx[J][I];
        ik = CIjAb->params->rowidx[I][K];
        ki = CIjAb->params->rowidx[K][I];
        jk = CIjAb->params->rowidx[J][K];
        kj = CIjAb->params->rowidx[K][J];

        T3_RHF_ic(W3, nirreps, I, Gi, J, Gj, K, Gk, CIjAb, WAbEi, WMbIj, 
          fIJ, fAB, occpi, occ_off, virtpi, vir_off, omega);

        if(newtrips) {
        /* generate lincombs of triples needed for singles and doubles */
        /* W(abc) = 2 T(abc) - T(cba) - T(acb) */
        /* V(abc) = T(abc) - T(cba) */
        for(Gab=0; Gab < nirreps; Gab++) {
          Gc = Gab ^ Gijk ^ GX3;
          for(ab=0; ab < WAbEi->params->coltot[Gab]; ab++) {
            A = WAbEi->params->colorb[Gab][ab][0];
            B = WAbEi->params->colorb[Gab][ab][1];
            Ga = WAbEi->params->rsym[A];
            Gb = Gab^Ga;
            a = A - vir_off[Ga];
            b = B - vir_off[Gb];
            Gcb = Gc ^ Gb;
            Gac = Ga ^ Gc;
            for(c=0; c < virtpi[Gc]; c++) {
              C = vir_off[Gc] + c;
              cb = WAbEi->params->colidx[C][B];
              ac = WAbEi->params->colidx[A][C];
              bc = WAbEi->params->colidx[B][C];
              W[Gab][ab][c] = 2*W3[Gab][ab][c]-W3[Gcb][cb][a]-W3[Gac][ac][b];
              V[Gab][ab][c] = W3[Gab][ab][c] - W3[Gcb][cb][a];
            }
          }
        }
        } /* newtrips */


        /* do (Wmnie*X3(ab,c) --> SIjAb) contractions that use W3(ab,c) first */
        if (do_doubles) {
          if(newtrips) {
          /* S_liba <-- -W_ijkabc W_jklc */
          /* S_ilab <-- -W_ijkabc W_jklc */
          for(Gl=0; Gl < nirreps; Gl++) {
            Gli = Gl ^ Gi;
            Gab = Gli ^ GS;
            Gc = Gjk ^ Gl ^ GW;

            lc = WMnIe->col_offset[Gjk][Gl];
            nrows = SIjAb_local.params->coltot[Gab];
            ncols = occpi[Gl];
            nlinks = virtpi[Gc];

            if(nrows && ncols && nlinks) {
              Z = dpd_block_matrix(nrows, ncols);

              C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W[Gab][0],
                  nlinks, &(WMnIe->matrix[Gjk][jk][lc]), nlinks,
                  0.0, Z[0], ncols);

              for(l=0; l < ncols; l++) {
                L = occ_off[Gl] + l;
                li = SIjAb_local.params->rowidx[L][I];
                il = SIjAb_local.params->rowidx[I][L];
                for(ab=0; ab < nrows; ab++) {
                  A = SIjAb_local.params->colorb[Gab][ab][0];
                  B = SIjAb_local.params->colorb[Gab][ab][1];
                  ba = SIjAb_local.params->colidx[B][A];
                  SIjAb_local.matrix[Gli][li][ba] -= Z[ab][l];
                  SIjAb_local.matrix[Gli][il][ab] -= Z[ab][l];
                }
              }
              dpd_free_block(Z, nrows, ncols);
            }
          }
          }  /* newtrips */
          else {

          /* S_liba <-- -2.0 * t_ijkabc W_jklc */
          /* S_ilab <-- -2.0 * t_ijkabc W_jklc */
          for(Gl=0; Gl < nirreps; Gl++) {
            Gli = Gl ^ Gi;
            Gab = Gli ^ GS;
            Gc = Gjk ^ Gl ^ GW;

            lc = WMnIe->col_offset[Gjk][Gl];
            nrows = SIjAb_local.params->coltot[Gab];
            ncols = occpi[Gl];
            nlinks = virtpi[Gc];

            if(nrows && ncols && nlinks) {
              Z = dpd_block_matrix(nrows, ncols);

              C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W3[Gab][0], nlinks,
                &(WMnIe->matrix[Gjk][jk][lc]), nlinks, 0.0, Z[0], ncols);

              for(l=0; l < ncols; l++) {
                L = occ_off[Gl] + l;
                li = SIjAb_local.params->rowidx[L][I];
                il = SIjAb_local.params->rowidx[I][L];
                for(ab=0; ab < nrows; ab++) {
                  A = SIjAb_local.params->colorb[Gab][ab][0];
                  B = SIjAb_local.params->colorb[Gab][ab][1];
                  ba = SIjAb_local.params->colidx[B][A];
                  SIjAb_local.matrix[Gli][li][ba] -= 2.0 * Z[ab][l];
                  SIjAb_local.matrix[Gli][il][ab] -= 2.0 * Z[ab][l];
                }
              }
              dpd_free_block(Z, nrows, ncols);
            }
          }

          /* S_lkba <-- + t_ijkabc W_jilc */
          /* S_klab <-- + t_ijkabc W_jilc */
          for(Gl=0; Gl < nirreps; Gl++) {
            Glk = Gl ^ Gk;
            Gab = Glk ^ GS;
            Gc = Gji ^ Gl ^ GW;

            lc = WMnIe->col_offset[Gji][Gl];
            nrows = SIjAb_local.params->coltot[Gab];
            ncols = occpi[Gl];
            nlinks = virtpi[Gc];

            if(nrows && ncols && nlinks) {
              Z = dpd_block_matrix(nrows, ncols);

              C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W3[Gab][0], nlinks,
                &(WMnIe->matrix[Gji][ji][lc]), nlinks, 0.0, Z[0], ncols);

              for(l=0; l < ncols; l++) {
                L = occ_off[Gl] + l;
                lk = SIjAb_local.params->rowidx[L][K];
                kl = SIjAb_local.params->rowidx[K][L];
                for(ab=0; ab < nrows; ab++) {
                  A = SIjAb_local.params->colorb[Gab][ab][0];
                  B = SIjAb_local.params->colorb[Gab][ab][1];
                  ba = SIjAb_local.params->colidx[B][A];
                  SIjAb_local.matrix[Glk][lk][ba] += Z[ab][l];
                  SIjAb_local.matrix[Glk][kl][ab] += Z[ab][l];
                }
              }
              dpd_free_block(Z, nrows, ncols);
            }
          }

          /* S_liba <-- + t_ijkabc W_kjlc */
          /* S_ilab <-- + t_ijkabc W_kjlc */
          for(Gl=0; Gl < nirreps; Gl++) {
            Gli = Gl ^ Gi;
            Gab = Gli ^ GS;
            Gc = Gjk ^ Gl ^ GW;

            lc = WMnIe->col_offset[Gkj][Gl];
            nrows = SIjAb_local.params->coltot[Gab];
            ncols = occpi[Gl];
            nlinks = virtpi[Gc];

            if(nrows && ncols && nlinks) {
              Z = dpd_block_matrix(nrows, ncols);

              C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W3[Gab][0], nlinks,
                &(WMnIe->matrix[Gkj][kj][lc]), nlinks, 0.0, Z[0], ncols);

              for(l=0; l < ncols; l++) {
                L = occ_off[Gl] + l;
                li = SIjAb_local.params->rowidx[L][I];
                il = SIjAb_local.params->rowidx[I][L];
                for(ab=0; ab < nrows; ab++) {
                  A = SIjAb_local.params->colorb[Gab][ab][0];
                  B = SIjAb_local.params->colorb[Gab][ab][1];
                  ba = SIjAb_local.params->colidx[B][A];
                  SIjAb_local.matrix[Gli][li][ba] += Z[ab][l];
                  SIjAb_local.matrix[Gli][il][ab] += Z[ab][l];
                }
              }
              dpd_free_block(Z, nrows, ncols);
            }
          }
        } /* no newtrips */

        } /* end Wmnie*X3 doubles contributions */

        if(newtrips) {
          /* sort W(ab,c) to W(a,bc) */
          for(Gab=0; Gab < nirreps; Gab++) {
            Gc = Gab ^ Gijk ^ GX3;
            for(ab=0; ab < WAbEi->params->coltot[Gab]; ab++ ){
              A = WAbEi->params->colorb[Gab][ab][0];
              B = WAbEi->params->colorb[Gab][ab][1];
              Ga = WAbEi->params->rsym[A];
              Gb = Gab^Ga;
              a = A - vir_off[Ga];
              for(c=0; c < virtpi[Gc]; c++) {
                C = vir_off[Gc] + c;
                bc = WAbEi->params->colidx[B][C];
                Wa[Ga][a][bc] = W[Gab][ab][c];
                Va[Ga][a][bc] = V[Gab][ab][c];
              }
            }
          }
        }
        else {
        /* sort W(ab,c) to W(a,bc) */
        for(Gab=0; Gab < nirreps; Gab++) {
          Gc = Gab ^ Gijk ^ GX3;
          for(ab=0; ab < WAbEi->params->coltot[Gab]; ab++ ){
            A = WAbEi->params->colorb[Gab][ab][0];
            B = WAbEi->params->colorb[Gab][ab][1];
            Ga = WAbEi->params->rsym[A];
            Gb = Gab^Ga;
            a = A - vir_off[Ga];
            for(c=0; c < virtpi[Gc]; c++) {
              C = vir_off[Gc] + c;
              bc = WAbEi->params->colidx[B][C];
              W3a[Ga][a][bc] = W3[Gab][ab][c];
            }
          }
        }
        }

        /*** X3(a,bc)*Dints --> SIA Contributions ***/
        if (do_singles) {
          /* Note: the do_singles code has not been used where Dints!=A1
          as this non-symmetric case does not arise for CC3 EOM energies */

          if(newtrips) {
            /* S_ia <-- t_ijkabc Djkbc */
            Ga = Gi ^ GS;
            Gbc = Gjk ^ GW;
            nrows = virtpi[Ga];
            ncols = Dints->params->coltot[Gbc];
            if(nrows && ncols)
              C_DGEMV('n', nrows, ncols, 1.0, Va[Ga][0], ncols,
                 &(Dints->matrix[Gjk][jk][0]), 1, 1.0,
                 SIA_local.matrix[Gi][i], 1);
          }
          else {
          /* S_ia <-- t_ijkabc Djkbc */
          Ga = Gi ^ GS;
          Gbc = Gjk ^ GW;
          nrows = virtpi[Ga];
          ncols = Dints->params->coltot[Gbc];

          if(nrows && ncols)
            C_DGEMV('n', nrows, ncols, 1.0, W3a[Ga][0], ncols,
              &(Dints->matrix[Gjk][jk][0]), 1, 1.0, SIA_local.matrix[Gi][i], 1);

          /* S_ka <-- tijkabc Djibc */
          Ga = Gk ^ GS;
          Gbc = Gji ^ GW;
          nrows = virtpi[Ga];
          ncols = Dints->params->coltot[Gbc];

          if(nrows && ncols)
            C_DGEMV('n', nrows, ncols, -1.0, W3a[Ga][0], ncols,
              &(Dints->matrix[Gji][ji][0]), 1, 1.0, SIA_local.matrix[Gk][k], 1);
          }
        }

        /*** X3(a,bc)*Wamef --> SIjAb Contributions ***/
        if (do_doubles) {

          if(newtrips) {
            /* S_IjAb <-- V_ijkabc F_kc */
            /* S_jIbA <-- V_ijkabc F_kc */
            Gc = Gk ^ GW;
            Gab = Gij ^ GS;
            nrows = SIjAb_local.params->coltot[Gij^GS];
            ncols = virtpi[Gc];

            if(nrows && ncols) {
              tvect = init_array(nrows);
              C_DGEMV('n', nrows, ncols, 1.0, V[Gab][0], ncols,
                      FME->matrix[Gk][k], 1, 0.0, tvect, 1);

              for(cnt=0; cnt<nrows; ++cnt) {
                A = SIjAb_local.params->colorb[Gab][cnt][0];
                B = SIjAb_local.params->colorb[Gab][cnt][1];
                ba = SIjAb_local.params->colidx[B][A];
                SIjAb_local.matrix[Gij][ij][cnt] += tvect[cnt];
                SIjAb_local.matrix[Gij][ji][ba] += tvect[cnt];
              }
              free(tvect);
            }
            /* S_ijad <-- W_ijkabc W_kdbc */
            /* S_jida <-- W_ijkabc W_kdbc */
            for(Gd=0; Gd < nirreps; Gd++) {                  
              Gkd = Gk ^ Gd;
              Ga = Gd ^ Gij ^ GS;

              nrows = virtpi[Ga];
              ncols = virtpi[Gd];
              nlinks = WmAEf->params->coltot[Gkd^GW];

              if(nrows && ncols && nlinks) {
                kd = WmAEf->row_offset[Gkd][K];

                Z = block_matrix(virtpi[Ga], virtpi[Gd]);
                C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, Wa[Ga][0],
                  nlinks, WmAEf->matrix[Gkd][kd], nlinks, 1.0, Z[0], ncols);

                for(a=0; a < virtpi[Ga]; a++) {
                  A = vir_off[Ga] + a;
                  for(d=0; d < virtpi[Gd]; d++) {
                    D = vir_off[Gd] + d;
                    ad = SIjAb_local.params->colidx[A][D];
                    da = SIjAb_local.params->colidx[D][A];
                    SIjAb_local.matrix[Gij][ij][ad] += Z[a][d];
                    SIjAb_local.matrix[Gij][ji][da] += Z[a][d];
                  }
                }

                free_block(Z);
              }
            }
          } /* newtrips */
          else {
          /* S_IjAb <-- t_ijkabc F_kc */
          /* S_jIbA <-- t_ijkabc F_kc */
          Gc = Gk ^ GW;
          Gab = Gij ^ GS;
          nrows = SIjAb_local.params->coltot[Gij^GS];
          ncols = virtpi[Gc];
              
          if(nrows && ncols) {
            tvect = init_array(nrows);
            C_DGEMV('n', nrows, ncols, 1.0, W3[Gab][0], ncols, FME->matrix[Gk][k], 1,
              0.0, tvect, 1);

            for(cnt=0; cnt<nrows; ++cnt) {
              A = SIjAb_local.params->colorb[Gab][cnt][0];
              B = SIjAb_local.params->colorb[Gab][cnt][1];
              ba = SIjAb_local.params->colidx[B][A];
              SIjAb_local.matrix[Gij][ij][cnt] += tvect[cnt];
              SIjAb_local.matrix[Gij][ji][ba] += tvect[cnt];
            }
            free(tvect);
          }

          /* S_KjAb <-- - t_ijkabc F_ic */
          /* S_jKbA <-- - t_ijkabc F_ic */
          Gc = Gi ^ GW;
          Gab = Gkj ^ GS;
          nrows = SIjAb_local.params->coltot[Gkj^GS];
          ncols = virtpi[Gc];
              
          if(nrows && ncols) {
            tvect = init_array(nrows);
            C_DGEMV('n', nrows, ncols, 1.0, W3[Gab][0], ncols, FME->matrix[Gi][i], 1,
              0.0, tvect, 1);

            for(cnt=0; cnt<nrows; ++cnt) {
              A = SIjAb_local.params->colorb[Gab][cnt][0];
              B = SIjAb_local.params->colorb[Gab][cnt][1];
              ba = SIjAb_local.params->colidx[B][A];
              SIjAb_local.matrix[Gkj][kj][cnt] -= tvect[cnt];
              SIjAb_local.matrix[Gkj][jk][ba] -= tvect[cnt];
            }
            free(tvect);
          }

          /* S_ijad <-- 2.0 * t_ijkabc W_kdbc */
          /* S_jida <-- 2.0 * t_ijkabc W_kdbc */
          for(Gd=0; Gd < nirreps; Gd++) {
            Gkd = Gk ^ Gd;
            Ga = Gd ^ Gij ^ GS;

            nrows = virtpi[Ga];
            ncols = virtpi[Gd];
            nlinks = WmAEf->params->coltot[Gkd^GW];

            if(nrows && ncols && nlinks) {
              kd = WmAEf->row_offset[Gkd][K];
               /*   WmAEf->matrix[Gkd] = dpd_block_matrix(virtpi[Gd], WmAEf->params->coltot[Gkd^GW]);
                  dpd_buf4_mat_irrep_rd_block(WmAEf, Gkd, kd, virtpi[Gd]); */

              Z = block_matrix(virtpi[Ga], virtpi[Gd]);
              C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W3a[Ga][0], nlinks,
                WmAEf->matrix[Gkd][kd], nlinks, 1.0, Z[0], ncols);

              for(a=0; a < virtpi[Ga]; a++) {
                A = vir_off[Ga] + a;
                for(d=0; d < virtpi[Gd]; d++) {
                  D = vir_off[Gd] + d;
                  ad = SIjAb_local.params->colidx[A][D];
                  da = SIjAb_local.params->colidx[D][A];
                  SIjAb_local.matrix[Gij][ij][ad] += 2.0 * Z[a][d];
                  SIjAb_local.matrix[Gij][ji][da] += 2.0 * Z[a][d];
                }
              }
              free_block(Z);
              /* dpd_free_block(WmAEf->matrix[Gkd], virtpi[Gd], WmAEf->params->coltot[Gkd^GW]); */
            }
          }

          /* S_kjad <-- - t_ijkabc W_idbc */
          /* S_jkda <-- - t_ijkabc W_idbc */
          for(Gd=0; Gd < nirreps; Gd++) {
            Gid = Gi ^ Gd;
            Ga = Gd ^ Gjk ^ GS;

            nrows = virtpi[Ga];
            ncols = virtpi[Gd];
            nlinks = WmAEf->params->coltot[Gid^GW];

            if(nrows && ncols && nlinks) {
              id = WmAEf->row_offset[Gid][I];
              /* WmAEf->matrix[Gid] = dpd_block_matrix(virtpi[Gd], WmAEf->params->coltot[Gid^GW]);
              dpd_buf4_mat_irrep_rd_block(WmAEf, Gid, id, virtpi[Gd]); */

              Z = block_matrix(virtpi[Ga], virtpi[Gd]);
              C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W3a[Ga][0], nlinks,
                WmAEf->matrix[Gid][id], nlinks, 1.0, Z[0], ncols);

              for(a=0; a < virtpi[Ga]; a++) {
                A = vir_off[Ga] + a;
                for(d=0; d < virtpi[Gd]; d++) {
                  D = vir_off[Gd] + d;
                  ad = SIjAb_local.params->colidx[A][D];
                  da = SIjAb_local.params->colidx[D][A];
                  SIjAb_local.matrix[Gkj][kj][ad] -= Z[a][d];
                  SIjAb_local.matrix[Gjk][jk][da] -= Z[a][d];
                }
              }

              free_block(Z);
              /* dpd_free_block(WmAEf->matrix[Gid], virtpi[Gd], WmAEf->params->coltot[Gid^GW]); */
            }
          }

          /* S_kida <-- - t_ijkabc W_jdbc */
          /* S_ikad <-- - t_ijkabc W_jdbc */
          for(Gd=0; Gd < nirreps; Gd++) {
            Gjd = Gj ^ Gd;
            Ga = Gd ^ Gik ^ GS;

            nrows = virtpi[Ga];
            ncols = virtpi[Gd];
            nlinks = WmAEf->params->coltot[Gjd^GW];

            if(nrows && ncols && nlinks) {
              jd = WmAEf->row_offset[Gjd][J];
              /* WmAEf->matrix[Gjd] = dpd_block_matrix(virtpi[Gd], WmAEf->params->coltot[Gjd^GW]);
              dpd_buf4_mat_irrep_rd_block(WmAEf, Gjd, jd, virtpi[Gd]); */

              Z = block_matrix(virtpi[Ga], virtpi[Gd]);
              C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W3a[Ga][0], nlinks,
                WmAEf->matrix[Gjd][jd], nlinks, 1.0, Z[0], ncols);

              for(a=0; a < virtpi[Ga]; a++) {
                A = vir_off[Ga] + a;
                for(d=0; d < virtpi[Gd]; d++) {
                  D = vir_off[Gd] + d;
                  ad = SIjAb_local.params->colidx[A][D];
                  da = SIjAb_local.params->colidx[D][A];
                  SIjAb_local.matrix[Gki][ki][da] -= Z[a][d];
                  SIjAb_local.matrix[Gik][ik][ad] -= Z[a][d];
                }
              }
              free_block(Z);
              /* dpd_free_block(WmAEf->matrix[Gjd], virtpi[Gd], WmAEf->params->coltot[Gjd^GW]); */
            }
          }
          } /* oldtrips */
        } /* end do_doubles */
      } /* k */
    } /* j */
  } /* i */

  for(Gab=0; Gab < nirreps; Gab++) {
    Gc = Gab ^ Gijk ^ GX3;
    dpd_free_block(W3[Gab], WAbEi->params->coltot[Gab], virtpi[Gc]);
    if(newtrips) {
      dpd_free_block(W[Gab], WAbEi->params->coltot[Gab], virtpi[Gc]);
      dpd_free_block(V[Gab], WAbEi->params->coltot[Gab], virtpi[Gc]);
    }
  }
  for(Ga=0; Ga < nirreps; Ga++) {
    Gbc = Ga ^ Gijk ^ GX3;
    if(newtrips) {
      dpd_free_block(Wa[Ga], virtpi[Ga], WAbEi->params->coltot[Gbc]);
      dpd_free_block(Va[Ga], virtpi[Ga], WAbEi->params->coltot[Gbc]);
    }
    else
      dpd_free_block(W3a[Ga], virtpi[Ga], WAbEi->params->coltot[Gbc]);
  }
  free(W3);
  free(W3a);
  free(W);
  free(V);
  free(Wa);
  free(Va);

  pthread_exit(NULL);
}

}
