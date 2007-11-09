/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <ccfiles.h>

extern "C" {

/*
  This function computes contributions to singles and doubles of
  matrix elements of triples:
    SIA   <-- <S|(Dints)           <T|(Wmbij,Wabei) CMNEF |0> |T> / (w-wt)
    SIjAb <-- <D|(FME,WmAEf,WMnIe) <T|(Wmbij,Wabei) CMNEF |0> |T> / (w-wt)
  Irrep variables are:
    GS <--            GW                  GWX3    ^ GC
                                  (           GX3            )
  These are used to make X3 quantity in T3_RHF:
    CIjAb, WAbEi, WMbIj, fIJ2, fAB2, omega
*/

void cc3_sigma_RHF(dpdbuf4 *CIjAb, dpdbuf4 *WAbEi, dpdbuf4 *WMbIj,
    int do_singles, dpdbuf4 *Dints, dpdfile2 *SIA,
    int do_doubles, dpdfile2 *FME, dpdbuf4 *WmAEf, dpdbuf4 *WMnIe,
    dpdbuf4 *SIjAb, int *occpi, int *occ_off, int *virtpi, int *vir_off,
    double omega, FILE *outfile)
{
  int h, nirreps;
  int Gi, Gj, Gk, Gl, Ga, Gb, Gc, Gd;
  int i, j, k, l, a, b, c, d;
  int I, J, K, L, A, B, C, D;
  int kj, jk, ji, ij, ik, ki;
  int Gkj, Gjk, Gji, Gij, Gik, Gki, Gkd;
  int Gijk, GS, GC, GWX3, GW, GX3, nrows,ncols;
  int ab, ba, ac, ca, bc, cb;
  int Gab, Gba, Gac, Gca, Gbc, Gcb, Gid, Gjd;
  int id, jd, kd, ad, bd, cd;
  int il, jl, kl, la, lb, lc, li, lk;
  int da, di, dj, dk;
  int Gad, Gdi, Gdj, Gdk, Glc, Gli, Glk,cnt;
  int nlinks;
  double value, F_val, t_val, E_val;
  double dijk, denom, *tvect, **Z;
  double value_ia, value_ka, denom_ia, denom_ka;
  dpdfile2 fIJ, fIJ2, fAB, fAB2, SIA_inc;
  dpdbuf4 SIjAb_inc, buf4_tmp;
  double ***W3, ***W3a;

  nirreps = CIjAb->params->nirreps;

  /* these are sent to T3 function */
  dpd_file2_init(&fIJ2, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB2, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_mat_init(&fIJ2);
  dpd_file2_mat_init(&fAB2);
  dpd_file2_mat_rd(&fIJ2);
  dpd_file2_mat_rd(&fAB2);

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
    dpd_file2_init(&SIA_inc, CC_TMP0, GS, 0, 1, "CC3 SIA");  /* T3->S1 increment */
    dpd_file2_mat_init(&SIA_inc);
  }
  if (do_doubles) {
    dpd_buf4_init(&SIjAb_inc, CC_TMP0, GS, 0, 5, 0, 5, 0, "CC3 SIjAb");
    dpd_buf4_scm(&SIjAb_inc, 0);
  }

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(CIjAb, h);
    dpd_buf4_mat_irrep_rd(CIjAb, h);
    dpd_buf4_mat_irrep_init(WMbIj, h);
    dpd_buf4_mat_irrep_rd(WMbIj, h);

    if (do_singles) {
      dpd_buf4_mat_irrep_init(Dints, h);
      dpd_buf4_mat_irrep_rd(Dints, h);
    }
    if (do_doubles) {
      dpd_buf4_mat_irrep_init(WMnIe, h);
      dpd_buf4_mat_irrep_rd(WMnIe, h);
    }
    dpd_buf4_mat_irrep_init(&SIjAb_inc, h);
  }

  W3 = (double ***) malloc(nirreps * sizeof(double **));
  W3a = (double ***) malloc(nirreps * sizeof(double **));

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      Gij = Gji = Gi ^ Gj; 
      for(Gk=0; Gk < nirreps; Gk++) {
        Gkj = Gjk = Gk ^ Gj;
        Gik = Gki = Gi ^ Gk;
        Gijk = Gi ^ Gj ^ Gk;
        
        /* allocate memory for all irrep blocks of (ab,c) */
        for(Gab=0; Gab < nirreps; Gab++) {
          Gc = Gab ^ Gijk ^ GX3;
          W3[Gab] = dpd_block_matrix(WAbEi->params->coltot[Gab], virtpi[Gc]);
        }
        for(Ga=0; Ga < nirreps; Ga++) {
          Gbc = Ga ^ Gijk ^ GX3;
          W3a[Ga] = dpd_block_matrix(virtpi[Ga], WAbEi->params->coltot[Gbc]);
        }
        
        for(i=0; i < occpi[Gi]; i++) {
          I = occ_off[Gi] + i;
          for(j=0; j < occpi[Gj]; j++) {
            J = occ_off[Gj] + j;
            for(k=0; k < occpi[Gk]; k++) {
              K = occ_off[Gk] + k;
              
              ij = CIjAb->params->rowidx[I][J];
              ji = CIjAb->params->rowidx[J][I];
              ik = CIjAb->params->rowidx[I][K];
              ki = CIjAb->params->rowidx[K][I];
              jk = CIjAb->params->rowidx[J][K];
              kj = CIjAb->params->rowidx[K][J];
#ifdef T3_TIMER_ON              
  timer_on("T3_RHF");
#endif
              T3_RHF(W3, nirreps, I, Gi, J, Gj, K, Gk, CIjAb, WAbEi, WMbIj, 
                     &fIJ2, &fAB2, occpi, occ_off, virtpi, vir_off, omega);
#ifdef T3_TIMER_ON              
  timer_off("T3_RHF");
#endif

              /* do (Wmnie*X3(ab,c) --> SIjAb) contractions that use W3(ab,c) first */
              if (do_doubles) {
#ifdef T3_TIMER_ON              
  timer_on("X3*Wmnie");
#endif
                /* S_liba <-- -2.0 * t_ijkabc W_jklc */
                /* S_ilab <-- -2.0 * t_ijkabc W_jklc */
                for(Gl=0; Gl < nirreps; Gl++) {
                  Gli = Gl ^ Gi;
                  Gab = Gli ^ GS;
                  Gc = Gjk ^ Gl ^ GW;

                  lc = WMnIe->col_offset[Gjk][Gl];
                  nrows = SIjAb_inc.params->coltot[Gab];
                  ncols = occpi[Gl];
                  nlinks = virtpi[Gc];

                  if(nrows && ncols && nlinks) {
                    Z = dpd_block_matrix(nrows, ncols);

                    C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W3[Gab][0], nlinks,
                        &(WMnIe->matrix[Gjk][jk][lc]), nlinks, 0.0, Z[0], ncols);

                    for(l=0; l < ncols; l++) {
                      L = occ_off[Gl] + l;
                      li = SIjAb_inc.params->rowidx[L][I];
                      il = SIjAb_inc.params->rowidx[I][L];
                      for(ab=0; ab < nrows; ab++) {
                        A = SIjAb_inc.params->colorb[Gab][ab][0];
                        B = SIjAb_inc.params->colorb[Gab][ab][1];
                        ba = SIjAb_inc.params->colidx[B][A];
                        SIjAb_inc.matrix[Gli][li][ba] -= 2.0 * Z[ab][l];
                        SIjAb_inc.matrix[Gli][il][ab] -= 2.0 * Z[ab][l];
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
                  nrows = SIjAb_inc.params->coltot[Gab];
                  ncols = occpi[Gl];
                  nlinks = virtpi[Gc];

                  if(nrows && ncols && nlinks) {
                    Z = dpd_block_matrix(nrows, ncols);

                    C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W3[Gab][0], nlinks,
                        &(WMnIe->matrix[Gji][ji][lc]), nlinks, 0.0, Z[0], ncols);

                    for(l=0; l < ncols; l++) {
                      L = occ_off[Gl] + l;
                      lk = SIjAb_inc.params->rowidx[L][K];
                      kl = SIjAb_inc.params->rowidx[K][L];
                      for(ab=0; ab < nrows; ab++) {
                        A = SIjAb_inc.params->colorb[Gab][ab][0];
                        B = SIjAb_inc.params->colorb[Gab][ab][1];
                        ba = SIjAb_inc.params->colidx[B][A];
                        SIjAb_inc.matrix[Glk][lk][ba] += Z[ab][l];
                        SIjAb_inc.matrix[Glk][kl][ab] += Z[ab][l];
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
                  nrows = SIjAb_inc.params->coltot[Gab];
                  ncols = occpi[Gl];
                  nlinks = virtpi[Gc];

                  if(nrows && ncols && nlinks) {
                    Z = dpd_block_matrix(nrows, ncols);

                    C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W3[Gab][0], nlinks,
                        &(WMnIe->matrix[Gkj][kj][lc]), nlinks, 0.0, Z[0], ncols);

                    for(l=0; l < ncols; l++) {
                      L = occ_off[Gl] + l;
                      li = SIjAb_inc.params->rowidx[L][I];
                      il = SIjAb_inc.params->rowidx[I][L];
                      for(ab=0; ab < nrows; ab++) {
                        A = SIjAb_inc.params->colorb[Gab][ab][0];
                        B = SIjAb_inc.params->colorb[Gab][ab][1];
                        ba = SIjAb_inc.params->colidx[B][A];
                        SIjAb_inc.matrix[Gli][li][ba] += Z[ab][l];
                        SIjAb_inc.matrix[Gli][il][ab] += Z[ab][l];
                      }
                    }
                    dpd_free_block(Z, nrows, ncols);
                  }
                }
#ifdef T3_TIMER_ON              
  timer_off("X3*Wmnie");
#endif
              } /* end Wmnie*X3 doubles contributions */

#ifdef T3_TIMER_ON              
  timer_on("X3_sort");
#endif
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
#ifdef T3_TIMER_ON              
  timer_off("X3_sort");
#endif

              /*** X3(a,bc)*Dints --> SIA Contributions ***/
              if (do_singles) {
                /* Note: the do_singles code has not been used where Dints!=A1
                   as this non-symmetric case does not arise for CC3 EOM energies */

                /* S_ia <-- t_ijkabc Djkbc */
                Ga = Gi ^ GS;
                Gbc = Gjk ^ GW;
                nrows = virtpi[Ga];
                ncols = Dints->params->coltot[Gbc];

                if(nrows && ncols)
                  C_DGEMV('n', nrows, ncols, 1.0, W3a[Ga][0], ncols,
                    &(Dints->matrix[Gjk][jk][0]), 1, 1.0, SIA_inc.matrix[Gi][i], 1);

                /* S_ka <-- tijkabc Djibc */
                Ga = Gk ^ GS;
                Gbc = Gji ^ GW;
                nrows = virtpi[Ga];
                ncols = Dints->params->coltot[Gbc];

                if(nrows && ncols)
                  C_DGEMV('n', nrows, ncols, -1.0, W3a[Ga][0], ncols,
                  &(Dints->matrix[Gji][ji][0]), 1, 1.0, SIA_inc.matrix[Gk][k], 1);
              }

              /*** X3(a,bc)*Wamef --> SIjAb Contributions ***/
              if (do_doubles) {
#ifdef T3_TIMER_ON              
  timer_on("X3*Wamef");
#endif
                /* S_IjAb <-- t_ijkabc F_kc */
                /* S_jIbA <-- t_ijkabc F_kc */
                Gc = Gk ^ GW;
                Gab = Gij ^ GS;
                nrows = SIjAb_inc.params->coltot[Gij^GS];
                ncols = virtpi[Gc];
                
                if(nrows && ncols) {
                  tvect = init_array(nrows);
                  C_DGEMV('n', nrows, ncols, 1.0, W3[Gab][0], ncols, FME->matrix[Gk][k], 1,
                    0.0, tvect, 1);

                  for(cnt=0; cnt<nrows; ++cnt) {
                    A = SIjAb_inc.params->colorb[Gab][cnt][0];
                    B = SIjAb_inc.params->colorb[Gab][cnt][1];
                    ba = SIjAb_inc.params->colidx[B][A];
                    SIjAb_inc.matrix[Gij][ij][cnt] += tvect[cnt];
                    SIjAb_inc.matrix[Gij][ji][ba] += tvect[cnt];
                  }
                  free(tvect);
                }

                /* S_KjAb <-- - t_ijkabc F_ic */
                /* S_jKbA <-- - t_ijkabc F_ic */
                Gc = Gi ^ GW;
                Gab = Gkj ^ GS;
                nrows = SIjAb_inc.params->coltot[Gkj^GS];
                ncols = virtpi[Gc];
                
                if(nrows && ncols) {
                  tvect = init_array(nrows);
                  C_DGEMV('n', nrows, ncols, 1.0, W3[Gab][0], ncols, FME->matrix[Gi][i], 1,
                      0.0, tvect, 1);

                  for(cnt=0; cnt<nrows; ++cnt) {
                    A = SIjAb_inc.params->colorb[Gab][cnt][0];
                    B = SIjAb_inc.params->colorb[Gab][cnt][1];
                    ba = SIjAb_inc.params->colidx[B][A];
                    SIjAb_inc.matrix[Gkj][kj][cnt] -= tvect[cnt];
                    SIjAb_inc.matrix[Gkj][jk][ba] -= tvect[cnt];
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
                    WmAEf->matrix[Gkd] = dpd_block_matrix(virtpi[Gd], WmAEf->params->coltot[Gkd^GW]);
                    dpd_buf4_mat_irrep_rd_block(WmAEf, Gkd, kd, virtpi[Gd]);

                    Z = block_matrix(virtpi[Ga], virtpi[Gd]);
                    C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W3a[Ga][0], nlinks,
                      WmAEf->matrix[Gkd][0], nlinks, 1.0, Z[0], ncols);

                    for(a=0; a < virtpi[Ga]; a++) {
                      A = vir_off[Ga] + a;
                      for(d=0; d < virtpi[Gd]; d++) {
                        D = vir_off[Gd] + d;
                        ad = SIjAb_inc.params->colidx[A][D];
                        da = SIjAb_inc.params->colidx[D][A];
                        SIjAb_inc.matrix[Gij][ij][ad] += 2.0 * Z[a][d];
                        SIjAb_inc.matrix[Gij][ji][da] += 2.0 * Z[a][d];
                      }
                    }

                    free_block(Z);
                    dpd_free_block(WmAEf->matrix[Gkd], virtpi[Gd], WmAEf->params->coltot[Gkd^GW]);
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
                    WmAEf->matrix[Gid] = dpd_block_matrix(virtpi[Gd], WmAEf->params->coltot[Gid^GW]);
                    dpd_buf4_mat_irrep_rd_block(WmAEf, Gid, id, virtpi[Gd]);

                    Z = block_matrix(virtpi[Ga], virtpi[Gd]);
                    C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W3a[Ga][0], nlinks,
                      WmAEf->matrix[Gid][0], nlinks, 1.0, Z[0], ncols);

                    for(a=0; a < virtpi[Ga]; a++) {
                      A = vir_off[Ga] + a;
                      for(d=0; d < virtpi[Gd]; d++) {
                        D = vir_off[Gd] + d;
                        ad = SIjAb_inc.params->colidx[A][D];
                        da = SIjAb_inc.params->colidx[D][A];
                        SIjAb_inc.matrix[Gkj][kj][ad] -= Z[a][d];
                        SIjAb_inc.matrix[Gjk][jk][da] -= Z[a][d];
                      }
                    }

                    free_block(Z);
                    dpd_free_block(WmAEf->matrix[Gid], virtpi[Gd], WmAEf->params->coltot[Gid^GW]);
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
                    WmAEf->matrix[Gjd] = dpd_block_matrix(virtpi[Gd], WmAEf->params->coltot[Gjd^GW]);
                    dpd_buf4_mat_irrep_rd_block(WmAEf, Gjd, jd, virtpi[Gd]);

                    Z = block_matrix(virtpi[Ga], virtpi[Gd]);
                    C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W3a[Ga][0], nlinks,
                      WmAEf->matrix[Gjd][0], nlinks, 1.0, Z[0], ncols);

                    for(a=0; a < virtpi[Ga]; a++) {
                      A = vir_off[Ga] + a;
                      for(d=0; d < virtpi[Gd]; d++) {
                        D = vir_off[Gd] + d;
                        ad = SIjAb_inc.params->colidx[A][D];
                        da = SIjAb_inc.params->colidx[D][A];
                        SIjAb_inc.matrix[Gki][ki][da] -= Z[a][d];
                        SIjAb_inc.matrix[Gik][ik][ad] -= Z[a][d];
                      }
                    }
                    free_block(Z);
                    dpd_free_block(WmAEf->matrix[Gjd], virtpi[Gd], WmAEf->params->coltot[Gjd^GW]);
                  }
                }
#ifdef T3_TIMER_ON              
timer_off("X3*Wamef");
#endif
              } /* end Wamef*X3->Sijab contributions */
            } /* k */
          } /* j */
         } /* i */

        for(Gab=0; Gab < nirreps; Gab++) {
          Gc = Gab ^ Gijk ^ GX3;
          dpd_free_block(W3[Gab], WAbEi->params->coltot[Gab], virtpi[Gc]);
        }
        for(Ga=0; Ga < nirreps; Ga++) {
          Gbc = Ga ^ Gijk ^ GX3;
          dpd_free_block(W3a[Ga], virtpi[Ga], WAbEi->params->coltot[Gbc]);
        }
      } /* Gk */
    } /* Gj */
  } /* Gi */

  /* close up files and update sigma vectors */
  dpd_file2_mat_close(&fIJ);
  dpd_file2_mat_close(&fAB);
  dpd_file2_mat_close(FME);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fAB);

  dpd_file2_mat_close(&fIJ2);
  dpd_file2_mat_close(&fAB2);
  dpd_file2_close(&fIJ2);
  dpd_file2_close(&fAB2);

  if (do_singles) {
    for(h=0; h < nirreps; h++)
      dpd_buf4_mat_irrep_close(Dints, h);
    dpd_file2_mat_wrt(&SIA_inc);
    dpd_file2_mat_close(&SIA_inc);
    dpd_file2_axpy(&SIA_inc, SIA, 1, 0);
    dpd_file2_close(&SIA_inc);
  }

  if (do_doubles) {
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_close(WMnIe, h);
      dpd_buf4_mat_irrep_wrt(&SIjAb_inc, h);
      dpd_buf4_mat_irrep_close(&SIjAb_inc, h);
    }
    dpd_buf4_axpy(&SIjAb_inc, SIjAb, 1);
    dpd_buf4_close(&SIjAb_inc);
  }
}

} /* extern "C" */
