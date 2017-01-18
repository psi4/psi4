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

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace cctriples {

extern void T3_UHF_AAA_abc(double ***W, double ***V, int disc, int nirreps,
    int A, int Ga, int B, int Gb, int C, int Gc, dpdbuf4 *C2, dpdbuf4 *F,
    dpdbuf4 *E, dpdfile2 *C1, dpdbuf4 *D, dpdfile2 *fIA, dpdfile2 *fIJ,
    dpdfile2 *fAB, int *occpi, int *occ_off, int *virtpi, int *vir_off,
    double omega);

extern void T3_UHF_AAB_abc(double ***W, double ***V, int disc, int nirreps,
    int I, int Gi, int J, int Gj, int K, int Gk, dpdbuf4 *T2AA, dpdbuf4 *T2AB,
    dpdbuf4 *T2BA, dpdbuf4 *FAA, dpdbuf4 *FAB, dpdbuf4 *FBA, dpdbuf4 *EAA,
    dpdbuf4 *EAB, dpdbuf4 *EBA, dpdfile2 *T1A, dpdfile2 *T1B, dpdbuf4 *DAA,
    dpdbuf4 *DAB, dpdfile2 *fIA, dpdfile2 *fia, dpdfile2 *fIJ, dpdfile2 *fij,
    dpdfile2 *fAB, dpdfile2 *fab, int *aoccpi, int *aocc_off, int *boccpi,
    int *bocc_off, int *avirtpi, int *avir_off, int *bvirtpi, int *bvir_off,
    double omega);

void test_abc_loops_AAA() {
  int i, j, k, a, b, c, ij, I, J, K;
  int A, B, C;
  int Ga, Gb, Gc, Gi, Gj, Gij, Gk;
  dpdfile2 fIJ, fAB, fIA, T1;

  double denom, ET;
  int nirreps = moinfo.nirreps;
  int *occpi = moinfo.aoccpi;
  int *virtpi = moinfo.avirtpi;
  int *occ_off = moinfo.aocc_off;
  int *vir_off = moinfo.avir_off;
  double ***WIJK = (double ***) malloc(nirreps * sizeof(double **));
  double ***VIJK = (double ***) malloc(nirreps * sizeof(double **));

  global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");

  dpdbuf4 EAAints;
  global_dpd_->buf4_init(&EAAints, PSIF_CC_EINTS, 0, 21, 0, 21, 2, 0, "E <AK||IJ> (AK, I>J)");
  dpdbuf4 FAAints;
  global_dpd_->buf4_init(&FAAints, PSIF_CC_FINTS, 0, 5, 20, 7, 20, 0, "F <BC||IA>");
  dpdbuf4 T2AA;
  global_dpd_->buf4_init(&T2AA, PSIF_CC_TAMPS, 0, 5, 0, 7, 2, 0, "tABIJ");
  dpdbuf4 DAAints;
  global_dpd_->buf4_init(&DAAints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");

  int Gabc;
  ET = 0.0;
  for (Ga=0; Ga < nirreps; ++Ga) {
    for (a=0; a<virtpi[Ga]; ++a) {
      A = vir_off[Ga] + a;
      for (Gb=0; Gb < nirreps; ++Gb) {
        for (b=0; b<virtpi[Gb]; ++b) {
          B = vir_off[Gb] + b;
          for (Gc=0; Gc < nirreps; ++Gc) {
            for (c=0; c < virtpi[Gc]; ++c) {
              C = vir_off[Gc] + c;
              Gabc = Ga ^ Gb ^ Gc;
              //Allocate the memory for connected and disconnected triples
              for (Gij=0; Gij < nirreps; ++Gij) {
                Gk = Gij ^ Gabc;
                WIJK[Gij] = global_dpd_->dpd_block_matrix(T2AA.params->coltot[Gij],
                    occpi[Gk]);
                VIJK[Gij] = global_dpd_->dpd_block_matrix(T2AA.params->coltot[Gij],
                    occpi[Gk]);
              }

              T3_UHF_AAA_abc(WIJK, VIJK, 1, nirreps, A, Ga, B, Gb, C, Gc,
                  &T2AA, &FAAints, &EAAints, &T1, &DAAints, &fIA, &fIJ, &fAB,
                  occpi, occ_off, virtpi, vir_off, 0.0);
              global_dpd_->file2_mat_init(&fIJ);
              global_dpd_->file2_mat_init(&fIA);
              global_dpd_->file2_mat_init(&fAB);
              global_dpd_->file2_mat_init(&T1);
              global_dpd_->file2_mat_rd(&fIJ);
              global_dpd_->file2_mat_rd(&fIA);
              global_dpd_->file2_mat_rd(&fAB);
              global_dpd_->file2_mat_rd(&T1);
              double dabc = 0.0;
              if (fAB.params->rowtot[Ga])
                dabc -= fAB.matrix[Ga][a][a];
              if (fAB.params->rowtot[Gb])
                dabc -= fAB.matrix[Gb][b][b];
              if (fAB.params->rowtot[Gc])
                dabc -= fAB.matrix[Gc][c][c];

              for (Gij=0; Gij < nirreps; ++Gij) {
                Gk = Gabc ^ Gij;
                for (ij=0; ij<T2AA.params->coltot[Gij]; ++ij) {
                  I = T2AA.params->colorb[Gij][ij][0];
                  Gi = T2AA.params->rsym[I];
                  i = I-occ_off[Gi];
                  J = T2AA.params->colorb[Gij][ij][1];
                  Gj = T2AA.params->ssym[J];
                  j = J-occ_off[Gj];
                  for (k=0; k<occpi[Gk]; ++k) {
                    denom = dabc;
                    K = k + occ_off[Gk];
                    if (fIJ.params->rowtot[Gi])
                      denom += fIJ.matrix[Gi][i][i];
                    if (fIJ.params->rowtot[Gj])
                      denom += fIJ.matrix[Gj][j][j];
                    if (fIJ.params->rowtot[Gk])
                      denom += fIJ.matrix[Gk][k][k];
                    ET += WIJK[Gij][ij][k] * (WIJK[Gij][ij][k]
                        +VIJK[Gij][ij][k]) * denom;
                  }
                }
              }
              global_dpd_->file2_mat_close(&fIJ);
              global_dpd_->file2_mat_close(&fAB);
              global_dpd_->file2_mat_close(&fIA);
              global_dpd_->file2_mat_close(&T1);

              //Deallocate the memory for connected and disconnected triples
              for (Gij=0; Gij < nirreps; ++Gij) {
                Gk = Gij ^ Gabc;
                global_dpd_->free_dpd_block(WIJK[Gij], T2AA.params->coltot[Gij], occpi[Gk]);
                global_dpd_->free_dpd_block(VIJK[Gij], T2AA.params->coltot[Gij], occpi[Gk]);
              }
            }
          }
        }
      }
    }
  }
  outfile->Printf( "    E(T) AAA from Andy's code = %20.15f\n", ET/36.0);
  global_dpd_->buf4_close(&EAAints);
  global_dpd_->buf4_close(&FAAints);
  global_dpd_->buf4_close(&T2AA);
  global_dpd_->file2_close(&T1);
  global_dpd_->file2_close(&fIJ);
  global_dpd_->file2_close(&fIA);
  global_dpd_->file2_close(&fAB);
}

void test_abc_loops_BBB() {
  int i, j, k, a, b, c, ij, I, J, K;
  int A, B, C;
  int Ga, Gb, Gc, Gi, Gj, Gij, Gk;
  dpdfile2 fij, fab, fia, T1;

  double denom, ET;
  int nirreps = moinfo.nirreps;
  int *occpi = moinfo.boccpi;
  int *virtpi = moinfo.bvirtpi;
  int *occ_off = moinfo.bocc_off;
  int *vir_off = moinfo.bvir_off;
  double ***WIJK = (double ***) malloc(nirreps * sizeof(double **));
  double ***VIJK = (double ***) malloc(nirreps * sizeof(double **));

  global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij");
  global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab");
  global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");

  dpdbuf4 EBBints;
  global_dpd_->buf4_init(&EBBints, PSIF_CC_EINTS, 0, 31, 10, 31, 12, 0,
      "E <ak||ij> (ak, i>j)");
  dpdbuf4 FBBints;
  global_dpd_->buf4_init(&FBBints, PSIF_CC_FINTS, 0, 15, 30, 17, 30, 0, "F <bc||ia>");
  dpdbuf4 T2BB;
  global_dpd_->buf4_init(&T2BB, PSIF_CC_TAMPS, 0, 15, 10, 17, 12, 0, "tabij");
  dpdbuf4 DBBints;
  global_dpd_->buf4_init(&DBBints, PSIF_CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");

  int Gabc;
  ET = 0.0;
  for (Ga=0; Ga < nirreps; ++Ga) {
    for (a=0; a<virtpi[Ga]; ++a) {
      A = vir_off[Ga] + a;
      for (Gb=0; Gb < nirreps; ++Gb) {
        for (b=0; b<virtpi[Gb]; ++b) {
          B = vir_off[Gb] + b;
          for (Gc=0; Gc < nirreps; ++Gc) {
            for (c=0; c < virtpi[Gc]; ++c) {
              C = vir_off[Gc] + c;
              Gabc = Ga ^ Gb ^ Gc;
              //Allocate the memory for connected and disconnected triples
              for (Gij=0; Gij < nirreps; ++Gij) {
                Gk = Gij ^ Gabc;
                WIJK[Gij] = global_dpd_->dpd_block_matrix(T2BB.params->coltot[Gij],
                    occpi[Gk]);
                VIJK[Gij] = global_dpd_->dpd_block_matrix(T2BB.params->coltot[Gij],
                    occpi[Gk]);
              }

              T3_UHF_AAA_abc(WIJK, VIJK, 1, nirreps, A, Ga, B, Gb, C, Gc,
                  &T2BB, &FBBints, &EBBints, &T1, &DBBints, &fia, &fij, &fab,
                  occpi, occ_off, virtpi, vir_off, 0.0);
              global_dpd_->file2_mat_init(&fij);
              global_dpd_->file2_mat_init(&fia);
              global_dpd_->file2_mat_init(&fab);
              global_dpd_->file2_mat_init(&T1);
              global_dpd_->file2_mat_rd(&fij);
              global_dpd_->file2_mat_rd(&fia);
              global_dpd_->file2_mat_rd(&fab);
              global_dpd_->file2_mat_rd(&T1);
              double dabc = 0.0;
              if (fab.params->rowtot[Ga])
                dabc -= fab.matrix[Ga][a][a];
              if (fab.params->rowtot[Gb])
                dabc -= fab.matrix[Gb][b][b];
              if (fab.params->rowtot[Gc])
                dabc -= fab.matrix[Gc][c][c];

              for (Gij=0; Gij < nirreps; ++Gij) {
                Gk = Gabc ^ Gij;
                for (ij=0; ij<T2BB.params->coltot[Gij]; ++ij) {
                  I = T2BB.params->colorb[Gij][ij][0];
                  Gi = T2BB.params->rsym[I];
                  i = I-occ_off[Gi];
                  J = T2BB.params->colorb[Gij][ij][1];
                  Gj = T2BB.params->ssym[J];
                  j = J-occ_off[Gj];
                  for (k=0; k<occpi[Gk]; ++k) {
                    denom = dabc;
                    K = k + occ_off[Gk];
                    if (fij.params->rowtot[Gi])
                      denom += fij.matrix[Gi][i][i];
                    if (fij.params->rowtot[Gj])
                      denom += fij.matrix[Gj][j][j];
                    if (fij.params->rowtot[Gk])
                      denom += fij.matrix[Gk][k][k];
                    ET += WIJK[Gij][ij][k] * (WIJK[Gij][ij][k]
                        +VIJK[Gij][ij][k]) * denom;
                  }
                }
              }
              global_dpd_->file2_mat_close(&fij);
              global_dpd_->file2_mat_close(&fab);
              global_dpd_->file2_mat_close(&fia);
              global_dpd_->file2_mat_close(&T1);

              //Deallocate the memory for connected and disconnected triples
              for (Gij=0; Gij < nirreps; ++Gij) {
                Gk = Gij ^ Gabc;
                global_dpd_->free_dpd_block(WIJK[Gij], T2BB.params->coltot[Gij], occpi[Gk]);
                global_dpd_->free_dpd_block(VIJK[Gij], T2BB.params->coltot[Gij], occpi[Gk]);
              }
            }
          }
        }
      }
    }
  }
  outfile->Printf( "    E(T) BBB from Andy's code = %20.15f\n", ET/36.0);
  global_dpd_->buf4_close(&EBBints);
  global_dpd_->buf4_close(&FBBints);
  global_dpd_->buf4_close(&T2BB);
  global_dpd_->file2_close(&T1);
  global_dpd_->file2_close(&fij);
  global_dpd_->file2_close(&fia);
  global_dpd_->file2_close(&fab);
}

void test_abc_loops_AAB() {
  int i, j, k, a, b, c, ij, I, J, K;
  int A, B, C;
  int Ga, Gb, Gc, Gi, Gj, Gij, Gk;
  dpdfile2 fij, fIJ, fab, fAB, fia, fIA, T1A, T1B;

  double denom, ET;
  int nirreps = moinfo.nirreps;
  int *aoccpi = moinfo.aoccpi;
  int *avirtpi = moinfo.avirtpi;
  int *aocc_off = moinfo.aocc_off;
  int *avir_off = moinfo.avir_off;
  int *boccpi = moinfo.boccpi;
  int *bvirtpi = moinfo.bvirtpi;
  int *bocc_off = moinfo.bocc_off;
  int *bvir_off = moinfo.bvir_off;
  double ***WIJk = (double ***) malloc(nirreps * sizeof(double **));
  double ***VIJk = (double ***) malloc(nirreps * sizeof(double **));

  global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij");
  global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab");
  global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
  global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");

  dpdbuf4 EAAints;
  dpdbuf4 EABints;
  dpdbuf4 EBAints;
  dpdbuf4 DAAints;
  dpdbuf4 DABints;
  dpdbuf4 FAAints;
  dpdbuf4 FABints;
  dpdbuf4 FBAints;
  dpdbuf4 T2AA;
  dpdbuf4 T2BA;
  dpdbuf4 T2AB;
  global_dpd_->buf4_init(&T2AA, PSIF_CC_TAMPS, 0, 5, 0, 7, 2, 0, "tABIJ");
  global_dpd_->buf4_init(&T2AB, PSIF_CC_TAMPS, 0, 28, 22, 28, 22, 0, "tAbIj");
  global_dpd_->buf4_init(&T2BA, PSIF_CC_TAMPS, 0, 29, 23, 29, 23, 0, "taBiJ");
  global_dpd_->buf4_init(&FAAints, PSIF_CC_FINTS, 0, 5, 20, 7, 20, 0, "F <BC||IA>");
  global_dpd_->buf4_init(&FBAints, PSIF_CC_FINTS, 0, 29, 27, 29, 27, 0, "F <bC|iA>");
  global_dpd_->buf4_init(&FABints, PSIF_CC_FINTS, 0, 28, 24, 28, 24, 0, "F <Bc|Ia>");
  global_dpd_->buf4_init(&EAAints, PSIF_CC_EINTS, 0, 21, 0, 21, 2, 0, "E <AK||IJ> (AK, I>J)");
  global_dpd_->buf4_init(&EABints, PSIF_CC_EINTS, 0, 25, 22, 25, 22, 0, "E <aK|Ij>");
  global_dpd_->buf4_init(&EBAints, PSIF_CC_EINTS, 0, 26, 23, 26, 23, 0, "E <Ak|iJ>");
  global_dpd_->buf4_init(&DAAints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
  global_dpd_->buf4_init(&DABints, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");

  int Gabc;
  ET = 0.0;
  for (Ga=0; Ga < nirreps; ++Ga) {
    for (a=0; a<avirtpi[Ga]; ++a) {
      A = avir_off[Ga] + a;
      for (Gb=0; Gb < nirreps; ++Gb) {
        for (b=0; b<avirtpi[Gb]; ++b) {
          B = avir_off[Gb] + b;
          for (Gc=0; Gc < nirreps; ++Gc) {
            for (c=0; c < bvirtpi[Gc]; ++c) {
              C = bvir_off[Gc] + c;
              Gabc = Ga ^ Gb ^ Gc;
              //Allocate the memory for connected and disconnected triples
              for (Gij=0; Gij < nirreps; ++Gij) {
                Gk = Gij ^ Gabc;
                WIJk[Gij] = global_dpd_->dpd_block_matrix(T2AA.params->coltot[Gij],
                    boccpi[Gk]);
                VIJk[Gij] = global_dpd_->dpd_block_matrix(T2AA.params->coltot[Gij],
                    boccpi[Gk]);
              }
              T3_UHF_AAB_abc(WIJk, VIJk, 1, nirreps, A, Ga, B, Gb, C, Gc,
                  &T2AA, &T2AB, &T2BA, &FAAints, &FABints, &FBAints, &EAAints,
                  &EABints, &EBAints, &T1A, &T1B, &DAAints, &DABints, &fIA,
                  &fia, &fIJ, &fij, &fAB, &fab, aoccpi, aocc_off, boccpi,
                  bocc_off, avirtpi, avir_off, bvirtpi, bvir_off, 0.0);
              global_dpd_->file2_mat_init(&fij);
              global_dpd_->file2_mat_init(&fIJ);
              global_dpd_->file2_mat_init(&fab);
              global_dpd_->file2_mat_init(&fAB);
              global_dpd_->file2_mat_rd(&fij);
              global_dpd_->file2_mat_rd(&fIJ);
              global_dpd_->file2_mat_rd(&fab);
              global_dpd_->file2_mat_rd(&fAB);

              double dabc = 0.0;
              if (fAB.params->rowtot[Ga])
                dabc -= fAB.matrix[Ga][a][a];
              if (fAB.params->rowtot[Gb])
                dabc -= fAB.matrix[Gb][b][b];
              if (fab.params->rowtot[Gc])
                dabc -= fab.matrix[Gc][c][c];

              for (Gij=0; Gij < nirreps; ++Gij) {
                Gk = Gabc ^ Gij;
                for (ij=0; ij<T2AA.params->coltot[Gij]; ++ij) {
                  I = T2AA.params->colorb[Gij][ij][0];
                  Gi = T2AA.params->rsym[I];
                  i = I-aocc_off[Gi];
                  J = T2AA.params->colorb[Gij][ij][1];
                  Gj = T2AA.params->ssym[J];
                  j = J-aocc_off[Gj];
                  for (k=0; k<boccpi[Gk]; ++k) {
                    denom = dabc;
                    K = k + bocc_off[Gk];
                    if (fIJ.params->rowtot[Gi])
                      denom += fIJ.matrix[Gi][i][i];
                    if (fIJ.params->rowtot[Gj])
                      denom += fIJ.matrix[Gj][j][j];
                    if (fij.params->rowtot[Gk])
                      denom += fij.matrix[Gk][k][k];
                    //							  outfile->Printf( "%d,%d,%d,%d,%d,%d = %20.12f\n",I,J,K,A,B,C,VIJk[Gij][ij][k]);


                    ET += WIJk[Gij][ij][k] * (WIJk[Gij][ij][k]
                        +VIJk[Gij][ij][k]) * denom;
                  }
                }
              }
              global_dpd_->file2_mat_close(&fij);
              global_dpd_->file2_mat_close(&fIJ);
              global_dpd_->file2_mat_close(&fab);
              global_dpd_->file2_mat_close(&fAB);

              //Deallocate the memory for connected and disconnected triples
              for (Gij=0; Gij < nirreps; ++Gij) {
                Gk = Gij ^ Gabc;
                global_dpd_->free_dpd_block(WIJk[Gij], T2AA.params->coltot[Gij], boccpi[Gk]);
                global_dpd_->free_dpd_block(VIJk[Gij], T2AA.params->coltot[Gij], boccpi[Gk]);
              }
            }
          }
        }
      }
    }
  }

  outfile->Printf( "    E(T) AAB from Andy's code = %20.15f\n", ET/4.0);
  global_dpd_->file2_close(&fij);
  global_dpd_->file2_close(&fIJ);
  global_dpd_->file2_close(&fab);
  global_dpd_->file2_close(&fAB);

  global_dpd_->buf4_close(&T2AA);
  global_dpd_->buf4_close(&T2AB);
  global_dpd_->buf4_close(&T2BA);
  global_dpd_->buf4_close(&EAAints);
  global_dpd_->buf4_close(&EABints);
  global_dpd_->buf4_close(&EBAints);
  global_dpd_->buf4_close(&FAAints);
  global_dpd_->buf4_close(&FABints);
  global_dpd_->buf4_close(&FBAints);
  global_dpd_->buf4_close(&DAAints);
  global_dpd_->buf4_close(&DABints);
}

void test_abc_loops_BBA() {
  int i, j, k, a, b, c, ij, I, J, K;
  int A, B, C;
  int Ga, Gb, Gc, Gi, Gj, Gij, Gk;
  dpdfile2 fij, fIJ, fab, fAB, fia, fIA, T1A, T1B;

  double denom, ET;
  int nirreps = moinfo.nirreps;
  int *aoccpi = moinfo.aoccpi;
  int *avirtpi = moinfo.avirtpi;
  int *aocc_off = moinfo.aocc_off;
  int *avir_off = moinfo.avir_off;
  int *boccpi = moinfo.boccpi;
  int *bvirtpi = moinfo.bvirtpi;
  int *bocc_off = moinfo.bocc_off;
  int *bvir_off = moinfo.bvir_off;
  double ***WijK = (double ***) malloc(nirreps * sizeof(double **));
  double ***VijK = (double ***) malloc(nirreps * sizeof(double **));

  global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij");
  global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab");
  global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
  global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");

  dpdbuf4 EBBints;
  dpdbuf4 EABints;
  dpdbuf4 EBAints;
  dpdbuf4 DBBints;
  dpdbuf4 DBAints;
  dpdbuf4 FBBints;
  dpdbuf4 FABints;
  dpdbuf4 FBAints;
  dpdbuf4 T2BB;
  dpdbuf4 T2BA;
  dpdbuf4 T2AB;
  global_dpd_->buf4_init(&T2BB, PSIF_CC_TAMPS, 0, 15, 10, 17, 12, 0, "tabij");
  global_dpd_->buf4_init(&T2AB, PSIF_CC_TAMPS, 0, 28, 22, 28, 22, 0, "tAbIj");
  global_dpd_->buf4_init(&T2BA, PSIF_CC_TAMPS, 0, 29, 23, 29, 23, 0, "taBiJ");
  global_dpd_->buf4_init(&FBBints, PSIF_CC_FINTS, 0, 15, 30, 17, 30, 0, "F <bc||ia>");
  global_dpd_->buf4_init(&FBAints, PSIF_CC_FINTS, 0, 29, 27, 29, 27, 0, "F <bC|iA>");
  global_dpd_->buf4_init(&FABints, PSIF_CC_FINTS, 0, 28, 24, 28, 24, 0, "F <Bc|Ia>");
  global_dpd_->buf4_init(&EBBints, PSIF_CC_EINTS, 0, 31, 10, 31, 12, 0,
      "E <ak||ij> (ak, i>j)");
  global_dpd_->buf4_init(&EABints, PSIF_CC_EINTS, 0, 25, 22, 25, 22, 0, "E <aK|Ij>");
  global_dpd_->buf4_init(&EBAints, PSIF_CC_EINTS, 0, 26, 23, 26, 23, 0, "E <Ak|iJ>");
  global_dpd_->buf4_init(&DBBints, PSIF_CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
  global_dpd_->buf4_init(&DBAints, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");

  int Gabc;
  ET = 0.0;
  for (Ga=0; Ga < nirreps; ++Ga) {
    for (a=0; a<bvirtpi[Ga]; ++a) {
      A = bvir_off[Ga] + a;
      for (Gb=0; Gb < nirreps; ++Gb) {
        for (b=0; b<bvirtpi[Gb]; ++b) {
          B = bvir_off[Gb] + b;
          for (Gc=0; Gc < nirreps; ++Gc) {
            for (c=0; c < avirtpi[Gc]; ++c) {
              C = avir_off[Gc] + c;
              Gabc = Ga ^ Gb ^ Gc;
              //Allocate the memory for connected and disconnected triples
              for (Gij=0; Gij < nirreps; ++Gij) {
                Gk = Gij ^ Gabc;
                WijK[Gij] = global_dpd_->dpd_block_matrix(T2BB.params->coltot[Gij],
                    aoccpi[Gk]);
                VijK[Gij] = global_dpd_->dpd_block_matrix(T2BB.params->coltot[Gij],
                    aoccpi[Gk]);
              }

              T3_UHF_AAB_abc(WijK, VijK, 1, nirreps, A, Ga, B, Gb, C, Gc,
                  &T2BB, &T2BA, &T2AB, &FBBints, &FBAints, &FABints, &EBBints,
                  &EBAints, &EABints, &T1B, &T1A, &DBBints, &DBAints, &fia,
                  &fIA, &fij, &fIJ, &fab, &fAB, boccpi, bocc_off, aoccpi,
                  aocc_off, bvirtpi, bvir_off, avirtpi, avir_off, 0.0);

              global_dpd_->file2_mat_init(&fij);
              global_dpd_->file2_mat_init(&fIJ);
              global_dpd_->file2_mat_init(&fab);
              global_dpd_->file2_mat_init(&fAB);
              global_dpd_->file2_mat_rd(&fij);
              global_dpd_->file2_mat_rd(&fIJ);
              global_dpd_->file2_mat_rd(&fab);
              global_dpd_->file2_mat_rd(&fAB);

              double dabc = 0.0;
              if (fab.params->rowtot[Ga])
                dabc -= fab.matrix[Ga][a][a];
              if (fab.params->rowtot[Gb])
                dabc -= fab.matrix[Gb][b][b];
              if (fAB.params->rowtot[Gc])
                dabc -= fAB.matrix[Gc][c][c];

              for (Gij=0; Gij < nirreps; ++Gij) {
                Gk = Gabc ^ Gij;
                for (ij=0; ij<T2BB.params->coltot[Gij]; ++ij) {
                  I = T2BB.params->colorb[Gij][ij][0];
                  Gi = T2BB.params->rsym[I];
                  i = I-bocc_off[Gi];
                  J = T2BB.params->colorb[Gij][ij][1];
                  Gj = T2BB.params->ssym[J];
                  j = J-bocc_off[Gj];
                  for (k=0; k<aoccpi[Gk]; ++k) {
                    denom = dabc;
                    K = k + aocc_off[Gk];
                    if (fij.params->rowtot[Gi])
                      denom += fij.matrix[Gi][i][i];
                    if (fij.params->rowtot[Gj])
                      denom += fij.matrix[Gj][j][j];
                    if (fIJ.params->rowtot[Gk])
                      denom += fIJ.matrix[Gk][k][k];
                    //							  outfile->Printf( "%d,%d,%d,%d,%d,%d = %20.12f\n",I,J,K,A,B,C,VIJk[Gij][ij][k]);


                    ET += WijK[Gij][ij][k] * (WijK[Gij][ij][k]
                        +VijK[Gij][ij][k]) * denom;
                  }
                }
              }
              global_dpd_->file2_mat_close(&fij);
              global_dpd_->file2_mat_close(&fIJ);
              global_dpd_->file2_mat_close(&fab);
              global_dpd_->file2_mat_close(&fAB);

              //Deallocate the memory for connected and disconnected triples
              for (Gij=0; Gij < nirreps; ++Gij) {
                Gk = Gij ^ Gabc;
                global_dpd_->free_dpd_block(WijK[Gij], T2BB.params->coltot[Gij], aoccpi[Gk]);
                global_dpd_->free_dpd_block(VijK[Gij], T2BB.params->coltot[Gij], aoccpi[Gk]);
              }
            }
          }
        }
      }
    }
  }

  outfile->Printf( "    E(T) BBA from Andy's code = %20.15f\n", ET/4.0);
  global_dpd_->file2_close(&fij);
  global_dpd_->file2_close(&fIJ);
  global_dpd_->file2_close(&fab);
  global_dpd_->file2_close(&fAB);

  global_dpd_->buf4_close(&T2BB);
  global_dpd_->buf4_close(&T2AB);
  global_dpd_->buf4_close(&T2BA);
  global_dpd_->buf4_close(&EBBints);
  global_dpd_->buf4_close(&EABints);
  global_dpd_->buf4_close(&EBAints);
  global_dpd_->buf4_close(&FBBints);
  global_dpd_->buf4_close(&FABints);
  global_dpd_->buf4_close(&FBAints);
  global_dpd_->buf4_close(&DBBints);
  global_dpd_->buf4_close(&DBAints);
}

}
} // namespace psi::CCTRIPLES
