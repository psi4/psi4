#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <psifiles.h>
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

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");

  dpdbuf4 EAAints;
  dpd_buf4_init(&EAAints, CC_EINTS, 0, 21, 0, 21, 2, 0, "E <AK||IJ> (AK, I>J)");
  dpdbuf4 FAAints;
  dpd_buf4_init(&FAAints, CC_FINTS, 0, 5, 20, 7, 20, 0, "F <BC||IA>");
  dpdbuf4 T2AA;
  dpd_buf4_init(&T2AA, CC_TAMPS, 0, 5, 0, 7, 2, 0, "tABIJ");
  dpdbuf4 DAAints;
  dpd_buf4_init(&DAAints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");

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
                WIJK[Gij] = dpd_block_matrix(T2AA.params->coltot[Gij],
                    occpi[Gk]);
                VIJK[Gij] = dpd_block_matrix(T2AA.params->coltot[Gij],
                    occpi[Gk]);
              }

              T3_UHF_AAA_abc(WIJK, VIJK, 1, nirreps, A, Ga, B, Gb, C, Gc,
                  &T2AA, &FAAints, &EAAints, &T1, &DAAints, &fIA, &fIJ, &fAB,
                  occpi, occ_off, virtpi, vir_off, 0.0);
              dpd_file2_mat_init(&fIJ);
              dpd_file2_mat_init(&fIA);
              dpd_file2_mat_init(&fAB);
              dpd_file2_mat_init(&T1);
              dpd_file2_mat_rd(&fIJ);
              dpd_file2_mat_rd(&fIA);
              dpd_file2_mat_rd(&fAB);
              dpd_file2_mat_rd(&T1);
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
              dpd_file2_mat_close(&fIJ);
              dpd_file2_mat_close(&fAB);
              dpd_file2_mat_close(&fIA);
              dpd_file2_mat_close(&T1);

              //Deallocate the memory for connected and disconnected triples
              for (Gij=0; Gij < nirreps; ++Gij) {
                Gk = Gij ^ Gabc;
                dpd_free_block(WIJK[Gij], T2AA.params->coltot[Gij], occpi[Gk]);
                dpd_free_block(VIJK[Gij], T2AA.params->coltot[Gij], occpi[Gk]);
              }
            }
          }
        }
      }
    }
  }
  fprintf(outfile, "\tE(T) AAA from Andy's code = %20.15f\n", ET/36.0);
  dpd_buf4_close(&EAAints);
  dpd_buf4_close(&FAAints);
  dpd_buf4_close(&T2AA);
  dpd_file2_close(&T1);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fIA);
  dpd_file2_close(&fAB);
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

  dpd_file2_init(&fij, CC_OEI, 0, 2, 2, "fij");
  dpd_file2_init(&fab, CC_OEI, 0, 3, 3, "fab");
  dpd_file2_init(&fia, CC_OEI, 0, 2, 3, "fia");
  dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");

  dpdbuf4 EBBints;
  dpd_buf4_init(&EBBints, CC_EINTS, 0, 31, 10, 31, 12, 0,
      "E <ak||ij> (ak, i>j)");
  dpdbuf4 FBBints;
  dpd_buf4_init(&FBBints, CC_FINTS, 0, 15, 30, 17, 30, 0, "F <bc||ia>");
  dpdbuf4 T2BB;
  dpd_buf4_init(&T2BB, CC_TAMPS, 0, 15, 10, 17, 12, 0, "tabij");
  dpdbuf4 DBBints;
  dpd_buf4_init(&DBBints, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");

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
                WIJK[Gij] = dpd_block_matrix(T2BB.params->coltot[Gij],
                    occpi[Gk]);
                VIJK[Gij] = dpd_block_matrix(T2BB.params->coltot[Gij],
                    occpi[Gk]);
              }

              T3_UHF_AAA_abc(WIJK, VIJK, 1, nirreps, A, Ga, B, Gb, C, Gc,
                  &T2BB, &FBBints, &EBBints, &T1, &DBBints, &fia, &fij, &fab,
                  occpi, occ_off, virtpi, vir_off, 0.0);
              dpd_file2_mat_init(&fij);
              dpd_file2_mat_init(&fia);
              dpd_file2_mat_init(&fab);
              dpd_file2_mat_init(&T1);
              dpd_file2_mat_rd(&fij);
              dpd_file2_mat_rd(&fia);
              dpd_file2_mat_rd(&fab);
              dpd_file2_mat_rd(&T1);
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
              dpd_file2_mat_close(&fij);
              dpd_file2_mat_close(&fab);
              dpd_file2_mat_close(&fia);
              dpd_file2_mat_close(&T1);

              //Deallocate the memory for connected and disconnected triples
              for (Gij=0; Gij < nirreps; ++Gij) {
                Gk = Gij ^ Gabc;
                dpd_free_block(WIJK[Gij], T2BB.params->coltot[Gij], occpi[Gk]);
                dpd_free_block(VIJK[Gij], T2BB.params->coltot[Gij], occpi[Gk]);
              }
            }
          }
        }
      }
    }
  }
  fprintf(outfile, "\tE(T) BBB from Andy's code = %20.15f\n", ET/36.0);
  dpd_buf4_close(&EBBints);
  dpd_buf4_close(&FBBints);
  dpd_buf4_close(&T2BB);
  dpd_file2_close(&T1);
  dpd_file2_close(&fij);
  dpd_file2_close(&fia);
  dpd_file2_close(&fab);
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

  dpd_file2_init(&fij, CC_OEI, 0, 2, 2, "fij");
  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fab, CC_OEI, 0, 3, 3, "fab");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_init(&fia, CC_OEI, 0, 2, 3, "fia");
  dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
  dpd_file2_init(&T1B, CC_OEI, 0, 2, 3, "tia");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");

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
  dpd_buf4_init(&T2AA, CC_TAMPS, 0, 5, 0, 7, 2, 0, "tABIJ");
  dpd_buf4_init(&T2AB, CC_TAMPS, 0, 28, 22, 28, 22, 0, "tAbIj");
  dpd_buf4_init(&T2BA, CC_TAMPS, 0, 29, 23, 29, 23, 0, "taBiJ");
  dpd_buf4_init(&FAAints, CC_FINTS, 0, 5, 20, 7, 20, 0, "F <BC||IA>");
  dpd_buf4_init(&FBAints, CC_FINTS, 0, 29, 27, 29, 27, 0, "F <bC|iA>");
  dpd_buf4_init(&FABints, CC_FINTS, 0, 28, 24, 28, 24, 0, "F <Bc|Ia>");
  dpd_buf4_init(&EAAints, CC_EINTS, 0, 21, 0, 21, 2, 0, "E <AK||IJ> (AK, I>J)");
  dpd_buf4_init(&EABints, CC_EINTS, 0, 25, 22, 25, 22, 0, "E <aK|Ij>");
  dpd_buf4_init(&EBAints, CC_EINTS, 0, 26, 23, 26, 23, 0, "E <Ak|iJ>");
  dpd_buf4_init(&DAAints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
  dpd_buf4_init(&DABints, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");

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
                WIJk[Gij] = dpd_block_matrix(T2AA.params->coltot[Gij],
                    boccpi[Gk]);
                VIJk[Gij] = dpd_block_matrix(T2AA.params->coltot[Gij],
                    boccpi[Gk]);
              }
              T3_UHF_AAB_abc(WIJk, VIJk, 1, nirreps, A, Ga, B, Gb, C, Gc,
                  &T2AA, &T2AB, &T2BA, &FAAints, &FABints, &FBAints, &EAAints,
                  &EABints, &EBAints, &T1A, &T1B, &DAAints, &DABints, &fIA,
                  &fia, &fIJ, &fij, &fAB, &fab, aoccpi, aocc_off, boccpi,
                  bocc_off, avirtpi, avir_off, bvirtpi, bvir_off, 0.0);
              dpd_file2_mat_init(&fij);
              dpd_file2_mat_init(&fIJ);
              dpd_file2_mat_init(&fab);
              dpd_file2_mat_init(&fAB);
              dpd_file2_mat_rd(&fij);
              dpd_file2_mat_rd(&fIJ);
              dpd_file2_mat_rd(&fab);
              dpd_file2_mat_rd(&fAB);

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
                    //							  fprintf(outfile, "%d,%d,%d,%d,%d,%d = %20.12f\n",I,J,K,A,B,C,VIJk[Gij][ij][k]);


                    ET += WIJk[Gij][ij][k] * (WIJk[Gij][ij][k]
                        +VIJk[Gij][ij][k]) * denom;
                  }
                }
              }
              dpd_file2_mat_close(&fij);
              dpd_file2_mat_close(&fIJ);
              dpd_file2_mat_close(&fab);
              dpd_file2_mat_close(&fAB);

              //Deallocate the memory for connected and disconnected triples
              for (Gij=0; Gij < nirreps; ++Gij) {
                Gk = Gij ^ Gabc;
                dpd_free_block(WIJk[Gij], T2AA.params->coltot[Gij], boccpi[Gk]);
                dpd_free_block(VIJk[Gij], T2AA.params->coltot[Gij], boccpi[Gk]);
              }
            }
          }
        }
      }
    }
  }

  fprintf(outfile, "\tE(T) AAB from Andy's code = %20.15f\n", ET/4.0);
  dpd_file2_close(&fij);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fab);
  dpd_file2_close(&fAB);

  dpd_buf4_close(&T2AA);
  dpd_buf4_close(&T2AB);
  dpd_buf4_close(&T2BA);
  dpd_buf4_close(&EAAints);
  dpd_buf4_close(&EABints);
  dpd_buf4_close(&EBAints);
  dpd_buf4_close(&FAAints);
  dpd_buf4_close(&FABints);
  dpd_buf4_close(&FBAints);
  dpd_buf4_close(&DAAints);
  dpd_buf4_close(&DABints);
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

  dpd_file2_init(&fij, CC_OEI, 0, 2, 2, "fij");
  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fab, CC_OEI, 0, 3, 3, "fab");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_init(&fia, CC_OEI, 0, 2, 3, "fia");
  dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
  dpd_file2_init(&T1B, CC_OEI, 0, 2, 3, "tia");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");

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
  dpd_buf4_init(&T2BB, CC_TAMPS, 0, 15, 10, 17, 12, 0, "tabij");
  dpd_buf4_init(&T2AB, CC_TAMPS, 0, 28, 22, 28, 22, 0, "tAbIj");
  dpd_buf4_init(&T2BA, CC_TAMPS, 0, 29, 23, 29, 23, 0, "taBiJ");
  dpd_buf4_init(&FBBints, CC_FINTS, 0, 15, 30, 17, 30, 0, "F <bc||ia>");
  dpd_buf4_init(&FBAints, CC_FINTS, 0, 29, 27, 29, 27, 0, "F <bC|iA>");
  dpd_buf4_init(&FABints, CC_FINTS, 0, 28, 24, 28, 24, 0, "F <Bc|Ia>");
  dpd_buf4_init(&EBBints, CC_EINTS, 0, 31, 10, 31, 12, 0,
      "E <ak||ij> (ak, i>j)");
  dpd_buf4_init(&EABints, CC_EINTS, 0, 25, 22, 25, 22, 0, "E <aK|Ij>");
  dpd_buf4_init(&EBAints, CC_EINTS, 0, 26, 23, 26, 23, 0, "E <Ak|iJ>");
  dpd_buf4_init(&DBBints, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
  dpd_buf4_init(&DBAints, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");

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
                WijK[Gij] = dpd_block_matrix(T2BB.params->coltot[Gij],
                    aoccpi[Gk]);
                VijK[Gij] = dpd_block_matrix(T2BB.params->coltot[Gij],
                    aoccpi[Gk]);
              }

              T3_UHF_AAB_abc(WijK, VijK, 1, nirreps, A, Ga, B, Gb, C, Gc,
                  &T2BB, &T2BA, &T2AB, &FBBints, &FBAints, &FABints, &EBBints,
                  &EBAints, &EABints, &T1B, &T1A, &DBBints, &DBAints, &fia,
                  &fIA, &fij, &fIJ, &fab, &fAB, boccpi, bocc_off, aoccpi,
                  aocc_off, bvirtpi, bvir_off, avirtpi, avir_off, 0.0);

              dpd_file2_mat_init(&fij);
              dpd_file2_mat_init(&fIJ);
              dpd_file2_mat_init(&fab);
              dpd_file2_mat_init(&fAB);
              dpd_file2_mat_rd(&fij);
              dpd_file2_mat_rd(&fIJ);
              dpd_file2_mat_rd(&fab);
              dpd_file2_mat_rd(&fAB);

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
                    //							  fprintf(outfile, "%d,%d,%d,%d,%d,%d = %20.12f\n",I,J,K,A,B,C,VIJk[Gij][ij][k]);


                    ET += WijK[Gij][ij][k] * (WijK[Gij][ij][k]
                        +VijK[Gij][ij][k]) * denom;
                  }
                }
              }
              dpd_file2_mat_close(&fij);
              dpd_file2_mat_close(&fIJ);
              dpd_file2_mat_close(&fab);
              dpd_file2_mat_close(&fAB);

              //Deallocate the memory for connected and disconnected triples
              for (Gij=0; Gij < nirreps; ++Gij) {
                Gk = Gij ^ Gabc;
                dpd_free_block(WijK[Gij], T2BB.params->coltot[Gij], aoccpi[Gk]);
                dpd_free_block(VijK[Gij], T2BB.params->coltot[Gij], aoccpi[Gk]);
              }
            }
          }
        }
      }
    }
  }

  fprintf(outfile, "\tE(T) BBA from Andy's code = %20.15f\n", ET/4.0);
  dpd_file2_close(&fij);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fab);
  dpd_file2_close(&fAB);

  dpd_buf4_close(&T2BB);
  dpd_buf4_close(&T2AB);
  dpd_buf4_close(&T2BA);
  dpd_buf4_close(&EBBints);
  dpd_buf4_close(&EABints);
  dpd_buf4_close(&EBAints);
  dpd_buf4_close(&FBBints);
  dpd_buf4_close(&FABints);
  dpd_buf4_close(&FBAints);
  dpd_buf4_close(&DBBints);
  dpd_buf4_close(&DBAints);
}

}
} // namespace psi::CCTRIPLES
