/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void x_Gijka_6_uhf(void);
void x_Gijka_7_uhf(void);
void x_Gijka_8_uhf(void);

/* This function computes the non-R0 parts of the 2pdm density matrix
   Gijka = 0.5 *(rho_kaij + rho_Gijka) */

void x_Gijka_uhf(void) { 
  int h, nirreps, i, j, k, a, I, J, K, A, Isym, Jsym, Ksym, Asym, row, col;
  int II,JJ,IIsym,JJsym;
  int L_irr, R_irr, G_irr;
  double value;
  dpdfile2 L1A, T1A, L1B, T1B, R1A, R1B, I1A, I1B;
  dpdbuf4 G, V, T, L, Z, Z1, Z2, Tau;
 
  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* term 1, rho_kaij += Lijae * Rke */
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 2, 20, 2, 20, 0, "L2R1_OOVO(pqsr)");
  dpd_buf4_copy(&Z, EOM_TMP0, "GIJKA");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 12, 30, 12, 30, 0, "L2R1_oovo(pqsr)");
  dpd_buf4_copy(&Z, EOM_TMP0, "Gijka");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 22, 24, 22, 24, 0, "L2R1_OovO(pqsr)");
  dpd_buf4_copy(&Z, EOM_TMP0, "GIjKa");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 23, 27, 23, 27, 0, "L2R1_OoVo(qpsr)");
  dpd_buf4_copy(&Z, EOM_TMP0, "GiJkA");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  dpd_buf4_scm(&Z, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  dpd_buf4_scm(&Z, -1.0);
  dpd_buf4_close(&Z);

  /* term 2, rho_ijka += Rijae * Lke */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 20, 2, 20, 0, "GIJKA");
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 2, 20, 2, 20, 0, "L1R2_OOVO(pqsr)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 12, 30, 12, 30, 0, "Gijka");
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 12, 30, 12, 30, 0, "L1R2_oovo(pqsr)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 22, 24, 22, 24, 0, "L1R2_OovO(pqsr)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 23, 27, 23, 27, 0, "L1R2_OoVo(qpsr)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  /* term 3, rho_ijka += 0.5 Rijef Lkmef tma  */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 20, 2, 20, 0, "GIJKA");
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 2, 0, 2, 2, 0, "R2L2_OOOO");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&Z, &T1A, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&T1A);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 12, 30, 12, 30, 0, "Gijka");
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 12, 10, 12, 12, 0, "R2L2_oooo");
  dpd_file2_init(&T1B, CC_OEI, 0, 2, 3, "tia");
  dpd_contract424(&Z, &T1B, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&T1B);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 22, 22, 22, 22, 0, "R2L2_OoOo");
  dpd_file2_init(&T1B, CC_OEI, 0, 2, 3, "tia");
  dpd_contract424(&Z, &T1B, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&T1B);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 23, 23, 23, 23, 0, "R2L2_oOoO");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&Z, &T1A, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&T1A);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  /* term 4, rho_ijka += 0.5 [Tijef + P(ij) Tie Tjf] Lkmef Rma */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 20, 2, 20, 0, "GIJKA");
  dpd_buf4_init(&Z, EOM_TMP, L_irr, 2, 0, 2, 2, 0, "Tau2L2_OOOO");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&Z, &R1A, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&R1A);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 12, 30, 12, 30, 0, "Gijka");
  dpd_buf4_init(&Z, EOM_TMP, L_irr, 12, 10, 12, 12, 0, "Tau2L2_oooo");
  dpd_file2_init(&R1B, CC_GR, R_irr, 2, 3, "Ria");
  dpd_contract424(&Z, &R1B, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&R1B);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  dpd_buf4_init(&Z, EOM_TMP, L_irr, 22, 22, 22, 22, 0, "Tau2L2_OoOo");
  dpd_file2_init(&R1B, CC_GR, R_irr, 2, 3, "Ria");
  dpd_contract424(&Z, &R1B, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&R1B);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  dpd_buf4_init(&Z, EOM_TMP, L_irr, 23, 23, 23, 23, 0, "Tau2L2_oOoO");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&Z, &R1A, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&R1A);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  
  /* term 5, rho_ijka += - (Lkmef Rmf) (Tijea - P(ij) Tie Tja) */
  if (!params.connect_xi) {
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 20, 2, 20, 0, "GIJKA");
    dpd_buf4_init(&Tau, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tauIJAB");
    dpd_file2_init(&I1A, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_contract244(&I1A, &Tau, &G, 1, 2, 1, -1.0, 1.0);
    dpd_file2_close(&I1A);
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 12, 30, 12, 30, 0, "Gijka");
    dpd_buf4_init(&Tau, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tauijab");
    dpd_file2_init(&I1B, EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    dpd_contract244(&I1B, &Tau, &G, 1, 2, 1, -1.0, 1.0);
    dpd_file2_close(&I1B);
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
    dpd_buf4_init(&Tau, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    dpd_file2_init(&I1A, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    dpd_contract244(&I1A, &Tau, &G, 1, 2, 1, -1.0, 1.0);
    dpd_file2_close(&I1A);
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&G);
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
    dpd_buf4_init(&Tau, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tauiJaB");
    dpd_file2_init(&I1B, EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    dpd_contract244(&I1B, &Tau, &G, 1, 2, 1, -1.0, 1.0);
    dpd_file2_close(&I1B);
    dpd_buf4_close(&Tau);
    dpd_buf4_close(&G);
  }

  x_Gijka_6_uhf();
  x_Gijka_7_uhf();

  /* term 8, +P(ij) Lkmfe rimae tjf */
  /* term 9, +P(ij) Lkmfe Timae Rjf, uses Z3, Z4 */
  x_Gijka_8_uhf();

  /* term 10, +P(IJ) LKMEF RJF TMA TIE */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 2, 0, 2, 0, "Z5(JI,KM)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 2, 20, 2, 20, 0, "L2R1_OOVO(pqsr)");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&V, &T1A, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_buf4_close(&V);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 0, 0, 2, 0, "Z5(JI,KM)");
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 20, 0, 20, 0, "Z5(JI,KA)");
  dpd_contract424(&Z, &T1A, &Z2, 3, 0, 0, 1.0, 0.0);
  dpd_file2_close(&T1A);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 20, 2, 20, 0, "GIJKA");
  dpd_buf4_axpy(&Z2, &G, -1.0);
  dpd_buf4_sort(&Z2, EOM_TMP1, qprs, 0, 20, "Z5(IJ,KA)");
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 0, 20, 0, 20, 0, "Z5(IJ,KA)");
  dpd_buf4_axpy(&Z2, &G, 1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&G);
  /* term 10, +P(ij) lkmef rjf tma tie */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 12, 10, 12, 0, "Z5(ji,km)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 12, 30, 12, 30, 0, "L2R1_oovo(pqsr)");
  dpd_file2_init(&T1B, CC_OEI, 0, 2, 3, "tia");
  dpd_contract424(&V, &T1B, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_buf4_close(&V);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 10, 10, 12, 0, "Z5(ji,km)");
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 30, 10, 30, 0, "Z5(ji,ka)");
  dpd_contract424(&Z, &T1B, &Z2, 3, 0, 0, 1.0, 0.0);
  dpd_file2_close(&T1B);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 30, 12, 30, 0, "Gijka");
  dpd_buf4_axpy(&Z2, &G, -1.0);
  dpd_buf4_sort(&Z2, EOM_TMP1, qprs, 10, 30, "Z5(ij,ka)");
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&Z2, EOM_TMP1, G_irr, 10, 30, 10, 30, 0, "Z5(ij,ka)");
  dpd_buf4_axpy(&Z2, &G, 1.0);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&G);
  /* term 10, GIjKa += LKmEf Rjf TIE tma */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 22, 22, 22, 22, 0, "Z(Ij,Km)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 22, 26, 22, 26, 0, "L2R1_OoVo");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1A, &V, &Z, 1, 2, 0, 1.0, 0.0);
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_init(&V, EOM_TMP, G_irr, 22, 24, 22, 24, 0, "L2R1_OovO(pqsr)");
  dpd_file2_init(&T1B, CC_OEI, 0, 2, 3, "tia");
  dpd_contract424(&V, &T1B, &Z, 3, 1, 1, 1.0, 1.0);
  dpd_file2_close(&T1B);
  dpd_buf4_close(&V);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  dpd_file2_init(&T1B, CC_OEI, 0, 2, 3, "tia");
  dpd_contract424(&Z, &T1B, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&T1B);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  /* term 10, GiJkA += P(ij) LkMeF RJF Tie tMA */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 23, 23, 23, 23, 0, "Z(iJ,kM)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 23, 25, 23, 25, 0, "L2R1_OovO(qprs)");
  dpd_file2_init(&T1B, CC_OEI, 0, 2, 3, "tia");
  dpd_contract244(&T1B, &V, &Z, 1, 2, 0, 1.0, 0.0);
  dpd_file2_close(&T1B);
  dpd_buf4_close(&V);
  dpd_buf4_init(&V, EOM_TMP, G_irr, 23, 27, 23, 27, 0, "L2R1_OoVo(qpsr)");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&V, &T1A, &Z, 3, 1, 1, 1.0, 1.0);
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&Z, &T1A, &G, 3, 0, 0, 1.0, 1.0);
  dpd_file2_close(&T1A);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  psio_close(EOM_TMP1, 0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);

  /* add 1/2 to ground-state parts of density */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 2, 20, 2, 20, 0, "GIJKA");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 2, 20, 2, 20, 0, "GIJKA");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 12, 30, 12, 30, 0, "Gijka");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 12, 30, 12, 30, 0, "Gijka");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);

  /* clear out temporary files */
  psio_close(EOM_TMP0, 0);
  psio_open(EOM_TMP0, PSIO_OPEN_NEW);

  return;
}



/* This function computes term 6,
   rho_ijka -= P(ij) Lke Rie Tja + P(ij) Lmkef Rmief Tja or
   rho_ijka -= P(ij) LR_OO(k,i) T(j,a) */

void x_Gijka_6_uhf(void) { 
  int h, nirreps, i, j, k, a, I, J, K, A, Isym, Jsym, Ksym, Asym, row, col;
  int II,JJ,IIsym,JJsym;
  int L_irr, R_irr, G_irr;
  dpdfile2 LRA, LRB, T1A, T1B;
  dpdbuf4 G;
 
  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* open one-electron files for the nasty terms */
  dpd_file2_init(&LRA, EOM_TMP, G_irr, 0, 0, "LR_OO");
  dpd_file2_init(&LRB, EOM_TMP, G_irr, 2, 2, "LR_oo");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&T1B, CC_OEI, 0, 2, 3, "tia");
  dpd_file2_mat_init(&T1A);   dpd_file2_mat_init(&T1B);
  dpd_file2_mat_init(&LRA);  dpd_file2_mat_init(&LRB);
  dpd_file2_mat_rd(&T1A);     dpd_file2_mat_rd(&T1B);
  dpd_file2_mat_rd(&LRA);    dpd_file2_mat_rd(&LRB);

  /* rho_IJKA += - LR_OO(K,I) T(J,A) + LR_OO(K,J) T(I,A) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 20, 2, 20, 0, "GIJKA");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LRA.params->colidx[i]; Isym = LRA.params->qsym[i];
      J = T1A.params->rowidx[j]; Jsym = T1A.params->psym[j];
      II = T1A.params->rowidx[i]; IIsym = T1A.params->psym[i];
      JJ = LRA.params->colidx[j]; JJsym = LRA.params->qsym[j];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        k = G.params->colorb[h^G_irr][col][0];
        a = G.params->colorb[h^G_irr][col][1];
        K = LRA.params->rowidx[k]; Ksym = LRA.params->psym[k];
        A = T1A.params->colidx[a]; Asym = T1A.params->qsym[a];
        if( ((Ksym^Isym)==G_irr) && (Jsym==Asym))
          G.matrix[h][row][col] -=
            LRA.matrix[Ksym][K][I] * T1A.matrix[Jsym][J][A];
        if( ((Ksym^JJsym)==G_irr) && (IIsym==Asym))
          G.matrix[h][row][col] +=
            LRA.matrix[Ksym][K][JJ] * T1A.matrix[IIsym][II][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  /* rho_ijka += - LR_oo(k,i) T(j,a) + LR_oo(k,j) T(i,a) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 30, 12, 30, 0, "Gijka");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LRB.params->colidx[i]; Isym = LRB.params->qsym[i];
      J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
      II = T1B.params->rowidx[i]; IIsym = T1B.params->psym[i];
      JJ = LRB.params->colidx[j]; JJsym = LRB.params->qsym[j];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        k = G.params->colorb[h^G_irr][col][0];
        a = G.params->colorb[h^G_irr][col][1];
        K = LRB.params->rowidx[k]; Ksym = LRB.params->psym[k];
        A = T1B.params->colidx[a]; Asym = T1B.params->qsym[a];
        if( ((Ksym^Isym)==G_irr) && (Jsym==Asym))
          G.matrix[h][row][col] -=
            LRB.matrix[Ksym][K][I] * T1B.matrix[Jsym][J][A];
        if( ((Ksym^JJsym)==G_irr) && (IIsym==Asym))
          G.matrix[h][row][col] +=
            LRB.matrix[Ksym][K][JJ] * T1B.matrix[IIsym][II][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  /* rho_IjKa += - LR_OO(K,I) T(j,a) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LRA.params->colidx[i]; Isym = LRA.params->qsym[i];
      J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        k = G.params->colorb[h^G_irr][col][0];
        a = G.params->colorb[h^G_irr][col][1];
        K = LRA.params->rowidx[k]; Ksym = LRA.params->psym[k];
        A = T1B.params->colidx[a]; Asym = T1B.params->qsym[a];
        if( ((Ksym^Isym)==G_irr) && (Jsym==Asym))
          G.matrix[h][row][col] -=
            LRA.matrix[Ksym][K][I] * T1B.matrix[Jsym][J][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  /* rho_iJkA += - LR_oo(k,i) T(J,A) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LRB.params->colidx[i]; Isym = LRB.params->qsym[i];
      J = T1A.params->rowidx[j]; Jsym = T1A.params->psym[j];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        k = G.params->colorb[h^G_irr][col][0];
        a = G.params->colorb[h^G_irr][col][1];
        K = LRB.params->rowidx[k]; Ksym = LRB.params->psym[k];
        A = T1A.params->colidx[a]; Asym = T1A.params->qsym[a];
        if( ((Ksym^Isym)==G_irr) && (Jsym==Asym))
          G.matrix[h][row][col] -=
            LRB.matrix[Ksym][K][I] * T1A.matrix[Jsym][J][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);
  dpd_file2_mat_close(&LRA);
  dpd_file2_mat_close(&LRB);
  dpd_file2_close(&LRA);
  dpd_file2_close(&LRB);

  dpd_file2_mat_close(&T1A);
  dpd_file2_mat_close(&T1B);
  dpd_file2_close(&T1A);
  dpd_file2_close(&T1B);

  return;
}



/* This function computes 
   Gijka -= P(ij) Lke Tie Rja + Lmkef Tmief Rja or
         -= P(ij) LT_OO(k,i) * R(j,a) */
void x_Gijka_7_uhf(void) { 
  int h, nirreps, i, j, k, a, I, J, K, A, Isym, Jsym, Ksym, Asym, row, col;
  int II,JJ,IIsym,JJsym;
  int L_irr, R_irr, G_irr;
  dpdfile2 R1A, R1B, LTA, LTB;
  dpdbuf4 G;
 
  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* open one-electron files for the nasty terms */
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_file2_init(&R1B, CC_GR, R_irr, 2, 3, "Ria");
  dpd_file2_init(&LTA, EOM_TMP, L_irr, 0, 0, "LT_OO");
  dpd_file2_init(&LTB, EOM_TMP, L_irr, 2, 2, "LT_oo");
  dpd_file2_mat_init(&R1A);
  dpd_file2_mat_init(&R1B);
  dpd_file2_mat_init(&LTA);
  dpd_file2_mat_init(&LTB);
  dpd_file2_mat_rd(&R1A);
  dpd_file2_mat_rd(&R1B);
  dpd_file2_mat_rd(&LTA);
  dpd_file2_mat_rd(&LTB);

  /* rho_IJKA += - LT_OO(K,I)   R(J,A) + LT_OO(K,J)   R(I,A) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 20, 2, 20, 0, "GIJKA");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LTA.params->colidx[i];  Isym = LTA.params->qsym[i];
      J = R1A.params->rowidx[j];  Jsym = R1A.params->psym[j];
      II = R1A.params->rowidx[i]; IIsym = R1A.params->psym[i];
      JJ = LTA.params->colidx[j]; JJsym = LTA.params->qsym[j];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        k = G.params->colorb[h^G_irr][col][0];
        a = G.params->colorb[h^G_irr][col][1];
        K = LTA.params->rowidx[k]; Ksym = LTA.params->psym[k];
        A = R1A.params->colidx[a]; Asym = R1A.params->qsym[a];
        if( ((Ksym^Isym)==L_irr) && ((Jsym^Asym)==R_irr) )
          G.matrix[h][row][col] -=
            LTA.matrix[Ksym][K][I] * R1A.matrix[Jsym][J][A];
        if( ((Ksym^JJsym)==L_irr) && ((IIsym^Asym)==R_irr) )
          G.matrix[h][row][col] +=
            LTA.matrix[Ksym][K][JJ] * R1A.matrix[IIsym][II][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  /* rho_ijka += - LT_oo(k,i)   R(j,a) + LT_oo(k,j)   R(i,a) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 30, 12, 30, 0, "Gijka");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LTB.params->colidx[i];  Isym = LTB.params->qsym[i];
      J = R1B.params->rowidx[j];  Jsym = R1B.params->psym[j];
      II = R1B.params->rowidx[i]; IIsym = R1B.params->psym[i];
      JJ = LTB.params->colidx[j]; JJsym = LTB.params->qsym[j];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        k = G.params->colorb[h^G_irr][col][0];
        a = G.params->colorb[h^G_irr][col][1];
        K = LTB.params->rowidx[k]; Ksym = LTB.params->psym[k];
        A = R1B.params->colidx[a]; Asym = R1B.params->qsym[a];
        if( ((Ksym^Isym)==L_irr) && ((Jsym^Asym)==R_irr) )
          G.matrix[h][row][col] -=
            LTB.matrix[Ksym][K][I] * R1B.matrix[Jsym][J][A];
        if( ((Ksym^JJsym)==L_irr) && ((IIsym^Asym)==R_irr) )
          G.matrix[h][row][col] +=
            LTB.matrix[Ksym][K][JJ] * R1B.matrix[IIsym][II][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  /* rho_IjKa += - LT_OO(K,I)   R(j,a) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LTA.params->colidx[i];  Isym = LTA.params->qsym[i];
      J = R1B.params->rowidx[j];  Jsym = R1B.params->psym[j];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        k = G.params->colorb[h^G_irr][col][0];
        a = G.params->colorb[h^G_irr][col][1];
        K = LTA.params->rowidx[k]; Ksym = LTA.params->psym[k];
        A = R1B.params->colidx[a]; Asym = R1B.params->qsym[a];
        if( ((Ksym^Isym)==L_irr) && ((Jsym^Asym)==R_irr) )
          G.matrix[h][row][col] -=
            LTA.matrix[Ksym][K][I] * R1B.matrix[Jsym][J][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  /* rho_iJkA += - LT_oo(k,i)   R(J,A) */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LTB.params->colidx[i];  Isym = LTB.params->qsym[i];
      J = R1A.params->rowidx[j];  Jsym = R1A.params->psym[j];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        k = G.params->colorb[h^G_irr][col][0];
        a = G.params->colorb[h^G_irr][col][1];
        K = LTB.params->rowidx[k]; Ksym = LTB.params->psym[k];
        A = R1A.params->colidx[a]; Asym = R1A.params->qsym[a];
        if( ((Ksym^Isym)==L_irr) && ((Jsym^Asym)==R_irr) )
          G.matrix[h][row][col] -=
            LTB.matrix[Ksym][K][I] * R1A.matrix[Jsym][J][A];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  dpd_file2_mat_close(&LTA);
  dpd_file2_mat_close(&LTB);
  dpd_file2_mat_close(&R1A);
  dpd_file2_mat_close(&R1B);

  dpd_file2_close(&LTA);
  dpd_file2_close(&LTB);
  dpd_file2_close(&R1A);
  dpd_file2_close(&R1B);

  return;
}



/* This function computes term 8 and term 9 of Gijka */
/* term 8, +P(ij) Lkmfe Rimae Tjf */
/* term 9, +P(ij) Lkmfe Timae Rjf */

void x_Gijka_8_uhf(void) { 
  int h, nirreps, i, j, k, a, I, J, K, A, Isym, Jsym, Ksym, Asym, row, col;
  int II,JJ,IIsym,JJsym;
  int L_irr, R_irr, G_irr;
  double value;
  dpdfile2 L1A, T1A, L1B, T1B, R1A, R1B, I1A, I1B;
  dpdbuf4 G, V, T, L, Z, Z1, Z2, Tau;

  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* term 8, +P(ij) Lkmfe Rimae Tjf */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 20, 0, 20, 0, 0, "Z(IA,KJ)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 20, 20, 20, 20, 0, "R2L2_OVOV");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, psrq, 0, 20, "Z(IJ,KA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 20, 0, 20, 0, "Z(IJ,KA)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 20, 2, 20, 0, "GIJKA");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 0, 20, "Z(JI,KA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 20, 0, 20, 0, "Z(JI,KA)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 30, 10, 30, 10, 0, "Z(ia,kj)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 30, 30, 30, 30, 0, "R2L2_ovov");
  dpd_file2_init(&T1A, CC_OEI, 0, 2, 3, "tia");
  dpd_contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, psrq, 10, 30, "Z(ij,ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 30, 10, 30, 0, "Z(ij,ka)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 30, 12, 30, 0, "Gijka");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 10, 30, "Z(ji,ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 30, 10, 30, 0, "Z(ji,ka)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  /* GIjKa += R2L2_OvOv(Ia,Kf) T(j,f) */
  /* GIjKa -= R2L2_OvOv(ja,KF) T(I,F) */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 24, 22, 24, 22, 0, "Z(Ia,Kj)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 24, 24, 24, 24, 0, "R2L2_OvOv");
  dpd_file2_init(&T1A, CC_OEI, 0, 2, 3, "tia");
  dpd_contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, psrq, 22, 24, "Z(Ij,Ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 22, 24, 22, 24, 0, "Z(Ij,Ka)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 30, 0, 30, 0, 0, "Z2(ja,KI)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 30, 20, 30, 20, 0, "R2L2_ovOV");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, psrq, 23, 24, "Z2(jI,Ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 23, 24, 23, 24, 0, "Z2(jI,Ka)");
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 22, 24, "Z2(Ij,Ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 22, 24, 22, 24, 0, "Z2(Ij,Ka)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  /* GiJkA += R2L2_oVoV(iA,kF) T(J,F) */
  /* GiJkA += R2L2_OVov(JA,kf) T(i,f) */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 27, 23, 27, 23, 0, "Z(iA,kJ)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 27, 27, 27, 27, 0, "R2L2_oVoV");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, psrq, 23, 27, "Z(iJ,kA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 23, 27, 23, 27, 0, "Z(iJ,kA)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 20, 10, 20, 10, 0, "Z2(JA,ki)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 20, 30, 20, 30, 0, "R2L2_OVov");
  dpd_file2_init(&T1A, CC_OEI, 0, 2, 3, "tia");
  dpd_contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, psrq, 22, 27, "Z2(Ji,kA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 22, 27, 22, 27, 0, "Z2(Ji,kA)");
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 23, 27, "Z2(iJ,kA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 23, 27, 23, 27, 0, "Z2(iJ,kA)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  psio_close(EOM_TMP1, 0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);

  /* term 9, +P(ij) Lkmfe Timae Rjf */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 20, 0, 20, 0, 0, "Z3(IA,KJ)");
  dpd_buf4_init(&V, EOM_TMP, L_irr, 20, 20, 20, 20, 0, "VIAJB");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&V, &R1A, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&R1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, psrq, 0, 20, "Z3(IJ,KA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 20, 0, 20, 0, "Z3(IJ,KA)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 0, 20, 2, 20, 0, "GIJKA");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 0, 20, "Z3(JI,KA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 0, 20, 0, 20, 0, "Z3(JI,KA)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 30, 10, 30, 10, 0, "Z3(ia,kj)");
  dpd_buf4_init(&V, EOM_TMP, L_irr, 30, 30, 30, 30, 0, "Viajb");
  dpd_file2_init(&R1B, CC_GR, R_irr, 2, 3, "Ria");
  dpd_contract424(&V, &R1B, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&R1B);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, psrq, 10, 30, "Z3(ij,ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 30, 10, 30, 0, "Z3(ij,ka)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 10, 30, 12, 30, 0, "Gijka");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 10, 30, "Z3(ji,ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 30, 10, 30, 0, "Z3(ji,ka)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  /* GIjKa += L2T2_OvOv(Ia,Kf) R(j,f) */
  /* GIjKa -= L2T2_OvOv(ja,KF) R(I,F) */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 24, 22, 24, 22, 0, "Z3(Ia,Kj)");
  dpd_buf4_init(&V, EOM_TMP, L_irr, 24, 24, 24, 24, 0, "VIaJb");
  dpd_file2_init(&R1B, CC_GR, R_irr, 2, 3, "Ria");
  dpd_contract424(&V, &R1B, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&R1B);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, psrq, 22, 24, "Z3(Ij,Ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 22, 24, 22, 24, 0, "Z3(Ij,Ka)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 30, 0, 30, 0, 0, "Z4(ja,KI)");
  dpd_buf4_init(&V, EOM_TMP, L_irr, 30, 20, 30, 20, 0, "ViaJB");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&V, &R1A, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&R1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, psrq, 23, 24, "Z4(jI,Ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 23, 24, 23, 24, 0, "Z4(jI,Ka)");
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 22, 24, "Z4(Ij,Ka)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 22, 24, 22, 24, 0, "Z4(Ij,Ka)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  /* GiJkA += L2T2_oVoV(iA,kF) R(J,F) */
  /* GiJkA += L2T2_OVov(JA,kf) R(i,f) */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 27, 23, 27, 23, 0, "Z3(iA,kJ)");
  dpd_buf4_init(&V, EOM_TMP, L_irr, 27, 27, 27, 27, 0, "ViAjB");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&V, &R1A, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&R1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, psrq, 23, 27, "Z3(iJ,kA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 23, 27, 23, 27, 0, "Z3(iJ,kA)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 20, 10, 20, 10, 0, "Z4(JA,ki)");
  dpd_buf4_init(&V, EOM_TMP, L_irr, 20, 30, 20, 30, 0, "VIAjb");
  dpd_file2_init(&R1B, CC_GR, R_irr, 2, 3, "Ria");
  dpd_contract424(&V, &R1B, &Z, 3, 1, 0, 1.0, 0.0); 
  dpd_file2_close(&R1B);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, psrq, 22, 27, "Z4(Ji,kA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 22, 27, 22, 27, 0, "Z4(Ji,kA)");
  dpd_buf4_sort(&Z, EOM_TMP1, qprs, 23, 27, "Z4(iJ,kA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 23, 27, 23, 27, 0, "Z4(iJ,kA)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  return;
}

}} // namespace psi::ccdensity
