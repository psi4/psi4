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

void x_Gciab_rohf(void);
void x_Gciab_6(void);
void x_Gciab_7(void);
void x_Gciab_8_rohf(void);
extern void x_Gciab_uhf(void);

/* This function computes the non-R0 parts of the 2pdm density matrix
   Gciab = 0.5 *(rho_abci + rho_ciab) */
void x_Gciab(void) {
  if (params.ref == 0 || params.ref == 1)
    x_Gciab_rohf();
  else
    x_Gciab_uhf();
}

void x_Gciab_rohf(void) { 
  int h, nirreps, i, j, k, a, I, J, K, A, Isym, Jsym, Ksym, Asym, row, col;
  int II,JJ,IIsym,JJsym;
  int L_irr, R_irr, G_irr;
  double value;
  dpdfile2 L1, T1, R1, I1;
  dpdbuf4 G, V, T, L, Z, Z2, R, Tau;
 
  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* term 1, rho_abci += Lmiab * Rmc */
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 7, 11, 7, 11, 0, "L2R1_VVOV(pqsr)");
  dpd_buf4_sort(&Z, EOM_TMP0, rspq, 11, 7, "GCIAB");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 7, 11, 7, 11, 0, "L2R1_vvov(pqsr)");
  dpd_buf4_sort(&Z, EOM_TMP0, rspq, 11, 7, "Gciab");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 5, 11, 5, 11, 0, "L2R1_VvoV(pqsr)");
  dpd_buf4_sort(&Z, EOM_TMP0, rspq, 11, 5, "GCiAb");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 5, 11, 5, 11, 0, "L2R1_VvOv(pqsr)");
  dpd_buf4_sort(&Z, EOM_TMP0, rsqp, 11, 5, "GcIaB");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  dpd_buf4_scm(&Z, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "Gciab");
  dpd_buf4_scm(&Z, -1.0);
  dpd_buf4_close(&Z);

  /* term 2, rho_ciab += Rmiab * Lmc */
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 7, 11, 7, 11, 0, "R2L1_VVOV(pqsr)");
  dpd_buf4_sort_axpy(&Z, EOM_TMP0, rspq, 11, 7, "GCIAB",-1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 7, 11, 7, 11, 0, "R2L1_vvov(pqsr)");
  dpd_buf4_sort_axpy(&Z, EOM_TMP0, rspq, 11, 7, "Gciab",-1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 5, 11, 5, 11, 0, "R2L1_VvoV(pqsr)");
  dpd_buf4_sort_axpy(&Z, EOM_TMP0, rspq, 11, 5, "GCiAb", 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 5, 11, 5, 11, 0, "R2L1_VvOv(pqsr)");
  dpd_buf4_sort_axpy(&Z, EOM_TMP0, rsqp, 11, 5, "GcIaB", 1.0);
  dpd_buf4_close(&Z);
  /* same in two-steps, perhaps suggestive of future out-of-core approach
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "Gciab");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 7, 11, 7, 0, "Gciab");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);
  */

  /* term 3, rho_CIAB -= 0.5 LMNCE RMNAB tIE */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 2, 11, 2, 0, "Z(CI,MN)");
  dpd_buf4_init(&L, CC_GL, L_irr, 2, 5, 2, 7, 0, "LIJAB");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&L, &T1, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&L);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  dpd_buf4_init(&R, CC_GR, R_irr, 2, 7, 2, 7, 0, "RIJAB");
  dpd_contract444(&Z, &R, &G, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&R);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 3, rho_ciab -= 0.5 Lmnce Rmnab tie */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 2, 11, 2, 0, "Z(ci,mn)");
  dpd_buf4_init(&L, CC_GL, L_irr, 2, 5, 2, 7, 0, "Lijab");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&L, &T1, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&L);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "Gciab");
  dpd_buf4_init(&R, CC_GR, R_irr, 2, 7, 2, 7, 0, "Rijab");
  dpd_contract444(&Z, &R, &G, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&R);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 3, rho_CiAb -= 0.5 LMnCe RMnAb tie */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 0, 11, 0, 0, "Z(Ci,Mn)");
  dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&L, &T1, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&L);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_init(&R, CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
  dpd_contract444(&Z, &R, &G, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&R);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 3, rho_cIaB -= 0.5 LmNcE RmNaB tIE */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 0, 11, 0, 0, "Z(cI,mN)");
  dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&L, &T1, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&L);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_init(&R, CC_GR, R_irr, 0, 5, 0, 5, 0, "RiJaB");
  dpd_contract444(&Z, &R, &G, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&R);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);

  /* term 4, rho_CIAB -= 0.5 LMNCE TauMNAB RIE */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 2, 11, 2, 0, "Z(CI,MN)");
  dpd_buf4_init(&L, CC_GL, L_irr, 2, 5, 2, 7, 0, "LIJAB");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&L, &R1, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&L);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  dpd_contract444(&Z, &T, &G, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&T);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 4, rho_ciab -= 0.5 Lmnce Taumnab Rie */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 2, 11, 2, 0, "Z(ci,mn)");
  dpd_buf4_init(&L, CC_GL, L_irr, 2, 5, 2, 7, 0, "Lijab");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&L, &R1, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&L);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "Gciab");
  dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
  dpd_contract444(&Z, &T, &G, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&T);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 4, rho_CiAb -= 0.5 LMnCe TauMnAb Rie */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 0, 11, 0, 0, "Z(Ci,Mn)");
  dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&L, &R1, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&L);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  dpd_contract444(&Z, &T, &G, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&T);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 4, rho_cIaB -= 0.5 LmNcE TaumNaB RIE */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 0, 11, 0, 0, "Z(cI,mN)");
  dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&L, &R1, &Z, 3, 1, 1, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&L);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauiJaB");
  dpd_contract444(&Z, &T, &G, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&T);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  /* term 5, rho_ciab += - (Lmnec Rme) (Tniab + P(ij) Tna Tib) */
  if (!params.connect_xi) {
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  dpd_buf4_init(&Tau, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tauIJAB");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
  dpd_contract244(&I1, &Tau, &G, 0, 0, 0, 1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&Tau);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "Gciab");
  dpd_buf4_init(&Tau, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tauijab");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
  dpd_contract244(&I1, &Tau, &G, 0, 0, 0, 1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&Tau);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_init(&Tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_OV");
  dpd_contract244(&I1, &Tau, &G, 0, 0, 0, 1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&Tau);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_init(&Tau, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauiJaB");
  dpd_file2_init(&I1, EOM_TMP, G_irr, 0, 1, "L2R1_ov");
  dpd_contract244(&I1, &Tau, &G, 0, 0, 0, 1.0, 1.0);
  dpd_file2_close(&I1);
  dpd_buf4_close(&Tau);
  dpd_buf4_close(&G);
  }

  /* +P(ab) LR_VV(c,a) t(i,b) */
  x_Gciab_6();

  /* +P(ab) LR_TT(c,a) R(i,b) */
  x_Gciab_7();

  /* -P(ab) Lmnce Rinae Tmb, term 8 */
  /* -P(ab) Lmnce Tinae Rmb, term 9 */
  x_Gciab_8_rohf();

  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);

  /* term 10, -P(AB) LNMCE TIE TNA RMB */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 2, 11, 2, 11, 0, "Z(NM,CI)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_buf4_init(&L, CC_GL, L_irr, 2, 5, 2, 7, 0, "LIJAB");
  dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_buf4_close(&L);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 0, 11, 2, 11, 0, "Z(NM,CI)");
  dpd_buf4_init(&Z2, EOM_TMP1, L_irr, 11, 11, 11, 11, 0, "Z(CI,AM)");
  dpd_contract244(&T1, &Z, &Z2, 0, 0, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(CI,AB)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&Z2, &R1, &Z, 3, 0, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "GCIAB");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 11, 5, "Z(CI,BA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(CI,BA)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  /* term 10, -P(ab) lmnce tie tna rmb */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 2, 11, 2, 11, 0, "Z(nm,ci)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_buf4_init(&L, CC_GL, L_irr, 2, 5, 2, 7, 0, "Lijab");
  dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_buf4_close(&L);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 0, 11, 2, 11, 0, "Z(nm,ci)");
  dpd_buf4_init(&Z2, EOM_TMP1, L_irr, 11, 11, 11, 11, 0, "Z(ci,am)");
  dpd_contract244(&T1, &Z, &Z2, 0, 0, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(ci,ab)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&Z2, &R1, &Z, 3, 0, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&Z2);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "Gciab");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 11, 5, "Z(ci,ba)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(ci,ba)");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  /* term 10, GCiAb -= P(AB) LNmCe Tie TNA Rmb */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(Nm,Ci)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_buf4_close(&L);
  dpd_file2_close(&T1);
  dpd_buf4_init(&Z2, EOM_TMP1, L_irr, 11, 11, 11, 11, 0, "Z(Ci,Am)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z, &Z2, 0, 0, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&Z2, &R1, &G, 3, 0, 0, -1.0, 1.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&G);

  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(nM,Ci)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJAb");
  dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_buf4_close(&L);
  dpd_file2_close(&T1);
  dpd_buf4_init(&Z2, EOM_TMP1, L_irr, 11, 11, 11, 11, 0, "Z(Ci,bM)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract244(&T1, &Z, &Z2, 0, 0, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(Ci,bA)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&Z2, &R1, &Z, 3, 0, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&Z2);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 11, 5, "Z(Ci,Ab)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(Ci,Ab)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);

  /* term 10, GcIaB - LnMcE TIE Tna RMB + LNmcE TIE TNB Rma */
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(nM,cI)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LiJaB");
  dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_buf4_close(&L);
  dpd_file2_close(&T1);
  dpd_buf4_init(&Z2, EOM_TMP1, L_irr, 11, 11, 11, 11, 0, "Z(cI,aM)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract244(&T1, &Z, &Z2, 0, 0, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&Z2, &R1, &G, 3, 0, 0, -1.0, 1.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&Z2);
  dpd_buf4_close(&G);

  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 0, 11, 0, 11, 0, "Z(Nm,cI)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_buf4_init(&L, CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjaB");
  dpd_contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  dpd_buf4_close(&L);
  dpd_file2_close(&T1);
  dpd_buf4_init(&Z2, EOM_TMP1, L_irr, 11, 11, 11, 11, 0, "Z(cI,Bm)");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1, &Z, &Z2, 0, 0, 1, 1.0, 0.0);
  dpd_file2_close(&T1);
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(cI,Ba)");
  dpd_file2_init(&R1, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&Z2, &R1, &Z, 3, 0, 0, 1.0, 0.0);
  dpd_file2_close(&R1);
  dpd_buf4_close(&Z2);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 11, 5, "Z(cI,aB)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(cI,aB)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);


  /* add 1/2 to ground-state parts of density */
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 7, 11, 7, 0, "Gciab");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 11, 7, 11, 7, 0, "Gciab");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_buf4_axpy(&G, &V, 0.5);
  dpd_buf4_close(&V);
  dpd_buf4_close(&G);

  /* clear out temporary files */
  psio_close(EOM_TMP0, 0);
  psio_open(EOM_TMP0, PSIO_OPEN_NEW);

  /*
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 11, 7, 11, 7, 0, "GCIAB");
  value = dpd_buf4_dot_self(&V);
  dpd_buf4_close(&V);
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 11, 7, 11, 7, 0, "Gciab");
  value += dpd_buf4_dot_self(&V);
  dpd_buf4_close(&V);
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  value += dpd_buf4_dot_self(&V);
  dpd_buf4_close(&V);
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  value += dpd_buf4_dot_self(&V);
  dpd_buf4_close(&V);
  fprintf(outfile,"<Gciab|Gciab> = %15.10lf\n",value);
  */

}



/* This function computes term 6,
   rho_ciab += P(ab) LR1_VV(c,a) T(i,b)
*/

void x_Gciab_6(void) { 
  int h, nirreps, c, i, a, b, C, I, A, B, Csym, Isym, Asym, Bsym, row, col;
  int AA, BB, AAsym, BBsym;
  int L_irr, R_irr, G_irr;
  dpdfile2 LR1A, LR1B, T1A, T1B;
  dpdbuf4 G;
 
  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* open one-electron files for the nasty terms */
  dpd_file2_init(&LR1A, EOM_TMP, G_irr, 1, 1, "LR_VV");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  if (params.ref == 0 || params.ref == 1) {
    dpd_file2_init(&LR1B, EOM_TMP, G_irr, 1, 1, "LR_vv");
    dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  }
  else {
    dpd_file2_init(&LR1B, EOM_TMP, G_irr, 3, 3, "LR_vv");
    dpd_file2_init(&T1B, CC_OEI, 0, 2, 3, "tia");
  }
  dpd_file2_mat_init(&T1A);   dpd_file2_mat_init(&T1B);
  dpd_file2_mat_init(&LR1A);   dpd_file2_mat_init(&LR1B);
  dpd_file2_mat_rd(&T1A);     dpd_file2_mat_rd(&T1B);
  dpd_file2_mat_rd(&LR1A);     dpd_file2_mat_rd(&LR1B);

  /* rho_CIAB += LR1_VV(C,A) T(I,B) - LR1_VV(C,B) T(I,A) */
  if (params.ref == 0 || params.ref == 1)
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "GCIAB");
  else
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 21, 5, 21, 7, 0, "GCIAB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      c = G.params->roworb[h][row][0];
      i = G.params->roworb[h][row][1];
      C = LR1A.params->rowidx[c]; Csym = LR1A.params->psym[c];
      I = T1A.params->rowidx[i]; Isym = T1A.params->psym[i];
      for(col=0; col < G.params->coltot[h]; col++) {
        a = G.params->colorb[h^G_irr][col][0];
        b = G.params->colorb[h^G_irr][col][1];
        A = LR1A.params->colidx[a]; Asym = LR1A.params->qsym[a];
        B = T1A.params->colidx[b]; Bsym = T1A.params->qsym[b];
        AA = T1A.params->colidx[a]; AAsym = T1A.params->qsym[a];
        BB = LR1A.params->colidx[b]; BBsym = LR1A.params->qsym[b];
        if( ((Csym^Asym)==G_irr) && (Isym==Bsym))
          G.matrix[h][row][col] +=
            LR1A.matrix[Csym][C][A] * T1A.matrix[Isym][I][B];
        if( ((Csym^BBsym)==G_irr) && (Isym==AAsym))
          G.matrix[h][row][col] -=
            LR1A.matrix[Csym][C][BB] * T1A.matrix[Isym][I][AA];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);
  /* rho_ciab += LR1_vv(c,a) T(i,b) - LR1_vv(c,b) T(i,a) */
  if (params.ref == 0 || params.ref == 1)
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "Gciab");
  else
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 31, 15, 31, 17, 0, "Gciab");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      c = G.params->roworb[h][row][0];
      i = G.params->roworb[h][row][1];
      C = LR1B.params->rowidx[c]; Csym = LR1B.params->psym[c];
      I = T1B.params->rowidx[i]; Isym = T1B.params->psym[i];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        a = G.params->colorb[h^G_irr][col][0];
        b = G.params->colorb[h^G_irr][col][1];
        A = LR1B.params->colidx[a]; Asym = LR1B.params->qsym[a];
        B = T1B.params->colidx[b]; Bsym = T1B.params->qsym[b];
        AA = T1B.params->colidx[a]; AAsym = T1B.params->qsym[a];
        BB = LR1B.params->colidx[b]; BBsym = LR1B.params->qsym[b];
        if( ((Csym^Asym)==G_irr) && (Isym==Bsym))
          G.matrix[h][row][col] +=
            LR1B.matrix[Csym][C][A] * T1B.matrix[Isym][I][B];
        if( ((Csym^BBsym)==G_irr) && (Isym==AAsym))
          G.matrix[h][row][col] -=
            LR1B.matrix[Csym][C][BB] * T1B.matrix[Isym][I][AA];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);
  /* rho_CiAb += LR1_VV(C,A) T(i,b) */
  if (params.ref == 0 || params.ref == 1)
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  else
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 26, 28, 26, 28, 0, "GCiAb");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      c = G.params->roworb[h][row][0];
      i = G.params->roworb[h][row][1];
      C = LR1A.params->rowidx[c]; Csym = LR1A.params->psym[c];
      I = T1B.params->rowidx[i]; Isym = T1B.params->psym[i];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        a = G.params->colorb[h^G_irr][col][0];
        b = G.params->colorb[h^G_irr][col][1];
        A = LR1A.params->colidx[a]; Asym = LR1A.params->qsym[a];
        B = T1B.params->colidx[b]; Bsym = T1B.params->qsym[b];
        if( ((Csym^Asym)==G_irr) && (Isym==Bsym))
          G.matrix[h][row][col] +=
            LR1A.matrix[Csym][C][A] * T1B.matrix[Isym][I][B];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);
  /* rho_cIaB += LR1_vv(c,a) T(I,B) */
  if (params.ref == 0 || params.ref == 1)
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  else
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 25, 29, 25, 29, 0, "GcIaB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      c = G.params->roworb[h][row][0];
      i = G.params->roworb[h][row][1];
      C = LR1B.params->rowidx[c]; Csym = LR1B.params->psym[c];
      I = T1A.params->rowidx[i]; Isym = T1A.params->psym[i];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        a = G.params->colorb[h^G_irr][col][0];
        b = G.params->colorb[h^G_irr][col][1];
        A = LR1B.params->colidx[a]; Asym = LR1B.params->qsym[a];
        B = T1A.params->colidx[b]; Bsym = T1A.params->qsym[b];
        if( ((Csym^Asym)==G_irr) && (Isym==Bsym))
          G.matrix[h][row][col] +=
            LR1B.matrix[Csym][C][A] * T1A.matrix[Isym][I][B];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  dpd_file2_mat_close(&LR1A);
  dpd_file2_mat_close(&LR1B);
  dpd_file2_close(&LR1A);
  dpd_file2_close(&LR1B);

  dpd_file2_mat_close(&T1A);
  dpd_file2_mat_close(&T1B);
  dpd_file2_close(&T1A);
  dpd_file2_close(&T1B);

  return;
}



/* This function computes term 7,
   rho_ciab += P(ab) LT_VV(c,a) R(i,b)
*/

void x_Gciab_7(void) { 
  int h, nirreps, c, i, a, b, C, I, A, B, Csym, Isym, Asym, Bsym, row, col;
  int AA, BB, AAsym, BBsym;
  int L_irr, R_irr, G_irr;
  dpdfile2 LT1A, LT1B, R1A, R1B;
  dpdbuf4 G;
 
  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* open one-electron files for the nasty terms */
  dpd_file2_init(&LT1A, EOM_TMP, L_irr, 1, 1, "LT_VV");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  if (params.ref == 0 || params.ref == 1) {
    dpd_file2_init(&LT1B, EOM_TMP, L_irr, 1, 1, "LT_vv");
    dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  }
  else {
    dpd_file2_init(&LT1B, EOM_TMP, L_irr, 3, 3, "LT_vv");
    dpd_file2_init(&R1B, CC_GR, R_irr, 2, 3, "Ria");
  }
  dpd_file2_mat_init(&R1A);   dpd_file2_mat_init(&R1B);
  dpd_file2_mat_init(&LT1A);   dpd_file2_mat_init(&LT1B);
  dpd_file2_mat_rd(&R1A);     dpd_file2_mat_rd(&R1B);
  dpd_file2_mat_rd(&LT1A);     dpd_file2_mat_rd(&LT1B);

  /* rho_CIAB += LT_VV(C,A) R(I,B) - LT_VV(C,B) R(I,A) */
  if (params.ref == 0 || params.ref == 1)
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "GCIAB");
  else
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 21, 5, 21, 7, 0, "GCIAB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      c = G.params->roworb[h][row][0];
      i = G.params->roworb[h][row][1];
      C = LT1A.params->rowidx[c]; Csym = LT1A.params->psym[c];
      I = R1A.params->rowidx[i]; Isym = R1A.params->psym[i];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        a = G.params->colorb[h^G_irr][col][0];
        b = G.params->colorb[h^G_irr][col][1];
        A = LT1A.params->colidx[a]; Asym = LT1A.params->qsym[a];
        B = R1A.params->colidx[b]; Bsym = R1A.params->qsym[b];
        AA = R1A.params->colidx[a]; AAsym = R1A.params->qsym[a];
        BB = LT1A.params->colidx[b]; BBsym = LT1A.params->qsym[b];
        if( ((Csym^Asym)==L_irr) && ((Isym^Bsym)==R_irr) )
          G.matrix[h][row][col] +=
            LT1A.matrix[Csym][C][A] * R1A.matrix[Isym][I][B];
        if( ((Csym^BBsym)==L_irr) && ((Isym^AAsym)==R_irr) )
          G.matrix[h][row][col] -=
            LT1A.matrix[Csym][C][BB] * R1A.matrix[Isym][I][AA];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);
  /* rho_ciab += LT_vv(c,a) R(i,b) - LT_vv(c,b) R(i,a) */
  if (params.ref == 0 || params.ref == 1)
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "Gciab");
  else
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 31, 15, 31, 17, 0, "Gciab");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      c = G.params->roworb[h][row][0];
      i = G.params->roworb[h][row][1];
      C = LT1B.params->rowidx[c]; Csym = LT1B.params->psym[c];
      I = R1B.params->rowidx[i]; Isym = R1B.params->psym[i];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        a = G.params->colorb[h^G_irr][col][0];
        b = G.params->colorb[h^G_irr][col][1];
        A = LT1B.params->colidx[a]; Asym = LT1B.params->qsym[a];
        B = R1B.params->colidx[b]; Bsym = R1B.params->qsym[b];
        AA = R1B.params->colidx[a]; AAsym = R1B.params->qsym[a];
        BB = LT1B.params->colidx[b]; BBsym = LT1B.params->qsym[b];
        if( ((Csym^Asym)==L_irr) && ((Isym^Bsym)==R_irr) )
          G.matrix[h][row][col] +=
            LT1B.matrix[Csym][C][A] * R1B.matrix[Isym][I][B];
        if( ((Csym^BBsym)==L_irr) && ((Isym^AAsym)==R_irr) )
          G.matrix[h][row][col] -=
            LT1B.matrix[Csym][C][BB] * R1B.matrix[Isym][I][AA];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);
  /* rho_CiAb += LT_VV(C,A) R(i,b) */
  if (params.ref == 0 || params.ref == 1)
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  else
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 26, 28, 26, 28, 0, "GCiAb");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      c = G.params->roworb[h][row][0];
      i = G.params->roworb[h][row][1];
      C = LT1A.params->rowidx[c]; Csym = LT1A.params->psym[c];
      I = R1B.params->rowidx[i]; Isym = R1B.params->psym[i];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        a = G.params->colorb[h^G_irr][col][0];
        b = G.params->colorb[h^G_irr][col][1];
        A = LT1A.params->colidx[a]; Asym = LT1A.params->qsym[a];
        B = R1B.params->colidx[b]; Bsym = R1B.params->qsym[b];
        if( ((Csym^Asym)==L_irr) && ((Isym^Bsym)==R_irr) )
          G.matrix[h][row][col] +=
            LT1A.matrix[Csym][C][A] * R1B.matrix[Isym][I][B];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);
  /* rho_cIaB += LT_vv(c,a) R(I,B) */
  if (params.ref == 0 || params.ref == 1)
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  else
    dpd_buf4_init(&G, EOM_TMP0, G_irr, 25, 29, 25, 29, 0, "GcIaB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&G, h);
    dpd_buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      c = G.params->roworb[h][row][0];
      i = G.params->roworb[h][row][1];
      C = LT1B.params->rowidx[c]; Csym = LT1B.params->psym[c];
      I = R1A.params->rowidx[i]; Isym = R1A.params->psym[i];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        a = G.params->colorb[h^G_irr][col][0];
        b = G.params->colorb[h^G_irr][col][1];
        A = LT1B.params->colidx[a]; Asym = LT1B.params->qsym[a];
        B = R1A.params->colidx[b]; Bsym = R1A.params->qsym[b];
        if( ((Csym^Asym)==L_irr) && ((Isym^Bsym)==R_irr) )
          G.matrix[h][row][col] +=
            LT1B.matrix[Csym][C][A] * R1A.matrix[Isym][I][B];
      }
    }
    dpd_buf4_mat_irrep_wrt(&G, h);
    dpd_buf4_mat_irrep_close(&G, h);
  }
  dpd_buf4_close(&G);

  dpd_file2_mat_close(&LT1A);
  dpd_file2_mat_close(&LT1B);
  dpd_file2_close(&LT1A);
  dpd_file2_close(&LT1B);

  dpd_file2_mat_close(&R1A);
  dpd_file2_mat_close(&R1B);
  dpd_file2_close(&R1A);
  dpd_file2_close(&R1B);

  return;
}



/* This function computes terms of excited Gciab
   term 8  +P(ab) Lmnce Rinae Tmb
   term 9, +P(AB) LMNCE TINAE RMB
*/

void x_Gciab_8_rohf(void) { 
  int h, nirreps, i, j, k, a, I, J, K, A, Isym, Jsym, Ksym, Asym, row, col;
  int II,JJ,IIsym,JJsym;
  int L_irr, R_irr, G_irr;
  double value;
  dpdfile2 L1A, T1A, L1B, T1B, R1A, R1B, I1A, I1B;
  dpdbuf4 G, V, T, L, Z, Z1, Z2, Tau;
 
  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* term 8, +P(AB) LMNCE RINAE TMB */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 5, 10, 5, 0, "Z(IA,BC)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1A, &V, &Z, 0, 2, 1, 1.0, 0.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, sprq, 11, 5, "Z(CI,BA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(CI,BA)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "GCIAB");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 11, 5, "Z(CI,AB)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(CI,AB)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 8, +P(ab) Lmnce Rinae Tmb */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 5, 10, 5, 0, "Z(ia,bc)");
  dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovov");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tia");
  dpd_contract244(&T1A, &V, &Z, 0, 2, 1, 1.0, 0.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, sprq, 11, 5, "Z(ci,ba)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(ci,ba)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "Gciab");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 11, 5, "Z(ci,ab)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(ci,ab)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  
  /* term 8, GCiAb -= LmNCe RiNAe tmb */
  dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_oVoV");
  dpd_buf4_sort(&V, EOM_TMP1, spqr, 11, 11, "Z(Ci,Am)");
  dpd_buf4_close(&V);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 11, 11, 11, 0, "Z(Ci,Am)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_contract424(&Z, &T1B, &G, 3, 0, 0, -1.0, 1.0); 
  dpd_file2_close(&T1B);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 8, GCiAb -= (LMnCe Rinbe + LMNCE RiNbE) TMA */
  dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovOV");
  dpd_buf4_sort(&V, EOM_TMP1, sprq, 11, 10, "Z(Ci,Mb)");
  dpd_buf4_close(&V);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z(Ci,Mb)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract244(&T1A, &Z, &G, 0, 2, 1, 1.0, 1.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  /* term 8, GcIaB -= LMncE RInaE tMB */
  dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  dpd_buf4_sort(&V, EOM_TMP1, spqr, 11, 11, "Z(cI,aM)");
  dpd_buf4_close(&V);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 11, 11, 11, 0, "Z(cI,aM)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract424(&Z, &T1A, &G, 3, 0, 0, -1.0, 1.0); 
  dpd_file2_close(&T1A);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 8, GcIaB -= (LmNcE RINBE + Lmnce RInBe) Tma */
  dpd_buf4_init(&V, EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
  dpd_buf4_sort(&V, EOM_TMP1, sprq, 11, 10, "Z(cI,mB)");
  dpd_buf4_close(&V);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 10, 11, 10, 0, "Z(cI,mB)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_file2_init(&T1B, CC_OEI, 0, 0, 1, "tia");
  dpd_contract244(&T1B, &Z, &G, 0, 2, 1, 1.0, 1.0); 
  dpd_file2_close(&T1B);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  psio_close(EOM_TMP1, 0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);

  /* term 9, +P(AB) LMNCE TINAE RMB */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 5, 10, 5, 0, "Z(IA,BC)");
  dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIAJB");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract244(&R1A, &V, &Z, 0, 2, 1, 1.0, 0.0); 
  dpd_file2_close(&R1A);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, sprq, 11, 5, "Z(CI,BA)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(CI,BA)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "GCIAB");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 11, 5, "Z(CI,AB)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(CI,AB)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 9, +P(ab) Lmnce Tinae Rmb */
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 10, 5, 10, 5, 0, "Z(ia,bc)");
  dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "Viajb");
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract244(&R1B, &V, &Z, 0, 2, 1, 1.0, 0.0); 
  dpd_file2_close(&R1B);
  dpd_buf4_close(&V);
  dpd_buf4_sort(&Z, EOM_TMP1, sprq, 11, 5, "Z(ci,ba)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(ci,ba)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 7, 0, "Gciab");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_sort(&Z, EOM_TMP1, pqsr, 11, 5, "Z(ci,ab)");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 11, 5, 11, 5, 0, "Z(ci,ab)");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  
  /* term 9, GCiAb -= LmNCe TiNAe Rmb */
  dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "ViAjB");
  dpd_buf4_sort(&V, EOM_TMP1, spqr, 11, 11, "Z(Ci,Am)");
  dpd_buf4_close(&V);
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 11, 11, 11, 0, "Z(Ci,Am)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract424(&Z, &R1B, &G, 3, 0, 0, -1.0, 1.0); 
  dpd_file2_close(&R1B);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 9, GCiAb -= (LMnCe Tinbe + LMNCE TiNbE) RMA */
  dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "ViaJB");
  dpd_buf4_sort(&V, EOM_TMP1, sprq, 11, 10, "Z(Ci,Mb)");
  dpd_buf4_close(&V);
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 10, 11, 10, 0, "Z(Ci,Mb)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GCiAb");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract244(&R1A, &Z, &G, 0, 2, 1, 1.0, 1.0); 
  dpd_file2_close(&R1A);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  /* term 9, GcIaB -= LMncE TInaE RMB */
  dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIaJb");
  dpd_buf4_sort(&V, EOM_TMP1, spqr, 11, 11, "Z(cI,aM)");
  dpd_buf4_close(&V);
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 11, 11, 11, 0, "Z(cI,aM)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_file2_init(&R1A, CC_GR, R_irr, 0, 1, "RIA");
  dpd_contract424(&Z, &R1A, &G, 3, 0, 0, -1.0, 1.0); 
  dpd_file2_close(&R1A);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);
  /* term 9, GcIaB -= (LmNcE TINBE + Lmnce TInBe) Rma */
  dpd_buf4_init(&V, EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIAjb");
  dpd_buf4_sort(&V, EOM_TMP1, sprq, 11, 10, "Z(cI,mB)");
  dpd_buf4_close(&V);
  dpd_buf4_init(&Z, EOM_TMP1, L_irr, 11, 10, 11, 10, 0, "Z(cI,mB)");
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 11, 5, 11, 5, 0, "GcIaB");
  dpd_file2_init(&R1B, CC_GR, R_irr, 0, 1, "Ria");
  dpd_contract244(&R1B, &Z, &G, 0, 2, 1, 1.0, 1.0); 
  dpd_file2_close(&R1B);
  dpd_buf4_close(&G);
  dpd_buf4_close(&Z);

  return;
}

}} // namespace psi::ccdensity
