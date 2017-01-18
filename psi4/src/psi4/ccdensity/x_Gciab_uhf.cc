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

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void x_Gciab_uhf(void);
extern void x_Gciab_6(void);
extern void x_Gciab_7(void);
void x_Gciab_8_uhf(void);

/* This function computes the non-R0 parts of the 2pdm density matrix
   Gciab = 0.5 *(rho_abci + rho_ciab) */

void x_Gciab_uhf(void) {
  int h, nirreps, i, j, k, a, I, J, K, A, Isym, Jsym, Ksym, Asym, row, col;
  int II,JJ,IIsym,JJsym;
  int L_irr, R_irr, G_irr;
  double value;
  dpdfile2 L1, T1, R1, I1;
  dpdbuf4 G, V, T, L, Z, Z2, R, Tau;

  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* term 1, rho_abci += Lmiab * Rmc */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 7, 21, 7, 21, 0, "L2R1_VVOV(pqsr)");
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP0, rspq, 21, 7, "GCIAB");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 17, 31, 17, 31, 0, "L2R1_vvov(pqsr)");
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP0, rspq, 31, 17, "Gciab");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 28, 26, 28, 26, 0, "L2R1_VvoV(pqsr)");
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP0, rspq, 26, 28, "GCiAb");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 28, 25, 28, 25, 0, "L2R1_VvOv(pqsr)");
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP0, rsqp, 25, 29, "GcIaB");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP0, G_irr, 21, 7, 21, 7, 0, "GCIAB");
  global_dpd_->buf4_scm(&Z, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP0, G_irr, 31, 17, 31, 17, 0, "Gciab");
  global_dpd_->buf4_scm(&Z, -1.0);
  global_dpd_->buf4_close(&Z);

  /* change to sort_axpy */
  /* term 2, rho_ciab += Rmiab * Lmc */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 7, 21, 7, 21, 0, "R2L1_VVOV(pqsr)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_TMP0, rspq, 21, 7, "GCIAB", -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 17, 31, 17, 31, 0, "R2L1_vvov(pqsr)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_TMP0, rspq, 31, 17, "Gciab", -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 28, 26, 28, 26, 0, "R2L1_VvoV(pqsr)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_TMP0, rspq, 26, 28, "GCiAb", 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 28, 25, 28, 25, 0, "R2L1_VvOv(pqsr)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_TMP0, rsqp, 25, 29, "GcIaB", 1.0);
  global_dpd_->buf4_close(&Z);
  /* two step code - perhaps suggestive of future out of core approach
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 7, 21, 7, 21, 0, "R2L1_VVOV(pqsr)");
  dpd_buf4_sort(&Z, EOM_TMP1, rspq, 21, 7, "GCIAB");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 17, 31, 17, 31, 0, "R2L1_vvov(pqsr)");
  dpd_buf4_sort(&Z, EOM_TMP1, rspq, 31, 17, "Gciab");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 28, 26, 28, 26, 0, "R2L1_VvoV(pqsr)");
  dpd_buf4_sort(&Z, EOM_TMP1, rspq, 26, 28, "GCiAb");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&Z, EOM_TMP, G_irr, 28, 25, 28, 25, 0, "R2L1_VvOv(pqsr)");
  dpd_buf4_sort(&Z, EOM_TMP1, rsqp, 25, 29, "GcIaB");
  dpd_buf4_close(&Z);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 21, 7, 21, 7, 0, "GCIAB");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 21, 7, 21, 7, 0, "GCIAB");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 31, 17, 31, 17, 0, "Gciab");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 31, 17, 31, 17, 0, "Gciab");
  dpd_buf4_axpy(&Z, &G, -1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 26, 28, 26, 28, 0, "GCiAb");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 26, 28, 26, 28, 0, "GCiAb");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  dpd_buf4_init(&G, EOM_TMP0, G_irr, 25, 29, 25, 29, 0, "GcIaB");
  dpd_buf4_init(&Z, EOM_TMP1, G_irr, 25, 29, 25, 29, 0, "GcIaB");
  dpd_buf4_axpy(&Z, &G, 1.0);
  dpd_buf4_close(&Z);
  dpd_buf4_close(&G);
  psio_close(EOM_TMP1,0);
  psio_open(EOM_TMP1, PSIO_OPEN_NEW);
  */

  /* term 3, rho_CIAB -= 0.5 LMNCE RMNAB tIE */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 21, 2, 21, 2, 0, "Z(CI,MN)");
  global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 2, 5, 2, 7, 0, "LIJAB");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&L, &T1, &Z, 3, 1, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&L);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 21, 7, 21, 7, 0, "GCIAB");
  global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 2, 7, 2, 7, 0, "RIJAB");
  global_dpd_->contract444(&Z, &R, &G, 0, 1, -1.0, 1.0);
  global_dpd_->buf4_close(&R);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);
  /* term 3, rho_ciab -= 0.5 Lmnce Rmnab tie */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 31, 12, 31, 12, 0, "Z(ci,mn)");
  global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 12, 15, 12, 17, 0, "Lijab");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&L, &T1, &Z, 3, 1, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&L);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 31, 17, 31, 17, 0, "Gciab");
  global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 12, 17, 12, 17, 0, "Rijab");
  global_dpd_->contract444(&Z, &R, &G, 0, 1, -1.0, 1.0);
  global_dpd_->buf4_close(&R);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);
  /* term 3, rho_CiAb -= 0.5 LMnCe RMnAb tie */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 26, 22, 26, 22, 0, "Z(Ci,Mn)");
  global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&L, &T1, &Z, 3, 1, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&L);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 26, 28, 26, 28, 0, "GCiAb");
  global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 22, 28, 22, 28, 0, "RIjAb");
  global_dpd_->contract444(&Z, &R, &G, 0, 1, -1.0, 1.0);
  global_dpd_->buf4_close(&R);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);
  /* term 3, rho_cIaB -= 0.5 LmNcE RmNaB tIE */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 25, 23, 25, 23, 0, "Z(cI,mN)");
  global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 23, 29, 23, 29, 0, "LiJaB");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&L, &T1, &Z, 3, 1, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&L);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 25, 29, 25, 29, 0, "GcIaB");
  global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 23, 29, 23, 29, 0, "RiJaB");
  global_dpd_->contract444(&Z, &R, &G, 0, 1, -1.0, 1.0);
  global_dpd_->buf4_close(&R);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);

  /* term 4, rho_CIAB -= 0.5 LMNCE TauMNAB RIE */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 21, 2, 21, 2, 0, "Z(CI,MN)");
  global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 2, 5, 2, 7, 0, "LIJAB");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&L, &R1, &Z, 3, 1, 1, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&L);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 21, 7, 21, 7, 0, "GCIAB");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  global_dpd_->contract444(&Z, &T, &G, 0, 1, -1.0, 1.0);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);
  /* term 4, rho_ciab -= 0.5 Lmnce Taumnab Rie */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 31, 12, 31, 12, 0, "Z(ci,mn)");
  global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 12, 15, 12, 17, 0, "Lijab");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract424(&L, &R1, &Z, 3, 1, 1, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&L);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 31, 17, 31, 17, 0, "Gciab");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
  global_dpd_->contract444(&Z, &T, &G, 0, 1, -1.0, 1.0);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);
  /* term 4, rho_CiAb -= 0.5 LMnCe TauMnAb Rie */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 26, 22, 26, 22, 0, "Z(Ci,Mn)");
  global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract424(&L, &R1, &Z, 3, 1, 1, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&L);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 26, 28, 26, 28, 0, "GCiAb");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
  global_dpd_->contract444(&Z, &T, &G, 0, 1, -1.0, 1.0);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);
  /* term 4, rho_cIaB -= 0.5 LmNcE TaumNaB RIE */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 25, 23, 25, 23, 0, "Z(cI,mN)");
  global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 23, 29, 23, 29, 0, "LiJaB");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&L, &R1, &Z, 3, 1, 1, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&L);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 25, 29, 25, 29, 0, "GcIaB");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tauiJaB");
  global_dpd_->contract444(&Z, &T, &G, 0, 1, -1.0, 1.0);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);

  /* term 5, rho_ciab += - (Lmnec Rme) (Tniab + P(ij) Tna Tib) */
  if (!params.connect_xi) {
    global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 21, 7, 21, 7, 0, "GCIAB");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->contract244(&I1, &Tau, &G, 0, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&I1);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 31, 17, 31, 17, 0, "Gciab");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tauijab");
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    global_dpd_->contract244(&I1, &Tau, &G, 0, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&I1);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 26, 28, 26, 28, 0, "GCiAb");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->contract244(&I1, &Tau, &G, 0, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&I1);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 25, 29, 25, 29, 0, "GcIaB");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tauiJaB");
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    global_dpd_->contract244(&I1, &Tau, &G, 0, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&I1);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&G);
  }

  /* +P(ab) LR_VV(c,a) t(i,b) */
  x_Gciab_6();

  /* +P(ab) LR_TT(c,a) R(i,b) */
  x_Gciab_7();

  /* -P(ab) Lmnce Rinae Tmb, term 8 */
  /* -P(ab) Lmnce Tinae Rmb, term 9 */
  x_Gciab_8_uhf();

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);

  /* term 10, -P(AB) LNMCE TIE TNA RMB */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 2, 21, 2, 21, 0, "Z(NM,CI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 2, 5, 2, 7, 0, "LIJAB");
  global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&L);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 0, 21, 2, 21, 0, "Z(NM,CI)");
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, L_irr, 21, 21, 21, 21, 0, "Z(CI,AM)");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 0, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 21, 5, 21, 5, 0, "Z(CI,AB)");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&Z2, &R1, &Z, 3, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 21, 5, 21, 7, 0, "GCIAB");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, pqsr, 21, 5, "Z(CI,BA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 21, 5, 21, 5, 0, "Z(CI,BA)");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  /* term 10, -P(ab) lmnce tie tna rmb */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 12, 31, 12, 31, 0, "Z(nm,ci)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 12, 15, 12, 17, 0, "Lijab");
  global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&L);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 10, 31, 12, 31, 0, "Z(nm,ci)");
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, L_irr, 31, 31, 31, 31, 0, "Z(ci,am)");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 0, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 31, 15, 31, 15, 0, "Z(ci,ab)");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract424(&Z2, &R1, &Z, 3, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 31, 15, 31, 17, 0, "Gciab");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, pqsr, 31, 15, "Z(ci,ba)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 31, 15, 31, 15, 0, "Z(ci,ba)");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  /* term 10, GCiAb -= P(AB) LNmCe Tie TNA Rmb */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 22, 26, 22, 26, 0, "Z(Nm,Ci)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 22, 28, 22, 28, 0, "LIjAb");
  global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&L);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, L_irr, 26, 26, 26, 26, 0, "Z(Ci,Am)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 0, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 26, 28, 26, 28, 0, "GCiAb");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract424(&Z2, &R1, &G, 3, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 23, 26, 23, 26, 0, "Z(nM,Ci)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 23, 28, 23, 28, 0, "LiJAb");
  global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&L);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, L_irr, 26, 25, 26, 25, 0, "Z(Ci,bM)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 0, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 26, 29, 26, 29, 0, "Z(Ci,bA)");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&Z2, &R1, &Z, 3, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, pqsr, 26, 28, "Z(Ci,Ab)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 26, 28, 26, 28, 0, "Z(Ci,Ab)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 26, 28, 26, 28, 0, "GCiAb");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  /* term 10, GcIaB - LnMcE TIE Tna RMB + LNmcE TIE TNB Rma */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 23, 25, 23, 25, 0, "Z(nM,cI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 23, 29, 23, 29, 0, "LiJaB");
  global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&L);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, L_irr, 25, 25, 25, 25, 0, "Z(cI,aM)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 0, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 25, 29, 25, 29, 0, "GcIaB");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&Z2, &R1, &G, 3, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 22, 25, 22, 25, 0, "Z(Nm,cI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->buf4_init(&L, PSIF_CC_GL, L_irr, 22, 29, 22, 29, 0, "LIjaB");
  global_dpd_->contract424(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&L);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, L_irr, 25, 26, 25, 26, 0, "Z(cI,Bm)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 0, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 25, 28, 25, 28, 0, "Z(cI,Ba)");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract424(&Z2, &R1, &Z, 3, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, pqsr, 25, 29, "Z(cI,aB)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 25, 29, 25, 29, 0, "Z(cI,aB)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 25, 29, 25, 29, 0, "GcIaB");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);


  /* add 1/2 to ground-state parts of density */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 21, 7, 21, 7, 0, "GCIAB");
  global_dpd_->buf4_init(&V, PSIF_CC_GAMMA, G_irr, 21, 7, 21, 7, 0, "GCIAB");
  global_dpd_->buf4_axpy(&G, &V, 0.5);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 31, 17, 31, 17, 0, "Gciab");
  global_dpd_->buf4_init(&V, PSIF_CC_GAMMA, G_irr, 31, 17, 31, 17, 0, "Gciab");
  global_dpd_->buf4_axpy(&G, &V, 0.5);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 26, 28, 26, 28, 0, "GCiAb");
  global_dpd_->buf4_init(&V, PSIF_CC_GAMMA, G_irr, 26, 28, 26, 28, 0, "GCiAb");
  global_dpd_->buf4_axpy(&G, &V, 0.5);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 25, 29, 25, 29, 0, "GcIaB");
  global_dpd_->buf4_init(&V, PSIF_CC_GAMMA, G_irr, 25, 29, 25, 29, 0, "GcIaB");
  global_dpd_->buf4_axpy(&G, &V, 0.5);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&G);

  /* clear out temporary files */
  psio_close(PSIF_EOM_TMP0, 0);
  psio_open(PSIF_EOM_TMP0, PSIO_OPEN_NEW);

  /*
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 21, 7, 21, 7, 0, "GCIAB");
  value = dpd_buf4_dot_self(&V);
  dpd_buf4_close(&V);
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 31, 17, 31, 17, 0, "Gciab");
  value += dpd_buf4_dot_self(&V);
  dpd_buf4_close(&V);
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 26, 28, 26, 28, 0, "GCiAb");
  value += dpd_buf4_dot_self(&V);
  dpd_buf4_close(&V);
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 25, 29, 25, 29, 0, "GcIaB");
  value += dpd_buf4_dot_self(&V);
  dpd_buf4_close(&V);
  outfile->Printf("<Gciab|Gciab> = %15.10lf\n",value);
  */

}



/* This function computes terms of excited Gciab
   term 8  +P(ab) Lmnce Rinae Tmb
   term 9, +P(AB) LMNCE TINAE RMB
*/

void x_Gciab_8_uhf(void) {
  int h, nirreps, i, j, k, a, I, J, K, A, Isym, Jsym, Ksym, Asym, row, col;
  int II,JJ,IIsym,JJsym;
  int L_irr, R_irr, G_irr;
  double value;
  dpdfile2 L1A, T1A, L1B, T1B, R1A, R1B, I1A, I1B;
  dpdbuf4 G, V, T, L, Z, Z1, Z2, Tau;

  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* term 8, +P(AB) LMNCE RINAE TMB */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 20, 5, 20, 5, 0, "Z(IA,BC)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 20, 20, 20, 20, 0, "R2L2_OVOV");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1A, &V, &Z, 0, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, sprq, 21, 5, "Z(CI,BA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 21, 5, 21, 5, 0, "Z(CI,BA)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 21, 5, 21, 7, 0, "GCIAB");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, pqsr, 21, 5, "Z(CI,AB)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 21, 5, 21, 5, 0, "Z(CI,AB)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);
  /* term 8, +P(ab) Lmnce Rinae Tmb */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 30, 15, 30, 15, 0, "Z(ia,bc)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 30, 30, 30, 30, 0, "R2L2_ovov");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract244(&T1A, &V, &Z, 0, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, sprq, 31, 15, "Z(ci,ba)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 31, 15, 31, 15, 0, "Z(ci,ba)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 31, 15, 31, 17, 0, "Gciab");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, pqsr, 31, 15, "Z(ci,ab)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 31, 15, 31, 15, 0, "Z(ci,ab)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);

  /* term 8, GCiAb -= LmNCe RiNAe tmb */
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 27, 27, 27, 27, 0, "R2L2_oVoV");
  global_dpd_->buf4_sort(&V, PSIF_EOM_TMP1, spqr, 26, 26, "Z(Ci,Am)");
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 26, 26, 26, 26, 0, "Z(Ci,Am)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 26, 28, 26, 28, 0, "GCiAb");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&Z, &T1B, &G, 3, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);
  /* term 8, GCiAb -= (LMnCe Rinbe + LMNCE RiNbE) TMA */
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 30, 20, 30, 20, 0, "R2L2_ovOV");
  global_dpd_->buf4_sort(&V, PSIF_EOM_TMP1, sprq, 26, 24, "Z(Ci,Mb)");
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 26, 24, 26, 24, 0, "Z(Ci,Mb)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 26, 28, 26, 28, 0, "GCiAb");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1A, &Z, &G, 0, 2, 1, 1.0, 1.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);

  /* term 8, GcIaB -= LMncE RInaE tMB */
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 24, 24, 24, 24, 0, "R2L2_OvOv");
  global_dpd_->buf4_sort(&V, PSIF_EOM_TMP1, spqr, 25, 25, "Z(cI,aM)");
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 25, 25, 25, 25, 0, "Z(cI,aM)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 25, 29, 25, 29, 0, "GcIaB");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&Z, &T1A, &G, 3, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);
  /* term 8, GcIaB -= (LmNcE RINBE + Lmnce RInBe) Tma */
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 20, 30, 20, 30, 0, "R2L2_OVov");
  global_dpd_->buf4_sort(&V, PSIF_EOM_TMP1, sprq, 25, 27, "Z(cI,mB)");
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 25, 27, 25, 27, 0, "Z(cI,mB)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 25, 29, 25, 29, 0, "GcIaB");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract244(&T1B, &Z, &G, 0, 2, 1, 1.0, 1.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);

  psio_close(PSIF_EOM_TMP1, 0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);

  /* term 9, +P(AB) LMNCE TINAE RMB */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 20, 5, 20, 5, 0, "Z(IA,BC)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 20, 20, 20, 20, 0, "VIAJB");
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract244(&R1A, &V, &Z, 0, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&R1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, sprq, 21, 5, "Z(CI,BA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 21, 5, 21, 5, 0, "Z(CI,BA)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 21, 5, 21, 7, 0, "GCIAB");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, pqsr, 21, 5, "Z(CI,AB)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 21, 5, 21, 5, 0, "Z(CI,AB)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);
  /* term 9, +P(ab) Lmnce Tinae Rmb */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 30, 15, 30, 15, 0, "Z(ia,bc)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 30, 30, 30, 30, 0, "Viajb");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract244(&R1B, &V, &Z, 0, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&R1B);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, sprq, 31, 15, "Z(ci,ba)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 31, 15, 31, 15, 0, "Z(ci,ba)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 31, 15, 31, 17, 0, "Gciab");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, pqsr, 31, 15, "Z(ci,ab)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 31, 15, 31, 15, 0, "Z(ci,ab)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);

  /* term 9, GCiAb -= LmNCe TiNAe Rmb */
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 27, 27, 27, 27, 0, "ViAjB");
  global_dpd_->buf4_sort(&V, PSIF_EOM_TMP1, spqr, 26, 26, "Z(Ci,Am)");
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 26, 26, 26, 26, 0, "Z(Ci,Am)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 26, 28, 26, 28, 0, "GCiAb");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract424(&Z, &R1B, &G, 3, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&R1B);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);
  /* term 9, GCiAb -= (LMnCe Tinbe + LMNCE TiNbE) RMA */
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 30, 20, 30, 20, 0, "ViaJB");
  global_dpd_->buf4_sort(&V, PSIF_EOM_TMP1, sprq, 26, 24, "Z(Ci,Mb)");
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 26, 24, 26, 24, 0, "Z(Ci,Mb)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 26, 28, 26, 28, 0, "GCiAb");
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract244(&R1A, &Z, &G, 0, 2, 1, 1.0, 1.0);
  global_dpd_->file2_close(&R1A);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);

  /* term 9, GcIaB -= LMncE TInaE RMB */
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 24, 24, 24, 24, 0, "VIaJb");
  global_dpd_->buf4_sort(&V, PSIF_EOM_TMP1, spqr, 25, 25, "Z(cI,aM)");
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 25, 25, 25, 25, 0, "Z(cI,aM)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 25, 29, 25, 29, 0, "GcIaB");
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&Z, &R1A, &G, 3, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&R1A);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);
  /* term 9, GcIaB -= (LmNcE TINBE + Lmnce TInBe) Rma */
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 20, 30, 20, 30, 0, "VIAjb");
  global_dpd_->buf4_sort(&V, PSIF_EOM_TMP1, sprq, 25, 27, "Z(cI,mB)");
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 25, 27, 25, 27, 0, "Z(cI,mB)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 25, 29, 25, 29, 0, "GcIaB");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract244(&R1B, &Z, &G, 0, 2, 1, 1.0, 1.0);
  global_dpd_->file2_close(&R1B);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);

  return;
}

}} // namespace psi::ccdensity
