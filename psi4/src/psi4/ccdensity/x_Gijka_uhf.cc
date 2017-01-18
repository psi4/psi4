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
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 2, 20, 2, 20, 0, "L2R1_OOVO(pqsr)");
  global_dpd_->buf4_copy(&Z, PSIF_EOM_TMP0, "GIJKA");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 12, 30, 12, 30, 0, "L2R1_oovo(pqsr)");
  global_dpd_->buf4_copy(&Z, PSIF_EOM_TMP0, "Gijka");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 22, 24, 22, 24, 0, "L2R1_OovO(pqsr)");
  global_dpd_->buf4_copy(&Z, PSIF_EOM_TMP0, "GIjKa");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 23, 27, 23, 27, 0, "L2R1_OoVo(qpsr)");
  global_dpd_->buf4_copy(&Z, PSIF_EOM_TMP0, "GiJkA");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  global_dpd_->buf4_scm(&Z, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  global_dpd_->buf4_scm(&Z, -1.0);
  global_dpd_->buf4_close(&Z);

  /* term 2, rho_ijka += Rijae * Lke */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 20, 2, 20, 0, "GIJKA");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 2, 20, 2, 20, 0, "L1R2_OOVO(pqsr)");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 12, 30, 12, 30, 0, "Gijka");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 12, 30, 12, 30, 0, "L1R2_oovo(pqsr)");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 22, 24, 22, 24, 0, "L1R2_OovO(pqsr)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 23, 27, 23, 27, 0, "L1R2_OoVo(qpsr)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  /* term 3, rho_ijka += 0.5 Rijef Lkmef tma  */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 20, 2, 20, 0, "GIJKA");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 2, 0, 2, 2, 0, "R2L2_OOOO");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&Z, &T1A, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 12, 30, 12, 30, 0, "Gijka");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 12, 10, 12, 12, 0, "R2L2_oooo");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&Z, &T1B, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 22, 22, 22, 22, 0, "R2L2_OoOo");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&Z, &T1B, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 23, 23, 23, 23, 0, "R2L2_oOoO");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&Z, &T1A, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  /* term 4, rho_ijka += 0.5 [Tijef + P(ij) Tie Tjf] Lkmef Rma */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 20, 2, 20, 0, "GIJKA");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, L_irr, 2, 0, 2, 2, 0, "Tau2L2_OOOO");
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&Z, &R1A, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&R1A);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 12, 30, 12, 30, 0, "Gijka");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, L_irr, 12, 10, 12, 12, 0, "Tau2L2_oooo");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract424(&Z, &R1B, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&R1B);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, L_irr, 22, 22, 22, 22, 0, "Tau2L2_OoOo");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract424(&Z, &R1B, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&R1B);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, L_irr, 23, 23, 23, 23, 0, "Tau2L2_oOoO");
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&Z, &R1A, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&R1A);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);


  /* term 5, rho_ijka += - (Lkmef Rmf) (Tijea - P(ij) Tie Tja) */
  if (!params.connect_xi) {
    global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 20, 2, 20, 0, "GIJKA");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tauIJAB");
    global_dpd_->file2_init(&I1A, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->contract244(&I1A, &Tau, &G, 1, 2, 1, -1.0, 1.0);
    global_dpd_->file2_close(&I1A);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 12, 30, 12, 30, 0, "Gijka");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tauijab");
    global_dpd_->file2_init(&I1B, PSIF_EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    global_dpd_->contract244(&I1B, &Tau, &G, 1, 2, 1, -1.0, 1.0);
    global_dpd_->file2_close(&I1B);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    global_dpd_->file2_init(&I1A, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->contract244(&I1A, &Tau, &G, 1, 2, 1, -1.0, 1.0);
    global_dpd_->file2_close(&I1A);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tauiJaB");
    global_dpd_->file2_init(&I1B, PSIF_EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    global_dpd_->contract244(&I1B, &Tau, &G, 1, 2, 1, -1.0, 1.0);
    global_dpd_->file2_close(&I1B);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&G);
  }

  x_Gijka_6_uhf();
  x_Gijka_7_uhf();

  /* term 8, +P(ij) Lkmfe rimae tjf */
  /* term 9, +P(ij) Lkmfe Timae Rjf, uses Z3, Z4 */
  x_Gijka_8_uhf();

  /* term 10, +P(IJ) LKMEF RJF TMA TIE */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 2, 0, 2, 0, "Z5(JI,KM)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 2, 20, 2, 20, 0, "L2R1_OOVO(pqsr)");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&V, &T1A, &Z, 3, 1, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 0, 0, 2, 0, "Z5(JI,KM)");
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 20, 0, 20, 0, "Z5(JI,KA)");
  global_dpd_->contract424(&Z, &T1A, &Z2, 3, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 20, 2, 20, 0, "GIJKA");
  global_dpd_->buf4_axpy(&Z2, &G, -1.0);
  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, qprs, 0, 20, "Z5(IJ,KA)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 20, 0, 20, 0, "Z5(IJ,KA)");
  global_dpd_->buf4_axpy(&Z2, &G, 1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&G);
  /* term 10, +P(ij) lkmef rjf tma tie */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 12, 10, 12, 0, "Z5(ji,km)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 12, 30, 12, 30, 0, "L2R1_oovo(pqsr)");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&V, &T1B, &Z, 3, 1, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 12, 0, "Z5(ji,km)");
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 10, 30, 10, 30, 0, "Z5(ji,ka)");
  global_dpd_->contract424(&Z, &T1B, &Z2, 3, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 30, 12, 30, 0, "Gijka");
  global_dpd_->buf4_axpy(&Z2, &G, -1.0);
  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, qprs, 10, 30, "Z5(ij,ka)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 10, 30, 10, 30, 0, "Z5(ij,ka)");
  global_dpd_->buf4_axpy(&Z2, &G, 1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&G);
  /* term 10, GIjKa += LKmEf Rjf TIE tma */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 22, 22, 22, 22, 0, "Z(Ij,Km)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 22, 26, 22, 26, 0, "L2R1_OoVo");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1A, &V, &Z, 1, 2, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 22, 24, 22, 24, 0, "L2R1_OovO(pqsr)");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&V, &T1B, &Z, 3, 1, 1, 1.0, 1.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&Z, &T1B, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  /* term 10, GiJkA += P(ij) LkMeF RJF Tie tMA */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 23, 23, 23, 23, 0, "Z(iJ,kM)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 23, 25, 23, 25, 0, "L2R1_OovO(qprs)");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract244(&T1B, &V, &Z, 1, 2, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 23, 27, 23, 27, 0, "L2R1_OoVo(qpsr)");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&V, &T1A, &Z, 3, 1, 1, 1.0, 1.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&Z, &T1A, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  psio_close(PSIF_EOM_TMP1, 0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);

  /* add 1/2 to ground-state parts of density */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 20, 2, 20, 0, "GIJKA");
  global_dpd_->buf4_init(&V, PSIF_CC_GAMMA, G_irr, 2, 20, 2, 20, 0, "GIJKA");
  global_dpd_->buf4_axpy(&G, &V, 0.5);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 12, 30, 12, 30, 0, "Gijka");
  global_dpd_->buf4_init(&V, PSIF_CC_GAMMA, G_irr, 12, 30, 12, 30, 0, "Gijka");
  global_dpd_->buf4_axpy(&G, &V, 0.5);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  global_dpd_->buf4_init(&V, PSIF_CC_GAMMA, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  global_dpd_->buf4_axpy(&G, &V, 0.5);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  global_dpd_->buf4_init(&V, PSIF_CC_GAMMA, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  global_dpd_->buf4_axpy(&G, &V, 0.5);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&G);

  /* clear out temporary files */
  psio_close(PSIF_EOM_TMP0, 0);
  psio_open(PSIF_EOM_TMP0, PSIO_OPEN_NEW);

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
  global_dpd_->file2_init(&LRA, PSIF_EOM_TMP, G_irr, 0, 0, "LR_OO");
  global_dpd_->file2_init(&LRB, PSIF_EOM_TMP, G_irr, 2, 2, "LR_oo");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->file2_mat_init(&T1A);   global_dpd_->file2_mat_init(&T1B);
  global_dpd_->file2_mat_init(&LRA);  global_dpd_->file2_mat_init(&LRB);
  global_dpd_->file2_mat_rd(&T1A);     global_dpd_->file2_mat_rd(&T1B);
  global_dpd_->file2_mat_rd(&LRA);    global_dpd_->file2_mat_rd(&LRB);

  /* rho_IJKA += - LR_OO(K,I) T(J,A) + LR_OO(K,J) T(I,A) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 20, 2, 20, 0, "GIJKA");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
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
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  /* rho_ijka += - LR_oo(k,i) T(j,a) + LR_oo(k,j) T(i,a) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 30, 12, 30, 0, "Gijka");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
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
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  /* rho_IjKa += - LR_OO(K,I) T(j,a) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
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
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  /* rho_iJkA += - LR_oo(k,i) T(J,A) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
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
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);
  global_dpd_->file2_mat_close(&LRA);
  global_dpd_->file2_mat_close(&LRB);
  global_dpd_->file2_close(&LRA);
  global_dpd_->file2_close(&LRB);

  global_dpd_->file2_mat_close(&T1A);
  global_dpd_->file2_mat_close(&T1B);
  global_dpd_->file2_close(&T1A);
  global_dpd_->file2_close(&T1B);

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
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->file2_init(&LTA, PSIF_EOM_TMP, L_irr, 0, 0, "LT_OO");
  global_dpd_->file2_init(&LTB, PSIF_EOM_TMP, L_irr, 2, 2, "LT_oo");
  global_dpd_->file2_mat_init(&R1A);
  global_dpd_->file2_mat_init(&R1B);
  global_dpd_->file2_mat_init(&LTA);
  global_dpd_->file2_mat_init(&LTB);
  global_dpd_->file2_mat_rd(&R1A);
  global_dpd_->file2_mat_rd(&R1B);
  global_dpd_->file2_mat_rd(&LTA);
  global_dpd_->file2_mat_rd(&LTB);

  /* rho_IJKA += - LT_OO(K,I)   R(J,A) + LT_OO(K,J)   R(I,A) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 20, 2, 20, 0, "GIJKA");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
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
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  /* rho_ijka += - LT_oo(k,i)   R(j,a) + LT_oo(k,j)   R(i,a) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 30, 12, 30, 0, "Gijka");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
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
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  /* rho_IjKa += - LT_OO(K,I)   R(j,a) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
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
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  /* rho_iJkA += - LT_oo(k,i)   R(J,A) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
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
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&LTA);
  global_dpd_->file2_mat_close(&LTB);
  global_dpd_->file2_mat_close(&R1A);
  global_dpd_->file2_mat_close(&R1B);

  global_dpd_->file2_close(&LTA);
  global_dpd_->file2_close(&LTB);
  global_dpd_->file2_close(&R1A);
  global_dpd_->file2_close(&R1B);

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
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 20, 0, 20, 0, 0, "Z(IA,KJ)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 20, 20, 20, 20, 0, "R2L2_OVOV");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 0, 20, "Z(IJ,KA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 20, 0, 20, 0, "Z(IJ,KA)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 20, 2, 20, 0, "GIJKA");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 0, 20, "Z(JI,KA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 20, 0, 20, 0, "Z(JI,KA)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 30, 10, 30, 10, 0, "Z(ia,kj)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 30, 30, 30, 30, 0, "R2L2_ovov");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 10, 30, "Z(ij,ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 30, 10, 30, 0, "Z(ij,ka)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 30, 12, 30, 0, "Gijka");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 10, 30, "Z(ji,ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 30, 10, 30, 0, "Z(ji,ka)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);

  /* GIjKa += R2L2_OvOv(Ia,Kf) T(j,f) */
  /* GIjKa -= R2L2_OvOv(ja,KF) T(I,F) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 24, 22, 24, 22, 0, "Z(Ia,Kj)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 24, 24, 24, 24, 0, "R2L2_OvOv");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 22, 24, "Z(Ij,Ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 22, 24, 22, 24, 0, "Z(Ij,Ka)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 30, 0, 30, 0, 0, "Z2(ja,KI)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 30, 20, 30, 20, 0, "R2L2_ovOV");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 23, 24, "Z2(jI,Ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 23, 24, 23, 24, 0, "Z2(jI,Ka)");
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 22, 24, "Z2(Ij,Ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 22, 24, 22, 24, 0, "Z2(Ij,Ka)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  /* GiJkA += R2L2_oVoV(iA,kF) T(J,F) */
  /* GiJkA += R2L2_OVov(JA,kf) T(i,f) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 27, 23, 27, 23, 0, "Z(iA,kJ)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 27, 27, 27, 27, 0, "R2L2_oVoV");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 23, 27, "Z(iJ,kA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 23, 27, 23, 27, 0, "Z(iJ,kA)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 20, 10, 20, 10, 0, "Z2(JA,ki)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 20, 30, 20, 30, 0, "R2L2_OVov");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 22, 27, "Z2(Ji,kA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 22, 27, 22, 27, 0, "Z2(Ji,kA)");
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 23, 27, "Z2(iJ,kA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 23, 27, 23, 27, 0, "Z2(iJ,kA)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  psio_close(PSIF_EOM_TMP1, 0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);

  /* term 9, +P(ij) Lkmfe Timae Rjf */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 20, 0, 20, 0, 0, "Z3(IA,KJ)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 20, 20, 20, 20, 0, "VIAJB");
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&V, &R1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 0, 20, "Z3(IJ,KA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 20, 0, 20, 0, "Z3(IJ,KA)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 20, 2, 20, 0, "GIJKA");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 0, 20, "Z3(JI,KA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 20, 0, 20, 0, "Z3(JI,KA)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 30, 10, 30, 10, 0, "Z3(ia,kj)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 30, 30, 30, 30, 0, "Viajb");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract424(&V, &R1B, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1B);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 10, 30, "Z3(ij,ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 30, 10, 30, 0, "Z3(ij,ka)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 30, 12, 30, 0, "Gijka");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 10, 30, "Z3(ji,ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 30, 10, 30, 0, "Z3(ji,ka)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);

  /* GIjKa += L2T2_OvOv(Ia,Kf) R(j,f) */
  /* GIjKa -= L2T2_OvOv(ja,KF) R(I,F) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 24, 22, 24, 22, 0, "Z3(Ia,Kj)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 24, 24, 24, 24, 0, "VIaJb");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract424(&V, &R1B, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1B);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 22, 24, "Z3(Ij,Ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 22, 24, 22, 24, 0, "Z3(Ij,Ka)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 30, 0, 30, 0, 0, "Z4(ja,KI)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 30, 20, 30, 20, 0, "ViaJB");
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&V, &R1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 23, 24, "Z4(jI,Ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 23, 24, 23, 24, 0, "Z4(jI,Ka)");
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 22, 24, "Z4(Ij,Ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 22, 24, 22, 24, 0, "Z4(Ij,Ka)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 24, 22, 24, 0, "GIjKa");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  /* GiJkA += L2T2_oVoV(iA,kF) R(J,F) */
  /* GiJkA += L2T2_OVov(JA,kf) R(i,f) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 27, 23, 27, 23, 0, "Z3(iA,kJ)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 27, 27, 27, 27, 0, "ViAjB");
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&V, &R1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 23, 27, "Z3(iJ,kA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 23, 27, 23, 27, 0, "Z3(iJ,kA)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 20, 10, 20, 10, 0, "Z4(JA,ki)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 20, 30, 20, 30, 0, "VIAjb");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract424(&V, &R1B, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1B);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 22, 27, "Z4(Ji,kA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 22, 27, 22, 27, 0, "Z4(Ji,kA)");
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 23, 27, "Z4(iJ,kA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 23, 27, 23, 27, 0, "Z4(iJ,kA)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 23, 27, 23, 27, 0, "GiJkA");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  return;
}

}} // namespace psi::ccdensity
