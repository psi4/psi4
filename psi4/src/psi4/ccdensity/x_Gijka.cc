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

void x_Gijka_rohf(void);
void x_Gijka_6_rohf(void);
void x_Gijka_7_rohf(void);
void x_Gijka_8_rohf(void);
extern void x_Gijka_uhf(void);

void x_Gijka(void) {
  if (params.ref == 0 || params.ref == 1)
    x_Gijka_rohf();
  else
    x_Gijka_uhf();
  return;
}

/* This function computes the non-R0 parts of the 2pdm density matrix
   Gijka = 0.5 *(rho_kaij + rho_Gijka) */

void x_Gijka_rohf(void) {
  int h, nirreps, i, j, k, a, I, J, K, A, Isym, Jsym, Ksym, Asym, row, col;
  int II,JJ,IIsym,JJsym;
  int L_irr, R_irr, G_irr;
  double value;
  dpdfile2 L1A, T1A, L1B, T1B, R1A, R1B, I1A, I1B;
  dpdbuf4 G, V, T, L, Z, Z1, Z2, Tau;

  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* term 1, rho_kaij += Lijae * Rke */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 2, 10, 2, 10, 0, "L2R1_OOVO(pqsr)");
  global_dpd_->buf4_copy(&Z, PSIF_EOM_TMP0, "GIJKA");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 2, 10, 2, 10, 0, "L2R1_oovo(pqsr)");
  global_dpd_->buf4_copy(&Z, PSIF_EOM_TMP0, "Gijka");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 0, 10, 0, 10, 0, "L2R1_OovO(pqsr)");
  global_dpd_->buf4_copy(&Z, PSIF_EOM_TMP0, "GIjKa");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 0, 10, 0, 10, 0, "L2R1_OoVo(qpsr)");
  global_dpd_->buf4_copy(&Z, PSIF_EOM_TMP0, "GiJkA");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  global_dpd_->buf4_scm(&Z, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  global_dpd_->buf4_scm(&Z, -1.0);
  global_dpd_->buf4_close(&Z);

  /* term 2, rho_ijka += Rijae * Lke */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "GIJKA");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 2, 10, 2, 10, 0, "L1R2_OOVO(pqsr)");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "Gijka");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 2, 10, 2, 10, 0, "L1R2_oovo(pqsr)");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 0, 10, 0, 10, 0, "L1R2_OovO(pqsr)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 0, 10, 0, 10, 0, "L1R2_OoVo(qpsr)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  /* term 3, rho_ijka += 0.5 Rijef Lkmef tma  */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "GIJKA");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 2, 0, 2, 2, 0, "R2L2_OOOO");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&Z, &T1A, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "Gijka");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 2, 0, 2, 2, 0, "R2L2_oooo");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract424(&Z, &T1B, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 0, 0, 0, 0, 0, "R2L2_OoOo");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract424(&Z, &T1B, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 0, 0, 0, 0, 0, "R2L2_oOoO");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&Z, &T1A, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  /* term 4, rho_ijka += 0.5 [Tijef + P(ij) Tie Tjf] Lkmef Rma */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "GIJKA");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, L_irr, 2, 0, 2, 2, 0, "Tau2L2_OOOO");
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&Z, &R1A, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&R1A);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "Gijka");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, L_irr, 2, 0, 2, 2, 0, "Tau2L2_oooo");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->contract424(&Z, &R1B, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&R1B);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, L_irr, 0, 0, 0, 0, 0, "Tau2L2_OoOo");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->contract424(&Z, &R1B, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&R1B);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, L_irr, 0, 0, 0, 0, 0, "Tau2L2_oOoO");
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&Z, &R1A, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&R1A);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  /* term 5, rho_ijka += - (Lkmef Rmf) (Tijea - P(ij) Tie Tja) */
  if (!params.connect_xi) {
    global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "GIJKA");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tauIJAB");
    global_dpd_->file2_init(&I1A, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->contract244(&I1A, &Tau, &G, 1, 2, 1, -1.0, 1.0);
    global_dpd_->file2_close(&I1A);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "Gijka");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tauijab");
    global_dpd_->file2_init(&I1B, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    global_dpd_->contract244(&I1B, &Tau, &G, 1, 2, 1, -1.0, 1.0);
    global_dpd_->file2_close(&I1B);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->file2_init(&I1A, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->contract244(&I1A, &Tau, &G, 1, 2, 1, -1.0, 1.0);
    global_dpd_->file2_close(&I1A);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
    global_dpd_->buf4_init(&Tau, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauiJaB");
    global_dpd_->file2_init(&I1B, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    global_dpd_->contract244(&I1B, &Tau, &G, 1, 2, 1, -1.0, 1.0);
    global_dpd_->file2_close(&I1B);
    global_dpd_->buf4_close(&Tau);
    global_dpd_->buf4_close(&G);
  }

  x_Gijka_6_rohf();
  x_Gijka_7_rohf();

  /* term 8, +P(ij) Lkmfe rimae tjf */
  /* term 9, +P(ij) Lkmfe Timae Rjf, uses Z3, Z4 */
  x_Gijka_8_rohf();

  /* term 10, +P(IJ) LKMEF RJF TMA TIE */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 2, 0, 2, 0, "Z5(JI,KM)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 2, 10, 2, 10, 0, "L2R1_OOVO(pqsr)");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&V, &T1A, &Z, 3, 1, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 0, 0, 2, 0, "Z5(JI,KM)");
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z5(JI,KA)");
  global_dpd_->contract424(&Z, &T1A, &Z2, 3, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 2, 10, 0, "GIJKA");
  global_dpd_->buf4_axpy(&Z2, &G, -1.0);
  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, qprs, 0, 10, "Z5(IJ,KA)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z5(IJ,KA)");
  global_dpd_->buf4_axpy(&Z2, &G, 1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&G);
  /* term 10, +P(ij) lkmef rjf tma tie */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 2, 0, 2, 0, "Z5(ji,km)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 2, 10, 2, 10, 0, "L2R1_oovo(pqsr)");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract424(&V, &T1B, &Z, 3, 1, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 0, 0, 2, 0, "Z5(ji,km)");
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z5(ji,ka)");
  global_dpd_->contract424(&Z, &T1B, &Z2, 3, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 2, 10, 0, "Gijka");
  global_dpd_->buf4_axpy(&Z2, &G, -1.0);
  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, qprs, 0, 10, "Z5(ij,ka)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z5(ij,ka)");
  global_dpd_->buf4_axpy(&Z2, &G, 1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&G);
  /* term 10, GIjKa += LKmEf Rjf TIE tma */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 0, 0, 0, 0, "Z(Ij,Km)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OoVo");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1A, &V, &Z, 1, 2, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 0, 10, 0, 10, 0, "L2R1_OovO(pqsr)");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract424(&V, &T1B, &Z, 3, 1, 1, 1.0, 1.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract424(&Z, &T1B, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  /* term 10, GiJkA += P(ij) LkMeF RJF Tie tMA */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 0, 0, 0, 0, "Z(iJ,kM)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OovO(qprs)");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract244(&T1B, &V, &Z, 1, 2, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 0, 10, 0, 10, 0, "L2R1_OoVo(qpsr)");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&V, &T1A, &Z, 3, 1, 1, 1.0, 1.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&Z, &T1A, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  psio_close(PSIF_EOM_TMP1, 0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);

  /* add 1/2 to ground-state parts of density */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "GIJKA");
  global_dpd_->buf4_init(&V, PSIF_CC_GAMMA, G_irr, 2, 10, 2, 10, 0, "GIJKA");
  global_dpd_->buf4_axpy(&G, &V, 0.5);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 10, 2, 10, 0, "Gijka");
  global_dpd_->buf4_init(&V, PSIF_CC_GAMMA, G_irr, 2, 10, 2, 10, 0, "Gijka");
  global_dpd_->buf4_axpy(&G, &V, 0.5);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  global_dpd_->buf4_init(&V, PSIF_CC_GAMMA, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  global_dpd_->buf4_axpy(&G, &V, 0.5);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  global_dpd_->buf4_init(&V, PSIF_CC_GAMMA, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  global_dpd_->buf4_axpy(&G, &V, 0.5);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&G);

  /* clear out temporary files */
  psio_close(PSIF_EOM_TMP0, 0);
  psio_open(PSIF_EOM_TMP0, PSIO_OPEN_NEW);

  return;
}



/* This function computes term 6,
   rho_ijka -= P(ij) lke rie tja or
   rho_ijka -= P(ij) LR1_OO(k,i) T(j,a) */

void x_Gijka_6_rohf(void) {
  int h, nirreps, i, j, k, a, I, J, K, A, Isym, Jsym, Ksym, Asym, row, col;
  int II,JJ,IIsym,JJsym;
  int L_irr, R_irr, G_irr;
  dpdfile2 LR1A, LR1B, T1A, T1B;
  dpdbuf4 G;

  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* open one-electron files for the nasty terms */
  global_dpd_->file2_init(&LR1A, PSIF_EOM_TMP, G_irr, 0, 0, "LR_OO");
  global_dpd_->file2_init(&LR1B, PSIF_EOM_TMP, G_irr, 0, 0, "LR_oo");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->file2_mat_init(&T1A);   global_dpd_->file2_mat_init(&T1B);
  global_dpd_->file2_mat_init(&LR1A);   global_dpd_->file2_mat_init(&LR1B);
  global_dpd_->file2_mat_rd(&T1A);     global_dpd_->file2_mat_rd(&T1B);
  global_dpd_->file2_mat_rd(&LR1A);     global_dpd_->file2_mat_rd(&LR1B);

  /* rho_IJKA += - LR1_OO(K,I) T(J,A) + LR1_OO(K,J) T(I,A) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 2, 10, 0, "GIJKA");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LR1A.params->colidx[i]; Isym = LR1A.params->qsym[i];
      J = T1A.params->rowidx[j]; Jsym = T1A.params->psym[j];
      II = T1A.params->rowidx[i]; IIsym = T1A.params->psym[i];
      JJ = LR1A.params->colidx[j]; JJsym = LR1A.params->qsym[j];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        k = G.params->colorb[h^G_irr][col][0];
        a = G.params->colorb[h^G_irr][col][1];
        K = LR1A.params->rowidx[k]; Ksym = LR1A.params->psym[k];
        A = T1A.params->colidx[a]; Asym = T1A.params->qsym[a];
        if( ((Ksym^Isym)==G_irr) && (Jsym==Asym))
          G.matrix[h][row][col] -=
            LR1A.matrix[Ksym][K][I] * T1A.matrix[Jsym][J][A];
        if( ((Ksym^JJsym)==G_irr) && (IIsym==Asym))
          G.matrix[h][row][col] +=
            LR1A.matrix[Ksym][K][JJ] * T1A.matrix[IIsym][II][A];
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  /* rho_ijka += - LR1_oo(k,i) T(j,a) + LR1_oo(k,j) T(i,a) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 2, 10, 0, "Gijka");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LR1B.params->colidx[i]; Isym = LR1B.params->qsym[i];
      J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
      II = T1B.params->rowidx[i]; IIsym = T1B.params->psym[i];
      JJ = LR1B.params->colidx[j]; JJsym = LR1B.params->qsym[j];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        k = G.params->colorb[h^G_irr][col][0];
        a = G.params->colorb[h^G_irr][col][1];
        K = LR1B.params->rowidx[k]; Ksym = LR1B.params->psym[k];
        A = T1B.params->colidx[a]; Asym = T1B.params->qsym[a];
        if( ((Ksym^Isym)==G_irr) && (Jsym==Asym))
          G.matrix[h][row][col] -=
            LR1B.matrix[Ksym][K][I] * T1B.matrix[Jsym][J][A];
        if( ((Ksym^JJsym)==G_irr) && (IIsym==Asym))
          G.matrix[h][row][col] +=
            LR1B.matrix[Ksym][K][JJ] * T1B.matrix[IIsym][II][A];
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  /* rho_IjKa += - LR1_OO(K,I) T(j,a) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LR1A.params->colidx[i]; Isym = LR1A.params->qsym[i];
      J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        k = G.params->colorb[h^G_irr][col][0];
        a = G.params->colorb[h^G_irr][col][1];
        K = LR1A.params->rowidx[k]; Ksym = LR1A.params->psym[k];
        A = T1B.params->colidx[a]; Asym = T1B.params->qsym[a];
        if( ((Ksym^Isym)==G_irr) && (Jsym==Asym))
          G.matrix[h][row][col] -=
            LR1A.matrix[Ksym][K][I] * T1B.matrix[Jsym][J][A];
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  /* rho_iJkA += - LR1_oo(k,i) T(J,A) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I = LR1B.params->colidx[i]; Isym = LR1B.params->qsym[i];
      J = T1A.params->rowidx[j]; Jsym = T1A.params->psym[j];
      for(col=0; col < G.params->coltot[h^G_irr]; col++) {
        k = G.params->colorb[h^G_irr][col][0];
        a = G.params->colorb[h^G_irr][col][1];
        K = LR1B.params->rowidx[k]; Ksym = LR1B.params->psym[k];
        A = T1A.params->colidx[a]; Asym = T1A.params->qsym[a];
        if( ((Ksym^Isym)==G_irr) && (Jsym==Asym))
          G.matrix[h][row][col] -=
            LR1B.matrix[Ksym][K][I] * T1A.matrix[Jsym][J][A];
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);
  global_dpd_->file2_mat_close(&LR1A);
  global_dpd_->file2_mat_close(&LR1B);
  global_dpd_->file2_close(&LR1A);
  global_dpd_->file2_close(&LR1B);

  global_dpd_->file2_mat_close(&T1A);
  global_dpd_->file2_mat_close(&T1B);
  global_dpd_->file2_close(&T1A);
  global_dpd_->file2_close(&T1B);

  return;
}



/* This function computes
   Gijka -= P(ij) LT_OO(k,i) * R(j,a) */

void x_Gijka_7_rohf(void) {
  int h, nirreps, i, j, k, a, I, J, K, A, Isym, Jsym, Ksym, Asym, row, col;
  int II,JJ,IIsym,JJsym;
  int L_irr, R_irr, G_irr;
  dpdfile2 R1A, R1B, LTA, LTB;
  dpdbuf4 G;

  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* open one-electron files for the nasty terms */
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->file2_init(&LTA, PSIF_EOM_TMP, L_irr, 0, 0, "LT_OO");
  global_dpd_->file2_init(&LTB, PSIF_EOM_TMP, L_irr, 0, 0, "LT_oo");
  global_dpd_->file2_mat_init(&R1A);
  global_dpd_->file2_mat_init(&R1B);
  global_dpd_->file2_mat_init(&LTA);
  global_dpd_->file2_mat_init(&LTB);
  global_dpd_->file2_mat_rd(&R1A);
  global_dpd_->file2_mat_rd(&R1B);
  global_dpd_->file2_mat_rd(&LTA);
  global_dpd_->file2_mat_rd(&LTB);

  /* rho_IJKA += - LT_OO(K,I)   R(J,A) + LT_OO(K,J)   R(I,A) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 2, 10, 0, "GIJKA");
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
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 2, 10, 0, "Gijka");
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
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
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
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
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
/* term 8, +P(ij) Lkmfe rimae tjf */
/* term 9, +P(ij) Lkmfe rimae tjf */

void x_Gijka_8_rohf(void) {
  int h, nirreps, i, j, k, a, I, J, K, A, Isym, Jsym, Ksym, Asym, row, col;
  int II,JJ,IIsym,JJsym;
  int L_irr, R_irr, G_irr;
  double value;
  dpdfile2 L1A, T1A, L1B, T1B, R1A, R1B, I1A, I1B;
  dpdbuf4 G, V, T, L, Z, Z1, Z2, Tau;


  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* term 8, +P(ij) Lkmfe rimae tjf */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z(IA,KJ)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 0, 10, "Z(IJ,KA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z(IJ,KA)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 2, 10, 0, "GIJKA");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 0, 10, "Z(JI,KA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z(JI,KA)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z(ia,kj)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovov");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 0, 10, "Z(ij,ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z(ij,ka)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 2, 10, 0, "Gijka");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 0, 10, "Z(ji,ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z(ji,ka)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);

  /* GIjKa += R2L2_OvOv(Ia,Kf) T(j,f) */
  /* GIjKa -= R2L2_OvOv(ja,KF) T(I,F) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z(Ia,Kj)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 0, 10, "Z(Ij,Ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z(Ij,Ka)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z2(ja,KI)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovOV");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 0, 10, "Z2(jI,Ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z2(jI,Ka)");
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 0, 10, "Z2(Ij,Ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z2(Ij,Ka)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  /* GiJkA += R2L2_oVoV(iA,kF) T(J,F) */
  /* GiJkA += R2L2_OVov(JA,kf) T(i,f) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z(iA,kJ)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_oVoV");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 0, 10, "Z(iJ,kA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z(iJ,kA)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z2(JA,ki)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tia");
  global_dpd_->contract424(&V, &T1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 0, 10, "Z2(Ji,kA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z2(Ji,kA)");
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 0, 10, "Z2(iJ,kA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z2(iJ,kA)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  psio_close(PSIF_EOM_TMP1, 0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);

  /* term 9, +P(ij) Lkmfe rimae tjf */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z3(IA,KJ)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIAJB");
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&V, &R1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 0, 10, "Z3(IJ,KA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z3(IJ,KA)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 2, 10, 0, "GIJKA");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 0, 10, "Z3(JI,KA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z3(JI,KA)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z3(ia,kj)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 10, 10, 10, 10, 0, "Viajb");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->contract424(&V, &R1B, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1B);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 0, 10, "Z3(ij,ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z3(ij,ka)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 2, 10, 0, "Gijka");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 0, 10, "Z3(ji,ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z3(ji,ka)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_close(&Z);

  /* GIjKa += R2L2_OvOv(Ia,Kf) T(j,f) */
  /* GIjKa -= R2L2_OvOv(ja,KF) T(I,F) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z3(Ia,Kj)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIaJb");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->contract424(&V, &R1B, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1B);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 0, 10, "Z3(Ij,Ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z3(Ij,Ka)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z4(ja,KI)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 10, 10, 10, 10, 0, "ViaJB");
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&V, &R1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 0, 10, "Z4(jI,Ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z4(jI,Ka)");
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 0, 10, "Z4(Ij,Ka)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z4(Ij,Ka)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GIjKa");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  /* GiJkA += R2L2_oVoV(iA,kF) T(J,F) */
  /* GiJkA += R2L2_OVov(JA,kf) T(i,f) */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z3(iA,kJ)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 10, 10, 10, 10, 0, "ViAjB");
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&V, &R1A, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1A);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 0, 10, "Z3(iJ,kA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z3(iJ,kA)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 0, 10, 0, 0, "Z4(JA,ki)");
  global_dpd_->buf4_init(&V, PSIF_EOM_TMP, L_irr, 10, 10, 10, 10, 0, "VIAjb");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->contract424(&V, &R1B, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1B);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 0, 10, "Z4(Ji,kA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z4(Ji,kA)");
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 0, 10, "Z4(iJ,kA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 10, 0, 10, 0, "Z4(iJ,kA)");
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 10, 0, 10, 0, "GiJkA");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  return;
}

}} // namespace psi::ccdensity
