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

/* computes non R0 parts of EOM CCSD Gijab */
void x_Gijab_uhf_2(void);
void x_Gijab_uhf_3(void);

void x_Gijab_uhf(void)
{
  int h, nirreps, II;
  int R_irr, L_irr, G_irr;
  double value, tval;
  dpdfile2 T1, L1, I1, T1A, T1B, Z1, R1;
  dpdbuf4 R, I, G, L, T, V, Z, Z2;

  nirreps = moinfo.nirreps;
  R_irr = params.R_irr; L_irr = params.L_irr; G_irr = params.G_irr;

  /* term 1,2: (L*R) * Tau(IJ,AB), see comments in xi1_connected */
  if (G_irr == 0) {
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->buf4_copy(&T, PSIF_EOM_TMP0, "GIJAB");
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    global_dpd_->buf4_copy(&T, PSIF_EOM_TMP0, "Gijab");
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    global_dpd_->buf4_copy(&T, PSIF_EOM_TMP0, "GIjAb");
    global_dpd_->buf4_close(&T);
  }
  else {
    global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 7, 2, 7, 0, "GIJAB");
    global_dpd_->buf4_scm(&G,0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 12, 17, 12, 17, 0, "Gijab");
    global_dpd_->buf4_scm(&G,0.0);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 28, 22, 28, 0, "GIjAb");
    global_dpd_->buf4_scm(&G,0.0);
    global_dpd_->buf4_close(&G);
  }

  /* -P(ij) LR_OO(M,I) Tau(MJ,AB); terms 4,5,8,9 */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z(IJ,A>B)");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 0, "LR_OO");
  global_dpd_->file2_init(&Z1, PSIF_EOM_TMP1, G_irr, 0, 0, "Z(N,I)");
  global_dpd_->file2_axpy(&I1, &Z1, 1.0, 0);
  global_dpd_->file2_close(&I1);
  /* -P(ij) L2R1_OV(M,F) T(I,F) Tau(MJ,AB); terms 35, 36 */
  if (!params.connect_xi) {
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&I1, &T1, &Z1, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&I1);
  }
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tauIJAB");
  global_dpd_->contract244(&Z1, &T, &Z, 0, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&Z1);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 7, 2, 7, 0, "GIJAB");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 0, 7, "Z(JI,A>B)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z(JI,A>B)");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  /* -P(ij) LR_oo(m,i) Tau(mj,ab); terms 4,5,8,9 */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 17, 10, 17, 0, "Z(ij,a>b)");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 2, 2, "LR_oo");
  global_dpd_->file2_init(&Z1, PSIF_EOM_TMP1, G_irr, 2, 2, "Z(n,i)");
  global_dpd_->file2_axpy(&I1, &Z1, 1.0, 0);
  global_dpd_->file2_close(&I1);
  /* -P(ij) L2R1_ov(m,f) T(i,f) Tau(mj,ab); terms 35, 36 */
  if (!params.connect_xi) {
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract222(&I1, &T1, &Z1, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&I1);
  }
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tauijab");
  global_dpd_->contract244(&Z1, &T, &Z, 0, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&Z1);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 17, 12, 17, 0, "Gijab");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 10, 17, "Z(ji,a>b)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 17, 10, 17, 0, "Z(ji,a>b)");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  psio_close(PSIF_EOM_TMP1, 0);
  psio_open(PSIF_EOM_TMP1,PSIO_OPEN_NEW);

  /* GIjAb += -P(Ij) LR_OO(M,I) Tau(Mj,Ab); terms 4,5,8,9 */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 28, 22, 28, 0, "GIjAb");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 0, "LR_OO");
  global_dpd_->file2_init(&Z1, PSIF_EOM_TMP1, G_irr, 0, 0, "Z(N,I)");
  global_dpd_->file2_axpy(&I1, &Z1, 1.0, 0);
  global_dpd_->file2_close(&I1);
  /* -P(Ij) L2R1_OV(M,F) T(I,F) Tau(Nm,Ab); terms 35, 36 */
  if (!params.connect_xi) {
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract222(&I1, &T1, &Z1, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&I1);
  }
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
  global_dpd_->contract244(&Z1, &T, &G, 0, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&Z1);
  global_dpd_->buf4_close(&T);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 23, 29, 23, 29, 0, "Z(jI,bA)");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 2, 2, "LR_oo");
  global_dpd_->file2_init(&Z1, PSIF_EOM_TMP1, G_irr, 2, 2, "Z(n,i)");
  global_dpd_->file2_axpy(&I1, &Z1, 1.0, 0);
  global_dpd_->file2_close(&I1);
  if (!params.connect_xi) {
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract222(&I1, &T1, &Z1, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_close(&I1);
  }
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tauiJaB");
  global_dpd_->contract244(&Z1, &T, &Z, 0, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&Z1);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qpsr, 22, 28, "Z(Ij,Ab)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 22, 28, 22, 28, 0, "Z(Ij,Ab)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&G);

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1,PSIO_OPEN_NEW);

  /* -P(ab) LR_VV(F,A) Tau(IJ,FB); terms 6,7,10,11 */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(I>J,AB)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tauIJAB");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 1, 1, "LR_VV");
  global_dpd_->contract244(&I1, &T, &Z, 0, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 5, 2, 7, 0, "GIJAB");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, pqsr, 2, 5, "Z(I>J,BA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(I>J,BA)");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  /* -P(ab) LR_vv(f,a) Tau(ij,fb); */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 12, 15, 12, 15, 0, "Z(i>j,ab)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tauijab");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 3, 3, "LR_vv");
  global_dpd_->contract244(&I1, &T, &Z, 0, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 12, 15, 12, 17, 0, "Gijab");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, pqsr, 12, 15, "Z(i>j,ba)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 12, 15, 12, 15, 0, "Z(i>j,ba)");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);
  /* GIjAb += - LR_VV(F,A) Tau(Ij,Fb) + LR_VV(f,a) Tau(Ij,fA) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 28, 22, 28, 0, "GIjAb");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 1, 1, "LR_VV");
  global_dpd_->contract244(&I1, &T, &G, 0, 2, 1, -1.0, 1.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 22, 29, 22, 29, 0, "Z(Ij,bA)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 29, 22, 29, 0, "tauIjbA");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 3, 3, "LR_vv");
  global_dpd_->contract244(&I1, &T, &Z, 0, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, pqsr, 22, 28, "Z(Ij,Ab)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 22, 28, 22, 28, 0, "Z(Ij,Ab)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&G);

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1,PSIO_OPEN_NEW);

  /* + 1/4 Lmnef Rmnab Tau_ijef, terms 13, 15 */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 7, 2, 7, 0, "GIJAB");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, L_irr, 2, 2, 2, 2, 0, "Tau2L2_OOOO");
  global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 2, 7, 2, 7, 0, "RIJAB");
  global_dpd_->contract444(&I, &R, &G, 0, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&R);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 12, 17, 12, 17, 0, "Gijab");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, L_irr, 12, 12, 12, 12, 0, "Tau2L2_oooo");
  global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 12, 17, 12, 17, 0, "Rijab");
  global_dpd_->contract444(&I, &R, &G, 0, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&R);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 28, 22, 28, 0, "GIjAb");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, L_irr, 22, 22, 22, 22, 0, "Tau2L2_OoOo");
  global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 22, 28, 22, 28, 0, "RIjAb");
  global_dpd_->contract444(&I, &R, &G, 0, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&R);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&G);

  /* - 0.5 P(ij) (Lmnfe Rie) Tjf (taumnab), terms 24, 26 */
  /* + 1/4 Lmnef Rijef Tau_mnab, terms 12, 14 */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 2, 0, 2, 0, "Z(IJ,M>N)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 2, 20, 2, 20, 0, "L2R1_OOVO(pqsr)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&I, &T1, &Z, 3, 1, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&I);
  /* add terms 12, 14 */
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 0, 2, 2, 2, 0, "R2L2_OOOO");
  global_dpd_->buf4_axpy(&I, &Z, -0.5);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z(IJ,A>B)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  global_dpd_->contract444(&Z, &T, &Z2, 0, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 7, 2, 7, 0, "GIJAB");
  global_dpd_->buf4_axpy(&Z2, &G, -1.0);
  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, qprs, 0, 7, "Z(JI,A>B)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z(JI,A>B)");
  global_dpd_->buf4_axpy(&Z2, &G, 1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 12, 10, 12, 0, "Z(ij,m>n)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 12, 30, 12, 30, 0, "L2R1_oovo(pqsr)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&I, &T1, &Z, 3, 1, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&I);
  /* add terms 12, 14 */
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 10, 12, 12, 12, 0, "R2L2_oooo");
  global_dpd_->buf4_axpy(&I, &Z, -0.5);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 10, 17, 10, 17, 0, "Z(ij,a>b)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
  global_dpd_->contract444(&Z, &T, &Z2, 0, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 17, 12, 17, 0, "Gijab");
  global_dpd_->buf4_axpy(&Z2, &G, -1.0);
  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, qprs, 10, 17, "Z(ji,a>b)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 10, 17, 10, 17, 0, "Z(ji,a>b)");
  global_dpd_->buf4_axpy(&Z2, &G, 1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 22, 22, 22, 22, 0, "Z(Ij,Mn)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 22, 24, 22, 24, 0, "L2R1_OovO(pqsr)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&I, &T1, &Z, 3, 1, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 22, 26, 22, 26, 0, "L2R1_OoVo");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &I, &Z, 1, 2, 0, 1.0, 1.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&I);
  /* add terms 12, 14 */
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 22, 22, 22, 22, 0, "R2L2_OoOo");
  global_dpd_->buf4_axpy(&I, &Z, 1.0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 28, 22, 28, 0, "GIjAb");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
  global_dpd_->contract444(&Z, &T, &G, 0, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1,PSIO_OPEN_NEW);

  /* + 0.5 P(AB) (tau_IJEF LMNEF) RMA TNB ; terms 25, 27 */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 2, 21, 2, 21, 0, "Z(I>J,AN)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, L_irr, 2, 0, 2, 2, 0, "Tau2L2_OOOO");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract244(&R1, &I, &Z, 0, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(I>J,AB)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&Z, &T1, &Z2, 3, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 5, 2, 7, 0, "GIJAB");
  global_dpd_->buf4_axpy(&Z2, &G, 1.0);
  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, pqsr, 2, 5, "Z(I>J,BA)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(I>J,BA)");
  global_dpd_->buf4_axpy(&Z2, &G, -1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&G);
  /* + 0.5 P(ab) (tau_ijef Lmnef) Rma Tnb ; terms 25, 27 */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 12, 31, 12, 31, 0, "Z(i>j,an)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, L_irr, 12, 10, 12, 12, 0, "Tau2L2_oooo");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract244(&R1, &I, &Z, 0, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 12, 15, 12, 15, 0, "Z(i>j,ab)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&Z, &T1, &Z2, 3, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 12, 15, 12, 17, 0, "Gijab");
  global_dpd_->buf4_axpy(&Z2, &G, 1.0);
  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, pqsr, 12, 15, "Z(i>j,ba)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 12, 15, 12, 15, 0, "Z(i>j,ba)");
  global_dpd_->buf4_axpy(&Z2, &G, -1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&G);
  /* + tau_IjEf LMnEf RMA Tnb ; terms 25, 27 */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 22, 26, 22, 26, 0, "Z(Ij,An)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, L_irr, 22, 22, 22, 22, 0, "Tau2L2_OoOo");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract244(&R1, &I, &Z, 0, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 28, 22, 28, 0, "GIjAb");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&Z, &T1, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  /* + tau_IjEf LNmEf Rmb TNA ; terms 25, 27 */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 22, 26, 22, 26, 0, "Z2(Ij,An)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, L_irr, 22, 22, 22, 22, 0, "Tau2L2_OoOo");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &I, &Z, 0, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&I);
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract424(&Z, &R1, &G, 3, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&G);

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1,PSIO_OPEN_NEW);

  /* terms combined to P(ij)P(ab) Z(i,a) T(j,b), 3,22,23,33 */
  x_Gijab_uhf_2();

  /* terms combined to P(ij)P(ab) Z(i,a) R(j,b), 18,32,34,19 */
  x_Gijab_uhf_3();

  /* -P(ij)(Lme Tie + 0.5 Lmnef Tinef) Rmjab, term 16, 30 */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z(IJ,A>B)");
  global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 0, 7, 2, 7, 0, "RIJAB");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, L_irr, 0, 0, "LT_OO");
  global_dpd_->contract244(&I1, &R, &Z, 0, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_close(&R);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 7, 2, 7, 0, "GIJAB");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 0, 7, "Z(JI,A>B)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z(JI,A>B)");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 17, 10, 17, 0, "Z(ij,a>b)");
  global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 10, 17, 12, 17, 0, "Rijab");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, L_irr, 2, 2, "LT_oo");
  global_dpd_->contract244(&I1, &R, &Z, 0, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_close(&R);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 17, 12, 17, 0, "Gijab");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 10, 17, "Z(ji,a>b)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 17, 10, 17, 0, "Z(ji,a>b)");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 28, 22, 28, 0, "GIjAb");
  global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 22, 28, 22, 28, 0, "RIjAb");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, L_irr, 0, 0, "LT_OO");
  global_dpd_->contract244(&I1, &R, &G, 0, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_close(&R);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 23, 28, 23, 28, 0, "Z(jI,Ab)");
  global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 23, 28, 23, 28, 0, "RiJAb");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, L_irr, 2, 2, "LT_oo");
  global_dpd_->contract244(&I1, &R, &Z, 0, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_close(&R);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 22, 28, "Z(Ij,Ab)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 22, 28, 22, 28, 0, "Z(Ij,Ab)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);

  /* -P(ab)(Lme Tmb + 0.5 Lmnfe Tmnfb) Rijae, term 17, 31 */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(I>J,AB)");
  global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 2, 5, 2, 7, 0, "RIJAB");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, L_irr, 1, 1, "LT_VV");
  global_dpd_->contract424(&R, &I1, &Z, 3, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_close(&R);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 5, 2, 7, 0, "GIJAB");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, pqsr, 2, 5, "Z(I>J,BA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(I>J,BA)");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 12, 15, 12, 15, 0, "Z(i>j,ab)");
  global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 12, 15, 12, 17, 0, "Rijab");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, L_irr, 3, 3, "LT_vv");
  global_dpd_->contract424(&R, &I1, &Z, 3, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_close(&R);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 12, 15, 12, 17, 0, "Gijab");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, pqsr, 12, 15, "Z(i>j,ba)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 12, 15, 12, 15, 0, "Z(i>j,ba)");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 28, 22, 28, 0, "GIjAb");
  global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 22, 28, 22, 28, 0, "RIjAb");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, L_irr, 3, 3, "LT_vv");
  global_dpd_->contract424(&R, &I1, &G, 3, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_close(&R);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 22, 29, 22, 29, 0, "Z(Ij,bA)");
  global_dpd_->buf4_init(&R, PSIF_CC_GR, R_irr, 22, 29, 22, 29, 0, "RIjaB");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, L_irr, 1, 1, "LT_VV");
  global_dpd_->contract424(&R, &I1, &Z, 3, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_close(&R);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, pqsr, 22, 28, "Z(Ij,Ab)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 22, 28, 22, 28, 0, "Z(Ij,Ab)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);

  /* -P(ab) lmnef rme tijfb tna, term 37 */
  if (!params.connect_xi) {
    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 2, 21, 2, 21, 0, "Z(I>J,AN)");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->contract424(&T, &I1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&I1);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z(I>J,AB)");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract424(&Z, &T1, &Z2, 3, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 5, 2, 7, 0, "GIJAB");
    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, 0, 2, 5, 2, 5, 0, "Z(I>J,AB)");
    global_dpd_->buf4_axpy(&Z, &G, -1.0);
    global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, pqsr, 2, 5, "Z(I>J,BA)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, 0, 2, 5, 2, 5, 0, "Z(I>J,BA)");
    global_dpd_->buf4_axpy(&Z, &G, 1.0);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 12, 31, 12, 31, 0, "Z(i>j,an)");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    global_dpd_->contract424(&T, &I1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&I1);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 12, 15, 12, 15, 0, "Z(i>j,ab)");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract424(&Z, &T1, &Z2, 3, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 12, 15, 12, 17, 0, "Gijab");
    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 12, 15, 12, 15, 0, "Z(i>j,ab)");
    global_dpd_->buf4_axpy(&Z, &G, -1.0);
    global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, pqsr, 12, 15, "Z(i>j,ba)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 12, 15, 12, 15, 0, "Z(i>j,ba)");
    global_dpd_->buf4_axpy(&Z, &G, 1.0);
    global_dpd_->buf4_close(&G);

    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 22, 26, 22, 26, 0, "Z(Ij,An)");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    global_dpd_->contract424(&T, &I1, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&I1);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 28, 22, 28, 0, "GIjAb");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->contract424(&Z, &T1, &G, 3, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 22, 24, 22, 24, 0, "Z(Ij,Nb)");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->contract244(&I1, &T, &Z, 1, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&I1);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 28, 22, 28, 0, "GIjAb");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract244(&T1, &Z, &G, 0, 2, 1, -1.0, 1.0);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&G);
  }

  /* compute Z(IA,JB) Z(ia,jb) and Z(IA,jb) for terms
     20, 28, 29, 21 then permute and add in */
  /* + P(ij) P(ab) (Rimae Lnmfe) Tnjfb, term 20 */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP0, G_irr, 20, 20, 20, 20, 0, "Z(IA,JB)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 20, 20, 20, 20, 0, "R2L2_OVOV");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  global_dpd_->contract444(&I, &T, &Z, 0, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 20, 30, 20, 30, 0, "R2L2_OVov");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
  global_dpd_->contract444(&I, &T, &Z, 0, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&T);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP0, G_irr, 30, 30, 30, 30, 0, "Z(ia,jb)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 30, 30, 30, 30, 0, "R2L2_ovov");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
  global_dpd_->contract444(&I, &T, &Z, 0, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 30, 20, 30, 20, 0, "R2L2_ovOV");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
  global_dpd_->contract444(&I, &T, &Z, 0, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&T);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP0, G_irr, 20, 30, 20, 30, 0, "Z(IA,jb)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 20, 20, 20, 20, 0, "R2L2_OVOV");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
  global_dpd_->contract444(&I, &T, &Z, 0, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 20, 30, 20, 30, 0, "R2L2_OVov");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
  global_dpd_->contract444(&I, &T, &Z, 0, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 30, 20, 30, 20, 0, "R2L2_ovOV");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  global_dpd_->contract444(&T, &I, &Z, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 30, 30, 30, 30, 0, "R2L2_ovov");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
  global_dpd_->contract444(&T, &I, &Z, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 24, 27, 24, 27, 0, "Z(Ib,jA)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 27, 27, 27, 27, 0, "R2L2_oVoV");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
  global_dpd_->contract444(&T, &I, &Z, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 24, 24, 24, 24, 0, "R2L2_OvOv");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
  global_dpd_->contract444(&I, &T, &Z, 0, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_TMP0, psrq, 20, 30, "Z(IA,jb)", 1.0);
  global_dpd_->buf4_close(&Z);

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1,PSIO_OPEN_NEW);

  /* - P(ij) P(ab) (Tjmbe Lnmfe) Tif Rna, term 28 */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 20, 0, 20, 0, 0, "Z(JB,NI)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, L_irr, 20, 20, 20, 20, 0, "VIAJB");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 21, 20, 21, 20, 0, "Z(AI,JB)");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract244(&R1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_TMP0, qprs, 20, 20, "Z(IA,JB)", -1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 30, 10, 30, 10, 0, "Z(jb,ni)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, L_irr, 30, 30, 30, 30, 0, "Viajb");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 31, 30, 31, 30, 0, "Z(ai,jb)");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract244(&R1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_TMP0, qprs, 30, 30, "Z(ia,jb)", -1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 30, 0, 30, 0, 0, "Z(jb,NI)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, L_irr, 30, 20, 30, 20, 0, "ViaJB");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 21, 30, 21, 30, 0, "Z(AI,jb)");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract244(&R1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_TMP0, qprs, 20, 30, "Z(IA,jb)", -1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 24, 22, 24, 22, 0, "Z(Ib,Nj)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, L_irr, 24, 24, 24, 24, 0, "VIaJb");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 24, 26, 24, 26, 0, "Z(Ib,Aj)");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract244(&R1, &Z, &Z2, 0, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_TMP0, prsq, 20, 30, "Z(IA,jb)", +1.0);
  global_dpd_->buf4_close(&Z2);

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1,PSIO_OPEN_NEW);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 27, 23, 27, 23, 0, "Z(jA,nI)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, L_irr, 27, 27, 27, 27, 0, "ViAjB");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 25, 27, 25, 27, 0, "Z(bI,jA)");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract244(&R1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, sqrp, 21, 30, "Z(AI,jb)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 21, 30, 21, 30, 0, "Z(AI,jb)");
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_TMP0, qprs, 20, 30, "Z(IA,jb)", +1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, L_irr, 20, 10, 20, 10, 0, "Z(IA,nj)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, L_irr, 20, 30, 20, 30, 0, "VIAjb");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 20, 31, 20, 31, 0, "Z(IA,bj)");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract244(&R1, &Z, &Z2, 0, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_TMP0, pqsr, 20, 30, "Z(IA,jb)", -1.0);
  global_dpd_->buf4_close(&Z2);

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1,PSIO_OPEN_NEW);

  /* - P(ij) P(ab) (Tjmbe Lnmfe) Rif Tna, term 29 */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 20, 0, 20, 0, 0, "Z(JB,NI)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, L_irr, 20, 20, 20, 20, 0, "VIAJB");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&I, &R1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 21, 20, 21, 20, 0, "Z(AI,JB)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_TMP0, qprs, 20, 20, "Z(IA,JB)", -1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 30, 10, 30, 10, 0, "Z(jb,ni)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, L_irr, 30, 30, 30, 30, 0, "Viajb");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract424(&I, &R1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 31, 30, 31, 30, 0, "Z(ai,jb)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_TMP0, qprs, 30, 30, "Z(ia,jb)", -1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 30, 0, 30, 0, 0, "Z(jb,NI)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, L_irr, 30, 20, 30, 20, 0, "ViaJB");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&I, &R1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 21, 30, 21, 30, 0, "Z(AI,jb)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_TMP0, qprs, 20, 30, "Z(IA,jb)", -1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 24, 22, 24, 22, 0, "Z(Ib,Nj)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, L_irr, 24, 24, 24, 24, 0, "VIaJb");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract424(&I, &R1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 24, 26, 24, 26, 0, "Z(Ib,Aj)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_TMP0, prsq, 20, 30, "Z(IA,jb)", +1.0);
  global_dpd_->buf4_close(&Z2);

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1,PSIO_OPEN_NEW);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 27, 23, 27, 23, 0, "Z(jA,nI)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, L_irr, 27, 27, 27, 27, 0, "ViAjB");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&I, &R1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 25, 27, 25, 27, 0, "Z(bI,jA)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, sqrp, 21, 30, "Z(AI,jb)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 21, 30, 21, 30, 0, "Z(AI,jb)");
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_TMP0, qprs, 20, 30, "Z(IA,jb)", +1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 20, 10, 20, 10, 0, "Z(IA,nj)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, L_irr, 20, 30, 20, 30, 0, "VIAjb");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract424(&I, &R1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 20, 31, 20, 31, 0, "Z(IA,bj)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_TMP0, pqsr, 20, 30, "Z(IA,jb)", -1.0);
  global_dpd_->buf4_close(&Z2);

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1,PSIO_OPEN_NEW);

  /* - P(ij) P(ab) (Rjmbe Lnmfe) Tif Tna, term 21 */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 20, 0, 20, 0, 0, "Z(JB,NI)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 20, 20, 20, 20, 0, "R2L2_OVOV");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 21, 20, 21, 20, 0, "Z(AI,JB)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_TMP0, qprs, 20, 20, "Z(IA,JB)", -1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 30, 10, 30, 10, 0, "Z(jb,ni)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 30, 30, 30, 30, 0, "R2L2_ovov");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 31, 30, 31, 30, 0, "Z(ai,jb)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_TMP0, qprs, 30, 30, "Z(ia,jb)", -1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 30, 0, 30, 0, 0, "Z(jb,NI)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 30, 20, 30, 20, 0, "R2L2_ovOV");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 21, 30, 21, 30, 0, "Z(AI,jb)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_TMP0, qprs, 20, 30, "Z(IA,jb)", -1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 24, 22, 24, 22, 0, "Z(Ib,Nj)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 24, 24, 24, 24, 0, "R2L2_OvOv");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 24, 26, 24, 26, 0, "Z(Ib,Aj)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_TMP0, prsq, 20, 30, "Z(IA,jb)", +1.0);
  global_dpd_->buf4_close(&Z2);

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1,PSIO_OPEN_NEW);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 27, 23, 27, 23, 0, "Z(jA,nI)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 27, 27, 27, 27, 0, "R2L2_oVoV");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 25, 27, 25, 27, 0, "Z(bI,jA)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 2, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, sqrp, 21, 30, "Z(AI,jb)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 21, 30, 21, 30, 0, "Z(AI,jb)");
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_TMP0, qprs, 20, 30, "Z(IA,jb)", +1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 20, 10, 20, 10, 0, "Z(IA,nj)");
  global_dpd_->buf4_init(&I, PSIF_EOM_TMP, G_irr, 20, 30, 20, 30, 0, "R2L2_OVov");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract424(&I, &T1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 20, 31, 20, 31, 0, "Z(IA,bj)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract244(&T1, &Z, &Z2, 0, 2, 1, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_TMP0, pqsr, 20, 30, "Z(IA,jb)", -1.0);
  global_dpd_->buf4_close(&Z2);

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1,PSIO_OPEN_NEW);


  /* Now permute Z(IA,JB) and Z(ia,jb) and add them in along with Z(IA,jb) */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP0, G_irr, 20, 20, 20, 20, 0, "Z(IA,JB)");
  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, prqs, 0, 5, "Z(IJ,AB)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 0, 5, 2, 7, 0, "GIJAB");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(IJ,AB)");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 0, 5, "Z(JI,AB)");
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, pqsr, 0, 5, "Z(IJ,BA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(JI,AB)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(IJ,BA)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 0, 5, "Z(JI,BA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z(JI,BA)");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP0, G_irr, 30, 30, 30, 30, 0, "Z(ia,jb)");
  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, prqs, 10, 15, "Z(ij,ab)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 10, 15, 12, 17, 0, "Gijab");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 15, 10, 15, 0, "Z(ij,ab)");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 10, 15, "Z(ji,ab)");
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, pqsr, 10, 15, "Z(ij,ba)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 15, 10, 15, 0, "Z(ji,ab)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 15, 10, 15, 0, "Z(ij,ba)");
  global_dpd_->buf4_axpy(&Z, &G, -1.0);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 10, 15, "Z(ji,ba)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 15, 10, 15, 0, "Z(ji,ba)");
  global_dpd_->buf4_axpy(&Z, &G, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP0, G_irr, 20, 30, 20, 30, 0, "Z(IA,jb)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_TMP0, prqs, 22, 28, "GIjAb", 1.0);
  global_dpd_->buf4_close(&Z);

  /* add to ground state parts */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 7, 2, 7, 0, "GIJAB");
  global_dpd_->buf4_init(&V, PSIF_CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "GIJAB");
  global_dpd_->buf4_axpy(&G, &V, 0.5);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 12, 17, 12, 17, 0, "Gijab");
  global_dpd_->buf4_init(&V, PSIF_CC_GAMMA, G_irr, 12, 17, 12, 17, 0, "Gijab");
  global_dpd_->buf4_axpy(&G, &V, 0.5);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&G);
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 28, 22, 28, 0, "GIjAb");
  global_dpd_->buf4_init(&V, PSIF_CC_GAMMA, G_irr, 22, 28, 22, 28, 0, "GIjAb");
  global_dpd_->buf4_axpy(&G, &V, 0.5);
  global_dpd_->buf4_close(&V);
  global_dpd_->buf4_close(&G);

  /*
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 2, 7, 2, 7, 0, "GIJAB");
  tval = dpd_buf4_dot_self(&V);
  dpd_buf4_close(&V);
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 12, 17, 12, 17, 0, "Gijab");
  tval += dpd_buf4_dot_self(&V);
  dpd_buf4_close(&V);
  dpd_buf4_init(&V, CC_GAMMA, G_irr, 22, 28, 22, 28, 0, "GIjAb");
  tval += dpd_buf4_dot_self(&V);
  dpd_buf4_close(&V);
  outfile->Printf("<Gijab|Gijab> = %15.10lf\n", tval);
  */


  return;
}

/* This function computes the following EOM Gijab terms
   Z(i,a) += Rimae Lme ; term 3
   Z(i,a) -= 0.5 Lmnef Tmnea Rif ; term 22
   Z(i,a) -= 0.5 Lmnef Timef Rna; term 23
   Z(i,a) += lmnef rme tinaf; term 33
   P(ij) P(ab) [ Z(i,a) * T(j,b) ]
 */

void x_Gijab_uhf_2(void) {
  int h, nirreps, row, col;
  int i,j,a,b;
  int I1, I2, I3, I4, J1, J2, J3, J4, A1, A2, A3, A4, B1, B2, B3, B4;
  int I1sym, I2sym, I3sym, I4sym, J1sym, J2sym, J3sym, J4sym;
  int A1sym, A2sym, A3sym, A4sym, B1sym, B2sym, B3sym, B4sym;
  int L_irr, R_irr, G_irr;
  dpdfile2 L1R2A, L1R2B, T1A, T1B, Z1A, Z1B, I1A, I1B, R1A, R1B;
  dpdbuf4 G, I, Z, Z2, T, L;

  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* Z(I,A) += L1R2_OV, term 3 */
  global_dpd_->file2_init(&L1R2A, PSIF_EOM_TMP, G_irr, 0, 1, "L1R2_OV");
  global_dpd_->file2_init(&L1R2B, PSIF_EOM_TMP, G_irr, 2, 3, "L1R2_ov");
  global_dpd_->file2_copy(&L1R2A, PSIF_EOM_TMP1, "ZIA");
  global_dpd_->file2_copy(&L1R2B, PSIF_EOM_TMP1, "Zia");
  global_dpd_->file2_close(&L1R2A);
  global_dpd_->file2_close(&L1R2B);

  global_dpd_->file2_init(&Z1A, PSIF_EOM_TMP1, G_irr, 0, 1, "ZIA");
  global_dpd_->file2_init(&Z1B, PSIF_EOM_TMP1, G_irr, 2, 3, "Zia");

  /* Z(I,A) += 0.5 (lmnef tmnea) rif, term 22 */
  global_dpd_->file2_init(&I1A, PSIF_EOM_TMP, L_irr, 1, 1, "LT2_VV");
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract222(&R1A, &I1A, &Z1A, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&R1A);
  global_dpd_->file2_close(&I1A);
  global_dpd_->file2_init(&I1B, PSIF_EOM_TMP, L_irr, 3, 3, "LT2_vv");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract222(&R1B, &I1B, &Z1B, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&R1B);
  global_dpd_->file2_close(&I1B);

  /* Z(i,a) -= 0.5 (timef lnmef) rna; term 23 */
  global_dpd_->file2_init(&I1A, PSIF_EOM_TMP, L_irr, 0, 0, "LT2_OO");
  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract222(&I1A, &R1A, &Z1A, 1, 1, -1.0, 1.0);
  global_dpd_->file2_close(&R1A);
  global_dpd_->file2_close(&I1A);
  global_dpd_->file2_init(&I1B, PSIF_EOM_TMP, L_irr, 2, 2, "LT2_oo");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->contract222(&I1B, &R1B, &Z1B, 1, 1, -1.0, 1.0);
  global_dpd_->file2_close(&R1B);
  global_dpd_->file2_close(&I1B);

  /* Z(i,a) += lmnef rme tinaf; term 33  */
  if (!params.connect_xi) {
    global_dpd_->file2_init(&I1A, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->file2_init(&I1B, PSIF_EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    global_dpd_->dot24(&I1A, &T, &Z1A, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->dot24(&I1B, &T, &Z1A, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    global_dpd_->dot24(&I1B, &T, &Z1B, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->dot24(&I1A, &T, &Z1B, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T);
    global_dpd_->file2_close(&I1A);
    global_dpd_->file2_close(&I1B);
  }

  /* open one-electron files for the nasty permutations */
  global_dpd_->file2_mat_init(&Z1A);   global_dpd_->file2_mat_init(&Z1B);
  global_dpd_->file2_mat_rd(&Z1A);     global_dpd_->file2_mat_rd(&Z1B);

  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->file2_mat_init(&T1A);   global_dpd_->file2_mat_init(&T1B);
  global_dpd_->file2_mat_rd(&T1A);     global_dpd_->file2_mat_rd(&T1B);

  /* + Z(I,A) T(J,B) */
  /* - Z(I,B) T(J,A) */
  /* + T(I,A) Z(J,B) */
  /* - T(I,B) Z(J,A) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 7, 2, 7, 0, "GIJAB");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I1 = Z1A.params->rowidx[i]; I1sym = Z1A.params->psym[i];
      I2 = I1; I2sym = I1sym;
      I3 = T1A.params->rowidx[i]; I3sym = T1A.params->psym[i];
      I4 = I3; I4sym = I3sym;
      J1 = T1A.params->rowidx[j]; J1sym = T1A.params->psym[j];
      J2 = J1; J2sym=J1sym;
      J3 = Z1A.params->rowidx[j]; J3sym = Z1A.params->psym[j];
      J4 = J3; J4sym=J3sym;
      for(col=0; col < G.params->coltot[h]; col++) {
        a = G.params->colorb[h][col][0];
        b = G.params->colorb[h][col][1];
        A1 = Z1A.params->colidx[a]; A1sym = Z1A.params->qsym[a];
        A4 = A1; A4sym = A1sym;
        A2 = T1A.params->colidx[a]; A2sym = T1A.params->qsym[a];
        A3 = A2; A3sym = A2sym;
        B1 = T1A.params->colidx[b]; B1sym = T1A.params->qsym[b];
        B4 = B1; B4sym = B1sym;
        B2 = Z1A.params->colidx[b]; B2sym = Z1A.params->qsym[b];
        B3 = B2; B3sym = B2sym;
        /* + Z(I,A) T(J,B) */
        if ( ((I1sym^A1sym)==G_irr) && (J1sym==B1sym) )
          G.matrix[h][row][col] +=
            Z1A.matrix[I1sym][I1][A1] * T1A.matrix[J1sym][J1][B1];
        /* - Z(I,B) T(J,A) */
        if ( ((I2sym^B2sym)==G_irr) && (J2sym==A2sym) )
          G.matrix[h][row][col] -=
            Z1A.matrix[I2sym][I2][B2] * T1A.matrix[J2sym][J2][A2];
        /* + T(I,A) Z(J,B) */
        if ( ((J3sym^B3sym)==G_irr) && (I3sym==A3sym) )
          G.matrix[h][row][col] +=
            Z1A.matrix[J3sym][J3][B3] * T1A.matrix[I3sym][I3][A3];
        /* - T(I,B) Z(J,A) */
        if ( ((J4sym^A4sym)==G_irr) && (I4sym==B4sym) )
          G.matrix[h][row][col] -=
            Z1A.matrix[J4sym][J4][A4] * T1A.matrix[I4sym][I4][B4];
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 12, 17, 12, 17, 0, "Gijab");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I1 = Z1B.params->rowidx[i]; I1sym = Z1B.params->psym[i];
      I2 = I1; I2sym = I1sym;
      I3 = T1B.params->rowidx[i]; I3sym = T1B.params->psym[i];
      I4 = I3; I4sym = I3sym;
      J1 = T1B.params->rowidx[j]; J1sym = T1B.params->psym[j];
      J2 = J1; J2sym=J1sym;
      J3 = Z1B.params->rowidx[j]; J3sym = Z1B.params->psym[j];
      J4 = J3; J4sym=J3sym;
      for(col=0; col < G.params->coltot[h]; col++) {
        a = G.params->colorb[h][col][0];
        b = G.params->colorb[h][col][1];
        A1 = Z1B.params->colidx[a]; A1sym = Z1B.params->qsym[a];
        A4 = A1; A4sym = A1sym;
        A2 = T1B.params->colidx[a]; A2sym = T1B.params->qsym[a];
        A3 = A2; A3sym = A2sym;
        B1 = T1B.params->colidx[b]; B1sym = T1B.params->qsym[b];
        B4 = B1; B4sym = B1sym;
        B2 = Z1B.params->colidx[b]; B2sym = Z1B.params->qsym[b];
        B3 = B2; B3sym = B2sym;
        /* + Z(i,a) T(j,b) */
        if ( ((I1sym^A1sym)==G_irr) && (J1sym==B1sym) )
          G.matrix[h][row][col] +=
            Z1B.matrix[I1sym][I1][A1] * T1B.matrix[J1sym][J1][B1];
        /* - Z(i,b) T(j,a) */
        if ( ((I2sym^B2sym)==G_irr) && (J2sym==A2sym) )
          G.matrix[h][row][col] -=
            Z1B.matrix[I2sym][I2][B2] * T1B.matrix[J2sym][J2][A2];
        /* + T(i,a) Z(j,b) */
        if ( ((J3sym^B3sym)==G_irr) && (I3sym==A3sym) )
          G.matrix[h][row][col] +=
            Z1B.matrix[J3sym][J3][B3] * T1B.matrix[I3sym][I3][A3];
        /* - T(i,b) Z(j,a) */
        if ( ((J4sym^A4sym)==G_irr) && (I4sym==B4sym) )
          G.matrix[h][row][col] -=
            Z1B.matrix[J4sym][J4][A4] * T1B.matrix[I4sym][I4][B4];
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  /* + Z(I,A) T(j,b) */
  /* + T(I,A) Z(j,b) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 28, 22, 28, 0, "GIjAb");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I1 = Z1A.params->rowidx[i]; I1sym = Z1A.params->psym[i];
      I2 = T1A.params->rowidx[i]; I2sym = T1A.params->psym[i];
      J1 = T1B.params->rowidx[j]; J1sym = T1B.params->psym[j];
      J2 = Z1B.params->rowidx[j]; J2sym = Z1B.params->psym[j];
      for(col=0; col < G.params->coltot[h]; col++) {
        a = G.params->colorb[h][col][0];
        b = G.params->colorb[h][col][1];
        A1 = Z1A.params->colidx[a]; A1sym = Z1A.params->qsym[a];
        A2 = T1A.params->colidx[a]; A2sym = T1A.params->qsym[a];
        B1 = T1B.params->colidx[b]; B1sym = T1B.params->qsym[b];
        B2 = Z1B.params->colidx[b]; B2sym = Z1B.params->qsym[b];
        /* + Z(I,A) T(j,b) */
        if ( ((I1sym^A1sym)==G_irr) && (J1sym==B1sym) )
          G.matrix[h][row][col] +=
            Z1A.matrix[I1sym][I1][A1] * T1B.matrix[J1sym][J1][B1];
        /* + T(I,A) Z(j,b) */
        if ( ((J2sym^B2sym)==G_irr) && (I2sym==A2sym) )
          G.matrix[h][row][col] +=
            T1A.matrix[I2sym][I2][A2] * Z1B.matrix[J2sym][J2][B2];
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&T1A); global_dpd_->file2_mat_close(&T1B);
  global_dpd_->file2_close(&T1A); global_dpd_->file2_close(&T1B);

  global_dpd_->file2_mat_close(&Z1A); global_dpd_->file2_mat_close(&Z1B);
  global_dpd_->file2_close(&Z1A); global_dpd_->file2_close(&Z1B);

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1,PSIO_OPEN_NEW);
  return;
}



/* This function computes the following EOM Gijab terms
   P(ij) P(ab) [ Z(i,a) * R(j,b) ]
   Z(i,a) += Timae Lme, term 18
   Z(i,a) -= 0.5 (Lmnef Tmnea) Tif, term 32
   Z(i,a) -= 0.5 (Lmnef Tmief) Tna; term 34
   Z(i,a) += Lme Tma Tie ; term 19
 */

void x_Gijab_uhf_3(void) {
  int h, nirreps, row, col;
  int i,j,a,b;
  int I1, I2, I3, I4, J1, J2, J3, J4, A1, A2, A3, A4, B1, B2, B3, B4;
  int I1sym, I2sym, I3sym, I4sym, J1sym, J2sym, J3sym, J4sym;
  int A1sym, A2sym, A3sym, A4sym, B1sym, B2sym, B3sym, B4sym;
  int L_irr, R_irr, G_irr;
  dpdfile2 L1T2A, L1T2B, T1A, T1B, Z1A, Z1B, I1A, I1B, R1A, R1B;
  dpdbuf4 G, I, Z, Z2, T, L;

  L_irr = params.L_irr; R_irr = params.R_irr; G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* Z(I,A) += L1T2_OV, term 18 */
  global_dpd_->file2_init(&L1T2A, PSIF_EOM_TMP, L_irr, 0, 1, "L1T2_OV");
  global_dpd_->file2_init(&L1T2B, PSIF_EOM_TMP, L_irr, 2, 3, "L1T2_ov");
  global_dpd_->file2_copy(&L1T2A, PSIF_EOM_TMP1, "ZIA");
  global_dpd_->file2_copy(&L1T2B, PSIF_EOM_TMP1, "Zia");
  global_dpd_->file2_close(&L1T2A);
  global_dpd_->file2_close(&L1T2B);

  global_dpd_->file2_init(&Z1A, PSIF_EOM_TMP1, L_irr, 0, 1, "ZIA");
  global_dpd_->file2_init(&Z1B, PSIF_EOM_TMP1, L_irr, 2, 3, "Zia");

  /* Z(I,A) -= 0.5 Lmnef Tmnea Tif, term 32 */
  global_dpd_->file2_init(&I1A, PSIF_EOM_TMP, L_irr, 1, 1, "LT2_VV");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract222(&T1A, &I1A, &Z1A, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->file2_close(&I1A);
  global_dpd_->file2_init(&I1B, PSIF_EOM_TMP, L_irr, 3, 3, "LT2_vv");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract222(&T1B, &I1B, &Z1B, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->file2_close(&I1B);

  /* Z(i,a) -= 0.5 (Lmnef Tmief ) Tna; term 34 */
  global_dpd_->file2_init(&I1A, PSIF_EOM_TMP, L_irr, 0, 0, "LT2_OO");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract222(&I1A, &T1A, &Z1A, 1, 1, -1.0, 1.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->file2_close(&I1A);
  global_dpd_->file2_init(&I1B, PSIF_EOM_TMP, L_irr, 2, 2, "LT2_oo");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract222(&I1B, &T1B, &Z1B, 1, 1, -1.0, 1.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->file2_close(&I1B);

  /* Z(i,a) += Lme Tma Tie ; term 19  */
  global_dpd_->file2_init(&I1A, PSIF_EOM_TMP, L_irr, 1, 1, "LT1_VV");
  global_dpd_->file2_init(&T1A, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract222(&T1A, &I1A, &Z1A, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&T1A);
  global_dpd_->file2_close(&I1A);
  global_dpd_->file2_init(&I1B, PSIF_EOM_TMP, L_irr, 3, 3, "LT1_vv");
  global_dpd_->file2_init(&T1B, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract222(&T1B, &I1B, &Z1B, 0, 1, -1.0, 1.0);
  global_dpd_->file2_close(&T1B);
  global_dpd_->file2_close(&I1B);

  /* open one-electron files for the nasty terms */
  global_dpd_->file2_mat_init(&Z1A);   global_dpd_->file2_mat_init(&Z1B);
  global_dpd_->file2_mat_rd(&Z1A);     global_dpd_->file2_mat_rd(&Z1B);

  global_dpd_->file2_init(&R1A, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->file2_init(&R1B, PSIF_CC_GR, R_irr, 2, 3, "Ria");
  global_dpd_->file2_mat_init(&R1A);   global_dpd_->file2_mat_init(&R1B);
  global_dpd_->file2_mat_rd(&R1A);     global_dpd_->file2_mat_rd(&R1B);

  /* + Z(I,A) R(J,B) */
  /* - Z(I,B) R(J,A) */
  /* + R(I,A) Z(J,B) */
  /* - R(I,B) Z(J,A) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 2, 7, 2, 7, 0, "GIJAB");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I1 = Z1A.params->rowidx[i]; I1sym = Z1A.params->psym[i];
      I2 = I1; I2sym = I1sym;
      I3 = R1A.params->rowidx[i]; I3sym = R1A.params->psym[i];
      I4 = I3; I4sym = I3sym;
      J1 = R1A.params->rowidx[j]; J1sym = R1A.params->psym[j];
      J2 = J1; J2sym=J1sym;
      J3 = Z1A.params->rowidx[j]; J3sym = Z1A.params->psym[j];
      J4 = J3; J4sym=J3sym;
      for(col=0; col < G.params->coltot[h]; col++) {
        a = G.params->colorb[h][col][0];
        b = G.params->colorb[h][col][1];
        A1 = Z1A.params->colidx[a]; A1sym = Z1A.params->qsym[a];
        A4 = A1; A4sym = A1sym;
        A2 = R1A.params->colidx[a]; A2sym = R1A.params->qsym[a];
        A3 = A2; A3sym = A2sym;
        B1 = R1A.params->colidx[b]; B1sym = R1A.params->qsym[b];
        B4 = B1; B4sym = B1sym;
        B2 = Z1A.params->colidx[b]; B2sym = Z1A.params->qsym[b];
        B3 = B2; B3sym = B2sym;
        /* + Z(I,A) R(J,B) */
        if ( ((I1sym^A1sym)==L_irr) && ((J1sym^B1sym)==R_irr) )
          G.matrix[h][row][col] +=
            Z1A.matrix[I1sym][I1][A1] * R1A.matrix[J1sym][J1][B1];
        /* - Z(I,B) R(J,A) */
        if ( ((I2sym^B2sym)==L_irr) && ((J2sym^A2sym)==R_irr) )
          G.matrix[h][row][col] -=
            Z1A.matrix[I2sym][I2][B2] * R1A.matrix[J2sym][J2][A2];
        /* + R(I,A) Z(J,B) */
        if ( ((J3sym^B3sym)==L_irr) && ((I3sym^A3sym)==R_irr) )
          G.matrix[h][row][col] +=
            Z1A.matrix[J3sym][J3][B3] * R1A.matrix[I3sym][I3][A3];
        /* - R(I,B) Z(J,A) */
        if ( ((J4sym^A4sym)==L_irr) && ((I4sym^B4sym)==R_irr) )
          G.matrix[h][row][col] -=
            Z1A.matrix[J4sym][J4][A4] * R1A.matrix[I4sym][I4][B4];
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 12, 17, 12, 17, 0, "Gijab");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I1 = Z1B.params->rowidx[i]; I1sym = Z1B.params->psym[i];
      I2 = I1; I2sym = I1sym;
      I3 = R1B.params->rowidx[i]; I3sym = R1B.params->psym[i];
      I4 = I3; I4sym = I3sym;
      J1 = R1B.params->rowidx[j]; J1sym = R1B.params->psym[j];
      J2 = J1; J2sym=J1sym;
      J3 = Z1B.params->rowidx[j]; J3sym = Z1B.params->psym[j];
      J4 = J3; J4sym=J3sym;
      for(col=0; col < G.params->coltot[h]; col++) {
        a = G.params->colorb[h][col][0];
        b = G.params->colorb[h][col][1];
        A1 = Z1B.params->colidx[a]; A1sym = Z1B.params->qsym[a];
        A4 = A1; A4sym = A1sym;
        A2 = R1B.params->colidx[a]; A2sym = R1B.params->qsym[a];
        A3 = A2; A3sym = A2sym;
        B1 = R1B.params->colidx[b]; B1sym = R1B.params->qsym[b];
        B4 = B1; B4sym = B1sym;
        B2 = Z1B.params->colidx[b]; B2sym = Z1B.params->qsym[b];
        B3 = B2; B3sym = B2sym;
        /* + Z(i,a) R(j,b) */
        if ( ((I1sym^A1sym)==L_irr) && ((J1sym^B1sym)==R_irr ) )
          G.matrix[h][row][col] +=
            Z1B.matrix[I1sym][I1][A1] * R1B.matrix[J1sym][J1][B1];
        /* - Z(i,b) R(j,a) */
        if ( ((I2sym^B2sym)==L_irr) && ((J2sym^A2sym)==R_irr ) )
          G.matrix[h][row][col] -=
            Z1B.matrix[I2sym][I2][B2] * R1B.matrix[J2sym][J2][A2];
        /* + R(i,a) Z(j,b) */
        if ( ((J3sym^B3sym)==L_irr) && ((I3sym^A3sym)==R_irr ) )
          G.matrix[h][row][col] +=
            Z1B.matrix[J3sym][J3][B3] * R1B.matrix[I3sym][I3][A3];
        /* - R(i,b) Z(j,a) */
        if ( ((J4sym^A4sym)==L_irr) && ((I4sym^B4sym)==R_irr ) )
          G.matrix[h][row][col] -=
            Z1B.matrix[J4sym][J4][A4] * R1B.matrix[I4sym][I4][B4];
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  /* + Z(I,A) R(j,b) */
  /* + R(I,A) Z(j,b) */
  global_dpd_->buf4_init(&G, PSIF_EOM_TMP0, G_irr, 22, 28, 22, 28, 0, "GIjAb");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&G, h);
    global_dpd_->buf4_mat_irrep_rd(&G, h);
    for(row=0; row < G.params->rowtot[h]; row++) {
      i = G.params->roworb[h][row][0];
      j = G.params->roworb[h][row][1];
      I1 = Z1A.params->rowidx[i]; I1sym = Z1A.params->psym[i];
      I2 = R1A.params->rowidx[i]; I2sym = R1A.params->psym[i];
      J1 = R1B.params->rowidx[j]; J1sym = R1B.params->psym[j];
      J2 = Z1B.params->rowidx[j]; J2sym = Z1B.params->psym[j];
      for(col=0; col < G.params->coltot[h]; col++) {
        a = G.params->colorb[h][col][0];
        b = G.params->colorb[h][col][1];
        A1 = Z1A.params->colidx[a]; A1sym = Z1A.params->qsym[a];
        A2 = R1A.params->colidx[a]; A2sym = R1A.params->qsym[a];
        B1 = R1B.params->colidx[b]; B1sym = R1B.params->qsym[b];
        B2 = Z1B.params->colidx[b]; B2sym = Z1B.params->qsym[b];
        /* + Z(I,A) R(j,b) */
        if ( ((I1sym^A1sym)==L_irr) && ((J1sym^B1sym)==R_irr ) )
          G.matrix[h][row][col] +=
            Z1A.matrix[I1sym][I1][A1] * R1B.matrix[J1sym][J1][B1];
        /* + R(I,A) Z(j,b) */
        if ( ((J2sym^B2sym)==L_irr) && ((I2sym^A2sym)==R_irr ) )
          G.matrix[h][row][col] +=
            R1A.matrix[I2sym][I2][A2] * Z1B.matrix[J2sym][J2][B2];
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&G, h);
    global_dpd_->buf4_mat_irrep_close(&G, h);
  }
  global_dpd_->buf4_close(&G);

  global_dpd_->file2_mat_close(&R1A); global_dpd_->file2_mat_close(&R1B);
  global_dpd_->file2_close(&R1A); global_dpd_->file2_close(&R1B);

  global_dpd_->file2_mat_close(&Z1A); global_dpd_->file2_mat_close(&Z1B);
  global_dpd_->file2_close(&Z1A); global_dpd_->file2_close(&Z1B);

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1,PSIO_OPEN_NEW);
  return;
}

}} // namespace psi::ccdensity
