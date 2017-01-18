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
#include "math.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

void x_xi2_4_rhf(void);
extern void x_xi2_14(void);
extern void x_xi_check(char *term_lbl);
extern double norm_C_rhf(dpdfile2 *CME, dpdbuf4 *CMnEf, dpdbuf4 *CMnfE);

/* compute xi_2 amplitudes for RHF wavefunctions for zeta equations */

void x_xi2_rhf(void)
{
  dpdfile2 L1, XIA, Xia, I1, R1, F1, Z1A, Z1B;
  int L_irr, R_irr, G_irr;
  double tval;
  dpdbuf4 D2, R2, L2, H2, I2, Z, Z2, XIJAB, Xijab, XIjAb;

  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;

#ifdef DEBUG_XI
  x_xi_check("reset");
#endif

  /* terms 1 and 5, Xijab += (Lme Rme + 0.25 Lmnef Rmnef) <ij||eb> */
  /* overlaps in params are assigned in x_xi1.c */
  /* see comments in xi1_connected.c */
  global_dpd_->buf4_init(&D2, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->buf4_scmcopy(&D2, PSIF_EOM_XI, "XIjAb", params.overlap1+params.overlap2);
  global_dpd_->buf4_close(&D2);
#ifdef DEBUG_XI
x_xi_check("terms 1 and 5");
#endif

  /* terms 2 and 9, Xijab -= P(ab) (Lma Rme + Lmnfa Rmnfe) <ij||eb */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "X (Ij,Ab)");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 1, 1, "LR_VV");
  global_dpd_->buf4_init(&D2, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->contract424(&D2, &I1, &Z2, 3, 1, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&D2);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->buf4_axpy(&Z2, &XIjAb, -1.0);
  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, qpsr, 0, 5, "X (jI,bA)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "X (jI,bA)");
  global_dpd_->buf4_axpy(&Z2, &XIjAb, -1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("terms 2 and 9");
#endif

  /* terms 3 and 10, Xijab -= P(ij) (Lie Rme + 0.5 Linef Rmnef) <mj||ab> */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "X (Ij,Ab)");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 0, "LR_OO");
  global_dpd_->buf4_init(&D2, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->contract244(&I1, &D2, &Z2, 1, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&D2);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->buf4_axpy(&Z2, &XIjAb, -1.0);
  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, qpsr, 0, 5, "X (jI,bA)");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "X (jI,bA)");
  global_dpd_->buf4_axpy(&Z2, &XIjAb, -1.0);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("terms 3 and 10");
#endif

  x_xi2_4_rhf();
#ifdef DEBUG_XI
x_xi_check("terms 4 and 6");
#endif

  /* term 7, Xijab += 0.25 Lmnab Rmnef <ij||ef> */
  global_dpd_->buf4_init(&D2, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 0, 5, 0, 5, 0, "RIjAb");
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, R_irr, 0, 0, 0, 0, 0, "Z (Ij,Mn)");
  global_dpd_->contract444(&D2, &R2, &Z2, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&R2);
  global_dpd_->buf4_close(&D2);
  global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  global_dpd_->contract444(&Z2, &L2, &XIjAb, 0, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("term 7");
#endif

  /* term 8, Xijab += 0.25 Rmnef Lijef <mn||ab> */
  global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 0, 0, 0, 0, 0, "R2L2_OoOo");
  global_dpd_->buf4_init(&D2, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->contract444(&I2, &D2, &XIjAb, 1, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&D2);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("term 8");
#endif

  /* term 11, Xijab -= 0.5 P(ab) Lijfb (Rmnef <mn||ea>) */
  /* term 17        -=     P(ab) Lijfb (Rmf Fma) */
  /* term 20        +=     P(ab) Lijfb (Rme Wfmae) */
  /* build 1-e intermediates to include term 17 */
      /* for term 11: */
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP_XI, R_irr, 1, 1, "RD_VV");
  global_dpd_->file2_copy(&I1, PSIF_EOM_TMP1, "X1 (F,A)");
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 1, 1, "X1 (F,A)");
      /* for term 20: */
  global_dpd_->file2_init(&Z1A, PSIF_EOM_TMP_XI, R_irr, 1, 1, "R1Wamef_VV");
  global_dpd_->file2_axpy(&Z1A, &I1, -1.0, 0);
  global_dpd_->file2_close(&Z1A);
      /* for term 17: */
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->contract222(&R1, &F1, &I1, 1, 1, 1.0, 1.0);
  global_dpd_->file2_close(&F1);
  global_dpd_->file2_close(&R1);
  global_dpd_->file2_close(&I1);

  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 1, 1, "X1 (F,A)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  global_dpd_->contract244(&I1, &L2, &Z2, 0, 2, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->buf4_axpy(&Z2, &XIjAb, -1.0);
  global_dpd_->buf4_close(&XIjAb);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qpsr, 0, 5, "XIjAb", -1.0);
  global_dpd_->buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("terms 11, 17 and 20");
#endif

  /* term 12, Xijab -= 0.5 P(ij) Lmjab (Rmnef <in||ef>) */
  /* term 16,       -=     P(ij) Lmjab (Rme Fie) */
  /* term 21,       -=     P(ij) Lmjab (Rne Winme) */
  /* make 1-electron intermediates to include terms 16 and 21 as well */
      /* for term 12: */
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
  global_dpd_->file2_copy(&I1, PSIF_EOM_TMP1, "X1 (M,I)");
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 0, 0, "X1 (M,I)");
       /* for term 21 */
  global_dpd_->file2_init(&Z1A, PSIF_EOM_TMP_XI, R_irr, 0, 0, "R1Wmnie_OO");
  global_dpd_->file2_axpy(&Z1A, &I1, 1.0, 1);
  global_dpd_->file2_close(&Z1A);
      /* for term 16: */
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->contract222(&R1, &F1, &I1, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&F1);
  global_dpd_->file2_close(&R1);
  global_dpd_->file2_close(&I1);

  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 0, 0, "X1 (M,I)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  global_dpd_->contract244(&I1, &L2, &Z2, 0, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->buf4_axpy(&Z2, &XIjAb, -1.0);
  global_dpd_->buf4_close(&XIjAb);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qpsr, 0, 5, "XIjAb", -1.0);
  global_dpd_->buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("term 12, 16 and 21");
#endif

  /* term 13 + 15, (0.25 Rmnef <mn||ef> + Rme Fme) Lijab */
  if ( (L_irr == 0) && (!params.connect_xi)) {
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "FME");
    tval = 2.0 * global_dpd_->file2_dot(&R1, &F1);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_close(&R1);
    tval += params.RD_overlap;

    global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->buf4_axpy(&L2, &XIjAb, tval);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("term 13 (ijab) and 15 (Fme)");
#endif
  }

  /* term 14, +P(ij) P(ab) Lmjeb Rme Fia */
  if (!params.connect_xi) {
    x_xi2_14();
#ifdef DEBUG_XI
x_xi_check("term 14 (Fme)");
#endif
  }

  /* term 22, +P(ij) (Limfe Rme) Wfjab */
  if (!params.connect_xi) {
    global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "X (Ij,Ab)");
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->contract244(&I1, &H2, &Z2, 1, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&I1);
    global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    global_dpd_->buf4_axpy(&Z2, &XIjAb, 1.0);
    global_dpd_->buf4_close(&XIjAb);
    global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qpsr, 0, 5, "XIjAb", 1.0);
    global_dpd_->buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("term 22 (Wamef)");
#endif
  }

  /* term 23, -P(ab) (Lmnea Rme) Wijnb */
  if (!params.connect_xi) {
    global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "X (Ij,Ab)");
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->contract244(&I1, &H2, &Z2, 0, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&I1);
    global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    global_dpd_->buf4_axpy(&Z2, &XIjAb, -1.0);
    global_dpd_->buf4_close(&XIjAb);
    global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qpsr, 0, 5, "XIjAb", -1.0);
    global_dpd_->buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("term 23 (Wmnie)");
#endif
  }

  /* term 25, Xijab += (Lnmab Rme) Wijne */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, R_irr, 0, 0, 0, 0, 0, "X (Ij,Nm)");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&H2, &R1, &Z2, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  global_dpd_->contract444(&Z2, &L2, &Z, 0, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->buf4_axpy(&Z, &XIjAb, 1.0);
  global_dpd_->buf4_close(&XIjAb);
  global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_XI, qpsr, 0, 5, "XIjAb", 1.0);
  global_dpd_->buf4_close(&Z);
#ifdef DEBUG_XI
x_xi_check("term 25 (Wmnie)");
#endif

  /* term 24, Xijab -= (Lijfe Rme) Wfmab */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "X (Ij,Ab)");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OoVo");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  global_dpd_->contract444(&I2, &H2, &Z2, 0, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->buf4_axpy(&Z2, &XIjAb, -1.0);
  global_dpd_->buf4_close(&XIjAb);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qpsr, 0, 5, "XIjAb", -1.0);
  global_dpd_->buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("term 24 (Wamef)");
#endif

  /* terms 18, 19: Xijab -= P(ij) P(ab) Linae (Rme Wmjnb + Rnf Wejbf) */
  /* construct Z(JB,NE) = RME WMJNB + RNF WEJBF */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 11, 10, 11, 10, 0, "Z (Ej,Nb)");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract244(&R1, &H2, &Z, 0, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  global_dpd_->contract244(&R1, &H2, &Z, 1, 2, 1, -1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 10, 10, "Z (jE,Nb)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (jE,Nb)");
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 10, 10, "Z (jb,NE)");
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (Je,Nb)");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&H2, &R1, &Z, 1, 0, 1, -1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 11, 11, 11, 11, 0, "Z (eJ,bN)");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&H2, &R1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_TMP1, qpsr, 10, 10, "Z (Je,Nb)", 1.0);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (Je,Nb)");
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 10, 10, "Z (Jb,Ne)");
  global_dpd_->buf4_close(&Z);

  /* XIjAb += - L(IA,NE) Z(jb,NE) - L(IA,ne) Z(jb,ne)  */
  /*          - L(jb,ne) Z(IA,ne) - L(jb,NE) Z(IA,NE)  transpose of line above */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "XIjAb (IA,jb)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "2LIjAb - LIjbA (IA,jb)");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (jb,NE)");
  global_dpd_->contract444(&L2, &Z, &Z2, 0, 0, -1.0, 0.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (Jb,Ne)");
  global_dpd_->contract444(&L2, &Z, &Z2, 0, 0, -1.0, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&L2);

  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, rspq, 10, 10, "XIjAb (jb,IA)");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "XIjAb (jb,IA)");
  global_dpd_->buf4_axpy(&Z, &Z2, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, prqs, 0, 5, "XIjAb", 1.0);
  global_dpd_->buf4_close(&Z2);

  /* XIjAb += + L(jA,Ne) Z(Ib,Ne) + L(Ib,nE) Z(jA,nE)  */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (Ib,jA)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LjAIb");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (Jb,Ne)");
  global_dpd_->contract444(&Z, &L2, &Z2, 0, 0, -1.0, 0.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&L2);

  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, rspq, 10, 10, "Z2 (jA,Ib)");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (jA,Ib)");
  global_dpd_->buf4_axpy(&Z, &Z2, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, prsq, 0, 5, "XIjAb", 1.0);
  global_dpd_->buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("terms 18, 19 (Wmnie, Wamef)");
#endif

  /* Write irrep of XI amplitudes to CC_INFO */
  psio_write_entry(PSIF_CC_INFO, "XI Irrep", (char *) &G_irr,sizeof(int));

  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  tval = 2.0 * global_dpd_->file2_dot_self(&XIA);
  global_dpd_->file2_close(&XIA);
  outfile->Printf("XI_IA amplitudes: Norm=%15.10lf, Dot=%15.10lf\n", sqrt(tval), tval );

  global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->buf4_sort(&XIjAb, PSIF_EOM_TMP1, pqsr, 0, 5, "XIjbA");
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "XIjbA");
  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  tval = norm_C_rhf(&XIA, &XIjAb, &Z2);
  global_dpd_->file2_close(&XIA);
  global_dpd_->buf4_close(&XIjAb);
  global_dpd_->buf4_close(&Z2);
  outfile->Printf("XI amplitudes   : Norm=%15.10lf, Dot=%15.10lf\n", sqrt(tval), tval );

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);
  return;
}



/* compute terms 4 and 6 of xi2 amplitudes */
/* Xijab += P(ij) P(ab) (Rme Lia + Rmnef Linaf) <mj||eb> */
void x_xi2_4_rhf(void)
{
  dpdfile2 RIA, LIA;
  int L_irr, R_irr, G_irr, nirreps;
  int I, A, M, E, i, a, m, e, h, row, col, Isym, Esym, Asym, Msym;
  dpdbuf4 D, R2, L2, H2, I2, Z, Z2, XIjAb;

  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* construct RL = Rme Lia + Rmnef Linaf */
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
  global_dpd_->buf4_copy(&I2, PSIF_EOM_TMP1, "RL_OVov");
  global_dpd_->buf4_close(&I2);

  /* RL_OVOV(me,ia) += Rme Lia */
  global_dpd_->file2_init(&RIA, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->file2_init(&LIA, PSIF_CC_GL, L_irr, 0, 1, "LIA");

  global_dpd_->file2_mat_init(&RIA);
  global_dpd_->file2_mat_init(&LIA);
  global_dpd_->file2_mat_rd(&RIA);
  global_dpd_->file2_mat_rd(&LIA);

  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_OVov");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&I2, h);
    global_dpd_->buf4_mat_irrep_rd(&I2, h);
    for(row=0; row < I2.params->rowtot[h]; row++) {
      m = I2.params->roworb[h][row][0];
      e = I2.params->roworb[h][row][1];
      M = RIA.params->rowidx[m]; Msym = RIA.params->psym[m];
      E = RIA.params->colidx[e]; Esym = RIA.params->qsym[e];
      for(col=0; col < I2.params->coltot[h^G_irr]; col++) {
        i = I2.params->colorb[h^G_irr][col][0];
        a = I2.params->colorb[h^G_irr][col][1];
        I = LIA.params->rowidx[i]; Isym = LIA.params->psym[i];
        A = LIA.params->colidx[a]; Asym = LIA.params->qsym[a];
        if( ((Msym^Esym)==R_irr) && ((Isym^Asym)==L_irr) )
          I2.matrix[h][row][col] += RIA.matrix[Msym][M][E] * LIA.matrix[Isym][I][A];
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&I2, h);
    global_dpd_->buf4_mat_irrep_close(&I2, h);
  }
  global_dpd_->buf4_close(&I2);

  /* XIjAb += RL_OVov(me,IA) * (2<Mj|Eb> - <Mj|Be>) for the ZOVOV alpha case
   and the reverse of this multiplication for the Zovov beta case */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "X2 (IA,jb)");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_OVov");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
  global_dpd_->contract444(&Z, &D, &Z2, 1, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);
  /* XIjAb += RL_OvOv(Me,Ia) * (<Mj|Eb> (ME,jb) and the reverse of this
   to finish the all beta and all alpha spin orbital components */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  global_dpd_->contract444(&D, &Z, &Z2, 1, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, prqs, 0, 5, "XIjAb");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->buf4_axpy(&Z2, &XIjAb, 1.0);
  global_dpd_->buf4_close(&XIjAb);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qpsr, 0, 5, "XIjAb", 1.0);
  global_dpd_->buf4_close(&Z2);

  /* Now do Z alpha beta parts */
  /* XIjAb += ZmEjA <mI|bE> (mE,Ib) + ZMeIb <Mj|Ae> (Me,jA) */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, 0, 10, 10, 10, 10, 0, "X2 (Ib,jA)");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  global_dpd_->contract444(&Z, &D, &Z2, 1, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&D);

  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, prsq, 0, 5, "XIjAb");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->buf4_axpy(&Z2, &XIjAb, 1.0);
  global_dpd_->buf4_close(&XIjAb);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qpsr, 0, 5, "XIjAb", 1.0);
  global_dpd_->buf4_close(&Z2);

  return;
}


}} // namespace psi::ccdensity
