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

void x_xi2_4_rohf(void);
void x_xi2_14(void);
void x_xi2_rohf(void);
extern void x_xi_check(char *term_lbl);
extern void x_xi2_uhf(void);
extern void x_xi2_rhf(void);

/* compute xi_2 amplitudes for zeta equations */

void x_xi2(void) {
  if (params.ref == 0)
    x_xi2_rhf();
  else if (params.ref == 1)
    x_xi2_rohf();
  else
    x_xi2_uhf();
  return;
}

void x_xi2_rohf(void)
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
  global_dpd_->buf4_init(&D2, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
  global_dpd_->buf4_scmcopy(&D2, PSIF_EOM_XI, "XIJAB", params.overlap1+params.overlap2);
  global_dpd_->buf4_scmcopy(&D2, PSIF_EOM_XI, "Xijab", params.overlap1+params.overlap2);
  global_dpd_->buf4_close(&D2);
  global_dpd_->buf4_init(&D2, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->buf4_scmcopy(&D2, PSIF_EOM_XI, "XIjAb", params.overlap1+params.overlap2);
  global_dpd_->buf4_close(&D2);
#ifdef DEBUG_XI
x_xi_check("terms 1 and 5");
#endif

  /* terms 2 and 9, Xijab -= P(ab) (Lma Rme + Lmnfa Rmnfe) <ij||eb */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, 0, 2, 5, 2, 5, 0, "Z (I>J,AB)");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 1, 1, "LR_VV");
  global_dpd_->buf4_init(&D2, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
  global_dpd_->contract244(&I1, &D2, &Z2, 1, 2, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&D2);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 2, 5, 2, 7, 0, "XIJAB");
  global_dpd_->buf4_axpy(&Z2, &XIJAB, -1.0);
  global_dpd_->buf4_close(&XIJAB);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, pqsr, 2, 7, "XIJAB", 1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, 0, 2, 5, 2, 5, 0, "Z (i>j,ab)");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 1, 1, "LR_vv");
  global_dpd_->buf4_init(&D2, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
  global_dpd_->contract244(&I1, &D2, &Z2, 1, 2, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&D2);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_init(&Xijab, PSIF_EOM_XI, G_irr, 2, 5, 2, 7, 0, "Xijab");
  global_dpd_->buf4_axpy(&Z2, &Xijab, -1.0);
  global_dpd_->buf4_close(&Xijab);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, pqsr, 2, 7, "Xijab", 1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 1, 1, "LR_VV");
  global_dpd_->buf4_init(&D2, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->contract244(&I1, &D2, &XIjAb, 1, 2, 1, -1.0, 1.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 1, 1, "LR_vv");
  global_dpd_->contract424(&D2, &I1, &XIjAb, 3, 1, 0, -1.0, 1.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_close(&D2);
  global_dpd_->buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("terms 2 and 9");
#endif

  /* terms 3 and 10, Xijab -= P(ij) (Lie Rme + 0.5 Linef Rmnef) <mj||ab> */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, 0, 0, 7, 0, 7, 0, "Z (IJ,A>B)");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 0, "LR_OO");
  global_dpd_->buf4_init(&D2, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
  global_dpd_->contract244(&I1, &D2, &Z2, 1, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&D2);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 0, 7, 2, 7, 0, "XIJAB");
  global_dpd_->buf4_axpy(&Z2, &XIJAB, -1.0);
  global_dpd_->buf4_close(&XIJAB);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qprs, 2, 7, "XIJAB", 1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, 0, 0, 7, 0, 7, 0, "Z (ij,a>b)");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 0, "LR_oo");
  global_dpd_->buf4_init(&D2, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
  global_dpd_->contract244(&I1, &D2, &Z2, 1, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&D2);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_init(&Xijab, PSIF_EOM_XI, G_irr, 0, 7, 2, 7, 0, "Xijab");
  global_dpd_->buf4_axpy(&Z2, &Xijab, -1.0);
  global_dpd_->buf4_close(&Xijab);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qprs, 2, 7, "Xijab", 1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 0, "LR_OO");
  global_dpd_->buf4_init(&D2, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->contract244(&I1, &D2, &XIjAb, 1, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 0, "LR_oo");
  global_dpd_->contract424(&D2, &I1, &XIjAb, 1, 1, 1, -1.0, 1.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_close(&D2);
  global_dpd_->buf4_close(&XIjAb);

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);
#ifdef DEBUG_XI
x_xi_check("terms 3 and 10");
#endif

  x_xi2_4_rohf();

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);
#ifdef DEBUG_XI
x_xi_check("terms 4 and 6");
#endif

  /* term 7, Xijab += 0.25 Lmnab Rmnef <ij||ef> */
  global_dpd_->buf4_init(&D2, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
  global_dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 2, 7, 2, 7, 0, "RIJAB");
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, R_irr, 2, 2, 2, 2, 0, "Z (I>J,M>N)");
  global_dpd_->contract444(&D2, &R2, &Z2, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&R2);
  global_dpd_->buf4_close(&D2);
  global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 2, 7, 2, 7, 0, "LIJAB");
  global_dpd_->contract444(&Z2, &L2, &XIJAB, 0, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&XIJAB);

  global_dpd_->buf4_init(&D2, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
  global_dpd_->buf4_init(&R2, PSIF_CC_GR, R_irr, 2, 7, 2, 7, 0, "Rijab");
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, R_irr, 2, 2, 2, 2, 0, "Z (i>j,m>n)");
  global_dpd_->contract444(&D2, &R2, &Z2, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&R2);
  global_dpd_->buf4_close(&D2);
  global_dpd_->buf4_init(&Xijab, PSIF_EOM_XI, G_irr, 2, 7, 2, 7, 0, "Xijab");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 2, 7, 2, 7, 0, "Lijab");
  global_dpd_->contract444(&Z2, &L2, &Xijab, 0, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&Xijab);

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

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);
#ifdef DEBUG_XI
x_xi_check("term 7");
#endif

  /* term 8, Xijab += 0.25 Rmnef Lijef <mn||ab> */
  global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 2, 2, 2, 2, 0, "R2L2_OOOO");
  global_dpd_->buf4_init(&D2, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
  global_dpd_->contract444(&I2, &D2, &XIJAB, 1, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&D2);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&XIJAB);
  global_dpd_->buf4_init(&Xijab, PSIF_EOM_XI, G_irr, 2, 7, 2, 7, 0, "Xijab");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 2, 2, 2, 2, 0, "R2L2_oooo");
  global_dpd_->buf4_init(&D2, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
  global_dpd_->contract444(&I2, &D2, &Xijab, 1, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&D2);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&Xijab);
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
  global_dpd_->file2_copy(&I1, PSIF_EOM_TMP1, "Z (F,A)");
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 1, 1, "Z (F,A)");
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

      /* for term 11: */
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP_XI, R_irr, 1, 1, "RD_vv");
  global_dpd_->file2_copy(&I1, PSIF_EOM_TMP1, "Z (f,a)");
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 1, 1, "Z (f,a)");
      /* for term 20: */
  global_dpd_->file2_init(&Z1B, PSIF_EOM_TMP_XI, R_irr, 1, 1, "R1Wamef_vv");
  global_dpd_->file2_axpy(&Z1B, &I1, -1.0, 0);
  global_dpd_->file2_close(&Z1B);
      /* for term 17: */
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "Fme");
  global_dpd_->contract222(&R1, &F1, &I1, 1, 1, 1.0, 1.0);
  global_dpd_->file2_close(&F1);
  global_dpd_->file2_close(&R1);
  global_dpd_->file2_close(&I1);

  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z (I>J,AB)");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 1, 1, "Z (F,A)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 2, 5, 2, 7, 0, "LIJAB");
  global_dpd_->contract244(&I1, &L2, &Z2, 0, 2, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 2, 5, 2, 7, 0, "XIJAB");
  global_dpd_->buf4_axpy(&Z2, &XIJAB, -1.0);
  global_dpd_->buf4_close(&XIJAB);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, pqsr, 2, 7, "XIJAB", 1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z (i>j,ab)");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 1, 1, "Z (f,a)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 2, 5, 2, 7, 0, "Lijab");
  global_dpd_->contract244(&I1, &L2, &Z2, 0, 2, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_init(&Xijab, PSIF_EOM_XI, G_irr, 2, 5, 2, 7, 0, "Xijab");
  global_dpd_->buf4_axpy(&Z2, &Xijab, -1.0);
  global_dpd_->buf4_close(&Xijab);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, pqsr, 2, 7, "Xijab", 1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 1, 1, "Z (F,A)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  global_dpd_->contract244(&I1, &L2, &XIjAb, 0, 2, 1, -1.0, 1.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 1, 1, "Z (f,a)");
  global_dpd_->contract424(&L2, &I1, &XIjAb, 3, 0, 0, -1.0, 1.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("terms 11, 17 and 20");
#endif

  /* term 12, Xijab -= 0.5 P(ij) Lmjab (Rmnef <in||ef>) */
  /* term 16,       -=     P(ij) Lmjab (Rme Fie) */
  /* term 21,       -=     P(ij) Lmjab (Rne Winme) */
  /* make 1-electron intermediates to include terms 16 and 21 as well */
      /* for term 12: */
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP_XI, R_irr, 0, 0, "RD_OO");
  global_dpd_->file2_copy(&I1, PSIF_EOM_TMP1, "Z (M,I)");
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 0, 0, "Z (M,I)");
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

      /* for term 12 */
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP_XI, R_irr, 0, 0, "RD_oo");
  global_dpd_->file2_copy(&I1, PSIF_EOM_TMP1, "Z (m,i)");
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 0, 0, "Z (m,i)");
      /* for term 21 */
  global_dpd_->file2_init(&Z1B, PSIF_EOM_TMP_XI, R_irr, 0, 0, "R1Wmnie_oo");
  global_dpd_->file2_axpy(&Z1B, &I1, 1.0, 1);
  global_dpd_->file2_close(&Z1B);
      /* for term 16 */
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "Fme");
  global_dpd_->contract222(&R1, &F1, &I1, 0, 0, 1.0, 1.0);
  global_dpd_->file2_close(&F1);
  global_dpd_->file2_close(&R1);
  global_dpd_->file2_close(&I1);

  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z (IJ,A>B)");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 0, 0, "Z (M,I)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 7, 2, 7, 0, "LIJAB");
  global_dpd_->contract244(&I1, &L2, &Z2, 0, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 0, 7, 2, 7, 0, "XIJAB");
  global_dpd_->buf4_axpy(&Z2, &XIJAB, -1.0);
  global_dpd_->buf4_close(&XIJAB);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qprs, 2, 7, "XIJAB", 1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z (ij,a>b)");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 0, 0, "Z (m,i)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 7, 2, 7, 0, "Lijab");
  global_dpd_->contract244(&I1, &L2, &Z2, 0, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 0, 7, 2, 7, 0, "Xijab");
  global_dpd_->buf4_axpy(&Z2, &XIJAB, -1.0);
  global_dpd_->buf4_close(&XIJAB);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qprs, 2, 7, "Xijab", 1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 0, 0, "Z (M,I)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  global_dpd_->contract244(&I1, &L2, &XIjAb, 0, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&I1);
  global_dpd_->file2_init(&I1, PSIF_EOM_TMP1, R_irr, 0, 0, "Z (m,i)");
  global_dpd_->contract424(&L2, &I1, &XIjAb, 1, 0, 1, -1.0, 1.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->file2_close(&I1);
  global_dpd_->buf4_close(&XIjAb);
#ifdef DEBUG_XI
x_xi_check("term 12, 16 and 21");
#endif

  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);

    /* term 13 + 15, (0.25 Rmnef <mn||ef> + Rme Fme) Lijab */
  if ( (L_irr == 0) && (!params.connect_xi)) {
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "FME");
    tval = global_dpd_->file2_dot(&R1, &F1);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_close(&R1);
    global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
    global_dpd_->file2_init(&F1, PSIF_CC_OEI, 0, 0, 1, "Fme");
    tval += global_dpd_->file2_dot(&R1, &F1);
    global_dpd_->file2_close(&F1);
    global_dpd_->file2_close(&R1);
    tval += params.RD_overlap;

    global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    global_dpd_->buf4_axpy(&L2, &XIJAB, tval);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&XIJAB);
    global_dpd_->buf4_init(&Xijab, PSIF_EOM_XI, G_irr, 2, 7, 2, 7, 0, "Xijab");
    global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 2, 7, 2, 7, 0, "Lijab");
    global_dpd_->buf4_axpy(&L2, &Xijab, tval);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&Xijab);
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
  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);

    /* term 22, +P(ij) (Limfe Rme) Wfjab */
  if (!params.connect_xi) {
    global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z (IJ,A>B)");
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "WAMEF");
    global_dpd_->contract244(&I1, &H2, &Z2, 1, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&I1);
    global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 0, 7, 2, 7, 0, "XIJAB");
    global_dpd_->buf4_axpy(&Z2, &XIJAB, 1.0);
    global_dpd_->buf4_close(&XIJAB);
    global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qprs, 2, 7, "XIJAB", -1.0);
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 7, 0, 7, 0, "Z (ij,a>b)");
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "Wamef");
    global_dpd_->contract244(&I1, &H2, &Z2, 1, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&I1);
    global_dpd_->buf4_init(&Xijab, PSIF_EOM_XI, G_irr, 0, 7, 2, 7, 0, "Xijab");
    global_dpd_->buf4_axpy(&Z2, &Xijab, 1.0);
    global_dpd_->buf4_close(&Xijab);
    global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qprs, 2, 7, "Xijab", -1.0);
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->contract244(&I1, &H2, &XIjAb, 1, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&I1);
    global_dpd_->buf4_close(&XIjAb);

    global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z (jI,bA)");
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
    global_dpd_->contract244(&I1, &H2, &Z2, 1, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&I1);
    global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qpsr, 0, 5, "XIjAb", 1.0);
    global_dpd_->buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("term 22 (Wamef)");
#endif
    psio_close(PSIF_EOM_TMP1,0);
    psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);
  }

    /* term 23, -P(ab) (Lmnea Rme) Wijnb */
  if (!params.connect_xi) {
    global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z (I>J,AB)");
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 2, 10, 2, 10, 0, "WMNIE");
    global_dpd_->contract244(&I1, &H2, &Z2, 0, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&I1);
    global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 2, 5, 2, 7, 0, "XIJAB");
    global_dpd_->buf4_axpy(&Z2, &XIJAB, -1.0);
    global_dpd_->buf4_close(&XIJAB);
    global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, pqsr, 2, 7, "XIJAB", 1.0);
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 2, 5, 2, 5, 0, "Z (i>j,ab)");
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 2, 10, 2, 10, 0, "Wmnie");
    global_dpd_->contract244(&I1, &H2, &Z2, 0, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&I1);
    global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 2, 5, 2, 7, 0, "Xijab");
    global_dpd_->buf4_axpy(&Z2, &XIJAB, -1.0);
    global_dpd_->buf4_close(&XIJAB);
    global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, pqsr, 2, 7, "Xijab", 1.0);
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
    global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->contract244(&I1, &H2, &XIjAb, 0, 2, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&I1);
    global_dpd_->buf4_close(&XIjAb);

    global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z (jI,Ab)");
    global_dpd_->file2_init(&I1, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WmNiE (mN,Ei)");
    global_dpd_->contract424(&H2, &I1, &Z2, 3, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&H2);
    global_dpd_->file2_close(&I1);
    global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qprs, 0, 5, "XIjAb", -1.0);
    global_dpd_->buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("term 23 (Wmnie)");
#endif
  }

  /* term 25, Xijab += (Lnmab Rme) Wijne */
  global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 7, 10, 7, 10, 0, "L2R1_VVOV");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 2, 10, 2, 10, 0, "WMNIE");
  global_dpd_->contract444(&H2, &I2, &XIJAB, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&XIJAB);

  global_dpd_->buf4_init(&Xijab, PSIF_EOM_XI, G_irr, 2, 7, 2, 7, 0, "Xijab");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 7, 10, 7, 10, 0, "L2R1_vvov");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 2, 10, 2, 10, 0, "Wmnie");
  global_dpd_->contract444(&H2, &I2, &Xijab, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&Xijab);

  global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 5, 10, 5, 10, 0, "L2R1_VvOv");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
  global_dpd_->contract444(&H2, &I2, &XIjAb, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&XIjAb);

  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z (jI,Ab)");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 5, 10, 5, 10, 0, "L2R1_VvoV");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WmNiE");
  global_dpd_->contract444(&H2, &I2, &Z2, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qprs, 0, 5, "XIjAb", 1.0);
  global_dpd_->buf4_close(&Z2);

#ifdef DEBUG_XI
x_xi_check("term 25 (Wmnie)");
#endif
  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);

    /* term 24, Xijab -= (Lijfe Rme) Wfmab */
  global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L2R1_OOVO");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "WAMEF");
  global_dpd_->contract444(&I2, &H2, &XIJAB, 0, 1, -1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&XIJAB);

  global_dpd_->buf4_init(&Xijab, PSIF_EOM_XI, G_irr, 2, 7, 2, 7, 0, "Xijab");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 2, 11, 2, 11, 0, "L2R1_oovo");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 7, 11, 7, 0, "Wamef");
  global_dpd_->contract444(&I2, &H2, &Xijab, 0, 1, -1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&Xijab);

  global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OoVo");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  global_dpd_->contract444(&I2, &H2, &XIjAb, 0, 1, -1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_close(&XIjAb);

  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z (Ij,bA)");
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 0, 11, 0, 11, 0, "L2R1_OovO");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
  global_dpd_->contract444(&I2, &H2, &Z2, 0, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, pqsr, 0, 5, "XIjAb", -1.0);
  global_dpd_->buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("term 24 (Wamef)");
#endif
  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);

  /* terms 18, 19: Xijab -= P(ij) P(ab) Linae (Rme Wmjnb + Rnf Wejbf) */
  /* construct Z(JB,NE) = RME WMJNB + RNF WEJBF */
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 11, 11, 11, 11, 0, "Z (EJ,BN)");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 11, 2, 11, 0, "WMNIE (M>N,EI)");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract244(&R1, &H2, &Z, 0, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 7, 0, "WAMEF");
  global_dpd_->contract424(&H2, &R1, &Z, 3, 1, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qpsr, 10, 10, "Z (JE,NB)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (JE,NB)");
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 10, 10, "Z (JB,NE)");
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 11, 10, 11, 10, 0, "Z (eJ,nB)");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WmNiE");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->contract244(&R1, &H2, &Z, 0, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
  global_dpd_->contract244(&R1, &H2, &Z, 1, 2, 1, -1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qprs, 10, 10, "Z (Je,nB)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (Je,nB)");
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 10, 10, "Z (JB,ne)");
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 11, 11, 11, 11, 0, "Z (ej,bn)");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 11, 2, 11, 0, "Wmnie (m>n,ei)");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->contract244(&R1, &H2, &Z, 0, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 7, 0, "Wamef");
  global_dpd_->contract424(&H2, &R1, &Z, 3, 1, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&H2);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, qpsr, 10, 10, "Z (je,nb)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (je,nb)");
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, psrq, 10, 10, "Z (jb,ne)");
  global_dpd_->buf4_close(&Z);

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

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 11, 10, 11, 0, "Z (Je,bN)");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->contract424(&H2, &R1, &Z, 1, 0, 1, -1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, prsq, 10, 10, "Z (Jb,Ne)");
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 11, 11, 11, 11, 0, "Z2 (eJ,bN)");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WaMeF");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&H2, &R1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_TMP1, qrsp, 10, 10, "Z (Jb,Ne)", 1.0);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 11, 10, 11, 0, "Z (jE,Bn)");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WmNiE (mN,Ei)");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->contract424(&H2, &R1, &Z, 1, 0, 1, -1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP1, prsq, 10, 10, "Z (jB,nE)");
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 11, 11, 11, 11, 0, "Z2 (Ej,Bn)");
  global_dpd_->buf4_init(&H2, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
  global_dpd_->file2_init(&R1, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->contract424(&H2, &R1, &Z, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&R1);
  global_dpd_->buf4_close(&H2);
  global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_TMP1, qrsp, 10, 10, "Z (jB,nE)", 1.0);
  global_dpd_->buf4_close(&Z);

  /* XIJAB -= P(IJ) P(AB) L(IA,NE) Z(NE,JB) */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (IA,JB)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAJB");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (JB,NE)");
  global_dpd_->contract444(&L2, &Z, &Z2, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (JB,ne)");
  global_dpd_->contract444(&L2, &Z, &Z2, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, prqs, 0, 5, "Z2 (IJ,AB)");
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z2 (IJ,AB)");
  global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 0, 5, 2, 7, 0, "XIJAB");
  global_dpd_->buf4_axpy(&Z2, &XIJAB, -1.0);
  global_dpd_->buf4_close(&XIJAB);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qprs, 2, 7, "XIJAB", 1.0);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, pqsr, 2, 7, "XIJAB", 1.0);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qpsr, 2, 7, "XIJAB", -1.0);
  global_dpd_->buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("terms 18, 19 (Wmnie, Wamef)");
#endif

  /* Xijab -= P(ij) P(ab) L(ia,ne) Z(ne,jb) */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (ia,jb)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "Liajb");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (jb,ne)");
  global_dpd_->contract444(&L2, &Z, &Z2, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LiaJB");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (jb,NE)");
  global_dpd_->contract444(&L2, &Z, &Z2, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, prqs, 0, 5, "Z2 (ij,ab)");
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
  global_dpd_->buf4_init(&Xijab, PSIF_EOM_XI, G_irr, 0, 5, 2, 7, 0, "Xijab");
  global_dpd_->buf4_axpy(&Z2, &Xijab, -1.0);
  global_dpd_->buf4_close(&Xijab);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qprs, 2, 7, "Xijab", 1.0);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, pqsr, 2, 7, "Xijab", 1.0);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qpsr, 2, 7, "Xijab", -1.0);
  global_dpd_->buf4_close(&Z2);

  /* XIjAb += - L(IA,NE) Z(jb,NE) - L(IA,ne) Z(jb,ne)  */
  /*          - L(jb,ne) Z(IA,ne) - L(jb,NE) Z(IA,NE)  */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (IA,jb)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAJB");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (jb,NE)");
  global_dpd_->contract444(&L2, &Z, &Z2, 0, 0, 1.0, 0.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIAjb");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (jb,ne)");
  global_dpd_->contract444(&L2, &Z, &Z2, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (JB,ne)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "Liajb");
  global_dpd_->contract444(&Z, &L2, &Z2, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (JB,NE)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LiaJB");
  global_dpd_->contract444(&Z, &L2, &Z2, 0, 0, 1.0, 1.0);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, prqs, 0, 5, "XIjAb", -1.0);
  global_dpd_->buf4_close(&Z2);

  /* XIjAb += + L(jA,Ne) Z(Ib,Ne) + L(Ib,nE) Z(jA,nE)  */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (Ib,jA)");
  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LjAIb");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (Jb,Ne)");
  global_dpd_->contract444(&Z, &L2, &Z2, 0, 0, -1.0, 0.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&L2);

  global_dpd_->buf4_init(&L2, PSIF_CC_GL, L_irr, 10, 10, 10, 10, 0, "LIbjA");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, R_irr, 10, 10, 10, 10, 0, "Z (jB,nE)");
  global_dpd_->contract444(&L2, &Z, &Z2, 0, 0, -1.0, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, prsq, 0, 5, "XIjAb", 1.0);
  global_dpd_->buf4_close(&Z2);
#ifdef DEBUG_XI
x_xi_check("terms 18, 19 (Wmnie, Wamef)");
#endif
  psio_close(PSIF_EOM_TMP1,0);
  psio_open(PSIF_EOM_TMP1, PSIO_OPEN_NEW);

  /* Write irrep of XI amplitudes to CC_INFO */
  psio_write_entry(PSIF_CC_INFO, "XI Irrep", (char *) &G_irr,sizeof(int));

  global_dpd_->file2_init(&XIA, PSIF_EOM_XI, G_irr, 0, 1, "XIA");
  tval = global_dpd_->file2_dot_self(&XIA);
  global_dpd_->file2_close(&XIA);
  outfile->Printf("XIA amplitudes: norm=%20.15lf dot=%20.15lf\n", sqrt(tval), tval );
  global_dpd_->file2_init(&Xia, PSIF_EOM_XI, G_irr, 0, 1, "Xia");
  tval += global_dpd_->file2_dot_self(&Xia);
  outfile->Printf("X1 amplitudes:  norm=%20.15lf dot=%20.15lf\n", sqrt(tval), tval );
  global_dpd_->file2_close(&Xia);
  global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 2, 7, 2, 7, 0, "XIJAB");
  tval += global_dpd_->buf4_dot_self(&XIJAB);
  global_dpd_->buf4_close(&XIJAB);
  global_dpd_->buf4_init(&Xijab, PSIF_EOM_XI, G_irr, 2, 7, 2, 7, 0, "Xijab");
  tval += global_dpd_->buf4_dot_self(&Xijab);
  global_dpd_->buf4_close(&Xijab);
  global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
  tval += global_dpd_->buf4_dot_self(&XIjAb);
  global_dpd_->buf4_close(&XIjAb);
  outfile->Printf("Norm of Xi: %20.15lf\n", sqrt(tval) );
  return;
}



/* compute terms 4 and 6 of xi2 amplitudes */
/* Xijab += P(ij) P(ab) (Rme Lia + Rmnef Linaf) <mj||eb> */
void x_xi2_4_rohf(void)
{
  dpdfile2 RIA, Ria, LIA, Lia;
  int L_irr, R_irr, G_irr, nirreps;
  int I, A, M, E, i, a, m, e, h, row, col, Isym, Esym, Asym, Msym;
  dpdbuf4 D, R2, L2, H2, I2, Z, Z2, XIJAB, Xijab, XIjAb;

  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  /* construct RL = Rme Lia + Rmnef Linaf */
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVOV");
  global_dpd_->buf4_copy(&I2, PSIF_EOM_TMP1, "RL_OVOV");
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovov");
  global_dpd_->buf4_copy(&I2, PSIF_EOM_TMP1, "RL_ovov");
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OVov");
  global_dpd_->buf4_copy(&I2, PSIF_EOM_TMP1, "RL_OVov");
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_ovOV");
  global_dpd_->buf4_copy(&I2, PSIF_EOM_TMP1, "RL_ovOV");
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_oVoV");
  global_dpd_->buf4_copy(&I2, PSIF_EOM_TMP1, "RL_oVoV");
  global_dpd_->buf4_close(&I2);
  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP, G_irr, 10, 10, 10, 10, 0, "R2L2_OvOv");
  global_dpd_->buf4_copy(&I2, PSIF_EOM_TMP1, "RL_OvOv");
  global_dpd_->buf4_close(&I2);

  /* RL_OVOV(me,ia) += Rme Lia */
  global_dpd_->file2_init(&RIA, PSIF_CC_GR, R_irr, 0, 1, "RIA");
  global_dpd_->file2_init(&LIA, PSIF_CC_GL, L_irr, 0, 1, "LIA");
  global_dpd_->file2_init(&Ria, PSIF_CC_GR, R_irr, 0, 1, "Ria");
  global_dpd_->file2_init(&Lia, PSIF_CC_GL, L_irr, 0, 1, "Lia");

  global_dpd_->file2_mat_init(&RIA);
  global_dpd_->file2_mat_init(&Ria);
  global_dpd_->file2_mat_init(&LIA);
  global_dpd_->file2_mat_init(&Lia);
  global_dpd_->file2_mat_rd(&RIA);
  global_dpd_->file2_mat_rd(&Ria);
  global_dpd_->file2_mat_rd(&LIA);
  global_dpd_->file2_mat_rd(&Lia);

  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_OVOV");
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

  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_ovov");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&I2, h);
    global_dpd_->buf4_mat_irrep_rd(&I2, h);
    for(row=0; row < I2.params->rowtot[h]; row++) {
      m = I2.params->roworb[h][row][0];
      e = I2.params->roworb[h][row][1];
      M = Ria.params->rowidx[m]; Msym = Ria.params->psym[m];
      E = Ria.params->colidx[e]; Esym = Ria.params->qsym[e];
      for(col=0; col < I2.params->coltot[h^G_irr]; col++) {
        i = I2.params->colorb[h^G_irr][col][0];
        a = I2.params->colorb[h^G_irr][col][1];
        I = Lia.params->rowidx[i]; Isym = Lia.params->psym[i];
        A = Lia.params->colidx[a]; Asym = Lia.params->qsym[a];
        if( ((Msym^Esym)==R_irr) && ((Isym^Asym)==L_irr) )
          I2.matrix[h][row][col] += Ria.matrix[Msym][M][E] * Lia.matrix[Isym][I][A];
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&I2, h);
    global_dpd_->buf4_mat_irrep_close(&I2, h);
  }
  global_dpd_->buf4_close(&I2);

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
        I = Lia.params->rowidx[i]; Isym = Lia.params->psym[i];
        A = Lia.params->colidx[a]; Asym = Lia.params->qsym[a];
        if( ((Msym^Esym)==R_irr) && ((Isym^Asym)==L_irr) )
          I2.matrix[h][row][col] += RIA.matrix[Msym][M][E] * Lia.matrix[Isym][I][A];
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&I2, h);
    global_dpd_->buf4_mat_irrep_close(&I2, h);
  }
  global_dpd_->buf4_close(&I2);

  global_dpd_->buf4_init(&I2, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_ovOV");
  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&I2, h);
    global_dpd_->buf4_mat_irrep_rd(&I2, h);
    for(row=0; row < I2.params->rowtot[h]; row++) {
      m = I2.params->roworb[h][row][0];
      e = I2.params->roworb[h][row][1];
      M = Ria.params->rowidx[m]; Msym = Ria.params->psym[m];
      E = Ria.params->colidx[e]; Esym = Ria.params->qsym[e];
      for(col=0; col < I2.params->coltot[h^G_irr]; col++) {
        i = I2.params->colorb[h^G_irr][col][0];
        a = I2.params->colorb[h^G_irr][col][1];
        I = LIA.params->rowidx[i]; Isym = LIA.params->psym[i];
        A = LIA.params->colidx[a]; Asym = LIA.params->qsym[a];
        if( ((Msym^Esym)==R_irr) && ((Isym^Asym)==L_irr) )
          I2.matrix[h][row][col] += Ria.matrix[Msym][M][E] * LIA.matrix[Isym][I][A];
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&I2, h);
    global_dpd_->buf4_mat_irrep_close(&I2, h);
  }
  global_dpd_->buf4_close(&I2);

  global_dpd_->file2_mat_close(&RIA);
  global_dpd_->file2_mat_close(&Ria);
  global_dpd_->file2_mat_close(&LIA);
  global_dpd_->file2_mat_close(&Lia);
  global_dpd_->file2_close(&RIA);
  global_dpd_->file2_close(&Ria);
  global_dpd_->file2_close(&LIA);
  global_dpd_->file2_close(&Lia);

  /* Z2(ME,IA) <MJ|EB>(ME,JB) + Z2(me,IA) <mJ|eB>(me,JB) */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (IA,JB)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_OVOV");
  global_dpd_->contract444(&Z, &D, &Z2, 1, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_ovOV");
  global_dpd_->contract444(&Z, &D, &Z2, 1, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, prqs, 0, 5, "Z2 (IJ,AB)");
  global_dpd_->buf4_close(&Z2);

  /* add in Z2(IJ,AB) permutations */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z2 (IJ,AB)");
  global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 0, 5, 2, 7, 0, "XIJAB");
  global_dpd_->buf4_axpy(&Z2, &XIJAB, 1.0);
  global_dpd_->buf4_close(&XIJAB);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qprs, 2, 7, "XIJAB", -1.0);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, pqsr, 2, 7, "XIJAB", -1.0);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qpsr, 2, 7, "XIJAB", 1.0);
  global_dpd_->buf4_close(&Z2);

  /* Z2(me,ia) <mj|eb>(me,jb) + Z2(ME,ia) <Mj|Eb>(ME,jb) */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (ia,jb)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_ovov");
  global_dpd_->contract444(&Z, &D, &Z2, 1, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_OVov");
  global_dpd_->contract444(&Z, &D, &Z2, 1, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_sort(&Z2, PSIF_EOM_TMP1, prqs, 0, 5, "Z2 (ij,ab)");
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
  global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 0, 5, 2, 7, 0, "Xijab");
  global_dpd_->buf4_axpy(&Z2, &XIJAB, 1.0);
  global_dpd_->buf4_close(&XIJAB);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qprs, 2, 7, "Xijab", -1.0);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, pqsr, 2, 7, "Xijab", -1.0);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, qpsr, 2, 7, "Xijab", 1.0);
  global_dpd_->buf4_close(&Z2);

  /* Z2(ME,IA) <Mj|Eb>(ME,jb) + Z2(me,IA) <mj|eb>(me,jb) */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (IA,jb)");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_OVOV");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  global_dpd_->contract444(&Z, &D, &Z2, 1, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_ovOV");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
  global_dpd_->contract444(&Z, &D, &Z2, 1, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&D);

  /* XIjAb += ZMjEb <MI|EA> + Zmjeb <mI|eA> */
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_OVov");
  global_dpd_->contract444(&D, &Z, &Z2, 1, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_ovov");
  global_dpd_->contract444(&D, &Z, &Z2, 1, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, prqs, 0, 5, "XIjAb", 1.0);
  global_dpd_->buf4_close(&Z2);

  /* XIjAb += ZmjEA <mI|bE> + ZMIeb <Mj|Ae> */
  global_dpd_->buf4_init(&Z2, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "Z2 (Ib,jA)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_oVoV");
  global_dpd_->contract444(&D, &Z, &Z2, 1, 1, 1.0, 0.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&D);

  global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 10, 10, 10, 0, "RL_OvOv");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  global_dpd_->contract444(&Z, &D, &Z2, 1, 1, 1.0, 1.0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_sort_axpy(&Z2, PSIF_EOM_XI, prsq, 0, 5, "XIjAb", 1.0);
  global_dpd_->buf4_close(&Z2);
  return;
}



/* compute term 14 of Xi2 amplitudes for all spin cases */
/* Xijab += P(ij) P(ab) (Lmjeb Rme) Fia */

void x_xi2_14(void)
{
  dpdfile2 IIA, Iia, FME, Fme;
  int L_irr, R_irr, G_irr, nirreps;
  int I, A, J, B, i, a, j, b, h, row, col, Isym, Asym, Jsym, Bsym;
  int II, AA, JJ, BB, ii, aa, jj, bb, IIsym, AAsym, JJsym, BBsym;
  dpdbuf4 Z, XIJAB, Xijab, XIjAb;

  L_irr = params.L_irr;
  R_irr = params.R_irr;
  G_irr = params.G_irr;
  nirreps = moinfo.nirreps;

  global_dpd_->file2_init(&IIA, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_OV");
  global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->file2_mat_init(&IIA);
  global_dpd_->file2_mat_init(&FME);
  global_dpd_->file2_mat_rd(&IIA);
  global_dpd_->file2_mat_rd(&FME);

  if (params.ref == 1) {
    global_dpd_->file2_init(&Iia, PSIF_EOM_TMP, G_irr, 0, 1, "L2R1_ov");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");
  }
  else if (params.ref == 2) {
    global_dpd_->file2_init(&Iia, PSIF_EOM_TMP, G_irr, 2, 3, "L2R1_ov");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
  }

  if (params.ref != 0) {
    global_dpd_->file2_mat_init(&Iia);
    global_dpd_->file2_mat_rd(&Iia);
    global_dpd_->file2_mat_init(&Fme);
    global_dpd_->file2_mat_rd(&Fme);
  }

  /* build Z(IJ,AB) = FME(I,A) L2R1_OV(J,B) */
  if (params.ref != 0) { /* ROHF or UHF for XIJAB */
    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z (IJ,AB)");
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&Z, h);
      for(row=0; row < Z.params->rowtot[h]; row++) {
        i = Z.params->roworb[h][row][0];
        j = Z.params->roworb[h][row][1];
        I = FME.params->rowidx[i]; Isym = FME.params->psym[i];
        J = IIA.params->rowidx[j]; Jsym = IIA.params->psym[j];
        for(col=0; col < Z.params->coltot[h^G_irr]; col++) {
          a = Z.params->colorb[h^G_irr][col][0];
          b = Z.params->colorb[h^G_irr][col][1];
          A = FME.params->colidx[a]; Asym = FME.params->qsym[a];
          B = IIA.params->colidx[b]; Bsym = IIA.params->qsym[b];
          if( (Isym==Asym) && ((Jsym^Bsym)==G_irr) )
            Z.matrix[h][row][col] += FME.matrix[Isym][I][A] * IIA.matrix[Jsym][J][B];
        }
      }
      global_dpd_->buf4_mat_irrep_wrt(&Z, h);
      global_dpd_->buf4_mat_irrep_close(&Z, h);
    }
    /* XIJAB += P(IJ) P(AB) Z(IJ,AB) */
    global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, G_irr, 0, 5, 2, 7, 0, "XIJAB");
    global_dpd_->buf4_axpy(&Z, &XIJAB, 1.0);
    global_dpd_->buf4_close(&XIJAB);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_XI, qprs, 2, 7, "XIJAB", -1.0);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_XI, pqsr, 2, 7, "XIJAB", -1.0);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_XI, qpsr, 2, 7, "XIJAB", 1.0);
    global_dpd_->buf4_close(&Z);
  }

  /* build Z(ij,ab) = Fme(i,a) L2R1_ov(j,b) */
  if (params.ref != 0) { /* ROHF or UHF for Xijab */
    if (params.ref == 1)
      global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 0, 5, 0, 5, 0, "Z (ij,ab)");
    else if (params.ref == 2)
      global_dpd_->buf4_init(&Z, PSIF_EOM_TMP1, G_irr, 10, 15, 10, 15, 0, "Z (ij,ab)");
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&Z, h);
      for(row=0; row < Z.params->rowtot[h]; row++) {
        i = Z.params->roworb[h][row][0];
        j = Z.params->roworb[h][row][1];
        I = Fme.params->rowidx[i]; Isym = Fme.params->psym[i];
        J = Iia.params->rowidx[j]; Jsym = Iia.params->psym[j];
        for(col=0; col < Z.params->coltot[h^G_irr]; col++) {
          a = Z.params->colorb[h^G_irr][col][0];
          b = Z.params->colorb[h^G_irr][col][1];
          A = Fme.params->colidx[a]; Asym = Fme.params->qsym[a];
          B = Iia.params->colidx[b]; Bsym = Iia.params->qsym[b];
          if( (Isym==Asym) && ((Jsym^Bsym)==G_irr) )
            Z.matrix[h][row][col] += Fme.matrix[Isym][I][A] * Iia.matrix[Jsym][J][B];
        }
      }
      global_dpd_->buf4_mat_irrep_wrt(&Z, h);
      global_dpd_->buf4_mat_irrep_close(&Z, h);
    }
    /* Xijab += P(ij) P(ab) Z(ij,ab) */
    if (params.ref == 1) {
      global_dpd_->buf4_init(&Xijab, PSIF_EOM_XI, G_irr, 0, 5, 2, 7, 0, "Xijab");
      global_dpd_->buf4_axpy(&Z, &Xijab, 1.0);
      global_dpd_->buf4_close(&Xijab);
      global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_XI, qprs, 2, 7, "Xijab", -1.0);
      global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_XI, pqsr, 2, 7, "Xijab", -1.0);
      global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_XI, qpsr, 2, 7, "Xijab", 1.0);
      global_dpd_->buf4_close(&Z);
    }
    else if (params.ref == 2) {
      global_dpd_->buf4_init(&Xijab, PSIF_EOM_XI, G_irr, 10, 15, 12, 17, 0, "Xijab");
      global_dpd_->buf4_axpy(&Z, &Xijab, 1.0);
      global_dpd_->buf4_close(&Xijab);
      global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_XI, qprs, 12, 17, "Xijab", -1.0);
      global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_XI, pqsr, 12, 17, "Xijab", -1.0);
      global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_XI, qpsr, 12, 17, "Xijab", 1.0);
      global_dpd_->buf4_close(&Z);
    }
  }

  /* XIjAb += FME(I,A) L2R1_ov(j,b) + IIA(I,A) F(j,b) */
  if (params.ref == 0) { /* RHF XIjAb */
    global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&XIjAb, h);
      global_dpd_->buf4_mat_irrep_rd(&XIjAb, h);
      for(row=0; row < XIjAb.params->rowtot[h]; row++) {
        i = XIjAb.params->roworb[h][row][0];
        j = XIjAb.params->roworb[h][row][1];
        I = FME.params->rowidx[i]; Isym = FME.params->psym[i];
        J = IIA.params->rowidx[j]; Jsym = IIA.params->psym[j];
        II = IIA.params->rowidx[i]; IIsym = IIA.params->psym[i];
        JJ = FME.params->rowidx[j]; JJsym = FME.params->psym[j];
        for(col=0; col < XIjAb.params->coltot[h^G_irr]; col++) {
          a = XIjAb.params->colorb[h^G_irr][col][0];
          b = XIjAb.params->colorb[h^G_irr][col][1];
          A = FME.params->colidx[a]; Asym = FME.params->qsym[a];
            B = IIA.params->colidx[b]; Bsym = IIA.params->qsym[b];
            AA = IIA.params->colidx[a]; AAsym = IIA.params->qsym[a];
            BB = FME.params->colidx[b]; BBsym = FME.params->qsym[b];
            if( (Isym==Asym) && ((Jsym^Bsym)==G_irr) )
              XIjAb.matrix[h][row][col] += FME.matrix[Isym][I][A] * IIA.matrix[Jsym][J][B];
            if( ((IIsym^AAsym)==G_irr) && (JJsym==BBsym) )
              XIjAb.matrix[h][row][col] += IIA.matrix[IIsym][II][AA] * FME.matrix[JJsym][JJ][BB];
        }
      }
      global_dpd_->buf4_mat_irrep_wrt(&XIjAb, h);
      global_dpd_->buf4_mat_irrep_close(&XIjAb, h);
    }
    global_dpd_->buf4_close(&XIjAb);
    global_dpd_->file2_mat_close(&IIA);
    global_dpd_->file2_close(&IIA);
    global_dpd_->file2_mat_close(&FME);
    global_dpd_->file2_close(&FME);
  }
  else { /* ROHF and UHF XIjAb */
    if (params.ref == 1)
      global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 0, 5, 0, 5, 0, "XIjAb");
    else if (params.ref == 2)
      global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, G_irr, 22, 28, 22, 28, 0, "XIjAb");
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&XIjAb, h);
      global_dpd_->buf4_mat_irrep_rd(&XIjAb, h);
      for(row=0; row < XIjAb.params->rowtot[h]; row++) {
        i = XIjAb.params->roworb[h][row][0];
        j = XIjAb.params->roworb[h][row][1];
        I = FME.params->rowidx[i]; Isym = FME.params->psym[i];
        J = Iia.params->rowidx[j]; Jsym = Iia.params->psym[j];
        II = IIA.params->rowidx[i]; IIsym = IIA.params->psym[i];
        JJ = Fme.params->rowidx[j]; JJsym = Fme.params->psym[j];
        for(col=0; col < XIjAb.params->coltot[h^G_irr]; col++) {
          a = XIjAb.params->colorb[h^G_irr][col][0];
          b = XIjAb.params->colorb[h^G_irr][col][1];
          A = FME.params->colidx[a]; Asym = FME.params->qsym[a];
            B = Iia.params->colidx[b]; Bsym = Iia.params->qsym[b];
            AA = IIA.params->colidx[a]; AAsym = IIA.params->qsym[a];
            BB = Fme.params->colidx[b]; BBsym = Fme.params->qsym[b];
            if( (Isym==Asym) && ((Jsym^Bsym)==G_irr) )
              XIjAb.matrix[h][row][col] += FME.matrix[Isym][I][A] * Iia.matrix[Jsym][J][B];
            if( ((IIsym^AAsym)==G_irr) && (JJsym==BBsym) )
              XIjAb.matrix[h][row][col] += IIA.matrix[IIsym][II][AA] * Fme.matrix[JJsym][JJ][BB];
        }
      }
      global_dpd_->buf4_mat_irrep_wrt(&XIjAb, h);
      global_dpd_->buf4_mat_irrep_close(&XIjAb, h);
    }
    global_dpd_->buf4_close(&XIjAb);
    global_dpd_->file2_mat_close(&IIA); global_dpd_->file2_mat_close(&Iia);
    global_dpd_->file2_close(&IIA); global_dpd_->file2_close(&Iia);
    global_dpd_->file2_mat_close(&FME); global_dpd_->file2_mat_close(&Fme);
    global_dpd_->file2_close(&FME); global_dpd_->file2_close(&Fme);
  }

  return;
}

}} // namespace psi::ccdensity
