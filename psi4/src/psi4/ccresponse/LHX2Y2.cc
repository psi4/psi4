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
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstring>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

double LHX2Y2(const char *pert_x, int irrep_x, double omega_x,
	      const char *pert_y, int irrep_y, double omega_y)
{
  dpdbuf4 X2, Y2, I, W, Z, Z1, Z2, W1, W2, L2;
  dpdfile2 z;
  char lbl[32];
  double polar;

  sprintf(lbl, "Z_%s_MbeJ (Me,Jb)", pert_y);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep_y, 10, 10, 10, 10, 0, lbl);
  sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_y, omega_y);
  global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_y, 10, 10, 10, 10, 0, lbl);
  global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
  global_dpd_->contract444(&I, &Y2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&Y2);
  global_dpd_->buf4_close(&Z);

  sprintf(lbl, "Z_%s_MbEj (ME,jb)", pert_y);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep_y, 10, 10, 10, 10, 0, lbl);

  sprintf(lbl, "X_%s_(2IAjb-jAIb) (%5.3f)", pert_y, omega_y);
  global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_y, 10, 10, 10, 10, 0, lbl);
  global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
  global_dpd_->contract444(&I, &Y2, &Z, 0, 1, 0.5, 0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&Y2);

  sprintf(lbl, "Z_%s_MbeJ (Me,Jb)", pert_y);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, irrep_y, 10, 10, 10, 10, 0, lbl);
  global_dpd_->buf4_axpy(&Z1, &Z, -0.5);
  sprintf(lbl, "Z_%s_(2MbEj+MbeJ) (ME,JB)", pert_y);
  global_dpd_->buf4_scmcopy(&Z, PSIF_CC_TMP0, lbl, 2);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep_y, 10, 10, 10, 10, 0, lbl);
  global_dpd_->buf4_axpy(&Z1, &Z, 1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&Z1);


  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb) I");
  sprintf(lbl, "X_%s_(2IAjb-IbjA) (%5.3f)", pert_x, omega_x);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_x, 10, 10, 10, 10, 0, lbl);
  sprintf(lbl, "Z_%s_(2MbEj+MbeJ) (ME,JB)", pert_y);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep_y, 10, 10, 10, 10, 0, lbl);
  global_dpd_->contract444(&X2, &Z, &Z1, 0, 1, 0.5, 0);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&X2);
  global_dpd_->buf4_close(&Z1);

  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb) Ib");
  sprintf(lbl, "Z_%s_MbeJ (Me,Jb)", pert_y);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep_y, 10, 10, 10, 10, 0, lbl);
  sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_x, omega_x);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_x, 10, 10, 10, 10, 0, lbl);
  global_dpd_->contract444(&X2, &Z, &Z1, 0, 1, 1, 0);
  global_dpd_->buf4_close(&X2);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&Z1);

  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb) I");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb) Ib");
  global_dpd_->buf4_axpy(&Z2, &Z1, 0.5);
  global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, psrq, 10, 10, "Z(IA,jb) III");
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&Z1);

  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb) I");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(IA,jb) III");
  global_dpd_->buf4_axpy(&Z2, &Z1, 1);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, prqs, 0, 5, "Z(Ij,Ab) I+III");
  global_dpd_->buf4_close(&Z1);

  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) I+III");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qpsr, 0, 5, "Z(Ij,Ab) II+IV");
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  global_dpd_->buf4_scm(&Z, 0);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) I+III");
  global_dpd_->buf4_axpy(&Z1, &Z, 1);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) II+IV");
  global_dpd_->buf4_axpy(&Z1, &Z, 1);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&Z);


  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  sprintf(lbl, "Z_%s_MnIj", pert_x);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, irrep_x, 0, 0, 0, 0, 0, lbl);
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_x, omega_x);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
  global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->contract444(&I, &X2, &Z1, 0, 0, 1, 0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&X2);
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
  global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  global_dpd_->contract444(&Z1, &Y2, &Z, 1, 1, 1, 1);
  global_dpd_->buf4_close(&Y2);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  sprintf(lbl, "Z_%s_MnIj", pert_y);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, irrep_y, 0, 0, 0, 0, 0, lbl);
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  global_dpd_->contract444(&I, &X2, &Z1, 0, 0, 1, 0);
  global_dpd_->buf4_close(&I);
  global_dpd_->buf4_close(&X2);
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_x, omega_x);
  global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
  global_dpd_->contract444(&Z1, &Y2, &Z, 1, 1, 1, 1);
  global_dpd_->buf4_close(&Y2);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&Z);


  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");

  sprintf(lbl, "Z_%s_AE", pert_y);
  global_dpd_->file2_init(&z, PSIF_CC_TMP0, irrep_y, 1, 1, lbl);
  global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
  global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  global_dpd_->contract442(&Y2, &I, &z, 3, 3, -1, 0);
  global_dpd_->buf4_close(&Y2);
  global_dpd_->buf4_close(&I);

  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_x, omega_x);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
  global_dpd_->contract424(&X2, &z, &Z1, 3, 1, 0, 1, 0);
  global_dpd_->buf4_close(&X2);
  global_dpd_->file2_close(&z);
  global_dpd_->buf4_axpy(&Z1, &Z, 1);
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
  global_dpd_->buf4_axpy(&Z1, &Z, 1);

  global_dpd_->buf4_close(&Z);


  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");

  sprintf(lbl, "Z_%s_AE", pert_x);
  global_dpd_->file2_init(&z, PSIF_CC_TMP0, irrep_x, 1, 1, lbl);
  global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_x, omega_x);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
  global_dpd_->contract442(&X2, &I, &z, 3, 3, -1, 0);
  global_dpd_->buf4_close(&X2);
  global_dpd_->buf4_close(&I);

  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
  global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  global_dpd_->contract424(&Y2, &z, &Z1, 3, 1, 0, 1, 0);
  global_dpd_->buf4_close(&Y2);
  global_dpd_->file2_close(&z);
  global_dpd_->buf4_axpy(&Z1, &Z, 1);
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
  global_dpd_->buf4_axpy(&Z1, &Z, 1);

  global_dpd_->buf4_close(&Z);


  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");

  sprintf(lbl, "Z_%s_MI", pert_y);
  global_dpd_->file2_init(&z, PSIF_CC_TMP0, irrep_y, 0, 0, lbl);
  global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
  global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  global_dpd_->contract442(&I, &Y2, &z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&Y2);
  global_dpd_->buf4_close(&I);

  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_x, omega_x);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
  global_dpd_->contract244(&z, &X2, &Z1, 0, 0, 0, 1, 0);
  global_dpd_->buf4_close(&X2);
  global_dpd_->file2_close(&z);
  global_dpd_->buf4_axpy(&Z1, &Z, -1);
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
  global_dpd_->buf4_axpy(&Z1, &Z, -1);

  global_dpd_->buf4_close(&Z);


  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");

  sprintf(lbl, "Z_%s_MI", pert_x);
  global_dpd_->file2_init(&z, PSIF_CC_TMP0, irrep_x, 0, 0, lbl);
  global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_x, omega_x);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
  global_dpd_->contract442(&I, &X2, &z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&X2);
  global_dpd_->buf4_close(&I);

  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
  global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
  global_dpd_->contract244(&z, &Y2, &Z1, 0, 0, 0, 1, 0);
  global_dpd_->buf4_close(&Y2);
  global_dpd_->file2_close(&z);
  global_dpd_->buf4_axpy(&Z1, &Z, -1);
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
  global_dpd_->buf4_axpy(&Z1, &Z, -1);

  global_dpd_->buf4_close(&Z);


  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
  polar = global_dpd_->buf4_dot(&L2, &Z);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_close(&Z);

  return polar;
}

}} // namespace psi::ccresponse
