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
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

double LCX(const char *pert_c, int irrep_c,
	   const char *pert_x, int irrep_x, double omega)
{
  double polar=0.0;
  dpdfile2 X1, mu1, z1, l1, mu, lt, xc;
  dpdbuf4 X2, mu2, z2, l2, Z;
  char lbl[32];

  /*** Mu * X1 ***/

  sprintf(lbl, "%s_IA", pert_c);
  global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep_c, 0, 1, lbl);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega);
  global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
  polar += 2.0 * global_dpd_->file2_dot(&mu1, &X1);
  global_dpd_->file2_close(&X1);
  global_dpd_->file2_close(&mu1);

  /*** L1 * MuBAR * X1 + L1 * MuBAR * X2 ***/

  global_dpd_->file2_init(&z1, PSIF_CC_TMP0, 0, 0, 1, "z_IA");

  sprintf(lbl, "%sBAR_MI", pert_c);
  global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep_c, 0, 0, lbl);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega);
  global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
  global_dpd_->contract222(&mu1, &X1, &z1, 1, 1, -1, 0);
  global_dpd_->file2_close(&X1);
  global_dpd_->file2_close(&mu1);

  sprintf(lbl, "%sBAR_AE", pert_c);
  global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep_c, 1, 1, lbl);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega);
  global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
  global_dpd_->contract222(&X1, &mu1, &z1, 0, 0, 1, 1);
  global_dpd_->file2_close(&X1);
  global_dpd_->file2_close(&mu1);

  sprintf(lbl, "%sBAR_ME", pert_c);
  global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep_c, 0, 1, lbl);
  sprintf(lbl, "X_%s_(2IjAb-IjbA) (%5.3f)", pert_x, omega);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
  global_dpd_->dot24(&mu1, &X2, &z1, 0, 0, 1, 1);
  global_dpd_->buf4_close(&X2);
  global_dpd_->file2_close(&mu1);

  global_dpd_->file2_init(&l1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1");
  polar += 2.0 * global_dpd_->file2_dot(&z1, &l1);
  global_dpd_->file2_close(&l1);

  global_dpd_->file2_close(&z1);

  /*** L2 * MuBAR * X1 + L2 * MuBAR * X2 ***/

  global_dpd_->file2_init(&xc, PSIF_CC_TMP0, 0, 0, 0, "XC_IJ");
  sprintf(lbl, "%s_IA", pert_c);
  global_dpd_->file2_init(&mu, PSIF_CC_OEI, irrep_c, 0, 1, lbl);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega);
  global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
  global_dpd_->contract222(&X1, &mu, &xc, 0, 0, 2.0, 0.0);
  global_dpd_->file2_close(&X1);
  global_dpd_->file2_close(&mu);
  global_dpd_->file2_close(&xc);

  global_dpd_->file2_init(&lt, PSIF_CC_OEI, 0, 0, 0, "Lt_IJ");
  global_dpd_->file2_init(&xc, PSIF_CC_TMP0, 0, 0, 0, "XC_IJ");
  polar += global_dpd_->file2_dot(&lt, &xc);
  global_dpd_->file2_close(&xc);
  global_dpd_->file2_close(&lt);

  global_dpd_->buf4_init(&z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");
  global_dpd_->buf4_scm(&z2, 0);
  global_dpd_->buf4_close(&z2);

  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega);
  global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);

  global_dpd_->buf4_init(&z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) Final");

  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
  sprintf(lbl, "%sBAR_MbIj", pert_c);
  global_dpd_->buf4_init(&mu2, PSIF_CC_LR, irrep_c, 10, 0, 10, 0, 0, lbl);
  global_dpd_->contract244(&X1, &mu2, &Z, 0, 0, 1, 1, 0);
  global_dpd_->buf4_close(&mu2);
  global_dpd_->buf4_axpy(&Z, &z2, -1);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
  global_dpd_->buf4_axpy(&Z, &z2, -1);
  global_dpd_->buf4_close(&Z);

  global_dpd_->file2_close(&X1);

  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_x, omega);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);

  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");

  sprintf(lbl, "%sBAR_AE", pert_c);
  global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep_c, 1, 1, lbl);
  global_dpd_->contract424(&X2, &mu1, &Z, 3, 1, 0, 1, 0);
  global_dpd_->file2_close(&mu1);
  global_dpd_->buf4_axpy(&Z, &z2, 1);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
  global_dpd_->buf4_axpy(&Z, &z2, 1);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");

  sprintf(lbl, "%sBAR_MI", pert_c);
  global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep_c, 0, 0, lbl);
  global_dpd_->contract244(&mu1, &X2, &Z, 0, 0, 0, 1, 0);
  global_dpd_->file2_close(&mu1);
  global_dpd_->buf4_axpy(&Z, &z2, -1);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, qpsr, 0, 5, "Z(jI,bA)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 0, 5, 0, 5, 0, "Z(jI,bA)");
  global_dpd_->buf4_axpy(&Z, &z2, -1);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_close(&X2);

  global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
  polar += global_dpd_->buf4_dot(&l2, &z2);
  global_dpd_->buf4_close(&l2);

  global_dpd_->buf4_close(&z2);

  if(params.sekino) {  /* disconnected piece for Sekino-Bartlett modelIII */
    /* L2 * MUBAR * X1 */
    sprintf(lbl, "%sZ_IA", pert_c);
    global_dpd_->file2_init(&z1, PSIF_CC_TMP0, irrep_c, 0, 1, lbl);
    sprintf(lbl, "%sBAR_IA", pert_c);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep_c, 0, 1, lbl);
    global_dpd_->buf4_init(&l2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->dot24(&mu1, &l2, &z1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&l2);
    global_dpd_->file2_close(&mu1);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
    polar += 2.0 * global_dpd_->file2_dot(&X1, &z1);
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&z1);
  }

  return polar;
}

}} // namespace psi::ccresponse
