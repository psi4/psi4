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

void denom2(dpdbuf4 *X2, double omega);
void local_filter_T2(dpdbuf4 *T2);

void cc2_X2_build(const char *pert, int irrep, double omega)
{
  dpdfile2 X1, z, F, t1;
  dpdbuf4 X2, X2new, Z, Z1, Z2, W, I;
  char lbl[32];

  sprintf(lbl, "%sBAR_IjAb", pert);
  global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_copy(&X2new, PSIF_CC_LR, lbl);
  global_dpd_->buf4_close(&X2new);
  global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

  /*** D-S ***/

  sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
  global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);

  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 10, 0, 10, 0, 0, "CC2 WMbIj");
  global_dpd_->contract244(&X1, &W, &Z, 0, 0, 1, 1, 0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_axpy(&Z, &X2new, -1);
  global_dpd_->buf4_close(&X2new);
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LR, qpsr, 0, 5, lbl, -1);
  global_dpd_->buf4_close(&Z);


  sprintf(lbl, "Z(Ab,Ij) %s", pert);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 5, 0, 5, 0, 0, lbl);
  global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 5, 11, 5, 11, 0, "CC2 WAbEi");
  global_dpd_->contract244(&X1, &W, &Z, 1, 2, 1, 1, 0);
  global_dpd_->buf4_close(&W);
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LR, rspq, 0, 5, lbl, 1);
  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, srqp, 0, 5, lbl);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->buf4_axpy(&Z, &X2new, 1);
  global_dpd_->buf4_close(&X2new);
  global_dpd_->buf4_close(&Z);

  global_dpd_->file2_close(&X1);

  /*** D-D ***/

  global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

  global_dpd_->buf4_axpy(&X2, &X2new, -omega);

  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "fAB");
  global_dpd_->contract424(&X2, &F, &Z, 3, 1, 0, 1, 0);
  global_dpd_->file2_close(&F);
  global_dpd_->buf4_axpy(&Z, &X2new, 1);
  sprintf(lbl, "Z(jI,bA) %s", pert);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qpsr, 0, 5, lbl);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->buf4_axpy(&Z, &X2new, 1);
  global_dpd_->buf4_close(&Z);

  sprintf(lbl, "Z(Ij,Ab) %s", pert);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->contract244(&F, &X2, &Z, 0, 0, 0, 1, 0);
  global_dpd_->file2_close(&F);
  global_dpd_->buf4_axpy(&Z, &X2new, -1);
  sprintf(lbl, "Z(jI,bA) %s", pert);
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qpsr, 0, 5, lbl);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->buf4_axpy(&Z, &X2new, -1);
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_close(&X2);

  /** Filter and apply denominator **/
  if(params.local) local_filter_T2(&X2new);
  else denom2(&X2new, omega);
  global_dpd_->buf4_close(&X2new);
}

}} // namespace psi::ccresponse
