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

void sort_X(const char *pert, int irrep, double omega)
{
  dpdbuf4 X;
  char lbl[32];

  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_init(&X, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  sprintf(lbl, "X_%s_IAjb (%5.3f)", pert, omega);
  global_dpd_->buf4_sort(&X, PSIF_CC_LR, prqs, 10, 10, lbl);
  sprintf(lbl, "X_%s_IbjA (%5.3f)", pert, omega);
  global_dpd_->buf4_sort(&X, PSIF_CC_LR, psqr, 10, 10, lbl);
  sprintf(lbl, "X_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
  global_dpd_->buf4_scmcopy(&X, PSIF_CC_LR, lbl, 2);
  global_dpd_->buf4_sort_axpy(&X, PSIF_CC_LR, pqsr, 0, 5, lbl, -1);
  global_dpd_->buf4_close(&X);

  sprintf(lbl, "X_%s_IAjb (%5.3f)", pert, omega);
  global_dpd_->buf4_init(&X, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
  sprintf(lbl, "X_%s_(2IAjb-IbjA) (%5.3f)", pert, omega);
  global_dpd_->buf4_scmcopy(&X, PSIF_CC_LR, lbl, 2);
  global_dpd_->buf4_sort_axpy(&X, PSIF_CC_LR, psrq, 10, 10, lbl, -1);
  global_dpd_->buf4_close(&X);

  sprintf(lbl, "X_%s_IAjb (%5.3f)", pert, omega);
  global_dpd_->buf4_init(&X, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
  sprintf(lbl, "X_%s_(2IAjb-jAIb) (%5.3f)", pert, omega);
  global_dpd_->buf4_scmcopy(&X, PSIF_CC_LR, lbl, 2);
  global_dpd_->buf4_sort_axpy(&X, PSIF_CC_LR, rqps, 10, 10, lbl, -1);
  global_dpd_->buf4_close(&X);

  if(params.ref == 0 && params.abcd == "NEW") {
    /* X(-)(ij,ab) (i>j, a>b) = X(ij,ab) - X(ij,ba) */
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X, PSIF_CC_LR, irrep, 4, 9, 0, 5, 1, lbl);
    sprintf(lbl, "X_%s_(-)(ij,ab) (%5.3f)", pert, omega);
    global_dpd_->buf4_copy(&X, PSIF_CC_LR, lbl);
    global_dpd_->buf4_close(&X);

    /* X(+)(ij,ab) (i>=j, a>=b) = X(ij,ab) + X(ij,ba) */
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    sprintf(lbl, "X_%s_(+)(ij,ab) (%5.3f)", pert, omega);
    global_dpd_->buf4_copy(&X, PSIF_CC_TMP0, lbl);
    global_dpd_->buf4_sort_axpy(&X, PSIF_CC_TMP0, pqsr, 0, 5, lbl, 1);
    global_dpd_->buf4_close(&X);
    global_dpd_->buf4_init(&X, PSIF_CC_TMP0, irrep, 3, 8, 0, 5, 0, lbl);
    global_dpd_->buf4_copy(&X, PSIF_CC_LR, lbl);
    global_dpd_->buf4_close(&X);
  }
}

}} // namespace psi::ccresponse
