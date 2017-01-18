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

void update_X(const char *pert, int irrep, double omega)
{
  dpdfile2 X1new, X1;
  dpdbuf4 X2new, X2;
  char lbl[32];

  sprintf(lbl, "New X_%s_IA (%5.3f)", pert, omega);
  global_dpd_->file2_init(&X1new, PSIF_CC_OEI, irrep, 0, 1, lbl);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
  global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
  global_dpd_->file2_axpy(&X1, &X1new, 1, 0);
  global_dpd_->file2_close(&X1);
  global_dpd_->file2_close(&X1new);

  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  global_dpd_->buf4_axpy(&X2, &X2new, 1);
  global_dpd_->buf4_close(&X2);
  global_dpd_->buf4_close(&X2new);
}

}} // namespace psi::ccresponse
