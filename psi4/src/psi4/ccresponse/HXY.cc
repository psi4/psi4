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

/*! \defgroup ccresponse ccresponse: Coupled-cluster response module */

#include <cstdio>
#include <cstring>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

double HXY(const char *pert_x, int irrep_x, double omega_x,
	   const char *pert_y, int irrep_y, double omega_y)
{
  double polar;
  dpdfile2 X1, Y1, z;
  dpdbuf4 I;
  char lbl[32];

  sprintf(lbl, "Z_%s_IA", pert_y);
  global_dpd_->file2_init(&z, PSIF_CC_TMP0, irrep_y, 0, 1, lbl);

  sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
  global_dpd_->file2_init(&Y1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
  global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  global_dpd_->dot24(&Y1, &I, &z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&I);
  global_dpd_->file2_close(&Y1);

  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
  global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
  polar = 2.0 * global_dpd_->file2_dot(&X1, &z);
  global_dpd_->file2_close(&X1);

  global_dpd_->file2_close(&z);

  return polar;
}

}} // namespace psi::ccresponse
