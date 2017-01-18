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
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

double converged(const char *pert, int irrep, double omega)
{
  dpdfile2 X1, X1new;
  dpdbuf4 X2, X2new;
  double rms=0.0, value;
  int row, col, h, nirreps;
  char lbl[32];

  nirreps = moinfo.nirreps;

  sprintf(lbl, "New X_%s_IA (%5.3f)", pert, omega);
  global_dpd_->file2_init(&X1new, PSIF_CC_OEI, irrep, 0, 1, lbl);
  global_dpd_->file2_mat_init(&X1new);
  global_dpd_->file2_mat_rd(&X1new);
  sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
  global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
  global_dpd_->file2_mat_init(&X1);
  global_dpd_->file2_mat_rd(&X1);

  for(h=0; h < nirreps; h++)
    for(row=0; row < X1.params->rowtot[h]; row++)
      for(col=0; col < X1.params->coltot[h^irrep]; col++) {
	value = X1new.matrix[h][row][col] - X1.matrix[h][row][col];
	rms += value * value;
      }
  global_dpd_->file2_mat_close(&X1new);
  global_dpd_->file2_close(&X1new);
  global_dpd_->file2_mat_close(&X1);
  global_dpd_->file2_close(&X1);

  sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_init(&X2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
  sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
  global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&X2new, h);
    global_dpd_->buf4_mat_irrep_rd(&X2new, h);
    global_dpd_->buf4_mat_irrep_init(&X2, h);
    global_dpd_->buf4_mat_irrep_rd(&X2, h);

    for(row=0; row < X2.params->rowtot[h]; row++)
      for(col=0; col < X2.params->coltot[h^irrep]; col++) {
	value = X2new.matrix[h][row][col] - X2.matrix[h][row][col];
	rms += value * value;
      }

    global_dpd_->buf4_mat_irrep_close(&X2new, h);
    global_dpd_->buf4_mat_irrep_close(&X2, h);
  }
  global_dpd_->buf4_close(&X2new);
  global_dpd_->buf4_close(&X2);

  return sqrt(rms);
}

}} // namespace psi::ccresponse
