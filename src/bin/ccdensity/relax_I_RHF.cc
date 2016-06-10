/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* relax_I_RHF(): Add the ROHF orbital-response contributions from
** the one-electron density matrix to the I(I,J) and I(I,A) blocks of
** the Lagrangian.  These terms arise from the first-order CPHF
** equations.  I *think* the following code is general enough to deal
** with both RHF and ROHF cases. */

void relax_I_RHF(void)
{
  dpdfile2 I, D, f;
  dpdbuf4 E;
  int h, nirreps, i, j, e, *occpi, *virtpi, *openpi;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi;
  virtpi = moinfo.virtpi;
  openpi = moinfo.openpi;

  /* I(I,A) = I'(I,A) - sum_M f(I,M) D(orb)(A,M) */
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I'IA");
  global_dpd_->file2_copy(&I, PSIF_CC_OEI, "I(I,A)");
  global_dpd_->file2_close(&I);
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 1, "I(I,A)");
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  global_dpd_->file2_init(&f, PSIF_CC_OEI, 0, 0, 0, "fIJ");
  global_dpd_->contract222(&f, &D, &I, 0, 0, -1.0, 1.0);
  global_dpd_->file2_close(&f);
  global_dpd_->file2_close(&D);
  global_dpd_->file2_close(&I);

  /* RHF Case: I(i,j) = I'(i,j) - D(orb)(e,c) [4 <ei|mj> - <ei|jm> - <ej|im>] */
  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I'IJ");
  global_dpd_->file2_copy(&I, PSIF_CC_OEI, "I(I,J)");
  global_dpd_->file2_close(&I);

  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  global_dpd_->buf4_scmcopy(&E, PSIF_CC_EINTS, "4 <ei|mj> - <ei|jm> - <ej|im>", 4);
  global_dpd_->buf4_sort_axpy(&E, PSIF_CC_EINTS, pqsr, 11, 0, "4 <ei|mj> - <ei|jm> - <ej|im>", -1);
  global_dpd_->buf4_sort_axpy(&E, PSIF_CC_EINTS, psqr, 11, 0, "4 <ei|mj> - <ei|jm> - <ej|im>", -1);
  global_dpd_->buf4_close(&E);

  global_dpd_->file2_init(&I, PSIF_CC_OEI, 0, 0, 0, "I(I,J)");
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "4 <ei|mj> - <ei|jm> - <ej|im>");
  global_dpd_->dot13(&D, &E, &I, 0, 0, -1.0, 1.0);
  global_dpd_->buf4_close(&E);
  global_dpd_->file2_close(&D);

  global_dpd_->file2_close(&I);
}

}} // namespace psi::ccdensity