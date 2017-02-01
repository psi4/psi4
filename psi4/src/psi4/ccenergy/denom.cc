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
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libdpd/dpd.h"
#include "Params.h"
#include "Local.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

void local_filter_T1(dpdfile2 *T1);
void dijabT2(void);

/* apply denominators to t1 and t2 */

void CCEnergyWavefunction::denom(void)
{
  dpdfile2 newtIA, dIA, tIA, newtia, dia, tia;

  if (params_.ref == 0) {
    global_dpd_->file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    global_dpd_->file2_copy(&newtIA, PSIF_CC_OEI, "New tIA Increment");
    global_dpd_->file2_close(&newtIA);

    global_dpd_->file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "New tIA Increment");
    if(params_.local && local_.filter_singles) {
      local_filter_T1(&newtIA);
    }
    else {
      global_dpd_->file2_init(&dIA, PSIF_CC_OEI, 0, 0, 1, "dIA");
      global_dpd_->file2_dirprd(&dIA, &newtIA);
      global_dpd_->file2_close(&dIA);
    }
    global_dpd_->file2_close(&newtIA);

    /* Add the new increment to the old tIA to get the New tIA */
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_copy(&tIA, PSIF_CC_OEI, "New tIA");
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "New tIA Increment");
    global_dpd_->file2_axpy(&tIA, &newtIA, 1, 0);
    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&newtIA);
  }
  else if (params_.ref == 1) {
    global_dpd_->file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    global_dpd_->file2_init(&dIA, PSIF_CC_OEI, 0, 0, 1, "dIA");
    global_dpd_->file2_dirprd(&dIA, &newtIA);
    global_dpd_->file2_close(&dIA);
    global_dpd_->file2_close(&newtIA);

    global_dpd_->file2_init(&newtia, PSIF_CC_OEI, 0, 0, 1, "New tia");
    global_dpd_->file2_init(&dia, PSIF_CC_OEI, 0, 0, 1, "dia");
    global_dpd_->file2_dirprd(&dia, &newtia);
    global_dpd_->file2_close(&dia);
    global_dpd_->file2_close(&newtia);
  }
  else if (params_.ref == 2) {
    global_dpd_->file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    global_dpd_->file2_init(&dIA, PSIF_CC_OEI, 0, 0, 1, "dIA");
    global_dpd_->file2_dirprd(&dIA, &newtIA);
    global_dpd_->file2_close(&dIA);
    global_dpd_->file2_close(&newtIA);

    global_dpd_->file2_init(&newtia, PSIF_CC_OEI, 0, 2, 3, "New tia");
    global_dpd_->file2_init(&dia, PSIF_CC_OEI, 0, 2, 3, "dia");
    global_dpd_->file2_dirprd(&dia, &newtia);
    global_dpd_->file2_close(&dia);
    global_dpd_->file2_close(&newtia);
  }

  dijabT2();

  return;
}
}} // namespace psi::ccenergy
