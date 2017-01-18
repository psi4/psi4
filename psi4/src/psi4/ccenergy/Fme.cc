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
#include "ccwave.h"

namespace psi { namespace ccenergy {

void CCEnergyWavefunction::Fme_build(void)
{
  dpdfile2 FME, Fme, fIA, fia, tIA, tia;
  dpdbuf4 D_anti, D;

  if(params_.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_copy(&fIA, PSIF_CC_OEI, "FME");
    global_dpd_->file2_close(&fIA);

    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");

    global_dpd_->buf4_init(&D_anti, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    global_dpd_->dot13(&tIA, &D_anti, &FME, 0, 0, 1.0, 1.0);
    global_dpd_->dot13(&tIA, &D, &FME, 0, 0, 1.0, 1.0);

    global_dpd_->file2_close(&tIA);
    global_dpd_->buf4_close(&D_anti);
    global_dpd_->buf4_close(&D);

    global_dpd_->file2_close(&FME);
  }
  else if(params_.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_copy(&fIA, PSIF_CC_OEI, "FME");
    global_dpd_->file2_close(&fIA);

    global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
    global_dpd_->file2_copy(&fia, PSIF_CC_OEI, "Fme");
    global_dpd_->file2_close(&fia);

    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");

    global_dpd_->buf4_init(&D_anti, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->dot13(&tIA, &D_anti, &FME, 0, 0, 1.0, 1.0);
    global_dpd_->dot13(&tia, &D, &FME, 0, 0, 1.0, 1.0);

    global_dpd_->dot13(&tia, &D_anti, &Fme, 0, 0, 1.0, 1.0);
    global_dpd_->dot13(&tIA, &D, &Fme, 0, 0, 1.0, 1.0);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);
    global_dpd_->buf4_close(&D_anti);
    global_dpd_->buf4_close(&D);

    global_dpd_->file2_close(&FME);
    global_dpd_->file2_close(&Fme);
  }
  else if(params_.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    global_dpd_->file2_copy(&fIA, PSIF_CC_OEI, "FME");
    global_dpd_->file2_close(&fIA);

    global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
    global_dpd_->file2_copy(&fia, PSIF_CC_OEI, "Fme");
    global_dpd_->file2_close(&fia);

    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    global_dpd_->contract422(&D, &tIA, &FME, 0, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    global_dpd_->contract422(&D, &tia, &FME, 0, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->contract422(&D, &tia, &Fme, 0, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    global_dpd_->contract422(&D, &tIA, &Fme, 0, 0, 1, 1);
    global_dpd_->buf4_close(&D);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    global_dpd_->file2_close(&FME);
    global_dpd_->file2_close(&Fme);

  }
}
}} // namespace psi::ccenergy
