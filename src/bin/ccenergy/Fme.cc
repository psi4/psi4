/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void Fme_build(void)
{
  dpdfile2 FME, Fme, fIA, fia, tIA, tia;
  dpdbuf4 D_anti, D;

  if(params.ref == 0) { /** RHF **/
    dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_->file2_copy(&fIA, PSIF_CC_OEI, "FME");
    dpd_->file2_close(&fIA);

    dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
  
    dpd_->buf4_init(&D_anti, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    dpd_->dot13(&tIA, &D_anti, &FME, 0, 0, 1.0, 1.0);
    dpd_->dot13(&tIA, &D, &FME, 0, 0, 1.0, 1.0);

    dpd_->file2_close(&tIA);
    dpd_->buf4_close(&D_anti);
    dpd_->buf4_close(&D);

    dpd_->file2_close(&FME);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_->file2_copy(&fIA, PSIF_CC_OEI, "FME");
    dpd_->file2_close(&fIA);

    dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
    dpd_->file2_copy(&fia, PSIF_CC_OEI, "Fme");
    dpd_->file2_close(&fia);
  
    dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");
  
    dpd_->buf4_init(&D_anti, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    dpd_->dot13(&tIA, &D_anti, &FME, 0, 0, 1.0, 1.0);
    dpd_->dot13(&tia, &D, &FME, 0, 0, 1.0, 1.0);

    dpd_->dot13(&tia, &D_anti, &Fme, 0, 0, 1.0, 1.0);
    dpd_->dot13(&tIA, &D, &Fme, 0, 0, 1.0, 1.0);

    dpd_->file2_close(&tIA);
    dpd_->file2_close(&tia);
    dpd_->buf4_close(&D_anti);
    dpd_->buf4_close(&D);

    dpd_->file2_close(&FME);
    dpd_->file2_close(&Fme);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
    dpd_->file2_copy(&fIA, PSIF_CC_OEI, "FME");
    dpd_->file2_close(&fIA);

    dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
    dpd_->file2_copy(&fia, PSIF_CC_OEI, "Fme");
    dpd_->file2_close(&fia);

    dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
  
    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    dpd_->contract422(&D, &tIA, &FME, 0, 0, 1, 1);
    dpd_->buf4_close(&D);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    dpd_->contract422(&D, &tia, &FME, 0, 0, 1, 1);
    dpd_->buf4_close(&D);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    dpd_->contract422(&D, &tia, &Fme, 0, 0, 1, 1);
    dpd_->buf4_close(&D);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    dpd_->contract422(&D, &tIA, &Fme, 0, 0, 1, 1);
    dpd_->buf4_close(&D);

    dpd_->file2_close(&tIA);
    dpd_->file2_close(&tia);

    dpd_->file2_close(&FME);
    dpd_->file2_close(&Fme);

  }
}
}} // namespace psi::ccenergy
