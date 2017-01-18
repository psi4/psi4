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
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libdpd/dpd.h"
#include "Params.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

void CCEnergyWavefunction::init_amps(void)
{
  dpdfile2 tIA, tia, fIA, fia, dIA, dia;
  dpdbuf4 tIJAB, tijab, tIjAb, D, dIJAB, dijab, dIjAb;

  if(params_.ref == 0) { /** RHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    if(!params_.restart || !psio_tocscan(PSIF_CC_OEI, "tIA"))
      global_dpd_->file2_scm(&tIA, 0);
    else outfile->Printf( "    Using old T1 amplitudes.\n");
    global_dpd_->file2_close(&tIA);

    if(!params_.restart || !psio_tocscan(PSIF_CC_TAMPS, "tIjAb")) {
      global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      global_dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "tIjAb");
      global_dpd_->buf4_close(&D);

      global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      if(params_.local) {
	local_filter_T2(&tIjAb);
      }
      else {
	global_dpd_->buf4_init(&dIjAb, PSIF_CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
	global_dpd_->buf4_dirprd(&dIjAb, &tIjAb);
	global_dpd_->buf4_close(&dIjAb);
      }
      global_dpd_->buf4_close(&tIjAb);
    }
    else outfile->Printf( "    Using old T2 amplitudes.\n\n");
  }
  else if(params_.ref == 1) { /*** ROHF ***/
    if(!params_.restart || !psio_tocscan(PSIF_CC_OEI, "tIA") ||
       !psio_tocscan(PSIF_CC_OEI, "tia")) {

      global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
      global_dpd_->file2_copy(&fIA, PSIF_CC_OEI, "tIA");
      global_dpd_->file2_close(&fIA);
      global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
      /*  dpd_oe_scm(&tIA, 0);  */
      global_dpd_->file2_close(&tIA);

      global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
      global_dpd_->file2_copy(&fia, PSIF_CC_OEI, "tia");
      global_dpd_->file2_close(&fia);
      global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
      /*  dpd_oe_scm(&tia, 0); */
      global_dpd_->file2_close(&tia);

      global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
      global_dpd_->file2_init(&dIA, PSIF_CC_OEI, 0, 0, 1, "dIA");
      global_dpd_->file2_dirprd(&dIA, &tIA);
      global_dpd_->file2_close(&tIA);
      global_dpd_->file2_close(&dIA);

      global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
      global_dpd_->file2_init(&dia, PSIF_CC_OEI, 0, 0, 1, "dia");
      global_dpd_->file2_dirprd(&dia, &tia);
      global_dpd_->file2_close(&tia);
      global_dpd_->file2_close(&dia);
    }
    else outfile->Printf( "    Using old T1 amplitudes.\n");

    if(!params_.restart || !psio_tocscan(PSIF_CC_TAMPS, "tIjAb") ||
       !psio_tocscan(PSIF_CC_TAMPS, "tIJAB") || !psio_tocscan(PSIF_CC_TAMPS, "tijab")) {

      global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
      global_dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "tIJAB");
      global_dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "tijab");
      global_dpd_->buf4_close(&D);

      global_dpd_->buf4_init(&dIJAB, PSIF_CC_DENOM, 0, 1, 6, 1, 6, 0, "dIJAB");
      global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
      global_dpd_->buf4_dirprd(&dIJAB, &tIJAB);
      global_dpd_->buf4_close(&tIJAB);
      global_dpd_->buf4_close(&dIJAB);

      global_dpd_->buf4_init(&dijab, PSIF_CC_DENOM, 0, 1, 6, 1, 6, 0, "dijab");
      global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
      global_dpd_->buf4_dirprd(&dijab, &tijab);
      global_dpd_->buf4_close(&tijab);
      global_dpd_->buf4_close(&dijab);

      global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
      global_dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "tIjAb");
      global_dpd_->buf4_close(&D);

      global_dpd_->buf4_init(&dIjAb, PSIF_CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
      global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      global_dpd_->buf4_dirprd(&dIjAb, &tIjAb);
      global_dpd_->buf4_close(&tIjAb);
      global_dpd_->buf4_close(&dIjAb);
    }
    else
      outfile->Printf( "    Using old T2 amplitudes.\n");
  }
  else if(params_.ref == 2) { /*** UHF ***/

    if(!params_.restart || !psio_tocscan(PSIF_CC_OEI, "tIA") ||
       !psio_tocscan(PSIF_CC_OEI, "tia")) {

      global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
      global_dpd_->file2_copy(&fIA, PSIF_CC_OEI, "tIA");
      global_dpd_->file2_close(&fIA);
      global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
      /*    dpd_file2_scm(&tIA, 0); */
      global_dpd_->file2_close(&tIA);

      global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
      global_dpd_->file2_copy(&fia, PSIF_CC_OEI, "tia");
      global_dpd_->file2_close(&fia);
      global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
      /*    dpd_file2_scm(&tia, 0); */
      global_dpd_->file2_close(&tia);

      global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
      global_dpd_->file2_init(&dIA, PSIF_CC_OEI, 0, 0, 1, "dIA");
      global_dpd_->file2_dirprd(&dIA, &tIA);
      global_dpd_->file2_close(&tIA);
      global_dpd_->file2_close(&dIA);

      global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
      global_dpd_->file2_init(&dia, PSIF_CC_OEI, 0, 2, 3, "dia");
      global_dpd_->file2_dirprd(&dia, &tia);
      global_dpd_->file2_close(&tia);
      global_dpd_->file2_close(&dia);
    }
    else outfile->Printf( "    Using old T1 amplitudes.\n");

    if(!params_.restart || !psio_tocscan(PSIF_CC_TAMPS, "tIjAb") ||
       !psio_tocscan(PSIF_CC_TAMPS, "tIJAB") || !psio_tocscan(PSIF_CC_TAMPS, "tijab")) {

      global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
      global_dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "tIJAB");
      global_dpd_->buf4_close(&D);
      global_dpd_->buf4_init(&dIJAB, PSIF_CC_DENOM, 0, 1, 6, 1, 6, 0, "dIJAB");
      global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
      global_dpd_->buf4_dirprd(&dIJAB, &tIJAB);
      global_dpd_->buf4_close(&tIJAB);
      global_dpd_->buf4_close(&dIJAB);

      global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
      global_dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "tijab");
      global_dpd_->buf4_close(&D);
      global_dpd_->buf4_init(&dIJAB, PSIF_CC_DENOM, 0, 11, 16, 11, 16, 0, "dijab");
      global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
      global_dpd_->buf4_dirprd(&dIJAB, &tIJAB);
      global_dpd_->buf4_close(&tIJAB);
      global_dpd_->buf4_close(&dIJAB);

      global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
      global_dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "tIjAb");
      global_dpd_->buf4_close(&D);
      global_dpd_->buf4_init(&dIJAB, PSIF_CC_DENOM, 0, 22, 28, 22, 28, 0, "dIjAb");
      global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
      global_dpd_->buf4_dirprd(&dIJAB, &tIJAB);
      global_dpd_->buf4_close(&tIJAB);
      global_dpd_->buf4_close(&dIJAB);
    }
    else
      outfile->Printf( "    Using old T2 amplitudes.\n");
  }
  else {  /*** RHF/ROHF ***/

  }
}
}} // namespace psi::ccenergy
