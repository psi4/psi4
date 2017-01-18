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
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

int CCEnergyWavefunction::converged(double ediff)
{
  int row,col,h,nirreps;
  double rms=0.0;
  dpdfile2 T1, T1old;
  dpdbuf4 T2, T2old;

  nirreps = moinfo_.nirreps;

  if(params_.ref == 0) { /** RHF **/

    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    global_dpd_->file2_mat_init(&T1);
    global_dpd_->file2_mat_rd(&T1);
    global_dpd_->file2_init(&T1old, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&T1old);
    global_dpd_->file2_mat_rd(&T1old);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1.params->rowtot[h]; row++)
	for(col=0; col < T1.params->coltot[h]; col++)
	  rms += (T1.matrix[h][row][col] - T1old.matrix[h][row][col]) *
	    (T1.matrix[h][row][col] - T1old.matrix[h][row][col]);

    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_mat_close(&T1old);
    global_dpd_->file2_close(&T1old);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    global_dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&T2, h);
      global_dpd_->buf4_mat_irrep_rd(&T2, h);
      global_dpd_->buf4_mat_irrep_init(&T2old, h);
      global_dpd_->buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      global_dpd_->buf4_mat_irrep_close(&T2, h);
      global_dpd_->buf4_mat_irrep_close(&T2old, h);
    }
    global_dpd_->buf4_close(&T2old);
    global_dpd_->buf4_close(&T2);

  }
  else if(params_.ref == 1) { /** ROHF **/

//    outfile->Printf("I am a ROHF Wavefunction\n");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    global_dpd_->file2_mat_init(&T1);
    global_dpd_->file2_mat_rd(&T1);
    global_dpd_->file2_init(&T1old, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&T1old);
    global_dpd_->file2_mat_rd(&T1old);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1.params->rowtot[h]; row++)
	for(col=0; col < T1.params->coltot[h]; col++)
	  rms += (T1.matrix[h][row][col] - T1old.matrix[h][row][col]) *
	    (T1.matrix[h][row][col] - T1old.matrix[h][row][col]);

    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_mat_close(&T1old);
    global_dpd_->file2_close(&T1old);

    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "New tia");
    global_dpd_->file2_mat_init(&T1);
    global_dpd_->file2_mat_rd(&T1);
    global_dpd_->file2_init(&T1old, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->file2_mat_init(&T1old);
    global_dpd_->file2_mat_rd(&T1old);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1.params->rowtot[h]; row++)
	for(col=0; col < T1.params->coltot[h]; col++)
	  rms += (T1.matrix[h][row][col] - T1old.matrix[h][row][col]) *
	    (T1.matrix[h][row][col] - T1old.matrix[h][row][col]);

    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_mat_close(&T1old);
    global_dpd_->file2_close(&T1old);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    global_dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&T2, h);
      global_dpd_->buf4_mat_irrep_rd(&T2, h);
      global_dpd_->buf4_mat_irrep_init(&T2old, h);
      global_dpd_->buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      global_dpd_->buf4_mat_irrep_close(&T2, h);
      global_dpd_->buf4_mat_irrep_close(&T2old, h);
    }
    global_dpd_->buf4_close(&T2old);
    global_dpd_->buf4_close(&T2);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tijab");
    global_dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&T2, h);
      global_dpd_->buf4_mat_irrep_rd(&T2, h);
      global_dpd_->buf4_mat_irrep_init(&T2old, h);
      global_dpd_->buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      global_dpd_->buf4_mat_irrep_close(&T2, h);
      global_dpd_->buf4_mat_irrep_close(&T2old, h);
    }
    global_dpd_->buf4_close(&T2old);
    global_dpd_->buf4_close(&T2);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    global_dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&T2, h);
      global_dpd_->buf4_mat_irrep_rd(&T2, h);
      global_dpd_->buf4_mat_irrep_init(&T2old, h);
      global_dpd_->buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      global_dpd_->buf4_mat_irrep_close(&T2, h);
      global_dpd_->buf4_mat_irrep_close(&T2old, h);
    }
    global_dpd_->buf4_close(&T2old);
    global_dpd_->buf4_close(&T2);
  }
  else if(params_.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    global_dpd_->file2_mat_init(&T1);
    global_dpd_->file2_mat_rd(&T1);
    global_dpd_->file2_init(&T1old, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&T1old);
    global_dpd_->file2_mat_rd(&T1old);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1.params->rowtot[h]; row++)
	for(col=0; col < T1.params->coltot[h]; col++)
	  rms += (T1.matrix[h][row][col] - T1old.matrix[h][row][col]) *
	    (T1.matrix[h][row][col] - T1old.matrix[h][row][col]);

    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_mat_close(&T1old);
    global_dpd_->file2_close(&T1old);

    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "New tia");
    global_dpd_->file2_mat_init(&T1);
    global_dpd_->file2_mat_rd(&T1);
    global_dpd_->file2_init(&T1old, PSIF_CC_OEI, 0, 2, 3, "tia");
    global_dpd_->file2_mat_init(&T1old);
    global_dpd_->file2_mat_rd(&T1old);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1.params->rowtot[h]; row++)
	for(col=0; col < T1.params->coltot[h]; col++)
	  rms += (T1.matrix[h][row][col] - T1old.matrix[h][row][col]) *
	    (T1.matrix[h][row][col] - T1old.matrix[h][row][col]);

    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);
    global_dpd_->file2_mat_close(&T1old);
    global_dpd_->file2_close(&T1old);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    global_dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&T2, h);
      global_dpd_->buf4_mat_irrep_rd(&T2, h);
      global_dpd_->buf4_mat_irrep_init(&T2old, h);
      global_dpd_->buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      global_dpd_->buf4_mat_irrep_close(&T2, h);
      global_dpd_->buf4_mat_irrep_close(&T2old, h);
    }
    global_dpd_->buf4_close(&T2old);
    global_dpd_->buf4_close(&T2);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
    global_dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&T2, h);
      global_dpd_->buf4_mat_irrep_rd(&T2, h);
      global_dpd_->buf4_mat_irrep_init(&T2old, h);
      global_dpd_->buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      global_dpd_->buf4_mat_irrep_close(&T2, h);
      global_dpd_->buf4_mat_irrep_close(&T2old, h);
    }
    global_dpd_->buf4_close(&T2old);
    global_dpd_->buf4_close(&T2);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    global_dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&T2, h);
      global_dpd_->buf4_mat_irrep_rd(&T2, h);
      global_dpd_->buf4_mat_irrep_init(&T2old, h);
      global_dpd_->buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      global_dpd_->buf4_mat_irrep_close(&T2, h);
      global_dpd_->buf4_mat_irrep_close(&T2old, h);
    }
    global_dpd_->buf4_close(&T2old);
    global_dpd_->buf4_close(&T2);
  }

  rms = sqrt(rms);
  moinfo_.conv = rms;

  if((rms < params_.convergence) && (fabs(ediff) < params_.e_convergence)) return 1;
  else return 0;
}
}} // namespace psi::ccenergy
