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
#include <cmath>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

int converged(double ediff)
{
  int row,col,h,nirreps;
  double rms=0.0;
  dpdfile2 T1, T1old;
  dpdbuf4 T2, T2old;

  nirreps = moinfo.nirreps;

  if(params.ref == 0) { /** RHF **/

    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    dpd_->file2_mat_init(&T1);
    dpd_->file2_mat_rd(&T1);
    dpd_->file2_init(&T1old, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_mat_init(&T1old);
    dpd_->file2_mat_rd(&T1old);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1.params->rowtot[h]; row++)
	for(col=0; col < T1.params->coltot[h]; col++)
	  rms += (T1.matrix[h][row][col] - T1old.matrix[h][row][col]) *
	    (T1.matrix[h][row][col] - T1old.matrix[h][row][col]);

    dpd_->file2_mat_close(&T1);
    dpd_->file2_close(&T1);
    dpd_->file2_mat_close(&T1old);
    dpd_->file2_close(&T1old);

    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_->buf4_mat_irrep_init(&T2, h);
      dpd_->buf4_mat_irrep_rd(&T2, h);
      dpd_->buf4_mat_irrep_init(&T2old, h);
      dpd_->buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      dpd_->buf4_mat_irrep_close(&T2, h);
      dpd_->buf4_mat_irrep_close(&T2old, h);
    }
    dpd_->buf4_close(&T2old);
    dpd_->buf4_close(&T2);

  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    dpd_->file2_mat_init(&T1);
    dpd_->file2_mat_rd(&T1);
    dpd_->file2_init(&T1old, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_mat_init(&T1old);
    dpd_->file2_mat_rd(&T1old);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1.params->rowtot[h]; row++)
	for(col=0; col < T1.params->coltot[h]; col++)
	  rms += (T1.matrix[h][row][col] - T1old.matrix[h][row][col]) *
	    (T1.matrix[h][row][col] - T1old.matrix[h][row][col]);

    dpd_->file2_mat_close(&T1);
    dpd_->file2_close(&T1);
    dpd_->file2_mat_close(&T1old);
    dpd_->file2_close(&T1old);

    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "New tia");
    dpd_->file2_mat_init(&T1);
    dpd_->file2_mat_rd(&T1);
    dpd_->file2_init(&T1old, PSIF_CC_OEI, 0, 0, 1, "tia");
    dpd_->file2_mat_init(&T1old);
    dpd_->file2_mat_rd(&T1old);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1.params->rowtot[h]; row++)
	for(col=0; col < T1.params->coltot[h]; col++)
	  rms += (T1.matrix[h][row][col] - T1old.matrix[h][row][col]) *
	    (T1.matrix[h][row][col] - T1old.matrix[h][row][col]);

    dpd_->file2_mat_close(&T1);
    dpd_->file2_close(&T1);
    dpd_->file2_mat_close(&T1old);
    dpd_->file2_close(&T1old);

    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
    for(h=0; h < nirreps; h++) {
      dpd_->buf4_mat_irrep_init(&T2, h);
      dpd_->buf4_mat_irrep_rd(&T2, h);
      dpd_->buf4_mat_irrep_init(&T2old, h);
      dpd_->buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      dpd_->buf4_mat_irrep_close(&T2, h);
      dpd_->buf4_mat_irrep_close(&T2old, h);
    }
    dpd_->buf4_close(&T2old);
    dpd_->buf4_close(&T2);

    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tijab");
    dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
    for(h=0; h < nirreps; h++) {
      dpd_->buf4_mat_irrep_init(&T2, h);
      dpd_->buf4_mat_irrep_rd(&T2, h);
      dpd_->buf4_mat_irrep_init(&T2old, h);
      dpd_->buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      dpd_->buf4_mat_irrep_close(&T2, h);
      dpd_->buf4_mat_irrep_close(&T2old, h);
    }
    dpd_->buf4_close(&T2old);
    dpd_->buf4_close(&T2);

    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_->buf4_mat_irrep_init(&T2, h);
      dpd_->buf4_mat_irrep_rd(&T2, h);
      dpd_->buf4_mat_irrep_init(&T2old, h);
      dpd_->buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      dpd_->buf4_mat_irrep_close(&T2, h);
      dpd_->buf4_mat_irrep_close(&T2old, h);
    }
    dpd_->buf4_close(&T2old);
    dpd_->buf4_close(&T2);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    dpd_->file2_mat_init(&T1);
    dpd_->file2_mat_rd(&T1);
    dpd_->file2_init(&T1old, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_mat_init(&T1old);
    dpd_->file2_mat_rd(&T1old);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1.params->rowtot[h]; row++)
	for(col=0; col < T1.params->coltot[h]; col++)
	  rms += (T1.matrix[h][row][col] - T1old.matrix[h][row][col]) *
	    (T1.matrix[h][row][col] - T1old.matrix[h][row][col]);

    dpd_->file2_mat_close(&T1);
    dpd_->file2_close(&T1);
    dpd_->file2_mat_close(&T1old);
    dpd_->file2_close(&T1old);

    dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "New tia");
    dpd_->file2_mat_init(&T1);
    dpd_->file2_mat_rd(&T1);
    dpd_->file2_init(&T1old, PSIF_CC_OEI, 0, 2, 3, "tia");
    dpd_->file2_mat_init(&T1old);
    dpd_->file2_mat_rd(&T1old);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1.params->rowtot[h]; row++)
	for(col=0; col < T1.params->coltot[h]; col++)
	  rms += (T1.matrix[h][row][col] - T1old.matrix[h][row][col]) *
	    (T1.matrix[h][row][col] - T1old.matrix[h][row][col]);

    dpd_->file2_mat_close(&T1);
    dpd_->file2_close(&T1);
    dpd_->file2_mat_close(&T1old);
    dpd_->file2_close(&T1old);

    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
    dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
    for(h=0; h < nirreps; h++) {
      dpd_->buf4_mat_irrep_init(&T2, h);
      dpd_->buf4_mat_irrep_rd(&T2, h);
      dpd_->buf4_mat_irrep_init(&T2old, h);
      dpd_->buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      dpd_->buf4_mat_irrep_close(&T2, h);
      dpd_->buf4_mat_irrep_close(&T2old, h);
    }
    dpd_->buf4_close(&T2old);
    dpd_->buf4_close(&T2);

    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
    dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
    for(h=0; h < nirreps; h++) {
      dpd_->buf4_mat_irrep_init(&T2, h);
      dpd_->buf4_mat_irrep_rd(&T2, h);
      dpd_->buf4_mat_irrep_init(&T2old, h);
      dpd_->buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      dpd_->buf4_mat_irrep_close(&T2, h);
      dpd_->buf4_mat_irrep_close(&T2old, h);
    }
    dpd_->buf4_close(&T2old);
    dpd_->buf4_close(&T2);

    dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    dpd_->buf4_init(&T2old, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_->buf4_mat_irrep_init(&T2, h);
      dpd_->buf4_mat_irrep_rd(&T2, h);
      dpd_->buf4_mat_irrep_init(&T2old, h);
      dpd_->buf4_mat_irrep_rd(&T2old, h);
      for(row=0; row < T2.params->rowtot[h]; row++)
	for(col=0; col < T2.params->coltot[h]; col++)
	  rms += (T2.matrix[h][row][col] - T2old.matrix[h][row][col]) *
	    (T2.matrix[h][row][col] - T2old.matrix[h][row][col]);
      dpd_->buf4_mat_irrep_close(&T2, h);
      dpd_->buf4_mat_irrep_close(&T2old, h);
    }
    dpd_->buf4_close(&T2old);
    dpd_->buf4_close(&T2);
  }

  rms = sqrt(rms);
  moinfo.conv = rms;

  if((rms < params.convergence) && (fabs(ediff) < params.e_convergence)) return 1;
  else return 0;
}
}} // namespace psi::ccenergy
