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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

int converged(int L_irr)
{
  int row,col,h,nirreps;
  double rms=0.0;
  dpdfile2 L1, L1old;
  dpdbuf4 L2, L2old;

  nirreps = moinfo.nirreps;

  dpd_file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "New LIA");
  dpd_file2_mat_init(&L1);
  dpd_file2_mat_rd(&L1);
  dpd_file2_init(&L1old, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
  dpd_file2_mat_init(&L1old);
  dpd_file2_mat_rd(&L1old);

  for(h=0; h < nirreps; h++)
    for(row=0; row < L1.params->rowtot[h]; row++)
      for(col=0; col < L1.params->coltot[h^L_irr]; col++)
	rms += (L1.matrix[h][row][col] - L1old.matrix[h][row][col]) *
	  (L1.matrix[h][row][col] - L1old.matrix[h][row][col]);

  dpd_file2_mat_close(&L1);
  dpd_file2_close(&L1);
  dpd_file2_mat_close(&L1old);
  dpd_file2_close(&L1old);

  if(params.ref == 0) rms *= 2.0;

  if(params.ref == 1) { /** ROHF **/

    dpd_file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "New Lia");
    dpd_file2_mat_init(&L1);
    dpd_file2_mat_rd(&L1);
    dpd_file2_init(&L1old, PSIF_CC_LAMBDA, L_irr, 0, 1, "Lia");
    dpd_file2_mat_init(&L1old);
    dpd_file2_mat_rd(&L1old);

  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 2, 3, "New Lia");
    dpd_file2_mat_init(&L1);
    dpd_file2_mat_rd(&L1);
    dpd_file2_init(&L1old, PSIF_CC_LAMBDA, L_irr, 2, 3, "Lia");
    dpd_file2_mat_init(&L1old);
    dpd_file2_mat_rd(&L1old);

  }

  if(params.ref == 1 || params.ref == 2) {
    for(h=0; h < nirreps; h++)
      for(row=0; row < L1.params->rowtot[h]; row++)
	for(col=0; col < L1.params->coltot[h^L_irr]; col++)
	  rms += (L1.matrix[h][row][col] - L1old.matrix[h][row][col]) *
	    (L1.matrix[h][row][col] - L1old.matrix[h][row][col]);

    dpd_file2_mat_close(&L1);
    dpd_file2_close(&L1);
    dpd_file2_mat_close(&L1old);
    dpd_file2_close(&L1old);
  }

  if(params.ref == 1 || params.ref == 2) { 
    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_init(&L2old, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2, h);
      dpd_buf4_mat_irrep_rd(&L2, h);
      dpd_buf4_mat_irrep_init(&L2old, h);
      dpd_buf4_mat_irrep_rd(&L2old, h);
      for(row=0; row < L2.params->rowtot[h]; row++)
	for(col=0; col < L2.params->coltot[h^L_irr]; col++)
	  rms += (L2.matrix[h][row][col] - L2old.matrix[h][row][col]) *
	    (L2.matrix[h][row][col] - L2old.matrix[h][row][col]);
      dpd_buf4_mat_irrep_close(&L2, h);
      dpd_buf4_mat_irrep_close(&L2old, h);
    }
    dpd_buf4_close(&L2old);
    dpd_buf4_close(&L2);
  }

  if(params.ref == 1) { /** ROHF **/
    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    dpd_buf4_init(&L2old, PSIF_CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    dpd_buf4_init(&L2old, PSIF_CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
  }

  if(params.ref == 1 || params.ref == 2) {
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2, h);
      dpd_buf4_mat_irrep_rd(&L2, h);
      dpd_buf4_mat_irrep_init(&L2old, h);
      dpd_buf4_mat_irrep_rd(&L2old, h);
      for(row=0; row < L2.params->rowtot[h]; row++)
	for(col=0; col < L2.params->coltot[h^L_irr]; col++)
	  rms += (L2.matrix[h][row][col] - L2old.matrix[h][row][col]) *
	    (L2.matrix[h][row][col] - L2old.matrix[h][row][col]);
      dpd_buf4_mat_irrep_close(&L2, h);
      dpd_buf4_mat_irrep_close(&L2old, h);
    }
    dpd_buf4_close(&L2old);
    dpd_buf4_close(&L2);
  }

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_init(&L2old, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
  }
  else if(params.ref == 2) { /** UHF **/
    dpd_buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_buf4_init(&L2old, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
  }

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&L2, h);
    dpd_buf4_mat_irrep_rd(&L2, h);
    dpd_buf4_mat_irrep_init(&L2old, h);
    dpd_buf4_mat_irrep_rd(&L2old, h);
    for(row=0; row < L2.params->rowtot[h]; row++)
      for(col=0; col < L2.params->coltot[h^L_irr]; col++)
	rms += (L2.matrix[h][row][col] - L2old.matrix[h][row][col]) *
	  (L2.matrix[h][row][col] - L2old.matrix[h][row][col]);
    dpd_buf4_mat_irrep_close(&L2, h);
    dpd_buf4_mat_irrep_close(&L2old, h);
  }
  dpd_buf4_close(&L2old);
  dpd_buf4_close(&L2);

  rms = sqrt(rms);
  moinfo.conv = rms;

  if(rms < params.convergence) return 1;
  else return 0;
}

}} // namespace psi::cclambda
