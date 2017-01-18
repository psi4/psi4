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
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

void denom1(dpdfile2 *X1, double omega)
{
  int nirreps, h, irrep;
  int i, a;
  int *occpi, *virtpi;
  dpdfile2 FAE, FMI;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi;
  virtpi = moinfo.virtpi;

  irrep = X1->my_irrep;

  if(params.wfn == "CC2" || params.wfn == "EOM_CC2") {
    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->file2_mat_init(&FMI);
    global_dpd_->file2_mat_rd(&FMI);

    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_mat_init(&FAE);
    global_dpd_->file2_mat_rd(&FAE);
  }
  else {
    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->file2_mat_init(&FAE);
    global_dpd_->file2_mat_rd(&FAE);

    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->file2_mat_init(&FMI);
    global_dpd_->file2_mat_rd(&FMI);
  }

  global_dpd_->file2_mat_init(X1);
  global_dpd_->file2_mat_rd(X1);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < occpi[h]; i++)
      for(a=0; a < virtpi[h^irrep]; a++)
	X1->matrix[h][i][a] /= (FMI.matrix[h][i][i] - FAE.matrix[h^irrep][a][a] + omega);
  }
  global_dpd_->file2_mat_wrt(X1);
  global_dpd_->file2_mat_close(X1);

  global_dpd_->file2_mat_close(&FAE);
  global_dpd_->file2_mat_close(&FMI);
  global_dpd_->file2_close(&FAE);
  global_dpd_->file2_close(&FMI);
}

void denom2(dpdbuf4 *X2, double omega)
{
  int nirreps, h, row, col, irrep;
  int i, j, I, J, a, b, A, B, isym, jsym, asym, bsym;
  dpdfile2 FAE, FMI;

  nirreps = moinfo.nirreps;
  irrep = X2->file.my_irrep;

  if(params.wfn == "CC2" || params.wfn == "EOM_CC2") {
    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->file2_mat_init(&FMI);
    global_dpd_->file2_mat_rd(&FMI);

    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_mat_init(&FAE);
    global_dpd_->file2_mat_rd(&FAE);
  }
  else {
    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->file2_mat_init(&FAE);
    global_dpd_->file2_mat_rd(&FAE);

    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->file2_mat_init(&FMI);
    global_dpd_->file2_mat_rd(&FMI);
  }

  for(h=0; h < nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(X2, h);
    global_dpd_->buf4_mat_irrep_rd(X2, h);

    for(row=0; row < X2->params->rowtot[h]; row++) {

      i = X2->params->roworb[h][row][0];
      j = X2->params->roworb[h][row][1];
      isym = X2->params->psym[i];
      jsym = X2->params->qsym[j];
      I = i - moinfo.occ_off[isym];
      J = j - moinfo.occ_off[jsym];

      for(col=0; col < X2->params->coltot[h^irrep]; col++) {

	a = X2->params->colorb[h^irrep][col][0];
	b = X2->params->colorb[h^irrep][col][1];
	asym = X2->params->rsym[a];
	bsym = X2->params->ssym[b];
	A = a - moinfo.vir_off[asym];
	B = b - moinfo.vir_off[bsym];

	X2->matrix[h][row][col] /= (FMI.matrix[isym][I][I] + FMI.matrix[jsym][J][J] -
				    FAE.matrix[asym][A][A] - FAE.matrix[bsym][B][B] + omega);

      }
    }

    global_dpd_->buf4_mat_irrep_wrt(X2, h);
    global_dpd_->buf4_mat_irrep_close(X2, h);
  }

  global_dpd_->file2_mat_close(&FAE);
  global_dpd_->file2_mat_close(&FMI);
  global_dpd_->file2_close(&FAE);
  global_dpd_->file2_close(&FMI);
}

}} // namespace psi::ccresponse
