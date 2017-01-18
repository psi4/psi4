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
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

/* lambda_residuals(): Computes the "lambda residual" contributions to
** the <0|L*(HBAR*X1*Y1)c|0> part of the response function.
**
** Type-I residuals are taken directly from cclambda: L2*D, L2*Wmnij,
** L2*Wabef, L2*Wamef, L2*Wmnie, and L2*three-body-terms.  These are
** already computed by the L2 code in cclambda (see cclambda/L2.c).
**
** Type-II residuals are computed once in lambda_residuals(): L1*Fme
** and L2*Wmbej.  These are computed explicitly here with a minimum
** number of contractions.
**
** TDC, 9/10/05
*/

void lambda_residuals(void)
{
  dpdbuf4 L2, Z, W;
  dpdfile2 L1, F;
  int h, i, e, m, a, I, E, M, A, Isym, Esym, Msym, Asym;
  int row, col;

  /* Generate spin-adapted Type-I Lambda-residual contributions to LHX1Y1 */
  global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "LHX1Y1 Residual I");
  global_dpd_->buf4_scmcopy(&L2, PSIF_CC_LAMPS, "LHX1Y1 I (2 Lijab - Lijba)", 2);
  global_dpd_->buf4_sort_axpy(&L2, PSIF_CC_LAMPS, pqsr, 0, 5, "LHX1Y1 I (2 Lijab - Lijba)", -1);
  global_dpd_->buf4_close(&L2);

  /* Generate spin-adapted Type-II Lambda-residual contributions to LHX1Y1 */
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "LHX1Y1 Residual II");

  global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1");
  global_dpd_->file2_mat_init(&L1);
  global_dpd_->file2_mat_rd(&L1);
  global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->file2_mat_init(&F);
  global_dpd_->file2_mat_rd(&F);

  for(h=0; h < moinfo.nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&Z, h);
    for(row=0; row < Z.params->rowtot[h]; row++) {
      i = Z.params->roworb[h][row][0];
      a = Z.params->roworb[h][row][1];
      for(col=0; col < Z.params->coltot[h]; col++) {
	m = Z.params->colorb[h][col][0];
	e = Z.params->colorb[h][col][1];
	I = L1.params->rowidx[i]; Isym = L1.params->psym[i];
	A = L1.params->colidx[a]; Asym = L1.params->qsym[a];
	M = F.params->rowidx[m]; Msym = F.params->psym[m];
	E = F.params->colidx[e]; Esym = F.params->qsym[e];
	if((Isym == Asym) && (Msym == Esym))
	  Z.matrix[h][row][col] = (L1.matrix[Isym][I][A] * F.matrix[Msym][M][E]);
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&Z, h);
    global_dpd_->buf4_mat_irrep_close(&Z, h);
  }
  global_dpd_->file2_mat_close(&F);
  global_dpd_->file2_close(&F);
  global_dpd_->file2_mat_close(&L1);
  global_dpd_->file2_close(&L1);

  global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
  global_dpd_->buf4_sort(&L2, PSIF_CC_TMP0, prqs, 10, 10, "2 Lijab - Lijba (ia,jb)");
  global_dpd_->buf4_sort(&L2, PSIF_CC_TMP0, psqr, 10, 10, "2 Lijab - Lijba (ib,ja)");
  global_dpd_->buf4_close(&L2);

  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj");
  global_dpd_->buf4_init(&L2, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "2 Lijab - Lijba (ia,jb)");
  global_dpd_->contract444(&L2, &W, &Z, 0, 0, 1, 1);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_close(&W);

  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
  global_dpd_->buf4_init(&L2, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "2 Lijab - Lijba (ib,ja)");
  global_dpd_->contract444(&L2, &W, &Z, 0, 0, -1, 1);
  global_dpd_->buf4_close(&L2);
  global_dpd_->buf4_close(&W);

  global_dpd_->buf4_sort(&Z, PSIF_CC_LAMPS, psrq, 10, 10, "LHX1Y1 Residual II");
  global_dpd_->buf4_close(&Z);
}

}} // namespace psi::ccresponse
