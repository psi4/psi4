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
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

/* LHX1Y1(): Computes the <0|L*(HBAR*X1*Y1)c |0> contributions to the
** linear response function.  Note that most of the work is actually
** already done in other codes because the contractions may all be
** formulated in as products of X1*Y1 with L*HBAR, the latter of which
** are referred to in this code as "lambda residuals".
**
** Type-I residuals are taken directly from cclambda: L2*D, L2*Wmnij,
** L2*Wabef, L2*Wamef, L2*Wmnie, and L2*three-body-terms.
**
** Type-II residuals are computed once in lambda_residuals(): L1*Fme
** and L2*Wmbej.
**
** TDC, 9/10/05
*/

void build_XY(const char *, int, double, const char *, int, double);

double LHX1Y1(const char *pert_x, int irrep_x, double omega_x,
	      const char *pert_y, int irrep_y, double omega_y)
{

  dpdfile2 F, X1, Y1, Zmi, Zae_1, Zae_2, Zfb, Znj, ZIA, L1, t1, z;
  dpdbuf4 Z1, Z2, I, tau, W1, W2, ZIjAb, L2, T2, W, Z;
  double polar, polar_I, polar_II;
  char lbl[32];
  int Gbm, Gfe, bm, b, m, Gb, Gm, Ge, Gf, B, M, fe, f, e, ef, nrows, ncols;
  double *X;

  build_XY(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y);

  /* Type-I L2 residual */
  global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "LHX1Y1 I (2 Lijab - Lijba)");
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X*Y(ij,ab)");
  polar_I = 2.0 * global_dpd_->buf4_dot(&L2, &Z1);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&L2);

  /* Type-II L2 residual */
  global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "LHX1Y1 Residual II");
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "(X*Y+Y*X)(ie,ma)");
  polar_II = -2.0 * global_dpd_->buf4_dot(&L2, &Z1);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&L2);

  return polar_I+polar_II;
}

/* build_XY(): Compute products of X1 and Y1 for the
** <0|L*(HBAR*X1*Y1)|0> part of the response function.
**
** For the Type-I residuals:
** X*Y(ij,ab) = X(i,a) * Y(j,b)
**
** For the Type-II residuals:
** (X*Y+Y*X)(ie,ma) = [X(i,e) * Y(m,a) + X(m,a) * Y(i,e)]
**
** TDC, 9/10/05
**/

void build_XY(const char *pert_x, int irrep_x, double omega_x,
              const char *pert_y, int irrep_y, double omega_y)
{
  int h, row, col, i, j, m, e, f, a, I, J, M, E, F, A, ij, ef;
  int Isym, Jsym, Msym, Esym, Fsym, Asym;
  int nirreps;
  dpdbuf4 Z, Z1;
  dpdfile2 X1, Y1;
  char lbl[32];

  nirreps = moinfo.nirreps;

  sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
  global_dpd_->file2_init(&Y1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
  global_dpd_->file2_mat_init(&Y1);
  global_dpd_->file2_mat_rd(&Y1);

  sprintf(lbl, "X_%s_IA (%5.3f)", pert_x, omega_x);
  global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
  global_dpd_->file2_mat_init(&X1);
  global_dpd_->file2_mat_rd(&X1);

  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X*Y(ij,ab)");
  global_dpd_->buf4_scm(&Z, 0.0);
  for(h=0; h< nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&Z, h);
    for(row=0; row < Z.params->rowtot[h]; row++) {
      i = Z.params->roworb[h][row][0];
      j = Z.params->roworb[h][row][1];
      I = X1.params->rowidx[i];
      J = Y1.params->rowidx[j];
      Isym = X1.params->psym[i];
      Jsym = Y1.params->psym[j];
      for(col=0; col < Z.params->coltot[h]; col++) {
	e = Z.params->colorb[h][col][0];
	f = Z.params->colorb[h][col][1];
	E = X1.params->colidx[e];
	F = Y1.params->colidx[f];
	Esym = X1.params->qsym[e];
	Fsym = Y1.params->qsym[f];
	if((Isym^Esym)==irrep_x && (Jsym^Fsym)==irrep_y)
	  Z.matrix[h][row][col] = (X1.matrix[Isym][I][E] * Y1.matrix[Jsym][J][F]);
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&Z, h);
    global_dpd_->buf4_mat_irrep_close(&Z, h);
  }
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "(X*Y+Y*X)(ie,ma)");
  global_dpd_->buf4_scm(&Z, 0.0);
  for(h=0; h< nirreps; h++) {
    global_dpd_->buf4_mat_irrep_init(&Z, h);
    for(row=0; row < Z.params->rowtot[h]; row++) {
      i = Z.params->roworb[h][row][0];
      e = Z.params->roworb[h][row][1];
      I = X1.params->rowidx[i];
      E = X1.params->colidx[e];
      Isym = X1.params->psym[i];
      Esym = X1.params->qsym[e];
      for(col=0; col < Z.params->coltot[h]; col++) {
	m = Z.params->colorb[h][col][0];
	a = Z.params->colorb[h][col][1];
	M = Y1.params->rowidx[m];
	A = Y1.params->colidx[a];
	Msym = Y1.params->psym[m];
	Asym = Y1.params->qsym[a];

	if(((Isym^Esym)==irrep_x) && ((Msym^Asym)==irrep_y))
	  Z.matrix[h][row][col] =
	    (X1.matrix[Isym][I][E] * Y1.matrix[Msym][M][A]) +
	    (Y1.matrix[Isym][I][E] * X1.matrix[Msym][M][A]);
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&Z, h);
    global_dpd_->buf4_mat_irrep_close(&Z, h);
  }
  global_dpd_->buf4_close(&Z);

  global_dpd_->file2_mat_close(&X1);
  global_dpd_->file2_close(&X1);
  global_dpd_->file2_mat_close(&Y1);
  global_dpd_->file2_close(&Y1);

}

}} // namespace psi::ccresponse
