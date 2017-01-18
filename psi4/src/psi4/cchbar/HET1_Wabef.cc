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
    \ingroup CCHBAR
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

/* HET1_Wabef(): Builds an ROHF-like [H,e^T1] Wabef for use in CC3
 * response codes.  Eventually this function will go away, but it will
 * be useful for getting the first versions of the CC3 codes up and
 * running.
 *
 * Wabef = <ab||ef> - P(ab) t_n^b <an||ef> + 1/2 P(ab) t_m^a <mn||ef> t_n^b
 *
 * The last two terms may be factored together to give:
 *
 * Wabef = <ab||ef> - P(ab) t_n^b Z_anef
 *
 * where Z_anef = <an||ef> - 1/2 t_m^a <mn||ef>
 *
 * To avoid memory bottlenecks, which could be severe in this case, I
 * compute the Z_anef contributions out-of-core with a major loop over
 * rows ab of the target Wabef.  With an inner loop over irreps of the
 * index n, I grab from disk all Z_anef for a given index a and
 * compute a DGEMV of this with the b-th column of T1.  This algorithm
 * requires only one quantity of size OV^2, and up to three others of
 * size V^2, depending on the spin case.
 *
 * Note that, although this code runs in ROHF mode, it will really
 * only be used in RHF-based codes, so I only produce the ABEF and
 * AbEf spin cases.
 *
 * TDC, 7/04
*/

void HET1_Wabef(void)
{
  dpdbuf4 Bints, Fints, Z, D, W;
  dpdfile2 T1;
  int Gab, nirreps, ab, ba, A, B, Ga, Gb, a, b;
  int Gn, Gan, an, Gbn, bn;
  int ef, fe, E, F;
  int nrows, ncols;
  double *ZEf, *ZfE;

  nirreps = moinfo.nirreps;

  if(params.ref == 1) {

    global_dpd_->buf4_init(&Bints, PSIF_CC_BINTS, 0, 7, 7, 5, 5, 1, "B <ab|cd>");
    global_dpd_->buf4_copy(&Bints, PSIF_CC3_HET1, "CC3 WABEF");
    global_dpd_->buf4_close(&Bints);

    global_dpd_->buf4_init(&Bints, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    global_dpd_->buf4_copy(&Bints, PSIF_CC3_HET1, "CC3 WAbEf");
    global_dpd_->buf4_close(&Bints);

    /* ZANEF = <AN||EF> - 1/2 t_M^A <MN||EF> */
    global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 11, 7, 11, 5, 1, "F <ai|bc>");
    global_dpd_->buf4_copy(&Fints, PSIF_CC_TMP0, "ZANEF (AN,E>F)");
    global_dpd_->buf4_close(&Fints);

    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 7, 11, 7, 0, "ZANEF (AN,E>F)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    global_dpd_->contract244(&T1, &D, &Z, 0, 0, 0, -0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&T1);

    /* WABEF <-- -P(AB) t_N^B ZANEF */
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 7, 7, 7, 7, 0, "CC3 WABEF");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 7, 11, 7, 0, "ZANEF (AN,E>F)");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&T1);
    global_dpd_->file2_mat_rd(&T1);
    for(Gab=0; Gab < nirreps; Gab++) {
      ncols = W.params->coltot[Gab];
      W.matrix[Gab] = global_dpd_->dpd_block_matrix(1, ncols);
      for(ab=0; ab < W.params->rowtot[Gab]; ab++) {
	A = W.params->roworb[Gab][ab][0];
	B = W.params->roworb[Gab][ab][1];
	Ga = W.params->psym[A];
	Gb = W.params->qsym[B];
	a = A - W.params->poff[Ga];
	b = B - W.params->qoff[Gb];

	global_dpd_->buf4_mat_irrep_rd_block(&W, Gab, ab, 1);

	for(Gn=0; Gn < nirreps; Gn++) {
	  nrows = Z.params->qpi[Gn];

	  if(Gn == Gb) {
	    Gan = Ga ^ Gn;
	    an = Z.row_offset[Gan][A];
	    Z.matrix[Gan] = global_dpd_->dpd_block_matrix(nrows, ncols);
	    global_dpd_->buf4_mat_irrep_rd_block(&Z, Gan, an, nrows);

	    if(nrows && ncols)
	      C_DGEMV('t', nrows, ncols, -1.0, Z.matrix[Gan][0], ncols,
		      &(T1.matrix[Gn][0][b]), T1.params->coltot[Gn], 1.0, W.matrix[Gab][0], 1);

	    global_dpd_->free_dpd_block(Z.matrix[Gan], nrows, ncols);
	  }

	  if(Gn == Ga) {
	    Gbn = Gb ^ Gn;
	    bn = Z.row_offset[Gbn][B];
	    Z.matrix[Gbn] = global_dpd_->dpd_block_matrix(nrows, ncols);
	    global_dpd_->buf4_mat_irrep_rd_block(&Z, Gbn, bn, nrows);

	    if(nrows && ncols)
	      C_DGEMV('t', nrows, ncols, 1.0, Z.matrix[Gbn][0], ncols,
		      &(T1.matrix[Gn][0][a]), T1.params->coltot[Gn], 1.0, W.matrix[Gab][0], 1);

	    global_dpd_->free_dpd_block(Z.matrix[Gbn], nrows, ncols);
	  }
	}

	global_dpd_->buf4_mat_irrep_wrt_block(&W, Gab, ab, 1);
      }
      global_dpd_->free_dpd_block(W.matrix[Gab], 1, ncols);
    }
    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /* ZAnEf = <An|Ef> - 1/2 t_MA <Mn|Ef> */
    global_dpd_->buf4_init(&Fints, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    global_dpd_->buf4_copy(&Fints, PSIF_CC_TMP0, "ZAnEf");
    global_dpd_->buf4_close(&Fints);

    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "ZAnEf");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract244(&T1, &D, &Z, 0, 0, 0, -0.5, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&T1);

    /* WAbEf <-- P(Ab) t_n^b ZAnEf */
    global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 5, 5, 5, 5, 0, "CC3 WAbEf");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 5, 11, 5, 0, "ZAnEf");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&T1);
    global_dpd_->file2_mat_rd(&T1);
    for(Gab=0; Gab < nirreps; Gab++) {
      ncols = W.params->coltot[Gab];
      W.matrix[Gab] = global_dpd_->dpd_block_matrix(1, ncols);
      ZEf = init_array(ncols);
      ZfE = init_array(ncols);
      for(ab=0; ab < W.params->rowtot[Gab]; ab++) {
	A = W.params->roworb[Gab][ab][0];
	B = W.params->roworb[Gab][ab][1];
	Ga = W.params->psym[A];
	Gb = W.params->qsym[B];
	a = A - W.params->poff[Ga];
	b = B - W.params->qoff[Gb];
	ba = W.params->rowidx[B][A];
	zero_arr(ZEf, ncols);  zero_arr(ZfE, ncols);
	for(Gn=0; Gn < nirreps; Gn++) {
	  nrows = Z.params->qpi[Gn];
	  if(Gn == Gb) {
	    Gan = Ga ^ Gn;
	    an = Z.row_offset[Gan][A];
	    Z.matrix[Gan] = global_dpd_->dpd_block_matrix(nrows, ncols);
	    global_dpd_->buf4_mat_irrep_rd_block(&Z, Gan, an, nrows);
	    if(nrows && ncols)
	      C_DGEMV('t', nrows, ncols, -1.0, Z.matrix[Gan][0], ncols,
		      &(T1.matrix[Gn][0][b]), T1.params->coltot[Gn], 1.0, ZEf, 1);
	    global_dpd_->free_dpd_block(Z.matrix[Gan], nrows, ncols);
	  }
	}
	/* Sort Ef->fE */
	for(ef=0; ef < W.params->coltot[Gab]; ef++) {
	  E = W.params->colorb[Gab][ef][0];
	  F = W.params->colorb[Gab][ef][1];
	  fe = W.params->colidx[F][E];
	  ZfE[fe] = ZEf[ef];
	}
	global_dpd_->buf4_mat_irrep_rd_block(&W, Gab, ab, 1);
	C_DAXPY(ncols, 1.0, ZEf, 1, W.matrix[Gab][0], 1);
	global_dpd_->buf4_mat_irrep_wrt_block(&W, Gab, ab, 1);
	global_dpd_->buf4_mat_irrep_rd_block(&W, Gab, ba, 1);
	C_DAXPY(ncols, 1.0, ZfE, 1, W.matrix[Gab][0], 1);
	global_dpd_->buf4_mat_irrep_wrt_block(&W, Gab, ba, 1);
      }
      free(ZEf);
      free(ZfE);
      global_dpd_->free_dpd_block(W.matrix[Gab], 1, ncols);
    }

    global_dpd_->file2_mat_close(&T1);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);
  }
}

}} // namespace psi::cchbar
