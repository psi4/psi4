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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

void cc3_l3l1(void)
{
  dpdfile2 L1, D1, L1new;
  dpdbuf4 Z, W;
  int nirreps, Gde, Gg, Gi, Ga;
  int de, ig, ag;
  int nrows, ncols, nlinks;

  nirreps = moinfo.nirreps;

  global_dpd_->buf4_init(&Z, PSIF_CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZIGDE");
  global_dpd_->buf4_sort(&Z, PSIF_CC3_MISC, rspq, 5, 10, "CC3 ZIGDE (DE,IG)");
  global_dpd_->buf4_close(&Z);

  global_dpd_->buf4_init(&Z, PSIF_CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZIgDe");
  global_dpd_->buf4_sort(&Z, PSIF_CC3_MISC, rspq, 5, 10, "CC3 ZIgDe (De,Ig)");
  global_dpd_->buf4_close(&Z);

  global_dpd_->file2_init(&L1, PSIF_CC3_MISC, 0, 0, 1, "CC3 LIA");
  global_dpd_->file2_mat_init(&L1);

  global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 5, 5, 7, 7, 0, "CC3 WABEF");
  global_dpd_->buf4_init(&Z, PSIF_CC3_MISC, 0, 5, 10, 5, 10, 0, "CC3 ZIGDE (DE,IG)");
  for(Gde=0; Gde < nirreps; Gde++) {
    if(Z.params->coltot[Gde] && W.params->coltot[Gde]) {
      Z.matrix[Gde] = global_dpd_->dpd_block_matrix(1, Z.params->coltot[Gde]);
      W.matrix[Gde] = global_dpd_->dpd_block_matrix(1, W.params->coltot[Gde]);
      for(de=0; de < Z.params->rowtot[Gde]; de++) {
	global_dpd_->buf4_mat_irrep_rd_block(&W, Gde, de, 1);
	global_dpd_->buf4_mat_irrep_rd_block(&Z, Gde, de, 1);

	for(Gg=0; Gg < nirreps; Gg++) {
	  Ga = Gi = Gg ^ Gde; /* totally symmetric */
	  nrows = L1.params->rowtot[Gi];
	  ncols = L1.params->coltot[Gi];
	  nlinks = Z.params->spi[Gg];

	  ig = Z.col_offset[Gde][Gi];
	  ag = W.col_offset[Gde][Ga];

	  if(nrows && ncols && nlinks)
	    C_DGEMM('n', 't', nrows, ncols, nlinks, 0.5, &(Z.matrix[Gde][0][ig]), nlinks,
		    &(W.matrix[Gde][0][ag]), nlinks, 1.0, L1.matrix[Gi][0], ncols);
	}
      }
      global_dpd_->free_dpd_block(Z.matrix[Gde], 1, Z.params->coltot[Gde]);
      global_dpd_->free_dpd_block(W.matrix[Gde], 1, W.params->coltot[Gde]);
    }
  }
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 5, 5, 5, 5, 0, "CC3 WAbEf");
  global_dpd_->buf4_init(&Z, PSIF_CC3_MISC, 0, 5, 10, 5, 10, 0, "CC3 ZIgDe (De,Ig)");
  for(Gde=0; Gde < nirreps; Gde++) {
    if(Z.params->coltot[Gde] && W.params->coltot[Gde]) {
      Z.matrix[Gde] = global_dpd_->dpd_block_matrix(1, Z.params->coltot[Gde]);
      W.matrix[Gde] = global_dpd_->dpd_block_matrix(1, W.params->coltot[Gde]);
      for(de=0; de < Z.params->rowtot[Gde]; de++) {
	global_dpd_->buf4_mat_irrep_rd_block(&W, Gde, de, 1);
	global_dpd_->buf4_mat_irrep_rd_block(&Z, Gde, de, 1);

	for(Gg=0; Gg < nirreps; Gg++) {
	  Ga = Gi = Gg ^ Gde; /* totally symmetric */
	  nrows = L1.params->rowtot[Gi];
	  ncols = L1.params->coltot[Gi];
	  nlinks = Z.params->spi[Gg];

	  ig = Z.col_offset[Gde][Gi];
	  ag = W.col_offset[Gde][Ga];

	  if(nrows && ncols && nlinks)
	    C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, &(Z.matrix[Gde][0][ig]), nlinks,
		    &(W.matrix[Gde][0][ag]), nlinks, 1.0, L1.matrix[Gi][0], ncols);
	}
      }
      global_dpd_->free_dpd_block(Z.matrix[Gde], 1, Z.params->coltot[Gde]);
      global_dpd_->free_dpd_block(W.matrix[Gde], 1, W.params->coltot[Gde]);
    }
  }
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  global_dpd_->file2_mat_wrt(&L1);
  global_dpd_->file2_mat_close(&L1);

  /* Wmbej --> L1 */

  global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMBEJ (ME,JB)");
  global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, psrq, 10, 10, "CC3 WMBEJ (MB,JE)");
  global_dpd_->buf4_close(&W);

  global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbEj (ME,jb)");
  global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, psrq, 10, 10, "CC3 WMbEj (Mb,jE)");
  global_dpd_->buf4_close(&W);

  global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbeJ (Me,Jb)");
  global_dpd_->buf4_sort(&W, PSIF_CC3_HET1, psrq, 10, 10, "CC3 WMbeJ (Mb,Je)");
  global_dpd_->buf4_close(&W);

  global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMBEJ (MB,JE)");
  global_dpd_->buf4_init(&Z, PSIF_CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZDMAE (MD,AE)");
  global_dpd_->contract442(&W, &Z, &L1, 0, 2, -1, 1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbEj (Mb,jE)");
  global_dpd_->buf4_init(&Z, PSIF_CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZDmAe (mD,Ae)");
  global_dpd_->contract442(&W, &Z, &L1, 0, 2, -1, 1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 10, 10, 10, 0, "CC3 WMbeJ (Mb,Je)");
  global_dpd_->buf4_init(&Z, PSIF_CC3_MISC, 0, 10, 5, 10, 5, 0, "CC3 ZdMAe (Md,Ae)");
  global_dpd_->contract442(&W, &Z, &L1, 0, 2, -1, 1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 11, 10, 11, 0, "CC3 WMBEJ (MB,EJ)");
  global_dpd_->buf4_init(&Z, PSIF_CC3_MISC, 0, 0, 10, 0, 10, 0, "CC3 ZIMLE");
  global_dpd_->contract442(&Z, &W, &L1, 0, 2, 1, 1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 11, 10, 11, 0, "CC3 WMbEj (Mb,Ej)");
  global_dpd_->buf4_init(&Z, PSIF_CC3_MISC, 0, 0, 10, 0, 10, 0, "CC3 ZImLe");
  global_dpd_->contract442(&Z, &W, &L1, 0, 2, 1, 1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 10, 11, 10, 11, 0, "CC3 WmBEj (mB,Ej)");
  global_dpd_->buf4_init(&Z, PSIF_CC3_MISC, 0, 0, 10, 0, 10, 0, "CC3 ZImlE");
  global_dpd_->contract442(&Z, &W, &L1, 0, 2, 1, 1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /* Wmnij -> L1 */

  global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 0, 2, 2, 0, "CC3 WMNIJ (M>N,I>J)");
  global_dpd_->buf4_init(&Z, PSIF_CC3_MISC, 0, 0, 11, 0, 11, 0, "CC3 ZLMAO");
  global_dpd_->contract442(&W, &Z, &L1, 0, 2, -0.5, 1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  global_dpd_->buf4_init(&W, PSIF_CC3_HET1, 0, 0, 0, 0, 0, 0, "CC3 WMnIj (Mn,Ij)");
  global_dpd_->buf4_init(&Z, PSIF_CC3_MISC, 0, 0, 11, 0, 11, 0, "CC3 ZLmAo");
  global_dpd_->contract442(&W, &Z, &L1, 0, 2, -1, 1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  global_dpd_->file2_init(&D1, PSIF_CC_DENOM, 0, 0, 1, "dIA");
  global_dpd_->file2_dirprd(&D1, &L1);
  global_dpd_->file2_close(&D1);
  global_dpd_->file2_init(&L1new, PSIF_CC_LAMBDA, 0, 0, 1, "New LIA");
  global_dpd_->file2_axpy(&L1, &L1new, 1, 0);
  global_dpd_->file2_copy(&L1new, PSIF_CC_LAMBDA, "New Lia");
  global_dpd_->file2_close(&L1new);
  global_dpd_->file2_close(&L1);

}

}} // namespace psi::cclambda
