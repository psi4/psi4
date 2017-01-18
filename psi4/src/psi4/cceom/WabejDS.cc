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
    \ingroup CCEOM
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cmath>
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/* This function computes the H-bar doubles-singles block contribution
   of Wabej to a Sigma vector stored at Sigma plus 'i' */

void WabejDS(int i, int C_irr) {
  dpdfile2 CME, Cme;
  dpdbuf4 SIJAB, Sijab, SIjAb;
  dpdbuf4 WABEI, Wabei, WAbEi, WaBeI, WM, WP, Z;
  char CME_lbl[32], Cme_lbl[32], SIJAB_lbl[32], Sijab_lbl[32], SIjAb_lbl[32];
  int Gej, Gab, Gij, Gj, Gi, Ge, nrows, length, E, e, I;
  dpdbuf4 W;

  if (params.eom_ref == 0) { /* RHF */
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WabejDS Z(Ij,Ab)");
    global_dpd_->buf4_scm(&Z, 0);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
/*     dpd_contract244(&CME, &WAbEi, &Z, 1, 0, 0, 1.0, 0.0); */
    global_dpd_->file2_mat_init(&CME);
    global_dpd_->file2_mat_rd(&CME);
    for(Gej=0; Gej < moinfo.nirreps; Gej++) {
      Gab = Gej ^ H_IRR;
      Gij = Gab ^ C_irr;

      global_dpd_->buf4_mat_irrep_init(&Z, Gij);
      global_dpd_->buf4_mat_irrep_shift13(&Z, Gij);

      for(Ge=0; Ge < moinfo.nirreps; Ge++) {
	Gj = Ge ^ Gej;
	Gi = Gj ^ Gij;

	nrows = moinfo.occpi[Gj];
	length = nrows * W.params->coltot[Gab];
	global_dpd_->buf4_mat_irrep_init_block(&W, Gej, nrows);

	for(E=0; E < moinfo.virtpi[Ge]; E++) {
	  e = moinfo.vir_off[Ge] + E;
	  global_dpd_->buf4_mat_irrep_rd_block(&W, Gej, W.row_offset[Gej][e], nrows);

	  for(I=0; I < moinfo.occpi[Gi]; I++) {
	    if(length)
	      C_DAXPY(length, CME.matrix[Gi][I][E], W.matrix[Gej][0], 1,
		      Z.shift.matrix[Gij][Gi][I], 1);
	  }
	}

	global_dpd_->buf4_mat_irrep_close_block(&W, Gej, nrows);
      }

      global_dpd_->buf4_mat_irrep_wrt(&Z, Gij);
      global_dpd_->buf4_mat_irrep_close(&Z, Gij);

    }
    global_dpd_->file2_mat_close(&CME);
    global_dpd_->file2_close(&CME);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_SIjAb, qpsr, 0, 5, SIjAb_lbl, 1);
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    global_dpd_->buf4_axpy(&Z, &SIjAb, 1.0);
    global_dpd_->buf4_close(&SIjAb);
    global_dpd_->buf4_close(&Z);
  }

  else if (params.eom_ref == 1) { /* ROHF */
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(Cme_lbl, "%s %d", "Cme", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIJAB += WABEJ * CIE - WABEI * CJE */
    global_dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WabejDS_P");
    global_dpd_->buf4_init(&WABEI, PSIF_CC_HBAR, H_IRR, 11, 7, 11, 7, 0, "WEIAB");
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->contract244(&CME, &WABEI, &WP, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&CME);
    global_dpd_->buf4_close(&WABEI);
    global_dpd_->buf4_sort(&WP, PSIF_EOM_TMP, qprs, 0, 7, "WabejDS_M");
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 0, 7, 2, 7, 0, SIJAB_lbl);
    global_dpd_->buf4_axpy(&WP, &SIJAB, 1.0);
    global_dpd_->buf4_close(&WP);
    global_dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WabejDS_M");
    global_dpd_->buf4_axpy(&WM, &SIJAB, -1.0);
    global_dpd_->buf4_close(&WM);
    global_dpd_->buf4_close(&SIJAB);

    /* Sijab += Wabej * Cie - Wabei * Cje */
    global_dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WabejDS_P");
    global_dpd_->buf4_init(&Wabei, PSIF_CC_HBAR, H_IRR, 11, 7, 11, 7, 0, "Weiab");
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, Cme_lbl);
    global_dpd_->contract244(&Cme, &Wabei, &WP, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&Cme);
    global_dpd_->buf4_close(&Wabei);
    global_dpd_->buf4_sort(&WP, PSIF_EOM_TMP, qprs, 0, 7, "WabejDS_M");
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 0, 7, 2, 7, 0, Sijab_lbl);
    global_dpd_->buf4_axpy(&WP, &Sijab, 1.0);
    global_dpd_->buf4_close(&WP);
    global_dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WabejDS_M");
    global_dpd_->buf4_axpy(&WM, &Sijab, -1.0);
    global_dpd_->buf4_close(&WM);
    global_dpd_->buf4_close(&Sijab);


    /* SIjAb += WAbEj * CIE - WAbeI * Cje */
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->buf4_init(&WAbEi, PSIF_CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WEiAb");
    global_dpd_->contract244(&CME, &WAbEi, &SIjAb, 1, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&CME);
    global_dpd_->buf4_close(&WAbEi);
    global_dpd_->buf4_init(&WaBeI, PSIF_CC_HBAR, H_IRR, 10, 5, 10, 5, 0, "WeIaB (Ie,Ab)");
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, Cme_lbl);
    global_dpd_->contract424(&WaBeI, &Cme, &SIjAb, 1, 1, 1, 1.0, 1.0);
    global_dpd_->file2_close(&Cme);
    global_dpd_->buf4_close(&WaBeI);
    global_dpd_->buf4_close(&SIjAb);
  }

  else if (params.eom_ref == 2) {
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(Cme_lbl, "%s %d", "Cme", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIJAB += WABEJ * CIE - WABEI * CJE */
    global_dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WabejDS_P");
    global_dpd_->buf4_init(&WABEI, PSIF_CC_HBAR, H_IRR, 21, 7, 21, 7, 0, "WEIAB");
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->contract244(&CME, &WABEI, &WP, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&CME);
    global_dpd_->buf4_close(&WABEI);
    global_dpd_->buf4_sort(&WP, PSIF_EOM_TMP, qprs, 0, 7, "WabejDS_M");
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 0, 7, 2, 7, 0, SIJAB_lbl);
    global_dpd_->buf4_axpy(&WP, &SIJAB, 1.0);
    global_dpd_->buf4_close(&WP);
    global_dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WabejDS_M");
    global_dpd_->buf4_axpy(&WM, &SIJAB, -1.0);
    global_dpd_->buf4_close(&WM);
    global_dpd_->buf4_close(&SIJAB);

    /* Sijab += Wabej * Cie - Wabei * Cje */
    global_dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 10, 17, 10, 17, 0, "WabejDS_PB");
    global_dpd_->buf4_init(&Wabei, PSIF_CC_HBAR, H_IRR, 31, 17, 31, 17, 0, "Weiab");
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, Cme_lbl);
    global_dpd_->contract244(&Cme, &Wabei, &WP, 1, 0, 0, 1.0, 0.0);
    global_dpd_->file2_close(&Cme);
    global_dpd_->buf4_close(&Wabei);
    global_dpd_->buf4_sort(&WP, PSIF_EOM_TMP, qprs, 10, 17, "WabejDS_MB");
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 10, 17, 12, 17, 0, Sijab_lbl);
    global_dpd_->buf4_axpy(&WP, &Sijab, 1.0);
    global_dpd_->buf4_close(&WP);
    global_dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 10, 17, 10, 17, 0, "WabejDS_MB");
    global_dpd_->buf4_axpy(&WM, &Sijab, -1.0);
    global_dpd_->buf4_close(&WM);
    global_dpd_->buf4_close(&Sijab);


    /* SIjAb += WAbEj * CIE - WAbeI * Cje */
    /* start here */
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, SIjAb_lbl);
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->buf4_init(&WAbEi, PSIF_CC_HBAR, H_IRR, 26, 28, 26, 28, 0, "WEiAb");
    global_dpd_->contract244(&CME, &WAbEi, &SIjAb, 1, 0, 0, 1.0, 1.0);
    global_dpd_->file2_close(&CME);
    global_dpd_->buf4_close(&WAbEi);
    global_dpd_->buf4_init(&WaBeI, PSIF_CC_HBAR, H_IRR, 24, 28, 24, 28, 0, "WeIaB (Ie,Ab)");
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, Cme_lbl);
    global_dpd_->contract424(&WaBeI, &Cme, &SIjAb, 1, 1, 1, 1.0, 1.0);
    global_dpd_->file2_close(&Cme);
    global_dpd_->buf4_close(&WaBeI);
    global_dpd_->buf4_close(&SIjAb);
  }

#ifdef EOM_DEBUG
  check_sum("WabejDS",i,C_irr);
#endif
  return;
}


}} // namespace psi::cceom
