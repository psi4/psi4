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
#include <cstdlib>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/* This function computes the H-bar doubles-singles block contribution
   of Wbmfe to a Sigma vector stored at Sigma plus 'i' */

void WbmfeDS(int i, int C_irr) {
  dpdfile2 CME, Cme, XBF, Xbf;
  dpdbuf4 SIJAB, Sijab, SIjAb;
  dpdbuf4 WAMEF, Wamef, WAmEf, WaMeF, WM, WP, W, Z;
  dpdbuf4 TIJAB, TIjAb, Tijab;
  char CME_lbl[32], Cme_lbl[32], SIJAB_lbl[32], Sijab_lbl[32], SIjAb_lbl[32];
  int Gbm, Gfe, bm, b, m, Gb, Gm, Ge, Gf, B, M, f, e, fe, ef, nrows, ncols;
  dpdfile2 C;
  double *X;

  if (params.eom_ref == 0) { /* RHF */
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* Form Xbf intermediates */
/*     dpd_file2_init(&XBF, EOM_TMP, C_irr, 1, 1, "XBF"); */
/*     dpd_file2_scm(&XBF, 0.0); */
/*     dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl); */
/*     dpd_buf4_init(&W, CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)"); */
/*     dpd_dot24(&CME, &W, &XBF, 0, 0, 1.0, 1.0); */
/*     dpd_buf4_close(&W); */
/*     dpd_file2_close(&CME); */

    /* OOC code below added 7/27/05, -TDC */
    /* X(b,f) = [ 2 Wbmfe - Wbmef ] * C(m,e) */
    global_dpd_->file2_init(&XBF, PSIF_EOM_TMP, C_irr, 1, 1, "XBF");
    global_dpd_->file2_scm(&XBF, 0.0);
    global_dpd_->file2_mat_init(&XBF);
    global_dpd_->file2_mat_rd(&XBF);
    global_dpd_->file2_init(&C, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->file2_mat_init(&C);
    global_dpd_->file2_mat_rd(&C);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WAmEf");
    for(Gbm=0; Gbm < moinfo.nirreps; Gbm++) {
      Gfe = Gbm ^ H_IRR;
      global_dpd_->buf4_mat_irrep_row_init(&W, Gbm);
      X = init_array(W.params->coltot[Gfe]);
      for(bm=0; bm < W.params->rowtot[Gbm]; bm++) {
	global_dpd_->buf4_mat_irrep_row_rd(&W, Gbm, bm);

	b = W.params->roworb[Gbm][bm][0];
	m = W.params->roworb[Gbm][bm][1];
	Gb = W.params->psym[b];
	Gm = Gbm ^ Gb;
	Ge = Gm ^ C_irr;
	Gf = Ge ^ Gfe;
	B = b - moinfo.vir_off[Gb];
	M = m - moinfo.occ_off[Gm];

	zero_arr(X, W.params->coltot[Gfe]);

	for(fe=0; fe < W.params->coltot[Gfe]; fe++) {
	  f = W.params->colorb[Gfe][fe][0];
	  e = W.params->colorb[Gfe][fe][1];
	  ef = W.params->colidx[e][f];
	  X[fe] = 2.0 * W.matrix[Gbm][0][fe] - W.matrix[Gbm][0][ef];
	}

	nrows = moinfo.virtpi[Gf];
	ncols = moinfo.virtpi[Ge];

	if(nrows && ncols)
	  C_DGEMV('n',nrows,ncols,1,&X[W.col_offset[Gfe][Gf]],ncols,
		  C.matrix[Gm][M],1,1,XBF.matrix[Gb][B],1);

      }
      free(X);
      global_dpd_->buf4_mat_irrep_row_close(&W, Gbm);
    }
    global_dpd_->buf4_close(&W);
    global_dpd_->file2_close(&C);
    global_dpd_->file2_mat_wrt(&XBF);
    global_dpd_->file2_mat_close(&XBF);

    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WbmfeDS Z(Ij,Ab)");
    global_dpd_->buf4_init(&TIjAb, PSIF_CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->contract424(&TIjAb, &XBF, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&TIjAb);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_EOM_SIjAb, qpsr, 0, 5, SIjAb_lbl, 1);
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    global_dpd_->buf4_axpy(&Z, &SIjAb, 1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&SIjAb);
    global_dpd_->file2_close(&XBF);
  }

  else if (params.eom_ref == 1) { /* ROHF */
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(Cme_lbl, "%s %d", "Cme", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* Form Xbf intermediates */
    /* XBF = CME * WBMFE + Cme * WBmFe */
    global_dpd_->file2_init(&XBF, PSIF_EOM_TMP, C_irr, 1, 1, "XBF");
    global_dpd_->file2_scm(&XBF, 0.0);
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->buf4_init(&WAMEF, PSIF_CC_HBAR, H_IRR, 11, 5, 11, 7, 0, "WAMEF");
    global_dpd_->dot24(&CME, &WAMEF, &XBF, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&WAMEF);
    global_dpd_->file2_close(&CME);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, Cme_lbl);
    global_dpd_->buf4_init(&WAmEf, PSIF_CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->dot24(&Cme, &WAmEf, &XBF, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&WAmEf);
    global_dpd_->file2_close(&Cme);
    global_dpd_->file2_close(&XBF);

    /* Xbf = Cme * Wbmfe + CME * WbMfE */
    global_dpd_->file2_init(&Xbf, PSIF_EOM_TMP, C_irr, 1, 1, "Xbf");
    global_dpd_->file2_scm(&Xbf, 0.0);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, Cme_lbl);
    global_dpd_->buf4_init(&Wamef, PSIF_CC_HBAR, H_IRR, 11, 5, 11, 7, 0, "Wamef");
    global_dpd_->dot24(&Cme, &Wamef, &Xbf, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&Wamef);
    global_dpd_->file2_close(&Cme);
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->buf4_init(&WaMeF, PSIF_CC_HBAR, H_IRR, 11, 5, 11, 5, 0, "WaMeF");
    global_dpd_->dot24(&CME, &WaMeF, &Xbf, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&WaMeF);
    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&Xbf);

    /* SIJAB += XBF * TIJAF - XAF * TIJBF */
    global_dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WbmfeDS_P");
    global_dpd_->file2_init(&XBF, PSIF_EOM_TMP, C_irr, 1, 1, "XBF");
    global_dpd_->buf4_init(&TIJAB, PSIF_CC_TAMPS, H_IRR, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->contract424(&TIJAB, &XBF, &WP, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&TIJAB);
    global_dpd_->file2_close(&XBF);
    global_dpd_->buf4_sort(&WP, PSIF_EOM_TMP, pqsr, 2, 5, "WbmfeDS_M");
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 2, 5, 2, 7, 0, SIJAB_lbl);
    global_dpd_->buf4_axpy(&WP, &SIJAB, 1.0);
    global_dpd_->buf4_close(&WP);
    global_dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WbmfeDS_M");
    global_dpd_->buf4_axpy(&WM, &SIJAB, -1.0);
    global_dpd_->buf4_close(&WM);
    global_dpd_->buf4_close(&SIJAB);

    /* Sijab += Xbf * Tijaf - Xaf * Tijbf */
    global_dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WbmfeDS_P");
    global_dpd_->file2_init(&Xbf, PSIF_EOM_TMP, C_irr, 1, 1, "Xbf");
    global_dpd_->buf4_init(&Tijab, PSIF_CC_TAMPS, H_IRR, 2, 5, 2, 7, 0, "tijab");
    global_dpd_->contract424(&Tijab, &Xbf, &WP, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Tijab);
    global_dpd_->file2_close(&Xbf);
    global_dpd_->buf4_sort(&WP, PSIF_EOM_TMP, pqsr, 2, 5, "WbmfeDS_M");
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 2, 5, 2, 7, 0, Sijab_lbl);
    global_dpd_->buf4_axpy(&WP, &Sijab, 1.0);
    global_dpd_->buf4_close(&WP);
    global_dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WbmfeDS_M");
    global_dpd_->buf4_axpy(&WM, &Sijab, -1.0);
    global_dpd_->buf4_close(&WM);
    global_dpd_->buf4_close(&Sijab);

    /* SIjAb += Xbf * tIjAf + XAF * TIjbF */
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    global_dpd_->buf4_init(&TIjAb, PSIF_CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->file2_init(&Xbf, PSIF_EOM_TMP, C_irr, 1, 1, "Xbf");
    global_dpd_->contract424(&TIjAb, &Xbf, &SIjAb, 3, 1, 0, 1.0, 1.0);
    global_dpd_->file2_close(&Xbf);
    global_dpd_->file2_init(&XBF, PSIF_EOM_TMP, C_irr, 1, 1, "XBF");
    global_dpd_->contract244(&XBF, &TIjAb, &SIjAb, 1, 2, 1, 1.0, 1.0);
    global_dpd_->file2_close(&XBF);
    global_dpd_->buf4_close(&TIjAb);
    global_dpd_->buf4_close(&SIjAb);
  }

  else { /* UHF */
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(Cme_lbl, "%s %d", "Cme", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* Form Xbf intermediates */
    /* XBF = CME * WBMFE + Cme * WBmFe */
    global_dpd_->file2_init(&XBF, PSIF_EOM_TMP, C_irr, 1, 1, "XBF");
    global_dpd_->file2_scm(&XBF, 0.0);
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->buf4_init(&WAMEF, PSIF_CC_HBAR, H_IRR, 21, 5, 21, 7, 0, "WAMEF");
    global_dpd_->dot24(&CME, &WAMEF, &XBF, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&WAMEF);
    global_dpd_->file2_close(&CME);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, Cme_lbl);
    global_dpd_->buf4_init(&WAmEf, PSIF_CC_HBAR, H_IRR, 26, 28, 26, 28, 0, "WAmEf");
    global_dpd_->dot24(&Cme, &WAmEf, &XBF, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&WAmEf);
    global_dpd_->file2_close(&Cme);
/*
outfile->Printf("XBF self dot %15.10lf\n", dpd_file2_dot_self(&XBF));
*/
    global_dpd_->file2_close(&XBF);

    /* Xbf = Cme * Wbmfe + CME * WbMfE */
    global_dpd_->file2_init(&Xbf, PSIF_EOM_TMP, C_irr, 3, 3, "Xbf");
    global_dpd_->file2_scm(&Xbf, 0.0);
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, Cme_lbl);
    global_dpd_->buf4_init(&Wamef, PSIF_CC_HBAR, H_IRR, 31, 15, 31, 17, 0, "Wamef");
    global_dpd_->dot24(&Cme, &Wamef, &Xbf, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&Wamef);
    global_dpd_->file2_close(&Cme);
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->buf4_init(&WaMeF, PSIF_CC_HBAR, H_IRR, 25, 29, 25, 29, 0, "WaMeF");
    global_dpd_->dot24(&CME, &WaMeF, &Xbf, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&WaMeF);
    global_dpd_->file2_close(&CME);
/*
outfile->Printf("Xbf self dot %15.10lf\n", dpd_file2_dot_self(&Xbf));
*/
    global_dpd_->file2_close(&Xbf);

    /* SIJAB += XBF * TIJAF - XAF * TIJBF */
    global_dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WbmfeDS_P");
    global_dpd_->file2_init(&XBF, PSIF_EOM_TMP, C_irr, 1, 1, "XBF");
    global_dpd_->buf4_init(&TIJAB, PSIF_CC_TAMPS, H_IRR, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->contract424(&TIJAB, &XBF, &WP, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&TIJAB);
    global_dpd_->file2_close(&XBF);
    global_dpd_->buf4_sort(&WP, PSIF_EOM_TMP, pqsr, 2, 5, "WbmfeDS_M");
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 2, 5, 2, 7, 0, SIJAB_lbl);
    global_dpd_->buf4_axpy(&WP, &SIJAB, 1.0);
    global_dpd_->buf4_close(&WP);
    global_dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WbmfeDS_M");
    global_dpd_->buf4_axpy(&WM, &SIJAB, -1.0);
    global_dpd_->buf4_close(&WM);
    global_dpd_->buf4_close(&SIJAB);

    /* Sijab += Xbf * Tijaf - Xaf * Tijbf */
    global_dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 12, 15, 12, 15, 0, "WbmfeDS_PB");
    global_dpd_->file2_init(&Xbf, PSIF_EOM_TMP, C_irr, 3, 3, "Xbf");
    global_dpd_->buf4_init(&Tijab, PSIF_CC_TAMPS, H_IRR, 12, 15, 12, 17, 0, "tijab");
    global_dpd_->contract424(&Tijab, &Xbf, &WP, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Tijab);
    global_dpd_->file2_close(&Xbf);
    global_dpd_->buf4_sort(&WP, PSIF_EOM_TMP, pqsr, 12, 15, "WbmfeDS_MB");
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 12, 15, 12, 17, 0, Sijab_lbl);
    global_dpd_->buf4_axpy(&WP, &Sijab, 1.0);
    global_dpd_->buf4_close(&WP);
    global_dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 12, 15, 12, 15, 0, "WbmfeDS_MB");
    global_dpd_->buf4_axpy(&WM, &Sijab, -1.0);
    global_dpd_->buf4_close(&WM);
    global_dpd_->buf4_close(&Sijab);

    /* SIjAb += Xbf * tIjAf + XAF * TIjbF */
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, SIjAb_lbl);
    global_dpd_->buf4_init(&TIjAb, PSIF_CC_TAMPS, H_IRR, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->file2_init(&Xbf, PSIF_EOM_TMP, C_irr, 3, 3, "Xbf");
    global_dpd_->contract424(&TIjAb, &Xbf, &SIjAb, 3, 1, 0, 1.0, 1.0);
    global_dpd_->file2_close(&Xbf);
    global_dpd_->file2_init(&XBF, PSIF_EOM_TMP, C_irr, 1, 1, "XBF");
    global_dpd_->contract244(&XBF, &TIjAb, &SIjAb, 1, 2, 1, 1.0, 1.0);
    global_dpd_->file2_close(&XBF);
    global_dpd_->buf4_close(&TIjAb);
    global_dpd_->buf4_close(&SIjAb);
  }

#ifdef EOM_DEBUG
  check_sum("WbmfeDS",i,C_irr);
#endif
  return;
}

}} // namespace psi::cceom
