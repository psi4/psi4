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
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/* This function computes the H-bar doubles-singles block contribution
   of Wnmje to a Sigma vector stored at Sigma plus 'i' */

void WnmjeDS(int i, int C_irr) {
  dpdfile2 CME, Cme, XNJ, Xnj;
  dpdbuf4 SIJAB, Sijab, SIjAb, Z;
  dpdbuf4 WMNIE, Wmnie, WMnIe, WmNiE, WM, WP, W;
  dpdbuf4 TIJAB, TIjAb, Tijab;
  char CME_lbl[32], Cme_lbl[32], SIJAB_lbl[32], Sijab_lbl[32], SIjAb_lbl[32];
  double tval;

  if (params.eom_ref == 0) { /* RHF */
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* Form XNJ intermediates */
    dpd_->file2_init(&XNJ, PSIF_EOM_TMP, C_irr, 0, 0, "XNJ");
    dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
    dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_->dot23(&CME, &W, &XNJ, 0, 0, 1.0, 0.0);
    dpd_->file2_close(&CME);
    dpd_->buf4_close(&W);

    dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WnmjeDS Z(Ij,Ab)");
    dpd_->buf4_init(&TIjAb, PSIF_CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tIjAb");
    dpd_->contract244(&XNJ, &TIjAb, &Z, 0, 0, 0, 1.0, 0.0);
    dpd_->buf4_close(&TIjAb);
    dpd_->file2_close(&XNJ);
    dpd_->buf4_sort(&Z, PSIF_EOM_TMP, qpsr, 0, 5, "WnmjeDS Z(jI,bA)"); 

    dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_->buf4_axpy(&Z, &SIjAb, -1.0);
    dpd_->buf4_close(&Z);
    dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WnmjeDS Z(jI,bA)");
    dpd_->buf4_axpy(&Z, &SIjAb, -1.0);
    dpd_->buf4_close(&Z);
    dpd_->buf4_close(&SIjAb);
  }

  else if (params.eom_ref == 1) { /* ROHF */
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(Cme_lbl, "%s %d", "Cme", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* Form XNJ intermediates */
    /* XNJ = CME * WNMJE + Cme * WNmJe */
    dpd_->file2_init(&XNJ, PSIF_EOM_TMP, C_irr, 0, 0, "XNJ");
    dpd_->file2_scm(&XNJ, 0.0);
    dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_->buf4_init(&WMNIE, PSIF_CC_HBAR, H_IRR, 0, 11, 2, 11, 0, "WMNIE (M>N,EI)");
    dpd_->dot23(&CME, &WMNIE, &XNJ, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&WMNIE);
    dpd_->file2_close(&CME);
    dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, Cme_lbl);
    dpd_->buf4_init(&WMnIe, PSIF_CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
    dpd_->dot23(&Cme, &WMnIe, &XNJ, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&WMnIe);
    dpd_->file2_close(&Cme);
    dpd_->file2_close(&XNJ);


    /* Xnj = Cme * Wnmje + CME * WnMjE */
    dpd_->file2_init(&Xnj, PSIF_EOM_TMP, C_irr, 0, 0, "Xnj");
    dpd_->file2_scm(&Xnj, 0.0);
    dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, Cme_lbl);
    dpd_->buf4_init(&Wmnie, PSIF_CC_HBAR, H_IRR, 0, 11, 2, 11, 0, "Wmnie (m>n,ei)");
    dpd_->dot23(&Cme, &Wmnie, &Xnj, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&Wmnie);
    dpd_->file2_close(&Cme);
    dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_->buf4_init(&WmNiE, PSIF_CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "WmNiE (mN,Ei)");
    dpd_->dot23(&CME, &WmNiE, &Xnj, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&WmNiE);
    dpd_->file2_close(&CME);
    dpd_->file2_close(&Xnj);

    /* SIJAB -= XNJ * TINAB + XNI * TJNAB */
    dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WnmjeDS_M");
    dpd_->buf4_init(&TIJAB, PSIF_CC_TAMPS, H_IRR, 0, 7, 2, 7, 0, "tIJAB");
    dpd_->file2_init(&XNJ, PSIF_EOM_TMP, C_irr, 0, 0, "XNJ");
    dpd_->contract424(&TIJAB, &XNJ, &WM, 1, 0, 1, 1.0, 0.0);
    dpd_->file2_close(&XNJ);
    dpd_->buf4_close(&TIJAB);
    dpd_->buf4_sort(&WM, PSIF_EOM_TMP, qprs, 0, 7, "WnmjeDS_P");
    dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 0, 7, 2, 7, 0, SIJAB_lbl);
    dpd_->buf4_axpy(&WM, &SIJAB, -1.0);
    dpd_->buf4_close(&WM);
    dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WnmjeDS_P");
    dpd_->buf4_axpy(&WP, &SIJAB, 1.0);
    dpd_->buf4_close(&WP);
    dpd_->buf4_close(&SIJAB);

    /* Sijab -= Xnj * Tinab + Xni * Tjnab */
    dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WnmjeDS_M");
    dpd_->buf4_init(&Tijab, PSIF_CC_TAMPS, H_IRR, 0, 7, 2, 7, 0, "tijab");
    dpd_->file2_init(&Xnj, PSIF_EOM_TMP, C_irr, 0, 0, "Xnj");
    dpd_->contract424(&Tijab, &Xnj, &WM, 1, 0, 1, 1.0, 0.0);
    dpd_->file2_close(&Xnj);
    dpd_->buf4_close(&Tijab);
    dpd_->buf4_sort(&WM, PSIF_EOM_TMP, qprs, 0, 7, "WnmjeDS_P");
    dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 0, 7, 2, 7, 0, Sijab_lbl);
    dpd_->buf4_axpy(&WM, &Sijab, -1.0);
    dpd_->buf4_close(&WM);
    dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WnmjeDS_P");
    dpd_->buf4_axpy(&WP, &Sijab, 1.0);
    dpd_->buf4_close(&WP);
    dpd_->buf4_close(&Sijab);

    /* SIjAb -= Xnj * tInAb + XNI * TjNAb */
    dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_->buf4_init(&TIjAb, PSIF_CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tIjAb");
    dpd_->file2_init(&Xnj, PSIF_EOM_TMP, C_irr, 0, 0, "Xnj");
    dpd_->contract424(&TIjAb, &Xnj, &SIjAb, 1, 0, 1, -1.0, 1.0);
    dpd_->file2_close(&Xnj);
    dpd_->file2_init(&XNJ, PSIF_EOM_TMP, C_irr, 0, 0, "XNJ");
    dpd_->contract244(&XNJ, &TIjAb, &SIjAb, 0, 0, 0, -1.0, 1.0);
    dpd_->file2_close(&XNJ);
    dpd_->buf4_close(&TIjAb);
    dpd_->buf4_close(&SIjAb);
  }

  else if (params.eom_ref == 2) {
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(Cme_lbl, "%s %d", "Cme", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* Form XNJ intermediates */
    /* XNJ = CME * WNMJE + Cme * WNmJe */
    dpd_->file2_init(&XNJ, PSIF_EOM_TMP, C_irr, 0, 0, "XNJ");
    dpd_->file2_scm(&XNJ, 0.0);
    dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_->buf4_init(&WMNIE, PSIF_CC_HBAR, H_IRR, 0, 21, 2, 21, 0, "WMNIE (M>N,EI)");
    dpd_->dot23(&CME, &WMNIE, &XNJ, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&WMNIE);
    dpd_->file2_close(&CME);
    dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, Cme_lbl);
    dpd_->buf4_init(&WMnIe, PSIF_CC_HBAR, H_IRR, 22, 25, 22, 25, 0, "WMnIe (Mn,eI)");
    dpd_->dot23(&Cme, &WMnIe, &XNJ, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&WMnIe);
    dpd_->file2_close(&Cme);
/*
    tval = dpd_file2_dot_self(&XNJ);
    fprintf(outfile,"XNJ self dot %15.10lf\n",tval);
*/
    dpd_->file2_close(&XNJ);


    /* Xnj = Cme * Wnmje + CME * WnMjE */
    dpd_->file2_init(&Xnj, PSIF_EOM_TMP, C_irr, 2, 2, "Xnj");
    dpd_->file2_scm(&Xnj, 0.0);
    dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, Cme_lbl);
    dpd_->buf4_init(&Wmnie, PSIF_CC_HBAR, H_IRR, 10, 31, 12, 31, 0, "Wmnie (m>n,ei)");
    dpd_->dot23(&Cme, &Wmnie, &Xnj, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&Wmnie);
    dpd_->file2_close(&Cme);
    dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_->buf4_init(&WmNiE, PSIF_CC_HBAR, H_IRR, 23, 26, 23, 26, 0, "WmNiE (mN,Ei)");
    dpd_->dot23(&CME, &WmNiE, &Xnj, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&WmNiE);
    dpd_->file2_close(&CME);

/*
    tval = dpd_file2_dot_self(&Xnj);
    fprintf(outfile,"Xnj self dot %15.10lf\n",tval);
*/
    dpd_->file2_close(&Xnj);

    /* SIJAB -= XNJ * TINAB + XNI * TJNAB */
    dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WnmjeDS_M");
    dpd_->buf4_init(&TIJAB, PSIF_CC_TAMPS, H_IRR, 0, 7, 2, 7, 0, "tIJAB");
    dpd_->file2_init(&XNJ, PSIF_EOM_TMP, C_irr, 0, 0, "XNJ");
    dpd_->contract424(&TIJAB, &XNJ, &WM, 1, 0, 1, 1.0, 0.0);
    dpd_->file2_close(&XNJ);
    dpd_->buf4_close(&TIJAB);
    dpd_->buf4_sort(&WM, PSIF_EOM_TMP, qprs, 0, 7, "WnmjeDS_P");
    dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 0, 7, 2, 7, 0, SIJAB_lbl);
    dpd_->buf4_axpy(&WM, &SIJAB, -1.0);
    dpd_->buf4_close(&WM);
    dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WnmjeDS_P");
    dpd_->buf4_axpy(&WP, &SIJAB, 1.0);
    dpd_->buf4_close(&WP);
    dpd_->buf4_close(&SIJAB);

    /* Sijab -= Xnj * Tinab + Xni * Tjnab */
    dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 10, 17, 10, 17, 0, "WnmjeDS_MB");
    dpd_->buf4_init(&Tijab, PSIF_CC_TAMPS, H_IRR, 10, 17, 12, 17, 0, "tijab");
    dpd_->file2_init(&Xnj, PSIF_EOM_TMP, C_irr, 2, 2, "Xnj");
    dpd_->contract424(&Tijab, &Xnj, &WM, 1, 0, 1, 1.0, 0.0);
    dpd_->file2_close(&Xnj);
    dpd_->buf4_close(&Tijab);
    dpd_->buf4_sort(&WM, PSIF_EOM_TMP, qprs, 10, 17, "WnmjeDS_PB");
    dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 10, 17, 12, 17, 0, Sijab_lbl);
    dpd_->buf4_axpy(&WM, &Sijab, -1.0);
    dpd_->buf4_close(&WM);
    dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 10, 17, 10, 17, 0, "WnmjeDS_PB");
    dpd_->buf4_axpy(&WP, &Sijab, 1.0);
    dpd_->buf4_close(&WP);
    dpd_->buf4_close(&Sijab);

    /* SIjAb -= Xnj * tInAb + XNI * TjNAb */
    dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, SIjAb_lbl);
    dpd_->buf4_init(&TIjAb, PSIF_CC_TAMPS, H_IRR, 22, 28, 22, 28, 0, "tIjAb");
    dpd_->file2_init(&Xnj, PSIF_EOM_TMP, C_irr, 2, 2, "Xnj");
    dpd_->contract424(&TIjAb, &Xnj, &SIjAb, 1, 0, 1, -1.0, 1.0);
    dpd_->file2_close(&Xnj);
    dpd_->file2_init(&XNJ, PSIF_EOM_TMP, C_irr, 0, 0, "XNJ");
    dpd_->contract244(&XNJ, &TIjAb, &SIjAb, 0, 0, 0, -1.0, 1.0);
    dpd_->file2_close(&XNJ);
    dpd_->buf4_close(&TIjAb);
    dpd_->buf4_close(&SIjAb);
  }

#ifdef EOM_DEBUG
  check_sum("WnmjeDS",i,C_irr);
#endif
  return;
}

}} // namespace psi::cceom
