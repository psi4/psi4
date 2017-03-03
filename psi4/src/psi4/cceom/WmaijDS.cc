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
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/* This function computes the H-bar doubles-singles block contribution
   of Wmaij to a Sigma vector stored at Sigma plus 'i' */

void WmaijDS(int i, int C_irr) {
  dpdfile2 CME, Cme;
  dpdbuf4 SIJAB, Sijab, SIjAb, Z;
  dpdbuf4 W, WP, WM, WMBIJ, Wmbij, WMbIj, WmBiJ;
  char CME_lbl[32], Cme_lbl[32], SIJAB_lbl[32], Sijab_lbl[32], SIjAb_lbl[32];

  if (params.eom_ref == 0) { /* RHF */
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIjAb += WmAIj * Cmb - WMbIj * CMA */
    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WmaijDS Z(Ij,Ab)");
    global_dpd_->buf4_init(&WMbIj, PSIF_CC_HBAR, H_IRR, 10, 0, 10, 0, 0, "WMbIj");
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->contract244(&CME, &WMbIj, &Z, 0, 0, 1, 1.0, 0.0);
    global_dpd_->file2_close(&CME);
    global_dpd_->buf4_close(&WMbIj);
    global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP, qpsr, 0, 5, "WmaijDS Z(jI,bA)");
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    global_dpd_->buf4_axpy(&Z, &SIjAb,  -1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WmaijDS Z(jI,bA)");
    global_dpd_->buf4_axpy(&Z, &SIjAb,  -1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&SIjAb);
  }

  else if (params.eom_ref == 1) { /* ROHF */
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(Cme_lbl, "%s %d", "Cme", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIJAB += WMAIJ * CMB - WMBIJ * CMA */
    global_dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WmaijDS_M");
    global_dpd_->buf4_init(&WMBIJ, PSIF_CC_HBAR, H_IRR, 10, 2, 10, 2, 0, "WMBIJ");
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->contract244(&CME, &WMBIJ, &WM, 0, 0, 1, 1.0, 0.0);
    global_dpd_->file2_close(&CME);
    global_dpd_->buf4_close(&WMBIJ);
    global_dpd_->buf4_sort(&WM, PSIF_EOM_TMP, pqsr, 2, 5, "WmaijDS_P");
    global_dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WmaijDS_P");
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 2, 5, 2, 7, 0, SIJAB_lbl);
    global_dpd_->buf4_axpy(&WM, &SIJAB, -1.0);
    global_dpd_->buf4_axpy(&WP, &SIJAB,  1.0);
    global_dpd_->buf4_close(&SIJAB);
    global_dpd_->buf4_close(&WM);
    global_dpd_->buf4_close(&WP);

    /* Sijab += Wmaij * Cmb - Wmbij * Cma */
    global_dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WmaijDS_M");
    global_dpd_->buf4_init(&Wmbij, PSIF_CC_HBAR, H_IRR, 10, 2, 10, 2, 0, "Wmbij");
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, Cme_lbl);
    global_dpd_->contract244(&Cme, &Wmbij, &WM, 0, 0, 1, 1.0, 0.0);
    global_dpd_->file2_close(&Cme);
    global_dpd_->buf4_close(&Wmbij);
    global_dpd_->buf4_sort(&WM, PSIF_EOM_TMP, pqsr, 2, 5, "WmaijDS_P");
    global_dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WmaijDS_P");
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 2, 5, 2, 7, 0, Sijab_lbl);
    global_dpd_->buf4_axpy(&WM, &Sijab, -1.0);
    global_dpd_->buf4_axpy(&WP, &Sijab,  1.0);
    global_dpd_->buf4_close(&Sijab);
    global_dpd_->buf4_close(&WM);
    global_dpd_->buf4_close(&WP);

    /* SIjAb += WmAIj * Cmb - WMbIj * CMA */
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    global_dpd_->buf4_init(&WmBiJ, PSIF_CC_HBAR, H_IRR, 11, 0, 11, 0, 0, "WmBiJ (Bm,Ji)");
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, Cme_lbl);
    global_dpd_->contract424(&WmBiJ, &Cme, &SIjAb, 1, 0, 0, -1.0, 0.0);
    global_dpd_->file2_close(&Cme);
    global_dpd_->buf4_close(&WmBiJ);

    global_dpd_->buf4_init(&WMbIj, PSIF_CC_HBAR, H_IRR, 10, 0, 10, 0, 0, "WMbIj");
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->contract244(&CME, &WMbIj, &SIjAb, 0, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&CME);
    global_dpd_->buf4_close(&WMbIj);
    global_dpd_->buf4_close(&SIjAb);
  }

  else if (params.eom_ref == 2) { /* UHF */
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(Cme_lbl, "%s %d", "Cme", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIJAB += WMAIJ * CMB - WMBIJ * CMA */
    global_dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WmaijDS_M");
    global_dpd_->buf4_init(&WMBIJ, PSIF_CC_HBAR, H_IRR, 20, 2, 20, 2, 0, "WMBIJ");
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->contract244(&CME, &WMBIJ, &WM, 0, 0, 1, 1.0, 0.0);
    global_dpd_->file2_close(&CME);
    global_dpd_->buf4_close(&WMBIJ);
    global_dpd_->buf4_sort(&WM, PSIF_EOM_TMP, pqsr, 2, 5, "WmaijDS_P");
    global_dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "WmaijDS_P");
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 2, 5, 2, 7, 0, SIJAB_lbl);
    global_dpd_->buf4_axpy(&WM, &SIJAB, -1.0);
    global_dpd_->buf4_axpy(&WP, &SIJAB,  1.0);
    global_dpd_->buf4_close(&SIJAB);
    global_dpd_->buf4_close(&WM);
    global_dpd_->buf4_close(&WP);

    /* Sijab += Wmaij * Cmb - Wmbij * Cma */
    global_dpd_->buf4_init(&WM, PSIF_EOM_TMP, C_irr, 12, 15, 12, 15, 0, "WmaijDS_MB");
    global_dpd_->buf4_init(&Wmbij, PSIF_CC_HBAR, H_IRR, 30, 12, 30, 12, 0, "Wmbij");
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, Cme_lbl);
    global_dpd_->contract244(&Cme, &Wmbij, &WM, 0, 0, 1, 1.0, 0.0);
    global_dpd_->file2_close(&Cme);
    global_dpd_->buf4_close(&Wmbij);
    global_dpd_->buf4_sort(&WM, PSIF_EOM_TMP, pqsr, 12, 15, "WmaijDS_PB");
    global_dpd_->buf4_init(&WP, PSIF_EOM_TMP, C_irr, 12, 15, 12, 15, 0, "WmaijDS_PB");
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 12, 15, 12, 17, 0, Sijab_lbl);
    global_dpd_->buf4_axpy(&WM, &Sijab, -1.0);
    global_dpd_->buf4_axpy(&WP, &Sijab,  1.0);
    global_dpd_->buf4_close(&Sijab);
    global_dpd_->buf4_close(&WM);
    global_dpd_->buf4_close(&WP);

    /* SIjAb += WmAIj * Cmb - WMbIj * CMA */
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, SIjAb_lbl);
    global_dpd_->buf4_init(&WmBiJ, PSIF_CC_HBAR, H_IRR, 26, 22, 26, 22, 0, "WmBiJ (Bm,Ji)");
    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, Cme_lbl);
    global_dpd_->contract424(&WmBiJ, &Cme, &SIjAb, 1, 0, 0, -1.0, 0.0);
    global_dpd_->file2_close(&Cme);
    global_dpd_->buf4_close(&WmBiJ);

    global_dpd_->buf4_init(&WMbIj, PSIF_CC_HBAR, H_IRR, 24, 22, 24, 22, 0, "WMbIj");
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
    global_dpd_->contract244(&CME, &WMbIj, &SIjAb, 0, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&CME);
    global_dpd_->buf4_close(&WMbIj);
    global_dpd_->buf4_close(&SIjAb);
  }


#ifdef EOM_DEBUG
  check_sum("WmaijDS",i,C_irr);
#endif

  return;
}

}} // namespace psi::cceom
