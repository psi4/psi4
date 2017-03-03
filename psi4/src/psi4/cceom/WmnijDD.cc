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

/* This function computes the H-bar doubles-doubles block contribution
   from Wmnij*Cmnab to a Sigma vector stored at Sigma plus 'i' */

void WmnijDD(int i, int C_irr) {
  dpdbuf4 SIJAB, Sijab, SIjAb;
  dpdbuf4 CMNEF, Cmnef, CMnEf;
  dpdbuf4 WMNIJ, Wmnij, WMnIj;
  char CMNEF_lbl[32], Cmnef_lbl[32], CMnEf_lbl[32];
  char SIJAB_lbl[32], Sijab_lbl[32], SIjAb_lbl[32];

  if (params.eom_ref == 0) { /* RHF */
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIjAb += WMnIj * CMnAb */
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, H_IRR, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->contract444(&WMnIj, &CMnEf, &SIjAb, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&WMnIj);
    global_dpd_->buf4_close(&SIjAb);
    global_dpd_->buf4_close(&CMnEf);
  }

  else if (params.eom_ref == 1) { /* ROHF */
    sprintf(CMNEF_lbl, "%s %d", "CMNEF", i);
    sprintf(Cmnef_lbl, "%s %d", "Cmnef", i);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIJAB += WMNIJ * CMNAB */
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, CMNEF_lbl);
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, SIJAB_lbl);
    global_dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, H_IRR, 2, 2, 2, 2, 0, "WMNIJ");
    global_dpd_->contract444(&WMNIJ, &CMNEF, &SIJAB, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&WMNIJ);
    global_dpd_->buf4_close(&SIJAB);
    global_dpd_->buf4_close(&CMNEF);

    /* Sijab += Wmnij * Cmnab */
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, Cmnef_lbl);
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 2, 7, 2, 7, 0, Sijab_lbl);
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, H_IRR, 2, 2, 2, 2, 0, "Wmnij");
    global_dpd_->contract444(&Wmnij, &Cmnef, &Sijab, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_close(&Sijab);
    global_dpd_->buf4_close(&Cmnef);

    /* SIjAb += WMnIj * CMnAb */
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, H_IRR, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->contract444(&WMnIj, &CMnEf, &SIjAb, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&WMnIj);
    global_dpd_->buf4_close(&SIjAb);
    global_dpd_->buf4_close(&CMnEf);
  }

  else if (params.eom_ref == 2) { /* UHF */
    sprintf(CMNEF_lbl, "%s %d", "CMNEF", i);
    sprintf(Cmnef_lbl, "%s %d", "Cmnef", i);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIJAB += WMNIJ * CMNAB */
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, CMNEF_lbl);
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, SIJAB_lbl);
    global_dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, H_IRR, 2, 2, 2, 2, 0, "WMNIJ");
    global_dpd_->contract444(&WMNIJ, &CMNEF, &SIJAB, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&WMNIJ);
    global_dpd_->buf4_close(&SIJAB);
    global_dpd_->buf4_close(&CMNEF);

    /* Sijab += Wmnij * Cmnab */
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 12, 17, 12, 17, 0, Cmnef_lbl);
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 12, 17, 12, 17, 0, Sijab_lbl);
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, H_IRR, 12, 12, 12, 12, 0, "Wmnij");
    global_dpd_->contract444(&Wmnij, &Cmnef, &Sijab, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_close(&Sijab);
    global_dpd_->buf4_close(&Cmnef);

    /* SIjAb += WMnIj * CMnAb */
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, CMnEf_lbl);
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, SIjAb_lbl);
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, H_IRR, 22, 22, 22, 22, 0, "WMnIj");
    global_dpd_->contract444(&WMnIj, &CMnEf, &SIjAb, 1, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&WMnIj);
    global_dpd_->buf4_close(&SIjAb);
    global_dpd_->buf4_close(&CMnEf);
  }

#ifdef EOM_DEBUG
  check_sum("WmnijDD",i,C_irr);
#endif
  return;
}

}} // namespace psi::cceom
