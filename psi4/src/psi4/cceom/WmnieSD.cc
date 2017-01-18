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
    \brief Computes the H-bar SD block contribution of Wmnie to a Sigma vector[i]. 
*/
#include <cstdio>
#include <cmath>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/* This function computes the H-bar singles-doubles block contribution
   of Wmnie to a Sigma vector stored at Sigma plus 'i' */

void WmnieSD(int i, int C_irr) {
  dpdfile2 SIA, Sia;
  dpdbuf4 CMNEF, Cmnef, CMnEf, CmNeF;
  dpdbuf4 WMNIE, Wmnie, WMnIe, WmNiE;
  char lbl[32];

  if (params.eom_ref == 0) { /* RHF */
    sprintf(lbl, "%s %d", "SIA", i);
    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "CMnEf", i);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&WMnIe, PSIF_CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "WMnIe - 2WnMIe (Mn,eI)");
    /* dpd_buf4_print(&WMnIe,outfile,1);
       dpd_buf4_print(&CMnEf,outfile,1);
       outfile->Printf(stdout,"starting Wmnie*CMNEF ->SIA\n");
       outfile->Printf("starting Wmnie*CMNEF ->SIA\n"); */
    global_dpd_->contract442(&WMnIe, &CMnEf, &SIA, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_close(&WMnIe);
    global_dpd_->file2_close(&SIA);
  }

  else if (params.eom_ref == 1) { /* ROHF */
    sprintf(lbl, "%s %d", "SIA", i);
    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "Sia", i);
    global_dpd_->file2_init(&Sia, PSIF_EOM_Sia, C_irr, 0, 1, lbl);

    /* SIA += 0.5 WMNIE * CMNAE + WMnIe * CMnAe */
    global_dpd_->buf4_init(&WMNIE, PSIF_CC_HBAR, H_IRR, 2, 11, 2, 11, 0, "WMNIE (M>N,EI)");
    sprintf(lbl, "%s %d", "CMNEF", i);
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 2, 5, 2, 7, 0, lbl);
    global_dpd_->contract442(&WMNIE, &CMNEF, &SIA, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&CMNEF);
    global_dpd_->buf4_close(&WMNIE);

    global_dpd_->buf4_init(&WMnIe, PSIF_CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CMnfE");
    global_dpd_->contract442(&WMnIe, &CMnEf, &SIA, 3, 3, -1.0, 1.0);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_close(&WMnIe);

    /* Sia += 0.5 Wmnie * Cmnae + Wmnie * Cmnae */
    global_dpd_->buf4_init(&Wmnie, PSIF_CC_HBAR, H_IRR, 2, 11, 2, 11, 0, "Wmnie (m>n,ei)");
    sprintf(lbl, "%s %d", "Cmnef", i);
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 2, 5, 2, 7, 0, lbl);
    global_dpd_->contract442(&Wmnie, &Cmnef, &Sia, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&Cmnef);
    global_dpd_->buf4_close(&Wmnie);

    global_dpd_->buf4_init(&WmNiE, PSIF_CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "WmNiE (mN,Ei)");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CnMEf");
    global_dpd_->contract442(&WmNiE, &CMnEf, &Sia, 3, 3, -1.0, 1.0);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_close(&WmNiE);

    global_dpd_->file2_close(&SIA);
    global_dpd_->file2_close(&Sia);
  }

  else { /* UHF */
    sprintf(lbl, "%s %d", "SIA", i);
    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "Sia", i);
    global_dpd_->file2_init(&Sia, PSIF_EOM_Sia, C_irr, 2, 3, lbl);

    /* SIA += 0.5 WMNIE * CMNAE + WMnIe * CMnAe */
    global_dpd_->buf4_init(&WMNIE, PSIF_CC_HBAR, H_IRR, 2, 21, 2, 21, 0, "WMNIE (M>N,EI)");
    sprintf(lbl, "%s %d", "CMNEF", i);
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 2, 5, 2, 7, 0, lbl);
    global_dpd_->contract442(&WMNIE, &CMNEF, &SIA, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&CMNEF);
    global_dpd_->buf4_close(&WMNIE);

    global_dpd_->buf4_init(&WMnIe, PSIF_CC_HBAR, H_IRR, 22, 25, 22, 25, 0, "WMnIe (Mn,eI)");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_TMP, C_irr, 22, 29, 22, 29, 0, "CMnfE");
    global_dpd_->contract442(&WMnIe, &CMnEf, &SIA, 3, 3, -1.0, 1.0);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_close(&WMnIe);

    /* Sia += 0.5 Wmnie * Cmnae + Wmnie * Cmnae */
    global_dpd_->buf4_init(&Wmnie, PSIF_CC_HBAR, H_IRR, 12, 31, 12, 31, 0, "Wmnie (m>n,ei)");
    sprintf(lbl, "%s %d", "Cmnef", i);
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 12, 15, 12, 17, 0, lbl);
    global_dpd_->contract442(&Wmnie, &Cmnef, &Sia, 3, 3, 1.0, 1.0);
    global_dpd_->buf4_close(&Cmnef);
    global_dpd_->buf4_close(&Wmnie);

    global_dpd_->buf4_init(&WmNiE, PSIF_CC_HBAR, H_IRR, 23, 26, 23, 26, 0, "WmNiE (mN,Ei)");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_TMP, C_irr, 23, 28, 23, 28, 0, "CnMEf");
    global_dpd_->contract442(&WmNiE, &CMnEf, &Sia, 3, 3, -1.0, 1.0);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_close(&WmNiE);

    global_dpd_->file2_close(&SIA);
    global_dpd_->file2_close(&Sia);
  }

#ifdef EOM_DEBUG
  check_sum("WmnieSD",i,C_irr);
#endif
  return;
}


}} // namespace psi::cceom
