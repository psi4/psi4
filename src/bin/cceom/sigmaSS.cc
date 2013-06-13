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

/* This function computes the H-bar singles-singles block contribution
   to a Sigma vector stored at Sigma plus 'i' */

void sigmaSS(int i, int C_irr) {
  dpdfile2 FMI, Fmi, FAE, Fae, Cjunk;
  dpdfile2 CME, Cme, SIA, Sia;
  dpdbuf4 W;
  char lbl[32];

  if (params.eom_ref == 0) { /* RHF */
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "CME", i);
    dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);

    /* SIA = FAE*CIE */
    dpd_->file2_init(&FAE, PSIF_CC_OEI, H_IRR, 1, 1, "FAE");
    dpd_->contract222(&CME, &FAE, &SIA, 0, 0, 1.0, 0.0);
    dpd_->file2_close(&FAE);

    /* SIA -= FMI*CMA */
    dpd_->file2_init(&FMI, PSIF_CC_OEI, H_IRR, 0, 0, "FMI");
    dpd_->contract222(&FMI, &CME, &SIA, 1, 1, -1.0, 1.0);
    dpd_->file2_close(&FMI);

    /* SIA += WMAEI*CME */
    dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "2 W(jb,ME) + W(Jb,Me)");
    dpd_->contract422(&W, &CME, &SIA, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&W);

    dpd_->file2_close(&CME);
    dpd_->file2_close(&SIA);
  }

  else if (params.eom_ref == 1) { /* ROHF */
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "CME", i);
    dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);

    sprintf(lbl, "%s %d", "Sia", i);
    dpd_->file2_init(&Sia, PSIF_EOM_Sia, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "Cme", i);
    dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, lbl);

    /* SIA = FAE*CIE */
    dpd_->file2_init(&FAE, PSIF_CC_OEI, H_IRR, 1, 1, "FAE");
    dpd_->contract222(&CME, &FAE, &SIA, 0, 0, 1.0, 0.0);
    dpd_->file2_close(&FAE);

    /* SIA -= FMI*CMA */
    dpd_->file2_init(&FMI, PSIF_CC_OEI, H_IRR, 0, 0, "FMI");
    dpd_->contract222(&FMI, &CME, &SIA, 1, 1, -1.0, 1.0);
    dpd_->file2_close(&FMI);

    /* SIA += WMAEI*CME */
    dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WMBEJ (JB,ME)");
    dpd_->contract422(&W, &CME, &SIA, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&W);

    /* SIA += WmAeI*Cme */
    dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WmBeJ (JB,me)");
    dpd_->contract422(&W, &Cme, &SIA, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&W);

    /* Sia = Fae*Cie */
    dpd_->file2_init(&Fae, PSIF_CC_OEI, H_IRR, 1, 1, "Fae");
    dpd_->contract222(&Cme, &Fae, &Sia, 0, 0, 1.0, 0.0);
    dpd_->file2_close(&Fae);

    /* Sia -= Fmi*Cma */
    dpd_->file2_init(&Fmi, PSIF_CC_OEI, H_IRR, 0, 0, "Fmi");
    dpd_->contract222(&Fmi, &Cme, &Sia, 1, 1, -1.0, 1.0);
    dpd_->file2_close(&Fmi);

    /* Sia += Wmaei*Cme */
    dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "Wmbej (jb,me)");
    dpd_->contract422(&W,&Cme,&Sia, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&W);

    /* Sia += WMaEi*CME */
    dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WMbEj (jb,ME)");
    dpd_->contract422(&W,&CME,&Sia, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&W);

    dpd_->file2_close(&Cme);
    dpd_->file2_close(&Sia);
    dpd_->file2_close(&CME);
    dpd_->file2_close(&SIA);
  }
  else { /* UHF */
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "CME", i);
    dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);

    sprintf(lbl, "%s %d", "Sia", i);
    dpd_->file2_init(&Sia, PSIF_EOM_Sia, C_irr, 2, 3, lbl);
    sprintf(lbl, "%s %d", "Cme", i);
    dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, lbl);

    /* SIA = FAE*CIE */
    dpd_->file2_init(&FAE, PSIF_CC_OEI, H_IRR, 1, 1, "FAE");
    dpd_->contract222(&CME, &FAE, &SIA, 0, 0, 1.0, 0.0);
    dpd_->file2_close(&FAE);

    /* SIA -= FMI*CMA */
    dpd_->file2_init(&FMI, PSIF_CC_OEI, H_IRR, 0, 0, "FMI");
    dpd_->contract222(&FMI, &CME, &SIA, 1, 1, -1.0, 1.0);
    dpd_->file2_close(&FMI);

    /* SIA += WMAEI*CME */
    dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 20, 20, 20, 20, 0, "WMBEJ (JB,ME)");
    dpd_->contract422(&W, &CME, &SIA, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&W);

    /* SIA += WmAeI*Cme */
    dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 20, 30, 20, 30, 0, "WmBeJ (JB,me)");
    dpd_->contract422(&W, &Cme, &SIA, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&W);

    /* Sia = Fae*Cie */
    dpd_->file2_init(&Fae, PSIF_CC_OEI, H_IRR, 3, 3, "Fae");
    dpd_->contract222(&Cme, &Fae, &Sia, 0, 0, 1.0, 0.0);
    dpd_->file2_close(&Fae);

    /* Sia -= Fmi*Cma */
    dpd_->file2_init(&Fmi, PSIF_CC_OEI, H_IRR, 2, 2, "Fmi");
    dpd_->contract222(&Fmi, &Cme, &Sia, 1, 1, -1.0, 1.0);
    dpd_->file2_close(&Fmi);

    /* Sia += Wmaei*Cme */
    dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 30, 30, 30, 30, 0, "Wmbej (jb,me)");
    dpd_->contract422(&W,&Cme,&Sia, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&W);

    /* Sia += WMaEi*CME */
    dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 30, 20, 30, 20, 0, "WMbEj (jb,ME)");
    dpd_->contract422(&W,&CME,&Sia, 0, 0, 1.0, 1.0);
    dpd_->buf4_close(&W);

    dpd_->file2_close(&Cme);
    dpd_->file2_close(&Sia);
    dpd_->file2_close(&CME);
    dpd_->file2_close(&SIA);
  }

#ifdef EOM_DEBUG
  check_sum("SigmaSS",i,C_irr);
#endif
  return;
}


}} // namespace psi::cceom
