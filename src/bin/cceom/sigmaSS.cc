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
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "CME", i);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);

    /* SIA = FAE*CIE */
    dpd_file2_init(&FAE, CC_OEI, H_IRR, 1, 1, "FAE");
    dpd_contract222(&CME, &FAE, &SIA, 0, 0, 1.0, 0.0);
    dpd_file2_close(&FAE);

    /* SIA -= FMI*CMA */
    dpd_file2_init(&FMI, CC_OEI, H_IRR, 0, 0, "FMI");
    dpd_contract222(&FMI, &CME, &SIA, 1, 1, -1.0, 1.0);
    dpd_file2_close(&FMI);

    /* SIA += WMAEI*CME */
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "2 W(jb,ME) + W(Jb,Me)");
    dpd_contract422(&W, &CME, &SIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);

    dpd_file2_close(&CME);
    dpd_file2_close(&SIA);
  }

  else if (params.eom_ref == 1) { /* ROHF */
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "CME", i);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);

    sprintf(lbl, "%s %d", "Sia", i);
    dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "Cme", i);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);

    /* SIA = FAE*CIE */
    dpd_file2_init(&FAE, CC_OEI, H_IRR, 1, 1, "FAE");
    dpd_contract222(&CME, &FAE, &SIA, 0, 0, 1.0, 0.0);
    dpd_file2_close(&FAE);

    /* SIA -= FMI*CMA */
    dpd_file2_init(&FMI, CC_OEI, H_IRR, 0, 0, "FMI");
    dpd_contract222(&FMI, &CME, &SIA, 1, 1, -1.0, 1.0);
    dpd_file2_close(&FMI);

    /* SIA += WMAEI*CME */
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WMBEJ (JB,ME)");
    dpd_contract422(&W, &CME, &SIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);

    /* SIA += WmAeI*Cme */
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WmBeJ (JB,me)");
    dpd_contract422(&W, &Cme, &SIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);

    /* Sia = Fae*Cie */
    dpd_file2_init(&Fae, CC_OEI, H_IRR, 1, 1, "Fae");
    dpd_contract222(&Cme, &Fae, &Sia, 0, 0, 1.0, 0.0);
    dpd_file2_close(&Fae);

    /* Sia -= Fmi*Cma */
    dpd_file2_init(&Fmi, CC_OEI, H_IRR, 0, 0, "Fmi");
    dpd_contract222(&Fmi, &Cme, &Sia, 1, 1, -1.0, 1.0);
    dpd_file2_close(&Fmi);

    /* Sia += Wmaei*Cme */
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "Wmbej (jb,me)");
    dpd_contract422(&W,&Cme,&Sia, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);

    /* Sia += WMaEi*CME */
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 10, 10, 10, 10, 0, "WMbEj (jb,ME)");
    dpd_contract422(&W,&CME,&Sia, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);

    dpd_file2_close(&Cme);
    dpd_file2_close(&Sia);
    dpd_file2_close(&CME);
    dpd_file2_close(&SIA);
  }
  else { /* UHF */
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "CME", i);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);

    sprintf(lbl, "%s %d", "Sia", i);
    dpd_file2_init(&Sia, EOM_Sia, C_irr, 2, 3, lbl);
    sprintf(lbl, "%s %d", "Cme", i);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, lbl);

    /* SIA = FAE*CIE */
    dpd_file2_init(&FAE, CC_OEI, H_IRR, 1, 1, "FAE");
    dpd_contract222(&CME, &FAE, &SIA, 0, 0, 1.0, 0.0);
    dpd_file2_close(&FAE);

    /* SIA -= FMI*CMA */
    dpd_file2_init(&FMI, CC_OEI, H_IRR, 0, 0, "FMI");
    dpd_contract222(&FMI, &CME, &SIA, 1, 1, -1.0, 1.0);
    dpd_file2_close(&FMI);

    /* SIA += WMAEI*CME */
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 20, 20, 20, 20, 0, "WMBEJ (JB,ME)");
    dpd_contract422(&W, &CME, &SIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);

    /* SIA += WmAeI*Cme */
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 20, 30, 20, 30, 0, "WmBeJ (JB,me)");
    dpd_contract422(&W, &Cme, &SIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);

    /* Sia = Fae*Cie */
    dpd_file2_init(&Fae, CC_OEI, H_IRR, 3, 3, "Fae");
    dpd_contract222(&Cme, &Fae, &Sia, 0, 0, 1.0, 0.0);
    dpd_file2_close(&Fae);

    /* Sia -= Fmi*Cma */
    dpd_file2_init(&Fmi, CC_OEI, H_IRR, 2, 2, "Fmi");
    dpd_contract222(&Fmi, &Cme, &Sia, 1, 1, -1.0, 1.0);
    dpd_file2_close(&Fmi);

    /* Sia += Wmaei*Cme */
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 30, 30, 30, 30, 0, "Wmbej (jb,me)");
    dpd_contract422(&W,&Cme,&Sia, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);

    /* Sia += WMaEi*CME */
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 30, 20, 30, 20, 0, "WMbEj (jb,ME)");
    dpd_contract422(&W,&CME,&Sia, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&W);

    dpd_file2_close(&Cme);
    dpd_file2_close(&Sia);
    dpd_file2_close(&CME);
    dpd_file2_close(&SIA);
  }

#ifdef EOM_DEBUG
  check_sum("SigmaSS",i,C_irr);
#endif
  return;
}


}} // namespace psi::cceom
