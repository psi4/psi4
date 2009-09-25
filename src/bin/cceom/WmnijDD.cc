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
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_init(&WMnIj, CC_HBAR, H_IRR, 0, 0, 0, 0, 0, "WMnIj");
    dpd_contract444(&WMnIj, &CMnEf, &SIjAb, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&WMnIj);
    dpd_buf4_close(&SIjAb);
    dpd_buf4_close(&CMnEf);
  }

  else if (params.eom_ref == 1) { /* ROHF */
    sprintf(CMNEF_lbl, "%s %d", "CMNEF", i);
    sprintf(Cmnef_lbl, "%s %d", "Cmnef", i);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIJAB += WMNIJ * CMNAB */
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, CMNEF_lbl);
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_init(&WMNIJ, CC_HBAR, H_IRR, 2, 2, 2, 2, 0, "WMNIJ");
    dpd_contract444(&WMNIJ, &CMNEF, &SIJAB, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&WMNIJ);
    dpd_buf4_close(&SIJAB);
    dpd_buf4_close(&CMNEF);

    /* Sijab += Wmnij * Cmnab */
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 7, 2, 7, 0, Cmnef_lbl);
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 7, 2, 7, 0, Sijab_lbl);
    dpd_buf4_init(&Wmnij, CC_HBAR, H_IRR, 2, 2, 2, 2, 0, "Wmnij");
    dpd_contract444(&Wmnij, &Cmnef, &Sijab, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&Wmnij);
    dpd_buf4_close(&Sijab);
    dpd_buf4_close(&Cmnef);

    /* SIjAb += WMnIj * CMnAb */
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_init(&WMnIj, CC_HBAR, H_IRR, 0, 0, 0, 0, 0, "WMnIj");
    dpd_contract444(&WMnIj, &CMnEf, &SIjAb, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&WMnIj);
    dpd_buf4_close(&SIjAb);
    dpd_buf4_close(&CMnEf);
  }

  else if (params.eom_ref == 2) { /* UHF */
    sprintf(CMNEF_lbl, "%s %d", "CMNEF", i);
    sprintf(Cmnef_lbl, "%s %d", "Cmnef", i);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIJAB += WMNIJ * CMNAB */
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 7, 2, 7, 0, CMNEF_lbl);
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 7, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_init(&WMNIJ, CC_HBAR, H_IRR, 2, 2, 2, 2, 0, "WMNIJ");
    dpd_contract444(&WMNIJ, &CMNEF, &SIJAB, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&WMNIJ);
    dpd_buf4_close(&SIJAB);
    dpd_buf4_close(&CMNEF);

    /* Sijab += Wmnij * Cmnab */
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 12, 17, 12, 17, 0, Cmnef_lbl);
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 12, 17, 12, 17, 0, Sijab_lbl);
    dpd_buf4_init(&Wmnij, CC_HBAR, H_IRR, 12, 12, 12, 12, 0, "Wmnij");
    dpd_contract444(&Wmnij, &Cmnef, &Sijab, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&Wmnij);
    dpd_buf4_close(&Sijab);
    dpd_buf4_close(&Cmnef);

    /* SIjAb += WMnIj * CMnAb */
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, CMnEf_lbl);
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, SIjAb_lbl);
    dpd_buf4_init(&WMnIj, CC_HBAR, H_IRR, 22, 22, 22, 22, 0, "WMnIj");
    dpd_contract444(&WMnIj, &CMnEf, &SIjAb, 1, 1, 1.0, 1.0);
    dpd_buf4_close(&WMnIj);
    dpd_buf4_close(&SIjAb);
    dpd_buf4_close(&CMnEf);
  }

#ifdef EOM_DEBUG
  check_sum("WmnijDD",i,C_irr);
#endif
  return;
}

}} // namespace psi::cceom
