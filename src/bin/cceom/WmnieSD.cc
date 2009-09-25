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

/* This function computes the H-bar singles-doubles block contribution
   of Wmnie to a Sigma vector stored at Sigma plus 'i' */

void WmnieSD(int i, int C_irr) {
  dpdfile2 SIA, Sia;
  dpdbuf4 CMNEF, Cmnef, CMnEf, CmNeF;
  dpdbuf4 WMNIE, Wmnie, WMnIe, WmNiE;
  char lbl[32];

  if (params.eom_ref == 0) { /* RHF */
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "CMnEf", i);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_init(&WMnIe, CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "WMnIe - 2WnMIe (Mn,eI)");
    /* dpd_buf4_print(&WMnIe,outfile,1);
       dpd_buf4_print(&CMnEf,outfile,1);
       fprintf(stdout,"starting Wmnie*CMNEF ->SIA\n");
       fprintf(outfile,"starting Wmnie*CMNEF ->SIA\n"); */
    dpd_contract442(&WMnIe, &CMnEf, &SIA, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_close(&WMnIe);
    dpd_file2_close(&SIA);
  }

  else if (params.eom_ref == 1) { /* ROHF */
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "Sia", i);
    dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, lbl);

    /* SIA += 0.5 WMNIE * CMNAE + WMnIe * CMnAe */
    dpd_buf4_init(&WMNIE, CC_HBAR, H_IRR, 2, 11, 2, 11, 0, "WMNIE (M>N,EI)");
    sprintf(lbl, "%s %d", "CMNEF", i);
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 5, 2, 7, 0, lbl);
    dpd_contract442(&WMNIE, &CMNEF, &SIA, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&CMNEF);
    dpd_buf4_close(&WMNIE);

    dpd_buf4_init(&WMnIe, CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
    dpd_buf4_init(&CMnEf, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CMnfE");
    dpd_contract442(&WMnIe, &CMnEf, &SIA, 3, 3, -1.0, 1.0);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_close(&WMnIe);

    /* Sia += 0.5 Wmnie * Cmnae + Wmnie * Cmnae */
    dpd_buf4_init(&Wmnie, CC_HBAR, H_IRR, 2, 11, 2, 11, 0, "Wmnie (m>n,ei)");
    sprintf(lbl, "%s %d", "Cmnef", i);
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 5, 2, 7, 0, lbl);
    dpd_contract442(&Wmnie, &Cmnef, &Sia, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&Cmnef);
    dpd_buf4_close(&Wmnie);

    dpd_buf4_init(&WmNiE, CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "WmNiE (mN,Ei)");
    dpd_buf4_init(&CMnEf, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CnMEf");
    dpd_contract442(&WmNiE, &CMnEf, &Sia, 3, 3, -1.0, 1.0);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_close(&WmNiE);

    dpd_file2_close(&SIA);
    dpd_file2_close(&Sia);
  }

  else { /* UHF */
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "Sia", i);
    dpd_file2_init(&Sia, EOM_Sia, C_irr, 2, 3, lbl);

    /* SIA += 0.5 WMNIE * CMNAE + WMnIe * CMnAe */
    dpd_buf4_init(&WMNIE, CC_HBAR, H_IRR, 2, 21, 2, 21, 0, "WMNIE (M>N,EI)");
    sprintf(lbl, "%s %d", "CMNEF", i);
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 5, 2, 7, 0, lbl);
    dpd_contract442(&WMNIE, &CMNEF, &SIA, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&CMNEF);
    dpd_buf4_close(&WMNIE);

    dpd_buf4_init(&WMnIe, CC_HBAR, H_IRR, 22, 25, 22, 25, 0, "WMnIe (Mn,eI)");
    dpd_buf4_init(&CMnEf, EOM_TMP, C_irr, 22, 29, 22, 29, 0, "CMnfE");
    dpd_contract442(&WMnIe, &CMnEf, &SIA, 3, 3, -1.0, 1.0);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_close(&WMnIe);

    /* Sia += 0.5 Wmnie * Cmnae + Wmnie * Cmnae */
    dpd_buf4_init(&Wmnie, CC_HBAR, H_IRR, 12, 31, 12, 31, 0, "Wmnie (m>n,ei)");
    sprintf(lbl, "%s %d", "Cmnef", i);
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 12, 15, 12, 17, 0, lbl);
    dpd_contract442(&Wmnie, &Cmnef, &Sia, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&Cmnef);
    dpd_buf4_close(&Wmnie);

    dpd_buf4_init(&WmNiE, CC_HBAR, H_IRR, 23, 26, 23, 26, 0, "WmNiE (mN,Ei)");
    dpd_buf4_init(&CMnEf, EOM_TMP, C_irr, 23, 28, 23, 28, 0, "CnMEf");
    dpd_contract442(&WmNiE, &CMnEf, &Sia, 3, 3, -1.0, 1.0);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_close(&WmNiE);

    dpd_file2_close(&SIA);
    dpd_file2_close(&Sia);
  }

#ifdef EOM_DEBUG
  check_sum("WmnieSD",i,C_irr);
#endif
  return;
}


}} // namespace psi::cceom
