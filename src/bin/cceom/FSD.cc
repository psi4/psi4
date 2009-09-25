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
   to a Sigma vector stored at Sigma plus 'i' */

void FSD(int i, int C_irr) {
  dpdfile2 Fme, FME;
  dpdfile2 SIA, Sia;
  dpdbuf4 CMNEF, Cmnef, CMnEf, CmNeF;
  char lbl[32];

  if (params.eom_ref == 0) {  /* RHF */
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    dpd_file2_init(&FME, CC_OEI, H_IRR, 0, 1, "FME");
    dpd_buf4_init(&CMnEf, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "2CMnEf - CMnfE");
    dpd_dot24(&FME,&CMnEf,&SIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&CMnEf);
    dpd_file2_close(&FME);
    dpd_file2_close(&SIA);
  }

  else if (params.eom_ref == 1) { /* ROHF */
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "Sia", i);
    dpd_file2_init(&Sia, EOM_Sia, C_irr, 0, 1, lbl);

    /* SIA += FME*CIMAE + Fme*CImAe */
    dpd_file2_init(&FME, CC_OEI, H_IRR, 0, 1, "FME");
    sprintf(lbl, "%s %d", "CMNEF", i);
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 0, 5, 2, 7, 0, lbl);
    dpd_dot24(&FME,&CMNEF,&SIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&CMNEF);
    dpd_file2_close(&FME);

    dpd_file2_init(&Fme, CC_OEI, H_IRR, 0, 1, "Fme");
    sprintf(lbl, "%s %d", "CMnEf", i);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    dpd_dot24(&Fme,&CMnEf,&SIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&CMnEf);
    dpd_file2_close(&Fme);

    /* Sia += Fme*Cimae + FME*CiMaE */
    dpd_file2_init(&Fme, CC_OEI, H_IRR, 0, 1, "Fme");
    sprintf(lbl, "%s %d", "Cmnef", i);
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 0, 5, 2, 7, 0, lbl);
    dpd_dot24(&Fme,&Cmnef,&Sia, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&Cmnef);
    dpd_file2_close(&Fme);

    dpd_file2_init(&FME, CC_OEI, H_IRR, 0, 1, "FME");
    dpd_buf4_init(&CmNeF, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CmNeF");
    dpd_dot24(&FME,&CmNeF,&Sia, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&CmNeF);
    dpd_file2_close(&FME);

    dpd_file2_close(&Sia);
    dpd_file2_close(&SIA);
  }

  else { /* UHF */

    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "Sia", i);
    dpd_file2_init(&Sia, EOM_Sia, C_irr, 2, 3, lbl);

    /* SIA += FME*CIMAE + Fme*CImAe */
    dpd_file2_init(&FME, CC_OEI, H_IRR, 0, 1, "FME");
    sprintf(lbl, "%s %d", "CMNEF", i);
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 0, 5, 2, 7, 0, lbl);
    dpd_dot24(&FME,&CMNEF,&SIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&CMNEF);
    dpd_file2_close(&FME);

    dpd_file2_init(&Fme, CC_OEI, H_IRR, 2, 3, "Fme");
    sprintf(lbl, "%s %d", "CMnEf", i);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, lbl);
    dpd_dot24(&Fme,&CMnEf,&SIA, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&CMnEf);
    dpd_file2_close(&Fme);

    /* Sia += Fme*Cimae + FME*CiMaE */
    dpd_file2_init(&Fme, CC_OEI, H_IRR, 2, 3, "Fme");
    sprintf(lbl, "%s %d", "Cmnef", i);
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 10, 15, 12, 17, 0, lbl);
    dpd_dot24(&Fme,&Cmnef,&Sia, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&Cmnef);
    dpd_file2_close(&Fme);

    dpd_file2_init(&FME, CC_OEI, H_IRR, 0, 1, "FME");
    dpd_buf4_init(&CmNeF, EOM_TMP, C_irr, 23, 29, 23, 29, 0, "CmNeF");
    dpd_dot24(&FME,&CmNeF,&Sia, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&CmNeF);
    dpd_file2_close(&FME);

    dpd_file2_close(&Sia);
    dpd_file2_close(&SIA);
  }

#ifdef EOM_DEBUG
  check_sum("FSD    ",i,C_irr);
#endif

  return;
}

}} // namespace psi::cceom
