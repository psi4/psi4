/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/

/*! \defgroup CCEOM cceom: Equation-of-Motion Coupled-Cluster */

#include <cstdio>
#include <cmath>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/* This function computes the H-bar doubles-doubles block contribution
   to a Sigma vector stored at Sigma plus 'i' */

void FDD(int i, int C_irr) {
  dpdfile2 FAE, Fae, FMI, Fmi;
  dpdbuf4 SIJAB, Sijab, SIjAb, FP, FM;
  dpdbuf4 CMNEF, Cmnef, CMnEf, X, X2, Z, Z2, Z3;
  char CMNEF_lbl[32], Cmnef_lbl[32], CMnEf_lbl[32], CmNeF_lbl[32];
  char SIJAB_lbl[32], Sijab_lbl[32], SIjAb_lbl[32];

  if (params.eom_ref == 0) { /* RHF */
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIjAb += Fbe * CIjAe - FAE * CIjbE */
    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "FDD_Fbe Z(Ij,Ab)");
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    dpd_file2_init(&FAE, CC_OEI, H_IRR, 1, 1, "FAE");
    dpd_contract424(&CMnEf, &FAE, &Z, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&FAE);
    dpd_buf4_close(&CMnEf);

    dpd_buf4_sort(&Z, EOM_TMP, qpsr, 0, 5, "FDD_Fbe Z(jI,bA)");
    dpd_buf4_init(&Z2, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "FDD_Fbe Z(jI,bA)");

    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_axpy(&Z, &SIjAb, 1.0);
    dpd_buf4_axpy(&Z2, &SIjAb, 1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&Z2);
    dpd_buf4_close(&SIjAb);

#ifdef EOM_DEBUG
    check_sum("FDD_Fbe",i,C_irr);
#endif

    /* SIjAb -= FMJ * CImAb - FMI * CjMAb */
    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "FDD_Fmj Z(Ij,Ab)");
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    dpd_file2_init(&FMI, CC_OEI, H_IRR, 0, 0, "FMI");
    dpd_contract244(&FMI, &CMnEf, &Z, 0, 0, 0, 1.0, 0.0);
    dpd_file2_close(&FMI);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_sort(&Z, EOM_TMP, qpsr, 0, 5, "FDD_Fmj Z(jI,bA)"); 
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_axpy(&Z, &SIjAb, -1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "FDD_Fmj Z(jI,bA)");
    dpd_buf4_axpy(&Z, &SIjAb, -1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&SIjAb);

#ifdef EOM_DEBUG
    check_sum("FDD_Fmj",i,C_irr);
#endif
  }

  else if (params.eom_ref == 1) { /* ROHF */
    sprintf(CMNEF_lbl, "%s %d", "CMNEF", i);
    sprintf(Cmnef_lbl, "%s %d", "Cmnef", i);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    sprintf(CmNeF_lbl, "%s %d", "CmNeF", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIJAB += FBE * CIJAE - FAE * CIJBE */
    dpd_buf4_init(&FP, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "FDD_FBEP");
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 5, 2, 7, 0, CMNEF_lbl);
    dpd_file2_init(&FAE, CC_OEI, H_IRR, 1, 1, "FAE");
    dpd_contract424(&CMNEF, &FAE, &FP, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&FAE);
    dpd_buf4_close(&CMNEF);
    dpd_buf4_sort(&FP, EOM_TMP, pqsr, 2, 5, "FDD_FBEM");
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 5, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_axpy(&FP, &SIJAB, 1.0);
    dpd_buf4_close(&FP);
    dpd_buf4_init(&FM, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "FDD_FBEM");
    dpd_buf4_axpy(&FM, &SIJAB, -1.0);
    dpd_buf4_close(&FM);
    dpd_buf4_close(&SIJAB);

    /* Sijab += Fbe * Cijae - Fae * Cijbe */
    dpd_buf4_init(&FP, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "FDD_FBEP");
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 2, 5, 2, 7, 0, Cmnef_lbl);
    dpd_file2_init(&Fae, CC_OEI, H_IRR, 1, 1, "Fae");
    dpd_contract424(&Cmnef, &Fae, &FP, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&Fae);
    dpd_buf4_close(&Cmnef);
    dpd_buf4_sort(&FP, EOM_TMP, pqsr, 2, 5, "FDD_FBEM");
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 2, 5, 2, 7, 0, Sijab_lbl);
    dpd_buf4_axpy(&FP, &Sijab, 1.0);
    dpd_buf4_close(&FP);
    dpd_buf4_init(&FM, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "FDD_FBEM");
    dpd_buf4_axpy(&FM, &Sijab, -1.0);
    dpd_buf4_close(&FM);
    dpd_buf4_close(&Sijab);

    /* SIjAb += Fbe * CIjAe - FAE * CIjbE */
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    dpd_file2_init(&Fae, CC_OEI, H_IRR, 1, 1, "Fae");
    dpd_contract424(&CMnEf, &Fae, &SIjAb, 3, 1, 0, 1.0, 1.0);
    dpd_file2_close(&Fae);
    dpd_file2_init(&FAE, CC_OEI, H_IRR, 1, 1, "FAE");
    dpd_contract244(&FAE, &CMnEf, &SIjAb, 1, 2, 1, 1.0, 1.0);
    dpd_file2_close(&FAE);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_close(&SIjAb);

#ifdef EOM_DEBUG
    check_sum("Fbe_FDD ",i,C_irr);
#endif

    /* SIJAB -= FMJ * CIMAB - FMI * CJMAB */
    dpd_buf4_init(&FM, EOM_TMP, C_irr, 0, 7, 0, 7, 0, "FDD_FMJM");
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 0, 7, 2, 7, 0, CMNEF_lbl);
    dpd_file2_init(&FMI, CC_OEI, H_IRR, 0, 0, "FMI");
    dpd_contract424(&CMNEF, &FMI, &FM, 1, 0, 1, 1.0, 0.0);
    dpd_file2_close(&FMI);
    dpd_buf4_close(&CMNEF);
    dpd_buf4_sort(&FM, EOM_TMP, qprs, 0, 7, "FDD_FMJP");
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 0, 7, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_axpy(&FM, &SIJAB, -1.0);
    dpd_buf4_close(&FM);
    dpd_buf4_init(&FP, EOM_TMP, C_irr, 0, 7, 0, 7, 0, "FDD_FMJP");
    dpd_buf4_axpy(&FP, &SIJAB, 1.0);
    dpd_buf4_close(&FP);
    dpd_buf4_close(&SIJAB);

    /* Sijab -= Fmj * Cimab - Fmi * Cjmab */
    dpd_buf4_init(&FM, EOM_TMP, C_irr, 0, 7, 0, 7, 0, "FDD_FMJM");
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 0, 7, 2, 7, 0, Cmnef_lbl);
    dpd_file2_init(&Fmi, CC_OEI, H_IRR, 0, 0, "Fmi");
    dpd_contract424(&Cmnef, &Fmi, &FM, 1, 0, 1, 1.0, 0.0);
    dpd_file2_close(&Fmi);
    dpd_buf4_close(&Cmnef);
    dpd_buf4_sort(&FM, EOM_TMP, qprs, 0, 7, "FDD_FMJP");
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 0, 7, 2, 7, 0, Sijab_lbl);
    dpd_buf4_axpy(&FM, &Sijab, -1.0);
    dpd_buf4_close(&FM);
    dpd_buf4_init(&FP, EOM_TMP, C_irr, 0, 7, 0, 7, 0, "FDD_FMJP");
    dpd_buf4_axpy(&FP, &Sijab, 1.0);
    dpd_buf4_close(&FP);
    dpd_buf4_close(&Sijab);

    /* SIjAb -= Fmj * CImAb - FMI * CjMAb */
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    dpd_file2_init(&Fmi, CC_OEI, H_IRR, 0, 0, "Fmi");
    dpd_contract424(&CMnEf, &Fmi, &SIjAb, 1, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Fmi);
    dpd_file2_init(&FMI, CC_OEI, H_IRR, 0, 0, "FMI");
    dpd_contract244(&FMI, &CMnEf, &SIjAb, 0, 0, 0, -1.0, 1.0);
    dpd_file2_close(&FMI);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_close(&SIjAb);

#ifdef EOM_DEBUG
    check_sum("Fmj_DD",i,C_irr);
#endif
  }

  else if (params.eom_ref == 2) { /* UHF */
    sprintf(CMNEF_lbl, "%s %d", "CMNEF", i);
    sprintf(Cmnef_lbl, "%s %d", "Cmnef", i);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    sprintf(CmNeF_lbl, "%s %d", "CmNeF", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* SIJAB += FBE * CIJAE - FAE * CIJBE */
    dpd_buf4_init(&FP, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "FDD_FBEP");
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 2, 5, 2, 7, 0, CMNEF_lbl);
    dpd_file2_init(&FAE, CC_OEI, H_IRR, 1, 1, "FAE");
    dpd_contract424(&CMNEF, &FAE, &FP, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&FAE);
    dpd_buf4_close(&CMNEF);
    dpd_buf4_sort(&FP, EOM_TMP, pqsr, 2, 5, "FDD_FBEM");
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 2, 5, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_axpy(&FP, &SIJAB, 1.0);
    dpd_buf4_close(&FP);
    dpd_buf4_init(&FM, EOM_TMP, C_irr, 2, 5, 2, 5, 0, "FDD_FBEM");
    dpd_buf4_axpy(&FM, &SIJAB, -1.0);
    dpd_buf4_close(&FM);
    dpd_buf4_close(&SIJAB);

    /* Sijab += Fbe * Cijae - Fae * Cijbe */
    dpd_buf4_init(&FP, EOM_TMP, C_irr, 12, 15, 12, 15, 0, "FDD_FbePB");
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 12, 15, 12, 17, 0, Cmnef_lbl);
    dpd_file2_init(&Fae, CC_OEI, H_IRR, 3, 3, "Fae");
    dpd_contract424(&Cmnef, &Fae, &FP, 3, 1, 0, 1.0, 0.0);
    dpd_file2_close(&Fae);
    dpd_buf4_close(&Cmnef);
    dpd_buf4_sort(&FP, EOM_TMP, pqsr, 12, 15, "FDD_FbeMB");
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 12, 15, 12, 17, 0, Sijab_lbl);
    dpd_buf4_axpy(&FP, &Sijab, 1.0);
    dpd_buf4_close(&FP);
    dpd_buf4_init(&FM, EOM_TMP, C_irr, 12, 15, 12, 15, 0, "FDD_FbeMB");
    dpd_buf4_axpy(&FM, &Sijab, -1.0);
    dpd_buf4_close(&FM);
    dpd_buf4_close(&Sijab);

    /* SIjAb += Fbe * CIjAe - FAE * CIjbE */
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, SIjAb_lbl);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, CMnEf_lbl);
    dpd_file2_init(&Fae, CC_OEI, H_IRR, 3, 3, "Fae");
    dpd_contract424(&CMnEf, &Fae, &SIjAb, 3, 1, 0, 1.0, 1.0);
    dpd_file2_close(&Fae);
    dpd_file2_init(&FAE, CC_OEI, H_IRR, 1, 1, "FAE");
    dpd_contract244(&FAE, &CMnEf, &SIjAb, 1, 2, 1, 1.0, 1.0);
    dpd_file2_close(&FAE);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_close(&SIjAb);

#ifdef EOM_DEBUG
    check_sum("Fbe_FDD ",i,C_irr);
#endif

    /* SIJAB -= FMJ * CIMAB - FMI * CJMAB */
    dpd_buf4_init(&FM, EOM_TMP, C_irr, 0, 7, 0, 7, 0, "FDD_FMJM");
    dpd_buf4_init(&CMNEF, EOM_CMNEF, C_irr, 0, 7, 2, 7, 0, CMNEF_lbl);
    dpd_file2_init(&FMI, CC_OEI, H_IRR, 0, 0, "FMI");
    dpd_contract424(&CMNEF, &FMI, &FM, 1, 0, 1, 1.0, 0.0);
    dpd_file2_close(&FMI);
    dpd_buf4_close(&CMNEF);
    dpd_buf4_sort(&FM, EOM_TMP, qprs, 0, 7, "FDD_FMJP");
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 0, 7, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_axpy(&FM, &SIJAB, -1.0);
    dpd_buf4_close(&FM);
    dpd_buf4_init(&FP, EOM_TMP, C_irr, 0, 7, 0, 7, 0, "FDD_FMJP");
    dpd_buf4_axpy(&FP, &SIJAB, 1.0);
    dpd_buf4_close(&FP);
    dpd_buf4_close(&SIJAB);

    /* Sijab -= Fmj * Cimab - Fmi * Cjmab */
    dpd_buf4_init(&FM, EOM_TMP, C_irr, 10, 17, 10, 17, 0, "FDD_FmjMB");
    dpd_buf4_init(&Cmnef, EOM_Cmnef, C_irr, 10, 17, 12, 17, 0, Cmnef_lbl);
    dpd_file2_init(&Fmi, CC_OEI, H_IRR, 2, 2, "Fmi");
    dpd_contract424(&Cmnef, &Fmi, &FM, 1, 0, 1, 1.0, 0.0);
    dpd_file2_close(&Fmi);
    dpd_buf4_close(&Cmnef);
    dpd_buf4_sort(&FM, EOM_TMP, qprs, 10, 17, "FDD_FmjPB");
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 10, 17, 12, 17, 0, Sijab_lbl);
    dpd_buf4_axpy(&FM, &Sijab, -1.0);
    dpd_buf4_close(&FM);
    dpd_buf4_init(&FP, EOM_TMP, C_irr, 10, 17, 10, 17, 0, "FDD_FmjPB");
    dpd_buf4_axpy(&FP, &Sijab, 1.0);
    dpd_buf4_close(&FP);
    dpd_buf4_close(&Sijab);

    /* SIjAb -= Fmj * CImAb - FMI * CjMAb */
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, SIjAb_lbl);
    dpd_buf4_init(&CMnEf, EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, CMnEf_lbl);
    dpd_file2_init(&Fmi, CC_OEI, H_IRR, 2, 2, "Fmi");
    dpd_contract424(&CMnEf, &Fmi, &SIjAb, 1, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Fmi);
    dpd_file2_init(&FMI, CC_OEI, H_IRR, 0, 0, "FMI");
    dpd_contract244(&FMI, &CMnEf, &SIjAb, 0, 0, 0, -1.0, 1.0);
    dpd_file2_close(&FMI);
    dpd_buf4_close(&CMnEf);
    dpd_buf4_close(&SIjAb);

#ifdef EOM_DEBUG
    check_sum("Fmj_DD",i,C_irr);
#endif
  }
  return;
}

}} // namespace psi::cceom
