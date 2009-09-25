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
    dpd_file2_init(&XNJ, EOM_TMP, C_irr, 0, 0, "XNJ");
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_dot23(&CME, &W, &XNJ, 0, 0, 1.0, 0.0);
    dpd_file2_close(&CME);
    dpd_buf4_close(&W);

    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WnmjeDS Z(Ij,Ab)");
    dpd_buf4_init(&TIjAb, CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract244(&XNJ, &TIjAb, &Z, 0, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&TIjAb);
    dpd_file2_close(&XNJ);
    dpd_buf4_sort(&Z, EOM_TMP, qpsr, 0, 5, "WnmjeDS Z(jI,bA)"); 

    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_axpy(&Z, &SIjAb, -1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_init(&Z, EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WnmjeDS Z(jI,bA)");
    dpd_buf4_axpy(&Z, &SIjAb, -1.0);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&SIjAb);
  }

  else if (params.eom_ref == 1) { /* ROHF */
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(Cme_lbl, "%s %d", "Cme", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* Form XNJ intermediates */
    /* XNJ = CME * WNMJE + Cme * WNmJe */
    dpd_file2_init(&XNJ, EOM_TMP, C_irr, 0, 0, "XNJ");
    dpd_file2_scm(&XNJ, 0.0);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_buf4_init(&WMNIE, CC_HBAR, H_IRR, 0, 11, 2, 11, 0, "WMNIE (M>N,EI)");
    dpd_dot23(&CME, &WMNIE, &XNJ, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&WMNIE);
    dpd_file2_close(&CME);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);
    dpd_buf4_init(&WMnIe, CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
    dpd_dot23(&Cme, &WMnIe, &XNJ, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&WMnIe);
    dpd_file2_close(&Cme);
    dpd_file2_close(&XNJ);


    /* Xnj = Cme * Wnmje + CME * WnMjE */
    dpd_file2_init(&Xnj, EOM_TMP, C_irr, 0, 0, "Xnj");
    dpd_file2_scm(&Xnj, 0.0);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, Cme_lbl);
    dpd_buf4_init(&Wmnie, CC_HBAR, H_IRR, 0, 11, 2, 11, 0, "Wmnie (m>n,ei)");
    dpd_dot23(&Cme, &Wmnie, &Xnj, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&Wmnie);
    dpd_file2_close(&Cme);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_buf4_init(&WmNiE, CC_HBAR, H_IRR, 0, 11, 0, 11, 0, "WmNiE (mN,Ei)");
    dpd_dot23(&CME, &WmNiE, &Xnj, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&WmNiE);
    dpd_file2_close(&CME);
    dpd_file2_close(&Xnj);

    /* SIJAB -= XNJ * TINAB + XNI * TJNAB */
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WnmjeDS_M");
    dpd_buf4_init(&TIJAB, CC_TAMPS, H_IRR, 0, 7, 2, 7, 0, "tIJAB");
    dpd_file2_init(&XNJ, EOM_TMP, C_irr, 0, 0, "XNJ");
    dpd_contract424(&TIJAB, &XNJ, &WM, 1, 0, 1, 1.0, 0.0);
    dpd_file2_close(&XNJ);
    dpd_buf4_close(&TIJAB);
    dpd_buf4_sort(&WM, EOM_TMP, qprs, 0, 7, "WnmjeDS_P");
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 0, 7, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_axpy(&WM, &SIJAB, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WnmjeDS_P");
    dpd_buf4_axpy(&WP, &SIJAB, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_close(&SIJAB);

    /* Sijab -= Xnj * Tinab + Xni * Tjnab */
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WnmjeDS_M");
    dpd_buf4_init(&Tijab, CC_TAMPS, H_IRR, 0, 7, 2, 7, 0, "tijab");
    dpd_file2_init(&Xnj, EOM_TMP, C_irr, 0, 0, "Xnj");
    dpd_contract424(&Tijab, &Xnj, &WM, 1, 0, 1, 1.0, 0.0);
    dpd_file2_close(&Xnj);
    dpd_buf4_close(&Tijab);
    dpd_buf4_sort(&WM, EOM_TMP, qprs, 0, 7, "WnmjeDS_P");
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 0, 7, 2, 7, 0, Sijab_lbl);
    dpd_buf4_axpy(&WM, &Sijab, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WnmjeDS_P");
    dpd_buf4_axpy(&WP, &Sijab, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_close(&Sijab);

    /* SIjAb -= Xnj * tInAb + XNI * TjNAb */
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    dpd_buf4_init(&TIjAb, CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tIjAb");
    dpd_file2_init(&Xnj, EOM_TMP, C_irr, 0, 0, "Xnj");
    dpd_contract424(&TIjAb, &Xnj, &SIjAb, 1, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Xnj);
    dpd_file2_init(&XNJ, EOM_TMP, C_irr, 0, 0, "XNJ");
    dpd_contract244(&XNJ, &TIjAb, &SIjAb, 0, 0, 0, -1.0, 1.0);
    dpd_file2_close(&XNJ);
    dpd_buf4_close(&TIjAb);
    dpd_buf4_close(&SIjAb);
  }

  else if (params.eom_ref == 2) {
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(Cme_lbl, "%s %d", "Cme", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* Form XNJ intermediates */
    /* XNJ = CME * WNMJE + Cme * WNmJe */
    dpd_file2_init(&XNJ, EOM_TMP, C_irr, 0, 0, "XNJ");
    dpd_file2_scm(&XNJ, 0.0);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_buf4_init(&WMNIE, CC_HBAR, H_IRR, 0, 21, 2, 21, 0, "WMNIE (M>N,EI)");
    dpd_dot23(&CME, &WMNIE, &XNJ, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&WMNIE);
    dpd_file2_close(&CME);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, Cme_lbl);
    dpd_buf4_init(&WMnIe, CC_HBAR, H_IRR, 22, 25, 22, 25, 0, "WMnIe (Mn,eI)");
    dpd_dot23(&Cme, &WMnIe, &XNJ, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&WMnIe);
    dpd_file2_close(&Cme);
/*
    tval = dpd_file2_dot_self(&XNJ);
    fprintf(outfile,"XNJ self dot %15.10lf\n",tval);
*/
    dpd_file2_close(&XNJ);


    /* Xnj = Cme * Wnmje + CME * WnMjE */
    dpd_file2_init(&Xnj, EOM_TMP, C_irr, 2, 2, "Xnj");
    dpd_file2_scm(&Xnj, 0.0);
    dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, Cme_lbl);
    dpd_buf4_init(&Wmnie, CC_HBAR, H_IRR, 10, 31, 12, 31, 0, "Wmnie (m>n,ei)");
    dpd_dot23(&Cme, &Wmnie, &Xnj, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&Wmnie);
    dpd_file2_close(&Cme);
    dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, CME_lbl);
    dpd_buf4_init(&WmNiE, CC_HBAR, H_IRR, 23, 26, 23, 26, 0, "WmNiE (mN,Ei)");
    dpd_dot23(&CME, &WmNiE, &Xnj, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&WmNiE);
    dpd_file2_close(&CME);

/*
    tval = dpd_file2_dot_self(&Xnj);
    fprintf(outfile,"Xnj self dot %15.10lf\n",tval);
*/
    dpd_file2_close(&Xnj);

    /* SIJAB -= XNJ * TINAB + XNI * TJNAB */
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WnmjeDS_M");
    dpd_buf4_init(&TIJAB, CC_TAMPS, H_IRR, 0, 7, 2, 7, 0, "tIJAB");
    dpd_file2_init(&XNJ, EOM_TMP, C_irr, 0, 0, "XNJ");
    dpd_contract424(&TIJAB, &XNJ, &WM, 1, 0, 1, 1.0, 0.0);
    dpd_file2_close(&XNJ);
    dpd_buf4_close(&TIJAB);
    dpd_buf4_sort(&WM, EOM_TMP, qprs, 0, 7, "WnmjeDS_P");
    dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 0, 7, 2, 7, 0, SIJAB_lbl);
    dpd_buf4_axpy(&WM, &SIJAB, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 0, 7, 0, 7, 0, "WnmjeDS_P");
    dpd_buf4_axpy(&WP, &SIJAB, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_close(&SIJAB);

    /* Sijab -= Xnj * Tinab + Xni * Tjnab */
    dpd_buf4_init(&WM, EOM_TMP, C_irr, 10, 17, 10, 17, 0, "WnmjeDS_MB");
    dpd_buf4_init(&Tijab, CC_TAMPS, H_IRR, 10, 17, 12, 17, 0, "tijab");
    dpd_file2_init(&Xnj, EOM_TMP, C_irr, 2, 2, "Xnj");
    dpd_contract424(&Tijab, &Xnj, &WM, 1, 0, 1, 1.0, 0.0);
    dpd_file2_close(&Xnj);
    dpd_buf4_close(&Tijab);
    dpd_buf4_sort(&WM, EOM_TMP, qprs, 10, 17, "WnmjeDS_PB");
    dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 10, 17, 12, 17, 0, Sijab_lbl);
    dpd_buf4_axpy(&WM, &Sijab, -1.0);
    dpd_buf4_close(&WM);
    dpd_buf4_init(&WP, EOM_TMP, C_irr, 10, 17, 10, 17, 0, "WnmjeDS_PB");
    dpd_buf4_axpy(&WP, &Sijab, 1.0);
    dpd_buf4_close(&WP);
    dpd_buf4_close(&Sijab);

    /* SIjAb -= Xnj * tInAb + XNI * TjNAb */
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, SIjAb_lbl);
    dpd_buf4_init(&TIjAb, CC_TAMPS, H_IRR, 22, 28, 22, 28, 0, "tIjAb");
    dpd_file2_init(&Xnj, EOM_TMP, C_irr, 2, 2, "Xnj");
    dpd_contract424(&TIjAb, &Xnj, &SIjAb, 1, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Xnj);
    dpd_file2_init(&XNJ, EOM_TMP, C_irr, 0, 0, "XNJ");
    dpd_contract244(&XNJ, &TIjAb, &SIjAb, 0, 0, 0, -1.0, 1.0);
    dpd_file2_close(&XNJ);
    dpd_buf4_close(&TIjAb);
    dpd_buf4_close(&SIjAb);
  }

#ifdef EOM_DEBUG
  check_sum("WnmjeDS",i,C_irr);
#endif
  return;
}

}} // namespace psi::cceom
