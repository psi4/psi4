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
   from -0.5*P(ab)*Wnmfe*Cmnea*tijfb and +0.5*Wnmfe*Cimfe*tjnab to a
   Sigma vector stored at Sigma plus 'i' */

void WmnefDD(int i, int C_irr) {
  dpdbuf4 C2, T2, S2, D;
  dpdfile2 X;
  dpdbuf4 SIJAB, Sijab, SIjAb, B;
  dpdbuf4 CMNEF, Cmnef, CMnEf, F, tau, Z;
  char CMNEF_lbl[32], Cmnef_lbl[32], CMnEf_lbl[32];
  char SIJAB_lbl[32], Sijab_lbl[32], SIjAb_lbl[32];
  int l,I,a,f,h,nC_irrs,*occpi,*virtpi,*openpi;

  nC_irrs = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi; openpi = moinfo.openpi;

  if (params.eom_ref == 0) { /* RHF */
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);
    /* XAF = CMNAE * WMNFE + CMnAe * WMnFe */
    global_dpd_->file2_init(&X, PSIF_EOM_TMP, C_irr, 1, 1, "XFA");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "2CMnEf - CMnfE");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract442(&D, &CMnEf, &X, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&CMnEf);

    /* SIjAb += -XFA * TIjFb - TIjAf * Xfb */
    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WmnefDD Z(Ij,Ab)");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->contract424(&T2, &X, &Z, 3, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&X);
    global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP, qpsr, 0, 5, "WmnefDD Z(jI,bA)");
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    global_dpd_->buf4_axpy(&Z, &SIjAb, -1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WmnefDD Z(jI,bA)");
    global_dpd_->buf4_axpy(&Z, &SIjAb, -1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&SIjAb);

#ifdef EOM_DEBUG
    check_sum("WmnefDD XAF",i,C_irr);
#endif

    /* XLI = WLMEF * CIMEF + WLmEf * CImEf */
    global_dpd_->file2_init(&X, PSIF_EOM_TMP, C_irr, 0, 0, "XLI");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "2CMnEf - CMnfE");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract442(&D, &CMnEf, &X, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&CMnEf);

    /* SIjAb += -XLI * TLjAb - Xli * TIlAb */
    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WmnefDD Z(Ij,Ab)");

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->contract244(&X, &T2, &Z, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&X);

    global_dpd_->buf4_sort(&Z, PSIF_EOM_TMP, qpsr, 0, 5, "WmnefDD Z(jI,bA)");

    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    global_dpd_->buf4_axpy(&Z, &SIjAb, -1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_init(&Z, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "WmnefDD Z(jI,bA)");
    global_dpd_->buf4_axpy(&Z, &SIjAb, -1.0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&SIjAb);
  }

  else if (params.eom_ref == 1) { /* ROHF */
    sprintf(CMNEF_lbl, "%s %d", "CMNEF", i);
    sprintf(Cmnef_lbl, "%s %d", "Cmnef", i);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* XAF = CMNAE * WMNFE + CMnAe * WMnFe */
    /* SIJAB -= P(ab) XAF * TIJFB */
    global_dpd_->file2_init(&X, PSIF_EOM_TMP, C_irr, 1, 1, "XFA");
    global_dpd_->file2_scm(&X, 0.0);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 2, 5, 2, 7, 0, CMNEF_lbl);
    global_dpd_->contract442(&D, &CMNEF, &X, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&CMNEF);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    global_dpd_->contract442(&D, &CMnEf, &X, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&S2, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "SIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, H_IRR, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->contract244(&X, &T2, &S2, 0, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&X);
    global_dpd_->buf4_sort(&S2, PSIF_EOM_TMP, pqsr, 2, 5, "SIJBA");
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 2, 5, 2, 7, 0, SIJAB_lbl);
    global_dpd_->buf4_axpy(&S2, &SIJAB, -1.0);
    global_dpd_->buf4_close(&S2);
    global_dpd_->buf4_init(&S2, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "SIJBA");
    global_dpd_->buf4_axpy(&S2, &SIJAB, 1.0);
    global_dpd_->buf4_close(&S2);
    global_dpd_->buf4_close(&SIJAB);

    /* Xaf = Cmnae * Wmnfe + CmNaE * WmNfE */
    /* Sijab -= P(ab) Xfa * Tijfb */
    global_dpd_->file2_init(&X, PSIF_EOM_TMP, C_irr, 1, 1, "Xfa");
    global_dpd_->file2_scm(&X, 0.0);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 2, 5, 2, 7, 0, Cmnef_lbl);
    global_dpd_->contract442(&D, &Cmnef, &X, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&Cmnef);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CmNeF");
    global_dpd_->contract442(&D, &CMnEf, &X, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&S2, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "Sijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, H_IRR, 2, 5, 2, 7, 0, "tijab");
    global_dpd_->contract244(&X, &T2, &S2, 0, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&X);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_sort(&S2, PSIF_EOM_TMP, pqsr, 2, 5, "Sijba");
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 2, 5, 2, 7, 0, Sijab_lbl);
    global_dpd_->buf4_axpy(&S2, &Sijab, -1.0);
    global_dpd_->buf4_close(&S2);
    global_dpd_->buf4_init(&S2, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "Sijba");
    global_dpd_->buf4_axpy(&S2, &Sijab, 1.0);
    global_dpd_->buf4_close(&S2);
    global_dpd_->buf4_close(&Sijab);

    /* SIjAb += -XFA * TIjFb - TIjAf * Xfb */
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    global_dpd_->file2_init(&X, PSIF_EOM_TMP, C_irr, 1, 1, "XFA");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->contract244(&X, &T2, &SIjAb, 0, 2, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&X);
    global_dpd_->file2_init(&X, PSIF_EOM_TMP, C_irr, 1, 1, "Xfa");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->contract424(&T2, &X, &SIjAb, 3, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&X);
    global_dpd_->buf4_close(&SIjAb);

#ifdef EOM_DEBUG
    check_sum("WmnefDD XAF",i,C_irr);
#endif

    /* XLI = WLMEF * CIMEF + WLmEf * CImEf */
    /* SIJAB += P(IJ) XLI * TLJAB */
    global_dpd_->file2_init(&X, PSIF_EOM_TMP, C_irr, 0, 0, "XLI");
    global_dpd_->file2_scm(&X, 0.0);
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 0, 7, 2, 7, 0, CMNEF_lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    global_dpd_->contract442(&D, &CMNEF, &X, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&CMNEF);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, CMnEf_lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract442(&D, &CMnEf, &X, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&CMnEf);

    /*
       dpd_file2_mat_init(&X);
       dpd_file2_mat_rd(&X);
       for(h=0; h < nC_irrs; h++) {
       for(l=0; l<occpi[h]; l++)
       for(I=(occpi[h]-openpi[h]); I<occpi[h]; I++)
       X.matrix[h][l][I] = 0.0;
       for(I=0; I<occpi[h]; I++)
       for(l=(occpi[h]-openpi[h]); l<occpi[h]; l++)
       X.matrix[h][l][I] = 0.0;
       }
       dpd_file2_mat_wrt(&X);
       dpd_file2_mat_close(&X);
     */

    global_dpd_->buf4_init(&S2, PSIF_EOM_TMP, C_irr, 0, 7, 0, 7, 0, "SIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, H_IRR, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->contract244(&X, &T2, &S2, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&X);
    global_dpd_->buf4_sort(&S2, PSIF_EOM_TMP, qprs, 0, 7, "SJIAB");
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 0, 7, 2, 7, 0, SIJAB_lbl);
    global_dpd_->buf4_axpy(&S2, &SIJAB, -1.0);
    global_dpd_->buf4_close(&S2);
    global_dpd_->buf4_init(&S2, PSIF_EOM_TMP, C_irr, 0, 7, 0, 7, 0, "SJIAB");
    global_dpd_->buf4_axpy(&S2, &SIJAB, 1.0);
    global_dpd_->buf4_close(&S2);
    global_dpd_->buf4_close(&SIJAB);


    /* Xli = Wlmef * Cimef + WlMeF * CiMeF */
    /* Sijab += P(ij) Xli * Tljab */
    global_dpd_->file2_init(&X, PSIF_EOM_TMP, C_irr, 0, 0, "Xli");
    global_dpd_->file2_scm(&X, 0.0);
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 0, 7, 2, 7, 0, Cmnef_lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    global_dpd_->contract442(&D, &Cmnef, &X, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Cmnef);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_TMP, C_irr, 0, 5, 0, 5, 0, "CmNeF");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract442(&D, &CMnEf, &X, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_init(&S2, PSIF_EOM_TMP, C_irr, 0, 7, 0, 7, 0, "Sijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, H_IRR, 0, 7, 2, 7, 0, "tijab");
    global_dpd_->contract244(&X, &T2, &S2, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&X);
    global_dpd_->buf4_sort(&S2, PSIF_EOM_TMP, qprs, 0, 7, "Sjiab");
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 0, 7, 2, 7, 0, Sijab_lbl);
    global_dpd_->buf4_axpy(&S2, &Sijab, -1.0);
    global_dpd_->buf4_close(&S2);
    global_dpd_->buf4_init(&S2, PSIF_EOM_TMP, C_irr, 0, 7, 0, 7, 0, "Sjiab");
    global_dpd_->buf4_axpy(&S2, &Sijab, 1.0);
    global_dpd_->buf4_close(&S2);
    global_dpd_->buf4_close(&Sijab);

    /* SIjAb += -XLI * TLjAb - Xli * TIlAb */
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, SIjAb_lbl);
    global_dpd_->file2_init(&X, PSIF_EOM_TMP, C_irr, 0, 0, "XLI");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, H_IRR, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->contract244(&X, &T2, &SIjAb, 0, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&X);
    global_dpd_->file2_init(&X, PSIF_EOM_TMP, C_irr, 0, 0, "Xli");
    global_dpd_->contract424(&T2, &X, &SIjAb, 1, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&X);
    global_dpd_->buf4_close(&SIjAb);
  }

  else if (params.eom_ref == 2) {
    sprintf(CMNEF_lbl, "%s %d", "CMNEF", i);
    sprintf(Cmnef_lbl, "%s %d", "Cmnef", i);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    sprintf(SIJAB_lbl, "%s %d", "SIJAB", i);
    sprintf(Sijab_lbl, "%s %d", "Sijab", i);
    sprintf(SIjAb_lbl, "%s %d", "SIjAb", i);

    /* XAF = CMNAE * WMNFE + CMnAe * WMnFe */
    /* SIJAB -= P(ab) XAF * TIJFB */
    global_dpd_->file2_init(&X, PSIF_EOM_TMP, C_irr, 1, 1, "XFA");
    global_dpd_->file2_scm(&X, 0.0);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 2, 5, 2, 7, 0, CMNEF_lbl);
    global_dpd_->contract442(&D, &CMNEF, &X, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&CMNEF);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, CMnEf_lbl);
    global_dpd_->contract442(&D, &CMnEf, &X, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&S2, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "SIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, H_IRR, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->contract244(&X, &T2, &S2, 0, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&X);
    global_dpd_->buf4_sort(&S2, PSIF_EOM_TMP, pqsr, 2, 5, "SIJBA");
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 2, 5, 2, 7, 0, SIJAB_lbl);
    global_dpd_->buf4_axpy(&S2, &SIJAB, -1.0);
    global_dpd_->buf4_close(&S2);
    global_dpd_->buf4_init(&S2, PSIF_EOM_TMP, C_irr, 2, 5, 2, 5, 0, "SIJBA");
    global_dpd_->buf4_axpy(&S2, &SIJAB, 1.0);
    global_dpd_->buf4_close(&S2);
    global_dpd_->buf4_close(&SIJAB);

    /* Xaf = Cmnae * Wmnfe + CmNaE * WmNfE */
    /* Sijab -= P(ab) Xfa * Tijfb */
    global_dpd_->file2_init(&X, PSIF_EOM_TMP, C_irr, 3, 3, "Xfa");
    global_dpd_->file2_scm(&X, 0.0);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 12, 15, 12, 17, 0, Cmnef_lbl);
    global_dpd_->contract442(&D, &Cmnef, &X, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&Cmnef);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 23, 29, 23, 29, 0, "D <iJ|aB>");
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_TMP, C_irr, 23, 29, 23, 29, 0, "CmNeF");
    global_dpd_->contract442(&D, &CMnEf, &X, 2, 2, 1.0, 1.0);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&S2, PSIF_EOM_TMP, C_irr, 12, 15, 12, 15, 0, "Sijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, H_IRR, 12, 15, 12, 17, 0, "tijab");
    global_dpd_->contract244(&X, &T2, &S2, 0, 2, 1, 1.0, 0.0);
    global_dpd_->file2_close(&X);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_sort(&S2, PSIF_EOM_TMP, pqsr, 12, 15, "Sijba");
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 12, 15, 12, 17, 0, Sijab_lbl);
    global_dpd_->buf4_axpy(&S2, &Sijab, -1.0);
    global_dpd_->buf4_close(&S2);
    global_dpd_->buf4_init(&S2, PSIF_EOM_TMP, C_irr, 12, 15, 12, 15, 0, "Sijba");
    global_dpd_->buf4_axpy(&S2, &Sijab, 1.0);
    global_dpd_->buf4_close(&S2);
    global_dpd_->buf4_close(&Sijab);

    /* SIjAb += -XFA * TIjFb - TIjAf * Xfb */
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, SIjAb_lbl);
    global_dpd_->file2_init(&X, PSIF_EOM_TMP, C_irr, 1, 1, "XFA");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, H_IRR, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->contract244(&X, &T2, &SIjAb, 0, 2, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&X);
    global_dpd_->file2_init(&X, PSIF_EOM_TMP, C_irr, 3, 3, "Xfa");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, H_IRR, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->contract424(&T2, &X, &SIjAb, 3, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&X);
    global_dpd_->buf4_close(&SIjAb);

#ifdef EOM_DEBUG
    check_sum("WmnefDD XAF",i,C_irr);
#endif

    /* XLI = WLMEF * CIMEF + WLmEf * CImEf */
    /* SIJAB += P(IJ) XLI * TLJAB */
    global_dpd_->file2_init(&X, PSIF_EOM_TMP, C_irr, 0, 0, "XLI");
    global_dpd_->file2_scm(&X, 0.0);
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 0, 7, 2, 7, 0, CMNEF_lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    global_dpd_->contract442(&D, &CMNEF, &X, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&CMNEF);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, CMnEf_lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->contract442(&D, &CMnEf, &X, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&CMnEf);

    global_dpd_->buf4_init(&S2, PSIF_EOM_TMP, C_irr, 0, 7, 0, 7, 0, "SIJAB");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, H_IRR, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->contract244(&X, &T2, &S2, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&X);
    global_dpd_->buf4_sort(&S2, PSIF_EOM_TMP, qprs, 0, 7, "SJIAB");
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 0, 7, 2, 7, 0, SIJAB_lbl);
    global_dpd_->buf4_axpy(&S2, &SIJAB, -1.0);
    global_dpd_->buf4_close(&S2);
    global_dpd_->buf4_init(&S2, PSIF_EOM_TMP, C_irr, 0, 7, 0, 7, 0, "SJIAB");
    global_dpd_->buf4_axpy(&S2, &SIJAB, 1.0);
    global_dpd_->buf4_close(&S2);
    global_dpd_->buf4_close(&SIJAB);


    /* Xli = Wlmef * Cimef + WlMeF * CiMeF */
    /* Sijab += P(ij) Xli * Tljab */
    global_dpd_->file2_init(&X, PSIF_EOM_TMP, C_irr, 2, 2, "Xli");
    global_dpd_->file2_scm(&X, 0.0);
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 10, 17, 12, 17, 0, Cmnef_lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    global_dpd_->contract442(&D, &Cmnef, &X, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Cmnef);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_TMP, C_irr, 23, 29, 23, 29, 0, "CmNeF");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, H_IRR, 23, 29, 23, 29, 0, "D <iJ|aB>");
    global_dpd_->contract442(&D, &CMnEf, &X, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_init(&S2, PSIF_EOM_TMP, C_irr, 10, 17, 10, 17, 0, "Sijab");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, H_IRR, 10, 17, 12, 17, 0, "tijab");
    global_dpd_->contract244(&X, &T2, &S2, 0, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&X);
    global_dpd_->buf4_sort(&S2, PSIF_EOM_TMP, qprs, 10, 17, "Sjiab");
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 10, 17, 12, 17, 0, Sijab_lbl);
    global_dpd_->buf4_axpy(&S2, &Sijab, -1.0);
    global_dpd_->buf4_close(&S2);
    global_dpd_->buf4_init(&S2, PSIF_EOM_TMP, C_irr, 10, 17, 10, 17, 0, "Sjiab");
    global_dpd_->buf4_axpy(&S2, &Sijab, 1.0);
    global_dpd_->buf4_close(&S2);
    global_dpd_->buf4_close(&Sijab);

    /* SIjAb += -XLI * TLjAb - Xli * TIlAb */
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, SIjAb_lbl);
    global_dpd_->file2_init(&X, PSIF_EOM_TMP, C_irr, 0, 0, "XLI");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, H_IRR, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->contract244(&X, &T2, &SIjAb, 0, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&X);
    global_dpd_->file2_init(&X, PSIF_EOM_TMP, C_irr, 2, 2, "Xli");
    global_dpd_->contract424(&T2, &X, &SIjAb, 1, 0, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&X);
    global_dpd_->buf4_close(&SIjAb);
  }

#ifdef EOM_DEBUG
  check_sum("WmnefDD XLI",i,C_irr);
#endif
  return;
}

}} // namespace psi::cceom
