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
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/* This function computes the extra contributions to sigma_1 and sigma_2
  for EOM_CC3 computations that are not normally present in a EOM_CCSD
  calculation */

/* The additional terms are:
 * Term 1:
 * <S| H    <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_1
 * <D| Hhat <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_2
 * Term 2:
 * <S| H    <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_1
 * <D| Hhat <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_2
 * Term 3:
 * <D| H'   <T| (Uhat T2)c   |0> |T> / (-wt) -> sigma_2
 *
 *  See Eqn. (83) in JCP, 103, 7429, 1995
 *  All three terms can be evaluated by the same functions in 
 *  cc3_sigma_UHF() and cc3_sigma_RHF() in libdpd given
 *  different matrix elements.
 *
 * - RHF case added by generalizing ccenergy code RAK 2006
 *
 * */

void sigmaCC3(int i, int C_irr, double omega) {
  dpdfile2 CME, Cme, SIA, Sia, FME, Fme;
  dpdbuf4 CMNEF, CMnEf, Cmnef, CmNeF, SIJAB, Sijab, SIjAb;
  dpdbuf4 WAMEF, WMNIE, WABEI, WMBIJ, DIJAB_anti, TIJAB;
  dpdbuf4 Wamef, Wmnie, Wabei, Wmbij, Dijab_anti, Tijab;
  dpdbuf4 WAmEf, WMnIe, WAbEi, WMbIj, DIjAb, TIjAb, tIjAb;
  dpdbuf4 WmAEf, WaMeF, WmNiE, WaBeI, WmBiJ, DiJaB, TiJaB, Dints;
  char lbl[32];

  if (params.eom_ref == 0) { /* RHF */
    sprintf(lbl, "%s %d", "SIA", i);
    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "SIjAb", i);
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);

    /*** alpha-alpha-beta term 1 ***/ 
    /* quantities to compute X3 */
    sprintf(lbl, "%s %d", "CMnEf", i);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&WAbEi, PSIF_CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAbEi (iE,bA)");
    global_dpd_->buf4_init(&WMbIj, PSIF_CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMbIj (Ij,Mb)");
    /* quantities to compute sigma */
    global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->buf4_init(&WmAEf, PSIF_CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAmEf (mA,Ef)");
    global_dpd_->buf4_init(&WMnIe, PSIF_CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMnIe (Mn,Ie)");
  
         /* * <S| H    <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_1
            * <D| Hhat <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_2 */

    if (params.t3_Ws_incore)
      global_dpd_->cc3_sigma_RHF_ic(&CMnEf, &WAbEi, &WMbIj, 1,  &Dints, &SIA, 
        1, &FME, &WmAEf, &WMnIe, &SIjAb, moinfo.occpi, moinfo.occ_off,
        moinfo.virtpi, moinfo.vir_off, omega, "outfile", params.nthreads,
        params.newtrips);
    else
      global_dpd_->cc3_sigma_RHF(&CMnEf, &WAbEi, &WMbIj, 1,  &Dints, &SIA, 
        1, &FME, &WmAEf, &WMnIe, &SIjAb, moinfo.occpi, moinfo.occ_off,
        moinfo.virtpi, moinfo.vir_off, omega, "outfile", params.newtrips);
  
    global_dpd_->buf4_close(&CMnEf);
    global_dpd_->buf4_close(&WAbEi);
    global_dpd_->buf4_close(&WMbIj);
    global_dpd_->buf4_close(&Dints);
    global_dpd_->file2_close(&FME);
    global_dpd_->buf4_close(&WmAEf);
    global_dpd_->buf4_close(&WMnIe);

#ifdef EOM_DEBUG
    dpd_file2_close(&SIA);
    dpd_buf4_close(&SIjAb);
    check_sum("<Psi|Hhat<T|(Uhat C2)c|0>|T>/(w-wt)", i, C_irr);
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "SIjAb", i);
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
#endif

    /* do alpha-alpha-beta term 2 */
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&WAbEi, PSIF_CC3_HC1ET1, C_irr, 10, 5, 10, 5, 0, "Ht_WAbEi (iE,bA)");
    global_dpd_->buf4_init(&WMbIj, PSIF_CC3_HC1ET1, C_irr, 0, 10, 0, 10, 0, "Ht_WMbIj (Ij,Mb)");
  
    global_dpd_->buf4_init(&Dints, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->buf4_init(&WmAEf, PSIF_CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAmEf (mA,Ef)");
    global_dpd_->buf4_init(&WMnIe, PSIF_CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMnIe (Mn,Ie)");
  
           /* * <S| H    <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_1
              * <D| Hhat <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_2 */
  
    if (params.t3_Ws_incore)
      global_dpd_->cc3_sigma_RHF_ic(&tIjAb, &WAbEi, &WMbIj, 1,  &Dints, &SIA, 
        1, &FME, &WmAEf, &WMnIe, &SIjAb, moinfo.occpi, moinfo.occ_off,
        moinfo.virtpi, moinfo.vir_off, omega, "outfile", params.nthreads,
        params.newtrips);
    else
      global_dpd_->cc3_sigma_RHF(&tIjAb, &WAbEi, &WMbIj, 1,  &Dints, &SIA,
         1, &FME, &WmAEf, &WMnIe, &SIjAb, moinfo.occpi, moinfo.occ_off,
         moinfo.virtpi, moinfo.vir_off, omega, "outfile", params.newtrips);
  
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_close(&WAbEi);
    global_dpd_->buf4_close(&WMbIj);
    global_dpd_->buf4_close(&Dints);
    global_dpd_->file2_close(&FME);
    global_dpd_->buf4_close(&WmAEf);
    global_dpd_->buf4_close(&WMnIe);
  
#ifdef EOM_DEBUG
    dpd_file2_close(&SIA);
    dpd_buf4_close(&SIjAb);
    check_sum("<Psi|Hhat<T|(Utilde T2)c|0>|T>/(w-wt)", i, C_irr);
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "SIjAb", i);
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
#endif

    /* alpha-alpha-beta term 3 */
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&WAbEi, PSIF_CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAbEi (iE,bA)");
    global_dpd_->buf4_init(&WMbIj, PSIF_CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMbIj (Ij,Mb)");
  
    global_dpd_->file2_init(&FME, PSIF_CC3_HC1, C_irr, 0, 1, "HC1 FME");
    global_dpd_->buf4_init(&WmAEf, PSIF_CC3_HC1, C_irr, 10, 5, 10, 5, 0, "HC1 WAmEf (mA,Ef)");
    global_dpd_->buf4_init(&WMnIe, PSIF_CC3_HC1, C_irr, 0, 10, 0, 10, 0, "HC1 WMnIe (Mn,Ie)");
  
           /* <D| H'   <T| (Uhat T2)c   |0> |T> / (-wt) -> sigma_2 */
  
    if (params.t3_Ws_incore)
      global_dpd_->cc3_sigma_RHF_ic(&tIjAb, &WAbEi, &WMbIj, 0, NULL, NULL, 
        1, &FME, &WmAEf, &WMnIe, &SIjAb, moinfo.occpi, moinfo.occ_off,
        moinfo.virtpi, moinfo.vir_off, 0.0, "outfile", params.nthreads,
        params.newtrips);
    else
    global_dpd_->cc3_sigma_RHF(&tIjAb, &WAbEi, &WMbIj, 0, NULL, NULL,
       1, &FME, &WmAEf, &WMnIe, &SIjAb, moinfo.occpi, moinfo.occ_off,
       moinfo.virtpi, moinfo.vir_off, 0.0, "outfile", params.newtrips);
  
    global_dpd_->buf4_close(&tIjAb); 
    global_dpd_->buf4_close(&WAbEi);
    global_dpd_->buf4_close(&WMbIj);
    global_dpd_->file2_close(&FME);
    global_dpd_->buf4_close(&WmAEf);
    global_dpd_->buf4_close(&WMnIe);

#ifdef EOM_DEBUG
    dpd_file2_close(&SIA);
    dpd_buf4_close(&SIjAb);
    check_sum("<Psi|H'<T|(Uhat T2)c|0>|T>/(w-wt)", i, C_irr);
    sprintf(lbl, "%s %d", "SIA", i);
    dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "SIjAb", i);
    dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 0, 5, 0, 5, 0, lbl);
#endif

    global_dpd_->file2_close(&SIA);
    global_dpd_->buf4_close(&SIjAb);
  }
  else if (params.eom_ref == 1) { /* ROHF */
  }
  else if (params.eom_ref == 2) { /* UHF */
    /* open all sigma (output) files */
    sprintf(lbl, "%s %d", "SIA", i);
    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
    sprintf(lbl, "%s %d", "Sia", i);
    global_dpd_->file2_init(&Sia, PSIF_EOM_Sia, C_irr, 2, 3, lbl);
    sprintf(lbl, "%s %d", "SIJAB", i);
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, C_irr, 0, 5, 2, 7, 0, lbl);
    sprintf(lbl, "%s %d", "Sijab", i);
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, C_irr, 10, 15, 12, 17, 0, lbl);
    sprintf(lbl, "%s %d", "SIjAb", i);
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, lbl);

    /*** alpha-alpha-alpha term 1 */

    sprintf(lbl, "%s %d", "CMNEF", i);
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 0, 5, 2, 7, 0, lbl);
    global_dpd_->buf4_init(&WABEI, PSIF_CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WABEI (IE,B>A)");
    global_dpd_->buf4_init(&WMBIJ, PSIF_CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMBIJ (I>J,MB)");
    global_dpd_->buf4_init(&DIJAB_anti, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->buf4_init(&WAMEF, PSIF_CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WAMEF (MA,F>E)");
    global_dpd_->buf4_init(&WMNIE, PSIF_CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMNIE (M>N,IE)");

         /* * <S| H    <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_1
            * <D| Hhat <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_2 */

    global_dpd_->cc3_sigma_UHF_AAA(&CMNEF, &WABEI, &WMBIJ, 1, &DIJAB_anti, &SIA,
        1, &FME, &WAMEF, &WMNIE, &SIJAB, moinfo.aoccpi, moinfo.aocc_off,
        moinfo.avirtpi, moinfo.avir_off, omega, "outfile");

    global_dpd_->buf4_close(&CMNEF);
    global_dpd_->buf4_close(&WABEI);
    global_dpd_->buf4_close(&WMBIJ);
    global_dpd_->buf4_close(&DIJAB_anti);
    global_dpd_->file2_close(&FME);
    global_dpd_->buf4_close(&WAMEF);
    global_dpd_->buf4_close(&WMNIE);

    /*** beta-beta-beta term 1 */

    sprintf(lbl, "%s %d", "Cmnef", i);
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 10, 15, 12, 17, 0, lbl);
    global_dpd_->buf4_init(&Wabei, PSIF_CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wabei (ie,b>a)");
    global_dpd_->buf4_init(&Wmbij, PSIF_CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmbij (i>j,mb)");
    global_dpd_->buf4_init(&Dijab_anti, PSIF_CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
    global_dpd_->buf4_init(&Wamef, PSIF_CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wamef (ma,f>e)");
    global_dpd_->buf4_init(&Wmnie, PSIF_CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmnie (m>n,ie)");

         /* * <S| H    <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_1
            * <D| Hhat <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_2 */

    global_dpd_->cc3_sigma_UHF_BBB(&Cmnef, &Wabei, &Wmbij, 1, &Dijab_anti, &Sia,
        1, &Fme, &Wamef, &Wmnie, &Sijab, moinfo.boccpi, moinfo.bocc_off,
        moinfo.bvirtpi, moinfo.bvir_off, omega, "outfile");

    global_dpd_->buf4_close(&Cmnef);
    global_dpd_->buf4_close(&Wabei);
    global_dpd_->buf4_close(&Wmbij);
    global_dpd_->buf4_close(&Dijab_anti);
    global_dpd_->file2_close(&Fme);
    global_dpd_->buf4_close(&Wamef);
    global_dpd_->buf4_close(&Wmnie);

    /*** alpha-alpha-beta term 1 */ 

    sprintf(lbl, "%s %d", "CMNEF", i);
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, C_irr, 0, 5, 2, 7, 0, lbl);
    sprintf(lbl, "%s %d", "CMnEf", i);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, lbl);
    global_dpd_->buf4_init(&CmNeF, PSIF_EOM_TMP, C_irr, 23, 29, 23, 29, 0, "CmNeF");

    global_dpd_->buf4_init(&WABEI, PSIF_CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WABEI (IE,B>A)");
    global_dpd_->buf4_init(&WaBeI, PSIF_CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaBeI (Ie,Ba)");
    global_dpd_->buf4_init(&WAbEi, PSIF_CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAbEi (iE,bA)");
    global_dpd_->buf4_init(&WMBIJ, PSIF_CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMBIJ (I>J,MB)");
    global_dpd_->buf4_init(&WMbIj, PSIF_CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMbIj (Ij,Mb)");
    global_dpd_->buf4_init(&WmBiJ, PSIF_CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmBiJ (iJ,mB)");

    global_dpd_->buf4_init(&DIJAB_anti, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
    global_dpd_->buf4_init(&DIjAb, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");

    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
    global_dpd_->buf4_init(&WAMEF, PSIF_CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WAMEF (MA,F>E)");
    global_dpd_->buf4_init(&WaMeF, PSIF_CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaMeF (Ma,Fe)");
    global_dpd_->buf4_init(&WAmEf, PSIF_CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAmEf (mA,fE)");
    global_dpd_->buf4_init(&WMNIE, PSIF_CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMNIE (M>N,IE)");
    global_dpd_->buf4_init(&WMnIe, PSIF_CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMnIe (Mn,Ie)");
    global_dpd_->buf4_init(&WmNiE, PSIF_CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmNiE (mN,iE)");

         /* * <S| H    <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_1
            * <D| Hhat <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_2 */

    global_dpd_->cc3_sigma_UHF_AAB(&CMNEF, &CMnEf, &CmNeF, &WABEI, &WaBeI, &WAbEi,
       &WMBIJ, &WMbIj, &WmBiJ, 1,  &DIJAB_anti, &DIjAb, &SIA, &Sia,
       1, &FME, &Fme, &WAMEF, &WaMeF, &WAmEf, &WMNIE, &WMnIe, &WmNiE,
       &SIJAB, &SIjAb, moinfo.aoccpi, moinfo.aocc_off, moinfo.boccpi,
       moinfo.bocc_off, moinfo.avirtpi, moinfo.avir_off, moinfo.bvirtpi,
       moinfo.bvir_off, omega, "outfile");

    global_dpd_->buf4_close(&CMNEF); global_dpd_->buf4_close(&CMnEf); global_dpd_->buf4_close(&CmNeF);
    global_dpd_->buf4_close(&WABEI); global_dpd_->buf4_close(&WaBeI); global_dpd_->buf4_close(&WAbEi);
    global_dpd_->buf4_close(&WMBIJ); global_dpd_->buf4_close(&WMbIj); global_dpd_->buf4_close(&WmBiJ);
    global_dpd_->buf4_close(&DIJAB_anti); global_dpd_->buf4_close(&DIjAb);
    global_dpd_->file2_close(&FME); global_dpd_->file2_close(&Fme);
    global_dpd_->buf4_close(&WAMEF); global_dpd_->buf4_close(&WaMeF); global_dpd_->buf4_close(&WAmEf);
    global_dpd_->buf4_close(&WMNIE); global_dpd_->buf4_close(&WMnIe); global_dpd_->buf4_close(&WmNiE);

    /*** beta-beta-alpha term 1 */ 

    sprintf(lbl, "%s %d", "Cmnef", i);
    global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, C_irr, 10, 15, 12, 17, 0, lbl);
    sprintf(lbl, "%s %d", "CMnEf", i);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, C_irr, 22, 28, 22, 28, 0, lbl);
    global_dpd_->buf4_init(&CmNeF, PSIF_EOM_TMP, C_irr, 23, 29, 23, 29, 0, "CmNeF");

    global_dpd_->buf4_init(&Wabei, PSIF_CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wabei (ie,b>a)");
    global_dpd_->buf4_init(&WaBeI, PSIF_CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaBeI (Ie,Ba)");
    global_dpd_->buf4_init(&WAbEi, PSIF_CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAbEi (iE,bA)");
    global_dpd_->buf4_init(&Wmbij, PSIF_CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmbij (i>j,mb)");
    global_dpd_->buf4_init(&WMbIj, PSIF_CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMbIj (Ij,Mb)");
    global_dpd_->buf4_init(&WmBiJ, PSIF_CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmBiJ (iJ,mB)");

    global_dpd_->buf4_init(&Dijab_anti, PSIF_CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
    global_dpd_->buf4_init(&DiJaB, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");

    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
    global_dpd_->buf4_init(&Wamef, PSIF_CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wamef (ma,f>e)");
    global_dpd_->buf4_init(&WaMeF, PSIF_CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaMeF (Ma,Fe)");
    global_dpd_->buf4_init(&WAmEf, PSIF_CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAmEf (mA,fE)");
    global_dpd_->buf4_init(&Wmnie, PSIF_CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmnie (m>n,ie)");
    global_dpd_->buf4_init(&WMnIe, PSIF_CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMnIe (Mn,Ie)");
    global_dpd_->buf4_init(&WmNiE, PSIF_CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmNiE (mN,iE)");

         /* * <S| H    <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_1
            * <D| Hhat <T| (Uhat C2)c   |0> |T> / (w-wt) -> sigma_2 */

    global_dpd_->cc3_sigma_UHF_BBA(&Cmnef, &CMnEf, &CmNeF, &Wabei, &WaBeI, &WAbEi,
      &Wmbij, &WMbIj, &WmBiJ, 1, &Dijab_anti, &DiJaB, &SIA, &Sia,
      1, &FME, &Fme, &Wamef, &WaMeF, &WAmEf, &Wmnie, &WMnIe, &WmNiE,
      &Sijab, &SIjAb, moinfo.aoccpi, moinfo.aocc_off, moinfo.boccpi,
      moinfo.bocc_off, moinfo.avirtpi, moinfo.avir_off, moinfo.bvirtpi,
      moinfo.bvir_off, omega, "outfile");

    global_dpd_->buf4_close(&Cmnef); global_dpd_->buf4_close(&CMnEf); global_dpd_->buf4_close(&CmNeF);
    global_dpd_->buf4_close(&Wabei); global_dpd_->buf4_close(&WaBeI); global_dpd_->buf4_close(&WAbEi);
    global_dpd_->buf4_close(&Wmbij); global_dpd_->buf4_close(&WMbIj); global_dpd_->buf4_close(&WmBiJ);
    global_dpd_->buf4_close(&Dijab_anti); global_dpd_->buf4_close(&DiJaB);
    global_dpd_->file2_close(&FME); global_dpd_->file2_close(&Fme);
    global_dpd_->buf4_close(&Wamef); global_dpd_->buf4_close(&WaMeF); global_dpd_->buf4_close(&WAmEf);
    global_dpd_->buf4_close(&Wmnie); global_dpd_->buf4_close(&WMnIe); global_dpd_->buf4_close(&WmNiE);

#ifdef EOM_DEBUG
  dpd_file2_close(&SIA);
  dpd_file2_close(&Sia);
  dpd_buf4_close(&SIJAB);
  dpd_buf4_close(&Sijab);
  dpd_buf4_close(&SIjAb);
  check_sum("<Psi|Hhat<T|(Uhat C2)c|0>|T>/(w-wt)", i, C_irr);
  sprintf(lbl, "%s %d", "SIA", i);
  dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
  sprintf(lbl, "%s %d", "Sia", i);
  dpd_file2_init(&Sia, EOM_Sia, C_irr, 2, 3, lbl);
  sprintf(lbl, "%s %d", "SIJAB", i);
  dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 0, 5, 2, 7, 0, lbl);
  sprintf(lbl, "%s %d", "Sijab", i);
  dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 10, 15, 12, 17, 0, lbl);
  sprintf(lbl, "%s %d", "SIjAb", i);
  dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, lbl);
#endif

    /*** alpha-alpha-alpha term 2 */

    global_dpd_->buf4_init(&TIJAB, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&WABEI, PSIF_CC3_HC1ET1, C_irr, 20, 5, 20, 7, 0, "Ht_WABEI (IE,B>A)");
    global_dpd_->buf4_init(&WMBIJ, PSIF_CC3_HC1ET1, C_irr, 0, 20, 2, 20, 0, "Ht_WMBIJ (I>J,MB)");
    global_dpd_->buf4_init(&DIJAB_anti, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->buf4_init(&WAMEF, PSIF_CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WAMEF (MA,F>E)");
    global_dpd_->buf4_init(&WMNIE, PSIF_CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMNIE (M>N,IE)");

         /* * <S| H    <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_1
            * <D| Hhat <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_2 */

    global_dpd_->cc3_sigma_UHF_AAA(&TIJAB, &WABEI, &WMBIJ, 1, &DIJAB_anti, &SIA,
        1, &FME, &WAMEF, &WMNIE, &SIJAB, moinfo.aoccpi, moinfo.aocc_off, 
        moinfo.avirtpi, moinfo.avir_off, omega, "outfile");

    global_dpd_->buf4_close(&TIJAB);
    global_dpd_->buf4_close(&WABEI);
    global_dpd_->buf4_close(&WMBIJ);
    global_dpd_->buf4_close(&DIJAB_anti);
    global_dpd_->file2_close(&FME);
    global_dpd_->buf4_close(&WAMEF);
    global_dpd_->buf4_close(&WMNIE);

    /*** beta-beta-beta term 2 */

    global_dpd_->buf4_init(&Tijab, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    global_dpd_->buf4_init(&Wabei, PSIF_CC3_HC1ET1, C_irr, 30, 15, 30, 17, 0, "Ht_Wabei (ie,b>a)");
    global_dpd_->buf4_init(&Wmbij, PSIF_CC3_HC1ET1, C_irr, 10, 30, 12, 30, 0, "Ht_Wmbij (i>j,mb)");
    global_dpd_->buf4_init(&Dijab_anti, PSIF_CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
    global_dpd_->buf4_init(&Wamef, PSIF_CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wamef (ma,f>e)");
    global_dpd_->buf4_init(&Wmnie, PSIF_CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmnie (m>n,ie)");

         /* * <S| H    <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_1
            * <D| Hhat <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_2 */

    global_dpd_->cc3_sigma_UHF_BBB(&Tijab, &Wabei, &Wmbij, 1, &Dijab_anti, &Sia,
        1, &Fme, &Wamef, &Wmnie, &Sijab, moinfo.boccpi, moinfo.bocc_off,
        moinfo.bvirtpi, moinfo.bvir_off, omega, "outfile");

    global_dpd_->buf4_close(&Tijab);
    global_dpd_->buf4_close(&Wabei);
    global_dpd_->buf4_close(&Wmbij);
    global_dpd_->buf4_close(&Dijab_anti);
    global_dpd_->file2_close(&Fme);
    global_dpd_->buf4_close(&Wamef);
    global_dpd_->buf4_close(&Wmnie);

    /*** do alpha-alpha-beta term 2 */

    global_dpd_->buf4_init(&TIJAB, PSIF_CC_TAMPS, 0,  0,  5,  2,  7, 0, "tIJAB");
    global_dpd_->buf4_init(&TIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->buf4_init(&TiJaB, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");

    global_dpd_->buf4_init(&WABEI, PSIF_CC3_HC1ET1, C_irr, 20, 5, 20, 7, 0, "Ht_WABEI (IE,B>A)");
    global_dpd_->buf4_init(&WaBeI, PSIF_CC3_HC1ET1, C_irr, 24, 28, 24, 28, 0, "Ht_WaBeI (Ie,Ba)");
    global_dpd_->buf4_init(&WAbEi, PSIF_CC3_HC1ET1, C_irr, 27, 29, 27, 29, 0, "Ht_WAbEi (iE,bA)");
    global_dpd_->buf4_init(&WMBIJ, PSIF_CC3_HC1ET1, C_irr, 0, 20, 2, 20, 0, "Ht_WMBIJ (I>J,MB)");
    global_dpd_->buf4_init(&WMbIj, PSIF_CC3_HC1ET1, C_irr, 22, 24, 22, 24, 0, "Ht_WMbIj (Ij,Mb)");
    global_dpd_->buf4_init(&WmBiJ, PSIF_CC3_HC1ET1, C_irr, 23, 27, 23, 27, 0, "Ht_WmBiJ (iJ,mB)");

    global_dpd_->buf4_init(&DIJAB_anti, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
    global_dpd_->buf4_init(&DIjAb, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");

    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
    global_dpd_->buf4_init(&WAMEF, PSIF_CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WAMEF (MA,F>E)");
    global_dpd_->buf4_init(&WaMeF, PSIF_CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaMeF (Ma,Fe)");
    global_dpd_->buf4_init(&WAmEf, PSIF_CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAmEf (mA,fE)");
    global_dpd_->buf4_init(&WMNIE, PSIF_CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMNIE (M>N,IE)");
    global_dpd_->buf4_init(&WMnIe, PSIF_CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMnIe (Mn,Ie)");
    global_dpd_->buf4_init(&WmNiE, PSIF_CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmNiE (mN,iE)");

         /* * <S| H    <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_1
            * <D| Hhat <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_2 */

    global_dpd_->cc3_sigma_UHF_AAB(&TIJAB, &TIjAb, &TiJaB, &WABEI, &WaBeI, &WAbEi,
       &WMBIJ, &WMbIj, &WmBiJ, 1,  &DIJAB_anti, &DIjAb, &SIA, &Sia,
       1, &FME, &Fme, &WAMEF, &WaMeF, &WAmEf, &WMNIE, &WMnIe, &WmNiE,
       &SIJAB, &SIjAb, moinfo.aoccpi, moinfo.aocc_off, moinfo.boccpi,
       moinfo.bocc_off, moinfo.avirtpi, moinfo.avir_off, moinfo.bvirtpi,
       moinfo.bvir_off, omega, "outfile");

    global_dpd_->buf4_close(&TIJAB); global_dpd_->buf4_close(&TIjAb); global_dpd_->buf4_close(&TiJaB);
    global_dpd_->buf4_close(&WABEI); global_dpd_->buf4_close(&WaBeI); global_dpd_->buf4_close(&WAbEi);
    global_dpd_->buf4_close(&WMBIJ); global_dpd_->buf4_close(&WMbIj); global_dpd_->buf4_close(&WmBiJ);
    global_dpd_->buf4_close(&DIJAB_anti); global_dpd_->buf4_close(&DIjAb);
    global_dpd_->file2_close(&FME); global_dpd_->file2_close(&Fme);
    global_dpd_->buf4_close(&WAMEF); global_dpd_->buf4_close(&WaMeF); global_dpd_->buf4_close(&WAmEf);
    global_dpd_->buf4_close(&WMNIE); global_dpd_->buf4_close(&WMnIe); global_dpd_->buf4_close(&WmNiE);

    /* beta-beta-alpha term 2 */ 

    global_dpd_->buf4_init(&Tijab, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    global_dpd_->buf4_init(&TIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->buf4_init(&TiJaB, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");

    global_dpd_->buf4_init(&Wabei, PSIF_CC3_HC1ET1, C_irr, 30, 15, 30, 17, 0, "Ht_Wabei (ie,b>a)");
    global_dpd_->buf4_init(&WaBeI, PSIF_CC3_HC1ET1, C_irr, 24, 28, 24, 28, 0, "Ht_WaBeI (Ie,Ba)");
    global_dpd_->buf4_init(&WAbEi, PSIF_CC3_HC1ET1, C_irr, 27, 29, 27, 29, 0, "Ht_WAbEi (iE,bA)");
    global_dpd_->buf4_init(&Wmbij, PSIF_CC3_HC1ET1, C_irr, 10, 30, 12, 30, 0, "Ht_Wmbij (i>j,mb)");
    global_dpd_->buf4_init(&WMbIj, PSIF_CC3_HC1ET1, C_irr, 22, 24, 22, 24, 0, "Ht_WMbIj (Ij,Mb)");
    global_dpd_->buf4_init(&WmBiJ, PSIF_CC3_HC1ET1, C_irr, 23, 27, 23, 27, 0, "Ht_WmBiJ (iJ,mB)");

    global_dpd_->buf4_init(&Dijab_anti, PSIF_CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
    global_dpd_->buf4_init(&DiJaB, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");

    global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
    global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
    global_dpd_->buf4_init(&Wamef, PSIF_CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wamef (ma,f>e)");
    global_dpd_->buf4_init(&WaMeF, PSIF_CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaMeF (Ma,Fe)");
    global_dpd_->buf4_init(&WAmEf, PSIF_CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAmEf (mA,fE)");
    global_dpd_->buf4_init(&Wmnie, PSIF_CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmnie (m>n,ie)");
    global_dpd_->buf4_init(&WMnIe, PSIF_CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMnIe (Mn,Ie)");
    global_dpd_->buf4_init(&WmNiE, PSIF_CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmNiE (mN,iE)");

         /* * <S| H    <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_1
            * <D| Hhat <T| (Utilde T2)c |0> |T> / (w-wt) -> sigma_2 */

    global_dpd_->cc3_sigma_UHF_BBA(&Tijab, &TIjAb, &TiJaB, &Wabei, &WaBeI, &WAbEi,
      &Wmbij, &WMbIj, &WmBiJ, 1, &Dijab_anti, &DiJaB, &SIA, &Sia,
      1, &FME, &Fme, &Wamef, &WaMeF, &WAmEf, &Wmnie, &WMnIe, &WmNiE,
      &Sijab, &SIjAb, moinfo.aoccpi, moinfo.aocc_off, moinfo.boccpi,
      moinfo.bocc_off, moinfo.avirtpi, moinfo.avir_off, moinfo.bvirtpi,
      moinfo.bvir_off, omega, "outfile");

    global_dpd_->buf4_close(&Tijab); global_dpd_->buf4_close(&TIjAb); global_dpd_->buf4_close(&TiJaB);
    global_dpd_->buf4_close(&Wabei); global_dpd_->buf4_close(&WaBeI); global_dpd_->buf4_close(&WAbEi);
    global_dpd_->buf4_close(&Wmbij); global_dpd_->buf4_close(&WMbIj); global_dpd_->buf4_close(&WmBiJ);
    global_dpd_->buf4_close(&Dijab_anti); global_dpd_->buf4_close(&DiJaB);
    global_dpd_->file2_close(&FME); global_dpd_->file2_close(&Fme);
    global_dpd_->buf4_close(&Wamef); global_dpd_->buf4_close(&WaMeF); global_dpd_->buf4_close(&WAmEf);
    global_dpd_->buf4_close(&Wmnie); global_dpd_->buf4_close(&WMnIe); global_dpd_->buf4_close(&WmNiE);

#ifdef EOM_DEBUG
  dpd_file2_close(&SIA);
  dpd_file2_close(&Sia);
  dpd_buf4_close(&SIJAB);
  dpd_buf4_close(&Sijab);
  dpd_buf4_close(&SIjAb);
  check_sum("<Psi|Hhat<T|(Utilde T2)c|0>|T>/(w-wt)", i, C_irr);
  sprintf(lbl, "%s %d", "SIA", i);
  dpd_file2_init(&SIA, EOM_SIA, C_irr, 0, 1, lbl);
  sprintf(lbl, "%s %d", "Sia", i);
  dpd_file2_init(&Sia, EOM_Sia, C_irr, 2, 3, lbl);
  sprintf(lbl, "%s %d", "SIJAB", i);
  dpd_buf4_init(&SIJAB, EOM_SIJAB, C_irr, 0, 5, 2, 7, 0, lbl);
  sprintf(lbl, "%s %d", "Sijab", i);
  dpd_buf4_init(&Sijab, EOM_Sijab, C_irr, 10, 15, 12, 17, 0, lbl);
  sprintf(lbl, "%s %d", "SIjAb", i);
  dpd_buf4_init(&SIjAb, EOM_SIjAb, C_irr, 22, 28, 22, 28, 0, lbl);
#endif

    /*** alpha-alpha-alpha term 3 */

    global_dpd_->buf4_init(&TIJAB, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&WABEI, PSIF_CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WABEI (IE,B>A)");
    global_dpd_->buf4_init(&WMBIJ, PSIF_CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMBIJ (I>J,MB)");
    global_dpd_->file2_init(&FME, PSIF_CC3_HC1, C_irr, 0, 1, "HC1 FME");
    global_dpd_->buf4_init(&WAMEF, PSIF_CC3_HC1, C_irr, 20, 5, 20, 7, 0, "HC1 WAMEF (MA,F>E)");
    global_dpd_->buf4_init(&WMNIE, PSIF_CC3_HC1, C_irr, 0, 20, 2, 20, 0, "HC1 WMNIE (M>N,IE)");

         /* <D| H'   <T| (Uhat T2)c   |0> |T> / (-wt) -> sigma_2 */

    global_dpd_->cc3_sigma_UHF_AAA(&TIJAB, &WABEI, &WMBIJ, 0, NULL, NULL,
        1, &FME, &WAMEF, &WMNIE, &SIJAB, moinfo.aoccpi, moinfo.aocc_off,
        moinfo.avirtpi, moinfo.avir_off, 0.0, "outfile");

    global_dpd_->buf4_close(&TIJAB);
    global_dpd_->buf4_close(&WABEI);
    global_dpd_->buf4_close(&WMBIJ);
    global_dpd_->file2_close(&FME);
    global_dpd_->buf4_close(&WAMEF);
    global_dpd_->buf4_close(&WMNIE);

    /*** beta-beta-beta term 3 */

    global_dpd_->buf4_init(&Tijab, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    global_dpd_->buf4_init(&Wabei, PSIF_CC3_HET1,0, 30, 15, 30, 17, 0, "CC3 Wabei (ie,b>a)");
    global_dpd_->buf4_init(&Wmbij, PSIF_CC3_HET1,0, 10, 30, 12, 30, 0, "CC3 Wmbij (i>j,mb)");
    global_dpd_->file2_init(&Fme, PSIF_CC3_HC1, C_irr, 2, 3, "HC1 Fme");
    global_dpd_->buf4_init(&Wamef, PSIF_CC3_HC1, C_irr, 30, 15, 30, 17, 0, "HC1 Wamef (ma,f>e)");
    global_dpd_->buf4_init(&Wmnie, PSIF_CC3_HC1, C_irr, 10, 30, 12, 30, 0, "HC1 Wmnie (m>n,ie)");

         /* <D| H'   <T| (Uhat T2)c   |0> |T> / (-wt) -> sigma_2 */

    global_dpd_->cc3_sigma_UHF_BBB(&Tijab, &Wabei, &Wmbij, 0, NULL, NULL,
        1, &Fme, &Wamef, &Wmnie, &Sijab, moinfo.boccpi, moinfo.bocc_off, 
        moinfo.bvirtpi, moinfo.bvir_off, 0.0, "outfile");

    global_dpd_->buf4_close(&Tijab);
    global_dpd_->buf4_close(&Wabei);
    global_dpd_->buf4_close(&Wmbij);
    global_dpd_->file2_close(&Fme);
    global_dpd_->buf4_close(&Wamef);
    global_dpd_->buf4_close(&Wmnie);

    /*** alpha-alpha-beta term 3 */

    global_dpd_->buf4_init(&TIJAB, PSIF_CC_TAMPS, 0,  0,  5,  2,  7, 0, "tIJAB");
    global_dpd_->buf4_init(&TIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->buf4_init(&TiJaB, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");

    global_dpd_->buf4_init(&WABEI, PSIF_CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WABEI (IE,B>A)");
    global_dpd_->buf4_init(&WaBeI, PSIF_CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaBeI (Ie,Ba)");
    global_dpd_->buf4_init(&WAbEi, PSIF_CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAbEi (iE,bA)");
    global_dpd_->buf4_init(&WMBIJ, PSIF_CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMBIJ (I>J,MB)");
    global_dpd_->buf4_init(&WMbIj, PSIF_CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMbIj (Ij,Mb)");
    global_dpd_->buf4_init(&WmBiJ, PSIF_CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmBiJ (iJ,mB)");

    global_dpd_->file2_init(&FME, PSIF_CC3_HC1, C_irr, 0, 1, "HC1 FME");
    global_dpd_->file2_init(&Fme, PSIF_CC3_HC1, C_irr, 2, 3, "HC1 Fme");
    global_dpd_->buf4_init(&WAMEF, PSIF_CC3_HC1, C_irr, 20, 5, 20, 7, 0, "HC1 WAMEF (MA,F>E)");
    global_dpd_->buf4_init(&WaMeF, PSIF_CC3_HC1, C_irr, 24, 28, 24, 28, 0, "HC1 WaMeF (Ma,Fe)");
    global_dpd_->buf4_init(&WAmEf, PSIF_CC3_HC1, C_irr, 27, 29, 27, 29, 0, "HC1 WAmEf (mA,fE)");
    global_dpd_->buf4_init(&WMNIE, PSIF_CC3_HC1, C_irr, 0, 20, 2, 20, 0, "HC1 WMNIE (M>N,IE)");
    global_dpd_->buf4_init(&WMnIe, PSIF_CC3_HC1, C_irr, 22, 24, 22, 24, 0, "HC1 WMnIe (Mn,Ie)");
    global_dpd_->buf4_init(&WmNiE, PSIF_CC3_HC1, C_irr, 23, 27, 23, 27, 0, "HC1 WmNiE (mN,iE)");

         /* <D| H'   <T| (Uhat T2)c   |0> |T> / (-wt) -> sigma_2 */

    global_dpd_->cc3_sigma_UHF_AAB(&TIJAB, &TIjAb, &TiJaB, &WABEI, &WaBeI, &WAbEi,
       &WMBIJ, &WMbIj, &WmBiJ, 0, NULL, NULL, NULL, NULL,
       1, &FME, &Fme, &WAMEF, &WaMeF, &WAmEf, &WMNIE, &WMnIe, &WmNiE,
       &SIJAB, &SIjAb, moinfo.aoccpi, moinfo.aocc_off, moinfo.boccpi,
       moinfo.bocc_off, moinfo.avirtpi, moinfo.avir_off, moinfo.bvirtpi,
       moinfo.bvir_off, 0.0, "outfile");

    global_dpd_->buf4_close(&TIJAB); global_dpd_->buf4_close(&TIjAb); global_dpd_->buf4_close(&TiJaB);
    global_dpd_->buf4_close(&WABEI); global_dpd_->buf4_close(&WaBeI); global_dpd_->buf4_close(&WAbEi);
    global_dpd_->buf4_close(&WMBIJ); global_dpd_->buf4_close(&WMbIj); global_dpd_->buf4_close(&WmBiJ);
    global_dpd_->file2_close(&FME); global_dpd_->file2_close(&Fme);
    global_dpd_->buf4_close(&WAMEF); global_dpd_->buf4_close(&WaMeF); global_dpd_->buf4_close(&WAmEf);
    global_dpd_->buf4_close(&WMNIE); global_dpd_->buf4_close(&WMnIe); global_dpd_->buf4_close(&WmNiE);

    /*** beta-beta-alpha term 3 */ 

    global_dpd_->buf4_init(&Tijab, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    global_dpd_->buf4_init(&TIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->buf4_init(&TiJaB, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");

    global_dpd_->buf4_init(&Wabei, PSIF_CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wabei (ie,b>a)");
    global_dpd_->buf4_init(&WaBeI, PSIF_CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaBeI (Ie,Ba)");
    global_dpd_->buf4_init(&WAbEi, PSIF_CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAbEi (iE,bA)");
    global_dpd_->buf4_init(&Wmbij, PSIF_CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmbij (i>j,mb)");
    global_dpd_->buf4_init(&WMbIj, PSIF_CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMbIj (Ij,Mb)");
    global_dpd_->buf4_init(&WmBiJ, PSIF_CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmBiJ (iJ,mB)");

    global_dpd_->file2_init(&FME, PSIF_CC3_HC1, C_irr, 0, 1, "HC1 FME");
    global_dpd_->file2_init(&Fme, PSIF_CC3_HC1, C_irr, 2, 3, "HC1 Fme");
    global_dpd_->buf4_init(&Wamef, PSIF_CC3_HC1, C_irr, 30, 15, 30, 17, 0, "HC1 Wamef (ma,f>e)");
    global_dpd_->buf4_init(&WaMeF, PSIF_CC3_HC1, C_irr, 24, 28, 24, 28, 0, "HC1 WaMeF (Ma,Fe)");
    global_dpd_->buf4_init(&WAmEf, PSIF_CC3_HC1, C_irr, 27, 29, 27, 29, 0, "HC1 WAmEf (mA,fE)");
    global_dpd_->buf4_init(&Wmnie, PSIF_CC3_HC1, C_irr, 10, 30, 12, 30, 0, "HC1 Wmnie (m>n,ie)");
    global_dpd_->buf4_init(&WMnIe, PSIF_CC3_HC1, C_irr, 22, 24, 22, 24, 0, "HC1 WMnIe (Mn,Ie)");
    global_dpd_->buf4_init(&WmNiE, PSIF_CC3_HC1, C_irr, 23, 27, 23, 27, 0, "HC1 WmNiE (mN,iE)");

         /* <D| H'   <T| (Uhat T2)c   |0> |T> / (-wt) -> sigma_2 */

    global_dpd_->cc3_sigma_UHF_BBA(&Tijab, &TIjAb, &TiJaB, &Wabei, &WaBeI, &WAbEi,
      &Wmbij, &WMbIj, &WmBiJ, 0, NULL, NULL, NULL, NULL,
      1, &FME, &Fme, &Wamef, &WaMeF, &WAmEf, &Wmnie, &WMnIe, &WmNiE,
      &Sijab, &SIjAb, moinfo.aoccpi, moinfo.aocc_off, moinfo.boccpi,
      moinfo.bocc_off, moinfo.avirtpi, moinfo.avir_off, moinfo.bvirtpi,
      moinfo.bvir_off, 0.0, "outfile");

    global_dpd_->buf4_close(&Tijab); global_dpd_->buf4_close(&TIjAb); global_dpd_->buf4_close(&TiJaB);
    global_dpd_->buf4_close(&Wabei); global_dpd_->buf4_close(&WaBeI); global_dpd_->buf4_close(&WAbEi);
    global_dpd_->buf4_close(&Wmbij); global_dpd_->buf4_close(&WMbIj); global_dpd_->buf4_close(&WmBiJ);
    global_dpd_->file2_close(&FME); global_dpd_->file2_close(&Fme);
    global_dpd_->buf4_close(&Wamef); global_dpd_->buf4_close(&WaMeF); global_dpd_->buf4_close(&WAmEf);
    global_dpd_->buf4_close(&Wmnie); global_dpd_->buf4_close(&WMnIe); global_dpd_->buf4_close(&WmNiE);

    global_dpd_->file2_close(&SIA);
    global_dpd_->file2_close(&Sia);
    global_dpd_->buf4_close(&SIJAB);
    global_dpd_->buf4_close(&Sijab);
    global_dpd_->buf4_close(&SIjAb);

#ifdef EOM_DEBUG
  check_sum("<Psi|H'<T|(Uhat T2)c|0>|T>/(w-wt)", i, C_irr);
#endif
  }
  return;
}

}} // namespace psi::cceom
